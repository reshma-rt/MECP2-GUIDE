#!/usr/bin/env python3
"""Explicit, approval-gated Tier-1 composite ranking. Standard library only.
No robustness, BLAST, folding, assembly, promotion, subprocess or dynamic code execution.
"""
import argparse
import contextlib
import csv
import ctypes
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import sys
from collections import Counter

VERSION = 'MECP2_EXPLICIT_RANKING_ADAPTER_V1'
WEIGHTS = [0.25, 0.30, 0.30, 0.10, 0.05]
MUTATIONS = {'143546', '143559', '143597', '2505626', '3237211', '372650', '4279604', '548706'}
COLS = {
    'a': 'Target_A_unpaired_probability',
    'a5': 'Target_A_5nt_mean_unpaired_probability',
    'b': 'Bystander_precision_score_0_1',
    's': 'Fullwindow_GuideAware_transcriptome_specificity_score_0_1',
    'd': 'Guide_target_duplex_MFE_kcal_mol',
    'g': 'Guide_MFE_kcal_mol',
    'tie': 'Fullwindow_editable_A_offtarget_hits_le3mm',
}
COMP = ['Stage2E_v2_accessibility_component_0_1', 'Stage2E_v2_bystander_component_0_1',
        'Stage2E_v2_specificity_component_0_1', 'Stage2E_v2_duplex_component_0_1',
        'Stage2E_v2_self_structure_component_0_1']
SCORE = 'Stage2E_v2_composite_score_0_100'
RANK = 'Stage2E_v2_rank_within_mutation'
WEIGHT_FIELDS = ['Stage2E_v2_weight_' + x for x in
                 ['accessibility', 'bystander', 'specificity', 'duplex', 'self_structure']]
NEW_FIELDS = COMP + [SCORE, RANK] + WEIGHT_FIELDS + ['Stage2E_v2_status']
SUMMARY_FIELDS = ['VariationID', 'Stage2E_v2_top_Guide_ID', 'Stage2E_v2_top_score_0_100',
                  'Editable_A_offtarget_hits_le3mm']
ROLES = {'adapter', 'schema', 'verification_report', 'ranking_source', 'robustness_source',
         'candidates', 'structure', 'bystander', 'features', 'feature_summary', 'feature_manifest'}
OUTPUTS = {'ranked.csv', 'top_guides.csv', 'mutation_summary.csv', 'normalization.json',
           'run_manifest.json', 'output_manifest.json', 'run_result.json'}

def digest(data):
    return hashlib.sha256(data).hexdigest()

def require(condition, message):
    if not condition:
        raise ValueError(message)

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='microseconds')

def strict_json(data):
    def unique(pairs):
        result = {}
        for k, v in pairs:
            require(k not in result, 'Duplicate JSON key: ' + k)
            result[k] = v
        return result
    return json.loads(data, object_pairs_hook=unique,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON: ' + x)))

def table(data, required=None):
    reader = csv.DictReader(io.StringIO(data.decode('utf-8-sig'), newline=''))
    fields = reader.fieldnames or []
    require(len(fields) == len(set(fields)), 'Duplicate/ambiguous CSV columns')
    require(all(fields), 'Empty header')
    if required is not None:
        require(fields == required, 'Exact schema mismatch; no aliases or fallback allowed')
    rows = list(reader)
    require(all(None not in r and all(v is not None for v in r.values()) for r in rows), 'Malformed CSV row')
    return rows, fields

def number(value, unit=False):
    x = float(value)
    require(math.isfinite(x), 'Nonfinite input')
    require(not unit or 0 <= x <= 1, 'Value outside [0,1]; no clamping')
    return x

def integer(value):
    require(isinstance(value, str) and bool(re.fullmatch(r'\d+', value)), 'Expected nonnegative integer string')
    return int(value)

def indexed(rows):
    require(all(r.get('Guide_ID') for r in rows), 'Missing Guide_ID')
    require(len(rows) == len({r['Guide_ID'] for r in rows}), 'Duplicate Guide_ID')
    return {r['Guide_ID']: r for r in rows}

def cohort(rows):
    ids = indexed(rows)
    require(set(ids) == {'T1G%04d' % i for i in range(1, 377)}, 'Incomplete/extra frozen guide cohort')
    counts = Counter(r['VariationID'] for r in rows)
    require(set(counts) == MUTATIONS and set(counts.values()) == {47}, 'Incomplete/mismatched mutation groups')
    return ids

def validate_records(candidates, structure, bystander, features):
    maps = [cohort(rows) for rows in (candidates, structure, bystander, features)]
    for gid in sorted(maps[0]):
        c, s, b, f = [m[gid] for m in maps]
        # Preserve full C lineage, and every frozen B value except recorded stage metadata.
        require(all(f[k] == v for k, v in b.items()), gid + ': Stage2C lineage mismatch')
        require(all(f[k] == v for k, v in s.items()
                    if k not in {'Dataset2_stage', 'Bystander_risk_model_status'}), gid + ': structure mismatch')
        for k in ['VariationID', 'Transcript', 'Guide_sequence_RNA_5to3', 'Target_segment_RNA_5to3',
                  'Target_A_transcript_position', 'Guide_length_nt', 'Target_A_position_in_guide_1based',
                  'Target_A_position_in_target_segment_1based', 'Rescue_class']:
            require(c[k] == f[k], gid + ': frozen candidate identity mismatch: ' + k)
        length = integer(f['Guide_length_nt'])
        gp = integer(f['Target_A_position_in_guide_1based'])
        qp = length - gp + 1
        require(length in {21, 25, 31} and 6 <= gp <= length - 5, gid + ': invalid guide-relative position')
        require(integer(f['Target_A_position_in_target_segment_1based']) == qp
                == integer(f['Fullwindow_catalytic_query_position_1based']), gid + ': coordinate conflict')
        guide, target = f['Guide_sequence_RNA_5to3'], f['Target_segment_RNA_5to3']
        require(len(guide) == len(target) == length and set(guide + target) <= set('ACGU'), gid + ': RNA sequence invalid')
        require(guide[gp-1] == 'C' and target[qp-1] == 'A', gid + ': catalytic base mismatch')
        require(f['Stage2D_v3_status'] == 'FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE', gid + ': upstream not complete')
        require(f['Fullwindow_intended_transcript_hit_le3mm'] == 'YES'
                and integer(f['Fullwindow_intended_best_mismatch_count']) == 0, gid + ': intended-hit QC failed')
        for k in ['a', 'a5', 'd', 'g']:
            number(s[COLS[k]], k in {'a', 'a5'})
        for k in ['b', 's']:
            number(f[COLS[k]], True)
        integer(f[COLS['tie']])
    return maps

def normalize(values, lower=False):
    lo, hi = min(values), max(values)
    scaled = [1.0 if hi == lo else (v - lo) / (hi - lo) for v in values]
    return ([1.0 - x for x in scaled] if lower else scaled), [lo, hi]

def rank_rows(features, structure):
    """Pure scoring function. Called on real data ONLY by the approved ranking stage."""
    s = indexed(structure)
    rows = [dict(r) for r in sorted(features, key=lambda r: r['Guide_ID'])]
    ds, drange = normalize([number(s[r['Guide_ID']][COLS['d']]) for r in rows], lower=True)
    gs, grange = normalize([number(s[r['Guide_ID']][COLS['g']]) for r in rows])
    for i, r in enumerate(rows):
        st = s[r['Guide_ID']]
        comps = [(number(st[COLS['a']], True) + number(st[COLS['a5']], True)) / 2,
                 number(r[COLS['b']], True), number(r[COLS['s']], True), ds[i], gs[i]]
        for field, value in zip(COMP, comps):
            r[field] = value
        # Explicit addition order matches inspected Stage2E v2.
        a, b, sp, d, g = comps
        r[SCORE] = (0.25*a + 0.30*b + 0.30*sp + 0.10*d + 0.05*g) * 100.0
        for field, value in zip(WEIGHT_FIELDS, WEIGHTS):
            r[field] = value
        r['Stage2E_v2_status'] = 'HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE'
    top = []
    for vid in sorted({r['VariationID'] for r in rows}):
        group = sorted([r for r in rows if r['VariationID'] == vid],
                       key=lambda r: (-r[SCORE], integer(r[COLS['tie']]), r['Guide_ID']))
        for rank, r in enumerate(group, 1):
            r[RANK] = rank
        top.append(group[0])
    rows.sort(key=lambda r: (r['VariationID'], r[RANK], r['Guide_ID']))
    summary = [{'VariationID': r['VariationID'], 'Stage2E_v2_top_Guide_ID': r['Guide_ID'],
                'Stage2E_v2_top_score_0_100': r[SCORE],
                'Editable_A_offtarget_hits_le3mm': r[COLS['tie']]} for r in top]
    return rows, top, summary, {'duplex_MFE_range': drange, 'guide_MFE_range': grange,
                               'cohort_size': len(rows), 'constant_duplex_component': 0.0,
                               'constant_self_component': 1.0}

@contextlib.contextmanager
def locked_input(path):
    # Windows share-read only: deny concurrent write/delete while processing.
    # This is cooperative process protection, not an OS sandbox or production deployment.
    require(path.is_file() and not path.is_symlink(), 'Input unavailable or symlink: ' + str(path))
    for ancestor in [path, *path.parents]:
        st = ancestor.lstat()
        require(not getattr(st, 'st_file_attributes', 0) & 0x400, 'Reparse path rejected: ' + str(ancestor))
    if os.name == 'nt':
        import msvcrt
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.CreateFileW.argtypes = [ctypes.c_wchar_p, ctypes.c_uint32, ctypes.c_uint32,
                                      ctypes.c_void_p, ctypes.c_uint32, ctypes.c_uint32, ctypes.c_void_p]
        kernel.CreateFileW.restype = ctypes.c_void_p
        handle = kernel.CreateFileW(str(path), 0x80000000, 1, None, 3, 0x80, None)
        require(handle != ctypes.c_void_p(-1).value, 'Cannot lock read-only input: ' + str(path))
        fd = msvcrt.open_osfhandle(handle, os.O_RDONLY | os.O_BINARY)
        stream = os.fdopen(fd, 'rb')
    else:
        stream = path.open('rb')
    with stream:
        yield stream

def write_json(path, obj):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')

def write_csv(path, rows, fields):
    with path.open('x', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='raise')
        w.writeheader()
        w.writerows(rows)

def approval_check(approval, config_hash, run, now=None):
    require(approval.get('approved') is True and bool(approval.get('approved_by')), 'Approval not granted/named')
    require(approval.get('stage') == 'composite-ranking', 'Wrong approved stage')
    require(approval.get('config_sha256') == config_hash and approval.get('run_directory') == str(run), 'Approval binding mismatch')
    require(approval.get('python_sha256') == digest(Path(sys.executable).read_bytes()), 'Approval interpreter hash mismatch')
    expiry = dt.datetime.fromisoformat(approval.get('expires_utc', '').replace('Z', '+00:00'))
    require(expiry.tzinfo is not None and expiry > (now or dt.datetime.now(dt.timezone.utc)), 'Approval expired/invalid')
    require(approval.get('robustness_approved') is False and approval.get('promotion_approved') is False,
            'Approval must restrict downstream stages')

def execute(stage, config_path, expected_hash, approval_path=None):
    require(stage in {'verify-inputs', 'composite-ranking'}, 'Stage is not implemented or permitted')
    config_path = Path(config_path).resolve(strict=True)
    started = utc()
    run = None
    with contextlib.ExitStack() as stack:
        config_stream = stack.enter_context(locked_input(config_path))
        config_bytes = config_stream.read()
        require(digest(config_bytes) == expected_hash, 'Stale configuration hash')
        cfg = strict_json(config_bytes)
        require(cfg['version'] == VERSION and cfg['columns'] == COLS and cfg['weights'] == WEIGHTS,
                'Unsupported contract/columns/weights')
        require(set(cfg['inputs']) == ROLES, 'Missing/extra input binding')
        require(cfg['stages'] == {'composite-ranking': 'PROPOSED', 'robustness': 'NOT_AUTHORIZED',
                                 'final-library': 'NOT_AUTHORIZED', 'promotion': 'NOT_AUTHORIZED'}, 'Stage separation violated')
        root = Path(cfg['staging_root'])
        require(root.is_absolute() and root.is_dir(), 'Unverified staging root')
        require(not root.is_symlink(), 'Staging root is symlink')
        for a in [root, *root.parents]:
            require(not getattr(a.lstat(), 'st_file_attributes', 0) & 0x400 and not a.is_symlink(), 'Unsafe staging path')
        require(bool(re.fullmatch(r'[A-Za-z0-9_-]+', cfg['run_id'])), 'Unsafe run ID')
        proposed = root / cfg['run_id']
        require(not proposed.exists(), 'Run directory already exists; no overwrite or reuse')
        for official in cfg['protected_roots']:
            off = Path(official).resolve()
            require(not proposed.resolve().is_relative_to(off), 'Prohibited official output path')
        streams, payloads, paths = {}, {}, {}
        for role, binding in cfg['inputs'].items():
            require(binding.get('path') and re.fullmatch('[a-f0-9]{64}', binding.get('sha256', '')), 'Unassigned input: ' + role)
            p = Path(binding['path'])
            if not p.is_absolute():
                require(role in {'adapter', 'schema', 'verification_report', 'ranking_source', 'robustness_source'}, 'Data path must be absolute')
                p = config_path.parent / p
            p = p.absolute()
            require(not p.is_relative_to(proposed), 'Input is under run output directory')
            streams[role] = stack.enter_context(locked_input(p))
            payloads[role] = streams[role].read()
            require(digest(payloads[role]) == binding['sha256'], 'Stale input hash: ' + role)
            paths[role] = str(p)
        require(len(set(paths.values())) == len(paths), 'Conflicting input paths')
        require(Path(paths['adapter']).resolve() == Path(__file__).resolve(), 'Executing adapter binding mismatch')
        schema = strict_json(payloads['schema'])
        parsed = {role: table(payloads[role], schema[role]) for role in
                  ['candidates', 'structure', 'bystander', 'features', 'feature_summary']}
        validate_records(*[parsed[r][0] for r in ['candidates', 'structure', 'bystander', 'features']])
        summary = parsed['feature_summary'][0]
        require(len(summary) == 8 and {r['VariationID'] for r in summary} == MUTATIONS, 'Feature summary incomplete')
        require(all(integer(r['Guide_count']) == integer(r['Valid_intended_fullwindow_le3mm']) == 47 for r in summary), 'Feature summary counts fail')
        fm = strict_json(payloads['feature_manifest'])
        # Manifest entries must explicitly bind the two supplied feature outputs.
        entries = fm.get('files', []) if isinstance(fm, dict) else []
        for role in ['features', 'feature_summary']:
            matches = [r for r in entries if r.get('sha256') == cfg['inputs'][role]['sha256']]
            require(len(matches) == 1, 'Feature manifest hash association missing/ambiguous: ' + role)
        qc = {'guides': 376, 'mutations': 8, 'guides_per_mutation': 47,
              'exact_column_bindings': COLS, 'scoring_performed': False}
        if stage == 'verify-inputs':
            # Read-only preflight; does not create/consume the proposed ranking directory.
            for role, stream in streams.items():
                stream.seek(0)
                require(digest(stream.read()) == cfg['inputs'][role]['sha256'], 'Preservation failed: ' + role)
            config_stream.seek(0)
            require(digest(config_stream.read()) == expected_hash, 'Config preservation failed')
            return {'status': 'VERIFIED_FOR_APPROVAL_REVIEW', 'input_qc': qc,
                    'config_sha256': expected_hash, 'proposed_run_directory': str(proposed),
                    'python_executable': sys.executable, 'python_version': sys.version,
                    'python_sha256': digest(Path(sys.executable).read_bytes()),
                    'production_host_enforcement_verified': False, 'exit_code': 0}
        require(approval_path is not None, 'Separate ranking approval required')
        approval_path = Path(approval_path).resolve(strict=True)
        approval_stream = stack.enter_context(locked_input(approval_path))
        ab = approval_stream.read()
        approval_check(strict_json(ab), expected_hash, proposed)
        try:
            proposed.mkdir(exist_ok=False)
            run = proposed
            write_json(run / 'run_manifest.json', {'version': VERSION, 'stage': stage, 'started_utc': started,
                       'config_sha256': expected_hash, 'config': cfg, 'input_paths': paths,
                       'approval_sha256': digest(ab), 'python_executable': sys.executable,
                       'python_version': sys.version, 'python_sha256': digest(Path(sys.executable).read_bytes()),
                       'historical_provenance': 'UNRESOLVED', 'production_host_enforcement_verified': False})
            ranked, top, summary, norm = rank_rows(parsed['features'][0], parsed['structure'][0])
            require(len(ranked) == 376 and len(top) == len(summary) == 8, 'Output row count failed')
            for vid in MUTATIONS:
                require(sorted(r[RANK] for r in ranked if r['VariationID'] == vid) == list(range(1, 48)), 'Rank completeness failed')
            fields = parsed['features'][1] + NEW_FIELDS
            write_csv(run / 'ranked.csv', ranked, fields)
            write_csv(run / 'top_guides.csv', top, fields)
            write_csv(run / 'mutation_summary.csv', summary, SUMMARY_FIELDS)
            write_json(run / 'normalization.json', norm)
            for name, expected_rows in [('ranked.csv', ranked), ('top_guides.csv', top), ('mutation_summary.csv', summary)]:
                actual, headers = table((run / name).read_bytes())
                expected = [{k: str(v) for k, v in r.items()} for r in expected_rows]
                require(actual == expected, 'Written output content/QC failure: ' + name)
            for role, stream in streams.items():
                stream.seek(0)
                require(digest(stream.read()) == cfg['inputs'][role]['sha256'], 'Preservation failed: ' + role)
            config_stream.seek(0)
            require(digest(config_stream.read()) == expected_hash, 'Config preservation failed')
            approval_stream.seek(0)
            require(digest(approval_stream.read()) == digest(ab), 'Approval preservation failed')
            entries = [{'path': p.name, 'bytes': p.stat().st_size, 'sha256': digest(p.read_bytes())}
                       for p in sorted(run.iterdir())]
            write_json(run / 'output_manifest.json', {'files': entries})
            result = {'status': 'COMPLETE', 'stage': stage, 'rows': 376, 'top_rows': 8,
                      'summary_rows': 8, 'preservation': 'PASS', 'finished_utc': utc(),
                      'output_manifest_sha256': digest((run / 'output_manifest.json').read_bytes()),
                      'downstream_stages_executed': [], 'production_host_enforcement_verified': False, 'exit_code': 0}
            require({p.name for p in run.iterdir()} == OUTPUTS - {'run_result.json'}, 'Unexpected/missing output file')
            write_json(run / 'run_result.json', result)
            return result
        except Exception as e:
            if run is not None and not (run / 'run_result.json').exists():
                write_json(run / 'run_result.json', {'status': 'FAILED', 'reason': str(e), 'exit_code': 4,
                                                    'partial_outputs_not_approved': True, 'finished_utc': utc()})
            raise

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['verify-inputs', 'composite-ranking'])
    p.add_argument('--config', required=True)
    p.add_argument('--config-sha256', required=True)
    p.add_argument('--approval')
    args = p.parse_args(argv)
    try:
        result = execute(args.stage, args.config, args.config_sha256, args.approval)
        print(json.dumps(result, indent=2))
        return 0
    except Exception as e:
        print(json.dumps({'status': 'FAILED_CLOSED', 'reason': str(e), 'stage': args.stage, 'exit_code': 4}, indent=2))
        return 4

if __name__ == '__main__':
    sys.exit(main())
