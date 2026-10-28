"""Synthetic fixtures ONLY. No real guide scores/winners/robustness are calculated."""
import ast
import contextlib
import copy
import csv
import datetime as dt
import io
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch
import mecp2_ranking_adapter as m

HERE = Path(__file__).resolve().parent

def fixtures():
    schema = json.loads((HERE/'schema_contract.json').read_text())
    sets = {k: [] for k in ['candidates', 'structure', 'bystander', 'features']}
    for i in range(376):
        gid = 'T1G%04d' % (i+1)
        vid = sorted(m.MUTATIONS)[i//47]
        common = dict(Guide_ID=gid, VariationID=vid, Transcript='SYNTHETIC_ONLY',
                      Rescue_class='SYNTHETIC_EXACT_RESTORATION', Guide_length_nt='21',
                      Target_A_transcript_position='1000', Target_A_position_in_guide_1based='6',
                      Target_A_position_in_target_segment_1based='16',
                      Guide_sequence_RNA_5to3='UUUUUCUUUUUUUUUUUUUUU', Target_segment_RNA_5to3='A'*21)
        s = {k: common.get(k, '') for k in schema['structure']}
        s.update({m.COLS['a']: str((i%9)/10), m.COLS['a5']: '0.4',
                  m.COLS['d']: str(-20-i%7), m.COLS['g']: str(-8+i%5)})
        b = {k: s.get(k, '') for k in schema['bystander']}
        b[m.COLS['b']] = '0.7'
        f = {k: b.get(k, '') for k in schema['features']}
        f.update({m.COLS['s']: str((i%11)/11), m.COLS['tie']: str(i%3),
                  'Stage2D_v3_status': 'FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE',
                  'Fullwindow_intended_transcript_hit_le3mm': 'YES', 'Fullwindow_intended_best_mismatch_count': '0',
                  'Fullwindow_catalytic_query_position_1based': '16'})
        c = {k: common.get(k, '') for k in schema['candidates']}
        for role, row in [('candidates', c), ('structure', s), ('bystander', b), ('features', f)]:
            sets[role].append(row)
    return sets, schema

class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.rows, self.schema = fixtures()

    def validate(self):
        return m.validate_records(*[self.rows[x] for x in ['candidates', 'structure', 'bystander', 'features']])

    def test_valid_cohort(self):
        self.validate()

    def test_duplicate_guide(self):
        self.rows['features'][1]['Guide_ID'] = self.rows['features'][0]['Guide_ID']
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.validate()

    def test_incomplete_cohort(self):
        self.rows['features'].pop()
        with self.assertRaisesRegex(ValueError, 'cohort'):
            self.validate()

    def test_incomplete_group(self):
        self.rows['features'][0]['VariationID'] = self.rows['features'][-1]['VariationID']
        with self.assertRaisesRegex(ValueError, 'mutation groups'):
            self.validate()

    def test_nonfinite_and_missing(self):
        for v in ['NaN', 'NAN', 'inf', '-Infinity', '', 'NA']:
            with self.subTest(v=v), self.assertRaises(ValueError):
                m.number(v)

    def test_no_clamping(self):
        for v in ['-0.001', '1.0000000001']:
            with self.assertRaises(ValueError): m.number(v, True)

    def test_missing_exact_column(self):
        with self.assertRaisesRegex(ValueError, 'schema'):
            m.table(b'Guide_ID,Specificity_component_0_1\nG,1\n', ['Guide_ID', m.COLS['s']])

    def test_duplicate_header(self):
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            m.table(b'Guide_ID,Guide_ID\nA,A\n')

    def test_malformed_row(self):
        for b in [b'A,B\n1\n', b'A\n1,2\n']:
            with self.assertRaisesRegex(ValueError, 'Malformed'): m.table(b)

    def test_duplicate_config_keys(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate JSON'): m.strict_json('{"x":1,"x":2}')

    def test_coordinate_conflict(self):
        self.rows['features'][0]['Fullwindow_catalytic_query_position_1based'] = '20'
        with self.assertRaisesRegex(ValueError, 'coordinate conflict'): self.validate()

    def test_upstream_incomplete(self):
        self.rows['features'][0]['Stage2D_v3_status'] = 'INCOMPLETE'
        with self.assertRaisesRegex(ValueError, 'upstream'): self.validate()

    def test_intended_qc(self):
        self.rows['features'][0]['Fullwindow_intended_best_mismatch_count'] = '1'
        with self.assertRaisesRegex(ValueError, 'intended-hit'): self.validate()

    def test_intended_is_yes_flag_not_integer_count(self):
        self.rows['features'][0]['Fullwindow_intended_transcript_hit_le3mm'] = '1'
        with self.assertRaisesRegex(ValueError, 'intended-hit'): self.validate()

    def test_structure_mismatch(self):
        self.rows['structure'][0][m.COLS['d']] = '-99'
        with self.assertRaisesRegex(ValueError, 'structure mismatch'): self.validate()

    def test_bystander_lineage(self):
        self.rows['features'][0][m.COLS['b']] = '0.1'
        with self.assertRaisesRegex(ValueError, 'lineage'): self.validate()

    def test_known_score(self):
        # All components [0.5,0.4,0.8,1,1] -> 63.5, synthetic only.
        f = self.rows['features'][:2]; s = self.rows['structure'][:2]
        for r in s: r.update({m.COLS['a']: '0.4', m.COLS['a5']: '0.6'})
        for r in f: r.update({m.COLS['b']: '0.4', m.COLS['s']: '0.8'})
        s[0].update({m.COLS['d']: '-30', m.COLS['g']: '-1'})
        s[1].update({m.COLS['d']: '-20', m.COLS['g']: '-9'})
        ranked, top, _, _ = m.rank_rows(f, s)
        self.assertAlmostEqual(top[0][m.SCORE], 63.5)
        self.assertEqual([r[m.RANK] for r in ranked], [1, 2])

    def test_constant_energy_convention(self):
        self.assertEqual(m.normalize([-5, -5], True)[0], [0, 0])
        self.assertEqual(m.normalize([-5, -5])[0], [1, 1])

    def test_tie_zero_count_then_id(self):
        f = self.rows['features'][:3]; s = self.rows['structure'][:3]
        for r in f: r.update({m.COLS['b']: '0.5', m.COLS['s']: '0.5'})
        for r in s: r.update({m.COLS['a']: '0.5', m.COLS['a5']: '0.5', m.COLS['d']: '-20', m.COLS['g']: '-5'})
        f[0][m.COLS['tie']] = '1'; f[1][m.COLS['tie']] = f[2][m.COLS['tie']] = '0'
        ranked, _, _, _ = m.rank_rows(f, s)
        self.assertEqual([r['Guide_ID'] for r in ranked], ['T1G0002', 'T1G0003', 'T1G0001'])

    def test_deterministic_shuffled_input(self):
        a = m.rank_rows(self.rows['features'], self.rows['structure'])
        random.Random(8).shuffle(self.rows['features'])
        random.Random(9).shuffle(self.rows['structure'])
        self.assertEqual(a, m.rank_rows(self.rows['features'], self.rows['structure']))

    def test_no_legacy_status_selection(self):
        for rows in self.rows.values():
            for r in rows: r['Transcriptome_offtarget_status'] = 'NOT_RUN'
        self.validate()

    def test_source_formula_parity_synthetic(self):
        # Execute the original E v2 source only against a fresh synthetic 376-row fixture.
        # Redirect only ROOT binding in its AST, retaining all scoring/ranking statements.
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); pipe = root/'02_tier1_pipeline'; pipe.mkdir()
            m.write_csv(pipe/'MECP2_ADAR_Dataset2_Tier1_stage2B_structure.csv', self.rows['structure'], self.schema['structure'])
            m.write_csv(pipe/'MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv', self.rows['features'], self.schema['features'])
            source = HERE/'source_snapshots/stage2E_v2_tier1_harmonized_ranking.py'
            tree = ast.parse(source.read_text())
            for node in tree.body:
                if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'ROOT' for t in node.targets):
                    node.value = ast.Call(func=ast.Name(id='Path', ctx=ast.Load()), args=[ast.Constant(str(root))], keywords=[])
            ast.fix_missing_locations(tree)
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(tree, str(source), 'exec'), {'__file__': str(source)})
            original, _ = m.table((pipe/'MECP2_ADAR_Dataset2_Tier1_stage2E_v2_ranked.csv').read_bytes())
            ours = m.rank_rows(self.rows['features'], self.rows['structure'])[0]
            self.assertEqual([r['Guide_ID'] for r in original], [r['Guide_ID'] for r in ours])
            for old, new in zip(original, ours):
                for col in m.COMP + [m.SCORE, m.RANK]: self.assertEqual(float(old[col]), new[col])

class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name); self.stage = self.root/'stage'; self.stage.mkdir()
        self.rows, schema = fixtures()
        self.cfg = {'version': m.VERSION, 'columns': dict(m.COLS), 'weights': list(m.WEIGHTS),
                    'staging_root': str(self.stage), 'run_id': 'SYNTHETIC_ONLY', 'protected_roots': [str(self.root/'official')],
                    'stages': {'composite-ranking':'PROPOSED', 'robustness':'NOT_AUTHORIZED',
                               'final-library':'NOT_AUTHORIZED', 'promotion':'NOT_AUTHORIZED'}, 'inputs': {}}
        paths = {'adapter': HERE/'mecp2_ranking_adapter.py', 'schema': HERE/'schema_contract.json',
                 'ranking_source': HERE/'source_snapshots/stage2E_v2_tier1_harmonized_ranking.py',
                 'robustness_source': HERE/'source_snapshots/stage2F_v2_tier1_harmonized_sensitivity.py',
                 'verification_report': HERE/'source_snapshots/MECP2_Feature_Output_Verification.md'}
        for role in self.rows:
            p = self.root/(role+'.csv'); m.write_csv(p, self.rows[role], schema[role]); paths[role] = p
        summary = [{'VariationID': v, 'Guide_count':'47', 'Valid_intended_fullwindow_le3mm':'47',
                    'Zero_editable_A_offtarget_guides':'0', 'Editable_A_offtarget_ge1_guides':'47'} for v in sorted(m.MUTATIONS)]
        p = self.root/'summary.csv'; m.write_csv(p, summary, schema['feature_summary']); paths['feature_summary'] = p
        p = self.root/'manifest.json'; m.write_json(p, {'files':[{'sha256': m.digest(paths[r].read_bytes())} for r in ['features','feature_summary']]}); paths['feature_manifest'] = p
        for role, p in paths.items(): self.cfg['inputs'][role] = {'path':str(p), 'sha256':m.digest(p.read_bytes())}
        self.cpath = self.root/'config.json'; self.save()
        self.ap = self.root/'approval.json'; self.approve()

    def save(self):
        self.cpath.write_text(json.dumps(self.cfg)); self.chash = m.digest(self.cpath.read_bytes())

    def approve(self):
        self.ap.write_text(json.dumps({'approved':True, 'approved_by':'SYNTHETIC_TEST_ONLY', 'stage':'composite-ranking',
            'config_sha256':self.chash, 'run_directory':str(self.stage/self.cfg['run_id']),
            'python_sha256':m.digest(Path(m.sys.executable).read_bytes()),
            'expires_utc':(dt.datetime.now(dt.timezone.utc)+dt.timedelta(hours=1)).isoformat(),
            'robustness_approved':False, 'promotion_approved':False}))

    def execute(self, stage='composite-ranking'):
        return m.execute(stage, self.cpath, self.chash, self.ap)

    def test_verify_does_not_score_or_write(self):
        with patch.object(m, 'rank_rows', side_effect=AssertionError('Scoring forbidden')):
            self.assertEqual(self.execute('verify-inputs')['exit_code'], 0)
        self.assertEqual(list(self.stage.iterdir()), [])

    def test_success_synthetic_run_and_manifest(self):
        before = {k:Path(v['path']).read_bytes() for k,v in self.cfg['inputs'].items()}
        result = self.execute(); self.assertEqual(result['status'], 'COMPLETE')
        run = self.stage/self.cfg['run_id']; self.assertEqual({p.name for p in run.iterdir()}, m.OUTPUTS)
        for entry in json.loads((run/'output_manifest.json').read_text())['files']:
            self.assertEqual(entry['sha256'], m.digest((run/entry['path']).read_bytes()))
        self.assertEqual(before, {k:Path(v['path']).read_bytes() for k,v in self.cfg['inputs'].items()})

    def test_stale_hash(self):
        Path(self.cfg['inputs']['features']['path']).write_text('STALE')
        with self.assertRaisesRegex(ValueError, 'Stale input'): self.execute()

    def test_existing_directory(self):
        (self.stage/self.cfg['run_id']).mkdir()
        with self.assertRaisesRegex(ValueError, 'already exists'): self.execute()

    def test_wrong_config_hash(self):
        with self.assertRaisesRegex(ValueError, 'configuration hash'):
            m.execute('verify-inputs', self.cpath, '0'*64)

    def test_unassigned_path(self):
        self.cfg['inputs']['features']['path'] = None; self.save()
        with self.assertRaisesRegex(ValueError, 'Unassigned'): self.execute()

    def test_ambiguous_column_binding(self):
        self.cfg['columns']['s'] = 'Specificity_component_0_1'; self.save()
        with self.assertRaisesRegex(ValueError, 'contract'): self.execute()

    def test_no_approval(self):
        with self.assertRaisesRegex(ValueError, 'approval'):
            m.execute('composite-ranking', self.cpath, self.chash)

    def test_unapproved_record(self):
        a = json.loads(self.ap.read_text()); a['approved'] = False; self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError, 'Approval not granted'): self.execute()

    def test_approval_bound_to_config(self):
        a = json.loads(self.ap.read_text()); a['config_sha256'] = '0'*64; self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError, 'binding mismatch'): self.execute()

    def test_expired_approval(self):
        a = json.loads(self.ap.read_text()); a['expires_utc'] = '2000-01-01T00:00:00Z'; self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError, 'expired'): self.execute()

    def test_wrong_interpreter_hash(self):
        a = json.loads(self.ap.read_text()); a['python_sha256'] = '0'*64; self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError, 'interpreter hash'): self.execute()

    def test_no_unsupported_stage(self):
        for stage in ['robustness', 'promotion', 'final-library', 'blast']:
            with self.subTest(stage=stage), self.assertRaisesRegex(ValueError, 'not implemented'):
                self.execute(stage)

    def test_prohibited_write_root(self):
        self.cfg['protected_roots'] = [str(self.stage)]; self.save()
        with self.assertRaisesRegex(ValueError, 'Prohibited'): self.execute()

    def test_output_partial_fails_closed(self):
        with patch.object(m, 'rank_rows', return_value=([],[],[],{})):
            with self.assertRaisesRegex(ValueError, 'row count'): self.execute()
        result = json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text())
        self.assertEqual((result['status'], result['exit_code']), ('FAILED',4))

    def test_no_output_overwrite(self):
        p = self.root/'exists.json'; p.write_text('PRESERVE')
        with self.assertRaises(FileExistsError): m.write_json(p, {})
        self.assertEqual(p.read_text(), 'PRESERVE')

    def test_preservation_failure_blocks_complete(self):
        original = m.rank_rows
        def changed(*args):
            output = original(*args)
            # POSIX synthetic race: modify an open read-only input using another descriptor.
            # Windows share-read locking would prevent this write itself.
            if os.name != 'nt':
                Path(self.cfg['inputs']['features']['path']).write_text('CHANGED DURING SYNTHETIC RUN')
            else:
                raise ValueError('Preservation failed: synthetic injected failure')
            return output
        with patch.object(m, 'rank_rows', side_effect=changed):
            with self.assertRaisesRegex(ValueError, 'Preservation failed'): self.execute()
        result = json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text())
        self.assertEqual((result['status'], result['exit_code']), ('FAILED',4))

    def test_output_write_failure_blocks_complete(self):
        with patch.object(m, 'write_csv', side_effect=OSError('SYNTHETIC write failure')):
            with self.assertRaisesRegex(OSError, 'write failure'): self.execute()
        result = json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text())
        self.assertEqual((result['status'], result['exit_code']), ('FAILED',4))

    def test_incomplete_written_table_fails_closed(self):
        original = m.write_csv
        def truncated(path, rows, fields):
            original(path, rows[:-1] if path.name == 'ranked.csv' else rows, fields)
        with patch.object(m, 'write_csv', side_effect=truncated):
            with self.assertRaisesRegex(ValueError, 'Written output content'): self.execute()
        result = json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text())
        self.assertEqual((result['status'], result['exit_code']), ('FAILED',4))

    def test_cli_failure_status_and_exit_agree(self):
        self.cfg['inputs']['features']['path'] = None; self.save()
        with contextlib.redirect_stdout(io.StringIO()) as out:
            code = m.main(['verify-inputs','--config',str(self.cpath),'--config-sha256',self.chash])
        self.assertEqual(code,4)
        self.assertEqual(json.loads(out.getvalue())['exit_code'],code)

    def test_symlink_input_rejected(self):
        if os.name == 'nt':
            self.skipTest('Symlink privilege is not assumed on Windows')
        p = self.root/'linked.csv'; p.symlink_to(Path(self.cfg['inputs']['features']['path']))
        self.cfg['inputs']['features']['path'] = str(p); self.save()
        with self.assertRaisesRegex(ValueError, 'symlink'): self.execute()

if __name__ == '__main__':
    unittest.main(verbosity=2)
