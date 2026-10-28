#!/usr/bin/env python3
"""Approval-gated Tier-1 robustness; no ranking regeneration or promotion."""
import argparse
import contextlib
import datetime as dt
import json
import math
from pathlib import Path
import re
import sys

# Python -I omits script directory; only the explicit colocated modules are imported.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import robustness_core as core
import robustness_io as util

VERSION = 'MECP2_EXPLICIT_ROBUSTNESS_ADAPTER_V1'
COUNT = 10000
SEED = 20260914
ROLES = {'adapter', 'core', 'io', 'schema', 'ranking_source', 'robustness_source',
         'ranked', 'top', 'summary', 'normalization', 'ranking_manifest', 'ranking_result', 'ranking_run_manifest'}
STAGES = {'robustness': 'PROPOSED', 'composite-ranking': 'NOT_AUTHORIZED',
          'final-library': 'NOT_AUTHORIZED', 'promotion': 'NOT_AUTHORIZED'}
EXTRA = ['Sensitivity_random_top_count', 'Sensitivity_random_top_fraction', 'Sensitivity_mean_rank',
         'Sensitivity_scenario_top_count_of_6', 'Sensitivity_stability_label', 'Stage2F_v2_status']
SCEN_FIELDS = ['Scenario', 'VariationID', 'Top_Guide_ID', 'Scenario_score_0_1',
               'Weight_accessibility', 'Weight_bystander', 'Weight_specificity', 'Weight_duplex', 'Weight_self_structure']
MUT_FIELDS = ['VariationID', 'Stage2E_v2_top_Guide_ID', 'Stage2F_v2_robust_top_Guide_ID',
              'Stage2E_top_remains_robust_top', 'Robust_top_random_top_fraction', 'Robust_top_mean_rank',
              'Robust_top_scenario_top_count_of_6', 'Robust_top_stability_label',
              'Robust_top_editable_A_offtarget_hits_le3mm', 'Robust_top_Stage2E_v2_score_0_100']
OUTPUT_FILES = {'sensitivity_all_guides.csv', 'robust_top_guides.csv', 'scenarios.csv', 'mutation_summary.csv',
                'simulation_qc.json', 'run_manifest.json', 'events.jsonl', 'output_manifest.json', 'run_result.json'}

def label(f):
    return 'HIGH_STABILITY' if f >= .50 else 'MODERATE_STABILITY' if f >= .20 else 'WEIGHT_SENSITIVE'

def validate_baseline(rows, top, summary, normalization):
    """Check frozen records only; no new cohort ranking or simulation is generated."""
    ids = util.cohort(rows)
    util.require(len(top) == len(summary) == 8, 'Incomplete baseline winner/summary files')
    util.indexed(top)
    vids = util.MUTATIONS
    util.require({r['VariationID'] for r in top} == vids and {r['VariationID'] for r in summary} == vids,
                 'Duplicate/missing baseline mutation')
    dvals = [util.number(r['Guide_target_duplex_MFE_kcal_mol']) for r in rows]
    gvals = [util.number(r['Guide_MFE_kcal_mol']) for r in rows]
    util.require(normalization['duplex_MFE_range'] == [min(dvals), max(dvals)]
                 and normalization['guide_MFE_range'] == [min(gvals), max(gvals)]
                 and normalization['cohort_size'] == 376, 'Frozen normalization mismatch')
    util.require(normalization['constant_duplex_component'] == 0.0
                 and normalization['constant_self_component'] == 1.0, 'Constant-energy convention mismatch')
    for r in rows:
        gid = r['Guide_ID']
        util.require(r['Stage2E_v2_status'] == 'HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE'
                     and r['Stage2D_v3_status'] == 'FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE', gid+': upstream incomplete')
        length = util.integer(r['Guide_length_nt']); gp = util.integer(r['Target_A_position_in_guide_1based'])
        qp = length-gp+1
        util.require(length in {21,25,31} and 6 <= gp <= length-5 and
                     util.integer(r['Target_A_position_in_target_segment_1based']) == qp ==
                     util.integer(r['Fullwindow_catalytic_query_position_1based']), gid+': coordinate conflict')
        guide, target = r['Guide_sequence_RNA_5to3'], r['Target_segment_RNA_5to3']
        util.require(len(guide) == len(target) == length and set(guide+target) <= set('ACGU')
                     and guide[gp-1] == 'C' and target[qp-1] == 'A', gid+': sequence/coordinate mismatch')
        util.require(r['Transcript'] and util.integer(r['Target_A_transcript_position']) > 0, gid+': missing target metadata')
        a = (util.number(r['Target_A_unpaired_probability'],True)+util.number(r['Target_A_5nt_mean_unpaired_probability'],True))/2
        d = util.number(r['Guide_target_duplex_MFE_kcal_mol']); g = util.number(r['Guide_MFE_kcal_mol'])
        dn = 0.0 if min(dvals) == max(dvals) else 1-(d-min(dvals))/(max(dvals)-min(dvals))
        gn = 1.0 if min(gvals) == max(gvals) else (g-min(gvals))/(max(gvals)-min(gvals))
        expected = [a,util.number(r['Bystander_precision_score_0_1'],True),
                    util.number(r['Fullwindow_GuideAware_transcriptome_specificity_score_0_1'],True),dn,gn]
        comps = [util.number(r[c],True) for c in core.COMPONENTS]
        util.require(all(abs(x-y) <= 1e-12 for x,y in zip(expected,comps)), gid+': baseline component mismatch')
        weights = list(core.SCENARIOS['baseline'])
        for name,w in zip(['accessibility','bystander','specificity','duplex','self_structure'],weights):
            util.require(util.number(r['Stage2E_v2_weight_'+name]) == w, 'Baseline weight mismatch')
        score = util.number(r[core.BASELINE_SCORE])
        util.require(0 <= score <= 100 and abs(score-100*sum(w*c for w,c in zip(weights,comps))) <= 1e-10,
                     gid+': baseline score mismatch')
        util.integer(r['Fullwindow_editable_A_offtarget_hits_le3mm'])
    for vid in sorted(vids):
        # Verify supplied ordinal order via adjacent sort keys; do not regenerate ranking outputs.
        group = [r for r in rows if r['VariationID'] == vid]
        for column in ['Transcript','cDNA_change','Protein_change','Target_A_transcript_position','Rescue_class']:
            util.require(len({r[column] for r in group}) == 1, 'Mutation/target metadata conflict: '+column)
        util.require([util.integer(r[core.BASELINE_RANK]) for r in group] == list(range(1,48)), 'Baseline ranks incomplete/order conflict')
        keys = [(-util.number(r[core.BASELINE_SCORE]),util.integer(r['Fullwindow_editable_A_offtarget_hits_le3mm']),r['Guide_ID']) for r in group]
        util.require(all(x <= y for x,y in zip(keys,keys[1:])), 'Baseline ordering/tie conflict')
        t = next(r for r in top if r['VariationID'] == vid)
        util.require(t == group[0], 'Baseline top record mismatch')
        sm = next(r for r in summary if r['VariationID'] == vid)
        util.require(sm['Stage2E_v2_top_Guide_ID'] == t['Guide_ID'] and
                     util.number(sm['Stage2E_v2_top_score_0_100']) == util.number(t[core.BASELINE_SCORE]) and
                     sm['Editable_A_offtarget_hits_le3mm'] == t['Fullwindow_editable_A_offtarget_hits_le3mm'], 'Baseline summary mismatch')
    return ids

def reconcile(original, outputs, n_sim=COUNT):
    """Independent output invariants; no random draw or second simulation."""
    all_rows, winners, scenarios, summaries, counts, sums = outputs
    ids = util.indexed(original); out = util.indexed(all_rows)
    util.require(set(out) == set(ids) and len(winners) == len(summaries) == 8 and len(scenarios) == 48, 'Incomplete simulation outputs')
    for gid,r in out.items():
        util.require(all(r[k] == v for k,v in ids[gid].items()), 'Frozen ranking cells changed')
        count = r['Sensitivity_random_top_count']; frac = r['Sensitivity_random_top_fraction']; mean = r['Sensitivity_mean_rank']
        util.require(isinstance(count,int) and 0 <= count <= n_sim and counts.get(gid,0) == count
                     and frac == count/n_sim, 'Winner count/frequency denominator mismatch')
        util.require(mean == sums[gid]/n_sim and 1 <= mean <= 47, 'Mean-rank mismatch')
        util.require(r['Sensitivity_stability_label'] == label(frac) and
                     r['Stage2F_v2_status'] == 'HARMONIZED_WEIGHT_SENSITIVITY_COMPLETE', 'Stability/status mismatch')
    util.require(len({(r['Scenario'],r['VariationID']) for r in scenarios}) == 48, 'Duplicate/missing scenario')
    for vid in sorted(util.MUTATIONS):
        group = [r for r in all_rows if r['VariationID'] == vid]
        util.require(len(group) == 47 and sum(r['Sensitivity_random_top_count'] for r in group) == n_sim,
                     'Mutation winner-count conservation failed')
        util.require(sum(sums[r['Guide_ID']] for r in group) == n_sim*47*48//2, 'Mutation rank-sum conservation failed')
        group_scen = [r for r in scenarios if r['VariationID'] == vid]
        util.require({r['Scenario'] for r in group_scen} == set(core.SCENARIOS), 'Scenario coverage failed')
        for sc in group_scen:
            w = core.SCENARIOS[sc['Scenario']]
            scored = [(sum(wi*util.number(r[c]) for wi,c in zip(w,core.COMPONENTS)),r['Guide_ID']) for r in group]
            best = min(scored,key=lambda pair:(-pair[0],pair[1]))
            util.require(sc['Top_Guide_ID'] == best[1] and sc['Scenario_score_0_1'] == best[0], 'Scenario tie/score mismatch')
            util.require([sc['Weight_'+c] for c in ['accessibility','bystander','specificity','duplex','self_structure']] == w,
                         'Scenario weights mismatch')
        for r in group:
            util.require(r['Sensitivity_scenario_top_count_of_6'] == sum(sc['Top_Guide_ID'] == r['Guide_ID'] for sc in group_scen),
                         'Scenario counts mismatch')
        best = min(group,key=lambda r:(-r['Sensitivity_random_top_fraction'],r['Sensitivity_mean_rank'],
                                      -util.number(r[core.BASELINE_SCORE]),r['Guide_ID']))
        win = [r for r in winners if r['VariationID'] == vid]
        util.require(len(win) == 1 and win[0] == best, 'Robust selection mismatch')
        sm = [r for r in summaries if r['VariationID'] == vid]
        baseline = next(r for r in group if util.integer(r[core.BASELINE_RANK]) == 1)
        expected = {'VariationID':vid,'Stage2E_v2_top_Guide_ID':baseline['Guide_ID'],
                    'Stage2F_v2_robust_top_Guide_ID':best['Guide_ID'],
                    'Stage2E_top_remains_robust_top':'YES' if baseline['Guide_ID']==best['Guide_ID'] else 'NO',
                    'Robust_top_random_top_fraction':best['Sensitivity_random_top_fraction'],
                    'Robust_top_mean_rank':best['Sensitivity_mean_rank'],
                    'Robust_top_scenario_top_count_of_6':best['Sensitivity_scenario_top_count_of_6'],
                    'Robust_top_stability_label':best['Sensitivity_stability_label'],
                    'Robust_top_editable_A_offtarget_hits_le3mm':best['Fullwindow_editable_A_offtarget_hits_le3mm'],
                    'Robust_top_Stage2E_v2_score_0_100':best[core.BASELINE_SCORE]}
        util.require(sm == [expected], 'Mutation summary mismatch')
    return {'guides':376,'mutations':8,'simulations':n_sim,'seed':SEED,'scenario_rows':48,
            'winner_and_rank_sum_conservation':'PASS','fixed_cells_unchanged':376*89,
            'note':'Computational weight sensitivity; not experimental validation'}

def approve(record, chash, run, python_hash):
    util.require(record.get('approved') is True and bool(record.get('approved_by')), 'Approval not granted/named')
    util.require(record.get('stage') == 'robustness' and record.get('config_sha256') == chash
                 and record.get('run_directory') == str(run) and record.get('python_sha256') == python_hash,
                 'Approval stage/hash/path mismatch')
    util.require(record.get('simulations') == COUNT and record.get('seed') == SEED
                 and record.get('final_library_approved') is False and record.get('promotion_approved') is False,
                 'Approval scope mismatch')
    expiry = dt.datetime.fromisoformat(record.get('expires_utc','').replace('Z','+00:00'))
    util.require(expiry.tzinfo is not None and expiry > dt.datetime.now(dt.timezone.utc), 'Approval expired/invalid')

def execute(stage, config_path, config_hash, approval_path=None):
    util.require(stage in {'verify-inputs','robustness'}, 'Stage prohibited/not implemented')
    config_path = Path(config_path).resolve(strict=True)
    with contextlib.ExitStack() as stack:
        cs = stack.enter_context(util.locked_input(config_path)); cb = cs.read()
        util.require(util.digest(cb) == config_hash, 'Stale configuration hash')
        cfg = util.strict_json(cb)
        util.require(cfg['version'] == VERSION and cfg['simulations'] == COUNT and cfg['seed'] == SEED
                     and cfg['components'] == core.COMPONENTS and cfg['scenarios'] == core.SCENARIOS
                     and cfg['stages'] == STAGES, 'Method/schema/stage contract mismatch')
        util.require(set(cfg['inputs']) == ROLES, 'Missing/extra input role')
        ph = util.digest(Path(sys.executable).read_bytes())
        util.require(ph == cfg['python_sha256'], 'Interpreter hash mismatch')
        root = Path(cfg['staging_root'])
        util.require(root.is_absolute() and root.is_dir(), 'Staging root unavailable')
        for a in [root,*root.parents]:
            util.require(not a.is_symlink() and not getattr(a.lstat(),'st_file_attributes',0)&0x400,'Unsafe staging path')
        util.require(bool(re.fullmatch('[A-Za-z0-9_-]+',cfg['run_id'])), 'Unsafe run ID')
        run = root/cfg['run_id']
        util.require(not run.exists(), 'Run directory already exists')
        for protected in cfg['protected_roots']:
            util.require(not run.resolve().is_relative_to(Path(protected).resolve()), 'Prohibited output root')
        streams, data, paths = {}, {}, {}
        for role,binding in cfg['inputs'].items():
            util.require(binding.get('path') and bool(re.fullmatch('[a-f0-9]{64}',binding.get('sha256','')))
                         and isinstance(binding.get('bytes'),int), 'Unassigned input: '+role)
            path = Path(binding['path'])
            if not path.is_absolute():
                util.require(role in {'adapter','core','io','schema','ranking_source','robustness_source'},'Data path must be absolute')
                path = config_path.parent/path
            path = path.absolute()
            util.require(not path.is_relative_to(run), 'Input overlaps outputs')
            streams[role] = stack.enter_context(util.locked_input(path)); data[role] = streams[role].read()
            util.require(len(data[role]) == binding['bytes'] and util.digest(data[role]) == binding['sha256'], 'Stale/size/hash mismatch: '+role)
            paths[role] = str(path)
        util.require(len(set(paths.values())) == len(paths), 'Duplicate input path')
        for role,module in [('adapter',sys.modules[__name__]),('core',core),('io',util)]:
            util.require(Path(paths[role]).resolve() == Path(module.__file__).resolve(), 'Imported code path mismatch')
        schema = util.strict_json(data['schema'])
        rows,fields = util.table(data['ranked'],schema['ranked'])
        top,_ = util.table(data['top'],schema['top']); summary,_ = util.table(data['summary'],schema['summary'])
        validate_baseline(rows,top,summary,util.strict_json(data['normalization']))
        manifest = util.strict_json(data['ranking_manifest']); result = util.strict_json(data['ranking_result'])
        rm = util.strict_json(data['ranking_run_manifest'])
        for role in ['ranked','top','summary','normalization','ranking_run_manifest']:
            name = {'top':'top_guides.csv','summary':'mutation_summary.csv','ranking_run_manifest':'run_manifest.json'}.get(role,role+'.csv' if role=='ranked' else role+'.json')
            hits = [e for e in manifest['files'] if e['path']==name]
            util.require(len(hits)==1 and hits[0]['sha256']==cfg['inputs'][role]['sha256']
                         and hits[0]['bytes']==cfg['inputs'][role]['bytes'], 'Ranking manifest association mismatch')
        util.require(result['status']=='COMPLETE' and result['exit_code']==0 and result['preservation']=='PASS'
                     and result['stage']=='composite-ranking' and result['downstream_stages_executed']==[]
                     and result['output_manifest_sha256']==util.digest(data['ranking_manifest']), 'Upstream ranking not complete')
        util.require(rm['stage']=='composite-ranking' and rm['config_sha256']==cfg['ranking_config_sha256']
                     and rm['python_sha256']==cfg['python_sha256'], 'Upstream run provenance mismatch')
        def preserve():
            for role,s in streams.items():
                s.seek(0)
                util.require(util.digest(s.read())==cfg['inputs'][role]['sha256'], 'Preservation failure: '+role)
            cs.seek(0); util.require(util.digest(cs.read())==config_hash,'Config preservation failure')
        if stage=='verify-inputs':
            preserve()
            return {'status':'VERIFIED_FOR_APPROVAL_REVIEW','stage':stage,'guides':376,'mutations':8,
                    'simulations_executed':0,'config_sha256':config_hash,'python_sha256':ph,'python_version':sys.version,
                    'proposed_run_directory':str(run),'production_host_enforcement_verified':False,'exit_code':0}
        util.require(approval_path is not None,'Separate robustness approval required')
        ap = Path(approval_path).resolve(strict=True)
        aps = stack.enter_context(util.locked_input(ap)); ab = aps.read()
        approve(util.strict_json(ab),config_hash,run,ph)
        run.mkdir(exist_ok=False)
        try:
            util.write_json(run/'run_manifest.json',{'stage':stage,'version':VERSION,'started_utc':util.utc(),
                            'config_sha256':config_hash,'config':cfg,'input_paths':paths,'approval_sha256':util.digest(ab),
                            'python_sha256':ph,'python_version':sys.version,'historical_provenance':'UNRESOLVED',
                            'production_host_enforcement_verified':False})
            with (run/'events.jsonl').open('x',encoding='utf-8') as log:
                log.write(json.dumps({'event':'START','utc':util.utc(),'simulations':COUNT,'seed':SEED})+'\n');log.flush()
                outputs = core.calculate(rows,n_sim=COUNT,seed=SEED)
                qc = reconcile(rows,outputs)
                log.write(json.dumps({'event':'SIMULATION_AND_QC_FINISHED','utc':util.utc()})+'\n')
            for name,records,headers in [('sensitivity_all_guides.csv',outputs[0],fields+EXTRA),
                                         ('robust_top_guides.csv',outputs[1],fields+EXTRA),
                                         ('scenarios.csv',outputs[2],SCEN_FIELDS),('mutation_summary.csv',outputs[3],MUT_FIELDS)]:
                util.write_csv(run/name,records,headers)
                written,wh = util.table((run/name).read_bytes(),headers)
                util.require(written==[{k:str(v) for k,v in r.items()} for r in records],'Incomplete/wrong written output: '+name)
            util.write_json(run/'simulation_qc.json',qc)
            preserve(); aps.seek(0); util.require(util.digest(aps.read())==util.digest(ab),'Approval preservation failure')
            util.require({p.name for p in run.iterdir()}==OUTPUT_FILES-{'output_manifest.json','run_result.json'},'Unexpected/missing output file')
            entries = [{'path':p.name,'bytes':p.stat().st_size,'sha256':util.digest(p.read_bytes())} for p in sorted(run.iterdir())]
            util.write_json(run/'output_manifest.json',{'files':entries,'excludes':['output_manifest.json','run_result.json']})
            status = {'status':'COMPLETE','stage':stage,'rows':376,'robust_top_rows':8,'scenario_rows':48,'summary_rows':8,
                      'simulations':COUNT,'seed':SEED,'preservation':'PASS','exit_code':0,'finished_utc':util.utc(),
                      'output_manifest_sha256':util.digest((run/'output_manifest.json').read_bytes()),
                      'downstream_stages_executed':[],'production_host_enforcement_verified':False}
            util.write_json(run/'run_result.json',status)
            return status
        except BaseException as e:
            if not (run/'run_result.json').exists():
                util.write_json(run/'run_result.json',{'status':'FAILED_OR_INTERRUPTED','stage':stage,'exit_code':4,
                                'reason':str(e),'partial_outputs_not_approved':True,'finished_utc':util.utc()})
            raise

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['verify-inputs','robustness']);p.add_argument('--config',required=True)
    p.add_argument('--config-sha256',required=True);p.add_argument('--approval')
    args=p.parse_args(argv)
    try:
        result=execute(args.stage,args.config,args.config_sha256,args.approval)
        print(json.dumps(result,indent=2));return 0
    except (Exception,KeyboardInterrupt) as e:
        print(json.dumps({'status':'FAILED_CLOSED','stage':args.stage,'reason':str(e),'exit_code':4},indent=2));return 4

if __name__=='__main__':
    sys.exit(main())
