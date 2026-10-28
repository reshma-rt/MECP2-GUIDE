#!/usr/bin/env python3
"""Proposal software: approval-gated Tier-1 assembly; no ranking/simulation code."""
import argparse
import contextlib
import datetime as dt
import json
import math
from pathlib import Path
import re
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import assembly_io as u

VERSION='MECP2_TIER1_LIBRARY_ASSEMBLY_PROPOSAL_V1'
POLICY='ROBUST_FREQUENCY_MEAN_RANK_BASELINE_SCORE_GUIDE_ID_V2'
ROLES={'adapter','io','schema','ranked','baseline_top','ranking_summary','normalization','ranking_manifest','ranking_result',
       'ranking_run_manifest','sensitivity','robust_top','robust_summary','scenarios','qc','robustness_manifest',
       'robustness_result','robustness_run_manifest'}
EXPECTED={'143546':'T1G0168','143559':'T1G0139','143597':'T1G0079','2505626':'T1G0327',
          '3237211':'T1G0271','372650':'T1G0360','4279604':'T1G0216','548706':'T1G0027'}
SCORE='Stage2E_v2_composite_score_0_100'
RANK='Stage2E_v2_rank_within_mutation'
ADDED=['Tier1_library_scope','Selection_policy','Selection_rationale','Baseline_Guide_ID',
       'Baseline_composite_score_0_100','Baseline_random_top_fraction','Baseline_mean_rank','Baseline_stability_label',
       'Robust_matches_baseline','Assembly_status','Biological_validation_status','Source_ranking_run_directory',
       'Source_robustness_run_directory','Source_ranked_sha256','Source_robust_top_sha256','Source_sensitivity_sha256',
       'Source_ranking_manifest_sha256','Source_robustness_manifest_sha256','Assembly_config_sha256','Assembly_run_id',
       'Historical_execution_provenance','Events_log_independent_verification','Production_host_enforcement_verified']
STAGES={'assemble-tier1':'PROPOSED','promotion':'NOT_AUTHORIZED','ranking':'NOT_AUTHORIZED','robustness':'NOT_AUTHORIZED'}
FILES={'tier1_candidate_library.csv','assembly_qc.json','run_manifest.json','output_manifest.json','run_result.json'}

def stability(f):return 'HIGH_STABILITY' if f>=.5 else 'MODERATE_STABILITY' if f>=.2 else 'WEIGHT_SENSITIVE'

def validate_tables(t,norm):
    """Read-only checks: no new scores, random weights, or assembled records."""
    ranked=u.cohort(t['ranked']);sens=u.cohort(t['sensitivity'])
    for role in ['baseline_top','robust_top']:
        u.indexed(t[role]);u.require(len(t[role])==8 and {r['VariationID'] for r in t[role]}==u.MUTATIONS,'Incomplete/duplicate winners: '+role)
    for role in ['ranking_summary','robust_summary']:
        u.require(len(t[role])==8 and {r['VariationID'] for r in t[role]}==u.MUTATIONS,'Incomplete/duplicate summary: '+role)
    u.require(len(t['scenarios'])==48 and len({(r['VariationID'],r['Scenario']) for r in t['scenarios']})==48,'Incomplete scenario coverage')
    scenario_names={'baseline','equal','specificity_heavy','bystander_heavy','accessibility_heavy','structure_heavy'}
    for gid,r in sens.items():
        u.require(r['Rescue_class']=='EXACT_RESTORATION','Out-of-scope rescue class')
        u.require(all(r[k]==v for k,v in ranked[gid].items()),'Ranking/sensitivity conflict: '+gid)
        for k,v in r.items():
            if any(x in k for x in ['_component_0_1','_probability','_MFE_kcal_mol','_score_0_1','_score_0_100','Sensitivity_mean_rank','Sensitivity_random_top_fraction']):
                u.number(v, k.endswith('_0_1') or k.endswith('_probability'))
        u.require(r['Stage2E_v2_status']=='HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE'
                  and r['Stage2F_v2_status']=='HARMONIZED_WEIGHT_SENSITIVITY_COMPLETE','Incomplete upstream row')
        length=u.integer(r['Guide_length_nt']);pos=u.integer(r['Target_A_position_in_guide_1based']);q=length-pos+1
        u.require(length in {21,25,31} and 6<=pos<=length-5 and q==u.integer(r['Fullwindow_catalytic_query_position_1based'])
                  ==u.integer(r['Target_A_position_in_target_segment_1based']),'Conflicting coordinate metadata')
        guide=r['Guide_sequence_RNA_5to3'];target=r['Target_segment_RNA_5to3']
        u.require(len(guide)==len(target)==length and set(guide+target)<=set('ACGU') and guide[pos-1]=='C' and target[q-1]=='A','Sequence mismatch')
        count=u.integer(r['Sensitivity_random_top_count']);frequency=u.number(r['Sensitivity_random_top_fraction'],True)
        mean=u.number(r['Sensitivity_mean_rank']);u.require(count<=10000 and frequency==count/10000 and 1<=mean<=47,'Invalid sensitivity metrics')
        u.require(abs(mean*10000-round(mean*10000))<1e-7 and r['Sensitivity_stability_label']==stability(frequency),'Invalid mean rank/label')
        u.integer(r['Sensitivity_scenario_top_count_of_6'])
    u.require(norm['cohort_size']==376,'Normalization cohort mismatch')
    for key in ['duplex_MFE_range','guide_MFE_range']:
        u.require(len(norm[key])==2 and all(math.isfinite(x) for x in norm[key]) and norm[key][0]<=norm[key][1],'Invalid normalization metadata')
    for vid in sorted(u.MUTATIONS):
        group=[r for r in sens.values() if r['VariationID']==vid]
        u.require(sorted(u.integer(r[RANK]) for r in group)==list(range(1,48)),'Incomplete baseline ranks')
        for col in ['Transcript','cDNA_change','Protein_change','Rescue_class','Target_A_transcript_position']:
            u.require(len({r[col] for r in group})==1,'Conflicting mutation assignment: '+col)
        u.require(sum(u.integer(r['Sensitivity_random_top_count']) for r in group)==10000,'Incomplete winner counts')
        u.require(sum(round(u.number(r['Sensitivity_mean_rank'])*10000) for r in group)==11280000,'Rank conservation failure')
        scenarios=[r for r in t['scenarios'] if r['VariationID']==vid]
        u.require({r['Scenario'] for r in scenarios}==scenario_names,'Incomplete mutation/scenario set')
        for sc in scenarios:
            u.require(sc['Top_Guide_ID'] in sens and sens[sc['Top_Guide_ID']]['VariationID']==vid,'Scenario candidate conflict')
            u.number(sc['Scenario_score_0_1'],True)
            for c in ['accessibility','bystander','specificity','duplex','self_structure']:u.number(sc['Weight_'+c],True)
        for r in group:
            u.require(u.integer(r['Sensitivity_scenario_top_count_of_6'])==sum(s['Top_Guide_ID']==r['Guide_ID'] for s in scenarios),'Scenario metric mismatch')
        # Verify the source-selected winner key; no simulation or reranking is performed.
        best=min(group,key=lambda r:(-u.number(r['Sensitivity_random_top_fraction']),u.number(r['Sensitivity_mean_rank']),-u.number(r[SCORE]),r['Guide_ID']))
        win=next(r for r in t['robust_top'] if r['VariationID']==vid)
        u.require(win==best,'Robust winner mismatch')
        base=next(r for r in group if r[RANK]=='1')
        u.require(next(r for r in t['baseline_top'] if r['VariationID']==vid)==ranked[base['Guide_ID']],'Baseline winner mismatch')
        rs=next(r for r in t['ranking_summary'] if r['VariationID']==vid)
        u.require(rs['Stage2E_v2_top_Guide_ID']==base['Guide_ID'] and u.number(rs['Stage2E_v2_top_score_0_100'])==u.number(base[SCORE])
                  and rs['Editable_A_offtarget_hits_le3mm']==base['Fullwindow_editable_A_offtarget_hits_le3mm'],'Ranking summary mismatch')
        sm=next(r for r in t['robust_summary'] if r['VariationID']==vid)
        expected={'VariationID':vid,'Stage2E_v2_top_Guide_ID':base['Guide_ID'],'Stage2F_v2_robust_top_Guide_ID':win['Guide_ID'],
                  'Stage2E_top_remains_robust_top':'YES' if base['Guide_ID']==win['Guide_ID'] else 'NO',
                  'Robust_top_random_top_fraction':win['Sensitivity_random_top_fraction'],'Robust_top_mean_rank':win['Sensitivity_mean_rank'],
                  'Robust_top_scenario_top_count_of_6':win['Sensitivity_scenario_top_count_of_6'],
                  'Robust_top_stability_label':win['Sensitivity_stability_label'],
                  'Robust_top_editable_A_offtarget_hits_le3mm':win['Fullwindow_editable_A_offtarget_hits_le3mm'],
                  'Robust_top_Stage2E_v2_score_0_100':win[SCORE]}
        u.require(sm==expected,'Robust summary mismatch')
    return {'status':'PASS','guides':376,'mutations':8,'proposed_entries':8,'library_created':False}

def assemble(t,cfg,chash):
    """Pure projection/join. Real invocation is behind explicit assembly approval."""
    sensitivity=u.indexed(t['sensitivity']);rows=[]
    baseline={r['VariationID']:r for r in t['baseline_top']}
    for win in sorted(t['robust_top'],key=lambda r:r['VariationID']):
        r=dict(win);b=baseline[r['VariationID']];bs=sensitivity[b['Guide_ID']]
        r.update(Tier1_library_scope='EXACT_RESTORATION_TIER1_8',Selection_policy=POLICY,
                 Selection_rationale='Maximum random win fraction; minimum mean rank; maximum baseline score; lexical Guide_ID',
                 Baseline_Guide_ID=b['Guide_ID'],Baseline_composite_score_0_100=b[SCORE],
                 Baseline_random_top_fraction=bs['Sensitivity_random_top_fraction'],Baseline_mean_rank=bs['Sensitivity_mean_rank'],
                 Baseline_stability_label=bs['Sensitivity_stability_label'],Robust_matches_baseline='YES' if b['Guide_ID']==r['Guide_ID'] else 'NO',
                 Assembly_status='STAGED_COMPUTATIONAL_SELECTION_NOT_PROMOTED',Biological_validation_status='NOT_ESTABLISHED',
                 Source_ranking_run_directory=cfg['ranking_run_directory'],Source_robustness_run_directory=cfg['robustness_run_directory'],
                 Source_ranked_sha256=cfg['inputs']['ranked']['sha256'],Source_robust_top_sha256=cfg['inputs']['robust_top']['sha256'],
                 Source_sensitivity_sha256=cfg['inputs']['sensitivity']['sha256'],
                 Source_ranking_manifest_sha256=cfg['inputs']['ranking_manifest']['sha256'],
                 Source_robustness_manifest_sha256=cfg['inputs']['robustness_manifest']['sha256'],
                 Assembly_config_sha256=chash,Assembly_run_id=cfg['run_id'],Historical_execution_provenance='UNRESOLVED',
                 Events_log_independent_verification='NOT_SUPPLIED',Production_host_enforcement_verified='false')
        rows.append(r)
    return rows

def check_manifest(manifest,roles,data,cfg,extra=None):
    entries=manifest['files'];names=[r['path'] for r in entries]
    u.require(len(names)==len(set(names)) and set(names)==set(roles.values())|set(extra or {}),'Incomplete/duplicate/unexpected manifest')
    for role,name in roles.items():
        e=next(e for e in entries if e['path']==name)
        u.require(e['sha256']==u.digest(data[role]) and e['bytes']==len(data[role]),'Manifest hash/size mismatch: '+role)
    for name,binding in (extra or {}).items():
        e=next(e for e in entries if e['path']==name)
        u.require(e['sha256']==binding['sha256'] and e['bytes']==binding['bytes'],'Uninspected artifact metadata mismatch')

def execute(stage,config_path,chash,approval_path=None):
    u.require(stage in {'verify-inputs','assemble-tier1'},'Stage not permitted')
    cp=Path(config_path).resolve(strict=True)
    with contextlib.ExitStack() as stack:
        cs=stack.enter_context(u.locked_input(cp));cb=cs.read();u.require(u.digest(cb)==chash,'Stale config hash');cfg=u.strict_json(cb)
        u.require(cfg['version']==VERSION and cfg['selection_policy']==POLICY and cfg['stages']==STAGES,'Invalid stage/policy contract')
        u.require(set(cfg['inputs'])==ROLES,'Missing/extra role')
        ph=u.digest(Path(sys.executable).read_bytes());u.require(ph==cfg['python_sha256'],'Interpreter hash mismatch')
        root=Path(cfg['staging_root']);u.require(root.is_absolute() and root.is_dir(),'Staging base missing')
        for a in [root,*root.parents]:u.require(not a.is_symlink() and not getattr(a.lstat(),'st_file_attributes',0)&0x400,'Unsafe staging path')
        u.require(bool(re.fullmatch('[A-Za-z0-9_-]+',cfg['run_id'])),'Unsafe run ID');run=root/cfg['run_id']
        u.require(not run.exists(),'Output directory already exists')
        for r in cfg['protected_roots']:u.require(not run.resolve().is_relative_to(Path(r).resolve()),'Prohibited output path')
        streams={};data={};paths={}
        for role,b in cfg['inputs'].items():
            u.require(b.get('path') and bool(re.fullmatch('[a-f0-9]{64}',b.get('sha256',''))) and isinstance(b.get('bytes'),int),'Unassigned binding: '+role)
            p=Path(b['path'])
            if not p.is_absolute():
                u.require(role in {'adapter','io','schema'},'Data path must be absolute');p=cp.parent/p
            p=p.absolute();u.require(not p.is_relative_to(run),'Input/output overlap');paths[role]=str(p)
            streams[role]=stack.enter_context(u.locked_input(p));data[role]=streams[role].read()
            u.require(len(data[role])==b['bytes'] and u.digest(data[role])==b['sha256'],'Stale/hash/size mismatch: '+role)
        u.require(len(set(paths.values()))==len(paths),'Duplicate path bindings')
        u.require(Path(paths['adapter']).resolve()==Path(__file__).resolve() and Path(paths['io']).resolve()==Path(u.__file__).resolve(),'Executing code mismatch')
        schema=u.strict_json(data['schema']);tables={role:u.table(data[role],schema[role])[0] for role in schema}
        qc=validate_tables(tables,u.strict_json(data['normalization']))
        rankroles={'ranked':'ranked.csv','baseline_top':'top_guides.csv','ranking_summary':'mutation_summary.csv',
                   'normalization':'normalization.json','ranking_run_manifest':'run_manifest.json'}
        robroles={'sensitivity':'sensitivity_all_guides.csv','robust_top':'robust_top_guides.csv','robust_summary':'mutation_summary.csv',
                  'scenarios':'scenarios.csv','qc':'simulation_qc.json','robustness_run_manifest':'run_manifest.json'}
        check_manifest(u.strict_json(data['ranking_manifest']),rankroles,data,cfg)
        check_manifest(u.strict_json(data['robustness_manifest']),robroles,data,cfg,cfg['uninspected_artifacts'])
        for prefix,stage_name in [('ranking','composite-ranking'),('robustness','robustness')]:
            r=u.strict_json(data[prefix+'_result']);rm=u.strict_json(data[prefix+'_run_manifest'])
            u.require(r['status']=='COMPLETE' and r['exit_code']==0 and r['preservation']=='PASS' and r['stage']==stage_name
                      and r['downstream_stages_executed']==[] and r['output_manifest_sha256']==u.digest(data[prefix+'_manifest']),'Upstream result incomplete')
            u.require(rm['stage']==stage_name and rm['config_sha256']==cfg[prefix+'_config_sha256'] and rm['python_sha256']==ph,'Upstream provenance mismatch')
        q=u.strict_json(data['qc']);u.require(q['guides']==376 and q['mutations']==8 and q['simulations']==10000 and q['seed']==20260914
                                             and q['fixed_cells_unchanged']==33464 and q['scenario_rows']==48 and q['winner_and_rank_sum_conservation']=='PASS','Upstream QC mismatch')
        u.require({r['VariationID']:r['Guide_ID'] for r in tables['robust_top']}==cfg['expected_winners'],'Expected winner set mismatch')
        def preserve():
            for role,s in streams.items():s.seek(0);u.require(u.digest(s.read())==cfg['inputs'][role]['sha256'],'Preservation failure: '+role)
            cs.seek(0);u.require(u.digest(cs.read())==chash,'Config preservation failure')
        if stage=='verify-inputs':
            preserve();return {**qc,'status':'VERIFIED_FOR_APPROVAL_REVIEW','stage':stage,'library_created':False,
                               'proposed_run_directory':str(run),'config_sha256':chash,'python_sha256':ph,
                               'events_log_independent_verification':'NOT_SUPPLIED','production_host_enforcement_verified':False,'exit_code':0}
        u.require(approval_path is not None,'Separate assembly approval required')
        aps=stack.enter_context(u.locked_input(Path(approval_path).resolve(strict=True)));ab=aps.read();ap=u.strict_json(ab)
        u.require(ap.get('approved') is True and bool(ap.get('approved_by')) and ap.get('stage')=='assemble-tier1'
                  and ap.get('config_sha256')==chash and ap.get('run_directory')==str(run) and ap.get('python_sha256')==ph
                  and ap.get('promotion_approved') is False,'Approval missing/wrong scope/hash/path')
        deadline=dt.datetime.fromisoformat(ap.get('expires_utc','').replace('Z','+00:00'))
        u.require(deadline.tzinfo is not None and deadline>dt.datetime.now(dt.timezone.utc),'Approval expired/invalid')
        run.mkdir(exist_ok=False)
        try:
            u.write_json(run/'run_manifest.json',{'stage':stage,'started_utc':u.utc(),'config':cfg,'config_sha256':chash,
                         'input_paths':paths,'approval_sha256':u.digest(ab),'python_sha256':ph,'python_version':sys.version,
                         'events_log_independent_verification':'NOT_SUPPLIED','production_host_enforcement_verified':False})
            rows=assemble(tables,cfg,chash);fields=schema['robust_top']+ADDED
            u.require(len(rows)==8 and len({r['Guide_ID'] for r in rows})==8 and {r['VariationID'] for r in rows}==u.MUTATIONS,'Assembly coverage/identity failure')
            for r in rows:
                source=next(x for x in tables['robust_top'] if x['VariationID']==r['VariationID'])
                u.require(all(r[c]==v for c,v in source.items()),'Source cells changed')
                u.require(r['Assembly_config_sha256']==chash and r['Source_robust_top_sha256']==cfg['inputs']['robust_top']['sha256'],'Provenance missing')
            u.write_csv(run/'tier1_candidate_library.csv',rows,fields)
            actual,_=u.table((run/'tier1_candidate_library.csv').read_bytes(),fields)
            u.require(actual==rows,'Incomplete/wrong written assembly')
            u.write_json(run/'assembly_qc.json',{'rows':8,'mutations':8,'inherited_source_cells_unchanged':760,'provenance':'PASS',
                         'events_log_independent_verification':'NOT_SUPPLIED','official_promotion':False})
            preserve();aps.seek(0);u.require(u.digest(aps.read())==u.digest(ab),'Approval preservation failure')
            u.require({p.name for p in run.iterdir()}==FILES-{'output_manifest.json','run_result.json'},'Unexpected/missing output')
            u.write_json(run/'output_manifest.json',{'files':[{'path':p.name,'bytes':p.stat().st_size,'sha256':u.digest(p.read_bytes())} for p in sorted(run.iterdir())]})
            result={'status':'COMPLETE','stage':stage,'rows':8,'mutations':8,'preservation':'PASS','exit_code':0,
                    'output_manifest_sha256':u.digest((run/'output_manifest.json').read_bytes()),'finished_utc':u.utc(),
                    'official_promotion':False,'downstream_stages_executed':[],'events_log_independent_verification':'NOT_SUPPLIED',
                    'production_host_enforcement_verified':False}
            u.write_json(run/'run_result.json',result);return result
        except BaseException as e:
            if not (run/'run_result.json').exists():u.write_json(run/'run_result.json',{'status':'FAILED_OR_INTERRUPTED','reason':str(e),'exit_code':4})
            raise

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['verify-inputs','assemble-tier1'])
    p.add_argument('--config',required=True);p.add_argument('--config-sha256',required=True);p.add_argument('--approval');a=p.parse_args(argv)
    try:print(json.dumps(execute(a.stage,a.config,a.config_sha256,a.approval),indent=2));return 0
    except (Exception,KeyboardInterrupt) as e:print(json.dumps({'status':'FAILED_CLOSED','reason':str(e),'exit_code':4},indent=2));return 4

if __name__=='__main__':sys.exit(main())
