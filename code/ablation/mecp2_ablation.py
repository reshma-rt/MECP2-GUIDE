#!/usr/bin/env python3
"""Supplementary comparison only; never replaces frozen rankings or candidates."""
import argparse
import contextlib
import datetime as dt
import json
from pathlib import Path
import re
import statistics
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import ablation_io as u

VERSION='MECP2_COMPONENT_ABLATION_V1'
NAMES=['accessibility','bystander','specificity','duplex','self_structure']
COLS=['Stage2E_v2_'+c+'_component_0_1' for c in NAMES]
BASE=[.25,.30,.30,.10,.05]
SCORE='Stage2E_v2_composite_score_0_100'
RANK='Stage2E_v2_rank_within_mutation'
COUNT='Fullwindow_editable_A_offtarget_hits_le3mm'
ROLES={'adapter','io','schema','ranked','robust_top','ranking_manifest','robustness_manifest','ranking_result','robustness_result'}
OUT={'comparison_rankings.csv','comparison_winners.csv','mutation_comparisons.csv','scenario_summary.csv',
     'scenario_definitions.json','run_manifest.json','output_manifest.json','run_result.json'}
R_FIELDS=['Scenario','VariationID','Guide_ID','Guide_length_nt','Comparison_score_0_100','Comparison_rank',
          'Frozen_baseline_score_0_100','Frozen_baseline_rank','Primary_score_tie_count','Is_frozen_baseline_winner','Is_frozen_robust_winner']+COLS
W_FIELDS=R_FIELDS+['Frozen_baseline_winner','Frozen_robust_winner','Matches_baseline','Matches_robust']
M_FIELDS=['Scenario','VariationID','Winner','Matches_baseline','Matches_robust','Winner_length_nt',
          'Ordinal_rank_spearman_vs_baseline','Top5_jaccard_vs_baseline','Primary_top_score_tie_count']
S_FIELDS=['Scenario','Mutation_count','Baseline_winner_changes','Robust_winner_changes','Winners_21nt','Winners_25nt','Winners_31nt',
          'Median_ordinal_rank_spearman','Mean_top5_jaccard']

def definitions():
    result={'baseline_control':list(BASE)}
    for i in range(1,5):
        w=[0.]*5;w[i]=1.;result['only_'+NAMES[i]]=w
    for i in range(1,5):
        result['without_'+NAMES[i]]=[0. if j==i else x/(1-BASE[i]) for j,x in enumerate(BASE)]
    return result

def left_score(row,w):
    value=0.
    for wi,col in zip(w,COLS):value+=wi*u.number(row[col],True)
    return value*100.

def validate(rows,robust):
    indexed=u.cohort(rows);u.indexed(robust)
    u.require(len(robust)==8 and {r['VariationID'] for r in robust}==u.MUTATIONS,'Incomplete robust selections')
    for r in robust:u.require(r['Guide_ID'] in indexed and all(r[k]==v for k,v in indexed[r['Guide_ID']].items()),'Robust/ranked lineage mismatch')
    for row in rows:
        for col in COLS:u.number(row[col],True)
        u.require(row['Rescue_class']=='EXACT_RESTORATION' and row['Stage2E_v2_status']=='HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE','Wrong cohort/upstream status')
        u.require(abs(left_score(row,BASE)-u.number(row[SCORE]))<1e-10,'Stored baseline score mismatch')
        u.integer(row[COUNT]);u.integer(row[RANK])
        length=u.integer(row['Guide_length_nt']);gp=u.integer(row['Target_A_position_in_guide_1based'])
        q=length-gp+1
        u.require(length in {21,25,31} and 6<=gp<=length-5 and q==u.integer(row['Target_A_position_in_target_segment_1based'])
                  ==u.integer(row['Fullwindow_catalytic_query_position_1based']),'Coordinate conflict')
        for name,w in zip(NAMES,BASE):u.require(u.number(row['Stage2E_v2_weight_'+name])==w,'Frozen weights mismatch')
    for vid in sorted(u.MUTATIONS):
        group=[r for r in rows if r['VariationID']==vid]
        u.require(sorted(u.integer(r[RANK]) for r in group)==list(range(1,48)),'Incomplete baseline ranks')
        ordered=sorted(group,key=lambda r:u.integer(r[RANK]))
        keys=[(-u.number(r[SCORE]),u.integer(r[COUNT]),r['Guide_ID']) for r in ordered]
        u.require(all(a<=b for a,b in zip(keys,keys[1:])),'Frozen baseline tie/order mismatch')
        u.require(len({r[COLS[0]] for r in group})==1,'Accessibility no longer constant: revise this analysis specification')
    return {'guides':376,'mutations':8,'scenario_count':9,'comparison_calculated':False}

def compare(rows,robust):
    """Call on real data ONLY after separate component-ablation approval."""
    robust_map={r['VariationID']:r['Guide_ID'] for r in robust}
    ranked_output=[];winner_output=[];comparisons=[];summaries=[]
    for scenario,w in definitions().items():
        for vid in sorted(u.MUTATIONS):
            group=[r for r in rows if r['VariationID']==vid]
            baseline=min(group,key=lambda r:u.integer(r[RANK]))['Guide_ID']
            scored=[(left_score(r,w),r) for r in group]
            scored.sort(key=lambda pair:(-pair[0],u.integer(pair[1][COUNT]),pair[1]['Guide_ID']))
            values=[x[0] for x in scored];records=[]
            for rank,(score,row) in enumerate(scored,1):
                record={'Scenario':scenario,'VariationID':vid,'Guide_ID':row['Guide_ID'],'Guide_length_nt':row['Guide_length_nt'],
                        'Comparison_score_0_100':score,'Comparison_rank':rank,'Frozen_baseline_score_0_100':row[SCORE],
                        'Frozen_baseline_rank':row[RANK],'Primary_score_tie_count':values.count(score),
                        'Is_frozen_baseline_winner':'YES' if row['Guide_ID']==baseline else 'NO',
                        'Is_frozen_robust_winner':'YES' if row['Guide_ID']==robust_map[vid] else 'NO'}
                record.update({c:row[c] for c in COLS});records.append(record)
            winner=dict(records[0]);winner.update(Frozen_baseline_winner=baseline,Frozen_robust_winner=robust_map[vid],
                          Matches_baseline=records[0]['Is_frozen_baseline_winner'],Matches_robust=records[0]['Is_frozen_robust_winner'])
            n=len(records);d2=sum((r['Comparison_rank']-u.integer(r['Frozen_baseline_rank']))**2 for r in records)
            rho=1.-6*d2/(n*(n*n-1))
            oldtop={r['Guide_ID'] for r in records if u.integer(r['Frozen_baseline_rank'])<=5};newtop={r['Guide_ID'] for r in records[:5]}
            jaccard=len(oldtop&newtop)/len(oldtop|newtop)
            comparisons.append({'Scenario':scenario,'VariationID':vid,'Winner':winner['Guide_ID'],
                                'Matches_baseline':winner['Matches_baseline'],'Matches_robust':winner['Matches_robust'],
                                'Winner_length_nt':winner['Guide_length_nt'],'Ordinal_rank_spearman_vs_baseline':rho,
                                'Top5_jaccard_vs_baseline':jaccard,'Primary_top_score_tie_count':winner['Primary_score_tie_count']})
            if scenario=='baseline_control':
                u.require(all(u.integer(r['Frozen_baseline_rank'])==r['Comparison_rank'] for r in records),'Baseline control rank mismatch')
            ranked_output.extend(records);winner_output.append(winner)
        sub=[r for r in comparisons if r['Scenario']==scenario]
        summaries.append({'Scenario':scenario,'Mutation_count':8,'Baseline_winner_changes':sum(r['Matches_baseline']=='NO' for r in sub),
                          'Robust_winner_changes':sum(r['Matches_robust']=='NO' for r in sub),
                          **{'Winners_'+str(n)+'nt':sum(int(r['Winner_length_nt'])==n for r in sub) for n in [21,25,31]},
                          'Median_ordinal_rank_spearman':statistics.median(r['Ordinal_rank_spearman_vs_baseline'] for r in sub),
                          'Mean_top5_jaccard':sum(r['Top5_jaccard_vs_baseline'] for r in sub)/8})
    return ranked_output,winner_output,comparisons,summaries

def check_outputs(outputs):
    ranked,winners,mutations,summaries=outputs
    u.require([len(x) for x in outputs]==[3384,72,72,9],'Incomplete comparison outputs')
    u.require(len({(r['Scenario'],r['VariationID'],r['Guide_ID']) for r in ranked})==3384,'Duplicate comparison rows')
    u.require(len({(r['Scenario'],r['VariationID']) for r in winners})==72,'Duplicate/missing winner')
    for scenario in definitions():
        u.require(sum(r['Scenario']==scenario for r in ranked)==376,'Scenario incomplete')
        for vid in u.MUTATIONS:
            group=[r for r in ranked if r['Scenario']==scenario and r['VariationID']==vid]
            u.require(sorted(r['Comparison_rank'] for r in group)==list(range(1,48)),'Ranks incomplete')
            u.require(all(0<=r['Comparison_score_0_100']<=100.000000000001 for r in group),'Invalid score')

def execute(stage,cp,chash,ap=None):
    u.require(stage in {'verify-inputs','component-ablation'},'Stage prohibited')
    cp=Path(cp).resolve(strict=True)
    with contextlib.ExitStack() as stack:
        cs=stack.enter_context(u.locked_input(cp));cb=cs.read();u.require(u.digest(cb)==chash,'Stale config hash');cfg=u.strict_json(cb)
        u.require(cfg['version']==VERSION and cfg['scenarios']==definitions() and cfg['tie_rule']=='score_desc_editable_count_asc_Guide_ID_asc'
                  and cfg['promotion_approved'] is False,'Method contract mismatch')
        u.require(set(cfg['inputs'])==ROLES,'Missing/extra binding')
        ph=u.digest(Path(sys.executable).read_bytes());u.require(ph==cfg['python_sha256'],'Interpreter mismatch')
        root=Path(cfg['staging_root']);u.require(root.is_absolute() and root.is_dir(),'Staging unavailable')
        for parent in [root,*root.parents]:u.require(not parent.is_symlink() and not getattr(parent.lstat(),'st_file_attributes',0)&0x400,'Unsafe staging root')
        u.require(bool(re.fullmatch('[A-Za-z0-9_-]+',cfg['run_id'])),'Invalid run ID');run=root/cfg['run_id']
        u.require(not run.exists(),'Run directory already exists')
        for protected in cfg['protected_roots']:u.require(not run.resolve().is_relative_to(Path(protected).resolve()),'Prohibited write root')
        streams={};data={};paths={}
        for role,b in cfg['inputs'].items():
            u.require(b.get('path') and bool(re.fullmatch('[a-f0-9]{64}',b.get('sha256',''))) and isinstance(b.get('bytes'),int),'Unassigned binding')
            p=Path(b['path'])
            if not p.is_absolute():
                u.require(role in {'adapter','io','schema'},'Data path must be absolute');p=cp.parent/p
            p=p.absolute();u.require(not p.is_relative_to(run),'Input/output overlap');paths[role]=str(p)
            streams[role]=stack.enter_context(u.locked_input(p));data[role]=streams[role].read()
            u.require(len(data[role])==b['bytes'] and u.digest(data[role])==b['sha256'],'Stale input/hash/size: '+role)
        u.require(len(set(paths.values()))==len(paths),'Duplicate input bindings')
        u.require(Path(paths['adapter']).resolve()==Path(__file__).resolve() and Path(paths['io']).resolve()==Path(u.__file__).resolve(),'Code path mismatch')
        schema=u.strict_json(data['schema']);rows,fields=u.table(data['ranked'],schema['ranked']);robust,_=u.table(data['robust_top'],schema['robust_top'])
        qc=validate(rows,robust)
        for prefix,role,name,expectedstage in [('ranking','ranked','ranked.csv','composite-ranking'),('robustness','robust_top','robust_top_guides.csv','robustness')]:
            manifest=u.strict_json(data[prefix+'_manifest']);names=[e['path'] for e in manifest['files']]
            u.require(len(names)==len(set(names)),'Duplicate upstream manifest entry')
            hit=[e for e in manifest['files'] if e['path']==name]
            u.require(len(hit)==1 and hit[0]['sha256']==cfg['inputs'][role]['sha256'] and hit[0]['bytes']==len(data[role]),'Upstream manifest mismatch')
            result=u.strict_json(data[prefix+'_result'])
            u.require(result['stage']==expectedstage and result['status']=='COMPLETE' and result['exit_code']==0 and result['preservation']=='PASS'
                      and result['output_manifest_sha256']==u.digest(data[prefix+'_manifest']),'Upstream result failed')
        def preserve():
            for role,s in streams.items():s.seek(0);u.require(u.digest(s.read())==cfg['inputs'][role]['sha256'],'Preservation failure')
            cs.seek(0);u.require(u.digest(cs.read())==chash,'Config preservation failure')
        if stage=='verify-inputs':
            preserve();return {**qc,'status':'VERIFIED_FOR_APPROVAL_REVIEW','config_sha256':chash,'python_sha256':ph,
                               'proposed_run_directory':str(run),'files_written':0,'production_host_enforcement_verified':False,'exit_code':0}
        u.require(ap is not None,'Separate component-ablation approval required')
        aps=stack.enter_context(u.locked_input(Path(ap).resolve(strict=True)));ab=aps.read();approval=u.strict_json(ab)
        u.require(approval.get('approved') is True and bool(approval.get('approved_by')) and approval.get('stage')==stage
                  and approval.get('config_sha256')==chash and approval.get('python_sha256')==ph and approval.get('run_directory')==str(run)
                  and approval.get('promotion_approved') is False,'Approval scope/hash/path mismatch')
        expiry=dt.datetime.fromisoformat(approval.get('expires_utc','').replace('Z','+00:00'))
        u.require(expiry.tzinfo is not None and expiry>dt.datetime.now(dt.timezone.utc),'Approval expired/invalid')
        run.mkdir(exist_ok=False)
        try:
            u.write_json(run/'run_manifest.json',{'stage':stage,'config':cfg,'config_sha256':chash,'input_paths':paths,'approval_sha256':u.digest(ab),
                         'python_sha256':ph,'python_version':sys.version,'started_utc':u.utc(),'production_host_enforcement_verified':False})
            outputs=compare(rows,robust);check_outputs(outputs)
            for name,records,headers in zip(['comparison_rankings.csv','comparison_winners.csv','mutation_comparisons.csv','scenario_summary.csv'],outputs,[R_FIELDS,W_FIELDS,M_FIELDS,S_FIELDS]):
                u.write_csv(run/name,records,headers);actual,_=u.table((run/name).read_bytes(),headers)
                u.require(actual==[{k:str(v) for k,v in r.items()} for r in records],'Incomplete/wrong written output')
            u.write_json(run/'scenario_definitions.json',{'scenarios':definitions(),'components':COLS,'tie_rule':cfg['tie_rule'],
                         'interpretation':'Decision sensitivity only; no experimental superiority or editing probability.'})
            preserve();aps.seek(0);u.require(u.digest(aps.read())==u.digest(ab),'Approval preservation failure')
            u.require({p.name for p in run.iterdir()}==OUT-{'output_manifest.json','run_result.json'},'Unexpected/missing output')
            u.write_json(run/'output_manifest.json',{'files':[{'path':p.name,'bytes':p.stat().st_size,'sha256':u.digest(p.read_bytes())} for p in sorted(run.iterdir())]})
            result={'status':'COMPLETE','stage':stage,'scenarios':9,'ranking_rows':3384,'winner_rows':72,'mutation_comparison_rows':72,
                    'preservation':'PASS','exit_code':0,'output_manifest_sha256':u.digest((run/'output_manifest.json').read_bytes()),
                    'finished_utc':u.utc(),'official_outputs_modified':False,'downstream_stages_executed':[],
                    'production_host_enforcement_verified':False}
            u.write_json(run/'run_result.json',result);return result
        except BaseException as e:
            if not (run/'run_result.json').exists():u.write_json(run/'run_result.json',{'status':'FAILED_OR_INTERRUPTED','reason':str(e),'exit_code':4})
            raise

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['verify-inputs','component-ablation'])
    p.add_argument('--config',required=True);p.add_argument('--config-sha256',required=True);p.add_argument('--approval');a=p.parse_args(argv)
    try:print(json.dumps(execute(a.stage,a.config,a.config_sha256,a.approval),indent=2));return 0
    except (Exception,KeyboardInterrupt) as e:print(json.dumps({'status':'FAILED_CLOSED','reason':str(e),'exit_code':4},indent=2));return 4

if __name__=='__main__':sys.exit(main())
