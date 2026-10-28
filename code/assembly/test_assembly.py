"""Synthetic projection and failure tests; no real library is assembled."""
import contextlib
import copy
import datetime as dt
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import mecp2_tier1_assembly as m

HERE=Path(__file__).resolve().parent

def fixtures():
    schema=json.loads((HERE/'input_schema.json').read_text());ranked=[];sens=[];top=[];wins=[];rs=[];sm=[];sc=[]
    for vi,vid in enumerate(sorted(m.u.MUTATIONS)):
        for i in range(47):
            r={k:'' for k in schema['ranked']};gid='T1G%04d'%(vi*47+i+1)
            for k in r:
                if any(x in k for x in ['_component_0_1','_probability','_MFE_kcal_mol','_score_0_1','_score_0_100']):r[k]='0.5'
            r.update(Guide_ID=gid,VariationID=vid,Rescue_class='EXACT_RESTORATION',Transcript='SYNTHETIC_ONLY',Guide_length_nt='21',Target_A_transcript_position='1000',
                     Target_A_position_in_guide_1based='6',Target_A_position_in_target_segment_1based='16',
                     Fullwindow_catalytic_query_position_1based='16',Guide_sequence_RNA_5to3='UUUUUCUUUUUUUUUUUUUUU',
                     Target_segment_RNA_5to3='A'*21,Stage2E_v2_status='HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE',
                     Fullwindow_editable_A_offtarget_hits_le3mm='0',Stage2E_v2_composite_score_0_100=str(100-i),
                     Stage2E_v2_rank_within_mutation=str(i+1))
            ranked.append(r);s=dict(r)
            s.update(Sensitivity_random_top_count='10000' if i==1 else '0',Sensitivity_random_top_fraction='1.0' if i==1 else '0.0',
                     Sensitivity_mean_rank=str(1 if i==1 else 2 if i==0 else i+1),Sensitivity_scenario_top_count_of_6='6' if i==1 else '0',
                     Sensitivity_stability_label='HIGH_STABILITY' if i==1 else 'WEIGHT_SENSITIVE',Stage2F_v2_status='HARMONIZED_WEIGHT_SENSITIVITY_COMPLETE')
            sens.append(s)
            if i==0:top.append(r)
            if i==1:wins.append(s)
        b=top[-1];w=wins[-1]
        rs.append({'VariationID':vid,'Stage2E_v2_top_Guide_ID':b['Guide_ID'],'Stage2E_v2_top_score_0_100':b[m.SCORE],
                   'Editable_A_offtarget_hits_le3mm':'0'})
        sm.append({'VariationID':vid,'Stage2E_v2_top_Guide_ID':b['Guide_ID'],'Stage2F_v2_robust_top_Guide_ID':w['Guide_ID'],
                   'Stage2E_top_remains_robust_top':'NO','Robust_top_random_top_fraction':'1.0','Robust_top_mean_rank':'1',
                   'Robust_top_scenario_top_count_of_6':'6','Robust_top_stability_label':'HIGH_STABILITY',
                   'Robust_top_editable_A_offtarget_hits_le3mm':'0','Robust_top_Stage2E_v2_score_0_100':w[m.SCORE]})
        for scenario in ['baseline','equal','specificity_heavy','bystander_heavy','accessibility_heavy','structure_heavy']:
            d={k:'0.2' for k in schema['scenarios']};d.update(VariationID=vid,Scenario=scenario,Top_Guide_ID=w['Guide_ID']);sc.append(d)
    return {'ranked':ranked,'baseline_top':top,'ranking_summary':rs,'sensitivity':sens,'robust_top':wins,'robust_summary':sm,'scenarios':sc},schema

class MappingTests(unittest.TestCase):
    def setUp(self):
        self.t,self.schema=fixtures();self.norm={'cohort_size':376,'duplex_MFE_range':[-10,0],'guide_MFE_range':[-5,0]}
        self.cfg={'ranking_run_directory':'SYNTHETIC_RANK','robustness_run_directory':'SYNTHETIC_ROB','run_id':'SYNTHETIC_ONLY',
                  'inputs':{k:{'sha256':'a'*64} for k in ['ranked','robust_top','sensitivity','ranking_manifest','robustness_manifest']}}
    def validate(self):return m.validate_tables(self.t,self.norm)
    def test_valid_input_contract(self):self.validate()
    def test_one_per_mutation_and_identity(self):
        rows=m.assemble(self.t,self.cfg,'c'*64)
        self.assertEqual(len(rows),8);self.assertEqual({r['VariationID'] for r in rows},m.u.MUTATIONS)
        self.assertEqual({r['Guide_ID'] for r in rows},{r['Guide_ID'] for r in self.t['robust_top']})
    def test_all_source_cells_retained(self):
        for row in m.assemble(self.t,self.cfg,'c'*64):
            source=next(r for r in self.t['robust_top'] if r['Guide_ID']==row['Guide_ID'])
            self.assertTrue(all(row[k]==v for k,v in source.items()))
    def test_baseline_separate_from_robust(self):
        for r in m.assemble(self.t,self.cfg,'c'*64):
            self.assertNotEqual(r['Guide_ID'],r['Baseline_Guide_ID']);self.assertEqual(r['Robust_matches_baseline'],'NO')
            self.assertEqual(r['Baseline_random_top_fraction'],'0.0');self.assertEqual(r['Sensitivity_random_top_fraction'],'1.0')
    def test_provenance_complete(self):
        for r in m.assemble(self.t,self.cfg,'c'*64):
            self.assertEqual(r['Assembly_config_sha256'],'c'*64);self.assertEqual(r['Source_robust_top_sha256'],'a'*64)
            self.assertEqual(r['Events_log_independent_verification'],'NOT_SUPPLIED');self.assertEqual(r['Biological_validation_status'],'NOT_ESTABLISHED')
    def test_deterministic_order(self):
        a=m.assemble(self.t,self.cfg,'c'*64);self.t['robust_top'].reverse()
        self.assertEqual(a,m.assemble(self.t,self.cfg,'c'*64))
    def test_missing_candidate(self):
        self.t['robust_top'].pop()
        with self.assertRaisesRegex(ValueError,'Incomplete'):self.validate()
    def test_duplicate_candidate(self):
        self.t['robust_top'][1]=self.t['robust_top'][0]
        with self.assertRaisesRegex(ValueError,'Duplicate'):self.validate()
    def test_missing_cohort_guide(self):
        self.t['sensitivity'].pop()
        with self.assertRaisesRegex(ValueError,'cohort'):self.validate()
    def test_mutation_conflict(self):
        self.t['robust_top'][0]=dict(self.t['robust_top'][0]);self.t['robust_top'][0]['VariationID']='999'
        with self.assertRaises(ValueError):self.validate()
    def test_winner_mismatch(self):
        self.t['robust_top'][0]=self.t['sensitivity'][0]
        with self.assertRaisesRegex(ValueError,'winner mismatch'):self.validate()
    def test_ranking_lineage_conflict(self):
        self.t['sensitivity'][0]['Transcript']='OTHER'
        with self.assertRaisesRegex(ValueError,'conflict'):self.validate()
    def test_nonfinite_values(self):
        for v in ['NaN','inf','-Infinity','']:
            t=copy.deepcopy(self.t);t['sensitivity'][0]['Sensitivity_mean_rank']=v
            with self.subTest(value=v),self.assertRaises(ValueError):m.validate_tables(t,self.norm)
    def test_nonfinite_inherited_feature(self):
        for role in ['ranked','sensitivity']:self.t[role][0]['Target_A_unpaired_probability']='NaN'
        with self.assertRaisesRegex(ValueError,'Nonfinite'):self.validate()
    def test_invalid_coordinate(self):
        for role in ['ranked','sensitivity']:self.t[role][0]['Fullwindow_catalytic_query_position_1based']='20'
        with self.assertRaisesRegex(ValueError,'coordinate'):self.validate()
    def test_frequency_denominator(self):
        self.t['sensitivity'][0]['Sensitivity_random_top_fraction']='.5'
        with self.assertRaisesRegex(ValueError,'metrics'):self.validate()
    def test_missing_or_ambiguous_fields(self):
        for b in [b'Guide_ID\nG\n',b'Guide_ID,Guide_ID\nG,G\n']:
            with self.assertRaises(ValueError):m.u.table(b,self.schema['ranked'])

class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        self.stage=self.root/'staging';self.stage.mkdir();self.t,self.schema=fixtures();self.paths={}
        self.cfg={'version':m.VERSION,'selection_policy':m.POLICY,'stages':dict(m.STAGES),'staging_root':str(self.stage),
                  'run_id':'SYNTHETIC_ASSEMBLY','protected_roots':[str(self.root/'official')],
                  'ranking_run_directory':'SYNTHETIC_RANK','robustness_run_directory':'SYNTHETIC_ROB',
                  'python_sha256':m.u.digest(Path(sys.executable).read_bytes()),'ranking_config_sha256':'1'*64,
                  'robustness_config_sha256':'2'*64,'expected_winners':{r['VariationID']:r['Guide_ID'] for r in self.t['robust_top']},
                  'uninspected_artifacts':{'events.jsonl':{'sha256':'f'*64,'bytes':187}},'inputs':{}}
        self.names={'ranked':'ranked.csv','baseline_top':'top_guides.csv','ranking_summary':'mutation_summary.csv','normalization':'normalization.json',
               'ranking_manifest':'output_manifest.json','ranking_result':'run_result.json','ranking_run_manifest':'run_manifest.json',
               'sensitivity':'sensitivity_all_guides.csv','robust_top':'robust_top_guides.csv','robust_summary':'mutation_summary.csv',
               'scenarios':'scenarios.csv','qc':'simulation_qc.json','robustness_manifest':'output_manifest.json',
               'robustness_result':'run_result.json','robustness_run_manifest':'run_manifest.json'}
        for role in self.names:
            d=self.root/('ranking' if role in ['ranked','baseline_top','ranking_summary','normalization'] or role.startswith('ranking_') else 'robustness');d.mkdir(exist_ok=True)
            self.paths[role]=d/self.names[role]
        for role in self.t:m.u.write_csv(self.paths[role],self.t[role],self.schema[role])
        m.u.write_json(self.paths['normalization'],{'cohort_size':376,'duplex_MFE_range':[-10,0],'guide_MFE_range':[-5,0]})
        m.u.write_json(self.paths['qc'],{'guides':376,'mutations':8,'simulations':10000,'seed':20260914,'fixed_cells_unchanged':33464,'scenario_rows':48,'winner_and_rank_sum_conservation':'PASS'})
        for prefix,stage in [('ranking','composite-ranking'),('robustness','robustness')]:
            m.u.write_json(self.paths[prefix+'_run_manifest'],{'stage':stage,'config_sha256':self.cfg[prefix+'_config_sha256'],'python_sha256':self.cfg['python_sha256']})
            roles=['ranked','baseline_top','ranking_summary','normalization','ranking_run_manifest'] if prefix=='ranking' else ['sensitivity','robust_top','robust_summary','scenarios','qc','robustness_run_manifest']
            entries=[{'path':self.names[r],'sha256':m.u.digest(self.paths[r].read_bytes()),'bytes':self.paths[r].stat().st_size} for r in roles]
            if prefix=='robustness':entries.append({'path':'events.jsonl','sha256':'f'*64,'bytes':187})
            m.u.write_json(self.paths[prefix+'_manifest'],{'files':entries})
            m.u.write_json(self.paths[prefix+'_result'],{'stage':stage,'status':'COMPLETE','exit_code':0,'preservation':'PASS','downstream_stages_executed':[],
                                                      'output_manifest_sha256':m.u.digest(self.paths[prefix+'_manifest'].read_bytes())})
        self.paths.update(adapter=HERE/'mecp2_tier1_assembly.py',io=HERE/'assembly_io.py',schema=HERE/'input_schema.json')
        for role,p in self.paths.items():self.cfg['inputs'][role]={'path':str(p),'bytes':p.stat().st_size,'sha256':m.u.digest(p.read_bytes())}
        self.cp=self.root/'config.json';self.save();self.ap=self.root/'approval.json';self.approve()
    def save(self):self.cp.write_text(json.dumps(self.cfg));self.ch=m.u.digest(self.cp.read_bytes())
    def approve(self):
        self.ap.write_text(json.dumps({'approved':True,'approved_by':'SYNTHETIC_TEST_ONLY','stage':'assemble-tier1','config_sha256':self.ch,
                      'run_directory':str(self.stage/self.cfg['run_id']),'python_sha256':self.cfg['python_sha256'],'promotion_approved':False,
                      'expires_utc':(dt.datetime.now(dt.timezone.utc)+dt.timedelta(hours=1)).isoformat()}))
    def execute(self,stage='assemble-tier1'):return m.execute(stage,self.cp,self.ch,self.ap)
    def test_preflight_does_not_assemble(self):
        with patch.object(m,'assemble',side_effect=AssertionError('No real assembly in preflight')):
            r=self.execute('verify-inputs');self.assertEqual(r['status'],'VERIFIED_FOR_APPROVAL_REVIEW');self.assertFalse(r['library_created'])
        self.assertEqual(list(self.stage.iterdir()),[])
    def test_complete_synthetic_assembly(self):
        before={r:p.read_bytes() for r,p in self.paths.items()};out=self.execute();self.assertEqual(out['status'],'COMPLETE')
        run=self.stage/self.cfg['run_id'];self.assertEqual({p.name for p in run.iterdir()},m.FILES)
        rows,fields=m.u.table((run/'tier1_candidate_library.csv').read_bytes());self.assertEqual((len(rows),len(fields)),(8,118))
        for e in json.loads((run/'output_manifest.json').read_text())['files']:self.assertEqual(m.u.digest((run/e['path']).read_bytes()),e['sha256'])
        self.assertEqual(before,{r:p.read_bytes() for r,p in self.paths.items()})
    def test_missing_approval(self):
        with self.assertRaisesRegex(ValueError,'approval'):m.execute('assemble-tier1',self.cp,self.ch)
    def test_unapproved_record(self):
        a=json.loads(self.ap.read_text());a['approved']=False;self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'Approval'):self.execute()
    def test_wrong_approval_scope(self):
        a=json.loads(self.ap.read_text());a['stage']='robustness';self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'scope'):self.execute()
    def test_expired_approval(self):
        a=json.loads(self.ap.read_text());a['expires_utc']='2000-01-01T00:00:00Z';self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'expired'):self.execute()
    def test_stale_input(self):
        self.paths['robust_top'].write_text('CHANGED')
        with self.assertRaisesRegex(ValueError,'Stale/hash'):self.execute()
    def test_missing_input(self):
        self.cfg['inputs']['robust_top']['path']=str(self.root/'ABSENT');self.save()
        with self.assertRaisesRegex(ValueError,'unavailable'):self.execute()
    def test_incomplete_manifest(self):
        p=self.paths['ranking_manifest'];d=json.loads(p.read_text());d['files'].pop();p.write_text(json.dumps(d))
        self.cfg['inputs']['ranking_manifest'].update(sha256=m.u.digest(p.read_bytes()),bytes=p.stat().st_size);self.save()
        with self.assertRaisesRegex(ValueError,'manifest'):self.execute()
    def test_manifest_duplicate(self):
        p=self.paths['robustness_manifest'];d=json.loads(p.read_text());d['files'].append(d['files'][0]);p.write_text(json.dumps(d))
        self.cfg['inputs']['robustness_manifest'].update(sha256=m.u.digest(p.read_bytes()),bytes=p.stat().st_size);self.save()
        with self.assertRaisesRegex(ValueError,'manifest'):self.execute()
    def test_wrong_expected_winner(self):
        self.cfg['expected_winners']['143546']='T1G0000';self.save()
        with self.assertRaisesRegex(ValueError,'Expected winner'):self.execute()
    def test_existing_directory(self):
        (self.stage/self.cfg['run_id']).mkdir()
        with self.assertRaisesRegex(ValueError,'already exists'):self.execute()
    def test_wrong_config_hash(self):
        with self.assertRaisesRegex(ValueError,'config hash'):m.execute('verify-inputs',self.cp,'0'*64)
    def test_interpreter_hash(self):
        self.cfg['python_sha256']='0'*64;self.save()
        with self.assertRaisesRegex(ValueError,'Interpreter'):self.execute()
    def test_prohibited_root(self):
        self.cfg['protected_roots']=[str(self.stage)];self.save()
        with self.assertRaisesRegex(ValueError,'Prohibited'):self.execute()
    def test_promotion_and_scientific_commands_absent(self):
        for s in ['promotion','ranking','robustness','blast']:
            with self.assertRaisesRegex(ValueError,'not permitted'):self.execute(s)
    def test_preservation_failure(self):
        original=m.assemble
        def changed(*a,**k):
            out=original(*a,**k)
            if os.name!='nt':self.paths['robust_top'].write_text('SYNTHETIC CHANGE')
            else:raise ValueError('Preservation failure injected')
            return out
        with patch.object(m,'assemble',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'Preservation failure'):self.execute()
    def test_incomplete_assembly(self):
        with patch.object(m,'assemble',return_value=[]):
            with self.assertRaisesRegex(ValueError,'coverage'):self.execute()
    def test_write_failure(self):
        with patch.object(m.u,'write_csv',side_effect=OSError('SYNTHETIC write failure')):
            with self.assertRaises(OSError):self.execute()
    def test_wrong_approval_hash(self):
        a=json.loads(self.ap.read_text());a['config_sha256']='0'*64;self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'hash/path'):self.execute()
    def test_partial_written_library(self):
        original=m.u.write_csv
        def truncated(path,rows,fields):original(path,rows[:-1],fields)
        with patch.object(m.u,'write_csv',side_effect=truncated):
            with self.assertRaisesRegex(ValueError,'Incomplete/wrong written'):self.execute()
        result=json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text())
        self.assertEqual((result['status'],result['exit_code']),('FAILED_OR_INTERRUPTED',4))
    def test_manifest_digest_mismatch(self):
        p=self.paths['ranking_manifest'];d=json.loads(p.read_text());d['files'][0]['sha256']='0'*64;p.write_text(json.dumps(d))
        self.cfg['inputs']['ranking_manifest'].update(sha256=m.u.digest(p.read_bytes()),bytes=p.stat().st_size);self.save()
        with self.assertRaisesRegex(ValueError,'Manifest hash'):self.execute()
    def test_status_exit_consistency(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:code=m.main(['verify-inputs','--config',str(self.cp),'--config-sha256','0'*64])
        self.assertEqual(code,4);self.assertEqual(json.loads(out.getvalue())['exit_code'],code)

if __name__=='__main__':unittest.main(verbosity=2)
