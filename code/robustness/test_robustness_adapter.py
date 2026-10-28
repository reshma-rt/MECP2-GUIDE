"""Synthetic-only unit/integration tests. No real ranked input is read."""
import ast
import contextlib
import copy
import datetime as dt
import io
import json
import os
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch
import mecp2_robustness_adapter as m
from synthetic_fixtures import fixtures

HERE=Path(__file__).resolve().parent

class MethodTests(unittest.TestCase):
    def setUp(self):
        self.rows,self.top,self.summary,self.norm,self.schema=fixtures()

    def validate(self):
        return m.validate_baseline(self.rows,self.top,self.summary,self.norm)

    def test_baseline_contract(self): self.validate()

    def test_seed_determinism(self):
        a=m.core.calculate(self.rows,n_sim=31,seed=20260914)
        b=m.core.calculate(self.rows,n_sim=31,seed=20260914)
        self.assertEqual(a,b)

    def test_exact_source_rng(self):
        tree=ast.parse((HERE/'source_snapshots/stage2F_v2_tier1_harmonized_sensitivity.py').read_text())
        fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='dirichlet_ones')
        scope={};exec(compile(ast.Module(body=[fn],type_ignores=[]),'source_rng','exec'),scope)
        a=random.Random(20260914);b=random.Random(20260914)
        for i in range(100):
            x=m.core.dirichlet_ones(5,a);self.assertEqual(x,scope['dirichlet_ones'](5,b))
            self.assertTrue(all(0 < v < 1 for v in x));self.assertAlmostEqual(sum(x),1)

    def test_scientific_body_ast_identical(self):
        original=ast.parse((HERE/'source_snapshots/stage2F_v2_tier1_harmonized_sensitivity.py').read_text())
        start=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='guide_data' for t in n.targets))
        end=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='extra_fields' for t in n.targets))
        engine=ast.parse((HERE/'robustness_core.py').read_text())
        fn=next(n for n in engine.body if isinstance(n,ast.FunctionDef) and n.name=='calculate')
        self.assertEqual(ast.dump(ast.Module(body=original.body[start:end],type_ignores=[]),include_attributes=False),
                         ast.dump(ast.Module(body=fn.body[2:-1],type_ignores=[]),include_attributes=False))

    def test_parity_with_original_full_source_synthetic(self):
        source=HERE/'source_snapshots/stage2F_v2_tier1_harmonized_sensitivity.py'
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);pipe=root/'02_tier1_pipeline';pipe.mkdir()
            m.util.write_csv(pipe/'MECP2_ADAR_Dataset2_Tier1_stage2E_v2_ranked.csv',self.rows,self.schema['ranked'])
            tree=ast.parse(source.read_text())
            for n in tree.body:
                if isinstance(n,ast.Assign):
                    names=[t.id for t in n.targets if isinstance(t,ast.Name)]
                    if 'ROOT' in names:n.value=ast.Call(func=ast.Name(id='Path',ctx=ast.Load()),args=[ast.Constant(str(root))],keywords=[])
                    if 'N_SIM' in names:n.value=ast.Constant(23)
            ast.fix_missing_locations(tree)
            with contextlib.redirect_stdout(io.StringIO()):exec(compile(tree,str(source),'exec'),{'__file__':str(source)})
            ours=m.core.calculate(self.rows,n_sim=23,seed=20260914)
            names=['sensitivity_all_guides','robust_top_guides','scenarios','mutation_summary']
            for name,records in zip(names,ours[:4]):
                old,_=m.util.table((pipe/('MECP2_ADAR_Dataset2_Tier1_stage2F_v2_'+name+'.csv')).read_bytes())
                self.assertEqual(old,[{k:str(v) for k,v in r.items()} for r in records])

    def test_denominators_and_mean_ranks(self):
        output=m.core.calculate(self.rows,n_sim=19)
        qc=m.reconcile(self.rows,output,n_sim=19)
        self.assertEqual(qc['winner_and_rank_sum_conservation'],'PASS')
        for r in output[0]:
            self.assertEqual(r['Sensitivity_random_top_fraction'],r['Sensitivity_random_top_count']/19)
            self.assertEqual(r['Sensitivity_mean_rank'],output[5][r['Guide_ID']]/19)

    def test_all_ties_use_guide_id_not_editable_count(self):
        for r in self.rows:
            for c in m.core.COMPONENTS:r[c]='0.5'
            r['Fullwindow_editable_A_offtarget_hits_le3mm']='0' if r['Guide_ID'].endswith('2') else '9'
        output=m.core.calculate(self.rows,n_sim=3)
        for vid in m.util.MUTATIONS:
            expected=min(r['Guide_ID'] for r in self.rows if r['VariationID']==vid)
            sc=[r for r in output[2] if r['VariationID']==vid]
            self.assertTrue(all(r['Top_Guide_ID']==expected for r in sc))
            ann=next(r for r in output[0] if r['Guide_ID']==expected)
            self.assertEqual(ann['Sensitivity_random_top_count'],3)
            self.assertEqual(ann['Sensitivity_mean_rank'],1)

    def test_threshold_edges_match_source(self):
        tree=ast.parse((HERE/'source_snapshots/stage2F_v2_tier1_harmonized_sensitivity.py').read_text())
        node=next(n for n in ast.walk(tree) if isinstance(n,ast.If) and ast.unparse(n.test)=='top_fraction >= 0.5')
        for f in [0,.199999999,.2,.499999999,.5,1]:
            scope={'top_fraction':f};exec(compile(ast.Module(body=[node],type_ignores=[]),'source_label','exec'),scope)
            self.assertEqual(m.label(f),scope['label'])

    def test_robust_winner_uses_baseline_before_id_but_display_sort_does_not(self):
        expected={}
        for vid in m.util.MUTATIONS:
            group=sorted([r for r in self.rows if r['VariationID']==vid],key=lambda r:r['Guide_ID'])
            for r in group:
                for c in m.core.COMPONENTS:r[c]='0'
                r[m.core.BASELINE_SCORE]='0'
            group[0][m.core.COMPONENTS[0]]='1';group[0][m.core.BASELINE_SCORE]='25'
            group[1][m.core.COMPONENTS[1]]='1';group[1][m.core.BASELINE_SCORE]='30'
            expected[vid]=(group[0]['Guide_ID'],group[1]['Guide_ID'])
        with patch.object(m.core,'dirichlet_ones',side_effect=[[1,0,0,0,0],[0,1,0,0,0]]):
            output=m.core.calculate(self.rows,n_sim=2)
        for vid,(display_gid,winner_gid) in expected.items():
            display=next(r for r in output[0] if r['VariationID']==vid)
            winner=next(r for r in output[1] if r['VariationID']==vid)
            self.assertEqual(display['Guide_ID'],display_gid)
            self.assertEqual(winner['Guide_ID'],winner_gid)

    def test_no_energy_renormalization_in_F(self):
        original=m.core.calculate(self.rows,n_sim=3)
        for r in self.rows:r['Guide_target_duplex_MFE_kcal_mol']='-5';r['Guide_MFE_kcal_mol']='-5'
        changed=m.core.calculate(self.rows,n_sim=3)
        # Scientific F uses supplied components, not raw energy columns.
        self.assertEqual([r['Guide_ID'] for r in original[1]],[r['Guide_ID'] for r in changed[1]])
        for x,y in zip(original[0],changed[0]):
            for c in m.EXTRA:self.assertEqual(x[c],y[c])

    def test_constant_energy_baseline_convention(self):
        for r in self.rows:
            r['Guide_target_duplex_MFE_kcal_mol']='-5';r['Guide_MFE_kcal_mol']='-5'
            r[m.core.COMPONENTS[3]]='0.0';r[m.core.COMPONENTS[4]]='1.0'
        # Compare arithmetic convention without generating an alternative real baseline.
        self.norm['duplex_MFE_range']=[-5,-5];self.norm['guide_MFE_range']=[-5,-5]
        # Existing synthetic scores/ranks intentionally stale: fail rather than repair them.
        with self.assertRaisesRegex(ValueError,'baseline score'):self.validate()

    def test_duplicate_ids(self):
        self.rows[1]['Guide_ID']=self.rows[0]['Guide_ID']
        with self.assertRaisesRegex(ValueError,'Duplicate'):self.validate()

    def test_missing_guide(self):
        self.rows.pop()
        with self.assertRaisesRegex(ValueError,'cohort'):self.validate()

    def test_wrong_mutation_group(self):
        self.rows[0]['VariationID']=self.rows[-1]['VariationID']
        with self.assertRaisesRegex(ValueError,'mutation groups'):self.validate()

    def test_nonfinite(self):
        for v in ['NaN','NAN','inf','-Infinity','','NA']:
            rows=copy.deepcopy(self.rows);rows[0][m.core.COMPONENTS[0]]=v
            with self.subTest(v=v),self.assertRaises(ValueError):m.validate_baseline(rows,self.top,self.summary,self.norm)

    def test_out_of_range(self):
        self.rows[0][m.core.COMPONENTS[0]]='1.0001'
        with self.assertRaisesRegex(ValueError,'outside'):self.validate()

    def test_conflicting_coordinate(self):
        self.rows[0]['Fullwindow_catalytic_query_position_1based']='20'
        with self.assertRaisesRegex(ValueError,'coordinate conflict'):self.validate()

    def test_mutation_target_metadata_conflict(self):
        self.rows[1]['Transcript']='CONFLICT'
        with self.assertRaisesRegex(ValueError,'metadata conflict'):self.validate()

    def test_valid_constant_energy_baseline(self):
        rows,top,summary,norm,_=fixtures(constant=True)
        m.validate_baseline(rows,top,summary,norm)
        result=m.core.calculate(rows,n_sim=2)
        m.reconcile(rows,result,n_sim=2)

    def test_conflicting_target_sequence(self):
        self.rows[0]['Target_segment_RNA_5to3']='U'*21
        with self.assertRaisesRegex(ValueError,'sequence'):self.validate()

    def test_conflicting_component(self):
        self.rows[0][m.core.COMPONENTS[2]]='0.123456'
        with self.assertRaisesRegex(ValueError,'component mismatch'):self.validate()

    def test_conflicting_baseline_score(self):
        self.rows[0][m.core.BASELINE_SCORE]='99'
        with self.assertRaisesRegex(ValueError,'baseline score'):self.validate()

    def test_conflicting_baseline_rank(self):
        self.rows[0][m.core.BASELINE_RANK]='2'
        with self.assertRaisesRegex(ValueError,'ranks'):self.validate()

    def test_conflicting_top_record(self):
        self.top[0]=dict(self.top[0]);self.top[0]['Transcript']='OTHER'
        with self.assertRaisesRegex(ValueError,'top record'):self.validate()

    def test_conflicting_summary(self):
        self.summary[0]['Stage2E_v2_top_score_0_100']='0'
        with self.assertRaisesRegex(ValueError,'summary'):self.validate()

    def test_incomplete_output(self):
        output=m.core.calculate(self.rows,n_sim=2);output[0].pop()
        with self.assertRaisesRegex(ValueError,'Incomplete'):m.reconcile(self.rows,output,n_sim=2)

    def test_bad_frequency(self):
        output=m.core.calculate(self.rows,n_sim=2);output[0][0]['Sensitivity_random_top_fraction']=.123
        with self.assertRaisesRegex(ValueError,'denominator'):m.reconcile(self.rows,output,n_sim=2)

    def test_bad_scenario(self):
        output=m.core.calculate(self.rows,n_sim=2);output[2][0]['Top_Guide_ID']='T1G0000'
        with self.assertRaisesRegex(ValueError,'Scenario'):m.reconcile(self.rows,output,n_sim=2)

    def test_missing_ambiguous_extra_columns(self):
        for headers in [self.schema['ranked'][:-1],self.schema['ranked']+['Unexpected'],self.schema['ranked']+['Guide_ID']]:
            raw=(','.join(headers)+'\n').encode()
            with self.assertRaises(ValueError):m.util.table(raw,self.schema['ranked'])

class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.stage=self.root/'stage';self.stage.mkdir()
        rows,top,summary,norm,schema=fixtures()
        paths={}
        for role,records in [('ranked',rows),('top',top),('summary',summary)]:
            name={'ranked':'ranked.csv','top':'top_guides.csv','summary':'mutation_summary.csv'}[role]
            paths[role]=self.root/name;m.util.write_csv(paths[role],records,schema[role])
        paths['normalization']=self.root/'normalization.json';m.util.write_json(paths['normalization'],norm)
        self.ph=m.util.digest(Path(sys_executable()).read_bytes())
        paths['ranking_run_manifest']=self.root/'run_manifest.json'
        m.util.write_json(paths['ranking_run_manifest'],{'stage':'composite-ranking','config_sha256':'a'*64,'python_sha256':self.ph})
        entries=[{'path':p.name,'bytes':p.stat().st_size,'sha256':m.util.digest(p.read_bytes())} for p in paths.values()]
        paths['ranking_manifest']=self.root/'output_manifest.json';m.util.write_json(paths['ranking_manifest'],{'files':entries})
        paths['ranking_result']=self.root/'run_result.json'
        m.util.write_json(paths['ranking_result'],{'status':'COMPLETE','stage':'composite-ranking','preservation':'PASS','exit_code':0,
                        'downstream_stages_executed':[],'output_manifest_sha256':m.util.digest(paths['ranking_manifest'].read_bytes())})
        paths.update(adapter=HERE/'mecp2_robustness_adapter.py',core=HERE/'robustness_core.py',io=HERE/'robustness_io.py',
                     schema=HERE/'schema_contract.json',ranking_source=HERE/'source_snapshots/stage2E_v2_tier1_harmonized_ranking.py',
                     robustness_source=HERE/'source_snapshots/stage2F_v2_tier1_harmonized_sensitivity.py')
        self.cfg={'version':m.VERSION,'simulations':10000,'seed':20260914,'components':list(m.core.COMPONENTS),
                  'scenarios':copy.deepcopy(m.core.SCENARIOS),'stages':dict(m.STAGES),'python_sha256':self.ph,
                  'ranking_config_sha256':'a'*64,'staging_root':str(self.stage),'run_id':'SYNTHETIC_ONLY',
                  'protected_roots':[str(self.root/'official')],
                  'inputs':{role:{'path':str(p),'sha256':m.util.digest(p.read_bytes()),'bytes':p.stat().st_size} for role,p in paths.items()}}
        self.cp=self.root/'config.json';self.save();self.ap=self.root/'approval.json';self.approve()

    def save(self):
        self.cp.write_text(json.dumps(self.cfg));self.ch=m.util.digest(self.cp.read_bytes())

    def approve(self):
        self.ap.write_text(json.dumps({'approved':True,'approved_by':'SYNTHETIC_TEST_ONLY','stage':'robustness','config_sha256':self.ch,
            'run_directory':str(self.stage/self.cfg['run_id']),'python_sha256':self.ph,'simulations':10000,'seed':20260914,
            'final_library_approved':False,'promotion_approved':False,
            'expires_utc':(dt.datetime.now(dt.timezone.utc)+dt.timedelta(hours=1)).isoformat()}))

    def execute(self,stage='robustness'):
        return m.execute(stage,self.cp,self.ch,self.ap)

    def test_read_only_preflight(self):
        with patch.object(m.core,'calculate',side_effect=AssertionError('No simulation in preflight')):
            self.assertEqual(self.execute('verify-inputs')['simulations_executed'],0)
        self.assertEqual(list(self.stage.iterdir()),[])

    def test_complete_10000_synthetic_only(self):
        before={role:Path(v['path']).read_bytes() for role,v in self.cfg['inputs'].items()}
        result=self.execute();self.assertEqual((result['status'],result['exit_code']),('COMPLETE',0))
        self.assertEqual(result['simulations'],10000)
        run=self.stage/self.cfg['run_id'];self.assertEqual({p.name for p in run.iterdir()},m.OUTPUT_FILES)
        for e in json.loads((run/'output_manifest.json').read_text())['files']:
            self.assertEqual(m.util.digest((run/e['path']).read_bytes()),e['sha256'])
        self.assertEqual(before,{role:Path(v['path']).read_bytes() for role,v in self.cfg['inputs'].items()})

    def test_missing_approval(self):
        with self.assertRaisesRegex(ValueError,'approval'):m.execute('robustness',self.cp,self.ch)

    def test_unapproved_record(self):
        a=json.loads(self.ap.read_text());a['approved']=False;self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'not granted'):self.execute()

    def test_wrong_approval_hash(self):
        a=json.loads(self.ap.read_text());a['config_sha256']='0'*64;self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'hash/path'):self.execute()

    def test_expired_approval(self):
        a=json.loads(self.ap.read_text());a['expires_utc']='2000-01-01T00:00:00Z';self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'expired'):self.execute()

    def test_wrong_scope(self):
        a=json.loads(self.ap.read_text());a['stage']='composite-ranking';self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'stage/hash'):self.execute()

    def test_existing_output_directory(self):
        (self.stage/self.cfg['run_id']).mkdir()
        with self.assertRaisesRegex(ValueError,'already exists'):self.execute()

    def test_stale_hash(self):
        Path(self.cfg['inputs']['ranked']['path']).write_text('CHANGED')
        with self.assertRaisesRegex(ValueError,'Stale/size/hash'):self.execute()

    def test_missing_file(self):
        self.cfg['inputs']['ranked']['path']=str(self.root/'ABSENT');self.save()
        with self.assertRaisesRegex(ValueError,'unavailable'):self.execute()

    def test_unset_binding(self):
        self.cfg['inputs']['ranked']['path']=None;self.save()
        with self.assertRaisesRegex(ValueError,'Unassigned'):self.execute()

    def test_wrong_config_hash(self):
        with self.assertRaisesRegex(ValueError,'configuration hash'):m.execute('verify-inputs',self.cp,'0'*64)

    def test_method_change_rejected(self):
        self.cfg['seed']=123;self.save()
        with self.assertRaisesRegex(ValueError,'Method/schema'):self.execute()

    def test_interpreter_hash(self):
        self.cfg['python_sha256']='0'*64;self.save()
        with self.assertRaisesRegex(ValueError,'Interpreter'):self.execute()

    def test_output_root_protected(self):
        self.cfg['protected_roots']=[str(self.stage)];self.save()
        with self.assertRaisesRegex(ValueError,'Prohibited'):self.execute()

    def test_other_stages_rejected(self):
        for s in ['composite-ranking','final-library','promotion','blast']:
            with self.subTest(stage=s),self.assertRaisesRegex(ValueError,'prohibited'):self.execute(s)

    def test_interruption_failure_status(self):
        with patch.object(m.core,'calculate',side_effect=KeyboardInterrupt('SYNTHETIC interruption')):
            with self.assertRaises(KeyboardInterrupt):self.execute()
        result=json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text())
        self.assertEqual((result['status'],result['exit_code']),('FAILED_OR_INTERRUPTED',4))

    def test_incomplete_computation(self):
        with patch.object(m.core,'calculate',return_value=([],[],[],[],{},{})):
            with self.assertRaisesRegex(ValueError,'Incomplete'):self.execute()

    def cheap_output(self,*args,**kwargs):
        return self.original(args[0],n_sim=3,seed=20260914)

    def test_preservation_failure(self):
        self.original=m.core.calculate
        def altered(*args,**kwargs):
            result=self.cheap_output(*args,**kwargs)
            if os.name!='nt':Path(self.cfg['inputs']['ranked']['path']).write_text('SYNTHETIC MUTATION')
            else:raise ValueError('Preservation failure injected')
            return result
        # A fixture-only short calculation and QC stub exercise boundary preservation, not scientific conclusions.
        with patch.object(m.core,'calculate',side_effect=altered),patch.object(m,'reconcile',return_value={'synthetic':True}):
            with self.assertRaisesRegex(ValueError,'Preservation failure'):self.execute()

    def test_incomplete_written_output(self):
        self.original=m.core.calculate;write=m.util.write_csv
        def truncate(path,rows,fields):write(path,rows[:-1],fields)
        with patch.object(m.core,'calculate',side_effect=self.cheap_output),patch.object(m,'reconcile',return_value={'synthetic':True}),patch.object(m.util,'write_csv',side_effect=truncate):
            with self.assertRaisesRegex(ValueError,'Incomplete/wrong written'):self.execute()

    def test_failed_write(self):
        self.original=m.core.calculate
        with patch.object(m.core,'calculate',side_effect=self.cheap_output),patch.object(m,'reconcile',return_value={'synthetic':True}),patch.object(m.util,'write_csv',side_effect=OSError('SYNTHETIC disk failure')):
            with self.assertRaises(OSError):self.execute()

    def test_status_exit_agree(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            code=m.main(['verify-inputs','--config',str(self.cp),'--config-sha256','0'*64])
        self.assertEqual(code,4);self.assertEqual(json.loads(out.getvalue())['exit_code'],code)

    def test_exclusive_write(self):
        p=self.root/'keep.json';p.write_text('KEEP')
        with self.assertRaises(FileExistsError):m.util.write_json(p,{})
        self.assertEqual(p.read_text(),'KEEP')

def sys_executable():
    import sys
    return sys.executable

if __name__=='__main__':unittest.main(verbosity=2)
