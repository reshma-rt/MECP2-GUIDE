"""Artificial fixtures only. Does not read real guide tables."""
import copy
import datetime as dt
import contextlib
import io
import json
import os
from pathlib import Path
import random
import sys
import tempfile
import unittest
from unittest.mock import patch
import mecp2_ablation as m

HERE=Path(__file__).resolve().parent

def fixtures():
    schema=json.loads((HERE/'schema.json').read_text());rows=[];robust=[]
    for vi,vid in enumerate(sorted(m.u.MUTATIONS)):
        group=[]
        for i in range(47):
            r={c:'' for c in schema['ranked']}
            r.update(Guide_ID='T1G%04d'%(vi*47+i+1),VariationID=vid,Rescue_class='EXACT_RESTORATION',Transcript='SYNTHETIC_ONLY',
                     Guide_length_nt='21',Target_A_position_in_guide_1based='6',Target_A_position_in_target_segment_1based='16',
                     Fullwindow_catalytic_query_position_1based='16',Stage2E_v2_status='HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE')
            for c,x in zip(m.COLS,[.2+vi*.05,(i%7)/8,(i%11)/12,i/46,(i%3)/4]):r[c]=str(x)
            for c,w in zip(m.NAMES,m.BASE):r['Stage2E_v2_weight_'+c]=str(w)
            r[m.COUNT]=str(i%3);r[m.SCORE]=str(m.left_score(r,m.BASE));group.append(r)
        group.sort(key=lambda r:(-float(r[m.SCORE]),int(r[m.COUNT]),r['Guide_ID']))
        for rank,r in enumerate(group,1):r[m.RANK]=str(rank)
        rows.extend(group);rr={c:group[1].get(c,'SYNTHETIC_ONLY') for c in schema['robust_top']};robust.append(rr)
    return rows,robust,schema

class MethodTests(unittest.TestCase):
    def setUp(self):self.rows,self.robust,self.schema=fixtures()
    def test_contract(self):m.validate(self.rows,self.robust)
    def test_nine_predefined_scenarios(self):
        d=m.definitions();self.assertEqual(len(d),9);self.assertNotIn('only_accessibility',d)
        for w in d.values():self.assertAlmostEqual(sum(w),1)
    def test_leave_one_out_weights(self):
        for i in range(1,5):
            w=m.definitions()['without_'+m.NAMES[i]];self.assertEqual(w[i],0)
            for j in range(5):
                if j!=i:self.assertEqual(w[j],m.BASE[j]/(1-m.BASE[i]))
    def test_exact_known_score(self):
        row=dict(zip(m.COLS,['.5','.4','.8','1','1']));self.assertAlmostEqual(m.left_score(row,m.BASE),63.5)
    def test_only_component_score(self):
        row=dict(zip(m.COLS,['.5','.4','.8','1','1']));self.assertEqual(m.left_score(row,m.definitions()['only_specificity']),80)
    def test_baseline_reproduction(self):
        out=m.compare(self.rows,self.robust);m.check_outputs(out)
        base=[r for r in out[0] if r['Scenario']=='baseline_control']
        self.assertTrue(all(r['Comparison_rank']==int(r['Frozen_baseline_rank']) for r in base))
        stats=next(r for r in out[3] if r['Scenario']=='baseline_control')
        self.assertEqual(stats['Baseline_winner_changes'],0);self.assertEqual(stats['Median_ordinal_rank_spearman'],1)
        self.assertEqual(stats['Mean_top5_jaccard'],1)
    def test_fixed_ties_and_tie_count(self):
        rows=copy.deepcopy(self.rows)
        for r in rows:r[m.COLS[2]]='0.5'
        # Keep baseline control consistent with synthetic changes, without touching real inputs.
        for vid in m.u.MUTATIONS:
            g=[r for r in rows if r['VariationID']==vid]
            for r in g:r[m.SCORE]=str(m.left_score(r,m.BASE))
            g.sort(key=lambda r:(-float(r[m.SCORE]),int(r[m.COUNT]),r['Guide_ID']))
            for i,r in enumerate(g,1):r[m.RANK]=str(i)
        out=m.compare(rows,self.robust)
        for v in m.u.MUTATIONS:
            g=[r for r in rows if r['VariationID']==v]
            expected=min(g,key=lambda r:(int(r[m.COUNT]),r['Guide_ID']))['Guide_ID']
            winner=next(r for r in out[1] if r['Scenario']=='only_specificity' and r['VariationID']==v)
            self.assertEqual(winner['Guide_ID'],expected);self.assertEqual(winner['Primary_score_tie_count'],47)
    def test_deterministic_input_order(self):
        a=m.compare(self.rows,self.robust);random.Random(5).shuffle(self.rows);self.robust.reverse()
        self.assertEqual(a,m.compare(self.rows,self.robust))
    def test_duplicate_guide(self):
        self.rows[1]['Guide_ID']=self.rows[0]['Guide_ID']
        with self.assertRaisesRegex(ValueError,'Duplicate'):m.validate(self.rows,self.robust)
    def test_missing_guide(self):
        self.rows.pop()
        with self.assertRaisesRegex(ValueError,'cohort'):m.validate(self.rows,self.robust)
    def test_nonfinite(self):
        for v in ['NaN','inf','-Infinity','']:
            with self.subTest(value=v),self.assertRaises(ValueError):m.u.number(v,True)
    def test_unit_range(self):
        with self.assertRaises(ValueError):m.u.number('1.1',True)
    def test_missing_duplicate_extra_columns(self):
        for data in [b'Guide_ID\nA\n',b'Guide_ID,Guide_ID\nA,A\n']:
            with self.assertRaises(ValueError):m.u.table(data,self.schema['ranked'])
    def test_coordinate_conflict(self):
        self.rows[0]['Fullwindow_catalytic_query_position_1based']='20'
        with self.assertRaisesRegex(ValueError,'Coordinate'):m.validate(self.rows,self.robust)
    def test_accessibility_must_be_constant(self):
        self.rows[-1][m.COLS[0]]='0.123';self.rows[-1][m.SCORE]=str(m.left_score(self.rows[-1],m.BASE))
        with self.assertRaises(ValueError):m.validate(self.rows,self.robust)
    def test_incomplete_output(self):
        out=m.compare(self.rows,self.robust);out[0].pop()
        with self.assertRaisesRegex(ValueError,'Incomplete'):m.check_outputs(out)

class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);self.stage=self.root/'stage';self.stage.mkdir()
        rows,robust,schema=fixtures();paths={}
        for role,data in [('ranked',rows),('robust_top',robust)]:
            paths[role]=self.root/(role+'.csv');m.u.write_csv(paths[role],data,schema[role])
        for prefix,role,name,stage in [('ranking','ranked','ranked.csv','composite-ranking'),('robustness','robust_top','robust_top_guides.csv','robustness')]:
            paths[prefix+'_manifest']=self.root/(prefix+'_manifest.json')
            m.u.write_json(paths[prefix+'_manifest'],{'files':[{'path':name,'bytes':paths[role].stat().st_size,'sha256':m.u.digest(paths[role].read_bytes())}]})
            paths[prefix+'_result']=self.root/(prefix+'_result.json')
            m.u.write_json(paths[prefix+'_result'],{'stage':stage,'status':'COMPLETE','exit_code':0,'preservation':'PASS',
                             'output_manifest_sha256':m.u.digest(paths[prefix+'_manifest'].read_bytes())})
        paths.update(adapter=HERE/'mecp2_ablation.py',io=HERE/'ablation_io.py',schema=HERE/'schema.json')
        self.paths=paths;self.cfg={'version':m.VERSION,'scenarios':m.definitions(),'tie_rule':'score_desc_editable_count_asc_Guide_ID_asc',
                 'promotion_approved':False,'python_sha256':m.u.digest(Path(sys.executable).read_bytes()),'staging_root':str(self.stage),
                 'run_id':'SYNTHETIC_ONLY','protected_roots':[str(self.root/'official')],
                 'inputs':{r:{'path':str(p),'bytes':p.stat().st_size,'sha256':m.u.digest(p.read_bytes())} for r,p in paths.items()}}
        self.cp=self.root/'config.json';self.save();self.ap=self.root/'approval.json';self.approve()
    def save(self):self.cp.write_text(json.dumps(self.cfg));self.ch=m.u.digest(self.cp.read_bytes())
    def approve(self):
        self.ap.write_text(json.dumps({'approved':True,'approved_by':'SYNTHETIC_TEST_ONLY','stage':'component-ablation','config_sha256':self.ch,
         'python_sha256':self.cfg['python_sha256'],'run_directory':str(self.stage/self.cfg['run_id']),'promotion_approved':False,
         'expires_utc':(dt.datetime.now(dt.timezone.utc)+dt.timedelta(hours=1)).isoformat()}))
    def execute(self,stage='component-ablation'):return m.execute(stage,self.cp,self.ch,self.ap)
    def test_preflight_never_compares(self):
        with patch.object(m,'compare',side_effect=AssertionError('No comparison in preflight')):
            self.assertFalse(self.execute('verify-inputs')['comparison_calculated'])
        self.assertEqual(list(self.stage.iterdir()),[])
    def test_success_and_preservation(self):
        before={r:p.read_bytes() for r,p in self.paths.items()};result=self.execute();self.assertEqual((result['status'],result['exit_code']),('COMPLETE',0))
        folder=self.stage/self.cfg['run_id'];self.assertEqual({p.name for p in folder.iterdir()},m.OUT)
        for e in json.loads((folder/'output_manifest.json').read_text())['files']:self.assertEqual(m.u.digest((folder/e['path']).read_bytes()),e['sha256'])
        self.assertEqual(before,{r:p.read_bytes() for r,p in self.paths.items()})
    def test_missing_approval(self):
        with self.assertRaisesRegex(ValueError,'approval'):m.execute('component-ablation',self.cp,self.ch)
    def test_wrong_approval(self):
        a=json.loads(self.ap.read_text());a['config_sha256']='0'*64;self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'Approval'):self.execute()
    def test_expired_approval(self):
        a=json.loads(self.ap.read_text());a['expires_utc']='2000-01-01T00:00:00Z';self.ap.write_text(json.dumps(a))
        with self.assertRaisesRegex(ValueError,'expired'):self.execute()
    def test_stale_hash(self):
        self.paths['ranked'].write_text('STALE')
        with self.assertRaisesRegex(ValueError,'Stale input'):self.execute()
    def test_wrong_config_hash(self):
        with self.assertRaisesRegex(ValueError,'config hash'):m.execute('verify-inputs',self.cp,'0'*64)
    def test_missing_binding(self):
        self.cfg['inputs']['ranked']['path']=None;self.save()
        with self.assertRaisesRegex(ValueError,'Unassigned'):self.execute()
    def test_existing_directory(self):
        (self.stage/self.cfg['run_id']).mkdir()
        with self.assertRaisesRegex(ValueError,'already exists'):self.execute()
    def test_protected_root(self):
        self.cfg['protected_roots']=[str(self.stage)];self.save()
        with self.assertRaisesRegex(ValueError,'Prohibited'):self.execute()
    def test_wrong_interpreter(self):
        self.cfg['python_sha256']='0'*64;self.save()
        with self.assertRaisesRegex(ValueError,'Interpreter'):self.execute()
    def test_interruption(self):
        with patch.object(m,'compare',side_effect=KeyboardInterrupt('Synthetic interruption')):
            with self.assertRaises(KeyboardInterrupt):self.execute()
        r=json.loads((self.stage/self.cfg['run_id']/'run_result.json').read_text());self.assertEqual(r['exit_code'],4)
    def test_incomplete_output(self):
        with patch.object(m,'compare',return_value=([],[],[],[])):
            with self.assertRaisesRegex(ValueError,'Incomplete'):self.execute()
    def test_preservation_failure(self):
        original=m.compare
        def alter(*a,**k):
            out=original(*a,**k)
            if os.name!='nt':self.paths['ranked'].write_text('Synthetic change')
            else:raise ValueError('Preservation failure injected')
            return out
        with patch.object(m,'compare',side_effect=alter):
            with self.assertRaisesRegex(ValueError,'Preservation failure'):self.execute()
    def test_status_exit(self):
        with contextlib.redirect_stdout(io.StringIO()) as out:code=m.main(['verify-inputs','--config',str(self.cp),'--config-sha256','0'*64])
        self.assertEqual(code,json.loads(out.getvalue())['exit_code']);self.assertEqual(code,4)
    def test_other_stages_refused(self):
        for stage in ['ranking','robustness','promotion','assemble-tier1']:
            with self.assertRaisesRegex(ValueError,'prohibited'):self.execute(stage)

if __name__=='__main__':unittest.main(verbosity=2)
