"""Focused synthetic integration tests. NEVER executes the scientific processor.
Feature-path tests use a fixed inert subprocess that writes synthetic artifacts.
"""
import ast,contextlib,csv,hashlib,importlib.util,io,json,os,shutil,subprocess,sys,tempfile,time,unittest,uuid
from pathlib import Path
from datetime import datetime,timezone,timedelta
import mecp2_launcher as L
import mecp2_processor_bridge as B
HERE=Path(__file__).resolve().parent
SOURCES=HERE/'test_sources'
EVIDENCE=HERE/'test_evidence'/(time.strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]) if __name__=='__main__' else None
if EVIDENCE is not None:EVIDENCE.mkdir(parents=True,exist_ok=False)
RESULTS=[]
def table(p,fields,rows,delimiter=','):
 with p.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter=delimiter);w.writeheader();w.writerows(rows)
def fixture(root):
 official=root/'official';official.mkdir();stage=root/'staging';stage.mkdir();inp=root/'inputs';inp.mkdir()
 ctx='A'*201;guide='U'*5+'C'+'U'*15;query=L.query_for(guide)
 design=[]
 for m in range(8):
  for i in range(47):design.append({'Guide_ID':f'T{m:02}{i:03}','VariationID':str(m),'Transcript':f'NM_{900000+m}.1','Guide_sequence_RNA_5to3':guide,'Target_segment_RNA_5to3':'A'*21,'Guide_length_nt':'21','Target_A_position_in_guide_1based':'6','Target_A_position_in_target_segment_1based':'16','Target_A_transcript_position':'20'})
 fields=list(design[0]);table(inp/'candidates.csv',fields,design)
 nums=['Target_A_unpaired_probability','Target_A_3nt_mean_unpaired_probability','Target_A_5nt_mean_unpaired_probability','Target_A_11nt_mean_unpaired_probability','Target_context_201nt_MFE_kcal_mol','Target_context_201nt_ensemble_free_energy_kcal_mol','Guide_MFE_kcal_mol','Target_segment_MFE_kcal_mol','Guide_target_duplex_MFE_kcal_mol']
 structures=[dict(r,**{n:('0.5' if 'probability' in n else '-1') for n in nums},Stage2B_structure_feature_status='VIENNARNA_FEATURES_COMPUTED') for r in design]
 for n in ['structure','stage2c']:table(inp/(n+'.csv'),list(structures[0]),structures)
 table(inp/'dataset1.csv',['VariationID','Mutant_targetA_context_201nt'],[{'VariationID':str(i),'Mutant_targetA_context_201nt':ctx if i<8 else ''} for i in range(141)])
 (inp/'queries.fasta').write_text(''.join('>'+r['Guide_ID']+'\n'+query+'\n' for r in design))
 (inp/'reference.fasta').write_text(''.join('>'+f'NM_{900000+m}.1'+'\n'+query+'\n' for m in range(8)))
 with (inp/'raw.tsv').open('w') as f:
  for r in design:f.write('\t'.join([r['Guide_ID'],r['Transcript'],'100','21','0','0','1','21','1','21','0','42','21','21',query,query])+'\n')
 cfg={'schema_version':1,'analysis_mode':'existing_raw_anchors','expected_guides':376,'expected_mutations':8,'max_mismatches':3,'optional_stage3d':False,'staging_root':str(stage),'official_project_root':str(official),'raw_format':{'header':False,'fields':L.RAW16},'bindings':{}}
 for role,file in [('dataset1','dataset1.csv'),('candidates','candidates.csv'),('structure','structure.csv'),('stage2c','stage2c.csv'),('queries','queries.fasta'),('raw_hits','raw.tsv'),('reference','reference.fasta')]:cfg['bindings'][role]={'path':str(inp/file),'sha256':L.sha(inp/file)}
 for role,p in [('processor',SOURCES/'processor.py'),('helper',SOURCES/'stage2d_reference_safeguard.py'),('launcher',HERE/'mecp2_launcher.py'),('bridge',HERE/'mecp2_processor_bridge.py'),('python',Path(sys.executable).resolve())]:cfg['bindings'][role]={'path':str(p),'sha256':L.sha(p)}
 p=root/'config.json';p.write_text(json.dumps(cfg,indent=2));return p,cfg

def create_inert_outputs(run):
 """Subprocess fixture output, not calculated editing/risk measurements."""
 req=L.read_json(run/'worker_request.json');cfg=L.read_json(req['config_path']);out=run/'outputs';fields,rows=L.full_table(req['inputs']['stage2c']);main=[];hits=[]
 for a in rows:
  r=dict(a);Lnt,p,q=L.positions(a);r.update({k:'0' for k in L.MAIN_EXTRA});r.update(Fullwindow_intended_transcript_hit_le3mm='YES',Fullwindow_intended_best_mismatch_count='0',Fullwindow_GuideAware_transcriptome_specificity_score_0_1='1',Fullwindow_GuideAware_relative_risk_category='LOWER_RELATIVE_RISK',Fullwindow_catalytic_query_position_1based=str(q),Fullwindow_search_definition='SYNTHETIC_FIXTURE',Fullwindow_offtarget_interpretation='SYNTHETIC_NOT_SCIENTIFIC',Stage2D_v3_status='FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE');main.append(r)
  hit={k:'' for k in L.HIT_FIELDS};hit.update(Guide_ID=a['Guide_ID'],VariationID=a['VariationID'],Subject_accession=a['Transcript'],Full_window_start_1based='1',Full_window_end_1based='21',Guide_length_nt='21',Query_transcript_equivalent_RNA=L.query_for(a['Guide_sequence_RNA_5to3']).replace('T','U'),Subject_full_window_RNA=L.query_for(a['Guide_sequence_RNA_5to3']).replace('T','U'),Full_window_mismatch_count='0',Anchor_qstart='1',Anchor_qend='21',Anchor_sstart='1',Anchor_send='21',Is_intended_transcript='YES',Is_editable_A_offtarget='NO',Catalytic_query_position_1based=str(q));hits.append(hit)
 table(out/L.MAIN_NAME,fields+L.MAIN_EXTRA,main)
 table(out/L.SUMMARY_NAME,L.SUMMARY_FIELDS,[dict(zip(L.SUMMARY_FIELDS,[str(m),'47','47','47','0'])) for m in range(8)])
 table(out/L.HITS_NAME,L.HIT_FIELDS,hits,'\t')
 guard={'run_id':'synthetic-helper','status':'COMPLETE','exit_code':0,'reference_coverage':'PASS','missing_subject_count':0,'intended_hit_qc':'PASS'}
 (out/L.STATUS_NAME).write_text(json.dumps(guard));table(out/L.GUARD_SUMMARY,list(guard),[{k:str(v) for k,v in guard.items()}])
 (run/'resolved_processor.py').write_text('# INERT FIXTURE: SCIENTIFIC PROCESSOR NOT EXECUTED\n');L.new_json(run/'adapter_changes.json',{'scope':'INERT_FIXTURE_NO_PROCESSOR_EXECUTION'})

class FocusedTests(unittest.TestCase):
 def setUp(self):
  self.root=EVIDENCE/self._testMethodName;self.root.mkdir();self.cp,self.cfg=fixture(self.root)
 def save(self):self.cp.write_text(json.dumps(self.cfg,indent=2))
 def edit(self,role,fn,rehash=True):
  entry=self.cfg['bindings'][role];p=Path(entry['path']);fn(p)
  if rehash:entry['sha256']=L.sha(p);self.save()
 def execute(self,stage='verify-inputs',name='check',**kw):
  stdout=io.StringIO();stderr=io.StringIO()
  with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):r=L.run_stage(str(self.cp),stage,name,**kw)
  (self.root/(name+'.stdout')).write_text(stdout.getvalue());(self.root/(name+'.stderr')).write_text(stderr.getvalue())
  if r.get('run_directory'):
   final=L.read_json(Path(r['run_directory'])/'run_result.json');self.assertEqual(r,final)
  self.assertEqual(r['exit_code'],L.CODES[r['status']]);return r
 def good_verification(self):
  v=self.execute();self.assertEqual(v['status'],'VERIFIED_INPUTS');return Path(v['run_directory'])/'run_result.json'
 def feature(self,fault=None):
  v=self.good_verification();target=Path(self.cfg['staging_root'])/'features';approval=self.root/'owner_approval_FIXTURE_ONLY.json'
  approval.write_text(json.dumps({'approved':True,'stage':'fullwindow-features','config_sha256':L.sha(self.cp),'run_directory':str(target),'verification_result_sha256':L.sha(v),'approver':'SYNTHETIC_TEST_AUTHORITY_NOT_REAL_APPROVAL','expires_utc':(datetime.now(timezone.utc)+timedelta(hours=1)).isoformat()}))
  def executor(command,cwd,stdout,stderr,check):
   # Actual inert subprocess with same path/log handoff, never actual bridge/main.
   child="import sys;sys.path.insert(0,sys.argv[1]);from test_launcher import create_inert_outputs;from pathlib import Path;create_inert_outputs(Path(sys.argv[2]))"
   proc=subprocess.run([sys.executable,'-I','-B','-c',child,str(HERE),cwd],stdout=stdout,stderr=stderr,check=False)
   if proc.returncode:return proc
   if fault:fault(Path(cwd))
   return subprocess.CompletedProcess(command,0)
  return self.execute('fullwindow-features','features',approval_path=str(approval),approval_hash=L.sha(approval),verification_path=str(v),verification_hash=L.sha(v),_executor=executor)
 def fail_input(self,expected):
  r=self.execute();self.assertEqual(r['status'],expected);self.assertFalse(r['accepted']);self.assertEqual(r['dispatches'],0)
 def test_valid_input_stage(self):
  r=self.execute();self.assertEqual(r['status'],'VERIFIED_INPUTS');self.assertEqual(r['input_qc']['guide_count'],376);self.assertEqual(r['dispatches'],0)
 def test_missing_raw_selection(self):self.cfg['bindings']['raw_hits']['path']=None;self.save();self.fail_input('FAILED_INPUT')
 def test_missing_raw_hash(self):self.cfg['bindings']['raw_hits']['sha256']=None;self.save();self.fail_input('FAILED_INPUT')
 def test_stale_hash(self):self.edit('raw_hits',lambda p:p.write_text(p.read_text()+'stale\n'),False);self.fail_input('FAILED_PRESERVATION')
 def test_missing_reference(self):Path(self.cfg['bindings']['reference']['path']).unlink();self.fail_input('BLOCKED_DEPENDENCY')
 def test_partial_reference(self):self.edit('reference',lambda p:p.write_text('\n'.join(p.read_text().splitlines()[:-2])+'\n'));self.fail_input('INCOMPLETE_REFERENCE')
 def test_partial_raw_queries(self):self.edit('raw_hits',lambda p:p.write_text('\n'.join(p.read_text().splitlines()[:-1])+'\n'));self.fail_input('FAILED_INPUT')
 def change_coordinate(self,value):
  def change(p):
   fields,rows=L.full_table(p);rows[0]['Target_A_position_in_target_segment_1based']=value;table(p,fields,rows)
  self.edit('candidates',change)
 def test_coordinate_conflict(self):self.change_coordinate('15');self.fail_input('FAILED_INPUT')
 def test_coordinate_missing(self):self.change_coordinate('');self.fail_input('FAILED_INPUT')
 def test_coordinate_fractional(self):self.change_coordinate('16.0');self.fail_input('FAILED_INPUT')
 def test_coordinate_out_of_bounds(self):self.change_coordinate('22');self.fail_input('FAILED_INPUT')
 def test_duplicate_coordinate_header(self):
  def change(p):
   rows=list(csv.reader(p.open()));rows[0].append('Target_A_position_in_guide_1based');[r.append('6') for r in rows[1:]]
   with p.open('w',newline='') as f:csv.writer(f).writerows(rows)
  self.edit('candidates',change);self.fail_input('FAILED_INPUT')
 def test_empty_candidate_table(self):self.edit('candidates',lambda p:p.write_text(p.read_text().splitlines()[0]+'\n'));self.fail_input('FAILED_INPUT')
 def test_unknown_raw_query(self):self.edit('raw_hits',lambda p:p.write_text(p.read_text().replace('T00000','UNKNOWN',1)));self.fail_input('FAILED_INPUT')
 def test_raw_query_sequence_mismatch(self):self.edit('queries',lambda p:p.write_text(p.read_text().replace('G','C',1)));self.fail_input('FAILED_INPUT')
 def test_existing_run_is_untouched(self):
  run=Path(self.cfg['staging_root'])/'check';run.mkdir();sentinel=run/'keep';sentinel.write_text('unchanged');self.fail_input('FAILED_OUTPUT');self.assertEqual(sentinel.read_text(),'unchanged');self.assertEqual({x.name for x in run.iterdir()},{'keep'})
 def test_prohibited_output_overlap(self):self.cfg['staging_root']=self.cfg['official_project_root'];self.save();self.fail_input('FAILED_INPUT')
 def test_symlink_output_parent(self):
  real=Path(self.cfg['staging_root']);link=self.root/'link';link.symlink_to(real,target_is_directory=True);self.cfg['staging_root']=str(link);self.save();self.fail_input('FAILED_INPUT')
 def test_forbidden_stage(self):
  for stage in ['raw-search','ranking','robustness','promote']:
   r=self.execute(stage,stage);self.assertEqual(r['status'],'BLOCKED_APPROVAL');self.assertEqual(r['dispatches'],0)
 def test_missing_feature_approval(self):
  r=self.execute('fullwindow-features');self.assertEqual(r['status'],'BLOCKED_APPROVAL');self.assertEqual(r['dispatches'],0)
 def test_complete_synthetic_feature_handoff(self):
  r=self.feature();self.assertEqual(r['status'],'COMPLETE');self.assertTrue(r['accepted']);self.assertEqual(r['output_qc']['main_rows'],376)
 def test_partial_main_output(self):
  def fault(run):
   p=run/'outputs'/L.MAIN_NAME;fields,rows=L.full_table(p);table(p,fields,rows[:-1])
  self.assertEqual(self.feature(fault)['status'],'FAILED_OUTPUT')
 def test_empty_main_output_false_success(self):
  def fault(run):
   p=run/'outputs'/L.MAIN_NAME;fields,_=L.full_table(p);table(p,fields,[])
  self.assertEqual(self.feature(fault)['status'],'FAILED_OUTPUT')
 def test_missing_output(self):self.assertEqual(self.feature(lambda r:(r/'outputs'/L.HITS_NAME).unlink())['status'],'FAILED_OUTPUT')
 def test_bad_output_schema(self):self.assertEqual(self.feature(lambda r:(r/'outputs'/L.MAIN_NAME).write_text('wrong\n'))['status'],'FAILED_OUTPUT')
 def test_status_exit_disagreement(self):
  def fault(run):
   p=run/'outputs'/L.STATUS_NAME;g=L.read_json(p);g['exit_code']=2;p.write_text(json.dumps(g))
  self.assertEqual(self.feature(fault)['status'],'FAILED_OUTPUT')
 def test_worker_failure(self):
  def fault(run):
   p=run/'outputs'/L.STATUS_NAME;g=L.read_json(p);g.update(status='FAILED_RUNTIME',exit_code=6);p.write_text(json.dumps(g));table(run/'outputs'/L.GUARD_SUMMARY,list(g),[{k:str(v) for k,v in g.items()}])
  # Nonzero actual inert subprocess at parent handoff.
  v=self.good_verification();a=self.root/'approval.json';a.write_text(json.dumps({'approved':True,'stage':'fullwindow-features','config_sha256':L.sha(self.cp),'run_directory':str(Path(self.cfg['staging_root'])/'features'),'verification_result_sha256':L.sha(v),'approver':'FIXTURE','expires_utc':(datetime.now(timezone.utc)+timedelta(hours=1)).isoformat()}))
  def worker(cmd,cwd,stdout,stderr,check):
   run=Path(cwd);create_inert_outputs(run);fault(run)
   return subprocess.run([sys.executable,'-I','-B','-c','raise SystemExit(6)'],cwd=cwd,stdout=stdout,stderr=stderr,check=False)
  r=self.execute('fullwindow-features','features',approval_path=str(a),approval_hash=L.sha(a),verification_path=str(v),verification_hash=L.sha(v),_executor=worker);self.assertEqual(r['status'],'FAILED_RUNTIME');self.assertEqual(r['worker_exit_code'],6)
 def test_missing_guard_status(self):self.assertEqual(self.feature(lambda r:(r/'outputs'/L.STATUS_NAME).unlink())['status'],'FAILED_OUTPUT')
 def test_preservation_failure(self):
  r=self.feature(lambda run:Path(self.cfg['bindings']['stage2c']['path']).write_text('SYNTHETIC_TAMPER'));self.assertEqual(r['status'],'FAILED_PRESERVATION');self.assertFalse(r['accepted'])
 def test_unexpected_run_file(self):self.assertEqual(self.feature(lambda run:(run/'unexpected').write_text('fixture'))['status'],'FAILED_OUTPUT')
 def test_adapter_scientific_body_unchanged(self):
  source=(SOURCES/'processor.py').read_text();request={'run_directory':str(self.root/'run'),'output_directory':str(self.root/'run/outputs'),'inputs':{k:v['path'] for k,v in self.cfg['bindings'].items()}}
  adapted,changes=B.adapt_source(source,request,self.cfg);compile(adapted,'ADAPTER_COMPILE_ONLY','exec')
  changed={x['name'] for x in changes['changes']}
  def scientific_tree(text):
   tr=ast.parse(text);tr.body=[n for n in tr.body if not(isinstance(n,ast.FunctionDef) and n.name in changed) and not(isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in changed)];return ast.dump(tr,include_attributes=False)
  self.assertEqual(scientific_tree(source),scientific_tree(adapted));self.assertFalse(changes['scientific_formula_or_threshold_changes'])
 def test_actual_helper_rejects_partial_reference(self):
  spec=importlib.util.spec_from_file_location('guard_fixture_only',SOURCES/'stage2d_reference_safeguard.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
  def work():
   subjects=g.validate_subject_manifest(['NM_1.1','NM_2.1'],lambda s:s);g.require_reference_coverage(subjects,{'NM_1.1':'A'});self.fail('Partial coverage passed')
  code=g.run_guarded(work,self.root/'guard.json',self.root/'guard.csv','SYNTHETIC_REFERENCE_GUARD_ONLY');self.assertEqual(code,2);self.assertEqual(L.read_json(self.root/'guard.json')['status'],'INCOMPLETE_REFERENCE')
 def test_actual_audit_write_boundary(self):
  run=self.root/'run';run.mkdir();protected=self.root/'protected.txt';protected.write_text('PRESERVE')
  script="import sys;sys.path.insert(0,sys.argv[1]);from pathlib import Path;import mecp2_processor_bridge as b;b.install_write_boundary(Path(sys.argv[2]));\ntry:\n Path(sys.argv[3]).write_text('BAD')\nexcept Exception:\n print('WRITE_BLOCKED')\nelse:\n raise SystemExit(9)\n"
  proc=subprocess.run([sys.executable,'-I','-B','-c',script,str(HERE),str(run),str(protected)],capture_output=True,text=True);(self.root/'write_boundary.stdout').write_text(proc.stdout);(self.root/'write_boundary.stderr').write_text(proc.stderr);self.assertEqual(proc.returncode,0);self.assertIn('WRITE_BLOCKED',proc.stdout);self.assertEqual(protected.read_text(),'PRESERVE')
 def test_partial_reference_worker_status_propagates(self):
  def fault(run):
   p=run/'outputs'/L.STATUS_NAME;g=L.read_json(p);g.update(status='INCOMPLETE_REFERENCE',exit_code=2,reference_coverage='INCOMPLETE',missing_subject_count=1,intended_hit_qc='NOT_RUN');p.write_text(json.dumps(g));table(run/'outputs'/L.GUARD_SUMMARY,list(g),[{k:str(v) for k,v in g.items()}])
  v=self.good_verification();a=self.root/'approval.json';a.write_text(json.dumps({'approved':True,'stage':'fullwindow-features','config_sha256':L.sha(self.cp),'run_directory':str(Path(self.cfg['staging_root'])/'features'),'verification_result_sha256':L.sha(v),'approver':'FIXTURE','expires_utc':(datetime.now(timezone.utc)+timedelta(hours=1)).isoformat()}))
  def worker(cmd,cwd,stdout,stderr,check):
   run=Path(cwd);create_inert_outputs(run);fault(run)
   return subprocess.run([sys.executable,'-I','-B','-c','raise SystemExit(2)'],cwd=cwd,stdout=stdout,stderr=stderr,check=False)
  r=self.execute('fullwindow-features','features',approval_path=str(a),approval_hash=L.sha(a),verification_path=str(v),verification_hash=L.sha(v),_executor=worker);self.assertEqual(r['status'],'INCOMPLETE_REFERENCE');self.assertEqual(r['exit_code'],2);self.assertFalse(r['accepted'])
 def test_invalid_approval_scope_and_expiry(self):
  v=self.good_verification()
  for name,delta in [('wrong_scope',{'stage':'ranking'}),('expired',{'expires_utc':'2000-01-01T00:00:00Z'})]:
   a=self.root/(name+'.json');record={'approved':True,'stage':'fullwindow-features','config_sha256':L.sha(self.cp),'run_directory':str(Path(self.cfg['staging_root'])/name),'verification_result_sha256':L.sha(v),'approver':'FIXTURE','expires_utc':(datetime.now(timezone.utc)+timedelta(hours=1)).isoformat()};record.update(delta);a.write_text(json.dumps(record))
   r=self.execute('fullwindow-features',name,approval_path=str(a),approval_hash=L.sha(a),verification_path=str(v),verification_hash=L.sha(v));self.assertEqual(r['status'],'BLOCKED_APPROVAL');self.assertEqual(r['dispatches'],0)
 def test_guide_lengths_and_coordinate_bounds(self):
  for length in [21,25,31]:
   for pos in [1,6,length-5,length]:
    r={'Guide_sequence_RNA_5to3':'U'*(pos-1)+'C'+'U'*(length-pos),'Guide_length_nt':str(length),'Target_A_position_in_guide_1based':str(pos),'Target_A_position_in_target_segment_1based':str(length-pos+1)}
    self.assertEqual(L.positions(r),(length,pos,length-pos+1))
 def test_duplicate_json_keys(self):
  self.cp.write_text('{"schema_version":1,"schema_version":1}')
  self.fail_input('FAILED_INPUT')
 def test_audit_boundary_forbids_nested_subprocess(self):
  run=self.root/'run';run.mkdir()
  script="import sys;sys.path.insert(0,sys.argv[1]);from pathlib import Path;import mecp2_processor_bridge as b;import subprocess;b.install_write_boundary(Path(sys.argv[2]));\ntry:\n subprocess.run([sys.executable,'-c','pass'])\nexcept PermissionError:\n print('EXEC_BLOCKED')\nelse:\n raise SystemExit(9)\n"
  proc=subprocess.run([sys.executable,'-I','-B','-c',script,str(HERE),str(run)],capture_output=True,text=True);self.assertEqual(proc.returncode,0,proc.stderr);self.assertIn('EXEC_BLOCKED',proc.stdout)
 def test_existing_output_is_not_overwritten_by_worker(self):
  run=self.root/'run';run.mkdir();(run/'outputs').mkdir();p=run/'outputs'/L.MAIN_NAME;p.write_text('PRESERVE')
  script="import sys;sys.path.insert(0,sys.argv[1]);from pathlib import Path;import mecp2_processor_bridge as b;b.install_write_boundary(Path(sys.argv[2]));\ntry:\n Path(sys.argv[3]).write_text('BAD')\nexcept PermissionError:\n print('OVERWRITE_BLOCKED')\nelse:\n raise SystemExit(9)\n"
  proc=subprocess.run([sys.executable,'-I','-B','-c',script,str(HERE),str(run),str(p)],capture_output=True,text=True);self.assertEqual(proc.returncode,0,proc.stderr);self.assertIn('OVERWRITE_BLOCKED',proc.stdout);self.assertEqual(p.read_text(),'PRESERVE')
 def test_actual_helper_persistence_under_write_boundary(self):
  # Exercises status RUNNING->failure replacement under the actual hook;
  # does not import or execute the processor scientific module.
  run=self.root/'run';run.mkdir();out=run/'outputs';out.mkdir()
  script="import sys;sys.path.insert(0,sys.argv[1]);from pathlib import Path;import importlib.util;import mecp2_processor_bridge as b;import mecp2_launcher as L;s=importlib.util.spec_from_file_location('fixture_guard',sys.argv[3]);g=importlib.util.module_from_spec(s);s.loader.exec_module(g);b.install_write_boundary(Path(sys.argv[2]));\ndef work():\n subs=g.validate_subject_manifest(['NM_1.1','NM_2.1'],lambda x:x);g.require_reference_coverage(subs,{'NM_1.1':'A'})\nraise SystemExit(g.run_guarded(work,Path(sys.argv[2])/'outputs'/L.STATUS_NAME,Path(sys.argv[2])/'outputs'/L.GUARD_SUMMARY,'SYNTHETIC_ONLY'))"
  proc=subprocess.run([sys.executable,'-I','-B','-c',script,str(HERE),str(run),str(SOURCES/'stage2d_reference_safeguard.py')],capture_output=True,text=True);(self.root/'guard_hook.stdout').write_text(proc.stdout);(self.root/'guard_hook.stderr').write_text(proc.stderr);self.assertEqual(proc.returncode,2,proc.stderr);self.assertEqual(L.read_json(out/L.STATUS_NAME)['status'],'INCOMPLETE_REFERENCE')
 def test_guard_inside_run_write_allowed(self):
  run=self.root/'run';run.mkdir();(run/'outputs').mkdir();p=run/'outputs'/L.MAIN_NAME
  script="import sys;sys.path.insert(0,sys.argv[1]);from pathlib import Path;import mecp2_processor_bridge as b;b.install_write_boundary(Path(sys.argv[2]));Path(sys.argv[3]).write_text('INERT')"
  proc=subprocess.run([sys.executable,'-I','-B','-c',script,str(HERE),str(run),str(p)],capture_output=True,text=True);self.assertEqual(proc.returncode,0,proc.stderr);self.assertEqual(p.read_text(),'INERT')

class RecordingResult(unittest.TextTestResult):
 def addSuccess(self,test):super().addSuccess(test);RESULTS.append({'test':test.id().split('.')[-1],'status':'PASS'})
 def addFailure(self,test,err):super().addFailure(test,err);RESULTS.append({'test':test.id().split('.')[-1],'status':'FAIL'})
 def addError(self,test,err):super().addError(test,err);RESULTS.append({'test':test.id().split('.')[-1],'status':'ERROR'})
if __name__=='__main__':
 runner=unittest.TextTestRunner(verbosity=2,resultclass=RecordingResult);r=runner.run(unittest.defaultTestLoader.loadTestsFromTestCase(FocusedTests))
 table(HERE/'TEST_RESULTS.csv',['test','status'],RESULTS)
 (HERE/'TEST_SUMMARY.json').write_text(json.dumps({'tests_run':r.testsRun,'passed':sum(x['status']=='PASS' for x in RESULTS),'failures':len(r.failures),'errors':len(r.errors),'evidence_directory':str(EVIDENCE),'scientific_processor_executions':0,'BLAST_executions':0,'Windows_host_enforcement_tested':False},indent=2)+'\n');sys.exit(0 if r.wasSuccessful() else 1)
