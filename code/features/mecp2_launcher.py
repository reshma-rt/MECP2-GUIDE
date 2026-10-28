"""Explicit, hash-bound Tier1 launcher. No BLAST/ranking/promotion entry points.
Python 3.10+, standard library. Scientific worker runs ONLY with a separate
stage-specific approval and a matching earlier input-verification receipt.
"""
import argparse,contextlib,csv,hashlib,json,math,os,re,stat,subprocess,sys,traceback
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

PROCESSOR_SHA='43b709217564efa0d218be6a03b1ec733dd350e4369978c0dfd185a5c6e43ca1'
HELPER_SHA='220d5e62c64a037609cebdec72660bf1a64cc6717f457efcde6ae4470690369e'
INPUT_ROLES=('dataset1','candidates','structure','stage2c','queries','raw_hits','reference','processor','helper')
CODE_ROLES=('launcher','bridge','python')
RAW16='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq'.split()
MAIN_NAME='MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv'
SUMMARY_NAME='MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_summary.csv'
HITS_NAME='stage2D_v3_reconstructed_hits.tsv'
STATUS_NAME='stage2D_v4_safeguard_run_status.json'
GUARD_SUMMARY='stage2D_v4_safeguard_run_summary.csv'
MAIN_EXTRA=['Fullwindow_intended_transcript_hit_le3mm','Fullwindow_intended_best_mismatch_count','Fullwindow_non_target_hits_le3mm',*[f'Fullwindow_non_target_{i}_mismatch_hits' for i in range(4)],'Fullwindow_editable_A_offtarget_hits_le3mm',*[f'Fullwindow_editable_A_offtarget_{i}mm' for i in range(4)],'Fullwindow_GuideAware_raw_offtarget_risk','Fullwindow_GuideAware_transcriptome_risk_score_0_1','Fullwindow_GuideAware_transcriptome_specificity_score_0_1','Fullwindow_GuideAware_relative_risk_category','Fullwindow_catalytic_query_position_1based','Fullwindow_search_definition','Fullwindow_offtarget_interpretation','Stage2D_v3_status']
SUMMARY_FIELDS='VariationID Guide_count Valid_intended_fullwindow_le3mm Zero_editable_A_offtarget_guides Editable_A_offtarget_ge1_guides'.split()
HIT_FIELDS='Guide_ID VariationID Subject_accession Full_window_start_1based Full_window_end_1based Guide_length_nt Query_transcript_equivalent_RNA Subject_full_window_RNA Full_window_mismatch_count Mismatch_positions_1based Anchor_qstart Anchor_qend Anchor_sstart Anchor_send Is_intended_transcript Is_editable_A_offtarget Catalytic_query_position_1based'.split()
CODES={'VERIFIED_INPUTS':0,'COMPLETE':0,'INCOMPLETE_REFERENCE':2,'BLOCKED_DEPENDENCY':3,'FAILED_INPUT':4,'FAILED_INTENDED_QC':5,'FAILED_RUNTIME':6,'INTERRUPTED':7,'BLOCKED_APPROVAL':8,'FAILED_OUTPUT':9,'FAILED_PRESERVATION':10}
class Stop(Exception):
 def __init__(self,status,message):self.status=status;super().__init__(message)
def need(condition,message,status='FAILED_INPUT'):
 if not condition:raise Stop(status,message)
def utc():return datetime.now(timezone.utc).isoformat()
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def no_duplicate_keys(pairs):
 d={}
 for k,v in pairs:
  if k in d:raise Stop('FAILED_INPUT','Duplicate JSON key: '+k)
  d[k]=v
 return d
def read_json(p):
 with open(p,encoding='utf-8-sig') as f:return json.load(f,object_pairs_hook=no_duplicate_keys)
def new_json(p,value):
 with open(p,'x',encoding='utf-8') as f:json.dump(value,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def plain_path(value,exists=True):
 need(isinstance(value,str) and bool(value),'Unassigned path binding')
 p=Path(value);need(p.is_absolute(),'Path must be absolute: '+value)
 for a in (p,*p.parents):
  if a.exists() or a.is_symlink():
   need(not a.is_symlink() and not(getattr(a.lstat(),'st_file_attributes',0)&getattr(stat,'FILE_ATTRIBUTE_REPARSE_POINT',1024)),'Symlink/reparse point refused: '+str(a))
 if exists:need(p.is_file(),'Required file unavailable: '+str(p),'BLOCKED_DEPENDENCY')
 return p.resolve(strict=exists)
def full_table(path,delimiter=','):
 with open(path,encoding='utf-8-sig',newline='') as f:
  rd=csv.DictReader(f,delimiter=delimiter);fields=rd.fieldnames or []
  need(fields and len(fields)==len(set(fields)),'Missing/duplicate column names: '+str(path))
  rows=list(rd)
  need(all(None not in r and all(v is not None for v in r.values()) for r in rows),'Malformed/short/extra CSV fields: '+str(path))
 return fields,rows
def number(s,label):
 try:v=float(s)
 except (TypeError,ValueError):raise Stop('FAILED_INPUT','Invalid numeric '+label)
 need(math.isfinite(v),'Nonfinite '+label);return v
def integer(s,label,minimum=0):
 need(isinstance(s,str) and re.fullmatch(r'[0-9]+',s) is not None,'Explicit integer required: '+label)
 v=int(s);need(v>=minimum,'Integer below bound: '+label);return v
def accession(s):
 token=str(s).strip().split()[0] if str(s).strip() else ''
 return next((x for x in token.split('|') if x.startswith(('NM_','NR_','XM_','XR_'))),token.strip('|'))
def query_for(guide):return guide.translate(str.maketrans('AUGC','TACG'))[::-1]
def positions(r):
 guide=r['Guide_sequence_RNA_5to3'];L=integer(r['Guide_length_nt'],'Guide_length_nt',1)
 need(L in (21,25,31) and len(guide)==L and set(guide)<=set('AUGC'),'Guide sequence/length mismatch')
 p=integer(r['Target_A_position_in_guide_1based'],'guide-relative target',1)
 need(p<=L,'Guide-relative position outside guide');q=L-p+1
 if 'Target_A_position_in_target_segment_1based' in r:
  need(integer(r['Target_A_position_in_target_segment_1based'],'target-segment position',1)==q,'Conflicting guide/target-segment coordinates')
 need(guide[p-1]=='C','Designated guide base is not C')
 return L,p,q

@contextlib.contextmanager
def read_locks(paths):
 """On Windows hold read-only handles denying write/delete sharing. Else no
 OS locking claim: hashes/reparse checks still apply. Locks do not defeat admin."""
 handles=[];kernel=None
 try:
  if os.name=='nt':
   import ctypes
   from ctypes import wintypes
   kernel=ctypes.WinDLL('kernel32',use_last_error=True)
   kernel.CreateFileW.argtypes=[wintypes.LPCWSTR,wintypes.DWORD,wintypes.DWORD,wintypes.LPVOID,wintypes.DWORD,wintypes.DWORD,wintypes.HANDLE];kernel.CreateFileW.restype=wintypes.HANDLE
   kernel.CloseHandle.argtypes=[wintypes.HANDLE];kernel.CloseHandle.restype=wintypes.BOOL
   for p in sorted(set(map(str,paths))):
    handle=kernel.CreateFileW(p,0x80000000,1,None,3,0x08000000,None)
    need(handle!=ctypes.c_void_p(-1).value,'Cannot acquire Windows read-only input lock: '+p,'BLOCKED_DEPENDENCY');handles.append(handle)
  yield
 finally:
  if kernel:
   for handle in handles:kernel.CloseHandle(handle)

def bind(config_path,cfg):
 need(cfg.get('schema_version')==1,'Unsupported/missing configuration version')
 need(cfg.get('analysis_mode')=='existing_raw_anchors','Only explicitly selected existing raw anchors are supported')
 need(cfg.get('expected_guides')==376 and cfg.get('expected_mutations')==8,'Frozen cohort contract must remain376/8')
 need(cfg.get('max_mismatches')==3 and cfg.get('optional_stage3d') is False,'Unsupported scientific scope change')
 need(cfg.get('raw_format')=={'header':False,'fields':RAW16},'Explicit16-column headerless raw format required')
 out={}
 for role in INPUT_ROLES+CODE_ROLES:
  entry=cfg.get('bindings',{}).get(role,{})
  p=plain_path(entry.get('path'));expected=entry.get('sha256')
  need(isinstance(expected,str) and re.fullmatch('[0-9a-f]{64}',expected) is not None,'Unassigned/invalid SHA256: '+role)
  need(sha(p)==expected,'Stale/hash-mismatched binding: '+role,'FAILED_PRESERVATION');out[role]=p
 need(cfg['bindings']['processor']['sha256']==PROCESSOR_SHA and cfg['bindings']['helper']['sha256']==HELPER_SHA,'Unexpected processor/helper source version')
 need(out['launcher']==Path(__file__).resolve(),'Running launcher differs from configured executable source')
 need(out['bridge']==Path(__file__).with_name('mecp2_processor_bridge.py').resolve(),'Bridge must be bound adjacent to launcher')
 need(out['python']==Path(sys.executable).resolve(),'Actual Python executable differs from binding')
 need(out['processor'].parent==out['helper'].parent,'Processor/helper staging directories differ')
 stage=plain_path(cfg.get('staging_root'),exists=False);need(stage.is_dir(),'Staging directory unavailable','BLOCKED_DEPENDENCY')
 official=plain_path(cfg.get('official_project_root'),exists=False)
 need(not stage.is_relative_to(official) and not official.is_relative_to(stage),'Staging/official project roots overlap')
 # No input or executable may be a future child of a run directory.
 out['config']=plain_path(str(config_path));return out,stage

def raw_rows(path,cfg):
 with open(path,encoding='utf-8-sig',newline='') as f:
  for line,values in enumerate(csv.reader(f,delimiter='\t'),1):
   need(len(values)==16,'Raw row must have exactly16 fields at line '+str(line))
   yield dict(zip(RAW16,values))
def validate_inputs(cfg,bound):
 required=['Guide_ID','VariationID','Transcript','Guide_sequence_RNA_5to3','Target_segment_RNA_5to3','Guide_length_nt','Target_A_position_in_guide_1based']
 tables={};fieldsets={}
 for role in ['candidates','structure','stage2c']:
  fields,rows=full_table(bound[role]);need(set(required)<=set(fields),'Required canonical design column missing in '+role)
  need(len(rows)==376,'Expected376 '+role+' records')
  idx={r['Guide_ID']:r for r in rows};need(len(idx)==376 and all(idx),'Missing/duplicate guide ID in '+role)
  need(len({r['VariationID'] for r in rows})==8,'Expected8 mutations in '+role)
  tables[role]=idx;fieldsets[role]=fields
 base=tables['candidates'];need(all(set(x)==set(base) for x in tables.values()),'Guide ID sets differ across stages')
 _,d1=full_table(bound['dataset1']);need(len(d1)==141,'Expected frozen141-row Dataset1')
 contexts={r['VariationID']:r.get('Mutant_targetA_context_201nt','') for r in d1}
 need(len(contexts)==141,'Duplicate Dataset1 IDs')
 for gid,a in base.items():
  L,p,q=positions(a);target=a['Target_segment_RNA_5to3'];ctx=contexts.get(a['VariationID'],'')
  need(len(ctx)==201 and set(ctx)<=set('AUGC') and ctx[100]=='A','Missing/inconsistent target context for '+gid)
  need(len(target)==L and set(target)<=set('AUGC') and target==ctx[101-q:101-q+L],'Target segment/context mismatch for '+gid)
  for role in ['structure','stage2c']:
   r=tables[role][gid];positions(r)
   for key in required+['Target_A_position_in_target_segment_1based']:
    if key in a or key in r:need(a.get(key)==r.get(key),'Cross-stage sequence/coordinate identity conflict '+gid+'/'+key)
  st=tables['structure'][gid]
  need(st.get('Stage2B_structure_feature_status')=='VIENNARNA_FEATURES_COMPUTED','Missing/failed stored structure result')
  for key in ['Target_A_unpaired_probability','Target_A_3nt_mean_unpaired_probability','Target_A_5nt_mean_unpaired_probability','Target_A_11nt_mean_unpaired_probability','Target_context_201nt_MFE_kcal_mol','Target_context_201nt_ensemble_free_energy_kcal_mol','Guide_MFE_kcal_mol','Target_segment_MFE_kcal_mol','Guide_target_duplex_MFE_kcal_mol']:
   v=number(st.get(key),key)
   if 'probability' in key:need(0<=v<=1,'Probability outside[0,1]')
   need(st.get(key)==tables['stage2c'][gid].get(key),'Stored structure feature mismatch across B/C')
 queries={};current=None
 with open(bound['queries'],encoding='utf-8-sig') as f:
  for line in f:
   line=line.strip()
   if line.startswith('>'):
    current=line[1:];need(current and current not in queries,'Missing/duplicate query ID');queries[current]=''
   elif line:
    need(current is not None,'Malformed query FASTA');queries[current]+=line
 need(set(queries)==set(base),'Query manifest guide IDs differ')
 need(all(queries[g]==query_for(r['Guide_sequence_RNA_5to3']) for g,r in base.items()),'Query sequence not actual-guide reverse complement')
 subjects={};qseen=set();raw_count=0
 for r in raw_rows(bound['raw_hits'],cfg):
  gid=r['qseqid'];need(gid in queries,'Unknown raw query ID');qseen.add(gid);raw_count+=1
  sid=accession(r['sseqid']);need(re.fullmatch(r'(?:NM|NR|XM|XR)_[0-9]+\.[0-9]+',sid) is not None,'Invalid/versionless subject accession')
  qstart=integer(r['qstart'],'qstart',1);qend=integer(r['qend'],'qend',1);ss=integer(r['sstart'],'sstart',1);se=integer(r['send'],'send',1)
  qlen=integer(r['qlen'],'qlen',1);slen=integer(r['slen'],'slen',1)
  need(qstart<=qend<=qlen==len(queries[gid]) and max(ss,se)<=slen,'Raw alignment coordinate/length conflict')
  need(sid not in subjects or subjects[sid]==slen,'Conflicting raw subject lengths');subjects[sid]=slen
  for col in ['pident','evalue','bitscore']:number(r[col],col)
  for col in ['length','mismatch','gapopen']:integer(r[col],col)
  need(0<=float(r['pident'])<=100,'Invalid percent identity')
  aligned=r['qseq'].replace('-','');need(aligned==queries[gid][qstart-1:qend],'Raw aligned query mismatch/stale anchors')
  need(len(r['qseq'])==len(r['sseq'])==int(r['length']),'Raw aligned sequence lengths differ')
 need(raw_count and qseen==set(base),'Empty/partial raw query coverage')
 wanted=set(subjects)|{accession(r['Transcript']) for r in base.values()};found={};seen=set();current=None;length=0
 def commit():
  if current is not None:
   need(length>0,'Empty reference sequence')
   if current in wanted:found[current]=length
 with open(bound['reference'],encoding='utf-8-sig') as f:
  for line in f:
   if line.startswith('>'):
    commit();current=accession(line[1:]);need(current and current not in seen,'Duplicate/malformed normalized FASTA identifier');seen.add(current);length=0
   elif line.strip():
    need(current is not None,'Reference data before first FASTA header');sequence=line.strip().upper();need(set(sequence)<=set('ACGTURYSWKMBDHVN'),'Malformed reference alphabet');length+=len(sequence)
  commit()
 missing=sorted(wanted-set(found));need(not missing,'Required reference subjects unavailable: '+','.join(missing[:20]),'INCOMPLETE_REFERENCE')
 need(all(found[s]==L for s,L in subjects.items()),'Raw slen/reference-length mismatch')
 return {'tables':tables,'stage2c_fields':fieldsets['stage2c'],'input_qc':{'guide_count':376,'mutation_count':8,'Dataset1_count':141,'raw_rows':raw_count,'required_subjects':len(wanted),'missing_subjects':0,'query_identity':'PASS','coordinates':'PASS','structure_context_consistency':'PASS','reference_id_coverage':'PASS','historical_raw_execution_provenance':'UNRESOLVED'}}

def validate_outputs(output,qc):
 """Check completeness/schema and reconcile declared loci/counts. No E/F scoring."""
 required=[MAIN_NAME,SUMMARY_NAME,HITS_NAME,STATUS_NAME,GUARD_SUMMARY]
 need({p.name for p in output.iterdir()}==set(required),'Unexpected/missing output files','FAILED_OUTPUT')
 try:
  fields,rows=full_table(output/MAIN_NAME);expected_fields=qc['stage2c_fields']+MAIN_EXTRA
  need(fields==expected_fields,'Main output schema differs','FAILED_OUTPUT')
  base=qc['tables']['stage2c'];main={r['Guide_ID']:r for r in rows}
  need(len(rows)==len(main)==376 and set(main)==set(base),'Partial/duplicate main guide output','FAILED_OUTPUT')
  sf,summary=full_table(output/SUMMARY_NAME);need(sf==SUMMARY_FIELDS and len(summary)==8,'Summary schema/row count','FAILED_OUTPUT')
  need(len({r['VariationID'] for r in summary})==8,'Duplicate summary mutations','FAILED_OUTPUT')
  hf,hits=full_table(output/HITS_NAME,'\t');need(hf==HIT_FIELDS and hits,'Missing/reconstructed-hit schema','FAILED_OUTPUT')
  counts=defaultdict(Counter);editable=defaultdict(Counter);intended=set();loci=set()
  for hit in hits:
   gid=hit['Guide_ID'];need(gid in base,'Unknown reconstructed guide','FAILED_OUTPUT');a=base[gid];L,p,q=positions(a)
   start=integer(hit['Full_window_start_1based'],'window start',1);end=integer(hit['Full_window_end_1based'],'window end',1)
   need(hit['VariationID']==a['VariationID'] and int(hit['Guide_length_nt'])==L and end-start+1==L,'Hit identity/window mismatch','FAILED_OUTPUT')
   locus=(gid,hit['Subject_accession'],start,end);need(locus not in loci,'Duplicate reconstructed locus','FAILED_OUTPUT');loci.add(locus)
   query=query_for(a['Guide_sequence_RNA_5to3']).replace('T','U');subject=hit['Subject_full_window_RNA'];mm=integer(hit['Full_window_mismatch_count'],'mismatch')
   need(len(subject)==L and hit['Query_transcript_equivalent_RNA']==query and mm<=3 and int(hit['Catalytic_query_position_1based'])==q,'Hit sequence/register inconsistency','FAILED_OUTPUT')
   diffs=[i+1 for i,(x,y) in enumerate(zip(query,subject)) if x!=y]
   need(len(diffs)==mm and hit['Mismatch_positions_1based']==','.join(map(str,diffs)),'Output mismatch annotation conflict','FAILED_OUTPUT')
   is_intended=accession(hit['Subject_accession'])==accession(a['Transcript']);is_editable=not is_intended and subject[q-1]=='A'
   need(hit['Is_intended_transcript']==('YES' if is_intended else 'NO') and hit['Is_editable_A_offtarget']==('YES' if is_editable else 'NO'),'Intended/editable annotation conflict','FAILED_OUTPUT')
   if is_intended and mm==0:intended.add(gid)
   if not is_intended:counts[gid][mm]+=1
   if is_editable:editable[gid][mm]+=1
  need(intended==set(base),'Intended-hit output coverage incomplete','FAILED_OUTPUT')
  for gid,r in main.items():
   for key in qc['stage2c_fields']:need(r[key]==base[gid][key],'Inherited input field changed','FAILED_OUTPUT')
   need(r['Stage2D_v3_status']=='FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE' and r['Fullwindow_intended_transcript_hit_le3mm']=='YES' and r['Fullwindow_intended_best_mismatch_count']=='0','Main intended/completion status','FAILED_OUTPUT')
   need(int(r['Fullwindow_catalytic_query_position_1based'])==positions(base[gid])[2],'Main catalytic register','FAILED_OUTPUT')
   for mm in range(4):
    need(integer(r[f'Fullwindow_non_target_{mm}_mismatch_hits'],'non-target count')==counts[gid][mm],'Hit/main mismatch-count discrepancy','FAILED_OUTPUT')
    need(integer(r[f'Fullwindow_editable_A_offtarget_{mm}mm'],'editable count')==editable[gid][mm],'Hit/main editable discrepancy','FAILED_OUTPUT')
   need(int(r['Fullwindow_non_target_hits_le3mm'])==sum(counts[gid].values()) and int(r['Fullwindow_editable_A_offtarget_hits_le3mm'])==sum(editable[gid].values()),'Hit/main totals differ','FAILED_OUTPUT')
   need(number(r['Fullwindow_GuideAware_raw_offtarget_risk'],'risk')>=0,'Negative raw risk','FAILED_OUTPUT')
   risk=number(r['Fullwindow_GuideAware_transcriptome_risk_score_0_1'],'risk component');specificity=number(r['Fullwindow_GuideAware_transcriptome_specificity_score_0_1'],'specificity')
   need(0<=risk<=1 and 0<=specificity<=1 and abs(risk+specificity-1)<1e-12,'Risk/specificity range/complement','FAILED_OUTPUT')
  bymut=defaultdict(list)
  for r in rows:bymut[r['VariationID']].append(r)
  for s in summary:
   group=bymut.get(s['VariationID'],[]);zero=sum(int(r['Fullwindow_editable_A_offtarget_hits_le3mm'])==0 for r in group)
   need(int(s['Guide_count'])==int(s['Valid_intended_fullwindow_le3mm'])==len(group)>0 and int(s['Zero_editable_A_offtarget_guides'])==zero and int(s['Editable_A_offtarget_ge1_guides'])==len(group)-zero,'Summary does not reconcile','FAILED_OUTPUT')
 except (Stop,ValueError,KeyError,TypeError,OSError) as exc:
  raise Stop('FAILED_OUTPUT',str(exc))
 return {'main_rows':376,'mutation_summary_rows':8,'reconstructed_rows':len(hits),'intended_exact_guides':len(intended),'reconciled':'PASS'}

def approval_for(stage,cfg_hash,run,approval_path,approval_hash,verification_path,verification_hash):
 need(approval_path and approval_hash and verification_path and verification_hash,'Feature calculation needs separate approval and verification receipt','BLOCKED_APPROVAL')
 ap=plain_path(approval_path);vp=plain_path(verification_path)
 need(sha(ap)==approval_hash and sha(vp)==verification_hash,'Approval/verification receipt hash mismatch','BLOCKED_APPROVAL')
 a=read_json(ap);v=read_json(vp)
 need(v.get('status')=='VERIFIED_INPUTS' and v.get('stage')=='verify-inputs' and v.get('accepted') is True and v.get('config_sha256')==cfg_hash and v.get('exit_code')==0,'Earlier input verification not successful for exact config','BLOCKED_APPROVAL')
 need(a.get('approved') is True and a.get('stage')==stage and a.get('config_sha256')==cfg_hash and a.get('run_directory')==str(run) and a.get('verification_result_sha256')==verification_hash,'Approval wrong scope/run/config','BLOCKED_APPROVAL')
 try:expires=datetime.fromisoformat(a['expires_utc'].replace('Z','+00:00'))
 except (KeyError,ValueError,TypeError):raise Stop('BLOCKED_APPROVAL','Missing/malformed approval expiry')
 need(expires.tzinfo is not None and expires>datetime.now(timezone.utc),'Approval expired','BLOCKED_APPROVAL')
 need(isinstance(a.get('approver'),str) and a['approver'].strip(),'Missing named owner approval record','BLOCKED_APPROVAL')
 return [ap,vp]

def run_stage(config_path,stage,run_id,approval_path=None,approval_hash=None,verification_path=None,verification_hash=None,*,_executor=None):
 """_executor is an internal test seam. CLI has no worker/command override."""
 result={'stage':stage,'run_id':run_id,'started_utc':utc(),'dispatches':0,'accepted':False,'downstream_stages_executed':[],'production_host_enforcement_verified':False};run=None;before={};cfg_hash=None
 try:
  need(stage in ('verify-inputs','fullwindow-features'),'Unsupported stage; raw search/ranking/promotion unavailable','BLOCKED_APPROVAL')
  need(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}',run_id or '') is not None,'Invalid/unassigned run ID')
  cp=plain_path(str(Path(config_path).absolute()));cfg_hash=sha(cp);cfg=read_json(cp);result['config_sha256']=cfg_hash
  bound,staging=bind(cp,cfg);run=staging/run_id
  need(not run.exists(),'Fresh run directory already exists; no overwrite','FAILED_OUTPUT')
  need(all(not p.is_relative_to(run) for p in bound.values()),'Input path lies within output run')
  extras=[]
  if stage=='fullwindow-features':extras=approval_for(stage,cfg_hash,run,approval_path,approval_hash,verification_path,verification_hash)
  before={str(p):sha(p) for p in [*bound.values(),*extras]}
  # Exclusive mkdir reserves this run ID; never adopts an existing run.
  try:run.mkdir()
  except FileExistsError:raise Stop('FAILED_OUTPUT','Run directory was created concurrently; refused')
  result['run_directory']=str(run)
  new_json(run/'run_manifest.json',{'schema_version':1,'stage':stage,'run_id':run_id,'started_utc':result['started_utc'],'config_sha256':cfg_hash,'bound_files':before,'python_version':sys.version,'command':sys.argv,'analysis_mode':cfg['analysis_mode'],'scope':'New controlled analysis of existing anchors; not historical reproduction','windows_read_locks_requested':os.name=='nt','staging_only':True})
  with read_locks(before):
   need(all(sha(p)==hh for p,hh in before.items()),'Input changed before lock/processing','FAILED_PRESERVATION')
   qc=validate_inputs(cfg,bound);new_json(run/'input_qc.json',qc['input_qc']);result['input_qc']=qc['input_qc']
   if stage=='verify-inputs':result['status']='VERIFIED_INPUTS'
   else:
    output=run/'outputs';output.mkdir()
    request={'config_path':str(cp),'config_sha256':cfg_hash,'run_directory':str(run),'inputs':{k:str(v) for k,v in bound.items()},'expected_hashes':{k:sha(v) for k,v in bound.items()},'output_directory':str(output),'scope':'FULLWINDOW_FEATURES_ONLY'}
    new_json(run/'worker_request.json',request)
    command=[str(bound['python']),'-I','-B',str(bound['bridge']),'--request',str(run/'worker_request.json'),'--request-sha256',sha(run/'worker_request.json')]
    result['worker_command']=command;new_json(run/'worker_command.json',command)
    with open(run/'worker.stdout.log','x',encoding='utf-8') as stdout,open(run/'worker.stderr.log','x',encoding='utf-8') as stderr:
     result['dispatches']=1
     worker=(_executor or subprocess.run)(command,cwd=str(run),stdout=stdout,stderr=stderr,check=False)
    result['worker_exit_code']=worker.returncode
    need((output/STATUS_NAME).is_file() and (output/GUARD_SUMMARY).is_file(),'Worker did not produce final safeguard records','FAILED_OUTPUT')
    guard=read_json(output/STATUS_NAME);_,guard_csv=full_table(output/GUARD_SUMMARY)
    need(len(guard_csv)==1 and all(str(guard[k])==guard_csv[0][k] for k in ['run_id','status','exit_code']),'Safeguard JSON/CSV disagreement','FAILED_OUTPUT')
    need(guard.get('exit_code')==worker.returncode and CODES.get(guard.get('status'))==worker.returncode,'Safeguard status/process exit disagreement','FAILED_OUTPUT')
    if worker.returncode:
     raise Stop(guard['status'] if guard['status'] in CODES else 'FAILED_RUNTIME','Worker failed; no downstream stage accepted')
    need(guard.get('status')=='COMPLETE' and guard.get('reference_coverage')=='PASS' and guard.get('missing_subject_count')==0 and guard.get('intended_hit_qc')=='PASS','Incomplete safeguard evidence','FAILED_OUTPUT')
    result['output_qc']=validate_outputs(output,qc);result['status']='COMPLETE'
   need(all(sha(p)==hh for p,hh in before.items()),'Protected input changed during stage','FAILED_PRESERVATION')
   result['preservation']='PASS';result['accepted']=True
  # Closed logs and immutable receipts are hashed; manifest excludes itself and final result.
  expected={'run_manifest.json','input_qc.json'}
  if stage=='fullwindow-features':expected|={'worker_request.json','worker_command.json','worker.stdout.log','worker.stderr.log','resolved_processor.py','adapter_changes.json',*('outputs/'+n for n in [MAIN_NAME,SUMMARY_NAME,HITS_NAME,STATUS_NAME,GUARD_SUMMARY])}
  actual={p.relative_to(run).as_posix() for p in run.rglob('*') if p.is_file()}
  need(actual==expected,'Unexpected/missing run artifacts','FAILED_OUTPUT')
  new_json(run/'output_manifest.json',{'files':[{'path':n,'sha256':sha(run/n),'bytes':(run/n).stat().st_size} for n in sorted(actual)],'excludes':['output_manifest.json','run_result.json']})
  result['output_manifest_sha256']=sha(run/'output_manifest.json')
 except Stop as exc:result.update(status=exc.status,accepted=False,reason=str(exc))
 except KeyboardInterrupt:result.update(status='INTERRUPTED',accepted=False,reason='Interrupted')
 except Exception as exc:result.update(status='FAILED_RUNTIME',accepted=False,reason=repr(exc));traceback.print_exc()
 # Detect changes even after a failed worker; never repair files automatically.
 if before:
  changed=[]
  for p,hh in before.items():
   try:
    if not Path(p).is_file() or sha(p)!=hh:changed.append(p)
   except OSError:changed.append(p)
  if changed:result.update(status='FAILED_PRESERVATION',accepted=False,reason='Protected inputs changed',changed_paths=changed)
 result.update(finished_utc=utc(),exit_code=CODES[result['status']])
 if run is not None and run.is_dir() and result.get('run_directory')==str(run):
  try:new_json(run/'run_result.json',result)
  except Exception as exc:result.update(status='FAILED_RUNTIME',exit_code=6,accepted=False,reason='Cannot persist final run result: '+repr(exc))
 print(json.dumps(result,indent=2));return result

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('stage',choices=['verify-inputs','fullwindow-features']);ap.add_argument('--config',required=True);ap.add_argument('--run-id',required=True);ap.add_argument('--approval');ap.add_argument('--approval-sha256');ap.add_argument('--verification-result');ap.add_argument('--verification-sha256');a=ap.parse_args()
 r=run_stage(a.config,a.stage,a.run_id,a.approval,a.approval_sha256,a.verification_result,a.verification_sha256);return r['exit_code']
if __name__=='__main__':sys.exit(main())
