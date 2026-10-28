"""Fixed adapter for the exact staged Stage2D source; no standalone BLAST.
Only invoked by the launcher for a separately approved fullwindow-features run.
Synthetic tests compile/inspect this adapter but never execute the processor.
"""
import argparse,ast,hashlib,importlib.util,json,os,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import mecp2_launcher as L

def adapt_source(source,request,cfg):
 """Change I/O bindings and explicit raw schema only; scientific body retained."""
 L.need(hashlib.sha256(source.encode()).hexdigest()==L.PROCESSOR_SHA,'Processor source version mismatch')
 tree=ast.parse(source);output=Path(request['output_directory']);inputs=request['inputs']
 assignments={'ROOT':str(Path(request['run_directory'])),'PIPE':str(output),'REF':str(Path(inputs['reference']).parent),'INPUT':inputs['stage2c'],'FASTA':inputs['reference'],'OUT_MAIN':str(output/L.MAIN_NAME),'OUT_SUM':str(output/L.SUMMARY_NAME),'OUT_HITS':str(output/L.HITS_NAME)}
 changed=[];counts={n:0 for n in [*assignments,'KNOWN_GUIDE_COLUMNS','detect_raw_hits','parse_blast_hits']}
 for i,node in enumerate(tree.body):
  if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
   name=node.targets[0].id
   if name in assignments:
    node.value=ast.Call(func=ast.Name(id='Path',ctx=ast.Load()),args=[ast.Constant(assignments[name])],keywords=[]);counts[name]+=1;changed.append({'name':name,'original_line':node.lineno,'operation':'explicit_path_binding'})
   elif name=='KNOWN_GUIDE_COLUMNS':
    node.value=ast.List(elts=[ast.Constant('Guide_sequence_RNA_5to3')],ctx=ast.Load());counts[name]+=1;changed.append({'name':name,'original_line':node.lineno,'operation':'canonical_validated_sequence_column'})
  elif isinstance(node,ast.FunctionDef) and node.name=='detect_raw_hits':
   replacement=ast.parse('def detect_raw_hits():\n    return Path('+repr(inputs['raw_hits'])+')').body[0];tree.body[i]=replacement;counts[node.name]+=1;changed.append({'name':node.name,'original_line':node.lineno,'operation':'remove_all_discovery_return_explicit_raw_path'})
  elif isinstance(node,ast.FunctionDef) and node.name=='parse_blast_hits':
   replacement=ast.parse('def parse_blast_hits(path):\n    if Path(path).resolve() != Path('+repr(inputs['raw_hits'])+').resolve():\n        raise ValueError("Unexpected raw input path")\n    return list(_explicit_raw_rows(path, _explicit_config)), list(_explicit_raw_fields)').body[0];tree.body[i]=replacement;counts[node.name]+=1;changed.append({'name':node.name,'original_line':node.lineno,'operation':'declared16_column_schema_no_silent_row_skipping'})
 L.need(all(v==1 for v in counts.values()),'Source interface changed; adapter refuses unexpected AST')
 tree=ast.fix_missing_locations(tree);text=ast.unparse(tree)+'\n'
 return text,{'processor_original_sha256':L.PROCESSOR_SHA,'adapted_sha256':hashlib.sha256(text.encode()).hexdigest(),'changes':changed,'scientific_formula_or_threshold_changes':False}

def install_write_boundary(run):
 """Whitelist worker outputs; not a sandbox against native code/admin bypass."""
 import re
 run=Path(run).resolve();output=run/'outputs';owned=set()
 regular={output/n for n in [L.MAIN_NAME,L.SUMMARY_NAME,L.HITS_NAME,L.STATUS_NAME,L.GUARD_SUMMARY]}
 status_dest={output/L.STATUS_NAME,output/L.GUARD_SUMMARY}
 def checked(value):
  if isinstance(value,int):raise PermissionError('Raw descriptor write-open unavailable')
  p=Path(os.fsdecode(value))
  if not p.is_absolute():p=Path.cwd()/p
  L.need(p.resolve().is_relative_to(run),'Prohibited write outside this fresh run: '+str(p),'FAILED_PRESERVATION')
  for ancestor in (p,*p.parents):
   if ancestor.exists() or ancestor.is_symlink():
    L.need(not ancestor.is_symlink() and not(getattr(ancestor.lstat(),'st_file_attributes',0)&1024),'Reparse write path refused','FAILED_PRESERVATION')
  return p.resolve()
 def allowed_file(p):
  if p in regular:return True
  return p.parent==output and any(re.fullmatch(re.escape(name)+r'\.[0-9a-f-]{36}\.tmp',p.name) for name in [L.STATUS_NAME,L.GUARD_SUMMARY])
 def audit(event,args):
  if event=='open':
   path,mode,flags=args
   if (isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)):
    p=checked(path)
    if not allowed_file(p) or p.exists():raise PermissionError('Unexpected/existing worker output refused: '+str(p))
    owned.add(p)
  elif event=='os.mkdir':
   p=checked(args[0])
   if p!=output or not p.is_dir():raise PermissionError('Worker directory creation unavailable')
  elif event in ('os.rename','os.replace'):
   if len(args)>2 and any(x not in (-1,None) for x in args[2:]):raise PermissionError('dir_fd indirection unavailable')
   src,dst=checked(args[0]),checked(args[1])
   if src not in owned or src.suffix!='.tmp' or dst not in status_dest or (dst.exists() and dst not in owned):raise PermissionError('Unowned/non-status replacement refused')
   owned.add(dst)
  elif event in ('os.remove','os.unlink','os.rmdir','os.chmod','os.utime','os.truncate','os.link','os.symlink','subprocess.Popen','os.system','os.exec','os.posix_spawn','socket.connect','socket.bind') or event.startswith('os.exec'):
   raise PermissionError('Forbidden worker operation: '+event)
 sys.addaudithook(audit)

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--request',required=True);ap.add_argument('--request-sha256',required=True);a=ap.parse_args()
 L.need(L.sha(a.request)==a.request_sha256,'Worker request hash mismatch','FAILED_PRESERVATION')
 request=L.read_json(a.request);L.need(request.get('scope')=='FULLWINDOW_FEATURES_ONLY','Invalid worker stage','BLOCKED_APPROVAL')
 run=Path(request['run_directory']).resolve();output=Path(request['output_directory']).resolve()
 L.need(Path(a.request).resolve()==run/'worker_request.json' and Path.cwd().resolve()==run and output==run/'outputs','Worker run paths differ')
 L.need(output.is_dir() and not any(output.iterdir()),'Worker output directory is not empty','FAILED_OUTPUT')
 cp=Path(request['config_path']);L.need(L.sha(cp)==request['config_sha256'],'Worker config stale','FAILED_PRESERVATION')
 cfg=L.read_json(cp);bound,_=L.bind(cp,cfg)
 L.need(all(str(bound[k])==request['inputs'][k] and L.sha(bound[k])==request['expected_hashes'][k] for k in bound),'Worker request/input binding mismatch','FAILED_PRESERVATION')
 text,changes=adapt_source(bound['processor'].read_text(encoding='utf-8'),request,cfg)
 # These are new run artifacts only; original/staged sources remain untouched.
 with open(run/'resolved_processor.py','x',encoding='utf-8') as f:f.write(text)
 L.new_json(run/'adapter_changes.json',changes)
 spec=importlib.util.spec_from_file_location('stage2d_reference_safeguard',bound['helper']);helper=importlib.util.module_from_spec(spec);sys.modules[spec.name]=helper;spec.loader.exec_module(helper)
 install_write_boundary(run)
 def work():
  env={'__name__':helper.WORKER_NAME,'__file__':str(run/'resolved_processor.py'),'_explicit_raw_rows':L.raw_rows,'_explicit_config':cfg,'_explicit_raw_fields':L.RAW16}
  exec(compile(text,str(run/'resolved_processor.py'),'exec'),env)
 return helper.run_guarded(work,output/L.STATUS_NAME,output/L.GUARD_SUMMARY,'TIER1_EXPLICIT_INPUT_FULLWINDOW_FEATURES')
if __name__=='__main__':
 try:sys.exit(main())
 except L.Stop as exc:print(json.dumps({'status':exc.status,'reason':str(exc)}),file=sys.stderr);sys.exit(L.CODES[exc.status])
