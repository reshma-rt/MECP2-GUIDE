"""Artificial test records only; never load scientific input files here."""
from pathlib import Path
import json
import robustness_core as core
import robustness_io as util

def fixtures(constant=False):
    schema=json.loads((Path(__file__).parent/'schema_contract.json').read_text())
    rows=[]
    for i in range(376):
        r={c:'' for c in schema['ranked']}
        r.update(Guide_ID='T1G%04d'%(i+1),VariationID=sorted(util.MUTATIONS)[i//47],Transcript='SYNTHETIC_ONLY',
                 Guide_length_nt='21',Target_A_transcript_position='1000',Target_A_position_in_guide_1based='6',
                 Target_A_position_in_target_segment_1based='16',Fullwindow_catalytic_query_position_1based='16',
                 Guide_sequence_RNA_5to3='UUUUUCUUUUUUUUUUUUUUU',Target_segment_RNA_5to3='A'*21,
                 Stage2E_v2_status='HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE',
                 Stage2D_v3_status='FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE',
                 Target_A_unpaired_probability=str((i%9)/10),Target_A_5nt_mean_unpaired_probability='.4',
                 Bystander_precision_score_0_1='.7',Fullwindow_GuideAware_transcriptome_specificity_score_0_1=str((i%11)/11),
                 Guide_target_duplex_MFE_kcal_mol=str(-20-i%7),Guide_MFE_kcal_mol=str(-8+i%5),
                 Fullwindow_editable_A_offtarget_hits_le3mm=str(i%3))
        comps=[((i%9)/10+.4)/2,.7,(i%11)/11,1-((-20-i%7)+26)/6,((-8+i%5)+8)/4]
        if constant:
            r['Guide_target_duplex_MFE_kcal_mol']='-5';r['Guide_MFE_kcal_mol']='-5'
            comps[3:]=[0.0,1.0]
        for c,x in zip(core.COMPONENTS,comps):r[c]=str(x)
        r[core.BASELINE_SCORE]=str(sum(w*x for w,x in zip(core.SCENARIOS['baseline'],comps))*100)
        for c,w in zip(['accessibility','bystander','specificity','duplex','self_structure'],core.SCENARIOS['baseline']):
            r['Stage2E_v2_weight_'+c]=str(w)
        rows.append(r)
    top=[]
    for v in sorted(util.MUTATIONS):
        g=sorted([r for r in rows if r['VariationID']==v],key=lambda r:(-float(r[core.BASELINE_SCORE]),int(r['Fullwindow_editable_A_offtarget_hits_le3mm']),r['Guide_ID']))
        for j,r in enumerate(g,1):r[core.BASELINE_RANK]=str(j)
        top.append(g[0])
    rows.sort(key=lambda r:(r['VariationID'],int(r[core.BASELINE_RANK]),r['Guide_ID']))
    summary=[{'VariationID':r['VariationID'],'Stage2E_v2_top_Guide_ID':r['Guide_ID'],
              'Stage2E_v2_top_score_0_100':r[core.BASELINE_SCORE],
              'Editable_A_offtarget_hits_le3mm':r['Fullwindow_editable_A_offtarget_hits_le3mm']} for r in top]
    norm={'duplex_MFE_range':[-26,-20],'guide_MFE_range':[-8,-4],'cohort_size':376,
          'constant_duplex_component':0.,'constant_self_component':1.}
    if constant:
        norm['duplex_MFE_range']=[-5,-5];norm['guide_MFE_range']=[-5,-5]
    return rows,top,summary,norm,schema
