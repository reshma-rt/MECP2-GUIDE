# MECP2 Phase3F completed feature outputs — verification result

25 October 2026, Asia/Calcutta. **Decision: the uploaded feature outputs pass the inspected integrity, identity, schema, annotation, arithmetic and raw-anchor reconciliation checks and are suitable as inputs to a separately approved ranking proposal.** No ranking, robustness, promotion, BLAST or production processor was run during this verification.

## Access and evidence boundary

I inspected the four uploaded copies in the current workspace. I did not open the user's Windows drive directly. The exact uploaded output_manifest.json SHA256 is `48bd98e1188b0d3f3cb340282a6fec641804d643227544d16cb96cb519580249`, matching the launcher result posted in the conversation. All three uploaded table bytes/sizes match their canonical manifest entries; `(1)` is only the upload filename suffix and does not indicate a data change.

The manifest lists13 artifacts. Three listed table artifacts were supplied and independently checked; the manifest itself was checked against the reported final result. The other10 listed run artifacts (including helper status/summary,actual generated processor,request and runtime logs) were not supplied as current-run files and cannot be independently declared inspected. Their hashes in the manifest and the posted COMPLETE/exit0 receipt remain evidence reported by the completed run. No new blocker is invented for preparing a ranking proposal, but this scope must remain explicit. The actual current Windows paths/state can be checked with the optional read-only command below.

| Uploaded filename | Rows | Columns | Bytes | Independently computed SHA256 |
| --- | ---: | ---: | ---: | --- |
| `MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget(1).csv` | 376 | 76 | 633552 | `36683b5b89853e8a49c044e4e78376ee5034c1dea219120b6e38da2074e585ab` |
| `MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_summary(1).csv` | 8 | 5 | 280 | `4675a61f80327e67502deb17512f4cfb915da2ff307a5d1908d9b04f6eab4724` |
| `stage2D_v3_reconstructed_hits.tsv` | 6782 | 17 | 840761 | `b02fec5cb7b3286b0634012193e26c103289c5721482fba625051ce72d630d1d` |

## Identity, completeness and consistency

- Exact guide-ID set T1G0001 through T1G0376 matches the recovered frozen candidate,B and C tables. No extra/missing/duplicate guide records. Eight unique VariationIDs,47 guides per variant; all exact-restoration annotations retained.
- Main schema76 columns =56 inherited Stage2C columns plus20 full-window fields; summary5 columns; reconstructed-hit table17 columns. No duplicate headers or short/extra-field rows. All **21,056 inherited Stage2C cells** are unchanged.
- The recovered B/C source tables were verified against their pinned original hashes. Raw B structure values match; two expected stage-metadata differences between B and C are Dataset2_stage and Bystander_risk_model_status. An initial over-broad checker treated one as a failure; it was corrected to compare the actual structural fields and retained C lineage, not to alter any input. These were checker assumptions, not repaired output data.
- All376 catalytic query positions satisfy q=L−guide_position+1 and agree with explicit target-segment positions. Designated guide bases are C. No absolute transcript coordinate was substituted. All guide/query identities,sequence lengths and mutation assignments agree.
- All6,782 retained loci have unique `(Guide_ID,subject accession,start,end)` keys, correct anchor-derived starts/inclusive lengths,Hamming mismatch counts/position lists,at most3 mismatches,intended-accession flags and editable-A annotations at the corrected register.
- Exactly376 intended rows (one per guide),all zero mismatches;6,406 non-target rows. All per-guide mismatch-bin,total and editable counts independently reconcile to stored hits. Current flags total27 rows across8 guides;368 guides have no detected editable-A row within this search definition.
- Raw-risk arithmetic checked against the supplied windows: sum of editable counts/(1+mm) plus0.25N0+0.15N1+0.08N2+0.04N3. Stored raw-risk range0.75–7.72 and min–max risk/specificity values reconcile to floating-point tolerance. This is verification of supplied values, not generation of replacement feature tables. All required probabilities/specificity/precision values are finite and in[0,1].

| VariationID | Guides | Exact intended guides | Zero detected editable-A guides | Guides with detected editable-A |
| --- | ---: | ---: | ---: | ---: |
| 143546 | 47 | 47 | 47 | 0 |
| 143559 | 47 | 47 | 47 | 0 |
| 143597 | 47 | 47 | 47 | 0 |
| 2505626 | 47 | 47 | 47 | 0 |
| 3237211 | 47 | 47 | 46 | 1 |
| 372650 | 47 | 47 | 43 | 4 |
| 4279604 | 47 | 47 | 47 | 0 |
| 548706 | 47 | 47 | 44 | 3 |

## Exact reconciliation of6,782 versus6,359

Locus key is `(Guide_ID,normalized subject accession,full_start,full_end)`. Projection follows the documented processor: plus subject orientation; `full_start=sstart−qstart+1`, `full_end=full_start+guide_length−1`; deduplicate that key,retaining the longest query-span anchor and the first on equal span. Retained output is then filtered to full-window Hamming distance≤3. The register is used afterward for editable-A annotation; it does not select window coordinates or mismatch eligibility.

| Set relationship | Loci |
| --- | ---: |
| Archived windows |6,359|
| New windows |6,782|
| Shared |6,344|
| New only |438|
| Archive only |15|
| Net difference |438−15=423|

**Every addition/removal is explained by the actual raw-anchor sets.** None of the438 new-only loci projects from the older raw file; all438 project from the selected guide-aware file. All15 archive-only loci project from the older raw file; none projects from the selected guide-aware file. All6,344 shared loci project from both. The438 additions affect51 guides (408 at3mm,30 at2mm); all15 removals are3mm sites affecting6 guides. All additions are non-target rows with current editable flagNO. Fourteen removed rows had historical editable flagYES; their absence under another anchor set is not proof of biological false positivity.

| Actual raw file | HSP rows | Fields | Archived exact recorded anchors found | New exact recorded anchors found | Actual raw SHA256 |
| --- | ---: | ---: | ---: | ---: | --- |
| `stage2D_blast_hits.tsv` | 1,230,305 | 14 | 6,359/6,359 | 6,317/6,782 | `2d2153b25a229c0a35e366ed53bd749c9190ab307805abf12eadef1151ee4712` |
| `stage2D_v2_blast_hits.tsv` | 659,348 | 16 | 6,317/6,359 | 6,782/6,782 | `d2ce7eabbc38b9ccc48b25ce67dd06b720bcfcd57f2358b002664ad2b5550798` |

The selected guide-aware raw file's34315 distinct subject accessions match the completed run's required-subject count. Its query FASTA SHA is `af17f1cb2dc65535cbf0a5bc92b8fcf73d24f348466a5b36108e4ee7db3fb44e`;376/376 queries equal actual-guide reverse complements. All659348 aligned raw qseq values match their paired query slices. All6,782 recorded new anchors match the selected file's longest-span/first-tie choice; their raw aligned subject substrings agree with the overlapping stored windows.

The older raw has1,230,305 HSP rows,67228 distinct subjects,and query FASTA SHA `3d46827703f0a5b0c31e0791ad0a593cb98de572831907bd9798b0effe71ac0e`:376/376 query sequences match mutant target segments. Its14-column schema ends in **qlen,slen**, confirmed from recovered stage2D_transcriptome_offtarget.py lines53–55 (source SHA f9360a024d8e17828e184ee7876dbf1de0980c736a13031bf7c6e9b73b84a50a). It does NOT contain aligned qseq/sseq, so no old aligned-query/subject-sequence validation is claimed. An initial comparison wrongly applied aligned-sequence indices to those numeric columns; source inspection resolved that checker error before the final reconciliation. The new16-column file includes qlen,slen,qseq,sseq and was parsed explicitly as such.

All6,359 archived recorded anchors match the older file's longest-span/first-tie choice. Only6,317 of those exact archived anchor tuples occur in the guide-aware file, but6,344 archived **loci** are projectable there:27 shared loci use a different selected HSP. Together with15 missing old loci,this explains the previously reported42 archived anchor tuples absent fromv2. Exact HSP identity and locus coverage must not be conflated.

This establishes current byte/query/anchor consistency. It strongly supports the older raw content for the archived outputs under the documented algorithm; it still does not authenticate the historical command,filename at execution or environment. No whole-reference reconstruction or unsearched transcriptome coverage was independently rerun here.

## Register effects versus anchor-set effects

On all6,344 shared loci,query and subject sequences,mismatch counts and mismatch positions are unchanged. Register differs on340 of these rows. Shared editable flags:6,264 NO→NO,53 YES→NO,21 NO→YES,6 YES→YES.

Of the341 historically wrong-register rows,340 remain shared and one is archive-only. All53 historically positive flags within that subset remain shared and areNO under the corrected register. This is a demonstrated change in computational annotation on stored windows—not53 experimentally confirmed false positives. The14 historically flagged archive-only rows are a separate anchor-set effect.

Across current versus archived stored guide tables:51 guides change total non-target counts;18 change editable counts;60 change raw risk;357 change normalized specificity. Among these357,297 have unchanged raw risk and change through cohort normalization. Old raw-risk extrema were0.75–13.75; new extrema are0.75–7.72. No composite score,rank,winner or robustness statistic was calculated. The overall comparison includes both anchor-source selection and register correction and must not be labelled a register-only sensitivity experiment.

## Ranking-proposal suitability and remaining boundaries

**PASS for a separately reviewed ranking proposal.** The required full-window specificity/bystander columns are present, numeric and complete; frozen structure values and exact guide keys are compatible. Use the explicit new full-window columns and same frozen376-guide cohort; preserve legacy E/F without renaming/converting them or assuming selected winners are unchanged.

The generic inherited Transcriptome_offtarget_status remains NOT_RUN and Ranking_status remains UNRANKED in all376 rows. These are unchanged source placeholders; the new Stage2D_v3_status and full-window fields record this run. Do not silently rewrite them or let downstream code select generic legacy status/features. A ranking adapter should bind the exact structure/full-window inputs and columns and fail on missing/ambiguous/nonfinite values. Existing legacy top tables contain multiple guides per variant; future comparisons must explicitly define the baseline winner rather than use dictionary last-row overwrite. These are requirements for the separate ranking proposal,not execution performed now.

Historical provenance,unsearched-locus sensitivity and experimental efficacy remain unestablished. Reference authentication from prior evidence is carried forward; only stored windows/local raw overlaps were checked in this task,not every base against a newly opened full reference FASTA. Production-host enforcement is not inferred from either this verification or the earlier isolated tests.

## Optional direct Windows confirmation — read only

The uploaded copies already match the posted manifest; no repeat production run is needed. To confirm the actual Windows paths still contain those bytes,run this read-only PowerShell block and paste its JSON only if you want that additional current-path evidence. It writes no files and performs no scientific processing.

```powershell
$ErrorActionPreference = 'Stop'
$run = '[LOCAL_PATH]
$manifestPath = Join-Path $run 'output_manifest.json'
$manifestHash = (Get-FileHash -LiteralPath $manifestPath -Algorithm SHA256).Hash.ToLowerInvariant()
if ($manifestHash -ne '48bd98e1188b0d3f3cb340282a6fec641804d643227544d16cb96cb519580249') { throw 'Manifest hash mismatch' }
$manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
$specs = @(
    @('outputs/MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv', ',', 376),
    @('outputs/MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_summary.csv', ',', 8),
    @('outputs/stage2D_v3_reconstructed_hits.tsv', "`t", 6782)
)
$report = @()
foreach ($spec in $specs) {
    $path = Join-Path $run $spec[0]
    $entries = @($manifest.files | Where-Object { $_.path -ceq $spec[0] })
    if ($entries.Count -ne 1) { throw "Ambiguous/missing manifest entry: $path" }
    $entry = $entries[0]
    $hash = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
    $bytes = (Get-Item -LiteralPath $path).Length
    $rows = @(Import-Csv -LiteralPath $path -Delimiter $spec[1])
    if ($hash -ne $entry.sha256 -or $bytes -ne $entry.bytes -or $rows.Count -ne $spec[2]) { throw "Local inspection mismatch: $path" }
    $report += [pscustomobject]@{path=$path;SHA256=$hash;bytes=$bytes;rows=$rows.Count;columns=$rows[0].PSObject.Properties.Name.Count;status='MATCH'}
}
[pscustomobject]@{scope='READ_ONLY_LOCAL_CONFIRMATION';manifest_SHA256=$manifestHash;files=$report;production_processing_performed=$false} | ConvertTo-Json -Depth 6
```

## Preserved evidence and reproducibility

Sources read: the four current uploads; exact candidate/B/C members from MECP2_Stage1_to_Stage3_Recovery_Bundle.zip; the three relevant Tier1 archived output members; two actual Tier1 raw TSVs and paired FASTAs from MECP2_Historical_BLAST_Evidence.zip; the recovered producer source. No unrelated nonsense-window processing,BLAST,ranking,robustness or promotion was performed. Existing files were not modified. Independent checks and intermediate locus-difference records are stored in this continuing Phase3F verification workspace. This standalone report is supplied instead of another audit-only ZIP.
