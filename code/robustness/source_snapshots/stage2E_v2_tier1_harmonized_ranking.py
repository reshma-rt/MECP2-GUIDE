#!/usr/bin/env python3
"""
MECP2 ADAR Dataset 2 — Tier 1
Stage 2E v2: Re-ranking after harmonized full-window off-target analysis

Purpose
-------
Recompute Tier 1 multi-objective ranking using the corrected Stage 2D v3
full-window transcriptome specificity score.

The guide set is unchanged (376 guides across 8 exact-restoration mutations).
Only the transcriptome-specificity component is replaced by the harmonized
full-window metric.

Inputs
------
02_tier1_pipeline/
    MECP2_ADAR_Dataset2_Tier1_stage2B_structure.csv
    MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv

The Stage 2D v3 file already carries Stage 2C bystander columns.

Scoring
-------
Accessibility component:
    mean(Target_A_unpaired_probability,
         Target_A_5nt_mean_unpaired_probability)

Bystander component:
    Bystander_precision_score_0_1

Transcriptome specificity:
    Fullwindow_GuideAware_transcriptome_specificity_score_0_1

Duplex component:
    min-max normalized duplex MFE, MORE NEGATIVE favored

Guide self-structure component:
    min-max normalized guide MFE, LESS NEGATIVE favored

Weights:
    0.25 accessibility
    0.30 bystander precision
    0.30 transcriptome specificity
    0.10 duplex
    0.05 guide self-structure

Outputs
-------
02_tier1_pipeline/
    MECP2_ADAR_Dataset2_Tier1_stage2E_v2_ranked.csv
    MECP2_ADAR_Dataset2_Tier1_stage2E_v2_top_guides.csv
    MECP2_ADAR_Dataset2_Tier1_stage2E_v2_summary.csv

Interpretation
--------------
Scores are computational, cohort-relative prioritization metrics and are not
experimental editing efficiencies.
"""

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "02_tier1_pipeline"

STRUCT = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2B_structure.csv"
OFFTARGET = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv"

OUT_RANKED = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2E_v2_ranked.csv"
OUT_TOP = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2E_v2_top_guides.csv"
OUT_SUM = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2E_v2_summary.csv"

EXPECTED_GUIDES = 376
EXPECTED_MUTATIONS = 8

W_ACCESS = 0.25
W_BYSTANDER = 0.30
W_SPECIFICITY = 0.30
W_DUPLEX = 0.10
W_SELF = 0.05


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        rows = list(rd)
        return rows, list(rd.fieldnames or [])


def write_csv(path, rows, fields):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)


def fnum(v):
    try:
        x = float(v)
        if math.isnan(x) or math.isinf(x):
            return None
        return x
    except Exception:
        return None


def find_col(fields, exact=None, tokens=None, excludes=None, required=False):
    exact = exact or []
    tokens = [x.lower() for x in (tokens or [])]
    excludes = [x.lower() for x in (excludes or [])]

    for c in exact:
        if c in fields:
            return c

    matches = []
    for f in fields:
        fl = f.lower()
        if tokens and not all(t in fl for t in tokens):
            continue
        if any(t in fl for t in excludes):
            continue
        matches.append(f)

    if matches:
        return sorted(matches, key=lambda x: (len(x), x.lower()))[0]

    if required:
        raise KeyError(
            f"Required column not found. exact={exact}, tokens={tokens}, excludes={excludes}"
        )
    return None


def minmax(values, higher_is_better=True):
    vals = [v for v in values if v is not None]
    if not vals:
        raise ValueError("No numeric values supplied for normalization.")

    lo, hi = min(vals), max(vals)

    out = []
    for v in values:
        if v is None:
            out.append(None)
            continue
        if hi == lo:
            score = 1.0
        else:
            score = (v - lo) / (hi - lo)
        if not higher_is_better:
            score = 1.0 - score
        out.append(score)
    return out, lo, hi


if not STRUCT.exists():
    raise FileNotFoundError(f"Missing Stage 2B structure file: {STRUCT}")
if not OFFTARGET.exists():
    raise FileNotFoundError(f"Missing harmonized Stage 2D v3 file: {OFFTARGET}")

struct_rows, struct_fields = read_csv(STRUCT)
off_rows, off_fields = read_csv(OFFTARGET)

if len(off_rows) != EXPECTED_GUIDES:
    raise ValueError(f"Expected {EXPECTED_GUIDES} Stage 2D v3 guides, found {len(off_rows)}.")

if len({r["Guide_ID"] for r in off_rows}) != EXPECTED_GUIDES:
    raise ValueError("Stage 2D v3 contains duplicate Guide_ID values.")

if len({r["VariationID"] for r in off_rows}) != EXPECTED_MUTATIONS:
    raise ValueError("Expected 8 Tier 1 mutations in Stage 2D v3.")

struct_by_id = {r["Guide_ID"]: r for r in struct_rows}

missing_structure = [r["Guide_ID"] for r in off_rows if r["Guide_ID"] not in struct_by_id]
if missing_structure:
    raise RuntimeError(
        f"{len(missing_structure)} Stage 2D guides are missing from Stage 2B structure file. "
        f"First: {missing_structure[:20]}"
    )

# Resolve structure columns.
target_a_col = find_col(
    struct_fields,
    exact=["Target_A_unpaired_probability"],
    tokens=["target", "a", "unpaired", "probability"],
    excludes=["5nt", "3nt", "11nt", "mean"],
    required=True,
)

five_col = find_col(
    struct_fields,
    exact=["Target_A_5nt_mean_unpaired_probability"],
    tokens=["5nt", "unpaired"],
    required=True,
)

duplex_col = find_col(
    struct_fields,
    exact=[
        "Guide_target_duplex_MFE_kcal_per_mol",
        "Guide_target_duplex_MFE",
        "Duplex_MFE_kcal_per_mol",
    ],
    tokens=["duplex", "mfe"],
    required=True,
)

guide_mfe_col = find_col(
    struct_fields,
    exact=["Guide_MFE_kcal_per_mol", "Guide_MFE"],
    tokens=["guide", "mfe"],
    excludes=["duplex", "target"],
    required=True,
)

# Resolve Stage 2C / Stage 2D v3 scoring columns.
bystander_col = find_col(
    off_fields,
    exact=["Bystander_precision_score_0_1"],
    tokens=["bystander", "precision"],
    required=True,
)

specificity_col = find_col(
    off_fields,
    exact=["Fullwindow_GuideAware_transcriptome_specificity_score_0_1"],
    tokens=["fullwindow", "specificity"],
    required=True,
)

print("STAGE 2E v2 — resolved scoring columns")
print("  Target-A unpaired       :", target_a_col)
print("  5-nt accessibility      :", five_col)
print("  Bystander precision     :", bystander_col)
print("  Full-window specificity :", specificity_col)
print("  Duplex MFE              :", duplex_col)
print("  Guide MFE               :", guide_mfe_col)

# Merge structure into Stage 2D v3 rows.
rows = []
for r in off_rows:
    x = dict(r)
    s = struct_by_id[r["Guide_ID"]]
    for c in struct_fields:
        if c not in x:
            x[c] = s.get(c, "")
    rows.append(x)

target_a_vals = []
five_vals = []
bystander_vals = []
specificity_vals = []
duplex_vals = []
guide_mfe_vals = []

for r in rows:
    s = struct_by_id[r["Guide_ID"]]

    target_a_vals.append(fnum(s.get(target_a_col)))
    five_vals.append(fnum(s.get(five_col)))
    bystander_vals.append(fnum(r.get(bystander_col)))
    specificity_vals.append(fnum(r.get(specificity_col)))
    duplex_vals.append(fnum(s.get(duplex_col)))
    guide_mfe_vals.append(fnum(s.get(guide_mfe_col)))

# Hard completeness QC.
component_arrays = {
    "target-A unpaired": target_a_vals,
    "5nt accessibility": five_vals,
    "bystander precision": bystander_vals,
    "full-window specificity": specificity_vals,
    "duplex MFE": duplex_vals,
    "guide MFE": guide_mfe_vals,
}
for name, arr in component_arrays.items():
    nmiss = sum(v is None for v in arr)
    if nmiss:
        raise RuntimeError(f"{name}: {nmiss} missing/non-numeric values.")

# Accessibility is the mean of two already 0-1 probabilities.
accessibility = [(a + b) / 2.0 for a, b in zip(target_a_vals, five_vals)]

# Structural normalization:
# duplex: more negative = better
# guide self-MFE: less negative = better
duplex_score, duplex_lo, duplex_hi = minmax(duplex_vals, higher_is_better=False)
self_score, self_lo, self_hi = minmax(guide_mfe_vals, higher_is_better=True)

for i, r in enumerate(rows):
    a = accessibility[i]
    b = bystander_vals[i]
    sp = specificity_vals[i]
    d = duplex_score[i]
    gs = self_score[i]

    composite_01 = (
        W_ACCESS * a
        + W_BYSTANDER * b
        + W_SPECIFICITY * sp
        + W_DUPLEX * d
        + W_SELF * gs
    )

    r["Stage2E_v2_accessibility_component_0_1"] = a
    r["Stage2E_v2_bystander_component_0_1"] = b
    r["Stage2E_v2_specificity_component_0_1"] = sp
    r["Stage2E_v2_duplex_component_0_1"] = d
    r["Stage2E_v2_self_structure_component_0_1"] = gs
    r["Stage2E_v2_composite_score_0_100"] = composite_01 * 100.0

    r["Stage2E_v2_weight_accessibility"] = W_ACCESS
    r["Stage2E_v2_weight_bystander"] = W_BYSTANDER
    r["Stage2E_v2_weight_specificity"] = W_SPECIFICITY
    r["Stage2E_v2_weight_duplex"] = W_DUPLEX
    r["Stage2E_v2_weight_self_structure"] = W_SELF
    r["Stage2E_v2_status"] = "HARMONIZED_FULLWINDOW_SPECIFICITY_RANKING_COMPLETE"

# Rank separately within each mutation.
by_mut = defaultdict(list)
for r in rows:
    by_mut[r["VariationID"]].append(r)

top_rows = []
for vid, grp in by_mut.items():
    grp.sort(
        key=lambda r: (
            -float(r["Stage2E_v2_composite_score_0_100"]),
            int(r.get("Fullwindow_editable_A_offtarget_hits_le3mm", 10**9) or 10**9),
            r["Guide_ID"],
        )
    )
    for rank, r in enumerate(grp, 1):
        r["Stage2E_v2_rank_within_mutation"] = rank
    top_rows.append(grp[0])

# Restore deterministic all-guide ordering.
rows.sort(key=lambda r: (str(r["VariationID"]), int(r["Stage2E_v2_rank_within_mutation"]), r["Guide_ID"]))
top_rows.sort(key=lambda r: str(r["VariationID"]))

# Compare to old Stage 2E top guides if present.
old_top_path = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2E_top_guides.csv"
old_top = {}
if old_top_path.exists():
    old_rows, _ = read_csv(old_top_path)
    old_top = {r["VariationID"]: r["Guide_ID"] for r in old_rows}

summary_rows = []
changed = 0

for r in top_rows:
    vid = r["VariationID"]
    old_gid = old_top.get(vid, "")
    new_gid = r["Guide_ID"]
    same = "" if not old_gid else ("YES" if old_gid == new_gid else "NO")
    if same == "NO":
        changed += 1

    summary_rows.append({
        "VariationID": vid,
        "Stage2E_v2_top_Guide_ID": new_gid,
        "Stage2E_v2_top_score_0_100": r["Stage2E_v2_composite_score_0_100"],
        "Editable_A_offtarget_hits_le3mm": r.get("Fullwindow_editable_A_offtarget_hits_le3mm", ""),
        "Old_Stage2E_top_Guide_ID": old_gid,
        "Old_top_retained": same,
    })

# Output fields.
new_fields = [
    "Stage2E_v2_accessibility_component_0_1",
    "Stage2E_v2_bystander_component_0_1",
    "Stage2E_v2_specificity_component_0_1",
    "Stage2E_v2_duplex_component_0_1",
    "Stage2E_v2_self_structure_component_0_1",
    "Stage2E_v2_composite_score_0_100",
    "Stage2E_v2_rank_within_mutation",
    "Stage2E_v2_weight_accessibility",
    "Stage2E_v2_weight_bystander",
    "Stage2E_v2_weight_specificity",
    "Stage2E_v2_weight_duplex",
    "Stage2E_v2_weight_self_structure",
    "Stage2E_v2_status",
]

out_fields = list(off_fields)
for c in struct_fields:
    if c not in out_fields:
        out_fields.append(c)
for c in new_fields:
    if c not in out_fields:
        out_fields.append(c)

write_csv(OUT_RANKED, rows, out_fields)
write_csv(OUT_TOP, top_rows, out_fields)

sum_fields = [
    "VariationID",
    "Stage2E_v2_top_Guide_ID",
    "Stage2E_v2_top_score_0_100",
    "Editable_A_offtarget_hits_le3mm",
    "Old_Stage2E_top_Guide_ID",
    "Old_top_retained",
]
write_csv(OUT_SUM, summary_rows, sum_fields)

print("\n" + "=" * 88)
print("STAGE 2E v2 COMPLETE — TIER 1 HARMONIZED RE-RANKING")
print("=" * 88)
print("Guides ranked                         :", len(rows))
print("Mutations                             :", len(by_mut))
print("Top guides                            :", len(top_rows))
print("Missing scoring components            : 0")
print("Duplex MFE normalization range        :", duplex_lo, "to", duplex_hi)
print("Guide MFE normalization range         :", self_lo, "to", self_hi)

if old_top:
    print("Old Stage 2E top guides retained      :", EXPECTED_MUTATIONS - changed, "/", EXPECTED_MUTATIONS)
    print("Top-guide changes after harmonization :", changed, "/", EXPECTED_MUTATIONS)
else:
    print("Old Stage 2E top comparison           : old top-guide file not found")

print("\nNew top guide per mutation:")
for r in top_rows:
    print(
        f"  {r['VariationID']:>10}  {r['Guide_ID']:<10} "
        f"score={float(r['Stage2E_v2_composite_score_0_100']):.3f}  "
        f"editable-A-offtargets={r.get('Fullwindow_editable_A_offtarget_hits_le3mm','')}"
    )

print("\nOutputs:")
print("  Ranked :", OUT_RANKED)
print("  Top    :", OUT_TOP)
print("  Summary:", OUT_SUM)

print(
    "\nInterpretation: Stage 2E v2 is a computational within-cohort ranking. "
    "The composite score is not an experimentally measured editing efficiency."
)
