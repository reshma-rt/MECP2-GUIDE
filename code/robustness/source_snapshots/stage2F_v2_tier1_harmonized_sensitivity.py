#!/usr/bin/env python3
"""
MECP2 ADAR Dataset 2 — Tier 1
Stage 2F v2: Sensitivity analysis after harmonized Stage 2E v2 ranking

Purpose
-------
Run the same weight-sensitivity framework used previously, but on the
harmonized Stage 2E v2 component scores that incorporate corrected
full-window transcriptome specificity.

Inputs
------
02_tier1_pipeline/
    MECP2_ADAR_Dataset2_Tier1_stage2E_v2_ranked.csv

Method
------
- 10,000 random Dirichlet(1,1,1,1,1) weight vectors
- Seed: 20260914
- Components:
    accessibility
    bystander precision
    full-window transcriptome specificity
    duplex
    guide self-structure
- Six fixed scenarios:
    baseline
    equal
    specificity-heavy
    bystander-heavy
    accessibility-heavy
    structure-heavy
- Stability label by random top fraction:
    >= 0.50  HIGH_STABILITY
    >= 0.20  MODERATE_STABILITY
    < 0.20   WEIGHT_SENSITIVE
- Robust top per mutation:
    highest random top fraction,
    then best mean rank,
    then highest Stage 2E v2 baseline score.

Outputs
-------
02_tier1_pipeline/
    MECP2_ADAR_Dataset2_Tier1_stage2F_v2_sensitivity_all_guides.csv
    MECP2_ADAR_Dataset2_Tier1_stage2F_v2_robust_top_guides.csv
    MECP2_ADAR_Dataset2_Tier1_stage2F_v2_scenarios.csv
    MECP2_ADAR_Dataset2_Tier1_stage2F_v2_mutation_summary.csv
"""

import csv
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "02_tier1_pipeline"

INPUT = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2E_v2_ranked.csv"

OUT_ALL = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2F_v2_sensitivity_all_guides.csv"
OUT_ROBUST = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2F_v2_robust_top_guides.csv"
OUT_SCEN = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2F_v2_scenarios.csv"
OUT_MUT = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2F_v2_mutation_summary.csv"

EXPECTED_GUIDES = 376
EXPECTED_MUTATIONS = 8
N_SIM = 10000
SEED = 20260914

COMPONENTS = [
    "Stage2E_v2_accessibility_component_0_1",
    "Stage2E_v2_bystander_component_0_1",
    "Stage2E_v2_specificity_component_0_1",
    "Stage2E_v2_duplex_component_0_1",
    "Stage2E_v2_self_structure_component_0_1",
]

BASELINE_SCORE = "Stage2E_v2_composite_score_0_100"
BASELINE_RANK = "Stage2E_v2_rank_within_mutation"

SCENARIOS = {
    "baseline":            [0.25, 0.30, 0.30, 0.10, 0.05],
    "equal":               [0.20, 0.20, 0.20, 0.20, 0.20],
    "specificity_heavy":   [0.15, 0.20, 0.45, 0.10, 0.10],
    "bystander_heavy":     [0.15, 0.45, 0.20, 0.10, 0.10],
    "accessibility_heavy": [0.45, 0.20, 0.20, 0.10, 0.05],
    "structure_heavy":     [0.15, 0.20, 0.20, 0.25, 0.20],
}


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
    x = float(v)
    if math.isnan(x) or math.isinf(x):
        raise ValueError(f"Non-finite numeric value: {v}")
    return x


def dirichlet_ones(k, rng):
    # Gamma(1,1) == exponential(1)
    vals = [rng.expovariate(1.0) for _ in range(k)]
    s = sum(vals)
    return [v / s for v in vals]


if not INPUT.exists():
    raise FileNotFoundError(f"Missing Stage 2E v2 ranked file: {INPUT}")

rows, fields = read_csv(INPUT)

if len(rows) != EXPECTED_GUIDES:
    raise ValueError(f"Expected {EXPECTED_GUIDES} guides, found {len(rows)}.")

if len({r["Guide_ID"] for r in rows}) != EXPECTED_GUIDES:
    raise ValueError("Duplicate Guide_ID values detected.")

if len({r["VariationID"] for r in rows}) != EXPECTED_MUTATIONS:
    raise ValueError(
        f"Expected {EXPECTED_MUTATIONS} mutations, found "
        f"{len({r['VariationID'] for r in rows})}."
    )

for c in COMPONENTS + [BASELINE_SCORE, BASELINE_RANK]:
    if c not in fields:
        raise KeyError(f"Required Stage 2E v2 column missing: {c}")

# Parse component matrix.
guide_data = {}
by_mut = defaultdict(list)

for r in rows:
    gid = r["Guide_ID"]
    comps = [fnum(r[c]) for c in COMPONENTS]
    if any(v < 0 or v > 1 for v in comps):
        raise ValueError(f"{gid}: component score outside [0,1]: {comps}")

    d = {
        "row": r,
        "components": comps,
        "baseline_score": fnum(r[BASELINE_SCORE]),
        "baseline_rank": int(float(r[BASELINE_RANK])),
    }
    guide_data[gid] = d
    by_mut[r["VariationID"]].append(gid)

# Hard QC: 47 guides per mutation.
for vid, gids in by_mut.items():
    if len(gids) != 47:
        raise RuntimeError(f"VariationID {vid}: expected 47 guides, found {len(gids)}")

# Accumulators.
top_counts = Counter()
rank_sums = Counter()
scenario_top_counts = Counter()
scenario_results = []

rng = random.Random(SEED)

# ---------------------------------------------------------------------
# 10,000 random weight simulations.
# ---------------------------------------------------------------------
for sim_i in range(N_SIM):
    w = dirichlet_ones(5, rng)

    for vid, gids in by_mut.items():
        scored = []
        for gid in gids:
            c = guide_data[gid]["components"]
            score = sum(wi * ci for wi, ci in zip(w, c))
            scored.append((score, gid))

        # Deterministic tie-break by Guide_ID.
        scored.sort(key=lambda x: (-x[0], x[1]))

        for rank, (_, gid) in enumerate(scored, 1):
            rank_sums[gid] += rank

        top_counts[scored[0][1]] += 1

# ---------------------------------------------------------------------
# Fixed scenarios.
# ---------------------------------------------------------------------
for sname, w in SCENARIOS.items():
    for vid, gids in by_mut.items():
        scored = []
        for gid in gids:
            c = guide_data[gid]["components"]
            score = sum(wi * ci for wi, ci in zip(w, c))
            scored.append((score, gid))

        scored.sort(key=lambda x: (-x[0], x[1]))
        top_gid = scored[0][1]
        scenario_top_counts[top_gid] += 1

        scenario_results.append({
            "Scenario": sname,
            "VariationID": vid,
            "Top_Guide_ID": top_gid,
            "Scenario_score_0_1": scored[0][0],
            "Weight_accessibility": w[0],
            "Weight_bystander": w[1],
            "Weight_specificity": w[2],
            "Weight_duplex": w[3],
            "Weight_self_structure": w[4],
        })

# ---------------------------------------------------------------------
# Annotate all guides.
# ---------------------------------------------------------------------
annotated = []

for r in rows:
    gid = r["Guide_ID"]
    top_fraction = top_counts[gid] / N_SIM
    mean_rank = rank_sums[gid] / N_SIM

    if top_fraction >= 0.50:
        label = "HIGH_STABILITY"
    elif top_fraction >= 0.20:
        label = "MODERATE_STABILITY"
    else:
        label = "WEIGHT_SENSITIVE"

    x = dict(r)
    x["Sensitivity_random_top_count"] = top_counts[gid]
    x["Sensitivity_random_top_fraction"] = top_fraction
    x["Sensitivity_mean_rank"] = mean_rank
    x["Sensitivity_scenario_top_count_of_6"] = scenario_top_counts[gid]
    x["Sensitivity_stability_label"] = label
    x["Stage2F_v2_status"] = "HARMONIZED_WEIGHT_SENSITIVITY_COMPLETE"
    annotated.append(x)

# ---------------------------------------------------------------------
# Pick robust top guide per mutation.
# ---------------------------------------------------------------------
ann_by_mut = defaultdict(list)
for r in annotated:
    ann_by_mut[r["VariationID"]].append(r)

robust_rows = []
mutation_summary = []

for vid, grp in ann_by_mut.items():
    grp_sorted = sorted(
        grp,
        key=lambda r: (
            -float(r["Sensitivity_random_top_fraction"]),
            float(r["Sensitivity_mean_rank"]),
            -float(r[BASELINE_SCORE]),
            r["Guide_ID"],
        )
    )

    robust = grp_sorted[0]
    robust_rows.append(robust)

    baseline = min(grp, key=lambda r: int(float(r[BASELINE_RANK])))

    same = "YES" if robust["Guide_ID"] == baseline["Guide_ID"] else "NO"

    mutation_summary.append({
        "VariationID": vid,
        "Stage2E_v2_top_Guide_ID": baseline["Guide_ID"],
        "Stage2F_v2_robust_top_Guide_ID": robust["Guide_ID"],
        "Stage2E_top_remains_robust_top": same,
        "Robust_top_random_top_fraction": robust["Sensitivity_random_top_fraction"],
        "Robust_top_mean_rank": robust["Sensitivity_mean_rank"],
        "Robust_top_scenario_top_count_of_6": robust["Sensitivity_scenario_top_count_of_6"],
        "Robust_top_stability_label": robust["Sensitivity_stability_label"],
        "Robust_top_editable_A_offtarget_hits_le3mm": robust.get(
            "Fullwindow_editable_A_offtarget_hits_le3mm", ""
        ),
        "Robust_top_Stage2E_v2_score_0_100": robust[BASELINE_SCORE],
    })

annotated.sort(
    key=lambda r: (
        str(r["VariationID"]),
        -float(r["Sensitivity_random_top_fraction"]),
        float(r["Sensitivity_mean_rank"]),
        r["Guide_ID"],
    )
)
robust_rows.sort(key=lambda r: str(r["VariationID"]))
mutation_summary.sort(key=lambda r: str(r["VariationID"]))

# ---------------------------------------------------------------------
# Write outputs.
# ---------------------------------------------------------------------
extra_fields = [
    "Sensitivity_random_top_count",
    "Sensitivity_random_top_fraction",
    "Sensitivity_mean_rank",
    "Sensitivity_scenario_top_count_of_6",
    "Sensitivity_stability_label",
    "Stage2F_v2_status",
]

out_fields = list(fields)
for c in extra_fields:
    if c not in out_fields:
        out_fields.append(c)

write_csv(OUT_ALL, annotated, out_fields)
write_csv(OUT_ROBUST, robust_rows, out_fields)

scenario_fields = [
    "Scenario",
    "VariationID",
    "Top_Guide_ID",
    "Scenario_score_0_1",
    "Weight_accessibility",
    "Weight_bystander",
    "Weight_specificity",
    "Weight_duplex",
    "Weight_self_structure",
]
write_csv(OUT_SCEN, scenario_results, scenario_fields)

mutation_fields = [
    "VariationID",
    "Stage2E_v2_top_Guide_ID",
    "Stage2F_v2_robust_top_Guide_ID",
    "Stage2E_top_remains_robust_top",
    "Robust_top_random_top_fraction",
    "Robust_top_mean_rank",
    "Robust_top_scenario_top_count_of_6",
    "Robust_top_stability_label",
    "Robust_top_editable_A_offtarget_hits_le3mm",
    "Robust_top_Stage2E_v2_score_0_100",
]
write_csv(OUT_MUT, mutation_summary, mutation_fields)

# ---------------------------------------------------------------------
# Console summary.
# ---------------------------------------------------------------------
labels = Counter(r["Sensitivity_stability_label"] for r in annotated)
robust_labels = Counter(r["Sensitivity_stability_label"] for r in robust_rows)
retained = sum(
    r["Stage2E_top_remains_robust_top"] == "YES"
    for r in mutation_summary
)

print("\n" + "=" * 88)
print("STAGE 2F v2 COMPLETE — TIER 1 HARMONIZED SENSITIVITY ANALYSIS")
print("=" * 88)
print("Guides analyzed                          :", len(annotated))
print("Mutations                                :", len(ann_by_mut))
print("Random weight simulations                :", N_SIM)
print("Random seed                              :", SEED)
print("Fixed scenarios                          :", len(SCENARIOS))
print()
print("All-guide stability labels:")
for k in ["HIGH_STABILITY", "MODERATE_STABILITY", "WEIGHT_SENSITIVE"]:
    print(f"  {k:<20}:", labels.get(k, 0))
print()
print("Mutation-level robust candidates:")
print("  Stage 2E v2 top remained robust top    :", retained, "/", EXPECTED_MUTATIONS)
for k in ["HIGH_STABILITY", "MODERATE_STABILITY", "WEIGHT_SENSITIVE"]:
    print(f"  Robust top {k:<20}:", robust_labels.get(k, 0))

print("\nRobust top guide per mutation:")
for r in robust_rows:
    print(
        f"  {r['VariationID']:>10}  {r['Guide_ID']:<10} "
        f"top_fraction={float(r['Sensitivity_random_top_fraction']):.4f}  "
        f"mean_rank={float(r['Sensitivity_mean_rank']):.3f}  "
        f"stability={r['Sensitivity_stability_label']}  "
        f"editable-A-offtargets={r.get('Fullwindow_editable_A_offtarget_hits_le3mm','')}"
    )

print("\nOutputs:")
print("  All guides       :", OUT_ALL)
print("  Robust top guides:", OUT_ROBUST)
print("  Scenarios        :", OUT_SCEN)
print("  Mutation summary :", OUT_MUT)

print(
    "\nInterpretation: primary Tier 1 candidates should now be taken from "
    "Stage 2F v2 robust-top guides. Stage 2E v2 rank-1 remains the baseline "
    "weighted candidate. Stability reflects ranking robustness to weight "
    "uncertainty, not experimental efficacy."
)
