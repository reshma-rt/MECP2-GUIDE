"""Scientific body extracted from pinned Stage2F v2 AST; no file I/O.
Synthetic tests may override n_sim; the real CLI fixes10000 and seed20260914.
"""
import math
import random
from collections import Counter, defaultdict
COMPONENTS = ['Stage2E_v2_accessibility_component_0_1', 'Stage2E_v2_bystander_component_0_1', 'Stage2E_v2_specificity_component_0_1', 'Stage2E_v2_duplex_component_0_1', 'Stage2E_v2_self_structure_component_0_1']
BASELINE_SCORE = 'Stage2E_v2_composite_score_0_100'
BASELINE_RANK = 'Stage2E_v2_rank_within_mutation'
SCENARIOS = {'baseline': [0.25, 0.3, 0.3, 0.1, 0.05], 'equal': [0.2, 0.2, 0.2, 0.2, 0.2], 'specificity_heavy': [0.15, 0.2, 0.45, 0.1, 0.1], 'bystander_heavy': [0.15, 0.45, 0.2, 0.1, 0.1], 'accessibility_heavy': [0.45, 0.2, 0.2, 0.1, 0.05], 'structure_heavy': [0.15, 0.2, 0.2, 0.25, 0.2]}

def fnum(v):
    x = float(v)
    if math.isnan(x) or math.isinf(x):
        raise ValueError(f'Non-finite numeric value: {v}')
    return x

def dirichlet_ones(k, rng):
    vals = [rng.expovariate(1.0) for _ in range(k)]
    s = sum(vals)
    return [v / s for v in vals]

def calculate(rows, n_sim=10000, seed=20260914):
    N_SIM = n_sim
    SEED = seed
    guide_data = {}
    by_mut = defaultdict(list)
    for r in rows:
        gid = r['Guide_ID']
        comps = [fnum(r[c]) for c in COMPONENTS]
        if any((v < 0 or v > 1 for v in comps)):
            raise ValueError(f'{gid}: component score outside [0,1]: {comps}')
        d = {'row': r, 'components': comps, 'baseline_score': fnum(r[BASELINE_SCORE]), 'baseline_rank': int(float(r[BASELINE_RANK]))}
        guide_data[gid] = d
        by_mut[r['VariationID']].append(gid)
    for vid, gids in by_mut.items():
        if len(gids) != 47:
            raise RuntimeError(f'VariationID {vid}: expected 47 guides, found {len(gids)}')
    top_counts = Counter()
    rank_sums = Counter()
    scenario_top_counts = Counter()
    scenario_results = []
    rng = random.Random(SEED)
    for sim_i in range(N_SIM):
        w = dirichlet_ones(5, rng)
        for vid, gids in by_mut.items():
            scored = []
            for gid in gids:
                c = guide_data[gid]['components']
                score = sum((wi * ci for wi, ci in zip(w, c)))
                scored.append((score, gid))
            scored.sort(key=lambda x: (-x[0], x[1]))
            for rank, (_, gid) in enumerate(scored, 1):
                rank_sums[gid] += rank
            top_counts[scored[0][1]] += 1
    for sname, w in SCENARIOS.items():
        for vid, gids in by_mut.items():
            scored = []
            for gid in gids:
                c = guide_data[gid]['components']
                score = sum((wi * ci for wi, ci in zip(w, c)))
                scored.append((score, gid))
            scored.sort(key=lambda x: (-x[0], x[1]))
            top_gid = scored[0][1]
            scenario_top_counts[top_gid] += 1
            scenario_results.append({'Scenario': sname, 'VariationID': vid, 'Top_Guide_ID': top_gid, 'Scenario_score_0_1': scored[0][0], 'Weight_accessibility': w[0], 'Weight_bystander': w[1], 'Weight_specificity': w[2], 'Weight_duplex': w[3], 'Weight_self_structure': w[4]})
    annotated = []
    for r in rows:
        gid = r['Guide_ID']
        top_fraction = top_counts[gid] / N_SIM
        mean_rank = rank_sums[gid] / N_SIM
        if top_fraction >= 0.5:
            label = 'HIGH_STABILITY'
        elif top_fraction >= 0.2:
            label = 'MODERATE_STABILITY'
        else:
            label = 'WEIGHT_SENSITIVE'
        x = dict(r)
        x['Sensitivity_random_top_count'] = top_counts[gid]
        x['Sensitivity_random_top_fraction'] = top_fraction
        x['Sensitivity_mean_rank'] = mean_rank
        x['Sensitivity_scenario_top_count_of_6'] = scenario_top_counts[gid]
        x['Sensitivity_stability_label'] = label
        x['Stage2F_v2_status'] = 'HARMONIZED_WEIGHT_SENSITIVITY_COMPLETE'
        annotated.append(x)
    ann_by_mut = defaultdict(list)
    for r in annotated:
        ann_by_mut[r['VariationID']].append(r)
    robust_rows = []
    mutation_summary = []
    for vid, grp in ann_by_mut.items():
        grp_sorted = sorted(grp, key=lambda r: (-float(r['Sensitivity_random_top_fraction']), float(r['Sensitivity_mean_rank']), -float(r[BASELINE_SCORE]), r['Guide_ID']))
        robust = grp_sorted[0]
        robust_rows.append(robust)
        baseline = min(grp, key=lambda r: int(float(r[BASELINE_RANK])))
        same = 'YES' if robust['Guide_ID'] == baseline['Guide_ID'] else 'NO'
        mutation_summary.append({'VariationID': vid, 'Stage2E_v2_top_Guide_ID': baseline['Guide_ID'], 'Stage2F_v2_robust_top_Guide_ID': robust['Guide_ID'], 'Stage2E_top_remains_robust_top': same, 'Robust_top_random_top_fraction': robust['Sensitivity_random_top_fraction'], 'Robust_top_mean_rank': robust['Sensitivity_mean_rank'], 'Robust_top_scenario_top_count_of_6': robust['Sensitivity_scenario_top_count_of_6'], 'Robust_top_stability_label': robust['Sensitivity_stability_label'], 'Robust_top_editable_A_offtarget_hits_le3mm': robust.get('Fullwindow_editable_A_offtarget_hits_le3mm', ''), 'Robust_top_Stage2E_v2_score_0_100': robust[BASELINE_SCORE]})
    annotated.sort(key=lambda r: (str(r['VariationID']), -float(r['Sensitivity_random_top_fraction']), float(r['Sensitivity_mean_rank']), r['Guide_ID']))
    robust_rows.sort(key=lambda r: str(r['VariationID']))
    mutation_summary.sort(key=lambda r: str(r['VariationID']))
    return (annotated, robust_rows, scenario_results, mutation_summary, dict(top_counts), dict(rank_sums))
