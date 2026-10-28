# Methods summary

Guide lengths: 21/25/31 nt; observed 11/15/21 designs per mutation (47 total). Coordinate q=L-g+1 from explicit guide-relative position. Candidate generation is frozen; original Stage2A implementation is not recovered.

Composite score: 100*(0.25 accessibility + 0.30 bystander precision + 0.30 full-window specificity + 0.10 duplex + 0.05 self-structure). Accessibility averages target-A and five-base mean unpaired probabilities and is constant within mutation. Duplex favors more negative MFE; self-structure favors less negative guide MFE, with cohort-relative min–max normalization. Baseline ties: editable-A count then Guide_ID. Full-window similarity is BLAST-anchor defined and uses <=3 mismatches; matching intended-accession loci are excluded from non-target counts.

Robustness: 10,000 prescribed random normalized-exponential weight vectors (Dirichlet ones), seed 20260914; shared weights per iteration across targets. Simulation ties use Guide_ID. Win-frequency thresholds: High >=0.50, Moderate >=0.20, otherwise weight-sensitive. Robust selection: descending win frequency, ascending mean rank, descending baseline score, then ID.

Ablation: baseline control, four component-only rules, and four leave-one-component-out rules. Remaining original weights are renormalized. No accessibility-only rule, because it is constant within group. Comparisons use score, editable-A count, then ID. These are internal decision rules, not predictive efficacy evaluation.

Eight selected records contain six distinct sequences; source IDs and per-variant registers must remain separate. No clinical metrics, biological efficacy or therapeutic claims are supported.
