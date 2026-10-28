#!/usr/bin/env python3
"""
MECP2 ADAR Dataset 2 — Tier 1
Stage 2D v3: FULL-WINDOW HARMONIZED TRANSCRIPTOME OFF-TARGET ANALYSIS

Purpose
-------
Recompute Tier 1 transcriptome specificity using the same corrected full-window
method used for Dataset 3 / Stage 3D v2.

This script DOES NOT change guide sequences or candidate enumeration.
It only replaces the older near-full BLAST-alignment interpretation with:

    BLAST hit anchor
        -> reconstruct the complete subject window
        -> compare full 21/25/31-nt window to guide-derived transcript-equivalent query
        -> retain full-window sites with <=3 mismatches
        -> separate intended transcript
        -> count potential editable-A off-targets

Important Tier 1 expectation
----------------------------
For Tier 1 exact-restoration guides, the guide-derived transcript-equivalent
query should reproduce the WT intended transcript sequence at the edited site.
Therefore all 376 guides are expected to have an intended full-window
0-mismatch hit.

This is the key QC difference from the nonsense-suppression cohort, where
1-2 intended mismatches were expected.

Inputs
------
02_tier1_pipeline/
    MECP2_ADAR_Dataset2_Tier1_stage2C_bystander.csv
    existing Tier 1 Stage 2D raw BLAST-hit TSV

03_reference/
    GRCh38.p14_RefSeq_rna.fna

The script auto-detects a likely Tier 1 raw BLAST hit file under
02_tier1_pipeline. If it cannot find one, it prints the TSV files present.

Outputs
-------
02_tier1_pipeline/
    MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv
    MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_summary.csv
    stage2D_v3_reconstructed_hits.tsv

Methodological note
-------------------
This full-window reconstruction corrects BLAST local-end trimming.
It still inherits candidate-locus discovery sensitivity from the original
BLAST search used to create the raw hits.

Interpretation
--------------
"Zero editable-A off-targets" means no candidate editable-A off-target sites
were detected within the defined BLAST-anchored, reconstructed full-window,
<=3-mismatch RefSeq transcriptome search space. It does NOT mean biological
absence of off-target editing.
"""

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPE = ROOT / "02_tier1_pipeline"
REF = ROOT / "03_reference"

INPUT = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2C_bystander.csv"
FASTA = REF / "GRCh38.p14_RefSeq_rna.fna"

OUT_MAIN = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_offtarget.csv"
OUT_SUM = PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2D_v3_fullwindow_summary.csv"
OUT_HITS = PIPE / "stage2D_v3_reconstructed_hits.tsv"

EXPECTED_GUIDES = 376
EXPECTED_MUTATIONS = 8
MAX_MM = 3

# Separate proposal only. The guard executes this unchanged processing module
# as a worker, catches failures and commits status only after successful return.
from stage2d_reference_safeguard import (
    WORKER_NAME, run_script, validate_subject_manifest, require_reference_coverage,
    IntendedHitQCError, mark_intended_qc_passed, assert_ready_for_results,
)
if __name__ == "__main__":
    raise SystemExit(run_script(
        __file__, PIPE / "stage2D_v4_safeguard_run_status.json",
        PIPE / "stage2D_v4_safeguard_run_summary.csv"))
if __name__ != WORKER_NAME:
    raise RuntimeError("Use the guarded script entry point; direct worker imports are prohibited.")

print("SCRIPT VERSION: TIER1_STAGE2D_V4_ACCESSION_FIX")


def read_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        return list(rd), list(rd.fieldnames or [])


def write_csv(path, rows, fields):
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)


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


def revcomp_rna(seq):
    seq = seq.upper().replace("T", "U")
    comp = str.maketrans("AUGC", "UACG")
    return seq.translate(comp)[::-1]


def rna_to_dna(seq):
    return seq.upper().replace("U", "T")


def dna_to_rna(seq):
    return seq.upper().replace("T", "U")


def clean_accession(s):
    """
    Normalize RefSeq accession identifiers from BLAST and FASTA.

    Examples:
        ref|NM_000016.6|   -> NM_000016.6
        NM_000016.6        -> NM_000016.6
        gi|...|ref|NM_...| -> NM_....
        >NM_000016.6 ...   -> NM_000016.6

    We intentionally preserve the version suffix because the transcriptome
    FASTA and intended transcript annotations are versioned.
    """
    x = str(s or "").strip()
    if not x:
        return ""

    if x.startswith(">"):
        x = x[1:].strip()

    # First token only; description follows whitespace in FASTA headers.
    x = x.split()[0]

    # BLAST DB created from RefSeq FASTA commonly returns:
    #   ref|NM_000016.6|
    # or longer pipe-delimited identifiers.
    parts = [p for p in x.split("|") if p]

    for p in parts:
        if p.startswith(("NM_", "NR_", "XM_", "XR_")):
            return p

    # If there was no RefSeq-looking pipe component, use the token itself.
    return x.strip("|")


def detect_raw_hits():
    candidates = [
        PIPE / "stage2D_blast_hits.tsv",
        PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2D_blast_hits.tsv",
        PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2D_v2_blast_hits.tsv",
        PIPE / "stage2D_v2_blast_hits.tsv",
        PIPE / "MECP2_ADAR_Dataset2_Tier1_stage2D_raw_hits.tsv",
    ]

    for p in candidates:
        if p.exists():
            return p

    # Search all TSVs, favor names containing stage2D + blast/hit.
    tsvs = sorted(PIPE.glob("*.tsv"))
    ranked = []
    for p in tsvs:
        n = p.name.lower()
        score = 0
        if "stage2d" in n:
            score += 10
        if "blast" in n:
            score += 5
        if "hit" in n:
            score += 3
        if "reconstructed" in n or "v3" in n:
            score -= 20
        ranked.append((score, p))

    ranked.sort(key=lambda x: (-x[0], x[1].name.lower()))
    if ranked and ranked[0][0] > 0:
        return ranked[0][1]

    print("\nCould not auto-detect the Tier 1 raw BLAST-hit TSV.")
    print("TSV files currently present in 02_tier1_pipeline:")
    for p in tsvs:
        print("  ", p.name)
    print(
        "\nRename/copy the old Tier 1 raw BLAST file to:\n"
        "  02_tier1_pipeline\\stage2D_blast_hits.tsv\n"
        "and rerun."
    )
    raise FileNotFoundError("Required Tier-1 raw BLAST-hit file unavailable.")


def parse_blast_hits(path):
    """
    Supports either:
    A) headered BLAST TSV
    B) standard custom outfmt rows with these columns:
       qseqid sseqid pident length mismatch gapopen
       qstart qend sstart send evalue bitscore qseq sseq
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        first = f.readline().rstrip("\n\r")

    parts = first.split("\t")
    lower = [x.strip().lower() for x in parts]

    known = {"qseqid", "sseqid", "qstart", "qend", "sstart", "send"}
    has_header = len(known.intersection(lower)) >= 4

    rows = []

    if has_header:
        with open(path, "r", encoding="utf-8-sig", newline="") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for r in rd:
                rows.append(r)
        return rows, list(rows[0].keys()) if rows else lower

    default_fields = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qseq", "sseq",
    ]

    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rd = csv.reader(f, delimiter="\t")
        for vals in rd:
            if not vals or len(vals) < 10:
                continue
            d = {}
            for i, v in enumerate(vals):
                if i < len(default_fields):
                    d[default_fields[i]] = v
                else:
                    d[f"extra_{i+1}"] = v
            rows.append(d)

    return rows, default_fields


def load_selected_fasta(path, wanted):
    """
    Stream the transcriptome FASTA and retain only accessions that occur
    in the BLAST hit table.
    """
    wanted = set(wanted)
    seqs = {}

    current = None
    chunks = []

    def commit():
        nonlocal current, chunks
        if current in wanted:
            seqs[current] = "".join(chunks).upper()

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if current is not None:
                    commit()
                header = line[1:].strip()
                current = clean_accession(header)
                chunks = []
            else:
                if current in wanted:
                    chunks.append(line.strip())

        if current is not None:
            commit()

    return seqs


# ---------------------------------------------------------------------
# Load guide table.
# ---------------------------------------------------------------------
if not INPUT.exists():
    raise FileNotFoundError(f"Tier 1 Stage 2C input not found: {INPUT}")

if not FASTA.exists():
    raise FileNotFoundError(f"RefSeq transcriptome FASTA not found: {FASTA}")

rows, fields = read_csv(INPUT)

if len(rows) != EXPECTED_GUIDES:
    raise ValueError(f"Expected {EXPECTED_GUIDES} Tier 1 guides, found {len(rows)}.")

if len({r["Guide_ID"] for r in rows}) != EXPECTED_GUIDES:
    raise ValueError("Duplicate Guide_ID values detected.")

if len({r["VariationID"] for r in rows}) != EXPECTED_MUTATIONS:
    raise ValueError(
        f"Expected {EXPECTED_MUTATIONS} Tier 1 mutations, "
        f"found {len({r['VariationID'] for r in rows})}."
    )


# ---------------------------------------------------------------------
# Resolve the ACTUAL guide-sequence column.
# IMPORTANT: do not call find_col(... tokens=None) here because that helper
# treats an empty token list as "match anything", which can select Guide_ID.
# ---------------------------------------------------------------------
KNOWN_GUIDE_COLUMNS = [
    "Guide_RNA",
    "Guide_sequence",
    "Guide_RNA_sequence",
    "Guide_sequence_RNA",
    "Guide_seq",
    "Guide",
    "gRNA_sequence",
    "Guide_with_AC_mismatch",
    "Guide_RNA_with_AC_mismatch",
    "Final_guide_RNA",
    "Candidate_guide_RNA",
]

guide_col = next((c for c in KNOWN_GUIDE_COLUMNS if c in fields), None)

def sequence_column_score(col):
    vals = [
        str(rr.get(col, "") or "").strip().upper().replace(" ", "")
        for rr in rows[:min(150, len(rows))]
    ]
    vals = [v for v in vals if v]
    if not vals:
        return None

    valid_alpha = sum(all(ch in "AUGCT" for ch in v) for v in vals) / len(vals)
    valid_len = sum(
        all(ch in "AUGCT" for ch in v) and len(v) in (21, 25, 31)
        for v in vals
    ) / len(vals)

    if valid_alpha < 0.95 or valid_len < 0.80:
        return None

    cl = col.lower()
    bonus = 0
    if "guide" in cl or "grna" in cl:
        bonus += 10
    if "sequence" in cl or "seq" in cl:
        bonus += 5
    if "target" in cl:
        bonus -= 2

    return (bonus, valid_len, valid_alpha)

# Even when a known-name column exists, validate that its CONTENT is actually
# 21/25/31-nt nucleotide sequence. This prevents Guide_ID-like mistakes.
if guide_col is not None and sequence_column_score(guide_col) is None:
    guide_col = None

if guide_col is None:
    candidates = []
    for c in fields:
        score = sequence_column_score(c)
        if score is not None:
            candidates.append((*score, c))

    if not candidates:
        print("\nERROR: Could not identify the guide-RNA sequence column.")
        print("Available Stage 2C columns:")
        for c in fields:
            sample = next(
                (str(rr.get(c, "") or "") for rr in rows[:5] if str(rr.get(c, "") or "")),
                ""
            )
            print(f"  {c:45s} sample={sample[:70]}")
        raise KeyError(
            "No column contains predominantly 21/25/31-nt A/U/G/C/T guide sequences."
        )

    candidates.sort(key=lambda x: (-x[0], -x[1], -x[2], len(x[3]), x[3].lower()))
    guide_col = candidates[0][3]
    print("Auto-detected guide sequence column:", guide_col)
    if len(candidates) > 1:
        print("Other sequence-like columns:", [x[3] for x in candidates[1:6]])

# Final hard guard: Guide_ID can never be used as sequence.
if guide_col.lower() in {"guide_id", "guideid", "id"}:
    raise RuntimeError(
        f"Internal column-resolution error: {guide_col} was selected as guide sequence."
    )

guide_len_col = find_col(
    fields,
    exact=["Guide_length_nt", "Guide_length"],
    tokens=["guide", "length"],
    required=False,
)

transcript_col = find_col(
    fields,
    exact=[
        "Transcript",
        "Target_transcript",
        "RefSeq_transcript",
        "Transcript_accession",
        "Target_transcript_accession",
        "RefSeq_transcript_accession",
    ],
    tokens=["transcript"],
    excludes=["offtarget", "specificity", "count"],
    required=False,
)

if not transcript_col:
    # Detect RefSeq transcript accession column from values such as NM_001110792.2.
    for c in fields:
        vals = [str(rr.get(c, "") or "").strip() for rr in rows[:50]]
        vals = [v for v in vals if v]
        if vals and sum(v.startswith(("NM_", "NR_", "XM_", "XR_")) for v in vals) / len(vals) >= 0.8:
            transcript_col = c
            print("Auto-detected transcript column:", transcript_col)
            break

if not transcript_col:
    print("\nCould not identify intended transcript column.")
    print("Available Stage 2C columns:")
    for c in fields:
        print("  ", c)
    raise KeyError("No transcript/accession column could be resolved.")

# Coordinate contract: use an explicitly guide-relative field only.
# Absolute transcript/CDS and generic position columns are never substitutes.
target_pos_col = "Target_A_position_in_guide_1based"
if fields.count(target_pos_col) != 1:
    raise KeyError(
        "Require exactly one Target_A_position_in_guide_1based column; "
        "legacy absolute/generic coordinates are not guide positions."
    )
segment_pos_col = "Target_A_position_in_target_segment_1based"
if fields.count(segment_pos_col) > 1:
    raise KeyError("Duplicate explicit target-segment coordinate column.")

print("Resolved Tier 1 columns:")
print("  Guide sequence :", guide_col)
print("  Guide length   :", guide_len_col or "derived from sequence")
print("  Transcript     :", transcript_col)
print("  Target position:", target_pos_col or "auto-resolved from position metadata")

# Use exact Stage 2A/2C design provenance if available.
guide_by_id = {}
for r in rows:
    gid = r["Guide_ID"]
    guide = str(r[guide_col]).strip().upper().replace("T", "U")
    if not guide:
        raise ValueError(f"Empty guide sequence for {gid}")

    # The transcript-equivalent query is the reverse complement of the
    # ACTUAL guide. Because guide C is opposite catalytic target A,
    # that position appears as G in the query.
    query_rna = revcomp_rna(guide)
    query_dna = rna_to_dna(query_rna)

    g_len = len(guide)
    if guide_len_col:
        stated = int(float(r[guide_len_col]))
        if stated != g_len:
            raise ValueError(
                f"{gid}: stated guide length {stated} != sequence length {g_len}"
            )

    guide_by_id[gid] = {
        "row": r,
        "guide_rna": guide,
        "query_rna": query_rna,
        "query_dna": query_dna,
        "qlen": g_len,
        "intended_transcript": clean_accession(r[transcript_col]),
    }


# ---------------------------------------------------------------------
# Load old raw BLAST anchors.
# ---------------------------------------------------------------------
RAW_HITS = detect_raw_hits()
print("Tier 1 raw BLAST hits detected:", RAW_HITS)

blast_rows, blast_fields = parse_blast_hits(RAW_HITS)

required_blast = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send"]
for c in required_blast:
    if c not in blast_fields:
        raise KeyError(f"Raw BLAST table missing required column: {c}")

print("Raw BLAST hits loaded:", len(blast_rows))

subjects = validate_subject_manifest(
    [r.get("sseqid") for r in blast_rows], clean_accession
)

# Hard guard against accidentally running the pre-fix script.
bad_subjects = [s for s in subjects[:1000] if s.startswith("ref|") or s.endswith("|")]
if bad_subjects:
    raise RuntimeError(
        "ACCESSION NORMALIZATION FIX IS NOT ACTIVE. "
        f"Example unnormalized subject: {bad_subjects[0]}"
    )

print("Unique transcript accessions needed:", len(subjects))
print("First normalized BLAST subject accessions:", subjects[:5])

seqs = load_selected_fasta(FASTA, subjects)

missing_subjects = require_reference_coverage(subjects, seqs)
print("Required transcript sequences loaded:", len(seqs))
print("BLAST subjects absent from FASTA dictionary:", len(missing_subjects))

if missing_subjects:
    print("First missing subjects:", missing_subjects[:20])


# ---------------------------------------------------------------------
# Reconstruct full subject windows from BLAST local anchors.
# ---------------------------------------------------------------------
reconstructed = defaultdict(dict)
total_reconstructed_anchor_loci = 0

for i, h in enumerate(blast_rows, 1):
    gid = str(h.get("qseqid", "")).strip()
    if gid not in guide_by_id:
        continue

    sid = clean_accession(h.get("sseqid", ""))
    if sid not in seqs:
        continue

    try:
        qstart = int(float(h["qstart"]))
        qend = int(float(h["qend"]))
        sstart = int(float(h["sstart"]))
        send = int(float(h["send"]))
    except Exception:
        continue

    # This harmonized analysis only uses plus-strand transcript hits,
    # matching the prior guide-aware Tier 1 approach.
    if send < sstart:
        continue

    qlen = guide_by_id[gid]["qlen"]

    # For a local alignment beginning at query position qstart and subject
    # coordinate sstart, infer where the complete query-length subject window
    # begins.
    full_start = sstart - (qstart - 1)
    full_end = full_start + qlen - 1

    if full_start < 1:
        continue

    sseq = seqs[sid]
    if full_end > len(sseq):
        continue

    subject_window_dna = sseq[full_start - 1:full_end].upper().replace("U", "T")
    if len(subject_window_dna) != qlen:
        continue

    query_dna = guide_by_id[gid]["query_dna"]

    mm_positions = [
        j + 1 for j, (a, b) in enumerate(zip(query_dna, subject_window_dna))
        if a != b
    ]
    mm = len(mm_positions)

    key = (sid, full_start, full_end)

    candidate = {
        "Guide_ID": gid,
        "VariationID": guide_by_id[gid]["row"]["VariationID"],
        "Subject_accession": sid,
        "Full_window_start_1based": full_start,
        "Full_window_end_1based": full_end,
        "Guide_length_nt": qlen,
        "Query_transcript_equivalent_RNA": guide_by_id[gid]["query_rna"],
        "Subject_full_window_RNA": dna_to_rna(subject_window_dna),
        "Full_window_mismatch_count": mm,
        "Mismatch_positions_1based": ",".join(map(str, mm_positions)),
        "Anchor_qstart": qstart,
        "Anchor_qend": qend,
        "Anchor_sstart": sstart,
        "Anchor_send": send,
    }

    # Deduplicate same reconstructed locus. If multiple BLAST local alignments
    # point to the same full window, retain the one with the longest query span.
    span = abs(qend - qstart) + 1
    prev = reconstructed[gid].get(key)
    if prev is None:
        candidate["_anchor_span"] = span
        reconstructed[gid][key] = candidate
        total_reconstructed_anchor_loci += 1
    elif span > prev["_anchor_span"]:
        candidate["_anchor_span"] = span
        reconstructed[gid][key] = candidate


# ---------------------------------------------------------------------
# Resolve the catalytic query coordinate from explicit guide geometry.
# A missing segment column is allowed; if supplied, it must agree exactly.
# ---------------------------------------------------------------------
def catalytic_query_pos_1based(r, guide):
    gid = r.get("Guide_ID", "UNKNOWN")
    guide_field = "Target_A_position_in_guide_1based"
    segment_field = "Target_A_position_in_target_segment_1based"

    def integer_position(name):
        value = r.get(name)
        text = "" if value is None else str(value).strip()
        if not text or any(ch not in "0123456789" for ch in text):
            raise ValueError(f"{gid}: {name} must be an explicit integer coordinate.")
        position = int(text)
        if not 1 <= position <= len(guide):
            raise ValueError(f"{gid}: {name} outside guide-length bounds.")
        return position

    guide_pos = integer_position(guide_field)
    query_pos = len(guide) - guide_pos + 1
    if segment_field in r:
        segment_pos = integer_position(segment_field)
        if segment_pos != query_pos:
            raise ValueError(
                f"{gid}: explicit guide and target-segment positions disagree."
            )
    if guide[guide_pos - 1].upper() != "C":
        raise ValueError(f"{gid}: designated guide base must be C for this A:C design.")
    return query_pos



# ---------------------------------------------------------------------
# Score each guide.
# ---------------------------------------------------------------------
main_rows = []
audit_rows = []
summary_by_mut = defaultdict(lambda: {
    "guides": 0,
    "intended_valid": 0,
    "zero_editable": 0,
    "editable_ge1": 0,
})

all_raw_risks = []

for gid, meta in guide_by_id.items():
    r = dict(meta["row"])
    intended = meta["intended_transcript"]
    qlen = meta["qlen"]

    cat_pos = catalytic_query_pos_1based(r, meta["guide_rna"])
    if not (1 <= cat_pos <= qlen):
        raise ValueError(f"{gid}: catalytic query position out of bounds: {cat_pos}")

    loci = list(reconstructed.get(gid, {}).values())

    valid = [x for x in loci if int(x["Full_window_mismatch_count"]) <= MAX_MM]

    intended_hits = [
        x for x in valid
        if clean_accession(x["Subject_accession"]) == intended
    ]

    non_target_hits = [
        x for x in valid
        if clean_accession(x["Subject_accession"]) != intended
    ]

    best_intended_mm = (
        min(int(x["Full_window_mismatch_count"]) for x in intended_hits)
        if intended_hits else None
    )

    # Editable-A off-target:
    # subject transcript actually has A at catalytic target position.
    editable_hits = []
    for x in non_target_hits:
        s = x["Subject_full_window_RNA"]
        if len(s) >= cat_pos and s[cat_pos - 1] == "A":
            editable_hits.append(x)

    # Counts by mismatch class.
    nt_counts = Counter(
        int(x["Full_window_mismatch_count"]) for x in non_target_hits
    )
    ed_counts = Counter(
        int(x["Full_window_mismatch_count"]) for x in editable_hits
    )

    # Same raw risk formula used in Dataset 3 / Stage 3D v2.
    editable_component = sum(
        1.0 / (1.0 + int(x["Full_window_mismatch_count"]))
        for x in editable_hits
    )

    similarity_component = (
        nt_counts.get(0, 0) * 0.25
        + nt_counts.get(1, 0) * 0.15
        + nt_counts.get(2, 0) * 0.08
        + nt_counts.get(3, 0) * 0.04
    )

    raw_risk = editable_component + similarity_component
    all_raw_risks.append(raw_risk)

    r["Fullwindow_intended_transcript_hit_le3mm"] = "YES" if intended_hits else "NO"
    r["Fullwindow_intended_best_mismatch_count"] = (
        best_intended_mm if best_intended_mm is not None else ""
    )
    r["Fullwindow_non_target_hits_le3mm"] = len(non_target_hits)
    r["Fullwindow_non_target_0_mismatch_hits"] = nt_counts.get(0, 0)
    r["Fullwindow_non_target_1_mismatch_hits"] = nt_counts.get(1, 0)
    r["Fullwindow_non_target_2_mismatch_hits"] = nt_counts.get(2, 0)
    r["Fullwindow_non_target_3_mismatch_hits"] = nt_counts.get(3, 0)

    r["Fullwindow_editable_A_offtarget_hits_le3mm"] = len(editable_hits)
    r["Fullwindow_editable_A_offtarget_0mm"] = ed_counts.get(0, 0)
    r["Fullwindow_editable_A_offtarget_1mm"] = ed_counts.get(1, 0)
    r["Fullwindow_editable_A_offtarget_2mm"] = ed_counts.get(2, 0)
    r["Fullwindow_editable_A_offtarget_3mm"] = ed_counts.get(3, 0)

    r["Fullwindow_GuideAware_raw_offtarget_risk"] = raw_risk
    r["Fullwindow_catalytic_query_position_1based"] = cat_pos
    r["Fullwindow_search_definition"] = (
        "BLAST_ANCHORED_RECONSTRUCTED_FULL_WINDOW_LE3_MISMATCH_REFSEQ_TRANSCRIPTOME"
    )
    r["Fullwindow_offtarget_interpretation"] = (
        "No detected editable-A off-target means none detected within the "
        "defined BLAST-anchored reconstructed full-window <=3-mismatch "
        "RefSeq transcriptome search space."
    )
    r["Stage2D_v3_status"] = "PENDING_COMPLETENESS_AND_INTENDED_QC"

    main_rows.append(r)

    sm = summary_by_mut[r["VariationID"]]
    sm["guides"] += 1
    if intended_hits:
        sm["intended_valid"] += 1
    if len(editable_hits) == 0:
        sm["zero_editable"] += 1
    else:
        sm["editable_ge1"] += 1

    for x in valid:
        ax = dict(x)
        ax.pop("_anchor_span", None)
        ax["Is_intended_transcript"] = (
            "YES" if clean_accession(x["Subject_accession"]) == intended else "NO"
        )
        ax["Is_editable_A_offtarget"] = (
            "YES"
            if (
                clean_accession(x["Subject_accession"]) != intended
                and x["Subject_full_window_RNA"][cat_pos - 1] == "A"
            )
            else "NO"
        )
        ax["Catalytic_query_position_1based"] = cat_pos
        audit_rows.append(ax)


# ---------------------------------------------------------------------
# Normalize raw risk within Tier 1 cohort, exactly as a cohort-relative score.
# ---------------------------------------------------------------------
rmin = min(all_raw_risks)
rmax = max(all_raw_risks)

for r in main_rows:
    raw = float(r["Fullwindow_GuideAware_raw_offtarget_risk"])
    if rmax > rmin:
        norm = (raw - rmin) / (rmax - rmin)
    else:
        norm = 0.0

    specificity = 1.0 - norm

    r["Fullwindow_GuideAware_transcriptome_risk_score_0_1"] = norm
    r["Fullwindow_GuideAware_transcriptome_specificity_score_0_1"] = specificity


# Quartile categories, descriptive only.
sorted_risks = sorted(
    float(r["Fullwindow_GuideAware_transcriptome_risk_score_0_1"])
    for r in main_rows
)

def percentile(vals, p):
    if not vals:
        return 0.0
    pos = (len(vals) - 1) * p
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    return vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)

q25 = percentile(sorted_risks, 0.25)
q75 = percentile(sorted_risks, 0.75)

for r in main_rows:
    v = float(r["Fullwindow_GuideAware_transcriptome_risk_score_0_1"])
    if v <= q25:
        cat = "LOWER_RELATIVE_RISK"
    elif v >= q75:
        cat = "HIGHER_RELATIVE_RISK"
    else:
        cat = "INTERMEDIATE_RELATIVE_RISK"
    r["Fullwindow_GuideAware_relative_risk_category"] = cat


# ---------------------------------------------------------------------
# Tier 1 hard QC.
# Exact-restoration guide query should exactly match intended WT transcript.
# ---------------------------------------------------------------------
intended_valid = sum(
    r["Fullwindow_intended_transcript_hit_le3mm"] == "YES"
    for r in main_rows
)

best_mm_counts = Counter(
    str(r["Fullwindow_intended_best_mismatch_count"])
    for r in main_rows
)

if intended_valid != EXPECTED_GUIDES:
    missing = [
        r["Guide_ID"] for r in main_rows
        if r["Fullwindow_intended_transcript_hit_le3mm"] != "YES"
    ]
    raise IntendedHitQCError(
        f"Tier 1 intended-target QC failed: only {intended_valid}/{EXPECTED_GUIDES} "
        f"guides have a valid reconstructed intended hit <=3 mismatches. "
        f"First missing: {missing[:20]}"
    )

# Tier 1 should be exact 0-mm at intended transcript after applying
# guide-derived transcript-equivalent query.
nonzero_intended = [
    (r["Guide_ID"], r["Fullwindow_intended_best_mismatch_count"])
    for r in main_rows
    if str(r["Fullwindow_intended_best_mismatch_count"]) != "0"
]

if nonzero_intended:
    raise IntendedHitQCError(
        "Tier 1 exact-restoration QC failed: expected every guide-derived "
        "transcript-equivalent query to have a 0-mismatch intended WT hit. "
        f"First non-zero cases: {nonzero_intended[:20]}"
    )


# ---------------------------------------------------------------------
# Write outputs.
# ---------------------------------------------------------------------
mark_intended_qc_passed()
assert_ready_for_results()
for checked_row in main_rows:
    checked_row["Stage2D_v3_status"] = "FULLWINDOW_HARMONIZED_OFFTARGET_ANALYSIS_COMPLETE"
new_fields = [
    "Fullwindow_intended_transcript_hit_le3mm",
    "Fullwindow_intended_best_mismatch_count",
    "Fullwindow_non_target_hits_le3mm",
    "Fullwindow_non_target_0_mismatch_hits",
    "Fullwindow_non_target_1_mismatch_hits",
    "Fullwindow_non_target_2_mismatch_hits",
    "Fullwindow_non_target_3_mismatch_hits",
    "Fullwindow_editable_A_offtarget_hits_le3mm",
    "Fullwindow_editable_A_offtarget_0mm",
    "Fullwindow_editable_A_offtarget_1mm",
    "Fullwindow_editable_A_offtarget_2mm",
    "Fullwindow_editable_A_offtarget_3mm",
    "Fullwindow_GuideAware_raw_offtarget_risk",
    "Fullwindow_GuideAware_transcriptome_risk_score_0_1",
    "Fullwindow_GuideAware_transcriptome_specificity_score_0_1",
    "Fullwindow_GuideAware_relative_risk_category",
    "Fullwindow_catalytic_query_position_1based",
    "Fullwindow_search_definition",
    "Fullwindow_offtarget_interpretation",
    "Stage2D_v3_status",
]

out_fields = list(fields)
for c in new_fields:
    if c not in out_fields:
        out_fields.append(c)

write_csv(OUT_MAIN, main_rows, out_fields)

summary_fields = [
    "VariationID",
    "Guide_count",
    "Valid_intended_fullwindow_le3mm",
    "Zero_editable_A_offtarget_guides",
    "Editable_A_offtarget_ge1_guides",
]

summary_rows = []
for vid, sm in sorted(summary_by_mut.items()):
    summary_rows.append({
        "VariationID": vid,
        "Guide_count": sm["guides"],
        "Valid_intended_fullwindow_le3mm": sm["intended_valid"],
        "Zero_editable_A_offtarget_guides": sm["zero_editable"],
        "Editable_A_offtarget_ge1_guides": sm["editable_ge1"],
    })

write_csv(OUT_SUM, summary_rows, summary_fields)

audit_fields = [
    "Guide_ID",
    "VariationID",
    "Subject_accession",
    "Full_window_start_1based",
    "Full_window_end_1based",
    "Guide_length_nt",
    "Query_transcript_equivalent_RNA",
    "Subject_full_window_RNA",
    "Full_window_mismatch_count",
    "Mismatch_positions_1based",
    "Anchor_qstart",
    "Anchor_qend",
    "Anchor_sstart",
    "Anchor_send",
    "Is_intended_transcript",
    "Is_editable_A_offtarget",
    "Catalytic_query_position_1based",
]

with open(OUT_HITS, "w", encoding="utf-8-sig", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=audit_fields, delimiter="\t", extrasaction="ignore")
    wr.writeheader()
    wr.writerows(audit_rows)


# ---------------------------------------------------------------------
# Console report.
# ---------------------------------------------------------------------
editable_ge1 = sum(
    int(r["Fullwindow_editable_A_offtarget_hits_le3mm"]) >= 1
    for r in main_rows
)
editable_zero = EXPECTED_GUIDES - editable_ge1

risk_categories = Counter(
    r["Fullwindow_GuideAware_relative_risk_category"]
    for r in main_rows
)

print("\n" + "=" * 90)
print("STAGE 2D v3 COMPLETE — TIER 1 FULL-WINDOW HARMONIZED OFF-TARGET ANALYSIS")
print("=" * 90)
print("Input guides                                  :", len(main_rows))
print("Mutations                                     :", len(summary_by_mut))
print("Guides with valid intended full-window <=3mm  :", intended_valid)
print("Guides without valid intended full-window     :", EXPECTED_GUIDES - intended_valid)
print("Best intended full-window mismatch counts     :", dict(best_mm_counts))
print("Guides with editable-A off-target >=1         :", editable_ge1)
print("Guides with zero editable-A off-targets       :", editable_zero)
print(
    "Corrected raw risk range                      :",
    min(all_raw_risks), "to", max(all_raw_risks)
)
print("Relative categories                           :", risk_categories)
print("Reconstructed BLAST-anchor loci               :", total_reconstructed_anchor_loci)
print("BLAST subjects missing from FASTA             :", len(missing_subjects))
print()
print("Main output:")
print(OUT_MAIN)
print("Summary:")
print(OUT_SUM)
print("Reconstructed-hit audit:")
print(OUT_HITS)
print()
print("INTENDED-TARGET QC PASSED:")
print("All 376 Tier 1 guides have an exact 0-mismatch reconstructed")
print("intended WT transcript full-window hit.")
print()
print(
    "Interpretation: specificity is a Tier 1 within-cohort computational "
    "metric. Zero detected editable-A off-targets means none detected "
    "within the defined BLAST-anchored reconstructed full-window "
    "<=3-mismatch RefSeq transcriptome search space."
)
