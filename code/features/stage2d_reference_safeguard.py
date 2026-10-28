"""Run-completeness reporting only. No BLAST, scoring or sequence algorithms.

The JSON is the authoritative final record. CSV and JSON carry the same run_id.
No claim about unsearched loci, reference identity or biological outcomes.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import runpy
import re
import sys
import traceback
import uuid

WORKER_NAME = "_mecp2_stage2d_guarded_worker"
EXIT_CODES = {"COMPLETE": 0, "INCOMPLETE_REFERENCE": 2, "BLOCKED_DEPENDENCY": 3,
              "FAILED_INPUT": 4, "FAILED_INTENDED_QC": 5, "FAILED_RUNTIME": 6,
              "INTERRUPTED": 7}
_active = None


class StopRun(RuntimeError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


class IntendedHitQCError(StopRun):
    def __init__(self, message):
        if _active is not None:
            _active["intended_hit_qc"] = "FAIL"
        super().__init__("FAILED_INTENDED_QC", message)


def state():
    if _active is None:
        raise RuntimeError("Safeguard entry point is required; do not import/run the worker directly.")
    return _active


def require_dependencies(required):
    """Check explicitly required files only. No discovery, download or DB launch."""
    if not isinstance(required, dict) or not required:
        raise ValueError("Required dependency manifest must be a nonempty mapping.")
    for name, path in required.items():
        if not isinstance(name, str) or not name or not isinstance(path, (str, Path)) or not str(path):
            raise ValueError("Malformed required dependency manifest.")
        if "database" in name.lower():
            state()["blast_database_requirement"] = "REQUIRED_BY_CALLER"
        if not Path(path).is_file():
            raise StopRun("BLOCKED_DEPENDENCY", f"Required dependency unavailable: {name}: {path}")


def validate_subject_manifest(raw_subjects, normalize):
    """Repeated raw HSP subjects are normal; validate each, then deduplicate."""
    if not isinstance(raw_subjects, (list, tuple)) or not raw_subjects:
        raise ValueError("Subject manifest must contain at least one raw-hit subject.")
    normalized = []
    for index, value in enumerate(raw_subjects, 1):
        if not isinstance(value, str) or not value.strip() or any(c in value for c in '\r\n\t'):
            raise ValueError(f"Missing/malformed raw subject identifier at row {index}.")
        accession = normalize(value)
        if not isinstance(accession, str) or re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.|:-]*', accession) is None:
            raise ValueError(f"Malformed normalized subject at row {index}.")
        normalized.append(accession)
    subjects = sorted(set(normalized))
    s = state()
    s["subject_manifest_rows"] = len(raw_subjects)
    s["expected_subject_count"] = len(subjects)
    s["subject_manifest_valid"] = True
    s["subject_manifest_sha256"] = hashlib.sha256(json.dumps(subjects).encode()).hexdigest()
    return subjects


def require_reference_coverage(subjects, sequences):
    """Fail before reconstruction or result writes on ANY missing subject."""
    s = state()
    if (not s["subject_manifest_valid"] or not subjects or len(subjects) != len(set(subjects))
            or len(subjects) != s["expected_subject_count"] or not isinstance(sequences, dict)):
        raise ValueError("Missing, empty or inconsistent validated subject manifest.")
    if hashlib.sha256(json.dumps(subjects).encode()).hexdigest() != s["subject_manifest_sha256"]:
        raise ValueError("Subject manifest changed after validation.")
    missing = [subject for subject in subjects if subject not in sequences or not sequences[subject]]
    s["missing_subjects"] = missing
    s["missing_subject_count"] = len(missing)
    s["available_subject_count"] = len(subjects) - len(missing)
    s["reference_coverage"] = "INCOMPLETE" if missing else "PASS"
    if missing:
        raise StopRun("INCOMPLETE_REFERENCE",
                      f"Required reference coverage incomplete: {len(missing)}/{len(subjects)} subjects unavailable. "
                      "Reconstruction and intended-hit QC were not run; no biological-invalidity conclusion.")
    return missing


def mark_intended_qc_passed():
    s = state()
    if s["reference_coverage"] != "PASS":
        raise StopRun("INCOMPLETE_REFERENCE", "Cannot accept QC without complete required subject coverage.")
    s["intended_hit_qc"] = "PASS"


def assert_ready_for_results():
    s = state()
    if s["reference_coverage"] != "PASS" or s["missing_subject_count"] != 0:
        raise StopRun("INCOMPLETE_REFERENCE", "Result completion denied: reference coverage not established.")
    if s["intended_hit_qc"] != "PASS":
        raise IntendedHitQCError("Result completion denied: intended-hit QC has not passed.")


def _persist(report, status_path, summary_path):
    """Prepare both files, replace CSV first, authoritative JSON last."""
    status_path, summary_path = Path(status_path), Path(summary_path)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    csv_tmp = summary_path.with_name(summary_path.name + "." + report["run_id"] + ".tmp")
    json_tmp = status_path.with_name(status_path.name + "." + report["run_id"] + ".tmp")
    with csv_tmp.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(report))
        writer.writeheader()
        writer.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (list, dict, bool)) else v
                         for k, v in report.items()})
    json_tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(csv_tmp, summary_path)
    os.replace(json_tmp, status_path)


def run_guarded(work, status_path, summary_path, stage):
    """Callable boundary used by real entry point and isolated synthetic tests."""
    global _active
    if _active is not None:
        raise RuntimeError("Nested safeguard runs are prohibited.")
    report = {"run_id": str(uuid.uuid4()), "stage": stage, "status": "RUNNING", "exit_code": -1,
              "reason": "Checks have not completed.", "subject_manifest_valid": False,
              "subject_manifest_rows": 0, "expected_subject_count": 0,
              "subject_manifest_sha256": "",
              "available_subject_count": 0, "missing_subject_count": 0, "missing_subjects": [],
              "reference_coverage": "NOT_CHECKED", "intended_hit_qc": "NOT_RUN",
              "blast_database_requirement": "NOT_REQUIRED_FOR_EXISTING_RAW_HIT_RESCORING",
              "scope": "Required subjects in supplied BLAST anchors only; reference identity not established."}
    _active = report
    try:
        _persist(report, status_path, summary_path)
        work()
        assert_ready_for_results()
        status, reason = "COMPLETE", "Required subject coverage and intended-hit QC passed; worker returned successfully."
    except StopRun as exc:
        status, reason = exc.status, str(exc)
    except FileNotFoundError as exc:
        status, reason = "BLOCKED_DEPENDENCY", str(exc)
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        status, reason = "FAILED_INPUT", f"{type(exc).__name__}: {exc}"
    except KeyboardInterrupt:
        status, reason = "INTERRUPTED", "Execution interrupted; not complete."
    except BaseException as exc:
        # Includes SystemExit(0): early returns must never bypass required checks.
        status, reason = "FAILED_RUNTIME", f"{type(exc).__name__}: {exc}"
        traceback.print_exc()
    try:
        report.update(status=status, exit_code=EXIT_CODES[status], reason=reason)
        _persist(report, status_path, summary_path)
        print(f"SAFEGUARD FINAL STATUS: {status}; exit={report['exit_code']}")
        return report["exit_code"]
    except BaseException:
        # Output failure must not yield exit 0; RUNNING is not a success record.
        traceback.print_exc()
        return EXIT_CODES["FAILED_RUNTIME"]
    finally:
        _active = None


def run_script(script, status_path, summary_path):
    return run_guarded(lambda: runpy.run_path(str(script), run_name=WORKER_NAME),
                       status_path, summary_path, Path(script).name)
