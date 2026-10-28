"""Run focused synthetic tests; write a fresh evidence directory only."""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import platform
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

class Results(unittest.TextTestResult):
    records = None
    def startTestRun(self):
        self.records = []
        super().startTestRun()
    def addSuccess(self, test):
        super().addSuccess(test); self.records.append((test.id(), 'PASS', ''))
    def addFailure(self, test, err):
        super().addFailure(test, err); self.records.append((test.id(), 'FAIL', str(err[1])))
    def addError(self, test, err):
        super().addError(test, err); self.records.append((test.id(), 'ERROR', str(err[1])))
    def addSkip(self, test, reason):
        super().addSkip(test, reason); self.records.append((test.id(), 'SKIP', reason))

def main():
    destination = HERE / ('test_evidence_' + dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
    destination.mkdir(exist_ok=False)
    suite = unittest.defaultTestLoader.discover(str(HERE), pattern='test_ablation.py')
    with (destination/'TEST_LOG.txt').open('x', encoding='utf-8') as log:
        result = unittest.TextTestRunner(stream=log, verbosity=2, resultclass=Results).run(suite)
    with (destination/'TEST_RESULTS.csv').open('x', encoding='utf-8', newline='') as f:
        w = csv.writer(f); w.writerow(['test', 'status', 'detail']); w.writerows(result.records)
    summary = {'tests': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
               'skipped': len(result.skipped), 'passed': sum(r[1] == 'PASS' for r in result.records),
               'exit_code': 0 if result.wasSuccessful() else 1,
               'python': sys.version, 'platform': platform.platform(),
               'synthetic_only': True, 'real_cohort_scoring': False, 'real_comparisons_executed': False, 'robustness_simulations': False,
               'production_host_enforcement_verified': False,
               'adapter_sha256': hashlib.sha256((HERE/'mecp2_ablation.py').read_bytes()).hexdigest(),
               'command': 'python -I -B run_tests.py', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat()}
    (destination/'TEST_SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'directory':str(destination), **summary}, indent=2))
    return summary['exit_code']

if __name__ == '__main__':
    sys.exit(main())
