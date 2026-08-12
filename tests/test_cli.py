from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    def test_reproduce_creates_manifest_and_all_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "result"
            environment = dict(os.environ)
            environment["PYTHONPATH"] = str(ROOT / "src")
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "provfold.cli",
                    "reproduce",
                    "--recipe",
                    str(ROOT / "examples/minimal/recipe.json"),
                    "--output-dir",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
                env=environment,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((output / "run_manifest.json").exists())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["software"], "provfold")
            self.assertGreaterEqual(len(manifest["outputs"]), 10)
            self.assertTrue((output / "aggregate/taxon_summaries.csv").exists())
            self.assertTrue((output / "compare/method_comparison.csv").exists())
            self.assertTrue((output / "sensitivity/threshold_surface.csv").exists())
            self.assertTrue((output / "omit_family/leave_one_family_out.csv").exists())


if __name__ == "__main__":
    unittest.main()
