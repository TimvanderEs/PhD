import contextlib
import gzip
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import ldref
import make_demo
import build_panel


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.path = make_demo.create(self.root / "input with spaces")
        self.config = ldref.load(self.path)

    def tearDown(self):
        self.temp.cleanup()

    def test_valid_inputs_and_plan(self):
        result = ldref.validate(self.config)
        self.assertEqual(result["traits"]["toy"]["overlap"], 220)
        self.assertTrue(result["warnings"])
        job = next(ldref.jobs(self.config, self.root / "output"))
        self.assertIn("--not-M-5-50", job["argv"])
        self.assertIn("input with spaces", job["argv"][3])

    def test_sample_size_error(self):
        self.config["traits"]["toy"]["max_n"] = 999
        with self.assertRaisesRegex(ValueError, "exceeds"):
            ldref.validate(self.config)

    def test_panel_order_error(self):
        path = self.path.parent / "panel_b/1.l2.ldscore.gz"
        with gzip.open(path, "rt") as handle:
            lines = handle.readlines()
        lines[1], lines[2] = lines[2], lines[1]
        with gzip.open(path, "wt") as handle:
            handle.writelines(lines)
        with self.assertRaisesRegex(ValueError, "ordered"):
            ldref.validate(self.config)

    def test_build_error(self):
        config = json.loads(self.path.read_text())
        config["traits"]["toy"]["build"] = "GRCh37"
        self.path.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError, "Build mismatch"):
            ldref.load(self.path)

    def test_paired_variance(self):
        result = ldref.paired(self.path.parent / "blocks_a.tsv", self.path.parent / "blocks_b.tsv", .2, .19)
        self.assertAlmostEqual(result["difference"], .01)
        self.assertAlmostEqual(result["paired_se"], math.sqrt(.75*.0002))

    def test_invalid_pairing(self):
        path = self.path.parent / "blocks_b.tsv"
        lines = path.read_text().splitlines()
        lines[1], lines[2] = lines[2], lines[1]
        path.write_text("\n".join(lines)+"\n")
        with self.assertRaisesRegex(ValueError, "order differs"):
            ldref.paired(self.path.parent / "blocks_a.tsv", path, .2, .19)

    def test_run_with_fake_engine_and_overwrite_protection(self):
        # Deliberately fake backend tests orchestration only, never LDSC science.
        engine = self.root / "fake_engine.py"
        engine.write_text("import sys, pathlib\np = sys.argv[sys.argv.index('--out')+1]\npathlib.Path(p+'.log').write_text('Total Observed scale h2: 0.2 (0.03)\\n')\n")
        config = json.loads(self.path.read_text())
        config["engine"] = dict(python=sys.executable, script=str(engine), version="fake-test-only")
        self.path.write_text(json.dumps(config))
        command = [sys.executable, str(Path(ldref.__file__)), "run", str(self.path), "--output", str(self.root / "results")]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = json.loads((self.root / "results/run.json").read_text())
        self.assertEqual(state["status"], "complete")
        self.assertTrue((self.root / "results/h2_summary.tsv").exists())
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)

    def test_failure_recorded(self):
        engine = self.root / "failed_engine.py"
        engine.write_text("raise SystemExit(9)\n")
        config = json.loads(self.path.read_text())
        config["engine"] = dict(python=sys.executable, script=str(engine), version="fake-test-only")
        self.path.write_text(json.dumps(config))
        command = [sys.executable, str(Path(ldref.__file__)), "run", str(self.path), "--output", str(self.root / "failed")]
        self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
        state = json.loads((self.root / "failed/run.json").read_text())
        self.assertEqual(state["status"], "failed")
        self.assertEqual(state["jobs"][0]["returncode"], 9)

    def test_archived_rg_parser(self):
        source = Path(__file__).resolve().parents[1] / "validation/ancestry_matched_sg100k/rg/Chinese/c1_sbp_avg__c8_bmi.chinese_matched.log"
        if not source.exists():
            self.skipTest("Archived Project4 log not included in standalone copy")
        output = self.root / "archive_summary"
        output.mkdir()
        prefix = output / "pair"
        Path(str(prefix)+".log").write_text(source.read_text())
        ldref.write_json(output / "run.json", dict(status="complete", jobs=[dict(analysis="rg", panel="Chinese", traits=["sbp", "bmi"], prefix=str(prefix))]))
        ldref.summarise(output)
        row = list(ldref.rows(output / "rg_summary.tsv"))[0]
        self.assertAlmostEqual(float(row["rg"]), .3137)
        self.assertAlmostEqual(float(row["se"]), .0753)

    def build_fixture(self):
        directory = self.root / "genotypes"
        directory.mkdir()
        for chrom in range(1, 23):
            (directory / f"{chrom}.fam").write_text("f a 0 0 0 -9\nf b 0 0 0 -9\nf c 0 0 0 -9\nf d 0 0 0 -9\n")
            (directory / f"{chrom}.bim").write_text(f"{chrom} rs1 0 100 A C\n{chrom} rs2 1 200 G T\n")
            (directory / f"{chrom}.bed").write_bytes(b"\x6c\x1b\x01\x00\x00")
        (directory / "pcs.txt").write_text("f a -0.5\nf b -0.5\nf c 0.5\nf d 0.5\n")
        config = dict(schema_version=1, build="GRCh38", ancestry="synthetic", source="fake genotype fixture", qc_complete=True,
                      bfile="genotypes/{chr}", covariates="genotypes/pcs.txt", window_cm=1,
                      engine=dict(python=sys.executable, script="absent_engine.py", version="test"))
        path = self.root / "build.json"
        path.write_text(json.dumps(config))
        return path

    def test_builder_plan(self):
        state = build_panel.prepare(self.build_fixture(), self.root / "build_output")
        self.assertEqual(len(state["commands"]), 22)
        self.assertIn("--cov", state["commands"][0])

    def test_builder_covariate_order(self):
        path = self.build_fixture()
        cov = self.root / "genotypes/pcs.txt"
        cov.write_text("f b -0.5\nf a -0.5\nf c 0.5\nf d 0.5\n")
        with self.assertRaisesRegex(ValueError, "order"):
            build_panel.prepare(path, self.root / "bad_build")

    def test_builder_rejects_nonorthonormal_pcs(self):
        path = self.build_fixture()
        (self.root / "genotypes/pcs.txt").write_text("f a -1\nf b -1\nf c 1\nf d 1\n")
        with self.assertRaisesRegex(ValueError, "orthonormal"):
            build_panel.prepare(path, self.root / "bad_build")


if __name__ == "__main__":
    unittest.main()
