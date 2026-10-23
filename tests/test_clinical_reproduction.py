"""A source-to-result run must isolate outputs without requiring an owner drive."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from reproduce_clinical import validate_layout, validate_stage_state


class ClinicalReproductionLayout(unittest.TestCase):
    def test_separate_portable_roots_are_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            result = validate_layout(base / "outputs", base / "inputs" / "CLVer.zip", base / "code")
            self.assertEqual(result, (base / "outputs").resolve())

    def test_overlapping_repository_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "code"
            for run in (repo, repo / "data", repo.parent):
                with self.subTest(run=run), self.assertRaises(ValueError):
                    validate_layout(run, repo=repo)

    def test_original_input_cannot_hide_in_generated_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            run = base / "outputs"
            with self.assertRaises(ValueError):
                validate_layout(run, run / "input.zip", base / "code")
            copied_source = run / "card003" / "raw" / "CLVerPublicDataset.zip"
            self.assertEqual(validate_layout(run, copied_source, base / "code"), run.resolve())

    def test_output_cannot_hide_inside_original_input_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            with self.assertRaises(ValueError):
                validate_layout(base / "inputs" / "outputs", base / "inputs" / "CLVer.zip", base / "code")

    def test_new_clinical_component_can_share_run_with_other_studies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "card001").mkdir()
            validate_stage_state(root, "complete")

    def test_source_preparation_preserves_existing_clinical_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            study = root / "card003"
            study.mkdir()
            source = study / "source_reproduction.json"
            source.write_text("preserved", encoding="utf-8")
            for stage in ("prepare", "complete"):
                with self.subTest(stage=stage), self.assertRaises(FileExistsError):
                    validate_stage_state(root, stage)
            self.assertEqual(source.read_text(encoding="utf-8"), "preserved")

    def test_analysis_refuses_an_existing_frozen_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            policy = root / "card003" / "models" / "frozen_policy.json"
            policy.parent.mkdir(parents=True)
            policy.write_text("preserved", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                validate_stage_state(root, "analyze")
            validate_stage_state(root, "verify")
            self.assertEqual(policy.read_text(encoding="utf-8"), "preserved")


if __name__ == "__main__":
    unittest.main()
