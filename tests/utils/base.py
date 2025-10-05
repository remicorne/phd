import os
import shutil
import unittest
from unittest.mock import patch
from matplotlib.testing.compare import compare_images
from pathlib import Path


class PlotterTestCase(unittest.TestCase):
    """
    Base test case for plotter functional tests.

    Centralizes:
    - stable Matplotlib/Seaborn/Numpy setup (via import side-effect)
    - project/dataset defaults
    - common patches for interactive prompts and file edits
    - cleanup and handy assertions
    """

    project_name = "tcb2_test_project"
    dataset_name = "hplc"
    ARTIFACT_ROOT = Path("test_artifacts") / "result_images"

    @classmethod
    def setUpClass(cls):
        cls.p_yes_meta = patch("module.core.Metadata.yes_or_no", return_value=True)
        cls.p_yes_ds = patch("module.core.ProjectDataset.yes_or_no", return_value=True)
        cls.p_input = patch(
            "module.core.ProjectDataset.input_escape",
            side_effect=[
                "module/example_project/tcb2_hplc_data.csv",
                "module/example_project/tcb2_fake_behavior_data.csv",
            ],
        )
        cls.p_edit = patch(
            "module.core.Metadata.ProjectMetadata.make_user_edit_excel",
            return_value=None,
        )

        cls.m_yes_meta = cls.p_yes_meta.start()
        cls.m_yes_ds = cls.p_yes_ds.start()
        cls.m_input = cls.p_input.start()
        cls.m_edit = cls.p_edit.start()

    @classmethod
    def tearDownClass(cls):
        cls.p_yes_meta.stop()
        cls.p_yes_ds.stop()
        cls.p_input.stop()
        cls.p_edit.stop()
        # FileSystem.delete_project(cls.project_name)

    def assert_image_similar(self, expected_path, actual_path, tol=2.0):
        self.assertTrue(
            os.path.exists(expected_path), f"Expected not found: {expected_path}"
        )
        self.assertTrue(os.path.exists(actual_path), f"Actual not found: {actual_path}")

        res = compare_images(expected_path, actual_path, tol=tol)
        if res is not None:
            # Build a per-test folder like test_artifacts/result_images/package_Class_test_method/
            safe_name = self.id().replace(".", "_")
            outdir = self.ARTIFACT_ROOT / safe_name
            outdir.mkdir(parents=True, exist_ok=True)

            # Copy the three images (expected, actual, diff) if present
            try:
                shutil.copyfile(expected_path, outdir / "expected.png")
            except Exception:
                pass
            try:
                shutil.copyfile(actual_path, outdir / "actual.png")
            except Exception:
                pass

            # Also write a tiny summary (RMS, message)
            summary = outdir / "summary.txt"
            with summary.open("w", encoding="utf-8") as fh:
                fh.write(f"RMS: {res.get('rms')}\n")
                fh.write(f"Message: {res.get('msg')}\n")
                fh.write(f"Expected: {expected_path}\nActual:   {actual_path}\n")

            # Fail with a pointer to the artifact folder
            self.fail(
                f"Image files did not match (RMS={res.get('rms')}). "
                f"See artifacts in: {outdir.as_posix()}"
            )
