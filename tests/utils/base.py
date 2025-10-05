import os
import unittest
from unittest.mock import patch
from tests.utils import _mpl_test_setup  # noqa: F401
from matplotlib.testing.compare import compare_images
from module.core.FileSystem import FileSystem


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
        FileSystem.delete_project(cls.project_name)

    def assert_image_similar(self, expected_path, actual_path, tol=2.0):
        """
        Assert two PNGs are visually similar using Matplotlib's standard comparator.

        Args:
            expected_path: baseline image path (PNG)
            actual_path:   test output image path (PNG)
            tol: RMS tolerance in pixel intensity (float).
                 0 means identical, larger tolerates minor AA/font rasterization diffs.
        """
        self.assertTrue(
            os.path.exists(expected_path), f"Expected not found: {expected_path}"
        )
        self.assertTrue(os.path.exists(actual_path), f"Actual not found: {actual_path}")

        res = compare_images(expected_path, actual_path, tol=tol)
        if res is not None:
            self.fail(res)
