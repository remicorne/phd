import pandas as pd
from tests.utils.base import PlotterTestCase
from module.core.ProjectDataset import ProjectDataset
from module.core.enums import DatasetColumn

EXPERIMENT = "agonist_antagonist"


class TestProjectDataset(PlotterTestCase):
    """Unit tests for module.core.ProjectDataset using example project data.

    Inherits PlotterTestCase to reuse stable environment and patched IO.
    """

    def make_ds(self) -> ProjectDataset:
        return ProjectDataset(project=self.project_name, filename=self.dataset_name)

    def test_stats_with_outliers(self):
        ds = self.make_ds().select(
            experiment=EXPERIMENT,
            compound="DA",
            region=["OF"],
        )
        group_stats = ds.group_statistics.iloc[0]
        pd.testing.assert_series_equal(
            group_stats,
            pd.Series(
                {
                    DatasetColumn.GROUP_ID: 1,
                    "compound": "DA",
                    "region": "OF",
                    "shapiro_F": 0.8205090372142911,
                    "shapiro_p": 0.017515024431137262,
                    "is_parametric": False,
                    "mean": 0.019056090909090908,
                    "std": 0.014390968671041886,
                    "sem": 0.0043390403137823115,
                    "values": [
                        0.012725,
                        0.056825,
                        0.02732,
                        0.019416,
                        0.009179,
                        0.009371,
                        0.021207,
                        0.016852,
                        0.015104,
                        0.00123,
                        0.020388,
                    ],
                },
            ),
            check_names=False,
        )
        quantitative_stats = ds.quantitative_statistics_table.copy()
        cols_without_result = [col for col in quantitative_stats if col != "result"]
        pd.testing.assert_frame_equal(
            quantitative_stats[cols_without_result],
            pd.DataFrame(
                [
                    {
                        "fully_significant": False,
                        "region": "OF",
                        "test": "two_way_anova",
                        "p_value": 0.505,
                        "result_string": "F(1, 36) = 0.453, p = 0.505 ",
                        "compound": "DA",
                        "is_significant": False,
                        "p_value_threshold": 0.05,
                        "experiment": "agonist_antagonist",
                    },
                    {
                        "fully_significant": False,
                        "region": "OF",
                        "test": "one_way_anova",
                        "p_value": 0.5452311344861052,
                        "result_string": "F(3, 36) = 0.722, p = 0.55 ",
                        "compound": "DA",
                        "is_significant": False,
                        "p_value_threshold": 0.05,
                        "experiment": "agonist_antagonist",
                    },
                    {
                        "fully_significant": False,
                        "region": "OF",
                        "test": "tukey",
                        "p_value": [[], []],
                        "result_string": pd.DataFrame(columns=["pairs", "p_values"]),
                        "compound": "DA",
                        "is_significant": False,
                        "p_value_threshold": 0.05,
                        "experiment": "agonist_antagonist",
                    },
                ],
            ),
            check_like=True,
        )

    def test_stats_without_grubbs_outliers(self):
        ds = self.make_ds().select(
            experiment=EXPERIMENT,
            compound="DA",
            region=["OF"],
            remove_outliers={"grubbs": "calculated"},
        )
        group_stats = ds.group_statistics.iloc[0]
        pd.testing.assert_series_equal(
            group_stats,
            pd.Series(
                {
                    DatasetColumn.GROUP_ID: 1,
                    "compound": "DA",
                    "region": "OF",
                    "shapiro_F": 0.9801003059192006,
                    "shapiro_p": 0.9657026230138452,
                    "is_parametric": True,
                    "mean": 0.015279200000000001,
                    "std": 0.007467412534175111,
                    "sem": 0.0023614031836083297,
                    "values": [
                        0.012725,
                        0.02732,
                        0.019416,
                        0.009179,
                        0.009371,
                        0.021207,
                        0.016852,
                        0.015104,
                        0.00123,
                        0.020388,
                    ],
                },
            ),
            check_names=False,
        )
        quantitative_stats = ds.quantitative_statistics_table.copy()
        cols_without_result = [col for col in quantitative_stats if col != "result"]
        pd.testing.assert_frame_equal(
            quantitative_stats[cols_without_result],
            pd.DataFrame(
                [
                    {
                        "p_value": 0.576,
                        "is_significant": False,
                        "result_string": "F(1, 33) = 0.318, p = 0.576 ",
                        "test": "two_way_anova",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                    {
                        "p_value": 0.17304417123590965,
                        "is_significant": False,
                        "result_string": "F(3, 33) = 1.76, p = 0.17 ",
                        "test": "one_way_anova",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                    {
                        "p_value": [[], []],
                        "is_significant": False,
                        "result_string": pd.DataFrame(columns=["pairs", "p_values"]),
                        "test": "tukey",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                ],
            ),
            check_like=True,
        )

    def test_stats_without_iqr_outliers(self):
        ds = self.make_ds().select(
            experiment=EXPERIMENT,
            compound="DA",
            region=["OF"],
            remove_outliers={"iqr": "calculated"},
        )
        group_stats = ds.group_statistics.iloc[0]
        pd.testing.assert_series_equal(
            group_stats,
            pd.Series(
                {
                    DatasetColumn.GROUP_ID: 1,
                    "compound": "DA",
                    "region": "OF",
                    "shapiro_F": 0.9801003059192006,
                    "shapiro_p": 0.9657026230138452,
                    "is_parametric": True,
                    "mean": 0.015279200000000001,
                    "std": 0.007467412534175111,
                    "sem": 0.0023614031836083297,
                    "values": [
                        0.012725,
                        0.02732,
                        0.019416,
                        0.009179,
                        0.009371,
                        0.021207,
                        0.016852,
                        0.015104,
                        0.00123,
                        0.020388,
                    ],
                },
            ),
            check_names=False,
        )
        quantitative_stats = ds.quantitative_statistics_table.copy()
        cols_without_result = [col for col in quantitative_stats if col != "result"]
        pd.testing.assert_frame_equal(
            quantitative_stats[cols_without_result],
            pd.DataFrame(
                [
                    {
                        "p_value": 0.576,
                        "is_significant": False,
                        "result_string": "F(1, 33) = 0.318, p = 0.576 ",
                        "test": "two_way_anova",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                    {
                        "p_value": 0.17304417123590965,
                        "is_significant": False,
                        "result_string": "F(3, 33) = 1.76, p = 0.17 ",
                        "test": "one_way_anova",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                    {
                        "p_value": [[], []],
                        "is_significant": False,
                        "result_string": pd.DataFrame(columns=["pairs", "p_values"]),
                        "test": "tukey",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                ],
            ),
            check_like=True,
        )
