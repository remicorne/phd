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
                    "shapiro_F": 0.8204995939256723,
                    "shapiro_p": 0.01751001431570624,
                    "is_parametric": False,
                    "mean": 0.019057272727272725,
                    "std": 0.014392273684799835,
                    "sem": 0.004339433790235163,
                    "values": [
                        0.01273,
                        0.00123,
                        0.0151,
                        0.01685,
                        0.02121,
                        0.02039,
                        0.00918,
                        0.01942,
                        0.02732,
                        0.05683,
                        0.00937,
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
                        "p_value": 0.5452438477964717,
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
                    "shapiro_F": 0.9800975750452369,
                    "shapiro_p": 0.9656881301230172,
                    "is_parametric": True,
                    "mean": 0.015279999999999998,
                    "std": 0.007467846632954732,
                    "sem": 0.0023615404576956397,
                    "values": [
                        0.01273,
                        0.00123,
                        0.0151,
                        0.01685,
                        0.02121,
                        0.02039,
                        0.00918,
                        0.01942,
                        0.02732,
                        0.00937,
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
                        "p_value": 0.577,
                        "is_significant": False,
                        "result_string": "F(1, 33) = 0.318, p = 0.577 ",
                        "test": "two_way_anova",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                    {
                        "p_value": 0.17306839435164126,
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
                    "shapiro_F": 0.9800975750452369,
                    "shapiro_p": 0.9656881301230172,
                    "is_parametric": True,
                    "mean": 0.015279999999999998,
                    "std": 0.007467846632954732,
                    "sem": 0.0023615404576956397,
                    "values": [
                        0.01273,
                        0.00123,
                        0.0151,
                        0.01685,
                        0.02121,
                        0.02039,
                        0.00918,
                        0.01942,
                        0.02732,
                        0.00937,
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
                        "p_value": 0.577,
                        "is_significant": False,
                        "result_string": "F(1, 33) = 0.318, p = 0.577 ",
                        "test": "two_way_anova",
                        "p_value_threshold": 0.05,
                        "compound": "DA",
                        "region": "OF",
                        "experiment": "agonist_antagonist",
                        "fully_significant": False,
                    },
                    {
                        "p_value": 0.17306839435164126,
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
