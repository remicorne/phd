import os, shutil
import unittest
from module.core.HPLC import HPLC
from module.core.Statistics import QuantitativeStatistic, AggregateStatistics
from module.core.Project import Project, ROOT
from module.core.Figure import (
    Histogram,
    Correlation,
    Correlogram,
    Network,
    Table,
    StatisticsTable,
)


class TestProject(unittest.TestCase):
    def setUp(self):
        # Setup initial conditions here if necessary

        self.project_name = "TEST"
        self.project_location = f"{ROOT}/{self.project_name}"
        self.raw_data_filename = "raw_data.xlsx"

        # If deletion idn't happen (quit debug session too early for ex)
        raw_data_location = os.path.join(os.getcwd(), self.raw_data_filename)
        if os.path.exists(raw_data_location):
            os.remove(raw_data_location)
        if os.path.exists(self.project_location):
            shutil.rmtree(self.project_location)

        os.mkdir(self.project_location)

        # Copy necessary test data to the project location
        test_data_dir = f"{os.getcwd()}/tests/test_project_data/"
        for filename in os.listdir(test_data_dir):
            test_data_file = os.path.join(test_data_dir, filename)
            project_test_location = os.path.join(self.project_location, filename)
            if filename == self.raw_data_filename:
                shutil.copy(test_data_file, os.getcwd())
            else:
                shutil.copy(test_data_file, project_test_location)

    def test_project(self):
        project = Project(self.project_name)
        assert project.data.columns.to_list() == [
            "mouse_id",
            "group_id",
            "value",
            "compound",
            "region",
            "treatment",
            "color",
            "independant_variables",
            "label",
            "is_outlier",
            "outlier_status",
        ]
        assert len(project.data) == 1152
        assert (
            len(
                project.data.select(is_outlier=True).select(
                    treatment=["vehicles", "MDL"]
                )
            )
            == 28
        )
        assert len(project.statistics.significant_results) == 12

    def test_stats(self):
        stats = QuantitativeStatistic.calculate(
            project=self.project_name,
            experiment="agonist antagonist",
            compound="DA/5HT",
            region="OF",
            p_value_threshold=0.05,
        )
        assert stats.is_significant.to_list() == [False, True, True]
        stats = QuantitativeStatistic.calculate(
            project=self.project_name,
            experiment="agonist antagonist",
            compound="DA/5HT",
            region="OF",
            p_value_threshold=0.01,
        )
        assert stats.is_significant.to_list() == [False, True, False]
        assert len(stats.select(fully_significant=True)) == 0

        agg_stats = AggregateStatistics(project=self.project_name)
        assert (
            len(
                agg_stats.select(
                    is_parametric=True, compound=lambda compound: "/" not in compound
                )
            )
            == 19
        )
        assert len(agg_stats.select(is_parametric=False)) == 24

    def test_figures(self):
        assert Histogram(
            project=self.project_name,
            experiment="agonist antagonist",
            compound="5HT/5HTP",
            region="OF",
            remove_outliers="calculated",
            from_scratch=True,
        ).statistics_table.is_significant.to_list() == [False, False, False]
        Correlogram(
            project=self.project_name,
            compound="5HT-DA",
            from_scratch=True,
            remove_outliers="calculated",
        )
        Correlation(
            project=self.project_name,
            experiment="agonist antagonist",
            treatment="vehicles",
            compound="DA",
            region=["OF", "CB"],
            from_scratch=True,
            remove_outliers="calculated",
        )
        Network(
            project=self.project_name,
            compound="DA, 5HT",
            region="thalamocortical_interaction",
            from_scratch=True,
            remove_outliers="calculated",
        )
        Table(
            project=self.project_name,
            treatment="vehicles",
            from_scratch=True,
            remove_outliers="calculated",
        )
        Histogram(
            project=self.project_name,
            compound="weight",
            region="OF",
            from_scratch=True,
            remove_outliers="calculated",
        )
        StatisticsTable(
            project=self.project_name,
            compound="neurotransmitters",
            experiment="agonist antagonist",
            from_scratch=True,
        )
        Histogram(
            project=self.project_name,
            experiment="agonist antagonist",
            from_scratch=True,
            remove_outliers="calculated",
            pool="treatment",
        )

    def test_utils(self):
        filtered = (
            HPLC(project=self.project_name)
            .select(region="thalamocortical_interaction", compound="monoamines")
            .select(region="cortex", compound="neurotransmitters")
        )
        assert set(filtered.region.unique()) == {"OF"}
        assert set(filtered.compound.unique()) == {"DA", "5HT"}

    def tearDown(self):
        # Cleanup code here
        print("DELETING TEST FILES")
        os.remove(os.path.join(os.getcwd(), self.raw_data_filename))
        shutil.rmtree(self.project_location)
        print("TEST FILES DELETED")


if __name__ == "__main__":
    unittest.main()
