from cyberlabrat.core.Statistics import QuantitativeStatistic, AggregateStatistics
import unittest
from cyberlabrat.core.Project import Project, ROOT
from cyberlabrat.core.Figure import Histogram, Correlation, Correlogram, Network, Table, StatisticsTable
import os, shutil


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
        # This runs the test
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
            == 18
        )
        stats = QuantitativeStatistic.calculate(
            project=self.project_name,
            experiment="agonist antagonist",
            compound="DA/5HT",
            region="OF",
            p_value_threshold=0.05,
        )
        assert stats.is_significant.to_list() == [True, True, True]
        stats = QuantitativeStatistic.calculate(
            project=self.project_name,
            experiment="agonist antagonist",
            compound="DA/5HT",
            region="OF",
            p_value_threshold=0.01,
        )
        assert stats.is_significant.to_list() == [False, True, False]
        stats = QuantitativeStatistic.calculate(
            project="TCB2", experiment="agonist antagonist", compound="5HT"
        )
        AggregateStatistics(project=self.project_name)
        assert len(AggregateStatistics(project=self.project_name).select(is_parametric=True, compound=lambda compound: '/' not in compound)) == 20
        assert len(AggregateStatistics(project=self.project_name).select(is_parametric=False)) == 4
        assert len(stats.select(fully_significant=True)) == 9
        assert len(project.statistics.significant_results) == 6

        assert Histogram(
            project="TCB2",
            experiment="dose response",
            compound="5HIAA/5HT",
            region="OF",
            from_scratch=True,
        ).statistics_table.is_significant.to_list() == [True, True]
        assert Histogram(
            project="TCB2",
            experiment="dose response",
            compound="5HIAA/5HT",
            from_scratch=True,
            remove_outliers="calculated",
        ).statistics_table.is_significant.all()
        legit_regions = [
            "OF",
            "PL",
            "aCC",
            "M1",
            "SJ",
            "S1L",
            "S1R",
            "AC",
            "V1",
            "A",
            "dH",
            "vH",
            "NAc",
            "VM",
            "DM",
            "VL",
            "DL",
            "MD",
            "VPL",
            "VPR",
            "DLG",
            "HY",
            "SC",
            "SN",
            "VTA",
            "DR",
            "MR",
            "CB",
        ]
        Correlogram(
            project="TCB2",
            compound="5HT-DA",
            region=legit_regions,
            from_scratch=True,
            remove_outliers="calculated",
        )
        Correlation(
            project="TCB2",
            experiment="agonist antagonist",
            treatment="vehicles",
            compound="DA",
            region=["VPL", "SN"],
            from_scratch=True,
            remove_outliers="calculated",
        )
        node_positions = {
            "OF": (5.9, 2.4),
            "PL": (5.6, 3.5),
            "aCC": (7.1, 4),
            "M1": (6.2, 3.8),
            "SJ": (7.5, 3.8),
            "S1L": (8.5, 3.8),
            "S1R": (8.2, 4.1),
            "AC": (9.9, 4.2),
            "V1": (10.9, 4),
            "MD": (10.3, 1.45),
            "VPR": (9.1, 1.45),
            "VPL": (9.4, 1.15),
            "DLG": (9.7, 1.9),
        }

        Network(
            project="TCB2",
            node_positions=node_positions,
            compound="5HT",
            region="thalamocortical_interaction",
            from_scratch=True,
            remove_outliers="calculated",
        )
        Table(
            project="TCB2",
            compound=["5HT", "5HIAA", "DA", "DOPAC", "HVA", "3MT", "NA"],
            treatment="vehicles",
            region=legit_regions,
            from_scratch=True,
            remove_outliers="calculated",
        )
        Histogram(project='TEST', experiment="weight", region='OF', from_scratch=True, remove_outliers="calculated")
        StatisticsTable(project="TCB2", compound=["DA", "5HT"], experiment ="agonist antagonist", from_scratch=True)
        Histogram(project='TCB2', 
                experiment='agonist antagonist',
                compound=['DA', 'NA'], 
                region=['OF', 'PL'],
                from_scratch=True, 
                remove_outliers='calculated',
                pool="treatment"
                )
        
    def tearDown(self):
        # Cleanup code here
        print("DELETING TEST FILES")
        os.remove(os.path.join(os.getcwd(), self.raw_data_filename))
        shutil.rmtree(self.project_location)
        print("TEST FILES DELETED")


if __name__ == "__main__":
    unittest.main()
