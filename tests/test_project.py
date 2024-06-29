from module.core.Statistics import Statistics
import unittest
from module.core.Project import Project, ROOT
from module.core.Figure import Histogram
import os, shutil

class TestProject(unittest.TestCase):
    def setUp(self):
        # Setup initial conditions here if necessary
        self.project_name = 'TEST'
        self.project_location = f"{ROOT}/{self.project_name}"
        self.raw_data_filename = 'raw_data.xlsx'
        
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
        assert project.data.columns.to_list() == ['mouse_id', 'group_id', 'value', 'compound', 'region', 'treatment', 'color', 'independant_variables', 'label', 'is_outlier', 'outlier_status']
        assert len(project.data) == 1152
        assert len(project.data.select(is_outlier=True).select(treatment=["vehicles", "MDL"])) == 18
        stats = Statistics(self.project_name).get_quantitative_stats(experiment='agonist antagonist', compound='DA/5HT', region='OF', p_value_threshold=0.05)
        assert stats.results.is_significant.to_list() == [True, True, True]
        stats = Statistics(self.project_name).get_quantitative_stats(experiment='agonist antagonist', compound='DA/5HT', region='OF', p_value_threshold=0.01)
        assert stats.results.is_significant.to_list() == [False, True, False]
        assert len(project.statistics.select(is_significant = True)) == 9
        assert len(project.statistics.significant_results) == 6
        assert len(project.statistics.select(fully_significant=True)) == 6
        assert len(project.statistics.insufficent_data) == 5
        assert Histogram("TCB2", experiment='dose response', compound="5HIAA/5HT", region="OF", from_scratch=True, handle_outliers=False).statistics.is_significant.to_list() == [True, True]
        assert Histogram('TCB2', experiment='dose response', compound="5HIAA/5HT", from_scratch=True, handle_outliers=False).statistics.is_significant.all()
        print('PROJECT INITIALIZATION TEST PASSED')

    def tearDown(self):
        # Cleanup code here
        print("DELETING TEST FILES")
        os.remove(os.path.join(os.getcwd(), self.raw_data_filename))
        shutil.rmtree(self.project_location)
        print('TEST FILES DELETED')

        
if __name__ == '__main__':
    unittest.main()