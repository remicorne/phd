import unittest
from module.core.Project import Project, ROOT
import os, shutil

class TestProject(unittest.TestCase):
    def setUp(self):
        # Setup initial conditions here if necessary
        self.project_name = 'TEST'
        self.project_location = f"{ROOT}/{self.project_name}"
        self.raw_data_filename = 'raw_data.xlsx'
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
        project.full_df.select({'is_outlier': True}).select({"treatment": ["vehicles", "MDL"]})
        project.statistics.df
        project.statistics.significant_tests
        project.statistics.significant_results
        project.statistics.insufficent_data
        print('PROJECT INITIALIZATION TEST PASSED')

    def tearDown(self):
        # Cleanup code here
        print("DELETING TEST FILES")
        os.remove(os.path.join(os.getcwd(), self.raw_data_filename))
        shutil.rmtree(self.project_location)
        print('TEST FILES DELETED')
        
        
if __name__ == '__main__':
    unittest.main()