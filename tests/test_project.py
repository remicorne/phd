import unittest
from module.core.Project import Project, ROOT
import os, shutil

class TestProject(unittest.TestCase):
    def test_project(self):
        
        project_name = 'TEST'
        project_location = f"{ROOT}/{project_name}"
        raw_data_filename = 'raw_data.xlsx'
        
        os.mkdir(project_location)
        
        test_data_dir= f"{os.getcwd()}/tests/test_project_data/"
        
        for filename in os.listdir(test_data_dir):
            test_data_file = os.path.join(test_data_dir, filename)
            project_test_location = os.path.join(project_location, filename)
            if filename == raw_data_filename :
                shutil.copy(test_data_file, os.getcwd())    
            else:
                shutil.copy(test_data_file, project_test_location)
        
        Project(project_name)
        
        print("DELETING TEST FILES")
        
        os.remove(os.path.join(os.getcwd(), raw_data_filename))
        shutil.rmtree(project_location)
        
        print('PROJECT INITIALISATION TEST PASSED')
        
        
        
        
if __name__ == '__main__':
    unittest.main()