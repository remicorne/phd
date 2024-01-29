
from dataclasses import dataclass
from module.core.Questions import Questions
from module.core.Metadata import Metadata
from module.core.HPLC import HPLC
from module.constants import ROOT
import os


@dataclass
class Project:
        
    name: str
    
    def __post_init__(self):
        self.location = f"{ROOT}/{self.name}"
        if not os.path.exists(self.location):
            if Questions.yes_or_no(f"INITIALIZE NEW PROJECT: '{self.name}' ?"):
                os.mkdir(self.location)
            else:
                print(f"UNKNOWN PROJECT: {self.name}")
                print(f"KNOW PROJECTS ARE: {os.listdir(ROOT)}")
                exit(1)
        self.hplc = HPLC(self.name)
        # self.metadata = Metadata(self.name)
            

        

            
        
        