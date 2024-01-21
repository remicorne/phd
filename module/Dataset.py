import os
import pickle
from dataclasses import dataclass

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
INPUT_DIR = f"{ROOT}/input"
OUTPUT_DIR = f"{ROOT}/output"
CACHE_DIR = f"{INPUT_DIR}/cache"


class Dataset:
    
    type: str = None
    project: str
    
    def __post_init__(self):
        self.pkl_path = f"{os.getcwd()}/input/{self.project}/{self.type}.pkl"
        
    
    def load(self):
        with open(self.pkl_path, "rb") as file:
            dataset = pickle.load(file)
            print(f'RETRIEVED "{self.type}" FROM "{self.project}" CACHE')
            return dataset

    def save(self, dataset):
        with open(self.pkl_path, "wb") as file:
            pickle.dump(dataset, file)
        print(f"CACHED {self.project}/{self.type}.pkl")
        
    
    

    @property
    def is_cached(self):
        return os.path.isfile(self.self.filepath)
