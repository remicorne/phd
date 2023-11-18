import shutil
import os
import pickle
from module.Filesystem import Filesystem

class Cache:

    def __init__(self, path, data_type) -> None:
        self.path = path
        self.data_type = data_type
        self.cache_filepath = f"{self.path}/{self.data_type}"

    
    def isCached(self):
        return os.path.isfile(f"{self.cache_filepath}.pkl")

    
    
    def get_or_add(self, builder_cb):
        if self.isCached():
            return self.getCache()
        # Build useing callback otherwise and cache result
        print(f'BUILDING "{self.data_type}"')
        to_cache = builder_cb()
        self.cache(to_cache)


    def cache_pkl(self, to_cache):
        with open(self.cache_filepath, "wb") as file:
            pickle.dump(to_cache, file)
        print(f"CACHED {self.cache_filepath}.pkl")

    def cache_xlsx(self, df):
        df.to_excel(f"{self.cache_filepath}.xlsx")
        print(f"CACHED {self.cache_filepath}.xlsx")


    def getCache(self):
        print(f'RETRIEVED "{self.data_type}" FROM CACHE')
        with open(self.cache_filepath, "rb") as file:
            return pickle.load(file)
    
    
    def resetCache(self):
        shutil.rmtree(self.path)
        os.mkdir(self.path)
        print("CACHE CLEARED")