import shutil
import os
import pickle
from module.Filesystem import Filesystem


class Cache:
    @staticmethod
    def is_cached(self, filepath):
        return os.path.isfile(filepath)

    @staticmethod
    def get_or_add(self, builder_cb):
        if self.isCached():
            return self.getCache()
        # Build useing callback otherwise and cache result
        print(f'BUILDING "{self.data_type}"')
        to_cache = builder_cb()
        self.cache(to_cache)

    @staticmethod
    def cache_pkl(self, to_cache, filepath):
        pkl_filepath = f"{filepath}.pkl"
        with open(pkl_filepath, "wb") as file:
            pickle.dump(to_cache, file)
        print(f"CACHED {pkl_filepath}")

    @staticmethod
    def cache_xlsx(self, df, filepath):
        xlsx_filepath = f"{filepath}.xlsx"
        df.to_excel(xlsx_filepath)
        print(f"CACHED {xlsx_filepath}")

    @staticmethod
    def get_pickle(self):
        print(f'RETRIEVED "{self.data_type}" FROM CACHE')
        with open(self.cache_filepath, "rb") as file:
            return pickle.load(file)

    @staticmethod
    def resetCache(self):
        shutil.rmtree(self.path)
        os.mkdir(self.path)
        print("CACHE CLEARED")
