import os
import json
from module.Questions import Questions


class Filesystem:
    ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
    INPUT = f"{ROOT}/input"
    OUTPUT = f"{ROOT}/output"

    @staticmethod
    def check_directory(path):
        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            print(f"CREATED {path} DIRECTORY")

    @staticmethod
    def check_or_add_raw_data(experiment_name, raw_data_name):
        experiment_directory = f"{Filesystem.INPUT}/{experiment_name}"
        raw_data_filename = f"{experiment_name}_{raw_data_name}.csv"
        while not os.path.exists(f"{experiment_directory}/{raw_data_filename}.csv"):
            add_file = Questions.askYesorNo(
                f"PLEASE ADD {raw_data_filename} TO {experiment_directory}. PRESS 'y' WHEN DONE or 'n' to exit"
            )
            if not add_file:
                return
        return

    @staticmethod
    def saveJSON(path, dict_to_save):
        with open(path, "w", encoding="utf8") as json_file:
            json.dump(dict_to_save, json_file)

    @staticmethod
    def getJSON(path):
        with open(path) as outfile:
            loaded_json = json.load(outfile)
        return loaded_json
