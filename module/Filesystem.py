import os
import json
from module.Configuration import Configuration
from module.Dataset import Dataset
from module.Questions import Questions


class Filesystem:
    ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
    INPUT = f"{ROOT}/input"
    OUTPUT = f"{ROOT}/output"

    @staticmethod
    def check_directory(filepath):
        if not os.path.exists(filepath):
            os.mkdir(filepath)
            print(f"CREATED {filepath} DIRECTORY")

    @staticmethod
    def check_input_directory():
        Filesystem.check_directory(Filesystem.INPUT)

    @staticmethod
    def check_output_directory():
        Filesystem.check_directory(Filesystem.OUTPUT)

    @staticmethod
    def check_config_file(experiment, config_file):
        path = f"{Filesystem.INPUT}/{experiment}/{config_file}.json"
        if not os.path.isfile(path):
            with open(path, "w", encoding="utf8") as json_file:
                json.dump({}, json_file)
            print(f"CREATED EMPTY {path} FILE, PLEASE FILL IT IN")

    @staticmethod
    def check_or_add_raw_data(experiment_name, raw_data_name):
        experiment_directory = f"{Filesystem.INPUT}/{experiment_name}"
        raw_data_filename = f"{experiment_name}_{raw_data_name}.csv"
        while not os.path.exists(f"{experiment_directory}/{raw_data_name}.csv"):
            if not (
                Questions.askYesorNo(
                    f"PLEASE ADD {raw_data_filename} TO {experiment_directory}. PRESS 'y' WHEN DONE or 'n' TO CONTINUE WITHOUT IT"
                )
            ):
                print(f"MISSING {raw_data_name} RAW DATA. CONTINUING WITHOUT IT")
                return
        print(f"FOUND {raw_data_name} RAW DATA")

    @staticmethod
    def initiate_experiment(experiment_name):
        Filesystem.check_directory(f"{Filesystem.INPUT}/{experiment_name}")
        Filesystem.check_directory(f"{Filesystem.OUTPUT}/{experiment_name}")
        for configuration in Configuration.CONFIGURATIONS:
            Filesystem.check_config_file(experiment_name, configuration)
        for dataset_type in Dataset.DATASETS_TYPES:
            Filesystem.check_or_add_raw_data(experiment_name, dataset_type)
