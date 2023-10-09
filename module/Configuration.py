from module.Filesystem import Filesystem


class Configuration:
    CONFIGURATIONS = [
        "treatment_mapping",
        "experimental_info",
        "compound_ratio_mapping",
        "region_subclassification",
        "compound_subclassification",
    ]

    def __init__(self, experiment):
        self.experiment = experiment
        self.treatment_mapping = self.get_config("treatment_mapping")
        self.experimental_info = self.get_config("experimental_info")
        self.compound_ratio_mapping = self.get_config("compound_ratio_mapping")
        self.region_subclassification = self.get_config("region_subclassification")
        self.compound_subclassification = self.get_config("compound_subclassification")
        
    def get_config(self, configuration):
        for configuration in self.CONFIGURATIONS:
            path = f"{Filesystem.INPUT}/{self.experiment}/{configuration}.json"
            try:
                return Filesystem.getJSON(path)
            except FileNotFoundError:
                Filesystem.saveJSON(path, {})
                print(f"CREATED EMPTY {path} FILE, PLEASE FILL IT IN")
        
    def apply_treatment_mapping(self, df):
        if not self.treatment_mapping:
            raise FileNotFoundError(
                f"MISSING TREATMENT MAPPING FOR {self.experiment}, PLEASE ADD IT TO {Filesystem.INPUT}/{self.experiment}/treatment_mapping.json"
            )
        new_columns = list(list(self.treatment_mapping.values())[0].keys())
        df[new_columns] = df.apply(
            lambda x: self.treatment_mapping[str(int(x["group_id"]))].values(),
            axis=1,
            result_type="expand",
        )  # Get alll the values and assign to corresponding columns
        
        return df.explode("experiments").rename(columns={"experiments": "experiment"})
        
        
    