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
        self.__file = None
        self.__data = None
