class Dataset:
    DATASETS_TYPES = ["HPLC", "HT"]
    name = None

    def __init__(self, experiment) -> None:
        self.experiment = experiment
