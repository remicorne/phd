from module.core.FileSystem import FileSystem
from module.core.HPLC import HPLC, TissueWeight
from module.core.Parameters import StandardizedParameter

projects = FileSystem.list_projects()

for project in projects:

    hplc = HPLC(project=project)
    tissue_weights = TissueWeight(project=project)

    for dataset in [hplc, tissue_weights]:
        df = dataset.df
        for parameter_type in ["compound", "region"]:
            df[parameter_type] = df[parameter_type].apply(
                lambda x: str(x.name) if isinstance(x, StandardizedParameter) else x
            )
        dataset.save(df)
