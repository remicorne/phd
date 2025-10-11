"""
Utility functions shared across plotting modules.
"""

from typing import Union

import pandas as pd

from module.core.ProjectDataset import ProjectDataset, MergedDatasets
from module.core.Metadata import ProjectMetadata


def get_dataset(project: str, request: dict) -> Union[ProjectDataset, MergedDatasets]:
    """
    Get a dataset or merged datasets based on the request.

    Args:
        project: Name of the project
        request: Dictionary containing dataset specifications and filters

    Returns:
        A ProjectDataset or MergedDatasets instance
    """
    datasets = request["datasets"]
    if len(datasets) == 1:
        dataset, selector = next(iter(request["datasets"].items()))
        dataset = ProjectDataset(project=project, filename=dataset).select(**selector)
    else:
        dataset = MergedDatasets(
            [
                ProjectDataset(project=project, filename=dataset).select(**selector)
                for dataset, selector in request["datasets"].items()
            ]
        )
    if "experiment" in request:
        dataset = dataset.select(experiment=request.get("experiment"))
    return dataset
