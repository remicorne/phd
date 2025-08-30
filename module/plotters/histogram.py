from module.plotters.utils import get_dataset
from module.core.Metadata import ProjectMetadata
from module.core.FileSystem import FileSystem
from module.core.Figure import Histogram, SummaryHistogram
from module.core.ProjectDataset import MergedDatasets


def histogram(project, request, filename=None, custom_params=None):
    dataset = get_dataset(project, request)
    custom_params = custom_params or {}
    x = hue = custom_params.get("x", "group_name")
    custom_params["palette"] = custom_params.get(
        "palette", dataset.get_palette("color")
    )
    ylabel = ", ".join(dataset.get_units())
    custom_params["ylabel"] = custom_params.get("ylabel", ylabel)
    custom_params["order"] = ProjectMetadata(project).groups.select(
        group_id=dataset.data.group_id.unique()
    )[x]
    if "experiment" in request:
        dataset.calculate_quantitative_statistics()
        statistic = dataset.statistics[0]
    else:
        statistic = []
    filename = filename or dataset.get_selection_string()
    filepath = FileSystem.get_location(
        project=project, figure_type="histogram", filename=filename
    )
    Histogram(
        filename,
        filepath,
        dataset.data,
        x,
        hue,
        statistic,
        custom_params=custom_params,
    )
    return dataset


def summary_histogram(
    project, request, filename=None, invert_hue=False, custom_params=None
):
    dataset = get_dataset(project, request)
    custom_params = custom_params or {}
    if isinstance(dataset, MergedDatasets):
        x = "measurement"
    else:
        multiple_measurement_columns = list(
            filter(
                lambda col: len(dataset.data[col].unique()) > 1,
                dataset.measurement_columns,
            )
        )
        if len(multiple_measurement_columns) > 1:
            dataset.to_generic()
            x = "measurement"
        elif len(multiple_measurement_columns) == 1:
            x = next(iter(multiple_measurement_columns))
        else:
            x = dataset.measurement_columns[0]

    hue = custom_params.get("hue", "group_name")
    if invert_hue:
        hue, x = x, hue
        custom_params["invert_hue"] = True
    else:
        custom_params["palette"] = custom_params.get(
            "palette", dataset.get_palette("color")
        )
    custom_params["significance_palette"] = custom_params.get(
        "significance_palette", dataset.get_palette("significance_symbol")
    )
    if "experiment" in request:
        dataset.calculate_quantitative_statistics()
        statistics = dataset.statistics
    else:
        statistics = []

    custom_params["order"] = list(dataset.data[x].cat.categories)
    custom_params["hue_order"] = list(dataset.data[hue].cat.categories)

    ylabel = ", ".join(dataset.get_units())
    custom_params["ylabel"] = custom_params.get("ylabel", ylabel)
    filename = filename or dataset.get_selection_string()
    filepath = FileSystem.get_location(
        project=project, figure_type="summary_histogram", filename=filename
    )
    SummaryHistogram(
        filename,
        filepath,
        dataset.data,
        x,
        hue,
        statistics,
        custom_params=custom_params,
    )
    return dataset
