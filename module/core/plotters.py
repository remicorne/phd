import os
from module.core.ProjectDataset import Dataset, MergedDatasets
from module.core.Figure import (
    Histogram,
    SummaryHistogram,
    Correlogram,
    NetworkFigure,
    NetworkDegreesFigure,
)
from module.core.FileSystem import FileSystem
from module.core.Metadata import GroupInformation
from module.core.Matrix import MatrixGroup, NetworkGroup
from module.core.Constants import ConstantRegistry
from module.core.questions import input_escape


def get_dataset(project, request):
    datasets = request["datasets"]
    if len(datasets) == 1:
        dataset, selector = next(iter(request["datasets"].items()))
        dataset = Dataset(project=project, filename=dataset).select(**selector)
    else:
        dataset = MergedDatasets(
            [
                Dataset(project=project, filename=dataset).select(**selector)
                for dataset, selector in request["datasets"].items()
            ]
        )
    return dataset.select(**request.get("selector", {}))


def histogram(project, request, custom_params=None):
    dataset = get_dataset(project, request)
    custom_params = custom_params or {}
    x = hue = custom_params.get("x", "group_name")
    custom_params["palette"] = custom_params.get(
        "palette", dataset.get_palette("color")
    )
    ylabel = ", ".join(dataset.get_units())
    custom_params["ylabel"] = custom_params.get("ylabel", ylabel)
    custom_params["x_order"] = GroupInformation(project).select(
        group_id=dataset.data.group_id.unique()
    )[x]
    title = dataset.get_selection_string()
    location = FileSystem.get_location(
        **{
            "project": project,
            "experiment": dataset.selector.get("experiment", "project"),
        }
    )
    filepath = os.path.join(location, "histogram", title)
    if "experiment" in dataset.selector:
        dataset.calculate_quantitative_statistics()
        statistic = dataset.statistics[0]
    else:
        statistic = []
    Histogram(
        title,
        filepath,
        dataset.data,
        x,
        hue,
        statistic,
        custom_params=custom_params,
    )
    return dataset


def summary_histogram(project, request, custom_params=dict()):
    dataset = get_dataset(project, request)
    custom_params = custom_params or {}
    if isinstance(dataset, MergedDatasets):
        x = "measurement"
    else:
        distinct_measurement_columns = list(
            filter(
                lambda col: len(dataset.data[col].unique()) > 1,
                dataset.measurement_columns,
            )
        )
        if len(distinct_measurement_columns) == 1:
            x = next(iter(distinct_measurement_columns))
        else:
            dataset.to_generic()
            x = "measurement"

    hue = custom_params.get("hue", "group_name")
    custom_params["palette"] = custom_params.get(
        "palette", dataset.get_palette("color")
    )
    custom_params["significance_palette"] = custom_params.get(
        "palette", dataset.get_palette("significance")
    )
    ylabel = ", ".join(dataset.get_units())
    custom_params["ylabel"] = custom_params.get("ylabel", ylabel)
    custom_params["order"] = dataset.data[x].unique()
    custom_params["hue_order"] = dataset.data.group_name.unique()
    title = dataset.get_selection_string()
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "summary_histogram", title)
    if "experiment" in dataset.selector:
        dataset.calculate_quantitative_statistics()
        statistics = dataset.statistics
    else:
        statistics = []
    SummaryHistogram(
        title,
        filepath,
        dataset.data,
        x,
        hue,
        statistics,
        custom_params=custom_params,
    )
    return dataset


def correlogram(project, request, between, custom_params=None):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
    )
    title = dataset.get_selection_string()
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "correlogram", title)
    Correlogram(title, filepath, matrices.matrices, custom_params=custom_params)
    return dataset


def network(project, request, between, custom_params=None):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
    )
    networks = NetworkGroup(matrices).networks
    title = dataset.get_selection_string()

    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "network", title)
    positions = ConstantRegistry.get_registry(name="region_classes_positions").get(
        dataset.selector.get("region")
    )
    NetworkFigure(
        title,
        filepath,
        networks,
        positions=positions,
        custom_params=custom_params,
    )
    return dataset


def network_degrees(project, request, between, custom_params=None):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
    )
    networks = NetworkGroup(matrices).networks
    title = input_escape("ENter figure title/filename")

    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "network_degrees", title)
    NetworkDegreesFigure(
        title,
        filepath,
        networks,
        custom_params=custom_params or {},
    )
    return dataset


def network_summary(project, request, between, measurement: str, custom_params=None):
    custom_params = custom_params or {}
    custom_params["plot_bar"] = False
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data, "group_name", dataset.measurement_columns, between=between
    )
    network_summary_df = (
        NetworkGroup(
            matrices,
        )
        .get_summary_df()
        .select(measurement=measurement)
    )
    x = hue = custom_params.get("x", "group_name")

    location = FileSystem.get_location(  # IMPROVE
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    custom_params["ylabel"] = "AU"
    custom_params["swarm_hue"] = next(iter(between.keys()))
    title = input_escape("ENter figure title/filename")
    filepath = os.path.join(location, "network_summary", title)
    Histogram(
        None,
        filepath,
        network_summary_df,
        x,
        hue,
        custom_params=custom_params,
    )
    return network_summary_df


def summary_network_summary(
    project, request, between, measurement: list[str] = None, custom_params=dict()
):
    custom_params = custom_params or {}
    custom_params["plot_bar"] = False
    custom_params["plot_swarm"] = True
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data, "group_name", dataset.measurement_columns, between=between
    )
    network_summary_df = (
        NetworkGroup(
            matrices,
        )
        .get_summary_df()
        .select(measurement=measurement)
    )
    title = input_escape("ENter figure title/filename")
    hue = custom_params.get("x", "group_name")
    x = "measurement"

    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    custom_params["swarm_hue"] = next(iter(between.keys()))
    filepath = os.path.join(location, title)

    custom_params["significance_palette"] = custom_params.get(
        "palette", dataset.get_palette("significance")
    )
    custom_params["ylabel"] = "AU"
    custom_params["hue_order"] = dataset.data.group_name.unique()
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "summary_network_summary", title)
    SummaryHistogram(
        title,
        filepath,
        network_summary_df,
        x,
        hue,
        custom_params=custom_params,
    )
    return network_summary_df
