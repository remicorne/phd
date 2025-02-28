import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from module.core.ProjectDataset import Dataset, MergedDatasets
from module.core.Figure import (
    Histogram,
    SummaryHistogram,
    Correlogram,
    NetworkFigure,
    NetworkDegreesFigure,
    Correlation,
)
from module.core.FileSystem import FileSystem
from module.core.Metadata import GroupInformation
from module.core.Matrix import MatrixGroup, NetworkGroup
from module.core.Constants import ConstantRegistry
from module.core.questions import input_escape
from module.core.Statistics import QuantitativeStatistic


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
    custom_params["order"] = GroupInformation(project).select(
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


def summary_histogram(project, request, invert_hue=False, custom_params=None):
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
        "significance_palette", dataset.get_palette("significance")
    )
    if "experiment" in dataset.selector:
        dataset.calculate_quantitative_statistics()
        statistics = dataset.statistics
    else:
        statistics = []

    custom_params["order"] = list(dataset.data[x].unique())
    custom_params["hue_order"] = list(dataset.data[hue].unique())

    ylabel = ", ".join(dataset.get_units())
    custom_params["ylabel"] = custom_params.get("ylabel", ylabel)
    title = dataset.get_selection_string()
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "summary_histogram", title)
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
        pvalue_threshold=custom_params.get("p_value_threshold", 0.05),
    )
    title = dataset.get_selection_string()
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "correlogram", title)
    Correlogram(title, filepath, matrices.matrices, custom_params=custom_params)
    return matrices


def network(project, request, between, custom_params=None):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=custom_params.get("p_value_threshold", 0.05),
    )
    networks = NetworkGroup(matrices).networks
    title = dataset.get_selection_string()

    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "network", title)
    region_class = request["datasets"].get("hplc", {}).get("region")
    positions = ConstantRegistry.get_registry(name="region_classes_positions").get(
        region_class
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
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=custom_params.get("p_value_threshold", 0.05),
    )
    network_summary_df = (
        NetworkGroup(
            matrices,
        )
        .get_summary_df()
        .select(measurement=measurement)
    )
    # TODO generalize stats logic + dataset logic for when mouse_id not there
    experiment_information = dataset.experiment_information.iloc[0, :]
    network_summary_df = GroupInformation(project).extend(network_summary_df)
    statistic = (
        QuantitativeStatistic(
            data=network_summary_df,
            group_column="group_name",
            independant_variables=experiment_information.independant_variables,
            is_paired=experiment_information.paired,
            is_parametric=experiment_information.parametric,
            p_value_threshold=0.05,
            delay_execution=False,
            metadata=dict(measurement=measurement),
        )
        if "experiment" in request.get("selector", {})
        else None
    )  # TODO fix by generalizing concept of dataset further DerivedDataset?
    x = hue = custom_params.get("x", "group_name")

    # colormapping by vehicle rank #REMI CLEAN ME
    if "palette" not in custom_params.keys():
        cmap = plt.cm.viridis
        vehicle_values = network_summary_df[
            network_summary_df[hue] == "vehicles"
        ]  # HARD CODE
        compound_value_dict = dict(
            zip(vehicle_values["compound"], vehicle_values["value"])
        )
        sorted_compound_value = {
            k: v
            for k, v in sorted(compound_value_dict.items(), key=lambda item: item[1])
        }
        colors = cmap(np.linspace(0, 1, len(sorted_compound_value)))
        custom_params["palette"] = {
            compound: color
            for compound, color in zip(sorted_compound_value.keys(), colors)
        }

    location = FileSystem.get_location(  # IMPROVE
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    custom_params["ylabel"] = measurement
    custom_params["swarm_hue"] = next(iter(between.keys()))
    title = input_escape("ENter figure title/filename")
    filepath = os.path.join(location, "network_summary", title)

    Histogram(
        None,
        filepath,
        network_summary_df,
        x,
        hue,
        statistic,
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


def correlation(project, x, y, grouper, custom_params=None):
    custom_params = custom_params or {}

    data = []
    for selector in [x, y]:
        dataset = Dataset(project=project, filename=selector.pop("dataset"))
        dataset.select(**selector, **grouper)
        data.append(dataset)

    common_mouse_ids = list(
        set(data[0].data[data[0].project_information.subject_column]).intersection(
            set(data[1].data[data[1].project_information.subject_column])
        )
    )
    x_data = (
        data[0]
        .data.set_index(data[0].project_information.subject_column)
        .loc[common_mouse_ids]
        .value.values
    )
    y_data = (
        data[1]
        .data.set_index(data[1].project_information.subject_column)
        .loc[common_mouse_ids]
        .value.values
    )

    x_label = data[0].get_selection_string() + " " + ", ".join(data[0].get_units())
    y_label = data[1].get_selection_string() + " " + ", ".join(data[1].get_units())

    title = next(iter(grouper.values()))
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "correlation", title)
    Correlation(
        title,
        filepath,
        {"data": x_data, "label": x_label},
        {"data": y_data, "label": y_label},
    )


def statistics_table(project, request):
    if len(request["datasets"]) > 1:
        raise NotImplementedError("Multiple datasets not supported")
    dataset = get_dataset(project, request)
    if "experiment" in dataset.selector:
        dataset.calculate_quantitative_statistics()
        statistics = dataset.statistics
    else:
        statistics = []

    stats_results = []
    for statistic in statistics:
        data = statistic.results
        data = data[data["test"] == statistic.statistical_test][
            ["test", *dataset.measurement_columns, "result_string"]
        ]
        stats_results.append(data)
    stats_results = dataset.sort_values(pd.concat(stats_results))

    index, *column = sorted(
        dataset.measurement_columns,
        key=lambda col: dataset.data[col].unique().size,
        reverse=True,
    )
    stats_table = stats_results.pivot_table(
        index=index,
        columns=["test", *column],
        values="result_string",
        aggfunc="first",
    )

    stats_table = stats_table.sort_index()

    title = dataset.get_selection_string()
    location = FileSystem.get_location(
        **{"project": project, "experiment": dataset.selector.get("experiment", "All")}
    )
    filepath = os.path.join(location, "statistics_table", title + ".xlsx")

    dirpath, _ = os.path.split(filepath)
    os.makedirs(dirpath, exist_ok=True)
    stats_table.to_excel(filepath)
    print(f"Saved {filepath}")
    return stats_table
