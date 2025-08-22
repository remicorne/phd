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
from module.core.Metadata import ProjectMetadata
from module.core.Matrix import MatrixGroup, NetworkGroup
from module.core.Registry import Registry
from module.core.Statistics import QuantitativeStatistic
from module.core.constants import DatasetColumn


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


def histogram(project, request, title=None, filename=None, custom_params=None):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
    filepath = FileSystem.get_location(
        project=project, figure_type="histogram", filename=filename or title
    )
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
    if "experiment" in dataset.selector:
        dataset.calculate_quantitative_statistics()
        statistic = dataset.statistics[0]
    else:
        statistic = []
    Histogram(
        title or filename,
        filepath,
        dataset.data,
        x,
        hue,
        statistic,
        custom_params=custom_params,
    )
    return dataset


def summary_histogram(
    project, request, title=None, filename=None, invert_hue=False, custom_params=None
):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
    filepath = FileSystem.get_location(
        project=project, figure_type="summary_histogram", filename=filename or title
    )
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
    if "experiment" in dataset.selector:
        dataset.calculate_quantitative_statistics()
        statistics = dataset.statistics
    else:
        statistics = []

    custom_params["order"] = list(dataset.data[x].cat.categories)
    custom_params["hue_order"] = list(dataset.data[hue].cat.categories)

    ylabel = ", ".join(dataset.get_units())
    custom_params["ylabel"] = custom_params.get("ylabel", ylabel)
    SummaryHistogram(
        title or filename,
        filepath,
        dataset.data,
        x,
        hue,
        statistics,
        custom_params=custom_params,
    )
    return dataset


def correlogram(project, request, between, title=None, filename=None, custom_params=None):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=custom_params.get("p_value_threshold", 0.05),
        fdr_correction=custom_params.get("fdr_correction", False),
    )
    filepath = FileSystem.get_location(
        project=project, figure_type="correlogram", filename=filename or title
    )
    Correlogram(title or filename, filepath, matrices.matrices, custom_params=custom_params)
    return matrices


def network(project, request, between, title=None, filename=None, layout=None, custom_params=None):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=custom_params.get("p_value_threshold", 0.05),
    )
    networks = NetworkGroup(matrices)
    filepath = FileSystem.get_location(
        project=project, figure_type="network", filename=filename or title
    )
    positions = Registry.get_registry(name="region_classes_positions").get(layout)
    if positions:
        if missing_positions := set(networks.nodes) - set(positions.keys()):
            raise ValueError(f"Missing positions for {missing_positions}")
    NetworkFigure(
        title or filename,
        filepath,
        networks,
        positions=positions,
        custom_params=custom_params,
    )
    return dataset


def network_degrees(project, request, between, title=None, filename=None, custom_params=None):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
    )
    networks = NetworkGroup(matrices).networks
    filepath = FileSystem.get_location(
        project=project, figure_type="network_degrees", filename=filename or title
    )
    NetworkDegreesFigure(
        title or filename,
        filepath,
        networks,
        custom_params=custom_params or {},
    )
    return dataset


def network_summary(project, request, between, measurement: str, title=None, filename=None, custom_params=None):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
    custom_params = custom_params or {}
    custom_params["plot_bar"] = custom_params.get("plot_bar", False)
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
    if "experiment" in request.get("selector", {}):
        statistic = QuantitativeStatistic(
            data=network_summary_df,
            group_column="group_name",
            independant_variables=dataset.selected_experiment.independant_variables,
            is_paired=dataset.selected_experiment.paired,
            is_parametric=dataset.selected_experiment.parametric,
            p_value_threshold=0.05,
            delay_execution=False,
            metadata=dict(measurement=measurement),
        )  # TODO fix by generalizing concept of dataset further DerivedDataset?
    else:
        statistic = None
    network_summary_df = ProjectMetadata(project).groups.extend(network_summary_df)
    x = hue = custom_params.get("x", "group_name")

    # colormapping by vehicle rank #REMI CLEAN ME
    if "palette" not in custom_params.keys():
        cmap = plt.cm.viridis
        vehicle_values = network_summary_df[
            network_summary_df[hue] == "vehicles"
        ]  # HARD CODE
        compound_value_dict = dict(
            zip(vehicle_values["compound"], vehicle_values[DatasetColumn.VALUE])
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

    custom_params["ylabel"] = measurement
    custom_params["swarm_hue"] = next(iter(between.keys()))
    filepath = FileSystem.get_location(
        project=project, figure_type="network_summary", filename=filename or title
    )

    Histogram(
        title or filename,
        filepath,
        network_summary_df,
        x,
        hue,
        statistic,
        custom_params=custom_params,
    )
    return network_summary_df


def summary_network_summary(
    project, request, between, title=None, filename=None, measurement: list[str] = None, custom_params=dict()
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
    hue = custom_params.get("x", "group_name")
    x = "measurement"

    custom_params["swarm_hue"] = next(iter(between.keys()))
    filepath = FileSystem.get_location(
        project=project, figure_type="summary_network_summary", filename=filename or title
    )

    custom_params["significance_palette"] = custom_params.get(
        "palette", dataset.get_palette("significance_symbol")
    )
    custom_params["ylabel"] = "AU"
    custom_params["hue_order"] = dataset.data.group_name.unique()
    SummaryHistogram(
        title or filename,
        filepath,
        network_summary_df,
        x,
        hue,
        custom_params=custom_params,
    )
    return network_summary_df


def correlation(project, x, y, grouper, title=None, filename=None, custom_params=None):
    custom_params = custom_params or {}
    data = []
    for selector in [x, y]:
        dataset = Dataset(project=project, filename=selector.pop("dataset"))
        dataset.select(**selector, **grouper)
        data.append(dataset)
    subject_column = DatasetColumn.SUBJECT_ID
    common_subject_ids = list(
        set(data[0].data[subject_column]).intersection(
            set(data[1].data[subject_column])
        )
    )
    x_data = data[0].data.set_index(subject_column).loc[common_subject_ids].value.values
    y_data = data[1].data.set_index(subject_column).loc[common_subject_ids].value.values

    x_label = data[0].get_selection_string() + " " + ", ".join(data[0].get_units())
    y_label = data[1].get_selection_string() + " " + ", ".join(data[1].get_units())

    Correlation(
        title or filename,
        FileSystem.get_location(
            project=project, figure_type="correlation", filename=filename or title
        ),
        {"data": x_data, "label": x_label},
        {"data": y_data, "label": y_label},
    )


def statistics_table(project, request, title=None, filename=None):
    if not (title or filename):
        raise Exception("Must specify parameter filename or title")
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
    stats_results = dataset.sort_values(pd.concat(stats_results), dataset.selection)

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

    filepath = FileSystem.get_location(
        project=project, figure_type="statistics_table", filename=(filename or title) + ".xlsx"
    )
    stats_table.to_excel(filepath)
    print(f"Saved {filepath}")
    return stats_table
