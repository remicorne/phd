import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from module.core.ProjectDataset import ProjectDataset, MergedDatasets
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
from module.core.enums import DatasetColumn


def get_dataset(project, request) -> ProjectDataset | MergedDatasets:
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


def correlogram(
    project,
    request,
    between,
    filename=None,
    pvalue_threshold=0.05,
    fdr_correction=None,
    density_thresholding=None,
    custom_params=None,
):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=pvalue_threshold,
        fdr_correction=fdr_correction,
        density_thresholding=density_thresholding,
    )
    filename = filename or dataset.get_selection_string()
    filepath = FileSystem.get_location(
        project=project, figure_type="correlogram", filename=filename
    )
    Correlogram(filename, filepath, matrices.matrices, custom_params=custom_params)
    return matrices


def network(
    project,
    request,
    between,
    filename=None,
    layout=None,
    pvalue_threshold=0.05,
    fdr_correction=None,
    density_thresholding=None,
    custom_params=None,
):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=pvalue_threshold,
        fdr_correction=fdr_correction,
        density_thresholding=density_thresholding,
    )
    networks = NetworkGroup(matrices)
    positions = Registry.get_registry(name="region_classes_positions").get(layout)
    if positions:
        if missing_positions := set(networks.nodes) - set(positions.keys()):
            raise ValueError(f"Missing positions for {missing_positions}")
    filename = filename or dataset.get_selection_string()
    filepath = FileSystem.get_location(
        project=project, figure_type="network", filename=filename
    )
    NetworkFigure(
        filename,
        filepath,
        networks,
        positions=positions,
        custom_params=custom_params,
    )
    return dataset


def network_degrees(
    project,
    request,
    between,
    filename=None,
    pvalue_threshold=0.05,
    fdr_correction=None,
    density_thresholding=None,
    custom_params=None,
):
    custom_params = custom_params or {}
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=pvalue_threshold,
        fdr_correction=fdr_correction,
        density_thresholding=density_thresholding,
    )
    networks = NetworkGroup(matrices).networks
    filename = filename or dataset.get_selection_string()
    filepath = FileSystem.get_location(
        project=project, figure_type="network_degrees", filename=filename
    )
    NetworkDegreesFigure(
        filename,
        filepath,
        networks,
        custom_params=custom_params or {},
    )
    return dataset


def network_summary(
    project,
    request,
    between,
    measurement: str,
    filename=None,
    pvalue_threshold=0.05,
    fdr_correction=None,
    density_thresholding=None,
    custom_params=None,
):
    custom_params = custom_params or {}
    custom_params["plot_bar"] = custom_params.get("plot_bar", False)
    if len(request["datasets"]) > 1:
        raise NotImplementedError("Multiple datasets not yet supported")
    dataset = get_dataset(project, request)
    matrices = MatrixGroup(
        dataset.data,
        "group_name",
        dataset.measurement_columns,
        between=between,
        pvalue_threshold=pvalue_threshold,
        fdr_correction=fdr_correction,
        density_thresholding=density_thresholding,
    )
    network_summary_df = (
        NetworkGroup(
            matrices,
        )
        .get_summary_df()
        .select(measurement=measurement)
    )
    # TODO generalize stats logic + dataset logic for when mouse_id not there
    if "experiment" in request:
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
        project=project, figure_type="network_summary", filename=filename
    )

    Histogram(
        filename,
        filepath,
        network_summary_df,
        x,
        hue,
        statistic,
        custom_params=custom_params,
    )
    return network_summary_df


def summary_network_summary(  # NOT FUNCTIONAL
    project,
    request,
    between,
    filename=None,
    measurement: list[str] = None,
    custom_params=None,
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

    custom_params["significance_palette"] = custom_params.get(
        "palette", dataset.get_palette("significance_symbol")
    )
    custom_params["ylabel"] = "AU"
    custom_params["hue_order"] = dataset.data.group_name.unique()
    filename = filename or dataset.get_selection_string()
    filepath = FileSystem.get_location(
        project=project, figure_type="summary_network_summary", filename=filename
    )
    SummaryHistogram(
        filename,
        filepath,
        network_summary_df,
        x,
        hue,
        custom_params=custom_params,
    )
    return network_summary_df


def correlation(project, x, y, grouper, filename=None, custom_params=None):
    custom_params = custom_params or {}
    data = []
    for selector in [x, y]:
        dataset = ProjectDataset(project=project, filename=selector.pop("dataset"))
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

    filename = next(iter(grouper.values()))
    filepath = FileSystem.get_location(
        project=project, figure_type="correlation", filename=filename
    )
    return Correlation(
        filename,
        filepath,
        {"data": x_data, "label": x_label},
        {"data": y_data, "label": y_label},
    )


def statistics_table(project, request, filename=None):
    if len(request["datasets"]) > 1:
        raise NotImplementedError("Multiple datasets not supported")
    if "experiment" not in request:
        request["experiment"] = "default"
        print(
            "No experiment specified, using default: unpaired, parametric, 1 independant variable"
        )

    dataset = get_dataset(project, request)
    dataset.calculate_quantitative_statistics()
    stats_results = []
    for statistic in dataset.statistics:
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
        project=project,
        figure_type="statistics_table",
        filename=(filename or dataset.get_selection_string()) + ".xlsx",
    )
    stats_table.to_excel(filepath)
    print(f"Saved {filepath}")
    return stats_table
