import numpy as np
import matplotlib.pyplot as plt

from module.plotters.utils import get_dataset
from module.core.Matrix import MatrixGroup, NetworkGroup, NetworkCharacteristic
from module.core.Figure import (
    Histogram,
    SummaryHistogram,
    Correlogram,
    NetworkFigure,
    NetworkDegreesFigure,
)
from module.core.Registry import Registry
from module.core.FileSystem import FileSystem
from module.core.Statistics import QuantitativeStatistic
from module.core.ProjectDataset import DatasetColumn
from module.core.Metadata import ProjectMetadata


def correlogram(
    project,
    request,
    between,
    filename=None,
    pvalue_threshold=0.05,
    fdr_correction=False,
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


def network(
    project,
    request,
    between,
    filename=None,
    layout=None,
    pvalue_threshold=0.05,
    fdr_correction=False,
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
    return networks.get_summary_df()


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


def network_summary(
    project,
    request,
    between,
    measurement: NetworkCharacteristic,
    filename=None,
    pvalue_threshold=0.05,
    fdr_correction=False,
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
    # Reintroduce group infos
    network_summary_df = ProjectMetadata(project).groups.extend(network_summary_df)
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
    return statistic.results


# TODO
# def summary_network_summary(  # NOT FUNCTIONAL
#     project,
#     request,
#     between,
#     filename=None,
#     measurement: list[str] = None,
#     custom_params=None,
# ):
#     custom_params = custom_params or {}
#     custom_params["plot_bar"] = False
#     custom_params["plot_swarm"] = True
#     dataset = get_dataset(project, request)
#     matrices = MatrixGroup(
#         dataset.data, "group_name", dataset.measurement_columns, between=between
#     )
#     network_summary_df = (
#         NetworkGroup(
#             matrices,
#         )
#         .get_summary_df()
#         .select(measurement=measurement)
#     )
#     hue = custom_params.get("x", "group_name")
#     x = "measurement"

#     custom_params["swarm_hue"] = next(iter(between.keys()))

#     custom_params["significance_palette"] = custom_params.get(
#         "palette", dataset.get_palette("significance_symbol")
#     )
#     custom_params["ylabel"] = "AU"
#     custom_params["hue_order"] = dataset.data.group_name.unique()
#     filename = filename or dataset.get_selection_string()
#     filepath = FileSystem.get_location(
#         project=project, figure_type="summary_network_summary", filename=filename
#     )
#     SummaryHistogram(
#         filename,
#         filepath,
#         network_summary_df,
#         x,
#         hue,
#         custom_params=custom_params,
#     )
#     return network_summary_df
