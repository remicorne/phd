from module.core.ProjectDataset import ProjectDataset
from module.core.ProjectDataset import DatasetColumn
from module.core.Figure import Correlation
from module.core.FileSystem import FileSystem


def correlation(project, x, y, filename=None, custom_params=None):
    custom_params = custom_params or {}
    data: list[ProjectDataset] = []
    for selector in [x, y]:
        dataset = ProjectDataset(project=project, filename=selector.pop("dataset"))
        dataset.select(**selector)
        data.append(dataset)
    subject_column = DatasetColumn.SUBJECT_ID
    common_subject_ids = list(
        set(data[0].data[subject_column]).intersection(
            set(data[1].data[subject_column])
        )
    )
    if not common_subject_ids:
        raise ValueError(
            "No subject ids in common for current selection, nothing to correlate"
        )
    x_data = data[0].data.set_index(subject_column).loc[common_subject_ids].value.values
    y_data = data[1].data.set_index(subject_column).loc[common_subject_ids].value.values

    x_label = data[0].get_selection_string() + " " + ", ".join(data[0].get_units())
    y_label = data[1].get_selection_string() + " " + ", ".join(data[1].get_units())

    filepath = FileSystem.get_location(
        project=project,
        figure_type="correlation",
        filename=filename or x_label + y_label,
    )
    Correlation(
        filename,
        filepath,
        {"data": x_data, "label": x_label},
        {"data": y_data, "label": y_label},
    )
