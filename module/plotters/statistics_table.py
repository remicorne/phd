import pandas as pd
from module.plotters.utils import get_dataset
from module.core.FileSystem import FileSystem


def statistics_table(project, request, filename=None):
    if len(request["datasets"]) > 1:
        raise NotImplementedError("Multiple datasets not supported")
    if "experiment" not in request:
        request["experiment"] = "default"
        print(
            "No experiment specified, using default: unpaired, parametric, 1 independant variable"
        )

    dataset = get_dataset(project, request)
    stats_results = []
    for statistic in dataset.quantitative_statistics:
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
