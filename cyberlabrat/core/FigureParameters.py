from functools import wraps
from inspect import signature
from pydantic import BaseModel, model_validator
from cyberlabrat.core.questions import yes_or_no
from cyberlabrat.core.FileSystem import FileSystem
from cyberlabrat.core.Metadata import ProjectMetadata


class PlotterValidator(BaseModel, extra="forbid"):
    @classmethod
    def validate(cls, func):
        sig = signature(func)

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            filtered = {
                k: v for k, v in bound.arguments.items() if k in cls.model_fields
            }
            cls.model_validate(filtered)
            return func(*args, **kwargs)

        return wrapper


class MeasurementSelector(PlotterValidator):
    measurement_selector: dict[str, list | str]


class ProjectFigureParameters(PlotterValidator):
    """HAHA"""

    project: str
    dataset_selector: dict[str, MeasurementSelector | None]
    experiment: str | None
    custom_params: dict | None

    @model_validator(mode="after")
    def check_if_init(self):
        if self.project not in FileSystem.list_projects():
            if not yes_or_no(f"Unknown project {self.project}, initialize?"):
                raise ValueError(f"Unknown project {self.project}")
        if unknown_datasets := set(self.dataset_selector.keys()) - set(
            ProjectMetadata(self.project).datasets.label
        ):
            raise ValueError(f"Unknown datasets {unknown_datasets}")
        return self


class HistogramFigureParameters(ProjectFigureParameters):
    """ProjectFigureParameters"""


class SummaryHistogramFigureParameters(HistogramFigureParameters):
    invert_hue: bool


if __name__ == "__main__":
    print(HistogramFigureParameters.__doc__)
