import ipywidgets as widgets
from IPython.display import display, clear_output
from module.core.Metadata import (
    TreatmentInformation,
    ExperimentInformation,
)
from module.core.FileSystem import FileSystem
from module.core.HPLC import HPLC
from module.core.Figure import Correlogram, Histogram, Correlation, Network, Table
from module.core.Constants import COMPOUNDS, COMPOUND_CLASSES, REGION_CLASSES, REGIONS


class OrderedSelectMultipleWithAll(widgets.SelectMultiple):
    def __init__(self, options, **kwargs):
        super().__init__(options=options, **kwargs)

        self.value_ordered = []

        self.observe(self.handle_change, names="value")

    def __setattr__(self, name: str, value: list) -> None:
        if name == "options":
            value = ("Select all", *value)
        super().__setattr__(name, value)

    def __getattribute__(self, name: str) -> tuple:
        return (
            [item for item in super().__getattribute__(name) if item != "Select all"]
            if name == "value"
            else super().__getattribute__(name)
        )

    def handle_change(self, change):
        new_values = set(change.new)
        old_values = set(change.old)

        if new_values == {"Select all"}:
            if "Select all" not in old_values:
                self.value = tuple(self.options)
            else:
                self.value = tuple()
        else:
            added = new_values - old_values
            new_value_ordered = []
            for item in self.value_ordered:
                if item in new_values:
                    new_value_ordered.append(item)
            self.value_ordered = [*new_value_ordered, *added]


figure_mapping = {
    "correlogram": Correlogram,
    "histogram": Histogram,
    "correlation": Correlation,
    "network": Network,
    "table": Table,
}
class_constants_mapping = {
    "compounds": COMPOUND_CLASSES,
    "regions": REGION_CLASSES,
}


def playground():
    projects = FileSystem().list_projects()
    project_dropdown = widgets.Dropdown(
        options=projects, value=projects[0], description="Project:"
    )
    figure_dropdown = widgets.Dropdown(
        options=figure_mapping.keys(),
        value=list(figure_mapping.keys())[0],
        description="Figure:",
    )
    region_ordered_select_multiple = OrderedSelectMultipleWithAll(
        description="Regions:",
        options=REGIONS.list,
        style={"description_width": "initial"},
    )
    region_class_checkbox = OrderedSelectMultipleWithAll(
        description="Classes:",
        options=REGION_CLASSES.list,
        style={"description_width": "initial"},
    )
    compound_ordered_select_multiple = OrderedSelectMultipleWithAll(
        description="Compounds:",
        options=COMPOUNDS.list,
        style={"description_width": "initial"},
    )
    compound_class_checkbox = OrderedSelectMultipleWithAll(
        description="Classes:",
        options=COMPOUND_CLASSES.list,
        style={"description_width": "initial"},
    )
    treatment_ordered_select_multiple = OrderedSelectMultipleWithAll(
        description="Treatments:", options=[], style={"description_width": "initial"}
    )
    experiment_dropdown = widgets.Dropdown(
        description="Experiment:", options=[], style={"description_width": "initial"}
    )
    outliers_dropdown = widgets.Dropdown(
        options=["calculated", "eliminated", False],
        value="calculated",
        description="Eliminate outliers:",
        style={"description_width": "initial"},
    )
    submit = widgets.Button(
        description="Submit", layout=widgets.Layout(width="auto", margin="10px")
    )
    
    def display_interface():
        display(figure_dropdown)
        print(
            "Press control + click to select/deselect multiple items, and shift + click for range select"
        )

        widgets_vbox = widgets.VBox(
            [
                widgets.HBox(
                    [
                        project_dropdown,
                        experiment_dropdown,
                        treatment_ordered_select_multiple,
                        outliers_dropdown,
                    ]
                ),
                widgets.HBox(
                    [
                        region_ordered_select_multiple,
                        region_class_checkbox,
                        compound_ordered_select_multiple,
                        compound_class_checkbox,
                    ]
                ),
                submit,
            ],
            layout=widgets.Layout(align_items="center", justify_content="space-around"),
        )
        display(widgets_vbox)

    display_interface()

    def apply_class_dropdownion(type, dropdown):
        single_select_multiple = {
            "compounds": compound_ordered_select_multiple,
            "regions": region_ordered_select_multiple,
        }[type]

        def apply_class_dropdown_callback(_):
            classes_constant = class_constants_mapping[type]
            single_select_multiple.value = ()
            for selected_class in dropdown.value_ordered:
                for element in classes_constant[selected_class]:
                    single_select_multiple.value = (*single_select_multiple.value_ordered, element)
            
        return apply_class_dropdown_callback

    region_class_checkbox.observe(apply_class_dropdownion("regions", region_class_checkbox), names="value")
    compound_class_checkbox.observe(apply_class_dropdownion("compounds", compound_class_checkbox), names="value")

    def update_project(_):
        project = project_dropdown.value
        print(f"PROJECT: {project}")
        data = HPLC(project).df
        treatment_information = TreatmentInformation(project)
        experiment_information = ExperimentInformation(project)

        compounds_constant = COMPOUNDS.list
        data_compounds = data.compound.unique()
        regions_constant = REGIONS.list
        data_regions = data.region.unique()

        region_ordered_select_multiple.options = [
            region for region in regions_constant if region in data_regions
        ]
        compound_ordered_select_multiple.options = [
            compound for compound in compounds_constant if compound in data_compounds
        ]
        treatment_ordered_select_multiple.options = (
            treatment_information.label.unique()
        )
        experiment_dropdown.options = ["None"] + list(
            experiment_information.label.unique()
        )
        experiment_dropdown.value = "None"

        def experiment_dropdown_updated(experiment_dropdown):
            if experiment_dropdown.new != "None":
                treatment_ordered_select_multiple.value = (
                    experiment_information.df.select(
                        experiment=experiment_dropdown.new
                    ).treatments
                )
            else:
                treatment_ordered_select_multiple.value = ()

        experiment_dropdown.observe(experiment_dropdown_updated, names="value")

        def treatment_dropdown_updated(treatment_ordered_select_multiple):
            if set(treatment_ordered_select_multiple.new) not in [
                set(treatments) for treatments in experiment_information.df.treatments
            ]:
                experiment_dropdown.value = "None"

        treatment_ordered_select_multiple.observe(
            treatment_dropdown_updated, names="value"
        )

    project_dropdown.observe(update_project, names="value")
    
    def display_figure(_):
        clear_output(wait=True)
        display_interface()
        print(f"Generating {figure_dropdown.value}")
        figure_mapping[figure_dropdown.value](
            project=project_dropdown.value,
            compound=compound_ordered_select_multiple.value_ordered,
            region=region_ordered_select_multiple.value,
            experiment=None if experiment_dropdown.value == "None" else experiment_dropdown.value,
            treatment=treatment_ordered_select_multiple.value_ordered,
            remove_outliers=outliers_dropdown.value,
            from_scratch=True,
        )

    submit.on_click(display_figure)

    if len(projects) == 1:
        update_project(None)
