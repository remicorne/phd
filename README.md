# Data Analysis Tool

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)

A Python-based toolkit for analyzing data, calculating statistics and generating figures designed to run in Jupyter notebooks.

The documentation covers only basic functionnalites currently, contact remi.corne@gmail.com or 01111996jasmine@gmail.com for suggestions or questions. I also recommend using free AI tools (in-editor agents) to explain the code.

I don't recommend extending the code as I've implemented some things in bordeline criminal ways to get them done fast. Future versions will try to avoid exotic design patterns.

## Getting Started

### Prerequisites
- Python 3.11
- Jupyter Notebook/Lab

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/remicorne/phd.git
   ```

2. Create and activate a virtual environment:
   - **Windows**:
     ```powershell
     python -m venv venv
     .\venv\Scripts\activate
     ```
   - **macOS/Linux**:
     ```bash
     python3 -m venv venv
     source venv/bin/activate
     ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Call a plotter function with the required parameters and the code will walk you through whatever might need to be done (initializing a project, adding a dataset, editing metadata etc.)

## Usage

### Example use cases
A sample of already published neurochemical data comes packaged with the code (module/example_project/tcb2_hplc_data.csv)

A jupyter notebook 'INTERFACE.ipynb' is provided to guide you through the process of initializing a project, adding a dataset, editing metadata etc.

Run a cell (any plotter) in INTERFACE.ipynb and follow the instructions:
1. User is prompted on whether to initialize a project if the project name is not recognized (y/n)
2. Project is initialized (folder created and metadata.xlsx created)
3. The metadata.xlsx file is opened in the default editor. User is requested to edit metadata.xlsx (necessary before adding datasets, see **Metadata** section). **Edit unnecesary when running the tutorial as the metadata file is already consistent with the example data.**
4. The user is prompted whether to add dataset (y/n) (see **Dataset** section)
5. User must input raw dataset filepath (use 'module/example_project/tcb2_hplc_data.csv' for example)
6. The dataset is loaded, validated and a copy is save in the project's folder
7. The figure is generated and saved in PROJECTS/{project_name}/figure_type/figure_name.{svg|png|xlsx}

### Jupyter Notebook
The primary way to use this code is through Jupyter notebooks (see INTERFACE.ipynb for examples)

### Standalone Python Script
Work as well for running as a standalone Python script but the multiprocessing may pose problems on window due to how child processes are spawned. You may have to use freeze support in this case:

```python
if __name__ == "__main__":
    import multiprocessing as mp
    mp.freeze_support()

    # Your plotting or analysis code here
    from module.core import your_analysis_function
    your_analysis_function(*parameters)
```

## Repository Structure
```
 phd/
 ├── module/                      # Core Python modules
 │   └── core/                    # Main analysis code
 │   └── json/                    # JSON files for project metadata
*├── PROJECTS/                    # Project data and results (automatically created when initialising a project or creating figures)
*│   └── some_project/            # User created project
*│       └── dataset.pkl          # A dataset of the project (added by user)
*│       └── metadata.xlsx        # Metadata of the project (auto generated + user edits)
*│       └── figure_type/         # Figures of of a given type (auto generated)
 ├── tests/                       # Test suite (run with pytest)
 └── requirements.txt             # Project dependencies
```
(\* gitignored)

## Project Management

The `PROJECTS/` directory is the central location for storing project data and generated figures. Here's how it works:

### Project Initialization
When you create a plot or analysis:
1. The system checks for an existing project with the specified name
2. If the project doesn't exist, you'll be prompted to initialize it
3. A new project folder is created with a `metadata.xlsx` template (more on metadata later)

### Project Structure
Each project follows this structure:
- `{dataset_name}.pkl`: Your project's dataset files (manually added)
- `metadata.xlsx`: Project metadata (auto-generated, requires user input)
- `{figure_type}/`: Directories for different figure types (auto-created)
  - `{figure_name}.{svg|png|xlsx}`: Generated figures with automatic naming

### Workflow
1. **Data Loading**:
   - If the requested dataset isn't found, you'll be prompted to provide the file path
   - The dataset is then loaded (more on datasets later)
   - The requested data is processed for figure generation/stats etc.. (more on requests later)

2. **Figure Generation**:
   - Figures are automatically saved in `PROJECTS/{project_name}/{figure_type}/`
   - Figure names and titles are auto-generated for consistency
   - This allows for easy exploration and comparison of different visualizations

## Datasets

Datasets represent collections of measurements collected using the same technique and sharing common identifying characteristics. They are designed with flexibility in mind while maintaining a consistent structure for analysis.

### Structure
All datasets must include these required columns:
- `subject_id`: Unique identifier for the subject being measured
- `value`: The measurement value
- `unit`: The measurement unit
- One or more measurement characteristic columns (e.g., 'behavior', 'compound', 'location')

### Examples

#### Behavioral Study Example
```
subject_id | behavior  | value | unit
-----------|-----------|-------|------
mouse_1    | freezing  | 12.5  | s
mouse_1    | jumping   | 5.2   | s
mouse_2    | freezing  | 8.7   | s
mouse_2    | jumping   | 7.1   | s
```

#### HPLC Study Example
```
subject_id | compound  | location    | value | unit
-----------|-----------|-------------|-------|------
rat_1      | serotonin | cortex      | 5.6   | nmol
rat_1      | dopamine  | cerebellum  | 0.8   | nmol
rat_2      | dopamine  | cortex      | 5.8   | nmol
rat_2      | serotonin | cerebellum  | 0.7   | nmol
```

### Important Notes
- Each row represents a single measurement for a specific subject and characteristic combination
- Subjects can have multiple rows (one per unique combination of measurement characteristic columns)
- The combination of `subject_id` and measurement characteristic columns must be unique
- It's recommended to handle or remove null values before importing, as they may cause inconsistencies depending on how they were implemented
- This structure enables flexible data selection and filtering during analysis

## Requests

Requests are used to specify what data to analyze and how to process it. They are passed as a dictionary to plotting functions via the ```request``` parameter.

### Basic Structure
```python
{
    "datasets": {
        "dataset_name": {
            "column_name": "value",  # Filter condition
            "another_column": ["value1", "value2"]  # Multiple values
            "experiment": "experiment_name",  # Uses metadata, only selection here, not used for stats
            "remove_outliers": {"test_name": "mode"}
        }
    },
    "experiment": "experiment_name"  # Both selects the data and is used as a flag that group comparison should be performed (histogram, summary histogram, network summary, stats table)
    # "experiment": "default" will keep all groups and treat them as unpaired, parametric, and with the group being the single independant variable. Used in cases where the experiments actual independant variables are not relevant to stats such as in network summary 
}

```

### Key Components

#### 1. Dataset Selection
- The `datasets` key contains one or more dataset specifications
- Each key in `datasets` should match a dataset name
- Within each dataset specification, key-value pairs act as filters:
  - Keys are column names
  - Values can be single values or lists of values to match
  - Only rows matching all specified conditions will be selected
- The `remove_outliers` key is optional and can be used to remove outliers from the dataset. It takes a dictionary with the following keys:
  - `test_name`: The name of the outlier test to use (currently only grubbs and iqr are supported)
  - `mode`: The mode of the outlier test (e.g. "calculated"). Previously a functionality enabling manual outlier selection existed and will be restored in future versions.
- The program also support the selection of ratios of measurements. For example with {"compound": "dopamine/serotonin"} the program will compute on the fly the ratio of dopamine/serotonin and handle it as any other value. Only "simple" ratios are currently supported, development of complex ratios (eg dopamine in cortex / serotonin in cerebellum) is in progress.

#### 2. Multiple Datasets
When multiple datasets are specified:
- Data is combined into a single dataset
- A new `dataset` column is added to identify the source
- Measurement characteristics are combined into tuples in the `measurement` column
  *Note: This behavior may change in future versions to improve handling of datasets with different structures*

#### 3. Experiment Specification
- The optional `experiment` key enables statistical comparisons
- When present, the system will compute statistical comparisons between groups defined in the experiment's metadata
- Requires an experiment name that matches an entry in the project's metadata

### Example
```python
# Request for behavioral data from control and treated groups
request = {
    "datasets": {
        "behavior_data": {
            "group": ["control", "treated"],
            "timepoint": "day7",
            "remove_outliers": {"grubbs": "calculated"}
        },
        "hplc_data": {
            "compound": ["dopamine", "serotonin"],
            "location": "cortex"
        }
    },
    "experiment": "lsd_dose_response"
}
```

I understand that this is a bit complex, but I think it's the most flexible way to specify what data to analyze. I'm working on a more user-friendly interface to make it easier to create requests.

## Metadata

Metadata orchestrates all aspects of a project, defining datasets, experiments, groups, and visualization parameters. The metadata is stored in an Excel file with multiple sheets, each serving a specific purpose. Users with knowledge in data modeling will notice an effort to make the metadata as close to a database schema as possible. Development of a proper database is in progress.

### 1. Datasets Sheet
Defines the structure of each dataset in the project.

| Column | Type | Description |
|--------|------|-------------|
| `label` | string | Name of the dataset |
| `measurement_columns` | CSV string | Column names that identify unique measurements (e.g., "behavior" or "compound,location") |

### 2. Experiments Sheet
Configures experimental designs and statistical parameters. Statistical test pipeline is automaticcaly selected based on these parameters.

| Column | Type | Description |
|--------|------|-------------|
| `label` | string | Name of the experiment |
| `group_ids` | CSV list | IDs of groups included in this experiment |
| `independant_variables` | bool (0/1) | Whether the experiment has independent variables |
| `paired` | bool (0/1) | Whether groups are paired |
| `parametric` | bool (0/1) | Whether to assume parametric distribution (future versions will auto-detect) |

### 3. Groups Sheet
Defines subject groups and their properties.

| Column | Type | Description |
|--------|------|-------------|
| `label` | string | Name of the group |
| `group_id` | int | Unique identifier for the group |
| `independant_variables` | string | Independent variables in the group |
| `subject_ids` | CSV list | List of subject IDs in this group |

### 4. Palette Sheet
Controls visualization appearance for different groups.

| Column | Type | Description |
|--------|------|-------------|
| `group_id` | int | Reference to group ID |
| `color` | string | Color for bars/points in visualizations |
| `significance_symbol` | string | Symbol used to denote statistical significance |

### 5. Statistics Sheet
Configures default statistical parameters.

| Column | Type | Description |
|--------|------|-------------|
| `max_outliers` | int | Maximum number of outliers to remove |
| `p_value_threshold` | float | Threshold for statistical significance |

### Example Metadata Structure
```yaml
datasets:
  - label: "behavior"
    measurement_columns: "behavior"

experiments:
  - label: "drug_study"
    group_ids: "1,2,3"
    independant_variables: 1
    paired: 0
    parametric: 1

groups:
  - label: "control"
    group_id: 1
    independant_variables: ""
    subject_ids: "1,2,3,4"

palette:
  - group_id: 1
    color: "#1f77b4"
    significance_symbol: "*"

statistics:
  max_outliers: 2
  p_value_threshold: 0.05
```


## Additionnal features (incomplete)

### Constants

There exists a number of files in module.json. These work with the classes in module.core.Registry and are simply mapping that may be used to store information that is more persistent than a single project. Obviously this is a very poor setup and should be part of project metadata, it was just simpler to do it this way. Proper implementation will wait for the database backend.

#### Measurements

These files (regions.json, compounds.json, measures.json) store the "true" names of various measurement characteristics. If a charcteristic exists within a dataset, the system will attempt to validate it using the correspondingly named contant file. This ensures consistent naming across datasets with similar data.

#### Classes

These files store groupings of measurement characteristics. They enable the user to more easily select recurring groups of regions all the while keeping the name of the group in the automatic figure titling and naming.

#### Positions

A special case of classes that include locations. Only used in the network figure to position nodes.


## Future Development

We're actively working on improvements to make this tool more powerful and user-friendly:

- **Multi project figures**: Enabling the combination of data collected in different experiments with automtic normalization of the data to increase statistical power.
- **Python Package**: Converting the codebase into a proper Python package for easier installation and distribution.
- **RESTful API**: Developing a comprehensive API to enable programmatic access to analysis functions.
- **Database Backend**: Implementing a database solution for better data management and querying capabilities.
- **Web Interface**: Creating an intuitive web-based interface to make the tool accessible to non-technical users.

These enhancements will maintain all current functionality while making the tool more robust and easier to use.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
