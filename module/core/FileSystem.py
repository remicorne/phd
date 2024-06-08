import os


class FileSystem:

    ROOT = os.getcwd()
    PROJECTS = f"{ROOT}/PROJECTS"
    CONSTANTS = f"{ROOT}/module/json"
    PATH_ELEMENT_ORDER = ["project", "experiment", "figure_type"]
    
    @staticmethod
    def to_path(**path_elements):
        path = os.path.join(FileSystem.PROJECTS)
        for element in FileSystem.PATH_ELEMENT_ORDER:
            if element in path_elements:
                path = os.path.join(path, path_elements[element])
        return path

    @staticmethod
    def list_projects():
        projects = os.listdir(FileSystem.PROJECTS)
        return [project.split("/")[-1] for project in projects]

    @staticmethod
    def list_experiments(project=None):
        if project:
            experiments = []
            for file in os.listdir(project):
                if os.path.isdir(file):
                    experiments.append(file.split("/")[-1])
            return experiments
        else:
            experiments = {}
            for project in FileSystem.list_projects():
                experiments[project] = FileSystem.list_experiments(project)
            return experiments

    @staticmethod
    def get_location(**path_elements):
        for i, element in enumerate(FileSystem.PATH_ELEMENT_ORDER):
            location = FileSystem.to_path(
                **{
                    element: path_elements[element]
                    for element in FileSystem.PATH_ELEMENT_ORDER[: i + 1]
                    if element in path_elements
                }
            )
            if not os.path.exists(location):
                os.mkdir(location)
                print(f"Created {path_elements[element]} {element} folder")
        return location
