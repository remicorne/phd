import os


class FileSystem:
    """
    Static class that handles the filesystem architecture.
    Automatically creates folders if they don't exist.
    Automatically extracts the project, experiment and figure_type from the path.

    """

    ROOT = os.path.dirname(os.path.abspath(__file__))
    PROJECTS = os.path.join(ROOT, "PROJECTS")
    CONSTANTS = os.path.join(ROOT, "json")

    PATH_ELEMENT_ORDER = ["project", "experiment", "figure_type"]
    
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
        location = FileSystem.PROJECTS
        for element in FileSystem.PATH_ELEMENT_ORDER:
            if path_elements.get(element):
                location = os.path.join(location, path_elements[element])
                if not os.path.exists(location):
                    os.mkdir(location)
                    print(f"Created {path_elements[element]} {element} folder")
        return location
