import os
import shutil


class FileSystem:
    """
    Static class that handles the filesystem architecture.
    Automatically creates folders if they don't exist.
    Automatically extracts the project, experiment and figure_type from the path.

    """

    ROOT = os.getcwd()
    PROJECTS = f"{ROOT}/PROJECTS"
    CONSTANTS = f"{ROOT}/module/json"
    PATH_ELEMENT_ORDER = ["project", "experiment", "figure_type"]

    @staticmethod
    def project_exists(project_name):
        return os.path.exists(os.path.join(FileSystem.PROJECTS, project_name))

    @staticmethod
    def list_datasets(project):
        if os.path.exists(f"{FileSystem.PROJECTS}/{project}"):
            files = os.listdir(f"{FileSystem.PROJECTS}/{project}")
            return [file.split("/")[-1] for file in files if file.endswith(".pkl")]
        return []

    @staticmethod
    def get_location(**path_elements):
        location = FileSystem.PROJECTS
        if not os.path.exists(location):
            os.mkdir(location)
        for element in FileSystem.PATH_ELEMENT_ORDER:
            if path_elements.get(element):
                location = os.path.join(location, path_elements[element])
                if not os.path.exists(location):
                    os.mkdir(location)
                    print(f"Created {path_elements[element]} {element} folder")
        if path_elements.get("filename"):
            location = os.path.join(location, path_elements["filename"])
        return location

    @staticmethod
    def delete_project(project):
        project_folder = os.path.join(FileSystem.PROJECTS, project)
        shutil.rmtree(project_folder)
