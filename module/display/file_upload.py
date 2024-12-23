import ipywidgets as widgets
from IPython.display import display
import pandas as pd
from io import BytesIO
import easygui


def file_upload(description, extensions, callback):
    upload = widgets.FileUpload(
        description=description, accept=",".join(extensions), multiple=False
    )
    display(upload)

    def handle_upload(change):

        if upload.value:
            try:
                uploaded_file = upload.value[0]
                content = uploaded_file["content"]
                callback(content, uploaded_file["name"])
                upload.close()
            except Exception as e:
                print(f"Error processing file: {e}")

    upload.observe(handle_upload, names="value")


def upload_dataframe(description, callback):

    def process_file(content, filename):
        if filename.endswith(".xlsx"):
            df = pd.read_excel(BytesIO(content))
        elif filename.endswith(".pkl"):
            df = pd.read_pickle(BytesIO(content))
        else:
            raise ValueError(f"Unsupported file type: {filename}")
        callback(df)

    file_upload(description, [".xlsx", ".pkl"], process_file)


def select_file():
    file_path = easygui.fileopenbox(
        title="Select a file",
        filetypes=["*.xlsx", "*.pkl"],
    )

    if not file_path:
        raise ValueError("No file selected")

    return file_path


def upload_dataset():
    file_path = select_file()
    if file_path.endswith(".xlsx"):
        df = pd.read_excel(file_path)
    elif file_path.endswith(".pkl"):
        df = pd.read_pickle(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_path}")

    return df
