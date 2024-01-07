# Check filesystem is set up for write operations
import itertools
import json
import os
import pickle
import shutil
import sys
import warnings
from module.constants import *
from PIL import Image
import re
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import threading
from IPython.display import Image, display


# This function saves dictionnaries, JSON is a dictionnary text format that you use to not have to reintroduce dictionnaries as variables
def saveJSON(path, dict_to_save):
    with open(path, "w", encoding="utf8") as json_file:
        json.dump(dict_to_save, json_file)


# This function gets JSON files and makes them into python dictionnaries
def getJSON(path):
    with open(path) as outfile:
        loaded_json = json.load(outfile)
    return loaded_json


def checkFileSystem(filepath):
    if not os.path.exists(filepath):
        os.mkdir(filepath)


# This checks that the filesystem has all the requisite folders (input, cache, etc..) and creates them if not


def initiateFileSystem():
    checkFileSystem(INPUT_DIR)
    checkFileSystem(CACHE_DIR)
    checkFileSystem(OUTPUT_DIR)


# This function deletes all cached files, it is used when you want to start from square one because all intermediary results will be cached


def resetCache():
    shutil.rmtree(CACHE_DIR)
    os.mkdir(CACHE_DIR)
    print("CACHE CLEARED")


# This function cahces (aka saves in a easily readable format) all dataframes used
def cache(filename, identifier, to_cache):
    filename = filename.split(".")[0]
    cache_subdir = f"{CACHE_DIR}/{filename}"
    checkFileSystem(cache_subdir)
    with open(f"{cache_subdir}/{identifier}.pkl", "wb") as file:
        pickle.dump(to_cache, file)
    print(f"CACHED {cache_subdir}/{identifier}.pkl")


# This function gets the dataframes that are cached


def getCache(filename, identifier):
    filename = filename.split(".")[0]
    print(f'RETRIEVED "{identifier}" FROM "{filename}" CACHE')
    with open(f"{CACHE_DIR}/{filename}/{identifier}.pkl", "rb") as file:
        return pickle.load(file)


# This checks if a particulat dataframe/dataset is cached, return boolean


def isCached(filename, identifier):
    filename = filename.split(".")[0]
    return os.path.isfile(f"{CACHE_DIR}/{filename}/{identifier}.pkl")


def save_figure(fig, filepath):
    def target():
        fig.savefig(f"{filepath}.svg")  # dpi also?
        fig.savefig(f"{filepath}.png", dpi=fig.dpi)
        print(f"SAVED {filepath}.svg")
        print(f"SAVED {filepath}.png")
    threading.Thread(target=target).start()



def dictToFilename(dict_to_stringify):
    result = str(dict_to_stringify)
    for replacement in [
        ["{", ""],
        ["}", ""],
        [",", "_"],
        [":", ":"],
        ["'", ""],
        ["[", ""],
        ["]", ""],
        [" ", ""],
    ]:
        # This syntaxt wil unpack the list as if I had written 'result.replace(replacement[0], replacement[1])'
        result = result.replace(*replacement)
    return result


def buildCorrelogramFilenmae(to_correlate, columns):
    return f"{to_correlate.replace('/', ':')}_{(',').join(columns)}"


def flatten(two_d_list):
    return list(itertools.chain.from_iterable(two_d_list))


def inputEscape(question):
    answer = input(question)
    if answer == "":
        if askYesorNo("EXIT PROGRAM?"):
            print("EXITING PROGRAM")
            sys.exit()
    return answer


def askMultipleChoice(question, choices):
    if len(choices) == 1:
        return list(choices)[0]
    choice_mapping = {f"{i}": choice for i, choice in enumerate(choices)}
    options = "\n".join([f"{i}: {choice}" for i, choice in choice_mapping.items()])
    choice = inputEscape(f"{question}\n{options}\n")
    while choice not in choice_mapping.keys():
        choice = inputEscape(f"""Invalid choice, possibilities are:\{options}\n""")
    return choice_mapping[choice]


def askSelectParameter(data, column):
    options = set(data[column])
    print(options)
    answers = (
        inputEscape(f"""Select {column}?\n{', '.join(options)}\n""")
        .replace(" ", "")
        .split(",")
    )
    for i, answer in enumerate(answers):
        while answer not in options:
            answer = inputEscape(
                f"Invalid choice: '{answer}', possibilities are:\n{', '.join(options)}\n"
            )
            answers[i] = answer
    return delistify(answers)


def askYesorNo(question):
    answer = inputEscape(f"""{question} (y/n)\n""").upper()
    while answer not in ["Y", "N"]:
        answer = inputEscape(f"""Invalid choice, possibilities are: (y/n)\n""").upper()
    return answer == "Y"


def maskDf(df, mask_conditions):
    complex_filter = True
    for column, value in mask_conditions.items():
        if isinstance(value, list):
            atomic_filter = df[column].isin(value)
        else:
            atomic_filter = df[column] == value
        complex_filter &= atomic_filter
    return complex_filter


def subselectDf(df, subselection):
    df = df[maskDf(df, subselection)]
    return df


# This is a decorator design pattern. Its basically a function that wraps another function and does some operations before
# Calling the wrapped function and returning the result. In this cas it extracts the names of the parameter for the stat function called
# This makes is possible to generically call stat functions and pass them every variable we have all the while not knowing exactly
# Which one(s) they need. If we don't do this, we will be rejected when passing unexpected arguments to the function
# The decorator will look at the function definition to get only the relevant ones from the kwargs variable (keyword args)
# We could also have used **kwargs in all the stat functions and done the parameter selection after that but the decorator saves codelines and also its cool
def select_params(stat_function):
    def wrapper(*args, **kwargs):
        warnings.filterwarnings(
            "error"
        )  # This is to have warning raise exceptions as if they were errors in order to capture and store them
        # This here selects the name of the arguments the function need (it uses internal python variable, rtm if necessary)
        stat_function_args = stat_function.__code__.co_varnames[
            : stat_function.__code__.co_argcount
        ]
        # Select the argument based on their name from kwargs, kwargs is all the keyword arguments passed to the function including the useless ones
        selected_args = [kwargs[arg_name] for arg_name in stat_function_args]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                result = stat_function(*selected_args)
                if w and issubclass(w[-1].category, Warning):
                    print("WARNING")
                    print(w[-1].message)
                return result
            except Exception as e:
                print(e)
                return ["ERROR", str(e)]

    return wrapper


IDENTIFIERS = {
    "histogram": 'f"{experiment}_for_{compound}_in_{region}"',
    "correlogram": 'f"{experiment}_{correlogram_type}_{buildCorrelogramFilenmae(to_correlate, columns)}"',
    "head_twitch_histogram": 'f"head_twitch_histogram_{experiment}_for_{vairable}"',
    "percentage_vehicles": 'f"percentage_vehicles_{experiment}_for_{compound}_in_{regions}"',
    "quantitative_summary": 'f"quantitative_summary_{experiment}_for_{to_plot}_in_{columns}"',
    "pca": 'f"pca_{experiment}_for_{compounds}_in_{regions}"',
    "network": 'f"graph_{experiment}_{correlogram_type}_{buildCorrelogramFilenmae(to_correlate, columns)}"',
}


def buildIdentifier(identifier_type, **kwargs):
    locals().update(kwargs)
    return eval(IDENTIFIERS[identifier_type])


# TODO there should be seperate get and add decorators
def figure_cache(identifier_type):
    def decorator(builder_func):
        def wrapper(*args, **kwargs):
            identifier = buildIdentifier(identifier_type, **kwargs)
            identifier = checkIdentifier(identifier)
            filepath = f"{OUTPUT_DIR}/{identifier}"
            from_scratch = (
                kwargs.get("from_scratch")
                if kwargs.get("from_scratch") is not None
                else askYesorNo("Recalculate figure even if previous version exists?")
            )
            filename = args[0]
            if from_scratch or not isCached(filename, identifier):
                result = builder_func(*args, **kwargs)
                plt.show()
                save_figure(result, filepath)
            else:
                
                display(Image(filename=f'{filepath}.png'))

        return wrapper

    return decorator


# checks for symbols that can be be save in name /
def checkIdentifier(identifier):
    invalid_chars_pattern = r'[\/:*?"<>|%\&#$@!=+,;\'`~]'
    sanitized_identifier = re.sub(invalid_chars_pattern, "_", identifier)
    if re.search(invalid_chars_pattern, identifier):
        print("Invalid characters in identifier, replacing with '_' ")
    return sanitized_identifier


def listify(var):
    return var if isinstance(var, list) else [var]


def delistify(var):
    return var[0] if isinstance(var, list) and len(var) == 1 else var


def parallel_process(executor, cases):
    """Executes the operations performed by {executor} once per case in cases in parralel

    Args:
        executor (function): functio
        cases (list(args)): different arguments to be used by executor

    Returns:
        [executor(case) for case in cases]
    """
    with Pool() as pool:
        results = pool.starmap(executor, cases)
    return results


###### Generic Plotters
def generate_figure(experiment_data):
    """
    Generic function to create subplots of correct dimentions
    input: experimental data listed by treatment, plotter function that takes single treatment data
    ~optional_experimental_info may be passed such that the plotter_cb may scal axis the same for instance
    output: plotted and saved figure at experimental level
    """
    # determin number of treatments to corrispond to number of subplots
    num_treatments = len(experiment_data)
    num_cols = min(int(np.sqrt(num_treatments)), 2)  # max of 2 columns
    num_rows = (num_treatments + num_cols - 1) // num_cols  # Compute the number of rows

    # define the base size and a scaling factor for the figure size
    base_size = 11
    scale_factor = 1

    # create subplots
    fig, axs = plt.subplots(
        num_rows,
        num_cols,
        figsize=(
            num_cols * base_size * scale_factor,
            num_rows * base_size * scale_factor,
        ),
    )
    axs = axs.flatten()
    fig.tight_layout()
    
    return fig, axs

