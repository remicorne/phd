#Check filesystem is set up for write operations
import itertools
import json
import os
import pickle
import shutil
import warnings

from module.constants import (CACHE_DIR,
                              INPUT_DIR,
                              OUTPUT_DIR)

# 
def saveMetadata(filename, treatment_mapping, experimental_info, region_subclassification):
    subcache_dir = f"{CACHE_DIR}/{filename.split('.')[0]}"
    checkFileSystem(subcache_dir)
    saveJSON(f"{subcache_dir}/treatment_mapping.json", treatment_mapping)
    print(f"TREATMENT MAPPING {treatment_mapping} SAVED TO {subcache_dir} SUBCACHE")
    saveJSON(f"{subcache_dir}/experimental_info.json", experimental_info)
    print(f"EXPERIMENTAL INFO {experimental_info} SAVED TO {subcache_dir} SUBCACHE")
    saveJSON(f"{subcache_dir}/region_subclassification.json", region_subclassification)
    print(f"REGION SUBCLASSIFICATION {region_subclassification} SAVED TO {subcache_dir} SUBCACHE")


# This function saves dictionnaries, JSON is a dictionnary text format that you use to not have to reintroduce dictionnaries as variables


def saveJSON(path, dict_to_save):
    with open(path, "w", encoding="utf8") as json_file:
        json.dump(dict_to_save, json_file)


    
#This function gets JSON files and makes them into python dictionnaries
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
    print(f"CREATED {cache_subdir}/{identifier}.pkl CACHE")


# This function gets the dataframes that are cached


def getCache(filename, identifier):
    filename = filename.split(".")[0]
    print(f'GETTING "{identifier}" FROM "{filename}" CACHE')
    with open(f"{CACHE_DIR}/{filename}/{identifier}.pkl", "rb") as file:
        return pickle.load(file)


# This checks if a particulat dataframe/dataset is cached, return boolean


def isCached(filename, identifier):
    filename = filename.split(".")[0]
    return os.path.isfile(f"{CACHE_DIR}/{filename}/{identifier}.pkl")


def saveFigure(fig, identifier, fig_type):
    output_subdir = f"{OUTPUT_DIR}/{fig_type}"
    checkFileSystem(output_subdir)
    fig.savefig(f"{output_subdir}/{identifier}.svg")  # dpi also?
    print(f"SAVED {output_subdir}/{identifier}.svg")
    # https://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
    fig.savefig(f"{output_subdir}/{identifier}.png", dpi=fig.dpi)
    print(f"SAVED {output_subdir}/{identifier}.png")


def saveCorrelogram(fig, identifier):
    saveFigure(fig, identifier, "correlograms")


def saveHistogram(fig, identifier):
    saveFigure(fig, identifier, "histograms")


def saveHeadTwitchHistogram(fig, identifier):
    saveFigure(fig, identifier, "head_twitch_histograms")


def saveQuantitativeSummaryFig(fig, identifier):
    saveFigure(fig, identifier, 'QuantitativeSummaryFigs')



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
    
def flatten(two_d_list):
    return list(itertools.chain.from_iterable(two_d_list))



def askMultipleChoice(question, choices):
    return choices[
        int(
            input(
                question
                + "\n"
                + "\n".join([f"{i}: {choice}" for i, choice in choices.items()])
                + "\n"
            )
        )
    ]



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
        try:
            # call the function with the args it needs
            return ["OK", stat_function(*selected_args)]
        except Warning as w:
            warnings.resetwarnings()
            return ["WARNING", stat_function(*selected_args)]
        except Exception as e:
            return ["ERROR", str(e)]

    return wrapper
