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


def saveFigure(fig, identifier, fig_type):
    output_subdir = f"{OUTPUT_DIR}/{fig_type}"
    checkFileSystem(output_subdir)
    fig.savefig(f"{output_subdir}/{identifier}.svg")  # dpi also?
    print(f"SAVED {output_subdir}/{identifier}.svg")
    # https://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
    fig.savefig(f"{output_subdir}/{identifier}.png", dpi=fig.dpi)
    # Image.open(f"{output_subdir}/{identifier}.png").show()
    print(f"SAVED {output_subdir}/{identifier}.png")


def saveCorrelogram(fig, identifier):
    saveFigure(fig, identifier, "correlograms")


def saveHistogram(fig, identifier):
    saveFigure(fig, identifier, "histograms")



def saveQuantitativeSummaryFig(fig, identifier):
    saveFigure(fig, identifier, "quantitative_summary")


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
        if input("EXIT PROGRAM? (y/n)") == "y":
            print("EXITING PROGRAM")
            sys.exit()
    return answer


def askMultipleChoice(question, choices):
    choice_mapping = {f"{i}": choice for i, choice in enumerate(choices)}
    options = "\n".join([f"{i}: {choice}" for i, choice in choice_mapping.items()])
    choice = inputEscape(f"{question}\n{options}\n")
    while choice not in choice_mapping.keys():
        choice = inputEscape(f"Invalid choice, possibilities are:\{options}\n")
    return choice_mapping[choice]


def askSelectParameter(data, column):
    options = set(data[column])
    print(options)
    answer = inputEscape(f"""Select {column}?\n{', '.join(options)}\n""")#.upper()# i cant access vH, dH ect
    while answer not in options:
        print(f".{answer}.")
        answer = inputEscape(
            f"Invalid choice, possibilities are:\n{', '.join(options).upper()}\n"
        ).upper()
    return answer


def askYesorNo(question):
    return inputEscape(f"{question} (y/n)\n").upper() == "Y"


def subselectDf(df, subselect):
    return df.query(
        "&".join([f"{column}=='{value}'" for column, value in subselect.items()])
    )


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
                return ["ERROR", str(e)]

    return wrapper


IDENTIFIERS = {
    "histogram": 'f"{experiment}_for_{compound}_in_{region}"',
    "correlogram": 'f"{experiment}_{correlogram_type}_{buildCorrelogramFilenmae(to_correlate, columns)}"',
    "head_twitch_histogram": 'f"head_twitch_histogram_{experiment}_for_{vairable}"', 
    "quantitative_summary":'f"quantitative_summary_{experiment}_for_{compound}_in_{region}"' #REMI: i tried to add this here but dont understand yet how it works
}


def buildIdentifier(identifier_type, **kwargs):
    locals().update(kwargs)
    return eval(IDENTIFIERS[identifier_type])


# TODO there should be seperate get and add decorators
def get_or_add(identifier_type):
    def decorator(builder_func):
        def wrapper(*args, **kwargs):
            identifier = buildIdentifier(identifier_type, **kwargs)
            identifier=checkIdentifier(identifier)
            from_scratch = (
                kwargs.get("from_scratch")
                if kwargs.get("from_scratch") is not None
                else askYesorNo("Recalculate figure even if previous version exists?")
            )
            filename = args[0]
            if from_scratch or not isCached(filename, identifier):
                result = builder_func(*args, **kwargs)
                cache(args[0], identifier, result)
                saveFigure(result, identifier, identifier_type)
            else:
                result = getCache(filename, identifier)

        return wrapper

    return decorator

#checks for symbols that can be be save in name / 
def checkIdentifier(identifier):
    invalid_chars_pattern = r'[\/:*?"<>|%\&#$@!=+,;\'`~]'
    sanitized_identifier = re.sub(invalid_chars_pattern, '_', identifier)
    if re.search(invalid_chars_pattern, identifier):
        print("Invalid characters in identifier, replacing with '_' ")
    return sanitized_identifier
