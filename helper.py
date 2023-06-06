#!/usr/bin/env python3
"""Generic helper functions, primarily for I/O tasks such as pickle and JSON file handlings."""

## IMPORTS ##
# External
import json
import matplotlib.pyplot as plt
import os
import pickle
import zipfile
from typing import Any, List
from zipfile import ZipFile


def json_load(path: str) -> Any:
    """Loads the given JSON file and returns it as json_data (a list
    or a dictionary).

    Arguments
    ----------
    * path: str ~ The path of the JSON file
    """
    with open(path) as f:
        json_data = json.load(f)
    return json_data


def json_write(path: str, json_data: Any) -> None:
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written
    * json_data: Any ~ The dictionary or list which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(json_data, indent=4)
    with open(path, "w+", encoding="utf-8") as f:
        f.write(json_output)


def json_zip_load(path: str) -> Any:
    """Loads the given zipped JSON file and returns it as json_data (a list
    or a dictionary).

    Arguments
    ----------
    * path: str ~ The path of the JSON file without ".zip" at the end
    """
    old_wd = os.getcwd()
    folder = os.path.dirname(path)
    filename = os.path.basename(path)
    os.chdir(folder)

    with ZipFile(filename+".zip", 'r') as zip:
        zip.extractall()

    with open(filename) as f:
        json_data = json.load(f)

    os.remove(filename)
    os.chdir(old_wd)

    return json_data


def json_zip_write(path: str, json_data: Any, zip_method: int=zipfile.ZIP_LZMA) -> None:
    """Writes a zipped JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path: str ~  The path of the JSON file that shall be written without ".zip" at the end
    * json_data: Any ~ The dictionary or list which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(json_data, indent=4)
    with open(path, "w+", encoding="utf-8") as f:
        f.write(json_output)

    old_wd = os.getcwd()
    folder = os.path.dirname(path)
    filename = os.path.basename(path)
    os.chdir(folder)
    with ZipFile(filename+".zip", 'w', zip_method) as zip:
        zip.write(filename)
    os.chdir(old_wd)

    os.remove(path)


def ensure_folder_existence(folder: str) -> None:
    """Checks if the given folder exists. If not, the folder is created.

    Argument
    ----------
    * folder: str ~ The folder whose existence shall be enforced.
    """
    if os.path.isdir(folder):
        return
    os.makedirs(folder)


def ensure_json_existence(path: str) -> None:
    if os.path.isfile(path):
        return
    with open(path, "w") as f:
        f.write("{}")


def get_files(path: str) -> List[str]:
    """Returns the names of the files in the given folder as a list of strings.

    Arguments
    ----------
    * path: str ~ The path to the folder of which the file names shall be returned
    """
    files: List[str] = []
    for (_, _, filenames) in os.walk(path):
        files.extend(filenames)
    return files


def pickle_load(path: str) -> Any:
    """Returns the value of the given pickle file.
    Arguments
    ----------
    * path: str ~ The path to the pickle file.
    """
    pickle_file = open(path, 'rb')
    pickled_object = pickle.load(pickle_file)
    pickle_file.close()
    return pickled_object


def pickle_write(path: str, pickled_object: Any) -> None:
    """Writes the given object as pickled file with the given path
    Arguments
    ----------
    * path: str ~ The path of the pickled file that shall be created
    * pickled_object: Any ~ The object which shall be saved in the pickle file
    """
    pickle_file = open(path, 'wb')
    pickle.dump(pickled_object, pickle_file)
    pickle_file.close()


def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder
