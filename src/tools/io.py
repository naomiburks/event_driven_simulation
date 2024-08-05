"""This module handles saving and loading from the filesystem."""
# pylint:disable=missing-function-docstring
import json
import csv
from os import listdir

BASE_PATH = "output"
SIMULATION_PATH = "simulations"
PLOT_PATH = "plots"
IGNORED_PREFIX = "."


def read_simulation(filename):
    file_path = f"{BASE_PATH}/{SIMULATION_PATH}/{filename}.json"
    contents = _read_json(file_path)
    # json converts all keys to strings, but we want some to be floats
    # so we have to convert them back.
    contents = _convert_json_keys(contents)
    return contents

def write_simulation(data, filename):
    file_path = f"{BASE_PATH}/{SIMULATION_PATH}/{filename}.json"
    return _write_json(data, file_path)

def get_simulations():
    dir_path = f"{BASE_PATH}/{SIMULATION_PATH}/"
    return _get_tracked_files(dir_path)

# save figure

def save_figure(fig, filename):
    fig.savefig(f"{BASE_PATH}/{PLOT_PATH}/{filename}.svg")








# Read and write JSONs and CSVs

def _read_json(file_path):
    with open(file_path, encoding="utf8") as json_file:
        data = json.load(json_file)
    return data


def _write_json(data, file_path):
    with open(file_path, 'w', encoding="utf8") as outfile:
        json.dump(data, outfile)


def _read_csv(file_path):
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',',
                               quotechar='"', quoting=csv.QUOTE_MINIMAL)
        return [row for row in csvreader]


def _write_csv(file_path, data):
    with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',',
                               quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerows(data)

def _get_tracked_files(dir_path):
    return list(filter(_is_tracked, listdir(dir_path)))

def _is_tracked(file_name):
    return not file_name[0] == IGNORED_PREFIX

def _convert_json_keys(json):
    """Returns a json with string keys converted to floats"""
    if isinstance(json, list):
        return [_convert_json_keys(item) for item in json]
    if isinstance(json, dict):
        converted_dict = {}
        for key, val in json.items():
            try:
                key = int(key)
            except ValueError:
                try: 
                    key = float(key)
                except ValueError:
                    pass
            
            converted_dict[key] = _convert_json_keys(val)
        return converted_dict
    return json



