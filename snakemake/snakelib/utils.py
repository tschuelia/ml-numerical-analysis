import json
import regex
from typing import Tuple
from .custom_types import *


def get_parameter_value(filename: FilePath, param_identifier: str) -> float:
    with open(filename) as f:
        data = json.load(f)

    if param_identifier not in data:
        raise ValueError(
            f"The given parameter identifier {param_identifier} is not stored in the given file {filename}."
        )

    return data[param_identifier]


def get_iqtree_run_param_value(line: str, param_identifier: str) -> float:
    iqtree_run_param_regex = regex.compile(
        rf"/*-{param_identifier}\s+(\d+(?:\.\d+)?(?:[e][-+]?\d+)?)/*"
    )
    value = regex.search(iqtree_run_param_regex, line).groups()[0]
    return float(value)


def get_iqtree_run_param_values_from_file(
    iqtree_log: FilePath, param_identifier: str
) -> TreeIndexed[float]:
    content = read_file_contents(iqtree_log)

    values = []

    for line in content:
        if line.startswith("Command:"):
            values.append(get_iqtree_run_param_value(line, param_identifier))

    return values


def get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f'The given line "{line}" does not contain the search string "{search_string}".'
    )


def get_single_value_from_file(input_file: FilePath, search_string: str) -> float:
    with open(input_file) as f:
        lines = f.readlines()

    for l in lines:
        if search_string in l:
            return get_value_from_line(l, search_string)

    raise ValueError(
        f'The given input file {input_file} does not contain the search string "{search_string}".'
    )


def get_multiple_values_from_file(
    input_file: FilePath, search_string: str
) -> List[float]:
    with open(input_file) as f:
        lines = f.readlines()

    values = []
    for l in lines:
        if search_string in l:
            values.append(get_value_from_line(l, search_string))

    return values


def get_all_iqtree_llhs(iqtree_file: FilePath) -> TreeIndexed[float]:
    STR = "BEST SCORE FOUND :"
    return get_multiple_values_from_file(iqtree_file, STR)


def get_best_iqtree_llh(iqtree_file: FilePath) -> float:
    all_llhs = get_all_iqtree_llhs(iqtree_file)
    return max(all_llhs)


def get_iqtree_wallclock_time(log_file: FilePath) -> TreeIndexed[float]:
    content = read_file_contents(log_file)

    all_times = []

    for line in content:
        if not "Total wall-clock time used:" in line:
            continue

        # correct line looks like this: "Total wall-clock time used: 0.530 sec (0h:0m:0s)"
        value = line.split(" ")[4]
        all_times.append(float(value))

    if not all_times:
        raise ValueError(
            f"The given input file {log_file} does not contain the wall clock time."
        )

    return all_times


def get_iqtree_treesearch_wallclock_time_entire_run(log_file: FilePath) -> float:
    return sum(get_iqtree_wallclock_time(log_file))


def read_file_contents(file_path: FilePath) -> List[str]:
    with open(file_path) as f:
        content = f.readlines()

    return [l.strip() for l in content]