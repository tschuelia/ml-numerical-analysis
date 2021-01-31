import json
import regex
from typing import Tuple
from custom_types import *


def get_parameter_value(filename: FilePath, param_identifier: str) -> float:
    with open(filename) as f:
        data = json.load(f)

    if param_identifier not in data:
        raise ValueError(
            f"The given parameter identifier {param_identifier} is not stored in the given file {filename}."
        )

    return data[param_identifier]


def get_raxml_run_param_value(line: str, param_identifier: str) -> float:
    raxml_run_param_regex = regex.compile(
        rf"/*--{param_identifier}\s+(\d+(?:\.\d+)?(?:[e][-+]?\d+)?)/*"
    )
    value = regex.search(raxml_run_param_regex, line).groups()[0]
    return float(value)


def get_raxml_run_param_values_from_file(
    raxml_log: FilePath, raxml_command: str, param_identifier: str
) -> float:
    with open(raxml_log) as f:
        lines = f.readlines()

    values = []

    for l in lines:
        if l.startswith(raxml_command):
            values.append(get_raxml_run_param_value(l, param_identifier))

    return values


def _get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f'The given line "{line}" does not contain the search string "{search_string}".'
    )


def _get_single_value_from_file(input_file: FilePath, search_string: str) -> float:
    with open(input_file) as f:
        lines = f.readlines()

    for l in lines:
        if search_string in l:
            return _get_value_from_line(l, search_string)

    raise ValueError(
        f'The given input file {input_file} does not contain the search string "{search_string}".'
    )


def _get_multiple_values_from_file(
    input_file: FilePath, search_string: str
) -> List[float]:
    with open(input_file) as f:
        lines = f.readlines()

    values = []
    for l in lines:
        if search_string in l:
            values.append(_get_value_from_line(l, search_string))

    return values


def get_all_raxml_seeds(raxml_file: FilePath) -> TreeIndexed[int]:
    STR = "random seed:"
    return _get_multiple_values_from_file(raxml_file, STR)


def get_all_raxml_llhs(raxml_file: FilePath) -> TreeIndexed[float]:
    STR = "Final LogLikelihood:"
    return _get_multiple_values_from_file(raxml_file, STR)


def get_best_raxml_llh(raxml_file: FilePath) -> float:
    all_llhs = get_all_raxml_llhs(raxml_file)
    return max(all_llhs)


def get_all_iqtree_llhs(iqtree_file: FilePath) -> TreeIndexed[float]:
    STR = "BEST SCORE FOUND :"
    return _get_multiple_values_from_file(iqtree_file, STR)


def get_best_iqtree_llh(iqtree_file: FilePath) -> float:
    all_llhs = get_all_iqtree_llhs(iqtree_file)
    return max(all_llhs)


def get_raxml_abs_rf_distance(log_file: FilePath) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return _get_single_value_from_file(log_file, STR)


def get_raxml_rel_rf_distance(log_file: FilePath) -> float:
    STR = "Average relative RF distance in this tree set:"
    return _get_single_value_from_file(log_file, STR)


def get_raxml_num_unique_topos(log_file: FilePath) -> int:
    STR = "Number of unique topologies in this tree set:"
    return _get_single_value_from_file(log_file, STR)


def get_raxml_treesearch_elapsed_time(log_file: FilePath) -> TreeIndexed[float]:
    with open(log_file) as f:
        content = f.readlines()

    all_times = []

    for line in content:
        if not "Elapsed time:" in line:
            continue

        # correct line looks like this: "Elapsed time: 1976.056 seconds"
        value = line.split(" ")[2]
        all_times.append(float(value))

    if not all_times:
        raise ValueError(
            f"The given input file {log_file} does not contain the elapsed time."
        )

    return all_times


def get_raxml_treesearch_elapsed_time_entire_run(log_file: FilePath) -> float:
    return sum(get_raxml_treesearch_elapsed_time(log_file))


def get_cleaned_rf_dist(raw_line: str) -> Tuple[int, int, float, float]:
    line_regex = regex.compile(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d+)\s*")
    tree_idx1, tree_idx2, plain_dist, normalized_dist = regex.search(
        line_regex, raw_line
    ).groups()
    return int(tree_idx1), int(tree_idx2), float(plain_dist), float(normalized_dist)
