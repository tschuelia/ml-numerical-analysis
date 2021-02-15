import regex
from typing import Tuple

from snakelib.utils import (
    get_multiple_values_from_file,
    get_single_value_from_file,
    read_file_contents,
)
from snakelib.custom_types import *


def get_raxml_run_param_value(line: str, param_identifier: str) -> float:
    raxml_run_param_regex = regex.compile(
        rf"/*--{param_identifier}\s+(\d+(?:\.\d+)?(?:[e][-+]?\d+)?)/*"
    )
    value = regex.search(raxml_run_param_regex, line).groups()[0]
    return float(value)


def get_raxml_run_param_values_from_file(
    raxml_log: FilePath, raxml_command: str, param_identifier: str
) -> TreeIndexed[float]:
    lines = read_file_contents(raxml_log)

    values = []

    for l in lines:
        if l.startswith(raxml_command):
            values.append(get_raxml_run_param_value(l, param_identifier))

    return values


def get_all_raxml_seeds(raxml_file: FilePath) -> TreeIndexed[int]:
    STR = "random seed:"
    return get_multiple_values_from_file(raxml_file, STR)


def get_all_raxml_llhs(raxml_file: FilePath) -> TreeIndexed[float]:
    STR = "Final LogLikelihood:"
    return get_multiple_values_from_file(raxml_file, STR)


def get_best_raxml_llh(raxml_file: FilePath) -> float:
    all_llhs = get_all_raxml_llhs(raxml_file)
    return max(all_llhs)


def get_raxml_abs_rf_distance(log_file: FilePath) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_rel_rf_distance(log_file: FilePath) -> float:
    STR = "Average relative RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_num_unique_topos(log_file: FilePath) -> int:
    STR = "Number of unique topologies in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_elapsed_time(log_file: FilePath) -> TreeIndexed[float]:
    content = read_file_contents(log_file)

    all_times = []

    for line in content:
        if not "Elapsed time:" in line:
            continue
        # two cases now:
        # either the run was cancelled an rescheduled
        if "restarts" in line:
            # line looks like this: "Elapsed time: 5562.869 seconds (this run) / 91413.668 seconds (total with restarts)"
            _, right = line.split("/")
            value = right.split(" ")[1]
            all_times.append(float(value))

        # ...or the run ran in one sitting...
        else:
            # line looks like this: "Elapsed time: 63514.086 seconds"
            value = line.split(" ")[2]
            all_times.append(float(value))

    if not all_times:
        raise ValueError(
            f"The given input file {log_file} does not contain the elapsed time."
        )

    return all_times


def get_raxml_treesearch_elapsed_time_entire_run(log_file: FilePath) -> float:
    return sum(get_raxml_elapsed_time(log_file))


def get_cleaned_rf_dist(raw_line: str) -> Tuple[int, int, float, float]:
    line_regex = regex.compile(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d+)\s*")
    tree_idx1, tree_idx2, plain_dist, normalized_dist = regex.search(
        line_regex, raw_line
    ).groups()
    return int(tree_idx1), int(tree_idx2), float(plain_dist), float(normalized_dist)


def read_rfdistances(
    rfdistances_file_path: FilePath,
) -> TreeTreeIndexed:
    rfdistances = read_file_contents(rfdistances_file_path)

    res = {}

    for line in rfdistances:
        idx1, idx2, plain, norm = get_cleaned_rf_dist(line)
        res[(idx1, idx2)] = res[(idx2, idx1)] = (plain, norm)

    return res