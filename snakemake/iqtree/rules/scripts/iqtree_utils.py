import regex

from snakelib.utils import (
    get_multiple_values_from_file,
    read_file_contents,
)
from snakelib.custom_types import *


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

    # we need to detect checkpoints, because IQ-Tree logs the command again after resuming
    encountered_checkpoint = False

    for i, line in enumerate(content):
        # check for checkpoint:
        if (line.startswith("*****") and
                i < len(content) - 1 and  # check that we can look one line ahead
                content[i + 1].startswith("CHECKPOINT")):  # check that next line says CHECKPOINT
            encountered_checkpoint = True

        if line.startswith("Command:"):
            if encountered_checkpoint:
                # If the command is after the checkpoint call, do not store the value
                encountered_checkpoint = False
                continue
            else:
                values.append(get_iqtree_run_param_value(line, param_identifier))

    return values


def get_all_iqtree_llhs(iqtree_file: FilePath) -> TreeIndexed[float]:
    STR = "BEST SCORE FOUND :"
    return get_multiple_values_from_file(iqtree_file, STR)


def get_best_iqtree_llh(iqtree_file: FilePath) -> float:
    all_llhs = get_all_iqtree_llhs(iqtree_file)
    return max(all_llhs)


def get_iqtree_cpu_time(log_file: FilePath) -> TreeIndexed[float]:
    content = read_file_contents(log_file)

    all_times = []

    for line in content:
        if not "Total CPU time used:" in line:
            continue

        # correct line looks like this: "Total CPU time used: 0.530 sec (0h:0m:0s)"
        value = line.split(" ")[4]
        all_times.append(float(value))

    if not all_times:
        raise ValueError(
            f"The given input file {log_file} does not contain the cpu time."
        )

    return all_times


def get_iqtree_treesearch_cpu_time_entire_run(log_file: FilePath) -> float:
    return sum(get_iqtree_cpu_time(log_file))
