import regex

from snakelib.utils import (
    read_file_contents,
)
from snakelib.custom_types import *

def get_llh_from_line(line: str) -> float:
    line_re = r"\s*Gamma\s*.*LogLk\s*=\s*([-+]?\d+(?:\.\d+)?(?:[e][-+]?\d+)?).*"
    line_regex = regex.compile(line_re)
    m = regex.match(line_regex, line)

    if not m:
        raise ValueError("The given line does not contain a llh that matches the defined regex.")

    value = m[1]
    return float(value)

def get_all_fasttree_llhs(fasttree_file: FilePath) -> TreeIndexed[float]:
    content = read_file_contents(fasttree_file)
    llhs = []
    for line in content:
        if not line.startswith("Gamma"):
            continue
        llhs.append(get_llh_from_line(line))
    return llhs


def get_best_fasttree_llh(fasttree_log: FilePath):
    all_llhs = get_all_fasttree_llhs(fasttree_log)
    return max(all_llhs)


def get_fasttree_runtimes(fasttree_log: FilePath) -> TreeIndexed[float]:
    content = read_file_contents(fasttree_log)
    line_re = r"Total time:\s*(\d+(?:\.\d+)?(?:[e][-+]?\d+)?).*"
    line_regex = regex.compile(line_re)

    runtimes = []
    for line in content:
        if not "Total time" in line:
            continue
        m = regex.match(line_regex, line)

        if not m:
            raise ValueError("The given line does not contain a runtime that matches the defined regex.")

        value = m[1]
        runtimes.append(float(value))

    return runtimes


def get_fasttree_treesearch_entire_run(fasttree_log: FilePath):
    return sum(get_fasttree_runtimes(fasttree_log))


def get_fasttree_run_param_values_from_file(fasttree_log: FilePath, param_identifier: str) -> TreeIndexed[float]:
    content = read_file_contents(fasttree_log)
    values = []
    for line in content:
        if not line.startswith(param_identifier):
            continue
        _, value = line.split(":")
        value = value.strip()
        values.append(float(value))
    return values