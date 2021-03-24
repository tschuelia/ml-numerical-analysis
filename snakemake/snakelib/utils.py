import json
from Bio import Phylo
from .custom_types import *


def get_parameter_value(filename: FilePath, param_identifier: str) -> float:
    with open(filename) as f:
        data = json.load(f)

    if param_identifier not in data:
        raise ValueError(
            f"The given parameter identifier {param_identifier} is not stored in the given file {filename}."
        )

    return data[param_identifier]


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


def read_file_contents(file_path: FilePath) -> List[str]:
    with open(file_path) as f:
        content = f.readlines()

    return [l.strip() for l in content]

def get_tree_object(newick_str: Newick) -> Phylo.Newick.Tree:
    trees = list(Phylo.NewickIO.Parser.from_string(newick_str).parse())
    return trees[0]


def get_number_of_taxa_for_tree(newick_str: Newick) -> int:
    tree = get_tree_object(newick_str)
    return tree.count_terminals()


def get_total_branch_length_for_tree(newick_str: Newick) -> float:
    tree = get_tree_object(newick_str)
    return tree.total_branch_length()


def get_average_branch_length_for_tree(newick_str: Newick) -> float:
    total_brlen = get_total_branch_length_for_tree(newick_str)
    num_taxa = get_number_of_taxa_for_tree(newick_str)
    return total_brlen / num_taxa