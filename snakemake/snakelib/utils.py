import dataclasses
import json
from Bio import Phylo
import numpy as np
from .custom_types import *
from tempfile import TemporaryDirectory
import subprocess
import pickle
import tqdm.contrib.concurrent

from .iqtree_statstest_parser import get_iqtree_results, get_default_entry


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
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )


def get_single_value_from_file(input_file: FilePath, search_string: str) -> float:
    with open(input_file) as f:
        lines = f.readlines()

    for l in lines:
        if search_string in l:
            return get_value_from_line(l, search_string)

    raise ValueError(
        f"The given input file {input_file} does not contain the search string '{search_string}'."
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


@dataclasses.dataclass
class NewickTree:
    newick_str: NewickString
    number_of_taxa: int
    total_branch_length: float
    average_branch_length: float
    stdev_branch_length: float
    all_branch_lengths: List[float]
    min_branch_length: float
    max_branch_length: float


def parse_newick_string(newick_string: NewickString) -> NewickTree:
    trees = list(Phylo.NewickIO.Parser.from_string(newick_string).parse())
    tree = trees[0]

    num_taxa = tree.count_terminals()
    total_brlen = tree.total_branch_length()
    all_brlens = [node.branch_length for node in tree.find_clades(branch_length=True)]

    return NewickTree(
        newick_str=newick_string.strip(),
        number_of_taxa=num_taxa,
        total_branch_length=total_brlen,
        average_branch_length=total_brlen / num_taxa,
        stdev_branch_length=np.std(all_brlens),
        all_branch_lengths=all_brlens,
        min_branch_length=min(all_brlens),
        max_branch_length=max(all_brlens)
    )


def cat_input_to_output(input_files, output_file):
    output = []

    for tree_file in input_files:
        with open(tree_file) as f:
            content = f.readlines()
            content = [l.strip() for l in content]
            content = [l for l in content if l]
            output += content

    with open(output_file, "w") as f:
        f.write("\n".join(output))


def get_raxml_abs_rf_distance(log_file: FilePath) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_rel_rf_distance(log_file: FilePath) -> float:
    STR = "Average relative RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxml_num_unique_topos(log_file: FilePath) -> int:
    STR = "Number of unique topologies in this tree set:"
    return get_single_value_from_file(log_file, STR)


def raxml_rfdist(
        tree_strings,
        raxml_command,
        log_path,
):
    with TemporaryDirectory() as tmp_dir:
        trees = tmp_dir + "/trees"

        with open(trees, "w") as f:
            f.write("\n".join(tree_strings))

        cmd = [
            raxml_command,
            "--rfdist",
            trees,
            "--nofiles",
            "--redo",
        ]

        out = subprocess.check_output(cmd, encoding='utf-8')

        with open(log_path, "w")as f:
            f.write(out)

        num_topos = get_raxml_num_unique_topos(log_path)
        abs_rfdist = get_raxml_abs_rf_distance(log_path)

        return num_topos, abs_rfdist


def get_rfdist_clusters(log_path, all_trees):
    content = read_file_contents(log_path)
    clusters = []

    for line in content:
        line = line.strip()
        if not line.startswith("["):
            continue
        # the line is a string representation of a list of ints
        # like this: [1, 2, 3, 4, ]
        tree_ids = eval(line)
        cluster = set()

        for id in tree_ids:
            cluster.add(all_trees[id])

        clusters.append(cluster)

    return clusters



def filter_tree_topologies(
        raxml_command: str,
        eval_trees: List[NewickString],
        filtered_trees_path: FilePath,
        clusters_path: FilePath,
        log_path: FilePath,
):
    """
    Helper method to filter out duplicate tree topologies in the eval_trees.
    """
    num_trees = len(eval_trees)

    if num_trees > 1:
        num_topos, abs_rfdist = raxml_rfdist(
            eval_trees,
            raxml_command,
            log_path,
        )

        if num_topos > 1:
            clusters = get_rfdist_clusters(log_path, eval_trees)

            assert len(clusters) == num_topos
        else:
            clusters = [set(eval_trees)]

    else:
        clusters = [set(eval_trees)]
        num_topos = 1

    # for each cluster: keep only one tree as representative of the cluster
    unique_trees = [next(iter(cluster)) for cluster in clusters]

    # sanity checks
    #assert len(clusters) == num_topos
    #assert len(unique_trees) == num_topos
    assert sum([len(s) for s in clusters]) <= num_trees

    open(filtered_trees_path, "w").write("\n".join(unique_trees))

    with open(clusters_path, "wb") as f:
        pickle.dump(clusters, f)


def get_iqtree_results_for_eval_tree_str(iqtree_results, eval_tree_str, clusters):
    # returns the results for this eval_tree_id as well as the cluster ID
    for i, cluster in enumerate(clusters):
        if eval_tree_str in cluster:
            return iqtree_results[i], i

    raise ValueError("This newick_string belongs to no cluster. newick_str: ", eval_tree_str[:100])


def run_pairwise_iqtree(arg):
    iqtree_command, tree, best_tree, msa, model = arg
    if tree == best_tree:
        return get_default_entry()
    with TemporaryDirectory() as tmpdir:
        fp = tmpdir + "/trees"
        best = tmpdir + "/best"
        log_path = fp + ".significance.iqtree"

        with open(fp, "w") as f:
            f.write(f"{best_tree.strip()}\n{tree.strip()}")

        with open(best, "w") as f:
            f.write(best_tree.strip())

        cmd = [
            iqtree_command,
            "-s",
            msa,
            "-m",
            model,
            "-pre",
            fp + ".significance",
            "-redo",
            "-z",
            fp,
            "-te",
            best,
            "-n",
            "0",
            "-zb",
            "1000",
            "-zw",
            "-au",
            "-nt",
            "2",
            "-seed",
            "0",
            "-safe"
        ]
        subprocess.check_output(cmd, encoding='utf-8')

        results = get_iqtree_results(log_path)
        assert len(results) == 2, "Length of pairwise iqtree comparison results does not match."

        return results[1]


def pairwise_iqtree(iqtree_command, filtered_trees, best_tree, msa, model, max_workers):
    argument_list = [(iqtree_command, tree, best_tree, msa, model) for tree in filtered_trees]
    all_results = tqdm.contrib.concurrent.thread_map(run_pairwise_iqtree, argument_list, max_workers=max_workers, total=len(filtered_trees))
    return all_results
