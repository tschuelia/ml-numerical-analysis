import dataclasses
import json
from Bio import Phylo
import numpy as np
from .custom_types import *
from tempfile import TemporaryDirectory
import subprocess
import regex
import pickle

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


def get_cleaned_rf_dist(raw_line: str) -> Tuple[int, int, float, float]:
    line_regex = regex.compile(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s*")
    tree_idx1, tree_idx2, plain_dist, normalized_dist = regex.search(
        line_regex, raw_line
    ).groups()
    return int(tree_idx1), int(tree_idx2), float(plain_dist), float(normalized_dist)


def read_rfdistances(
        rfdistances_file_path,
):
    with open(rfdistances_file_path) as f:
        rfdistances = f.readlines()

    rfdistances = [l.strip() for l in rfdistances if l]

    # based on the number of lines in the rfdistances file
    # we can compute the number of trees that were compared
    # for n trees there are n(n-1)/2 unique pairs
    # so based on this, we can compute n as
    # [1 + sqrt(1 + 8 * len(rfdistances))] / 2
    size = (1 + np.sqrt(1 + 8 * len(rfdistances))) / 2
    assert int(size) == size # this needs to be an int, otherwise we computed something wrong
    size = int(size)
    abs_res = np.zeros((size, size))

    for line in rfdistances:
        # TODO: das hier als numpy matrix -> OOM
        idx1, idx2, plain, norm = get_cleaned_rf_dist(line)
        abs_res[idx1][idx2] = plain
        abs_res[idx2][idx1] = plain

    return abs_res


def raxml_rfdist(
        tree_strings,
        raxml_command,
        get_num_topos=False,
        get_rfdist=False,
        get_pairwise_rfdists=False
):
    with TemporaryDirectory() as tmp_dir:
        trees = tmp_dir + "/trees"

        with open(trees, "w") as f:
            f.write("\n".join(tree_strings))

        cmd = [
            raxml_command,
            "--rfdist",
            trees,
            "--redo"
        ]

        subprocess.check_output(cmd, encoding='utf-8', stderr=subprocess.STDOUT)

        log_pth = trees + ".raxml.log"
        rfdist_pth = trees + ".raxml.rfDistances"

        num_topos = get_raxml_num_unique_topos(log_pth) if get_num_topos else None
        abs_rfdist = get_raxml_abs_rf_distance(log_pth) if get_rfdist else None
        rel_rfdist = get_raxml_rel_rf_distance(log_pth) if get_rfdist else None
        abs_pairwise = read_rfdistances(rfdist_pth) if get_pairwise_rfdists else None

    return num_topos, abs_rfdist, rel_rfdist, abs_pairwise


def get_rfdist_clusters(rfdistances, trees):
    # instead of indexing via the eval_tree.id
    # use the newick string
    # this is slower but more robust
    clusters = []
    size_x, size_y = rfdistances.shape
    for t1 in range(size_x):
        for t2 in range(size_y):
            if t1 >= t2:
                continue
            dist = rfdistances[t1][t2]
            tree1 = trees[t1].strip()
            tree2 = trees[t2].strip()
            seen_t1 = False
            seen_t2 = False
            for s in clusters:
                if tree1 in s:
                    seen_t1 = True
                    if dist == 0:
                        s.add(tree2)
                        seen_t2 = True

                if tree2 in s:
                    seen_t2 = True
                    if dist == 0:
                        s.add(tree1)
                        seen_t1 = True

            if not seen_t1:
                if dist == 0:
                    clusters.append({tree1, tree2})
                    seen_t1 = True
                    seen_t2 = True
                else:
                    clusters.append({tree1})
                    seen_t1 = True

            if not seen_t2:
                clusters.append({tree2})
                seen_t2 = True

    # remove duplicates
    removed_duplicates = []
    for s in clusters:
        if s not in removed_duplicates:
            removed_duplicates.append(s)

    # check if sets are disjoint
    union = set().union(*removed_duplicates)
    n = sum(len(s) for s in removed_duplicates)
    assert n == len(union)

    return removed_duplicates


def filter_tree_topologies(
        raxml_command: str,
        eval_trees: List[NewickString],
        filtered_trees_path: FilePath,
        clusters_path: FilePath
):
    """
    Helper method to filter out duplicate tree topologies in the eval_trees.
    """
    num_trees = len(eval_trees)

    if num_trees > 1:
        num_topos, _, _, abs_pairwise = raxml_rfdist(
            eval_trees,
            raxml_command,
            get_num_topos=True,
            get_rfdist=False,
            get_pairwise_rfdists=True
        )

        if num_topos == 1:
            clusters = [set(eval_trees)]

        else:
            clusters = get_rfdist_clusters(abs_pairwise, eval_trees)

    else:
        clusters = [set(eval_trees)]
        num_topos = 1

    # for each cluster: keep only one tree as representative of the cluster
    unique_trees = [next(iter(cluster)) for cluster in clusters]

    # sanity checks
    assert len(clusters) == num_topos
    assert len(unique_trees) == num_topos
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
