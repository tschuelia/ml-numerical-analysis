from .utils import get_all_raxml_llhs, get_all_iqtree_llhs
from .custom_types import FilePath


def save_best_tree(
    trees_file: FilePath, logs_file: FilePath, output_file: FilePath, program_name: str
) -> None:
    with open(trees_file) as f:
        all_trees = f.readlines()

    if program_name == "raxml":
        all_llhs = get_all_raxml_llhs(logs_file)
    elif program_name == "iqtree":
        all_llhs = get_all_iqtree_llhs(logs_file)
    else:
        raise ValueError(
            "Either the specified tree inference program is not supported or you did not specify any."
        )

    idx_best = all_llhs.index(max(all_llhs))

    with open(output_file, "w") as f:
        f.write(all_trees[idx_best])
