import json
import pandas as pd

def save_best_tree_and_log(tree_file, all_llhs, all_runtimes, out_tree, out_log):
    with open(tree_file) as f:
        all_trees = f.readlines()

    idx_best = all_llhs.index(max(all_llhs))

    with open(out_tree, "w") as f:
        f.write(all_trees[idx_best])

    data = {
        "llh": all_llhs[idx_best],
        "runtime": all_runtimes[idx_best]
    }

    with open(out_log, "w") as f:
        json.dump(data, f)


def save_overall_best_tree(tree_files, log_files, out_file):
    trees = []
    logs = []

    for tf in tree_files:
        with open(tf) as f:
            trees.append(f.readlines()[0])

    for l in log_files:
        with open(l) as f:
            logs.append(json.load(f))

    df = pd.DataFrame(logs)
    df["newick"] = trees
    df = df.sort_values(by=["llh", "runtime"], ascending=[False, True])

    with open(out_file, "w") as f:
        f.write(df.head(1)["newick"].item())