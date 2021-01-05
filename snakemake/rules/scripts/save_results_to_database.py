from database import db, Run, Raxml_Tree, IQ_Tree, RFDistance
from utils import *
from iqtree_parser import get_iqtree_results
import regex


def create_Run(
    params_file, raxml_treesearch_log_file, iqtree_test_log_file, rfDistances_log_file
):
    blmin = get_parameter_value(params_file, "blmin")
    blmax = get_parameter_value(params_file, "blmax")
    raxml_llh = get_raxml_llh(raxml_treesearch_log_file)
    iqtree_llh = get_iqtree_llh(iqtree_test_log_file)

    average_absolute_rf_distance = get_raxml_abs_rf_distance(rfDistances_log_file)
    average_relative_rf_distance = get_raxml_rel_rf_distance(rfDistances_log_file)
    unique_topos = get_raxml_num_unique_topos(rfDistances_log_file)

    # store run in database
    return Run.create(
        blmin=blmin,
        blmax=blmax,
        final_raxml_llh=raxml_llh,
        best_iqtree_llh=iqtree_llh,
        average_absolute_rf_distance=average_absolute_rf_distance,
        average_relative_rf_distance=average_relative_rf_distance,
        unique_topos=unique_topos,
    )


def create_IQ_Trees(run, iqtree_trees_file, iqtree_results_file):
    with open(iqtree_trees_file) as f:
        iqtrees = f.readlines()

    iqtrees = [t.strip() for t in iqtrees]

    iqtree_results = get_iqtree_results(iqtree_results_file)

    for tree_res in iqtree_results:
        num = tree_res["tree_id"]
        tree = iqtrees[num - 1]

        iqtree = IQ_Tree.create(
            run=run,
            iq_tree=tree,
            llh=tree_res["logL"],
            deltaL=tree_res["deltaL"],
        )

        tests = tree_res["tests"]

        if "bp-RELL" in tests:
            iqtree.bpRell = tests["bp-RELL"]["score"]
            iqtree.bpRell_significant = tests["bp-RELL"]["significant"]

        if "p-KH" in tests:
            iqtree.pKH = tests["p-KH"]["score"]
            iqtree.pKH_significant = tests["p-KH"]["significant"]

        if "p-SH" in tests:
            iqtree.pSH = tests["p-SH"]["score"]
            iqtree.pSH_significant = tests["p-SH"]["significant"]

        if "p-WKH" in tests:
            iqtree.pWKH = tests["p-WKH"]["score"]
            iqtree.pWKH_significant = tests["p-WKH"]["significant"]

        if "p-WSH" in tests:
            iqtree.pWSH = tests["p-WSH"]["score"]
            iqtree.pWSH_significant = tests["p-WSH"]["significant"]

        if "c-ELW" in tests:
            iqtree.cELW = tests["c-ELW"]["score"]
            iqtree.cELW_significant = tests["c-ELW"]["significant"]

        if "p-AU" in tests:
            iqtree.pAU = tests["p-AU"]["score"]
            iqtree.pAU_significant = tests["p-AU"]["significant"]

        iqtree.save()


def get_cleaned_rf_dist(raw_line):
    line_regex = regex.compile(r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+.\d+)\s*")
    tree_idx1, tree_idx2, plain_dist, normalized_dist = regex.search(
        line_regex, line
    ).groups()
    return int(tree_idx1), int(tree_idx2), float(plain_dist), float(normalized_dist)


db.init(snakemake.output.database)
db.connect()
db.create_tables([Run, Raxml_Tree, IQ_Tree, RFDistance])

params_files = snakemake.input.params_file
best_tree_files = snakemake.input.best_tree_raxml
all_trees_raxml_files = snakemake.input.all_trees_raxml
iqtree_results_files = snakemake.input.iqtree_results
iqtree_trees_files = snakemake.input.iqtree_trees
raxml_treesearch_log_files = snakemake.input.raxml_treesearch_log
iqtree_test_log_files = snakemake.input.iqtree_test_log
rfDistances_log_files = snakemake.input.rfDistances_log
rfDistances = snakemake.input.rfDistances

rfDistances_best_trees = snakemake.input.rfDistances_best_trees
best_trees_collected = snakemake.input.best_trees_collected

num_runs = len(best_tree_files)

assert (
    len(best_tree_files)
    == len(all_trees_raxml_files)
    == len(iqtree_results_files)
    == len(raxml_treesearch_log_files)
    == len(iqtree_test_log_files)
    == len(rfDistances_log_files)
    == len(rfDistances)
)

best_tree_objects = []

for i in range(num_runs):

    run = create_Run(
        params_files[i],
        raxml_treesearch_log_files[i],
        iqtree_test_log_files[i],
        rfDistances_log_files[i],
    )

    # store Trees
    # read raxml_best_tree
    with open(best_tree_files[i]) as f:
        # read only first line as second line is an empty line
        raxml_best_tree = f.readline()

    raxml_best_tree = raxml_best_tree.strip()

    # read raxml_all_trees
    with open(all_trees_raxml_files[i]) as f:
        raxml_all_trees = f.readlines()

    raxml_all_trees = [t.strip() for t in raxml_all_trees]

    tree_objects = []

    for tree in raxml_all_trees:
        is_best = tree == raxml_best_tree
        t = Raxml_Tree.create(
            run=run,
            raxml_tree=tree,
            is_best=is_best,
        )
        tree_objects.append(t)
        if is_best:
            best_tree_objects.append(t)

    with open(rfDistances[i]) as f:
        rf_dist_lines = f.readlines()

    # TODO: see http://docs.peewee-orm.com/en/latest/peewee/querying.html for correct implementation
    for line in rf_dist_lines:
        tree_idx1, tree_idx2, plain_dist, normalized_dist = get_cleaned_rf_dist(line)

        RFDistance.create(
            tree1=tree_objects[tree_idx1],
            tree2=tree_objects[tree_idx2],
            plain_rf_distance=plain_dist,
            normalized_rf_distance=normalized_dist,
        )

    create_IQ_Trees(run, iqtree_trees_files[i], iqtree_results_files[i])

with open(rfDistances_best_trees[0]) as f:
    rf_dist_lines = f.readlines()

with open(best_trees_collected[0]) as f:
    best_trees = f.readlines()


for line in rf_dist_lines:
    tree_idx1, tree_idx2, plain_dist, normalized_dist = get_cleaned_rf_dist(line)

    tree_str1 = best_trees[tree_idx1].strip()
    tree_str2 = best_trees[tree_idx2].strip()

    tree1 = None
    tree2 = None

    # find the matching trees
    for tree in best_tree_objects:
        if tree.raxml_tree == tree_str1:
            tree1 = tree
        elif tree.raxml_tree == tree_str2:
            tree2 = tree

    if not tree1 or not tree2:
        raise ValueError(
            f"Search for best trees with indices {tree_idx1} and {tree_idx2} failed: no corresponding best tree found."
        )

    RFDistance.create(
        tree1=tree1,
        tree2=tree2,
        plain_rf_distance=plain_dist,
        normalized_rf_distance=normalized_dist,
    )
