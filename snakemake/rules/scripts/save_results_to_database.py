from database import db, Run, Raxml_Tree, IQ_Tree
from utils import get_iqtree_llh, get_raxml_llh, get_parameter_value
from iqtree_parser import get_iqtree_results

db.init(snakemake.output.database)
db.connect()
db.create_tables([Run, Raxml_Tree, IQ_Tree])

params_files = snakemake.input.params_file
best_trees_files = snakemake.input.all_trees_raxml
all_trees_raxml_files = snakemake.input.all_trees_raxml
iqtree_results_files = snakemake.input.iqtree_results
iqtree_trees_files = snakemake.input.iqtree_trees
raxml_treesearch_log_files = snakemake.input.raxml_treesearch_log
iqtree_test_log_files = snakemake.input.iqtree_test_log

num_runs = len(best_trees_files)

assert (
    len(best_trees_files)
    == len(all_trees_raxml_files)
    == len(iqtree_results_files)
    == len(raxml_treesearch_log_files)
    == len(iqtree_test_log_files)
)
for i in range(num_runs):
    blmin = get_parameter_value(params_files[i], "blmin")
    blmax = get_parameter_value(params_files[i], "blmax")
    raxml_llh = get_raxml_llh(raxml_treesearch_log_files[i])
    iqtree_llh = get_iqtree_llh(iqtree_test_log_files[i])

    # store run in database
    run = Run.create(
        blmin=blmin,
        blmax=blmax,
        final_raxml_llh=raxml_llh,
        best_iqtree_llh=iqtree_llh,
    )

    # store Trees
    # read raxml_best_tree
    with open(best_trees_files[i]) as f:
        raxml_best_tree = f.read()

    raxml_best_tree.strip("\n").strip()

    Raxml_Tree.create(
        run=run,
        raxml_tree=raxml_best_tree,
        is_best=True,
    )

    # read raxml_all_trees
    with open(all_trees_raxml_files[i]) as f:
        raxml_all_trees = f.readlines()

    raxml_all_trees = [t.strip("\n").strip() for t in raxml_all_trees]

    for tree in raxml_all_trees:
        Raxml_Tree.create(
            run=run,
            raxml_tree=tree,
            is_best=False,
        )

    # read iqtrees
    with open(iqtree_trees_files[i]) as f:
        iqtrees = f.readlines()

    iqtrees = [t.strip("\n").strip() for t in iqtrees]

    iqtree_results = get_iqtree_results(iqtree_results_files[i])

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