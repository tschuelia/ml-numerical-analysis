import sqlite3
import pandas as pd
from MasterThesis.snakemake.snakelib.iqtree_statstest_parser import get_iqtree_results
import subprocess
import os
import shutil


dataset_name = "1288"
alignment = f"/home/ju/Masterarbeit/test-Datasets/DNA-Data/{dataset_name}/{dataset_name}.phy"
model = "GTR+FO+G4"

programs = ["raxml", "iqtree", "fasttree"]
table_names = {
    "raxml": "raxmlevaltree",
    "iqtree": "iqtreeevaltree",
    "fasttree": "fasttreetreesearchtree"
}
base_dir = "./tmp/"


for program in programs:
    print("processing trees for program ", program)

    try:
        shutil.rmtree(f"tmp/{program}")
    except:
        pass
    os.mkdir(f"tmp/{program}")
    os.mkdir(f"tmp/{program}/iqtree_statstest_evaltrees")
    os.mkdir(f"tmp/{program}/iqtree_statstest_evaltrees/iqtree")

    prefix = "./tmp/" + program + "/iqtree_statstest_evaltrees/iqtree/statstest_evaltrees"
    db = base_dir + program + "_results.sqlite3"
    best_tree_path = "./tmp/" + program + "/iqtree_statstest_evaltrees/best_overall_eval_tree.tree"
    all_trees_path = "./tmp/" + program + "/iqtree_statstest_evaltrees/all_eval_trees.trees"

    con = sqlite3.connect(db)
    cur = con.cursor()

    eval_trees = pd.read_sql_query(f"SELECT * FROM {table_names[program]}", con)
    eval_trees = eval_trees.sort_values(by=["llh", "compute_time"], ascending=[False, True])

    best_tree = eval_trees.head(1)["newick_tree"].item()
    best_tree_id = eval_trees.head(1)["id"].item()

    with open(best_tree_path, "w") as f:
        f.write(best_tree)

    trees = list(eval_trees["newick_tree"])

    with open(all_trees_path, "w") as f:
        f.write("\n".join(trees))

    cmd = []
    cmd.append("iqtree")
    cmd.append("-s")
    cmd.append(alignment)
    cmd.append("-m")
    cmd.append(model)
    cmd.append("-pre")
    cmd.append(prefix)
    cmd.append("-redo")
    cmd.append("-z")
    cmd.append(all_trees_path)
    cmd.append("-te")
    cmd.append(best_tree_path)
    cmd.append("-n")
    cmd.append("0")
    cmd.append("-zb")
    cmd.append("10000")
    cmd.append("-zw")
    cmd.append("-au")
    cmd.append("-nt")
    cmd.append(str(4))

    subprocess.check_output(cmd, encoding='utf-8')

    # parse results
    results = get_iqtree_results(prefix + ".iqtree")

    # add EvalTreeStatsTest table to database
    table_name = f"{program}evaltreestatstest"
    columns = "(reference_tree_id, tree_id, deltaL, bpRell, bpRell_significant, pKH, pKH_significant, pSH, pSH_significant, \
     pWKH, pWKH_significant, pWSH, pWSH_significant, cELW, cELW_significant, pAU, pAU_significant)"

    cur.execute(
        f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            reference_tree_id INTEGER NOT NULL,
            tree_id INTEGER NOT NULL,
            deltaL REAL,
            bpRell REAL,
            bpRell_significant INTEGER,
            pKH REAL,
            pKH_significant INTEGER,
            pSH REAL,
            pSH_significant INTEGER,
            pWKH REAL,
            pWKH_significant INTEGER,
            pWSH REAL,
            pWSH_significant INTEGER,
            cELW REAL,
            cELW_significant INTEGER,
            pAU REAL,
            pAU_significant INTEGER
        )
        """
    )
    # make sure the table is empty before saving the new results
    cur.execute(
        f"""DELETE FROM {table_name}"""
    )

    for i, entry in enumerate(results):
        tree_values = {}
        reference_tree_id = best_tree_id
        tree_id = eval_trees.iloc[i]["id"].item()
        tests = entry["tests"]
        deltaL = entry["deltaL"]

        if "bp-RELL" in tests:
            bpRell = tests["bp-RELL"]["score"]
            bpRell_significant = tests["bp-RELL"]["significant"]
        else:
            bpRell = None
            bpRell_significant = None

        if "p-KH" in tests:
            pKH = tests["p-KH"]["score"]
            pKH_significant = tests["p-KH"]["significant"]
        else:
            pKH = None
            pKH_significant = None

        if "p-SH" in tests:
            pSH = tests["p-SH"]["score"]
            pSH_significant = tests["p-SH"]["significant"]
        else:
            pSH = None
            pSH_significant = None

        if "p-WKH" in tests:
            pWKH = tests["p-WKH"]["score"]
            pWKH_significant = tests["p-WKH"]["significant"]
        else:
            pWKH = None
            pWKH_significant = None

        if "p-WSH" in tests:
            pWSH = tests["p-WSH"]["score"]
            pWSH_significant = tests["p-WSH"]["significant"]
        else:
            pWSH = None
            pWSH_significant = None

        if "c-ELW" in tests:
            cELW = tests["c-ELW"]["score"]
            cELW_significant = tests["c-ELW"]["significant"]
        else:
            cELW = None
            cELW_significant = None

        if "p-AU" in tests:
            pAU = tests["p-AU"]["score"]
            pAU_significant = tests["p-AU"]["significant"]
        else:
            pAU = None
            pAU_significant = None

        cur.execute(
            f"""
            INSERT INTO {table_name} {columns}
            VALUES
            {reference_tree_id, tree_id, deltaL, bpRell, bpRell_significant, pKH, pKH_significant, pSH, pSH_significant,
            pWKH, pWKH_significant, pWSH, pWSH_significant, cELW, cELW_significant, pAU, pAU_significant}
            """
        )

    con.commit()
    con.close()
    print("done for ", program)

print("all done")

# import sqlite3
# import pandas as pd
# db = "/Users/julia/Desktop/Masterarbeit/DatenCluster/cascade/1288/results_23Feb21/iqtree_results_all_eval_settings.sqlite3"
# db_new = "/Users/julia/Desktop/Masterarbeit/DatenCluster/cascade/1288/results_23Feb21/iqtree_results.sqlite3"
# con = sqlite3.connect(db)
# con_new = sqlite3.connect(db_new)
#
# df = pd.read_sql_query("SELECT * FROM iqtreeevaltree", con)
# print("before ", len(df.index))
# df = df.loc[
#     (df.eval_blmin == 1e-10) &
#     (df.eval_blmax == 100) &
#     (df.eval_lh_eps == 0.001)
# ]
# print("after ", len(df.index))
# df.to_sql("iqtreeevaltree", con_new)
#
# df = pd.read_sql_query("SELECT * FROM iqtree", con)
# df.to_sql("iqtree", con_new)
#
# df = pd.read_sql_query("SELECT * FROM iqtreetreesearchtree", con)
# df.to_sql("iqtreetreesearchtree", con_new)