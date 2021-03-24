from snakelib.utils import (
    get_number_of_taxa_for_tree,
    get_total_branch_length_for_tree,
    get_average_branch_length_for_tree
)

import sqlite3

program = "fasttree"

db_path = f"/Users/julia/Desktop/Masterarbeit/DatenCluster/cascade/2000/results_18Mar21/{program}_results.sqlite3"
table_names = [f"{program}treesearchtree",]# f"{program}evaltree"]
columns = [
    ("number_of_taxa", "INTEGER"),
    ("total_branch_length", "REAL"),
    ("average_branch_length", "REAL")
]

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

con = sqlite3.connect(db_path)
con.row_factory = dict_factory

cur = con.cursor()

for table in table_names:
    for name, type in columns:
        cur.execute(f"ALTER TABLE {table} ADD COLUMN {name} {type}")


    update_data = []

    cur.execute(f"SELECT * FROM {table}")
    for row in cur:
        id = row["id"]
        newick_str = row["newick_tree"]
        num_taxa = get_number_of_taxa_for_tree(newick_str)
        total_brlen = get_total_branch_length_for_tree(newick_str)
        avg_brlen = get_average_branch_length_for_tree(newick_str)
        update_data.append((id, num_taxa, total_brlen, avg_brlen))
        print("computing for id in table ", id, table)

    for id, num_taxa, total_brlen, avg_brlen in update_data:
        cur.execute(f"UPDATE {table} SET number_of_taxa = {num_taxa} WHERE id = {id}")
        cur.execute(f"UPDATE {table} SET total_branch_length = {total_brlen} WHERE id = {id}")
        cur.execute(f"UPDATE {table} SET average_branch_length = {avg_brlen} WHERE id = {id}")
        print("updating for id in table ", id, table)


con.commit()
con.close()
