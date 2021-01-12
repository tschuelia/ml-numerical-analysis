from database import IQ_Tree, Raxml_Tree, RFDistance, Run, db
from iqtree_parser import get_iqtree_results
from utils import (
    get_cleaned_rf_dist,
    get_iqtree_llh,
    get_parameter_value,
    get_raxml_abs_rf_distance,
    get_raxml_llh,
    get_raxml_num_unique_topos,
    get_raxml_rel_rf_distance,
)
from typing import List


def create_Run(
    params_file: str,
    raxml_treesearch_log_file: str,
    iqtree_test_log_file: str,
    rfDistances_log_file: str,
) -> Run:
    """
    Creates an Run object in the database and returns it.

    Args:
        params_file: Path to a file containing the parameter settings
        raxml_treesearch_log_file: Path to a file containing the max likelihood
            as found by the raxml run (raxml treesearch log most likely)
        iqtree_test_log_file: Path to a file containing the max likelihood as
            reevaluated by iqtree (iqtree log most likely)
        rfDistances_log_file: Path to a file containing the rf distance
            run summary
    Returns:
        The created run database object

    Raises:
    """
    blmin = get_parameter_value(params_file, "blmin")
    blmax = get_parameter_value(params_file, "blmax")
    raxml_llh = get_raxml_llh(raxml_treesearch_log_file)
    iqtree_llh = get_iqtree_llh(iqtree_test_log_file)

    average_absolute_rf_distance = get_raxml_abs_rf_distance(
        rfDistances_log_file
    )
    average_relative_rf_distance = get_raxml_rel_rf_distance(
        rfDistances_log_file
    )
    unique_topos = get_raxml_num_unique_topos(rfDistances_log_file)

    return Run.create(
        blmin=blmin,
        blmax=blmax,
        final_raxml_llh=raxml_llh,
        best_iqtree_llh=iqtree_llh,
        average_absolute_rf_distance=average_absolute_rf_distance,
        average_relative_rf_distance=average_relative_rf_distance,
        unique_topos=unique_topos,
    )


def create_IQ_Trees(
    run: Run, iqtree_trees_file: str, iqtree_results_file: str
) -> None:
    """
    Creates IQ_Tree objects based on the given files.

    Args:
        run: The database Run object the iqtrees are associated with
        iqtree_trees_file: Path to the file containing all iqtree strings
        iqtree_results_file: Path to the file containing the results of the
            performed statistical tests
    """
    with open(iqtree_trees_file) as f:
        iqtrees = f.readlines()

    iqtrees = [t.strip() for t in iqtrees]

    iqtree_results = get_iqtree_results(iqtree_results_file)

    # instead of saving each tree seperately in the database (slow!)
    # we collect the information...
    data_source = []

    for tree_res in iqtree_results:
        values = {}

        num = tree_res["tree_id"]
        tree = iqtrees[num - 1]

        values["run"] = run
        values["iq_tree"] = tree
        values["llh"] = tree_res["logL"]
        values["deltaL"] = tree_res["deltaL"]

        tests = tree_res["tests"]

        if "bp-RELL" in tests:
            values["bpRell"] = tests["bp-RELL"]["score"]
            values["bpRell_significant"] = tests["bp-RELL"]["significant"]

        if "p-KH" in tests:
            values["pKH"] = tests["p-KH"]["score"]
            values["pKH_significant"] = tests["p-KH"]["significant"]

        if "p-SH" in tests:
            values["pSH"] = tests["p-SH"]["score"]
            values["pSH_significant"] = tests["p-SH"]["significant"]

        if "p-WKH" in tests:
            values["pWKH"] = tests["p-WKH"]["score"]
            values["pWKH_significant"] = tests["p-WKH"]["significant"]

        if "p-WSH" in tests:
            values["pWSH"] = tests["p-WSH"]["score"]
            values["pWSH_significant"] = tests["p-WSH"]["significant"]

        if "c-ELW" in tests:
            values["cELW"] = tests["c-ELW"]["score"]
            values["cELW_significant"] = tests["c-ELW"]["significant"]

        if "p-AU" in tests:
            values["pAU"] = tests["p-AU"]["score"]
            values["pAU_significant"] = tests["p-AU"]["significant"]

        data_source.append(values)

    # ... and then bulk insert all objects at once
    with db.atomic():
        IQ_Tree.insert_many(data_source).execute()


def create_Raxml_Trees(
    run: Run, best_tree_file: str, all_trees_file: str
) -> List[Raxml_Tree]:
    """
    Creates Raxml_Tree objects based on the given files.

    Args:
        run: The database Run object the trees are associated with
        best_tree_file: Path to the file containing the best found raxml tree
        all_trees_file: Path to the file containing all found raxml trees
    Returns:
        A list of the created Raxml_Tree database objects
    """
    # read best tree
    with open(best_tree_file) as f:
        # read only first line as second line is an empty line
        best_tree_str = f.readline().strip()

    # read all trees
    with open(all_trees_file) as f:
        all_trees = f.readlines()

    all_trees = [t.strip() for t in all_trees]

    tree_objects = []
    for tree_str in all_trees:
        # tree is the best tree if the newick string matches exactly
        # the best_tree newick string
        is_best = tree_str == best_tree_str
        t = Raxml_Tree.create(run=run, raxml_tree=tree_str, is_best=is_best)
        tree_objects.append(t)

    return tree_objects


def create_RF_Distances(
    tree_objects: List[Raxml_Tree], rfDistances_file: str
) -> None:
    """
    Creates RF_Distance objects in the database.
    Each RF_Distance is respective two trees and contains the plain distance
        and the normalized distance.

    Args:
        tree_objects: List of Raxml_Tree objects. Required to associate the
            distances to the correct trees.
        rfDistances_file: Path to the file containing all pairwise rf distances

    Returns:
        None

    Raises:
        ValueError: If the rfDistance_file contains tree indices that exceed
            the number of trees given in the tree_objects list.
    """
    with open(rfDistances_file) as f:
        rf_dist_lines = f.readlines()

    data_source = []

    for line in rf_dist_lines:
        values = {}

        (
            tree_idx1,
            tree_idx2,
            plain_dist,
            normalized_dist,
        ) = get_cleaned_rf_dist(line)

        # sanity check whether the tree indices do not exceed the number of tree_objects
        if tree_idx1 >= len(tree_objects) or tree_idx2 >= len(tree_objects):
            raise ValueError(
                f"The tree_object list only contains {len(tree_objects)} trees.\
                    Tried to access trees {tree_idx1}, {tree_idx2}.")

        values["tree1"] = tree_objects[tree_idx1]
        values["tree2"] = tree_objects[tree_idx2]
        values["plain_rf_distance"] = plain_dist
        values["normalized_rf_distance"] = normalized_dist

        data_source.append(values)

    # bulk insert all rf distances
    with db.atomic():
        RFDistance.insert_many(data_source).execute()


def create_RF_Distances_best_trees(
        rfDistances_file: str,
        best_trees_file: str,
        best_tree_objects: List[Raxml_Tree]
) -> None:
    """
    Creates RFDistance objects for the rf distances between the best trees
    for each run setting.

    Args:
        rfDistances_file: Path to the file containing the pairwise rf
            distances for the best trees
        best_trees_file: Path to the file containing all best trees in
            newick format
        best_tree_objects: List of Raxml_Tree database objects
            marked as best tree.

    Returns:
        None

    Raises:
        ValueError:
            - If the rfDistance_file contains tree indices that exceed
                the number of trees given in the best_tree_objects list.
            - If for one tree string in best_trees_file exists no
                matching tree in best_tree_objects
    """
   with open(rfDistances_file) as f:
        rf_dist_lines = f.readlines()

    with open(best_trees_file) as f:
        best_trees = f.readlines()

    """
    For each line in the rfdistance log file:
      1. get the tree indices and the distances
      2. sanity check whether the tree indices do not exceed the number of best trees
      3. using the tree indices: get the newick tree string from the best_trees_file
      4. find the corresponding best_tree database object by comparing the newick tree strings
      5. sanity check whether trees were found
      6. store the data for later bulk insert
    """
    data_source = []

    for line in rf_dist_lines:

        values = {}

        (
            tree_idx1,
            tree_idx2,
            plain_dist,
            normalized_dist,
        ) = get_cleaned_rf_dist(line)

        if tree_idx1 >= len(best_tree_objects) or tree_idx2 >= len(best_tree_objects):
            raise ValueError(
                f"The tree_object list only contains {len(best_tree_objects)} trees.\
                    Tried to access trees {tree_idx1}, {tree_idx2}.")

        tree_str1 = best_trees[tree_idx1].strip()
        tree_str2 = best_trees[tree_idx2].strip()

        tree1 = None
        tree2 = None

        """ 
        In order to reference the correct tree objects in the RFDistance object we need to find the 
        best_tree object for each tree string.
        To do this we compare the newick tree string for all best_trees 
        with the newick string we read from the best_trees_file.
        """
        for tree in best_tree_objects:
            if tree.raxml_tree == tree_str1:
                tree1 = tree
            elif tree.raxml_tree == tree_str2:
                tree2 = tree

        if not tree1 or not tree2:
            raise ValueError(
                f"Search for best trees with indices {tree_idx1} and {tree_idx2} failed:\
                    no corresponding best tree found.")

        values["tree1"] = tree1
        values["tree2"] = tree2
        values["plain_rf_distance"] = plain_dist
        values["normalized_rf_distance"] = normalized_dist

        data_source.append(values)
    
    # bulk insert all rf distances
    with db.atomic():
        RFDistance.insert_many(data_source).execute()


# initialize empty database
db.init(snakemake.output.database)
db.connect()
db.create_tables([Run, Raxml_Tree, IQ_Tree, RFDistance])

# snakemake.input.[something] is a list of filepaths
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

# sanity check whether we got all files
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

    tree_objects = create_Raxml_Trees(
        run, best_tree_files[i], all_trees_raxml_files[i]
    )

    for t in tree_objects:
        if t.is_best:
            best_tree_objects.append(t)

    create_RF_Distances(tree_objects, rfDistances[i])

    create_IQ_Trees(run, iqtree_trees_files[i], iqtree_results_files[i])

create_RF_Distances_best_trees(
    rfDistances_best_trees[0],
    best_trees_collected[0],
    best_tree_objects)
