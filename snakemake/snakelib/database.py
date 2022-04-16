import peewee as P
import uuid

from .program_parser import Program
from .custom_types import TreeIndexed, Dict

raxml_db = P.SqliteDatabase(None)
iqtree_db = P.SqliteDatabase(None)
fasttree_db = P.SqliteDatabase(None)


class BaseProgram(P.Model):
    # parameter values
    uuid = P.UUIDField()
    lh_epsilon = P.FloatField(null=True)

    num_pars_trees = P.IntegerField(null=True)
    num_rand_trees = P.IntegerField(null=True)
    best_treesearch_llh = P.FloatField()
    best_evaluation_llh = P.FloatField(null=True)
    treesearch_total_time = P.FloatField()


class Raxmlng(BaseProgram):
    class Meta:
        database = raxml_db


class Iqtree(BaseProgram):
    class Meta:
        database = iqtree_db


class Fasttree(BaseProgram):
    class Meta:
        database = fasttree_db


class BaseTree(P.Model):
    uuid = P.UUIDField()
    llh = P.FloatField()
    compute_time = P.FloatField()
    newick_tree = P.CharField()
    is_best = P.BooleanField(null=True)
    number_of_taxa = P.IntegerField(null=True)
    total_branch_length = P.FloatField(null=True)
    average_branch_length = P.FloatField(null=True)


class TreesearchTree(BaseTree):
    program = P.ForeignKeyField(BaseProgram)
    program_uuid = P.UUIDField()
    seed = P.FloatField()


class RaxmlTreesearchTree(TreesearchTree):
    class Meta:
        database = raxml_db


class IqtreeTreesearchTree(TreesearchTree):
    class Meta:
        database = iqtree_db


class FasttreeTreesearchTree(TreesearchTree):
    class Meta:
        database = fasttree_db


class EvalTree(BaseTree):
    start_tree = P.ForeignKeyField(TreesearchTree)
    start_tree_uuid = P.UUIDField()

class RaxmlEvalTree(EvalTree):
    class Meta:
        database = raxml_db


class RaxmlEvalStartingTree(EvalTree):
    class Meta:
        database = raxml_db


class IqtreeEvalTree(EvalTree):
    class Meta:
        database = iqtree_db


class BaseTreeStatsTest(P.Model):
    # id of the respective eval tree
    tree = P.ForeignKeyField(EvalTree)
    tree_uuid = P.UUIDField()
    plausible = P.BooleanField()
    # id of the cluster (unique topologies cluster; this is useful for debugging)
    cluster_id = P.IntegerField()

    bpRell = P.FloatField(null=True)  # bootstrap proportion using RELL method.
    # True denotes the 95% confidence sets.
    # False denotes significant exclusion.
    bpRell_significant = P.BooleanField(null=True)

    # p-value of one sided Kishino-Hasegawa test (1989).
    pKH = P.FloatField(null=True)
    pKH_significant = P.BooleanField(null=True)

    # p-value of Shimodaira-Hasegawa test (2000).
    pSH = P.FloatField(null=True)
    pSH_significant = P.BooleanField(null=True)

    pWKH = P.FloatField(null=True)  # p-value of weighted KH test.
    pWKH_significant = P.BooleanField(null=True)

    pWSH = P.FloatField(null=True)  # p-value of weighted SH test.
    pWSH_significant = P.BooleanField(null=True)

    cELW = P.FloatField(null=True)  # Expected Likelihood Weight
    cELW_significant = P.BooleanField(null=True)

    # p-value of approximately unbiased (AU) test.
    pAU = P.FloatField(null=True)
    pAU_significant = P.BooleanField(null=True)


class RaxmlEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = raxml_db


class RaxmlEvalAndStartingTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = raxml_db


class RaxmlPairwiseEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        pass
    database = raxml_db


class RaxmlPairwiseEvalAndStartingTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = raxml_db


class IqtreeEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = iqtree_db


class IQTreePairwiseEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = iqtree_db




class FasttreeEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = fasttree_db


def insert_program_data(program: Program, database_table: BaseProgram) -> BaseProgram:
    return database_table.create(
        uuid=uuid.uuid4().hex,
        lh_epsilo=program.lh_epsilon,

        num_pars_trees=program.num_pars_trees,
        num_rand_trees=program.num_rand_trees,
        best_treesearch_llh=program.best_treesearch_llh,
        best_evaluation_llh=program.best_evaluation_llh,
        treesearch_total_time=program.treesearch_total_time,
    )


def insert_treesarch_data(
        program: Program,
        program_database_object: BaseProgram,
        database_table: TreesearchTree
) -> (TreeIndexed[TreesearchTree], TreesearchTree):
    """
    This function inserts all treesearch trees corresponding to the given Program object

    Arguments:
        program: Program object that contains the entire data (runs, treesearch, eval)
        database_table: Database table where the data should be inserted into
        program_database_object: BaseProgram database object the trees should be references to (= belong to)

    Returns:
        The list of database objects of the inserted treesearch trees and the tree marked as best tree
    """
    treesarch_objects = []
    best_treesearch_object = None
    for tree_idx in range(program.get_number_of_trees()):
        is_best = (
                program.tree_for_index_is_best(tree_idx)
                and not best_treesearch_object
        )

        tree = database_table.create(
            uuid=uuid.uuid4().hex,
            llh=program.treesearch_llhs[tree_idx],
            compute_time=program.treesearch_compute_times[tree_idx],
            newick_tree=program.treesearch_trees[tree_idx].newick_str,
            is_best=is_best,
            number_of_taxa=program.treesearch_trees[tree_idx].number_of_taxa,
            total_branch_length=program.treesearch_trees[tree_idx].total_branch_length,
            average_branch_length=program.treesearch_trees[tree_idx].average_branch_length,
            program=program_database_object,
            program_uuid=program_database_object.uuid,
            seed=program.treeseach_seeds[tree_idx],
        )
        treesarch_objects.append(tree)

        if is_best:
            best_treesearch_object = tree

    assert best_treesearch_object
    return treesarch_objects, best_treesearch_object


def insert_eval_data(
        program: Program,
        start_tree_database_objects: TreeIndexed[TreesearchTree],
        database_table: EvalTree
) -> (TreeIndexed[EvalTree], EvalTree):
    """
    This function inserts all eval trees corresponding to the given Program object

    Arguments:
        program: Program object that contains the entire data (runs, treesearch, eval)
        start_tree_database_object: Reference to the corresponding treesearch tree
        database_table: Database table where the data should be inserted into

    Returns:
        The list of database objects of the inserted eval trees and the tree marked as best tree
    """
    best_eval_tree = None
    eval_trees = []
    for eval_tree_idx in range(program.get_number_of_eval_trees()):
        is_best = program.eval_tree_for_index_is_best(eval_tree_idx)

        # fmt: off
        eval_tree = database_table.create(
            uuid                        = uuid.uuid4().hex,
            start_tree                  = start_tree_database_objects[eval_tree_idx],
            start_tree_uuid             = start_tree_database_objects[eval_tree_idx].uuid,
            llh                         = program.eval_llhs[eval_tree_idx],
            newick_tree                 = program.eval_trees[eval_tree_idx].newick_str,
            compute_time                = program.eval_compute_times[eval_tree_idx],
            is_best                     = is_best,
            number_of_taxa              = program.eval_trees[eval_tree_idx].number_of_taxa,
            total_branch_length         = program.eval_trees[eval_tree_idx].total_branch_length,
            average_branch_length       = program.eval_trees[eval_tree_idx].average_branch_length,
        )
        # fmt: on

        eval_trees.append(eval_tree)

        if is_best:
            best_eval_tree = eval_tree

    return eval_trees, best_eval_tree


def insert_starting_eval_data(
        program: Program,
        start_tree_database_objects: TreeIndexed[TreesearchTree],
        database_table: EvalTree
) -> (TreeIndexed[EvalTree], EvalTree):
    """
    This function inserts all eval trees corresponding to the given Program object

    Arguments:
        program: Program object that contains the entire data (runs, treesearch, eval, starting_eval)
        start_tree_database_objects: Reference to the corresponding treesearch tree
        database_table: Database table where the data should be inserted into

    Returns:
        The list of database objects of the inserted starting eval trees and the tree marked as best tree
    """

    eval_trees = []
    for eval_tree_idx in range(program.get_number_of_eval_trees()):

        # fmt: off
        eval_tree = database_table.create(
            uuid            = uuid.uuid4().hex,
            start_tree      = start_tree_database_objects[eval_tree_idx],
            start_tree_uuid = start_tree_database_objects[eval_tree_idx].uuid,
            llh             = program.starting_eval_llhs[eval_tree_idx],
            newick_tree     = program.starting_eval_trees[eval_tree_idx].newick_str,
            compute_time    = program.starting_eval_compute_times[eval_tree_idx],
        )
        # fmt: on

        eval_trees.append(eval_tree)

    return eval_trees


def insert_statstest_data(
        eval_trees: TreeIndexed[EvalTree],
        statstest_results: TreeIndexed[Dict],
        cluster_ids: TreeIndexed[int],
        database_table: BaseTreeStatsTest,
):
    for i, eval_tree in enumerate(eval_trees):
        tests = statstest_results[i]["tests"]
        statstest_values = {}
        statstest_values["tree"] = eval_tree
        statstest_values["tree_uuid"] = eval_tree.uuid
        statstest_values["cluster_id"] = cluster_ids[i]
        statstest_values["plausible"] = statstest_results[i]["plausible"]

        if "bp-RELL" in tests:
            statstest_values["bpRell"] = tests["bp-RELL"]["score"]
            statstest_values["bpRell_significant"] = tests["bp-RELL"]["significant"]

        if "p-KH" in tests:
            statstest_values["pKH"] = tests["p-KH"]["score"]
            statstest_values["pKH_significant"] = tests["p-KH"]["significant"]

        if "p-SH" in tests:
            statstest_values["pSH"] = tests["p-SH"]["score"]
            statstest_values["pSH_significant"] = tests["p-SH"]["significant"]

        if "p-WKH" in tests:
            statstest_values["pWKH"] = tests["p-WKH"]["score"]
            statstest_values["pWKH_significant"] = tests["p-WKH"]["significant"]

        if "p-WSH" in tests:
            statstest_values["pWSH"] = tests["p-WSH"]["score"]
            statstest_values["pWSH_significant"] = tests["p-WSH"]["significant"]

        if "c-ELW" in tests:
            statstest_values["cELW"] = tests["c-ELW"]["score"]
            statstest_values["cELW_significant"] = tests["c-ELW"]["significant"]

        if "p-AU" in tests:
            statstest_values["pAU"] = tests["p-AU"]["score"]
            statstest_values["pAU_significant"] = tests["p-AU"]["significant"]

        database_table.create(**statstest_values)


