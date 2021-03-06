import peewee as P

from .program_parser import Program
from .custom_types import TreeIndexed, Dict

raxml_db = P.SqliteDatabase(None)
iqtree_db = P.SqliteDatabase(None)
fasttree_db = P.SqliteDatabase(None)


class BaseProgram(P.Model):
    # parameter values
    blmin = P.FloatField(null=True)
    blmax = P.FloatField(null=True)
    lh_epsilon = P.FloatField(null=True)
    model_epsilon = P.FloatField(null=True)
    branch_length_smoothing = P.IntegerField(null=True)
    spr_lh_epsilon = P.FloatField(null=True)
    bfgs_factor = P.FloatField(null=True)

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
    llh = P.FloatField()
    compute_time = P.FloatField()
    newick_tree = P.CharField()
    is_best = P.BooleanField()
    number_of_taxa = P.IntegerField(null=True)
    total_branch_length = P.FloatField(null=True)
    average_branch_length = P.FloatField(null=True)


class TreesearchTree(BaseTree):
    program = P.ForeignKeyField(BaseProgram)
    seed = P.FloatField()
    starting_tree_type = P.CharField(choices=[("parsimony", "parsimony"), ("random", "random")], default="parsimony")


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
    eval_blmin = P.FloatField(null=True)
    eval_blmax = P.FloatField(null=True)
    eval_lh_epsilon = P.FloatField(null=True)
    eval_model_epsilon = P.FloatField(null=True)
    eval_raxml_brlen_smoothings = P.IntegerField(null=True)
    eval_bfgs_factor = P.FloatField(null=True)


class RaxmlEvalTree(EvalTree):
    class Meta:
        database = raxml_db


class IqtreeEvalTree(EvalTree):
    class Meta:
        database = iqtree_db


class BaseTreeStatsTest(P.Model):
    # id of the respective eval tree
    tree_id = P.IntegerField()
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


class IqtreeEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = iqtree_db


class FasttreeEvalTreeStatsTest(BaseTreeStatsTest):
    class Meta:
        database = fasttree_db


def insert_program_data(program: Program, database_table: BaseProgram) -> BaseProgram:
    return database_table.create(
        blmin=program.blmin,
        blmax=program.blmax,
        lh_epsilon=program.lh_epsilon,
        model_epsilon=program.model_epsilon,
        branch_length_smoothing=program.branch_length_smoothing,
        spr_lh_epsilon=program.spr_lh_epsilon,
        bfgs_factor=program.bfgs_factor,

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
            llh=program.treesearch_llhs[tree_idx],
            compute_time=program.treesearch_compute_times[tree_idx],
            newick_tree=program.treesearch_trees[tree_idx].newick_str,
            is_best=is_best,
            number_of_taxa=program.treesearch_trees[tree_idx].number_of_taxa,
            total_branch_length=program.treesearch_trees[tree_idx].total_branch_length,
            average_branch_length=program.treesearch_trees[tree_idx].average_branch_length,
            program=program_database_object,
            seed=program.treeseach_seeds[tree_idx],
            starting_tree_type=program.starting_tree_types[tree_idx] if program.starting_tree_types else "parsimony"
        )
        treesarch_objects.append(tree)

        if is_best:
            best_treesearch_object = tree

    assert best_treesearch_object
    return treesarch_objects, best_treesearch_object


def insert_eval_data(
        program: Program,
        start_tree_database_objects: TreeIndexed[TreesearchTree],
        database_table: EvalTree) -> (TreeIndexed[EvalTree], EvalTree):
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
            start_tree                  = start_tree_database_objects[eval_tree_idx],
            llh                         = program.eval_llhs[eval_tree_idx],
            newick_tree                 = program.eval_trees[eval_tree_idx].newick_str,
            compute_time                = program.eval_compute_times[eval_tree_idx],
            is_best                     = is_best,
            number_of_taxa              = program.eval_trees[eval_tree_idx].number_of_taxa,
            total_branch_length         = program.eval_trees[eval_tree_idx].total_branch_length,
            average_branch_length       = program.eval_trees[eval_tree_idx].average_branch_length,
            eval_blmin                  = program.eval_blmins[eval_tree_idx] if program.eval_blmins else None,
            eval_blmax                  = program.eval_blmaxs[eval_tree_idx] if program.eval_blmaxs else None,
            eval_lh_epsilon             = program.eval_lh_epsilons[eval_tree_idx] if program.eval_lh_epsilons else None,
            eval_model_epsilon          = program.eval_model_param_epsilons[eval_tree_idx] if program.eval_model_param_epsilons else None,
            eval_raxml_brlen_smoothings = program.eval_raxml_brlen_smoothings[eval_tree_idx] if program.eval_raxml_brlen_smoothings else None,
            eval_bfgs_factor            = program.eval_bfgs_factors[eval_tree_idx] if program.eval_bfgs_factors else None
        )
        # fmt: on

        eval_trees.append(eval_tree)

        if is_best:
            best_eval_tree = eval_tree

    return eval_trees, best_eval_tree


def insert_statstest_data(
        eval_trees: TreeIndexed[EvalTree],
        statstest_results: TreeIndexed[Dict],
        cluster_ids: TreeIndexed[int],
        database_table: BaseTreeStatsTest,
):
    for i, eval_tree in enumerate(eval_trees):
        tests = statstest_results[i]["tests"]
        statstest_values = {}
        statstest_values["tree_id"] = eval_tree
        statstest_values["cluster_id"] = cluster_ids[i]

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


