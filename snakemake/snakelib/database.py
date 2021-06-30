import peewee as P

raxml_db = P.SqliteDatabase(None)
iqtree_db = P.SqliteDatabase(None)
fasttree_db = P.SqliteDatabase(None)


class BaseProgram(P.Model):
    # parameter values
    blmin = P.FloatField(null=True)
    blmax = P.FloatField(null=True)
    lh_eps = P.FloatField(null=True)

    num_pars_trees = P.IntegerField()
    num_rand_trees = P.IntegerField()
    best_treesearch_llh = P.FloatField()
    best_evaluation_llh = P.FloatField(null=True)
    treesearch_total_time = P.FloatField()


class Raxmlng(BaseProgram):
    raxml_param_epsilon = P.FloatField()
    branch_length_smoothing = P.IntegerField()
    spr_lh_epsilon = P.FloatField()
    bfgs_factor = P.FloatField()

    avg_abs_rfdist_treesearch = P.FloatField()
    avg_rel_rfdist_treesearch = P.FloatField()
    num_unique_topos_treesearch = P.IntegerField()

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
    eval_blmin = P.FloatField()
    eval_blmax = P.FloatField()
    eval_lh_eps = P.FloatField()


class RaxmlEvalTree(EvalTree):
    eval_raxml_param_epsilon = P.FloatField()
    eval_raxml_brlen_smoothings = P.IntegerField()
    eval_spr_lh_epsilon = P.FloatField()
    eval_bgfs_factor = P.FloatField()

    class Meta:
        database = raxml_db


class IqtreeEvalTree(EvalTree):
    class Meta:
        database = iqtree_db


class BaseRFDistance(P.Model):
    plain_rfdist = P.FloatField()
    normalized_rfdist = P.FloatField()


class RFDistTreesearchTree(BaseRFDistance):
    tree1 = P.ForeignKeyField(TreesearchTree)
    tree2 = P.ForeignKeyField(TreesearchTree)

    class Meta:
        database = raxml_db


class RFDistEvalTree(BaseRFDistance):
    tree1 = P.ForeignKeyField(EvalTree)
    tree2 = P.ForeignKeyField(EvalTree)

    class Meta:
        database = raxml_db

class BaseTreeStatsTest(P.Model):
    # all tests are performed respective this tree
    reference_tree_id = P.IntegerField()
    # tree results
    tree_id = P.IntegerField()
    iqtree_llh = P.FloatField(null=True)
    deltaL = P.FloatField(
        null=True
    )  # llh difference to max llh in the set according to iqtree run

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