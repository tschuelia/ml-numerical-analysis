import peewee as P

db = P.SqliteDatabase(None)


class BaseModel(P.Model):
    class Meta:
        database = db


class Run(BaseModel):
    blmin = P.FloatField()
    blmax = P.FloatField()


class BaseProgram(BaseModel):
    run = P.ForeignKeyField(Run)
    num_pars_trees = P.IntegerField()
    num_rand_trees = P.IntegerField()
    best_treesearch_llh = P.FloatField()
    best_evaluation_llh = P.FloatField()
    treesearch_total_time = P.FloatField()


class Raxmlng(BaseProgram):
    avg_abs_rfdist_treesearch = P.FloatField()
    avg_rel_rfdist_treesearch = P.FloatField()
    num_unique_topos_treesearch = P.IntegerField()


class Iqtree(BaseProgram):
    pass


class BaseTree(BaseModel):
    llh = P.FloatField()
    compute_time = P.FloatField()
    newick_tree = P.CharField()
    is_best = P.BooleanField()


class TreesearchTree(BaseTree):
    program = P.ForeignKeyField(BaseProgram)
    seed = P.FloatField()


class RaxmlTreesearchTree(TreesearchTree):
    iqtree_llh = P.FloatField()

    # the following are the results of the iqtree run
    # achieved by comparing with best raxml treesearch tree
    deltaL = (
        P.FloatField()
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


class IqtreeTreesearchTree(TreesearchTree):
    pass


class EvalTree(BaseTree):
    start_tree = P.ForeignKeyField(TreesearchTree)
    eval_blmin = P.FloatField()
    eval_blmax = P.FloatField()


class RaxmlEvalTree(EvalTree):
    pass


class IqtreeEvalTree(EvalTree):
    pass


class BaseRFDistance(BaseModel):
    plain_rfdist = P.FloatField()
    normalized_rfdist = P.FloatField()


class RFDistTreesearchTree(BaseRFDistance):
    tree1 = P.ForeignKeyField(TreesearchTree)
    tree2 = P.ForeignKeyField(TreesearchTree)


class RFDistEvalTree(BaseRFDistance):
    tree1 = P.ForeignKeyField(EvalTree)
    tree2 = P.ForeignKeyField(EvalTree)