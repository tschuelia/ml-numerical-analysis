import peewee as P

db = P.SqliteDatabase(None)


class Run(P.Model):
    blmin = P.FloatField()
    blmax = P.FloatField()
    average_absolute_rf_distance = P.FloatField()
    average_relative_rf_distance = P.FloatField()
    num_unique_topos = P.IntegerField()

    class Meta:
        database = db


class Tree(P.Model):
    run = P.ForeignKeyField(Run)
    raxml_tree = P.CharField()
    iqtree_tree = P.CharField()
    raxml_llh = P.FloatField()
    iqtree_llh = P.FloatField()
    is_best = P.BooleanField()

    # the following are the results of the iqtree run
    # achieved by comparing with best raxml tree
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

    class Meta:
        database = db


class RFDistance(P.Model):
    tree1 = P.ForeignKeyField(Tree)
    tree2 = P.ForeignKeyField(Tree)
    plain_rf_distance = P.FloatField()
    normalized_rf_distance = P.FloatField()

    class Meta:
        database = db
