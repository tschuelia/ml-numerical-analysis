from peewee import *

db = SqliteDatabase(None)


class Run(Model):
    blmin = FloatField()
    blmax = FloatField()
    final_raxml_llh = FloatField()
    best_iqtree_llh = FloatField()
    average_absolute_rf_distance = FloatField()
    average_relative_rf_distance = FloatField()
    unique_topos = IntegerField()

    class Meta:
        database = db


class Raxml_Tree(Model):
    run = ForeignKeyField(Run)
    raxml_tree = CharField()
    is_best = BooleanField()

    class Meta:
        database = db


class IQ_Tree(Model):
    run = ForeignKeyField(Run)
    iq_tree = CharField()  # tree string of iq_tree

    # the following are the results of the iqtree run
    # achieved by comparing with best raxml tree

    llh = FloatField()  # llh as recalculated by iqtree run
    deltaL = (
        FloatField()
    )  # llh difference to max llh in the set according to iqtree run

    bpRell = FloatField(null=True)  # bootstrap proportion using RELL method.
    # True denotes the 95% confidence sets.
    # False denotes significant exclusion.
    bpRell_significant = BooleanField(null=True)

    pKH = FloatField(null=True)  # p-value of one sided Kishino-Hasegawa test (1989).
    pKH_significant = BooleanField(null=True)

    pSH = FloatField(null=True)  # p-value of Shimodaira-Hasegawa test (2000).
    pSH_significant = BooleanField(null=True)

    pWKH = FloatField(null=True)  # p-value of weighted KH test.
    pWKH_significant = BooleanField(null=True)

    pWSH = FloatField(null=True)  # p-value of weighted SH test.
    pWSH_significant = BooleanField(null=True)

    cELW = FloatField(null=True)  # Expected Likelihood Weight
    cELW_significant = BooleanField(null=True)

    pAU = FloatField(null=True)  # p-value of approximately unbiased (AU) test.
    pAU_significant = BooleanField(null=True)

    class Meta:
        database = db


class RFDistance(Model):
    tree1 = ForeignKeyField(Raxml_Tree)
    tree2 = ForeignKeyField(Raxml_Tree)
    plain_rf_distance = FloatField()
    normalized_rf_distance = FloatField()

    class Meta:
        database = db