import pandas as pd
import numpy as np
from skbio.diversity import alpha


metrics = [
    (alpha.chao1, "chao1"),
    (alpha.shannon, "shannon"),
    (alpha.observed_otus, "observed"),
    (alpha.ace, "ace"),
    (alpha.berger_parker_d, "berger-parker"),
    (alpha.brillouin_d, "brillouin"),
    (alpha.dominance, "dominance"),
    (alpha.doubles, "doubles"),
    (alpha.enspie, "enspie"),
#    (alpha.esty_ci, "esty_ci"), # Disabled because of runtime error
    (alpha.fisher_alpha, "fisher alpha"),
    (alpha.gini_index, "gini"),
    (alpha.goods_coverage, "goods coverage"),
    (alpha.heip_e, "heip evenness"),
#    (alpha.kempton_taylor_q, "kempton_taylor_q"), # Disabled because of runtime error
#    (alpha.lladser_pe, "lladser"), # Disabled because of runtime error
    (alpha.margalef, "margalef"),
    (alpha.mcintosh_d, "mcintosh_d"),
    (alpha.mcintosh_e, "mcintosh_e"),
    (alpha.menhinick, "menhinick"),
#    (alpha.michaelis_menten_fit, "michaelis_menten_fit"), # Disabled because of runtime error
    (alpha.pielou_e, "pielou_e"),
    (alpha.robbins, "robbins"),
    (alpha.simpson, "simpson"),
    (alpha.simpson_e, "simpson_e"),
    (alpha.singles, "singles"),
    (alpha.strong, "strong"),
    ]


print("sample_name", *(list(zip(*metrics))[1]), "read_count", sep=",")

table = pd.read_csv("abundances.csv", index_col=0)
for index, row in table.iterrows():
    v = list(row)
    print(index, *map(lambda f: f(v), list(zip(*metrics))[0]), sum(v), sep=",")
