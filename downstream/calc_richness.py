#!/bin/python3

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.diversity import alpha


def rarefy(v, f):
    X = []
    Y = []
    step = 5
    sizes = [sum(v) * f // 100 for f in range(step, 101, step)]
    arr = []
    for i, x in enumerate(v):
        arr += [i] * x
    np.random.shuffle(arr)

    cur = [0] * len(v)

    j = 0
    for i, x in enumerate(arr):
        cur[x] += 1
        if i + 1 == sizes[j]:
            X.append(sum(cur))
            Y.append(f(cur))
            j += 1
    return X, Y


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

fig, ax = plt.subplots()
ax.grid()
ax.ticklabel_format(axis='x', style='plain')
ax.xaxis.set_tick_params(rotation=90)
ax.set_xlabel("Read count")
ax.set_ylabel("Observed species")

table = pd.read_csv("abundances.csv", index_col=0)
for index, row in table.iterrows():
    v = list(row)
    print(index, *map(lambda f: f(v), list(zip(*metrics))[0]), sum(v), sep=",")


    x, y = rarefy(v, alpha.observed_otus)
    ax.plot(x, y, linewidth=0.12)
    ax.text(x[-1], y[-1] - 2, index, fontsize=2)

ax.set_xlim(0, None)
ax.set_ylim(0, None)

plt.savefig("rarefaction.svg", bbox_inches='tight')
