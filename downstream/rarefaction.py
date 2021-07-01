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


fig, ax = plt.subplots()
ax.grid()
ax.ticklabel_format(axis='x', style='plain')
ax.xaxis.set_tick_params(rotation=90)
ax.set_xlabel("Read count")
ax.set_ylabel("Observed species")


table = pd.read_csv("abundances.csv", index_col=0)
for index, row in table.iterrows():
    x, y = rarefy(list(row), alpha.observed_otus)
    ax.plot(x, y, linewidth=0.12)
    ax.text(x[-1], y[-1] - 2, index, fontsize=2)

ax.set_xlim(0, None)
ax.set_ylim(0, None)

plt.savefig("rarefaction.svg", bbox_inches='tight')
