#!/bin/python3

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def shannon(v):
    s = sum(v)
    return sum(-(x/s) * math.log(x/s) if x > 0 else 0 for x in v)


def observed(v):
    return sum(1 if x > 0 else 0 for x in v)


def chao1(v):
    obs = observed(v)
    n1 = sum(1 if x == 1 else 0 for x in v)
    n2 = sum(1 if x == 2 else 0 for x in v)

    return obs + n1 * (n1 - 1) / (2 * (n2 + 1))


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


print("sample_name", "chao1", "shannon", "observed", sep=",")

fig, ax = plt.subplots()
ax.grid()
ax.ticklabel_format(axis='x', style='plain')
ax.xaxis.set_tick_params(rotation=90)
ax.set_xlabel("Read count")
ax.set_ylabel("Observed species")

table = pd.read_csv("abundances.csv", index_col=0)
for index, row in table.iterrows():
    v = list(row)
    print(index, chao1(v), shannon(v), observed(v), sep=",")

    x, y = rarefy(v, observed)
    ax.plot(x, y, linewidth=0.12)
    ax.text(x[-1], y[-1] - 2, index, fontsize=2)

ax.set_xlim(0, None)
ax.set_ylim(0, None)

plt.savefig("rarefaction.svg", bbox_inches='tight')
