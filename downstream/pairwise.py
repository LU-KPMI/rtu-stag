import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance


def bray_curtis(u, v):
    assert(len(u) == len(v))
    c = sum(min(x, y) for x, y in zip(u, v))
    return 1 - 2 * c / (sum(u) + sum(v))

genus_to_id = {}
genus_count = 0


def read_abundance_vector(filename):
    v = [0] * genus_count
    with open(filename) as f:
        first = True
        for line in f:
            if not first:
                data = line.rstrip().split('\t')
                genus = data[0]
                cnt = float(data[5])
                v[genus_to_id[genus]] = cnt
            first = False
    return v


for filename in os.listdir("./bracken_1"):
    if not filename.endswith(".bracken"):
        continue
    with open(os.path.join("bracken_1", filename)) as f:
        first = True
        for line in f:
            if not first:
                genus = line.rstrip().split('\t')[0]
                if not genus in genus_to_id:
                    genus_to_id[genus] = genus_count
                    genus_count += 1
            first = False
    

sample_abundances = []

for t in ["STD2", "STD3", "Control"]:
    with open("patient_richness.csv") as f:
        first = True
        for line in f:
            if not first:
                data = line.rstrip().split(',')
                if data[1] != t:
                    continue
                sample_name_before = data[0] + "_" + data[1] + "_" + "T1"
                sample_name_after = data[0] + "_" + data[1] + "_" + "T2"

                before = read_abundance_vector(os.path.join("bracken_1", sample_name_before + ".bracken"))
                after = read_abundance_vector(os.path.join("bracken_1", sample_name_after + ".bracken"))

                sample_abundances.append([data[0], data[1], before, after])
            first = False

n = len(sample_abundances)

d_before, d_after = np.zeros(shape=(n, n)), np.zeros(shape=(n, n))


for i in range(n):
    for j in range(n):
        d_before[i][j] = distance.jensenshannon(sample_abundances[i][2], sample_abundances[j][2])
        d_after[i][j]  = distance.jensenshannon(sample_abundances[i][3], sample_abundances[j][3])

avg_d_before = {}
avg_d_after = {}

for i in range(n):
    for j in range(i):
        if sample_abundances[i][1] == sample_abundances[j][1]:
            t = sample_abundances[i][1]
            if t not in avg_d_before:
                avg_d_before[t] = 0
                avg_d_after[t] = 0
            avg_d_before[t] += d_before[i][j]
            avg_d_after[t] += d_after[i][j]

for t in avg_d_before:
    c = 0
    for s in sample_abundances:
        if s[1] == t:
            c += 1
    avg_d_before[t] /= c * (c - 1) / 2
    avg_d_after[t] /= c * (c - 1) / 2


print('avg_d_before =', avg_d_before)
print('avg_d_after =', avg_d_after)

fig, ax = plt.subplots(2, 1, figsize=(7, 14))

g1 = sns.heatmap(d_before, linewidth = 0.5, ax=ax[0], square=True)
g1.set(xlabel = None, ylabel = None)
g1.set(xticklabels = [], yticklabels = [])
g1.tick_params(left=False, bottom=False)
g1.set(title="Before")

g2 = sns.heatmap(d_after, linewidth = 0.5, ax=ax[1], square=True)
g2.set(xlabel = None, ylabel = None)
g2.set(xticklabels = [], yticklabels = [])
g2.tick_params(left=False, bottom=False)
g2.set(title="After")

plt.tight_layout()
plt.savefig('heatmap.png')
