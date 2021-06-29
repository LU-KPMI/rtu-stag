import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster.hierarchy import dendrogram, linkage


def bray_curtis(u, v):
    assert(len(u) == len(v))
    c = sum(min(x, y) for x, y in zip(u, v))
    return 1 - 2 * c / (sum(u) + sum(v))

data = pd.read_csv("./abundances.csv", index_col=0)
metadata = pd.read_csv("./sample_metadata.csv", index_col=0)

keep = set(metadata[metadata["time"] == "T1"]["patient_id"]).intersection(set(metadata[metadata["time"] == "T2"]["patient_id"]))
metadata = metadata[metadata["patient_id"].isin(keep)]
data = data[data.index.isin(metadata.index)]

sample_abundances = []

for treatment in ["STD2", "STD3", "Control"]:
    for patient_id in set(metadata[metadata["treatment"] == treatment]["patient_id"]):
        sample_name_before = patient_id + "_" + treatment + "_" + "T1"
        sample_name_after = patient_id + "_" + treatment + "_" + "T2"
        before = np.array(data.loc[sample_name_before])
        after = np.array(data.loc[sample_name_after])
        sample_abundances.append([patient_id, treatment, before, after])



n = len(sample_abundances)

d_before, d_after = np.zeros(shape=(n, n)), np.zeros(shape=(n, n))

for i in range(n):
    for j in range(n):
        d_before[i][j] = distance.jensenshannon(sample_abundances[i][2], sample_abundances[j][2])
        d_after[i][j]  = distance.jensenshannon(sample_abundances[i][3], sample_abundances[j][3])

#avg_d_before = {}
#avg_d_after = {}
#
#for i in range(n):
#    for j in range(i):
#        if sample_abundances[i][1] == sample_abundances[j][1]:
#            t = sample_abundances[i][1]
#            if t not in avg_d_before:
#                avg_d_before[t] = 0
#                avg_d_after[t] = 0
#            avg_d_before[t] += d_before[i][j]
#            avg_d_after[t] += d_after[i][j]
#
#for t in avg_d_before:
#    c = 0
#    for s in sample_abundances:
#        if s[1] == t:
#            c += 1
#    avg_d_before[t] /= c * (c - 1) / 2
#    avg_d_after[t] /= c * (c - 1) / 2
#
#
#print('avg_d_before =', avg_d_before)
#print('avg_d_after =', avg_d_after)

fig, axs = plt.subplots(2, 1, figsize=(7, 14))

g1 = sns.heatmap(d_before, linewidth = 0.5, ax=axs[0], square=True)
g1.set(xlabel = None, ylabel = None)
g1.set(xticklabels = [], yticklabels = [])
g1.tick_params(left=False, bottom=False)
g1.set(title="Before")

g2 = sns.heatmap(d_after, linewidth = 0.5, ax=axs[1], square=True)
g2.set(xlabel = None, ylabel = None)
g2.set(xticklabels = [], yticklabels = [])
g2.tick_params(left=False, bottom=False)
g2.set(title="After")

plt.tight_layout()
plt.savefig('heatmap.png')


fig, axs = plt.subplots(2, 1, figsize=(11, 6))
linkage_matrix_before = linkage(distance.squareform(d_before), "single")
linkage_matrix_after = linkage(distance.squareform(d_after), "single")
dendrogram(linkage_matrix_before, labels=list(zip(*sample_abundances))[1], ax=axs[0])
dendrogram(linkage_matrix_after, labels=list(zip(*sample_abundances))[1], ax=axs[1])
axs[0].set_title("before")
axs[1].set_title("after")
plt.tight_layout()
plt.savefig("dendrograms.svg")

