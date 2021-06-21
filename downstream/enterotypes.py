import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_significant_categories(data):
    significant = []

    for index, row in data.iterrows():
        v = np.array(row)
        v = v / sum(v)

        while True:
            s = 0
            for i in significant:
                s += v[i]
            if s > 0.3:
                break
            i_take = -1
            for i in range(len(row)):
                if i in significant:
                    continue
                if i_take == -1 or v[i] > v[i_take]:
                    i_take = i
            significant.append(i_take)

    category_names = [data.columns[i] for i in significant]
#    return ['Bacteroides', 'Prevotella']
    return category_names


def get_relative_abundances(data, category_names):
    d = []
    for index, row in data.iterrows():
        v = [row[name] for name in category_names]
        v.append(sum(np.array(row)) - sum(v))
        v = np.array(v)
        v = v / sum(np.array(row))

        d.append(v)
    return np.array(d)

data = pd.read_csv("./abundances.csv", index_col=0)
category_names = get_significant_categories(data)
d = get_relative_abundances(data, category_names)
d_cum = d.cumsum(axis=1)

category_names.append("Other")
labels = list(data.index)
category_colors = plt.get_cmap('tab20')(np.linspace(0.05, 0.95, len(category_names)))

fig, ax = plt.subplots(figsize=(70, 100))
for i, (name, color) in enumerate(zip(category_names, category_colors)):
    widths = d[:, i]
    starts = d_cum[:, i] - widths
    rects = ax.barh(labels, widths, left=starts, height=0.5, label=name, color=color)
    r, g, b, _ = color
    text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
#    ax.bar_label(rects, label_type='center', color=text_color)
ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')

plt.savefig("enterotypes.svg")
