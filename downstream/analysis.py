import functools
import math
import numpy as np
import os
import operator
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from skbio.diversity import alpha
from scipy.spatial import distance
from scipy.cluster.hierarchy import dendrogram, linkage


def get_non_none_attributes(v, attr):
    f = operator.attrgetter(attr)
    return [f(x) for x in v if f(x) != None]


def parse_bracken(filename):
    with open(filename) as f:
        abundance = {}
        first = True
        for line in f:
            if first:
                first = False
                continue
            data = line.rstrip().split('\t')
            genus, cnt = data[0], int(data[5])
            abundance[genus] = cnt
        return abundance


class Sample:
    def __init__(self, id):
        self.id = id
        self.abundance = parse_bracken(os.path.join("bracken_output", id + ".bracken"))

    def read_count(self):
        return sum(self.abundance.values()) # TODO: Speedup and fix


class Patient:
    def __init__(self, id, sample_id_before, sample_id_after, treatment):
        self.id = id
        self.sample_before = Sample(sample_id_before) if type(sample_id_before) == str else None
        self.sample_after = Sample(sample_id_after) if type(sample_id_after) == str else None
        self.treatment = treatment


def get_all_genuses(samples):
    all_genuses = set()
    for sample in samples:
        for g in sample.abundance:
            all_genuses.add(g)
    return sorted(all_genuses)


def vectorize(sample, otus):
    return [sample.abundance[otu] if otu in sample.abundance else 0 for otu in otus]


def print_abundance_table(samples, output_file):
    all_genuses = get_all_genuses(samples)
    with open(output_file, "w") as f:
        print("sample_id", *all_genuses, sep = ",", file = f)
        for sample in samples:
            print(sample.id, *vectorize(sample, all_genuses), sep = ",", file = f)


def all_samples(patients):
    samples = []
    for patient in patients:
        if patient.sample_before != None:
            samples.append(patient.sample_before)
        if patient.sample_after != None:
            samples.append(patient.sample_after)
    return samples


def get_patients_from_table(filename):
    patients = []
    table = pd.read_csv(filename, index_col = "patient_id")
    for index, row in table.iterrows():
        patients.append(Patient(index, row["sample_id_before"], row["sample_id_after"], row["treatment"]))
    return patients


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


def get_significant_categories(samples):
    significant = []

    for sample in samples:
        total_reads = sample.read_count()

        while True:
            s = 0
            for i in significant:
                s += sample.abundance[i] if i in sample.abundance else 0
            if s > 0.3 * total_reads:
                break
            i_take = None
            for i in sample.abundance:
                if i in significant:
                    continue
                if i_take == None or sample.abundance[i] > sample.abundance[i_take]:
                    i_take = i
            significant.append(i_take)

#    return ['Bacteroides', 'Prevotella']
    return significant


def get_relative_abundances(samples, category_names):
    d = []
    for sample in samples:
        v = vectorize(sample, category_names)
        v.append(sample.read_count() - sum(v))
        v = np.array(v)
        v = v / sum(v)
        d.append(v)
    return np.array(d)


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
    (alpha.fisher_alpha, "fisher alpha"),
    (alpha.gini_index, "gini"),
    (alpha.goods_coverage, "goods coverage"),
    (alpha.heip_e, "heip evenness"),
    (alpha.margalef, "margalef"),
    (alpha.mcintosh_d, "mcintosh_d"),
    (alpha.mcintosh_e, "mcintosh_e"),
    (alpha.menhinick, "menhinick"),
    (alpha.pielou_e, "pielou_e"),
    (alpha.robbins, "robbins"),
    (alpha.simpson, "simpson"),
    (alpha.simpson_e, "simpson_e"),
    (alpha.singles, "singles"),
    (alpha.strong, "strong"),
]


def order_by_treatment(patients):
    P = []
    for treatment in ["STD2", "STD3", "Control"]:
        for p in patients:
            if p.treatment == treatment:
                P.append(p)
    return P


def get_dist_matrix(m, f):
    n = len(m)
    d = np.zeros(shape=(n, n))
    for i in range(n):
        for j in range(n):
            d[i][j] = f(m[i], m[j])
    return d


def draw_violin_plots(patients_both, metric, metric_name, pdf):
    fig, axs = plt.subplots(3, 2, figsize=(10, 20))

    max_y = max([metric(list(s.abundance.values())) for s in all_samples(patients_both)])
    if max_y < 10:
        max_y = math.ceil(max_y)
    else:
        max_y = math.ceil(max_y / 100) * 100

    fig.suptitle(metric_name)

    axs[0][0].set_title("Before")
    axs[0][1].set_title("After")

    for i, treatment in enumerate(["Control", "STD2", "STD3"]):
        axs[i][0].violinplot([metric(list(p.sample_before.abundance.values())) for p in patients_both if p.treatment == treatment])
        axs[i][1].violinplot([metric(list(p.sample_after.abundance.values())) for p in patients_both if p.treatment == treatment])

        for ax in axs[i]:
            ax.set_ylim(0, max_y)
        axs[i][0].set_ylabel(treatment)

    pdf.savefig(fig)
    plt.close(fig)


def draw_rarefaction(samples, f, pdf):
    fig, ax = plt.subplots()
    ax.grid()
    ax.ticklabel_format(axis='x', style='plain')
    ax.xaxis.set_tick_params(rotation=90)
    ax.set_xlabel("Read count")
    ax.set_ylabel("Observed species")

    for s in samples:
        x, y = rarefy(list(s.abundance.values()), f)
        ax.plot(x, y, linewidth=0.12)
        ax.text(x[-1], y[-1] - 2, s.id, fontsize=2)

    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


def draw_heatmap(matrix, title, ax):
    g = sns.heatmap(matrix, linewidth = 0.5, ax=ax, square=True)
    g.set(xlabel = None, ylabel = None)
    g.set(xticklabels = [], yticklabels = [])
    g.tick_params(left=False, bottom=False)
    g.set(title=title)
    return g


def draw_before_after_heatmaps(d_before, d_after, pdf):
    fig, axs = plt.subplots(2, 1, figsize=(7, 14))
    draw_heatmap(d_before, "Before", axs[0])
    draw_heatmap(d_after, "After", axs[1])
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def draw_before_after_dendrograms(d_before, d_after, labels, pdf):
    fig, axs = plt.subplots(2, 1, figsize=(11, 6))
    linkage_matrix_before = linkage(distance.squareform(d_before), "single")
    linkage_matrix_after = linkage(distance.squareform(d_after), "single")
    dendrogram(linkage_matrix_before, labels=labels, ax=axs[0])
    dendrogram(linkage_matrix_after, labels=labels, ax=axs[1])
    axs[0].set_title("before")
    axs[1].set_title("after")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def draw_enterotypes(samples, pdf):
    category_names = get_significant_categories(samples)
    d = get_relative_abundances(samples, category_names)
    d_cum = d.cumsum(axis=1)

    category_names.append("Other")
    labels = [s.id for s in samples]
    category_colors = plt.get_cmap('tab20')(np.linspace(0.05, 0.95, len(category_names)))

    fig, ax = plt.subplots(figsize=(70, 100))
    for i, (name, color) in enumerate(zip(category_names, category_colors)):
        widths = d[:, i]
        starts = d_cum[:, i] - widths
        rects = ax.barh(labels, widths, left=starts, height=0.5, label=name, color=color)
        r, g, b, _ = color
        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
#       ax.bar_label(rects, label_type='center', color=text_color)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')

    pdf.savefig(fig)
    plt.close(fig)


def draw_scatter_plot(title, v, f_x, f_y, xlabel, ylabel, pdf):
    plt.scatter([f_x(s) for s in v], [f_y(s) for s in v], 1)
    plt.title(title)
    plt.xlabel("read count")
    plt.ylabel("observed")
    plt.xlim(0, 1.2e8)
    plt.ylim(0, 1050)
    pdf.savefig()
    plt.close()


if __name__ == "__main__":
    patients = get_patients_from_table("metadata.csv")
    patients_both = order_by_treatment([p for p in patients if p.sample_before != None and p.sample_after != None])

    all_genuses = get_all_genuses(all_samples(patients))

    d_before = get_dist_matrix([vectorize(p.sample_before, all_genuses) for p in patients_both], distance.jensenshannon)
    d_after = get_dist_matrix([vectorize(p.sample_after, all_genuses) for p in patients_both], distance.jensenshannon)

    print_abundance_table(all_samples(patients), "abundances.csv")

    with PdfPages("graphs.pdf") as pdf:
        draw = functools.partial(draw_scatter_plot,
                f_x = lambda t : t.read_count(),
                f_y = lambda t : alpha.observed_otus(list(t.abundance.values())),
                xlabel = "read count",
                ylabel = "observed",
                pdf = pdf)

        draw("Visi", all_samples(patients))
        for treatment in ["Control", "STD2", "STD3"]:
            get = functools.partial(get_non_none_attributes,
                    v = [p for p in patients if p.treatment == treatment])
            draw(treatment + " & before", get(attr = 'sample_before'))
            draw(treatment + " & after",  get(attr = 'sample_after'))

        for metric, metric_name in metrics:
            draw_violin_plots(patients_both, metric, metric_name, pdf)

        draw_rarefaction(all_samples(patients), alpha.observed_otus, pdf)
        draw_before_after_heatmaps(d_before, d_after, pdf)
        draw_before_after_dendrograms(d_before, d_after, [p.id for p in patients_both], pdf)
        draw_enterotypes(all_samples(patients), pdf)
