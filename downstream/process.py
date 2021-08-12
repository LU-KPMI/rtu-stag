import os
import numpy as np
import pandas as pd
from skbio.diversity import alpha
import matplotlib.pyplot as plt


def parse_bracken(filename):
    if not os.path.isfile(filename):
        return None
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


def parse_amrplusplus(filename):
    if not os.path.isfile(filename):
        return None
    with open(filename) as f:
        resistome = {}
        first = True
        for line in f:
            if first:
                first = False
                continue
            data = line.rstrip().split('\t')
            genus, cnt = data[1], int(data[2])
            resistome[genus] = cnt
        return resistome


class Sample:
    def __init__(self, id):
        self.id = id
        self.abundance = parse_bracken(os.path.join("bracken_output", id + ".bracken"))
        self.resistome = parse_amrplusplus(os.path.join("amrplusplus_report", id + ".tsv"))

    def vectorize(self, keys, attr):
        m = getattr(self, attr)
        return [m[x] if x in m else 0 for x in keys]

    def calc_enterotype(self, significant_otus):
        if self.abundance == None:
            return None
        total = sum(self.abundance.values())
        max_otu = None
        for otu in significant_otus:
            if otu in self.abundance:
                if max_otu == None or self.abundance[otu] > self.abundance[max_otu]:
                    max_otu = otu
        if self.abundance[max_otu] > total / 5:
            if max_otu == "Bacteroides":
                return "A"
            elif max_otu == "Prevotella":
                return "B"
            else:
                return "C"
        return "D"

    def apply_diversity_metric(self, metric):
        if self.abundance == None:
            return None
        return metric(list(self.abundance.values()))


class Patient:
    def __init__(self, id, sample_id_before, sample_id_after, treatment):
        self.id = id
        self.sample_before = Sample(sample_id_before) if type(sample_id_before) == str else None
        self.sample_after = Sample(sample_id_after) if type(sample_id_after) == str else None
        self.treatment = treatment


def read_patients_from_csv(filename):
    patients = []
    table = pd.read_csv(filename, index_col = "patient_id")
    for index, row in table.iterrows():
        patients.append(Patient(index, row["sample_id_before"], row["sample_id_after"], row["treatment"]))
    return patients


def get_all_samples(patients):
    samples = []
    for p in patients:
        if p.sample_before != None:
            samples.append(p.sample_before)
        if p.sample_after != None:
            samples.append(p.sample_after)
    return samples


def get_all_sample_keys(samples, attr):
    all_keys = set()
    for sample in samples:
        if getattr(sample, attr) != None:
            for k in getattr(sample, attr):
                all_keys.add(k)
    return sorted(all_keys)


def print_attr_table(samples, output_file, attr):
    keys = get_all_sample_keys(samples, attr)
    with open(output_file, "w") as f:
        print("sample_id", *keys, sep = ",", file = f)
        for sample in samples:
            if getattr(sample, attr) != None:
                print(sample.id, *sample.vectorize(keys, attr), sep = ",", file = f)


def print_abundance_table(samples, output_file):
    print_attr_table(samples, output_file, "abundance")


def print_resistome_table(samples, output_file):
    print_attr_table(samples, output_file, "resistome")


def get_significant_keys(samples, attr, offset):
    significant = []

    for sample in samples:
        m = getattr(sample, attr)
        if m == None:
            continue
        total_reads = sum(m.values())

        while True:
            s = sum(sample.vectorize(significant, attr))
            if s > offset * total_reads:
                break
            i_take = None
            for i in m:
                if i in significant:
                    continue
                if i_take == None or m[i] > m[i_take]:
                    i_take = i
            significant.append(i_take)

#    return ['Bacteroides', 'Prevotella']
    return significant


def get_significant_otus(samples):
    return get_significant_keys(samples, "abundance", 0.3)


def get_significant_resistance_mechanisms(samples):
    return get_significant_keys(samples, "resistome", 0.7)


def get_patient_table(patients):
    def skip_if_none(f, x):
        return f(x) if x != None else None
    def get_metric_applier(metric):
        return lambda sample: sample.apply_diversity_metric(metric)
    def append_both(name, func):
        columns.append((name + "_before", lambda patient: skip_if_none(func, patient.sample_before)))
        columns.append((name + "_after", lambda patient: skip_if_none(func, patient.sample_after)))

    metrics = [
        (alpha.chao1, "chao1"),
        (alpha.shannon, "shannon"),
        (alpha.observed_otus, "observed"),
        (alpha.ace, "ace"),
        (alpha.berger_parker_d, "berger_parker"),
    ]

    significant_otus = get_significant_otus(get_all_samples(patients))

    columns = []

    for metric, metric_name in metrics:
        append_both(metric_name, get_metric_applier(metric))
    append_both("enterotype", lambda sample: sample.calc_enterotype(significant_otus))
    append_both("read_count", lambda sample: sum(sample.abundance.values()))
    columns.append(("treatment", lambda patient: patient.treatment))

    all_data = []

    for patient in patients:
        print(patient.id)
        all_data.append([f(patient) for _, f in columns])

    df = pd.DataFrame(columns=[name for name, _ in columns], index=[p.id for p in patients], data = all_data)
    df.index.name = "patient_id"
    return df


def get_relative_abundances(samples, category_names, attr):
    d = []
    for sample in samples:
        if sample == None or getattr(sample, attr) == None:
            v = [0] * (len(category_names) + 1)
        else:
            m = getattr(sample, attr)
            v = sample.vectorize(category_names, attr)
            v.append(sum(m.values()) - sum(v))
            v = np.array(v)
            v = v / sum(v)
        d.append(v)
    return np.array(d)


def draw_evenness(patients, suptitle, significant_categories_getter, attr, filename):
    fig, axs = plt.subplots(2, 3, figsize=(60, 30))
    fig.suptitle(suptitle, fontsize=40)
    axs = axs.transpose()

    significant_categories = significant_categories_getter(get_all_samples(patients))
    category_names = significant_categories + ["Other"]
    category_colors = plt.get_cmap('tab20')(np.linspace(0.05, 0.95, len(category_names)))

    def draw_chart(d, labels, ax):
        d_cum = d.cumsum(axis=1)
        for i, (name, color) in enumerate(zip(category_names, category_colors)):
            widths = d[:, i]
            starts = d_cum[:, i] - widths
            ax.barh(labels, widths, left=starts, height=0.5, label=name, color=color)

    axs[0][0].set_ylabel("Before", fontsize=40)
    axs[0][1].set_ylabel("After", fontsize=40)

    for treatment, ax_col in zip(["Control", "STD2", "STD3"], axs):
        subdata = [p for p in patients if p.treatment == treatment]

        labels = [p.id for p in subdata]

        ax_col[0].set_title(treatment, fontsize=40)

        draw_chart(get_relative_abundances([p.sample_before for p in subdata], significant_categories, attr),
                labels,
                ax_col[0])
        draw_chart(get_relative_abundances([p.sample_after for p in subdata], significant_categories, attr),
                labels,
                ax_col[1])

    axs[0][0].legend(ncol=len(category_names), bbox_to_anchor=(0, 1.2), loc='lower left', fontsize='small')

    plt.savefig(filename)
    plt.close(fig)


def draw_enterotypes(patients, filename):
    draw_evenness(patients, "Relative abundance", get_significant_otus, "abundance", filename)


def draw_resistomes(patients, filename):
    draw_evenness(patients, "Resistome mechanisms", get_significant_resistance_mechanisms, "resistome", filename)


def draw_helicobacter_abundance(patients, filename):
    fig, axs = plt.subplots(3, figsize=(8, 24))
    for treatment, ax in zip(["Control", "STD2", "STD3"], axs):
        X = []
        Y = []
        data = [p for p in patients if p.treatment == treatment]
        for p in data:
            if p.sample_before != None and p.sample_after != None:
                x = p.sample_before.abundance["Helicobacter"] if "Helicobacter" in p.sample_before.abundance else 0
                x /= sum(p.sample_before.abundance.values())
                X.append(x)
                y = p.sample_after.abundance["Helicobacter"] if "Helicobacter" in p.sample_after.abundance else 0
                y /= sum(p.sample_after.abundance.values())
                Y.append(y)
        ax.scatter(X, Y)
        ax.set_title(treatment)
        ax.set_xlabel("abundance before")
        ax.set_ylabel("abundance after")
        ax.set_xlim([0, 0.0005])
        ax.set_ylim([0, 0.0005])
    plt.savefig(filename)


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


def draw_rarefaction(samples, filename):
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.grid()
    ax.ticklabel_format(axis='x', style='plain')
    ax.xaxis.set_tick_params(rotation=90)
    ax.set_xlabel("Read count")
    ax.set_ylabel("Observed species")

    for s in samples:
        x, y = rarefy(list(s.abundance.values()), alpha.observed_otus)
        ax.plot(x, y, linewidth=0.12)
        ax.text(x[-1], y[-1] - 2, s.id, fontsize=2)

    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    patients = read_patients_from_csv("metadata.csv")

    print_abundance_table(get_all_samples(patients), "abundances.csv")
    print_resistome_table(get_all_samples(patients), "resistomes.csv")
    get_patient_table(patients).to_csv("patients.csv")
    draw_enterotypes(patients, "enterotypes.svg")
    draw_resistomes(patients, "resistomes.svg")
    draw_helicobacter_abundance(patients, "helicobacter.svg")
    draw_rarefaction(get_all_samples(patients), "rarefaction.svg")
