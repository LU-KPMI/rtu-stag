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

    def read_count(self):
        return sum(self.abundance.values())
    
    def read_count_resistome(self):
        return sum(self.resistome.values())


class Patient:
    def __init__(self, id, sample_id_before, sample_id_after, treatment):
        self.id = id
        self.sample_before = Sample(sample_id_before) if type(sample_id_before) == str else None
        self.sample_after = Sample(sample_id_after) if type(sample_id_after) == str else None
        self.treatment = treatment


def get_patients_from_table(filename):
    patients = []
    table = pd.read_csv(filename, index_col = "patient_id")
    for index, row in table.iterrows():
        patients.append(Patient(index, row["sample_id_before"], row["sample_id_after"], row["treatment"]))
    return patients


def get_all_otus(samples):
    all_otus = set()
    for sample in samples:
        for otu in sample.abundance:
            all_otus.add(otu)
    return sorted(all_otus)


def get_all_mechanisms(samples):
    all_mechanisms = set()
    for sample in samples:
        if sample.resistome != None:
            for otu in sample.resistome:
                all_mechanisms.add(otu)
    return sorted(all_mechanisms)


def get_all_samples(patients):
    samples = []
    for p in patients:
        if p.sample_before != None:
            samples.append(p.sample_before)
        if p.sample_after != None:
            samples.append(p.sample_after)
    return samples


def vectorize_abundance(sample, otus):
    return [sample.abundance[otu] if otu in sample.abundance else 0 for otu in otus]


def vectorize_resistome(sample, mechanisms):
    return [sample.resistome[mechanism] if mechanism in sample.resistome else 0 for mechanism in mechanisms]


def print_abundance_table(samples, output_file):
    all_otus = get_all_otus(samples)
    with open(output_file, "w") as f:
        print("sample_id", *all_otus, sep = ",", file = f)
        for sample in samples:
            print(sample.id, *vectorize_abundance(sample, all_otus), sep = ",", file = f)


def print_resistome_table(samples, output_file):
    all_mechanisms = get_all_mechanisms(samples)
    with open(output_file, "w") as f:
        print("sample_id", *all_mechanisms, sep = ",", file = f)
        for sample in samples:
            if sample.resistome != None:
                print(sample.id, *vectorize_resistome(sample, all_mechanisms), sep = ",", file = f)


def get_significant_otus(samples):
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


def get_significant_mechanisms(samples):
    significant = []

    for sample in samples:
        if sample.resistome == None:
            continue
        total_reads = sample.read_count_resistome()

        while True:
            s = 0
            for i in significant:
                s += sample.resistome[i] if i in sample.resistome else 0
            if s > 0.7 * total_reads:
                break
            i_take = None
            for i in sample.resistome:
                if i in significant:
                    continue
                if i_take == None or sample.resistome[i] > sample.resistome[i_take]:
                    i_take = i
            significant.append(i_take)

    return significant


def calc_enterotype(sample, significant_otus):
    total = sample.read_count()
    max_otu = None
    for otu in significant_otus:
        if otu in sample.abundance:
            if max_otu == None or sample.abundance[otu] > sample.abundance[max_otu]:
                max_otu = otu
    if sample.abundance[max_otu] > total / 5:
        if max_otu == "Bacteroides":
            return "A"
        elif max_otu == "Prevotella":
            return "B"
        else:
            return "C"
    return "D"


def get_patient_table(patients): # TODO: pretify
    metrics = [
        (alpha.chao1, "chao1"),
        (alpha.shannon, "shannon"),
        (alpha.observed_otus, "observed"),
        (alpha.ace, "ace"),
        (alpha.berger_parker_d, "berger_parker"),
    ]

    significant_otus = get_significant_otus(get_all_samples(patients))

    column_names = []
    for _, metric in metrics:
        column_names.append(metric + "_before")
        column_names.append(metric + "_after")
    column_names.append("enterotype_before")
    column_names.append("enterotype_after")
    column_names.append("read_count_before")
    column_names.append("read_count_after")
    column_names.append("treatment")

    all_data = []

    for patient in patients:
        print(patient.id)
        values = []
        for metric, _ in metrics:
            values.append(metric(list(patient.sample_before.abundance.values())) if patient.sample_before != None else None)
            values.append(metric(list(patient.sample_after.abundance.values())) if patient.sample_after != None else None)
        values.append(calc_enterotype(patient.sample_before, significant_otus) if patient.sample_before != None else None)
        values.append(calc_enterotype(patient.sample_after, significant_otus) if patient.sample_after != None else None)
        values.append(patient.sample_before.read_count() if patient.sample_before != None else None)
        values.append(patient.sample_after.read_count() if patient.sample_after != None else None)
        values.append(patient.treatment)
        all_data.append(values)

    df = pd.DataFrame(columns=column_names, index=[p.id for p in patients], data = all_data)
    df.index.name = "patient_id"
    return df


def get_relative_abundances(samples, category_names):
    d = []
    for sample in samples:
        v = vectorize_abundance(sample, category_names)
        v.append(sample.read_count() - sum(v))
        v = np.array(v)
        v = v / sum(v)
        d.append(v)
    return np.array(d)


def draw_enterotypes(samples, filename):
    category_names = get_significant_otus(samples)
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

    plt.savefig(filename)
    plt.close(fig)


def get_relative_abundances2(samples, category_names):
    d = []
    for sample in samples:
        if sample == None or sample.resistome == None:
            v = [0] * (len(category_names) + 1)
        else:
            v = vectorize_resistome(sample, category_names)
            v.append(sample.read_count_resistome() - sum(v))
            v = np.array(v)
            v = v / sum(v)
        d.append(v)
    return np.array(d)


def draw_chart(d, labels, category_names, category_colors, ax):
    d_cum = d.cumsum(axis=1)
    for i, (name, color) in enumerate(zip(category_names, category_colors)):
        widths = d[:, i]
        starts = d_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5, label=name, color=color)


def draw_resistomes(patients, filename):
    fig, axs = plt.subplots(2, 3, figsize=(60, 30))
    fig.suptitle("Resistome mechanisms", fontsize=40)
    axs = axs.transpose()

    significant_mechanisms = get_significant_mechanisms(get_all_samples(patients))
    category_names = significant_mechanisms + ["Other"]
    axs[0][0].set_ylabel("Before", fontsize=40)
    axs[0][1].set_ylabel("After", fontsize=40)

    for treatment, ax_col in zip(["Control", "STD2", "STD3"], axs):
        subdata = [p for p in patients if p.treatment == treatment]

        labels = [p.id for p in subdata]
        category_colors = plt.get_cmap('tab20')(np.linspace(0.05, 0.95, len(category_names)))

        ax_col[0].set_title(treatment, fontsize=40)

        draw_chart(get_relative_abundances2([p.sample_before for p in subdata], significant_mechanisms),
                labels,
                category_names,
                category_colors,
                ax_col[0])
        draw_chart(get_relative_abundances2([p.sample_after for p in subdata], significant_mechanisms),
                labels,
                category_names,
                category_colors,
                ax_col[1])

    axs[0][0].legend(ncol=len(category_names), bbox_to_anchor=(0, 1.2), loc='lower left', fontsize='small')

    plt.savefig(filename)
    plt.close(fig)


def print_helicobacter_abundance(patients):
    fig, axs = plt.subplots(3, figsize=(8, 24))
    for treatment, ax in zip(["Control", "STD2", "STD3"], axs):
        X = []
        Y = []
        print("Treatment:", treatment)
        data = [p for p in patients if p.treatment == treatment]
        for p in data:
            if p.sample_before != None and p.sample_after != None:
                x = p.sample_before.abundance["Helicobacter"] if "Helicobacter" in p.sample_before.abundance else 0
                x /= p.sample_before.read_count()
                X.append(x)
                y = p.sample_after.abundance["Helicobacter"] if "Helicobacter" in p.sample_after.abundance else 0
                y /= p.sample_after.read_count()
                Y.append(y)


            if p.sample_before != None:
                x = p.sample_before.abundance["Helicobacter"] if "Helicobacter" in p.sample_before.abundance else 0
                x /= p.sample_before.read_count()
                x = "{:.2e}".format(x)
            else:
                x = "N/A"

            if p.sample_after != None:
                y = p.sample_after.abundance["Helicobacter"] if "Helicobacter" in p.sample_after.abundance else 0
                y /= p.sample_after.read_count()
                y = "{:.2e}".format(y)
            else:
                y = "N/A"
            print(x, y)
        ax.scatter(X, Y)
        ax.set_title(treatment)
        ax.set_xlabel("abundance before")
        ax.set_ylabel("abundance after")
        ax.set_xlim([0, 0.0005])
        ax.set_ylim([0, 0.0005])
    plt.savefig("helicobacter.svg")

if __name__ == "__main__":
    patients = get_patients_from_table("metadata.csv")
#    print_abundance_table(get_all_samples(patients), "abundances.csv")
#    print_resistome_table(get_all_samples(patients), "resistomes.csv")
#    get_patient_table(patients).to_csv("patients.csv")
#    draw_enterotypes(get_all_samples(patients), "sample_enterotypes.svg")
#    draw_resistomes(patients, "sample_resistomes.svg")
    print_helicobacter_abundance(patients)

