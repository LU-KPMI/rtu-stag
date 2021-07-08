import os
import numpy as np
import pandas as pd
from skbio.diversity import alpha
import matplotlib.pyplot as plt


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
        return sum(self.abundance.values())


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


def get_all_samples(patients):
    samples = []
    for p in patients:
        if p.sample_before != None:
            samples.append(p.sample_before)
        if p.sample_after != None:
            samples.append(p.sample_after)
    return samples


def vectorize(sample, otus):
    return [sample.abundance[otu] if otu in sample.abundance else 0 for otu in otus]


def print_abundance_table(samples, output_file):
    all_otus = get_all_otus(samples)
    with open(output_file, "w") as f:
        print("sample_id", *all_otus, sep = ",", file = f)
        for sample in samples:
            print(sample.id, *vectorize(sample, all_otus), sep = ",", file = f)


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
        v = vectorize(sample, category_names)
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


if __name__ == "__main__":
    patients = get_patients_from_table("metadata.csv")
    print_abundance_table(get_all_samples(patients), "abundances.csv")
    get_patient_table(patients).to_csv("patients.csv")
    draw_enterotypes(get_all_samples(patients), "sample_enterotypes.svg")
