import itertools
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def subsetting(data):
    numerical_before = ["shannon_before", "observed_before", "read_count_before"]
    categorical_before = ["enterotype_before", "treatment"]
    numerical_after = ["shannon_after", "observed_after", "read_count_after"]
    categorical_after = ["enterotype_after"]

    categories = [sorted(data[c].dropna().unique()) for c in categorical_before]
    with PdfPages("outputs/pictures/subsetting.pdf") as pdf:
        for target in numerical_after + categorical_after:
            fig, axs = plt.subplots(10, 10, figsize=(50, 50))
            fig.suptitle(target, size=60)
            axs = [x for ax in axs for x in ax]
            i = 0
            for vals in itertools.product(*categories):
                df = data
                for c, v in zip(categorical_before, vals):
                    df = df[df[c] == v]
                for n in numerical_before:
                    df2 = df[[n, target]].dropna()
                    axs[i].scatter(x=df2[n], y=df2[target])
                    axs[i].set_title(" + ".join(vals) + " : " + n)
                    i += 1
            pdf.savefig(fig)
            plt.close(fig)


numerical = ["chao1_before", "chao1_after", "shannon_before", "shannon_after", "observed_before", "observed_after", "ace_before", "ace_after", "berger_parker_before", "berger_parker_after", "read_count_before", "read_count_after"]
categorical = ["enterotype_before", "enterotype_after", "treatment"]


def draw_pairplots(data):
    g = sns.pairplot(data, diag_kind="kde", palette=["#F4F1DE", "#AA3F22"])
    g.map_lower(sns.kdeplot, levels=4, color=".4")
    plt.savefig("outputs/pictures/pairplots.svg")
    plt.close()


def draw_categorical_count(data):
    fig, ax = plt.subplots(1, 3, figsize=(20, 10))
    for variable, subplot in zip(categorical, ax.flatten()):
        sns.countplot(data[variable], ax=subplot)
        for label in subplot.get_xticklabels():
            label.set_rotation(90)
    plt.savefig("outputs/pictures/categorical_count.svg")
    plt.close()


def draw_violin_plots(data):
    metrics = ["chao1", "shannon", "observed", "ace", "berger_parker"]
    fig, axs = plt.subplots(nrows=1, ncols=len(metrics), figsize=(25, 5))
    fig.suptitle("Alpha diversity")
    for ax, metric in zip(axs, metrics):
        data_sub = data[["treatment", metric + "_before", metric + "_after"]].dropna()

        max_y = max(max(data_sub[metric + "_before"]), max(data_sub[metric + "_after"]))
        if max_y < 10:
            max_y = math.ceil(max_y)
        else:
            max_y = math.ceil(max_y / 100) * 100

        d = []
        labels = []
        for time in [metric + "_before", metric + "_after"]:
            for treatment in ["STD2", "STD3", "Control"]:
                d.append(list(data_sub[data_sub["treatment"] == treatment][time]))
                labels.append(time + " " + treatment)

        parts = ax.violinplot(d)
        ax.set_xticks(range(1, 7))
        ax.set_xticklabels(labels, rotation=90)
        ax.set_title(metric)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(['red', 'green', 'blue'][i % 3])

    plt.savefig("outputs/pictures/violin_plots.svg", bbox_inches="tight")
    plt.close(fig)


def draw_numerical_vs_categorical(data):
    fig, axs = plt.subplots(len(numerical), len(categorical), figsize=(50, 100))
    for i, numvar in enumerate(numerical):
        for j, catvar in enumerate(categorical):
            sns.boxplot(x=catvar, y=numvar, data=data, ax=axs[i][j])
    plt.savefig("outputs/pictures/numerical_vs_categorical.svg", bbox_inches="tight")
    plt.close(fig)


def draw_correlation(data):
    fig, ax = plt.subplots(figsize=(30, 30))
    corr_df = data[numerical]
    corrMatrix = corr_df.corr()
    sns.heatmap(corrMatrix, annot=True)
    sns.set(rc={'figure.figsize': (30, 30)})
    plt.savefig("outputs/pictures/corr.svg")
    plt.close(fig)


if __name__ == "__main__":
    data = pd.read_csv("outputs/tables/patients.csv", index_col="patient_id")

    subsetting(data)
    draw_pairplots(data)
    draw_categorical_count(data)
    draw_violin_plots(data)
    draw_numerical_vs_categorical(data)
    draw_correlation(data)
