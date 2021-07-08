import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages



if __name__ == "__main__":
    data = pd.read_csv("patients.csv", index_col="patient_id")

    g = sns.pairplot(data, diag_kind="kde", palette=["#F4F1DE","#AA3F22"])
    g.map_lower(sns.kdeplot, levels=4, color=".4")
    plt.savefig("output/pairplots.svg")
    plt.close()



    numerical = ["chao1_before", "chao1_after", "shannon_before", "shannon_after", "observed_before", "observed_after", "ace_before", "ace_after", "berger_parker_before", "berger_parker_after", "read_count_before", "read_count_after"]
    categorical = ["enterotype_before", "enterotype_after", "treatment"]



    fig, ax = plt.subplots(1, 3, figsize=(20, 10))
    for variable, subplot in zip(categorical, ax.flatten()):
        sns.countplot(data[variable], ax=subplot)
        for label in subplot.get_xticklabels():
            label.set_rotation(90)
    plt.savefig("output/categorical_count.svg")
    plt.close()



    with PdfPages("output/violin_plots.pdf") as pdf:
        for metric in ["chao1", "shannon", "observed", "ace", "berger_parker"]:
            data_sub = data[["treatment", metric + "_before", metric + "_after"]].dropna()

            fig, axs = plt.subplots(3, 2, figsize=(10, 20))
            max_y = max(max(data_sub[metric + "_before"]), max(data_sub[metric + "_after"]))
            if max_y < 10:
                max_y = math.ceil(max_y)
            else:
                max_y = math.ceil(max_y / 100) * 100

            fig.suptitle(metric)

            axs[0][0].set_title("Before")
            axs[0][1].set_title("After")

            for i, treatment in enumerate(["Control", "STD2", "STD3"]):
                axs[i][0].violinplot(data_sub[metric + "_before"])
                axs[i][1].violinplot(data_sub[metric + "_after"])

                for ax in axs[i]:
                    ax.set_ylim(0, max_y)
                axs[i][0].set_ylabel(treatment)

            pdf.savefig(fig)
            plt.close(fig)



    fig, axs = plt.subplots(len(numerical), len(categorical), figsize=(50, 100))
    for i, numvar in enumerate(numerical):
        for j, catvar in enumerate(categorical):
            sns.boxplot(x=catvar, y=numvar, data=data, ax=axs[i][j])
    plt.savefig("output/numerical_vs_categorical.svg", bbox_inches="tight")
    plt.close(fig)



    fig, ax = plt.subplots(figsize=(30, 30))
    corr_df=data[numerical]
    corrMatrix = corr_df.corr()
    sns.heatmap(corrMatrix, annot=True)
    sns.set(rc={'figure.figsize':(30,30)})
    plt.savefig("output/corr.svg")
    plt.close(fig)
