import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def add_to_pdf(pdf, title, data):
    plt.scatter(data["read_count"], data["observed"], 1)
    plt.title(title)
    plt.xlabel("read count")
    plt.ylabel("observed")
    plt.xlim(0, 1.2e8)
    plt.ylim(0, 1050)
    pdf.savefig()
    plt.close()

data = pd.read_csv("./richness.csv", index_col="sample_name")
metadata = pd.read_csv("./sample_metadata.csv", index_col="sample_name")

data = data.merge(metadata, left_index=True, right_index=True)

with PdfPages("graphs.pdf") as pdf:
    add_to_pdf(pdf, "Visi", data)

    for treatment in ["Control", "STD2", "STD3"]:
        for time in ["T1", "T2"]:
            add_to_pdf(pdf, str(treatment) + " & " + str(time), data[(data["treatment"] == treatment) & (data["time"] == time)])

    data_before = data[data["time"] == "T1"].set_index("patient_id").drop("time", 'columns')
    data_after = data[data["time"] == "T2"].set_index("patient_id").drop("time", 'columns')

    keep = data_before.index.intersection(data_after.index)

    print("Patients missing data after therapy:", list(data_before[~data_before.index.isin(keep)].index))
    print("Patients missing data before therapy:", list(data_after[~data_after.index.isin(keep)].index))

    data_before = data_before[data_before.index.isin(keep)]
    data_after = data_after[data_after.index.isin(keep)]


    for metric in ["chao1", "shannon", "observed"]:
        fig, axs = plt.subplots(3, 2)

        fig.suptitle(metric)

        axs[0][0].set_title("Before")
        axs[0][1].set_title("After")

        for i, treatment in enumerate(["Control", "STD2", "STD3"]):
            axs[i][0].violinplot(data_before[data_before["treatment"] == treatment][metric])
            axs[i][1].violinplot(data_after[data_after["treatment"] == treatment][metric])

            for ax in axs[i]:
                ax.set_ylim(0, 5 if metric == "shannon" else 1500)
            axs[i][0].set_ylabel(treatment)

        pdf.savefig(fig)
