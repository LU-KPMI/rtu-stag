import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math

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
metadata = pd.read_csv("./metadata.csv", index_col="patient_id")

with PdfPages("graphs.pdf") as pdf:
    add_to_pdf(pdf, "Visi", data)

    for treatment in ["Control", "STD2", "STD3"]:
        for time, row_name in [("T1", "sample_id_before"), ("T2", "sample_id_after")]:
            sample_ids = set(metadata[metadata["treatment"] == treatment][row_name].dropna())
            add_to_pdf(pdf, str(treatment) + " & " + str(time), data[data.index.isin(sample_ids)])

    metadata.dropna(inplace=True)
    data_before = metadata.\
                    merge(data, left_on = "sample_id_before", right_index = True).\
                    drop(["sample_id_before", "sample_id_after"], "columns")
    data_after  = metadata.\
                    merge(data, left_on = "sample_id_after", right_index = True).\
                    drop(["sample_id_before", "sample_id_after"], "columns")

    metrics = list(data_before.columns)
    metrics.remove("treatment")
    metrics.remove("read_count")

    for metric in metrics:
        fig, axs = plt.subplots(3, 2, figsize=(10, 20))

        max_y = max(max(data_before[metric]), max(data_after[metric]))
        if max_y < 10:
            max_y = math.ceil(max_y)
        else:
            max_y = math.ceil(max_y / 100) * 100

        fig.suptitle(metric)

        axs[0][0].set_title("Before")
        axs[0][1].set_title("After")

        for i, treatment in enumerate(["Control", "STD2", "STD3"]):
            axs[i][0].violinplot(data_before[data_before["treatment"] == treatment][metric])
            axs[i][1].violinplot(data_after[data_after["treatment"] == treatment][metric])

            for ax in axs[i]:
                ax.set_ylim(0, max_y)
            axs[i][0].set_ylabel(treatment)

        pdf.savefig(fig)
        plt.close(fig)
