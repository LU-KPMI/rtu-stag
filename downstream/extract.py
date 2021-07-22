# TODO: This file shouldn't exist in a well planned program:
#       * Bracken should be done before downstream analysis;
#       * There should already be a file containing patient_id, sample_id_before, sample_id_after
#         and various metadata instead of hacking it together from sample ids.
import os
import shutil
import subprocess


group_outputs_dir = os.path.join(os.getenv("GROUP"), "outputs")
bracken_db_path = os.path.join(os.getenv("GROUP"), "databases", "full_ref_bafp", "database150mers.kmer_distrib")


def get_all_sample_ids():
    return sorted(list(filter(lambda x : x.startswith("LV"), os.listdir(group_outputs_dir))))


def print_metadata(output_filename):
    m = {}
    for sample_id in get_all_sample_ids():
        patient_id, treatment, time = sample_id.split('_')
        if patient_id not in m:
            m[patient_id] = [treatment, "", ""]
        if time == "T1":
            m[patient_id][1] = sample_id
        elif time == "T2":
            m[patient_id][2] = sample_id
        else:
            raise RuntimeError("Wrong time in sample name")
    with open(output_filename, "w") as f:
        print("patient_id", "treatment", "sample_id_before", "sample_id_after", sep = ',', file = f)
        for patient_id in m:
            print(patient_id, *(m[patient_id]), sep = ',', file = f)


def get_kraken_filename(sample_id):
    sample_subdir = os.path.join(group_outputs_dir, sample_id)
    for trial in os.listdir(sample_subdir):
        kraken_dir = os.path.join(sample_subdir, trial, "kraken2")
        if os.path.isdir(kraken_dir):
            for s in ["1", "2", sample_id]:
                kraken_filename = os.path.join(kraken_dir, s + ".kreport")
                if os.path.isfile(kraken_filename):
                    return kraken_filename
            raise AssertionError("Can't find .kreport file in kraken2 subdirectory")
    return None


def get_amrplusplus_filename(sample_id):
    sample_subdir = os.path.join(group_outputs_dir, sample_id)
    for trial in os.listdir(sample_subdir):
        amrplusplus_dir = os.path.join(sample_subdir, trial, "amrplusplus", "RunResistome")
        if os.path.isdir(amrplusplus_dir):
            for f in os.listdir(amrplusplus_dir):
                if f.endswith(".mechanism.tsv"):
                    return os.path.join(amrplusplus_dir, f)
            raise AssertionError("Can't find .mechanism.tsv file im amrplusplus subdirectory")
    return None


def run_bracken(kraken_filename, output_filename, level):
    proc = subprocess.run(["../../Bracken/src/est_abundance.py",
        "-i", kraken_filename,
        "-k", bracken_db_path,
        "-o", output_filename,
        "-t", "1",
        "-l", level,
        "--out-report", "./discard"])
    if proc.returncode != 0:
        raise RuntimeError("Bracken failed")
    os.remove("./discard")


def create_bracken_reports():
    os.mkdir("bracken_output")
    for sample_id in get_all_sample_ids():
        kraken_filename = get_kraken_filename(sample_id)
        if kraken_filename == None:
            print(sample_id, " missing kraken2, skipping sample")
            continue

        out_name = os.path.join("bracken_output", str(sample_id) + ".bracken")
        run_bracken(kraken_filename, out_name, "G")


def extract_amrplusplus_reports():
    os.mkdir("amrplusplus_report")
    for sample_id in get_all_sample_ids():
        amrplusplus_filename = get_amrplusplus_filename(sample_id)
        if amrplusplus_filename == None:
            print(sample_id, " missing amrplusplus, skipping sample")
            continue
        shutil.copyfile(amrplusplus_filename, os.path.join("amrplusplus_report", str(sample_id) + ".tsv"))


if __name__ == "__main__":
    print_metadata("metadata.csv")
    create_bracken_reports()
    extract_amrplusplus_reports()