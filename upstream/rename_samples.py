import os
from glob import glob
import subprocess
import pandas as pd

SOURCE_DIR = "/mnt/home/groups/lu_kpmi/raw_mgi_data"
OUTPUT_DIR = "/mnt/home/groups/lu_kpmi/renamed_samples"
SAMPLE_LIST_DIR = "/mnt/home/groups/lu_kpmi/sample_list_data"


def get_list(s):
    if '-' in s:
        a, b = map(int, s.split('-'))
        return range(a, b + 1)
    else:
        return [int(x) for x in s.split(',')]


def get_old_name(row):
    return row[6] + "_L" + f"{int(row[7]):02d}"


if __name__ == "__main__":
    processed = set(os.listdir(OUTPUT_DIR))

    for filename in os.listdir(SAMPLE_LIST_DIR):
        df = pd.read_csv(os.path.join(SAMPLE_LIST_DIR, filename),
                         header=None,
                         delimiter=';' if "deleted" in filename else ',')
        for _, row in df.iterrows():
            old_name, new_name = get_old_name(row), row[10]
            print("Processing", old_name, "->", new_name)
            if new_name + "_1.fq.gz" in processed:
                print("Skipping...")
                continue
            nums = get_list(row[5])
            for d in [1, 2]:
                files = []
                for t in nums:
                    orig = old_name + "_" + str(t) + "_" + str(d) + ".fq.gz"
                    files.append(glob(os.path.join(SOURCE_DIR, "**", orig),
                                      recursive=True)[0])
                print(files)
                dest = os.path.join(OUTPUT_DIR,
                                    new_name + "_" + str(d) + ".fq.gz")
                print("to", dest)
                with open(dest, 'w') as out:
                    subprocess.run(["cat"] + files, stdout=out)
