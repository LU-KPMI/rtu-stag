import argparse
import pathlib
import os
import shutil


def copy_subdir(sample_dir, subdir, dest_dir):
    for run in os.listdir(sample_dir):
        if subdir in os.listdir(os.path.join(sample_dir, run)):
            shutil.copytree(os.path.join(sample_dir, run, subdir),
                            os.path.join(dest_dir, subdir))
            return
    print("Warning: Cannot find subdirectory " + subdir + " in " + sample_dir)


def extract_sample(sample_dir, dest_dir):
    if not os.path.exists(sample_dir):
        print("Error: sample " + sample_dir + " doesn't exist")
        exit(0)
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    os.mkdir(dest_dir)
    copy_subdir(sample_dir, "kraken2", dest_dir)
    copy_subdir(sample_dir, "amrplusplus", dest_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Populates directory with samples. Sample list should be in file samples.txt')
    parser.add_argument("dest_dir",
                        type=pathlib.Path,
                        help="Destination directory")
    args = parser.parse_args()

    dest_dir = args.dest_dir

    with open(os.path.join(dest_dir, "samples.txt")) as sample_list:
        for sample in sample_list:
            sample = sample.rstrip()
            extract_sample(os.path.join("/home/groups/lu_kpmi/outputs/",
                                        sample),
                           os.path.join(dest_dir, sample))
