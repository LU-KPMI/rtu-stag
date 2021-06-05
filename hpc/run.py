#!/bin/python3

import os
import re
import sys
import argparse, textwrap
import subprocess

group_dir = "/home/groups/lu_kpmi"


def list_all_samples():
    samples = []

    files = os.listdir(os.path.join(group_dir, "renamed_samples"))
    files = [re.sub('_[1-2].fq.gz', '', f) for f in files]
    sample_names = sorted(set(files))

    samples += [(s, group_dir + '/renamed_samples/' + s + '_1.fq.gz', group_dir + '/renamed_samples/' + s + '_2.fq.gz') for s in sample_names]

    for d in os.listdir(os.path.join(group_dir, 'raw_mgi_data/MMHP_Latvia_20210514')):
        if d.endswith('.md5'):
            continue
        read_1 = ''
        read_2 = ''
        for f in os.listdir(os.path.join(group_dir, 'raw_mgi_data/MMHP_Latvia_20210514', d)):
            if f.endswith('_1.fq.gz'):
                read_1 = os.path.join(group_dir, 'raw_mgi_data/MMHP_Latvia_20210514', d, f)
            elif f.endswith('_2.fq.gz'):
                read_2 = os.path.join(group_dir, 'raw_mgi_data/MMHP_Latvia_20210514', d, f)
        assert(read_1 != '' and read_2 != '')
        samples.append((d, read_1, read_2))

    return samples


def sample_has_subfolder(sample_name, subfolder_name):
    if not os.path.isdir(os.path.join(group_dir, 'outputs', sample_name)):
        return False
    for d in os.listdir(os.path.join(group_dir, 'outputs', sample_name)):
        if os.path.isdir(os.path.join(group_dir, 'outputs', sample_name, d, subfolder_name)):
            return True
    return False


def sample_has_kraken(sample_name):
    return sample_has_subfolder(sample_name, 'kraken2')


def sample_has_amrplusplus(sample_name):
    return sample_has_subfolder(sample_name, 'amrplusplus')


def sample_has_groot(sample_name):
    return sample_has_subfolder(sample_name, 'groot')


def plus_or_minus(b):
    return '+' if b else '-'


def print_sample_stats(sample_name):
    if not os.path.isdir(os.path.join(group_dir, 'outputs', sample_name)):
        print('{sample: <25}                     n/a'.format(sample=sample_name))
    else:
        kraken = plus_or_minus(sample_has_kraken(sample_name))
        amrplusplus = plus_or_minus(sample_has_amrplusplus(sample_name))
        groot = plus_or_minus(sample_has_groot(sample_name))
        print('{sample: <25} {kraken2: >11} {amrplusplus: >11} {groot: >11}'.format(sample=sample_name, kraken2=kraken, amrplusplus=amrplusplus, groot=groot))


def print_stats():
    print('{sample: <25} {kraken2: >11} {amrplusplus: >11} {groot: >11}'.format(sample='Sample name', kraken2='kraken2', amrplusplus='amrplusplus', groot='groot'))
    l = list_all_samples()
    for s in l:
        print_sample_stats(s[0])
    print('Total: {}'.format(len(l)))


def diagnostics():
    sample_names = [s[0] for s in list_all_samples()]
    for f in os.listdir(os.path.join(group_dir, 'outputs')):
        if f not in sample_names:
            print("Folder {f} doesn't correspond to any sample!".format(f=f))


def print_all_samples():
    for s in list_all_samples():
        print(s[0])


def list_mode(args):
    print_all_samples()


def diagnostics_mode(args):
    diagnostics()


def filter_mode(args):
    for s in sys.stdin:
        s = s.rstrip()
        if args.kraken and sample_has_kraken(s):
            continue
        if args.amrplusplus and sample_has_amrplusplus(s):
            continue
        if args.groot and sample_has_groot(s):
            continue
        print(s)


def stats_mode(args):
    print_stats()


def run_mode(args):
    all_samples = list_all_samples()
    unprocessed = []

    with open(args.filename) as file:
        for sample_name in file:
            sample_name = sample_name.rstrip()
            read_1 = ''
            read_2 = ''
            for s in all_samples:
                if s[0] == sample_name:
                    read_1 = s[1]
                    read_2 = s[2]
                    break

            assert(read_1 != '' and read_2 != '')

            qsub_params = '{name} {read_1} {read_2} {work_path} {taxon_db_path} \
{human_ref_path} {resistome_path} {output_path} {enable_qc_reads} {enable_host_removal} \
{enable_kraken2} {enable_groot} {enable_amrplusplus} {bracken_threshold}'.format(
                            name=sample_name,
                            read_1=read_1,
                            read_2=read_2,
                            work_path=os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')),
                            taxon_db_path="/home/groups/lu_kpmi/databases/full_ref_bafp",
                            human_ref_path="/home/groups/lu_kpmi/databases/human_reference",
                            resistome_path="/home/groups/lu_kpmi/databases/groot_db/arg-annot_index",
                            output_path="/home/groups/lu_kpmi/outputs",
                            enable_qc_reads=True,
                            enable_host_removal=True,
                            enable_kraken2=args.kraken,
                            enable_groot=args.groot,
                            enable_amrplusplus=str(args.amrplusplus),
                            bracken_threshold=10)

            command = subprocess.run(['qsub', 'subscripts/sub.run.sh', '-F', qsub_params])

            if command.returncode == 0:
                print(sample_name, "in queue")
            else:
                print(sample_name, "couldn't run")
                unprocessed.append(sample_name)

    with open(args.filename, 'w') as file:
        file.write('\n'.join(unprocessed))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Main interface to sample processing.',
            formatter_class=argparse.RawTextHelpFormatter,
            epilog=textwrap.dedent('''
            Example usage:
                Do amrplusplus analysis for samples which don't have amrplusplus analysis yet:

                ./run.py list | ./run.py filter --amrplusplus > list
                ./run.py run --disable-kraken --disable-groot list

                Stats:

                ./run.py stats | less

                Diagnostics:

                ./run.py diagnostics
            '''))
    parser.set_defaults(func=lambda x : exit(0))
    subparsers = parser.add_subparsers()

    parser_list = subparsers.add_parser('list')
    parser_list.set_defaults(func=list_mode)

    pareser_diagnostics = subparsers.add_parser('diagnostics')
    pareser_diagnostics.set_defaults(func=diagnostics_mode)

    parser_filter = subparsers.add_parser('filter')
    parser_filter.set_defaults(func=filter_mode)
    parser_filter.add_argument('--kraken', action='store_true')
    parser_filter.add_argument('--amrplusplus', action='store_true')
    parser_filter.add_argument('--groot', action='store_true')

    parser_stats = subparsers.add_parser('stats')
    parser_stats.set_defaults(func=stats_mode)

    parser_run = subparsers.add_parser('run')
    parser_run.set_defaults(func=run_mode)
    parser_run.add_argument('filename')
    parser_run.add_argument('--disable-kraken', dest='kraken', action='store_false')
    parser_run.add_argument('--disable-amrplusplus', dest='amrplusplus', action='store_false')
    parser_run.add_argument('--disable-groot', dest='groot', action='store_false')

    args = parser.parse_args()
    args.func(args)

