import os

if __name__ == "__main__":
    sample_to_abundance = {}
    all_genuses = set([])
    for filename in os.listdir("bracken"):
        if not filename.endswith(".bracken"):
            continue
        with open(os.path.join("bracken", filename)) as f:
            abundance = {}
            first = True
            for line in f:
                if first:
                    first = False
                    continue
                data = line.rstrip().split('\t')
                genus, cnt = data[0], int(data[5])
                all_genuses.add(genus)
                abundance[genus] = cnt
            sample_to_abundance[filename.rstrip(".bracken")] = abundance

    all_genuses = sorted(all_genuses)

    print("", *all_genuses, sep=',')

    for sample in sample_to_abundance:
        abundance = sample_to_abundance[sample]
        v = [abundance[g] if g in abundance else 0 for g in all_genuses]
        print(sample, *v, sep=",")

