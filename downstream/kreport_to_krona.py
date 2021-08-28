import sys

full_name = []
level_map = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5}

if __name__ == "__main__":
    for line in sys.stdin:
        entries = line.rstrip().split('\t')
        if entries[3] not in level_map:
            continue
        level = level_map[entries[3]]

        while len(full_name) > level:
            full_name.pop()

        full_name.append(entries[5].lstrip().rstrip())
        if level == 5:
            print(entries[1] + "\t" + "\t".join(full_name))
