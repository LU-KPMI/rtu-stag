#!/bin/python3

import os
import re
import random
import math


def shannon(v):
    s = sum(v)
    return sum(-(x/s) * math.log(x/s) if x > 0 else 0 for x in v)


def observed(v):
    return sum(1 if x > 0 else 0 for x in v)


def chao1(v):
    obs = observed(v)
    n1 = sum(1 if x == 1 else 0 for x in v)
    n2 = sum(1 if x == 2 else 0 for x in v)

    return obs + n1 * (n1 - 1) / (2 * (n2 + 1))


def rarefy(v, f):
    res = []
    step = 5
    sizes = [sum(v) * f // 100 for f in range(step, 101, step)]
    cur = [0] * len(v)
    for s in sizes:
        x = []
        for i in range(len(v)):
            x += [i] * (v[i] - cur[i])
        sample = random.sample(x, s - sum(cur))
        for t in sample:
            cur[t] += 1
        res.append([sum(cur) / sum(v), f(cur)])
    return res


def read_bracken_file(filename):
    v = []
    with open(filename) as f:
        first = True
        for line in f:
            if not first:
                data = line.rstrip().split('\t')
                n = int(data[5])
                v.append(n)
            first = False
    return v


print("sample_name", "chao1", "shannon", "observed", sep=",")

for filename in os.listdir("./bracken_1/"):
    if not filename.endswith(".bracken"):
        continue
    sample_name = re.sub(".bracken", "", filename)

    v = read_bracken_file(os.path.join("bracken_1", filename))
    print(sample_name, chao1(v), shannon(v), observed(v), sep=",")

