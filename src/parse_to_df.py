import re
import argparse
import numpy as np
import pandas as pd

p = "([A-Za-z0-9_\(\)\/\s]+)\.([A-Za-z0-9_\.]+):(\d+)-(\d+) ([+|-]).*"
pattern = re.compile(p)
columns = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]


def find_indices(lst, condition):
    return [i for i, elem in enumerate(lst) if condition(elem)]


def parse_to_df(file_name):
    with open(file_name) as f:
        lines = f.readlines()

    last_line = len(lines) - 1
    while lines[last_line] == '\n': last_line -= 1

    n_at_end = len(lines) - 1 - last_line
    for _ in range(1 - n_at_end): lines.append('\n')

    bs = np.split(lines, find_indices(lines, lambda x: x[0] == ">"))
    temp = []

    for i, b in enumerate(bs):
        if len(b) == 0: continue
        b_i = int(b[0][1:])

        for oc in b[1:-1]:
            m = pattern.match(oc)
            temp.append([b_i, m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)])

    return pd.DataFrame(temp, columns=columns)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert .infercars format to .tsv')
    parser.add_argument('-i', '--input', type=str, help='Path to input file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    args = parser.parse_args()

    in_file = args.input
    df = parse_to_df(in_file)
    df.to_csv(args.output, sep='\t', index=None)

