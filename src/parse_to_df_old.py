import numpy as np
import re
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


def export_df(df, file_name):
    with open(file_name, 'w') as f:
        for block, block_df in df.groupby('block'):
            print(f'>{block}', file=f)
            for i, row in block_df.iterrows():
                print(f'{row["species"]}.{row["chr"]}:{row["chr_beg"]}-{row["chr_end"]} {row["orientation"]}', file=f)
            print(file=f)


def filter_unique_gene(df):
    allowed_blocks = set()
    all_sp = len(df['species'].unique())
    for block, df_block in df.groupby('block'):
        if len(df_block) == len(df_block['species'].unique()) == all_sp:
            allowed_blocks.add(block)

    return df[df.apply(lambda x: x['block'] in allowed_blocks, axis=1)]


def dist_between_blocks(df):
    ds = []
    for sp, df_sp in df.groupby('species'):
        df_sp = df_sp.sort_values(by=['chr_beg'])
        ds += (start_ - end_ for start_, end_ in zip(df_sp['chr_beg'][1:], df_sp['chr_end']))
    return ds

all_blocks = '/home/yulia/Downloads/great_again/summer_ipaH_story/data/blocks_coords.infercars'
common_blocks = '/home/yulia/Downloads/great_again/summer_ipaH_story/data/blocks_coords_unique_gene.infercars'

all_blocks = parse_to_df(all_blocks)
common_blocks = parse_to_df(common_blocks)

all_blocks.to_csv('/home/yulia/Downloads/great_again/summer_ipaH_story/results/all_blocks.tsv', sep='\t', index=None)
common_blocks.to_csv('/home/yulia/Downloads/great_again/summer_ipaH_story/results/common_blocks.tsv', sep='\t', index=None)

# in_file = '/home/yulia/Downloads/great_again/summer_ipaH_story/blocks_coords.infercars'
# blocks = [1360, 1362]
# df = parse_to_df(in_file)
# sorted_df = df.sort_values(by=['species', 'chr_beg'])

# df.to_csv('blocks_df.tsv', sep='\t')
# sorted_df.to_csv('sorted_df.tsv', sep='\t')
# print(df[(df.block.isin(blocks))])
# print(df[(df.block.isin(blocks))])
# print(df)
# print(sorted_df)
# sorted_df = sorted_df.reset_index()
# idxs = list(sorted_df[sorted_df.block.isin(blocks)].index)
# idxs_min = list(map(lambda x: x - 1, idxs))
# idxs_plus = list(map(lambda x: x + 1, idxs))
# all_idxs = idxs + idxs_min + idxs_plus
# ipah_blocks_and_neighbours = sorted_df.loc[[x for x in all_idxs]].sort_values(by=['species', 'chr_beg'])
# ipah_blocks_and_neighbours.to_csv('ipah_blocks_and_neighbours.tsv', sep='\t')
#
# from collections import Counter

# print(Counter(ipah_blocks_and_neighbours.block))
# counter_dict = {1360: 51, 1362: 35, 348: 31, 532: 30, 805: 22, 1366: 18, 764: 14, 531: 14, 1363: 13, 1710: 11, 615: 3,
#                 618: 3, 541: 2, 339: 2, 1907: 2, 1359: 2, 1165: 1, 1357: 1, 436: 1, 2979: 1, 1358: 1}
# data = {'block': list(counter_dict.keys()), 'frequency': list(counter_dict.values())}
# print(pd.DataFrame.from_dict(data))
# pd.DataFrame.from_dict(data).to_csv('neighbourgs_frequencies.tsv', sep='\t')
# print(len(set(ipah_blocks_and_neighbours.species))) # всего 35 видов, у кого есть блоки 1360 или 1362

# in_file = '/home/yulia/Downloads/great_again/summer_ipaH_story/blocks_coords/blocks_coords_unique_gene.infercars'

# df = parse_to_df(in_file)
# sorted_df = df.sort_values(by=['species', 'chr_beg'])
# df.to_csv('blocks_coords_unique_gene.tsv', sep='\t')
# sorted_df.to_csv('blocks_coords_unique_gene_sorted.tsv', sep='\t')
#
# ipah_blocks = pd.read_csv('1360_1362_coords.tsv', sep='\t', index_col=0)

# print(sorted_df)
# print(ipah_blocks)

# 16_07_2020
# frames = (sorted_df, ipah_blocks)
# result = pd.concat(frames).sort_values(by=['species', 'chr_beg'])
# result.to_csv('blocks_coords_unique_gene_and_ipah_blocks_sorted.tsv', sep='\t')
#
# result = result.reset_index()
# idxs = list(result[result.block.isin(blocks)].index)
# idxs_min = list(map(lambda x: x - 1, idxs))
# idxs_plus = list(map(lambda x: x + 1, idxs))
# all_idxs = idxs + idxs_min + idxs_plus
# ipah_blocks_and_neighbours = result.loc[[x for x in all_idxs]].sort_values(by=['species', 'chr_beg'])
# ipah_blocks_and_neighbours.to_csv('ipah_blocks_and_common_neighbours.tsv', sep='\t')
#
# from collections import Counter
#
# print(Counter(ipah_blocks_and_neighbours.block))
#
# counter_dict = {1360: 51, 1362: 35, 532: 30, 358: 30, 531: 29, 347: 29, 857: 19, 704: 17, 623: 3, 678: 3, 541: 2, 701: 2, 339: 2, 487: 1, 529: 1, 436: 1, 462: 1, 434: 1, 1082: 1}
# data = {'block': list(counter_dict.keys()), 'frequency': list(counter_dict.values())}
# pd.DataFrame.from_dict(data).to_csv('common_neighbourgs_frequencies.tsv', sep='\t')
# print(len(set(ipah_blocks_and_neighbours.species))) # всего 35 видов, у кого есть блоки 1360 или 1362


