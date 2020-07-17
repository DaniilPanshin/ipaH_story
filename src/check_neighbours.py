import pandas as pd
from collections import Counter

all_blocks = '/home/yulia/Downloads/great_again/summer_ipaH_story/results/all_blocks.tsv'
common_blocks = '/home/yulia/Downloads/great_again/summer_ipaH_story/results/common_blocks.tsv'
total_stats = '/home/yulia/Downloads/great_again/summer_ipaH_story/data/total_stats.csv'
outf_prefix = '/home/yulia/Downloads/great_again/summer_ipaH_story/results/'


def get_ids_for_lineage(dict_with_ids_ans_names, lineage):
    lineage_ids = []
    for k in dict_with_ids_ans_names:
        if lineage in k:
            lineage_ids.append(dict_with_ids_ans_names[k][:-2])
    return lineage_ids


def select_coord_by_species(df, ids):
    frames = []
    for acc in ids:
        frames.append(df[df['species'] == acc])
    result = pd.concat(frames).sort_values(by=['species', 'chr_beg'])
    return result


def take_blocks_and_neighbours(df, blocks):
    """
    return df with coordinates for each block from list and coordinates of block neighbours
    """
    df = df.reset_index()
    idxs = list(df[df['block'].isin(blocks)].index)
    idxs_min = list(map(lambda x: x - 1, idxs))
    idxs_plus = list(map(lambda x: x + 1, idxs))
    all_idxs = idxs + idxs_min + idxs_plus
    blocks_and_neighbours = df.loc[[x for x in all_idxs]]
    return blocks_and_neighbours.sort_values(by=['species', 'chr_beg'])


def check_neighours_frequencies(df):
    block_freq = dict(Counter(df.block))
    data = {'block': list(block_freq.keys()), 'frequency': list(block_freq.values())}
    return pd.DataFrame.from_dict(data)


def write_to_files(lineage, dict_with_ids_and_names, df_with_blocks, interesting_blocks, out_prefix, out_suffix):
    accessions = get_ids_for_lineage(dict_with_ids_and_names, lineage)
    species_df = select_coord_by_species(df_with_blocks, accessions)
    ipah_and_neighbours_coords = take_blocks_and_neighbours(species_df, interesting_blocks)
    ipah_and_neighbours_coords.to_csv(f"{out_prefix}{sp}/ipah_blocks_and_{out_suffix}_neighbours_coords.tsv",
                                      sep='\t', index=None)
    ipah_and_neighbours_frequencies = check_neighours_frequencies(ipah_and_neighbours_coords).\
        sort_values(by='frequency', ascending=False)
    ipah_and_neighbours_frequencies.to_csv(f"{out_prefix}{sp}/ipah_blocks_and_{out_suffix}_neighbours_frequencies.tsv",
                                           sep='\t', index=None)


# read files
total_stats = pd.read_csv(total_stats)
all_blocks = pd.read_csv(all_blocks, sep='\t')
common_blocks = pd.read_csv(common_blocks, sep='\t')

# take part of df with shigella info
shigella_df = total_stats[total_stats['Organism'].str.contains('Shigella')]

# make dict with shigella names and ids
ids_and_names = dict(zip(shigella_df['Organism'], shigella_df['AssemblyID']))

lineages = ('flexneri', 'boydii', 'dysenteriae', 'sonnei')
ipah_blocks = [1360, 1362]
ipah_coords = all_blocks[all_blocks['block'].isin(ipah_blocks)]
common_blocks_and_ipah = pd.concat([common_blocks, ipah_coords]).sort_values(by=['species', 'chr_beg'])

for sp in lineages:
    write_to_files(sp, ids_and_names, all_blocks, ipah_blocks, outf_prefix, 'all')
    write_to_files(sp, ids_and_names, common_blocks_and_ipah, ipah_blocks, outf_prefix, 'common')