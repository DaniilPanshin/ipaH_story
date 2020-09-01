import argparse
import numpy as np
import pandas as pd
from io import StringIO

# dir = 'example_data/SibeliaZ/hg38_mm10/far_params_12/5000/'

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert Ragout output to infercars format')
    parser.add_argument('-dir', '--dir', type=str, help='Path to ragout_out directory')
    args = parser.parse_args()

    dir = args.dir
    file = 'blocks_coords.txt'
    file_inf = 'blocks_coords.infercars'
    sep = '-' * 80
    split_by_underscore = False

    lines = open(dir + file).read().split('\n')[:-1]
    ls = np.array(lines)
    bs = np.split(ls, np.where(ls == sep)[0])

    # names of chromosomes
    df_names = pd.read_csv(StringIO('\n'.join(bs[0])), sep='\t')
    chr_names = {}
    for index, row in df_names.iterrows():
        chr_names[row['Seq_id']] = row['Description']

    # blocks data
    with open(dir + file_inf, 'w') as f:
        for b in bs[1:-1]:
            block = b[1].split('#')[1]
            df_block = pd.read_csv(StringIO('\n'.join(b[2:])), sep='\t')

            print(f'>{block}', file=f)
            for index, row in df_block.iterrows():
                chr_name, start, end, strand = chr_names[row['Seq_id']], row['Start'], row['End'], row['Strand']
                if split_by_underscore:
                    splitted = chr_name.split('_', 1)
                    chr = splitted[0] + '.' + splitted[1].split('.')[0]
                else:
                    chr = chr_name
                if not 'alt' in chr:
                    print(chr + ':' + (f'{start}-{end}' if strand == '+' else f'{end}-{start}') + ' ' + strand, file=f)
            print(file=f)

