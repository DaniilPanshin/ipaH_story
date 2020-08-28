import glob
import argparse
from Bio import SeqIO
import pandas as pd


def make_stats(path):
    with open('total_stats.tsv', 'w') as out_f:
        out_f.write(f"AssemblyID\tAssemblyName\tType\tGenBankID\tOrganism\tLength(bp)\tPlasmid\n")
        for pth in path:
            records = [el for el in SeqIO.parse(pth, 'genbank')]

            pth = pth.split('/')[-1]
            print(pth)
            acc_id_list = pth.split('_')

            acc_id = f"{acc_id_list[0]}_{acc_id_list[1]}"
            ass_name = acc_id_list[2]

            for rec in records:
                organism = rec.features[0].qualifiers['organism'][0]
                if 'strain' in rec.features[0].qualifiers:
                    if rec.features[0].qualifiers['strain'][0] == str(organism.split(' ')[-1]):
                        strain = ''
                    else:
                        strain = rec.features[0].qualifiers['strain'][0]
                else:
                    strain = ''
                animal = f"{organism} {strain}"
                if 'plasmid' in rec.description:
                    out_f.write(
                        f"{acc_id}\t{ass_name}\tPlasmid\t{rec.id}\t{animal}\t{len(rec)}\t{'yes' if len(records) > 1 else 'no'}\n")
                else:
                    out_f.write(
                        f"{acc_id}\t{ass_name}\tChromosome\t{rec.id}\t{animal}\t{len(rec)}\t{'yes' if len(records) > 1 else 'no'}\n")


def sort_table_by_cromosome(df):
    df = df.sort_values(by=['Type', 'AssemblyID'])
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make total_stats.tsv table for list of GenBank file')
    parser.add_argument('-gb', '--genbank', type=str, nargs='*', help='Input GenBank files')
    args = parser.parse_args()

    files = args.genbank
    make_stats(files)
    out_df = sort_table_by_cromosome(pd.read_csv('total_stats.tsv', sep='\t'))
    out_df.to_csv('total_stats.tsv', sep='\t', index=None)
