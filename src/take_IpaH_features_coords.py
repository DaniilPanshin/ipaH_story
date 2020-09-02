import os
import argparse
import pandas as pd
from Bio import SeqIO


def check_CDS_and_genes(f, key_words):
    if f.type == 'CDS' or f.type == 'gene':
        for k in f.qualifiers:
            for word in key_words:
                if word in f.qualifiers[k][0]:
                    return True


def check_q_values(f, key_words):
    for val in f.qualifiers.values():
        for word in key_words:
            if word in val[0]:
                return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Write all IpaH features to table with coordinates')
    parser.add_argument('-i', '--input', type=str, nargs='*', help='Path to .gbff files')
    parser.add_argument('-o', '--output', type=str, default='IpaH_features.tsv',
                        help='Path to output file. By default file will be written to current working directory under "IpaH_features.tsv" name')
    args = parser.parse_args()

    files = args.input
    key_words = ('IpaH', 'ipaH', 'invasion plasmid antigen')
    restriction_words = ('IpaJ', 'ipaJ')
    df = pd.DataFrame(columns=['AssemblyID', 'GenBankID', 'ChrType', 'FeatureName', 'FeatureStart', 'FeatureEnd'])

    count = 0
    for file in files:
        count += 1
        assembly_id = '_'.join(file.split('/')[-1].split('_')[0:2]).split('.')[0]
        recs = SeqIO.parse(file, 'genbank')
        for rec in recs:
            if 'plasmid' in rec.description:
                chr_type = 'Plasmid'
            else:
                chr_type = 'Chromosome'
            gb_id = rec.name
            print(f"Processing... {count}/{len(files)}\t{assembly_id}")

            chromosomal_features = []
            for f in rec.features:
                if check_CDS_and_genes(f, key_words) or check_q_values(f, key_words):
                    chromosomal_features.append(f)

            # print(len(set(chromosomal_features)))
            for feat in set(chromosomal_features):
                f_start = feat.location.start
                f_end = feat.location.end

                if feat.type == 'CDS' or feat.type == 'gene':
                    if 'gene' in feat.qualifiers:
                        f_name = feat.qualifiers['gene'][0]
                    elif 'product' in feat.qualifiers:
                        f_name = feat.qualifiers['product'][0]
                else:
                    f_name = feat.qualifiers['note'][0]

                if f_name not in restriction_words:
                    df = df.append({'AssemblyID': assembly_id,
                                    'GenBankID': gb_id,
                                    'FeatureName': f_name,
                                    'ChrType': chr_type,
                                    'FeatureStart': f_start,
                                    'FeatureEnd': f_end},
                                   ignore_index=True)

    df = df.sort_values(by=['AssemblyID', 'ChrType', 'FeatureStart'])
    df = df.drop_duplicates().reset_index(drop=True)
    df.to_csv(args.output, sep='\t')
    print(f"Output is written to {os.getcwd()}/{args.output} file :>")
