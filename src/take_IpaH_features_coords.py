import os
import argparse
from Bio import SeqIO


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Write all IpaH features to table with coordinates')
    parser.add_argument('-i', '--input', type=str, nargs='*', help='Path to .gbff files')
    args = parser.parse_args()

    files = args.input

    with open('IpaH_features.tsv', 'w') as out_f:
        out_f.write(f"AssemblyID\tGenBankID\tFeatureName\tFeatureStart\tFeatureEnd\n")
        for file in files:
            assembly_id = '_'.join(file.split('/')[-1].split('_')[0:2]).split('.')[0]
            recs = SeqIO.parse(file, 'genbank')
            for rec in recs:
                if 'plasmid' in rec.description:
                    pass
                else:
                    gb_id = rec.name
                    print(f"Processing...\t{gb_id}")
                    for feature in rec.features:
                        for val in feature.qualifiers.values():
                            if 'IpaH' in val[0] or 'ipaH' in val[0] or 'invasion plasmid antigen' in val[0]:
                                f_start = feature.location.start
                                f_end = feature.location.end
                                out_f.write(f"{assembly_id}\t{gb_id}\t{val[0]}\t{f_start}\t{f_end}\n")

    print(f"Output is written to {os.getcwd()}/IpaH_features.tsv file :>")