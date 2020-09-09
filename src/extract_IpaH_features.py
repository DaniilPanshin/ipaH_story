import argparse
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


def define_f_name(f):
    if feat.type == 'CDS' or feat.type == 'gene':
        if 'gene' in feat.qualifiers:
            f_name = feat.qualifiers['gene'][0]
            return f_name
        elif 'product' in feat.qualifiers:
            f_name = feat.qualifiers['product'][0]
            return f_name
    else:
        f_name = feat.qualifiers['note'][0]
        return f_name


def filter_restrictions(f_name, restriction_words):
    result = set()
    for word in restriction_words:
        if word in f_name:
            result.add(True)
        else:
            result.add(False)
    if True in result:
        return False
    else:
        return True


if __name__ == "__main__":

    # GCA_000006925 is a reference Shigella

    parser = argparse.ArgumentParser(description='Write all IpaH features to table with coordinates')
    parser.add_argument('-i', '--input', type=str, nargs='*', help='Path to .gbff files')
    parser.add_argument('-o', '--output', type=str, default='IpaH_features.fa',
                        help='Path to output file. By default file will be written to current working directory under "IpaH_features.tsv" name')
    args = parser.parse_args()

    files = args.input
    key_words = ('IpaH', 'ipaH', 'invasion plasmid antigen')
    restriction_words = ('IpaJ', 'ipaJ', 'SpaS')

    with open(args.output, 'w') as out_f:
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

                features = []
                for f in rec.features:
                    if check_CDS_and_genes(f, key_words) or check_q_values(f, key_words):
                        features.append(f)

                for feat in features:
                    if feat.type == 'CDS':
                        f_name = define_f_name(f)
                        if filter_restrictions(f_name, restriction_words):
                            new_feat = feat.extract(rec)
                            new_feat.id = f_name
                            new_feat.description = chr_type
                            SeqIO.write(new_feat, out_f, 'fasta')