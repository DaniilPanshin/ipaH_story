import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqFeature import SeqFeature, FeatureLocation


def check_coods(start, end):
    if start > end:
        return end, start
    else:
        return start, end


header_list = ['qseqid', 'sseqid', 'bitscore', 'pident', 'qcovhsp', 'length', 'evalue',
               'qstart', 'qend', 'sstart', 'send', 'sframe']

df = pd.read_csv('tblastn_result.tsv', sep='\t', header=None, comment='#', names=header_list)
df_filter_evalue = df[df.evalue == 0.0]

for gene in set(df_filter_evalue.qseqid):
    if gene == 'invasion':
        pass
    else:
        # print(gene)
        df_qseqid = df_filter_evalue[df_filter_evalue.qseqid == gene]

        plt.figure(figsize=(30, 10))

        plt.subplot(2, 2, 1)
        plt_bitscore = df_qseqid.bitscore.hist(bins=100)
        plt_bitscore.set_title("bitscore")

        plt.subplot(2, 2, 2)
        plt_pident = df_qseqid.pident.hist(bins=100)
        plt_pident.set_title("pident")

        plt.subplot(2, 2, 3)
        plt_qcovhsp = df_qseqid.qcovhsp.hist(bins=100)
        plt_qcovhsp.set_title("qcovhsp")

        plt.subplot(2, 2, 4)
        plt_qcovhsp = df_qseqid.length.hist(bins=100)
        plt_qcovhsp.set_title("length")

        plt.suptitle(f"Distributions for {gene}")

        # plt.show()
        plt.savefig(f"{gene}.png", fromat="png")
