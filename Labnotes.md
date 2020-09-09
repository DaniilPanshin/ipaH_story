### 0. Data preprocessing

#### Make table of available genomes

<pre><code>python /src/make_total_stats.py --genbank /path_to/genome_assemblies/gbff/*.gbff</code></pre> or <pre><code>python /src/make_total_stats.py -gb /path_to/genome_assemblies/gbff/*.gbff</code></pre>

Output will be written to *total_stats.tsv* file

#### Make table of IpaH/invasion plasmid antigen features

<pre><code>python /src/take_IpaH_features_coords.py -i /path_to/genome_assemblies/gbff/*.gbff</code></pre>

Output will be written to *IpaH_features.tsv* file

#### Extract all IpaH genes from reference *Shigella flexneri 2a str. 301* genome (assembly GCA_000006925.2)

<pre><code>python /src/extract_IpaH_features.py -i /path_to/genome_assemblies/gbff/GCA_000006925.2_ASM692v2_genomic.gbff -o /data/IpaH_from_reference.fa</code></pre>

This will produce output with 13 IpaH genes.

### 1. Run SibeliaZ

<pre><code>sibeliaz -a 2000 -k 15 -b 300 -n -t 2 -o /results/sibeliaz_out genome_assemblies/fna/merged.fna</code></pre>

### 2. Run Ragout

<pre><code>ragout-maf2synteny -s fine.txt -o /results/ragout_out/ /results/sibeliaz_out/blocks_coords.gff -b 1000</code></pre>

*fine.txt* file contains:

<pre><code>30 150
100 500
500 1500
</code></pre>

### 3. Convert Ragout output to .infercars format

<pre><code>python /src/ragout_to_infercars.py -dir /results/ragout_out/1000/</code></pre>

*blocks_coords.infercars* file will appear in *ragout_out* folder.

### 4. Convert .invercars to .tsv

<pre><code>python /src/parse_to_df.py -i /results/ragout_out/1000/blocks_coords.infercars -o /results/ragout_out/1000/blocks_coords.tsv</code></pre>

### 1. Covert infercars format to tsv

**parse_to_df.py** script converts _data/blocks_coords.infercars_, _data/blocks_coords_unique_gene.infercars_ files to _results/all_blocks.tsv_ and _results/common_blocks.tsv_ files correspondingly.

### 2. Check neighbour blocks of IpaH genes

**check_neighbours.py** script looks to all and common blocks for each lineage (_flexneri_, _boydii_, _dysenteriae_ and _sonnei_) from _results/all_blocks.tsv_ and _results/common_blocks.tsv_ files and found blocks, which are neighbours to IpaH genes blocks. The coordinates of neighbour blocks are in _results/lineage_name/*coords.tsv_, the frequencies of blocks are in _results/lineage_name/*frequencies.tsv_.

