### 1. Covert infercars format to tsv

**parse_to_df.py** script converts _data/blocks_coords.infercars_, _data/blocks_coords_unique_gene.infercars_ files to _results/all_blocks.tsv_ and _results/common_blocks.tsv_ files correspondingly.

### 2. Check neighbour blocks of IpaH genes

**check_neighbours.py** script looks to all and common blocks for each lineage (_flexneri_, _boydii_, _dysenteriae_ and _sonnei_) from _results/all_blocks.tsv_ and _results/common_blocks.tsv_ files and found blocks, which are neighbours to IpaH genes blocks. The coordinates of neighbour blocks are in _results/lineage_name/*coords.tsv_, the frequencies of blocks are in _results/lineage_name/*frequencies.tsv_.
