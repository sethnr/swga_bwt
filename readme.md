scripts for counting swga tags across target genome
uses BWT for indexing of fastas

swga_bwt_idx = build index of chrs from fasta (currently uncompressed), split into blocks
> python swga_bwt_idx.py FASTA.fa block_size

swga_bwt_cnt = count occurrences of tags in index blocks
patterns will be read from first block of file, split on whitespace
> python swga_bwt_cnt.py index_name pattern_file
