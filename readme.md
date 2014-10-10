scripts for selecting and designing tags for selective WGA (swga)
counting swga tags across target genome, and using BWT for indexing of fastas
SWGA technique developed by Leichty / Brisson [http://www.genetics.org/content/early/2014/08/05/genetics.114.165498.abstract | genetics 2014]


requires numpy, hpy5

swga_bwt_idx = build index of chrs from fasta (currently uncompressed), split into blocks
> python swga_bwt_idx.py FASTA.fa block_size

swga_bwt_cnt = count occurrences of tags in index blocks
patterns will be read from first block of file, split on whitespace
> python swga_bwt_cnt.py index_name pattern_file

to quickly run all in index, do something like:
> cd index_dir
> ls *BWT* | perl -pe 's/\.BWT\.npy//gi' | xargs -Irepl python ../swga_bwt_cnt.py repl ../PATTERNS_FILE.txt ../out/repl.out
