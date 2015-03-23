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

0: make index dir:
mkdir ./idx/

1: run indexing for target genome:
  python swga_01_idx_all.py -f <FASTA.fa> -b <BLOCKSIZE (3000)> -p <PERCENT (100)> -I <INDEXDIR [./idx/]>
  python swga_01_idx_all.py -f pvivax_ref.fasta -b 3000

2: run indexing for subset of N percent of background genome(s)
  python swga_01_idx_all.py -f hsapiens_ref.fasta -b 3000 -p 2
  python swga_01_idx_all.py -f agambiae_ref.fasta -b 3000 -p 10

3: call primers from target genome (exhaustive)
### nb parallel code not currently working
  python swga_02_primers.py -t ~/refs/Pf3D7_v3.fasta  -b BACKGROUND1 [-b BACKGROUND2 ]  \
     -c <MAX_TM> -n <primer min length> -x <primer max length> \
	 -T <threads [1]>

4: take only those with highest genome counts (will roll into program later)
sort -k2gr candidates.txt | head -n 200 > top200candidates.txt

5: find those with highest whole-genome ratios
python swga_bwt_ratio.py top200candidates.tx

6: show across genome
#backgrounds currently hard-coded
python swga_bwt_cnt.py INDEX PATTERNS BACKGROUND1:BACKGROUND2

