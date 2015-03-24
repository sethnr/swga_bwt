scripts for selecting and designing tags for selective WGA (swga)
counting swga tags across target genome, and using BWT for indexing of fastas
SWGA technique developed by Leichty / Brisson [http://www.genetics.org/content/early/2014/08/05/genetics.114.165498.abstract | genetics 2014]


requires numpy, hpy5

0: make index dir:
mkdir ./idx/

1.1: run indexing for target genome:
  python swga_01_idx_all.py -f <FASTA.fa> -b <BLOCKSIZE (3000)> -p <PERCENT (100)> -I <INDEXDIR [./idx/]>
  python swga_01_idx_all.py -f pvivax_ref.fasta -b 3000

1.2: run indexing for subset of N percent of background genome(s)
  python swga_01_idx_all.py -f hsapiens_ref.fasta -b 3000 -p 2
  python swga_01_idx_all.py -f agambiae_ref.fasta -b 3000 -p 10

2.1: call primers from target genome (exhaustive)
  python swga_02_primers.py -t ~/refs/Pf3D7_v3.fasta  -b BACKGROUND1 [-b BACKGROUND2 ]  \
     -c <MAX_TM> -n <primer min length> -x <primer max length> \
	 -T <threads [1]>
  #4 threads good, 8 better
  python swga_02_get_primers.py -t PlasmoDB-13.0_PvivaxSal1_Genome.CHR01.fasta \
  -b idx/Homo_sapiens_GRCh38.0.01.3000 -b idx/Pf3D7.0.005.3000 \
  -B 100 -T 8 -o candidates.txt

2.2: take only those with highest genome counts (will roll into program later)
  sort -k2gr candidates.txt | head -n 200 > top200candidates.txt

3: fill match index for primers
3.1: initialise match index
  python swga_03_build_match_idx.py -p candidates.txt \
    -t pvivax_ref.IDX.hdf5  -b hsapiens_ref.IDX.hsf5 -b agambiae_ref.IDX.hdf5 \
    -I -i matchIdx.hdf5
3.2: split candidates
    mkdir candidates
    split -n 100 candidates.txt -prefix candidates/c
3.3: run fill match index (on cluster)
  bsub 
  python swga_03_build_match_idx.py -p candidates.txt \
    -t pvivax_ref.IDX.hdf5  -b hsapiens_ref.IDX.hsf5 -b agambiae_ref.IDX.hdf5 \
    -A -i matchIdx.hdf5
   
4: chose primers with broadest representation in matcharray
  python swga_04_assess_plexes.py -p <primer_candidates> -i <match matrix (from previous step)> -c <no of primer combinations to assess [1000]> \
  -N <best N combinations 100> -r <min ratio target:background> -M <don't calculate max plex (exhaustive on large datasets)>
  -F forward weighting (bias of primer selection towards high-ratio samples) -c min no of counts in matrix -D decrease count and ratio for 1/D for D steps when chosing suboptimal
   
  python swga_04_find_plex.py -p top200candidates.txt -i matchIdx.hdf5 -c 1000 \
  -N 100 -k 20 -r 200 -c 100 -M 


5: show across genome
#fix later:
python swga_bwt_cnt.py INDEX PATTERNS BACKGROUND1:BACKGROUND2
swga_bwt_cnt = count occurrences of tags in index blocks
patterns will be read from first block of file, split on whitespace
> python swga_bwt_cnt.py index_name pattern_file

to quickly run all in index, do something like:
> cd index_dir
> ls *BWT* | perl -pe 's/\.BWT\.npy//gi' | xargs -Irepl python ../swga_bwt_cnt.py repl ../PATTERNS_FILE.txt ../out/repl.out


