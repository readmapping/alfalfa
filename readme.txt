ALFALFA v0.5
-------------------
  
Installing: 
	ALFALFA comes with 2 distributions: 1 compiled under Ubuntu Linux and the other under Cygwin-Windows
  
Usage:  
	./alfalfa  COMMAND [options] -x <reference-file> [-U <query-file>] [-1 <upstream mates> -2 <downstream mates>] [-S ouput-file]
  
Command should be one of the following: 
  index         only build the index for given <reference-file>, used for time calculations
  aln           map the reads from <query-file> to the index build for <reference-file>
  		if no output file is specified, output is written to <query-file>.sam


I/O OPTIONS
-----------
-x (string)                reference sequence in mult-fasta
-1                         query file with first mates (fasta or fastq)
-2                         query file with second mates (fasta or fastq)
-U                         query file with unpaired reads (fasta or fastq)
-S                         output file name (will be sam)

PERFORMANCE OPTIONS 
-------------------
-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index, 
			   value should be lower than -L parameter [1]
-p/--threads (int)    	   number of threads [1]

ALIGNMENT OPTIONS 
-----------------
-d/--errors (double)       percentage of errors allowed according to the edit distance [0.08]
-L/--seedminlength (int)   minimum length of the seeds used [depending on errorPercent and read length [20]
-k/--alignments (int)      expected number of alignments required per strand per read [100]
-T/--trials (int)          maximum number of times alignment is attempted before we give up [10]
-C/--mincoverage (int)     minimum percent of bases of read the seeds have to cover [25]
--tryharder                enable: 'try harder': when no seeds have been found, search using less stringent parameters
--seedthreads (int)        number of threads for calculating the seeds [1]
--nofw                     do not compute forward matches
--norc                     do not compute reverse complement matches
-n/--wildcards             treat Ns as wildcard characters
--softclipping             allow soft clipping at the beginning and end of an alignment

DYNAMIC PROGRAMMING OPTIONS 
---------------------------
-m/--match (int)           match bonus [0]
-u/--mismatch (int)        mismatch penalty [-2]
-o/--gapopen (int)         gap open penalty (set 0 for non-affine gap-penalties) [0]
-e/--gapextend (int)       gap extension penalty [-2]

PAIRED END OPTIONS 
------------------
-I/--minins (int)          minimum insert size [0]
-X/--maxins (int)          maximum insert sier [500]
--fr/--rf/--ff             orientation of the mates: fr means forward upstream mate 1 
                           and reverse complement downstream mate 2, or vise versa. rf and ff are similar [fr]
--no-mixed                 never search single mate alignments
--no-discordant            never discordant alignments
--dovetail                 allow reads to dovetail (changing up- and downstream of reads)
--no-contain               disallow a mate to be fully contained in the other
--no-overlap               disallow a mate to overlap with the other
--paired-mode (int)        choose algorithm to calculate paired-end reads

MISC OPTIONS 
------------
--verbose                  enable verbose mode (not by default)
-h/--help                  print this statement