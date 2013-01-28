ALFALFA v0.6.2
-------------------
  
Installing: 
	ALFALFA comes with 2 distributions: 1 compiled under Ubuntu Linux and the other under Cygwin-Windows
  
Usage:  
	./alfalfa  COMMAND [options] -x <reference-file> [-U <query-file>] [-1 <upstream mates> -2 <downstream mates>] [-S ouput-file]
  
Command should be one of the following: " << endl;
    index                       build the index for given <reference-file> and save to disk
                                this command is not necessairy for mapping, as aln can first construct the index
    aln                         map the reads to the index build for <reference-file>

call alfalfa COMMAND --help [or -h] for more detailed information

index COMMAND
-------------
-------------
./alfalfa  index [options] -x <reference-file>

OPTIONS
-------
-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index.
			   Note that the value needs to be lower than -L parameter in the ALN command[1].
-p/--prefix (string/path)  prefix of the index names [reference-file name]
--save (0 or 1)            save index to disk or not [1]"


aln COMMAND
-------------
-------------

I/O OPTIONS
-----------
-x (string)                reference sequence in mult-fasta
-i/--index (string)        prefix or path of the index to load, if not set, index will first be calculated
-1                         query file with first mates (fasta or fastq)
-2                         query file with second mates (fasta or fastq)
-U                         query file with unpaired reads (fasta or fastq)
-S                         output file name (will be sam) [referenceName.sam]
--save (0 or 1)            if --index is not set, save index to disk [0]
-p                         prefix of index that will be saved [reference sequence name]

PERFORMANCE OPTIONS 
-------------------
-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index if it is not yet constructed [1].
-q/--threads (int)    	   number of threads [1]

ALIGNMENT OPTIONS 
-----------------
-d/--errors (double)       percentage of errors allowed according to the edit distance [0.08]
-L/--seedminlength (int)   minimum length of the seeds used [depending on errorPercent and read length [20]
-k/--alignments (int)      expected number of alignments required per strand per read [100]
-T/--trials (int)          maximum number of times alignment is attempted before we give up [10]
-C/--mincoverage (int)     minimum percent of bases of read the seeds have to cover [25]
--tryharder                enable: 'try harder': when no seeds have been found, search using less stringent parameters
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
--paired-mode (int)        choose algorithm to calculate paired-end reads 1,2,3,4 [1]

MISC OPTIONS 
------------
--verbose                  enable verbose mode (not by default)
-h/--help                  print this statement