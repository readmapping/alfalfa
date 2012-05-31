ALFALFA v0.3
-------------------
  
Installing: 
	ALFALFA comes with 2 distributions: 1 compiled under Ubuntu Linux and the other under Cygwin-Windows
  
Usage:  
	Linux	./alfalfa  COMMAND [options] <reference-file> <query-file> [ouput-file]
	Windows	alfalfa  COMMAND [options] <reference-file> <query-file> [ouput-file]
  
Command should be one of the following: 
  index         only build the index for given <reference-file>, used for time calculations
  aln           map the reads from <query-file> to the index build for <reference-file>
  		if no output file is specified, output is written to <query-file>.sam
  
PERFORMANCE OPTIONS (default options between square brackets)
  -s (int)         the sparsity factor of the sparse suffix array index, values between 1 and 4 are preferred [1]
  -q (int)         number of threads (untested) [1]
  
I/O OPTIONS 
  -v               enable verbose mode (not by default)

ALIGNMENT OPTIONS 
  -l (int)         minimum length of the seeds used [depending on errorPercent and read length, min 20]
  -k (int)         expected number of alignments required per strand per read [1]
  -T (int)         maximum number of times alignment is attempted before we give up [10]
  -C (int)         minimum percent of bases of read the seeds have to cover [25]
  -H (int)         enable: 'try harder': when no seeds have been found, search using less stringent parameters
  -t (int)         number of threads for calculating the seeds [1]
  -d (double)      percentage of errors allowed according to the edit distance [0.08]
  -f               do not compute forward matches
  -r               do not compute reverse complement matches
  -n               treat Ns as wildcard characters
  -c               allow soft clipping at the beginning and end of an alignment

DYNAMIC PROGRAMMING OPTIONS 
  -m (int)         match bonus [2]
  -u (int)         mismatch penalty [-1]
  -o (int)         gap open penalty (set 0 for non-affine gap-penalties) [-2]
  -e (int)         gap extension penalty [-1]