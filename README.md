alfalfa
=======

ALFALFA: fast and accurate long read mapper

url: https://github.ugent.be/pages/ComputationalBiology/alfalfa/

Rapid evolutions in sequencing technology force read mappers into flexible adaptation to longer reads, changing error models, memory barriers and novel applications. The long read mapper ALFALFA achieves high performance in accurately mapping long (>500bp) single-end and paired-end reads to gigabase-scale reference genomes, while remaining competitive for mapping shorter (>100bp) reads. Its seed-and-extend workflow is underpinned by fast retrieval of super-maximal exact matches from an enhanced sparse suffix array, with flexible parameter tuning to balance performance, memory footprint and accuracy.

# Downloads

[zip containing current version](https://github.ugent.be/ComputationalBiology/alfalfa/archive/master.zip)
[ALFALFA v0.7.2 zip](https://github.ugent.be/ComputationalBiology/alfalfa/archive/v0.7.2.zip)
[ALFALFA v0.7.2 tarball](https://github.ugent.be/ComputationalBiology/alfalfa/archive/v0.7.2.tar.gz)

# Installation

## Clone repository

1.  optional: install git
2.  type `git clone git://github.ugent.be/ComputationalBiology/alfalfa.git`
3.  change into the folder `alfalfa`
4.  type `make`

## Current version zip file

1. unpack using uncompression program such as [7-zip](www.7-zip.org/‎) or using `unzip alfalfa-master.zip`
2. change into the folder `alfalfa-master`
3. type `make`

##specific version tarball  or zip-file

1. unpack
  * zip: unpack using uncompression program such as [7-zip](www.7-zip.org/‎) or using `unzip alfalfa-version.zip`
  * tar: `tar xvzf alfalfa-version.tar.gz`
2. change into the folder `alfalfa-version`
3. type `make`

# Usage

Usage and command line parameter information can be found on the [web page](https://github.ugent.be/pages/ComputationalBiology/alfalfa/).

# Version history

* Version 0.7.2 (June 13, 2013)  
  1. Improved online help
  2. Reordered source code in folders and updated makefile
* Version 0.7.1 (May 23, 2013)  
  1. Changed parameter names and defaults.  
  2. Removed deprecated code.  
  3. Improved paired-end mode 4.  
* Version 0.7.0 (Apr 19, 2013)  
  1. General improvements to the code.  
  2. Small bug fixes.  
  3. Added more debug information to verbose mode.   
  4. Evaluate command can now calculate edit distance of alignments using reference genome and CIGAR string.  
  5. Added local alignment option.  
  6. Improved SMEM calculation.  
* Version 0.6.4 (Mar 04, 2013)  
  1. Introduction of essaMEM features.  
  2. Introduction of a k-mer array to speed up seed-finding.  
  3. Bug fixes.  
  4. Fixed some typos.  
  5. General improvements to the code and removal of deprecated code.  
* Version 0.6.3 (Feb 21, 2013)  
  1. General improvements to the code.  
  2. Added verbose mode for debugging purposes.  
* Version 0.6.2 (Jan 23, 2013)  
  1. Fixed bugs introduced by previous version.  
  2. Added wgsim subcommand to evaluate command.  
* Version 0.6.1 (Jan 14, 2013)  
  1. Artifacts indicating mate origin at the end of paired-end reads (/0/1/2)  are now removed from read name.  
  2. Bug fixes in evaluate command.  
* Version 0.6.0 (Nov 16, 2012)  
  1. Added option to store index to disk.  
  2. Bug fixes in evaluate command.  
* Version 0.6.0 (Nov 16, 2012)  
  1. Added option to store index to disk.  
  2. Bug fixes in evaluate command.  
* Version 0.5.4 (Oct 09, 2012)  
  1. Added evaluate command to check accuracy of mapping on simulated reads.  
* Version 0.5.3 (Sep 27, 2012)  
  1. Bug fixes for paired-end mapping.  
  2. Memory leak fixed that was introduced in previous version.  
  3. Memory usage decreased for index.  
* Version 0.5.2 (Sep 24, 2012)  
  1. Bug fixes for paired-end mapping.  
  2. Decreased memory usage of dynamic programming matrix.  
* Version 0.5.1 (Sep 13, 2012)  
  1. Bug fixes for paired-end mapping.  
  2. Bug fixes for reference genomes larger than 2Gbp.  
* Version 0.5.0 (Sep 03, 2012)  
  1. Major reordering of source code. 
  2. Many bug fixes.  
  3. Added new behavior for single-end read alignment that finds seeds for both strands first, followed by prioritzed extension of candidate regions, instead of aligning both srands separately.  
  4. Changed command line parameters and default parameter values.  
  5. Added support for paired-end reads.  
  6. Added many parameters for paired-end read mapping.  
  7. Added different algorithms for paired-end read mapping.  
  8. Added banded dynamic programming with non-symmetrical band.  
* Version 0.3.5 (Jul 03, 2012)  
  1. Dynamic programming matrix is now semi-static, decreasing computation time.  
  2. Trace calculation of dynamic programming without trace matrix.  
* Version 0.3.4 (Jun 29, 2012)  
  1. Changed from normal dynamic programming to banded dynamic programming.  
  2. Score matrix added for dynamic programming (internal).  
* Version 0.3.3 (Jun 26, 2012)  
  1. Enhanced sparse suffix array with sparse child array.  
* Version 0.3.2 (Jun 08, 2012)  
  1. Multithreading support added.  
  2. Minor improvements in the code.  
* Version 0.3.1 (May 31, 2012)  
  1. Optimized dynamic programming when no affine gap penalties are used.  
  2. Changes to default parameters.  
* Version 0.3.0 (May 31, 2012)  
  1. First Github release.  
* Version 0.2.0 (March 19, 2012)  
  1. Allow multiple alignments to be returned per read. Added parameter to set the maximum number of alignments returned per read.  
  2. Added algorithm to select multiple candidate regions.  
  3. Added some heuristics to speed up mapping time.  
* Version 0.1.0 (March 15, 2012)  
  1. Changed MAM seeds to SMEM seeds.  
  2. Added parameter to limit number of seeds per query offset.
  3. Bugfixes for local alignment, mapping quality calculation and CIGAR string calculation.  
  4. Dynamic programming at boundaries now takes into account current edit distance of alignment in the gaps between seeds.  
* Version 0.0.0 (March 01, 2012)  
  1. First release of ALFALFA.  

# Links  

[essaMEM project](https://github.ugent.be/ComputationalBiology/essaMEM)

# Contact  

Don't hesitate to contact us at michael.vyverman[at]ugent.be if you have any further questions or suggestions. 
