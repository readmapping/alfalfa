alfalfa
=======

ALFALFA: fast and accurate long read mapper

Rapid evolutions in sequencing technology force read mappers into flexible adaptation to longer reads, changing error models, memory barriers and novel applications. The long read mapper ALFALFA achieves high performance in accurately mapping long (>500bp) single-end and paired-end reads to gigabase-scale reference genomes, while remaining competitive for mapping shorter (>100bp) reads. Its seed-and-extend workflow is underpinned by fast retrieval of super-maximal exact matches from an enhanced sparse suffix array, with flexible parameter tuning to balance performance, memory footprint and accuracy.

# Downloads

[zip containing current version](https://github.com/readmapping/alfalfa/archive/master.zip)
[ALFALFA v0.8 zip](https://github.com/readmapping/alfalfa/archive/v0.8.zip)
[ALFALFA v0.8 tarball](https://github.com/readmapping/alfalfa/archive/v0.8.tar.gz)

# Installation

## Clone repository

1.  optional: install git
2.  type `git clone git://github.com/readmapping/alfalfa.git`
3.  change into the folder `alfalfa`
4.  type `make`

## Current version zip file

1. unpack using uncompression program such as [7-zip](www.7-zip.org/‎) or using `unzip alfalfa-master.zip`
2. change into the folder `alfalfa-master`
3.  type `make`

##specific version tarball  or zip-file

1. unpack
  * zip: unpack using uncompression program such as [7-zip](www.7-zip.org/‎) or using `unzip alfalfa-version.zip`
  * tar: `tar xvzf alfalfa-version.tar.gz`
2. change into the folder `alfalfa-version`
3. type `make`

# Usage

Usage and command line parameter information can be found on the [web page](http://alfalfa.ugent.be).