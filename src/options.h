/*
 * Copyright 2012, Michael Vyverman <michael.vyverman@ugent.be>
 *
 * This file is part of ALFALFA.
 *
 * ALFALFA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ALFALFA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ALFALFA.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPTIONS_H
#define	OPTIONS_H

#include <string.h>
#include <stdint.h>

#include "optionparser.h"

enum mum_t { MUM, MAM, MEM, SMEM, SMAM, SMEMP };
static const int MINOPTIONCOUNT = 2;
enum command_t {INDEX, MATCHES, ALN};
enum orientation_t { PAIR_FR, PAIR_RF, PAIR_FF};

//match options
struct align_opt {
    int match;
    int mismatch;
    int openGap;
    int extendGap;
    int8_t mat[25];
    bool scoresSetByUser;
    bool global;
    double errorPercent;
    int minMemLength;
    bool fixedMinLength;
    int maxSeedCandidates;
    int maxSmemMemStart;
    int minMemRegLength;
    int maxRegionMems;
    int alignmentCount;
    int maxTrial;
    int minTrial;//
    bool maxTrialBest;
    int maxAlnCount;
    double minCoverage;//deprecated
    bool alwaysUnique;
    bool rescue;
    bool noSparseQ;
    bool fullDP;
    bool noFW; 
    bool noRC;
    //print Values: 0 = None, 1 = Basic, 2 = performance, 3 = read + pairing in PE, 4 = alignment
    // 5 = region, 6 = mems, 7 = DP info
    int print;
    FILE * debugFile;
    mum_t memType;
    int sparseMult;
    //next parameters are currently only used for local alignment
    int fixedBandSize;
    
    void fill_m();
    
};

struct paired_opt {
    int minInsert;
    int maxInsert;
    bool mixed;
    bool discordant;
    bool dovetail;
    bool contain;
    bool overlap;
    int mode;
    bool pairedRescue;
    orientation_t orientation;
};

struct mapOptions_t{//commentary + sort + constructor
    mapOptions_t(){ initOptions(); }
    void initOptions();
    void printOptions();
    
    //command option
    command_t command;
    //performance options
    int K;//sparsity factor
    int query_threads;
    //I/O options
    std::string unpairedQ;
    std::string pair1;
    std::string pair2;
    std::string ref_fasta;
    std::string index_prefix;
    bool saveIndex;
    bool hasChild;
    bool hasSuflink;
    bool hasKmer;
    std::string indexLocation;
    std::string outputName;
    //Sequence options
    bool nucleotidesOnly;
    //alignment options
    align_opt alnOptions;    
    //paired end options
    paired_opt pairedOpt;
};

struct Arg: public option::Arg
{
  static option::ArgStatus Required(const option::Option& option, bool)
  {
    return option.arg == 0 ? option::ARG_ILLEGAL : option::ARG_OK;
  }
  static option::ArgStatus Empty(const option::Option& option, bool)
  {
    return (option.arg == 0 || option.arg[0] == 0) ? option::ARG_OK : option::ARG_IGNORE;
  }
};

//Reads from - to ; trim ; phred quals
enum optionIndex {
    ARG_UNKNOWN,          //not found
    ARG_REFERENCE,
    ARG_SPARSENESS,
    ARG_PREFIX,
    ARG_NO_CHILD,
    ARG_SUFLINK,
    ARG_NO_KMER,
    ARG_INDEX,
    ARG_SAVE,
    ARG_SINGLE,
    ARG_MATES1,
    ARG_MATES2,
    ARG_OUTPUT,
    ARG_ALIGNMENTS,
    ARG_NO_FORWARD,
    ARG_NO_REVERSE,
    ARG_EDIT_DISTANCE,
    ARG_THREADS,
    ARG_SEED,
    ARG_MIN_LENGTH,
    ARG_MAX_SEEDS,
    ARG_NO_RESCUE,
    ARG_NO_SPARSEQ,
    ARG_FULLDP,
    ARG_MAX_SMEM,
    ARG_MIN_MEML,
    ARG_MAX_MEM,
    ARG_MAX_FAILURES,
    ARG_MAX_F_BEST,
    ARG_MIN_ALIGN,
    ARG_MIN_COVERAGE,
    ARG_LOCAL,
    ARG_ALWAYS_UNIQUE,
    ARG_BANDWIDTH,
    ARG_MATCH,
    ARG_MISMATCH,
    ARG_GAP_OPEN,
    ARG_GAP_EXTEND,
    ARG_MIN_INSERT,
    ARG_MAX_INSERT,
    ARG_ORIENTATION,     //--orientation
    ARG_NO_MIXED,        //--no-mixed
    ARG_NO_DISCORDANT,   //--no-discordant
    ARG_DOVETAIL,        //--dovetail
    ARG_NO_CONTAIN,      //--no-contain
    ARG_NO_OVERLAP,      //--no-overlap
    ARG_PAIRED_MODE,       //--paired-mode
    ARG_PAIRED_RESCUE,
    ARG_VERBOSE,
    ARG_DEBUG,
    ARG_HELP
};

const option::Descriptor alignUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          0},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nUsage: alfalfa align [option...]\nThe align command is used for mapping and aligning a read set onto a reference genome. As this process can be customized through a long list of options, we have grouped them into several categories."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nI/O options"},
    {ARG_REFERENCE,     0, "r", "reference",     Arg::Required,      "-r/--reference \t(file, automatic in index).\vSpecifies the location of a file that contains the reference genome in multi-fasta format."},
    {ARG_INDEX,         0, "i", "index",         Arg::Required,      "-i/--index \t(string).\vSpecifies the prefix used to name all generated index files. If this option is not set explicitly, an index will be computed from the reference genome according to the settings of the options that also apply to the index command."},
    {ARG_SAVE,          0, "",  "save",          Arg::None,          "--save \t.\vSpecifies that if an index is constructed by the align command itself, it will be stored to disk. This option is ignored if the index is loaded from disk (option -i)."},
    {ARG_SINGLE,        0, "0", "single",        Arg::Required,      "-0/--single \t(file).\vSpecifies the location of a file that contains single-end reads. Both fasta and fastQ formats are accepted. If both single-end and paired-end reads are specified, single-end reads are processed first."},
    {ARG_MATES1,        0, "1", "mates1",        Arg::Required,      "-1/--mates1 \t(file).\vSpecifies the location of a file that contains the first mates of paired-end reads. Both fasta and fastQ formats are accepted."},
    {ARG_MATES2,        0, "2", "mates2",        Arg::Required,      "-2/--mates2 \t(file).\vSpecifies the location of a file that contains the second mates of paired-end reads. Both fasta and fastQ formats are accepted."},
    {ARG_OUTPUT,        0, "o", "output",        Arg::Required,      "-o/--output \t(file, filename passed to the -r option with additional .sam extension).\vSpecifies the location of the generated SAM output file containing the results of read mapping and alignment."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nAlignment options"},
    {ARG_ALIGNMENTS,    0, "a", "alignments",    Arg::Required,      "-a/--alignments \t(int, 1).\vSpecifies the maximum number of alignments reported per read."},
    {ARG_NO_FORWARD,    0, "",  "no-forward",    Arg::None,          "--no-forward \t.\vDo not compute alignments on the forward strand."},
    {ARG_NO_REVERSE,    0, "",  "no-reverse",    Arg::None,          "--no-reverse \t.\vDo not compute alignments on the reverse complement strand."},
    {ARG_EDIT_DISTANCE, 0, "e", "edit-distance", Arg::Required,      "-e/--edit-distance \t(float, 0.08).\vRepresents the maximum percentage of differences allowed in accepting alignments and used in combination with the dynamic programming score function to calculate the minimum alignment score."},
    {ARG_NO_RESCUE,     0, "",  "no-rescue",     Arg::None,          "--no-rescue \t.\vDisables rescue procedures that are normally initiated when ALFALFA does not find seeds and/or alignments with the current parameters."},
    {ARG_THREADS,       0, "t", "threads",       Arg::Required,      "-t/--threads \t(int, 1).\vNumber of threads used during read mapping. Using more than one thread results in reporting read alignments in a different order compared to the order in which they are read from the input file(s)."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nSeed options"},
    {ARG_SEED,          0, "",  "seed",          Arg::Required,      "--seed \t(MEM | SMEM | PSMEM, SMEM).\vSpecifies the type of seeds used for read mapping. Possible values are MEM for maximal exact matches, SMEM for super-maximal exact matches, and PSMEM for a combination of both. The use of SMEMs generally boosts performance without having a negative impact on accuracy compared to the use of MEMs. On the other hand, there are usually many more MEMs than SMEMs, in general resulting in a higher number of candidate genomic regions. The latter might be useful if reporting more candidate mapping locations is preferred."},
    {ARG_MIN_LENGTH,    0, "l", "min-length",    Arg::Required,      "-l/--min-length \t(int, auto).\vSpecifies the minimum seed length. This value must be greater than the sparseness value used to build the index (option -s). By default, the value of this option is computed automatically using the following procedure. A value of 40 is used for reads shorter than 1kbp. The value is incremented by 20 for every 500bp above 1kbp, with the total increment being divided by the maximum percentage of errors allowed in accepting alignments (option -e)."},
    {ARG_MAX_SEEDS,     0, "m", "max-seeds",     Arg::Required,      "-m/--max-seeds \t(int, 10000).\vSpecifies the maximum number of seeds per length and offset in the read sequence. The value passed to this option is multiplied by the automatically computed skip factor that determines sparse matching of sampled suffixes from the read sequence. As a result, the actual number of seeds per starting position in the read might still vary. Higher values of this option result in higher numbers of seeds, increasing in turn the number of candidate genomic regions."},
    {ARG_MAX_SMEM,      0, "",  "max-smems",     Arg::Required,      "--max-smems \t(int, 10).\vSpecifies the maximum number of SMEMs per offset in the read sequence to allow MEM-finding. This only applies for PSMEM seeds."},
    {ARG_MAX_MEM,       0, "",  "max-mems",      Arg::Required,      "--max-mems \t(int, 20).\vSpecifies the maximum number of MEMs per offset in the read that can be used for candidate region identification."},
    {ARG_MIN_MEML,      0, "",  "min-mem-length",Arg::Required,      "--min-mem-length \t(int, 50).\vSpecifies minimum length need to have to be used for candidate region identification."},
    {ARG_NO_SPARSEQ,    0, "",  "no-sparseness", Arg::None,          "--no-sparseness \t.\vDisables the use of sparseness in the read sequence during seed-finding. This option is increases runtime and can have a slight effect on accuracy."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nExtend options"},
    {ARG_MAX_FAILURES,  0, "f", "max-failures",  Arg::Required,      "-f/--max-failures \t(int, 10).\vSpecifies the maximum number of successive candidate regions that are investigated without success before ALFALFA stops extending the candidate regions of a read. Extension can be restarted only if the remaining candidate regions contain unique seeds."},
    {ARG_MAX_F_BEST,    0, "",  "reset--failures",Arg::None,          "--reset-failures \t.\v If set, the counter of successive candidate regions that are investigated without success is reset if a feasible alignment is found. By default the counter is only reset if a new best alignment is found."},
    {ARG_MIN_ALIGN,     0, "",  "max-alignments",Arg::Required,      "--max-alignments \t(int, 5000).\vSpecifies the maximum number of alignments calculated per read. This value should be higher than the number of reported alignments (option -a). Decreasing this value can increase performance of the algorithm, at the cost of a lower accuracy and worse mapping quality estimation."},
    {ARG_MIN_COVERAGE,  0, "c", "min-coverage",  Arg::Required,      "-c/--min-coverage \t(float, 0.25).\vSpecifies the minimum percentage of the read length that candidate region containing a single seed need to cover before extension of the candidate region is taken into consideration."},
    {ARG_ALWAYS_UNIQUE, 0, "",  "skip-unique",   Arg::Required,      "--skip-unique \t.\vBy default, ALFALFA extends all candidate regions containing unique seeds. If this flag is set, this criterium is not taken into account when deciding to extend a candidate region."},
    {ARG_FULLDP,        0, "",  "full-dp",       Arg::None,          "--full-dp \t.\vBy default, ALFALFA uses a chain-guided alignment to retrieve the CIGAR alignment. If this parameter is set, banded dynamic programming is performed instead. The use of this parameter can greatly increase runtime, but can lead to more optimal alignments in some cases."},
    {ARG_LOCAL,         0, "",  "local",         Arg::None,          "--local \t.\vBy default, ALFALFA uses global alignment during the last phase of the mapping process. Global alignment in essence is end-to-end alignment, as it entirely covers the read but only covers the reference genome in part. Local alignment is used during the last phase of the mapping process if the --local option is set, which may result in soft clipping of the read."},
    {ARG_BANDWIDTH,     0, "b", "bandwidth",     Arg::Required,      "-b/--bandwidth \t(int, 100).\vSpecifies the maximum bandwidth that is used by the banded alignment algorithm. The bandwidth used is automatically inferred from the specification of the maximum percentage of errors allowed in accepting alignments (option -e), but is bounded by this parameter."},
    {ARG_MATCH,         0, "M", "match",         Arg::Required,      "-M/--match \t(int, 1).\vSpecifies the positive score assigned to matches in the dynamic programming extension phase."},
    {ARG_MISMATCH,      0, "U", "mismatch",      Arg::Required,      "-U/--mismatch \t(int, -4).\vSpecifies the penalty assigned to mismatches in the dynamic programming extension phase."},
    {ARG_GAP_OPEN,      0, "O", "gap-open",      Arg::Required,      "-O/--gap-open \t(int, -6).\vSpecifies the penalty O for opening a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E. The use of affine gap penalties can be disabled by setting this value to zero."},
    {ARG_GAP_EXTEND,    0, "E", "gap-extend",    Arg::Required,      "-E/--gap-extend \t(int, -1).\vSpecifies the penalty E for extending a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nPaired-end mapping options"},
    {ARG_MIN_INSERT,    0, "I", "min-insert",    Arg::Required,      "-I/--min-insert \t(int, 0).\vSpecifies the minimum insert size."},
    {ARG_MAX_INSERT,    0, "X", "max-insert",    Arg::Required,      "-X/--max-insert \t(int, 1000).\vSpecifies the maximum insert size."},
    {ARG_ORIENTATION,   0, "",  "orientation",   Arg::Required,      "--orientation \t(fr | rf | ff, fr).\vSpecifies the orientation of mates. fr means a forward upstream first mate and reverse complemented downstream second mate or vice versa. rf means a reverse complemented upstream first mate and forward downstream second mate or vice versa. ff means a forward upstream first mate and forward downstream second mate or vice versa. Note that these definitions are literally taken over from Bowtie 2."},
    {ARG_NO_MIXED,      0, "",  "no-mixed",      Arg::None,          "--no-mixed \t.\vDisables searching for unpaired alignments."},
    {ARG_NO_DISCORDANT, 0, "",  "no-discordant", Arg::None,          "--no-discordant \t.\vDisables searching for discordant alignments."},
    {ARG_DOVETAIL,      0, "",  "dovetail",      Arg::None,          "--dovetail \t.\vAllows switching between upstream and downstream mates in the definition of their orientation (option --orientation)."},
    {ARG_NO_CONTAIN,    0, "",  "no-contain",    Arg::None,          "--no-contain \t.\vDisallows concordant mates to be fully contained within each other."},
    {ARG_NO_OVERLAP,    0, "",  "no-overlap",    Arg::None,          "--no-overlap \t.\vDisallows concordant mates to overlap each other."},
    {ARG_PAIRED_MODE,   0, "",  "paired-mode",   Arg::Required,      "--paired-mode \t(1 | 2 | 3 | 4 | 5 | 6, 1).\vSpecifies the algorithm used to align paired-end reads. The possible algorithms are discussed in detail in the methods section. Algorithms 1 and 2 do not use information from candidate regions. Algorithms 3 and 4 prioritize extension of candidate regions over both reads. Algorithms 5 and 6 filter the list of candidate regions using the paired-end restraints. Algorithms with an odd number pair mapped reads afterward alignment. Algorithms with an even number perform dynamic programming across a window defined by the insert size restrictions to search for a bridging alignment reaching the other mate."},
    {ARG_PAIRED_RESCUE, 0, "",  "paired-rescue", Arg::None,          "--paired-rescue \t.\vEnable an automatical rescue procedue if no concordant alignment was found using the current parameter settings."},
    {0,0,0,0,0,0}
};

const option::Descriptor basicUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "Usage: alfalfa <command> [<subcommand>] [option...]\nCommand should be index, align or evaluate\nSubcommand is only required for the evaluate command\n\ncommands:\n"},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "index \tis used to construct the data structures for indexing a given reference genome."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "align \tis used for mapping and aligning a read set onto a reference genome."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "evaluate \tis used for evaluating the accuracy of simulated reads and summarizing statistics from the SAM-formatted alignments reported by a read mapper."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\ncall alfalfa <command> -h/--help for more detailed information on the specific commands\n"},
    {0,0,0,0,0,0}
};

const option::Descriptor indexUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "Usage: alfalfa index [option...]\nindex is used to construct the data structures for indexing a given reference genome.\n\noptions \n"},
    {ARG_REFERENCE,     0, "r", "reference",     Arg::Required,      "-r/--reference \t(file).\vSpecifies the location of a file that contains the reference genome in multi-fasta format."},
    {ARG_SPARSENESS,    0, "s", "sparseness",    Arg::Required,      "-s/--sparseness \t(int, 12).\vSpecifies the sparseness of the index structure as a way to control part of the speed-memory trade-off."},
    {ARG_PREFIX,        0, "p", "prefix",        Arg::Required,      "-p/--prefix \t(string, filename passed to the -r option).\vSpecifies the prefix that will be used to name all generated index files. The same prefix has to be passed to the -i option of the align command to load the index structure when mapping reads."},
    {ARG_NO_CHILD,      0, "",  "no-child",      Arg::None,          "--no-child \t.\vBy default, a sparse child array is constructed and stored in an index file with extension .child. The construction of this sparse child array is skipped when the --no-child option is set. This data structure speeds up seed-finding at the cost of (4/s) bytes per base in the reference genome. As the data structure provides a major speed-up, it is advised to have it constructed."},
    {ARG_SUFLINK,       0, "",  "suflink",       Arg::None,          "--suflink \t.\vSuffix link support is disabled by default. Suffix link support is enabled when the --suflink option is set, resulting in an index file with extension .isa to be generated. This data structure speeds up seed-finding at the cost of (4/s) bytes per base. It is only useful when sparseness is less than four and minimum seed length is very low (less than 10), because it conflicts with skipping suffixes in matching the read. In practice, this is rarely the case."},
    {ARG_NO_KMER,       0, "",  "no-kmer",       Arg::None,          "--no-kmer \t.\vBy default, a 10-mer lookup table is constructed that contains the suffix array interval positions to depth 10 in the virtual suffix tree. It is stored in an index file with extension .kmer and required only 8MB of memory. The construction of this lookup table is skipped when the --no-kmer option is set. The lookup table stores intervals for sequences of length 10 that only contain {A,C,G,T}. This data structure speeds up seed-finding if the minimum seed length is greater than 10."},
    {0,0,0,0,0,0}
};

const option::Descriptor miscUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nMiscellaneous  options"},
    {ARG_VERBOSE,       0, "v", "verbose",       Arg::Required,      "-v/--verbose \t(int, 0).\vTurns on lots of progress reporting about the alignment process. Higher numbers give more verbose output. Information is printed to standard error and is useful for debugging purposes. The default value 0 disables progress reporting. The maximum verbosity level is 7."},
    {ARG_HELP,          0, "h", "help",          Arg::None,          "-h/--help \t.\vPrints to standard error the version number, usage description and an overview of the options that can be used to customize the software package."},
    {ARG_DEBUG,         0, "",  "debug",         Arg::Required,      "--debug \t(file).\vSpecifies file to print specific debug information to that can be further examined. This info is structured to allow usage by other (inhouse) debugging tools and is currently experimental."},
    {0,0,0,0,0,0}
};

const option::Descriptor fullUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          0},
    {ARG_REFERENCE,     0, "r", "reference",     Arg::Required,      "-r/--reference \t(file, automatic from index).\vSpecifies the location of a file that contains the reference genome in multi-fasta format."},
    {ARG_SPARSENESS,    0, "s", "sparseness",    Arg::Required,      "-s/--sparseness \t(int, 12).\vSpecifies the sparseness of the index structure as a way to control part of the speed-memory trade-off."},
    {ARG_PREFIX,        0, "p", "prefix",        Arg::Required,      "-p/--prefix \t(string, filename passed to the -r option).\vSpecifies the prefix that will be used to name all generated index files. The same prefix has to be passed to the -i option of the align command to load the index structure when mapping reads."},
    {ARG_NO_CHILD,      0, "",  "no-child",      Arg::None,          "--no-child \t.\vBy default, a sparse child array is constructed and stored in an index file with extension .child. The construction of this sparse child array is skipped when the --no-child option is set. This data structure speeds up seed-finding at the cost of (4/s) bytes per base in the reference genome. As the data structure provides a major speed-up, it is advised to have it constructed."},
    {ARG_SUFLINK,       0, "",  "suflink",       Arg::None,          "--suflink \t.\vSuffix link support is disabled by default. Suffix link support is enabled when the --suflink option is set, resulting in an index file with extension .isa to be generated. This data structure speeds up seed-finding at the cost of (4/s) bytes per base. It is only useful when sparseness is less than four and minimum seed length is very low (less than 10), because it conflicts with skipping suffixes in matching the read. In practice, this is rarely the case."},
    {ARG_NO_KMER,       0, "",  "no-kmer",       Arg::None,          "--no-kmer \t.\vBy default, a 10-mer lookup table is constructed that contains the suffix array interval positions to depth 10 in the virtual suffix tree. It is stored in an index file with extension .kmer and required only 8MB of memory. The construction of this lookup table is skipped when the --no-kmer option is set. The lookup table stores intervals for sequences of length 10 that only contain {A,C,G,T}. This data structure speeds up seed-finding if the minimum seed length is greater than 10."},
    {ARG_REFERENCE,     0, "r", "reference",     Arg::Required,      "-r/--reference \t(file).\vSpecifies the location of a file that contains the reference genome in multi-fasta format."},
    {ARG_INDEX,         0, "i", "index",         Arg::Required,      "-i/--index \t(string).\vSpecifies the prefix used to name all generated index files. If this option is not set explicitly, an index will be computed from the reference genome according to the settings of the options that also apply to the index command."},
    {ARG_SAVE,          0, "",  "save",          Arg::None,          "--save \t.\vSpecifies that if an index is constructed by the align command itself, it will be stored to disk. This option is ignored if the index is loaded from disk (option -i)."},
    {ARG_SINGLE,        0, "0", "single",        Arg::Required,      "-0/--single \t(file).\vSpecifies the location of a file that contains single-end reads. Both fasta and fastQ formats are accepted. If both single-end and paired-end reads are specified, single-end reads are processed first."},
    {ARG_MATES1,        0, "1", "mates1",        Arg::Required,      "-1/--mates1 \t(file).\vSpecifies the location of a file that contains the first mates of paired-end reads. Both fasta and fastQ formats are accepted."},
    {ARG_MATES2,        0, "2", "mates2",        Arg::Required,      "-2/--mates2 \t(file).\vSpecifies the location of a file that contains the second mates of paired-end reads. Both fasta and fastQ formats are accepted."},
    {ARG_OUTPUT,        0, "o", "output",        Arg::Required,      "-o/--output \t(file, filename passed to the -r option with additional .sam extension).\vSpecifies the location of the generated SAM output file containing the results of read mapping and alignment."},
    {ARG_ALIGNMENTS,    0, "a", "alignments",    Arg::Required,      "-a/--alignments \t(int, 1).\vSpecifies the maximum number of alignments reported per read."},
    {ARG_NO_FORWARD,    0, "",  "no-forward",    Arg::None,          "--no-forward \t.\vDo not compute alignments on the forward strand."},
    {ARG_NO_REVERSE,    0, "",  "no-reverse",    Arg::None,          "--no-reverse \t.\vDo not compute alignments on the reverse complement strand."},
    {ARG_EDIT_DISTANCE, 0, "e", "edit-distance", Arg::Required,      "-e/--edit-distance \t(float, 0.08).\vRepresents the maximum percentage of differences allowed in accepting alignments and used in combination with the dynamic programming score function to calculate the minimum alignment score."},
    {ARG_THREADS,       0, "t", "threads",       Arg::Required,      "-t/--threads \t(int, 1).\vNumber of threads used during read mapping. Using more than one thread results in reporting read alignments in a different order compared to the order in which they are read from the input file(s)."},
    {ARG_SEED,          0, "",  "seed",          Arg::Required,      "--seed \t(MEM | SMEM | PSMEM, SMEM).\vSpecifies the type of seeds used for read mapping. Possible values are MEM for maximal exact matches, SMEM for super-maximal exact matches, and PSMEM for a combination of both. The use of SMEMs generally boosts performance without having a negative impact on accuracy compared to the use of MEMs. On the other hand, there are usually many more MEMs than SMEMs, in general resulting in a higher number of candidate genomic regions. The latter might be useful if reporting more candidate mapping locations is preferred."},
    {ARG_MIN_LENGTH,    0, "l", "min-length",    Arg::Required,      "-l/--min-length \t(int, auto).\vSpecifies the minimum seed length. This value must be greater than the sparseness value used to build the index (option -s). By default, the value of this option is computed automatically using the following procedure. A value of 40 is used for reads shorter than 1kbp. The value is incremented by 20 for every 500bp above 1kbp, with the total increment being divided by the maximum percentage of errors allowed in accepting alignments (option -e)."},
    {ARG_MAX_SEEDS,     0, "m", "max-seeds",     Arg::Required,      "-m/--max-seeds \t(int, 10000).\vSpecifies the maximum number of seeds per length and offset in the read sequence. The value passed to this option is multiplied by the automatically computed skip factor that determines sparse matching of sampled suffixes from the read sequence. As a result, the actual number of seeds per starting position in the read might still vary. Higher values of this option result in higher numbers of seeds, increasing in turn the number of candidate genomic regions."},
    {ARG_MAX_SMEM,      0, "",  "max-smems",     Arg::Required,      "--max-smems \t(int, 10).\vSpecifies the maximum number of SMEMs per offset in the read sequence to allow MEM-finding. This only applies for PSMEM seeds."},
    {ARG_MAX_MEM,       0, "",  "max-mems",      Arg::Required,      "--max-mems \t(int, 20).\vSpecifies the maximum number of MEMs per offset in the read that can be used for candidate region identification."},
    {ARG_MIN_MEML,      0, "",  "min-mem-length",Arg::Required,      "--min-mem-length \t(int, 50).\vSpecifies minimum length need to have to be used for candidate region identification."},
    {ARG_NO_RESCUE,     0, "",  "no-rescue",     Arg::None,          "--no-rescue \t.\vDisables rescue procedures that are normally initiated when ALFALFA does not find seeds and/or alignments with the current parameters."},
    {ARG_NO_SPARSEQ,    0, "",  "no-sparseness", Arg::None,          "--no-sparseness \t.\vDisables the use of sparseness in the read sequence during seed-finding. This option is increases runtime and can have a slight effect on accuracy."},
    {ARG_FULLDP,        0, "",  "full-dp",       Arg::None,          "--full-dp \t.\vBy default, ALFALFA uses a chain-guided alignment to retrieve the CIGAR alignment. If this parameter is set, banded dynamic programming is performed instead. The use of this parameter can greatly increase runtime, but can lead to more optimal alignments in some cases."},
    {ARG_MAX_FAILURES,  0, "f", "max-failures",  Arg::Required,      "-f/--max-failures \t(int, 10).\vSpecifies the maximum number of successive candidate regions that are investigated without success before ALFALFA stops extending the candidate regions of a read. Extension can be restarted only if the remaining candidate regions contain unique seeds."},
    {ARG_MAX_F_BEST,    0, "",  "reset--failures",Arg::None,         "--reset-failures \t.\v If set, the counter of successive candidate regions that are investigated without success is reset if a feasible alignment is found. By default the counter is only reset if a new best alignment is found."},
    {ARG_MIN_ALIGN,     0, "",  "max-alignments",Arg::Required,      "--max-alignments \t(int, 5000).\vSpecifies the maximum number of alignments calculated per read. This value should be higher than the number of reported alignments (option -a). Decreasing this value can increase performance of the algorithm, at the cost of a lower accuracy and worse mapping quality estimation."},
    {ARG_MIN_COVERAGE,  0, "c", "min-coverage",  Arg::Required,      "-c/--min-coverage \t(float, 0.25).\vSpecifies the minimum percentage of the read length that candidate region containing a single seed need to cover before extension of the candidate region is taken into consideration."},
    {ARG_ALWAYS_UNIQUE, 0, "",  "skip-unique",   Arg::Required,      "--skip-unique \t.\vBy default, ALFALFA extends all candidate regions containing unique seeds. If this flag is set, this criterium is not taken into account when deciding to extend a candidate region."},
    {ARG_LOCAL,         0, "",  "local",         Arg::None,          "--local \t.\vBy default, ALFALFA uses global alignment during the last phase of the mapping process. Global alignment in essence is end-to-end alignment, as it entirely covers the read but only covers the reference genome in part. Local alignment is used during the last phase of the mapping process if the --local option is set, which may result in soft clipping of the read."},
    {ARG_BANDWIDTH,     0, "b", "bandwidth",     Arg::Required,      "-b/--bandwidth \t(int, 100).\vSpecifies the maximum bandwidth that is used by the banded alignment algorithm. The bandwidth used is automatically inferred from the specification of the maximum percentage of errors allowed in accepting alignments (option -e), but is bounded by this parameter."},
    {ARG_MATCH,         0, "M", "match",         Arg::Required,      "-M/--match \t(int, 1).\vSpecifies the positive score assigned to matches in the dynamic programming extension phase."},
    {ARG_MISMATCH,      0, "U", "mismatch",      Arg::Required,      "-U/--mismatch \t(int, -4).\vSpecifies the penalty assigned to mismatches in the dynamic programming extension phase."},
    {ARG_GAP_OPEN,      0, "O", "gap-open",      Arg::Required,      "-O/--gap-open \t(int, -6).\vSpecifies the penalty O for opening a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E. The use of affine gap penalties can be disabled by setting this value to zero."},
    {ARG_GAP_EXTEND,    0, "E", "gap-extend",    Arg::Required,      "-E/--gap-extend \t(int, -1).\vSpecifies the penalty E for extending a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E."},
    {ARG_MIN_INSERT,    0, "I", "min-insert",    Arg::Required,      "-I/--min-insert \t(int, 0).\vSpecifies the minimum insert size."},
    {ARG_MAX_INSERT,    0, "X", "max-insert",    Arg::Required,      "-X/--max-insert \t(int, 1000).\vSpecifies the maximum insert size."},
    {ARG_ORIENTATION,   0, "",  "orientation",   Arg::Required,      "--orientation \t(fr | rf | ff, fr).\vSpecifies the orientation of mates. fr means a forward upstream first mate and reverse complemented downstream second mate or vice versa. rf means a reverse complemented upstream first mate and forward downstream second mate or vice versa. ff means a forward upstream first mate and forward downstream second mate or vice versa. Note that these definitions are literally taken over from Bowtie 2."},
    {ARG_NO_MIXED,      0, "",  "no-mixed",      Arg::None,          "--no-mixed \t.\vDisables searching for unpaired alignments."},
    {ARG_NO_DISCORDANT, 0, "",  "no-discordant", Arg::None,          "--no-discordant \t.\vDisables searching for discordant alignments."},
    {ARG_DOVETAIL,      0, "",  "dovetail",      Arg::None,          "--dovetail \t.\vAllows switching between upstream and downstream mates in the definition of their orientation (option --orientation)."},
    {ARG_NO_CONTAIN,    0, "",  "no-contain",    Arg::None,          "--no-contain \t.\vDisallows concordant mates to be fully contained within each other."},
    {ARG_NO_OVERLAP,    0, "",  "no-overlap",    Arg::None,          "--no-overlap \t.\vDisallows concordant mates to overlap each other."},
    {ARG_PAIRED_MODE,   0, "",  "paired-mode",   Arg::Required,      "--paired-mode \t(1 | 2 | 3 | 4 | 5 | 6, 1).\vSpecifies the algorithm used to align paired-end reads. The possible algorithms are discussed in detail in the methods section. Algorithms 1 and 2 do not use information from candidate regions. Algorithms 3 and 4 prioritize extension of candidate regions over both reads. Algorithms 5 and 6 filter the list of candidate regions using the paired-end restraints. Algorithms with an odd number pair mapped reads afterward alignment. Algorithms with an even number perform dynamic programming across a window defined by the insert size restrictions to search for a bridging alignment reaching the other mate."},
    {ARG_PAIRED_RESCUE, 0, "",  "paired-rescue", Arg::None,          "--paired-rescue \t.\vEnable an automatical rescue procedue if no concordant alignment was found using the current parameter settings."},
    {ARG_VERBOSE,       0, "v", "verbose",       Arg::Required,      "-v/--verbose \t(int, 0).\vTurns on lots of progress reporting about the alignment process. Higher numbers give more verbose output. Information is printed to standard error and is useful for debugging purposes. The default value 0 disables progress reporting. The maximum verbosity level is 7."},
    {ARG_DEBUG,         0, "",  "debug",         Arg::Required,      "--debug \t(file).\vSpecifies file to print specific debug information to that can be further examined. This info is structured to allow usage by other (inhouse) debugging tools and is currently experimental."},
    {ARG_HELP,          0, "h", "help",          Arg::None,          "-h/--help \t.\vPrints to standard error the version number, usage description and an overview of the options that can be used to customize the software package."},
    {0,0,0,0,0,0}
};

//functions
mum_t parseSeedType(std::string type);
orientation_t parseOrientation(std::string type);
void processParameters(int argc, char* argv[], mapOptions_t& opt);

#endif	/* OPTIONS_H */

