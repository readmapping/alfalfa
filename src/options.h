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

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "dp.h"
#include <iostream>
#include "optionparser.h"

using option::Option;
using option::Descriptor;
using option::Parser;
using option::Stats;
using option::ArgStatus;

enum mum_t { MUM, MAM, MEM, SMAM };
static const int MINOPTIONCOUNT = 2;
enum command_t {INDEX, MATCHES, ALN};
enum orientation_t { PAIR_FR, PAIR_RF, PAIR_FF};

//match options
struct align_opt {
    dp_scores scores;
    bool scoresSetByUser;
    bool noClipping;
    double errorPercent;
    int minMemLength;
    bool fixedMinLength;
    int maxSeedCandidates;
    int alignmentCount;
    int maxTrial;
    double minCoverage;
    bool tryHarder;
    bool unique;
    bool noFW; 
    bool noRC;
    int print;
    mum_t memType;
    int sparseMult;
    //next parameters are currently only used for local alignment
    int fixedBandSize;
    int minScoreFixed;
    double minScoreLength;
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
    orientation_t orientation;
};

struct mapOptions_t{//commentary + sort + constructor
    mapOptions_t(){ initOptions(); }
    void initOptions(){
        K = 12; query_threads = 1;
        nucleotidesOnly = false;
        alnOptions.minMemLength = 40;
        alnOptions.fixedMinLength = false;
        alnOptions.memType = SMAM;
        alnOptions.sparseMult = 1;
        alnOptions.noFW = alnOptions.noRC = false;
        alnOptions.errorPercent = 0.08;
        alnOptions.print = 0;
        alnOptions.scores.match = 0;//(2,-1,-2,-1);
        alnOptions.scores.mismatch = -2;
        alnOptions.scores.openGap = 0;
        alnOptions.scores.extendGap = -2;
        alnOptions.noClipping = true;
        alnOptions.alignmentCount = 4;
        alnOptions.maxTrial = 5;
        alnOptions.minCoverage = 0.10;
        alnOptions.tryHarder = true;
        alnOptions.unique = false;
        alnOptions.maxSeedCandidates = 1;
        alnOptions.scoresSetByUser = false;
        alnOptions.fixedBandSize = 33;
        alnOptions.minScoreFixed = 37;
        alnOptions.minScoreLength = 5.5;
        pairedOpt.contain = true;
        pairedOpt.discordant = true;
        pairedOpt.dovetail = false;
        pairedOpt.maxInsert = 1000;//return to default!//
        pairedOpt.minInsert = 0;//return to default!//
        pairedOpt.mixed = true;
        pairedOpt.orientation = PAIR_FR;
        pairedOpt.overlap = true;
        pairedOpt.mode = 1;
        unpairedQ = pair1 = pair2 = ref_fasta = outputName = "";
        command = ALN;
        saveIndex = false;
        hasChild = true;
        hasSuflink = false;
        hasKmer = true;
        index_prefix = indexLocation = "";
    }
    void printOptions(){
        cerr << "ALFALFA was executed with following options: " << endl;
        cerr << "reference sequence\t\t" << ref_fasta << endl;
        if(!outputName.empty())
            cerr << "output file\t\t" << outputName << endl;
        if(!indexLocation.empty())
            cerr << "index location\t\t\t" << indexLocation << endl;
        if(!index_prefix.empty())
            cerr << "index prefix\t\t\t" << index_prefix << endl;
        cerr << "save index \t\t\t" << (saveIndex ? "yes" : "no") << endl;
        if(!unpairedQ.empty())
            cerr << "single-end reads\t\t\t" << unpairedQ << endl;
        if(!pair1.empty() && !pair2.empty()){
            cerr << "mate 1   \t\t\t" << pair1 << endl;
            cerr << "mate 2   \t\t\t" << pair2 << endl;
        }
        cerr << "sparseness \t\t\t" << K << endl;
        cerr << "threads  \t\t\t" << query_threads << endl;
        cerr << "edit distance   \t\t\t" << alnOptions.errorPercent << endl;
        cerr << "min seed length\t\t\t" << alnOptions.minMemLength << endl;
        cerr << "auto seed length\t\t\t" << (alnOptions.fixedMinLength ? "no" : "yes") << endl;
        cerr << "max seed candidates\t\t\t" << alnOptions.maxSeedCandidates << endl;
        cerr << "type of seeds\t\t\t" << alnOptions.memType << endl;
        cerr << "save seed finding\t\t\t" << (alnOptions.tryHarder ? "yes" : "no") << endl;
        cerr << "#alignments\t\t\t" << alnOptions.alignmentCount << endl;
        cerr << "#try and error\t\t\t" << alnOptions.maxTrial << endl;
        cerr << "min query coverage\t\t" << alnOptions.minCoverage << endl;
        cerr << "use local alignment?\t\t\t" << (alnOptions.noClipping ? "no" : "yes") << endl;
        if(!alnOptions.noClipping){
            cerr << "use band size\t\t\t" << alnOptions.fixedBandSize << endl;
            cerr << "use minimum score\t\t\t max(" << alnOptions.minScoreFixed << ",log(readlength)*" << alnOptions.minScoreLength << ")" << endl;
        }
        cerr << "verbosity level\t\t\t" << alnOptions.print << endl;
        cerr << "calculate forward\t\t\t" << (alnOptions.noFW ? "no" : "yes") << endl;
        cerr << "calculate reverse\t\t\t" << (alnOptions.noRC ? "no" : "yes") << endl;
        cerr << "match score\t\t\t" << alnOptions.scores.match << endl;
        cerr << "mismatch score\t\t\t" << alnOptions.scores.mismatch << endl;
        cerr << "gap open score\t\t\t" << alnOptions.scores.openGap << endl;
        cerr << "gap extend score\t\t\t" << alnOptions.scores.extendGap << endl;
        if(!pair1.empty() && !pair2.empty()){
            cerr << "paired options: " << endl;
            cerr << "paired mode\t\t\t" << pairedOpt.mode << endl;
            cerr << "type pairing\t\t\t" << pairedOpt.orientation << endl;
            cerr << "min insert size\t\t\t" << pairedOpt.minInsert << endl;
            cerr << "max insert size\t\t\t" << pairedOpt.maxInsert << endl;
            cerr << "allow overlap\t\t\t" << (pairedOpt.overlap ? "yes" : "no") << endl;
            cerr << "allow contained\t\t\t" << (pairedOpt.contain ? "yes" : "no") << endl;
            cerr << "allow discordant\t\t\t" << (pairedOpt.discordant ? "yes" : "no") << endl;
            cerr << "allow dovetail\t\t\t" << (pairedOpt.dovetail ? "yes" : "no") << endl;
            cerr << "allow mixed\t\t\t" << (pairedOpt.mixed ? "yes" : "no") << endl;
        }
    }
    //command option
    command_t command;
    //performance options
    int K;//sparsity factor
    int query_threads;
    //I/O options
    string unpairedQ;
    string pair1;
    string pair2;
    string ref_fasta;
    string index_prefix;
    bool saveIndex;
    bool hasChild;
    bool hasSuflink;
    bool hasKmer;
    string indexLocation;
    string outputName;
    //Sequence options
    bool nucleotidesOnly;
    //alignment options
    align_opt alnOptions;    
    //paired end options
    paired_opt pairedOpt;
};

struct Arg: public option::Arg
{
  static ArgStatus Required(const Option& option, bool)
  {
    return option.arg == 0 ? option::ARG_ILLEGAL : option::ARG_OK;
  }
  static ArgStatus Empty(const Option& option, bool)
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
    ARG_MAX_FAILURES,
    ARG_MIN_COVERAGE,
    ARG_LOCAL,
    ARG_BANDWIDTH,
    ARG_ALPHA,
    ARG_BETA,
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
    ARG_VERBOSE,
    ARG_HELP
};

const option::Descriptor alignUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          0},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nUsage: alfalfa align [option...]\nThe align command is used for mapping and aligning a read set onto a reference genome. As this process can be customized through a long list of options, we have grouped them into several categories."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nI/O options"},
    {ARG_REFERENCE,     0, "r", "reference",     Arg::Required,      "-r/--reference \t(file).\vSpecifies the location of a file that contains the reference genome in multi-fasta format."},
    {ARG_INDEX,         0, "i", "index",         Arg::Required,      "-i/--index \t(string).\vSpecifies the prefix used to name all generated index files. If this option is not set explicitly, an index will be computed from the reference genome according to the settings of the options that also apply to the index command."},
    {ARG_SAVE,          0, "",  "save",          Arg::None,          "--save \t.\vSpecifies that if an index is constructed by the align command itself, it will be stored to disk. This option is ignored if the index is loaded from disk (option -i)."},
    {ARG_SINGLE,        0, "0", "single",        Arg::Required,      "-0/--single \t(file).\vSpecifies the location of a file that contains single-end reads. Both fasta and fastQ formats are accepted. If both single-end and paired-end reads are specified, single-end reads are processed first."},
    {ARG_MATES1,        0, "1", "mates1",        Arg::Required,      "-1/--mates1 \t(file).\vSpecifies the location of a file that contains the first mates of paired-end reads. Both fasta and fastQ formats are accepted."},
    {ARG_MATES2,        0, "2", "mates2",        Arg::Required,      "-2/--mates2 \t(file).\vSpecifies the location of a file that contains the second mates of paired-end reads. Both fasta and fastQ formats are accepted."},
    {ARG_OUTPUT,        0, "o", "output",        Arg::Required,      "-o/--output \t(file, filename passed to the -r option with additional .sam extension).\vSpecifies the location of the generated SAM output file containing the results of read mapping and alignment."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nAlignment options"},
    {ARG_ALIGNMENTS,    0, "a", "alignments",    Arg::Required,      "-a/--alignments \t(int, 4).\vSpecifies the maximum number of alignments reported per read."},
    {ARG_NO_FORWARD,    0, "",  "no-forward",    Arg::None,          "--no-forward \t.\vDo not compute alignments on the forward strand."},
    {ARG_NO_REVERSE,    0, "",  "no-reverse",    Arg::None,          "--no-reverse \t.\vDo not compute alignments on the reverse complement strand."},
    {ARG_EDIT_DISTANCE, 0, "e", "edit-distance", Arg::Required,      "-e/--edit-distance \t(float, 0.08).\vSpecifies the maximum percentage of errors allowed in accepting alignments. For global alignment, this value is taken as a fixed boundary to decide upon accepting alignments. For local alignment, a scoring function is used instead."},
    {ARG_THREADS,       0, "t", "threads",       Arg::Required,      "-t/--threads \t(int, 1).\vNumber of threads used during read mapping. Using more than one thread results in reporting read alignments in a different order compared to the order in which they are read from the input file(s)."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nSeed options"},
    {ARG_SEED,          0, "",  "seed",          Arg::Required,      "--seed \t(MEM or SMEM, SMEM).\vSpecifies the type of seeds used for read mapping. Possible values are MEM for maximal exact matches and SMEM for super-maximal exact matches. The use of SMEMs generally boosts performance without having a negative impact on accuracy compared to the use of MEMs. On the other hand, there are usually many more MEMs than SMEMs, in general resulting in a higher number of candidate genomic regions. The latter might be useful if reporting more candidate mapping locations is preferred."},
    {ARG_MIN_LENGTH,    0, "l", "min-length",    Arg::Required,      "-l/--min-length \t(int, auto).\vSpecifies the minimum seed length. This value must be greater than the sparseness value used to build the index (option -s). By default, the value of this option is computed automatically using the following procedure. A value of 40 is used for reads shorter than 1kbp. The value is incremented by 20 for every 500bp above 1kbp, with the total increment being divided by the maximum percentage of errors allowed in accepting alignments (option -e)."},
    {ARG_MAX_SEEDS,     0, "m", "max-seeds",     Arg::Required,      "-m/--max-seeds \t(int, 1).\vSpecifies the maximum number of seeds that will be selected per starting position in the read sequence. The value passed to this option is multiplied by the automatically computed skip factor that determines sparse matching of sampled suffixes from the read sequence. As a result, the actual number of seeds per starting position in the read might still vary. Higher values of this option result in higher numbers of seeds, increasing in turn the number of candidate genomic regions."},
    {ARG_NO_RESCUE,     0, "",  "no-rescue",     Arg::None,          "--no-rescue \t.\vIf ALFALFA finds no seeds having the minimum seed length specified (option -l), it attempts to find gradually shorter seeds. This behaviour is disabled when specifying the --no-rescue option, which is handy in rare cases where the rescue procedure results in a significant performance drop."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nExtend options"},
    {ARG_MAX_FAILURES,  0, "f", "max-failures",  Arg::Required,      "-f/--max-failures \t(int, 5).\vSpecifies the maximum number of successive candidate regions that are investigated without success before ALFALFA stops extending the seeds of a read."},
    {ARG_MIN_COVERAGE,  0, "c", "min-coverage",  Arg::Required,      "-c/--min-coverage \t(float, 0.10).\vSpecifies the minimum percentage of the read length that seeds within a candidate region need to cover before extension of the candidate region is taken into consideration."},
    {ARG_LOCAL,         0, "",  "local",         Arg::None,          "--local \t.\vBy default, ALFALFA uses global alignment during the last phase of the mapping process. Global alignment in essence is end-to-end alignment, as it entirely covers the read but only covers the reference genome in part. Local alignment is used during the last phase of the mapping process if the --local option is set, which may result in soft clipping of the read. This option determines both the type of dynamic programming and the default scoring function that will be used."},
    {ARG_BANDWIDTH,     0, "b", "bandwidth",     Arg::Required,      "-b/--bandwidth \t(int, 33).\vSpecifies the fixed bandwidth that is used by the banded local alignment algorithm. This option is ignored if the --local option is not set, as the bandwidth used by the global alignment algorithm is automatically inferred from the specification of the maximum percentage of errors allowed in accepting alignments (option -e)."},
    {ARG_ALPHA,         0, "A", "alpha",         Arg::Required,      "-A/--alpha \t(int, 37).\vSpecifies the alpha constant used in calculating the threshold of the scoring function that decides upon accepting alignments when local alignment is used (option --local). For read length L and constants alpha and beta, the minimum score is computed as s*max(alpha, beta*log(L)). The value s is the score given to a match in the alignment."},
    {ARG_BETA,          0, "B", "beta",          Arg::Required,      "-B/--beta \t(float, 5.5).\vSpecifies the beta constant used in calculating the threshold of the scoring function that decides upon accepting alignments when local alignment is used(option --local). For read length L and constants alpha and beta, the minimum score is computed as s*max(alpha, beta*log(L)). The value s is the score given to a match in the alignment."},
    {ARG_MATCH,         0, "M", "match",         Arg::Required,      "-M/--match \t(int, 0 for global alignment and 1 for local alignment).\vSpecifies the positive score assigned to matches in the dynamic programming extension phase."},
    {ARG_MISMATCH,      0, "U", "mismatch",      Arg::Required,      "-U/--mismatch \t(int, -2 for global alignment and -3 for local alignment).\vSpecifies the penalty assigned to mismatches in the dynamic programming extension phase."},
    {ARG_GAP_OPEN,      0, "O", "gap-open",      Arg::Required,      "-O/--gap-open \t(int, 0 for global alignment and -5 for local alignment).\vSpecifies the penalty O for opening a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E. The use of affine gap penalties can be disabled by setting this value to zero."},
    {ARG_GAP_EXTEND,    0, "E", "gap-extend",    Arg::Required,      "-E/--gap-extend \t(int, -2 for global alignment and -2 for local alignment).\vSpecifies the penalty E for extending a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E."},
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          "\nPaired-end mapping options"},
    {ARG_MIN_INSERT,    0, "I", "min-insert",    Arg::Required,      "-I/--min-insert \t(int, 0).\vSpecifies the minimum insert size."},
    {ARG_MAX_INSERT,    0, "X", "max-insert",    Arg::Required,      "-X/--max-insert \t(int, 1000).\vSpecifies the maximum insert size."},
    {ARG_ORIENTATION,   0, "",  "orientation",   Arg::Required,      "--orientation \t(fr | rf | ff, fr).\vSpecifies the orientation of mates. fr means a forward upstream first mate and reverse complemented downstream second mate or vice versa. rf means a reverse complemented upstream first mate and forward downstream second mate or vice versa. ff means a forward upstream first mate and forward downstream second mate or vice versa. Note that these definitions are literally taken over from Bowtie 2."},
    {ARG_NO_MIXED,      0, "",  "no-mixed",      Arg::None,          "--no-mixed \t.\vDisables searching for unpaired alignments."},
    {ARG_NO_DISCORDANT, 0, "",  "no-discordant", Arg::None,          "--no-discordant \t.\vDisables searching for discordant alignments."},
    {ARG_DOVETAIL,      0, "",  "dovetail",      Arg::None,          "--dovetail \t.\vAllows switching between upstream and downstream mates in the definition of their orientation (option --orientation)."},
    {ARG_NO_CONTAIN,    0, "",  "no-contain",    Arg::None,          "--no-contain \t.\vDisallows concordant mates to be fully contained within each other."},
    {ARG_NO_OVERLAP,    0, "",  "no-overlap",    Arg::None,          "--no-overlap \t.\vDisallows concordant mates to overlap each other."},
    {ARG_PAIRED_MODE,   0, "",  "paired-mode",   Arg::Required,      "--paired-mode \t(1 | 2 | 3 | 4, 1).\vSpecifies the algorithm used to align paired-end reads. The four possible algorithms are discussed in detail in the Online Methods. Algorithm 1 maps both mates independently and pairs the mapped reads afterwards. For every alignment of one mate, algorithm 2 performs full dynamic programming across a window defined by the insert size restrictions (options --min-insert and --max-insert) to search for a bridging alignment reaching the other mate. Algorithm 3 independently searches candidate regions for both mates, pairs them and only performs the extension phase for paired candidate regions. The hybrid algorithm 4 proceeds as algorithm 2 but uses the normal extension phase for all candidate regions across the window defined by the insert size restrictions (options --min-insert and --max-insert), instead of full dynamic programming."},
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
    {ARG_VERBOSE,       0, "v", "verbose",       Arg::Required,      "-v/--verbose \t(int, 0).\vTurns on lots of progress reporting about the alignment process. Higher numbers give more verbose output. Information is printed to standard error and is useful for debugging purposes. The default value 0 disables progress reporting. The maximum verbosity level is 2."},
    {ARG_HELP,          0, "h", "help",          Arg::None,          "-h/--help \t.\vPrints to standard error the version number, usage description and an overview of the options that can be used to customize the software package."},
    {0,0,0,0,0,0}
};

const option::Descriptor fullUsage[] = 
{
    {ARG_UNKNOWN,       0, "",  "",              Arg::None,          0},
    {ARG_REFERENCE,     0, "r", "reference",     Arg::Required,      "-r/--reference \t(file).\vSpecifies the location of a file that contains the reference genome in multi-fasta format."},
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
    {ARG_ALIGNMENTS,    0, "a", "alignments",    Arg::Required,      "-a/--alignments \t(int, 4).\vSpecifies the maximum number of alignments reported per read."},
    {ARG_NO_FORWARD,    0, "",  "no-forward",    Arg::None,          "--no-forward \t.\vDo not compute alignments on the forward strand."},
    {ARG_NO_REVERSE,    0, "",  "no-reverse",    Arg::None,          "--no-reverse \t.\vDo not compute alignments on the reverse complement strand."},
    {ARG_EDIT_DISTANCE, 0, "e", "edit-distance", Arg::Required,      "-e/--edit-distance \t(float, 0.08).\vSpecifies the maximum percentage of errors allowed in accepting alignments. For global alignment, this value is taken as a fixed boundary to decide upon accepting alignments. For local alignment, a scoring function is used instead."},
    {ARG_THREADS,       0, "t", "threads",       Arg::Required,      "-t/--threads \t(int, 1).\vNumber of threads used during read mapping. Using more than one thread results in reporting read alignments in a different order compared to the order in which they are read from the input file(s)."},
    {ARG_SEED,          0, "",  "seed",          Arg::Required,      "--seed \t(MEM or SMEM, SMEM).\vSpecifies the type of seeds used for read mapping. Possible values are MEM for maximal exact matches and SMEM for super-maximal exact matches. The use of SMEMs generally boosts performance without having a negative impact on accuracy compared to the use of MEMs. On the other hand, there are usually many more MEMs than SMEMs, in general resulting in a higher number of candidate genomic regions. The latter might be useful if reporting more candidate mapping locations is preferred."},
    {ARG_MIN_LENGTH,    0, "l", "min-length",    Arg::Required,      "-l/--min-length \t(int, auto).\vSpecifies the minimum seed length. This value must be greater than the sparseness value used to build the index (option -s). By default, the value of this option is computed automatically using the following procedure. A value of 40 is used for reads shorter than 1kbp. The value is incremented by 20 for every 500bp above 1kbp, with the total increment being divided by the maximum percentage of errors allowed in accepting alignments (option -e)."},
    {ARG_MAX_SEEDS,     0, "m", "max-seeds",     Arg::Required,      "-m/--max-seeds \t(int, 1).\vSpecifies the maximum number of seeds that will be selected per starting position in the read sequence. The value passed to this option is multiplied by the automatically computed skip factor that determines sparse matching of sampled suffixes from the read sequence. As a result, the actual number of seeds per starting position in the read might still vary. Higher values of this option result in higher numbers of seeds, increasing in turn the number of candidate genomic regions."},
    {ARG_NO_RESCUE,     0, "",  "no-rescue",     Arg::None,          "--no-rescue \t.\vIf ALFALFA finds no seeds having the minimum seed length specified (option -l), it attempts to find gradually shorter seeds. This behaviour is disabled when specifying the --no-rescue option, which is handy in rare cases where the rescue procedure results in a significant performance drop."},
    {ARG_MAX_FAILURES,  0, "f", "max-failures",  Arg::Required,      "-f/--max-failures \t(int, 5).\vSpecifies the maximum number of successive candidate regions that are investigated without success before ALFALFA stops extending the seeds of a read."},
    {ARG_MIN_COVERAGE,  0, "c", "min-coverage",  Arg::Required,      "-c/--min-coverage \t(float, 0.10).\vSpecifies the minimum percentage of the read length that seeds within a candidate region need to cover before extension of the candidate region is taken into consideration."},
    {ARG_LOCAL,         0, "",  "local",         Arg::None,          "--local \t.\vBy default, ALFALFA uses global alignment during the last phase of the mapping process. Global alignment in essence is end-to-end alignment, as it entirely covers the read but only covers the reference genome in part. Local alignment is used during the last phase of the mapping process if the --local option is set, which may result in soft clipping of the read. This option determines both the type of dynamic programming and the default scoring function that will be used."},
    {ARG_BANDWIDTH,     0, "b", "bandwidth",     Arg::Required,      "-b/--bandwidth \t(int, 33).\vSpecifies the fixed bandwidth that is used by the banded local alignment algorithm. This option is ignored if the --local option is not set, as the bandwidth used by the global alignment algorithm is automatically inferred from the specification of the maximum percentage of errors allowed in accepting alignments (option -e)."},
    {ARG_ALPHA,         0, "A", "alpha",         Arg::Required,      "-A/--alpha \t(int, 37).\vSpecifies the alpha constant used in calculating the threshold of the scoring function that decides upon accepting alignments when local alignment is used (option --local). For read length L and constants alpha and beta, the minimum score is computed as s*max(alpha, beta*log(L)). The value s is the score given to a match in the alignment."},
    {ARG_BETA,          0, "B", "beta",          Arg::Required,      "-B/--beta \t(float, 5.5).\vSpecifies the beta constant used in calculating the threshold of the scoring function that decides upon accepting alignments when local alignment is used(option --local). For read length L and constants alpha and beta, the minimum score is computed as s*max(alpha, beta*log(L)). The value s is the score given to a match in the alignment."},
    {ARG_MATCH,         0, "M", "match",         Arg::Required,      "-M/--match \t(int, 0 for global alignment and 1 for local alignment).\vSpecifies the positive score assigned to matches in the dynamic programming extension phase."},
    {ARG_MISMATCH,      0, "U", "mismatch",      Arg::Required,      "-U/--mismatch \t(int, -2 for global alignment and -3 for local alignment).\vSpecifies the penalty assigned to mismatches in the dynamic programming extension phase."},
    {ARG_GAP_OPEN,      0, "O", "gap-open",      Arg::Required,      "-O/--gap-open \t(int, 0 for global alignment and -5 for local alignment).\vSpecifies the penalty O for opening a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E. The use of affine gap penalties can be disabled by setting this value to zero."},
    {ARG_GAP_EXTEND,    0, "E", "gap-extend",    Arg::Required,      "-E/--gap-extend \t(int, -2 for global alignment and -2 for local alignment).\vSpecifies the penalty E for extending a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E."},
    {ARG_MIN_INSERT,    0, "I", "min-insert",    Arg::Required,      "-I/--min-insert \t(int, 0).\vSpecifies the minimum insert size."},
    {ARG_MAX_INSERT,    0, "X", "max-insert",    Arg::Required,      "-X/--max-insert \t(int, 1000).\vSpecifies the maximum insert size."},
    {ARG_ORIENTATION,   0, "",  "orientation",   Arg::Required,      "--orientation \t(fr | rf | ff, fr).\vSpecifies the orientation of mates. fr means a forward upstream first mate and reverse complemented downstream second mate or vice versa. rf means a reverse complemented upstream first mate and forward downstream second mate or vice versa. ff means a forward upstream first mate and forward downstream second mate or vice versa. Note that these definitions are literally taken over from Bowtie 2."},
    {ARG_NO_MIXED,      0, "",  "no-mixed",      Arg::None,          "--no-mixed \t.\vDisables searching for unpaired alignments."},
    {ARG_NO_DISCORDANT, 0, "",  "no-discordant", Arg::None,          "--no-discordant \t.\vDisables searching for discordant alignments."},
    {ARG_DOVETAIL,      0, "",  "dovetail",      Arg::None,          "--dovetail \t.\vAllows switching between upstream and downstream mates in the definition of their orientation (option --orientation)."},
    {ARG_NO_CONTAIN,    0, "",  "no-contain",    Arg::None,          "--no-contain \t.\vDisallows concordant mates to be fully contained within each other."},
    {ARG_NO_OVERLAP,    0, "",  "no-overlap",    Arg::None,          "--no-overlap \t.\vDisallows concordant mates to overlap each other."},
    {ARG_PAIRED_MODE,   0, "",  "paired-mode",   Arg::Required,      "--paired-mode \t(1 | 2 | 3 | 4, 1).\vSpecifies the algorithm used to align paired-end reads. The four possible algorithms are discussed in detail in the Online Methods. Algorithm 1 maps both mates independently and pairs the mapped reads afterwards. For every alignment of one mate, algorithm 2 performs full dynamic programming across a window defined by the insert size restrictions (options --min-insert and --max-insert) to search for a bridging alignment reaching the other mate. Algorithm 3 independently searches candidate regions for both mates, pairs them and only performs the extension phase for paired candidate regions. The hybrid algorithm 4 proceeds as algorithm 2 but uses the normal extension phase for all candidate regions across the window defined by the insert size restrictions (options --min-insert and --max-insert), instead of full dynamic programming."},
    {ARG_VERBOSE,       0, "v", "verbose",       Arg::Required,      "-v/--verbose \t(int, 0).\vTurns on lots of progress reporting about the alignment process. Higher numbers give more verbose output. Information is printed to standard error and is useful for debugging purposes. The default value 0 disables progress reporting. The maximum verbosity level is 2."},
    {ARG_HELP,          0, "h", "help",          Arg::None,          "-h/--help \t.\vPrints to standard error the version number, usage description and an overview of the options that can be used to customize the software package."},
    {0,0,0,0,0,0}
};

static mum_t parseSeedType(string type){
    if(type.compare("MEM")==0)
        return MEM;
    else if(type.compare("MAM")==0)
        return MAM;
    else if(type.compare("MUM")==0)
        return MUM;
    else if(type.compare("SMEM")==0)
        return SMAM;
    else{
        cerr << "WARNING: unrecognized seed type. SMEMs are used instead." << endl;
        return SMAM;
    }
}

static orientation_t parseOrientation(string type){
    if(type.compare("fr")==0)
        return PAIR_FR;
    else if(type.compare("rf")==0)
        return PAIR_RF;
    else if(type.compare("ff")==0)
        return PAIR_FF;
    else{
        cerr << "WARNING: unrecognized orientation. fr is used instead." << endl;
        return PAIR_FR;
    }
}


inline void processParameters(int argc, char* argv[], mapOptions_t& opt){
    //parse command
    if(strcmp(argv[1], "index") == 0){ 
        opt.command = INDEX;
        opt.saveIndex = true;
    }
    else if(strcmp(argv[1], "align") == 0){ 
        opt.command = ALN;
    }
    else{
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[0]);
        option::printUsage(std::cerr, basicUsage);
    }
    cerr << "parsing options: ..." << endl;
    argc-=2; argv++; argv++;
    option::Stats  stats(fullUsage, argc, argv);
    //misc option parsing
    option::Option* options = new option::Option[stats.options_max];
    option::Option* buffer  = new option::Option[stats.buffer_max];
    option::Parser parser(fullUsage, argc, argv, options, buffer);
    if (parser.error()){
        cerr << "error during parsing of command line options" << endl;
        exit(1);
    }
    if ( options[ARG_HELP] ){
        if(opt.command == INDEX){
            option::printUsage(std::cerr, indexUsage);
        }
        else{
            option::printUsage(std::cerr, alignUsage);
        }
        option::printUsage(std::cerr, miscUsage);
        exit(0);
    }
    if( options[ARG_VERBOSE])
        opt.alnOptions.print = atoi(options[ARG_VERBOSE].arg);
    //index parsing options
    for (int i = 0; i < parser.optionsCount(); ++i) {
        Option& option = buffer[i];
        switch(option.index()) {
            case ARG_REFERENCE: opt.ref_fasta = option.arg; break;
            case ARG_SPARSENESS: opt.K = atoi(option.arg); break;
            case ARG_PREFIX: opt.index_prefix = option.arg; break;
            case ARG_NO_CHILD: opt.hasChild = false; break;
            case ARG_SUFLINK: opt.hasSuflink = true; break;
            case ARG_NO_KMER: opt.hasKmer = false; break;
            default: break;
            throw 1;
        }
    }
    //check needed options and correctness
    if (!options[ARG_REFERENCE] ){
        cerr << "ERROR: no reference file given." << endl;
        exit(1);
    }
    if(opt.command == INDEX){
        if(opt.index_prefix.empty()){
            opt.index_prefix = opt.ref_fasta;
        }
        if(opt.hasChild && opt.hasSuflink)
            fprintf(stderr, "WARNING: having both child array and suffix link support is marginally useful and requires more memory! \n");
        if(!opt.hasChild && opt.K>=3)
            fprintf(stderr, "WARNING: no child array means much slower search times! \n");
    }
    //mapping options parsing
    if(opt.command == ALN){
        for (int i = 0; i < parser.optionsCount(); ++i) {
            Option& option = buffer[i];
            switch(option.index()) {
                case ARG_INDEX: opt.indexLocation = option.arg; break;
                case ARG_SAVE: opt.saveIndex = true; break;
                case ARG_SINGLE: opt.unpairedQ = option.arg; break;
                case ARG_MATES1: opt.pair1 = option.arg; break;
                case ARG_MATES2: opt.pair2 = option.arg; break;
                case ARG_OUTPUT: opt.outputName = option.arg; break;
                case ARG_ALIGNMENTS: opt.alnOptions.alignmentCount = atoi(option.arg); break;
                case ARG_NO_FORWARD: opt.alnOptions.noFW = 1; break;
                case ARG_NO_REVERSE: opt.alnOptions.noRC = 1; break;
                case ARG_EDIT_DISTANCE: opt.alnOptions.errorPercent = atof(option.arg); break;
                case ARG_THREADS: opt.query_threads = atoi(option.arg); break;
                case ARG_SEED: opt.alnOptions.memType = parseSeedType(option.arg); break;
                case ARG_MIN_LENGTH: opt.alnOptions.minMemLength = atoi(option.arg); opt.alnOptions.fixedMinLength = true; break;
                case ARG_MAX_SEEDS: opt.alnOptions.maxSeedCandidates = atoi(option.arg); break;
                case ARG_NO_RESCUE: opt.alnOptions.tryHarder = false; break;
                case ARG_MAX_FAILURES: opt.alnOptions.maxTrial = atoi(option.arg); break;
                case ARG_MIN_COVERAGE: opt.alnOptions.minCoverage = atof(option.arg); break;
                case ARG_LOCAL: opt.alnOptions.noClipping = false; break;
                case ARG_BANDWIDTH: opt.alnOptions.fixedBandSize = atoi(option.arg); break;
                case ARG_ALPHA: opt.alnOptions.minScoreFixed = atoi(option.arg); break;
                case ARG_BETA: opt.alnOptions.minScoreLength = atof(option.arg); break;
                case ARG_MATCH: opt.alnOptions.scores.match = atoi(option.arg); opt.alnOptions.scoresSetByUser = true; break;
                case ARG_MISMATCH: opt.alnOptions.scores.mismatch = atoi(option.arg); opt.alnOptions.scoresSetByUser = true; break;
                case ARG_GAP_OPEN: opt.alnOptions.scores.openGap = atoi(option.arg); opt.alnOptions.scoresSetByUser = true; break;
                case ARG_GAP_EXTEND: opt.alnOptions.scores.extendGap = atoi(option.arg); opt.alnOptions.scoresSetByUser = true; break;
                case ARG_MIN_INSERT: opt.pairedOpt.minInsert = atoi(option.arg); break;
                case ARG_MAX_INSERT: opt.pairedOpt.maxInsert = atoi(option.arg); break;
                case ARG_ORIENTATION: opt.pairedOpt.orientation = parseOrientation(option.arg); break;
                case ARG_NO_MIXED: opt.pairedOpt.mixed = false; break;
                case ARG_NO_DISCORDANT: opt.pairedOpt.discordant = false; break;
                case ARG_DOVETAIL: opt.pairedOpt.dovetail = true; break;
                case ARG_NO_CONTAIN: opt.pairedOpt.contain = false; break;
                case ARG_NO_OVERLAP: opt.pairedOpt.overlap = false; break;
                case ARG_PAIRED_MODE: opt.pairedOpt.mode = atoi(option.arg); break;
                default: break;
                throw 1;
            }
        }
        if(opt.indexLocation.empty() && opt.index_prefix.empty()){
            opt.index_prefix = opt.ref_fasta;
        }
        if(opt.hasChild && opt.hasSuflink)
            fprintf(stderr, "WARNING: having both child array and suffix link support is marginally useful and requires more memory! \n");
        if(!opt.hasChild && opt.K>=3)
            fprintf(stderr, "WARNING: no child array means much slower search times! \n");
        if(!opt.saveIndex && opt.indexLocation.empty())
            fprintf(stderr, "WARNING: index will not be saved to disk! \n");
        if(opt.unpairedQ.empty() && (opt.pair1.empty() || opt.pair2.empty())){
            fprintf(stderr, "ERROR: ALFALFA requires query files. Specify query files with -0 or -1,-2 options.\n");
            exit(1);
        }
        if(opt.outputName.empty())
            opt.outputName = opt.ref_fasta.substr().append(".sam");
        if(opt.alnOptions.errorPercent < 0.0 || opt.alnOptions.errorPercent > 1.0){
            fprintf(stderr, "ERROR: uncorrect value for edit distance. Value should lie in the interval [0,1[\n");
            exit(1);
        }
        if(opt.alnOptions.minCoverage < 0.0 || opt.alnOptions.minCoverage > 1.0){
            fprintf(stderr, "ERROR: uncorrect value for minimum coverage. Value should lie in the interval [0,1[\n");
            exit(1);
        }
        if(!opt.alnOptions.noClipping && !opt.alnOptions.scoresSetByUser){
            //local alignment and scores not set: BWA-SW defaults
            opt.alnOptions.scores.match = 1;
            opt.alnOptions.scores.mismatch = -3;
            opt.alnOptions.scores.openGap = -5;
            opt.alnOptions.scores.extendGap = -2;
        }
        if(opt.pairedOpt.minInsert > opt.pairedOpt.maxInsert){
            fprintf(stderr, "ERROR: minimum insert size is larger than maximum insert size\n");
            exit(1);
        }
        if(opt.pairedOpt.mode < 1 || opt.pairedOpt.mode > 4){
            fprintf(stderr, "ERROR: chosen paired-end mode not supported.\n");
            exit(1);
        }
    }
    if(opt.alnOptions.print > 0)
        opt.printOptions();
    delete[] options;
    delete[] buffer;
}

#endif	/* OPTIONS_H */

