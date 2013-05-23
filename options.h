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

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "dp.h"

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

static const char * short_options = "r:s:p:i:0:1:2:o:a:e:t:l:m:f:c:b:A:B:M:U:O:E:I:X:v:h";

//Reads from - to ; trim ; phred quals
enum {
    ARG_ZERO_OPT = 255,          //not found
    ARG_TRY_HARDER,      //--no-rescue
    ARG_NOFW,            //--no-forward
    ARG_NORC,            //--no-reverse
    ARG_CLIP,            //--local
    ARG_ORIENTATION,     //--orientation
    ARG_NO_MIXED,        //--no-mixed
    ARG_NO_DISCORDANT,   //--no-discordant
    ARG_DOVETAIL,        //--dovetail
    ARG_NO_CONTAIN,      //--no-contain
    ARG_NO_OVERLAP,      //--no-overlap
    ARG_PAIR_MODE,       //--paired-mode
    ARG_SAVE_INDEX,      //--save
    ARG_SEEDTYPE,        //--seed
    ARG_CHILD,           //--no-child
    ARG_SUFLINK,         //--suflink
    ARG_KMER,            //--no-kmer
};

static struct option long_options[] = {
    {(char*)"reference",        required_argument, 0,            'r'},
    {(char*)"sparseness",       required_argument, 0,            's'},
    {(char*)"prefix",           required_argument, 0,            'p'},
    {(char*)"no-child",         no_argument,       0,            ARG_CHILD},
    {(char*)"suflink",          no_argument,       0,            ARG_SUFLINK},
    {(char*)"no-kmer",          no_argument,       0,            ARG_KMER},
    {(char*)"index",            required_argument, 0,            'i'},
    {(char*)"save",             no_argument,       0,            ARG_SAVE_INDEX},
    {(char*)"single",           required_argument, 0,            '0'},
    {(char*)"mates1",           required_argument, 0,            '1'},
    {(char*)"mates2",           required_argument, 0,            '2'},
    {(char*)"output",           required_argument, 0,            'o'},
    {(char*)"alignments",       required_argument, 0,            'a'},
    {(char*)"no-forward",       no_argument,       0,            ARG_NOFW},
    {(char*)"no-reverse",       no_argument,       0,            ARG_NORC},
    {(char*)"edit-distance",    required_argument, 0,            'e'},
    {(char*)"threads",          required_argument, 0,            't'},
    {(char*)"seed",             required_argument, 0,            ARG_SEEDTYPE},
    {(char*)"min-length",       required_argument, 0,            'l'},
    {(char*)"max-seeds",        required_argument, 0,            'm'},
    {(char*)"no-rescue",        no_argument,       0,            ARG_TRY_HARDER},
    {(char*)"max-failures",     required_argument, 0,            'f'},
    {(char*)"min-coverage",     required_argument, 0,            'c'},
    {(char*)"local",            no_argument,       0,            ARG_CLIP},
    {(char*)"bandwidth",        required_argument, 0,            'b'},
    {(char*)"alpha",            required_argument, 0,            'A'},
    {(char*)"beta",             required_argument, 0,            'B'},
    {(char*)"match",            required_argument, 0,            'M'},
    {(char*)"mismatch",         required_argument, 0,            'U'},
    {(char*)"gap-open",         required_argument, 0,            'O'},
    {(char*)"gap-extend",       required_argument, 0,            'E'},
    {(char*)"min-insert",       required_argument, 0,            'I'},
    {(char*)"max-insert",       required_argument, 0,            'X'},
    {(char*)"orientation",      required_argument, 0,            ARG_ORIENTATION},
    {(char*)"no-mixed",         no_argument,       0,            ARG_NO_MIXED},
    {(char*)"no-discordant",    no_argument,       0,            ARG_NO_DISCORDANT},
    {(char*)"dovetail",         no_argument,       0,            ARG_DOVETAIL},
    {(char*)"no-contain",       no_argument,       0,            ARG_NO_CONTAIN},
    {(char*)"no-overlap",       no_argument,       0,            ARG_NO_OVERLAP},
    {(char*)"paired-mode",      required_argument, 0,            ARG_PAIR_MODE},
    {(char*)"verbose",          required_argument, 0,            'v'},
    {(char*)"help",             no_argument,       0,            'h'},
    {(char*)0, 0, 0, 0} // terminator
};

static void usage() {
  cerr << "Usage: " << "alfalfa <command> [<subcommand>] [option...]" << endl;
  cerr << "Command should be index, align or evaluate" << endl;
  cerr << "Subcommand is only required for the evaluate command" << endl;
  cerr << endl;
  cerr << "commands:" << endl;
  cerr << "index is used to construct the data structures for indexing a given reference genome" << endl;
  cerr << "align is used for mapping and aligning a read set onto a reference genome" << endl;
  cerr << "evaluate is used for evaluating the accuracy of simulated reads and summarizing statistics from the SAM-formatted alignments reported by a read mapper" << endl;
  cerr << endl;
  cerr << "call " << "alfalfa <command> -h/--help for more detailed information on the specific commands" << endl;
  exit(1);
}

static void usageIndex(const string prog) {
  cerr << "Usage: alfalfa index [option...]" << endl;
  cerr << "index is used to construct the data structures for indexing a given reference genome" << endl;
  cerr << endl;
  cerr << "options " << endl;
  cerr << "-r/--reference (file). Specifies the location of a file that contains the reference genome in multi-fasta format." << endl;
  cerr << "-s/--sparseness (int, 12). Specifies the sparseness of the index structure as a way to control part of the speed-memory trade-off." << endl;
  cerr << "-p/--prefix (string, filename passed to the -r option). Specifies the prefix that will be used to name all generated index files. The same prefix has to be passed to the -i option of the align command to load the index structure when mapping reads." << endl;
  cerr << "--no-child. By default, a sparse child array is constructed and stored in an index file with extension .child. The construction of this sparse child array is skipped when the --no-child option is set. This data structure speeds up seed-finding at the cost of (4/s) bytes per base in the reference genome. As the data structure provides a major speed-up, it is advised to have it constructed." << endl;
  cerr << "--suflink. Suffix link support is disabled by default. Suffix link support is enabled when the --suflink option is set, resulting in an index file with extension .isa to be generated. This data structure speeds up seed-finding at the cost of (4/s) bytes per base. It is only useful when sparseness is less than four and minimum seed length is very low (less than 10), because it conflicts with skipping suffixes in matching the read. In practice, this is rarely the case." << endl;
  cerr << "--no-kmer. By default, a 10-mer lookup table is constructed that contains the suffix array interval positions to depth 10 in the virtual suffix tree. It is stored in an index file with extension .kmer and required only 8MB of memory. The construction of this lookup table is skipped when the --no-kmer option is set. The lookup table stores intervals for sequences of length 10 that only contain {A,C,G,T}. This data structure speeds up seed-finding if the minimum seed length is greater than 10." << endl;
  cerr << "-h/--help. Prints to standard error the version number, usage description and an overview of the options that can be used to customize the software package." << endl;
  exit(1);
}

static void usageAln(const string prog) {
  cerr << "Usage: alfalfa align [option...]" << endl;
  cerr << "index is used for mapping and aligning a read set onto a reference genome" << endl;
  cerr << endl;
  cerr << "options" << endl;
  cerr << "I/O options" << endl;
  cerr << "r/--reference (file). Specifies the location of a file that contains the reference genome in multi-fasta format." << endl;
  cerr << "i/--index (string). Specifies the prefix used to name all generated index files. If this option is not set explicitly, an index will be computed from the reference genome according to the settings of the options that also apply to the index command." << endl;
  cerr << "-save. Specifies that if an index is constructed by the align command itself, it will be stored to disk. This option is ignored if the index is loaded from disk (option -i)." << endl;
  cerr << "0/--single (file). Specifies the location of a file that contains single-end reads. Both fasta and fastQ formats are accepted. If both single-end and paired-end reads are specified, single-end reads are processed first." << endl;
  cerr << "1/--mates1 (file). Specifies the location of a file that contains the first mates of paired-end reads. Both fasta and fastQ formats are accepted." << endl;
  cerr << "2/--mates2 (file). Specifies the location of a file that contains the second mates of paired-end reads. Both fasta and fastQ formats are accepted." << endl;
  cerr << "o/--output (file, filename passed to the -r option with additional .sam extension). Specifies the location of the generated SAM output file containing the results of read mapping and alignment." << endl;
  cerr << endl;
  cerr << "Alignment options" << endl;
  cerr << "a/--alignments (int, 4). Specifies the maximum number of alignments reported per read." << endl;
  cerr << "-no-forward. Do not compute alignments on the forward strand." << endl;
  cerr << "-no-reverse. Do not compute alignments on the reverse complement strand." << endl;
  cerr << "e/--edit-distance (float, 0.08). Specifies the maximum percentage of errors allowed in accepting alignments. For global alignment, this value is taken as a fixed boundary to decide upon accepting alignments. For local alignment, a scoring function is used instead." << endl;
  cerr << "t/--threads (int, 1). Number of threads used during read mapping. Using more than one thread results in reporting read alignments in a different order compared to the order in which they are read from the input file(s)." << endl;
  cerr << endl;
  cerr << "Seed options" << endl;
  cerr << "--seed (MEM or SMEM, SMEM). Specifies the type of seeds used for read mapping. Possible values are MEM for maximal exact matches and SMEM for super-maximal exact matches. The use of SMEMs generally boosts performance without having a negative impact on accuracy compared to the use of MEMs. On the other hand, there are usually many more MEMs than SMEMs, in general resulting in a higher number of candidate genomic regions. The latter might be useful if reporting more candidate mapping locations is preferred." << endl;
  cerr << "l/--min-length (int, auto). Specifies the minimum seed length. This value must be greater than the sparseness value used to build the index (option -s). By default, the value of this option is computed automatically using the following procedure. A value of 40 is used for reads shorter than 1kbp. The value is incremented by 20 for every 500bp above 1kbp, with the total increment being divided by the maximum percentage of errors allowed in accepting alignments (option -e)." << endl;
  cerr << "m/--max-seeds (int, 1). Specifies the maximum number of seeds that will be selected per starting position in the read sequence. The value passed to this option is multiplied by the automatically computed skip factor that determines sparse matching of sampled suffixes from the read sequence. As a result, the actual number of seeds per starting position in the read might still vary. Higher values of this option result in higher numbers of seeds, increasing in turn the number of candidate genomic regions." << endl;
  cerr << "-no-rescue. If ALFALFA finds no seeds having the minimum seed length specified (option -l), it attempts to find gradually shorter seeds. This behaviour is disabled when specifying the --no-rescue option, which is handy in rare cases where the rescue procedure results in a significant performance drop." << endl;
  cerr << endl;
  cerr << "Extend options" << endl;
  cerr << "f/--max-failures (int, 5). Specifies the maximum number of successive candidate regions that are investigated without success before ALFALFA stops extending the seeds of a read." << endl;
  cerr << "c/--min-coverage (float, 0.10). Specifies the minimum percentage of the read length that seeds within a candidate region need to cover before extension of the candidate region is taken into consideration." << endl;
  cerr << "-local. By default, ALFALFA uses global alignment during the last phase of the mapping process. Global alignment in essence is end-to-end alignment, as it entirely covers the read but only covers the reference genome in part. Local alignment is used during the last phase of the mapping process if the --local option is set, which may result in soft clipping of the read. This option determines both the type of dynamic programming and the default scoring function that will be used." << endl;
  cerr << "b/--bandwidth (int, 33). Specifies the fixed bandwidth that is used by the banded local alignment algorithm. This option is ignored if the --local option is not set, as the bandwidth used by the global alignment algorithm is automatically inferred from the specification of the maximum percentage of errors allowed in accepting alignments (option -e)." << endl;
  cerr << "A/--alpha (int, 37). Specifies the alpha constant used in calculating the threshold of the scoring function that decides upon accepting alignments when local alignment is used (option --local). For read length L and constants alpha and beta, the minimum score is computed as s*max(alpha, beta*log(L)). The value s is the score given to a match in the alignment." << endl;
  cerr << "B/--beta (float, 5.5). Specifies the beta constant used in calculating the threshold of the scoring function that decides upon accepting alignments when local alignment is used(option --local). For read length L and constants alpha and beta, the minimum score is computed as s*max(alpha, beta*log(L)). The value s is the score given to a match in the alignment." << endl;
  cerr << "M/--match (int, 0 for global alignment and 1 for local alignment). Specifies the positive score assigned to matches in the dynamic programming extension phase." << endl;
  cerr << "U/--mismatch (int, -2 for global alignment and -3 for local alignment). Specifies the penalty assigned to mismatches in the dynamic programming extension phase." << endl;
  cerr << "O/--gap-open (int, 0 for global alignment and -5 for local alignment). Specifies the penalty O for opening a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E. The use of affine gap penalties can be disabled by setting this value to zero." << endl;
  cerr << "E/--gap-extend (int, -2 for global alignment and -2 for local alignment). Specifies the penalty E for extending a gap (insertion or deletion) in the dynamic programming extension phase. The total penalty for a gap of length L equals O + L*E." << endl;
  cerr << endl;
  cerr << "Paired-end mapping options" << endl;
  cerr << "I/--min-insert (int, 0). Specifies the minimum insert size." << endl;
  cerr << "X/--max-insert (int, 1000). Specifies the maximum insert size." << endl;
  cerr << "-orientation (fr | rf | ff, fr). Specifies the orientation of mates. fr means a forward upstream first mate and reverse complemented downstream second mate or vice versa. rf means a reverse complemented upstream first mate and forward downstream second mate or vice versa. ff means a forward upstream first mate and forward downstream second mate or vice versa. Note that these definitions are literally taken over from Bowtie 2." << endl;
  cerr << "-no-mixed. Disables searching for unpaired alignments." << endl;
  cerr << "-no-discordant. Disables searching for discordant alignments." << endl;
  cerr << "-dovetail. Allows switching between upstream and downstream mates in the definition of their orientation (option --orientation)." << endl;
  cerr << "-no-contain. Disallows concordant mates to be fully contained within each other." << endl;
  cerr << "-no-overlap. Disallows concordant mates to overlap each other." << endl;
  cerr << "-paired-mode (1 | 2 | 3 | 4, 1). Specifies the algorithm used to align paired-end reads. The four possible algorithms are discussed in detail in the Online Methods. Algorithm 1 maps both mates independently and pairs the mapped reads afterwards. For every alignment of one mate, algorithm 2 performs full dynamic programming across a window defined by the insert size restrictions (options --min-insert and --max-insert) to search for a bridging alignment reaching the other mate. Algorithm 3 independently searches candidate regions for both mates, pairs them and only performs the extension phase for paired candidate regions. The hybrid algorithm 4 proceeds as algorithm 2 but uses the normal extension phase for all candidate regions across the window defined by the insert size restrictions (options --min-insert and --max-insert), instead of full dynamic programming." << endl;
  cerr << endl;
  cerr << "Miscellaneous  options" << endl;
  cerr << "v/--verbose (int, 0). Turns on lots of progress reporting about the alignment process. Higher numbers give more verbose output. Information is printed to standard error and is useful for debugging purposes. The default value 0 disables progress reporting. The maximum verbosity level is 2." << endl;
  cerr << "h/--help. Prints to standard error the version number, usage description and an overview of the options that can be used to customize the software package." << endl;
  exit(1);
}

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

//move process Command to main, 
inline void processParameters(int argc, char* argv[], mapOptions_t& opt, const string program){
    if(argc < MINOPTIONCOUNT)
        usage();
    else{
        //parse the command
        if(strcmp(argv[1], "index") == 0){ 
            opt.command = INDEX;
            opt.saveIndex = true;
        }
        else if(strcmp(argv[1], "align") == 0) opt.command = ALN;
        else{
            fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
            usage();
        }
        cerr << "COMMAND: " << argv[1] << endl;
        cerr << "parsing options: ..." << endl;
        int option_index = 0;
        int c;
        while ((c = getopt_long(
        argc-1, argv+1,
        short_options, long_options, &option_index)) != -1) {//memType cannot be chosen
            switch (c) {
                case 'r': opt.ref_fasta = optarg; break;
                case 's': opt.K = atoi(optarg); break;
                case 'p': opt.index_prefix = optarg; break;
                case ARG_CHILD: opt.hasChild = false; break;
                case ARG_SUFLINK: opt.hasSuflink = true; break;
                case ARG_KMER: opt.hasKmer = false; break;
                case 'i': opt.indexLocation = optarg; break;
                case ARG_SAVE_INDEX: opt.saveIndex = true; break;
                case '0': opt.unpairedQ = optarg; break;
                case '1': opt.pair1 = optarg; break;
                case '2': opt.pair2 = optarg; break;
                case 'o': opt.outputName = optarg; break;
                case 'a': opt.alnOptions.alignmentCount = atoi(optarg); break;
                case ARG_NOFW: opt.alnOptions.noFW = 1; break;
                case ARG_NORC: opt.alnOptions.noRC = 1; break;
                case 'e': opt.alnOptions.errorPercent = atof(optarg); break;
                case 't': opt.query_threads = atoi(optarg); break;
                case ARG_SEEDTYPE: opt.alnOptions.memType = parseSeedType(optarg); break;
                case 'l': opt.alnOptions.minMemLength = atoi(optarg); opt.alnOptions.fixedMinLength = true; break;
                case 'm': opt.alnOptions.maxSeedCandidates = atoi(optarg); break;
                case ARG_TRY_HARDER: opt.alnOptions.tryHarder = false; break;
                case 'f': opt.alnOptions.maxTrial = atoi(optarg); break;
                case 'c': opt.alnOptions.minCoverage = atof(optarg); break;
                case ARG_CLIP: opt.alnOptions.noClipping = false; break;
                case 'b': opt.alnOptions.fixedBandSize = atoi(optarg); break;
                case 'A': opt.alnOptions.minScoreFixed = atoi(optarg); break;
                case 'B': opt.alnOptions.minScoreLength = atof(optarg); break;
                case 'M': opt.alnOptions.scores.match = atoi(optarg); opt.alnOptions.scoresSetByUser = true; break;
                case 'U': opt.alnOptions.scores.mismatch = atoi(optarg); opt.alnOptions.scoresSetByUser = true; break;
                case 'O': opt.alnOptions.scores.openGap = atoi(optarg); opt.alnOptions.scoresSetByUser = true; break;
                case 'E': opt.alnOptions.scores.extendGap = atoi(optarg); opt.alnOptions.scoresSetByUser = true; break;
                case 'I': opt.pairedOpt.minInsert = atoi(optarg); break;
                case 'X': opt.pairedOpt.maxInsert = atoi(optarg); break;
                case ARG_ORIENTATION: opt.pairedOpt.orientation = parseOrientation(optarg); break;
                case ARG_NO_MIXED: opt.pairedOpt.mixed = false; break;
                case ARG_NO_DISCORDANT: opt.pairedOpt.discordant = false; break;
                case ARG_DOVETAIL: opt.pairedOpt.dovetail = true; break;
                case ARG_NO_CONTAIN: opt.pairedOpt.contain = false; break;
                case ARG_NO_OVERLAP: opt.pairedOpt.overlap = false; break;
                case ARG_PAIR_MODE: opt.pairedOpt.mode = atoi(optarg); break;
                case 'v': opt.alnOptions.print = atoi(optarg); break;
                case 'h': opt.command == INDEX ? usageIndex(program) : usageAln(program); break;
                case -1: /* Done with options. */
                break;
                case 0: if (long_options[option_index].flag != 0) break;
                default: opt.command == INDEX ? usageIndex(program) : usageAln(program);
                throw 1;
            }
        }
        //setting defaults and checking necessary parameters. Checking correctness parameters + giving warnings
        if(opt.ref_fasta.empty()){
            cerr << "ERROR: no reference file given. Specify reference file with -r option." << endl;
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
        if(opt.command == ALN){
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
    }
}

#endif	/* OPTIONS_H */

