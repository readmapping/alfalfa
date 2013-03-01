/* 
 * File:   options.h
 * Author: mvyvermn
 *
 * Created on 2 augustus 2012, 16:19
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
    bool noClipping;
    double errorPercent;
    int minMemLength;
    int maxSeedCandidates;
    int alignmentCount;
    int maxTrial;
    int minCoverage;
    bool tryHarder;
    bool unique;
    bool noFW; 
    bool noRC;
    int print;
    mum_t memType;
    int sparseMult;
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
        K = 1; query_threads = 1;
        nucleotidesOnly = false;
        alnOptions.minMemLength = 20;
        alnOptions.memType = MEM;
        alnOptions.sparseMult = 1;
        alnOptions.noFW = alnOptions.noRC = false;
        alnOptions.errorPercent = 0.08;
        alnOptions.print = 0;
        alnOptions.scores.match = 0;//(2,-1,-2,-1);
        alnOptions.scores.mismatch = -2;
        alnOptions.scores.openGap = 0;
        alnOptions.scores.extendGap = -2;
        alnOptions.noClipping = true;
        alnOptions.alignmentCount = 100;
        alnOptions.maxTrial = 10;
        alnOptions.minCoverage = 25;
        alnOptions.tryHarder = false;
        alnOptions.unique = false;
        alnOptions.maxSeedCandidates = -1;
        pairedOpt.contain = true;
        pairedOpt.discordant = true;
        pairedOpt.dovetail = false;
        pairedOpt.maxInsert = 500;//return to default!//
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
            cerr << "unpaired reads\t\t\t" << unpairedQ << endl;
        if(!pair1.empty() && !pair2.empty()){
            cerr << "mate 1   \t\t\t" << pair1 << endl;
            cerr << "mate 2   \t\t\t" << pair2 << endl;
        }
        cerr << "sparseness \t\t\t" << K << endl;
        cerr << "threads  \t\t\t" << query_threads << endl;
        cerr << "errors   \t\t\t" << alnOptions.errorPercent << endl;
        cerr << "min seed length\t\t\t" << alnOptions.minMemLength << endl;
        cerr << "max seed candidates\t\t\t" << alnOptions.maxSeedCandidates << endl;
        cerr << "type of seeds\t\t\t" << alnOptions.memType << endl;
        cerr << "save seed finding\t\t\t" << (alnOptions.tryHarder ? "yes" : "no") << endl;
        cerr << "use other characters too\t\t" << (nucleotidesOnly ? "no" : "yes") << endl;
        cerr << "#alignments\t\t\t" << alnOptions.alignmentCount << endl;
        cerr << "#try and error\t\t\t" << alnOptions.maxTrial << endl;
        cerr << "min query coverage\t\t" << alnOptions.minCoverage << endl;
        cerr << "use clipping\t\t\t" << (alnOptions.noClipping ? "no" : "yes") << endl;
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
    string indexLocation;
    string outputName;
    //Sequence options
    bool nucleotidesOnly;
    //alignment options
    align_opt alnOptions;    
    //paired end options
    paired_opt pairedOpt;
};

static const char * short_options = "i:s:k:L:np:q:d:m:u:o:e:C:T:hx:U:1:2:S:I:X:M:";

//Reads from - to ; trim ; phred quals
enum {
    ARG_ZERO_OPT = 255,          //not found
    ARG_TRY_HARDER,      //--tryharder
    ARG_NOFW,            //--noFw
    ARG_NORC,            //--noRc
    ARG_CLIP,            //--softclipping
    ARG_VERBOSE,         //--verbose
    ARG_FR,              //--fr
    ARG_RF,              //--rf
    ARG_FF,              //--ff
    ARG_NO_MIXED,        //--no-mixed
    ARG_NO_DISCORDANT,   //--no-discordant
    ARG_DOVETAIL,        //--dovetail
    ARG_NO_CONTAIN,      //--no-contain
    ARG_NO_OVERLAP,      //--no-overlap
    ARG_PAIR_MODE,       //--paired-mode
    ARG_SAVE_INDEX,      //--save
    ARG_SEEDTYPE,         //--memtype
    ARG_CHILD,         //--memtype
    ARG_SUFLINK,         //--memtype
    ARG_VVERBOSE         //--vverbose
};

static struct option long_options[] = {
    {(char*)"sparsityfactor",   required_argument, 0,            's'},
    {(char*)"threads",          required_argument, 0,            'q'},
    {(char*)"verbose",          no_argument,       0,            ARG_VERBOSE},
    {(char*)"vverbose",         no_argument,       0,            ARG_VVERBOSE},
    {(char*)"seedminlength",    required_argument, 0,            'L'},
    {(char*)"alignments",       required_argument, 0,            'k'},
    {(char*)"trials",           required_argument, 0,            'T'},
    {(char*)"mincoverage",      required_argument, 0,            'C'},
    {(char*)"errors",           required_argument, 0,            'd'},
    {(char*)"tryharder",        no_argument,       0,            ARG_TRY_HARDER},
    {(char*)"nofw",             no_argument,       0,            ARG_NOFW},
    {(char*)"norc",             no_argument,       0,            ARG_NORC},
    {(char*)"wildcards",        no_argument,       0,            'N'},
    {(char*)"softclipping",     no_argument,       0,            ARG_CLIP},
    {(char*)"match",            required_argument, 0,            'm'},
    {(char*)"mismatch",         required_argument, 0,            'u'},
    {(char*)"gapopen",          required_argument, 0,            'o'},
    {(char*)"gapextend",        required_argument, 0,            'e'},
    {(char*)"minins",           required_argument, 0,            'I'},
    {(char*)"maxins",           required_argument, 0,            'X'},
    {(char*)"help",             no_argument,       0,            'h'},
    {(char*)"fr",               no_argument,       0,            ARG_FR},
    {(char*)"rf",               no_argument,       0,            ARG_RF},
    {(char*)"ff",               no_argument,       0,            ARG_FF},
    {(char*)"no-mixed",         no_argument,       0,            ARG_NO_MIXED},
    {(char*)"no-discordant",    no_argument,       0,            ARG_NO_DISCORDANT},
    {(char*)"dovetail",         no_argument,       0,            ARG_DOVETAIL},
    {(char*)"no-contain",       no_argument,       0,            ARG_NO_CONTAIN},
    {(char*)"no-overlap",       no_argument,       0,            ARG_NO_OVERLAP},
    {(char*)"paired-mode",      required_argument, 0,            ARG_PAIR_MODE},
    {(char*)"index",            required_argument, 0,            'i'},
    {(char*)"prefix",           required_argument, 0,            'p'},
    {(char*)"seedtype",         required_argument, 0,            ARG_SEEDTYPE},
    {(char*)"save",             required_argument, 0,            ARG_SAVE_INDEX},
    {(char*)"child",            required_argument, 0,            ARG_CHILD},
    {(char*)"suflink",          required_argument, 0,            ARG_SUFLINK},
    {(char*)"seedmaxcand",      required_argument, 0,            'M'},
    {(char*)0, 0, 0, 0} // terminator
};

static void usage(const string prog) {
  cerr << "Usage: " << prog << " COMMAND [options]" << endl;
  cerr << "Command should be one of the following: " << endl;
  cerr << "index                      build the index for given <reference-file> and save to disk" << endl;
  cerr << "                           this command is not necessairy for mapping, as aln can first construct the index" << endl;
  cerr << "aln                        map the reads to the index build for <reference-file>" << endl;
  cerr << "check                      contains several commands for summarizing the accuracy of an output SAM file" << endl;
  cerr << endl;
  cerr << "call " << prog << " COMMAND --help [or -h] for more detailed information" << endl;
  exit(1);
}

static void usageIndex(const string prog) {
  cerr << "Usage: " << prog << " index [options] -x <reference-file>" << endl;
  cerr << "build the index for given <reference-file> and save to disk" << endl;
  cerr << endl;
  cerr << "OPTIONS " << endl;
  cerr << "-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index. "
          << "Note that the value needs to be lower than -L parameter in the ALN command[1]." << endl;
  cerr << "-p/--prefix (string/path)  prefix of the index names [reference-file name]" << endl;
  cerr << "--save (0 or 1)            save index to disk or not [1]" << endl;
  cerr << "--child (0 or 1)           use sparse child array (useful for s>=3 and for MEMs)[1]" << endl;
  cerr << "--suflink (0 or 1)         use suffix links (useful for s<=3 for MAM and SMAM)[1]" << endl;
  exit(1);
}

static void usageAln(const string prog) {
  cerr << "Usage: " << prog << " aln [options]" << endl;
  cerr << "the options are ordered by functionality: " << endl;
  cerr << endl;
  cerr << "I/O OPTIONS " << endl;
  cerr << "-x (string)                reference sequence in multi-fasta" << endl;
  cerr << "-i/--index (string)        prefix or path of the index to load. If not set, index will first be calculated" << endl;
  cerr << "-1                         query file with first mates (fasta or fastq)" << endl;
  cerr << "-2                         query file with second mates (fasta or fastq)" << endl;
  cerr << "-U                         query file with unpaired reads (fasta or fastq)" << endl;
  cerr << "-S                         output file name (will be sam) [referenceName.sam]" << endl;
  cerr << "--save (0 or 1)            if --index is not set, save index to disk [0]" << endl;
  cerr << "-p                         prefix of index that will be saved [reference sequence name]" << endl;
  cerr << endl;
  cerr << "PERFORMANCE OPTIONS " << endl;
  cerr << "-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index if it is not yet constructed [1]." << endl;
  cerr << "-q/--threads (int)    number of threads [1]" << endl;
  cerr << endl;
  cerr << "ALIGNMENT OPTIONS " << endl;
  cerr << "-d/--errors (double)       percentage of errors allowed according to the edit distance [0.08]" << endl;
  cerr << "-L/--seedminlength (int)   minimum length of the seeds used [INT_MAX because seedmaxcand is used]" << endl;
  cerr << "-M/--seedmaxcand (int)     max number of right max matches per query index [max(-k,20)]" << endl;
  cerr << "--seedtype (string)        type of seeds used, choice of MEM,MAM,MUM,SMAM [SMAM]." << endl;
  cerr << "-k/--alignments (int)      expected number of alignments required per strand per read [50]" << endl;
  cerr << "-T/--trials (int)          maximum number of times alignment is attempted before we give up [10]" << endl;
  cerr << "-C/--mincoverage (int)     minimum percent of bases of read the seeds have to cover [25]" << endl;
  cerr << "--tryharder                enable: 'try harder': when no seeds have been found, search using less stringent parameters" << endl;
  cerr << "--nofw                     do not compute forward matches" << endl;
  cerr << "--norc                     do not compute reverse complement matches" << endl;
  cerr << "-n/--wildcards             treat Ns as wildcard characters" << endl;
  cerr << "--softclipping          allow soft clipping at the beginning and end of an alignment" << endl;
  cerr << endl;
  cerr << "DYNAMIC PROGRAMMING OPTIONS " << endl;
  cerr << "-m/--match (int)           match bonus [0]" << endl;
  cerr << "-u/--mismatch (int)        mismatch penalty [-2]" << endl;
  cerr << "-o/--gapopen (int)         gap open penalty (set 0 for non-affine gap-penalties) [0]" << endl;
  cerr << "-e/--gapextend (int)       gap extension penalty [-2]" << endl;
  cerr << endl;
  cerr << "PAIRED END OPTIONS " << endl;
  cerr << "-I/--minins (int)          minimum insert size [0]" << endl;
  cerr << "-X/--maxins (int)          maximum insert sier [500]" << endl;
  cerr << "--fr/--rf/--ff             orientation of the mates: fr means forward upstream mate 1" << endl; 
  cerr << "                           and reverse complement downstream mate 2, or vise versa. rf and ff are similar [fr]" << endl;
  cerr << "--no-mixed                 never search single mate alignments" << endl;
  cerr << "--no-discordant            never discordant alignments" << endl;
  cerr << "--dovetail                 allow reads to dovetail (changing up- and downstream of reads)" << endl;
  cerr << "--no-contain               disallow a mate to be fully contained in the other" << endl;
  cerr << "--no-overlap               disallow a mate to overlap with the other" << endl;
  cerr << "--paired-mode (int)        choose algorithm to calculate paired-end reads" << endl;
  cerr << endl;
  cerr << "MISC OPTIONS " << endl;
  cerr << "--verbose                  enable verbose mode (not by default)" << endl;
  cerr << "--vverbose                 enable very verbose mode (not by default)" << endl;
  cerr << "-h/--help                  print this statement" << endl;
  exit(1);
}

static mum_t parseSeedType(string type){
    if(type.compare("MEM")==0)
        return MEM;
    else if(type.compare("MAM")==0)
        return MAM;
    else if(type.compare("MUM")==0)
        return MUM;
    else if(type.compare("SMAM")==0)
        return SMAM;
    else
        return SMAM;
}

//move process Command to main, 
inline void processParameters(int argc, char* argv[], mapOptions_t& opt, const string program){
    if(argc < MINOPTIONCOUNT)
        usage(program);
    else{
        //parse the command
        if(strcmp(argv[1], "index") == 0){ 
            opt.command = INDEX;
            opt.saveIndex = true;
        }
        else if(strcmp(argv[1], "mem") == 0){
        opt.command = MATCHES;
            printf("This command is not yet supported\n");
            exit(1);
        }
        else if(strcmp(argv[1], "aln") == 0) opt.command = ALN;
        else{
            fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
            usage(program);
        }
        cerr << "COMMAND: " << argv[1] << endl;
        cerr << "parsing options: ..." << endl;
        int option_index = 0;
        int c;
        while ((c = getopt_long(
        argc-1, argv+1,
        short_options, long_options, &option_index)) != -1) {//memType cannot be chosen
            switch (c) {
                case 'x': opt.ref_fasta = optarg; break;
                case 'U': opt.unpairedQ = optarg; break;
                case '1': opt.pair1 = optarg; break;
                case '2': opt.pair2 = optarg; break;
                case 'S': opt.outputName = optarg; break;
                case 's': opt.K = atoi(optarg); break;
                case 'k': opt.alnOptions.alignmentCount = atoi(optarg); break;
                case 'L': opt.alnOptions.minMemLength = atoi(optarg); break;
                case 'M': opt.alnOptions.maxSeedCandidates = atoi(optarg); break;
                case ARG_NOFW: opt.alnOptions.noFW = 1; break;
                case ARG_NORC: opt.alnOptions.noRC = 1; break;
                case 'n': opt.nucleotidesOnly = 1; break;
                case ARG_VERBOSE: opt.alnOptions.print = 1; break;
                case ARG_VVERBOSE: opt.alnOptions.print = 2; break;
                case ARG_CLIP: opt.alnOptions.noClipping = false; break;
                case 'q': opt.query_threads = atoi(optarg); break;
                case 'm': opt.alnOptions.scores.match = atoi(optarg); break;
                case 'u': opt.alnOptions.scores.mismatch = atoi(optarg); break;
                case 'o': opt.alnOptions.scores.openGap = atoi(optarg); break;
                case 'e': opt.alnOptions.scores.extendGap = atoi(optarg); break;
                case 'd': opt.alnOptions.errorPercent = atof(optarg); break;
                case 'T': opt.alnOptions.maxTrial = atoi(optarg); break;
                case 'C': opt.alnOptions.minCoverage = atoi(optarg); break;
                case ARG_TRY_HARDER: opt.alnOptions.tryHarder = true; break;
                case 'h': opt.command == INDEX ? usageIndex(program) : usageAln(program); break;
                case 'I': opt.pairedOpt.minInsert = atoi(optarg); break;
                case 'X': opt.pairedOpt.maxInsert = atoi(optarg); break;
                case ARG_FR: opt.pairedOpt.orientation = PAIR_FR; break;
                case ARG_RF: opt.pairedOpt.orientation = PAIR_RF; break;
                case ARG_FF: opt.pairedOpt.orientation = PAIR_FF; break;
                case ARG_NO_MIXED: opt.pairedOpt.mixed = false; break;
                case ARG_NO_DISCORDANT: opt.pairedOpt.discordant = false; break;
                case ARG_DOVETAIL: opt.pairedOpt.dovetail = true; break;
                case ARG_NO_CONTAIN: opt.pairedOpt.contain = false; break;
                case ARG_NO_OVERLAP: opt.pairedOpt.overlap = false; break;
                case ARG_PAIR_MODE: opt.pairedOpt.mode = atoi(optarg); break;
                case 'i': opt.indexLocation = optarg; break;
                case 'p': opt.index_prefix = optarg; break;
                case ARG_SAVE_INDEX: opt.saveIndex = atoi(optarg); break;
                case ARG_CHILD: opt.hasChild = atoi(optarg); break;
                case ARG_SUFLINK: opt.hasSuflink = atoi(optarg); break;
                case ARG_SEEDTYPE: opt.alnOptions.memType = parseSeedType(optarg); break;
                case -1: /* Done with options. */
                break;
                case 0: if (long_options[option_index].flag != 0) break;
                default: opt.command == INDEX ? usageIndex(program) : usageAln(program);
                throw 1;
            }
        }
        if(opt.alnOptions.maxSeedCandidates == -1){
           opt.alnOptions.maxSeedCandidates = max(opt.alnOptions.alignmentCount,5);
        }
        if(opt.alnOptions.alignmentCount < 1){
            opt.alnOptions.alignmentCount = 1;
            opt.alnOptions.unique = true;
        }
        //add the reference query and output files
        if(opt.ref_fasta.empty()){
            fprintf(stderr, "ALFALFA requires input reference file \n");
            exit(1);
        }
        if(opt.command == INDEX){
            if(opt.saveIndex == false)
                fprintf(stderr, "Warning, index will not be saved to memory! \n");
            if(opt.hasChild && opt.hasSuflink)
                fprintf(stderr, "Warning, having both child array and suffix link support is marginally useful and requires more memory! \n");
            if(!opt.hasChild && opt.K>=3)
                fprintf(stderr, "Warning no child array means much slower search times! \n");
        }
        if(opt.command == ALN && (opt.unpairedQ.empty() && (opt.pair1.empty() || opt.pair2.empty()))){
            fprintf(stderr, "ALFALFA requires query files (1 unpaired and/or 2 mate files \n");
            exit(1);
        }
        if(opt.saveIndex && opt.index_prefix.empty())
                opt.index_prefix = opt.ref_fasta;
        if(opt.K > 1 && (opt.alnOptions.memType == MAM || opt.alnOptions.memType == MUM)){
                fprintf(stderr, "MAM and MUM seeds are not supported with s>1\n");
                exit(1);
        }
        if(opt.alnOptions.print > 0)
            opt.printOptions();
    }
}

#endif	/* OPTIONS_H */

