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
#include "dp.h"

//match options
struct align_opt {
    dp_scores scores;
    bool noClipping;
    double errorPercent;
    int minMemLength;
    int numThreads;
    int alignmentCount;
    int maxTrial;
    int minCoverage;
    bool tryHarder;
    bool fixedMinLength;
    bool unique;
    bool noFW; 
    bool noRC;
};

enum orientation_t { PAIR_FR, PAIR_RF, PAIR_FF};

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

enum mum_t { MUM, MAM, MEM, SMAM };
static const int MINOPTIONCOUNT = 2;
enum command_t {INDEX, MATCHES, ALN};

struct mapOptions_t{//commentary + sort + constructor
    mapOptions_t(){ initOptions(); }
    void initOptions(){
        K = 1; query_threads = 1;
        _4column = false;
        nucleotidesOnly = false;
        alnOptions.minMemLength = 20;
        memType = SMAM;
        alnOptions.noFW = alnOptions.noRC = false;
        alnOptions.errorPercent = 0.08;
        verbose = false;
        alnOptions.scores.match = 0;//(2,-1,-2,-1);
        alnOptions.scores.mismatch = -2;
        alnOptions.scores.openGap = 0;
        alnOptions.scores.extendGap = -2;
        alnOptions.noClipping = true;
        alnOptions.numThreads = 1;
        alnOptions.alignmentCount = 100;
        alnOptions.maxTrial = 10;
        alnOptions.minCoverage = 25;
        alnOptions.tryHarder = false;
        alnOptions.fixedMinLength = false;
        alnOptions.unique = false;
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
    string outputName;
    bool _4column;//for MEM output
    bool verbose;
    //Sequence options
    bool nucleotidesOnly;
    //MAM options
    enum mum_t memType;
    //alignment options
    align_opt alnOptions;    
    //paired end options
    paired_opt pairedOpt;
};

static const char * short_options = "s:k:L:np:d:m:u:o:e:C:T:hx:U:1:2:S:I:X:";

//Reads from - to ; trim ; phred quals
enum {
    ARG_ZERO_OPT = 255,          //not found
    ARG_TRY_HARDER,      //--tryharder
    ARG_SEED_THREADS,    //--seedthreads
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
    ARG_PAIR_MODE        //--paired-mode
};

static struct option long_options[] = {
    {(char*)"sparsityfactor",   required_argument, 0,            's'},
    {(char*)"threads",          required_argument, 0,            'p'},
    {(char*)"verbose",          no_argument,       0,            ARG_VERBOSE},
    {(char*)"seedminlength",    required_argument, 0,            'L'},
    {(char*)"alignments",       required_argument, 0,            'k'},
    {(char*)"trials",           required_argument, 0,            'T'},
    {(char*)"mincoverage",      required_argument, 0,            'C'},
    {(char*)"errors",           required_argument, 0,            'd'},
    {(char*)"tryharder",        no_argument,       0,            ARG_TRY_HARDER},
    {(char*)"seedthreads",      required_argument, 0,            ARG_SEED_THREADS},
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
    {(char*)0, 0, 0, 0} // terminator
};

static void usage(const string prog) {
  cerr << "Usage: " << prog << " COMMAND [options]" << endl;
  cerr << "Command should be one of the following: " << endl;
  cerr << "index                      only build the index for given <reference-file>, used for time calculations" << endl;
  cerr << "                           this command is not necessairy for mapping, as aln first constructs the index too" << endl;
  cerr << "aln                        map the reads to the index build for <reference-file>" << endl;
  cerr << "check                      contains several commands for summarizing the accuracy of an output SAM file" << endl;
  cerr << endl;
  cerr << "call " << prog << " COMMAND --help [or -h] for more detailed information" << endl;
  exit(1);
}

static void usageIndex(const string prog) {
  cerr << "Usage: " << prog << " index -x <reference-file>" << endl;
  cerr << "Only build the index for given <reference-file>, used for time calculations." << endl;
  cerr << "This command is not necessairy for mapping, as aln first constructs the index too" << endl;
  exit(1);
}

static void usageAln(const string prog) {
  cerr << "Usage: " << prog << " aln [options]" << endl;
  cerr << "the options for are ordered by functionality: " << endl;
  cerr << endl;
  cerr << "I/O OPTIONS " << endl;
  cerr << "-x (string)                reference sequence in mult-fasta" << endl;
  cerr << "-1                         query file with first mates (fasta or fastq)" << endl;
  cerr << "-2                         query file with second mates (fasta or fastq)" << endl;
  cerr << "-U                         query file with unpaired reads (fasta or fastq)" << endl;
  cerr << "-S                         output file name (will be sam) [referenceName.sam]" << endl;
  cerr << endl;
  cerr << "PERFORMANCE OPTIONS " << endl;
  cerr << "-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index, value needs to be lower than -L parameter [1]" << endl;
  cerr << "-p/--threads (int)    number of threads [1]" << endl;
  cerr << endl;
  cerr << "ALIGNMENT OPTIONS " << endl;
  cerr << "-d/--errors (double)       percentage of errors allowed according to the edit distance [0.08]" << endl;
  cerr << "-L/--seedminlength (int)   minimum length of the seeds used [depending on errorPercent and read length, [min 20]" << endl;
  cerr << "-k/--alignments (int)      expected number of alignments required per strand per read [50]" << endl;
  cerr << "-T/--trials (int)          maximum number of times alignment is attempted before we give up [10]" << endl;
  cerr << "-C/--mincoverage (int)     minimum percent of bases of read the seeds have to cover [25]" << endl;
  cerr << "--tryharder                enable: 'try harder': when no seeds have been found, search using less stringent parameters" << endl;
  cerr << "--seedthreads (int)        number of threads for calculating the seeds [1]" << endl;
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
  cerr << "--verbose               enable verbose mode (not by default)" << endl;
  cerr << "-h/--help                  print this statement" << endl;
  exit(1);
}

//move process Command to main, 
static void processParameters(int argc, char* argv[], mapOptions_t& opt, const string program){
    if(argc < MINOPTIONCOUNT)
        usage(program);
    else{
        //parse the command
        if(strcmp(argv[1], "index") == 0) opt.command = INDEX;
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
                case 'L': opt.alnOptions.minMemLength = atoi(optarg); opt.alnOptions.fixedMinLength = true; break;
                case ARG_NOFW: opt.alnOptions.noFW = 1; break;
                case ARG_NORC: opt.alnOptions.noRC = 1; break;
                case 'n': opt.nucleotidesOnly = 1; break;
                case ARG_VERBOSE: opt.verbose = 1; break;
                case ARG_CLIP: opt.alnOptions.noClipping = false; break;
                case 't': opt.alnOptions.numThreads = atoi(optarg); break;
                case 'p': opt.query_threads = atoi(optarg); break;
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
                case -1: /* Done with options. */
                break;
                case 0: if (long_options[option_index].flag != 0) break;
                default: opt.command == INDEX ? usageIndex(program) : usageAln(program);
                throw 1;
            }
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
        if(opt.command == ALN && (opt.unpairedQ.empty() && (opt.pair1.empty() || opt.pair2.empty()))){
            fprintf(stderr, "ALFALFA requires query files (1 unpaired and/or 2 mate files \n");
            exit(1);
        }
    }
}

#endif	/* OPTIONS_H */

