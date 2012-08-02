/* 
 * File:   options.h
 * Author: mvyvermn
 *
 * Created on 2 augustus 2012, 16:19
 */

#ifndef OPTIONS_H
#define	OPTIONS_H

#include <getopt.h>
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
    orientation_t orientation;
};

enum mum_t { MUM, MAM, MEM, SMAM };
enum command_t {INDEX, MATCHES, ALN};

struct mapOptions_t{//commentary + sort + constructor
    mapOptions_t(){ initOptions(); }
    void initOptions(){
        K = 1; query_threads = 1;
        _4column = false;
        nucleotidesOnly = false;
        alnOptions.minMemLength = 20;
        memType = SMAM;
        noFW = noRC = false;
        alnOptions.errorPercent = 0.08;
        verbose = false;
        alnOptions.scores.match = 0;//(2,-1,-2,-1);
        alnOptions.scores.mismatch = -2;
        alnOptions.scores.openGap = 0;
        alnOptions.scores.extendGap = -2;
        alnOptions.noClipping = true;
        alnOptions.numThreads = 1;
        alnOptions.alignmentCount = 50;
        alnOptions.maxTrial = 10;
        alnOptions.minCoverage = 25;
        alnOptions.tryHarder = false;
        alnOptions.fixedMinLength = false;
        alnOptions.unique = false;
        pairedOpt.contain = true;
        pairedOpt.discordant = true;
        pairedOpt.dovetail = false;
        pairedOpt.maxInsert = 500;
        pairedOpt.minInsert = 0;
        pairedOpt.mixed = true;
        pairedOpt.orientation = PAIR_FR;
        pairedOpt.overlap = true;
    }
    //performance options
    int K;//sparsity factor
    int query_threads;
    //I/O options
    bool _4column;//for MEM output
    bool verbose;
    //Sequence options
    bool nucleotidesOnly;
    //MAM options
    enum mum_t memType;
    //alignment options
    align_opt alnOptions;
    bool noFW; bool noRC;
    //paired end options
    paired_opt pairedOpt;
};

static const char * short_options = "s:k:l:nq:d:m:u:o:e:C:T:hx:U:1:2:S:I:X:";

//Reads from - to ; trim ; phred quals
enum {
    ARG_NO_OPT,          //not found
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
    ARG_NO_OVERLAP       //--no-overlap
};

static struct option long_options[] = {
    {(char*)"sparsityfactor",   required_argument, 0,            's'},
    {(char*)"threads",          required_argument, 0,            'q'},
    {(char*)"verbose",          no_argument,       0,            ARG_VERBOSE},
    {(char*)"seedminlength",    required_argument, 0,            'l'},
    {(char*)"alignments",       required_argument, 0,            'k'},
    {(char*)"trials",           required_argument, 0,            'T'},
    {(char*)"mincoverage",      required_argument, 0,            'C'},
    {(char*)"errors",           required_argument, 0,            'd'},
    {(char*)"tryharder",        no_argument,       0,            ARG_TRY_HARDER},
    {(char*)"seedthreads",      required_argument, 0,            ARG_SEED_THREADS},
    {(char*)"noFw",             no_argument,       0,            ARG_NOFW},
    {(char*)"noRc",             no_argument,       0,            ARG_NORC},
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
    {(char*)0, 0, 0, 0} // terminator
};

#endif	/* OPTIONS_H */

