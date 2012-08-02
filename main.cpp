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
 * 
 * Parts of the code are based on "A Practical Algorithm for Finding 
 * Maximal Exact Matches in Large Sequence Data Sets Using Sparse Suffix Arrays"
 * By Khan et al. Copyright (c) 2009 
 * You should have received a copy of the copyright notice with this code.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "sparseSA.h"
#include "fasta.h"
#include "utils.h"

#include <getopt.h>
#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <string.h>
#include <bitset>

#include <cctype> // std::tolower(), uppercase/lowercase conversion

#include "dp.h"

// NOTE use of special characters ~, `, and $ !!!!!!!!

using namespace std;

void usage(string prog);

enum mum_t { MUM, MAM, MEM, SMAM };
enum command_t {INDEX, MATCHES, ALN};
//mapper options

int MINOPTIONCOUNT = 2;
static const string PROG = "ALFALFA";
static const string SAM_VERSION = "1.4";
static const string PROG_VERSION = "0.3.5";
static string NAN = "*";

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
//output struct
struct samOutput{//construct
    string header;
    vector<string> refdescr;
    vector<long> startpos;
    void createHeader(int argc, char* argv[], long refLength){
        stringstream ss;
        ss << "@HD\tVN:" << SAM_VERSION << endl;
        for(int i=0; i < refdescr.size()-1; i++){
            long length = startpos[i+1]-startpos[i]-1;
            ss << "@SQ\tSN:" << refdescr[i] << "\tLN:" << length << endl;//startpos contains refSequences with '`' in between.
        }
        long length = refLength-startpos[refdescr.size()-1];
        ss << "@SQ\tSN:" << refdescr[refdescr.size()-1] << "\tLN:" << length << endl;
        ss << "@PG\tID:" << PROG << "\tCL:";
        for(int i=0; i<argc; i++){
            if(i > 0)
                ss << " ";
            ss << argv[i];//add number to the stream
        }
        ss << "\tVN:" << PROG_VERSION << endl;
        header = ss.str();
    }
};

//FIELDS
static sparseSA *sa;
static fastqInputReader *queryReader;
static FILE * outfile;
static pthread_mutex_t writeLock_;

static fastqInputReader *mate1Reader;
static fastqInputReader *mate2Reader;

//I/O NAMES
static string unpairedQ = "";
static string pair1 = "";
static string pair2 = "";
static string ref_fasta = "";
static string outputName = "";

struct query_arg {
    query_arg(): skip0(0), skip(0), opt(0){}
  int skip0;
  int skip;
  mapOptions_t * opt;
  pthread_mutex_t *readLock;
  pthread_mutex_t *writeLock;
};

void *query_thread(void *arg_) {

  query_arg *arg = (query_arg *)arg_;
  bool print = arg->opt->verbose;
  long seq_cnt = 0;
  bool hasRead = true;
  while(hasRead){
      read_t read;
      pthread_mutex_lock(arg->readLock);
      hasRead = queryReader->nextRead(read.qname,read.sequence,read.qual);
      pthread_mutex_unlock(arg->readLock);
      if(hasRead){
          seq_cnt++;
          if(!arg->opt->noFW)
                sa->inexactMatch(read, arg->opt->alnOptions, true, print);
          if(!arg->opt->noRC)
                sa->inexactMatch(read, arg->opt->alnOptions, false, print);
          read.postprocess(arg->opt->alnOptions.scores);
          //From global to local pos and write to stringstream
          stringstream * ss = new stringstream;
          if(read.alignments.empty())
              *ss << read.qname << "\t4\t*\t0\t0\t*\t*\t0\t0\t" 
                      << read.sequence << "\t" << read.qual << endl;
          else{
                long globPos;
                string revCompl = read.sequence;
                //TODO: do this during thread-work, perhaps has worse locality for startpos and refdescr tables
                Utils::reverse_complement(revCompl, false);
                string qualRC = read.qual;
                reverse(qualRC.begin(),qualRC.end());
                for(int k = 0; k < read.alignments.size(); k++){
                        alignment_t & a = read.alignments[k];
                        globPos = a.pos;
                        long descIndex;
                        sa->from_set(globPos, descIndex, a.pos);
                        a.rname = sa->descr[descIndex];
                        *ss << read.qname << "\t" << a.flag.to_ulong() << "\t"
                                << a.rname << "\t" << a.pos << "\t" << 
                                a.mapq << "\t" << a.cigar << "\t" << 
                                a.rnext << "\t" << a.pnext << "\t" << 
                                a.tLength << "\t" << (a.flag.test(4) ? revCompl : read.sequence) <<
                                "\t" << (a.flag.test(4) ? qualRC : read.qual) << "\tAS:i:" << 
                                a.alignmentScore << "\tNM:i:" << a.editDist << "\tX0:Z:" << 
                                a.NMtag << endl;
                }
          }
          pthread_mutex_lock(arg->writeLock);
          fprintf(outfile,"%s",ss->str().c_str());
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences mapped: %ld\n", seq_cnt);
  
  pthread_exit(NULL);
}

void *paired_thread(void *arg_) {

  query_arg *arg = (query_arg *)arg_;
  bool print = arg->opt->verbose;
  long seq_cnt = 0;
  bool hasRead = true;
//  while(hasRead){
//      read_t read;
//      pthread_mutex_lock(arg->readLock);
//      hasRead = queryReader->nextRead(read.qname,read.sequence,read.qual);
//      pthread_mutex_unlock(arg->readLock);
//      if(hasRead){
//          seq_cnt++;
//          if(!arg->opt->noFW)
//                sa->inexactMatch(read, arg->opt->alnOptions, true, print);
//          if(!arg->opt->noRC)
//                sa->inexactMatch(read, arg->opt->alnOptions, false, print);
//          read.postprocess(arg->opt->alnOptions.scores);
//          //From global to local pos and write to stringstream
//          stringstream * ss = new stringstream;
//          if(read.alignments.empty())
//              *ss << read.qname << "\t4\t*\t0\t0\t*\t*\t0\t0\t" 
//                      << read.sequence << "\t" << read.qual << endl;
//          else{
//                long globPos;
//                string revCompl = read.sequence;
//                //TODO: do this during thread-work, perhaps has worse locality for startpos and refdescr tables
//                Utils::reverse_complement(revCompl, false);
//                string qualRC = read.qual;
//                reverse(qualRC.begin(),qualRC.end());
//                for(int k = 0; k < read.alignments.size(); k++){
//                        alignment_t & a = read.alignments[k];
//                        globPos = a.pos;
//                        long descIndex;
//                        sa->from_set(globPos, descIndex, a.pos);
//                        a.rname = sa->descr[descIndex];
//                        *ss << read.qname << "\t" << a.flag.to_ulong() << "\t"
//                                << a.rname << "\t" << a.pos << "\t" << 
//                                a.mapq << "\t" << a.cigar << "\t" << 
//                                a.rnext << "\t" << a.pnext << "\t" << 
//                                a.tLength << "\t" << (a.flag.test(4) ? revCompl : read.sequence) <<
//                                "\t" << (a.flag.test(4) ? qualRC : read.qual) << "\tAS:i:" << 
//                                a.alignmentScore << "\tNM:i:" << a.editDist << "\tX0:Z:" << 
//                                a.NMtag << endl;
//                }
//          }
//          pthread_mutex_lock(arg->writeLock);
//          fprintf(outfile,"%s",ss->str().c_str());
//          pthread_mutex_unlock(arg->writeLock);
//          delete ss;
//      }
//  }
  printf("sequences mapped: %ld\n", seq_cnt);
  
  pthread_exit(NULL);
}

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

void usage(const string prog) {
  cerr << "Usage: " << prog << " COMMAND [options] -x <reference-file> -1 <m1> -2 <m2> -U <unpaired> [-S ouput-file]" << endl;
  cerr << "Command should be one of the following: " << endl;
  cerr << "index                      only build the index for given <reference-file>, used for time calculations" << endl;
  cerr << "aln                        map the reads to the index build for <reference-file>" << endl;
  cerr << "the options are ordered by functionality: " << endl;
  cerr << endl;
  cerr << "I/O OPTIONS " << endl;
  cerr << "-x (string)                reference sequence in mult-fasta" << endl;
  cerr << "-1                         query file with first mates (fasta or fastq)" << endl;
  cerr << "-2                         query file with second mates (fasta or fastq)" << endl;
  cerr << "-U                         query file with unpaired reads (fasta or fastq)" << endl;
  cerr << "-S                         output file name (will be sam)" << endl;
  cerr << endl;
  cerr << "PERFORMANCE OPTIONS " << endl;
  cerr << "-s/--sparsityfactor (int)  the sparsity factor of the sparse suffix array index, values between 1 and 4 are preferred [1]" << endl;
  cerr << "-q/--threads (int)    number of threads [1]" << endl;
  cerr << endl;
  cerr << "ALIGNMENT OPTIONS " << endl;
  cerr << "-d/--errors (double)       percentage of errors allowed according to the edit distance [0.08]" << endl;
  cerr << "-l/--seedminlength (int)   minimum length of the seeds used [depending on errorPercent and read length, min 20]" << endl;
  cerr << "-k/--alignments (int)      expected number of alignments required per strand per read [50]" << endl;
  cerr << "-T/--trials (int)          maximum number of times alignment is attempted before we give up [10]" << endl;
  cerr << "-C/--mincoverage (int)     minimum percent of bases of read the seeds have to cover [25]" << endl;
  cerr << "--tryharder                enable: 'try harder': when no seeds have been found, search using less stringent parameters" << endl;
  cerr << "--seedthreads (int)        number of threads for calculating the seeds [1]" << endl;
  cerr << "--noFw                     do not compute forward matches" << endl;
  cerr << "--noRc                     do not compute reverse complement matches" << endl;
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
  cerr << endl;
  cerr << "MISC OPTIONS " << endl;
  cerr << "--verbose               enable verbose mode (not by default)" << endl;
  cerr << "-h/--help                  print this statement" << endl;
  exit(1);
}

int main(int argc, char* argv[]){
    if(argc < MINOPTIONCOUNT)
        usage(PROG);
    else{
        mapOptions_t opt;
        samOutput output;
        opt.initOptions();
        //parse the command
        enum command_t command;
        if(strcmp(argv[1], "index") == 0) command = INDEX;
        else if(strcmp(argv[1], "mem") == 0){
            command = MATCHES;
             printf("This command is not yet supported\n");
             exit(1);
        }
        else if(strcmp(argv[1], "aln") == 0) command = ALN;
        else{
            fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
            exit(1);
        }
        cerr << "@PG\tID:" << PROG << "\tVN:" << PROG_VERSION << endl;
        cerr << "COMMAND: " << argv[1] << endl;
        cerr << "parsing options: ..." << endl;
        int option_index = 0;
        int c;
        while ((c = getopt_long(
			argc-1, argv+1,
			short_options, long_options, &option_index)) != -1) {//memType cannot be chosen
            switch (c) {
                case 'x': ref_fasta = optarg; break;
                case 'U': unpairedQ = optarg; break;
                case '1': pair1 = optarg; break;
                case '2': pair2 = optarg; break;
                case 'S': outputName = optarg; break;
                case 's': opt.K = atoi(optarg); break;
                case 'k': opt.alnOptions.alignmentCount = atoi(optarg); break;
                case 'l': opt.alnOptions.minMemLength = atoi(optarg); opt.alnOptions.fixedMinLength = true; break;
                case ARG_NOFW: opt.noFW = 1; break;
                case ARG_NORC: opt.noRC = 1; break;
                case 'n': opt.nucleotidesOnly = 1; break;
                case ARG_VERBOSE: opt.verbose = 1; break;
                case ARG_CLIP: opt.alnOptions.noClipping = false; break;
                case 't': opt.alnOptions.numThreads = atoi(optarg); break;
                case 'q': opt.query_threads = atoi(optarg); break;
                case 'm': opt.alnOptions.scores.match = atoi(optarg); break;
                case 'u': opt.alnOptions.scores.mismatch = atoi(optarg); break;
                case 'o': opt.alnOptions.scores.openGap = atoi(optarg); break;
                case 'e': opt.alnOptions.scores.extendGap = atoi(optarg); break;
                case 'd': opt.alnOptions.errorPercent = atof(optarg); break;
                case 'T': opt.alnOptions.maxTrial = atoi(optarg); break;
                case 'C': opt.alnOptions.minCoverage = atoi(optarg); break;
                case ARG_TRY_HARDER: opt.alnOptions.tryHarder = true; break;
                case 'h': usage(PROG); break;
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
                case -1: /* Done with options. */
				break;
                case 0: if (long_options[option_index].flag != 0) break;
                default: usage(PROG);
				throw 1;
            }
        }
        opt.alnOptions.scores.updateScoreMatrixDna();
        initDPMatrix(2048, opt.alnOptions.scores.openGap != 0);
        if(opt.alnOptions.alignmentCount < 1){
            opt.alnOptions.alignmentCount = 1;
            opt.alnOptions.unique = true;
        }
        //add the reference query and output files
        if(ref_fasta.empty()){
            fprintf(stderr, "ALFALFA requires input reference file \n");
            exit(1);
        }
        //query name + read initialization
        if(command == ALN && (unpairedQ.empty() && (pair1.empty() || pair2.empty()))){
            fprintf(stderr, "ALFALFA requires query files (1 unpaired and/or 2 mate files \n");
            exit(1);
        }
        cerr << "parsing options: done" << endl;
        cerr << "loading ref sequences: ..." << endl;
        string ref;
        load_fasta(ref_fasta, ref, output.refdescr, output.startpos);
        cerr << "loading ref sequences: done" << endl;
        cerr << "# ref sequence: " << ref_fasta.data() << endl;
        if(command >= INDEX){
            //build index
            cerr << "building index with s = " << opt.K  << " ... "<< endl;
            clock_t start = clock();
            sa = new sparseSA(ref, output.refdescr, output.startpos, opt._4column, opt.K);
            clock_t end = clock();
            double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
            cerr << "building index: done" << endl;
            cerr << "time for building index structure: " << cpu_time << endl;

            if(command == INDEX) delete sa;
        }
        if(command == ALN){
            if(outputName.empty()) 
                outputName = ref_fasta.substr().append(".sam");
            //Print SAM Header, require refdescre, startpos and argc/argv
            output.createHeader(argc, argv, ref.length());
            outfile = fopen( outputName.c_str(), "w" );
            fprintf(outfile,"%s",output.header.c_str());
            //FIRST: Unpaired reads
            if(!unpairedQ.empty()){
                queryReader = new fastqInputReader(unpairedQ, opt.nucleotidesOnly);
                pthread_mutex_init(&writeLock_, NULL);
                pthread_attr_t attr;  pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                vector<query_arg> args(opt.query_threads);
                vector<pthread_t> thread_ids(opt.query_threads);
                cerr << "Mapping unpaired reads to the index using " << opt.query_threads << " threads ..." << endl;
                clock_t start = clock();
                // Initialize additional thread data.
                for(int i = 0; i < opt.query_threads; i++) {
                    args[i].skip = opt.query_threads;
                    args[i].skip0 = i;
                    args[i].opt = & opt;
                    args[i].readLock = &queryReader->readLock_;
                    args[i].writeLock = &writeLock_;
                }
                // Create joinable threads to find MEMs.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_create(&thread_ids[i], &attr, query_thread, (void *)&args[i]);
                // Wait for all threads to terminate.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_join(thread_ids[i], NULL);
                clock_t end = clock();
                double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
                cerr << "mapping unpaired: done" << endl;
                cerr << "time for mapping: " << cpu_time << endl;
                delete queryReader;
            }
            //SECOND: Paired reads
            if(!pair1.empty() && !pair2.empty()){
                mate1Reader = new fastqInputReader(pair1, opt.nucleotidesOnly);
                mate2Reader = new fastqInputReader(pair2, opt.nucleotidesOnly);
                pthread_attr_t attr;  pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                vector<query_arg> args(opt.query_threads);
                vector<pthread_t> thread_ids(opt.query_threads);
                cerr << "Mapping paired reads to the index using " << opt.query_threads << " threads ..." << endl;
                clock_t start = clock();
                // Initialize additional thread data.
                for(int i = 0; i < opt.query_threads; i++) {
                    args[i].skip = opt.query_threads;
                    args[i].skip0 = i;
                    args[i].opt = & opt;
                    args[i].readLock = &mate1Reader->readLock_;
                    args[i].writeLock = &writeLock_;
                }
                // Create joinable threads to find MEMs.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_create(&thread_ids[i], &attr, paired_thread, (void *)&args[i]);
                // Wait for all threads to terminate.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_join(thread_ids[i], NULL);
                clock_t end = clock();
                double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
                cerr << "mapping paired: done" << endl;
                cerr << "time for mapping: " << cpu_time << endl;
                delete mate1Reader;
                delete mate2Reader;
            }
            else if(!pair1.empty() || !pair2.empty()){
                fprintf(stderr, "ALFALFA requires 2 separate mate pairs files \n");
                exit(1);
            }
            fclose( outfile );
            delete sa;
            deleteDPMatrix(opt.alnOptions.scores.openGap != 0);
            cerr << "FINISHED" << endl;
        }
    }
    return 0;
}