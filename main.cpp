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

void getlijn(ifstream & input, string & lijn){
    //write windows version of getline
    getline(input, lijn, '\n'); // Load one line at a time.
    if(lijn.length() > 0 && lijn[lijn.length()-1]=='\r')
        lijn.erase(--lijn.end());
}

enum mum_t { MUM, MAM, MEM, SMAM };
enum command_t {INDEX, MATCHES, ALN};
//mapper options

int MINOPTIONCOUNT = 2;
static const string PROG = "ALFALFA";
static const string SAM_VERSION = "1.4";
static const string PROG_VERSION = "0.3.1";
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
    }
    //performance options
    int K;//sparsity factor
    int query_threads;
    //I/O options
    bool _4column;//for MEM output
    bool verbose;
    string query_fast;
    string ref_fasta;
    string outputName;
    //Sequence options
    bool nucleotidesOnly;
    //MAM options
    enum mum_t memType;
    //alignment options
    align_opt alnOptions;
    bool noFW; bool noRC;
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

struct query_arg {
    query_arg(): alignments(0), skip0(0), skip(0), opt(0){}
  int skip0;
  int skip;
  mapOptions_t * opt;
  vector<read_t> alignments;
};

void *query_thread(void *arg_) {

  query_arg *arg = (query_arg *)arg_;

  string meta, line = NAN;
  ifstream data(arg->opt->query_fast.c_str());
  bool print = arg->opt->verbose;

  long seq_cnt = 0;

  if(!data.is_open()) { cerr << "unable to open " << arg->opt->query_fast << endl; exit(1); }

  bool fastq = true;
  string *P = new string;
  string *qual = new string;

  while(!data.eof() && fastq) {//first read which serves as swich between fasta and fastq
    getlijn(data, line); // Load one line at a time.
    if(line.length() == 0) continue;
    if(line[0] == '@') {
      long start = 1, end = line.length() - 1;
      trim(line, start, end);
      meta = line.substr(start,end-start+1);
      getlijn(data, line); //sequence line
      start = 0; end = line.length() - 1;
      trim(line, start,end);
      for(long i = start; i <= end; i++) {
        char c = std::tolower(line[i]);
        if(arg->opt->nucleotidesOnly) {
            switch(c) {
                case 'a': case 't': case 'g': case 'c': break;
                default:
                    c = '~';
            }
        }
        *P += c;
      }
      getlijn(data, line); //'+' line
      getlijn(data, line); //qual line
      if(seq_cnt % arg->skip == arg->skip0) {
          *qual = line;
          read_t read(meta,*P,*qual);
          if(!arg->opt->noFW)
            sa->inexactMatch(read, arg->opt->alnOptions, true, print);
          if(!arg->opt->noRC)
            sa->inexactMatch(read, arg->opt->alnOptions, false, print);
          read.postprocess(arg->opt->alnOptions.scores);
          arg->alignments.push_back(read);
      }
      seq_cnt++;
      delete qual; qual = new string;
      delete P; P = new string; meta = "";
      break;
    }
    else if(line[0] == '>'){
        long start = 1, end = line.length() - 1;
        trim(line, start, end);
        meta = line.substr(start,end-start+1);
        fastq = false;
    }
  }
  if(fastq){//fastq
    while(!data.eof()) {// Load one line at a time: cycle through meta, sequence, '+' and quality
        getlijn(data, line); // new meta
        if(line.length() == 0) continue;
        long start = 1, end = line.length()-1;
        trim(line, start , end);
        meta = line.substr(start,end-start+1);
        getlijn(data, line); //sequence line
        start = 0; end = line.length() - 1;
        trim(line, start,end);
        for(long i = start; i <= end; i++) {
            char c = std::tolower(line[i]);
            if(arg->opt->nucleotidesOnly) {
                switch(c) {
                    case 'a': case 't': case 'g': case 'c': break;
                    default:
                        c = '~';
                }
            }
            *P += c;
        }
        getlijn(data, line); //'+' line
        getlijn(data, line); //qual line
        if(seq_cnt % arg->skip == arg->skip0) {
            *qual = line;
            read_t read(meta,*P,*qual);
            if(!arg->opt->noFW)
                sa->inexactMatch(read, arg->opt->alnOptions, true, print);
            if(!arg->opt->noRC)
                sa->inexactMatch(read, arg->opt->alnOptions, false, print);
            read.postprocess(arg->opt->alnOptions.scores);
            arg->alignments.push_back(read);
        }
        seq_cnt++;
        delete qual; qual = new string;
        delete P; P = new string; meta = "";
    }
  }
  else{//fasta
        while(!data.eof()) {
        getlijn(data, line); // Load one line at a time.
        if(line.length() == 0) continue;
        long start = 0, end = line.length() - 1;
        // Meta tag line and start of a new sequence.
        // Collect meta data.
        if(line[0] == '>') {
          if(meta != "") {
            if(seq_cnt % arg->skip == arg->skip0) {
              read_t read(meta,*P,NAN);
              if(!arg->opt->noFW)
                sa->inexactMatch(read, arg->opt->alnOptions, true, print);
              if(!arg->opt->noRC)
                sa->inexactMatch(read, arg->opt->alnOptions, false, print);
              read.postprocess(arg->opt->alnOptions.scores);
              arg->alignments.push_back(read);
            }
            seq_cnt++;
            delete P; P = new string; meta = "";
          }
          start = 1;
          trim(line, start, end);
          meta = line.substr(start,end-start+1);
        }
        else { // Collect sequence data.
          trim(line, start,end);
          for(long i = start; i <= end; i++) {
            char c = std::tolower(line[i]);
            if(arg->opt->nucleotidesOnly) {
              switch(c) {
              case 'a': case 't': case 'g': case 'c': break;
              default:
                c = '~';
              }
            }
            *P += c;
          }
        }
      }
      // Handle very last sequence.
      if(meta != "") {
        if(seq_cnt % arg->skip == arg->skip0) {
          read_t read(meta,*P,NAN);
          if(!arg->opt->noFW)
            sa->inexactMatch(read, arg->opt->alnOptions, true, print);
          if(!arg->opt->noRC)
            sa->inexactMatch(read, arg->opt->alnOptions, false, print);
          read.postprocess(arg->opt->alnOptions.scores);
          arg->alignments.push_back(read);
        }
      }
  }
  delete P;
  delete qual;
  pthread_exit(NULL);
}

static const char * short_options = "s:k:l:frnt:q:d:vm:u:o:e:cC:T:Hh";
static struct option long_options[] = {
    {(char*)"sparsityfactor",   required_argument, 0,            's'},
    {(char*)"threads",          required_argument, 0,            'q'},
    {(char*)"verbose",          no_argument,       0,            'v'},
    {(char*)"seedminlength",    required_argument, 0,            'l'},
    {(char*)"alignments",       required_argument, 0,            'k'},
    {(char*)"trials",           required_argument, 0,            'T'},
    {(char*)"mincoverage",      required_argument, 0,            'C'},
    {(char*)"errors",           required_argument, 0,            'd'},
    {(char*)"tryharder",        no_argument,       0,            'H'},
    {(char*)"seedthreads",      required_argument, 0,            't'},
    {(char*)"noFw",             no_argument,       0,            'f'},
    {(char*)"noRc",             no_argument,       0,            'r'},
    {(char*)"wildcards",        no_argument,       0,            'N'},
    {(char*)"softclipping",     no_argument,       0,            'c'},
    {(char*)"match",            required_argument, 0,            'm'},
    {(char*)"mismatch",         required_argument, 0,            'u'},
    {(char*)"gapopen",          required_argument, 0,            'o'},
    {(char*)"gapextend",        required_argument, 0,            'e'},
    {(char*)"help",             no_argument,       0,            'h'},
    {(char*)0, 0, 0, 0} // terminator
};

void usage(string prog) {
  cerr << "Usage: " << prog << " COMMAND [options] <reference-file> <query-file> [ouput-file]" << endl;
  cerr << "Command should be one of the following: " << endl;
  cerr << "index                      only build the index for given <reference-file>, used for time calculations" << endl;
  cerr << "aln                        map the reads from <query-file> to the index build for <reference-file>" << endl;
  cerr << "the options are ordered by functionality: " << endl;
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
  cerr << "-H/--tryharder             enable: 'try harder': when no seeds have been found, search using less stringent parameters" << endl;
  cerr << "-t/--seedthreads (int)     number of threads for calculating the seeds [1]" << endl;
  cerr << "-f/--noFw                  do not compute forward matches" << endl;
  cerr << "-r/--noRc                  do not compute reverse complement matches" << endl;
  cerr << "-n/--wildcards             treat Ns as wildcard characters" << endl;
  cerr << "-c/--softclipping          allow soft clipping at the beginning and end of an alignment" << endl;
  cerr << endl;
  cerr << "DYNAMIC PROGRAMMING OPTIONS " << endl;
  cerr << "-m/--match (int)           match bonus [0]" << endl;
  cerr << "-u/--mismatch (int)        mismatch penalty [-2]" << endl;
  cerr << "-o/--gapopen (int)         gap open penalty (set 0 for non-affine gap-penalties) [0]" << endl;
  cerr << "-e/--gapextend (int)       gap extension penalty [-2]" << endl;
  cerr << endl;
  cerr << "MISC OPTIONS " << endl;
  cerr << "-v/--verbose               enable verbose mode (not by default)" << endl;
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
		case 's': opt.K = atoi(optarg); break;
                case 'k': opt.alnOptions.alignmentCount = atoi(optarg); break;
                case 'l': opt.alnOptions.minMemLength = atoi(optarg); opt.alnOptions.fixedMinLength = true; break;
                case 'f': opt.noFW = 1; break;
		case 'r': opt.noRC = 1; break;
                case 'n': opt.nucleotidesOnly = 1; break;
                case 'v': opt.verbose = 1; break;
                case 'c': opt.alnOptions.noClipping = false; break;
                case 't': opt.alnOptions.numThreads = atoi(optarg); break;
                case 'q': opt.query_threads = atoi(optarg); break;
                case 'm': opt.alnOptions.scores.match = atoi(optarg); break;
                case 'u': opt.alnOptions.scores.mismatch = atoi(optarg); break;
                case 'o': opt.alnOptions.scores.openGap = atoi(optarg); break;
                case 'e': opt.alnOptions.scores.extendGap = atoi(optarg); break;
		case 'd': opt.alnOptions.errorPercent = atof(optarg); break;
                case 'T': opt.alnOptions.maxTrial = atoi(optarg); break;
                case 'C': opt.alnOptions.minCoverage = atoi(optarg); break;
                case 'H': opt.alnOptions.tryHarder = true; break;
                case 'h': usage(PROG); break;
                case -1: /* Done with options. */
				break;
                case 0: if (long_options[option_index].flag != 0) break;
                default: usage(PROG);
				throw 1;
            }
	}
        if(opt.alnOptions.alignmentCount < 1){
            opt.alnOptions.alignmentCount = 1;
            opt.alnOptions.unique = true;
        }
        //add the reference query and output files
        if(optind + 1 < argc) opt.ref_fasta = argv[optind+1];
        else{
            fprintf(stderr, "ALFALFA requires input reference file \n");
            exit(1);
        }
        cerr << "parsing options: done" << endl;
        cerr << "loading ref sequences: ..." << endl;
        string ref;
        load_fasta(opt.ref_fasta, ref, output.refdescr, output.startpos);
        cerr << "loading ref sequences: done" << endl;
        cerr << "# ref sequence: " << opt.ref_fasta.data() << endl;
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
            //Print SAM Header, require refdescre, startpos and argc/argv
            opt.query_fast = optind+2 < argc ? argv[optind+2] : "";
            opt.outputName = optind+3 < argc ? argv[optind+3] : opt.query_fast.substr().append(".sam");
            output.createHeader(argc, argv, ref.length());
            pthread_attr_t attr;  pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            vector<query_arg> args(opt.query_threads);
            vector<pthread_t> thread_ids(opt.query_threads);
            cerr << "Mapping reads to the index using " << opt.query_threads << " threads ..." << endl;
            clock_t start = clock();
            // Initialize additional thread data.
            for(int i = 0; i < opt.query_threads; i++) {
                args[i].skip = opt.query_threads;
                args[i].skip0 = i;
                args[i].opt = & opt;
            }
            // Create joinable threads to find MEMs.
            for(int i = 0; i < opt.query_threads; i++)
                pthread_create(&thread_ids[i], &attr, query_thread, (void *)&args[i]);
            // Wait for all threads to terminate.
            for(int i = 0; i < opt.query_threads; i++)
                pthread_join(thread_ids[i], NULL);
            clock_t end = clock();
            double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
            cerr << "mapping: done" << endl;
            cerr << "time for mapping: " << cpu_time << endl;
            delete sa;
            cerr << "generating SAM and writing to " << opt.outputName << endl;
            FILE * outfile = fopen( opt.outputName.c_str(), "w" );  // open "shoppingList
            fprintf(outfile,"%s",output.header.c_str());
            //From global to local pos and write to stringstream
            for(int i = 0; i < opt.query_threads; i++) {
                for(int j = 0; j < args[i].alignments.size(); j++){
                    read_t & read = args[i].alignments[j];
                    if(read.alignments.empty())
                        fprintf(outfile,"%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\n",
                                read.qname.c_str(),4,"*",0,0,"*","*",0,0,read.sequence.c_str(),read.qual.c_str());
                    else{
                        long globPos;
                        long it;
                        string revCompl = read.sequence;
                        reverse_complement(revCompl, false);
                        string qualRC = read.qual;
                        reverse(qualRC.begin(),qualRC.end());
                        for(int k = 0; k < read.alignments.size(); k++){
                            alignment_t & a = read.alignments[k];
                            globPos = a.pos;
                            it = 1;
                            while(it < output.startpos.size() && globPos > output.startpos[it]){
                                it++;
                            }
                            a.pos = globPos - output.startpos[it-1];
                            a.rname = output.refdescr[it-1];
                            fprintf(outfile,"%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\tAS:i:%d\tNM:i:%d\tX0:Z:%s\n",
                        read.qname.c_str(),(int)a.flag.to_ulong(),a.rname.c_str(),a.pos,a.mapq,a.cigar.c_str(),
                        a.rnext.c_str(),a.pnext,a.tLength, a.flag.test(4) ? revCompl.c_str() : read.sequence.c_str(), a.flag.test(4) ? qualRC.c_str() : read.qual.c_str(),
                        a.alignmentScore,a.editDist,a.NMtag.c_str());
                        }
                    }
                }
            }
            fclose( outfile );
            cerr << "FINISHED" << endl;
        }
    }
  return 0;
}