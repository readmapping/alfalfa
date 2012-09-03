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

#include "fasta.h"
#include "mapper.h"

#include <time.h>
#include <stdio.h>
#include <sys/time.h>
#include <sstream>
#include <string>
#include <string.h>
#include <bitset>

#include <cctype> // std::tolower(), uppercase/lowercase conversion

// NOTE use of special characters ~, `, and $ !!!!!!!!

using namespace std;

//mapper options
static const string PROG = "ALFALFA";
static const string SAM_VERSION = "1.4";
static const string PROG_VERSION = "0.5.0";
static string NAN = "*";

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

struct query_arg {
    query_arg(): skip0(0), skip(0), opt(0){}
  int skip0;
  int skip;
  mapOptions_t * opt;
  dynProg * dp;
  pthread_mutex_t *readLock;
  pthread_mutex_t *writeLock;
};

void *unpaired_thread(void *arg_) {

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
          read.init(arg->opt->nucleotidesOnly);
          unpairedMatch(*sa, *arg->dp, read, arg->opt->alnOptions, print);
          read.postprocess(arg->opt->alnOptions.scores, *sa);
          //Ouput
          pthread_mutex_lock(arg->writeLock);
          stringstream * ss = new stringstream;
          if(read.alignments.empty())
              fprintf(outfile,"%s",read.emptyAlingment(false,false).c_str());
          else
              for(int k = 0; k < read.alignmentCount(); k++)
                fprintf(outfile,"%s",read.printAlignment(k).c_str());
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences mapped: %ld\n", seq_cnt);
  delete arg->dp;
  pthread_exit(NULL);
}

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
          if(!arg->opt->alnOptions.noFW)//to inner funtion
                inexactMatch(*sa, *arg->dp, read, arg->opt->alnOptions, true, print);
          if(!arg->opt->alnOptions.noRC){
              read.init(arg->opt->nucleotidesOnly);
              inexactMatch(*sa, *arg->dp, read, arg->opt->alnOptions, false, print);
          }
          read.postprocess(arg->opt->alnOptions.scores, *sa);
          //Ouput
          pthread_mutex_lock(arg->writeLock);
          stringstream * ss = new stringstream;
          if(read.alignments.empty())
              fprintf(outfile,"%s",read.emptyAlingment(false,false).c_str());
          else
              for(int k = 0; k < read.alignmentCount(); k++)
                fprintf(outfile,"%s",read.printAlignment(k).c_str());
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences mapped: %ld\n", seq_cnt);
  delete arg->dp;
  pthread_exit(NULL);
}

void *paired_thread1(void *arg_) {

  query_arg *arg = (query_arg *)arg_;
  bool print = arg->opt->verbose;
  long seq_cnt = 0;
  bool hasRead = true;
  while(hasRead){
      read_t mate1;
      read_t mate2;
      pthread_mutex_lock(arg->readLock);
      hasRead = mate1Reader->nextRead(mate1.qname,mate1.sequence,mate1.qual);
      hasRead = mate2Reader->nextRead(mate2.qname,mate2.sequence,mate2.qual);
      pthread_mutex_unlock(arg->readLock);
      if(hasRead){
          seq_cnt++;
          mate1.init(arg->opt->nucleotidesOnly);
          mate2.init(arg->opt->nucleotidesOnly);
          pairedMatch(*sa, *arg->dp, mate1, mate2, arg->opt->alnOptions, arg->opt->pairedOpt, print);
          mate1.postprocess(arg->opt->alnOptions.scores, *sa);
          mate2.postprocess(arg->opt->alnOptions.scores, *sa);
          //Ouput
          pthread_mutex_lock(arg->writeLock);
          stringstream * ss = new stringstream;
          if(mate1.alignments.empty())
              fprintf(outfile,"%s",mate1.emptyAlingment(true,mate2.alignments.empty()).c_str());
          else{
              for(int k = 0; k < mate1.alignmentCount(); k++)
                  if(mate1.alignments[k].flag.test(0))
                    fprintf(outfile,"%s",mate1.printAlignment(k).c_str());
          }
          if(mate2.alignments.empty())
              fprintf(outfile,"%s",mate2.emptyAlingment(true,mate1.alignments.empty()).c_str());
          else{
              for(int k = 0; k < mate2.alignmentCount(); k++)
                  if(mate2.alignments[k].flag.test(0))
                    fprintf(outfile,"%s",mate2.printAlignment(k).c_str());
          }
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences mapped: %ld\n", seq_cnt);
  delete arg->dp;
  pthread_exit(NULL);
}

int main(int argc, char* argv[]){
    cerr << "@PG\tID:" << PROG << "\tVN:" << PROG_VERSION << endl;
    mapOptions_t opt;
    opt.initOptions();
    processParameters(argc, argv, opt, PROG);
    //query name + read initialization
    cerr << "parsing options: done" << endl;
    cerr << "loading ref sequences: ..." << endl;
    string ref;
    samOutput output;
    load_fasta(opt.ref_fasta, ref, output.refdescr, output.startpos);
    
    cerr << "loading ref sequences: done" << endl;
    cerr << "# ref sequence: " << opt.ref_fasta.data() << endl;
    if(opt.command >= INDEX){
        //build index
        cerr << "building index with s = " << opt.K  << " ... "<< endl;
        clock_t start = clock();
        sa = new sparseSA(ref, output.refdescr, output.startpos, opt._4column, opt.K);
        clock_t end = clock();
        double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
        cerr << "building index: done" << endl;
        cerr << "time for building index structure: " << cpu_time << endl;

        if(opt.command == INDEX) delete sa;
    }
    if(opt.command == ALN){
        if(opt.outputName.empty()) 
            opt.outputName = opt.ref_fasta.substr().append(".sam");
        //Print SAM Header, require refdescre, startpos and argc/argv
        output.createHeader(argc, argv, ref.length());
        outfile = fopen( opt.outputName.c_str(), "w" );
        fprintf(outfile,"%s",output.header.c_str());
        //FIRST: Unpaired reads
        if(!opt.unpairedQ.empty()){
            queryReader = new fastqInputReader(opt.unpairedQ, opt.nucleotidesOnly);
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
                args[i].dp = new dynProg(2048, opt.alnOptions.scores.openGap!=0, opt.alnOptions.scores);
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
        if(!opt.pair1.empty() && !opt.pair2.empty()){
            mate1Reader = new fastqInputReader(opt.pair1, opt.nucleotidesOnly);
            mate2Reader = new fastqInputReader(opt.pair2, opt.nucleotidesOnly);
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
                args[i].dp = new dynProg(2048, opt.alnOptions.scores.openGap!=0, opt.alnOptions.scores);
            }
            // Create joinable threads to find MEMs.
            for(int i = 0; i < opt.query_threads; i++)
                pthread_create(&thread_ids[i], &attr, paired_thread1, (void *)&args[i]);
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
        else if(!opt.pair1.empty() || !opt.pair2.empty()){
            fprintf(stderr, "ALFALFA requires 2 separate mate pairs files \n");
            exit(1);
        }
        fclose( outfile );
        delete sa;
        cerr << "FINISHED" << endl;
    }
    return 0;
}