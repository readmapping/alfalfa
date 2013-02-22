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
#include "performanceUtils.h"

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
static const string SAM_VERSION = "1.3";
static const string PROG_VERSION = "0.6.3.1";
static const string NOT_AVAILABLE = "*";
static const long INIT_DP_DIMENSION = 2048;

//output struct
struct samOutput{//construct
    string header;
    vector<string> refdescr;
    vector<long> startpos;
    void createHeader(int argc, char* argv[], long refLength){
        stringstream ss;
        ss << "@HD\tVN:" << SAM_VERSION << endl;
        for(size_t i=0; i < refdescr.size()-1; i++){
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

//currently not used
void *unpaired_thread(void *arg_) {

  query_arg *arg = (query_arg *)arg_;
  long seq_cnt = 0;
  long seq_mapped = 0;
  long alignments_printed = 0;
  bool hasRead = true;
  while(hasRead){
      read_t read;
      pthread_mutex_lock(arg->readLock);
      hasRead = queryReader->nextRead(read.qname,read.sequence,read.qual);
      pthread_mutex_unlock(arg->readLock);
      if(hasRead){
          seq_cnt++;
          if(seq_cnt%10000==0)
              cerr << ".";
          read.init(arg->opt->nucleotidesOnly);
          unpairedMatch(*sa, *arg->dp, read, arg->opt->alnOptions);
          read.postprocess(arg->opt->alnOptions.scores, *sa);
          //Ouput
          pthread_mutex_lock(arg->writeLock);
          stringstream * ss = new stringstream;
          if(read.alignments.empty()){
              fprintf(outfile,"%s",read.emptyAlingment(false,false,true).c_str());
          }
          else{
              seq_mapped++;
              alignments_printed += read.alignmentCount();
              for(int k = 0; k < read.alignmentCount(); k++)
                fprintf(outfile,"%s",read.printUnpairedAlignment(k).c_str());
          }
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences read by thread %d: %ld\n",arg->skip0, seq_cnt);
  printf("sequences mapped by thread %d: %ld\n",arg->skip0, seq_mapped);
  printf("alignments printed by thread %d: %ld\n",arg->skip0, alignments_printed);
  delete arg->dp;
  pthread_exit(NULL);
}

void *query_thread(void *arg_) {

  query_arg *arg = (query_arg *)arg_;
  long seq_cnt = 0;
  long seq_mapped = 0;
  long alignments_printed = 0;
  bool hasRead = true;
  while(hasRead){
      read_t read;
      pthread_mutex_lock(arg->readLock);
      hasRead = queryReader->nextRead(read.qname,read.sequence,read.qual);
      pthread_mutex_unlock(arg->readLock);
      if(hasRead){
          seq_cnt++;
          if(seq_cnt%10000==0)
              cerr << ".";
          if(!arg->opt->alnOptions.noFW){//to inner funtion
              if(arg->opt->alnOptions.print>0)
                cerr << "match " << read.qname << " forward" << endl;
              inexactMatch(*sa, *arg->dp, read, arg->opt->alnOptions, true);
          }
          if(!arg->opt->alnOptions.noRC){
              read.init(arg->opt->nucleotidesOnly);
              if(arg->opt->alnOptions.print>0)
                cerr << "match " << read.qname << " backward" << endl;
              inexactMatch(*sa, *arg->dp, read, arg->opt->alnOptions, false);
          }
          if(arg->opt->alnOptions.print>0) cerr << "post process read alignments" << endl;
          read.postprocess(arg->opt->alnOptions.scores, *sa);
          if(arg->opt->alnOptions.print>0) cerr << "read " << read.qname << " has " << read.alignmentCount() << " alignments" << endl;
          //Ouput
          pthread_mutex_lock(arg->writeLock);
          stringstream * ss = new stringstream;
          if(read.alignments.empty())
              fprintf(outfile,"%s",read.emptyAlingment(false,false,true).c_str());
          else{
              seq_mapped++;
              alignments_printed += read.alignmentCount();
              for(int k = 0; k < read.alignmentCount(); k++)
                fprintf(outfile,"%s",read.printUnpairedAlignment(k).c_str());
          }
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences read by thread %d: %ld\n",arg->skip0, seq_cnt);
  printf("sequences mapped by thread %d: %ld\n",arg->skip0, seq_mapped);
  printf("alignments printed by thread %d: %ld\n",arg->skip0, alignments_printed);
  delete arg->dp;
  pthread_exit(NULL);
}

void *paired_thread1(void *arg_) {

  query_arg *arg = (query_arg *)arg_;
  long seq_cnt = 0;
  long seq_mapped1 = 0;
  long alignments_printed1 = 0;
  long seq_mapped2 = 0;
  long alignments_printed2 = 0;
  bool hasRead = true;
  while(hasRead){
      read_t mate1;
      read_t mate2;
      pthread_mutex_lock(arg->readLock);
      hasRead = mate1Reader->nextRead(mate1.qname,mate1.sequence,mate1.qual);
      hasRead = mate2Reader->nextRead(mate2.qname,mate2.sequence,mate2.qual);
      pthread_mutex_unlock(arg->readLock);
      if(hasRead){
          if(mate1.qname.length() > 2 && mate1.qname[mate1.qname.length()-2]=='/')
              mate1.qname.erase(mate1.qname.length()-2);
          if(mate2.qname.length() > 2 && mate2.qname[mate2.qname.length()-2]=='/')
              mate2.qname.erase(mate2.qname.length()-2);
          seq_cnt++;
          if(seq_cnt%10000==0)
              cerr << ".";
          mate1.init(arg->opt->nucleotidesOnly);
          mate2.init(arg->opt->nucleotidesOnly);
          if(arg->opt->alnOptions.print>0) cerr << "match " << mate1.qname << endl;
          pairedMatch(*sa, *arg->dp, mate1, mate2, arg->opt->alnOptions, arg->opt->pairedOpt);
          if(arg->opt->alnOptions.print>0) cerr << "postprocess both mates of " << mate1.qname << endl;
          mate1.postprocess(arg->opt->alnOptions.scores, *sa);
          mate2.postprocess(arg->opt->alnOptions.scores, *sa);
          if(arg->opt->alnOptions.print>0) cerr << mate1.qname << " has " << mate1.alignmentCount() << " alignments" << endl;
          if(arg->opt->alnOptions.print>0) cerr << mate2.qname << " has " << mate2.alignmentCount() << " alignments" << endl;
          //Ouput
          pthread_mutex_lock(arg->writeLock);
          stringstream * ss = new stringstream;
          if(mate1.alignments.empty())
              fprintf(outfile,"%s",mate1.emptyAlingment(true,mate2.alignments.empty(),true).c_str());
          else{
              seq_mapped1++;
              alignments_printed1+= mate1.pairedAlignmentCount;
              for(int k = 0; k < mate1.alignmentCount(); k++)
                  if(mate1.alignments[k]->flag.test(0))
                    fprintf(outfile,"%s",mate1.printPairedAlignments(k).c_str());
          }
          if(mate2.alignments.empty())
              fprintf(outfile,"%s",mate2.emptyAlingment(true,mate1.alignments.empty(),false).c_str());
          else{
              seq_mapped2++;
              alignments_printed2+= mate2.pairedAlignmentCount;
              for(int k = 0; k < mate2.alignmentCount(); k++)
                  if(mate2.alignments[k]->flag.test(0))
                    fprintf(outfile,"%s",mate2.printPairedAlignments(k).c_str());
          }
          pthread_mutex_unlock(arg->writeLock);
          delete ss;
      }
  }
  printf("sequences read by thread %d: %ld\n", arg->skip0,seq_cnt);
  printf("sequences mapped (mate1) by thread %d: %ld\n", arg->skip0, seq_mapped1);
  printf("sequences mapped (mate2) by thread %d: %ld\n", arg->skip0, seq_mapped2);
  printf("alignments printed (mate1) by thread %d: %ld\n", arg->skip0, alignments_printed1);
  printf("alignments printed (mate2) by thread %d: %ld\n", arg->skip0, alignments_printed2);
  delete arg->dp;
  pthread_exit(NULL);
}

int main(int argc, char* argv[]){
    cerr << "@PG\tID:" << PROG << "\tVN:" << PROG_VERSION << endl;
    if(argc < 2 ){
        usage(PROG);
        exit(1);
    }
    if(strcmp(argv[1], "check") == 0){
        samCheckOptions_t opt;
        opt.initOptions();
        processCheckParameters(argc-1, argv+1, opt, PROG );
        cerr << "parsing options: done" << endl;
        if(opt.subcommand == ORACLE)
            checkOracle(opt);
        else if(opt.subcommand == SUMMARY)
            checkSummary(opt);
        else if(opt.subcommand == COMPARE)
            checkCompare(opt);
        else if(opt.subcommand == WGSIM)
            checkWgsim(opt);
    }
    else{
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
        if(opt.command == INDEX){
            //build index
            cerr << "building index with s = " << opt.K  << " ... "<< endl;
            clock_t start = clock();
            sa = new sparseSA(ref, output.refdescr, output.startpos, opt._4column, opt.K, opt.hasSuflink, opt.hasChild);
            sa->construct();
            cerr << "building index: done" << endl;
            if(opt.saveIndex){
                cerr << "saving index to disk" << " ... "<< endl;
                sa->save(opt.index_prefix);
                cerr << "saving index to disk: done" << endl;
            }
            clock_t end = clock();
            double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
            cerr << "time for building index structure: " << cpu_time << endl;

            delete sa;
        }
        else if(opt.command == ALN){
            //build or load index
            sa = new sparseSA(ref, output.refdescr, output.startpos, opt._4column, opt.K, opt.hasSuflink, opt.hasChild);
            if(opt.indexLocation.empty() || !sa->load(opt.indexLocation)){
                cerr << "building index with s = " << opt.K  << " ... "<< endl;
                //adjust data structure to seed choice and sparseness
                if((opt.alnOptions.memType == MAM || opt.alnOptions.memType == MUM) && opt.K>3){
                    opt.hasChild = sa->hasChild = false;
                    opt.hasSuflink = sa->hasSufLink = true;
                    cerr << "because of MAM or MUM seeds, ESSA has suffix links and no child array" << endl;
                }
                else if(opt.alnOptions.memType == SMAM){
                    if(opt.K >= 3){
                        opt.hasChild = sa->hasChild = true;
                        opt.hasSuflink = sa->hasSufLink = false;
                        cerr << "because of SMAM seeds and s= " << opt.K << ", ESSA has child array and no suffix links" << endl;
                    }
                    else{
                        opt.hasChild = sa->hasChild = true;
                        opt.hasSuflink = sa->hasSufLink = false;
                        cerr << "because of SMAM seeds and s= " << opt.K << ", ESSA has suffix links and no child array" << endl;
                    }
                }
                else if(opt.alnOptions.memType == MEM){
                    opt.hasChild = sa->hasChild = true;
                    opt.hasSuflink = sa->hasSufLink = false;
                    cerr << "because of MEM seeds, ESSA has child array and no suffix links" << endl;
                }
                clock_t start = clock();
                sa->construct();
                cerr << "building index: done" << endl;
                if(opt.saveIndex){  
                    cerr << "saving index to disk" << " ... "<< endl;
                    sa->save(opt.index_prefix);
                    cerr << "saving index to disk: done" << endl;
                }
                clock_t end = clock();
                double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
                cerr << "time for building index structure: " << cpu_time << endl;
            }
            //calculate skip parameter for MEMs
            if(opt.alnOptions.memType == MEM){
                if (opt.K >= 4) sa->sparseMult = (int) (opt.alnOptions.minMemLength - 10) / opt.K;
                else sa->sparseMult = (int) (opt.alnOptions.minMemLength - 12) / opt.K;
                opt.alnOptions.sparseMult = sa->sparseMult;
                if(opt.alnOptions.print) cerr << "skip factor was set to " << opt.alnOptions.sparseMult << endl;
            }
            if(opt.alnOptions.print>0){
                cerr << "index has sparseness " << opt.K << endl;
                cerr << "index uses suffix links? " << (sa->hasSufLink ? "yes" : "no") << endl;
                cerr << "index uses child array? " << (sa->hasChild ? "yes" : "no") << endl;
                cerr << "skip factor is set to " << sa->sparseMult << endl;
                cerr << "INDEX SIZE IN BYTES: " << sa->index_size_in_bytes() << endl;
            }
            if(opt.outputName.empty()){
                opt.outputName = opt.ref_fasta.substr().append(".sam");
                if(opt.alnOptions.print>0)
                    cerr << "output name changed to " << opt.outputName << endl;
            }
            //Print SAM Header, require refdescre, startpos and argc/argv
            output.createHeader(argc, argv, ref.length());
            outfile = fopen( opt.outputName.c_str(), "w" );
            fprintf(outfile,"%s",output.header.c_str());
            if(opt.alnOptions.print>0) cerr << "header printed to output" << endl;
            if(opt.alnOptions.print>0) cerr << "dp matrices will be created for every thread with initial " << INIT_DP_DIMENSION << "x" << INIT_DP_DIMENSION << " dimension" << endl;
            //FIRST: Unpaired reads
            if(!opt.unpairedQ.empty()){
                queryReader = new fastqInputReader(opt.unpairedQ, opt.nucleotidesOnly);
                pthread_mutex_init(&writeLock_, NULL);
                pthread_attr_t attr;  pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                vector<query_arg> args(opt.query_threads);
                vector<pthread_t> thread_ids(opt.query_threads);
                cerr << "Mapping unpaired reads to the index using " << opt.query_threads << " threads ..." << endl;
                cerr << "Progress (each dot represents 10.000 reads processed: ";
                clock_t start = clock();
                // Initialize additional thread data.
                for(int i = 0; i < opt.query_threads; i++) {
                    args[i].skip = opt.query_threads;
                    args[i].skip0 = i;
                    args[i].opt = & opt;
                    args[i].readLock = &queryReader->readLock_;
                    args[i].writeLock = &writeLock_;
                    args[i].dp = new dynProg(INIT_DP_DIMENSION, opt.alnOptions.scores.openGap!=0, opt.alnOptions.scores);
                }
                // Create joinable threads to find MEMs.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_create(&thread_ids[i], &attr, query_thread, (void *)&args[i]);
                // Wait for all threads to terminate.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_join(thread_ids[i], NULL);
                clock_t end = clock();
                double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
                cerr << endl;
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
                cerr << "Progress (each dot represents 10.000 reads processed: ";
                clock_t start = clock();
                // Initialize additional thread data.
                for(int i = 0; i < opt.query_threads; i++) {
                    args[i].skip = opt.query_threads;
                    args[i].skip0 = i;
                    args[i].opt = & opt;
                    args[i].readLock = &mate1Reader->readLock_;
                    args[i].writeLock = &writeLock_;
                    args[i].dp = new dynProg(INIT_DP_DIMENSION, opt.alnOptions.scores.openGap!=0, opt.alnOptions.scores);
                }
                // Create joinable threads to find MEMs.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_create(&thread_ids[i], &attr, paired_thread1, (void *)&args[i]);
                // Wait for all threads to terminate.
                for(int i = 0; i < opt.query_threads; i++)
                    pthread_join(thread_ids[i], NULL);
                clock_t end = clock();
                double cpu_time = (double)( end - start ) /CLOCKS_PER_SEC;
                cerr << endl;
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
    }
    return 0;
}

