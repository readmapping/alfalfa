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

// NOTE use of special characters ~, `, and $ !!!!!!!!

//TODO: in makefile add gcc -DDEBUG
//use ifdef instead of if with command line parameter #ifdef DEBUG #endif

using namespace std;

//mapper options
static const string PROG = "alfalfa";
static const string SAM_VERSION = "1.3";
static const string PROG_VERSION = "0.8";
static const string NOT_AVAILABLE = "*";

//FIELDS
static sparseSA *sa;
static fastqInputReader *queryReader;
static FILE * outfile;
static pthread_mutex_t writeLock_;

static fastqInputReader *mate1Reader;
static fastqInputReader *mate2Reader;

struct query_arg {

    query_arg() : skip0(0), skip(0), opt(0) {
    }
    int skip0;
    int skip;
    mapOptions_t * opt;
    pthread_mutex_t *readLock;
    pthread_mutex_t *writeLock;
};

//currently not used

void *unpaired_thread(void *arg_) {

    query_arg *arg = (query_arg *) arg_;
    long seq_cnt = 0;
    long seq_mapped = 0;
    long alignments_printed = 0;
    bool hasRead = true;
#ifndef NDEBUG
    if (arg->opt->alnOptions.debugFile != NULL)
        fprintf(arg->opt->alnOptions.debugFile, "#ALN\n");
#endif
    while (hasRead) {
        read_t read;
        pthread_mutex_lock(arg->readLock);
        hasRead = queryReader->nextRead(read.qname, read.sequence, read.qual);
        pthread_mutex_unlock(arg->readLock);
        if (hasRead) {
            if (read.qname.length() > 2 && read.qname[read.qname.length() - 2] == '/')
                read.qname.erase(read.qname.length() - 2);
            seq_cnt++;
            if (seq_cnt % 10000 == 0)
                cerr << ".";
            read.init(arg->opt->nucleotidesOnly);
            if (arg->opt->alnOptions.print >= 3) cerr << "match " << read.qname << " with length " << read.sequence.length() << endl;
            if (arg->opt->alnOptions.print >= 2) cerr << "#READ" << endl;
#ifndef NDEBUG
            if (arg->opt->alnOptions.debugFile != NULL) {
                fprintf(arg->opt->alnOptions.debugFile, ">%s\n", read.qname.c_str());
                fprintf(arg->opt->alnOptions.debugFile, "%lu\n", read.sequence.size());
            }
#endif
            unpairedMatch(*sa, read, arg->opt->alnOptions);
            read.postprocess(*sa, arg->opt->alnOptions, false);
            //Ouput
            pthread_mutex_lock(arg->writeLock);
            if (read.alignments.empty())
                fprintf(outfile, "%s", read.emptyAlingment(false, false, true).c_str());
            else {
                seq_mapped++;
                alignments_printed += min(read.alignmentCount(), arg->opt->alnOptions.alignmentCount);
                for (int k = 0; k < min(read.alignmentCount(), arg->opt->alnOptions.alignmentCount); k++)
                    fprintf(outfile, "%s", read.printUnpairedAlignment(k).c_str());
            }
            pthread_mutex_unlock(arg->writeLock);
        }
    }
    cerr << endl;
    printf("sequences read by thread %d: %ld\n", arg->skip0, seq_cnt);
    printf("sequences mapped by thread %d: %ld\n", arg->skip0, seq_mapped);
    printf("alignments printed by thread %d: %ld\n", arg->skip0, alignments_printed);
    pthread_exit(NULL);
}

void *paired_thread1(void *arg_) {

    query_arg *arg = (query_arg *) arg_;
    long seq_cnt = 0;
    long seq_mapped1 = 0;
    long alignments_printed1 = 0;
    long seq_mapped2 = 0;
    long alignments_printed2 = 0;
    bool hasRead = true;
    while (hasRead) {
        read_t mate1;
        read_t mate2;
        pthread_mutex_lock(arg->readLock);
        hasRead = mate1Reader->nextRead(mate1.qname, mate1.sequence, mate1.qual);
        hasRead = mate2Reader->nextRead(mate2.qname, mate2.sequence, mate2.qual);
        pthread_mutex_unlock(arg->readLock);
        if (hasRead) {
            if (mate1.qname.length() > 2 && mate1.qname[mate1.qname.length() - 2] == '/')
                mate1.qname.erase(mate1.qname.length() - 2);
            if (mate2.qname.length() > 2 && mate2.qname[mate2.qname.length() - 2] == '/')
                mate2.qname.erase(mate2.qname.length() - 2);
            seq_cnt++;
            if (seq_cnt % 10000 == 0)
                cerr << ".";
            mate1.init(arg->opt->nucleotidesOnly);
            mate2.init(arg->opt->nucleotidesOnly);
            if (arg->opt->alnOptions.print >= 3) cerr << "match " << mate1.qname << " with lengths " << mate1.sequence.length() << " and " << mate2.sequence.length() << endl;
            if (arg->opt->alnOptions.print >= 2) cerr << "#READ" << endl;
            pairedMatch(*sa, mate1, mate2, arg->opt->alnOptions, arg->opt->pairedOpt);
            mate1.postprocess(*sa, arg->opt->alnOptions, true);
            mate2.postprocess(*sa, arg->opt->alnOptions, true);
            //Ouput
            pthread_mutex_lock(arg->writeLock);
            if (mate1.alignments.empty())
                fprintf(outfile, "%s", mate1.emptyAlingment(true, mate2.alignments.empty(), true).c_str());
            else {
                seq_mapped1++;
                int mate1AlnCount;
                if (arg->opt->pairedOpt.mixed)
                    mate1AlnCount = min(mate1.alignmentCount(), arg->opt->alnOptions.alignmentCount);
                else
                    mate1AlnCount = min(mate1.pairedAlignmentCount, arg->opt->alnOptions.alignmentCount);
                alignments_printed1 += mate1AlnCount;
                for (int k = 0; k < mate1AlnCount; k++)
                    fprintf(outfile, "%s", mate1.printPairedAlignments(k, true).c_str());
            }
            if (mate2.alignments.empty())
                fprintf(outfile, "%s", mate2.emptyAlingment(true, mate1.alignments.empty(), false).c_str());
            else {
                seq_mapped2++;
                int mate2AlnCount;
                if (arg->opt->pairedOpt.mixed)
                    mate2AlnCount = min(mate2.alignmentCount(), arg->opt->alnOptions.alignmentCount);
                else
                    mate2AlnCount = min(mate2.pairedAlignmentCount, arg->opt->alnOptions.alignmentCount);
                alignments_printed2 += mate2AlnCount;
                for (int k = 0; k < mate2AlnCount; k++)
                    fprintf(outfile, "%s", mate2.printPairedAlignments(k, false).c_str());
            }
            pthread_mutex_unlock(arg->writeLock);
        }
    }
    cerr << endl;
    printf("sequences read by thread %d: %ld\n", arg->skip0, seq_cnt);
    printf("sequences mapped (mate1) by thread %d: %ld\n", arg->skip0, seq_mapped1);
    printf("sequences mapped (mate2) by thread %d: %ld\n", arg->skip0, seq_mapped2);
    printf("alignments printed (mate1) by thread %d: %ld\n", arg->skip0, alignments_printed1);
    printf("alignments printed (mate2) by thread %d: %ld\n", arg->skip0, alignments_printed2);
    pthread_exit(NULL);
}

string createHeader(int argc, char* argv[], long refLength, sparseSA *sa) {
    stringstream ss;
    ss << "@HD\tVN:" << SAM_VERSION << endl;
    for (size_t i = 0; i < sa->descr.size() - 1; i++) {
        long length = sa->startpos[i + 1] - sa->startpos[i] - 1;
        ss << "@SQ\tSN:" << sa->descr[i] << "\tLN:" << length << endl; //startpos contains refSequences with '`' in between.
    }
    long length = refLength - sa->startpos[sa->descr.size() - 1];
    ss << "@SQ\tSN:" << sa->descr[sa->descr.size() - 1] << "\tLN:" << length << endl;
    ss << "@PG\tID:" << PROG << "\tCL:";
    for (int i = 0; i < argc; i++) {
        if (i > 0)
            ss << " ";
        ss << argv[i]; //add number to the stream
    }
    ss << "\tVN:" << PROG_VERSION << endl;
    return ss.str();
}

int main(int argc, char* argv[]) {
    cerr << "@PG\tID:" << PROG << "\tVN:" << PROG_VERSION << endl;
    //handle command structure and general help
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        option::printUsage(std::cerr, basicUsage);
        exit(1);
    }
    if (strcmp(argv[1], "evaluate") == 0) {
        samCheckOptions_t opt;
        opt.initOptions();
        processCheckParameters(argc - 1, argv + 1, opt);
        cerr << "parsing options: done" << endl;
        if (opt.subcommand == ORACLE)
            checkOracle(opt);
        else if (opt.subcommand == SUMMARY)
            checkSummary(opt);
        else if (opt.subcommand == WGSIM)
            checkWgsim(opt);
    } else {
        mapOptions_t opt;
        opt.initOptions();
        processParameters(argc, argv, opt);
        //query name + read initialization
        cerr << "parsing options: done" << endl;
        sa = new sparseSA();
        if (opt.command == INDEX) {
            cerr << "loading ref sequences: ..." << endl;
            sa->loadRef(opt.ref_fasta);
            cerr << "loading ref sequences: done" << endl;
            cerr << "# ref sequence: " << opt.ref_fasta.data() << endl;
            //build index
            cerr << "building index with s = " << opt.K << " ... " << endl;
            clock_t start = clock();
            sa->init(opt.K, opt.hasSuflink, opt.hasChild, opt.hasKmer);
            sa->construct();
            cerr << "building index: done" << endl;
            if (opt.saveIndex) {
                cerr << "saving index to disk" << " ... " << endl;
                sa->save(opt.index_prefix);
                cerr << "saving index to disk: done" << endl;
            }
            clock_t end = clock();
            double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
            cerr << "time for building index structure: " << cpu_time << endl;
            cerr << "INDEX SIZE IN BYTES: " << sa->index_size_in_bytes() << endl;
            delete sa;
        } else if (opt.command == ALN) {
            //build or load index
            if (opt.indexLocation.empty() || !sa->load(opt.indexLocation)) {
                sa->loadRef(opt.ref_fasta);
                sa->init(opt.K, opt.hasSuflink, opt.hasChild, opt.hasKmer); //because S needs to be changed
                //TODO: set choice for kmer: situation yes or no
                cerr << "building index with s = " << opt.K << " ... " << endl;
                //adjust data structure to seed choice and sparseness
                clock_t start = clock();
                sa->construct();
                cerr << "building index: done" << endl;
                if (opt.saveIndex) {
                    cerr << "saving index to disk" << " ... " << endl;
                    sa->save(opt.index_prefix);
                    cerr << "saving index to disk: done" << endl;
                }
                clock_t end = clock();
                double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
                cerr << "time for building index structure: " << cpu_time << endl;
            } else {
                opt.K = sa->K;
                opt.hasChild = sa->hasChild;
                opt.hasSuflink = sa->hasSufLink;
                opt.hasKmer = sa->hasKmer;
            }
            if (opt.alnOptions.print >= 1) {
                cerr << "index has sparseness " << opt.K << endl;
                cerr << "index uses suffix links? " << (sa->hasSufLink ? "yes" : "no") << endl;
                cerr << "index uses child array? " << (sa->hasChild ? "yes" : "no") << endl;
                cerr << "INDEX SIZE IN BYTES: " << sa->index_size_in_bytes() << endl;
            }
#ifndef NDEBUG
            if (opt.alnOptions.debugFile != NULL) {
                fprintf(opt.alnOptions.debugFile, "#REF %s\n", opt.ref_fasta.c_str());
                fprintf(opt.alnOptions.debugFile, "%lu\n", sa->S.length());
                long sequenceCount = sa->descr.size();
                fprintf(opt.alnOptions.debugFile, "%ld\n", sequenceCount);
                if (sequenceCount == 1)
                    fprintf(opt.alnOptions.debugFile, "%d %s %lu %ld\n", 0, sa->descr[0].c_str(), sa->S.length(), sa->startpos[0]);
                else {
                    for (long i = 0; i < sequenceCount - 1; i++)
                        fprintf(opt.alnOptions.debugFile, "%ld %s %ld %ld\n", i, sa->descr[i].c_str(), sa->startpos[i + 1] - sa->startpos[i], sa->startpos[i]);
                    fprintf(opt.alnOptions.debugFile, "%ld %s %ld %ld\n", sequenceCount - 1, sa->descr[sequenceCount - 1].c_str(), (long) (sa->S.length() - sa->startpos[sequenceCount - 1]), sa->startpos[sequenceCount - 1]);
                }
            }
#endif
            //Print SAM Header, require refdescre, startpos and argc/argv
            outfile = fopen(opt.outputName.c_str(), "w");
            string header = createHeader(argc, argv, sa->S.length(), sa);
            fprintf(outfile, "%s", header.c_str());
            if (opt.alnOptions.print >= 1) cerr << "header printed to output" << endl;
            //FIRST: Unpaired reads
            if (!opt.unpairedQ.empty()) {
                if (opt.alnOptions.print >= 2) cerr << "#SE" << endl;
                queryReader = new fastqInputReader(opt.unpairedQ, opt.nucleotidesOnly);
                pthread_mutex_init(&writeLock_, NULL);
                pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                vector<query_arg> args(opt.query_threads);
                vector<pthread_t> thread_ids(opt.query_threads);
                cerr << "Mapping unpaired reads to the index using " << opt.query_threads << " threads ..." << endl;
                cerr << "Progress (each dot represents 10.000 reads processed): " << endl;
                clock_t start = clock();
                // Initialize additional thread data.
                for (int i = 0; i < opt.query_threads; i++) {
                    args[i].skip = opt.query_threads;
                    args[i].skip0 = i;
                    args[i].opt = &opt;
                    args[i].readLock = &queryReader->readLock_;
                    args[i].writeLock = &writeLock_;
                }
                // Create joinable threads to find MEMs.
                for (int i = 0; i < opt.query_threads; i++)
                    pthread_create(&thread_ids[i], &attr, unpaired_thread, (void *) &args[i]);
                // Wait for all threads to terminate.
                for (int i = 0; i < opt.query_threads; i++)
                    pthread_join(thread_ids[i], NULL);
                clock_t end = clock();
                double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
                cerr << endl;
                cerr << "mapping unpaired: done" << endl;
                cerr << "time for mapping: " << cpu_time << endl;
                delete queryReader;
            }
            //SECOND: Paired reads
            if (!opt.pair1.empty() && !opt.pair2.empty()) {
                if (opt.alnOptions.print >= 2) cerr << "#PE" << endl;
                mate1Reader = new fastqInputReader(opt.pair1, opt.nucleotidesOnly);
                mate2Reader = new fastqInputReader(opt.pair2, opt.nucleotidesOnly);
                pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
                vector<query_arg> args(opt.query_threads);
                vector<pthread_t> thread_ids(opt.query_threads);
                cerr << "Mapping paired reads to the index using " << opt.query_threads << " threads ..." << endl;
                cerr << "Progress (each dot represents 10.000 reads processed): ";
                clock_t start = clock();
                // Initialize additional thread data.
                for (int i = 0; i < opt.query_threads; i++) {
                    args[i].skip = opt.query_threads;
                    args[i].skip0 = i;
                    args[i].opt = &opt;
                    args[i].readLock = &mate1Reader->readLock_;
                    args[i].writeLock = &writeLock_;
                }
                // Create joinable threads to find MEMs.
                for (int i = 0; i < opt.query_threads; i++)
                    pthread_create(&thread_ids[i], &attr, paired_thread1, (void *) &args[i]);
                // Wait for all threads to terminate.
                for (int i = 0; i < opt.query_threads; i++)
                    pthread_join(thread_ids[i], NULL);
                clock_t end = clock();
                double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
                cerr << endl;
                cerr << "mapping paired: done" << endl;
                cerr << "time for mapping: " << cpu_time << endl;
                delete mate1Reader;
                delete mate2Reader;
            } else if (!opt.pair1.empty() || !opt.pair2.empty()) {
                fprintf(stderr, "ALFALFA requires 2 separate mate pairs files \n");
                exit(1);
            }
            fclose(outfile);
#ifndef NDEBUG
            if (opt.alnOptions.debugFile != NULL) fclose(opt.alnOptions.debugFile);
#endif
            delete sa;
            cerr << "FINISHED" << endl;
        }
    }
    return 0;
}

