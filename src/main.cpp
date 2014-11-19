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

#include "kthread.h"

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
static const string PROG_VERSION = "0.8.1";
static const string NOT_AVAILABLE = "*";
static const int BUFFER_SIZE = 10000000;

//FIELDS
static sparseSA *sa;
static fastqInputReader *queryReader;
static FILE * outfile;

static fastqInputReader *mate1Reader;
static fastqInputReader *mate2Reader;

int readBatch(int bpCount, vector<read_t> & reads){
    int curCount = 0;
    bool hasRead = true;
    while(curCount < bpCount && hasRead){
        read_t read;
        hasRead = queryReader->nextRead(read.qname, read.sequence, read.qual);
        if(hasRead){
            curCount += read.sequence.length();
            reads.push_back(read);
        }
    }
    return reads.size();
}

struct unpaired_data {
    unpaired_data(const mapOptions_t * opt_, const sparseSA * ssa_, 
        vector<read_t> * reads_) : opt(opt_), ssa(ssa_), reads(reads_) {
    }
    const mapOptions_t * opt;
    const sparseSA * ssa;
    vector<read_t> * reads; 
};

static void map_unpaired(void * unpairedData, int i, int threadid)
{
    unpaired_data * data = (unpaired_data*) unpairedData;
    read_t & read = data->reads->at(i);
    if (read.qname.length() > 2 && read.qname[read.qname.length() - 2] == '/')
        read.qname.erase(read.qname.length() - 2);
    read.init(data->opt->nucleotidesOnly);
    if (data->opt->alnOptions.print >= 3) cerr << "match " << read.qname << " with length " << read.sequence.length() << endl;
    if (data->opt->alnOptions.print >= 2) cerr << "#READ" << endl;
#ifndef NDEBUG
    if (data->opt->alnOptions.debugFile != NULL && data->opt->query_threads == 1) {
        fprintf(data->opt->alnOptions.debugFile, ">%s\n", read.qname.c_str());
        fprintf(data->opt->alnOptions.debugFile, "%lu\n", read.sequence.size());
    }
#endif
    unpairedMatch(*data->ssa, read, data->opt->alnOptions);
    read.postprocess(*data->ssa, data->opt->alnOptions, false);
}

void unpaired_thread(const mapOptions_t * opt) {
    long total_seq_cnt = 0;
    long seq_mapped = 0;
    long alignments_printed = 0;
#ifndef NDEBUG
    if (opt->alnOptions.debugFile != NULL)
        fprintf(opt->alnOptions.debugFile, "#ALN\n");
#endif
    //init fields
    vector<read_t> reads;
    int bpCount = opt->query_threads * BUFFER_SIZE;
    //read reads
    while(readBatch(bpCount, reads) > 0){// have reads
        //count reads, init fields
        int seq_cnt = reads.size();
        total_seq_cnt += seq_cnt;
        cerr << "Progress " << seq_cnt << " reads ..." << endl;
        clock_t start = clock();
        unpaired_data unpairedData(opt, sa, &reads);
        //parallel loop
        kt_for(opt->query_threads, seq_cnt, map_unpaired, &unpairedData);
        //output
        for(int i = 0; i < reads.size(); i++){
            if (reads[i].alignments.empty())
                fprintf(outfile, "%s", reads[i].emptyAlingment(false, false, true).c_str());
            else{
                seq_mapped++;
                alignments_printed += min(reads[i].alignmentCount(), opt->alnOptions.alignmentCount);
                for (int k = 0; k < min(reads[i].alignmentCount(), opt->alnOptions.alignmentCount); k++)
                    fprintf(outfile, "%s", reads[i].printUnpairedAlignment(k).c_str());
            }
        }
        reads.clear();
        clock_t end = clock();
        double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
        cerr << "Progress " << seq_cnt << " reads: done in " << cpu_time << " seconds" << endl;
    }
    //clean up
    reads.clear();
    //report statistics
    printf("sequences read: %ld\n", total_seq_cnt);
    printf("sequences mapped: %ld\n", seq_mapped);
    printf("alignments printed: %ld\n", alignments_printed);
}

int readBatchPaired(int bpCount, vector<read_t> & mates1, vector<read_t> & mates2){
    int curCount = 0;
    bool hasRead = true;
    while(curCount < bpCount && hasRead){
        read_t mate1;
        read_t mate2;
        hasRead = mate1Reader->nextRead(mate1.qname, mate1.sequence, mate1.qual);
        hasRead = mate2Reader->nextRead(mate2.qname, mate2.sequence, mate2.qual);
        if(hasRead){
            curCount += mate1.sequence.length();
            mates1.push_back(mate1);
            mates2.push_back(mate2);
        }
    }
    return mates1.size();
}

struct paired_data {
    paired_data(const mapOptions_t * opt_, const sparseSA * ssa_, 
        vector<read_t> * mates1_, vector<read_t> * mates2_) : opt(opt_), 
        ssa(ssa_), mates1(mates1_), mates2(mates2_) {
    }
    const mapOptions_t * opt;
    const sparseSA * ssa;
    vector<read_t> * mates1; 
    vector<read_t> * mates2; 
};

static void map_paired(void * pairedData, int i, int threadid)
{
    paired_data * data = (paired_data*) pairedData;
    read_t & mate1 = data->mates1->at(i);
    read_t & mate2 = data->mates2->at(i);
    if (mate1.qname.length() > 2 && mate1.qname[mate1.qname.length() - 2] == '/')
        mate1.qname.erase(mate1.qname.length() - 2);
    if (mate2.qname.length() > 2 && mate2.qname[mate2.qname.length() - 2] == '/')
        mate2.qname.erase(mate2.qname.length() - 2);
    mate1.init(data->opt->nucleotidesOnly);
    mate2.init(data->opt->nucleotidesOnly);
    if (data->opt->alnOptions.print >= 3) cerr << "match " << mate1.qname << " with lengths " << mate1.sequence.length() << " and " << mate2.sequence.length() << endl;
    if (data->opt->alnOptions.print >= 2) cerr << "#READ" << endl;
    pairedMatch(*data->ssa, mate1, mate2, data->opt->alnOptions, data->opt->pairedOpt);
    mate1.postprocess(*data->ssa, data->opt->alnOptions, true);
    mate2.postprocess(*data->ssa, data->opt->alnOptions, true);
}

void paired_thread1(const mapOptions_t * opt) {
    long total_seq_cnt = 0;
    long seq_mapped1 = 0;
    long alignments_printed1 = 0;
    long seq_mapped2 = 0;
    long alignments_printed2 = 0;
#ifndef NDEBUG
    if (opt->alnOptions.debugFile != NULL)
        fprintf(opt->alnOptions.debugFile, "#ALN\n");
#endif
    //init fields
    vector<read_t> mates1;
    vector<read_t> mates2;
    int bpCount = opt->query_threads * BUFFER_SIZE;
    //read reads
    while (readBatchPaired(bpCount, mates1, mates2) > 0) {
        //count reads, init fields
        int seq_cnt = mates1.size();
        total_seq_cnt += seq_cnt;
        cerr << "Progress " << seq_cnt << " reads ..." << endl;
        clock_t start = clock();
        paired_data pairedData(opt, sa, &mates1, &mates2);
        //parallel loop
        kt_for(opt->query_threads, seq_cnt, map_paired, &pairedData);
        //output
        for(int i = 0; i < mates1.size(); i++){
            if (mates1[i].alignments.empty())
                fprintf(outfile, "%s", mates1[i].emptyAlingment(true, mates2[i].alignments.empty(), true).c_str());
            else{
                seq_mapped1++;
                int mate1AlnCount;
                if (opt->pairedOpt.mixed)
                    mate1AlnCount = min(mates1[i].alignmentCount(), opt->alnOptions.alignmentCount);
                else
                    mate1AlnCount = min(mates1[i].pairedAlignmentCount, opt->alnOptions.alignmentCount);
                alignments_printed1 += mate1AlnCount;
                for (int k = 0; k < mate1AlnCount; k++)
                    fprintf(outfile, "%s", mates1[i].printPairedAlignments(k, true).c_str());
            }
            if (mates2[i].alignments.empty())
                fprintf(outfile, "%s", mates2[i].emptyAlingment(true, mates1[i].alignments.empty(), false).c_str());
            else {
                seq_mapped2++;
                int mate2AlnCount;
                if (opt->pairedOpt.mixed)
                    mate2AlnCount = min(mates2[i].alignmentCount(), opt->alnOptions.alignmentCount);
                else
                    mate2AlnCount = min(mates2[i].pairedAlignmentCount, opt->alnOptions.alignmentCount);
                alignments_printed2 += mate2AlnCount;
                for (int k = 0; k < mate2AlnCount; k++)
                    fprintf(outfile, "%s", mates2[i].printPairedAlignments(k, false).c_str());
            }
        }
        mates1.clear(); mates2.clear();
        clock_t end = clock();
        double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
        cerr << "Progress " << seq_cnt << " reads: done in " << cpu_time << " seconds" << endl;
    }
    //clean up
    mates1.clear();
    mates2.clear();
    //report statistics
    printf("sequences read: %ld\n", total_seq_cnt);
    printf("sequences mapped (mate1): %ld\n", seq_mapped1);
    printf("sequences mapped (mate2): %ld\n", seq_mapped2);
    printf("alignments printed (mate1): %ld\n", alignments_printed1);
    printf("alignments printed (mate2): %ld\n", alignments_printed2);
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
                cerr << "Mapping unpaired reads to the index using " << opt.query_threads << " threads ..." << endl;
                queryReader = new fastqInputReader(opt.unpairedQ, opt.nucleotidesOnly);
                clock_t start = clock();
                unpaired_thread(&opt);
                clock_t end = clock();
                double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
                cerr << endl;
                cerr << "mapping unpaired: done in " << cpu_time << " seconds." << endl;
                delete queryReader;
            }
            //SECOND: Paired reads
            if (!opt.pair1.empty() && !opt.pair2.empty()) {
                if (opt.alnOptions.print >= 2) cerr << "#PE" << endl;
                cerr << "Mapping paired reads to the index using " << opt.query_threads << " threads ..." << endl;
                mate1Reader = new fastqInputReader(opt.pair1, opt.nucleotidesOnly);
                mate2Reader = new fastqInputReader(opt.pair2, opt.nucleotidesOnly);
                clock_t start = clock();
                paired_thread1(&opt);
                clock_t end = clock();
                double cpu_time = (double) (end - start) / CLOCKS_PER_SEC;
                cerr << endl;
                cerr << "mapping paired: done in " << cpu_time << " seconds" << endl;
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

