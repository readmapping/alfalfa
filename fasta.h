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
 * Part of the code originates from "A Practical Algorithm for Finding 
 * Maximal Exact Matches in Large Sequence Data Sets Using Sparse Suffix Arrays"
 * By Khan et al. Copyright (c) 2009 
 * You should have received a copy of the copyright notice with this code.
 */

#ifndef FASTA_H
#define	FASTA_H

#include <string>
#include <vector>
#include <pthread.h>
#include <iostream>
#include <fstream>

using namespace std;

void reverse_complement(string &seq_rc, bool nucleotides_only);
void trim(string &line, long &start, long &end);
void load_fasta(string filename, string &S, vector<string> &descr, vector<long> &startpos);

enum filetype_t { UNKNOWN, FASTA, FASTQ };

/**
 * Global reader class for the sequence reads
 * @param fileName
 * @param fileT
 */
struct fastqInputReader {
    fastqInputReader();
    fastqInputReader(const string& fileName, bool nucOnly, filetype_t fileT = UNKNOWN);
    ~fastqInputReader();
    void open(const string& fileName, bool nucOnly, filetype_t fileT = UNKNOWN);
    
    void determineType();
    //returns 0 if succeeded
    bool nextRead(string& meta, string& sequence, string& qualities);
    bool nextReadFastQ(string& meta, string& sequence, string& qualities);
    bool nextReadFastA(string& meta, string& sequence);
    
    bool nucleotidesOnly;
    ifstream data;
    filetype_t fileType;
    string filename;
    pthread_mutex_t readLock_;
};

#endif	/* FASTA_H */

