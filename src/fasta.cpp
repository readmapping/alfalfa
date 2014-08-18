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

#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>

#include "fasta.h"

using namespace std;

// Trim a string, giving start and end in trimmed version.
// NOTE: Assumes line.length() > 0!!!!
void trim(string &line, long &start, long &end) {
  // Trim leading spaces.
  for(long i = start; i < (int)line.length(); i++) {
    if(line[i] != ' ') { start = i; break; }
  }
  // Trim trailing spaces.
  for(long i = line.length() - 1; i >= 0; i--) {
    if(line[i] != ' ') { end = i; break; }
    if(i == 0) break;
  }
}


// Concatenate new sequences to set, keep track of lengths.
// NOTE: Concatenation using the '`' character to separate strings!
void load_fasta(string filename, string &S, vector<string> &descr, vector<long> &startpos) {
  string meta, line;
  long length = 0;

  // Everything starts at zero.
  startpos.push_back(0);

  ifstream data(filename.c_str());

  if(!data.is_open()) { cerr << "unable to open " << filename << endl; exit(1); }

  while(!data.eof()) {
    getline(data, line, '\n'); // Load one line at a time.
    if(line.length() == 0) continue;
    if(line[line.length()-1]=='\r')
        line.erase(--line.end());

    long start = 0, end = line.length() - 1;

    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
      // Save previous sequence and meta data.
      if(length > 0) {
	descr.push_back(meta);
	S += '`'; // ` character used to separate strings
	startpos.push_back(S.length());
	//lengths.push_back(length+1);
      }
      // Reset parser state.
      start = 1; meta = ""; length = 0;
    }
    trim(line, start, end);
    // Collect meta data.
    if(line[0] == '>') {
      for(long i = start; i <= end; i++) { if(line[i] == ' ') break; meta += line[i]; }
    }
    else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) {
	S += std::tolower(line[i]);
      }
    }
  }
  if(length > 0) {
    descr.push_back(meta);
  }
  cerr << "# full reference length = " << S.length() << endl;
  for(long i = 0; i < (long)descr.size(); i++) {
    cerr << "# " << descr[i] << " " << startpos[i] << endl;
  }
}

fastqInputReader::fastqInputReader(): nucleotidesOnly(false), filename("") {
    pthread_mutex_init(&readLock_, NULL);
}

fastqInputReader::fastqInputReader(const string& fName, bool nucOnly): 
        nucleotidesOnly(nucOnly), filename(fName) {
    pthread_mutex_init(&readLock_, NULL);
    fp = gzopen(fName.c_str(),"r");
    if(fp == NULL){
        fprintf(stderr, "ERROR: file %s could not be read.", fName.c_str());
        exit(1);
    }
    seq = kseq_init(fp);
}

void fastqInputReader::open(const string& fileName, bool nucOnly){
    nucleotidesOnly = nucOnly;
    filename = fileName;
    fp = gzopen(fileName.c_str(),"r");
    if(fp == NULL){
        fprintf(stderr, "ERROR: file %s could not be read.", fileName.c_str());
        exit(1);
    }
    seq = kseq_init(fp);
}

fastqInputReader::~fastqInputReader(){
    kseq_destroy(seq);
    gzclose(fp);
    pthread_mutex_destroy(&readLock_);
}

bool fastqInputReader::nextRead(string& meta, string& sequence, string& qualities){
    l = kseq_read(seq);
    if(l < 0)
        return false;
    else{
        meta = seq->name.s;
        sequence = seq->seq.s;
        for(size_t i = 0; i <= sequence.length(); i++) {
            char c = std::tolower(sequence[i]);
            if(nucleotidesOnly) {
                switch(c) {
                    case 'a': case 't': case 'g': case 'c': break;
                    default:
                        c = '~';
                }
            }
            sequence[i] = c;
        }
        
        if (seq->qual.l)
            qualities = seq->qual.s;
        
        return true;
    }
}
