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

#ifndef sparseSA_h__
#define sparseSA_h__

#include <vector>
#include <string>
#include <limits>
#include <limits.h>
#include <algorithm>//sort + lower_bound
#include <stdio.h>//cerr
#include <iostream>

static const long KMERSIZE = 10;
static const long TABLESIZE = 1 << (2*KMERSIZE);//10 is k-mer size
static const unsigned int BITADD[256] = { UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//0-9
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//10-19
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//20-29
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//30-39
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//40-49
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//50-59
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 0,        UINT_MAX, 1,        UINT_MAX, UINT_MAX,//60-69 65:A, 67:C
                                         UINT_MAX, 2,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//70-79 71:G
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 3,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//80-89 84:T
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 0,        UINT_MAX, 1,       //90-99 97:a, 99: c
                                         UINT_MAX, UINT_MAX, UINT_MAX, 2,        UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//100-109 103:g
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, 3,        UINT_MAX, UINT_MAX, UINT_MAX,//110-119 116:t
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//120-129
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//130-139
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//140-149
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//150-159
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//160-169
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//170-179
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//180-189
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//190-199
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//200-209
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//210-219
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//220-229
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//230-239
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX,//240-249
                                         UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX };//250-255

// Stores the LCP array in an unsigned char (0-255).  Values larger
// than or equal to 255 are stored in a sorted array.
// Simulates a vector<int> LCP
struct vec_uchar {
  struct item_t{
    item_t(){}
    item_t(size_t i, int v) { idx = i; val = v; }
    size_t idx; int val;
    bool operator < (item_t t) const { return idx < t.idx;  }
  };
  std::vector<unsigned char> vec;  // LCP values from 0-65534
  std::vector<item_t> M;
  void resize(size_t N){ vec.resize(N); }
  // Vector X[i] notation to get LCP values.
  int operator[] (size_t idx) const {
    if(vec[idx] == std::numeric_limits<unsigned char>::max())
        return std::lower_bound(M.begin(), M.end(), item_t(idx,0))->val;
    else
        return vec[idx];
  }
  // Actually set LCP values, distingushes large and small LCP
  // values.
  void set(size_t idx, int v){
    if(v >= std::numeric_limits<unsigned char>::max()) {
        vec.at(idx) = std::numeric_limits<unsigned char>::max();
        M.push_back(item_t(idx, v));
    }
    else { vec.at(idx) = (unsigned char)v; }
  }
  
  // Once all the values are set, call init. This will assure the
  // values >= 255 are sorted by index for fast retrieval.
  void init();
  
  long index_size_in_bytes();
};

// Match find by findMEM.
struct match_t {
  match_t(): ref(0), query(0), len(0) {}
  match_t(long r, int q, int l): ref(r), query(q), len(l) {}
  match_t(const match_t & o): ref(o.ref), query(o.query), len(o.len) {}
  long ref; // position in reference sequence
  int query; // position in query
  int len; // length of match
};

struct saTuple_t {
    saTuple_t(): left(0), right(0) {}
    saTuple_t(unsigned int l, unsigned int r): left(l), right(r){}
    unsigned int left;
    unsigned int right;
};

// depth : [start...end]
struct interval_t {
  interval_t(): start(1), end(0), depth(-1) {}
  interval_t(long s, long e, int d): start(s), end(e), depth(d) {}
  interval_t(const interval_t & other): start(other.start), end(other.end), depth(other.depth) {}
  void reset(long e) { start = 0; end = e; depth = 0; }
  long start, end; 
  int depth;
  inline long size() const { return end - start + 1; }
};

// depth : [start...end]
struct intervTuple_t {
  intervTuple_t(): mli(), xmi(), prefix(0) {}
  intervTuple_t(const interval_t &m, const interval_t &x, int pref): mli(m), xmi(x), prefix(pref) {}
  intervTuple_t(const intervTuple_t & other): mli(other.mli), xmi(other.xmi), prefix(other.prefix) {}
  interval_t mli, xmi;
  int prefix;
  int value() const { return xmi.depth; }
  int xmiEnd() const { return prefix + xmi.depth; }
  bool max() const { return mli.depth < 0; }
};

struct sparseSA {
  std::vector<std::string> descr; // Descriptions of concatenated sequences.
  std::vector<long> startpos; // Lengths of concatenated sequences.
  long maxdescrlen; // Maximum length of the sequence description, used for formatting.

  long N; //!< Length of the sequence.
  long logN; // ceil(log(N))
  long NKm1; // N/K - 1
  std::string S; //!< Reference to sequence data.
  std::vector<unsigned int> SA;  // Suffix array.
  std::vector<int> ISA;  // Inverse suffix array.
  vec_uchar LCP; // Simulates a vector<int> LCP.
  std::vector<int> CHILD; //child table

  long K; // suffix sampling, K = 1 every suffix, K = 2 every other suffix, K = 3, every 3rd sffix
  bool hasChild;
  bool hasSufLink;
  int sparseMult;
  
  //fields for lookup table of sa intervals to a certain small depth
  bool hasKmer;
  std::vector<saTuple_t> KMR;
  
  //CONSTRUCT
  // Constructor builds sparse suffix array.
  sparseSA(){
      sparseMult = 1;
  }
  //ACCESS INFO FUNCTIONS
  long index_size_in_bytes();
  // Maps a hit in the concatenated sequence set to a position in that sequence.
  void from_set(long hit, long &seq, long &seqpos) const{
    // Use binary search to locate index of sequence and position
    // within sequence.
    std::vector<long>::const_iterator it = std::upper_bound(startpos.begin(), startpos.end(), hit);
    seq = std::distance(startpos.begin(), it) - 1;
    it--;
    seqpos = hit - *it;
  }
  // Get the index of the start and the end of the chromosome in the concatenated sequence
  void getChromBounds(const long hit, long& chrStart, long& chrEnd) const{
    long seq, seqpos;
    from_set(hit, seq, seqpos);
    chrStart = startpos[seq];
    if(startpos.size() > 1 && seq < (long)(startpos.size() - 1))
        chrEnd = startpos[seq+1] - 1; // Minus one for the separation character
    else
        chrEnd =  N - 1; // Minus one for the separation character
  }
  
  //CREATION
  void loadRef(const std::string & fileName);
  void init(long K_, bool suflink_, bool child_, bool kmer_);
  // Modified Kasai et all for LCP computation.
  void computeLCP();
  //Modified Abouelhoda et all for CHILD Computation.
  void computeChild();
  //build look-up table for sa intervals of kmers up to depth 10
  void computeKmer();
  // Radix sort required to construct transformed text for sparse SA construction.
  void radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h);
  //save index to files
  void save(const std::string &prefix);
  //load index from file
  bool load(const std::string &prefix);
  //construct
  void construct();

  //SA-SEARCH METHODS
  inline bool top_down_faster(char c, long i, long &start, long &end) const;
  inline bool top_down_child(char c, interval_t &cur) const;

  //TRAVERSE FUNCTIONS
  // Traverse pattern P starting from a given prefix and interval
  // until mismatch or min_len characters reached.
  inline void traverse(const std::string &P, int prefix, interval_t &cur, int min_len, int maxCount) const;
  inline void traverse_faster(const std::string &P,const int prefix, interval_t &cur, int min_len, int maxCount) const;
  inline void traverse_keep(const std::string &P,const int prefix, std::vector<interval_t> &cur, int min_len, int maxCount) const;
  // Simulate a suffix link.
  inline bool suffixlink(interval_t &m) const;
  // Expand ISA/LCP interval. Used to simulate suffix links.
  inline bool expand_link(interval_t &link) const {
    long thresh = 2 * link.depth * logN, exp = 0; // Threshold link expansion.
    long start = link.start;
    long end = link.end;
    while(LCP[start] >= link.depth) {
      exp++;
      if(exp >= thresh) return false;
      start--;
    }
    while(end < NKm1 && LCP[end+1] >= link.depth) {
      exp++;
      if(exp >= thresh) return false;
      end++;
    }
    link.start = start; link.end = end;
    return true;
  }

  //MEM-finding functions
  
  // Given a position i in S, finds a left maximal match of minimum
  // length within K steps.  
  inline void find_LmaximalStrict(const std::string &P, int prefix, long i, int len, std::vector<match_t> &matches, int min_len, int border) const;

  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.    
  void collectAllMEMs(const std::string &P, int prefix, const interval_t mli, 
  interval_t xmi, std::vector<match_t> &matches, int min_len, int maxLeft) const;

  void collectAllSMAMs(const std::string &P, int prefix, interval_t xmi, 
  std::vector<match_t> &matches, int min_len, int maxCount, int maxLeft) const;
  
  void collectXmiAndMli(const std::string &P, intervTuple_t interval, 
  std::vector<match_t> &matches, bool calcMli, int nextEndStop, 
  int nextEndStopPrefix, int min_len, int maxCount, int maxLeft) const;
  
  void collectIntervals(const std::string &P, int prefix, std::vector<interval_t> & intervals, 
  std::vector<match_t> &matches, bool calcMli, int nextEndStop, int min_len, int maxCount, int maxLeft) const;

  // Find all MEMs given a prefix pattern offset k.
  void findMEM(int k, const std::string &P, std::vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const;

  // Find all MEMs given a prefix pattern offset k.
  void findSMAM(int k, const std::string &P, std::vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const;

  // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
  // sequence S. as computed by MUMmer version 2 by Salzberg
  // et. al. Note this is a "one-sided" query. It "streams" the query
  // P throught he index.  Consequently, repeats can occur in the
  // pattern P.
  // Find Maximal Exact Matches (MEMs)
  void MEM(const std::string &P, std::vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const;

  // Find Maximal Exact Matches (MEMs)
  void SMAM(const std::string &P, std::vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const;
  
  // Find Maximal Exact Matches (MEMs)
  void SMEM(const std::string &P, std::vector<match_t> &matches, int min_len, int maxCount, int sparseQ, int max_smem_cont) const;
  
  // Find Maximal Exact Matches (MEMs)
  void SMEMP(const std::string &P, std::vector<match_t> &matches, int min_len, int maxCount, int sparseQ, int max_smem_cont) const;

};


#endif // __sparseSA_hpp__

