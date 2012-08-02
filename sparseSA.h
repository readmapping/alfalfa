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
#include <iostream>
#include <algorithm>
#include <limits>
#include <sstream>
#include <assert.h>

#include "dp.h"
#include "utils.h"
#include <bitset>

using namespace std;

// Stores the LCP array in an unsigned char (0-255).  Values larger
// than or equal to 255 are stored in a sorted array.
// Simulates a vector<int> LCP;
struct vec_uchar {
  struct item_t{
    item_t(size_t i, int v) { idx = i; val = v; }
    size_t idx; int val;
    bool operator < (item_t t) const { return idx < t.idx;  }
  };
  vector<unsigned char> vec;  // LCP values from 0-65534
  vector<item_t> M;
  void resize(size_t N) { vec.resize(N); }
  // Vector X[i] notation to get LCP values.
  int operator[] (size_t idx) {
    if(vec[idx] == numeric_limits<unsigned char>::max())
      return lower_bound(M.begin(), M.end(), item_t(idx,0))->val;
    else
      return vec[idx];
  }
  // Actually set LCP values, distingushes large and small LCP
  // values.
  void set(size_t idx, int v) {
    if(v >= numeric_limits<unsigned char>::max()) {
      vec.at(idx) = numeric_limits<unsigned char>::max();
      M.push_back(item_t(idx, v));
    }
    else { vec.at(idx) = (unsigned char)v; }
  }
  // Once all the values are set, call init. This will assure the
  // values >= 255 are sorted by index for fast retrieval.
  void init() { sort(M.begin(), M.end()); cerr << "M.size()=" << M.size() << endl; }
};

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

// Match find by findMEM.
struct match_t {
  match_t(): ref(0), query(0), len(0) {}
  match_t(long r, long q, long l): ref(r), query(q), len(l) {}
  long ref; // position in reference sequence
  long query; // position in query
  long len; // length of match
};

// interval in match results + bases covering the result
struct lis_t {
  lis_t(): begin(0), end(0), len(0) {}
  lis_t(int b, int e, int l): begin(b), end(e), len(l) {}
  int begin; // position in reference sequence
  int end; // position in query
  int len; // length of match
};

struct alignment_t {
  alignment_t(): pos(0), cigar("*"), flag(0), rname("*"), mapq(0), tLength(0),
  rnext("*"), pnext(0), editDist(0), alignmentScore(0), cigarChars(0), cigarLengths(0), NMtag("*") {}
  string cigar;//TODO: remove this fields, only used when printed
  string NMtag;//TODO: remove this fields, only used when printed
  vector<char> cigarChars;//Change these to fixed-length values
  vector<int> cigarLengths;//Change these to fixed-length values (make sure to increase them when necessary): make own string and vector classes
  string rname;//leave out, only for printing
  long pos; // position in reference sequence
  bitset<11> flag;
  int mapq;
  //TODO: paired-end mapping
  string rnext;
  long pnext;
  int editDist;
  int alignmentScore;
  int tLength;
  void createCigar(bool newVersion){
      stringstream ss;
      assert(cigarChars.size()==cigarLengths.size());
      if(newVersion)
          for(int i = 0; i < cigarChars.size(); i++)
              ss << cigarLengths[i] << cigarChars[i];
      else{
          int i = 0;
          while(i < cigarChars.size()){
              if(cigarChars[i]=='=' || cigarChars[i]=='X'){
                    int tempLength = cigarLengths[i];
                    while(i < cigarChars.size()-1 && (cigarChars[i+1]=='=' || cigarChars[i+1]=='X')){
                        tempLength += cigarLengths[i+1];
                        i++;
                    }
                    ss << tempLength << 'M';
              }
              else{
                ss << cigarLengths[i] << cigarChars[i];
              }
              i++;
          }
      }
      cigar = ss.str();
  }
  void setFieldsFromCigar(const dp_scores & scores){//TODO: change to use special vectors
      stringstream sNM;
      stringstream sCig;
      assert(cigarChars.size()==cigarLengths.size());
      assert(cigarChars.size()>0);
      int i = 0;
      editDist = 0;
      while(i < cigarChars.size()){
          if(cigarChars[i]=='=' || cigarChars[i]=='X'){
                int tempLength = cigarLengths[i];
                sNM << cigarLengths[i] << cigarChars[i];
                editDist += (cigarChars[i]=='X' ? cigarLengths[i] : 0);
                alignmentScore += (cigarChars[i]=='X' ? scores.mismatch*cigarLengths[i] : scores.match*cigarLengths[i]);
                while(i < cigarChars.size()-1 && (cigarChars[i+1]=='=' || cigarChars[i+1]=='X')){
                    i++;
                    tempLength += cigarLengths[i];
                    sNM << cigarLengths[i] << cigarChars[i];
                    editDist += (cigarChars[i]=='X' ? cigarLengths[i] : 0);
                    alignmentScore += (cigarChars[i]=='X' ? scores.mismatch*cigarLengths[i] : scores.match*cigarLengths[i]);
                }
                sCig << tempLength << 'M';
          }
          else{
            sCig << cigarLengths[i] << cigarChars[i];
            sNM << cigarLengths[i] << cigarChars[i];
            if(cigarChars[i]=='D' || cigarChars[i]=='I'){
                editDist += cigarLengths[i];
                alignmentScore += scores.openGap + scores.extendGap*cigarLengths[i];
            }
          }
          i++;
      }
      cigar = sCig.str();
      NMtag = sNM.str();
  }
};

struct read_t {
    read_t(): qname(""),sequence(""),qual("*"), alignments(0) {}
    read_t(string name, string &seq, string &qw): qname(name),sequence(seq),qual(qw), alignments(0) {}
    void postprocess(const dp_scores & scores){
        int maxScore = scores.mismatch*sequence.length();
        assert(maxScore < 0);//works only for score<0 for now (to do: add sign switch to allow positive min scores)
        int secBestScore = scores.mismatch*sequence.length();
        for(int j = 0; j < alignments.size(); j++){
             alignments[j].setFieldsFromCigar(scores);
             if(alignments[j].alignmentScore > maxScore){
                 secBestScore = maxScore;
                 maxScore = alignments[j].alignmentScore;
             }
             else if(alignments[j].alignmentScore < maxScore && alignments[j].alignmentScore > secBestScore)
                 secBestScore = alignments[j].alignmentScore;
        }
        int mapq = (secBestScore == scores.mismatch*sequence.length()) ? 255 :
            250*(maxScore-secBestScore)/(maxScore-scores.mismatch*sequence.length()) ;
        assert(mapq <= 255 && mapq >= 0);
        for(int j = 0; j < alignments.size(); j++){
             if(alignments[j].alignmentScore == maxScore)
                 alignments[j].mapq = mapq;
             else
                 alignments[j].flag.set(8,true);
        }
    }
    void printAlignments(){
        if(alignments.empty()){
            printf("%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\n",qname.c_str(),4,"*",0,0,"*","*",0,0,sequence.c_str(),qual.c_str());
        }
        else{
            alignment_t & a = alignments[0];
            string revCompl = sequence;
            Utils::reverse_complement(revCompl, false);
            string qualRC = qual;
            reverse(qualRC.begin(),qualRC.end());
            printf("%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\tAS:i:%d\tNH:i:%ld\tNM:i:%d\tX0:Z:%s\n",
                        qname.c_str(),(int)a.flag.to_ulong(),a.rname.c_str(),a.pos,a.mapq,a.cigar.c_str(),
                        a.rnext.c_str(),a.pnext,a.tLength, a.flag.test(4) ? revCompl.c_str() : sequence.c_str(), a.flag.test(4) ? qualRC.c_str() : qual.c_str(),
                        a.alignmentScore,alignments.size(),a.editDist,a.NMtag.c_str());
            for(int i=1; i < alignments.size(); i++){
                a = alignments[i];
                printf("%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\tAS:i:%d\tNM:i:%d\tX0:Z:%s\n",
                        qname.c_str(),(int)a.flag.to_ulong(),a.rname.c_str(),a.pos,a.mapq,a.cigar.c_str(),
                        a.rnext.c_str(),a.pnext,a.tLength,a.flag.test(4) ? revCompl.c_str() : sequence.c_str(),
                        a.flag.test(4) ? qualRC.c_str() : qual.c_str(),a.alignmentScore,a.editDist,a.NMtag.c_str());
            }
        }
    }
    string qname;//TODO should be reference
    vector<alignment_t> alignments;
    string sequence;//TODO should be reference
    string qual;//TODO should be reference
};

// depth : [start...end]
struct interval_t {
  interval_t(): start(1), end(0), depth(-1) {}
  interval_t(long s, long e, long d): start(s), end(e), depth(d) {}
  void reset(long e) { start = 0; end = e; depth = 0; }
  long depth, start, end;
  long size() { return end - start + 1; }
};

struct sparseSA {
  vector<string> &descr; // Descriptions of concatenated sequences.
  vector<long> &startpos; // Lengths of concatenated sequences.
  long maxdescrlen; // Maximum length of the sequence description, used for formatting.
  bool _4column; // Use 4 column output format.

  long N; //!< Length of the sequence.
  long logN; // ceil(log(N))
  long NKm1; // N/K - 1
  string &S; //!< Reference to sequence data.
  vector<unsigned int> SA;  // Suffix array.
  vector<int> ISA;  // Inverse suffix array.
  vec_uchar LCP; // Simulates a vector<int> LCP.
  vector<int> CHILD; //child table

  long K; // suffix sampling, K = 1 every suffix, K = 2 every other suffix, K = 3, every 3rd sffix
  bool hasChild;
  bool hasSufLink;

  // Maps a hit in the concatenated sequence set to a position in that sequence.
  void from_set(long hit, long &seq, long &seqpos) {
    // Use binary search to locate index of sequence and position
    // within sequence.
    vector<long>::iterator it = upper_bound(startpos.begin(), startpos.end(), hit);
    seq = distance(startpos.begin(), it) - 1;
    it--;
    seqpos = hit - *it;
  }

  // Constructor builds sparse suffix array.
  sparseSA(string &S_, vector<string> &descr_, vector<long> &startpos_, bool __4column, long K_);

  // Modified Kasai et all for LCP computation.
  void computeLCP();
  //Modified Abouelhoda et all for CHILD Computation.
  void computeChild();

  // Radix sort required to construct transformed text for sparse SA construction.
  void radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h);

  // Prints match to cout.
  void print_match(const match_t m);
  void print_match(const match_t m, vector<match_t> &buf); // buffered version
  void print_match(const string meta, vector<match_t> &buf, bool rc); // buffered version

  // Binary search for left boundry of interval.
  inline long bsearch_left(char c, long i, long s, long e);
  // Binary search for right boundry of interval.
  inline long bsearch_right(char c, long i, long s, long e);

  // Simple suffix array search.
  inline bool search(const string &P, long &start, long &end);

  // Simple top down traversal of a suffix array.
  inline bool top_down(char c, long i, long &start, long &end);
  inline bool top_down_faster(char c, long i, long &start, long &end);
  inline bool top_down_child(char c, interval_t &cur);

  // Traverse pattern P starting from a given prefix and interval
  // until mismatch or min_len characters reached.
  inline void traverse(const string &P, long prefix, interval_t &cur, int min_len);
  inline void traverse_faster(const string &P,const long prefix, interval_t &cur, int min_len);

  // Simulate a suffix link.
  inline bool suffixlink(interval_t &m);

  // Expand ISA/LCP interval. Used to simulate suffix links.
  inline bool expand_link(interval_t &link) {
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

  // Given a position i in S, finds a left maximal match of minimum
  // length within K steps.
  inline void find_Lmaximal(const string &P, long prefix, long i, long len, vector<match_t> &matches, int min_len, bool print);

  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.
  void collectMEMs(const string &P, long prefix, const interval_t mli, 
  interval_t xmi, vector<match_t> &matches, int min_len, bool print);

  void collectSMAMs(const string &P, long prefix, const interval_t mli, 
  interval_t xmi, vector<match_t> &matches, int min_len, int maxCount, bool print);

  // Find all MEMs given a prefix pattern offset k.
  void findMEM(long k, const string &P, vector<match_t> &matches, int min_len, bool print);

  // Find all MEMs given a prefix pattern offset k.
  void findSMAM(long k, const string &P, vector<match_t> &matches, int min_len, int maxCount, bool print);

  // NOTE: min_len must be > 1
  void findMAM(const string &P, vector<match_t> &matches, int min_len, bool print);
  inline bool is_leftmaximal(const string &P, long p1, long p2);

  // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
  // sequence S. as computed by MUMmer version 2 by Salzberg
  // et. al. Note this is a "one-sided" query. It "streams" the query
  // P throught he index.  Consequently, repeats can occur in the
  // pattern P.
  void MAM(const string &P, vector<match_t> &matches, int min_len, bool print) {
    if(K != 1) return;  // Only valid for full suffix array.
    findMAM(P, matches, min_len, print);
  }

  // Find Maximal Exact Matches (MEMs)
  void MEM(const string &P, vector<match_t> &matches, int min_len, bool print, int num_threads = 1);

  // Find Maximal Exact Matches (MEMs)
  void SMAM(const string &P, vector<match_t> &matches, int min_len, int maxCount, bool print, int num_threads = 1);

  // Maximal Unique Match (MUM)
  void MUM(const string &P, vector<match_t> &unique, int min_len, bool print);

  //post process MEMs, MAMs, etc...
  void postProcess(vector<match_t> &matches);

  void inexactMatch(read_t& read,const align_opt & alnOptions, bool fwStrand, bool print);
  //TODO: calculate global position in above function
};


#endif // __sparseSA_hpp__

