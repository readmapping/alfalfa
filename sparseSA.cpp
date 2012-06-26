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

#include <stdio.h>
#include <math.h>
#include <pthread.h>


#include <bitset>
#include <limits.h>
#include <stack>

#include "sparseSA.h"
#include "dp.h"
#include "utils.h"

// LS suffix sorter (integer alphabet).
extern "C" { void suffixsort(int *x, int *p, int n, int k, int l); };

pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;

sparseSA::sparseSA(string &S_, vector<string> &descr_, vector<long> &startpos_, bool __4column, long K_) :
  descr(descr_), startpos(startpos_), S(S_) {
  _4column = __4column;

  // Get maximum query sequence description length.
  maxdescrlen = 0;
  for(long i = 0; i < (long)descr.size(); i++) {
    if(maxdescrlen < (long)descr[i].length())  maxdescrlen = descr[i].length();
  }
  K = K_;

  // Increase string length so divisible by K.
  // Don't forget to count $ termination character.
  if(S.length() % K != 0) {
    long appendK = K - S.length() % K ;
    for(long i = 0; i < appendK; i++) S += '$';
  }
  // Make sure last K-sampled characeter is this special character as well!!
  for(long i = 0; i < K; i++) S += '$'; // Append "special" end character. Note: It must be lexicographically less.
  N = S.length();

  if(K > 1) {
    long bucketNr = 1;
    int *intSA = new int[N/K+1];  for(int i = 0; i < N/K; i++) intSA[i] = i; // Init SA.
    int* t_new = new int[N/K+1];
    long* BucketBegin = new long[256]; // array to save current bucket beginnings
    radixStep(t_new, intSA, bucketNr, BucketBegin, 0, N/K-1, 0); // start radix sort
    t_new[N/K] = 0; // Terminate new integer string.
    delete[] BucketBegin;

    // Suffix sort integer text.
    cerr << "# suffixsort()" << endl;
    suffixsort(t_new, intSA, N/K, bucketNr, 0);
    cerr << "# DONE suffixsort()" << endl;

    delete[] t_new;

    // Translate suffix array.
    SA.resize(N/K);
    for (long i=0; i<N/K; i++) SA[i] = (unsigned int)intSA[i+1] * K;
    delete[] intSA;

    // Build ISA using sparse SA.
    ISA.resize(N/K);
    for(long i = 0; i < N/K; i++) { ISA[SA[i]/K] = i; }
  }
  else {
    SA.resize(N);
    ISA.resize(N);
    int char2int[UCHAR_MAX+1]; // Map from char to integer alphabet.

    // Zero char2int mapping.
    for (int i=0; i<=UCHAR_MAX; i++) char2int[i]=0;

    // Determine which characters are used in the string S.
    for (long i = 0; i < N; i++) char2int[(int)S[i]]=1;

    // Count the size of the alphabet.
    int alphasz = 0;
    for(int i=0; i <= UCHAR_MAX; i++) {
      if (char2int[i]) char2int[i]=alphasz++;
      else char2int[i] = -1;
    }

    // Remap the alphabet.
    for(long i = 0; i < N; i++) ISA[i] = (int)S[i];
    for (long i = 0; i < N; i++) ISA[i]=char2int[ISA[i]] + 1;
    // First "character" equals 1 because of above plus one, l=1 in suffixsort().
    int alphalast = alphasz + 1;

    // Use LS algorithm to construct the suffix array.
    int *SAint = (int*)(&SA[0]);
    suffixsort(&ISA[0], SAint , N-1, alphalast, 1);
  }

  cerr << "N=" << N << endl;

  // Adjust to "sampled" size.
  logN = (long)ceil(log(N/K) / log(2.0));

  LCP.resize(N/K);
  cerr << "N/K=" << N/K << endl;
  // Use algorithm by Kasai et al to construct LCP array.
  computeLCP();  // SA + ISA -> LCP
    LCP.init();
    if(K >= 4){
        hasChild = true;
        hasSufLink = false;
        ISA.clear();
        CHILD.resize(N/K);
        //Use algorithm by Abouelhoda et al to construct CHILD array
        computeChild();
    }
    else{
        hasChild = false;
        hasSufLink = true;
    }

    NKm1 = N/K-1;
}

// Uses the algorithm of Kasai et al 2001 which was described in
// Manzini 2004 to compute the LCP array. Modified to handle sparse
// suffix arrays and inverse sparse suffix arrays.
void sparseSA::computeLCP() {
  long h=0;
  for(long i = 0; i < N; i+=K) {
    long m = ISA[i/K];
    if(m==0) LCP.set(m, 0); // LCP[m]=0;
    else {
      long j = SA[m-1];
      while(i+h < N && j+h < N && S[i+h] == S[j+h])  h++;
      LCP.set(m, h); //LCP[m] = h;
    }
    h = max(0L, h - K);
  }
}

// Child array construction algorithm
void sparseSA::computeChild() {
    for(int i = 0; i < N/K; i++){
        CHILD[i] = -1;
    }
        //Compute up and down values
        int lastIndex  = -1;
        stack<int,vector<int> > stapelUD;
        stapelUD.push(0);
        for(int i = 1; i < N/K; i++){
            while(LCP[i] < LCP[stapelUD.top()]){
                lastIndex = stapelUD.top();
                stapelUD.pop();
                if(LCP[i] <= LCP[stapelUD.top()] && LCP[stapelUD.top()] != LCP[lastIndex]){
                    CHILD[stapelUD.top()] = lastIndex;
                }
            }
            //now LCP[i] >= LCP[top] holds
            if(lastIndex != -1){
                CHILD[i-1] = lastIndex;
                lastIndex = -1;
            }
            stapelUD.push(i);
        }
        while(0 < LCP[stapelUD.top()]){//last row (fix for last character of sequence not being unique
            lastIndex = stapelUD.top();
            stapelUD.pop();
            if(0 <= LCP[stapelUD.top()] && LCP[stapelUD.top()] != LCP[lastIndex]){
                CHILD[stapelUD.top()] = lastIndex;
            }
        }
        //Compute Next L-index values
        stack<int,vector<int> > stapelNL;
        stapelNL.push(0);
        for(int i = 1; i < N/K; i++){
            while(LCP[i] < LCP[stapelNL.top()])
                stapelNL.pop();
            lastIndex = stapelNL.top();
            if(LCP[i] == LCP[lastIndex]){
                stapelNL.pop();
                CHILD[lastIndex] = i;
            }
            stapelNL.push(i);
        }
}


// Implements a variant of American flag sort (McIlroy radix sort).
// Recurse until big-K size prefixes are sorted. Adapted from the C++
// source code for the wordSA implementation from the following paper:
// Ferragina and Fischer. Suffix Arrays on Words. CPM 2007.
void sparseSA::radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h) {
  if(h >= K) return;
  // first pass: count
  vector<long> Sigma(256, 0); // Sigma counts occurring characters in bucket
  for (long i = l; i <= r; i++) Sigma[ S[ SA[i]*K + h ] ]++; // count characters
  BucketBegin[0] = l; for (long i = 1; i < 256; i++) { BucketBegin[i] = Sigma[i-1] + BucketBegin[i-1]; } // accumulate

  // second pass: move (this variant does *not* need an additional array!)
  unsigned char currentKey = 0;    // character of current bucket
  long end = l-1+Sigma[currentKey]; // end of current bucket
  long pos = l;                     // 'pos' is current position in bucket
  while (1) {
    if (pos > end) { // Reached the end of the bucket.
      if (currentKey == 255) break; // Last character?
      currentKey++; // Advance to next characer.
      pos = BucketBegin[currentKey]; // Next bucket start.
      end += Sigma[currentKey]; // Next bucket end.
    }
    else {
      // American flag sort of McIlroy et al. 1993. BucketBegin keeps
      // track of current position where to add to bucket set.
      int tmp = SA[ BucketBegin[ S[ SA[pos]*K + h ] ] ];
      SA[ BucketBegin[ S[ SA[pos]*K + h] ]++ ] = SA[pos];  // Move bucket beginning to the right, and replace
      SA[ pos ] = tmp; // Save value at bucket beginning.
      if (S[ SA[pos]*K + h ] == currentKey) pos++; // Advance to next position if the right character.
    }
  }
  // recursively refine buckets and calculate new text:
  long beg = l; end = l-1;
  for (long i = 1; i < 256; i++) { // step through Sigma to find bucket borders
    end += Sigma[i];
    if (beg <= end) {
      if(h == K-1) {
	for (long j = beg; j <= end; j++) {
	  t_new[ SA[j] ] = bucketNr; // set new text
	}
	bucketNr++;
      }
      else {
	radixStep(t_new, SA, bucketNr, BucketBegin, beg, end, h+1); // recursive refinement
      }
      beg = end + 1; // advance to next bucket
    }
  }
}

// Binary search for left boundry of interval.
long sparseSA::bsearch_left(char c, long i, long s, long e) {
  if(c == S[SA[s]+i]) return s;
  long l = s, r = e;
  while (r - l > 1) {
    long m = (l+r) / 2;
    if (c <= S[SA[m] + i]) r = m;
    else l = m;
  }
  return r;
}

// Binary search for right boundry of interval.
long sparseSA::bsearch_right(char c, long i, long s, long e) {
  if(c == S[SA[e]+i]) return e;
  long l = s, r = e;
  while (r - l > 1) {
    long m = (l+r) / 2;
    if (c < S[SA[m] + i]) r = m;
    else l = m;
  }
  return l;
}


// Simple top down traversal of a suffix array.
bool sparseSA::top_down(char c, long i, long &start, long &end) {
  if(c < S[SA[start]+i]) return false;
  if(c > S[SA[end]+i]) return false;
  long l = bsearch_left(c, i, start, end);
  long l2 = bsearch_right(c, i, start, end);
  start = l; end = l2;
  return l <= l2;
}

// Top down traversal of the suffix array to match a pattern.  NOTE:
// NO childtab as in the enhanced suffix array (ESA).
bool sparseSA::search(const string &P, long &start, long &end) {
  start = 0; end = N - 1;
  long i = 0;
  while(i < (long)P.length()) {
    if(top_down(P[i], i, start, end) == false) {
      return false;
    }
    i++;
  }
  return true;
}


// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
void sparseSA::traverse(const string &P, long prefix, interval_t &cur, int min_len) {
  if(cur.depth >= min_len) return;

  while(prefix+cur.depth < (long)P.length()) {
    long start = cur.start; long end = cur.end;
    // If we reach a mismatch, stop.
    if(top_down_faster(P[prefix+cur.depth], cur.depth, start, end) == false) return;

    // Advance to next interval.
    cur.depth += 1; cur.start = start; cur.end = end;

    // If we reach min_len, stop.
    if(cur.depth == min_len) return;
  }
}

// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
// Uses the child table for faster traversal
void sparseSA::traverse_faster(const string &P,const long prefix, interval_t &cur, int min_len){
        if(cur.depth >= min_len) return;
        int c = prefix + cur.depth;
        bool intervalFound = c < P.length();
        if(intervalFound && cur.size() > 1)
            intervalFound = top_down_child(P[c], cur);
        else if(intervalFound)
            intervalFound = P[c] == S[SA[cur.start]+cur.depth];
        bool mismatchFound = false;
        while(intervalFound && !mismatchFound &&
                c < P.length() && cur.depth < min_len){
            c++;
            cur.depth++;
            if(cur.start != cur.end){
                int childLCP;
                //calculate LCP of child node, which is now cur. the LCP value
                //of the parent is currently c - prefix
                if(cur.start < CHILD[cur.end] && CHILD[cur.end] <= cur.end)
                    childLCP = LCP[CHILD[cur.end]];
                else
                    childLCP = LCP[CHILD[cur.start]];
                int minimum = min(childLCP,min_len);
                //match along branch
                while(!mismatchFound && c < P.length() && cur.depth < minimum){
                    mismatchFound = S[SA[cur.start]+cur.depth] != P[c];
                    c++;
                    cur.depth += !mismatchFound;
                }
                intervalFound = c < P.length() && !mismatchFound &&
                        cur.depth < min_len && top_down_child(P[c], cur);
            }
            else{
                while(!mismatchFound && c < P.length() && cur.depth < min_len){
                    mismatchFound = SA[cur.start]+cur.depth >= S.length() ||
                            S[SA[cur.start]+cur.depth] != P[c];
                    c++;
                    cur.depth += !mismatchFound;
                }
            }
        }
}

//finds the child interval of cur that starts with character c
//updates left and right bounds of cur to child interval if found, or returns
//cur if not found (also returns true/false if found or not)
bool sparseSA::top_down_child(char c, interval_t &cur){
    long left = cur.start;
    long right = CHILD[cur.end];
    if(cur.start >= right || right > cur.end)
        right = CHILD[cur.start];
    //now left and right point to first child
    if(S[SA[cur.start]+cur.depth] == c){
        cur.end = right-1;
        return true;
    }
    left = right;
    //while has next L-index
    while(CHILD[right] > right && LCP[right] == LCP[CHILD[right]]){
        right = CHILD[right];
        if(S[SA[left]+cur.depth] == c){
            cur.start = left; cur.end = right - 1;
            return true;
        }
        left = right;
    }
    //last interval
    if(S[SA[left]+cur.depth] == c){
            cur.start = left;
            return true;
    }
    return false;
}

// Given SA interval apply binary search to match character c at
// position i in the search string. Adapted from the C++ source code
// for the wordSA implementation from the following paper: Ferragina
// and Fischer. Suffix Arrays on Words. CPM 2007.
bool sparseSA::top_down_faster(char c, long i, long &start, long &end) {
  long l, r, m, r2=end, l2=start, vgl;
  bool found = false;
  long cmp_with_first = (long)c - (long)S[SA[start]+i];
  long cmp_with_last = (long)c - (long)S[SA[end]+i];
  if(cmp_with_first < 0) {
    l = start+1; l2 = start; // pattern doesn't occur!
  }
  else if(cmp_with_last > 0) {
    l = end+1; l2 = end;
    // pattern doesn't occur!
  }
  else {
    // search for left border:
    l = start; r = end;
    if (cmp_with_first == 0) {
      found = true; r2 = r;
    }
    else {
      while (r - l > 1) {
	m = (l+r) / 2;
	vgl = (long)c - (long)S[SA[m] + i];
	if (vgl <= 0) {
	  if (!found && vgl == 0) {
	    found = true;
	    l2 = m; r2 = r; // search interval for right border
	  }
	  r = m;
	}
	else l = m;
      }
      l = r;
    }
    // search for right border (in the range [l2:r2])
    if (!found) {
      l2 = l - 1; // pattern not found => right border to the left of 'l'
    }
    if (cmp_with_last == 0) {
      l2 = end; // right border is the end of the array
    }
    else {
      while (r2 - l2 > 1) {
	m = (l2 + r2) / 2;
	vgl = (long)c - (long)S[SA[m] + i];
	if (vgl < 0) r2 = m;
	else l2 = m;
      }
    }
  }
  start = l;
  end = l2;
  return l <= l2;
}


// Suffix link simulation using ISA/LCP heuristic.
bool sparseSA::suffixlink(interval_t &m) {
  m.depth -= K;
  if( m.depth <= 0) return false;
  m.start = ISA[SA[m.start] / K + 1];
  m.end = ISA[SA[m.end] / K + 1];
  return expand_link(m);
}

// For a given offset in the prefix k, find all MEMs.
void sparseSA::findMEM(long k, const string &P, vector<match_t> &matches, int min_len, bool print) {
  if(k < 0 || k >= K) { cerr << "Invalid k." << endl; return; }
  // Offset all intervals at different start points.
  long prefix = k;
  interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval

  // Right-most match used to terminate search.
  int min_lenK = min_len - (K-1);

  while( prefix <= (long)P.length() - (K-k)) {
    traverse(P, prefix, mli, min_lenK);    // Traverse until minimum length matched.
    if(mli.depth > xmi.depth) xmi = mli;
    if(mli.depth <= 1) { mli.reset(N/K-1); xmi.reset(N/K-1); prefix+=K; continue; }

    if(mli.depth >= min_lenK) {
      traverse(P, prefix, xmi, P.length()); // Traverse until mismatch.
      collectMEMs(P, prefix, mli, xmi, matches, min_len, print); // Using LCP info to find MEM length.
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix+=K;
      if( suffixlink(mli) == false ) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      suffixlink(xmi);
    }
    else {
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix+=K;
      if( suffixlink(mli) == false ) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      xmi = mli;
    }
  }
#ifndef NDEBUG
  if(print) print_match(match_t(), matches);   // Clear buffered matches.
#endif
}


// Use LCP information to locate right maximal matches. Test each for
// left maximality.
void sparseSA::collectMEMs(const string &P, long prefix, const interval_t mli,
        interval_t xmi, vector<match_t> &matches, int min_len, bool print) {
  // All of the suffixes in xmi's interval are right maximal.
  for(long i = xmi.start; i <= xmi.end; i++) find_Lmaximal(P, prefix, SA[i], xmi.depth, matches, min_len, print);

  if(mli.start == xmi.start && mli.end == xmi.end) return;

  while(xmi.depth >= mli.depth) {
    // Attempt to "unmatch" xmi using LCP information.
    if(xmi.end+1 < N/K) xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
    else xmi.depth = LCP[xmi.start];

    // If unmatched XMI is > matched depth from mli, then examine rmems.
    if(xmi.depth >= mli.depth) {
      // Scan RMEMs to the left, check their left maximality..
      while(LCP[xmi.start] >= xmi.depth) {
	xmi.start--;
	find_Lmaximal(P, prefix, SA[xmi.start], xmi.depth, matches, min_len, print);
      }
      // Find RMEMs to the right, check their left maximality.
      while(xmi.end+1 < N/K && LCP[xmi.end+1] >= xmi.depth) {
	xmi.end++;
	find_Lmaximal(P, prefix, SA[xmi.end], xmi.depth, matches, min_len, print);
      }
    }
  }
}


// Finds left maximal matches given a right maximal match at position i.
void sparseSA::find_Lmaximal(const string &P, long prefix, long i,
        long len, vector<match_t> &matches, int min_len, bool print) {
  // Advance to the left up to K steps.
  for(long k = 0; k < K; k++) {
    // If we reach the end and the match is long enough, print.
    if(prefix == 0 || i == 0) {
      if(len >= min_len) {
#ifndef NDEBUG
	if(print) print_match(match_t(i, prefix, len), matches);
#endif
	matches.push_back(match_t(i, prefix, len));
      }
      return; // Reached mismatch, done.
    }
    else if(P[prefix-1] != S[i-1]){
      // If we reached a mismatch, print the match if it is long enough.
      if(len >= min_len) {
#ifndef NDEBUG
	if(print) print_match(match_t(i, prefix, len), matches);
#endif
	matches.push_back(match_t(i, prefix, len));
      }
      return; // Reached mismatch, done.
    }
    prefix--; i--; len++; // Continue matching.
  }
}


// Print results in format used by MUMmer v3.  Prints results
// 1-indexed, instead of 0-indexed.
void sparseSA::print_match(const match_t m) {
  if(_4column == false) {
    printf("%8ld  %8ld  %8ld\n", m.ref + 1, m.query + 1, m.len);
  }
  else {
    long refseq=0, refpos=0;
    from_set(m.ref, refseq, refpos); // from_set is slow!!!
    // printf works faster than count... why? I don't know!!
    printf("  %s", descr[refseq].c_str());
    for(long j = 0; j < maxdescrlen - (long)descr[refseq].length() + 1; j++) putchar(' ');
    printf(" %8ld  %8ld  %8ld\n", refpos + 1L, m.query + 1L, m.len);
  }
}

// This version of print match places m_new in a buffer. The buffer is
// flushed if m_new.len <= 0 or it reaches 1000 entries.  Buffering
// limits the number of locks on cout.
void sparseSA::print_match(const match_t m_new, vector<match_t> &buf) {
  if(m_new.len > 0)  buf.push_back(m_new);
  if(buf.size() > 1000 || m_new.len <= 0) {
    pthread_mutex_lock(&cout_mutex);
    for(long i = 0; i < (long)buf.size(); i++) print_match(buf[i]);
    pthread_mutex_unlock(&cout_mutex);
    buf.clear();
  }
}

void sparseSA::print_match(const string meta, vector<match_t> &buf, bool rc) {
  pthread_mutex_lock(&cout_mutex);
  if(!rc) printf("> %s\n", meta.c_str());
  else printf("> %s Reverse\n", meta.c_str());
  for(long i = 0; i < (long)buf.size(); i++) print_match(buf[i]);
  pthread_mutex_unlock(&cout_mutex);
  buf.clear();
}

// Finds maximal almost-unique matches (MAMs) These can repeat in the
// given query pattern P, but occur uniquely in the indexed reference S.
void sparseSA::findMAM(const string &P, vector<match_t> &matches, int min_len, bool print) {
  interval_t cur(0, N-1, 0);
  long prefix = 0;
  while(prefix < (long)P.length()) {
    // Traverse SA top down until mismatch or full string is matched.
    traverse(P, prefix, cur, P.length());
    if(cur.depth <= 1) { cur.depth = 0; cur.start = 0; cur.end = N-1; prefix++; continue; }
    if(cur.depth >= min_len) {//cur.size() == 1 &&
        for(int i = 0; i < cur.size(); i++){
          if(is_leftmaximal(P, prefix, SA[cur.start+i])) {
            // Yes, it's a MAM.
            match_t m; m.ref = SA[cur.start+i]; m.query = prefix; m.len = cur.depth;
#ifndef NDEBUG
            if(print) print_match(m);
#endif
            matches.push_back(m);
          }
        }
    }
    do {
      cur.depth = cur.depth-1;
      cur.start = ISA[SA[cur.start] + 1];
      cur.end = ISA[SA[cur.end] + 1];
      prefix++;
      if( cur.depth == 0 || expand_link(cur) == false ) { cur.depth = 0; cur.start = 0; cur.end = N-1; break; }
    } while(cur.depth > 0 && cur.size() == 1);
  }
}

// Returns true if the position p1 in the query pattern and p2 in the
// reference is left maximal.
bool sparseSA::is_leftmaximal(const string &P, long p1, long p2) {
  if(p1 == 0 || p2 == 0) return true;
  else return P[p1-1] != S[p2-1];
}


struct by_ref { bool operator() (const match_t &a, const match_t &b) const { if(a.ref == b.ref) return a.len > b.len; else return a.ref < b.ref; }  };

// Maximal Unique Match (MUM)
void sparseSA::MUM(const string &P, vector<match_t> &unique, int min_len, bool print) {
  // Find unique MEMs.
  vector<match_t> matches;
  MAM(P, matches, min_len, false);

  // Adapted from Stephan Kurtz's code in cleanMUMcand.c in MUMMer v3.20.
  long currentright, dbright = 0;
  bool ignorecurrent, ignoreprevious = false;
  sort(matches.begin(), matches.end(), by_ref());
  for(long i = 0; i < (long)matches.size(); i++) {
    ignorecurrent = false;
    currentright = matches[i].ref + matches[i].len - 1;
    if(dbright > currentright)
      ignorecurrent = true;
    else {
      if(dbright == currentright) {
	ignorecurrent = true;
	if(!ignoreprevious && matches[i-1].ref == matches[i].ref)
	  ignoreprevious = true;
      }
      else {
	dbright = currentright;
      }
    }
    if(i > 0 && !ignoreprevious) {
#ifndef NDEBUG
      if(print)	print_match(matches[i-1]);
#endif
      unique.push_back(matches[i-1]);
    }
    ignoreprevious = ignorecurrent;
  }
  if(!ignoreprevious) {
    if(matches.size() > 0) {
#ifndef NDEBUG
      if(print) print_match(matches[matches.size()-1]);
#endif
      unique.push_back(matches[matches.size()-1]);
    }
  }

}

struct thread_data {
  vector<long> Kvalues; // Values of K this thread should process.
  sparseSA *sa; // Suffix array + aux informaton
  int min_len; // Minimum length of match.
  bool print;//verbose or not
  int maxCount;//max number of mems per query position
  const string *P; // Query string.
};

void *MEMthread(void *arg) {
  thread_data *data = (thread_data*)arg;
  vector<long> &K = data->Kvalues;
  sparseSA *sa = data->sa;

  // Find MEMs for all assigned offsets to this thread.

  vector<match_t> matches; // place-holder
  matches.reserve(2000);   // TODO: Use this as a buffer again!!!!!!

  for(long k = 0; k < (long)K.size(); k++) {  sa->findMEM(K[k], *(data->P), matches, data->min_len, true); }

  pthread_exit(NULL);
}

void sparseSA::MEM(const string &P, vector<match_t> &matches, int min_len, bool print, int num_threads) {
  if(min_len < K) return;
  if(num_threads == 1) {
    for(int k = 0; k < K; k++) { findMEM(k, P, matches, min_len, print); }
  }
  else if(num_threads > 1) {
    vector<pthread_t> thread_ids(num_threads);
    vector<thread_data> data(num_threads);

    // Make sure all num_threads are joinable.
    pthread_attr_t attr;  pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Distribute K-values evenly between num_threads.
    int t = 0;
    for(int k = 0; k < K; k++) {
      data[t].Kvalues.push_back(k);
      t++;
      if(t == num_threads) t = 0;
    }
    // Initialize additional thread data.
    for(int i = 0; i < num_threads; i++) { data[i].sa = this; data[i].min_len = min_len;  data[i].P = &P; }
    // Create joinable threads to find MEMs.
    for(int i = 0; i < num_threads; i++) pthread_create(&thread_ids[i], &attr, MEMthread, (void *)&data[i]);
    // Wait for all threads to terminate.
    for(int i = 0; i < num_threads; i++) pthread_join(thread_ids[i], NULL);
  }
}
// Use LCP information to locate right maximal matches. Test each for
// left maximality.
void sparseSA::collectSMAMs(const string &P, long prefix,
        const interval_t mli, interval_t xmi, vector<match_t> &matches, int min_len, int maxCount, bool print) {
  // All of the suffixes in xmi's interval are right maximal.
  //if(xmi.size() > maxCount ) return;// --> many long matches is ok, do not have to be unique!!!
  int upperLimit = xmi.size() < maxCount ? xmi.end : xmi.start + maxCount-1;
  for(long i = xmi.start; i <= upperLimit; i++) find_Lmaximal(P, prefix, SA[i], xmi.depth, matches, min_len, print);
}

// For a given offset in the prefix k, find all MEMs.
void sparseSA::findSMAM(long k, const string &P, vector<match_t> &matches, int min_len, int maxCount, bool print) {
  if(k < 0 || k >= K) { cerr << "Invalid k." << endl; return; }
  // Offset all intervals at different start points.
  long prefix = k;
  interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval

  // Right-most match used to terminate search.
  int min_lenK = min_len - (K-1);

  while( prefix <= (long)P.length() - (K-k)) {
    if(hasChild) 
        traverse_faster(P, prefix, xmi, P.length());
    else
        traverse(P, prefix, xmi, P.length());// Traverse until mismatch.
    if(xmi.depth <= 1) { xmi.reset(N/K-1); prefix+=K; continue; }
    if(xmi.depth >= min_lenK) {
        collectSMAMs(P, prefix, mli, xmi, matches, min_len, maxCount, print); // Using LCP info to find MEM length.
        // When using ISA/LCP trick, depth = depth - K. prefix += K.
        prefix+=K;
        if( !hasSufLink || suffixlink(xmi) == false ) { xmi.reset(N/K-1); continue; }
    }
    else {
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix+=K;
      if( !hasSufLink || suffixlink(xmi) == false ) { xmi.reset(N/K-1); continue; }
    }
  }
#ifndef NDEBUG
  if(print) print_match(match_t(), matches);   // Clear buffered matches.
#endif
}

void *SMAMthread(void *arg) {
  thread_data *data = (thread_data*)arg;
  vector<long> &K = data->Kvalues;
  sparseSA *sa = data->sa;

  // Find MEMs for all assigned offsets to this thread.

  vector<match_t> matches; // place-holder
  matches.reserve(2000);   // TODO: Use this as a buffer again!!!!!! Use this everywhere, but static initialize is even better

  for(long k = 0; k < (long)K.size(); k++) {  sa->findSMAM(K[k], *(data->P), matches, data->min_len, data->maxCount, data->print); }

  pthread_exit(NULL);
}

void sparseSA::SMAM(const string &P, vector<match_t> &matches, int min_len, int maxCount, bool print, int num_threads) {
  if(min_len < K) return;
  if(num_threads == 1) {
    for(int k = 0; k < K; k++) { findSMAM(k, P, matches, min_len, maxCount, print); }
  }
  else if(num_threads > 1) {
    vector<pthread_t> thread_ids(num_threads);
    vector<thread_data> data(num_threads);

    // Make sure all num_threads are joinable.
    pthread_attr_t attr;  pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Distribute K-values evenly between num_threads.
    int t = 0;
    for(int k = 0; k < K; k++) {
      data[t].Kvalues.push_back(k);
      t++;
      if(t == num_threads) t = 0;
    }
    // Initialize additional thread data.
    for(int i = 0; i < num_threads; i++) { data[i].sa = this; data[i].min_len = min_len;  data[i].P = &P; data[i].print = print; data[i].maxCount = maxCount; }
    // Create joinable threads to find MEMs.
    for(int i = 0; i < num_threads; i++) pthread_create(&thread_ids[i], &attr, SMAMthread, (void *)&data[i]);
    // Wait for all threads to terminate.
    for(int i = 0; i < num_threads; i++) pthread_join(thread_ids[i], NULL);
  }

}

bool compMatches(const match_t i, const match_t j){
    return (i.ref < j.ref || (i.ref == j.ref && i.query < j.query));
    //return (i.query < j.query || (i.query == j.query && i.ref < j.ref));
}

bool compMatchesQuery(const match_t i,const match_t j){
    //return (i.ref < j.ref || (i.ref == j.ref && i.query < j.query));
    return (i.query < j.query || (i.query == j.query && i.ref < j.ref));
}

bool compIntervals(const lis_t i,const lis_t j){
    return (i.len > j.len || (i.len == j.len && i.begin < j.begin));
}

void sparseSA::postProcess(vector<match_t> &matches){
    sort(matches.begin(),matches.end(), compMatches);
}

void sparseSA::inexactMatch(read_t & read,const align_opt & alnOptions, bool fwStrand, bool print){
#ifndef NDEBUG
    int loops = 0;
    int dps = 0;
    int dpsizeSum = 0;
    if(print) printf("read %s of length %d, strand %s\n",read.qname.c_str(),read.sequence.length(), fwStrand ? "FORWARD" : "REVERSE");
#endif
    string P = read.sequence;
    int Plength = P.length();
    int editDist = (int)(alnOptions.errorPercent*Plength)+1;
    if(!fwStrand)
        Utils::reverse_complement(P,false);
    bool clipping = !alnOptions.noClipping;
    int min_len = alnOptions.minMemLength;
    if(!alnOptions.fixedMinLength && Plength/editDist > min_len)
        min_len = Plength/editDist;
    const dp_scores & scores = alnOptions.scores;
    outputType outputT = ERRORSTRING;
    dp_output output;
    vector<match_t> matches;
    //calc seeds
    SMAM(P, matches, min_len, alnOptions.alignmentCount, false);
    if(alnOptions.tryHarder){//TODO: change try-harder to recalculate only after forward + reverse has been tried
        //easy solution: try reverse and if found: skip
        if(matches.empty())
            SMAM(P, matches, min_len, 1000, false);
        if(matches.empty())
            SMAM(P, matches, 20, 1000, false);
        if(matches.empty())
            SMAM(P, matches, Plength/editDist, 1000, false);
    }
#ifndef NDEBUG
    if(print) for(long index = 0; index < (long)matches.size(); index++) print_match(matches[index]);
    int matchesFound = matches.size();
    int matchesOverCount = 0;
    int matchCount[P.length()];
    for( int i = 0; i < P.length(); i++)
        matchCount[i] = 0;
    for( int i = 0; i < matches.size(); i++)
        matchCount[matches[i].query]++;
    for(int i = 0; i < P.length(); i++)
        if(matchCount[i] > alnOptions.alignmentCount)
            matchesOverCount+=matchCount[i]-alnOptions.alignmentCount;
#endif
    //sort matches
    if(matches.size()>0){
        postProcess(matches);//sort matches
        /////////////////////////
        //FIND CANDIDATE REGIONS
        /////////////////////////
        vector<lis_t> lisIntervals;
        int begin = 0;
        int end = 0;
        int len = 0;
        while(begin < matches.size()){
            int refEnd = matches[begin].ref - matches[begin].query + Plength + editDist;
            while(end < matches.size() && matches[end].ref + matches[end].len <= refEnd){
                len += matches[end].len;
                end++;
            }
            lisIntervals.push_back(lis_t(begin,end-1,len));
            begin = end;
            len = 0;
        }
        //sort candidate regions for likelyhood of an alignment
        sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
        //for every interval, try to align
        int alnCount = 0;
        int lisIndex = 0;
        //trial parameter for performance reasons
        int trial = 0;
        //////////////////////
        //MAIN ALIGNMENT LOOP
        //////////////////////
        while(alnCount < alnOptions.alignmentCount &&
                lisIndex < lisIntervals.size() &&
                lisIntervals[lisIndex].len > (Plength*alnOptions.minCoverage)/100 &&
                trial < alnOptions.maxTrial){
            //sort matches in this interval according to query position
            begin = lisIntervals[lisIndex].begin;
            end = lisIntervals[lisIndex].end;
#ifndef NDEBUG
            loops++;
            int dpSumInLoop = 0;
if(print){
    printf("region in reference to test: %ld + %d length\n",matches[begin].ref-matches[begin].query,Plength+editDist);
    printf("Number of mems is %d and number of bases covered by MEMs is %d\n",end-begin+1,lisIntervals[lisIndex].len);
    printf("This is interval %d of %d tested and there have been %d alignments found so far\n",lisIndex, lisIntervals.size(),alnCount);
    printf("filtered\n");
        for(long index = begin; index < end; index++) print_match(matches[index]);
}
#endif
            //sort this candidate region by query position
            sort(matches.begin()+begin,matches.begin()+end+1, compMatchesQuery);
            alignment_t alignment;
            int curEditDist = 0;
            ///////////////////
            //FIRST SEED
            //////////////////
            match_t firstSeed = matches[begin];
#ifndef NDEBUG
    if( print) printf("first seed with ref %ld, query %ld and length %ld\n",firstSeed.ref,firstSeed.query,firstSeed.len);
#endif
            int refstrLB = firstSeed.ref;
            int refstrRB = refstrLB + firstSeed.len -1;
            int queryLB = firstSeed.query;
            int queryRB = queryLB + firstSeed.len-1;
            output.clear();
            //no mem starting at 0 + fill in the starting position in the ref sequence
            if(firstSeed.query > 0 && firstSeed.ref > 0){
                //Alignment now!!! with beginQ-E in ref to begin match and beginQ to match
                int alignmentBoundLeft = refstrLB-queryLB-editDist;
                if(alignmentBoundLeft < 0)
                    alignmentBoundLeft = 0;
                boundaries grenzen(alignmentBoundLeft,refstrLB-1,0,queryLB-1);
                dp_type types;
                types.freeRefB = true;
                types.freeQueryB = clipping;
                dp( S, P, grenzen, scores, types, outputT, output, false);
                if(clipping && grenzen.queryB> 0){
                    alignment.cigarChars.push_back('S');
                    alignment.cigarLengths.push_back(grenzen.queryB);
                }
                alignment.cigarChars.insert(alignment.cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
                alignment.cigarLengths.insert(alignment.cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
                alignment.pos = grenzen.refB+1;
                curEditDist += output.editDist;
#ifndef NDEBUG
                    dps++;
                    dpsizeSum += ((refstrLB-alignmentBoundLeft)*(queryLB));
                    dpSumInLoop += ((refstrLB-alignmentBoundLeft)*(queryLB));
    if(print)printf("primary dp executed with editDist %d\n",output.editDist);
#endif
                output.clear();
            }//fill in the starting position in the ref sequence
            else {
                alignment.pos = firstSeed.ref + 1;
            }
            if (firstSeed.ref == 0 && firstSeed.query > 0) {
                alignment.cigarChars.push_back(clipping ? 'S' : 'I');
                alignment.cigarLengths.push_back(firstSeed.query);
                curEditDist += clipping ? 0 : firstSeed.query;
            }
            alignment.cigarChars.push_back('=');
            alignment.cigarLengths.push_back(firstSeed.len);
            //extend seed to the right, filling gaps between mems
            bool foundNext = true;
            int memsIndex = begin+1;
            //////////////////////
            //NEXT SEEDS LOOP
            //////////////////////
            while (memsIndex <= end && foundNext && curEditDist < editDist) {
                foundNext = false;
                //search matchingstats and calc cost
                int indexIncrease = 0;
                int minDistMem = memsIndex - 1;
                int minDist = P.length();
                int minQDist = P.length();
                int minRefDist = P.length();
                //temporary distance for this mem
                int dist = minDist;
                while (indexIncrease + memsIndex <= end && dist <= minDist) {
                    match_t match = matches[memsIndex + indexIncrease];
                    //check distance: k is enlarged and compare minRef-refstrRB
                    int qDist = match.query - queryRB;
                    int refDist = match.ref - refstrRB;
                    if (qDist < 0 || refDist < 0) {
                        if (qDist >= refDist) {
                            qDist -= refDist; //qDist insertions
                            if (qDist + queryRB >= P.length() || 1 - refDist >= match.len)//and refDist +1 matches lost
                                qDist = minDist + 1;
                        } else {
                            refDist -= qDist;
                            if (refDist + refstrRB >= S.length() || 1 - qDist >= match.len)
                                refDist = minDist + 1;
                        }
                    }
                    if (qDist < refDist)
                        dist = refDist;
                    else
                        dist = qDist;
                    if (dist <= minDist) {
                        minDist = dist;
                        minDistMem = indexIncrease + memsIndex;
                        minQDist = qDist;
                        minRefDist = refDist;
                        foundNext = true;
                    }
                    indexIncrease++;
                }
                if (foundNext) {
                    match_t match = matches[minDistMem];
#ifndef NDEBUG
                    if(print) printf("next match with ref %ld, query %ld and length %ld\n", match.ref, match.query, match.len);
#endif
                    //add indels and mutations to sol
                    if (minQDist <= 0) {
                        //minRefDist deletions
                        //1+|qDist| matches lost of next mem
                        alignment.cigarChars.push_back('D');
                        alignment.cigarChars.push_back('=');
                        alignment.cigarLengths.push_back(minRefDist);
                        alignment.cigarLengths.push_back(match.len - 1 + minQDist);
                        curEditDist += Utils::contains(S, refstrRB + 1, refstrRB + minRefDist - 1, '`') ? editDist + 1 : minRefDist; //boundary of ref sequences is passed
                        //TODO: only use of contains: add inline here
#ifndef NDEBUG
                       if(print) printf("extra match required deletion of %d chars, and edit increase of %d\n", minRefDist, Utils::contains(S, refstrRB + 1, refstrRB + minRefDist - 1, '`') ? editDist + 1 : minRefDist);
#endif
                    } else if (minRefDist <= 0) {
                        alignment.cigarChars.push_back('I');
                        alignment.cigarChars.push_back('=');
                        alignment.cigarLengths.push_back(minQDist);
                        alignment.cigarLengths.push_back(match.len - 1 + minRefDist);
                        curEditDist += minQDist;
#ifndef NDEBUG
     if(print)                   printf("extra match required insertion of %d chars, and edit increase of %d\n", minQDist, minQDist);
#endif
                    } else {//both distances are positive and not equal to (1,1)
                        boundaries grenzen(refstrRB + 1, refstrRB + minRefDist - 1, queryRB + 1, queryRB + minQDist - 1);
                        dp_type types;
                        dp(S, P, grenzen, scores, types, outputT, output, false);
                        alignment.cigarChars.insert(alignment.cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
                        alignment.cigarLengths.insert(alignment.cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
                        alignment.cigarChars.push_back('=');
                        alignment.cigarLengths.push_back(match.len);
                        curEditDist += output.editDist;
#ifndef NDEBUG
                       if(print) printf("extra match required dp, resulting in extra edit dist of %d\n", output.editDist);
                        dps++;
                        dpsizeSum += ((minRefDist)*(minQDist));
                        dpSumInLoop += ((minRefDist)*(minQDist));
#endif
                        output.clear();
                    }
                    //set new positions
                    memsIndex = minDistMem + 1;
                    //add matches
                    refstrRB = match.ref + match.len - 1;
                    queryRB = match.query + match.len - 1;
                }
            }
            //possible at the end characters
            //mapInterval += Integer.toString(queryRB+1)+"]";
            if (queryRB + 1 < P.length() && refstrRB < S.length() - 1 && curEditDist <= editDist) {
                int refAlignRB = refstrRB + (P.length() - queryRB) + 1 + editDist - curEditDist;
                if (refAlignRB >= S.length())
                    refAlignRB = S.length() - 1;
                boundaries grenzen(refstrRB + 1, refAlignRB, queryRB + 1, P.length() - 1);
                dp_type types;
                types.freeRefE = true;
                types.freeQueryE = clipping;
                dp(S, P, grenzen, scores, types, outputT, output, false);
                alignment.cigarChars.insert(alignment.cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
                alignment.cigarLengths.insert(alignment.cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
                if (clipping && grenzen.queryE < P.length() - 1) {
                    alignment.cigarChars.push_back('S');
                    alignment.cigarLengths.push_back(P.length() - 1 - grenzen.queryE);
                }
                curEditDist += output.editDist;
#ifndef NDEBUG
               if(print) printf("final dp required, resulting in extra edit dist of %d\n", output.editDist);
                dps++;
                dpsizeSum += ((refAlignRB-refstrRB)*(P.length()-queryRB-1));
                dpSumInLoop += ((refAlignRB-refstrRB)*(P.length()-queryRB-1));
#endif
                output.clear();
            } else if (queryRB + 1 < P.length() && curEditDist < editDist) {
                alignment.cigarChars.push_back(clipping ? 'S' : 'I');
                alignment.cigarLengths.push_back(P.length() - queryRB - 1);
                curEditDist += clipping ? 0 : P.length() - 1 - queryRB;
#ifndef NDEBUG
              if(print)  printf("final insertion required of %d chars\n", P.length() - queryRB - 1);
#endif
            }
#ifndef NDEBUG
            if(print){
    printf("alignment of %d edit dist found, while max is %d\n", curEditDist, editDist);
    alignment.setFieldsFromCigar(scores);
    printf("alignment cigar is %s\n", alignment.NMtag.c_str());
            }
#endif
    //TODO: check for possible optimizations (static initialisations
    //TODO: reorder matches according to best hits
            if(curEditDist <= editDist){
                if(!fwStrand){
                    alignment.flag.set(4,true);
                    if(alnOptions.unique && !read.alignments.empty() && curEditDist < read.alignments[0].editDist)
                        read.alignments.pop_back();
                }
                alnCount++;
                alignment.editDist = curEditDist;
                read.alignments.push_back(alignment);
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
    }
#ifndef NDEBUG
       if(print){       printf("%d;%d;%d;%d\n",loops,dps,dpsizeSum,matchesFound);
                        printf("%d\n",matchesOverCount);
       }
#endif
}


