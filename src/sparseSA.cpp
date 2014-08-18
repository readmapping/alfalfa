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

#include <stack>
#include <fstream> //printf
#include <math.h>//ceil
#include <stdlib.h>//atoi

#include "zlib.h"
#include "kseq.h"
#include "sparseSA.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

// LS suffix sorter (integer alphabet).
extern "C" { void suffixsort(int *x, int *p, int n, int k, int l); }

struct Less_intervTuple_t {
  bool operator()(const intervTuple_t& lhs, const intervTuple_t& rhs) const
  {
    return lhs.xmi.depth < rhs.xmi.depth;
  }
};

bool compIntervTuple(const intervTuple_t & i,const intervTuple_t & j){
    return (i.prefix+i.xmi.depth < j.prefix+j.xmi.depth || (i.prefix+i.xmi.depth == j.prefix+j.xmi.depth && i.prefix < j.prefix));
}

bool compIntervTupleRev(const intervTuple_t & i,const intervTuple_t & j){
    return (i.prefix+i.xmi.depth > j.prefix+j.xmi.depth || (i.prefix+i.xmi.depth == j.prefix+j.xmi.depth && i.prefix < j.prefix));
}

//////////////////////////
//struct vec_uchar
//////////////////////////

void vec_uchar::init(){ 
    sort(M.begin(), M.end()); 
    cerr << "number of lcp values >= 255: " << M.size() << endl; 
    vector<item_t>(M).swap(M);
}

long vec_uchar::index_size_in_bytes(){
    long indexSize = 0L;
    indexSize += sizeof(vec) + vec.capacity()*sizeof(unsigned char);
    indexSize += sizeof(M) + M.capacity()*(sizeof(size_t)+sizeof(int));
    return indexSize;
}

//////////////////////////
//struct sparseSA
//////////////////////////

long sparseSA::index_size_in_bytes(){
      long indexSize = 0L;
      indexSize += sizeof(sparseMult);
      indexSize += sizeof(hasSufLink);
      indexSize += sizeof(hasChild);
      indexSize += sizeof(K);
      indexSize += sizeof(NKm1);
      indexSize += sizeof(logN);
      indexSize += sizeof(N);
      indexSize += sizeof(maxdescrlen);
      indexSize += sizeof(descr);
      for(size_t i = 0; i < descr.size(); i++){
          indexSize += descr[i].capacity();
      }
      indexSize += sizeof(startpos) + startpos.capacity()*sizeof(long);
      indexSize += S.capacity();
      indexSize += sizeof(SA) + SA.capacity()*sizeof(unsigned int);
      indexSize += sizeof(ISA) + ISA.capacity()*sizeof(int);
      indexSize += sizeof(CHILD) + CHILD.capacity()*sizeof(int);
      indexSize += LCP.index_size_in_bytes();
      indexSize += sizeof(hasKmer);
      indexSize += sizeof(KMR) + KMR.capacity()*sizeof(saTuple_t);
      return indexSize;
}

void sparseSA::init(long K_, bool _hasSufLink, bool _hasChild, bool _hasKmer){
  hasSufLink = _hasSufLink;
  hasChild = _hasChild;
  hasKmer = _hasKmer;
  K = K_;

  // Increase string length so divisible by K.
  // Don't forget to count $ termination character.
  S.reserve(S.length()+2*K);
  if(S.length() % K != 0) {
    long appendK = K - S.length() % K ;
    for(long i = 0; i < appendK; i++) S += '$';
  }
  // Make sure last K-sampled characeter is this special character as well!!
  for(long i = 0; i < K; i++) S += '$'; // Append "special" end character. Note: It must be lexicographically less.
  N = S.length();
  std::string(S.data(), S.size()).swap(S);

  // Adjust to "sampled" size.
  logN = (long)ceil(log(N/K) / log(2.0));
  NKm1 = N/K-1;
}

void sparseSA::loadRef(const string & fileName){
    gzFile fp;
    kseq_t *seq;
    //first iteration for lengths (reserve in string)
    size_t Slength = 0;
    int l;
    fp = gzopen(fileName.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        Slength += seq->seq.l+1;
    }
    gzclose(fp);
    S.reserve(Slength);
    fp = gzopen(fileName.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
            printf("name: %s\n", seq->name.s);
            if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
            //to lower case
            for(size_t i = 0; i < seq->seq.l; i++){
                seq->seq.s[i] = std::tolower(seq->seq.s[i]);
            }
            string descrName(seq->name.s, seq->name.l);
            if(S.length() > 0)
                S += '`';
            startpos.push_back(S.length());
            S.append(seq->seq.s,seq->seq.l);
            descr.push_back(descrName);
    }
    kseq_destroy(seq);
    gzclose(fp);
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
    for(long i = 0; i < N/K; i++){
        CHILD[i] = -1;
    }
    //Compute up and down values
    long lastIndex  = -1;
    stack<long,vector<long> > stapelUD;
    stapelUD.push(0);
    for(long i = 1; i < N/K; i++){
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
    stack<long,vector<long> > stapelNL;
    stapelNL.push(0);
    for(long i = 1; i < N/K; i++){
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

// Look-up table construction algorithm
void sparseSA::computeKmer() {
    stack<interval_t> intervalStack;
    stack<unsigned int> indexStack;
    
    interval_t curInterval(0,N/K-1,0);
    unsigned int curIndex = 0;
    unsigned int newIndex = 0;
    
    intervalStack.push(interval_t(0,N/K-1,0));
    indexStack.push(curIndex);
    
    while(!intervalStack.empty()){
        curInterval = intervalStack.top(); intervalStack.pop();
        curIndex = indexStack.top(); indexStack.pop();
        if(KMERSIZE == (long) (curInterval.depth)){
            if(curIndex < TABLESIZE){
                KMR[curIndex].left = curInterval.start;
                KMR[curIndex].right = curInterval.end;
            }
        }
        else{
            if(hasChild){//similar to function traverse_faster
                //walk up to depth KMERSIZE or new child
                long curLCP; //max LCP of this interval
                if (curInterval.start == curInterval.end)//BUGFIX: was omitted first and caused a bug
                    curLCP = N - SA[curInterval.start];
                else if (curInterval.start < CHILD[curInterval.end] && CHILD[curInterval.end] <= curInterval.end)
                    curLCP = LCP[CHILD[curInterval.end]];
                else
                    curLCP = LCP[CHILD[curInterval.start]];
                long minimum = min(curLCP,KMERSIZE);
                newIndex = curIndex;
                while(curInterval.depth < (long) minimum){
                    newIndex = (newIndex << 2) | BITADD[(unsigned char) S[SA[curInterval.start]+curInterval.depth]];
                    curInterval.depth ++;
                }
                if(curInterval.depth == KMERSIZE){//reached KMERSIZE in the middle of an edge
                    if(newIndex < TABLESIZE){
                        KMR[newIndex].left = curInterval.start;
                        KMR[newIndex].right = curInterval.end;
                    }
                }
                else{//find child intervals
                    long left = curInterval.start;
                    long right = CHILD[curInterval.end];
                    curIndex = newIndex;
                    if(curInterval.start >= right || right > curInterval.end)
                        right = CHILD[curInterval.start];
                    //now left and right point to first child
                    newIndex = (curIndex << 2) | BITADD[(unsigned char) S[SA[left]+curInterval.depth]];
                    if(newIndex < TABLESIZE){
                        intervalStack.push(interval_t(left,right-1,curInterval.depth+1));
                        indexStack.push(newIndex);
                    }
                    left = right;
                    //while has next L-index
                    while(CHILD[right] > right && LCP[right] == LCP[CHILD[right]]){
                        right = CHILD[right];
                        newIndex = (curIndex << 2) | BITADD[ (unsigned char) S[SA[left]+curInterval.depth]];
                        if(newIndex < TABLESIZE){
                            intervalStack.push(interval_t(left,right-1,curInterval.depth+1));
                            indexStack.push(newIndex);
                        }
                        left = right;
                    }
                    //last interval
                    newIndex = (curIndex << 2) | BITADD[ (unsigned char) S[SA[left]+curInterval.depth]];
                    if(newIndex < TABLESIZE){
                        intervalStack.push(interval_t(left,curInterval.end,curInterval.depth+1));
                        indexStack.push(newIndex);
                    }
                }
            }
            else{//similar to function traverse
                long start = curInterval.start; long end = curInterval.end;
                while(start <= curInterval.end){//bug: did not calculate correctly
                    newIndex = (curIndex << 2) | BITADD[ (unsigned char) S[SA[start]+curInterval.depth]];
                    start = curInterval.start;
                    end = curInterval.end;
                    top_down_faster(S[SA[start]+curInterval.depth], curInterval.depth, start, end);
                    if(newIndex < TABLESIZE){
                        intervalStack.push(interval_t(start,end,curInterval.depth+1));
                        indexStack.push(newIndex);
                    }
                    // Advance to next interval.
                    start = end+1; end = curInterval.end;
                }
            }
        }
    }
    
}

//TODO: add error handling and messages
void sparseSA::save(const string &prefix){
    string basic = prefix;
    string aux = basic + ".aux";
    string sa = basic + ".sa";
    string lcp = basic + ".lcp";
    string ref = basic + ".ref";
    //print auxiliary information
    ofstream aux_s (aux.c_str());
    if (aux_s.is_open()){
        aux_s << N << endl;
        aux_s << K << endl;
        aux_s << logN << endl;
        aux_s << NKm1 << endl;
        aux_s << (hasSufLink ? 1 : 0) << endl;
        aux_s << (hasChild ? 1 : 0) << endl;
        aux_s << (hasKmer ? 1 : 0) << endl;
        aux_s << startpos.size() << endl;
        for(size_t i = 0; i < startpos.size(); i++)
            aux_s << startpos[i] << endl;
        for(size_t i = 0; i < descr.size(); i++)
            aux_s << descr[i] << endl;
        aux_s.close();
    }
    else cerr << "Unable to open " << aux.c_str() << endl;
    //print S
    ofstream ref_s (ref.c_str());
    if (ref_s.is_open()){
        ref_s << S << endl;
    }
    else cerr << "Unable to open " << ref.c_str() << endl;
    //print sa
    ofstream sa_s (sa.c_str(), ios::binary);
    unsigned int sizeSA = SA.size();
    sa_s.write((const char*)&sizeSA,sizeof(sizeSA));
    sa_s.write((const char*)&SA[0],sizeSA*sizeof(unsigned int));
    sa_s.close();
    //print LCP
    ofstream lcp_s (lcp.c_str(), ios::binary);
    unsigned int sizeLCP = LCP.vec.size();
    unsigned int sizeM = LCP.M.size();
    lcp_s.write((const char*)&sizeLCP,sizeof(sizeLCP));
    lcp_s.write((const char*)&sizeM,sizeof(sizeM));
    lcp_s.write((const char*)&LCP.vec[0],sizeLCP*sizeof(unsigned char));
    lcp_s.write((const char*)&LCP.M[0],sizeM*sizeof(vec_uchar::item_t));
    lcp_s.close();
    //print ISA if nec
    if(hasSufLink){
        string isa = basic + ".isa";
        ofstream isa_s (isa.c_str(), ios::binary);
        unsigned int sizeISA = ISA.size();
        isa_s.write((const char*)&sizeISA,sizeof(sizeISA));
        isa_s.write((const char*)&ISA[0],sizeISA*sizeof(int));
        isa_s.close();
    }
    //print child if nec
    if(hasChild){
        string child = basic + ".child";
        ofstream child_s (child.c_str(), ios::binary);
        unsigned int sizeCHILD = CHILD.size();
        child_s.write((const char*)&sizeCHILD,sizeof(sizeCHILD));
        child_s.write((const char*)&CHILD[0],sizeCHILD*sizeof(int));
        child_s.close();
    }
    //print child if nec
    if(hasKmer){
        string kmer = basic + ".kmer";
        ofstream kmer_s (kmer.c_str(), ios::binary);
        unsigned int sizeKMR = KMR.size();
        kmer_s.write((const char*)&sizeKMR,sizeof(sizeKMR));
        kmer_s.write((const char*)&KMR[0],sizeKMR*sizeof(saTuple_t));
        kmer_s.close();
    }
}

bool sparseSA::load(const string &prefix){
    cerr << "atempting to load index " << prefix  << " ... "<< endl;
    string basic = prefix;
    string aux = basic + ".aux";
    string sa = basic + ".sa";
    string lcp = basic + ".lcp";
    string ref = basic + ".ref";
    ifstream aux_s (aux.c_str());
    if(!aux_s.is_open()){
        cerr << "unable to open " << prefix << endl;
        return false;
    }
    string line;
    getline(aux_s,line); N = atol(line.c_str());
    getline(aux_s,line); K = atol(line.c_str());
    getline(aux_s,line); logN = atol(line.c_str());
    getline(aux_s,line); NKm1 = atol(line.c_str());
    getline(aux_s,line); hasSufLink = atoi(line.c_str()) == 1 ? true : false;
    getline(aux_s,line); hasChild = atoi(line.c_str()) == 1 ? true : false;
    getline(aux_s,line); hasKmer = atoi(line.c_str()) == 1 ? true : false;
    getline(aux_s,line); 
    int chromNumber = atoi(line.c_str());
    for(int i = 0; i < chromNumber; i++){
        getline(aux_s,line); 
        startpos.push_back(atol(line.c_str()));
    }
    for(int i = 0; i < chromNumber; i++){
        getline(aux_s,line); 
        descr.push_back(line);
    }
    //read auxiliary information
    aux_s.close();
    //read ref
    ifstream ref_s (ref.c_str());
    getline(ref_s,S);
    ref_s.close();
    //read sa
    ifstream sa_s (sa.c_str(), ios::binary);
    unsigned int sizeSA;
    sa_s.read((char*)&sizeSA,sizeof(sizeSA));
    SA.resize(sizeSA);
    sa_s.read((char*)&SA[0],sizeSA*sizeof(unsigned int));
    sa_s.close();
    //read LCP
    ifstream lcp_s (lcp.c_str(), ios::binary);
    unsigned int sizeLCP;
    unsigned int sizeM;
    lcp_s.read((char*)&sizeLCP,sizeof(sizeLCP));
    lcp_s.read((char*)&sizeM,sizeof(sizeM));
    LCP.vec.resize(sizeLCP);
    LCP.M.resize(sizeM);
    lcp_s.read((char*)&LCP.vec[0],sizeLCP*sizeof(unsigned char));
    lcp_s.read((char*)&LCP.M[0],sizeM*sizeof(vec_uchar::item_t));
    lcp_s.close();
    //read ISA if nec
    if(hasSufLink){
        string isa = basic + ".isa";
        ifstream isa_s (isa.c_str(), ios::binary);
        unsigned int sizeISA;
        isa_s.read((char*)&sizeISA,sizeof(sizeISA));
        ISA.resize(sizeISA);
        isa_s.read((char*)&ISA[0],sizeISA*sizeof(int));
        isa_s.close();
    }
    //read child if nec
    if(hasChild){
        string child = basic + ".child";
        ifstream child_s (child.c_str(), ios::binary);
        unsigned int sizeCHILD;
        child_s.read((char*)&sizeCHILD,sizeof(sizeCHILD));
        CHILD.resize(sizeCHILD);
        child_s.read((char*)&CHILD[0],sizeCHILD*sizeof(int));
        child_s.close();
    }
    //read kmer table if nec
    if(hasKmer){
        string kmer = basic + ".kmer";
        ifstream kmer_s (kmer.c_str(), ios::binary);
        unsigned int sizeKMR;
        kmer_s.read((char*)&sizeKMR,sizeof(sizeKMR));
        KMR.resize(sizeKMR);
        kmer_s.read((char*)&KMR[0],sizeKMR*sizeof(saTuple_t));
        kmer_s.close();
    }
    cerr << "index loaded succesful" << endl;
    return true;
}

void sparseSA::construct(){
    if(K > 1) {
        long bucketNr = 1;
        int *intSA = new int[N/K+1];  
        for(int i = 0; i < N/K; i++) intSA[i] = i; // Init SA.
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

    cerr << "size concatenated sequence = " << N << endl;

    LCP.resize(N/K);
    cerr << "number of suffix array values = " << N/K << endl;
    // Use algorithm by Kasai et al to construct LCP array.
    computeLCP();  // SA + ISA -> LCP
    LCP.init();
    if(!hasSufLink){
        {
          vector<int> tmp;
          ISA.swap(tmp);
        }
    }
    if(hasChild){
        CHILD.resize(N/K);
        //Use algorithm by Abouelhoda et al to construct CHILD array
        computeChild();
    }
    if(hasKmer){
        KMR.resize(TABLESIZE, saTuple_t());
        computeKmer();
    }

    NKm1 = N/K-1;
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
      int tmp = SA[ BucketBegin[ (unsigned char) S[ SA[pos]*K + h ] ] ];
      SA[ BucketBegin[ (unsigned char) S[ SA[pos]*K + h] ]++ ] = SA[pos];  // Move bucket beginning to the right, and replace
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

// Given SA interval apply binary search to match character c at
// position i in the search string. Adapted from the C++ source code
// for the wordSA implementation from the following paper: Ferragina
// and Fischer. Suffix Arrays on Words. CPM 2007.
bool sparseSA::top_down_faster(char c, long i, long &start, long &end) const {
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

//finds the child interval of cur that starts with character c
//updates left and right bounds of cur to child interval if found, or returns
//cur if not found (also returns true/false if found or not)
bool sparseSA::top_down_child(char c, interval_t &cur) const{
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

// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
void sparseSA::traverse(const string &P, int prefix, interval_t &cur, int min_len, int maxBranching) const {
  if(hasKmer && cur.depth == 0 && min_len >= KMERSIZE){//free match first bases
    unsigned int index = 0;
    for(long i = 0; i < KMERSIZE; i++)
        index = (index << 2 ) | BITADD[ (unsigned char) P[prefix + i]];
    if(index < TABLESIZE && KMR[index].right>0){
        cur.depth = KMERSIZE;
        cur.start = KMR[index].left;
        cur.end = KMR[index].right;
    }
    else{
        return;//this results in no found seeds where the first KMERSIZE bases contain a non-ACGT character
    }
  }
  if(cur.depth >= min_len || cur.size() <= maxBranching) return;

  while(prefix+cur.depth < (int)P.length()) {
    long start = cur.start; long end = cur.end;
    // If we reach a mismatch, stop.
    if(top_down_faster(P[prefix+cur.depth], cur.depth, start, end) == false) return;

    // Advance to next interval.
    cur.depth += 1; cur.start = start; cur.end = end;

    // If we reach min_len, stop.
    if(cur.depth == min_len || cur.size() <= maxBranching) return;
  }
}

// Traverse pattern P starting from a given prefix and interval
// until mismatch or min_len characters reached.
// Uses the child table for faster traversal
void sparseSA::traverse_faster(const string &P,const int prefix, interval_t &cur, int min_len, int maxBranching) const{
        if(hasKmer && cur.depth == 0 && min_len >= KMERSIZE){//free match first bases
            unsigned int index = 0;
            for(long i = 0; i < KMERSIZE; i++)
                index = (index << 2 ) | BITADD[ (unsigned char) P[prefix + i]];
            if(index < TABLESIZE && KMR[index].right>0){
                cur.depth = KMERSIZE;
                cur.start = KMR[index].left;
                cur.end = KMR[index].right;
            }
            else{
                return;//this results in no found seeds where the first KMERSIZE bases contain a non-ACGT character
            }
        }
        if(cur.depth >= (int) min_len || cur.size() <= maxBranching) return;
        size_t c = prefix + cur.depth;
        bool intervalFound = c < P.length();
        int curLCP; //check if this is correct for root interval (unlikely case)
        if (cur.start < CHILD[cur.end] && CHILD[cur.end] <= cur.end)
            curLCP = LCP[CHILD[cur.end]];
        else
            curLCP = LCP[CHILD[cur.start]];
        if (intervalFound && cur.size() > 1 && curLCP == cur.depth)
            intervalFound = top_down_child(P[c], cur);
        else if(intervalFound)
            intervalFound = P[c] == S[((long)SA[cur.start])+(long)cur.depth];
        bool mismatchFound = false;
        while(intervalFound && !mismatchFound && c < P.length() && 
                (cur.depth < min_len && cur.size() > maxBranching)){
            c++;
            cur.depth++;
            if(cur.start != cur.end){
                int childLCP;
                //calculate LCP of child node, which is now cur. the LCP value
                //of the parent is currently c - prefix
                if(cur.start < (long) CHILD[cur.end] && (long) CHILD[cur.end] <= cur.end)
                    childLCP = LCP[CHILD[cur.end]];
                else
                    childLCP = LCP[CHILD[cur.start]];
                int minimum = min(childLCP,min_len);
                //match along branch, note: no cur.size() restriction as this does not change 
                while(!mismatchFound && c < P.length() && cur.depth < minimum){
                    mismatchFound = S[((long)SA[cur.start])+(long)cur.depth] != P[c];
                    c++;
                    cur.depth += !mismatchFound;
                }
                intervalFound = c < P.length() && !mismatchFound &&
                        cur.depth < min_len && top_down_child(P[c], cur);
            }
            else{
                while(!mismatchFound && c < P.length() && cur.depth < min_len){
                    mismatchFound = (long) SA[cur.start]+cur.depth >= N ||
                            S[((long)SA[cur.start])+cur.depth] != P[c];
                    c++;
                    cur.depth += !mismatchFound;
                }
            }
        }
}

void sparseSA::traverse_keep(const string &P,const int prefix, vector<interval_t> &current, int min_len, int maxBranching) const{
    interval_t cur = current.back();
    if(cur.size() > maxBranching ) current.pop_back();
    if(cur.depth >= min_len || cur.size() <= -1) return;
    int c = prefix + cur.depth;
    bool intervalFound = c < (int) P.length();
    int curLCP; //check if this is correct for root interval (unlikely case)
    if (cur.start < CHILD[cur.end] && CHILD[cur.end] <= cur.end)
        curLCP = LCP[CHILD[cur.end]];
    else
        curLCP = LCP[CHILD[cur.start]];
    if (intervalFound && cur.size() > 1 && curLCP == cur.depth){
        intervalFound = top_down_child(P[c], cur);
        if(intervalFound && cur.size() <= maxBranching) current.push_back(cur);
    }
    else if(intervalFound)
        intervalFound = P[c] == S[((long)SA[cur.start])+cur.depth];
    bool mismatchFound = false;
    while(intervalFound && !mismatchFound &&
            c < (int) P.length() && (cur.depth < min_len && cur.size() > -1)){
        c++;
        cur.depth++;
        if(cur.start != cur.end){
            int childLCP;
            //calculate LCP of child node, which is now cur. the LCP value
            //of the parent is currently c - prefix
            if(cur.start < (long) CHILD[cur.end] && (long) CHILD[cur.end] <= cur.end)
                childLCP = LCP[CHILD[cur.end]];
            else
                childLCP = LCP[CHILD[cur.start]];
            int minimum = min(childLCP,min_len);
            //match along branch, note: no cur.size() restriction as this does not change 
            while(!mismatchFound && c < (int) P.length() && cur.depth < minimum){
                mismatchFound = S[((long)SA[cur.start])+cur.depth] != P[c];
                c++;
                cur.depth += !mismatchFound;
            }
            intervalFound = c < (int) P.length() && !mismatchFound &&
                    cur.depth < min_len && top_down_child(P[c], cur);
            if(intervalFound && cur.size() <= maxBranching) current.push_back(cur);
        }
        else{
            while(!mismatchFound && c < (int) P.length() && cur.depth < min_len){
                mismatchFound = (long) SA[cur.start]+cur.depth >= N ||
                        S[((long)SA[cur.start])+cur.depth] != P[c];
                c++;
                cur.depth += !mismatchFound;
            }
        }
    }
    if(!current.empty() && current.back().start == cur.start && current.back().end == cur.end && current.back().depth < cur.depth)
        current.back().depth = cur.depth;
}


// Suffix link simulation using ISA/LCP heuristic.
bool sparseSA::suffixlink(interval_t &m) const {
  m.depth -= K;
  if( m.depth <= 0) return false;
  m.start = ISA[SA[m.start] / K + 1];
  m.end = ISA[SA[m.end] / K + 1];
  return expand_link(m);
}

// Finds left maximal matches given a right maximal match at position i.
void sparseSA::find_LmaximalStrict(const string &P, int prefix, long i,
        int len, vector<match_t> &matches, int min_len, int border) const {
  // Advance to the left up to K steps.
  for(int k = 0; k < border; k++) {
    // If we reach the end and the match is long enough, print.
    if(prefix == 0 || i == 0 || P[prefix-1] != S[i-1]) {
        if(len >= min_len){
            matches.push_back(match_t(i, prefix, len));
        }
        return; // Reached mismatch, done.
    }
    prefix--; i--; len++; // Continue matching.
  }
}

void sparseSA::collectAllMEMs(const string &P, int prefix, const interval_t mli,
        interval_t xmi, vector<match_t> &matches, int min_len, int maxLeft) const {
  // All of the suffixes in xmi's interval are right maximal.
  long xmiEnd = xmi.end;
  for(long i = xmi.start; i <= xmiEnd; i++) find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);

  if(mli.start == xmi.start && mli.end == xmi.end){ return;}
  while(xmi.depth >= mli.depth) {
    // Attempt to "unmatch" xmi using LCP information.
    if(xmi.end+1 < N/K) xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
    else xmi.depth = LCP[xmi.start];

    // If unmatched XMI is > matched depth from mli, then examine rmems.
    if(xmi.depth >= mli.depth) {
      // Scan RMEMs to the left, check their left maximality..
      while(LCP[xmi.start] >= xmi.depth) {
	xmi.start--;
	find_LmaximalStrict(P, prefix, SA[xmi.start], xmi.depth, matches, min_len, maxLeft);
      }
      // Find RMEMs to the right, check their left maximality.
      while(xmi.end+1 < N/K && LCP[xmi.end+1] >= xmi.depth) {
	xmi.end++;
	find_LmaximalStrict(P, prefix, SA[xmi.end], xmi.depth, matches, min_len, maxLeft);
      }
    }
  }
}

void sparseSA::collectAllSMAMs(const string &P, int prefix,
        interval_t xmi, vector<match_t> &matches, int min_len, int maxCount, int maxLeft) const {
  // All of the suffixes in xmi's interval are right maximal.
  if(xmi.size() <= maxCount)
        for(long i = xmi.start; i <= xmi.end; i++) 
            find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
}

//Methode 3 collectmems
void sparseSA::collectXmiAndMli(const string &P, intervTuple_t interval, 
        vector<match_t> &matches, bool calcMli, int nextEndStop, 
        int nextEndStopPrefix, int min_len, int maxCount, int maxLeft) const {
    int prefix = interval.prefix;
    interval_t xmi = interval.xmi;
    interval_t mli = interval.mli;
    if(xmi.depth + prefix < min_len) return;
    //always calculate xmi MEMs
    for(long i = xmi.start; i <= xmi.end; i++)
        find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
    //if we need to calculate mli:
    if(!calcMli || (mli.start == xmi.start && mli.end == xmi.end)) return;
    //calc all mli until new end is smaller or equal to nextEndStop or ...
    if(xmi.end+1 < N/K) xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
    else xmi.depth = LCP[xmi.start];
    //iterate over all intervals
    while(xmi.size() < maxCount && prefix+xmi.depth > nextEndStop && xmi.depth >= mli.depth && xmi.depth + prefix >= min_len){
        long prevStart = xmi.start;
        long prevStop = xmi.end;
        // Scan RMEMs to the left, check their left maximality..
        long minStart = max(0L,xmi.start-(maxCount-xmi.size()));
        while(xmi.start >= minStart && LCP[xmi.start] >= xmi.depth)
            xmi.start--;
        // Find RMEMs to the right, check their left maximality.
        long maxEnd = min(N/K,xmi.end+1+(maxCount-xmi.size()));
        while(xmi.end+1 < maxEnd && LCP[xmi.end+1] >= xmi.depth)
            xmi.end++;
        if(xmi.size() < maxCount){
            for(long i = xmi.start; i < prevStart; i++)
                find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
            for(long i = prevStop+1; i <= xmi.end; i++)
                find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
            // Attempt to "unmatch" xmi using LCP information.
            if(xmi.end+1 < N/K) xmi.depth = max(LCP[xmi.start], LCP[xmi.end+1]);
            else xmi.depth = LCP[xmi.start];
        }
    }
    //last possible mli
    if(xmi.size() < maxCount && prefix+xmi.depth == nextEndStop && xmi.depth >= mli.depth && xmi.depth + prefix >= min_len && prefix <= nextEndStopPrefix+K){
        long prevStart = xmi.start;
        long prevStop = xmi.end;
        // Scan RMEMs to the left, check their left maximality..
        long minStart = max(0L,xmi.start-(maxCount-xmi.size()));
        while(xmi.start >= minStart && LCP[xmi.start] >= xmi.depth)
            xmi.start--;
        // Find RMEMs to the right, check their left maximality.
        long maxEnd = min(N/K,xmi.end+(maxCount-xmi.size()));
        while(xmi.end+1 < maxEnd && LCP[xmi.end+1] >= xmi.depth)
            xmi.end++;
        if(xmi.size() < maxCount){
            for(long i = xmi.start; i < prevStart; i++)
                find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
            for(long i = prevStop+1; i <= xmi.end; i++)
                find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
        }
    }
}

void sparseSA::collectIntervals(const string &P, int prefix, vector<interval_t> & intervals, 
  vector<match_t> &matches, bool calcMli, int nextEndStop, int min_len, int maxCount, int maxLeft) const{
    interval_t xmi = intervals.back();
    if(xmi.depth + prefix < min_len) return;
    //always calculate xmi MEMs
    for(long i = xmi.start; i <= xmi.end; i++)
        find_LmaximalStrict(P, prefix, SA[i], xmi.depth, matches, min_len, maxLeft);
    //if we need to calculate mli:
    if(!calcMli || intervals.size()==1) return;
    //iterate over all intervals
    int index = intervals.size()-2;
    
    while(index >= 0 && intervals[index].size() <= (maxCount)/K && prefix+intervals[index].depth > max((int)(nextEndStop-K),min_len)){
        interval_t mli = intervals[index];
        long prevStart = intervals[index+1].start;
        long prevStop = intervals[index+1].end;
        for(long i = mli.start; i < prevStart; i++)
            find_LmaximalStrict(P, prefix, SA[i], mli.depth, matches, min_len, maxLeft);
        for(long i = prevStop+1; i <= xmi.end; i++)
            find_LmaximalStrict(P, prefix, SA[i], mli.depth, matches, min_len, maxLeft);
        index--;
    }
}

// For a given offset in the prefix k, find all MEMs.
void sparseSA::findMEM(int k, const string &P, vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const {
  if(k < 0 || k >= K) { cerr << "Invalid k." << endl; return; }
  // Offset all intervals at different start points.
  int prefix = k;
  interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval

  // Right-most match used to terminate search.
  int min_lenK = min_len - (sparseQ * K - 1);

  while( prefix <= (int)P.length() - (min_lenK)) {//bugfix (was (K-k) before)
    if (hasChild)
        traverse_faster(P, prefix, mli, min_lenK, -1); // Traverse until minimum length matched.
    else
        traverse(P, prefix, mli, min_lenK, -1);    // Traverse until minimum length matched.
    if(mli.depth > xmi.depth) xmi = mli;
    if(mli.depth <= 1) { mli.reset(N/K-1); xmi.reset(N/K-1); prefix += sparseQ*K; continue; }

    if(mli.depth >= min_lenK) {
      if (hasChild)
          traverse_faster(P, prefix, xmi, P.length(), -1); // Traverse until mismatch.
      else
          traverse(P, prefix, xmi, P.length(), -1); // Traverse until mismatch.
      collectAllMEMs(P, prefix, mli, xmi, matches, min_len, sparseQ*K); // Using LCP info to find MEM length.
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix += sparseQ*K;
      if (!hasSufLink) {
        mli.reset(N / K - 1);
        xmi.reset(N / K - 1);
        continue;
      }
      else{
        int i = 0;
        bool succes = true;
        while (i < sparseQ && (succes = suffixlink(mli))) {
            suffixlink(xmi);
            i++;
        }
        if (!succes) {
            mli.reset(N / K - 1);
            xmi.reset(N / K - 1);
            continue;
        }
      }
    }
    else {
        // When using ISA/LCP trick, depth = depth - K. prefix += K.
        prefix += sparseQ*K;
        if (!hasSufLink) {
            mli.reset(N / K - 1);
            xmi.reset(N / K - 1);
            continue;
        } else {
            int i = 0;
            bool succes = true;
            while (i < sparseQ && (succes = suffixlink(mli))) {
                i++;
            }
            if (!succes) {
                mli.reset(N / K - 1);
                xmi.reset(N / K - 1);
                continue;
            }
            xmi = mli;
        }
    }
  }
}

// For a given offset in the prefix k, find all MEMs.
void sparseSA::findSMAM(int k, const string &P, vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const {
  if(k < 0 || k >= K) { cerr << "Invalid k." << endl; return; }
  // Offset all intervals at different start points.
  int prefix = k;
  interval_t xmi(0,N/K-1,0); // max match interval

  // Right-most match used to terminate search.
  int min_lenK = min_len - (sparseQ*K-1);
  int minMin_len = sparseQ > 1 ? min_len-K : min_len;
  while( prefix <= (int)P.length() - (min_lenK)) {//bugfix (was (K-k) before)
    if(hasChild) 
        traverse_faster(P, prefix, xmi, P.length(), -1);
    else
        traverse(P, prefix, xmi, P.length(), -1);// Traverse until mismatch.
    if(xmi.depth <= 1) { xmi.reset(N/K-1); prefix += sparseQ*K; continue; }
    if(xmi.depth >= min_lenK) {//less restriction seems to give best results
        collectAllSMAMs(P, prefix, xmi, matches, minMin_len, sparseQ*maxCount, sparseQ*K); // Using LCP info to find MEM length.
        // When using ISA/LCP trick, depth = depth - K. prefix += K.
        prefix += sparseQ*K;
        if (!hasSufLink) {
            xmi.reset(N / K - 1);
            continue;
        } else {
            int i = 0;
            bool succes = true;
            while (i < sparseQ && (succes = suffixlink(xmi)))
                i++;
            if (!succes) {
                xmi.reset(N / K - 1);
                continue;
            }
        }
    }
    else {
      // When using ISA/LCP trick, depth = depth - K. prefix += K.
      prefix += sparseQ*K;
      if (!hasSufLink) {
            xmi.reset(N / K - 1);
            continue;
      } else {
        int i = 0;
        bool succes = true;
        while (i < sparseQ && (succes = suffixlink(xmi)))
            i++;
        if (!succes) {
            xmi.reset(N / K - 1);
            continue;
        }
      }
    }
  }
}

void sparseSA::MEM(const string &P, vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const {
  if(min_len < K) return;
  for(int k = 0; k < K; k++) { findMEM(k, P, matches, min_len, maxCount, sparseQ); }
}

void sparseSA::SMAM(const string &P, vector<match_t> &matches, int min_len, int maxCount, int sparseQ) const {
  if(min_len < K) return;
  for(long k = 0; k < K; k++) { findSMAM(k, P, matches, min_len, maxCount, sparseQ); }
}

//Methode 3: achter naar voor en rekening houden met grootte van xmi, start, ...
void sparseSA::SMEM(const string &P, vector<match_t> &matches, int min_len, int maxCount, int sparseQ, int max_smem_cont) const {
  if(min_len < K) return;
  // Offset all intervals at different start points.
  int min_lenK = min_len - (sparseQ*K-1);
  interval_t mli(0,N/K-1,0); // min length interval
  interval_t xmi(0,N/K-1,0); // max match interval
  vector<intervTuple_t> intervals;
  intervals.reserve(P.length() - (min_lenK));
  for(int prefix = 0; prefix <= (int)P.length() - (min_lenK); prefix++) {
      // Right-most match used to terminate search.
      traverse_faster(P, prefix, mli, min_lenK, -1); // Traverse until minimum length matched.
      if(mli.depth > xmi.depth) xmi = mli;
      if(mli.depth <= 1) { mli.reset(N/K-1); xmi.reset(N/K-1); continue; }
      if(mli.depth >= min_lenK) {
          traverse_faster(P, prefix, xmi, P.length(), -1); // Traverse until mismatch.
          //check if we want to add this prefix to priority queue
          if(xmi.depth + prefix >= min_len){
              intervals.push_back(intervTuple_t(mli, xmi, prefix));
          }
      }
      mli.reset(N / K - 1);
      xmi.reset(N / K - 1);
      if(prefix % K == K-1) prefix += (sparseQ-1)*K;
  }
  sort(intervals.begin(),intervals.end(),compIntervTupleRev);
  vector<int> beginLong;
  vector<int> smemCount;
  beginLong.assign(P.length()+1,-1);
  smemCount.assign(P.length()+1,-1);
  //fill tables: for every RSMEM ending: make note if this position has a high number of RSMEMs.
  size_t i = 0;
  while(i < intervals.size()){
      long curSize = 0L;
      long curEnd = intervals[i].xmiEnd();
      bool foundMaxSmemCont = false;
      while(i < intervals.size() && intervals[i].xmiEnd() == curEnd){
          if(!foundMaxSmemCont || (intervals[i].prefix <= beginLong[curEnd] + K))
                curSize += intervals[i].xmi.size();
          if(!foundMaxSmemCont && curSize >= max_smem_cont){
              beginLong[intervals[i].xmiEnd()] = intervals[i].prefix;
              foundMaxSmemCont = true;
          }
          i++;
      }
      smemCount[intervals[i-1].xmiEnd()] = curSize;
  }
  i = 0;
  int curMinPos = P.length();
  int minMin_len = sparseQ > 1 ? min_len-K : min_len;
  while(i < intervals.size()){
      int curEnd = intervals[i].xmiEnd();
      //check if we need to skip all rsmems for this end position
      //do this only if there is a frequent rsmem with higher ending position and lower starting position
      bool skipEnd = smemCount[curEnd] >= maxCount || intervals[i].prefix > curMinPos+K;
      if(!skipEnd && smemCount[curEnd] >= max_smem_cont)
          curMinPos = min(curMinPos, beginLong[curEnd]);
      //iterate over this end pos
      while(i < intervals.size() && intervals[i].xmiEnd() == curEnd){
          //check if we need to collect SMEMs
          if(!skipEnd && intervals[i].prefix <= curMinPos+K && 
                  (smemCount[curEnd] < max_smem_cont || intervals[i].prefix <= beginLong[curEnd] + K )){
              //collect all smems in xmi
              //calculate stop criterium for mli
              //if xmi occurs frequent: do not calc mli
              if(smemCount[curEnd] >= max_smem_cont){// && (curPrefix + K <= beginLong[curEnd] || intervals[i].prefix > curPrefix + K)){
                  collectXmiAndMli(P, intervals[i], matches, false, 0L, 0L, minMin_len, maxCount, K*sparseQ);
              }
              else{
                  //if xmi does not occur frequent: calc mli until mem length reaches frequent smem
                  int nextEndStop = curEnd-1;
                  int nextEndPrefix = curEnd;
                  while(nextEndStop>=0 && smemCount[nextEndStop] < max_smem_cont)
                      nextEndStop--;
                  if(nextEndStop>=0)
                      nextEndPrefix = beginLong[nextEndStop];
                  collectXmiAndMli(P, intervals[i], matches, true, nextEndStop, nextEndPrefix, minMin_len, maxCount, K*sparseQ);
              }
          }
          i++;
      }
  }
}

//Methode 3: achter naar voor en rekening houden met grootte van xmi, start, ...
void sparseSA::SMEMP(const string &P, vector<match_t> &matches, int min_len, int maxCount, int sparseQ, int max_smem_cont) const {
  if(min_len < K) return;
  // Offset all intervals at different start points.
  int min_lenK = min_len - (sparseQ*K-1);
  interval_t mli(0,N/K-1,0); // min length interval
  vector<interval_t> intervals;
  int prevEnd = 0;
  int prevEndPrefix = 0;
  long prevSize = -1L;
  int minMin_len = sparseQ > 1 ? min_len-K : min_len;
  for(int prefix = 0; prefix <= (int)P.length() - (min_lenK); prefix++) {
      // Right-most match used to terminate search.
      intervals.clear();
      traverse_faster(P, prefix, mli, min_lenK, -1); // Traverse until minimum length matched.
      if(mli.depth <= 1) { mli.reset(N/K-1); continue; }
      if(mli.depth >= min_lenK) {
          intervals.push_back(mli);
          traverse_keep(P, prefix, intervals, P.length(), maxCount); // Traverse until mismatch.
          //decide if this interval should be extended
          if(!intervals.empty()){
              bool extendXmi = (prefix+intervals.back().depth) > prevEnd || 
              prevSize <= max_smem_cont || prefix <= prevEndPrefix+K || 
              intervals.back().size() <= max_smem_cont;
              bool extendMli = extendXmi && intervals.back().size() <= max_smem_cont && 
              ((prefix+intervals.back().depth) > prevEnd || prevSize <= max_smem_cont);
              //collect
              if(extendXmi)
                collectIntervals(P, prefix, intervals, matches, extendMli, prevEnd, minMin_len, maxCount, K*sparseQ);
              //update prevEnd and prevSize:
              if(prevEnd < prefix+intervals.back().depth){
                  prevEnd = prefix+intervals.back().depth;
                  prevSize = intervals.back().size();
                  prevEndPrefix = prevSize > max_smem_cont ? prefix : -1;
              }
              else{
                  prevSize += intervals.back().size();
                  if(prevEndPrefix < 0 && prevSize > max_smem_cont) prevEndPrefix = prefix;
              }
          }
      }
      mli.reset(N / K - 1);
      if(prefix % K == K-1) prefix += (sparseQ-1)*K;
  }
}
