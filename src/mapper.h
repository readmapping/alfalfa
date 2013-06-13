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
 */

#ifndef MAPPER_H
#define	MAPPER_H

#include <string>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <bitset>

#include "sparseSA.h"
#include "utils.h"
#include "options.h"

using namespace std;

struct mate_t {
    mate_t(): flag(0), rnext("*"), pnext(0), pnextGlob(0), tLength(0), concordant(false) {}
    mate_t(const mate_t & o): flag(o.flag.to_ulong()), rnext(o.rnext), pnext(o.pnext), 
    pnextGlob(o.pnextGlob), tLength(o.tLength), concordant(o.concordant) {}
    bitset<11> flag;
    string rnext;
    long pnext;
    long pnextGlob;
    int tLength;
    bool concordant;
};

struct alignment_t {
  alignment_t(): globPos(0), pos(0), cigar("*"), flag(0), rname("*"), mapq(0), editDist(0), 
  alignmentScore(0), cigarChars(0), cigarLengths(0), NMtag("*"), refLength(0), mateInfo(0) {}
  alignment_t(const alignment_t & o): globPos(o.globPos), pos(o.pos), cigar(o.cigar), flag(o.flag.to_ulong()), 
  rname(o.rname), mapq(o.mapq), editDist(o.editDist), alignmentScore(o.alignmentScore), 
  cigarChars(o.cigarChars), cigarLengths(o.cigarLengths), NMtag(o.NMtag), 
  refLength(o.refLength), mateInfo(o.mateInfo) {}
  long globPos; // position in index (concat of all reference sequences)
  long pos; // position in reference sequence
  string cigar;//TODO: remove this fields, only used when printed
  bitset<11> flag;
  string rname;//leave out, only for printing
  int mapq;
  int editDist;
  int alignmentScore;
  vector<char> cigarChars;//Change these to fixed-length values
  vector<int> cigarLengths;//Change these to fixed-length values (make sure to increase them when necessary): make own string and vector classes
  string NMtag;//TODO: remove this fields, only used when printed
  int refLength;//length in reference sequence spanned by this alignment
  //paired-end fields
  vector<mate_t> mateInfo;
  //functions
  bool paired(){
      return !mateInfo.empty();
  }
  bool concordant() const{
      return !mateInfo.empty() && mateInfo[0].concordant;
  }
  void setLocalPos(const sparseSA& sa){
      if(rname.empty() || rname == "*"){
        long descIndex;
        sa.from_set(globPos, descIndex, pos);
        rname = sa.descr[descIndex];
      }
  }
  int pairedCount(){
      return mateInfo.size();
  }
  void createCigar(bool newVersion){
      stringstream ss;
      assert(cigarChars.size()==cigarLengths.size());
      if(newVersion)
          for(size_t i = 0; i < cigarChars.size(); i++)
              ss << cigarLengths[i] << cigarChars[i];
      else{
          size_t i = 0;
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
      size_t i = 0;
      while(i < cigarChars.size()){
          if(cigarChars[i]=='=' || cigarChars[i]=='X'){
                int tempLength = cigarLengths[i];
                sNM << cigarLengths[i] << cigarChars[i];
                alignmentScore += (cigarChars[i]=='X' ? scores.mismatch*cigarLengths[i] : scores.match*cigarLengths[i]);
                while(i < cigarChars.size()-1 && (cigarChars[i+1]=='=' || cigarChars[i+1]=='X')){
                    i++;
                    tempLength += cigarLengths[i];
                    sNM << cigarLengths[i] << cigarChars[i];
                    alignmentScore += (cigarChars[i]=='X' ? scores.mismatch*cigarLengths[i] : scores.match*cigarLengths[i]);
                }
                sCig << tempLength << 'M';
          }
          else{
            sCig << cigarLengths[i] << cigarChars[i];
            sNM << cigarLengths[i] << cigarChars[i];
            if(cigarChars[i]=='D' || cigarChars[i]=='I'){
                alignmentScore += scores.openGap + scores.extendGap*cigarLengths[i];
            }
          }
          i++;
      }
      cigar = sCig.str();
      NMtag = sNM.str();
  }
  void addMate(const alignment_t * o, bool concordant, bool upstream){
      mate_t mate;
      mate.flag.set(0,true);//positions for mate to be important: 0,1,5,6,7
      mate.flag.set(1,true);
      mate.flag.set(5,o->flag.test(4));
      mate.flag.set(6,upstream);
      mate.flag.set(7,!upstream);
      mate.concordant = concordant;
      mate.pnext = o->pos;
      mate.pnextGlob = o->globPos;
      if(rname == o->rname){
        mate.rnext = "=";
        mate.tLength = max(pos+(long)refLength-1L,o->pos+(long)o->refLength-1L) - min(pos,o->pos) + 1;
        long firstPos = flag.test(4) ? pos+(long)refLength-1L : pos;
        long secondPos = o->flag.test(4) ? o->pos+(long)o->refLength-1L : o->pos;
        if(firstPos > secondPos) mate.tLength *= -1;
      }
      else{
        mate.rnext = o->rname;
        mate.tLength = 0;
      }
      mateInfo.push_back(mate);
  }
};

// interval in match results + bases covering the result
struct lis_t {
  lis_t(): begin(0), end(0), len(0), fw(1), alignment(NULL), extended(0) {}
  lis_t(vector<match_t> * matches, int b, int e, int l, bool fw_): begin(b), end(e), len(l), fw(fw_), alignment(NULL), extended(0), matches(matches) {}
  int begin; // position in reference sequence
  int end; // position in query
  int len; // length of match
  bool fw;
  alignment_t * alignment;
  bool extended;
  vector<match_t> * matches;
};

struct read_t {
    read_t(): qname(""),sequence(""),qual("*"),rcSequence(""),rQual(""), alignments(0),pairedAlignmentCount(0) {}
    read_t(string name, string &seq, string &qw, bool nucleotidesOnly): qname(name),
    sequence(seq),qual(qw), alignments(0), pairedAlignmentCount(0){
        rcSequence = sequence;//copy
        Utils::reverse_complement(rcSequence, nucleotidesOnly);
        rQual = qual;
        reverse(rQual.begin(),rQual.end());
    }
    ~read_t() {
        for(size_t i=0; i < alignments.size(); i++){
            delete alignments[i];
        }
    }
    void init(bool nucleotidesOnly){
        rcSequence = sequence;//copy
        Utils::reverse_complement(rcSequence, nucleotidesOnly);
        rQual = qual;
        reverse(rQual.begin(),rQual.end());
    }
    void postprocess(const dp_scores & scores, const sparseSA& sa){
        int maxScore = scores.mismatch*sequence.length();
        assert(maxScore < 0);//works only for score<0 for now (to do: add sign switch to allow positive min scores)
        int secBestScore = scores.mismatch*sequence.length();
        for(size_t j = 0; j < alignments.size(); j++){
             alignments[j]->setFieldsFromCigar(scores);
             if(alignments[j]->alignmentScore > maxScore){
                 secBestScore = maxScore;
                 maxScore = alignments[j]->alignmentScore;
             }
             else if(alignments[j]->alignmentScore < maxScore && alignments[j]->alignmentScore > secBestScore)
                 secBestScore = alignments[j]->alignmentScore;
        }
        int mapq = (secBestScore == scores.mismatch*sequence.length()) ? 255 :
            250*(maxScore-secBestScore)/(maxScore-scores.mismatch*sequence.length()) ;//signed-unsigned comparison: fix?
        assert(mapq <= 255 && mapq >= 0);
        for(size_t j = 0; j < alignments.size(); j++){
            alignments[j]->setLocalPos(sa);
             if(alignments[j]->alignmentScore == maxScore)
                 alignments[j]->mapq = mapq;
             else
                 alignments[j]->flag.set(8,true);
        }
    }
    string emptyAlingment(bool paired, bool mateFailed, bool upstream){
        stringstream ss;
        int flag = 4;
        if(paired) flag += 1;
        if(mateFailed) flag += 8;
        if(paired && upstream) flag += 64;
        if(paired && !upstream) flag += 128;
        ss << qname << "\t" << flag << "\t*\t0\t0\t*\t*\t0\t0\t" << sequence << "\t" << qual << endl;
        return ss.str();
    }
    int alignmentCount(){//check the correct use of this
        return alignments.size();
    }
    void removeLastAlignment(){
        delete alignments[alignments.size()-1];
        alignments.pop_back();
    }
    string printUnpairedAlignment(int i){
        stringstream ss;
        alignment_t * a = alignments[i];
        ss << qname << "\t" << a->flag.to_ulong() << "\t"
            << a->rname << "\t" << a->pos << "\t" << 
            a->mapq << "\t" << a->cigar << "\t" << 
            "*" << "\t" << 0 << "\t" << 
            0 << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
            "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" << 
            a->alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" << 
            a->NMtag << endl;
        return ss.str();
    }
    string printPairedAlignments(int i){
        stringstream ss;
        alignment_t * a = alignments[i];
        if(a->paired()){
            for(int j = 0; j < a->pairedCount(); j++){
                //flag has to be changed
                ss << qname << "\t" << (a->flag.to_ulong() |a->mateInfo[j].flag.to_ulong())  << "\t"
                << a->rname << "\t" << a->pos << "\t" << 
                a->mapq << "\t" << a->cigar << "\t" << 
                a->mateInfo[j].rnext << "\t" << a->mateInfo[j].pnext << "\t" << 
                a->mateInfo[j].tLength << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
                "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" << 
                a->alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" << 
                a->NMtag << endl;
            }
            return ss.str();
        }
        else{
            return printUnpairedAlignment(i);
        }
    }
    string qname;//TODO should be reference
    string sequence;//TODO should be reference
    string qual;//TODO should be reference
    string rcSequence;
    string rQual;
    vector<alignment_t *> alignments;
    int pairedAlignmentCount;
};

extern void unpairedMatch(const sparseSA& sa, dynProg& dp_, read_t & read,const align_opt & alnOptions);
extern void inexactMatch(const sparseSA& sa, dynProg& dp_, read_t& read, const align_opt & alnOptions, bool fwStrand);
//TODO: calculate global position in above function

//PAIRED END FUNCTIONS
//extern bool isConcordant(const alignment_t& mate1, const alignment_t& mate2, const paired_opt& options);
extern void pairedMatch(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt);

#endif	/* MAPPER_H */

