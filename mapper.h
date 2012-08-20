/* 
 * File:   mapper.h
 * Author: mvyvermn
 *
 * Created on 2 augustus 2012, 16:18
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


struct alignment_t {
  alignment_t(): pnextGlob(0), concordant(false), globPos(0), pos(0), cigar("*"), flag(0), rname("*"), mapq(0), tLength(0),
  rnext("*"), pnext(0), editDist(0), alignmentScore(0), cigarChars(0), cigarLengths(0), NMtag("*"), refLength(0) {}
  alignment_t(const alignment_t & o): concordant(false), globPos(o.globPos), pos(o.pos), cigar(o.cigar), flag(o.flag.to_ulong()), rname(o.rname), mapq(o.mapq), tLength(o.tLength),
  rnext(o.rnext), pnext(o.pnext), pnextGlob(o.pnextGlob), editDist(o.editDist), alignmentScore(o.alignmentScore), cigarChars(o.cigarChars), cigarLengths(o.cigarLengths), NMtag(o.NMtag), refLength(o.refLength) {}
  string cigar;//TODO: remove this fields, only used when printed
  string NMtag;//TODO: remove this fields, only used when printed
  vector<char> cigarChars;//Change these to fixed-length values
  vector<int> cigarLengths;//Change these to fixed-length values (make sure to increase them when necessary): make own string and vector classes
  string rname;//leave out, only for printing
  long globPos; // position in index (concat of all reference sequences)
  long pos; // position in reference sequence
  bitset<11> flag;
  int mapq;
  //TODO: paired-end mapping
  string rnext;
  long pnext;
  long pnextGlob;
  int editDist;
  int alignmentScore;
  int tLength;
  int refLength;//length in reference sequence spanned by this alignment
  bool concordant;
  void setLocalPos(const sparseSA& sa){
      if(rname.empty() || rname == "*"){
        long descIndex;
        sa.from_set(globPos, descIndex, pos);
        rname = sa.descr[descIndex];
      }
  }
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
};

// interval in match results + bases covering the result
struct lis_t {
  lis_t(): begin(0), end(0), len(0), fw(1), alignment(NULL) {}
  lis_t(vector<match_t> * matches, int b, int e, int l, bool fw): matches(matches), begin(b), end(e), len(l), fw(fw) { alignment = NULL;}
  int begin; // position in reference sequence
  int end; // position in query
  int len; // length of match
  bool fw;
  vector<match_t> * matches;
  alignment_t * alignment;
};

struct read_t {
    read_t(): qname(""),sequence(""),qual("*"),rcSequence(""),rQual(""), alignments(0) {}
    read_t(string name, string &seq, string &qw, bool nucleotidesOnly): qname(name),sequence(seq),qual(qw), alignments(0){
        rcSequence = sequence;//copy
        Utils::reverse_complement(rcSequence, nucleotidesOnly);
        rQual = qual;
        reverse(rQual.begin(),rQual.end());
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
            alignments[j].setLocalPos(sa);
             if(alignments[j].alignmentScore == maxScore)
                 alignments[j].mapq = mapq;
             else
                 alignments[j].flag.set(8,true);
        }
    }
    string emptyAlingment(bool paired, bool mateFailed){
        stringstream * ss = new stringstream;
        int flag = 4;
        if(paired) flag += 1;
        if(mateFailed) flag += 8;
        *ss << qname << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << sequence << "\t" << qual << endl;
        return ss->str();
    }
    int alignmentCount(){
        return alignments.size();
    }
    string printAlignment(int i){
        stringstream * ss = new stringstream;
        alignment_t & a = alignments[i];
        *ss << qname << "\t" << a.flag.to_ulong() << "\t"
            << a.rname << "\t" << a.pos << "\t" << 
            a.mapq << "\t" << a.cigar << "\t" << 
            a.rnext << "\t" << a.pnext << "\t" << 
            a.tLength << "\t" << (a.flag.test(4) ? rcSequence : sequence) <<
            "\t" << (a.flag.test(4) ? rQual : qual) << "\tAS:i:" << 
            a.alignmentScore << "\tNM:i:" << a.editDist << "\tX0:Z:" << 
            a.NMtag << endl;
        return ss->str();
    }
    
    void printAlignments(){
        if(alignments.empty()){
            printf("%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\n",qname.c_str(),4,"*",0,0,"*","*",0,0,sequence.c_str(),qual.c_str());
        }
        else{
            alignment_t & a = alignments[0];
            printf("%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\tAS:i:%d\tNH:i:%ld\tNM:i:%d\tX0:Z:%s\n",
                        qname.c_str(),(int)a.flag.to_ulong(),a.rname.c_str(),a.pos,a.mapq,a.cigar.c_str(),
                        a.rnext.c_str(),a.pnext,a.tLength, a.flag.test(4) ? rcSequence.c_str() : sequence.c_str(), a.flag.test(4) ? rQual.c_str() : qual.c_str(),
                        a.alignmentScore,alignments.size(),a.editDist,a.NMtag.c_str());
            for(int i=1; i < alignments.size(); i++){
                a = alignments[i];
                printf("%s\t%d\t%s\t%ld\t%d\t%s\t%s\t%ld\t%d\t%s\t%s\tAS:i:%d\tNM:i:%d\tX0:Z:%s\n",
                        qname.c_str(),(int)a.flag.to_ulong(),a.rname.c_str(),a.pos,a.mapq,a.cigar.c_str(),
                        a.rnext.c_str(),a.pnext,a.tLength,a.flag.test(4) ? rcSequence.c_str() : sequence.c_str(),
                        a.flag.test(4) ? rQual.c_str() : qual.c_str(),a.alignmentScore,a.editDist,a.NMtag.c_str());
            }
        }
    }
    string qname;//TODO should be reference
    vector<alignment_t> alignments;
    string sequence;//TODO should be reference
    string rcSequence;
    string qual;//TODO should be reference
    string rQual;
};

extern void unpairedMatch(const sparseSA& sa, dynProg& dp_, read_t & read,const align_opt & alnOptions, bool print);
extern void inexactMatch(const sparseSA& sa, dynProg& dp_, read_t& read, const align_opt & alnOptions, bool fwStrand, bool print);
//TODO: calculate global position in above function

//PAIRED END FUNCTIONS
extern bool isConcordant(const alignment_t& mate1, const alignment_t& mate2, const paired_opt& options);
extern void pairedMatch(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt, int mode, bool print);

#endif	/* MAPPER_H */
