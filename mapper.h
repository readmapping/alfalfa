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

// interval in match results + bases covering the result
struct lis_t {
  lis_t(): begin(0), end(0), len(0), fw(1) {}
  lis_t(vector<match_t> * matches, int b, int e, int l, bool fw): matches(matches), begin(b), end(e), len(l), fw(fw) {}
  int begin; // position in reference sequence
  int end; // position in query
  int len; // length of match
  bool fw;
  vector<match_t> * matches;
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
  void setGlobalPos(const sparseSA& sa){
      if(rname.empty() || rname == "*"){
        long globPos = pos;
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
            alignments[j].setGlobalPos(sa);
             if(alignments[j].alignmentScore == maxScore)
                 alignments[j].mapq = mapq;
             else
                 alignments[j].flag.set(8,true);
        }
    }
    string emptyAlingment(){
        stringstream * ss = new stringstream;
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

extern void unpairedMatch(const sparseSA& sa, read_t & read,const align_opt & alnOptions, bool print);
extern void inexactMatch(const sparseSA& sa, read_t& read, const align_opt & alnOptions, bool fwStrand, bool print);
//TODO: calculate global position in above function

#endif	/* MAPPER_H */

