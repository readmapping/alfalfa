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
#include <bitset>
#include <vector>
#include <stdint.h>
#include <algorithm> //min,max
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include <iostream>

#include "sparseSA.h"
#include "options.h"
#include "utils.h"

static bool compMatchesQuery(const match_t & i,const match_t & j){
    return (i.query < j.query || (i.query == j.query && i.ref < j.ref));
}

// interval in match results + bases covering the result
struct lis_t {
  lis_t(): score(0), fw(1), refBegin(0L), refEnd(0L), extended(0) {}
  lis_t(bool fw_): score(0), fw(fw_), refBegin(0L), refEnd(0L), extended(0) {}
  lis_t(const lis_t & o): score(o.score), fw(o.fw), refBegin(o.refBegin), refEnd(o.refEnd), extended(o.extended), matches(o.matches){}
  ~lis_t(){ matches.clear(); }
  int score; // score of the region
  bool fw;
  long refBegin;
  long refEnd;
  bool extended;
  std::vector<match_t> matches;
  
  void addSeed(const match_t & seed){matches.push_back(seed); }
  
  void addSeeds(const std::vector<match_t> & seeds, size_t begin, size_t end){
      for(size_t i = begin; i <= end; i++)
          matches.push_back(seeds[i]);
  }

  void sortByReadPos(){
      std::sort(matches.begin(),matches.end(), compMatchesQuery);
  }
  
};

struct mate_t {
    mate_t(): flag(0), rnext("*"), pnext(0), pnextGlob(0), tLength(0), concordant(false), alignmentScore(0) {}
    mate_t(const mate_t & o): flag(o.flag.to_ulong()), rnext(o.rnext), pnext(o.pnext), 
    pnextGlob(o.pnextGlob), tLength(o.tLength), concordant(o.concordant), alignmentScore(o.alignmentScore) {}
    std::bitset<11> flag;
    std::string rnext;
    long pnext;
    long pnextGlob;
    int tLength;
    bool concordant;
    int alignmentScore;
};

struct alignment_t {
  alignment_t(): refB(0), refE(0), qB(0), qE(0), pos(0), cigar("*"), flag(0), rname("*"), mapq(0), editDist(0), 
  alignmentScore(0), NMtag("*"), refLength(0), mate(), chain(0) {}
  alignment_t(const alignment_t & o): refB(o.refB), refE(o.refE), qB(o.qB), qE(o.qE), pos(o.pos), cigar(o.cigar), flag(o.flag.to_ulong()), 
  rname(o.rname), mapq(o.mapq), editDist(o.editDist), alignmentScore(o.alignmentScore),
  NMtag(o.NMtag), refLength(o.refLength), mate(o.mate), chain(o.chain) {}
  long refB, refE; // position in index (concat of all reference sequences)
  int qB, qE;
  long pos; // position in reference sequence
  std::string cigar;//TODO: remove this fields, only used when printed
  std::bitset<11> flag;
  std::string rname;//leave out, only for printing
  int mapq;
  int editDist;
  int alignmentScore;
  std::string NMtag;//TODO: remove this fields, only used when printed
  int refLength;//length in reference sequence spanned by this alignment
  //paired-end fields
  mate_t mate;
  std::vector<match_t> chain;
  //functions
  bool paired() const { return flag.test(0); }
  bool concordant() const{ return paired() && mate.concordant; }
  
  void setLocalPos(const sparseSA& sa){
    if(rname.empty() || rname == "*"){
        long descIndex;
        refB++;
        sa.from_set(refB, descIndex, pos);
        refB--;
        rname = sa.descr[descIndex];
    }
  }
  
  void setUnPaired(bool upstream){
    flag.set(0,true);
    flag.set(6,upstream);
    flag.set(7,!upstream);
    flag.set(3,true);
  }
  
  void setFieldsKsw(int n_cigar, uint32_t *cigar, const uint8_t * query, const uint8_t * rseq, int beginGap, int endGap);
  void constructCIGAR(const std::string& S, const uint8_t * readP, int Plen, int fixedBand, const int8_t *mat, int match, int openGap, int extendGap, int print);
  void setMate(const alignment_t * o, bool concordant, bool upstream);
};

struct read_t {
    read_t(): qname(""),sequence(""),qual("*"),rcSequence(""),rQual(""), alignments(0),pairedAlignmentCount(0) {}
    read_t(std::string name, std::string &seq, std::string &qw, bool nucleotidesOnly): qname(name),
    sequence(seq),qual(qw), alignments(0), pairedAlignmentCount(0){
        rcSequence = sequence;//copy
        Utils::reverse_complement(rcSequence, nucleotidesOnly);
        rQual = qual;
        std::reverse(rQual.begin(),rQual.end());
    }
    ~read_t() {
        if(readP != NULL){
            std::free(readP); readP = NULL;
            std::free(rcReadP); rcReadP = NULL;
        }
        for(size_t i=0; i < alignments.size(); i++){
            delete alignments[i];
        }
    }
    
    std::string qname;//TODO should be reference
    std::string sequence;//TODO should be reference
    std::string qual;//TODO should be reference
    std::string rcSequence;
    std::string rQual;
    uint8_t *readP;
    uint8_t *rcReadP;
    std::vector<alignment_t *> alignments;
    int pairedAlignmentCount;
    
    void init(bool nucleotidesOnly){
        rcSequence = sequence;//copy
        Utils::reverse_complement(rcSequence, nucleotidesOnly);
        rQual = qual;
        std::reverse(rQual.begin(),rQual.end());
        readP = (uint8_t *)std::malloc(sequence.length());
        rcReadP = (uint8_t *)std::malloc(sequence.length());
        for (size_t i = 0; i < sequence.length(); ++i){ // convert to the nt4 encoding
            readP[i] = Utils::ORDVALUE[(size_t)sequence[i]];
            rcReadP[i] = Utils::ORDVALUE[(size_t)rcSequence[i]];
        }
    }
    
    void postprocess(const sparseSA& sa, const align_opt & alnOptions, bool paired);
    
    std::string emptyAlingment(bool paired, bool mateFailed, bool upstream){
        std::stringstream ss;
        int flag = 4;
        if(paired) flag += 1;
        if(mateFailed) flag += 8;
        if(paired && upstream) flag += 64;
        if(paired && !upstream) flag += 128;
        ss << qname << "\t" << flag << "\t*\t0\t0\t*\t*\t0\t0\t" << sequence << "\t" << qual << std::endl;
        return ss.str();
    }

    std::string printUnpairedAlignment(int i){
        std::stringstream ss;
        alignment_t * a = alignments[i];
        ss << qname << "\t" << a->flag.to_ulong() << "\t"
            << a->rname << "\t" << a->pos << "\t" << 
            a->mapq << "\t" << a->cigar << "\t" << 
            "*" << "\t" << 0 << "\t" << 
            0 << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
            "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" << 
            a->alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" << 
            a->NMtag << std::endl;
        return ss.str();
    }

    //(a->flag.test(8) ? "*" :(a->flag.test(4) ? rQual : qual))
    std::string printPairedAlignments(int i, bool firstMate){
        std::stringstream ss;
        alignment_t * a = alignments[i];
        if(a->paired()){
            //flag has to be changed
            ss << qname << "\t" << (a->flag.to_ulong() |a->mate.flag.to_ulong())  << "\t"
            << a->rname << "\t" << a->pos << "\t" << 
            a->mapq << "\t" << a->cigar << "\t" << 
            (a->rname == a->mate.rnext ? "=" : a->mate.rnext) << "\t" << a->mate.pnext << "\t" << 
            a->mate.tLength << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
            "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" << 
            a->alignmentScore+a->mate.alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" << 
            a->NMtag << std::endl;
            return ss.str();
        }
        else{
            a->setUnPaired(firstMate);
            return printUnpairedAlignment(i);
        }
    }
    
    int alignmentCount(){ return alignments.size(); }
    
};

void unpairedMatch(const sparseSA& sa, read_t & read, const align_opt & alnOptions);

//PAIRED END FUNCTIONS
void pairedMatch(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt);

#endif	/* MAPPER_H */

