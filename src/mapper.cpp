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

#include <set>//pairedmatch3
#include <list>//chain
#include <assert.h>

#include "mapper.h"
#include "kbtree.h"
#include "ksw.h"
#include "dp.h" //eventually remove

using namespace std;

static const int MIN_SEED_LENGTH = 20;
static const int MAX_MAPQ = 60;
static const long MAX_OVERLAP = -10;

//////////////////////////
//struct alignment_t
//////////////////////////

int bandSizeEstimate(int fixedBandSize, int qLen, int refLen, int matchScore, int gapOpen, int gapExtend, int alignmentScore){
    int minBandSize = max(qLen-refLen,refLen-qLen)+3;
    int scoreBandSize = (min(qLen,refLen)*matchScore - alignmentScore + gapOpen)/(-1*gapExtend) + 1;
    int maxGapBandSize = minBandSize-3 + ((qLen>>1)*matchScore+gapOpen)/(-1*gapExtend) + 1;
    int band = min(maxGapBandSize,scoreBandSize);
    band = min(fixedBandSize,band);
    band = max(minBandSize,band);
    return band;
}

string recoverCIGAR3from4(string cigar4){
    stringstream ss(cigar4);
    stringstream cigar3;
    int i;
    char symbol;
    int matchLength = 0;
    while (ss >> i){
        ss.get(symbol);
        if(symbol == 'X' || symbol == '=')
            matchLength += i;
        else{
            //finish previous series of match/mismatch
            if(matchLength > 0){
                cigar3 << matchLength << 'M';
                matchLength = 0;
            }
            cigar3 << i << symbol;
        }
    }
    //finish possible last series of matches
    if(matchLength > 0)
        cigar3 << matchLength << 'M';
    return cigar3.str();
}

//recover cigar string (2 versions) and edit distance
int recoverCIGAR4(int n_cigar, uint32_t *BamCigar, const uint8_t * query, const uint8_t * rseq, stringstream & sNM){
    int i, k, x, y, n_mm = 0, n_gap = 0;
    for (k = 0, x = y = 0; k < n_cigar; ++k) {
        int op  = BamCigar[k]&0xf;
        int len = BamCigar[k]>>4;
        if (op == 0) { // match
            int matchLen = 0;
            int mmLen = 0;
            for(i = 0; i < len; i++){
                if(query[x + i] != rseq[y + i]){
                    ++n_mm; ++mmLen;
                    if(matchLen){
                        sNM << matchLen << "=";
                        matchLen = 0;
                    }
                }
                else{
                    ++matchLen;
                    if(mmLen){
                        sNM << mmLen << "X";
                        mmLen = 0;
                    }
                }
            }
            if(matchLen)
                sNM << matchLen << "=";
            else
                sNM << mmLen << "X";
            x += len; y += len;
        } 
        else{
            sNM << len << "MIDSH"[op];
            if (op == 1) x += len, n_gap += len;
            else if (op == 2) y += len, n_gap += len;
        }
    }
    return n_mm + n_gap;
}

// compute NM
void alignment_t::setFieldsKsw(int n_cigar, uint32_t *BamCigar, const uint8_t * query, const uint8_t * rseq, int beginGap, int endGap){
    int i, k, x, y, n_mm = 0, n_gap = 0;
    stringstream sNM;
    stringstream sCig;
    if(beginGap){
        sCig << beginGap << "S";
        sNM << beginGap << "S";
    }    
    for (k = 0, x = y = 0; k < n_cigar; ++k) {
        int op  = BamCigar[k]&0xf;
        int len = BamCigar[k]>>4;
        sCig << len << "MIDSH"[op];
        if (op == 0) { // match
            int matchLen = 0;
            int mmLen = 0;
            for(i = 0; i < len; i++){
                if(query[x + i] != rseq[y + i]){
                    ++n_mm; ++mmLen;
                    if(matchLen){
                        sNM << matchLen << "=";
                        matchLen = 0;
                    }
                }
                else{
                    ++matchLen;
                    if(mmLen){
                        sNM << mmLen << "X";
                        mmLen = 0;
                    }
                }
            }
            if(matchLen)
                sNM << matchLen << "=";
            else
                sNM << mmLen << "X";
            x += len; y += len;
        } 
        else{
            sNM << len << "MIDSH"[op];
            if (op == 1) x += len, n_gap += len;
            else if (op == 2) y += len, n_gap += len;
        }
    }
    if(endGap){
        sCig << endGap << "S";
        sNM << endGap << "S";
    }  
    editDist = n_mm + n_gap;
    cigar = sCig.str();
    NMtag = sNM.str();
}

void alignment_t::constructCIGAR(const string& S, const uint8_t * readP, 
        int Plen, int fixedBand, const int8_t *mat, int match, int openGap, int extendGap, int print){
    //use chain to recreate succesfull alignment, but now use other method that recovers alignment
    stringstream sNM;
    editDist = 0;
    int costEx = 0; //performanceCounter
    //start with beginGap:
    if(qB){
        if(print >= 4) cerr << "clipping before start alignment: " << qB << endl;
        sNM << qB << "S";
    }
    assert(chain.size()>0);
    int qLen = qE - qB + 1;
    int refLen = (int)(refE - refB + 1L);
    int band;
    uint8_t *target;
    target = (uint8_t *)malloc(refLen);
    for (int i = 0; i < refLen; ++i) // convert to the nt4 encoding
        target[i] = Utils::ORDVALUE[(size_t)S[refB+i]];
    //before first element of chain
    if(chain[0].query > qB){
        qLen = chain[0].query - qB;
        refLen = (int)(chain[0].ref - refB);
        int n_cigar;
        uint32_t * cigar;
        band = bandSizeEstimate(fixedBand, qLen, refLen, match, openGap, extendGap, alignmentScore);
        costEx += refLen*band;
        if(print >=7 ) cerr << "dp called with dimension " << (qLen) << "x" << (refLen) << " and band " << band << endl;
        ksw_global(qLen, readP+qB, refLen, target, 5, mat, -1*openGap, -1*extendGap, band, &n_cigar, &cigar);
        editDist += recoverCIGAR4(n_cigar, cigar, readP+qB, target, sNM);
        if(print >=4 ) cerr << "current cigar " << sNM.str() << endl;
        free(cigar);        
    }
    sNM << chain[0].len << '=';
    //traverse chain
    for(size_t i = 1; i < chain.size(); i++){
        int qDist = chain[i].query - (chain[i-1].query + chain[i-1].len - 1);
        int refDist = (int) (chain[i].ref - (chain[i-1].ref + chain[i-1].len - 1L));
        if (qDist < 0 || refDist < 0) {
            if (qDist >= refDist)
                qDist -= refDist; //qDist insertions
            else
                refDist -= qDist;
        }
        if (qDist <= 0) {//minRefDist deletions //1+|qDist| matches lost of next mem
            sNM << refDist << 'D' << (chain[i].len - 1 + qDist) << '=';
            editDist += refDist;
        } else if (refDist <= 0) {
            sNM << qDist << 'I' << (chain[i].len - 1 + refDist) << '=';
            editDist += qDist;
        } else if (qDist == 1) {
            sNM << (refDist-1) << 'D' << chain[i].len << '=';
            editDist += refDist-1;
        } else if (refDist == 1) {
            sNM << (qDist-1) << 'I' << chain[i].len << '=';
            editDist += qDist-1;
        } else if (refDist == 2 && qDist==2) {
            sNM << 1 << 'X' << chain[i].len << '=';
            editDist += 1;
        } else {//both distances are positive and not equal to (1,1), otherwise seeds would not be maximal!
            qLen = chain[i].query - 1 - (chain[i-1].query + chain[i-1].len) + 1;
            refLen = (int)(chain[i].ref - 1 - (chain[i-1].ref + chain[i-1].len) + 1);
            int n_cigar;
            uint32_t * cigar;
            band = bandSizeEstimate(fixedBand, qLen, refLen, match, openGap, extendGap, alignmentScore);
            costEx += refLen*band;
            if(print >= 7) cerr << "dp called with dimension " << (qLen) << "x" << (refLen) << " and band " << band << endl;
            ksw_global(qLen, readP+chain[i-1].query + chain[i-1].len, refLen, target+(chain[i-1].ref + chain[i-1].len-refB), 5, mat, -1*openGap, -1*extendGap, band, &n_cigar, &cigar);
            editDist += recoverCIGAR4(n_cigar, cigar, readP+chain[i-1].query + chain[i-1].len, target + (chain[i-1].ref + chain[i-1].len-refB), sNM);
            sNM << chain[i].len << '=';
            free(cigar);
        }
        if(print>= 4) cerr << "current cigar " << sNM.str() << endl;
    }
    //after last element of chain
    int lastIndex = chain.size()-1;
    if(chain[lastIndex].query+chain[lastIndex].len-1 < qE){
        qLen = qE - (chain[lastIndex].query + chain[lastIndex].len) + 1;
        refLen = (int)(refE - (chain[lastIndex].ref + chain[lastIndex].len) + 1);
        int n_cigar;
        uint32_t * cigar;
        band = bandSizeEstimate(fixedBand, qLen, refLen, match, openGap, extendGap, alignmentScore);
        if(print>=7) cerr << "dp called with dimension " << (qLen) << "x" << (refLen) << " and band " << band << endl;
        costEx += refLen*band;
        ksw_global(qLen, readP+chain[lastIndex].query + chain[lastIndex].len, refLen, target+(chain[lastIndex].ref + chain[lastIndex].len-refB), 5, mat, -1*openGap, -1*extendGap, band, &n_cigar, &cigar);
        editDist += recoverCIGAR4(n_cigar, cigar, readP+chain[lastIndex].query + chain[lastIndex].len, target+(chain[lastIndex].ref + chain[lastIndex].len-refB), sNM);
        if(print>=4) cerr << "current cigar " << sNM.str() << endl;
        free(cigar);        
    }    
    //end with possible endGap
    if(Plen-1-qE){
        if(print>=4) cerr << "clipping after end alignment: " << (Plen-1-qE) << endl;
        sNM << (Plen-1-qE) << "S";
    }
    //final results:
    free(target);    
    NMtag = sNM.str();
    cigar = recoverCIGAR3from4(NMtag);
    if(print>=2) cerr << "#COSTPOSTPROC " << costEx << endl;
}

void alignment_t::setMate(const alignment_t * o, bool concordant, bool upstream){
    mate.flag.set(0,true);//positions for mate to be important: 0,1,5,6,7
    mate.flag.set(1,true);
    mate.flag.set(5,o->flag.test(4));
    mate.flag.set(6,upstream);
    mate.flag.set(7,!upstream);
    mate.concordant = concordant;
    mate.pnext = o->pos;
    mate.pnextGlob = o->refB;
    mate.tLength = max(refE,o->refE) - min(refB,o->refB) + 1;
    long firstPos = flag.test(4) ? refE : refB;
    long secondPos = o->flag.test(4) ? o->refE : o->refB;
    if(firstPos > secondPos) mate.tLength *= -1;
    mate.rnext = rname == o->rname ? "=" : o->rname;
    if(concordant)
        mate.alignmentScore = o->alignmentScore;
}

//////////////////////////
//struct read_t
//////////////////////////

bool compAlignments(const alignment_t * i,const alignment_t * j){
    return (i->alignmentScore > j->alignmentScore);
}

bool compPairedAlignments(const alignment_t * i,const alignment_t * j){
    return (i->alignmentScore+i->mate.alignmentScore > j->alignmentScore+j->mate.alignmentScore);
}

void fullDP(const string & S, const uint8_t * readP, int Plen, alignment_t * alignment, int fixedBand, 
        const int8_t *mat, int match, int openGap, int extendGap, int print){
    int qLen = alignment->qE - alignment->qB + 1;
    int refLen = (int)(alignment->refE - alignment->refB + 1L);
    int band = bandSizeEstimate(fixedBand, qLen, refLen, match, openGap, extendGap, alignment->alignmentScore);
    //calculate CIGAR STRING
    uint8_t *target;
    target = (uint8_t *)malloc(refLen);
    for (int i = 0; i < refLen; ++i) // convert to the nt4 encoding
        target[i] = Utils::ORDVALUE[(size_t)S[alignment->refB+i]];
    int n_cigar;
    uint32_t * cigar;
    if(print>=2) cerr << "#COSTPOSTPROC " << refLen*band << endl;
    alignment->alignmentScore = ksw_global(qLen, readP+alignment->qB, refLen, target, 5, mat, -1*openGap, -1*extendGap, band, &n_cigar, &cigar);
    alignment->setFieldsKsw(n_cigar, cigar, readP+alignment->qB, target, alignment->qB, Plen-1-alignment->qE);
    free(target);
    free(cigar);
}

void read_t::postprocess(const sparseSA& sa, const align_opt & alnOptions, bool paired){
    //mapq independent from match/mismatch
    if(alnOptions.print >= 3) cerr << "postprocessing read ... " << endl;
    if(!alignments.empty()){
        //FIRST SORT
        if(!paired)
            sort(alignments.begin(),alignments.end(), compAlignments);
        else
            sort(alignments.begin(),alignments.end(), compPairedAlignments);
        int chainCount = 0;
        int maxDP = min((int)alignments.size(),alnOptions.alignmentCount);
        if(alnOptions.fullDP)
            maxDP = min((int)alignments.size(),max(4,alnOptions.alignmentCount));//TODO set parameter
        for(int j = 0; j < maxDP; j++){
            alignments[j]->setLocalPos(sa);
            if(!alignments[j]->chain.empty() && !alnOptions.fullDP){
                chainCount++;
                alignments[j]->constructCIGAR(sa.S, alignments[j]->flag.test(4) ? rcReadP : readP, 
                    sequence.length(), alnOptions.fixedBandSize, alnOptions.mat, 
                    alnOptions.match, alnOptions.openGap, alnOptions.extendGap, alnOptions.print);
            }
            else{
                //full dp version
                fullDP(sa.S, alignments[j]->flag.test(4) ? rcReadP : readP, sequence.size(), 
             alignments[j], alnOptions.fixedBandSize, alnOptions.mat, 
             alnOptions.match, alnOptions.openGap, alnOptions.extendGap, alnOptions.print);
            }
            if(alnOptions.print >= 3) cerr << "edit dist: " << alignments[j]->editDist << " and cigar " << alignments[j]->cigar << endl;
        }
        if(alnOptions.print >= 2){ 
            cerr << "#POSTPROCCHAIN " << chainCount << endl;
            cerr << "#POSTPROCFULL " << (min((int)alignments.size(),alnOptions.alignmentCount)-chainCount) << endl;
        }
        //SECOND SORT
        if(!paired)
            sort(alignments.begin(),alignments.end(), compAlignments);
        else
            sort(alignments.begin(),alignments.end(), compPairedAlignments);
        int maxScore = alignments[0]->alignmentScore + alignments[0]->mate.alignmentScore;
        int secBestScore = alignments.size() > 1 ? alignments[1]->alignmentScore + alignments[1]->mate.alignmentScore : alnOptions.mismatch*sequence.length();
        int mapq = (secBestScore == (int)(alnOptions.mismatch*sequence.length())) ? MAX_MAPQ :
            (int)(250*((double)(maxScore-secBestScore) / maxScore) + 0.499);
        mapq = maxScore == secBestScore ? 0 : mapq;
        mapq = min(max(0,mapq),MAX_MAPQ);
        if(alnOptions.print >= 2){ 
            cerr << "#MAPQ " << mapq << endl;
            cerr << "#POSTPROC " << min((int)alignments.size(),alnOptions.alignmentCount) << endl;
        }
        for(int j = 0; j < min((int)alignments.size(),alnOptions.alignmentCount); j++){
            if(alnOptions.print >= 3) cerr << "constructing cigar for alignment on " << alignments[j]->rname << " and pos " << alignments[j]->pos << endl;
            if(alignments[j]->alignmentScore + alignments[j]->mate.alignmentScore == maxScore)
                alignments[j]->mapq = mapq;
            else
                alignments[j]->flag.set(8,true);
        }
        
    }
}

/////////////////////
//DEBUG FUNCTIONS
////////////////////
#ifndef NDEBUG

string createCigar(vector<char> & chars, vector<int> & lengths){
    stringstream sNM;
    assert(chars.size()==lengths.size());
    assert(chars.size()>0);
    size_t i = 0;
    while(i < chars.size()){
        if(chars[i]=='=' || chars[i]=='X'){
            int tempLength = lengths[i];
            sNM << lengths[i] << chars[i];
            while(i < chars.size()-1 && (chars[i+1]=='=' || chars[i+1]=='X')){
                i++;
                tempLength += lengths[i];
                sNM << lengths[i] << chars[i];
            }
        }
        else
            sNM << lengths[i] << chars[i];
        i++;
    }
    return sNM.str();
}

static string previousWgsimEntry(string & line, char delimeter, int& endPos){
    int tabPos = line.rfind(delimeter,endPos);
    string substring = line.substr(tabPos+1, endPos-tabPos);
    endPos = tabPos-1;
    return substring;
}

dp_output dpInRegion(const string & S, const string & P, long & refB, long refE, int qB, int qE){
    dp_output output;
    boundaries grenzen(refB,refE,qB,qE);
    dp_type types;
    types.freeRefB = true;
    types.freeQueryB = false;
    types.freeQueryE = false;
    types.freeRefE = true;
    long INIT_DP_DIMENSION = 2048;
    dynProg * dp_;
    dp_scores scores = dp_scores(0, -2, 0, -2);
    dp_ = new dynProg(INIT_DP_DIMENSION, false, scores);
    dp_->scores.updateScoreMatrixDna();
    dp_->dpBasic( S, P, grenzen, types, ERRORSTRING, output);
    delete dp_;
    refB = grenzen.refB;
    return output;
}

void printSimRegion(const sparseSA& sa, const read_t & read, const align_opt & alnOptions, int editDist){
    int Plength = read.sequence.length();
    //wgsim read name to dp simulated position
    string qname = read.qname;
    if(qname.length() > 2 && qname[qname.length()-2]=='/')
        qname.erase(qname.length()-2);
    int namePos=qname.size()-1;
    previousWgsimEntry(qname, '_', namePos);//counter
    int revEdit = atoi(previousWgsimEntry(qname, ':', namePos).c_str());
    revEdit += atoi(previousWgsimEntry(qname, ':', namePos).c_str());
    revEdit += atoi(previousWgsimEntry(qname, '_', namePos).c_str());
    int fwdEdit = atoi(previousWgsimEntry(qname, ':', namePos).c_str());
    fwdEdit += atoi(previousWgsimEntry(qname, ':', namePos).c_str());
    fwdEdit += atoi(previousWgsimEntry(qname, '_', namePos).c_str());
    long revPos = atoi(previousWgsimEntry(qname, '_', namePos).c_str());
    long fwdPos = atoi(previousWgsimEntry(qname, '_', namePos).c_str());
    string realchrom = qname.substr(0, namePos+1);
    size_t chrIndex = 0;
    while(chrIndex < sa.descr.size() && sa.descr[chrIndex].compare(realchrom)!=0)
        chrIndex++;
    if(chrIndex < sa.descr.size()){
        //calculate alignment on simulated region
        long chrStart = sa.startpos[chrIndex];
        long chrEnd = (chrIndex+1 < sa.descr.size() ? sa.startpos[chrIndex+1]-1 : sa.N);
        long alignmentBoundLeft = max(min(fwdPos,(long)(revPos-Plength))+chrStart-3*editDist,chrStart);
        long alignmentBoundRight = min(max(revPos,(long)(fwdPos+Plength))+chrStart+3*editDist,chrEnd);
        //perform alignment twice: 1 time for both strands
        long refBfw = alignmentBoundLeft;
        long refBrc = alignmentBoundLeft;
        dp_output outputFw = dpInRegion(sa.S, read.sequence, refBfw, alignmentBoundRight, 0, Plength-1);
        dp_output outputRc = dpInRegion(sa.S, read.rcSequence, refBrc, alignmentBoundRight, 0, Plength-1);
        alignment_t * alignmentFull = new alignment_t();
        string strand = "";
        long bestEdit = 0L;
        string cigarstring = "";
        if(outputFw.editDist < outputRc.editDist){
            strand = "f";
            alignmentFull->refB = refBfw+1L;
            bestEdit = outputFw.editDist;
            cigarstring = createCigar(outputFw.cigarChars, outputFw.cigarLengths);
        }
        else{
            strand = "r";
            alignmentFull->refB = refBrc+1L;
            bestEdit = outputRc.editDist;
            cigarstring = createCigar(outputRc.cigarChars, outputRc.cigarLengths);
        }
        alignmentFull->setLocalPos(sa);
        fprintf(alnOptions.debugFile,"#SIM %s %ld %s %ld %ld %s\n", strand.c_str(), alignmentFull->refB, realchrom.c_str(), alignmentFull->pos, bestEdit, cigarstring.c_str());
        delete alignmentFull;
    }
    else{//bad code or different reference
        cerr << "problem in debug code: simulated chromosome is not found in index. Either non-Wgsim reads are used or a bug has been found" << endl;
        exit(1);
    }
}

void printLisData(const sparseSA& sa, const lis_t & cluster, const align_opt & alnOptions, int index, bool extended, bool unique, bool trial, int Plength, int editDist){
    const vector<match_t> & matchVector = cluster.matches;
    long windowBegin, windowEnd, windowOffset, windowReach,seqIndex; 
    string chromosome;
    windowBegin = matchVector[0].ref - matchVector[0].query - editDist;
    for(size_t i= 0; i < matchVector.size(); i++)
        windowBegin = min(windowBegin,matchVector[0].ref-matchVector[0].query-editDist);
    windowEnd = matchVector[0].ref-matchVector[0].query+Plength+editDist;
    for(size_t i= 0; i < matchVector.size(); i++)
        windowEnd = max(windowEnd,(long)(matchVector[0].ref-matchVector[0].query+Plength+editDist));
    sa.from_set(windowBegin,seqIndex,windowOffset);
    sa.from_set(windowEnd,seqIndex,windowReach);
    chromosome = sa.descr[seqIndex];
    fprintf(alnOptions.debugFile,"$%d %d %s %lu %d %d %d %ld %ld %s %ld %ld\n", index, 
            extended ? 1 : 0, 
            cluster.fw ? "f" : "r", 
            matchVector.size(), cluster.score, unique ? 1 : 0, trial ? 1 : 0, 
            windowBegin, windowEnd, chromosome.c_str(), windowOffset, windowReach);
    for(size_t i = 0; i < matchVector.size(); i++)
        fprintf(alnOptions.debugFile,"%d %ld %d\n", matchVector[i].query, matchVector[i].ref, matchVector[i].len);
}

#endif

////////////////////////////
//CALCULATE SEED FUNCTIONS
///////////////////////////

inline void executeSeedType(const sparseSA& sa, const string& P, int min_len, mum_t memType, int maxBranchWidth, int sparseSkip, vector<match_t>& matches, int maxSmemMemStart){
    if(memType == SMAM)
        sa.SMAM(P, matches, min_len, maxBranchWidth, sparseSkip);
    else if(memType == SMEM)
        sa.SMEM(P, matches, min_len, maxBranchWidth, sparseSkip, maxSmemMemStart);
    else if(memType == MEM)
        sa.MEM(P, matches, min_len, maxBranchWidth, sparseSkip);
    else if(memType == SMEMP)
        sa.SMEMP(P, matches, min_len, maxBranchWidth, sparseSkip, maxSmemMemStart);
}

void calculateSeedsBothStrands(const sparseSA& sa, const string& P, const string& Prc, int min_len, const align_opt & alnOptions, vector<match_t>& matches, vector<match_t>& matchesRC){
    int maxBranchWidth = alnOptions.maxSeedCandidates;
    mum_t memType = alnOptions.memType;
    int sparseSkip = 1;
    if (sa.K >= 4) sparseSkip = max((min_len - 10) / (int)sa.K,1);
    else sparseSkip = max((min_len - 12) / (int)sa.K,1);
    if(alnOptions.noSparseQ) sparseSkip = 1;
    if(alnOptions.print >= 6) cerr << "min_len " << min_len << " and sparseSkip " << sparseSkip << endl;
    if(!alnOptions.noFW)
        executeSeedType(sa, P, min_len, memType, maxBranchWidth, sparseSkip, matches, alnOptions.maxSmemMemStart);
    if(!alnOptions.noRC)
        executeSeedType(sa, Prc, min_len, memType, maxBranchWidth, sparseSkip, matchesRC, alnOptions.maxSmemMemStart);
    if(alnOptions.rescue && matches.empty() && matchesRC.empty()){
        if(alnOptions.print >= 2) cerr << "#MEMRESCUE" << endl;
        int minLength = max(MIN_SEED_LENGTH, (int)sa.K+10);//requires some difference with sparseness
        if(!alnOptions.noFW)
            executeSeedType(sa, P, (minLength+min_len)/2, memType, 100000, 1, matches, alnOptions.maxSmemMemStart);
        if(!alnOptions.noRC)
            executeSeedType(sa, Prc, (minLength+min_len)/2, memType, 100000, 1, matchesRC, alnOptions.maxSmemMemStart);
        if(matches.empty() && matchesRC.empty()){
            if(!alnOptions.noFW)
                executeSeedType(sa, P, minLength, memType, 100000, 1, matches, alnOptions.maxSmemMemStart);
            if(!alnOptions.noRC)
                executeSeedType(sa, Prc, minLength, memType, 100000, 1, matchesRC, alnOptions.maxSmemMemStart);
        }
    }
    if(alnOptions.print >= 2){ 
        cerr << "#MEMF " << matches.size() << endl;
        cerr << "#MEMR " << matchesRC.size() << endl;
        if(alnOptions.print >= 6){
            cerr << "seeds fw: " << endl;
            for(size_t i = 0; i < matches.size(); i++)
                cerr << "q,r,l " << matches[i].query << "," << matches[i].ref << "," << matches[i].len << endl;
            cerr << "seeds rev: " << endl;
            for(size_t i = 0; i < matchesRC.size(); i++)
                cerr << "q,r,l " << matchesRC[i].query << "," << matchesRC[i].ref << "," << matchesRC[i].len << endl;
        }
    }
}

bool compMatches(const match_t & i, const match_t & j){
    return (i.ref < j.ref || (i.ref == j.ref && i.query < j.query));
}

bool compMatchesSMEM(const match_t & i,const match_t & j){
    return (i.query < j.query || (i.query == j.query && i.len > j.len));
}

void identifyAndFilter(vector<match_t>& matches, bool fw, int qLength, int editDist, vector<lis_t>& lisIntervals, const align_opt & alnOptions){
    if(!matches.empty()){
        size_t emptySize = lisIntervals.size();
        //FIND SMEMS AND REGIONS
        size_t maxSmemOcc = alnOptions.maxSeedCandidates;
        long minMemLength = alnOptions.minMemRegLength;
        size_t maxPossOccIncrease = alnOptions.maxRegionMems;
        //sort according to query offset and len
        sort(matches.begin(),matches.end(), compMatchesSMEM);
        ///////////////////////////////////////////////////
        //TRIAGE: smems and mems to fill candidate regions
        ///////////////////////////////////////////////////
        vector<match_t> smems;
        vector<match_t> reserves;
        smems.reserve(matches.size());
        reserves.reserve(matches.size());
        long previousStart = matches[0].query;
        long previousEnd = matches[0].query + matches[0].len;
        long previousREnd = -1L;
        size_t rangeStart = 0;
        size_t currentSmemRange = 0;
        size_t index = 0;
        ///////////////////////////
        //IDENTIFY SMEMS
        //////////////////////////
        if(alnOptions.print >= 5) cerr << "number of MEMs before calculating regions: " << matches.size() << endl;
        while(index < matches.size()){
            rangeStart = index;
            previousStart = matches[index].query;
            previousEnd = matches[index].query + matches[index].len;
            index++;
            //same mem/smem
            while(index < matches.size() && matches[index].query == previousStart && matches[index].query+matches[index].len == previousEnd)
                index++;
            //region of same query position seeds: insert this region in triage
            assert(index>=rangeStart);
            if(previousEnd > previousREnd)
                currentSmemRange = index-rangeStart;
            else
                currentSmemRange += index-rangeStart;
            if ((previousEnd > previousREnd && currentSmemRange <= maxSmemOcc) || 
                    (currentSmemRange <= maxPossOccIncrease 
                    && previousEnd-previousStart >= minMemLength)){
                        for(size_t j = rangeStart; j < index; j++)
                            smems.push_back(matches[j]);
            }
            else {
                for(size_t j = rangeStart; j < index; j++)
                    reserves.push_back(matches[j]);
            }
            //new values of fields
            previousREnd = previousEnd;
        }
        if(alnOptions.print >= 5){ 
            cerr << "number of SMEMs: " << smems.size() << endl;
            cerr << "number of reserves: " << reserves.size() << endl;
        }
        ///////////////////////////
        //CREATE CANDIDATE REGIONS
        //////////////////////////
        sort(smems.begin(),smems.end(), compMatches);
        int begin, end, longest;
        begin = end = longest = 0;
        long refEnd, refBegin;
        while(begin < (int) smems.size()){
            //left most --> define right most boundary
            refEnd = smems[begin].ref - smems[begin].query + qLength + editDist;
            longest = begin;
            end = begin+1;
            //find right most boundary
            while(end < (int) smems.size() && smems[end].ref + smems[end].len <= refEnd){
                if(smems[longest].len < smems[end].len){
                    longest = end;
                    refEnd = smems[longest].ref - smems[longest].query + qLength + editDist;
                }
                end++;
            }
            //first to left
            refBegin = smems[longest].ref - smems[longest].query - editDist;
            //begin = longest-1;
            begin--;
            int beginFix = begin;//bugfix: merge region contains some seeds twice
            while(begin >= 0 && smems[begin].ref >= refBegin)
                begin--;
            //add interval to list and calculate coverage
            if(emptySize==lisIntervals.size() || lisIntervals.back().refEnd < refBegin){
                //no overlap, calc coverage and push
                lisIntervals.push_back(lis_t(fw));
                lisIntervals.back().addSeeds(smems,begin+1,end-1);
                lisIntervals.back().refBegin = refBegin;
                lisIntervals.back().refEnd = refEnd;
                if(alnOptions.print >= 5) cerr << "region found with " << (end-begin-1) << " seeds" << endl;
            }
            else{//overlap with previous: merge two regions
                lisIntervals.back().addSeeds(smems,beginFix+1,end-1);
                lisIntervals.back().refBegin = min(lisIntervals.back().refBegin,refBegin);
                lisIntervals.back().refEnd = min(lisIntervals.back().refEnd,refEnd);
                if(alnOptions.print >= 5) cerr << "previous region extended with " << (end-beginFix-1) << " seeds" << endl;
            }
            //start new interval
            begin = end;
        }
        if(alnOptions.print >= 5) cerr << (lisIntervals.size() - emptySize) << " regions found" << endl;
        ///////////////////////////
        //FILL CANDIDATE REGIONS
        //////////////////////////
        //region fill: merge reserves in regions fill both matches and lisIntervals
        int reservesAdded = 0;
        if(emptySize < lisIntervals.size()){
            sort(reserves.begin(),reserves.end(), compMatches);
            size_t curReserve = 0;
            size_t curInterval = emptySize;
            while(curReserve < reserves.size() && curInterval < lisIntervals.size()){
                if(reserves[curReserve].ref < lisIntervals[curInterval].refBegin)
                    curReserve++;
                else if(lisIntervals[curInterval].refEnd < reserves[curReserve].ref)
                    curInterval++;
                else{
                    lisIntervals[curInterval].addSeed(reserves[curReserve]);
                    reservesAdded++;
                    curReserve++;
                }
            }
        }
        if(alnOptions.print >= 5) cerr << reservesAdded << " extra seeds added to region" << endl;
    }
}

int coveredBases(lis_t& interval){
    //sort this candidate region by query position
    interval.sortByReadPos();
    vector<match_t> & seeds = interval.matches;
    int trueCover = seeds[0].len;
    int prevEnd = seeds[0].query + seeds[0].len-1;
    for(size_t i = 1; i < seeds.size(); i++ ){
        int start = max(seeds[i].query,prevEnd+1);
        int stop = max(prevEnd,seeds[i].query + seeds[i].len -1);
        prevEnd = stop;
        trueCover += max(0,stop-start+1);
    }
    return trueCover;
}

int nextChain(const vector<match_t>& matches, list<size_t>& chain, size_t anchorPos, double errorPercent, vector<bool>& used, size_t qLength, int & coverage){
    int curEditDist = 0;
    coverage = 0;
    coverage += (int)matches[anchorPos].len;
    chain.clear();
    chain.push_back(anchorPos);
    //boundaries
    const match_t & firstSeed = matches[anchorPos];
    long refstrLB = firstSeed.ref;
    long queryLB = firstSeed.query;
    long refstrRB = refstrLB + firstSeed.len -1L;
    long queryRB = queryLB + firstSeed.len-1L;
    
    //extend chain to the right
    int memsIndex = anchorPos+1;
    while (memsIndex < (int) matches.size()) {
            const match_t & match = matches[memsIndex];
            //if maxoverlap read, reference, skew
            long qDist = match.query - (queryRB + 1);
            long refDist = match.ref - (refstrRB + 1);
            long qDistLong = match.query+match.len - queryLB;
            long refDistLong = match.ref+match.len - refstrLB;
            if(qDist>=MAX_OVERLAP && refDist>=MAX_OVERLAP && abs(qDist-refDist)<=(long)(errorPercent*max(qDistLong,refDistLong))){
                chain.push_back(memsIndex);
                curEditDist += max(1L,abs(qDist-refDist));
                refstrLB = match.ref;
                refstrRB = match.ref + match.len - 1L;
                queryLB = match.query;
                queryRB = match.query + match.len - 1L;
                used[memsIndex] = true;
                coverage += min(qDist,refDist) < 0 ? (int)(matches[memsIndex].len+min(qDist,refDist)) : (int)matches[memsIndex].len;
            }
        memsIndex++;
    }
    curEditDist += (int)(queryRB + 1) < (int) qLength ? 1 : 0;
    //extend chain to the left
    memsIndex = anchorPos-1;
    refstrLB = firstSeed.ref;
    queryLB = firstSeed.query;
    refstrRB = refstrLB + firstSeed.len -1L;
    queryRB = queryLB + firstSeed.len-1L;
    while (memsIndex >= 0) {
        const match_t & match = matches[memsIndex];
        long qDist = queryLB - (match.query+match.len);
        long refDist = refstrLB - (match.ref+match.len);
        long qDistLong = queryRB - (match.query)+1;
        long refDistLong = refstrRB - (match.ref)+1;
        if(qDist>=MAX_OVERLAP && refDist>=MAX_OVERLAP && abs(qDist-refDist)<=(long)(errorPercent*max(qDistLong,refDistLong))){//errorPercent*(qLength/4)
            chain.push_front(memsIndex);
            curEditDist += max(1L,abs(qDist-refDist));
            refstrLB = match.ref;
            refstrRB = match.ref + match.len -1L;
            queryLB = match.query;
            queryRB = match.query + match.len-1L;
            used[memsIndex] = true;
            coverage += min(qDist,refDist) < 0 ? (int)(matches[memsIndex].len+min(qDist,refDist)) : (int)matches[memsIndex].len;
        }
        memsIndex--;
    }
    curEditDist += queryLB > 0 ? 1 : 0;
    return curEditDist;
}

inline void leftExtension(const string& S, const string& P, int & tmpQB, long & tmpRB, 
        long chrStart, const align_opt & alnOptions, int & curScore, int bandSize, int editDist){
    if(tmpQB>0 && tmpRB > chrStart){
        int zdrop = 100;
        if(alnOptions.print >= 7) cerr << "dp required before first seed" << endl;
        int qle, tle, gtle, gscore, max_off;
        uint8_t *target;
        long alignmentBoundLeft = max(tmpRB-tmpQB-(editDist),chrStart);
        int tLen = tmpRB-alignmentBoundLeft;
        target = (uint8_t *)malloc(tLen);
        for (int j = 0; j < tLen; ++j) // convert to the nt4 encoding
            target[j] = Utils::ORDVALUE[(size_t)S[tmpRB-j-1]];
        uint8_t *query;
        query = (uint8_t *)malloc(tmpQB);
        for (int j = 0; j < tmpQB; ++j)
            query[j] = Utils::ORDVALUE[(size_t)P[tmpQB-j-1]];
        if(alnOptions.print >= 7) cerr << "dp called with dimension " << (tmpQB) << "x" << (tLen) << " and band " << min(bandSize,tmpQB) << endl;
        if(alnOptions.print >= 2) cerr << "#DP " << (min(bandSize,tmpQB)*tLen) << endl;
        curScore = ksw_extend(tmpQB, query, tLen, target, 5, alnOptions.mat, 
            -1*alnOptions.openGap, -1*alnOptions.extendGap, 
            min(bandSize,tmpQB), zdrop, curScore, &qle, &tle, &gtle, &gscore, &max_off);
        if(alnOptions.print >= 7) cerr << "curScore is " << curScore << endl;
        if(alnOptions.global || gscore >= curScore + (alnOptions.mismatch-alnOptions.match)*3){
            curScore = gscore;
            tmpQB = 0;
            tmpRB = tmpRB-gtle;
            if(alnOptions.print >= 7) cerr << "global score is " << gscore << " with refbegin " << tmpRB << endl;
        }
        else{
            tmpQB = tmpQB-qle;
            tmpRB = tmpRB-tle;
            if(alnOptions.print >= 7) cerr << "local begin pos is " << tmpQB << " with refBegin " << tmpRB << endl;
        }
        free(query);
        free(target);
    }
    else if(tmpQB>0 && tmpRB == chrStart && alnOptions.global){
        if(alnOptions.print >= 7) cerr << "begin of refseq: " << tmpQB << " insert" << endl;
        curScore += alnOptions.openGap + alnOptions.extendGap*tmpQB;
    }
}

inline void rightExtension(const string& S, const string& P, int & tmpQE, long & tmpRE, 
        long chrEnd, const align_opt & alnOptions, int & curScore, int bandSize, int editDist){
    if(tmpQE < (int)(P.length()-1) && tmpRE < chrEnd){
        int zdrop = 100;
        if(alnOptions.print >= 7) cerr << "end DP required" << endl;
        int qle, tle, gtle, gscore, max_off;
        uint8_t *target;
        long alignmentBoundRight = min((long)(tmpRE+(P.length()-1-tmpQE)+editDist),chrEnd);
        int tLen = alignmentBoundRight-tmpRE;
        target = (uint8_t *)malloc(tLen);
        for (int j = 0; j < tLen; ++j) // convert to the nt4 encoding
            target[j] = Utils::ORDVALUE[(size_t)S[tmpRE+1+j]];
        uint8_t *query;
        int qLen = P.length()-1-tmpQE;
        query = (uint8_t *)malloc(qLen);
        for (int j = 0; j < qLen; ++j)
            query[j] = Utils::ORDVALUE[(size_t)P[tmpQE+1+j]];
        if(alnOptions.print >= 7) cerr << "dp called with dimension " << (qLen) << "x" << (tLen) << " and band " << min(bandSize,qLen) << endl;
        if(alnOptions.print >= 2) cerr << "#DP " << (min(bandSize,qLen)*tLen) << endl;
        curScore = ksw_extend(qLen, query, tLen, target, 5, alnOptions.mat, 
                -1*alnOptions.openGap, -1*alnOptions.extendGap, 
                min(bandSize,qLen), zdrop, curScore, &qle, &tle, &gtle, &gscore, &max_off);
        if(alnOptions.print >= 7) cerr << "curScore is " << curScore << endl;
        if(alnOptions.global || gscore >= curScore + (alnOptions.mismatch-alnOptions.match)*3){
            curScore = gscore;
            tmpQE = P.length()-1;
            tmpRE = tmpRE+gtle;
            if(alnOptions.print >= 7) cerr << "global score is " << gscore << " with refEnd " << tmpRE << endl;
        }
        else{
            tmpQE = tmpQE+qle;
            tmpRE = tmpRE+tle;
            if(alnOptions.print >= 7) cerr << "local end pos is " << tmpQE << " with refEnd " << tmpRE << endl;
        }
        free(query);
        free(target);
    }
    else if(tmpQE < (int)(P.length()-1) && tmpRE == chrEnd && alnOptions.global){
        if(alnOptions.print >= 7) cerr << "end of refseq: " << P.length()-1-tmpQE << " insert" << endl;
        curScore += alnOptions.openGap + alnOptions.extendGap*((int)(P.length()-1-tmpQE));
    }
}

int extendAlignmentChain(const string& S, const string& P, 
        const vector<match_t>& matches, int editDist, int minScore, int minCoverage, 
        const align_opt & alnOptions, long chrStart, long chrEnd, read_t & read, bool fw){
    int maxAlnScore = min(-1, minScore-10);
    if(alnOptions.print >= 5){
        cerr << "min Score is " << minScore << endl;
        cerr << "seeds are: " << endl;
        for(size_t i = 0; i < matches.size(); i++)
            cerr << matches[i].query << "," << matches[i].ref << "," << matches[i].len << endl;
    }
    list<size_t> chain;
    //seeds in order of 
    vector<long> seedOrder;
    seedOrder.reserve(matches.size());
    vector<bool> used;
    used.assign(matches.size(),false);
    for(size_t i = 0; i < matches.size(); i++)
        seedOrder.push_back((long)matches[i].len<<32 | i);
    sort(seedOrder.begin(), seedOrder.end());
    //for all seeds:
    vector<boundaries> alnBounds;
    vector<int> alnScores;
    vector<long> alnBegins;
    for(int i = seedOrder.size()-1; i >=0; i--){
        size_t anchorPos = (uint32_t)seedOrder[i];
        match_t anchor = matches[anchorPos];
        if(alnOptions.print >= 5) cerr << "anchor " << anchor.query << "," << anchor.ref << "," << anchor.len << endl;
        //check if current anchor is not already in alignment
        bool extend = !used[anchorPos];
        size_t j = 0;
        while(j < alnBounds.size() && extend){
            bool notContained = anchor.query < alnBounds[j].queryB || anchor.query+anchor.len >= alnBounds[j].queryE || anchor.ref < alnBounds[j].refB || anchor.ref+anchor.len >= alnBounds[j].refE;
            if(!notContained){
                int qDist = anchor.query - alnBounds[j].queryB;
                int refDist = anchor.ref - alnBounds[j].refB;
                int max_gapL = min(alnOptions.fixedBandSize, (int)((min(qDist,refDist)*alnOptions.match + alnOptions.openGap)/alnOptions.extendGap + 1.));
                extend = max(qDist,refDist) >= max_gapL;
                if(extend){
                    qDist = alnBounds[j].queryE - (anchor.query+anchor.len);
                    refDist = alnBounds[j].refE - (anchor.ref+anchor.len);
                    max_gapL = min(alnOptions.fixedBandSize, (int)((min(qDist,refDist)*alnOptions.match + alnOptions.openGap)/alnOptions.extendGap + 1.));
                    extend = max(qDist,refDist) >= max_gapL;
                }
            }
            j++;
        }
        if(extend){
            //CREATE CHAIN AROUND THIS ANCHOR
            int curScore = 0;
            int minEdit = nextChain(matches, chain, anchorPos, alnOptions.errorPercent, used, P.length(), curScore);
            int chainScore = curScore;//for debug purpose
            if(alnOptions.print >= 2){
                cerr << "#CHAIN " << chain.size() << endl;
                cerr << "#CHAINSCORE " << curScore << endl;
                cerr << "#CHAINSKEW " << minEdit << endl;
                cerr << "#CHAINMINSCORE " << minCoverage << endl;
                if(alnOptions.print >= 5)
                        cerr << "extend chain of size " << chain.size() << " elements and minEdit,curScore =  " << minEdit << "," << curScore << endl;
            }
            //ITERATE THROUGH CHAIN
            if(minEdit <= editDist && (chain.size() > 1 || anchor.len >= minCoverage )){
                std::list<size_t>::const_iterator iterator = chain.begin(), itEnd = chain.end();
                int bandSize = min(alnOptions.fixedBandSize, (int)(editDist-minEdit+1));
                anchor = matches[*iterator];
                iterator++;
                int tmpQB = anchor.query;
                int tmpQE = anchor.query+anchor.len-1;
                long tmpRB = anchor.ref;
                long tmpRE = anchor.ref + anchor.len-1L;
                if(alnOptions.print >= 7) cerr << "first seed (q/r/l): " << tmpQB << "," << tmpRB << "," << anchor.len << endl;
                //INSIDE CHAIN
                while(iterator != itEnd){
                    match_t seed = matches[*iterator];
                    if(alnOptions.print >= 6) cerr << "current seed (q/r/l): " << seed.query << "," << seed.ref << "," << seed.len << endl;
                    long qDist = seed.query - tmpQE;
                    long refDist = seed.ref - tmpRE;
                    if (qDist < 0 || refDist < 0) {
                        if (qDist >= refDist)
                            qDist -= refDist; //qDist insertions
                        else
                            refDist -= qDist;
                    }
                    if (qDist <= 0) {//minRefDist deletions
                        //1+|qDist| matches lost of next mem
                        curScore += alnOptions.openGap + alnOptions.extendGap*refDist;
                        if(alnOptions.print >= 7) cerr << "added deletion of size " << refDist << endl;
                    } else if (refDist <= 0) {
                        curScore += alnOptions.openGap + alnOptions.extendGap*qDist;
                        if(alnOptions.print >= 7) cerr << "added insertion of size " << qDist << endl;
                    } else if (qDist == 1) {
                        curScore += alnOptions.openGap + alnOptions.extendGap*(refDist-1);
                        if(alnOptions.print >= 7) cerr << "added deletion of size " << refDist-1 << endl;
                    } else if (refDist == 1) {
                        curScore += alnOptions.openGap + alnOptions.extendGap*(qDist-1);
                        if(alnOptions.print >= 7) cerr << "added insertion of size " << qDist-1 << endl;
                    } else if (refDist == 2 && qDist==2) {
                        curScore += alnOptions.mismatch;
                        if(alnOptions.print >= 7) cerr << "added mutation" << endl;
                    } else {//both distances are positive and not equal to (1,1), otherwise seeds would not be maximal!
                        uint8_t *target;
                        int tLen = seed.ref - tmpRE - 1;
                        target = (uint8_t *)malloc(tLen);
                        for (int j = 0; j < tLen; ++j) // convert to the nt4 encoding
                            target[j] = Utils::ORDVALUE[(size_t)S[tmpRE+1+j]];
                        uint8_t *query;
                        int qLen = seed.query - tmpQE - 1;
                        query = (uint8_t *)malloc(qLen);
                        for (int j = 0; j < qLen; ++j)
                            query[j] = Utils::ORDVALUE[(size_t)P[tmpQE+1+j]];
                        int innerBand = min(alnOptions.fixedBandSize, (int)(editDist-minEdit+abs(qLen-tLen)+3));
                        innerBand = max(innerBand,abs(qLen-tLen)+(int)(abs(qLen-tLen)*alnOptions.errorPercent));
                        if(alnOptions.print >= 2) cerr << "#DP " << (tLen*innerBand) << endl;
                        if(alnOptions.print >= 7) cerr << "dp called with dimension " << (seed.query-tmpQE+1) << "x" << (seed.ref-tmpRE+1) << " and band " << innerBand << endl;
                        curScore += ksw_nw(qLen, query, tLen, target, 5, alnOptions.mat, 
                                -1*alnOptions.openGap, -1*alnOptions.extendGap, innerBand);
                        free(query);
                        free(target);
                    }
                    if(alnOptions.print >= 6) cerr << "curScore is " << curScore << endl;
                    //update bounds
                    tmpRE = seed.ref + seed.len - 1L;
                    tmpQE = seed.query + seed.len - 1;            
                    ++iterator;
                }
                //lef-right extensions
                if(tmpQB <= (int)(P.length()-1)-tmpQE){
                    leftExtension(S, P, tmpQB, tmpRB, chrStart, alnOptions, curScore, bandSize, editDist-minEdit+1);
                    rightExtension(S, P, tmpQE, tmpRE, chrEnd, alnOptions, curScore, bandSize, editDist);
                }
                else{
                    rightExtension(S, P, tmpQE, tmpRE, chrEnd, alnOptions, curScore, bandSize, editDist);
                    leftExtension(S, P, tmpQB, tmpRB, chrStart, alnOptions, curScore, bandSize, editDist-minEdit+1);
                }
                if(alnOptions.print >= 5) cerr << "current alignments startpos would be " << tmpRB << endl;
                maxAlnScore = max(maxAlnScore, curScore);
                if(curScore >= minScore){
                    if(alnOptions.print >= 2) cerr << "#SUCCES " << curScore << endl;
                    bool equal = false;
                    for(size_t i = 0; i < alnScores.size(); i++){
                        if(alnScores[i]>= curScore && abs(alnBegins[i]-tmpRB) < editDist)
                            equal = true;
                    }
                    if(!equal){
                        alignment_t * alignmentFound = new alignment_t();
                        alignmentFound->refB = tmpRB;
                        alignmentFound->refE = tmpRE;
                        alignmentFound->qB = tmpQB;
                        alignmentFound->qE = tmpQE;
                        alignmentFound->alignmentScore = curScore;
                        if(!fw){
                            alignmentFound->flag.set(4,true);
                        }
                        alignmentFound->refLength = (tmpRE-tmpRB+1);
                        for(std::list<size_t>::const_iterator iterator = chain.begin(); iterator != chain.end(); iterator++){
                            match_t seed = matches[*iterator];
                            alignmentFound->chain.push_back(seed);
                        }
                        read.alignments.push_back(alignmentFound);
                        alnScores.push_back(curScore);
                        alnBegins.push_back(tmpRB);
                    }
                }
                else if(alnOptions.print >= 2)
                    cerr << "#FAIL " << curScore << endl;
                if(alnOptions.print >= 5) cerr << "final score is " << curScore << endl;
                alnBounds.push_back(boundaries(tmpRB, tmpRE, tmpQB, tmpQE));
#ifndef NDEBUG
                if(alnOptions.debugFile != NULL){
                    //print alignment
                    fprintf(alnOptions.debugFile,"@CHA %lu %d %d\n", chain.size(), chainScore, anchor.len);
                    for( list<size_t>::const_iterator iterator = chain.begin(); iterator != chain.end(); iterator++){
                        fprintf(alnOptions.debugFile,"%d %ld %d\n", matches[*iterator].query, matches[*iterator].ref, matches[*iterator].len);
                    }
                    if(curScore >= (int) (P.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist)){
                        alignment_t * alignmentFound = new alignment_t();
                        alignmentFound->refB = tmpRB;
                        alignmentFound->refE = tmpRE;
                        alignmentFound->qB = tmpQB;
                        alignmentFound->qE = tmpQE;
                        alignmentFound->alignmentScore = curScore;
                        for(std::list<size_t>::const_iterator iterator = chain.begin(); iterator != chain.end(); iterator++){
                            match_t seed = matches[*iterator];
                            alignmentFound->chain.push_back(seed);
                        }
                        uint8_t *query;
                        query = (uint8_t *)malloc(P.length());
                        for (size_t i = 0; i < P.length(); ++i) // convert to the nt4 encoding
                            query[i] = Utils::ORDVALUE[(size_t)P[i]];
                        alignmentFound->constructCIGAR(S, query, P.length(), alnOptions.fixedBandSize, alnOptions.mat, 
                                alnOptions.match, alnOptions.openGap, alnOptions.extendGap, 0);
                        free(query);
                        fprintf(alnOptions.debugFile,"@SUC %d %ld %s %d\n", alignmentFound->editDist, tmpRB, alignmentFound->NMtag.c_str(), curScore);
                        delete alignmentFound;
                    }
                    else{//print fail reason
                        fprintf(alnOptions.debugFile,"@ERR %ld %d\n", tmpRB, curScore);
                    }
                    //print best possible alignment
                    long refB = matches[chain.front()].ref - matches[chain.front()].query - 3*editDist;
                    long refE = refB+P.length() + 6*editDist;
                    dp_output output = dpInRegion(S, P, refB, refE, 0, P.length());
                    string cigarstring = createCigar(output.cigarChars, output.cigarLengths);
                    fprintf(alnOptions.debugFile,"@BEST %d %ld %s %d\n", output.editDist, refB, cigarstring.c_str(), output.dpScore);
                }
#endif
            }
        }
    }
    int scoreReturned = -1;
    if(maxAlnScore >= minScore)
        scoreReturned = maxAlnScore;
    return scoreReturned;
}

typedef struct {
    long q, l;
} q_cov_t;

#define cov_cmp(a, b) (((b).q < (a).q) - ((a).q < (b).q))
KBTREE_INIT(qCov, q_cov_t, cov_cmp)

bool hasUnique(const lis_t & lis, kbtree_t(qCov) *tree){
    const vector<match_t> & seeds = lis.matches;
    bool unique = false;
    for(size_t i = 0; i < seeds.size(); i++){
        q_cov_t seed;
        seed.q = seeds[i].query;
        seed.l = seeds[i].len;
        if (kb_size(tree)) {
            q_cov_t *lower, *upper;
            kb_interval(qCov, tree, seed, &lower, &upper); // find the closest read pos that is already covered
            if(!lower || lower->q+lower->l < seed.q+seed.l){
                if(lower && lower->q == seed.q) kb_del(qCov, tree, *lower);
                kb_put(qCov, tree, seed);
                unique = true;
            }
        }
        else{
            kb_put(qCov, tree, seed);
            unique = true;
        }
    }
    return unique;
}

bool compIntervals(const lis_t & i,const lis_t & j){
    return (i.score > j.score);
}

void processLisIntervals(const sparseSA& sa, read_t & read, const align_opt & alnOptions, 
        const vector<lis_t> & lisIntervals, int editDist, int maxAlnCount, int minCoverage, int minScore){
    //for every interval, try to align
    size_t lisIndex = 0;
    //trial parameter for performance reasons
    int trialF = 0;
    int trialR = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    //tree containing information for uniqueness option
    kbtree_t(qCov) *treeF;
    kbtree_t(qCov) *treeR;
    treeF = kb_init(qCov, KB_DEFAULT_SIZE);
    treeR = kb_init(qCov, KB_DEFAULT_SIZE);
    int maxAlnScore = 0;
    while((int)(read.alignments.size()) < maxAlnCount && lisIndex < lisIntervals.size()){
        bool unique = hasUnique(lisIntervals[lisIndex], lisIntervals[lisIndex].fw ? treeF : treeR);
        bool extend = unique;
        int trial = lisIntervals[lisIndex].fw ? trialF : trialR;
        bool extensionCriterium = (trial < alnOptions.maxTrial || extend) && (lisIntervals[lisIndex].matches.size() > 1 || lisIntervals[lisIndex].score >= minCoverage);
#ifndef NDEBUG
        if(alnOptions.debugFile != NULL){
            int Plength = read.sequence.length();
            printLisData(sa, lisIntervals[lisIndex], alnOptions, lisIndex, extensionCriterium, unique, trial < alnOptions.maxTrial, Plength, editDist);
        }
#endif
        if(alnOptions.print >=2 ){
            cerr << "#LISI " << lisIndex << endl;
            cerr << "#UNIQUE " << unique << endl;
            cerr << "#TRIAL " << trial << endl;
            cerr << "#EXTEND " << extensionCriterium << endl;
            cerr << "#LISMEMS " << lisIntervals[lisIndex].matches.size() << endl;
            cerr << "#LISSCORE " << lisIntervals[lisIndex].score << endl;
        }
        if(extensionCriterium){
            //sort matches in this interval according to query position
            const lis_t & cluster = lisIntervals[lisIndex];
            const vector<match_t> & matchVector = cluster.matches;
            long chrStart, chrEnd;
            sa.getChromBounds(matchVector[0].ref, chrStart, chrEnd);
            if(alnOptions.print >=5) cerr << "try cluster " << lisIndex << " with " << matchVector.size() << " seeds on " << (cluster.fw ? "forward" : "reverse complemented") << " strand." << endl;
            int returnedAlnScore = extendAlignmentChain(sa.S, cluster.fw ? read.sequence : read.rcSequence, matchVector, editDist, minScore, minCoverage, alnOptions, chrStart, chrEnd, read, cluster.fw);
            //deal with new alignments
            if(cluster.fw) trialF++;
            else trialR++;
            if(returnedAlnScore >= 0 && (!alnOptions.maxTrialBest || maxAlnScore < returnedAlnScore)){
                if(cluster.fw) trialF = 0;
                else trialR = 0;
            }
            maxAlnScore = max(maxAlnScore, returnedAlnScore);
        }
        lisIndex++;
    }
#ifndef NDEBUG
    if(alnOptions.debugFile != NULL){
        int Plength = read.sequence.length();
        fprintf(alnOptions.debugFile,"$E\n");
        size_t lisAfterIndex =  lisIndex;
        while( lisAfterIndex < lisIntervals.size()){
            printLisData(sa, lisIntervals[lisAfterIndex], alnOptions, lisAfterIndex, false, false, false, Plength, editDist);
            lisAfterIndex++;
        }
        fprintf(alnOptions.debugFile,"#END\n");
    }
#endif         
    kb_destroy(qCov, treeF);
    kb_destroy(qCov, treeR);
    if(alnOptions.print >=5){
        cerr << "stopped matching because: " << endl; 
        if((int)(read.alignments.size()) >= maxAlnCount) cerr << " enough alignments were found."<< endl;
        else if(lisIndex >= lisIntervals.size()) cerr << " no more clusters are available."<< endl;
        else if(lisIntervals[lisIndex].score <= minCoverage) cerr << " no new clusters reach minimum query coverage."<< endl;
        else if(min(trialF,trialR) >= alnOptions.maxTrial) cerr << " max number of extensions without result have been reached."<< endl;
        if(lisIndex < lisIntervals.size()){
            cerr << "remaining clusters map approx to:" << endl;
            for(size_t i = lisIndex; i < lisIntervals.size(); i++){
                //sort this candidate region by query position
                const lis_t & cluster = lisIntervals[i];
                const vector<match_t> & matchVector = cluster.matches;
                long chrStart, chrEnd;
                cerr << "cluster " << i << " with " << matchVector.size() << " seeds and " << cluster.score << " bases covered" << endl;
                sa.getChromBounds(matchVector[0].ref, chrStart, chrEnd);
            }
        }
    }
}

void filterLISwithAlignments(read_t & read, vector<lis_t> & lisIntervals){
    for(size_t i = 0; i < read.alignments.size(); i++){
        alignment_t * a = read.alignments[i];
        vector< lis_t >::iterator it = lisIntervals.begin();
        while(it != lisIntervals.end()) {
            //If alignment overlaps with LIS
            if(a->flag.test(4) != it->fw && a->refE >= it->refBegin && a->refB <= it->refEnd) 
                it = lisIntervals.erase(it);
            else 
                ++it;
        }
    }
}

void calculateCandRegions(const sparseSA& sa, const read_t & read, const align_opt & alnOptions, vector<lis_t> & lisIntervals, double errorPercent){
    const string & P = read.sequence;
    const string & Prc = read.rcSequence;
    int Plength = P.length();
    int editDist = (int)(errorPercent*Plength)+1;
    int min_len = alnOptions.minMemLength;
    if(!alnOptions.fixedMinLength && Plength >= 1000){
        min_len += (20*((int)(Plength/500)))/((int)100*errorPercent);
    }
    if(alnOptions.fixedMinLength && Plength/editDist > min_len){
        min_len = Plength/editDist;
    }
    ////////////////////////////
    //FIND SEEDS
    ///////////////////////////
    vector<match_t> matches;
    vector<match_t> matchesRC;
    //calc seeds
    if(alnOptions.print >= 4) cerr << "calculate seeds" << endl;
    calculateSeedsBothStrands(sa, P, Prc, min_len, alnOptions, matches, matchesRC);
#ifndef NDEBUG
    if(alnOptions.debugFile != NULL){
        printSimRegion(sa, read, alnOptions, editDist);
        fprintf(alnOptions.debugFile,"#MEM %lu %lu\n", matches.size(), matchesRC.size());
        fprintf(alnOptions.debugFile,"#MEMF\n");
        for(size_t i = 0; i < matches.size(); i++)
            fprintf(alnOptions.debugFile,"%d %ld %d\n", matches[i].query, matches[i].ref, matches[i].len);
        fprintf(alnOptions.debugFile,"#MEMR\n");
        for(size_t i = 0; i < matchesRC.size(); i++)
            fprintf(alnOptions.debugFile,"%d %ld %d\n", matchesRC[i].query, matchesRC[i].ref, matchesRC[i].len);
    }
#endif
    if(matches.size()>0 || matchesRC.size()>0){
        /////////////////////////
        //FIND CANDIDATE REGIONS
        /////////////////////////
        if(alnOptions.print >= 4) cerr << "calculate clusters" << endl;
        identifyAndFilter(matches, true, Plength, editDist, lisIntervals, alnOptions);
        size_t lisForward = lisIntervals.size();
        if(alnOptions.print >= 2) cerr << "#LISF " << lisIntervals.size() << endl;
        identifyAndFilter(matchesRC, false, Plength, editDist, lisIntervals, alnOptions);
        if(alnOptions.print >= 2) cerr << "#LISR " << (lisIntervals.size()-lisForward) << endl;
        //fair calculation of len: coverage
        for(size_t i = 0; i < lisIntervals.size(); i++)
            lisIntervals[i].score = coveredBases(lisIntervals[i]);
        sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
    }
}

void unpairedMatch(const sparseSA& sa, read_t & read,const align_opt & alnOptions){
    /////////////////////////
    //FIND CANDIDATE REGIONS
    /////////////////////////
    vector<lis_t> lisIntervals;
    calculateCandRegions(sa, read, alnOptions, lisIntervals, alnOptions.errorPercent);
    if(!lisIntervals.empty()){
        //FILTER LIS WITH EXISTING ALIGNMENTS
        if(!read.alignments.empty())
            filterLISwithAlignments(read, lisIntervals);
        int Plength = read.sequence.length();
        int editDist = (int)(alnOptions.errorPercent*Plength)+1;
        int minScore = (int) (Plength*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist);
        processLisIntervals(sa, read, alnOptions, lisIntervals, editDist, alnOptions.maxAlnCount, Plength*alnOptions.minCoverage, minScore);
        if(read.alignments.empty() && alnOptions.rescue){
            if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
            minScore = max(0,(int)(Plength*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength));
            int maxAlignmentCount = 20;
            int minCoverage = Plength*0.01;
            processLisIntervals(sa, read, alnOptions, lisIntervals, editDist, maxAlignmentCount, minCoverage, minScore);
        }
        if(read.alignments.empty() && alnOptions.rescue){
            if(alnOptions.print >= 2) cerr << "#RESCUE2" << endl;
            lisIntervals.clear();
            calculateCandRegions(sa, read, alnOptions, lisIntervals, 2.0*alnOptions.errorPercent);
            if(!lisIntervals.empty()){
                minScore = max(0,(int)(Plength*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength));
                int maxAlignmentCount = 20;
                int minCoverage = Plength*0.01;
                processLisIntervals(sa, read, alnOptions, lisIntervals, editDist, maxAlignmentCount, minCoverage, minScore);
            }
        }
    }
}

////////////////////
//PAIRED END METHODS
////////////////////

/**
 * returns if two alignments are concordantly paired or not
 */
bool isConcordant(const alignment_t* mate1, const alignment_t* mate2, const paired_opt& options, int print){
    assert(mate1->rname.compare("*") != 0);
    assert(mate2->rname.compare("*") != 0);
    if(print >= 4) cerr << "check concordancy" << endl;
    bool firstLeft = true;
    //Check same chromosome
    if(mate1->rname != mate2->rname){
        if(print >= 4) cerr << "different chromosomes" << endl;
        return false;
    }
    //Check strand
    if(options.orientation == PAIR_FF){
        if(mate1->flag.test(4)!=mate2->flag.test(4)){
            if(print >= 4) cerr << "different orientation" << endl;
            return false;
        }
        else
            firstLeft = !mate1->flag.test(4);//if rc must be right
    }
    else if(options.orientation == PAIR_FR){
        if(mate1->flag.test(4)==mate2->flag.test(4)){
            if(print >= 4) cerr << "different orientation" << endl;
            return false;
        }
        else
            firstLeft = !mate1->flag.test(4);//if rc must be right
    }
    else{
        if(mate1->flag.test(4)==mate2->flag.test(4)){
            if(print >= 4) cerr << "different orientation" << endl;
            return false;
        }
        else
            firstLeft = mate1->flag.test(4);//if rc must be left
    }
    if(firstLeft && mate2->refE <= mate1->refB){
        if(print >= 4) cerr << "bad order of reads" << endl;
        return false;
    }
    else if(!firstLeft && mate1->refE <= mate2->refB){
        if(print >= 4) cerr << "bad order of reads" << endl;
        return false;
    }
    //check fragment insert size
    long begin = min(mate1->refB,mate2->refB);
    long end = max(mate1->refE,mate2->refE);
    long size = end - begin +1L;
    if(size < options.minInsert || size > options.maxInsert){
        if(print >= 2) cerr << "#INSERTERR " << size << endl;
        if(print >= 4) cerr << "bad insert: " << size << " is not in interval " << options.minInsert << ", " << options.maxInsert << endl;
        return false;
    }
    long begin1 = mate1->pos;
    long end1 = begin1 + (long)mate1->refLength-1L;
    long begin2 = mate2->pos;
    long end2 = begin2 + (long)mate2->refLength-1L;
    //check contain
    bool contain = (begin1 >= begin2 && end1 <= end2) || (begin2 >= begin1 && end2 <= end1);
    if(contain && !options.contain){
        if(print >= 4) cerr << "is contained" << endl;
        return false;
    }
    //check overlap
    bool overlap = contain || (begin2 >= begin1 && begin2 <= end1) || (end2 >= begin1 && end2 <= end1);
    if(!options.overlap && overlap){
        if(print >= 4) cerr << "is overlapping" << endl;
        return false;
    }
    //check dovetail
    if((firstLeft && (begin2 < begin1)) || (!firstLeft && (begin1 < begin2)) ){
        if(!overlap){
            if(print >= 4) cerr << "is dovetail" << endl;
            return false;
        }
        else if(!options.dovetail){
            if(print >= 4) cerr << "is dovetail" << endl;
            return false;
        }
    }
    if(print >= 2) cerr << "#CONCORDANT" << endl;
    return true;
}

alignment_t * unpairedCopy(const alignment_t * a){
  alignment_t * newAln = new alignment_t();
  newAln->refB = a->refB;
  newAln->refE = a->refE;
  newAln->qB = a->qB;
  newAln->qE = a->qE;
  newAln->pos = a->pos;
  newAln->rname = a->rname;
  newAln->refLength = a->refLength;
  newAln->alignmentScore = a->alignmentScore;
  newAln->flag.set(4, a->flag.test(4));
  for(size_t i = 0; i < a->chain.size(); i++)
      newAln->chain.push_back(a->chain[i]);
  return newAln;
}

void setPaired(alignment_t* mate1, alignment_t* mate2, bool concordant){
    assert(!mate1->flag.test(0));
    assert(!mate2->flag.test(0));
    mate1->flag.set(0, true);
    mate2->flag.set(0, true);
    mate1->setMate(mate2, concordant, true);
    mate2->setMate(mate1, concordant, false);
}

void setPaired(const sparseSA& sa, read_t & mate1, read_t & mate2, int i, int j, bool concordant){
    mate1.alignments[i]->setLocalPos(sa);
    mate2.alignments[j]->setLocalPos(sa);
    int mate1index = i;
    int mate2index = j;
    if(mate1.alignments[mate1index]->paired()){
        alignment_t * newAln = unpairedCopy(mate1.alignments[mate1index]);
        mate1.alignments.push_back(newAln);
        mate1index = mate1.alignments.size()-1;
    }
    if(mate2.alignments[mate2index]->paired()){
        alignment_t * newAln = unpairedCopy(mate2.alignments[mate2index]);
        mate2.alignments.push_back(newAln);
        mate2index = mate2.alignments.size()-1;
    }
    mate1.pairedAlignmentCount++;
    mate2.pairedAlignmentCount++;
    setPaired(mate1.alignments[mate1index], mate2.alignments[mate2index], concordant);
}

void calculateDiscordantAlignments(const sparseSA& sa, read_t & mate1, read_t & mate2){
    if(!mate1.alignments.empty() && !mate2.alignments.empty()){
        sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignments);
        sort(mate2.alignments.begin(),mate2.alignments.end(), compAlignments);
        //FIND MAX SCORE OF BOTH MATES
        int maxScoreFirst = mate1.alignments[0]->alignmentScore;
        int maxScoreSecond = mate2.alignments[0]->alignmentScore;
        //DISCORDANT ONLY BETWEEN UNPAIRED MAX SCORE ALIGNMENTS
        size_t mate1AlnCount = mate1.alignments.size();
        size_t mate2AlnCount = mate2.alignments.size();
        size_t i = 0;
        while(i < mate1AlnCount && mate1.alignments[i]->alignmentScore == maxScoreFirst){
            if(!mate1.alignments[i]->concordant()){
                size_t j = i;
                while(j < mate2AlnCount && mate2.alignments[j]->alignmentScore == maxScoreSecond){
                    if(!mate2.alignments[j]->concordant())
                        setPaired(sa, mate1,mate2,i,j,false);
                    j++;
                }
            }
            i++;
        }
    }
}

bool dpWindow(long mateRefB, long mateRefE, bool mateIsFirst, bool mateIsFw, 
        int Plength, const paired_opt& options, long & refB, long & refE, 
        bool& otherFW){
    //FIND VALUE OF OTHER ALNS BEGIN, ORIENTATION and IF ALN IS BEFORE OR AFTER THIS ONE
    long refBmate = mateIsFw ? mateRefB : mateRefE;
    bool rightWindow;
    if(options.orientation == PAIR_FF){
        rightWindow = mateIsFirst ? mateIsFw : !mateIsFw;
        otherFW = mateIsFw;
    }
    else if(options.orientation == PAIR_FR){
        rightWindow = mateIsFw;
        otherFW = !mateIsFw;
    }
    else{
        rightWindow = !mateIsFw;
        otherFW = !mateIsFw;
    }
    //FIND BASIC WINDOW VALUES
    if(rightWindow){
        refB = refBmate + options.minInsert;
        refE = refBmate + options.maxInsert;
    }
    else{
        refB = refBmate - options.maxInsert;
        refE = refBmate - options.minInsert;
    }
    if(otherFW)
        refE += Plength;
    else
        refB -= Plength;
    //ADJUST TO OVERLAP AND DOVETAIL
    if(!options.overlap){
        if(rightWindow)
            refB = max(refB, mateRefE);
        else
            refE = min(refE, mateRefB);
    }
    if(!options.dovetail){
        if(rightWindow)
            refB = max(refB, mateRefB);
        else
            refE = min(refE, mateRefE);
    }
    //RETURN IF WINDOW IS LARGE ENOUGH
    return refE > refB;
}

struct alnScoreInfo_t {
  bool firstMate;
  int position; // alignment in second mate
  int alignmentScore;
  alnScoreInfo_t(): firstMate(1), position(0), alignmentScore(0) {}
  alnScoreInfo_t(bool fMate, int pos, int alnScore): firstMate(fMate), position(pos), alignmentScore(alnScore) {}
};

bool compAlnScorePaired(const alnScoreInfo_t & i,const alnScoreInfo_t & j){//pairedMatch 1
    return i.alignmentScore > j.alignmentScore;
}

//FIND ALIGNMENT IN MATE THAT FITS WINDOW refB, refE, fw and has highest alnScore
//Optimize to not work with vector, but search tree
int bestAlignmentInWindow(const read_t & mate, const vector<alnScoreInfo_t> & mateAlns, 
        bool firstMate, long refB, long refE, bool fw){
    int bestAln = -1;
    int bestAlnScore = -1*mate.sequence.size();
    for(size_t i = 0; i < mateAlns.size(); i++){
        if(firstMate != mateAlns[i].firstMate){
            alignment_t * mateAln = mate.alignments[mateAlns[i].position];
            if(mateAln->refB >= refB && mateAln->refE <= refE && mateAln->flag.test(4) != fw && mateAln->alignmentScore > bestAlnScore){
                bestAln = mateAlns[i].position;
                bestAlnScore = mateAln->alignmentScore;
            }
        }
    }
    return bestAln;
}

void pairAlignmentToWindow(const sparseSA& sa, read_t & mate1, read_t & mate2, 
        const align_opt & alnOptions, const paired_opt & pairedOpt, vector<alnScoreInfo_t> & alnMateAlns, 
        vector<alnScoreInfo_t> & windowMateAlns, int & alnCount, int maxAlnCount){
    size_t i = 0;
    int trial = 0;
    while(i < alnMateAlns.size() && alnCount < maxAlnCount && trial < alnOptions.maxTrial){
        alignment_t * aln = alnMateAlns[i].firstMate ? mate1.alignments[alnMateAlns[i].position] : mate2.alignments[alnMateAlns[i].position]; 
        if(!aln->paired()){
            if(alnOptions.print >= 4) cerr << "find window for aligment of mate " << (alnMateAlns[i].firstMate ? "1" : "0") << " on strand " << (aln->flag.test(4) ? "R" : "F") << " and position " << aln->refB << endl;
            trial++;
            //FIND WINDOW
            long refB, refE;
            bool fw;
            read_t & windowMate = alnMateAlns[i].firstMate ? mate2 : mate1;
            int Plength = windowMate.sequence.length();
            dpWindow(aln->refB, aln->refE, alnMateAlns[i].firstMate, !aln->flag.test(4), 
                    Plength, pairedOpt, refB, refE, fw);
            if(alnOptions.print >= 4) cerr << "window is located on " << (fw ? "F" : "R") << " strand and position has interval " << refB << "," << refE << endl;
            //CHECK IF AN ALIGNMENT ALREADY EXISTS IN THIS WINDOW
            int windowMatePos = bestAlignmentInWindow(windowMate, windowMateAlns, alnMateAlns[i].firstMate, refB, refE, fw);
            if(alnOptions.print >= 4) cerr << "an alignment of the mate was found as the " << windowMatePos << " 'th alignment of the mate" << endl;
            if(alnOptions.print >= 2) cerr << "#WINDOWLIST " << (windowMatePos>=0) << endl;
            if(windowMatePos < 0){
                //TRY TO ALIGN IN WINDOW    
                int tLen = refE-refB+1;
                uint8_t * target = (uint8_t *)malloc(tLen);
                for (int j = 0; j < tLen; ++j) // convert to the nt4 encoding
                    target[j] = Utils::ORDVALUE[(size_t)sa.S[refB+j]];
                int xtra = KSW_XSTART | (Plength * alnOptions.match < 250? KSW_XBYTE : 0);
                kswr_t alnResult = ksw_align(Plength, fw ? windowMate.readP : windowMate.rcReadP, tLen, 
                        target, 5, alnOptions.mat, alnOptions.openGap, alnOptions.extendGap, xtra, 0);
                free(target);
                //if succes, set paired, trial back to 0
                int editDist = (int)(alnOptions.errorPercent*Plength)+1;
                int minScore = (int) (Plength*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist);
                if(alnOptions.print >= 4) cerr << "an alignment was produced with score " << alnResult.score << ". Min score is " << minScore << endl;
                if(alnResult.score > minScore){
                    alignment_t * newAln = new alignment_t();
                    newAln->alignmentScore = alnResult.score;
                    newAln->qB = alnResult.qb;
                    newAln->qE = alnResult.qe;
                    newAln->refB = refB + alnResult.tb;
                    newAln->refE = refB + alnResult.te;
                    newAln->refLength = alnResult.te - alnResult.tb + 1;
                    newAln->flag.set(4, !fw);
                    windowMate.alignments.push_back(newAln);
                    windowMatePos = windowMate.alignments.size()-1;
                    windowMateAlns.push_back(alnScoreInfo_t(!alnMateAlns[i].firstMate,windowMatePos,alnResult.score));
                }
                if(alnOptions.print >= 2) cerr << "#WINDOWALIGN " << (windowMatePos >= 0) << endl;
            }
            if(windowMatePos >= 0){
                //PAIR ALIGNMENTS, SUCCES
                if(alnOptions.print >= 4) cerr << "found concordancy between 2 alns." << endl;
                if(alnMateAlns[i].firstMate)
                    setPaired(sa, mate1, mate2, alnMateAlns[i].position,windowMatePos, true);
                else
                    setPaired(sa, mate1, mate2,windowMatePos, alnMateAlns[i].position, true);
                trial = 0;
                alnCount++;
            }
        }
        i++;
    }
}

void pairAlignmentsOrdered(const sparseSA& sa, read_t & mate1, read_t & mate2, 
        const align_opt & alnOptions, const paired_opt & pairedOpt, 
        const vector<alnScoreInfo_t> & pairedOrder, int & alnCount, int maxPairedAlignments){
    if(alnOptions.print >= 2) cerr << "#COUPLEALN " << pairedOrder.size() << endl;
    size_t i = 0;
    //TRY AND PAIR THE ALIGNMENTS
    while(i < pairedOrder.size() && alnCount < maxPairedAlignments){
        size_t j = i+1;
        bool paired = pairedOrder[i].firstMate ? mate1.alignments[pairedOrder[i].position]->paired() : mate2.alignments[pairedOrder[i].position]->paired();
        if(!paired){
            while(j < pairedOrder.size() && alnCount < maxPairedAlignments){
                if(pairedOrder[i].firstMate != pairedOrder[j].firstMate){
                    int mate1Index = pairedOrder[i].firstMate ? pairedOrder[i].position : pairedOrder[j].position;
                    int mate2Index = pairedOrder[i].firstMate ? pairedOrder[j].position : pairedOrder[i].position;
                    alignment_t * aln1 = mate1.alignments[mate1Index];
                    aln1->setLocalPos(sa);
                    alignment_t * aln2 = mate2.alignments[mate2Index];
                    aln2->setLocalPos(sa);
                    if(isConcordant(aln1, aln2, pairedOpt, alnOptions.print)){
                        if(alnOptions.print >= 4) cerr << "found concordancy between aln " << mate1Index << " of mate 1 and aln " << mate2Index << " of mate 2" << endl;
                        setPaired(sa, mate1, mate2,mate1Index,mate2Index, true);
                        alnCount++;
                    }
                }
                j++;
            }
        }
        i++;
    }
}

//aln both and pair
void pairedMatch1(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    //Calculate the mappings for read 1 and 2
    if(alnOptions.print >= 3) cerr << "matching first mate ..." << endl;
    unpairedMatch(sa, mate1, alnOptions);
    if(alnOptions.print >= 3) cerr << "first mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    if(alnOptions.print >= 3) cerr << "matching second mate ..." << endl;
    unpairedMatch(sa, mate2, alnOptions);
    if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
    int alnCount1 = mate1.alignmentCount();
    int alnCount2 = mate2.alignmentCount();
    int alnCount = 0;
    int maxPairedAlignments = alnOptions.alignmentCount+1;
    //ORDER ALIGNMENTS ACCORDING TO SCORE
    vector<alnScoreInfo_t> pairedOrder;
    pairedOrder.reserve(alnCount1+alnCount2);
    for(size_t i = 0; i < mate1.alignments.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(true, i,mate1.alignments[i]->alignmentScore ));
    for(size_t i = 0; i < mate2.alignments.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(false, i,mate2.alignments[i]->alignmentScore ));
    sort(pairedOrder.begin(),pairedOrder.end(),compAlnScorePaired);
    //TRY AND PAIR THE ALIGNMENTS
    if(alnCount1>0 && alnCount2>0){
        if(alnOptions.print >= 3) cerr << "search concordant alignments of " << pairedOrder.size() << " alignments" << endl;
        pairAlignmentsOrdered(sa, mate1, mate2, alnOptions, pairedOpt, 
        pairedOrder, alnCount, maxPairedAlignments);
        if(alnOptions.print >= 3) cerr << "found " << alnCount << " concordant alignments" << endl;
    }
    //RESCUE PROCEDURE IF NO CONCORDANT ALIGNMENTS WERE FOUND
    if(pairedOpt.pairedRescue && alnCount < maxPairedAlignments && max(alnCount1,alnCount2)>0){
        if(alnOptions.print >= 2) cerr << "#PAIREDRESCUE" << endl;
        if(alnOptions.print >= 3) cerr << "try to recover concordant alignments using windows " << endl;
        vector<alnScoreInfo_t> newAlnFromWindow;
        pairAlignmentToWindow(sa, mate1, mate2, alnOptions, pairedOpt, pairedOrder, 
        newAlnFromWindow, alnCount, maxPairedAlignments);
        if(alnOptions.print >= 3) cerr << "rescue procedure resulted in " << alnCount << " concordant alignments" << endl;
    }
    //FIND DISCORDANT ALIGNMENTS AMONGST HIGHEST SCORING UNPAIRED ALIGNMENTS
    int pairedAlignments = mate1.pairedAlignmentCount;
    if(alnOptions.print >= 2) cerr << "#CONCORDANTCOUNT " << pairedAlignments << endl;
    if(pairedOpt.discordant && alnCount < maxPairedAlignments)
        calculateDiscordantAlignments(sa, mate1, mate2);
    if(alnOptions.print >= 2) cerr << "#DISCORDANTCOUNT " << mate1.pairedAlignmentCount-pairedAlignments << endl;
}

bool compIntervalsPointer(const lis_t * i,const lis_t * j){//pairedMatch 4
    return (i->score > j->score);
}

bool compIntervalsRef(const lis_t & i,const lis_t & j){//pairedMatch 4
    return compMatches(i.matches[0], j.matches[0]);
}

bool isConcordantLIS(const lis_t& mate1, const lis_t& mate2, int editDist1, int editDist2, const paired_opt& options, int print){
    if(print >= 5) cerr << "check concordancy between LIS" << endl;
    //Check strand
    bool firstLeft;
    bool mate1FW = mate1.fw;
    bool mate2FW = mate2.fw;
    if((options.orientation == PAIR_RF || options.orientation == PAIR_FR) && mate1FW == mate2FW){
        if(print >= 5) cerr << "bad orientation" << endl;
        return false;
    }
    else if(options.orientation == PAIR_FF && mate1FW != mate2FW ){
        if(print >= 5) cerr << "bad orientation" << endl;
        return false;
    }
    //check left-right
    if(options.orientation == PAIR_FF || options.orientation == PAIR_FR)
        firstLeft = mate1FW;//if rc must be right
    else
        firstLeft = !mate1FW;//if rc must be left
    if(firstLeft && mate2.refEnd <= mate1.refBegin){
        if(print >= 5) cerr << "bad order of reads" << endl;
        return false;
    }
    else if(!firstLeft && mate1.refEnd <= mate2.refBegin){
        if(print >= 5) cerr << "bad order of reads" << endl;
        return false;
    }
    //check coordinates versus insert size
    long begin = min(mate1.refBegin,mate2.refBegin);
    long end = max(mate1.refEnd,mate2.refEnd);
    long size = end - begin +1L;
    if(size+(editDist1+editDist2) < options.minInsert || size-(editDist1+editDist2) > options.maxInsert){
        if(print >= 2) cerr << "#INSERTERR " << size+(editDist1+editDist2) << endl;
        return false;
    }
    //check contain, dovetail and overlap
    long begin1 = mate1.refBegin;
    long end1 = mate1.refEnd;
    long begin2 = mate2.refBegin;
    long end2 = mate1.refEnd;
    //check contain
    bool contain = (begin1-editDist1 >= begin2+editDist2 && end1+editDist1 <= end2-editDist2) 
    || (begin2-editDist2 >= begin1+editDist1 && end2+editDist2 <= end1-editDist1);
    if(contain && !options.contain){
        if(print >= 5) cerr << "contain" << endl;
        return false;
    }
    //check overlap
    bool overlap = contain || !(end2-editDist2 < begin1+editDist1 || begin2+editDist2 > end1-editDist1);
    if(!options.overlap && overlap){
        if(print >= 5) cerr << "overlap" << endl;
        return false;
    }
    //check dovetail
    if((firstLeft && (begin2+editDist2 < begin1-editDist1)) || (!firstLeft && (begin1+editDist1 < begin2-editDist2)) ){
        if(!overlap){
            if(print >= 5) cerr << "dovetail" << endl;
            return false;
        }
        else if(!options.dovetail){
            if(print >= 5) cerr << "dovetail" << endl;
            return false;
        }
    }
    if(print >= 5) cerr << "concordant pair found" << endl;
    return true;
}

//Calculate LIS of both mates, match together 2 LIS, if match, calculate both alignments
void pairedMatch3(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    int editDist1 = (int)(alnOptions.errorPercent*mate1.sequence.length())+1;
    int editDist2 = (int)(alnOptions.errorPercent*mate2.sequence.length())+1;
    //FIND CANDIDATE REGIONS FOR BOTH MATES
    vector<lis_t>  lisIntervals1;
    vector<lis_t>  lisIntervals2;
    if(alnOptions.print >= 3) cerr << "calculate LIS of both mates"  << endl;
    calculateCandRegions(sa, mate1, alnOptions, lisIntervals1, alnOptions.errorPercent);
    calculateCandRegions(sa, mate2, alnOptions, lisIntervals2, alnOptions.errorPercent);
    //FILTER LIS: ONLY THOSE THAT PAIR   
    vector<lis_t>  lisIntervals1Filtered;
    vector<lis_t>  lisIntervals2Filtered;
    if(alnOptions.print >= 3) cerr << "filter LIS of both mates using pairing"  << endl;
    for(size_t i = 0; i < lisIntervals1.size(); i++){
        for(size_t j = 0; j < lisIntervals2.size(); j++){
            if(!(lisIntervals1[i].extended && lisIntervals2[j].extended) && isConcordantLIS(lisIntervals1[i], lisIntervals2[j], editDist1, editDist2, pairedOpt, alnOptions.print)){
                if(!lisIntervals1[i].extended){
                    lisIntervals1Filtered.push_back(lisIntervals1[i]);
                    lisIntervals1[i].extended = true;
                }
                if(!lisIntervals2[j].extended){
                    lisIntervals2Filtered.push_back(lisIntervals2[j]);
                    lisIntervals2[j].extended = true;
                }
            }
        }
    }
    if(alnOptions.print >= 3) cerr << "number of paired LIS in mate1 " << lisIntervals1Filtered.size()  << endl;
    if(alnOptions.print >= 3) cerr << "number of paired LIS in mate2 " << lisIntervals2Filtered.size()  << endl;
    //RESCUE IF NO PAIRED LIS WERE FOUND
    if(lisIntervals1Filtered.empty() && lisIntervals2Filtered.empty()){
        if(alnOptions.print >= 2) cerr << "#EMPTYFILTEREDLIST"  << endl;
        pairedMatch1(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else{
        //FIND ALIGNMENTS FOR FILTERED LIST OF ALIGNMENTS
        //BE NICE AND GIVE THE OPTION OF LOWER MINSCORE
        int minScore1 = (int) (mate1.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*(editDist1*1.5));
        int minScore2 = (int) (mate2.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*(editDist2*1.5));
        processLisIntervals(sa, mate1, alnOptions, lisIntervals1Filtered, editDist1, alnOptions.maxAlnCount, 0, minScore1);
        processLisIntervals(sa, mate2, alnOptions, lisIntervals2Filtered, editDist2, alnOptions.maxAlnCount, 0, minScore2);
        if(alnOptions.print >= 3) cerr << "first  mate resulted in " << mate1.alignments.size() << " alignments" << endl;
        if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
        //RESCUE
        if(mate1.alignments.empty() && alnOptions.rescue){
            if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
            minScore1 = max(0,(int)(mate1.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*mate1.sequence.length()));
            processLisIntervals(sa, mate1, alnOptions, lisIntervals1Filtered, editDist1, alnOptions.maxAlnCount, 0, minScore1);
        }
        if(mate2.alignments.empty() && alnOptions.rescue){
            if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
            minScore2 = max(0,(int)(mate2.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*mate2.sequence.length()));
            processLisIntervals(sa, mate2, alnOptions, lisIntervals2Filtered, editDist2, alnOptions.maxAlnCount, 0, minScore2);
        }
        //TRY AND PAIR THE ALIGNMENTS
        int alnCount = 0;
        int maxPairedAlignments = alnOptions.alignmentCount+1;
        int alnCount1 = mate1.alignmentCount();
        int alnCount2 = mate2.alignmentCount();
        //ORDER ALIGNMENTS ACCORDING TO SCORE
        vector<alnScoreInfo_t> pairedOrder;
        pairedOrder.reserve(alnCount1+alnCount2);
        for(size_t i = 0; i < mate1.alignments.size(); i++)
            pairedOrder.push_back(alnScoreInfo_t(true, i,mate1.alignments[i]->alignmentScore ));
        for(size_t i = 0; i < mate2.alignments.size(); i++)
            pairedOrder.push_back(alnScoreInfo_t(false, i,mate2.alignments[i]->alignmentScore ));
        sort(pairedOrder.begin(),pairedOrder.end(),compAlnScorePaired);
        //TRY AND PAIR THE ALIGNMENTS
        if(alnCount1>0 && alnCount2>0){
            if(alnOptions.print >= 3) cerr << "search concordant alignments" << endl;
            pairAlignmentsOrdered(sa, mate1, mate2, alnOptions, pairedOpt, 
            pairedOrder, alnCount, maxPairedAlignments);
        }
        //FIND DISCORDANT ALIGNMENTS AMONGST HIGHEST SCORING UNPAIRED ALIGNMENTS
        int pairedAlignments = mate1.pairedAlignmentCount;
        if(alnOptions.print >= 2) cerr << "#CONCORDANTCOUNT " << pairedAlignments << endl;
        if(pairedOpt.discordant && alnCount < maxPairedAlignments){
            if(alnOptions.print >= 3) cerr << "search discordant alignments" << endl;
            calculateDiscordantAlignments(sa, mate1, mate2);
        }
        if(alnOptions.print >= 2) cerr << "#DISCORDANTCOUNT " << mate1.pairedAlignmentCount-pairedAlignments << endl;
    }
}

struct LIScouple_t {
  int firstPosition;
  int secondPosition; // alignment in second mate
  int score;
  LIScouple_t(): firstPosition(0), secondPosition(0), score(0) {}
  LIScouple_t(int fPos, int sPos, int sc): firstPosition(fPos), secondPosition(sPos), score(sc) {}
};

bool compLIScouple(const LIScouple_t & i,const LIScouple_t & j){//pairedMatch 1
    return i.score > j.score;
}

inline void processLisIntervalsSimultaneous(const sparseSA& sa, read_t & mate1, read_t & mate2, 
        const align_opt & alnOptions, const paired_opt & pairedOpt, const vector<lis_t> & lisIntervals1, 
        const vector<lis_t> & lisIntervals2, const vector<LIScouple_t> & pairedOrder, 
        int maxAlnCount, double minCoveragePercent, int minScore1, int minScore2){
    int minCoverage1 = mate1.sequence.length()*minCoveragePercent;
    int minCoverage2 = mate2.sequence.length()*minCoveragePercent;
    int editDist1 = (int)(alnOptions.errorPercent*mate1.sequence.length())+1;
    int editDist2 = (int)(alnOptions.errorPercent*mate2.sequence.length())+1;
    //for every interval, try to align
    size_t lisIndex = 0;
    //trial parameter for performance reasons
    int trial = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    //tree containing information for uniqueness option
    int maxAlnScore = 0;
    int alnCount = 0;
    while(alnCount < maxAlnCount && lisIndex < pairedOrder.size()){
        const lis_t & cluster1 = lisIntervals1[pairedOrder[lisIndex].firstPosition];
        const lis_t & cluster2 = lisIntervals2[pairedOrder[lisIndex].secondPosition];
        bool extensionCriterium = (trial < alnOptions.maxTrial);        
        if(alnOptions.print >=2 ){
            cerr << "#LISI " << lisIndex << endl;
            cerr << "#TRIAL " << trial << endl;
            cerr << "#EXTEND " << extensionCriterium << endl;
            cerr << "#LISMEMS " << (cluster1.matches.size()+cluster2.matches.size()) << endl;
            cerr << "#LISSCORE " << pairedOrder[lisIndex].score << endl;
        }
        if(extensionCriterium){
            const vector<match_t> & matchVector1 = cluster1.matches;
            const vector<match_t> & matchVector2 = cluster2.matches;
            long chrStart1, chrEnd1, chrStart2, chrEnd2;
            sa.getChromBounds(matchVector1[0].ref, chrStart1, chrEnd1);
            sa.getChromBounds(matchVector2[0].ref, chrStart2, chrEnd2);
            int returnedAlnScore, returnedAlnScore1, returnedAlnScore2;
            size_t baseAln1 = mate1.alignments.size();
            size_t baseAln2 = mate2.alignments.size();
            //extend cluster1
            if(alnOptions.print >=5) cerr << "try cluster " << pairedOrder[lisIndex].firstPosition << " with " << matchVector1.size() << " seeds on " << (cluster1.fw ? "forward" : "reverse complemented") << " strand." << endl;
            returnedAlnScore1 = extendAlignmentChain(sa.S, cluster1.fw ? mate1.sequence : mate1.rcSequence, matchVector1, editDist1, minScore1, minCoverage1, alnOptions, chrStart1, chrEnd1, mate1, cluster1.fw);
            //extend cluster2
            if(alnOptions.print >=5) cerr << "try cluster " << pairedOrder[lisIndex].secondPosition << " with " << matchVector2.size() << " seeds on " << (cluster2.fw ? "forward" : "reverse complemented") << " strand." << endl;
            returnedAlnScore2 = extendAlignmentChain(sa.S, cluster2.fw ? mate2.sequence : mate2.rcSequence, matchVector2, editDist2, minScore2, minCoverage2, alnOptions, chrStart2, chrEnd2, mate2, cluster2.fw);
            //rescue cluster1 or 2
            if(baseAln1 == mate1.alignments.size() && baseAln2 < mate2.alignments.size() && alnOptions.rescue){
                if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
                minScore1 = max(0,(int)(mate1.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*mate1.sequence.length()));
                returnedAlnScore1 = extendAlignmentChain(sa.S, cluster1.fw ? mate1.sequence : mate1.rcSequence, matchVector1, editDist1, minScore1, minCoverage1, alnOptions, chrStart1, chrEnd1, mate1, cluster1.fw);
            }
            else if(baseAln1 < mate1.alignments.size() && baseAln2 == mate2.alignments.size() && alnOptions.rescue){
                if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
                minScore2 = max(0,(int)(mate2.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*mate2.sequence.length()));
                returnedAlnScore2 = extendAlignmentChain(sa.S, cluster2.fw ? mate2.sequence : mate2.rcSequence, matchVector2, editDist2, minScore2, minCoverage2, alnOptions, chrStart2, chrEnd2, mate2, cluster2.fw);
            }
            //pair if possible
            if(baseAln1 < mate1.alignments.size() && baseAln2 < mate2.alignments.size()){
                size_t mate1MaxSize = mate1.alignments.size();
                size_t mate2MaxSize = mate2.alignments.size();
                for(size_t i = baseAln1; i < mate1MaxSize; i++){
                    for(size_t j = baseAln2; j < mate2MaxSize; j++){
                        alignment_t * aln1 = mate1.alignments[i];
                        aln1->setLocalPos(sa);
                        alignment_t * aln2 = mate2.alignments[j];
                        aln2->setLocalPos(sa);
                        if(isConcordant(aln1, aln2, pairedOpt, alnOptions.print)){
                            if(alnOptions.print >= 4) cerr << "found concordancy between aln " << i << " of mate 1 and aln " << j << " of mate 2" << endl;
                            setPaired(sa, mate1, mate2,i,j, true);
                            alnCount++;
                        }
                    }
                }
            }
            //deal with new alignments
            trial++;
            returnedAlnScore = returnedAlnScore1 + returnedAlnScore2;
            if(returnedAlnScore >= 0 && (!alnOptions.maxTrialBest || maxAlnScore < returnedAlnScore))
                trial = 0;
            maxAlnScore = max(maxAlnScore, returnedAlnScore);
        }
        lisIndex++;
    }  
}

void pairedMatch6(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    int editDist1 = (int)(alnOptions.errorPercent*mate1.sequence.length())+1;
    int editDist2 = (int)(alnOptions.errorPercent*mate2.sequence.length())+1;
    //FIND CANDIDATE REGIONS FOR BOTH MATES
    vector<lis_t>  lisIntervals1;
    vector<lis_t>  lisIntervals2;
    if(alnOptions.print >= 3) cerr << "calculate LIS of both mates"  << endl;
    calculateCandRegions(sa, mate1, alnOptions, lisIntervals1, alnOptions.errorPercent);
    calculateCandRegions(sa, mate2, alnOptions, lisIntervals2, alnOptions.errorPercent);
    //FILTER LIS: ONLY THOSE THAT PAIR : PUT IN LIST OF LISCOUPLE 
    vector<LIScouple_t>  lisIntervalsFiltered;
    if(alnOptions.print >= 3) cerr << "filter LIS of both mates using pairing"  << endl;
    for(size_t i = 0; i < lisIntervals1.size(); i++){
        for(size_t j = 0; j < lisIntervals2.size(); j++){
            if(!(lisIntervals1[i].extended && lisIntervals2[j].extended) && isConcordantLIS(lisIntervals1[i], lisIntervals2[j], editDist1, editDist2, pairedOpt, alnOptions.print))
                lisIntervalsFiltered.push_back(LIScouple_t(i,j,lisIntervals1[i].score+lisIntervals2[j].score));
        }
    }
    if(alnOptions.print >= 3) cerr << "number of paired LIS in filtered list " << lisIntervalsFiltered.size()  << endl;
    //RESCUE IF NO PAIRED LIS WERE FOUND
    if(lisIntervalsFiltered.empty()){
        if(alnOptions.print >= 2) cerr << "#EMPTYFILTEREDLIST"  << endl;
        pairedMatch1(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else{
        sort(lisIntervalsFiltered.begin(),lisIntervalsFiltered.end(),compLIScouple);
        //FIND ALIGNMENTS FOR FILTERED LIST OF ALIGNMENTS
        //BE NICE AND GIVE THE OPTION OF LOWER MINSCORE
        int minScore1 = (int) (mate1.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*(editDist1*1.5));
        int minScore2 = (int) (mate2.sequence.length()*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*(editDist2*1.5));
        processLisIntervalsSimultaneous(sa, mate1, mate2, alnOptions, pairedOpt, lisIntervals1, lisIntervals2, lisIntervalsFiltered, 
                alnOptions.maxAlnCount, 0, minScore1, minScore2);
        if(alnOptions.print >= 3) cerr << "first  mate resulted in " << mate1.alignments.size() << " alignments" << endl;
        if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
        //FIND DISCORDANT ALIGNMENTS AMONGST HIGHEST SCORING UNPAIRED ALIGNMENTS
        int maxPairedAlignments = alnOptions.alignmentCount+1;
        int pairedAlignments = mate1.pairedAlignmentCount;
        if(alnOptions.print >= 2) cerr << "#CONCORDANTCOUNT " << pairedAlignments << endl;
        if(pairedOpt.discordant && pairedAlignments < maxPairedAlignments){
            if(alnOptions.print >= 3) cerr << "search discordant alignments" << endl;
            calculateDiscordantAlignments(sa, mate1, mate2);
        }
        if(alnOptions.print >= 2) cerr << "#DISCORDANTCOUNT " << mate1.pairedAlignmentCount-pairedAlignments << endl;
    }
}

inline void processLisIntervalsPaired(const sparseSA& sa, read_t & mate1, read_t & mate2, 
        const align_opt & alnOptions, const vector<lis_t> & lisIntervals1, 
        const vector<lis_t> & lisIntervals2, const vector<alnScoreInfo_t> & pairedOrder, 
        int maxAlnCount, double minCoveragePercent, int minScore1, int minScore2){
    int minCoverage1 = mate1.sequence.length()*minCoveragePercent;
    int minCoverage2 = mate2.sequence.length()*minCoveragePercent;
    int editDist1 = (int)(alnOptions.errorPercent*mate1.sequence.length())+1;
    int editDist2 = (int)(alnOptions.errorPercent*mate2.sequence.length())+1;
    //for every interval, try to align
    size_t lisIndex = 0;
    //trial parameter for performance reasons
    int trialF = 0;
    int trialR = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    //tree containing information for uniqueness option
    kbtree_t(qCov) *treeF;
    kbtree_t(qCov) *treeR;
    treeF = kb_init(qCov, KB_DEFAULT_SIZE);
    treeR = kb_init(qCov, KB_DEFAULT_SIZE);
    int maxAlnScore = 0;
    int alnCount = 0;
    while(alnCount < maxAlnCount && lisIndex < pairedOrder.size()){
        const lis_t & cluster = pairedOrder[lisIndex].firstMate ? lisIntervals1[pairedOrder[lisIndex].position] : lisIntervals2[pairedOrder[lisIndex].position];
        bool unique = hasUnique(cluster, cluster.fw ? treeF : treeR);
        bool extend = unique;
        int trial = cluster.fw ? trialF : trialR;
        bool highCoverage = pairedOrder[lisIndex].firstMate ? cluster.score >= minCoverage1 : cluster.score >= minCoverage2;
        bool extensionCriterium = (trial < alnOptions.maxTrial || extend) && (cluster.matches.size() > 1 || highCoverage);        
        if(alnOptions.print >=2 ){
            cerr << "#LISI " << lisIndex << endl;
            cerr << "#UNIQUE " << unique << endl;
            cerr << "#TRIAL " << trial << endl;
            cerr << "#EXTEND " << extensionCriterium << endl;
            cerr << "#LISMEMS " << cluster.matches.size() << endl;
            cerr << "#LISSCORE " << cluster.score << endl;
        }
        if(extensionCriterium){
            //sort matches in this interval according to query position
            const vector<match_t> & matchVector = cluster.matches;
            long chrStart, chrEnd;
            sa.getChromBounds(matchVector[0].ref, chrStart, chrEnd);
            int returnedAlnScore;
            if(alnOptions.print >=5) cerr << "try cluster " << lisIndex << " with " << matchVector.size() << " seeds on " << (cluster.fw ? "forward" : "reverse complemented") << " strand." << endl;
            if(pairedOrder[lisIndex].firstMate)
                returnedAlnScore = extendAlignmentChain(sa.S, cluster.fw ? mate1.sequence : mate1.rcSequence, matchVector, editDist1, minScore1, minCoverage1, alnOptions, chrStart, chrEnd, mate1, cluster.fw);
            else
                returnedAlnScore = extendAlignmentChain(sa.S, cluster.fw ? mate2.sequence : mate2.rcSequence, matchVector, editDist2, minScore2, minCoverage2, alnOptions, chrStart, chrEnd, mate2, cluster.fw);
            //deal with new alignments
            if(cluster.fw) trialF++;
            else trialR++;
            if(returnedAlnScore >= 0 && (!alnOptions.maxTrialBest || maxAlnScore < returnedAlnScore)){
                if(cluster.fw) trialF = 0;
                else trialR = 0;
            }
            alnCount = mate1.alignments.size() + mate2.alignments.size();
            maxAlnScore = max(maxAlnScore, returnedAlnScore);
        }
        lisIndex++;
    }  
    kb_destroy(qCov, treeF);
    kb_destroy(qCov, treeR);
}

//Calculate LIS of both mates, align 1, match LIS to other, align 2
void pairedMatch4(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    //FIND CANDIDATE REGIONS FOR BOTH MATES
    vector<lis_t>  lisIntervals1;
    vector<lis_t>  lisIntervals2;
    if(alnOptions.print >= 3) cerr << "calculate LIS of both mates"  << endl;
    calculateCandRegions(sa, mate1, alnOptions, lisIntervals1, alnOptions.errorPercent);
    calculateCandRegions(sa, mate2, alnOptions, lisIntervals2, alnOptions.errorPercent);
    //SORT ALL LIS
    vector<alnScoreInfo_t> pairedOrder;
    pairedOrder.reserve(lisIntervals1.size()+lisIntervals2.size());
    for(size_t i = 0; i < lisIntervals1.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(true, i,lisIntervals1[i].score));
    for(size_t i = 0; i < lisIntervals2.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(false, i,lisIntervals2[i].score));
    sort(pairedOrder.begin(),pairedOrder.end(),compAlnScorePaired);
    //TRY TO EXTEND LIS
    int Plength1 = mate1.sequence.length();
    int editDist1 = (int)(alnOptions.errorPercent*Plength1)+1;
    int minScore1 = (int) (Plength1*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist1);
    int Plength2 = mate2.sequence.length();
    int editDist2 = (int)(alnOptions.errorPercent*Plength2)+1;
    int minScore2 = (int) (Plength2*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist2);
    if(alnOptions.print >= 3) cerr << "process LIS in sorted order of both paires"  << endl;
    processLisIntervalsPaired(sa, mate1, mate2, alnOptions, lisIntervals1, lisIntervals2, 
    pairedOrder, alnOptions.maxAlnCount, alnOptions.minCoverage, minScore1, minScore2);
    if(mate1.alignments.empty() && !mate2.alignments.empty() && alnOptions.rescue){
        if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
        minScore1 = max(0,(int)(Plength1*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength1));
        int maxAlignmentCount = 20;
        int minCoverage = Plength1*0.01;
        processLisIntervals(sa, mate1, alnOptions, lisIntervals1, editDist1, maxAlignmentCount, minCoverage, minScore1);
    }
    else if(mate2.alignments.empty() && !mate1.alignments.empty() && alnOptions.rescue){
        if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
        minScore2 = max(0,(int)(Plength2*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength2));
        int maxAlignmentCount = 20;
        int minCoverage = Plength1*0.01;
        processLisIntervals(sa, mate2, alnOptions, lisIntervals2, editDist2, maxAlignmentCount, minCoverage, minScore2);
    }
    if(alnOptions.print >= 3) cerr << "first  mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
    if(mate1.alignments.empty() && mate2.alignments.empty() && pairedOpt.pairedRescue){
        minScore1 = max(0,(int)(Plength1*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength1));
        minScore2 = max(0,(int)(Plength2*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength2));
        processLisIntervalsPaired(sa, mate1, mate2, alnOptions, lisIntervals1, lisIntervals2, 
    pairedOrder, 20, 0.01, minScore1, minScore2);
    }
    if(alnOptions.print >= 3) cerr << "first mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
    int alnCount1 = mate1.alignmentCount();
    int alnCount2 = mate2.alignmentCount();
    int alnCount = 0;
    int maxPairedAlignments = alnOptions.alignmentCount+1;
    //ORDER ALIGNMENTS ACCORDING TO SCORE
    pairedOrder.clear();
    pairedOrder.reserve(alnCount1+alnCount2);
    for(size_t i = 0; i < mate1.alignments.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(true, i,mate1.alignments[i]->alignmentScore ));
    for(size_t i = 0; i < mate2.alignments.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(false, i,mate2.alignments[i]->alignmentScore ));
    sort(pairedOrder.begin(),pairedOrder.end(),compAlnScorePaired);
    //TRY AND PAIR THE ALIGNMENTS
    if(alnCount1>0 && alnCount2>0){
        if(alnOptions.print >= 3) cerr << "search concordant alignments" << endl;
        pairAlignmentsOrdered(sa, mate1, mate2, alnOptions, pairedOpt, 
        pairedOrder, alnCount, maxPairedAlignments);
    }
    //RESCUE PROCEDURE IF NO CONCORDANT ALIGNMENTS WERE FOUND
    if(pairedOpt.pairedRescue && alnCount < maxPairedAlignments && max(alnCount1,alnCount2)>0){
        if(alnOptions.print >= 2) cerr << "#PAIREDRESCUE" << endl;
        vector<alnScoreInfo_t> newAlnFromWindow;
        pairAlignmentToWindow(sa, mate1, mate2, alnOptions, pairedOpt, pairedOrder, 
        newAlnFromWindow, alnCount, maxPairedAlignments);
    }
    //FIND DISCORDANT ALIGNMENTS AMONGST HIGHEST SCORING UNPAIRED ALIGNMENTS
    int pairedAlignments = mate1.pairedAlignmentCount;
    if(alnOptions.print >= 2) cerr << "#CONCORDANTCOUNT " << pairedAlignments << endl;
    if(pairedOpt.discordant && alnCount < maxPairedAlignments){
        if(alnOptions.print >= 3) cerr << "search discordant alignments" << endl;
        calculateDiscordantAlignments(sa, mate1, mate2);
    }
    if(alnOptions.print >= 2) cerr << "#DISCORDANTCOUNT " << mate1.pairedAlignmentCount-pairedAlignments << endl;
}

inline int bestLISInWindow(vector<lis_t> & lisIntervalsWindow, int editDist, long refB, long refE, bool fw){
    int bestAln = -1;
    int bestAlnScore = -10000;//set to INT_MIN
    for(size_t i = 0; i < lisIntervalsWindow.size(); i++){
        if(lisIntervalsWindow[i].refBegin+editDist >= refB && 
                lisIntervalsWindow[i].refEnd-editDist <= refE && 
                lisIntervalsWindow[i].fw == fw && lisIntervalsWindow[i].score > bestAlnScore){
            bestAln = i;
            bestAlnScore = lisIntervalsWindow[i].score;
        }
    }
    return bestAln;
}

inline void findOrAlignMateToAlignment(const sparseSA& sa, read_t & alignmentMate, read_t & windowMate, 
        bool alnMateFirst, const align_opt & alnOptions, const paired_opt & pairedOpt, 
        vector<lis_t> & lisIntervalsWindow, int minScore, vector<alnScoreInfo_t> & windowAlns,
        int alignmentIndex){
    alignment_t * aln = alignmentMate.alignments[alignmentIndex];
    //FIND WINDOW
    long refB, refE;
    bool fw;
    int Plength = windowMate.sequence.length();
    dpWindow(aln->refB, aln->refE, alnMateFirst, !aln->flag.test(4), Plength, pairedOpt, refB, refE, fw);
    //CHECK IF AN ALIGNMENT ALREADY EXISTS IN THIS WINDOW
    if(alnOptions.print >= 4) cerr << "find window for aligment of mate " << (alnMateFirst ? "1" : "2") << " on strand " << (aln->flag.test(4) ? "R" : "F") << " and position " << aln->refB << endl;
    int windowMatePos = bestAlignmentInWindow(windowMate, windowAlns, alnMateFirst, refB, refE, fw);
    if(alnOptions.print >= 2) cerr << "#WINDOWLIST " << (windowMatePos>=0) << endl;
    int windowCount = 1;
    if(windowMatePos < 0){
        if(alnOptions.print >= 4) cerr << "found no existing concordant alignment" << endl;
        //TRY TO FIND LIS THAT IS CONCORDANT AND EXTEND THAT
        int editDist = (int)(alnOptions.errorPercent*Plength)+1;
        int lisPos = bestLISInWindow(lisIntervalsWindow, editDist, refB, refE, fw);
        if(lisPos >= 0){
            if(alnOptions.print >= 4) cerr << "found existing concordant LIS" << endl;
            //TRY TO EXTEND THIS LIS
            const vector<match_t> & matchVector = lisIntervalsWindow[lisPos].matches;
            long chrStart, chrEnd;
            sa.getChromBounds(matchVector[0].ref, chrStart, chrEnd);
            size_t extendCount = windowMate.alignments.size();
            extendAlignmentChain(sa.S, lisIntervalsWindow[lisPos].fw ? windowMate.sequence : windowMate.rcSequence, 
                    matchVector, editDist, minScore, Plength*0.01, alnOptions, chrStart, chrEnd, windowMate, lisIntervalsWindow[lisPos].fw);
            if(windowMate.alignments.size() > extendCount){
                lisIntervalsWindow[lisPos].extended = true;
                windowMatePos = extendCount;
                windowCount = windowMate.alignments.size() - extendCount;
            }
        }
        if(alnOptions.print >= 2) cerr << "#WINDOWALIGN " << (windowMatePos >= 0) << endl;
    }
    if(windowMatePos < 0 && pairedOpt.pairedRescue){
        if(alnOptions.print >= 4) cerr << "found no existing concordant LIS or alignment" << endl;
        //TRY TO ALIGN IN WINDOW
        int tLen = refE-refB+1;
        uint8_t * target = (uint8_t *)malloc(tLen);
        for (int j = 0; j < tLen; ++j) // convert to the nt4 encoding
            target[j] = Utils::ORDVALUE[(size_t)sa.S[refB+j]];
        int xtra = KSW_XSTART | (Plength * alnOptions.match < 250? KSW_XBYTE : 0);
        kswr_t alnResult = ksw_align(Plength, fw ? windowMate.readP : windowMate.rcReadP, tLen, 
                target, 5, alnOptions.mat, alnOptions.openGap, alnOptions.extendGap, xtra, 0);
        free(target);
        if(alnResult.score > minScore){
            alignment_t * newAln = new alignment_t();
            newAln->alignmentScore = alnResult.score;
            newAln->qB = alnResult.qb;
            newAln->qE = alnResult.qe;
            newAln->refB = refB + alnResult.tb;
            newAln->refE = refB + alnResult.te;
            newAln->refLength = alnResult.te - alnResult.tb + 1;
            newAln->flag.set(4, !fw);
            windowMate.alignments.push_back(newAln);
            windowMatePos = windowMate.alignments.size()-1;
        }
    }
    if(windowMatePos >= 0){
        if(alnOptions.print >= 4) cerr << "found concordant alignment" << endl;
        for(int i = 0; i < windowCount; i++){
            windowAlns.push_back(alnScoreInfo_t(!alnMateFirst,windowMatePos+i,windowMate.alignments[windowMatePos+i]->alignmentScore));
            //PAIR ALIGNMENTS, SUCCES
            if(alnMateFirst)
                setPaired(sa, alignmentMate, windowMate, alignmentIndex, windowMatePos+i, true);
            else
                setPaired(sa, windowMate, alignmentMate, windowMatePos+i, alignmentIndex, true);
        }
    }
}

inline void processLisIntervalsAndPair(const sparseSA& sa, read_t & mate1, read_t & mate2, 
        const align_opt & alnOptions, const paired_opt & pairedOpt, vector<lis_t> & lisIntervals1, 
        vector<lis_t> & lisIntervals2, const vector<alnScoreInfo_t> & pairedOrder, 
        int maxAlnCount, double minCoveragePercent, int minScore1, int minScore2){
    int minCoverage1 = mate1.sequence.length()*minCoveragePercent;
    int minCoverage2 = mate2.sequence.length()*minCoveragePercent;
    int editDist1 = (int)(alnOptions.errorPercent*mate1.sequence.length())+1;
    int editDist2 = (int)(alnOptions.errorPercent*mate2.sequence.length())+1;
    //for every interval, try to align
    size_t lisIndex = 0;
    //trial parameter for performance reasons
    int trialF = 0;
    int trialR = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    //tree containing information for uniqueness option
    kbtree_t(qCov) *treeF;
    kbtree_t(qCov) *treeR;
    treeF = kb_init(qCov, KB_DEFAULT_SIZE);
    treeR = kb_init(qCov, KB_DEFAULT_SIZE);
    int maxAlnScore = 0;
    int alnCount = 0;
    vector<alnScoreInfo_t>  mate1Alns;
    vector<alnScoreInfo_t>  mate2Alns;
    while(alnCount < maxAlnCount && lisIndex < pairedOrder.size()){
        lis_t & cluster = pairedOrder[lisIndex].firstMate ? lisIntervals1[pairedOrder[lisIndex].position] : lisIntervals2[pairedOrder[lisIndex].position];
        bool unique = hasUnique(cluster, cluster.fw ? treeF : treeR);
        bool extend = unique;
        int trial = cluster.fw ? trialF : trialR;
        bool highCoverage = pairedOrder[lisIndex].firstMate ? cluster.score >= minCoverage1 : cluster.score >= minCoverage2;
        bool extensionCriterium = (trial < alnOptions.maxTrial || extend) && (cluster.matches.size() > 1 || highCoverage);        
        if(alnOptions.print >=2 ){
            cerr << "#LISI " << lisIndex << endl;
            cerr << "#UNIQUE " << unique << endl;
            cerr << "#TRIAL " << trial << endl;
            cerr << "#EXTEND " << (extensionCriterium && !cluster.extended) << endl;
            cerr << "#LISMEMS " << cluster.matches.size() << endl;
            cerr << "#LISSCORE " << cluster.score << endl;
        }
        if(extensionCriterium && !cluster.extended){
            //sort matches in this interval according to query position
            const vector<match_t> & matchVector = cluster.matches;
            long chrStart, chrEnd;
            sa.getChromBounds(matchVector[0].ref, chrStart, chrEnd);
            int returnedAlnScore;
            size_t extendCount;
            if(pairedOrder[lisIndex].firstMate){
                if(alnOptions.print >= 5) cerr << "extend cluster on first mate" << endl;
                extendCount = mate1.alignments.size();
                returnedAlnScore = extendAlignmentChain(sa.S, cluster.fw ? mate1.sequence : mate1.rcSequence, matchVector, editDist1, minScore1, minCoverage1, alnOptions, chrStart, chrEnd, mate1, cluster.fw);
                //NEW PART FOR THIS METHOD: TRY TO FIND A MATE FOR THE ALIGNMENTS THAT WERE FOUND HERE
                if(alnOptions.print >= 5) cerr << "alignment of lis in mate 1 resulted in " << (mate1.alignments.size()-extendCount) << " alignments" << endl;
                size_t lastAln = mate1.alignments.size();
                for(size_t i = extendCount; i < lastAln; i++){
                    mate1Alns.push_back(alnScoreInfo_t(true, extendCount, mate1.alignments[extendCount]->alignmentScore));
                    findOrAlignMateToAlignment(sa, mate1, mate2, true, alnOptions, 
                            pairedOpt, lisIntervals2, minScore2, mate2Alns, extendCount);
                }
            }
            else{
                if(alnOptions.print >= 5) cerr << "extend cluster on second mate" << endl;
                extendCount = mate2.alignments.size();
                returnedAlnScore = extendAlignmentChain(sa.S, cluster.fw ? mate2.sequence : mate2.rcSequence, matchVector, editDist2, minScore2, minCoverage2, alnOptions, chrStart, chrEnd, mate2, cluster.fw);
                if(alnOptions.print >= 5) cerr << "alignment of lis in mate 2 resulted in " << (mate1.alignments.size()-extendCount) << " alignments" << endl;
                size_t lastAln = mate2.alignments.size();
                for(size_t i = extendCount; i < lastAln; i++){
                    mate1Alns.push_back(alnScoreInfo_t(false, extendCount, mate2.alignments[extendCount]->alignmentScore));
                    findOrAlignMateToAlignment(sa, mate2, mate1, false, alnOptions, 
                            pairedOpt, lisIntervals1, minScore1, mate1Alns, extendCount);
                }
            }
            //deal with new alignments
            cluster.extended = true;
            if(cluster.fw) trialF++;
            else trialR++;
            if(returnedAlnScore >= 0 && (!alnOptions.maxTrialBest || maxAlnScore < returnedAlnScore)){
                if(cluster.fw) trialF = 0;
                else trialR = 0;
            }
            alnCount = mate1.alignments.size() + mate2.alignments.size();
            maxAlnScore = max(maxAlnScore, returnedAlnScore);
        }
        lisIndex++;
    }  
    kb_destroy(qCov, treeF);
    kb_destroy(qCov, treeR);
}

//Calculate LIS of both mates, align 1, match LIS to other, align 2
void pairedMatch5(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    //FIND CANDIDATE REGIONS FOR BOTH MATES
    vector<lis_t>  lisIntervals1;
    vector<lis_t>  lisIntervals2;
    if(alnOptions.print >= 3) cerr << "calculate LIS of both mates"  << endl;
    calculateCandRegions(sa, mate1, alnOptions, lisIntervals1, alnOptions.errorPercent);
    calculateCandRegions(sa, mate2, alnOptions, lisIntervals2, alnOptions.errorPercent);
    //SORT ALL LIS
    if(alnOptions.print >= 3) cerr << "first mate resulted in " << lisIntervals1.size() << " intervals" << endl;
    if(alnOptions.print >= 3) cerr << "second mate resulted in " << lisIntervals2.size() << " intervals" << endl;
    vector<alnScoreInfo_t> pairedOrder;
    pairedOrder.reserve(lisIntervals1.size()+lisIntervals2.size());
    for(size_t i = 0; i < lisIntervals1.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(true, i,lisIntervals1[i].score));
    for(size_t i = 0; i < lisIntervals2.size(); i++)
        pairedOrder.push_back(alnScoreInfo_t(false, i,lisIntervals2[i].score));
    sort(pairedOrder.begin(),pairedOrder.end(),compAlnScorePaired);
    //TRY TO EXTEND LIS
    int Plength1 = mate1.sequence.length();
    int editDist1 = (int)(alnOptions.errorPercent*Plength1)+1;
    int minScore1 = (int) (Plength1*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist1);
    int Plength2 = mate2.sequence.length();
    int editDist2 = (int)(alnOptions.errorPercent*Plength2)+1;
    int minScore2 = (int) (Plength2*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*editDist2);
    processLisIntervalsAndPair(sa, mate1, mate2, alnOptions, pairedOpt, lisIntervals1, lisIntervals2, 
    pairedOrder, alnOptions.maxAlnCount, alnOptions.minCoverage, minScore1, minScore2);
    if(mate1.alignments.empty() && !mate2.alignments.empty() && alnOptions.rescue){
        if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
        minScore1 = max(0,(int)(Plength1*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength1));
        int maxAlignmentCount = 20;
        int minCoverage = Plength1*0.01;
        processLisIntervals(sa, mate1, alnOptions, lisIntervals1, editDist1, maxAlignmentCount, minCoverage, minScore1);
    }
    else if(mate2.alignments.empty() && !mate1.alignments.empty() && alnOptions.rescue){
        if(alnOptions.print >= 2) cerr << "#RESCUE" << endl;
        minScore2 = max(0,(int)(Plength2*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength2));
        int maxAlignmentCount = 20;
        int minCoverage = Plength1*0.01;
        processLisIntervals(sa, mate2, alnOptions, lisIntervals2, editDist2, maxAlignmentCount, minCoverage, minScore2);
    }
    if(alnOptions.print >= 3) cerr << "first  mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
    if(mate1.alignments.empty() && mate2.alignments.empty() && pairedOpt.pairedRescue){
        minScore1 = max(0,(int)(Plength1*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength1));
        minScore2 = max(0,(int)(Plength2*alnOptions.match + (alnOptions.mismatch-alnOptions.match)*0.2*Plength2));
        processLisIntervalsAndPair(sa, mate1, mate2, alnOptions, pairedOpt, lisIntervals1, lisIntervals2, 
    pairedOrder, 20, 0.01, minScore1, minScore2);
    }
    if(alnOptions.print >= 3) cerr << "first mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
    //FIND DISCORDANT ALIGNMENTS AMONGST HIGHEST SCORING UNPAIRED ALIGNMENTS
    int maxPairedAlignments = alnOptions.alignmentCount+1;
    int pairedAlignments = mate1.pairedAlignmentCount;
    if(alnOptions.print >= 2) cerr << "#CONCORDANTCOUNT " << pairedAlignments << endl;
    int alnCount = max(mate1.pairedAlignmentCount, mate2.pairedAlignmentCount);
    if(pairedOpt.discordant && alnCount < maxPairedAlignments)
        calculateDiscordantAlignments(sa, mate1, mate2);
    if(alnOptions.print >= 2) cerr << "#DISCORDANTCOUNT " << mate1.pairedAlignmentCount-pairedAlignments << endl;
}

//Calculate Mates 1 at a time, without calculating seeds for the other
//Optimalization: use 'Bailing' method: to enable/disable the calculation for the other mate
void pairedMatch2(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    //ASSYMETRIC ALIGNMENT STRATEGY: ALIGN 1 MATE + RESCUE
    if(alnOptions.print >= 3) cerr << "matching first mate ..." << endl;
    unpairedMatch(sa, mate1, alnOptions);
    if(alnOptions.print >= 3) cerr << "first mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignments);
    int alnCount = 0;
    int maxPairedAlignments = alnOptions.alignmentCount+1;
    //TRAVERSE ALIGNMENTS FOUND FOR MATE 1 AND FIND CONCORDANT ALIGNMENTS IN WINDOW
    vector<alnScoreInfo_t> mate1Alns;
    for(int j = 0; j < (int) mate1.alignments.size(); j++)
        mate1Alns.push_back(alnScoreInfo_t(true, j, mate1.alignments[j]->alignmentScore));
    sort(mate1Alns.begin(),mate1Alns.end(), compAlnScorePaired);
    vector<alnScoreInfo_t> mate2Alns;
    pairAlignmentToWindow(sa, mate1, mate2, alnOptions, pairedOpt, mate1Alns, mate2Alns, alnCount, maxPairedAlignments);
    if(alnOptions.print >= 3) cerr << "part one has resulted in " << alnCount << " concordant alignments" << endl;
    //OTHER DIRECTION STARTING FROM MATE 2
    if(alnCount < maxPairedAlignments){
        //FILTERED UNPAIRED MATCHING OF MATE 2
        if(alnOptions.print >= 3) cerr << "matching second mate ..." << endl;
        unpairedMatch(sa, mate2, alnOptions);
        if(alnOptions.print >= 3) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
        //TRAVERSE ALIGNMENTS FOUND FOR MATE 2 AND FIND CONCORDANT ALIGNMENTS IN WINDOW
        mate2Alns.clear();
        for(int j = 0; j < (int) mate2.alignments.size(); j++)
            if(!mate2.alignments[j]->paired())
                mate2Alns.push_back(alnScoreInfo_t(false, j, mate2.alignments[j]->alignmentScore));
        sort(mate2Alns.begin(),mate2Alns.end(), compAlnScorePaired);
        pairAlignmentToWindow(sa, mate1, mate2, alnOptions, pairedOpt, mate2Alns, mate1Alns, alnCount, maxPairedAlignments);
        if(alnOptions.print >= 3) cerr << "part two has resulted in " << alnCount << " concordant alignments" << endl;
    }
    //RESCUE: PAIR UNPAIRED ALIGNMENTS
    if(pairedOpt.pairedRescue && alnCount < maxPairedAlignments && 
            max((int)mate1.alignments.size(),(int)mate2.alignments.size())>0){
        if(alnOptions.print >= 2) cerr << "#PAIREDRESCUE" << endl;
        vector<alnScoreInfo_t> pairedOrder;
        for(size_t i = 0; i < mate1.alignments.size(); i++)
            pairedOrder.push_back(alnScoreInfo_t(true, i,mate1.alignments[i]->alignmentScore ));
        for(size_t i = 0; i < mate2.alignments.size(); i++)
            pairedOrder.push_back(alnScoreInfo_t(false, i,mate2.alignments[i]->alignmentScore ));
        sort(pairedOrder.begin(),pairedOrder.end(),compAlnScorePaired);
        pairAlignmentsOrdered(sa, mate1, mate2, alnOptions, pairedOpt, pairedOrder, alnCount, maxPairedAlignments);
    }
    //INCLUDE DISCORDANT ALIGNMENTS
    int pairedAlignments = mate1.pairedAlignmentCount;
    if(alnOptions.print >= 2) cerr << "#CONCORDANTCOUNT " << pairedAlignments << endl;
    if(pairedOpt.discordant && alnCount < maxPairedAlignments)
        calculateDiscordantAlignments(sa, mate1, mate2);
    if(alnOptions.print >= 2) cerr << "#DISCORDANTCOUNT " << mate1.pairedAlignmentCount-pairedAlignments << endl;
}

//MODE 3 AND 4 NOT WORKING: [3] crash, [4] not pairing!!!
void pairedMatch(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    if(pairedOpt.mode==1){
        pairedMatch1(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==2){
        pairedMatch2(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==3){
        pairedMatch3(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==4){
        pairedMatch4(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==5){
        pairedMatch5(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==6){
        pairedMatch6(sa, mate1, mate2, alnOptions, pairedOpt);
    }
    else{
        pairedMatch1(sa, mate1, mate2, alnOptions, pairedOpt);
    }      
}

