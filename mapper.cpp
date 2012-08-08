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

#include <algorithm>

#include "mapper.h"

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

bool compAlignmentOffset(const alignment_t i,const alignment_t j){
    return (i.pos < j.pos);
}

bool compAlignmentScore(const alignment_t i,const alignment_t j){
    return (i.alignmentScore > j.alignmentScore);
}

bool compPairScore(const alignment_t i,const alignment_t j){
    return (i.flag.test(0) > j.flag.test(0) || ( i.flag.test(0) == j.flag.test(0) && i.flag.test(1) > j.flag.test(1)) ||
    (i.flag.test(0) == j.flag.test(0) && i.flag.test(1) == j.flag.test(1) && i.alignmentScore > j.alignmentScore));
}

void postProcess(vector<match_t> &matches){
    sort(matches.begin(),matches.end(), compMatches);
}

void calculateSeeds(const sparseSA& sa, const string& P, int min_len, int alignmentCount, vector<match_t>& matches, bool tryHarder){
    sa.SMAM(P, matches, min_len, alignmentCount, false);
    if(tryHarder){//TODO: change try-harder to recalculate only after forward + reverse has been tried
        //easy solution: try reverse and if found: skip
        if(matches.empty())
            sa.SMAM(P, matches, min_len, 1000, false);
        if(matches.empty())
            sa.SMAM(P, matches, 20, 1000, false);
    }
    postProcess(matches);
}

void calculateLISintervals(vector<match_t>& matches, bool fw, int qLength, int editDist, vector<lis_t>& lisIntervals){
    if(matches.size() > 0){
        int begin = 0;
        int end = 0;
        int len = 0;
        while(begin < matches.size()){
            int refEnd = matches[begin].ref - matches[begin].query + qLength + editDist;
            while(end < matches.size() && matches[end].ref + matches[end].len <= refEnd){
                len += matches[end].len;
                end++;
            }
            lisIntervals.push_back(lis_t(&matches,begin,end-1,len,fw));
            begin = end;
            len = 0;
        }
    }
}

void calculateLISintervalsFair(vector<match_t>& matches, bool fw, int qLength, int editDist, vector<lis_t>& lisIntervals){
    int begin = 0;
    int end = 1;
    int len = 0;
    while(begin < matches.size()){
        len = matches[begin].len;
        int refEnd = matches[begin].ref - matches[begin].query + qLength + editDist;
        while(end < matches.size() && matches[end].ref + matches[end].len <= refEnd){
            int leftb = max(matches[end-1].ref+matches[end-1].len-1,matches[end].ref);
            int rightb = min(matches[end-1].ref+matches[end-1].len-1,matches[end].ref+matches[end].len-1);
            len += max(0,rightb-leftb+1);
            end++;
        }
        lisIntervals.push_back(lis_t(&matches,fw,begin,end-1,len));
        begin = end;
        end++;
    }
}

bool extendAlignment(const string& S, const string& P, alignment_t& alignment, vector<match_t>& matches, int begin, int end, int editDist, const align_opt & alnOptions){
    dp_output output;
    bool clipping = !alnOptions.noClipping;
    int curEditDist = 0;
    ///////////////////
    //FIRST SEED
    //////////////////
    match_t firstSeed = matches[begin];
    int refstrLB = firstSeed.ref;
    int refstrRB = refstrLB + firstSeed.len -1;
    int queryLB = firstSeed.query;
    int queryRB = queryLB + firstSeed.len-1;
    if(firstSeed.query > 0 && firstSeed.ref > 0){
        //Alignment now!!! with beginQ-E in ref to begin match and beginQ to match
        int alignmentBoundLeft = max(refstrLB-queryLB-editDist,0);
        boundaries grenzen(alignmentBoundLeft,refstrLB-1,0,queryLB-1);
        dp_type types;
        types.freeRefB = true;
        types.freeQueryB = clipping;
        dpBandStatic( S, P, grenzen, alnOptions.scores, types, ERRORSTRING, output, editDist-curEditDist, false);
        if(clipping && grenzen.queryB> 0){
            alignment.cigarChars.push_back('S');
            alignment.cigarLengths.push_back(grenzen.queryB);
        }
        alignment.cigarChars.insert(alignment.cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
        alignment.cigarLengths.insert(alignment.cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
        alignment.pos = grenzen.refB+1;
        curEditDist += output.editDist;
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
            } else if (minRefDist <= 0) {
                alignment.cigarChars.push_back('I');
                alignment.cigarChars.push_back('=');
                alignment.cigarLengths.push_back(minQDist);
                alignment.cigarLengths.push_back(match.len - 1 + minRefDist);
                curEditDist += minQDist;
            } else {//both distances are positive and not equal to (1,1)
                boundaries grenzen(refstrRB + 1, refstrRB + minRefDist - 1, queryRB + 1, queryRB + minQDist - 1);
                dp_type types;
                dpBandStatic( S, P, grenzen, alnOptions.scores, types, ERRORSTRING, output, editDist-curEditDist, false);
                alignment.cigarChars.insert(alignment.cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
                alignment.cigarLengths.insert(alignment.cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
                alignment.cigarChars.push_back('=');
                alignment.cigarLengths.push_back(match.len);
                curEditDist += output.editDist;
                output.clear();
            }
            //set new positions
            memsIndex = minDistMem + 1;
            //add matches
            refstrRB = match.ref + match.len - 1;
            queryRB = match.query + match.len - 1;
        }
    }
    /////////////
    //FINAL DP?
    ////////////
    //possible at the end characters
    //mapInterval += Integer.toString(queryRB+1)+"]";
    alignment.refLength = queryRB-alignment.pos+1;
    if (queryRB + 1 < P.length() && refstrRB < S.length() - 1 && curEditDist <= editDist) {
        int refAlignRB = refstrRB + (P.length() - queryRB) + 1 + editDist - curEditDist;
        if (refAlignRB >= S.length())
            refAlignRB = S.length() - 1;
        boundaries grenzen(refstrRB + 1, refAlignRB, queryRB + 1, P.length() - 1);
        dp_type types;
        types.freeRefE = true;
        types.freeQueryE = clipping;
        dpBandStatic( S, P, grenzen, alnOptions.scores, types, ERRORSTRING, output, editDist-curEditDist+1, false);
        if(grenzen.queryE > queryRB){
            int addToLength = grenzen.queryE-queryRB;
            if(output.cigarChars[output.cigarChars.size()-1] == 'I')
                addToLength -= output.cigarLengths[output.cigarLengths.size()-1];
            alignment.refLength += addToLength;
        }
        alignment.cigarChars.insert(alignment.cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
        alignment.cigarLengths.insert(alignment.cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
        if (clipping && grenzen.queryE < P.length() - 1) {
            alignment.cigarChars.push_back('S');
            alignment.cigarLengths.push_back(P.length() - 1 - grenzen.queryE);
        }
        curEditDist += output.editDist;
        output.clear();
    } else if (queryRB + 1 < P.length() && curEditDist < editDist) {
        alignment.cigarChars.push_back(clipping ? 'S' : 'I');
        alignment.cigarLengths.push_back(P.length() - queryRB - 1);
        curEditDist += clipping ? 0 : P.length() - 1 - queryRB;
    }
    //TODO: check for possible optimizations (static initialisations
    //TODO: reorder matches according to best hits
    if(curEditDist <= editDist){
        alignment.editDist = curEditDist;
        return true;
    }
    else
        return false;
}

void inexactMatch(const sparseSA& sa, read_t & read,const align_opt & alnOptions, bool fwStrand, bool print){
    string P = read.sequence;
    int Plength = P.length();
    int editDist = (int)(alnOptions.errorPercent*Plength)+1;
    if(!fwStrand)
        P = read.rcSequence;
    int min_len = alnOptions.minMemLength;
    if(!alnOptions.fixedMinLength && Plength/editDist > min_len)
        min_len = Plength/editDist;
    vector<match_t> matches;
    //calc seeds
    calculateSeeds(sa, P, min_len, alnOptions.alignmentCount/2, matches, alnOptions.tryHarder);
    //sort matches
    if(matches.size()>0){
        /////////////////////////
        //FIND CANDIDATE REGIONS
        /////////////////////////
        vector<lis_t> lisIntervals;
        calculateLISintervals(matches, fwStrand, Plength, editDist, lisIntervals);
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
        while(alnCount < alnOptions.alignmentCount/2 &&
                lisIndex < lisIntervals.size() &&
                lisIntervals[lisIndex].len > (Plength*alnOptions.minCoverage)/100 &&
                trial < alnOptions.maxTrial){
            //sort matches in this interval according to query position
            int begin = lisIntervals[lisIndex].begin;
            int end = lisIntervals[lisIndex].end;
            //sort this candidate region by query position
            sort(matches.begin()+begin,matches.begin()+end+1, compMatchesQuery);
            alignment_t alignment;
            bool extended = extendAlignment(sa.S, P, alignment, matches, begin, end, editDist, alnOptions);
            if(extended){
                if(!fwStrand){
                    alignment.flag.set(4,true);
                    if(alnOptions.unique && !read.alignments.empty() && alignment.editDist < read.alignments[0].editDist)
                        read.alignments.pop_back();
                }
                alnCount++;
                read.alignments.push_back(alignment);
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
    }
}

void unpairedMatch(const sparseSA& sa, read_t & read,const align_opt & alnOptions, bool print){
    string P = read.sequence;
    string Prc = read.rcSequence;
    int Plength = P.length();
    int editDist = (int)(alnOptions.errorPercent*Plength)+1;
    int min_len = alnOptions.minMemLength;
    if(!alnOptions.fixedMinLength && Plength/editDist > min_len)
        min_len = Plength/editDist;
    vector<match_t> matches;
    vector<match_t> matchesRC;
    //calc seeds
    if(!alnOptions.noFW){
        calculateSeeds(sa, P, min_len, alnOptions.alignmentCount, matches, alnOptions.tryHarder);
    }
    if(!alnOptions.noRC){
        calculateSeeds(sa, Prc, min_len, alnOptions.alignmentCount, matchesRC, alnOptions.tryHarder);
    }
    //sort matches
    if(matches.size()>0 || matchesRC.size()>0){
        /////////////////////////
        //FIND CANDIDATE REGIONS
        /////////////////////////
        vector<lis_t> lisIntervals;
        calculateLISintervals(matches, true, Plength, editDist, lisIntervals);
        calculateLISintervals(matchesRC, false, Plength, editDist, lisIntervals);
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
            int begin = lisIntervals[lisIndex].begin;
            int end = lisIntervals[lisIndex].end;
            vector<match_t>* matchVector = lisIntervals[lisIndex].matches;
            //sort this candidate region by query position
            sort(matchVector->begin()+begin,matchVector->begin()+end+1, compMatchesQuery);
            alignment_t alignment;
            bool extended = extendAlignment(sa.S, lisIntervals[lisIndex].fw ? P : Prc, alignment, *matchVector, begin, end, editDist, alnOptions);
            if(extended){
                if(!lisIntervals[lisIndex].fw){
                    alignment.flag.set(4,true);
//                    if(alnOptions.unique && !read.alignments.empty() && alignment.editDist < read.alignments[0].editDist)
//                        read.alignments.pop_back();
                }
                alnCount++;
                read.alignments.push_back(alignment);
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
    }
}

//PAIRED END OPTIONS
bool isConcordant(const alignment_t& mate1, const alignment_t& mate2, const paired_opt& options){
    bool firstLeft = true;
    //Check same chromosome
    if(mate1.rname != mate2.rname)
        return false;
    //Check strand
    if(options.orientation == PAIR_FF){
        if(mate1.flag.test(4)!=mate1.flag.test(4))
            return false;
        else
            firstLeft = !mate1.flag.test(4);//if rc must be right
    }
    else if(options.orientation == PAIR_FR){
        if(mate1.flag.test(4)==mate1.flag.test(4))
            return false;
        else
            firstLeft = !mate1.flag.test(4);//if rc must be right
    }
    else{
        if(mate1.flag.test(4)==mate1.flag.test(4))
            return false;
        else
            firstLeft = mate1.flag.test(4);//if rc must be left
    }
    //check fragment insert size
    long begin = min(mate1.pos,mate2.pos);
    long end = max(mate1.pos+mate1.refLength-1,mate2.pos+mate2.refLength-1);
    long size = end - begin +1;
    if(size < options.minInsert || size > options.maxInsert)
        return false;
    long begin1 = mate1.pos;
    long end1 = begin1 + mate1.refLength-1;
    long begin2 = mate2.pos;
    long end2 = begin2 + mate2.refLength-1;
    //chack contain
    bool contain = (begin1 >= begin2 && end1 <= end2) || (begin2 >= begin1 && end2 <= end1);
    if(contain && !options.contain){
        return false;
    }
    //check overlap
    bool overlap = contain || (begin2 >= begin1 && begin2 <= end1) || (end2 >= begin1 && end2 <= end1);
    if(!options.overlap && overlap){
        return false;
    }
    //check dovetail
    if((firstLeft && (begin2 < begin1)) || (!firstLeft && (begin1 < begin2)) ){
        if(!overlap)
            return false;
        else if(!options.dovetail)
            return false;
    }
    return true;
}
struct pair_t {
  int mate1; // alignment in first mate
  int mate2; // alignment in second mate
  pair_t(int first, int second): mate1(first), mate2(second) {}
};

void setPaired(alignment_t& toset, alignment_t& mate, bool firstMate){
    toset.flag.set(0,true);
    toset.flag.set(1,true);
    toset.flag.set(5,mate.flag.test(4));
    toset.flag.set(6,firstMate);
    toset.flag.set(7,!firstMate);
    toset.pnext = mate.pos;
    if(toset.rname == mate.rname){
        toset.rnext = "=";
        toset.tLength = max(toset.pos+toset.refLength-1,mate.pos+mate.refLength-1) - min(toset.pos,mate.pos) + 1;
        long firstPos = toset.flag.test(4) ? toset.pos+toset.refLength-1 : toset.pos;
        long secondPos = mate.flag.test(4) ? mate.pos+mate.refLength-1 : mate.pos;
        if(firstPos > secondPos) toset.tLength *= -1;
    }
    else{
        toset.rnext = mate.rnext;
        toset.tLength = 0;
    }
}

void setUnPaired(alignment_t& toset){
    toset.flag.set(0,true);
}

//Only works when only requiring concordant pairs OR combi concordant/unpaired
void pairedMatch1(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt, bool print){
    //Calculate the mappings for read 1 and 2
    unpairedMatch(sa, mate1, alnOptions, print);
    unpairedMatch(sa, mate2, alnOptions, print);
    int alnCount1 = mate1.alignmentCount();
    int alnCount2 = mate2.alignmentCount();
    if(alnCount1>0 && alnCount2>0){
        int concordant = 0;
        int discordant = 0;
        vector<pair_t> discordantAln;
        //sort according to alignmentScore
        sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignmentScore);
        sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignmentScore);
        int maxScoreFirst = mate1.alignments[0].alignmentScore;
        int maxScoreSecond = mate2.alignments[0].alignmentScore;
        for(int j = 0; j < alnCount2; j++){
            mate2.alignments[j].setGlobalPos(sa);
        }
        int i = 0;
        int alnCount = 0;
        //search concordant alignments
        while(i < alnCount1 && alnCount < alnOptions.alignmentCount){
            alignment_t & aln1 = mate1.alignments[i];
            aln1.setGlobalPos(sa);
            int j = 0;
            while(j < alnCount2 && alnCount < alnOptions.alignmentCount){
                alignment_t & aln2 = mate2.alignments[j];
                if(isConcordant(aln1, aln2, pairedOpt)){
                    //Set the concordants and add another if necessary
                    if(!aln1.flag.test(0)){
                        setPaired(aln1, aln2, true);
                    }
                    else{
                        alignment_t newAln1(aln1);
                        setPaired(newAln1, aln2, true);
                        mate1.alignments.push_back(newAln1);
                    }
                    if(!aln2.flag.test(0)){
                        setPaired(aln2, aln1, false);
                    }
                    else{
                        alignment_t newAln2(aln2);
                        setPaired(newAln2, aln1, false);
                        mate2.alignments.push_back(newAln2);
                    }
                    concordant++;
                    alnCount++;
                }
                else if(pairedOpt.discordant && maxScoreFirst == aln1.alignmentScore && maxScoreSecond == aln2.alignmentScore){
                    discordantAln.push_back(pair_t(i,j));
                    discordant++;
                }
                j++;
            }
            i++;
        }
        //search discordant alignments
        if(pairedOpt.discordant){
            i = 0;
            while(i < discordantAln.size() && alnCount < alnOptions.alignmentCount){
                alignment_t & aln1 = mate1.alignments[discordantAln[i].mate1];
                alignment_t & aln2 = mate2.alignments[discordantAln[i].mate2];
                if(!aln1.flag.test(0)){
                    setPaired(aln1, aln2, true);
                }
                else{
                    alignment_t newAln1(aln1);
                    setPaired(newAln1, aln2, true);
                    mate1.alignments.push_back(newAln1);
                }
                if(!aln2.flag.test(0)){
                    setPaired(aln2, aln1, false);
                }
                else{
                    alignment_t newAln2(aln2);
                    setPaired(newAln2, aln1, false);
                    mate2.alignments.push_back(newAln2);
                }
                i++;
                alnCount++;
            }
        }
        //sort alignments such that first concordant, than discordant, than unpaired
        if(pairedOpt.mixed){
            i = 0;
            int alnCountFirst = alnCount;
            while(i < alnCount1 && alnCountFirst < alnOptions.alignmentCount){
                if(!mate1.alignments[i].flag.test(0))
                    setUnPaired(mate1.alignments[i]);
                i++;
                alnCountFirst++;
            }
            i = 0;
            int alnCountSecond = alnCount;
            while(i < alnCount2 && alnCountSecond < alnOptions.alignmentCount){
                if(!mate2.alignments[i].flag.test(0))
                    setUnPaired(mate2.alignments[i]);
                i++;
                alnCountSecond++;
            }
        }
//        sort(mate1.alignments.begin(),mate1.alignments.end(), compPairScore);
//        sort(mate1.alignments.begin(),mate1.alignments.end(), compPairScore);
    }
    else if(alnCount1>0 || alnCount2>0){//bail1 + set to unpaired??? flag 3
        //IMPLEMENT BAIL ME: PARTIAL MODE2 ON THE ONE THAT IS NON_ZERO
    }
    else{//bail2
        //IMPLEMENT BAIL ME: FULL MODE2
    }
}

void pairedMatch(const sparseSA& sa, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt, int mode, bool print){
    if(mode==1){
        pairedMatch1(sa, mate1, mate2, alnOptions, pairedOpt, print);
    }
    else{
        //mode unknown
    }
}