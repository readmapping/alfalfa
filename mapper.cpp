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
#include <set>

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

bool compIntervalsRef(const lis_t i,const lis_t j){
    return compMatches(i.matches->at(i.begin), j.matches->at(j.begin));
}

bool compAlignmentScore(const alignment_t i,const alignment_t j){
    return (i.alignmentScore > j.alignmentScore || (i.alignmentScore == j.alignmentScore && i.globPos < j.globPos));
}

bool compAlignmentAndScore(const lis_t i,const lis_t j){
    return (i.extended&&!j.extended )|| 
            ((i.extended && j.extended ) && 
            (i.alignment.alignmentScore > j.alignment.alignmentScore 
            || (i.alignment.alignmentScore == j.alignment.alignmentScore && 
            (i.alignment.globPos < j.alignment.globPos))));
}

bool compPairScore(const alignment_t i,const alignment_t j){
    if(i.concordant() > j.concordant()) return true;
    else if(i.concordant() == j.concordant()){
        return (i.flag.test(1) > j.flag.test(1) || ( i.flag.test(1) == j.flag.test(1) && i.flag.test(0) > j.flag.test(0)) ||
    (i.flag.test(0) == j.flag.test(0) && i.flag.test(1) == j.flag.test(1) && i.alignmentScore > j.alignmentScore));
    }
    else{
        return false;
    }
}

void postProcess(vector<match_t> &matches){
    sort(matches.begin(),matches.end(), compMatches);
}

void calculateSeeds(const sparseSA& sa, const string& P, int min_len, int alignmentCount, vector<match_t>& matches, bool tryHarder){
    int maxBranchWidth = max(alignmentCount,20);
    sa.SMAM(P, matches, min_len, maxBranchWidth, false);
    if(tryHarder){//TODO: change try-harder to recalculate only after forward + reverse has been tried
        //easy solution: try reverse and if found: skip
        if(matches.empty())
            sa.SMAM(P, matches, min_len, 1000, false);
        if(matches.empty())
            sa.SMAM(P, matches, 20, 1000, false);
    }
    postProcess(matches);
}

void calculateLISintervals(vector<match_t>& matches, bool fw, long qLength, long editDist, vector<lis_t>& lisIntervals){
    if(matches.size() > 0){
        int begin = 0;
        int end = 0;
        long len = 0;
        while(begin < matches.size()){
            long refEnd = matches[begin].ref - matches[begin].query + qLength + editDist;
            while(end < matches.size() && matches[end].ref + matches[end].len <= refEnd){
                len += matches[end].len;
                end++;
            }
            lisIntervals.push_back(lis_t(&matches,begin,end-1,len,fw));
            begin = end;
            len = 0L;
        }
    }
}

void calculateLISintervalsFair(vector<match_t>& matches, bool fw, long qLength, long editDist, vector<lis_t>& lisIntervals){
    int begin = 0;
    int end = 1;
    long len = 0;
    while(begin < matches.size()){
        len = matches[begin].len;
        long refEnd = matches[begin].ref - matches[begin].query + qLength + editDist;
        while(end < matches.size() && matches[end].ref + matches[end].len <= refEnd){
            long leftb = max(matches[end-1].ref+matches[end-1].len-1,matches[end].ref);
            long rightb = min(matches[end-1].ref+matches[end-1].len-1,matches[end].ref+matches[end].len-1);
            len += max(0L,rightb-leftb+1L);
            end++;
        }
        lisIntervals.push_back(lis_t(&matches,fw,begin,end-1,len));
        begin = end;
        end++;
    }
}

bool extendAlignment(dynProg& dp_, const string& S, const string& P, 
        alignment_t& alignment, vector<match_t>& matches, 
        int begin, int end, long editDist, const align_opt & alnOptions){
    dp_output output;
    bool clipping = !alnOptions.noClipping;
    long curEditDist = 0;
    ///////////////////
    //FIRST SEED
    //////////////////
    match_t firstSeed = matches[begin];
    long refstrLB = firstSeed.ref;
    long refstrRB = refstrLB + firstSeed.len -1L;
    long queryLB = firstSeed.query;
    long queryRB = queryLB + firstSeed.len-1L;
    if(firstSeed.query > 0L && firstSeed.ref > 0L){
        //Alignment now!!! with beginQ-E in ref to begin match and beginQ to match
        long alignmentBoundLeft = max(refstrLB-queryLB-editDist,0L);
        boundaries grenzen(alignmentBoundLeft,refstrLB-1,0,queryLB-1);
        dp_type types;
        types.freeRefB = true;
        types.freeQueryB = clipping;
        dp_.dpBandStatic( S, P, grenzen, types, ERRORSTRING, output, editDist-curEditDist, false);
        queryLB = grenzen.queryB;
        if(output.cigarChars[output.cigarChars[0]] == 'I')
                queryLB += (long)output.cigarLengths[0];
        if(clipping && grenzen.queryB> 0){
            alignment.cigarChars.push_back('S');
            alignment.cigarLengths.push_back(grenzen.queryB);
        }
        alignment.cigarChars.insert(alignment.cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
        alignment.cigarLengths.insert(alignment.cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
        alignment.globPos = grenzen.refB+1L;
        alignment.alignmentScore += output.dpScore;
        curEditDist += output.editDist;
        output.clear();
    }//fill in the starting position in the ref sequence
    else {
        alignment.globPos = firstSeed.ref + 1L;
    }
    if (firstSeed.ref == 0L && firstSeed.query > 0L) {
        alignment.cigarChars.push_back(clipping ? 'S' : 'I');
        alignment.cigarLengths.push_back(firstSeed.query);
        curEditDist += clipping ? 0 : firstSeed.query;
        alignment.alignmentScore += (clipping ? 0 : dp_.scores.openGap + dp_.scores.extendGap*firstSeed.query);
    }
    alignment.cigarChars.push_back('=');
    alignment.cigarLengths.push_back(firstSeed.len);
    alignment.alignmentScore += dp_.scores.match*firstSeed.len;
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
        long minDist = (long)P.length();
        long minQDist = (long)P.length();
        long minRefDist = (long)P.length();
        //temporary distance for this mem
        long dist = minDist;
        while (indexIncrease + memsIndex <= end && dist <= minDist) {
            match_t match = matches[memsIndex + indexIncrease];
            //check distance: k is enlarged and compare minRef-refstrRB
            long qDist = match.query - queryRB;
            long refDist = match.ref - refstrRB;
            if (qDist < 0 || refDist < 0) {
                if (qDist >= refDist) {
                    qDist -= refDist; //qDist insertions
                    if (qDist + queryRB >= (long)P.length() || 1 - refDist >= match.len)//and refDist +1 matches lost
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
                alignment.alignmentScore += dp_.scores.match*(match.len - 1 + minQDist) + dp_.scores.openGap + dp_.scores.extendGap*minRefDist;
                curEditDist += Utils::contains(S, refstrRB + 1L, refstrRB + minRefDist - 1L, '`') ? editDist + 1 : minRefDist; //boundary of ref sequences is passed
                //TODO: only use of contains: add inline here
            } else if (minRefDist <= 0) {
                alignment.cigarChars.push_back('I');
                alignment.cigarChars.push_back('=');
                alignment.cigarLengths.push_back(minQDist);
                alignment.cigarLengths.push_back(match.len - 1 + minRefDist);
                alignment.alignmentScore += dp_.scores.match*(match.len - 1 + minRefDist) + dp_.scores.openGap + dp_.scores.extendGap*minQDist;
                curEditDist += minQDist;
            } else {//both distances are positive and not equal to (1,1)
                boundaries grenzen(refstrRB + 1L, refstrRB + minRefDist - 1L, queryRB + 1, queryRB + minQDist - 1);
                dp_type types;
                dp_.dpBandStatic( S, P, grenzen, types, ERRORSTRING, output, editDist-curEditDist, false);
                alignment.cigarChars.insert(alignment.cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
                alignment.cigarLengths.insert(alignment.cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
                alignment.cigarChars.push_back('=');
                alignment.cigarLengths.push_back(match.len);
                alignment.alignmentScore += dp_.scores.match*match.len + output.dpScore;
                curEditDist += output.editDist;
                output.clear();
            }
            //set new positions
            memsIndex = minDistMem + 1;
            //add matches
            refstrRB = match.ref + match.len - 1L;
            queryRB = match.query + match.len - 1L;
        }
    }
    /////////////
    //FINAL DP?
    ////////////
    //possible at the end characters
    //mapInterval += Integer.toString(queryRB+1)+"]";
    alignment.refLength = queryRB-queryLB;
    if (queryRB + 1 < P.length() && refstrRB < (long)S.length() - 1L && curEditDist <= editDist) {
        long refAlignRB = refstrRB + ((long)P.length() - queryRB) + 1L + editDist - curEditDist;
        if (refAlignRB >= (long)S.length())
            refAlignRB = (long)S.length() - 1L;
        boundaries grenzen(refstrRB + 1L, refAlignRB, queryRB + 1, P.length() - 1);
        dp_type types;
        types.freeRefE = true;
        types.freeQueryE = clipping;
        dp_.dpBandStatic( S, P, grenzen, types, ERRORSTRING, output, editDist-curEditDist+1, false);
        if(grenzen.queryE > queryRB){
            int addToLength = grenzen.queryE-queryRB;
            if(output.cigarChars[output.cigarChars.size()-1] == 'I')
                addToLength -= output.cigarLengths[output.cigarLengths.size()-1];
            alignment.refLength += addToLength;
        }
        alignment.cigarChars.insert(alignment.cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
        alignment.cigarLengths.insert(alignment.cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
        alignment.alignmentScore += output.dpScore;
        if (clipping && grenzen.queryE < P.length() - 1) {
            alignment.cigarChars.push_back('S');
            alignment.cigarLengths.push_back(P.length() - 1 - grenzen.queryE);
        }
        curEditDist += output.editDist;
        output.clear();
    } else if (queryRB + 1 < P.length() && curEditDist < editDist) {
        alignment.cigarChars.push_back(clipping ? 'S' : 'I');
        alignment.cigarLengths.push_back(P.length() - queryRB - 1);
        alignment.alignmentScore += (clipping ? 0 : ((P.length() - 1 - queryRB)*dp_.scores.extendGap + dp_.scores.openGap));
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

void inexactMatch(const sparseSA& sa, dynProg& dp_, read_t & read,const align_opt & alnOptions, bool fwStrand, bool print){
    string P = read.sequence;
    long Plength = (long) P.length();
    long editDist = (long)(alnOptions.errorPercent*Plength)+1;
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
        while(alnCount < max(1,alnOptions.alignmentCount/2) &&
                lisIndex < lisIntervals.size() &&
                lisIntervals[lisIndex].len > (Plength*alnOptions.minCoverage)/100 &&
                trial < alnOptions.maxTrial){
            //sort matches in this interval according to query position
            int begin = lisIntervals[lisIndex].begin;
            int end = lisIntervals[lisIndex].end;
            //sort this candidate region by query position
            sort(matches.begin()+begin,matches.begin()+end+1, compMatchesQuery);
            alignment_t alignment;
            bool extended = extendAlignment(dp_, sa.S, P, alignment, matches, begin, end, editDist, alnOptions);
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

void unpairedMatch(const sparseSA& sa, dynProg& dp_, read_t & read,const align_opt & alnOptions, bool print){
    string P = read.sequence;
    string Prc = read.rcSequence;
    long Plength = (long) P.length();
    long editDist = (long)(alnOptions.errorPercent*Plength)+1;
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
            bool extended = extendAlignment(dp_, sa.S, lisIntervals[lisIndex].fw ? P : Prc, alignment, *matchVector, begin, end, editDist, alnOptions);
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
        if(mate1.flag.test(4)!=mate2.flag.test(4))
            return false;
        else
            firstLeft = !mate1.flag.test(4);//if rc must be right
    }
    else if(options.orientation == PAIR_FR){
        if(mate1.flag.test(4)==mate2.flag.test(4))
            return false;
        else
            firstLeft = !mate1.flag.test(4);//if rc must be right
    }
    else{
        if(mate1.flag.test(4)==mate2.flag.test(4))
            return false;
        else
            firstLeft = mate1.flag.test(4);//if rc must be left
    }
    //check fragment insert size
    long begin = min(mate1.pos,mate2.pos);
    long end = max(mate1.pos+(long)mate1.refLength-1L,mate2.pos+(long)mate2.refLength-1L);
    long size = end - begin +1L;
    if(size < options.minInsert || size > options.maxInsert)
        return false;
    long begin1 = mate1.pos;
    long end1 = begin1 + (long)mate1.refLength-1L;
    long begin2 = mate2.pos;
    long end2 = begin2 + (long)mate2.refLength-1L;
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
  pair_t(): mate1(0), mate2(0) {}
  pair_t(int first, int second): mate1(first), mate2(second) {}
};

void setUnPaired(alignment_t& toset, read_t& read, bool upstream){
    toset.flag.set(0,true);
    toset.flag.set(6,upstream);
    toset.flag.set(7,!upstream);
    toset.flag.set(3,true);
    read.pairedAlignmentCount++;
}

void setPaired(alignment_t& mate1, alignment_t& mate2, read_t& upstream, read_t& downstream, bool concordant){
    mate1.flag.set(0, true);
    mate2.flag.set(0, true);
    mate1.addMate(mate2, concordant, true);
    mate2.addMate(mate1, concordant, false);
    upstream.pairedAlignmentCount++;
    downstream.pairedAlignmentCount++;
}

//aln both and pair
void pairedMatch1(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt, bool print){
    //Calculate the mappings for read 1 and 2
    unpairedMatch(sa, dp_, mate1, alnOptions, print);
    unpairedMatch(sa, dp_, mate2, alnOptions, print);
    int alnCount1 = mate1.alignmentCount();
    int alnCount2 = mate2.alignmentCount();
    int i = 0;
    int alnCount = 0;
    if(alnCount1>0 && alnCount2>0){
        int concordant = 0;
        int discordant = 0;
        vector<pair_t> discordantAln;
        //sort according to alignmentScore
        sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignmentScore);
        sort(mate2.alignments.begin(),mate2.alignments.end(), compAlignmentScore);
        int maxScoreFirst = mate1.alignments[0].alignmentScore;
        int maxScoreSecond = mate2.alignments[0].alignmentScore;
        for(int j = 0; j < alnCount2; j++){
            mate2.alignments[j].setLocalPos(sa);
        }
        //search concordant alignments
        while(i < alnCount1 && alnCount < alnOptions.alignmentCount){
            alignment_t & aln1 = mate1.alignments[i];
            aln1.setLocalPos(sa);
            int j = 0;
            while(j < alnCount2 && alnCount < alnOptions.alignmentCount){
                alignment_t & aln2 = mate2.alignments[j];
                if(isConcordant(aln1, aln2, pairedOpt)){
                    //Set the concordants and add another if necessary
                    setPaired(aln1,aln2,mate1,mate2,true);
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
                setPaired(aln1,aln2,mate1,mate2,false);
                i++;
                alnCount++;
            }
        }
    }
    //sort alignments such that first concordant, than discordant, than unpaired
    if(pairedOpt.mixed){
        i = 0;
        int alnCountFirst = alnCount;
        while(i < alnCount1 && alnCountFirst < alnOptions.alignmentCount){
            if(!mate1.alignments[i].flag.test(0)){
                setUnPaired(mate1.alignments[i],mate1, true);
                alnCountFirst++;
            }
            i++;
        }
        i = 0;
        int alnCountSecond = alnCount;
        while(i < alnCount2 && alnCountSecond < alnOptions.alignmentCount){
            if(!mate2.alignments[i].flag.test(0)){
                setUnPaired(mate2.alignments[i],mate2, false);
                alnCountSecond++;
            }
            i++;
        }
        mate2.pairedAlignmentCount = alnCountSecond;
    }
    else if(alnCount1>0 || alnCount2>0){//bail1 + set to unpaired??? flag 3
        //IMPLEMENT BAIL ME: PARTIAL MODE2 ON THE ONE THAT IS NON_ZERO
    }
    else{//bail2
        //IMPLEMENT BAIL ME: FULL MODE2
    }
}

bool alignFromLIS(const sparseSA& sa, dynProg& dp_, read_t& read, lis_t & lis, long editDist, const align_opt & alnOptions){
    if(lis.extended)
        return true;
    int begin = lis.begin;
    int end = lis.end;
    vector<match_t>* matchVector = lis.matches;
    //sort this candidate region by query position
    sort(matchVector->begin()+begin,matchVector->begin()+end+1, compMatchesQuery);
    lis.extended = extendAlignment(dp_, sa.S, lis.fw ? read.sequence : read.rcSequence, lis.alignment, *matchVector, begin, end, editDist, alnOptions);
    if(lis.extended){
        if(!lis.fw){
            lis.alignment.flag.set(4,true);
        }
        lis.alignment.setLocalPos(sa);
    }
    return lis.extended;
}

bool isPosConcordant(const lis_t& mate1, const lis_t& mate2, long editDist1, long editDist2, long M1length, long M2length, const paired_opt& options){
    bool firstLeft = true;
    bool mate1FW = mate1.fw;
    bool mate2FW = mate2.fw;
    //Check strand
    if(options.orientation == PAIR_FF || options.orientation == PAIR_FR){
        firstLeft = mate1FW;//if rc must be right
    }
    else{
        firstLeft = !mate1FW;//if rc must be left
    }
    //find possible offset
    long posLeft1 = mate1.matches->at(mate1.begin).ref-mate1.matches->at(mate1.begin).query;
    long posLeft2 = mate2.matches->at(mate2.begin).ref-mate2.matches->at(mate2.begin).query;
    long posRight1 = posLeft1 + M1length;
    long posRight2 = posLeft2 + M2length;
    //check fragment insert size
    long minsize = max(posRight2+editDist2,posRight1+editDist1) - min(posLeft2+editDist2,posLeft1+editDist1) +1;
    if(minsize < options.minInsert)
        return false;
    long maxsize = max(posRight2-editDist2,posRight1-editDist1) - min(posLeft2-editDist2,posLeft1-editDist1) +1;
    if(maxsize > options.maxInsert)
        return false;
    //check flow direction is correct
    long anchor1 = mate1FW ? posLeft1 : posRight1;
    long anchor2 = mate2FW ? posLeft2 : posRight2;
    if(firstLeft){
        if(anchor1-editDist1 > anchor2 + editDist2) return false;
    }
    else{
        if(anchor2 - editDist2 > anchor1 + editDist1) return false;
    }
    //check contain
    bool contain = ((posLeft1-editDist1 >= posLeft2+editDist2 && posRight1+editDist1 <= posRight2-editDist2) 
    || (posLeft2-editDist2 >= posLeft1+editDist1 && posRight2+editDist2 <= posRight1-editDist1));
    if(contain && !options.contain){
        return false;
    }
    //check overlap
    bool overlap = contain || !(posRight2-editDist2 < posLeft1+editDist1 || posLeft2+editDist2 > posRight1-editDist1);
    if(!options.overlap && overlap){
        return false;
    }
    //check dovetail
    if((firstLeft && (posLeft2+editDist2 < posLeft1-editDist1)) || (!firstLeft && (posLeft1+editDist1 < posLeft2-editDist2)) ){
        if(!overlap)
            return false;
        else if(!options.dovetail)
            return false;
    }
    return true;
}

void coupleLIS(const sparseSA& sa, 
        dynProg & dp_,
        read_t & mate1, 
        read_t & mate2,
        vector<lis_t>& lisIntervalsFM1, 
        vector<lis_t>& lisIntervalsFM2,
        int& concordant, 
        int& discordant,
        const align_opt & alnOptions,
        const paired_opt & pairedOpt){
    long M1length = (long)mate1.sequence.size();
    long M2length = (long)mate2.sequence.size();
    long editDistM1 = (long)(alnOptions.errorPercent*M1length)+1;
    long editDistM2 = (long)(alnOptions.errorPercent*M2length)+1;
    int state = 0;
    int i = 0;
    int j = 0;
    int begin = j;
    while(i < lisIntervalsFM1.size() && concordant < alnOptions.alignmentCount){
        state = 0;
        begin = j;
        lis_t & lis1 = lisIntervalsFM1[i];
        while(j < lisIntervalsFM2.size() && concordant < alnOptions.alignmentCount && state < 2){
            lis_t & lis2 = lisIntervalsFM2[j];
            if(isPosConcordant(lis1, lis2, editDistM1, editDistM2, M1length, M2length, pairedOpt)){//concordant (pos)
                if(state==0){ state++; begin = j;}//first pos where a concordant alignment was found
                //Align if necessary, 
                bool firstAligned = alignFromLIS(sa, dp_, mate1, lis1, editDistM1, alnOptions);
                bool secondAligned = alignFromLIS(sa, dp_, mate2, lis2, editDistM2, alnOptions);
                //If both align: check pairing and set conc/disc: add to pairedInfo
                if(firstAligned && secondAligned){
                    if(isConcordant(lis1.alignment, lis2.alignment, pairedOpt)){//concordant
                        concordant++;
                        lis1.len = 0; lis2.len = 0;
                        setPaired(lis1.alignment, lis2.alignment,mate1,mate2,true);
                    }
                    else if(pairedOpt.discordant){
                        lis1.len = 0; lis2.len = 0;
                        discordant++;
                        setPaired(lis1.alignment,lis2.alignment,mate1,mate2, false);
                    }
                }
            }
            else if(state==1)//discordant, move on
                state++;
            j++;
        }
        j = begin;
        i++;
    }
}

void unpairedMatchFromLIS(const sparseSA& sa, dynProg& dp_, read_t & read, vector<lis_t> & lisIntervals, const align_opt & alnOptions, int & alnCount){
    long editDist = (long)(alnOptions.errorPercent*read.sequence.size())+1;
    //sort candidate regions for likelyhood of an alignment
    sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
    //for every interval, try to align
    int lisIndex = 0;
    //trial parameter for performance reasons
    int trial = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    while(alnCount < alnOptions.alignmentCount &&
            lisIndex < lisIntervals.size() &&
            lisIntervals[lisIndex].len > (read.sequence.size()*alnOptions.minCoverage)/100 &&
            trial < alnOptions.maxTrial){
        //sort matches in this interval according to query position
        int begin = lisIntervals[lisIndex].begin;
        int end = lisIntervals[lisIndex].end;
        vector<match_t>* matchVector = lisIntervals[lisIndex].matches;
        //sort this candidate region by query position
        sort(matchVector->begin()+begin,matchVector->begin()+end+1, compMatchesQuery);
        bool extended = extendAlignment(dp_, sa.S, lisIntervals[lisIndex].fw ? read.sequence : read.rcSequence, 
                lisIntervals[lisIndex].alignment, *matchVector, begin, end, editDist, alnOptions);
        if(extended){
            if(!lisIntervals[lisIndex].fw){
                lisIntervals[lisIndex].alignment.flag.set(4,true);
            }
            alnCount++;
            lisIntervals[lisIndex].alignment.setLocalPos(sa);
            lisIntervals[lisIndex].extended = true;
            trial = 0;
        }
        lisIndex++;
        trial++;
    }
}

//Calculate LIS of both mates, match together 2 LIS, if match, calculate both alignments
void pairedMatch3(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    int concordant = 0;
    int discordant = 0;
    int min_len = alnOptions.minMemLength;
    if((alnOptions.noFW || alnOptions.noRC))
        return;//no possible alignments
    //fields
    bool mate1FWfirst = pairedOpt.orientation == PAIR_FR || pairedOpt.orientation == PAIR_FF;
    bool mate2FWfirst = pairedOpt.orientation == PAIR_RF || pairedOpt.orientation == PAIR_FF;
    vector<match_t> firstmatchesM1;
    vector<match_t> firstmatchesM2;
    vector<match_t> secondmatchesM1;
    vector<match_t> secondmatchesM2;
    vector<lis_t> lisIntervalsFM1;
    vector<lis_t> lisIntervalsFM2;
    vector<lis_t> lisIntervalsSM1;
    vector<lis_t> lisIntervalsSM2;
    long M1length = (long) mate1.sequence.size();
    long M2length = (long) mate2.sequence.size();
    long editDistM1 = (long)(alnOptions.errorPercent*M1length)+1;
    long editDistM2 = (long)(alnOptions.errorPercent*M2length)+1;
    //Calc seeds for first direction (put this in separate function for each direction
    calculateSeeds(sa, mate1FWfirst ? mate1.sequence : mate1.rcSequence, min_len, alnOptions.alignmentCount, firstmatchesM1, alnOptions.tryHarder);
    calculateSeeds(sa, mate2FWfirst ? mate2.sequence : mate2.rcSequence, min_len, alnOptions.alignmentCount, firstmatchesM2, alnOptions.tryHarder);
    if(firstmatchesM1.size()>0 && firstmatchesM2.size()>0){
        //Calc LIS intervals
        //lis intervals are intervals in seed matches that are ordered according to reference offset
        calculateLISintervals(firstmatchesM1, mate1FWfirst, M1length, editDistM1, lisIntervalsFM1);
        calculateLISintervals(firstmatchesM2, mate2FWfirst, M2length, editDistM2, lisIntervalsFM2);
        //find concordant intervals
        coupleLIS(sa, dp_, mate1, mate2, lisIntervalsFM1, lisIntervalsFM2, concordant, discordant, alnOptions, pairedOpt);
    }
    //Do the same for other direction if necessary
    if(concordant < alnOptions.alignmentCount){
        calculateSeeds(sa, mate1FWfirst ? mate1.rcSequence : mate1.sequence, min_len, alnOptions.alignmentCount, secondmatchesM1, alnOptions.tryHarder);
        calculateSeeds(sa, mate2FWfirst ? mate2.rcSequence : mate2.sequence, min_len, alnOptions.alignmentCount, secondmatchesM2, alnOptions.tryHarder);
        if(secondmatchesM1.size()>0 && secondmatchesM2.size()>0){
            //Calc LIS intervals
            //lis intervals are intervals in seed matches that are ordered according to reference offset
            calculateLISintervals(secondmatchesM1, !mate1FWfirst, M1length, editDistM1, lisIntervalsSM1);
            calculateLISintervals(secondmatchesM2, !mate2FWfirst, M2length, editDistM2, lisIntervalsSM2);
            //find concordant intervals
            coupleLIS(sa, dp_, mate1, mate2, lisIntervalsSM1, lisIntervalsSM2, concordant, discordant, alnOptions, pairedOpt);
        }
    }
    concordant = concordant + discordant;
    int maxScoreFirst;
    int maxScoreSecond;
    //preprocessing for discordant and mixed: calculate more alignments
    if(concordant < alnOptions.alignmentCount && (pairedOpt.discordant || pairedOpt.mixed)){
        //concatenate the lists per mate
        lisIntervalsFM1.insert(lisIntervalsFM1.end(),lisIntervalsSM1.begin(),lisIntervalsSM1.end());
        lisIntervalsFM2.insert(lisIntervalsFM2.end(),lisIntervalsSM2.begin(),lisIntervalsSM2.end());
        //try to align some more LIS intervals
        unpairedMatchFromLIS(sa, dp_, mate1, lisIntervalsFM1, alnOptions, concordant);
        unpairedMatchFromLIS(sa, dp_, mate2, lisIntervalsFM2, alnOptions, concordant);
        //sort according to alignmentScore and globalPos
        sort(lisIntervalsFM1.begin(),lisIntervalsFM1.end(), compAlignmentAndScore);
        sort(lisIntervalsFM2.begin(),lisIntervalsFM2.end(), compAlignmentAndScore);
        maxScoreFirst = lisIntervalsFM1[0].alignment.alignmentScore;
        maxScoreSecond = lisIntervalsFM2[0].alignment.alignmentScore;
    }
    //Set discordant to LIS, count something else, and make sure we can retrieve the right ones
    if(pairedOpt.discordant && concordant < alnOptions.alignmentCount ){
        //extra discordant searches depending on max score
        int i=0;
        while(i < lisIntervalsFM1.size() && 
                concordant < alnOptions.alignmentCount && 
                lisIntervalsFM1[i].extended &&
                lisIntervalsFM1[i].alignment.alignmentScore == maxScoreFirst){
            set<long> globPosMate;
            if(lisIntervalsFM1[i].alignment.paired())
                for(int idx=0; idx < lisIntervalsFM1[i].alignment.pairedCount(); idx++)
                    globPosMate.insert(lisIntervalsFM1[i].alignment.mateInfo[idx].pnextGlob);
            int j = 0;
            while(j < lisIntervalsFM2.size() && 
                    concordant < alnOptions.alignmentCount &&
                    lisIntervalsFM2[j].extended &&
                    lisIntervalsFM2[j].alignment.alignmentScore == maxScoreSecond){
                if(globPosMate.count(lisIntervalsFM2[j].alignment.globPos)==0){
                    //found extra discordant
                    concordant++;
                    setPaired(lisIntervalsFM1[i].alignment,lisIntervalsFM2[j].alignment,mate1,mate2, false);
                }
                j++;
            }
            i++;
        }
    }
    if(pairedOpt.mixed && concordant < alnOptions.alignmentCount ){
        //extra unpaired alignments
        int i = 0;
        int alnCountFirst = concordant;
        while(i < lisIntervalsFM1.size() && lisIntervalsFM1[i].extended && alnCountFirst < alnOptions.alignmentCount){
            if(!lisIntervalsFM1[i].alignment.paired()){
                setUnPaired(lisIntervalsFM1[i].alignment,mate1, true);
                alnCountFirst++;
            }
            i++;
        }
        i = 0;
        int alnCountSecond = concordant;
        while(i < lisIntervalsFM2.size() && lisIntervalsFM2[i].extended && alnCountSecond < alnOptions.alignmentCount){
            if(!lisIntervalsFM2[i].alignment.paired()){
                setUnPaired(lisIntervalsFM2[i].alignment,mate2, false);
                alnCountSecond++;
            }
            i++;
        }
    }
    for(int i=0; i < lisIntervalsFM1.size(); i++ )
        if(lisIntervalsFM1[i].extended)
            mate1.alignments.push_back(lisIntervalsFM1[i].alignment);
    for(int i=0; i < lisIntervalsFM2.size(); i++ )
        if(lisIntervalsFM2[i].extended)
            mate2.alignments.push_back(lisIntervalsFM2[i].alignment);
}

bool isConcordantAlnToLis(const alignment_t& mate1, const lis_t& mate2, bool mate1isFirst, long editDist2, long M2length, const paired_opt& options){
    bool firstLeft = true;
    bool mate1FW = !mate1.flag.test(4);
    bool mate2FW = mate2.fw;
    //Check strand
    if(options.orientation == PAIR_FF || options.orientation == PAIR_FR){
        firstLeft = mate1isFirst ? mate1FW : mate2FW;//if rc must be right
    }
    else{
        firstLeft = mate1isFirst ? !mate1FW : !mate2FW;//if rc must be left
    }
    //find possible offset
    long posLeft1 = mate1.globPos;
    long posLeft2 = mate2.matches->at(mate2.begin).ref-mate2.matches->at(mate2.begin).query;
    long posRight1 = posLeft1 + (long) mate1.refLength;
    long posRight2 = posLeft2 + M2length;
    //check fragment insert size
    long minsize = max(posRight2+editDist2,posRight1) - min(posLeft2+editDist2,posLeft1) +1L;
    if(minsize < options.minInsert)
        return false;
    long maxsize = max(posRight2-editDist2,posRight1) - min(posLeft2-editDist2,posLeft1) +1L;
    if(maxsize > options.maxInsert)
        return false;
    //check flow direction is correct
    long anchor1 = mate1FW ? posLeft1 : posRight1;
    long anchor2 = mate2FW ? posLeft2 : posRight2;
    if(firstLeft){
        if(anchor1 > anchor2 + editDist2) return false;
    }
    else{
        if(anchor2 - editDist2 > anchor1) return false;
    }
    //check contain
    bool contain = ((posLeft1 >= posLeft2+editDist2 && posRight1 <= posRight2-editDist2) 
    || (posLeft2-editDist2 >= posLeft1 && posRight2+editDist2 <= posRight1));
    if(contain && !options.contain){
        return false;
    }
    //check overlap
    bool overlap = contain || !(posRight2-editDist2 < posLeft1 || posLeft2+editDist2 > posRight1);
    if(!options.overlap && overlap){
        return false;
    }
    //check dovetail
    if((firstLeft && (posLeft2+editDist2 < posLeft1)) || (!firstLeft && (posLeft1 < posLeft2-editDist2)) ){
        if(!overlap)
            return false;
        else if(!options.dovetail)
            return false;
    }
    return true;
}

void matchStrandOfPair(const sparseSA& sa, 
        dynProg & dp_,
        read_t & mate1, 
        read_t & mate2,
        vector<lis_t>& lisIntervalsFM1, 
        vector<lis_t>& lisIntervalsFM2,
        int& concordant, 
        int& discordant,
        bool mate1isFirst,
        const align_opt & alnOptions,
        const paired_opt & pairedOpt){
    long M1length = (long) mate1.sequence.size();
    long M2length = (long) mate2.sequence.size();
    long editDistM1 = (long)(alnOptions.errorPercent*M1length)+1;
    long editDistM2 = (long)(alnOptions.errorPercent*M2length)+1;
    //Matching strand and mate: sort acocrding to score
    sort(lisIntervalsFM1.begin(), lisIntervalsFM1.end(), compIntervals);
    //Other strand: for pairing
    sort(lisIntervalsFM2.begin(), lisIntervalsFM2.end(), compIntervalsRef);
    int lisIndex = 0;
    //trial parameter for performance reasons
    int trial = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    while(concordant < alnOptions.alignmentCount && 
            lisIndex < lisIntervalsFM1.size() && 
            lisIntervalsFM1[lisIndex].len > (M1length*alnOptions.minCoverage)/100 && 
            trial < alnOptions.maxTrial){
        //Alignment of mate1
        if(!lisIntervalsFM1[lisIndex].extended){
            //sort matches in this interval according to query position
            int begin = lisIntervalsFM1[lisIndex].begin;
            int end = lisIntervalsFM1[lisIndex].end;
            vector<match_t>* matchVector = lisIntervalsFM1[lisIndex].matches;
            //sort this candidate region by query position
            sort(matchVector->begin()+begin,matchVector->begin()+end+1, compMatchesQuery);
            lisIntervalsFM1[lisIndex].extended = extendAlignment(dp_, sa.S, lisIntervalsFM1[lisIndex].fw ? mate1.sequence : mate1.rcSequence, 
                    lisIntervalsFM1[lisIndex].alignment, *matchVector, begin, end, editDistM1, alnOptions);
            if(lisIntervalsFM1[lisIndex].extended){
                if(!lisIntervalsFM1[lisIndex].fw)
                    lisIntervalsFM1[lisIndex].alignment.flag.set(4,true);
                lisIntervalsFM1[lisIndex].alignment.setLocalPos(sa);
            }
        }
        //If an alignment, continue to find match for other mate
        if(lisIntervalsFM1[lisIndex].extended){
            trial = 0;
            lis_t & lis1 = lisIntervalsFM1[lisIndex];
            int state = 0;
            int j = 0;
            while(j < lisIntervalsFM2.size() && concordant < alnOptions.alignmentCount && state < 2){
                lis_t & lis2 = lisIntervalsFM2[j];
                if((lis2.extended && isConcordant(mate1isFirst ? lis1.alignment : lis2.alignment, 
                        mate1isFirst ? lis2.alignment : lis1.alignment, pairedOpt)) ){
                    state = 1;
                    if(mate1isFirst){//if not: repeated find
                        concordant++;
                        setPaired(lis1.alignment,lis2.alignment,mate1,mate2, true);
                    }
                }
                else if(!lis2.extended && isConcordantAlnToLis(lis1.alignment, lis2, mate1isFirst, editDistM2, M2length, pairedOpt)){
                    state=1;
                    if(alignFromLIS(sa, dp_, mate2, lis2, editDistM2, alnOptions)){
                        if(isConcordant(mate1isFirst ? lis1.alignment : lis2.alignment, 
                                mate1isFirst ? lis2.alignment : lis1.alignment, pairedOpt)){//concordant
                            concordant++;
                            mate1isFirst ? setPaired(lis1.alignment,lis2.alignment,mate1,mate2 , true) : 
                                setPaired(lis2.alignment,lis1.alignment,mate2,mate1 , true);
                        }
                        else if(pairedOpt.discordant){
                            discordant++;
                            mate1isFirst ? setPaired(lis1.alignment,lis2.alignment,mate1,mate2 , false) : 
                                setPaired(lis2.alignment,lis1.alignment,mate2,mate1 , false);
                        }
                    }
                }
                else if(state==1)//discordant, move on
                    state++;
                j++;
            }
        }
        lisIndex++;
        trial++;
    }
}

//Calculate LIS of both mates, align 1, match LIS to other, align 2
void pairedMatch4(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    int concordant = 0;
    int discordant = 0;
    int min_len = alnOptions.minMemLength;
    if((alnOptions.noFW || alnOptions.noRC))
        return;//no possible alignments
    //fields
    bool mate1FWfirst = pairedOpt.orientation == PAIR_FR || pairedOpt.orientation == PAIR_FF;
    bool mate2FWfirst = pairedOpt.orientation == PAIR_RF || pairedOpt.orientation == PAIR_FF;
    vector<match_t> firstmatchesM1;
    vector<match_t> firstmatchesM2;
    vector<match_t> secondmatchesM1;
    vector<match_t> secondmatchesM2;
    vector<lis_t> lisIntervalsFM1;
    vector<lis_t> lisIntervalsFM2;
    vector<lis_t> lisIntervalsSM1;
    vector<lis_t> lisIntervalsSM2;
    long M1length = (long) mate1.sequence.size();
    long M2length = (long) mate2.sequence.size();
    long editDistM1 = (long)(alnOptions.errorPercent*M1length)+1;
    long editDistM2 = (long)(alnOptions.errorPercent*M2length)+1;
    //Calc seeds for first direction (put this in separate function for each direction
    calculateSeeds(sa, mate1FWfirst ? mate1.sequence : mate1.rcSequence, min_len, alnOptions.alignmentCount, firstmatchesM1, alnOptions.tryHarder);
    calculateSeeds(sa, mate2FWfirst ? mate2.sequence : mate2.rcSequence, min_len, alnOptions.alignmentCount, firstmatchesM2, alnOptions.tryHarder);
    if(firstmatchesM1.size()>0 && firstmatchesM2.size()>0){
        //Calc LIS intervals
        //lis intervals are intervals in seed matches that are ordered according to reference offset
        calculateLISintervals(firstmatchesM1, mate1FWfirst, M1length, editDistM1, lisIntervalsFM1);
        calculateLISintervals(firstmatchesM2, mate2FWfirst, M2length, editDistM2, lisIntervalsFM2);
        //find concordant intervals
        matchStrandOfPair(sa, dp_, mate1, mate2, lisIntervalsFM1, lisIntervalsFM2, concordant, discordant, true, alnOptions, pairedOpt);
    }
    //Do the same for other direction if necessary
    if(concordant < alnOptions.alignmentCount){
        calculateSeeds(sa, mate1FWfirst ? mate1.rcSequence : mate1.sequence, min_len, alnOptions.alignmentCount, secondmatchesM1, alnOptions.tryHarder);
        calculateSeeds(sa, mate2FWfirst ? mate2.rcSequence : mate2.sequence, min_len, alnOptions.alignmentCount, secondmatchesM2, alnOptions.tryHarder);
        if(secondmatchesM1.size()>0 && secondmatchesM2.size()>0){
            //Calc LIS intervals
            //lis intervals are intervals in seed matches that are ordered according to reference offset
            calculateLISintervals(secondmatchesM1, !mate1FWfirst, M1length, editDistM1, lisIntervalsSM1);
            calculateLISintervals(secondmatchesM2, !mate2FWfirst, M2length, editDistM2, lisIntervalsSM2);
            //find concordant intervals
            matchStrandOfPair(sa, dp_, mate1, mate2, lisIntervalsSM1, lisIntervalsSM2, concordant, discordant, true, alnOptions, pairedOpt);
        }
    }
    //Do the same for first direction of other mate
    if(concordant < alnOptions.alignmentCount){
        if(firstmatchesM1.size()>0 && firstmatchesM2.size()>0){
            //find concordant intervals
            matchStrandOfPair(sa, dp_, mate2, mate1, lisIntervalsFM2, lisIntervalsFM1, concordant, discordant, false, alnOptions, pairedOpt);
        }
    }
    //Do the same for second direction of other mate
    if(concordant < alnOptions.alignmentCount){
        if(secondmatchesM1.size()>0 && secondmatchesM2.size()>0){
            //find concordant intervals
            matchStrandOfPair(sa, dp_, mate2, mate1, lisIntervalsSM2, lisIntervalsSM1, concordant, discordant, false, alnOptions, pairedOpt);
        }
    }
    concordant = concordant + discordant;
    int maxScoreFirst;
    int maxScoreSecond;
    //preprocessing for discordant and mixed: calculate more alignments
    if(concordant < alnOptions.alignmentCount && (pairedOpt.discordant || pairedOpt.mixed)){
        //concatenate the lists per mate
        lisIntervalsFM1.insert(lisIntervalsFM1.end(),lisIntervalsSM1.begin(),lisIntervalsSM1.end());
        lisIntervalsFM2.insert(lisIntervalsFM2.end(),lisIntervalsSM2.begin(),lisIntervalsSM2.end());
        //sort according to alignmentScore and globalPos
        sort(lisIntervalsFM1.begin(),lisIntervalsFM1.end(), compAlignmentAndScore);
        sort(lisIntervalsFM2.begin(),lisIntervalsFM2.end(), compAlignmentAndScore);
        maxScoreFirst = lisIntervalsFM1[0].alignment.alignmentScore;
        maxScoreSecond = lisIntervalsFM2[0].alignment.alignmentScore;
    }
    if(pairedOpt.discordant && concordant < alnOptions.alignmentCount ){
        //extra discordant searches depending on max score
        int i=0;
        while(i < lisIntervalsFM1.size() && 
                concordant < alnOptions.alignmentCount && 
                lisIntervalsFM1[i].extended &&
                lisIntervalsFM1[i].alignment.alignmentScore == maxScoreFirst){
            set<long> globPosMate;
            if(lisIntervalsFM1[i].alignment.paired())
                for(int idx=0; idx < lisIntervalsFM1[i].alignment.pairedCount(); idx++)
                    globPosMate.insert(lisIntervalsFM1[i].alignment.mateInfo[idx].pnextGlob);
            int j = 0;
            while(j < lisIntervalsFM2.size() && 
                    concordant < alnOptions.alignmentCount &&
                    lisIntervalsFM2[j].extended &&
                    lisIntervalsFM2[j].alignment.alignmentScore == maxScoreSecond){
                if(globPosMate.count(lisIntervalsFM2[j].alignment.globPos)==0){
                    //found extra discordant
                    concordant++;
                    setPaired(lisIntervalsFM1[i].alignment,lisIntervalsFM2[j].alignment,mate1,mate2, false);
                }
                j++;
            }
            i++;
        }
    }
    if(pairedOpt.mixed && concordant < alnOptions.alignmentCount ){
        //extra unpaired alignments
        int i = 0;
        int alnCountFirst = concordant;
        while(i < lisIntervalsFM1.size() && lisIntervalsFM1[i].extended && alnCountFirst < alnOptions.alignmentCount){
            if(!lisIntervalsFM1[i].alignment.paired()){
                setUnPaired(lisIntervalsFM1[i].alignment,mate1, true);
                alnCountFirst++;
            }
            i++;
        }
        i = 0;
        int alnCountSecond = concordant;
        while(i < lisIntervalsFM2.size() && lisIntervalsFM2[i].extended && alnCountSecond < alnOptions.alignmentCount){
            if(!lisIntervalsFM2[i].alignment.paired()){
                setUnPaired(lisIntervalsFM2[i].alignment,mate2, false);
                alnCountSecond++;
            }
            i++;
        }
    }
    for(int i=0; i < lisIntervalsFM1.size(); i++ )
        if(lisIntervalsFM1[i].extended)
            mate1.alignments.push_back(lisIntervalsFM1[i].alignment);
    for(int i=0; i < lisIntervalsFM2.size(); i++ )
        if(lisIntervalsFM2[i].extended)
            mate2.alignments.push_back(lisIntervalsFM2[i].alignment);
}

bool dpWindow(const alignment_t& mate1, bool mate1isFirst, long editDist2, long M2length, const paired_opt& options,
        boundaries& grenzen, bool& otherFW, int& bandLeft, int& bandRight){
    bool firstLeft = true;
    bool mate1FW = !mate1.flag.test(4);
    //Check strand
    if(options.orientation == PAIR_FF){
        firstLeft = mate1isFirst ? mate1FW : !mate1FW;
        otherFW = mate1FW;
    }
    else if(options.orientation == PAIR_FR){
        firstLeft = mate1isFirst ? mate1FW : !mate1FW;
        otherFW = !mate1FW;
    }
    else{
        firstLeft = mate1isFirst ? !mate1FW : mate1FW;//if rc must be left
        otherFW = !mate1FW;
    }
    grenzen.queryB = 0;
    grenzen.queryE = M2length-1;
    long leftMate = mate1.globPos;
    long rightMate = leftMate + (long) mate1.refLength-1L;
    //Both first boundaries: maxfrag
    grenzen.refB = rightMate - (long)options.maxInsert + 1L;
    grenzen.refE = leftMate + (long) options.maxInsert - 1L;
    if(firstLeft){//Window to the right of already aligned mate
		long bandEnd  = leftMate + (long) options.minInsert -1L;
        if(!options.overlap){
            grenzen.refB = max(grenzen.refB, rightMate+1);
        }
        if(!options.dovetail)
            grenzen.refB = max(grenzen.refB, leftMate);
        //orientation
        grenzen.refB = max(grenzen.refB, (mate1FW ? leftMate : rightMate) + 1 - (otherFW ? 0 : M2length+editDist2));
        if(bandEnd < grenzen.refB)
            bandEnd = grenzen.refB;
        else{
            grenzen.refB = max(grenzen.refB, bandEnd - M2length+1 - editDist2);
        }
        bandRight = grenzen.refE - M2length + editDist2 + 1;
        bandEnd = min(grenzen.refB - editDist2, bandEnd - M2length -editDist2);
        bandLeft = abs(bandEnd - grenzen.refB +1);
    }
    else{//Window to the left of already aligned mate
        long bandEnd = rightMate - options.minInsert +1;
        if(!options.overlap){
            grenzen.refE = min(grenzen.refE, leftMate-1);
        }
        if(!options.dovetail)
            grenzen.refE = min(grenzen.refE, rightMate);
        //orientation
        grenzen.refE = min(grenzen.refE, (mate1FW ? leftMate : rightMate) - 1 + (otherFW ? M2length+editDist2 : 0));
        if(bandEnd > grenzen.refE)
            bandEnd = grenzen.refE;
        else{
            grenzen.refE = min(grenzen.refE, bandEnd + M2length-1 + editDist2);
        }
        bandLeft = editDist2;
        bandEnd = min(bandEnd + editDist2, grenzen.refE - M2length + editDist2);
        bandRight = bandEnd - grenzen.refB +1;
    }
    return (grenzen.refE > grenzen.refB && bandRight > 0 && bandRight < grenzen.refE-grenzen.refB+1);
}

void pairedBowtie2(const sparseSA& sa, 
        dynProg& dp_,
        read_t & mate1, 
        read_t & mate2,
        bool alignFirstMate,
        bool forward,
        int& concordant, 
        const align_opt & alnOptions,
        const paired_opt & pairedOpt){
    //First: matching one read in one direction unpaired
    read_t & base = alignFirstMate ? mate1 : mate2;
    read_t & other = alignFirstMate ? mate2 : mate1;
    string & P = forward ? base.sequence : base.rcSequence;
    long Plength = (long) P.size();
    long Olength = (long) other.sequence.size();
    long editDistBase = (long)(alnOptions.errorPercent*Plength)+1;
    long editDistOther = (long)(alnOptions.errorPercent*Olength)+1;
    int min_len = alnOptions.minMemLength;
    if(!alnOptions.fixedMinLength && Plength/editDistBase > min_len)
        min_len = Plength/editDistBase;
    vector<match_t> matches;
    //calc seeds
    calculateSeeds(sa, P, min_len, alnOptions.alignmentCount, matches, alnOptions.tryHarder);
    //sort matches
    if(matches.size()>0){
        vector<lis_t> lisIntervals;
        calculateLISintervals(matches, forward, Plength, editDistBase, lisIntervals);
        //sort candidate regions for likelyhood of an alignment
        sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
        int lisIndex = 0;
        int trial = 0;//paired succes times
        while(concordant < alnOptions.alignmentCount &&
                lisIndex < lisIntervals.size() &&
                lisIntervals[lisIndex].len > (Plength*alnOptions.minCoverage)/100 &&
                trial < alnOptions.maxTrial){
            //try to find out if we already aligned this LIS-interval when pairing the other mate
            int i = 0;
            while(i < other.alignmentCount() && isConcordantAlnToLis(other.alignments[i], lisIntervals[lisIndex], !alignFirstMate, editDistBase, Plength, pairedOpt))
                i++;
            if(alignFirstMate || i == other.alignmentCount()){
                //have to extend alignment
                int begin = lisIntervals[lisIndex].begin;
                int end = lisIntervals[lisIndex].end;
                sort(matches.begin()+begin,matches.begin()+end+1, compMatchesQuery);
                alignment_t alignment;
                bool extended = extendAlignment(dp_, sa.S, P, alignment, matches, begin, end, editDistBase, alnOptions);
                if(extended){
                    alignment.setLocalPos(sa);
                    if(!forward)
                        alignment.flag.set(4,true);
                    base.alignments.push_back(alignment);
                    //try to pair
                    dp_output output;
                    bool clipping = !alnOptions.noClipping;
                    dp_type types;
                    types.freeQueryB = types.freeQueryE = clipping;
                    types.freeRefB = types.freeRefE = true;
                    boundaries grenzen;
                    int bandLeft = 0;
                    int bandRight = 0;
                    bool otherFW = true;
                    //find dp window
                    if(dpWindow(alignment, alignFirstMate, editDistOther, Olength, pairedOpt, 
                            grenzen, otherFW, bandLeft, bandRight))
                        dp_.dpBandFull( sa.S, otherFW ? other.sequence : other.rcSequence, grenzen, 
                            types, ERRORSTRING, output, bandLeft, bandRight, false);
                    else
                        output.editDist = 2*editDistOther;
                    if(output.editDist <= editDistOther){
                        //SUCCES!!!!!!!
                        alignment_t mate;
                        mate.cigarChars.insert(mate.cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
                        mate.cigarLengths.insert(mate.cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
                        mate.globPos = grenzen.refB+1;
                        mate.refLength = grenzen.refE-grenzen.refB+1;
                        mate.editDist = output.editDist;
                        mate.alignmentScore = output.dpScore;
                        if(!otherFW)
                            mate.flag.set(4,true);
                        mate.setLocalPos(sa);
                        other.alignments.push_back(mate);
                        //set paired: alignments are the ones last pushed
                        setPaired(mate1.alignments[mate1.alignmentCount()-1],mate2.alignments[mate2.alignmentCount()-1],mate1,mate2, true);
                        concordant++;
                        trial = 0;
                    }
                }
            }
            else{
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
    }
}

//Calculate Mates 1 at a time, without calculating seeds for the other
//Optimalization: use 'Bailing' method: to enable/disable the calculation for the other mate
void pairedMatch2(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    int concordant = 0;
    if((alnOptions.noFW || alnOptions.noRC))
        return;//no possible alignments
    //fields
    bool mate1FWfirst = pairedOpt.orientation == PAIR_FR || pairedOpt.orientation == PAIR_FF;
    bool mate2FWfirst = pairedOpt.orientation == PAIR_RF || pairedOpt.orientation == PAIR_FF;
    //concordant
    pairedBowtie2(sa, dp_, mate1, mate2, true, mate1FWfirst, concordant, alnOptions, pairedOpt);
    if(concordant < alnOptions.alignmentCount)
        pairedBowtie2(sa, dp_, mate1, mate2, true, !mate1FWfirst, concordant, alnOptions, pairedOpt);
    if(concordant < alnOptions.alignmentCount)
        pairedBowtie2(sa, dp_, mate1, mate2, false, mate2FWfirst, concordant, alnOptions, pairedOpt);
    if(concordant < alnOptions.alignmentCount)
        pairedBowtie2(sa, dp_, mate1, mate2, false, !mate2FWfirst, concordant, alnOptions, pairedOpt);
    //
    int maxScoreFirst;
    int maxScoreSecond;
    //preprocessing for discordant and mixed: calculate more alignments
    if(concordant < alnOptions.alignmentCount && (pairedOpt.discordant || pairedOpt.mixed)){
        //sort according to alignmentScore and globalPos
        sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignmentScore);
        sort(mate2.alignments.begin(),mate2.alignments.end(), compAlignmentScore);
        if(mate1.alignmentCount()>0 && mate2.alignmentCount() > 0){
            maxScoreFirst = mate1.alignments[0].alignmentScore;
            maxScoreSecond = mate2.alignments[0].alignmentScore;
        }
    }
    if(mate1.alignmentCount()>0 && mate2.alignmentCount() > 0 && pairedOpt.discordant && concordant < alnOptions.alignmentCount ){
        //extra discordant searches depending on max score
        int i=0;
        while(i < mate1.alignmentCount() && 
                concordant < alnOptions.alignmentCount && 
                mate1.alignments[i].alignmentScore == maxScoreFirst){
            int j = 0;
            while(j < mate2.alignmentCount() && 
                    concordant < alnOptions.alignmentCount &&
                    mate2.alignments[j].alignmentScore == maxScoreSecond){
                if((!mate1.alignments[i].paired() || !mate2.alignments[j].paired()) || 
                        mate1.alignments[i].mateInfo[0].pnextGlob!=mate2.alignments[j].globPos){
                    //found extra discordant
                    concordant++;
                    setPaired(mate1.alignments[i],mate2.alignments[j],mate1,mate2, false);
                }
                j++;
            }
            i++;
        }
    }
    if(pairedOpt.mixed && concordant < alnOptions.alignmentCount ){
        //extra unpaired alignments
        int i = 0;
        int alnCountFirst = concordant;
        while(i < mate1.alignmentCount() && alnCountFirst < alnOptions.alignmentCount){
            if(!mate1.alignments[i].paired()){
                setUnPaired(mate1.alignments[i],mate1, true);
                alnCountFirst++;
            }
            i++;
        }
        i = 0;
        int alnCountSecond = concordant;
        while(i < mate2.alignmentCount() && alnCountSecond < alnOptions.alignmentCount){
            if(!mate2.alignments[i].paired()){
                setUnPaired(mate2.alignments[i],mate2, false);
                alnCountSecond++;
            }
            i++;
        }
    }
}

//MODE 3 AND 4 NOT WORKING: [3] crash, [4] not pairing!!!
void pairedMatch(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt, bool print){
    if(pairedOpt.mode==1){
        pairedMatch1(sa, dp_, mate1, mate2, alnOptions, pairedOpt, print);
    }
    else if(pairedOpt.mode==2){
        pairedMatch2(sa, dp_, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==3){
        pairedMatch3(sa, dp_, mate1, mate2, alnOptions, pairedOpt);
    }
    else if(pairedOpt.mode==4){
        pairedMatch4(sa, dp_, mate1, mate2, alnOptions, pairedOpt);
    }
    else{
        //mode unknown
    }
}