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

bool compMatches(const match_t & i, const match_t & j){
    return (i.ref < j.ref || (i.ref == j.ref && i.query < j.query));
    //return (i.query < j.query || (i.query == j.query && i.ref < j.ref));
}

bool compMatchesQuery(const match_t & i,const match_t & j){
    //return (i.ref < j.ref || (i.ref == j.ref && i.query < j.query));
    return (i.query < j.query || (i.query == j.query && i.ref < j.ref));
}

bool compIntervals(const lis_t & i,const lis_t & j){
    return (i.len > j.len || (i.len == j.len && i.begin < j.begin));
}

bool compIntervalsRef(const lis_t & i,const lis_t & j){
    return compMatches(i.matches->at(i.begin), j.matches->at(j.begin));
}

bool compAlignmentScore(const alignment_t * i,const alignment_t * j){
    return (i->alignmentScore > j->alignmentScore || (i->alignmentScore == j->alignmentScore && i->globPos < j->globPos));
}

bool compAlignmentAndScore(const lis_t & i,const lis_t & j){
    return (i.extended&&!j.extended )|| 
            ((i.extended && j.extended ) && 
            (i.alignment->alignmentScore > j.alignment->alignmentScore 
            || (i.alignment->alignmentScore == j.alignment->alignmentScore && 
            (i.alignment->globPos < j.alignment->globPos))));
}

bool compPairScore(const alignment_t * i,const alignment_t * j){
    if(i->concordant() > j->concordant()) return true;
    else if(i->concordant() == j->concordant()){
        return (i->flag.test(1) > j->flag.test(1) || ( i->flag.test(1) == j->flag.test(1) && i->flag.test(0) > j->flag.test(0)) ||
    (i->flag.test(0) == j->flag.test(0) && i->flag.test(1) == j->flag.test(1) && i->alignmentScore > j->alignmentScore));
    }
    else{
        return false;
    }
}

void postProcess(vector<match_t> &matches){
    sort(matches.begin(),matches.end(), compMatches);
}

void calculateSeeds(const sparseSA& sa, const string& P, int min_len, int maxBranchWidth, vector<match_t>& matches, bool tryHarder, mum_t memType){
    if(memType == SMAM)
        sa.SMAM(P, matches, min_len, maxBranchWidth);
    else if(memType == MEM)
        sa.MEM(P, matches, min_len, maxBranchWidth);
    else if(memType == MAM)
        sa.MAM(P, matches, min_len, maxBranchWidth);
    else if(memType == MUM)
        sa.MUM(P, matches, min_len, maxBranchWidth);
    if(tryHarder){//TODO: change try-harder to recalculate only after forward + reverse has been tried
        //easy solution: try reverse and if found: skip
        if(matches.empty())
            sa.SMAM(P, matches, min_len, 1000);
        if(matches.empty())
            sa.SMAM(P, matches, 20, 1000);
    }
    postProcess(matches);
}

void calculateLISintervals(vector<match_t>& matches, bool fw, long qLength, long editDist, vector<lis_t>& lisIntervals){
    if(matches.size() > 0){
        size_t begin = 0;
        size_t end = 0;
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
    size_t begin = 0;
    size_t end = 1;
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

alignment_t * extendAlignment(dynProg& dp_, const string& S, const string& P, 
        vector<match_t>& matches, int begin, int end, long editDist, 
        const align_opt & alnOptions, long chrStart, long chrEnd){
    int print = alnOptions.print;//print lvl 1: only display steps and results, lvl2: display seeds
    alignment_t * alignment = new alignment_t();
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
    if(print) cerr << "max edit distance is " << editDist << endl;
    if(print>1) cerr << "first seed (q/r/l): " << firstSeed.query << "," << firstSeed.ref << "," << firstSeed.len << endl;
    if(firstSeed.query > 0L && firstSeed.ref > chrStart){
        if(print) cerr << "dp required before first seed";
        //Alignment now!!! with beginQ-E in ref to begin match and beginQ to match
        long alignmentBoundLeft = max(refstrLB-queryLB-editDist,chrStart);
        boundaries grenzen(alignmentBoundLeft,refstrLB-1,0,queryLB-1);
        dp_type types;
        types.freeRefB = true;
        types.freeQueryB = clipping;
        if(print) cerr << "dp called with dimension " << (grenzen.queryE-grenzen.queryB+1) << "x" << (grenzen.refE-grenzen.refB+1) << endl;
        dp_.dpBandStatic( S, P, grenzen, types, ERRORSTRING, output, editDist-curEditDist, print>1);
        if(curEditDist + output.editDist <= editDist){//required if dp_ returns fail or too high editDist (output may not be initialized
            queryLB = grenzen.queryB;
            if(output.cigarChars[0] == 'I')
                    queryLB += (long)output.cigarLengths[0];
            if(clipping && grenzen.queryB> 0){
                alignment->cigarChars.push_back('S');
                alignment->cigarLengths.push_back(grenzen.queryB);
            }
            alignment->cigarChars.insert(alignment->cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
            alignment->cigarLengths.insert(alignment->cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
            alignment->globPos = grenzen.refB+1L;
            alignment->alignmentScore += output.dpScore;
        }
        curEditDist += output.editDist;
        if(print) cerr << "begin DP returned with edit distance " << output.editDist << endl;
        output.clear();
    }//fill in the starting position in the ref sequence
    else {
        alignment->globPos = firstSeed.ref + 1L;
    }
    if (firstSeed.ref == chrStart && firstSeed.query > 0L) {
        alignment->cigarChars.push_back(clipping ? 'S' : 'I');
        alignment->cigarLengths.push_back(firstSeed.query);
        curEditDist += clipping ? 0 : firstSeed.query;
        alignment->alignmentScore += (clipping ? 0 : dp_.scores.openGap + dp_.scores.extendGap*firstSeed.query);
    }
    if(print) cerr << "current alignments startpos would be " << alignment->globPos << endl;
    alignment->cigarChars.push_back('=');
    alignment->cigarLengths.push_back(firstSeed.len);
    alignment->alignmentScore += dp_.scores.match*firstSeed.len;
    //extend seed to the right, filling gaps between mems
    bool foundNext = true;
    int memsIndex = begin+1;
    //////////////////////
    //NEXT SEEDS LOOP
    //////////////////////
    if(print) cerr << "enter chain loop" << endl;
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
            if(print>1) cerr << "try " << match.query << "," << match.ref << "," << match.len << endl;
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
                    if (refDist + refstrRB > chrEnd || 1 - qDist >= match.len)
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
        if(print>1){
            cerr << "previous seed number " << (memsIndex-1-begin) << " was " << matches[memsIndex].query << "," << matches[memsIndex].ref << "," << matches[memsIndex].len << endl;
            if(foundNext){
                for(int i = memsIndex; i < minDistMem; i++)
                    cerr << "seed skipped: " << matches[i].query << "," << matches[i].ref << "," << matches[i].len << endl;
            }
            else
                cerr << "found no new seed" << endl;
        }
        if (foundNext) {
            match_t match = matches[minDistMem];
            //add indels and mutations to sol
            if (minQDist <= 0) {
                //minRefDist deletions
                //1+|qDist| matches lost of next mem
                alignment->cigarChars.push_back('D');
                alignment->cigarChars.push_back('=');
                alignment->cigarLengths.push_back(minRefDist);
                alignment->cigarLengths.push_back(match.len - 1 + minQDist);
                alignment->alignmentScore += dp_.scores.match*(match.len - 1 + minQDist) + dp_.scores.openGap + dp_.scores.extendGap*minRefDist;
                curEditDist += minRefDist; //boundary of ref sequences is passed
                if(print) cerr << "added deletion of size " << minRefDist << endl;
            } else if (minRefDist <= 0) {
                alignment->cigarChars.push_back('I');
                alignment->cigarChars.push_back('=');
                alignment->cigarLengths.push_back(minQDist);
                alignment->cigarLengths.push_back(match.len - 1 + minRefDist);
                alignment->alignmentScore += dp_.scores.match*(match.len - 1 + minRefDist) + dp_.scores.openGap + dp_.scores.extendGap*minQDist;
                curEditDist += minQDist;
                if(print) cerr << "added insertion of size " << minQDist << endl;
            } else {//both distances are positive and not equal to (1,1)
                boundaries grenzen(refstrRB + 1L, refstrRB + minRefDist - 1L, queryRB + 1, queryRB + minQDist - 1);
                dp_type types;
                if(print) cerr << "dp called with dimension " << (grenzen.queryE-grenzen.queryB+1) << "x" << (grenzen.refE-grenzen.refB+1) << endl;
                dp_.dpBandStatic( S, P, grenzen, types, ERRORSTRING, output, editDist-curEditDist, print>1);
                if(curEditDist + output.editDist <= editDist){//required if dp_ returns fail or too high editDist (output may not be initialized
                    alignment->cigarChars.insert(alignment->cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
                    alignment->cigarLengths.insert(alignment->cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
                    alignment->cigarChars.push_back('=');
                    alignment->cigarLengths.push_back(match.len);
                    alignment->alignmentScore += dp_.scores.match*match.len + output.dpScore;
                }
                curEditDist += output.editDist;
                if(print) cerr << "DP returned with edit distance " << output.editDist << endl;
                output.clear();
            }
            if(print) cerr << "new seed number " << (minDistMem-begin) << " was " << match.query << "," << match.ref << "," << match.len << endl;
            //set new positions
            memsIndex = minDistMem + 1;
            //add matches
            refstrRB = match.ref + match.len - 1L;
            queryRB = match.query + match.len - 1L;
        }
    }
    if(print){
        cerr << "stopped chain because ";
        if(memsIndex > end) cerr << "all seeds in cluster were tried." << endl;
        if(!foundNext) cerr << "no seed found within limits of the previous seed" << endl;     
        if(curEditDist >= editDist) cerr << "max edit distance was reached" << endl;     
    }
    /////////////
    //FINAL DP?
    ////////////
    //possible at the end characters
    //mapInterval += Integer.toString(queryRB+1)+"]";
    alignment->refLength = queryRB-queryLB;
    if (queryRB + 1 < (long)P.length() && refstrRB < chrEnd && curEditDist <= editDist) {
        long refAlignRB = refstrRB + ((long)P.length() - queryRB) + 1L + editDist - curEditDist;
        if (refAlignRB > chrEnd)
            refAlignRB = chrEnd;
        boundaries grenzen(refstrRB + 1L, refAlignRB, queryRB + 1, P.length() - 1);
        dp_type types;
        types.freeRefE = true;
        types.freeQueryE = clipping;
        if(print) cerr << "END dp called with dimension " << (grenzen.queryE-grenzen.queryB+1) << "x" << (grenzen.refE-grenzen.refB+1) << endl;
        dp_.dpBandStatic( S, P, grenzen, types, ERRORSTRING, output, editDist-curEditDist+1, print>1);
        if(curEditDist + output.editDist <= editDist){//required if dp_ returns fail or too high editDist (output may not be initialized
            if(grenzen.queryE > queryRB){
                int addToLength = grenzen.queryE-queryRB;
                if(output.cigarChars[output.cigarChars.size()-1] == 'I')
                    addToLength -= output.cigarLengths[output.cigarLengths.size()-1];
                alignment->refLength += addToLength;
            }
            alignment->cigarChars.insert(alignment->cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
            alignment->cigarLengths.insert(alignment->cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
            alignment->alignmentScore += output.dpScore;
            if (clipping && grenzen.queryE < (long)P.length() - 1) {
                alignment->cigarChars.push_back('S');
                alignment->cigarLengths.push_back(P.length() - 1 - grenzen.queryE);
            }
        }
        if(print) cerr << "end DP returned with edit distance " << output.editDist << endl;
        curEditDist += output.editDist;
        output.clear();
    } else if (queryRB + 1 < (long)P.length() && curEditDist < editDist) {
        alignment->cigarChars.push_back(clipping ? 'S' : 'I');
        alignment->cigarLengths.push_back(P.length() - queryRB - 1);
        alignment->alignmentScore += (clipping ? 0 : ((P.length() - 1 - queryRB)*dp_.scores.extendGap + dp_.scores.openGap));
        curEditDist += clipping ? 0 : P.length() - 1 - queryRB;
    }
    //TODO: check for possible optimizations (static initialisations
    //TODO: reorder matches according to best hits
    if(curEditDist <= editDist){
        if(print) cerr << "extension returned an alignment with edit distance " << curEditDist << endl;
        alignment->editDist = curEditDist;
        return alignment;
    }
    else{
        if(print) cerr << "extension failed because edit distance " << curEditDist << ">" << editDist << endl;
        delete alignment;
        return NULL;
    }
}

void inexactMatch(const sparseSA& sa, dynProg& dp_, read_t & read,const align_opt & alnOptions, bool fwStrand){
    bool print = (alnOptions.print>0);
    string P = read.sequence;
    long Plength = (long) P.length();
    long editDist = (long)(alnOptions.errorPercent*Plength)+1;
    if(!fwStrand)
        P = read.rcSequence;
    int min_len = alnOptions.minMemLength;
    if(Plength/editDist > min_len){
        min_len = Plength/editDist;
        if(print) cerr << "min length adjusted to " << min_len << endl;
    }
    vector<match_t> matches;
    //calc seeds
    if(print) cerr << "calculate seeds" << endl;
    calculateSeeds(sa, P, min_len, alnOptions.maxSeedCandidates, matches, alnOptions.tryHarder, alnOptions.memType);
    if(print) cerr << "found " << matches.size() << " seeds"<< endl;
    //sort matches
    if(matches.size()>0){
        /////////////////////////
        //FIND CANDIDATE REGIONS
        /////////////////////////
        vector<lis_t> lisIntervals;
        if(print) cerr << "calculate clusters" << endl;
        calculateLISintervals(matches, fwStrand, Plength, editDist, lisIntervals);
        if(print) cerr << "found " << lisIntervals.size() << " clusters"<< endl;
        //sort candidate regions for likelyhood of an alignment
        sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
        //for every interval, try to align
        int alnCount = 0;
        size_t lisIndex = 0;
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
            long chrStart, chrEnd;
            sa.getChromBounds(matches[begin].ref, chrStart, chrEnd);
            if(print) cerr << "try cluster " << lisIndex << " with " << (end-begin+1) << " seeds." << endl;
            alignment_t * alignment = extendAlignment(dp_, sa.S, P, matches, begin, end, editDist, alnOptions, chrStart, chrEnd);
            if(alignment!=NULL){
                if(!fwStrand){
                    alignment->flag.set(4,true);
                    if(alnOptions.unique && !read.alignments.empty() && alignment->editDist < read.alignments[0]->editDist)
                        read.removeLastAlignment();
                }
                alnCount++;
                read.alignments.push_back(alignment);
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
        if(print){ 
            cerr << "stopped matching because "; 
            if(alnCount >= max(1,alnOptions.alignmentCount/2)) cerr << " enough alignments were found."<< endl;
            else if(lisIndex >= lisIntervals.size()) cerr << " no more clusters are available."<< endl;
            else if(lisIntervals[lisIndex].len <= (Plength*alnOptions.minCoverage)/100) cerr << " no new clusters reach minimum query coverage."<< endl;
            else if(trial >= alnOptions.maxTrial) cerr << " max number of extensions without result have been reached."<< endl;
        }
    }
}

void unpairedMatch(const sparseSA& sa, dynProg& dp_, read_t & read,const align_opt & alnOptions){
    bool print = (alnOptions.print>0);
    string P = read.sequence;
    string Prc = read.rcSequence;
    long Plength = (long) P.length();
    long editDist = (long)(alnOptions.errorPercent*Plength)+1;
    int min_len = alnOptions.minMemLength;
    if(Plength/editDist > min_len){
        min_len = Plength/editDist;
        if(print) cerr << "min length adjusted to " << min_len << endl;
    }
    vector<match_t> matches;
    vector<match_t> matchesRC;
    //calc seeds
    if(print) cerr << "calculate seeds" << endl;
    if(!alnOptions.noFW){
        calculateSeeds(sa, P, min_len, alnOptions.maxSeedCandidates, matches, alnOptions.tryHarder, alnOptions.memType);
    }
    if(print) cerr << "found " << matches.size() << " seeds on forward direction"<< endl;
    if(!alnOptions.noRC){
        calculateSeeds(sa, Prc, min_len, alnOptions.maxSeedCandidates, matchesRC, alnOptions.tryHarder, alnOptions.memType);
    }
    if(print) cerr << "found " << matchesRC.size() << " seeds on reverse direction"<< endl;
    //sort matches
    if(matches.size()>0 || matchesRC.size()>0){
        /////////////////////////
        //FIND CANDIDATE REGIONS
        /////////////////////////
        vector<lis_t> lisIntervals;
        if(print) cerr << "calculate clusters" << endl;
        calculateLISintervals(matches, true, Plength, editDist, lisIntervals);
        if(print) cerr << "found " << lisIntervals.size() << " clusters forward"<< endl;
        calculateLISintervals(matchesRC, false, Plength, editDist, lisIntervals);
        if(print) cerr << "found " << lisIntervals.size() << " clusters TOTAL"<< endl;
        //sort candidate regions for likelyhood of an alignment
        sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
        //for every interval, try to align
        int alnCount = 0;
        size_t lisIndex = 0;
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
            long chrStart, chrEnd;
            sa.getChromBounds(matchVector->at(begin).ref, chrStart, chrEnd);
            if(print) cerr << "try cluster " << lisIndex << " with " << (end-begin+1) << " seeds." << endl;
            alignment_t * alignment = extendAlignment(dp_, sa.S, lisIntervals[lisIndex].fw ? P : Prc, *matchVector, begin, end, editDist, alnOptions, chrStart, chrEnd);
            if(alignment!=NULL){
                if(!lisIntervals[lisIndex].fw)
                    alignment->flag.set(4,true);
                alnCount++;
                read.alignments.push_back(alignment);
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
        if(print){ 
            cerr << "stopped matching because "; 
            if(alnCount >= max(1,alnOptions.alignmentCount)) cerr << " enough alignments were found."<< endl;
            else if(lisIndex >= lisIntervals.size()) cerr << " no more clusters are available."<< endl;
            else if(lisIntervals[lisIndex].len <= (Plength*alnOptions.minCoverage)/100) cerr << " no new clusters reach minimum query coverage."<< endl;
            else if(trial >= alnOptions.maxTrial) cerr << " max number of extensions without result have been reached."<< endl;
        }
    }
}

//PAIRED END OPTIONS
bool isConcordant(const alignment_t* mate1, const alignment_t* mate2, const paired_opt& options){
    bool firstLeft = true;
    //Check same chromosome
    if(mate1->rname != mate2->rname)
        return false;
    //Check strand
    if(options.orientation == PAIR_FF){
        if(mate1->flag.test(4)!=mate2->flag.test(4))
            return false;
        else
            firstLeft = !mate1->flag.test(4);//if rc must be right
    }
    else if(options.orientation == PAIR_FR){
        if(mate1->flag.test(4)==mate2->flag.test(4))
            return false;
        else
            firstLeft = !mate1->flag.test(4);//if rc must be right
    }
    else{
        if(mate1->flag.test(4)==mate2->flag.test(4))
            return false;
        else
            firstLeft = mate1->flag.test(4);//if rc must be left
    }
    //check fragment insert size
    long begin = min(mate1->pos,mate2->pos);
    long end = max(mate1->pos+(long)mate1->refLength-1L,mate2->pos+(long)mate2->refLength-1L);
    long size = end - begin +1L;
    if(size < options.minInsert || size > options.maxInsert)
        return false;
    long begin1 = mate1->pos;
    long end1 = begin1 + (long)mate1->refLength-1L;
    long begin2 = mate2->pos;
    long end2 = begin2 + (long)mate2->refLength-1L;
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

void setUnPaired(alignment_t * toset, read_t& read, bool upstream){
    toset->flag.set(0,true);
    toset->flag.set(6,upstream);
    toset->flag.set(7,!upstream);
    toset->flag.set(3,true);
    read.pairedAlignmentCount++;
}

void setPaired(alignment_t* mate1, alignment_t* mate2, read_t& upstream, read_t& downstream, bool concordant){
    mate1->flag.set(0, true);
    mate2->flag.set(0, true);
    mate1->addMate(mate2, concordant, true);
    mate2->addMate(mate1, concordant, false);
    upstream.pairedAlignmentCount++;
    downstream.pairedAlignmentCount++;
}

//aln both and pair
void pairedMatch1(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    //Calculate the mappings for read 1 and 2
    bool print = alnOptions.print>0;
    if(print) cerr << "matching first mate ..." << endl;
    unpairedMatch(sa, dp_, mate1, alnOptions);
    if(print) cerr << "first mate resulted in " << mate1.alignments.size() << " alignments" << endl;
    if(print) cerr << "matching second mate ..." << endl;
    unpairedMatch(sa, dp_, mate2, alnOptions);
    if(print) cerr << "second mate resulted in " << mate2.alignments.size() << " alignments" << endl;
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
        int maxScoreFirst = mate1.alignments[0]->alignmentScore;
        if(print) cerr << "max score mate 1 is " << maxScoreFirst << endl;
        int maxScoreSecond = mate2.alignments[0]->alignmentScore;
        if(print) cerr << "max score mate 2 is " << maxScoreSecond << endl;
        for(int j = 0; j < alnCount2; j++){
            mate2.alignments[j]->setLocalPos(sa);
        }
        //search concordant alignments
        if(print) cerr << "search concordant alignments" << endl;
        while(i < alnCount1 && alnCount < alnOptions.alignmentCount){
            alignment_t * aln1 = mate1.alignments[i];
            aln1->setLocalPos(sa);
            int j = 0;
            while(j < alnCount2 && alnCount < alnOptions.alignmentCount){
                alignment_t * aln2 = mate2.alignments[j];
                if(isConcordant(aln1, aln2, pairedOpt)){
                    if(print) cerr << "found concordancy between aln " << i << " of mate 1 and aln " << j << " of mate 2" << endl;
                    //Set the concordants and add another if necessary
                    setPaired(aln1,aln2,mate1,mate2,true);
                    concordant++;
                    alnCount++;
                }
                else if(pairedOpt.discordant && maxScoreFirst == aln1->alignmentScore && maxScoreSecond == aln2->alignmentScore){
                    if(print) cerr << "found discordant pair between aln " << i << " of mate 1 and aln " << j << " of mate 2" << endl;
                    discordantAln.push_back(pair_t(i,j));
                    discordant++;
                }
                j++;
            }
            i++;
        }
        //search discordant alignments
        if(pairedOpt.discordant){
            size_t j = 0;
            while(j < discordantAln.size() && alnCount < alnOptions.alignmentCount){
                alignment_t * aln1 = mate1.alignments[discordantAln[j].mate1];
                alignment_t * aln2 = mate2.alignments[discordantAln[j].mate2];
                setPaired(aln1,aln2,mate1,mate2,false);
                j++;
                alnCount++;
            }
        }
    }
    //sort alignments such that first concordant, than discordant, than unpaired
    if(pairedOpt.mixed){
        i = 0;
        int alnCountFirst = alnCount;
        while(i < alnCount1 && alnCountFirst < alnOptions.alignmentCount){
            if(!mate1.alignments[i]->flag.test(0)){
                if(print) cerr << "found unpaired aln " << i << " of mate 1" << endl;
                setUnPaired(mate1.alignments[i],mate1, true);
                alnCountFirst++;
            }
            i++;
        }
        mate1.pairedAlignmentCount = alnCountFirst;
        i = 0;
        int alnCountSecond = alnCount;
        while(i < alnCount2 && alnCountSecond < alnOptions.alignmentCount){
            if(!mate2.alignments[i]->flag.test(0)){
                if(print) cerr << "found unpaired aln " << i << " of mate 2" << endl;
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
    long chrStart, chrEnd;
    sa.getChromBounds(matchVector->at(begin).ref, chrStart, chrEnd);
    alignment_t * alignment = extendAlignment(dp_, sa.S, lis.fw ? read.sequence : read.rcSequence, *matchVector, begin, end, editDist, alnOptions, chrStart, chrEnd);
    if(alignment!=NULL){
        lis.extended = true;
        if(!lis.fw){
            alignment->flag.set(4,true);
        }
        alignment->setLocalPos(sa);
        lis.alignment = alignment;
        if(alnOptions.print) cerr << "alignment found" << endl;
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
    bool print = alnOptions.print > 0;
    long M1length = (long)mate1.sequence.size();
    long M2length = (long)mate2.sequence.size();
    long editDistM1 = (long)(alnOptions.errorPercent*M1length)+1;
    long editDistM2 = (long)(alnOptions.errorPercent*M2length)+1;
    int state = 0;
    size_t i = 0;
    size_t j = 0;
    size_t begin = j;
    while(i < lisIntervalsFM1.size() && concordant < alnOptions.alignmentCount){
        state = 0;
        begin = j;
        lis_t & lis1 = lisIntervalsFM1[i];
        while(j < lisIntervalsFM2.size() && concordant < alnOptions.alignmentCount && state < 2){
            lis_t & lis2 = lisIntervalsFM2[j];
            if(print) cerr << "check concordancy of cluster " << i << " and " << j << endl;
            if(isPosConcordant(lis1, lis2, editDistM1, editDistM2, M1length, M2length, pairedOpt)){//concordant (pos)
                if(print) cerr << "cluster " << i << " and " << j << " could be concordant: check if they have good aln" << endl;
                if(state==0){ state++; begin = j;}//first pos where a concordant alignment was found
                //Align if necessary, 
                bool firstAligned = alignFromLIS(sa, dp_, mate1, lis1, editDistM1, alnOptions);
                if(print) cerr << "mate 1 aligned? " << firstAligned << endl;
                bool secondAligned = alignFromLIS(sa, dp_, mate2, lis2, editDistM2, alnOptions);
                if(print) cerr << "mate 2 aligned? " << secondAligned << endl;
                //If both align: check pairing and set conc/disc: add to pairedInfo
                if(firstAligned && secondAligned){
                    if(print) cerr << "final check for concordancy" << endl;
                    if(isConcordant(lis1.alignment, lis2.alignment, pairedOpt)){//concordant
                        concordant++;
                        if(print) cerr << "concordant" << endl;
                        lis1.len = 0; lis2.len = 0;
                        setPaired(lis1.alignment, lis2.alignment,mate1,mate2,true);
                    }
                    else if(pairedOpt.discordant){
                        lis1.len = 0; lis2.len = 0;
                        discordant++;
                        if(print) cerr << "discordant" << endl;
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
    bool print = alnOptions.print>0;
    long editDist = (long)(alnOptions.errorPercent*read.sequence.size())+1;
    //sort candidate regions for likelyhood of an alignment
    sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
    //for every interval, try to align
    size_t lisIndex = 0;
    //trial parameter for performance reasons
    int trial = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    if(print) cerr << "start extending clusters" << endl;
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
        long chrStart, chrEnd;
        sa.getChromBounds(matchVector->at(begin).ref, chrStart, chrEnd);
        if(print) cerr << "try cluster " << lisIndex << " of " << lisIntervals.size() << " with " << (end-begin+1) << " seeds" << endl;
        alignment_t * alignment = extendAlignment(dp_, sa.S, lisIntervals[lisIndex].fw ? read.sequence : read.rcSequence, 
                *matchVector, begin, end, editDist, alnOptions, chrStart, chrEnd);
        if(alignment!=NULL){
            if(print) cerr << "extension succes" << endl;
            if(!lisIntervals[lisIndex].fw){
                alignment->flag.set(4,true);
            }
            alnCount++;
            alignment->setLocalPos(sa);
            lisIntervals[lisIndex].extended = true;
            lisIntervals[lisIndex].alignment = alignment;
            trial = 0;
        }
        lisIndex++;
        trial++;
    }
    if(print){
        cerr << "stopped extending clusters because ";
        if(alnCount >= alnOptions.alignmentCount) cerr << " enough alignments were found."<< endl;
        else if(lisIndex >= lisIntervals.size()) cerr << " no more clusters are available."<< endl;
        else if(lisIntervals[lisIndex].len <= (read.sequence.size()*alnOptions.minCoverage)/100) cerr << " no new clusters reach minimum query coverage."<< endl;
        else if(trial >= alnOptions.maxTrial) cerr << " max number of extensions without result have been reached."<< endl;
    }
}

//Calculate LIS of both mates, match together 2 LIS, if match, calculate both alignments
void pairedMatch3(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    bool print = alnOptions.print>0;
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
    if(print) cerr << "calculate seeds for mate 1 and 2 in first direction" << endl;
    calculateSeeds(sa, mate1FWfirst ? mate1.sequence : mate1.rcSequence, min_len, alnOptions.maxSeedCandidates, firstmatchesM1, alnOptions.tryHarder, alnOptions.memType);
    calculateSeeds(sa, mate2FWfirst ? mate2.sequence : mate2.rcSequence, min_len, alnOptions.maxSeedCandidates, firstmatchesM2, alnOptions.tryHarder, alnOptions.memType);
    if(print) cerr << "found " << firstmatchesM1.size() << " seeds for mate 1 and " << firstmatchesM2.size() << " seeds for mate 2" << endl;
    if(firstmatchesM1.size()>0 && firstmatchesM2.size()>0){
        //Calc LIS intervals
        //lis intervals are intervals in seed matches that are ordered according to reference offset
        calculateLISintervals(firstmatchesM1, mate1FWfirst, M1length, editDistM1, lisIntervalsFM1);
        calculateLISintervals(firstmatchesM2, mate2FWfirst, M2length, editDistM2, lisIntervalsFM2);
        //find concordant intervals
        if(print) cerr << "try to connect " << lisIntervalsFM1.size() << " clusters of mate 1 with " << lisIntervalsFM2.size() << " clusters of mate 2" << endl;
        coupleLIS(sa, dp_, mate1, mate2, lisIntervalsFM1, lisIntervalsFM2, concordant, discordant, alnOptions, pairedOpt);
    }
    if(print) cerr << "found " << concordant << " concordant alignments" << endl;
    //Do the same for other direction if necessary
    if(concordant < alnOptions.alignmentCount){
        if(print) cerr << "calculate seeds for mate 1 and 2 in second direction" << endl;
        calculateSeeds(sa, mate1FWfirst ? mate1.rcSequence : mate1.sequence, min_len, alnOptions.maxSeedCandidates, secondmatchesM1, alnOptions.tryHarder, alnOptions.memType);
        calculateSeeds(sa, mate2FWfirst ? mate2.rcSequence : mate2.sequence, min_len, alnOptions.maxSeedCandidates, secondmatchesM2, alnOptions.tryHarder, alnOptions.memType);
        if(print) cerr << "found " << secondmatchesM1.size() << " seeds for mate 1 and " << secondmatchesM2.size() << " seeds for mate 2" << endl;
        if(secondmatchesM1.size()>0 && secondmatchesM2.size()>0){
            //Calc LIS intervals
            //lis intervals are intervals in seed matches that are ordered according to reference offset
            calculateLISintervals(secondmatchesM1, !mate1FWfirst, M1length, editDistM1, lisIntervalsSM1);
            calculateLISintervals(secondmatchesM2, !mate2FWfirst, M2length, editDistM2, lisIntervalsSM2);
            //find concordant intervals
            if(print) cerr << "try to connect " << lisIntervalsSM1.size() << " clusters of mate 1 with " << lisIntervalsSM2.size() << " clusters of mate 2" << endl;
            coupleLIS(sa, dp_, mate1, mate2, lisIntervalsSM1, lisIntervalsSM2, concordant, discordant, alnOptions, pairedOpt);
        }
    }
    concordant = concordant + discordant;
    int maxScoreFirst;
    int maxScoreSecond;
    if(print) cerr << "found " << concordant << " concordant alignments so far" << endl;
    //preprocessing for discordant and mixed: calculate more alignments
    if(concordant < alnOptions.alignmentCount && (pairedOpt.discordant || pairedOpt.mixed)){
        cerr << "try to find discordant and/or unpaired alignments" << endl;
        //concatenate the lists per mate
        lisIntervalsFM1.insert(lisIntervalsFM1.end(),lisIntervalsSM1.begin(),lisIntervalsSM1.end());
        lisIntervalsFM2.insert(lisIntervalsFM2.end(),lisIntervalsSM2.begin(),lisIntervalsSM2.end());
        //try to align some more LIS intervals
        cerr << "extend more clusters to full alignments" << endl;
        unpairedMatchFromLIS(sa, dp_, mate1, lisIntervalsFM1, alnOptions, concordant);
        unpairedMatchFromLIS(sa, dp_, mate2, lisIntervalsFM2, alnOptions, concordant);
        //sort according to alignmentScore and globalPos
        cerr << "sort clusters according to alignment score" << endl;
        sort(lisIntervalsFM1.begin(),lisIntervalsFM1.end(), compAlignmentAndScore);
        sort(lisIntervalsFM2.begin(),lisIntervalsFM2.end(), compAlignmentAndScore);
        maxScoreFirst = lisIntervalsFM1[0].alignment->alignmentScore;
        maxScoreSecond = lisIntervalsFM2[0].alignment->alignmentScore;
            //Set discordant to LIS, count something else, and make sure we can retrieve the right ones
        if(pairedOpt.discordant){
            //extra discordant searches depending on max score
            size_t i=0;
            while(i < lisIntervalsFM1.size() && 
                    concordant < alnOptions.alignmentCount && 
                    lisIntervalsFM1[i].extended &&
                    lisIntervalsFM1[i].alignment->alignmentScore == maxScoreFirst){
                set<long> globPosMate;
                if(lisIntervalsFM1[i].alignment->paired())
                    for(int idx=0; idx < lisIntervalsFM1[i].alignment->pairedCount(); idx++)
                        globPosMate.insert(lisIntervalsFM1[i].alignment->mateInfo[idx].pnextGlob);
                size_t j = 0;
                while(j < lisIntervalsFM2.size() && 
                        concordant < alnOptions.alignmentCount &&
                        lisIntervalsFM2[j].extended &&
                        lisIntervalsFM2[j].alignment->alignmentScore == maxScoreSecond){
                    if(globPosMate.count(lisIntervalsFM2[j].alignment->globPos)==0){
                        //found extra discordant
                        concordant++;
                        cerr << "discordant pair found between aln " << i << " of mate 1 and aln " << j << " of mate 2" << endl;
                        setPaired(lisIntervalsFM1[i].alignment,lisIntervalsFM2[j].alignment,mate1,mate2, false);
                    }
                    j++;
                }
                i++;
            }
        }
        if(pairedOpt.mixed && concordant < alnOptions.alignmentCount ){
            //extra unpaired alignments
            cerr << "search for unpaired aln" << endl;
            size_t i = 0;
            int alnCountFirst = concordant;
            while(i < lisIntervalsFM1.size() && lisIntervalsFM1[i].extended && alnCountFirst < alnOptions.alignmentCount){
                if(!lisIntervalsFM1[i].alignment->paired()){
                    cerr << "unpaired aln " << i << " of mate 1 found" << endl;
                    setUnPaired(lisIntervalsFM1[i].alignment,mate1, true);
                    alnCountFirst++;
                }
                i++;
            }
            i = 0;
            int alnCountSecond = concordant;
            while(i < lisIntervalsFM2.size() && lisIntervalsFM2[i].extended && alnCountSecond < alnOptions.alignmentCount){
                if(!lisIntervalsFM2[i].alignment->paired()){
                    cerr << "unpaired aln " << i << " of mate 1 found" << endl;
                    setUnPaired(lisIntervalsFM2[i].alignment,mate2, false);
                    alnCountSecond++;
                }
                i++;
            }
        }
    }
    for(size_t i=0; i < lisIntervalsFM1.size(); i++ )
        if(lisIntervalsFM1[i].extended)
            mate1.alignments.push_back(lisIntervalsFM1[i].alignment);
    for(size_t i=0; i < lisIntervalsFM2.size(); i++ )
        if(lisIntervalsFM2[i].extended)
            mate2.alignments.push_back(lisIntervalsFM2[i].alignment);
}

bool isConcordantAlnToLis(const alignment_t* mate1, const lis_t& mate2, bool mate1isFirst, long editDist2, long M2length, const paired_opt& options){
    bool firstLeft = true;
    bool mate1FW = !mate1->flag.test(4);
    bool mate2FW = mate2.fw;
    //Check strand
    if(options.orientation == PAIR_FF || options.orientation == PAIR_FR){
        firstLeft = mate1isFirst ? mate1FW : mate2FW;//if rc must be right
    }
    else{
        firstLeft = mate1isFirst ? !mate1FW : !mate2FW;//if rc must be left
    }
    //find possible offset
    long posLeft1 = mate1->globPos;
    long posLeft2 = mate2.matches->at(mate2.begin).ref-mate2.matches->at(mate2.begin).query;
    long posRight1 = posLeft1 + (long) mate1->refLength;
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
    bool print = alnOptions.print>0;
    long M1length = (long) mate1.sequence.size();
    long M2length = (long) mate2.sequence.size();
    long editDistM1 = (long)(alnOptions.errorPercent*M1length)+1;
    long editDistM2 = (long)(alnOptions.errorPercent*M2length)+1;
    //Matching strand and mate: sort acocrding to score
    sort(lisIntervalsFM1.begin(), lisIntervalsFM1.end(), compIntervals);
    //Other strand: for pairing
    sort(lisIntervalsFM2.begin(), lisIntervalsFM2.end(), compIntervalsRef);
    size_t lisIndex = 0;
    //trial parameter for performance reasons
    int trial = 0;
    //////////////////////
    //MAIN ALIGNMENT LOOP
    //////////////////////
    if(print) cerr << "start extension loop" << endl;
    while(concordant < alnOptions.alignmentCount && 
            lisIndex < lisIntervalsFM1.size() && 
            lisIntervalsFM1[lisIndex].len > (M1length*alnOptions.minCoverage)/100 && 
            trial < alnOptions.maxTrial){
        //Alignment of mate1
        if(print) cerr << "try cluster " << lisIndex << endl;
        if(!lisIntervalsFM1[lisIndex].extended){
            if(print) cerr << "try to extend to alignment" << endl;
            //sort matches in this interval according to query position
            int begin = lisIntervalsFM1[lisIndex].begin;
            int end = lisIntervalsFM1[lisIndex].end;
            vector<match_t>* matchVector = lisIntervalsFM1[lisIndex].matches;
            //sort this candidate region by query position
            sort(matchVector->begin()+begin,matchVector->begin()+end+1, compMatchesQuery);
            long chrStart, chrEnd;
            sa.getChromBounds(matchVector->at(begin).ref, chrStart, chrEnd);
            alignment_t * alignment = extendAlignment(dp_, sa.S, lisIntervalsFM1[lisIndex].fw ? mate1.sequence : mate1.rcSequence, 
                    *matchVector, begin, end, editDistM1, alnOptions, chrStart, chrEnd);
            if(alignment!=NULL){
                if(print) cerr << "cluster succesfully extended to an alignment" << endl;
                lisIntervalsFM1[lisIndex].extended = true;
                if(!lisIntervalsFM1[lisIndex].fw)
                    alignment->flag.set(4,true);
                alignment->setLocalPos(sa);
                lisIntervalsFM1[lisIndex].alignment = alignment;
            }
        }
        else{
            if(print) cerr << "this cluster was already extended" << endl;
        }
        //If an alignment, continue to find match for other mate
        if(lisIntervalsFM1[lisIndex].extended){
            if(print) cerr << "try to find concordant cluster for this clusters alignment" << endl;
            trial = 0;
            lis_t & lis1 = lisIntervalsFM1[lisIndex];
            int state = 0;
            int j = 0;
            while(j < lisIntervalsFM2.size() && concordant < alnOptions.alignmentCount && state < 2){
                lis_t & lis2 = lisIntervalsFM2[j];
                if(print) cerr << "try cluster " << j << " which is extended (0/1) " << lis2.extended << endl;
                if((lis2.extended && isConcordant(mate1isFirst ? lis1.alignment : lis2.alignment, 
                        mate1isFirst ? lis2.alignment : lis1.alignment, pairedOpt)) ){
                    state = 1;
                    if(mate1isFirst){//if not: repeated find
                        concordant++;
                        setPaired(lis1.alignment,lis2.alignment,mate1,mate2, true);
                        if(print) cerr << "concordant alignment found" << endl;
                    }
                }
                else if(!lis2.extended && isConcordantAlnToLis(lis1.alignment, lis2, mate1isFirst, editDistM2, M2length, pairedOpt)){
                    state=1;
                    if(alignFromLIS(sa, dp_, mate2, lis2, editDistM2, alnOptions)){
                        if(isConcordant(mate1isFirst ? lis1.alignment : lis2.alignment, 
                                mate1isFirst ? lis2.alignment : lis1.alignment, pairedOpt)){//concordant
                            concordant++;
                            if(print) cerr << "concordant alignment found" << endl;
                            mate1isFirst ? setPaired(lis1.alignment,lis2.alignment,mate1,mate2 , true) : 
                                setPaired(lis2.alignment,lis1.alignment,mate2,mate1 , true);
                        }
                        else if(pairedOpt.discordant){
                            discordant++;
                            if(print) cerr << "discordant alignment found" << endl;
                            mate1isFirst ? setPaired(lis1.alignment,lis2.alignment,mate1,mate2 , false) : 
                                setPaired(lis2.alignment,lis1.alignment,mate2,mate1 , false);
                        }
                    }
                }
                else if(state==1)//discordant, move on
                    state++;
                j++;
            }
            if(print){
                cerr << "stopped searching for matching alignment because ";
                if(concordant >= alnOptions.alignmentCount) cerr << " enough alignments were found."<< endl;
                else if(j >= lisIntervalsFM2.size()) cerr << " no more clusters are available."<< endl;
                else if(state == 2) cerr << " concordant aln were found and remaining clusters are discordant"<< endl;
            }
        }
        lisIndex++;
        trial++;
    }
    if(print){
        cerr << "stopped extending clusters because ";
        if(concordant >= alnOptions.alignmentCount) cerr << " enough alignments were found."<< endl;
        else if(lisIndex >= lisIntervalsFM1.size()) cerr << " no more clusters are available."<< endl;
        else if(lisIntervalsFM1[lisIndex].len <= (M1length*alnOptions.minCoverage)/100) cerr << " no new clusters reach minimum query coverage."<< endl;
        else if(trial >= alnOptions.maxTrial) cerr << " max number of extensions without result have been reached."<< endl;
    }
}

//Calculate LIS of both mates, align 1, match LIS to other, align 2
void pairedMatch4(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    bool print = alnOptions.print>0;
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
    if(print) cerr << "calculate seeds for mate 1 and 2 in first direction" << endl;
    calculateSeeds(sa, mate1FWfirst ? mate1.sequence : mate1.rcSequence, min_len, alnOptions.maxSeedCandidates, firstmatchesM1, alnOptions.tryHarder, alnOptions.memType);
    calculateSeeds(sa, mate2FWfirst ? mate2.sequence : mate2.rcSequence, min_len, alnOptions.maxSeedCandidates, firstmatchesM2, alnOptions.tryHarder, alnOptions.memType);
    if(print) cerr << "found " << firstmatchesM1.size() << " seeds for mate 1 and " << firstmatchesM2.size() << " seeds for mate 2" << endl;
    if(firstmatchesM1.size()>0 && firstmatchesM2.size()>0){
        //Calc LIS intervals
        //lis intervals are intervals in seed matches that are ordered according to reference offset
        calculateLISintervals(firstmatchesM1, mate1FWfirst, M1length, editDistM1, lisIntervalsFM1);
        calculateLISintervals(firstmatchesM2, mate2FWfirst, M2length, editDistM2, lisIntervalsFM2);
        //find concordant intervals
        if(print) cerr << "try to expand " << lisIntervalsFM1.size() << " clusters of mate 1 and connect them to " << lisIntervalsFM2.size() << " clusters of mate 2" << endl;
        matchStrandOfPair(sa, dp_, mate1, mate2, lisIntervalsFM1, lisIntervalsFM2, concordant, discordant, true, alnOptions, pairedOpt);
        if(print) cerr << "found " << concordant << " concordant alignments" << endl;
    }
    //Do the same for other direction if necessary
    if(concordant < alnOptions.alignmentCount){
        if(print) cerr << "calculate seeds for mate 1 and 2 in second direction" << endl;
        calculateSeeds(sa, mate1FWfirst ? mate1.rcSequence : mate1.sequence, min_len, alnOptions.maxSeedCandidates, secondmatchesM1, alnOptions.tryHarder, alnOptions.memType);
        calculateSeeds(sa, mate2FWfirst ? mate2.rcSequence : mate2.sequence, min_len, alnOptions.maxSeedCandidates, secondmatchesM2, alnOptions.tryHarder, alnOptions.memType);
        if(print) cerr << "found " << secondmatchesM1.size() << " seeds for mate 1 and " << secondmatchesM2.size() << " seeds for mate 2" << endl;
        if(secondmatchesM1.size()>0 && secondmatchesM2.size()>0){
            //Calc LIS intervals
            //lis intervals are intervals in seed matches that are ordered according to reference offset
            calculateLISintervals(secondmatchesM1, !mate1FWfirst, M1length, editDistM1, lisIntervalsSM1);
            calculateLISintervals(secondmatchesM2, !mate2FWfirst, M2length, editDistM2, lisIntervalsSM2);
            if(print) cerr << "try to expand " << lisIntervalsSM1.size() << " clusters of mate 1 (other dir) and connect them to " << lisIntervalsSM2.size() << " clusters of mate 2" << endl;
            //find concordant intervals
            matchStrandOfPair(sa, dp_, mate1, mate2, lisIntervalsSM1, lisIntervalsSM2, concordant, discordant, true, alnOptions, pairedOpt);
            if(print) cerr << "found " << concordant << " concordant alignments" << endl;
        }
    }
    //Do the same for first direction of other mate
    if(concordant < alnOptions.alignmentCount){
        if(firstmatchesM1.size()>0 && firstmatchesM2.size()>0){
            //find concordant intervals
            if(print) cerr << "try to expand " << lisIntervalsFM2.size() << " clusters of mate 2 and connect them to " << lisIntervalsFM1.size() << " clusters of mate 1" << endl;
            matchStrandOfPair(sa, dp_, mate2, mate1, lisIntervalsFM2, lisIntervalsFM1, concordant, discordant, false, alnOptions, pairedOpt);
            if(print) cerr << "found " << concordant << " concordant alignments" << endl;
        }
    }
    //Do the same for second direction of other mate
    if(concordant < alnOptions.alignmentCount){
        if(secondmatchesM1.size()>0 && secondmatchesM2.size()>0){
            //find concordant intervals
            if(print) cerr << "try to expand " << lisIntervalsSM2.size() << " clusters of mate 2 and connect them to " << lisIntervalsSM1.size() << " clusters of mate 1" << endl;
            matchStrandOfPair(sa, dp_, mate2, mate1, lisIntervalsSM2, lisIntervalsSM1, concordant, discordant, false, alnOptions, pairedOpt);
            if(print) cerr << "found " << concordant << " concordant alignments" << endl;
        }
    }
    concordant = concordant + discordant;
    int maxScoreFirst;
    int maxScoreSecond;
    //preprocessing for discordant and mixed: calculate more alignments
    if(concordant < alnOptions.alignmentCount && (pairedOpt.discordant || pairedOpt.mixed)){
        if(print) cerr << "search for discordant and/or unpaired alignments" << endl;
        //concatenate the lists per mate
        lisIntervalsFM1.insert(lisIntervalsFM1.end(),lisIntervalsSM1.begin(),lisIntervalsSM1.end());
        lisIntervalsFM2.insert(lisIntervalsFM2.end(),lisIntervalsSM2.begin(),lisIntervalsSM2.end());
        //sort according to alignmentScore and globalPos
        if(print) cerr << "sort aln according to alignment score and global pos" << endl;
        sort(lisIntervalsFM1.begin(),lisIntervalsFM1.end(), compAlignmentAndScore);
        sort(lisIntervalsFM2.begin(),lisIntervalsFM2.end(), compAlignmentAndScore);
        maxScoreFirst = lisIntervalsFM1[0].alignment->alignmentScore;
        maxScoreSecond = lisIntervalsFM2[0].alignment->alignmentScore;
        if(pairedOpt.discordant){
            //extra discordant searches depending on max score
            size_t i=0;
            while(i < lisIntervalsFM1.size() && 
                    concordant < alnOptions.alignmentCount && 
                    lisIntervalsFM1[i].extended &&
                    lisIntervalsFM1[i].alignment->alignmentScore == maxScoreFirst){
                set<long> globPosMate;
                if(lisIntervalsFM1[i].alignment->paired())
                    for(int idx=0; idx < lisIntervalsFM1[i].alignment->pairedCount(); idx++)
                        globPosMate.insert(lisIntervalsFM1[i].alignment->mateInfo[idx].pnextGlob);
                size_t j = 0;
                while(j < lisIntervalsFM2.size() && 
                        concordant < alnOptions.alignmentCount &&
                        lisIntervalsFM2[j].extended &&
                        lisIntervalsFM2[j].alignment->alignmentScore == maxScoreSecond){
                    if(globPosMate.count(lisIntervalsFM2[j].alignment->globPos)==0){
                        //found extra discordant
                        concordant++;
                        if(print) cerr << "found discordant aln between aln " << i << " of mate 1 and aln " << j << " of mate 2" << endl;
                        setPaired(lisIntervalsFM1[i].alignment,lisIntervalsFM2[j].alignment,mate1,mate2, false);
                    }
                    j++;
                }
                i++;
            }
        }
        if(pairedOpt.mixed && concordant < alnOptions.alignmentCount ){
            //extra unpaired alignments
            size_t i = 0;
            int alnCountFirst = concordant;
            while(i < lisIntervalsFM1.size() && lisIntervalsFM1[i].extended && alnCountFirst < alnOptions.alignmentCount){
                if(!lisIntervalsFM1[i].alignment->paired()){
                    if(print) cerr << "found unpaired aln of mate 1: " << i << endl;
                    setUnPaired(lisIntervalsFM1[i].alignment,mate1, true);
                    alnCountFirst++;
                }
                i++;
            }
            i = 0;
            int alnCountSecond = concordant;
            while(i < lisIntervalsFM2.size() && lisIntervalsFM2[i].extended && alnCountSecond < alnOptions.alignmentCount){
                if(!lisIntervalsFM2[i].alignment->paired()){
                    if(print) cerr << "found unpaired aln of mate 1: " << i << endl;
                    setUnPaired(lisIntervalsFM2[i].alignment,mate2, false);
                    alnCountSecond++;
                }
                i++;
            }
        }
    }
    if(print) cerr << "add alignments to output data structures" << endl;
    for(size_t i=0; i < lisIntervalsFM1.size(); i++ )
        if(lisIntervalsFM1[i].extended)
            mate1.alignments.push_back(lisIntervalsFM1[i].alignment);
    for(size_t i=0; i < lisIntervalsFM2.size(); i++ )
        if(lisIntervalsFM2[i].extended)
            mate2.alignments.push_back(lisIntervalsFM2[i].alignment);
}

bool dpWindow(const alignment_t * mate1, bool mate1isFirst, long editDist2, long M2length, const paired_opt& options,
        boundaries& grenzen, bool& otherFW, int& bandLeft, int& bandRight){
    bool firstLeft = true;
    bool mate1FW = !mate1->flag.test(4);
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
    long leftMate = mate1->globPos;
    long rightMate = leftMate + (long) mate1->refLength-1L;
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
    bool print = alnOptions.print>0;
    read_t & base = alignFirstMate ? mate1 : mate2;
    read_t & other = alignFirstMate ? mate2 : mate1;
    string & P = forward ? base.sequence : base.rcSequence;
    long Plength = (long) P.size();
    long Olength = (long) other.sequence.size();
    long editDistBase = (long)(alnOptions.errorPercent*Plength)+1;
    long editDistOther = (long)(alnOptions.errorPercent*Olength)+1;
    int min_len = alnOptions.minMemLength;
    if(Plength/editDistBase > min_len)
        min_len = Plength/editDistBase;
    vector<match_t> matches;
    //calc seeds
    calculateSeeds(sa, P, min_len, alnOptions.maxSeedCandidates, matches, alnOptions.tryHarder, alnOptions.memType);
    //sort matches
    if(print) cerr << "found " << matches.size() << " seeds" << endl;
    if(matches.size()>0){
        vector<lis_t> lisIntervals;
        calculateLISintervals(matches, forward, Plength, editDistBase, lisIntervals);
        //sort candidate regions for likelyhood of an alignment
        sort(lisIntervals.begin(),lisIntervals.end(), compIntervals);
        size_t lisIndex = 0;
        int trial = 0;//paired succes times
        if(print) cerr << "found " << lisIntervals.size() << " clusters" << endl;
        if(print) cerr << "try to extend clusters ..." << endl;
        while(concordant < alnOptions.alignmentCount &&
                lisIndex < lisIntervals.size() &&
                lisIntervals[lisIndex].len > (Plength*alnOptions.minCoverage)/100 &&
                trial < alnOptions.maxTrial){
            if(print) cerr << "try cluster " << lisIndex << endl;
            //try to find out if we already aligned this LIS-interval when pairing the other mate
            int i = 0;
            while(i < other.alignmentCount() && isConcordantAlnToLis(other.alignments[i], lisIntervals[lisIndex], !alignFirstMate, editDistBase, Plength, pairedOpt))
                i++;
            if(alignFirstMate || i == other.alignmentCount()){
                //have to extend alignment
                int begin = lisIntervals[lisIndex].begin;
                int end = lisIntervals[lisIndex].end;
                sort(matches.begin()+begin,matches.begin()+end+1, compMatchesQuery);
                long chrStart, chrEnd;
                sa.getChromBounds(matches[begin].ref, chrStart, chrEnd);
                if(print) cerr << "try to extend cluster " << lisIndex << endl;
                alignment_t * alignment = extendAlignment(dp_, sa.S, P, matches, begin, end, editDistBase, alnOptions, chrStart, chrEnd);
                if(alignment!=NULL){
                    if(print) cerr << "cluster " << lisIndex << " yielded a good alignment" << endl;
                    alignment->setLocalPos(sa);
                    if(!forward)
                        alignment->flag.set(4,true);
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
                    if(print) cerr << "try to find alignment for other mate using DP in insertregion " << endl;
                    if(dpWindow(alignment, alignFirstMate, editDistOther, Olength, pairedOpt, 
                            grenzen, otherFW, bandLeft, bandRight))
                        dp_.dpBandFull( sa.S, otherFW ? other.sequence : other.rcSequence, grenzen, 
                            types, ERRORSTRING, output, bandLeft, bandRight, alnOptions.print>1);
                    else{
                        if(print) cerr << "unable to establish DP window" << endl;
                        output.editDist = 2*editDistOther;
                    }
                    if(output.editDist <= editDistOther){
                        //SUCCES!!!!!!!
                        if(print) cerr << "found concordant alignment" << endl;
                        alignment_t * mate = new alignment_t();
                        mate->cigarChars.insert(mate->cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
                        mate->cigarLengths.insert(mate->cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
                        mate->globPos = grenzen.refB+1;
                        mate->refLength = grenzen.refE-grenzen.refB+1;
                        mate->editDist = output.editDist;
                        mate->alignmentScore = output.dpScore;
                        if(!otherFW)
                            mate->flag.set(4,true);
                        mate->setLocalPos(sa);
                        other.alignments.push_back(mate);
                        //set paired: alignments are the ones last pushed
                        setPaired(mate1.alignments[mate1.alignmentCount()-1],mate2.alignments[mate2.alignmentCount()-1],mate1,mate2, true);
                        concordant++;
                        trial = 0;
                    }
                }
            }
            else{
                if(print) cerr << "cluset " << lisIndex << " was already extended in previous steps" << endl;
                trial = 0;
            }
            lisIndex++;
            trial++;
        }
        if(print){
            cerr << "stopped extending clusters because ";
            if(concordant >= alnOptions.alignmentCount) cerr << " enough alignments were found."<< endl;
            else if(lisIndex >= lisIntervals.size()) cerr << " no more clusters are available."<< endl;
            else if(lisIntervals[lisIndex].len <= (Plength*alnOptions.minCoverage)/100) cerr << " no new clusters reach minimum query coverage."<< endl;
            else if(trial >= alnOptions.maxTrial) cerr << " max number of extensions without result have been reached."<< endl;
        }
    }
}

//Calculate Mates 1 at a time, without calculating seeds for the other
//Optimalization: use 'Bailing' method: to enable/disable the calculation for the other mate
void pairedMatch2(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    bool print = alnOptions.print>0;
    int concordant = 0;
    if((alnOptions.noFW || alnOptions.noRC))
        return;//no possible alignments
    //fields
    bool mate1FWfirst = pairedOpt.orientation == PAIR_FR || pairedOpt.orientation == PAIR_FF;
    bool mate2FWfirst = pairedOpt.orientation == PAIR_RF || pairedOpt.orientation == PAIR_FF;
    //concordant
    if(print) cerr << "match mate 1 fwd as anchor for mate 2" << endl;
    pairedBowtie2(sa, dp_, mate1, mate2, true, mate1FWfirst, concordant, alnOptions, pairedOpt);
    if(print) cerr << "got " << concordant << " concordant alignments" << endl;
    if(concordant < alnOptions.alignmentCount){
        if(print) cerr << "match mate 1 rev as anchor for mate 2" << endl;
        pairedBowtie2(sa, dp_, mate1, mate2, true, !mate1FWfirst, concordant, alnOptions, pairedOpt);
        if(print) cerr << "got " << concordant << " concordant alignments" << endl;
    }
    if(concordant < alnOptions.alignmentCount){
        if(print) cerr << "match mate 2 fwd as anchor for mate 1" << endl;
        pairedBowtie2(sa, dp_, mate1, mate2, false, mate2FWfirst, concordant, alnOptions, pairedOpt);
        if(print) cerr << "got " << concordant << " concordant alignments" << endl;
    }
    if(concordant < alnOptions.alignmentCount){
        if(print) cerr << "match mate 2 rev as anchor for mate 1" << endl;
        pairedBowtie2(sa, dp_, mate1, mate2, false, !mate2FWfirst, concordant, alnOptions, pairedOpt);
        if(print) cerr << "got " << concordant << " concordant alignments" << endl;
    }
    //
    int maxScoreFirst;
    int maxScoreSecond;
    //preprocessing for discordant and mixed: calculate more alignments
    if(concordant < alnOptions.alignmentCount && (pairedOpt.discordant || pairedOpt.mixed)){
        if(print) cerr << "search for discordant and/or unpaired alignments" << endl;
        //sort according to alignmentScore and globalPos
        if(print) cerr << "sort aln according to alignment score" << endl;
        sort(mate1.alignments.begin(),mate1.alignments.end(), compAlignmentScore);
        sort(mate2.alignments.begin(),mate2.alignments.end(), compAlignmentScore);
        if(mate1.alignmentCount()>0 && mate2.alignmentCount() > 0){
            maxScoreFirst = mate1.alignments[0]->alignmentScore;
            maxScoreSecond = mate2.alignments[0]->alignmentScore;
            if(pairedOpt.discordant){
                //extra discordant searches depending on max score
                int i=0;
                while(i < mate1.alignmentCount() && 
                        concordant < alnOptions.alignmentCount && 
                        mate1.alignments[i]->alignmentScore == maxScoreFirst){
                    int j = 0;
                    while(j < mate2.alignmentCount() && 
                            concordant < alnOptions.alignmentCount &&
                            mate2.alignments[j]->alignmentScore == maxScoreSecond){
                        if((!mate1.alignments[i]->paired() || !mate2.alignments[j]->paired()) || 
                                mate1.alignments[i]->mateInfo[0].pnextGlob!=mate2.alignments[j]->globPos){
                            //found extra discordant
                            concordant++;
                            setPaired(mate1.alignments[i],mate2.alignments[j],mate1,mate2, false);
                            if(print) cerr << "found discordant aln between aln " << i << " of mate 1 and aln " << j << " of mate 2" << endl;
                        }
                        j++;
                    }
                    i++;
                }
            }
        }
    }
    if(pairedOpt.mixed && concordant < alnOptions.alignmentCount ){
        //extra unpaired alignments
        int i = 0;
        int alnCountFirst = concordant;
        while(i < mate1.alignmentCount() && alnCountFirst < alnOptions.alignmentCount){
            if(!mate1.alignments[i]->paired()){
                if(print) cerr << "found unpaired aln for mate 1: " << i << endl;
                setUnPaired(mate1.alignments[i],mate1, true);
                alnCountFirst++;
            }
            i++;
        }
        i = 0;
        int alnCountSecond = concordant;
        while(i < mate2.alignmentCount() && alnCountSecond < alnOptions.alignmentCount){
            if(!mate2.alignments[i]->paired()){
                if(print) cerr << "found unpaired aln for mate 1: " << i << endl;
                setUnPaired(mate2.alignments[i],mate2, false);
                alnCountSecond++;
            }
            i++;
        }
    }
}

//MODE 3 AND 4 NOT WORKING: [3] crash, [4] not pairing!!!
void pairedMatch(const sparseSA& sa, dynProg& dp_, read_t & mate1, read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt){
    if(pairedOpt.mode==1){
        pairedMatch1(sa, dp_, mate1, mate2, alnOptions, pairedOpt);
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

