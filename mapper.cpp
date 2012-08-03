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

void postProcess(vector<match_t> &matches){
    sort(matches.begin(),matches.end(), compMatches);
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
    if (queryRB + 1 < P.length() && refstrRB < S.length() - 1 && curEditDist <= editDist) {
        int refAlignRB = refstrRB + (P.length() - queryRB) + 1 + editDist - curEditDist;
        if (refAlignRB >= S.length())
            refAlignRB = S.length() - 1;
        boundaries grenzen(refstrRB + 1, refAlignRB, queryRB + 1, P.length() - 1);
        dp_type types;
        types.freeRefE = true;
        types.freeQueryE = clipping;
        dpBandStatic( S, P, grenzen, alnOptions.scores, types, ERRORSTRING, output, editDist-curEditDist+1, false);
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
        Utils::reverse_complement(P,false);
    int min_len = alnOptions.minMemLength;
    if(!alnOptions.fixedMinLength && Plength/editDist > min_len)
        min_len = Plength/editDist;
    vector<match_t> matches;
    //calc seeds
    sa.SMAM(P, matches, min_len, alnOptions.alignmentCount, false);
    if(alnOptions.tryHarder){//TODO: change try-harder to recalculate only after forward + reverse has been tried
        //easy solution: try reverse and if found: skip
        if(matches.empty())
            sa.SMAM(P, matches, min_len, 1000, false);
        if(matches.empty())
            sa.SMAM(P, matches, 20, 1000, false);
        if(matches.empty())
            sa.SMAM(P, matches, Plength/editDist, 1000, false);
    }
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
