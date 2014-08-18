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

#ifndef DP_H
#define	DP_H

#include <string>
#include <vector>
#include <sstream>

#define DEBUG 0

struct boundaries{
    boundaries(): refB(0), refE(0), queryB(0), queryE(0) {}
    boundaries(long refBegin, long refEnd, int queryBegin, int queryEnd): refB(refBegin), refE(refEnd), queryB(queryBegin), queryE(queryEnd){}
    long refB;
    long refE;
    int queryB;
    int queryE;
};

struct dp_scores{
    dp_scores(): match(0), mismatch(-1), openGap(0), extendGap(-1) {}
    dp_scores(int m, int mm, int oG, int eG): match(m),mismatch(mm), openGap(oG), extendGap(eG) {}
    int match;
    int mismatch;
    int openGap;
    int extendGap;
    int scoreMatrix[127][127];
    void updateScoreMatrixDna();
};

struct dp_type{
    dp_type(): freeRefB(false),freeRefE(false),freeQueryB(false),freeQueryE(false),local(false) {}
    dp_type(bool refB, bool refE, bool qB, bool qE, bool loc):freeRefB(refB),freeRefE(refE),freeQueryB(qB),freeQueryE(qE),local(loc) {}
    bool freeRefB;
    bool freeRefE;
    bool freeQueryB;
    bool freeQueryE;
    bool local;
};

struct dp_output{
    dp_output(): editDist(0),dpScore(0), cigarChars(0), cigarLengths(0) {}
    int editDist;
    int dpScore;
    std::vector<char> cigarChars;//static, fixed length field?
    std::vector<int> cigarLengths;
    void clear(){cigarChars.clear(); cigarLengths.clear(); editDist = dpScore = 0; };
};

enum outputType{
    DPSCORE,
    SCORE,
    ALIGNMENT,
    ERRORSTRING,
    ALL
};

struct dynProg {

//These fields should maybe be positioned somewhere else for multi-threading!!!
int ** M;
int ** UP;
int ** LEFT;
int DP_L2;
int DP_L1;
int L1;
int L2;
int bandSize;
int bandLeft;
int bandRight;
bool banded;
dp_scores& scores;

// Constructor builds sparse suffix array.
dynProg(int dimension, bool affine, dp_scores& scoreFunction);
~dynProg();

void initDPMatrix(int dimensionL2, int dimensionL1, bool affine);
void resizeDPMatrix(int dimensionL2, int dimensionL1, bool affine);
void deleteDPMatrix(bool affine);

int updateMatrix(const dp_type& type, bool backward, bool print);
int dpFillOptStatic(const std::string& ref,const std::string& query, bool forward, 
        const boundaries& offset, const dp_type& type);
int dpFillOptStaticBackward(const std::string& ref,const std::string& query,
        const boundaries& offset);
int dpFillStatic(const std::string& ref,const std::string& query, bool forward,
        const boundaries& offset, const dp_type& type);
int dpFillStaticBackward(const std::string& ref,const std::string& query,
        const boundaries& offset);
int findTraceBackPosStatic(int* const i, int* const j,const dp_type& type);
int dpTraceBackFull(int& i, int& j, dp_output& output, std::stringstream & ss, 
        const boundaries& offset, bool forward, const std::string& ref, const std::string& query);
int dpTraceBackOptFull(int& i, int& j, const dp_type& type, dp_output& output, 
        std::stringstream & ss, const boundaries& offset, bool forward, const std::string& ref, const std::string& query);
int dpTraceBackStatic(int& i, int& j, dp_output& output, std::stringstream & ss, 
        const boundaries& offset, bool forward, const std::string& ref, const std::string& query);
int dpTraceBackOptStatic(int& i, int& j, const dp_type& type, dp_output& output, 
        std::stringstream & ss, const boundaries& offset, bool forward, const std::string& ref, const std::string& query);

void  print_matrices(const std::string& ref, const std::string& query, const boundaries& offset, bool gapMatrices );
void  print_seq( const std::string& ref,const std::string& query, boundaries& offset );

int dpBandStatic( const std::string&, const std::string&, boundaries&, 
        const dp_type&, const outputType&, dp_output&, int minbandSize, int maxbandSize, bool print = false);

int dpBandFull( const std::string&, const std::string&, boundaries&, 
        const dp_type&, const outputType&, dp_output&, int bandLeft, int bandRight, bool print = false);

int dpTraceBackScore(int& i, int& j, const boundaries& offset, dp_output& output, 
        const std::string& ref, const std::string& query);

int dpBandScore( const std::string&, const std::string&, boundaries&, 
        const dp_type&, const outputType&, dp_output&, int bandLeft, int bandRight, bool print = false);

int dpBasic( const std::string&, const std::string&, boundaries&,
        const dp_type&, const outputType&, dp_output&);

};

#endif	/* DP_H */

