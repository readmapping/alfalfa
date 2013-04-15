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

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

#define DEBUG 0

using namespace std;

static const unsigned int ORDVALUE[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//0-9
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//10-19
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//20-29
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//30-39
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//40-49
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//50-59
                                         4, 4, 4, 4, 4, 0, 4, 1, 4, 4,//60-69 65:A, 67:C
                                         4, 2, 4, 4, 4, 4, 4, 4, 4, 4,//70-79 71:G
                                         4, 4, 4, 4, 3, 4, 4, 4, 4, 4,//80-89 84:T
                                         4, 4, 4, 4, 4, 4, 4, 0, 4, 1,//90-99 97:a, 99: c
                                         4, 4, 4, 2, 4, 4, 4, 4, 4, 4,//100-109 103:g
                                         4, 4, 4, 4, 4, 4, 3, 4, 4, 4,//110-119 116:t
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//120-129
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//130-139
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//140-149
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//150-159
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//160-169
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//170-179
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//180-189
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//190-199
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//200-209
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//210-219
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//220-229
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//230-239
                                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4,//240-249
                                         4, 4, 4, 4, 4, 4 };//250-255


struct boundaries{
    boundaries(): refB(0), refE(0), queryB(0), queryE(0) {}
    boundaries(long refBegin, long refEnd, long queryBegin, long queryEnd): refB(refBegin), refE(refEnd), queryB(queryBegin), queryE(queryEnd){}
    long refB;
    long refE;
    long queryB;
    long queryE;
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
    dp_output(): traceRef(""),traceQuery(""),editDist(0),dpScore(0), cigarChars(0), cigarLengths(0) {}
    string traceRef;//remove???
    string traceQuery;
    int editDist;
    int dpScore;
    vector<char> cigarChars;//static, fixed length field?
    vector<int> cigarLengths;
    void clear(){cigarChars.clear(); cigarLengths.clear(); editDist = dpScore = 0; traceRef = traceQuery = ""; };
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

int updateMatrix(const dp_type& type, bool print);
int dpFillOptStatic(const string& ref,const string& query, bool forward, 
        const boundaries& offset, const dp_type& type);
int dpFillOptStaticBackward(const string& ref,const string& query, bool forward, 
        const boundaries& offset, const dp_type& type);
int dpFillStatic(const string& ref,const string& query, bool forward,
        const boundaries& offset, const dp_type& type);
int dpFillStaticBackward(const string& ref,const string& query, bool forward,
        const boundaries& offset, const dp_type& type);
int findTraceBackPosStatic(bool forward, int* const i, int* const j,const dp_type& type);
int dpTraceBackStatic(int& i, int& j, const dp_type& type, dp_output& output, 
        stringstream & ss, const boundaries& offset, bool forward, const string& ref, const string& query);
int dpTraceBackOptStatic(int& i, int& j, const dp_type& type, dp_output& output, 
        stringstream & ss, const boundaries& offset, bool forward, const string& ref, const string& query);

void  print_matrices(const string& ref, const string& query, const boundaries& offset, bool gapMatrices );
void  print_seq( const string& ref,const string& query, boundaries& offset );

//int dp( const string&, const string&, boundaries&, const dp_type&, 
//        const outputType&, dp_output&, bool print = false);

//int dpBand( const string&, const string&, boundaries&, const dp_type&, 
//        const outputType&, dp_output&, int bandSize, bool print = false);

int dpBandStatic( const string&, const string&, boundaries&, 
        const dp_type&, const outputType&, dp_output&, int minbandSize, int maxbandSize, bool print = false);

int dpMyers( const string&, const string&, boundaries&, 
        const dp_type&, const outputType&, dp_output&, int minbandSize, int maxbandSize, bool print = false);

int dpBandFull( const string&, const string&, boundaries&, 
        const dp_type&, const outputType&, dp_output&, int bandLeft, int bandRight, bool print = false);

};

#endif	/* DP_H */

