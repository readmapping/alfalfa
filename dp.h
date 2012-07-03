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

struct boundaries{
    boundaries(): refB(0), queryB(0), refE(0), queryE(0) {}
    boundaries(int refBegin, int refEnd, int queryBegin, int queryEnd): refB(refBegin), refE(refEnd), queryB(queryBegin), queryE(queryEnd){}
    int refB;
    int refE;
    int queryB;
    int queryE;
};

struct dp_scores{
    dp_scores(): match(0), mismatch(-1), extendGap(-1), openGap(0) {}
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

//These fields should maybe be positioned somewhere else for multi-threading!!!
static int ** M;
static int ** UP;
static int ** LEFT;
static int DP_DIM;
extern void initDPMatrix(int dimension, bool affine);
extern void resizeDPMatrix(int dimension, bool affine);
extern void deleteDPMatrix(bool affine);

extern int dp( const string&, const string&, boundaries&, const dp_scores&, 
        const dp_type&, const outputType&, dp_output&, bool print = false);

extern int dpBand( const string&, const string&, boundaries&, const dp_scores&, 
        const dp_type&, const outputType&, dp_output&, int bandSize, bool print = false);

extern int dpBandStatic( const string&, const string&, boundaries&, const dp_scores&, 
        const dp_type&, const outputType&, dp_output&, int bandSize, bool print = false);

#endif	/* DP_H */

