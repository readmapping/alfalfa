/*
 * File:   dp.h
 * Author: mvyvermn
 *
 * Created on 14 december 2011, 12:25
 */

#ifndef DP_H
#define	DP_H


/*
 * nw.h for program nw.
 *
 */

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
    string traceRef;
    string traceQuery;
    int editDist;
    int dpScore;
    vector<char> cigarChars;
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

extern int dp( const string&, const string&, boundaries&, const dp_scores&, 
        const dp_type&, const outputType&, dp_output&, bool print = false);

#endif	/* DP_H */

