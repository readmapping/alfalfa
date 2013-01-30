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

/*-------------------------------------------
 *
 *        Dynamic programming
 *
 -------------------------------------------*/

#include <sstream>
#include <vector>
#include <assert.h>

#include "dp.h"

using namespace std;

dynProg::dynProg(int dimension, bool affine, dp_scores& scoreFunction): 
    DP_L2(0), DP_L1(0), L1(0), L2(0), bandSize(0), bandLeft(0), bandRight(0), banded(false), scores(scoreFunction) {
    scores.updateScoreMatrixDna();
    initDPMatrix(dimension, dimension, affine);
}

dynProg::~dynProg(){
    deleteDPMatrix(scores.openGap!=0);
}

void dynProg::initDPMatrix(int dimensionL2, int dimensionL1, bool affine){
    DP_L2 = dimensionL2;
    DP_L1 = dimensionL1;
    M = new int * [DP_L2];
    for(int i = 0; i < DP_L2; i++){
        M[i] = new int[DP_L1];
    }
    if(affine){
        UP = new int * [DP_L1];
        LEFT = new int * [DP_L1];
        for(int i = 0; i < DP_L2; i++){
            UP[i] = new int[DP_L1];
            LEFT[i] = new int[DP_L1];
        }
    }
}

void dynProg::resizeDPMatrix(int dimensionL2, int dimensionL1, bool affine){
    deleteDPMatrix(affine);
    initDPMatrix(dimensionL2, dimensionL1, affine);
}

void dynProg::deleteDPMatrix(bool affine){
    for(int i = 0; i < DP_L2; i++)
        delete[] M[ i ];
    delete[] M;
    if(affine){
        for(int i = 0; i < DP_L2; i++){
            delete[] UP[ i ];
            delete[] LEFT[ i ];
        }
        delete[] UP;
        delete[] LEFT;
    }
}

void dp_scores::updateScoreMatrixDna(){//presumes match/mismatch has been added correctly
    for(int i = 0; i < 127; i++)
        for(int j = 0; j < 127; j++)
            //set to value that is small enough to be used with large reads and not cause overflow with current score
            scoreMatrix[i][j] = -200000;
    for(int i = 'A'; i <= 'Z'; i++)
        for(int j = 'A'; j <= 'Z'; j++)
            scoreMatrix[i][j] = i==j ? match : mismatch;
    for(int i = 'a'; i <= 'z'; i++)
        for(int j = 'a'; j <= 'z'; j++)
            scoreMatrix[i][j] = i==j ? match : mismatch;
}

//int initializeMatrix(dp_matrices& matrices,const dp_scores& scores,const dp_type& type){
//    //TODO: change this to mallocs and initialize once per thread + realloc if nec.
//    assert(matrices.L1 >=0 && matrices.L2 >= 0);
//    matrices.M = new int * [matrices.L2+1];
//    matrices.traceBack = new char * [matrices.L2+1];
//    for(int i = 0; i <= matrices.L2; i++){
//        matrices.M[i] = new int[matrices.L1+1];
//        matrices.traceBack[i] = new char[matrices.L1+1];
//    }
//    //init matrices
//    matrices.M[0][0] = 0;
//    matrices.traceBack[0][0] = '0';//n is where traceBack stops
//    for(int i = 1; i <= matrices.L2; i++)
//        matrices.traceBack[i][0] = '0';
//    for(int i = 1; i <= matrices.L1; i++)
//        matrices.traceBack[0][i] = '0';
//    //assert(!banded || (!type.freeQueryB && !type.freeRefB));
//    if(type.freeQueryB)
//        for(int i = 1; i <= matrices.L2; i++)
//            matrices.M[i][0] = 0;
//    else
//        for(int i = 1; i <= matrices.L2; i++)
//            matrices.M[i][0] = scores.openGap + i*scores.extendGap;
//    if(type.freeRefB)
//        for(int i = 1; i <= matrices.L1; i++)
//            matrices.M[0][i] = 0;
//    else
//        for(int i = 1; i <= matrices.L1; i++)
//            matrices.M[0][i] = scores.openGap + i*scores.extendGap;
//}

int dynProg::updateMatrix(const dp_type& type){
    //TODO: change this to mallocs and initialize once per thread + realloc if nec.
    assert(L1 >=0 && L2 >= 0);
    int newDimL2 = DP_L2;
    while(L2 > newDimL2){
        newDimL2 = 2*newDimL2;
    }
    int newDimL1 = DP_L1;
    while(L1 > newDimL1){
        newDimL1 = 2*newDimL1;
    }
    if(newDimL2 > DP_L2 || newDimL1 > DP_L1)
        resizeDPMatrix(newDimL2, newDimL1, scores.openGap != 0);
    //assert(!banded || (!type.freeQueryB && !type.freeRefB));
    M[0][0] = 0;
    if(type.freeQueryB)//no need for reinitialization if same values
        for(int i = 1; i <= L2; i++)
            M[i][0] = 0;
    else
        for(int i = 1; i <= L2; i++)
            M[i][0] = scores.openGap + i*scores.extendGap;
    if(type.freeRefB)
        for(int i = 1; i <= L1; i++)
            M[0][i] = 0;
    else
        for(int i = 1; i <= L1; i++)
            M[0][i] = scores.openGap + i*scores.extendGap;
    if(scores.openGap != 0){//affine gap penalties
        UP[0][0] = LEFT[0][0] = 0;
        for(int i = 1; i <= L1; i++){
            UP[0][i] = M[0][i] + scores.openGap;
            LEFT[0][i] = M[0][i] + scores.openGap;
        }
        for(int i = 1; i <= L2; i++){
            UP[i][0] = M[i][0] + scores.openGap;
            LEFT[i][0] = M[i][0] + scores.openGap;
        }
    }
    return 0;//return is not really used, maybe for future use
}

//only called for affine gap penalties
//int initializeMatrices(dp_matrices& matrices,const dp_scores& scores,const dp_type& type){
//    initializeMatrix(matrices, scores, type);
//    matrices.gapUp = new int * [matrices.L2+1];
//    matrices.gapLeft = new int * [matrices.L2+1];
//    for(int i = 0; i <= matrices.L2; i++){
//        matrices.gapUp[i] = new int[matrices.L1+1];
//        matrices.gapLeft[i] = new int[matrices.L1+1];
//    }
//    matrices.gapLeft[0][0] = matrices.gapUp[0][0] = 0;
//    for(int i = 1; i <= matrices.L1; i++){
//        matrices.gapUp[0][i] = matrices.M[0][i] + scores.openGap;
//        matrices.gapLeft[0][i] = matrices.M[0][i] + scores.openGap;
//    }
//    for(int i = 1; i <= matrices.L2; i++){
//        matrices.gapUp[i][0] = matrices.M[i][0] + scores.openGap;
//        matrices.gapLeft[i][0] = matrices.M[i][0] + scores.openGap;
//    }
//}

//int deleteMatrices(dp_matrices& matrices,const dp_scores& scores,const outputType& oType){
//    for(int  i = 0; i <= matrices.L2; i++ ){
//        delete[] matrices.M[ i ];
//    }
//    delete[] matrices.M;
//    if(oType != DPSCORE){
//        for(int  i = 0; i <= matrices.L2; i++ ){
//            delete[] matrices.traceBack[ i ];
//        }
//        delete[] matrices.traceBack;
//    }
//    if(scores.openGap!=0){
//        for(int  i = 0; i <= matrices.L2; i++ ){
//            delete[] matrices.gapUp[ i ];
//            delete[] matrices.gapLeft[ i ];
//        }
//        delete[] matrices.gapUp;
//        delete[] matrices.gapLeft;
//    }
//    else{
//        matrices.gapUp = NULL;
//        matrices.gapLeft = NULL;
//    }
//    return 0;
//}

//int deleteMatrix(dp_matrices& matrices){
//    for(int  i = 0; i <= matrices.L2; i++ ){
//        delete[] matrices.M[ i ];
//        delete[] matrices.traceBack[ i ];
//    }
//    delete[] matrices.M;
//    delete[] matrices.traceBack;
//    matrices.gapUp = NULL;
//    matrices.gapLeft = NULL;
//    return 0;
//}

//not yet adjusted for different band sizes
//int findTraceBackPos(const dp_matrices& matrices, int* const i, int* const j,const dp_type& type){
//    //find beginPosition for traceBack and find dpScore
//    int maximum = matrices.M[*i][*j];
//    if((type.freeQueryE && type.freeRefE) || type.local){//local
//        if(!matrices.banded){
//            for(int idx = matrices.L2; idx >= 0; idx--){
//                for(int idx2 = matrices.L1; idx2 >= 0; idx2--){
//                    if(matrices.M[idx][idx2] > maximum){
//                        maximum = matrices.M[idx][idx2];
//                        *i = idx;
//                        *j = idx2;
//                    }
//                }
//            }
//        }
//        else{
//            for(int idx = matrices.L2; idx >= 0; idx--){
//                for(int idx2 = min(idx+matrices.bandSize,matrices.L1); idx2 >= max(0,idx-matrices.bandSize); idx2--){
//                    if(matrices.M[idx][idx2] > maximum){
//                        maximum = matrices.M[idx][idx2];
//                        *i = idx;
//                        *j = idx2;
//                    }
//                }
//            }
//        }
//    }
//    else if(type.freeQueryE){//banded update
//        for(int idx = min(matrices.L2, matrices.L1 + matrices.bandSize); 
//                idx >= min(matrices.L2, max(0,matrices.L1 - matrices.bandSize)); 
//                idx--){
//            if(matrices.M[idx][*j] > maximum){
//                maximum = matrices.M[idx][*j];
//                *i = idx;
//            }
//        }
//    }
//    else if(type.freeRefE){
//        for(int idx = min(matrices.L1, matrices.L2 + matrices.bandSize); 
//                idx >= min(matrices.L1, max(0,matrices.L2 - matrices.bandSize)); 
//                idx--){
//            if(matrices.M[*i][idx] > maximum){
//                maximum = matrices.M[*i][idx];
//                *j = idx;
//            }
//        }
//    }
//    return 0;
//}

//--> possible bug ? do not take into account forward or reverse
//it was ok for non-local dp or when gaps are not allowed in both begin and end, which was the only case untill now
int dynProg::findTraceBackPosStatic(bool forward, int* const i, int* const j,const dp_type& type){
    //find beginPosition for traceBack and find dpScore
    int maximum = M[*i][*j];
    int bandCorrection = forward ? 0 : L1 - L2;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    if((type.freeQueryE && type.freeRefE) || type.local){//local
        if(!banded){
            for(int idx = L2; idx >= 0; idx--){
                for(int idx2 = L1; idx2 >= 0; idx2--){
                    if(M[idx][idx2] > maximum){
                        maximum = M[idx][idx2];
                        *i = idx;
                        *j = idx2;
                    }
                }
            }
        }
        else{
            for(int idx = L2; idx >= 0; idx--){
                for(int idx2 = min(bandCorrection+idx+bandR,L1); idx2 >= max(0,bandCorrection+idx-bandL); idx2--){
                    if(M[idx][idx2] > maximum){
                        maximum = M[idx][idx2];
                        *i = idx;
                        *j = idx2;
                    }
                }
            }
        }
    }
    else if(type.freeQueryE){//banded update
        for(int idx = min(L2, bandCorrection + L1 + bandR); 
                idx >= min(L2, max(0, bandCorrection + L1 - bandL)); 
                idx--){
            if(M[idx][*j] > maximum){
                maximum = M[idx][*j];
                *i = idx;
            }
        }
    }
    else if(type.freeRefE){
        for(int idx = min(L1, bandCorrection + L2 + bandR); 
                idx >= min(L1, max(0, bandCorrection + L2 - bandL)); 
                idx--){
            if(M[*i][idx] > maximum){
                maximum = M[*i][idx];
                *j = idx;
            }
        }
    }
    return 0;
}

inline int  maximum( int f1, int f2, int f3, char * ptr )
{
        int  max = 0 ;

        if( f2 >= f1 && f2 >= f3 )
        {
                max = f2 ;
                *ptr = '\\' ;
        }
        else if( f1 > f3 )
        {
                max = f1 ;
                *ptr = '|' ;
        }
        else
        {
                max = f3 ;
                *ptr = '-' ;
        }

        return  max ;
}

//not yet adjusted for different band sizes
//int dpFill(dp_matrices& matrices,const string& ref,const string& query, bool forward,
//        const boundaries& offset, const dp_scores& scores, 
//        const dp_type& type){
//    int        d = 0;
//    int        fU, fD, fL ;
//    char       ptr;
//
//    for(int i = 1; i <= matrices.L2; i++ ){
//        int idx2 = forward ? i-matrices.bandSize : matrices.L1-matrices.L2 + i - matrices.bandSize;
//        int colRB = forward ? i+matrices.bandSize : matrices.L1-matrices.L2 + i + matrices.bandSize;
//        if(idx2>0){//left value undefined
//            ptr = '\\';
//            d = scores.scoreMatrix[ref[offset.refB+idx2-1]][query[offset.queryB+i-1]];
//            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ] + d;
//            matrices.gapUp[i][idx2] = max(matrices.M[ i-1 ][ idx2 ] + scores.extendGap + scores.openGap, matrices.gapUp[i-1][idx2]+scores.extendGap);
//            if(matrices.M[ i ][ idx2 ] < matrices.gapUp[i][idx2]){
//                matrices.M[ i ][ idx2 ] = matrices.gapUp[i][idx2];
//                ptr = '|';
//            }
//            //fill in generic value for gapLeft, such that it is defined, but not chosen in next step
//            matrices.gapLeft[i][idx2] = matrices.M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
//            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
//                matrices.M[i][idx2] = 0;
//                ptr = '0';
//            }
//            idx2++;
//            matrices.traceBack[i][idx2] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
//        }
//        for(int j = max( 1, idx2); j <= min(matrices.L1, colRB-1); j++ ){
//            //here insert substitution matrix!
//            d = scores.scoreMatrix[ref[offset.refB+j-1]][query[offset.queryB+i-1]];
//            fU = max(matrices.M[ i-1 ][ j ] + scores.extendGap + scores.openGap, matrices.gapUp[i-1][j]+scores.extendGap);
//            fD = matrices.M[ i-1 ][ j-1 ] + d;
//            fL = max(matrices.M[ i ][ j-1 ] + scores.extendGap + scores.openGap, matrices.gapLeft[i][j-1]+scores.extendGap);
//            matrices.gapUp[i][j] = fU;
//            matrices.gapLeft[i][j] = fL;
//            matrices.M[ i ][ j ] = maximum( fU, fD, fL, &ptr );
//            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][j]<0){
//                matrices.M[i][j] = 0;
//                matrices.traceBack[i][j] = '0';
//            }
//            else
//                matrices.traceBack[i][j] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
//        }
//        if(colRB <= matrices.L1){
//            idx2 = colRB;
//            ptr = '\\';
//            d = scores.scoreMatrix[ref[offset.refB+idx2-1]][query[offset.queryB+i-1]];
//            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ] + d;
//            matrices.gapLeft[i][idx2] = max(matrices.M[ i ][ idx2-1 ] + scores.extendGap + scores.openGap, matrices.gapLeft[i][idx2-1]+scores.extendGap);
//            if(matrices.M[ i ][ idx2 ] < matrices.gapLeft[i][idx2]){
//                matrices.M[ i ][ idx2 ] = matrices.gapLeft[i][idx2];
//                ptr = '-';
//            }
//            //fill in generic value for gapUp, such that it is defined, but not chosen in next step
//            matrices.gapUp[i][idx2] = matrices.M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
//            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
//                matrices.M[i][idx2] = 0;
//                ptr = '0';
//            }
//            matrices.traceBack[i][idx2] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
//        }
//    }
//    return 0;
//}

int dynProg::dpFillStatic(const string& ref,const string& query, bool forward,
        const boundaries& offset, const dp_type& type){
    int        d = 0;
    int        local = type.local || (type.freeQueryB && type.freeRefB) ? 0 : 1;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    for(long i = 1; i <= L2; i++ ){
        long idx2 = forward ? i-bandL : L1 - L2 + i - bandL;
        long colRB = forward ? i+bandR : L1 - L2 + i + bandR;
        if(idx2>0){//left value undefined
            d = scores.scoreMatrix[size_t (ref[offset.refB+idx2-1])][size_t (query[offset.queryB+i-1])];
            M[ i ][ idx2 ] = M[ i-1 ][ idx2-1 ] + d;
            UP[i][idx2] = max(M[ i-1 ][ idx2 ] + scores.extendGap + scores.openGap, UP[i-1][idx2]+scores.extendGap);
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ] , UP[i][idx2]);
            //fill in generic value for gapLeft, such that it is defined, but not chosen in next step
            LEFT[i][idx2] = M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ], 0 + local*M[ i ][ idx2 ]);
            idx2++;
        }
        for(long j = max( 1L, idx2); j <= min((long)L1, colRB-1); j++ ){
            //here insert substitution matrix!
            d = scores.scoreMatrix[size_t(ref[offset.refB+j-1])][size_t(query[offset.queryB+i-1])];
            UP[i][j] = max(M[ i-1 ][ j ] + scores.extendGap + scores.openGap, UP[i-1][j]+scores.extendGap);
            M[i][j] = M[ i-1 ][ j-1 ] + d;
            LEFT[i][j] = max(M[ i ][ j-1 ] + scores.extendGap + scores.openGap, LEFT[i][j-1]+scores.extendGap);
            M[i][j] = max(M[i][j], LEFT[i][j]);
            M[i][j] = max(M[i][j], UP[i][j]);
            M[ i ][ j ] = max(M[ i ][ j ], 0 + local*M[ i ][ j ]);
        }
        if(colRB <= L1){
            idx2 = colRB;
            d = scores.scoreMatrix[size_t(ref[offset.refB+idx2-1])][size_t(query[offset.queryB+i-1])];
            M[ i ][ idx2 ] = M[ i-1 ][ idx2-1 ] + d;
            LEFT[i][idx2] = max(M[ i ][ idx2-1 ] + scores.extendGap + scores.openGap, LEFT[i][idx2-1]+scores.extendGap);
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ], LEFT[i][idx2]);
            //fill in generic value for gapUp, such that it is defined, but not chosen in next step
            UP[i][idx2] = M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ], 0 + local*M[ i ][ idx2 ]);
        }
    }
    return 0;
}

//not yet adjusted for different band sizes
//int dpFillOpt(dp_matrices& matrices,const string& ref,const string& query, bool forward, 
//        const boundaries& offset, const dp_scores& scores, const dp_type& type){
//    //TODO: when changed to mallocs: use row-order traversal
//    int        d = 0;
//    char       ptr;
//    
//    for(int i = 1; i <= matrices.L2; i++ ){
//        int idx2 = forward ? i-matrices.bandSize : matrices.L1-matrices.L2 + i - matrices.bandSize;
//        int colRB = forward ? i+matrices.bandSize : matrices.L1-matrices.L2 + i + matrices.bandSize;
//        if(idx2>0){
//            ptr = '\\';
//            d = scores.scoreMatrix[ref[offset.refB+idx2-1]][query[offset.queryB+i-1]];
//            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ] + d;
//            if(matrices.M[ i ][ idx2 ] < matrices.M[ i-1 ][ idx2 ] + scores.extendGap){
//                matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2 ] + scores.extendGap;
//                ptr = '|';
//            }
//            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
//                matrices.M[i][idx2] = 0;
//                ptr = '0';
//            }
//            matrices.traceBack[i][idx2] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
//            idx2++;
//        }
//        for(int j = max( 1, idx2); j <= min(matrices.L1, colRB-1); j++ ){
//            d = scores.scoreMatrix[ref[offset.refB+j-1]][query[offset.queryB+i-1]];
//            matrices.M[ i ][ j ] = maximum( matrices.M[ i-1 ][ j ] + scores.extendGap,
//                                matrices.M[ i-1 ][ j-1 ] + d,
//                                matrices.M[ i ][ j-1 ] + scores.extendGap,
//                                &ptr );
//            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][j]<0){//TODO: put this if up front for sw programming
//                matrices.M[i][j] = 0;
//                matrices.traceBack[i][j] = '0';
//            }
//            else
//                matrices.traceBack[i][j] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
//            //TODO: remove if by using table on previous match/mismatch table and storing the pointer in the max function
//        }
//        if(colRB <= matrices.L1){
//            idx2 = colRB;
//            ptr = '\\';
//            d = scores.scoreMatrix[ref[offset.refB+idx2-1]][query[offset.queryB+i-1]];
//            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ] + d;
//            if(matrices.M[ i ][ idx2 ] < matrices.M[ i ][ idx2-1 ] + scores.extendGap){
//                matrices.M[ i ][ idx2 ] = matrices.M[ i ][ idx2-1 ] + scores.extendGap;
//                ptr = '-';
//            }
//            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
//                matrices.M[i][idx2] = 0;
//                ptr = '0';
//            }
//            matrices.traceBack[i][idx2] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
//        }
//    }
//    return 0;
//}

int dynProg::dpFillOptStatic(const string& ref,const string& query, bool forward, 
        const boundaries& offset, const dp_type& type){
    //TODO: when changed to mallocs: use row-order traversal
    int        d = 0;
    int        local = type.local || (type.freeQueryB && type.freeRefB) ? 0 : 1;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    for(long i = 1; i <= L2; i++ ){
        long idx2 = forward ? i-bandL : L1 - L2 + i - bandL;
        long colRB = forward ? i+bandR : L1 - L2 + i + bandR;
        if(idx2>0){
            d = scores.scoreMatrix[size_t(ref[offset.refB+idx2-1])][size_t(query[offset.queryB+i-1])];
            M[ i ][ idx2 ] = max(M[ i-1 ][ idx2-1 ] + d, M[ i-1 ][ idx2 ] + scores.extendGap);
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ], 0 + local*M[ i ][ idx2 ]);
            idx2++;
        }
        for(long j = max( 1L, idx2); j <= min((long)L1, colRB-1); j++ ){
            d = scores.scoreMatrix[size_t(ref[offset.refB+j-1])][size_t(query[offset.queryB+i-1])];
            M[ i ][ j ] = max(M[ i-1 ][ j-1 ] + d, M[ i ][ j-1 ] + scores.extendGap);
            M[ i ][ j ] = max(M[ i ][ j ], M[ i-1 ][ j ] + scores.extendGap);
            M[ i ][ j ] = max(M[ i ][ j ], 0 + local*M[ i ][ j ]);
        }
        if(colRB <= L1){
            idx2 = colRB;
            d = scores.scoreMatrix[size_t(ref[offset.refB+idx2-1])][size_t(query[offset.queryB+i-1])];
            M[ i ][ idx2 ] = max(M[ i-1 ][ idx2-1 ] + d, M[ i ][ idx2-1 ] + scores.extendGap);
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ], 0 + local*M[ i ][ idx2 ]);
        }
    }
    return 0;
}

//not yet adjusted for different band sizes
//int dpTraceBack(const dp_matrices& matrices, int& i, int& j, dp_output& output, std::stringstream & ss){
//    //TODO: after mallocs: change this traversal to row-order
//    //change stringstream to char array of length dimensions array
//    while( matrices.traceBack[ i ][ j ] != '0' ){
//        switch( matrices.traceBack[ i ][ j ] ){
//            case '|' :      ss << 'I';
//                            output.editDist++;
//                            i-- ;
//                            break ;
//
//            case '\\':      ss << '=';
//                            i-- ;  j-- ;
//                            break ;
//
//            case ':':       ss << 'X';
//                            output.editDist++;
//                            i-- ;  j-- ;
//                            break ;
//
//            case '-' :      ss << 'D';
//                            output.editDist++;
//                            j-- ;
//                            break;
//            default :       cerr << "Error in dynamic programming" << endl;
//                            break;
//        }
//    }
//    return 0;
//}

int dynProg::dpTraceBackStatic(int& i, int& j, const dp_type& type, dp_output& output, 
        std::stringstream & ss, const boundaries& offset, bool forward, const string& ref, const string& query){
    //TODO: after mallocs: change this traversal to row-order
    //change stringstream to char array of length dimensions array
    int qBound = 0; int rBound = 0;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    if(type.freeQueryB)
        qBound = L2;
    if(type.freeRefB)
        rBound = L1;
    bool local = type.local || (type.freeQueryB && type.freeRefB);
    int bandCorrection = forward ? 0 : L1 - L2;
    while( (i > qBound || j > rBound) && (!local || M[ i ][ j ] > 0)){
        if( j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.match && ref[offset.refB+j-1] == query[offset.queryB+i-1]){
            ss << '=';
            i-- ;  j-- ;
        }
        else if(j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.mismatch){
            ss << 'X';
            output.editDist++;
            i-- ;  j-- ;
        }
        else if(j > 0 && j > i+bandCorrection - bandL && M[ i ][ j ] == LEFT[ i ][ j ]){
            ss << 'D';
            output.editDist++;
            j-- ;
        }
        else if(i > 0 && j < i+bandCorrection + bandR && M[ i ][ j ] == UP[ i ][ j ]){
            ss << 'I';
            output.editDist++;
            i-- ;
        }
        else{
            cerr << "Error in dynamic programming" << endl;
            cerr << "No valid trace position found" << endl;
            cerr << "Position " << i << ", " << j << endl;
            cerr << "Scores: Current " << M[ i ][ j ] << endl;
            if(j > 0 && i > 0){
                cerr << "Scores: Diag "    << M[ i-1 ][ j-1 ] << endl;
                cerr << "Characters "    << ref[offset.refB+j-1] << "  " << query[offset.queryB+i-1] << endl; 
            }
            cerr << "BandSize " << max(max(bandSize, bandL), bandRight) << endl;
            if(j > 0)
                cerr << "Scores: Left "    << M[ i ][ j-1 ] << endl;
            if(i > 0)
                cerr << "Scores: Up "    << M[ i-1 ][ j ] << endl;
            exit(1);
        }
    }
    return 0;
}

int dynProg::dpTraceBackOptStatic(int& i, int& j, const dp_type& type, dp_output& output, 
        std::stringstream & ss, const boundaries& offset, bool forward, const string& ref, const string& query){
    //TODO: after mallocs: change this traversal to row-order
    //change stringstream to char array of length dimensions array
    int qBound = 0; int rBound = 0;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    if(type.freeQueryB)
        qBound = L2;
    if(type.freeRefB)
        rBound = L1;
    bool local = type.local || (type.freeQueryB && type.freeRefB);
    int bandCorrection = forward ? 0 : L1 - L2;
    while( (i > qBound || j > rBound) && (!local || M[ i ][ j ] > 0)){
        if( j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.match && ref[offset.refB+j-1] == query[offset.queryB+i-1]){
            ss << '=';
            i-- ;  j-- ;
        }
        else if(j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.mismatch){
            ss << 'X';
            output.editDist++;
            i-- ;  j-- ;
        }
        else if(j > 0 && j > i+bandCorrection - bandL && M[ i ][ j ] == M[ i ][ j-1 ] + scores.extendGap){
            ss << 'D';
            output.editDist++;
            j-- ;
        }
        else if(i > 0 && j < i+bandCorrection + bandR && M[ i ][ j ] == M[ i-1 ][ j ] + scores.extendGap){
            ss << 'I';
            output.editDist++;
            i-- ;
        }
        else{
            cerr << "Error in dynamic programming" << endl;
            cerr << "No valid trace position found" << endl;
            cerr << "Position " << i << ", " << j << endl;
            cerr << "Scores: Current " << M[ i ][ j ] << endl;
            if(j > 0 && i > 0){
                cerr << "Scores: Diag "    << M[ i-1 ][ j-1 ] << endl;
                cerr << "Characters "    << ref[offset.refB+j-1] << "  " << query[offset.queryB+i-1] << endl; 
            }
            cerr << "BandSize " << max(max(bandSize, bandL), bandR) << endl;
            if(j > 0)
                cerr << "Scores: Left "    << M[ i ][ j-1 ] << endl;
            if(i > 0)
                cerr << "Scores: Up "    << M[ i-1 ][ j ] << endl;
            exit(1);
        }
    }
    return 0;
}

void  dynProg::print_matrices(const string& ref, 
        const string& query, const boundaries& offset, bool gapMatrices )
{
        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << ref[ offset.refB+j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << query[ offset.queryB+i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << M[ i ][ j ] << " ";
                }
                cout << endl;
        }
        cout << endl;
        if(gapMatrices){
            cout << "        ";
            for( int j = 0; j < L1; j++ )
            {
                    cout << ref[ offset.refB+j ] << "   ";
            }
            cout << "\n  ";

            for( int i = 0; i <= L2; i++ )
            {
                    if( i > 0 )
                    {
                            cout << query[ offset.queryB+i-1 ] << " ";
                    }
                    for( int j = 0; j <= L1; j++ )
                    {
                            cout.width( 3 );
                            cout << LEFT[ i ][ j ] << " ";
                    }
                    cout << endl;
            }
            cout << endl;
            cout << "        ";
            for( int j = 0; j < L1; j++ )
            {
                    cout << ref[ offset.refB+j ] << "   ";
            }
            cout << "\n  ";

            for( int i = 0; i <= L2; i++ )
            {
                    if( i > 0 )
                    {
                            cout << query[ offset.queryB+i-1 ] << " ";
                    }
                    for( int j = 0; j <= L1; j++ )
                    {
                            cout.width( 3 );
                            cout << UP[ i ][ j ] << " ";
                    }
                    cout << endl;
            }
            cout << endl;
        }
        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << ref[ offset.refB+j ] << "   ";
        }
        cout << "\n  ";

//        for( int i = 0; i <= L2; i++ )
//        {
//                if( i > 0 )
//                {
//                        cout << query[ offset.queryB+i-1 ] << " ";
//                }
//                for( int j = 0; j <= L1; j++ )
//                {
//                        cout.width( 3 );
//                        cout << traceBack[ i ][ j ] << " ";
//                }
//                cout << endl;
//        }
        cout << endl;
}

void  dynProg::print_seq( const string& ref,const string& query, boundaries& offset )
{
        cout << "" << endl << ref.substr(offset.refB,offset.refE-offset.refB+1) << endl;
        cout << query.substr(offset.queryB,offset.queryE-offset.queryB+1) << endl << "" << endl;
}

////(extra parameters offset1 and offset2, length1 and length2 in struct boundaries)
////extra: substitution matrix!
//int dp( const string&     ref,
//        const string&     query,
//        boundaries& offset,
//        const dp_scores&  scores,
//        const dp_type&    type,
//        const outputType& oType,
//        dp_output&  output,
//        bool print)
//{
//    if(print) print_seq(ref, query, offset);
//    dp_matrices matrices;
//    matrices.L1 = offset.refE-offset.refB+1;
//    matrices.L2 = offset.queryE-offset.queryB+1;
//    matrices.bandSize = max(matrices.L1, matrices.L2) + 1;
//    matrices.banded = false;
//    assert(matrices.L1 >=0 && matrices.L2 >= 0);
//    if(scores.openGap != 0){
//        //build matrices
//        initializeMatrices(matrices, scores, type);
//        //Check for affine gap: if else
//        //dp
//        dpFill(matrices, ref, query, true, offset, scores, type);
//        if(print) print_matrices(matrices, ref, query, offset, scores.openGap!=0);
//        //traceback and output
//        int        i = matrices.L2, j = matrices.L1;
//        //find beginPosition for traceBack and find dpScore
//        findTraceBackPos(matrices,&i,&j,type);
//        output.dpScore = matrices.M[i][j];
//        //Next is for output only
//        if(oType>DPSCORE){
//            //alignment output
//            //Alter boundaries to return the new positions when free gaps were introduced.
//            if((type.local || type.freeQueryE) && i < matrices.L2)
//                offset.queryE = offset.queryB+i-1;
//            if((type.local || type.freeRefE) && j < matrices.L1)
//                offset.refE = offset.refB+j-1;
//            //traceBack
//            std::stringstream ss;
//            dpTraceBack(matrices, i, j, output, ss);
//            //Edit Dist ouput
//            if( i > 0 && !type.local && !type.freeQueryB)
//                output.editDist += i;
//            if(j > 0 && !type.local && !type.freeRefB)
//                output.editDist += j;
//            //Errorstring output
//            if(oType == ERRORSTRING || oType == ALL){
//                string errorString = ss.str();
//                if(i > 0 && !type.local && !type.freeQueryB)
//                    errorString.append(i,'I');
//                else if(j > 0 && !type.local && !type.freeRefB)
//                    errorString.append(j, 'D');
//                reverse( errorString.begin(), errorString.end() );
//                //create output: errorString
//                int iter = 0;
//                int temp;
//                while(iter < errorString.length()){
//                    temp = iter;
//                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
//                        iter++;
//                    iter++;
//                    output.cigarChars.push_back(errorString[temp]);
//                    output.cigarLengths.push_back(iter-temp);
//                }
//            }
//        }
//        //change boundaries for local and freeGap types to return changed free boundaries.
//        if((type.local || type.freeQueryB) && i > 0)
//            offset.queryB = offset.queryB+i;
//        if((type.local || type.freeRefB) && j > 0)
//            offset.refB = offset.refB+j;
//        deleteMatrices(matrices, scores, oType);
//    }
//    else{
//        //build matrices
//        initializeMatrix(matrices, scores, type);
//        //dp
//        dpFillOpt(matrices, ref, query, true, offset, scores, type);
//        //traceback and output
//        int        i = matrices.L2, j = matrices.L1;
//        //find beginPosition for traceBack and find dpScore
//        findTraceBackPos(matrices,&i,&j,type);
//        output.dpScore = matrices.M[i][j];
//        //Next is for output only
//        if(oType>DPSCORE){
//            //Alter boundaries to return the new positions when free gaps were introduced.
//            if((type.local || type.freeQueryE) && i < matrices.L2)
//                offset.queryE = offset.queryB+i-1;
//            if((type.local || type.freeRefE) && j < matrices.L1)
//                offset.refE = offset.refB+j-1;
//            //traceBack
//            std::stringstream ss;
//            dpTraceBack(matrices, i, j, output, ss);
//            //Edit Dist ouput
//            if( i > 0 && !type.local && !type.freeQueryB)
//                output.editDist += i;
//            if(j > 0 && !type.local && !type.freeRefB)
//                output.editDist += j;
//            //Errorstring output
//            if(oType == ERRORSTRING || oType == ALL){
//                string errorString = ss.str();
//                if(i > 0 && !type.local && !type.freeQueryB)
//                    errorString.append(i,'I');
//                else if(j > 0 && !type.local && !type.freeRefB)
//                    errorString.append(j, 'D');
//                reverse( errorString.begin(), errorString.end() );//TODO: leave this out, but use different iterator
//                //create output: errorString
//                int iter = 0;
//                int temp;
//                while(iter < errorString.length()){
//                    temp = iter;
//                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
//                        iter++;
//                    iter++;
//                    output.cigarChars.push_back(errorString[temp]);
//                    output.cigarLengths.push_back(iter-temp);
//                }
//            }
//        }
//        //change boundaries for local and freeGap types to return changed free boundaries.
//        if((type.local || type.freeQueryB) && i > 0)
//            offset.queryB = offset.queryB+i;
//        if((type.local || type.freeRefB) && j > 0)
//            offset.refB = offset.refB+j;
//        deleteMatrix(matrices);
//    }
//
//    return 0;
//}

//(extra parameters offset1 and offset2, length1 and length2 in struct boundaries)
//extra: substitution matrix!
//Banded alignment only for band sizes > 0 and non-local alignment!!!
//int dpBand( const string&     ref,
//        const string&     query,
//        boundaries& offset,
//        const dp_scores&  scores,
//        const dp_type&    type,
//        const outputType& oType,
//        dp_output&  output,
//        int bandSize,
//        bool print)
//{
//    dp_matrices matrices;
//    matrices.L1 = offset.refE-offset.refB+1;
//    matrices.L2 = offset.queryE-offset.queryB+1;
//    matrices.bandSize = bandSize;
//    matrices.banded = true;
//    assert(matrices.L1 >=0 && matrices.L2 >= 0 && bandSize >0);
//    assert(!dp_type.local);
//    assert((!dp_type.freeQueryB && !dp_type.freeRefB) || (!dp_type.freeQueryE && !dp_type.freeRefE));
//    //IF no free start at begin:
//    //possible end of alignment should be accessible for traceback: 
//    //range of allowed traceback startpositions should have a non-empty intersection
//    //with the band, otherwise: give bad alignment and return
//    //ELSE IF free start at begin in cols or rows: (resulting in no free ending):
//    //check range of starting point (0,0) lies in range of band if rows or columns were not free
//    if((type.freeQueryB && !type.freeRefB && matrices.L2 + bandSize < matrices.L1) ||
//            (type.freeRefB && !type.freeQueryB && matrices.L2 - bandSize > matrices.L1) ||
//            (!(type.freeQueryB || type.freeRefB) && !type.freeRefE && matrices.L2 + bandSize < matrices.L1) ||
//            (!(type.freeQueryB || type.freeRefB) && !type.freeQueryE && matrices.L2 - bandSize > matrices.L1)){
//        output.dpScore = scores.mismatch*query.length();
//        output.editDist = query.length();
//        return 1;
//    }
//    
//    if(scores.openGap != 0){
//        //build matrices
//        initializeMatrices(matrices, scores, type);
//        //Check for affine gap: if else
//        //dp
//        dpFill(matrices, ref, query, !type.freeQueryB && !type.freeRefB, offset, scores, type);
//        if(print) print_matrices(matrices, ref, query, offset, scores.openGap!=0);
//        //traceback and output
//        int        i = min(matrices.L2,matrices.L1+matrices.bandSize), j = min(matrices.L1,matrices.L2+matrices.bandSize);
//        //find beginPosition for traceBack and find dpScore
//        findTraceBackPos(matrices,&i,&j,type);
//        output.dpScore = matrices.M[i][j];
//        //Next is for output only
//        if(oType>DPSCORE){
//            //Alter boundaries to return the new positions when free gaps were introduced.
//            if(type.freeQueryE && i < matrices.L2)
//                offset.queryE = offset.queryB+i-1;
//            if(type.freeRefE && j < matrices.L1)
//                offset.refE = offset.refB+j-1;
//            //traceBack
//            std::stringstream ss;
//            dpTraceBack(matrices, i, j, output, ss);
//            //Edit Dist ouput
//            if( i > 0 && !type.freeQueryB)
//                output.editDist += i;
//            if(j > 0 && !type.freeRefB)
//                output.editDist += j;
//            //Errorstring output
//            if(oType == ERRORSTRING || oType == ALL){
//                string errorString = ss.str();
//                if(i > 0 && !type.freeQueryB)
//                    errorString.append(i,'I');
//                else if(j > 0 && !type.freeRefB)
//                    errorString.append(j, 'D');
//                reverse( errorString.begin(), errorString.end() );
//                //create output: errorString
//                int iter = 0;
//                int temp;
//                while(iter < errorString.length()){
//                    temp = iter;
//                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
//                        iter++;
//                    iter++;
//                    output.cigarChars.push_back(errorString[temp]);
//                    output.cigarLengths.push_back(iter-temp);
//                }
//            }
//        }
//        //change boundaries for local and freeGap types to return changed free boundaries.
//        if(type.freeQueryB && i > 0)
//            offset.queryB = offset.queryB+i;
//        if(type.freeRefB && j > 0)
//            offset.refB = offset.refB+j;
//        deleteMatrices(matrices, scores, oType);
//    }
//    else{
//        //build matrices
//        initializeMatrix(matrices, scores, type);
//        //dp
//        dpFillOpt(matrices, ref, query, !type.freeQueryB && !type.freeRefB, offset, scores, type);
//        //traceback and output
//        int        i = min(matrices.L2,matrices.L1+matrices.bandSize), j = min(matrices.L1,matrices.L2+matrices.bandSize);
//        //find beginPosition for traceBack and find dpScore
//        findTraceBackPos(matrices,&i,&j,type);
//        output.dpScore = matrices.M[i][j];
//        //Next is for output only
//        if(oType>DPSCORE){
//            //Alter boundaries to return the new positions when free gaps were introduced.
//            if((type.local || type.freeQueryE) && i < matrices.L2)
//                offset.queryE = offset.queryB+i-1;
//            if((type.local || type.freeRefE) && j < matrices.L1)
//                offset.refE = offset.refB+j-1;
//            //traceBack
//            std::stringstream ss;
//            dpTraceBack(matrices, i, j, output, ss);
//            //Edit Dist ouput
//            if( i > 0 && !type.local && !type.freeQueryB)
//                output.editDist += i;
//            if(j > 0 && !type.local && !type.freeRefB)
//                output.editDist += j;
//            //Errorstring output
//            if(oType == ERRORSTRING || oType == ALL){
//                string errorString = ss.str();
//                if(i > 0 && !type.local && !type.freeQueryB)
//                    errorString.append(i,'I');
//                else if(j > 0 && !type.local && !type.freeRefB)
//                    errorString.append(j, 'D');
//                reverse( errorString.begin(), errorString.end() );//TODO: leave this out, but use different iterator
//                //create output: errorString
//                int iter = 0;
//                int temp;
//                while(iter < errorString.length()){
//                    temp = iter;
//                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
//                        iter++;
//                    iter++;
//                    output.cigarChars.push_back(errorString[temp]);
//                    output.cigarLengths.push_back(iter-temp);
//                }
//            }
//        }
//        //change boundaries for local and freeGap types to return changed free boundaries.
//        if((type.local || type.freeQueryB) && i > 0)
//            offset.queryB = offset.queryB+i;
//        if((type.local || type.freeRefB) && j > 0)
//            offset.refB = offset.refB+j;
//        deleteMatrix(matrices);
//    }
//
//    return 0;
//}
//(extra parameters offset1 and offset2, length1 and length2 in struct boundaries)
//extra: substitution matrix!
//Banded alignment only for band sizes > 0 and non-local alignment!!!


int dynProg::dpBandStatic( const string&     ref,
        const string&     query,
        boundaries& offset,
        const dp_type&    type,
        const outputType& oType,
        dp_output&  output,
        int bandS,
        bool print)
{
    L1 = offset.refE-offset.refB+1;
    L2 = offset.queryE-offset.queryB+1;
    bandSize = bandS;
    banded = true;
    assert(L1 >=0 && L2 >= 0 && bandSize >0);
    assert(!type.local);
    assert((!type.freeQueryB && !type.freeRefB) || (!type.freeQueryE && !type.freeRefE));
    //IF no free start at begin:
    //possible end of alignment should be accessible for traceback: 
    //range of allowed traceback startpositions should have a non-empty intersection
    //with the band, otherwise: give bad alignment and return
    //ELSE IF free start at begin in cols or rows: (resulting in no free ending):
    //check range of starting point (0,0) lies in range of band if rows or columns were not free
    if((type.freeQueryB && !type.freeRefB && L2 + bandSize < L1) ||
            (type.freeRefB && !type.freeQueryB && L2 - bandSize > L1) ||
            (!(type.freeQueryB || type.freeRefB) && !type.freeRefE && L2 + bandSize < L1) ||
            (!(type.freeQueryB || type.freeRefB) && !type.freeQueryE && L2 - bandSize > L1)){
        output.dpScore = scores.mismatch*query.length();
        output.editDist = query.length();
        return 1;
    }
    //build matrices
    updateMatrix(type);
    bool affine = scores.openGap != 0;
    //dp
    if(affine)
        dpFillStatic(ref, query, !type.freeQueryB && !type.freeRefB, offset, type);
    else
        dpFillOptStatic(ref, query, !type.freeQueryB && !type.freeRefB, offset, type);
    //traceback and output
    int        i = min(L2,L1+bandSize), j = min(L1,L2+bandSize);
    //find beginPosition for traceBack and find dpScore
    findTraceBackPosStatic((!type.freeQueryB && !type.freeRefB), &i,&j,type);
    output.dpScore = M[i][j];
    //Next is for output only
    if(oType>DPSCORE){
        //Alter boundaries to return the new positions when free gaps were introduced.
        if(type.freeQueryE && i < L2)
            offset.queryE = offset.queryB+i-1;
        if(type.freeRefE && j < L1)
            offset.refE = offset.refB+j-1;
        //traceBack
        std::stringstream ss;
        if(affine)
            dpTraceBackStatic(i, j, type, output, ss, offset, !type.freeQueryB && !type.freeRefB, ref, query);
        else
            dpTraceBackOptStatic(i, j, type, output, ss, offset, !type.freeQueryB && !type.freeRefB, ref, query);
        //Edit Dist ouput
        if( i > 0 && !type.freeQueryB)
            output.editDist += i;
        if(j > 0 && !type.freeRefB)
            output.editDist += j;
        //Errorstring output
        if(oType == ERRORSTRING || oType == ALL){
            string errorString = ss.str();
            if(i > 0 && !type.freeQueryB)
                errorString.append(i,'I');
            else if(j > 0 && !type.freeRefB)
                errorString.append(j, 'D');
            reverse( errorString.begin(), errorString.end() );
            //create output: errorString
            int iter = 0;
            int temp;
            int errorStringLength = errorString.length();
            while(iter < errorStringLength){
                temp = iter;
                while(iter+1 < errorStringLength && errorString[iter+1]==errorString[temp])
                    iter++;
                iter++;
                output.cigarChars.push_back(errorString[temp]);
                output.cigarLengths.push_back(iter-temp);
            }
        }
    }
    //change boundaries for local and freeGap types to return changed free boundaries.
    if(type.freeQueryB && i > 0)
        offset.queryB = offset.queryB+i;
    if(type.freeRefB && j > 0)
        offset.refB = offset.refB+j;

    return 0;
}

int dynProg::dpBandFull( const string& ref, 
        const string&     query,
        boundaries& offset,
        const dp_type&    type,
        const outputType& oType,
        dp_output&  output, 
        int bandL, 
        int bandR, 
        bool print){
    L1 = offset.refE-offset.refB+1;
    L2 = offset.queryE-offset.queryB+1;
    bandSize = 0;
    banded = true;
    bandLeft = bandL;
    bandRight = bandR;
    assert(L1 >=0 && L2 >= 0 && bandLeft >0 && bandRight > 0);
    //build matrices
    updateMatrix(type);
    bool affine = scores.openGap != 0;
    //dp
    if(affine)
        dpFillStatic( ref, query, true, offset, type);
    else
        dpFillOptStatic( ref, query, true, offset, type);
    //traceback and output
    int        i = min(L2,L1+bandLeft), j = min(L1,L2+bandRight);
    //find beginPosition for traceBack and find dpScore
    findTraceBackPosStatic(true, &i,&j,type);
    output.dpScore = M[i][j];
    //Next is for output only
    if(oType>DPSCORE){
        //Alter boundaries to return the new positions when free gaps were introduced.
        if(type.freeQueryE && i < L2)
            offset.queryE = offset.queryB+i-1;
        if(type.freeRefE && j < L1)
            offset.refE = offset.refB+j-1;
        //traceBack
        std::stringstream ss;
        if(affine)
            dpTraceBackStatic( i, j, type, output, ss, offset, true, ref, query);
        else
            dpTraceBackOptStatic( i, j, type, output, ss, offset, true, ref, query);
        //Edit Dist ouput
        if( i > 0 && !type.freeQueryB)
            output.editDist += i;
        if(j > 0 && !type.freeRefB)
            output.editDist += j;
        //Errorstring output
        if(oType == ERRORSTRING || oType == ALL){
            string errorString = ss.str();
            if(i > 0 && !type.freeQueryB)
                errorString.append(i,'I');
            else if(j > 0 && !type.freeRefB)
                errorString.append(j, 'D');
            reverse( errorString.begin(), errorString.end() );
            //create output: errorString
            int iter = 0;
            int temp;
            int errorStringLength = errorString.length();
            while(iter < errorStringLength){
                temp = iter;
                while(iter+1 < errorStringLength && errorString[iter+1]==errorString[temp])
                    iter++;
                iter++;
                output.cigarChars.push_back(errorString[temp]);
                output.cigarLengths.push_back(iter-temp);
            }
        }
    }
    //change boundaries for local and freeGap types to return changed free boundaries.
    if(type.freeQueryB && i > 0)
        offset.queryB = offset.queryB+i;
    if(type.freeRefB && j > 0)
        offset.refB = offset.refB+j;

    return 0;
}

