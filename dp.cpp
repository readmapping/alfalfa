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
#include <limits>
#include <limits.h>

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

int dynProg::updateMatrix(const dp_type& type, bool print){
    //TODO: change this to mallocs and initialize once per thread + realloc if nec.
    assert(L1 >=0 && L2 >= 0);
    int newDimL2 = DP_L2;
    while(L2 >= newDimL2){
        newDimL2 = 2*newDimL2;
    }
    int newDimL1 = DP_L1;
    while(L1 >= newDimL1){
        newDimL1 = 2*newDimL1;
    }
    if(newDimL2 > DP_L2 || newDimL1 > DP_L1){
        if(print) cerr << "resize DP matrix from " << DP_L2 << "x" << DP_L1 << " to " << newDimL2 << "x" << newDimL1 << endl;
        resizeDPMatrix(newDimL2, newDimL1, scores.openGap != 0);
        if(print) cerr << "resize succesful" << endl;
    }
    M[0][0] = 0;
        for(int i = 1; i <= L2; i++)
            M[i][0] = scores.openGap + i*scores.extendGap;
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

//--> possible bug ? do not take into account forward or reverse
//it was ok for non-local dp or when gaps are not allowed in both begin and end, which was the only case untill now
int dynProg::findTraceBackPosStatic(bool forward, int* const i, int* const j,const dp_type& type){
    //find beginPosition for traceBack and find dpScore
    int maximum = M[*i][*j];
    int bandCorrection = 0;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    if((type.freeQueryE && type.freeRefE) || (type.freeQueryB && type.freeRefB) || type.local){//local
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
    else if(type.freeQueryE || type.freeQueryB){//banded update
        for(int idx = min(L2, bandCorrection + L1 + bandR); 
                idx >= min(L2, max(0, bandCorrection + L1 - bandL)); 
                idx--){
            if(M[idx][*j] > maximum){
                maximum = M[idx][*j];
                *i = idx;
            }
        }
    }
    else if(type.freeRefE || type.freeRefB){
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

int dynProg::dpFillStatic(const string& ref,const string& query, bool forward,
        const boundaries& offset, const dp_type& type){
    if(!forward)
        return dpFillStaticBackward(ref, query, forward, offset, type);
    int        d = 0;
    int        local = type.local || (type.freeQueryB && type.freeRefB) ? 0 : 1;//TODO: remove because this will never happen (speed up)
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

int dynProg::dpFillStaticBackward(const string& ref,const string& query, bool forward,
        const boundaries& offset, const dp_type& type){
    int        d = 0;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    for(long i = 1; i <= L2; i++ ){
        //always backward: traverse from right to left!
        long idx2 = i-bandL;
        long colRB = i+bandR;
        if(idx2>0){//left value undefined
            d = scores.scoreMatrix[size_t (ref[offset.refE-idx2+1])][size_t (query[offset.queryE-i+1])];
            M[ i ][ idx2 ] = M[ i-1 ][ idx2-1 ] + d;
            UP[i][idx2] = max(M[ i-1 ][ idx2 ] + scores.extendGap + scores.openGap, UP[i-1][idx2]+scores.extendGap);
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ] , UP[i][idx2]);
            //fill in generic value for gapLeft, such that it is defined, but not chosen in next step
            LEFT[i][idx2] = M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
            idx2++;
        }
        for(long j = max( 1L, idx2); j <= min((long)L1, colRB-1); j++ ){
            //here insert substitution matrix!
            d = scores.scoreMatrix[size_t(ref[offset.refE-j+1])][size_t(query[offset.queryE-i+1])];
            UP[i][j] = max(M[ i-1 ][ j ] + scores.extendGap + scores.openGap, UP[i-1][j]+scores.extendGap);
            M[i][j] = M[ i-1 ][ j-1 ] + d;
            LEFT[i][j] = max(M[ i ][ j-1 ] + scores.extendGap + scores.openGap, LEFT[i][j-1]+scores.extendGap);
            M[i][j] = max(M[i][j], LEFT[i][j]);
            M[i][j] = max(M[i][j], UP[i][j]);
        }
        if(colRB <= L1){
            idx2 = colRB;
            d = scores.scoreMatrix[size_t(ref[offset.refE-idx2+1])][size_t(query[offset.queryE-i+1])];
            M[ i ][ idx2 ] = M[ i-1 ][ idx2-1 ] + d;
            LEFT[i][idx2] = max(M[ i ][ idx2-1 ] + scores.extendGap + scores.openGap, LEFT[i][idx2-1]+scores.extendGap);
            M[ i ][ idx2 ] = max(M[ i ][ idx2 ], LEFT[i][idx2]);
            //fill in generic value for gapUp, such that it is defined, but not chosen in next step
            UP[i][idx2] = M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
        }
    }
    return 0;
}


int dynProg::dpFillOptStatic(const string& ref,const string& query, bool forward, 
        const boundaries& offset, const dp_type& type){
    if(!forward)
        return dpFillOptStaticBackward(ref, query, forward, offset, type);
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

int dynProg::dpFillOptStaticBackward(const string& ref,const string& query, bool forward, 
        const boundaries& offset, const dp_type& type){
    //TODO: when changed to mallocs: use row-order traversal
    int        d = 0;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    for(long i = 1; i <= L2; i++ ){
        long idx2 = i-bandL;
        long colRB = i+bandR;
        if(idx2>0){
            d = scores.scoreMatrix[size_t(ref[offset.refE-idx2+1])][size_t(query[offset.queryE-i+1])];
            M[ i ][ idx2 ] = max(M[ i-1 ][ idx2-1 ] + d, M[ i-1 ][ idx2 ] + scores.extendGap);
            idx2++;
        }
        for(long j = max( 1L, idx2); j <= min((long)L1, colRB-1); j++ ){
            d = scores.scoreMatrix[size_t(ref[offset.refE-j+1])][size_t(query[offset.queryE-i+1])];
            M[ i ][ j ] = max(M[ i-1 ][ j-1 ] + d, M[ i ][ j-1 ] + scores.extendGap);
            M[ i ][ j ] = max(M[ i ][ j ], M[ i-1 ][ j ] + scores.extendGap);
        }
        if(colRB <= L1){
            idx2 = colRB;
            d = scores.scoreMatrix[size_t(ref[offset.refE-idx2+1])][size_t(query[offset.queryE-i+1])];
            M[ i ][ idx2 ] = max(M[ i-1 ][ idx2-1 ] + d, M[ i ][ idx2-1 ] + scores.extendGap);
        }
    }
    return 0;
}

int dynProg::dpTraceBackStatic(int& i, int& j, const dp_type& type, dp_output& output, 
        std::stringstream & ss, const boundaries& offset, bool forward, const string& ref, const string& query){
    //TODO: after mallocs: change this traversal to row-order
    //change stringstream to char array of length dimensions array
    int qBound = 0; int rBound = 0;
    int bandL = max(bandLeft, bandSize);
    int bandR = max(bandRight, bandSize);
    bool local = type.local || (type.freeQueryB && type.freeRefB);
    int bandCorrection = 0;
    while( (i > qBound || j > rBound)){
        bool equal = false;
        if(j > 0 && i > 0)
            equal = forward ? ref[offset.refB+j-1] == query[offset.queryB+i-1] : ref[offset.refE-j+1] == query[offset.queryE-i+1];
        if( j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.match && equal){
            ss << '=';
            i-- ;  j-- ;
        }
        else if(j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.mismatch){
            ss << 'X';
            output.editDist++;
            i-- ;  j-- ;
        }
        else if(j > 0 && j > i+bandCorrection - bandL && (M[ i ][ j ] == LEFT[ i ][ j ] || i==0)){
            ss << 'D';
            output.editDist++;
            j-- ;
        }
        else if(i > 0 && j < i+bandCorrection + bandR && (M[ i ][ j ] == UP[ i ][ j ] || j==0)){
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
    bool local = type.local || (type.freeQueryB && type.freeRefB);
    int bandCorrection = 0;
    while( (i > qBound || j > rBound || local)){
        bool equal = false;
        if(j > 0 && i > 0)
            equal = forward ? ref[offset.refB+j-1] == query[offset.queryB+i-1] : ref[offset.refE-j+1] == query[offset.queryE-i+1];
        if( j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.match && equal){
            ss << '=';
            i-- ;  j-- ;
        }
        else if(j > 0 && i > 0 && M[ i ][ j ] == M[ i-1 ][ j-1 ] + scores.mismatch){
            ss << 'X';
            output.editDist++;
            i-- ;  j-- ;
        }
        else if(j > 0 && j > i+bandCorrection - bandL && (M[ i ][ j ] == M[ i ][ j-1 ] + scores.extendGap || i==0)){
            ss << 'D';
            output.editDist++;
            j-- ;
        }
        else if(i > 0 && j < i+bandCorrection + bandR && (M[ i ][ j ] == M[ i-1 ][ j ] + scores.extendGap || j==0)){
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

        cout << endl;
}

void  dynProg::print_seq( const string& ref,const string& query, boundaries& offset )
{
        cout << "" << endl << ref.substr(offset.refB,offset.refE-offset.refB+1) << endl;
        cout << query.substr(offset.queryB,offset.queryE-offset.queryB+1) << endl << "" << endl;
}


int dynProg::dpBandStatic( const string&     ref,
        const string&     query,
        boundaries& offset,
        const dp_type&    type,
        const outputType& oType,
        dp_output&  output,
        int minbandS,
        int maxbandS,
        bool print)
{
    L1 = offset.refE-offset.refB+1;
    L2 = offset.queryE-offset.queryB+1;
    bandSize = minbandS;
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
        if(max(L2,L1) - min(L2,L1) < maxbandS){
            bandSize = max(L2,L1) - min(L2,L1);
            if(print) cerr << "adjusted bandSize to " << bandSize << endl;
        }
        else{
            if(print) cerr << "possible end of alignment would lay beyond the dimensions of the matrix" << endl;
            output.dpScore = scores.mismatch*query.length();
            output.editDist = query.length();
            return 1;
        }
    }
    bool backward = (type.freeQueryB || type.freeRefB) && !type.local;
    //build matrices
    updateMatrix(type, print);
    bool affine = scores.openGap != 0;
    if(print)
        cerr << "fill DP matrix" << endl;
    //dp
    if(affine)
        dpFillStatic(ref, query, !backward, offset, type);
    else
        dpFillOptStatic(ref, query, !backward, offset, type);
    if(print){
        //print_matrices(ref, query, offset, false);
    }
    //traceback and output
    int        i = min(L2,L1+bandSize), j = min(L1,L2+bandSize);
    if(print)
        cerr << "find traceback position" << endl;
    //find beginPosition for traceBack and find dpScore
    findTraceBackPosStatic(!backward, &i,&j,type);
    if(print)
        cerr << "tracebackpos is " << i << ", " << j << endl;
    output.dpScore = M[i][j];
    if(print)
        cerr << "dp score is " << output.dpScore << endl;
    //Next is for output only
    if(oType>DPSCORE){
        //Alter boundaries to return the new positions when free gaps were introduced.
        if(type.freeQueryE && i < L2)
            offset.queryE = offset.queryB+i-1;
        if(type.freeRefE && j < L1)
            offset.refE = offset.refB+j-1;
        if(type.freeQueryB && i < L2)
            offset.queryB = offset.queryE-i+1;
        if(type.freeRefB && j < L1)
            offset.refB = offset.refE-j+1;
        //traceBack
        if(print)
                cerr << "search trace" << endl;
        std::stringstream ss;
        if(affine)
            dpTraceBackStatic(i, j, type, output, ss, offset, !backward, ref, query);
        else
            dpTraceBackOptStatic(i, j, type, output, ss, offset, !backward, ref, query);
        if(print)
                cerr << "trace ended at position " << i << ", " << j << endl;
        assert(j ==0 && i==0);
        //Errorstring output
        if(oType == ERRORSTRING || oType == ALL){
            if(print)
                cerr << "write errorstring to cigar vector" << endl;
            string errorString = ss.str();
            if(!backward)
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
    if(print)
        cerr << "DP finished" << endl;

    return 0;
}

//small pattern: no banded if band is as large as query

int dynProg::dpMyers( const string&     ref,
        const string&     query,
        boundaries& offset,
        const dp_type&    type,
        const outputType& oType,
        dp_output&  output,
        int minbandS,
        int maxbandS,
        bool print)
{
    L1 = offset.refE-offset.refB+1;
    L2 = offset.queryE-offset.queryB+1;
    bandSize = minbandS;
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
        if(max(L2,L1) - min(L2,L1) < maxbandS){
            bandSize = max(L2,L1) - min(L2,L1);
            if(print) cerr << "adjusted bandSize to " << bandSize << endl;
        }
        else{
            if(print) cerr << "possible end of alignment would lay beyond the dimensions of the matrix" << endl;
            output.dpScore = scores.mismatch*query.length();
            output.editDist = query.length();
            return 1;
        }
    }
    bool backward = (type.freeQueryB || type.freeRefB) && !type.local;
    //build matrices
    
    //MYERS STARTS HERE
    // Use size of unsigned int as blocksize for bit-vectors.
    const unsigned int BLOCK_SIZE = 32;//ulong size
    unsigned int len_x = abs(L1);
    unsigned int len_y = abs(L2);
    
    //PARAMETERS MASAI PAPER
    //m = L2
    //n = L1
    unsigned int w = 2*bandSize+1;
    unsigned int c = bandSize;
    
    unsigned int pos = 0;
    int score = (-1)*len_y;//Set to c-parameter
    unsigned int patternAlphabetSize = 5;
    unsigned int blockCount = (len_y + BLOCK_SIZE - 1) / BLOCK_SIZE;
    unsigned int scoreMask = 1 << ((len_y + BLOCK_SIZE -1)% BLOCK_SIZE);	// the mask with a bit set at the position of the last active cell
    vector<unsigned int> VP;
    VP.resize(blockCount, UINT_MAX);
    vector<unsigned int> VN;
    VN.resize(blockCount, 0);
    vector<unsigned> bitMask;
    bitMask.resize(patternAlphabetSize * blockCount, 0);
    
    // encoding the letters as bit-vectors
    for (unsigned int j = 0; j < len_y; j++){
            bitMask[blockCount * ORDVALUE[query[offset.queryB+j]] + j/BLOCK_SIZE] |= 1 << (j%BLOCK_SIZE);
            bitMask[blockCount * 4 + j/BLOCK_SIZE] &= ~(1 << (j%BLOCK_SIZE));
    }
    
    // compute score
    unsigned int X, D0, HN, HP;
    if(blockCount == 1){
        while (pos < len_x) {
            X = bitMask[ORDVALUE[ref[offset.refB+pos]]] | VN[0];

            D0 = ((VP[0] + (X & VP[0])) ^ VP[0]) | X;
            HN = VP[0] & D0;
            HP = VN[0] | ~(VP[0] | D0);

            // customized to compute edit distance
            X = (HP << 1) | 1;
            VN[0] = X & D0;
            VP[0] = (HN << 1) | ~(X | D0);

            if (HP & scoreMask)
                    score--;
            else if (HN & scoreMask)
                    score++;

            ++pos;
        }
    }// end compute score - short pattern
    else{
        unsigned int temp, shift, currentBlock;
        unsigned int carryD0, carryHP, carryHN;
        while (pos < len_x) {
            // set vars
            carryD0 = carryHP = carryHN = 0;
            shift = blockCount * ORDVALUE[ref[offset.refB+pos]];
            
            // computing first the top most block
            X = bitMask[shift] | VN[0];
	
            temp = VP[0] + (X & VP[0]);
            carryD0 = temp < VP[0];

            D0 = (temp ^ VP[0]) | X;
            HN = VP[0] & D0;
            HP = VN[0] | ~(VP[0] | D0);

            // customized to compute edit distance
            X = (HP << 1) | 1;
            carryHP = HP >> (BLOCK_SIZE - 1);

            VN[0] = X & D0;

            temp = (HN << 1);
            carryHN = HN >> (BLOCK_SIZE - 1);

            VP[0] = temp | ~(X | D0);

            // computing the necessary blocks, carries between blocks following one another are stored
            for (currentBlock = 1; currentBlock < blockCount; currentBlock++) {
                    X = bitMask[shift + currentBlock] | VN[currentBlock];

                    temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;

                    carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);

                    D0 = (temp ^ VP[currentBlock]) | X;
                    HN = VP[currentBlock] & D0;
                    HP = VN[currentBlock] | ~(VP[currentBlock] | D0);

                    X = (HP << 1) | carryHP;
                    carryHP = HP >> (BLOCK_SIZE-1);

                    VN[currentBlock] = X & D0;

                    temp = (HN << 1) | carryHN;
                    carryHN = HN >> (BLOCK_SIZE - 1);

                    VP[currentBlock] = temp | ~(X | D0);
            }

            // update score with the HP and HN values of the last block the last block
            if (HP & scoreMask)
                    score--;
            else if (HN & scoreMask)
                    score++;
            ++pos;
        }

    }// end compute score - long pattern
    output.editDist = score;
    return score;
    
    updateMatrix(type, print);
    bool affine = scores.openGap != 0;
    if(print)
        cerr << "fill DP matrix" << endl;
    //dp
    if(affine)
        dpFillStatic(ref, query, !backward, offset, type);
    else
        dpFillOptStatic(ref, query, !backward, offset, type);
    if(print){
        //print_matrices(ref, query, offset, false);
    }
    //traceback and output
    int        i = min(L2,L1+bandSize), j = min(L1,L2+bandSize);
    if(print)
        cerr << "find traceback position" << endl;
    //find beginPosition for traceBack and find dpScore
    findTraceBackPosStatic(!backward, &i,&j,type);
    if(print)
        cerr << "tracebackpos is " << i << ", " << j << endl;
    output.dpScore = M[i][j];
    if(print)
        cerr << "dp score is " << output.dpScore << endl;
    //Next is for output only
    if(oType>DPSCORE){
        //Alter boundaries to return the new positions when free gaps were introduced.
        if(type.freeQueryE && i < L2)
            offset.queryE = offset.queryB+i-1;
        if(type.freeRefE && j < L1)
            offset.refE = offset.refB+j-1;
        if(type.freeQueryB && i < L2)
            offset.queryB = offset.queryE-i+1;
        if(type.freeRefB && j < L1)
            offset.refB = offset.refE-j+1;
        //traceBack
        if(print)
                cerr << "search trace" << endl;
        std::stringstream ss;
        if(affine)
            dpTraceBackStatic(i, j, type, output, ss, offset, !backward, ref, query);
        else
            dpTraceBackOptStatic(i, j, type, output, ss, offset, !backward, ref, query);
        if(print)
                cerr << "trace ended at position " << i << ", " << j << endl;
        assert(j ==0 && i==0);
        //Errorstring output
        if(oType == ERRORSTRING || oType == ALL){
            if(print)
                cerr << "write errorstring to cigar vector" << endl;
            string errorString = ss.str();
            if(!backward)
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
    if(print)
        cerr << "DP finished" << endl;

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
    updateMatrix(type, print);
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

