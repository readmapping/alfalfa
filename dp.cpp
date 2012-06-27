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

struct dp_matrices{
    dp_matrices(): L1(0), L2(0), bandSize(0), banded(false) {}
    int ** M;
    char ** traceBack;
    int ** gapUp;
    int ** gapLeft;
    int L1;
    int L2;
    int bandSize;
    bool banded;
};

bool insideBand(int i, int j, int band){
    return i-j <= band || i-j >= -1*band;
}

int initializeMatrix(dp_matrices& matrices,const dp_scores& scores,const dp_type& type){
    //TODO: change this to mallocs and initialize once per thread + realloc if nec.
    assert(matrices.L1 >=0 && matrices.L2 >= 0);
    matrices.M = new int * [matrices.L2+1];
    matrices.traceBack = new char * [matrices.L2+1];
    for(int i = 0; i <= matrices.L2; i++){
        matrices.M[i] = new int[matrices.L1+1];
        matrices.traceBack[i] = new char[matrices.L1+1];
    }
    //init matrices
    matrices.M[0][0] = 0;
    matrices.traceBack[0][0] = '0';//n is where traceBack stops
    int L2Bound = min(matrices.bandSize, matrices.L2);
    int L1Bound = min(matrices.bandSize, matrices.L1);
    if(matrices.bandSize < matrices.L2 || matrices.bandSize < matrices.L1)
        matrices.banded = true;
    else
        matrices.banded = false;
    for(int i = 1; i <= L2Bound; i++)
        matrices.traceBack[i][0] = '0';
    for(int i = 1; i <= L1Bound; i++)
        matrices.traceBack[0][i] = '0';
    assert(!banded || (!type.freeQueryB && !type.freeRefB));
    if(type.freeQueryB)
        for(int i = 1; i <= L2Bound; i++)
            matrices.M[i][0] = 0;
    else
        for(int i = 1; i <= L2Bound; i++)
            matrices.M[i][0] = scores.openGap + i*scores.extendGap;
    if(type.freeRefB)
        for(int i = 1; i <= L1Bound; i++)
            matrices.M[0][i] = 0;
    else
        for(int i = 1; i <= L1Bound; i++)
            matrices.M[0][i] = scores.openGap + i*scores.extendGap;
}

//only called for affine gap penalties
int initializeMatrices(dp_matrices& matrices,const dp_scores& scores,const dp_type& type){
    initializeMatrix(matrices, scores, type);
    matrices.gapUp = new int * [matrices.L2+1];
    matrices.gapLeft = new int * [matrices.L2+1];
    int L2Bound = min(matrices.bandSize, matrices.L2);
    int L1Bound = min(matrices.bandSize, matrices.L1);
    for(int i = 0; i <= matrices.L2; i++){
        matrices.gapUp[i] = new int[matrices.L1+1];
        matrices.gapLeft[i] = new int[matrices.L1+1];
    }
    matrices.gapLeft[0][0] = matrices.gapUp[0][0] = 0;
    for(int i = 1; i <= L1Bound; i++){
        matrices.gapUp[0][i] = matrices.M[0][i] + scores.openGap;
        matrices.gapLeft[0][i] = matrices.M[0][i] + scores.openGap;
    }
    for(int i = 1; i <= L2Bound; i++){
        matrices.gapUp[i][0] = matrices.M[i][0] + scores.openGap;
        matrices.gapLeft[i][0] = matrices.M[i][0] + scores.openGap;
    }
}

int deleteMatrices(dp_matrices& matrices,const dp_scores& scores,const outputType& oType){
    for(int  i = 0; i <= matrices.L2; i++ ){
        delete[] matrices.M[ i ];
    }
    delete[] matrices.M;
    if(oType != DPSCORE){
        for(int  i = 0; i <= matrices.L2; i++ ){
            delete[] matrices.traceBack[ i ];
        }
        delete[] matrices.traceBack;
    }
    if(scores.openGap!=0){
        for(int  i = 0; i <= matrices.L2; i++ ){
            delete[] matrices.gapUp[ i ];
            delete[] matrices.gapLeft[ i ];
        }
        delete[] matrices.gapUp;
        delete[] matrices.gapLeft;
    }
    else{
        matrices.gapUp = NULL;
        matrices.gapLeft = NULL;
    }
    return 0;
}

int deleteMatrix(dp_matrices& matrices){
    for(int  i = 0; i <= matrices.L2; i++ ){
        delete[] matrices.M[ i ];
        delete[] matrices.traceBack[ i ];
    }
    delete[] matrices.M;
    delete[] matrices.traceBack;
    matrices.gapUp = NULL;
    matrices.gapLeft = NULL;
    return 0;
}

int findTraceBackPos(const dp_matrices& matrices, int* const i, int* const j,const dp_type& type){
    //find beginPosition for traceBack and find dpScore
    int maximum = matrices.M[*i][*j];
    if(type.freeQueryE && type.freeRefE || type.local){//local
        if(!matrices.banded){
            for(int idx = matrices.L2; idx >= 0; idx--){
                for(int idx2 = matrices.L1; idx2 >= 0; idx2--){
                    if(matrices.M[idx][idx2] > maximum){
                        maximum = matrices.M[idx][idx2];
                        *i = idx;
                        *j = idx2;
                    }
                }
            }
        }
        else{
            for(int idx = matrices.L2; idx >= 0; idx--){
                for(int idx2 = matrices.L1; idx2 >= 0; idx2--){
                    if(insideBand(idx, idx2, matrices.bandSize) && matrices.M[idx][idx2] > maximum){
                        maximum = matrices.M[idx][idx2];
                        *i = idx;
                        *j = idx2;
                    }
                }
            }
        }
    }
    else if(type.freeQueryE){//banded update
        for(int idx = min(matrices.L2, matrices.L1 + matrices.bandSize); 
                idx >= min(matrices.L2, max(0,matrices.L1 - matrices.bandSize)); 
                idx--){
            if(matrices.M[idx][*j] > maximum){
                maximum = matrices.M[idx][*j];
                *i = idx;
            }
        }
    }
    else if(type.freeRefE){
        for(int idx = min(matrices.L1, matrices.L2 + matrices.bandSize); 
                idx >= min(matrices.L1, max(0,matrices.L2 - matrices.bandSize)); 
                idx--){
            if(matrices.M[*i][idx] > maximum){
                maximum = matrices.M[*i][idx];
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

int dpFill(dp_matrices& matrices,const string& ref,const string& query, bool forward,
        const boundaries& offset, const dp_scores& scores, 
        const dp_type& type){
    int        d = 0;
    int        fU, fD, fL ;
    char       ptr;

    for(int i = 1; i <= matrices.L2; i++ ){
        int idx2 = i-matrices.bandSize;
        if(idx2>0){//left value undefined
            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ];
            if(ref[offset.refB+idx2-1]==query[offset.queryB+i-1]){//TODO: remove if by using table
                matrices.M[ i ][ idx2 ] = scores.match;
                ptr = '\\';
            }
            else{
                matrices.M[ i ][ idx2 ] = ref[offset.refB+idx2-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
                ptr = ':';
                //remove this if with a table entry: init table!!!
            }
            matrices.gapUp[i][idx2] = max(matrices.M[ i-1 ][ idx2 ] + scores.extendGap + scores.openGap, matrices.gapUp[i-1][idx2]+scores.extendGap);
            if(matrices.M[ i ][ idx2 ] < matrices.gapUp[i][idx2]){
                matrices.M[ i ][ idx2 ] = matrices.gapUp[i][idx2];
                ptr = '|';
            }
            //fill in generic value for gapLeft, such that it is defined, but not chosen in next step
            matrices.gapLeft[i][idx2] = matrices.M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
                matrices.M[i][idx2] = 0;
                matrices.traceBack[i][idx2] = '0';
            }
            idx2++;
        }
        for(int j = max( 1, idx2); j <= min(matrices.L1, i+matrices.bandSize-1); j++ ){
            //here insert substitution matrix!
            if(ref[offset.refB+j-1]==query[offset.queryB+i-1])
                d = scores.match;
            else{
                d = ref[offset.refB+j-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
            }
            fU = max(matrices.M[ i-1 ][ j ] + scores.extendGap + scores.openGap, matrices.gapUp[i-1][j]+scores.extendGap);
            fD = matrices.M[ i-1 ][ j-1 ] + d;
            fL = max(matrices.M[ i ][ j-1 ] + scores.extendGap + scores.openGap, matrices.gapLeft[i][j-1]+scores.extendGap);
            matrices.gapUp[i][j] = fU;
            matrices.gapLeft[i][j] = fL;
            matrices.M[ i ][ j ] = maximum( fU, fD, fL, &ptr );
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][j]<0){
                matrices.M[i][j] = 0;
                matrices.traceBack[i][j] = '0';
            }
            else
                matrices.traceBack[i][j] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
        }
        if(i+matrices.bandSize <= matrices.L1){
            idx2 = i + matrices.bandSize;
            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ];
            if(ref[offset.refB+idx2-1]==query[offset.queryB+i-1]){//TODO: remove if by using table
                matrices.M[ i ][ idx2 ] += scores.match;
                ptr = '\\';
            }
            else{
                matrices.M[ i ][ idx2 ] += ref[offset.refB+idx2-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
                ptr = ':';
                //remove this if with a table entry: init table!!!
            }
            matrices.gapLeft[i][idx2] = max(matrices.M[ i ][ idx2-1 ] + scores.extendGap + scores.openGap, matrices.gapLeft[i][idx2-1]+scores.extendGap);
            if(matrices.M[ i ][ idx2 ] < matrices.gapLeft[i][idx2]){
                matrices.M[ i ][ idx2 ] = matrices.gapLeft[i][idx2];
                ptr = '-';
            }
            //fill in generic value for gapUp, such that it is defined, but not chosen in next step
            matrices.gapUp[i][idx2] = matrices.M[ i ][ idx2 ] + scores.extendGap + 2*scores.openGap;
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
                matrices.M[i][idx2] = 0;
                matrices.traceBack[i][idx2] = '0';
            }
        }
    }
    return 0;
}

int dpFillOpt(dp_matrices& matrices,const string& ref,const string& query, bool forward, 
        const boundaries& offset, const dp_scores& scores, const dp_type& type){
    //TODO: when changed to mallocs: use row-order traversal
    int        d = 0;
    char       ptr;
    
    if(forward)
    for(int i = 1; i <= matrices.L2; i++ ){
        int idx2 = i-matrices.bandSize;
        if(idx2>0){
            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ];
            if(ref[offset.refB+idx2-1]==query[offset.queryB+i-1]){//TODO: remove if by using table
                matrices.M[ i ][ idx2 ] += scores.match;
                ptr = '\\';
            }
            else{
                matrices.M[ i ][ idx2 ] += ref[offset.refB+idx2-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
                ptr = ':';
                //remove this if with a table entry: init table!!!
            }
            if(matrices.M[ i ][ idx2 ] < matrices.M[ i-1 ][ idx2 ] + scores.extendGap){
                matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2 ] + scores.extendGap;
                ptr = '|';
            }
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
                matrices.M[i][idx2] = 0;
                matrices.traceBack[i][idx2] = '0';
            }
            idx2++;
        }
        for(int j = max( 1, idx2); j <= min(matrices.L1, i+matrices.bandSize-1); j++ ){
            //here insert substitution matrix!
            if(ref[offset.refB+j-1]==query[offset.queryB+i-1])//TODO: remove if by using table
                d = scores.match;
            else{
                d = ref[offset.refB+j-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
                //remove this if with a table entry: init table!!!
            }
            matrices.M[ i ][ j ] = maximum( matrices.M[ i-1 ][ j ] + scores.extendGap,
                                matrices.M[ i-1 ][ j-1 ] + d,
                                matrices.M[ i ][ j-1 ] + scores.extendGap,
                                &ptr );
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][j]<0){//TODO: put this if up front for sw programming
                matrices.M[i][j] = 0;
                matrices.traceBack[i][j] = '0';
            }
            else
                matrices.traceBack[i][j] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
            //TODO: remove if by using table on previous match/mismatch table and storing the pointer in the max function
        }
        if(i+matrices.bandSize <= matrices.L1){
            idx2 = i + matrices.bandSize;
            matrices.M[ i ][ idx2 ] = matrices.M[ i-1 ][ idx2-1 ];
            if(ref[offset.refB+idx2-1]==query[offset.queryB+i-1]){//TODO: remove if by using table
                matrices.M[ i ][ idx2 ] += scores.match;
                ptr = '\\';
            }
            else{
                matrices.M[ i ][ idx2 ] += ref[offset.refB+idx2-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
                ptr = ':';
                //remove this if with a table entry: init table!!!
            }
            if(matrices.M[ i ][ idx2 ] < matrices.M[ i ][ idx2-1 ] + scores.extendGap){
                matrices.M[ i ][ idx2 ] = matrices.M[ i ][ idx2-1 ] + scores.extendGap;
                ptr = '-';
            }
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][idx2]<0){//TODO: put this if up front for sw programming
                matrices.M[i][idx2] = 0;
                matrices.traceBack[i][idx2] = '0';
            }
        }
    }
    return 0;
}

int dpTraceBack(const dp_matrices& matrices, int& i, int& j, dp_output& output, std::stringstream & ss){
    //TODO: after mallocs: change this traversal to row-order
    //change stringstream to char array of length dimensions array
    while( matrices.traceBack[ i ][ j ] != '0' ){
        switch( matrices.traceBack[ i ][ j ] ){
            case '|' :      ss << 'I';
                            output.editDist++;
                            i-- ;
                            break ;

            case '\\':      ss << '=';
                            i-- ;  j-- ;
                            break ;

            case ':':       ss << 'X';
                            output.editDist++;
                            i-- ;  j-- ;
                            break ;

            case '-' :      ss << 'D';
                            output.editDist++;
                            j-- ;
        }
    }
    return 0;
}

void  print_matrices(const dp_matrices& matrices,const string& ref, 
        const string& query, const boundaries& offset, bool gapMatrices )
{
        cout << "        ";
        for( int j = 0; j < matrices.L1; j++ )
        {
                cout << ref[ offset.refB+j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= matrices.L2; i++ )
        {
                if( i > 0 )
                {
                        cout << query[ offset.queryB+i-1 ] << " ";
                }
                for( int j = 0; j <= matrices.L1; j++ )
                {
                        cout.width( 3 );
                        cout << matrices.M[ i ][ j ] << " ";
                }
                cout << endl;
        }
        cout << endl;
        if(gapMatrices){
            cout << "        ";
            for( int j = 0; j < matrices.L1; j++ )
            {
                    cout << ref[ offset.refB+j ] << "   ";
            }
            cout << "\n  ";

            for( int i = 0; i <= matrices.L2; i++ )
            {
                    if( i > 0 )
                    {
                            cout << query[ offset.queryB+i-1 ] << " ";
                    }
                    for( int j = 0; j <= matrices.L1; j++ )
                    {
                            cout.width( 3 );
                            cout << matrices.gapLeft[ i ][ j ] << " ";
                    }
                    cout << endl;
            }
            cout << endl;
            cout << "        ";
            for( int j = 0; j < matrices.L1; j++ )
            {
                    cout << ref[ offset.refB+j ] << "   ";
            }
            cout << "\n  ";

            for( int i = 0; i <= matrices.L2; i++ )
            {
                    if( i > 0 )
                    {
                            cout << query[ offset.queryB+i-1 ] << " ";
                    }
                    for( int j = 0; j <= matrices.L1; j++ )
                    {
                            cout.width( 3 );
                            cout << matrices.gapUp[ i ][ j ] << " ";
                    }
                    cout << endl;
            }
            cout << endl;
        }
        cout << "        ";
        for( int j = 0; j < matrices.L1; j++ )
        {
                cout << ref[ offset.refB+j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= matrices.L2; i++ )
        {
                if( i > 0 )
                {
                        cout << query[ offset.queryB+i-1 ] << " ";
                }
                for( int j = 0; j <= matrices.L1; j++ )
                {
                        cout.width( 3 );
                        cout << matrices.traceBack[ i ][ j ] << " ";
                }
                cout << endl;
        }
        cout << endl;
}

void  print_seq( const string& ref,const string& query, boundaries& offset )
{
        cout << "" << endl << ref.substr(offset.refB,offset.refE-offset.refB+1) << endl;
        cout << query.substr(offset.queryB,offset.queryE-offset.queryB+1) << endl << "" << endl;
}

//(extra parameters offset1 and offset2, length1 and length2 in struct boundaries)
//extra: substitution matrix!
int dp( const string&     ref,
        const string&     query,
        boundaries& offset,
        const dp_scores&  scores,
        const dp_type&    type,
        const outputType& oType,
        dp_output&  output,
        bool print)
{
    if(print) print_seq(ref, query, offset);
    dp_matrices matrices;
    matrices.L1 = offset.refE-offset.refB+1;
    matrices.L2 = offset.queryE-offset.queryB+1;
    matrices.bandSize = max(matrices.L1, matrices.L2) + 1;
    assert(matrices.L1 >=0 && matrices.L2 >= 0);
    if(scores.openGap != 0){
        //build matrices
        initializeMatrices(matrices, scores, type);
        //Check for affine gap: if else
        //dp
        dpFill(matrices, ref, query, true, offset, scores, type);
        if(print) print_matrices(matrices, ref, query, offset, scores.openGap!=0);
        //traceback and output
        int        i = matrices.L2, j = matrices.L1;
        //find beginPosition for traceBack and find dpScore
        findTraceBackPos(matrices,&i,&j,type);
        output.dpScore = matrices.M[i][j];
        //Next is for output only
        if(oType>DPSCORE){
            //alignment output
            //Alter boundaries to return the new positions when free gaps were introduced.
            if((type.local || type.freeQueryE) && i < matrices.L2)
                offset.queryE = offset.queryB+i-1;
            if((type.local || type.freeRefE) && j < matrices.L1)
                offset.refE = offset.refB+j-1;
            //traceBack
            std::stringstream ss;
            dpTraceBack(matrices, i, j, output, ss);
            //Edit Dist ouput
            if( i > 0 && !type.local && !type.freeQueryB)
                output.editDist += i;
            if(j > 0 && !type.local && !type.freeRefB)
                output.editDist += j;
            //Errorstring output
            if(oType == ERRORSTRING || oType == ALL){
                string errorString = ss.str();
                if(i > 0 && !type.local && !type.freeQueryB)
                    errorString.append(i,'I');
                else if(j > 0 && !type.local && !type.freeRefB)
                    errorString.append(j, 'D');
                reverse( errorString.begin(), errorString.end() );
                //create output: errorString
                int iter = 0;
                int temp;
                while(iter < errorString.length()){
                    temp = iter;
                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
                        iter++;
                    iter++;
                    output.cigarChars.push_back(errorString[temp]);
                    output.cigarLengths.push_back(iter-temp);
                }
            }
        }
        //change boundaries for local and freeGap types to return changed free boundaries.
        if((type.local || type.freeQueryB) && i > 0)
            offset.queryB = offset.queryB+i;
        if((type.local || type.freeRefB) && j > 0)
            offset.refB = offset.refB+j;
        deleteMatrices(matrices, scores, oType);
    }
    else{
        //build matrices
        initializeMatrix(matrices, scores, type);
        //dp
        dpFillOpt(matrices, ref, query, true, offset, scores, type);
        //traceback and output
        int        i = matrices.L2, j = matrices.L1;
        //find beginPosition for traceBack and find dpScore
        findTraceBackPos(matrices,&i,&j,type);
        output.dpScore = matrices.M[i][j];
        //Next is for output only
        if(oType>DPSCORE){
            //Alter boundaries to return the new positions when free gaps were introduced.
            if((type.local || type.freeQueryE) && i < matrices.L2)
                offset.queryE = offset.queryB+i-1;
            if((type.local || type.freeRefE) && j < matrices.L1)
                offset.refE = offset.refB+j-1;
            //traceBack
            std::stringstream ss;
            dpTraceBack(matrices, i, j, output, ss);
            //Edit Dist ouput
            if( i > 0 && !type.local && !type.freeQueryB)
                output.editDist += i;
            if(j > 0 && !type.local && !type.freeRefB)
                output.editDist += j;
            //Errorstring output
            if(oType == ERRORSTRING || oType == ALL){
                string errorString = ss.str();
                if(i > 0 && !type.local && !type.freeQueryB)
                    errorString.append(i,'I');
                else if(j > 0 && !type.local && !type.freeRefB)
                    errorString.append(j, 'D');
                reverse( errorString.begin(), errorString.end() );//TODO: leave this out, but use different iterator
                //create output: errorString
                int iter = 0;
                int temp;
                while(iter < errorString.length()){
                    temp = iter;
                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
                        iter++;
                    iter++;
                    output.cigarChars.push_back(errorString[temp]);
                    output.cigarLengths.push_back(iter-temp);
                }
            }
        }
        //change boundaries for local and freeGap types to return changed free boundaries.
        if((type.local || type.freeQueryB) && i > 0)
            offset.queryB = offset.queryB+i;
        if((type.local || type.freeRefB) && j > 0)
            offset.refB = offset.refB+j;
        deleteMatrix(matrices);
    }

    return 0;
}

//(extra parameters offset1 and offset2, length1 and length2 in struct boundaries)
//extra: substitution matrix!
//Banded alignment only for band sizes > 0 and non-local alignment!!!
int dpBand( const string&     ref,
        const string&     query,
        boundaries& offset,
        const dp_scores&  scores,
        const dp_type&    type,
        const outputType& oType,
        dp_output&  output,
        int bandSize,
        bool print)
{
    dp_matrices matrices;
    matrices.L1 = offset.refE-offset.refB+1;
    matrices.L2 = offset.queryE-offset.queryB+1;
    matrices.bandSize = bandSize;
    assert(matrices.L1 >=0 && matrices.L2 >= 0 && bandSize >0);
    assert(!dp_type.local);
    assert((!dp_type.freeQueryB && !dp_type.freeRefB) || (!dp_type.freeQueryE && !dp_type.freeRefE));
    //possible end of alignment should be accessible for traceback: 
    //range of allowed traceback startpositions should have a non-empty intersection
    //with the band, otherwise: give bad alignment and return
    if((!type.freeRefE && matrices.L2 + bandSize < matrices.L1) ||
            (!type.freeQueryE && matrices.L1 + bandSize < matrices.L2)){
        output.dpScore = scores.mismatch*query.length();
        output.editDist = query.length();
        return 1;
    }
    
    if(scores.openGap != 0){
        //build matrices
        initializeMatrices(matrices, scores, type);
        //Check for affine gap: if else
        //dp
        dpFill(matrices, ref, query, true, offset, scores, type);
        if(print) print_matrices(matrices, ref, query, offset, scores.openGap!=0);
        //traceback and output
        int        i = matrices.L2, j = matrices.L1;
        //find beginPosition for traceBack and find dpScore
        findTraceBackPos(matrices,&i,&j,type);
        output.dpScore = matrices.M[i][j];
        //Next is for output only
        if(oType>DPSCORE){
            //Alter boundaries to return the new positions when free gaps were introduced.
            if(type.freeQueryE && i < matrices.L2)
                offset.queryE = offset.queryB+i-1;
            if(type.freeRefE && j < matrices.L1)
                offset.refE = offset.refB+j-1;
            //traceBack
            std::stringstream ss;
            dpTraceBack(matrices, i, j, output, ss);
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
                while(iter < errorString.length()){
                    temp = iter;
                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
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
        deleteMatrices(matrices, scores, oType);
    }
    else{
        //build matrices
        initializeMatrix(matrices, scores, type);
        //dp
        dpFillOpt(matrices, ref, query, true, offset, scores, type);
        //traceback and output
        int        i = matrices.L2, j = matrices.L1;
        //find beginPosition for traceBack and find dpScore
        findTraceBackPos(matrices,&i,&j,type);
        output.dpScore = matrices.M[i][j];
        //Next is for output only
        if(oType>DPSCORE){
            //Alter boundaries to return the new positions when free gaps were introduced.
            if((type.local || type.freeQueryE) && i < matrices.L2)
                offset.queryE = offset.queryB+i-1;
            if((type.local || type.freeRefE) && j < matrices.L1)
                offset.refE = offset.refB+j-1;
            //traceBack
            std::stringstream ss;
            dpTraceBack(matrices, i, j, output, ss);
            //Edit Dist ouput
            if( i > 0 && !type.local && !type.freeQueryB)
                output.editDist += i;
            if(j > 0 && !type.local && !type.freeRefB)
                output.editDist += j;
            //Errorstring output
            if(oType == ERRORSTRING || oType == ALL){
                string errorString = ss.str();
                if(i > 0 && !type.local && !type.freeQueryB)
                    errorString.append(i,'I');
                else if(j > 0 && !type.local && !type.freeRefB)
                    errorString.append(j, 'D');
                reverse( errorString.begin(), errorString.end() );//TODO: leave this out, but use different iterator
                //create output: errorString
                int iter = 0;
                int temp;
                while(iter < errorString.length()){
                    temp = iter;
                    while(iter+1 < errorString.length() && errorString[iter+1]==errorString[temp])
                        iter++;
                    iter++;
                    output.cigarChars.push_back(errorString[temp]);
                    output.cigarLengths.push_back(iter-temp);
                }
            }
        }
        //change boundaries for local and freeGap types to return changed free boundaries.
        if((type.local || type.freeQueryB) && i > 0)
            offset.queryB = offset.queryB+i;
        if((type.local || type.freeRefB) && j > 0)
            offset.refB = offset.refB+j;
        deleteMatrix(matrices);
    }

    return 0;
}
