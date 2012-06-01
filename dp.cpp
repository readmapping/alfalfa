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
    int ** M;
    char ** traceBack;
    int ** gapUp;
    int ** gapLeft;
    int L1;
    int L2;
};

int initializeMatrices(dp_matrices& matrices,const dp_scores& scores,const dp_type& type,const outputType& oType){
    assert(matrices.L1 >=0 && matrices.L2 >= 0);
    matrices.M = new int * [matrices.L2+1];
    for(int i = 0; i <= matrices.L2; i++)
        matrices.M[i] = new int[matrices.L1+1];
    //init matrices
    matrices.M[0][0] = 0;
    if(type.freeQueryB)
        for(int i = 1; i <= matrices.L2; i++)
            matrices.M[i][0] = 0;
    else
        for(int i = 1; i <= matrices.L2; i++)
            matrices.M[i][0] = scores.openGap+i*scores.extendGap;
    if(type.freeRefB)
        for(int i = 1; i <= matrices.L1; i++)
            matrices.M[0][i] = 0;
    else
        for(int i = 1; i <= matrices.L1; i++)
            matrices.M[0][i] = scores.openGap+i*scores.extendGap;
    if(oType != DPSCORE){
        matrices.traceBack = new char * [matrices.L2+1];
        for(int i = 0; i <= matrices.L2; i++)
            matrices.traceBack[i] = new char[matrices.L1+1];
        matrices.traceBack[0][0] = '0';//n is where traceBack stops
        for(int i = 1; i <= matrices.L2; i++)
            matrices.traceBack[i][0] = '0';
        for(int i = 1; i <= matrices.L1; i++)
            matrices.traceBack[0][i] = '0';
    }
    if(scores.openGap!=0){//affine gap penalty
        matrices.gapUp = new int * [matrices.L2+1];
        matrices.gapLeft = new int * [matrices.L2+1];
        for(int i = 0; i <= matrices.L2; i++){
            matrices.gapUp[i] = new int[matrices.L1+1];
            matrices.gapLeft[i] = new int[matrices.L1+1];
        }
        matrices.gapLeft[0][0] = matrices.gapUp[0][0] = 0;
        for(int i = 1; i <= matrices.L1; i++){
            matrices.gapUp[0][i] = matrices.M[0][i];
            matrices.gapLeft[0][i] = matrices.M[0][i];
        }
        for(int i = 1; i <= matrices.L2; i++){
            matrices.gapUp[i][0] = matrices.M[i][0];
            matrices.gapLeft[i][0] = matrices.M[i][0];
        }
    }
    else{
        matrices.gapUp = matrices.M;//now it should be a pointer to the same memory
        matrices.gapLeft = matrices.M;
    }
}

int initializeMatrix(dp_matrices& matrices,const dp_scores& scores,const dp_type& type){
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
    for(int i = 1; i <= matrices.L2; i++)
        matrices.traceBack[i][0] = '0';
    for(int i = 1; i <= matrices.L1; i++)
        matrices.traceBack[0][i] = '0';
    if(type.freeQueryB)
        for(int i = 1; i <= matrices.L2; i++)
            matrices.M[i][0] = 0;
    else
        for(int i = 1; i <= matrices.L2; i++)
            matrices.M[i][0] = i*scores.extendGap;
    if(type.freeRefB)
        for(int i = 1; i <= matrices.L1; i++)
            matrices.M[0][i] = 0;
    else
        for(int i = 1; i <= matrices.L1; i++)
            matrices.M[0][i] = i*scores.extendGap;
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
    int max = matrices.M[*i][*j];
    if(type.freeQueryE && type.freeRefE || type.local){//local
        for(int idx = matrices.L2; idx >= 0; idx--){
            for(int idx2 = matrices.L1; idx2 >= 0; idx2--){
                if(matrices.M[idx][idx2] > max){
                    max = matrices.M[idx][idx2];
                    *i = idx;
                    *j = idx2;
                }
            }
        }
    }
    else if(type.freeQueryE){
        for(int idx = matrices.L2; idx >= 0; idx--){
            if(matrices.M[idx][*j] > max){
                max = matrices.M[idx][*j];
                *i = idx;
            }
        }
    }
    else if(type.freeRefE){
        for(int idx = *j; idx >= 0; idx--){
            if(matrices.M[*i][idx] > max){
                max = matrices.M[*i][idx];
                *j = idx;
            }
        }
    }
    return 0;
}

int  maximum( int f1, int f2, int f3, char * ptr )
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

int dpFill(dp_matrices& matrices,const string& ref,const string& query,
        const boundaries& offset, const dp_scores& scores, 
        const dp_type& type,const outputType& oType){
    int        d = 0;
    int        fU, fD, fL ;
    char       ptr;

    for(int i = 1; i <= matrices.L2; i++ ){
        for(int j = 1; j <= matrices.L1; j++ ){
            //here insert substitution matrix!
            if(ref[offset.refB+j-1]==query[offset.queryB+i-1])
                d = scores.match;
            else{
                d = ref[offset.refB+j-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
            }
            fU = maximum(matrices.M[ i-1 ][ j ] + scores.extendGap + scores.openGap, matrices.gapUp[i-1][j]+scores.extendGap, matrices.gapUp[i-1][j]+scores.extendGap , &ptr);
            fD = matrices.M[ i-1 ][ j-1 ] + d;
            fL = maximum(matrices.M[ i ][ j-1 ] + scores.extendGap + scores.openGap, matrices.gapLeft[i][j-1]+scores.extendGap, matrices.gapLeft[i][j-1]+scores.extendGap , &ptr);
            matrices.gapUp[i][j] = fU;
            matrices.gapLeft[i][j] = fL;
            matrices.M[ i ][ j ] = maximum( fU, fD, fL, &ptr );
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][j]<0){
                matrices.M[i][j] = 0;
                matrices.traceBack[i][j] = '0';
            }
            else if(oType != DPSCORE){
                if(ptr=='\\' && d == scores.mismatch)
                    matrices.traceBack[i][j] = ':';
                else
                    matrices.traceBack[ i ][ j ] =  ptr ;
            }
        }
    }
    return 0;
}

int dpFillOpt(dp_matrices& matrices,const string& ref,const string& query, 
        const boundaries& offset, const dp_scores& scores, const dp_type& type){
    int        d = 0;
    char       ptr;

    for(int i = 1; i <= matrices.L2; i++ ){
        for(int j = 1; j <= matrices.L1; j++ ){
            //here insert substitution matrix!
            if(ref[offset.refB+j-1]==query[offset.queryB+i-1])
                d = scores.match;
            else{
                d = ref[offset.refB+j-1] == '`' ? scores.mismatch*query.length() : scores.mismatch;// if ref boundary is passed: huge penalty (A bit hacky)
            }
            matrices.M[ i ][ j ] = maximum( matrices.M[ i-1 ][ j ] + scores.extendGap,
                                matrices.M[ i-1 ][ j-1 ] + d,
                                matrices.M[ i ][ j-1 ] + scores.extendGap,
                                &ptr );
            if((type.local || (type.freeQueryB && type.freeRefB)) && matrices.M[i][j]<0){
                matrices.M[i][j] = 0;
                matrices.traceBack[i][j] = '0';
            }
            else
                matrices.traceBack[i][j] = (ptr=='\\' && d == scores.mismatch) ? ':' : ptr;
        }
    }
    return 0;
}

//int scoreOnly(char ref, char q, char err, dp_output& output){
//    return 0;
//}
//
//int alignmentOutput(char ref, char q, char err, dp_output& output){
//    output.traceRef += ref;
//    output.traceQuery += q;
//    return 0;
//}
//
//int errStringOutput(char ref, char q, char err, dp_output& output){
//    output.errorString += err;
//    return 0;
//}
//
//int allOutput(char ref, char q, char err, dp_output& output){
//    output.traceRef += ref;
//    output.traceQuery += q;
//    output.errorString += err;
//    return 0;
//}

//int dpTraceBack(dp_matrices& matrices, int& i, int& j, dp_output& output, string& ref, string& query, boundaries& offset, int(*outputSwitch)(char,char,char,dp_output&)){
//    while( matrices.traceBack[ i ][ j ] != '0' ){
//        switch( matrices.traceBack[ i ][ j ] ){
//            case '|' :      outputSwitch('-', query[offset.queryB+ i-1 ], 'I', output);
//                            output.editDist++;
//                            i-- ;
//                            break ;
//
//            case '\\':      outputSwitch(ref[offset.refB + j-1 ], query[offset.queryB+ i-1 ], '=', output);
//                            i-- ;  j-- ;
//                            break ;
//
//            case ':':       outputSwitch(ref[offset.refB + j-1 ], query[offset.queryB+ i-1 ], 'X', output);
//                            output.editDist++;
//                            i-- ;  j-- ;
//                            break ;
//
//            case '-' :      outputSwitch(ref[offset.refB +  j-1 ], '-', 'D', output);
//                            output.editDist++;
//                            j-- ;
//        }
//    }
//    return 0;
//}

int dpTraceBack(const dp_matrices& matrices, int& i, int& j, dp_output& output, std::stringstream & ss){
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
    assert(matrices.L1 >=0 && matrices.L2 >= 0);
    if(scores.openGap != 0){
        //build matrices
        initializeMatrices(matrices, scores, type, oType);
        //Check for affine gap: if else
        //dp
        dpFill(matrices, ref, query, offset, scores, type, oType);
        if(print) print_matrices(matrices, ref, query, offset, scores.openGap!=0);
        //traceback and output
        int        i = matrices.L2, j = matrices.L1;
        //find beginPosition for traceBack and find dpScore
        findTraceBackPos(matrices,&i,&j,type);
        output.dpScore = matrices.M[i][j];
        //Next is for output only
        if(oType>DPSCORE){
            //alignment output
    //        if(oType == ALIGNMENT || oType == ALL){
    //            output.traceRef = "";
    //            output.traceQuery = "";
    //            if(type.freeQueryE && !type.local && i < matrices.L2){
    //                output.traceRef.append(matrices.L2-i,'-');
    //                output.traceQuery.append(query.rbegin()+query.length()-offset.queryE-1,query.rbegin()+query.length()-offset.queryB-1);
    //            }
    //            else if(type.freeRefE && !type.local && j < matrices.L1){
    //                output.traceRef.append(ref.rbegin()+ref.length()-offset.refE-1,ref.rbegin()+ref.length()-offset.refB-1);
    //                output.traceQuery.append(matrices.L1-j,'-');
    //            }
    //        }
            //Alter boundaries to return the new positions when free gaps were introduced.
            if((type.local || type.freeQueryE) && i < matrices.L2)
                offset.queryE = offset.queryB+i-1;
            if((type.local || type.freeRefE) && j < matrices.L1)
                offset.refE = offset.refB+j-1;
            //traceBack
            std::stringstream ss;
            dpTraceBack(matrices, i, j, output, ss);
    //        switch(oType){
    //            case SCORE : dpTraceBack(matrices, i, j, output, ref, query, offset, scoreOnly); break;
    //            case ALIGNMENT : dpTraceBack(matrices, i, j, output, ref, query, offset, alignmentOutput); break;
    //            case ERRORSTRING : dpTraceBack(matrices, i, j, output, ref, query, offset, errStringOutput); break;
    //            case ALL : dpTraceBack(matrices, i, j, output, ref, query, offset, allOutput);
    //        }
            //Edit Dist ouput
            if( i > 0 && !type.local && !type.freeQueryB)
                output.editDist += i;
            if(j > 0 && !type.local && !type.freeRefB)
                output.editDist += j;
            //Alignment output
    //        if(oType == ALIGNMENT || oType == ALL){
    //            if(i > 0 && !type.local){
    //                output.traceRef.append(i,'-');
    //                output.traceQuery.append(query, offset.queryB,i);
    //            }
    //            else if(j > 0 && !type.local){
    //                output.traceRef.append(ref,offset.refB,j);
    //                output.traceQuery.append(j, '-');
    //            }
    //            reverse( output.traceRef.begin(), output.traceRef.end() );
    //            reverse( output.traceQuery.begin(), output.traceQuery.end() );
    //        }
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
        dpFillOpt(matrices, ref, query, offset, scores, type);
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
        deleteMatrix(matrices);
    }

    return 0;
}


