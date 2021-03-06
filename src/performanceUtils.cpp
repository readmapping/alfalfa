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

#include <iostream>

#include <cmath>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>

#include "performanceUtils.h"
#include "fasta.h"

using namespace std;

//////////////////////////
//struct oracleLine_t
//////////////////////////

void oracleLine_t::addAlignment(int score){
    tempReadResult = max(tempReadResult, score);
    alignmentCnt++;
    if(score == 2)
        alnGood++;
    else if(score == 1)
        alnOk++;
    else
        alnBad++;
}
void oracleLine_t::finishRead(){
    if(tempReadResult==2)
        readGood++;
    else if(tempReadResult ==1)
        readOk++;
    else
        readBad++;
    tempReadResult = 0;
}
string oracleLine_t::printLine(long readcnt){
    stringstream ss;
    ss << alnGood << "\t " << ((alnGood*100)/(alignmentCnt > 0 ? alignmentCnt : 1)) << "\t ";
    ss << alnOk << "\t " << ((alnOk*100)/(alignmentCnt > 0 ? alignmentCnt : 1)) << "\t "; 
    ss << alnBad << "\t " << ((alnBad*100)/(alignmentCnt > 0 ? alignmentCnt : 1)) << "\t ";
    ss << readGood << "\t " << ((readGood*100)/(readcnt > 0 ? readcnt : 1)) << "\t ";
    ss << readOk << "\t " << ((readOk*100)/(readcnt > 0 ? readcnt : 1)) << "\t ";
    ss << readBad << "\t " << ((readBad*100)/(readcnt > 0 ? readcnt : 1)) << endl;
    return ss.str();
}

//////////////////////////
//struct summaryLine_t
//////////////////////////

void summaryLine_t::addAlignment(int flag, bool mateFound){
    int score = 0;
    //set fields
    tempAlnCount++;
    alignmentCnt++;
    if(!(flag & (1<<2))){
        mappedAln++;
        score++;
        if(mateFound){
            mappedPaired++;
            score++;
            if(flag & (1<<1)){
                mappedPairedCorrect++;
                score++;
            }
        }
        if(flag & (1<<8)){
            secondAln++;
        }
    }
    tempReadResult = max(tempReadResult, score);
}
    
void summaryLine_t::finishRead(){
    if(tempReadResult==3)
        readPairedCorrect++;
    else if(tempReadResult ==2)
        readPaired++;
    else if(tempReadResult ==1)
        readMapped++;
    tempReadResult = 0;
    alnPerRead[min(tempAlnCount,101)]++;
    tempAlnCount = 0;
}

string summaryLine_t::printLine(long readcnt, bool paired){
    stringstream ss;
    int upperbound = 101;
    while(alnPerRead[upperbound]==0)
        upperbound--;
    if(paired){
        ss << "~~~\t num_mapped\t percent_mapped\t num_paired\t percent_paired\t /"
                "num_paired_correct\t percent_paired_correct\t num_unmapped\t /"
                "percent_unmapped\t num_second_aligned\t count" << endl;
        ss << "alignments\t " << mappedAln << "\t " << (mappedAln*100)/alignmentCnt;
        ss << "\t " << mappedPaired << "\t " << (mappedPaired*100)/alignmentCnt;
        ss << "\t " << mappedPairedCorrect << "\t " << (mappedPairedCorrect*100)/alignmentCnt;
        ss << "\t " << (alignmentCnt-mappedAln) << "\t " << ((alignmentCnt-mappedAln)*100)/alignmentCnt;
        ss << "\t " << secondAln << "\t " << alignmentCnt << endl;
        ss << "reads\t " << readMapped << "\t " << (readMapped*100)/readcnt;
        ss << "\t " << readPaired << "\t " << (readPaired*100)/readcnt;
        ss << "\t " << readPairedCorrect << "\t " << (readPairedCorrect*100)/readcnt;
        ss << "\t " << (readcnt-readMapped-readPaired-readPairedCorrect) << "\t " << ((readcnt-readMapped-readPaired-readPairedCorrect)*100)/readcnt;
        ss << "\t " << "N.A." << "\t " << readcnt << endl;
    }
    else{
        ss << "~~~\t num_mapped\t percent_mapped\t num_unmapped\t /"
                "percent_unmapped\t num_second_aligned\t count" << endl;
        ss << "alignments\t " << mappedAln << "\t " << (mappedAln*100)/(alignmentCnt > 0 ? alignmentCnt: 1 );
        ss << "\t " << (alignmentCnt-mappedAln) << "\t " << ((alignmentCnt-mappedAln)*100)/(alignmentCnt > 0 ? alignmentCnt: 1);
        ss << "\t " << secondAln << "\t " << alignmentCnt << endl;
        ss << "reads\t " << readMapped << "\t " << (readMapped*100)/(readcnt>0?readcnt:1);
        ss << "\t " << (readcnt-readMapped) << "\t " << ((readcnt-readMapped)*100)/(readcnt>0?readcnt:1);
        ss << "\t " << "N.A." << "\t " << readcnt << endl;
    }
    ss << endl;
    ss << "number of reads with x alignments" << endl;
    ss << 0 << "\t " << (readcnt-readMapped-readPaired-readPairedCorrect) << endl;
    for(int i = 1; i<= upperbound; i++)
        ss << i << "\t " << alnPerRead[i] << endl;
    return ss.str();
}
    
//////////////////////////
//static functions
//////////////////////////

template<class T>
string printVector(vector<T> vec){
    stringstream ss;
    if(!vec.empty()){
        ss << vec[0];
        for(size_t i=1; i < vec.size(); i++)
            ss << "," << vec[i];
    }
    return ss.str();
}

void processCommaSepListInt(string list, vector<int> & options, string name, int lowerBound, int upperBound){
    int beginPos = 0;
    int commaPos = list.find(',',beginPos);
    while(commaPos>=0){
        int value = atoi(list.substr(beginPos,commaPos-beginPos).c_str());
        if(value >= lowerBound && value <= upperBound)
            options.push_back(value);
        else{
            fprintf(stderr, "Value %d for parameter %s does not fall within logical bounds. This value will be ignored.\n", value, name.c_str());
        }
        beginPos = commaPos+1;
        commaPos = list.find(',',beginPos);
    }
    int value = atoi(list.substr(beginPos).c_str());
    if(value >= lowerBound && value <= upperBound)
        options.push_back(value);
    else{
        fprintf(stderr, "Value %d for parameter %s does not fall within logical bounds. This value will be ignored.\n", value, name.c_str());
    }
}

void processCheckParameters(int argc, char* argv[], samCheckOptions_t& opt){
    //parse the subcommand
    if(argc < 2 || strcmp(argv[1], "-h")==0 || strcmp(argv[1], "--help")==0){
        option::printUsage(std::cerr, evaluateUsage);
        exit(1);
    }
    string subcommand = argv[1];
    if(strcmp(argv[1], "sam") == 0) opt.subcommand = ORACLE;
    else if(strcmp(argv[1], "summary") == 0) opt.subcommand = SUMMARY;
    else if(strcmp(argv[1], "wgsim") == 0) opt.subcommand = WGSIM;
    else{
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        option::printUsage(std::cerr, evaluateUsage);
        exit(1);
    }
    cerr << "COMMAND: evaluate " << argv[1] << endl;
    argc-=2; argv++; argv++;
    //summary options
    if(opt.subcommand == SUMMARY){
        option::Stats  stats(summaryUsage, argc, argv);
        option::Option* options = new option::Option[stats.options_max];
        option::Option* buffer  = new option::Option[stats.buffer_max];
        option::Parser parser(summaryUsage, argc, argv, options, buffer);
        if (parser.error()){
            cerr << "error during parsing of command line options" << endl;
            exit(1);
        }
        if(parser.optionsCount() == 0){
            option::printUsage(std::cerr, summaryUsage); exit(0);
        }
        for (int i = 0; i < parser.optionsCount(); ++i) {
            option::Option& option = buffer[i];
            switch(option.index()) {
                case ARG_INPUT_SAM: opt.querySam = option.arg; break;
                case ARG_OUTPUT_EVAL: opt.outputFile = option.arg; break;
                case ARG_QUALITY: processCommaSepListInt(option.arg, opt.qualityValues, "quality-values", 0, 255); break;
                case ARG_READS: opt.numReads = atoi(option.arg); break;
                case ARG_PAIRED: opt.paired = true; break;
                case ARG_PRINT: opt.print = true; break;
                case ARG_HELP_EVAL: option::printUsage(std::cerr, summaryUsage); exit(0);
                default: cerr << "unknown option: " << option.name << endl; break;
                throw 1;
            }
        }
        delete[] options;
        delete[] buffer;
    }
    else if(opt.subcommand == ORACLE){
        option::Stats  stats(samUsage, argc, argv);
        option::Option* options = new option::Option[stats.options_max];
        option::Option* buffer  = new option::Option[stats.buffer_max];
        option::Parser parser(samUsage, argc, argv, options, buffer);
        if (parser.error()){
            cerr << "error during parsing of command line options" << endl;
            exit(1);
        }
        if(parser.optionsCount() == 0){
            option::printUsage(std::cerr, samUsage); exit(0);
        }
        for (int i = 0; i < parser.optionsCount(); ++i) {
            option::Option& option = buffer[i];
            switch(option.index()) {
                case ARG_INPUT_SAM: opt.querySam = option.arg; break;
                case ARG_OUTPUT_EVAL: opt.outputFile = option.arg; break;
                case ARG_QUALITY: processCommaSepListInt(option.arg, opt.qualityValues, "quality-values", 0, 255); break;
                case ARG_REFERENCE_EVAL: opt.refString = option.arg; break;
                case ARG_WINDOW: processCommaSepListInt(option.arg, opt.correctRange, "ranges", 0, 1000000); break;
                case ARG_INPUT_EDIT: opt.qtag = option.arg; break;
                case ARG_REFERENCE_SAM: opt.oracleSam = option.arg; break;
                case ARG_REFERENCE_EDIT: opt.otag = option.arg; break;
                case ARG_PRINT: opt.print = true; break;
                case ARG_HELP_EVAL: option::printUsage(std::cerr, samUsage); exit(0);
                default: cerr << "unknown option: " << option.name << endl; break;
                throw 1;
            }
        }
        delete[] options;
        delete[] buffer;
    }
    else if(opt.subcommand == WGSIM){
        option::Stats  stats(wgsimUsage, argc, argv);
        option::Option* options = new option::Option[stats.options_max];
        option::Option* buffer  = new option::Option[stats.buffer_max];
        option::Parser parser(wgsimUsage, argc, argv, options, buffer);
        if (parser.error()){
            cerr << "error during parsing of command line options" << endl;
            exit(1);
        }
        if(0==parser.optionsCount()){
            option::printUsage(std::cerr, wgsimUsage); exit(0);
        }
        for (int i = 0; i < parser.optionsCount(); ++i) {
            option::Option& option = buffer[i];
            switch(option.index()) {
                case ARG_INPUT_SAM: opt.querySam = option.arg; break;
                case ARG_OUTPUT_EVAL: opt.outputFile = option.arg; break;
                case ARG_QUALITY: processCommaSepListInt(option.arg, opt.qualityValues, "quality-values", 0, 255); break;
                case ARG_REFERENCE_EVAL: opt.refString = option.arg; break;
                case ARG_WINDOW: processCommaSepListInt(option.arg, opt.correctRange, "ranges", 0, 1000000); break;
                case ARG_INPUT_EDIT: opt.qtag = option.arg; break;
                case ARG_PRINT: opt.print = true; break;
                case ARG_HELP_EVAL: option::printUsage(std::cerr, wgsimUsage); exit(0);
                default: cerr << "unknown option: " << option.name << endl; break;
                throw 1;
            }
        }
        delete[] options;
        delete[] buffer;
    }
    if(opt.querySam.empty()){
        fprintf(stderr, "ERROR: The query SAM file has not been specified using the -i option.\n"); 
        exit(1);
    }
    if(opt.subcommand == SUMMARY){
        if(opt.numReads==0)
            fprintf(stderr, "WARNING: Number of reads was not specified. Number of reads in input sam will be used.\n");
    }
    if(opt.subcommand == ORACLE){
        if(opt.refString.empty()){
            fprintf(stderr, "ERROR: The reference genome has not been specified using the -r option.\n"); 
            exit(1);
        }
        if(opt.oracleSam.empty()){
            fprintf(stderr, "ERROR: The reference SAM file has not been specified using the --reference-sam option.\n"); 
            exit(1);
        }
    }
    if(opt.subcommand == WGSIM){
        if(opt.refString.empty()){
            fprintf(stderr, "ERROR: The reference genome has not been specified using the -r option.\n"); 
            exit(1);
        }
    }
    if(opt.subcommand != SUMMARY && opt.correctRange.empty())
        opt.correctRange.push_back(50);
}

string nextField(const string & line, char delimeter, int& beginPos){
    int tabPos = line.find(delimeter,beginPos);
    string substring = line.substr(beginPos, tabPos-beginPos);
    beginPos = tabPos+1;
    return substring;
}

string previousField(string & line, char delimeter, int& endPos){
    int tabPos = line.rfind(delimeter,endPos);
    string substring = line.substr(tabPos+1, endPos-tabPos);
    endPos = tabPos-1;
    return substring;
}

void getlijn2(ifstream & input, string & lijn){
    //write windows version of getline
    getline(input, lijn, '\n'); // Load one line at a time.
    if(lijn.length() > 0 && lijn[lijn.length()-1]=='\r')
        lijn.erase(--lijn.end());
}

bool readRecord(string & line, samRecord_t & record, string tag){
    int tabPos = 0;
    char delimeter = '\t';
    //qname
    record.qname = nextField(line, delimeter, tabPos);
    if(record.qname.length() > 2 && record.qname[record.qname.length()-2]=='/')
        record.qname.erase(record.qname.length()-2);
    //flag
    record.flag = bitset<11>((ulong) atoi(nextField(line, delimeter, tabPos).c_str()));
    //rname
    record.rname = nextField(line, delimeter, tabPos);
    //pos
    record.pos = atoi(nextField(line, delimeter, tabPos).c_str());
    //mapq
    record.mapq = atoi(nextField(line, delimeter, tabPos).c_str());
    //skip CIGAR
    string cigar = nextField(line, delimeter, tabPos);
    int cigarPos = 0;
    while(cigar[cigarPos] >= 48 && cigar[cigarPos] <= 57)//integer
        cigarPos++;
    //now cigar[cigarPos] is character
    if(cigar[cigarPos] == 'S' || cigar[cigarPos] == 'H'){
        record.pos -= atoi(cigar.substr(0,cigarPos).c_str());
    }
    //rnext
    record.rnext = nextField(line, delimeter, tabPos);
    //pnext
    record.pnext = atoi(nextField(line, delimeter, tabPos).c_str());
    //edit distance NM tag
    //skip fields to make sure 
    tabPos = line.find(delimeter,tabPos)+1;
    tabPos = line.find(delimeter,tabPos)+1;
    tabPos = line.find(delimeter,tabPos)+1;
    //now passed qual
    string editTag = tag+":i:";
    int NMpos = line.find(editTag,tabPos);
    if(NMpos == -1 && !record.flag.test(2))
            record.edit = -1;
    else if(NMpos>=tabPos){
        tabPos = line.find(delimeter,NMpos);
        if(tabPos == -1)
            tabPos = line.length();
        record.edit = atoi(line.substr(NMpos+5,tabPos-NMpos-5).c_str());
    }
    return true;
}

int unpairedCompareSam(samRecord_t& query, samRecord_t& oracle, int range){
    //correct if within X bases of origin
    if(query.flag.test(2)){
        return 0;//unmapped
    }
    if(query.flag.test(4)==oracle.flag.test(4) && 
            query.rname.compare(oracle.rname)==0 && 
            query.pos >= oracle.pos-range && query.pos <= oracle.pos + range){
        return 2;
    }
    else if(query.edit <= oracle.edit){
        return 1;
    }
    //ok if edit distance is fine
    else{
        return 0;
    }
    //else incorrect
}

int pairedCompareSam(samRecord_t& query, samRecord_t& oracle, int range){
    //correct if within X bases of origin
    if(query.flag.test(2)){
        return 0;//unmapped
    }
    if((query.flag.to_ulong() % 128) == (oracle.flag.to_ulong() % 128) &&
            query.rname.compare(oracle.rname)==0 && 
            query.pos >= oracle.pos-range && query.pos <= oracle.pos + range && 
            query.rnext.compare(oracle.rnext) == 0 ){// && query.pnext >= oracle.pnext-range && query.pnext <= oracle.pnext+range
        return 2;
    }
    else if(query.edit <= oracle.edit){
        return 1;
    }
    //ok if edit distance is fine
    else{
        return 0;
    }
    //else incorrect
}

int calculatedDPScore(const string& S, vector<string> & refdescr, vector<long> & startpos, const string& queryLine){
    //read boundaries and sequence from SAM line
    int tabPos = 0;
    char delimeter = '\t';
    //qname
    nextField(queryLine, delimeter, tabPos);
    //flag
    nextField(queryLine, delimeter, tabPos);
    //rname
    string rname = nextField(queryLine, delimeter, tabPos);
    long chromStart = 0;
    size_t index = 0;
    while(index < refdescr.size() && refdescr[index].compare(rname)!=0)
        index++;
    if(index < refdescr.size()){
        chromStart = startpos[index];
    }
    else{
        cerr << "query found in sequence " << rname << ". But this sequence is not found in reference string" << endl;
        cerr << queryLine << endl;
        cerr << "length of line: " << queryLine.length() << endl;
        exit(1);
    }
    //pos
    long refStart = atoi(nextField(queryLine, delimeter, tabPos).c_str()) -1;//-1 because SAM is 1-based
    refStart += chromStart;
    //mapq
    nextField(queryLine, delimeter, tabPos).c_str();
    //skip CIGAR
    string cigar = nextField(queryLine, delimeter, tabPos);
    //rnext
    nextField(queryLine, delimeter, tabPos);
    //pnext
    nextField(queryLine, delimeter, tabPos);
    //tlen
    nextField(queryLine, delimeter, tabPos);
    //seq
    string P = nextField(queryLine, delimeter, tabPos);
    for(size_t i = 0; i < P.length(); i++){
        P[i] = std::tolower(P[i]);
    }
    int cigarPos = 0;
    long refEnd = refStart;
    long qStart = 0;
    int editDist = 0;
    while(cigarPos < (int)cigar.length()){
        int beginPos = cigarPos;
        while(cigar[cigarPos] >= 48 && cigar[cigarPos] <= 57)//integer
                cigarPos++;
        //now cigar[cigarPos] is character
        long partsize = atoi(cigar.substr(beginPos,cigarPos-beginPos).c_str());
        if(cigar[cigarPos] == 'D' || cigar[cigarPos] == 'X' || cigar[cigarPos] == 'I'){
            editDist += partsize;
            qStart += partsize;
            refEnd += partsize;
        }
        else if(cigar[cigarPos] == 'S' || cigar[cigarPos] == 'H'){
            qStart += partsize;
        }
        else if(cigar[cigarPos] == 'M'){
            
            for(int i= 0; i < partsize; i++){
                if(P[qStart+i]!=S[refEnd+i])
                    editDist++;
            }
            qStart += partsize;
            refEnd += partsize;
        }
        cigarPos++;
    }
    return editDist;
}

void checkOracle(samCheckOptions_t & opt){
    //initialize fields and data for calculating edit distance
    string ref;
    vector<string> refdescr;
    vector<long> startpos;
    load_fasta(opt.refString, ref, refdescr, startpos);
    //fields: input
    samRecord_t oracleUp;
    samRecord_t oracleDown;
    vector<samRecord_t> mapped;
    ifstream oracleFile(opt.oracleSam.c_str());
    ifstream queryFile(opt.querySam.c_str());
    string oracleLine = "";
    string queryLine = "";
    streampos oraclePos = oracleFile.tellg();
    if(oracleFile.good())
        getlijn2(oracleFile, oracleLine);
    if(queryFile.good())
        getlijn2(queryFile,queryLine);
    //fields: output
    long readcnt = 0;
    size_t rangeCount = opt.correctRange.size();
    size_t qualValCount = opt.qualityValues.size();
    oracleLine_t ** results = new oracleLine_t * [rangeCount];
    for(size_t i = 0; i < rangeCount; i++){
        results[i] = new oracleLine_t[qualValCount+2];
    }
    
    //init results
    for(size_t i=0; i<rangeCount; i++)
        for(size_t j=0; j<qualValCount+2; j++)
            results[i][j].init();
    
    //iterate over both SORTED sam files simultaneously
    cerr << "summary: checking accuracy of SAM file, using oracle SAM file containing the original simulated positions" << endl;
    //print header of both files
    cerr << "header of oracle SAM file: " << endl;
    while(oracleFile.good() && !oracleLine.empty() && oracleLine[0]=='@'){
        cerr << oracleLine << endl;
        oraclePos = oracleFile.tellg();
        getlijn2(oracleFile,oracleLine);
    }
    cerr << "header of query SAM file: " << endl;
    while(queryFile.good() && !queryLine.empty() && queryLine[0]=='@'){
        cerr << queryLine << endl;
        getlijn2(queryFile,queryLine);
    }
    cerr << endl;
    cerr << "the ranges for which an alignment is considered valid are: " << printVector(opt.correctRange) << endl;
    cerr << "extra lines containing only alignments with minimal quality values of: " << printVector(opt.qualityValues) << endl;    
    cerr << "checking accuracy of alignments: ..." << endl;
    //if at least one record can be found
    bool paired = false;
    if(opt.print)
        cerr << "summary of reads: " << endl;
    if(!oracleLine.empty() && oracleLine[0]!='@' && !queryLine.empty() && queryLine[0]!='@'){
        //check if it are paired reads
        samRecord_t tempQuery;
        samRecord_t tempOracle;
        readRecord(oracleLine, tempOracle,opt.otag);
        readRecord(queryLine, tempQuery,opt.qtag);
        paired = tempOracle.flag.test(0);
        oracleFile.seekg(oraclePos);
        if(paired){
            while(oracleFile.good()){//next reads
                getlijn2(oracleFile, oracleLine);
                readRecord(oracleLine, tempOracle,opt.otag);
                if(!tempOracle.qname.empty()){
                    tempOracle.flag.test(6) ? oracleUp = tempOracle : oracleDown = tempOracle;
                    getlijn2(oracleFile, oracleLine);
                    readRecord(oracleLine, tempOracle,opt.otag);
                    tempOracle.flag.test(6) ? oracleUp = tempOracle : oracleDown = tempOracle;
                    mapped.clear();
                    //iterate over the SAM file. To work query should always be more advanced than oracle
                    while(queryFile.good() && !queryLine.empty() && tempQuery.qname.compare(oracleUp.qname)==0){
                        //fill mappedUp
                        if(tempQuery.edit<0 && tempQuery.rname.compare("*")!=0){
                            tempQuery.edit = calculatedDPScore(ref, refdescr, startpos, queryLine);
                        }
                        mapped.push_back(tempQuery);
                        getlijn2(queryFile, queryLine);
                        readRecord(queryLine, tempQuery, opt.qtag); 
                    }
                    if(!queryLine.empty() && tempQuery.qname.compare(oracleUp.qname)==0){//add last query of file
                        if(tempQuery.edit<0 && tempQuery.rname.compare("*")!=0){
                            tempQuery.edit = calculatedDPScore(ref, refdescr, startpos, queryLine);
                        }
                        mapped.push_back(tempQuery);
                    }
                    //process this read's alignments (paired and unpaired)
                    for(size_t i=0; i < mapped.size(); i++){
                        samRecord_t& query = mapped[i];
                        if(!query.flag.test(2)){
                            for(size_t j=0; j < rangeCount; j++){
                                int match = pairedCompareSam(query, query.flag.test(6) ? oracleUp : oracleDown, opt.correctRange[j]);
                                //all
                                results[j][0].addAlignment(match);
                                //unique
                                if(mapped.size()==2) results[j][1].addAlignment(match);
                                //quality values
                                for(size_t k=0; k < qualValCount; k++)
                                    if(query.mapq >= opt.qualityValues[k])
                                        results[j][2+k].addAlignment(match);
                            }
                        }
                    }
                    //sumarize for read
                    for(size_t j=0; j < rangeCount; j++){
                        if(opt.print)
                            cerr << oracleUp.qname << "\t" << (results[j][0].tempReadResult==0 ? 0 : 1) << endl;
                        results[j][0].finishRead();
                        results[j][1].finishRead();
                        for(size_t k=0; k < qualValCount; k++)
                            results[j][2+k].finishRead();
                    }
                    readcnt++;//update readcnt
                }
            }
        }
        else{
            while(oracleFile.good()){//next reads
                getlijn2(oracleFile, oracleLine);
                readRecord(oracleLine, oracleUp,opt.otag);
                if(!oracleUp.qname.empty()){
                    mapped.clear();
                    //iterate over the SAM file. To work query should always be more advanced than oracle
                    while(queryFile.good() && !queryLine.empty() && tempQuery.qname.compare(oracleUp.qname)==0){
                        //fill mappedUp
                        if(tempQuery.edit<0 && tempQuery.rname.compare("*")!=0){
                            tempQuery.edit = calculatedDPScore(ref, refdescr, startpos, queryLine);
                        }
                        mapped.push_back(tempQuery);
                        getlijn2(queryFile, queryLine);
                        readRecord(queryLine, tempQuery, opt.qtag); 
                    }
                    if(!queryLine.empty() && tempQuery.qname.compare(oracleUp.qname)==0){//add last query of file
                        if(tempQuery.edit<0 && tempQuery.rname.compare("*")!=0){
                            tempQuery.edit = calculatedDPScore(ref, refdescr, startpos, queryLine);
                        }
                        mapped.push_back(tempQuery);
                    }
                    //process this read's alignments
                    for(size_t i=0; i < mapped.size(); i++){
                        if(!mapped[i].flag.test(2)){
                            for(size_t j=0; j < rangeCount; j++){
                                int match = unpairedCompareSam(mapped[i], oracleUp, opt.correctRange[j]);
                                //all
                                results[j][0].addAlignment(match);
                                //unique
                                if(mapped.size()==1) results[j][1].addAlignment(match);
                                //quality values
                                for(size_t k=0; k < qualValCount; k++)
                                    if(mapped[i].mapq >= opt.qualityValues[k])
                                        results[j][2+k].addAlignment(match);
                            }
                        }
                    }
                    //sumarize for read
                    for(size_t j=0; j < rangeCount; j++){
                        if(opt.print)
                            cerr << oracleUp.qname << "\t" << (results[j][0].tempReadResult==0 ? 0 : 1) << endl;
                        results[j][0].finishRead();
                        results[j][1].finishRead();
                        for(size_t k=0; k < qualValCount; k++)
                            results[j][2+k].finishRead();
                    }
                    readcnt++;//update readcnt
                }
            }
        }
    }
    cerr << "checking accuracy of alignments: done" << endl;
    queryFile.close();
    oracleFile.close();
    cerr << endl;
    cerr << "printing results" << endl;
    cerr << "Warning, these results only make sense if both the oracle and query file are sorted according to \
    query name and if the edit distance field is filled correctly." << endl;
    if(!opt.outputFile.empty()){
        ofstream outfile ( opt.outputFile.c_str() );
        outfile << "##ALFALFA check alignment accuracy using an oracle file for simulated reads" << endl;
        outfile << "##" << endl;
        outfile << "##Warning: these results only make sense if both the oracle and query file are sorted according to \
        query name and if the edit distance field is filled correctly." << endl;
        outfile << "##" << endl;
        outfile << "##SAM file containing alignments: " << opt.querySam << endl;
        outfile << "##SAM file containing simulated origins: " << opt.oracleSam << endl;
        outfile << "##" << endl;
        outfile << "##An alignment is considered correct if it falls within a certain range from the simulated origin" << endl;
        outfile << "##The ranges that are considered are: " << printVector(opt.correctRange) << endl;
        outfile << "##" << endl;
        outfile << "##Results will be presented for all alignments in the file, unique alignments and for the \
        fraction of reads for which the mapping quality is at least: " << printVector(opt.qualityValues)  << endl;
        outfile << "##" << endl;
        outfile << "##The first 6 columns represent number and percentage of alignments that are correct, ok \
        (not close to origin/not paired correctly, but with a smaller or equal edit distance) and incorrect alignments." << endl;
        outfile << "##The second 6 columns represent these numbers and percentages for the reads. A read is considered correctly mapped\
        if at least one correct alignment has been found (and likewise for ok reads)." << endl;
        outfile << "##" << endl;
        outfile << "##" << endl;
        outfile << "##Results" << endl;
        outfile << "##" << endl;
        outfile << "##Number of reads: " << readcnt << endl;
        outfile << "##Reads are paired?: " << paired << endl;
        outfile << "##" << endl;
        for(size_t i=0; i < rangeCount; i++){
            outfile << "##Results for allowed range of " << opt.correctRange[i] << endl;
            outfile << "#\t\t\t alignments\t\t\t\t\t\t reads" << endl;
            outfile << "#resultType\t num_correct\t percent_correct\t num_ok\t percent_ok\t /"
                    "num_incorrect\t percent_incorrect\t num_correct\t percent_correct\t /"
                    "num_ok\t percent_ok\t num_incorrect\t percent_incorrect" << endl;
            outfile << "all\t " << results[i][0].printLine(readcnt);
            outfile << "unique aligned\t " << results[i][1].printLine(readcnt);
            for(size_t j=0; j< qualValCount; j++){
                outfile << "min mapq of " << opt.qualityValues[j] << "\t " << results[i][2+j].printLine(readcnt);
            }
            outfile << "##" << endl;
            outfile << "##" << endl;
        }
        outfile.close();
    }
    else{
        cout << "##ALFALFA check alignment accuracy using an oracle file for simulated reads" << endl;
        cout << "##" << endl;
        cout << "##Warning: these results only make sense if both the oracle and query file are sorted according to \
        query name and if the edit distance field is filled correctly." << endl;
        cout << "##" << endl;
        cout << "##SAM file containing alignments: " << opt.querySam << endl;
        cout << "##SAM file containing simulated origins: " << opt.oracleSam << endl;
        cout << "##" << endl;
        cout << "##An alignment is considered correct if it falls within a certain range from the simulated origin" << endl;
        cout << "##The ranges that are considered are: " << printVector(opt.correctRange) << endl;
        cout << "##" << endl;
        cout << "##Results will be presented for all alignments in the file, unique alignments and for the \
        fraction of reads for which the mapping quality is at least: " << printVector(opt.qualityValues)  << endl;
        cout << "##" << endl;
        cout << "##The first 6 columns represent number and percentage of alignments that are correct, ok \
        (not close to origin/not paired correctly, but with a smaller or equal edit distance) and incorrect alignments." << endl;
        cout << "##The second 6 columns represent these numbers and percentages for the reads. A read is considered correctly mapped\
        if at least one correct alignment has been found (and likewise for ok reads)." << endl;
        cout << "##" << endl;
        cout << "##" << endl;
        cout << "##Results" << endl;
        cout << "##" << endl;
        cout << "##Number of reads: " << readcnt << endl;
        cout << "##Reads are paired?: " << paired << endl;
        cout << "##" << endl;
        for(size_t i=0; i < rangeCount; i++){
            cout << "##Results for allowed range of " << opt.correctRange[i] << endl;
            cout << "#\t\t\t alignments\t\t\t\t\t\t reads" << endl;
            cout << "#resultType\t num_correct\t percent_correct\t num_ok\t percent_ok\t /"
                    "num_incorrect\t percent_incorrect\t num_correct\t percent_correct\t /"
                    "num_ok\t percent_ok\t num_incorrect\t percent_incorrect" << endl;
            cout << "all\t " << results[i][0].printLine(readcnt);
            cout << "unique aligned\t " << results[i][1].printLine(readcnt);
            for(size_t j=0; j< qualValCount; j++){
                cout << "min mapq of " << opt.qualityValues[j] << "\t " << results[i][2+j].printLine(readcnt);
            }
            cout << "##" << endl;
            cout << "##" << endl;
        }        
    }
    for(size_t i = 0; i < rangeCount; i++)
        delete[] results[ i ];
    delete[] results;
}

void checkWgsim(samCheckOptions_t & opt){
    //initialize fields and data for calculating edit distance
    string ref;
    vector<string> refdescr;
    vector<long> startpos;
    load_fasta(opt.refString, ref, refdescr, startpos);
    //should provide sorted query file to count query and check missing queries
    if(opt.print)
        cerr << "Warning current version does not take into account unmapped reads." << endl;
    //fields: input
    ifstream queryFile(opt.querySam.c_str());
    string queryLine = "";
    if(queryFile.good())
        getlijn2(queryFile,queryLine);
    //fields: output
    long readcnt = 0;
    size_t rangeCount = opt.correctRange.size();
    size_t qualValCount = opt.qualityValues.size();
    oracleLine_t ** results = new oracleLine_t * [rangeCount];
    for(size_t i = 0; i < rangeCount; i++){
        results[i] = new oracleLine_t[qualValCount+2];
    }

    //init results
    for(size_t i=0; i<rangeCount; i++)
        for(size_t j=0; j<qualValCount+2; j++)
            results[i][j].init();
    
    //iterate over both SORTED sam files simultaneously
    cerr << "summary: checking accuracy of SAM file for reads produced by wgsim" << endl;
    //print header of SAM file
    cerr << "header of SAM file: " << endl;
    while(queryFile.good() && !queryLine.empty() && queryLine[0]=='@'){
        cerr << queryLine << endl;
        getlijn2(queryFile,queryLine);
    }
    cerr << endl;
    cerr << "the ranges for which an alignment is considered valid are: " << printVector(opt.correctRange) << endl;
    cerr << "extra lines containing only alignments with minimal quality values of: " << printVector(opt.qualityValues) << endl;    
    cerr << "checking accuracy of alignments: ..." << endl;
    //if at least one record can be found
    string prevQname = "";
    string qname = "";
    int curAln=0;
    char delimeter = '\t';
    if(opt.print)
        cerr << "reads that failed to align correctly: " << endl;
    if(queryFile.good() && !queryLine.empty()){
        do{
            //read in all fields for this line
            int tabPos = 0;
            //qname
            qname = nextField(queryLine, delimeter, tabPos);
            if(qname.length() > 2 && qname[qname.length()-2]=='/')
                qname.erase(qname.length()-2);
            int namePos=qname.size()-1;
            previousField(qname, '_', namePos);//counter
            int revEdit = atoi(previousField(qname, ':', namePos).c_str());
            revEdit += atoi(previousField(qname, ':', namePos).c_str());
            revEdit += atoi(previousField(qname, '_', namePos).c_str());
            int fwdEdit = atoi(previousField(qname, ':', namePos).c_str());
            fwdEdit += atoi(previousField(qname, ':', namePos).c_str());
            fwdEdit += atoi(previousField(qname, '_', namePos).c_str());
            long revPos = atoi(previousField(qname, '_', namePos).c_str());
            long fwdPos = atoi(previousField(qname, '_', namePos).c_str());
            string realchrom = qname.substr(0, namePos+1);
            //flag
            bitset<11> flag = bitset<11>((ulong) atoi(nextField(queryLine, delimeter, tabPos).c_str()));
            //rname
            string chrom = nextField(queryLine, delimeter, tabPos);
            //pos
            long leftpos = atoi(nextField(queryLine, delimeter, tabPos).c_str());
            long rightpos = leftpos;
            //mapq
            int mapq = atoi(nextField(queryLine, delimeter, tabPos).c_str());
            //CIGAR
            string cigar = nextField(queryLine, delimeter, tabPos);
            //determine left and right pos using cigar
            int cigarPos = 0;
            int cigarlength = cigar.size();
            while(cigarPos < cigarlength){
                int beginPos = cigarPos;
                while(cigar[cigarPos] >= 48 && cigar[cigarPos] <= 57)//integer
                    cigarPos++;
                //now cigar[cigarPos] is character
                if(cigar[cigarPos] == 'M' || cigar[cigarPos] == 'N' || cigar[cigarPos] == 'D')
                    rightpos += atoi(cigar.substr(beginPos,cigarPos-beginPos).c_str());
                cigarPos++;
            }
            rightpos--;
            //correct left and right pos for clipping
            long leftpos0 = leftpos;
            long rightpos0 = rightpos;
            cigarPos = 0;
            while(cigar[cigarPos] >= 48 && cigar[cigarPos] <= 57)//integer
                cigarPos++;
                //now cigar[cigarPos] is character
            if(cigar[cigarPos] == 'S' || cigar[cigarPos] == 'H'){
                leftpos -= atoi(cigar.substr(0,cigarPos).c_str());
                rightpos0 += atoi(cigar.substr(0,cigarPos).c_str());
            }
            if(cigar[cigar.size()-1] == 'S' || cigar[cigar.size()-1] == 'H'){
                cigarPos = cigar.size()-2;
                while(cigarPos >= 0 && cigar[cigarPos] >= 48 && cigar[cigarPos] <= 57)//integer
                    cigarPos--;
                rightpos += atoi(cigar.substr(cigarPos+1,cigar.size()-1-cigarPos).c_str());
                leftpos0 -= atoi(cigar.substr(cigarPos+1,cigar.size()-1-cigarPos).c_str());
            }
            //rnext
            nextField(queryLine, delimeter, tabPos);
            //pnext
            nextField(queryLine, delimeter, tabPos);
            //edit distance NM tag
            int edit = 0;
            //skip fields to make sure 
            tabPos = queryLine.find(delimeter,tabPos)+1;
            tabPos = queryLine.find(delimeter,tabPos)+1;
            tabPos = queryLine.find(delimeter,tabPos)+1;
            //now passed qual
            string editTag = opt.qtag+":i:";
            int NMpos = queryLine.find(editTag,tabPos);
            if(NMpos == -1 && !flag.test(2)){
                edit = calculatedDPScore(ref, refdescr, startpos, queryLine);
            }
            else if(NMpos>=tabPos){
                tabPos = queryLine.find(delimeter,NMpos);
                if(tabPos == -1)
                    tabPos = queryLine.length();
                edit = atoi(queryLine.substr(NMpos+5,tabPos-NMpos-5).c_str());
            }
            //if new query, finish previous
            if(prevQname != "" && prevQname.compare(qname) != 0){
                //sumarize for read
                for(size_t j=0; j < rangeCount; j++){
                    if(curAln <= 2)
                        results[j][1].addAlignment(results[j][0].tempReadResult);
                    if(opt.print && results[j][0].tempReadResult==0)
                        cerr << prevQname << endl;
                    results[j][0].finishRead();
                    results[j][1].finishRead();
                    for(size_t k=0; k < qualValCount; k++)
                        results[j][2+k].finishRead();
                }
                curAln=0;
                readcnt++;//update readcnt
            }
            prevQname = qname;
            //add alignment
            if(!flag.test(2)){
                curAln++;
                for(size_t j=0; j < rangeCount; j++){
                    //decide if current alignment is a match
                    int match = 0;
                    if(!flag.test(4)){//fwd read
                        if(realchrom.compare(chrom)==0 && 
                                    (abs((int)(fwdPos-leftpos)) <= opt.correctRange[j] || abs((int)(fwdPos-leftpos0)) <= opt.correctRange[j])){
                            match = 2;
                        }
                        else if(edit <= fwdEdit)
                            match = 1;
                    }
                    else{//reverse read
                        if(realchrom.compare(chrom)==0 && 
                                    (abs((int)(revPos-rightpos)) <= opt.correctRange[j] || abs((int)(revPos-rightpos0)) <= opt.correctRange[j])){
                            match = 2;
                        }
                        else if(edit <= revEdit)
                            match = 1;
                    }
                    //all
                    results[j][0].addAlignment(match);
                    //quality values
                    for(size_t k=0; k < qualValCount; k++)
                        if(mapq >= opt.qualityValues[k])
                            results[j][2+k].addAlignment(match);
                }
            }
            //read next line
            getlijn2(queryFile, queryLine);
        } while(queryFile.good() && !queryLine.empty());
        //if new query, finish previous
        if(prevQname != ""){
            //sumarize for read
            for(size_t j=0; j < rangeCount; j++){
                if(opt.print && results[j][0].tempReadResult==0)
                    cerr << prevQname << endl;
                results[j][0].finishRead();
                results[j][1].finishRead();
                for(size_t k=0; k < qualValCount; k++)
                    results[j][2+k].finishRead();
            }
            readcnt++;//update readcnt
        }
    }
    cerr << "checking accuracy of alignments: done" << endl;
    queryFile.close();
    cerr << endl;
    cerr << "printing results" << endl;
    cerr << "Warning, these results only make sense if the SAM file is sorted according to \
    query name and if the edit distance field is filled correctly." << endl;
    if(!opt.outputFile.empty()){
        ofstream outfile ( opt.outputFile.c_str() );
        outfile << "##ALFALFA check alignment accuracy using an oracle file for simulated reads" << endl;
        outfile << "##" << endl;
        outfile << "##Warning, these results only make sense if the SAM file is sorted according to \
        query name and if the edit distance field is filled correctly." << endl;
        outfile << "##" << endl;
        outfile << "##SAM file containing alignments: " << opt.querySam << endl;
        outfile << "##" << endl;
        outfile << "##An alignment is considered correct if it falls within a certain range from the simulated origin" << endl;
        outfile << "##The ranges that are considered are: " << printVector(opt.correctRange) << endl;
        outfile << "##" << endl;
        outfile << "##Results will be presented for all alignments in the file, unique alignments and for the \
        fraction of reads for which the mapping quality is at least: " << printVector(opt.qualityValues)  << endl;
        outfile << "##" << endl;
        outfile << "##The first 6 columns represent number and percentage of alignments that are correct, ok \
        (not close to origin/not paired correctly, but with a smaller or equal edit distance) and incorrect alignments." << endl;
        outfile << "##The second 6 columns represent these numbers and percentages for the reads. A read is considered correctly mapped\
        if at least one correct alignment has been found (and likewise for ok reads)." << endl;
        outfile << "##" << endl;
        outfile << "##" << endl;
        outfile << "##Results" << endl;
        outfile << "##" << endl;
        outfile << "##Number of reads: " << readcnt << endl;
        outfile << "##Reads are paired?: " << "true" << endl;
        outfile << "##" << endl;
        for(size_t i=0; i < rangeCount; i++){
            outfile << "##Results for allowed range of " << opt.correctRange[i] << endl;
            outfile << "#\t\t\t alignments\t\t\t\t\t\t reads" << endl;
            outfile << "#resultType\t num_correct\t percent_correct\t num_ok\t percent_ok\t /"
                    "num_incorrect\t percent_incorrect\t num_correct\t percent_correct\t /"
                    "num_ok\t percent_ok\t num_incorrect\t percent_incorrect" << endl;
            outfile << "all\t " << results[i][0].printLine(readcnt);
            outfile << "unique aligned\t " << results[i][1].printLine(readcnt);
            for(size_t j=0; j< qualValCount; j++){
                outfile << "min mapq of " << opt.qualityValues[j] << "\t " << results[i][2+j].printLine(readcnt);
            }
            outfile << "##" << endl;
            outfile << "##" << endl;
        }
        outfile.close();
    }
    else{
        cout << "##ALFALFA check alignment accuracy using an oracle file for simulated reads" << endl;
        cout << "##" << endl;
        cout << "##Warning, these results only make sense if the SAM file is sorted according to \
        query name and if the edit distance field is filled correctly." << endl;
        cout << "##" << endl;
        cout << "##SAM file containing alignments: " << opt.querySam << endl;
        cout << "##" << endl;
        cout << "##An alignment is considered correct if it falls within a certain range from the simulated origin" << endl;
        cout << "##The ranges that are considered are: " << printVector(opt.correctRange) << endl;
        cout << "##" << endl;
        cout << "##Results will be presented for all alignments in the file, unique alignments and for the \
        fraction of reads for which the mapping quality is at least: " << printVector(opt.qualityValues)  << endl;
        cout << "##" << endl;
        cout << "##The first 6 columns represent number and percentage of alignments that are correct, ok \
        (not close to origin/not paired correctly, but with a smaller or equal edit distance) and incorrect alignments." << endl;
        cout << "##The second 6 columns represent these numbers and percentages for the reads. A read is considered correctly mapped\
        if at least one correct alignment has been found (and likewise for ok reads)." << endl;
        cout << "##" << endl;
        cout << "##" << endl;
        cout << "##Results" << endl;
        cout << "##" << endl;
        cout << "##Number of reads: " << readcnt << endl;
        cout << "##Reads are paired?: " << "true" << endl;
        cout << "##" << endl;
        for(size_t i=0; i < rangeCount; i++){
            cout << "##Results for allowed range of " << opt.correctRange[i] << endl;
            cout << "#\t\t\t alignments\t\t\t\t\t\t reads" << endl;
            cout << "#resultType\t num_correct\t percent_correct\t num_ok\t percent_ok\t /"
                    "num_incorrect\t percent_incorrect\t num_correct\t percent_correct\t /"
                    "num_ok\t percent_ok\t num_incorrect\t percent_incorrect" << endl;
            cout << "all\t " << results[i][0].printLine(readcnt);
            cout << "unique aligned\t " << results[i][1].printLine(readcnt);
            for(size_t j=0; j< qualValCount; j++){
                cout << "min mapq of " << opt.qualityValues[j] << "\t " << results[i][2+j].printLine(readcnt);
            }
            cout << "##" << endl;
            cout << "##" << endl;
        }        
    }
    for(size_t i = 0; i < rangeCount; i++)
        delete[] results[ i ];
    delete[] results;
}

void checkSummary(samCheckOptions_t & opt){
    //fields: input
    ifstream queryFile(opt.querySam.c_str());
    string queryLine = "";
    streampos queryPos = queryFile.tellg();
    long readCount = 0;
    if(queryFile.good())
        getlijn2(queryFile,queryLine);//first line
    //fields: output
    int qualValCount = opt.qualityValues.size();
    summaryLine_t * results = new summaryLine_t[qualValCount+2];
    cerr << "summary: summary of alignments found (including pairing info and mapped/mapping quality" << endl;
    //print header of both files
    cerr << "header of SAM file: " << endl;
    while(queryFile.good() && !queryLine.empty() && queryLine[0]=='@'){
        cerr << queryLine << endl;
        queryPos = queryFile.tellg();
        getlijn2(queryFile,queryLine);
    }
    cerr << endl;
    cerr << "extra lines containing only alignments with minimal quality values of: " << printVector(opt.qualityValues) << endl;    
    cerr << "compiling summary: ..." << endl;
    if(!queryLine.empty() && queryLine[0]!='@'){
        string qnamePrev = "";
        int flagPrev = 0;
        string rnextPrev = "";
        string qname = "";
        int flag = 0;
        int mapq = 0;
        string rnext = "";
        int tabPos = 0;
        char delimeter = '\t';
        do{
            //set prev to current
            qnamePrev = qname;
            flagPrev = flag;
            rnextPrev = rnext;
            //translate current line to information
            tabPos = 0;
            //qname
            qname = nextField(queryLine, delimeter, tabPos);
            if(qname.length() > 2 && qname[qname.length()-2]=='/')
              qname.erase(qname.length()-2);
            //flag
            flag = atoi(nextField(queryLine, delimeter, tabPos).c_str());
            //skip fields rname, pos
            tabPos = queryLine.find(delimeter,tabPos)+1;//behind rname
            tabPos = queryLine.find(delimeter,tabPos)+1;//behind pos
            //mapq
            mapq = atoi(nextField(queryLine, delimeter, tabPos).c_str());
            //skip CIGAR
            tabPos = queryLine.find(delimeter,tabPos)+1;//behind cigar
            rnext = nextField(queryLine, delimeter, tabPos);
            //fill in the fields
            if(qname.compare(qnamePrev)!=0 && !qnamePrev.empty()){
                readCount++;
                if(results[0].tempAlnCount==1 || (opt.paired && results[0].tempAlnCount==2))
                    results[1].addAlignment(flagPrev, rnextPrev.compare("*")!=0);
                results[0].finishRead();
                results[1].finishRead();
                for(int i=0; i < qualValCount; i++)
                    results[2+i].finishRead();
            }
            if(flag!=4){
                results[0].addAlignment(flag, rnext.compare("*")!=0);
                for(int i=0; i < qualValCount; i++)
                    if(mapq >= opt.qualityValues[i])
                        results[2+i].addAlignment(flag, rnext.compare("*")!=0);
            }
            getlijn2(queryFile,queryLine);//next line
        }while(queryFile.good());
        if(!qnamePrev.empty()){
            readCount++;
            if(results[0].tempAlnCount==1 || (opt.paired && results[0].tempAlnCount==2))
                results[1].addAlignment(flagPrev, rnextPrev.compare("*")!=0);
            results[0].finishRead();
            results[1].finishRead();
            for(int i=0; i < qualValCount; i++)
                results[2+i].finishRead();
        }
    }
    cerr << "compiling summary: done" << endl;
    queryFile.close();
    cerr << endl;
    cerr << "printing results" << endl;
    cerr << "Warning, these results only make sense if the query file is sorted according to query name." << endl;
    long maxReads = opt.numReads;
    if(maxReads==0){
        maxReads = readCount;
        cerr << "Warning, these number of reads is taken from the result and might differ from the number of reads that were given to the aligner." << endl;
    }
    if(!opt.outputFile.empty()){
        ofstream outfile ( opt.outputFile.c_str() );
        outfile << "##ALFALFA summary of alignment accuracy" << endl;
        outfile << "##" << endl;
        outfile << "##Warning: these results only make sense if the query file is sorted according to query name." << endl;
        outfile << "##" << endl;
        outfile << "##SAM file containing alignments: " << opt.querySam << endl;
        outfile << "##Results will be presented for all alignments in the file, unique alignments and for the \
        fraction of reads for which the mapping quality is at least: " << printVector(opt.qualityValues)  << endl;
        outfile << "##For every category, the fraction of alignments/reads mapped (and paired correctly) will be shown." << endl;
        outfile << "##Furthermore, a list is given of number of alignments per read." << endl;
        outfile << "##" << endl;
        outfile << "##" << endl;
        outfile << "##Results" << endl;
        outfile << "##" << endl;
        outfile << "##Number of reads: " << maxReads << endl;
        outfile << "##Reads are paired?: " << opt.paired << endl;
        outfile << "##" << endl;
        outfile << "##All alignments:" << endl;
        outfile << results[0].printLine(maxReads, opt.paired);
        outfile << "##" << endl;
        outfile << "##Unique alignments:" << endl;
        outfile << results[1].printLine(maxReads, opt.paired);
        for(int j=0; j< qualValCount; j++){
            outfile << "##Alignments with min Q-value " << opt.qualityValues[j] << endl;
            outfile << results[2+j].printLine(maxReads, opt.paired);
        }
        outfile.close();
    }
    else{
        cout << "##ALFALFA summary of alignment accuracy" << endl;
        cout << "##" << endl;
        cout << "##Warning: these results only make sense if the query file is sorted according to query name." << endl;
        cout << "##" << endl;
        cout << "##SAM file containing alignments: " << opt.querySam << endl;
        cout << "##Results will be presented for all alignments in the file, unique alignments and for the \
        fraction of reads for which the mapping quality is at least: " << printVector(opt.qualityValues)  << endl;
        cout << "##For every category, the fraction of alignments/reads mapped (and paired correctly) will be shown." << endl;
        cout << "##Furthermore, a list is given of number of alignments per read." << endl;
        cout << "##" << endl;
        cout << "##" << endl;
        cout << "##Results" << endl;
        cout << "##" << endl;
        cout << "##Number of reads: " << maxReads << endl;
        cout << "##Reads are paired?: " << opt.paired << endl;
        cout << "##" << endl;
        cout << "##All alignments:" << endl;
        cout << results[0].printLine(maxReads, opt.paired);
        cout << "##" << endl;
        cout << "##Unique alignments:" << endl;
        cout << results[1].printLine(maxReads, opt.paired);
        for(int j=0; j< qualValCount; j++){
            cout << "##Alignments with min Q-value " << opt.qualityValues[j] << endl;
            cout << results[2+j].printLine(maxReads, opt.paired);
        }
    }
    delete [] results;
}

