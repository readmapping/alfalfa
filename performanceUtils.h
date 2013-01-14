/* 
 * File:   performanceUtils.h
 * Author: mvyvermn
 *
 * Created on 27 september 2012, 14:17
 */

#ifndef PERFORMANCEUTILS_H
#define	PERFORMANCEUTILS_H

#include <getopt.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define SET_BIT(var,pos) ((var) | (1<<(pos)))

//should be located somewhere else
template<class T>
string printVector(vector<T> vec){
    stringstream ss;
    if(!vec.empty()){
        ss << vec[0];
        for(int i=1; i < vec.size(); i++)
            ss << "," << vec[i];
    }
    return ss.str();
}

enum checkCommand_t {ORACLE, SUMMARY, COMPARE};

struct samCheckOptions_t {
    samCheckOptions_t(){ initOptions(); }
    void initOptions(){
        subcommand = SUMMARY;
        outputFile = "";
        oracleSam = "";
        paired = false;
        querySam = "";
        qtag = "NM";
        otag = "NM";
        numReads=0;
    }
    checkCommand_t subcommand;
    string outputFile;
    string oracleSam;
    bool paired;
    string querySam;
    vector<string> compareFiles;
    vector<int> qualityValues;
    vector<int> correctRange;
    string qtag;
    string otag;
    long numReads;
    //other options to change or reduce the output
};

struct samRecord_t {
    samRecord_t(): qname(""), flag(0), rname(""), pos(0), mapq(0), rnext(""), pnext(0), edit(0){}
    samRecord_t(const samRecord_t& o): qname(o.qname), rname(o.rname), pos(o.pos), 
    mapq(o.mapq), rnext(o.rnext), pnext(o.pnext), edit(o.edit), flag(o.flag.to_ulong()) {}
    string qname;
    bitset<11> flag;
    string rname;
    long pos;
    int mapq;
    string rnext;
    long pnext;
    int edit;
};

//o: output filename
static const char * short_optionsC = "o:S:I:2:3:4:Q:r:hC:";

//Numbers for options without short option
enum checkOptionNumber {
    ARG_NO_OPT = 255,          //not found
    ARG_QTAG,      //--qtag
    ARG_OTAG,    //--otag
    ARG_PAIRED,     //--paired
    ARG_NUM_READS   //--num-reads
};

static struct option long_optionsC[] = {
    {(char*)"output",           required_argument, 0,            'o'},
    {(char*)"oracle",           required_argument, 0,            'S'},
    {(char*)"input-sam",        required_argument, 0,            'I'},
    {(char*)"compare-sam",        required_argument, 0,          'C'},
    {(char*)"quality-values",   required_argument, 0,            'Q'},
    {(char*)"ranges",           required_argument, 0,            'r'},
    {(char*)"help",             no_argument, 0,                  'h'},
    {(char*)"qtag",             required_argument, 0,            ARG_QTAG},
    {(char*)"otag",             required_argument, 0,            ARG_OTAG},
    {(char*)"paired",           no_argument, 0,                  ARG_PAIRED},
    {(char*)"num-reads",        required_argument, 0,            ARG_NUM_READS},
    {(char*)0, 0, 0, 0} // terminator
};

static void usageCheck(const string prog) {
  cerr << "Usage: " << prog << " check SUBCOMMAND [options]" << endl;
  cerr << "SUBCOMMAND should be one of the following: " << endl;
  cerr << "oracle                      compare the output SAM with an oracle file containing original simulated mapping positions" << endl;
  cerr << "summary                     give a summary of number of mapped reads and alignments according to mapping quality" << endl;
  cerr << "compare                     compare the result of up to 32 SAM results, showing the number of alignments/reads mapped that are found by a combination of the SAM files" << endl;
  cerr << endl;
  cerr << "COMMON OPTIONS" << endl;
  cerr << "-h/--help                  print this statement" << endl;
  cerr << "-o/--output                the filename to write the output to [std-out]" << endl;
  cerr << endl;
  cerr << "ORACLE OPTIONS" << endl;
  cerr << "-S/--oracle (string)       reference SAM containing the simulated read coordinates" << endl;
  cerr << "-I/--input-sam             the SAM file to check" << endl;
  cerr << "-Q/--quality-values        comma separated list of quality values, for which the statistics will be displayed. [0 included always]" << endl;
  cerr << "-r/--ranges                comma separated list of ranges to the simulated origin, for which the statistics will be displayed. [50]" << endl;
  cerr << "                           for example: range 50 means that an alignment found within 50bp of the simulated origin is considered correct." << endl;
  cerr << "--qtag (string)            the tag that stores the edit distance in the --input-sam file [NM]" << endl;
  cerr << "--otag (string)            the tag that stores the edit distance in the --oracle file [NM]" << endl;
  cerr << endl;
  cerr << "SUMMARY OPTIONS " << endl;
  cerr << "-I/--input-sam             the SAM file to check" << endl;
  cerr << "--paired                   the reads were aligned in a paired fashion [false]" << endl;
  cerr << "--num-reads                the number of reads that were (tried to be) aligned" << endl;
  cerr << "-Q/--quality-values        comma separated list of quality values, for which the statistics will be displayed. [0 included always]" << endl;
  cerr << endl;
  cerr << "COMPARE OPTIONS " << endl;
  cerr << "--paired                   the reads were aligned in a paired fashion" << endl;
  cerr << "--num-reads                the number of reads that were (tried to be) aligned" << endl;
  cerr << "-C/--compare-sam           comma-separated list of SAM files to compare (at least two). All SAM files should be ordered according to read name" << endl;
  cerr << "-r/--ranges                comma separated list of ranges to the simulated origin, for which the statistics will be displayed. [50]" << endl;
  cerr << "                           for example: range 50 means that an alignment found within 50bp of the simulated origin is considered correct." << endl;
  exit(1);
}

static void processCommaSepListInt(string list, vector<int> & options, string name, int lowerBound, int upperBound){
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

static void processCommaSepListString(string list, vector<string> & options){
    int beginPos = 0;
    int commaPos = list.find(',',beginPos);
    while(commaPos>=0){
        options.push_back(list.substr(beginPos,commaPos-beginPos));
        beginPos = commaPos+1;
        commaPos = list.find(',',beginPos);
    }
}

static void processParameters(int argc, char* argv[], samCheckOptions_t& opt, const string program){
    //parse the subcommand
    string subcommand = argv[1];
    if(strcmp(argv[1], "oracle") == 0) opt.subcommand = ORACLE;
    else if(strcmp(argv[1], "summary") == 0) opt.subcommand = SUMMARY;
    else if(strcmp(argv[1], "compare") == 0) opt.subcommand = COMPARE;
    else{
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        usageCheck(program);
        exit(1);
    }
    cerr << "COMMAND: check " << argv[1] << endl;
    int option_index = 0;
    int c;
    while ((c = getopt_long(
    argc-1, argv+1,
    short_optionsC, long_optionsC, &option_index)) != -1) {//memType cannot be chosen
        switch (c) {//"o:S:U:D:1:2:3:4:Q:r:h"//report parameters that will be ignored, or failed
            case 'o': opt.outputFile = optarg; break;
            case 'S': if(opt.subcommand == ORACLE) opt.oracleSam = optarg;
            else fprintf(stderr, "oracle SAM file not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case ARG_PAIRED: if(opt.subcommand != ORACLE) opt.paired = true;
            else fprintf(stderr, "paired/ not paired not necessairy for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case ARG_NUM_READS: if(opt.subcommand != ORACLE) opt.numReads = atoi(optarg);
            else fprintf(stderr, "number of reads not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case 'I': if(opt.subcommand != COMPARE) opt.querySam = optarg;
            else fprintf(stderr, "separate SAM file not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case 'C': if(opt.subcommand == COMPARE) processCommaSepListString(optarg, opt.compareFiles);
            else fprintf(stderr, "more SAM files not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case 'h': usageCheck(program); break;
            case 'Q': if(opt.subcommand != COMPARE) processCommaSepListInt(optarg, opt.qualityValues, "quality-values", 0, 255);
            else fprintf(stderr, "quality value cut-offs not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str());
            break;
            case 'r': if(opt.subcommand != SUMMARY) processCommaSepListInt(optarg, opt.correctRange, "ranges", 0, 1000000);
            else fprintf(stderr, "quality value cut-offs not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case ARG_QTAG: if(opt.subcommand == ORACLE) opt.qtag = optarg;
            else fprintf(stderr, "tag for edit distance not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str());
            break;
            case ARG_OTAG: if(opt.subcommand == ORACLE) opt.otag = optarg;
            else fprintf(stderr, "tag for edit distance not needed for '%s' mode. This parameter will be ignored.\n", subcommand.c_str()); 
            break;
            case -1: /* Done with options. */
            break;
            case 0: if (long_options[option_index].flag != 0) break;
            default: usageCheck(program);
            throw 1;
        }
    }
    if(opt.subcommand != COMPARE && opt.querySam.empty()){
        fprintf(stderr, "The query SAM file has not been specified using the -I parameter: terminate.\n"); exit(1);
    }
    if(opt.subcommand == COMPARE && (opt.compareFiles.empty() || opt.compareFiles.size() < 2)){
        fprintf(stderr, "Not enough SAM files specified to compare.\n"); exit(1);
    }
    if(opt.subcommand == ORACLE && opt.oracleSam.empty()){
        fprintf(stderr, "The oracle SAM file has not been specified for the oracle mode: terminate.\n"); exit(1);
    }
    else{
        if(opt.subcommand != ORACLE && opt.numReads==0){
            fprintf(stderr, "Number of reads was not specified.\n"); exit(1);
        }
    }
    if(opt.subcommand != SUMMARY && opt.correctRange.empty())
        opt.correctRange.push_back(50);
}

static string nextField(string & line, char delimeter, int& beginPos){
    int tabPos = line.find(delimeter,beginPos);
    string substring = line.substr(beginPos, tabPos-beginPos);
    beginPos = tabPos+1;
    return substring;
}

//TODO: temp method copied from fasta.cpp, should be placed elsewhere
void getlijn2(ifstream & input, string & lijn){
    //write windows version of getline
    getline(input, lijn, '\n'); // Load one line at a time.
    if(lijn.length() > 0 && lijn[lijn.length()-1]=='\r')
        lijn.erase(--lijn.end());
}

static bool readRecord(string & line, samRecord_t & record, string tag){
    int tabPos = 0;
    char delimeter = '\t';
    //qname
    record.qname = nextField(line, delimeter, tabPos);
    //flag
    record.flag = bitset<11>((ulong) atoi(nextField(line, delimeter, tabPos).c_str()));
    //rname
    record.rname = nextField(line, delimeter, tabPos);
    //pos
    record.pos = atoi(nextField(line, delimeter, tabPos).c_str());
    //mapq
    record.mapq = atoi(nextField(line, delimeter, tabPos).c_str());
    //skip CIGAR
    nextField(line, delimeter, tabPos);
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
        return false;
    else if(NMpos>=tabPos){
        tabPos = line.find(delimeter,NMpos);
        if(tabPos == -1)
            tabPos = line.length();
        record.edit = atoi(line.substr(NMpos+5,tabPos-NMpos-5).c_str());
    }
    return true;
}

struct oracleLine_t {
    oracleLine_t(): alignmentCnt(0),alnGood(0),alnOk(0),alnBad(0),readGood(0),readOk(0),readBad(0), tempReadResult(0){}
    long alignmentCnt;
    long alnGood;
    long alnOk;
    long alnBad;
    long readGood;
    long readOk;
    long readBad;
    int tempReadResult;
    void addAlignment(int score){
        tempReadResult = max(tempReadResult, score);
        alignmentCnt++;
        if(score == 2)
            alnGood++;
        else if(score == 1)
            alnOk++;
        else
            alnBad++;
    }
    void finishRead(){
        if(tempReadResult==2)
            readGood++;
        else if(tempReadResult ==1)
            readOk++;
        else
            readBad++;
        tempReadResult = 0;
    }
    string printLine(long readcnt){
        stringstream ss;
        ss << alnGood << "\t " << ((alnGood*100)/alignmentCnt) << "\t ";
        ss << alnOk << "\t " << ((alnOk*100)/alignmentCnt) << "\t "; 
        ss << alnBad << "\t " << ((alnBad*100)/alignmentCnt) << "\t ";
        ss << readGood << "\t " << ((readGood*100)/readcnt) << "\t ";
        ss << readOk << "\t " << ((readOk*100)/readcnt) << "\t ";
        ss << readBad << "\t " << ((readBad*100)/readcnt) << endl;
        return ss.str();
    }
};

static int unpairedCompareSam(samRecord_t& query, samRecord_t& oracle, int range){
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

static int pairedCompareSam(samRecord_t& query, samRecord_t& oracle, int range){
    //correct if within X bases of origin
    if(query.flag.test(2)){
        return 0;//unmapped
    }
    if((query.flag.to_ulong() % 128) == (oracle.flag.to_ulong() % 128) &&
            query.rname.compare(oracle.rname)==0 && 
            query.pos >= oracle.pos-range && query.pos <= oracle.pos + range && 
            query.rnext.compare(oracle.rnext) == 0 &&
            query.pnext >= oracle.pnext-range && query.pnext <= oracle.pnext+range){
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

static void checkOracle(samCheckOptions_t & opt){
    //fields: input
    samRecord_t oracleUp;
    samRecord_t oracleDown;
    vector<samRecord_t> mapped;
    ifstream oracleFile(opt.oracleSam.c_str());
    ifstream queryFile(opt.querySam.c_str());
    string oracleLine = "";
    string queryLine = "";
    streampos oraclePos = oracleFile.tellg();
    if(!oracleFile.eof())
        getlijn2(oracleFile, oracleLine);
    if(!queryFile.eof())
        getlijn2(queryFile,queryLine);
    //fields: output
    ofstream outfile ( opt.outputFile.c_str() );
    long readcnt = 0;
    int rangeCount = opt.correctRange.size();
    int qualValCount = opt.qualityValues.size();
    oracleLine_t results[rangeCount][qualValCount+2];
    //iterate over both SORTED sam files simultaneously
    cerr << "summary: checking accuracy of SAM file, using oracle SAM file containing the original simulated positions" << endl;
    //print header of both files
    cerr << "header of oracle SAM file: " << endl;
    while(!oracleFile.eof() && !oracleLine.empty() && oracleLine[0]=='@'){
        cerr << oracleLine << endl;
        oraclePos = oracleFile.tellg();
        getlijn2(oracleFile,oracleLine);
    }
    cerr << "header of query SAM file: " << endl;
    while(!queryFile.eof() && !queryLine.empty() && queryLine[0]=='@'){
        cerr << queryLine << endl;
        getlijn2(queryFile,queryLine);
    }
    cerr << endl;
    cerr << "the ranges for which an alignment is considered valid are: " << printVector(opt.correctRange) << endl;
    cerr << "extra lines containing only alignments with minimal quality values of: " << printVector(opt.qualityValues) << endl;    
    cerr << "checking accuracy of alignments: ..." << endl;
    //if at least one record can be found
    bool paired = false;
    if(!oracleLine.empty() && oracleLine[0]!='@' && !queryLine.empty() && queryLine[0]!='@'){
        //check if it are paired reads
        samRecord_t tempQuery;
        samRecord_t tempOracle;
        readRecord(oracleLine, tempOracle,opt.otag);
        readRecord(queryLine, tempQuery,opt.qtag);
        paired = tempOracle.flag.test(0);
        oracleFile.seekg(oraclePos);
        if(paired){
            while(!oracleFile.eof()){//next reads
                getlijn2(oracleFile, oracleLine);
                readRecord(oracleLine, tempOracle,opt.otag);
                tempOracle.flag.test(6) ? oracleUp = tempOracle : oracleDown = tempOracle;
                getlijn2(oracleFile, oracleLine);
                readRecord(oracleLine, tempOracle,opt.otag);
                tempOracle.flag.test(6) ? oracleUp = tempOracle : oracleDown = tempOracle;
                mapped.clear();
                //iterate over the SAM file. To work query should always be more advanced than oracle
                while(!queryFile.eof() && tempQuery.qname.compare(oracleUp.qname)==0){
                    //fill mappedUp
                    mapped.push_back(tempQuery);
                    getlijn2(queryFile, queryLine);
                    readRecord(queryLine, tempQuery, opt.qtag); 
                }
                if(tempQuery.qname.compare(oracleUp.qname)==0){//add last query of file
                    mapped.push_back(tempQuery);
                }
                //process this read's alignments (paired and unpaired)
                for(int i=0; i < mapped.size(); i++){
                    samRecord_t& query = mapped[i];
                    for(int j=0; j < rangeCount; j++){
                        int match = pairedCompareSam(query, query.flag.test(6) ? oracleUp : oracleDown, opt.correctRange[j]);
                        //all
                        results[j][0].addAlignment(match);
                        //unique
                        if(mapped.size()==1) results[j][1].addAlignment(match);
                        //quality values
                        for(int k=0; k < qualValCount; k++)
                            if(query.mapq >= opt.qualityValues[k])
                                results[j][2+k].addAlignment(match);
                    }
                }
                //sumarize for read
                for(int j=0; j < rangeCount; j++){
                    results[j][0].finishRead();
                    results[j][1].finishRead();
                    for(int k=0; k < qualValCount; k++)
                        results[j][2+k].finishRead();
                }
                readcnt++;//update readcnt
            }
        }
        else{
            while(!oracleFile.eof()){//next reads
                getlijn2(oracleFile, oracleLine);
                readRecord(oracleLine, oracleUp,opt.otag);
                mapped.clear();
                //iterate over the SAM file. To work query should always be more advanced than oracle
                while(!queryFile.eof() && tempQuery.qname.compare(oracleUp.qname)==0){
                    //fill mappedUp
                    mapped.push_back(tempQuery);
                    getlijn2(queryFile, queryLine);
                    readRecord(queryLine, tempQuery, opt.qtag); 
                }
                if(tempQuery.qname.compare(oracleUp.qname)==0){//add last query of file
                    mapped.push_back(tempQuery);
                }
                //process this read's alignments
                for(int i=0; i < mapped.size(); i++){
                    for(int j=0; j < rangeCount; j++){
                        int match = unpairedCompareSam(mapped[i], oracleUp, opt.correctRange[j]);
                        //all
                        results[j][0].addAlignment(match);
                        //unique
                        if(mapped.size()==1) results[j][1].addAlignment(match);
                        //quality values
                        for(int k=0; k < qualValCount; k++)
                            if(mapped[i].mapq >= opt.qualityValues[k])
                                results[j][2+k].addAlignment(match);
                    }
                }
                //sumarize for read
                for(int j=0; j < rangeCount; j++){
                    results[j][0].finishRead();
                    results[j][1].finishRead();
                    for(int k=0; k < qualValCount; k++)
                        results[j][2+k].finishRead();
                }
                readcnt++;//update readcnt
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
    for(int i=0; i < rangeCount; i++){
        outfile << "##Results for allowed range of " << opt.correctRange[i] << endl;
        outfile << "#\t\t\t alignments\t\t\t\t\t\t reads" << endl;
        outfile << "#resultType\t num_correct\t percent_correct\t num_ok\t percent_ok\t /"
                "num_incorrect\t percent_incorrect\t num_correct\t percent_correct\t /"
                "num_ok\t percent_ok\t num_incorrect\t percent_incorrect" << endl;
        outfile << "all\t " << results[i][0].printLine(readcnt);
        outfile << "unique aligned\t " << results[i][1].printLine(readcnt);
        for(int j=0; j< qualValCount; j++){
            outfile << "min mapq of " << opt.qualityValues[j] << "\t " << results[i][2+j].printLine(readcnt);
        }
        outfile << "##" << endl;
        outfile << "##" << endl;
    }
    outfile.close();
}

struct summaryLine_t {
    summaryLine_t(): alignmentCnt(0),readMapped(0),readPaired(0),readPairedCorrect(0),
    secondAln(0),mappedPairedCorrect(0), mappedPaired(0),mappedAln(0), tempReadResult(0), 
    tempAlnCount(0) {for(int i=0; i < 102; i++){alnPerRead[i]=0;}}
    long readMapped;
    long readPaired;
    long readPairedCorrect;
    long secondAln;
    long mappedPairedCorrect;
    long mappedPaired;
    long mappedAln;
    long alignmentCnt;
    int alnPerRead[102];//0-100 + 1 meaning 100+
    int tempReadResult;
    int tempAlnCount;
    void addAlignment(int flag, bool mateFound){
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
    void finishRead(){
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
    string printLine(long readcnt, bool paired){
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
            ss << "\t " << "N.A." << "\t " << "N.A.";
            ss << "\t " << secondAln << "\t " << alignmentCnt << endl;
            ss << "reads\t " << readMapped << "\t " << (readMapped*100)/readcnt;
            ss << "\t " << readPaired << "\t " << (readPaired*100)/readcnt;
            ss << "\t " << readPairedCorrect << "\t " << (readPairedCorrect*100)/readcnt;
            ss << "\t " << (readcnt-readMapped) << "\t " << ((readcnt-readMapped)*100)/readcnt;
            ss << "\t " << "N.A." << "\t " << readcnt << endl;
        }
        else{
            ss << "~~~\t num_mapped\t percent_mapped\t num_unmapped\t /"
                    "percent_unmapped\t num_second_aligned\t count" << endl;
            ss << "alignments\t " << mappedAln << "\t " << (mappedAln*100)/alignmentCnt;
            ss << "\t " << "N.A." << "\t " << "N.A.";
            ss << "\t " << secondAln << "\t " << alignmentCnt << endl;
            ss << "reads\t " << readMapped << "\t " << (readMapped*100)/readcnt;
            ss << "\t " << (readcnt-readMapped) << "\t " << ((readcnt-readMapped)*100)/readcnt;
            ss << "\t " << "N.A." << "\t " << readcnt << endl;
        }
        ss << endl;
        ss << "number of reads with x alignments" << endl;
        ss << 0 << "\t " << (readcnt-readMapped) << endl;
        for(int i = 1; i<= upperbound; i++)
            ss << i << "\t " << alnPerRead[i] << endl;
        return ss.str();
    }
};

static void checkSummary(samCheckOptions_t & opt){
    //fields: input
    ifstream queryFile(opt.querySam.c_str());
    string queryLine = "";
    streampos queryPos = queryFile.tellg();
    if(!queryFile.eof())
        getlijn2(queryFile,queryLine);//first line
    //fields: output
    int qualValCount = opt.qualityValues.size();
    summaryLine_t results[qualValCount+2];
    cerr << "summary: summary of alignments found (including pairing info and mapped/mapping quality" << endl;
    //print header of both files
    cerr << "header of SAM file: " << endl;
    while(!queryFile.eof() && !queryLine.empty() && queryLine[0]=='@'){
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
                if(results[0].tempAlnCount==1)
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
        }while(!queryFile.eof());
    }
    cerr << "compiling summary: done" << endl;
    queryFile.close();
    cerr << endl;
    cerr << "printing results" << endl;
    cerr << "Warning, these results only make sense if the query file is sorted according to query name." << endl;
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
        outfile << "##Number of reads: " << opt.numReads << endl;
        outfile << "##Reads are paired?: " << opt.paired << endl;
        outfile << "##" << endl;
        outfile << "##All alignments:" << endl;
        outfile << results[0].printLine(opt.numReads, opt.paired);
        outfile << "##" << endl;
        outfile << "##Unique alignments:" << endl;
        outfile << results[1].printLine(opt.numReads, opt.paired);
        for(int j=0; j< qualValCount; j++){
            outfile << "##Alignments with min Q-value " << opt.qualityValues[j] << endl;
            outfile << results[2+j].printLine(opt.numReads, opt.paired);
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
        cout << "##Number of reads: " << opt.numReads << endl;
        cout << "##Reads are paired?: " << opt.paired << endl;
        cout << "##" << endl;
        cout << "##All alignments:" << endl;
        cout << results[0].printLine(opt.numReads, opt.paired);
        cout << "##" << endl;
        cout << "##Unique alignments:" << endl;
        cout << results[1].printLine(opt.numReads, opt.paired);
        for(int j=0; j< qualValCount; j++){
            cout << "##Alignments with min Q-value " << opt.qualityValues[j] << endl;
            cout << results[2+j].printLine(opt.numReads, opt.paired);
        }
    }
}

long maxRange(samRecord_t & first, samRecord_t & second){
    long result = -1;
    if(first.rname.compare(second.rname)!=0 || 
            first.rnext.compare(second.rnext) !=0 || 
            (first.flag.to_ulong() % 256) != (second.flag.to_ulong() % 256))
        return result;
    else{
        result = max(first.pnext, second.pnext) - min(first.pnext, second.pnext);
        result = max(result, max(first.pos, second.pos) - min(first.pos, second.pos));
    }
    return result;
}

static void checkCompare(samCheckOptions_t & opt){
    //fields: input
    int mapperCount = opt.compareFiles.size();
    ifstream inputFiles[mapperCount];
    string inputLines[mapperCount];
    for(int i=0; i < mapperCount; i++){
        inputFiles[i].open(opt.querySam.c_str());
        inputLines[i] = "";
        if(!inputFiles[i].eof())
            getlijn2(inputFiles[i],inputLines[i]);//first line
    }    
    //fields: output
    ofstream outfile ( opt.outputFile.c_str() );
    int rangeValCount = opt.correctRange.size();
    int combinations = (1<<mapperCount);
    long resultsAln[rangeValCount][mapperCount][mapperCount];
    long resultsRead[combinations];
    
    //fields temp
    vector< vector<samRecord_t> > records(mapperCount, vector<samRecord_t>(0) );
    bitset<32> owned;
    long readCnt = 0;
    
    cerr << "summary: comparison of alignment results: which alignments are shared among SAM outputs and which are not shared." << endl;
    //print header of all Files
    cerr << "header of SAM file: PG lines" << endl;
    for(int i = 0; i < mapperCount; i++){
        cerr << "File " << i << "with name " << opt.compareFiles[i] << endl;
        while(!inputFiles[i].eof() && !inputLines[i].empty() && inputLines[i][0]=='@'){
            if(inputLines[i].substr(0, 3).compare("@PG")==0)
                cerr << inputLines[i] << endl;
            getlijn2(inputFiles[i],inputLines[i]);
        }
    }
    cerr << endl;
    cerr << "different values for the range in which alignments are considered equal: " << printVector(opt.correctRange) << endl;    
    cerr << "compiling comparison: ..." << endl;
    string tag = "NM";
    bitset<32> hasRead;
    hasRead.set();
    samRecord_t currentRec;
    bitset<32> toCompare;
    
    //fill for first reads
    for(int i = 0; i < mapperCount; i++){
        toCompare[i] = false;
        records[i].clear();
        bool sameRead = true;
        readRecord(inputLines[i], currentRec, tag);
        records[i].push_back(currentRec);
        if(!inputFiles[i].eof()){
            getlijn2(inputFiles[i],inputLines[i]);
        }
        else{
            inputLines[i]="";
        }
        while(!inputLines[i].empty() && sameRead){
            readRecord(inputLines[i], currentRec, tag);
            if(records[i][0].qname.compare(currentRec.qname)==0){
                records[i].push_back(currentRec);
                if(!inputFiles[i].eof()){
                    getlijn2(inputFiles[i],inputLines[i]);
                }
                else{
                    inputLines[i]="";
                }
            }
            else{
                sameRead = false;
            }
        }
    }
    while(hasRead.any()){
        //reset some parameters
        toCompare.reset();
        //search lexicographical smallest reads.
        string smallest="~";//Find a constant string that will be largest!!!
        for(int i = 0; i < mapperCount; i++){
            if(hasRead.test(i) && smallest.compare(records[i][0].qname) > 0){
                smallest = records[i][0].qname;
            }
        }
        //set the files which need to be checked
        owned.reset();
        for(int i = 0; i < mapperCount; i++){
            if(hasRead.test(i) && smallest.compare(records[i][0].qname) == 0){
                toCompare.set(i, true);
                for(int j=0; j < records[i].size(); j++){
                    if(records[i][j].flag.test(2)){
                        owned.set(i, true);
                        records[i][j].mapq = 0;//VERY DIRTY: SAVE OTHER INFO IN FIELD
                    }
                }
            }
        }
        //Set the read solution
        resultsRead[owned.to_ulong()]++;
        //set the alignment solution
        if(toCompare.any()){
            readCnt++;
            for(int i=0; i < mapperCount; i++){
                if(toCompare.test(i)){
                    for(int j=i+1; j < mapperCount; j++){
                        if(toCompare.test(j)){
                            for(int k = 0; k < records[i].size(); k++){
                                samRecord_t & first = records[i][k];
                                if(!first.flag.test(2)){
                                    for(int m = 0; m < records[j].size(); m++){
                                        samRecord_t & second = records[j][m];
                                        if(!second.flag.test(2)){
                                            long maxDist = maxRange(first, second);
                                            if(maxDist >= 0){// < 0 means a difference unrelated to pos
                                                for(int n = 0; n < rangeValCount; n++){
                                                    if(maxDist <= opt.correctRange[n]){
                                                        //add a combo
                                                        resultsAln[n][i][j]++;
                                                        //set not unique for these 2 alignments and this range
                                                        first.mapq = SET_BIT(first.mapq,n);
                                                        second.mapq = SET_BIT(second.mapq,n);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //Fill in the alignments unique to this output, for all ranges and all alignments
                    for(int k = 0; k < records[i].size(); k++){
                        samRecord_t & first = records[i][k];
                        if(!first.flag.test(2)){
                            for(int n = 0; n < rangeValCount; n++){
                                if(CHECK_BIT(first.mapq,n)){
                                    //set unique
                                    resultsAln[n][i][i]++;
                                }
                            }
                        }
                    }
                }
            }
        }
        //update those records that were checked.
        for(int i = 0; i < mapperCount; i++){
            if(toCompare.test(i)){
                records[i].clear();
                if(!inputLines[i].empty()){
                    bool sameRead = true;
                    readRecord(inputLines[i], currentRec, tag);
                    records[i].push_back(currentRec);
                    while(!inputFiles[i].eof() && !inputLines[i].empty() && sameRead){
                        getlijn2(inputFiles[i],inputLines[i]);
                        readRecord(inputLines[i], currentRec, tag);
                        if(currentRec.qname.compare(records[i][0].qname)==0){
                            records[i].push_back(currentRec);
                        }
                        else
                            sameRead = false;
                    }
                    if(inputFiles[i].eof() && sameRead)
                        inputLines[i].clear();
                }
                else{
                    hasRead.set(i, false);
                }
            }
        }
    }
    //extra reads that no aligner found:
    resultsRead[0] += opt.numReads-readCnt;
    cerr << "compiling comparison: done" << endl;
    for(int i=0; i < mapperCount; i++){
        inputFiles[i].close();
    }
    cerr << endl;
    cerr << "printing results" << endl;
    cerr << "Warning, these results only make sense if the all files are sorted according to read name and contain at least one alignment." << endl;
    outfile << "##ALFALFA comparison between alignments" << endl;
    outfile << "##" << endl;
    outfile << "##Warning: these results only make sense if all files are sorted according to read name and contain at least one alignment." << endl;
    outfile << "##" << endl;
    outfile << "##SAM files: " << endl;
    for(int i=0; i < opt.compareFiles.size(); i++)
        outfile << "##mapper " << i << " equals file " << opt.compareFiles[i] << endl;
    outfile << "##Distances to check: " << printVector(opt.correctRange)  << endl;
    outfile << "##Results show the number of reads that are shared among alignment files and \
    for every input distance, the number of alignments which are considered equal" << endl;
    outfile << "##" << endl;
    outfile << "##" << endl;
    outfile << "##Results" << endl;
    outfile << "##" << endl;
    outfile << "##Number of reads: " << opt.numReads << endl;
    outfile << "##Reads are paired?: " << opt.paired << endl;
    outfile << "##" << endl;
    outfile << "Reads aligned by none: " << resultsRead[0] << endl;
    outfile << "##Reads uniquely aligned by mappers:" << endl;
    outfile << "mapper\t number_reads" << endl;
    for(int i=0; i < mapperCount; i++){
        outfile << i << "\t" << resultsRead[(1<<i)] << endl;
    }
    outfile << "##Reads aligned by multiple mappers:" << endl;
    outfile << "mappers\t number_reads" << endl;
    for(int i=3; i< combinations; i++){
        int currentCombi = i;
        int mapperIndex=0;
        while(currentCombi > 0 && (currentCombi & 1)==0){
            mapperIndex++;
            currentCombi = currentCombi >> 1;
        }
        if(currentCombi > 1){//at least 2 mappers, not power of 2
            outfile << mapperIndex;
            mapperIndex++;
            currentCombi = currentCombi >> 1;
            while(currentCombi > 0){
                if((currentCombi & 1) == 1)
                    outfile << "," << mapperIndex;
                mapperIndex++;
                currentCombi = currentCombi >> 1;
            }
            outfile << "\t" << resultsRead[i] << endl;
        }
    }
    outfile << "##" << endl;
    outfile << "##For every distance: the pairwise shared alignments" << endl;
    outfile << "##" << endl;
    for(int i=0; i< rangeValCount; i++){
        outfile << "##Distance: " << opt.correctRange[i] << endl;
        outfile << "Mapper";
        for(int j=0; j < mapperCount; j++){
            outfile << "\t" << j;
        }
        outfile << endl;
        for(int j=0; j < mapperCount; j++){//rows
            outfile << j;
            for(int k=0; k < mapperCount; k++){//columns
                outfile << "\t";
                if(k >= j){
                    outfile << resultsAln[i][j][k];
                }
            }
            outfile << endl;
        }
        outfile << "##" << endl;
        outfile << "##" << endl;
    }
    outfile.close();
}

#endif	/* PERFORMANCEUTILS_H */

