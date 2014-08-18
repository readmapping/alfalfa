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
 *
 * Part of the code originates from "A Practical Algorithm for Finding
 * Maximal Exact Matches in Large Sequence Data Sets Using Sparse Suffix Arrays"
 * By Khan et al. Copyright (c) 2009
 * You should have received a copy of the copyright notice with this code.
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "options.h"

using namespace std;

void align_opt::fill_m() {
    int i, j, k;
    for (i = k = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j)
            mat[k++] = i == j ? match : mismatch;
        mat[k++] = -1; // ambiguous base
    }
    for (j = 0; j < 5; ++j) mat[k++] = -1;
}

void mapOptions_t::initOptions() {
    K = 12;
    query_threads = 1;
    nucleotidesOnly = false;
    alnOptions.minMemLength = 40;
    alnOptions.fixedMinLength = false;
    alnOptions.memType = SMAM;
    alnOptions.sparseMult = 1;
    alnOptions.noFW = alnOptions.noRC = false;
    alnOptions.errorPercent = 0.08;
    alnOptions.print = 0;
    alnOptions.debugFile = NULL;
    alnOptions.match = 1; //(2,-1,-2,-1);
    alnOptions.mismatch = -4;
    alnOptions.openGap = -6;
    alnOptions.extendGap = -1;
    alnOptions.global = true;
    alnOptions.alignmentCount = 1;
    alnOptions.maxTrial = 10;
    alnOptions.minTrial = 5000;
    alnOptions.maxTrialBest = true;
    alnOptions.alwaysUnique = true;
    alnOptions.minCoverage = 0.25;
    alnOptions.rescue = true;
    alnOptions.noSparseQ = false;
    alnOptions.fullDP = false;
    alnOptions.maxSeedCandidates = 10000;
    alnOptions.maxSmemMemStart = 10;
    alnOptions.maxRegionMems = 20;
    alnOptions.minMemRegLength = 50; //SET TO AUTOMATIC VALUE
    alnOptions.scoresSetByUser = false;
    alnOptions.fixedBandSize = 100;
    pairedOpt.contain = true;
    pairedOpt.discordant = true;
    pairedOpt.dovetail = false;
    pairedOpt.maxInsert = 1000; //return to default!//
    pairedOpt.minInsert = 0; //return to default!//
    pairedOpt.mixed = true;
    pairedOpt.orientation = PAIR_FR;
    pairedOpt.overlap = true;
    pairedOpt.mode = 1;
    pairedOpt.pairedRescue = false;
    unpairedQ = pair1 = pair2 = ref_fasta = outputName = "";
    command = ALN;
    saveIndex = false;
    hasChild = true;
    hasSuflink = false;
    hasKmer = true;
    index_prefix = indexLocation = "";
}

void mapOptions_t::printOptions() {
    cerr << "ALFALFA was executed with following options: " << endl;
    cerr << "reference sequence\t\t" << ref_fasta << endl;
    if (!outputName.empty())
        cerr << "output file\t\t" << outputName << endl;
    if (!indexLocation.empty())
        cerr << "index location\t\t\t" << indexLocation << endl;
    if (!index_prefix.empty())
        cerr << "index prefix\t\t\t" << index_prefix << endl;
    cerr << "save index \t\t\t" << (saveIndex ? "yes" : "no") << endl;
    if (!unpairedQ.empty())
        cerr << "single-end reads\t\t\t" << unpairedQ << endl;
    if (!pair1.empty() && !pair2.empty()) {
        cerr << "mate 1   \t\t\t" << pair1 << endl;
        cerr << "mate 2   \t\t\t" << pair2 << endl;
    }
    cerr << "sparseness \t\t\t" << K << endl;
    cerr << "threads  \t\t\t" << query_threads << endl;
    cerr << "edit distance   \t\t\t" << alnOptions.errorPercent << endl;
    cerr << "min seed length\t\t\t" << alnOptions.minMemLength << endl;
    cerr << "auto seed length\t\t\t" << (alnOptions.fixedMinLength ? "no" : "yes") << endl;
    cerr << "max seed candidates\t\t\t" << alnOptions.maxSeedCandidates << endl;
    cerr << "max SMEM occurrence before mems are calculated " << alnOptions.maxSmemMemStart << endl;
    cerr << "min MEM length to consider for region selection " << alnOptions.minMemRegLength << endl;
    cerr << "max number of MEMs per SMEM for region selection " << alnOptions.maxRegionMems << endl;
    cerr << "type of seeds\t\t\t" << alnOptions.memType << endl;
    cerr << "rescue procedures on?\t\t\t" << (alnOptions.rescue ? "yes" : "no") << endl;
    cerr << "full dynamic programming for post-processing\t" << (alnOptions.fullDP ? "yes" : "no") << endl;
    cerr << "#alignments\t\t\t" << alnOptions.alignmentCount << endl;
    cerr << "#try and error\t\t\t" << alnOptions.maxTrial << endl;
    cerr << "min #alignments to decide mapq\t\t" << alnOptions.minTrial << endl;
    cerr << "min #alignments checks from best alignment?\t\t" << (alnOptions.maxTrialBest ? "yes" : "no") << endl;
    cerr << "always consider regions containing unique SMEMs\t\t" << (alnOptions.alwaysUnique ? "yes" : "no") << endl;
    cerr << "min query coverage\t\t" << alnOptions.minCoverage << endl;
    cerr << "use local alignment?\t\t\t" << (alnOptions.global ? "no" : "yes") << endl;
    if (!alnOptions.global) {
        cerr << "use band size\t\t\t" << alnOptions.fixedBandSize << endl;
    }
    cerr << "verbosity level\t\t\t" << alnOptions.print << endl;
    cerr << "debug info to file?" << (alnOptions.debugFile == NULL ? "no" : "yes") << endl;
    cerr << "calculate forward\t\t\t" << (alnOptions.noFW ? "no" : "yes") << endl;
    cerr << "calculate reverse\t\t\t" << (alnOptions.noRC ? "no" : "yes") << endl;
    cerr << "match score\t\t\t" << alnOptions.match << endl;
    cerr << "mismatch score\t\t\t" << alnOptions.mismatch << endl;
    cerr << "gap open score\t\t\t" << alnOptions.openGap << endl;
    cerr << "gap extend score\t\t\t" << alnOptions.extendGap << endl;
    if (!pair1.empty() && !pair2.empty()) {
        cerr << "paired options: " << endl;
        cerr << "paired mode\t\t\t" << pairedOpt.mode << endl;
        cerr << "type pairing\t\t\t" << pairedOpt.orientation << endl;
        cerr << "min insert size\t\t\t" << pairedOpt.minInsert << endl;
        cerr << "max insert size\t\t\t" << pairedOpt.maxInsert << endl;
        cerr << "allow overlap\t\t\t" << (pairedOpt.overlap ? "yes" : "no") << endl;
        cerr << "allow contained\t\t\t" << (pairedOpt.contain ? "yes" : "no") << endl;
        cerr << "allow discordant\t\t\t" << (pairedOpt.discordant ? "yes" : "no") << endl;
        cerr << "allow dovetail\t\t\t" << (pairedOpt.dovetail ? "yes" : "no") << endl;
        cerr << "allow mixed\t\t\t" << (pairedOpt.mixed ? "yes" : "no") << endl;
        cerr << "paired rescue?\t\t\t" << (pairedOpt.pairedRescue ? "yes" : "no") << endl;
    }
}

mum_t parseSeedType(string type) {
    if (type.compare("MEM") == 0)
        return MEM;
    else if (type.compare("MAM") == 0)
        return MAM;
    else if (type.compare("MUM") == 0)
        return MUM;
    else if (type.compare("SMEM") == 0)
        return SMAM;
    else if (type.compare("SMAM") == 0)
        return SMAM;
    else if (type.compare("PSMEM") == 0)
        return SMEM;
    else {
        cerr << "WARNING: unrecognized seed type. SMEMs are used instead." << endl;
        return SMAM;
    }
}

orientation_t parseOrientation(string type) {
    if (type.compare("fr") == 0)
        return PAIR_FR;
    else if (type.compare("rf") == 0)
        return PAIR_RF;
    else if (type.compare("ff") == 0)
        return PAIR_FF;
    else {
        cerr << "WARNING: unrecognized orientation. fr is used instead." << endl;
        return PAIR_FR;
    }
}

void processParameters(int argc, char* argv[], mapOptions_t& opt) {
    //parse command
    if (strcmp(argv[1], "index") == 0) {
        opt.command = INDEX;
        opt.saveIndex = true;
    } else if (strcmp(argv[1], "align") == 0) {
        opt.command = ALN;
    } else {
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[0]);
        option::printUsage(std::cerr, basicUsage);
    }

    cerr << "parsing options: ..." << endl;
    argc -= 2;
    argv++;
    argv++;
    option::Stats stats(fullUsage, argc, argv);
    //misc option parsing
    option::Option* options = new option::Option[stats.options_max];
    option::Option* buffer = new option::Option[stats.buffer_max];
    option::Parser parser(fullUsage, argc, argv, options, buffer);
    if (parser.error()) {
        cerr << "error during parsing of command line options" << endl;
        exit(1);
    }
    if (options[ARG_HELP]) {
        if (opt.command == INDEX) {
            option::printUsage(std::cerr, indexUsage);
        } else {
            option::printUsage(std::cerr, alignUsage);
        }
        option::printUsage(std::cerr, miscUsage);
        exit(0);
    }
    if (options[ARG_VERBOSE])
        opt.alnOptions.print = atoi(options[ARG_VERBOSE].arg);
    if (options[ARG_DEBUG])
        opt.alnOptions.debugFile = fopen(options[ARG_DEBUG].arg, "w");
    //index parsing options
    for (int i = 0; i < parser.optionsCount(); ++i) {
        option::Option& option = buffer[i];
        switch (option.index()) {
            case ARG_REFERENCE: opt.ref_fasta = option.arg;
                break;
            case ARG_SPARSENESS: opt.K = atoi(option.arg);
                break;
            case ARG_PREFIX: opt.index_prefix = option.arg;
                break;
            case ARG_NO_CHILD: opt.hasChild = false;
                break;
            case ARG_SUFLINK: opt.hasSuflink = true;
                break;
            case ARG_NO_KMER: opt.hasKmer = false;
                break;
            default: break;
                throw 1;
        }
    }
    //check needed options and correctness
    if (opt.command == INDEX && !options[ARG_REFERENCE]) {
        cerr << "ERROR: no reference file given." << endl;
        exit(1);
    }
    if (opt.command == INDEX) {
        if (opt.index_prefix.empty()) {
            opt.index_prefix = opt.ref_fasta;
        }
        if (opt.hasChild && opt.hasSuflink)
            fprintf(stderr, "WARNING: having both child array and suffix link support is marginally useful and requires more memory! \n");
        if (!opt.hasChild && opt.K >= 3)
            fprintf(stderr, "WARNING: no child array means much slower search times! \n");
    }
    if (opt.command == ALN && !(options[ARG_REFERENCE] || options[ARG_INDEX])) {
        cerr << "ERROR: no reference or index name given." << endl;
        exit(1);
    }
    //mapping options parsing
    if (opt.command == ALN) {
        for (int i = 0; i < parser.optionsCount(); ++i) {
            option::Option& option = buffer[i];
            switch (option.index()) {
                case ARG_INDEX: opt.indexLocation = option.arg;
                    break;
                case ARG_SAVE: opt.saveIndex = true;
                    break;
                case ARG_SINGLE: opt.unpairedQ = option.arg;
                    break;
                case ARG_MATES1: opt.pair1 = option.arg;
                    break;
                case ARG_MATES2: opt.pair2 = option.arg;
                    break;
                case ARG_OUTPUT: opt.outputName = option.arg;
                    break;
                case ARG_ALIGNMENTS: opt.alnOptions.alignmentCount = atoi(option.arg);
                    break;
                case ARG_NO_FORWARD: opt.alnOptions.noFW = 1;
                    break;
                case ARG_NO_REVERSE: opt.alnOptions.noRC = 1;
                    break;
                case ARG_EDIT_DISTANCE: opt.alnOptions.errorPercent = atof(option.arg) + 0.01;
                    break;
                case ARG_THREADS: opt.query_threads = atoi(option.arg);
                    break;
                case ARG_SEED: opt.alnOptions.memType = parseSeedType(option.arg);
                    break;
                case ARG_MIN_LENGTH: opt.alnOptions.minMemLength = atoi(option.arg);
                    opt.alnOptions.fixedMinLength = true;
                    break;
                case ARG_MAX_SEEDS: opt.alnOptions.maxSeedCandidates = atoi(option.arg);
                    break;
                case ARG_MAX_SMEM: opt.alnOptions.maxSmemMemStart = atoi(option.arg);
                    break;
                case ARG_MAX_MEM: opt.alnOptions.maxRegionMems = atoi(option.arg);
                    break;
                case ARG_MIN_MEML: opt.alnOptions.minMemRegLength = atoi(option.arg);
                    break;
                case ARG_NO_RESCUE: opt.alnOptions.rescue = false;
                    break;
                case ARG_NO_SPARSEQ: opt.alnOptions.noSparseQ = true;
                    break;
                case ARG_FULLDP: opt.alnOptions.fullDP = true;
                    break;
                case ARG_MAX_FAILURES: opt.alnOptions.maxTrial = atoi(option.arg);
                    break;
                case ARG_MAX_F_BEST: opt.alnOptions.maxTrialBest = false;
                    break;
                case ARG_MIN_ALIGN: opt.alnOptions.minTrial = atoi(option.arg);
                    break;
                case ARG_MIN_COVERAGE: opt.alnOptions.minCoverage = atof(option.arg);
                    break;
                case ARG_ALWAYS_UNIQUE: opt.alnOptions.alwaysUnique = false;
                    break;
                case ARG_LOCAL: opt.alnOptions.global = false;
                    break;
                case ARG_BANDWIDTH: opt.alnOptions.fixedBandSize = atoi(option.arg);
                    break;
                case ARG_MATCH: opt.alnOptions.match = atoi(option.arg);
                    opt.alnOptions.scoresSetByUser = true;
                    break;
                case ARG_MISMATCH: opt.alnOptions.mismatch = atoi(option.arg);
                    opt.alnOptions.scoresSetByUser = true;
                    break;
                case ARG_GAP_OPEN: opt.alnOptions.openGap = atoi(option.arg);
                    opt.alnOptions.scoresSetByUser = true;
                    break;
                case ARG_GAP_EXTEND: opt.alnOptions.extendGap = atoi(option.arg);
                    opt.alnOptions.scoresSetByUser = true;
                    break;
                case ARG_MIN_INSERT: opt.pairedOpt.minInsert = atoi(option.arg);
                    break;
                case ARG_MAX_INSERT: opt.pairedOpt.maxInsert = atoi(option.arg);
                    break;
                case ARG_ORIENTATION: opt.pairedOpt.orientation = parseOrientation(option.arg);
                    break;
                case ARG_NO_MIXED: opt.pairedOpt.mixed = false;
                    break;
                case ARG_NO_DISCORDANT: opt.pairedOpt.discordant = false;
                    break;
                case ARG_DOVETAIL: opt.pairedOpt.dovetail = true;
                    break;
                case ARG_NO_CONTAIN: opt.pairedOpt.contain = false;
                    break;
                case ARG_NO_OVERLAP: opt.pairedOpt.overlap = false;
                    break;
                case ARG_PAIRED_MODE: opt.pairedOpt.mode = atoi(option.arg);
                    break;
                case ARG_PAIRED_RESCUE: opt.pairedOpt.pairedRescue = true;
                    break;
                default: break;
                    throw 1;
            }
        }
        if (opt.indexLocation.empty() && opt.index_prefix.empty()) {
            opt.index_prefix = opt.ref_fasta;
        }
        if (opt.hasChild && opt.hasSuflink)
            fprintf(stderr, "WARNING: having both child array and suffix link support is marginally useful and requires more memory! \n");
        if (!opt.hasChild && opt.K >= 3)
            fprintf(stderr, "WARNING: no child array means much slower search times! \n");
        if (!opt.saveIndex && opt.indexLocation.empty())
            fprintf(stderr, "WARNING: index will not be saved to disk! \n");
        if (opt.unpairedQ.empty() && (opt.pair1.empty() || opt.pair2.empty())) {
            fprintf(stderr, "ERROR: ALFALFA requires query files. Specify query files with -0 or -1,-2 options.\n");
            exit(1);
        }
        if (opt.outputName.empty())
            opt.outputName = opt.ref_fasta.substr().append(".sam");
        if (opt.alnOptions.errorPercent < 0.0 || opt.alnOptions.errorPercent > 1.0) {
            fprintf(stderr, "ERROR: uncorrect value for edit distance. Value should lie in the interval [0,1[\n");
            exit(1);
        }
        if (opt.alnOptions.minCoverage < 0.0 || opt.alnOptions.minCoverage > 1.0) {
            fprintf(stderr, "ERROR: uncorrect value for minimum coverage. Value should lie in the interval [0,1[\n");
            exit(1);
        }
        if (opt.pairedOpt.minInsert > opt.pairedOpt.maxInsert) {
            fprintf(stderr, "ERROR: minimum insert size is larger than maximum insert size\n");
            exit(1);
        }
        if (opt.pairedOpt.mode < 1 || opt.pairedOpt.mode > 6) {
            fprintf(stderr, "ERROR: chosen paired-end mode not supported.\n");
            exit(1);
        }
        opt.alnOptions.maxAlnCount = max(opt.alnOptions.alignmentCount, opt.alnOptions.minTrial);
        opt.alnOptions.fill_m();
    }
    if (opt.alnOptions.print >= 1)
        opt.printOptions();
    delete[] options;
    delete[] buffer;
}

