/*
 * File:   fasta.h
 * Author: mvyvermn
 *
 * Created on 31 augustus 2011, 15:30
 *
 */

#ifndef FASTA_H
#define	FASTA_H

#include <string>
#include <vector>

using namespace std;

void reverse_complement(string &seq_rc, bool nucleotides_only);
void trim(string &line, long &start, long &end);
void load_fasta(string filename, string &S, vector<string> &descr, vector<long> &startpos);


#endif	/* FASTA_H */

