#ifndef SA_H
#define SA_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <map>
#include "Sequence_FastQ.h"
#include "kmer.h"

class SAcomp {

	private : 
		Sequence_FastA* const& seqA;
		std::vector <size_t> SA;
		std::vector <size_t> LCP;

	public :
		SAcomp();
		SAcomp(Sequence_FastA* const& seq); 
		struct SATest;
		void buildLCP();
		size_t length();
		vector <size_t> const& getSA() const;
		void print() const;
		int dichoM (string W);
		int dicho (Sequence_FastA* const& W);
		Sequence_FastA* getIseq(size_t i);
		Sequence_FastA* getIfacteur(size_t pos, size_t lg);
		vector<int> rechercheAllM (string W); 
		vector<int> rechercheAll (Sequence_FastA* const& W); 
		size_t iLength(size_t i); 
		Sequence_FastA* getOccurences(size_t lg);
		size_t getDistincts(size_t lg);
		Sequence_FastA* getMot (size_t k, size_t pos);
		map< int, vector<int>> rechercheBrins(Sequence_FastA* const& W);
		map <int,size_t> searchAllKmers(vector<Kmer*> & kmers);
		int* rechercheK (int pos, Sequence_FastA* const& kmer);
		vector<int> rechercheMMP(Sequence_FastA* const& kmer);
		void  getMappingKmer (map <int, size_t> & res, vector<Kmer*> const & kmers, ofstream & log_file); 
		size_t prefcommSA(size_t i, Sequence_FastA* const& seq);
		size_t SprefcommSA(size_t i, string W);
		size_t prefcommSuff(size_t i, size_t j);



};
//size_t prefcommLCP(size_t i, size_t j);
size_t Sprefcomm(size_t i, Sequence_FastA const& seq, string W);
size_t prefcomm(size_t i, Sequence_FastA* const& seq1, Sequence_FastA* const& seq2);
string char2string (char c);
void printLocalise(unsigned char i,unsigned char j );


#endif