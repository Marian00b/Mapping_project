#ifndef K_H
#define K_H

#include <map>
#include "suffixe.h"
#include <cmath>

class Kmer {

	private : 
		vector<Sequence_FastA*> kmers; 
		size_t longueur;
		//vector <size_t> ordre;
		map <size_t, vector<int>> exact_match; 
		map <size_t, vector<int>> non_match ; 
		int localisation; 
		bool isLocalise_exact ;
		bool isLocalise_NonExact ;

	public :
		Kmer(); 
		Kmer(Sequence_FastA* const& seq, size_t lg);
		//Kmer(vector<Sequence_FastA*> const& k); 
		void print(); 
		bool isKmer (Sequence_FastA* const&  seq) const;
		//void insertKmer (Sequence_FastA* const& seq, size_t nb);
		Sequence_FastA* const& getKmer(size_t i);
		size_t length() const;
		~Kmer();
		size_t nombre() const;
		void setMatch(size_t i, int v);
		void setNonMatch(size_t i, int v);
		void setNewNonMatch(size_t i, vector<int> v);
		vector<int> & getPositions(size_t i);
		int isSuite (map <int, size_t> const& res, int ecart);
		bool getIsLocalise_exact ();
		bool getIsLocalise_NonExact ();
		int getLocalisation();

};


vector<Kmer*> createAllKmers (vector<Sequence_FastQ*> const sequences, size_t longueur);
bool sameBrin (int nb1, int nb2);
void find_closer(int pos, vector<int> & v);
map< Sequence_FastA*, int> getNonLocalise(vector<Kmer*> Reads);
//void insertKmer (vector<size_t> & ordre, vector<Sequence_FastA*> & non_localise, Sequence_FastA* const& seq, size_t nb);

#endif