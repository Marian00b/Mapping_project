#ifndef SEQQ_H
#define SEQQ_H

#include "Sequence_FastA.h"
#include "FakeSeq.h"

class Sequence_FastQ : public Sequence_FastA
{
	private:
		char * encodage_Q;

	public:
		Sequence_FastQ();
		Sequence_FastQ(string titre, size_t lg);
		Sequence_FastQ(string titre, const char * seq, const char * seq_Q);
	    Sequence_FastQ(string titre, size_t lg, unsigned char * seq, char * seq_Q);
	    Sequence_FastQ(const Sequence_FastA* const s, const char * seq_Q);
	    Sequence_FastQ* complement() const ;
		Sequence_FastQ* reverse() const ;
		Sequence_FastQ* subSequence(int begin, int end = 0) const;
		Sequence_FastQ(Sequence_FastQ* const& s);
		void supprimePrefixe(size_t pref);
		void supprimeSuffixe(size_t suff); 
		void supprimeSousSequence(size_t begin, size_t end); 
		void supprimePolyA (size_t polyA);
		size_t supprimeNextBadQuality(size_t lg_min, size_t q_min);
		void supprimeBadQuality (size_t lg_min, size_t q_min);
	  	void print() const ;
	  	bool operator==(Sequence_FastQ const &s);
	  	Sequence_FastQ &operator=(Sequence_FastQ const &s);
		~Sequence_FastQ();
		void printQualite() const;
		char* getQualite() const;
};

Sequence_FastQ* createIntituleQ (string nomFicFasta, size_t position);
Sequence_FastQ* createSequenceQ (string nomFicFasta,  FakeSeq * seq);
vector <size_t> extractPosSeqQ (string nomFicFasta);
int next_sequence (string nomFicFasta, int position);
vector <Sequence_FastQ*> createAllQ (string nomFicFasta, vector <size_t> positions, size_t debut, size_t fin, bool isSeq, vector <FakeSeq*> & seq);
unsigned char *  joinEncodage(Sequence_FastQ* const& s1, Sequence_FastQ* const& s2, Sequence_FastQ* const& s0, size_t milieu);
Sequence_FastQ* concatAllSeqQ (vector <Sequence_FastQ*> sequences);

#endif
