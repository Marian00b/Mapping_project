#ifndef SEQA_H
#define SEQA_H

#include "Sequence_FastX.h"
#include "FakeSeq.h"

class Sequence_FastA : public Sequence_FastX
{
	protected:
		unsigned char * encodage;

	public:
		Sequence_FastA();
		Sequence_FastA(string titre, size_t lg);
		Sequence_FastA(string titre, const char * const seq);
		Sequence_FastA(string titre, size_t lg, unsigned char * seq);
		Sequence_FastA(const Sequence_FastA* const &s);
		Sequence_FastA(Sequence_FastA const &s);
		Sequence_FastA* complement() const ;
		Sequence_FastA* reverse() const ;
		Sequence_FastA* subSequence(int begin, int end = 0) const;
		string operator[](size_t i) const; 
		Sequence_FastA &operator=(Sequence_FastA const &s) ;
		bool operator==(Sequence_FastA const &s);
		void setSymbolAt(size_t i, char c);
	 	virtual void print() const;	
		unsigned char * getEncodage() const; 
	    ~Sequence_FastA();
	    friend ostream& operator<<(ostream &os, Sequence_FastA const& es);

};

string decode(char num);
size_t donneOctet (size_t pos, size_t lg);
size_t donnePosOctet (size_t pos); 
unsigned char encode(char nuc) ;
bool estValide(char nuc);
unsigned char donneNext (size_t octet, size_t position,  unsigned char * const& encodage);
Sequence_FastA* createIntitule (string nomFicFasta, size_t position, char separator);
Sequence_FastA* createSequence (string nomFicFasta, FakeSeq * seq);
vector <Sequence_FastA*> createAll (string nomFicFasta, vector <size_t> positions, size_t debut, size_t fin, bool isSeq, vector <FakeSeq*> seq = vector<FakeSeq*>());
vector <size_t> extractPosSeq (string nomFicFasta);
unsigned char * remplissageOctet(unsigned char * const& encodage, int pos_octet, size_t octet, int lg);
Sequence_FastA* concatSeq (Sequence_FastA* s, Sequence_FastA* iseq);
Sequence_FastA* concatAllSeq (vector <Sequence_FastA*> sequences);


#endif
