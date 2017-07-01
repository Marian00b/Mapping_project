#ifndef SEQF_H
#define SEQF_H

#include "Sequence_FastX.h"

class FakeSeq : public Sequence_FastX
{
	private:
		size_t pos_encodage;
		size_t pos_encodageQ; 

	public:
		FakeSeq(string titre = "", size_t lg = 0, size_t seq = 0, size_t seq_Q = 0);
		FakeSeq(const FakeSeq* const& s);
	  	virtual void print() const;
	  	bool operator==(FakeSeq const& s);
	  	FakeSeq &operator=(FakeSeq const& s);
	  	size_t getPos_encodage() const;
	  	size_t getPos_encodageQ() const;
		~FakeSeq();
};

FakeSeq* createFakeSeqQ (string nomFicFasta, size_t position);
FakeSeq* createFakeSeqA (string nomFicFasta, size_t position); 
vector <FakeSeq*> createAllFakeSeq (string nomFicFasta, vector <size_t> positions, size_t debut, size_t fin, bool isSeqQ);

#endif 
