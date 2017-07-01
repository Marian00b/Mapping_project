#include "Sequence_FastQ.h"

Sequence_FastQ::Sequence_FastQ(): Sequence_FastA(), encodage_Q(NULL) {}

Sequence_FastQ::Sequence_FastQ(string titre, size_t lg): Sequence_FastA(titre, lg),encodage_Q(NULL) {}

Sequence_FastQ::Sequence_FastQ(string titre, const char * seq, const char * seq_Q): Sequence_FastA(titre, seq) {
	encodage_Q = new  char [longueur];

	for (size_t i=0; i < longueur; i++) {
		encodage_Q[i] = seq_Q[i];
	}

}

Sequence_FastQ::Sequence_FastQ(const Sequence_FastA* const s, const char * seq_Q): 
Sequence_FastA(s), encodage_Q(new char[longueur]) {
	for (size_t i=0; i < longueur; i++) {
		encodage_Q[i] = seq_Q[i];
	}
}

Sequence_FastQ::Sequence_FastQ(string titre, size_t lg, unsigned char * seq, char * seq_Q): 
 Sequence_FastA(titre, lg, seq), encodage_Q(new char[longueur])  {

 	for (size_t i=0; i < longueur; i++) {
		encodage_Q[i] = seq_Q[i];
	}
}

Sequence_FastQ::Sequence_FastQ(Sequence_FastQ* const& s)
:Sequence_FastA(s) {
	encodage_Q = new char[longueur];
	for (size_t i= 0; i < longueur; i++){
		encodage_Q[i] = s->encodage_Q[i];
	}
}


Sequence_FastQ::~Sequence_FastQ(){
	if (encodage_Q != NULL) 
		delete[] encodage_Q;
}

void Sequence_FastQ::print() const{

	Sequence_FastA::print(); 

	if (encodage_Q != NULL) {

		cout<<endl<<"Sequence qualité : "<<endl;

		for (size_t i=0; i < longueur; i++) {
			cout<<encodage_Q[i];
		}
		cout<<endl;
	}
}

Sequence_FastQ* Sequence_FastQ::complement() const {
 	Sequence_FastA* s = Sequence_FastA::complement();
	Sequence_FastQ* sQ= new Sequence_FastQ(s, encodage_Q);

	delete s; 
	return  sQ;
}

Sequence_FastQ* Sequence_FastQ::reverse() const {
	Sequence_FastA* s = Sequence_FastA::reverse();
	char * seq_Q = new char [longueur];
	for (size_t i = 0; i < longueur; i++ ){
		seq_Q[i] = encodage_Q[longueur-i-1];
	}
	Sequence_FastQ* sQ= new Sequence_FastQ(s, seq_Q);

	delete[] seq_Q;
	delete s; 
	return  sQ;
}


Sequence_FastQ* Sequence_FastQ::subSequence(int begin, int end) const {
	if (begin >=  int(longueur)) {
		cout<<"Sous sequence impossible."<<endl;
		return new Sequence_FastQ(intitule, longueur, encodage, encodage_Q);
	}
	else {
		if (begin < 0) {
			begin = 0;
		}
		Sequence_FastA* s = Sequence_FastA::subSequence(begin, end);

		char * sub_Q = new char[s->length()];
		int i = begin;

		for (size_t j = 0; j < s->length(); j++){
			sub_Q[j] = encodage_Q[i];
			i++;
		}

		Sequence_FastQ* sQ= new Sequence_FastQ(s, sub_Q);

		delete s;
		delete[] sub_Q;
		return sQ;
	}
}

void Sequence_FastQ::supprimePrefixe(size_t pref) {
 	Sequence_FastQ* sQ = this->subSequence(pref,longueur-1);
	longueur = sQ->length();
	delete[] encodage;
	encodage = sQ->encodage;
	delete[] encodage_Q;
	encodage_Q = sQ->encodage_Q;
}

void Sequence_FastQ::supprimeSuffixe(size_t suff) {
	if (suff > longueur) {
		suff = longueur;
	}
	Sequence_FastQ* sQ = this->subSequence(0,longueur-suff-1);
 	longueur = sQ->length();
 	delete[] encodage;
 	encodage = sQ->encodage; 
 	delete[] encodage_Q;
 	encodage_Q = sQ->encodage_Q;
} 

void Sequence_FastQ::supprimeSousSequence(size_t begin, size_t  end){
	// rajouter si inverse!
	// if (begin < 0) {
	// 	begin = 0;
	// }
	if (end >= longueur){
		end = longueur-1; 
	}

	if (end == longueur-1) {
		this->supprimeSuffixe(end-begin+1);
	}
	else {
		if (begin == 0) {
			this->supprimePrefixe(end+1);
		}
		else {

			Sequence_FastQ* sP = this->subSequence(0,begin-1);	
	
			Sequence_FastQ* sS = this->subSequence(end+1,longueur-1);
			delete[] encodage;
			encodage = joinEncodage(sP, sS, this, end-begin+2); 

			longueur = sP->length() + sS->length();

			char * seq_Q = new  char [longueur];

			size_t i_seq = 0; 
			size_t i_enco = 0; 
			while (i_seq <  longueur) {
				if (i_enco == sP->length()) {
					i_enco += end-begin+1;
				}
				seq_Q[i_seq] = encodage_Q[i_enco];
				i_enco++;
				i_seq++;
			} 

			delete[] encodage_Q; 

			encodage_Q = new char [longueur]; 

			for (size_t i = 0; i < longueur; i++) {
				encodage_Q[i] = seq_Q[i];
			}

			delete[] seq_Q;

			delete sS;
			delete sP;
		}
	}
}

// parametre longueur polyA
void Sequence_FastQ::supprimePolyA (size_t polyA) {
	size_t i = longueur - 1;
	while ((*this)[i] == "A") {
		i--;
	}

	if (longueur - 1 - i >= polyA) {
		supprimeSuffixe(longueur-i-1);
	}
}

size_t Sequence_FastQ::supprimeNextBadQuality(size_t lg_min, size_t q_min) {
	int bad_q = -1;
	size_t i = 0;
	bool hasBeenFound = false; 

	// attention -33! 
	while (i < longueur && !hasBeenFound) {
		if (int(encodage_Q[i]) <= int(q_min+33)) {
			if (bad_q == -1) {
				bad_q = i;  
			}
		}
		else {
			if (bad_q != -1 && i-bad_q > lg_min) {
				supprimeSousSequence(bad_q, i-1);
				hasBeenFound = true; 
			}
			bad_q = -1;
		}
		i++;
	}

	if (!hasBeenFound && bad_q != -1 && i-bad_q > lg_min) {
		supprimeSousSequence(bad_q, i-1);
	}
	return i == longueur ? -1 : i ;
}

// parametre qualite et longueur
void Sequence_FastQ::supprimeBadQuality (size_t lg_min, size_t q_min) {
	int i = 0;
	while (i != -1) {
		i = supprimeNextBadQuality(lg_min, q_min);
	}
}



