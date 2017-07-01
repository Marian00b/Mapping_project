#include "Sequence_FastQ.h"

unsigned char *  joinEncodage(Sequence_FastQ* const& s1, Sequence_FastQ* const& s2, Sequence_FastQ* const& s0, size_t milieu) {
	size_t longueur = s1->length() + s2->length();
	size_t final_octet = donneOctet(s1->length()-1, s1->length());
	int final_position = donnePosOctet(s1->length()-1);
	size_t i_tab; 
	int i = 6;
	unsigned char next = 0; 
	unsigned char * new_encodage = new unsigned char [longueur/4 + (longueur%4 != 0)]; 

	for (i_tab = 0; i_tab < final_octet; i_tab++){
		new_encodage[i_tab] = s1->getEncodage()[i_tab];
	}

	// remplissage de l'octet particulier avec le premier encodage
	while (i >= final_position) {
		next<<= 2; 
		next = next | donneNext(final_octet, i ,s1->getEncodage());
		i-=2;
	}

	final_position = 6; 

	//remplissage avec le second
	while( i != -2 ){
		next<<= 2; 
		next = next | donneNext(final_octet, final_position ,s2->getEncodage());
		final_position-=2;
		i-=2;
	}

	new_encodage[final_octet] = next;

	if (final_position == -2) {
		final_octet +=1;
		final_position = 6;
	}

	// recuperation position dans encodage

	size_t position = s1->length() + milieu + (3-final_position/2) - 1;


	unsigned char * enco2 = remplissageOctet(s0->getEncodage(), donnePosOctet(position), donneOctet(position, s0->length()),s2->length());

	for (size_t j = i_tab+1 ; j < (longueur/4 + (longueur%4 != 0)); j++) {
		new_encodage[j] = enco2[j-(i_tab+1)];
	} 

	delete[] enco2;

	return new_encodage;
}

Sequence_FastQ* createIntituleQ (string nomFicFasta, size_t position) {
	Sequence_FastA* sA = createIntitule(nomFicFasta, position, ':');

	Sequence_FastQ* sQ = new Sequence_FastQ(sA->getIntitule(), sA->length());

	delete sA;
	return sQ;

}

Sequence_FastQ* createSequenceQ (string nomFicFasta, FakeSeq* seq) {
	size_t longueur = 0;
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);
	string sequence = ""; 
	string sequence_Q = ""; 
	char caractere;

	IfFicFasta.seekg(seq->getPos_encodage()); 

	while(IfFicFasta.get(caractere) and caractere != '+') {
		if (caractere != '\n') {
			sequence+=caractere;
		}
	}

	IfFicFasta.seekg(seq->getPos_encodageQ()); 

	while(IfFicFasta.get(caractere) and longueur != seq->length()) {
		if (caractere != '\n') {
			sequence_Q+=caractere;
			longueur++; 
		}
	}

	IfFicFasta.close();
	return new Sequence_FastQ(seq->getIntitule(), sequence.c_str(), sequence_Q.c_str());
}

int next_sequence (string nomFicFasta, int position) {
	size_t longueur = 0;
	size_t i = 0; 
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	char caractere; 
	// leture de l'entête

	IfFicFasta.seekg(position); 

	IfFicFasta.ignore(numeric_limits<streamsize>::max(), '\n' );

	// lecture de la séquence 
	while(IfFicFasta.get(caractere) and caractere != '+') {
		if (caractere != '\n') {
			longueur+=1;
		}
	}

	IfFicFasta.ignore(numeric_limits<streamsize>::max(), '\n' );

	while(IfFicFasta.get(caractere) and i != longueur) {
		if (caractere != '\n') {
			i++;
		}
	}

	// verifie si il y a une sequence ensuite ou la fin du fichier
	IfFicFasta.peek() != '@' ? position = -1 : position = IfFicFasta.tellg(); 

	IfFicFasta.close();

	return position; 
			
}


vector <size_t> extractPosSeqQ (string nomFicFasta){
	vector <size_t> positions;
	int position = 0;
	
	// verifier que 1er est @
	positions.push_back(position+1);

	while(position != -1) {
		position = next_sequence(nomFicFasta, position);
		if (position != -1)
			positions.push_back(position+1);
		
	}

	return positions;
}

vector <Sequence_FastQ*> createAllQ (string nomFicFasta, vector <size_t> positions, size_t debut, size_t fin, bool isSeq,  vector <FakeSeq*> & seq){
	vector <Sequence_FastQ*> allseq; 
	size_t pos;
	for (size_t i = debut; i < fin; i++){
		if (debut)
			pos = debut-i;
		else 
			pos = i; 
		
		if (!isSeq) 
			allseq.push_back(createIntituleQ(nomFicFasta,positions[pos]));
		else 
			allseq.push_back(createSequenceQ(nomFicFasta,seq[pos]));
	}

	return allseq;
}

bool Sequence_FastQ::operator==(Sequence_FastQ const &s) {
	bool isEqual = (*this).Sequence_FastA::operator==(s); 
	
	if (isEqual) {
		for (size_t i = 0; i < longueur; i++) {
			if (this->encodage_Q[i] != s.encodage_Q[i]) {
				isEqual = false; 
			}
		}
	} else {
		isEqual = false; 
	}
	return isEqual; 
}

Sequence_FastQ &Sequence_FastQ::operator=(Sequence_FastQ const &s) {
	(*this).Sequence_FastA::operator=(s); 
	delete[] encodage_Q; 
	encodage_Q = NULL;
	if (longueur) {
		encodage_Q = new char (longueur); 
		for (size_t i = 0; i < longueur; i++){
			encodage_Q[i] = s.encodage_Q[i];
		}
	}
	return (*this);
}

// Sequence_FastA* concatAllSeq (vector <Sequence_FastQ*> sequences) {
// 	Sequence_FastA * s = sequences[0] ; 
// 	for (size_t i = 1; i < sequences.size(); i++) {
// 		s = concatSeq(s, sequences[i]);
// 		delete sequences[i];
// 	}
// 	return s; 
// }

void Sequence_FastQ::printQualite() const{
	if (encodage_Q != NULL) 
		for (size_t i=0; i < longueur; i++) {
			cout<<encodage_Q[i];
		}
}


char* Sequence_FastQ::getQualite() const{
	return encodage_Q;
}


Sequence_FastQ* concatAllSeqQ (vector <Sequence_FastQ*> sequences){
	//Sequence_FastA * s = concatAllSeq(sequences);
	Sequence_FastA * s = sequences[0] ; 

	for (size_t i = 1; i < sequences.size(); i++) {
		s = concatSeq(s, sequences[i]);
	}

	char * newQ = new char [s->length()];

	size_t pos = 0; 

	for (size_t i = 0; i < sequences.size(); i++ ){
		for (size_t j = 0; j < sequences[i]->length(); j++) {
			newQ[pos] = sequences[i]->getQualite()[j];
			pos++;
		}
	}

	Sequence_FastQ* sQ= new Sequence_FastQ(s, newQ);

	delete[] newQ;
	delete s; 
	return  sQ;
}
