#include "FakeSeq.h"

FakeSeq::FakeSeq(string titre, size_t lg, size_t seq, size_t seq_Q ): Sequence_FastX(titre,lg),pos_encodage(seq),pos_encodageQ(seq_Q){}

FakeSeq::FakeSeq(const FakeSeq* const &s):Sequence_FastX(s), pos_encodage(s->pos_encodage), pos_encodageQ(s->pos_encodageQ){}

FakeSeq::~FakeSeq(){}

size_t FakeSeq::getPos_encodage() const {
	return pos_encodage;
}
size_t FakeSeq::getPos_encodageQ() const {
	return pos_encodageQ;
}

void FakeSeq::print() const{
	Sequence_FastX::print(); 

	cout<<"Début encodage à : "<<pos_encodage<<endl;
	cout<<"Début encodage qualité à : "<<pos_encodageQ<<endl;
}

bool FakeSeq::operator==(FakeSeq const &s) {
	bool isEqual = (*this).Sequence_FastX::operator==(s); 
	return (isEqual && (this->pos_encodage == s.pos_encodage) && (this->pos_encodageQ == s.pos_encodageQ)) ;
}


FakeSeq* createFakeSeqQ (string nomFicFasta, size_t position) {
	size_t pos_encodage, pos_encodageQ;
	string intitule = ""; 
	size_t longueur = 0;
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	char caractere;
	// leture de l'entête 
	vector <string> analyse = analyseEntete(nomFicFasta, position, ':');

	IfFicFasta.seekg(position); 

	// lecture intitule
	while(IfFicFasta.get(caractere) and caractere != '\n') {
		intitule += caractere;
	}

	// debut de la sequence
	pos_encodage = IfFicFasta.tellg(); 

	// lecture de la séquence 
	while(IfFicFasta.get(caractere) and caractere != '+') {
		if (caractere != '\n') {
			longueur+=1;
		}
	}

	//lecture de la deuxieme entête 

	bool isEqual = compareEntete(analyse, analyseEntete(nomFicFasta, IfFicFasta.tellg(), ':'));

	if (!isEqual) {
		cout<<"Attention les deux entêtes ne sont pas identiques!"<<endl;
	}

	IfFicFasta.ignore(numeric_limits<streamsize>::max(), '\n' ); 

	// debut de la sequence qualite
	pos_encodageQ = IfFicFasta.tellg();

	IfFicFasta.close();

	return new FakeSeq(intitule, longueur, pos_encodage, pos_encodageQ);
}

FakeSeq* createFakeSeqA (string nomFicFasta, size_t position){
	string intitule = ""; 
	size_t longueur = 0;
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	char caractere;

	vector <string> analyse = analyseEntete(nomFicFasta, position-1, '|');
	intitule = analyse.back();

	IfFicFasta.seekg(position); 

	IfFicFasta.ignore(numeric_limits<streamsize>::max(), '\n' ); 

	position = IfFicFasta.tellg();

	while(IfFicFasta.get(caractere) and caractere != '>') {
		if (caractere != '\n') {
			longueur+=1;
		}
	}

	IfFicFasta.close();

	return new FakeSeq(intitule, longueur, position);
}


vector <FakeSeq*> createAllFakeSeq (string nomFicFasta, vector <size_t> positions, size_t debut, size_t fin, bool isSeqQ){
	vector <FakeSeq*> allseq; 
	for (size_t i = debut; i < fin; i++){
		if (!isSeqQ) 
			allseq.push_back(createFakeSeqA(nomFicFasta,positions[i]));
		else 
			allseq.push_back(createFakeSeqQ(nomFicFasta,positions[i]));
	}
	return allseq;
}

FakeSeq &FakeSeq::operator=(FakeSeq const& s) {
	intitule = s.intitule;
	longueur = s.longueur;
	pos_encodage = s.pos_encodage;
	pos_encodageQ = s.pos_encodageQ;
	return (*this);
}


