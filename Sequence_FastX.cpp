#include "Sequence_FastX.h"


Sequence_FastX::Sequence_FastX():intitule(),longueur(0) {}

// Destructeur
Sequence_FastX::~Sequence_FastX(){}

Sequence_FastX::Sequence_FastX(string titre, size_t lg):intitule(titre),longueur(lg){}

// Constructeur par copie
Sequence_FastX::Sequence_FastX(const Sequence_FastX* const &s)
:intitule(s->intitule),longueur(s->longueur){}

string Sequence_FastX::getIntitule() const {
	return intitule;
}

void Sequence_FastX::setIntitule(string s){
	intitule =  s;
}

size_t Sequence_FastX::length() const {
	return longueur;
}

void Sequence_FastX::print() const{

	cout<<endl<<"Titre : "<<intitule<<endl;
	cout<<endl<<"Longueur : "<<longueur<<endl;	
}	


void printInfoEntente (vector <string> analyse) {

	if (analyse.size() == 3 && analyse[0] == "ref") {
		cout<<"Numero accession : "<<analyse[1]<<endl;
		cout<<"Locus : "<<analyse[2]<<endl;	
	} else {
		if (analyse.size() == 5 && analyse[0] == "gi") {
			cout<<"gi-number : "<<analyse[1]<<endl;
			cout<<"Numero accession : "<<analyse[3]<<endl;
			cout<<"Locus : "<<analyse[4]<<endl;	
		}
		else {
			cout<<"Voici les informations concernant l'entête:"<<endl;
			for (size_t i = 0; i < analyse.size(); i++) {
				cout<<"\t --> "<<analyse[i]<<endl;
			}
		}
	}

}

// separateur peut etre '|' ou ':' 
vector <string> analyseEntete (string nomFicFasta, size_t position, char separator) {
	string intitule = ""; 
	vector <string> analyse; 
	char caractere;
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);
	// Se place a la position du debut de l'entete a analyser
	IfFicFasta.seekg(position); 

	// recupere chaque partie de l'entete separee par le separateur
	while(IfFicFasta.get(caractere) and caractere != '\n') {
		if (caractere == separator ){
			analyse.push_back(intitule);
			intitule = "";
		}
		else {
			intitule += caractere;
		}
	}
	analyse.push_back(intitule);
	IfFicFasta.close(); 

	return analyse;

}

// return true si sequence FastA 
bool isFastA (string nomFicFasta) {
	int Qpos = 1, Apos = -1;
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	// Positions du premier @ 
	IfFicFasta.ignore(numeric_limits<streamsize>::max(), '@' ); 
	Qpos = IfFicFasta.tellg(); 

	// Et du premier >
	IfFicFasta.ignore(numeric_limits<streamsize>::max(), '>' ); 
	Apos = IfFicFasta.tellg(); 
	
	return Qpos < Apos ? false : true; 

	IfFicFasta.close();

}

// retunr true si le fichier a pu etre ouvert 
bool checkFile (string nomFicFasta) {
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	if (!IfFicFasta) {
		cout<<endl<<"Impossible d'ouvrir le fichier: "<<nomFicFasta<<endl<<endl;
		return 0;
	}

	IfFicFasta.close();
	return 1; 
}


bool compareEntete (vector <string> a1, vector <string> a2) {
	if (a1.size() == a2.size()) {
		for (size_t i = 0; i < a1.size(); i++) {
			if (a1[i] != a2[i]) {
				return false;
			}
		}
	}
	return true;
}

bool Sequence_FastX::operator==(Sequence_FastX const &s) {
 	return (this->longueur == s.longueur) ? true : false; 

}

