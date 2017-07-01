#include "Sequence_FastA.h"


bool estValide(char nuc)
{
	return (nuc=='A') || (nuc=='a')
		|| (nuc=='C') || (nuc=='c')
		|| (nuc=='G') || (nuc=='g')
		|| (nuc=='T') || (nuc=='t');
}

unsigned char encode(char nuc) 
{
	unsigned char num; 
	switch (nuc) {
		case 'c' : 
		case 'C' :
			num = 1; 
			break;
		case 'g' : 
		case 'G' : 
			num = 2;
			break;
		case 't': 
		case 'T': 
			num = 3;
			break;
		default : 
			num = 0; 

	}
	return num; 
}

size_t donneOctet (size_t pos, size_t lg) {

	return pos/4 < lg ? pos/4 : -1; 
}

size_t donnePosOctet (size_t pos) {
	return (3-(pos%4)) * 2; 
}

string decode(char num){
	string nuc; 
	switch (num) {
		case 1 :
			nuc = "C"; 
			break;
		case 2 : 
			nuc = "G";
			break;
		case 3 : 
			nuc = "T";
			break;
		case 0 : 
			nuc = "A";
	}
	return nuc;
}

unsigned char donneNext (size_t octet, size_t position, unsigned char * const& encodage) {
	// contient l'octet courant
	unsigned char tempo = encodage[octet]; 
	// va contenir les deux bits de l'octet en une position donnee
	unsigned char tmp = 0;

	position = position/2;
	size_t cpt = 0;

	// tant qu'on a pas atteint la position on decale de 2
	while (cpt <= position) {
		tmp = 0; 
		cpt++;
		// &3 pour le %4 -> donne les deux derniers bits
		tmp = tmp | (tempo&3); 
		tempo >>= 2;
	}
	
	return tmp; 
}


// Fonction qui va permettre de copier un tableau d'encodage à partir d'un octet et d'une position d'un octet
unsigned char * remplissageOctet(unsigned char * const& encodage, int pos_octet, size_t octet, int lg){
		size_t longueur_new = lg/4 + (lg%4 != 0);
		unsigned char * new_tab = new unsigned char [longueur_new];
		// compteur pour savoir où se positionner dans le tableau
		size_t i_tab = 0; 
		// compteur pour savoir quand un octet est complet 
		size_t cpt = 0; 
		// compteur pour savoir quand on ajoute toutes les lettres
		int i = 0;
		// va contenir l'octet a inserer, construit avec donneNext
		unsigned char next = 0; 

		// tant qu'on a pas dépassé la taille de l'encodage a copier 
		while (i <= lg) {
			// lorsqu'on a rempli un octet on l'ajoute dans le nouveau tableau
			if (cpt == 4) {
				new_tab[i_tab] = next;
				cpt = 0;
				next = 0; 
				i_tab++;
			} 
			else {
				// Pour le dernier octet 
				if (i== lg) {
					//int j = 0 ;
					while (cpt != 4){
					//while (j <= pos_octet/2) {
						next <<=2;
						//j++;
						cpt++;
					} 
					new_tab[i_tab] = next;
				}
			}  
			// lorsque l'on a terminé de parcourir un octet de l'encodage a copier, on change d'octet et de position
			// -2 car a la fin on fait -2 systématiquement
			if (pos_octet == -2) {
				pos_octet = 6;
				octet+=1;
			}
			next <<= 2; 
			// on va chercher les deux prochains bits de l'octet
			next = next | donneNext(octet, pos_octet, encodage);
			cpt++; 
			i++;
			pos_octet-=2;
		}
		return new_tab;
}



Sequence_FastA* createIntitule (string nomFicFasta, size_t position, char separator){
	string intitule = ""; 
	size_t longueur = 0 ; 
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	char caractere;

	vector <string> analyse = analyseEntete(nomFicFasta, position, separator);

	char end = '>';

	IfFicFasta.seekg(position);  

	while(IfFicFasta.get(caractere) and caractere != '\n') {
		position = IfFicFasta.tellg();
		if ((caractere != '>' || caractere != '@') && (separator == ':')) {
			intitule += caractere;
		}
	}

	if (separator == ':') {
		 end = '+';
	} else {
		intitule = analyse.back();
	}

	IfFicFasta.clear();
	IfFicFasta.seekg(position);

	while(IfFicFasta.get(caractere) and caractere != end) {
		if (caractere != '\n') {
			longueur+=1; 
		}
	}

	IfFicFasta.close();
	return new Sequence_FastA(intitule, longueur);
}

Sequence_FastA* createSequence (string nomFicFasta, FakeSeq * seq){
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);
	string sequence = ""; 
	char caractere;

	IfFicFasta.seekg(seq->getPos_encodage()); 

	while(IfFicFasta.get(caractere) and caractere != '>') {
		if (caractere != '\n') {
			sequence+=caractere;
		}
	}

	IfFicFasta.close();
	return new Sequence_FastA(seq->getIntitule(), sequence.c_str());
}

vector <size_t> extractPosSeq (string nomFicFasta){
	vector <size_t> positions;
	ifstream IfFicFasta (nomFicFasta.c_str(), ios::in);

	char caractere; 

	while(IfFicFasta.get(caractere)) {
		if (caractere == '>') {
			positions.push_back(IfFicFasta.tellg());
		} 	
	}
	IfFicFasta.close();
	return positions;
}

vector <Sequence_FastA*> createAll (string nomFicFasta, vector <size_t> positions, size_t debut, size_t fin, bool isSeq, vector <FakeSeq*> seq){
	vector <Sequence_FastA*> allseq; 
	size_t pos;
	
	for (size_t i = debut; i < fin; i++){
		if (debut)
			pos = debut-i;
		else 
			pos = i; 
		
		if (!isSeq) 
			allseq.push_back(createIntitule(nomFicFasta,positions[pos], '|'));
		else 
			allseq.push_back(createSequence(nomFicFasta,seq[pos]));
	}

	return allseq;
}

 ostream& operator<<(ostream &os, Sequence_FastA const& es) {
	for (size_t i = 0; i< es.length(); i++) {
		os<<es[i];
	}
	return os; 
}

bool Sequence_FastA::operator==(Sequence_FastA const &s) {
	bool isEqual = (*this).Sequence_FastX::operator==(s); 
	if (isEqual){
		for (size_t i = 0; i < longueur/4 + (longueur%4 != 0); i++) {
			if (this->encodage[i] != s.encodage[i]) {
				isEqual = false; 
			}
		}
	}
	
	return isEqual; 
}

// Sequence_FastA &Sequence_FastA::operator=(Sequence_FastA* const &s){
// 	cout<<"Test"<<endl;
// 	intitule = s->intitule;
// 	size_t longueur_tableau = longueur/4+((longueur&3)>0); 
// 	if (longueur != s->longueur){
// 		if (longueur) {
// 			delete[] encodage; 
// 			encodage = NULL;
// 		}
// 		longueur = s->longueur; 
// 		if (longueur) {
// 			encodage = new unsigned char (longueur_tableau); 
// 		}
// 	}
// 	if (longueur) {
// 		for (size_t i = 0; i < longueur_tableau; i++){
// 			encodage[i] = s->encodage[i];
// 		}
// 	}
// 	return (*this);	
// }




Sequence_FastA &Sequence_FastA::operator=(Sequence_FastA const &s){
	cout<<"Test"<<endl;
	intitule = s.intitule;
	size_t longueur_tableau = longueur/4+((longueur&3)>0); 
	if (longueur != s.longueur){
		if (longueur) {
			delete[] encodage; 
			encodage = NULL;
		}
		longueur = s.longueur; 
		if (longueur) {
			encodage = new unsigned char (longueur_tableau); 
		}
	}
	if (longueur) {
		for (size_t i = 0; i < longueur_tableau; i++){
			encodage[i] = s.encodage[i];
		}
	}
	return (*this);	
}


Sequence_FastA* concatSeq (Sequence_FastA* s, Sequence_FastA* iseq){
	size_t longueur = s->length() + iseq->length(); 
	size_t lg = s->length()/4 + (s->length()%4 != 0);

	unsigned char * big_seq = new unsigned char [longueur/4 + (longueur%4 != 0)]; 


	for (size_t j = 0; j < lg - 1 ; j++) {
		big_seq[j] = s->getEncodage()[j];
	}

	unsigned char tmp = s->getEncodage()[lg-1]; 
	unsigned char iTmp = 0; 
	size_t pos = 6;

	// construction de l'octet intermediaire 
	if (s->length() % 4) {
		for (size_t i = 0; i < (4 - (s->length()%4)) ; i++){
			iTmp = iTmp<<2;
			iTmp = iTmp | donneNext(0, pos,iseq->getEncodage());
			pos -=2; 
		}

		tmp = tmp | iTmp;
	}

	big_seq[lg-1] = tmp;

	size_t longueur_new; 

	// ce qui a deja ete pris pour completer l'octet
	if (iseq->length() > (4-(s->length()%4)) && s->length() % 4)
		longueur_new = iseq->length()-(4-(s->length()%4)); 
	else 
		longueur_new = iseq->length();

	//unsigned char * test = new unsigned char [longueur_new];

	unsigned char * test = remplissageOctet(iseq->getEncodage(), pos, 0, longueur_new);

	pos = 0; 

	for (size_t j = lg; j < longueur/4 + (longueur%4 != 0) ; j++){
		big_seq[j] = test[pos];
		pos++; 
	}

	delete[] test ;
	return new Sequence_FastA("BIG_SEQ", longueur, big_seq) ;  
} 


Sequence_FastA* concatAllSeq (vector <Sequence_FastA*> sequences){
	Sequence_FastA * s = sequences[0]; 
	for (size_t i = 1; i < sequences.size(); i++) {
		s = concatSeq(s, sequences[i]);
	}
	return s; 
}

