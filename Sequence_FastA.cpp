#include "Sequence_FastA.h"

Sequence_FastA::Sequence_FastA(): Sequence_FastX(),encodage(NULL){}

Sequence_FastA::Sequence_FastA(string titre, size_t lg): Sequence_FastX(titre,lg),encodage(NULL){}

Sequence_FastA::Sequence_FastA(string titre, size_t lg, unsigned char * seq):Sequence_FastX(titre,lg),encodage(seq){}

Sequence_FastA::Sequence_FastA(const Sequence_FastA* const &s)
:Sequence_FastX(s),encodage(new unsigned char[longueur/4 + (longueur%4 != 0)]){

	for (size_t i= 0; i < longueur/4 + (longueur%4 != 0); i++){
		encodage[i] = s->encodage[i];
	}
}

Sequence_FastA::Sequence_FastA(Sequence_FastA const &s)
:Sequence_FastX(s),encodage(new unsigned char[longueur/4 + (longueur%4 != 0)]){

	for (size_t i= 0; i < longueur/4 + (longueur%4 != 0); i++){
		encodage[i] = s.encodage[i];
	}
}

Sequence_FastA::~Sequence_FastA(){
	if (encodage != NULL) 
		delete[] encodage;
}

unsigned char * Sequence_FastA::getEncodage() const {
	return encodage;
}

// idée, donner le flux et la position
Sequence_FastA::Sequence_FastA(string titre, const char * const seq):Sequence_FastX(titre, 0),encodage(NULL)
{
	if (seq != NULL)
	{	
		// Recuperation de la longueur de la sequence et d'informations sur les caracteres invalides 
		bool erreur = false; 
		while(seq[longueur] != '\0')
		{

			if (!estValide(seq[longueur]))
			{
				erreur = true;
			}
			
			longueur++;
		}

		if (erreur) {
				//cout<<"Il y a un ou plusieurs caractere-s invalide-s!"<<"Transformation en A."<<endl; 
		}

		size_t longueur_tableau = longueur/4 + (longueur%4 != 0); 

		encodage = new unsigned char[longueur_tableau];

		// initialisation du tableau a 0 
		for (size_t i = 0; i < longueur_tableau; i++) {
			encodage[i] = 0; 
		}

		// remplissage du tableau 
		unsigned char tmp = 0;

		for (size_t i = 0; i < longueur; i++) {
			// tant que l'octet n'est pas rempli 
			if (!(donnePosOctet(i))) {
				tmp = tmp << 2;
				tmp = tmp | encode(seq[i]); 
				// ajout dans le tableau et reinitilisation de tmp 
				encodage[donneOctet(i, longueur_tableau)] = tmp; 
				tmp = 0; 
			}
			tmp = tmp << 2;
			tmp = tmp | encode(seq[i]); 
		}
		// ajout du dernier octet avec decalage des valeurs vers la droite pour obtenir 000.. vers les bits de poids faible
		// quand l'octet n'est pas rempli 
		if (longueur%4) 
			encodage [donneOctet(longueur-1, longueur_tableau)] = tmp << ((4-(longueur%4)) * 2);

		

	}
}

Sequence_FastA* Sequence_FastA::complement() const {
	size_t longueur_tableau = longueur/4 + (longueur%4 != 0); 
	unsigned char * sComplement = new unsigned char[longueur_tableau];
	unsigned char tmp; 

	for (size_t i=0; i<longueur_tableau; i++){
		tmp = encodage[i];
		// pour chaque case du tableau, bit complementaires 
		sComplement[i] = ~tmp; 
	}
	return new Sequence_FastA(intitule, longueur, sComplement);
}


Sequence_FastA* Sequence_FastA::reverse() const {
	size_t longueur_tableau = longueur/4 + (longueur%4 != 0); 
	unsigned char * reverse = new unsigned char[longueur_tableau];
	unsigned char  * treverse =  new unsigned char[longueur_tableau];
	unsigned char tmp = 0;
	unsigned char tempo = 0; 
	size_t cpt = 0;

	for (size_t i=0; i<longueur_tableau; i++) {
		tempo =  encodage[longueur_tableau-i-1]; 

		// jusqu'à qu'on ait rempli un octet
		while (cpt < 4) {
			cpt++; 
			tmp <<= 2;
			// recupérée puis décalée vers la gauche, pour inverser
			tmp = tmp | (tempo&3);
			tempo >>= 2;
		}
		cpt=0;
		reverse[i] = tmp;
		tmp = 0;
	}

	if (longueur%4){
		treverse = remplissageOctet(reverse, (longueur%4)*2-2, 0, longueur);
		return new Sequence_FastA(intitule, longueur, treverse);
	}
	else 
		return new Sequence_FastA(intitule, longueur, reverse);
}

// Extremite incluses
Sequence_FastA* Sequence_FastA::subSequence(int begin, int end) const {
	if (begin < 0 || begin >= int(longueur)) {
		begin = 0;
	}

	if (end >= int(longueur)) {
		end = longueur - 1;
	}

	size_t longueur_tableau = longueur/4 + (longueur%4 != 0); 

	// pour savoir où se positionner dans la Sequence_FastA encodee
	size_t octet = donneOctet(begin, longueur_tableau);
	int pos_octet = donnePosOctet(begin);

	return new Sequence_FastA(intitule, end-begin+1, remplissageOctet(encodage, pos_octet, octet, end-begin+1));
}

string Sequence_FastA::operator[](size_t i) const{
	return decode(donneNext(donneOctet(i, longueur/4 + (longueur%4 != 0)), donnePosOctet(i), encodage));
}

void Sequence_FastA::setSymbolAt(size_t i, char c){
	int pos_octet = donnePosOctet(i); 
	size_t octet = donneOctet(i, longueur/4 + (longueur%4 != 0));

	unsigned char new_symbol = encode(c);
	// position maximale octet
	int j = 6;
	unsigned char tmp = 0;

	// tant qu'on a pas parcouru tout l'octet
	while (j != -2) {
		tmp <<= 2;
		// si on atteint les 2 bits qu'on veut changer
		if (j == pos_octet){
			tmp = tmp | new_symbol;
			j-=2;
		}
		else { // sinon on recopie ceux qui étaient déjà présents
			tmp = tmp | donneNext(octet, j, encodage);
			j-=2;
		}
	}
	encodage[octet] = tmp;

}

void Sequence_FastA::print() const{

	Sequence_FastX::print();

	if (encodage != NULL) {
		cout<<endl<<"Sequence : "<<endl;	
		for (size_t i = 0; i < longueur; i++){
			cout<<(*this)[i];
		}
		cout<<endl;
	}
}	
