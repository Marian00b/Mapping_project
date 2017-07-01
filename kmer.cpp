#include "kmer.h"

using namespace std;

Kmer::Kmer(){
	longueur = 0; 
	kmers = vector <Sequence_FastA*> ();
	exact_match = map <size_t, vector<int>> (); 
	non_match = map <size_t, vector<int>> (); 
	isLocalise_exact = false;
	isLocalise_NonExact = false;
	localisation = 0;
}


Kmer::~Kmer() {
	if (kmers.size()) {
		for (size_t i = 0; i < kmers.size(); i++){
			delete kmers[i];
		}
	}
}

int Kmer::getLocalisation(){
	return localisation;
}

bool Kmer::getIsLocalise_exact (){
	return isLocalise_exact ;
}

bool Kmer::getIsLocalise_NonExact (){
	return isLocalise_NonExact ;
}


void Kmer::setMatch(size_t i, int v){
	if (exact_match.find(i) == exact_match.end())
		exact_match[i] = vector<int> (); 

	exact_match[i].push_back(v);
}

void Kmer::setNewNonMatch(size_t i, vector<int> v){
	non_match[i] = v; 
}


void Kmer::setNonMatch(size_t i, int v){
	if (non_match.find(i) == non_match.end())
		non_match[i] = vector<int> (); 

	non_match[i].push_back(v);
}


vector<int> & Kmer::getPositions(size_t i){
	return exact_match[i].size() ? exact_match[i] : non_match[i];
}

// est deja present
bool Kmer::isKmer (Sequence_FastA* const& seq) const {
	size_t i = 0; 
	while (i < kmers.size()) {
		if ((*seq) == (*kmers[i])) {
			return true;
		}
		i++;
	}
	return false;
}

// longueur du kmer
size_t Kmer::length() const {
	return longueur;
}

// nombre de kmer
size_t Kmer::nombre() const{
	return kmers.size();
}

Kmer::Kmer(Sequence_FastA* const& seq, size_t lg): longueur (lg){
	localisation = 0;
	isLocalise_NonExact = false;
	isLocalise_exact  = false;
	non_match = map <size_t, vector<int>> (); 
	exact_match = map <size_t, vector<int>> (); 
	kmers = vector <Sequence_FastA*> ();
	Sequence_FastA* k;
	size_t cpt = 0;
	string intitule = seq->getIntitule();

	if (seq->length() >= longueur ) {
		for (size_t i = 0; i < seq->length()-longueur+1; i++) {
			k = seq->subSequence(i, i+longueur-1);

				k->setIntitule(intitule + to_string(i) );
			if (!isKmer(k)) {
				kmers.push_back(k);
				cpt++;
			}
		}
	}
	else {
		//cout<<"Attention un Kmer vide!"<<endl;
		kmers.push_back(seq);
	}
}

void Kmer::print(){
	cout<<"Kmer de longueur "<<longueur<<" :"<<endl;

	if (kmers.size() > 0) {
		for (size_t i = 0; i < kmers.size(); i++) {
			cout<<kmers[i]->getIntitule()<<" "<<(*kmers[i])<<" "<<endl;
		}
	}
	else 
		cout<<"Aucun Kmers!";
	cout<<endl;
} 

vector<Kmer*>  createAllKmers (vector<Sequence_FastQ*> const sequences, size_t longueur) {
	vector<Kmer*> kmers; 
	string intitule = ""; 
	for (size_t i = 0; i < sequences.size(); i ++){
		intitule = "Read:" + to_string(i) + ":";
		sequences[i]->setIntitule(intitule+sequences[i]->getIntitule());
		kmers.push_back(new Kmer(sequences[i], longueur)); 
	}

	return kmers;
}

Sequence_FastA* const& Kmer::getKmer(size_t i ){
	return kmers[i];
}



	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	********************************************* ANALYSE KMER **************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

// est sur le meme brin 
bool sameBrin (int nb1, int nb2){
	return ((nb1 < 0 && nb2 < 0) || (nb1 >= 0 && nb2 >= 0)) ? true : false; 
}

// inverse le tableau en fonction de la position actuelle
void find_closer(int pos, vector<int> & v){
	int closer = 0;
	if (v.size()) {
		for (size_t i = 1; i < v.size(); i++){
			if ( abs((v[i] - pos)) < abs((v[closer] - pos))) {
				closer = i;
			}
		}

		pos = v[closer];

		v.erase (v.begin()+closer);

		v.insert(v.begin(),pos);
	}

}


// peut-être que le const ne va pas fonctionner
int Kmer::isSuite (map <int, size_t> const& res, int ecart){
	int suite = 0;
	bool isSuite; 
	int localisation_max;
	int last_of_suite; // va contenir la position qu'on a retenu 
	size_t k,j;
	int max = 0;
	size_t gap = 0;

	for (size_t i = 0; i < nombre() -1 ; i++) {
		j = 0;
		isSuite = false;

		// comparaison des positions jusqu'a trouver une suite de position
		while (!isSuite && j < getPositions(i).size()){
			k = 0;
				
			// on met la liste des position dans l'ordre
			if (getPositions(i+1).size() > 1)
				find_closer(getPositions(i)[j], getPositions(i+1));

			// avec les autres positions jusqu'a trouver une suite de position 
			while (!isSuite && k < getPositions(i+1).size()){

				// same brin et ecart entre 0 et ecart 
				if (sameBrin(getPositions(i+1)[k],getPositions(i)[j]) && 
					(getPositions(i+1))[k] - (getPositions(i))[j] <= ecart && (getPositions(i+1))[k] - (getPositions(i))[j] >= 0) {

						isSuite = true;

					// si c'est pas 0 ou qu'il y a un precedent, si l'actuel est le suivant du precedent
					if (!i || (getPositions(i-1).size() && gap == 0 &&  getPositions(i)[j] == last_of_suite)){
						
						// max si debut de suite
						if (suite == 0) {
							localisation_max = (getPositions(i))[j];
						}
						suite++;
						gap = 0;

						// memorisation du precedent
						last_of_suite  = (getPositions(i+1))[k];
					}
					else {
						//cas d'un gap entre deux kmer eloignes, gap est incremente quand on sort des 2 while
						if ( getPositions(i-gap).size() && sameBrin(getPositions(i-gap)[k],getPositions(i)[j])  
							&& ((getPositions(i))[j] - (getPositions(i-gap))[k]<= (ecart+gap)) ){
							if (suite == 0) {
								localisation_max = (getPositions(i-gap))[k];
							}
							suite+=2;
							gap = 0;
							last_of_suite  = (getPositions(i+1))[k];
						}
						else {
							// recherche du max quand on termine une suite
							if (suite > max) {
								max = suite; 
								localisation = localisation_max;

							}
							last_of_suite  = (getPositions(i+1))[k];	
							localisation_max = (getPositions(i))[j];
							gap = 0;
							suite = 0;
							suite++;
						}
					}
					
				}
				k++;
			}
			j++;

		}
		if (suite > max) {
			max = suite; 
			localisation = localisation_max;
		}
		if (!isSuite) {
			gap++;
		}	
	}

	// - 1 si un seul kmer ? suite 1 2 = 1 mais suite 1 = 0; 
	if (max >= (nombre()) * 0.8) {
		isLocalise_exact  = true;
	}
	else {
		if (max >= (nombre()) * 0.3) 
			isLocalise_NonExact  = true;
	}
	return max;
}


// map< Sequence_FastA*, int> getNonLocalise(vector<Kmer*> Reads) {
// 	vector <size_t> ordre;
// 	map< Sequence_FastA*, int> non_localise; 

// 	for (size_t i = 0; i < Reads.size(); i ++){
// 		for (size_t j = 0; j < Reads[i]->nombre(); j++){

// 			if (!(Reads[i]->getPositions(j)).size()){

// 				if (non_localise.find(Reads[i]->getKmer(j)) != non_localise.end()){
// 					non_localise[Reads[i]->getKmer(j)]= non_localise[Reads[i]->getKmer(j)] + 1;
// 				}
// 				else 
// 					non_localise[Reads[i]->getKmer(j)]=1;


	
// 			}



// 		}	
// 	}

// 	return non_localise;

// }

