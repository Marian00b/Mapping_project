#include "SA.h"
#include <algorithm> 
#include <numeric>
#include <sstream>
#include <cmath>

using namespace std; 

vector<size_t> const& SAcomp::getSA() const {
	return SA;
}

struct SAcomp::SATest {
	Sequence_FastA* const& seqA;

	SATest(Sequence_FastA* const& seq): seqA(seq){}

	bool operator()(size_t A,size_t B)const{
		bool aInfB = false; 

		while (A < seqA->length() && B < seqA->length() && (*seqA)[A] == (*seqA)[B]) {
			A++;
			B++;
		}

		if (B != seqA->length()) {
			if (A == seqA->length()) 
					aInfB = true; 
			else {
				if ((*seqA)[A] < (*seqA)[B]) 
					aInfB = true;
			}
		}
		return aInfB;
	}
};

SAcomp::SAcomp(Sequence_FastA* const& seq): seqA (seq) {

	SA = vector <size_t>();
	SA.resize(seqA->length());
	iota(SA.begin(), SA.end(), 0);
	make_heap(SA.begin(), SA.end(), SATest(seqA));
	sort_heap (SA.begin(), SA.end(), SATest(seqA));
	buildLCP();
}


void SAcomp::print() const{
	for (size_t i = 0; i < SA.size(); i++) {
		cout<<SA[i]<<" ";
	}
	cout<<endl;
}


void SAcomp::buildLCP(){
	for (size_t i = 0; i < SA.size()-1; i++){
		LCP.push_back(prefcommSuff(SA[i], SA[i+1]));
	}
}

// Entre deux suffixes d'une meme sequence.
size_t SAcomp::prefcommSuff(size_t i, size_t j){
	size_t cpt = 0;
	while(i < seqA->length() && j < seqA->length() && (*seqA)[i] == (*seqA)[j]) {
		cpt++;
		i++;
		j++;
	} 
	return cpt;
}


Sequence_FastA* SAcomp::getIseq(size_t i){
	Sequence_FastA* s ;
	if (i < SA.size())
		s = seqA->subSequence(SA[i], seqA->length()-1);
	else 
		s = seqA->subSequence(SA[0], seqA->length()-1) ;
	return s; 
}

// connaitre la longueur de la sequence
size_t SAcomp::length(){
	return SA.size();
}

// connaitre la longueur du suffixe i de la SA
size_t SAcomp::iLength(size_t i) {
	if (i < SA.size()) 
		return length() - SA[i];
	else 
		return length();
}


	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	******************************************* AUTRE ***********************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

Sequence_FastA* SAcomp::getMot (size_t k, size_t pos){
	if (pos < seqA->length()) {
		size_t longueur = seqA->length() - pos; 
		size_t i = 0; 

		while (iLength(i) != longueur) {
			i++;
		}

		if (SA[i]+k-1 > seqA->length()-1) {
			k = 1;
		}
		return seqA->subSequence(SA[i], SA[i]+k-1); 
	}
	else 
		return new Sequence_FastA();
}




// recuperer le ieme facteur de longueur k de la sequence
// un facteur est un prefixe d'un suffixe
// i dans le sens alphabetique
Sequence_FastA* SAcomp::getIfacteur(size_t pos, size_t lg) {
	//parcour de la SA
	for (size_t i = 0; i < SA.size(); i++){
		if (iLength(i) >= lg){
			if (!pos) {
				if (iLength(i)>=lg)
					return getIseq(i)->subSequence(0, lg-1);
				else 
					return getIseq(i);
			}
			else {
				pos--;
			}
		}
	}
	return new Sequence_FastA();
}

// Recuperer le facteur de longueur k qui a le plus d'occurrences dans la sequence
// Recuperation du 1er facteur seulement si plusieurs.
Sequence_FastA* SAcomp::getOccurences(size_t lg) {
	size_t cpt = 0; 
	Sequence_FastA * s; 
	size_t i = 0; 
	size_t max = 0; 

	while (i < SA.size()){
		if (iLength(i) >= lg){ 
			cpt = 1; 
			while (LCP[i] >= lg) {
				i++;
				cpt++;
			}
			if (cpt > max){
				max = cpt;
				s = getIseq(i-cpt+1)->subSequence(0, lg-1);
			}
		}
		if (cpt <=1)
			i++;
	}

	if (max != 0)
		return s; 
	else 
		return new Sequence_FastA();
}	

size_t SAcomp::getDistincts(size_t lg) {
	size_t i = 0; 
	size_t nb = 0; 

	if (lg) {
		while (i < SA.size()){
			if (iLength(i) >= lg){ 
				nb +=1; 
				while (i < SA.size()-1 && LCP[i] >= lg) {
					i++;
				}
				i++;
			}
			else 
				i++;
		}
	}
	return nb;

}	


	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	******************************************* STRNG ***********************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/
string char2string (char c){
	stringstream ss;
	string s;
	ss << c;
	ss >> s;
	return s;
}


//Entre un suffixe d'une sequence et une string
size_t Sprefcomm(size_t i, Sequence_FastA* const& seq, string W){
	stringstream ss;
	size_t cpt = 0;
	size_t j  = 0; 
	while(i < seq->length() && j < W.length() && (*seq)[i] == char2string(W[j])) {
		cpt++;
		i++;
		j++;
	} 
	return cpt;
}


//Entre un suffixe d'une sequence et une string
size_t SAcomp::SprefcommSA(size_t i, string W){
	stringstream ss;
	size_t j  = 0; 
	while(i < seqA->length() && j < W.length() && (*seqA)[i] == char2string(W[j])) {
		i++;
		j++;
	} 
	return j;
}




int SAcomp::dichoM (string W) {
	int d = -1;
	int f = SA.size();
	int m, l;

	Sequence_FastA* dm;
	
	while (d+1<f) {
		m = (d+f)/2;
		l = SprefcommSA(SA[m], W);
		
		dm = getIseq(m);
		if ( l == W.length()) {
			if (l == dm->length()) {
				return m;
			}
			else {
				f = m;
			}
		}
		else {
			if (l < dm->length() && l < W.length() && char2string(W[l]) < (*dm)[l]) {
				f = m;
			}
			else {
				d = m;
			}
		}	
	}

	if (m == SA.size()-1){
		delete dm;
		return -1;
	}
	else {
		if (!m)
			return m;
		else {
			if (SprefcommSA(SA[d], W) >= SprefcommSA(SA[f], W)){
				delete dm;
				return d;
			}
			else {
				delete dm;
				return f;
			}
		}
	}
}


// Rechercher toutes les occurences d'un facteur

vector<int> SAcomp::rechercheAllM (string W) {

	int pos_dicho = dichoM(W); 
	vector<int> positions; 

	if (pos_dicho == -1) 
		return positions;

	size_t lg = SprefcommSA(SA[pos_dicho], W);

	if (lg == W.length()) {
		positions.push_back(SA[pos_dicho]);

		while (pos_dicho < LCP.size() && LCP[pos_dicho] >= lg) {
			pos_dicho++;
			positions.push_back(SA[pos_dicho]);
		}
	}

	return positions;
}

	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	******************************************* ENCDE ***********************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/



// Entre deux sequences
size_t prefcomm(size_t i, Sequence_FastA* const& seq1, Sequence_FastA* const& seq2){
	size_t cpt = 0;

	while(i < seq1->length() && cpt < seq2->length() && (*seq1)[i] == (*seq2)[cpt]) {
		cpt++;
		i++;
	} 
	return cpt;
}


// Entre deux sequences
size_t SAcomp::prefcommSA(size_t i, Sequence_FastA* const& seq){
	size_t cpt = 0;

	while(i < seqA->length() && cpt < seq->length() && (*seqA)[i] == (*seq)[cpt]) {
		cpt++;
		i++;
	} 
	return cpt;
}

int SAcomp::dicho (Sequence_FastA* const& W) {
	
	int d = -1;
	int f = SA.size();
	int m, l;
	Sequence_FastA* dm;
	size_t longueurW = W->length();

	while (d+1 < f) {
		m = (d+f)/2;
		l = prefcommSA(SA[m], W);
		dm = getIseq(m);

		if (l == longueurW) {
			if (l == dm->length()) {
				delete dm;
				return m;
			}
			else {
				f = m;
			}
		}
		else {
			if (l < dm->length() && (*W)[l] < (*dm)[l]) {
				f = m;
			}
			else {
				d = m;
			}
		}	
	}
	//retourne où le mot devrait se positionner
	if (m == SA.size()-1){
		delete dm;
		return -1;
	}
	else {
		if (!m)
			return m;
		else {
			if (prefcommSA(SA[d], W) >= prefcommSA(SA[f], W)){
				delete dm;
				return d;
			}
			else {
				delete dm;
				return f;
			}
		}
	}

}


// Rechercher toutes les occurences d'un facteur
// amélioration possible avec getOccurences.
vector<int> SAcomp::rechercheAll (Sequence_FastA* const& W) {
	int pos_dicho = dicho(W); 
	vector<int> positions; 

	if (pos_dicho == -1) 
		return positions;

	size_t lg = prefcommSA(SA[pos_dicho], W);

	if (lg == W->length()) {
		positions.push_back(SA[pos_dicho]);

		while (pos_dicho < LCP.size() && LCP[pos_dicho] >= lg) {
			pos_dicho++;
			positions.push_back(SA[pos_dicho]);
		}
	}

	return positions;
}

//Rechercher un mot et son complementaire inverse
map< int, vector<int>> SAcomp::rechercheBrins(Sequence_FastA* const& W){
	vector <int> positif = rechercheAll(W);
	vector <int> negatif = rechercheAll((W->complement())->reverse());

	map <int, vector<int>> posneg; 
	posneg[0] = negatif; 
	posneg[1] = positif; 


	return posneg;

}


