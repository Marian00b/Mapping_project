#include "SA.h" 


vector<int> SAcomp::rechercheMMP(Sequence_FastA* const& kmer){

	map< int, vector<int>> resultat; 
	vector<int> retour;
	size_t i = kmer->length()-2;

	// cherche au moins 1 position pour un suffixe
	while (resultat[0].size() == 0 && resultat[1].size() == 0 && i > 0){
		resultat = rechercheBrins(kmer->subSequence(0,i));
		i--;
	}

	//ajout des positions neg et pos pour le retour
	for (map<int,vector<int>>::iterator it=resultat.begin(); it!=resultat.end(); ++it) {
		for (size_t k = 0; k < it->second.size(); k++){
			if (!it->first) 
				retour.push_back(- it->second[k]); 
			else 
				retour.push_back(it->second[k]);
		}
	}

	return retour;
}


int* SAcomp::rechercheK (int pos, Sequence_FastA* const& kmer){
	int * res = new int[2];
	//va contenir gap ou 0 si pas trouve dans les 5 prochaines positions
	res[1] = -1;
	//size_t i = 0;
	size_t cpt = 0;
	res[0] = pos+1; 

	if (pos < 0) {
		pos = abs(pos);
		--pos;
		Sequence_FastA* k = (kmer->complement())->reverse();


		// on trouve premiere lettre du kmer si négatif
		while (((*k)[0] != (*seqA)[(pos)]) && cpt < 5){
			//pos++;
			pos--;
			cpt++;
		}

		if (cpt != 5)
			res[1] = cpt; 

	}
	else {
		pos += kmer->length();

		//on cherche la derniere lettre 
		while ((*kmer)[kmer->length()-1] != (*seqA)[pos] && cpt < 5){
			pos++;
			cpt++;
		}

		if (cpt != 5)
			res[1] = cpt; 
	}

	return res; 
}



map <int,size_t> SAcomp::searchAllKmers(vector<Kmer*> & reads){

	map < int, vector<int> > pos ;
	map <int,size_t> cov; 
	int value; 

	for (size_t i = 0; i < reads.size(); i++) {

		for (size_t j = 0; j < reads[i]->nombre(); j++) {

			pos = rechercheBrins(reads[i]->getKmer(j));

			for (map<int,vector<int>>::iterator it=pos.begin(); it!=pos.end(); ++it) {
				for (size_t k = 0; k < it->second.size(); k++){
					// pour négatif on veut afficher quoi : debut ou fin?
					if (!it->first) 
						value = - (it->second[k] + reads[i]->nombre() - 1) ; 
					else 
						value = it->second[k];

					if (cov.find(value) == cov.end())
						cov[value] = 1; 
					else 
						cov[value] = cov[value] + 1;

					reads[i]->setMatch(j, value);
				}
			}

		}
		//cout<<" ]"<<endl;

	}

	return cov; 
}

void SAcomp::getMappingKmer (map <int, size_t> & res, vector<Kmer*> const& reads, ofstream & log_file) {

	vector<int> pos, pos2, precedent;
	int value; 
	int * recherche; 

	// afficher pour chaque read, les positions de chacun de ses k-mers, et leur couverture	

	for (size_t i = 0; i < reads.size(); i++) {
		log_file<<"Read "<<i<<" :"<<endl;

		for (size_t j = 0; j < reads[i]->nombre(); j++) {
			log_file<<"\tKmer "<<j<<"\t";
			pos = ((reads[i])->getPositions(j)); 

			// Si le kmer n'a pas mappé 
			if (!pos.size()){
				// si c'est celui en premiere positions
				if (!j){
					(reads[i])->setNewNonMatch(j,(rechercheMMP(reads[i]->getKmer(j))));
				}
				else {

					// position de son kmer mappé d'avant
					// fonction si jamais pas avant 
					precedent = ((reads[i])->getPositions(j-1));

					// if (!precedent.size())
					// 	precedent = ((reads[i])->getPositions(j-1));

					// pour toutes les positions où le precedent a mappe 
					for (size_t k = 0; k< precedent.size(); k++) {
						recherche = rechercheK(precedent[k], reads[i]->getKmer(j));
						// affiche avec le gap pour coincider avec suivant
						if (recherche[1]>=0){
							value = recherche[0]+recherche[1]; 
							(reads[i])->setNonMatch(j,value);

							// calcul de la couverture
							if (res.find(value) == res.end()){
								res[value] = 1; 
							}
							else 
								res[value] = res[value] + 1;
						}
					}
				}

				// affichage des positions du kmer si non exact
				pos2 = ((reads[i])->getPositions(j));

				for (size_t k = 0; k < pos2.size(); k++){
					log_file<<"N"<<res[pos2[k]]<<"\t";
					log_file<<pos2[k]<<"\t";
				} 
				
				
			}

			// affichage si exact
			for (size_t k = 0; k < pos.size(); k++){
				if (!k)
					log_file<<res[pos[k]]<<"\t";
				log_file<<pos[k]<<"\t";

			} 
			log_file<<endl;

		}
	}

}


void printLocalise(unsigned char i,unsigned char j ){
	if (i|j){
		if (i)
			cout<<"Exact";
		else 
			cout<<"Non_Exact";
	}
	else 
		cout<<"Non_Localise";
}
