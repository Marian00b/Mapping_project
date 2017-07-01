#include "Sequence_FastQ.h"
#include "FakeSeq.h"
#include "SA.h"
#include "suffixe.h"
#include "optionparser.h"
#include "kmer.h"
#include <sstream>
#include <time.h>   
#include <cmath>

	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	********************************************* LEAN MEAN *****************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

enum  optionIndex { UNKNOWN, HELP, FASTQ , SEQUENCE, CONCATENE, PRINT, PRINTS, SEARCH, ENCODE, OCCURENCE, DISTINCT, FACTEUR, MOT};

const option::Descriptor usage[] =
{
	{UNKNOWN, 0,"" , ""    ,option::Arg::None, "\nOPTIONS:\n" },
	{HELP,    0,"" , "help",option::Arg::None, "\t--help  \tPrint usage and exit.\n" },
	{SEQUENCE, 0,"n", "num",option::Arg::Optional, "\t--num, -n  \tNuméro des séquences à utiliser avec bord supérieur non inclus.\n"
													"\t\tDéfaut : seulement la première. \n" },
	{CONCATENE, 0,"c", "concat",option::Arg::None, "\t--concat, -t  \tConcatenation des séquences.\n" },
	{PRINT, 0,"p", "print",option::Arg::None, "\t--print, -p  \tImprime le tableau ordonné des suffixes.\n" },
	{PRINTS, 0,"i", "prints",option::Arg::None, "\t--prints, -i  \tImprime la table des suffixes.\n" },
	{SEARCH, 0,"s", "search",option::Arg::None, "\t--search, -s  \tCherche une ou plusieurs séquences.\n" },
	{ENCODE,    0,"e" , "enc",option::Arg::None, "\t--enc, -e \tSi la séquence doit etre encodee.\n" },
	{OCCURENCE,    0,"o" , "occ",option::Arg::Optional, "\t--occ, -o \tRecuperer les facteurs de longueur k qui a le plus d'occurrences dans la sequence.\n"
 																"\t\tRecuperation du 1er facteur seulement si plusieurs.\n" },
 	{DISTINCT,    0,"d" , "dist",option::Arg::Optional, "\t--dist, -d \tRecuperer le nombre de facteurs distincts de longueur k\n" },															
 	{FACTEUR,    0,"f" , "fact",option::Arg::Optional, "\t--fact, -f \tRecuperer le ieme facteur de longueur k de la sequence\n"
 																	"\t\tExemple : -fi,k\n"
 																	"\t\t(Defaut : 1,1)\n"},	
 	{MOT,    0,"m" , "mot",option::Arg::Optional, "\t--mot, -m \tRecuperer le mot de longueur K  à la position i\n" 
 																"\t\tExemple : -mi,k\n"
 																"\t\t(Defaut : 1,1)"},														
 	
	{0,0,0,0,0,0}
};

int main(int argc, char const *argv[])
{

	clock_t t;
	t = clock();


	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	********************************************* ARGUMENTS *****************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	option::Stats  stats(usage, argc, argv);
	vector<option::Option> options(stats.options_max);// non POD var size arrays are not C++!
	vector<option::Option> buffer(stats.buffer_max);
	option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

	if (parse.error())
	return 1;
	
	if (options[HELP] || argc == 0) {
		option::printUsage(std::cout, usage);
		return 0;
	}

	string texte = parse.nonOption(0);
	string mot;
	SAcomp * sa; 
	vector <FakeSeq*> sequencesF;
	vector <Sequence_FastA*> sequences;
	vector <Sequence_FastA*> sequencesQ;

	int beg = 0, end = 0, occ = 1, dis = 1, fk =1, fi = 1, mk = 1, mi = 1; 

	for (option::Option* opt = options[SEQUENCE]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s; 
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			beg = stoi(s.substr(0, pos));
			end = stoi(s.substr(pos+1));
		}
	}

	for (option::Option* opt = options[FACTEUR]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s;
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			fi = stoi(s.substr(0, pos));
			fk = stoi(s.substr(pos+1));
		}
	}

	for (option::Option* opt = options[MOT]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s;
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			mi = stoi(s.substr(0, pos));
			mk = stoi(s.substr(pos+1));
		}
	}

	if (options[SEARCH]) {
		mot = parse.nonOption(1);
	}

	
	for (option::Option* opt = options[OCCURENCE]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; 
			ss << opt->arg; ss >> occ;
		}
	}

	for (option::Option* opt = options[DISTINCT]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; 
			ss << opt->arg; ss >> dis;
		}
	}




	if (checkFile(texte)) {

	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	******************************************* TEXTE ***********************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

		vector <size_t> positions = extractPosSeq(texte);

		if (positions.size() == 0) {
			return (0);
		}

		if (!end) 
			end = 1;

		sequencesF = createAllFakeSeq(texte, positions, beg, end, false);


		sequences = createAll(texte, positions, beg , end, true, sequencesF);

		for (size_t i = 0; i < sequencesF.size(); i++){
			delete sequencesF[i];
		}


		if (options[CONCATENE] && sequences.size() > 1){
			Sequence_FastA* c = concatAllSeq(sequences);

			for (size_t i = 0; i < sequences.size(); i++){
				delete sequences[i];
			}

			sequences.clear();
			sequences.push_back(c);
		}

		if (options[PRINT]) {

			Suffixe suff = Suffixe(sequences[0]);
			suff.trier();
			suff.print();
		}

		if (options[SEARCH]||options[OCCURENCE]||options[DISTINCT]||options[FACTEUR]||options[MOT]||options[PRINTS]) {
			sa = new SAcomp(sequences[0]);
			cout<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;
		}

		if (options[PRINTS])
			sa->print();

		if(options[OCCURENCE]){
			cout<<*(sa->getOccurences(occ))<<endl;
		}	

		if(options[DISTINCT]){
			cout<<sa->getDistincts(dis)<<endl;
		}	

		if(options[FACTEUR]){
			cout<<(*sa->getIfacteur(fi,fk))<<endl;
		}

		if(options[MOT]){
			cout<<(*sa->getMot(mk,mi))<<endl;
		}


	}

	if (options[SEARCH] && checkFile(mot)) {
	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	******************************************* MOT *************************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

		if (options[ENCODE]){
			vector <size_t> positions = extractPosSeq(mot);

			if (positions.size() == 0) {
				return (0);
			}

			vector <FakeSeq*> sequencesF = createAllFakeSeq(mot, positions, 0, positions.size(), false);

			sequencesQ = createAll(mot, positions, 0, positions.size(), true, sequencesF);

			for (size_t i = 0; i < sequencesF.size(); i++){
				delete sequencesF[i];
			}

			for (size_t i = 0; i < sequencesQ.size(); i ++ ) {

				cout<<(*sequencesQ[i])<<" trouvé aux positions suivantes : "<<endl;

				map < int, vector<int> > pos = sa->rechercheBrins(sequencesQ[i]);

			 	for (map<int,vector<int>>::iterator it=pos.begin(); it!=pos.end(); ++it) {
			 		if (it->second.size()) {
			 			if (!it->first)
			 				cout<<"Brin négatif => ";
			 			else 
			 				cout<<"Brin positif => ";
			 		}
					for (size_t i = 0; i < it->second.size(); i++){
							if (!it->first)
								cout<<(-it->second[i])<<" "; 
							else 
								cout<<it->second[i]<<" ";
					}
					if (it->second.size())
						cout<<endl;
				}

			}

			for (size_t i = 0; i < sequencesQ.size(); i++){
				delete sequencesQ[i];
			}


		}

		else {

			ifstream Ifmot (mot.c_str(), ios::in);

			string mot = "";
			char carac ;

			vector <int> pos;

			while (Ifmot.get(carac)){
				if (carac != '\n')
					mot+= carac; 
			}

			cout<<mot<<" trouvé aux positions: [ ";

			pos = sa->rechercheAllM(mot);

			for (size_t i = 0; i < pos.size(); i++){
				cout<<pos[i]<<" ";
			}
			cout<<"]"<<endl;
			mot = "";

			Ifmot.close();

		}


	}

	for (size_t i = 0; i < sequences.size(); i++){
		delete sequences[i];
	}

	cout<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;

	return 0;
}
