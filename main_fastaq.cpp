#include "Sequence_FastQ.h"
#include "FakeSeq.h"
#include "SA.h"
#include "suffixe.h"
#include "optionparser.h"
#include <sstream>
#include <time.h>   

	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	********************************************* LEAN MEAN *****************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

enum  optionIndex { UNKNOWN, HELP, FASTQ , SEQUENCE, SUBSEQUENCE, REVERSE, COMPLEMENT, COUNT, ENTETE, CONCATENE, INFO,
					PREFIXE, SUFFIXE, POLYA, QUALITE};

const option::Descriptor usage[] =
{
	{UNKNOWN, 0,"" , ""    ,option::Arg::None, "\nOPTIONS:\n" },
	{HELP,    0,"" , "help",option::Arg::None, "\t--help  \tPrint usage and exit.\n" },
	{FASTQ, 0,"q", "fastq",option::Arg::None, "\t--fastq, -q  \tSi il s'agit d'un FastQ.\n" },
	{SEQUENCE, 0,"n", "num",option::Arg::Optional, "\t--num, -n  \tNuméro des séquences à utiliser avec bord supérieur non inclus.\n" 
																"\t\tExemple usage : -n1,2 \n"},
	{SUBSEQUENCE, 0,"s", "seq",option::Arg::Optional, "\t--seq, -s  \tSous-sequence avec extrémités incluses.\n"
																	"\t\tExemple usage : -s1,2 \n"},
	{REVERSE, 0,"r", "rev",option::Arg::None, "\t--rev, -r  \tInverser la/les séquences.\n" },
	{COMPLEMENT, 0,"c", "cmp",option::Arg::None, "\t--cmp, -c  \tSequence complementaire.\n" },
	{COUNT, 0,"w", "count",option::Arg::None, "\t--count, -w  \tSeulement compter les séquences.\n" },
	{ENTETE, 0,"e", "entete",option::Arg::None, "\t--entete, -e  \tAnalyse les entetes.\n" },
	{CONCATENE, 0,"t", "concat",option::Arg::None, "\t--concat, -t  \tConcatenation des séquences.\n" },
	{INFO, 0,"i", "info",option::Arg::None, "\t--info, -i  \tRecupère toutes les informations sur les séquences.\n" },
	{PREFIXE, 0,"p", "pref",option::Arg::Optional, "\t--pref, -p  \tSupprime préfixe d'une longueur donnée.\n" },
	{SUFFIXE, 0,"f", "suff",option::Arg::Optional, "\t--suff, -f  \tSupprime suffixe d'une longueur donnée.\n" },
	{POLYA, 0,"a", "polya",option::Arg::Optional, "\t--polya, -f  \tSupprime queue polyA d'une longueur minimale donnée. (Defaut: 4)\n" },
	{QUALITE, 0,"u", "qual",option::Arg::Optional, "\t--qual, -u  \tSupprime les parties de mauvaise qualité selon une qualité\n"
																"\t\tminimale et une longueur de sequence minimale.\n"
																"\t\t(Defaut: Lg = 15 et Qu = 20)\n"
																"\t\tExemple usage : -u10,20 avec 10 la longueur et 20 la qualité\n"},


	{0,0,0,0,0,0}
};

int main(int argc, char const *argv[])
{
	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	********************************************* ARGUMENTS *****************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

	//clock_t t;
	//t = clock();

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

	string nomFicFasta = parse.nonOption(0);
	bool isFastQ = !isFastA(nomFicFasta);
	int beg = 0, end = 0, s_beg =0, s_end = 0, pref = 0, suff = 0, polya = 4, lg = 15, qu = 20;

	if (options[FASTQ])
		isFastQ = true;

	for (option::Option* opt = options[SEQUENCE]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s; 
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			beg = stoi(s.substr(0, pos));
			end = stoi(s.substr(pos+1));
		}
	}


	for (option::Option* opt = options[SUBSEQUENCE]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s; 
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			s_beg = stoi(s.substr(0, pos));
			s_end = stoi(s.substr(pos+1));
		}
	}


	for (option::Option* opt = options[QUALITE]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s; 
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			lg = stoi(s.substr(0, pos));
			qu = stoi(s.substr(pos+1));
		}
	}


	for (option::Option* opt = options[PREFIXE]; opt; opt = opt->next()) {

		if (opt->arg){
			stringstream ss; 
			ss << opt->arg; ss >> pref;
		}
	}

	for (option::Option* opt = options[SUFFIXE]; opt; opt = opt->next()) {

		if (opt->arg){
			stringstream ss; 
			ss << opt->arg; ss >> suff;
		}
	}

	for (option::Option* opt = options[POLYA]; opt; opt = opt->next()) {

		if (opt->arg){
			stringstream ss; 
			ss << opt->arg; ss >> polya;
		}
	}


	if (checkFile(nomFicFasta)) {


	/************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************
	******************************************* OPTION FASTA ****************************************
	*************************************************************************************************
	*************************************************************************************************	
	*************************************************************************************************/

		if(!isFastQ) {
			vector <size_t> positions = extractPosSeq(nomFicFasta);

			if (positions.size() == 0) {
				return (0);
			}

			cout<<endl<<"Il y a "<<positions.size()<<" séquences présentes dans le fichier!"<<endl;

			if (options[ENTETE]) {

				for (size_t i = 0; i < positions.size(); i++){
					cout<<endl<<"Séquence "<<i<<" =>"<<endl;
					printInfoEntente(analyseEntete(nomFicFasta,positions[i], '|'));
				}
			}

			if (!end) 
				end = positions.size();

			vector <FakeSeq*> sequencesF = createAllFakeSeq(nomFicFasta, positions, beg, end, false);

			if (options[INFO]) {
				for (size_t i = 0; i < sequencesF.size(); i++){
					sequencesF[i]->print();
				}
			}

			if (!options[COUNT]) {

				vector <Sequence_FastA*> sequences = createAll(nomFicFasta, positions, beg , end, true, sequencesF);

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


				if (options[SUBSEQUENCE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						Sequence_FastA* s = sequences[i]->subSequence(s_beg,s_end);
						delete sequences[i];
						sequences[i] = s;
					}
				}

				if (options[REVERSE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						Sequence_FastA* s = sequences[i]->reverse();
						delete sequences[i];
						sequences[i] = s;
					}
				}
			
				if (options[COMPLEMENT]) {
					for (size_t i = 0; i < sequences.size(); i++){
						Sequence_FastA* s = sequences[i]->complement();
						delete sequences[i];
						sequences[i] = s; 
					}
				}

				for (size_t i = 0; i < sequences.size(); i++){
					cout<<endl<<(*sequences[i])<<endl;
				}

				for (size_t i = 0; i < sequences.size(); i++){
					delete sequences[i];
				}

			}


	}

		/************************************************************************************************
		*************************************************************************************************	
		*************************************************************************************************
		******************************************* OPTION FASTQ ****************************************
		*************************************************************************************************
		*************************************************************************************************	
		*************************************************************************************************/

		else {

			vector <size_t> positions = extractPosSeqQ(nomFicFasta);

			if (positions.size() == 0) {
				return 0;
			}



			cout<<endl<<"Il y a "<<positions.size()<<" séquences présentes dans le fichier!"<<endl;

			if (options[ENTETE]) {

				for (size_t i = 0; i < positions.size(); i++){
					cout<<endl<<"Séquence "<<i<<" =>"<<endl;
					printInfoEntente(analyseEntete(nomFicFasta,positions[i], ':'));
				}
			}

			if (!end) 
				end = positions.size();

			vector <FakeSeq*> sequencesF = createAllFakeSeq(nomFicFasta, positions, beg, end, true);

			if (options[INFO]) {
				for (size_t i = 0; i < sequencesF.size(); i++){
					sequencesF[i]->print();
				}
			}

			if (!options[COUNT]) {

				vector <Sequence_FastQ*> sequences = createAllQ(nomFicFasta, positions, beg , end, true, sequencesF);

				for (size_t i = 0; i < sequencesF.size(); i++){
					delete sequencesF[i];
				}


				// if (options[CONCATENE] && sequences.size() > 1){
				// 	Sequence_FastQ* c = concatAllSeqQ(sequences);

				// 	for (size_t i = 0; i < sequences.size(); i++){
				// 		delete sequences[i];
				// 	}

				// 	sequences.clear();
				// 	sequences.push_back(c);
				// }

				if (options[POLYA]) {
					for (size_t i = 0; i < sequences.size(); i++){
						sequences[i]->supprimePolyA(polya);
					}
				}

				if (options[PREFIXE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						sequences[i]->supprimePrefixe(pref);
					}
				}

				if (options[SUFFIXE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						sequences[i]->supprimeSuffixe(suff);
					}
				}

				if (options[SUBSEQUENCE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						Sequence_FastQ* s = sequences[i]->subSequence(s_beg,s_end);;
						delete sequences[i];
						sequences[i] = s;
					}
				}

				if (options[QUALITE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						sequences[i]->supprimeBadQuality(lg,qu);
					}
				}

				if (options[CONCATENE] && sequences.size() > 1){
					Sequence_FastQ* c = concatAllSeqQ(sequences);

					for (size_t i = 0; i < sequences.size(); i++){
						delete sequences[i];
					}

					sequences.clear();
					sequences.push_back(c);
				}


				if (options[REVERSE]) {
					for (size_t i = 0; i < sequences.size(); i++){
						Sequence_FastQ* s= sequences[i]->reverse();
						delete sequences[i];
						sequences[i] = s; 
					}
				}
			
				if (options[COMPLEMENT]) {
					for (size_t i = 0; i < sequences.size(); i++){
						Sequence_FastQ* s = sequences[i]->complement();
						delete sequences[i];
						sequences[i] = s; 
					}
				}

				for (size_t i = 0; i < sequences.size(); i++){
					cout<<endl<<(*sequences[i])<<endl;
					sequences[i]->printQualite();
					cout<<endl;
				}

				for (size_t i = 0; i < sequences.size(); i++){
					delete sequences[i];
				}
			}


		}
	}

		/************************************************************************************************
		*************************************************************************************************	
		*************************************************************************************************
		******************************************* DELETE***********************************************
		*************************************************************************************************
		*************************************************************************************************	
		*************************************************************************************************/


	//cout<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;



	return 0;
}
