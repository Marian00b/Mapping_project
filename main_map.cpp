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

enum  optionIndex { UNKNOWN, HELP, SEQUENCE, CONCATENE, PREFIXE, SUFFIXE, POLYA, QUALITE, KMER, GAP};

const option::Descriptor usage[] =
{
	{UNKNOWN, 0,"" , ""    ,option::Arg::None, "\nOPTIONS:\n" },
	{HELP,    0,"" , "help",option::Arg::None, "\t--help  \tPrint usage and exit.\n" },
	{SEQUENCE, 0,"n", "num",option::Arg::Optional, "\t--num, -n  \tNuméro des lectures à utiliser avec bord supérieur non inclus.\n" 
																"\t\tExemple usage : -n1,2 \n"},
	{CONCATENE, 0,"t", "concat",option::Arg::None, "\t--concat, -t  \tConcatenation des séquences du génome.\n" },
	{PREFIXE, 0,"p", "pref",option::Arg::Optional, "\t--pref, -p  \tSupprime préfixe d'une longueur donnée pour les Lectures.\n" },
	{SUFFIXE, 0,"f", "suff",option::Arg::Optional, "\t--suff, -f  \tSupprime suffixe d'une longueur donnée pour les lectures.\n" },
	{POLYA, 0,"a", "polya",option::Arg::Optional, "\t--polya, -f  \tSupprime queue polyA d'une longueur minimale donnée pour les lectures. (Defaut: 4)\n" },
	{QUALITE, 0,"q", "qual",option::Arg::Optional, "\t--qual, -u  \tSupprime les parties de mauvaise qualité des lectures selon une qualité\n"
																"\t\tminimale et une longueur de sequence minimale.\n"
																"\t\t(Defaut: Lg = 15 et Qu = 20)\n"
																"\t\tExemple usage : -u10,20 avec 10 la longueur et 20 la qualité\n"},
	{KMER, 0,"k", "kmer",option::Arg::Optional, "\t--kmer, -k  \tTaille des kmers.\n" },
	{GAP, 0,"g", "gap",option::Arg::Optional, "\t--gap, -g  \tTaille du gap pour la recherche partielle.\n" },														


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

	clock_t t;
	t = clock();

	std::ofstream log_file("log_file.txt");

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
	string mot = parse.nonOption(1);
	bool isFastQ = !isFastA(nomFicFasta);
	int beg = 0, end = 0, pref = 0, suff = 0, polya = 4, lg = 15, qu = 20, k = 22, gap = 5;
	SAcomp * sa; 
	vector <FakeSeq*> sequencesF;
	vector <Sequence_FastA*> sequences;

	for (option::Option* opt = options[SEQUENCE]; opt; opt = opt->next()) {

		if (opt->arg){

			stringstream ss; string s; 
			ss << opt->arg; ss >> s;

			size_t pos = s.find(",");  
			beg = stoi(s.substr(0, pos));
			end = stoi(s.substr(pos+1));
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

	for (option::Option* opt = options[KMER]; opt; opt = opt->next()) {

		if (opt->arg){
			stringstream ss; 
			ss << opt->arg; ss >>k;
		}
	}

	for (option::Option* opt = options[GAP]; opt; opt = opt->next()) {

		if (opt->arg){
			stringstream ss; 
			ss << opt->arg; ss >> gap;
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

			size_t s_end; 

			if (options[CONCATENE]) 
				s_end = positions.size();
			else 
				s_end = 1; 				

			vector <FakeSeq*> sequencesF = createAllFakeSeq(nomFicFasta, positions, 0, s_end, false);


			sequences = createAll(nomFicFasta, positions, 0 , s_end, true, sequencesF);

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

			sa = new SAcomp(sequences[0]);

			log_file<<"Génome indexé et LCP crée"<<endl;
			log_file<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;

		}
	}

		/************************************************************************************************
		*************************************************************************************************	
		*************************************************************************************************
		******************************************* OPTION FASTQ ****************************************
		*************************************************************************************************
		*************************************************************************************************	
		*************************************************************************************************/

	if (checkFile(mot)) {

		vector <size_t> positionsQ = extractPosSeqQ(mot);

		if (positionsQ.size() == 0) {
			return 0;
		}

		if (!end) 
			end = positionsQ.size();

		vector <FakeSeq*> sequencesFQ = createAllFakeSeq(mot, positionsQ, beg, end, true);

		vector <Sequence_FastQ*> sequencesQ = createAllQ(mot, positionsQ, beg , end, true, sequencesFQ);

		for (size_t i = 0; i < sequencesFQ.size(); i++){
			delete sequencesFQ[i];
		}

		if (options[POLYA]) {
			for (size_t i = 0; i < sequencesQ.size(); i++){
				sequencesQ[i]->supprimePolyA(polya);
			}
		}

		if (options[PREFIXE]) {
			for (size_t i = 0; i < sequencesQ.size(); i++){
				sequencesQ[i]->supprimePrefixe(pref);
			}
		}

		if (options[SUFFIXE]) {
			for (size_t i = 0; i < sequencesQ.size(); i++){
				sequencesQ[i]->supprimeSuffixe(suff);
			}
		}

		if (options[QUALITE]) {
			for (size_t i = 0; i < sequencesQ.size(); i++){
				sequencesQ[i]->supprimeBadQuality(lg,qu);
			}
		}

		vector<Kmer*> reads = createAllKmers(sequencesQ,k);

		log_file<<"Kmers construits!"<<endl;
		log_file<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;

		map <int, size_t> res; 

		res =sa->searchAllKmers(reads);

		log_file<<"Kmers mappés!"<<endl;
		log_file<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;

		sa->getMappingKmer(res, reads, log_file);


		for (size_t i = 0; i <reads.size(); i++){
			cout<<sequencesQ[i]->getIntitule()<<"\t\t";
			cout<<"\t"<<reads[i]->isSuite(res,gap);
			cout<<"\t";
			printLocalise(reads[i]->getIsLocalise_exact(), reads[i]->getIsLocalise_NonExact()); 
			if (reads[i]->getIsLocalise_exact() || reads[i]->getIsLocalise_NonExact()) 
				cout<<"\t"<<abs(reads[i]->getLocalisation())<<"\t"<<(reads[i]->getLocalisation()>=0 ? "+" : "-")<<"\t";
			else 
				cout<<"\t"<<"."<<"\t"<<"."<<"\t";

			cout<<(*sequencesQ[i])<<"\t";
			sequencesQ[i]->printQualite();
			cout<<endl;
		}

		for (size_t i = 0; i < sequencesQ.size(); i++){
			delete sequencesQ[i];
		}

		for (size_t i = 0; i < reads.size(); i++){
			delete reads[i];
		}

	}

	for (size_t i = 0; i < sequences.size(); i++){
		delete sequences[i];
	}

	log_file<<((float)(clock() -t))/CLOCKS_PER_SEC<<endl;

	log_file.close();

	return 0;
}
