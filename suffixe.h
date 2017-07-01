#ifndef SU_H
#define SU_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include "Sequence_FastQ.h"

class Suffixe {

	private : 
		std::string s; 
		Sequence_FastA seqA;
		std::vector <std::string> suffixes; 
		std::vector <Sequence_FastA*> suffixesA; 
		std::vector <size_t> SA;

	public :
		Suffixe();
		Suffixe(std::string seq); 
		Suffixe(Sequence_FastA const& seq); 
		void print();
		void trier();
		//~Suffixe(){

};

bool cmp_str(const std::string &a, const std::string &b);
bool cmp(const Sequence_FastA &a, const  Sequence_FastA &b);

#endif