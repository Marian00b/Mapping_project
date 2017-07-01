#ifndef SEQX_H
#define SEQX_H

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <limits>

using namespace std;

class Sequence_FastX 
{
	protected:
		string intitule;
		size_t longueur;

	public:
		Sequence_FastX();
		Sequence_FastX(string titre, size_t lg = 0) ;
		Sequence_FastX(const Sequence_FastX* const &s);
		bool operator==(Sequence_FastX const &s);
		Sequence_FastX* operator=(Sequence_FastX const &s);
	 	virtual void print() const = 0;	
	    string getIntitule() const ;
	    void setIntitule(string s);
 	    size_t length() const; 
		virtual ~Sequence_FastX();

};

void printInfoEntente (vector <string> analyse);
vector <string> analyseEntete (string nomFicFasta, size_t position, char separator);
bool isFastA (string nomFicFasta);
bool checkFile (string nomFicFasta);
bool compareEntete (vector <string> a1, vector <string> a2);

#endif 
