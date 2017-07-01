#include "suffixe.h"
#include <algorithm> 
#include <numeric>

using namespace std; 

Suffixe::Suffixe():s(),suffixes() {}

Suffixe::Suffixe(string seq) : s(seq) {

	for (size_t i = 0; i < s.length(); i++) {
		suffixes.push_back(seq.substr(i));
	}
}

Suffixe::Suffixe(Sequence_FastA const& seq): Suffixe()  {
	//tableau ordonné des suffixes
	for (size_t i = 0; i < seq.length(); i++) {
		suffixesA.push_back(seq.subSequence(i, seq.length()-1));
	}
	

}

bool cmp(const Sequence_FastA &a, const Sequence_FastA &b) {
	size_t i = 0; 
	bool aInfB = false; 

	while (i < a.length() && i < b.length() && a[i] == b[i]) {
		i++;
	}

	if (i != b.length()) {
		if (i == a.length()) 
				aInfB = true; 
		else {
			if (a[i] < b[i]) 
				aInfB = true;
		}
	}

	return aInfB;
}


bool cmp_str(const string &a, const string &b) {
	size_t i = 0; 
	bool aInfB = false; 

	while (i < a.length() && i < b.length() && toupper(a[i]) == toupper(b[i])) {
		i++;
	}

	if (i != b.length()) {
		if (i == a.length()) {
			aInfB = true; 
		} 
		else {
			if (toupper(a[i]) < toupper(b[i])) {
				aInfB = true;
			}
		}
	}
	return aInfB;
}


void Suffixe::print() {
	if (!s.empty()) {
		for (size_t i = 0; i < suffixes.size(); i++) {
			cout<<suffixes[i]<<endl;
		}
	}
	else {
		//char test = 0; 
		for (size_t i = 0; i < suffixesA.size(); i++) {
		 	cout<<i<<" "<<(*suffixesA[i])<<endl;
			//if (i) 
			// 	test = test|!(cmp((*suffixesA[i-1]), (*suffixesA[i])));
		}
	}
	//cout<<int(test)<<endl;


}

void Suffixe::trier() {
	if (!s.empty()) 
		sort(suffixes.begin(), suffixes.end(), cmp_str);
	else {
		//sort(suffixesA.begin(), suffixesA.end(), cmp);
		make_heap(suffixesA.begin(), suffixesA.end(), cmp);
		sort_heap(suffixesA.begin(), suffixesA.end(), cmp);
	}
}


