#include "alignement.h"
#include "R.h"
#include "algorithm"
#include "utility"

using namespace std;

Alignement::Alignement(const int& iseq, const int&jseq, const double& pid, const int& liseq, const int& ljseq, int& maxlen) {
	  this->iseq = iseq;
	  this->jseq = jseq;
	  this->pid = pid;
	  this->liseq = liseq;
	  this->ljseq = ljseq;
	  this->maxlen = maxlen;
	  //int maxij = std::max(liseq, ljseq);
	  /*
	  this->alignements = new int*[2];
	  for (int i = 0; i < 2; i++) {
		     alignements[i] = new int[maxlen+maxlen];
		  }
	*/
}
Alignement::~Alignement() {
}
/*
void Alignement::setAlign(stack<pair<int,int> >* alignpairs, double * frequences, double *freqconj, const int&alphasize) {
	int i = 0;
	std::pair<int,int> pr1;
	this->size=alignpairs->size();
	REprintf("Alignement size = %d\n", size);
	while(!alignpairs->empty()) {
		pr1 = alignpairs->top();
		REprintf("pr1.first = %d\n pr1.second = %d\n", pr1.first, pr1.second);
		alignements[0][i] = pr1.first;
		alignements[1][i] = pr1.second;
		alignpairs->pop();
		i++;
	}


}



void Alignement::printAlign() {
	int i;
	for(i=0;i<size;i++) {
		Rprintf("[0][%d] = [%d] | [1][%d] = [%d]\n", i, alignements[0][i], i, alignements[1][i]);
	}
}

void Alignement::setSize(const int& alignsize) {
	this->size = alignsize;
}
*/
