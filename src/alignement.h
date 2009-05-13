#ifndef _ALIGNEMENT_INCLUDED_
#define _ALIGNEMENT_INCLUDED_
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <algorithm>
#include <utility>
#include <stack>

#define TMRMATRIXINDEXC(ligne, colone,len) (ligne)+(colone)*(len)
using namespace std;
class Alignement {


protected:
	int iseq;
	int jseq;
	int liseq;
	int ljseq;
	int maxlen;
	int size;
	double pid;
//	int ** alignements;
public:
	Alignement(const int& iseq, const int&jseq, const double& pid, const int& liseq, const int& ljseq, int& maxlen);
	virtual ~Alignement();
//	void setAlign(stack<pair<int,int> >* alignpairs, double * frequences, double * freqconj, const int& alphasize);
	void printAlign();
//	void setSize(const int& size);
	inline int getliseq() {
	  return liseq;

	}
	inline int getljseq() {
	  return ljseq;
	}
	inline int getiseq() {
	  return iseq;
	}
	inline int getjseq() {
	  return jseq;
	}
	inline double getpid() {
	  return pid;
	}
	inline void setpid(const double& pid) {
	  this->pid=pid;
	}
//	double getPID();
};
#endif
