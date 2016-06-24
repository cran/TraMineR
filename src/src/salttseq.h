#ifndef _SALTTSEQ_INCLUDED_
#define _SALTTSEQ_INCLUDED_
#include <R.h>
#include <Rinternals.h>
#ifdef length
#undef length
#endif
#include <Rmath.h>
#include "alignement.h"
#include <stack>
#include <set>
#define TMRMATRIXINDEXC(ligne, colone,len) (ligne)+(colone)*(len)
class Salttseq {

protected:
	int norm;
	int nseq;
	int maxlen;
	int * sequences;
	int * slen;
	double * indel;
	int alphasize;
	double * scost;
	double * distmatrix;
	double * fmat;
	double * tbmat;
	double * salttcost;
	int fmatsize;
	double maxscost;
	std::stack<Alignement>* stackAlign;
	//std::stack<Alignement>* stackAlign2;
	double* frequences;
    double* freqconj;
    double* dayhoffcosts;
    double pid;
    double totalstates;
    double totalalignements;
    double alphafois;
    int log_odd_mode;



public:
	//Salttseq(const int& norm, const int& nseq, const int&maxlen, int * sequences, int * slen, double & indel, int & alphasize, double * costs, double * distmatrix);

	Salttseq(const int& norm, const int& nseq, int * slen, const int& maxlen, double * indel, int & alphasize, double * costs, double * distmatrix, std::stack<Alignement>* stackAlign, int *sequences, double *salttcost, double * fmat, double * tbmat, const int& fmatsize, const double& pid, const int& logoddmode);
	virtual ~Salttseq();
	int computeDistances(const int& withsaltt);
	double* getdist();
	double levenshtein(const int &is, const int &js, int &m, int &n);
	double normalizeDistance(const double& rawdist, const double& maxdist, const int& l1, const int& l2, const int&norm);
	double getPID(const int& is, const int &js, const int &m, const int &n, Alignement& align);
	void printfreqtables();
	void computeFreq();
	void computeConjfreq(const int &is, const int &js, const int &m, const int &n, const double& currentpid);
	void printsalttmatrix();
	void computeDayhoffcosts();
	void moveStack(std::stack<Alignement>* oldstack, std::stack<Alignement>* newstack);
	void computeConjfreq2(const int &is, const int &js, const int &m, const int &n, const double& currentpid);
	void printtraceback(int m, int n);
/*	inline double getmaxGap() {

		return(this->maxGap);
	}
	*/
};
#endif
