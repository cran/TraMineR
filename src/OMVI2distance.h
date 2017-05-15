#ifndef OMVI2distanceCALCULATOR_H
#define OMVI2distanceCALCULATOR_H
#include "OMdistance.h"


#define OMVI2TOL 0.0000001
class OMVI2distance: public OMdistance{
	double timecost;
	double localcost;
	int * opmat;
	int *prev_j_state;
	int *prev_i_state;

  public:
	OMVI2distance(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS);
	OMVI2distance(OMVI2distance *dc);
	virtual DistanceCalculator* copy(){return new OMVI2distance(this);}
	virtual void setParameters(SEXP params);
    virtual ~OMVI2distance();
    virtual double distance(const int&is, const int& js);
	double getFirstIndel(const int&state, const int& previous, const int& i);
	double getIndel(const int&state, const int& fmat_pos);
};

inline double OMVI2distance::getFirstIndel(const int& state, const int& previous, const int& i) {
  if (i > 1) {
    return timecost + localcost * scost[MINDICE(previous, state, alphasize)];
  } else {
    return timecost + localcost * maxscost;
  }
}

inline double OMVI2distance::getIndel(const int&state, const int& fmat_pos){
	int previ=this->prev_i_state[fmat_pos], prevj=this->prev_j_state[fmat_pos];
	TMRLOG(7,"state =%d, previ=%d, prevj=%d\n", state, previ, prevj);
	int curstate_indice = state*alphasize;
	if(previ==prevj){
		return timecost+localcost*(scost[curstate_indice+previ]);
	}
	else{
		return timecost+localcost*(fmin2(scost[curstate_indice+previ], scost[curstate_indice+prevj]));
	}
}

#endif
