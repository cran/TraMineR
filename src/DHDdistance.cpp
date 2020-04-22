#include "DHDdistance.h"

double DHDdistance::distance(const int&is, const int& js){


    int m=slen[is];
    int n=slen[js];
    // Computing min length
    int minimum = m;
    if (n<m) minimum = n;
    double cost=0;
    for (int i=0;i<minimum;i++) {
        cost += scost[ARINDICE(sequences[MINDICE(is,i,nseq)], sequences[MINDICE(js,i,nseq)], i, alphasize)];
    }
    TMRLOG(5, "DHD distance");
	
    //return normalizeDistance(cost, maxdist, ml, nl);
	// GR maxdist is the sum of maximum scost per position, i.e the length weighted by max cost
    return normalizeDistance(cost, maxdist, maxdist, maxdist);
}
