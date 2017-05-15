#include "LCPdistance.h"

double LCPdistance::distance(const int&is, const int& js){

    int m=slen[is];
    int n=slen[js];
    // Computing min length
    int minimum = m;
    if (n<m) minimum = n;
    int i;
    if (sign>0) {
        i=0;
        while (sequences[MINDICE(is,i,nseq)]==sequences[MINDICE(js,i,nseq)] && i<minimum) {
            i++;
        }
    } else {
        i=1;
        while (sequences[MINDICE(is,(m-i),nseq)]==sequences[MINDICE(js,(n-i),nseq)] && i<=minimum) {
            i++;
        }
        i--;
    }
    return normalizeDistance((double)n+(double)m-2.0*(double)i, (double)n+(double)m, m, n);
}
