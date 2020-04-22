#include "TraMineR.h"



// Getting a element by its names


extern "C" {

    SEXP tmrChisq(SEXP ChiTableS, SEXP tdimS, SEXP margeS) {
		TMRLOG(5, "tmrChisq\n");
        SEXP distS;
		int n = INTEGER(tdimS)[0];
		int n1 = n-1;
		PROTECT(distS=allocVector(REALSXP, n*(n-1)/2));
		double * dist=REAL(distS);
		int col =INTEGER(tdimS)[1];
		double * chitable=REAL(ChiTableS);
		double * marge=REAL(margeS);

		for(int i=0; i< n1; i++){
			int base_index= TMRDISTINDEX(i+1, 1, n);
			for(int j = i+1; j < n; j++){
				
				double cmpres=0;
				for(int c=0; c < col;c++){
					double diff=chitable[MINDICE(i, c, n)]-chitable[MINDICE(j, c, n)];
					cmpres += diff*diff/marge[c];
				}
				dist[base_index+j] = sqrt(cmpres);
			}
		}
        UNPROTECT(1);
        return distS;
    }    

    SEXP tmrChisqRef(SEXP ChiTableS, SEXP tdimS, SEXP margeS, SEXP refid) {
		TMRLOG(5, "tmrChisq\n");
        SEXP distS;
		int n = INTEGER(tdimS)[0];
		int rseq= INTEGER(refid)[0]-1;
		//int n1 = n-1;
		PROTECT(distS=allocVector(REALSXP, n));
		double * dist=REAL(distS);
		int col =INTEGER(tdimS)[1];
		double * chitable=REAL(ChiTableS);
		double * marge=REAL(margeS);

		for(int i=0; i< n; i++){
			//int base_index= TMRDISTINDEX(i+1, 1, n);
			//for(int j = i+1; j < n; j++){
				
			double cmpres=0;
			if(i!=rseq) {
				for(int c=0; c < col;c++){
					double diff=chitable[MINDICE(i, c, n)]-chitable[MINDICE(rseq, c, n)];
					cmpres += diff*diff/marge[c];
				}
			}
			dist[i] = sqrt(cmpres);
		}
        UNPROTECT(1);
        return distS;
    }    
}
