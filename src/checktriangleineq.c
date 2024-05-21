#include "TraMineR.h"

SEXP checktriangleineq(SEXP mat, SEXP matsize, SEXP tolS) {
    int n=INTEGER(matsize)[0];
	double tol = REAL(tolS)[0];
    SEXP ans;
    
    double *matrix=REAL(mat);
//  Rprintf("mlen = %i \n", mlen);
//  Rprintf("ilen = %i \n", ilen);
    int i, j, z, i_indiv, j_indiv;
    double d;
    for (i=0;i<n;i++) {
        i_indiv=i*n;
        matrix[i+i*n]=0;
        for (j=i+1;j<n;j++) {
			d = matrix[j+i_indiv];
			j_indiv=j*n;
			for (z=0; z<n; z++) {
				if (d-(matrix[i_indiv+ z] + matrix[j_indiv+z]) >= tol) {
					PROTECT(ans = Rf_allocVector(INTSXP, 3));
					INTEGER(ans)[0] =i+1;
					INTEGER(ans)[1] =j+1;
					INTEGER(ans)[2] =z+1;
					UNPROTECT(1);
					return ans;
				}

			}
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)],TMRMATRIXINDEX(indiv[i],indiv[j],mlen));
        }
    }
//  Rprintf("Sum = %f\n",(*result));
   return R_NilValue;
//   Rprintf("Inertia = %f\n",(*result));
}
