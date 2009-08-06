#include "TraMineR.h"

SEXP dist2matrix(SEXP diss, SEXP diss_size) {
    int n=INTEGER(diss_size)[0];
    SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP, n, n));
    double *matrix=REAL(ans);
    double *dmat=REAL(diss);
//  Rprintf("mlen = %i \n", mlen);
//  Rprintf("ilen = %i \n", ilen);
    int i, j, i_indiv, base_indice;
    double r;
    for (i=0;i<n;i++) {
        i_indiv=i+1;
        base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 -i_indiv-1;
        matrix[i+i*n]=0;
        for (j=i+1;j<n;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)],TMRMATRIXINDEX(indiv[i],indiv[j],mlen));
            r=dmat[base_indice+j+1];
            matrix[j+i*n]=r;
            matrix[i+j*n]=r;
        }
    }
//  Rprintf("Sum = %f\n",(*result));
    UNPROTECT(1);
    return ans;
//   Rprintf("Inertia = %f\n",(*result));
}
