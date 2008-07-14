#include<R.h>
#include <Rinternals.h>

#define TMRMATRIXINDEX(ligne, colone,len) (ligne-1)+(colone-1)*len

/** distmatrix: un vecteur de double représentant la matrice (concaténation des colonnes à la suite les une des autres
L'indexation est assurée par (ligne-1)+(colone-1)*len

    indiv: vecteur numéric des indexs des individus qui forme le groupe (commence à 1)
    matrixsize: taille (nbligne ou colonne de la matrice carrée)
    groupesize: taille du groupe
*/

SEXP tmrsubmatrixinertia(SEXP distmatrix, SEXP individuals) {
    int mlen=nrows(distmatrix);
    int ilen=length(individuals);
    int * indiv=INTEGER(individuals);
    double *dmat=REAL(distmatrix);
  //Rprintf("mlen = %i \n", mlen);
  //Rprintf("ilen = %i \n", ilen);
    int i, j;
    double result=0;
    for (i=0;i<ilen;i++) {
        for (j=i+1;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)],TMRMATRIXINDEX(indiv[i],indiv[j],mlen));
            result+=dmat[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)];
        }
    }

	//Rprintf("Sum = %f\n",result);
    if(ilen>0)result/=(double)ilen;
    //Rprintf("Inertia = %f\n",result);
    return ScalarReal(result);

}
SEXP tmrinertiacontrib(SEXP distmatrix, SEXP individuals) {
    int mlen=nrows(distmatrix);
    int ilen=length(individuals);
    int * indiv=INTEGER(individuals);
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, ilen));
    double *result=REAL(ans);
    double *dmat=REAL(distmatrix);
//  Rprintf("mlen = %i \n", mlen);
//  Rprintf("ilen = %i \n", ilen);
    int i, j;
    for (i=0;i<ilen;i++) {
    	result[i]=0;
        for (j=0;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)],TMRMATRIXINDEX(indiv[i],indiv[j],mlen));
            result[i]+=dmat[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)];
        }
        if(ilen>0)result[i]/=(double)ilen;
    }
//  Rprintf("Sum = %f\n",(*result));
    UNPROTECT(1);
    return ans;
//   Rprintf("Inertia = %f\n",(*result));
}

SEXP tmrinterinertia(SEXP distmatrix, SEXP grp1,SEXP grp2){
    int mlen=nrows(distmatrix);
    int ilen1=length(grp1);
    int ilen2=length(grp2);
    int * indiv1=INTEGER(grp1);
    int * indiv2=INTEGER(grp2);
    double *dmat=REAL(distmatrix);
	int i, j;
//	Rprintf("mlen = %i \n", mlen);
//	Rprintf("ilen1 = %i \n", ilen1);
//	Rprintf("ilen2 = %i \n", ilen2);

    double result=0;
    for (i=0;i<ilen1;i++) {
        for (j=0;j<ilen2;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[TMRMATRIXINDEX(indiv[i],indiv[j],mlen)],TMRMATRIXINDEX(indiv[i],indiv[j],mlen));
            result+=dmat[TMRMATRIXINDEX(indiv1[i],indiv2[j],mlen)];
        }
    }

	//Rprintf("Sum = %f\n",result);
    //if(ilen>0)result/=(double)ilen;
    //Rprintf("Inertia = %f\n",result);
    return ScalarReal(result);
//   Rprintf("Inertia = %f\n",(*result));
}
