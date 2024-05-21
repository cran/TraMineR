#include "TraMineR.h"

/** diss: Le vecteur de diss
L'indexation est donnée par
n*(i-1) - i*(i-1)/2 + j-i
=(i-1)*(n-i/2)- i +j

    indiv: vecteur numeric des indexs (sorted!) des individus qui forme le groupe (commence a 0)...............
    matrixsize: taille (nbligne ou colonne de la matrice carree)
    groupesize: taille du groupe
*/
SEXP tmrsubmatrixinertiadiss(SEXP diss, SEXP diss_size, SEXP individuals) {
    int n=INTEGER(diss_size)[0];
    int ilen=Rf_length(individuals);
    int * indiv=INTEGER(individuals);
    double *dmat=REAL(diss);
    //Rprintf("mlen = %i \n", mlen);
    //Rprintf("ilen = %i \n", ilen);
    int i, j, i_indiv, base_indice;
    double result=0;
    for (i=0;i<ilen;i++) {
        i_indiv=indiv[i];
        base_indice = (i_indiv-1)*(n-i_indiv/2)-i_indiv-1;
        for (j=i+1;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[MINDICER(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[MINDICER(indiv[i],indiv[j],mlen)],MINDICER(indiv[i],indiv[j],mlen));
            result+=dmat[base_indice+indiv[j]];
        }
    }

    //Rprintf("Sum = %f\n",result);
    if (ilen>0)result/=(double)ilen;
    //Rprintf("Inertia = %f\n",result);
    return Rf_ScalarReal(result);

}

/** distmatrix: un vecteur de double représentant la matrice (concaténation des colonnes à la suite les une des autres
L'indexation est assurée par (ligne-1)+(colone-1)*len
n*(i-1) - i*(i-1)/2 + j-i

    indiv: vecteur numérique des index des individus qui forme le groupe (commence à 1)
    matrixsize: taille (nbligne ou colonne de la matrice carrée)
    groupesize: taille du groupe
*/
SEXP tmrsubmatrixinertiaCindividuals(SEXP distmatrix, SEXP individuals) {
    int mlen=Rf_nrows(distmatrix);
    int ilen=Rf_length(individuals);
    int * indiv=INTEGER(individuals);
    double *dmat=REAL(distmatrix);
    int i, j, base_index;
    double result=0;
    for (i=0;i<ilen;i++) {
        base_index=indiv[i]*mlen;
        for (j=i+1;j<ilen;j++) {
            result+=dmat[indiv[j]+base_index];
        }
    }
    if (ilen>0)result/=(double)ilen;
    return Rf_ScalarReal(result);

}
SEXP tmrsubmatrixinertia(SEXP distmatrix, SEXP individuals) {
    int mlen=Rf_nrows(distmatrix);
    int ilen=Rf_length(individuals);
    int * indiv=INTEGER(individuals);
    double *dmat=REAL(distmatrix);
    //Rprintf("mlen = %i \n", mlen);
    //Rprintf("ilen = %i \n", ilen);
    int i, j;
    int i_index=0;
    double result=0;
    for (i=0;i<ilen;i++) {
        i_index=(-1)+(indiv[i]-1)*mlen;
        for (j=i+1;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[MINDICER(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[MINDICER(indiv[i],indiv[j],mlen)],MINDICER(indiv[i],indiv[j],mlen));
            result+=dmat[i_index+indiv[j]];
        }
    }

    //Rprintf("Sum = %f\n",result);
    if (ilen>0)result/=(double)ilen;
    //Rprintf("Inertia = %f\n",result);
    return Rf_ScalarReal(result);

}
SEXP tmrinertiacontrib(SEXP distmatrix, SEXP individuals) {
    int mlen=Rf_nrows(distmatrix);
    int ilen=Rf_length(individuals);
    int * indiv=INTEGER(individuals);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(REALSXP, ilen));
    double *result=REAL(ans);
    double *dmat=REAL(distmatrix);
//  Rprintf("mlen = %i \n", mlen);
//  Rprintf("ilen = %i \n", ilen);
    int i_index=0;
    int i, j;
    double r;
    for (i=0;i<ilen;i++) {
        result[i]=0;
    }
    for (i=0;i<ilen;i++) {
        i_index=(-1)+(indiv[i]-1)*mlen;
        for (j=i+1;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[MINDICER(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[MINDICER(indiv[i],indiv[j],mlen)],MINDICER(indiv[i],indiv[j],mlen));
            r=dmat[i_index+indiv[j]];
            result[i]+=r;
            result[j]+=r;
        }
        if (ilen>0)result[i]/=(double)ilen;
    }
//  Rprintf("Sum = %f\n",(*result));
    UNPROTECT(1);
    return ans;
//   Rprintf("Inertia = %f\n",(*result));
}

SEXP tmrinertiacontribext(SEXP distmatrix, SEXP individuals, SEXP extindivS) {
    int mlen=Rf_nrows(distmatrix);
    int ilen=Rf_length(individuals);
    int ilenExt=Rf_length(extindivS);
    int * indiv=INTEGER(individuals);
    int * indivExt=INTEGER(extindivS);
    int totlen=ilen+ilenExt;
    SEXP ans;
    PROTECT(ans = Rf_allocVector(REALSXP, totlen));
    double *result=REAL(ans);
    double *dmat=REAL(distmatrix);

//  Rprintf("mlen = %i \n", mlen);
//  Rprintf("ilen = %i \n", ilen);
    int i_index=0;
    int r_index=0;
    int i, j;
    double r;
    for (i=0;i<totlen;i++) {
        result[i]=0;
    }
    for (i=0;i<ilen;i++) {
        i_index=(-1)+(indiv[i]-1)*mlen;
        for (j=i+1;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[MINDICER(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[MINDICER(indiv[i],indiv[j],mlen)],MINDICER(indiv[i],indiv[j],mlen));
            r=dmat[i_index+indiv[j]];
            result[i]+=r;
            result[j]+=r;
        }
        result[i]/=(double)ilen;
    }
    for (i=0;i<ilenExt;i++) {
        i_index=(-1)+(indivExt[i]-1)*mlen;
        r_index=i+ilen;
        for (j=0;j<ilen;j++) {
            result[r_index]+=dmat[i_index+indiv[j]];
        }
        result[r_index]/=(double)ilen;
    }
//  Rprintf("Sum = %f\n",(*result));
    UNPROTECT(1);
    return ans;
//   Rprintf("Inertia = %f\n",(*result));
}


SEXP tmrinertiacontribdiss(SEXP diss, SEXP diss_size, SEXP individuals) {
    int ilen=Rf_length(individuals);
    int * indiv=INTEGER(individuals);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(REALSXP, ilen));
    double *result=REAL(ans);
    double *dmat=REAL(diss);
//  Rprintf("mlen = %i \n", mlen);
//  Rprintf("ilen = %i \n", ilen);
    int i, j, i_indiv, base_indice;
    int n=INTEGER(diss_size)[0];
    double r;
    for (i=0;i<ilen;i++) {
        result[i]=0;
    }
    for (i=0;i<ilen;i++) {
        i_indiv=indiv[i];
        //base_indice = (i_indiv-1)*(n-i_indiv/2)-i_indiv-1;
        //base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 + i_indiv;
        base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 -i_indiv-1;
        for (j=i+1;j<ilen;j++) {
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[MINDICER(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[MINDICER(indiv[i],indiv[j],mlen)],MINDICER(indiv[i],indiv[j],mlen));
            r=dmat[base_indice+indiv[j]];
            result[i]+=r;
            result[j]+=r;
        }
    }
    if (ilen>0) {
        r=(double)ilen;
        for (i=0;i<ilen;i++) {
            result[i]/=r;
        }
    }
//  Rprintf("Sum = %f\n",(*result));
    UNPROTECT(1);
    return ans;
//   Rprintf("Inertia = %f\n",(*result));
}


SEXP tmrinterinertia(SEXP distmatrix, SEXP grp1,SEXP grp2) {
    int mlen=Rf_nrows(distmatrix);
    int ilen1=Rf_length(grp1);
    int ilen2=Rf_length(grp2);
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
//      Rprintf("index coord(%i,%i)=%f\n",indiv[i],indiv[j], distmatrix[MINDICER(indiv[i],indiv[j],mlen)]);
//      Rprintf("cindex coord(%i,%i)=%f, (%i)\n",i,j, distmatrix[MINDICER(indiv[i],indiv[j],mlen)],MINDICER(indiv[i],indiv[j],mlen));
            result+=dmat[MINDICER(indiv1[i],indiv2[j],mlen)];
        }
    }

    //Rprintf("Sum = %f\n",result);
    //if(ilen>0)result/=(double)ilen;
    //Rprintf("Inertia = %f\n",result);
    return Rf_ScalarReal(result);
//   Rprintf("Inertia = %f\n",(*result));
}
