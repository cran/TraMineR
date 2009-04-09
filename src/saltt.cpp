#include<R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
//#include <math.h>

/**

*/

#define TMRMATRIXINDEXC(ligne, colone,len) (ligne)+(colone)*(len)
//#define TMRMIN(a,b) ((a)<(b))?a:b
static R_INLINE double normalizeDistance(const double& rawdist, const double& maxdist, const int& l1, const int& l2, const int&norm) {
    if (rawdist==0)return 0;
    switch (norm) {
    case 0:
        return rawdist;
    case 1:
        if (l1>l2)return rawdist/((double)l1);
        else if (l2>0) return rawdist/((double)l2);
        return 0;
    case 2:
        if (l1*l2==0) {
            if (l1!=l2)return 1;
            return 0;
        }
        return 1-((maxdist-rawdist)/(2*R_pow(((double)l1*l2),0.5)));
    case 3:
        if (maxdist==0)return 1;
        return rawdist/maxdist;
    }
}


extern "C" {

    SEXP saltt(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP indelS, SEXP alphasizeS, SEXP costsS, SEXP normS) {
        //Objet R
        SEXP ans, Fmat;
        //Indices, avec s pour s�quences
        int i, j, is, js;
        //longueur des s�quences m=i, n=j
        int m, n;
        //Couts de subsistutions
        double cost;
        //normalisation?
        int norm=INTEGER(normS)[0];
        //Nombre de s�quence
        int nseq=INTEGER(seqdim)[0];
        //nb colonnes des s�quences
        int maxlen=INTEGER(seqdim)[1];
        //Matrice des s�quences
        int* sequences=INTEGER(Ssequences);
        //Tailles des s�quences
        int* slen=INTEGER(lenS);
        //indel
        double indel=REAL(indelS)[0];
        //nb �tats
        int alphasize=INTEGER(alphasizeS)[0];
        //Matrice des co�ts de substitutions
        double* scost=REAL(costsS);

        double maxscost=0;

        //Alocation du vecteur de distance
        //REprintf("Final seq %d\n",finalnseq);
        PROTECT(ans = allocVector(REALSXP, nseq*nseq));
        double * distmatrix = REAL(ans);
        //Taille de la matrice F de levenshtein
        int fmatsize=0;
        double *fmat=NULL;
        fmatsize=maxlen+1;
        PROTECT(Fmat = allocVector(REALSXP, (fmatsize*fmatsize)));
        fmat=REAL(Fmat);
        for (i=0;i<alphasize;i++) {
          for(j=i; j<alphasize;j++){
            if (scost[TMRMATRIXINDEXC(i,j,alphasize)]>maxscost) {
              maxscost=scost[TMRMATRIXINDEXC(i,j,alphasize)];
            }
          }
        }
        maxscost=fmin(maxscost,2*indel);
        //Initialisation, peut �tre fait qu'une fois
        for (i=0;i<fmatsize;i++) {
          fmat[TMRMATRIXINDEXC(i,0,fmatsize)]=fmat[TMRMATRIXINDEXC(0,i,fmatsize)]=i*indel;
        }
        //Cout pour les diff�rentes possibilit�s
        double minimum=0, j_indel=0, sub=0;//, lenmax=0;
        //�tats compar�s
        int i_state, j_state, prefix;
        double maxpossiblecost=0;
        //starting store index
        //int i_start, j_start, i_end, j_end, i_index, j_index, base_index;
        //Pour chaque s�quence i
        for (is=0;is<nseq;is++) {
            //toutes les distances intra-groupes=0
            distmatrix[TMRMATRIXINDEXC(is,is,nseq)]=0;
            for (js=is+1;js<nseq;js++) {
                double cmpres=0;
                 ///Calcul des distances

                    //On passe les prefix commun
                    i=1;
                    j=1;
                    m=slen[is]+1;
                    n=slen[js]+1;
                    prefix=0;
                    while (i<m&&j<n&&sequences[TMRMATRIXINDEXC(is,i-1,nseq)]==sequences[TMRMATRIXINDEXC(js,j-1,nseq)]) {
                        i++;
                        j++;
                        prefix++;
                    }
                    //+1 pour correspondre � la matrice F
                    while (i<m) {
                        j=prefix+1;
                        while (j<n) {
                            i_state=sequences[TMRMATRIXINDEXC(is,i-1,nseq)];
                            j_state=sequences[TMRMATRIXINDEXC(js,j-1,nseq)];
                            if (i_state == j_state) {
                                cost = 0;
                            } else {
                                cost = scost[TMRMATRIXINDEXC(i_state,j_state,alphasize)];
                                //      				Rprintf("costs = %d %d, %d => %f \n",TMRMATRIXINDEXC(i_state,j_state,alphasize),i_state,j_state,cost);
                            }
                            minimum=fmat[TMRMATRIXINDEXC(i-prefix,j-1-prefix,fmatsize)]+ indel;

                            j_indel=fmat[TMRMATRIXINDEXC(i-1-prefix,j-prefix,fmatsize)]+ indel;
                            if (j_indel<minimum)minimum=j_indel;
                            sub=fmat[TMRMATRIXINDEXC(i-1-prefix,j-1-prefix,fmatsize)]+ cost;
                            if (sub<minimum)minimum=sub;
                            fmat[TMRMATRIXINDEXC(i-prefix,j-prefix,fmatsize)]=minimum;
                            j++;
                        }
                        i++;
                    }//Fmat build
                    m--;
                    n--;
                    //Warning! m and n decreased!!!!!
                    maxpossiblecost=abs(n-m)*indel+maxscost*fmin((double)m,(double)n);
                    cmpres=normalizeDistance(fmat[TMRMATRIXINDEXC(m-prefix,n-prefix,fmatsize)], maxpossiblecost, m, n, norm);

                
                //Same for j
                distmatrix[TMRMATRIXINDEXC(js,is,nseq)]=distmatrix[TMRMATRIXINDEXC(is,js,nseq)]=cmpres;


                //result[TMRMATRIXINDEXC(is,js,nseq)]=result[TMRMATRIXINDEXC(js,is,nseq)]=cmpres;
            }//end js
        }
        
        UNPROTECT(2);
        return ans;

    }

}
