#include<R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
//#include <math.h>

/**

*/

#define TMRMATRIXINDEXC(ligne, colone,len) (ligne)+(colone)*len
//#define TMRMIN(a,b) ((a)<(b))?a:b
#define TMRDISTINDEX(i,j,n) (n*(i-1) - i*(i-1)/2 + j-i-1)

static R_INLINE int distIndex(const int &i,const int &j,const int &n) {
    if (i<j)return TMRDISTINDEX(i,j,n);
    else return TMRDISTINDEX(j,i,n);
}

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
	return rawdist;
}

static R_INLINE void setDistance(const int &is,const int &js,const int* magicIndex, const int * magicSeq, const int& finalnseq, SEXP& ans, const double& cmpres) {
    int j_start=magicIndex[js];
    int j_end=magicIndex[js+1];
    int i_start=magicIndex[is];
    int i_end=magicIndex[is+1];
    int i_index, j_index, i, j, base_index;
    double *result=REAL(ans);
    for (i=i_start;i<i_end;i++) {
        i_index=magicSeq[i];
        for (j=j_start;j<j_end;j++) {
            j_index=magicSeq[j];
            if (i_index!=j_index) {
                base_index=distIndex(i_index,j_index,finalnseq);
                result[base_index]=cmpres;
            }
        }
    }
}


extern "C" {

    SEXP cstringdistance(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP indelS, SEXP alphasizeS, SEXP costsS, SEXP normS, SEXP magicIndexS, SEXP magicSeqS, SEXP disttypeS) {
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

        int disttype=INTEGER(disttypeS)[0];
        int* magicIndex=INTEGER(magicIndexS);
        int* magicSeq=INTEGER(magicSeqS);
        int finalnseq=length(magicSeqS);

        //Alocation du vecteur de distance
        //REprintf("Final seq %d\n",finalnseq);
        PROTECT(ans = allocVector(REALSXP, (finalnseq*(finalnseq-1)/2)));

        //Taille de la matrice F de levenshtein
        int fmatsize=0;
        double *fmat=NULL;
        if (disttype==1) {
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
            maxscost=fmin2(maxscost,2*indel);
            //Initialisation, peut �tre fait qu'une fois
            for (i=0;i<fmatsize;i++) {
                fmat[TMRMATRIXINDEXC(i,0,fmatsize)]=fmat[TMRMATRIXINDEXC(0,i,fmatsize)]=i*indel;
            }
        }
        //Cout pour les diff�rentes possibilit�s
        double minimum=0, j_indel=0, sub=0;//, lenmax=0;
        //�tats compar�s
        int i_state, j_state, prefix, sign;
        if (disttype==2) {
            sign=1;
        } else if (disttype==3) {
            sign=-1;
        }
        double maxpossiblecost=0;
        //starting store index
        //int i_start, j_start, i_end, j_end, i_index, j_index, base_index;
        //Pour chaque s�quence i

        for (is=0;is<nseq;is++) {
            //toutes les distances intra-groupes=0
            setDistance(is,is,magicIndex,magicSeq, finalnseq, ans, 0);
            for (js=is+1;js<nseq;js++) {
                double cmpres=0;
                if (disttype==1) { //optimal matching
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
                    maxpossiblecost=abs(n-m)*indel+maxscost*fmin2((double)m,(double)n);
                    cmpres=normalizeDistance(fmat[TMRMATRIXINDEXC(m-prefix,n-prefix,fmatsize)], maxpossiblecost, m, n, norm);

                } else if (disttype>1) { //BEGIN LCP
                    m=slen[is];
                    n=slen[js];
                    // Computing min length
                    if (n<m) minimum = n;
                    else minimum = m;

                    if (disttype==2) {
                        i=0;
                        while (sequences[TMRMATRIXINDEXC(is,i,nseq)]==sequences[TMRMATRIXINDEXC(js,i,nseq)] && i<minimum) {
                            i++;
                        }
                    } else {
                        i=1;
                        while (sequences[TMRMATRIXINDEXC(is,(m-i),nseq)]==sequences[TMRMATRIXINDEXC(js,(n-i),nseq)] && i<=minimum) {
                            i++;
                        }
                        i--;
                    }
                    cmpres=normalizeDistance((double)n+(double)m-2.0*(double)i, (double)n+(double)m, m, n, norm);
                }//End LCP
//		  Rprintf("cmpres = %d %d => %f \n",(1+is),(1+js), cmpres);
                //return Fmat;

                //Same for j
                setDistance(is,js,magicIndex,magicSeq, finalnseq, ans, cmpres);


                //result[TMRMATRIXINDEXC(is,js,nseq)]=result[TMRMATRIXINDEXC(js,is,nseq)]=cmpres;
            }//end js
        }
        if (disttype==1) {
            UNPROTECT(1);
        }
        UNPROTECT(1);
        return ans;

    }

}
