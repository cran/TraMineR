#include "TraMineR.h"

/**

*/

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
	case 4:
		if (maxdist==0)return 1;
        return (2*rawdist)/(rawdist+maxdist);
		
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
static R_INLINE double OMdistanceOpt(int* slen,const int &is,const int &js, const int&nseq, int* sequences, const int &alphasize, double * scost, double * fmat, const int& fmatsize, const double& indel, const double& maxscost, const int&norm) {

    //On passe les prefix commun
    double minimum=0, j_indel=0, sub=0;//, lenmax=0;
    //etats comparés
    int i_state, j_state;
    double maxpossiblecost;
    int i=1;
    int j=1;
    int m=slen[is];
    int n=slen[js];
	//int minlen = imin2(m, n);
	int mSuf = m+1, nSuf = n+1;
    int prefix=0;
	//Skipping common prefix
	//TMRLOG(6,"Skipping common prefix\n");
    while (i<mSuf && j<nSuf && sequences[MINDICE(is,i-1,nseq)]==sequences[MINDICE(js,j-1,nseq)]) {
        i++;
        j++;
        prefix++;
    }
	//Skipping common suffix
	//TMRLOG(6,"Skipping common suffix\n");
	while (mSuf>i && nSuf>j && sequences[MINDICE(is,(mSuf-2),nseq)]==sequences[MINDICE(js,(nSuf-2),nseq)]) {
		mSuf--;
		nSuf--;
    }
	//TMRLOG(6,"Skipping common suffix\n");
	
	//TMRLOG(5,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d, prefix=%d, fmatsize=%d\n", m, n, mSuf, nSuf, i, j, prefix, fmatsize);
    //+1 pour correspondre a la matrice F
	int fmat_ij_prefix=0;
	int i_state_indice=0;
    while (j<nSuf) {
        i=prefix+1;
		fmat_ij_prefix=1 + ((j-prefix)*fmatsize);
		j_state=sequences[MINDICE(js,j-1,nseq)];
		i_state_indice=is+prefix*nseq;
        while (i<mSuf) {
            //i_state=sequences[MINDICE(is,i-1,nseq)];
			//TMRLOG(6,"Getting i state\n");
            i_state=sequences[i_state_indice];
			//////////////////////////////
            //Computing current indel cost
			//////////////////////////////
			//fmat_ij_prefix=((i-prefix)+(j-prefix)*(fmatsize));
			//minimum=fmat[MINDICE(i-prefix,j-1-prefix,fmatsize)]+ indel;
			//TMRLOG(6,"fmat_ij_prefix =%d,th =%d \n", fmat_ij_prefix, (MINDICE(i-prefix,j-prefix,fmatsize)));
            minimum=fmat[fmat_ij_prefix-fmatsize];
            //j_indel=fmat[MINDICE(i-1-prefix,j-prefix,fmatsize)]+ indel;
            j_indel=fmat[fmat_ij_prefix-1];
            if (j_indel<minimum)minimum=j_indel;
			minimum+=indel;
			
			//////////////////////////////
            //Computing current indel cost
			//////////////////////////////
			//sub=fmat[MINDICE(i-1-prefix,j-1-prefix,fmatsize)]+ cost;
			//TMRLOG(6,"Substitution cost\n");
			if (i_state == j_state) {
                sub=fmat[fmat_ij_prefix-1-fmatsize];
            } else {
				//TMRLOG(6,"Sub cost\n");
                sub=fmat[fmat_ij_prefix-1-fmatsize]+ scost[MINDICE(i_state,j_state,alphasize)];
            }
            //sub=fmat[fmat_ij_prefix-1-fmatsize]+ cost;
            if (sub<minimum) {
				fmat[fmat_ij_prefix]=sub;
			} else {
				fmat[fmat_ij_prefix]=minimum;
			}
            //fmat[MINDICE(i-prefix,j-prefix,fmatsize)]=minimum;
            i++;
			fmat_ij_prefix++;
			i_state_indice+=nseq;
        }
		j++;
    }//Fmat build
	//Max possible cost
    maxpossiblecost=abs(n-m)*indel+maxscost*fmin2((double)m,(double)n);
	
	//TMRLOG(6,"End of dist compute\n");
    return normalizeDistance(fmat[MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)], maxpossiblecost, m, n, norm);
}

static R_INLINE double OMdistance(int* slen,const int &is,const int &js, const int&nseq, int* sequences, const int &alphasize, double * scost, double * fmat, const int& fmatsize, const double& indel, const double& maxscost, const int&norm) {

    //On passe les prefix commun
    double minimum=0, j_indel=0, sub=0;//, lenmax=0;
    //etats comparés
    int i_state, j_state;
    double cost, maxpossiblecost;
    int i=1;
    int j=1;
    int m=slen[is]+1;
    int n=slen[js]+1;
    int prefix=0;
    while (i<m&&j<n&&sequences[MINDICE(is,i-1,nseq)]==sequences[MINDICE(js,j-1,nseq)]) {
        i++;
        j++;
        prefix++;
    }
    //+1 pour correspondre � la matrice F
    while (i<m) {
        j=prefix+1;
        while (j<n) {
            i_state=sequences[MINDICE(is,i-1,nseq)];
            j_state=sequences[MINDICE(js,j-1,nseq)];
            if (i_state == j_state) {
                cost = 0;
            } else {
                cost = scost[MINDICE(i_state,j_state,alphasize)];
                //      				Rprintf("costs = %d %d, %d => %f \n",MINDICE(i_state,j_state,alphasize),i_state,j_state,cost);
            }
            minimum=fmat[MINDICE(i-prefix,j-1-prefix,fmatsize)]+ indel;

            j_indel=fmat[MINDICE(i-1-prefix,j-prefix,fmatsize)]+ indel;
            if (j_indel<minimum)minimum=j_indel;
            sub=fmat[MINDICE(i-1-prefix,j-1-prefix,fmatsize)]+ cost;
            if (sub<minimum)minimum=sub;
            fmat[MINDICE(i-prefix,j-prefix,fmatsize)]=minimum;
            j++;
        }
        i++;
    }//Fmat build
    m--;
    n--;
    //Warning! m and n decreased!!!!!
    maxpossiblecost=abs(n-m)*indel+maxscost*fmin2((double)m,(double)n);
    return normalizeDistance(fmat[MINDICE(m-prefix,n-prefix,fmatsize)], maxpossiblecost, m, n, norm);
}
static R_INLINE double LCPdistance(int* slen,const int &is,const int &js, const int&nseq, int* sequences, const int&norm, const int& sign) {

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
    return normalizeDistance((double)n+(double)m-2.0*(double)i, (double)n+(double)m, m, n, norm);
}

static R_INLINE double DHDdistance(int* slen,const int &is,const int &js, const int&nseq, int* sequences, const int&norm,  const int &alphasize, double * scost, const double& maxscost) {

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
    return normalizeDistance(cost, maxscost, m, n, norm);
}


extern "C" {

    SEXP cstringdistance(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP indelS, SEXP alphasizeS, SEXP costsS, SEXP normS, SEXP magicIndexS, SEXP magicSeqS, SEXP disttypeS) {
        //Objet R, matrice des distances (objet dist)
        SEXP ans;//, Fmat;
        //Indices, avec s pour séquences
        int i, j, is, js;
        //longueur des séquences m=i, n=j
        //int m, n;
        //normalisation?
        int norm=INTEGER(normS)[0];
        //Nombre de séquence
        int nseq=INTEGER(seqdim)[0];
        //nb colonnes des séquence
        int maxlen=INTEGER(seqdim)[1];
        //Matrice des séquence
        int* sequences=INTEGER(Ssequences);
        //Tailles des séquence
        int* slen=INTEGER(lenS);
        //indel
        double indel=REAL(indelS)[0];
        //nb etats
        int alphasize=INTEGER(alphasizeS)[0];
        //Matrice des couts de substitutions
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
        if (disttype<=1) {
            fmatsize=maxlen+1;
            //PROTECT(Fmat = allocVector(REALSXP, (fmatsize*fmatsize)));
            //fmat= new double[fmatsize*fmatsize];
			fmat= (double*) R_alloc(fmatsize*fmatsize,sizeof(double));
            for (i=0;i<alphasize;i++) {
                for (j=i; j<alphasize;j++) {
                    if (scost[MINDICE(i,j,alphasize)]>maxscost) {
                        maxscost=scost[MINDICE(i,j,alphasize)];
                    }
                }
            }
            maxscost=fmin2(maxscost,2*indel);
			if(norm==4) { //Maximum cost is thus defined as in Yujian and Bo (2007)
				maxscost=2*indel;
			}
            //Initialisation, peut etre fait qu'une fois
            for (i=0;i<fmatsize;i++) {
                fmat[MINDICE(i,0,fmatsize)]=fmat[MINDICE(0,i,fmatsize)]=i*indel;
            }
        }
        if (disttype == 4) {
            maxscost = indel;
        }
        //Cout pour les diff�rentes possibilit�s
        //double minimum=0, j_indel=0, sub=0;//, lenmax=0;
        //etats comparés
        int sign = 1;
        if (disttype==3) {
            sign=-1;
        }
        //starting store index
        //int i_start, j_start, i_end, j_end, i_index, j_index, base_index;
        //Pour chaque s�quence i

        for (is=0;is<nseq;is++) {
            //toutes les distances intra-groupes=0
			R_CheckUserInterrupt();
            setDistance(is,is,magicIndex,magicSeq, finalnseq, ans, 0);
            for (js=is+1;js<nseq;js++) {
                double cmpres=0;
				if (disttype==0) {
					cmpres=OMdistanceOpt(slen,is,js, nseq, sequences, alphasize, scost, fmat, fmatsize, indel,  maxscost, norm);
				} else if (disttype==1) { //optimal matching
                    ///Calcul des distances
                    cmpres=OMdistance(slen,is,js, nseq, sequences, alphasize, scost, fmat, fmatsize, indel,  maxscost, norm);

                } else if (disttype==4) { //DHD
                    cmpres=DHDdistance(slen,is,js,nseq, sequences, norm, alphasize, scost, maxscost);
                } else if (disttype>1) { //BEGIN LCP
                    cmpres=LCPdistance(slen,is,js, nseq, sequences, norm, sign);
                }//End LCP
                TMRLOG(5,"cmpres = %d %d => %f \n",(1+is),(1+js), cmpres);
                //return Fmat;

                //Same for j
                setDistance(is,js,magicIndex,magicSeq, finalnseq, ans, cmpres);


                //result[MINDICE(is,js,nseq)]=result[MINDICE(js,is,nseq)]=cmpres;
            }//end js
        }
        UNPROTECT(1);
        return ans;

    }

}

