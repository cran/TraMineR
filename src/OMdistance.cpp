#include "OMdistance.h"
OMdistance::OMdistance(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS)
	:DistanceCalculator(normS, Ssequences, seqdim, lenS), 
	scost(NULL), alphasize(0), indel(0), maxscost(0){
		fmatsize=maxlen+1;
		this->fmat = new double[fmatsize*fmatsize];
}
OMdistance::OMdistance(OMdistance *dc)
	:DistanceCalculator(dc), 
	scost(dc->scost), alphasize(dc->alphasize), indel(dc->indel), maxscost(dc->maxscost){
	fmatsize=maxlen+1;
	this->fmat = new double[fmatsize*fmatsize];
	for (int i=0;i<fmatsize;i++) {
		fmat[MINDICE(i,0,fmatsize)]=fmat[MINDICE(0,i,fmatsize)]=i*indel;
    }
}

void OMdistance::setParameters(SEXP params){
	scost = REAL(getListElement(params, "scost"));
	alphasize = INTEGER(getListElement(params, "alphasize"))[0];
	TMRLOGMATRIX(10,  scost,alphasize, alphasize, alphasize);
	indel = REAL(getListElement(params, "indel"))[0];
	if(norm==4) { //Maximum cost is thus defined as in Yujian and Bo (2007)
		maxscost=2*indel;
	} else {
		for (int i=0;i<alphasize;i++) {
			for (int j=i+1; j<alphasize;j++) {
				if (scost[MINDICE(i,j,alphasize)]>maxscost) {
					maxscost=scost[MINDICE(i,j,alphasize)];
				}
			}
		}
		maxscost=fmin2(maxscost,2*indel);
	}
	//fmat= (double*) R_alloc(fmatsize*fmatsize,sizeof(double));
            //Initialisation, peut etre fait qu'une fois
    for (int i=0;i<fmatsize;i++) {
		fmat[MINDICE(i,0,fmatsize)]=fmat[MINDICE(0,i,fmatsize)]=i*indel;
    }
}

OMdistance::~OMdistance(){
	delete [] fmat;
}
double OMdistance::distance(const int&is, const int& js){
     //On passe les prefix commun
    double minimum=0, j_indel=0, sub=0;//, lenmax=0;
    //etats comparés
    int i_state, j_state;
    double maxpossiblecost;
    int i=1;
    int j=1;
    int m=slen[is];
    int n=slen[js];
	double ml, nl;
	//int minlen = imin2(m, n);
	int mSuf = m+1, nSuf = n+1;
    int prefix=0;
	//Skipping common prefix
	TMRLOG(6,"Skipping common prefix\n");
    while (i<mSuf && j<nSuf && sequences[MINDICE(is,i-1,nseq)]==sequences[MINDICE(js,j-1,nseq)]) {
        i++;
        j++;
        prefix++;
    }
	//Skipping common suffix
	TMRLOG(6,"Skipping common suffix\n");
	while (mSuf>i && nSuf>j && sequences[MINDICE(is,(mSuf-2),nseq)]==sequences[MINDICE(js,(nSuf-2),nseq)]) {
		mSuf--;
		nSuf--;
    }
	TMRLOG(6,"Skipping common suffix\n");
	
	TMRLOG(5,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d, prefix=%d, fmatsize=%d\n", m, n, mSuf, nSuf, i, j, prefix, fmatsize);
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
	
	TMRLOG(6,"End of dist compute index %d val %g\n", MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize), fmat[MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)]);
	if(MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)<0||MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)>fmatsize*fmatsize){
		TMRLOG(4,"End of dist compute index %d val %g\n", MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize), fmat[MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)]);
		TMRLOG(4,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d, prefix=%d, fmatsize=%d\n", m, n, mSuf, nSuf, i, j, prefix, fmatsize);
		TMRLOG(4,"is =%d, js=%d\n", is, js);
	}
	TMRLOGMATRIX(10,  fmat, mSuf-1-prefix,nSuf-1-prefix,fmatsize);
	nl = double(n) * indel;
	ml = double(m) * indel;
    return normalizeDistance(fmat[MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)], maxpossiblecost, ml, nl);
}

