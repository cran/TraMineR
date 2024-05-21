#include "OMPerdistanceII.h"
OMPerdistanceII::OMPerdistanceII(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS)
	:OMdistance(normS, Ssequences, seqdim, lenS){
}
OMPerdistanceII::OMPerdistanceII(OMPerdistanceII *dc)
	:OMdistance(dc), timecost(dc->timecost), seqdur(dc->seqdur), indellist(dc->indellist), seqlen(dc->seqlen){
}

void OMPerdistanceII::setParameters(SEXP params){
	OMdistance::setParameters(params);
	this->timecost = REAL(getListElement(params, "timecost"))[0];
	this->seqdur=REAL(getListElement(params, "seqdur"));
	this->indellist=REAL(getListElement(params, "indels"));
	this->seqlen=INTEGER(getListElement(params, "seqlength"));
	this->tokdeplist=REAL(getListElement(params, "tokdepcoeff"));
}

OMPerdistanceII::~OMPerdistanceII(){
}

double OMPerdistanceII::distance(const int&is, const int& js){
     //On passe les prefix commun
    double minimum=0, j_indel=0, sub=0;//, lenmax=0;
    //etats comparés
    int i_state, j_state;
    double maxpossiblecost;
    int i=1;
    int j=1;
    int m=slen[is];
    int n=slen[js];
	int nn=seqlen[is];
	int mm=seqlen[js];
	double nl, ml;
	//int minlen = imin2(m, n);
	int mSuf = m+1, nSuf = n+1;
    int prefix=0;
	
	//Skipping common prefix
	TMRLOG(6,"Skipping common prefix\n");
    //Skipping common suffix
	TMRLOG(6,"Skipping common suffix\n");
	TMRLOG(6,"Skipping common suffix\n");
	
	for(int ii=prefix+1; ii<mSuf; ii++){
		i_state=sequences[MINDICE(is, ii-1, nseq)];
		fmat[MINDICE(ii-prefix,0,fmatsize)] = fmat[MINDICE(ii-prefix-1,0,fmatsize)]+
				getIndel(MINDICE(is, ii-1, nseq), i_state);
		//mm = mm + seqdur[MINDICE(is, ii-1, nseq)];
	}

	for(int ii=prefix+1; ii<nSuf; ii++){
		j_state=sequences[MINDICE(js, ii-1, nseq)];
		fmat[MINDICE(0,ii-prefix,fmatsize)] = fmat[MINDICE(0,ii-prefix-1,fmatsize)]+
				getIndel(MINDICE(js, ii-1, nseq), j_state); 
		//nn = nn + seqdur[MINDICE(js, ii-1, nseq)];
	}
		TMRLOGMATRIX(10,  fmat, mSuf-prefix, nSuf-prefix, fmatsize);
	TMRLOG(5,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d, prefix=%d, fmatsize=%d\n", m, n, mSuf, nSuf, i, j, prefix, fmatsize);
    //+1 pour correspondre a la matrice F
	int i_state_indice=0;
	int j_state_indice=0;
	int fmat_ij_prefix=0;
    while (j<nSuf) {
        i=prefix+1;
		fmat_ij_prefix=1 + ((j-prefix)*fmatsize);
		j_state_indice = MINDICE(js,j-1,nseq);
		j_state=sequences[j_state_indice];
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
			//Adding i_indel
            minimum=fmat[fmat_ij_prefix-fmatsize]+getIndel(j_state_indice, j_state);
            //j_indel=fmat[MINDICE(i-1-prefix,j-prefix,fmatsize)]+ indel;
			//add j_indel
            j_indel=fmat[fmat_ij_prefix-1]+getIndel(i_state_indice, i_state);
            if (j_indel<minimum)minimum=j_indel;
			
			//////////////////////////////
            //Computing current indel cost
			//////////////////////////////
			//sub=fmat[MINDICE(i-1-prefix,j-1-prefix,fmatsize)]+ cost;
			//TMRLOG(6,"Substitution cost\n");
				//TMRLOG(6,"Sub cost\n");
            sub=fmat[fmat_ij_prefix-1-fmatsize]+ getSubCost(i_state, j_state, i_state_indice, j_state_indice);
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

    maxpossiblecost=abs(nn-mm)*indel+maxscost*fmin2((double)mm,(double)nn);
	
	TMRLOG(6,"End of dist compute index %d val %g\n", MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize), fmat[MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)]);
	if(MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)<0||MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)>fmatsize*fmatsize){
		TMRLOG(4,"End of dist compute index %d val %g\n", MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize), fmat[MINDICE(mSuf-1-prefix,nSuf-1-prefix,fmatsize)]);
		TMRLOG(4,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d, prefix=%d, fmatsize=%d\n", m, n, mSuf, nSuf, i, j, prefix, fmatsize);
		TMRLOG(4,"is =%d, js=%d\n", is, js);
	}
	TMRLOGMATRIX(10,  fmat, mSuf-prefix, nSuf-prefix, fmatsize);

	nl = double(nn) * indel;
	ml = double(mm) * indel;
    return normalizeDistance(fmat[MINDICE(mSuf-1-prefix, nSuf-1-prefix, fmatsize)], maxpossiblecost, ml, nl);
}

