#include "OMVI2distance.h"

OMVI2distance::OMVI2distance(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS)
	:OMdistance(normS, Ssequences, seqdim, lenS){
		TMRLOG(3,"fmatsize =%d\n", fmatsize);
		this->opmat = new int[fmatsize*fmatsize];
		this->prev_j_state = new int[fmatsize*fmatsize];
		this->prev_i_state = new int[fmatsize*fmatsize];
}

OMVI2distance::OMVI2distance(OMVI2distance *dc)
  :OMdistance(dc), timecost(dc->timecost), localcost(dc->localcost) {
    TMRLOG(3,"fmatsize =%d\n", fmatsize);
    this->opmat = new int[fmatsize*fmatsize];
    this->prev_j_state = new int[fmatsize*fmatsize];
    this->prev_i_state = new int[fmatsize*fmatsize];
}

void OMVI2distance::setParameters(SEXP params) {
  OMdistance::setParameters(params);
  this->timecost = REAL(getListElement(params, "timecost"))[0] * maxscost;
  this->localcost = REAL(getListElement(params, "localcost"))[0];
  TMRLOG(7, "timecost=%g, localcost=%g\n", this->timecost, this->localcost);
}

OMVI2distance::~OMVI2distance(){
	delete [] this->opmat;
	delete [] prev_j_state;
	delete [] prev_i_state;
}


double OMVI2distance::distance(const int&is, const int& js){

     //On passe les prefix commun
    double minimum=0, j_indel=0, sub=0;//, lenmax=0;
    //etats comparés
    double maxpossiblecost;
    int i=1;
    int j=1;
    int m=slen[is];
    int n=slen[js];
	//int minlen = imin2(m, n);
	int mSuf = m+1, nSuf = n+1;
	//Skipping common prefix
	/* TMRLOG(6,"Skipping common prefix\n");
    while (i<m && j<n && sequences[MINDICE(is,i-1,nseq)]==sequences[MINDICE(js,j-1,nseq)]) {
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
	 */
	TMRLOG(5,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d, fmatsize=%d\n", m, n, mSuf, nSuf, i, j, fmatsize);
    //+1 pour correspondre a la matrice F
	int fmat_ij_prefix=0;
	int i_state_indice=0;
	int  j_state, i_state;
	fmat[0] = 0;

	this->prev_j_state[0]=sequences[MINDICE(js, 0, nseq)];
	this->prev_i_state[0]=sequences[MINDICE(is, 0, nseq)];
	//Opmat operation -1=subst; 1 = insert state from i in j; 2 = insert state from j in i
	opmat[0] =0;
	int prev= 0;
	for(int ii=1; ii<mSuf; ii++){
		i_state=sequences[MINDICE(is, ii-1, nseq)];
		fmat_ij_prefix=MINDICE(ii,0,fmatsize);
		prev = imax2(fmat_ij_prefix-1, 0);
		opmat[fmat_ij_prefix]=1;
		TMRLOG(6,"i_state =%d, fmat_ij_prefix=%d, prev=%d, pi=%d, pj=%d\n", i_state, fmat_ij_prefix, prev, this->prev_i_state[prev], this->prev_j_state[prev]);
		fmat[fmat_ij_prefix] = fmat[prev]+
				getFirstIndel(i_state, this->prev_j_state[prev], ii);

		this->prev_j_state[fmat_ij_prefix] = i_state; //i was inserted in j
		this->prev_i_state[fmat_ij_prefix] = i_state; //i was inserted in j
	}
	TMRLOGMATRIX(10,  fmat, mSuf, nSuf, fmatsize);

	for(int ii=1; ii<nSuf; ii++){
		j_state=sequences[MINDICE(js, ii-1, nseq)];
		fmat_ij_prefix=MINDICE(0, ii, fmatsize);
		prev = imax2(fmat_ij_prefix-fmatsize, 0);
		opmat[fmat_ij_prefix]=2;
		TMRLOG(6,"j_state =%d, fmat_ij_prefix=%d, prev=%d, pi=%d, pj=%d\n", j_state, fmat_ij_prefix, prev, this->prev_i_state[prev], this->prev_j_state[prev]);
		fmat[fmat_ij_prefix] = fmat[prev]+
				getFirstIndel(j_state, this->prev_i_state[prev], ii);
		this->prev_j_state[fmat_ij_prefix]=j_state; //j was inserted in i
		this->prev_i_state[fmat_ij_prefix]=j_state; //j was inserted in i
		//prev_istate=i_state;
	}
	TMRLOG(5,"Fmat initialized\n");
	TMRLOGMATRIX(10,  fmat, mSuf, nSuf, fmatsize);
	/* int * i_states_table= new int[fmatsize*fmatsize];
	int * j_states_table= new int[fmatsize*fmatsize];
	double * jindel_table= new double[fmatsize*fmatsize];
	double * iindel_table= new double[fmatsize*fmatsize];
	double * sub_table= new double[fmatsize*fmatsize]; */
	int diagonal = fmatsize+1;
    while (j<nSuf) {
        i=1;
		fmat_ij_prefix=1 + (j*fmatsize);
		j_state=sequences[MINDICE(js, j-1, nseq)];
		//j_indel_val = getIndel(j_state, sequences[MINDICE(js, (j==1?1:(j-2)),nseq)], sequences[MINDICE(js,(j==n?(n-2):j), nseq)]);
		//next_jstate= sequences[MINDICE(js,(j==n?(n-2):j), nseq)];
		i_state_indice= is;
		//prev_istate =sequences[MINDICE(is, firststate,nseq)];
        while (i<mSuf) {
            //i_state=sequences[MINDICE(is,i-1,nseq)];
			//TMRLOG(6,"Getting i state\n");
            i_state=sequences[i_state_indice];
			//////////////////////////////
            //Computing current indel cost
			//////////////////////////////
			//i_states_table[fmat_ij_prefix]=i_state;
			//j_states_table[fmat_ij_prefix]=j_state;
			//fmat_ij_prefix=((i-prefix)+(j-prefix)*(fmatsize));
			//minimum=fmat[MINDICE(i-prefix,j-1-prefix,fmatsize)]+ indel;
			//TMRLOG(6,"fmat_ij_prefix =%d,th =%d \n", fmat_ij_prefix, (MINDICE(i-prefix,j-prefix,fmatsize)));
			TMRLOG(6,"j_state =%d, fmat_ij_prefix=%d, prev=%d, pi=%d, pj=%d\n", j_state, fmat_ij_prefix, prev, this->prev_i_state[fmat_ij_prefix-1], this->prev_j_state[fmat_ij_prefix-1]);
			minimum=fmat[fmat_ij_prefix-1]+getIndel(i_state, fmat_ij_prefix-1);
			//iindel_table[fmat_ij_prefix] = minimum;
			TMRLOG(6,"j_state =%d, fmat_ij_prefix=%d, prev=%d, pi=%d, pj=%d\n", j_state, fmat_ij_prefix, prev, this->prev_i_state[fmat_ij_prefix-fmatsize], this->prev_j_state[fmat_ij_prefix-fmatsize]);
			j_indel=fmat[fmat_ij_prefix-fmatsize]+getIndel(j_state, fmat_ij_prefix-fmatsize);
			//jindel_table[fmat_ij_prefix]=j_indel;
            //j_indel=fmat[fmat_ij_prefix-fmatsize]+j_indel_val;
			double diffii =fabs(minimum - j_indel);
			if(fabs(minimum - j_indel)<OMVI2TOL){
				opmat[fmat_ij_prefix]=-2;
				prev_i_state[fmat_ij_prefix] = i_state;
				prev_j_state[fmat_ij_prefix] = j_state;
			}else if(minimum > j_indel){
				minimum=j_indel;
				//j state inserted in i
				opmat[fmat_ij_prefix]=2;
				prev_i_state[fmat_ij_prefix] = j_state;
				prev_j_state[fmat_ij_prefix] = j_state;
			} else {
				//i state inserted in j
				opmat[fmat_ij_prefix]=1;
				prev_i_state[fmat_ij_prefix] = i_state;
				prev_j_state[fmat_ij_prefix] = i_state;
			}

			//////////////////////////////
            //Computing current indel cost
			//////////////////////////////
			//sub=fmat[MINDICE(i-1-prefix,j-1-prefix,fmatsize)]+ cost;
			//TMRLOG(6,"Substitution cost\n");
			if (i_state == j_state) {
                sub=fmat[fmat_ij_prefix-diagonal];
            } else {
				//TMRLOG(6,"Sub cost\n");
				int previous_ij =fmat_ij_prefix-diagonal;
				if(opmat[previous_ij]<0 && prev_i_state[previous_ij]==i_state && prev_j_state[previous_ij] == j_state){
					sub=fmat[previous_ij] + 2 * timecost * scost[MINDICE(i_state, j_state, alphasize)];
				}
                else {
					sub= fmat[previous_ij] + scost[MINDICE(i_state, j_state, alphasize)];
				}
            }
				//sub_table[fmat_ij_prefix]=fabs(sub-minimum);
            //sub=fmat[fmat_ij_prefix-1-fmatsize]+ cost;
            if (sub-minimum < OMVI2TOL) {
				fmat[fmat_ij_prefix]=sub;
				//substitution
				opmat[fmat_ij_prefix]=-1;
				prev_i_state[fmat_ij_prefix] = i_state;
				prev_j_state[fmat_ij_prefix] = j_state;
			} else {
				fmat[fmat_ij_prefix]=minimum;
			}
			TMRLOG(5, "fmat_ij_prefix=%d, diffii=%g, diffsm=%g, sub=%g, minimum=%g, choosing=%d", fmat_ij_prefix, diffii, fabs(sub-minimum), sub, minimum, opmat[fmat_ij_prefix]);
            //fmat[MINDICE(i-prefix,j-prefix,fmatsize)]=minimum;
            i++;
			fmat_ij_prefix++;
			i_state_indice+=nseq;

        }
		j++;
    }//Fmat build
	//Max possible cost
    maxpossiblecost=abs(n-m)*indel+maxscost*fmin2((double)m,(double)n);

	TMRLOG(2,"End of dist compute index %d val %g\n", MINDICE(mSuf-1, nSuf-1, fmatsize), fmat[MINDICE(mSuf-1, nSuf-1, fmatsize)]);
	if(MINDICE(mSuf-1, nSuf-1, fmatsize)<0 || MINDICE(mSuf-1, nSuf-1, fmatsize) > fmatsize*fmatsize){
		TMRLOG(4,"End of dist compute index %d val %g\n", MINDICE(mSuf-1, nSuf-1, fmatsize), fmat[MINDICE(mSuf-1, nSuf-1, fmatsize)]);
		TMRLOG(4,"m =%d, n=%d, mSuf=%d, nSuf=%d i=%d, j=%d,  fmatsize=%d\n", m, n, mSuf, nSuf, i, j,  fmatsize);
		TMRLOG(4,"is =%d, js=%d\n", is, js);
	}
	TMRLOGMATRIX(2,  fmat, mSuf, nSuf, fmatsize);
	TMRLOGMATRIXINT(2,  prev_i_state, mSuf, nSuf, fmatsize);
	TMRLOGMATRIXINT(2,  prev_j_state, mSuf, nSuf, fmatsize);
	TMRLOGMATRIXINT(2,  opmat, mSuf, nSuf, fmatsize);
	//TMRLOGMATRIXINT(2,  i_states_table, mSuf-prefix, nSuf-prefix, fmatsize);
	//TMRLOGMATRIXINT(2,  j_states_table, mSuf-prefix, nSuf-prefix, fmatsize);
	//TMRLOGMATRIX(2,  jindel_table, mSuf-prefix, nSuf-prefix, fmatsize);
	//TMRLOGMATRIX(2,  iindel_table, mSuf-prefix, nSuf-prefix, fmatsize);
	//TMRLOGMATRIX(2,  sub_table, mSuf-prefix, nSuf-prefix, fmatsize);

	// delete[] i_states_table;
	// delete[] j_states_table;
	// delete[] jindel_table;
	// delete[] iindel_table;
	// delete[] sub_table;

    return normalizeDistance(fmat[MINDICE(mSuf-1, nSuf-1, fmatsize)], maxpossiblecost, m, n);
}
