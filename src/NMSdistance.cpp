#include "NMSdistance.h"

	
void SUBSEQdistance::setParameters(SEXP params){
	kweights = REAL(getListElement(params, "kweights"));
	distMethod = INTEGER(getListElement(params, "distMethod"))[0];
	distTransform = INTEGER(getListElement(params, "distTransform"))[0];
	for (int is=0;is<nseq;is++) {
		this->resetKvect();
		if(slen[is]>0){
			this->computeattr(is,is);
		}
		for(int i =0; i<maxlen; i++){
			selfmatvect[MINDICE(is,i,nseq)]=kvect[i];
		}
	}
}

double SUBSEQdistance::distance(const int&is, const int& js){
	/*int m=slen[is];
    int n=slen[js];
	int minimum = m;
    if (n<m) {
		minimum = n;
	}*/
	int minimum = maxlen;
	this->resetKvect();
	if(slen[is]>0&&slen[js]>0){
		this->computeattr(is, js);
	}
	if(this->distMethod==1){
		double s=0, kval=0, ktot=0;
		for(int i =0; i<minimum; i++){
			if(kweights[i]!=0){
			//Take care of potentially big numbers
				kval=kvect[i]/(sqrt(selfmatvect[MINDICE(is, i, nseq)]));
				kval/=sqrt(selfmatvect[MINDICE(js, i, nseq)]);
				kval*=kweights[i];
				TMRLOG(5,"kval(%d)=%g\n", i, kvect[i]);
				ktot+=kweights[i];
				s+=kval;
			}
		}
		TMRLOG(5,"s tot:%f\n",s);
		return (1.0-s/ktot);
	} else {
		double Aval=0, Ival=0, Jval=0, ktot=0, dist, maxdist;
		for(int i =0; i<minimum; i++){
			if(kweights[i]!=0){
				TMRLOG(5,"kval(%d)=%g\n", i, kvect[i]);
				Aval+=kweights[i]*kvect[i];
				Ival+=kweights[i]*selfmatvect[MINDICE(is, i, nseq)];
				TMRLOG(5,"ival(%d)=%g\n", i, selfmatvect[MINDICE(is, i, nseq)]);
				Jval+=kweights[i]*selfmatvect[MINDICE(js, i, nseq)];
				TMRLOG(5,"Jval(%d)=%g\n", i, selfmatvect[MINDICE(js, i, nseq)]);
				ktot+=kweights[i];
			}
		}
		if(distTransform){
			Ival = log1p(Ival);
			Jval = log1p(Jval);
			Aval = log1p(Aval);
		}
		dist = Ival+Jval- 2*Aval;
		maxdist=Ival+Jval;
		TMRLOG(2,"Coord (is=%d, js=%d) dist %g, maxdist %g, Aval %g, Ival %g, Jval%g\n", is, js, dist, maxdist, Aval, Ival, Jval);
		if(this->norm==4){
			return normalizeDistance(sqrt(dist), sqrt(maxdist), Jval, Ival);
		}
		else{
			return normalizeDistance(dist, maxdist, Jval, Ival);
		}
	}
	return 0;
}

NMSdistance::NMSdistance(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS)
	:SUBSEQdistance(normS, Ssequences, seqdim, lenS){
		TMRLOG(5, "Requesting zmat\n");
		zmatsize=maxlen*maxlen;
		zmat= new int[zmatsize*2];
		TMRLOG(5, "Requesting fmat\n");
		hmat= new double[maxlen*maxlen];
		TMRLOG(5, "Requesting vmat\n");
		vmat=  new double[maxlen*maxlen];
	}

NMSdistance::NMSdistance(NMSdistance *dc)
	:SUBSEQdistance(dc){
		TMRLOG(5, "Requesting zmat\n");
		zmatsize=maxlen*maxlen;
		zmat= new int[zmatsize*2];
		TMRLOG(5, "Requesting fmat\n");
		hmat= new double[maxlen*maxlen];
		TMRLOG(5, "Requesting vmat\n");
		vmat=  new double[maxlen*maxlen];
	}

NMSdistance::~NMSdistance(){
	delete[] zmat;
	delete[] hmat;
	delete[] vmat;
}



void  NMSdistance::computeattr(const int&is, const int& js){
	int m=slen[is];
    int n=slen[js];
	int ksize = imin2(m, n);
	int i, j, ij, zsize=0, k=0, zi, zj, zindex, base_index, maxlen1=maxlen+1, iseq;
	double htot=0;
	//building Zmat & hmat
	TMRLOG(5,"\n\nCOMPUTING BETWEEN %d %d\n", is, js);
	for(i=0; i<m; i++){
		iseq=sequences[MINDICE(is,i,nseq)];
		for(j=0; j<n; j++){
			if(iseq==sequences[MINDICE(js,j,nseq)]){
				zmat[MINDICE(zsize,0,zmatsize)] =i;
				zmat[MINDICE(zsize,1,zmatsize)] =j;
				zsize++;
			}
		}
	}
	TMRLOG(5,"\n\nZmat\n");
	TMRLOGMATRIX(5, zmat, zsize,2, zmatsize);
	/// We want to ensure that last column and row are=0
	//Last value of vmat
	zindex= MINDICE(m-1,n-1,maxlen);
	// Column
	ij = MINDICE(m-1,0,maxlen);
	while(ij<=zindex) {
		vmat[ij]=0;
		ij+=maxlen;
	}// Column
	ij = MINDICE(0,n-1,maxlen);
	while(ij<zindex) {
		vmat[ij]=0;
		ij++;
	}
	
	zindex=0;
	zi=MINDICE(zindex,0,zmatsize);
	zj=MINDICE(zindex,1,zmatsize);
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			ij = MINDICE(i,j,maxlen);
			if(zindex<zsize && zmat[zi]==i && zmat[zj] == j){
				hmat[ij] = 1;
				zindex++;
				zi=MINDICE(zindex,0,zmatsize);
				zj=MINDICE(zindex,1,zmatsize);
				htot += 1;
			}else{
				hmat[ij] = 0;
			}
			vmat[ij]=0;
		}
	}
	TMRLOG(5,"\n\nH%d (htot=%g):\n", k, htot);
	TMRLOGMATRIX(5, hmat, m,n, maxlen);
	TMRLOG(5,"\n\nSetting kvect %d to %g\n", k, htot);
	this->kvect[k] = htot;
	k++;
	TMRLOG(5,"\n\nPreEntering %d\n", k);
	if(m>1&&n>1){
		while(k<ksize && htot>0){
			TMRLOG(5,"\n\nEntering %d\n", k);
			//Building hmat
			if(htot==DBL_MAX){
				error(" [!] Number of subsequences is getting too big");
			} 
			
			//Building vmat
			//Building from reverse to go faster
			for(j=n-2; j>=0; j--){
				base_index = j*maxlen;
				ij = base_index+m-2;
				while(ij>=base_index){
					vmat[ij]=vmat[ij+1]+vmat[ij+maxlen]-vmat[ij+maxlen1]+ hmat[ij+maxlen1];
					ij--;
				}
			}
			if(vmat[0]==0){
				this->kvect[k] = 0;
				break;
			}
			TMRLOG(5,"\n\nV%d:\n", k);
			TMRLOGMATRIX(5, vmat, m,n, maxlen);
			htot=0;
			for(zindex=0; zindex<zsize; zindex++){
				// MINDICE(ligne, colone,len) ((ligne)+(colone)*(len))
				ij=MINDICE(zmat[zindex], zmat[zindex+zmatsize], maxlen);
				hmat[ij]=vmat[ij];
				htot+=vmat[ij];
			}
			TMRLOG(5,"\n\nH%d (htot=%g):\n", k, htot);
			TMRLOGMATRIX(5, hmat, m,n, maxlen);
			kvect[k] = htot;
			k++;
		}
	}
	for(;k<ksize;k++){
		kvect[k] = 0;
	}
	return;
	
	
}

