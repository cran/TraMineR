#include <stack>
#include "salttseq.h"
#include <cmath>
#include "TraMineR.h"
//#include <math.h>

#define TMRMATRIXINDEXC(ligne, colone,len) (ligne)+(colone)*(len)
//#define TMRMIN(a,b) ((a)<(b))?a:b



static R_INLINE double computeDelta(double* previousmatrix, double * salttmatrix, const int& alphasize) {
	int i, alphafois=alphasize*alphasize;
	double delta=0;
	double diff;
	for(i=0;i<alphafois;i++) {
		diff = (double)previousmatrix[i] - (double)salttmatrix[i];
		delta = delta+(diff*diff);
	}
	delta=sqrt(delta/alphafois);
	return delta;
}
static R_INLINE double plog2(const double &value) {
	double res;
	res = std::log(value)/std::log((double)2);
	return res;
}

extern "C" {

SEXP saltt(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP indelS,
		SEXP alphasizeS, SEXP costsS, SEXP normS, SEXP optimS, SEXP pidS, SEXP maxpassS, SEXP retS, SEXP logoddS) {

	SEXP ans, scostmat;
	int nseq = INTEGER(seqdim)[0];
	double pid = REAL(pidS)[0];
	PROTECT(ans = allocVector(REALSXP, nseq * nseq));
	double * distmatrix = REAL(ans);
	//nb �tats
	int alphasize = INTEGER(alphasizeS)[0];
	int optim = INTEGER(optimS)[0];
	int ret = INTEGER(retS)[0];
	int logoddmode = INTEGER(logoddS)[0];
	PROTECT(scostmat = allocVector(REALSXP, (alphasize * alphasize)));
	double * salttcost = REAL(scostmat);

	int i, is, js, liseq, ljseq, fmatsize;
	//longueur des s�quences m=i, n=j

	//normalisation?
	int norm = INTEGER(normS)[0];
	//Nombre de s�quence
	//nb colonnes des s�quences
	int maxlen = INTEGER(seqdim)[1];
	//Matrice des s�quences
	int* sequences = INTEGER(Ssequences);
	//Tailles des s�quences
	int* slen = INTEGER(lenS);
	//indel
	double indel = REAL(indelS)[0];
	double * fmat;
	double * tbmat;
	int maxpass= INTEGER(maxpassS)[0];
	//Matrice des co�ts de substitutions
	double* scost = REAL(costsS);

	double *previousmatrix = new double[alphasize * alphasize];
	double delta = 1;

	std::stack<Alignement>* stackAligne;
	stackAligne = new std::stack<Alignement>;

	SEXP Fmat, TBmat;

	fmatsize=maxlen+1;
    PROTECT(Fmat = allocVector(REALSXP, (fmatsize*fmatsize)));
    fmat=REAL(Fmat);

    PROTECT(TBmat = allocVector(REALSXP, (fmatsize*fmatsize)));
    tbmat=REAL(TBmat);



	for (is = 0; is < nseq; is++) {
		liseq = slen[is] + 1;
		for (js = is+1; js < nseq; js++) {
			ljseq = slen[js] + 1;
			//if(js!=is) {
				Alignement align = Alignement(is, js, 0, liseq, ljseq, maxlen);
				stackAligne->push(align);
			//}
		}
	}

	//REprintf("stack size = %d", stackAligne->size());

	Salttseq seq1 = Salttseq(norm, nseq, slen, maxlen, &indel, alphasize, scost,
			distmatrix, stackAligne, sequences, salttcost, fmat, tbmat, fmatsize, pid, logoddmode);
	int pass = 0;

	if (optim == 1) {
		while (delta > 0.0001) {
			if(pass==maxpass) {
				REprintf("\nWe have reached the maximum number of iteration (%d). We stop here.\n", pass);
				break;
			}



			for (i = 0; i < alphasize * alphasize; i++) {
				previousmatrix[i] = salttcost[i];
			}

			if(seq1.computeDistances(1)==-1) {
				REprintf("\n************************\nWarning\n************************\n");
				REprintf("No alignement has a PID>%f, must stop here\n", pid);
				delta = 0;
			}
			/*REprintf("New cost matrix, after pass #%d \n", pass);
			seq1.printsalttmatrix();
			REprintf("indel = %f\n", indel);
			*/
			delta = computeDelta(previousmatrix, salttcost, alphasize);
			REprintf("\ndelta = %f\n", delta);
			pass++;
		}

		REprintf("Final optimal matching on all sequences\n");
		for (is = 0; is < nseq; is++) {
				liseq = slen[is] + 1;
				for (js = is + 1; js < nseq; js++) {
					ljseq = slen[js] + 1;
					Alignement align = Alignement(is, js, 0, liseq, ljseq, maxlen);
					stackAligne->push(align);
				}
			}

		seq1.computeDistances(0);
	} else {
		seq1.computeDistances(0);
	}
	UNPROTECT(4);
	delete previousmatrix;
	if(ret==1) {
		return ans;
	}
	else {
		return scostmat;
	}
}

SEXP henikoff(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP alphasizeS, SEXP logoddS) {
	SEXP scostmat;
	REprintf("AAA");
	int nseq = INTEGER(seqdim)[0];
	int alphasize = INTEGER(alphasizeS)[0];
	//int logoddmode = INTEGER(logoddS)[0];
	PROTECT(scostmat = allocVector(REALSXP, (alphasize * alphasize)));
	double * scosts = REAL(scostmat);
	//Nombre de s�quence
	//nb colonnes des s�quences
	int maxlen = INTEGER(seqdim)[1];
	//Matrice des s�quences
	int* sequences = INTEGER(Ssequences);
	//Tailles des s�quences
	//int* slen = INTEGER(lenS);
	//indel
	int i,j,k, statea, stateb;
	double ** fmat = new double*[alphasize];
	double ** qijtable = new double*[alphasize];
	double ** eijtable = new double*[alphasize];

	for(i=0;i<alphasize;i++) {
		fmat[i] = new double[alphasize];

		eijtable[i] = new double[alphasize];
	}
	for(i=0;i<alphasize;i++) {
		qijtable[i] = new double[alphasize];
	}
	double * pri = new double[alphasize];

	for(i=0;i<alphasize;i++) {
		pri[i]=0;
		for(j=0;j<alphasize;j++) {
			fmat[i][j] = fmat[j][i] = qijtable[i][j] = qijtable[j][i] = scosts[TMRMATRIXINDEXC(i,j,alphasize)] = scosts[TMRMATRIXINDEXC(j,i,alphasize)] = 0;
		}
	}


	//Matrice des co�ts de substitutions


	double maxn = ((nseq*(nseq-1))/2)*maxlen;

	for(i=0;i<maxlen;i++) {
		for(j=0;j<nseq-1;j++) {
			for(k=j+1;k<nseq;k++) {

				statea = sequences[TMRMATRIXINDEXC(j,i,nseq)];
				stateb = sequences[TMRMATRIXINDEXC(k,i,nseq)];
				if(statea!=stateb) {
					fmat[statea][stateb]+=1;
					fmat[stateb][statea]+=1;
				}
				else {
					fmat[statea][stateb]+=1;
				}
			}
		}
	}


	for(i=0;i<alphasize;i++) {
		//REprintf("Init ok\n");
		for(j=0;j<alphasize;j++) {
			qijtable[i][j] = qijtable[j][i] = (fmat[i][j])/maxn;
		}
	}
	int a;
	double sumqij;

	for(a=0;a<alphasize;a++) {
		sumqij=0;
		for(i=0;i<alphasize-1;i++) {
			for(j=i+1;j<alphasize;j++) {
				sumqij+=qijtable[i][j];
			}
		}
		pri[a]=qijtable[a][a]+(sumqij/2);
		REprintf("pri[%d] = %f\n", a, pri[a]);
	}

	for(i=0;i<alphasize;i++) {
		for(j=i;j<alphasize;j++) {
			if(i==j) {
				eijtable[i][j] = eijtable[j][i] = pri[i]*pri[j];
			}
			else {
				eijtable[i][j] = eijtable[j][i] = pri[i]*pri[j]*2;
			}
		}
	}

	REprintf("Last loop\n");
	for(i=0;i<alphasize;i++) {
		for(j=i+1;j<alphasize;j++) {
			scosts[TMRMATRIXINDEXC(i,j,alphasize)] = scosts[TMRMATRIXINDEXC(j,i,alphasize)] = plog2(qijtable[i][j]/eijtable[i][j]);
			//- ((log2(qijtable[i][i]/eijtable[i][i])+log2(qijtable[j][j]/eijtable[j][j]))/2);
			REprintf("scost[%d][%d] = %f\n",i,j,scosts[TMRMATRIXINDEXC(i,j,alphasize)]);
		}
	}




	for(i=0;i<alphasize;i++) {
		delete[] fmat[i];
		delete[] qijtable[i];
		delete[] eijtable[i];
 	}
	delete[] fmat;
	delete[] pri;
	delete[] qijtable;
	delete[] eijtable;

	REprintf("End");
	UNPROTECT(1);
	return(scostmat);

	}
}
