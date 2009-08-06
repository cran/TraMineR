#include "salttseq.h"
#include "alignement.h"
#include "TraMineR.h"
#include <stack>
#include <set>
#include <algorithm>
#include <utility>

using namespace std;

Salttseq::Salttseq(const int& inorm, const int& inseq, int * islen, const int&imaxlen, double*  iindel, int  &ialphasize, double * iscost, double * idistmatrix, std::stack<Alignement>* istackAlign, int * isequences, double * isalttcost, double * ifmat, double * itbmat, const int& ifmatsize, const double& ipid, const int& ilogoddmode) {
	this->norm=inorm;
	this->nseq=inseq;
	this->maxlen=imaxlen;
	this->indel=iindel;
	this->alphasize=ialphasize;
	this->scost=iscost;
	this->stackAlign=istackAlign;
	this->frequences = new double[alphasize];
	this->alphafois = alphasize*alphasize;
	this->freqconj = new double[(int)alphafois];
	this->sequences=isequences;
	this->slen=islen;
	this->salttcost=isalttcost;
	this->fmat = ifmat;
	this->tbmat = itbmat;
	//	this->stackAlign2 = new std::stack<Alignement>;
	this->fmatsize=ifmatsize;
	this->pid=ipid;
	this->totalalignements=0;
	this->log_odd_mode = ilogoddmode;

	int i,j;
	maxscost=0;
	this->distmatrix=idistmatrix;


    for (i=0;i<alphasize;i++) {
      for(j=i; j<alphasize;j++){
        if (scost[TMRMATRIXINDEXC(i,j,alphasize)]>maxscost) {
          maxscost=scost[TMRMATRIXINDEXC(i,j,alphasize)];
        }
      }
    }
    maxscost=fmin2(maxscost,2*(*indel));
    //Initialisation, peut �tre fait qu'une fois
    for (i=0;i<fmatsize;i++) {
      fmat[TMRMATRIXINDEXC(i,0,fmatsize)]=fmat[TMRMATRIXINDEXC(0,i,fmatsize)]=i*(*indel);
      tbmat[TMRMATRIXINDEXC(i,0,fmatsize)]=2;
      tbmat[TMRMATRIXINDEXC(0,i,fmatsize)]=1;

    }
    tbmat[TMRMATRIXINDEXC(0,0,fmatsize)]=0;

    // initialisation de la matrice des couts
    // the initial costs are copied into salttcost. From now on we will modify only salttcost.
    for (i=0;i<alphafois;i++) {
    	salttcost[i] = scost[i];
    }


 }
Salttseq::~Salttseq() {
	delete this->frequences;
	delete this->freqconj;
	delete this->stackAlign;
	//delete this->stackAlign2;
}
double* Salttseq::getdist() {
	return(this->distmatrix);
}
// withsaltt = 1 if we want to compute optimal costs
// withsaltt = 0 if we just want an optimal matching distance matrix
int Salttseq::computeDistances(const int& withsaltt) {
	int i;
	double currentpid;

	// initialization of the frequencies tables
	for (i=0;i < alphasize;i++) {
		frequences[i] = 0.001;
	}

	// to avoid log = -inf, we set an initial value of 0.0001 instead of 0
	for (i=0;i< alphafois; i++) {
		freqconj[i] = 0.001;
	}


	int is, js, m, n;

	// initialization of the stack that will contain the new set of alignement (> pid)
	std::stack<Alignement>* stackAlign2;
	stackAlign2 = new std::stack<Alignement>;

	// initialization of the levenshtein and traceback matrices
    for (i=0;i<fmatsize;i++) {
      fmat[TMRMATRIXINDEXC(i,0,fmatsize)]=fmat[TMRMATRIXINDEXC(0,i,fmatsize)]=i*(*indel);
      tbmat[TMRMATRIXINDEXC(i,0,fmatsize)]=2;
      tbmat[TMRMATRIXINDEXC(0,i,fmatsize)]=1;

    }
    tbmat[TMRMATRIXINDEXC(0,0,fmatsize)]=0;


	//double cmpres=0;
	// will contain the total number of states in the set of alignements (in order to compute f_a)
	this->totalstates=0;
	// will contain the total number of state alignements (in order to compute f_a,b)
	this->totalalignements=0;
	//REprintf("stack size = %d, totalstates = %f, totalalignements = %f\n", stackAlign->size(), totalstates, totalalignements);

	// dépilage des alignements
	while(!stackAlign->empty()) {
		Alignement align = stackAlign->top();
		m = align.getliseq();
		n = align.getljseq();
		is = align.getiseq();
		js = align.getjseq();
		//Distance bewteen a sequence and itself is 0
		distmatrix[TMRMATRIXINDEXC(is,is,nseq)]=0;
		distmatrix[TMRMATRIXINDEXC(js,js,nseq)]=0;

		//Distance between is and js is computed by levenshtein()
		distmatrix[TMRMATRIXINDEXC(js,is,nseq)]=distmatrix[TMRMATRIXINDEXC(is,js,nseq)]=levenshtein(is,js,m,n);

		stackAlign->pop();

		// percentage identity of the alignement
		currentpid = getPID(is, js, m, n, align);
		//REprintf("align pid = %f\n", currentpid);

		// we select alignement with a higher pid than the threshold (pid)
		if(currentpid > pid) {
					align.setpid(currentpid);
					// state frequencies are computed and weighted by currentpid
					computeConjfreq(is,js,m,n,currentpid);
					// we fill the new stack
					stackAlign2->push(align);
			}
	}


	moveStack(stackAlign2, stackAlign);
	delete stackAlign2;
	//delete stackAlign;
	//this->stackAlign = stackAlign2;



	// if the new set of alignements is not empty, we compute the new cost matrix
	if(stackAlign->size() > 0 && withsaltt==1) {
		//frequencies are now computed in the computeConjfreq() method
		computeDayhoffcosts();
		return(1);
	}
	// if it is empty, or if we don't want the cost matrix, we stop here
	return(-1);
}

void Salttseq::computeDayhoffcosts() {
	double dayhoff;
	// table that will contain the dayhoff cost for the same state (Dayhoff(a,a))
	double *dayhoffaa=new double[alphasize*alphasize];
	double newindel=0;
	int i, j;
	double r1, r2;
	//printfreqtables();
	if(log_odd_mode==0) {
		for(i=0;i<alphasize;i++) {
			for(j=0; j<alphasize;j++) {

			    r1=(freqconj[TMRMATRIXINDEXC(i,j,alphasize)])/(freqconj[TMRMATRIXINDEXC(i,i,alphasize)]);
			    r2=(freqconj[TMRMATRIXINDEXC(j,i,alphasize)])/(freqconj[TMRMATRIXINDEXC(j,j,alphasize)]);
			    dayhoff=( ( log(r1) ) + ( log(r2) ) )/2;
				dayhoffaa[TMRMATRIXINDEXC(i, j, alphasize)] = (-1)*dayhoff;

			}
		}
	}


	else if(log_odd_mode==1) {
		for(int i=0;i<alphasize;i++) {
			for(int j=0; j<alphasize;j++) {
				r1 = (freqconj[TMRMATRIXINDEXC(i,j,alphasize)]+freqconj[TMRMATRIXINDEXC(j,i,alphasize)]+1)/(totalstates+1);
				r2 = ((frequences[i])/(totalstates))*((frequences[j])/(totalstates))*2;
				dayhoff=log(r1/r2);
				dayhoffaa[TMRMATRIXINDEXC(i, j, alphasize)] = (-1)*dayhoff;
			}
		}
	}

	else if(log_odd_mode==2) {
		double p_a, p_b, p_ab, p_aSb, p_bSa;
		for(int i=0;i<alphasize;i++) {
			for(int j=0; j<alphasize;j++) {

				p_a = frequences[i]/totalstates;
				p_b = frequences[j]/totalstates;
				p_ab = (freqconj[TMRMATRIXINDEXC(i,j,alphasize)] + freqconj[TMRMATRIXINDEXC(j,i,alphasize)])/totalstates;
				p_aSb = p_ab/p_b;
				p_bSa = p_ab/p_a;
				r1 = (p_aSb + p_bSa)/2;
				dayhoff = (-1)*log(r1/(p_a*p_b));
				dayhoffaa[TMRMATRIXINDEXC(i,j,alphasize)] = dayhoff;

			}
		}
	}

    // normalize dayhoff
    for(int i=0;i<alphasize;i++) {
    	for(int j=0; j<alphasize;j++) {
    		dayhoff = dayhoffaa[TMRMATRIXINDEXC(i, j, alphasize)] - ((dayhoffaa[TMRMATRIXINDEXC(i, i, alphasize)] + dayhoffaa[TMRMATRIXINDEXC(j, j, alphasize)])/2);
    		salttcost[TMRMATRIXINDEXC(i, j, alphasize)] = dayhoff;
    		//salttcost[TMRMATRIXINDEXC(i, j, alphasize)]=dayhoff = dayhoffaa[TMRMATRIXINDEXC(i, j, alphasize)];
    		if(j>i) newindel += (dayhoff);
		}
    }

	//*indel = newindel/((((double)alphasize*alphasize)-(double)alphasize)*0.5);
	*indel = newindel/(double)((alphasize*(alphasize-1))*0.5);
    //REprintf("new indel = %f\n", *indel);
	delete dayhoffaa;
}


// needleman-wunsch algorithm. return the levenshtein distance between is and js
double Salttseq::levenshtein(const int &is, const int &js, int &m, int &n) {
	int i, j, i_state, j_state, tbdirection;
	double cost, in, del, sub, minimum, maxpossiblecost=0, cmpres;

	for(i=1;i<m;i++) {
		i_state=sequences[TMRMATRIXINDEXC(is,i-1,nseq)];

              for(j=1;j<n;j++) {

                         j_state=sequences[TMRMATRIXINDEXC(js,j-1,nseq)];
                         if (i_state == j_state) {
                             cost = 0;
                         } else {
                             cost = salttcost[TMRMATRIXINDEXC(i_state,j_state,alphasize)];
                         }

                         in = fmat[TMRMATRIXINDEXC(i,j-1,fmatsize)]+ *indel;
                         del = fmat[TMRMATRIXINDEXC(i-1,j,fmatsize)]+ *indel;
                         sub = fmat[TMRMATRIXINDEXC(i-1,j-1,fmatsize)]+ cost;

                         minimum=in;
                         // in
                         //tbmat[TMRMATRIXINDEXC(i,j,fmatsize)]
                         tbdirection=1;
                         if(del < minimum) {
                         	minimum = del;
                         	// del
                         	tbdirection=2;
                         }
                         if(sub < minimum) {
                         	minimum = sub;
                         	// sub
                         	tbdirection=3;
                         }
                         fmat[TMRMATRIXINDEXC(i,j,fmatsize)]=minimum;
						 tbmat[TMRMATRIXINDEXC(i,j,fmatsize)]=tbdirection;
                 	}
                }
	m--; n--;
	maxpossiblecost=abs(n-m)*(*indel)+maxscost*fmin2((double)m,(double)n);
	cmpres=normalizeDistance(fmat[TMRMATRIXINDEXC(m,n,fmatsize)], maxpossiblecost, m, n, norm);
	//printtraceback(m,n);
	//cmpres=fmat[TMRMATRIXINDEXC(m-1,n-1,fmatsize)];
	return(cmpres);
}


// compute the percentage identity between two sequences using the traceback matrix
double Salttseq::getPID(const int &is, const int &js, const int&m, const int &n, Alignement &align) {
	int i,j,maxmn, i_state, j_state, tb_state, matches=0;
	double cpid=0;
	i = m-1;
	j = n-1;
	while(i>0 || j>0) {
		i_state = sequences[TMRMATRIXINDEXC(is,i-1,nseq)];
		j_state = sequences[TMRMATRIXINDEXC(js,j-1,nseq)];
		tb_state = tbmat[TMRMATRIXINDEXC(i,j,fmatsize)];
		if(tb_state==3) {
			if(i_state==j_state) {
				matches++;
			}
			 i--;
			 j--;
		}
		else if(tb_state==2) {
 			i--;
		}
		else if(tb_state==1) {
			j--;
		}
	}
	maxmn = max(m,n);

	maxmn--;
	cpid=(double)matches/(double)maxmn;
	//REprintf("\n matches = %d, max = %d, PID = %f\n", matches, maxmn, pid);
	return(cpid);
}

void Salttseq::computeConjfreq(const int &is, const int &js, const int&m, const int &n, const double &currentpid) {
	int i_state, j_state, tb_state;
	int i = m-1;
	int j = n-1;
	//stack<pair<int,int> >* alignpairs = new stack<pair<int,int> >();
	while(i>0 || j>0) {
		pair<int, int> pr1;
		i_state = sequences[TMRMATRIXINDEXC(is,i-1,nseq)];
		j_state = sequences[TMRMATRIXINDEXC(js,j-1,nseq)];
		tb_state = tbmat[TMRMATRIXINDEXC(i,j,fmatsize)];

		if(tb_state==3 || tb_state==0) {
			//pr1.first = i_state;
			//pr1.second = j_state;
			freqconj[TMRMATRIXINDEXC(i_state, j_state, alphasize)] +=  currentpid;
			frequences[i_state] += currentpid;
			frequences[j_state] += currentpid;
			totalstates+=2*currentpid;
			i--;
			j--;
		}
		else if(tb_state==2) {
			//pr1.first = i_state;
			//pr1.second = -1;
			frequences[i_state] += currentpid;
			totalstates+=currentpid;
 			i--;
		}
		else if(tb_state==1) {
			//pr1.first = -1;
			//pr1.second = j_state;
			frequences[j_state] += currentpid;
			totalstates+=currentpid;
			j--;
		}
		//alignpairs->push(pr1);
	}
	//align.setAlign(alignpairs, frequences, freqconj);
	//align.printAlign();
	//REprintf("alignpairs.size = %d\n", alignpairs->size());
}

void Salttseq::moveStack(std::stack<Alignement>* oldstack, std::stack<Alignement>* newstack) {
	if(!newstack->empty()) {
		REprintf("warning: new stack is not empty, stopping operation\n");
	}
	else {
		while(!oldstack->empty()) {
			newstack->push(oldstack->top());
			oldstack->pop();
		}
		//REprintf("new stack has now %d elements\n", newstack->size());
	}
}


double Salttseq::normalizeDistance(const double& rawdist, const double& maxdist, const int& l1, const int& l2, const int&inorm) {
  this->norm=inorm;
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


// print the frequency and conjoint frequencies tables
void Salttseq::printfreqtables() {
	int i, j;
		REprintf("frequences\n");
	for (i = 0; i < alphasize; i++) {
		REprintf("[%d]:%f\t", i, frequences[i]);
	}

	REprintf("\n");
	REprintf("\nfrquences conjointes\n");
	for (i=0;i < alphasize; i++) {
		REprintf("\n");
		for (j=0;j< alphasize; j++) {
			REprintf("[%d][%d]=%f | ", i, j, freqconj[TMRMATRIXINDEXC(i,j,alphasize)]);
		}
	}
	REprintf("\n");
}

// print the cost matrix
void Salttseq::printsalttmatrix() {
	int i, j;
	for (i=0;i<alphasize;i++) {
		for (j=0;j<alphasize;j++) {
			REprintf("[%d][%d]=%f \t", i, j, salttcost[TMRMATRIXINDEXC(i,j,alphasize)]);
		}
		REprintf("\n");
	}
}

void Salttseq::printtraceback(int m, int n) {
	int i,j;
	for(i=0;i<m;i++) {
		Rprintf("{");
		for(j=0;j<n;j++) {
			Rprintf("[%f]", tbmat[TMRMATRIXINDEXC(i,j,fmatsize)]);
		}
		Rprintf("}\n");
	}
}
