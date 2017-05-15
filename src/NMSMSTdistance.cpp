#include "NMSMSTdistance.h"
NMSMSTdistance::NMSMSTdistance(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS):SUBSEQdistance(normS, Ssequences, seqdim, lenS), seqdur(NULL){
	rowsize=maxlen+1;
	int matrixsize=rowsize*rowsize;
	e1=new double[matrixsize];
	e=new double[matrixsize];
	t1=new double[matrixsize];
	t=new double[matrixsize];
}
NMSMSTdistance::NMSMSTdistance(NMSMSTdistance *dc):SUBSEQdistance(dc), seqdur(dc->seqdur){
	rowsize=maxlen+1;
	int matrixsize=rowsize*rowsize;
	e1=new double[matrixsize];
	e=new double[matrixsize];
	t1=new double[matrixsize];
	t=new double[matrixsize];
}
NMSMSTdistance::~NMSMSTdistance(){
	delete [] e1;
	delete [] e;
	delete [] t1;
	delete [] t;

}
void NMSMSTdistance::setParameters(SEXP params){
	this->seqdur=REAL(getListElement(params, "seqdur"));
	SUBSEQdistance::setParameters(params);
}
void NMSMSTdistance::computeattr(const int&is, const int& js){
	int mrows=slen[is]+1;
	int ncols=slen[js]+1;
	int m=mrows-1,n=ncols-1, ij=0;
	double st=0.0,sum=0.0,temp=0.0,temps=0.0;
	double sumt,tempt,temptt;
	int seqi=-1;
	double ti;
	int k=0;
	for (int i=0;i<m;i++) {
		seqi=sequences[MINDICE(is,i,nseq)];
		ti=seqdur[MINDICE(is,i,nseq)];
		for (int j=0;j<n;j++) {
			ij = MINDICE(i,j,rowsize);
			if(seqi==sequences[MINDICE(js,j,nseq)]) {
				e1[ij]=1;
				e[ij]=1;
				t1[ij]=fmin2(ti, seqdur[MINDICE(js,j,nseq)]);
				t[ij]=t1[ij];
				st+=t[ij];
				if(st==DBL_MAX){
					error(" [!] Number of subsequences is getting too big");
				} 
			}
			else{
				e1[ij]=0.0;
				e[ij]=0.0;
				t1[ij]=0.0;
				t[ij]=0.0;
			}
		}
	}
	for(int i=0; i<m;i++){
		ij = MINDICE(i,ncols-1,rowsize);
		e1[ij]=0.0;
		e[ij]=0.0;
		t1[ij]=0.0;
		t[ij]=0.0;
	}
	for(int i=0; i<ncols;i++) {
		ij = MINDICE(mrows-1, i,rowsize);
		e1[ij]=0.0;
		e[ij]=0.0;
		t1[ij]=0.0;
		t[ij]=0.0;
	}
	TMRLOG(5,"\n\ne1\n");
	TMRLOGMATRIX(5, e1, mrows,ncols, rowsize);
	TMRLOG(5,"\n\ne\n");
	TMRLOGMATRIX(5, e, mrows,ncols, rowsize);
	TMRLOG(5,"\n\nt\n");
	TMRLOGMATRIX(5, t, mrows,ncols, rowsize);
	TMRLOG(5,"\n\nt1\n");
	TMRLOGMATRIX(5, t1, mrows,ncols, rowsize);
	this->kvect[k]=st;
	if(st==0.0){
		return;
	}
	while(mrows !=0 && ncols != 0) {
		k++;
		for (int irow=0;irow<mrows;irow++) {
			
			sum=0.0;
			sumt=0.0;
			for (int jcol=ncols-1;jcol>=0;jcol--) {
				ij = MINDICE(irow,jcol,rowsize);
				temp=sum;
				tempt=sumt;
				sum+=e[ij];
				sumt+=t[ij];
				e[ij]=temp;
				t[ij]=tempt;
			}
		}
		temps=0.0;
		temptt=0.0;
		for (int jcol=0;jcol<ncols;jcol++) {
			sum=0;
			sumt=0.0;
			for (int irow=mrows-1;irow>=0;irow--) {
				ij = MINDICE(irow, jcol, rowsize);
				temp=sum;
				tempt=sumt;
				sum+=e[ij];
				sumt+=t[ij];
				e[ij]=temp*e1[ij];
				t[ij]=e1[ij]*(e[ij]*t1[ij]+tempt);
				temps+=e[ij];
				temptt+=t[ij];
			}
		}
		if(temps==0.0){
			return;
		}
		this->kvect[k]=temptt;
		if(temptt==DBL_MAX){
			error(" [!] Number of subsequences is getting too big");
		} 
		mrows--;
		ncols--;
	}
}



// static R_INLINE void matrixfill(double *mat, const double & val, const int &row, const int & col, const int &matsize){
	// for(int i=0; i<row; i++){
		// for(int j=0; j<col; j++){
			// mat[MINDICE(i,j,matsize)] = val;
		// }
	// }
// }

// static R_INLINE void computevmatrix(double *vmat, const int &i, const int& j, const int &maxlen, double *hmat, const double &total){
	// if(i<0 || j<0) {
		// return;
	// }
	// ij= MINDICE(i,j,maxlen);
	// vmat[ij] = total;
	// computevmatrix(vmat, i-1, j, maxlen, hmat, total);
	// computevmatrix(vmat, i, j-1, maxlen, hmat, total);
	
	
	// vmat[ij] = computevmatrix(vmat, i+1, j+1, maxlen, hmat);
	// return vmat[ij]+ hmat[ij];
	// for(i=0; i<m; i++){
			// for(j=0; j<n; j++){
				
				// vmat[ij]=0;
				// for(a=i+1; a<m; a++){
					// for(b=j+1; b<n; b++){
						// vmat[ij] += hmat[MINDICE(a,b,maxlen)];
					// }
				// }
			// }
		// }
// }
	/***********************************************
	* Calculates the minimal shared time (MST) spent
	* in all matching subsequences from strings 
	* seq[iseq] and seq[jseq]
	***********************************************/
/* static R_INLINE double NMSMSTdistance(int* slen,const int &is,const int &js, const int&nseq, int* sequences, const int &maxlen, 
								double * hmat, double * vmat, double * zmat, double * selfmatvect, double *ijvect ) {

	int mrows=slen[is]+1;
	int ncols=slen[js]+1;
	int m=mrows-1,n=ncols-1;
	int[,] e1=new int[mrows,ncols];
	int[,] e=new int[mrows,ncols];
	double[,] t1=new double[mrows,ncols];
	double[,] t=new double[mrows,ncols]; 
	double st=0.0,sum=0.0,temp=0.0,temps=0.0;
	double sumt,tempt,temptt;
	int seqi=-1;
	double ti;
	for (int i=0;i<m;i++) {
		seqi=DataSpec.seq[iseq][i];
		ti=DataSpec.quant[iseq][i];
		for (int j=0;j<n;j++) {
			if(seqi==DataSpec.seq[jseq][j]) {
				e1[i,j]=1;
				e[i,j]=1;
				t1[i,j]=Math.Min(ti,DataSpec.quant[jseq][j]);
				t[i,j]=t1[i,j];
				st+=t[i,j];
			}
		}
	}
	if(st==0.0){
		return st;
	}
	while(mrows !=0 && ncols != 0) {
		for (int irow=0;irow<mrows;irow++) {
			sum=0.0;
			sumt=0.0;
			for (int jcol=ncols-1;jcol>=0;jcol--) {
				temp=sum;
				tempt=sumt;
				sum+=e[irow,jcol];
				sumt+=t[irow,jcol];
				e[irow,jcol]=(int)temp;
				t[irow,jcol]=tempt;
			}
		}
		temps=0.0;
		temptt=0.0;
		for (int jcol=0;jcol<ncols;jcol++) {
			sum=0;
			sumt=0.0;
			for (int irow=mrows-1;irow>=0;irow--) {
				temp=sum;
				tempt=sumt;
				sum+=e[irow,jcol];
				sumt+=t[irow,jcol];
				e[irow,jcol]=(int)temp*e1[irow,jcol];
				t[irow,jcol]=e1[irow,jcol]*(e[irow,jcol]*t1[irow,jcol]+tempt);
				temps+=e[irow,jcol];
				temptt+=t[irow,jcol];
			}
		}
		if(temps==0.0){
			return st;
		}
		st+=temptt;
		mrows--;
		ncols--;
	} 
	return st;
}
 */
 
