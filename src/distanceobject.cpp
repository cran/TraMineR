/* #include "distanceobject.h"

DistanceObject::DistanceObject(SEXP magicIndexS, SEXP magicSeqS){
  this->magicIndex=INTEGER(magicIndexS);
  this->magicSeq=INTEGER(magicSeqS);
  this->finalnseq=length(magicSeqS);
  PROTECT(ans = allocVector(REALSXP, (finalnseq*(finalnseq-1)/2)));
  result=REAL(ans);
}
DistanceObject::~DistanceObject(){
  UNPROTECT(1);
}

void DistanceObject::setDistance(const int &is,const int &js, const double& cmpres){
      int j_start=magicIndex[js];
		  int j_end=magicIndex[js+1];
		  int i_start=magicIndex[is];
		  int i_end=magicIndex[is+1];
		  int i_index, j_index, i, j, base_index;
		  for(i=i_start;i<i_end;i++){
		    i_index=magicSeq[i];
		    //n*(i-1) - i*(i-1)/2 + j-i
		    for(j=j_start;j<j_end;j++){
		      j_index=magicSeq[j];
		      if(i_index!=j_index){
		        base_index=distIndex(i_index,j_index,finalnseq);
//		        REprintf("Unique (%d,%d) => (%d,%d)(%d) => %f \n",is,js,i_index,j_index,(base_index),cmpres);
            result[base_index]=cmpres;
          }		      
        }
      }
}
 */
 
