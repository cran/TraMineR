#include<R.h>
#include "eventseq.h"
#include "prefixtree.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "eventdictionary.h"
//#include <math.h>

/**
	tmrsequence build a sequence obect and return an external pointer to that object
*/


#define TMRMATRIXINDEXC(ligne, colone,len) (ligne)+(colone)*len

extern "C" {
  

    SEXP tmrseqedist(SEXP seqs, SEXP Sidcost, SEXP vparam, SEXP Snorm, SEXP Sinterval) {
        double V=REAL(vparam)[0];
        int interval=INTEGER(Sinterval)[0];
        int norm=INTEGER(Snorm)[0];
        double* idcost=REAL(Sidcost);
        Sequence *s =NULL, *s2=NULL;
        int ns=length(seqs);
        SEXP ans;
        SEXP seq;
        PROTECT(ans = allocMatrix(REALSXP, ns, ns));
        double *matrix=REAL(ans);
        int maxevent=0, event=0;
        SequenceEventNode * sen=NULL, *sen2=NULL;
        for (int j=0;j<ns;j++) {
            seq=VECTOR_ELT(seqs,j);
            ASSIGN_TMRSEQ_TYPE(s,seq);
            if(s->hasEvent()){
              sen=s->getEvent();
              event=0;
              while(sen!=NULL){
                event++;
                sen=sen->getNext();
              }
              if(event>maxevent)maxevent=event;
            }
        }
        maxevent++;
        double* Fmat=new double[(maxevent)*(maxevent)];
        int ei, ej;
        double sen1cost, sen2cost, kcost, maxcost;
        double age1, age2, tmpage1;
        double scost=0, dist=0;
        for (int i=0;i<ns;i++) {
            seq=VECTOR_ELT(seqs,i);
            ASSIGN_TMRSEQ_TYPE(s,seq);
            matrix[TMRMATRIXINDEXC(i,i,ns)]=0;
            //Rprintf("Processing ");
            //sub->print();
            for (int j=i+1;j<ns;j++) {
                seq=VECTOR_ELT(seqs,j);
                ASSIGN_TMRSEQ_TYPE(s2,seq);
                sen=s->getEvent();
                ei=1;
                Fmat[0]=0;
                maxcost=0;
                while(sen!=NULL){
                  Fmat[TMRMATRIXINDEXC(ei,0,maxevent)]=Fmat[TMRMATRIXINDEXC(ei-1,0,maxevent)]+idcost[sen->getType()-1];
                  maxcost+=idcost[sen->getType()-1];
                  sen=sen->getNext();
                  ei++;
                }
                sen=s2->getEvent();
                ej=1;
                while(sen!=NULL){
                  Fmat[TMRMATRIXINDEXC(0,ej,maxevent)]=Fmat[TMRMATRIXINDEXC(0,ej-1,maxevent)]+idcost[sen->getType()-1];
                  maxcost+=idcost[sen->getType()-1];
                  sen=sen->getNext();
                  ej++;
                }
                ei=1;
                ej=1;
                
                sen=s->getEvent();
                age1=0;
                age2=0;
                dist=0;
                while(sen!=NULL){
                  sen2=s2->getEvent();
                  ej=1;
                  age1+=sen->getGap();
                  tmpage1=age1;
                  age2=0;
                  while(sen2!=NULL){
                    age2+=sen2->getGap();
                    
                    sen1cost=Fmat[TMRMATRIXINDEXC(ei-1,ej,maxevent)]+idcost[sen->getType()-1];
                    sen2cost=Fmat[TMRMATRIXINDEXC(ei,ej-1,maxevent)]+idcost[sen2->getType()-1];
                    scost=idcost[sen->getType()-1]+idcost[sen2->getType()-1];
                    if(sen->getType()==sen2->getType()){
                      kcost=V*fabs(tmpage1-age2);
                      if(kcost>scost)kcost=scost;
                    }else kcost=scost;
                    kcost+=Fmat[TMRMATRIXINDEXC(ei-1,ej-1,maxevent)];
                    
                    if(kcost<sen1cost&&kcost<sen2cost){
                        if(interval)age2=tmpage1;                      
                    }else if(sen1cost<sen2cost){
                        if(interval)age2=tmpage1;                      
                      kcost=sen1cost;
                    }else{
                      if(interval)age2=tmpage1;
                      kcost=sen2cost;
                    }
                   
                    Fmat[TMRMATRIXINDEXC(ei,ej,maxevent)]=kcost;
                    sen2=sen2->getNext();
                    ej++;
                  }
                  sen=sen->getNext();
                  ei++;
                }
                dist=Fmat[TMRMATRIXINDEXC(ei-1,ej-1,maxevent)];
                if(norm){
                  if(maxcost>0){
                    dist/=maxcost;
                  }else if(dist>0){
                    dist=1;
                  }
                }
                matrix[TMRMATRIXINDEXC(i,j,ns)]=dist;
                matrix[TMRMATRIXINDEXC(j,i,ns)]=dist;
                                                                 
            }
        }
        UNPROTECT(1);
        return ans;
    }

}

