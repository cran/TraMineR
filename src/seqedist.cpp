#include<R.h>
#include "eventseq.h"
#include "prefixtree.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "eventdictionary.h"
#include "TraMineR.h"
//#include <math.h>

	//put indel costs into costs. costs[1] => indel1 costs[2] => indel2  costs[3] => subst cost
	void getIndelSubstitutionCost(SequenceEventNode * sen1, SequenceEventNode *sen2, double* idcost, double& indel1, double& indel2, double& subst){
		SequenceEventNode * s1=sen1;
		SequenceEventNode * s2=sen2;
		indel1=indel2=subst=0;
		TMRLOG(5, "Treating S1");
		do{
			indel1+= idcost[s1->getType()-1];
			TMRLOG(5, "Treating %d", s1->getType());
			s1 =s1->getNextWithoutGap();
			
			
		}while(s1!=NULL);
		TMRLOG(5, "Treating S2");
		do{
			indel2+= idcost[s2->getType()-1];
			TMRLOG(5, "Treating S2 %d", s2->getType());
			s2 =s2->getNextWithoutGap();
		}while(s2!=NULL);
		
		s1=sen1;
		s2=sen2;
		do{
			if(s1->getType()==s2->getType()){
				s1= s1->getNextWithoutGap();
				s2= s2->getNextWithoutGap();
			}
			else if(s1->getType()<s2->getType()){
				subst+= idcost[s1->getType()-1];
				s1= s1->getNextWithoutGap();
			}else{
				subst+= idcost[s2->getType()-1];
				s2= s2->getNextWithoutGap();
			}
			if(s1==NULL){
				while(s2!=NULL){
					subst+= idcost[s2->getType()-1];
					s2= s2->getNextWithoutGap();
				}
				return;
			}
			if(s2==NULL){
				do{
					subst+= idcost[s1->getType()-1];
					s1= s1->getNextWithoutGap();
				}while(s1!=NULL);
				return;
			}
		}while(true);
	
	}

/**
	tmrsequence build a sequence obect and return an external pointer to that object
*/

extern "C" {
	

    SEXP tmrseqedist(SEXP seqs, SEXP Sidcost, SEXP vparam, SEXP Snorm, SEXP Sinterval) {
        double V=REAL(vparam)[0];
        int interval=INTEGER(Sinterval)[0];
        int norm=INTEGER(Snorm)[0];
        double* idcost=REAL(Sidcost);
        Sequence *s =NULL, *s2=NULL;
        int ns=Rf_length(seqs);
        SEXP ans;
        SEXP seq;
        PROTECT(ans = Rf_allocMatrix(REALSXP, ns, ns));
        double *matrix=REAL(ans);
        int maxevent=0, event=0;
        SequenceEventNode * sen=NULL, *sen2=NULL;
		
		//Calcul de la taille de la matrice, nb max event est une borne supérieur
        for (int j=0;j<ns;j++) {
            seq=VECTOR_ELT(seqs,j);
            ASSIGN_TMRSEQ_TYPE(s,seq);
            if(s->hasEvent()){
              sen=s->getEvent();
              event=0;
              while(sen!=NULL){
				if(event==0 || sen->getGap()>0){
					event++;
				}
                sen=sen->getNext();
              }
              if(event>maxevent)maxevent=event;
            }
        }
        maxevent++;
		TMRLOG(5, "Max Event %d", maxevent);
        double* Fmat=new double[(maxevent)*(maxevent)];
        int ei, ej;
        double sen1cost, sen2cost, kcost, maxcost;
        double age1, age2, tmpage1, tmpage2;
        double scost=0, scosttot=0, dist=0;
        for (int i=0;i<ns;i++) {
            seq=VECTOR_ELT(seqs,i);
            ASSIGN_TMRSEQ_TYPE(s,seq);
			//distance on diagonal => 0
            matrix[MINDICE(i,i,ns)]=0;
            //Rprintf("Processing ");
            //sub->print();
            for (int j=i+1;j<ns;j++) {
                seq=VECTOR_ELT(seqs,j);
                ASSIGN_TMRSEQ_TYPE(s2,seq);
                sen=s->getEvent();
				//Initializing the Fmatrix
				TMRLOG(5, "Initializing the Fmatrix");
                ei=1;
                Fmat[0]=0;
                maxcost=0;
				TMRLOG(5, "S1");
                while(sen!=NULL){
					sen2=sen;
					sen1cost=0;
					do{
						sen1cost+=idcost[sen2->getType()-1];
						TMRLOG(5, "S1 %d", sen2->getType());
						sen2=sen2->getNextWithoutGap();
					}while(sen2!=NULL);
					TMRLOG(5, "GAP %f", sen1cost);
					maxcost+=sen1cost;
					Fmat[MINDICE(ei,0,maxevent)]=Fmat[MINDICE(ei-1,0,maxevent)]+sen1cost;
					sen=sen->getNextWithGap();
					ei++;
                }
				
                sen=s2->getEvent();
                ej=1;
				TMRLOG(5, "S2");
                while(sen!=NULL){
					sen2=sen;
					sen1cost=0;
					do{
						sen1cost+=idcost[sen2->getType()-1];
						TMRLOG(5, "%d", sen2->getType());
						sen2=sen2->getNextWithoutGap();
					}while(sen2!=NULL);
					TMRLOG(5, "GAP %f", sen1cost);
					Fmat[MINDICE(0,ej,maxevent)]=Fmat[MINDICE(0,ej-1,maxevent)]+sen1cost;
					maxcost+=sen1cost;
					sen=sen->getNextWithGap();
					ej++;
                }
				
                ei=1;
                ej=1;
                TMRLOG(5, "Maximum cost %f", maxcost);
				//return Rf_ScalarInteger(1);
                sen=s->getEvent();
                age1=0;
                age2=0;
                dist=0;
                while(sen!=NULL){
					sen2=s2->getEvent();
					ej=1;
					
					if(interval==2) {
						age1=sen->getGap();
					}else if(interval==3){
						if(sen->getNextWithGap()!=NULL){
							tmpage1=sen->getNextWithGap()->getGap();
							age1+=tmpage1;
						}else{
							tmpage1=s->getObsTime()-age1;
						}
					} else {
						age1+=sen->getGap();
					}
					tmpage1=age1;
					age2=0;
					while(sen2!=NULL){
						if(interval==2) {
							age2=sen2->getGap();
						}else if(interval==3){
							if(sen2->getNextWithGap()!=NULL){
								tmpage2=sen2->getNextWithGap()->getGap();
								age2+=tmpage2;
							}else{
								tmpage2=s->getObsTime()-age2;
							}
						} else {
							age2+=sen2->getGap();
						}
						TMRLOG(5, "S1 %d, S2 %d", sen->getType(), sen2->getType());
						getIndelSubstitutionCost(sen, sen2, idcost, sen1cost, sen2cost, scost);
						scosttot=sen1cost+sen2cost;
						TMRLOG(5, "Indel1 %f, Indel2 %f, Scost %f, Scosttot %f ", sen1cost, sen2cost, scost, scosttot);
						sen1cost=Fmat[MINDICE(ei-1,ej,maxevent)]+sen1cost;
						sen2cost=Fmat[MINDICE(ei,ej-1,maxevent)]+sen2cost;
						
						//computing real kcost
						kcost=V*fabs(tmpage1-age2)+scost;
						if(kcost>scosttot) kcost=scosttot;
						kcost+=Fmat[MINDICE(ei-1,ej-1,maxevent)];
						//If interval, equalize ages
						if(interval==1)age2=tmpage1;
						
						if(kcost>sen1cost){
							kcost = sen1cost;
						}
						if(kcost>sen2cost){
							kcost = sen2cost;
						}
						Fmat[MINDICE(ei,ej,maxevent)]=kcost;
						sen2=sen2->getNextWithGap();
						ej++;
					}
					sen=sen->getNextWithGap();
					ei++;
                }
                dist=Fmat[MINDICE(ei-1,ej-1,maxevent)];
                if(norm>0){
					if(norm==1){
						if(maxcost>0){
							dist/=maxcost;
						}else if(dist>0){
							dist=1;
						}
					}else if(norm==2) {
						dist = (2*dist)/(maxcost+dist);
					}
                }
                matrix[MINDICE(i,j,ns)]=dist;
                matrix[MINDICE(j,i,ns)]=dist;
            }
        }
        UNPROTECT(1);
        return ans;
    }

}

