#include "TraMineR.h"
#include "distanceobject.h"

#include "LCPdistance.h"
#include "OMdistance.h"
#include "OMVIdistance.h"
#include "OMVI2distance.h"
#include "OMPerdistance.h"
#include "TWEDdistance.h"
#include "OMvdistance.h"
#include "NMSdistance.h"
#include "DHDdistance.h"
#include "NMSMSTdistance.h"
#include "NMSMSTSoftdistance.h"
#include "NMSMSTSoftdistanceII.h"
#include "NMSDURSoftdistance.h"

/**

*/

// Getting a element by its names

DistanceCalculator* getDistanceCalculatorObject(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP paramS, SEXP normS, SEXP disttypeS){
	int disttype=INTEGER(disttypeS)[0];
	TMRLOG(5, "Choosing distance type\n");
		DistanceCalculator* ds= NULL;
		if(disttype==1){
								//Base for DistanceCalculator
			ds = new OMdistance(normS, Ssequences, seqdim, lenS);
		}else if (disttype==2){
			//Setting one for LCP
			ds = new LCPdistance(normS, Ssequences, seqdim, lenS);
		}else if(disttype==4){
			ds = new DHDdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==5){
			ds = new NMSdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==6){
			ds = new NMSMSTdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==7){
			ds = new OMVIdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==8){
			ds = new OMPerdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==9){
			ds = new OMVI2distance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==10){
			ds = new OMvdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==11){
			ds = new NMSMSTSoftdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==12){
			ds = new NMSMSTSoftdistanceII(normS, Ssequences, seqdim, lenS);
		} else if(disttype==13){
			ds = new NMSDURSoftdistance(normS, Ssequences, seqdim, lenS);
		} else if(disttype==14){
			ds = new TWEDdistance(normS, Ssequences, seqdim, lenS);
		} else {
			error("Unsupported distance type");
		}
		TMRLOG(5, "Initparameters\n");
		ds->setParameters(paramS);
	return ds;
}


extern "C" {

  SEXP cstringdistance(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP paramS, SEXP normS, SEXP magicIndexS, SEXP magicSeqS, SEXP disttypeS) {
    TMRLOG(5, "Starting cstringdistance\n");
    //Objet R, matrice des distances (objet dist)
    //Indices, avec s pour séquences
    int is, js;
    //longueur des séquences m=i, n=j
    //int m, n;
    DistanceObject* distObj = new DistanceObject(magicIndexS, magicSeqS);

    int nseq= INTEGER(seqdim)[0];
    TMRLOG(5, "Choosing distance type\n");
    DistanceCalculator* ds= getDistanceCalculatorObject(Ssequences, seqdim, lenS, paramS, normS, disttypeS);

    // Ensure correct memory management
    SEXP workers;

    PROTECT(workers = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(workers,0, distanceObjectFactory(distObj));
    SET_VECTOR_ELT(workers,1, DistanceCalculator::distanceCalculatorFactory(ds));
    TMRLOG(5, "Initparameters finished\n");

    //starting store index
    //int i_start, j_start, i_end, j_end, i_index, j_index, base_index;
    //Pour chaque séquence i
    // double perc=0.05;
    // double totcompute=(double)nseq*(nseq-1);
    // double currentperc=0;
    // REprintf(" [>] Progress (#=2.5%%): ");
    double cmpres=0;
    for (is=0;is<nseq;is++) {
      //toutes les distances intra-groupes=0
      R_CheckUserInterrupt();
      distObj->setDistance(is,is, 0);
      for (js=is+1;js<nseq;js++) {
        cmpres = ds->distance(is,js);
        TMRLOG(5,"cmpres = %d %d => %f \n",(1+is),(1+js), cmpres);
        //return Fmat;
        //Same for j
        distObj->setDistance(is,js, cmpres);
        //result[MINDICE(is,js,nseq)]=result[MINDICE(js,is,nseq)]=cmpres;
      }//end js
      // currentperc+=(double)(nseq-is-1);
      // while(currentperc/totcompute>perc){
      // REprintf("#");
      // perc+=0.025;
      // }
    }
    SEXP ans = distObj->getDistObject();
    UNPROTECT(2);
    return ans;
  }

	SEXP cstringrefseqdistance(SEXP Ssequences, SEXP seqdim, SEXP lenS, SEXP paramS, SEXP normS, SEXP disttypeS, SEXP refseqS) {
    TMRLOG(5, "Starting cstringdistancerefseq\n");
    //Objet R, matrice des distances (objet dist)
    //Indices, avec s pour séquences

    //longueur des séquences m=i, n=j
    //int m, n;

    int nseq= INTEGER(seqdim)[0];
    int rseq= INTEGER(refseqS)[0]-1;

    TMRLOG(5, "Choosing distance type\n");
    DistanceCalculator* ds= getDistanceCalculatorObject(Ssequences, seqdim, lenS, paramS, normS, disttypeS);

    // Ensure correct memory management
    SEXP workers, ans;
    PROTECT(ans = allocVector(REALSXP, nseq));

    PROTECT(workers = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(workers,0, DistanceCalculator::distanceCalculatorFactory(ds));
    TMRLOG(5, "Initparameters finished\n");
    double * distances=REAL(ans);
    //starting store index
    //int i_start, j_start, i_end, j_end, i_index, j_index, base_index;
    //Pour chaque séquence i
    // double perc=0.05;
    // double totcompute=(double)nseq*(nseq-1);
    // double currentperc=0;
    // REprintf(" [>] Progress (#=2.5%%): ");
    double cmpres=0;
    for (int is=0;is<nseq;is++) {
      //toutes les distances intra-groupes=0
      R_CheckUserInterrupt();
      if(is==rseq){
        cmpres=0;
      } else {
        cmpres = ds->distance(is,rseq);
      }
      distances[is]=cmpres;
      TMRLOG(5,"cmpres = %d %d => %f \n",(1+is),(1+rseq), cmpres);
    }
    UNPROTECT(2);
    return ans;
    }
}
