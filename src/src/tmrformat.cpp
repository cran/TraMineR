#include "tmrformat.h"
	SEXP formatSymb=NULL;
	void TMRNumberFormatInit(){
		PROTECT(formatSymb=lang2(findFun(install("format"), R_GlobalEnv), R_NilValue));
	}
	SEXP TMRNumberFormat(const double &num){
		if(formatSymb==NULL){
			error(" [!!!!] TMRNumberFormat not initialized.\n");
		}
		SETCADR(formatSymb, ScalarReal(num));
		return eval(formatSymb, R_GlobalEnv);
	}
	void TMRNumberFormatClean(){
		formatSymb=NULL;
		UNPROTECT(1);
	}
