#include "eventdictionary.h"

EventDictionary::EventDictionary(SEXP flist) :numseq(0){
    for (int i = 0; i < length(flist); i++) {
        this->insert(std::make_pair(i+1,std::string(CHAR(STRING_ELT(flist, i)))));// add to respect R indice system
    }
    //ctor
}

EventDictionary::~EventDictionary() {
    //dtor
}

bool EventDictionary::codeExists(const int &code) const {
    return this->find(code)==this->end();
}
int EventDictionary::sprint(char * buffer, const char* start, const int&code)const {
    const_iterator it=this->find(code);
if (it!=this->end()) {
return sprintf(buffer,"%s%s",start,it->second.c_str());
}
return sprintf(buffer,"%s%i",start,code);
}

SEXP EventDictionary::getDictionary()const{
	SEXP ret;
	int s=this->size();
	//REprintf((char*)"size %i\n",s);
	PROTECT(ret = allocVector(STRSXP, s));
	for(const_iterator it=this->begin();it!=this->end();it++){
		//REprintf((char*)"Code %i=%s\n",it->first,it->second.c_str());
		if(it->first<=s){
			SET_STRING_ELT(ret, it->first-1, mkChar(it->second.c_str()));
			//REprintf((char*)"Code %i=%s\n",it->first,it->second.c_str());
		}
	}
    UNPROTECT(1);
    return ret;

}

void EventSet::add(SEXP elist){
	int * code=INTEGER(elist);
	for (int i = 0; i < length(elist); i++) {
        this->insert(code[i]);// add to respect R indice system
    }
}
