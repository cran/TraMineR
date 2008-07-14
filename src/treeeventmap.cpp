#include "treeeventmap.h"
#include "treeeventnode.h"

void TreeEventMap::clearAllPointers() {
    TreeEventMapIterator it;
    for (it = this->begin();it != this->end();it++) {
        delete it->second;
    }
    this->clear();
}
void TreeEventMap::clearSupport() {
    TreeEventMapIterator it;
    for (it = this->begin();it != this->end();it++) {
        it->second->clearSupport();
    }
}
void TreeEventMap::simplifyTreeMap(const int &minSup) {
    TreeEventMapIterator it=this->begin(), it2;
    TreeEventNode * n=NULL;
    while (it!=this->end()) {
        if (it->second->getSupport()<minSup) {
            it2=it;
            it2++;
            n=it->second;
            //REprintf("Deleting brother (Simplify) %i\n",TreeEventNode::getNodeCount());
            delete n; //Delete (free)
            this->erase(it);
            it = it2;
        } else {
            //REprintf("(SimplifySub) %i=>%i\n",this->type,it->second->type);
            it->second->simplifyTree(minSup);
            it++;
            //REprintf("(SimplifySub Finished) %i=>%i\n",this->type,it->second->type);
        }
    }

}

//Give an interpretation of the tree:
//- denote a brother
// underscore for child
void TreeEventMap::print(const int & prof, const bool& isbrother) {
    TreeEventMapIterator it;
    //print each brother branches
    for (it = this->begin();it != this->end();it++) {
        it->second->print(prof,isbrother);
    }
}

void TreeEventMap::getSubsequences(SEXP result,int * support, Sequence *s2, int *index,const double &step, SEXP classname,EventDictionary * ed) {
    SEXP tmpseq;
    TreeEventMapIterator it;
    Sequence * s=NULL;
    bool copy=(s2!=NULL);
    for (it = this->begin();it != this->end();it++) {
        if (copy)s=s2->copy();
        else s=new Sequence(-1,ed);
        s->addEvent(it->second->getType(),step);
        tmpseq = makeTMRSequence(s, classname);
        SET_VECTOR_ELT(result,*index,tmpseq);
        support[*index]=it->second->getSupport();
        *index=(*index)+1;
        it->second->getSubsequences(result,support,s,index,step,classname,ed);
    }
}
int TreeEventMap::countSubsequence(int minSup) {
    TreeEventMapIterator it;
    int count=0;
    for (it = this->begin();it != this->end();it++) {
        count+=it->second->countSubsequence(minSup);
    }
    return count;
}

