#include "prefixtree.h"
#include "eventseq.h"
#include "constraint.h"
//#include "subsequence.h"

PrefixTree::PrefixTree() {
    //ctor
}
//dtor
PrefixTree::~PrefixTree() {
    /*TreeEventMapIterator it;
    for (it = child.begin();it != child.end();it++) {
        delete it->second;
    }*/
    this->child.clearAllPointers();
}

//void PrefixTree::addSequence(Sequence *s,const double &maxGap,const double &windowSize, const double & ageMin, const double & ageMax,const double & ageMaxEnd, const int& k) {
void PrefixTree::addSequence(Sequence *s,Constraint * cst, const int& k) {

//subsequences (actually first symbol)
    if (!s->hasEvent())return;
    SequenceEventNode * e=s->getEvent();
    //Iterator to search for brother and child
    TreeEventMapIterator it;
    TreeEventNode *ten=NULL;
    double age=0;
    while (e!=NULL) {
        age+=e->getGap();
        //Only start when ageMin is reached
        if (age>cst->getageMaxBegin()) break;
        if (age>=cst->getageMinBegin()) {
            it=this->child.find(e->getType());
            if (it!=this->child.end()) {
                it->second->addSequenceInternal(s,e,cst,0,age,k, 2);
            } else if (k==1) { //Build new node only when k==1
                ten=new TreeEventNode(e->getType());
                this->child[e->getType()]=ten;
                //Rprintf("Adding event %i\n",ten->getType());
                ten->addSequenceInternal(s,e,cst,0,age,k, 2);
            }
        }
        e=e->getNext();//Get Next element
        //Rprintf("Adding starting at event %i",ten->getType());
    }
}



void PrefixTree::simplifyTree(double minSup) {

    this->child.simplifyTreeMap(minSup);
}

//Give an overview of this tree (param�tre prof==profondeur, interne)
void PrefixTree::print() {
    this->child.print(0,true);
}



int PrefixTree::countSubsequence(double minSup) {
    return this->child.countSubsequence(minSup);
}


void PrefixTree::getSubsequences(SEXP result,double * support, int *index, SEXP classname,EventDictionary * ed) {
    this->child.getSubsequences(result,support,NULL,index,0,classname,ed);
}




