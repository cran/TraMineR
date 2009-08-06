#include "treeeventnode.h"
#include "treeeventmap.h"
#include "constraint.h"

int TreeEventNode::nodeCount=0;
int TreeEventNode::getNodeCount() {
    return nodeCount;
}

void TreeEventNode::getSubsequences(SEXP result,int *  isupport, Sequence *s, int *index,const double &step, SEXP classname,EventDictionary * ed) {
    this->brother.getSubsequences(result,isupport,s,index,step, classname,ed);
    this->child.getSubsequences(result,isupport,s,index,step+1, classname,ed);

}


//Type of this event
/**
 *TreeEventNode
 Prefix Tree base node
*/
//Ctor
TreeEventNode::TreeEventNode(const int& t):type(t),support(0),lastID(-1) {
    nodeCount++;
}
//Dtor
TreeEventNode::~TreeEventNode() {
    TreeEventMapIterator it;

    nodeCount--;

    this->brother.clearAllPointers();
    this->child.clearAllPointers();
    //REprintf("Node deleted (DTOR finished)%i\n",nodeCount);
}



//Main function to build the tree
void TreeEventNode::addSequenceInternal(Sequence *s, SequenceEventNode * en, Constraint * cst, const double &gapConsumed, const double& currentAge, const int& k, const int&currentK) {
    //If we already reached this point with the specified sequence, don't increment
    /*if (this->lastID!=s->getIDpers()) {
        this->support++;//
        this->lastID=s->getIDpers();
    }*/

	// If we count several by sequences, or the last element added was from another sequence, we had this element to support.
	if (cst->getcountMethod()==2 || this->lastID!=s->getIDpers()) {
        this->support++;//
        this->lastID=s->getIDpers();
	}
    if (!en->hasNext())return;
    if (currentK>k)return;
    //Current Gap for next subsequence search
    double currentGap=0;
    //next subsequences
    SequenceEventNode *n=en;
    //Iterator to search for brother and child
    TreeEventMapIterator it;
    TreeEventNode *ten=NULL;
    while (n->hasNext()) {
        n=n->getNext();//Get Next element
        currentGap+=n->getGap(); //Increment current gap

        //3 terminations conditions
        if (	gapConsumed+currentGap>cst->getwindowSize() //current window Size too big
                ||currentGap>cst->getmaxGap()	//current gap too big
                ||currentGap+currentAge>cst->getageMaxEnd() //current age too big
           ) {
            break;
        }
        //IF the gap is >0, then we have children
        if (currentGap>0) {
            it=this->child.find(n->getType());//Search for child
            if (it!=this->child.end()) { //if exist
                ten=it->second;
            } else if (k==currentK) { //Build new child
                ten=new TreeEventNode(n->getType());
                this->child[n->getType()]=ten;
            } else {
                ten=NULL;
            }
        } else { //currentGap==0
            it=this->brother.find(n->getType());
            if (it!=this->brother.end()) {
                ten =it->second;
            } else if (k==currentK) { //Build new brother
                ten=new TreeEventNode(n->getType());
                this->brother[n->getType()]=ten	;
            } else {
                ten=NULL;
            }
        }
        if (ten!=NULL)ten->addSequenceInternal(s,n,cst,gapConsumed+currentGap, currentAge+currentGap,k,currentK+1);//Increment support and continue

    }//end while(hasNext())
}

void TreeEventNode::simplifyTree(int minSup) {
    this->brother.simplifyTreeMap(minSup);
    this->child.simplifyTreeMap(minSup);
}

void TreeEventNode::print(const int & prof, const bool& isbrother) {
    for (int i=0;i<prof;i++) {
        Rprintf((char*)"   ");
    }
    if (isbrother) {
        Rprintf((char*)"|--(%i:%i)[%i]\n",this->type,this->support,this->lastID);
    } else {
        Rprintf((char*)"|__(%i:%i)[%i]\n",this->type,this->support,this->lastID);
    }
    this->brother.print(prof+1,isbrother);
    //print each child branches
    this->child.print(prof+1,isbrother);

}

int TreeEventNode::countSubsequence(int minSup) {

    return 1+this->brother.countSubsequence(minSup)+this->child.countSubsequence(minSup);
}
void TreeEventNode::clearSupport() {
    this->support=0;
    this->child.clearSupport();
    this->brother.clearSupport();
}
