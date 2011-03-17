#include "eventseq.h"
#include<R.h>
#include <Rinternals.h>
using namespace std;

/** Sequence finalizer, used by R to free memory
*/
void finalizeSequence(SEXP ptr) {

    Sequence *s;
    ASSIGN_TMRSEQ_TYPE(s,ptr);
    delete s;
}
/** CLASS SequenceEventNode
  Contain one event in an indiviudal sequence, reponsible for deleting next event in sequence
*/

//recursive add
void SequenceEventNode::addEvent(const int &e,const double &t) {
    if (this->hasNext()){
    	if(this->next->greaterThan(e, t-this->gap)){
			SequenceEventNode * s= new SequenceEventNode(e,t-this->gap);
			this->next->gap-=t-this->gap;
			s->setNext(this->next);
			this->next=s;
    	}else this->next->addEvent(e,t-this->gap);
    }
    else {
        this->next=new SequenceEventNode(e,t-this->gap);
    }
}

void Sequence::addEvent(const int &e,const double &t) {
    if (this->hasEvent()){
    	if(this->event->greaterThan(e,t)){
			this->event->setGap(this->event->getGap()-t);
    		SequenceEventNode * s=new SequenceEventNode(e,t);
    		s->setNext(this->event);
    		this->event=s;
    	}else{
			this->event->addEvent(e,t);
    	}
    } else {
        this->event=new SequenceEventNode(e,t);
    }
}
SequenceEventNode * SequenceEventNode::copy() {
    SequenceEventNode *s=new SequenceEventNode(this->type, this->gap);
    if (this->hasNext())s->next=this->next->copy();
    return s;
}
Sequence * Sequence::copy() {
    Sequence *s=new Sequence(this->idpers,this->dict);
    if (this->hasEvent())s->event=this->event->copy();
    return s;
}
///CLASS Sequence
///Represent an individual sequence
///CTor
Sequence::Sequence(const int&id, EventDictionary* ed):dict(ed), obsTime(-1), weight(1) { //personnal time 0, type =0 (root)
//    this->ns=NULL;
	this->dict->addSequence();
    this->idpers=id;
    this->event=NULL;
}
///Dtor clear all events and next sequences.
Sequence::~Sequence() {
    //if (ns!=NULL)delete ns;
    if (this->event!=NULL) delete event;
    this->dict->removeSequence();
    if(this->dict->shouldDelete())delete this->dict;
}

int Sequence::sprint(char * buffer) {
    int n2=0;
    int n=0;
    //if(!this->isGeneric())n = sprintf(buffer, (char*)"[%i] ",this->idpers);
    //Rprintf((char*)"Current buffer %s\n",buffer);
    if (n==-1)return -1;
    if (this->hasEvent()) {
        n2=this->event->sprint(buffer, n,true, !this->isGeneric(), (*this->dict),this->obsTime);
        if (n2==-1)return -1;
    }
    return n+n2;

}
void Sequence::print() {
    char buffer[TMR_STRING_BUFFER_SIZE];
    buffer[0]='\0';
    int r=this->sprint(buffer);
    //Rprintf((char *)"%s %i",buffer,r);
    if (r>0)	Rprintf((char *)"%s\n",buffer);
    else Rprintf((char *)"Error %i\n",r);
}
int SequenceEventNode::sprint(char * buffer, int index, const bool& start, const bool &printGap, const EventDictionary& ed, const double & remainingTime) {
    int tmp=0;
    if (start) {
        if (this->gap>0&&printGap) {
            tmp=sprintf(&buffer[index],(char*)"%.2f-",this->gap);
            if (tmp==-1)return -1;
            index+=tmp;
            tmp=0;
            tmp=ed.sprint(&buffer[index],"(",this->type);
        } else {
            tmp=ed.sprint(&buffer[index],"(",this->type);
        }


    } else if (this->gap>0) {
        if (printGap) {
            tmp=sprintf(&buffer[index],(char*)")-%.2f-",this->gap);
            if (tmp==-1)return -1;
            index+=tmp;
            tmp=0;
            tmp=ed.sprint(&buffer[index],"(",this->type);
        } else {
            tmp=ed.sprint(&buffer[index],")-(",this->type);
        }
    } else {
        tmp=ed.sprint(&buffer[index],",",this->type);
    }
    if (tmp==-1)return -1; //An error occured
    index+=tmp; //Increment index
    //Print next elements
    if (this->hasNext()) {
        tmp=this->next->sprint(buffer, index, false, printGap,ed, remainingTime-this->gap);
    } else {
    	if(remainingTime>0){
    		tmp=sprintf(&buffer[index],(char*)")-%.2f", remainingTime-this->gap);
    	}
    	else{
    		tmp=sprintf(&buffer[index],(char*)")");
    	}

    }
    //Rprintf((char*)"Current buffer %s => %i\n",buffer, index);
    if (tmp==-1)return -1; //An error occured
    else return tmp+index;
}



double Sequence::first_occurence(Sequence * s, const double &maxGap, const double& windowSize, const double & ageMin, const double & ageMax, const double & ageMaxEnd) {
    if (!this->hasEvent()||!s->hasEvent()) return -1;
    double age=0;
    SequenceEventNode *sen=s->getEvent();
    while (sen!=NULL) {
        age+=sen->getGap();
        if (age>ageMax)return -1;
        if (age>=ageMin&&this->event->checkType(sen)&&this->event->count(sen,maxGap, windowSize,ageMaxEnd,0,age)>0) {
            return age;
        }
        sen=sen->getNext();
    }
    return -1;
}
int Sequence::count(Sequence * s, const double &maxGap, const double& windowSize, const double & ageMin, const double & ageMax, const double & ageMaxEnd) {
    if (!this->hasEvent()||!s->hasEvent()) return 0;
    SequenceEventNode * sen=s->getEvent();
    double age=0;
    int c=0;
    while (sen!=NULL) {
        age+=sen->getGap();
        if (age>ageMax)break;
        if (age>=ageMin&&this->event->checkType(sen)) {
            c+=this->event->count(sen, maxGap, windowSize,ageMaxEnd, 0,age);
        }
        sen=sen->getNext();
    }
    return c;
}
//Internal method, check, how many starting from this point already checked (we are looking for the rest of subsequence
int SequenceEventNode::count(SequenceEventNode * s, const double &maxGap, const double& windowSize, const double & ageMaxEnd, const double& gapConsumed, const double& currentAge) {
    int c=0;
    if (!this->hasNext())return 1;
    //Rprintf("Checked %i and %i %f\n", this->type,s->type, this->next->gap);
    SequenceEventNode * sen=s->getNext();
    if (this->next->gap==0) {
        //Rprintf("Equality %i and %i\n", this->next->type,sen->type);
        while (sen!=NULL&&sen->gap==0) {
            //Rprintf("Equality %i and %i\n", this->next->type,sen->type);
            if (this->next->checkType(sen)) {
                c+=this->next->count(sen, maxGap, windowSize, ageMaxEnd,gapConsumed,currentAge);
            }
            sen=sen->getNext();
        }
    } else {
        while (sen!=NULL&&sen->gap==0) {//Looking for next gap
            sen=sen->getNext();
            //Rprintf("Skipping %i\n", sen->type);
        }
        if (sen==NULL)return 0; //we didn't match
        double g=0;
        while (sen!=NULL) {
            //	Rprintf("Step %i and %i\n", this->next->type,sen->type);
            g+=sen->gap;
            //Time constraints, we don't need to look deeper
            if (g>maxGap||(g+gapConsumed)>windowSize)return c;
            if ((currentAge+g)>ageMaxEnd)return c;
            if (this->next->checkType(sen)) {
                c+=this->next->count(sen, maxGap, windowSize, ageMaxEnd,g+gapConsumed,g+currentAge);
            }
            sen=sen->getNext();
        }
    }
    return c;
}
bool Sequence::contain(const EventSet& es, const bool& exclude){
	if (!this->hasEvent()) return false;
    SequenceEventNode * sen=this->getEvent();
    while (sen!=NULL) {
		if(es.contain(sen->getType())){
			if(!exclude) return true;
		}else if(exclude)return false;
        sen=sen->getNext();
    }
    return exclude;
}
/*bool Sequence::notContain(const EventSet& es){
	if (!this->hasEvent()) return true;
    SequenceEventNode * sen=this->getEvent();
    while (sen!=NULL) {
		if(es.contain(sen->getType()))return false;
        sen=sen->getNext();
    }
    return true;
}*/
