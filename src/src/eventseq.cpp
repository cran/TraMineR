#include "eventseq.h"
#include <sstream>

#include "tmrformat.h"
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

string Sequence::sprint() {
	ostringstream oss;
	//oss.precision(2);
    //if(!this->isGeneric())n = sprintf(buffer, (char*)"[%i] ",this->idpers);
    //Rprintf((char*)"Current buffer %s\n",buffer);
    if (this->hasEvent()) {
        this->event->sprint(oss, true, !this->isGeneric(), (*this->dict), this->obsTime);
    }
    return oss.str();

}
void Sequence::print() {
	TMRNumberFormatInit();
    string r=this->sprint();
    //Rprintf((char *)"%s %i",buffer,r);
    REprintf((char *)"%s\n",r.c_str());
	TMRNumberFormatClean();
}
void SequenceEventNode::sprint(ostringstream &oss, const bool& start, const bool &printGap, const EventDictionary& ed, const double & remainingTime) {
    if (start) {
        if (this->gap>0&&printGap) {
			SEXP gg;
			PROTECT(gg=asChar(TMRNumberFormat(gap)));
			oss << CHAR(gg) << "-(" << ed.find(this->type)->second;
			UNPROTECT(1);
            //tmp=sprintf(&buffer[index],(char*)"%.2f-",this->gap);
        } else {
            //tmp=ed.sprint(&buffer[index],"(",this->type);
			oss << "(" << ed.find(this->type)->second;
        }


    } else if (this->gap>0) {
        if (printGap) {
			SEXP gg;
			PROTECT(gg=asChar(TMRNumberFormat(gap)));
			oss << ")-" << CHAR(gg)<< "-(" << ed.find(this->type)->second;
			UNPROTECT(1);
        } else {
			oss << ")-(" << ed.find(this->type)->second;
			
        }
    } else {
		oss << "," << ed.find(this->type)->second;
    }
    if (this->hasNext()) {
        this->next->sprint(oss, false, printGap, ed, remainingTime-this->gap);
    } else {
    	if(remainingTime>0){
			SEXP gg;
			PROTECT(gg=asChar(TMRNumberFormat(remainingTime-this->gap)));
			oss << ")-" << CHAR(gg);
			UNPROTECT(1);
    	}
    	else{
    		oss << ")";
    	}

    }
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

// count the number of subsequences in assigned sequence * s

int Sequence::count(Sequence * s, const double &maxGap,
		    const double &windowSize, const double &ageMin, 
		    const double &ageMax, const double &ageMaxEnd,
		    const int &cMethod) 
{
  if (!this->hasEvent()||!s->hasEvent()) return 0;
 // define SequenceEventNode from first event in entered sequence s 
  SequenceEventNode * sen=s->getEvent();
  // define age of first event
  double age = 0;
  int type = 0;
  int c = 0;
  int cs;
  
  switch (cMethod)
    {
    case 1: // COBJ method --------------------------------- //
      {
	sen = s->getEvent();
	age = 0;
	while (sen!=NULL)
	  {
	    // search in sequence for the first event of the subsequence
	    age+=sen->getGap();
	    if (age>ageMax) break; // stop if age is exceeded
	    if (age>=ageMin && this->event->checkType(sen))  // if event found
	      { 
		cs = this->event->count(sen,maxGap,windowSize,
					ageMaxEnd,0,age);
		if (cs>0) c = 1; // count 1 if subsequence occures in sequence
	      }
	    sen = sen->getNext();
	  }
      }
      break;
    case 2: // CDIST_0 method ------------------------------ //
      {
	sen = s->getEvent();
	age = 0;
	while (sen!=NULL)
	  {
	    // search in sequence for the first event of the subsequence
	    age+=sen->getGap();
	    if (age>ageMax) break; // stop if age is exceeded
	    if (age>=ageMin && this->event->checkType(sen))  // if event found
	      { 
		c+=this->event->count(sen,maxGap,windowSize,ageMaxEnd,0,age);
	      }
	    sen = sen->getNext();
	  }
      }
      break;
    case 3: // CWIN method --------------------------------- //
      {
	double startWin, stopWin;
	int lWin;
	if (ageMin==-DBL_MAX) 
	  {
	    sen = s->getEvent();
	    startWin = sen->getGap()-windowSize;
	  } else {
	  startWin = ageMin;
	}
	if (ageMaxEnd==DBL_MAX)
	  {
	    sen = s->getEvent();
	    age = 0;
	    while (sen!=NULL)
	      {
		age+=sen->getGap();
		sen = sen->getNext();
	      }
	    stopWin = age+windowSize;
	    TMRLOG(2,"Set endMaxTime to: %.1f\n",stopWin);
	  } else {
	  stopWin = ageMaxEnd;
	}
	lWin = round(stopWin-startWin-windowSize+1);
	double *tWin = new double[lWin];
	TMRLOG(3,"Windows: \n");
	for (int i=0;i<lWin;i++)
	  {
	    tWin[i] = startWin+i;
	    if ((tWin[i]+windowSize)<ageMaxEnd)
	      {
	    	TMRLOG(2,"[%.1f, %.1f]\n",tWin[i],tWin[i]+windowSize);
	      } else {
	      TMRLOG(2,"[%.1f, %.1f]\n",tWin[i],ageMaxEnd);
	    }
	  }
	int *Win =new int[lWin]; // initialize the one step
	std::fill_n(Win,lWin,0); // fill the array with zeros
	sen = s->getEvent();
	age = 0;
	while (sen!=NULL)
	  {
	    // search in sequence for the first event of the subsequence
	    age+=sen->getGap();
	    if (age>ageMax) break; // stop if age is exceeded
	    if (age>=ageMin && this->event->checkType(sen))  // if event found
	      { 
		c+= this->event->count3(sen,maxGap,windowSize,
					ageMaxEnd,0,age,Win,tWin,lWin);
	      }
	    sen = sen->getNext();
	  }
	// Special count mechanism for CWIN method: the number of occurences
	// are saved in vector WIN. A Win[i]=1 means that there is at 
	// least one subsequence within [tWin[i],tWin[i]+windowSpan)
	c = 0;
	for (int i=0;i<lWin;i++)
	  {
	    c = c+Win[i];
	  }
	delete [] tWin;
	delete [] Win;
      }
      break;

    case 4: // CMINWIN method ------------------------------ //
      {
	double minWin = DBL_MAX; // required in case of countMethod = 4
	double minWinOld = minWin;
	sen = s->getEvent();
	age = 0;
	while (sen!=NULL)
	  {
	    // search in sequence for the first event of the subsequence
	    age+=sen->getGap();
	    if (age>ageMax) break; // stop if age is exceeded
	    if (age>=ageMin && this->event->checkType(sen))  // if event found
	      { 
		TMRLOG(3,"CMinWin %.1f",maxGap);
		cs=this->event->count4(sen,maxGap,windowSize,
				       ageMaxEnd,0,age,minWin);
		// add one if the new found subsequence has the equal minimum
		// time window as the previous one(s)
		if (minWinOld==minWin && cs>0) c++;
		// if the new found subsequence has a smaller time window
		// as the previous, reset the counter
		if (minWinOld>minWin)
		  {
		    c = 1;
		    minWinOld = minWin;
		  } 
	      }
	    sen = sen->getNext();
	  }
      }
      break;

    case 5: // CDIST method -------------------------------- //
      {
	int lSen = 0; // number of elements in sequence
	while (sen!=NULL)
	  {
	    lSen++;
	    sen=sen->getNext();
	  }
	double *tSen = new double[lSen]; // event times
	int *typeSen = new int[lSen]; // event types
	sen = s->getEvent();
	for (int i=0;i<lSen;i++)
	  {
	    age+=sen->getGap();
	    tSen[i] = age;
	    type = sen->getType();
	    typeSen[i] = type;
	    sen=sen->getNext();
	  }
	int *flagSen = new int[lSen]; // flag previously counted events
	std::fill_n(flagSen,lSen,0); // fill the array with zeros
	sen = s->getEvent();
	age = 0;
	while (sen!=NULL)
	  {
	    // search in sequence for the first event of the subsequence
	    age+=sen->getGap();
	    if (age>ageMax) break; // stop if age is exceeded
	    if (age>=ageMin && this->event->checkType(sen))  // if event found
	      { 
		c+=this->event->count5(sen,maxGap,windowSize,ageMaxEnd,0,age,
				       typeSen,tSen,lSen,flagSen);
	      }
	    sen = sen->getNext();
	  }
	delete [] tSen;
	delete [] typeSen;
	delete [] flagSen;
      }
      break;
    case 6: // age of appearance --------------------------- // 
      {
	c+=this->first_occurence(s,maxGap,windowSize, 
				 ageMin,ageMax,ageMaxEnd);
	sen = NULL;
      }
      break;
    default: // -------------------------------------------- //
      {
	sen = s->getEvent();
	age = 0;
	while (sen!=NULL)
	  {
	    // search in sequence for the first event of the subsequence
	    age+=sen->getGap();
	    if (age>ageMax) break; // stop if age is exceeded
	    if (age>=ageMin && this->event->checkType(sen))  // if event found
	      { 
		cs = this->event->count(sen,maxGap,windowSize,
					ageMaxEnd,0,age);
		if (cs>0) c = 1; // count 1 if subsequence occures in sequence
	      }
	    sen = sen->getNext();
	  }
      }
      break;
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
// --------------------------------------------------------- //
// Internal method for CWIN counting method
// Reto Buergin, June 2011
// --------------------------------------------------------- //

/**
   Win : Out
   tWin : Out
   lWin : Our
 */

int SequenceEventNode::count3(SequenceEventNode *s, 
			      const double &maxGap, 
			      const double &windowSize, 
			      const double &ageMaxEnd, 
			      const double &gapConsumed, 
			      const double &currentAge,
			      int *Win, double *tWin, const int &lWin) 
{
  int c = 0; // counter
  if (!this->hasNext()) // end of sequence
    {
      TMRLOG(3,"Found a complete subsequence in window: ");
      TMRLOG(3,"[%.1f, %.1f]\n",currentAge-gapConsumed,currentAge);
      for (int i=0;i<lWin;i++)
      	{
	  if ((tWin[i]+windowSize)<=ageMaxEnd)
	    {
	      TMRLOG(3,"Check interval [%.1f, %.1f]:\t",
		     tWin[i],tWin[i]+windowSize);
	      if (currentAge-gapConsumed>=tWin[i] &&
		  currentAge<=tWin[i]+windowSize &&
		  Win[i]==0)
		{
		  TMRLOG(3,"TRUE\n"); 
		  Win[i] = 1;
		} else TMRLOG(3,"FALSE\n");
	    } else {
	    TMRLOG(3,"Check interval [%.1f, %.1f]:\t",
		   tWin[i],ageMaxEnd);
	    if (currentAge-gapConsumed>=tWin[i] &&
		currentAge<=ageMaxEnd &&
		  Win[i]==0)
	      {
		TMRLOG(3,"TRUE\n"); 
		Win[i] = 1;
	      } else TMRLOG(3,"FALSE\n");
	  }
      	}
      return 1;
    }
  SequenceEventNode * sen=s->getNext();
  // if next event in subsequence takes place at the same time
  if (this->next->gap==0) 
    {
      while (sen!=NULL&&sen->gap==0) 
	{
	  if (this->next->checkType(sen)) 
	    {
	      c+=this->next->count3(sen,maxGap,windowSize, 
				    ageMaxEnd,gapConsumed,currentAge,
				    Win,tWin,lWin);
	    }
	  sen=sen->getNext();
	}
    } else {
    while (sen!=NULL&&sen->gap==0) 
      {//Looking for next gap
	sen=sen->getNext();
	//Rprintf("Skipping %i\n", sen->type);
      }
    if (sen==NULL) return 0; //we didn't match
    double g = 0;
    // if next event in subsequence has a gap to the current
    while (sen!=NULL) 
      {
	g+=sen->gap; // get the gap to the current event
	// time constraints, we don't need to look deeper
	if (g>maxGap||(g+gapConsumed)>windowSize) return c;
	if ((currentAge+g)>ageMaxEnd) return c;
	if (this->next->checkType(sen)) 
	  {
	    c+=this->next->count3(sen,maxGap,windowSize,ageMaxEnd,
				  g+gapConsumed,g+currentAge,
				  Win,tWin,lWin);
	  }
	sen=sen->getNext();
      }
  }
  return c;
}

// --------------------------------------------------------- //
// Internal method for CMinWin counting method
// Reto Buergin, June 2011
// --------------------------------------------------------- //

/**
   minWin : Out
 */

int SequenceEventNode::count4(SequenceEventNode *s, 
			      const double &maxGap, 
			      const double &windowSize, 
			      const double &ageMaxEnd, 
			      const double &gapConsumed, 
			      const double &currentAge,
			      double &minWin) 
{ 
  int c = 0; // counter
  if (!this->hasNext())
    {
      return 1;
    }
  SequenceEventNode *sen = s->getNext();
  // if next event in subsequence takes place at the same time
  if (this->next->gap==0) 
    {
      while (sen!=NULL&&sen->gap==0) 
	{
	  if (this->next->checkType(sen))
	    { // found an event of subsequence in sequence
	      c+=this->next->count4(sen,maxGap,windowSize, 
				    ageMaxEnd,gapConsumed,currentAge,minWin);
	      if (!this->next->hasNext())
		{ // found the complete subsequence in sequence
		  if (gapConsumed<minWin)
		    { // renew minWin if required
		      if (minWin<DBL_MAX) 
			{
			  TMRLOG(3,"Delete current minWin\n");
			}
		      minWin = gapConsumed;
		      TMRLOG(3,"Found new minWin: [%.1f, %.1f] (length %.1f)\n",currentAge-gapConsumed,currentAge,minWin);
		    }
		}
	    }
	sen=sen->getNext();
	}
    } else {
    while (sen!=NULL&&sen->gap==0) 
      { // looking for next gap
	sen=sen->getNext();
      }
    if (sen==NULL) return 0; //we didn't match
    double g = 0;
    // if next event in subsequence has a gap to the current
    while (sen!=NULL) 
      {
	g+=sen->gap; // get the gap to the current event
	// leave if consumed gap larger than current minWin
	if (g+gapConsumed>minWin) return c;
	// time constraints, we don't need to look deeper
	if (g>maxGap||(g+gapConsumed)>windowSize) return c;
	if ((currentAge+g)>ageMaxEnd) return c;
	if (this->next->checkType(sen)) 
	  {  // found an event of subsequence in sequence
	    c+=this->next->count4(sen, maxGap, windowSize, ageMaxEnd,
				  g+gapConsumed,g+currentAge,minWin);
	    if (!this->next->hasNext())
	      { // found the complete subsequence in sequence
		if (g+gapConsumed<minWin)
		  { // renew minWin if required
		    if (minWin<DBL_MAX) TMRLOG(3,"Delete current minWin\n");
		    minWin = g+gapConsumed;
		    TMRLOG(3,"Found new minWin: [%.1f, %.1f] (length %.1f)\n",currentAge-gapConsumed,currentAge+g,minWin);
		    TMRLOG(3,"New minWin: %f\n",minWin);
		  }
	      }
	  }
	sen=sen->getNext();
      }
  }
  return c;
}

// --------------------------------------------------------- //
// Internal method for CDIST counting method
// Reto Buergin, June 2011
// --------------------------------------------------------- //

/**
   typeSen : Out
   tSen : Out
   flagSen : Out
 */

int SequenceEventNode::count5(SequenceEventNode *s, 
			      const double &maxGap, 
			      const double &windowSize, 
			      const double &ageMaxEnd, 
			      const double &gapConsumed, 
			      const double &currentAge,
			      const int *typeSen, const double *tSen,
			      const int &lSen, int *flagSen) 
{
  int c = 0; // counter
  int ind, indc;
  ind = 0;
  indc = 0;
  while (ind==0 && indc<lSen)
    {
      if (tSen[indc]==currentAge &&
	  typeSen[indc]==this->getType())
	{
	  ind = indc;
	}
      indc++;
    }
  if (flagSen[ind]!=0)
    {
      TMRLOG(3,"Event (%i,%.1f) was already flagged, return\n",
	     typeSen[ind],tSen[ind]);
      return 0;
    }
  TMRLOG(3,"Flag event (%i,%.1f)\n",typeSen[ind],tSen[ind]);
  flagSen[ind] = 1;
  if (!this->hasNext()) 
    {
      for (int i=0;i<lSen;i++)
	{
	  if (flagSen[i]==1) 
	    {
	      flagSen[i] = 2;
	      TMRLOG(3,"Flag event (%i,%.1f) definitly\n",
		     typeSen[i],tSen[i]);
	    }
	}
      TMRLOG(3,"Leave count5() function at beginning\n");
      return 1;
    }
  SequenceEventNode * sen=s->getNext();
  // if next event in subsequence takes place at the same time
  if (this->next->gap==0) 
    {
      while (sen!=NULL&&sen->gap==0) 
	{
	  if (this->next->checkType(sen)) 
	    {
	      TMRLOG(3,"Found an event within sequence\n");
 	      ind = 0;
	      indc = 0;
	      while ((ind==0) && (indc<lSen))
		{
		  if (tSen[indc]==currentAge &&
		      typeSen[indc]==this->next->getType())
		    {
		      ind = indc;
		    }
		  indc++;
		}
	      if (flagSen[ind]==0)
		{
		  c+=this->next->count5(sen,maxGap,windowSize, 
					ageMaxEnd,gapConsumed,currentAge,
					typeSen,tSen,lSen,flagSen);
		  if (flagSen[ind]==2)
		    {
		      TMRLOG(3,"Finish count5() function\n");
		      return c;
		    }
		} else {
		TMRLOG(3,"Event (%i,%.1f) was already flagged\n",
		       typeSen[ind],tSen[ind]);
	      }
	    }
	  sen=sen->getNext();
	}
    } else { // skipping to next event in sequence with time gap
    while (sen!=NULL&&sen->gap==0) 
      { //Looking for next gap
	sen=sen->getNext();
	//Rprintf("Skipping %i\n", sen->type);
      }
    if (sen==NULL) // if no success, erase temporarily flagged events
      {
	for (int i=0;i<lSen;i++)
	  {
	    if (flagSen[i]==1) flagSen[i] = 0;
	  }
	return c; //we didn't match
      }
    double g = 0;
    // if next event in subsequence has a gap to the current
    while (sen!=NULL) 
      {
	g+=sen->gap; // get the gap to the current event
	// time constraints, we don't need to look deeper
	if (g>maxGap||(g+gapConsumed)>windowSize ||
	    (currentAge+g)>ageMaxEnd)
	  {
	    for (int i=0;i<lSen;i++) // delete all temporary flagged events
	      {
		if (flagSen[i]==1) flagSen[i] = 0;
	      }
	    return c;
	  }
	if (this->next->checkType(sen)) 
	  {
	    TMRLOG(3,"Found an event within sequence\n");
	    ind = 0;
	    indc = 0;
	    while (ind==0 && indc<lSen)
	      {
	    	if (tSen[indc]==currentAge+g &&
		    typeSen[indc]==this->next->getType())
	    	  {
	    	    ind = indc;
	    	  }
	    	indc++;
	      }
	    if (flagSen[ind]==0)
	      {
		c+=this->next->count5(sen, maxGap, windowSize, ageMaxEnd,
				      g+gapConsumed,g+currentAge,
				      typeSen,tSen,lSen,flagSen);
		if (flagSen[ind]==2)
		  {
		    TMRLOG(3,"Finish count5() function\n");
		    return c;
		  }
	      } else {
	      TMRLOG(3,"Event (%i,%.1f) was already flagged\n",
		     typeSen[ind],tSen[ind]);
	      }
	  }
	sen=sen->getNext();
      }
  }
  for (int i=0;i<lSen;i++) // delete all temporary flagged events
    {
      if (flagSen[i]==1) flagSen[i] = 0;
    }
  TMRLOG(3,"Leave count5 function at the end\n");
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
