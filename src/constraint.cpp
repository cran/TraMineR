#include "constraint.h"
/*
 * constraint.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: nmuller
 */

Constraint::Constraint(const double& mg, const double& ws, const double& aminb, const double& amaxb, const double& amaxe, const double& cmethod) {
    // if (wSize==-1)wSize=DBL_MAX;
   //  if (mGap==-1)mGap=DBL_MAX;
   //  if (aMax==-1)aMax=DBL_MAX;
    // if (aMaxEnd==-1)aMaxEnd=DBL_MAX;
	if (mg==-1) {
		this->maxGap=DBL_MAX;
	}
	else {
		this->maxGap=mg;
	}
	if (ws==-1) {
		this->windowSize=DBL_MAX;
	}
	else {
		this->windowSize=ws;
	}

	this->ageMinBegin=aminb;

	if (amaxb==-1) {
		this->ageMaxBegin=DBL_MAX;
	}
	else {
		this->ageMaxBegin=amaxb;
	}
	if (amaxe==-1) {
		this->ageMaxEnd=DBL_MAX;
	}
	else {
		this->ageMaxEnd=amaxe;
	}
	if (cmethod==-1) {
		this->countMethod=1;
	}
	else {
		this->countMethod=cmethod;
	}

//	this->maxGap = mg;
//	this->windowSize = ws;
//	this->ageMinBegin = aminb;
//	this->ageMaxBegin = amaxb;
//	this->ageMaxEnd = amaxe;
}



