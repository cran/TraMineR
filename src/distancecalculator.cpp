#include "distancecalculator.h"

void DistanceCalculator::finalizeDistanceCalculator(SEXP ptr){
	DistanceCalculator * sdo;
	sdo= static_cast<DistanceCalculator *>(R_ExternalPtrAddr(ptr));
	delete sdo;
}
