#ifndef TMRFORMAT_H_INCLUDED
#define TMRFORMAT_H_INCLUDED

#include "TraMineR.h"

extern SEXP formatSymb;
void TMRNumberFormatInit();
SEXP TMRNumberFormat(const double &num);
void TMRNumberFormatClean();


#endif // TMRFORMAT_H_INCLUDED
