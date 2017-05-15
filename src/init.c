#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .C calls */
extern void cLCP(int *, int * , double *, int *);
extern void cLCS(int *, int *, double *, int *);

/* .Call calls */
extern SEXP checktriangleineq(SEXP, SEXP, SEXP);
extern SEXP cstringdistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cstringrefseqdistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dist2matrix(SEXP, SEXP);
extern SEXP getTraMineRDebugLevel();
extern SEXP setTraMineRDebugLevel(SEXP);
extern SEXP tmrChisq(SEXP, SEXP, SEXP);
extern SEXP tmreventinseq(SEXP, SEXP);
extern SEXP tmrinertiacontribext(SEXP, SEXP, SEXP);
extern SEXP tmrseqedist(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrsequence(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrsequencecontainevent(SEXP, SEXP, SEXP);
extern SEXP tmrsequencegetdictionary(SEXP);
extern SEXP tmrsequencegetid(SEXP);
extern SEXP tmrsequencegetlength(SEXP);
extern SEXP tmrsequencegetweight(SEXP);
extern SEXP tmrsequencesetlength(SEXP, SEXP);
extern SEXP tmrsequencesetweight(SEXP, SEXP);
extern SEXP tmrsequenceseveral(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrsequencestring(SEXP);
extern SEXP tmrsubmatrixinertia(SEXP, SEXP);
extern SEXP tmrWeightedInertiaContrib(SEXP, SEXP, SEXP);
extern SEXP tmrWeightedInertiaContribExt(SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrWeightedInertiaDist(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrWeightedInterInertia(SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrfindsubsequences(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrmatrixsubseqinseq(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP tmrseqetotse(SEXP);

static const R_CMethodDef CEntries[] = {
  {"cLCP", (DL_FUNC) &cLCP, 4},
  {"cLCS", (DL_FUNC) &cLCS, 4},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"checktriangleineq",            (DL_FUNC) &checktriangleineq,             3},
  {"cstringdistance",              (DL_FUNC) &cstringdistance,               8},
  {"cstringrefseqdistance",        (DL_FUNC) &cstringrefseqdistance,         7},
  {"dist2matrix",                  (DL_FUNC) &dist2matrix,                   2},
  {"getTraMineRDebugLevel",        (DL_FUNC) &getTraMineRDebugLevel,         0},
  {"setTraMineRDebugLevel",        (DL_FUNC) &setTraMineRDebugLevel,         1},
  {"tmrChisq",                     (DL_FUNC) &tmrChisq,                      3},
  {"tmreventinseq",                (DL_FUNC) &tmreventinseq,                 2},
  {"tmrinertiacontribext",         (DL_FUNC) &tmrinertiacontribext,          3},
  {"tmrseqedist",                  (DL_FUNC) &tmrseqedist,                   5},
  {"tmrsequence",                  (DL_FUNC) &tmrsequence,                   5},
  {"tmrsequencecontainevent",      (DL_FUNC) &tmrsequencecontainevent,       3},
  {"tmrsequencegetdictionary",     (DL_FUNC) &tmrsequencegetdictionary,      1},
  {"tmrsequencegetid",             (DL_FUNC) &tmrsequencegetid,              1},
  {"tmrsequencegetlength",         (DL_FUNC) &tmrsequencegetlength,          1},
  {"tmrsequencegetweight",         (DL_FUNC) &tmrsequencegetweight,          1},
  {"tmrsequencesetlength",         (DL_FUNC) &tmrsequencesetlength,          2},
  {"tmrsequencesetweight",         (DL_FUNC) &tmrsequencesetweight,          2},
  {"tmrsequenceseveral",           (DL_FUNC) &tmrsequenceseveral,            6},
  {"tmrsequencestring",            (DL_FUNC) &tmrsequencestring,             1},
  {"tmrsubmatrixinertia",          (DL_FUNC) &tmrsubmatrixinertia,           2},
  {"tmrWeightedInertiaContrib",    (DL_FUNC) &tmrWeightedInertiaContrib,     3},
  {"tmrWeightedInertiaContribExt", (DL_FUNC) &tmrWeightedInertiaContribExt,  4},
  {"tmrWeightedInertiaDist",       (DL_FUNC) &tmrWeightedInertiaDist,        6},
  {"tmrWeightedInterInertia",      (DL_FUNC) &tmrWeightedInterInertia,       4},
  {"tmrfindsubsequences",          (DL_FUNC) &tmrfindsubsequences,          10},
  {"tmrmatrixsubseqinseq",         (DL_FUNC) &tmrmatrixsubseqinseq,          8},
  {"tmrseqetotse",                 (DL_FUNC) &tmrseqetotse,                  1},
  {NULL, NULL, 0}
};

void R_init_TraMineR(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
