#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void loopC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loopU(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
//extern void mstepU(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mstepU (int *NN, int *pp, int *GG, double *z,
            double *sigmar, double *invsigmar, double *mu, double *mmtol,
            int *mmmax, double *x, char **covtype, double *PXgood);
extern void RdCN(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RestepC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RestepU(void *, void *, void *, void *, void *, void *, void *, void *);
extern void RllikelihoodC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RllikelihoodU(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"loopC",         (DL_FUNC) &loopC,         25},
  {"loopU",         (DL_FUNC) &loopU,         20},
  {"mstepU",        (DL_FUNC) &mstepU,        12},
  {"RdCN",          (DL_FUNC) &RdCN,           9},
  {"RestepC",       (DL_FUNC) &RestepC,       10},
  {"RestepU",       (DL_FUNC) &RestepU,        8},
  {"RllikelihoodC", (DL_FUNC) &RllikelihoodC, 10},
  {"RllikelihoodU", (DL_FUNC) &RllikelihoodU,  8},
  {NULL, NULL, 0}
};

void R_init_ContaminatedMixt(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
