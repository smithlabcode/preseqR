/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern void c_PS2CF(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"c_PS2CF", (DL_FUNC) &c_PS2CF, 12},
    {NULL, NULL, 0}
};

void R_init_preseqR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
