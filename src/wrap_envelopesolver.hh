#ifndef WRAP_ENVSOLVER_H
#define WRAP_ENVSOLVER_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_envelopesolver{
    void initEnvSolver(void);
    PyObject* getEnvSolverType(const char* name);
  }

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_ENVSOLVER
