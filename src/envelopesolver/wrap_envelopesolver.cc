#include "Python.h"
#include "orbit_mpi.hh"

#include "pyORBIT_Object.hh"

#include "envelopesolver.hh"
#include "wrap_envelopesolver.hh"

namespace wrap_envelopesolver
{
    void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C"
{
#endif

    static PyObject* wrap_BasicSolver(PyObject *self, PyObject *args)
    {
        PyObject* pyBunch;
        double EmitX, EmitY, SpaceCharge, Dispersion, dLen;

        if(!PyArg_ParseTuple(	args, "Oddddd:BasicSolver", &pyBunch, &EmitX, &EmitY, &SpaceCharge, &Dispersion, &dLen))
        {
            error("EnvelopeSolver - BasicSolver - cannot parse arguments!");
        }

        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
        envelope_solver::BasicSolver(cpp_bunch, EmitX, EmitY, SpaceCharge, Dispersion, dLen);
        Py_INCREF(Py_None);
        return Py_None;
    }

    static PyObject* wrap_LinacTiltSolver(PyObject *self, PyObject *args)
    {
        PyObject* pyBunch;
        double SpaceCharge, dLen;

        if(!PyArg_ParseTuple(	args, "Odd:LinacTiltSolver", &pyBunch, &SpaceCharge, &dLen))
        {
            error("EnvelopeSolver - LinacTiltSolver - cannot parse arguments!");
        }

        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object *) pyBunch)->cpp_obj;
        envelope_solver::LinacTiltSolver(cpp_bunch, SpaceCharge, dLen);
        Py_INCREF(Py_None);
        return Py_None;
    }


    static PyMethodDef envsolverMethods[] =
    {
      {"BasicSolver",         wrap_BasicSolver,           METH_VARARGS, "Basic Solver - no tilt, self consistance distro"},
      {"LinacTiltSolver",     wrap_LinacTiltSolver,       METH_VARARGS, "LINAC Tilting Envelope Solver - with tilting for self consistant distro"},
      { NULL, NULL }
  };

  void initEnvSolver(void)
  {
      PyObject *m, *d;
      m = Py_InitModule((char*)"envelope_solver", envsolverMethods);
      d = PyModule_GetDict(m);
  }

	PyObject* getEnvSolverType(const char* name){
		PyObject* mod = PyImport_ImportModule("envelope_solver");
		PyObject* pyType = PyObject_GetAttrString(mod,name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
	}


#ifdef __cplusplus
}
#endif

}
