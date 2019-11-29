#ifndef ENVSOLVER_H
#define ENVSOLVER_H

#include "Bunch.hh"

namespace envelope_solver
{
    void BasicSolver(Bunch* EnvBunch, double EmitX, double EmitY, double SpaceCharge, double Dispersion, double dLen);

    void LinacTiltSolver(Bunch* EnvBunch, double SpaceCharge, double dLen);
}

#endif // ENVSOLVER_H
