#include "OrbitConst.hh"
#include "Bunch.hh"
#include "SyncPart.hh"

#include "envelopesolver.hh"

#include <complex>

namespace envelope_solver
{
  void BasicSolver(Bunch* EnvBunch, double EmitX, double EmitY, double SpaceCharge, double Dispersion, double dLen)
  {
    double a, b, temp, SigE, Ap, Bp, **hold;

    //[part. index][x, xp, y, yp, z, zp]
    hold = EnvBunch->coordArr();

    a = std::abs(hold[0][0]);
    b = std::abs(hold[0][2]);
    SigE = hold[0][5];

    temp = sqrt( a*a + 4 * Dispersion*Dispersion * SigE*SigE) + b;

    Ap = ( ( EmitX*EmitX ) / ( a*a*a ) + 2*SpaceCharge / temp ) * dLen;
    Bp = ( ( EmitY*EmitY) / ( b*b*b ) + 2*SpaceCharge / temp ) * dLen;

    hold[0][1] += Ap;
    hold[0][3] += Bp;
  }


  void LinacTiltSolver(Bunch* EnvBunch, double SpaceCharge, double dLen)
  {
    double a,b,e,f, **hold;
    double phi, cosP, sinP;

    double C_x, C_y;

    hold = EnvBunch->coordArr();

    a = hold[0][0];
    e = hold[0][2];

    b = hold[1][0];
    f = hold[1][2];

    if(a*a+b*b-e*e-f*f == 0){
      phi = 3.14 / 4;
    }else{
      phi = -0.5 * atan( 2*(a*e+b*f) / (a*a+b*b-e*e-f*f) );
    }


    cosP = cos(phi);
    sinP = sin(phi);

    C_x = sqrt(  pow((a*f-b*e),2) / ( (e*e+f*f)*pow(cosP,2) + (a*a+b*b)*sinP*sinP +  2*(a*e+b*f)*cosP*sinP )  );
    C_y = sqrt(  pow((a*f-b*e),2) / ( (a*a+b*b)*pow(cosP,2) + (e*e+f*f)*sinP*sinP -  2*(a*e+b*f)*cosP*sinP )  );

    //These kick values are not strong/large enough

    // Particle 0: [ a, ap, e, ep]
    hold[0][1] += (  (2*SpaceCharge / (C_x + C_y)) * ( (a*pow(cosP,2)-e*sinP*cosP)/C_x + (a*pow(sinP,2) + e*sinP*cosP)/C_y )  ) * dLen;
    hold[0][3] += (  (2*SpaceCharge / (C_x + C_y)) * ( (e*pow(sinP,2)-a*sinP*cosP)/C_x + (e*pow(cosP,2) + a*sinP*cosP)/C_y )  ) * dLen;

    // Particle 1: [ b, bp, f, fp]
    hold[1][1] += (  (2*SpaceCharge / (C_x + C_y)) * ( (b*pow(cosP,2)-f*sinP*cosP)/C_x + (b*pow(sinP,2) + f*sinP*cosP)/C_y )  ) * dLen;
    hold[1][3] += (  (2*SpaceCharge / (C_x + C_y)) * ( (f*pow(sinP,2)-b*sinP*cosP)/C_x + (f*pow(cosP,2) + b*sinP*cosP)/C_y )  ) * dLen;

  }

}//end of namespace
