#ifndef MAINLOOPH
#define MAINLOOPH
#include "tonberry.h"

void mainLoop(
    Mesh &mesh,
    double tend, double dtmax, double dx, double gamma, double cfl
);

double getSlope(TB_ARRAY U, int i, double omega);

double calc_flux_rho(double rho, double u);

double calc_flux_mom(double rho, double u, double p);

double calc_flux_E(double u, double E, double p);
#endif
