#include "tonberry.h"
#include "flux.h"

#include <math.h>

Flux getFluxHLLC(
    double uL, double rhoL, double pL,
    double uR, double rhoR, double pR,
    double gamma) {

    // Pressure estimate from PVRS solver
    double aL = sqrt((gamma*pL)/rhoL);
    double aR = sqrt((gamma*pR)/rhoR);

    double rho_bar = 0.5*(rhoL + rhoR);
    double a_bar = 0.5*(aL + aR);

    double p_guess = 0.5*(pL + pR) - 0.5*(uR - uL)*(rho_bar*a_bar);


    double qL;
    if (p_guess <= pL) {
        qL = 1.0;
    } else {
        qL = sqrt(1.0 + (gamma + 1.0)/(2.0*gamma)*(p_guess/pL - 1.0));
    }

    double qR;
    if (p_guess <= pR) {
        qR = 1.0;
    } else {
        qR = sqrt(1.0 + (gamma + 1.0)/(2.0*gamma)*(p_guess/pR - 1.0));
    }

    // Estimate wave speeds
    double SL = uL - aL*qL;
    double SR = uR + aR*qR;
    double Sstar = pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR);
    Sstar /= (rhoL*(SL - uL) - rhoR*(SR - uR));

    // Get energy and momenta on boundaries
    double eL = pL/((gamma - 1.0)*rhoL);
    double EL = rhoL*(0.5*uL*uL + eL);
    double momL = uL*rhoL;

    double eR = pR/((gamma - 1.0)*rhoR);
    double ER = rhoR*(0.5*uR*uR + eR);
    double momR = uR*rhoR;

    // Funcion return
    Flux flux;

    // Get flux
    if (0 <= SL) {
        // Flux = F(UL)
        flux.rho = rhoL*uL;
        flux.mom = rhoL*uL*uL + pL;
        flux.E = uL*(EL + pL);
    } else if (SL < 0 && 0 <= Sstar) {
        // Flux = F(*L)
        double pLR = 0.5*(pL + pR + rhoL*(SL - uL)*(Sstar - uL) + rhoR*(SR - uR)*(Sstar - uR));
        double d = SL - Sstar;

        // Set F(UL)
        double fL_rho = rhoL*uL;
        double fL_mom = rhoL*uL*uL + pL;
        double fL_E = uL*(EL + pL);

        // Find F(*L)
        flux.rho = (Sstar*(SL*rhoL - fL_rho))/d;
        flux.mom = (Sstar*(SL*momL - fL_mom) + SL*pLR)/d;
        flux.E = (Sstar*(SL*EL - fL_E) + SL*pLR*Sstar)/d;
    } else if (Sstar < 0 && 0 <= SR) {
        // Flux = F(*R)
        double pLR = 0.5*(pL + pR + rhoL*(SL - uL)*(Sstar - uL) + rhoR*(SR - uR)*(Sstar - uR));
        double d = SR - Sstar;

        // Set F(UR)
        double fR_rho = rhoR*uR;
        double fR_mom = rhoR*uR*uR + pR;
        double fR_E = uR*(ER + pR);

        // Find F(*R)
        flux.rho = (Sstar*(SR*rhoR - fR_rho))/d;
        flux.mom = (Sstar*(SR*momR - fR_mom) + SR*pLR)/d;
        flux.E = (Sstar*(SR*ER - fR_E) + SR*pLR*Sstar)/d;
    } else {
        // Flux = F(UR)
        flux.rho = rhoR*uR;
        flux.mom = rhoR*uR*uR + pR;
        flux.E = uR*(ER + pR);
    }

    return flux;
}
