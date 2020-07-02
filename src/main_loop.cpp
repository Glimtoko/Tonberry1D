#include "flux.h"
#include "main_loop.h"
#include "comms.h"

#include <iostream>
#include <math.h>

#include <boost/mpi.hpp>

namespace mpi = boost::mpi;

void mainLoop(
    Mesh &mesh,
    double tend, double dtmax, double dx, double gamma, double cfl
){
    // MPI
    mpi::communicator world;

    // Number of non-ghost cells in mesh
    int nCellsGhosts = mesh.ncellsPlusGhosts;

    // Boundary indices
    int gL1 = 1;
    int L = 2;
    int gR2 = nCellsGhosts - 1;
    int gR1 = nCellsGhosts - 2;

    // Storage for fluxes
    std::vector<Flux> flux;
    Flux zeroFlux;

    // Left and Right boundary reconstructions
    TB_ARRAY rhoL;
    TB_ARRAY rhoR;
    TB_ARRAY momL;
    TB_ARRAY momR;
    TB_ARRAY EL;
    TB_ARRAY ER;

    // Slopes
    TB_ARRAY deltaiRho;
    TB_ARRAY deltaimom;
    TB_ARRAY deltaiE;


    // Allocate storage
    flux.assign(nCellsGhosts, zeroFlux);
    rhoL.assign(nCellsGhosts, 0.0);
    rhoR.assign(nCellsGhosts, 0.0);
    momL.assign(nCellsGhosts, 0.0);
    momR.assign(nCellsGhosts, 0.0);
    EL.assign(nCellsGhosts, 0.0);
    ER.assign(nCellsGhosts, 0.0);
    deltaiRho.assign(nCellsGhosts, 0.0);
    deltaimom.assign(nCellsGhosts, 0.0);
    deltaiE.assign(nCellsGhosts, 0.0);

    // Omega parameter sets bias of slope calculation
    double omega = 0.0;

    // Main loop
    double t = 0.0;
    int step = 1;
    double dt = 0.0;
    std::cout << "Stating main loop" << std::endl;
    for( ; ; ) {
        // Parallel update
        parallelUpdate(mesh);

        // Get estimate for maximum wave speed
        double S = 0;
        for (int i=0; i<nCellsGhosts; i++) {
            double u = mesh.mom[i]/mesh.rho[i];
            double p = (mesh.gamma - 1.0)*(mesh.E[i] - 0.5*mesh.rho[i]*u*u);
            double a = sqrt((gamma*p)/mesh.rho[i]);
            S = std::max(S, a + u);
        }

        // Calculate timestep
        dt = std::min(dtmax, cfl*dx/S);
        dt = mpi::all_reduce(world, dt, mpi::minimum<double>());

        if (world.rank() == 0) {
            std::cout << "Step: " << step;
            std::cout << ", time = " << t;
            std::cout << ", dt = " << dt << std::endl;
        }

        // Calculate slopes
        for (int i=gL1; i<gR2; i++) {
            deltaiRho[i] = getSlope(mesh.rho, i, omega);
            deltaimom[i] = getSlope(mesh.mom, i, omega);
            deltaiE[i] = getSlope(mesh.E, i, omega);
        }

        // Data reconstruction
        for (int i=gL1; i<gR2; i++) {
            rhoL[i] = mesh.rho[i] - 0.5*deltaiRho[i];
            rhoR[i] = mesh.rho[i] + 0.5*deltaiRho[i];

            momL[i] = mesh.mom[i] - 0.5*deltaimom[i];
            momR[i] = mesh.mom[i] + 0.5*deltaimom[i];

            EL[i] = mesh.E[i] - 0.5*deltaiE[i];
            ER[i] = mesh.E[i] + 0.5*deltaiE[i];
        }

        // Data evolution to half timestep
        double f = 0.5*dt/dx;
        for (int i=gL1; i<gR2; i++) {
            // Reconstruct u and p on boundaries using conserved quantities
            double uL = momL[i]/rhoL[i];
            double uR = momR[i]/rhoR[i];

            double pL = (gamma - 1.0)*(EL[i] - 0.5*rhoL[i]*uL*uL);
            double pR = (gamma - 1.0)*(ER[i] - 0.5*rhoR[i]*uR*uR);

            // Calculate flux difference across cell
            double dF_rho = f*(calc_flux_rho(rhoL[i], uL) -
                               calc_flux_rho(rhoR[i], uR));

            double dF_mom = f*(calc_flux_mom(rhoL[i], uL, pL) -
                               calc_flux_mom(rhoR[i], uR, pR));

            double dF_E = f*(calc_flux_E(uL, EL[i], pL) -
                             calc_flux_E(uR, ER[i], pR));

            // Use flux differences to update data
            rhoL[i] += dF_rho;
            rhoR[i] += dF_rho;

            momL[i] += dF_mom;
            momR[i] += dF_mom;

            EL[i] += dF_E;
            ER[i] += dF_E;
        }

        // Get intercell fluxes from Riemann solver
        for (int i=gL1; i<gR1; i++) {
            // N.b. use of uL_cell and uR_cell in p calculations is correct
            double rhoL_cell = rhoR[i];
            double uL_cell = momR[i]/rhoR[i];
            double pL_cell = (gamma - 1.0)*(ER[i] - 0.5*rhoR[i]*uL_cell*uL_cell);

            double rhoR_cell = rhoL[i+1];
            double uR_cell = momL[i+1]/rhoL[i+1];
            double pR_cell = (gamma - 1.0)*(EL[i+1] - 0.5*rhoL[i+1]*uR_cell*uR_cell);

            flux[i] = getFluxHLLC(
                uL_cell, rhoL_cell, pL_cell,
                uR_cell, rhoR_cell, pR_cell,
                gamma);

        }

        // Update conserved and primitive quantities
        f = dt/dx;
        for (int i=L; i<gR1; i++) {
            // Conservative update
            mesh.rho[i] += f*(flux[i-1].rho - flux[i].rho);
            mesh.mom[i] += f*(flux[i-1].mom - flux[i].mom);
            mesh.E[i] += f*(flux[i-1].E - flux[i].E);

            // Test for negative density and pressure
            if (mesh.rho[i] < 0.0) {
                std::cout << "Negative density in cell " << i << std::endl;
                std::cout << "Density = " << mesh.rho[i] << std::endl;
                std::cout << "Energy = " << mesh.E[i] << std::endl;
                std::cout << "Momentum = " << mesh.mom[i] << std::endl;
            }
        }

        if (t > tend) {
            break;
        }

        t += dt;
        step += 1;
    }
}


double getSlope(TB_ARRAY U, int i, double omega) {
    // Calculate slope of vector U in cell i
    double di1 = U[i] - U[i-1];
    double di2 = U[i+1] - U[i];
    double diU = 0.5*(1.0 + omega)*di1 + 0.5*(1.0 - omega)*di2;

    // Slope limiter - Van Leer
    if (di2 == 0) {
        diU = 0.0;
    } else {
        double r = di1/di2;
        if (r <= 0.0) {
            diU = 0.0;
        } else {
            double xiR = 2.0/(1.0 - omega + (1 + omega)*r);
            double xi = std::min(2*r/(1+r), xiR);
            diU = diU * xi;
        }
    }

    return diU;
}


double calc_flux_rho(double rho, double u) {
    return rho*u;
}


double calc_flux_mom(double rho, double u, double p) {
    return rho*u*u + p;
}


double calc_flux_E(double u, double E, double p) {
    return u*(E + p);
}
