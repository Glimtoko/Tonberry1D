#include "tonberry.h"
#include "setup.h"

Mesh setup(double length, double x0, int ncells, double gamma) {
    // Mesh object
    Mesh mesh;

    // Initial quick hack - only support Sod
    double uL = 0.0;
    double uR = 0.0;
    double rhoL = 1.0;
    double rhoR = 0.125;
    double pL = 1.0;
    double pR = 0.1;

    mesh.ncells = ncells;
    mesh.ncellsPlusGhosts = ncells + 4;
    mesh.dx = length/ncells;

    // Set indices
    int end = mesh.ncellsPlusGhosts - 1;
    int gL2 = 0;
    int gL1 = 1;
    int L = 2;
    int gR2 = end;
    int gR1 = end - 1;
    int R = end - 2;

    // Allocate storage - allow for ghosts
    mesh.x.assign(ncells+4, 0.0);
    mesh.rho.assign(ncells+4, 0.0);
    mesh.mom.assign(ncells+4, 0.0);
    mesh.E.assign(ncells+4, 0.0);
    mesh.p.assign(ncells+4, 0.0);
    mesh.u.assign(ncells+4, 0.0);

    // Set x - cell centre positions
    for (int i=1; i<ncells+2; i++) {
        mesh.x[i] = (i - 0.5)*mesh.dx;
    }
    mesh.x[gL1] = -mesh.dx;
    mesh.x[gL2] = -1.5*mesh.dx;
    mesh.x[gR1] = mesh.x[R] + mesh.dx;
    mesh.x[gR2] = mesh.x[R] + 1.5*mesh.dx;


    // Set initial fields
    for (int i=L; i<gR1; i++) {
        double xupper = mesh.x[i] + mesh.dx/2.0;
        if (xupper <= x0) {
            mesh.rho[i] = rhoL;
            mesh.mom[i] = rhoL*uL;
            mesh.p[i] = pL;
            mesh.u[i] = uL;
        } else {
            mesh.rho[i] = rhoR;
            mesh.mom[i] = rhoR*uR;
            mesh.p[i] = pR;
            mesh.u[i] = uR;
        }
        double e = mesh.p[i]/((gamma - 1.0)*mesh.rho[i]);
        mesh.E[i] = mesh.rho[i]*(0.5*mesh.u[i]*mesh.u[i] + e);
    }

    // Set boundaries
    setBoundaries(mesh);

    return mesh;
}


void setBoundaries(Mesh &mesh) {
    // Set indices
    int end = mesh.ncells + 4 - 1;
    int gL2 = 0;
    int gL1 = 1;
    int L = 2;
    int gR2 = end;
    int gR1 = end - 1;
    int R = end - 2;

    // Set boundary values
    mesh.rho[gL1] = mesh.rho[L];
    mesh.E[gL1] = mesh.E[L];
    mesh.mom[gL1] = mesh.mom[L];
    mesh.p[gL1] = mesh.p[L];
    mesh.u[gL1] = mesh.u[L];

    mesh.rho[gL2] = mesh.rho[L+1];
    mesh.E[gL2] = mesh.E[L+1];
    mesh.mom[gL2] = mesh.mom[L+1];
    mesh.p[gL2] = mesh.p[L+1];
    mesh.u[gL2] = mesh.u[L+1];

    mesh.rho[gR1] = mesh.rho[R];
    mesh.E[gR1] = mesh.E[R];
    mesh.mom[gR1] = mesh.mom[R];
    mesh.p[gR1] = mesh.p[R];
    mesh.u[gR1] = mesh.u[R];

    mesh.rho[gR2] = mesh.rho[R-1];
    mesh.E[gR2] = mesh.E[R-1];
    mesh.mom[gR2] = mesh.mom[R-1];
    mesh.p[gR2] = mesh.p[R-1];
    mesh.u[gR2] = mesh.u[R-1];
}
