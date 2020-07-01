#include "tonberry.h"
#include "setup.h"
#include "main_loop.h"
#include "flux.h"

#include <iostream>
#include <fstream>

int main() {
    // Hard code some parameters for initial testing
    int ncells = 100;
    double dtmax = 0.1;
    double tend = 0.25;
    double gamma = 1.4;
    double length = 1.0;
    double x0 = 0.5;
    double cfl = 0.6;

    Mesh mesh = setup(length, x0, ncells, gamma);

    // Main loop
    mainLoop(
        mesh.rho, mesh.p, mesh.u, mesh.E, mesh.mom,
        tend, dtmax, mesh.dx, gamma, cfl
    );


    // Output
    std::ofstream outFile;
    outFile.open("sod_c.dat");
    for (int i=0; i<ncells+4; i++) {
        double ein = mesh.E[i]/mesh.rho[i] - 0.5*mesh.u[i]*mesh.u[i];
        outFile << i << " ";
        outFile << mesh.x[i] << " ";
        outFile << mesh.rho[i] << " ";
        outFile << mesh.p[i] << " ";
        outFile << mesh.u[i] << " ";
        outFile << ein;
        outFile << std::endl;
    }
    outFile.close();

    return 0;
}
