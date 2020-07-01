#include "tonberry.h"
#include "setup.h"
#include "main_loop.h"
#include "flux.h"

#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: tonberry1d <input file>" << std::endl;
        return -1;
    }

    Problem problem = readProblemDetails(argv[1]);

    Mesh mesh = setup(problem);

    // Main loop
    mainLoop(
        mesh.rho, mesh.p, mesh.u, mesh.E, mesh.mom,
        problem.tend, problem.dtmax, mesh.dx, problem.gamma, problem.cfl
    );


    // Output
    std::ofstream outFile;
    outFile.open("sod_c.dat");
    for (int i=0; i<mesh.ncells+4; i++) {
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
