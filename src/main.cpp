#include "tonberry.h"
#include "setup.h"
#include "main_loop.h"
#include "flux.h"

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <iostream>
#include <fstream>
#include <string>

namespace mpi = boost::mpi;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: tonberry1d <input file>" << std::endl;
        return -1;
    }

    mpi::environment env;
    mpi::communicator world;
    std::cout << "I am process " << world.rank() << " of " << world.size()
                << "." << std::endl;

    Problem problem = readProblemDetails(argv[1], world.rank());

    Mesh mesh = setup(problem, world.rank(), world.size());

//     for (int i=0; i<mesh.ncellsPlusGhosts; i++) {
//         std::cout << world.rank() << " " << i-2 << " " << mesh.rho[i] << std::endl;
//     }
//
//     return -2;

    // Main loop
    mainLoop(
        mesh,
        problem.tend, problem.dtmax, mesh.dx, problem.gamma, problem.cfl
    );

    // Output
    std::string fname = "sod_c_" + std::to_string(world.rank()) + ".dat";
    std::ofstream outFile;
    outFile.open(fname);
    for (int i=2; i<mesh.ncells+2; i++) {
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
