#include "tonberry.h"
#include "setup.h"
#include "comms.h"

#include <boost/program_options.hpp>

#include <iostream>
#include <fstream>

namespace bpo = boost::program_options;

Problem readProblemDetails(char *fname, int rank) {
    bpo::options_description optionList;

    optionList.add_options()
        ("problem.ncells",bpo::value<int>()->default_value(100), "Number of cells")
        ("problem.tend",bpo::value<double>()->default_value(2.5), "End time")
        ("problem.dtmax",bpo::value<double>()->default_value(0.1), "Max timestep")
        ("problem.gamma",bpo::value<double>()->default_value(1.0), "Gamma (EoS)")
        ("problem.cfl",bpo::value<double>()->default_value(0.6), "CFL")

        ("domain.length", bpo::value<double>()->default_value(1.0), "Length of domain")
        ("domain.x0", bpo::value<double>()->default_value(0.5), "Position of interface")

        ("domain.rhoL", bpo::value<double>()->default_value(1.0), "Left-side density")
        ("domain.pL", bpo::value<double>()->default_value(1.0), "Left-side pressure")
        ("domain.uL", bpo::value<double>()->default_value(0.0), "Left-side velocity")

        ("domain.rhoR", bpo::value<double>()->default_value(0.1), "Right-side density")
        ("domain.pR", bpo::value<double>()->default_value(0.1), "Right-side pressure")
        ("domain.uR", bpo::value<double>()->default_value(0.0), "Right-side velocity")
        ;

    // Open an input file
    std::ifstream inFile;
    inFile.open(fname);

    // Parse file
    bpo::variables_map vm;
    bpo::store(bpo::parse_config_file(inFile, optionList), vm);

    inFile.close();

    // Create Problem struct
    Problem problem;

    problem.ncells = vm["problem.ncells"].as<int>();
    problem.tend = vm["problem.tend"].as<double>();
    problem.gamma = vm["problem.gamma"].as<double>();
    problem.dtmax = vm["problem.dtmax"].as<double>();
    problem.cfl = vm["problem.cfl"].as<double>();

    problem.length = vm["domain.length"].as<double>();
    problem.x0 = vm["domain.x0"].as<double>();

    problem.rhoL = vm["domain.rhoL"].as<double>();
    problem.pL = vm["domain.pL"].as<double>();
    problem.uL = vm["domain.uL"].as<double>();

    problem.rhoR = vm["domain.rhoR"].as<double>();
    problem.pR = vm["domain.pR"].as<double>();
    problem.uR = vm["domain.uR"].as<double>();

    if (rank == 0) {
        std::cout << "Problem input:" << std::endl;
        std::cout << "  Length: " << problem.length << std::endl
                << "  x0: " << problem.x0 << std::endl
                << "  ncells: " << problem.ncells << std::endl << std::endl

                << "  rhoL: " << problem.rhoL << std::endl
                << "  pL:" << problem.pL << std::endl
                << "  uL:" << problem.uL << std::endl << std::endl

                << "  rhoR: " << problem.rhoR << std::endl
                << "  pR:" << problem.pR << std::endl
                << "  uR:" << problem.uR << std::endl << std::endl

                << "  tend:" << problem.tend << std::endl
                << "  dtmax:" << problem.dtmax << std::endl
                << "  gamma:" << problem.gamma << std::endl
                << "  CFL:" << problem.cfl << std::endl
                ;
    }

    return problem;
}

Mesh setup(Problem problem, int rank, int size) {
    // Mesh object
    Mesh mesh;

    int ncellsPerProc = problem.ncells / size;

    if (size * ncellsPerProc != problem.ncells) {
        if (rank == 0) {
            std::cout << "ERROR: Number of cells does not divide ";
            std::cout << "by number of processors!";
            std::cout << std::endl;
        }
    }

    mesh.globalNCells = problem.ncells;
    mesh.ncells = ncellsPerProc;
    mesh.ncellsPlusGhosts = mesh.ncells + 4;
    mesh.dx = problem.length/problem.ncells;

    // Need to store gamma in mesh for post-comms reconstruction
    mesh.gamma = problem.gamma;

    // Set indices
    int end = mesh.ncellsPlusGhosts - 1;
    int gL2 = 0;
    int gL1 = 1;
    int L = 2;
    int gR2 = end;
    int gR1 = end - 1;
    int R = end - 2;

    // Allocate storage - allow for ghosts
    mesh.x.assign(problem.ncells+4, 0.0);
    mesh.rho.assign(problem.ncells+4, 0.0);
    mesh.mom.assign(problem.ncells+4, 0.0);
    mesh.E.assign(problem.ncells+4, 0.0);

    // Global index of left-hand  and right-hand cells
    mesh.globalL = ncellsPerProc*rank;
    mesh.globalR = ncellsPerProc*(rank+1) - 1;

    // Set x - cell centre positions
    for (int i=L; i<gR1; i++) {
        int gi = i - 1 + mesh.globalL;
        mesh.x[i] = (gi - 0.5)*mesh.dx;
    }
    mesh.x[gL1] = mesh.x[L] - mesh.dx;
    mesh.x[gL2] = mesh.x[L] - 2*mesh.dx;
    mesh.x[gR1] = mesh.x[R] + mesh.dx;
    mesh.x[gR2] = mesh.x[R] + 2*mesh.dx;


    // Set initial fields
    for (int i=L; i<gR1; i++) {
        double xupper = mesh.x[i] + mesh.dx/2.0;
        if (xupper <= problem.x0) {
            mesh.rho[i] = problem.rhoL;
            mesh.mom[i] = problem.rhoL*problem.uL;

            double e = problem.pL/((problem.gamma - 1.0)*mesh.rho[i]);
            mesh.E[i] = mesh.rho[i]*(0.5*problem.uL*problem.uL + e);
        } else {
            mesh.rho[i] = problem.rhoR;
            mesh.mom[i] = problem.rhoR*problem.uR;

            double e = problem.pR/((problem.gamma - 1.0)*mesh.rho[i]);
            mesh.E[i] = mesh.rho[i]*(0.5*problem.uR*problem.uR + e);
        }
    }

    // Set boundaries
    setBoundaries(mesh);

    // Parallel update
    parallelUpdate(mesh);

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
    if (mesh.globalL == 0) {
        mesh.rho[gL1] = mesh.rho[L];
        mesh.E[gL1] = mesh.E[L];
        mesh.mom[gL1] = mesh.mom[L];

        mesh.rho[gL2] = mesh.rho[L+1];
        mesh.E[gL2] = mesh.E[L+1];
        mesh.mom[gL2] = mesh.mom[L+1];
    }

    if (mesh.globalR == mesh.globalNCells - 1) {
        mesh.rho[gR1] = mesh.rho[R];
        mesh.E[gR1] = mesh.E[R];
        mesh.mom[gR1] = mesh.mom[R];

        mesh.rho[gR2] = mesh.rho[R-1];
        mesh.E[gR2] = mesh.E[R-1];
        mesh.mom[gR2] = mesh.mom[R-1];
    }
}
