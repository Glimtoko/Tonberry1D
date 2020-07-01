#ifndef TONBERRYH
#define TONBERRYH
#pragma once

#include <vector>
#include <iostream>

typedef std::vector<double> TB_ARRAY;

struct Problem {
    double length;
    double x0;
    int ncells;
    double cfl;

    double rhoL;
    double pL;
    double uL;

    double rhoR;
    double pR;
    double uR;

    double tend;
    double dtmax;
    double gamma;
};

struct Mesh {
    TB_ARRAY x;
    TB_ARRAY rho;
    TB_ARRAY mom;
    TB_ARRAY E;
    TB_ARRAY p;
    TB_ARRAY u;

    int ncells;
    int ncellsPlusGhosts;
    double dx;
};

struct Flux {
    double rho;
    double mom;
    double E;
};

// std::ostream& operator<<(std::ostream& os, const Problem& problem) {
//     return os /*<< "Length: " << problem.length << std::endl
//               << "x0: " << problem.x0 << std::endl
//               << "ncells: " << problem.ncells << std::endl << std::endl
//               << "rhoL: " << problem.rhoL << std::endl
//               << "pL:" << problem.pL << std::endl
//               << "uL:" << problem.uL << std::endl
//               ;*/
// }

#endif
