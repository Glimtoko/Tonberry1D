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

    int globalNCells;
    int ncells;
    int ncellsPlusGhosts;
    double dx;

    int globalL;
    int globalR;
};

struct Flux {
    double rho;
    double mom;
    double E;
};

#endif
