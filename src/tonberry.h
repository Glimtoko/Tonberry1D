#ifndef TONBERRYH
#define TONBERRYH

#include <vector>

typedef std::vector<double> TB_ARRAY;

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
#endif
