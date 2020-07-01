#include "tonberry.h"

Problem readProblemDetails(char *fname, int rank);
Mesh setup(Problem problem, int rank, int size);
void setBoundaries(Mesh &mesh);
