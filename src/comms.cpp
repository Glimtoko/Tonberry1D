#include "tonberry.h"

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;
void parallelUpdate(Mesh mesh){
    mpi::environment env;
    mpi::communicator world;


    int rank = world.rank();
    int size = world.size();

    int nL = rank - 1;
    int nR = rank + 1;

    int iL = 0;
    int iR = mesh.ncellsPlusGhosts - 3;



    // Pre-post receives
    if (rank == 0) {
        // Packages
        TB_ARRAY rhoRecvRPack(2, 0.0);
        TB_ARRAY pRecvRPack(2, 0.0);
        TB_ARRAY uRecvRPack(2, 0.0);

        mpi::request rReqs[3];
        rReqs[0] = world.irecv(nR, rank*100+1, rhoRecvRPack);
        rReqs[1] = world.irecv(nR, rank*100+2, pRecvRPack);
        rReqs[2] = world.irecv(nR, rank*100+3, uRecvRPack);
    } else if (rank == size - 1) {
        // Packages
        TB_ARRAY rhoRecvLPack(2, 0.0);
        TB_ARRAY pRecvLPack(2, 0.0);
        TB_ARRAY uRecvLPack(2, 0.0);

        mpi::request rReqs[3];
        rReqs[0] = world.irecv(nL, rank*100+11, rhoRecvLPack);
        rReqs[1] = world.irecv(nL, rank*100+12, pRecvLPack);
        rReqs[2] = world.irecv(nL, rank*100+13, uRecvLPack);
    } else {
        // Packages
        TB_ARRAY rhoRecvRPack(2, 0.0);
        TB_ARRAY pRecvRPack(2, 0.0);
        TB_ARRAY uRecvRPack(2, 0.0);
        TB_ARRAY rhoRecvLPack(2, 0.0);
        TB_ARRAY pRecvLPack(2, 0.0);
        TB_ARRAY uRecvLPack(2, 0.0);

        mpi::request rReqs[6];
        rReqs[0] = world.irecv(nR, rank*100+1, rhoRecvRPack);
        rReqs[1] = world.irecv(nR, rank*100+2, pRecvRPack);
        rReqs[2] = world.irecv(nR, rank*100+3, uRecvRPack);
        rReqs[3] = world.irecv(nL, rank*100+11, rhoRecvLPack);
        rReqs[4] = world.irecv(nL, rank*100+12, pRecvLPack);
        rReqs[5] = world.irecv(nL, rank*100+13, uRecvLPack);
    }

    // Post sends
    if (rank == 0) {
        // Packages
        TB_ARRAY rhoSendRPack(2, 0.0);
        TB_ARRAY pSendRPack(2, 0.0);
        TB_ARRAY uSendRPack(2, 0.0);

        rhoSendRPack[0] = mesh.rho[iL];
        rhoSendRPack[1] = mesh.rho[iL + 1];
        pSendRPack[0] = mesh.p[iL];
        pSendRPack[1] = mesh.p[iL + 1];
        uSendRPack[0] = mesh.u[iL];
        uSendRPack[1] = mesh.u[iL + 1];

        mpi::request sReqs[3];
        sReqs[0] = world.isend(nR, nR*100+11, rhoSendRPack);
        sReqs[1] = world.isend(nR, nR*100+12, pSendRPack);
        sReqs[2] = world.isend(nR, nR*100+13, uSendRPack);


    }
}
