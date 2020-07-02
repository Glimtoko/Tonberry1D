#include "tonberry.h"

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

namespace mpi = boost::mpi;
void parallelUpdate(Mesh &mesh){
    mpi::environment env;
    mpi::communicator world;


    int rank = world.rank();
    int size = world.size();

    if (size == 1) return;

    int nL = rank - 1;
    int nR = rank + 1;

    int Lg1 = 0;
    int Lg2 = 1;
    int L = 2;
    int R = mesh.ncellsPlusGhosts - 3;
    int Rg1 = mesh.ncellsPlusGhosts - 2;
    int Rg2 = mesh.ncellsPlusGhosts - 1;


    // Receive buffers
    TB_ARRAY rhoRecvRPack(2, 0.0);
    TB_ARRAY momRecvRPack(2, 0.0);
    TB_ARRAY ERecvRPack(2, 0.0);
    TB_ARRAY rhoRecvLPack(2, 0.0);
    TB_ARRAY momRecvLPack(2, 0.0);
    TB_ARRAY ERecvLPack(2, 0.0);


    // Pre-post receives
    int nreqs;
    mpi::request rReqs[6];
    if (rank == 0) {
        nreqs =  3;
        rReqs[0] = world.irecv(nR, rank*100+1, rhoRecvRPack);
        rReqs[1] = world.irecv(nR, rank*100+2, momRecvRPack);
        rReqs[2] = world.irecv(nR, rank*100+3, ERecvRPack);
    } else if (rank == size - 1) {
        nreqs =  3;
        rReqs[0] = world.irecv(nL, rank*100+11, rhoRecvLPack);
        rReqs[1] = world.irecv(nL, rank*100+12, momRecvLPack);
        rReqs[2] = world.irecv(nL, rank*100+13, ERecvLPack);
    } else {
        nreqs =  6;
        rReqs[0] = world.irecv(nR, rank*100+1, rhoRecvRPack);
        rReqs[1] = world.irecv(nR, rank*100+2, momRecvRPack);
        rReqs[2] = world.irecv(nR, rank*100+3, ERecvRPack);
        rReqs[3] = world.irecv(nL, rank*100+11, rhoRecvLPack);
        rReqs[4] = world.irecv(nL, rank*100+12, momRecvLPack);
        rReqs[5] = world.irecv(nL, rank*100+13, ERecvLPack);
    }

    // Post sends
    if (rank == 0) {
        // Packages
        TB_ARRAY rhoSendRPack(2, 0.0);
        TB_ARRAY momSendRPack(2, 0.0);
        TB_ARRAY ESendRPack(2, 0.0);

        // Pack data
        rhoSendRPack[0] = mesh.rho[R - 1];
        rhoSendRPack[1] = mesh.rho[R];
        momSendRPack[0] = mesh.mom[R - 1];
        momSendRPack[1] = mesh.mom[R];
        ESendRPack[0] = mesh.E[R - 1];
        ESendRPack[1] = mesh.E[R];

        mpi::request sReqs[3];
        sReqs[0] = world.isend(nR, nR*100+11, rhoSendRPack);
        sReqs[1] = world.isend(nR, nR*100+12, momSendRPack);
        sReqs[2] = world.isend(nR, nR*100+13, ESendRPack);
    } else if (rank == size - 1) {
        // Packages
        TB_ARRAY rhoSendLPack(2, 0.0);
        TB_ARRAY momSendLPack(2, 0.0);
        TB_ARRAY ESendLPack(2, 0.0);

        // Pack data
        rhoSendLPack[0] = mesh.rho[L];
        rhoSendLPack[1] = mesh.rho[L + 1];
        momSendLPack[0] = mesh.mom[L];
        momSendLPack[1] = mesh.mom[L + 1];
        ESendLPack[0] = mesh.E[L];
        ESendLPack[1] = mesh.E[L + 1];

        mpi::request sReqs[3];
        nreqs =  3;
        sReqs[0] = world.isend(nL, nL*100+1, rhoSendLPack);
        sReqs[1] = world.isend(nL, nL*100+2, momSendLPack);
        sReqs[2] = world.isend(nL, nL*100+3, ESendLPack);
    } else {
        // Packages
        TB_ARRAY rhoSendRPack(2, 0.0);
        TB_ARRAY momSendRPack(2, 0.0);
        TB_ARRAY ESendRPack(2, 0.0);
        TB_ARRAY rhoSendLPack(2, 0.0);
        TB_ARRAY momSendLPack(2, 0.0);
        TB_ARRAY ESendLPack(2, 0.0);

        // Pack data
        rhoSendRPack[0] = mesh.rho[R - 1];
        rhoSendRPack[1] = mesh.rho[R];
        momSendRPack[0] = mesh.mom[R - 1];
        momSendRPack[1] = mesh.mom[R];
        ESendRPack[0] = mesh.E[R - 1];
        ESendRPack[1] = mesh.E[R];

        rhoSendLPack[0] = mesh.rho[L];
        rhoSendLPack[1] = mesh.rho[L + 1];
        momSendLPack[0] = mesh.mom[L];
        momSendLPack[1] = mesh.mom[L + 1];
        ESendLPack[0] = mesh.E[L];
        ESendLPack[1] = mesh.E[L + 1];

        mpi::request sReqs[6];
        nreqs =  6;
        sReqs[0] = world.isend(nR, nR*100+11, rhoSendRPack);
        sReqs[1] = world.isend(nR, nR*100+12, momSendRPack);
        sReqs[2] = world.isend(nR, nR*100+13, ESendRPack);
        sReqs[3] = world.isend(nL, nL*100+1, rhoSendLPack);
        sReqs[4] = world.isend(nL, nL*100+2, momSendLPack);
        sReqs[5] = world.isend(nL, nL*100+3, ESendLPack);
    }

    mpi::wait_all(rReqs, rReqs+nreqs);

    // Unpack received data
    if (rank == 0) {
        mesh.rho[Rg1] = rhoRecvRPack[0];
        mesh.rho[Rg2] = rhoRecvRPack[1];

        mesh.mom[Rg1] = momRecvRPack[0];
        mesh.mom[Rg2] = momRecvRPack[1];

        mesh.E[Rg1] = ERecvRPack[0];
        mesh.E[Rg2] = ERecvRPack[1];
    } else if (rank == size - 1) {
        mesh.rho[Lg1] = rhoRecvLPack[0];
        mesh.rho[Lg2] = rhoRecvLPack[1];

        mesh.mom[Lg1] = momRecvLPack[0];
        mesh.mom[Lg2] = momRecvLPack[1];

        mesh.E[Lg1] = ERecvLPack[0];
        mesh.E[Lg2] = ERecvLPack[1];
    } else {
        mesh.rho[Rg1] = rhoRecvRPack[0];
        mesh.rho[Rg2] = rhoRecvRPack[1];

        mesh.mom[Rg1] = momRecvRPack[0];
        mesh.mom[Rg2] = momRecvRPack[1];

        mesh.E[Rg1] = ERecvRPack[0];
        mesh.E[Rg2] = ERecvRPack[1];

        mesh.rho[Lg1] = rhoRecvLPack[0];
        mesh.rho[Lg2] = rhoRecvLPack[1];

        mesh.mom[Lg1] = momRecvLPack[0];
        mesh.mom[Lg2] = momRecvLPack[1];

        mesh.E[Lg1] = ERecvLPack[0];
        mesh.E[Lg2] = ERecvLPack[1];
    }
}
