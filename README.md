# halobench

## Simple MPI benchmark of various methods for doing 3D halo swapping

The code sets up a 3D array of processes and swaps halo data up and
down in all three directions. It uses a number of different methods -
the aim is to see which is fastest:

*    pairwise swapping in turn using MPI_Sendrecv
*    pairwise swapping in turn but breaking the deadlock by splitting the processes into odd and even (here called red and black) sets: red processes send first and receive second, black processes do the opposite; this approach requires an even number of processes in each dimension
*    pairwise swapping using Isend / Recv / Wait
*    pairwise swapping using Irecv / Send / Wait
*    pairwise swapping using Irecv / Isend / Waitall (x2)
*    non-blocking comms in all six directions at once using Irecv / Isend / Waitall (x12)
*    persistent communications
*    neighbourhood collectives

Note that the neighbourhood collectives may fail the consistency test
on small numbers of processes due to an unfortunate bug in the current
implementation!