module haloswap

  use haloparams
  
  implicit none

  integer, parameter :: nrequest = 2  ! receive and send

  integer, parameter :: rrequest = 1
  integer, parameter :: srequest = 2

  integer, parameter :: defsendtag = 1
  integer, parameter :: defrecvtag = 1

  integer, private, dimension(MPI_STATUS_SIZE) :: status

contains

subroutine commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

  integer :: cartcomm, size, rank

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  integer :: dimsize, idim

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cartdim_get(cartcomm, dimsize, ierr)

  if (dimsize /= ndim) then

     if (rank == 0) write(*,*) "commdata: comm dimension ", dimsize, &
          " not equal to ndim = ", ndim
     call MPI_Finalize(ierr)
     stop

  end if

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

  !
  ! Find neighbours for this process
  !

  do idim = 1, ndim

     call MPI_Cart_shift(cartcomm, idim-1, 1, &
          neighbour(dndir, idim), &
          neighbour(updir, idim), ierr)

  end do

end subroutine commdata


subroutine halosendrecv(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idim, idir

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

!  do idim = 1, ndim
!    write(*,*) 'rank ', rank, ', neighbours in dim ', idim, ' are ', &
!                neighbour(dndir, idim), neighbour(updir, idim)
!  end do

! Halo swap

  do irep = 1, nrep

     do idim = 1, ndim
        do idir = 1, ndir

           call MPI_Sendrecv(sendbuf(1, idir, idim), n, MPI_DOUBLE_PRECISION, &
                             neighbour(idir, idim), defsendtag, &
                             recvbuf(1, ndir-idir+1, idim), n, MPI_DOUBLE_PRECISION, &
                             neighbour(ndir-idir+1, idim), defrecvtag, &
                             cartcomm, status, ierr)
        end do
     end do

  end do

end subroutine halosendrecv


subroutine haloredblack(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idir, idim, oddeven

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

  oddeven = mod(sum(coords), 2)

! Halo swap

  do irep = 1, nrep

     do idim = 1, ndim
        do idir = 1, ndir

           if (oddeven == 0) then

! Even processes send down and recv down, send up and recv up

              call MPI_Send(sendbuf(1, idir, idim), &
                            n, MPI_DOUBLE_PRECISION, &
                            neighbour(idir, idim), defsendtag, &
                           cartcomm, ierr)

              call MPI_Recv(recvbuf(1, idir, idim), &
                            n, MPI_DOUBLE_PRECISION, &
                            neighbour(idir, idim), defrecvtag, &
                            cartcomm, status, ierr)

           else
           
! Odd processes recv up and send up, recv down and send down

              call MPI_Recv(recvbuf(1, ndir-idir+1, idim), &
                            n, MPI_DOUBLE_PRECISION, &
                            neighbour(ndir-idir+1, idim), defrecvtag, &
                            cartcomm, status, ierr)

              call MPI_Send(sendbuf(1, ndir-idir+1, idim), &
                            n, MPI_DOUBLE_PRECISION, &
                            neighbour(ndir-idir+1, idim), defsendtag, &
                            cartcomm, ierr)
           end if

        end do
     end do

  end do

end subroutine haloredblack


subroutine haloisendrecvwait(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idim, idir
  integer :: request

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

! Halo swap

  do irep = 1, nrep

     do idim = 1, ndim
        do idir = 1, ndir

           call MPI_Isend(sendbuf(1,idir,idim), &
                          n, MPI_DOUBLE_PRECISION, &
                          neighbour(idir, idim), defsendtag, &
                          cartcomm, request, ierr)

           call MPI_Recv(recvbuf(1, ndir-idir+1,idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(ndir-idir+1, idim), defrecvtag, &
                         cartcomm, status, ierr)

           call MPI_Wait(request, status, ierr)

        end do
     end do

  end do

end subroutine haloisendrecvwait

subroutine haloirecvsendwait(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idim, idir
  integer :: request

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

! Halo swap

  do irep = 1, nrep

     do idim = 1, ndim
        do idir = 1, ndir

           call MPI_Irecv(recvbuf(1,ndir-idir+1,idim), &
                          n, MPI_DOUBLE_PRECISION, &
                          neighbour(ndir-idir+1, idim), defrecvtag, &
                          cartcomm, request, ierr)

           call MPI_Send(sendbuf(1,idir,idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(idir, idim), defsendtag, &
                         cartcomm, ierr)

           call MPI_Wait(request, status, ierr)

        end do
     end do

  end do

end subroutine haloirecvsendwait

subroutine haloirecvisendwaitpair(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idim, idir

  integer, dimension(MPI_STATUS_SIZE, nrequest) :: statuses
  integer :: requests(nrequest)

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

! Halo swap

  do irep = 1, nrep

     do idim = 1, ndim
        do idir = 1, ndir

           call MPI_Irecv(recvbuf(1,ndir-idir+1,idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(ndir-idir+1, idim), defrecvtag, &
                         cartcomm, requests(1), ierr)

           call MPI_Isend(sendbuf(1,idir,idim), &
                          n, MPI_DOUBLE_PRECISION, &
                          neighbour(idir, idim), defsendtag, &
                          cartcomm, requests(2), ierr)

           call MPI_Waitall(nrequest, requests, statuses, ierr)

        end do
     end do

  end do

end subroutine haloirecvisendwaitpair


subroutine haloirecvisendwaitall(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idim, idir, irequest, reqid

  integer, dimension(MPI_STATUS_SIZE, ndir*ndim*nrequest) :: statuses
  integer, dimension(ndir*ndim*nrequest) :: requests

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims


  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

  ! Halo swap

  do irep = 1, nrep

     reqid = 0

     do irequest = 1, nrequest
        do idim = 1, ndim
           do idir = 1, ndir

              reqid = reqid + 1

              if (irequest == rrequest) then

                 call MPI_Irecv(recvbuf(1,ndir-idir+1,idim), &
                                n, MPI_DOUBLE_PRECISION, &
                                neighbour(ndir-idir+1, idim), defrecvtag, &
                                cartcomm, requests(reqid), ierr)

              else

                 call MPI_Isend(sendbuf(1,idir,idim), &
                                n, MPI_DOUBLE_PRECISION, &
                                neighbour(idir, idim), defsendtag, &
                                cartcomm, requests(reqid), ierr)

              end if

           end do
        end do
     end do

     call MPI_Waitall(ndir*ndim*nrequest, requests, &
                      statuses, ierr)

  end do

end subroutine haloirecvisendwaitall

subroutine halopersist(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep, idim, idir, irequest, reqid

  integer, dimension(MPI_STATUS_SIZE, ndir*ndim*nrequest) :: statuses
  integer, dimension(ndir*ndim*nrequest) :: requests

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  ! Need tags as message ordering not maintained
  
  integer :: sendtag, recvtag

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

  ! Initialise persistent comms

  reqid = 0

  do irequest = 1, nrequest

     do idim = 1, ndim
        do idir = 1, ndir

           reqid = reqid + 1

           ! Need to create matching unique send and receive tags
           ! as persistent comms do not respect message ordering

           sendtag = (idim-1)*ndim + (idir-1)*ndir + 1
           recvtag = sendtag

!           write(*,*) 'rank, irequest, idim, idir, sendtag, recvtag = ', &
!                rank, irequest, idim, idir, sendtag, recvtag

           if (irequest == rrequest) then

              call MPI_Recv_init(recvbuf(1,ndir-idir+1,idim), &
                                 n, MPI_DOUBLE_PRECISION, &
                                 neighbour(ndir-idir+1, idim), recvtag, &
                                 cartcomm, requests(reqid), ierr)

           else
              
              call MPI_Send_init(sendbuf(1,idir,idim), &
                                 n, MPI_DOUBLE_PRECISION, &
                                 neighbour(idir, idim), sendtag, &
                                 cartcomm, requests(reqid), ierr)

           end if

        end do
     end do
  end do
  
! now execute halo swaps multiple times

  do irep = 1, nrep

     call MPI_Startall(ndir*ndim*nrequest, requests, ierr)

     call MPI_Waitall(ndir*ndim*nrequest, requests, &
                      statuses, ierr)

  end do

! Release the persistent handles

  do reqid = 1, ndir*ndim*nrequest

     call MPI_Request_free(requests(reqid), ierr)

  end do

end subroutine halopersist

subroutine haloneighboralltoall(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf

  integer :: size, rank
  
  integer :: irep

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

  ! Execute halo swaps multiple times

  do irep = 1, nrep

     call MPI_Neighbor_alltoall(sendbuf, n, MPI_DOUBLE_PRECISION, &
                                recvbuf, n, MPI_DOUBLE_PRECISION, &
                                cartcomm, ierr)
  end do

end subroutine haloneighboralltoall

subroutine inithalodata(rank, sendbuf, recvbuf, nbuf)

  double precision, dimension(nbuf, ndir, ndim) :: sendbuf, recvbuf

  integer :: nbuf, rank

  integer :: idim, idir, ibuf

  ! set receive buffer to have illegal values

 recvbuf(:,:,:) = -1
  
! Set send buffer to have unique values

  do idim = 1, ndim
     do idir = 1, ndir
        do ibuf = 1, nbuf

           sendbuf(ibuf, idir, idim) =   rank*(nbuf*ndir*ndim) &
                                       + (idim-1)*nbuf*ndir  &
                                       + (idir-1)*nbuf + ibuf

!           write(*,*) 'rank = ', rank, ', x(', ibuf, ', ', idir, ', ', idim, &
!                ') = ', sendbuf(ibuf,idir,idim)

        end do
     end do
  end do

end subroutine inithalodata

subroutine checkrecvdata(recvbuf, nbuf, cartcomm)

  integer :: nbuf, cartcomm
  double precision, dimension(nbuf, ndir, ndim) :: recvbuf

  integer :: size, rank
  integer :: idim, idir
  
  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  double precision, dimension(nbuf, ndir, ndim) :: sendbuf, nullbuf

  ! Check that the correct data is sent and received

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)
  
  do idim = 1, ndim
     do idir = 1, ndir

        ! Compute the send halo data for this neighbour

        call inithalodata(neighbour(idir, idim), sendbuf, nullbuf, nbuf)

        if(any(recvbuf(:, idir, idim) /= sendbuf(:, ndir-idir+1, idim))) then
           write(*,*) "checkrecvdata: error for rank, idim, idir = ", &
                rank, idim, idir
!                "r/s = ", recvbuf(:, idir, idim), sendbuf(:, ndir-idir+1, idim) 
        end if

     end do
  end do

end subroutine checkrecvdata

end module haloswap
