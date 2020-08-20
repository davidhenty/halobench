module haloswap

  use mpi
  implicit none

  integer, parameter :: ndim = 3
  integer, parameter :: nneigh = 2

  integer, parameter, private :: neighdn = 1
  integer, parameter, private :: neighup = 2

  integer, parameter, private :: nrequest = 2
  integer, parameter, private :: rrequest = 1
  integer, parameter, private :: srequest = 2

contains

subroutine halosendrecv(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n
  double precision, dimension(n, nneigh, ndim) :: sendbuf, recvbuf

  integer :: i, idim, irequest
  integer :: cartcomm, ierr, rank, size
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer, dimension(nneigh, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  integer :: sendtag = 1
  integer :: recvtag = 1

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

!
! Find neighbours for this process
!

  do idim = 1, ndim

     call MPI_Cart_shift(cartcomm, idim-1, 1, &
                         neighbour(neighdn, idim), &
                         neighbour(neighup, idim), ierr)

  end do

!  do idim = 1, ndim
!     write(*,*) 'rank ', rank, ', neighbours in dim ', idim, ' are ', &
!                neighbour(neighdn, idim), neighbour(neighup, idim)
!  end do

! Halo swap

  do i = 1, nrep

     do idim = 1, ndim

        call MPI_Sendrecv(sendbuf(1, neighup, idim), n, MPI_DOUBLE_PRECISION, &
                          neighbour(neighup, idim), sendtag, &
                          recvbuf(1, neighdn, idim), n, MPI_DOUBLE_PRECISION, &
                          neighbour(neighdn, idim), recvtag, &
                          cartcomm, status, ierr)

        call MPI_Sendrecv(sendbuf(1, neighdn, idim), n, MPI_DOUBLE_PRECISION, &
                          neighbour(neighdn, idim), sendtag, &
                          recvbuf(1, neighup, idim), n, MPI_DOUBLE_PRECISION, &
                          neighbour(neighup, idim), recvtag, &
                          cartcomm, status, ierr)

     end do

! do some fake "work"

     sendbuf(:,:,:) = recvbuf(:,:,:)

  end do

end subroutine halosendrecv

subroutine haloredblack(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n
  double precision, dimension(n, nneigh, ndim) :: sendbuf, recvbuf

  integer :: i, idim, irequest, irb, oddeven
  integer :: cartcomm, ierr, rank, size
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer, dimension(nneigh, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  integer :: sendtag = 1
  integer :: recvtag = 1

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

  oddeven = mod(sum(coords), 2)

!
! Find neighbours for this process
!

  do idim = 1, ndim

     call MPI_Cart_shift(cartcomm, idim-1, 1, &
                         neighbour(neighdn, idim), &
                         neighbour(neighup, idim), ierr)

  end do

!  do idim = 1, ndim
!     write(*,*) 'rank ', rank, ', neighbours in dim ', idim, ' are ', &
!                neighbour(neighdn, idim), neighbour(neighup, idim)
!  end do

! Halo swap

  do i = 1, nrep

     do idim = 1, ndim

        if (oddeven == 0) then

! Even processes send up and recv up, send down and recv down

!              write(*,*) 'Rank ', rank, ' Send to ', neighbour(neighup, idim)

           call MPI_Send(sendbuf(1, neighup, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighup, idim), sendtag, &
                         cartcomm, ierr)

!              write(*,*) 'Rank ', rank, ' rec from ', neighbour(neighup, idim)

           call MPI_Recv(recvbuf(1, neighup, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighup, idim), recvtag, &
                         cartcomm, status, ierr)

!              write(*,*) 'Rank ', rank, ' Send to ', neighbour(neighdn, idim)

           call MPI_Send(sendbuf(1, neighdn, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighdn, idim), sendtag, &
                         cartcomm, ierr)

!              write(*,*) 'Rank ', rank, ' rec from ', neighbour(neighdn, idim)

           call MPI_Recv(recvbuf(1, neighdn, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighdn, idim), recvtag, &
                         cartcomm, status, ierr)
        
        else

! Odd processes recv down and send down, recv up and send  up

!              write(*,*) 'Rank ', rank, ' rec from ', neighbour(neighdn, idim)

           call MPI_Recv(recvbuf(1, neighdn, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighdn, idim), recvtag, &
                         cartcomm, status, ierr)

!              write(*,*) 'Rank ', rank, ' Send to ', neighbour(neighdn, idim)

           call MPI_Send(sendbuf(1, neighdn, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighdn, idim), sendtag, &
                         cartcomm, ierr)

!              write(*,*) 'Rank ', rank, ' rec from ', neighbour(neighup, idim)

           call MPI_Recv(recvbuf(1, neighup, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighup, idim), recvtag, &
                         cartcomm, status, ierr)

!              write(*,*) 'Rank ', rank, ' Send to ', neighbour(neighup, idim)

           call MPI_Send(sendbuf(1, neighup, idim), &
                         n, MPI_DOUBLE_PRECISION, &
                         neighbour(neighup, idim), sendtag, &
                         cartcomm, ierr)

        end if
     end do

! do some fake "work"

     sendbuf(:,:,:) = recvbuf(:,:,:)

  end do

end subroutine haloredblack

subroutine halowaitall(nrep, sendbuf, recvbuf, n, cartcomm)

  integer :: nrep, n
  double precision, dimension(n, nneigh, ndim) :: sendbuf, recvbuf

  integer :: i, idim, ineigh, irequest, reqid, irep
  integer :: cartcomm, ierr, rank, size
  integer, dimension(MPI_STATUS_SIZE, nneigh*ndim*nrequest) :: statuses
  integer, dimension(nneigh*ndim*nrequest) :: requests

  integer, dimension(nneigh, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  integer :: sendtag = 1
  integer :: recvtag = 1

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

!
! Find neighbours for this process
!

  do idim = 1, ndim

     call MPI_Cart_shift(cartcomm, idim-1, 1, &
                         neighbour(neighdn, idim), &
                         neighbour(neighup, idim), ierr)

  end do

! Halo swap

  do irep = 1, nrep

     reqid = 0

     do irequest = 1, nrequest

        do idim = 1, ndim

           do ineigh = 1, nneigh

              reqid = reqid + 1

              if (irequest == rrequest) then

                 call MPI_Irecv(recvbuf(1,ineigh,idim), &
                                n, MPI_DOUBLE_PRECISION, &
                                neighbour(ineigh, idim), recvtag, &
                                cartcomm, requests(reqid), ierr)

              else

                 call MPI_Isend(sendbuf(1,ineigh,idim), &
                                n, MPI_DOUBLE_PRECISION, &
                                neighbour(ineigh, idim), sendtag, &
                                cartcomm, requests(reqid), ierr)

              end if

           end do
        end do
     end do

     call MPI_Waitall(nneigh*ndim*nrequest, requests, &
                      statuses, ierr)

! do some fake "work"

     sendbuf(:,:,:) = recvbuf(:,:,:)

  end do

end subroutine halowaitall

end module haloswap
