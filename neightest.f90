program neightest

  use iso_fortran_env
  use mpi
  
  implicit none

  integer, parameter :: ndir = 2    ! down and up
  integer, parameter :: ndim = 2
  integer, parameter :: nbuf = ndim*ndir  

  integer, dimension(ndim) :: dims, coords
  logical, dimension(ndim) :: periods
  logical :: reorder = .false.

  integer, dimension(nbuf) :: sendbuf, recvbuf

  integer :: rank, size, comm, cartcomm, ierr, ibuf

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

! Set processor grid with periodic bcs

  dims(:) = 0
  periods(:) = .true.

  call MPI_Dims_create(size, ndim, dims, ierr)

  if (rank == 0) then
     write(*,*) 'Process grid is (', dims(1), ' x ', dims(2), ')'
     write(*,*)
  end if

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)
  call MPI_Cart_coords(cartcomm, rank, ndim, coords, ierr)
  
  do ibuf = 1, nbuf
     sendbuf(ibuf) = rank*nbuf + ibuf
  end do

  recvbuf(:) = -1

  write(*,*) "Before: rank; coords; sendbuf = ", rank, "; (", coords(:), "); ", sendbuf(:)

  flush(output_unit)
  call MPI_Barrier(Comm, ierr)
  
  call MPI_Neighbor_alltoall(sendbuf, 1, MPI_INTEGER, &
                             recvbuf, 1, MPI_INTEGER, &
                             cartcomm, ierr)

  write(*,*) "After:  rank; coords; recvbuf = ", rank, "; (", coords(:), "); ", recvbuf(:)

  flush(output_unit)
  call MPI_Barrier(Comm, ierr)

  call MPI_Finalize(ierr)

end program neightest
