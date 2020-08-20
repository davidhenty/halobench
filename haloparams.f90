module haloparams

  use mpi

  integer, parameter :: ndim = 3
  integer, parameter :: ndir = 2    ! down and up

  integer, parameter :: dndir = 1
  integer, parameter :: updir = 2   ! ordering needs to be like this for neighourhood collectives

  integer :: ierr

end module haloparams
