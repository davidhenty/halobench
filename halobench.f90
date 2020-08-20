program benchio

  use benchclock
  use haloswap

  implicit none

! Set local buffer size

  integer, parameter :: n = 10000

  integer :: p1, p2, p3, nrep, i, ineigh, idim

  double precision, dimension(n, nneigh, ndim) :: sendbuf, recvbuf

  integer :: rank, size, ierr, comm, cartcomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: mib = 1024*1024

  logical :: reorder = .false.
  logical, dimension(ndim) :: periods = [.true., .true., .true.]

  double precision :: t0, t1, time, iorate, mibdata

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  dims(:) = 0

! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

  p1 = dims(1)
  p2 = dims(2)
  p3 = dims(3)

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  mibdata = float(dblesize*n*nneigh*ndim)/float(mib)

  nrep = 1000

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple haloswap benchmark'
     write(*,*) '-------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*) 'Process grid is (', p1, ', ', p2, ', ', p3, ')'
     write(*,*) 'Buffer size is ', n
     write(*,*)
     write(*,*) 'Total amount of halo data = ', mibdata, ' MiB'
     write(*,*)
     write(*,*) 'Number of repetitions = ', nrep
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ', usecs'
  end if
  
  dims(1) = p1
  dims(2) = p2
  dims(3) = p3

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

! Set data to illegal values

  recvbuf(:,:,:) = -1
  sendbuf(:,:,:) = -1
  
! Set halo data core to have unique values

  call MPI_Cart_coords(cartcomm, rank, ndim, coords, ierr)
  
  do idim = 1, ndim
     do ineigh = 1, nneigh
        do i = 1, n

           sendbuf(i, ineigh, idim) =   rank*(n*nneigh*ndim) &
                                      + (idim-1)*n*nneigh  &
                                      + (ineigh-1)*n + i

!           write(*,*) 'rank = ', rank, ', x(', i, ', ', ineigh, ', ', idim, &
!                ') = ', sendbuf(i,ineigh,idim)

        end do
     end do
  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Sendrecv'
     write(*,*) '--------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call halosendrecv(nrep, sendbuf, recvbuf, n, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Sendrecv time = ', time, ', IO rate = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Waitall'
     write(*,*) '-------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call halowaitall(nrep, sendbuf, recvbuf, n, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Waitall time =  ', time, ', IO rate = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Redblack'
     write(*,*) '--------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call haloredblack(nrep, sendbuf, recvbuf, n, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Redblack time =  ', time, ', IO rate = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Finished'
     write(*,*) '--------'
     write(*,*)
  end if

  call MPI_Finalize(ierr)
  
end program benchio
