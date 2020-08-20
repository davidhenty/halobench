program halobench

  use haloparams
  use benchclock
  use haloswap

  implicit none

  integer, parameter :: nhalotest = 8
  integer :: ihalotest
  logical :: halocalled
  
  !
  ! Set main paramemers: size of each halo and number of repetitions
  !

  integer, parameter :: nbuf = 10000
  integer, parameter :: nrep = 1000

  double precision, dimension(nbuf, ndir, ndim) :: sendbuf, recvbuf

  integer :: rank, size, comm, cartcomm, dblesize
  integer, dimension(ndim) :: dims

  integer, parameter :: mib = 1024*1024

  logical :: reorder = .false.

! Periodic boundaries

  logical, dimension(ndim) :: periods = [.true., .true., .true.]

  double precision :: t0, t1, time, iorate, mibdata

  integer :: ichar
  integer, parameter :: maxstring = 64
  character*(maxstring) :: halotest(nhalotest), divstring

  halotest(1) = "Sendrecv"
  halotest(2) = "Redblack"
  halotest(3) = "Isend / Recv / Wait"
  halotest(4) = "Irecv / Send / Wait"
  halotest(5) = "Irecv / Isend / Wait (pairwise)"
  halotest(6) = "Irecv / Isend / Waitall"
  halotest(7) = "Persistent / Startall / Waitall"
  halotest(8) = "Neighbourhood Collective"

  do ichar = 1, maxstring
     divstring(ichar:ichar) = "-"
  end do

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  dims(:) = 0

! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  mibdata = float(dblesize*nbuf*ndir*ndim)/float(mib)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple haloswap benchmark'
     write(*,*) '-------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*) 'Process grid is (', dims(1), ', ', dims(2), ', ', dims(3), ')'
     write(*,*) 'Each halo contains', nbuf, ' doubles'
     write(*,*)
     write(*,*) 'Total amount of halo data = ', mibdata, ' MiB per process'
     write(*,*)
     write(*,*) 'Number of repetitions = ', nrep
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ' microseconds'
  end if

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

  do ihalotest = 1, nhalotest

     halocalled = .true.

     if (rank == 0) then
        write(*,*)
        write(*,*) trim(halotest(ihalotest))
        write(*,*) divstring(1:len(trim(halotest(ihalotest))))
     end if


     call inithalodata(rank, sendbuf, recvbuf, nbuf)

     call MPI_Barrier(comm, ierr)
     t0 = benchtime()

     select case(ihalotest)

     case(1)
        call halosendrecv(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case(2)

        if (any(mod(dims(:),2) /= 0)) then

           if (rank == 0) then
              write(*,*) "Redblack requires all process dimensions to be even"
              write(*,*) "Skipping this test"
           end if

           halocalled = .false.

        else
           call haloredblack(nrep, sendbuf, recvbuf, nbuf, cartcomm)
        end if

     case(3)
        call haloisendrecvwait(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case(4)
        call haloirecvsendwait(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case(5)
        call haloirecvisendwaitpair(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case(6)
        call haloirecvisendwaitall(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case(7)
        call halopersist(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case(8)
        call haloneighboralltoall(nrep, sendbuf, recvbuf, nbuf, cartcomm)

     case default

        if (rank == 0) then
           write(*,*)
           write(*,*) "Illegal value for ihalotest: ", ihalotest
           write(*,*)
        end if

        halocalled = .false.

     end select

     if (halocalled) then

        call MPI_Barrier(comm, ierr)
        t1 = benchtime()

        time = t1 - t0
        iorate = dble(nrep)*mibdata/time

        if (rank == 0) then
           write(*,*) "Secs = ", time, ", bwidth = ", iorate, " MiB/s"
        end if

        call checkrecvdata(recvbuf, nbuf, cartcomm)

     end if

  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) "Finished"
     write(*,*) "--------"
     write(*,*)
  end if

  call MPI_Finalize(ierr)
  
end program halobench
