subroutine mpiexdel
  use mpi
  use GlobalVariables
  implicit none

  integer :: i, p, k, scount, rcount
   INTEGER, ALLOCATABLE :: del_send_count(:), del_recv_count(:)
   INTEGER, ALLOCATABLE :: del_send_pos(:),   del_recv_pos(:)
   REAL(DP), ALLOCATABLE :: del_send_data(:,:,:), del_recv_data(:,:,:)
   
  allocate(del_send_count(1:nprocs), del_recv_count(1:nprocs))
  allocate(del_send_pos(1:nprocs),   del_recv_pos(1:nprocs))
  allocate(del_send_data(max_pghosts, nprocs, nvar))
  allocate(del_recv_data(max_pghosts, nprocs, nvar))

!!$acc serial present(del_send_pos,del_recv_pos,del_send_count,del_recv_count)
  del_send_pos   = 0
  del_recv_pos   = 0
  del_send_count = 0
  del_recv_count = 0
!!$acc end serial
  !---------------- PACK ----------------

!$acc update host(delugx,delugy, delugz,delvgx, delvgy, delvgz,delwgx, delwgy, delwgz, &
!$acc deltgx, deltgy, deltgz,delupx, delupy, &
!$acc delupz,delvpx, delvpy, delvpz, delwpx, delwpy, delwpz,deltpx, deltpy, deltpz, &
!$acc delpgx, delpgy, delpgz,delppx, delppy, delppz,delphigx, delphigy, delphigz, &
!$acc delTkgx,delTkgy,delTkgz,delTegx,delTegy,delTegz)

!!$acc parallel loop private(i,p) present(delugx, delugy, delugz,delvgx, delvgy, delvgz, &
!!$acc delwgx, delwgy, delwgz,deltgx, deltgy, deltgz,delupx, delupy, delupz, &
!!$acc delvpx, delvpy, delvpz, delwpx, delwpy, delwpz,deltpx, deltpy, deltpz, &
!!$acc delpgx, delpgy, delpgz,delppx, delppy, delppz,delphigx, delphigy, delphigz, &
!!$acc pghosts,neigh_proc(:),del_send_pos(:),pp_elem(:))  
  do i = 1, pghosts
    p = neigh_proc(i)                 ! neighbor rank (0..nprocs-1)
    del_send_pos(p+1) = del_send_pos(p+1) + 1

  del_send_data(del_send_pos(p+1), p+1,  1) = delugx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  2) = delugy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  3) = delugz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  4) = delvgx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  5) = delvgy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  6) = delvgz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  7) = delwgx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  8) = delwgy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1,  9) = delwgz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 10) = deltgx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 11) = deltgy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 12) = deltgz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 13) = delupx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 14) = delupy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 15) = delupz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 16) = delvpx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 17) = delvpy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 18) = delvpz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 19) = delwpx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 20) = delwpy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 21) = delwpz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 22) = deltpx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 23) = deltpy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 24) = deltpz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 25) = delpgx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 26) = delpgy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 27) = delpgz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 28) = delppx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 29) = delppy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 30) = delppz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 31) = delphigx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 32) = delphigy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 33) = delphigz(pp_elem(i))
#ifdef KE_TURB
  del_send_data(del_send_pos(p+1), p+1, 34) = delTkgx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 35) = delTkgy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 36) = delTkgz(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 37) = delTegx(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 38) = delTegy(pp_elem(i))
  del_send_data(del_send_pos(p+1), p+1, 39) = delTegz(pp_elem(i))
#endif
  end do
!!$acc end parallel loop

!!$acc serial present(del_send_count,del_send_pos)
  del_send_count = del_send_pos
!!$acc end serial

!!$acc update host(del_send_count,del_send_data)
  !---------------- EXCHANGE COUNTS ----------------
  do p = 0, nprocs-1
    if (p /= myid) then
      scount = del_send_count(p+1)
      call MPI_sendrecv(scount, 1, MPI_INTEGER, p, 353, &
                        rcount, 1, MPI_INTEGER, p, 353, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      del_recv_count(p+1) = rcount
    end if
  end do
!print*, "counts", scount, rcount
  !---------------- EXCHANGE DATA (one message per neighbor) ----------------
  do p = 0, nprocs-1
    if (p /= myid) then
      if (del_send_count(p+1) > 0 .or. del_recv_count(p+1) > 0) then
              do k = 1,nvar
        call MPI_sendrecv( del_send_data(1,p+1,k), del_send_count(p+1), MPI_DOUBLE_PRECISION, p, 454, &
                           del_recv_data(1,p+1,k), del_recv_count(p+1), MPI_DOUBLE_PRECISION, p, 454, &
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr )
                   end do
      end if
    end if
  end do
!!$acc update device(del_recv_data)  

!!$acc serial present(del_recv_pos)
  del_recv_pos = 0
!!$acc end serial  
  !---------------- UNPACK ----------------
!!$acc parallel loop private(i,p) present(dest_ghost,del_recv_pos,neigh_proc,pghosts,delugx, &
!!$acc delugy, delugz,delvgx, delvgy, delvgz, &
!!$acc delwgx, delwgy, delwgz,deltgx, deltgy, deltgz,delupx, delupy, delupz, &
!!$acc delvpx, delvpy, delvpz, delwpx, delwpy, delwpz,deltpx, deltpy, deltpz, &
!!$acc delpgx, delpgy, delpgz,delppx, delppy, delppz,delphigx, delphigy, delphigz)

  do i = 1, pghosts
    p = neigh_proc(i)
    del_recv_pos(p+1) = del_recv_pos(p+1) + 1

    delugx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 1)
    delugy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 2)
    delugz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 3)

    delvgx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 4)
    delvgy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 5)
    delvgz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 6)

    delwgx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 7)
    delwgy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 8)
    delwgz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1, 9)

    deltgx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,10)
    deltgy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,11)
    deltgz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,12)

    delupx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,13)
    delupy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,14)
    delupz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,15)

    delvpx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,16)
    delvpy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,17)
    delvpz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,18)

    delwpx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,19)
    delwpy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,20)
    delwpz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,21)

    deltpx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,22)
    deltpy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,23)
    deltpz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,24)

    delpgx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,25)
    delpgy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,26)
    delpgz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,27)

    delppx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,28)
    delppy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,29)
    delppz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,30)

    delphigx(dest_ghost(i)) = del_recv_data(del_recv_pos(p+1), p+1,31)
    delphigy(dest_ghost(i)) = del_recv_data(del_recv_pos(p+1), p+1,32)
    delphigz(dest_ghost(i)) = del_recv_data(del_recv_pos(p+1), p+1,33)
#ifdef KE_TURB    
    delTkgx(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,34)
    delTkgy(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,35)
    delTkgz(dest_ghost(i))   = del_recv_data(del_recv_pos(p+1), p+1,36)

    delTegx(dest_ghost(i)) = del_recv_data(del_recv_pos(p+1), p+1,37)
    delTegy(dest_ghost(i)) = del_recv_data(del_recv_pos(p+1), p+1,38)
    delTegz(dest_ghost(i)) = del_recv_data(del_recv_pos(p+1), p+1,39)
#endif   
  end do
  
!$acc update device(delugx,delugy, delugz,delvgx, delvgy, delvgz,delwgx, delwgy, delwgz,deltgx, deltgy, deltgz,delupx, &
!$acc delupy, delupz, delvpx, delvpy, delvpz, delwpx, delwpy, delwpz,deltpx, deltpy, deltpz, &
!$acc delpgx, delpgy, delpgz,delppx, delppy, delppz,delphigx, delphigy, delphigz, &
!$acc delTkgx,delTkgy,delTkgz,delTegx,delTegy,delTegz)

deallocate(del_send_data, del_recv_data, del_send_pos, del_recv_pos, del_send_count, del_recv_count)
return
end subroutine mpiexdel

