subroutine mpiexcu
    USE mpi
    USE GlobalVariables
    IMPLICIT NONE
 
   integer :: i,k,l
   integer :: max_send, max_recv,p,scount,rcount
   integer, allocatable :: cud_send_pos(:), cud_recv_pos(:)
   integer, allocatable :: cud_send_count(:), cud_recv_count(:)
   real(dp) ,allocatable :: cud_recv_data(:,:,:)
   real(dp), allocatable :: cud_send_data(:,:,:)

allocate(cud_send_count(nprocs),cud_recv_count(nprocs))
allocate(cud_send_data(max_pghosts, nprocs,neq))
allocate(cud_recv_data(max_pghosts, nprocs,neq))
allocate(cud_send_pos(nprocs), cud_recv_pos(nprocs))

 cud_send_pos   = 0
 cud_recv_pos   = 0
 cud_send_count = 0
 cud_recv_count = 0

! PACK (same as geometry)
do i = 1, pghosts
  p = neigh_proc(i)
  cud_send_pos(p+1) = cud_send_pos(p+1) + 1
  do k = 1, neq
    cud_send_data(cud_send_pos(p+1), p+1, k) = cu(pp_elem(i), k)
  end do
end do
 cud_send_count = cud_send_pos

! EXCHANGE counts (same)
do p = 0, nprocs-1 
    if (p /= myid) then 
    scount = cud_send_count(p+1) 
    call MPI_Sendrecv(scount, 1, MPI_INTEGER, p, 111, rcount, 1, MPI_INTEGER, p, 111, &
            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) 
    cud_recv_count(p+1) = rcount 
    end if 
end do

! EXCHANGE loop: loop p, then loop k = 1..nvar and call Sendrecv exactly as geometry did:
do p = 0, nprocs-1
  if (p /= myid) then
    if (cud_send_count(p+1) > 0 .or. cud_recv_count(p+1) > 0) then
       do k = 1, neq
         call MPI_Sendrecv(cud_send_data(1,p+1,k), cud_send_count(p+1), MPI_DOUBLE_PRECISION, p, 222, &
                           cud_recv_data(1,p+1,k), cud_recv_count(p+1), MPI_DOUBLE_PRECISION, p, 222, &
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       end do
    end if
  end if
end do

! UNPACK (exactly mirror geometry)
 cud_recv_pos = 0
do i = 1, pghosts
  p = neigh_proc(i)
  cud_recv_pos(p+1) = cud_recv_pos(p+1) + 1
  do k = 1, neq
    cu(dest_ghost(i), k) = cud_recv_data(cud_recv_pos(p+1), p+1, k)
  end do
end do

deallocate( cud_send_count, cud_recv_count,cud_send_data,cud_recv_data,cud_send_pos,cud_recv_pos)
return
end subroutine mpiexcu
