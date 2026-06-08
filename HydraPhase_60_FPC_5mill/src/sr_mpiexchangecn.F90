subroutine mpiexcn
    USE mpi
    USE GlobalVariables
    IMPLICIT NONE
 
   integer :: i,k,l
   integer :: max_send, max_recv,p,scount,rcount
   integer, allocatable :: cnd_send_pos(:), cnd_recv_pos(:)
   integer, allocatable :: cnd_send_count(:), cnd_recv_count(:)
   real(dp) ,allocatable :: cnd_recv_data(:,:,:)
   real(dp), allocatable :: cnd_send_data(:,:,:)

allocate(cnd_send_count(nprocs),cnd_recv_count(nprocs))
allocate(cnd_send_data(max_pghosts, nprocs,neq))
allocate(cnd_recv_data(max_pghosts, nprocs,neq))
allocate(cnd_send_pos(nprocs), cnd_recv_pos(nprocs))

!if(myid==0)print*,"cnd_send_data",cnd_send_data(:,:,:)
!It will exchange cn data

 cnd_send_pos   = 0
 cnd_recv_pos   = 0
 cnd_send_count = 0
 cnd_recv_count = 0
 cnd_send_data  = 0
 cnd_recv_data = 0


!!$acc update host(total) 
! PACK (same as geometry)
!$acc update host(cn(:,:))
do i = 1, pghosts
  p = neigh_proc(i)
  cnd_send_pos(p+1) = cnd_send_pos(p+1) + 1
  do k = 1, neq
    !   if(myid==0.and.total==1)print*,"bef",i,k,pp_elem(i),cn(pp_elem(i), k)
      cnd_send_data(cnd_send_pos(p+1), p+1, k) = cn(pp_elem(i), k)
!    if(myid==0.and.k==11)print*,"A",cn(pp_elem(i), k),pp_elem(i)
  end do
end do

!stop

 cnd_send_count = cnd_send_pos


! EXCHANGE counts (same)
do p = 0, nprocs-1
    if (p /= myid) then
    scount = cnd_send_count(p+1)  
    call MPI_Sendrecv(scount, 1, MPI_INTEGER, p, 111, rcount, 1, MPI_INTEGER, p, 111, &
            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    cnd_recv_count(p+1) = rcount
    end if
end do


! EXCHANGE loop: loop p, then loop k = 1..nvar and call Sendrecv exactly as geometry did:
do p = 0, nprocs-1
  if (p /= myid) then
    if (cnd_send_count(p+1) > 0 .or. cnd_recv_count(p+1) > 0) then
       do k = 1, neq
         call MPI_Sendrecv(cnd_send_data(1,p+1,k), cnd_send_count(p+1), MPI_DOUBLE_PRECISION, p, 222, &
                           cnd_recv_data(1,p+1,k), cnd_recv_count(p+1), MPI_DOUBLE_PRECISION, p, 222, &
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       end do
    end if
  end if
end do
!!$acc update device(cnd_recv_data)


! UNPACK (exactly mirror geometry)
!!$acc serial present(cnd_recv_pos(:))
 cnd_recv_pos = 0
!!$acc end serial

!!$acc parallel loop present(pghosts,neigh_proc(:),cnd_recv_pos(:),dest_ghost(:),cnd_recv_data(:,:,:),cn(:,:)) &
!!$acc private(i,p,k) 
 do i = 1, pghosts
  p = neigh_proc(i)
  cnd_recv_pos(p+1) = cnd_recv_pos(p+1) + 1
!!$acc loop seq  
  do k = 1, neq
 !   if(myid==0)print*,"bef",i,k,dest_ghost(i),cn(dest_ghost(i), k)
    cn(dest_ghost(i), k) = cnd_recv_data(cnd_recv_pos(p+1), p+1, k)
 !   if(myid==0)print*,"aft",i,k,dest_ghost(i),cn(dest_ghost(i), k)
  end do
end do


!$acc update device(cn(nepg:ntot,:))

deallocate( cnd_send_count, cnd_recv_count,cnd_send_data,cnd_recv_data,cnd_send_pos,cnd_recv_pos)
return
end subroutine mpiexcn
