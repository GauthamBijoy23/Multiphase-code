SUBROUTINE geometry
    USE mpi
    USE GlobalVariables
    USE MatrixOps
    USE run_tioga, ONLY : iblankcells, nc_t
    IMPLICIT NONE
    
 !--------------------------------------------------------------------
! Local indices / counters
!--------------------------------------------------------------------
INTEGER :: i, j, k               ! generic loop indices
INTEGER :: iside                 ! face / side index
INTEGER :: n, n1, n2, n3, n4     ! generic integer helpers
INTEGER :: nel                   ! number of elements
INTEGER :: neigh                 ! neighbour id

!--------------------------------------------------------------------
! Problem / partition sizes
!--------------------------------------------------------------------
INTEGER :: np    ! Parent cell id 
INTEGER :: ng    ! Ghost cells id 
INTEGER :: ns    ! Side id 
INTEGER :: nt    ! Type for boundary condition 

!--------------------------------------------------------------------
! Geometric / metric quantities (double precision)
!--------------------------------------------------------------------
REAL(dp) :: dxyA, dxyB, dxyC
REAL(dp) :: dyzA, dyzB, dyzC
REAL(dp) :: dzxA, dzxB, dzxC

REAL(dp) :: v6        ! volume * 6 (or other scalar)
REAL(dp) :: xcface    ! x-coord of face centre
REAL(dp) :: ycface    ! y-coord of face centre
REAL(dp) :: zcface    ! z-coord of face centre

REAL(dp) :: dl12, dl13, dl14
REAL(dp) :: dl23, dl42, dl43
REAL(dp) :: dlmin      ! minimum edge length

REAL(dp) :: x1, y1, z1
REAL(dp) :: x2, y2, z2
REAL(dp) :: x3, y3, z3
REAL(dp) :: x4, y4, z4
REAL(dp) :: xcenter,ycenter,zcenter 

!--------------------------------------------------------------------
! File / I/O names and directories
!--------------------------------------------------------------------
CHARACTER(LEN=256) :: fname_grid    ! mesh file name
CHARACTER(LEN=256) :: fname_bc      ! boundary condition file
CHARACTER(LEN=256) :: fname_procbc  ! proc bc file
CHARACTER(LEN=256) :: geometry_dir  ! directory for geometry files

!--------------------------------------------------------------------
! Miscellaneous integers / flags
!--------------------------------------------------------------------
INTEGER :: p, dummy
INTEGER :: max_pghost    ! maximum number of physical ghosts

!--------------------------------------------------------------------
! MPI / communication bookkeeping
!--------------------------------------------------------------------
INTEGER :: r, s, v, rv
INTEGER :: status(MPI_STATUS_SIZE)
INTEGER :: scount, rcount

!--------------------------------------------------------------------
! File / directory checks
!--------------------------------------------------------------------
LOGICAL :: dir_exists

  INTEGER, ALLOCATABLE :: g_send_count(:), g_recv_count(:)
  INTEGER, ALLOCATABLE :: g_send_pos(:),   g_recv_pos(:)
  REAL(DP), ALLOCATABLE :: g_send_data(:,:,:), g_recv_data(:,:,:)

    !Build filenames for this rank
    geometry_dir = 'geometry_files/'
    
if (myid == 0) then
 inquire(file=geometry_dir, exist=dir_exists)
 if (.not. dir_exists) call system("mkdir -p " // trim(geometry_dir))
end if   
 
! processor-specific grid file
write(fname_grid,'(A,"grid_proc",I4.4,".dat")')  trim(geometry_dir), myid
open(unit=10, file=fname_grid, status="old")

! processor-specific BC file
write(fname_bc,'(A,"bc_proc",I4.4,".in")') trim(geometry_dir), myid
open(unit=11, file=fname_bc, status="old")

! processor-specific proc BC file
write(fname_procbc,'(A,"proc_bc",I4.4,".in")') trim(geometry_dir), myid
open(unit=12, file=fname_procbc, status="old")
          
    ! Reading nodes and neles 
    ! Rescale if needed 
    read(10,*)nodes,neles
    do i = 1, nodes
    read(10,*)n,x(n),y(n),z(n),btag(n) !(\mod3)
!        x(n)=x(n)*0.001d0
!        y(n)=y(n)*0.001d0
!        z(n)=z(n)*0.001d0
    enddo
    
     xmax = -1.0D30
     ymax = -1.0D30
     zmax = -1.0D30
     ymin =  1.0D30
     zmin =  1.0D30
 !    n1 = 0, n2 = 0, n3 = 0, n4 = 0 
     do n = 1, nodes
       xmax = max(xmax, x(n))
       ymax = max(ymax, y(n))
       zmax = max(zmax, z(n))
       ymin = min(ymin, y(n))
       zmin = min(zmin, z(n))
     end do
        
    do i = 1, neles
      read(10, *) nel, nod(nel, 1), nod(nel, 2), nod(nel, 3), nod(nel, 4), nc1(nel), nc2(nel), nc3(nel), nc4(nel)
    end do

open(unit=15,file='checking_nodes/iblankcellslist.dat',status='replace')
    do i = 1, neles              ! Replacing iblank cell read from files with iblankcells from tioga (after interpolation)
       iblank(i)=iblankcells(i)
       !if(iblankcells(i)==1) then
       write(15,*)i,iblankcells(i)
       !endif
    enddo
close(15)

    read(11, *) nghosts 
    do i = 1, nghosts
      read(11, *) nparent(i), nghost(i), ntype(i), nside(i)
    end do
    
    read(12,*) pghosts, max_pghost  ! number of partition ghost elements 
    do i = 1, pghosts 
      read(12,*) pp_elem(i), dummy, neigh_proc(i), dummy, dummy,dummy,dummy, dest_ghost(i)
    end do

    close(10)
    close(11)
    close(12)

!  print*, "pghost" , dest_ghost(:)

    do i = 1, neles
      do iside = 1, 4
        select case (iside)
          case (1)
                    n1 = nod(i, 2); n2 = nod(i, 3); n3 = nod(i, 4); neigh = nc1(i)
                case (2)
                    n1 = nod(i, 4); n2 = nod(i, 3); n3 = nod(i, 1); neigh = nc2(i)
                case (3)
                    n1 = nod(i, 2); n2 = nod(i, 4); n3 = nod(i, 1); neigh = nc3(i)
                case (4)
                    n1 = nod(i, 3); n2 = nod(i, 2); n3 = nod(i, 1); neigh = nc4(i)
            end select
            
            dxyA = (x(n1) - x(n2)) * (y(n1) + y(n2))
            dxyB = (x(n2) - x(n3)) * (y(n2) + y(n3))
            dxyC = (x(n3) - x(n1)) * (y(n3) + y(n1))
            
            dyzA = (y(n1) - y(n2)) * (z(n1) + z(n2))
            dyzB = (y(n2) - y(n3)) * (z(n2) + z(n3))
            dyzC = (y(n3) - y(n1)) * (z(n3) + z(n1))
            
            dzxA = (z(n1) - z(n2)) * (x(n1) + x(n2))
            dzxB = (z(n2) - z(n3)) * (x(n2) + x(n3))
            dzxC = (z(n3) - z(n1)) * (x(n3) + x(n1))
            
            select case (iside)
                case (1)
                    sc1x(i) = 0.5d0 * (dyzA + dyzB + dyzC)
                    sc1y(i) = 0.5d0 * (dzxA + dzxB + dzxC)
                    sc1z(i) = 0.5d0 * (dxyA + dxyB + dxyC)
                case (2)
                    sc2x(i) = 0.5d0 * (dyzA + dyzB + dyzC)
                    sc2y(i) = 0.5d0 * (dzxA + dzxB + dzxC)
                    sc2z(i) = 0.5d0 * (dxyA + dxyB + dxyC)
                case (3)
                    sc3x(i) = 0.5d0 * (dyzA + dyzB + dyzC)
                    sc3y(i) = 0.5d0 * (dzxA + dzxB + dzxC)
                    sc3z(i) = 0.5d0 * (dxyA + dxyB + dxyC)
                case (4)
                    sc4x(i) = 0.5d0 * (dyzA + dyzB + dyzC)
                    sc4y(i) = 0.5d0 * (dzxA + dzxB + dzxC)
                    sc4z(i) = 0.5d0 * (dxyA + dxyB + dxyC)
            end select
        end do
!        if(myid.eq.0) print*, myid, "scx",i, sc1x(i), sc2x(i), sc3x(i)
      
        n1 = nod(i,1)
        n2 = nod(i,2)
        n3 = nod(i,3)
        n4 = nod(i,4)

        x1 = x(n1); x2 = x(n2); x3 = x(n3); x4 = x(n4)
        y1 = y(n1); y2 = y(n2); y3 = y(n3); y4 = y(n4)
        z1 = z(n1); z2 = z(n2); z3 = z(n3); z4 = z(n4)

        ! Compute cell center
        xcel(i) = 0.25D0 * (x1 + x2 + x3 + x4)
        ycel(i) = 0.25D0 * (y1 + y2 + y3 + y4)
        zcel(i) = 0.25D0 * (z1 + z2 + z3 + z4)
        !if(myid==2.and.i==2) print*, "geom xcel", xcel(nc4(i))       
        ! Compute xmax value for the outflow boundary condition
      
        ! Compute minimum length for timestep
        dl12 = DSQRT((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        dl13 = DSQRT((x1-x3)**2 + (y1-y3)**2 + (z1-z3)**2)
        dl14 = DSQRT((x1-x4)**2 + (y1-y4)**2 + (z1-z4)**2)
        dl23 = DSQRT((x2-x3)**2 + (y2-y3)**2 + (z2-z3)**2)
        dl42 = DSQRT((x4-x2)**2 + (y4-y2)**2 + (z4-z2)**2)
        dl43 = DSQRT((x4-x3)**2 + (y4-y3)**2 + (z4-z3)**2)

        dl(i) = MIN(dl12, dl13, dl14, dl23, dl42, dl43)

        dlmin = MIN(dlmin, dl(i))
        
        X1=X(N1); X2=X(N2); X3=X(N3); X4=X(N4)
        Y1=Y(N1); Y2=Y(N2); Y3=Y(N3); Y4=Y(N4)
        Z1=Z(N1); Z2=Z(N2); Z3=Z(N3); Z4=Z(N4)

        ! Compute volume using determinant method
        v6 = DET(x2, x3, x4, y2, y3, y4, z2, z3, z4) - &
             DET(x1, x3, x4, y1, y3, y4, z1, z3, z4) + &
             DET(x1, x2, x4, y1, y2, y4, z1, z2, z4) - &
             DET(x1, x2, x3, y1, y2, y3, z1, z2, z3)

        vol(i) = ABS(v6) / 6.0D0
        
                
!if(myid==0)print '(A,I6,4F12.6)', "i,xcel,ycel,zcel,vol ", i, xcel(i), ycel(i), zcel(i), vol(i)

        

        IF (vol(i) <= 0.0D0) PRINT *, "Negative volume at element:", i, vol(i)
    END DO
     !if (myid.eq.2) print*, "nside list", nside
    !========================
    ! Ghost cell calculations
    !========================
    !print*, "nghosts", myid, nghosts
	DO i = 1, nghosts
        np = nparent(i)
        ng = nghost(i)
        ns = nside(i)
        nt = ntype(i)
        
        !if (myid.eq.1) print*, "ns nt", i,np,ng,ns, nt

        ! Determine ghost cell face
        SELECT CASE (ns)
            CASE (1)
                n1 = nod(np,2); n2 = nod(np,3); n3 = nod(np,4)
            CASE (2)
                n1 = nod(np,4); n2 = nod(np,3); n3 = nod(np,1)
            CASE (3)
                n1 = nod(np,2); n2 = nod(np,4); n3 = nod(np,1)
            CASE (4)
                n1 = nod(np,3); n2 = nod(np,2); n3 = nod(np,1)
        END SELECT

        x1 = x(n1); x2 = x(n2); x3 = x(n3)
        y1 = y(n1); y2 = y(n2); y3 = y(n3)
        z1 = z(n1); z2 = z(n2); z3 = z(n3)

        xcface = (x1 + x2 + x3) / 3.0D0
        ycface = (y1 + y2 + y3) / 3.0D0
        zcface = (z1 + z2 + z3) / 3.0D0

        ! Compute ghost cell center
        xcel(ng) = 2.0D0 * xcface - xcel(np)
        ycel(ng) = 2.0D0 * ycface - ycel(np)
        zcel(ng) = 2.0D0 * zcface - zcel(np)

        ! Assign volume to ghost cell
        vol(ng) = vol(np)
        
        !!!--- for neighbour's neighbour calculation -----
         nc1(ng)=ng
         nc2(ng)=ng
         nc3(ng)=ng
         nc4(ng)=ng
        !PRINT*, "BEF"
        
! if(myid==0)print '(A,I6,4F12.6)', "ng,xcel,ycel,zcel,vol ", ng, xcel(ng), ycel(ng), zcel(ng), vol(ng)
    END DO

!================================
!Calulate d for wave transmissive
!================================    
! Here need to do it for all boundary
do i = 1, nghosts

  np = NPARENT(i)
  ng = NGHOST(i) 
  nt = ntype(i)
  
   n1 = nod(np, 1)
   n2 = nod(np, 2)
   n3 = nod(np, 3)
   n4 = nod(np, 4)
   
    
if(nt == 3) then  
    x1 = x(n1)
    x2 = x(n2)
    x3 = x(n3)
    x4 = x(n4)
   
    
    xcenter=(x1+x2+x3+x4)*0.25d0      
               
    dd(np) = dabs(xmax-xcenter)  

elseif (nt == 51) then     
 
        y1 = y(n1)
        y2 = y(n2)
        y3 = y(n3)
        y4 = y(n4)
  
        ycenter=(y1+y2+y3+y4)*0.25d0    
       
       if(ycenter.gt.0)  then
        dd(np) = dabs(ymax-ycenter)
       else
        dd(np) =  dabs(ymin-ycenter)
       endif
 
elseif (nt == 53) then   
      z1 = z(n1)
        z2 = z(n2)
        z3 = z(n3)
        z4 = z(n4)
  
        zcenter=(z1+z2+z3+z4)*0.25d0    
       
       if(zcenter.gt.0.1) then 
       dd(np) = dabs(zmax-zcenter)
       else
       dd(np) =  dabs(zmin-zcenter)
       endif          
end if     
end do

allocate(g_send_count(1:nprocs), g_recv_count(1:nprocs))
allocate(g_send_data(max_pghosts, nprocs,4))
allocate(g_recv_data(max_pghosts, nprocs,4))
allocate(g_send_pos(1:nprocs), g_recv_pos(1:nprocs))

g_send_pos   = 0
g_recv_pos   = 0
g_send_count = 0
g_recv_count = 0

!-----------------------------------------------
! 1. PACK: fill send buffers
!-----------------------------------------------
do i = 1, pghosts
    p = neigh_proc(i)          ! destination rank for this ghost
    g_send_pos(p+1) = g_send_pos(p+1) + 1
    g_send_data( g_send_pos(p+1),p+1,1) = xcel(pp_elem(i))
    g_send_data( g_send_pos(p+1),p+1,2) = ycel(pp_elem(i))
    g_send_data( g_send_pos(p+1),p+1,3) = zcel(pp_elem(i))
    g_send_data( g_send_pos(p+1),p+1,4) = vol(pp_elem(i))

  ! if(myid.eq.0)print *, "PACK >> Rank", myid, "-> sending to", p, &
  !           "ghost_id(i)=", i, "pp_elem(i)=", pp_elem(i), &
  !           " value=", xcel(pp_elem(i)), " g_send_pos(p)=", g_send_pos(p+1)
end do

g_send_count = g_send_pos
!print *, "PACK SUMMARY >> Rank", myid, " g_send_count=", g_send_count

!print*, "send buffer packed on rank", myid, " counts=", g_send_count


!-----------------------------------------------
! 2. Exchange counts first (so each rank knows how much to expect)
!-----------------------------------------------

do p = 0, nprocs-1
    if (p /= myid) then
    scount = g_send_count(p+1)
!    print*, "scount", scount
        call MPI_Sendrecv(scount, 1, MPI_INTEGER, p, 111, &
                          rcount, 1, MPI_INTEGER, p, 111, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                          g_recv_count(p+1) = rcount
!                                  print *, "COUNT EXCH >> Rank", myid, " partner=", p, &
!                 " g_send_count=", g_send_count(p+1), " g_recv_count=", g_recv_count(p+1)
    end if
end do

!-----------------------------------------------
! 3. Exchange actual data
!-----------------------------------------------
do p = 0, nprocs-1
    if (p /= myid) then
        if (g_send_count(p+1) > 0 .or. g_recv_count(p+1) > 0) then
                do k = 1,4
       !  print *, "DATA EXCH START >> Rank", myid, " with partner", p, &
        !             " g_send_count=", g_send_count(p), " g_recv_count=", g_recv_count(p)
            call MPI_Sendrecv(g_send_data(1,p+1,k), g_send_count(p+1), MPI_DOUBLE_PRECISION, p, 222, &
                              g_recv_data(1,p+1,k), g_recv_count(p+1), MPI_DOUBLE_PRECISION, p, 222, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                end do
       end if
    end if
end do

!if(myid.eq.1)print*,"recieved data", g_recv_data
!-----------------------------------------------
! 4. UNPACK: put received values into ghosts
!-----------------------------------------------
g_recv_pos = 0
do i = 1, pghosts
    p = neigh_proc(i)                ! which rank owns this ghost
    g_recv_pos(p+1) = g_recv_pos(p+1) + 1
    xcel(dest_ghost(i)) = g_recv_data( g_recv_pos(p+1),p+1,1)
    ycel(dest_ghost(i)) = g_recv_data( g_recv_pos(p+1),p+1,2)
    zcel(dest_ghost(i)) = g_recv_data( g_recv_pos(p+1),p+1,3)
    vol(dest_ghost(i)) = g_recv_data( g_recv_pos(p+1),p+1,4)

    !if(myid.eq.1.and.i.lt.5) print*," what we recevied in myid  ", myid,"from", p, dest_ghost(i), xcel(dest_ghost(i))

end do

g_send_pos   = 0
g_recv_pos   = 0
g_send_count = 0
g_recv_count = 0
deallocate(g_send_pos,g_recv_pos,g_send_data,g_recv_data,g_send_count,g_recv_count)

!do ng = 1,ntot
! print '(A,2I6,4F12.6)', "ng,xcel,ycel,zcel,vol ", ng, myid,xcel(ng), ycel(ng), zcel(ng), vol(ng)
!end do

!$acc update device(nc1(:), nc2(:), nc3(:), nc4(:), nparent(:), nghost(:), ntype(:), dl(:), &
!$acc               xcel(:), ycel(:), zcel(:), nside(:), sc1x(:), sc2x(:), sc3x(:), sc4x(:), &
!$acc               sc1y(:), sc2y(:), sc3y(:), sc4y(:), sc1z(:), sc2z(:), sc3z(:), sc4z(:), &
!$acc  vol(:), nodes, neles, nghosts, itest(:), x(:), y(:), z(:), nod(:,:),xmax,ymax,zmax,ymin,zmin, &
!$acc  pp_elem(:), neigh_proc(:), dest_ghost(:),dd(:))   
 
end subroutine geometry

