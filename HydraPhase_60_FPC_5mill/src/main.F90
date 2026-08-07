PROGRAM Main
  USE mpi
  USE GlobalVariables
  USE FluxModule
  USE ISO_C_BINDING
  USE run_tioga, ONLY : tioga_init_conn, tioga_solutions, tioga_fin, iblankcells, nc_t,iblanknodes,cellvals
  implicit none

  external Read_flow
  external geometry

!WRITE WALL_PRESSURE_TRACE NOW AFTER CALL PRIMITIVE, RATHER THAN CALL RK - WRITES PG AFTER INTERPOLATION (mod9)

!===========================
! Variable Declarations
!===========================
! Integers
INTEGER :: i, j, k, l, ie, r, s, v, rv, p, alloc_anu,irank, dummy, writevar,acn


! Reals (double precision)
REAL(DP) :: qgi, qpi, bbb, ccc, aaa, phitotal, phigmin, eg, ep, romix, pmix, functiong, psig
REAL(DP) :: diffu_global, diffv_global, diffw_global, diffp_global, tkgmax,tkgmin,tegmax,tegmin,amug 
REAL(DP) :: epsilonmax,epsilonmin

!Declarations for IDW reconstruction
REAL(DP) :: cx, cy, cz, dist, w
INTEGER  :: m, n, inode
LOGICAL, ALLOCATABLE :: has_nan(:), has_valid(:)
LOGICAL,SAVE :: trace_file_started = .false.

!Variable for Y_p
REAL(DP) :: yfgdt,yogdt,ypgdt,cvg

!Variable for N_p
REAL(DP) :: rpcutoff,volrpcutoff

! Allocatable arrays
REAL(DP), ALLOCATABLE :: cn2(:,:)
INTEGER,  ALLOCATABLE :: reqs(:)

! Character strings
CHARACTER(LEN=7)   :: file1, file3, string
CHARACTER(LEN=20)  :: file2, file4
CHARACTER(LEN=256) :: fname_grid, fname_bc, fname_procbc, fname_restart
CHARACTER(LEN=256) :: restart_dir, tecplot_dir, fname_tec, geometry_dir

character(len=1024) :: varList
character(len=256) :: tioga_input_dir

! Logical flags
LOGICAL :: dir_exists

!-------------------Variables to write wall pressure trace----------------(mod9)

LOGICAL, SAVE :: wall_trace_started = .false.
INTEGER, SAVE :: wall_cells(100000), n_wallcells
REAL(DP), SAVE :: wall_theta(100000)
INTEGER :: eid, ios2, wcc
REAL(DP) :: cell_radius, theta

!--------------------Cp calculation variables (mod9)------------------
REAL(DP), PARAMETER :: p_inf   = 101325.0d0
REAL(DP), PARAMETER :: rho_inf = 1.293218976d0
REAL(DP), PARAMETER :: v_inf   = 33.11969203d0
REAL(DP) :: cp

!=======================================
!Createing and cheking necessary folders
!=======================================
    
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
   
   call get_environment_variable('TIOGA_INPUT_DIR', tioga_input_dir) ! Getting tioga input directory

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
  
   open(unit=21,file='itter.in')!..............do later
   open(unit=23,file='Velocity.dat',access='append')!.........do later 
  
   write(fname_restart,"(A,I4.4,A)")'restart',myid,'.in'

   read(10,*)nodes,neles
   
   read(11, *)nghosts
      
   read(12,*) pghosts,max_pghosts ! number of partition ghost elements 
    
   close(10)
 !  close(11)
   close(12)
  
  ! Define folder names
  restart_dir = 'restart_files/'
  tecplot_dir = 'tecplot_files/'
  
  
  ! Rank 0 creates the directories
  if (myid == 0) then
   inquire(file=restart_dir, exist=dir_exists)
   if (.not. dir_exists) call system("mkdir -p " // trim(restart_dir))

   inquire(file=tecplot_dir, exist=dir_exists)
   if (.not. dir_exists) call system("mkdir -p " // trim(tecplot_dir))
endif
    
  ! Make sure all ranks wait
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  
NTOT = neles + nghosts + pghosts
NEPG = neles + nghosts + 1

NEQ  = 10
NVAR = 33

CALL Read_flow

CALL tioga_init_conn(trim(tioga_input_dir))

!==================================
! Base equations: phases
!==================================
#ifdef FIVE_PHASE
   ALLOCATE(uq(ntot), vq(ntot), wq(ntot), phiq(ntot))
   ALLOCATE(ur(ntot), vr(ntot), wr(ntot), phir(ntot))
   ALLOCATE(us(ntot), vs(ntot), ws(ntot), phis(ntot))
   neq       = 22
   alloc_anu = 5
#else
   ALLOCATE(uq(1), vq(1), wq(1), phiq(1))
   ALLOCATE(ur(1), vr(1), wr(1), phir(1))
   ALLOCATE(us(1), vs(1), ws(1), phis(1))
   neq       = 10
   alloc_anu = 3
#endif


!=====================================
! Allocate turbulence fields
!=====================================
#ifdef KE_TURB
   ALLOCATE(vistkg(neles), visteg(neles))
   ALLOCATE(tkg(ntot), teg(ntot))
   ALLOCATE(delTkgx(ntot), delTkgy(ntot), delTkgz(ntot))
   ALLOCATE(delTegx(ntot), delTegy(ntot), delTegz(ntot))
#else
   ALLOCATE(vistkg(1), visteg(1))
   ALLOCATE(tkg(1), teg(1))
   ALLOCATE(delTkgx(1), delTkgy(1), delTkgz(1))
   ALLOCATE(delTegx(1), delTegy(1), delTegz(1))
#endif

!====================
! Condensation model
!====================
#ifdef NP_P
   ALLOCATE(enp(ntot))
#else
   ALLOCATE(enp(1))
#endif

!=================
! Species transport
!=================
#ifdef YP_P
   ALLOCATE(yfg(ntot), yog(ntot), ypg(ntot))
   ALLOCATE(visyfg(neles), visyog(neles), visypg(neles))
   alloc_anu = 4
#else
   ALLOCATE(yfg(1), yog(1), ypg(1))
   ALLOCATE(visyfg(1), visyog(1), visypg(1))
#endif

!=====================================
! Turbulence model (k-epsilon)
!=====================================
#ifdef KE_TURB
   nkg = neq + 1
   neg = neq + 2
   neq = neq + 2

   nvar = nvar + 6
#endif

!====================
! Condensation equation
!====================
#ifdef NP_P
   npp = neq + 1
   neq = neq + 1
#endif

!====================
! Species equations
!====================
#ifdef YP_P
   nyfg = neq + 1
   nyog = neq + 2
   nypg = neq + 3
   neq  = neq + 3
#endif

!===================
!Allocating variable
!===================
  ALLOCATE(ANU(ntot,alloc_anu))
  ALLOCATE(x(nodes), y(nodes), z(nodes), xcel(ntot), ycel(ntot), zcel(ntot))
  ALLOCATE(nod(ntot,4),nc1(ntot),nc2(ntot),nc3(ntot))
  ALLOCATE(nc4(ntot),vol(ntot),dl(neles))
  ALLOCATE(sc1x(neles),sc1y(neles),sc1z(neles),sc2x(neles),sc2y(neles),sc2z(neles),sc3x(neles))
  ALLOCATE(sc3y(neles),sc3z(neles),sc4x(neles),sc4y(neles),sc4z(neles))
  ALLOCATE(itest(nodes),ug(ntot),vg(ntot),wg(ntot),pg(ntot))
  ALLOCATE(pp(ntot),tg(ntot),phig(ntot))
  ALLOCATE(up(ntot),vp(ntot),wp(ntot),phip(ntot))
  ALLOCATE(tp(ntot),egtotal(ntot),eptotal(ntot))
  ALLOCATE(rop(ntot),rog(ntot))
  ALLOCATE(ag(ntot),ap(ntot),cn(ntot,neq),rhs(neles,neq),d24(neles,neq))
  ALLOCATE(d2q(ntot,neq),cn1(ntot,neq),cu(ntot,neq))
  ALLOCATE(visug(neles),visvg(neles),viswg(neles),vistg(neles),phi(ntot))
  ALLOCATE(visup(neles),visvp(neles),viswp(neles),vistp(neles),delrogy(ntot),delrogx(ntot),delphigy(ntot))
  ALLOCATE(delphigx(ntot),delpgz(ntot),delvgy(ntot),delvgx(ntot),deltgy(ntot),deltgx(ntot),delphigz(ntot))
  ALLOCATE(delwgz(ntot),delugz(ntot),delvgz(ntot),delTgz(ntot),delwgx(ntot),delugx(ntot),delrogz(ntot),delwgy(ntot))
  ALLOCATE(delugy(ntot),delpgy(ntot),delupz(ntot),delppz(ntot),delppx(ntot),delppy(ntot),delwpx(ntot),delwpy(ntot))
  ALLOCATE(delupx(ntot),delupy(ntot),delvpy(ntot),delvpx(ntot),deltpy(ntot),deltpx(ntot),delPgx(ntot),delvpz(ntot))
  ALLOCATE(delTpz(ntot),delropy(ntot),delropx(ntot),delwpz(ntot),delropz(ntot))
  ALLOCATE(nparent(nghosts), nghost(nghosts), ntype(nghosts), nside(nghosts))
  ALLOCATE(cn2(ntot,neq))
  ALLOCATE(pp_elem(pghosts),global_pp(pghosts), neigh_proc(pghosts), neigh_side(pghosts))
  ALLOCATE(dest_ghost(pghosts),loc_pghost(pghosts), neigh_glob(pghosts))
  ALLOCATE(cu_send_pos(nprocs), cu_recv_pos(nprocs))
  ALLOCATE(dd(ntot))
!============================================ 
!Overset grid identity and Node+Cell q vals
!============================================
  ALLOCATE(btag(ntot))
  ALLOCATE(iblank(ntot))
  ALLOCATE(qcell(neles*neq))
  ALLOCATE(q_node(nodes*neq))
  ALLOCATE(wsum(nodes))
!=============================================
!Pre caclulating the boundary cell information
!=============================================   
    do i = 1, nghosts
      read(11, *) nparent(i), nghost(i), ntype(i), nside(i)
    end do
    close(11)

    n_bc2  = 0; n_bc3  = 0; n_bc5  = 0
    n_bc51 = 0; n_bc53 = 0

do i = 1, nghosts
  select case (ntype(i))
  case (2);  n_bc2  = n_bc2  + 1
  case (3);  n_bc3  = n_bc3  + 1
  case (5);  n_bc5  = n_bc5  + 1
  case (51); n_bc51 = n_bc51 + 1
  case (53); n_bc53 = n_bc53 + 1
  end select
end do  

ALLOCATE(bc2_list(n_bc2), bc3_list(n_bc3), bc5_list(n_bc5), bc51_list(n_bc51), bc53_list(n_bc53))

n_bc2  = 0; n_bc3  = 0; n_bc5  = 0
n_bc51 = 0; n_bc53 = 0

do i = 1, nghosts
  select case (ntype(i))
  case (2)
    n_bc2 = n_bc2 + 1
    bc2_list(n_bc2) = i
  case (3)
    n_bc3 = n_bc3 + 1
    bc3_list(n_bc3) = i
  case (5)
    n_bc5 = n_bc5 + 1
    bc5_list(n_bc5) = i
  case (51)
    n_bc51 = n_bc51 + 1
    bc51_list(n_bc51) = i
  case (53)
    n_bc53 = n_bc53 + 1
    bc53_list(n_bc53) = i
  end select
end do



  
!=========================
!Createing varibles in GPU
!=========================   
 
!$acc data create( &
!$acc ag(:),ap(:),anu(:,:),as, &
!$acc cn(:,:),cpp,co2,co4,c23,cmug,cpg,cfl,cn1(:,:),cu, &
!$acc dl(:),dtmin,diap,d24(:,:),d2q(:,:),dt, &
!$acc delphigx(:),delropy(:),delrogz(:),delphigy(:), &
!$acc delrogx(:),delppy(:),delphigz(:),delpgx(:),deltpy(:),deltgz(:), &
!$acc delwpx(:),delvpy(:),delvgz(:),delupx(:),deltgy(:),delropz(:), &
!$acc delwgx(:),delvgy(:),delupz(:),delugx(:),delwpz(:), &
!$acc deltpz(:),delugy(:),delvgx(:),delvpz(:),delwgy(:),deltgx(:), &
!$acc delugz(:),delupy(:),delvpx(:),delwgz(:),delwpy(:),deltpx(:), &
!$acc delpgy(:),delppx(:),delppz(:),delrogy(:),delropx(:),delpgz(:), &
!$acc diffp,diffq,diffr,diffu,diffv,diffw,diffx,diffy,diffz,difft, &
!$acc epslnmax,epslnmin,total, &
!$acc foursigmabyd, &
!$acc gamag,gamap,gravity, &
!$acc itest(:),dd(:), &
!$acc myorder,egtotal,eptotal, &
!$acc neles,nod(:,:),nc1(:),nc2(:),nc3(:), &
!$acc nc4(:),nparent(:),nghost(:),ntype(:),neq, & 
!$acc nodes,neles,nghosts,nside(:), &
!$acc phig(:),phip(:),pg(:),pp(:),pinf,pmax,pmin,prg,phi(:),phipzero,prp, &
!$acc rhs(:,:),rog(:),rrg,rop(:), &
!$acc sc1x(:),sc2x(:),sc3x(:),sc1y(:), & 
!$acc sc2y(:),sc3y(:),sc4y(:),sc4x(:), &
!$acc sc3z(:),sc2z(:),sc1z(:),sc4z(:), &
!$acc tg(:),tp(:),tgmax,tpmax,tgmin,tpmin,time, ts,&
!$acc ug(:),up(:),ugmin,vgmin,wgmin, &
!$acc vp(:),vg(:),vol(:),vistp(:),visup(:),visvp(:),visug(:),Vel_mag, &
!$acc visvg(:),viswp(:),viswg(:),vistg(:), &
!$acc wp(:),wg(:), & 
!$acc x(:),xcel(:),xmax,ymax,zmax,ymin,zmin, &
!$acc y(:),ycel(:), &
!$acc z(:),zcel(:),nepg,ntot,nprocs,pghosts,pp_elem(:), neigh_proc(:), dest_ghost(:),&
!$acc uq(:),vq(:),wq(:),ur(:),vr(:),wr(:),us(:),vs(:),ws(:),phiq(:),phir(:),phis(:),& 
!$acc rhoq, rhor, rhos, radq, radr, rads, volq, volr, vols, qmass, rmass, smass, &
!$acc vistkg(:),visteg(:),tkg(:),teg(:),delTkgx(:),delTkgy(:),delTkgz(:),delTegx(:),delTegy(:),delTegz(:), &
!$acc TOLKG,TOLEG,nkg,neg,istart, &
!$acc max_pghosts,nvar, &
!$acc npp,enpcutoff,radp,volrp,enp(:), &
!$acc nyfg,nyog,nypg,tolsp,cpfg,cpog,cppg,amolfg,amolog,amolpg,ypg_inf,yfg_inf,yog_inf, &
!$acc yfg(:),yog(:),ypg(:),visyfg(:),visyog(:),visypg(:),amolg, &
!$acc cp_t,cp_o,cp_p,amol_f,amol_o,amol_p,sratio,hfp,preexp,ebyr, &
!$acc bc2_list(:), bc3_list(:), bc5_list(:), bc51_list(:), bc53_list(:), &
!$acc n_bc2, n_bc3, n_bc5, n_bc51, n_bc53)

!$acc update device(nepg,ntot,nprocs,pghosts,neq)
           
!============================
! Initialize global variables
!============================
!make these constant variables capital for distinguish 
  
  CO2 = 1.00 / 16.00
  CO4 = 1.00 / 1024.00
  C23 = 2.00 / 3.00
  CMUG = 11.848e-8 * DSQRT(amolg)
  RRG = 287.00
  GRAVITY = 0.00
! amup = 0.0010
  PRP = 10.00
!  akp = amup * cpp / prp
  FOURSIGMABYD = 4.00 * surfacetension / diap
!  PHIPZERO = epslnmin
  ITTER = 0
  TIME = 0.00
  DTMIN = 1.0E8
  
#ifdef FIVE_PHASE
  volq=(4.0d0/3.0d0)*3.1415926d0*radq**3
  volr=(4.0d0/3.0d0)*3.1415926d0*radr**3
  vols=(4.0d0/3.0d0)*3.1415926d0*rads**3
  qmass=rhoq*volq
  rmass=rhor*volr
  smass=rhos*vols
#endif  

#ifdef NP_P
  radp=1.0e-3
  rpcutoff=1.0e-6
  volrp=(4.0d0/3.0d0)*3.1415926d0*radp**3
  volrpcutoff=(4.0d0/3.0d0)*3.1415926d0*rpcutoff**3
  enpcutoff=epslnmin/volrpcutoff
#endif  

#ifdef YP_P
  tolsp=1.0e-6
  cpfg=2000.0d0
  cpog=1460.88d0
  cppg=1275.0d0
  amolfg=54.0d0
  amolog=27.8d0
  amolpg=26.8d0
  ypg_inf=tolsp
  yfg_inf=tolsp
  yog_inf=1.0d0-tolsp
  cpg=yfg_inf*cpfg+yog_inf*cpog+ypg_inf*cppg
  rrg=8314.0d0*(yfg_inf/amolfg+yog_inf/amolog+ypg_inf/amolpg)
  amolg=8314.0d0/rrg
  cvg=cpg-rrg
  gamag=cpg/cvg
  CMUG = 11.848e-8 * DSQRT(amolg)
#endif
  
!==============================  
! Call Geometry Subroutine once
!==============================  
CALL geometry

  rogmax = pmax / (rrg * tgmax)
  ropmax = (pmax + pinf) * gamap / ((gamap - 1.0D0) * cpp * tpmax)
  
  egtotalmax = rrg * tgmax / (gamag - 1.0D0) + 0.5D0 * (ugmax**2 + vgmax**2 + wgmax**2)
  eptotalmax = (cpp * tpmax) / gamap + pinf / ropmax + 0.5D0 * (upmax**2 + vpmax**2 + wpmax**2)
  
  rogmin = pmin / (rrg * tgmin)
  ropmin = (pmin + pinf) * gamap / ((gamap - 1.0D0) * cpp * tpmin)
  
!  print*,ropmin,pmin,pinf,gamap,gamap,cpp,tpmin
  egtotalmin = rrg * tgmin / (gamag - 1.0D0) + 0.5D0 * (ugmin**2 + vgmin**2 + wgmin**2)
  eptotalmin = (cpp * tpmin) / gamap + pinf / ropmin + 0.5D0 * (upmin**2 + vpmin**2 + wpmin**2)

  phigmin = epslninf
  phigmax = 1.0D0 - epslninf

  amug =  as*dsqrt(Tgmin)/(1.0d0+ ts/tgmin) 

#ifdef KE_TURB
    epsilonmax = 1.0D-10
    epsilonmin = 1.0D-12

    TKGMAX=max(0.666D0*(0.1D0*UGMAX)**2.0D0,epsilonmax)
    TKGMIN=max(0.666D0*(0.1D0*UGMIN)**2.0D0,epsilonmin)
    
    TEGMAX=0.0845D0*ROGMAX*TKGMAX*TKGMAX/(AMUG*0.1D0)
    TEGMIN=0.0845D0*ROGMIN*TKGMIN*TKGMIN/(AMUG*0.1D0)
	    
    TOLKG=DMIN1(TKGMAX,TKGMIN)*0.001D0  ! Gloabal
    TOLEG=DMIN1(TEGMAX,TEGMIN)*0.001D0 !Global
#endif 


!================================
! Read restart data if istart > 0
!================================  
if(istart.gt.0)then
       
        write(fname_restart,"(A,I4.4,A)") trim(restart_dir)//'restart', myid, '.in'
        open(unit=13, file=fname_restart, status='unknown', action='read')
        read(13,*)time
        read(13,*)((cu(i,k),i=1,neles),k=1,neq)
        read(13,*)(ug(nghost(i)),i=1,nghosts)
        read(13,*)(vg(nghost(i)),i=1,nghosts)
        read(13,*)(wg(nghost(i)),i=1,nghosts)
        read(13,*)(pg(nghost(i)),i=1,nghosts)
        read(13,*)(tg(nghost(i)),i=1,nghosts)
        read(13,*)(rog(nghost(i)),i=1,nghosts)
        read(21,*)itter
        close(13)
        close(21)
else   

!==================
!Initialize domain
!==================    
CALL INIT_COND   
    
endif

   Total = 0
   vel_mag = 0.0d0        
call mpiexcu

  ! Initialize viscosity variables
  DO i = 1, neles
   visug(i) = 0.0D0
   visvg(i) = 0.0D0
   viswg(i) = 0.0D0
   vistg(i) = 0.0D0
   visup(i) = 0.0D0
   visvp(i) = 0.0D0
   viswp(i) = 0.0D0
   vistp(i) = 0.0D0
#ifdef KE_TURB
   vistkg(i) = 0.0D0
   visteg(i) = 0.0D0
#endif  
#ifdef YP_P
   visyfg(i)= 0.0D0
   visyog(i)= 0.0D0
   visypg(i)= 0.0D0
#endif  
  END DO

 d24 = 0.0D0
   
 DO i = 1, ntot
   DO k = 1, neq
    cn(i,k) = cu(i,k)
   END DO
  END DO
  
!=============================
!Update global variable in GPU
!=============================
    
!$acc update device(cu(:,:),cn(:,:),d24(:,:),gravity,cpp,foursigmabyd,diap,pmax, &
!$acc rrg,tgmax,pinf,gamap,tpmax,dtmin,pmin,tgmin,tpmin,gamag,co2,co4,c23,cmug,  &
!$acc time,cpg,prg,epslnmax,epslnmin,cfl,myorder,ag(:),ap(:),phi(:),ug(:),dt,cn1(:,:), &
!$acc ugmin,vgmin,wgmin,as,total,ts,vg(:),pg(:),anu(:,:),tp(:),up(:),wp(:),wg(:),tg(:), &
!$acc rog(:),phig(:),d2q(:,:),vp(:),phip(:),pp(:),rop(:),phipzero,uq(:),vq(:),wq(:),ur(:), &
!$acc vr(:),wr(:),us(:),vs(:),ws(:),phiq(:),phir(:),phis(:),rhoq, rhor, rhos, radq, radr,rads,Vel_mag, &
!$acc volq, volr, vols, qmass, rmass, smass,prp,vistkg(:),visteg(:),tkg(:),teg(:),delTkgx(:), &
!$acc delTkgy(:),delTkgz(:),delTegx(:),delTegy(:),delTegz(:),TOLKG,TOLEG,nkg,neg,istart, &
!$acc max_pghosts,nvar,npp,enpcutoff,radp,volrp,enp(:),nyfg,nyog,nypg,tolsp,cpfg,cpog,cppg, &
!$acc amolfg,amolog,amolpg,ypg_inf,yfg_inf,yog_inf,yfg(:),yog(:),ypg(:),visyfg(:),visyog(:),visypg(:),&
!$acc bc2_list(:), bc3_list(:), bc5_list(:), bc51_list(:), bc53_list(:), &
!$acc n_bc2, n_bc3, n_bc5, n_bc51, n_bc53)

!!!$acc cnd_send_pos,cnd_recv_pos,cnd_send_count,cnd_recv_count,cnd_send_data,cnd_recv_data, &
!!!$acc del_send_pos,del_recv_pos,del_send_count,del_recv_count,del_send_data,del_recv_data, &
!========================================
!           Main iteration loop
!========================================
OPEN(unit=99, file='checking_nodes/active_cells.dat', status='replace', action='write') !-----writing active cells (mod3)
do i = 1, neles
    if(iblank(i) == 1) WRITE(99, *) i
end do
CLOSE(99)
  DO iter = 1, niter
   DO iters = 1, itersub
    DO iterss = 1,1
!print*,myid," Time loop start"
    CALL primitive
   
    call mpiexcn
    
   
!$acc serial present(total) 
        dtmin=1.0E8
        total = total + 1
!$acc end serial
     
!$acc parallel loop present(ap(:),rop(:),up(:),vp(:),pp(:),dtmin,pg(:),rog(:),ag(:),wp(:),wg(:),vg(:), &
!$acc ug(:),dl(:),gamag,gamap,pinf,neles,nyfg,nyog,nypg,cpg,cpfg,cpog,cppg,amolfg,amolog,amolpg,rrg,amolg) &
!$acc reduction(min:dtmin) &
!$acc private(qpi,qgi,i,yfgdt,yogdt,ypgdt,cvg)      
    DO i = 1, neles
    if(iblank(i)/=1)cycle
#ifdef YP_P
     yfgdt=cu(i,nyfg)/cu(i,1)
     yogdt=cu(i,nyog)/cu(i,1)
     ypgdt=cu(i,nypg)/cu(i,1)
     cpg=yfgdt*cpfg+yogdt*cpog+ypgdt*cppg
     rrg=8314.0d0*(yfgdt/amolfg+yogdt/amolog+ypgdt/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
!if(myid==1.and.i.lt.10)     print*,"gamag",yfgdt,yogdt,ypgdt
#endif	
     ag(i) = SQRT((gamag * pg(i)) / rog(i))
     ap(i) = SQRT((gamap * (pp(i) + pinf)) / rop(i))
     QGI = SQRT(ug(i)**2 + vg(i)**2 + wg(i)**2)
     QPI = SQRT(up(i)**2 + vp(i)**2 + wp(i)**2)
     dtmin = MIN(dtmin, dl(i) / (MAX(ag(i), ap(i)) + MAX(QGI, QPI)))
    END DO

!$acc end parallel loop 
!!$acc update host(ag,ap,dtmin)
!$acc update host(dtmin)

     call MPI_Allreduce(dtmin, dtmin_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
     
     ! Update timestep
     dt   = dtmin_global * cfl
     time = time + dt

!$acc update device(time,dt)

!======================
!Call smooth subroutine
!======================
 CALL SMOOTH

#ifdef NP_P
    CALL mpiexcn
#endif


#ifdef FIVE_PHASE
 CALL DISSIPATIONQRS
#endif


!==================
!Call RK subroutine
!==================
 CALL RK

!--------------------Writing wall pressure trace (mod9)--------------

!if (.not. wall_trace_started) then
!  n_wallcells = 0
!  open(unit=98,file='checking_nodes/wall_radius.dat',status='old')
!  do
!    read(98,*,iostat=ios2) eid, cell_radius, theta
!    if (ios2/=0) exit
!    n_wallcells = n_wallcells + 1
!    wall_cells(n_wallcells) = eid+262314   ! offset must be added if combined grid used (offset= no.bgcells)
!    wall_theta(n_wallcells) = theta
!  end do
!  write(*,*)"wallcells: ",n_wallcells
!  close(98)
!  open(unit=99,file='checking_nodes/wall_pressure_trace.dat',status='replace')
!  write(99,'(A)') '# iter  time  cell_id  pg  cp  theta'
!  wall_trace_started = .true.
!else
!  open(unit=99,file='checking_nodes/wall_pressure_trace.dat',position='append')
!end if
!
!do wcc = 1, n_wallcells
!  cp = (pg(wall_cells(wcc)) - p_inf) / (0.5d0 * rho_inf * v_inf**2)     ! Calculate and write cp (mod9)
!  write(99,'(I8,1X,F20.16,1X,I8,1X,F20.10,1X,F20.16,1X,F20.16)') &
!       iter, time, wall_cells(wcc), pg(wall_cells(wcc)), cp , wall_theta(wcc)
!end do
!close(99)

!========================================================= Calculate inverse-distance weighted node q values (mod4)
!=========================================================
acn=0      !Active node counter for test (mod8)
if(.not.allocated(has_nan)) allocate(has_nan(nodes))
if(.not.allocated(has_valid)) allocate(has_valid(nodes))
has_nan = .false.
has_valid = .false.
q_node = 0.0
wsum   = 0.0
DO ie = 1, neles
  if (iblank(ie) /= 1) CYCLE
  cx = SUM(x(nod(ie,1:4))) / 4.0
  cy = SUM(y(nod(ie,1:4))) / 4.0
  cz = SUM(z(nod(ie,1:4))) / 4.0
  DO m = 1, 4
    inode = nod(ie,m)
    if (iblanknodes(inode) /= 1) CYCLE
    acn=acn+1
    dist = SQRT((x(inode)-cx)**2 + (y(inode)-cy)**2 + (z(inode)-cz)**2)
    w = 1.0 / dist
    IF (cu(ie,1) /= cu(ie,1)) THEN
      has_nan(inode) = .true.
    ELSE
      has_valid(inode) = .true.
      DO k = 1, neq
        q_node((inode-1)*neq+k) = q_node((inode-1)*neq+k) + w*cu(ie,k)
      END DO
      wsum(inode) = wsum(inode) + w
    END IF
  END DO
END DO
write(*,*)"ACN= ",acn
open(unit=78,file='nan_field_cells.dat',status='replace')  !(mod7)
do ie=1,neles
  if (iblank(ie)==1 .and. cu(ie,1)/=cu(ie,1)) write(78,*) ie
end do
close(78)

open(unit=69,file='checking_nodes/wsum.dat',status='replace')
open(unit=70,file='checking_nodes/cu.dat',status='replace')
open(unit=71,file='checking_nodes/q_node.dat',status='replace')

m = 0
do n = 1, nodes
  if (wsum(n) == 0.0) m = m + 1
  write(69,*)wsum(n)
enddo

! check cu for NaN
m = 0
do ie = 1, neles
  if (isnan(cu(ie,1))) m = m + 1
  write(70,*) (cu(ie,k), k=1,neq)
enddo
write(*,*) 'writing cu,wsum,qnode'

! check q_node for NaN
m = 0
do n = 1, nodes
  if (isnan(q_node((n-1)*neq+1))) m = m + 1
 write(71,*) n, (q_node((n-1)*neq+k), k=1,neq)
! write(71,*)q_node((n-1)*neq+1)
enddo
close(69)
close(70)
close(71)

CALL tioga_solutions(q_node,neq,nodes)

do i = 1, neles                                               !copying only fringe cell info into cu (after tioga interp)(mod8)
  if (iblank(i) == -1) cu(i,1:neq) = cellvals(i,1:neq)
end do

open(unit=96,file='checking_nodes/cu2.dat',status='replace')  !writing cu after tioga interp. to check (mod8)
do i = 1, neles
  write(96,*) (cu(i,k),k=1,neq)
end do
close(96)

!$acc serial        
        diffp=0.0d0
        diffu=0.0d0
        diffv=0.0d0
        diffw=0.0d0
!$acc end serial
       
       diffu_global=0.0d0
       diffv_global=0.0d0
       diffw_global=0.0d0
       diffp_global=0.0d0


!=========================
!Call primitive subroutine
!=========================
CALL PRIMITIVE

!--------------------Writing wall pressure trace (mod9)--------------

if (.not. wall_trace_started) then
  n_wallcells = 0
  open(unit=98,file='checking_nodes/wall_radius.dat',status='old')
  do
    read(98,*,iostat=ios2) eid, cell_radius, theta
    if (ios2/=0) exit
    n_wallcells = n_wallcells + 1
    wall_cells(n_wallcells) = eid+262314   ! offset must be added if combined grid used (offset= no.bgcells)
    wall_theta(n_wallcells) = theta
  end do
  write(*,*)"wallcells: ",n_wallcells
  close(98)
  open(unit=99,file='checking_nodes/wall_pressure_trace.dat',status='replace')
  write(99,'(A)') '# iter  time  cell_id  pg  cp  theta'
  wall_trace_started = .true.
else
  open(unit=99,file='checking_nodes/wall_pressure_trace.dat',position='append')
end if

do wcc = 1, n_wallcells
  cp = (pg(wall_cells(wcc)) - p_inf) / (0.5d0 * rho_inf * v_inf**2)     ! Calculate and write cp (mod9)
  write(99,'(I8,1X,F20.16,1X,I8,1X,F20.10,1X,F20.16,1X,F20.16)') &
       iter, time, wall_cells(wcc), pg(wall_cells(wcc)), cp , wall_theta(wcc)
end do
close(99)

!================
!Call mpiexchange 
!================
 CALL MPIEXCN

! Update differences
!$acc parallel loop private(i) &
!$acc present(cu(:,:),cn(:,:),neles,diffp,diffu,diffv,diffw) 
do i = 1, neles
    diffp = max(diffp, abs(CN(i, 1) - Cu(i, 1)))
    diffu = max(diffu, abs(CN(i, 2) - Cu(i, 2)))
    diffv = max(diffv, abs(CN(i, 3) - Cu(i, 3)))
    diffw = max(diffw, abs(CN(i, 4) - Cu(i, 4)))
!    print*,"diffp",diffp
end do
!$acc end parallel loop 

!================
!Call mpiexchange 
!================
    call mpiexcn

!$acc parallel loop present(cu(:,:),cn(:,:),ntot,neq)private(i,j) 
  do i = 1,ntot
    do j = 1,neq
    Cu(i, j) = CN(i, j)
    end do 
   end do 
!$acc end parallel loop 


!if(total==2) then 
!do i = 1,neles
!if(myid==1) print*,"cn11", i, cn(i,11)
!enddo
!STOP 
!endif
   
!print*, "diffu",myid,diffu
 
!   write(23,72)time,sqrt(ug(6513446)**2+vg(6513446)**2+wg(6513446)**2)
!72 format(5f16.8)
 
!$acc serial present(ug,vg,wg,vel_mag)
vel_mag = sqrt(ug(1)**2+vg(1)**2+wg(1)**2)  !(mod2) 896880 -> 1 
!$acc end serial

 END DO 
 
!$acc update host(time,diffu,diffv,diffw,diffp,vel_mag)
call MPI_Allreduce(diffu, diffu_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(diffv, diffv_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(diffw, diffw_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_Allreduce(diffp, diffp_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)


if(myid==0)PRINT '(I6, 1PE20.10, 5F24.16)', iters, time, diffu_global, diffv_global, diffw_global, diffp_global
if(myid==0) write(23,*)time,vel_mag  

!1742938
END DO   !-----------Sub-iteration loop over

!$acc update host(cu(:,:),pg(:),ug(:),vg(:),wg(:),tg(:),rog(:))

!     file1='restart'
!     itter=itter+1
!        write(unit=string,fmt=' (i7.7) ')itter
!        file2=file1//string//'.in'
!        open(unit=41,file=file2)
!	rewind 41
!	write(41,*)time
!        write(41,*)((cu(i,k),i=1,neles),k=1,10)
!        write(41,*)(ug(nghost(i)),i=1,nghosts)
!        write(41,*)(vg(nghost(i)),i=1,nghosts)
!        write(41,*)(wg(nghost(i)),i=1,nghosts)
!        write(41,*)(pg(nghost(i)),i=1,nghosts)
!        write(41,*)(tg(nghost(i)),i=1,nghosts)
!        write(41,*)(rog(nghost(i)),i=1,nghosts)
!        close(41)
                   
!        do i = 1, nodes
!        itest(i)=0
!        cn2(i,1)=0.0d0
!        cn2(i,2)=0.0d0
!        cn2(i,3)=0.0d0
!        cn2(i,4)=0.0d0
!        cn2(i,5)=0.0d0
!        cn2(i,6)=0.0d0
!        cn2(i,7)=0.0d0
!        cn2(i,8)=0.0d0
!        cn2(i,9)=0.0d0
!        cn2(i,10)=0.0d0
!        enddo

!        do i = 1, neles
!        do j = 1, 4
!        k=nod(i,j)
!        itest(k)=itest(k)+1
!        do l = 1, 10
!        cn2(k,l)=cn2(k,l)+cu(i,l)
!        enddo
!        enddo
!        enddo

!        do i = 1, nodes
!        do j = 1, 10
!        cn2(i,j) = cn2(i,j) / DBLE(itest(i))
!        enddo
!        enddo



!        file3='tecplot'
!        write(unit=string,fmt=' (i7.7) ')itter
!        file4=file3//string//'.dat'
!        open(unit=43,file=file4)
!        write(43,*) 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "PG", "PP", "T", "UP", "VP", "WP", "PHIG", "TP", "ROG", "ROP"'
!        write(43,*) 'ZONE F=FEPOINT, ET=TETRAHEDRON, N =', nodes, ', E =', neles
!        write(43,*)'SOLUTIONTIME=',time


        
!        do i = 1,nodes

!        ug(i) = cn2(i, 2) / cn2(i, 1)
!        vg(i) = cn2(i, 3) / cn2(i, 1)
!        wg(i) = cn2(i, 4) / cn2(i, 1)
!        up(i) = cn2(i, 7) / cn2(i, 6)
!        vp(i) = cn2(i, 8) / cn2(i, 6)
!        wp(i) = cn2(i, 9) / cn2(i, 6)
!   eg = cn2(i, 5) / cn2(i, 1) - 0.5d0 * (ug(i)**2 + vg(i)**2 + wg(i)**2)
!   ep = cn2(i, 10) / cn2(i, 6) - 0.5d0 * (up(i)**2 + vp(i)**2 + wp(i)**2)
!   bbb = (foursigmabyd + pinf * gamap) - (cn2(i, 1) * (gamag - 1.0d0) * eg) - (cn2(i, 6) * (gamap - 1.0d0) * ep)
!   ccc = (cn2(i, 1) * (gamag - 1.0d0) * eg * (foursigmabyd + pinf * gamap))
!   pg(i) = (-bbb + dsqrt(bbb**2 + 4*ccc)) * 0.5d0
!   pp(i) = pg(i) + foursigmabyd

!   phig(i) = (cn2(i, 1) * (gamag - 1.0d0) * eg) / pg(i)
!   phip(i) = 1.0d0 - phig(i)

!   rog(i) = cn2(i, 1) / phig(i)
!   rop(i) = cn2(i, 6) / phip(i)
!   tg(i) = pg(i) / (rog(i) * rrg)
!   tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0d0) * rop(i) * cpp)
!------------------------------------------------------------------- 

!   IF (phig(i) < epslnmin) THEN
!    phig(i) = epslnmin
!    phip(i) = 1.0d0 - phig(i)
!    ug(i) = up(i)
!    vg(i) = vp(i)
!    wg(i) = wp(i)
!    tg(i) = tp(i)
!    pg(i) = pp(i)
!   END IF

!   IF (phip(i) < epslnmin) THEN
!    phip(i) = epslnmin
!    phig(i) = 1.0d0 - phip(i)
!    up(i) = ug(i)
!    vp(i) = vg(i)
!    wp(i) = wg(i)
!    tp(i) = tg(i)
!    pp(i) = pg(i)
!   END IF

!   IF (phig(i) >= epslnmin .AND. phig(i) <= epslnmax) THEN
!    psig = (phig(i) - epslnmin) / (epslnmax - epslnmin)
!    functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
!    ug(i) = functiong * ug(i) + (1.0d0 - functiong) * up(i)
!    vg(i) = functiong * vg(i) + (1.0d0 - functiong) * vp(i)
!    wg(i) = functiong * wg(i) + (1.0d0 - functiong) * wp(i)
!    tg(i) = functiong * tg(i) + (1.0d0 - functiong) * tp(i)
!    pg(i) = functiong * pg(i) + (1.0d0 - functiong) * pp(i)
!   END IF

!   IF (phip(i) >= epslnmin .AND. phip(i) <= epslnmax) THEN
!    psig = (phip(i) - epslnmin) / (epslnmax - epslnmin)
!    functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
!    up(i) = functiong * up(i) + (1.0d0 - functiong) * ug(i)
!    vp(i) = functiong * vp(i) + (1.0d0 - functiong) * vg(i)
!    wp(i) = functiong * wp(i) + (1.0d0 - functiong) * wg(i)
!    tp(i) = functiong * tp(i) + (1.0d0 - functiong) * tg(i)
!    pp(i) = functiong * pp(i) + (1.0d0 - functiong) * pg(i)
!   END IF

    ! Compute densities and mixed density
!   rog(i) = pg(i) / (rrg * tg(i))
!   rop(i) = (pp(i) + pinf) * gamap / ((gamap - 1.0d0) * cpp * tp(i))
   !rom(i) = phig(i) * rog(i) + phip(i) * rop(i)

    ! Write output
!   WRITE(43, *) x(i), y(i), z(i), ug(i), vg(i), wg(i), pg(i), &
!               pp(i), tg(i), up(i), vp(i), wp(i), phig(i), &
!                 tp(i), rog(i), rop(i)
!  END DO

  ! Loop over elements for nodal data
!  DO ie = 1, neles
!   WRITE(43, *) nod(ie, 1), nod(ie, 2), nod(ie, 3), nod(ie, 4)
!   if(myid==2)print*, ie, nod(ie,1), nod(ie, 2), nod(ie, 3), nod(ie, 4)
!  END DO   
 
                
! close(43)
!       ----------------------  
   
  end do
!=======================
!END calcualtions in GPU 
!======================= 
!$acc end data 

!==================================================
!Writing restart file seperatly for each processors
!==================================================
  rewind 13
  write(fname_restart,"(A,I4.4,A)") trim(restart_dir)//'restart', myid, '.in'
  open(unit=13, file=fname_restart, status='unknown', action='write')
  !write(fname_restart, "(A, I4.4,A)")'restart', myid,'.in'
  !open(unit=13, file=fname_restart, status='unknown',action='write')
  write(13,*)time
  write(13,*)((cu(i,k),i=1,neles),k=1,neq)
  write(13,*)(ug(nghost(i)),i=1,nghosts)
  write(13,*)(vg(nghost(i)),i=1,nghosts)
  write(13,*)(wg(nghost(i)),i=1,nghosts)
  write(13,*)(pg(nghost(i)),i=1,nghosts)
  write(13,*)(tg(nghost(i)),i=1,nghosts)
  write(13,*)(rog(nghost(i)),i=1,nghosts)
  close(13)
  rewind 21
  write(21,*) itter
  close(21)


! Process each element and write data for Tecplot
  do i = 1, neles
   ug(i) = cu(i, 2) / cu(i, 1)
   vg(i) = cu(i, 3) / cu(i, 1)
   wg(i) = cu(i, 4) / cu(i, 1)
   up(i) = cu(i, 7) / cu(i, 6)
   vp(i) = cu(i, 8) / cu(i, 6)
   wp(i) = cu(i, 9) / cu(i, 6)
   eg = cu(i, 5) / cu(i, 1) - 0.5d0 * (ug(i)**2 + vg(i)**2 + wg(i)**2)
   ep = cu(i, 10) / cu(i, 6) - 0.5d0 * (up(i)**2 + vp(i)**2 + wp(i)**2)
   
#ifdef KE_TURB
   tkg(i) = cu(i, nkg) / cu(i, 1)
   teg(i) = cu(i, neg) / cu(i, 1)
#endif      

#ifdef NP_P
    enp(i) = cu(i,npp)
#endif

#ifdef YP_P
      yfg(i) = cu(i, nyfg) / cu(i, 1)
      yog(i) = cu(i, nyog) / cu(i, 1)
      ypg(i) = cu(i, nypg) / cu(i, 1)
	
     cpg=yfg(i)*cpfg+yog(i)*cpog+ypg(i)*cppg
     rrg=8314.0d0*(yfg(i)/amolfg+yog(i)/amolog+ypg(i)/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
#endif  
   
#ifdef FIVE_PHASE
    phiq(i)=dmin1(dmax1(cu(i,11)/rhoq,epslnmin),1.0d0-epslnmin)
    phir(i)=dmin1(dmax1(cu(i,15)/rhor,epslnmin),1.0d0-epslnmin)
    phis(i)=dmin1(dmax1(cu(i,19)/rhos,epslnmin),1.0d0-epslnmin)
    	
   aaa=dmin1(dmax1(1.0d0-phiq(i)-phir(i)-phis(i),epslnmin),1.0d0-epslnmin) 	
   bbb=(aaa*pinf*gamap)-(cu(i,1)*(gamag-1.0d0)*eg)-(cu(i,6)*(gamap-1.0d0)*ep)
   ccc=(cu(i,1)*(gamag-1.0d0)*eg*(gamap*pinf))
   pg(i)=(-bbb+dsqrt(bbb**2+4.0d0*aaa*ccc))/(2.0d0*aaa)
   pp(i) = pg(i) + foursigmabyd
 
   ! Compute densities and temperatures      
   rog(i) = pg(i) / ((gamag - 1.0d0)*eg)
   rop(i) = (pp(i) + pinf*gamap)  / ((gamap - 1.0d0)*ep)
   
   tg(i) = pg(i) / (rog(i) * rrg)
   tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0_dp) * rop(i) * cpp)
  
   ! Compute phase fractions
   phig(i) = dmin1(dmax1(cu(i,1)/rog(i),epslnmin),(1-epslnmin))
   phip(i) = dmin1(dmax1(cu(i,6)/rop(i),epslnmin),(1-epslnmin)) 
   
   !Jameson
   phitotal=phig(i)+phip(i)+phiq(i)+phir(i)+phis(i)
   phig(i)=phig(i)/phitotal
   phip(i)=phip(i)/phitotal
   phiq(i)=phiq(i)/phitotal
   phir(i)=phir(i)/phitotal
   phis(i)=phis(i)/phitotal
   
    uq(i)=cu(i,12)/cu(i,11)
    vq(i)=cu(i,13)/cu(i,11)
    wq(i)=cu(i,14)/cu(i,11)
    ur(i)=cu(i,16)/cu(i,15)
    vr(i)=cu(i,17)/cu(i,15)
    wr(i)=cu(i,18)/cu(i,15)
    us(i)=cu(i,20)/cu(i,19)
    vs(i)=cu(i,21)/cu(i,19)
    ws(i)=cu(i,22)/cu(i,19)
    
#else
   
   !For 2 phase
   bbb = (foursigmabyd + pinf * gamap) - (cu(i, 1) * (gamag - 1.0d0) * eg) - (cu(i, 6) * (gamap - 1.0d0) * ep)
   ccc = (cu(i, 1) * (gamag - 1.0d0) * eg * (foursigmabyd + pinf * gamap))
   pg(i) = (-bbb + dsqrt(bbb**2 + 4*ccc)) * 0.5d0
   pp(i) = pg(i) + foursigmabyd
   phig(i) = (cu(i, 1) * (gamag - 1.0d0) * eg) / pg(i)
   phip(i) = 1.0d0 - phig(i)

   rog(i) = cu(i, 1) / phig(i)
   rop(i) = cu(i, 6) / phip(i)
   tg(i) = pg(i) / (rog(i) * rrg)
   tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0d0) * rop(i) * cpp)
 
#endif   
    
   if (dabs(ycel(i) - 0.065d0) .lt. 0.02 .and. dabs(zcel(i) - 0.065d0) .lt. 0.02) then
    romix = phig(i) * rog(i) + (1.0d0 - phig(i)) * rop(i)
    pmix = phig(i) * pg(i) + (1.0d0 - phig(i)) * pp(i)

   endif
  end do

! CPU Updates - Calculation loop
!  do i = 1, nodes
!   itest(i) = 0
!   do j= 1,neq
!   cn(i, j) = 0.0d0
!    end do
!   enddo
   

! Processing element data
!  do i = 1, neles
!   do j = 1, 4
!    k = nod(i, j)
!    itest(k) = itest(k) + 1
!    do l = 1, neq
!     cn(k, l) = cn(k, l) + cu(i, l)
!    end do
!   end do
!  end do

! Final data processing
!  do i = 1, nodes
!   do j = 1, neq
!    cn(i, j) = cn(i, j) / itest(i)
!   end do
!  end do

! Write the header for the Tecplot file

write(fname_tec,"(A,I4.4,A)") trim(tecplot_dir)//'tecplot', myid, '.dat'
open(unit=25, file=fname_tec, status='unknown', action='write')

!---------------------------------------------------------
! Base variables (always present)
!---------------------------------------------------------
varList = 'VARIABLES = "X", "Y", "Z", "U","V", "W", "PG", "PP", "TG", ' // &
          '"UP", "VP", "WP", "PHIG", "TP", "ROG", "ROP"'
writevar = 16          

#ifdef FIVE_PHASE
  varList = trim(varList) // ', "PHIP", "PHIQ", "PHIR", "PHIS"'
  writevar = writevar + 4
#endif

#ifdef KE_TURB
  varList = trim(varList) // ', "TKG", "TEG"'
  writevar = writevar + 2
#endif

#ifdef NP_P
  varList = trim(varList) // ', "ENP"'
  writevar = writevar + 1
#endif

#ifdef YP_P
  varList = trim(varList) // ', "YFG", "YOG", "YPG"'
  writevar = writevar + 3
#endif

!=========================================================
! Construct Tecplot VARIABLE list based on physics options
!=========================================================

! ================================
! HEADER
! ================================

write(25,*) trim(varList)

write(25,'(A,I0,A,I0,A,I0,A)') &
 'ZONE ZONETYPE=FETETRAHEDRON, NODES=', nodes, &
 ', ELEMENTS=', neles, &
', DATAPACKING=BLOCK, VARLOCATION=([4-', writevar, ']=CELLCENTERED)'
do ie = 1, neles
  ug(ie) = cu(ie,2) / cu(ie,1)
  vg(ie) = cu(ie,3) / cu(ie,1)
  wg(ie) = cu(ie,4) / cu(ie,1)

  up(ie) = cu(ie,7) / cu(ie,6)
  vp(ie) = cu(ie,8) / cu(ie,6)
  wp(ie) = cu(ie,9) / cu(ie,6)

  eg = cu(ie,5) / cu(ie,1) - 0.5d0*(ug(ie)**2+vg(ie)**2+wg(ie)**2)
  ep = cu(ie,10)/ cu(ie,6) - 0.5d0*(up(ie)**2+vp(ie)**2+wp(ie)**2)

#ifdef KE_TURB
  tkg(ie) = cu(ie,nkg)/cu(ie,1)
  teg(ie) = cu(ie,neg)/cu(ie,1)
#endif

#ifdef NP_P
  enp(ie) = cu(ie,npp)
#endif

#ifdef YP_P
  yfg(ie)=cu(ie,nyfg)/cu(ie,1)
  yog(ie)=cu(ie,nyog)/cu(ie,1)
  ypg(ie)=cu(ie,nypg)/cu(ie,1)
#endif
  
#ifdef FIVE_PHASE
!   if(myflux.eq.3) then
   ! For 5 phase
    phiq(ie)=dmin1(dmax1(cu(ie,11)/rhoq,epslnmin),1.0d0-epslnmin)
    phir(ie)=dmin1(dmax1(cu(ie,15)/rhor,epslnmin),1.0d0-epslnmin)
    phis(ie)=dmin1(dmax1(cu(ie,19)/rhos,epslnmin),1.0d0-epslnmin)
    	
   aaa=dmin1(dmax1(1.0d0-phiq(ie)-phir(ie)-phis(ie),epslnmin),1.0d0-epslnmin) 	
   bbb=(aaa*pinf*gamap)-(cu(i,1)*(gamag-1.0d0)*eg)-(cu(i,6)*(gamap-1.0d0)*ep)
   ccc=(cu(i,1)*(gamag-1.0d0)*eg*(gamap*pinf))
   pg(ie)=(-bbb+dsqrt(bbb**2+4.0d0*aaa*ccc))/(2.0d0*aaa)
   pp(ie) = pg(ie) + foursigmabyd
 
!=====================================================================

   ! Compute densities and temperatures      
   rog(ie) = pg(ie) / ((gamag - 1.0d0)*eg)
   rop(ie) = (pp(ie) + pinf*gamap)  / ((gamap - 1.0d0)*ep)
   
   tg(ie) = pg(ie) / (rog(ie) * rrg)
   tp(ie) = ((pp(ie) + pinf) * gamap) / ((gamap - 1.0_dp) * rop(ie) * cpp)
  
   ! Compute phase fractions
   phig(ie) = dmin1(dmax1(cu(ie,1)/rog(ie),epslnmin),(1-epslnmin))
   phip(ie) = dmin1(dmax1(cu(ie,6)/rop(ie),epslnmin),(1-epslnmin)) 
   
   !Jameson
   phitotal=phig(ie)+phip(ie)+phiq(ie)+phir(ie)+phis(ie)
   phig(ie)=phig(ie)/phitotal
   phip(ie)=phip(ie)/phitotal
   phiq(ie)=phiq(ie)/phitotal
   phir(ie)=phir(ie)/phitotal
   phis(ie)=phis(ie)/phitotal
   
   uq(ie)=cu(ie,12)/cu(ie,11)
   vq(ie)=cu(ie,13)/cu(ie,11)
   wq(ie)=cu(ie,14)/cu(ie,11)
   ur(ie)=cu(ie,16)/cu(ie,15)
   vr(ie)=cu(ie,17)/cu(ie,15)
   wr(ie)=cu(ie,18)/cu(ie,15)
   us(ie)=cu(ie,20)/cu(ie,19)
   vs(ie)=cu(ie,21)/cu(ie,19)
   ws(ie)=cu(ie,22)/cu(ie,19)
#else
!For 2 phase
      
   bbb = (foursigmabyd + pinf * gamap) - (cu(ie, 1) * (gamag - 1.0d0) * eg) - (cu(ie, 6) * (gamap - 1.0d0) * ep)
   ccc = (cu(ie, 1) * (gamag - 1.0d0) * eg * (foursigmabyd + pinf * gamap))
   pg(ie) = (-bbb + dsqrt(bbb**2 + 4*ccc)) * 0.5d0
   pp(ie) = pg(ie) + foursigmabyd

   phig(ie) = (cu(ie, 1) * (gamag - 1.0d0) * eg) / pg(ie)
   phip(ie) = 1.0d0 - phig(ie)

   rog(ie) = cu(ie, 1) / phig(ie)
   rop(ie) = cu(ie, 6) / phip(ie)
   tg(ie) = pg(ie) / (rog(ie) * rrg)
   tp(ie) = ((pp(ie) + pinf) * gamap) / ((gamap - 1.0d0) * rop(ie) * cpp)
#endif

!============================================================= 

   IF (phig(ie) < epslnmin) THEN
    phig(ie) = epslnmin
    phip(ie) = 1.0d0 - phig(ie)
    ug(ie) = up(ie)
    vg(ie) = vp(ie)
    wg(ie) = wp(ie)
    tg(ie) = tp(ie)
    pg(ie) = pp(ie)
   END IF

   IF (phip(ie) < epslnmin) THEN
    phip(ie) = epslnmin
    phig(ie) = 1.0d0 - phip(ie)
    up(ie) = ug(ie)
    vp(ie) = vg(ie)
    wp(ie) = wg(ie)
    tp(ie) = tg(ie)
    pp(ie) = pg(ie)
   END IF

   IF (phig(ie) >= epslnmin .AND. phig(ie) <= epslnmax) THEN
    psig = (phig(ie) - epslnmin) / (epslnmax - epslnmin)
    functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
    ug(ie) = functiong * ug(ie) + (1.0d0 - functiong) * up(ie)
    vg(ie) = functiong * vg(ie) + (1.0d0 - functiong) * vp(ie)
    wg(ie) = functiong * wg(ie) + (1.0d0 - functiong) * wp(ie)
    tg(ie) = functiong * tg(ie) + (1.0d0 - functiong) * tp(ie)
    pg(ie) = functiong * pg(ie) + (1.0d0 - functiong) * pp(ie)
   END IF

   IF (phip(ie) >= epslnmin .AND. phip(ie) <= epslnmax) THEN
    psig = (phip(ie) - epslnmin) / (epslnmax - epslnmin)
    functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
    up(ie) = functiong * up(ie) + (1.0d0 - functiong) * ug(ie)
    vp(ie) = functiong * vp(ie) + (1.0d0 - functiong) * vg(ie)
    wp(ie) = functiong * wp(ie) + (1.0d0 - functiong) * wg(ie)
    tp(ie) = functiong * tp(ie) + (1.0d0 - functiong) * tg(ie)
    pp(ie) = functiong * pp(ie) + (1.0d0 - functiong) * pg(ie)
   END IF

    ! Compute densities and mixed density
   rog(ie) = pg(ie) / (rrg * tg(ie))
   rop(ie) = (pp(ie) + pinf) * gamap / ((gamap - 1.0d0) * cpp * tp(ie))

end do  


! ================================
! NODAL COORDINATES
! ================================
do i = 1, nodes
  write(25,'(ES16.8)') x(i)
end do

do i = 1, nodes
  write(25,'(ES16.8)') y(i)
end do

do i = 1, nodes
  write(25,'(ES16.8)') z(i)
end do

! ================================
! CELL-CENTERED DATA
! ================================
do ie = 1, neles
  write(25,'(ES16.8)') ug(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') vg(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') wg(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') pg(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') pp(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') tg(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') up(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') vp(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') wp(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') phig(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') tp(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') rog(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') rop(ie)
end do

!---------------------------------------------------------
! Five-phase variables
!---------------------------------------------------------
#ifdef FIVE_PHASE
do ie = 1, neles
  write(25,'(ES16.8)') phip(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') phiq(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') phir(ie)
end do

do ie = 1, neles
  write(25,'(ES16.8)') phis(ie)
end do
#endif

!---------------------------------------------------------
! k–epsilon turbulence variables
!---------------------------------------------------------
#ifdef KE_TURB
do ie = 1, neles
  write(25,'(ES16.8)') tkg(ie)
end do
do ie = 1, neles
  write(25,'(ES16.8)') teg(ie)
end do
#endif

!---------------------------------------------------------
! Nanoparticle number density
!---------------------------------------------------------
#ifdef NP_P
do ie = 1, neles
  write(25,'(ES16.8)') enp(ie)
end do
#endif

!---------------------------------------------------------
! Species mass fractions
!---------------------------------------------------------
#ifdef YP_P
do ie = 1, neles
  write(25,'(ES16.8)') yfg(ie)
end do
do ie = 1, neles
  write(25,'(ES16.8)') yog(ie)
end do
do ie = 1, neles
  write(25,'(ES16.8)') ypg(ie)
end do
#endif

  ! Loop over elements for nodal data
DO ie = 1, neles
   WRITE(25, *) nod(ie, 1), nod(ie, 2), nod(ie, 3), nod(ie, 4)
END DO   
CALL tioga_fin()
 CLOSE(25)
 CLOSE(9)
 CLOSE(13)
 
  print *, 'Tecplot file written: ', trim(fname_tec)

CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
CALL DEALLOCATE_Variables()
CALL MPI_FINALIZE(ierr)   

END PROGRAM Main
