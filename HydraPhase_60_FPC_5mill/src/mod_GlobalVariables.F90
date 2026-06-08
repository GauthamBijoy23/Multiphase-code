MODULE GlobalVariables
implicit none
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)
  
  INTEGER :: NN,NB,NODES,NELES,NGHOSTS,NTOT,MYFLUX,MYORDER,ISTART,TOTAL,NVAR,NITER,ITTER,ITER
  INTEGER :: ITERS,ITERSS,ITERSUB,NEPG,NPARTITIONS,MYID,NPROCS,IERR,PGHOSTS,MAX_PGHOSTS
  INTEGER, ALLOCATABLE :: nparent(:),nghost(:),ntype(:),nside(:),nc1(:),nc2(:),nc3(:),nc4(:),nod(:,:)
  INTEGER, ALLOCATABLE :: g_id(:),iblank(:)
  REAL(DP), ALLOCATABLE :: x(:),y(:),z(:),xcel(:),ycel(:)
  REAL(DP), ALLOCATABLE :: zcel(:),vol(:),dl(:)
  REAL(DP), ALLOCATABLE :: sc1x(:),sc1y(:),sc1z(:),sc2x(:),sc2y(:),sc2z(:),sc3x(:),sc3y(:),sc3z(:),sc4x(:),sc4y(:)
  REAL(DP), ALLOCATABLE :: sc4z(:),itest(:)
  REAL(DP), ALLOCATABLE :: ug(:),vg(:),wg(:),pg(:),pp(:),tg(:),phig(:),up(:),vp(:),wp(:),phip(:),tp(:),rop(:),rog(:)
  REAL(DP), ALLOCATABLE :: ag(:),ap(:)
  REAL(DP), ALLOCATABLE :: cn(:,:),rhs(:,:),d24(:,:),anu(:,:),d2q(:,:),cn1(:,:),cu(:,:)
  REAL(DP), ALLOCATABLE :: visug(:),visvg(:),viswg(:),vistg(:),phi(:)
  REAL(DP), ALLOCATABLE :: visup(:),visvp(:),viswp(:),vistp(:)
  REAL(DP), ALLOCATABLE :: delrogy(:),delrogx(:),delphigy(:)
  REAL(DP), ALLOCATABLE :: delphigx(:),delpgz(:),delvgy(:),delvgx(:),deltgy(:),deltgx(:),delphigz(:)
  REAL(DP), ALLOCATABLE :: delwgz(:),delugz(:),delvgz(:),delTgz(:),delwgx(:),delugx(:),delrogz(:),delwgy(:)
  REAL(DP), ALLOCATABLE :: delugy(:),delpgy(:),delupz(:),delppz(:),delppx(:),delppy(:),delwpx(:),delwpy(:)
  REAL(DP), ALLOCATABLE :: delupx(:),delupy(:),delvpy(:),delvpx(:),deltpy(:),deltpx(:),delPgx(:),delvpz(:)
  REAL(DP), ALLOCATABLE :: delTpz(:),delropy(:),delropx(:),delwpz(:),delropz(:)
  REAL(DP) :: pinf,rrg,gamag,gamap,cpg,cpp,epslnmax,epslnmin,pmax,pmin,tgmax,tgmin
  REAL(DP) :: tpmax,tpmin,ugmax,vgmax,wgmax,ropmax,upmax,vpmax,wpmax,upmin,vpmin
  REAL(DP) :: rogmin,egtotalmin,ugmin,vgmin,wgmin,phigmax,cmug,prg,wpmin,xmax,ymax,zmax,ymin,zmin
  REAL(DP) :: gravity,c23,diap,foursigmabyd,time,dt,phipzero,as,eptotalmax,eptotalmin, dtmin_global
  REAL(DP) ::  epslninf, cfl,surfacetension, prp,dtmin,egtotalmax,rogmax,ropmin
  REAL(DP) :: diffp,diffq,diffr,diffu,diffv,diffw,diffx,diffy,diffz,difft,co2,co4,ts
  REAL(DP) :: amolg,Vel_mag

  INTEGER, ALLOCATABLE:: pp_elem(:), global_pp(:), neigh_proc(:), neigh_side(:), loc_pghost(:), neigh_glob(:),dest_ghost(:)
  
  INTEGER, ALLOCATABLE :: send_count(:), recv_count(:)
  INTEGER, ALLOCATABLE :: send_pos(:),   recv_pos(:)

  REAL(DP), ALLOCATABLE :: sendbuf(:,:), recvbuf(:,:)
  REAL(DP), ALLOCATABLE :: send_data(:,:,:), recv_data(:,:,:)
   
  INTEGER, ALLOCATABLE :: cu_send_pos(:), cu_recv_pos(:)
  INTEGER, ALLOCATABLE :: cu_send_count(:), cu_recv_count(:)
  REAL(DP), ALLOCATABLE :: cu_recv_data(:,:,:)

  REAL(DP), ALLOCATABLE :: uq(:),vq(:),wq(:),ur(:),vr(:),wr(:),us(:),vs(:),ws(:),phiq(:),phir(:),phis(:),dd(:)
  REAL(DP) :: rhoq, rhor, rhos, radq, radr, rads, volq, volr, vols, qmass, rmass, smass
  REAL(DP) :: phiginit,phipinit,phiqinit, phirinit, phisinit
   
  INTEGER :: nkg,neg,neq
  REAL(DP) :: TOLKG,TOLEG
  REAL(DP), ALLOCATABLE :: vistkg(:),visteg(:),tkg(:),teg(:),egtotal(:),eptotal(:)
  REAL(DP), ALLOCATABLE :: delTkgx(:),delTkgy(:),delTkgz(:),delTegx(:),delTegy(:),delTegz(:)
  
  
!Variable for NP
  INTEGER :: npp
  REAL(DP) :: enpcutoff,radp,volrp
  REAL(DP), ALLOCATABLE :: enp(:)
  
!Variable for Species transport  
  INTEGER :: nyfg,nyog,nypg
  REAL(DP) :: tolsp,cpfg,cpog,cppg,amolfg,amolog,amolpg,ypg_inf,yfg_inf,yog_inf
  REAL(DP), ALLOCATABLE :: yfg(:),yog(:),ypg(:),visyfg(:),visyog(:),visypg(:)
   
  REAL(DP)::cp_t,cp_o,cp_p,amol_f,amol_o,amol_p
  REAL(DP)::sratio,hfp,preexp,ebyr
       
  INTEGER, ALLOCATABLE :: bc2_list(:), bc3_list(:), bc5_list(:)
  INTEGER, ALLOCATABLE :: bc51_list(:), bc53_list(:)
  INTEGER :: n_bc2, n_bc3, n_bc5, n_bc51, n_bc53
   
  END MODULE GlobalVariables
