!===========================================================================!
!          Runge kutte 4 stage 3 Order (ESSPRK)         	            !
!===========================================================================!
		 	
SUBROUTINE RK
    USE GlobalVariables
    USE FluxModule 
    IMPLICIT NONE
    external geometry
    external primitive
    external smooth
    
    INTEGER :: i, j
    REAL(DP) :: dtbyvol

!=======================
! Runge-Kutta stage one
!=======================
CALL VISCOUSGP

CALL CallFluxSolver()

#ifdef NP_P
    CALL mpiexcn
#endif

!$acc parallel loop present(cn1(:,:),cn(:,:),d24(:,:),cu(:,:),rhs(:,:),neq, &
!$acc vol(:),dt,neles,neg,toleg,tolkg,nkg,npp,enpcutoff,nyfg,nyog,nypg,tolsp) private(i,dtbyvol,j)  
do i = 1, neles

if (iblank(i)/=1)cycle !-------(mod3) iblank check

    dtbyvol = dt / vol(i)
!$acc loop seq
    do j = 1, neq
        CN(i, j) = Cu(i, j) - (11.0D0 / 20.0D0) * dtbyvol * (RHS(i, j) - D24(i, j))
        CN1(i, j) = CN(i, j)
    end do
#ifdef KE_TURB
    if(cn(i,nkg)/cn(i,1).lt.tolkg)cn(i,nkg)=cn(i,1)*tolkg
    if(cn(i,neg)/cn(i,1).lt.toleg)cn(i,neg)=cn(i,1)*toleg
#endif	
#ifdef NP_P
    if(cn(i,npp).lt.enpcutoff)cn(i,npp)=enpcutoff
#endif
#ifdef YP_P
        if(cn(i,nyfg)/cn(i,1).lt.tolsp)cn(i,nyfg)=cn(i,1)*tolsp
	if(cn(i,nyog)/cn(i,1).lt.tolsp)cn(i,nyog)=cn(i,1)*tolsp
	if(cn(i,nypg)/cn(i,1).lt.tolsp)cn(i,nypg)=cn(i,1)*tolsp
	if(cn(i,nyfg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyfg)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nyog)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyog)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nypg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nypg)=cn(i,1)*(1.0d0-tolsp)
#endif
end do
!$acc end parallel loop

CALL mpiexcn


!=======================
! Runge-Kutta stage two
!=======================
CALL CallFluxSolver()

!$acc parallel loop present(cn(:,:),d24(:,:),cu(:,:),rhs(:,:),neq, &
!$acc vol(:),dt,neles,neg,toleg,tolkg,nkg,npp,enpcutoff,nyfg,nyog,nypg,tolsp) &
!$acc private(i,dtbyvol,j) 
do i = 1, neles
if (iblank(i)/=1)cycle !-------(mod3) iblank check
    dtbyvol = dt / vol(i)
!$acc loop seq
    do j = 1, neq
       CN(i, j) = (3.0D0/8.0D0) * Cu(i, j) + (5.0D0/8.0D0) * (CN(i, j) - (11.0D0/20.0D0) * dtbyvol * (RHS(i, j) - D24(i, j)))
    end do
#ifdef KE_TURB
    if(cn(i,nkg)/cn(i,1).lt.tolkg)cn(i,nkg)=cn(i,1)*tolkg
    if(cn(i,neg)/cn(i,1).lt.toleg)cn(i,neg)=cn(i,1)*toleg
#endif
#ifdef NP_P
    if(cn(i,npp).lt.enpcutoff)cn(i,npp)=enpcutoff
#endif
#ifdef YP_P
       if(cn(i,nyfg)/cn(i,1).lt.tolsp)cn(i,nyfg)=cn(i,1)*tolsp
	if(cn(i,nyog)/cn(i,1).lt.tolsp)cn(i,nyog)=cn(i,1)*tolsp
	if(cn(i,nypg)/cn(i,1).lt.tolsp)cn(i,nypg)=cn(i,1)*tolsp
	if(cn(i,nyfg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyfg)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nyog)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyog)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nypg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nypg)=cn(i,1)*(1.0d0-tolsp)
#endif
end do
!$acc end parallel loop  

CALL mpiexcn

!======================== 
! Runge-Kutta stage three
!========================
CALL CallFluxSolver()
!$acc parallel loop present(cu(:,:),d24(:,:),cn(:,:),neq, &
!$acc rhs(:,:),neles,dt,vol(:),neg,toleg,tolkg,nkg,npp,enpcutoff,nyfg,nyog,nypg,tolsp) &
!$acc private(j,dtbyvol,i)
do i = 1, neles
if (iblank(i)/=1)cycle !-------(mod3) iblank check
    dtbyvol = dt / vol(i)
!$acc loop seq     
    do j = 1, neq
        CN(i, j) = (4.0D0/9.0D0) * Cu(i, j) + (5.0D0/9.0D0) * (CN(i, j) - (11.0D0/20.0D0) * dtbyvol * (RHS(i, j)- D24(i,j)))
    end do
#ifdef KE_TURB
    if(cn(i,nkg)/cn(i,1).lt.tolkg)cn(i,nkg)=cn(i,1)*tolkg
    if(cn(i,neg)/cn(i,1).lt.toleg)cn(i,neg)=cn(i,1)*toleg
#endif
#ifdef NP_P
    if(cn(i,npp).lt.enpcutoff)cn(i,npp)=enpcutoff
#endif
#ifdef YP_P
        if(cn(i,nyfg)/cn(i,1).lt.tolsp)cn(i,nyfg)=cn(i,1)*tolsp
	if(cn(i,nyog)/cn(i,1).lt.tolsp)cn(i,nyog)=cn(i,1)*tolsp
	if(cn(i,nypg)/cn(i,1).lt.tolsp)cn(i,nypg)=cn(i,1)*tolsp
	if(cn(i,nyfg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyfg)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nyog)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyog)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nypg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nypg)=cn(i,1)*(1.0d0-tolsp)
#endif

end do
!$acc end parallel loop

CALL mpiexcn

!=======================
! Runge-Kutta stage four
!=======================
CALL CallFluxSolver()

!$acc parallel loop present(cu(:,:),d24(:,:),cn(:,:),neq, &
!$acc cn1(:,:),vol(:),rhs(:,:),neles,dt,neg,toleg,tolkg,nkg,npp,enpcutoff,nyfg,nyog,nypg,tolsp)  &        
!$acc private(dtbyvol,i,j)
do i = 1, neles
if (iblank(i)/=1)cycle !-------(mod3) iblank check
    dtbyvol = dt / vol(i)
!$acc loop seq 
    do j = 1, neq
        CN(i, j) = (111.0D0/1331.0D0) * Cu(i, j) + (260.0D0/1331.0D0) * CN1(i, j) + (960.0D0/1331.0D0) * (CN(i, j) - &
                   (11.0D0 / 20.0D0) * dtbyvol * (RHS(i, j) - D24(i, j)))
    end do
#ifdef KE_TURB
    if(cn(i,nkg)/cn(i,1).lt.tolkg)cn(i,nkg)=cn(i,1)*tolkg
    if(cn(i,neg)/cn(i,1).lt.toleg)cn(i,neg)=cn(i,1)*toleg
#endif    
#ifdef NP_P
    if(cn(i,npp).lt.enpcutoff)cn(i,npp)=enpcutoff
#endif
#ifdef YP_P
        if(cn(i,nyfg)/cn(i,1).lt.tolsp)cn(i,nyfg)=cn(i,1)*tolsp
	if(cn(i,nyog)/cn(i,1).lt.tolsp)cn(i,nyog)=cn(i,1)*tolsp
	if(cn(i,nypg)/cn(i,1).lt.tolsp)cn(i,nypg)=cn(i,1)*tolsp
	if(cn(i,nyfg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyfg)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nyog)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nyog)=cn(i,1)*(1.0d0-tolsp)
	if(cn(i,nypg)/cn(i,1).gt.(1.0d0-tolsp))cn(i,nypg)=cn(i,1)*(1.0d0-tolsp)
#endif    
end do
!$acc end parallel loop 

CALL mpiexcn

return
end subroutine RK

