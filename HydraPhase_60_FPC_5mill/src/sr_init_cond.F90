SUBROUTINE Init_Cond
USE GlobalVariables
IMPLICIT NONE
 
INTEGER :: i    
REAL(DP) :: phigmin
REAL(DP) :: TKGMAX,TKGMIN,TEGMAX,TEGMIN,amug 
REAL(DP) :: Epsilonmax,epsilonmin

!This loop need for istart more than 0
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
#endif

!============ 
! Shock Tube
!============
#ifdef ISHOCKTUBE
!====================
!LEFT SIDE SHOCK TUBE
!====================

do i = 1, neles
   if(xcel(i).le.5.0d0)then
   cu(i,1) = phigmax * rogmax
   cu(i,2) = phigmax * rogmax * ugmax
   cu(i,3) = phigmax * rogmax * vgmax
   cu(i,4) = phigmax * rogmax * wgmax
   cu(i,5) = phigmax * rogmax * egtotalmax
   cu(i,6) = (1.0D0 - phigmax) * ropmax
   cu(i,7) = (1.0D0 - phigmax) * ropmax * upmax
   cu(i,8) = (1.0D0 - phigmax) * ropmax * vpmax
   cu(i,9) = (1.0D0 - phigmax) * ropmax * wpmax
   cu(i,10) = (1.0D0 - phigmax) * ropmax * eptotalmax
#ifdef KE_TURB
   cu(i,11) = phigmax * rogmax * tkgmax
   cu(i,12) = phigmax * rogmax * tkgmax
#endif    
   else
   cu(i,1) = phigmin * rogmin
   cu(i,2) = phigmin * rogmin * ugmin
   cu(i,3) = phigmin * rogmin * vgmin
   cu(i,4) = phigmin * rogmin * wgmin
   cu(i,5) = phigmin * rogmin * egtotalmin
   cu(i,6) = (1.0D0 - phigmin) * ropmin
   cu(i,7) = (1.0D0 - phigmin) * ropmin * upmin
   cu(i,8) = (1.0D0 - phigmin) * ropmin * vpmin
   cu(i,9) = (1.0D0 - phigmin) * ropmin * wpmin
   cu(i,10) = (1.0D0 - phigmin) * ropmin * eptotalmin
#ifdef KE_TURB
   cu(i,11) = phigmin * rogmin * tkgmin
   cu(i,12) = phigmin * rogmin * tegmin
#endif    
   endif
enddo

! do i = 1, neles
!   if(xcel(i).le.5.0d0)then
   
!   phig(i) = phigmax 
!   phip(i) = 1-phig(i) 
!#ifdef FIVE_PHASE
!   phiq(i) = 4.92D-8
!   phir(i) = epslninf
!   phis(i) = epslninf
!   phig(i) = 1- phip(i) - phiq(i) - phir(i) - phis(i) 
!#endif
!   cu(i,1) = phig(i) * rogmax
!   cu(i,2) = phig(i) * rogmax * ugmax
!   cu(i,3) = phig(i) * rogmax * vgmax
!   cu(i,4) = phig(i) * rogmax * wgmax
!   cu(i,5) = phig(i) * rogmax * egtotalmax
!   cu(i,6) = phip(i) * ropmax
!   cu(i,7) = phip(i) * ropmax * upmax
!   cu(i,8) = phip(i) * ropmax * vpmax
!   cu(i,9) = phip(i) * ropmax * wpmax
!   cu(i,10) = phip(i) * ropmax * eptotalmax
!#ifdef FIVE_PHASE
!   cu(i,11) = phiq(i)*rhoq
!   cu(i,12) = phiq(i)*rhoq * upmax
!   cu(i,13) = phiq(i)*rhoq * vpmax
!   cu(i,14) = phiq(i)*rhoq * wpmax
!   cu(i,15) = phir(i)*rhor
!   cu(i,16) = phir(i)*rhor * upmax
!   cu(i,17) = phir(i)*rhor * vpmax
!   cu(i,18) = phir(i)*rhor * wpmax
!   cu(i,19) = phis(i)*rhos
!   cu(i,20) = phis(i)*rhos * upmax
!   cu(i,21) = phis(i)*rhos * vpmax
!   cu(i,22) = phis(i)*rhos * wpmax
!#endif
 
!#ifdef KE_TURB
!   cu(i,nkg) = phig(i) * rogmax * tkgmax
!   cu(i,neg) = phig(i) * rogmax * tegmax
!#endif 

!#ifdef NP_P
!   cu(i,npp) = (1.0D0 - phigmax)/volrp
!#endif 
!#ifdef YP_P
!   cu(i,nyfg) = phigmax * rogmax * yfg_inf
!   cu(i,nyog) = phigmax * rogmax * yog_inf
!   cu(i,nypg) = phigmax * rogmax * ypg_inf
!#endif   


!RIGHT SIDE SHOICK TUBE
!   else
   
!   phig(i) = phigmin 
!   phip(i) = 1- phig(i) 
!#ifdef FIVE_PHASE
!   phiq(i) = 4.92D-4 
!   phir(i) = epslninf
!   phis(i) = epslninf 
!   phig(i) = 1 -  phip(i) -  phiq(i) -  phir(i) -  phis(i) 
!#endif
!   cu(i,1) = phig(i) * rogmin
!   cu(i,2) = phig(i) * rogmin * ugmin
!   cu(i,3) = phig(i) * rogmin * vgmin
!   cu(i,4) = phig(i) * rogmin * wgmin
!   cu(i,5) = phig(i) * rogmin * egtotalmin
!   cu(i,6) = phip(i) * ropmin
!   cu(i,7) = phip(i) * ropmin * upmin
!   cu(i,8) = phip(i)* ropmin * vpmin
!   cu(i,9) = phip(i) * ropmin * wpmin
!   cu(i,10) = phip(i) * ropmin * eptotalmin
!#ifdef FIVE_PHASE
!   cu(i,11) = phiq(i)*rhoq
!   cu(i,12) = phiq(i)*rhoq * upmin
!   cu(i,13) = phiq(i)*rhoq * vpmin
!   cu(i,14) = phiq(i)*rhoq * wpmin
!   cu(i,15) = phir(i)*rhor
!   cu(i,16) = phir(i)*rhor * upmin
!   cu(i,17) = phir(i)*rhor * vpmin
!   cu(i,18) = phir(i)*rhor * wpmin
!   cu(i,19) = phis(i)*rhos
!   cu(i,20) = phis(i)*rhos * upmin
!   cu(i,21) = phis(i)*rhos * vpmin
!   cu(i,22) = phis(i)*rhos * wpmin
!#endif   
!#ifdef KE_TURB
!   cu(i,nkg) = phig(i) * rogmin * tkgmin
!   cu(i,neg) = phig(i) * rogmin * tegmin
!#endif  

!#ifdef NP_P
!   cu(i,npp) = (1.0D0 - phigmin)/volrp
!#endif 
!#ifdef YP_P
!   cu(i,nyfg) = phigmin * rogmin * yfg_inf
!   cu(i,nyog) = phigmin * rogmin * yog_inf
!   cu(i,nypg) = phigmin * rogmin * ypg_inf 
!#endif    
!    endif
!enddo


#else    
!===============================================================================
!Other than shock tube problem 
!===============================================================================


  DO i = 1, neles
   phip(i) = epslninf
   phig(i) = 1 - epslninf
#ifdef FIVE_PHASE
   phiq(i) = epslninf
   phir(i) = epslninf
   phis(i) = epslninf
   phig(i) = 1- phip(i) - phiq(i) - phir(i) - phis(i)   
#endif
   cu(i,1) = phig(i) * rogmin
   cu(i,2) = phig(i) * rogmin * ugmin
   cu(i,3) = phig(i) * rogmin * vgmin
   cu(i,4) = phig(i) * rogmin * wgmin
   cu(i,5) = phig(i) * rogmin * egtotalmin
   cu(i,6) = phip(i) * ropmin
   cu(i,7) = phip(i) * ropmin * upmin
   cu(i,8) = phip(i) * ropmin * vpmin
   cu(i,9) = phip(i) * ropmin * wpmin
   cu(i,10) = phip(i) * ropmin * eptotalmin
#ifdef FIVE_PHASE
   cu(i,11) = phiqinit*rhoq
   cu(i,12) = phiqinit*rhoq * upmin
   cu(i,13) = phiqinit*rhoq * vpmin
   cu(i,14) = phiqinit*rhoq * wpmin
   cu(i,15) = phirinit*rhor
   cu(i,16) = phirinit*rhor * upmin
   cu(i,17) = phirinit*rhor * vpmin
   cu(i,18) = phirinit*rhor * wpmin
   cu(i,19) = phisinit*rhos
   cu(i,20) = phisinit*rhos * upmin
   cu(i,21) = phisinit*rhos * vpmin
   cu(i,22) = phisinit*rhos * wpmin
#endif

#ifdef YP_P
   cu(i,nyfg)= 0.99
   cu(i,nyog)= 0.01
   cu(i,nypg)= 0.0 
#endif   

#ifdef NP_P
   cu(i,npp) = (1.0D0 - phigmax)/volrp
#endif 

#ifdef KE_TURB
   cu(i,nkg) = phig(i) * rogmax * tkgmax
   cu(i,neg) = phig(i) * rogmax * tegmax
#endif

End do

!if(myid==0.and.i.lt.100) print*,"cu(i,nkg)",i, nkg, cu(i,nkg),phig(i),rogmax, tkgmax 
#endif   


END SUBROUTINE Init_Cond

