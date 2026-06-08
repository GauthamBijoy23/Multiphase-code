SUBROUTINE primitiveplusbc
    USE GlobalVariables
    IMPLICIT NONE
    external geometry
    
    integer :: i,k,ng,np,jnt,ns,n1,n2,n3,n4,j,nt,startl,endl,ii
    REAL(DP):: functiong,bbb,psig,ccc,ep,eg,x1,x2,x4,x3,xcenter
    REAL(DP) :: d,wt,wu,alphat,deltat,alphau,fract,fracu,aaa,amug,cvg
    REAL(DP) :: y1,y2,y3,y4,z1,z2,z3,z4,ycenter,zcenter
    logical :: active


!$acc parallel loop present(ag(:),tg(:),rog(:),phig(:),phip(:),wp(:), &
!$acc wg(:),ug(:),up(:),ap(:),tp(:),vp(:),vg(:),cn(:,:),pg(:),pp(:),  &
!$acc rop(:),neles,cpp,gamag,gamap,rrg,epslnmax,pinf,tkg(:),teg(:),neg,nkg,egtotal,eptotal, &
!$acc foursigmabyd,epslnmin,nghosts,ntot,rhoq,rhos,rhor,phiq(:),phis(:),phir(:),enp(:),npp, &
!$acc yfg(:),yog(:),ypg(:),cpg,cpfg,cpog,cppg,rrg,amolfg,amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg) &
!$acc  private(functiong,bbb,psig,ccc,ep,eg,i,startl,endl,aaa,cvg,active)


do i = 1, ntot

   active = ( i <= neles ) .or. ( i > neles + nghosts )
   if (.not. active) cycle
!do k = 1,2

!        if (k .eq. 1) then
!           startl = 1
!           endl   = neles
!        else
!           startl = neles + nghosts + 1
!           endl   = ntot
!        end if

!  do i = startl, endl
   ug(i) = cn(i,2) / cn(i,1)
   vg(i) = cn(i,3) / cn(i,1)
   wg(i) = cn(i,4) / cn(i,1)
   up(i) = cn(i,7) / cn(i,6)
   vp(i) = cn(i,8) / cn(i,6)
   wp(i) = cn(i,9) / cn(i,6)
   eg = cn(i,5) / cn(i,1) - 0.5_dp * (ug(i)**2 + vg(i)**2 + wg(i)**2)
   ep = cn(i,10) / cn(i,6) - 0.5_dp * (up(i)**2 + vp(i)**2 + wp(i)**2)

#ifdef KE_TURB
    tkg(i)=cn(i,nkg)/cn(i,1)
    teg(i)=cn(i,neg)/cn(i,1)
#endif

#ifdef NP_P
     enp(i) = cn(i,npp)
#endif

#ifdef YP_P
     yfg(i) = cn(i,nyfg)/cn(i,1)
     yog(i) = cn(i,nyog)/cn(i,1)
     ypg(i) = cn(i,nypg)/cn(i,1)
	 
     cpg=yfg(i)*cpfg+yog(i)*cpog+ypg(i)*cppg
     rrg=8314.0d0*(yfg(i)/amolfg+yog(i)/amolog+ypg(i)/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
#endif 

#ifdef FIVE_PHASE
!    if(myflux == 3) then
    phiq(i)=dmin1(dmax1(cn(i,11)/rhoq,epslnmin),1.0d0-epslnmin)
    phir(i)=dmin1(dmax1(cn(i,15)/rhor,epslnmin),1.0d0-epslnmin)
    phis(i)=dmin1(dmax1(cn(i,19)/rhos,epslnmin),1.0d0-epslnmin)
    	
   aaa=dmin1(dmax1(1.0d0-phiq(i)-phir(i)-phis(i),epslnmin),1.0d0-epslnmin) 	
   bbb=(aaa*pinf*gamap)-(cn(i,1)*(gamag-1.0d0)*eg)-(cn(i,6)*(gamap-1.0d0)*ep)
   ccc=(cn(i,1)*(gamag-1.0d0)*eg*(gamap*pinf))
   pg(i)=(-bbb+dsqrt(bbb**2+4.0d0*aaa*ccc))/(2.0d0*aaa)
   pp(i) = pg(i) + foursigmabyd
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Compute densities and temperatures      
   rog(i) = pg(i) / ((gamag - 1.0d0)*eg)
   rop(i) = (pp(i) + pinf*gamap)  / ((gamap - 1.0d0)*ep)
   
   tg(i) = pg(i) / (rog(i) * rrg)
   tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0d0) * rop(i) * cpp)
  
   ! Compute phase fractions
   phig(i) = dmin1(dmax1(cn(i,1)/rog(i),epslnmin),(1-epslnmin))
   phip(i) = dmin1(dmax1(cn(i,6)/rop(i),epslnmin),(1-epslnmin)) 

!For 2 phase
#else

   bbb = (foursigmabyd + pinf * gamap) - (cn(i,1) * (gamag - 1.0_dp) * eg) - &
            (cn(i,6) * (gamap - 1.0_dp) * ep)
   ccc = cn(i,1) * (gamag - 1.0_dp) * eg * (foursigmabyd + pinf * gamap)
   pg(i) = (-bbb + sqrt(bbb**2 + 4 * ccc)) * 0.5_dp
   pp(i) = pg(i) + foursigmabyd
   phig(i) = (cn(i,1) * (gamag - 1.0_dp) * eg) / pg(i)
   phip(i) = 1.0_dp - phig(i)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   rog(i) = cn(i,1) / phig(i)
   rop(i) = cn(i,6) / phip(i)
   tg(i) = pg(i) / (rog(i) * rrg)
   tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0_dp) * rop(i) * cpp)
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ag(i) = sqrt((gamag * pg(i)) / rog(i))
   ap(i) = sqrt((gamap * (pp(i) + pinf)) / rop(i))

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Ensure minimum values
   if (phig(i) < epslnmin) then
    phig(i) = epslnmin
    phip(i) = 1.0_dp - phig(i)
    ug(i) = up(i)
    vg(i) = vp(i)
    wg(i) = wp(i)
    tg(i) = tp(i)
    pg(i) = pp(i)
   end if

   if (phip(i) < epslnmin) then
    phip(i) = epslnmin
    phig(i) = 1.0_dp - phip(i)
    up(i) = ug(i)
    vp(i) = vg(i)
    wp(i) = wg(i)
    tp(i) = tg(i)
    pp(i) = pg(i)
   end if

   ! Smooth transition between phases
   if (phig(i) >= epslnmin .and. phig(i) <= epslnmax) then
    psig = (phig(i) - epslnmin) / (epslnmax - epslnmin)
    functiong = -psig * psig * (2.0_dp * psig - 3.0_dp)
    ug(i) = functiong * ug(i) + (1.0_dp - functiong) * up(i)
    vg(i) = functiong * vg(i) + (1.0_dp - functiong) * vp(i)
    wg(i) = functiong * wg(i) + (1.0_dp - functiong) * wp(i)
    tg(i) = functiong * tg(i) + (1.0_dp - functiong) * tp(i)
    pg(i) = functiong * pg(i) + (1.0_dp - functiong) * pp(i)
   end if

   if (phip(i) >= epslnmin .and. phip(i) <= epslnmax) then
    psig = (phip(i) - epslnmin) / (epslnmax - epslnmin)
    functiong = -psig * psig * (2.0_dp * psig - 3.0_dp)
    up(i) = functiong * up(i) + (1.0_dp - functiong) * ug(i)
    vp(i) = functiong * vp(i) + (1.0_dp - functiong) * vg(i)
    wp(i) = functiong * wp(i) + (1.0_dp - functiong) * wg(i)
    tp(i) = functiong * tp(i) + (1.0_dp - functiong) * tg(i)
    pp(i) = functiong * pp(i) + (1.0_dp - functiong) * pg(i)
   end if

   rog(i) = pg(i) / (rrg * tg(i))
   rop(i) = (pp(i) + pinf) * gamap / ((gamap - 1.0_dp) * cpp * tp(i))
   eg = pg(i) / ((gamag - 1.0_dp) * rog(i))
   ep = (pp(i) + gamap * pinf) / ((gamap - 1.0_dp) * rop(i))
   egtotal(i) = eg + 0.5_dp * (ug(i)**2 + vg(i)**2 + wg(i)**2)
   eptotal(i) = ep + 0.5_dp * (up(i)**2 + vp(i)**2 + wp(i)**2)
   
   ! Update cn matrix
   cn(i,1) = phig(i) * rog(i)
   cn(i,2) = phig(i) * rog(i) * ug(i)
   cn(i,3) = phig(i) * rog(i) * vg(i)
   cn(i,4) = phig(i) * rog(i) * wg(i)
   cn(i,5) = phig(i) * rog(i) * egtotal(i)
   cn(i,6) = phip(i) * rop(i)
   cn(i,7) = phip(i) * rop(i) * up(i)
   cn(i,8) = phip(i) * rop(i) * vp(i)
   cn(i,9) = phip(i) * rop(i) * wp(i)
   cn(i,10) = phip(i) * rop(i) * eptotal(i)
#ifdef KE_TURB
        cn(i,nkg) = phig(i) * rog(i) * tkg(i)
        cn(i,neg) = phig(i) * rog(i) * teg(i)
#endif 

#ifdef NP_P
        cn(i,npp) = enp(i)
#endif

#ifdef YP_P
        cn(i,nyfg) = phig(i) * rog(i) * yfg(i)
	cn(i,nyog) = phig(i) * rog(i) * yog(i)
	cn(i,nypg) = phig(i) * rog(i) * ypg(i)
#endif
  !end do
end do
!$acc end parallel loop 

!=======================
!Boundary condition loop
!=======================

!=========================
!Subsonic Inflow Condition
!=========================
!$acc parallel loop present(bc2_list(:),n_bc2,phip(:),phipzero,phig(:),ug(:),vg(:),wg(:),ugmin,vgmin,wgmin, &
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,egtotal,eptotal, &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),nparent(:),nghost(:))  &
!$acc private(ii,np,ng,eg,ep,i)

do ii = 1, n_bc2
  i = bc2_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i)
  
#ifdef FIVE_PHASE

    phiq(ng)= epslnmin
    phir(ng)= epslnmin
    phis(ng)= epslnmin
    phip(ng) = phipzero
    phig(ng) = 1.0d0 - phip(ng) - phiq(ng) - phir(ng)- phis(ng)
    
#else   
    phip(ng) = phipzero
    phig(ng) = 1.0d0 - phip(ng)
    
#endif 

#ifdef NP_P
    enp(ng)=enp(np)   
#endif

#ifdef YP_P
     yfg(ng)=0.5d0
     yog(ng)=0.5d0
     ypg(ng)=dmax1(1.0d0-yfg(ng)-yog(ng),tolsp)
	
     cpg=yfg(ng)*cpfg+yog(ng)*cpog+ypg(ng)*cppg
     rrg=8314.0d0*(yfg(ng)/amolfg+yog(ng)/amolog+ypg(ng)/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
     CMUG = 11.848e-8 * DSQRT(amolg)
#endif  
     
    ug(ng) = ugmin
    vg(ng) = vgmin
    wg(ng) = wgmin
    tg(ng) = tgmin
    pg(ng) = pg(np)
    rog(ng) = pg(ng) / (rrg * tg(ng))
    
    up(ng) = 0.99d0 * ugmin
    vp(ng) = 0.99d0 * vgmin
    wp(ng) = 0.99d0 * wgmin
    tp(ng) = tgmin
    pp(ng) = pp(np)
    rop(ng) = gamap * (pg(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))
    
    eg = pg(ng) / (rog(ng) * (gamag - 1.0d0))
    ep = (pp(ng) + gamap * pinf) / (rop(ng) * (gamap - 1.0d0))
    egtotal(ng) = eg + 0.5d0 * (ug(ng)**2 + vg(ng)**2 + wg(ng)**2)
    eptotal(ng) = ep + 0.5d0 * (up(ng)**2 + vp(ng)**2 + wp(ng)**2)
  
#ifdef KE_TURB
     amug =  as*dsqrt(tg(ng))/(1.0d0+ ts/tg(ng))
     tkg(ng)=0.666d0*0.01d0*(ug(ng)*ug(ng)+vg(ng)*vg(ng))
     teg(ng)=0.0845d0*rog(ng)*tkg(ng)**2.0d0/(amug*0.1d0)
#endif	     
end do
!$acc end parallel loop

!=========================
!Subsonic Outflow Condition
!=========================
!$acc parallel loop present(bc3_list(:),n_bc3,phip(:),phipzero,phig(:),ug(:),vg(:),wg(:),ugmin,vgmin,wgmin, &
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,egtotal,eptotal,dd(:), &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),x(:),nod(:,:),nparent(:),nghost(:),xmax,total,pmin,dt)  &
!$acc private(i,ii,np,ng,eg,ep,n1,n2,n3,n4,x1,x2,x4,x3,ep,eg, &
!$acc fract,d,alphat,wt,fracu,alphau,xcenter,wu,amug,cvg)

do ii = 1, n_bc3
  i = bc3_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i)

#ifdef NP_P
     enp(ng)=enp(np)
#endif        


#ifdef YP_P
     yfg(ng)=yfg(np)
     yog(ng)=yog(np)
     ypg(ng)=ypg(np)
	
     cpg=yfg(ng)*cpfg+yog(ng)*cpog+ypg(ng)*cppg
     rrg=8314.0d0*(yfg(ng)/amolfg+yog(ng)/amolog+ypg(ng)/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
     CMUG = 11.848e-8 * DSQRT(amolg)
#endif
         phig(ng)=phig(np)
         phip(ng)=phip(np)

         wT = (ug(np))
         wU = (ug(np)) + sqrt(gamag*(pg(np))/(rog(np)))
         
         alphaT = wT*dt/dd(np)  
         alphaU = wU*dt/dd(np)
        
         fracT = 1/(1+ alphaT)
         fracU = 1/(1+ alphaU) 
         
         if(total.eq.1) then
           ug(ng)=ug(np)
           vg(ng)=vg(np)
           wg(ng)=wg(np)
           tg(ng)=tg(np)  
           pg(ng)=pmin
         endif
         
         if(ug(np).gt.ug(ng))then                                   
           ug(ng)=ug(ng)* fracU + (1- fracU)* (2.0d0*ug(np) - ug(ng))
           vg(ng)=vg(ng)* fracT + (1- fracT)* (2.0d0*vg(np) - vg(ng))
           wg(ng)=wg(ng)* fracT + (1- fracT)* (2.0d0*wg(np) - wg(ng))
           tg(ng)=tg(ng)* fracT + (1- fracT)* (2.0d0*tg(np) - tg(ng))
        else
           ug(ng)=ug(np)
           vg(ng)=vg(np)
           wg(ng)=wg(np)
           tg(ng)=tg(np)
        endif   
            pg(ng)=pmin
           
         rog(ng)=pg(ng)/(rrg*tg(ng))
                 
         vp(ng)=vp(np)
         up(ng)=up(np)
         wp(ng)=wp(np)
	 pp(ng)=pmin
	 tp(ng)=tp(np)   
    	 rop(ng)=gamap*(pp(ng)+pinf)/((gamap-1.0d0)*cpp*tp(ng))
    eg = rrg * tg(ng) / (gamag - 1.0d0)
    ep = cpp * tp(ng) / gamap + pinf / rop(ng)
    egtotal(ng) = eg + 0.5d0 * (ug(ng)**2 + vg(ng)**2 + wg(ng)**2)
    eptotal(ng) = ep + 0.5d0 * (up(ng)**2 + vp(ng)**2 + wp(ng)**2)
   
#ifdef KE_TURB
        tkg(ng)=tkg(np)
	teg(ng)=teg(np)
#endif 
end do      
!$acc end parallel loop


!=========================
!Wall Boundary Condition
!=========================
!$acc parallel loop present(bc5_list(:),n_bc5,phip(:),phipzero,phig(:),ug(:),vg(:),wg(:),ugmin,vgmin,wgmin, &
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,egtotal,eptotal, &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),nparent(:),nghost(:))  &
!$acc private(ii,np,ng,eg,ep,i)

do ii = 1, n_bc5
  i = bc5_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i)
   phig(ng) = phig(np)
    ug(ng) = -ug(np)
    vg(ng) = -vg(np)
    wg(ng) = -wg(np)
    tg(ng) = tg(np)
    
    phip(ng) = phip(np)
    up(ng) = -up(np)
    vp(ng) = -vp(np)
    wp(ng) = -wp(np)
    tp(ng) = tp(np)
    pg(ng) = pg(np)
    pp(ng) = pp(np)
    
#ifdef NP_P
     enp(ng)=enp(np)
#endif

#ifdef YP_P
      yfg(ng)=yfg(np)
      yog(ng)=yog(np)
      ypg(ng)=ypg(np)
	
     cpg=yfg(ng)*cpfg+yog(ng)*cpog+ypg(ng)*cppg
     rrg=8314.0d0*(yfg(ng)/amolfg+yog(ng)/amolog+ypg(ng)/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
     CMUG = 11.848e-8 * DSQRT(amolg)
#endif	
       

    rog(ng) = pg(ng) / (rrg * tg(ng))
    rop(ng) = gamap * (pp(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))
    eg = rrg * tg(ng) / (gamag - 1.0d0)
    ep = cpp * tp(ng) / gamap + pinf / rop(ng)
    egtotal(ng) = eg + 0.5d0 * (ug(ng)**2 + vg(ng)**2 + wg(ng)**2)
    eptotal(ng) = ep + 0.5d0 * (up(ng)**2 + vp(ng)**2 + wp(ng)**2)
    
#ifdef KE_TURB
    tkg(ng)=2.0d0*tolkg-tkg(np)
    teg(ng)=teg(np)
#endif  
end do
!$acc end parallel loop


!=================================
!Top and botoom boundary Condition
!=================================
!$acc parallel loop present(bc51_list(:),n_bc51,phip(:),phipzero,phig(:),ug(:),vg(:),wg(:),ugmin,vgmin,wgmin, &
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,egtotal,eptotal,dd(:), &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),y(:),nod(:,:),nparent(:),nghost(:),ymax,ymin,total,pmin,dt)&
!$acc private(i,ii,np,ng,eg,ep,n1,n2,n3,n4,y1,y2,y4,y3,ep,eg, &
!$acc fract,d,alphat,wt,fracu,alphau,ycenter,wu,amug,cvg)
do ii = 1, n_bc51
  i = bc51_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i)

         wT = (vg(np))
         wU = (vg(np)) + sqrt(gamag*(pg(np))/(rog(np)))
         
         alphaT = wT*dt/dd(np)  
         alphaU = wU*dt/dd(np)
        
         fracT = 1/(1+ alphaT)
         fracU = 1/(1+ alphaU) 
         
!         if(total.eq.1) then
!           ug(ng)=ug(np)
!           vg(ng)=-vg(np)
!           wg(ng)=wg(np)
!           tg(ng)=tg(np)  
         !  pg(ng)=pmin
!         endif
         
         if(vg(np).gt.vg(ng))then                                   
           vg(ng)=-vg(np)!* fracU + (1- fracU)* (2.0d0*vg(np) - vg(ng))
           ug(ng)=ug(ng)* fracT + (1- fracT)* (2.0d0*ug(np) - ug(ng))
           wg(ng)=wg(ng)* fracT + (1- fracT)* (2.0d0*wg(np) - wg(ng))
           tg(ng)=tg(ng)* fracT + (1- fracT)* (2.0d0*tg(np) - tg(ng))
        else
           ug(ng)=ug(np)
           vg(ng)=-vg(np)
           wg(ng)=wg(np)
           tg(ng)=tg(np)
        endif   
          pg(ng)=pg(np)
 !WT      
 
 !       phig(ng) = phig(np)
 !       ug(ng) = ug(np)
 !       vg(ng) = -vg(np)
 !       wg(ng) = wg(np)
 !       tg(ng) = tg(np)
    phig(ng) = phig(np)
 !   ug(ng) = ug(np)
 !   vg(ng) = -vg(np)
 !   wg(ng) = wg(np)
 !   tg(ng) = tg(np)
    phip(ng) = phip(np)
    
    up(ng) = up(np)
    vp(ng) = -vp(np)
    wp(ng) = wp(np) 
    tp(ng) = tp(np)
    pg(ng) = pg(np)
    pp(ng) = pp(np)
    
#ifdef NP_P
    enp(ng)=enp(np)
#endif   

#ifdef YP_P
     yfg(ng)=yfg(np)
     yog(ng)=yog(np)
     ypg(ng)=ypg(np)
	
     cpg=yfg(ng)*cpfg+yog(ng)*cpog+ypg(ng)*cppg
     rrg=8314.0d0*(yfg(ng)/amolfg+yog(ng)/amolog+ypg(ng)/amolpg)
     amolg=8314.0d0/rrg
     cvg=cpg-rrg
     gamag=cpg/cvg
     CMUG = 11.848e-8 * DSQRT(amolg)
#endif
           
    rog(ng) = pg(ng) / (rrg * tg(ng))
    rop(ng) = gamap * (pp(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))
    eg = rrg * tg(ng) / (gamag - 1.0d0)
    ep = cpp * tp(ng) / gamap + pinf / rop(ng)
    egtotal(ng) = eg + 0.5d0 * (ug(ng)**2 + vg(ng)**2 + wg(ng)**2)
    eptotal(ng) = ep + 0.5d0 * (up(ng)**2 + vp(ng)**2 + wp(ng)**2)
#ifdef KE_TURB
         tkg(ng)=tkg(np)
	 teg(ng)=teg(np)
#endif	
end do 
!$acc end parallel loop
   
!=================================
!Front and back boundary Condition
!=================================
!$acc parallel loop present(bc53_list(:),n_bc53,phip(:),phipzero,phig(:),ug(:),vg(:),wg(:),ugmin,vgmin,wgmin, &
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,egtotal,eptotal,dd(:), &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),z(:),nod(:,:),nparent(:),nghost(:),zmax,zmin,total,pmin,dt)&
!$acc private(i,ii,np,ng,eg,ep,n1,n2,n3,n4,z1,z2,z4,z3,ep,eg, &
!$acc fract,d,alphat,wt,fracu,alphau,zcenter,wu,amug,cvg)
do ii = 1, n_bc53
  i = bc53_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i) 
       
         wT = (wg(np))
         wU = (wg(np)) + sqrt(gamag*(pg(np))/(rog(np)))
         
         alphaT = wT*dt/dd(np)  
         alphaU = wU*dt/dd(np)
        
         fracT = 1/(1+ alphaT)
         fracU = 1/(1+ alphaU) 
         
      !   if(total.eq.1) then
      !     ug(ng)=ug(np)
      !     vg(ng)=vg(np)
      !     wg(ng)=-wg(np)
      !     tg(ng)=tg(np)  
         !  pg(ng)=pmin
      !   endif
         
         if(wg(np).gt.wg(ng))then                                   
           wg(ng)=-wg(np)!* fracU + (1- fracU)* (2.0d0*wg(np) - wg(ng))
           ug(ng)=ug(ng)* fracT + (1- fracT)* (2.0d0*ug(np) - ug(ng))
           vg(ng)=vg(ng)* fracT + (1- fracT)* (2.0d0*vg(np) - vg(ng))
           tg(ng)=tg(ng)* fracT + (1- fracT)* (2.0d0*tg(np) - tg(ng))
        else
           ug(ng)=ug(np)
           vg(ng)=vg(np)
           wg(ng)=-wg(np)
           tg(ng)=tg(np)
        endif   
          pg(ng)=pg(np)
 !WT  
      
   !     phig(ng) = phig(np)
   !     ug(ng) = ug(np)
   !     vg(ng) = vg(np)
   !     wg(ng) = -wg(np)
   !     tg(ng) = tg(np)
   
    phig(ng) = phig(np)
 !   ug(ng) = ug(np)
 !   vg(ng) = vg(np)
 !   wg(ng) = -wg(np)
 !   tg(ng) = tg(np)
    phip(ng) = phip(np)
    
    up(ng) = up(np)
    vp(ng) = vp(np)
    wp(ng) = -wp(np) 
    tp(ng) = tp(np)
    pg(ng) = pg(np)
    pp(ng) = pp(np)
    
#ifdef NP_P
    enp(ng)=enp(np)
#endif

    rog(ng) = pg(ng) / (rrg * tg(ng))
    rop(ng) = gamap * (pp(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))
    eg = rrg * tg(ng) / (gamag - 1.0d0)
    ep = cpp * tp(ng) / gamap + pinf / rop(ng)
    egtotal(ng) = eg + 0.5d0 * (ug(ng)**2 + vg(ng)**2 + wg(ng)**2)
    eptotal(ng) = ep + 0.5d0 * (up(ng)**2 + vp(ng)**2 + wp(ng)**2)
#ifdef KE_TURB
        tkg(ng)=tkg(np)
	teg(ng)=teg(np)
#endif	
end do 
!$acc end parallel loop
   



!$acc parallel loop present(nghosts,rog(:),vg(:),phig(:),nparent(:),vp(:),nghost(:),cn(:,:),up(:), &
!$acc rop(:),pp(:),ug(:),wp(:),wg(:),phip(:),ntype(:),nod(:,:),neq,ymax,ymin,zmin,zmax,  &
!$acc eptotal,egtotal,tkg(:),teg(:),neg,nkg,enp(:),npp,yfg(:),yog(:),ypg(:),nyfg,nyog,nypg) &
!$acc private(ng,np,i)

  do i = 1, nghosts
   np=nparent(i)
   ng=nghost(i)
    
   cn(ng,1) = phig(ng) * rog(ng)
   cn(ng,2) = phig(ng) * rog(ng) * ug(ng)
   cn(ng,3) = phig(ng) * rog(ng) * vg(ng)
   cn(ng,4) = phig(ng) * rog(ng) * wg(ng)
   cn(ng,5) = phig(ng) * rog(ng) * egtotal(ng)
   cn(ng,6) = phip(ng) * rop(ng)
   cn(ng,7) = phip(ng) * rop(ng) * up(ng)
   cn(ng,8) = phip(ng) * rop(ng) * vp(ng)
   cn(ng,9) = phip(ng) * rop(ng) * wp(ng)
   cn(ng,10) = phip(ng) * rop(ng) * eptotal(ng)
#ifdef KE_TURB
        cn(ng,nkg) = phig(ng) * rog(ng) * tkg(ng)
        cn(ng,neg) = phig(ng) * rog(ng) * teg(ng)
#endif

#ifdef NP_P
        cn(ng,npp) = enp(ng)
#endif

#ifdef YP_P
        cn(ng,nyfg) = phig(ng) * rog(ng) * yfg(ng)
	cn(ng,nyog) = phig(ng) * rog(ng) * yog(ng)
	cn(ng,nypg) = phig(ng) * rog(ng) * ypg(ng)
#endif

  end do
!$acc end parallel loop  

    return
end subroutine primitiveplusbc

