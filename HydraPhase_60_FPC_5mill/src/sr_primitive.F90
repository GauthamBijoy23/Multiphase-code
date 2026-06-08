SUBROUTINE primitive
    USE GlobalVariables
    IMPLICIT NONE
    external geometry
    
    INTEGER :: i,k,startl,endl
    REAL(DP) :: eg,ep,bbb,ccc,psig,functiong,aaa,cvg
    logical :: active
    
!$acc parallel loop & 
!$acc present(tg(:),rog(:),phig(:),phip(:),wp(:),wg(:),nkg,neg, &
!$acc ug(:),up(:),tp(:),vp(:),vg(:),pg(:),pp(:),rop(:),egtotal,eptotal, &
!$acc cpp,gamag,gamap,rrg,epslnmax,epslnmin,foursigmabyd,tkg(:),teg(:), &
!$acc pinf,neles,cn(:,:),nghosts,ntot,rhoq,rhos,rhor,phiq(:),phir(:),phis(:),enp(:),npp, &
!$acc yfg(:),yog(:),ypg(:),cpg,cpfg,cpog,cppg,rrg,amolfg,amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg)  &
!$acc private(functiong,bbb,psig,ccc,ep,eg,i,startl,endl,aaa,cvg,active) 

do i = 1, ntot

   active = ( i <= neles ) .or. ( i > neles + nghosts )
   if (.not. active) cycle

        ug(i) = cn(i,2) / cn(i,1)
        vg(i) = cn(i,3) / cn(i,1)
        wg(i) = cn(i,4) / cn(i,1)
        up(i) = cn(i,7) / cn(i,6)
        vp(i) = cn(i,8) / cn(i,6)
        wp(i) = cn(i,9) / cn(i,6)
           
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

!if(myid==1.and.i==5) print*,"prim",cu(i,1),cu(i,nkg),cn(i,1),cn(i,nkg)   
        ! Compute internal energies
        eg = cn(i,5) / cn(i,1) - 0.5d0 * (ug(i)**2 + vg(i)**2 + wg(i)**2)
        ep = cn(i,10) / cn(i,6) - 0.5d0 * (up(i)**2 + vp(i)**2 + wp(i)**2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef FIVE_PHASE
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
   tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0_dp) * rop(i) * cpp)
  
   ! Compute phase fractions
   phig(i) = dmin1(dmax1(cn(i,1)/rog(i),epslnmin),(1-epslnmin))
   phip(i) = dmin1(dmax1(cn(i,6)/rop(i),epslnmin),(1-epslnmin)) 

!For 2 phase
#else
     
        
       ! Compute pressures
        bbb = (foursigmabyd + pinf * gamap) - (cn(i,1) * (gamag - 1.0d0) * eg) - (cn(i,6) * (gamap - 1.0d0) * ep)
        ccc = (cn(i,1) * (gamag - 1.0d0) * eg * (foursigmabyd + pinf * gamap))
        pg(i) = (-bbb + sqrt(bbb**2 + 4.0d0 * ccc)) * 0.5d0
        pp(i) = pg(i) + foursigmabyd

        ! Compute phase fractions
        phig(i) = (cn(i,1) * (gamag - 1.0d0) * eg) / pg(i)
        phip(i) = 1.0d0 - phig(i)

        ! Compute densities and temperatures
        rog(i) = cn(i,1) / phig(i)
        rop(i) = cn(i,6) / phip(i)
        tg(i) = pg(i) / (rog(i) * rrg)
        tp(i) = ((pp(i) + pinf) * gamap) / ((gamap - 1.0d0) * rop(i) * cpp)
        
#endif

        ! Phase fraction correction
        if (phig(i) < epslnmin) then
            phig(i) = epslnmin
            phip(i) = 1.0d0 - phig(i)
            ug(i) = up(i)
            vg(i) = vp(i)
            wg(i) = wp(i)
            tg(i) = tp(i)
            pg(i) = pp(i)
        end if

        if (phip(i) < epslnmin) then
            phip(i) = epslnmin
            phig(i) = 1.0d0 - phip(i)
            up(i) = ug(i)
            vp(i) = vg(i)
            wp(i) = wg(i)
            tp(i) = tg(i)
            pp(i) = pg(i)
        end if

        ! Apply smoothing function
        if (phig(i) >= epslnmin .and. phig(i) <= epslnmax) then
            psig = (phig(i) - epslnmin) / (epslnmax - epslnmin)
            functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
            ug(i) = functiong * ug(i) + (1.0d0 - functiong) * up(i)
            vg(i) = functiong * vg(i) + (1.0d0 - functiong) * vp(i)
            wg(i) = functiong * wg(i) + (1.0d0 - functiong) * wp(i)
            tg(i) = functiong * tg(i) + (1.0d0 - functiong) * tp(i)
            pg(i) = functiong * pg(i) + (1.0d0 - functiong) * pp(i)
        end if

        if (phip(i) >= epslnmin .and. phip(i) <= epslnmax) then
            psig = (phip(i) - epslnmin) / (epslnmax - epslnmin)
            functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
            up(i) = functiong * up(i) + (1.0d0 - functiong) * ug(i)
            vp(i) = functiong * vp(i) + (1.0d0 - functiong) * vg(i)
            wp(i) = functiong * wp(i) + (1.0d0 - functiong) * wg(i)
            tp(i) = functiong * tp(i) + (1.0d0 - functiong) * tg(i)
            pp(i) = functiong * pp(i) + (1.0d0 - functiong) * pg(i)
        end if

        ! Update densities and energies
        rog(i) = pg(i) / (rrg * tg(i))
        rop(i) = (pp(i) + pinf) * gamap / ((gamap - 1.0d0) * cpp * tp(i))
        eg = pg(i) / ((gamag - 1.0d0) * rog(i))
        ep = (pp(i) + gamap * pinf) / ((gamap - 1.0d0) * rop(i))
        egtotal(i) = eg + 0.5d0 * (ug(i)**2 + vg(i)**2 + wg(i)**2)
        eptotal(i) = ep + 0.5d0 * (up(i)**2 + vp(i)**2 + wp(i)**2)

        ! Store updated values
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
!     end do
end do
!$acc end parallel loop   

return
end subroutine primitive

