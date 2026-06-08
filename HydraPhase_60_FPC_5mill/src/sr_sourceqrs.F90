!=======================================================================! 
!About this subroutine							!
!This subroutine is for source term only when 5 phases will be there 	!
!This subroutine will be called from jamesonqrs only			!		
!=======================================================================!
subroutine SOURCEQRS
      USE GlobalVariables
      USE MatrixOps
      IMPLICIT NONE

   INTEGER :: i
   REAL(DP) :: voli, reynq,relvq, cdq, dragq, reynr, relvr, cdr, dragr
   REAL(DP) :: reyns, relvs, cds, drags, amup,xi,xai,xbi,xci,xdi,yi,yai
   REAL(DP) :: ybi,yci,ydi,zi,zai,zbi,zci,zdi,upi,upai,upbi,upci,updi
   REAL(DP) :: vpi,vpai,vpbi,vpci,vpdi,wpi,wpai,wpbi,wpci,wpdi
   REAL(DP) :: wt1,wt2,wt3,wt4,a1,a2,a3,b1,b2,b3,cup1,cup2,cup3
   REAL(DP) :: cvp1,cvp2,cvp3,cwp1,cwp2,cwp3,dnr,dupx,dupy
   REAL(DP) :: dupz,dvpx,dvpy,dvpz,dwpx,dwpy,dwpz
   REAL(DP) :: gpart1,gpart2,gamapdot,amug,akp
   
!$acc parallel loop present(neles,us(:),vs(:),wp(:),vol(:),vp(:),rhs(:,:),rhoq,rhor, &
!$acc up(:),phis(:),wq(:),vq(:),uq(:),phiq(:),ws(:),wr(:),vr(:),ur(:),phir(:),gravity,rhos, &
!$acc rads,radq,radr,vols,volr,volq,tg(:),cpg,prg,as,ts,wg(:),vg(:),ug(:), &
!$acc nc4(:),zcel(:),ycel(:),xcel(:),nc3(:),nc2(:),nc1(:),phip(:)) &
!$acc private(reyns,relvs,voli,cds,drags,reynq,relvq,cdq,dragq, &
!$acc reynr,relvr,cdr,dragr,amup,xi,xai,xbi,xci,xdi,yi,yai,ybi,yci,ydi, &
!$acc zi,zai,zbi,zci,zdi,upi,upai,upbi,upci,updi, &
!$acc vpi,vpai,vpbi,vpci,vpdi,wpi,wpai,wpbi,wpci, &
!$acc wpdi,wt1,wt2,wt3,wt4,a1,a2,a3,b1,b2,b3, &
!$acc cup1,cup2,cup3,cvp1,cvp2,cvp3,cwp1,cwp2,cwp3,amug, &
!$acc dnr,dupx,dupy,dupz,dvpx,dvpy,dvpz,dwpx,dwpy,dwpz,gpart1,gamapdot,& 
!$acc gpart2)

do i=1,neles
	 
    voli=vol(i)
	 
    amup = 0.001d0 
    amug =  as*dsqrt(tg(i))/(1.0d0+ ts/tg(i))

#ifndef AMUP_CONSTANT 
        xi = xcel(i)
        xai = xcel(nc1(i)) - xi
        xbi = xcel(nc2(i)) - xi
        xci = xcel(nc3(i)) - xi
        xdi = xcel(nc4(i)) - xi

        yi = ycel(i)
        yai = ycel(nc1(i)) - yi
        ybi = ycel(nc2(i)) - yi
        yci = ycel(nc3(i)) - yi
        ydi = ycel(nc4(i)) - yi

        zi = zcel(i)
        zai = zcel(nc1(i)) - zi
        zbi = zcel(nc2(i)) - zi
        zci = zcel(nc3(i)) - zi
        zdi = zcel(nc4(i)) - zi

        upi  = up(i)
        upai = up(nc1(i)) - upi
        upbi = up(nc2(i)) - upi
        upci = up(nc3(i)) - upi
        updi = up(nc4(i)) - upi
		
	vpi  = vp(i)
        vpai = vp(nc1(i)) - vpi
        vpbi = vp(nc2(i)) - vpi
        vpci = vp(nc3(i)) - vpi
        vpdi = vp(nc4(i)) - vpi
		
	wpi  = wp(i)
        wpai = wp(nc1(i)) - wpi
        wpbi = wp(nc2(i)) - wpi
        wpci = wp(nc3(i)) - wpi
        wpdi = wp(nc4(i)) - wpi

        wt1 = SQRT((xai)**2 + (yai)**2 + (zai)**2)
        wt2 = SQRT((xbi)**2 + (ybi)**2 + (zbi)**2)
        wt3 = SQRT((xci)**2 + (yci)**2 + (zci)**2)
        wt4 = SQRT((xdi)**2 + (ydi)**2 + (zdi)**2)
		
	a1 = (xai)**2 / wt1 + (xbi)**2 / wt2 + (xci)**2 / wt3 + (xdi)**2 / wt4
        a2 = (yai)**2 / wt1 + (ybi)**2 / wt2 + (yci)**2 / wt3 + (ydi)**2 / wt4
        a3 = (zai)**2 / wt1 + (zbi)**2 / wt2 + (zci)**2 / wt3 + (zdi)**2 / wt4
        
        b1 = (xai) * (yai) / wt1 + (xbi) * (ybi) / wt2 + (xci) * (yci) / wt3 + (xdi) * (ydi) / wt4
        b2 = (xai) * (zai) / wt1 + (xbi) * (zbi) / wt2 + (xci) * (zci) / wt3 + (xdi) * (zdi) / wt4
        b3 = (yai) * (zai) / wt1 + (ybi) * (zbi) / wt2 + (yci) * (zci) / wt3 + (ydi) * (zdi) / wt4

        cup1 = (upai) * (xai) / wt1 + (upbi) * (xbi) / wt2 + (upci) * (xci) / wt3 + (updi) * (xdi) / wt4
        cup2 = (upai) * (yai) / wt1 + (upbi) * (ybi) / wt2 + (upci) * (yci) / wt3 + (updi) * (ydi) / wt4
        cup3 = (upai) * (zai) / wt1 + (upbi) * (zbi) / wt2 + (upci) * (zci) / wt3 + (updi) * (zdi) / wt4
		
	cvp1 = (vpai) * (xai) / wt1 + (vpbi) * (xbi) / wt2 + (vpci) * (xci) / wt3 + (vpdi) * (xdi) / wt4
        cvp2 = (vpai) * (yai) / wt1 + (vpbi) * (ybi) / wt2 + (vpci) * (yci) / wt3 + (vpdi) * (ydi) / wt4
        cvp3 = (vpai) * (zai) / wt1 + (vpbi) * (zbi) / wt2 + (vpci) * (zci) / wt3 + (vpdi) * (zdi) / wt4
		
	cwp1 = (wpai) * (xai) / wt1 + (wpbi) * (xbi) / wt2 + (wpci) * (xci) / wt3 + (wpdi) * (xdi) / wt4
        cwp2 = (wpai) * (yai) / wt1 + (wpbi) * (ybi) / wt2 + (wpci) * (yci) / wt3 + (wpdi) * (ydi) / wt4
        cwp3 = (wpai) * (zai) / wt1 + (wpbi) * (zbi) / wt2 + (wpci) * (zci) / wt3 + (wpdi) * (zdi) / wt4

        dnr = det(a1, b1, b2, b1, a2, b3, b2, b3, a3)
        dupx = det(cup1, cup2, cup3, b1, a2, b3, b2, b3, a3) / dnr
        dupy = det(a1, b1, b2, cup1, cup2, cup3, b2, b3, a3) / dnr
        dupz = det(a1, b1, b2, b1, a2, b3, cup1, cup2, cup3) / dnr
		
	dvpx = det(cvp1, cvp2, cvp3, b1, a2, b3, b2, b3, a3) / dnr
        dvpy = det(a1, b1, b2, cvp1, cvp2, cvp3, b2, b3, a3) / dnr
        dvpz = det(a1, b1, b2, b1, a2, b3, cvp1, cvp2, cvp3) / dnr
		
	dwpx = det(cwp1, cwp2, cwp3, b1, a2, b3, b2, b3, a3) / dnr
        dwpy = det(a1, b1, b2, cwp1, cwp2, cwp3, b2, b3, a3) / dnr
        dwpz = det(a1, b1, b2, b1, a2, b3, cwp1, cwp2, cwp3) / dnr

!--------------------------------------------------------------------------------------------------------------
        if(phip(i).gt.0.01d0)then
        gpart1=2.0d0*(dupx**2+dvpy**2+dwpz**2)
	gpart2=(dupy+dvpx)**2+(dupz+dwpx)**2+(dvpz+dwpy)**2
        gamapdot=min(max(sqrt(gpart1+gpart2),0.0001d0),10.0d0)	 

           if(gamapdot.lt.0.01d0)then
	     amup=2658.83d0*((2.0d0-0.814d0)+(0.814d0-1.0d0)*(gamapdot/0.01d0))+10.7d0*(2.0d0-(gamapdot/0.01d0))/0.01d0 
	   else
	     amup=10.7d0/gamapdot+2658.83d0*(gamapdot/0.01d0)**(0.814d0-1.0d0)
	   endif
        
 !       if(amup.lt.100) then 
!	    print*, "Warning low amup detected in cell id", i,", and its value is", gamapdot, amup    
!	 STOP
!	 endif  
	
! print*, "amup", i, amup, gamapdot	
	else
	   amup=amug
	endif
#endif	
!--------------------------------------------------------------------------------------------------------------	 
	 
	REYNQ = dmax1(rhoq*2.0d0*radq*dsqrt((up(i)-uq(i))**2+(vp(i)-vq(i))**2+(wp(i)-wq(i))**2)/amup,0.001d0)
	RELVQ = dsqrt((up(i)-uq(i))**2+(vp(i)-vq(i))**2+(wp(i)-wq(i))**2)

          IF(REYNQ.LT.0.1d0)THEN
	      	CDQ = 24.0d0/REYNQ 
          ELSEIF(REYNQ.LT.1000.0d0)THEN
	       	CDQ = 0.48d0 + 28.0d0/(REYNQ**0.85d0) 
          ELSE
 	      	CDQ = 0.48D0
          ENDIF

	  DRAGQ=0.5D0*CDQ*rhoq*3.1415926D0*RADQ*RADQ*RELVQ*(phiq(i)/volq)
	  
	  rhs(i,12)=rhs(i,12)-dragq*(up(i)-uq(i))*voli
	  rhs(i,13)=rhs(i,13)-dragq*(vp(i)-vq(i))*voli-phiq(i)*rhoq*gravity*voli
	  rhs(i,14)=rhs(i,14)-dragq*(wp(i)-wq(i))*voli
	  
	  rhs(i,7)=rhs(i,7)+dragq*(up(i)-uq(i))*voli
	  rhs(i,8)=rhs(i,8)+dragq*(vp(i)-vq(i))*voli
	  rhs(i,9)=rhs(i,9)+dragq*(wp(i)-wq(i))*voli
	  
!   -------------------------------------------------------------------------------
        REYNR = dmax1(rhor*2.0d0*radr*dsqrt((up(i)-ur(i))**2+(vp(i)-vr(i))**2+(wp(i)-wr(i))**2)/amup,0.001d0)
        RELVR = dsqrt((up(i)-ur(i))**2+(vp(i)-vr(i))**2+(wp(i)-wr(i))**2)

          IF(REYNR.LT.0.1d0)THEN
	      	CDR = 24.0d0/REYNR 
          ELSEIF(REYNR.LT.1000.0d0)THEN
	       	CDR = 0.48d0 + 28.0d0/(REYNR**0.85d0) 
          ELSE
 	      	CDR = 0.48D0
          ENDIF

	  DRAGR=0.5D0*CDR*rhor*3.1415926D0*RADR*RADR*RELVR*(phir(i)/volr)
	  
	  rhs(i,16)=rhs(i,16)-dragr*(up(i)-ur(i))*voli
	  rhs(i,17)=rhs(i,17)-dragr*(vp(i)-vr(i))*voli-phir(i)*rhor*gravity*voli
	  rhs(i,18)=rhs(i,18)-dragr*(wp(i)-wr(i))*voli
	  
	  rhs(i,7)=rhs(i,7)+dragr*(up(i)-ur(i))*voli
	  rhs(i,8)=rhs(i,8)+dragr*(vp(i)-vr(i))*voli	
	  rhs(i,9)=rhs(i,9)+dragr*(wp(i)-wr(i))*voli	

!   -------------------------------------------------------------------------------
        
        REYNS = dmax1(rhos*2.0d0*rads*dsqrt((up(i)-us(i))**2+(vp(i)-vs(i))**2+(wp(i)-ws(i))**2)/amup,0.001d0)
        RELVS = dsqrt((up(i)-us(i))**2+(vp(i)-vs(i))**2+(wp(i)-ws(i))**2)

          IF(REYNS.LT.0.1d0)THEN
	      	CDS = 24.0d0/REYNS 
          ELSEIF(REYNS.LT.1000.0d0)THEN
	       	CDS = 0.48d0 + 28.0d0/(REYNS**0.85d0) 
          ELSE
 	      	CDS = 0.48D0
          ENDIF

	  DRAGS=0.5D0*CDS*rhos*3.1415926D0*RADS*RADS*RELVS*(phis(i)/vols)
	  
	  rhs(i,20)=rhs(i,20)-drags*(up(i)-us(i))*voli
	  rhs(i,21)=rhs(i,21)-drags*(vp(i)-vs(i))*voli-phis(i)*rhos*gravity*voli
	  rhs(i,22)=rhs(i,22)-drags*(wp(i)-ws(i))*voli
	  
	  rhs(i,7)=rhs(i,7)+drags*(up(i)-us(i))*voli
	  rhs(i,8)=rhs(i,8)+drags*(vp(i)-vs(i))*voli	  
	  rhs(i,9)=rhs(i,9)+drags*(wp(i)-ws(i))*voli	  
	 enddo
!$acc end parallel loop	 

end subroutine sourceqrs
