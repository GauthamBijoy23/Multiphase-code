!=======================================================================! 
!About this subroutine							!
!This subroutine is for source term only when 2 phases will be there 	!
!This subroutine will be called from hllcbatten and hllzien		!		
!=======================================================================!

SUBROUTINE source
    USE GlobalVariables
    USE MatrixOps
    IMPLICIT NONE
    external geometry
    external primitive
    external smooth
    
    INTEGER :: i
    REAL(DP) :: voli, fdgx, fdgy, fdgz, qsource, cd, akg,amug
    REAL(DP) :: drag, enusseltp, delpstar, pint, uint, vint, wint
    REAL(DP) :: xi, xai, xbi, xci, xdi, yi, yai, ybi, yci, ydi, zi, zai, zbi, zci, zdi
    REAL(DP) :: phigi, phigai, phigbi, phigci, phigdi
    REAL(DP) :: wt1, wt2, wt3, wt4, dnr, dphigx, dphigy, dphigz
    REAL(DP) :: a1, a2, a3, b1, b2, b3, c1, c2, c3
    REAL(DP) :: amutg,sss1,sss2,sss3,sssg,pkgrng,etagrng,c1estarg,amutgbyr,stkg,steg
    REAL(DP) :: yfghat,yoghat,ypghat,yghatmin,grankappag,grangamag,omegadotfg,omegadotog,omegadotpg
        
!$acc parallel loop present(neles,viswg(:),visug(:),visvp(:),vistp(:), &
!$acc phip(:),ycel(:),zcel(:),xcel(:),nc2(:),wp(:),dl(:),rog(:),pg(:), &
!$acc rhs(:,:),ug(:),tp(:),vol(:),up(:),wg(:),vp(:),rop(:),viswp(:), &
!$acc tg(:),phig(:),nc4(:),nc1(:),nc3(:),visup(:),vg(:),visvg(:), &
!$acc vistg(:),prg,cmug,cpg,gravity,as,ts,&
!$acc tkg(:),vistkg(:),visteg(:),delugx(:),delugz(:),delugy(:),delvgz(:),delvgy(:), &
!$acc delwgx(:),cn(:,nkg),delvgx(:),delwgy(:),delwgz(:),teg(:),nkg,neg, &
!$acc  yfg(:),yog(:),ypg(:),sratio,nyfg,nypg,nyog) &
!$acc private(vint,enusseltp,fdgx,fdgy,c3,cd,pint,a3,dphigz,drag, &
!$acc phigi,qsource,amug,uint,voli,xdi,wt4,xai,xbi,xci,ydi,xi,yai,ybi,yci, &
!$acc zi,zdi,yi,zai,zbi,zci,phigdi,phigai,phigbi,phigci,wint,wt1,wt2,akg, &
!$acc wt3,a1,a2,b1,b2,b3,c1,c2,delpstar,dnr,dphigx,dphigy,fdgz, &
!$acc amutg,sss1,sss2,sss3,sssg,pkgrng,etagrng,c1estarg,amutgbyr,stkg,steg, &
!$acc yfghat,yoghat,ypghat,yghatmin,grankappag,grangamag,omegadotfg,omegadotog,omegadotpg)    

!This is because of jameson that is not calling derivative

    DO i = 1, neles

IF(iblank(i)/=1)CYCLE !-----mod3 checking for iblank

        voli = vol(i)
        fdgx = 0.0d0
        fdgy = 0.0d0
        fdgz = 0.0d0
        qsource = 0.0d0
        cd = 0.0d0
        amug =  as*dsqrt(tg(i))/(1.0d0+ ts/tg(i))  ! Ensure 'as' is properly declared elsewhere
        akg = amug * cpg / prg

        drag = 0.44d0 * phig(i) * (phig(i) - 1.0d0) * rog(i) * SQRT((ug(i) - up(i))**2 + &
               (vg(i) - vp(i))**2 + (wg(i) - wp(i))**2) * dl(i)
        fdgx = drag * (ug(i) - up(i))
        fdgy = drag * (vg(i) - vp(i))
        fdgz = drag * (wg(i) - wp(i))

!========================================================================
! if activating NP_P; choose a drag model as per the problem to be solved
!========================================================================
 
        enusseltp = 2.0d0
        qsource = 6.0d0 * phip(i) * amug * cpg * enusseltp * (tp(i) - tg(i)) / (prg * dl(i)**2)

        delpstar = 2.0d0 * (phig(i) * rog(i) * phip(i) * rop(i)) / (phig(i) * rop(i) + phip(i) * rog(i)) * &
                   ((ug(i) - up(i))**2 + (vg(i) - vp(i))**2 + (wg(i) - wp(i))**2)
        pint = pg(i) - MIN(delpstar, 0.01d0 * pg(i))

        uint = (phig(i) * rog(i) * ug(i) + (1.0 - phig(i)) * rop(i) * up(i)) / (phig(i) * rog(i) + (1.0 - phig(i)) * rop(i))
        vint = (phig(i) * rog(i) * vg(i) + (1.0 - phig(i)) * rop(i) * vp(i)) / (phig(i) * rog(i) + (1.0 - phig(i)) * rop(i))
        wint = (phig(i) * rog(i) * wg(i) + (1.0 - phig(i)) * rop(i) * wp(i)) / (phig(i) * rog(i) + (1.0 - phig(i)) * rop(i))

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

        phigi = phig(i)
        phigai = phig(nc1(i)) - phigi
        phigbi = phig(nc2(i)) - phigi
        phigci = phig(nc3(i)) - phigi
        phigdi = phig(nc4(i)) - phigi

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

        c1 = (phigai) * (xai) / wt1 + (phigbi) * (xbi) / wt2 + (phigci) * (xci) / wt3 + (phigdi) * (xdi) / wt4
        c2 = (phigai) * (yai) / wt1 + (phigbi) * (ybi) / wt2 + (phigci) * (yci) / wt3 + (phigdi) * (ydi) / wt4
        c3 = (phigai) * (zai) / wt1 + (phigbi) * (zbi) / wt2 + (phigci) * (zci) / wt3 + (phigdi) * (zdi) / wt4

        dnr = det(a1, b1, b2, b1, a2, b3, b2, b3, a3)
        dphigx = det(c1, c2, c3, b1, a2, b3, b2, b3, a3) / dnr
        dphigy = det(a1, b1, b2, c1, c2, c3, b2, b3, a3) / dnr
        dphigz = det(a1, b1, b2, b1, a2, b3, c1, c2, c3) / dnr
        
#ifdef KE_TURB
           amug =  as*dsqrt(tg(i))/(1.0d0+ ts/tg(i))
	   amutg=0.0845d0*rog(i)*tkg(i)*tkg(i)/teg(i)
		
	   sss1=2.0d0*(delugx(i)**2+delvgy(i)**2+delwgz(i)**2)
	   sss2=(delugy(i)+delvgx(i))**2+(delugz(i)+delwgx(i))**2+(delvgz(i)+delwgy(i))**2
	   sss3=(2.0d0/3.0d0)*(delugx(i)+delvgy(i)+delwgz(i))**2
	   sssg=sss1+sss2-sss3
	   pkgrng=amutg*sssg
		
        if(pkgrng.gt.20.0d0*cn(i,nkg))pkgrng=20.0d0*cn(i,nkg)
        if(pkgrng.lt.0.05d0*cn(i,nkg))pkgrng=0.05d0*cn(i,nkg)
		
	   etagrng=dsqrt(pkgrng/amutg)*tkg(i)/teg(i)
	   c1estarg=1.42d0-(etagrng*(1.0d0-etagrng/4.377d0))/(1.0d0+0.012d0*etagrng**3)
       	   amutgbyr=amutg/dmax1(0.01d0*amug,amutg)
	   stkg=phig(i)*(pkgrng-rog(i)*teg(i))*amutgbyr
           steg=phig(i)*((c1estarg*pkgrng-1.68d0*rog(i)*teg(i))*teg(i)/tkg(i))*amutgbyr
#endif

        rhs(i,2) = rhs(i,2) - (fdgx + pint * dphigx + phig(i) * rog(i) * gravity) * voli - visug(i)
        rhs(i,3) = rhs(i,3) - (fdgy + pint * dphigy) * voli - visvg(i)
        rhs(i,4) = rhs(i,4) - (fdgz + pint * dphigz) * voli - viswg(i)
        rhs(i,5) = rhs(i,5) - (pint * (uint * dphigx + vint * dphigy + wint * dphigz) + &
           uint * fdgx + vint * fdgy + wint * fdgz + qsource) * voli - vistg(i) 

        rhs(i,7) = rhs(i,7) + (fdgx + pint * dphigx - phip(i) * rop(i) * gravity) *voli - visup(i)        

        rhs(i,8) = rhs(i,8) + (fdgy + pint * dphigy) * voli - visvp(i)

        rhs(i,9) = rhs(i,9) + (fdgz + pint * dphigz) * voli - viswp(i)

        rhs(i,10) = rhs(i,10) + (pint * (uint * dphigx + vint * dphigy + wint * dphigz) + &
            uint * fdgx + vint * fdgy + wint * fdgz + qsource) * voli - vistp(i)
            
#ifdef KE_TURB
        rhs(i,nkg) = rhs(i,nkg) - (stkg * voli) - vistkg(i)
	rhs(i,neg) = rhs(i,neg) - (steg * voli) - visteg(i)
#endif
	
#ifdef IREACT
          yfghat=yfg(i)
	  yoghat=yog(i)/(1.0d0+sratio)
	  ypghat=ypg(i)/sratio
	  yghatmin=dmin1(yfghat,yoghat)
	  grankappag=(yghatmin+ypghat)**2/((yfghat+ypghat)*(yoghat+ypghat))
	  grangamag=21.0d0*(amug*teg(i))/(cn(i,1)*tkg(i)*tkg(i))**0.25d0
	  omegadotfg=-11.1d0*grankappag*teg(i)/tkg(i)*dmin1(yghatmin,ypghat,(yghatmin+ypghat)*grangamag)
	  omegadotog=sratio*omegadotfg
	  omegadotpg=-(1.0d0+sratio)+omegadotfg !Ask annie
	  
	  rhs(i,nyfg)=rhs(i,nyfg)-omegadotfg*voli
	  rhs(i,nyog)=rhs(i,nyog)-omegadotog*voli
	  rhs(i,nypg)=rhs(i,nypg)-omegadotpg*voli
#endif	   
            
END DO
!$acc end parallel loop 

end subroutine source

