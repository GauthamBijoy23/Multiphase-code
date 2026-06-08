SUBROUTINE derivative
    USE GlobalVariables
    USE MatrixOps
    IMPLICIT NONE
    external geometry
    external primitive
    external smooth
    
    INTEGER :: i
    REAL(DP) :: zi,yi,xi,wt4,wt3,wt2,wt1
    REAL(DP) :: wpi,wpdi,wpci,wpbi,wpai,wgi,wgdi,wgci,wgbi,wgai,vpi,vpdi,vpci,vpbi,vpai
    REAL(DP) :: vgi,vgdi,vgci,vgbi,vgai,upi,updi,upci,upbi,upai,ugi,ugdi,ugci,ugbi,ugai
    REAL(DP) :: tpi,tpdi,tpci,tpbi,tpai,tgi,tgdi,tgci,tgbi,tgai,ropi,ropdi,ropci,ropai
    REAL(DP) :: ropbi,rogi,rogdi,rogci,rogbi,rogai,ppi,ppdi,ppci,ppbi,ppai,phigi,phigdi
    REAL(DP) :: phigci,phigbi,phigai,pgi,pgdi,pgci,pgbi,pgai,a1,a2,a3,b1,b2,b3
    REAL(DP) :: cpg1,cpg2,cpg3,cphig1,cphig2,cphig3,cpp1,cpp2,cpp3,crog1,crog2
    REAL(DP) :: crog3,crop1,crop2,crop3,ctg1,ctg2,ctg3,ctp1,ctp2,ctp3
    REAL(DP) :: cug1,cug2,cug3,cup1,cup2,cup3,cvg1,cvg2,cvg3,cvp1,cvp2,cvp3
    REAL(DP) :: cwg1,cwg2,cwg3,cwp1,cwp2,cwp3,dnr,xai,xbi,xci,xdi,yai,ybi,yci,ydi
    REAL(DP) :: zai,zbi,zci,zdi
    REAL(DP) :: tkgi,tkgai,tkgbi,tkgci,tkgdi,tegi,tegai,tegbi,tegci,tegdi
    REAL(DP) :: cTkg1,cTkg2,cTkg3,cTkg4,cTeg1,cTeg2,cTeg3,cTeg4
    
!$acc parallel loop present(neles,delphigx(:),delropy(:),delrogz(:), &
!$acc delrogx(:),delppy(:),delphigz(:),delpgx(:),deltpy(:),deltgz(:), &
!$acc delwpx(:),delvpy(:),delvgz(:),delupx(:),deltgy(:),delropz(:), &
!$acc delwgx(:),delvgy(:),delupz(:),delugx(:),phig(:),nc4(:),rog(:), &
!$acc tg(:),wg(:),vg(:),ug(:),pg(:),pp(:),rop(:),vp(:),up(:),tp(:), &
!$acc zcel(:),ycel(:),xcel(:),nc2(:),wp(:),delwpz(:),nc1(:),nc3(:), &
!$acc deltpz(:),delugy(:),delvgx(:),delvpz(:),delwgy(:),deltgx(:), &
!$acc delugz(:),delupy(:),delvpx(:),delwgz(:),delwpy(:),deltpx(:), &
!$acc delpgy(:),delppx(:),delppz(:),delrogy(:),delropx(:),delpgz(:), & 
!$acc delphigy(:),delTkgx(:),delTkgy(:),delTkgz(:),delTegx(:),delTegy(:),delTegz(:), &
!$acc tkg(:),teg(:)) &
!$acc private(a1, a2, a3, b1, b2, b3, cpg1, cpg2, cpg3, cpp1, cpp2, cpp3, crop1, crop2, crop3, &
!$acc crog1, crog2, crog3, cphig1, cphig2, cphig3, ctp1, ctp2, ctp3, ctg1, ctg2, ctg3, &
!$acc cup1, cup2, cup3, cvg1, cvg2, cvg3, cvp1, cvp2, cvp3, cug1, cug2, cug3, &
!$acc cwg1, cwg2, cwg3, cwp1, cwp2, cwp3, dnr, pgi, pgdi, pgai, pgbi, pgci, phigai, phigbi, phigci, &
!$acc phigdi, phigi, ppi, ppdi, ppai, ppbi, ppci, rogi, rogdi, rogai, rogbi, rogci, ropi, ropdi, ropai, ropbi, ropci, &
!$acc tgi, tgai, tgbi, tgci, tgdi, tpai, tpbi, tpci, tpdi, tpi, ugai, ugbi, ugci, ugdi, ugi, upai, upbi, upci, updi, upi, &
!$acc vgi, vgai, vgbi, vgci, vgdi, vpai, vpbi, vpci, vpdi, vpi, wgi, wgai, wgbi, wgci, wgdi, wpai, wpbi, wpci, wpdi, wpi, &
!$acc wt1, wt2, wt3, wt4, xai, xbi, xci, xdi, yai, ybi, yci, ydi, zai, zbi, zci, zdi, xi, yi, zi , & 
!$acc tkgi,tkgai,tkgbi,tkgci,tkgdi,tegi,tegai,tegbi,tegci,tegdi,cTkg1,cTkg2,cTkg3,cTkg4,cTeg1,cTeg2,cTeg3,cTeg4)

  do i = 1, neles

IF (iblank(i) /= 1) CYCLE  !------ONLY PROCEED IF IBLANK = 1 ----------------(mod3)

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

    !if (myid.eq.0.and.k.eq.2) print*,"xi", i, myid, xi, xai, xbi, xci, xdi

    ugi = ug(i)
    ugai = ug(nc1(i)) - ugi	
    ugbi = ug(nc2(i)) - ugi			
    ugci = ug(nc3(i)) - ugi
    ugdi = ug(nc4(i)) - ugi

    vgi = vg(i)
    vgai = vg(nc1(i)) - vgi	
    vgbi = vg(nc2(i)) - vgi			
    vgci = vg(nc3(i)) - vgi
    vgdi = vg(nc4(i)) - vgi

    wgi = wg(i)
    wgai = wg(nc1(i)) - wgi	
    wgbi = wg(nc2(i)) - wgi			
    wgci = wg(nc3(i)) - wgi
    wgdi = wg(nc4(i)) - wgi

    Tgi = Tg(i)
    Tgai = Tg(nc1(i)) - Tgi
    Tgbi = Tg(nc2(i)) - Tgi
    Tgci = Tg(nc3(i)) - Tgi
    Tgdi = Tg(nc4(i)) - Tgi

    rogi = rog(i)
    rogai = rog(nc1(i)) - rogi
    rogbi = rog(nc2(i)) - rogi
    rogci = rog(nc3(i)) - rogi
    rogdi = rog(nc4(i)) - rogi
!   if (myid.eq.0) print*,"ugai", i, ugai, vgai, wgai, tgai

    Phigi = Phig(i)
    Phigai = Phig(nc1(i)) - Phigi
    Phigbi = Phig(nc2(i)) - Phigi
    Phigci = Phig(nc3(i)) - Phigi
    Phigdi = Phig(nc4(i)) - Phigi

    upi = up(i)
    upai = up(nc1(i)) - upi
    upbi = up(nc2(i)) - upi
    upci = up(nc3(i)) - upi
    updi = up(nc4(i)) - upi

    vpi = vp(i)
    vpai = vp(nc1(i)) - vpi
    vpbi = vp(nc2(i)) - vpi
    vpci = vp(nc3(i)) - vpi
    vpdi = vp(nc4(i)) - vpi

    wpi = wp(i)
    wpai = wp(nc1(i)) - wpi
    wpbi = wp(nc2(i)) - wpi
    wpci = wp(nc3(i)) - wpi
    wpdi = wp(nc4(i)) - wpi

    Tpi = Tp(i)
    Tpai = Tp(nc1(i)) - Tpi
    Tpbi = Tp(nc2(i)) - Tpi
    Tpci = Tp(nc3(i)) - Tpi
    Tpdi = Tp(nc4(i)) - Tpi

    ropi = rop(i)
    ropai = rop(nc1(i)) - ropi
    ropbi = rop(nc2(i)) - ropi
    ropci = rop(nc3(i)) - ropi
    ropdi = rop(nc4(i)) - ropi

    Pgi = Pg(i)
    Pgai = Pg(nc1(i)) - Pgi
    Pgbi = Pg(nc2(i)) - Pgi
    Pgci = Pg(nc3(i)) - Pgi
    Pgdi = Pg(nc4(i)) - Pgi

    Ppi = Pp(i)
    Ppai = Pp(nc1(i)) - Ppi
    Ppbi = Pp(nc2(i)) - Ppi
    Ppci = Pp(nc3(i)) - Ppi
    Ppdi = Pp(nc4(i)) - Ppi

#ifdef KE_TURB
        tkgi  = tkg(i)
        tkgai = tkg(nc1(i)) - tkgi
        tkgbi = tkg(nc2(i)) - tkgi
        tkgci = tkg(nc3(i)) - tkgi
        tkgdi = tkg(nc4(i)) - tkgi
		
	tegi  = teg(i)
        tegai = teg(nc1(i)) - tegi
        tegbi = teg(nc2(i)) - tegi
        tegci = teg(nc3(i)) - tegi
        tegdi = teg(nc4(i)) - tegi
#endif	

! Distance-based weighted derivative computation
    wt1 = (xai)**2 + (yai)**2 + (zai)**2
    wt2 = (xbi)**2 + (ybi)**2 + (zbi)**2
    wt3 = (xci)**2 + (yci)**2 + (zci)**2
    wt4 = (xdi)**2 + (ydi)**2 + (zdi)**2

    a1 = (xai)**2 / wt1 + (xbi)**2 / wt2 + (xci)**2 / wt3 + (xdi)**2 / wt4
    a2 = (yai)**2 / wt1 + (ybi)**2 / wt2 + (yci)**2 / wt3 + (ydi)**2 / wt4
    a3 = (zai)**2 / wt1 + (zbi)**2 / wt2 + (zci)**2 / wt3 + (zdi)**2 / wt4

    !if (myid.eq.1) print*, "xai",i, xdi, ydi, zdi

    b1 = (xai) * (yai) / wt1 + (xbi) * (ybi) / wt2 + (xci) * (yci) / wt3 + (xdi) * (ydi) / wt4
    b2 = (xai) * (zai) / wt1 + (xbi) * (zbi) / wt2 + (xci) * (zci) / wt3 + (xdi) * (zdi) / wt4
    b3 = (zai) * (yai) / wt1 + (zbi) * (ybi) / wt2 + (zci) * (yci) / wt3 + (zdi) * (ydi) / wt4

! dux duy computation
    cug1 = (ugai) * (xai) / wt1 + (ugbi) * (xbi) / wt2 + (ugci) * (xci) / wt3 + (ugdi) * (xdi) / wt4
    cug2 = (ugai) * (yai) / wt1 + (ugbi) * (ybi) / wt2 + (ugci) * (yci) / wt3 + (ugdi) * (ydi) / wt4
    cug3 = (ugai) * (zai) / wt1 + (ugbi) * (zbi) / wt2 + (ugci) * (zci) / wt3 + (ugdi) * (zdi) / wt4

! dvx dvy computation
     cvg1 = (vgai) * (xai) / wt1 + (vgbi) * (xbi) / wt2 + (vgci) * (xci) / wt3 + (vgdi) * (xdi) / wt4
     cvg2 = (vgai) * (yai) / wt1 + (vgbi) * (ybi) / wt2 + (vgci) * (yci) / wt3 + (vgdi) * (ydi) / wt4
     cvg3 = (vgai) * (zai) / wt1 + (vgbi) * (zbi) / wt2 + (vgci) * (zci) / wt3 + (vgdi) * (zdi) / wt4

! dwx dwy computation
     cwg1 = (wgai) * (xai) / wt1 + (wgbi) * (xbi) / wt2 + (wgci) * (xci) / wt3 + (wgdi) * (xdi) / wt4
     cwg2 = (wgai) * (yai) / wt1 + (wgbi) * (ybi) / wt2 + (wgci) * (yci) / wt3 + (wgdi) * (ydi) / wt4
     cwg3 = (wgai) * (zai) / wt1 + (wgbi) * (zbi) / wt2 + (wgci) * (zci) / wt3 + (wgdi) * (zdi) / wt4

! dTx dTy computation
     cTg1 = (Tgai) * (xai) / wt1 + (Tgbi) * (xbi) / wt2 + (Tgci) * (xci) / wt3 + (Tgdi) * (xdi) / wt4
     cTg2 = (Tgai) * (yai) / wt1 + (Tgbi) * (ybi) / wt2 + (Tgci) * (yci) / wt3 + (Tgdi) * (ydi) / wt4
     cTg3 = (Tgai) * (zai) / wt1 + (Tgbi) * (zbi) / wt2 + (Tgci) * (zci) / wt3 + (Tgdi) * (zdi) / wt4

! drox droy computation
     crog1 = (rogai) * (xai) / wt1 + (rogbi) * (xbi) / wt2 + (rogci) * (xci) / wt3 + (rogdi) * (xdi) / wt4
     crog2 = (rogai) * (yai) / wt1 + (rogbi) * (ybi) / wt2 + (rogci) * (yci) / wt3 + (rogdi) * (ydi) / wt4
     crog3 = (rogai) * (zai) / wt1 + (rogbi) * (zbi) / wt2 + (rogci) * (zci) / wt3 + (rogdi) * (zdi) / wt4

! dpx dpy computation
     cphig1 = (Phigai) * (xai) / wt1 + (Phigbi) * (xbi) / wt2 + (Phigci) * (xci) / wt3 + (Phigdi) * (xdi) / wt4
     cphig2 = (Phigai) * (yai) / wt1 + (Phigbi) * (ybi) / wt2 + (Phigci) * (yci) / wt3 + (Phigdi) * (ydi) / wt4
     cphig3 = (Phigai) * (zai) / wt1 + (Phigbi) * (zbi) / wt2 + (Phigci) * (zci) / wt3 + (Phigdi) * (zdi) / wt4

! dux duy computation
    cup1 = (upai) * (xai) / wt1 + (upbi) * (xbi) / wt2 + (upci) * (xci) / wt3 + (updi) * (xdi) / wt4
    cup2 = (upai) * (yai) / wt1 + (upbi) * (ybi) / wt2 + (upci) * (yci) / wt3 + (updi) * (ydi) / wt4
    cup3 = (upai) * (zai) / wt1 + (upbi) * (zbi) / wt2 + (upci) * (zci) / wt3 + (updi) * (zdi) / wt4

! dvx dvy computation
   cvp1 = (vpai) * (xai) / wt1 + (vpbi) * (xbi) / wt2 + (vpci) * (xci) / wt3 + (vpdi) * (xdi) / wt4
   cvp2 = (vpai) * (yai) / wt1 + (vpbi) * (ybi) / wt2 + (vpci) * (yci) / wt3 + (vpdi) * (ydi) / wt4
   cvp3 = (vpai) * (zai) / wt1 + (vpbi) * (zbi) / wt2 + (vpci) * (zci) / wt3 + (vpdi) * (zdi) / wt4

! dwx dwy computation
   cwp1 = (wpai) * (xai) / wt1 + (wpbi) * (xbi) / wt2 + (wpci) * (xci) / wt3 + (wpdi) * (xdi) / wt4
   cwp2 = (wpai) * (yai) / wt1 + (wpbi) * (ybi) / wt2 + (wpci) * (yci) / wt3 + (wpdi) * (ydi) / wt4
   cwp3 = (wpai) * (zai) / wt1 + (wpbi) * (zbi) / wt2 + (wpci) * (zci) / wt3 + (wpdi) * (zdi) / wt4

! dTx dTy computation
   cTp1 = (Tpai) * (xai) / wt1 + (Tpbi) * (xbi) / wt2 + (Tpci) * (xci) / wt3 + (Tpdi) * (xdi) / wt4
   cTp2 = (Tpai) * (yai) / wt1 + (Tpbi) * (ybi) / wt2 + (Tpci) * (yci) / wt3 + (Tpdi) * (ydi)
   cTp3 = (Tpai) * (zai) / wt1 + (Tpbi) * (zbi) / wt2 + (Tpci) * (zci) / wt3 + (Tpdi) * (zdi) / wt4

! drox droy computation 
    crop1 = (ropai) * (xai) / wt1 + (ropbi) * (xbi) / wt2 + (ropci) * (xci) / wt3 + (ropdi) * (xdi) / wt4
    crop2 = (ropai) * (yai) / wt1 + (ropbi) * (ybi) / wt2 + (ropci) * (yci) / wt3 + (ropdi) * (ydi) / wt4
    crop3 = (ropai) * (zai) / wt1 + (ropbi) * (zbi) / wt2 + (ropci) * (zci) / wt3 + (ropdi) * (zdi) / wt4
	
	! cpg1, cpg2, cpg3 computation
   cpg1 = (Pgai) * (xai) / wt1 + (Pgbi) * (xbi) / wt2 + (Pgci) * (xci) / wt3 + (Pgdi) * (xdi) / wt4
   cpg2 = (Pgai) * (yai) / wt1 + (Pgbi) * (ybi) / wt2 + (Pgci) * (yci) / wt3 + (Pgdi) * (ydi) / wt4
   cpg3 = (Pgai) * (zai) / wt1 + (Pgbi) * (zbi) / wt2 + (Pgci) * (zci) / wt3 + (Pgdi) * (zdi) / wt4

! dppx, dppy computation
   cpp1 = (Ppai) * (xai) / wt1 + (Ppbi) * (xbi) / wt2 + (Ppci) * (xci) / wt3 + (Ppdi) * (xdi) / wt4
   cpp2 = (Ppai) * (yai) / wt1 + (Ppbi) * (ybi) / wt2 + (Ppci) * (yci) / wt3 + (Ppdi) * (ydi) / wt4
   cpp3 = (Ppai) * (zai) / wt1 + (Ppbi) * (zbi) / wt2 + (Ppci) * (zci) / wt3 + (Ppdi) * (zdi) / wt4   
   
#ifdef KE_TURB
! dTkx dTky computation
     cTkg1 = (Tkgai) * (xai) / wt1 + (Tkgbi) * (xbi) / wt2 + (Tkgci) * (xci) / wt3 + (Tkgdi) * (xdi) / wt4
     cTkg2 = (Tkgai) * (yai) / wt1 + (Tkgbi) * (ybi) / wt2 + (Tkgci) * (yci) / wt3 + (Tkgdi) * (ydi) / wt4
     cTkg3 = (Tkgai) * (zai) / wt1 + (Tkgbi) * (zbi) / wt2 + (Tkgci) * (zci) / wt3 + (Tkgdi) * (zdi) / wt4
	 
! dTex dTey computation
     cTeg1 = (Tegai) * (xai) / wt1 + (Tegbi) * (xbi) / wt2 + (Tegci) * (xci) / wt3 + (Tegdi) * (xdi) / wt4
     cTeg2 = (Tegai) * (yai) / wt1 + (Tegbi) * (ybi) / wt2 + (Tegci) * (yci) / wt3 + (Tegdi) * (ydi) / wt4
     cTeg3 = (Tegai) * (zai) / wt1 + (Tegbi) * (zbi) / wt2 + (Tegci) * (zci) / wt3 + (Tegdi) * (zdi) / wt4	 
#endif      
         
    dnr=det(a1,b1,b2,b1,a2,b3,b2,b3,a3)
    !if(myid.eq.1) print*, "a1", i, a1, a2, a3
	 
! Computation for delugx, delugy, delugz, etc.
   delugx(i) = det(cug1, cug2, cug3, b1, a2, b3, b2, b3, a3) / dnr
   delugy(i) = det(a1, b1, b2, cug1, cug2, cug3, b2, b3, a3) / dnr
   delugz(i) = det(a1, b1, b2, b1, a2, b3, cug1, cug2, cug3) / dnr

   delvgx(i) = det(cvg1, cvg2, cvg3, b1, a2, b3, b2, b3, a3) / dnr
   delvgy(i) = det(a1, b1, b2, cvg1, cvg2, cvg3, b2, b3, a3) / dnr
   delvgz(i) = det(a1, b1, b2, b1, a2, b3, cvg1, cvg2, cvg3) / dnr

   delwgx(i) = det(cwg1, cwg2, cwg3, b1, a2, b3, b2, b3, a3) / dnr
   delwgy(i) = det(a1, b1, b2, cwg1, cwg2, cwg3, b2, b3, a3) / dnr
   delwgz(i) = det(a1, b1, b2, b1, a2, b3, cwg1, cwg2, cwg3) / dnr

   delTgx(i) = det(cTg1, cTg2, cTg3, b1, a2, b3, b2, b3, a3) / dnr
   delTgy(i) = det(a1, b1, b2, cTg1, cTg2, cTg3, b2, b3, a3) / dnr
   delTgz(i) = det(a1, b1, b2, b1, a2, b3, cTg1, cTg2, cTg3) / dnr

   delupx(i) = det(cup1, cup2, cup3, b1, a2, b3, b2, b3, a3) / dnr
   delupy(i) = det(a1, b1, b2, cup1, cup2, cup3, b2, b3, a3) / dnr
   delupz(i) = det(a1, b1, b2, b1, a2, b3, cup1, cup2, cup3) / dnr

   delvpx(i) = det(cvp1, cvp2, cvp3, b1, a2, b3, b2, b3, a3) / dnr
   delvpy(i) = det(a1, b1, b2, cvp1, cvp2, cvp3, b2, b3, a3) / dnr
   delvpz(i) = det(a1, b1, b2, b1, a2, b3, cvp1, cvp2, cvp3) / dnr

   delwpx(i) = det(cwp1, cwp2, cwp3, b1, a2, b3, b2, b3, a3) / dnr
   delwpy(i) = det(a1, b1, b2, cwp1, cwp2, cwp3, b2, b3, a3) / dnr
   delwpz(i) = det(a1, b1, b2, b1, a2, b3, cwp1, cwp2, cwp3) / dnr

   delTpx(i) = det(cTp1, cTp2, cTp3, b1, a2, b3, b2, b3, a3) / dnr
   delTpy(i) = det(a1, b1, b2, cTp1, cTp2, cTp3, b2, b3, a3) / dnr
   delTpz(i) = det(a1, b1, b2, b1, a2, b3, cTp1, cTp2, cTp3) / dnr

   delPgx(i) = det(cPg1, cPg2, cPg3, b1, a2, b3, b2, b3, a3) / dnr
   delPgy(i) = det(a1, b1, b2, cPg1, cPg2, cPg3, b2, b3, a3) / dnr
   delPgz(i) = det(a1, b1, b2, b1, a2, b3, cPg1, cPg2, cPg3) / dnr

   delPpx(i) = det(cPp1, cPp2, cPp3, b1, a2, b3, b2, b3, a3) / dnr
   delPpy(i) = det(a1, b1, b2, cPp1, cPp2, cPp3, b2, b3, a3) / dnr
   delPpz(i) = det(a1, b1, b2, b1, a2, b3, cPp1, cPp2, cPp3) / dnr

   delrogx(i) = det(crog1, crog2, crog3, b1, a2, b3, b2, b3, a3) / dnr
   delrogy(i) = det(a1, b1, b2, crog1, crog2, crog3, b2, b3, a3) / dnr
   delrogz(i) = det(a1, b1, b2, b1, a2, b3, crog1, crog2, crog3) / dnr
    
    delropx(i) = det(crop1, crop2, crop3, b1, a2, b3, b2, b3, a3) / dnr 
    delropy(i) = det(a1, b1, b2, crop1, crop2, crop3, b2, b3, a3) / dnr
    delropz(i) = det(a1, b1, b2, b1, a2, b3, crop1, crop2, crop3) / dnr

    delphigx(i) = det(cphig1, cphig2, cphig3, b1, a2, b3, b2, b3, a3) / dnr
    delphigy(i) = det(a1, b1, b2, cphig1, cphig2, cphig3, b2, b3, a3) / dnr
    delphigz(i) = det(a1, b1, b2, b1, a2, b3, cphig1, cphig2, cphig3) / dnr

   !if (myid==1) write(*,'(A,1X,I8,3(1X,F25.15))') "del-der", i, delvgx(i), delrogy(i), delTgz(i)
#ifdef KE_TURB
    delTkgx(i) = det(cTkg1, cTkg2, cTkg3, b1, a2, b3, b2, b3, a3) / dnr
    delTkgy(i) = det(a1, b1, b2, cTkg1, cTkg2, cTkg3, b2, b3, a3) / dnr
    delTkgz(i) = det(a1, b1, b2, b1, a2, b3, cTkg1, cTkg2, cTkg3) / dnr
	
    delTegx(i) = det(cTeg1, cTeg2, cTeg3, b1, a2, b3, b2, b3, a3) / dnr
    delTegy(i) = det(a1, b1, b2, cTeg1, cTeg2, cTeg3, b2, b3, a3) / dnr
    delTegz(i) = det(a1, b1, b2, b1, a2, b3, cTeg1, cTeg2, cTeg3) / dnr
#endif	   

enddo
!$acc end parallel loop	

CALL mpiexdel

end subroutine derivative

