SUBROUTINE hllcbatten
    USE GlobalVariables
    IMPLICIT NONE
    external geometry
    external primitive
    external smooth

    integer :: i,neigh,ng,np,ns,nt,n1,n2,n3,n4,j,iface,k

    real(dp) :: eg, ep, bbb, ccc,psig, functiong
    real(dp) :: zfn,zf0,zf,yfn,yf0,yf,xfn,xf0,xf
    REAL(DP) :: xcenter,x4,x3,x2,x1,wpr,wpl,wgr,wgl,vpr,vpl,vgr,vgl,upr,upl,ugr,ugl
    REAL(DP) :: tpr,tpl,tgr,tgl,sstar,sr,sl,scz,scy,scx,ror,roqr,roql,ropr
    REAL(DP) :: rol,ropl,rogl,qr,qpr,qpl,ql,rogr,qgr,qgl,prr,prl,ppstar,ppr,ppl
    REAL(DP) :: pgl,pgr,pgstar,phigl,phigr,phimin,phipl,phipr,hpr,hpl,hgr,hgl
    REAL(DP) :: fstara8,fstara7,fstara6,fstara5,fstara4,fstara3,fstara2,fstara1
    REAL(DP) :: fr1,fr2,fr3,fr4,fr5,fr6,fr7,fr8,fr9,FSTARA10,FSTARA9,fr10
    REAL(DP) :: fl1,fl2,fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10,epr,epl,enz,eny,enx
    REAL(DP) :: egr,egl,agr,apl,apr,area,cnl1,cnl10,cnl2,cnl3,cnl4
    REAL(DP) :: cnl7,cnl8,cnl9,agl,CNSTARR9,CNSTARR8,CNSTARR7
    REAL(DP) :: CNSTARR6,CNSTARR5,CNSTARR4,CNSTARR3,CNSTARR2,CNSTARR1,dnr
    REAL(DP) :: CNL6,CNR1,CNR10,CNR3,CNR4,CNR5,CNR6
    REAL(DP) :: CNR7,cnr2,cnr8,cnl5
    REAL(DP) :: sgl,spl,sgr,spr,sgstarnr,sgstardr,sgstar,spstarnr,spstardr,spstar,dnrg,dnrp 
    REAL(DP) :: CNSTARL1,CNSTARL10,CNSTARL2,CNSTARL4,CNSTARL5,CNSTARL6,CNSTARL7
    REAL(DP) :: CNSTARL8,CNSTARL9,CNSTARL3,CNSTARR10,cnr9
    REAL(DP) :: epsilon=1.0E-8
    
  call primitiveplusbc

  if (myorder == 2) then
   call derivative
   call venkatakrishnan
  end if

!$acc parallel loop present(pp(:),phig(:),nod(:,:),tg(:),wg(:),vg(:),ug(:),phi(:),sc4z(:), &
!$acc vp(:),up(:),tp(:),ycel(:),xcel(:),wp(:),sc4y(:),sc4x(:),sc3z(:),nc3(:),& 
!$acc nc2(:),sc2z(:),sc3x(:),sc3y(:),nc1(:),sc1z(:),sc2x(:),sc2y(:),nc4(:),delwpz(:), &
!$acc rhs(:,:),sc1x(:),sc1y(:),delupy(:),delupx(:),delugz(:),delvpy(:),delvpx(:), &
!$acc delvgz(:),delwpy(:),delwpx(:),delwgz(:),deltpy(:),deltpx(:),deltgz(:),delppy(:), & 
!$acc delppx(:),delphigz(:),delpgx(:),delpgy(:),delppz(:),deltgx(:),deltgy(:),delvpz(:), & 
!$acc delwgx(:),delwgy(:),delupz(:),delvgx(:),delvgy(:),deltpz(:),delugx(:),delugy(:), &
!$acc delpgz(:),delphigx(:),delphigy(:),zcel(:),z(:),y(:),x(:),pg(:),rrg,pinf,epslnmax, &
!$acc gamag,gamap,epslnmin,myorder,neles,cpp,ntype(:),nghosts) &  
!$acc  private(scx,scy,wpr,xfn,yfn,apr,egr,enx,eny,tpr,upr,vpr,spstarnr,pgstar,ugr,vgr, &
!$acc wgr,tgr,neigh,ugl,vgl,wgl,tgl,phigl,upl,vpl,wpl,tpl,pgl,ppl,hpr,phigr, &
!$acc xf,yf,zf,xf0,yf0,zfn,zf0,phimin,phipl,epr,ror,qpr,phipr,cnstarr9,rol,rogl, & 
!$acc ropl,dnrp,egl,hgl,qgl,agl,enz,hgr,qgr,agr,epl,hpl,qpl,apl,scz,sgstarnr, &
!$acc sgl,spl,rogr,ropr,roqr,roql,sgstardr,sgstar,sgr,spstardr,spstar,spr, &
!$acc pgr,ppr,area,cnl10,cnl2,cnl3,cnl4,cnl5,cnl6,cnl7,cnl8,cnl1,cnl9,cnr10, & 
!$acc cnr2,cnr3,cnr4,cnr5,cnr6,cnr7,cnr8,cnr1,cnr9,cnstarl10,cnstarl2,cnstarl3, &
!$acc cnstarl4,dnrg,cnstarl5,cnstarl6,cnstarl7,cnstarl8,cnstarl1,cnstarl9,cnstarr10,& 
!$acc cnstarr2,cnstarr3,cnstarr4,cnstarr5,cnstarr6,cnstarr7,cnstarr8,cnstarr1, &
!$acc fl10,fl2,fl3,fl4,fl5,fl6,fl7,fl8,fl1,fl9,fr10,fr2,fr3,fr4,fr5,fr6,fr7,fr8, &
!$acc fr1,fr9,fstara10,fstara2,fstara3,fstara4,fstara5,fstara6,fstara7,fstara8, &
!$acc fstara1,functiong,fstara9,psig,ppstar,n3,n2,n1,iface,nt,sgstarnr,sgstardr,sgstar, &
!$acc spstarnr,spstardr,spstar,rogr,ropr,roqr,roql,dnrg,dnrp)

  do i = 1, neles

IF (iblank(i) /= 1) CYCLE  !------ONLY PROCEED IF IBLANK = 1 ----------------(mod3)
    ! ---------------------------------------------------------------
   rhs(i, 1) = 0.0d0
   rhs(i, 2) = 0.0d0
   rhs(i, 3) = 0.0d0
   rhs(i, 4) = 0.0d0
   rhs(i, 5) = 0.0d0
   rhs(i, 6) = 0.0d0
   rhs(i, 7) = 0.0d0
   rhs(i, 8) = 0.0d0
   rhs(i, 9) = 0.0d0
   rhs(i, 10) = 0.0d0
   ! ---------------------------------------------------------------   
  
!$acc loop seq     
   do iface = 1, 4

    if (iface == 1) then
     n1 = nod(i, 2)
     n2 = nod(i, 4)
     n3 = nod(i, 3)
     neigh = NC1(i)
     SCX = SC1X(i)
     SCY = SC1Y(i)
     SCZ = SC1Z(i)

    elseif (iface == 2) then
     n1 = nod(i, 3)
     n2 = nod(i, 4)
     n3 = nod(i, 1)
     neigh = NC2(i)
     SCX = SC2X(i)
     SCY = SC2Y(i)
     SCZ = SC2Z(i)
       
    elseif (iface == 3) then
     n1 = nod(i, 4)
     n2 = nod(i, 2)
     n3 = nod(i, 1)
     neigh = NC3(i)
     SCX = SC3X(i)
     SCY = SC3Y(i)
     SCZ = SC3Z(i)

    elseif (iface == 4) then
     n1 = nod(i, 2)
     n2 = nod(i, 3)
     n3 = nod(i, 1)
     neigh = NC4(i)
     SCX = SC4X(i)
     SCY = SC4Y(i)
     SCZ = SC4Z(i)
    
    end if
 
    XF = (X(n1) + X(n2) + X(n3)) * 0.3333D0
    YF = (Y(n1) + Y(n2) + Y(n3)) * 0.3333D0
    ZF = (Z(n1) + Z(n2) + Z(n3)) * 0.3333D0

    AREA = DSQRT(SCX * SCX + SCY * SCY + SCZ * SCZ)

    ENX = SCX / AREA
    ENY = SCY / AREA
    ENZ = SCZ / AREA
              
    ! -------------------------------------------------------------
    nt = 0
    if ((neigh > neles).and.(neigh <= (neles+nghosts))) nt = ntype(neigh - neles) !for boundary ghosts 
    ! -------------------------------------------------------------
    !if (myid==1) print*, i, iface, nt 

    ! Store neighbour as right cells; before if loop; since uR, vR, HR, pR required in rhs
    ! ---------------------------------------------------------------------------------------
    rogL = rog(i)
    ugL = ug(i)
    vgL = vg(i)
    wgL = wg(i)
    tgL = tg(i)
    phigL = phig(i)
    ropL = rop(i)
    upL = up(i)
    vpL = vp(i)
    wpL = wp(i)
    tpL = tp(i)
    pgL = pg(i)
    ppL = pp(i)

    rogR = rog(neigh)
    ugR = ug(neigh)
    vgR = vg(neigh)
    wgR = wg(neigh)
    tgR = tg(neigh)
    phigR = phig(neigh)
    ropR = rop(neigh)
    upR = up(neigh)
    vpR = vp(neigh)
    wpR = wp(neigh)
    tpR = tp(neigh)
    pgR = pg(neigh)
    ppR = pp(neigh)
   !if (myid==2.and.i.eq.11) print*, i, rogl, ugl, vgl, wgl,tgl,ppl
!   if (myid == 1) write(*,'(2(1X,I8),4(1X,F25.15))') i, iface, ugR, upR, vpR, ppR

!-------------upr is correct here----------------------------------
    if (myorder == 2) then
    
          if ((neigh <= neles).or.(neigh > neles+nghosts)) then
    
      phimin = min(phi(i), phi(neigh))
      xf0 = xf - xcel(i)
      yf0 = yf - ycel(i)
      zf0 = zf - zcel(i)
      xfn = xf - xcel(neigh)
      yfn = yf - ycel(neigh)
      zfn = zf - zcel(neigh)

!      if(i.eq.6.and.myid.eq.1) print*, "phi", phimin
    ! -----------------------------------------------------------------------------------------------
    
      Phigl = Phigl + phimin * (xf0 * delPhigx(i) + yf0 * delPhigy(i) + zf0 * delPhigz(i))
      ugl   = ugl + phimin * (xf0 * delugx(i) + yf0 * delugy(i) + zf0 * delugz(i))
      vgl   = vgl + phimin * (xf0 * delvgx(i) + yf0 * delvgy(i) + zf0 * delvgz(i))
      wgl   = wgl + phimin * (xf0 * delwgx(i) + yf0 * delwgy(i) + zf0 * delwgz(i))
      Tgl   = Tgl + phimin * (xf0 * delTgx(i) + yf0 * delTgy(i) + zf0 * delTgz(i))
      Pgl   = Pgl + phimin * (xf0 * delPgx(i) + yf0 * delPgy(i) + zf0 * delPgz(i))
  
      !if (myid == 1) print*,  "delphi",i, delPhigz(i)

      Phigr = Phigr + phimin * (xfn * delPhigx(neigh) + yfn * delPhigy(neigh) + zfn * delPhigz(neigh))
      ugr = ugr + phimin * (xfn * delugx(neigh) + yfn * delugy(neigh) + zfn * delugz(neigh))
      vgr = vgr + phimin * (xfn * delvgx(neigh) + yfn * delvgy(neigh) + zfn * delvgz(neigh))
      wgr = wgr + phimin * (xfn * delwgx(neigh) + yfn * delwgy(neigh) + zfn * delwgz(neigh))
      Tgr = Tgr + phimin * (xfn * delTgx(neigh) + yfn * delTgy(neigh) + zfn * delTgz(neigh))        
      Pgr = Pgr + phimin * (xfn * delPgx(neigh) + yfn * delPgy(neigh) + zfn * delPgz(neigh)) 
     !if (myid==2.and.i.eq.11) print*, i, phigl, ugl, vgl, wgl,tgl,pgl
!   if (myid==2.and.i.eq.11) print*,"del from bat",  i,iface, delPhigx(neigh), delPhigy(neigh), delPhigz(neigh)

     endif
     
      if ((neigh <= neles) .or. (neigh > (neles+nghosts))) then
      
          xf0=xf-xcel(i)
          yf0=yf-ycel(i)
          zf0=zf-zcel(i)
          xfn=xf-xcel(neigh)
          yfn=yf-ycel(neigh)
          zfn=zf-zcel(neigh)
      
      upl = upl + phimin * (xf0 * delupx(i) + yf0 * delupy(i) + zf0 * delupz(i))
      vpl = vpl + phimin * (xf0 * delvpx(i) + yf0 * delvpy(i) + zf0 * delvpz(i))
      wpl = wpl + phimin * (xf0 * delwpx(i) + yf0 * delwpy(i) + zf0 * delwpz(i))
      Tpl = Tpl + phimin * (xf0 * delTpx(i) + yf0 * delTpy(i) + zf0 * delTpz(i))
      Ppl = Ppl + phimin * (xf0 * delPpx(i) + yf0 * delPpy(i) + zf0 * delPpz(i))  
    
 !     if(myid==1.and.i.eq.6) print*, upr, "before", iface
      upr = upr + phimin * (xfn * delupx(neigh) + yfn * delupy(neigh) + zfn * delupz(neigh))
      vpr = vpr + phimin * (xfn * delvpx(neigh) + yfn * delvpy(neigh) + zfn * delvpz(neigh))
      wpr = wpr + phimin * (xfn * delwpx(neigh) + yfn * delwpy(neigh) + zfn * delwpz(neigh))
      Tpr = Tpr + phimin * (xfn * delTpx(neigh) + yfn * delTpy(neigh) + zfn * delTpz(neigh))
      Ppr = Ppr + phimin * (xfn * delPpx(neigh) + yfn * delPpy(neigh) + zfn * delPpz(neigh))

     end if
    end if
    ! if (myid==1) write(*,'(2(1X,I8),4(1X,F25.15))') i, iface, ugR, upR, vpR, ppR

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PhipL = 1.0d0 - PhigL
    PhipR = 1.0d0 - PhigR

    if (phigL < epslnmin) then
            phigL = epslnmin
            phipL = 1.0d0 - phigL
            ugL = upL
            vgL = vpL
            wgL = wpL
            tgL = tpL
            pgL = ppL
    end if

    if (phipL < epslnmin) then
            phipL = epslnmin
            phigL = 1.0d0 - phipL
            upL = ugL
            vpL = vgL
            wpL = wgL
            tpL = tgL
            ppL = pgL
    end if

    if (phigL >= epslnmin .and. phigL <= epslnmax) then
            psig = (phigL - epslnmin) / (epslnmax - epslnmin)
            functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
            ugL = functiong * ugL + (1.0d0 - functiong) * upL  
            vgL = functiong * vgL + (1.0d0 - functiong) * vpL 
            wgL = functiong * wgL + (1.0d0 - functiong) * wpL 
            tgL = functiong * tgL + (1.0d0 - functiong) * tpL
            pgL = functiong * pgL + (1.0d0 - functiong) * ppL
    end if

    if (phipL >= epslnmin .and. phipL <= epslnmax) then
            psig = (phipL - epslnmin) / (epslnmax - epslnmin)
            functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
            upL = functiong * upL + (1.0d0 - functiong) * ugL  
            vpL = functiong * vpL + (1.0d0 - functiong) * vgL 
            wpL = functiong * wpL + (1.0d0 - functiong) * wgL 
            tpL = functiong * tpL + (1.0d0 - functiong) * tgL
            ppL = functiong * ppL + (1.0d0 - functiong) * pgL
    end if
   
    
    rogL = pgL / (rrg * tgL)
    ropL = gamap * (ppL + pinf) / ((gamap - 1.0d0) * cpp * tpL)

    if (phigR < epslnmin) then
            phigR = epslnmin
            phipR = 1.0d0 - phigR
            ugR = upR
            vgR = vpR
            wgR = wpR
            tgR = tpR
            pgR = ppR
    end if

    if (phipR < epslnmin) then
            phipR = epslnmin
            phigR = 1.0d0 - phipR
            upR = ugR
            vpR = vgR
            wpR = wgR
            tpR = tgR
            ppR = pgR
    end if
!   if(myid ==1) write(*,'(2(1X,I8),4(1X,F25.15))') i, iface, ugR, upR, vpR, ppR
 

    if (phigR >= epslnmin .and. phigR <= epslnmax) then
            psig = (phigR - epslnmin) / (epslnmax - epslnmin)
            functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
            ugR = functiong * ugR + (1.0d0 - functiong) * upR  
            vgR = functiong * vgR + (1.0d0 - functiong) * vpR 
            wgR = functiong * wgR + (1.0d0 - functiong) * wpR 
            tgR = functiong * tgR + (1.0d0 - functiong) * tpR
            pgR = functiong * pgR + (1.0d0 - functiong) * ppR
    end if

    if (phipR >= epslnmin .and. phipR <= epslnmax) then
            psig = (phipR - epslnmin) / (epslnmax - epslnmin)
            functiong = -psig * psig * (2.0d0 * psig - 3.0d0)
            upR = functiong * upR + (1.0d0 - functiong) * ugR  
            vpR = functiong * vpR + (1.0d0 - functiong) * vgR 
            wpR = functiong * wpR + (1.0d0 - functiong) * wgR 
            tpR = functiong * tpR + (1.0d0 - functiong) * tgR
            ppR = functiong * ppR + (1.0d0 - functiong) * pgR
    end if

    rogR = pgR / (rrg * tgR)
    ropR = gamap * (ppR + pinf) / ((gamap - 1.0d0) * cpp * tpR)
    !if(myid ==1) write(*,'(2(1X,I8),4(1X,F25.15))') i, iface, ugR, upR, vpR, ppR
    !print*, functiong


! ------------------------------------------------------------------------------------------------        
     IF(NT.EQ.5)THEN
!!$acc atomic update
     RHS(I, 2) = RHS(I, 2) + PHIGL * PGL * ENX * AREA
!!$acc atomic
     RHS(I, 3) = RHS(I, 3) + PHIGL * PGL * ENY * AREA
!!$acc atomic
     RHS(I, 4) = RHS(I, 4) + PHIGL * PGL * ENZ * AREA
!!$acc atomic
     RHS(I, 7) = RHS(I, 7) + PHIPL * PPL * ENX * AREA
!!$acc atomic
     RHS(I, 8) = RHS(I, 8) + PHIPL * PPL * ENY * AREA
!!$acc atomic update
     RHS(I, 9) = RHS(I, 9) + PHIPL * PPL * ENZ * AREA
    ELSE
        EGL=RRG*TGL/(GAMAG-1.0d0)+0.5d0*(UGL**2+VGL**2+WGL**2)
        HGL=EGL+PGL/ROGL
        QGL=UGL*ENX+VGL*ENY+WGL*ENZ
        AGL=DSQRT(GAMAG*PGL/ROGL)

        EGR=RRG*TGR/(GAMAG-1.0d0)+0.5d0*(UGR**2+VGR**2+WGR**2)
        HGR=EGR+PGR/ROGR
        QGR=UGR*ENX+VGR*ENY+WGR*ENZ
        AGR=DSQRT(GAMAG*PGR/ROGR)

        EPL=(CPP*TPL)/GAMAP+PINF/ROPL+0.5D0*(UPL**2+VPL**2+WPL**2) 
        HPL=EPL+PPL/ROPL
        QPL=UPL*ENX+VPL*ENY+WPL*ENZ
        APL=DSQRT(GAMAP*(PPL+PINF)/ROPL)

        EPR=(CPP*TPR)/GAMAP+PINF/ROPR+0.5D0*(UPR**2+VPR**2+WPR**2) 
        HPR=EPR+PPR/ROPR
        QPR=UPR*ENX+VPR*ENY+WPR*ENZ
        APR=DSQRT(GAMAP*(PPR+PINF)/ROPR)
    
!	-------------------------------------------------------------------------------------------------			  
        SGL=DMIN1(QGL-AGL,QGR-AGR,epsilon)
        SPL=DMIN1(QPL-APL,QPR-APR,epsilon)
        SGR=DMAX1(QGL+AGL,QGR+AGR,epsilon)
        SPR=DMAX1(QPL+APL,QPR+APR,epsilon)

               
        ROL=PHIGL*ROGL+PHIPL*ROPL
        ROQL=PHIGL*ROGL*QGL+PHIPL*ROPL*QPL
        QL=ROQL/ROL
	PRL=PHIGL*PGL+PHIPL*PPL       
       
        ROR=PHIGR*ROGR+PHIPR*ROPR
        ROQR=PHIGR*ROGR*QGR+PHIPR*ROPR*QPR
        QR=ROQR/ROR
        PRR=PHIGR*PGR+PHIPR*PPR
!	---------------------------------------------------------------
        SGSTARNR=PGR-PGL+ROGL*QGL*(SGL-QGL)-ROGR*QGR*(SGR-QGR)
        SGSTARDR=ROGL*(SGL-QGL)-ROGR*(SGR-QGR)
        SGSTAR=SGSTARNR/(SGSTARDR+1.0E-8)
!	---------------------------------------------------------------
        SPSTARNR=PPR-PPL+ROPL*QPL*(SPL-QPL)-ROPR*QPR*(SPR-QPR)
        SPSTARDR=ROPL*(SPL-QPL)-ROPR*(SPR-QPR)
        SPSTAR=SPSTARNR/(SPSTARDR+1.0E-8)
!	-----------------------------------------------------------	
        PGSTAR=DMAX1(PGL+ROGL*(QGL-SGL)*(QGL-SGSTAR),1.0E-6*DMIN1(PGL,PGR))
        PPSTAR=DMAX1(PPL+ROPL*(QPL-SPL)*(QPL-SPSTAR),1.0E-6*DMIN1(PPL,PPR))
  !      if(myid==1.and.i==1.and.iface==2)print*,"1st",PPSTAR,PPL,ROPL,QPL,SPL,QPL,SPSTAR,PPR
        !	-----------------------------------------------------------     

        CNL1 =PHIGL*ROGL
        CNL2 =PHIGL*ROGL*UGL
        CNL3 =PHIGL*ROGL*VGL
        CNL4 =PHIGL*ROGL*WGL
        CNL5 =PHIGL*ROGL*EGL
        CNL6 =PHIPL*ROPL
        CNL7 =PHIPL*ROPL*UPL
        CNL8 =PHIPL*ROPL*VPL
        CNL9 =PHIPL*ROPL*WPL
        CNL10=PHIPL*ROPL*EPL
                    
        CNR1 =PHIGR*ROGR
        CNR2 =PHIGR*ROGR*UGR
        CNR3 =PHIGR*ROGR*VGR
        CNR4 =PHIGR*ROGR*WGR
        CNR5 =PHIGR*ROGR*EGR
        CNR6 =PHIPR*ROPR
        CNR7 =PHIPR*ROPR*UPR
        CNR8 =PHIPR*ROPR*VPR
        CNR9 =PHIPR*ROPR*WPR
        CNR10=PHIPR*ROPR*EPR
        !if (myid==2.and.i.eq.11) print*, i,iface, phigr, phigl
     	DNRG=(SGSTAR-SGL)+1.0E-6
        CNSTARL1= (QGL-SGL)*CNL1/DNRG
        CNSTARL2=((QGL-SGL)*CNL2+PHIGL*(PGL*ENX-PGSTAR*ENX))/DNRG
        CNSTARL3=((QGL-SGL)*CNL3+PHIGL*(PGL*ENY-PGSTAR*ENY))/DNRG
	CNSTARL4=((QGL-SGL)*CNL4+PHIGL*(PGL*ENZ-PGSTAR*ENZ))/DNRG
        CNSTARL5=((QGL-SGL)*CNL5+PHIGL*(PGL*QGL-PGSTAR*SGSTAR))/DNRG
        DNRP=(SPSTAR-SPL)+1.0E-6
        CNSTARL6= (QPL-SPL)*CNL6/DNRP
        CNSTARL7=((QPL-SPL)*CNL7+PHIPL*(PPL*ENX-PPSTAR*ENX))/DNRP
        CNSTARL8=((QPL-SPL)*CNL8+PHIPL*(PPL*ENY-PPSTAR*ENY))/DNRP
	CNSTARL9=((QPL-SPL)*CNL9+PHIPL*(PPL*ENZ-PPSTAR*ENZ))/DNRP
!	if(myid==1.and.i==1.and.iface==2)print*,"this",QPL,SPL,CNL9,PHIPL,PPL,ENZ,PPSTAR,ENZ,DNRP
        CNSTARL10=((QPL-SPL)*CNL10+PHIPL*(PPL*QPL-PPSTAR*SPSTAR))/DNRP

! ---------------------------------------------------------------

  
        DNRG=(SGSTAR-SGR)+1.0E-6
        CNSTARR1= (QGR-SGR)*CNR1/DNRG
        CNSTARR2=((QGR-SGR)*CNR2+PHIGR*(PGR*ENX-PGSTAR*ENX))/DNRG
        CNSTARR3=((QGR-SGR)*CNR3+PHIGR*(PGR*ENY-PGSTAR*ENY))/DNRG
	CNSTARR4=((QGR-SGR)*CNR4+PHIGR*(PGR*ENZ-PGSTAR*ENZ))/DNRG
        CNSTARR5=((QGR-SGR)*CNR5+PHIGR*(PGR*QGR-PGSTAR*SGSTAR))/DNRG
        DNRP=(SPSTAR-SPR)+1.0E-6
        CNSTARR6= (QPR-SPR)*CNR6/DNRP
        CNSTARR7=((QPR-SPR)*CNR7+PHIPR*(PPR*ENX-PPSTAR*ENX))/DNRP
        CNSTARR8=((QPR-SPR)*CNR8+PHIPR*(PPR*ENY-PPSTAR*ENY))/DNRP
	CNSTARR9=((QPR-SPR)*CNR9+PHIPR*(PPR*ENZ-PPSTAR*ENZ))/DNRP
        CNSTARR10=((QPR-SPR)*CNR10+PHIPR*(PPR*QPR-PPSTAR*SPSTAR))/DNRP
!       ---------------------------------------------------------------

        FL1 =PHIGL*ROGL*QGL
        FL2 =PHIGL*ROGL*QGL*UGL+PHIGL*PGL*ENX
        FL3 =PHIGL*ROGL*QGL*VGL+PHIGL*PGL*ENY
        FL4 =PHIGL*ROGL*QGL*WGL+PHIGL*PGL*ENZ
        FL5 =PHIGL*ROGL*QGL*HGL
        FL6 =PHIPL*ROPL*QPL
        FL7 =PHIPL*ROPL*QPL*UPL+PHIPL*PPL*ENX
        FL8 =PHIPL*ROPL*QPL*VPL+PHIPL*PPL*ENY
        FL9 =PHIPL*ROPL*QPL*WPL+PHIPL*PPL*ENZ
        FL10=PHIPL*ROPL*QPL*HPL
                    
       
        FR1 =PHIGR*ROGR*QGR
        FR2 =PHIGR*ROGR*QGR*UGR+PHIGR*PGR*ENX
        FR3 =PHIGR*ROGR*QGR*VGR+PHIGR*PGR*ENY
        FR4 =PHIGR*ROGR*QGR*WGR+PHIGR*PGR*ENZ
        FR5 =PHIGR*ROGR*QGR*HGR
        FR6 =PHIPR*ROPR*QPR
        FR7 =PHIPR*ROPR*QPR*UPR+PHIPR*PPR*ENX
        FR8 =PHIPR*ROPR*QPR*VPR+PHIPR*PPR*ENY
        FR9 =PHIPR*ROPR*QPR*WPR+PHIPR*PPR*ENZ
        FR10=PHIPR*ROPR*QPR*HPR

     IF(SGL.GT.0.0D0)THEN

                FSTARA1=FL1*AREA
                FSTARA2=FL2*AREA
                FSTARA3=FL3*AREA
                FSTARA4=FL4*AREA
		FSTARA5=FL5*AREA

         ELSEIF(SGL.LE.0.0D0.AND.SGSTAR.GT.0.0D0)THEN

                FSTARA1=(FL1+SGL*(CNSTARL1-CNL1))*AREA
                FSTARA2=(FL2+SGL*(CNSTARL2-CNL2))*AREA
                FSTARA3=(FL3+SGL*(CNSTARL3-CNL3))*AREA
                FSTARA4=(FL4+SGL*(CNSTARL4-CNL4))*AREA
		FSTARA5=(FL5+SGL*(CNSTARL5-CNL5))*AREA

         ELSEIF(SGSTAR.LE.0.0D0.AND.SGR.GE.0.0D0)THEN

                FSTARA1=(FR1+SGR*(CNSTARR1-CNR1))*AREA
                FSTARA2=(FR2+SGR*(CNSTARR2-CNR2))*AREA
                FSTARA3=(FR3+SGR*(CNSTARR3-CNR3))*AREA
                FSTARA4=(FR4+SGR*(CNSTARR4-CNR4))*AREA
		FSTARA5=(FR5+SGR*(CNSTARR5-CNR5))*AREA

         ELSEIF(SGR.LT.0.0D0)THEN

                FSTARA1=FR1*AREA
                FSTARA2=FR2*AREA
                FSTARA3=FR3*AREA
                FSTARA4=FR4*AREA
		FSTARA5=FR5*AREA

         ENDIF

         IF(SPL.GT.0.0D0)THEN
		
		FSTARA6=FL6*AREA
                FSTARA7=FL7*AREA
                FSTARA8=FL8*AREA
                FSTARA9=FL9*AREA
		FSTARA10=FL10*AREA

         ELSEIF(SPL.LE.0.0D0.AND.SPSTAR.GT.0.0D0)THEN

		FSTARA6=(FL6+SPL*(CNSTARL6-CNL6))*AREA
                FSTARA7=(FL7+SPL*(CNSTARL7-CNL7))*AREA
                FSTARA8=(FL8+SPL*(CNSTARL8-CNL8))*AREA
                FSTARA9=(FL9+SPL*(CNSTARL9-CNL9))*AREA
		FSTARA10=(FL10+SPL*(CNSTARL10-CNL10))*AREA
 !              if(myid==1.and.i==1)print*,"l2", FL9,SPL,CNSTARL9,CNL9,AREA
         ELSEIF(SPSTAR.LE.0.0D0.AND.SPR.GE.0.0D0)THEN

		FSTARA6=(FR6+SPR*(CNSTARR6-CNR6))*AREA
                FSTARA7=(FR7+SPR*(CNSTARR7-CNR7))*AREA
                FSTARA8=(FR8+SPR*(CNSTARR8-CNR8))*AREA
                FSTARA9=(FR9+SPR*(CNSTARR9-CNR9))*AREA
		FSTARA10=(FR10+SPR*(CNSTARR10-CNR10))*AREA

         ELSEIF(SPR.LT.0.0D0)THEN

		FSTARA6=FR6*AREA
                FSTARA7=FR7*AREA
                FSTARA8=FR8*AREA
                FSTARA9=FR9*AREA
		FSTARA10=FR10*AREA
         ENDIF
!if (myid.eq.2.and.i.eq.11) print*,"FSTARTA", I, iface, FSTARA1, FSTARA3, FSTARA4, FSTARA5
!  if (myid.eq.2.and.i.eq.11) print*,"pre RHS", I, iface, RHS(I,1), RHS(I,3), RHS(I,4), RHS(I,5)

!!$acc atomic
               RHS(I,1) =RHS(I,1)+FSTARA1
!!$acc atomic                
                RHS(I,2) =RHS(I,2)+FSTARA2
!!$acc atomic
                RHS(I,3) =RHS(I,3)+FSTARA3
!!$acc atomic
                RHS(I,4) =RHS(I,4)+FSTARA4
!!$acc atomic
                RHS(I,5) =RHS(I,5)+FSTARA5
!!$acc atomic
                RHS(I,6) =RHS(I,6)+FSTARA6
!!$acc atomic
                RHS(I,7) =RHS(I,7)+FSTARA7
!!$acc atomic 
                RHS(I,8) =RHS(I,8)+FSTARA8
!!$acc atomic update
                RHS(I,9) =RHS(I,9)+FSTARA9
!!$acc atomic
                RHS(I,10)=RHS(I,10)+FSTARA10
!     if(myid==1.and.i==1)print*, RHS(I, 9), iface,FSTARA9  
    endif
   enddo
  enddo

!$acc end parallel loop 

CALL source

return
end subroutine hllcbatten

