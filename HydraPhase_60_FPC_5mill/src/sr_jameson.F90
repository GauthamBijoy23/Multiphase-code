SUBROUTINE jameson
    USE GlobalVariables
    IMPLICIT NONE
    external geometry
    external primitive
    external smooth

    INTEGER  :: i,neigh,ng,np,ns,nt,n1,n2,n3,n4,j,iface
    REAL(dp) :: eg,ep,zf,yf,xf,wpr,wpl,wgr,wgl,vpr,vpl,vgr,vgl,upr,upl,ugr,ugl,tpr,tpl,tgr,tgl,scz
    REAL(DP) :: scy,scx,ropr,ropl,rogl,qpr,qpl,rogr,qgr,qgl,prr,prl,ppr,ppl,pgl,pgr,phigl,phigr
    REAL(DP) :: phipl,phipr,hpr,hpl,hgr,hgl,epr,epl,enz,eny,enx,egr,egl,agr,apl,apr,area,agl
    REAL(DP) :: phighalf,roghalf,pghalf,ughalf,vghalf,wghalf,hghalf,qghalf,phiphalf,rophalf
    REAL(DP) :: pphalf,uphalf,vphalf,wphalf,hphalf,qphalf,rpgf,rugf,rvgf,rwgf,rtgf,rppf
    REAL(DP) :: rupf,rvpf,rwpf,rtpf,tkgL,tegL,tkgR,tegR,tkghalf,teghalf,rtkgf,rtegf
!For NP
    REAL(DP) :: enpL,enpR,flnp,frnp,fltk,frtk,flte,frte,fl1,fl2,fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10
    REAL(DP) :: fr1,fr2,fr3,fr4,fr5,fr6,fr7,fr8,fr9,fr10   
    
    REAL(DP) :: yfgL,yfgR,yogL,ypgL,ypgR,yogR,flyfg,flyog,flypg,fryfg,fryog,frypg
    REAL(DP) :: cpgL,rrgL,amolgL,cvgL,gamagL,cpgR,rrgR,amolgR,cvgR,gamagR
   
   CALL primitiveplusbc

#ifdef NP_P
    CALL mpiexcn
#endif
   
!#ifdef KE_TURB
!CALL mpiexdel
!#endif

!$acc parallel loop present(pp(:),phig(:),nod(:,:),tg(:),wg(:),vg(:),ug(:),phi(:),sc4z(:),vp(:),up(:),tp(:),ycel(:), &
!$acc xcel(:),wp(:),sc4y(:),sc4x(:),sc3z(:),nc3(:),nc2(:),sc2z(:),sc3x(:),sc3y(:),nc1(:),sc1z(:),sc2x(:),sc2y(:),nc4(:), &
!$acc rhs(:,:),sc1x(:),sc1y(:),zcel(:),z(:),y(:),x(:),pg(:),rrg,pinf,epslnmax,gamag,gamap,epslnmin,myorder,neles,cpp, &
!$acc ntype(:),rop(:),rog(:),tkg(:),teg(:),nkg,neg,npp,enp(:),nyfg,nyog,nypg,yfg(:),yog(:),ypg(:),nghosts) &  
!$acc private(scx,scy,wpr,apr,egr,enx,eny,tpr,upr,vpr,ugr,vgr,wgr,tgr,neigh,ugl,vgl,wgl,tgl,phigl,upl,vpl,wpl,tpl,pgl, &
!$acc ppl,hpr,phigr,xf,yf,zf,phipl,epr,qpr,phipr,rogl,ropl,egl,hgl,qgl,agl,enz,hgr,qgr,agr,epl,hpl,qpl,apl,scz,rogr, &
!$acc ropr,pgr,ppr,area,n3,n2,n1,iface,nt,pghalf,qphalf,hphalf,rwpf,uphalf,vphalf,pphalf,phighalf,roghalf,phiphalf, &
!$acc ughalf,vghalf,wphalf,wghalf,hghalf,qghalf,rophalf,rtpf,rupf,rvpf,rppf,rpgf,rugf,rvgf,rwgf,rtgf,tkgL,tegL,tkgR, &
!$acc tegR,TKGHALF,TEGHALF,RTKGF,RTEGF,enpL,enpR,FLNP,FRNP,fl1,fl2,fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10,fr1,fr2,fr3,fr4, &
!$acc fr5,fr6,fr7,fr8,fr9,fr10,yfgL,yogR,ypgL,yfgR,yogR,ypgR,frypg,yogl,flypg,flyfg,fryfg,flyog,fryog,fltk,frtk,flte,frte)

  do i = 1, neles
if(iblank(i) /=1)cycle !=======(mod3) Iblank check

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
#ifdef KE_TURB
   rhs(i, nkg) = 0.0d0
   rhs(i, neg) = 0.0d0
#endif  
#ifdef NP_P
   rhs(i, npp) = 0.0d0
#endif 
#ifdef YP_P
   rhs(i, nyfg) = 0.0d0
   rhs(i, nyog) = 0.0d0
   rhs(i, nypg) = 0.0d0
#endif   
  
!$acc loop seq  
!!$acc loop unroll   
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
              
!=====================================================================        
! Store neighbour as right cells; before if loop; since uR, vR, HR, pR required in rhs
!=====================================================================


    nt = 0
    if ((neigh > neles).and.(neigh <= (neles+nghosts))) nt = ntype(neigh - neles) !for boundary ghosts 
!    nt = 0
!    if (neigh > neles) nt = ntype(neigh - neles) 
! Neighbour is greater than the cell then boundary ghost only 

    rogl = rog(i)
    ugl = ug(i)
    vgl = vg(i)
    wgl = wg(i)
    tgl = tg(i)
    phigl = phig(i)
    ropl = rop(i)
    upl = up(i)
    vpl = vp(i)
    wpl = wp(i)
    tpl = tp(i)
    pgl = pg(i)
    ppl = pp(i)

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
#ifdef KE_TURB
    tkgL = tkg(i)
    tegL = teg(i)
    tkgR = tkg(neigh)
    tegR = teg(neigh)
#endif

#ifdef NP_P
    enpL = enp(i)
    enpR = enp(neigh)
#endif

#ifdef YP_P
    yfgL = yfg(i)
    yogL = yog(i)
    ypgL = ypg(i)
    yfgR = yfg(neigh)
    yogR = yog(neigh)
    ypgR = ypg(neigh)
#endif

!=====================================================================
    PhipL = 1.0d0 - PhigL
    PhipR = 1.0d0 - PhigR
!=====================================================================

     IF(NT.EQ.5)THEN
     RHS(I, 2) = RHS(I, 2) + PHIGL * PGL * ENX * AREA
     RHS(I, 3) = RHS(I, 3) + PHIGL * PGL * ENY * AREA
     RHS(I, 4) = RHS(I, 4) + PHIGL * PGL * ENZ * AREA
     RHS(I, 7) = RHS(I, 7) + PHIPL * PPL * ENX * AREA
     RHS(I, 8) = RHS(I, 8) + PHIPL * PPL * ENY * AREA
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
    
!=====================================================================

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
       	

#ifdef KE_TURB
                FLTK=PHIGL*ROGL*QGL*TKGL
		FRTK=PHIGR*ROGR*QGR*TKGR
		
		FLTE=PHIGL*ROGL*QGL*TEGL
		FRTE=PHIGR*ROGR*QGR*TEGR
#endif
#ifdef NP_P
                FLNP=ENPL*QPL
		FRNP=ENPR*QPR
#endif
#ifdef YP_P
                FLYFG=PHIGL*ROGL*QGL*YFGL
		FRYFG=PHIGR*ROGR*QGR*YFGR
		
		FLYOG=PHIGL*ROGL*QGL*YOGL
		FRYOG=PHIGR*ROGR*QGR*YOGR
		
		FLYPG=PHIGL*ROGL*QGL*YPGL
		FRYPG=PHIGR*ROGR*QGR*YPGR
#endif

                RHS(I,1) =RHS(I,1)+((FL1+FR1)*0.5D0*AREA)
                RHS(I,2) =RHS(I,2)+((FL2+FR2)*0.5D0*AREA)
                RHS(I,3) =RHS(I,3)+((FL3+FR3)*0.5D0*AREA)
                RHS(I,4) =RHS(I,4)+((FL4+FR4)*0.5D0*AREA)
                RHS(I,5) =RHS(I,5)+((FL5+FR5)*0.5D0*AREA)
                RHS(I,6) =RHS(I,6)+((FL6+FR6)*0.5D0*AREA)
                RHS(I,7) =RHS(I,7)+((FL7+FR7)*0.5D0*AREA)
                RHS(I,8) =RHS(I,8)+((FL8+FR8)*0.5D0*AREA)
                RHS(I,9) =RHS(I,9)+((FL9+FR9)*0.5D0*AREA)
                RHS(I,10)=RHS(I,10)+((FL10+FR10)*0.5D0*AREA)
				
#ifdef KE_TURB
                RHS(I,NKG) =RHS(I,NKG)+((FLTK+FRTK)*0.5D0*AREA)
                RHS(I,NEG) =RHS(I,NEG)+((FLTE+FRTE)*0.5D0*AREA)				
#endif 
  
#ifdef NP_P
                RHS(I,NPP) =RHS(I,NPP)+((FLNP+FRNP)*0.5D0*AREA)
#endif

#ifdef YP_P
                RHS(I,NYFG) =RHS(I,NYFG)+((FLYFG+FRYFG)*0.5D0*AREA)
                RHS(I,NYOG) =RHS(I,NYOG)+((FLYOG+FRYOG)*0.5D0*AREA)				
                RHS(I,NYPG) =RHS(I,NYPG)+((FLYPG+FRYPG)*0.5D0*AREA)
#endif            
    endif
   enddo
  enddo
!$acc end parallel loop 

CALL source

return
end subroutine jameson

