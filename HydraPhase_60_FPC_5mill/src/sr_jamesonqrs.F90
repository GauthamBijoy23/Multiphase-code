subroutine jamesonqrs
    USE GlobalVariables
    IMPLICIT NONE
    external geometry

  INTEGER :: i, j, iface, n1, n2, n3, neigh, nt, ns
  INTEGER :: np, ng
  REAL(DP) :: phitotal, area, enx, eny, enz
  REAL(DP) :: xf, yf, zf
  REAL(DP) :: roqL, uqL, vqL, wqL, phiqL, roqR, uqR, vqR, wqR, phiqR
  REAL(DP) :: rorL, urL, vrL, wrL, phirL, rorR, urR, vrR, wrR, phirR
  REAL(DP) :: rosL, usL, vsL, wsL, phisL, rosR, usR, vsR, wsR, phisR
  REAL(DP) :: qqnl, qqnr, phiqhalf, roqhalf, uqhalf, vqhalf, wqhalf, qqnhalf
  REAL(DP) :: rpqf, ruqf, rvqf, rwqf
  REAL(DP) :: qrnl, qrnr, phirhalf, rorhalf, urhalf, vrhalf, wrhalf, qrnhalf
  REAL(DP) :: rprf, rurf, rvrf, rwrf
  REAL(DP) :: qsnl, qsnr, phishalf, roshalf, ushalf, vshalf, wshalf, qsnhalf
  REAL(DP) :: rpsf, rusf, rvsf, rwsf
  REAL(DP) :: scx, scy, scz


!$acc parallel loop present(neles,phip(:),phiq(:),phir(:),phig(:),cn(:,:),phis(:), & 
!$acc us(:),vs(:),uq(:),vq(:),wq(:),ur(:),vr(:),ws(:),wr(:),epslnmin,rhoq,rhos,rhor) &
!$acc private(phitotal,i)

    do i=1,neles
   
IF (iblank(i) /= 1) CYCLE  !------ONLY PROCEED IF IBLANK = 1 ----------------(mod3)

	!if(myid==1)Print*,"phiq",i, cn(i,11), rhoq, epslnmin
   	
	phiq(i)=dmin1(dmax1(cn(i,11)/rhoq,epslnmin),1.0d0-epslnmin)
	phir(i)=dmin1(dmax1(cn(i,15)/rhor,epslnmin),1.0d0-epslnmin)
    	phis(i)=dmin1(dmax1(cn(i,19)/rhos,epslnmin),1.0d0-epslnmin)
	phitotal=phig(i)+phip(i)+phiq(i)+phir(i)+phis(i)
	phig(i)=phig(i)/phitotal
	phip(i)=phip(i)/phitotal
	phiq(i)=phiq(i)/phitotal
	phir(i)=phir(i)/phitotal
	phis(i)=phis(i)/phitotal
	
	!if(myid==1)Print*,i,phis(i)
	   
!	---------------------------------------------------------------		
	 uq(i)=cn(i,12)/cn(i,11)
	 vq(i)=cn(i,13)/cn(i,11)
	 wq(i)=cn(i,14)/cn(i,11)
	 ur(i)=cn(i,16)/cn(i,15)
	 vr(i)=cn(i,17)/cn(i,15)
	 wr(i)=cn(i,18)/cn(i,15)
	 us(i)=cn(i,20)/cn(i,19)
	 vs(i)=cn(i,21)/cn(i,19)
	 ws(i)=cn(i,22)/cn(i,19)
		
	cn(i,11)=phiq(i)*rhoq
        cn(i,12)=phiq(i)*rhoq*uq(i)
        cn(i,13)=phiq(i)*rhoq*vq(i)
        cn(i,14)=phiq(i)*rhoq*wq(i)

	cn(i,15)=phir(i)*rhor
        cn(i,16)=phir(i)*rhor*ur(i)
        cn(i,17)=phir(i)*rhor*vr(i)
        cn(i,18)=phir(i)*rhor*wr(i)

	cn(i,19)=phis(i)*rhos
        cn(i,20)=phis(i)*rhos*us(i)
        cn(i,21)=phis(i)*rhos*vs(i)
        cn(i,22)=phis(i)*rhos*ws(i)
     enddo
!$acc end parallel loop     
!if(myid==1) Print*,"loop 1 jamesonqrs" 


!$acc parallel loop present(nghosts,ws(:),wr(:),vr(:),ur(:),wq(:),vq(:),uq(:),wp(:), &
!$acc vs(:),vp(:),us(:),up(:),phis(:),ntype(:),phir(:),phiq(:),phip(:),phig(:), &
!$acc nparent(:),cn(:,:),nghost(:),epslnmin,rhoq,rhos,rhor) private(ng,j,nt,np)
 	
      do i=1,nghosts
	np=nparent(i)
        ng=nghost(i)
        nt=ntype(i)
        ns=nside(i)

!$acc loop seq		
	do j = 11, 22
         cn(ng,j)=cn(np,j)
        enddo
         
        if(nt.eq.2)then
        phip(ng)=epslnmin
        phiq(ng)=epslnmin
        phir(ng)=epslnmin
        phis(ng)=epslnmin
        phig(ng)=1.0d0-phip(ng)-phiq(ng)-phir(ng)-phis(ng)

	uq(ng)=up(ng)
	vq(ng)=vp(ng)
	wq(ng)=wp(ng)
	ur(ng)=up(ng)
	vr(ng)=vp(ng)
	wr(ng)=wp(ng)
	us(ng)=up(ng)
	vs(ng)=vp(ng)
	ws(ng)=wp(ng)
!	-------------------------------------------------------------------------
        elseif(nt.eq.3)then
 
        phip(ng)=epslnmin
        phiq(ng)=phiq(np)
        phir(ng)=phir(np)
        phis(ng)=phis(np)
        phig(ng)=1.0d0-phip(ng)-phiq(ng)-phir(ng)-phis(ng)

	uq(ng)=uq(np)
	vq(ng)=vq(np)
	wq(ng)=wq(np)
	ur(ng)=ur(np)
	vr(ng)=vr(np)
	wr(ng)=wr(np)
	us(ng)=us(np)
	vs(ng)=vs(np)
	ws(ng)=ws(np)

        elseif(nt.eq.5)then
         
        phiq(ng)=phiq(np)
        phir(ng)=phir(np)
        phis(ng)=phis(np)
	uq(ng)=-uq(np)
	vq(ng)=-vq(np)
	wq(ng)=-wq(np)
	ur(ng)=-ur(np)
	vr(ng)=-vr(np)
	wr(ng)=-wr(np)
	us(ng)=-us(np)
	vs(ng)=-vs(np)
	ws(ng)=-ws(np)

        endif

	cn(ng,11)=phiq(ng)*rhoq
        cn(ng,12)=phiq(ng)*rhoq*uq(ng)
        cn(ng,13)=phiq(ng)*rhoq*vq(ng)
        cn(ng,14)=phiq(ng)*rhoq*wq(ng)

	cn(ng,15)=phir(ng)*rhor
        cn(ng,16)=phir(ng)*rhor*ur(ng)
        cn(ng,17)=phir(ng)*rhor*vr(ng)
        cn(ng,18)=phir(ng)*rhor*wr(ng)

	cn(ng,19)=phis(ng)*rhos
        cn(ng,20)=phis(ng)*rhos*us(ng)
        cn(ng,21)=phis(ng)*rhos*vs(ng)
        cn(ng,22)=phis(ng)*rhos*ws(ng)
     enddo
!$acc end parallel loop     

!if(myid==1) Print*,"loop 2 jamson" 

!$acc parallel loop present(neles,phis(:),sc4y(:),sc4x(:),sc3z(:),nc3(:), &
!$acc nc2(:),sc2z(:),sc3x(:),sc3y(:),nc1(:),sc1z(:),sc2x(:),sc2y(:),rhs(:,11:22), &
!$acc sc1x(:),sc1y(:),nc4(:),sc4z(:),us(:),vs(:),ntype(:),uq(:),vq(:),wq(:),phiq(:), &
!$acc ur(:),vr(:),ws(:),wr(:),phir(:),rhoq,rhor,rhos) & 
!$acc private(enz,n1,n2,n3,rwsf,scx,scy,area,enx,eny,roqhalf,uqhalf, &
!$acc vqhalf,wqhalf,phiqhalf,roql,uql,vql,wql,phiql,rorhalf,urhalf,vrhalf, &
!$acc wrhalf,phirhalf,rorl,url,vrl,wrl,phirl,roshalf,ushalf,vshalf,wshalf, & 
!$acc phishalf,rosl,usl,vsl,wsr,wsl,phisl,qqnhalf,qqnl,qsnr,scz,usr,vsr,phisr, &
!$acc rosr,rpsf,rusf,rvsf,qrnhalf,qrnl,phiqr,roqr,uqr,vqr,wqr,qqnr,rpqf,ruqf, &
!$acc rvqf,rwqf,qsnhalf,qsnl,phirr,rorr,urr,vrr,wrr,qrnr,rprf,rurf,rvrf,rwrf,nt,neigh)
		
    do i=1,neles
IF (iblank(i) /=1)CYCLE !--------------------- skip if iblank /=1 (mod3)
	          rhs(i,11)=0.0d0
		  rhs(i,12)=0.0d0
		  rhs(i,13)=0.0d0
		  rhs(i,14)=0.0d0
		  rhs(i,15)=0.0d0
		  rhs(i,16)=0.0d0
		  rhs(i,17)=0.0d0
		  rhs(i,18)=0.0d0
		  rhs(i,19)=0.0d0
		  rhs(i,20)=0.0d0
		  rhs(i,21)=0.0d0
		  rhs(i,22)=0.0d0
!$acc loop seq			  
	do iface=1,4

	   IF(IFACE.EQ.1)THEN
             n1 = nod(I,2)
             n2 = nod(I,4)
             n3 = nod(I,3)
             neigh = NC1(I)
             SCX = SC1X(I)
             SCY = SC1Y(I)
             SCZ = SC1Z(I)
           ELSEIF(IFACE.EQ.2)THEN
             n1 = nod(I,3)
             n2 = nod(I,4)
             n3 = nod(I,1)
             neigh = NC2(I)
             SCX = SC2X(I)
             SCY = SC2Y(I)
             SCZ = SC2Z(I)
           ELSEIF(IFACE.EQ.3)THEN
             n1 = nod(I,4)
             n2 = nod(I,2)
             n3 = nod(I,1)
             neigh = NC3(I)
             SCX = SC3X(I)
             SCY = SC3Y(I)
             SCZ = SC3Z(I)
           ELSEIF(IFACE.EQ.4)THEN
             n1 = nod(I,2)
             n2 = nod(I,3)
             n3 = nod(I,1)
             neigh = NC4(I)
             SCX = SC4X(I)
             SCY = SC4Y(I)
             SCZ = SC4Z(I)
          ENDIF
 
            XF  = (X(N1)+X(N2)+X(N3))*0.3333D0
            YF  = (Y(N1)+Y(N2)+Y(N3))*0.3333D0
            ZF  = (Z(N1)+Z(N2)+Z(N3))*0.3333D0

            AREA=DSQRT(SCX*SCX+SCY*SCY+SCZ*SCZ)

            ENX = SCX/AREA
            ENY = SCY/AREA
            ENZ = SCZ/AREA
!	-----------------------------------------------------------------------------
	NT=0
	IF(NEIGH.GT.NELES)NT=NTYPE(NEIGH-NElES)
     
	    roqL=rhoq
            uqL=uq(i)
	    vqL=vq(i)
	    wqL=wq(i)
	    phiqL=phiq(i)

	    roqR=rhoq
	    uqR=uq(neigh)
	    vqR=vq(neigh)
	    wqR=wq(neigh)
	    phiqR=phiq(neigh)
!	-------------------------------------------------------------------------------------------------
	   
	   rorL=rhor
	   urL=ur(i)
	   vrL=vr(i)
	   wrL=wr(i)
	   phirL=phir(i)

          rorR=rhor
	  urR=ur(neigh)
	  vrR=vr(neigh)
	  wrR=wr(neigh)
	  phirR=phir(neigh)
!	-------------------------------------------------------------------------------------------------
	
	  rosL=rhos
	  usL=us(i)
	  vsL=vs(i)
	  wsL=ws(i)
	  phisL=phis(i)

	  rosR=rhos
	  usR=us(neigh)
	  vsR=vs(neigh)
	  wsR=ws(neigh)
	  phisR=phis(neigh)
				 
      IF(NT.NE.5)THEN
				 
	QQNL=UQL*ENX+VQL*ENY+WQL*ENZ
	QQNR=UQR*ENX+VQR*ENY+WQR*ENZ
	PHIQHALF=(PHIQL+PHIQR)*0.5D0
        ROQHALF=(ROQL+ROQR)*0.5D0
        UQHALF=(UQL+UQR)*0.5D0
        VQHALF=(VQL+VQR)*0.5D0
        WQHALF=(WQL+WQR)*0.5D0
        QQNHALF=(QQNL+QQNR)*0.5D0

        RPQF= (PHIQHALF*ROQHALF*QQNHALF)*AREA
        RUQF= (PHIQHALF*ROQHALF*QQNHALF*UQHALF)*AREA
        RVQF= (PHIQHALF*ROQHALF*QQNHALF*VQHALF)*AREA
        RWQF= (PHIQHALF*ROQHALF*QQNHALF*WQHALF)*AREA		
!	---------------------------------------------------------	
	QRNL=URL*ENX+VRL*ENY+WRL*ENZ
	QRNR=URR*ENX+VRR*ENY+WRR*ENZ
	PHIRHALF=(PHIRL+PHIRR)*0.5D0
        RORHALF=(RORL+RORR)*0.5D0
        URHALF=(URL+URR)*0.5D0
        VRHALF=(VRL+VRR)*0.5D0
        WRHALF=(WRL+WRR)*0.5D0
        QRNHALF=(QRNL+QRNR)*0.5D0

        RPRF= (PHIRHALF*RORHALF*QRNHALF)*AREA
        RURF= (PHIRHALF*RORHALF*QRNHALF*URHALF)*AREA
        RVRF= (PHIRHALF*RORHALF*QRNHALF*VRHALF)*AREA
        RWRF= (PHIRHALF*RORHALF*QRNHALF*WRHALF)*AREA
!	---------------------------------------------------------	
	QSNL=USL*ENX+VSL*ENY+WSL*ENZ
	QSNR=USR*ENX+VSR*ENY+WSR*ENZ
	PHISHALF=(PHISL+PHISR)*0.5D0
        ROSHALF=(ROSL+ROSR)*0.5D0
        USHALF=(USL+USR)*0.5D0
        VSHALF=(VSL+VSR)*0.5D0
        WSHALF=(WSL+WSR)*0.5D0
        QSNHALF=(QSNL+QSNR)*0.5D0

        RPSF= (PHISHALF*ROSHALF*QSNHALF)*AREA
        RUSF= (PHISHALF*ROSHALF*QSNHALF*USHALF)*AREA
        RVSF= (PHISHALF*ROSHALF*QSNHALF*VSHALF)*AREA
        RWSF= (PHISHALF*ROSHALF*QSNHALF*WSHALF)*AREA
!	---------------------------------------------------------

        RHS(I,11)=RHS(I,11)+RPQF
        RHS(I,12)=RHS(I,12)+RUQF
        RHS(I,13)=RHS(I,13)+RVQF
        RHS(I,14)=RHS(I,14)+RWQF
        RHS(I,15)=RHS(I,15)+RPRF
        RHS(I,16)=RHS(I,16)+RURF
        RHS(I,17)=RHS(I,17)+RVRF
        RHS(I,18)=RHS(I,18)+RWRF
        RHS(I,19)=RHS(I,19)+RPSF
        RHS(I,20)=RHS(I,20)+RUSF
        RHS(I,21)=RHS(I,21)+RVSF		
        RHS(I,22)=RHS(I,22)+RWSF		
       ENDIF
		
      enddo		
     enddo
!$acc end parallel loop     
!if(myid==1) Print*,"loop 3 jamsonqrs" 	    
    call sourceqrs

end subroutine jamesonqrs
