SUBROUTINE viscousgp
    USE GlobalVariables
    USE MatrixOps
    IMPLICIT NONE
    external geometry
    external primitive
    external smooth
    
    INTEGER  :: i,np,ng,n1,n2,n3,n4,m1,m2,m3,m4,m5,m6,m7,m8,iface,neigh,ii
    REAL(DP) :: xf,yf,zf,x1,x2,x3,x4,upf,ugf,vpf,vgf,wpf,wgf,tzzp,tzzg,tyzp,tyzg
    REAL(DP) :: tyyp,tyyg,txzp,txzg,txyp,txyg,txxp,txxg,tpf,tgf,scx,scy,scz,ropf,rogf
    REAL(DP) :: qzzp,qzzg,qyyp,qyyg,qxxp,qxxg,phipf,phigf,dwpx,dwpy,dwpz,dwgx,dwgy,dwgz
    REAL(DP) :: dvpx,dvpy,dvpz,dvgx,dvgy,dvgz,dupx,dupy,dupz,dugx,dugy,dugz
    REAL(DP) :: dtpx,dtpy,dtpz,dtgx,dtgy,dtgz,dnr,delvp,delvg,bwp1,bwp2,bwp3
    REAL(DP) :: bwg1,bwg2,bwg3,bvp1,bvp2,bvp3,bvg1,bvg2,bvg3,bup1,bup2,bup3
    REAL(DP) :: bug1,bug2,bug3,btg1,btg2,btg3,btp1,btp2,btp3,akg,amug
    REAL(DP) :: a11,a12,a13,a22,a23,a33,xcelm1f,xcelm2f,xcelm3f,xcelm4f,xcelm5f
    REAL(DP) :: ycelm1f,ycelm2f,ycelm3f,ycelm4f,ycelm5f,zcelm1f,zcelm2f,zcelm3f,zcelm4f,zcelm5f 
    REAL(DP) :: ugm1f,ugm2f,ugm3f,ugm4f,ugm5f,vgm1f,vgm2f,vgm3f,vgm4f,vgm5f
    REAL(DP) :: wgm1f,wgm2f,wgm3f,wgm4f,wgm5f,upm1f,upm2f,upm3f,upm4f,upm5f
    REAL(DP) :: wpm1f,wpm2f,wpm3f,wpm4f,wpm5f,tpm1f,tpm2f,tpm3f,tpm4f,tpm5f
    REAL(DP) :: tgm1f,tgm2f,tgm3f,tgm4f,tgm5f,vpm1f,vpm2f,vpm3f,vpm4f,vpm5f
    REAL(DP) :: xcelm6f,xcelm7f,xcelm8f,ycelm6f,ycelm7f,ycelm8f,zcelm6f,zcelm7f,zcelm8f
    REAL(DP) :: ugm6f,ugm7f,ugm8f,vgm6f,vgm7f,vgm8f,wgm6f,wgm7f,wgm8f,tgm6f,tgm7f,tgm8f
    REAL(DP) :: upm6f,upm7f,upm8f,vpm6f,vpm7f,vpm8f,wpm6f,wpm7f,wpm8f,tpm6f,tpm7f,tpm8f
    REAL(DP) :: wt,wu,alphat,deltat,alphau,fract,fracu,akp,amup,gpart1,gpart2,gamapdot
    REAL(DP) :: btkg1,btkg2,btkg3,bteg1,bteg2,bteg3,tkgf,tegf,tkgm1f,tkgm2f,tkgm3f,tkgm4f,tkgm5f,tegm1f,tegm2f
    REAL(DP) :: tegm3f,tegm4f,tegm5f,dtkgx,dtkgy,dtkgz,dtegx,dtegy,dtegz,amutg,aktg,fvtkg,gvtkg,hvtkg
    REAL(DP) :: fvteg,gvteg,hvteg,tkgm6f,tkgm7f,tkgm8f,tegm6f,tegm7f,tegm8f,yfgf,yogf,ypgf,yfgm1f,yfgm2f
    REAL(DP) :: yfgm3f,yfgm4f,yfgm5f,yogm1f,yogm2f,yogm3f,yogm4f,yogm5f,ypgm1f,ypgm2f,ypgm3f,ypgm4f,ypgm5f
    REAL(DP) :: byfg1,byfg2,byfg3,byog1,byog2,byog3,bypg1,bypg2,bypg3,dyfgx,dyfgy,dyfgz,dyogx,dyogy,dyogz
    REAL(DP) :: dypgx,dypgy,dypgz,cvg,yfgm6f,yfgm7f,yfgm8f,yogm6f,yogm7f,yogm8f,ypgm6f,ypgm7f,ypgm8f

!======================
!Initializing variables
!======================
!$acc parallel loop present(vistp(:),visup(:),visvp(:),visug(:), &
!$acc visvg(:),viswp(:),viswg(:),vistg(:),NELES,visteg(:),vistkg(:),visyfg(:),visyog(:),visypg(:))  
   do i = 1, neles
if(iblank(i)/=1)cycle !-----------(mod3)
    visug(i) = 0.0d0
    visvg(i) = 0.0d0
    viswg(i) = 0.0d0
    vistg(i) = 0.0d0
    visup(i) = 0.0d0
    visvp(i) = 0.0d0
    viswp(i) = 0.0d0
    vistp(i) = 0.0d0
    
#ifdef KE_TURB
    vistkg(i) = 0.0d0
    visteg(i) = 0.0d0
#endif

#ifdef YP_P
   visyfg(i)= 0.0D0
   visyog(i)= 0.0D0
   visypg(i)= 0.0D0
#endif	
  end do
!$acc end parallel loop  


!=====================================================================================================================================
!                                                                             Boundary condition loop
!=====================================================================================================================================

!=========================
!Subsonic Inflow Condition
!=========================
!$acc parallel loop present(bc2_list(:),n_bc2,phip(:),phipzero,phig(:),ug(:),vg(:),wg(:),ugmin,vgmin,wgmin, &
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp, &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),nparent(:),nghost(:))  &
!$acc private(ii,np,ng,i)

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
                    
#ifdef NP_P
    enp(ng)=enp(np)
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
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,dd(:), &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),nparent(:),nghost(:),xmax,total,pmin,dt)  &
!$acc private(i,ii,np,ng,fract,alphat,wt,fracu,alphau,wu,amug,cvg)

do ii = 1, n_bc3
  i = bc3_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i) 
  
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

#ifdef NP_P
    enp(ng)=enp(np)
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
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp, &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),nparent(:),nghost(:))  &
!$acc private(ii,np,ng,i)

do ii = 1, n_bc5
  i = bc5_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i)
  
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
        rog(ng) = pg(ng) / (rrg * tg(ng))
        rop(ng) = gamap * (pp(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))

#ifdef NP_P
    enp(ng)=enp(np)
#endif
        
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
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,dd(:), &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),y(:),nod(:,:),nparent(:),nghost(:),ymax,ymin,total,pmin,dt)&
!$acc private(i,ii,np,ng,fract,alphat,wt,fracu,alphau,wu,amug,cvg)
do ii = 1, n_bc51
  i = bc51_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i)
  
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
!WT        

       
         wT = (vg(np))
         wU = (vg(np)) + sqrt(gamag*(pg(np))/(rog(np)))
         
         alphaT = wT*dt/dd(np)  
         alphaU = wU*dt/dd(np)
        
         fracT = 1/(1+ alphaT)
         fracU = 1/(1+ alphaU) 
         
     !    if(total.eq.1) then
     !      ug(ng)=ug(np)
     !      vg(ng)=-vg(np)
     !      wg(ng)=wg(np)
     !      tg(ng)=tg(np)  
         !  pg(ng)=pmin
      !   endif
         
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
 
        phig(ng) = phig(np)
 !       ug(ng) = ug(np)
 !       vg(ng) = -vg(np)
 !       wg(ng) = wg(np)
 !       tg(ng) = tg(np)

        phip(ng) = phip(np)
        up(ng) = up(np)
        vp(ng) = -vp(np)
        wp(ng) = wp(np) 
        tp(ng) = tp(np)
        pg(ng) = pg(np)
        pp(ng) = pp(np)
        rog(ng) = pg(ng) / (rrg * tg(ng))
        rop(ng) = gamap * (pp(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))

#ifdef NP_P
    enp(ng)=enp(np)
#endif
        
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
!$acc tg(:),tgmin,pg(:),rog(:),rrg,up(:),vp(:),wp(:),tp(:),pp(:),rop(:),gamap,pinf,cpp,dd(:), &
!$acc phiq(:),epslnmin,phir(:),phis(:),phip(:),enp(:),yfg(:),yog(:),ypg(:),tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg, &
!$acc amolog,amolpg,amolg,cpg,gamag,nyfg,nyog,nypg,tkg(:),teg(:),z(:),nod(:,:),nparent(:),nghost(:),zmax,zmin,total,pmin,dt)&
!$acc private(i,ii,np,ng,fract,alphat,wt,fracu,alphau,wu,amug,cvg)
do ii = 1, n_bc53
  i = bc53_list(ii)

  np = NPARENT(i)
  ng = NGHOST(i) 
  
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
!WT        

       
         wT = (wg(np))
         wU = (wg(np)) + sqrt(gamag*(pg(np))/(rog(np)))
         
         alphaT = wT*dt/dd(np)  
         alphaU = wU*dt/dd(np)
        
         fracT = 1/(1+ alphaT)
         fracU = 1/(1+ alphaU) 
         
     !    if(total.eq.1) then
     !      ug(ng)=ug(np)
     !      vg(ng)=vg(np)
     !      wg(ng)=-wg(np)
     !      tg(ng)=tg(np)  
         !  pg(ng)=pmin
     !    endif
         
         if(wg(np).gt.wg(ng))then                                   
           wg(ng)=-wg(ng)!* fracU + (1- fracU)* (2.0d0*wg(np) - wg(ng))
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
      
        phig(ng) = phig(np)
   !     ug(ng) = ug(np)
   !     vg(ng) = vg(np)
   !     wg(ng) = -wg(np)
   !     tg(ng) = tg(np)

        phip(ng) = phip(np)
        up(ng) = up(np)
        vp(ng) = vp(np)
        wp(ng) = -wp(np) 
        tp(ng) = tp(np)
        pg(ng) = pg(np)
        pp(ng) = pp(np)
        rog(ng) = pg(ng) / (rrg * tg(ng))
        rop(ng) = gamap * (pp(ng) + pinf) / ((gamap - 1.0d0) * cpp * tp(ng))

#ifdef NP_P
    enp(ng)=enp(np)
#endif

#ifdef KE_TURB
        tkg(ng)=tkg(np)
	teg(ng)=teg(np)
#endif	
 end do 
!$acc end parallel loop    
  
     

!$acc parallel loop present(neles,nc4(:),c23,prg,cpg,cpp,prp, &
!$acc viswg(:),visug(:),visvp(:),vistp(:),y(:),phig(:),wg(:),ug(:), &
!$acc nod(:,:),vp(:),tp(:),sc4y(:),sc3z(:),sc4x(:),sc3y(:),sc2z(:), &
!$acc sc3x(:),sc2y(:),sc1z(:),sc2x(:),xcel(:),sc1y(:),phip(:),sc1x(:), &
!$acc wp(:),ycel(:),up(:),sc4z(:),viswp(:),tg(:),x(:),zcel(:),z(:), &
!$acc visup(:),vg(:),visvg(:),vistg(:),nc2(:),nc1(:),nc3(:),as,ts, &
!$acc tolkg,tkg(:),teg(:),rog(:),vistkg(:),visteg(:), &
!$acc yfg,yog,ypg,tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg,amolog,amolpg,amolg,gamag,cmug, &
!$acc yfg(:),yog(:),ypg(:),yfg,yog,ypg,tolsp,cpg,cpfg,cpog,cppg,rrg,amolfg,amolog,amolpg, &
!$acc amolg,gamag,cmug,visyfg,visyog,visypg,enp)  &
!$acc private(n1,n2,n3,n4,ropf,scx,scy,wpf,xf,zf,yf,m8,qzzp,tzzp,upf, &
!$acc vpf,scz,neigh,phipf,rogf,ugf,vgf,wgf,tgf,phigf,dwpz,m1,m2, &
!$acc m3,m4,akg, a13,a23,a11,a12,a22,btp3,bug1,bug2,bup3,bvg1,bvg2,bvp3, &
!$acc bwg1,bwg2,a33,btg1,btg2,bug3,bup1,bup2,bvg3,bvp1,bvp2,bwg3,bwp2,bwp1, &
!$acc delvp,dtpz,dugx,dugy,dupz,dvgx,dvgy,dvpz, btg3,btp1,btp2,&
!$acc dwgx,dwgy,dnr,dtgx,dugz,dupx,dupy,dvgz,dvpx,dvpy,dwgz,dwpx,dwpy, &
!$acc dtgz,dtpx,dtpy,bwp3,tpf,txzp,tyzp,txxp,txyp,tyyp,qxxp,qyyp,delvg, &
!$acc amug,dtgy,txxg,tyyg,tzzg,txyg,txzg,tyzg,qxxg,qyyg,qzzg,m5,m6,m7,xcelm1f, &
!$acc xcelm2f,xcelm3f,xcelm4f,xcelm5f,xcelm6f,xcelm7f,xcelm8f, &
!$acc ycelm1f,ycelm2f,ycelm3f,ycelm4f,ycelm5f, &
!$acc zcelm1f,zcelm2f,zcelm3f,zcelm4f,zcelm5f,ugm1f,ugm2f,ugm3f,ugm4f,ugm5f, &
!$acc wgm1f,wgm2f,wgm3f,wgm4f,wgm5f,upm1f,upm2f,upm3f,upm4f,upm5f, &
!$acc wpm1f,wpm2f,wpm3f,wpm4f,wpm5f,tpm1f,tpm2f,tpm3f,tpm4f,tpm5f, &
!$acc tgm1f,tgm2f,tgm3f,tgm4f,tgm5f,vpm1f,vpm2f,vpm3f,vpm4f,vpm5f, &
!$acc ugm6f,ugm7f,ugm8f,vgm6f,vgm7f,vgm8f,wgm6f,wgm7f,wgm8f,tgm6f,tgm7f,tgm8f, &
!$acc upm6f,upm7f,upm8f,vpm6f,vpm7f,vpm8f,wpm6f,wpm7f,wpm8f,tpm6f,tpm7f,tpm8f,&
!$acc vgm1f,vgm2f,vgm3f,vgm4f,vgm5f,ycelm6f,ycelm7f,ycelm8f,zcelm6f,zcelm8f,zcelm7f,&
!$acc akp,amup,gpart1,gpart2,gamapdot,btkg1,btkg2,btkg3,bteg1,bteg2,bteg3,tkgf,tegf, &
!$acc tkgm1f,tkgm2f,tkgm3f,tkgm4f,tkgm5f,tegm1f,tegm2f,tegm3f,tegm4f,tegm5f,dtkgx,dtkgy, &
!$acc dtkgz,dtegx,dtegy,dtegz,amutg,aktg,fvtkg,gvtkg,hvtkg,fvteg,gvteg,hvteg,tkgm6f, &
!$acc tkgm7f,tkgm8f,tegm6f,tegm7f,tegm8f,cvg,&
!$acc yfgf,yogf,ypgf,yfgm1f,yfgm2f,yfgm3f,yfgm4f,yfgm5f,yogm1f,yogm2f,yogm3f,yogm4f,yogm5f, &
!$acc ypgm1f,ypgm2f,ypgm3f,ypgm4f,ypgm5f, &
!$acc byfg1,byfg2,byfg3,byog1,byog2,byog3,bypg1,bypg2,bypg3,&
!$acc dyfgx,dyfgy,dyfgz,dyogx,dyogy,dyogz,dypgx,dypgy,dypgz,yfgm6f,yfgm7f,yfgm8f,yogm6f,yogm7f,yogm8f, &
!$acc ypgm6f,ypgm7f,ypgm8f)	
 
  
    do i = 1, neles
if(iblank(i)/=1)cycle !-------(mod3)

    n1 = nod(i,1)
    n2 = nod(i,2)
    n3 = nod(i,3)
    n4 = nod(i,4)
        
!$acc loop seq
        do iface = 1, 4 
      select case (iface)
        case (1)
          scx = sc1x(i)
          scy = sc1y(i)
          scz = sc1z(i)
          xf = (x(n2) + x(n3) + x(n4)) / 3.0d0
          yf = (y(n2) + y(n3) + y(n4)) / 3.0d0
          zf = (z(n2) + z(n3) + z(n4)) / 3.0d0
          neigh = nc1(i)
        case (2)
          scx = sc2x(i)
          scy = sc2y(i)
          scz = sc2z(i)
          xf = (x(n1) + x(n3) + x(n4)) / 3.0d0
          yf = (y(n1) + y(n3) + y(n4)) / 3.0d0
          zf = (z(n1) + z(n3) + z(n4)) / 3.0d0
          neigh = nc2(i)
        case (3)
          scx = sc3x(i)
          scy = sc3y(i)
          scz = sc3z(i)
          xf = (x(n1) + x(n2) + x(n4)) / 3.0d0
          yf = (y(n1) + y(n2) + y(n4)) / 3.0d0
          zf = (z(n1) + z(n2) + z(n4)) / 3.0d0
          neigh = nc3(i)
        case (4)
          scx = sc4x(i)
          scy = sc4y(i)
          scz = sc4z(i)
          xf = (x(n1) + x(n2) + x(n3)) / 3.0d0
          yf = (y(n1) + y(n2) + y(n3)) / 3.0d0
          zf = (z(n1) + z(n2) + z(n3)) / 3.0d0
          neigh = nc4(i)
      end select

#ifdef YP_P
      yfgf = (yfg(i)+yfg(neigh)) * 0.5d0
      yogf = (yog(i)+yog(neigh)) * 0.5d0
      ypgf = (ypg(i)+ypg(neigh)) * 0.5d0
	  
      cpg=yfgf*cpfg+yogf*cpog+ypgf*cppg
      rrg=8314.0d0*(yfgf/amolfg+yogf/amolog+ypgf/amolpg)
      amolg=8314.0d0/rrg
      cvg=cpg-rrg
      gamag=cpg/cvg
      CMUG = 11.848e-8 * DSQRT(amolg)
#endif

!      print*, scx, scy, scz, xf,yf,zf
!if(myid.eq.1.and.total.eq.2.and.neigh.gt.(neles+nghosts))  print*, "neighs", neigh
      rogf = (rog(i) + rog(neigh)) * 0.5d0
      ugf = (ug(i) + ug(neigh)) * 0.5d0
      vgf = (vg(i) + vg(neigh)) * 0.5d0
      wgf = (wg(i) + wg(neigh)) * 0.5d0
      tgf = (tg(i) + tg(neigh)) * 0.5d0
      phigf = (phig(i) + phig(neigh)) * 0.5d0
      amug =  as*dsqrt(tgf)/(1.0d0+ ts/tgf)
      akg = amug * cpg / prg

      ropf = (rop(i) + rop(neigh)) * 0.5d0
      upf = (up(i) + up(neigh)) * 0.5d0
      vpf = (vp(i) + vp(neigh)) * 0.5d0
      wpf = (wp(i) + wp(neigh)) * 0.5d0
      tpf = (tp(i) + tp(neigh)) * 0.5d0
      phipf = (phip(i) + phip(neigh)) * 0.5d0
      amup = 0.001d0
      akp = amup*cpp/prp

#ifdef KE_TURB
      tkgf = dmax1((tkg(i) + tkg(neigh)) * 0.5d0 , tolkg)
      tegf = (teg(i) + teg(neigh)) * 0.5d0
#endif   



 ! Check condition
  !if ( ((neigh > neles) .and. (neigh < neles+nghosts)) .or. (neigh > neles+nghosts) ) then !for bound ghosts and proc ghosts

!  if (neigh > neles.and.neigh<=neles+nghosts) then !Applies for bound ghosts only
if (neigh > neles)  then     ! Applies for ng and proc ghost both
           
  m1 = nc1(i);  m2 = nc2(i); m3 = nc3(i); m4 = nc4(i);  m5 = i
!if(myid.eq.1.and.total.eq.2.and.neigh.gt.(neles+nghosts))  print*,ug(neigh),tg(neigh),phig(neigh),phip(neigh) 

    xcelm1f = xcel(m1)-xf; xcelm2f = xcel(m2)-xf; xcelm3f = xcel(m3)-xf; xcelm4f = xcel(m4)-xf 
    xcelm5f = xcel(m5)-xf; 
    
    ycelm1f = ycel(m1)-yf; ycelm2f = ycel(m2)-yf; ycelm3f = ycel(m3)-yf; ycelm4f = ycel(m4)-yf; ycelm5f = ycel(m5)-yf
    zcelm1f = zcel(m1)-zf; zcelm2f = zcel(m2)-zf; zcelm3f = zcel(m3)-zf; zcelm4f = zcel(m4)-zf; zcelm5f = zcel(m5)-zf
    ugm1f = ug(m1)-ugf; ugm2f= ug(m2)-ugf; ugm3f= ug(m3)-ugf; ugm4f= ug(m4)-ugf;  ugm5f= ug(m5)-ugf
    vgm1f = vg(m1)-vgf; vgm2f= vg(m2)-vgf; vgm3f= vg(m3)-vgf; vgm4f= vg(m4)-vgf;  vgm5f= vg(m5)-vgf
    wgm1f = wg(m1)-wgf; wgm2f= wg(m2)-wgf; wgm3f= wg(m3)-wgf; wgm4f= wg(m4)-wgf;  wgm5f= wg(m5)-wgf
    tgm1f = tg(m1)-tgf; tgm2f = tg(m2)-tgf; tgm3f = tg(m3)-tgf; tgm4f = tg(m4)-tgf; tgm5f = tg(m5)-tgf
    upm1f = up(m1)-upf; upm2f= up(m2)-upf; upm3f= up(m3)-upf; upm4f= up(m4)-upf;  upm5f= up(m5)-upf
    vpm1f = vp(m1)-vpf; vpm2f= vp(m2)-vpf; vpm3f= vp(m3)-vpf; vpm4f= vp(m4)-vpf;  vpm5f= vp(m5)-vpf
    wpm1f = wp(m1)-wpf; wpm2f= wp(m2)-wpf; wpm3f= wp(m3)-wpf; wpm4f= wp(m4)-wpf;  wpm5f= wp(m5)-wpf
    tpm1f = tp(m1)-tpf; tpm2f = tp(m2)-tpf; tpm3f = tp(m3)-tpf; tpm4f = tp(m4)-tpf; tpm5f = tp(m5)-tpf

#ifdef KE_TURB
    tkgm1f = tkg(m1)-tkgf; tkgm2f = tkg(m2)-tkgf; tkgm3f = tkg(m3)-tkgf; tkgm4f = tkg(m4)-tkgf; tkgm5f = tkg(m5)-tkgf
    tegm1f = teg(m1)-tegf; tegm2f = teg(m2)-tegf; tegm3f = teg(m3)-tegf; tegm4f = teg(m4)-tegf; tegm5f = teg(m5)-tegf
#endif

#ifdef YP_P
        yfgm1f = yfg(m1)-yfgf; yfgm2f = yfg(m2)-yfgf; yfgm3f = yfg(m3)-yfgf; yfgm4f = yfg(m4)-yfgf; yfgm5f = yfg(m5)-yfgf
	yogm1f = yog(m1)-yogf; yogm2f = yog(m2)-yogf; yogm3f = yog(m3)-yogf; yogm4f = yog(m4)-yogf; yogm5f = yog(m5)-yogf
	ypgm1f = ypg(m1)-ypgf; ypgm2f = ypg(m2)-ypgf; ypgm3f = ypg(m3)-ypgf; ypgm4f = ypg(m4)-ypgf; ypgm5f = ypg(m5)-ypgf
#endif
    
    
     a11 = xcelm1f**2 + xcelm2f**2 + xcelm3f**2 + xcelm4f**2 + xcelm5f**2
     a22 = ycelm1f**2 + ycelm2f**2 + ycelm3f**2 + ycelm4f**2 + ycelm5f**2
     a33 = zcelm1f**2 + zcelm2f**2 + zcelm3f**2 + zcelm4f**2 + zcelm5f**2
            
     a12 = xcelm1f*(ycelm1f) + xcelm2f*(ycelm2f) + xcelm3f*(ycelm3f) + xcelm4f*(ycelm4f) + (xcelm5f)*(ycelm5f)
     a13 = xcelm1f*(zcelm1f) + xcelm2f*(zcelm2f) + xcelm3f*(zcelm3f) + xcelm4f*(zcelm4f) + (xcelm5f)*(zcelm5f)
     a23 = (ycelm1f)*(zcelm1f) + (ycelm2f)*(zcelm2f) + (ycelm3f)*(zcelm3f) + (ycelm4f)*(zcelm4f) + (ycelm5f)*(zcelm5f)
                                   
     bug1 = (ugm1f)*xcelm1f + (ugm2f)*xcelm2f + (ugm3f)*xcelm3f + (ugm4f)*xcelm4f + (ugm5f)*(xcelm5f) 
     bug2 = (ugm1f)*(ycelm1f) + (ugm2f)*(ycelm2f) + (ugm3f)*(ycelm3f) + (ugm4f)*(ycelm4f) + (ugm5f)*(ycelm5f) 
     bug3 = (ugm1f)*(zcelm1f) + (ugm2f)*(zcelm2f) + (ugm3f)*(zcelm3f) + (ugm4f)*(zcelm4f) + (ugm5f)*(zcelm5f) 
            
     bvg1 = (vgm1f)*xcelm1f + (vgm2f)*xcelm2f + (vgm3f)*xcelm3f + (vgm4f)*xcelm4f + (vgm5f)*(xcelm5f)
     bvg2 = (vgm1f)*(ycelm1f) + (vgm2f)*(ycelm2f) + (vgm3f)*(ycelm3f) + (vgm4f)*(ycelm4f) + (vgm5f)*(ycelm5f)
     bvg3 = (vgm1f)*(zcelm1f) + (vgm2f)*(zcelm2f) + (vgm3f)*(zcelm3f) + (vgm4f)*(zcelm4f) + (vgm5f)*(zcelm5f)
                 
     bwg1 = (wgm1f)*xcelm1f + (wgm2f)*xcelm2f + (wgm3f)*xcelm3f + (wgm4f)*xcelm4f + (wgm5f)*(xcelm5f)
     bwg2 = (wgm1f)*(ycelm1f) + (wgm2f)*(ycelm2f) + (wgm3f)*(ycelm3f) + (wgm4f)*(ycelm4f) + (wgm5f)*(ycelm5f)
     bwg3 = (wgm1f)*(zcelm1f) + (wgm2f)*(zcelm2f) + (wgm3f)*(zcelm3f) + (wgm4f)*(zcelm4f) + (wgm5f)*(zcelm5f)
            
     btg1 = (tgm1f)*xcelm1f + (tgm2f)*xcelm2f + (tgm3f)*xcelm3f + (tgm4f)*xcelm4f + (tgm5f)*(xcelm5f)
     btg2 = (tgm1f)*(ycelm1f) + (tgm2f)*(ycelm2f) + (tgm3f)*(ycelm3f) + (tgm4f)*(ycelm4f) + (tgm5f)*(ycelm5f)
     btg3 = (tgm1f)*(zcelm1f) + (tgm2f)*(zcelm2f) + (tgm3f)*(zcelm3f) + (tgm4f)*(zcelm4f) + (tgm5f)*(zcelm5f)

     bup1 = (upm1f)*xcelm1f + (upm2f)*xcelm2f + (upm3f)*xcelm3f + (upm4f)*xcelm4f + (upm5f)*(xcelm5f) 
     bup2 = (upm1f)*(ycelm1f) + (upm2f)*(ycelm2f) + (upm3f)*(ycelm3f) + (upm4f)*(ycelm4f) + (upm5f)*(ycelm5f) 
     bup3 = (upm1f)*(zcelm1f) + (upm2f)*(zcelm2f) + (upm3f)*(zcelm3f) + (upm4f)*(zcelm4f) + (upm5f)*(zcelm5f) 

     bvp1 = (vpm1f)*xcelm1f + (vpm2f)*xcelm2f + (vpm3f)*xcelm3f + (vpm4f)*xcelm4f + (vpm5f)*(xcelm5f)
     bvp2 = (vpm1f)*(ycelm1f) + (vpm2f)*(ycelm2f) + (vpm3f)*(ycelm3f) + (vpm4f)*(ycelm4f) + (vpm5f)*(ycelm5f)
     bvp3 = (vpm1f)*(zcelm1f) + (vpm2f)*(zcelm2f) + (vpm3f)*(zcelm3f) + (vpm4f)*(zcelm4f) + (vpm5f)*(zcelm5f)
     
     bwp1 = (wpm1f)*xcelm1f + (wpm2f)*xcelm2f + (wpm3f)*xcelm3f + (wpm4f)*xcelm4f + (wpm5f)*(xcelm5f)
     bwp2 = (wpm1f)*(ycelm1f) + (wpm2f)*(ycelm2f) + (wpm3f)*(ycelm3f) + (wpm4f)*(ycelm4f) + (wpm5f)*(ycelm5f)
     bwp3 = (wpm1f)*(zcelm1f) + (wpm2f)*(zcelm2f) + (wpm3f)*(zcelm3f) + (wpm4f)*(zcelm4f) + (wpm5f)*(zcelm5f)
            
     btp1 = (tpm1f)*xcelm1f + (tpm2f)*xcelm2f + (tpm3f)*xcelm3f + (tpm4f)*xcelm4f + (tpm5f)*(xcelm5f)
     btp2 = (tpm1f)*(ycelm1f) + (tpm2f)*(ycelm2f) + (tpm3f)*(ycelm3f) + (tpm4f)*(ycelm4f) + (tpm5f)*(ycelm5f)
     btp3 = (tpm1f)*(zcelm1f) + (tpm2f)*(zcelm2f) + (tpm3f)*(zcelm3f) + (tpm4f)*(zcelm4f) + (tpm5f)*(zcelm5f)

 
#ifdef KE_TURB
     btkg1 = (tkgm1f)*(xcelm1f) + (tkgm2f)*(xcelm2f) + (tkgm3f)*(xcelm3f) + (tkgm4f)*(xcelm4f) + (tkgm5f)*(xcelm5f)
     btkg2 = (tkgm1f)*(ycelm1f) + (tkgm2f)*(ycelm2f) + (tkgm3f)*(ycelm3f) + (tkgm4f)*(ycelm4f) + (tkgm5f)*(ycelm5f)
     btkg3 = (tkgm1f)*(zcelm1f) + (tkgm2f)*(zcelm2f) + (tkgm3f)*(zcelm3f) + (tkgm4f)*(zcelm4f) + (tkgm5f)*(zcelm5f)
	 
     bteg1 = (tegm1f)*(xcelm1f) + (tegm2f)*(xcelm2f) + (tegm3f)*(xcelm3f) + (tegm4f)*(xcelm4f) + (tegm5f)*(xcelm5f)
     bteg2 = (tegm1f)*(ycelm1f) + (tegm2f)*(ycelm2f) + (tegm3f)*(ycelm3f) + (tegm4f)*(ycelm4f) + (tegm5f)*(ycelm5f)
     bteg3 = (tegm1f)*(zcelm1f) + (tegm2f)*(zcelm2f) + (tegm3f)*(zcelm3f) + (tegm4f)*(zcelm4f) + (tegm5f)*(zcelm5f)
#endif    

#ifdef YP_P
     byfg1 = (yfgm1f)*(xcelm1f) + (yfgm2f)*(xcelm2f) + (yfgm3f)*(xcelm3f) + (yfgm4f)*(xcelm4f) + (yfgm5f)*(xcelm5f)
     byfg2 = (yfgm1f)*(ycelm1f) + (yfgm2f)*(ycelm2f) + (yfgm3f)*(ycelm3f) + (yfgm4f)*(ycelm4f) + (yfgm5f)*(ycelm5f)
     byfg3 = (yfgm1f)*(zcelm1f) + (yfgm2f)*(zcelm2f) + (yfgm3f)*(zcelm3f) + (yfgm4f)*(zcelm4f) + (yfgm5f)*(zcelm5f)
	 
     byog1 = (yogm1f)*(xcelm1f) + (yogm2f)*(xcelm2f) + (yogm3f)*(xcelm3f) + (yogm4f)*(xcelm4f) + (yogm5f)*(xcelm5f)
     byog2 = (yogm1f)*(ycelm1f) + (yogm2f)*(ycelm2f) + (yogm3f)*(ycelm3f) + (yogm4f)*(ycelm4f) + (yogm5f)*(ycelm5f)
     byog3 = (yogm1f)*(zcelm1f) + (yogm2f)*(zcelm2f) + (yogm3f)*(zcelm3f) + (yogm4f)*(zcelm4f) + (yogm5f)*(zcelm5f)
	 
     bypg1 = (ypgm1f)*(xcelm1f) + (ypgm2f)*(xcelm2f) + (ypgm3f)*(xcelm3f) + (ypgm4f)*(xcelm4f) + (ypgm5f)*(xcelm5f)
     bypg2 = (ypgm1f)*(ycelm1f) + (ypgm2f)*(ycelm2f) + (ypgm3f)*(ycelm3f) + (ypgm4f)*(ycelm4f) + (ypgm5f)*(ycelm5f)
     bypg3 = (ypgm1f)*(zcelm1f) + (ypgm2f)*(zcelm2f) + (ypgm3f)*(zcelm3f) + (ypgm4f)*(zcelm4f) + (ypgm5f)*(zcelm5f)
#endif

     ! Compute determinant of the matrix
  dnr = det(a11, a12, a13, a12, a22, a23, a13, a23, a33)

  dugx = det(bug1, bug2, bug3, a12, a22, a23, a13, a23, a33) / dnr
  dugy = det(a11, a12, a13, bug1, bug2, bug3, a13, a23, a33) / dnr
  dugz = det(a11, a12, a13, a12, a22, a23, bug1, bug2, bug3) / dnr

  dvgx = det(bvg1, bvg2, bvg3, a12, a22, a23, a13, a23, a33) / dnr
  dvgy = det(a11, a12, a13, bvg1, bvg2, bvg3, a13, a23, a33) / dnr
  dvgz = det(a11, a12, a13, a12, a22, a23, bvg1, bvg2, bvg3) / dnr

  dwgx = det(bwg1, bwg2, bwg3, a12, a22, a23, a13, a23, a33) / dnr
  dwgy = det(a11, a12, a13, bwg1, bwg2, bwg3, a13, a23, a33) / dnr
  dwgz = det(a11, a12, a13, a12, a22, a23, bwg1, bwg2, bwg3) / dnr

  dtgx = det(btg1, btg2, btg3, a12, a22, a23, a13, a23, a33) / dnr
  dtgy = det(a11, a12, a13, btg1, btg2, btg3, a13, a23, a33) / dnr
  dtgz = det(a11, a12, a13, a12, a22, a23, btg1, btg2, btg3) / dnr

  dupx = det(bup1, bup2, bup3, a12, a22, a23, a13, a23, a33) / dnr
  dupy = det(a11, a12, a13, bup1, bup2, bup3, a13, a23, a33) / dnr
  dupz = det(a11, a12, a13, a12, a22, a23, bup1, bup2, bup3) / dnr

  dvpx = det(bvp1, bvp2, bvp3, a12, a22, a23, a13, a23, a33) / dnr
  dvpy = det(a11, a12, a13, bvp1, bvp2, bvp3, a13, a23, a33) / dnr
  dvpz = det(a11, a12, a13, a12, a22, a23, bvp1, bvp2, bvp3) / dnr

  dwpx = det(bwp1, bwp2, bwp3, a12, a22, a23, a13, a23, a33) / dnr
  dwpy = det(a11, a12, a13, bwp1, bwp2, bwp3, a13, a23, a33) / dnr
  dwpz = det(a11, a12, a13, a12, a22, a23, bwp1, bwp2, bwp3) / dnr

  dtpx = det(btp1, btp2, btp3, a12, a22, a23, a13, a23, a33) / dnr
  dtpy = det(a11, a12, a13, btp1, btp2, btp3, a13, a23, a33) / dnr
  dtpz = det(a11, a12, a13, a12, a22, a23, btp1, btp2, btp3) / dnr

#ifdef KE_TURB
  dtkgx = det(btkg1, btkg2, btkg3, a12, a22, a23, a13, a23, a33) / dnr
  dtkgy = det(a11, a12, a13, btkg1, btkg2, btkg3, a13, a23, a33) / dnr
  dtkgz = det(a11, a12, a13, a12, a22, a23, btkg1, btkg2, btkg3) / dnr
  
  dtegx = det(bteg1, bteg2, bteg3, a12, a22, a23, a13, a23, a33) / dnr
  dtegy = det(a11, a12, a13, bteg1, bteg2, bteg3, a13, a23, a33) / dnr
  dtegz = det(a11, a12, a13, a12, a22, a23, bteg1, bteg2, bteg3) / dnr
#endif

#ifdef YP_P
  dyfgx = det(byfg1, byfg2, byfg3, a12, a22, a23, a13, a23, a33) / dnr
  dyfgy = det(a11, a12, a13, byfg1, byfg2, byfg3, a13, a23, a33) / dnr
  dyfgz = det(a11, a12, a13, a12, a22, a23, byfg1, byfg2, byfg3) / dnr
  
  dyogx = det(byog1, byog2, byog3, a12, a22, a23, a13, a23, a33) / dnr
  dyogy = det(a11, a12, a13, byog1, byog2, byog3, a13, a23, a33) / dnr
  dyogz = det(a11, a12, a13, a12, a22, a23, byog1, byog2, byog3) / dnr
  
  dypgx = det(bypg1, bypg2, bypg3, a12, a22, a23, a13, a23, a33) / dnr
  dypgy = det(a11, a12, a13, bypg1, bypg2, bypg3, a13, a23, a33) / dnr
  dypgz = det(a11, a12, a13, a12, a22, a23, bypg1, bypg2, bypg3) / dnr
#endif

!--------------------------------------------------------------------------------------------------------------
 
#ifndef AMUP_CONSTANT 
        if(phipf.gt.0.01d0)then
        gpart1=2.0d0*(dupx**2+dvpy**2+dwpz**2)
	gpart2=(dupy+dvpx)**2+(dupz+dwpx)**2+(dvpz+dwpy)**2
        gamapdot=min(max(sqrt(gpart1+gpart2),0.0001d0),10.0d0)	 
         
           if(gamapdot.lt.0.01d0)then
	     amup=2658.83d0*((2.0d0-0.814d0)+(0.814d0-1.0d0)*(gamapdot/0.01d0))+10.7d0*(2.0d0-(gamapdot/0.01d0))/0.01d0 
	   else
	     amup=10.7d0/gamapdot+2658.83d0*(gamapdot/0.01d0)**(0.814d0-1.0d0)
	   endif

	else
	   amup=0.001d0 !amug
	endif
		
	 akp=amup*cpp/prp
#endif

!--------------------------------------------------------------------------------------------------------------	 
#ifdef KE_TURB
       amutg=0.0845d0*rogf*tkgf*tkgf/tegf
       aktg=amutg*cpg/0.9d0
       amug=amug+amutg
       akg=akg+aktg
#endif

  delvg = dugx + dvgy + dwgz
  
  txxg = amug*(2.0d0*dugx - c23*delvg)
  tyyg = amug*(2.0d0*dvgy - c23*delvg)
  tzzg = amug*(2.0d0*dwgz - c23*delvg)
  txyg = amug*(dvgx + dugy)
  txzg = amug*(dwgx + dugz)
  tyzg = amug*(dvgz + dwgy)

  qxxg = ugf*txxg + vgf*txyg + wgf*txzg + akg*dtgx
  qyyg = ugf*txyg + vgf*tyyg + wgf*tyzg + akg*dtgy
  qzzg = ugf*txzg + vgf*tyzg + wgf*tzzg + akg*dtgz

  delvp = dupx + dvpy + dwpz
  txxp = amup*(2.0d0*dupx - c23*delvp)
  tyyp = amup*(2.0d0*dvpy - c23*delvp)
  tzzp = amup*(2.0d0*dwpz - c23*delvp)
  txyp = amup*(dvpx + dupy)
  txzp = amup*(dwpx + dupz)
  tyzp = amup*(dvpz + dwpy)

  qxxp = upf*txxp + vpf*txyp + wpf*txzp + akp*dtpx
  qyyp = upf*txyp + vpf*tyyp + wpf*tyzp + akp*dtpy
  qzzp = upf*txzp + vpf*tyzp + wpf*tzzp + akp*dtpz
  
#ifdef KE_TURB
    fvtkg=1.39d0*amug*dtkgx
    gvtkg=1.39d0*amug*dtkgy
    hvtkg=1.39d0*amug*dtkgz

    fvteg=1.39d0*amug*dtegx
    gvteg=1.39d0*amug*dtegy
    hvteg=1.39d0*amug*dtegz
#endif


  visug(i) = visug(i) + phigf*(txxg*scx + txyg*scy + txzg*scz)
  visvg(i) = visvg(i) + phigf*(txyg*scx + tyyg*scy + tyzg*scz)
  viswg(i) = viswg(i) + phigf*(txzg*scx + tyzg*scy + tzzg*scz)
  vistg(i) = vistg(i) + phigf*(qxxg*scx + qyyg*scy + qzzg*scz)
  visup(i) = visup(i) + phipf*(txxp*scx + txyp*scy + txzp*scz)
  visvp(i) = visvp(i) + phipf*(txyp*scx + tyyp*scy + tyzp*scz)
  viswp(i) = viswp(i) + phipf*(txzp*scx + tyzp*scy + tzzp*scz)
  vistp(i) = vistp(i) + phipf*(qxxp*scx + qyyp*scy + qzzp*scz)

#ifdef KE_TURB
  vistkg(i) = vistkg(i) + phigf*(fvtkg*scx + gvtkg*scy + hvtkg*scz)
  visteg(i) = visteg(i) + phigf*(fvteg*scx + gvteg*scy + hvteg*scz)
#endif


#ifdef YP_P
  visyfg(i) = visyfg(i) + phigf*amug*(dyfgx*scx + dyfgy*scy + dyfgz*scz)
  visyog(i) = visyog(i) + phigf*amug*(dyogx*scx + dyogy*scy + dyogz*scz)
  visypg(i) = visypg(i) + phigf*amug*(dypgx*scx + dypgy*scy + dypgz*scz)
#endif
!do i = 1, neles
!   if (myid.eq.1.and.total == 2) then
!      write(*,'(A,1X,I8,2(1X,F25.15))') "vis", i, dnr, dvgx
!   end if
!end do

else
     ! If a processor ghost cell appears here, the code will attempt to access 
     ! cell information that does not exist, leading to a memory error.
     ! Assign index values
    m1 = nc1(i);  m2 = nc2(i); m3 = nc3(i); m4 = nc4(i)
    m5 = nc1(neigh); m6 = nc2(neigh);  m7 = nc3(neigh);  m8 = nc4(neigh)
    
    xcelm1f = xcel(m1)-xf; xcelm2f = xcel(m2)-xf; xcelm3f = xcel(m3)-xf; xcelm4f = xcel(m4)-xf 
    xcelm5f = xcel(m5)-xf; xcelm6f = xcel(m6)-xf; xcelm7f = xcel(m7)-xf; xcelm8f = xcel(m8)-xf
    
    ycelm1f = ycel(m1)-yf; ycelm2f = ycel(m2)-yf; ycelm3f = ycel(m3)-yf; ycelm4f = ycel(m4)-yf 
    ycelm5f = ycel(m5)-yf; ycelm6f = ycel(m6)-yf; ycelm7f = ycel(m7)-yf; ycelm8f = ycel(m8)-yf 
    
    zcelm1f = zcel(m1)-zf; zcelm2f = zcel(m2)-zf; zcelm3f = zcel(m3)-zf; zcelm4f = zcel(m4)-zf
    zcelm5f = zcel(m5)-zf; zcelm6f = zcel(m6)-zf; zcelm7f = zcel(m7)-zf; zcelm8f = zcel(m8)-zf
    
    ugm1f = ug(m1)-ugf; ugm2f= ug(m2)-ugf; ugm3f= ug(m3)-ugf; ugm4f= ug(m4)-ugf
    ugm5f= ug(m5)-ugf; ugm6f= ug(m6)-ugf; ugm7f= ug(m7)-ugf; ugm8f= ug(m8)-ugf

    vgm1f = vg(m1)-vgf; vgm2f= vg(m2)-vgf; vgm3f= vg(m3)-vgf; vgm4f= vg(m4)-vgf
    vgm5f= vg(m5)-vgf;  vgm6f= vg(m6)-vgf;  vgm7f= vg(m7)-vgf;  vgm8f= vg(m8)-vgf 
    
    wgm1f = wg(m1)-wgf; wgm2f= wg(m2)-wgf; wgm3f= wg(m3)-wgf; wgm4f= wg(m4)-wgf
    wgm5f= wg(m5)-wgf; wgm6f= wg(m6)-wgf; wgm7f= wg(m7)-wgf; wgm8f= wg(m8)-wgf
    
    tgm1f = tg(m1)-tgf; tgm2f = tg(m2)-tgf; tgm3f = tg(m3)-tgf; tgm4f = tg(m4)-tgf
    tgm5f = tg(m5)-tgf; tgm6f = tg(m6)-tgf; tgm7f = tg(m7)-tgf; tgm8f = tg(m8)-tgf
    
    upm1f = up(m1)-upf; upm2f= up(m2)-upf; upm3f= up(m3)-upf; upm4f= up(m4)-upf
    upm5f= up(m5)-upf; upm6f= up(m6)-upf; upm7f= up(m7)-upf; upm8f= up(m8)-upf
     
    vpm1f = vp(m1)-vpf; vpm2f= vp(m2)-vpf; vpm3f= vp(m3)-vpf; vpm4f= vp(m4)-vpf
    vpm5f= vp(m5)-vpf; vpm6f= vp(m6)-vpf; vpm7f= vp(m7)-vpf; vpm8f= vp(m8)-vpf
    
    wpm1f = wp(m1)-wpf; wpm2f= wp(m2)-wpf; wpm3f= wp(m3)-wpf; wpm4f= wp(m4)-wpf
    wpm5f= wp(m5)-wpf; wpm6f= wp(m6)-wpf; wpm7f= wp(m7)-wpf; wpm8f= wp(m8)-wpf
     
    tpm1f = tp(m1)-tpf; tpm2f = tp(m2)-tpf; tpm3f = tp(m3)-tpf; tpm4f = tp(m4)-tpf 
    tpm5f = tp(m5)-tpf; tpm6f = tp(m6)-tpf; tpm7f = tp(m7)-tpf; tpm8f = tp(m8)-tpf 
   
#ifdef KE_TURB
    tkgm1f = tkg(m1)-tkgf; tkgm2f = tkg(m2)-tkgf; tkgm3f = tkg(m3)-tkgf; tkgm4f = tkg(m4)-tkgf
    tkgm5f = tkg(m5)-tkgf; tkgm6f = tkg(m6)-tkgf; tkgm7f = tkg(m7)-tkgf; tkgm8f = tkg(m8)-tkgf
	
    tegm1f = teg(m1)-tegf; tegm2f = teg(m2)-tegf; tegm3f = teg(m3)-tegf; tegm4f = teg(m4)-tegf
    tegm5f = teg(m5)-tegf; tegm6f = teg(m6)-tegf; tegm7f = teg(m7)-tegf; tegm8f = teg(m8)-tegf
#endif   

#ifdef YP_P
    yfgm1f = yfg(m1)-yfgf; yfgm2f = yfg(m2)-yfgf; yfgm3f = yfg(m3)-yfgf; yfgm4f = yfg(m4)-yfgf
    yfgm5f = yfg(m5)-yfgf; yfgm6f = yfg(m6)-yfgf; yfgm7f = yfg(m7)-yfgf; yfgm8f = yfg(m8)-yfgf
	
    yogm1f = yog(m1)-yogf; yogm2f = yog(m2)-yogf; yogm3f = yog(m3)-yogf; yogm4f = yog(m4)-yogf
    yogm5f = yog(m5)-yogf; yogm6f = yog(m6)-yogf; yogm7f = yog(m7)-yogf; yogm8f = yog(m8)-yogf
	
     ypgm1f = ypg(m1)-ypgf; ypgm2f = ypg(m2)-ypgf; ypgm3f = ypg(m3)-ypgf; ypgm4f = ypg(m4)-ypgf
    ypgm5f = ypg(m5)-ypgf; ypgm6f = ypg(m6)-ypgf; ypgm7f = ypg(m7)-ypgf; ypgm8f = ypg(m8)-ypgf
#endif  
    
     a11 = xcelm1f**2 + xcelm2f**2 + xcelm3f**2 + xcelm4f**2 + (xcelm5f)**2 + (xcelm6f)**2 + xcelm7f**2 + xcelm8f**2
     a22 = (ycelm1f)**2 + (ycelm2f)**2 + (ycelm3f)**2 + (ycelm4f)**2 + (ycelm5f)**2 + (ycelm6f)**2 + (ycelm7f)**2 + (ycelm8f)**2
     a33 = (zcelm1f)**2 + (zcelm2f)**2 + (zcelm3f)**2 + (zcelm4f)**2 + (zcelm5f)**2 + (zcelm6f)**2 + (zcelm7f)**2 + (zcelm8f)**2
            
     a12 = xcelm1f*(ycelm1f) + xcelm2f*(ycelm2f) + xcelm3f*(ycelm3f) + xcelm4f* (ycelm4f) + (xcelm5f)*(ycelm5f) + &
           (xcelm6f)*(ycelm6f) + xcelm7f* (ycelm7f) + xcelm8f*(ycelm8f) 
            
     a13 = xcelm1f*(zcelm1f) + xcelm2f*(zcelm2f) + xcelm3f*(zcelm3f) + xcelm4f*(zcelm4f) + (xcelm5f)*(zcelm5f) + &
           (xcelm6f)*(zcelm6f) + xcelm7f*(zcelm7f) + xcelm8f*(zcelm8f) 

     a23 = (ycelm1f)*(zcelm1f) + (ycelm2f)*(zcelm2f) + (ycelm3f)*(zcelm3f) + (ycelm4f)*(zcelm4f) + (ycelm5f)*(zcelm5f) + &
           (ycelm6f)*(zcelm6f) + (ycelm7f)*(zcelm7f) + (ycelm8f)*(zcelm8f) 
            
     bug1 = (ugm1f)*xcelm1f + (ugm2f)*xcelm2f + (ugm3f)*xcelm3f + (ugm4f)*xcelm4f + (ugm5f)*(xcelm5f) + (ugm6f)*(xcelm6f) + &
            (ugm7f)*xcelm7f + (ugm8f)*xcelm8f  
     bug2 = (ugm1f)*(ycelm1f) + (ugm2f)*(ycelm2f) + (ugm3f)*(ycelm3f) + (ugm4f)*(ycelm4f) + (ugm5f)*(ycelm5f) + (ugm6f)*(ycelm6f)+&
            (ugm7f)*(ycelm7f) + (ugm8f)*(ycelm8f)
     bug3 = (ugm1f)*(zcelm1f) + (ugm2f)*(zcelm2f) + (ugm3f)*(zcelm3f) + (ugm4f)*(zcelm4f) + (ugm5f)*(zcelm5f) + (ugm6f)*(zcelm6f)+&
            (ugm7f)*(zcelm7f) + (ugm8f)*(zcelm8f)

     bvg1 = (vgm1f)*xcelm1f + (vgm2f)*xcelm2f + (vgm3f)*xcelm3f + (vgm4f)*xcelm4f + (vgm5f)*(xcelm5f) + (vgm6f)*(xcelm6f) + &
            (vgm7f)*xcelm7f + (vgm8f)*xcelm8f  
     bvg2 = (vgm1f)*(ycelm1f) + (vgm2f)*(ycelm2f) + (vgm3f)*(ycelm3f) + (vgm4f)*(ycelm4f) + (vgm5f)*(ycelm5f) + (vgm6f)*(ycelm6f)+&
            (vgm7f)*(ycelm7f) + (vgm8f)*(ycelm8f)
     bvg3 = (vgm1f)*(zcelm1f) + (vgm2f)*(zcelm2f) + (vgm3f)*(zcelm3f) + (vgm4f)*(zcelm4f) + (vgm5f)*(zcelm5f) + (vgm6f)*(zcelm6f)+&
            (vgm7f)*(zcelm7f) + (vgm8f)*(zcelm8f)

     bwg1 = (wgm1f)*xcelm1f + (wgm2f)*xcelm2f + (wgm3f)*xcelm3f + (wgm4f)*xcelm4f + (wgm5f)*(xcelm5f) + (wgm6f)*(xcelm6f) + &
            (wgm7f)*xcelm7f + (wgm8f)*xcelm8f  
     bwg2 = (wgm1f)*(ycelm1f) + (wgm2f)*(ycelm2f) + (wgm3f)*(ycelm3f) + (wgm4f)*(ycelm4f) + (wgm5f)*(ycelm5f) + (wgm6f)*(ycelm6f)+&
            (wgm7f)*(ycelm7f) + (wgm8f)*(ycelm8f)
     bwg3 = (wgm1f)*(zcelm1f) + (wgm2f)*(zcelm2f) + (wgm3f)*(zcelm3f) + (wgm4f)*(zcelm4f) + (wgm5f)*(zcelm5f) + (wgm6f)*(zcelm6f)+&
            (wgm7f)*(zcelm7f) + (wgm8f)*(zcelm8f)
            
     btg1 = (tgm1f)*xcelm1f + (tgm2f)*xcelm2f + (tgm3f)*xcelm3f + (tgm4f)*xcelm4f + (tgm5f)*(xcelm5f) + (tgm6f)*(xcelm6f) + &
            (tgm7f)*xcelm7f + (tgm8f)*xcelm8f
     btg2 = (tgm1f)*(ycelm1f) + (tgm2f)*(ycelm2f) + (tgm3f)*(ycelm3f) + (tgm4f)*(ycelm4f) + (tgm5f)*(ycelm5f) + (tgm6f)*(ycelm6f)+&
            (tgm7f)*(ycelm7f) + (tgm8f)*(ycelm8f)
     btg3 = (tgm1f)*(zcelm1f) + (tgm2f)*(zcelm2f) + (tgm3f)*(zcelm3f) + (tgm4f)*(zcelm4f) + (tgm5f)*(zcelm5f) + (tgm6f)*(zcelm6f)+&
            (tgm7f)*(zcelm7f) + (tgm8f)*(zcelm8f)

     bup1 = (upm1f)*xcelm1f + (upm2f)*xcelm2f + (upm3f)*xcelm3f + (upm4f)*xcelm4f + (upm5f)*(xcelm5f) + (upm6f)*(xcelm6f) + &
            (upm7f)*xcelm7f + (upm8f)*xcelm8f  
     bup2 = (upm1f)*(ycelm1f) + (upm2f)*(ycelm2f) + (upm3f)*(ycelm3f) + (upm4f)*(ycelm4f) + (upm5f)*(ycelm5f) + (upm6f)*(ycelm6f)+&
            (upm7f)*(ycelm7f) + (upm8f)*(ycelm8f)
     bup3 = (upm1f)*(zcelm1f) + (upm2f)*(zcelm2f) + (upm3f)*(zcelm3f) + (upm4f)*(zcelm4f) + (upm5f)*(zcelm5f) + (upm6f)*(zcelm6f)+&
            (upm7f)*(zcelm7f) + (upm8f)*(zcelm8f)

     bvp1 = (vpm1f)*xcelm1f + (vpm2f)*xcelm2f + (vpm3f)*xcelm3f + (vpm4f)*xcelm4f + (vpm5f)*(xcelm5f) + (vpm6f)*(xcelm6f) + &
            (vpm7f)*xcelm7f + (vpm8f)*xcelm8f  
     bvp2 = (vpm1f)*(ycelm1f) + (vpm2f)*(ycelm2f) + (vpm3f)*(ycelm3f) + (vpm4f)*(ycelm4f) + (vpm5f)*(ycelm5f) + (vpm6f)*(ycelm6f)+&
            (vpm7f)*(ycelm7f) + (vpm8f)*(ycelm8f)
     bvp3 = (vpm1f)*(zcelm1f) + (vpm2f)*(zcelm2f) + (vpm3f)*(zcelm3f) + (vpm4f)*(zcelm4f) + (vpm5f)*(zcelm5f) + (vpm6f)*(zcelm6f)+&
            (vpm7f)*(zcelm7f) + (vpm8f)*(zcelm8f)

     bwp1 = (wpm1f)*xcelm1f + (wpm2f)*xcelm2f + (wpm3f)*xcelm3f + (wpm4f)*xcelm4f + (wpm5f)*(xcelm5f) + (wpm6f)*(xcelm6f) + &
            (wpm7f)*xcelm7f + (wpm8f)*xcelm8f  
     bwp2 = (wpm1f)*(ycelm1f) + (wpm2f)*(ycelm2f) + (wpm3f)*(ycelm3f) + (wpm4f)*(ycelm4f) + (wpm5f)*(ycelm5f) + (wpm6f)*(ycelm6f)+&
            (wpm7f)*(ycelm7f) + (wpm8f)*(ycelm8f)
     bwp3 = (wpm1f)*(zcelm1f) + (wpm2f)*(zcelm2f) + (wpm3f)*(zcelm3f) + (wpm4f)*(zcelm4f) + (wpm5f)*(zcelm5f) + (wpm6f)*(zcelm6f)+&
            (wpm7f)*(zcelm7f) + (wpm8f)*(zcelm8f)
            
     btp1 = tpm1f*xcelm1f + tpm2f*xcelm2f + tpm3f*xcelm3f + tpm4f*xcelm4f + tpm5f*xcelm5f + tpm6f*xcelm6f + tpm7f*xcelm7f + &
            tpm8f*xcelm8f
     btp2 = tpm1f*ycelm1f + tpm2f*ycelm2f + tpm3f*ycelm3f + tpm4f*ycelm4f + tpm5f*ycelm5f + tpm6f*ycelm6f + tpm7f*ycelm7f + &
            tpm8f*ycelm8f
     btp3 = tpm1f*zcelm1f + tpm2f*zcelm2f + tpm3f*zcelm3f + tpm4f*zcelm4f + tpm5f*zcelm5f + tpm6f*zcelm6f + tpm7f*zcelm7f + &
            tpm8f*zcelm8f
            
#ifdef KE_TURB
     btkg1 = (tkgm1f)*xcelm1f + (tkgm2f)*xcelm2f + (tkgm3f)*xcelm3f + (tkgm4f)*xcelm4f + (tkgm5f)*(xcelm5f) + &
             (tkgm6f)*(xcelm6f) + (tkgm7f)*xcelm7f + (tkgm8f)*xcelm8f
     btkg2 = (tkgm1f)*(ycelm1f) + (tkgm2f)*(ycelm2f) + (tkgm3f)*(ycelm3f) + (tkgm4f)*(ycelm4f) + (tkgm5f)*(ycelm5f) + &
             (tkgm6f)*(ycelm6f)+  (tkgm7f)*(ycelm7f) + (tkgm8f)*(ycelm8f)
     btkg3 = (tkgm1f)*(zcelm1f) + (tkgm2f)*(zcelm2f) + (tkgm3f)*(zcelm3f) + (tkgm4f)*(zcelm4f) + (tkgm5f)*(zcelm5f) + &
             (tkgm6f)*(zcelm6f)+ (tkgm7f)*(zcelm7f) + (tkgm8f)*(zcelm8f)
			
     bteg1 = (tegm1f)*xcelm1f + (tegm2f)*xcelm2f + (tegm3f)*xcelm3f + (tegm4f)*xcelm4f + (tegm5f)*(xcelm5f) + &
             (tegm6f)*(xcelm6f) + (tegm7f)*xcelm7f + (tegm8f)*xcelm8f
     bteg2 = (tegm1f)*(ycelm1f) + (tegm2f)*(ycelm2f) + (tegm3f)*(ycelm3f) + (tegm4f)*(ycelm4f) + (tegm5f)*(ycelm5f) + &
             (tegm6f)*(ycelm6f) + (tegm7f)*(ycelm7f) + (tegm8f)*(ycelm8f)
     bteg3 = (tegm1f)*(zcelm1f) + (tegm2f)*(zcelm2f) + (tegm3f)*(zcelm3f) + (tegm4f)*(zcelm4f) + (tegm5f)*(zcelm5f) + &
             (tegm6f)*(zcelm6f) +  (tegm7f)*(zcelm7f) + (tegm8f)*(zcelm8f)		
#endif   

#ifdef YP_P
     byfg1 = (yfgm1f)*(xcelm1f) + (yfgm2f)*(xcelm2f) + (yfgm3f)*(xcelm3f) + (yfgm4f)*(xcelm4f) + (yfgm5f)*(xcelm5f) + &
             (yfgm6f)*(xcelm6f) + (yfgm7f)*(xcelm7f) + (yfgm8f)*(xcelm8f)
     byfg2 = (yfgm1f)*(ycelm1f) + (yfgm2f)*(ycelm2f) + (yfgm3f)*(ycelm3f) + (yfgm4f)*(ycelm4f) + (yfgm5f)*(ycelm5f) + &
             (yfgm6f)*(ycelm6f) + (yfgm7f)*(ycelm7f) + (yfgm8f)*(ycelm8f)
     byfg3 = (yfgm1f)*(zcelm1f) + (yfgm2f)*(zcelm2f) + (yfgm3f)*(zcelm3f) + (yfgm4f)*(zcelm4f) + (yfgm5f)*(zcelm5f) + &
             (yfgm6f)*(zcelm6f) + (yfgm7f)*(zcelm7f) + (yfgm8f)*(zcelm8f)
			 
    byog1 = (yogm1f)*(xcelm1f) + (yogm2f)*(xcelm2f) + (yogm3f)*(xcelm3f) + (yogm4f)*(xcelm4f) + (yogm5f)*(xcelm5f) + &
             (yogm6f)*(xcelm6f) + (yogm7f)*(xcelm7f) + (yogm8f)*(xcelm8f)
     byog2 = (yogm1f)*(ycelm1f) + (yogm2f)*(ycelm2f) + (yogm3f)*(ycelm3f) + (yogm4f)*(ycelm4f) + (yogm5f)*(ycelm5f) + &
             (yogm6f)*(ycelm6f) + (yogm7f)*(ycelm7f) + (yogm8f)*(ycelm8f)
     byog3 = (yogm1f)*(zcelm1f) + (yogm2f)*(zcelm2f) + (yogm3f)*(zcelm3f) + (yogm4f)*(zcelm4f) + (yogm5f)*(zcelm5f) + &
             (yogm6f)*(zcelm6f) + (yogm7f)*(zcelm7f) + (yogm8f)*(zcelm8f)

     bypg1 = (ypgm1f)*(xcelm1f) + (ypgm2f)*(xcelm2f) + (ypgm3f)*(xcelm3f) + (ypgm4f)*(xcelm4f) + (ypgm5f)*(xcelm5f) + &
             (ypgm6f)*(xcelm6f) + (ypgm7f)*(xcelm7f) + (ypgm8f)*(xcelm8f)
     bypg2 = (ypgm1f)*(ycelm1f) + (ypgm2f)*(ycelm2f) + (ypgm3f)*(ycelm3f) + (ypgm4f)*(ycelm4f) + (ypgm5f)*(ycelm5f) + &
             (ypgm6f)*(ycelm6f) + (ypgm7f)*(ycelm7f) + (ypgm8f)*(ycelm8f)
     bypg3 = (ypgm1f)*(zcelm1f) + (ypgm2f)*(zcelm2f) + (ypgm3f)*(zcelm3f) + (ypgm4f)*(zcelm4f) + (ypgm5f)*(zcelm5f) + &
             (ypgm6f)*(zcelm6f) + (ypgm7f)*(zcelm7f) + (ypgm8f)*(zcelm8f)			 
#endif          

    ! Compute values
    dnr=det(a11,a12,a13,a12,a22,a23,a13,a23,a33)
    
    ! Calculate dugx, dugy, dugz, etc., using the determinant function
    dugx = det(bug1, bug2, bug3, a12, a22, a23, a13, a23, a33) / dnr
    dugy = det(a11, a12, a13, bug1, bug2, bug3, a13, a23, a33) / dnr
    dugz = det(a11, a12, a13, a12, a22, a23, bug1, bug2, bug3) / dnr

    dvgx = det(bvg1, bvg2, bvg3, a12, a22, a23, a13, a23, a33) / dnr
    dvgy = det(a11, a12, a13, bvg1, bvg2, bvg3, a13, a23, a33) / dnr
    dvgz = det(a11, a12, a13, a12, a22, a23, bvg1, bvg2, bvg3) / dnr

    dwgx = det(bwg1, bwg2, bwg3, a12, a22, a23, a13, a23, a33) / dnr
    dwgy = det(a11, a12, a13, bwg1, bwg2, bwg3, a13, a23, a33) / dnr
    dwgz = det(a11, a12, a13, a12, a22, a23, bwg1, bwg2, bwg3) / dnr

    dtgx = det(btg1, btg2, btg3, a12, a22, a23, a13, a23, a33) / dnr
    dtgy = det(a11, a12, a13, btg1, btg2, btg3, a13, a23, a33) / dnr
    dtgz = det(a11, a12, a13, a12, a22, a23, btg1, btg2, btg3) / dnr

    dupx = det(bup1, bup2, bup3, a12, a22, a23, a13, a23, a33) / dnr
    dupy = det(a11, a12, a13, bup1, bup2, bup3, a13, a23, a33) / dnr
    dupz = det(a11, a12, a13, a12, a22, a23, bup1, bup2, bup3) / dnr

    dvpx = det(bvp1, bvp2, bvp3, a12, a22, a23, a13, a23, a33) / dnr
    dvpy = det(a11, a12, a13, bvp1, bvp2, bvp3, a13, a23, a33) / dnr
    dvpz = det(a11, a12, a13, a12, a22, a23, bvp1, bvp2, bvp3) / dnr

    dwpx = det(bwp1, bwp2, bwp3, a12, a22, a23, a13, a23, a33) / dnr
    dwpy = det(a11, a12, a13, bwp1, bwp2, bwp3, a13, a23, a33) / dnr
    dwpz = det(a11, a12, a13, a12, a22, a23, bwp1, bwp2, bwp3) / dnr

    dtpx = det(btp1, btp2, btp3, a12, a22, a23, a13, a23, a33) / dnr
    dtpy = det(a11, a12, a13, btp1, btp2, btp3, a13, a23, a33) / dnr
    dtpz = det(a11, a12, a13, a12, a22, a23, btp1, btp2, btp3) / dnr

#ifdef KE_TURB
    dtkgx = det(btkg1, btkg2, btkg3, a12, a22, a23, a13, a23, a33) / dnr
    dtkgy = det(a11, a12, a13, btkg1, btkg2, btkg3, a13, a23, a33) / dnr
    dtkgz = det(a11, a12, a13, a12, a22, a23, btkg1, btkg2, btkg3) / dnr
	
    dtegx = det(bteg1, bteg2, bteg3, a12, a22, a23, a13, a23, a33) / dnr
    dtegy = det(a11, a12, a13, bteg1, bteg2, bteg3, a13, a23, a33) / dnr
    dtegz = det(a11, a12, a13, a12, a22, a23, bteg1, bteg2, bteg3) / dnr
#endif

#ifdef YP_P
    dyfgx = det(byfg1, byfg2, byfg3, a12, a22, a23, a13, a23, a33) / dnr
    dyfgy = det(a11, a12, a13, byfg1, byfg2, byfg3, a13, a23, a33) / dnr
    dyfgz = det(a11, a12, a13, a12, a22, a23, byfg1, byfg2, byfg3) / dnr
	
	dyogx = det(byog1, byog2, byog3, a12, a22, a23, a13, a23, a33) / dnr
    dyogy = det(a11, a12, a13, byog1, byog2, byog3, a13, a23, a33) / dnr
    dyogz = det(a11, a12, a13, a12, a22, a23, byog1, byog2, byog3) / dnr
	
	dypgx = det(bypg1, bypg2, bypg3, a12, a22, a23, a13, a23, a33) / dnr
    dypgy = det(a11, a12, a13, bypg1, bypg2, bypg3, a13, a23, a33) / dnr
    dypgz = det(a11, a12, a13, a12, a22, a23, bypg1, bypg2, bypg3) / dnr
#endif


!--------------------------------------------------------------------------------------------------------------
#ifndef AMUP_CONSTANT 
        if(phipf.gt.0.01d0)then
        gpart1=2.0d0*(dupx**2+dvpy**2+dwpz**2)
	gpart2=(dupy+dvpx)**2+(dupz+dwpx)**2+(dvpz+dwpy)**2
        gamapdot=min(max(sqrt(gpart1+gpart2),0.0001d0),10.0d0)	 
         
           if(gamapdot.lt.0.01d0)then
	     amup=2658.83d0*((2.0d0-0.814d0)+(0.814d0-1.0d0)*(gamapdot/0.01d0))+10.7d0*(2.0d0-(gamapdot/0.01d0))/0.01d0 
	   else
	     amup=10.7d0/gamapdot+2658.83d0*(gamapdot/0.01d0)**(0.814d0-1.0d0)
	   endif

	else
	   amup=0.001d0 !amug
	endif
		
		akp=amup*cpp/prp
#endif

!--------------------------------------------------------------------------------------------------------------	 
#ifdef KE_TURB
            amutg=0.0845d0*rogf*tkgf*tkgf/tegf
	    aktg=amutg*cpg/0.9d0
	    amug=amug+amutg
	    akg=akg+aktg
#endif

    ! Calculate delvg and stress terms
    delvg = dugx + dvgy + dwgz

    txxg = amug * (2.0d0 * dugx - c23 * delvg)
    tyyg = amug * (2.0d0 * dvgy - c23 * delvg)
    tzzg = amug * (2.0d0 * dwgz - c23 * delvg)
    txyg = amug * (dvgx + dugy)
    txzg = amug * (dwgx + dugz)
    tyzg = amug * (dvgz + dwgy)

    qxxg = ugf * txxg + vgf * txyg + wgf * txzg + akg * dtgx
    qyyg = ugf * txyg + vgf * tyyg + wgf * tyzg + akg * dtgy
    qzzg = ugf * txzg + vgf * tyzg + wgf * tzzg + akg * dtgz

    delvp = dupx + dvpy + dwpz
    txxp = amup * (2.0d0 * dupx - c23 * delvp)
    tyyp = amup * (2.0d0 * dvpy - c23 * delvp)
    tzzp = amup * (2.0d0 * dwpz - c23 * delvp)
    txyp = amup * (dvpx + dupy)
    txzp = amup * (dwpx + dupz)
    tyzp = amup * (dvpz + dwpy)

    qxxp = upf * txxp + vpf * txyp + wpf * txzp + akp * dtpx
    qyyp = upf * txyp + vpf * tyyp + wpf * tyzp + akp * dtpy
    qzzp = upf * txzp + vpf * tyzp + wpf * tzzp + akp * dtpz
    
#ifdef KE_TURB
    fvtkg=1.39d0*amug*dtkgx
    gvtkg=1.39d0*amug*dtkgy
    hvtkg=1.39d0*amug*dtkgz

    fvteg=1.39d0*amug*dtegx
    gvteg=1.39d0*amug*dtegy
    hvteg=1.39d0*amug*dtegz
#endif    

    ! Update viscosities
    visug(i) = visug(i) + phigf * (txxg * scx + txyg * scy + txzg * scz)
    visvg(i) = visvg(i) + phigf * (txyg * scx + tyyg * scy + tyzg * scz)
    viswg(i) = viswg(i) + phigf * (txzg * scx + tyzg * scy + tzzg * scz)
    vistg(i) = vistg(i) + phigf * (qxxg * scx + qyyg * scy + qzzg * scz)
    visup(i) = visup(i) + phipf * (txxp * scx + txyp * scy + txzp * scz)
    visvp(i) = visvp(i) + phipf * (txyp * scx + tyyp * scy + tyzp * scz)
    viswp(i) = viswp(i) + phipf * (txzp * scx + tyzp * scy + tzzp * scz)
    vistp(i) = vistp(i) + phipf * (qxxp * scx + qyyp * scy + qzzp * scz)

#ifdef KE_TURB
    vistkg(i) = vistkg(i) + phigf*(fvtkg*scx + gvtkg*scy + hvtkg*scz)
    visteg(i) = visteg(i) + phigf*(fvteg*scx + gvteg*scy + hvteg*scz)
#endif

#ifdef YP_P
  visyfg(i) = visyfg(i) + phigf*amug*(dyfgx*scx + dyfgy*scy + dyfgz*scz)
  visyog(i) = visyog(i) + phigf*amug*(dyogx*scx + dyogy*scy + dyogz*scz)
  visypg(i) = visypg(i) + phigf*amug*(dypgx*scx + dypgy*scy + dypgz*scz)
#endif

     endif
 end do
end do
!$acc end parallel loop    
end subroutine viscousgp

