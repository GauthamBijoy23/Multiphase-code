SUBROUTINE venkatakrishnan
    USE GlobalVariables
    USE MatrixOps
    IMPLICIT NONE
    external geometry
    external primitive
   
    INTEGER :: i,ne1,ne2,ne3,ne4,n1,n2,n3,n4,k,startl,endl
    REAL(DP) :: epssq,x0,y0,z0,xf1,xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4
    REAL(DP) :: xf10,xf20,xf30,xf40,yf10,yf20,yf30,yf40,zf10,zf20,zf30,zf40
    REAL(DP) :: delrogf1,delrogf2,delrogf3,delrogf4,delropf1,delropf2,delropf3,delropf4
    REAL(DP) :: delPgf1,delPgf2,delPgf3,delPgf4,delPpf1,delPpf2,delPpf3,delPpf4,phiPg1
    REAL(DP) :: phiPg2,phiPg3,phiPg4,venkatPpmax,venkatPpmin,phiPp1,phiPp2,phiPp3,phiPp4
    REAL(DP) :: venkatrogmax,venkatrogmin,venkatropmax,venkatropmin,phirog1,phirog2
    REAL(DP) :: phirog3,phirog4,phirop1,phirop2,phirop3,phirop4,venkatpgmax,venkatpgmin
    
!$acc parallel loop present(neles,pg(:),rog(:),pp(:),phi(:),nod(:,:), &
!$acc delpgz(:),delppx(:),delppy(:),delpgx(:),delpgy(:),delrogz(:), &
!$acc delropx(:),delropy(:),delppz(:),delrogx(:),delrogy(:),ycel(:), &
!$acc xcel(:),rop(:),nc4(:),nc2(:),dl(:),y(:),delropz(:),x(:),z(:), &
!$acc zcel(:),nc1(:),nc3(:)) &
!$acc private(delropf4,venkatropmin,n4,ne1,ne2,ne3,epssq,n1, &
!$acc n2,n3,x0,y0,z0,xf1, & 
!$acc xf2,xf3,xf4,yf1,yf2,yf3,yf4,zf1,zf2,zf3,zf4,delppf4, &
!$acc delrogf1,delrogf2,delrogf3,delrogf4,delropf1,delropf2,delropf3, &
!$acc delpgf1,delpgf2,delpgf3,delpgf4,delppf1,delppf2,delppf3,phirop4, &
!$acc venkatpgmax,ne4,phipg1,phipg2,phipg3,venkatpgmin,venkatppmax, &
!$acc phipg4,phipp1,phipp2,xf10,xf20,xf30,xf40,yf10,yf20,yf30,yf40,zf10,zf20,zf40,zf30, &
!$acc phipp3,venkatppmin,venkatrogmax,venkatrogmin,venkatropmax, &
!$acc phipp4,phirog1,phirog2,phirog3,phirog4,phirop1,phirop2,phirop3)

   do i=1,neles
if (iblank(i)/=1)cycle !-------(mod3) iblank check	 
	epssq=(0.3d0*dl(i))**3.0d0
	x0=xcel(i)
	y0=ycel(i)
	z0=zcel(i)
		
!Neighbours		
	ne1=nc1(i)
	ne2=nc2(i)
	ne3=nc3(i)
	ne4=nc4(i)
			
! Nodes are required to find del x and del y
	n1=nod(i,1)
	n2=nod(i,2)
	n3=nod(i,3)
	n4=nod(i,4)

!	 Nodes required for calculation of del r ie del x ;del y w.r.t element 0
         xf1=(x(n2)+x(n4)+x(n3))*0.3333d0    
         xf2=(x(n3)+x(n4)+x(n1))*0.3333d0    
         xf3=(x(n4)+x(n2)+x(n1))*0.3333d0    
         xf4=(x(n2)+x(n3)+x(n1))*0.3333d0    
         yf1=(y(n2)+y(n4)+y(n3))*0.3333d0    
         yf2=(y(n3)+y(n4)+y(n1))*0.3333d0    
         yf3=(y(n4)+y(n2)+y(n1))*0.3333d0    
         yf4=(y(n2)+y(n3)+y(n1))*0.3333d0    
         zf1=(z(n2)+z(n4)+z(n3))*0.3333d0    
         zf2=(z(n3)+z(n4)+z(n1))*0.3333d0    
         zf3=(z(n4)+z(n2)+z(n1))*0.3333d0    
         zf4=(z(n2)+z(n3)+z(n1))*0.3333d0    
			
 	 xf10=xf1-x0
	 xf20=xf2-x0
	 xf30=xf3-x0
	 xf40=xf4-x0
		
	 yf10=yf1-y0
	 yf20=yf2-y0
	 yf30=yf3-y0
	 yf40=yf4-y0
	 zf10=zf1-z0
	 zf20=zf2-z0
	 zf30=zf3-z0
	 zf40=zf4-z0

!del rog	is difference between face centre and element centre	
	delrogf1=xf10*delrogx(i)+yf10*delrogy(i)+zf10*delrogz(i)
	delrogf2=xf20*delrogx(i)+yf20*delrogy(i)+zf20*delrogz(i)
	delrogf3=xf30*delrogx(i)+yf30*delrogy(i)+zf30*delrogz(i)
	delrogf4=xf40*delrogx(i)+yf40*delrogy(i)+zf40*delrogz(i)
	
!	del rop	is difference between face centre and element centre	
	delropf1=xf10*delropx(i)+yf10*delropy(i)+zf10*delropz(i)
	delropf2=xf20*delropx(i)+yf20*delropy(i)+zf20*delropz(i)
	delropf3=xf30*delropx(i)+yf30*delropy(i)+zf30*delropz(i)
	delropf4=xf40*delropx(i)+yf40*delropy(i)+zf40*delropz(i)

!	del Pg is difference between face centre and element centre		
	delPgf1=xf10*delPgx(i)+yf10*delPgy(i)+zf10*delPgz(i)
	delPgf2=xf20*delPgx(i)+yf20*delPgy(i)+zf20*delPgz(i)
	delPgf3=xf30*delPgx(i)+yf30*delPgy(i)+zf30*delPgz(i)
	delPgf4=xf40*delPgx(i)+yf40*delPgy(i)+zf40*delPgz(i)
		
!	del Pp is difference between face centre and element centre		
	delPpf1=xf10*delPpx(i)+yf10*delPpy(i)+zf10*delPpz(i)
	delPpf2=xf20*delPpx(i)+yf20*delPpy(i)+zf20*delPpz(i)
	delPpf3=xf30*delPpx(i)+yf30*delPpy(i)+zf30*delPpz(i)
	delPpf4=xf40*delPpx(i)+yf40*delPpy(i)+zf40*delPpz(i)
			
	venkatPgmax=dmax1(Pg(i),Pg(ne1),Pg(ne2),Pg(ne3),Pg(ne4))
	venkatPgmin=dmin1(Pg(i),Pg(ne1),Pg(ne2),Pg(ne3),Pg(ne4))

!       Definition delPg	
		
	if(delPgf1.gt.0.0d0)then
		phiPg1=fvenkat(venkatPgmax-Pg(i),delPgf1,epssq)
	
	elseif(delPgf1.lt.0.0d0)then
		phiPg1=fvenkat(venkatPgmin-Pg(i),delPgf1,epssq)
	
	else
		phiPg1=1.0d0
	endif
	
	if(delPgf2.gt.0.0d0)then
		phiPg2=fvenkat(venkatPgmax-Pg(i),delPgf2,epssq)
			
	elseif(delPgf2.lt.0.0d0)then
		phiPg2=fvenkat(venkatPgmin-Pg(i),delPgf2,epssq)
		
	else
		phiPg2=1.0d0
	endif
			
	if(delPgf3.gt.0.0d0)then
		phiPg3=fvenkat(venkatPgmax-Pg(i),delPgf3,epssq)
			
	elseif(delPgf3.lt.0.0d0)then
		phiPg3=fvenkat(venkatPgmin-Pg(i),delPgf3,epssq)		
	else
	       phiPg3=1.0d0			
	endif
			
	if(delPgf4.gt.0.0d0)then
		phiPg4=fvenkat(venkatPgmax-Pg(i),delPgf4,epssq)
			
	elseif(delPgf4.lt.0.0d0)then
		phiPg4=fvenkat(venkatPgmin-Pg(i),delPgf4,epssq)		
	else
		phiPg4=1.0d0		
	endif
			
	venkatPpmax=dmax1(Pp(i),Pp(ne1),Pp(ne2),Pp(ne3),Pp(ne4))
	venkatPpmin=dmin1(Pp(i),Pp(ne1),Pp(ne2),Pp(ne3),Pp(ne4))

!	Definition delPp			
		
	if(delPpf1.gt.0.0d0)then	
		phiPp1=fvenkat(venkatPpmax-Pp(i),delPpf1,epssq)
		
	elseif(delPpf1.lt.0.0d0)then
		phiPp1=fvenkat(venkatPpmin-Pp(i),delPpf1,epssq)
	else
		phiPp1=1.0d0
	endif
			
	if(delPpf2.gt.0.0d0)then
		
		phiPp2=fvenkat(venkatPpmax-Pp(i),delPpf2,epssq)
			
	elseif(delPpf2.lt.0.0d0)then
		phiPp2=fvenkat(venkatPpmin-Pp(i),delPpf2,epssq)
				
	else
		phiPp2=1.0d0
	endif
			
	if(delPpf3.gt.0.0d0)then
			
		phiPp3=fvenkat(venkatPpmax-Pp(i),delPpf3,epssq)
			
	elseif(delPpf3.lt.0.0d0)then
			
		phiPp3=fvenkat(venkatPpmin-Pp(i),delPpf3,epssq)
				
	else
		phiPp3=1.0d0
	endif
!			
	if(delPpf4.gt.0.0d0)then
			
		phiPp4=fvenkat(venkatPpmax-Pp(i),delPpf4,epssq)
			
	elseif(delPpf4.lt.0.0d0)then
		phiPp4=fvenkat(venkatPpmin-Pp(i),delPpf4,epssq)
			
	else
		phiPp4=1.0d0
	endif

	venkatrogmax=dmax1(rog(i),rog(ne1),rog(ne2),rog(ne3),rog(ne4))
	venkatrogmin=dmin1(rog(i),rog(ne1),rog(ne2),rog(ne3),rog(ne4))
	venkatropmax=dmax1(rop(i),rop(ne1),rop(ne2),rop(ne3),rop(ne4))
	venkatropmin=dmin1(rop(i),rop(ne1),rop(ne2),rop(ne3),rop(ne4))
			
!       Definition delP			
	if(delrogf1.gt.0.0d0)then
	   phirog1=fvenkat(venkatrogmax-rog(i),delrogf1,epssq)
	elseif(delrogf1.lt.0.0d0)then
	   phirog1=fvenkat(venkatrogmin-rog(i),delrogf1,epssq)
	else
	   phirog1=1.0d0
	endif
	
	if(delrogf2.gt.0.0d0)then
	
	   phirog2=fvenkat(venkatrogmax-rog(i),delrogf2,epssq)
		
	elseif(delrogf2.lt.0.0d0)then
			
	   phirog2=fvenkat(venkatrogmin-rog(i),delrogf2,epssq)
	else
	   phirog2=1.0d0
	endif
	 
	if(delrogf3.gt.0.0d0)then
	   phirog3=fvenkat(venkatrogmax-rog(i),delrogf3,epssq)
		
	elseif(delrogf3.lt.0.0d0)then
	   phirog3=fvenkat(venkatrogmin-rog(i),delrogf3,epssq)
	else
	   phirog3=1.0d0
	endif
	
       if(delrogf4.gt.0.0d0)then
	  phirog4=fvenkat(venkatrogmax-rog(i),delrogf4,epssq)
       
       elseif(delrogf4.lt.0.0d0)then
	  phirog4=fvenkat(venkatrogmin-rog(i),delrogf4,epssq)
       else
	  phirog4=1.0d0
       endif
				
!	Definition delP			
	if(delropf1.gt.0.0d0)then
	    phirop1=fvenkat(venkatropmax-rop(i),delropf1,epssq)
	elseif(delropf1.lt.0.0d0)then
	    phirop1=fvenkat(venkatropmin-rop(i),delropf1,epssq)
	else
	    phirop1=1.0d0
	endif
	
	if(delropf2.gt.0.0d0)then
	   phirop2=fvenkat(venkatropmax-rop(i),delropf2,epssq)
	elseif(delropf2.lt.0.0d0)then
	   phirop2=fvenkat(venkatropmin-rop(i),delropf2,epssq)
        else
	   phirop2=1.0d0
	endif
	
	if(delropf3.gt.0.0d0)then
	   phirop3=fvenkat(venkatropmax-rop(i),delropf3,epssq)
	elseif(delropf3.lt.0.0d0)then
	    phirop3=fvenkat(venkatropmin-rop(i),delropf3,epssq)
        else
	    phirop3=1.0d0
	endif
	
	if(delropf4.gt.0.0d0)then
	  phirop4=fvenkat(venkatropmax-rop(i),delropf4,epssq)
	elseif(delropf4.lt.0.0d0)then
	  phirop4=fvenkat(venkatropmin-rop(i),delropf4,epssq)
	else
	   phirop4=1.0d0
	endif
	   phi(i)=dmin1(phipg1,phipg2,phipg3,phipg4)
	   phi(i)=dmin1(phi(i),phipp1,phipp2,phipp3,phipp4)
	   phi(i)=dmin1(phi(i),phirog1,phirog2,phirog3,phirog4)
	   phi(i)=dmin1(phi(i),phirop1,phirop2,phirop3,phirop4)
 
	  enddo
!$acc end parallel loop 


!!$acc parallel loop private(i)&
!!$acc present(neles,phi(:)

!do k = 1,2
!        if (k .eq. 1) then
!           startl = 1
!           endl   = neles
!        else
!           startl = neles + nghosts+ 1
!           endl   = ntot
!        end if

!        do i = startl, endl 
!          phi(i) = 1.0d0
!        enddo  
!end do        
!!$acc end parallel loop 

        
!$acc parallel loop private(i) &
!$acc present(nghosts,phi(:),nparent(:),nghost(:))        
       do i = 1, nghosts
	phi(nghost(i))=phi(nparent(i))
       enddo
!$acc end parallel loop
   return
end subroutine venkatakrishnan

