subroutine dissipationqrs
    USE GlobalVariables
    IMPLICIT NONE

INTEGER :: i, k                     ! Loop indices

REAL(DP) :: anun, anud, anupp         ! Intermediate smoothness indicators
REAL(DP) :: d2, d4                    ! Second and fourth-order diffusion terms
REAL(DP) :: d24nc1, d24nc2, d24nc3, d24nc4
REAL(DP) :: e2nc1, e2nc2, e2nc3, e2nc4
REAL(DP) :: e4nc1, e4nc2, e4nc3, e4nc4
REAL(DP) :: d2nc1, d2nc2, d2nc3, d2nc4
REAL(DP) :: d4nc1, d4nc2, d4nc3, d4nc4

!$acc parallel loop present(nghost(:),nparent(:),pp(:),nghosts,cn(:,:)) &
!$acc private(i,k)
      DO I = 1, NGHOSTS
         PP(NGHOST(I))=PP(NPARENT(I))
!$acc loop seq         
         DO K = 11, 22
            CN(NGHOST(I),K)=CN(NPARENT(I),K)
         ENDDO
      ENDDO
!$acc end parallel loop      
		
!$acc parallel loop present(neles,cn(:,:),anu(:,:),d2q(:,:),nc1(:),nc2(:),pp(:),nc4(:), &
!$acc nc3(:)) private(anud,anupp,anun)

      DO I = 1, NELES
		
	ANUN=DABS(PP(NC1(I))-PP(I))+DABS(PP(NC2(I))-PP(I))+DABS(PP(NC3(I))-PP(I))+DABS(PP(NC4(I))-PP(I))
        ANUD=DABS(PP(NC1(I))+PP(I))+DABS(PP(NC2(I))+PP(I))+DABS(PP(NC3(I))+PP(I))+DABS(PP(NC4(I))+PP(I))
        ANUPP=ANUN/ANUD
		
	ANUN=DABS(CN(NC1(I),6)-CN(I,6))+DABS(CN(NC2(I),6)-CN(I,6))+DABS(CN(NC3(I),6)-CN(I,6))+DABS(CN(NC4(I),6)-CN(I,6))
        ANUD=DABS(CN(NC1(I),6)+CN(I,6))+DABS(CN(NC2(I),6)+CN(I,6))+DABS(CN(NC3(I),6)+CN(I,6))+DABS(CN(NC4(I),6)+CN(I,6))
        ANUPP=ANUPP+ANUN/(ANUD+1.0E-16)
		
	ANUN=DABS(CN(NC1(I),11)-CN(I,11))+DABS(CN(NC2(I),11)-CN(I,11))+DABS(CN(NC3(I),11)-CN(I,11))+DABS(CN(NC4(I),11)-CN(I,11))
        ANUD=DABS(CN(NC1(I),11)+CN(I,11))+DABS(CN(NC2(I),11)+CN(I,11))+DABS(CN(NC3(I),11)+CN(I,11))+DABS(CN(NC4(I),11)+CN(I,11))
        ANU(I,1)=ANUPP+ANUN/(ANUD+1.0E-16)
		
	ANUN=DABS(CN(NC1(I),15)-CN(I,15))+DABS(CN(NC2(I),15)-CN(I,15))+DABS(CN(NC3(I),15)-CN(I,15))+DABS(CN(NC4(I),15)-CN(I,15))
        ANUD=DABS(CN(NC1(I),15)+CN(I,15))+DABS(CN(NC2(I),15)+CN(I,15))+DABS(CN(NC3(I),15)+CN(I,15))+DABS(CN(NC4(I),15)+CN(I,15))
        ANU(I,2)=ANUPP+ANUN/(ANUD+1.0E-16)
		
	ANUN=DABS(CN(NC1(I),19)-CN(I,19))+DABS(CN(NC2(I),19)-CN(I,19))+DABS(CN(NC3(I),19)-CN(I,19))+DABS(CN(NC4(I),19)-CN(I,19))
        ANUD=DABS(CN(NC1(I),19)+CN(I,19))+DABS(CN(NC2(I),19)+CN(I,19))+DABS(CN(NC3(I),19)+CN(I,19))+DABS(CN(NC4(I),19)+CN(I,19))
        ANU(I,3)=ANUPP+ANUN/(ANUD+1.0E-16)

!$acc loop seq  		
	DO K = 11, 22
           D2Q(I,K)=CN(NC1(I),K)-CN(I,K)+CN(NC2(I),K)-CN(I,K)+CN(NC3(I),K)-CN(I,K)+CN(NC4(I),K)-CN(I,K)
        ENDDO
     ENDDO
!$acc end parallel loop

!$acc parallel loop present(d2q(:,:),anu(:,:),nparent(:),nghost(:),nghosts) private(i,k)
      DO I = 1, NGHOSTS
!$acc loop seq       
         DO K = 11, 22
           D2Q(NGHOST(I),K)=D2Q(NPARENT(I),K)
        ENDDO
          ANU(NGHOST(I),1)=ANU(NPARENT(I),1)
          ANU(NGHOST(I),2)=ANU(NPARENT(I),2)
          ANU(NGHOST(I),3)=ANU(NPARENT(I),3)
     ENDDO
!$acc end parallel loop     

!$acc parallel loop present(neles,anu(:,:),d24(:,:),cn(:,:),nc3(:),nc2(:), &
!$acc nc1(:),vol(:),nc4(:),d2q(:,:),co4,dt,co2) &
!$acc private(d4nc4,d24nc1,d24nc2,d24nc3,e2nc1,e2nc2,e2nc3,e2nc4,e4nc1,e4nc2, &
!$acc e4nc3,e4nc4,d24nc4,d2nc1,d2nc2,d2nc3,d2nc4,d4nc1,d4nc2,d4nc3,d4,d2,i,k)
		
     DO I = 1, NELES

          D24NC1=0.5D0*(VOL(I)/DT+VOL(NC1(I))/DT)
          D24NC2=0.5D0*(VOL(I)/DT+VOL(NC2(I))/DT)
          D24NC3=0.5D0*(VOL(I)/DT+VOL(NC3(I))/DT)
          D24NC4=0.5D0*(VOL(I)/DT+VOL(NC4(I))/DT)

          E2NC1=CO2*DMAX1(ANU(I,1),ANU(NC1(I),1))
          E2NC2=CO2*DMAX1(ANU(I,1),ANU(NC2(I),1))
          E2NC3=CO2*DMAX1(ANU(I,1),ANU(NC3(I),1))
          E2NC4=CO2*DMAX1(ANU(I,1),ANU(NC4(I),1))
          E4NC1=DMAX1(0.0D0,(CO4-E2NC1))
          E4NC2=DMAX1(0.0D0,(CO4-E2NC2))
          E4NC3=DMAX1(0.0D0,(CO4-E2NC3))
          E4NC4=DMAX1(0.0D0,(CO4-E2NC4))
          D2NC1=E2NC1*D24NC1
          D2NC2=E2NC2*D24NC2
          D2NC3=E2NC3*D24NC3
          D2NC4=E2NC4*D24NC4
          D4NC1=E4NC1*D24NC1
          D4NC2=E4NC2*D24NC2
          D4NC3=E4NC3*D24NC3
          D4NC4=E4NC4*D24NC4
!$acc loop seq       
        DO K = 11, 14
        D2 = D2NC1*(CN(NC1(I),K)-CN(I,K)) + D2NC2*(CN(NC2(I),K)-CN(I,K)) + D2NC3*(CN(NC3(I),K)-CN(I,K)) + &
             D2NC4*(CN(NC4(I),K)-CN(I,K))
        D4 = D4NC1*(D2Q(NC1(I),K)-D2Q(I,K)) + D4NC2*(D2Q(NC2(I),K)-D2Q(I,K)) + D4NC3*(D2Q(NC3(I),K) - &
             D2Q(I,K))+D4NC4*(D2Q(NC4(I),K) - D2Q(I,K))
        D24(I,K)=D2-D4
        ENDDO
		
	E2NC1=CO2*DMAX1(ANU(I,2),ANU(NC1(I),2))
        E2NC2=CO2*DMAX1(ANU(I,2),ANU(NC2(I),2))
        E2NC3=CO2*DMAX1(ANU(I,2),ANU(NC3(I),2))
        E2NC4=CO2*DMAX1(ANU(I,2),ANU(NC4(I),2))
        E4NC1=DMAX1(0.0D0,(CO4-E2NC1))
        E4NC2=DMAX1(0.0D0,(CO4-E2NC2))
        E4NC3=DMAX1(0.0D0,(CO4-E2NC3))
        E4NC4=DMAX1(0.0D0,(CO4-E2NC4))
        D2NC1=E2NC1*D24NC1
        D2NC2=E2NC2*D24NC2
        D2NC3=E2NC3*D24NC3
        D2NC4=E2NC4*D24NC4
        D4NC1=E4NC1*D24NC1
        D4NC2=E4NC2*D24NC2
        D4NC3=E4NC3*D24NC3
        D4NC4=E4NC4*D24NC4
!$acc loop seq          
        DO K = 15, 18
        D2 = D2NC1*(CN(NC1(I),K)-CN(I,K)) + D2NC2*(CN(NC2(I),K)-CN(I,K)) + D2NC3*(CN(NC3(I),K)-CN(I,K)) + &
             D2NC4*(CN(NC4(I),K)-CN(I,K))
        D4 = D4NC1*(D2Q(NC1(I),K)-D2Q(I,K)) + D4NC2*(D2Q(NC2(I),K)-D2Q(I,K)) + D4NC3*(D2Q(NC3(I),K)-D2Q(I,K)) + &
             D4NC4*(D2Q(NC4(I),K) - D2Q(I,K))
        D24(I,K)=D2-D4
        ENDDO
		
		
	E2NC1=CO2*DMAX1(ANU(I,3),ANU(NC1(I),3))
        E2NC2=CO2*DMAX1(ANU(I,3),ANU(NC2(I),3))
        E2NC3=CO2*DMAX1(ANU(I,3),ANU(NC3(I),3))
        E2NC4=CO2*DMAX1(ANU(I,3),ANU(NC4(I),3))
        E4NC1=DMAX1(0.0D0,(CO4-E2NC1))
        E4NC2=DMAX1(0.0D0,(CO4-E2NC2))
        E4NC3=DMAX1(0.0D0,(CO4-E2NC3))
        E4NC4=DMAX1(0.0D0,(CO4-E2NC4))
        D2NC1=E2NC1*D24NC1
        D2NC2=E2NC2*D24NC2
        D2NC3=E2NC3*D24NC3
        D2NC4=E2NC4*D24NC4
        D4NC1=E4NC1*D24NC1
        D4NC2=E4NC2*D24NC2
        D4NC3=E4NC3*D24NC3
        D4NC4=E4NC4*D24NC4
!$acc loop seq          
        DO K = 19, 22
        D2 = D2NC1*(CN(NC1(I),K)-CN(I,K)) + D2NC2*(CN(NC2(I),K)-CN(I,K)) + D2NC3*(CN(NC3(I),K)-CN(I,K)) + &
             D2NC4*(CN(NC4(I),K)-CN(I,K))
        D4 = D4NC1*(D2Q(NC1(I),K)-D2Q(I,K)) + D4NC2*(D2Q(NC2(I),K)-D2Q(I,K)) + D4NC3*(D2Q(NC3(I),K)-D2Q(I,K)) + & 
             D4NC4*(D2Q(NC4(I),K)-D2Q(I,K))
        D24(I,K)=D2-D4
        ENDDO
      ENDDO  
!$acc end parallel loop		
end subroutine dissipationqrs
