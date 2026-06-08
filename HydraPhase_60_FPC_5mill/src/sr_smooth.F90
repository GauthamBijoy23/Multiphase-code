SUBROUTINE smooth
    USE GlobalVariables
    IMPLICIT NONE
    external geometry
    external primitive
   
    INTEGER :: i, k
    REAL(DP) :: ANUD,ANUN,ANUPG,ANUPP,D2,D24NC1,D24NC2,D24NC3,D24NC4
    REAL(DP) :: D4NC1,D2NC1,D2NC2,D2NC3,D2NC4,D4
    REAL(DP) :: D4NC2,D4NC3,D4NC4,E2NC1,E2NC2,E2NC3,E2NC4,E4NC1,E4NC2,E4NC3,E4NC4 

!=====================
! Loop through NGHOSTS
!=====================
!$acc parallel loop present(pp(:),pg(:),cn(:,:), &
!$acc nparent(:),nghost(:),nghosts,neg,nkg,npp,nyfg,nyog,nypg)  &
!$acc private(k,i)   

DO I = 1, NGHOSTS
    PG(NGHOST(I)) = PG(NPARENT(I))
    PP(NGHOST(I)) = PP(NPARENT(I))
!$acc loop seq   
    DO K = 1, 10 
        CN(NGHOST(I), K) = CN(NPARENT(I), K)
    END DO     

#ifdef KE_TURB
    CN(NGHOST(I), NKG) = CN(NPARENT(I), NKG)
    CN(NGHOST(I), NEG) = CN(NPARENT(I), NEG)
#endif

#ifdef NP_P
    CN(NGHOST(I), NPP) = CN(NPARENT(I), NPP)
#endif

#ifdef YP_P
       CN(NGHOST(I), NYFG) = CN(NPARENT(I), NYFG)
	CN(NGHOST(I), NYOG) = CN(NPARENT(I), NYOG)
	CN(NGHOST(I), NYPG) = CN(NPARENT(I), NYPG)
#endif
END DO
!$acc end parallel loop 


!=====================
! Loop through NELES
!=====================
!$acc parallel loop present(neles,cn(:,:),anu(:,:),d2q(:,:), &
!$acc nc1(:),nc2(:),nc4(:),nc3(:),pp(:),pg(:),neg,nkg,npp,nyfg,nyog,nypg) &      
!$acc private(anud,anun,anupp,anupg,k,i)

DO I = 1, NELES
    ANUN = ABS(PG(NC1(I)) - PG(I)) + ABS(PG(NC2(I)) - PG(I)) + ABS(PG(NC3(I)) - PG(I)) + ABS(PG(NC4(I)) - PG(I))
    ANUD = ABS(PG(NC1(I)) + PG(I)) + ABS(PG(NC2(I)) + PG(I)) + ABS(PG(NC3(I)) + PG(I)) + ABS(PG(NC4(I)) + PG(I))
    ANUPG = ANUN / ANUD

    ANUN = ABS(CN(NC1(I), 5) - CN(I, 5)) + ABS(CN(NC2(I), 5) - CN(I, 5)) + &
           ABS(CN(NC3(I), 5) - CN(I, 5)) + ABS(CN(NC4(I), 5) - CN(I, 5))
    ANUD = ABS(CN(NC1(I), 5) + CN(I, 5)) + ABS(CN(NC2(I), 5) + CN(I, 5)) + &
           ABS(CN(NC3(I), 5) + CN(I, 5)) + ABS(CN(NC4(I), 5) + CN(I, 5))

    ANU(I, 1) = ANUPG + ANUN / (ANUD + 1.0E-6)
  
    ANUN = ABS(CN(NC1(I), 1) - CN(I, 1)) + ABS(CN(NC2(I), 1) - CN(I, 1)) + &
           ABS(CN(NC3(I), 1) - CN(I, 1)) + ABS(CN(NC4(I), 1) - CN(I, 1))
    ANUD = ABS(CN(NC1(I), 1) + CN(I, 1)) + ABS(CN(NC2(I), 1) + CN(I, 1)) + &
           ABS(CN(NC3(I), 1) + CN(I, 1)) + ABS(CN(NC4(I), 1) + CN(I, 1))

    ANU(I, 1) = ANU(I, 1) + ANUN / (ANUD + 1.0E-6)

    ANUN = ABS(PP(NC1(I)) - PP(I)) + ABS(PP(NC2(I)) - PP(I)) + &
           ABS(PP(NC3(I)) - PP(I)) + ABS(PP(NC4(I)) - PP(I))
    ANUD = ABS(PP(NC1(I)) + PP(I)) + ABS(PP(NC2(I)) + PP(I)) + &
           ABS(PP(NC3(I)) + PP(I)) + ABS(PP(NC4(I)) + PP(I))

    ANUPP = ANUN / ANUD

    ANUN = ABS(CN(NC1(I), 6) - CN(I, 6)) + ABS(CN(NC2(I), 6) - CN(I, 6)) + &
           ABS(CN(NC3(I), 6) - CN(I, 6)) + ABS(CN(NC4(I), 6) - CN(I, 6))
    ANUD = ABS(CN(NC1(I), 6) + CN(I, 6)) + ABS(CN(NC2(I), 6) + CN(I, 6)) + &
           ABS(CN(NC3(I), 6) + CN(I, 6)) + ABS(CN(NC4(I), 6) + CN(I, 6))

    ANU(I, 2) = ANUPP + ANUN / (ANUD + 1.0E-6)

    ANUN = ABS(CN(NC1(I), 10) - CN(I, 10)) + ABS(CN(NC2(I), 10) - CN(I, 10)) + &
           ABS(CN(NC3(I), 10) - CN(I, 10)) + ABS(CN(NC4(I), 10) - CN(I, 10))
    ANUD = ABS(CN(NC1(I), 10) + CN(I, 10)) + ABS(CN(NC2(I), 10) + CN(I, 10)) + &
           ABS(CN(NC3(I), 10) + CN(I, 10)) + ABS(CN(NC4(I), 10) + CN(I, 10))

    ANU(I, 2) = ANU(I, 2) + ANUN / (ANUD + 1.0E-6)
    
!why this is coming bebore ke model 	
#ifdef NP_P
    ANUN = ABS(DLOG10(CN(NC1(I), npp)) - DLOG10(CN(I, npp))) + ABS(DLOG10(CN(NC2(I), npp)) - DLOG10(CN(I, npp))) + &
           ABS(DLOG10(CN(NC3(I), npp)) - DLOG10(CN(I, npp))) + ABS(DLOG10(CN(NC4(I), npp)) - DLOG10(CN(I, npp)))
    ANUD = ABS(DLOG10(CN(NC1(I), npp)) + DLOG10(CN(I, npp))) + ABS(DLOG10(CN(NC2(I), npp)) + DLOG10(CN(I, npp))) + &
           ABS(DLOG10(CN(NC3(I), npp)) + DLOG10(CN(I, npp))) + ABS(DLOG10(CN(NC4(I), npp)) + DLOG10(CN(I, npp)))

    ANU(I, 2) = ANU(I, 2) + ANUN / (ANUD + 1.0E-6)
#endif	    
    
    
#ifdef KE_TURB	
    ANUN = ABS(CN(NC1(I), NKG) - CN(I, NKG)) + ABS(CN(NC2(I), NKG) - CN(I, NKG)) + &
           ABS(CN(NC3(I), NKG) - CN(I, NKG)) + ABS(CN(NC4(I), NKG) - CN(I, NKG))
    ANUD = ABS(CN(NC1(I), NKG) + CN(I, NKG)) + ABS(CN(NC2(I), NKG) + CN(I, NKG)) + &
           ABS(CN(NC3(I), NKG) + CN(I, NKG)) + ABS(CN(NC4(I), NKG) + CN(I, NKG))


    ANU(I, 3) = ANUPG + ANUN / (ANUD + 1.0E-6)
  
    ANUN = ABS(CN(NC1(I), NEG) - CN(I, NEG)) + ABS(CN(NC2(I), NEG) - CN(I, NEG)) + &
           ABS(CN(NC3(I), NEG) - CN(I, NEG)) + ABS(CN(NC4(I), NEG) - CN(I, NEG))
    ANUD = ABS(CN(NC1(I), NEG) + CN(I, NEG)) + ABS(CN(NC2(I), NEG) + CN(I, NEG)) + &
           ABS(CN(NC3(I), NEG) + CN(I, NEG)) + ABS(CN(NC4(I), NEG) + CN(I, NEG))

    ANU(I, 3) = ANU(I, 3) + ANUN / (ANUD + 1.0E-6)
#endif	

#ifdef YP_P
    ANUN = ABS(CN(NC1(I), NYFG) - CN(I, NYFG)) + ABS(CN(NC2(I), NYFG) - CN(I, NYFG)) + &
           ABS(CN(NC3(I), NYFG) - CN(I, NYFG)) + ABS(CN(NC4(I), NYFG) - CN(I, NYFG))
    ANUD = ABS(CN(NC1(I), NYFG) + CN(I, NYFG)) + ABS(CN(NC2(I), NYFG) + CN(I, NYFG)) + &
           ABS(CN(NC3(I), NYFG) + CN(I, NYFG)) + ABS(CN(NC4(I), NYFG) + CN(I, NYFG))


    ANU(I, 4) = ANUPG + ANUN / (ANUD + 1.0E-6)
#endif

!$acc loop seq 
    DO K = 1, 10
        D2Q(I, K) = CN(NC1(I), K) - CN(I, K) + CN(NC2(I), K) - CN(I, K) + &
                    CN(NC3(I), K) - CN(I, K) + CN(NC4(I), K) - CN(I, K)
    END DO
#ifdef KE_TURB	
	D2Q(I, NKG) = CN(NC1(I), NKG) - CN(I, NKG) + CN(NC2(I), NKG) - CN(I, NKG) + &
                      CN(NC3(I), NKG) - CN(I, NKG) + CN(NC4(I), NKG) - CN(I, NKG)
					  
	D2Q(I, NEG) = CN(NC1(I), NEG) - CN(I, NEG) + CN(NC2(I), NEG) - CN(I, NEG) + &
                      CN(NC3(I), NEG) - CN(I, NEG) + CN(NC4(I), NEG) - CN(I, NEG)			  
#endif
#ifdef NP_P
     D2Q(I, NPP) = CN(NC1(I), NPP) - CN(I, NPP) + CN(NC2(I), NPP) - CN(I, NPP) + &
                   CN(NC3(I), NPP) - CN(I, NPP) + CN(NC4(I), NPP) - CN(I, NPP)
#endif
#ifdef YP_P
     D2Q(I, NYFG) = CN(NC1(I), NYFG) - CN(I, NYFG) + CN(NC2(I), NYFG) - CN(I, NYFG) + &
                    CN(NC3(I), NYFG) - CN(I, NYFG) + CN(NC4(I), NYFG) - CN(I, NYFG)
	 D2Q(I, NYOG) = CN(NC1(I), NYOG) - CN(I, NYOG) + CN(NC2(I), NYOG) - CN(I, NYOG) + &
                    CN(NC3(I), NYOG) - CN(I, NYOG) + CN(NC4(I), NYOG) - CN(I, NYOG)
     D2Q(I, NYPG) = CN(NC1(I), NYPG) - CN(I, NYPG) + CN(NC2(I), NYPG) - CN(I, NYPG) + &
                    CN(NC3(I), NYPG) - CN(I, NYPG) + CN(NC4(I), NYPG) - CN(I, NYPG)					
#endif
END DO
!$acc end parallel loop


!=====================
! Loop through NGHOSTS
!===================== 
!$acc parallel loop present(d2q(:,:),anu(:,:), &
!$acc nparent(:),nghost(:),nghosts,neg,nkg,npp,nyfg,nyog,nypg) &
!$acc private(i,k) 
DO I = 1, NGHOSTS
!$acc loop seq 
    DO K = 1, 10
        D2Q(NGHOST(I), K) = D2Q(NPARENT(I), K)
    END DO
    ANU(NGHOST(I), 1) = ANU(NPARENT(I), 1)
    ANU(NGHOST(I), 2) = ANU(NPARENT(I), 2)
    
#ifdef KE_TURB
        D2Q(NGHOST(I), NKG) = D2Q(NPARENT(I), NKG)
	D2Q(NGHOST(I), NEG) = D2Q(NPARENT(I), NEG)
	ANU(NGHOST(I), 3) = ANU(NPARENT(I), 3)	
#endif
#ifdef NP_P
    D2Q(NGHOST(I), NPP) = D2Q(NPARENT(I), NPP)
#endif

#ifdef YP_P
    D2Q(NGHOST(I), NYFG) = D2Q(NPARENT(I), NYFG)
	D2Q(NGHOST(I), NYOG) = D2Q(NPARENT(I), NYOG)
	D2Q(NGHOST(I), NYPG) = D2Q(NPARENT(I), NYPG)
	ANU(NGHOST(I), 4) = ANU(NPARENT(I), 4)	
#endif	    
END DO
!$acc end parallel loop 

!==============================================
! Loop through NELES for the final calculations
!==============================================
!$acc parallel loop present(neles,co2,co4,dt,nc3(:),nc1(:), &
!$acc d2q(:,:),vol(:),nc4(:),nc2(:),d24(:,:), &
!$acc anu(:,:),cn(:,:),neg,nkg,npp,nypg,nyfg) &
!$acc private(e2nc1,e2nc2,e2nc3,e2nc4,e4nc1,e4nc2,e4nc3,d4nc4, &
!$acc d24nc1,d24nc2,d24nc3,d24nc4,d2nc1,d2nc2,d2nc3,d2nc4,d4nc1, &
!$acc d4nc2,d4nc3,e4nc4,d4,d2,i,k)

DO I = 1, NELES
    E2NC1 = CO2 * MAX(ANU(I, 1), ANU(NC1(I), 1))
    E2NC2 = CO2 * MAX(ANU(I, 1), ANU(NC2(I), 1))
    E2NC3 = CO2 * MAX(ANU(I, 1), ANU(NC3(I), 1))
    E2NC4 = CO2 * MAX(ANU(I, 1), ANU(NC4(I), 1))
    E4NC1 = MAX(0.0D0, CO4 - E2NC1)
    E4NC2 = MAX(0.0D0, CO4 - E2NC2)
    E4NC3 = MAX(0.0D0, CO4 - E2NC3)
    E4NC4 = MAX(0.0D0, CO4 - E2NC4)

    D24NC1 = 0.5D0 * (VOL(I) / DT + VOL(NC1(I)) / DT)
    D24NC2 = 0.5D0 * (VOL(I) / DT + VOL(NC2(I)) / DT)
    D24NC3 = 0.5D0 * (VOL(I) / DT + VOL(NC3(I)) / DT)
    D24NC4 = 0.5D0 * (VOL(I) / DT + VOL(NC4(I)) / DT)

    D2NC1 = E2NC1 * D24NC1
    D2NC2 = E2NC2 * D24NC2
    D2NC3 = E2NC3 * D24NC3
    D2NC4 = E2NC4 * D24NC4
    D4NC1 = E4NC1 * D24NC1
    D4NC2 = E4NC2 * D24NC2
    D4NC3 = E4NC3 * D24NC3
    D4NC4 = E4NC4 * D24NC4

!$acc loop seq     
    ! Calculate D2 and D4 terms
    DO K = 1, 5
        D2 = D2NC1 * (CN(NC1(I), K) - CN(I, K)) + &
             D2NC2 * (CN(NC2(I), K) - CN(I, K)) + &
             D2NC3 * (CN(NC3(I), K) - CN(I, K)) + &
             D2NC4 * (CN(NC4(I), K) - CN(I, K))

        D4 = D4NC1 * (D2Q(NC1(I), K) - D2Q(I, K)) + &
             D4NC2 * (D2Q(NC2(I), K) - D2Q(I, K)) + &
             D4NC3 * (D2Q(NC3(I), K) - D2Q(I, K)) + &
             D4NC4 * (D2Q(NC4(I), K) - D2Q(I, K))

        D24(I, K) = D2 - D4
    END DO
    
    E2NC1 = CO2 * MAX(ANU(I, 2), ANU(NC1(I), 2))
    E2NC2 = CO2 * MAX(ANU(I, 2), ANU(NC2(I), 2))
    E2NC3 = CO2 * MAX(ANU(I, 2), ANU(NC3(I), 2))
    E2NC4 = CO2 * MAX(ANU(I, 2), ANU(NC4(I), 2))
    E4NC1 = MAX(0.0D0, CO4 - E2NC1)
    E4NC2 = MAX(0.0D0, CO4 - E2NC2)
    E4NC3 = MAX(0.0D0, CO4 - E2NC3)
    E4NC4 = MAX(0.0D0, CO4 - E2NC4)

    D24NC1 = 0.5D0 * (VOL(I) / DT + VOL(NC1(I)) / DT)
    D24NC2 = 0.5D0 * (VOL(I) / DT + VOL(NC2(I)) / DT)
    D24NC3 = 0.5D0 * (VOL(I) / DT + VOL(NC3(I)) / DT)
    D24NC4 = 0.5D0 * (VOL(I) / DT + VOL(NC4(I)) / DT)

    D2NC1 = E2NC1 * D24NC1
    D2NC2 = E2NC2 * D24NC2
    D2NC3 = E2NC3 * D24NC3
    D2NC4 = E2NC4 * D24NC4
    D4NC1 = E4NC1 * D24NC1
    D4NC2 = E4NC2 * D24NC2
    D4NC3 = E4NC3 * D24NC3
    D4NC4 = E4NC4 * D24NC4

!$acc loop seq      
    ! Calculate D2 and D4 terms for the second set of variables
    DO K = 6, 10
        D2 = D2NC1 * (CN(NC1(I), K) - CN(I, K)) + &
             D2NC2 * (CN(NC2(I), K) - CN(I, K)) + &
             D2NC3 * (CN(NC3(I), K) - CN(I, K)) + &
             D2NC4 * (CN(NC4(I), K) - CN(I, K))

        D4 = D4NC1 * (D2Q(NC1(I), K) - D2Q(I, K)) + &
             D4NC2 * (D2Q(NC2(I), K) - D2Q(I, K)) + &
             D4NC3 * (D2Q(NC3(I), K) - D2Q(I, K)) + &
             D4NC4 * (D2Q(NC4(I), K) - D2Q(I, K))

        D24(I, K) = D2 - D4
       ! if(total.eq.2.and.myid.eq.1) print*, i,D24(I,K)
    END DO  
#ifdef NP_P
        D2 = D2NC1 * (CN(NC1(I), NPP) - CN(I, NPP)) + &
             D2NC2 * (CN(NC2(I), NPP) - CN(I, NPP)) + &
             D2NC3 * (CN(NC3(I), NPP) - CN(I, NPP)) + &
             D2NC4 * (CN(NC4(I), NPP) - CN(I, NPP))

        D4 = D4NC1 * (D2Q(NC1(I), NPP) - D2Q(I, NPP)) + &
             D4NC2 * (D2Q(NC2(I), NPP) - D2Q(I, NPP)) + &
             D4NC3 * (D2Q(NC3(I), NPP) - D2Q(I, NPP)) + &
             D4NC4 * (D2Q(NC4(I), NPP) - D2Q(I, NPP))

        D24(I, NPP) = D2 - D4
#endif

#ifdef KE_TURB
    E2NC1 = CO2 * MAX(ANU(I, 3), ANU(NC1(I), 3))
    E2NC2 = CO2 * MAX(ANU(I, 3), ANU(NC2(I), 3))
    E2NC3 = CO2 * MAX(ANU(I, 3), ANU(NC3(I), 3))
    E2NC4 = CO2 * MAX(ANU(I, 3), ANU(NC4(I), 3))
    E4NC1 = MAX(0.0D0, CO4 - E2NC1)
    E4NC2 = MAX(0.0D0, CO4 - E2NC2)
    E4NC3 = MAX(0.0D0, CO4 - E2NC3)
    E4NC4 = MAX(0.0D0, CO4 - E2NC4)

    D24NC1 = 0.5D0 * (VOL(I) / DT + VOL(NC1(I)) / DT)
    D24NC2 = 0.5D0 * (VOL(I) / DT + VOL(NC2(I)) / DT)
    D24NC3 = 0.5D0 * (VOL(I) / DT + VOL(NC3(I)) / DT)
    D24NC4 = 0.5D0 * (VOL(I) / DT + VOL(NC4(I)) / DT)

    D2NC1 = E2NC1 * D24NC1
    D2NC2 = E2NC2 * D24NC2
    D2NC3 = E2NC3 * D24NC3
    D2NC4 = E2NC4 * D24NC4
    D4NC1 = E4NC1 * D24NC1
    D4NC2 = E4NC2 * D24NC2
    D4NC3 = E4NC3 * D24NC3
    D4NC4 = E4NC4 * D24NC4

!$acc loop seq     
    ! Calculate D2 and D4 terms
    DO K = NKG, NEG
        D2 = D2NC1 * (CN(NC1(I), K) - CN(I, K)) + &
             D2NC2 * (CN(NC2(I), K) - CN(I, K)) + &
             D2NC3 * (CN(NC3(I), K) - CN(I, K)) + &
             D2NC4 * (CN(NC4(I), K) - CN(I, K))

        D4 = D4NC1 * (D2Q(NC1(I), K) - D2Q(I, K)) + &
             D4NC2 * (D2Q(NC2(I), K) - D2Q(I, K)) + &
             D4NC3 * (D2Q(NC3(I), K) - D2Q(I, K)) + &
             D4NC4 * (D2Q(NC4(I), K) - D2Q(I, K))

        D24(I, K) = D2 - D4
    END DO 
#endif  

#ifdef YP_P
    E2NC1 = CO2 * MAX(ANU(I, 4), ANU(NC1(I), 4))
    E2NC2 = CO2 * MAX(ANU(I, 4), ANU(NC2(I), 4))
    E2NC3 = CO2 * MAX(ANU(I, 4), ANU(NC3(I), 4))
    E2NC4 = CO2 * MAX(ANU(I, 4), ANU(NC4(I), 4))
    E4NC1 = MAX(0.0D0, CO4 - E2NC1)
    E4NC2 = MAX(0.0D0, CO4 - E2NC2)
    E4NC3 = MAX(0.0D0, CO4 - E2NC3)
    E4NC4 = MAX(0.0D0, CO4 - E2NC4)

    D24NC1 = 0.5D0 * (VOL(I) / DT + VOL(NC1(I)) / DT)
    D24NC2 = 0.5D0 * (VOL(I) / DT + VOL(NC2(I)) / DT)
    D24NC3 = 0.5D0 * (VOL(I) / DT + VOL(NC3(I)) / DT)
    D24NC4 = 0.5D0 * (VOL(I) / DT + VOL(NC4(I)) / DT)

    D2NC1 = E2NC1 * D24NC1
    D2NC2 = E2NC2 * D24NC2
    D2NC3 = E2NC3 * D24NC3
    D2NC4 = E2NC4 * D24NC4
    D4NC1 = E4NC1 * D24NC1
    D4NC2 = E4NC2 * D24NC2
    D4NC3 = E4NC3 * D24NC3
    D4NC4 = E4NC4 * D24NC4

!$acc loop seq     
    ! Calculate D2 and D4 terms
    DO K = NYFG, NYPG
        D2 = D2NC1 * (CN(NC1(I), K) - CN(I, K)) + &
             D2NC2 * (CN(NC2(I), K) - CN(I, K)) + &
             D2NC3 * (CN(NC3(I), K) - CN(I, K)) + &
             D2NC4 * (CN(NC4(I), K) - CN(I, K))

        D4 = D4NC1 * (D2Q(NC1(I), K) - D2Q(I, K)) + &
             D4NC2 * (D2Q(NC2(I), K) - D2Q(I, K)) + &
             D4NC3 * (D2Q(NC3(I), K) - D2Q(I, K)) + &
             D4NC4 * (D2Q(NC4(I), K) - D2Q(I, K))

        D24(I, K) = D2 - D4
    END DO 
#endif              
END DO
!$acc end parallel loop

   return
end subroutine smooth

