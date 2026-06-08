!============================================!
! Subroutine Name: Read_Flow                 !
!Read all the flow properties                !
!============================================! 				

SUBROUTINE Read_flow
    USE GlobalVariables
    
    CHARACTER(LEN=256) :: line
    INTEGER :: ios,unitnum,comment_pos
    
  unitnum = 9
  OPEN(unitnum, FILE="flow.in", STATUS="OLD", ACTION="READ", IOSTAT=ios)
  IF (ios /= 0) THEN
    PRINT *, "Error opening flow.in, IOSTAT=", ios
    STOP
  ENDIF

  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT  ! Exit loop at end of file or error
    IF (LEN_TRIM(line) == 0) CYCLE  ! Skip empty lines
    IF (line(1:1) == "#") CYCLE  ! Skip comment lines

    ! Find the position of '#' (start of a comment)
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)  ! Remove comment

    ! Read first two numbers in the cleaned line
    READ(line, *, IOSTAT=ios) PMAX,PMIN
    IF (ios == 0) EXIT  
  ENDDO

  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) tgmax, tgmin
    IF (ios == 0) EXIT
  ENDDO

  ! Repeat the above pattern for all variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) ugmax, ugmin
    IF (ios == 0) EXIT
  ENDDO
  
    DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) vgmax, vgmin
    IF (ios == 0) EXIT
  ENDDO
  
    DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) wgmax, wgmin
    IF (ios == 0) EXIT
  ENDDO
  
    ! Repeat the above pattern for all variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) upmax, upmin
    IF (ios == 0) EXIT
  ENDDO

  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) vpmax, vpmin
    IF (ios == 0) EXIT
  ENDDO

  ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) wpmax, wpmin
    IF (ios == 0) EXIT
  ENDDO
  
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) tpmax, tpmin
    IF (ios == 0) EXIT
  ENDDO
 

  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) pinf, gamap
    IF (ios == 0) EXIT
  ENDDO

  ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios)  gamag, cpp 
    IF (ios == 0) EXIT
  ENDDO

   ! Repeat the above pattern for all variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) cpg, amolg
    IF (ios == 0) EXIT
  ENDDO

  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) prg, epslnmax  
    IF (ios == 0) EXIT
  ENDDO

  ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) epslnmin, epslninf, phipzero  
    IF (ios == 0) EXIT
  ENDDO

     DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) cfl, niter, itersub
    IF (ios == 0) EXIT
  ENDDO
   
     ! Repeat the above pattern for all variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) istart, diap
    IF (ios == 0) EXIT
  ENDDO

  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) surfacetension, myflux
    IF (ios == 0) EXIT
  ENDDO

  ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) myorder, as, ts  
    IF (ios == 0) EXIT
  ENDDO
  
    ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) rhoq, rhor, rhos  
    IF (ios == 0) EXIT
  ENDDO


  ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios) radq, radr, rads
    IF (ios == 0) EXIT
  ENDDO
  
  
  ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios)  phiginit,phipinit,phiqinit, phirinit, phisinit
    IF (ios == 0) EXIT
  ENDDO
  
    ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios)  cp_t,cp_o,cp_p,amol_f,amol_o,amol_p
    IF (ios == 0) EXIT
  ENDDO
  
      ! Repeat for all remaining variables
  DO
    READ(unitnum, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
    comment_pos = INDEX(line, "#")
    IF (comment_pos > 0) line = line(1:comment_pos - 1)
    READ(line, *, IOSTAT=ios)  tolsp,sratio,hfp,preexp,ebyr
    IF (ios == 0) EXIT
  ENDDO


if(myid.eq.0) PRINT *, "Successfully read input parameters:", "phisinit", phisinit
!  PRINT *, "pmax=", pmax, " pmin=", pmin
!  PRINT *, "tgmax=", tgmax, " tgmin=", tgmin
!  PRINT *, "ugmax=", ugmax, " ugmin=", ugmin
!  PRINT *, "vgmax=", vgmax, " vgmin=", vgmin
!  PRINT *, "wgmax=", wgmax, " wgmin=", wgmin
!  PRINT *, "upmax=", upmax, " upmin=", upmin
!  PRINT *, "vpmax=", vpmax, " vpmin=", vpmin
!  PRINT *, "wpmax=", wpmax, " wpmin=", wpmin
!  PRINT *, "tpmax=", tpmax, " tpmin=", tpmin
!  PRINT *, "pinf=", pinf, " gamap=", gamap
!  PRINT *, "gamag=", gamag, " cpp=", cpp
!  PRINT *, "cpg=", cpg, " amolg=", amolg
!  PRINT *, "prg=", prg, " epslnmax=", epslnmax
!  PRINT *, "epslnmin=", epslnmin, " epslninf=", epslninf
!  PRINT *, "cfl=", cfl, " niter=", niter, " itersub=", itersub
!  PRINT *, "istart=", istart, " diap=", diap
!  PRINT *, "surfacetension=", surfacetension, " myflux=", myflux
 ! PRINT *, "myorder=", myorder, " as=", as, " ts=", ts

  END SUBROUTINE Read_flow
