PROGRAM merge_grids
    IMPLICIT NONE

    ! --- scalars ---
    INTEGER :: bg_nodes, bg_neles, bg_nghosts, bg_pghosts, bg_maxpg
    INTEGER :: ov_nodes, ov_neles, ov_nghosts, ov_pghosts, ov_maxpg
    INTEGER :: tot_nodes, tot_neles, tot_nghosts, tot_pghosts, tot_maxpg

    ! --- temp variables for reading one line at a time ---
    INTEGER :: i, n
    INTEGER :: nel, n1, n2, n3, n4, c1, c2, c3, c4
    INTEGER :: np, ng, nt, ns
    INTEGER :: pp, d1, d2, d3, d4, d5, d6, dst, g, ib, ibc
    REAL(8) :: x, y, z

    CHARACTER(LEN=256) :: indir, outdir

    indir  = '/home/gautham_linux/Multiphase_code/HydraPhase_60_FPC_5mill/inputfiles/'
    outdir = '/home/gautham_linux/Multiphase_code/HydraPhase_60_FPC_5mill/geometry_files/'

    ! =========================================================
    ! STEP 1 — read headers to get counts
    ! =========================================================

    OPEN(10, file=trim(indir)//'grid_bg.dat',  status='old')
    OPEN(11, file=trim(indir)//'grid_ovr.dat', status='old')
    OPEN(12, file=trim(indir)//'bc_bg.in',    status='old')
    OPEN(13, file=trim(indir)//'bc_ovr.in',   status='old')
    OPEN(14, file=trim(indir)//'proc_bg.in',  status='old')
    OPEN(15, file=trim(indir)//'proc_ovr.in', status='old')

    READ(10,*) bg_nodes, bg_neles
    READ(11,*) ov_nodes, ov_neles
    READ(12,*) bg_nghosts
    READ(13,*) ov_nghosts
    READ(14,*) bg_pghosts, bg_maxpg
    READ(15,*) ov_pghosts, ov_maxpg

    tot_nodes   = bg_nodes   + ov_nodes
    tot_neles   = bg_neles   + ov_neles
    tot_nghosts = bg_nghosts + ov_nghosts
    tot_pghosts = bg_pghosts + ov_pghosts
    tot_maxpg   = MAX(bg_maxpg, ov_maxpg)

    PRINT*, 'BG  : nodes=', bg_nodes,  ' neles=', bg_neles,  ' nghosts=', bg_nghosts,  ' pghosts=', bg_pghosts
    PRINT*, 'OVR : nodes=', ov_nodes,  ' neles=', ov_neles,  ' nghosts=', ov_nghosts,  ' pghosts=', ov_pghosts
    PRINT*, 'OUT : nodes=', tot_nodes, ' neles=', tot_neles, ' nghosts=', tot_nghosts, ' pghosts=', tot_pghosts

    ! =========================================================
    ! STEP 2 — write merged grid file
    ! =========================================================

    OPEN(20, file=trim(outdir)//'grid_proc0000.dat', status='replace')

    WRITE(20,*) tot_nodes, tot_neles

    ! background nodes — no offset
    DO i = 1, bg_nodes
        READ(10,*) n, x, y, z, g, ib
        WRITE(20,*) n, x, y, z, g, ib
    END DO

    ! overset nodes — offset node index by bg_nodes
    DO i = 1, ov_nodes
        READ(11,*) n, x, y, z, g, ib
        WRITE(20,*) n + bg_nodes, x, y, z, g, ib
    END DO

    ! background elements — no offset
    DO i = 1, bg_neles
        READ(10,*) nel, n1, n2, n3, n4, c1, c2, c3, c4,ibc
        WRITE(20,*) nel, n1, n2, n3, n4, c1, c2, c3, c4
    END DO

    ! overset elements — offset nel and nod by bg_neles/bg_nodes, nc by bg_neles
    DO i = 1, ov_neles
        READ(11,*) nel, n1, n2, n3, n4, c1, c2, c3, c4
        WRITE(20,*) nel + bg_neles,        &
                    n1  + bg_nodes,        &
                    n2  + bg_nodes,        &
                    n3  + bg_nodes,        &
                    n4  + bg_nodes,        &
                    c1  + bg_neles,        &
                    c2  + bg_neles,        &
                    c3  + bg_neles,        &
                    c4  + bg_neles
    END DO

    CLOSE(20)

    ! =========================================================
    ! STEP 3 — write merged bc file
    ! =========================================================

    OPEN(21, file=trim(outdir)//'bc_proc0000.in', status='replace')

    WRITE(21,*) tot_nghosts

    ! background ghosts — offset nghost past tot_neles (bg_neles+ov_neles)
    ! so they don't collide with overset's real cells. nparent unchanged
    ! since bg real cell ids (1..bg_neles) are already correct.
    DO i = 1, bg_nghosts
        READ(12,*) np, ng, nt, ns
        WRITE(21,*) np, ng + ov_neles, nt, ns
    END DO

    ! overset ghosts — offset nparent by bg_neles, and offset nghost past
    ! tot_neles AND past all background ghosts, so they land right after
    ! the background ghost block ends (no collision).
    DO i = 1, ov_nghosts
        READ(13,*) np, ng, nt, ns
        WRITE(21,*) np + bg_neles, ng + bg_neles + bg_nghosts, nt, ns
    END DO

    CLOSE(21)

    ! =========================================================
    ! STEP 4 — write merged proc_bc file
    ! =========================================================

    OPEN(22, file=trim(outdir)//'proc_bc0000.in', status='replace')

    WRITE(22,*) tot_pghosts, tot_maxpg

    ! background pghosts — no offset
    DO i = 1, bg_pghosts
        READ(14,*) pp, d1, d2, d3, d4, d5, d6, dst
        WRITE(22,*) pp, d1, d2, d3, d4, d5, d6, dst
    END DO

    ! overset pghosts — offset pp_elem by bg_neles, dest_ghost by bg_neles+bg_nghosts
    DO i = 1, ov_pghosts
        READ(15,*) pp, d1, d2, d3, d4, d5, d6, dst
        WRITE(22,*) pp  + bg_neles,              &
                    d1, d2, d3, d4, d5, d6,      &
                    dst + bg_neles + bg_nghosts
    END DO

    CLOSE(22)

    CLOSE(10); CLOSE(11); CLOSE(12)
    CLOSE(13); CLOSE(14); CLOSE(15)

    PRINT*, 'Done. Merged files written to ', trim(outdir)

END PROGRAM merge_grids
