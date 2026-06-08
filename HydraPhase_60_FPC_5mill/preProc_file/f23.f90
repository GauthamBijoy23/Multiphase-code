!===========================================================
! Pre-processing code for tetrahedral mesh partitioning
!===========================================================
! Here I am writing in proc_bc 7th colum where writing cell id of the ghost cell of the neighbour processor 
! Now processor id will start from 0
Module metis_var
    ! Global arrays and variables used for mesh + METIS partitioning

    ! Node coordinates
!==============================================================
! Coordinates
!==============================================================
real, allocatable :: x(:), y(:), z(:)                ! Original node coordinates

!==============================================================
! Mesh topology: connectivity
!==============================================================
integer, allocatable :: nod(:,:)                     ! nod(e,1:4) = node IDs of element e
integer, allocatable :: nc1(:), nc2(:), nc3(:), nc4(:) ! Neighbor IDs across 4 faces of tetra

!==============================================================
! Ghost cell information (from METIS / domain decomposition)
!==============================================================
integer, allocatable :: nparent(:)                   ! Parent element index
integer, allocatable :: nghost(:)                    ! Ghost element index
integer, allocatable :: nside(:)                     ! Face/side index of ghost connection
integer, allocatable :: ntype(:)                     ! Type (-1 means from another processor)

!==============================================================
! Partitioning information
!==============================================================
integer, allocatable :: elem_count(:)                ! No. of elements per partition
integer, allocatable :: nod_count(:)                 ! No. of nodes per partition
integer, allocatable :: part(:)                      ! Part ID for each element

!==============================================================
! Processor local/global indexing
!==============================================================
integer, allocatable :: elem_g_to_l(:,:)             ! Element global→local mapping
integer, allocatable :: nod_g_to_l(:,:)              ! Node global→local mapping
integer, allocatable :: nod_present(:,:)             ! Flag: is node present in partition?
integer, allocatable :: proc_elem(:,:,:)             ! Element IDs in each processor
integer, allocatable :: proc_elem_l_to_g(:,:)        ! Local→global element IDs
integer, allocatable :: proc_nod(:,:)                ! Node IDs in each processor

!==============================================================
! Ghost tracking for each partition
!==============================================================
integer, allocatable :: ghost_count(:)               ! Number of ghost cells per partition
integer, allocatable :: ghost_local(:,:)             ! Ghost local indexing per partition
integer, allocatable :: np_local(:,:)                ! Local np mapping
integer, allocatable :: ng_local(:,:)                ! Local ng mapping
integer, allocatable :: ns_local(:,:)                ! Local ns mapping
integer, allocatable :: nt_local(:,:)                ! Local nt mapping

integer, allocatable :: gmap(:,:)                    ! gmap(local_elem, side) => local ghost id
integer, allocatable :: pghost_map(:,:,:)            ! Partition ghost map
integer, allocatable :: pghost_count(:)              ! Partition ghost count

integer, allocatable :: partghost_id(:,:,:)          ! Dimensions: (max_neigh, nfaces, npartitions)
integer, allocatable :: tmp(:,:)                     ! Temporary storage
integer, allocatable :: max_pghost(:)                ! Max ghost per partition
integer, allocatable :: column2(:), column8(:)

!==============================================================
! Scalars
!==============================================================
integer :: nodes, neles                              ! Total nodes, total elements
integer :: nghosts                                   ! Total number of ghost cells
integer :: numtotal                                  ! Total = neles + nghosts
integer :: npartitions                               ! Number of partitions
integer :: partition                                 ! Loop index for partition
integer :: blocksize

integer :: n_boundary, n_ghost, n_total              ! Local counters
integer :: e_loc, np, ng, g_loc
integer :: neigh_id, idx, nbndry
integer :: k, parent_local, side_k
integer :: m, m1, m2, m3, m4, n1, n2, n3, n4
integer :: e_global, neigh_global, neigh_local_in_neigh
integer :: maxghost
integer :: f_neigh
integer :: g, count, ppar, pgho, ttype, side
integer :: nt, neigh_proc
integer :: n, i, j, nel, ns, neigh, iface, nqu, blk_parent
integer :: e, ee, f, ff                             ! Element/face indices for loops
integer :: p, q
integer :: mmax_pghost
integer, allocatable :: idx_sort(:)
integer :: start, finish, this_proc

!==============================================================
! File handling
!==============================================================
character(len=256)  :: fname,mname
character(len=256)  :: filename
character(len=20)   :: numstr
integer :: ios, i1

integer, allocatable :: proc(:,:),idxx(:),tmpl(:,:)
integer :: col3_value

End Module

!===========================================================
Program metis_check
    Use metis_var
    Implicit none

    !===========================================
    ! Step 1: Read basic mesh and ghost info
    !===========================================

    ! Open files from the folder 'geometry_files'
    open(unit=10, file='geometry_files/grid.dat', status='old', action='read')
    open(unit=11, file='geometry_files/bc.in', status='old', action='read')


    read(10,*) nodes, neles          ! Total nodes and total elements
    read(11,*) nghosts               ! Number of ghost cells
    numtotal = neles + nghosts       ! Total cells including ghosts

    print*, 'Total number of ghost and elements Numtotal = ', numtotal

    ! Allocate memory for mesh data
    Allocate (x(nodes), y(nodes), z(nodes))                ! Node coordinates

    Allocate (nod(numtotal,4))                             ! Element node connectivity
    Allocate (nc1(numtotal), nc2(numtotal), nc3(numtotal), nc4(numtotal)) ! Neighbor indices

    ! Allocate ghost cell info arrays
    Allocate (nparent(nghosts), nghost(nghosts), ntype(nghosts), nside(nghosts))

    ! Allocate partition info
    Allocate (part(numtotal))


    !=============================
    ! Step 2: Read mesh data
    !=============================

    ! Read node coordinates
    Do i = 1, nodes
        read(10,*) n, x(n), y(n), z(n)
    End Do

    ! Read element connectivity and neighbors
    Do i = 1, neles
        read(10,*) nel, nod(nel,1), nod(nel,2), nod(nel,3), nod(nel,4), &
                   nc1(nel), nc2(nel), nc3(nel), nc4(nel)
    End Do


    ! Read ghost cell info (boundary conditions from domain decomposition)
    Do i = 1, nghosts
        read(11,*) nparent(i), nghost(i), ntype(i), nside(i)
    End Do

    close(10)
    close(11)

    !=================================================
    ! Step 3: Read partition information (from METIS)
    !=================================================

    ! --- Find the actual file name ---
    ! (Here I assume only one such file exists in the directory)
    call execute_command_line('ls metis_files/mesh_metis.dat.epart.* > tmpfile', wait=.true.)
    open(unit=99, file='tmpfile', status='old', action='read')
    read(99,'(A)',iostat=ios) filename
    close(99)
    call execute_command_line('rm -f tmpfile')
    
     if (ios /= 0) then
        print *, 'Error: no partition file found.'
        stop
    end if
    
    ! --- Extract the number after the last dot ---
    i1 = index(filename, '.', back=.true.) + 1
    numstr = filename(i1:)
    read(numstr,*) npartitions

    print *, 'Detected number of partitions = ', npartitions
    print *, 'Opening file: ', trim(filename)

    !npartitions = 3
    Allocate (elem_count(npartitions), nod_count(npartitions))
    Allocate (nod_g_to_l(nodes,npartitions))
    Allocate (elem_g_to_l(neles,npartitions))
    Allocate (nod_present(nodes,npartitions))
    Allocate (proc_elem_l_to_g(neles,npartitions))
    Allocate(max_pghost(npartitions))

    elem_count = 0

     ! --- Now open the correct file ---
    open(unit=12, file=trim(filename), status='old', action='read')
!    open(unit=12,file='mesh_metis.dat.epart.3')
    Do i = 1, neles
        Read(12,*) part(i)      ! Read partition ID from METIS output
        Do partition = 1, npartitions
            If (part(i) == partition-1) Then
                elem_count(partition) = elem_count(partition) + 1
                elem_g_to_l(i,partition) = elem_count(partition)
                proc_elem_l_to_g(elem_count(partition), partition) = i
            End If
        End Do
    End Do
    close(12)
    print*, proc_elem_l_to_g(46,1)

! Print maximum number of elements in a partition
!    print*,"Max Element",maxval(elem_count)
    print*, " Element count", elem_count(:) 

    Allocate (proc_elem(4, maxval(elem_count), npartitions)) !4 for tetra
    Allocate (proc_nod(nodes, npartitions))

    ! Finding nodes in each block
    elem_count = 0
    nod_present = 0
    Do i = 1, neles
        blk_parent = part(i) + 1
        elem_count(blk_parent) = elem_count(blk_parent) + 1
        Do j = 1, 4
            nqu = nod(i,j)
            proc_elem(j, elem_count(blk_parent), blk_parent) = nqu
            nod_present(nqu, blk_parent) = 1
        End Do
    End Do

    ! Global to local mapping for nodes
    nod_count = 0
    Do partition = 1, npartitions
        Do i = 1, nodes
            If (nod_present(i, partition) == 1) Then
                nod_count(partition) = nod_count(partition) + 1
                nod_g_to_l(i, partition) = nod_count(partition)
                proc_nod(nod_count(partition), partition) = i
            End If
        End Do
    End Do
!    print*,maxval(nod_count)

    !--------Establishing global to local link of nodes
        nod_count = 0
        Do partition = 1, npartitions
                Do i = 1,nodes
                        If(nod_present(i,partition) .eq. 1) Then
                                nod_count(partition) = nod_count(partition) + 1
                                nod_g_to_l(i,partition) = nod_count(partition)
                        Endif
                Enddo
        Enddo


    !=============================
    ! Step 4: Process ghost cells
    !=============================

    Allocate(ghost_count(npartitions))
    Allocate(ghost_local(maxval(elem_count), npartitions))
    Allocate(np_local(nghosts,npartitions), ng_local(nghosts,npartitions), ns_local(nghosts,npartitions))
    Allocate(nt_local(nghosts,npartitions))    ! <-- add ntype storage

    ghost_count = 0
    Do partition = 1, npartitions
        Do i = 1, nghosts
            np = nparent(i)
            ng = nghost(i)
            ns = nside(i)
            blk_parent = part(np)+1
            If (blk_parent == partition) Then
                ghost_count(partition) = ghost_count(partition) + 1
                ghost_local(elem_g_to_l(np,partition), partition) = elem_count(partition) + ghost_count(partition)
                np_local(ghost_count(partition),partition) = elem_g_to_l(np,partition)
                ng_local(ghost_count(partition),partition) =  ghost_local(elem_g_to_l(np,partition), partition)
                ns_local(ghost_count(partition),partition) = ns
                nt_local(ghost_count(partition),partition) = ntype(i)   ! <-- store ntype

                !if (partition.eq.1) print*, np, i, ns, ntype(i)
            End If
        End Do
    End Do
        !============================
    ! Step 5: Write bc_proc files
    !============================

    do partition = 1, npartitions
 !    print*, "ghost_count(", partition, ")" , ghost_count(partition)
    enddo

    Do partition = 1, npartitions
        write(fname, '(A,I4.4,A)') 'geometry_files/bc_proc', partition-1, '.in'
     !   write(fname,'(A,I0,A,I4.4,A)')  'geometry_files/bc_proc', partition-1, '.in'

!        write(fname,'(A,I4.4,A)') 'bc_proc', partition-1, '.in'
        open(unit=11,file=fname,status='replace',action='write')
        write(11,*) ghost_count(partition)

        Do i = 1, ghost_count(partition)
            write(11,'(4i10)') np_local(i,partition), ng_local(i,partition), nt_local(i,partition), ns_local(i,partition)
        End Do
        close(11)
!        print*, 'Partition', partition, 'ghosts written =', ghost_count(partition)
    End Do

!=============================================
!Writing glob to loc mapping 
!============================================
!=================================================
! Write local→global mapping for each partition
!=================================================


mname = 'geometry_files/map_all.in'
open(unit=21, file=mname, status='replace', action='write')
!open(unit=21, file='geometry_files/map_all.in', status='replace', action='write')

! Header
write(21,'(A)', advance='no') 'LocalID'
Do partition = 1, npartitions
    write(21,'(A,I4.4)', advance='no') '   GlobalId_Partition', partition-1
End Do
write(21,*)

! Loop over local IDs (assuming local IDs run up to max elem_count)
Do i = 1, maxval(elem_count)
    write(21,'(I10)', advance='no') i
    Do partition = 1, npartitions
        if (i <= elem_count(partition)) then
            write(21,'(I20)', advance='no') proc_elem_l_to_g(i, partition)
        else
            ! If this partition has fewer elements, just leave blank (or write 0)
            write(21,'(I20)', advance='no') 0
        end if
    End Do
    write(21,*)
End Do

close(21)
    !===================================
    ! Writing grid.dat for each process
    !===================================
    ! maxlocal = elem_count(partition)
    Allocate(partghost_id(neles, 4, npartitions)) !##CHANGED THIS
    partghost_id(:,:,:) = 0
    Allocate(gmap(maxval(elem_count), 4)) !##CHANGED THIS
    Allocate(pghost_count(npartitions))
    Allocate(pghost_map(neles,npartitions,4))   ! map global element → local pghost id (0 if not used)
    pghost_map = 0
    pghost_count = 0
    gmap = 0

Do partition = 1, npartitions

    pghost_count(partition) = 0
    ! Also zero only the portion of pghost_map for this partition (optional but safe)
    pghost_map(:, partition, :) = 0
    ! Zero gmap up to elem_count(partition)
    gmap(1:elem_count(partition), :) = 0

    write(fname, '(A,I4.4,A)') 'geometry_files/grid_proc', partition-1, '.dat'
    open(unit=13,file=fname,status='replace',action='write')

    ! Header: (#local nodes, #local elements)
    write(13,*) nod_count(partition), elem_count(partition)

    ! -------- Nodes (local list only) --------
    ! proc_nod(ln,partition) holds the global node id for local node ln.
    Do i = 1, nod_count(partition)
        n = proc_nod(i, partition)      ! global node id
        write(13,'(i10,3e23.15)') i, x(n), y(n), z(n)
    End Do

    ! ghost_count(partition) entries store:
       !   np_local(k,partition) = parent_local (local index of parent)
       !   ns_local(k,partition) = side (1..4)
       !   ng_local(k,partition) = local ghost id (we filled this earlier)
       do k = 1, ghost_count(partition)
         parent_local = np_local(k, partition)      ! local index of parent element
           side_k       = ns_local(k, partition)     ! side index 1..4
         if (parent_local >= 1 .and. parent_local <= elem_count(partition) .and. &
            side_k >= 1 .and. side_k <= 4) then
            gmap(parent_local, side_k) = ng_local(k, partition)
         end if
       end do

    ! -------- Elements (local) + local neighbor ids --------
    ! proc_elem_l_to_g(el,partition) gives global element id for local element el.
    Do i = 1, elem_count(partition)
        nel = proc_elem_l_to_g(i, partition)   ! global element id

        ! local node ids for this local element i (convert from global nodes in proc_elem)
        ! proc_elem(:,i,partition) stores the 4 GLOBAL node IDs of this element;
        ! nod_g_to_l maps them to local node IDs.
        n1 = nod_g_to_l( proc_elem(1,i,partition), partition )
        n2 = nod_g_to_l( proc_elem(2,i,partition), partition )
        n3 = nod_g_to_l( proc_elem(3,i,partition), partition )
        n4 = nod_g_to_l( proc_elem(4,i,partition), partition )


        ! neighbor local ids; 0 means outside this partition (either other proc or physical boundary)
        ! Face 1
        neigh = nc1(nel)
        if (neigh >= 1 .and. neigh <= neles) then
            if (part(neigh)+1 == partition) then
                m1 = elem_g_to_l(neigh, partition)
            else
                    if (pghost_map(neigh, partition, 1) == 0) then
                       pghost_count(partition) = pghost_count(partition) + 1
                    pghost_map(neigh,partition,1) = elem_count(partition) + ghost_count(partition) + pghost_count(partition)
                    end if
                    m1=pghost_map(neigh,partition,1)
 !                    print*, "neigh" , neigh
                    partghost_id(neigh, 1, partition) = m1
             end if
        else
            m1 =  gmap(i,1)
        end if

        ! Face 2
        neigh = nc2(nel)
        if (neigh >= 1 .and. neigh <= neles) then
            if (part(neigh)+1 == partition) then
                m2 = elem_g_to_l(neigh, partition)
            else
                    if (pghost_map(neigh,partition,2) == 0) then
                       pghost_count(partition) = pghost_count(partition) + 1
                    pghost_map(neigh,partition,2) = elem_count(partition) + ghost_count(partition) + pghost_count(partition)
                    end if
                    m2=pghost_map(neigh,partition,2)
                    partghost_id(neigh, 2, partition) = m2

            end if
        else
            m2 = gmap(i,2)
        end if

        ! Face 3
        neigh = nc3(nel)
        if (neigh >= 1 .and. neigh <= neles) then
            if (part(neigh)+1 == partition) then
                m3 = elem_g_to_l(neigh, partition)
            else
                    if (pghost_map(neigh,partition,3) == 0) then
                       pghost_count(partition) = pghost_count(partition) + 1
                    pghost_map(neigh,partition,3) = elem_count(partition) + ghost_count(partition) + pghost_count(partition)
                    end if
                    m3=pghost_map(neigh,partition,3)
                    partghost_id(neigh, 3, partition) = m3

            end if
        else
            m3 = gmap(i,3)
        end if

        ! Face 4
        neigh = nc4(nel)
        if (neigh >= 1 .and. neigh <= neles) then
            if (part(neigh)+1 == partition) then
                m4 = elem_g_to_l(neigh, partition)
            else
                    if (pghost_map(neigh,partition,4) == 0) then
                       pghost_count(partition) = pghost_count(partition) + 1
                    pghost_map(neigh,partition,4) = elem_count(partition) + ghost_count(partition) + pghost_count(partition)
                    end if
                    m4=pghost_map(neigh,partition,4)
                    partghost_id(neigh, 4, partition) = m4

            end if
        else
            m4 = gmap(i,4)
        !    print*, "gmap", i, " ", gmap(i,4)
        end if

        ! Write: local elem id, 4 local node ids, then 4 local neighbor ids
        write(13,'(9i10)') i, n1, n2, n3, n4, m1, m2, m3, m4
    End Do
    End do

!===========================================================================================================================
! ----------------------------
! Write proc_bc files (clean)
! ----------------------------

!calculate maximum p ghosts
do partition = 1, npartitions

  ! count boundary faces for this partition (only once)
  nbndry = 0
  do e_loc = 1, elem_count(partition)
    e_global = proc_elem_l_to_g(e_loc, partition)
    do iface = 1, 4
      select case(iface)
      case(1); neigh_global = nc1(e_global)
      case(2); neigh_global = nc2(e_global)
      case(3); neigh_global = nc3(e_global)
      case(4); neigh_global = nc4(e_global)
      end select
      if (neigh_global >= 1 .and. neigh_global <= neles) then
        if (part(neigh_global)+1 /= partition) then
          nbndry = nbndry + 1
        end if
      end if
    end do
  end do
           max_pghost(partition) = nbndry 
end do  

mmax_pghost = maxval(max_pghost(:))


do partition = 1, npartitions

  ! count boundary faces for this partition (only once)
  nbndry = 0
  do e_loc = 1, elem_count(partition)
    e_global = proc_elem_l_to_g(e_loc, partition)
    do iface = 1, 4
      select case(iface)
      case(1); neigh_global = nc1(e_global)
      case(2); neigh_global = nc2(e_global)
      case(3); neigh_global = nc3(e_global)
      case(4); neigh_global = nc4(e_global)
      end select
      if (neigh_global >= 1 .and. neigh_global <= neles) then
        if (part(neigh_global)+1 /= partition) then
          nbndry = nbndry + 1
        end if
      end if
    end do
  end do

  if (nbndry == 0) then
    ! Still write a small file with zero rows (keeps tools simple)
    write(fname,'("geometry_files/proc_bc",I4.4,".in")') partition-1
    open(unit=20, file=fname, status='replace', action='write')
    write(20,'(2I10)') 0, 0
    close(20)
    cycle
  end if

  ! allocate the tmp array: 8 rows x nbndry columns
  allocate(tmp(8, nbndry))
  tmp = 0
  idx = 0

  ! fill rows
  do e_loc = 1, elem_count(partition)
    e_global = proc_elem_l_to_g(e_loc, partition)

    do iface = 1, 4
      select case(iface)
      case(1); neigh_global = nc1(e_global)
      case(2); neigh_global = nc2(e_global)
      case(3); neigh_global = nc3(e_global)
      case(4); neigh_global = nc4(e_global)
      end select

      if (neigh_global >= 1 .and. neigh_global <= neles) then
        neigh_proc = part(neigh_global) + 1
        if (neigh_proc /= partition) then
          ! We have a boundary to another partition -> output row
          idx = idx + 1

          ! column 1: parent local id in this partition
          tmp(1, idx) = e_loc

          ! column 2: neighbor's local id in neighbor partition
          neigh_local_in_neigh = elem_g_to_l(neigh_global, neigh_proc)
          tmp(2, idx) = neigh_local_in_neigh

          ! column 3: neighbor proc (0-based)
          tmp(3, idx) = neigh_proc - 1

          ! column 4: face index on parent (1..4)
          tmp(4, idx) = iface

          ! column 5: local ghost id in THIS partition that represents 'neigh' (where THIS partition will store received values)
          ! This must have been set when writing grid_proc: partghost_id(neigh_global, face, partition)
          tmp(5, idx) = partghost_id(neigh_global, iface, partition)

          ! column 6: neighbor global id (for debug)
          tmp(6, idx) = neigh_global

          ! Now compute f_neigh = the face index in neigh_global that references e_global
          f_neigh = 0
          if (nc1(neigh_global) == e_global) then
             f_neigh = 1
          else if (nc2(neigh_global) == e_global) then
             f_neigh = 2
          else if (nc3(neigh_global) == e_global) then
             f_neigh = 3
          else if (nc4(neigh_global) == e_global) then
             f_neigh = 4
          else
             ! fallback: if neighbor's nc* does not directly reference e_global, try node-triplet match
             ! This is rarely needed but robust.
             ! Build neighbor face-node triplets and compare with the parent face nodes
             ! (Implement if you have node connectivity available. For now, set f_neigh=1 and warn.)
             f_neigh = 1
             print*, 'WARNING: could not find direct reciprocal face (neigh=', neigh_global, 'parent=', e_global, ')&
                     - using f_neigh=1'
          end if

          ! column 7: the ghost id in neighbor partition where THEY will receive parent value
          ! This must match the local ghost printed in neighbor's grid_*.dat
          tmp(7, idx) = partghost_id(e_global, f_neigh, neigh_proc)

          ! column 8: where THIS partition will receive data FROM neighbor (local ghost id) -> same as col5
          tmp(8, idx) = tmp(5, idx)

        end if
      end if
    end do
  end do


  ! final sanity check: idx must equal nbndry
  if (idx /= nbndry) then
     print*, 'Internal error: counted and filled nbndry mismatch in partition', partition, ' counted=', nbndry, ' filled=', idx
  end if

  

   !REORDERING COLUMN 8 TO RECIEVE IN SAME ORDER AS SEND ORDER
  ! Now reorder col(8) within each proc group by ascending col(2)
start = 1
do while (start <= nbndry)
   this_proc = tmp(3, start)

   ! find contiguous block for this_proc
   finish = start
   do while (finish < nbndry .and. tmp(3, finish+1) == this_proc)
      finish = finish + 1
   end do

   blocksize = finish - start + 1

   if (blocksize > 1) then
      ! Make a copy of col(2) and col(8) for this block
      allocate(column2(blocksize), column8(blocksize))
      do i = 1, blocksize
         column2(i) = tmp(2, start+i-1)
         column8(i) = tmp(8, start+i-1)
      end do

      ! Sort indices of block2 ascending
      allocate(idx_sort(blocksize))
      do i = 1, blocksize
         idx_sort(i) = i
      end do
     call sort_indices_by_values(column2, idx_sort, blocksize)

      ! Reorder col(8) according to sorted col(2)
      do i = 1, blocksize
         tmp(8, start+i-1) = column8(idx_sort(i))
      end do

      deallocate(column2, column8, idx_sort)
   end if

   start = finish + 1
end do

! Optionally sort by neighbor proc (column 3) for better grouping
  ! You can call your existing sort_by_column(tmp, 3, nbndry) if available.
   call sort_by_column(tmp, 3, nbndry)
   
!--------reorder done----------------
  ! Write proc_bc file
  write(fname,'("geometry_files/proc_bc",I4.4,".in")') partition-1
  open(unit=20, file=fname, status='replace', action='write')
  ! write header: number of rows and max_pghost (optional second int)
  ! We'll compute max_pghost as the maximum local ghost id used within this partition
  maxghost = 0
  do i = 1, nbndry
    if (tmp(5,i) > maxghost) maxghost = tmp(5,i)
    if (tmp(7,i) > maxghost) maxghost = max(maxghost, tmp(7,i))
  end do
  write(20,'(2I10)') nbndry, mmax_pghost

  do i = 1, nbndry
    write(20,'(8I10)') tmp(1,i), tmp(2,i), tmp(3,i), tmp(4,i), tmp(5,i), tmp(6,i), tmp(7,i), tmp(8,i)
  end do
  close(20)

  deallocate(tmp)

end do  ! partition loop

!!--------Pain loop-----------------------
 
allocate(proc(mmax_pghost, npartitions))
allocate(idxx(npartitions))
idxx(:) = 0
! Loop over proc_bc files
do partition = 1, npartitions
    write(fname,'("geometry_files/proc_bc",I4.4,".in")') partition-1
    open(unit=10, file=fname, status='old', action='read')

    ! Read first line: bndry_total, max_pghost
    read(10,*) nbndry, mmax_pghost

    ! Allocate temp array to read all 7 columns
    allocate(tmpl(8,nbndry))

    ! Read all entries
    do i = 1, nbndry
        read(10,*) tmpl(:,i)
    end do

    close(10)

    ! Loop through entries and separate based on 3rd column
    do i = 1, nbndry
       col3_value = tmpl(3,i) + 1   ! processor id +1 because Fortran arrays start at 1
      idxx(col3_value) = idxx(col3_value) + 1
      proc(idxx(col3_value), col3_value) = tmpl(7,i)   
    end do
   

 !if (partition.eq.1)   print*,"proc0",tmpl(8,:)
  ! print*,"proc1",proc(:,1)
    
   
    deallocate(tmpl)
end do


!==============================
! Separate loop to add 8th column
!==============================
do partition = 1, npartitions
    write(fname,'("geometry_files/proc_bc",I4.4,".in")') partition-1
    open(unit=10, file=fname, status='old', action='read')

    ! Read first line: bndry_total, max_pghost
    read(10,*) nbndry, mmax_pghost

    ! Allocate temp array to read 7 columns + 1 for 8th column
    allocate(tmpl(8,nbndry))

    ! Read first 7 columns and initialize 8th to 0
    do i = 1, nbndry
        read(10,*) tmpl(1:7,i)
        tmpl(8,i) = 0
    end do

    close(10)

    ! Fill 8th column using previously filled proc array
    do i = 1, nbndry
        col3_value = tmpl(3,i) + 1           ! processor id +1 for Fortran indexing
        tmpl(8,i) = proc(i, partition)
    end do

    ! Write back updated proc_bc file with 8 columns
    open(unit=10, file=fname, status='replace', action='write')
    write(10,'(2I10)') nbndry, mmax_pghost
    do i = 1, nbndry
        write(10,'(8I10)') tmpl(:,i)
    end do
    close(10)

    deallocate(tmpl)
end do




    !=====================   
    ! Writing tplot files!
    !=====================
    Do partition = 1, npartitions
        write(fname, '(A,I4.4,A)') 'tecplot_files/tplot_fin_proc', partition-1, '.dat'
        open(unit=12,file=fname)
        write(12,*) 'VARIABLES = "X","Y","Z","Part"'
        write(12,*) 'ZONE F=FEPOINT,ET=TETRAHEDRON,N =', nod_count(partition), ',E =', elem_count(partition)

        Do i = 1, nodes
            If (nod_present(i, partition) == 1) &
                write(12,'(3e23.15,i10)') x(i), y(i), z(i), partition
        End Do

        Do i = 1, elem_count(partition)
            write(12,'(4i10)') &
                nod_g_to_l(proc_elem(1,i,partition),partition), &
                nod_g_to_l(proc_elem(2,i,partition),partition), &
                nod_g_to_l(proc_elem(3,i,partition),partition), &
                nod_g_to_l(proc_elem(4,i,partition),partition)
        End Do
        close(12)
    End Do


! ----------------------------
! End of proc_bc writer
! ----------------------------
contains
    subroutine sort_by_column(arr, col, n)
    integer, intent(in) :: col, n
    integer, intent(inout) :: arr(:,:)
    integer :: i, j, k
    integer :: tmp_row(size(arr,1))

    do i = 1, n-1
        k = i
        do j = i+1, n
            if (arr(col,j) < arr(col,k)) k = j
        end do
        if (k /= i) then
            tmp_row = arr(:,i)
            arr(:,i) = arr(:,k)
            arr(:,k) = tmp_row
        end if
    end do
end subroutine
subroutine sort_indices_by_values(arr, idx, n)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: arr(n)
  integer, intent(inout) :: idx(n)
  integer :: i, j, tmp_idx, tmp_val

  ! simple insertion sort on indices
  do i = 2, n
     tmp_idx = idx(i)
     tmp_val = arr(tmp_idx)
     j = i - 1
     do while (j >= 1 .and. arr(idx(j)) > tmp_val)
        idx(j+1) = idx(j)
        j = j - 1
     end do
     idx(j+1) = tmp_idx
  end do
end subroutine sort_indices_by_values


End Program
