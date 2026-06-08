program compare_active_cells

    implicit none

    integer :: N,E
    integer :: i, ios

    integer :: elem
    integer :: n1,n2,n3,n4
    integer :: nb1,nb2,nb3,nb4
    integer :: iblank

    integer :: active_elem
    integer :: active_count
    integer :: iblank_count

    character(len=512) :: line

    !----------------------------------
    ! Open files
    !----------------------------------

    open(10,file='grid_proc0000.dat',status='old')
    open(20,file='active_cells.dat',status='old')

    read(10,*) N,E

    !----------------------------------
    ! Skip node section
    !----------------------------------

    do i=1,N
        read(10,'(A)') line
    end do

    !----------------------------------
    ! Read first active element
    !----------------------------------

    read(20,*,iostat=ios) active_elem

    if (ios /= 0) then
        print *, "active_cells.dat empty"
        stop
    end if

    active_count = 0
    iblank_count = 0

    !----------------------------------
    ! Compare
    !----------------------------------

    do i=1,E

        read(10,*) elem,n1,n2,n3,n4, &
                    nb1,nb2,nb3,nb4,iblank

        if (iblank /= 1) cycle

        iblank_count = iblank_count + 1

        if (elem /= active_elem) then

            print *
            print *, "Mismatch found"
            print *, "Grid iblank=1 element :", elem
            print *, "Expected active cell  :", active_elem
            stop

        end if

        active_count = active_count + 1

        read(20,*,iostat=ios) active_elem

    end do

    !----------------------------------
    ! Check extra active cells
    !----------------------------------

    if (ios == 0) then
        print *
        print *, "Extra entries exist in active_cells.dat"
        print *, "First extra cell:", active_elem
        stop
    end if

    print *
    print *, "MATCH SUCCESSFUL"
    print *, "iblank=1 elements checked :", iblank_count
    print *, "active cells checked      :", active_count

    close(10)
    close(20)

end program compare_active_cells
