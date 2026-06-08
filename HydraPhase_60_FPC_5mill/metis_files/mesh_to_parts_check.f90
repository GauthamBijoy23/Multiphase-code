program meshcheck
  implicit none
  integer :: i,j,n,nodes,neles,nel
  integer :: infile,ofile
  double precision, allocatable, dimension(:,:) :: c
  double precision, allocatable, dimension(:) :: nv
  integer, allocatable, dimension(:,:) :: node,neigh
  integer, allocatable,dimension(:) :: partition_node,part_element,num
  integer :: npartitions, ios, idx
  character(len=256) :: filename
  character(len=20)  :: numstr
  character(len=50) :: geometry_dir
  character(len=256) :: fname_gridfile, tmpfile
  character(len=256) :: metis_dir, fname_metis
   
character(len=256) :: ofile_name



  !!! mesh file reading
  
   geometry_dir = 'geometry_files/'
   metis_dir = 'metis_files/'
  
   write(fname_gridfile,'(A,"grid.dat")') trim(geometry_dir)
   infile = 52
   open(unit=infile, file=fname_gridfile, status='old', action='read')

!  infile=52
!  open(unit=infile,file='grid.dat',&
!       status='old',action='read')
  
  read(infile,*)nodes,neles

  allocate(c(nodes,3))
  allocate(node(neles,4))
  allocate(neigh(neles,4))
  allocate(nv(nodes))
  allocate(num(nodes))

  do i=1,nodes
    read(infile,*)n,c(n,1),c(n,2),c(n,3)
  end do
  
  do i=1,neles
    read(infile,*)nel,node(nel,1),node(nel,2),node(nel,3),&
         node(nel,4),&!,node(nel,5),node(nel,6),node(nel,7),&
         neigh(nel,1),neigh(nel,2),neigh(nel,3),&
         neigh(nel,4)!,neigh(nel,5),neigh(nel,6)
  end do

  close(infile)

  !!!read partition data
  allocate(partition_node(nodes))
  allocate(part_element(neles))
  
 
  ! Create a temporary file path
  tmpfile = trim(metis_dir)//'tmpfile'

! List npart files inside metis_dir
call execute_command_line('ls ' // trim(metis_dir) // 'mesh_metis.dat.npart.* > ' // trim(tmpfile), wait=.true.)

! Open the temporary file and read the first npart file
open(unit=99, file=tmpfile, status='old', action='read')
read(99,'(A)',iostat=ios) filename
 close(99)

! Remove temporary file
call execute_command_line('rm -f ' // trim(tmpfile))

! Check if read was successful
if (ios /= 0) then
    print *, 'Error: no mesh_metis.dat.npart.* file found in ', trim(metis_dir)
    stop
end if

 
  ! --- Extract the number after the last dot ---
  idx = index(filename, '.', back=.true.) + 1
  numstr = filename(idx:)
  read(numstr,*) npartitions

  print *, 'Detected partitions = ', npartitions

! --- Open and read npart file ---


call execute_command_line('ls ' // trim(metis_dir) // 'mesh_metis.dat.npart.* | xargs -n 1 basename > tmpfile', wait=.true.)

open(unit=99, file='tmpfile', status='old', action='read')
read(99,'(A)') filename
close(99)
call execute_command_line('rm -f tmpfile')

! Prepend metis_dir once
filename = trim(metis_dir)//trim(filename)
open(unit=infile, file=filename, status='old', action='read')


call execute_command_line('ls ' // trim(metis_dir) // 'mesh_metis.dat.epart.* | xargs -n 1 basename > tmpfile', wait=.true.)

open(unit=99, file='tmpfile', status='old', action='read')
read(99,'(A)') filename
close(99)
call execute_command_line('rm -f tmpfile')

! Prepend metis_dir once
filename = trim(metis_dir)//trim(filename)
open(unit=infile, file=filename, status='old', action='read')

  print *, 'Read partition data successfully.'

  call nodal_value()

  !!! output in tecplot format
    
!  ofile_name = trim(metis_dir)//'mesh_check.dat'
!  open(unit=ofile, file=ofile_name, status='replace', action='write')
    
!  ofile = 53
  !open(unit=ofile, file="mesh_check.dat", &
  !      status='replace', action="write")
!  write(ofile_name,*)'TITLE="mlater result from mesh_check"'
!  write(ofile_name,*)'VARIABLES = "X", "Y","Z","node","element"'
!  write(ofile_name,*)"ZONE N=",nodes, ' ,E=',neles,&
!    'DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'

 ! do i=1,nodes
 !   write(ofile_name,'(5(ES20.13,2X))')c(i,1),c(i,2),c(i,3),real(partition_node(i),8),nv(i)
 ! end do
  
 ! do i=1,neles
!    write(ofile_name,"(4(I8,2X))")node(i,1),node(i,2),node(i,3),node(i,4)!,&
!                             !node(i,7),node(i,8),node(i,4),node(i,3)
!  end do
!  close(ofile_name)
  
  
 
! Define folder and file name
ofile_name = trim(metis_dir)//'mesh_check.dat'

! Assign a unit number
ofile = 53

! Open file using unit number and filename
open(unit=ofile, file=ofile_name, status='replace', action='write')

! Write to the file using the **unit number**, not the filename
write(ofile,*) 'TITLE="mlater result from mesh_check"'
write(ofile,*) 'VARIABLES = "X", "Y","Z","node","element"'
write(ofile,*) "ZONE N=", nodes, ' ,E=', neles, &
               ' DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'

do i = 1, nodes
    write(ofile,'(5(ES20.13,2X))') c(i,1), c(i,2), c(i,3), real(partition_node(i),8), nv(i)
end do

do i = 1, neles
    write(ofile,'(4(I8,2X))') node(i,1), node(i,2), node(i,3), node(i,4)
end do

! Close the file
 close(ofile)

  
    
contains
  subroutine nodal_value
    implicit none
    integer :: i,j,k
    integer :: ns,npi
    integer :: n(4)
    do i=1,nodes
      nv(i)=0.0d0
      num(i)=0
    end do
    
    do i=1,neles
      do j=1,4
        k = node(i,j)
        nv(k) = nv(k)+part_element(i)
        num(k) = num(k)+1
      end do
    end do

    do i=1,nodes
      nv(i) = nv(i)/real(num(i),4)
    end do 
  end subroutine nodal_value


end program meshcheck
