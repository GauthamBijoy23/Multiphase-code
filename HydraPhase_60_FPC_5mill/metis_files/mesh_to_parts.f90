program main
  implicit none
  integer :: neles,nodes,nghosts
  integer :: infile,ofile
  integer :: i,j,n,nel,ele_type,ncon
  double precision,allocatable,dimension(:,:) :: c
  integer,allocatable,dimension(:,:) :: node,neigh
  character(len=50) :: geometry_dir
  character(len=256) :: fname_gridfile
  character(len=256) :: metis_dir, fname_metis
  logical :: dir_exists
  
  metis_dir = 'metis_files/'

  ! Create folder if it does not exist
  inquire(file=metis_dir, exist=dir_exists)
  if (.not. dir_exists) call system("mkdir -p " // trim(metis_dir))
  
   geometry_dir = 'geometry_files/'
  
   write(fname_gridfile,'(A,"grid.dat")') trim(geometry_dir)
   infile = 52
   open(unit=infile, file=fname_gridfile, status='old', action='read')
  
  read(infile,*)nodes,neles

  allocate(c(nodes,3))
  allocate(node(neles,4))
  allocate(neigh(neles,4))

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

  !!!!metis mesh output!!!!
  ele_type = 0
  ncon=0
  write(fname_metis,'(A,"mesh_metis.dat")') trim(metis_dir)
  ofile = 53
  open(unit=ofile, file=fname_metis, action='write')

 ! open(unit=ofile,file='mesh_metis.dat',action='write')
  
  write(ofile,'((I8,2X))')neles

  do i=1,neles
    write(ofile,'(4(I8,2X))')node(i,1),node(i,2),node(i,3),&
                node(i,4)!,node(i,3),node(i,7),node(i,8),node(i,4)
  end do

  close(ofile)

end program main
