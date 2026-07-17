program diag_tool
implicit none
integer :: i,k,ios
real(8) :: q(10)
integer :: cnt

! ==== SUBROUTINE CALLS (uncomment/comment to toggle) ====
!call check_nan('qvals_postinterp.dat')
!call check_nan('q_node_tioga_check.dat')
!call check_nan('q_node.dat')
call check_cu('q_node_tioga_check.dat',0,0)
call check_cu('qvals_postinterp.dat',0,0)
call check_cu('q_node.dat',0,0)
call check_wsum('wsum.dat', 0, 0)   ! args: fname, write_nan(1/0), write_zero(1/0)
call check_cu('cu.dat', 0, 0)
call check_iblank_counts()
!call iblank_nodes([1,0,-1])   ! pass any subset e.g. [1], [0,-1], [1,0,-1]
!call check_overlap()
! ==========================================================

contains

subroutine check_nan(fname)
character(len=*), intent(in) :: fname
integer :: i,k,ios,cnt
real(8) :: q(10)
cnt=0
open(unit=50,file=fname,status='old')
do
  read(50,*,iostat=ios) i,(q(k),k=1,10)
  if (ios/=0) exit
  if (any(q/=q)) cnt=cnt+1
end do
close(50)
print*, 'nan in ',trim(fname),' = ',cnt
end subroutine check_nan

subroutine check_wsum(fname, write_nan, write_zero)
character(len=*), intent(in) :: fname
integer, intent(in) :: write_nan, write_zero
integer :: ios,nan_cnt,zero_cnt,n
real(8) :: w
character(len=128) :: nan_out, zero_out
nan_cnt=0; zero_cnt=0
n=0
nan_out = trim(fname(1:index(fname,'.',back=.true.)-1))//'_nan.dat'
zero_out = trim(fname(1:index(fname,'.',back=.true.)-1))//'_zeros.dat'
if (write_nan==1) open(unit=61,file=trim(nan_out),status='replace')
if (write_zero==1) open(unit=62,file=trim(zero_out),status='replace')
open(unit=51,file=fname,status='old')
do
  n=n+1
  read(51,*,iostat=ios) w
  if (ios/=0) exit
  if (w/=w) then
    nan_cnt=nan_cnt+1
    if (write_nan==1) write(61,*) n
  end if
  if (w==0.0) then
    zero_cnt=zero_cnt+1
    if (write_zero==1) write(62,*) n
  end if
end do
close(51)
if (write_nan==1) close(61)
if (write_zero==1) close(62)
print*, 'nan in ',trim(fname),' = ',nan_cnt
print*, 'zero in ',trim(fname),' = ',zero_cnt
end subroutine check_wsum

subroutine check_cu(fname, write_nan, write_zero)
character(len=*), intent(in) :: fname
integer, intent(in) :: write_nan, write_zero
integer :: k,ios,nan_cnt,zero_cnt,n
real(8) :: q(10)
character(len=128) :: nan_out, zero_out
nan_cnt=0; zero_cnt=0
n=0
nan_out = trim(fname(1:index(fname,'.',back=.true.)-1))//'_nan.dat'
zero_out = trim(fname(1:index(fname,'.',back=.true.)-1))//'_zeros.dat'
if (write_nan==1) open(unit=63,file=trim(nan_out),status='replace')
if (write_zero==1) open(unit=64,file=trim(zero_out),status='replace')
open(unit=52,file=fname,status='old')
do
  n=n+1
  read(52,*,iostat=ios) (q(k),k=1,10)
  if (ios/=0) exit
  if (any(q/=q)) then
    nan_cnt=nan_cnt+1
    if (write_nan==1) write(63,*) n
  end if
  if (any(q==0.0)) then
    zero_cnt=zero_cnt+1
    if (write_zero==1) write(64,*) n
  end if
end do
close(52)
if (write_nan==1) close(63)
if (write_zero==1) close(64)
print*, 'nan in ',trim(fname),' = ',nan_cnt
print*, 'zero in ',trim(fname),' = ',zero_cnt
end subroutine check_cu

subroutine iblank_nodes(modes)
integer, intent(in) :: modes(:)
integer :: N,E,i,n1,gid,ibn,cid,ib,ios,mm
integer :: v1,v2,v3,v4,c1,c2,c3,c4
integer, allocatable :: nod(:,:)
logical, allocatable :: used_field(:), used_hole(:), used_fringe(:)
logical :: want_field, want_hole, want_fringe
real(8) :: xx,yy,zz

want_field  = any(modes==1)
want_hole   = any(modes==0)
want_fringe = any(modes==-1)

open(unit=10,file='grid_proc0000.dat',status='old')
read(10,*) N,E
allocate(nod(E,4))
allocate(used_field(N), used_hole(N), used_fringe(N))
used_field=.false.; used_hole=.false.; used_fringe=.false.

do i=1,N
  read(10,*) n1,xx,yy,zz,gid,ibn
end do

do i=1,E
  read(10,*) cid,v1,v2,v3,v4,c1,c2,c3,c4
  nod(cid,1)=v1; nod(cid,2)=v2; nod(cid,3)=v3; nod(cid,4)=v4
end do
close(10)

open(unit=11,file='iblankcellslist.dat',status='old')
do
  read(11,*,iostat=ios) cid,ib
  if (ios/=0) exit
  if (ib==1 .and. want_field) then
    used_field(nod(cid,1))=.true.; used_field(nod(cid,2))=.true.
    used_field(nod(cid,3))=.true.; used_field(nod(cid,4))=.true.
  else if (ib==0 .and. want_hole) then
    used_hole(nod(cid,1))=.true.; used_hole(nod(cid,2))=.true.
    used_hole(nod(cid,3))=.true.; used_hole(nod(cid,4))=.true.
  else if (ib==-1 .and. want_fringe) then
    used_fringe(nod(cid,1))=.true.; used_fringe(nod(cid,2))=.true.
    used_fringe(nod(cid,3))=.true.; used_fringe(nod(cid,4))=.true.
  end if
end do
close(11)

if (want_field) then
  open(unit=12,file='iblanknodeslist.dat',status='replace')
  do i=1,N
    if (used_field(i)) write(12,*) i
  end do
  close(12)
end if

if (want_hole) then
  open(unit=13,file='holenodes.dat',status='replace')
  do i=1,N
    if (used_hole(i)) write(13,*) i
  end do
  close(13)
end if

if (want_fringe) then
  open(unit=14,file='fringenodes.dat',status='replace')
  do i=1,N
    if (used_fringe(i)) write(14,*) i
  end do
  close(14)
end if

deallocate(nod, used_field, used_hole, used_fringe)
end subroutine iblank_nodes

subroutine check_overlap()
integer :: N,i,ios,n1
logical, allocatable :: field(:), hole(:), fringe(:)
integer :: fh,ff,hf,fhf,total_unique

N = 200000   ! set to max possible node id, adjust if needed
allocate(field(N), hole(N), fringe(N))
field=.false.; hole=.false.; fringe=.false.

open(unit=20,file='iblanknodeslist.dat',status='old')
do
  read(20,*,iostat=ios) n1
  if (ios/=0) exit
  field(n1)=.true.
end do
close(20)

open(unit=21,file='holenodes.dat',status='old')
do
  read(21,*,iostat=ios) n1
  if (ios/=0) exit
  hole(n1)=.true.
end do
close(21)

open(unit=22,file='fringenodes.dat',status='old')
do
  read(22,*,iostat=ios) n1
  if (ios/=0) exit
  fringe(n1)=.true.
end do
close(22)

fh=0; ff=0; hf=0; fhf=0; total_unique=0
do i=1,N
  if (field(i) .and. hole(i)) fh=fh+1
  if (field(i) .and. fringe(i)) ff=ff+1
  if (hole(i) .and. fringe(i)) hf=hf+1
  if (field(i) .and. hole(i) .and. fringe(i)) fhf=fhf+1
  if (field(i) .or. hole(i) .or. fringe(i)) total_unique=total_unique+1
end do

print*, 'field-hole overlap    = ', fh
print*, 'field-fringe overlap  = ', ff
print*, 'hole-fringe overlap   = ', hf
print*, 'all three overlap     = ', fhf
print*, 'total unique nodes    = ', total_unique

deallocate(field, hole, fringe)
end subroutine check_overlap

subroutine check_iblank_counts()
integer :: cid,ib,ios,field_cnt,hole_cnt,fringe_cnt,other_cnt
field_cnt=0; hole_cnt=0; fringe_cnt=0; other_cnt=0
open(unit=30,file='iblankcellslist.dat',status='old')
do
  read(30,*,iostat=ios) cid,ib
  if (ios/=0) exit
  if (ib==1) then
    field_cnt=field_cnt+1
  else if (ib==0) then
    hole_cnt=hole_cnt+1
  else if (ib==-1) then
    fringe_cnt=fringe_cnt+1
  else
    other_cnt=other_cnt+1
  end if
end do
close(30)
print*, 'field cells (iblank=1)  = ', field_cnt
print*, 'hole cells (iblank=0)   = ', hole_cnt
print*, 'fringe cells (iblank=-1)= ', fringe_cnt
if (other_cnt>0) print*, 'other/unexpected iblank = ', other_cnt
end subroutine check_iblank_counts

end program diag_tool
