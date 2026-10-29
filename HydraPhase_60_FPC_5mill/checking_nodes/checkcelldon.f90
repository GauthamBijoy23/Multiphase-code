program check_cell_donors
  implicit none
  character(len=*), parameter :: gridfile = 'grid_proc0000.dat'
  character(len=*), parameter :: prefile  = 'qcell_preinterp_0000.dat'
  character(len=*), parameter :: mapfile  = 'donorinfo_cell.dat'
  real(8), parameter :: tol = 1.d-8
  integer :: nn, ne, i, j, e, m, u, ios, r1, r2, r3, r4, off, rcell, dcell
  integer :: nrow, nout, nnf, nrf, nbad, gid, ibl
  integer, allocatable :: nod(:,:), ibc(:)
  real(8), allocatable :: x(:,:)
  real(8) :: w, xp(3), v(3,4), a1(3), a2(3), a3(3), bb(3), det, l(4), lmin, lworst
  character(len=512) :: line
integer, parameter :: neq = 10
real(8), allocatable :: qpre(:,:), qpost(:,:)
real(8) :: dq, dqmax
integer :: nq
character(len=*), parameter :: postfile = 'qcell_postinterp_0000.dat'

  open(newunit=u, file=gridfile, status='old')
  read(u,*) nn, ne
  allocate(x(3,nn), nod(4,ne), ibc(ne), qpre(neq,ne), qpost(neq,ne))
  do i = 1, nn
    read(u,*) j, x(1:3,j), gid, ibl
  end do
  do i = 1, ne
    read(u,*) j, nod(1:4,j)
  end do
  close(u)

open(newunit=u, file=prefile, status='old')
do i = 1, ne
  read(u,*) e, ibc(e), qpre(1:neq,e)
end do
close(u)

open(newunit=u, file=postfile, status='old')
do i = 1, ne
  read(u,*) e, ibl, qpost(1:neq,e)
end do
close(u)

nq = 0; dqmax = 0.d0
  nrow = 0; nout = 0; nnf = 0; nrf = 0; nbad = 0; lworst = 1.d0
  open(newunit=u, file=mapfile, status='old')
  do
    read(u,'(A)', iostat=ios) line
    if (ios /= 0) exit
    read(line,*, iostat=ios) r1, r2, r3, r4, off, w, rcell, dcell
    if (ios /= 0) cycle                       ! skips header lines
    nrow = nrow + 1
    if (rcell < 1 .or. rcell > ne .or. dcell < 1 .or. dcell > ne) then
      nbad = nbad + 1; cycle
    end if
    xp = 0.d0
    do m = 1, 4
      xp = xp + 0.25d0*x(:, nod(m,rcell))
      v(:,m) = x(:, nod(m,dcell))
    end do
    a1 = v(:,1)-v(:,4); a2 = v(:,2)-v(:,4); a3 = v(:,3)-v(:,4); bb = xp - v(:,4)
    det  = tp(a1,a2,a3)
    l(1) = tp(bb,a2,a3)/det
    l(2) = tp(a1,bb,a3)/det
    l(3) = tp(a1,a2,bb)/det
    l(4) = 1.d0 - l(1) - l(2) - l(3)
    lmin = minval(l); lworst = min(lworst, lmin)
    if (lmin < -tol) then
      nout = nout + 1
      if (nout <= 10) print*, 'outside: rec', rcell, 'don', dcell, 'lmin', lmin
    end if
    if (ibc(dcell) /= 1)  nnf = nnf + 1
    if (ibc(rcell) /= -1) nrf = nrf + 1
dq = maxval(abs(qpre(:,dcell) - qpost(:,rcell)) / (1.d0 + abs(qpre(:,dcell))))
dqmax = max(dqmax, dq)
if (dq > 1.d-12) then
  nq = nq + 1
  if (nq <= 10) print*, 'q mismatch: rec', rcell, 'don', dcell, 'max rel diff', dq
end if
  end do
  close(u)

  print*, 'rows checked          :', nrow
  print*, 'ids out of range      :', nbad
  print*, 'receptors outside donor:', nout
  print*, 'worst lambda          :', lworst
  print*, 'non-field donors      :', nnf
  print*, 'non-fringe receptors  :', nrf
  print*, 'q(1:10) donor /= receptor:', nq
  print*, 'max rel q difference     :', dqmax
contains
  real(8) function tp(a,b,c)
    real(8) :: a(3), b(3), c(3)
    tp = a(1)*(b(2)*c(3)-b(3)*c(2)) - a(2)*(b(1)*c(3)-b(3)*c(1)) + a(3)*(b(1)*c(2)-b(2)*c(1))
  end function
end program
