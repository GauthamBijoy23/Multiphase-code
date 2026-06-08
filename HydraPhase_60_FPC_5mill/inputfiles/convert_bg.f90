program convert

implicit none

integer, parameter :: dp = kind(1.0d0)

integer :: N,E,i,dummy
integer :: n1,n2,n3,n4
integer :: nb1,nb2,nb3,nb4

real(dp) :: x,y,z

open(10,file='grid_bg.dat')
open(20,file='flow00000.dat')

read(10,*) N,E

write(20,'(A)') "X,Y,Z,Btag"
write(20,'(A)') "Grid to Tioga"
write(20,'("N=",I0," E=",I0)') N,E

do i=1,N
    read(10,*) dummy,x,y,z
    write(20,'(3(ES25.16,1X),I1)') x,y,z,1
end do

do i=1,E
    read(10,*) dummy,n1,n2,n3,n4,nb1,nb2,nb3,nb4
    write(20,'(4(I10,1X))') n1,n2,n3,n4
end do

write(20,'(I1)') 0
write(20,'(I1)') 0

close(10)
close(20)

end program convert
