program householder
use houseHolderQR
IMPLICIT NONE

real, dimension(6,3)::A
real, dimension(6,1)::b
real, dimension(3,1)::x,var
real :: r
A(1,1) = -1
A(1,2) = 1
A(1,3) = 2
A(2,1) = -2
A(2,2) = 3
A(2,3) = 1
A(3,1) = 4
A(3,2) = 0
A(3,3) = 3
!A(4,1) = -1
!A(4,2) = 1
!A(4,3) = 0
!A(5,1) = -1
!A(5,2) = 0
!A(5,3) = 1
!A(6,1) = 0
!A(6,2) = -1
!A(6,3) = 1
b(1,1) = 1237
b(2,1) = 1941
b(3,1) = 2417
b(4,1) = 711
b(5,1) = 1177
b(6,1) = 475

!x(1,1) =0 
!x(2,1) = 0
!x(3,1) = 0

!call shhQR(A,b)
call variance(A,var)
call backSubstitution(A,b,x,r)

 
print*, var

end program householder
