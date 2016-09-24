module houseHolderQR
implicit none

contains
! compute the vector norm2 of a vector
real function vectorNorm2(v)
  integer, parameter :: dk = selected_real_kind(15, 307)
  real (kind = dk), dimension(:), intent(in):: v
  vectorNorm2 = sqrt(sum(v**2)) ! compute vector norm2 using euclidean distance formula
end function vectorNorm2

! transform a matrix using House Holder transformation
subroutine shhQR(A,b)
  integer, parameter :: dk = selected_real_kind(15, 307)
  real (kind = dk), dimension(:,:), intent(out):: A, b ! output matrices, dimension need to be allocated before calling this subroutine
  real (kind = dk), dimension(:,:), allocatable:: e, u, v 
  integer:: m, n, k, i, j
  real (kind = dk):: alpha, beta, gmma, bgmma, one = 1 

  m = size(A(:,1)) ! # of rows in matrix A
  n = size(A(1,:)) ! # of columns in matrix A
  allocate(e(m,1), u(m,1), v(m,1)) ! allocate dimensions of e, u, and v

  ! house holder transformation
  do k = 1, min(n, m-1) ! loop over columns
  alpha = -sign(one,A(k,k))*vectorNorm2(A(k:,k)) ! compute a scalar multiple that preserves the norm
  do i = 1,m ! create basis vectors with 1 at i=k, and 0 elsewhere
    if (i==k) then
      e(i,1) = 1
    else
      e(i,1) = 0
    end if

    if (i>=k) then ! create vector u with 0 for i<k and A(i,k) for i>=k
      u(i,1) = A(i,k)
    else
      u(i,1) = 0
    end if
  end do

  v(:,1) = u(:,1) - alpha*e(:,1) ! compute v for the present column
  beta = dot_product(v(:,1),v(:,1)) ! compute transpose(v)*v

  if (beta == 0) then ! continue with next k with beta = 0
    exit
  end if

  do j = k,n ! loop over the remaining submatrix
    gmma = dot_product(v(:,1),A(:,j))
    A(:,j) = A(:,j) - (2*gmma/beta)*v(:,1) ! householder transformation for each element of the submatrix
  end do

  bgmma = dot_product(v(:,1),b(:,1)) 
  b(:,1) = b(:,1) - (2*bgmma/beta)*v(:,1) ! householder transformation for the right hand side matrix
  end do
  return
end subroutine shhQR

! subroutine to solve Ax = b (A being upper tringular) by backward substitution
subroutine backSubstitution(A,b,x,rSq)
  integer, parameter :: dk = selected_real_kind(15, 307)
  real (kind = dk), dimension(:,:), intent(in):: A ! input upper triangular matrix A
  real (kind = dk), dimension(:,:), intent(out)::b ! input right hand sand matrix b
  real (kind = dk), dimension(:,:), intent(out)::x ! output solution matrix x
  real (kind = dk):: rSq 
  integer:: i, m, n, j

  m = size(A(:,1)) ! # of rows in A
  n = size(A(1,:)) ! # of columns in A
  rSq = vectorNorm2(b(n+1:,1))**2 ! residual squared 

  do j = n, 1, -1 ! loop over rows from bottom
    if (A(j,j) == 0) then ! continue with next n if the diagonal element is 0
      exit
    end if
    x(j,1) = b(j,1)/A(j,j) ! solve for jth x
    do i = 1, j-1 ! loop over remaing submatrix
      b(i,1) = b(i,1) - A(i,j)*x(j,1) ! update submatrix b
    end do
  end do
  return
end subroutine backSubstitution

end module houseHolderQR
