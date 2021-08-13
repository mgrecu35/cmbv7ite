! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
subroutine minvert(n, A, Ainv)
  real, intent(in) :: A(n,n)
  real,intent(out) :: Ainv(n,n)

  real, dimension(n) :: work  ! work array for LAPACK
  integer, dimension(n) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external SGETRF
  external SGETRI
  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  

  ! SGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.

  call SGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! SGETRI computes the inverse of a matrix using the LU factorization
  ! computed by SGETRF.
  call SGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end subroutine minvert
