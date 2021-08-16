subroutine oe_step(nobs,nvar,x,xa,sa_i,y,F,sy_i,K,xnext,sx,sx_i,prod2,xmin,xmax)
  !This subroutine performs the Gauss-Newton step as described by Rodgers (2000)
  implicit none

  integer :: nobs, nvar, a, b
  real, intent(in) :: x(nvar), xa(nvar)
  real, intent(out):: xnext(nvar)
  real, intent(out):: sx(nvar,nvar), sx_i(nvar,nvar), prod2(nvar,nvar)
  real, intent(in), optional :: xmin(nvar), xmax(nvar)
  real, intent(in) :: y(nobs), F(nobs)
  real, intent(in) :: sy_i(nobs,nobs),sa_i(nvar,nvar)
  real, intent(in) :: K(nobs,nvar)

  real :: prod1(nvar,nobs),prod3(nvar),prod4(nvar),prod5(nvar), K_t(nvar,nobs)
  real :: sum1(nvar), sum2(nobs), sum3(nvar), sum4(nobs)
  real*8 :: t1, dt
  
  !t1=secnds(0.0)
  K_t = transpose(K)
  prod1 = matmul(K_t,sy_i)
  prod2 = matmul(prod1,K)

  sx_i = sa_i+prod2
!   if(nvar .gt. 10) then
!     do a=1,nvar
!       do b=1,nvar
!         if(isNaN(sa_i(a,b))) print*, 1, a, b
!         if(isNaN(sx_i(a,b))) print*, 2, a, b
!         if(isNaN(prod2(a,b))) print*, 3, a, b
!       end do
!     end do
!   endif
  call minvert(nvar,sx_i,sx)

  
  sum1 = xa-x
  
  sum2 = y-F
  prod3 = matmul(sa_i,sum1)
  prod4 = matmul(prod1,sum2)
  sum3=prod4+prod3
  prod5 = matmul(sx,sum3)
  xnext = x+prod5
  
  !dt=secnds(t1)
  !if(nvar .gt. 10) print*, dt
          
  !If limits are present, ensure that variables are within hard limits
  if(present(xmax)) then
    do a=1,nvar	  
      if(xnext(a) .gt. xmax(a)) xnext(a) = xmax(a)
    end do
  endif
  if(present(xmin)) then
    do a=1,nvar	  
      if (xnext(a) .lt. xmin(a)) xnext(a) = xmin(a)
    end do
  endif
  do a=1,nvar
     sx(a,a)=max(sx(a,a),1e-9)
  enddo
end subroutine oe_step

subroutine oe_step_mirs(nobs,nvar,x,K,Sa,Sy,y,F,dx,xnext,xmin,xmax)
  !MIRS version of OE step. should be equivalent of above, but faster when nvar >> nobs. Have not verified results.
  implicit none

  integer :: nobs, nvar, a, b
  real, intent(in) :: x(nvar), dx(nvar,1)
  real, intent(out):: xnext(nvar)
  real, intent(in), optional :: xmin(nvar), xmax(nvar)
  real, intent(in) :: y(nobs), F(nobs)
  real, intent(in) :: sy(nobs,nobs),sa(nvar,nvar)
  real, intent(in) :: K(nobs,nvar)
  real :: K_t(nvar,nobs)

  
  real :: prod1(nobs,nvar), prod2(nobs,nobs), prod3(nobs,nobs), prod4(nvar,nobs), prod5(nvar,nobs), prod6(nobs,1)
  real :: sum1(nobs,nobs), sum2(nobs), sum3(nobs)
  real :: xstep(nvar)
  real*8 :: t1, dt
  
  t1=secnds(0.0)
  K_t=transpose(K)
  !dt=secnds(t1)
  !print*, dt
  
  !t1=secnds(0.0)
  prod1=matmul(K,Sa)
  !dt=secnds(t1)
  !print*, dt
  
  !t1=secnds(0.0)
  prod2=matmul(prod1,K_t)
  !dt=secnds(t1)
  !print*, dt
  
  !t1=secnds(0.0)
  sum1=prod2+Sy
  !dt=secnds(t1)
  !print*, dt
  
  !t1=secnds(0.0)
  call minvert(nobs,sum1,prod3)
  !dt=secnds(t1)
  !print*, dt
  
  !t1=secnds(0.0)
  prod4 = matmul(Sa,K_t)
  !dt=secnds(t1)
  !print*, dt
  
  !t1=secnds(0.0)
  prod5 = matmul(prod4,prod3)
  
  
  
  sum2=y-f
  prod6=matmul(K,dx)
  
  sum3=sum2+prod6(:,1)
  xstep=matmul(prod5,sum3)
  
  xnext=x+0.5*xstep
  dt=secnds(t1)
  print*, dt
  
  !If limits are present, ensure that variables are within hard limits
  if(present(xmax)) then
    do a=1,nvar	  
      if(xnext(a) .gt. xmax(a)) xnext(a) = xmax(a)
    end do
  endif
  if(present(xmin)) then
    do a=1,nvar	  
      if (xnext(a) .lt. xmin(a)) xnext(a) = xmin(a)
    end do
  endif

end subroutine oe_step_mirs

subroutine oe_diagnostics(nvar,nobs,x,xa,sa_i,F,y,sx,sy_i,prod2,A,X2)

  implicit none
  
  integer :: nvar,nobs
  
  real, intent(in) :: x(nvar), xa(nvar)
  real, intent(in) :: sa_i(nvar,nvar), sx(nvar,nvar), prod2(nvar,nvar)
  real, intent(in) :: y(nobs), F(nobs)
  real, intent(in) :: sy_i(nobs,nobs)
  real, intent(out) :: A(nvar,nvar)
  real, intent(out) :: X2
  
  real :: sum4(nobs,1), sum5(nvar,1), sum4_t(1,nobs), sum5_t(1,nvar)
  real :: prod8(1,nobs), prod9(1,1), prod10(1,nvar), prod11(1,1)
  A = matmul(sx,prod2) 
  sum4(:,1)=F-y
  sum5(:,1)=x-xa
  sum4_t=transpose(sum4)
  sum5_t=transpose(sum5)
  prod8=matmul(sum4_t,sy_i)
  prod9=matmul(prod8,sum4)
  prod10=matmul(sum5_t,sa_i)
  prod11=matmul(prod10,sum5)

  X2 = prod9(1,1)+prod11(1,1)

end subroutine oe_diagnostics
