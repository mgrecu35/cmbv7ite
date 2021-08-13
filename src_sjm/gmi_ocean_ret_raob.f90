module ocean_ret_raob

save

public  ::    gmi_ocean_ret_noclw, gmi_ocean_ret_clw1d, gmi_ocean_ret_clw2d

contains

SUBROUTINE gmi_ocean_ret_noclw(gdata, ret, relAz_in) !retrieval variables

  use profile_def
  use LUT_def
  use missingMod

  implicit none
  
  include 'constants_sjm.inc'
  include 'parametersMERRA.inc'
  include 'parametersGMI.inc'
  
  real*4  r4_normal_01
  character *8 :: fmt
  
  type(profile), intent(inout) :: gdata
  type(profile) :: gdata_p
  type(ret_data) :: ret
  
  real, intent(in), optional :: relAz_in
  
  real :: pc_std(1+2*MERRA_NLEV), prof_avg(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV)
  real :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  
  
  logical :: ch_mask(13)
  integer :: gmi_ch(13)
  real :: tb_out(13)
  integer :: nvar, neof, nobs, seed
!   !parameter (nvar=3)
! 
  real :: Tavg, Pavg, Vavg
  real :: emis, ebar, refl
  real :: tpw(100), tskin(100), tpw_mean, tskin_mean, tpw_rms, tskin_rms
  
  !variables for relative azimuth calc
  real :: x0, y0, x1, y1, dotprod, norm, du, dv, theta_w
!   ! Matrices to be used internally in optimal estimation expression
! 
  integer :: i, j, a, b, n, end_flag, max_iter, nlev
  parameter (max_iter=5)
  real, dimension(:), allocatable :: x, xa, xprime, xnext
  real, dimension(:), allocatable :: xmin, xmax, xstd, xbar
  !real, dimension(:), allocatable :: xall, xall_std, xall_bar
  real, dimension(:,:), allocatable :: sa, sa_i
  real, dimension(:,:), allocatable :: sy,sy_i
  real, dimension(:,:), allocatable :: K, K_t
  real, dimension(:,:), allocatable :: sx, sx_i, xdiff
  real, dimension(:,:), allocatable :: A_Matrix
  real, dimension(:), allocatable :: F, Fprime
  real, dimension(:), allocatable :: y
  real, dimension(:,:), allocatable :: prod6, prod2
  real :: Chi_Squared, prod7(1,1)

  seed=1983
  !First get the atmospheric mean profile and EOFS for this tskin value
  
  !if sfc pressure not present, extrapolate from tskin and sfc elev assuming 1013mb at sfc
  if(gdata%anc_psfc .gt. 0.) gdata%psfc = gdata%anc_psfc
  if(gdata%psfc .le. 0.) gdata%psfc = 1013.*exp(-1000.*gdata%sfc_elev*9.81/(float(gdata%tskin_index)*287.))
  !print*, gdata%tskin_index, gdata%sfc_elev, gdata%psfc
  
  prof_W = LUT%prof_W(1,:)
  prof_std = LUT%prof_std(1,:)
  pc_std = LUT%pc_std(1,:)
  prof_eofs = LUT%prof_eofs(1,:,:)
  prof_avg = 0.
  neof=5
  nvar=neof+3
  !stop
  !print*, nvar
  gmi_ch = 0
  !count channels with valid data
  nobs = 0
  ch_mask = .false.
  do i=1,13
!     print*, i, gdata%gmi_tb(i), gdata%eia, gdata%eia
    if(gdata%gmi_tb(i) .gt. 50. .and. gdata%gmi_tb(i) .lt. 350.) then
      nobs=nobs+1
      ch_mask(i) = .true.
      gmi_ch(nobs) = i
    endif
  end do
  ret%sim_tb0 = missing_r4
  ret%sim_tb = missing_r4
  ret%emis = missing_r4
  ret%semis = missing_r4
  ret%aemis = missing_r4
  ret%w10m_init = missing_r4
  ret%w10m_ret = missing_r4
  ret%sw10m_init = missing_r4
  ret%sw10m_ret = missing_r4
  ret%aw10m = missing_r4
  ret%relAz_init = missing_r4
  ret%relAz_ret = missing_r4
  ret%srelAz_init = missing_r4
  ret%srelAz_ret = missing_r4
  ret%arelAz = missing_r4
  ret%tpw_init = missing_r4
  ret%tpw_ret = missing_r4
  ret%stpw_init = missing_r4
  ret%stpw_ret = missing_r4
  ret%tskin_init = missing_r4
  ret%tskin_ret = missing_r4
  ret%stskin_init = missing_r4
  ret%stskin_ret = missing_r4
  ret%clwp_init = missing_r4
  ret%clwp_ret = missing_r4
  ret%clwp_init = missing_r4
  ret%sclwp_ret = missing_r4
  ret%aclwp = missing_r4
  ret%wc_ret(:) =missing_r4
  ret%swc(:) = missing_r4
  ret%awc(:) = 0.
  ret%rlwp_ret = 0.
  ret%iwp_ret = 0.
  ret%prate_sfc = 0.
  ret%chi_squared = missing_r4
  ret%obs_err=missing_r4
  ret%x_err=missing_r4
  ret%PIAKu = missing_r4
  ret%PIAKa = missing_r4
  ret%niter = 0
  ret%specularity=1.
  if(nobs .eq. 0) then
    
    return
  endif
  
  nvar=neof+3 !for wind speed and direction and SST, which is decoupled from the atmosphere
  allocate(x(nvar), xa(nvar), xprime(nvar), xnext(nvar), xdiff(1,nvar))
  !allocate(xall(5+1*MERRA_NLEV), xall_std(5+1*MERRA_NLEV), xall_bar(5+1*MERRA_NLEV))
  allocate(xmin(nvar), xmax(nvar), xstd(nvar), xbar(nvar))
  allocate(sa(nvar,nvar), sa_i(nvar,nvar))
  allocate(sy(nobs,nobs),sy_i(nobs,nobs))
  allocate(K(nobs,nvar), K_t(nvar,nobs))
  allocate(sx(nvar,nvar), sx_i(nvar,nvar))
  allocate(A_Matrix(nvar,nvar))
  allocate(F(nobs), Fprime(nobs))!, F_a(9), Fout(9)
  allocate(y(nobs))
  allocate(prod2(nvar,nvar),prod6(1,nvar))
  
  !initialize output

!-----Set bounds on state vector elements - converting all to normal distribution, so max at 5 standard deviations
  xmin = -5.
  xmax = 5.
  
      
!-----Declare the state covariance matrix
  sa(:,:)=0.
  do a=1,nvar
    sa(a,a)=1.
  end do
      
  !Standard Deviation of Retrieval Parameters
  xstd(1:neof) = pc_std(1:neof) !from raob error analysis
  !print*, xstd(1:neof)
  !stop
  !Mean value for retrieval parameters
  xbar(1:neof) = 0. !EOFs
  
  if(gdata%anc_u10m .ne. -99. .and. gdata%anc_v10m .ne. -99.) then
    xbar(neof+1) = sqrt(gdata%anc_u10m**2+gdata%anc_v10m**2)
    xstd(neof+1) = max(2.5,0.2*xbar(neof+1))
    !calculate relative azimuth
    !unfold longitude
    x0 = gdata%sclon
    if(gdata%sclon .gt. 90. .and. gdata%lon .lt. -90.) x0 = x0-360.
    if(gdata%sclon .lt. -90. .and. gdata%lon .gt. 90.) x0 = x0+360.
    x1 = gdata%lon
    y0 = gdata%sclat
    y1 = gdata%lat
    !print*, x0-x1,y0-y1
    dotprod = (x0-x1)*gdata%anc_u10m+(y0-y1)*gdata%anc_v10m
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(gdata%anc_u10m**2+gdata%anc_v10m**2)
    !print*, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      xbar(neof+2) = acos(dotprod/norm)/dtor 
      !set relAz_init to observed value (for evaluation)
      ret%relAz_init = acos(dotprod/norm)/dtor 
    else
      xbar(neof+2) = 0.
    endif

    du = 0.7*xstd(neof+1)/xbar(neof+1)*gdata%anc_v10m
    dv = -0.7*xstd(neof+1)/xbar(neof+1)*gdata%anc_u10m
    dotprod = (x0-x1)*(gdata%anc_u10m+du)+(y0-y1)*(gdata%anc_v10m+dv)
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt((gdata%anc_u10m+du)**2+(gdata%anc_v10m+dv)**2)
    !print*, du, dv, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      xstd(neof+2) =  abs(xbar(neof+2)-acos(dotprod/norm)/dtor)
      if(xstd(neof+2) .lt. 15.) xstd(neof+2) = 15.
    else
      xstd(neof+2) = 30.!abs(xbar(nvar)-0.)
    endif
  else
    xbar(neof+1) = 7.2!20. !wind
    xstd(neof+1) = 5.!10.
    xbar(neof+2) = 90.
    xstd(neof+2) = 900.
  endif
  
  !set relAz_init to observed value for evaluation
  if(gdata%iws .ge. 0. .and. gdata%iwd .ge. 0.) then
    x0 = gdata%sclon
    if(gdata%sclon .gt. 90. .and. gdata%lon .lt. -90.) x0 = x0-360.
    if(gdata%sclon .lt. -90. .and. gdata%lon .gt. 90.) x0 = x0+360.
    x1 = gdata%lon
    y0 = gdata%sclat
    y1 = gdata%lat
    !print*, x0-x1,y0-y1
  
    du = gdata%iws*cos(dtor*(270-gdata%iwd))
    dv = gdata%iws*sin(dtor*(270-gdata%iwd))
  
    dotprod = (x0-x1)*du+(y0-y1)*dv
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(du**2+dv**2)
    !print*, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      ret%aemis(13) = acos(dotprod/norm)/dtor 

    endif
  endif
  
  if(present(relAz_in)) then
    xbar(neof+2) = relAz_in
    xstd(neof+2) = 90.
  endif
  !ret%relAz_init = xbar(neof+2)
  !xbar(neof+2) = relAz_in
  !xstd(neof+2) = 45.
  !print*, gdata%anc_u10m, gdata%anc_v10m
  !print*, xbar(neof+1:nvar), xstd(neof+1:nvar)
  
  !limits for wind and clw
  xmin(neof+1) = (0.-xbar(neof+1))/xstd(neof+1)
  xmax(neof+1) = (50.-xbar(neof+1))/xstd(neof+1)
  xmin(neof+2) = (-90.-xbar(neof+2))/xstd(neof+2)
  xmax(neof+2) = (270.-xbar(neof+2))/xstd(neof+2)
  if(present(relAz_in)) then
    xmin(neof+2) = (-45.)/xstd(neof+2)
    xmax(neof+2) = (45.)/xstd(neof+2)
  endif

  if(gdata%anc_tskin .gt. 0.) then
    xbar(nvar) = gdata%anc_tskin
    xstd(nvar) = 1.
    xmin(nvar) = -3.
    xmax(nvar) = 3.
  else
    xbar(nvar) = 290.!gdata%tskin = prof_avg(1)
    xstd(nvar) = 10.
    xmin(nvar) = -2.
    xmax(nvar) = 2.
  endif
  !print*, xbar
  !print*, xstd
  !stop
  !set xall for tpw and tskin error comps.
  !xall = 0.
  !xall_bar = 0.
  !xall_std = 1.
  
  
  !Define variance in obs. No covariance in the observations
  !Eventually, will need to run a lot of orbits to estimate the true error covaraince matrix for pixels with 0 rain fraction and 0 land fraction.
  sy(:,:) = 0.
  !print*, nobs, gmi_ch
  do a=1,nobs
    do b=1,nobs
      !sy(a,b) = GMI_clr_errcov(gmi_ch(a),gmi_ch(b))
      sy(a,b) = LUT%eof_covmat_water(neof,gmi_ch(a),gmi_ch(b))
    end do
    sy(a,a) = sy(a,a)+GMI_NEDT(gmi_ch(a))**2
    !sy(a,a) = (0.0+GMI_NEDT(gmi_ch(a)))**2  !sensor NEDT + calibration stability
    !if(a .gt. 9) sy(a,a) = (3.+GMI_NEDT(GMI_nf(gmi_ch(a))))**2
    !print*, a, sqrt(sy(a,a))
  end do  
  !stop
!   !Add in non-uniformity variability using 85 as a proxy
!   sy(8,8) = (sqrt(sy(8,8))+abs(obs%TbL(8)-obs%TbH(1)))**2
!   sy(9,9) = (sqrt(sy(9,9))+abs(obs%TbL(9)-obs%TbH(2)))**2
  !print '(9F8.2)', sy
  !stop
  call minvert(nobs,sy,sy_i)
      
  !Initialize profile parameters from ancillary data or climatological means (based on tskin)
  nlev = MERRA_NLEV
  do while(gdata%psfc .lt. MERRA_PRESSLEV(MERRA_NLEV-nlev+1))
    nlev=nlev-1
  end do
         
  !Intialize a priori
  xa=0.

  !Initialize observation vector 
  y = gdata%gmi_tb(gmi_ch(1:nobs))

  call minvert(nvar,sa,sa_i)
! c-----GENERATE FIRST GUESS----- (currently setting first guess to apriori)

  x = xa

! c-----NEWTONIAN ITERATION-----
  end_flag=0
  n=0
  do while(end_flag .eq. 0)
    n=n+1 
    !print*, 'Iteration ', n
    
    call get_state_noclw_rh(nvar,neof,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)
    
    
    if(n .eq. 1) then ! get intial tpw and tskin uncertainty
      if(gdata%anc_tskin .gt. 0.) then
        call calc_tpw(gdata,nlev,ret%tpw_init)
        !ret%tskin_init = gdata%anc_tskin
      !else
        !ret%tpw_init = tpw_mean
        !ret%tskin_init = tskin_mean
      endif
      ret%tskin_init = x(nvar)*xstd(nvar)+xbar(nvar)
      !ret%stpw_init = tpw_rms
      !ret%stskin_init = tskin_rms
      ret%stskin_init = xstd(nvar)
      ret%w10m_init = gdata%w10m
      ret%sw10m_init = xstd(neof+2)
      !ret%relAz_init = gdata%relAz
      ret%srelAz_init = xstd(nvar)
      ret%clwp_init = 0.!exp(x(nvar)*xstd(nvar)+xbar(nvar))
      ret%sclwp_init = 0.!xstd(nvar)
      !print '(6F8.2)', ret%w10m_init, ret%clwp_init, ret%tpw_init, ret%tskin_init, tpw_mean, tskin_mean
      !print '(6F8.2)', ret%sw10m_init, ret%sclwp_init, ret%stpw_init, ret%stskin_init, tpw_rms, tskin_rms
    endif
    !print '(20F8.3)', x(1:neof), gdata%tskin, gdata%w10m, gdata%relAz
    !Simulate Tbs for this column
    !print*, 
    !print*, gdata%temp_lev
    !print*, gdata%qv_lev
    !print*, gdata%hgt_lev
    !print '(13F8.2)', y
    call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
    F = tb_out(gmi_ch(1:nobs))
    !print '(13F8.2)', F
    !stop
    if(n .eq. 1) then
      ret%sim_tb0(gmi_ch(1:nobs)) = tb_out(gmi_ch(1:nobs))
    endif
    
    !Calculate Jacobian
    
    do i=1,nvar
      !create perturbed atmosphere
      gdata_p = gdata
      xprime = x
      xprime(i) = x(i)+0.01
      call get_state_noclw_rh(nvar,neof,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata_p)
      call rtm_emission_csu(gdata_p,nlev,ch_mask,tb_out)
      Fprime = tb_out(gmi_ch(1:nobs))
      !print '(13F8.2)', Fprime
      do j=1,nobs
        K(j,i) = (Fprime(j)-F(j))/(xprime(i)-x(i))
      end do
      !print '(13F8.4)', K(:,i)
    end do
    
    call oe_step(nobs,nvar,x,xa,sa_i,y,F,sy_i,K,xnext,sx,sx_i,prod2,xmin,xmax)
    !Calculate convergence criterion (normalized difference between steps)
    xdiff(1,:)=x-xnext
    prod6 = matmul(xdiff,sx_i)
    prod7 = matmul(prod6,transpose(xdiff))

    !Reset x for next step
    x = 0.25*x+0.75*xnext !"memory" factor decreases steps to convergence in oscillations
    !Check for exit conditions
    if (prod7(1,1) .lt. 0.2*(nvar+nobs)) end_flag=1
    if (n .ge. max_iter) end_flag=1
    
    if (end_flag .eq. 1) then 
      
      !Evaluate A matrix and chi-squared
      call oe_diagnostics(nvar,nobs,x,xa,sa_i,F,y,sx,sy_i,prod2,A_Matrix,Chi_Squared)
      !write(fmt, '(A1,I2,A5)') '(',neof+3,'F8.4)'
      !print*, 'Sigma X:'
      !print fmt, sx
      !print*, 'A Matrix:'
      !print fmt, A_Matrix
      !print*, 'Chi-Squared:',Chi_Squared
      
      
      !run forward model and calculate Tbs
      call get_state_noclw_rh(nvar,neof,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)
      !print '(13F8.2)', y
      !simulate all Tbs, even if no used in retrieval
      !ch_mask(:) = .true.
      call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
      
      call calcPIA(gdata,nlev,ret%PIAKu,ret%PIAKa)
      !call calc_tpw(gdata,nlev,ret%tpw_ret)
      !print*, nlev
      call calc_tpw(gdata,nlev,tpw_mean)
      ret%tpw_ret = tpw_mean
      !print*, ret%tpw_init, ret%tpw_ret
      !print '(13F8.2)', F
      ret%sim_tb = tb_out
      !stop
      !diagnostic output
      do i=1,10
        xprime = x
        gdata_p=gdata
        do j=1,neof
          if(j .le. neof) xprime(j) = xprime(j)+sqrt(sx(j,j))*r4_normal_01(seed)
          !if(j .gt. neof) xall(j) = r4_normal_01(seed)
          !print*, j, xall(j)
        end do
        !stop
        call get_state_noclw_rh(nvar,neof,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata_p)
        call calc_tpw(gdata_p,nlev,tpw(i))
        tskin(i) = gdata_p%tskin
        !print*, tpw(i), tskin(i)
      end do
      
      tpw_mean = sum(tpw(1:10))/10.
      !print*, tpw_mean
!       !tskin_mean = sum(tskin)/100.
      tpw_rms = 0.
!       !tskin_rms = 0.
      do i=1,10
        tpw_rms = tpw_rms + (tpw(i)-tpw_mean)**2
!         !tskin_rms = tskin_rms + (tskin(i)-tskin_mean)**2
      end do
      tpw_rms = sqrt(tpw_rms/10.)
      !tskin_rms = sqrt(tskin_rms/100.)
      !ret%tpw_ret = tpw_mean
      !ret%tskin_ret = tskin_mean
      ret%tskin_ret = x(nvar)*xstd(nvar)+xbar(nvar)
      ret%stpw_ret = tpw_rms
      !ret%stskin_ret = tskin_rms
      ret%stskin_ret = xstd(nvar)*sqrt(sx(nvar,nvar))
      ret%atskin_ret = A_matrix(nvar,nvar)
      ret%w10m_ret = gdata%w10m
      ret%sw10m_ret = xstd(neof+1)*sqrt(sx(neof+1,neof+1))
      ret%aw10m = A_Matrix(neof+1,neof+1)
      ret%relAz_ret = gdata%relAz
      if(ret%relAz_ret .lt. 0.) ret%relAz_ret = -1.*ret%relAz_ret
      if(ret%relAz_ret .gt. 180.) ret%relAz_ret = 180.-ret%relAz_ret
      ret%srelAz_ret = xstd(neof+2)*sqrt(sx(neof+2,neof+2))
      ret%arelAz = A_Matrix(neof+2,neof+2)
      ret%clwp_init = gdata%anc_clwp
      ret%clwp_ret = 0.!exp(x(nvar)*xstd(nvar)+xbar(nvar))
      ret%sclwp_ret = 0.!xstd(nvar)*sqrt(sx(nvar,nvar))
      ret%aclwp = 0.!A_Matrix(nvar,nvar)
      !print*, ret%w10m_ret, ret%clwp_ret, ret%tpw_ret, ret%tskin_ret
      !print*, ret%sw10m_ret, ret%sclwp_ret, ret%stpw_ret, ret%stskin_ret
      !print*, ret%aw10m, ret%aclwp
      ret%niter = n
      ret%chi_squared = Chi_Squared/(nvar+nobs)
      ret%obs_err = 0.
      do i=1,nobs
        ret%obs_err = ret%obs_err+(F(i)-y(i))**2/(sy(i,i))
      end do
      ret%obs_err = ret%obs_err/nobs
      ret%x_err = 0.
      do i=1,nvar
        ret%x_err = ret%x_err+(x(i)-xa(i))**2/(sa(i,i))
      end do
      ret%x_err = ret%x_err/nvar
      ret%emis = missing_r4
      ret%aemis = missing_r4
      ret%semis = missing_r4
      do i=1,13
        if(ch_mask(i)) ret%emis(i) = gdata%emis(i)
      end do
      
      ret%err_eof = missing_r4
      ret%serr_eof = missing_r4
      ret%aerr_eof = missing_r4
      do i=1,min(10,neof)
        ret%err_eof(i) = x(i)*xstd(i)+xbar(i)
        ret%serr_eof(i) = sqrt(Sx(i,i))
        ret%aerr_eof(i) = A_Matrix(i,i)
      end do
      !stop
      do i =1,MERRA_NLEV
      	if(MERRA_PRESSLEV(i) .eq. gdata%press_lev(1)) exit
      end do
      ret%press_lev = MERRA_PRESSLEV
      ret%temp_lev(1:i) = gdata%temp_lev(1)
      ret%temp_lev(i:MERRA_NLEV) = gdata%temp_lev(1:MERRA_NLEV-i+1)
      ret%hgt_lev(i:MERRA_NLEV) = gdata%hgt_lev(1:MERRA_NLEV-i+1)
      ret%qv_lev(1:i) = gdata%qv_lev(1)
      ret%qv_lev(i:MERRA_NLEV) = gdata%qv_lev(1:MERRA_NLEV-i+1)
      ret%cloud_water(i:MERRA_NLEV) = gdata%cloud_water(1:MERRA_NLEV-i+1)
      ret%t2m_ret = gdata%t2m

      exit
    endif	!end final iteration
  end do	!end Newtonian loop
  
  deallocate(x, xa, xprime, xnext, xdiff)
  !deallocate(xall, xall_std, xall_bar)
  deallocate(xmin, xmax, xstd, xbar)
  deallocate(sa, sa_i)
  deallocate(sy,sy_i)
  deallocate(K, K_t)
  deallocate(sx, sx_i)
  deallocate(A_Matrix)
  deallocate(F, Fprime)!, F_a(9), Fout(9)
  deallocate(y)
  deallocate(prod2,prod6)

end subroutine gmi_ocean_ret_noclw

SUBROUTINE gmi_ocean_ret_clw1d(gdata, ret, relAz_in) !retrieval variables

  use profile_def
  use LUT_def
  use missingMod
  
  implicit none
  
  include 'constants_sjm.inc'
  include 'parametersMERRA.inc'
  include 'parametersGMI.inc'
  
  real*4  r4_normal_01
  character *8 :: fmt
  
  type(profile), intent(inout) :: gdata
  type(profile) :: gdata_p
  type(ret_data) :: ret
  
  real, intent(in), optional :: relAz_in
  
  real :: pc_std(1+2*MERRA_NLEV), prof_avg(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV)
  real :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  
  
  logical :: ch_mask(13)
  integer :: gmi_ch(13)
  real :: tb_out(13)
  integer :: nvar, neof, nobs, seed
!   !parameter (nvar=3)
! 
  real :: Tavg, Pavg, Vavg
  real :: emis, ebar, refl
  real :: tpw(100), tskin(100), tpw_mean, tskin_mean, tpw_rms, tskin_rms
  
  !variables for relative azimuth calc
  real :: x0, y0, x1, y1, dotprod, norm, du, dv, theta_w
!   ! Matrices to be used internally in optimal estimation expression
! 
  integer :: i, j, a, b, n, end_flag, max_iter, nlev
  parameter (max_iter=5)
  real, dimension(:), allocatable :: x, xa, xprime, xnext
  real, dimension(:), allocatable :: xmin, xmax, xstd, xbar
  !real, dimension(:), allocatable :: xall, xall_std, xall_bar
  real, dimension(:,:), allocatable :: sa, sa_i
  real, dimension(:,:), allocatable :: sy,sy_i
  real, dimension(:,:), allocatable :: K, K_t
  real, dimension(:,:), allocatable :: sx, sx_i, xdiff
  real, dimension(:,:), allocatable :: A_Matrix
  real, dimension(:), allocatable :: F, Fprime
  real, dimension(:), allocatable :: y
  real, dimension(:,:), allocatable :: prod6, prod2
  real :: Chi_Squared, prod7(1,1)

  seed=1983
  !First get the atmospheric mean profile and EOFS for this tskin value
  
  !if sfc pressure not present, extrapolate from tskin and sfc elev assuming 1013mb at sfc
  if(gdata%anc_psfc .gt. 0.) gdata%psfc = gdata%anc_psfc
  if(gdata%psfc .le. 0.) gdata%psfc = 1013.*exp(-1000.*gdata%sfc_elev*9.81/(float(gdata%tskin_index)*287.))
  !print*, gdata%tskin_index, gdata%sfc_elev, gdata%psfc
  
  prof_W = LUT%prof_W(1,:)
  prof_std = LUT%prof_std(1,:)
  pc_std = LUT%pc_std(1,:)
  prof_eofs = LUT%prof_eofs(1,:,:)
  prof_avg = 0.
  neof=5
  nvar=neof+3
  !stop
  !print*, nvar
  gmi_ch = 0
  !count channels with valid data
  nobs = 0
  ch_mask = .false.
  do i=1,13
!     print*, i, gdata%gmi_tb(i), gdata%eia, gdata%eia
    if(gdata%gmi_tb(i) .gt. 50. .and. gdata%gmi_tb(i) .lt. 350.) then
      nobs=nobs+1
      ch_mask(i) = .true.
      gmi_ch(nobs) = i
    endif
  end do
  ret%sim_tb0 = missing_r4
  ret%sim_tb = missing_r4
  ret%emis = missing_r4
  ret%semis = missing_r4
  ret%aemis = missing_r4
  ret%w10m_init = missing_r4
  ret%w10m_ret = missing_r4
  ret%sw10m_init = missing_r4
  ret%sw10m_ret = missing_r4
  ret%aw10m = missing_r4
  ret%relAz_init = missing_r4
  ret%relAz_ret = missing_r4
  ret%srelAz_init = missing_r4
  ret%srelAz_ret = missing_r4
  ret%arelAz = missing_r4
  ret%tpw_init = missing_r4
  ret%tpw_ret = missing_r4
  ret%stpw_init = missing_r4
  ret%stpw_ret = missing_r4
  ret%tskin_init = missing_r4
  ret%tskin_ret = missing_r4
  ret%stskin_init = missing_r4
  ret%stskin_ret = missing_r4
  ret%clwp_init = missing_r4
  ret%clwp_ret = missing_r4
  ret%clwp_init = missing_r4
  ret%sclwp_ret = missing_r4
  ret%aclwp = missing_r4
  ret%wc_ret(:) = 0.
  ret%swc(:) = missing_r4
  ret%awc(:) = 0.
  ret%rlwp_ret = 0.
  ret%iwp_ret = 0.
  ret%prate_sfc = 0.
  ret%chi_squared = missing_r4
  ret%obs_err=missing_r4
  ret%x_err=missing_r4
  ret%PIAKu = missing_r4
  ret%PIAKa = missing_r4
  ret%niter = 0
  ret%specularity=1.
  if(nobs .eq. 0) then
    
    return
  endif
  
  nvar=neof+4 !for wind speed and direction and SST, which is decoupled from the atmosphere
  allocate(x(nvar), xa(nvar), xprime(nvar), xnext(nvar), xdiff(1,nvar))
  !allocate(xall(5+1*MERRA_NLEV), xall_std(5+1*MERRA_NLEV), xall_bar(5+1*MERRA_NLEV))
  allocate(xmin(nvar), xmax(nvar), xstd(nvar), xbar(nvar))
  allocate(sa(nvar,nvar), sa_i(nvar,nvar))
  allocate(sy(nobs,nobs),sy_i(nobs,nobs))
  allocate(K(nobs,nvar), K_t(nvar,nobs))
  allocate(sx(nvar,nvar), sx_i(nvar,nvar))
  allocate(A_Matrix(nvar,nvar))
  allocate(F(nobs), Fprime(nobs))!, F_a(9), Fout(9)
  allocate(y(nobs))
  allocate(prod2(nvar,nvar),prod6(1,nvar))
  
  !initialize output

!-----Set bounds on state vector elements - converting all to normal distribution, so max at 5 standard deviations
  xmin = -5.
  xmax = 5.
  
      
!-----Declare the state covariance matrix
  sa(:,:)=0.
  do a=1,nvar
    sa(a,a)=1.
  end do
      
  !Standard Deviation of Retrieval Parameters
  xstd(1:neof) = pc_std(1:neof) !from raob error analysis
  !print*, xstd(1:neof)
  !stop
  !Mean value for retrieval parameters
  xbar(1:neof) = 0. !EOFs
  
  if(gdata%anc_u10m .ne. -99. .and. gdata%anc_v10m .ne. -99.) then
    xbar(neof+1) = sqrt(gdata%anc_u10m**2+gdata%anc_v10m**2)
    xstd(neof+1) = max(2.5,0.2*xbar(neof+1))
    !calculate relative azimuth
    !unfold longitude
    x0 = gdata%sclon
    if(gdata%sclon .gt. 90. .and. gdata%lon .lt. -90.) x0 = x0-360.
    if(gdata%sclon .lt. -90. .and. gdata%lon .gt. 90.) x0 = x0+360.
    x1 = gdata%lon
    y0 = gdata%sclat
    y1 = gdata%lat
    !print*, x0-x1,y0-y1
    dotprod = (x0-x1)*gdata%anc_u10m+(y0-y1)*gdata%anc_v10m
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(gdata%anc_u10m**2+gdata%anc_v10m**2)
    !print*, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      xbar(neof+2) = acos(dotprod/norm)/dtor 
      !set relAz_init to observed value (for evaluation)
      ret%relAz_init = acos(dotprod/norm)/dtor 
    else
      xbar(neof+2) = 0.
    endif

    du = 0.7*xstd(neof+1)/xbar(neof+1)*gdata%anc_v10m
    dv = -0.7*xstd(neof+1)/xbar(neof+1)*gdata%anc_u10m
    dotprod = (x0-x1)*(gdata%anc_u10m+du)+(y0-y1)*(gdata%anc_v10m+dv)
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt((gdata%anc_u10m+du)**2+(gdata%anc_v10m+dv)**2)
    !print*, du, dv, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      xstd(neof+2) =  abs(xbar(neof+2)-acos(dotprod/norm)/dtor)
      if(xstd(neof+2) .lt. 15.) xstd(neof+2) = 15.
    else
      xstd(neof+2) = 30.!abs(xbar(nvar)-0.)
    endif
  else
    xbar(neof+1) = 7.2!20. !wind
    xstd(neof+1) = 5.!10.
    xbar(neof+2) = 90.
    xstd(neof+2) = 900.
  endif
  
  !set relAz_init to observed value for evaluation
  if(gdata%iws .ge. 0. .and. gdata%iwd .ge. 0.) then
    x0 = gdata%sclon
    if(gdata%sclon .gt. 90. .and. gdata%lon .lt. -90.) x0 = x0-360.
    if(gdata%sclon .lt. -90. .and. gdata%lon .gt. 90.) x0 = x0+360.
    x1 = gdata%lon
    y0 = gdata%sclat
    y1 = gdata%lat
    !print*, x0-x1,y0-y1
  
    du = gdata%iws*cos(dtor*(270-gdata%iwd))
    dv = gdata%iws*sin(dtor*(270-gdata%iwd))
  
    dotprod = (x0-x1)*du+(y0-y1)*dv
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(du**2+dv**2)
    !print*, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      ret%aemis(13) = acos(dotprod/norm)/dtor 

    endif
  endif
  
  if(present(relAz_in)) then
    xbar(neof+2) = relAz_in
    xstd(neof+2) = 90.
  endif
  
  !limits for wind and clw
  xmin(neof+1) = (0.-xbar(neof+1))/xstd(neof+1)
  xmax(neof+1) = (50.-xbar(neof+1))/xstd(neof+1)
  xmin(neof+2) = (-90.-xbar(neof+2))/xstd(neof+2)
  xmax(neof+2) = (270.-xbar(neof+2))/xstd(neof+2)
  if(present(relAz_in)) then
    xmin(neof+2) = (-45.)/xstd(neof+2)
    xmax(neof+2) = (45.)/xstd(neof+2)
  endif

  if(gdata%anc_tskin .gt. 0.) then
    xbar(neof+3) = gdata%anc_tskin
    xstd(neof+3) = 1.
    xmin(neof+3) = -3.
    xmax(neof+3) = 3.
  else
    xbar(neof+3) = 290.!gdata%tskin = prof_avg(1)
    xstd(neof+3) = 10.
    xmin(neof+3) = -2.
    xmax(neof+3) = 2.
  endif
  
  !cloud water
  if(gdata%anc_clwp .gt. 0.01) then
    xbar(neof+4) = alog(gdata%anc_clwp)
    xstd(neof+4) = 2
    xmin(neof+4) = alog(0.001)
    xmax(neof+4) = alog(5.)
  else
    xbar(neof+4) = alog(0.01)
    xstd(neof+4) = 3
    xmin(neof+4) = alog(0.001)
    xmax(neof+4) = alog(5.)
  endif
  !print*, gdata%anc_clwp, xbar(neof+4), xmin(neof+4), xmax(neof+4)
  !print*, xbar
  !print*, xstd
  !stop
  !set xall for tpw and tskin error comps.
  !xall = 0.
  !xall_bar = 0.
  !xall_std = 1.
  
  
  !Define variance in obs. No covariance in the observations
  !Eventually, will need to run a lot of orbits to estimate the true error covaraince matrix for pixels with 0 rain fraction and 0 land fraction.
  sy(:,:) = 0.
  do a=1,nobs
    do b=1,nobs
      !sy(a,b) = GMI_clr_errcov(gmi_ch(a),gmi_ch(b))
      sy(a,b) = LUT%eof_covmat_water(neof,gmi_ch(a),gmi_ch(b))
    end do
    sy(a,a) = sy(a,a)+GMI_NEDT(gmi_ch(a))**2
    !sy(a,a) = (0.0+GMI_NEDT(gmi_ch(a)))**2  !sensor NEDT + calibration stability
    !if(a .gt. 9) sy(a,a) = (3.+GMI_NEDT(GMI_nf(gmi_ch(a))))**2
    !print*, a, sqrt(sy(a,a))
  end do  
  !stop
!   !Add in non-uniformity variability using 85 as a proxy
!   sy(8,8) = (sqrt(sy(8,8))+abs(obs%TbL(8)-obs%TbH(1)))**2
!   sy(9,9) = (sqrt(sy(9,9))+abs(obs%TbL(9)-obs%TbH(2)))**2
  !print '(9F8.2)', sy
  !stop
  call minvert(nobs,sy,sy_i)
      
  !Initialize profile parameters from ancillary data or climatological means (based on tskin)
  nlev = MERRA_NLEV
  do while(gdata%psfc .lt. MERRA_PRESSLEV(MERRA_NLEV-nlev+1))
    nlev=nlev-1
  end do
         
  !Intialize a priori
  xa=0.

  !Initialize observation vector 
  y = gdata%gmi_tb(gmi_ch(1:nobs))

  call minvert(nvar,sa,sa_i)
! c-----GENERATE FIRST GUESS----- (currently setting first guess to apriori)

  x = xa

! c-----NEWTONIAN ITERATION-----
  end_flag=0
  n=0
  do while(end_flag .eq. 0)
    n=n+1 
    !print*, 'Iteration ', n
    !print*, x
    !print*, xmax
    call get_state_clw1d_rh(nvar,neof,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)
    
    
    if(n .eq. 1) then ! get intial tpw and tskin uncertainty
      if(gdata%anc_tskin .gt. 0.) then
        call calc_tpw(gdata,nlev,ret%tpw_init)
        !ret%tskin_init = gdata%anc_tskin
      !else
        !ret%tpw_init = tpw_mean
        !ret%tskin_init = tskin_mean
      endif
      ret%tskin_init = x(neof+3)*xstd(neof+3)+xbar(neof+3)
      ret%stpw_init = tpw_rms
      !ret%stskin_init = tskin_rms
      ret%stskin_init = xstd(neof+3)
      ret%w10m_init = gdata%w10m
      ret%sw10m_init = xstd(neof+1)
      !ret%relAz_init = gdata%relAz
      ret%srelAz_init = xstd(neof+2)
      ret%clwp_init = exp(x(neof+4)*xstd(neof+4)+xbar(neof+4))
      ret%sclwp_init = xstd(neof+4)
      !print '(6F8.2)', ret%w10m_init, ret%clwp_init, ret%tpw_init, ret%tskin_init, tpw_mean, tskin_mean
      !print '(6F8.2)', ret%sw10m_init, ret%sclwp_init, ret%stpw_init, ret%stskin_init, tpw_rms, tskin_rms
    endif
    !print '(20F8.3)', x(1:neof), gdata%tskin, gdata%w10m, gdata%relAz, exp(x(neof+4)*xstd(neof+4)+xbar(neof+4))
    !Simulate Tbs for this column
    !print*, 
    !print*, gdata%temp_lev
    !print*, gdata%qv_lev
    !print*, gdata%hgt_lev
    !print '(13F8.2)', y
    call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
    F = tb_out(gmi_ch(1:nobs))
    !print '(13F8.2)', F
    
    if(n .eq. 1) then
      ret%sim_tb0(gmi_ch(1:nobs)) = tb_out(gmi_ch(1:nobs))
    endif
    
    !Calculate Jacobian
    
    do i=1,nvar
      !create perturbed atmosphere
      gdata_p = gdata
      xprime = x
      xprime(i) = x(i)+0.01
      call get_state_clw1d_rh(nvar,neof,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata_p)
      call rtm_emission_csu(gdata_p,nlev,ch_mask,tb_out)
      Fprime = tb_out(gmi_ch(1:nobs))
      !print '(13F8.2)', Fprime
      do j=1,nobs
        K(j,i) = (Fprime(j)-F(j))/(xprime(i)-x(i))
      end do
      !print '(13F8.4)', K(:,i)
    end do
    
    call oe_step(nobs,nvar,x,xa,sa_i,y,F,sy_i,K,xnext,sx,sx_i,prod2,xmin,xmax)
    !Calculate convergence criterion (normalized difference between steps)
    xdiff(1,:)=x-xnext
    prod6 = matmul(xdiff,sx_i)
    prod7 = matmul(prod6,transpose(xdiff))

    !Reset x for next step
    x = 0.25*x+0.75*xnext !"memory" factor decreases steps to convergence in oscillations
    !Check for exit conditions
    if (prod7(1,1) .lt. 0.2*(nvar+nobs)) end_flag=1
    if (n .ge. max_iter) end_flag=1
    
    if (end_flag .eq. 1) then 
      
      !Evaluate A matrix and chi-squared
      call oe_diagnostics(nvar,nobs,x,xa,sa_i,F,y,sx,sy_i,prod2,A_Matrix,Chi_Squared)
      !write(fmt, '(A1,I2,A5)') '(',neof+4,'F8.4)'
      !print*, 'Sigma X:'
      !print fmt, sx
      !print*, 'A Matrix:'
      !print fmt, A_Matrix
      !print*, 'Chi-Squared:',Chi_Squared
      
      
      !run forward model and calculate Tbs
      call get_state_clw1d_rh(nvar,neof,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)
      !print '(13F8.2)', y
      !simulate all Tbs, even if no used in retrieval
      !ch_mask(:) = .true.
      call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
      
      call calcPIA(gdata,nlev,ret%PIAKu,ret%PIAKa)
      !call calc_tpw(gdata,nlev,ret%tpw_ret)
      !print*, nlev
      call calc_tpw(gdata,nlev,tpw_mean)
      ret%tpw_ret = tpw_mean
      
      !print*, ret%tpw_init, ret%tpw_ret
      !print '(13F8.2)', F
      ret%sim_tb = tb_out
      
      !diagnostic output
      do i=1,10
        xprime = x
        gdata_p=gdata
        do j=1,neof
           if(sx(j,j)<1e-5) then
              print*, 'sx=',sx(j,j)
              sx(j,j)=1e-5
           endif
           if(j .le. neof) xprime(j) = xprime(j)+sqrt(sx(j,j))*r4_normal_01(seed)
           !if(j .gt. neof) xall(j) = r4_normal_01(seed)
           !print*, j, xall(j)
        end do
        !stop
        !call check_nan(gdata_p%hgt_lev,gdata_p%temp_lev)
        call get_state_clw1d_rh(nvar,neof,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata_p)
        call calc_tpw(gdata_p,nlev,tpw(i))
        tskin(i) = gdata_p%tskin
        !print*, tpw(i), tskin(i)
      end do
      
      tpw_mean = sum(tpw(1:10))/10.
      !print*, tpw_mean
!       !tskin_mean = sum(tskin)/100.
      tpw_rms = 0.
!       !tskin_rms = 0.
      do i=1,10
        tpw_rms = tpw_rms + (tpw(i)-tpw_mean)**2
!         !tskin_rms = tskin_rms + (tskin(i)-tskin_mean)**2
      end do
      tpw_rms = sqrt(tpw_rms/10.)
      !tskin_rms = sqrt(tskin_rms/100.)
      !ret%tpw_ret = tpw_mean
      !ret%tskin_ret = tskin_mean
      ret%tskin_ret = x(neof+3)*xstd(neof+3)+xbar(neof+3)
      ret%stpw_ret = tpw_rms
      !ret%stskin_ret = tskin_rms
      ret%stskin_ret = xstd(neof+3)*sqrt(sx(neof+3,neof+3))
      ret%atskin_ret = A_matrix(neof+3,neof+3)
      ret%w10m_ret = gdata%w10m
      ret%sw10m_ret = xstd(neof+1)*sqrt(sx(neof+1,neof+1))
      ret%aw10m = A_Matrix(neof+1,neof+1)
      ret%relAz_ret = gdata%relAz
      if(ret%relAz_ret .lt. 0.) ret%relAz_ret = -1.*ret%relAz_ret
      if(ret%relAz_ret .gt. 180.) ret%relAz_ret = 180.-ret%relAz_ret
      ret%srelAz_ret = xstd(neof+2)*sqrt(sx(neof+2,neof+2))
      ret%arelAz = A_Matrix(neof+2,neof+2)
      ret%clwp_init = gdata%anc_clwp
      ret%clwp_ret = exp(x(neof+4)*xstd(neof+4)+xbar(neof+4))
      ret%sclwp_ret = xstd(neof+4)*sqrt(sx(neof+4,neof+4))
      ret%aclwp = A_Matrix(neof+4,neof+4)
      !print*, ret%w10m_ret, ret%clwp_ret, ret%tpw_ret, ret%tskin_ret
      !print*, ret%sw10m_ret, ret%sclwp_ret, ret%stpw_ret, ret%stskin_ret
      !print*, ret%aw10m, ret%aclwp
      ret%niter = n
      ret%chi_squared = Chi_Squared/(nvar+nobs)
      ret%obs_err = 0.
      do i=1,nobs
        ret%obs_err = ret%obs_err+(F(i)-y(i))**2/(sy(i,i))
      end do
      ret%obs_err = ret%obs_err/nobs
      ret%x_err = 0.
      do i=1,nvar
        ret%x_err = ret%x_err+(x(i)-xa(i))**2/(sa(i,i))
      end do
      ret%x_err = ret%x_err/nvar
      ret%emis = missing_r4
      ret%aemis = missing_r4
      ret%semis = missing_r4
      do i=1,13
        if(ch_mask(i)) ret%emis(i) = gdata%emis(i)
      end do
      
      ret%err_eof = missing_r4
      ret%serr_eof = missing_r4
      ret%aerr_eof = missing_r4
      do i=1,min(10,neof)
        ret%err_eof(i) = x(i)*xstd(i)+xbar(i)
        ret%serr_eof(i) = sqrt(Sx(i,i))
        ret%aerr_eof(i) = A_Matrix(i,i)
      end do
      !stop
      do i =1,MERRA_NLEV
      	if(MERRA_PRESSLEV(i) .eq. gdata%press_lev(1)) exit
      end do
      ret%press_lev = MERRA_PRESSLEV
      ret%temp_lev(1:i) = gdata%temp_lev(1)
      ret%temp_lev(i:MERRA_NLEV) = gdata%temp_lev(1:MERRA_NLEV-i+1)
      ret%hgt_lev(i:MERRA_NLEV) = gdata%hgt_lev(1:MERRA_NLEV-i+1)
      ret%qv_lev(1:i) = gdata%qv_lev(1)
      ret%qv_lev(i:MERRA_NLEV) = gdata%qv_lev(1:MERRA_NLEV-i+1)
      ret%cloud_water(i:MERRA_NLEV) = gdata%cloud_water(1:MERRA_NLEV-i+1)
      ret%t2m_ret = gdata%t2m

      exit
    endif	!end final iteration
  end do	!end Newtonian loop
  
  deallocate(x, xa, xprime, xnext, xdiff)
  !deallocate(xall, xall_std, xall_bar)
  deallocate(xmin, xmax, xstd, xbar)
  deallocate(sa, sa_i)
  deallocate(sy,sy_i)
  deallocate(K, K_t)
  deallocate(sx, sx_i)
  deallocate(A_Matrix)
  deallocate(F, Fprime)!, F_a(9), Fout(9)
  deallocate(y)
  deallocate(prod2,prod6)

end subroutine gmi_ocean_ret_clw1d

SUBROUTINE gmi_ocean_ret_clw2d(gdata,tb2d, gmi_gain, ret, relAz_in) !retrieval variables

  use profile_def
  use LUT_def
  use missingMod

  implicit none
  
  include 'constants_sjm.inc'
  include 'parametersMERRA.inc'
  include 'parametersGMI.inc'
  
  real*4  r4_normal_01
  character *8 :: fmt
  
  type(profile), intent(inout) :: gdata
  type(profile) :: gdata_p
  type(ret_data) :: ret
  
  real, intent(in) :: tb2d(13,9), gmi_gain(13,9)
  real, intent(in), optional :: relAz_in
  
  real :: pc_std(1+2*MERRA_NLEV), prof_avg(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV)
  real :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  real :: clwprofs(9,MERRA_NLEV)
  
  logical :: ch_mask(13), ch_mask2(13)
  integer :: gmi_ch(13+18), y2locs(13*9), clwlocs(9)
  real :: tb_out(13), tb_out_2d(13,9), tb_out_2dp(13,9)
  integer :: neof, nclw, nvar, nobs, nobs1, nobs2, seed
!   !parameter (nvar=3)
! 
  real :: Tavg, Pavg, Vavg
  real :: emis, ebar, refl
  real :: tpw(100), tskin(100), tpw_mean, tskin_mean, tpw_rms, tskin_rms
  
  !variables for relative azimuth calc
  real :: x0, y0, x1, y1, dotprod, norm, du, dv, theta_w
!   ! Matrices to be used internally in optimal estimation expression
! 
  integer :: i, j, a, b, n, end_flag, max_iter, nlev
  parameter (max_iter=5)
  real, dimension(:), allocatable :: x, xa, xprime, xnext
  real, dimension(:), allocatable :: xmin, xmax, xstd, xbar
  real, dimension(:), allocatable :: xall, xall_std, xall_bar
  real, dimension(:,:), allocatable :: sa, sa_i
  real, dimension(:,:), allocatable :: sy,sy_i
  real, dimension(:,:), allocatable :: K, K_t
  real, dimension(:,:), allocatable :: sx, sx_i, xdiff
  real, dimension(:,:), allocatable :: A_Matrix
  real, dimension(:), allocatable :: F, Fprime
  real, dimension(:), allocatable :: y
  real, dimension(:,:), allocatable :: prod6, prod2
  real :: Chi_Squared, prod7(1,1)

  data ch_mask2 /.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.false.,.false.,.false.,.false./
  
  seed=1983
  !First get the atmospheric mean profile and EOFS for this tskin value
  
  !if sfc pressure not present, extrapolate from tskin and sfc elev assuming 1013mb at sfc
  if(gdata%anc_psfc .gt. 0.) gdata%psfc = gdata%anc_psfc
  if(gdata%psfc .le. 0.) gdata%psfc = 1013.*exp(-1000.*gdata%sfc_elev*9.81/(float(gdata%tskin_index)*287.))
  !print*, gdata%tskin_index, gdata%sfc_elev, gdata%psfc
  
  prof_W = LUT%prof_W(1,:)
  prof_std = LUT%prof_std(1,:)
  pc_std = LUT%pc_std(1,:)
  prof_eofs = LUT%prof_eofs(1,:,:)
  prof_avg = 0.
  neof=5
  nvar=neof+3
  !stop
  !print*, nvar
  gmi_ch = 0
  !count channels with valid data
  nobs1 = 0 !1D Tb vector
  nobs2 = 0 !2D Tbs
  nclw = 0
  ch_mask = .false.
  do i=1,13
!     print*, i, gdata%gmi_tb(i), gdata%eia, gdata%eia
    if(gdata%gmi_tb(i) .gt. 50. .and. gdata%gmi_tb(i) .lt. 350.) then
      nobs1=nobs1+1
      ch_mask(i) = .true.
      gmi_ch(nobs1) = i
    endif
  end do
  !2D Tb vector
  do j = 1,9
    do i=8,9 !only using 89V & 89H for cloud now; later test 166 V&H
      !print*, tb2d(i,j), gmi_gain(1,j)
      if(ch_mask(i) .and. tb2d(i,j) .gt. 50. .and. tb2d(i,j) .lt. 350. .and. j .ne. 5 .and. gmi_gain(1,j) .gt. 0.) then
        nobs2=nobs2+1
        gmi_ch(nobs1+nobs2) = i
        y2locs(nobs2) = j
      endif
    end do
    !if(gmi_gain(1,j) .gt. 0.) then
      nclw = nclw+1
      clwlocs(nclw) = j
    !endif
  end do
  nobs=nobs1+nobs2
  nvar = neof+nclw+3
  !print*, nobs1, nobs2, nobs
  !print*, nclw, neof, nvar
  !print*, clwlocs
  !stop
  ret%sim_tb0 = missing_r4
  ret%sim_tb = missing_r4
  ret%emis = missing_r4
  ret%semis = missing_r4
  ret%aemis = missing_r4
  ret%w10m_init = missing_r4
  ret%w10m_ret = missing_r4
  ret%sw10m_init = missing_r4
  ret%sw10m_ret = missing_r4
  ret%aw10m = missing_r4
  ret%relAz_init = missing_r4
  ret%relAz_ret = missing_r4
  ret%srelAz_init = missing_r4
  ret%srelAz_ret = missing_r4
  ret%arelAz = missing_r4
  ret%tpw_init = missing_r4
  ret%tpw_ret = missing_r4
  ret%stpw_init = missing_r4
  ret%stpw_ret = missing_r4
  ret%tskin_init = missing_r4
  ret%tskin_ret = missing_r4
  ret%stskin_init = missing_r4
  ret%stskin_ret = missing_r4
  ret%clwp_init = missing_r4
  ret%clwp_ret = missing_r4
  ret%clwp_init = missing_r4
  ret%sclwp_ret = missing_r4
  ret%aclwp = missing_r4
  ret%wc_ret(:) = 0.
  ret%swc(:) = missing_r4
  ret%awc(:) = 0.
  ret%rlwp_ret = 0.
  ret%iwp_ret = 0.
  ret%prate_sfc = 0.
  ret%chi_squared = missing_r4
  ret%obs_err=missing_r4
  ret%x_err=missing_r4
  ret%PIAKu = missing_r4
  ret%PIAKa = missing_r4
  ret%niter = 0
  ret%specularity=1.
  if(nobs .eq. 0) then
    
    return
  endif
  
  allocate(x(nvar), xa(nvar), xprime(nvar), xnext(nvar), xdiff(1,nvar))
  allocate(xall(5+1*MERRA_NLEV), xall_std(5+1*MERRA_NLEV), xall_bar(5+1*MERRA_NLEV))
  allocate(xmin(nvar), xmax(nvar), xstd(nvar), xbar(nvar))
  allocate(sa(nvar,nvar), sa_i(nvar,nvar))
  allocate(sy(nobs,nobs),sy_i(nobs,nobs))
  allocate(K(nobs,nvar), K_t(nvar,nobs))
  allocate(sx(nvar,nvar), sx_i(nvar,nvar))
  allocate(A_Matrix(nvar,nvar))
  allocate(F(nobs), Fprime(nobs))!, F_a(9), Fout(9)
  allocate(y(nobs))
  allocate(prod2(nvar,nvar),prod6(1,nvar))
  
  !initialize output

!-----Set bounds on state vector elements - converting all to normal distribution, so max at 5 standard deviations
  xmin = -5.
  xmax = 5.
  
      
!-----Declare the state covariance matrix
  sa(:,:)=0.
  do a=1,nvar
    sa(a,a)=1.
  end do
      
  !Standard Deviation of Retrieval Parameters
  xstd(1:neof) = pc_std(1:neof) !from raob error analysis
  !print*, xstd(1:neof)
  !stop
  !Mean value for retrieval parameters
  xbar(1:neof) = 0. !EOFs
  
  if(gdata%anc_u10m .ne. -99. .and. gdata%anc_v10m .ne. -99.) then
    xbar(neof+1) = sqrt(gdata%anc_u10m**2+gdata%anc_v10m**2)
    xstd(neof+1) = max(2.5,0.2*xbar(neof+1))
    !calculate relative azimuth
    !unfold longitude
    x0 = gdata%sclon
    if(gdata%sclon .gt. 90. .and. gdata%lon .lt. -90.) x0 = x0-360.
    if(gdata%sclon .lt. -90. .and. gdata%lon .gt. 90.) x0 = x0+360.
    x1 = gdata%lon
    y0 = gdata%sclat
    y1 = gdata%lat
    !print*, x0-x1,y0-y1
    dotprod = (x0-x1)*gdata%anc_u10m+(y0-y1)*gdata%anc_v10m
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(gdata%anc_u10m**2+gdata%anc_v10m**2)
    !print*, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      xbar(neof+2) = acos(dotprod/norm)/dtor 
      !set relAz_init to observed value (for evaluation)
      ret%relAz_init = acos(dotprod/norm)/dtor 
    else
      xbar(neof+2) = 0.
    endif

    du = 0.7*xstd(neof+1)/xbar(neof+1)*gdata%anc_v10m
    dv = -0.7*xstd(neof+1)/xbar(neof+1)*gdata%anc_u10m
    dotprod = (x0-x1)*(gdata%anc_u10m+du)+(y0-y1)*(gdata%anc_v10m+dv)
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt((gdata%anc_u10m+du)**2+(gdata%anc_v10m+dv)**2)
    !print*, du, dv, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      xstd(neof+2) =  abs(xbar(neof+2)-acos(dotprod/norm)/dtor)
      if(xstd(neof+2) .lt. 15.) xstd(neof+2) = 15.
    else
      xstd(neof+2) = 30.!abs(xbar(nvar)-0.)
    endif
  else
    xbar(neof+1) = 7.2!20. !wind
    xstd(neof+1) = 5.!10.
    xbar(neof+2) = 90.
    xstd(neof+2) = 900.
  endif
  
  !set relAz_init to observed value for evaluation
  if(gdata%iws .ge. 0. .and. gdata%iwd .ge. 0.) then
    x0 = gdata%sclon
    if(gdata%sclon .gt. 90. .and. gdata%lon .lt. -90.) x0 = x0-360.
    if(gdata%sclon .lt. -90. .and. gdata%lon .gt. 90.) x0 = x0+360.
    x1 = gdata%lon
    y0 = gdata%sclat
    y1 = gdata%lat
    !print*, x0-x1,y0-y1
  
    du = gdata%iws*cos(dtor*(270-gdata%iwd))
    dv = gdata%iws*sin(dtor*(270-gdata%iwd))
  
    dotprod = (x0-x1)*du+(y0-y1)*dv
    norm = sqrt((x0-x1)**2+(y0-y1)**2)*sqrt(du**2+dv**2)
    !print*, dotprod, norm
    if(norm .le. 0.) norm = 0.01*sqrt((x0-x1)**2+(y0-y1)**2)
    if(abs(dotprod/norm) .lt. 1.) then
      ret%aemis(13) = acos(dotprod/norm)/dtor 

    endif
  endif
  
  if(present(relAz_in)) then
    xbar(neof+2) = relAz_in
    xstd(neof+2) = 90.
  endif
  
  !limits for wind and clw
  xmin(neof+1) = (0.-xbar(neof+1))/xstd(neof+1)
  xmax(neof+1) = (50.-xbar(neof+1))/xstd(neof+1)
  xmin(neof+2) = (-90.-xbar(neof+2))/xstd(neof+2)
  xmax(neof+2) = (270.-xbar(neof+2))/xstd(neof+2)
  if(present(relAz_in)) then
    xmin(neof+2) = (-45.)/xstd(neof+2)
    xmax(neof+2) = (45.)/xstd(neof+2)
  endif

  if(gdata%anc_tskin .gt. 0.) then
    xbar(neof+3) = gdata%anc_tskin
    xstd(neof+3) = 1.
    xmin(neof+3) = -3.
    xmax(neof+3) = 3.
  else
    xbar(neof+3) = 290.!gdata%tskin = prof_avg(1)
    xstd(neof+3) = 10.
    xmin(neof+3) = -2.
    xmax(neof+3) = 2.
  endif
  
  !cloud water
  if(gdata%anc_clwp .gt. 0.01) then
    xbar(neof+4:nvar) = alog(gdata%anc_clwp)
    xstd(neof+4:nvar) = 2
    xmin(neof+4:nvar) = alog(0.001)
    xmax(neof+4:nvar) = alog(5.)
  else
    xbar(neof+4:nvar) = alog(0.01)
    xstd(neof+4:nvar) = 3
    xmin(neof+4:nvar) = alog(0.001)
    xmax(neof+4:nvar) = alog(5.)
  endif
  !print*, gdata%anc_clwp, xbar(neof+4), xmin(neof+4), xmax(neof+4)
  !print*, xbar
  !print*, xstd
  !stop
  !set xall for tpw and tskin error comps.
  xall = 0.
  xall_bar = 0.
  xall_std = 1.
  
  
  !Define variance in obs. No covariance in the observations
  !Eventually, will need to run a lot of orbits to estimate the true error covaraince matrix for pixels with 0 rain fraction and 0 land fraction.
  sy(:,:) = 0.
  do a=1,nobs1
    do b=1,nobs1
      !sy(a,b) = GMI_clr_errcov(gmi_ch(a),gmi_ch(b))
      sy(a,b) = LUT%eof_covmat_water(neof,gmi_ch(a),gmi_ch(b))
    end do
    sy(a,a) = sy(a,a)+GMI_NEDT(gmi_ch(a))**2
    !sy(a,a) = (0.0+GMI_NEDT(gmi_ch(a)))**2  !sensor NEDT + calibration stability
    !if(a .gt. 9) sy(a,a) = (3.+GMI_NEDT(GMI_nf(gmi_ch(a))))**2
    !print*, a, sqrt(sy(a,a))
  end do  
  do a=nobs1+1,nobs
    do b=nobs1+1,nobs
      !sy(a,b) = GMI_clr_errcov(gmi_ch(a),gmi_ch(b))
      sy(a,b) = LUT%eof_covmat_water(neof,gmi_ch(a),gmi_ch(b))
    end do
    sy(a,a) = sy(a,a)+GMI_NEDT(gmi_ch(a))**2
    !sy(a,a) = (0.0+GMI_NEDT(gmi_ch(a)))**2  !sensor NEDT + calibration stability
    !if(a .gt. 9) sy(a,a) = (3.+GMI_NEDT(GMI_nf(gmi_ch(a))))**2
    !print*, a, sqrt(sy(a,a))
  end do 
  !print '(29F6.2)', sy
  !stop
  call minvert(nobs,sy,sy_i)
      
  !Initialize profile parameters from ancillary data or climatological means (based on tskin)
  nlev = MERRA_NLEV
  do while(gdata%psfc .lt. MERRA_PRESSLEV(MERRA_NLEV-nlev+1))
    nlev=nlev-1
  end do
         
  !Intialize a priori
  xa=0.

  !Initialize observation vector 
  y(1:nobs1) = gdata%gmi_tb(gmi_ch(1:nobs1))
  !print*, gmi_ch(nobs1+1:nobs)
  !print*, y2locs(1:nobs2)
  do a=1,nobs2
    y(nobs1+a) = tb2d(gmi_ch(nobs1+a),y2locs(a))
  end do
  !print*, y
  !stop
  call minvert(nvar,sa,sa_i)
! c-----GENERATE FIRST GUESS----- (currently setting first guess to apriori)

  x = xa

! c-----NEWTONIAN ITERATION-----
  end_flag=0
  n=0
  do while(end_flag .eq. 0)
    n=n+1 
    !print*, 'Iteration ', n
    !print*, x
    !print*, xmax
    call get_state_clw2d_rh(nvar,neof,nclw,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,clwprofs,gdata)
    
    if(n .eq. 1) then ! get intial tpw and tskin uncertainty
      if(gdata%anc_tskin .gt. 0.) then
        call calc_tpw(gdata,nlev,ret%tpw_init)
        !ret%tskin_init = gdata%anc_tskin
      else
        ret%tpw_init = tpw_mean
        !ret%tskin_init = tskin_mean
      endif
      ret%tskin_init = x(neof+3)*xstd(neof+3)+xbar(neof+3)
      ret%stpw_init = tpw_rms
      !ret%stskin_init = tskin_rms
      ret%stskin_init = xstd(neof+3)
      ret%w10m_init = gdata%w10m
      ret%sw10m_init = xstd(neof+1)
      !ret%relAz_init = gdata%relAz
      ret%srelAz_init = xstd(neof+2)
      ret%clwp_init = exp(x(neof+8)*xstd(neof+8)+xbar(neof+8))
      ret%sclwp_init = xstd(neof+8)
      !print '(6F8.2)', ret%w10m_init, ret%clwp_init, ret%tpw_init, ret%tskin_init, tpw_mean, tskin_mean
      !print '(6F8.2)', ret%sw10m_init, ret%sclwp_init, ret%stpw_init, ret%stskin_init, tpw_rms, tskin_rms
    endif
    !print '(20F8.3)', x(1:neof), gdata%tskin, gdata%w10m, gdata%relAz, exp(x(neof+4)*xstd(neof+4)+xbar(neof+4))
    !Simulate Tbs for this column
    
    !print*, 
    !print*, gdata%temp_lev
    !print*, gdata%qv_lev
    !print*, gdata%hgt_lev
    !print '(29F8.2)', y
    do a = 1, nclw
      gdata%cloud_water = clwprofs(a,:)
      call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
      tb_out_2d(:,a) = tb_out
    end do
    do a=1, nobs1
      F(a) = sum(gmi_gain(gmi_ch(a),1:nclw)*tb_out_2d(gmi_ch(a),1:nclw))
    end do
    do a=1, nobs2
      F(nobs1+a) = tb_out_2d(gmi_ch(nobs1+a),y2locs(a))
    end do
    !print '(29F8.2)', F
   
    if(n .eq. 1) then
      ret%sim_tb0(gmi_ch(1:nobs1)) = tb_out(gmi_ch(1:nobs1))
    endif
    
    !Calculate Jacobian
    
    do i=1,neof+3
      !create perturbed atmosphere
      gdata_p = gdata
      xprime = x
      xprime(i) = x(i)+0.01
      call get_state_clw2d_rh(nvar,neof,nclw,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,clwprofs,gdata_p)
      do a=1,nclw
        gdata_p%cloud_water = clwprofs(clwlocs(a),:)
        call rtm_emission_csu(gdata_p,nlev,ch_mask,tb_out)
        tb_out_2dp(:,a) = tb_out
      end do
      do a=1, nobs1
        Fprime(a) = sum(gmi_gain(gmi_ch(a),1:nclw)*tb_out_2dp(gmi_ch(a),1:nclw))
      end do
      do a=1, nobs2
        Fprime(nobs1+a) = tb_out_2dp(gmi_ch(nobs1+a),y2locs(a))
      end do
      !print '(13F8.2)', Fprime
      do j=1,nobs
        K(j,i) = (Fprime(j)-F(j))/(xprime(i)-x(i))
      end do
      !print '(29F8.4)', K(:,i)
    end do
    do i=neof+4,nvar
      !create perturbed atmosphere
      gdata_p = gdata
      tb_out_2dp = tb_out_2d
      xprime = x
      xprime(i) = x(i)+0.01
      call get_state_clw2d_rh(nvar,neof,nclw,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,clwprofs,gdata_p)
      a = i-neof-3
      gdata_p%cloud_water = clwprofs(clwlocs(a),:)
      if(a .eq. 5) then
        call rtm_emission_csu(gdata_p,nlev,ch_mask,tb_out)
      else
        call rtm_emission_csu(gdata_p,nlev,ch_mask2,tb_out)
      endif
      tb_out_2dp(:,a) = tb_out
      
      do a=1, nobs1
        Fprime(a) = sum(gmi_gain(gmi_ch(a),1:nclw)*tb_out_2dp(gmi_ch(a),1:nclw))
      end do
      do a=1, nobs2
        Fprime(nobs1+a) = tb_out_2dp(gmi_ch(nobs1+a),y2locs(a))
      end do
      !Fprime = tb_out(gmi_ch(1:nobs))
      !print '(13F8.2)', Fprime
      do j=1,nobs
        K(j,i) = (Fprime(j)-F(j))/(xprime(i)-x(i))
      end do
      !print '(30F8.4)', K(:,i)
    end do
    !stop
    call oe_step(nobs,nvar,x,xa,sa_i,y,F,sy_i,K,xnext,sx,sx_i,prod2,xmin,xmax)
    !Calculate convergence criterion (normalized difference between steps)
    xdiff(1,:)=x-xnext
    prod6 = matmul(xdiff,sx_i)
    prod7 = matmul(prod6,transpose(xdiff))

    !Reset x for next step
    x = 0.25*x+0.75*xnext !"memory" factor decreases steps to convergence in oscillations
    !Check for exit conditions
    if (prod7(1,1) .lt. 0.2*(nvar+nobs)) end_flag=1
    if (n .ge. max_iter) end_flag=1
    
    if (end_flag .eq. 1) then 
      
      !Evaluate A matrix and chi-squared
      call oe_diagnostics(nvar,nobs,x,xa,sa_i,F,y,sx,sy_i,prod2,A_Matrix,Chi_Squared)
      !write(fmt, '(A1,I2,A5)') '(',neof+4,'F8.4)'
      !print*, 'Sigma X:'
      !print fmt, sx
      !print*, 'A Matrix:'
      !print fmt, A_Matrix
      !print*, 'Chi-Squared:',Chi_Squared
      
      
      !run forward model and calculate Tbs
      call get_state_clw2d_rh(nvar,neof,nclw,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,clwprofs,gdata)
      !print '(13F8.2)', y
      !simulate all Tbs, even if no used in retrieval
      !ch_mask(:) = .true.
      do a = 1, nclw
        gdata%cloud_water = clwprofs(a,:)
        call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
        tb_out_2d(:,a) = tb_out
      end do
      do a=1, nobs1
        F(a) = sum(gmi_gain(gmi_ch(a),1:nclw)*tb_out_2d(gmi_ch(a),1:nclw))
      end do
      ret%sim_tb(gmi_ch(1:nobs1)) = F(1:nobs1)
      
      gdata%cloud_water = clwprofs(5,:)
      call calcPIA(gdata,nlev,ret%PIAKu,ret%PIAKa)
      !call calc_tpw(gdata,nlev,ret%tpw_ret)
      !print*, nlev
      call calc_tpw(gdata,nlev,tpw_mean)
      ret%tpw_ret = tpw_mean
      
      !print*, ret%tpw_init, ret%tpw_ret
      !print '(13F8.2)', F
      !ret%sim_tb = tb_out
      
      !diagnostic output
!       do i=1,10
!         xall(1:nvar-3) = x(1:nvar-3)
!         gdata_p=gdata
!         do j=1,2+2*MERRA_NLEV
!           if(j .le. nvar-3) xall(j) = xall(j)+sqrt(sx(j,j))*r4_normal_01(seed)
!           if(j .gt. nvar-3) xall(j) = r4_normal_01(seed)
!           !print*, j, xall(j)
!         end do
!         !stop
!         call get_state_noclw(5+2*MERRA_NLEV,nlev,xall,xall_bar,xall_std,prof_avg,prof_eofs,gdata_p)
!         call calc_tpw(gdata_p,nlev,tpw(i))
!         !tskin(i) = gdata_p%tskin
!         print*, tpw(i), tskin(i)
!       end do
!       tpw_mean = sum(tpw(1:10))/10.
!       print*, tpw_mean
!       !tskin_mean = sum(tskin)/100.
!       tpw_rms = 0.
!       !tskin_rms = 0.
!       do i=1,10
!         tpw_rms = tpw_rms + (tpw(i)-tpw_mean)**2
!         !tskin_rms = tskin_rms + (tskin(i)-tskin_mean)**2
!       end do
!       tpw_rms = sqrt(tpw_rms/100.)
      !tskin_rms = sqrt(tskin_rms/100.)
      !ret%tpw_ret = tpw_mean
      !ret%tskin_ret = tskin_mean
      ret%tskin_ret = x(neof+3)*xstd(neof+3)+xbar(neof+3)
      ret%stpw_ret = missing_r4!tpw_rms
      !ret%stskin_ret = tskin_rms
      ret%stskin_ret = xstd(neof+3)*sqrt(sx(neof+3,neof+3))
      ret%atskin_ret = A_matrix(neof+3,neof+3)
      ret%w10m_ret = gdata%w10m
      ret%sw10m_ret = xstd(neof+1)*sqrt(sx(neof+1,neof+1))
      ret%aw10m = A_Matrix(neof+1,neof+1)
      ret%relAz_ret = gdata%relAz
      if(ret%relAz_ret .lt. 0.) ret%relAz_ret = -1.*ret%relAz_ret
      if(ret%relAz_ret .gt. 180.) ret%relAz_ret = 180.-ret%relAz_ret
      ret%srelAz_ret = xstd(neof+2)*sqrt(sx(neof+2,neof+2))
      ret%arelAz = A_Matrix(neof+2,neof+2)
      ret%clwp_init = gdata%anc_clwp
      ret%clwp_ret = exp(x(neof+8)*xstd(neof+8)+xbar(neof+8))
      ret%sclwp_ret = xstd(neof+8)*sqrt(sx(neof+8,neof+8))
      ret%aclwp = A_Matrix(neof+8,neof+8)
      !print*, ret%w10m_ret, ret%clwp_ret, ret%tpw_ret, ret%tskin_ret
      !print*, ret%sw10m_ret, ret%sclwp_ret, ret%stpw_ret, ret%stskin_ret
      !print*, ret%aw10m, ret%aclwp
      ret%niter = n
      ret%chi_squared = Chi_Squared/(nvar+nobs)
      ret%obs_err = 0.
      do i=1,nobs
        ret%obs_err = ret%obs_err+(F(i)-y(i))**2/(sy(i,i))
      end do
      ret%obs_err = ret%obs_err/nobs
      ret%x_err = 0.
      do i=1,nvar
        ret%x_err = ret%x_err+(x(i)-xa(i))**2/(sa(i,i))
      end do
      ret%x_err = ret%x_err/nvar
      ret%emis = missing_r4
      ret%aemis = missing_r4
      ret%semis = missing_r4
      do i=1,13
        if(ch_mask(i)) ret%emis(i) = gdata%emis(i)
      end do
      
      ret%err_eof = missing_r4
      ret%serr_eof = missing_r4
      ret%aerr_eof = missing_r4
      do i=1,min(10,neof)
        ret%err_eof(i) = x(i)*xstd(i)+xbar(i)
        ret%serr_eof(i) = sqrt(Sx(i,i))
        ret%aerr_eof(i) = A_Matrix(i,i)
      end do
      !stop
      do i =1,MERRA_NLEV
      	if(MERRA_PRESSLEV(i) .eq. gdata%press_lev(1)) exit
      end do
      ret%press_lev = MERRA_PRESSLEV
      ret%temp_lev(1:i) = gdata%temp_lev(1)
      ret%temp_lev(i:MERRA_NLEV) = gdata%temp_lev(1:MERRA_NLEV-i+1)
      ret%hgt_lev(i:MERRA_NLEV) = gdata%hgt_lev(1:MERRA_NLEV-i+1)
      ret%qv_lev(1:i) = gdata%qv_lev(1)
      ret%qv_lev(i:MERRA_NLEV) = gdata%qv_lev(1:MERRA_NLEV-i+1)
      ret%cloud_water(i:MERRA_NLEV) = gdata%cloud_water(1:MERRA_NLEV-i+1)
      ret%t2m_ret = gdata%t2m

      exit
    endif	!end final iteration
  end do	!end Newtonian loop
  
  deallocate(x, xa, xprime, xnext, xdiff)
  deallocate(xall, xall_std, xall_bar)
  deallocate(xmin, xmax, xstd, xbar)
  deallocate(sa, sa_i)
  deallocate(sy,sy_i)
  deallocate(K, K_t)
  deallocate(sx, sx_i)
  deallocate(A_Matrix)
  deallocate(F, Fprime)!, F_a(9), Fout(9)
  deallocate(y)
  deallocate(prod2,prod6)

end subroutine gmi_ocean_ret_clw2d

end module ocean_ret_raob
