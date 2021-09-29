SUBROUTINE gmi_land_ret(gdata, ret) !retrieval variables

  use profile_def
  use LUT_def
  use missingMod
  
  implicit none
  
  include 'constants_sjm.inc'
  include 'parametersMERRA.inc'
  include 'parametersGMI.inc'
  
  real*4  r4_normal_01
  
  type(profile), intent(inout) :: gdata
  type(profile) :: gdata_p
  type(ret_data) :: ret
  
  
  
  real :: pc_std(1+2*MERRA_NLEV), prof_std(1+2*MERRA_NLEV), prof_W(1+2*MERRA_NLEV), prof_avg(1+2*MERRA_NLEV)
  real :: prof_eofs(1+2*MERRA_NLEV,1+2*MERRA_NLEV)
  real :: emis_corrmat(LUT%sfc_nch,LUT%sfc_nch)
  
  character*7 :: fmt

  integer :: neof, nvar, nobs, seed, nemis
  logical :: ch_mask(13)
  integer :: gmi_ch(13), emis_ch(11)
!   !parameter (nvar=3)
! 
  real :: tb_out(13)
  real :: tpw(100), tskin(100), tpw_mean, tskin_mean, tpw_rms, tskin_rms
  real :: emis0(LUT%sfc_nch), estd0(LUT%sfc_nch)
!   ! Matrices to be used internally in optimal estimation expression
! 
  integer :: i, j, a, b, n, end_flag, max_iter, nlev
  parameter (max_iter=9)
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
  
  if(gdata%anc_psfc .gt. 0.) gdata%psfc = gdata%anc_psfc
  !print*, gdata%tskin_index, gdata%sfc_elev, gdata%psfc, gdata%lat, gdata%lon
  !stop
  !if sfc pressure not present, extrapolate from tskin and sfc elev assuming 1013mb at sfc
  if(gdata%psfc .le. 0.) gdata%psfc = 1013.*exp(-1000.*gdata%sfc_elev*9.81/(float(gdata%tskin_index)*287.))
  !print*, gdata%sfc_elev, gdata%psfc
  !read in the EOFs for this sfc type, sfc pressure, and tskin
  !call get_prof_eofs(gdata, prof_W, prof_avg, prof_std, prof_eofs, tpw_mean, tpw_rms, tskin_mean, tskin_rms)
  prof_W = LUT%prof_W(2,:)
  prof_std = LUT%prof_std(2,:)
  pc_std = LUT%pc_std(2,:)
  prof_eofs = LUT%prof_eofs(2,:,:)
  prof_avg = 0.
   
  neof = 5
  !print*, nvar
  
!   do i=1,neof
!    print*, i, pc_std(i)
!    print '(28F6.2)', prof_eofs(i,:)*pc_std(i)*prof_std
!   end do
!   stop
  gmi_ch = 0
  !count channels with valid data
  nobs = 0
  ch_mask = .false.
  do i=1,13
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
  
  !count number of emissivity variables to retrieve based on channel availability
  nemis=0
  do i=1,4 !10-18
    if(ch_mask(i)) then
      nemis=nemis+1
      emis_ch(nemis) = i
    endif
  end do
  do i=6,9 !36-89
    if(ch_mask(i)) then
      nemis=nemis+1
      emis_ch(nemis) = i
    endif
  end do
  !special case - need 18V and 36V to retrieve 23V
  if(ch_mask(3) .eq. .false. .and. ch_mask(6) .eq. .false.) then
    print*, 'Need 18 and 36 to retrieve 23'
    return
  endif
  
  if(ch_mask(10) .or. ch_mask(12) .or. ch_mask(13)) then
    nemis=nemis+1 !166V, 183 all use same emis
    emis_ch(nemis) = 10
  endif
  if(ch_mask(11)) then
    nemis=nemis+1 !166H
    emis_ch(nemis) = 11
  endif
  
  nvar=neof+nemis+1!special variable that controls weight of 18v and 36v to determine 23v
  !print*, nemis, neof, nvar, nobs, gdata%sfc_type
  
 
  allocate(x(nvar), xa(nvar), xprime(nvar), xnext(nvar), xdiff(1,nvar))
  !allocate(xall(nemis+1+2*MERRA_NLEV), xall_std(nemis+1+2*MERRA_NLEV), xall_bar(nemis+1+2*MERRA_NLEV))
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
  !To do: replace emissivity elements with emissivity correlation matrix
  !call get_emis_corrmat(gdata%sfc_type, emis_corrmat)
  
  

  !Standard Deviation of Retrieval Parameters
  xstd(1:neof) = pc_std(1:neof) !from raob error analysis
  !if(gdata%sfc_type .ge. 13) xstd(1:nvar-nemis) = 0.1 !even stronger constraint over coast)
  !xstd(nvar-nemis+1:nvar) = 0.5 !very broad emissivity allowance

  !Mean value for retrieval parameters
  xbar(1:nvar-(nemis+1)) = 0. !standardize emissivity
  !xbar(nvar-nemis+1:nvar) = 0.8 !emissivity
  
  !print*, gdata%sfc_type
  !print*, gdata%lat, gdata%lon, gdata%sfc_type
  call get_gmi_emis_std(gdata%lat, gdata%lon, gdata%sfc_type, emis0, estd0)
  do  i=1,nemis
    !print*, emis_ch(i), GMI_default_emis_std(emis_ch(i),gdata%sfc_type-1)
    !print*, i, emis_ch(i), emis0(emis_ch(i)), estd0(emis_ch(i))
    if(emis0(emis_ch(i)) .le. 0) then
      emis0(emis_ch(i)) = 0.75
      estd0(emis_ch(i)) = 0.25
    endif
    xbar(nvar-(nemis+1)+i) = emis0(emis_ch(i))!GMI_V4_emis(emis_ch(i),gdata%sfc_type-1)
    xstd(nvar-(nemis+1)+i) = estd0(emis_ch(i))!0.25!4.*GMI_V4_emis_std(emis_ch(i),gdata%sfc_type-1)
  end do
  
  !limits for emis
  xmin(nvar-(nemis+1)+1:nvar-1) = (0.25-xbar(nvar-(nemis+1)+1:nvar))/xstd(nvar-(nemis+1)+1:nvar) !lowest value for water
  xmax(nvar-(nemis+1)+1:nvar-1) = (1.5-xbar(nvar-(nemis+1)+1:nvar))/xstd(nvar-(nemis+1)+1:nvar) !higher than 1 for snow cover/RFI

  !emis weight variable for 23v
  xbar(nvar) = 0.3
  xstd(nvar) = 0.2
  xmin(nvar) = (0.0-xbar(nvar))/xstd(nvar)
  xmax(nvar) = (1.0-xbar(nvar))/xstd(nvar)
  !cloud water variable
  !xbar(nvar) = alog(0.03) !CLWP (retrieved in log space)  (kg/m^2)
  !xstd(nvar) = 2.
  !limits for wind and clw
  !xmin(nvar) = (alog(0.001)-xbar(nvar))/xstd(nvar)
  !xmax(nvar) = (alog(2.0)-xbar(nvar))/xstd(nvar)

  !set xall for tpw and tskin error comps.
  !xall = 0.
  !xall_bar = 0.
  !xall_std = 1.
  
  !Set observation + model error covariance matrix (Sy).
  !This is based on the number of EOFs
  !Add sensor noise to diagonal elements
  
  !Eventually, will need to run a lot of orbits to estimate the true error covaraince matrix for pixels with 0 rain fraction and 0 land fraction.
  sy(:,:) = 0.

  do a=1,nobs
    do b=1,nobs
      sy(a,b) = LUT%eof_covmat_land(neof,gmi_ch(a),gmi_ch(b))
    end do
    
    sy(a,a) = sy(a,a)+GMI_NEDT(gmi_ch(a))**2  !CRTM_TMI_RMS(a))**2!sensor error + model_error
    !if(a .gt. 9) sy(a,a) = (3.+GMI_NEDT(GMI_nf(gmi_ch(a))))**2
    !print*, a, sqrt(sy(a,a))
    !print '(13F8.4)', sy(a,:)
  end do  
  !stop

  call minvert(nobs,sy,sy_i)
  
  !Initialize profile parameters from ancillary data or climatological means (based on tskin)
  nlev = MERRA_NLEV
  do while(gdata%psfc .lt. MERRA_PRESSLEV(MERRA_NLEV-nlev+1))
    nlev=nlev-1
  end do
  !nlev = 0
  !do while(gdata%anc_temp_lev(i) .gt. 0 .and. gdata%anc_qv_lev(i) .ge. 0)
  !  nlev=nlev+1
  !end do
  !nlev = MERRA_NLEV-i+1
  !print*, nlev
  !stop
  !Intialize a priori
  xa=0.

  !Initialize observation vector 
  do i=1,nobs
    y(i) = gdata%gmi_tb(gmi_ch(i))
  end do

  call minvert(nvar,sa,sa_i)
! c-----GENERATE FIRST GUESS----- (currently setting first guess to apriori)

  x = xa

! c-----NEWTONIAN ITERATION-----
  end_flag=0
  n=0
  do while(end_flag .eq. 0)
    n=n+1 
    !print*, 'Iteration ', n
    !print*, nvar, nemis
    !print*, x(1:neof)
    !stop
    call get_state_land_rh(nvar,neof,nemis,emis_ch,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)
    !print '(13F8.4)', gdata%emis
    
    if(n .eq. 1) then ! get intial tpw and tskin uncertainty - this should be built into the LUTs
      if(gdata%anc_tskin .gt. 0.) then
        call calc_tpw(gdata,nlev,ret%tpw_init)
        ret%tskin_init = gdata%anc_tskin
      else
        ret%tpw_init = tpw_mean
        ret%tskin_init = tskin_mean
      endif
      ret%stpw_init = tpw_rms
      ret%stskin_init = tskin_rms
      !print*, ret%tpw_init, ret%tskin_init, ret%stpw_init, ret%stskin_init
      ret%w10m_init = -9999.9
      ret%sw10m_init = -9999.9
      ret%clwp_init = -9999.9
      ret%sclwp_init = -9999.9
      !print*, ret%tpw_init, ret%tskin_init
      !print*, ret%stpw_init, ret%stskin_init
    endif
    !print*, gdata%emis
    !Simulate Tbs for this column
    !print '(13F8.2)', y
    !call rtm_emission_sc(gdata,nlev,ch_mask,tb_out)
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
      xprime(i) = x(i)+0.1
      call get_state_land_rh(nvar,neof,nemis,emis_ch,nlev,xprime,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata_p)
      call rtm_emission_csu(gdata_p,nlev,ch_mask,tb_out)
      Fprime = tb_out(gmi_ch(1:nobs))
      !print*, Fprime
      do j=1,nobs
        K(j,i) = (Fprime(j)-F(j))/(xprime(i)-x(i))
      end do
      !print '(13F8.4)', K(:,i)
    end do
    !stop
    !return
    call oe_step(nobs,nvar,x,xa,sa_i,y,F,sy_i,K,xnext,sx,sx_i,prod2,xmin,xmax)
    
    !Calculate convergence criterion (normalized difference between steps)
    xdiff(1,:)=x-xnext
    prod6 = matmul(xdiff,sx_i)
    prod7 = matmul(prod6,transpose(xdiff))

    !Reset x for next step
    x = 0.25*x+0.75*xnext !25% "memory" factor decreases steps to convergence in oscillations
    !Check for exit conditions
    if (prod7(1,1) .lt. 0.2*(nvar+nobs)) end_flag=1
    if (n .gt. max_iter) end_flag=1
    !print*, prod7(1,1), 0.2*(nvar+nobs)
    if (end_flag .eq. 1) then 
      !Evaluate A matrix and chi-squared
      call oe_diagnostics(nvar,nobs,x,xa,sa_i,F,y,sx,sy_i,prod2,A_Matrix,Chi_Squared)
      !print*, 'Sigma X:'
      !write(fmt, '(A1,I1.1,A5)') '(',neof,'F8.4)'
      !print fmt, sx(1:neof,1:neof)
      !print*, 'A Matrix:'
      !print fmt, A_Matrix(1:neof,1:neof)
      !print*, 'Chi-Squared:',Chi_Squared
     
      
      !run forward model and calculate Tbs
      call get_state_land_rh(nvar,neof,nemis,emis_ch,nlev,x,xbar,xstd,prof_W,prof_avg,prof_std,prof_eofs,gdata)
      !print '(13F8.4)', gdata%emis
      !print '(13F8.2)', y
      call rtm_emission_csu(gdata,nlev,ch_mask,tb_out)
      F = tb_out(gmi_ch(1:nobs))
      !print '(13F8.2)', F
      call calc_tpw(gdata,nlev,tpw_mean)
      ret%tpw_ret = tpw_mean
      !print*, ret%tpw_init, ret%tpw_ret 
      !stop
      call calcPIA(gdata,nlev,ret%PIAKu,ret%PIAKa)
      !print*, ret%PIAKu, ret%PIAKa
      !ret%PIAKu = -99.
      !ret%PIAKa = -99.

      ret%sim_tb(gmi_ch(1:nobs)) = tb_out(gmi_ch(1:nobs))
      
      !diagnostic output
      do i=1,10
        xprime = x
        gdata_p=gdata
        do j=1,neof
          xprime(j) = xprime(j)+sqrt(sx(j,j))*r4_normal_01(seed)
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
      ret%stpw_ret = sqrt(tpw_rms/10.)

      ret%tskin_ret = gdata%anc_tskin!tskin_mean

      ret%clwp_init = gdata%anc_clwp

      ret%specularity=x(nvar)*xstd(nvar)+xbar(nvar)
      !print*, ret%tpw_ret, ret%tskin_ret
      !print*, ret%stpw_ret, ret%stskin_ret
      !print*, ret%specularity
      ret%emis = gdata%emis
      
      do i=1,min(10,neof)
        ret%err_eof(i) = x(i)*xstd(i)+xbar(i)
        ret%serr_eof(i) = sqrt(Sx(i,i))
        ret%aerr_eof(i) = A_Matrix(i,i)
      end do
      do i=1,nemis
        ret%semis(emis_ch(i)) = sqrt(sx(nvar-(nemis+1)+i,nvar-(nemis+1)+i))*xstd(nvar-(nemis+1)+i)
        ret%aemis(emis_ch(i)) = A_Matrix(nvar-(nemis+1)+i,nvar-(nemis+1)+i) 
      end do
      
      ret%semis(5) = sqrt(sx(nvar,nvar))!0.286517*ret%semis(6)+(1.-0.286517)*ret%semis(3) !interpolate 23V
      ret%aemis(5) = A_Matrix(nvar,nvar)!0.286517*ret%aemis(6)+(1.-0.286517)*ret%aemis(3) !interpolate 23V
      ret%semis(12:13) = ret%semis(10) !183 +/-3 and 183 +/-7 same as 166V
      ret%aemis(12:13) = ret%aemis(10)
      !print '(13F8.4)', ret%emis
      !print '(13F8.4)', ret%semis
      !print '(13F8.4)', ret%aemis
      
      ret%niter = n
      ret%chi_squared = Chi_Squared/(nvar+nobs)
      !print*, ret%chi_squared
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

end

