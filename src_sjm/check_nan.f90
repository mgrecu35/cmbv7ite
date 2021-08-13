subroutine check_nan(hgt_lev,temp_lev, istatus)
real :: hgt_lev(27),temp_lev(27)
integer :: i
istatus=0
do i=1,27
   if(isnan(hgt_lev(i))) then
      !print*, 'bingo'
      istatus=1
   end if
   if(isnan(temp_lev(i))) then
      !print*, 'bingo'
      istatus=1
   end if
enddo

end subroutine check_nan
