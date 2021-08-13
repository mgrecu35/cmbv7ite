subroutine check_nan(hgt_lev,temp_lev)
real :: hgt_lev(27),temp_lev(27)
integer :: i
do i=1,27
   if(isnan(hgt_lev(i))) then
      print*, 'bingo'
   end if
   if(isnan(temp_lev(i))) then
      print*, 'bingo'
   end if
enddo

end subroutine check_nan
