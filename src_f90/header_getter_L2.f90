!=============================================================================
!
!  NAME : header_getter    FORTRAN
!
!  PURPOSE: This is a quickie utility to retrieve the orbit number and date
!           number from an HDF5 formatted PPS file.  This information can
!           also be retrieved using Toolkit utilities; however, this is not
!           always convenient.  Note that if you don't give it an HDF5 file, 
!           it will give "H5F_OPEN" message and stop.
!
!  INPUT PARAMS:
!        input_file - HDF formatted input file
!
!  OUTPUT VALUE:
!        orbit_number  - orbit number in integer format
!        date_number   - date number in integer format, yyyymmdd
!        algorithm_ID  - type file (e.g. 1CGMI, 2AKu, etc)
!        error_code    - error output, integer; anything less than 0 ==> bad
!                        and you get a message
!
!  AUTHOR:  S.McLaughlin   03/18/2013
!
!  HISTORY: 
!    Date       Author        Change activity
!  03/18/2013 S.McLaughlin -  initial version
!  04/03/2013 S.McLaughlin -  added algorithm_ID
!  12/06/2013 S.McLaughlin -  adapted for use in L2 code (file name passing)
!
!=============================================================================

SUBROUTINE header_getter_L2(len_file, filename, orbit_number, date_number,     &
                            algorithm_ID, error_code)

USE tkio_f90  ! module contains necessary modules

IMPLICIT NONE

!...Calling sequence variables
    INTEGER*4, intent(in)   :: len_file
    CHARACTER(len=len_file) :: filename
    INTEGER*4, intent(out)  :: orbit_number
    INTEGER*4, intent(out)  :: date_number
    CHARACTER(len=10), intent(out) :: algorithm_ID
    INTEGER*4, intent(out)  :: error_code

!...Local computational variables
    INTEGER*4 :: place1, place3, place2, place_length
    INTEGER*4 :: year, month, day
    INTEGER*4 :: ret
    CHARACTER(len=10) :: dummy_jobid
    CHARACTER(len=255) :: date_string

    type (TKINFO), pointer  :: granuleHandleL2 => NULL()
   

    ret = TKgetAlgorithmVersion(dummy_jobid, filename, algorithm_ID)

    print*, algorithm_ID

    allocate(granuleHandleL2)
    ret = TKopen(filename,algorithm_ID,TKREAD, "HDF5", dummy_jobid,granuleHandleL2,1)
 
    ret = TKgetMetaInt(granuleHandleL2, 'FileHeader','GranuleNumber',orbit_number)
    print *, 'Granule Number ', orbit_number

    ret = TKgetMetaString(granuleHandleL2, 'FileHeader', 'StartGranuleDateTime',date_string)
    print *, 'StartGranuleDateTime', date_string

    READ(UNIT=date_string,FMT=53) year, month, day
    date_number = 10000*year + 100*month + day
53  FORMAT(i4,1x,i2,1x,i2)

    return
    end

    














