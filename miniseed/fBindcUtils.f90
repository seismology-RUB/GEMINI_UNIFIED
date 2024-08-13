!-------------------------------------------------------------
!  Utilities for calling C from Fortran
!
 module fBindcUtils
       use iso_c_binding
	implicit none
!
 contains
!--------------------------------------------------------------------------------------
!  convert a string to a c_char array
!  cut string if array is too small
!  append c_null_char if array is too large
!
	subroutine convertStringCharFBindcUtils(string,array)
	character (len=*) :: string
	character (kind = c_char), dimension(:) :: array
	integer :: n,i
!
	n = min(size(array)-1,len_trim(string))    ! leave space for \0 at end of array
	if (size(array) < len_trim(string)+1) then
		print *,'Warning'
		print *,'convertStringCharFBindcUtils: character array too short to take ',trim(string)
		print *,'Some characters will be lost.'
	endif
	array = c_null_char                        ! fill with \0
	do i=1,n
		array(i) = string(i:i)
	enddo
	end subroutine convertStringCharFBindcUtils
!-----------------------------------------------------------------------
!  convert a c_char array to a string
!  fill string up to its length
!
	subroutine convertCharStringFBindcUtils(array,string,ls,warn)
	character (kind = c_char), dimension(:) :: array
	character (len=*) :: string
	integer :: i,ls
	logical :: warn
!
!  fill string with blanks
!
	do i=1,len(string); string(i:i) = ' '; enddo
!
	if (warn .and. ls > len(string)) then
		print *,'Warning'
		print *,'convertCharStringFBindcUtils: character array too large for string'
		print *,'Some characters will be lost.'
		print *,ls,len(string)
	endif
	do i = 1,ls
		string(i:i) = array(i)
	enddo
	end subroutine convertCharStringFBindcUtils
!
 end module
