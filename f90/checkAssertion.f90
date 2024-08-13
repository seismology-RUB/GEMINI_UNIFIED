! ======================================================================
!  Check assertions set up in unit testing
! ======================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------
!> \brief Module for checking assertions made in unit testing
!------------------------------------------------------------------
module checkAssertion
    use errorMessage
    use string
    implicit none
    interface genericCheckAssertion
       module procedure realCheckAssertion
       module procedure integerCheckAssertion
    end interface
contains
    !---------------------------------------------------------------
    ! check if two real values agree
    !
    subroutine realCheckAssertion(description,val,refval,eps,res)
        character (len=*) :: description
        real :: val,refval,eps
        logical :: res
        !
        if (abs(val-refval) < eps) then
           res = .true.
        else
           res = .false.
           print *,'The assertion: ',trim(description),' failed!'
           print *,'Tested value: ',val,', reference value: ',refval
        endif
    end subroutine realCheckAssertion
    !---------------------------------------------------------------
    ! check if two real values agree
    !
    subroutine integerCheckAssertion(description,val,refval,res)
        character (len=*) :: description
        integer :: val,refval
        logical :: res
        !
        if (val == refval) then
           res = .true.
        else
           res = .false.
           print *,'The assertion: ',trim(description),' failed!'
           print *,'Tested value: ',val,', reference value: ',refval
        endif
    end subroutine integerCheckAssertion
end module checkAssertion
