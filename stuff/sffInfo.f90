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
!--------------------------------------------------------
!  module to describe a SFF info line
!--------------------------------------------------------
 module sffInfo
    implicit none
    private constructSFFInfo
    interface new
        module procedure constructSFFInfo
        module procedure setSFFInfo
    end interface
    type sff_info
        private
        character (len=1) :: cs        ! coordinate system
        real :: c1, c2, c3             ! location of info
        integer :: nstack              ! number of stacks
    end type sff_info
!
 contains
!-----------------------------------------------------------
!  constructor
!
    subroutine constructSFFInfo(this,lu)
    type (sff_info) :: this
    integer :: lu,ierr
!
    call sff_RInfo(lu,this%cs,this%c1,this%c2,this%c3,this%nstack,ierr)
    end subroutine constructSFFInfo
!---------------------------------------------------------
!  second constructor
!
    subroutine setSFFInfo(this,cs,c1,c2,c3,nstack)
    type (sff_info) :: this
    character (len=1) :: cs
    real :: c1,c2,c3
    integer :: nstack
!
    this = sff_info(cs,c1,c2,c3,nstack)
    end subroutine setSFFInfo
!-----------------------------------------------------------
!  return receiver location
!
    subroutine locationSFFInfo(this,x,y,z)
    type (sff_info) :: this
    real :: x,y,z
    x = this%c1; y = this%c2; z = this%c3
    end subroutine locationSFFInfo
!-----------------------------------------------------------
!  print info line
!
    subroutine printSFFInfo(this)
    type (sff_info) :: this
    print *,this%cs,this%c1,this%c2,this%c3,this%nstack
    end subroutine printSFFInfo
!---------------------------------------------------------
!  write info line to file
!
    subroutine writeSFFInfo(this,lu)
    type (sff_info) :: this
    integer :: lu
!
    call sff_WInfo(lu,this%cs,this%c1,this%c2,this%c3,this%nstack)
    end subroutine writeSFFInfo
!---------------------------------------------------------
 end module sffInfo
