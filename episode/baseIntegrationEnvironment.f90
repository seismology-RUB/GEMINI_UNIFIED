! ===========================================================================
!  Base class for handling integration environment
! ===========================================================================
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
!-----------------------------------------------------------------------------
!  Defines an abstract base class for an integration environment to be extended
!  by user-provided concrete realizations. These must provide procedures that
!  return the location of the nodes and discontinuities inside the open integration
!  range. Abstract interfaces for these procedures are defined here. It also 
!  provides procedures to set or increase the "layer" index, i.e.
!  the index of the range where the derivatives are continuous. The latter two
!  are implemented here. Material parameters are not dealt with here. This
!  should be done in the concrete implementations of the ODE procedures.
!------------------------------------------------------------------------
module baseIntegrationEnvironment
    use errorMessage
    implicit none
    type, abstract :: base_integration_environment
       integer :: nl                                     ! index of current layer between discontinuitites
    contains
        procedure :: setInterval
        procedure :: nextInterval
        procedure (getDiscontinuitiesIntegrationEnvironment), deferred :: getDiscontinuities
        procedure (getNodesIntegrationEnvironment), deferred :: getNodes
    end type base_integration_environment
!---------------------------------------------------
!  Interface for procedure returning discontinuities 
!  inside open integration range
!  x1:     start of integration interval (in)
!  x2:     end of integration interval   (in)
!  jl:     global index of first disco in interval (out)
!  jr:     global index of last disco in interval  (out)
!  xd:     array with location of discontinuities inside open interval (out)
!  errmsg: error message object (inout)
!
! Note: this procedure must return the discontinuities in
! - ascending order  if x2 > x1: ---x1---xd(1)---------------xd(n)----x2---->
! - descending order if x2 < x1: ---x2---xd(n)---------------xd(1)----x1---->
!
    abstract interface
       subroutine getDiscontinuitiesIntegrationEnvironment(this,x1,x2,jl,jr,xd,errmsg)
       import base_integration_environment
       import error_message
       class (base_integration_environment) :: this
       double precision :: x1,x2
       integer :: jl,jr
       double precision, dimension(:), allocatable :: xd
       type (error_message) :: errmsg
       end subroutine getDiscontinuitiesIntegrationEnvironment
    end interface
!---------------------------------------------------
!  Interface for procedure returning nodes
!  inside closed integration range
!  x1:     start of integration interval (in)
!  x2:     end of integration interval   (in)
!  jl:     global index of first node in interval (out)
!  jr:     global index of last node in interval  (out)
!  xn:     array with location of nodes inside open interval (out)
!  errmsg: error message object (inout)
!
! Note: this procedure must return the nodes in 
! - ascending order  if x2 > x1: ---x1---xn(1)---------------xn(n)----x2---->
! - descending order if x2 < x1: ---x2---xn(n)---------------xn(1)----x1---->
! Nodes may coincide with discontinuities
!
    abstract interface
       subroutine getNodesIntegrationEnvironment(this,x1,x2,jl,jr,xn,errmsg)
       import base_integration_environment
       import error_message
       class (base_integration_environment) :: this
       double precision :: x1,x2
       integer :: jl,jr
       double precision, dimension(:), allocatable :: xn
       type (error_message) :: errmsg
       end subroutine getNodesIntegrationEnvironment
    end interface
!
contains
!---------------------------------------------------------
!> \brief goto next interval, discontinuity was reached
!
    subroutine nextInterval(this,idir)
    class (base_integration_environment) :: this
    integer :: idir
    this%nl = this%nl+idir
    end subroutine nextInterval
!---------------------------------------------------------
!> \brief goto next interval, discontinuity was reached
!
    subroutine setInterval(this,nl)
    class (base_integration_environment) :: this
    integer :: nl
    this%nl = nl
    end subroutine setInterval
end module baseIntegrationEnvironment
