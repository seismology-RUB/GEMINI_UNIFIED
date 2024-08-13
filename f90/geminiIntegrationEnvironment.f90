! ======================================================================
!  Integration environment for GEMINI integrations
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
!-------------------------------------------------------------------------------
!  Implements an integration environment for use with GEMINI integrations.
!  Assumes a given earth model specifying discontinuities and material
!  parameters and externally defined radial nodes for solution storage.
!  Provides location of nodes and discontinuities to integration routines.
!  Always knows current layer index and is authorized to change it.
!  Access to material parameters is done via the ODE objects where they
!  are needed.
!------------------------------------------------------------------------
module geminiIntegrationEnvironment
    use baseIntegrationEnvironment
    use nodeEarthmodel
    use externalRadialNodes
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocGemini; end interface
    type, extends(base_integration_environment) ::  gemini_integration_environment
       double precision :: omre,omim                ! real and imaginary part of angular frequency
       double precision :: dll1                     ! harmonic degree (l*(l+1))
       type (node_earthmodel) :: nem                ! earth model object
       type (external_radial_nodes) :: exnod        ! external nodes object
    contains
        procedure :: createGemini
        procedure :: deallocGemini
        procedure :: getDiscontinuities => getDiscontinuitiesGemini
        procedure :: getNodes => getNodesGemini
    end type
!
contains
!--------------------------------------------------------
!> \brief Create integration environment
!! wn:     current wavenumber (in)
!! omre:   real part of current angular frequency (in)
!! omim:   imaginary part of current angular frequency (in)
!! nem:    node earth model object (in)
!! exnod:  external nodes object (in)
!! errmsg: error message object (in)
!
    subroutine createGemini(this,wn,omre,omim,nem,exnod,errmsg)
    class (gemini_integration_environment) :: this
    double precision :: omre,omim
    double precision :: wn
    type (node_earthmodel) :: nem
    type (external_radial_nodes) :: exnod
    type (error_message) :: errmsg
!
!  nl is set later dynamically, here only dummy value
!
    call addTrace(errmsg,'createGemini')
    this%nl = 1                            ! set layer index to 1 before more is known
    this%dll1 = (wn*(.rearth.nem))**2      ! l*(l+1) = (k*R)**2
    this%omre = omre; this%omim = omim
    this%nem = nem
    this%exnod = exnod
    end subroutine createGemini
!---------------------------------------------------------
!> \brief Deallocate (nothing to do)
!
    subroutine deallocGemini(this)
    class (gemini_integration_environment) :: this
    this%nl = 1                        ! dummy to avoid compiler warning
    end subroutine deallocGemini
!--------------------------------------------------------------------------
!> \brief Get array of discontinuities within open integration range
!! - ascending order  if x2 > x1: ---x1---xd(1)---------------xd(n)----x2---->
!! - descending order if x2 < x1: ---x2---xd(n)---------------xd(1)----x1---->
!! exclude discontinuities that coincide with integration limits
!! x1:     begin of integration range (in)
!! x2:     end of integration range (in)
!! jl:     index of disco right of or on xl = min(x1,x2)
!! jr:     index of disco left  of or on xr = max(x1,x2)
!! xd:     array of location of discos ordered as described above,
!!         not allocated if there are no discos in the range (out)
!! errmsg: error message objects (inout)
!
    subroutine getDiscontinuitiesGemini(this,x1,x2,jl,jr,xd,errmsg)
    class (gemini_integration_environment) :: this
    double precision :: x1,x2
    double precision, dimension(:), allocatable :: xd
    type (error_message) :: errmsg
    character(len=24) :: myname = 'getDiscontinuitiesGemini'
    double precision, dimension(:), pointer :: rb
    double precision :: xl,xr
    integer :: jl,jr,nd
!
    call addTrace(errmsg,myname)
    if (x1 > x2) then; xl = x2; xr = x1; else; xl = x1; xr = x2; endif
!
    rb => getRbArrayNodeEarthmodel(this%nem)
!
!  Handle both cases:
!
    jl = locate(xl,size(rb),rb)+1                          ! index of rb left of or on xl
    if (jl > size(rb)) then
       call add(errmsg,1,'all discos left of integration interval',myname)
       jr = jl-1                                           ! make jr-jl < 0, signifying no discos in interval
       return
    else
       if (dabs(1.d0-xl/rb(jl)) < epsilon(1.d0)) jl = jl+1     ! xl is on rb(jl), take next
       if (jl >= size(rb)) then
          call add(errmsg,1,'all discos left of integration interval',myname)
          jr = jl-1                                           ! make jr-jl < 0
          return
       endif
    endif
    jr = locate(xr,size(rb),rb)                            ! index of rb left of xr
    if (jr == 0) then
       call add(errmsg,1,'all discos right of integration interval',myname)
       jr = 0; jl = 1                                       ! make jr-jl < 0
       return
    endif
!
!  assign discos within open integration interval
!
    nd = jr-jl+1
    if (nd > 0) then
       allocate(xd(nd))
       if (x2 >= x1 ) then
          xd = rb(jl:jr)
       else
          xd = rb(jr:jl:-1)                               ! reverse order if x2 < x1
       endif
    endif
    end subroutine getDiscontinuitiesGemini
!----------------------------------------------------------------
!> \brief Get array storage nodes within open integration range
!! - ascending order  if x2 > x1: ---x1---xn(1)---------------xn(n)----x2---->
!! - descending order if x2 < x1: ---x2---xn(n)---------------xn(1)----x1---->
!! include nodes that coincide with integration limits
!! x1:     begin of integration range (in)
!! x2:     end of integration range (in)
!! jl:     index of node right of or on xl = min(x1,x2)
!! jr:     index of node left  of or on xr = max(x1,x2)
!! xn:     array of location of nodes ordered as described above,
!!         not allocated if there are no nodes in closed range (out)
!! errmsg: error message objects (inout)
!
    subroutine getNodesGemini(this,x1,x2,jl,jr,xn,errmsg)
    class (gemini_integration_environment) :: this
    double precision :: x1,x2
    type (error_message) :: errmsg
    double precision, dimension(:), allocatable :: xn
    double precision, dimension(:), pointer :: rnod
    double precision :: xl,xr
    integer :: jl,jr,nnod
    character(len=14) :: myname = 'getNodesGemini'
!
    call addTrace(errmsg,myname)
    if (x1 > x2) then; xl = x2; xr = x1; else; xl = x1; xr = x2; endif
!
    rnod => getPointerDoubleRadiiExternalRadialNodes(this%exnod)
!
!  Handle both cases:
!
    jl = locate(xl,size(rnod),rnod)+1                      ! index of rnod right or on xl
    if (jl > size(rnod)) then
       call add(errmsg,1,'all nodes left of integration interval',myname)
       jr = jl-1                                           ! make jr-jl < 0, signifying no nodes in interval
       nullify(rnod); return
    endif
    jr = locate(xr,size(rnod),rnod)                             ! index of rnod left of xr
    if (jr < size(rnod)) then
       if (dabs(1.d0-xr/rnod(jr+1)) < epsilon(1.d0)) jr = jr+1   ! node jr+1 is on xr, take it
    endif
    if (jr == 0) then
       call add(errmsg,1,'all nodes right of integration interval',myname)
       jr = 0; jl = 1                                       ! make jr-jl < 0
       nullify(rnod); return
    endif
!
!  assign nodes within integration interval
!
    nnod = jr-jl+1
    if (nnod > 0) then
       allocate(xn(nnod))
       if (x2 >= x1) then
          xn = rnod(jl:jr)
       else
          xn = rnod(jr:jl:-1)                    ! reverse order if x2 < x1
       endif
    endif
    nullify(rnod)
    end subroutine getNodesGemini
!
end module geminiIntegrationEnvironment
