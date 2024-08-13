! ============================================================================
!  Compute minors of toroidal system using EPISODE propagator
! ============================================================================
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
!---------------------------------------------------------------------
!  Drives integration of toroidal minors.
!  There are upwards minors (stored between starting radius and a source node)
!  and downwards minors (stored between surface and a source node). If a range of source nodes
!  is specified upward integration is done from starting radius to uppermost node, and
!  downward integration is done from surface to bottommost node or to the
!  starting radius. Hence, both type of minors are available between lowermost and
!  uppermost source node.
!  If only one node is a source node, the upward integration
!  is done from the starting radius to the source node, and downward integration
!  is done from the surface to the source node or the starting radius.
!  If a source node is below the starting radius, it will not
!  be reached by either integration. This is physically reasonable as it cannot excite
!  seismic motion. Minors and determinants there are explicitly set to zero.

!  For toroidal motion there are no separate minor 
!  equations but instead the toroidal system itself is integrated
!  Determinants are stored at the source nodes only where both type
!  of minors are available.
!--------------------------------------------------------------------
module toroidalMinors
    use propagateOde
    use toroidalOdeSystem
    use geminiIntegrationEnvironment
    use complexBulirschStep
    use initialValues
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocToroidalMinors; end interface
    interface operator (.done.); module procedure doneToroidalMinors; end interface
    integer, parameter :: secrecy_toroidal_minors = 3                   ! print screen output if secrecy <= this value
    type toroidal_minors
       type (gemini_integration_environment) :: intenv                  ! gemini as integration environment
       double complex, dimension(:,:), allocatable :: u_tor             ! upward toros at nodes (2,jnod)
       double complex, dimension(:,:), allocatable :: d_tor             ! downward toros at nodes (2,jnod)
       double complex, dimension(:), allocatable :: det                 ! determinant at source nodes
       logical :: done                                                  ! basis solutions done
    end type
!
contains
!--------------------------------------------------------------------------
!> \brief Compute toroidal minors
!!  Performs integration from starting radius upwards to uppermost of S/R nodes
!!  Performs integration from surface downwards to lowermost of S/R-nodes
!!  jsl:        index of bottommost source node (in)
!!  jsr:        index of topmost source node (in)
!!  gem_intenv: gemini integration environment object (in)
!!  eps:        desired accuracy (in)
!!  secrecy:    anti-debuggung level (in)
!
  subroutine computeToroidalMinors(this,jsl,jsr,gem_intenv,eps,secrecy,errmsg)
     type (toroidal_minors) :: this
     integer :: jsl,jsr
     type (gemini_integration_environment) :: gem_intenv
     double precision :: eps
     integer :: secrecy
     type (error_message) :: errmsg
     type (toroidal_ode_system) :: sode                               ! toriodal ode system
     type (complex_bulirsch_step) :: cbs                              ! complex bs-step as step engine
     type (propagate_ode) :: prop
     double complex, dimension(2) :: ystart
     double complex, dimension(:,:), allocatable :: zyn
     integer :: jlu,jru,jld,jrd,nl,top,j,js,nnod
     double precision :: rstart,rsl,rmin,rmax,htry,rsurf
     character(len=21) :: myname = 'computeToroidalMinors'
   !
     call addTrace(errmsg,myname)
   !
     nnod = .nnod.(gem_intenv%exnod)
     if (jsl < 1 .or. jsl > nnod) then
        call add(errmsg,2,'Invalid index for first source node',myname)
        return
     endif
     if (jsr < 1 .or. jsr > nnod) then
        call add(errmsg,2,'Invalid index for last source node',myname)
        return
     endif
   !
     rstart = startRadiusInitialValues(gem_intenv,eps,'tor',errmsg)
     if (.level.errmsg == 2) return
     this%intenv = gem_intenv    ! copy sufficient here, this%intenv will be updated
     this%done = .false.
   !
   !  set up ode system
   !  SODE's intenv member must point to this%intenv because it  
   !  (member nl) will be updated during integration (and not gem_intenv)
   !
     call sode%createToroidal(this%intenv,errmsg)
     if (.level.errmsg == 2) return
   !
   !  allocate space for minors
   !
     allocate(this%u_tor(2,nnod))
     allocate(this%d_tor(2,nnod))
     this%u_tor = dcmplx(0.d0,0.d0)
     this%d_tor = dcmplx(0.d0,0.d0)
     allocate(this%det(nnod))
     this%det = dcmplx(0.d0,0.d0)
   !
   !  check for vanishing wavenumber, minors stay zero
   !  no integration needed
   !    
     if (gem_intenv%dll1 < epsilon(1.d0)) then
        this%done = .true.
        return
     endif
   !
   !  Lower limit of downward integration from surface. It is
   !  max(rsl,rstart).
   !
     rsl = getDoubleRadiusSelectedExternalRadialNodes(this%intenv%exnod,jsl)
     rmin = max(rstart,rsl)
   !
   !  Upper limit of upward integration from starting radius.
   !  Later modified if ocean is present
   !
     rmax = getDoubleRadiusSelectedExternalRadialNodes(this%intenv%exnod,jsr)
   !
   !  if rstart >= rmax, then minors at all source nodes vanish.
   !  Toroidal motion for this value of wavenumber and frequency is not excited,
   !  neither by any source node. Nothing to do. Solution is zero.
   !
     if (rstart >= rmax) then
        if (secrecy <= secrecy_toroidal_minors) then
           print *,trim(myname),': starting radius above source nodes. Solution remains zero'
        endif
        this%done = .true.
        return
     endif
   !
   !  Establish propagator object for upwards and downwards integrations
   !  Use complex Bulirsch stepping (cbs)
   !
     call createPropagateOde(prop,sode,this%intenv,cbs,eps,secrecy,.false.,0.d0)
     if (.errlevel.prop == 2) then
        call printErrmsgPropagateOde(prop)
        call add(errmsg,2,'Problems with creating toroidal propagator',myname)
        call dealloc(prop)
        return
     endif
   !
   !  Downwards integration from surface to rmin
   !
     call toroidalSurfaceInitialValues(this%intenv,ystart,htry,rsurf,nl,errmsg)     ! ystart,htry,rsurf,nl on output, tests for ocean
     if (.level.errmsg == 2) return
     if (secrecy <= secrecy_toroidal_minors) then
        print *,'toroidalMinors: downwards integration started'
        print *,'rmin = ',rmin,' rsurf = ',rsurf
     endif
     if (rsurf < rmin) then                                                         ! all sources in ocean, toroidal minors all zero
        if (secrecy <= secrecy_toroidal_minors) then
           print *,'computeToroidalMinors: toroidal minors stay all zero because all sources are in ocean'
        endif
        this%done = .true.
        call dealloc(prop)
        call sode%deallocToroidal()
        if (secrecy <= secrecy_toroidal_minors) call print(errmsg)
        return
     endif
   !
   !  there are sources in the solid mantle
   !
     call this%intenv%setInterval(nl)
     call doPropagateOde(prop,rsurf,rmin,htry,real(ystart),imag(ystart))
     if (.errlevel.prop == 2) then
        call printErrmsgPropagateOde(prop)
        call add(errmsg,2,'Problems with doing toroidal propagation',myname)
        call dealloc(prop)
        return
     endif
   !
   !  read out values into d_tor
   !  solution is in descending order of nodes here but we store
   !  solution into d_tor in ascending order of nodes
   !  only set values at available nodes, others stay zero
   !
     call getIndexLimitsNodesPropagateOde(prop,jld,jrd)
     if (jrd-jld < 0) then
        call add(errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
        call dealloc(prop)
        return
     endif
     allocate(zyn(2,jrd-jld+1))
     call getComplexSolutionAtNodesPropagateOde(prop,zyn)
     this%d_tor(1:2,jrd:jld:-1) = zyn(1:2,1:jrd-jld+1)         ! inverse order here
     if (secrecy <= secrecy_toroidal_minors) then
        print *,'computeToroidalMinors: downward minors:'
        do j = jld,jrd
           write(6,'(i6,4e15.3)') j,this%d_tor(1:2,j)
        enddo
     endif
     deallocate(zyn)
   !
   !  Upwards integration of toroidal minors from starting radius to rmax,
   !  at least one of which we have assured to be above rstart
   !
     if (secrecy <= secrecy_toroidal_minors) then
        print *,'toroidalMinors: upwards integration started'
        print *,'rstart = ',rstart,' rmax = ',rmax
     endif
     rmax = min(rmax,rsurf)                                                  ! limit range for upwards integration to rsurf
     call getLayerIndexNodeEarthmodel(this%intenv%nem,rstart,nl,top)
     if (top == 1) nl = nl+1
     call this%intenv%setInterval(nl)
     call toroidalBottomInitialValues(this%intenv,rstart,ystart,htry,errmsg) ! ystart and htry on output
     if (.level.errmsg == 2) return
     call doPropagateOde(prop,rstart,rmax,htry,real(ystart),imag(ystart))
     if (.errlevel.prop == 2) then
        call printErrmsgPropagateOde(prop)
        call add(errmsg,2,'Problems with doing toroidal propagator',myname)
        call dealloc(prop)
        return
     endif
   !
   !  read out values into upward toroidal solution
   !  only set values at available nodes, others stay zero
   !
     call getIndexLimitsNodesPropagateOde(prop,jlu,jru)
     if (jru-jlu < 0) then
        call add(errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
        call dealloc(prop)
        return
     endif
     allocate(zyn(2,jru-jlu+1))
     call getComplexSolutionAtNodesPropagateOde(prop,zyn)
     this%u_tor(1:2,jlu:jru) = zyn(1:2,1:jru-jlu+1)
     call dealloc(prop)
     deallocate(zyn)
     if (secrecy <= secrecy_toroidal_minors) then
        print *,'toroidalMinors: upwards integration from starting radius done'
        do j = jlu,jru
           write(6,'(i6,4e15.3)') j,this%u_tor(1:2,j)
        enddo
     endif
   !
   !  form determinants at S/R nodes
   !
     forall (js = jsl:jsr) this%det(js) = this%u_tor(1,js)*this%d_tor(2,js)-this%d_tor(1,js)*this%u_tor(2,js)
   !
     this%done = .true.
   !
   !  dealloc propagator and sode
   !
     call dealloc(prop)
     call sode%deallocToroidal()
   !
     if (secrecy <= secrecy_toroidal_minors) call print(errmsg)
  end subroutine computeToroidalMinors
!----------------------------------------------------------------------
!> \brief Deallocate object
!
  subroutine deallocToroidalMinors(this)
     type (toroidal_minors) :: this
     if (allocated(this%u_tor)) deallocate(this%u_tor)
     if (allocated(this%d_tor)) deallocate(this%d_tor)
     if (allocated(this%det)) deallocate(this%det)
  end subroutine deallocToroidalMinors
!---------------------------------------------------------------
!> \brief Minor computation done successfully ?
!
  function doneToroidalMinors(this) result(res)
     type (toroidal_minors), intent(in) :: this
     logical :: res
     res = this%done
  end function doneToroidalMinors
!
end module toroidalMinors
