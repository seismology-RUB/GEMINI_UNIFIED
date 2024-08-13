! ============================================================================
!  Drives integration of toroidal system using EPISODE propagator
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
!   Organizes and carries out integration of Green function
!   in the omega-el domain for the toroidal case. There is no minors
!   integration step required here. The toroidal system itself is integrated
!   upwards from the starting radius to the S/R node and downwards from the
!   surface to the S/R node. A determinant at the S/R node is calculated.
!   It and the solution stored at the nodes may later be used to calculate
!   the Green function for given jump vectors at the S/R node.
!--------------------------------------------------------------------
module toroidalIntegration
    use propagateOde
    use toroidalOdeSystem
    use geminiIntegrationEnvironment
    use complexBulirschStep
    use initialValues
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocToroidalIntegration; end interface
    interface operator (.errlevel.); module procedure getErrlevelToroidalIntegration; end interface
    integer, parameter :: secrecy_toroidal_integration = 3              ! print screen output if secrecy <= this value
    type toroidal_integration
       private
       double precision :: rstart                                       ! starting radius at bottom
       integer :: jnodsr                                                ! index of S/R node
       type (gemini_integration_environment) :: intenv                  ! gemini integration environment
       type (toroidal_ode_system) :: sode                               ! toroidal ode system
       double precision :: eps                                          ! desired accuracy
       integer :: secrecy                                               ! secrecy level
       type (complex_bulirsch_step) :: cbs                              ! complex bs-step as step engine
       double complex, dimension(:,:), allocatable :: u_tor             ! upward toros at nodes (2,jnod)
       double complex, dimension(:,:), allocatable :: d_tor             ! downward toros at nodes (2,jnod)
       double complex :: det                                            ! determinant at S/R node (jnodsr)
       double complex, dimension(:,:), pointer :: gf => null()          ! gf(2,nodes), excess from outside needed
       logical :: derivflag                                             ! derivatives desired
       logical :: done                                                  ! basis solutions done
       type (error_message) :: errmsg                                   ! error message owned by object
    end type
!
contains
!--------------------------------------------------------------------------
!> \brief Create object
!!  rstart:     starting radius (in)
!!  jnodsr:     index of S/R node (in)
!!  gem_intenv: gemini integration environment object (in)
!!  eps:        desired accuracy (in)
!!  derivflag:  compute derivatives of Green functions? (in)
!!  secrecy:    anti-debuggung level (in)
!
    subroutine createToroidalIntegration(this,rstart,jnodsr,gem_intenv,eps,derivflag,secrecy)
    type (toroidal_integration) :: this
    integer :: jnodsr
    type (gemini_integration_environment) :: gem_intenv
    double precision :: rstart,eps
    logical :: derivflag
    integer :: secrecy
    integer :: nnod
    character(len=22) :: myname = 'createToroidalIntegration'
!
    call new(this%errmsg,myname)
    nnod = .nnod.(gem_intenv%exnod)
    if (jnodsr < 1 .or. jnodsr > nnod) then
       call add(this%errmsg,2,'Invalid index for S/R node',myname)
       return
    endif
!
    this%rstart = rstart
    this%jnodsr = jnodsr
    this%secrecy = secrecy
    this%intenv = gem_intenv    ! copy sufficient here, this%intenv will be updated
    this%eps = eps
    this%derivflag = derivflag
    this%done = .false.
!
!  set up ode system
!  SODE's intenv member must point to this%intenv because it  
!  (member nl) will be updated during integration (and not gem_intenv)
!
    call this%sode%createToroidal(this%intenv,this%errmsg)
    if (.level.(this%errmsg) == 2) return
!
!  allocate space for minors and dets
!
    allocate(this%u_tor(2,nnod))
    allocate(this%d_tor(2,nnod))
    this%u_tor = dcmplx(0.d0,0.d0)
    this%d_tor = dcmplx(0.d0,0.d0)
    if (derivflag) then
       allocate(this%gf(3,nnod))
    else
       allocate(this%gf(2,nnod))
    endif
    this%gf = dcmplx(0.d0,0.d0)
    end subroutine createToroidalIntegration
!----------------------------------------------------------------------
!> \brief Deallocate object
!
    subroutine deallocToroidalIntegration(this)
    type (toroidal_integration) :: this
    if (allocated(this%u_tor)) deallocate(this%u_tor)
    if (allocated(this%d_tor)) deallocate(this%d_tor)
    if (associated(this%gf)) deallocate(this%gf)
    call this%sode%deallocToroidal()
    if (this%secrecy <= secrecy_toroidal_integration) call print(this%errmsg)
    call dealloc(this%errmsg)
    end subroutine deallocToroidalIntegration
!----------------------------------------------------------------------
!> \brief Do toroidal integration
!
    subroutine doToroidalIntegration(this)
    type (toroidal_integration) :: this
    type (propagate_ode) :: prop
    double complex, dimension(2) :: ystart
    double complex, dimension(:,:), allocatable :: zyn
    integer :: jlu,jru,jld,jrd,nl,top,j
    double precision :: rearth,rnodsr,htry
    character(len=21) :: myname = 'doToroidalIntegration'
!
    call addTrace(this%errmsg,myname)
!
!  Check whether S/R node is below starting radius.
!
    rnodsr = getDoubleRadiusSelectedExternalRadialNodes(this%intenv%exnod,this%jnodsr)
    if (this%rstart >= rnodsr) then
       if (this%secrecy <= secrecy_toroidal_integration) then
          print *,'toroidalIntegration: Starting radius above S/R-node. DSV-solutions set to zero'
       endif
       return
    endif   
!
!  Establish propagator object for upwards and downwards integrations
!
    call createPropagateOde(prop,this%sode,this%intenv,this%cbs,this%eps,this%secrecy)
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call add(this%errmsg,2,'Problems with creating toroidal propagator',myname)
       call dealloc(prop)
       return
    endif
!
!  Downwards integration from surface to S/R node
!  We have assured that it is above rstart
!
    if (this%secrecy <= secrecy_toroidal_integration) then
       print *,'toroidalIntegration: downwards integration started, rstart = ',this%rstart
    endif
    nl = getNlayNodeEarthmodel(this%intenv%nem)
    call this%intenv%setInterval(nl)
    call toridalSurfaceInitialValues(this%intenv,ystart,htry,this%errmsg)   ! ystart and htry on output
    if (.level.(this%errmsg) == 2) return
    rearth = getEarthRadiusNodeEarthmodel(this%intenv%nem)                   ! could be different from top node
    call doPropagateOde(prop,rearth,rnodsr,htry,real(ystart),imag(ystart))
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call add(this%errmsg,2,'Problems with doing toroidal propagation',myname)
       call dealloc(prop)
       return
    endif
!
!  read out values into d_tor (descending order of nodes here)
!  only set values at available nodes, others stay zero
!
    call getIndexLimitsNodesPropagateOde(prop,jld,jrd)
    if (jrd-jld < 0) then
       call add(this%errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
       call dealloc(prop)
       return
    endif
    allocate(zyn(2,jrd-jld+1))
    call getComplexSolutionAtNodesPropagateOde(prop,zyn)
    this%d_tor(1:2,jrd:jld:-1) = zyn(1:2,1:jrd-jld+1)         ! inverse order here
    if (this%secrecy <= secrecy_toroidal_integration) then
       print *,'doToroidalIntegration: downward minors:'
       do j = jld,jrd
          print *,j,this%d_tor(1:2,j)
       enddo
    endif
    deallocate(zyn)
!
!  Upwards integration of toroidal solution from starting radius to S/R node,
!  which we have assured to be above rstart
!
    if (this%secrecy <= secrecy_toroidal_integration) then
       print *,'toroidalIntegration: upwards integration started'
       print *,'rstart = ',this%rstart
    endif
    if (.level.(this%errmsg) == 2) return
    call getLayerIndexNodeEarthmodel(this%intenv%nem,this%rstart,nl,top)
    if (top == 1) nl = nl+1
    call this%intenv%setInterval(nl)
    call toroidalBottomInitialValues(this%intenv,this%rstart,ystart,htry,this%errmsg) ! ystart and htry on output
    if (.level.(this%errmsg) == 2) return
    call doPropagateOde(prop,this%rstart,rnodsr,htry,real(ystart),imag(ystart))
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call add(this%errmsg,2,'Problems with doing toroidal propagator',myname)
       call dealloc(prop)
       return
    endif
!
!  read out values into upward toroidal solution
!  only set values at available nodes, others stay zero
!
    call getIndexLimitsNodesPropagateOde(prop,jlu,jru)
    if (jru-jlu < 0) then
       call add(this%errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
       call dealloc(prop)
       return
    endif
    allocate(zyn(2,jru-jlu+1))
    call getComplexSolutionAtNodesPropagateOde(prop,zyn)
    this%u_tor(1:2,jlu:jru) = zyn(1:2,1:jru-jlu+1)
    call dealloc(prop)
    deallocate(zyn)
    if (this%secrecy <= secrecy_toroidal_integration) then
       print *,'toroidalIntegration: upwards integration from starting radius done'
    endif
!
!  dealloc propagator 
!
    call dealloc(prop)
!
!  form determinant at S/R nodes
!
    this%det = this%u_tor(1,this%jnodsr)*this%d_tor(2,this%jnodsr) &
              -this%d_tor(1,this%jnodsr)*this%u_tor(2,this%jnodsr)
!
    this%done = .true.
!
    if (this%secrecy <= secrecy_toroidal_integration) then
       print *,'toroidalIntegration: done'
    endif
    end subroutine doToroidalIntegration
!---------------------------------------------------------------------
! \brief Compute GF for given jump vector
!! The GF jumps at the S/R point. Since the S/R point is a node
!! we should assign two values to this node. We do not do that but
!! instead follow the convention that a potential receiver or source
!! at S/R always gets the upper value (i.e. is assumed to sit at the
!! upper one of the two virtual nodes at S/R). Consequently, we store
!! the GF below S/R only beginning with the node below S/R.
!
    subroutine computeGFToroidalIntegration(this,zsr)
    type (toroidal_integration) :: this
    double complex, dimension(:) :: zsr
    integer :: jsr,nnod,n
    character(len=28) :: myname = 'computeGFToroidalIntegration'
!
    call addTrace(this%errmsg,myname)
    if (size(zsr) /= 2) then
       call add(this%errmsg,2,'Invalid dimension of jump vector',myname)
       return
    endif
    jsr = this%jnodsr
    nnod = .nnod.(this%intenv%exnod)
!
!  GF at nodes below S/R point, exclude jnodsr
!
    forall (n = 1:jsr-1) this%gf(:,n) = &
         (-zsr(1)*this%d_tor(2,jsr)+zsr(2)*this%d_tor(1,jsr))*this%u_tor(:,n)/this%det
!
!  GF at nodes above S/R point, include jnodsr
!
    forall (n = jsr:nnod) this%gf(:,n) = &
         (-zsr(1)*this%u_tor(2,jsr)+zsr(2)*this%u_tor(1,jsr))*this%d_tor(:,n)/this%det
!
!  compute radial derivatives if requested
!
    if (this%derivflag) call derivativesGFToroidalIntegration(this)
    end subroutine computeGFToroidalIntegration
!---------------------------------------------------------------------
!  \brief Compute radial derivatives of Green function at nodes
!
    subroutine derivativesGFToroidalIntegration(this)
    type (toroidal_integration) :: this
    integer :: nnod,n,top,nl,ios
    double precision, dimension(:), pointer :: rnod
    double complex, dimension(2,2) :: a
    double complex, dimension(2) :: gh
!
    rnod => getPointerDoubleRadiiExternalRadialNodes(this%intenv%exnod)
    nnod = size(rnod)
    do n = 1,nnod
       call getLayerIndexNodeEarthmodel(this%intenv%nem,rnod(n),nl,top)
       call this%intenv%setInterval(nl)
       call this%sode%getSysmatComplex(rnod(n),a,ios)
       gh = matmul(a,this%gf(1:2,n))
       this%gf(3,n) = gh(1)             ! displacement derivative only, dW/dr = gh(1)
    enddo
    end subroutine derivativesGFToroidalIntegration
!-------------------------------------------------------------------
!> \brief Get Green function at nodes
!
    function getGFToroidalIntegration(this) result(res)
    type (toroidal_integration) :: this
    double complex, dimension(:,:), pointer :: res
    res => this%gf
    end function getGFToroidalIntegration
!-------------------------------------------------------------------------------
!> \brief print error message
!
    subroutine printErrmsgToroidalIntegration(this)
    type (toroidal_integration) :: this
    call print(this%errmsg)
    end subroutine printErrmsgToroidalIntegration
!-------------------------------------------------------------------------------
!> \brief get error level
!
    function getErrlevelToroidalIntegration(this) result(res)
    type (toroidal_integration), intent(in) :: this
    integer :: res
    res = .level.(this%errmsg)
    end function getErrlevelToroidalIntegration
!
end module toroidalIntegration
