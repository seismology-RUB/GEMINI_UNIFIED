! ======================================================================
!  Drives integration of minors using EPISODE propagator
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
! Drives integration of minors.
! There are upwards minors (from starting radius to nodes)
! and downwards minors (from surface to nodes). External nodes should be
! chosen sufficiently dense to achieve desired sampling of solution.
! External nodes in the halfspace are not allowed. 
! This is checked before we come here.
!--------------------------------------------------------------------
module minorIntegration
    use propagateOde
    use minorsOdeSystem
    use geminiIntegrationEnvironment
    use complexBulirschStep
    use errorMessage
    use initialValues
    implicit none
    interface dealloc; module procedure deallocMinorIntegration; end interface
    interface operator (.done.); module procedure doneMinorIntegration; end interface
    interface operator (.errlevel.); module procedure getErrlevelMinorIntegration; end interface
    integer, parameter :: secrecy_minor_integration = 3                 ! print info if secrecy <= this value
    type minor_integration
       double precision :: rstart                                       ! starting radius at bottom
       integer :: jnodsr                                                ! index of source node
       type (gemini_integration_environment) :: intenv                  ! gemini as integration environment
       type (minors_ode_system) :: sode                                 ! minors as ode system
       double precision :: eps                                          ! desired accuracy
       integer :: secrecy                                               ! secrecy level
       type (complex_bulirsch_step) :: cbs                              ! complex bs-step as step engine
       double complex, dimension(:,:), allocatable :: u_minors          ! upward minors at nodes (6,jnod)
       double complex, dimension(:,:), allocatable :: d_minors          ! downward minors at nodes (6,jnod)
       double complex :: det                                            ! determinants at S/R node
       double complex, dimension(:,:,:), allocatable :: u_mtil          ! Upwards M-tilde (4x4) matrix at nodes
       double complex, dimension(:,:,:), allocatable :: d_mtil          ! Downwards M-tilde (4x4) matrix at nodes
       logical :: done                                                  ! integration done?
       type (error_message) :: errmsg                                   ! object owned error message
    end type
!
contains
!--------------------------------------------------------------------------
!> \brief Create object
!!  rstart:     starting radius (in)
!!  jnodsr:     index of S/R node (in)
!!  gem_intenv: integration environment (in)
!!  eps:        desired accuracy (in)
!!  secrecy:    anti-debugging level (in)
!
    subroutine createMinorIntegration(this,rstart,jnodsr,gem_intenv,eps,secrecy)
    type (minor_integration) :: this
    integer :: jnodsr
    type (gemini_integration_environment) :: gem_intenv
    double precision :: rstart,eps
    integer :: secrecy
    integer :: nnod
    character(len=22) :: myname = 'createMinorIntegration'
!
    call new(this%errmsg,myname)
    nnod = .nnod.(gem_intenv%exnod)
    if (jnodsr < 1 .or. jnodsr > nnod) then
       call add(this%errmsg,2,'Invalid index for S/R node',myname)
       return
    endif
    this%rstart = rstart
    this%jnodsr = jnodsr
    this%secrecy = secrecy
    this%intenv = gem_intenv    ! copy sufficient here, this%intenv will be updated
    this%eps = eps
    this%done = .false.
!
!  set up ode system
!  this%sode must point to this%intenv because the latter  
!  (member nl) will be updated during integration (and not gem_intenv)
!
    call this%sode%createMinors(this%intenv,this%errmsg)
    if (.level.(this%errmsg) == 2) return
!
!  allocate space for minors and dets
!
    nnod = .nnod.(this%intenv%exnod)
    allocate(this%u_minors(6,nnod))
    allocate(this%d_minors(6,nnod))
    this%u_minors = dcmplx(0.d0,0.d0)
    this%d_minors = dcmplx(0.d0,0.d0)
    allocate(this%u_mtil(4,4,nnod))
    allocate(this%d_mtil(4,4,nnod))
    this%u_mtil = dcmplx(0.d0,0.d0)
    this%d_mtil = dcmplx(0.d0,0.d0)
end subroutine createMinorIntegration
!----------------------------------------------------------------------
!> \brief Deallocate object
!
    subroutine deallocMinorIntegration(this)
    type (minor_integration) :: this
    if (allocated(this%u_minors)) deallocate(this%u_minors)
    if (allocated(this%d_minors)) deallocate(this%d_minors)
    if (allocated(this%u_mtil)) deallocate(this%u_mtil)
    if (allocated(this%d_mtil)) deallocate(this%d_mtil)
    call this%sode%deallocMinors()
    if (this%secrecy <= secrecy_minor_integration) call print(this%errmsg)
    call dealloc(this%errmsg)
    end subroutine deallocMinorIntegration
!----------------------------------------------------------------------
!> \brief Do minors integration
!
    subroutine doMinorIntegration(this)
    type (minor_integration) :: this
    type (propagate_ode) :: prop
    double complex, dimension(5) :: ystart
    double complex, dimension(:,:), allocatable :: zyn
    integer, dimension(4,4) :: idx
    integer :: jlu,jru,jld,jrd,nl,top,i,j,jsr
    double precision :: rearth,htry,rnodsr
    character(len=18) :: myname = 'doMinorIntegration'
!
    call addTrace(this%errmsg,myname)
!
!  Check whether S/R node is below starting radius.
!
    rnodsr = getDoubleRadiusSelectedExternalRadialNodes(this%intenv%exnod,this%jnodsr)
    if (this%rstart >= rnodsr) then
       if (this%secrecy <= secrecy_minor_integration) then
          print *,'minorIntegration: Starting radius above S/R-node. DSV-solutions set to zero'
       endif
       return
    endif   
!
!  Establish propagator object for upwards, and downwards integrations
!  The already created minorsOdeSystem object is passed
!  to the polymorphic base ODE system of propagateOde.
!  Also, the geminiIntegrationEnvironment object is passed to the
!  polymorphic integration environment of propagateOde.
!  The up to now only declared Bulirsch-Store integration object is passed
!  to the polymorphic base integration engine of propagateOde. It will
!  be created when running doIntegrationStep.
!
    call createPropagateOde(prop,this%sode,this%intenv,this%cbs,this%eps,this%secrecy)
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call dealloc(prop)
       call add(this%errmsg,2,'Problems creating propagator object',myname)
       return
    endif
!
!  Upwards integration of minors from starting radius to S/R node
!
    if (this%secrecy <= secrecy_minor_integration) then
       print *,'minorIntegration: Upwards integration started'
       print *,'minorIntegration: rstart = ',this%rstart
    endif
    call getLayerIndexNodeEarthmodel(this%intenv%nem,this%rstart,nl,top)
    if (top == 1) nl = nl+1
    call this%intenv%setInterval(nl)
    call minorsBottomInitialValues(this%intenv,this%rstart,ystart,htry,this%errmsg)  ! ystart and htry on output
    if (.level.(this%errmsg) == 2) return
    call doPropagateOde(prop,this%rstart,rnodsr,htry,real(ystart),imag(ystart))
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call add(this%errmsg,2,'Problems doing upward propagation',myname)
       call dealloc(prop)
       return
    endif
!
!  read out values into upward minors
!  only set minor values at available nodes, others stay zero
!
    call getIndexLimitsNodesPropagateOde(prop,jlu,jru)
    if (jru-jlu < 0) then
       call add(this%errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
       call dealloc(prop)
       return
    endif
    allocate(zyn(5,jru-jlu+1))
    call getComplexSolutionAtNodesPropagateOde(prop,zyn)
    this%u_minors(1:5,jlu:jru) = zyn
    this%u_minors(6,jlu:jru) = -this%u_minors(1,jlu:jru)/this%intenv%dll1
    deallocate(zyn)
    if (this%secrecy <= secrecy_minor_integration) then
       print *,'Upward minors:'
       do j = jlu,jru
          print *,j,this%u_minors(1:5,j)
       enddo
       print *,'minorIntegration: Upwards minor integration done'
    endif
!
!  Downwards integration from surface to S/R node
!
    if (this%secrecy <= secrecy_minor_integration) then
       print *,'minorIntegration: Downwards integration started'
    endif
    nl = getNlayNodeEarthmodel(this%intenv%nem)       ! layer index of surface
    call this%intenv%setInterval(nl)
    call minorsSurfaceInitialValues(this%intenv,ystart,htry,this%errmsg)  ! ystart and htry on output
    if (.level.(this%errmsg) == 2) return
    rearth = getEarthRadiusNodeEarthmodel(this%intenv%nem)
    call doPropagateOde(prop,rearth,rnodsr,htry,real(ystart),imag(ystart))
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call add(this%errmsg,2,'Problems doing downwards propagation',myname)
       call dealloc(prop)
       return
    endif
!
!  read out values into downward minors (descending order of nodes here)
!  only set minor values at available nodes, others stay zero
!
    call getIndexLimitsNodesPropagateOde(prop,jld,jrd)
    if (jrd-jld < 0) then
       call add(this%errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
       call dealloc(prop)
       return
    endif
    allocate(zyn(5,jrd-jld+1))
    call getComplexSolutionAtNodesPropagateOde(prop,zyn)
    this%d_minors(1:5,jrd:jld:-1) = zyn(1:5,1:jrd-jld+1)   ! descending order here
    this%d_minors(6,jld:jrd) = -this%d_minors(1,jld:jrd)/this%intenv%dll1
    deallocate(zyn)
    if (this%secrecy <= secrecy_minor_integration) then
       print *,'Downward minors:'
       do j = jld,jrd
          print *,j,this%d_minors(1:5,j)
       enddo
       print *,'minorIntegration: Downwards minor integration done'
    endif
!
!  dealloc propagator
!
    call dealloc(prop)
!
!  form determinant at S/R node
!
    jsr = this%jnodsr
    this%det =  this%u_minors(1,jsr)*this%d_minors(6,jsr) + this%u_minors(6,jsr)*this%d_minors(1,jsr) &
              - this%u_minors(2,jsr)*this%d_minors(5,jsr) + this%u_minors(3,jsr)*this%d_minors(4,jsr) &
              + this%u_minors(4,jsr)*this%d_minors(3,jsr) - this%u_minors(5,jsr)*this%d_minors(2,jsr)
!
!  calculate mtil Matrix
!  Evaluate index-array for mtil/dmtil-matrices. Only 6 values
!  according to the number of minors are necessary.
!
    idx(1,2) = 1
    idx(1,3) = 2
    idx(1,4) = 3
    idx(2,3) = 4
    idx(2,4) = 5
    idx(3,4) = 6
!
!  upwards and downwards M_tilde matrix
!  M_tilde matrices have been zeroed when creating the object
!
    do i=1,4
       do j=i+1,4
          this%u_mtil(i,j,jlu:jru)  = this%u_minors(idx(i,j),jlu:jru)
          this%u_mtil(j,i,jlu:jru)  = -this%u_mtil(i,j,jlu:jru)
          this%d_mtil(i,j,jld:jrd)  = this%d_minors(idx(i,j),jld:jrd)
          this%d_mtil(j,i,jld:jrd)  = -this%d_mtil(i,j,jld:jrd)
       enddo
    enddo
    if (this%secrecy <= secrecy_minor_integration) then
       print *,'Mtilde up:'
       do j = jlu,jru
          print *,j,this%u_mtil(1,2,j),this%u_mtil(1,3,j),this%u_mtil(1,4,j),this%u_mtil(2,3,j),this%u_mtil(2,4,j)
       enddo
       print *,'Mtilde down:'
       do j = jld,jrd
          print *,j,this%d_mtil(1,2,j),this%d_mtil(1,3,j),this%d_mtil(1,4,j),this%d_mtil(2,3,j),this%d_mtil(2,4,j)
       enddo
    endif
!
    this%done = .true.
    if (this%secrecy <= secrecy_minor_integration) then
       print *,'minorIntegration: done'
    endif
    end subroutine doMinorIntegration
!---------------------------------------------------------------
!> \brief Minor computation done successfully ?
!
    function doneMinorIntegration(this) result(res)
    type (minor_integration), intent(in) :: this
    logical :: res
    res = this%done
    end function doneMinorIntegration
!-------------------------------------------------------------------------------
!> \brief print error message
!
    subroutine printErrmsgMinorIntegration(this)
    type (minor_integration) :: this
    call print(this%errmsg)
    end subroutine printErrmsgMinorIntegration
!-------------------------------------------------------------------------------
!> \brief get error level
!
    function getErrlevelMinorIntegration(this) result(res)
    type (minor_integration), intent(in) :: this
    integer :: res
    res = .level.(this%errmsg)
    end function getErrlevelMinorIntegration
end module minorIntegration
