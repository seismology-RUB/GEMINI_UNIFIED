! ======================================================================
!  Drives integration of spheroidal minors using EPISODE propagator
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
!   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
! Drives integration of spheroidal minors.
! There are upwards minors (stored between starting radius and a source node)
! and downwards minors (stored between surface and a source node). If a range of source nodes
! is specified upward integration is done from the starting radius to the uppermost node,
! and downward integration is done from surface to bottommost node or to the
! starting radius. Hence, both type of minors are available between lowermost and
! uppermost source node. If only one source node is specified, the upward integration
! will be done from the starting radius to the source node, and downward integration
! will be done from the surface to the source node or the starting radius.
! If a source node is below the starting radius, it will not
! be reached by either integration. This is physically reasonable as it cannot excite
! seismic motion. Minors and determinants there are explicitly set to zero then.
! Determinants computed from the minors are only calculated within the source node range.
!---------------------------------------------------------------------------------------
module spheroidalMinors
    use propagateOde
    use minorsOdeSystem
    use fluidOdeSystem
    use radialOdeSystem
    use geminiIntegrationEnvironment
    use complexBulirschStep
    use errorMessage
    use initialValues
    implicit none
    interface dealloc; module procedure deallocSpheroidalMinors; end interface
    interface operator (.done.); module procedure doneSpheroidalMinors; end interface
    integer, parameter :: secrecy_spheroidal_minors = 3                 ! print info if secrecy <= this value
    type spheroidal_minors
       type (gemini_integration_environment) :: intenv                  ! gemini as integration environment
       double complex, dimension(:,:), allocatable :: u_minors          ! upward minors at nodes (6,jnod)
       double complex, dimension(:,:), allocatable :: d_minors          ! downward minors at nodes (6,jnod)
       double complex, dimension(:), allocatable :: det                 ! determinants at nodes
       double complex, dimension(:,:,:), allocatable :: u_mtil          ! Upwards minors' matrix (4x4) at nodes
       double complex, dimension(:,:,:), allocatable :: d_mtil          ! Downwards minors' matrix (4x4) at nodes
       double complex, dimension(:,:,:), allocatable :: u_meps          ! Upwards epsilon*minors' matrix (4x4) at nodes
       double complex, dimension(:,:,:), allocatable :: d_meps          ! Downwards epsilon*minors' matrix (4x4) at nodes
       double precision, dimension(:), allocatable :: u_lgynorm         ! Upwards log-normalization factors at nodes 
       double precision, dimension(:), allocatable :: d_lgynorm         ! Downwards log-normalization factors at nodes 
       logical :: done                                                  ! integration done?
    end type spheroidal_minors
!
contains
!--------------------------------------------------------------------------
!> \brief Create object
!!  jsl:        first index of source node range (in)
!!  jsr:        last  index of source node range (in)
!!  gem_intenv: integration environment (in)
!!  eps:        desired accuracy (in)
!!  secrecy:    anti-debugging level (in)
!
  subroutine computeSpheroidalMinors(this,jsl,jsr,gem_intenv,eps,secrecy,errmsg)
     type (spheroidal_minors) :: this
     integer :: jsl,jsr
     type (gemini_integration_environment) :: gem_intenv
     double precision :: eps
     integer :: secrecy
     type (error_message) :: errmsg
     type (radial_ode_system) :: ros                                  ! radial ode system
     type (fluid_ode_system) :: fos                                   ! fluid as ode system
     type (minors_ode_system) :: mos                                  ! minors as ode system
     type (complex_bulirsch_step) :: cbs                              ! complex bs-step as step engine
     integer :: nnod,n
     character(len=23) :: myname = 'computeSpheroidalMinors'
     type (propagate_ode) :: prop
     double complex, dimension(:), allocatable :: ystart,yend         ! start and end values, dimension depends on nvar
     double complex, dimension(:,:), allocatable :: zyn
     integer, dimension(:), pointer :: state_exnod
     integer, dimension(4,4) :: idx
     integer :: jlu,jru,jld,jrd,nl,top,i,j,nvar
     double precision, dimension(:), allocatable :: lgy
     double precision :: rstart,rearth,htry,rsl,rmin,rmax,rbeg,rend
     double precision :: vf,lgnorm,lgnormend,lgn
     logical :: liquid_flag_cur,liquid_flag_old
   !
     call addTrace(errmsg,myname)
   !
     nnod = .nnod.(gem_intenv%exnod)
     if (jsl < 1 .or. jsl > nnod) then
        call add(errmsg,2,'Invalid index for first source node index',myname)
        return
     endif
     if (jsr < 1 .or. jsr > nnod) then
        call add(errmsg,2,'Invalid index for last source node index',myname)
        return
     endif
   !
     state_exnod => getPointerStateExternalRadialNodes(gem_intenv%exnod)
     this%intenv = gem_intenv    ! copy sufficient here, this%intenv will be updated
     this%done = .false.
   !
   !  allocate space for minors and mtils and dets
   !  explicitly zero here
   !
     allocate(this%u_minors(6,nnod))
     allocate(this%d_minors(6,nnod))
     this%u_minors = dcmplx(0.d0,0.d0)
     this%d_minors = dcmplx(0.d0,0.d0)
     allocate(this%u_mtil(4,4,nnod))
     allocate(this%d_mtil(4,4,nnod))
     this%u_mtil = dcmplx(0.d0,0.d0)
     this%d_mtil = dcmplx(0.d0,0.d0)
     allocate(this%u_meps(4,4,nnod))
     allocate(this%d_meps(4,4,nnod))
     this%u_meps = dcmplx(0.d0,0.d0)
     this%d_meps = dcmplx(0.d0,0.d0)
     allocate(this%u_lgynorm(nnod))
     allocate(this%d_lgynorm(nnod))
     this%u_lgynorm = 0.d0
     this%d_lgynorm = 0.d0
     allocate(this%det(nnod))
     this%det = dcmplx(0.d0,0.d0)
   !
   !  Relation between minor matrix and linearly indexed minors. Only 6 values
   !  according to the number of minors are necessary.
   !
     idx(1,2) = 1
     idx(1,3) = 2
     idx(1,4) = 3
     idx(2,3) = 4
     idx(2,4) = 5
     idx(3,4) = 6
   !
   !  Lower limit of downward integration from surface. It is
   !  max(rsl,rstart).
   !
     rstart = startRadiusInitialValues(gem_intenv,eps,'sph',errmsg)
     if (.level.errmsg == 2) return
     rsl = getDoubleRadiusSelectedExternalRadialNodes(this%intenv%exnod,jsl)
     rmin = max(rstart,rsl)
   !
   !  Upper limit of upward integration from starting radius.
   !  rmax = rsr
   !
     rmax = getDoubleRadiusSelectedExternalRadialNodes(this%intenv%exnod,jsr)
   !
   !  if rstart >= rmax, then minors at all source nodes vanish.
   !  Spheroidal motion for this value of wavenumber and frequency is not excited,
   !  Nothing to do. Solution stays zero.
   !
     if (rstart >= rmax) then
        if (secrecy <= secrecy_spheroidal_minors) then
           print *,trim(myname),': starting radius above all source nodes. Solution remains zero'
        endif
        this%done = .true.
        return
     endif
     rearth = .rearth.this%intenv%nem
   !--------------------------------------------------------------------------------------------
   !  Upwards (towards source) integration of minors from starting radius to rmax
   !
     if (secrecy <= secrecy_spheroidal_minors) then
        print *,'spheroidalMinors: Upwards integration started'
        print *,'spheroidalMinors: rstart = ',rstart,', rmax = ',rmax
     endif
   !
   !  Set the smallest possible value for the end of the current integration interval
   !  and start loop over possible inner core--outer core--mantle-ocean sequence
   !  
     rend = rstart
     rbeg = rstart
     liquid_flag_old = .false.                               ! here dummy, set at end of while loop
     lgnorm = 0.d0                                           ! initialize current log normalization factor
     do while(rend < rmax)
        call getLayerIndexNodeEarthmodel(this%intenv%nem,rbeg,nl,top)
        if (top == 1) nl = nl+1                              ! move to next layer when rbeg is at top of a layer
        liquid_flag_cur = (this%intenv%nem).isfluid.nl       ! is layer fluid?
        call this%intenv%setInterval(nl)
        rend = findAbovePhaseChangeNodeEarthmodel(this%intenv%nem,nl)     ! top radius of layer just below next phase change
        rend = min(rend,rmax)
        if (secrecy <= secrecy_spheroidal_minors) then
           print *,'Integrating from ',rbeg,' to ',rend,'; State: ',liquid_flag_cur,'; Layer: ',nl
        endif
      !
      !  set up ode system
      !  sode must point to this%intenv because the latter  
      !  (member nl) will be updated during integration (and not gem_intenv)
      !
      !  Establish propagator object for integration
      !  The odeSystem object is passed
      !  to the polymorphic base ODE system of propagateOde.
      !  Also, the geminiIntegrationEnvironment object is passed to the
      !  polymorphic integration environment of propagateOde.
      !  The up to now only declared Bulirsch-Stoer integration object is passed
      !  to the polymorphic base integration engine of propagateOde. It will
      !  be created when running doIntegrationStep.
      !
        if (liquid_flag_cur) then                              ! also valid for l=0
           call fos%createFluid(this%intenv,errmsg)
           if (.level.errmsg == 2) return
           nvar = fos%getNvar()
           allocate(ystart(nvar))
           call createPropagateOde(prop,fos,this%intenv,cbs,eps,secrecy,.true.,lgnorm)
           if (.errlevel.prop == 2) then
              call printErrmsgPropagateOde(prop)
              call dealloc(prop)
              call add(errmsg,2,'Problems creating propagator object',myname)
              return
           endif
        else
           if (gem_intenv%dll1 < epsilon(1.d0)) then                 ! l = 0
              call ros%createRadial(this%intenv,errmsg)
              if (.level.errmsg == 2) return
              nvar = ros%getNvar()
              allocate(ystart(nvar))
              call createPropagateOde(prop,ros,this%intenv,cbs,eps,secrecy,.true.,lgnorm)
          else
              call mos%createMinors(this%intenv,errmsg)
              if (.level.errmsg == 2) return
              nvar = mos%getNvar()
              allocate(ystart(nvar))
              call createPropagateOde(prop,mos,this%intenv,cbs,eps,secrecy,.true.,lgnorm)
           endif
           if (.errlevel.prop == 2) then
              call printErrmsgPropagateOde(prop)
              call dealloc(prop)
              call add(errmsg,2,'Problems creating propagator object',myname)
              return
           endif
        endif
     !
     !  if we are at rstart
     !  get initial values, minorsBottomInitialValues handles solid and fluid case
     !  else we are at a phase transition and need to apply boundary conditions
     !   
        if (rbeg == rstart) then
           call minorsBottomInitialValues(this%intenv,rbeg,ystart,htry,errmsg)  ! ystart and htry on output
           if (.level.errmsg == 2) return
        else
           if (liquid_flag_cur .and. (.not. liquid_flag_old)) then          ! solid to liquid transition
              if (gem_intenv%dll1 < epsilon(1.d0)) then                     ! l = 0
                 ystart(1) = yend(1)
                 ystart(2) = yend(2)
              else
                 ystart(1) = yend(3)
                 ystart(2) = yend(5)
              endif
           else if ((.not.liquid_flag_cur) .and. liquid_flag_old) then      ! liquid to solid transition
              if (gem_intenv%dll1 < epsilon(1.d0)) then                     ! l = 0
                 ystart(1) = yend(1)
                 ystart(2) = yend(2)
              else
                 ystart(4) = yend(2)
                 ystart(2) = yend(1)
                 ystart(1) = 0.d0
                 ystart(3) = 0.d0
                 ystart(5) = 0.d0
              endif
           else
              call add(errmsg,2,'Error with solid-fluid states! liquid_flag_cur and liquid_flag_old are equal!',myname)
              return
           endif
           deallocate(yend)               ! yend no longer needed
        endif
     !
     !  perform integration, using propagator created before
     !
        call doPropagateOde(prop,rbeg,rend,htry,real(ystart),imag(ystart))
        if (.errlevel.prop == 2) then
           call printErrmsgPropagateOde(prop)
           call add(errmsg,2,'Problems doing upward propagation',myname)
           call dealloc(prop)
           return
        endif
     !
     !  read out values into upward minors
     !  account for differing number of variables in minor and fluid integration, respectively
     !  only set minor values at available nodes, others stay zero
     !  store minors for liquid case into first two columns of this%u   
     !
        call getIndexLimitsNodesPropagateOde(prop,jlu,jru)
        if (jru-jlu .ge. 0) then
           allocate(zyn(nvar,jru-jlu+1))
           allocate(lgy(jru-jlu+1))
           call getComplexSolutionAtNodesPropagateOde(prop,zyn)
           call getLogNormalizationAtNodesPropagateOde(prop,lgy)
           this%u_lgynorm(jlu:jru) = lgy
           this%u_minors(1:nvar,jlu:jru) = zyn
           this%u_minors(6,jlu:jru) = -this%u_minors(1,jlu:jru)/this%intenv%dll1
           deallocate(zyn,lgy)
     !
     !  calculate minors' matrix, only needed for solid regions and l > 0
     !  also calculate epsilon times minors' matrix (eq. 3)
     !  minor matrices have been zeroed when creating the object
     !
           if ((.not. liquid_flag_cur) .and. gem_intenv%dll1 > epsilon(1.d0)) then
              do i=1,4
                 do j=i+1,4
                    vf = -1.d0+2.d0*mod(i+j,2)     ! -1 for even (i+j), +1 for odd (i+j)
                    this%u_mtil(i,j,jlu:jru)  = this%u_minors(idx(i,j),jlu:jru)
                    this%u_mtil(j,i,jlu:jru)  = -this%u_mtil(i,j,jlu:jru)
                    this%u_meps(i,j,jlu:jru)  = vf*this%u_minors(7-idx(i,j),jlu:jru)
                    this%u_meps(j,i,jlu:jru)  = -this%u_meps(i,j,jlu:jru)
                 enddo
              enddo
              if (secrecy <= secrecy_spheroidal_minors) then
                 print *,'Mtilde up:'
                 do j = jlu,jru
                    write(6,'(i6,10e15.3)') j,this%u_mtil(1,2,j),this%u_mtil(1,3,j),&
                         this%u_mtil(1,4,j),this%u_mtil(2,3,j),this%u_mtil(2,4,j)
                 enddo
              endif
           endif
        endif
        if (secrecy <= secrecy_spheroidal_minors) then
           print *,'Upward minors:'
           do j = jlu,jru
              write(6,'(i6,10e15.3)') j,this%u_minors(1:nvar,j)
           enddo
        endif
     !
     !  update begin of next interval as end of last one
     !  overwrite start value with value at end of interval   
     !
        liquid_flag_old = liquid_flag_cur
        rbeg = rend
        allocate(yend(nvar))                               ! new yend with right dimension
        call getComplexSolutionPropagateOde(prop,yend)
        call getLogNormalizationPropagateOde(prop,lgnormend)
        lgnorm = lgnormend                                    ! new start value for normalization
      !
      !  deallocate propagator and ODE
      !
        call dealloc(prop)
        if (liquid_flag_cur) then
           call fos%deallocFluid()
        else
           if (gem_intenv%dll1 > epsilon(1.d0)) then
              call mos%deallocMinors()
           else
              call ros%deallocRadial()
           endif
        endif
        deallocate(ystart)
     enddo                                      ! end of lop over liquid-solid sequence
     if (allocated(yend)) deallocate(yend)
   !
   !  normalize all upward solutions to maximum reached normalization factor
   !
     do n = 1,nnod
        lgn = this%u_lgynorm(n)-lgnormend
        this%u_minors(:,n) = this%u_minors(:,n)*10**lgn
        this%u_mtil(:,:,n) = this%u_mtil(:,:,n)*10**lgn
        this%u_meps(:,:,n) = this%u_meps(:,:,n)*10**lgn
     enddo
     if (secrecy <= secrecy_spheroidal_minors) then
        print *,'spheroidalMinors: Upwards minor integration done'
     endif
   !----------------------------------------------------------------------------------------
   !  Downwards (towards source) integration from surface to rmin
   !
     if (secrecy <= secrecy_spheroidal_minors) then
        print *,'spheroidalMinors: Downwards integration started'
     endif
   !
   !  Set the largest possible value for the end of the current integration interval
   !  and start loop over possible ocean-crust-mantle-outer core-inner core sequence
   !
     rbeg = rearth
     rend = rearth
     liquid_flag_old = .false.
     lgnorm = 0.d0
     do while(rend > rmin)
        call getLayerIndexNodeEarthmodel(this%intenv%nem,rbeg,nl,top)
        liquid_flag_cur = (this%intenv%nem).isfluid.nl
        call this%intenv%setInterval(nl)
        rend = findBelowPhaseChangeNodeEarthmodel(this%intenv%nem,nl)  ! bottom radius of layer just above next phase change
        rend = max(rend,rmin)
        if (secrecy <= secrecy_spheroidal_minors) then
           print *,'Integrating from ',rbeg,' to ',rend,'; State: ',liquid_flag_cur,'; Layer: ',nl
        endif
     !
     !  set up ode system
     !  sode must point to this%intenv because the latter  
     !  (member nl) will be updated during integration (and not gem_intenv)
     !
        if (liquid_flag_cur) then                        ! also covers l = 0 case
           call fos%createFluid(this%intenv,errmsg)
           if (.level.errmsg == 2) return
           nvar = fos%getNvar()
           allocate(ystart(nvar))
           call createPropagateOde(prop,fos,this%intenv,cbs,eps,secrecy,.true.,lgnorm)
           if (.errlevel.prop == 2) then
              call printErrmsgPropagateOde(prop)
              call dealloc(prop)
              call add(errmsg,2,'Problems creating propagator object',myname)
              return
           endif
        else
           if (gem_intenv%dll1 < epsilon(1.d0)) then                 ! l = 0
              call ros%createRadial(this%intenv,errmsg)
              if (.level.errmsg == 2) return
              nvar = ros%getNvar()
              allocate(ystart(nvar))
              call createPropagateOde(prop,ros,this%intenv,cbs,eps,secrecy,.true.,lgnorm)
          else
              call mos%createMinors(this%intenv,errmsg)
              if (.level.errmsg == 2) return
              nvar = mos%getNvar()
              allocate(ystart(nvar))
              call createPropagateOde(prop,mos,this%intenv,cbs,eps,secrecy,.true.,lgnorm)
           endif
           if (.errlevel.prop == 2) then
              call printErrmsgPropagateOde(prop)
              call dealloc(prop)
              call add(errmsg,2,'Problems creating propagator object',myname)
              return
           endif
        endif
     !
     ! we are at the surface, minorsSurfaceInitialValues handles both solid and liquid case 
     !
        if (rbeg == rearth) then        
           call minorsSurfaceInitialValues(this%intenv,ystart,htry,errmsg)  ! ystart and htry on output
           if (.level.errmsg == 2) return
     !
     !  we are at some transition between solid and liquid
     !
        else
           if (liquid_flag_cur .and. (.not. liquid_flag_old)) then          ! solid to liquid transition
              if (gem_intenv%dll1 < epsilon(1.d0)) then                     ! l = 0
                 ystart(1) = yend(1)
                 ystart(2) = yend(2)
              else
                 ystart(1) = yend(3)
                 ystart(2) = yend(5)
              endif
           else if ((.not. liquid_flag_cur) .and. liquid_flag_old) then     ! liquid to solid transition
              if (gem_intenv%dll1 < epsilon(1.d0)) then                     ! l = 0
                 ystart(1) = yend(1)
                 ystart(2) = yend(2)
              else
                 ystart(4) = yend(2)
                 ystart(2) = yend(1)
                 ystart(1) = 0.d0
                 ystart(3) = 0.d0
                 ystart(5) = 0.d0
              endif
           else
              call add(errmsg,2,'Error with solid-fluid states! liquid_flag_cur and liquid_flag_old are equal!',myname)
              return
           endif
           deallocate(yend)               ! yend no longer needed
        endif
     !
     !  perform integration, using propagator created before
     !
        call doPropagateOde(prop,rbeg,rend,htry,real(ystart),imag(ystart))
        if (.errlevel.prop == 2) then
           call printErrmsgPropagateOde(prop)
           call add(errmsg,2,'Problems doing downwards propagation',myname)
           call dealloc(prop)
           return
        endif
     !
     !  read out values into downward minors
     !  solution is in descending order of nodes here
     !  but we store it into d_minors in ascending order of nodes
     !
        call getIndexLimitsNodesPropagateOde(prop,jld,jrd)
        if (jrd-jld < 0) then
           call add(errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
           call dealloc(prop)
           return
        endif
        allocate(zyn(nvar,jrd-jld+1))
        allocate(lgy(jrd-jld+1))
        call getComplexSolutionAtNodesPropagateOde(prop,zyn)
        call getLogNormalizationAtNodesPropagateOde(prop,lgy)
        this%d_lgynorm(jrd:jld:-1) = lgy(1:jrd-jld+1)
        this%d_minors(1:nvar,jrd:jld:-1) = zyn(1:nvar,1:jrd-jld+1)             ! descending order here
        this%d_minors(6,jld:jrd) = -this%d_minors(1,jld:jrd)/this%intenv%dll1
        deallocate(zyn,lgy)
     !
     !  calculate minors' matrix, only needed for solid regions and l > 0
     !  also calculate epsilon times minors' matrix (eq. 3)
     !  minor matrices have been zeroed when creating the object
     !
        if ((.not. liquid_flag_cur) .and. gem_intenv%dll1 > epsilon(1.d0)) then
           do i=1,4
              do j=i+1,4
                 vf = -1.d0+2.d0*mod(i+j,2)     ! -1 for even (i+j), +1 for odd (i+j)
                 this%d_mtil(i,j,jld:jrd)  = this%d_minors(idx(i,j),jld:jrd)
                 this%d_mtil(j,i,jld:jrd)  = -this%d_mtil(i,j,jld:jrd)
                 this%d_meps(i,j,jld:jrd)  = vf*this%d_minors(7-idx(i,j),jld:jrd)
                 this%d_meps(j,i,jld:jrd)  = -this%d_meps(i,j,jld:jrd)
              enddo
           enddo
           if (secrecy <= secrecy_spheroidal_minors) then
              print *,'Mtilde down:'
              do j = jld,jrd
                 write(6,'(i6,10e15.3)') j,this%d_mtil(1,2,j),this%d_mtil(1,3,j),this%d_mtil(1,4,j),&
                      this%d_mtil(2,3,j),this%d_mtil(2,4,j)
              enddo
           endif
        endif
        if (secrecy <= secrecy_spheroidal_minors) then
           print *,'Downward minors:'
           do j = jld,jrd
              write(6,'(i6,10e15.3)') j,this%d_minors(1:nvar,j)
           enddo
        endif
     !
     !  update begin of next interval as end of last one
     !  overwrite start value with value at end of interval   
     !
        liquid_flag_old = liquid_flag_cur
        rbeg = rend
        allocate(yend(nvar))                               ! new yend with right dimension
        call getComplexSolutionPropagateOde(prop,yend)
        call getLogNormalizationPropagateOde(prop,lgnormend)
        lgnorm = lgnormend
      !
      !  deallocate propagator and ODE
      !
        call dealloc(prop)
        if (liquid_flag_cur) then
           call fos%deallocFluid()
        else
           if (gem_intenv%dll1 > epsilon(1.d0)) then
              call mos%deallocMinors()
           else
              call ros%deallocRadial()
           endif
        endif
        deallocate(ystart)
     enddo                               ! end of lop over liquid-solid sequence
     if (allocated(yend)) deallocate(yend)
   !
   !  normalize all downward solutions to maximum reached normalization factor
   !
     do n = 1,nnod
        lgn = this%d_lgynorm(n)-lgnormend
        this%d_minors(:,n) = this%d_minors(:,n)*10**lgn
        this%d_mtil(:,:,n) = this%d_mtil(:,:,n)*10**lgn
        this%d_meps(:,:,n) = this%d_meps(:,:,n)*10**lgn
     enddo
     if (secrecy <= secrecy_spheroidal_minors) then
        print *,'spheroidalMinors: Downwards minor integration done'
     endif
   !
   !-------------------------------------------------------------------------------------------------
   !  form determinant at source nodes, if a source node was not reached
   !  everything stays zero
   !
     do j = jsl,jsr
        if (state_exnod(j) == 1 .or. gem_intenv%dll1 < epsilon(1.d0)) then    ! liquid case or l=0
           this%det(j) = this%u_minors(1,j)*this%d_minors(2,j)-this%d_minors(1,j)*this%u_minors(2,j)
        else
           this%det(j) =  &
             this%u_minors(1,j)*this%d_minors(6,j) + this%u_minors(6,j)*this%d_minors(1,j) &
             - this%u_minors(2,j)*this%d_minors(5,j) + this%u_minors(3,j)*this%d_minors(4,j) &
             + this%u_minors(4,j)*this%d_minors(3,j) - this%u_minors(5,j)*this%d_minors(2,j)
        endif
     enddo
     if (secrecy <= secrecy_spheroidal_minors) then
        print *,'Spheroidal determinants:'
        do j = jsl,jsr
           write(6,'(i6,2e15.3)') j,this%det(j)
        enddo
     endif
   !
     this%done = .true.
     if (secrecy <= secrecy_spheroidal_minors) call print(errmsg)
  end subroutine computeSpheroidalMinors
!----------------------------------------------------------------------
!> \brief Deallocate object
!
  subroutine deallocSpheroidalMinors(this)
     type (spheroidal_minors) :: this
     if (allocated(this%u_minors)) deallocate(this%u_minors)
     if (allocated(this%d_minors)) deallocate(this%d_minors)
     if (allocated(this%u_mtil)) deallocate(this%u_mtil)
     if (allocated(this%d_mtil)) deallocate(this%d_mtil)
     if (allocated(this%u_meps)) deallocate(this%u_meps)
     if (allocated(this%d_meps)) deallocate(this%d_meps)
     if (allocated(this%u_lgynorm)) deallocate(this%u_lgynorm)
     if (allocated(this%d_lgynorm)) deallocate(this%d_lgynorm)
     if (allocated(this%det)) deallocate(this%det)
  end subroutine deallocSpheroidalMinors
!---------------------------------------------------------------
!> \brief Minor computation done successfully ?
!
  function doneSpheroidalMinors(this) result(res)
     type (spheroidal_minors), intent(in) :: this
     logical :: res
     res = this%done
  end function doneSpheroidalMinors
end module spheroidalMinors
