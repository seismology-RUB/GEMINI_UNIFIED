! ==============================================================================
!  Compute travel time table
! ==============================================================================
!----------------------------------------------------------------------------
!   Copyright 2020 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!  Main program to compute a travel time table for varying turning point radius
!  and receiver radius. Produces a table with values for travel time and epicentral
!  distance for each combination of turning point radius and receiver radius.
!  For fixed turning point radius, Delta(re) can be considered as a ray with
!  associated travel times T(re). For fixed receiver radius, T(p) versus Delta(p)
!  is a travel time curve. This is put into variable raytable((delta,T),re,rt).
!  In addition, receiver radii, turning point radii and slownesses are written to
!  a HDF file.
!  When the travel timeand distance from source to receiver is needed, one can add the times
!  and distances from source to turning point and from turning point to receiver.
!-----------------------------------------------------------------------------
program computeTravelTimeTable
    use hdf5
    use argumentParser
    use nodeEarthmodel
    use externalRadialNodes
    use rayIntegrationEnvironment
    use rayOdeSystem
    use propagateOde
    use realBulirschStep
    use splineEarthmodelCoefficients
    use inputParameter
    use errorMessage
    use mpiSupport
    use hdfWrapper
    use anyRankRealArray
    use anyRankIntegerArray
    implicit none
    type (argument_parser) :: ap
    type (input_parameter) :: inpar
    type (node_earthmodel) :: nem
    type (external_radial_nodes) :: exnod
    type (ray_integration_environment) :: ray_intenv
    type (ray_ode_system) :: sode                               ! ray ode system
    type (real_bulirsch_step) :: rbs                            ! real bs-step as step engine
    type (propagate_ode) :: prop                                ! propagator
    type (error_message) :: errmsg
    type (mpi_support) :: mpisup
    type (any_rank_real_array) :: arra
    integer(hid_t) :: fid,xferprp
    integer :: nnod,ntp,jlu,jru,nl,top
    integer :: j,i,n
    integer :: myrank,numtasks,ierr,secrecy
    integer, dimension(:), pointer :: state,layer
    real, dimension(:), allocatable, target :: d
    real, dimension(:,:,:), allocatable :: ray_table
    double precision, dimension(:), pointer :: rnod
    double precision, dimension(:), allocatable :: rtp,slowness
    double precision, dimension(:,:), allocatable :: yn
    double precision, dimension(2) :: ystart
    double precision :: eps,rbot,drtp,v,rtmax
    character (len=3) :: crank
    character(len=5) :: raytype
    character(len=max_length_string) :: parfile,errstr,attmode,raytablefile,rt_type
    character (len=80), dimension(8) :: para_keywords
    character (len=80), dimension(4) :: attn_keywords
    character(len=34) :: myname = 'computeTravelTimeTable'
    integer, parameter :: secrecy_ray = 5                         ! print screen output if secrecy <= this value
    data para_keywords/'ATTENUATION_MODE', 'TURNING_NODES_DR', 'TURNING_NODES_RMAX', 'RAY_TABLE_TYPE',&
                   &   'EARTH_MODEL', 'ACCURACY', 'FILE_RAY_TABLE', 'EXTERNAL_NODES_FROM_FILE'/
    data attn_keywords/'ELASTIC','ATTENUATION_ONLY','DISPERSION_ONLY','ATTENUATION_AND_DISPERSION'/
!-----------------------------------------------------------------------------
!  initialise MPI
!
    call new(mpisup)
    myrank = .myrank.mpisup
    write(crank,'(i3.3)') myrank
    numtasks = .numtasks.mpisup
!----------------------------------------------------------------------------
    call init(ap,myname,'Travel time tables for spherically symmetric earth model')
    call addPosarg(ap,'parfile','sval','ray parameter file')
    call addOption(ap,'-S',.false.,'Compute S travel time table, Default is P')
    call addOption(ap,'-s',.true.,'Secrecy level','ival','6')
    call parse(ap)
    parfile = ap.sval.'parfile'
    secrecy = ap.ival.'-s'
    if (ap.optset.'-S') then
       raytype = 'S'
    else
       raytype = 'P'
    end if
    if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); stop; endif
    if (myrank == 0) call document(ap)
    call dealloc(ap)
!-----------------------------------------------------------------------------
    call new(errmsg,myname)
!-----------------------------------------------------------------------------
    call createKeywordsInputParameter(inpar,para_keywords)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 10
    if (myrank == 0) call printInputParameter(inpar)
!-----------------------------------------------------------------------------
!  define variables derived from input data
!
    eps = inpar.dval.'ACCURACY'
    raytablefile = inpar.sval.'FILE_RAY_TABLE'
!-----------------------------------------------------------------------------
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
    call setXferprpCollectiveHDFWrapper(xferprp,errmsg)        ! set data transfer prop to MPIO collective
    if (.level.errmsg == 2) goto 10
    call createFileHDFWrapper(trim(raytablefile),fid,errmsg)
    if (.level.errmsg == 2) goto 10
!-----------------------------------------------------------------------------
!  check validity of attenuation mode
!
    attmode = inpar.sval.'ATTENUATION_MODE'
    if (.not.(attmode.equal.attn_keywords(1)) .and. .not.(attmode.equal.attn_keywords(2)) .and. &
        .not.(attmode.equal.attn_keywords(3)) .and. .not.(attmode.equal.attn_keywords(4))) then
       call add(errmsg,2,'Invalid attenuation mode',myname)
       goto 10
    endif
!-----------------------------------------------------------------------------
    call createNodeEarthmodel(nem,1,inpar.sval.'EARTH_MODEL',errmsg)
    if (.level.errmsg == 2) goto 10
    if (myrank == 0) call printNodeEarthmodel(nem)
!       
    if ((inpar.ival.'EXTERNAL_NODES_FROM_FILE') == 0) then
       call createExternalRadialNodes(exnod,1,parfile,nem,errmsg)
    else
       call createFromRadiiExternalRadialNodes(exnod,1,parfile,nem,errmsg)
    endif
    if (.level.errmsg == 2) goto 10
!------------------------------------------------------------------------------
!  make sure that lowest external node is above halfspace or inner homogeneous sphere
!
    call getDoubleRadiusBottomNodeExternalRadialNodes(exnod,rbot)
    if (rbot < .rhs.nem) then
       errstr = 'External nodes within inner homogeneous sphere (halfspace) are not allowed.'
       call add(errmsg,2,errstr,myname)
       errstr = 'Choose smaller radius for inner homogeneous sphere in your earth model!'
       call add(errmsg,2,errstr,myname)
       goto 10
    endif
!------------------------------------------------------------------------------
!  external node information
!
    nnod = .nnod.exnod
    rnod => getPointerDoubleRadiiExternalRadialNodes(exnod)
    state => getPointerStateExternalRadialNodes(exnod)
    layer => getPointerLayerExternalRadialNodes(exnod)
!----------------------------------------------------------------------------
!  define turning point nodes (from CMB to surface)
!
    drtp = inpar.dval.'TURNING_NODES_DR'
    rtmax = inpar.dval.'TURNING_NODES_RMAX'
    ntp = int(rtmax-.roc.nem)/drtp+1
    rt_type = inpar.sval.'RAY_TABLE_TYPE'
    allocate(rtp(ntp),slowness(ntp))
    do i = 1,ntp
       rtp(i) = .roc.nem+(i-1)*drtp
    enddo
!
!  HDF output of rnod and rtp
!
    call writeStringAttributeHDFWrapper(fid,'rayTableType',trim(rt_type),errmsg)
    if (.level.errmsg == 2) goto 10
    allocate(d(nnod)); d = real(rnod)
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(fid,'receiverRadii',arra,errmsg)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
    
    allocate(d(ntp)); d = real(rtp)
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(fid,'turningPointRadii',arra,errmsg)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
!
!  spline node earth model, do not do that before createExternalRadialNodes
!  because the latter also does a spline interpolation to get naterial parameters
!  and then deallocates again.
!
    call computeSplineEarthmodelCoefficients(nem,.fref.nem,attmode,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  set up ode system
!  SODE's intenv member must point to this%intenv because it  
!  (member nl) will be updated during integration (and not gem_intenv)
!
    call ray_intenv%createRayIntegrationEnvironment(nem,exnod,raytype,errmsg)      ! ray parameter is set later
    if (.level.errmsg == 2) goto 10
    call sode%createRayOdeSystem(ray_intenv,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  Establish propagator object for all ray integrations
!  Use real Bulirsch stepping (rbs)
!
    call createPropagateOde(prop,sode,ray_intenv,rbs,eps,secrecy,.false.,0.d0)
    if (.errlevel.prop == 2) then
       call printErrmsgPropagateOde(prop)
       call add(errmsg,2,'Problems with creating toroidal propagator',myname)
       call dealloc(prop)
       goto 10
    endif
!
!  start integration of ray equations for each turning point
!
    allocate(ray_table(2,nnod,ntp))                          ! Delta and T at node radii for all turning radii (slownesses)
    ray_table = 0.0
    do n = 1,ntp
       call getLayerIndexNodeEarthmodel(nem,rtp(n),nl,top)
       if (top == 1) nl = nl+1                                ! if tprn sits on layer boundary, take next higher layer
       if (raytype == 'P') then
          v = getPVelocitySplineEarthmodel(nl,rtp(n))
       else if (raytype == 'S') then
          v = getSVelocitySplineEarthmodel(nl,rtp(n))
       else
          call add(errmsg,2,'unknown ray type',myname)
          goto 10
       endif
       slowness(n) = rtp(n)/v
       call ray_intenv%setSlownessRayIntegrationEnvironment(slowness(n))
       call ray_intenv%setInterval(nl)
       ystart = [datan(sqrt(2.d0*eps)),slowness(n)*datan(sqrt(2.d0*eps))]                   ! start at rtp*(1+eps)
       call doPropagateOde(prop,rtp(n)*(1.d0+eps),rnod(nnod),drtp,ystart,[0.d0,0.d0])
       if (.errlevel.prop == 2) then
          call printErrmsgPropagateOde(prop)
          call add(errmsg,2,'Problems with doing toroidal propagation',myname)
          call dealloc(prop)
          return
       endif
   !
   !  read out values
   !
       call getIndexLimitsNodesPropagateOde(prop,jlu,jru)
       if (jru-jlu < 0) then
          call add(errmsg,2,'Either problems with node indices or no nodes in integration range',myname)
          call dealloc(prop)
          return
       endif
       allocate(yn(2,jru-jlu+1))
       call getRealSolutionAtNodesPropagateOde(prop,yn)
       if (secrecy <= secrecy_ray) then
          print *,'Ray for turning point at: ',rtp(n)
          do j = jlu,jru
             write(6,'(i6,f12.3,2f12.5)') j,rnod(j),yn(1:2,j-jlu+1)
          enddo
       endif
       ray_table(:,jlu:jru,n) = yn(:,1:jru-jlu+1)                   ! store ray in ray table
       deallocate(yn)
    enddo
!
!  write ray table and slowness to HDF
!
    allocate(d(ntp)); d = real(slowness)
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(fid,'rayParameters',arra,errmsg)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
    
    call arra%assoc3d(ray_table)
    call writeArrayHDFWrapper(fid,'rayTable',arra,errmsg)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc()   
!
!  deallocation and clean up
!
    call h5pclose_f(xferprp,ierr)                              ! close xfer property list
    if (ierr < 0) goto 1
    call h5fclose_f(fid,ierr)                                  ! close gfk-all file
    if (ierr < 0) goto 1 
    call h5close_f(ierr)                                       ! close Fortran interface
    if (ierr < 0) goto 1
!
    call dealloc(prop)
    call sode%deallocRayOdeSystem()
    deallocate(ray_table,rtp,slowness)
    call dealloc(nem)
    call deallocSplineEarthmodelCoefficients
    call dealloc(exnod)
    call dealloc(mpisup)
!
!  treat error messages
!
 1  if (ierr < 0) then
       call add(errmsg,2,'HDF-problem',myname)
       call print(errmsg)
       call abort(mpisup)
    endif
10  if (.level.errmsg == 2) then
       call print(errmsg)
       call abort(mpisup)
    endif
end program 
