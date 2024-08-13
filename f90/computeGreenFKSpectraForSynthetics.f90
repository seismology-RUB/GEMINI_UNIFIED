! ==============================================================================
!  Compute Green functions in the frequency-wavenumber domain
! ==============================================================================
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
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!  Main program to compute Green displacement-stress vectors in the fk-domain
!  for ensuing calculation of synthetic seismograms. Assumes receiver positions
!  at one fixed external radial node (typically surface, borehole, ocean, ocean bottom)
!  and computes Green functions for potential sources at all depth nodes.
!  The computation first determines unit jump Green functions
!  at all external nodes for a source at the receiver position and then applies the
!  reciprocity relation to obtain the unit jump Green function at the receiver position
!  for a source at any external node.
!-----------------------------------------------------------------------------
program computeGreenFKSpectraForSynthetics
        use hdf5
    use argumentParser
    use nodeEarthmodel
    use externalRadialNodes
    use geminiIntegrationEnvironment
    use splineEarthmodelCoefficients
    use fluidOdeSystem
    use toroidalMinors
    use spheroidalMinors
    use unitJumpGreenFunctions
    use sourceTerms
    use reciprocity
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
    type (gemini_integration_environment) :: gem_intenv
    type (fluid_ode_system) :: flode
    type (toroidal_minors) :: toromin
    type (spheroidal_minors) :: spheromin
    type (unit_jump_green_functions) :: ujgf
    type (error_message) :: errmsg,errmsg2,flerrmsg
    type (mpi_support) :: mpisup
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    integer, dimension(:), allocatable, target :: nwn
    integer, dimension(4) :: dsvmask_sph
    integer, dimension(2) :: dsvmask_tor
    integer(hid_t) :: fid,xferprp,grpid,dsp,dsetr,dseti
    integer(hsize_t), dimension(3) :: dimsall,count,offset,dimsslab   ! Dimensions of gfk-all, count and offset of wn-slab, dims of wnslab
    integer :: dsvstep,nwnmax,nkfmax,nf2,nf1,nf,nnod,jre,npart
    integer :: njtor,njsph,ncsph,nctor
    integer :: iwn,j,jf,i,n,iwa,isp,k,jfe
    integer :: myrank,numtasks,global,ierr,istyp,nsp,secrecy
    integer, dimension(:), pointer :: state,layer
    integer, dimension(:), allocatable, target :: id
    real, dimension(:), allocatable, target :: d
    double precision, dimension(:), pointer :: rnod
    double precision :: eps,df,dwn,f,fmin,frac,ommax,omre,sigma,p1lim,p2lim,wn,wn2,wnmax,re,rbot
    real, dimension(:,:,:), allocatable :: gfre,gfim
    double complex, dimension(:,:), allocatable :: gfsph,gftor,bgsph,bgtor,zsph,ztor
    double complex, dimension(:), allocatable :: vcomp,vrcomp
    character (len=3) :: crank
    character(len=max_length_string) :: parfile,errstr,attmode,sourcetype,identity,dsvbasename
    character (len=80), dimension(15) :: para_keywords
    character (len=80), dimension(4) :: attn_keywords
    character (len=80), dimension(2) :: stype_keywords
    character(len=34) :: myname = 'computeGreenFKSpectraForSynthetics'
    integer, parameter :: secrecy_gfk = 5                                ! print screen output if secrecy <= this value
    data para_keywords/'TLEN','FMAX','ATTENUATION_MODE', &
                   &   'SLOWNESS_LIMIT_1','SLOWNESS_LIMIT_2','STRESSFLAG','XLEN',&
                   &   'WAVENUMBER_MARGIN_FRACTION','EARTH_MODEL','DSVBASENAME','SOURCE_TYPE',&
                   &   'REC_NODE_RADIUS','ACCURACY','GLOBAL','EXTERNAL_NODES_FROM_FILE'/
    data attn_keywords/'ELASTIC','ATTENUATION_ONLY','DISPERSION_ONLY','ATTENUATION_AND_DISPERSION'/
    data stype_keywords/'FORCE','MOMENT'/
!-----------------------------------------------------------------------------
!  initialise MPI
!
    call new(mpisup)
    myrank = .myrank.mpisup
    write(crank,'(i3.3)') myrank
    numtasks = .numtasks.mpisup
!----------------------------------------------------------------------------
    call init(ap,myname,'Green displacement stress vectors in the fk-domain')
    call addPosarg(ap,'parfile','sval','gfk parameter file')
    call addOption(ap,'-s',.true.,'Secrecy level','ival','6')
    call addOption(ap,'-fmin',.true.,'lower frequency limit','dval','0.d0')
    call parse(ap)
    parfile = ap.sval.'parfile'
    secrecy = ap.ival.'-s'
    fmin = ap.dval.'-fmin'
    if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); stop; endif
    if (myrank == 0) call document(ap)
    call dealloc(ap)
!-----------------------------------------------------------------------------
    call new(errmsg,myname)
!-----------------------------------------------------------------------------
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
    call setXferprpCollectiveHDFWrapper(xferprp,errmsg)                             ! set data transfer prop to MPIO collective
    if (.level.errmsg == 2) goto 10
!-----------------------------------------------------------------------------
    call createKeywordsInputParameter(inpar,para_keywords)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 10
    if (myrank == 0) call printInputParameter(inpar)
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
!  check validity of source types, set istyp
!
    sourcetype = inpar.sval.'SOURCE_TYPE'
    if (sourcetype.equal.stype_keywords(1)) then
       istyp = 0
    else if (sourcetype.equal.stype_keywords(2)) then
       istyp = 1
    else
       istyp = -1
       call add(errmsg,2,'Invalid source type',myname)
       call print(errmsg); call abort(mpisup)
    endif
!-----------------------------------------------------------------------------
    call createNodeEarthmodel(nem,1,inpar.sval.'EARTH_MODEL', errmsg)
    if (.level.errmsg == 2) goto 10
    if (myrank == 0) call printNodeEarthmodel(nem)
!       
    if ((inpar.ival.'EXTERNAL_NODES_FROM_FILE') == 0) then
       call createExternalRadialNodes(exnod,1,parfile,nem,errmsg)
    else
       call createFromRadiiExternalRadialNodes(exnod,1,parfile,nem,errmsg)
    endif
    if (.level.errmsg == 2) goto 10
    if (myrank == 0) call printExternalRadialNodes(exnod)
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
!  set and print radius of receiver node
!  and external node information
!
    re = inpar.dval.'REC_NODE_RADIUS'
    jre = getNodeIdxFromRadiusExternalRadialNodes(exnod,re)
    re = getDoubleRadiusSelectedExternalRadialNodes(exnod,jre)    ! overwrite desired value by closest node
    if (myrank == 0) print *,'Radius of receiver: ',re
    nnod = .nnod.exnod
    rnod => getPointerDoubleRadiiExternalRadialNodes(exnod)
    state => getPointerStateExternalRadialNodes(exnod)
    layer => getPointerLayerExternalRadialNodes(exnod)
!-----------------------------------------------------------------------------
!  define variables derived from input data
!
    dsvbasename = inpar.sval.'DSVBASENAME'
    global = inpar.ival.'GLOBAL'
    eps = inpar.dval.'ACCURACY'
    df = 1.d0/(inpar.dval.'TLEN')
    nf2 = nint((inpar.dval.'FMAX')/df)+1
    nf1 = max(nint(fmin/df),2)
    nf = nf2-nf1+1
    ommax = mc_two_pid*(nf2-1)*df
    sigma = 5.d0*df                                           ! Do not change the 5.0 here. It is the best value!
    dwn = mc_two_pid/(inpar.dval.'XLEN')                      ! only used for regional applications
    p1lim = inpar.dval.'SLOWNESS_LIMIT_1'
    p2lim = inpar.dval.'SLOWNESS_LIMIT_2'
    frac = inpar.dval.'WAVENUMBER_MARGIN_FRACTION'
    wnmax = 0.5d0*ommax*(p1lim+p2lim)/(1.-frac)
    nwnmax = int(wnmax/dwn)+1
    if (global == 1) then
       nwnmax = nint(wnmax*(.rearth.nem)-0.5)+1               ! exact: kR = sqrt(l*(l+1)), for l > 10: kR = l+0.5, including l=0
    endif
!------------------------------------------------------------------------------------
!  set masks of required DSV-spectra: Spheroidal: (U,R,V,S), Toroidal: (W,T)
!
    dsvmask_sph = (/ 1,0,1,0 /)
    dsvmask_tor = (/ 1,0 /)
    if (inpar.lval.'STRESSFLAG') then; dsvstep = 1; else; dsvstep = 2; endif
    if (dsvstep == 1) then
       dsvmask_sph = (/ 1,1,1,1 /)
       dsvmask_tor = (/ 1,1 /)
    endif
!------------------------------------------------------------------------------------------
!  prepare wavenumber frequency count
!  precalculate total number of values in fk-spectrum 
!
    allocate(nwn(nf2))       ! frequencies up to nf2
    nwn = 0                  ! zero nwn-array, also used to count l-values in global case
    nkfmax=0                 ! total count of wavenumbers over all frequencies of this task
    do jf = nf1,nf2
       f = (jf-1)*df
       omre = 2.d0*mc_pi*f
       if (f <= 0.5*(nf2-1)*df) then                 ! two different cutoff lines for fk-spectrum
          wn2 = omre*p1lim+frac*wnmax
       else
          wn2 = 0.5*ommax*(p1lim+p2lim)+(omre-ommax)*p2lim+frac*wnmax
       endif
       if (global == 1) then
          nwn(jf) = nint(wn2*(.rearth.nem)-0.5)+1    ! including l=0
       else
          nwn(jf) = int(wn2/dwn)+1
       endif
       nkfmax = nkfmax+nwn(jf)                       ! count total amount of wavenumbers/l-values
       if (myrank == 0) write(6,'(f15.3,e15.3,i8)') f,wn2,nwn(jf)
    enddo
!-----------------------------------------------------------------------------
!  First write GFK meta data to file. They are identical for all of
!  the following GFK data files. Produce identity string written into all files,
!  to allow checking later for compatibilty.    
!
    call fdate(identity)
    call createFileParallelAccessHDFWrapper(trim(dsvbasename)//'.'//trim(sourcetype),fid,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  write identity code as attribute to file
!
    call writeStringAttributeHDFWrapper(fid,'identity',trim(identity),errmsg)
    if (.level.errmsg == 2) goto 10
!
! info about earth model and attenuation mode as meta data
!       
    call writeStringHDFWrapper(fid,'earthModelName',inpar.sval.'EARTH_MODEL',errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call writeStringHDFWrapper(fid,'attenuationMode',inpar.sval.'ATTENUATION_MODE',errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
!
!  write meta information to gfk-file, data sets of root group
!  Data set in root group: reals
!
    d = (/ real(.rearth.nem),real(re),real(sigma),real(df),real(dwn),real(.fref.nem)/)
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(fid,'reals',arra,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
!
!  Data set in root group: full external radial node instance
!
    call h5gcreate_f(fid,'externalNodes',grpid,ierr)
    if (ierr < 0) then; print *,'h5gcreate meta'; goto 1; endif
    call writeHDFExternalRadialNodes(exnod,grpid,errmsg)
    call h5gclose_f(grpid,ierr)
    if (ierr < 0) then; print *,'h5gclose externalNodes'; goto 1; endif
!  
!  Data set in root group: nwn(jf)-array
!
    call aria%assoc1d(nwn)
    call writeArrayHDFWrapper(fid,'numberOfWavenumbersForFrequency',aria,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call aria%deassoc()
!    
!  Data set in root group: integers
!  derivflag = 0, dsvstep = 2
!
    id = (/ nf1,nf2,nkfmax,nwnmax,2,0,numtasks,global /)
    call aria%assoc1d(id)
    call writeArrayHDFWrapper(fid,'integers',aria,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call aria%deassoc(); deallocate(id)
!---------------------------------------------------------------------------------
!  ncsph: number of components of spheroidal Green function
!  nctor: number of components of toroidal Green function
!
    if (state(jre) == 1) then                 ! Receiver in ocean
       ncsph = 2; nctor = 0                   ! U,R,(V = -1/(om**2*r*ro)*R,S=0),(W=0,T=0)
       dsvmask_tor = (/ 0,0 /)                ! no toroidal motion at receiver in ocean
    else                                      ! Receiver in solid
       ncsph = 4; nctor = 2                   ! U,R,V,S and W,T
    endif
!---------------------------------------------------------------------------------
!  njsph: max number of spheroidal source vectors
!  njtor: max number of toroidal source vectors
!  values could be smaller if source is in ocean (njsph=2,njtor=0), checked later    
!
    if (sourcetype.equal.'FORCE') then        ! Force
       njsph = 2; njtor = 1                   ! (0,z_21,0,0) and (0,0,0,z_42) and (0,z_61)
    else                                      ! Moment
       njsph = 4; njtor = 2                   ! (z_11,z_21,0,z_41),(0,z_22,0,z_42),(0,0,z_33,0),(0,0,0,z_44) and (z_11,0),(0,z_22)
    endif
!
!  total number of gfk-spectra
!
    nsp = sum(dsvmask_sph)*njsph+sum(dsvmask_tor)*njtor
!
    allocate(gfsph(ncsph,njsph))                 ! temp storage for gfsph and gftor for given wn,f and node
    allocate(gftor(max(nctor,1),njtor))          ! avoid compiler warnings if nctor = 0
    allocate(vcomp(njsph),vrcomp(njsph))         ! storage for V-vomponent and derivative if receiver in ocean
!-------------------------------------------------------------------------------
!  write data specific integers as attribute
!
    id = (/ istyp,ncsph,nctor,njsph,njtor /)
    call aria%assoc1d(id)   
    call writeArrayAttributeHDFWrapper(fid,'dataSpecificIntegers',aria,errmsg)
    if (.level.errmsg == 2) goto 10
    call aria%deassoc(); deallocate(id)
!
!  dimensions of the entire gfk dataset
!  and independent transfer mode    
!
    dimsall = (/nkfmax,nsp,nnod/)
    if (.level.errmsg == 2) then; call print(errmsg); call abort(mpisup); endif
    call h5screate_simple_f(size(dimsall),dimsall,dsp,ierr)
    if (ierr < 0) then; print *,'h5screate_simple dimsall'; goto 1; endif
    call h5dcreate_f(fid,"GreenFKSpectraReal",H5T_NATIVE_REAL,dsp,dsetr,ierr)
    if (ierr < 0) then; print *,'h5dcreate_f GreenFKSpectraReal'; goto 1; endif
    call h5dcreate_f(fid,"GreenFKSpectraImag",H5T_NATIVE_REAL,dsp,dseti,ierr)
    if (ierr < 0) then; print *,'h5dcreate_f GreenFKSpectraImag'; goto 1; endif
!
!  switch to independent data transfer mode
!
    call h5pclose_f(xferprp,ierr)
    call setXferprpIndependentHDFWrapper(xferprp,errmsg)
    if (.level.errmsg == 2) goto 10
!-------------------------------------------------------------------------------------------
!  Frequency loop, parallelized
!  Parallelize into processes doing either
!  jf = nf1 +(rank+(k-1)*numtasks)     for 1 <= k <= nf/(2*numtasks)
!  jf = nf2 -(rank+(k-1)*numtasks)     for 1 <= k <= nf/(2*numtasks)
!
    if (myrank == 0) then
       print *,'Start frequency loop'
       print *,'Min frequency index: ',nf1
       print *,'Max frequency index: ',nf2
       print *,'Number of frequencies: ',nf
       print *,'Frequency step: ',df
       print *,'Total number of values in fk-spectrum: ',nkfmax
       print *,'Number of Green FK spectra: ',nsp
       print *,'Number of external nodes: ',nnod
    endif
!
!  Example of distribution of 11 frequencies on 4 processes
!
!           1 2 3 4 5 6 7 8 9 10 11    -> frequencies
! proc 0: |     x x               x |
! proc 1: |   x     x          x    |
! proc 2: | x         x     x       |
! proc 3: |             x x         |
!
!  number of partitions of frequency range
!
    if (mod(nf,numtasks) > 0) then
       npart = nf/numtasks+1
    else
       npart = nf/numtasks
    endif
!
!  loop over partitions starting from upper end of frequency range
!
    jfe = nf2
    do k = 1,npart
       if (mod(k,2) > 0) then            ! Odd partition: proc 0 takes highest frequency of partition
          jf = jfe-myrank
       else                              ! Even partition: proc 0 takes lowest frequency of partition
          jf = jfe+myrank-numtasks+1
       endif
       jfe = jfe-numtasks                ! reduce highest frequency index by numtasks for next partition
       if (jf < nf1) exit                ! stop loop here, nothing left for this process
       f = (jf-1)*df
       omre = 2.d0*mc_pi*f
       call new(flerrmsg,myname)
       call computeSplineEarthmodelCoefficients(nem,f,attmode,flerrmsg)
       if (.level.flerrmsg == 2) then; call print(flerrmsg); call abort(mpisup); endif
       if (secrecy <= secrecy_gfk) then
          print *,'Frequency: ',f
          print *,'Max. wavenumber : ',wn2
          print *,'Number of wavenumbers: ',nwn(jf)
       endif
   !--------------------------------------------------------------------------
   !  allocate space for Green functions per frequency and initialize to zero    
   !
       allocate(gfre(nwn(jf),nsp,nnod),gfim(nwn(jf),nsp,nnod))
       gfre = 0.d0; gfim = 0.d0
   !--------------------------------------------------------------------------
   !  Wavenumber loop
   !  if regional caluclation: starts at 2, zeros are written to gfk file for wn=0
   !  if global starts at 1, l=0 included
   !
       iwa = 2
       if (global == 1) iwa = 1
       do iwn = iwa,nwn(jf)
!       do iwn = 1415,1415
          wn = (iwn-1)*dwn
          if (global == 1) wn = sqrt(dble((iwn-1)*iwn))/(.rearth.nem)
          if (secrecy <= secrecy_gfk) then
             if (global == 1) print *,'Harmonic degree: ',iwn-1
             if (global == 0) print *,'Wavenumber: ',wn
          endif
          call new(errmsg2,myname)
          call gem_intenv%createGemini(wn,omre,-sigma,nem,exnod,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !-------------------------------------------------------------------------------
      !  compute toroidal minors
      !
          call computeToroidalMinors(toromin,jre,jre,gem_intenv,eps,secrecy,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !-------------------------------------------------------------------------------
      !  compute spheroidal minors
      !
          call computeSpheroidalMinors(spheromin,jre,jre,gem_intenv,eps,secrecy,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !--------------------------------------------------------------------------------
      !  compute unit jump Green functions
      !
          call computeUnitJumpGreenFunctions(ujgf,toromin,spheromin,jre,gem_intenv,eps,secrecy,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !---------------------------------------------------------------------------------
      !  for V component if receiver is in ocean
      !
          call flode%createFluid(gem_intenv,errmsg2) 
          if (.level.errmsg2 == 2) goto 20
      !---------------------------------------------------------------------------------
      !  apply reprocity relations and multiply with source terms
      !
          if (secrecy <= secrecy_gfk) print *,'Compute gfsph using reciprocity'
          do n = 1,nnod
             gfsph = 0.d0; gftor = 0.d0; vcomp = 0.d0; vrcomp = 0.d0
             call spheroidalReciprocity(gem_intenv%dll1,re,rnod(n),state(n),state(jre),ujgf%ujgfsph(:,:,n),bgsph,errmsg2)
             if (.level.errmsg2 == 2) goto 20
             if (nctor > 0) call toroidalReciprocity(re,rnod(n),ujgf%ujgftor(:,:,n),bgtor)
             call oneNodeSourceTerms(sourcetype,gem_intenv%dll1,gem_intenv%omre,rnod(n),state(n),layer(n),zsph,ztor)
         !
         !  relations are valid for both force and moment source
         !  ncsph is either 2 or 4 depending on whether receiver is in ocean or solid
         !  njsph is either 2 or 4 depending on source type
         !  if njsph=4 and source in ocean, then gfsph(:,3:4) = 0    
         !
             if (state(n) == 1) then                                                 ! force/moment in fluid (zsph(2,2) both)
                gfsph(1:ncsph,1:2) = matmul(bgsph(1:ncsph,1:2),zsph)                      ! ncsph = 2 or 4, njsph = 2 / 4
             else                                                                     ! force/moment in solid (zsph(4,2) / zsph(4,4))
                gfsph(1:ncsph,1:njsph) = matmul(bgsph(1:ncsph,1:4),zsph)                  ! ncsph is either 2 or 4, njsph = 2 / 4
                if (nctor > 0) then
                   gftor(1:nctor,1:njtor) = matmul(bgtor(1:nctor,1:2),ztor)               ! nctor = 2; njtor = 1 / 2 (ztor(2,1) / ztor(2,2))
                endif
             endif
             if (state(jre) == 1) call flode%getHorcompComplexFluid(rnod(n),gfsph(1:ncsph,1:njsph),vcomp,vrcomp)
         !
         !  store into a single array: gf(iwn,isp,n)
         !
             isp = 0
             do i = 1,ncsph,dsvstep                          ! regular spheroidal components
                do j = 1,njsph
                   isp = isp+1
                   gfre(iwn,isp,n) = real(gfsph(i,j))
                   gfim(iwn,isp,n) = aimag(gfsph(i,j))
                enddo
             enddo
             if (state(jre) == 1) then                        ! also store V-component if receiver in ocean
                do j = 1,njsph
                   isp = isp+1
                   gfre(iwn,isp,n) = real(vcomp(j))
                   gfim(iwn,isp,n) = aimag(vcomp(j))
                enddo
             endif
             do i = 1,nctor,dsvstep                           ! regular toroidal components
                do j = 1,njtor
                   isp = isp+1
                   gfre(iwn,isp,n) = real(gftor(i,j))
                   gfim(iwn,isp,n) = aimag(gftor(i,j))
                enddo
             enddo
             deallocate(bgsph,zsph)
             if (allocated(bgtor)) deallocate(bgtor)
             if (allocated(ztor)) deallocate(ztor)
          enddo                                               ! end node loop
      !---------------------------------------------------------------------------------
      !  deallocation
      !
          call dealloc(toromin); call dealloc(spheromin)
          call dealloc(ujgf); call deallocFluid(flode)
          call dealloc(errmsg2); call dealloc(gem_intenv)          
       enddo                                                  ! end wavenumber loop
    ! ------------------------------------------------------------------------------------
    !  write wavenumber slab to file
    !  use independent mode for write to HDF
    !  to avoid waiting of processes for each other
    !  collective mode makes IO very slow   
    !
       dimsslab = (/nwn(jf),nsp,nnod/)
       offset = (/sum(nwn(1:jf-1)),0,0/)
       count = (/nwn(jf),nsp,nnod/)
       call arra%assoc3d(gfre)
       call writeArrayHDFWrapper(fid,'GreenFKSpectraReal',arra,errmsg,xferprp,&
            ds = dsetr,offset = offset,count = count)
       call arra%deassoc()
    !
       call arra%assoc3d(gfim)
       call writeArrayHDFWrapper(fid,'GreenFKSpectraImag',arra,errmsg,xferprp,&
            ds = dseti,offset = offset,count = count)
       call arra%deassoc()
    !
    !  clean up
    !
       deallocate(gfre,gfim)
       call deallocSplineEarthmodelCoefficients
       if (myrank == 0) write(6,'(i5,i6,$)') jf,nwn(jf)
       call dealloc(flerrmsg)
    enddo                                                     ! end frequency loop
!--------------------------------------------------------------------------------
!  clean up
!
    deallocate(gfsph,gftor,vcomp,vrcomp)
    call h5sclose_f(dsp,ierr)
    if (ierr < 0) goto 1
    call h5dclose_f(dsetr,ierr); call h5dclose_f(dseti,ierr)
    if (ierr < 0) goto 1
    call h5pclose_f(xferprp,ierr)
    if (ierr < 0) goto 1
    call h5fclose_f(fid,ierr)                                                           ! close gfk-all file
    if (ierr < 0) goto 1
    call h5close_f(ierr)                                                                ! close Fortran interface
    if (ierr < 0) goto 1
!
!  treat HFD error
!
 1  if (ierr < 0) then
       call add(errmsg,2,'HDF-problem',myname)
       call print(errmsg)
       call abort(mpisup)
    endif
!
!  deallocation
!
    deallocate(nwn)
    call dealloc(exnod); call dealloc(nem)
    call dealloc(inpar)
    call dealloc(mpisup)
    if (secrecy <= secrecy_gfk) call print(errmsg)
!
!  treat error messages
!
10  if (.level.errmsg == 2) then
       call print(errmsg)
       call abort(mpisup)
    endif
20  if (.level.errmsg2 == 2) then
       call print(errmsg2)
       call abort(mpisup)
    endif
 end program computeGreenFKSpectraForSynthetics
