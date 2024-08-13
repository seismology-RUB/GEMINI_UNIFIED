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
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  Main program to compute Green displacement-stress vectors in the fk-domain
!  for ensuing calculation of ASKI waveform sensitivity kernels. 
!  Computes Green functions at external nodes for a range of source nodes
!  used later for computing the incident wavefield.
!  In addition, computes a special Green function for a single force at one specific node,
!  used later for computing back-propagated wavefield.
!  Writes a separate output file for each source node.
!-----------------------------------------------------------------------------
program computeGreenFKSpectraForASKI
       use hdf5
    use argumentParser
    use nodeEarthmodel
    use externalRadialNodes
    use geminiIntegrationEnvironment
    use splineEarthmodelCoefficients
    use toroidalMinors
    use spheroidalMinors
    use toroidalOdeSystem
    use fluidOdeSystem
    use solidOdeSystem
    use unitJumpGreenFunctions
    use sourceTerms
    use reciprocity
    use dsvDerivatives
    use inputParameter
    use errorMessage
    use hdfWrapper
    use anyRankRealArray
    use anyRankIntegerArray
    use mpiSupport
    implicit none
    type (argument_parser) :: ap
    type (input_parameter) :: inpar
    type (node_earthmodel) :: nem
    type (external_radial_nodes) :: exnod
    type (gemini_integration_environment) :: gem_intenv
    type (toroidal_minors) :: toromin
    type (toroidal_ode_system) :: torode
    type (spheroidal_minors) :: spheromin
    type (solid_ode_system) :: sphode
    type (fluid_ode_system) :: flode
    type (unit_jump_green_functions) :: ujgf
    type (error_message) :: errmsg,errmsg2,flerrmsg
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    type (mpi_support) :: mpisup
    integer(hid_t), dimension(:), allocatable :: fid,dsetr,dseti,dsp
    integer(hid_t) :: xferprp,grpid,fidmeta
    integer(hsize_t), dimension(3) :: dimsall,count,offset,dimsslab
    integer, dimension(:), allocatable :: nctor,njsph,njtor,nsp,istyp
    integer, dimension(:), allocatable, target :: nwn,id
    integer, dimension(6) :: dsvmask_sph
    integer, dimension(3) :: dsvmask_tor
    integer :: iwn,j,jf,i,n,js,jsh,isp,iwa,jre,jfe,k,ios
    integer :: nwnmax,nkfmax,nf2,nf1,nf,nspmax,npart
    integer :: isourcetype,ierr,secrecy,nnod,nsnod,ncsph,global
    integer :: myrank,numtasks
    integer, dimension(:), allocatable :: jsnod
    integer, dimension(:), pointer :: state,layer
    logical :: jre_in_jsnod
    real, dimension(:), allocatable, target :: d
    real, dimension(:,:,:,:), allocatable :: gfre,gfim        !(nwn(jf),nsp,nnod,jsmax-jsmin+2)
    double precision, dimension(:), pointer :: rnod
    double precision, dimension(:), allocatable :: rsnod
    double precision :: eps,df,dwn,f,fmin,frac,ommax,omre,sigma,re
    double precision :: p1lim,p2lim,wn,wn2,wnmax,rbot
    double complex, dimension(:,:), allocatable :: zsph,ztor
    double complex, dimension(:,:), allocatable :: gfsph,gftor
    double complex, dimension(:,:), allocatable :: gfds,gfdt
    character (len=3) :: cjs
    character (len=30) :: identity
    character(len=max_length_string) :: parfile,errstr,attmode,sourcetype,styp
    character (len=80), dimension(16) :: para_keywords
    character (len=80), dimension(4) :: attn_keywords
    character (len=80), dimension(2) :: stype_keywords
    character(len=28) :: myname = 'computeGreenFKSpectraForASKI'
    integer, parameter :: secrecy_gfk = 5                                ! print screen output if secrecy <= this value
    data para_keywords/'TLEN','FMAX','ATTENUATION_MODE','SOURCE_TYPE', &
                   &   'SLOWNESS_LIMIT_1','SLOWNESS_LIMIT_2','XLEN','GLOBAL',&
                   &   'WAVENUMBER_MARGIN_FRACTION','EARTH_MODEL','DSVBASENAME',&
                   &   'ACCURACY','EXTERNAL_NODES_FROM_FILE','SOURCE_NODE_RADII','NSNOD','REC_NODE_RADIUS'/
    data attn_keywords/'ELASTIC','ATTENUATION_ONLY','DISPERSION_ONLY','ATTENUATION_AND_DISPERSION'/
    data stype_keywords/'FORCE','MOMENT'/
!-----------------------------------------------------------------------------
    call new(mpisup)
    myrank = .myrank.mpisup
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
!  Create external radial nodes either from blocking and stepping information
!  or by reading desired radii from file. In both cases, nodes on internal
!  interfaces are split into with a few meters distance from the interface.
!  Hence, there are no nodes directly on any interface.
!  In addition, truly elastic constants at the reference frequency are computed
!  at the nodes.    
!  Caution: uses splineEarthModelCoefficients and deallocs them.
!  Do not call the routines while an integration is running.    
!
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
!  set and print radius of receiver node for back propagation
!  and external node information
!
    re = inpar.dval.'REC_NODE_RADIUS'
    jre = getNodeIdxFromRadiusExternalRadialNodes(exnod,re)
    re = getDoubleRadiusSelectedExternalRadialNodes(exnod,jre)    ! overwrite desired value by closest node
    if (myrank == 0) print *,'Radius of receiver: ',re,' Node index of receiver: ',jre
    nnod = .nnod.exnod
    rnod => getPointerDoubleRadiiExternalRadialNodes(exnod)
    state => getPointerStateExternalRadialNodes(exnod)
    layer => getPointerLayerExternalRadialNodes(exnod)
!-----------------------------------------------------------------------------
!  check validity of source types, set istyp
!
    sourcetype = inpar.sval.'SOURCE_TYPE'
    if (sourcetype.equal.stype_keywords(1)) then
       isourcetype = 0
    else if (sourcetype.equal.stype_keywords(2)) then
       isourcetype = 1
    else
       isourcetype = -1
       call add(errmsg,2,'Invalid source type',myname)
       goto 10
    endif
!---------------------------------------------------------------------------
!  set source nodes, source node indices
!
    nsnod = inpar.ival.'NSNOD'
    allocate(rsnod(nsnod),jsnod(nsnod))
    rsnod = dvec(inpar,'SOURCE_NODE_RADII',nsnod,ios)
    if (ios /= 0) then
       call add(errmsg,2,'problems reading out SOURCE_NODE_RADII',myname)
       goto 10
    endif
    jre_in_jsnod = .false.
    do js = 1,nsnod
       jsnod(js) = getNodeIdxFromRadiusExternalRadialNodes(exnod,rsnod(js))
       if (jsnod(js) == jre) jre_in_jsnod = .true.        ! check that jre is among the jsnod
       rsnod(js) = rnod(jsnod(js))                        ! redefine rs-nodes at radial nodes
       if (myrank == 0) then
          print *,"Source nodes: ",jsnod(js),rsnod(js)
       endif
    enddo
    if (.not. jre_in_jsnod) then
       call add(errmsg,2,'receiver node must also be a source node',myname)
       goto 10
    endif
!-----------------------------------------------------------------------------
!  define variables derived from input data
!
    global = inpar.ival.'GLOBAL'
    eps = inpar.dval.'ACCURACY'
    df = 1.d0/(inpar.dval.'TLEN')
    nf2 = nint((inpar.dval.'FMAX')/df)+1
    nf1 = max(nint(fmin/df),1)+1
    nf = nf2-nf1+1
    ommax = mc_two_pid*(nf2-1)*df
!    ommax = mc_two_pid*nint(0.15d0/df)*df
    sigma = 5.d0*df                                      ! Do not change the 5.0 here. It is the best value!
    dwn = mc_two_pid/(inpar.dval.'XLEN')                 ! Only used for regional applications
    p1lim = inpar.dval.'SLOWNESS_LIMIT_1'
    p2lim = inpar.dval.'SLOWNESS_LIMIT_2'
    frac = inpar.dval.'WAVENUMBER_MARGIN_FRACTION'
    wnmax = 0.5d0*ommax*(p1lim+p2lim)/(1.-frac)
    nwnmax = int(wnmax/dwn)+1
    if (global == 1) then
       nwnmax = nint(wnmax*(.rearth.nem)-0.5)+1          ! exact: kR = sqrt(l*(l+1)), for l > 10: kR = l+0.5, including l=0
    endif
!------------------------------------------------------------------------------------
!  set masks of required DSV-spectra: Spheroidal: (U,R,V,S,UP,VP), Toroidal: (W,T,WP)
!
    ncsph = 4                                            ! always 4 because receiver may be in solid
    dsvmask_sph = (/ 1,0,1,0,1,1 /)
    dsvmask_tor = (/ 1,0,1 /)
!------------------------------------------------------------------------------------------
!  prepare wavenumber frequency count
!  precalculate total number of values in fk-spectrum 
!
    allocate(nwn(nf2))       ! frequencies up to nf2
    nwn = 0                  ! zero nwn-array, also used to count l-values in global case
    nkfmax = 0                 ! total count of wavenumbers over all frequencies of this task
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
    enddo
!-----------------------------------------------------------------------------
!  First write GFK meta data to file. They are identical for all of
!  the following GFK data files. Produce identity string written into all files,
!  to allow checking later for compatibilty.    
!
    call fdate(identity)
    call createFileParallelAccessHDFWrapper((inpar.sval.'DSVBASENAME')+'.meta',fidmeta,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  write identity code as attribute to file
!
    call writeStringAttributeHDFWrapper(fidmeta,'identity',trim(identity),errmsg)
    if (.level.errmsg == 2) goto 10
!
! info about earth model and attenuation mode as meta data
!       
    call writeStringHDFWrapper(fidmeta,'earthModelName',inpar.sval.'EARTH_MODEL',errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call writeStringHDFWrapper(fidmeta,'attenuationMode',inpar.sval.'ATTENUATION_MODE',errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
!
!  write meta information to gfk-file, data sets of root group
!  Data set in root group: reals
!
    d = (/ real(.rearth.nem),real(re),real(sigma),real(df),real(dwn),real(.fref.nem)/)
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(fidmeta,'reals',arra,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
!
!  Data set in root group: full external radial node instance
!
    call h5gcreate_f(fidmeta,'externalNodes',grpid,ierr)
    if (ierr < 0) then; print *,'h5gcreate meta'; goto 1; endif
    call writeHDFExternalRadialNodes(exnod,grpid,errmsg)
    call h5gclose_f(grpid,ierr)
    if (ierr < 0) then; print *,'h5gclose externalNodes'; goto 1; endif
!  
!  Data set in root group: source node radii
!
    allocate(d(nsnod)); d = rsnod
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(fidmeta,'sourceNodeRadii',arra,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
!  
!  Data set in root group: nwn(jf)-array
!
    call aria%assoc1d(nwn)
    call writeArrayHDFWrapper(fidmeta,'numberOfWavenumbersForFrequency',aria,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call aria%deassoc()
!    
!  Data set in root group: integers
!  derivflag = 1, dsvstep = 2
!
    id = (/ nf1,nf2,nkfmax,nwnmax,2,1,numtasks,global /)
    call aria%assoc1d(id)
    call writeArrayHDFWrapper(fidmeta,'integers',aria,errmsg,xferprp)
    if (.level.errmsg == 2) goto 10
    call aria%deassoc(); deallocate(id)
!
!  close file
!
    call h5fclose_f(fidmeta,ierr)
    if (ierr < 0) then; print *,'h5fclose meta'; goto 1; endif
!-----------------------------------------------------------------------------
!  Open nsnod HDF-files
!  There are nsnod+1 source nodes (receiver node included)
!
    allocate(fid(nsnod+1),dsetr(nsnod+1),dseti(nsnod+1),dsp(nsnod+1))
    allocate(njsph(nsnod+1),nctor(nsnod+1),njtor(nsnod+1),nsp(nsnod+1),istyp(nsnod+1))
!
!  Loop over source nodes
!
    do js = 1,nsnod+1
       if (js == nsnod+1) then              ! set source type
          styp = 'FORCE'
          istyp(js) = 0
          jsh = jre
          write(cjs,'(i3.3)') js
       else
          styp = sourcetype
          istyp(js) = isourcetype
          jsh = jsnod(js)                    ! true source node index
          write(cjs,'(i3.3)') js
       endif
   !
   !  set njsph,nctor and njtor according to styp and state of node
   !
       if (styp.equal.'FORCE') then
          if (state(jsh) == 1) then              ! force in ocean
             njsph(js) = 2
             nctor(js) = 0; njtor(js) = 0
          else                                   ! force in solid
             njsph(js) = 2
             nctor(js) = 2; njtor(js) = 1             
          endif
       else
          if (state(jsh) == 1) then              ! moment in ocean
             njsph(js) = 2
             nctor(js) = 0; njtor(js) = 0
          else                                   ! moment in solid
             njsph(js) = 4
             nctor(js) = 2; njtor(js) = 2
          endif
       endif
       nsp(js) = sum(dsvmask_sph)*njsph(js)+sum(dsvmask_tor)*njtor(js)
   !---------------------------------------------------------------------------
   !  open HDF file for soure node js
   !
       call createFileParallelAccessHDFWrapper((inpar.sval.'DSVBASENAME')+'.'+styp+'.'+cjs,fid(js),errmsg)
       if (.level.errmsg == 2) goto 10
   !
   !  write identity code as attribute to file
   !
       call writeStringAttributeHDFWrapper(fid(js),'identity',trim(identity),errmsg)
       if (.level.errmsg == 2) goto 10
   !
   !  write data specific integers as attribute
   !
       id = (/ istyp(js),ncsph,nctor(js),njsph(js),njtor(js) /)
       call aria%assoc1d(id)   
       call writeArrayAttributeHDFWrapper(fid(js),'dataSpecificIntegers',aria,errmsg)
       if (.level.errmsg == 2) goto 10
       call aria%deassoc(); deallocate(id)
   !
   !  create dataspaces and datasets for the gfk-all output file in 3D (gfk3D(nkfmax,nsp,nnod)) 
   !
       dimsall = (/nkfmax,nsp(js),nnod/)
       call h5screate_simple_f(3,dimsall,dsp(js),ierr)                    ! create a simple data space for gfk-all data
       if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
       call h5dcreate_f(fid(js),"GreenFKSpectraReal",H5T_NATIVE_REAL,dsp(js),dsetr(js),ierr)
       if (ierr < 0) then; print *,'h5dcreate_f GreenFKSpectraReal'; goto 1; endif
       call h5dcreate_f(fid(js),"GreenFKSpectraImag",H5T_NATIVE_REAL,dsp(js),dseti(js),ierr)
       if (ierr < 0) then; print *,'h5dcreate_f GreenFKSpectraImag'; goto 1; endif
    enddo                                                                 ! end source node loop
    nspmax = maxval(nsp)
!
!  switch to independent data transfer mode
!
    call h5pclose_f(xferprp,ierr)
    call setXferprpIndependentHDFWrapper(xferprp,errmsg)
    if (.level.errmsg == 2) goto 10
!-------------------------------------------------------------------------------------------
!  Frequency loop, parallelized
!
    if (myrank == 0) then
       print *,'Start frequency loop'
       print *,'Min frequency index: ',nf1
       print *,'Max frequency index: ',nf2
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
   !-----------------------------------------------------------------------------  
   !  allocate space for Green functions per frequency and initialize to zero
   !
       allocate(gfre(nwn(jf),nspmax,nnod,nsnod+1),gfim(nwn(jf),nspmax,nnod,nsnod+1))
       gfre = 0.d0; gfim = 0.d0
   !--------------------------------------------------------------------------
   !  Wavenumber loop, starts at 2, zeros are written to gfk file for iwn=1
   !
       iwa = 2
       if (global == 1) iwa = 1
       do iwn = iwa,nwn(jf)
          wn = (iwn-1)*dwn
          if (global == 1) wn = sqrt(dble((iwn-1)*iwn))/(.rearth.nem)
          if (secrecy <= secrecy_gfk) print *,'Wavenumber: ',wn,iwn-1
          call new(errmsg2,myname)
          call gem_intenv%createGemini(wn,omre,-sigma,nem,exnod,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !-------------------------------------------------------------------------------
      !  compute toroidal minors
      !
          call computeToroidalMinors(toromin,minval(jsnod),maxval(jsnod),gem_intenv,eps,secrecy,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !-------------------------------------------------------------------------------
      !  compute spheroidal minors
      !
          call computeSpheroidalMinors(spheromin,minval(jsnod),maxval(jsnod),gem_intenv,eps,secrecy,errmsg2)
          if (.level.errmsg2 == 2) goto 20
      !-------------------------------------------------------------------------------
      !  prepare Green function computation
      !
          call sphode%createSolid(gem_intenv,errmsg2)        ! for derivatives
          if (.level.errmsg2 == 2) goto 20
          call torode%createToroidal(gem_intenv,errmsg2)     ! for derivatives
          if (.level.errmsg2 == 2) goto 20
          call flode%createFluid(gem_intenv,errmsg2)         ! for derivatives
          if (.level.errmsg2 == 2) goto 20
      !-------------------------------------------------------------------------------
      !  Loop over source nodes (including force at receiver position)
      !
          do js = 1,nsnod+1
             if (js == nsnod+1) then     ! force at receiver position (node index = jre)
                styp = 'FORCE'
                jsh = jre
             else
                styp = sourcetype        ! source at external node js
                jsh = jsnod(js)
             endif
         !
         !  compute unit jump Green functions
         !
             call computeUnitJumpGreenFunctions(ujgf,toromin,spheromin,jsh,gem_intenv,eps,secrecy,errmsg2)
             if (.level.errmsg2 == 2) goto 20
         !
         !  calculate source terms
         !
             call oneNodeSourceTerms(styp,gem_intenv%dll1,gem_intenv%omre,rnod(jsh),state(jsh),layer(jsh),zsph,ztor)
         !
         !  multiply unit jump Green functions with source terms
         !  relations are valid for both force and moment source
         !  ncsph is always 4
         !  njsph is either 2 or 4 depending on source type
         !  number of unit jumps is 2 for source in ocean and 4 for source in solid (0 and 2 for toroidal)    
         !
             allocate(gfsph(ncsph,njsph(js)),gfds(ncsph,njsph(js)))
             allocate(gftor(max(nctor(js),1),max(njtor(js),1)),gfdt(max(nctor(js),1),max(njtor(js),1)))    ! avoid compiler warning
             do n = 1,nnod
                gfsph = 0.d0; gftor = 0.d0
                gfds = 0.d0; gfdt = 0.d0
                if (state(jsh) == 1) then                                                ! force/moment in the ocean (zsph(2,2),njsph=2)
                   if (state(n) == 1) then
                      gfsph(1:2,:) = matmul(ujgf%ujgfsph(1:2,1:2,n),zsph)                ! receiver in ocean
                      call computeDsvDerivatives(gfsph(1:2,:),rnod(n),gem_intenv,flode,gfds(1:2,:),errmsg2)
                      if (.level.errmsg2 == 2) goto 20
                      call flode%getHorcompComplexFluid(rnod(n),gfsph(1:2,:),gfsph(3,:),gfds(3,:))
                   else
                      gfsph(:,:) = matmul(ujgf%ujgfsph(:,1:2,n),zsph)                    ! receiver in solid
                      call computeDsvDerivatives(gfsph,rnod(n),gem_intenv,sphode,gfds,errmsg2)
                      if (.level.errmsg2 == 2) goto 20
                   endif
                else                                                                    ! force/moment in solid (zsph(4,2) / zsph(4,4))
                   if (state(n) == 1) then                                              !                       (ztor(2,1) / ztor(2,2))
                      gfsph(1:2,:) = matmul(ujgf%ujgfsph(1:2,:,n),zsph)                  ! receiver in ocean
                      call computeDsvDerivatives(gfsph(1:2,:),rnod(n),gem_intenv,flode,gfds(1:2,:),errmsg2)
                      if (.level.errmsg2 == 2) goto 20
                      call flode%getHorcompComplexFluid(rnod(n),gfsph(1:2,:),gfsph(3,:),gfds(3,:))
                   else
                      gfsph(:,:) = matmul(ujgf%ujgfsph(:,:,n),zsph)                      ! receiver in solid
                      gftor(:,:) = matmul(ujgf%ujgftor(:,:,n),ztor)
                      call computeDsvDerivatives(gfsph,rnod(n),gem_intenv,sphode,gfds,errmsg2)
                      if (.level.errmsg2 == 2) goto 20
                      call computeDsvDerivatives(gftor,rnod(n),gem_intenv,torode,gfdt,errmsg2)
                      if (.level.errmsg2 == 2) goto 20
                   endif
                endif
          !
          !  store into a single array: gf(iwn,isp,n,js)
          !
                isp = 0
                do i = 1,ncsph,2                            ! U,V
                   do j = 1,njsph(js)
                      isp = isp+1
                      gfre(iwn,isp,n,js) = real(gfsph(i,j))
                      gfim(iwn,isp,n,js) = aimag(gfsph(i,j))
                   enddo
                enddo
                do i = 1,ncsph,2                            ! dU/dr, dV/dr                           
                   do j = 1,njsph(js)
                      isp = isp+1
                      gfre(iwn,isp,n,js) = real(gfds(i,j))
                      gfim(iwn,isp,n,js) = aimag(gfds(i,j))
                   enddo
                enddo
                do j = 1,njtor(js)                           ! W
                   isp = isp+1
                   gfre(iwn,isp,n,js) = real(gftor(1,j))
                   gfim(iwn,isp,n,js) = aimag(gftor(1,j))
                enddo
                do j = 1,njtor(js)                           ! dW/dr
                   isp = isp+1
                   gfre(iwn,isp,n,js) = real(gfdt(1,j))
                   gfim(iwn,isp,n,js) = aimag(gfdt(1,j))
                enddo
             enddo                                 ! end of receiver node loop
             deallocate(gfsph,gftor,gfds,gfdt)
             deallocate(zsph)
             if (allocated(ztor)) deallocate(ztor)
             call dealloc(ujgf)
          enddo                            !  end of source node loop
      !-------------------------------------------------------------------------------    
      !  deallocation
      !
          call dealloc(toromin); call dealloc(spheromin)
          call sphode%deallocSolid(); call torode%deallocToroidal(); call flode%deallocFluid()
          call dealloc(errmsg2)
          call dealloc(gem_intenv)
       enddo                         !  end of wavenumber loop
   ! ------------------------------------------------------------------------------------
   !  write wavenumber slab to files
   !
       do js = 1,nsnod+1
          dimsslab = (/nwn(jf),nsp(js),nnod/)
          offset = (/sum(nwn(1:jf-1)),0,0/)
          count = (/nwn(jf),nsp(js),nnod/)
          call arra%assoc3d(gfre(:,1:nsp(js),:,js))
          call writeArrayHDFWrapper(fid(js),'GreenFKSpectraReal',arra,errmsg,xferprp,&
               ds = dsetr(js),offset = offset,count = count)
          if (.level.errmsg == 2) goto 10
          call arra%deassoc()
    !
          call arra%assoc3d(gfim(:,1:nsp(js),:,js))
          call writeArrayHDFWrapper(fid(js),'GreenFKSpectraImag',arra,errmsg,xferprp,&
               ds = dseti(js),offset = offset,count = count)
          if (.level.errmsg == 2) goto 10
          call arra%deassoc()
       enddo
   !
   !  clean up
   !
       deallocate(gfre,gfim)
       call deallocSplineEarthmodelCoefficients
       if (myrank == 0) write(6,'(i5,i6,$)') jf,nwn(jf)
       call dealloc(flerrmsg)
    enddo                                          !  End of frequency loop
!--------------------------------------------------------------------------------
!  clean up
!
    call h5pclose_f(xferprp,ierr)                                                        ! close xfer property list
    if (ierr < 0) goto 1
    do js = 1,nsnod+1
       call h5dclose_f(dsetr(js),ierr)                                                   ! close gfk-all dataset for real part
       if (ierr < 0) goto 1
       call h5dclose_f(dseti(js),ierr)                                                   ! close gfk-all dataset for imag part
       if (ierr < 0) goto 1
       call h5fclose_f(fid(js),ierr)                                                     ! close gfk-all file
       if (ierr < 0) goto 1 
       call h5sclose_f(dsp(js),ierr)                                                     ! close gfk-all dataspace
       if (ierr < 0) goto 1
    enddo
    call h5close_f(ierr)                                                                 ! close Fortran interface
    if (ierr < 0) goto 1
    deallocate(fid,dsetr,dseti,dsp)
    deallocate(njsph,nctor,njtor,nsp)
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
    deallocate(nwn,rsnod,jsnod)
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
 end program computeGreenFKSpectraForASKI
