! ===============================================================================
!  Main program for the computation of synthetic seismograms for SPECFEM coupling
! ===============================================================================
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
!---------------------------------------------------------------------
!  Compute synthetic seismograms from Green FK spectra for ASKI
!  at boundary point of SPECFEM computazional box.
!  User provides boundary point coordinates, event file, filter
!  specifications for each event.
!  Seismograms (displacements and normal stresses) are written to HDF
!  at each boundary point.
!-----------------------------------------------------------------------------
program computeSyntheticsForSpecfemCoupling
    use hdf5
    use mathConstants
    use argumentParser
    use greenFKSpectra
    use displacementSpectra              ! NCOMP defined here
    use singleFrequencyDisplacement
    use harmonicDegreeSums
    use legendreFunctions
    use readEventStationFile
    use asciiDataIO
    use asciiSynseisIO
    use frequencyTime
    use eventFilter
    use errorMessage
    use inputParameter
    use fileUnitHandler
    use string
    use synseisHDF
    use boundaryPoints
    use travelTimes
    use mpiSupport
    implicit none
    type (mpi_support) :: mpisup
    type (argument_parser) :: ap
    type (green_fk_spectra) :: gfk
    type (error_message) :: errmsg,errmsg1,errmsg2
    type (input_parameter) :: inpar
    type (seismic_event) :: event
    type (seismic_event_list) :: event_list
    type (legendre_functions), dimension(:), allocatable :: alf
    type (legendre_functions) :: alfb
    type (boundary_points) :: bp
    type (raytable_travel_times) :: raytable
    type (single_boundary_point) :: cbp
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    real, dimension(:), allocatable :: d,urs,stf
    real, dimension(:,:), pointer :: gfkspr => null(),gfkspi => null()
    real, dimension(:,:), allocatable :: seis,uspre,uspim
    real :: dtstf
    double complex, dimension(:), allocatable :: evfilter,tsfil
    double complex, dimension(:,:), allocatable :: zsp,tensionsp
    double complex, dimension(:,:,:), allocatable :: strainsp
    double precision, dimension(:), pointer :: evfilspecs,rnod
    double precision, dimension(:), allocatable :: fomt,phi
    double precision, dimension(:), pointer :: uniquedis
    double precision :: rs,re,tanf,tanfstf,dt,dtp,tred,ttmin,ttmax,tsbuf,tend,twinlen,tscatlen,tlen
    double precision :: phic,phis,thetas,thetac,delta,xi,delmin,delmax,taperlen
    integer, dimension(:), allocatable :: point2dis
    integer, dimension(:), allocatable :: id
    integer(hid_t) :: fidseis,evgrid,fidmeta,fidgfk,xferprp
    integer(hsize_t), dimension(:), allocatable :: dimsall
    integer(hsize_t), dimension(:), allocatable :: dsp,dset,procid,dspcoor,dsetcoor,dsetre,dsetim
    integer(hsize_t), dimension(:), allocatable :: offset,countvt
    integer :: jeprev,je,js,ic,jf,nsamp,idis,nd,off,nstrue,nstf
    integer :: ierr,ibp,myrank,numtasks,ip,iproc,nproc
    logical :: compute_spectra, autoTimeshift,convolve_stf,single_event
    character(len=max_length_string) :: eventfile,boundaryfile,raytablefile, &
              parfile,synseispath,gfkmetafile,gfkbasename,evfiltype,dsvbasename, &
              source_type,evsel,stfpath,dummy1,dummy2,stffile,synfile,winlenfile, &
              scatlenfile,path_measured_seis,timing_file
    character(len=max_length_string) :: errstr
 !  character (len=10) :: time
    character (len=3) :: cjs
    character (len=6) :: cip
    character (len=27) :: myname = 'computeSyntheticsForSpecfemCoupling'
    character (len=80), dimension(14) :: par_keys
!
!  keywords for input parameters
!
    data par_keys/'PATH_SYNTHETICS','FILE_EVENT_LIST','EVENT_FILTER_TYPE','EVENT_FILTER_SPECS',&
         'EXTERNAL_NODES_HDF_FILE','DSVBASENAME','SOURCE_TYPE','COMPUTE_SPECTRA',&
         'FILE_RAY_TABLE','TIME_SHIFT_BUFFER','AUTOMATIC_TIMESHIFT','PHASE_WINDOW_LENGTH',&
         'END_TAPER_LENGTH','PATH_MEASURED_SEIS'/
!--------------------------------------------------------------------------------------------------
    call init(ap,myname,'Compute synthetic seismograms for SPECFEM coupling using Green FK-spectra')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call addOption(ap,'-stf',.false.,'Convolve with source time function. Assume name of file is stf.txt')
    call addOption(ap,'-single',.true.,'Run for specified event','sval','None')
    call parse(ap)
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    parfile = ap.sval.'parfile'
    convolve_stf = ap.optset.'-stf'
    evsel = ap.sval.'-single'
    single_event = ap.optset.'-single'
!-----------------------------------------------------------------------------
    call new(mpisup)
    myrank = .myrank.mpisup
    numtasks = .numtasks.mpisup
    if (myrank == 0) call document(ap)
    call dealloc(ap)
!
    call new(errmsg,myname)
!-------------------------------------------------------------
!  read input parameters from parfile
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 1
    source_type = inpar.sval.'SOURCE_TYPE'
    dsvbasename = inpar.sval.'DSVBASENAME'
    gfkmetafile = dsvbasename+".meta"
    gfkbasename = dsvbasename+"."+source_type
    synseispath = inpar.sval.'PATH_SYNTHETICS'
    eventfile = inpar.sval.'FILE_EVENT_LIST'
    boundaryfile = inpar.sval.'EXTERNAL_NODES_HDF_FILE'
    evfiltype = inpar.sval.'EVENT_FILTER_TYPE'
    compute_spectra = inpar.lval.'COMPUTE_SPECTRA'
    evfilspecs => dvecp(inpar,'EVENT_FILTER_SPECS',4,ierr)
    raytablefile = inpar.sval.'FILE_RAY_TABLE'
    tsbuf = inpar.dval.'TIME_SHIFT_BUFFER'
    taperlen = inpar.dval.'END_TAPER_LENGTH'
    autoTimeshift = inpar.lval.'AUTOMATIC_TIMESHIFT'
    twinlen = inpar.dval.'PHASE_WINDOW_LENGTH'
    path_measured_seis = inpar.sval.'PATH_MEASURED_SEIS'
    stfpath = path_measured_seis
    if (ierr /= 0) then
       call add(errmsg,2,'Input for event filter specs cannot be read',myname)
       goto 10
    endif
    if (compute_spectra .and. convolve_stf) then
       call add(errmsg,2,'Source wavelet convolution not compatible with compute spectra',myname)
       goto 10
    endif
    if (myrank == 0) call printInputParameter(inpar)
    call dealloc(inpar)
!----------------------------------------------------------------
!  open HDF environment    
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read event file
!  use ASKI convention for event file
!
    call createEventListFromEventFile(eventfile,1,'ASKI_events',event_list,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  check whether csys in event file is S
!
    if (.csys.event_list /= 'S') then
       call add(errmsg,2,'CSYS should be spherical',myname)
       goto 10
    endif
!
!  read SPECFEM boundary points, get box center    
!
    call readBoundaryPoints(bp,trim(boundaryfile),errmsg)
    if (.level.errmsg == 2) goto 1
    call getBoxCenterBoundaryPoints(bp,thetac,phic)
    nd = getSidesNumTotalBoundaryPoints(bp)
    if (mod(nd,numtasks) /= 0) then
       print *,'NBPVERT = ',nd,' NUMTASKS = ',numtasks
       call add(errmsg,2,'numtasks should be a divisor of nbpvert', myname)
       goto 10
    endif
    nd = getBottomNumTotalBoundaryPoints(bp)
    if (mod(nd,numtasks) /= 0) then
       print *,'NBOT = ',nd,' NUMTASKS = ',numtasks
       call add(errmsg,2,'numtasks should be a divisor of nbot', myname)
       goto 10
    endif
!
!  read GFK meta data, equal for all events and source depths
!
    call openFileParallelAccessHDFWrapper(gfkmetafile,fidmeta,errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'identity',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'attenuationMode',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'reals',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'integers',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'numberOfWavenumbersForFrequency',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'externalNodes',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'sourceNodeRadii',errmsg)
    if (.level.errmsg == 2) goto 10
    call h5fclose_f(fidmeta,ierr)
    if (ierr < 0) then; errstr = 'failed to close GFK meta file'; goto 1; endif
!
!  read ray table file
!
    if (numtasks > 1) then
       call readTableTravelTimes(raytable,raytablefile,errmsg,parallel=.true.)
    else
       call readTableTravelTimes(raytable,raytablefile,errmsg,parallel=.false.)
    end if
    if (.level.errmsg == 2) goto 10
!
!  Calculate event filter, only if seismograms are computed
!
    if (.not. compute_spectra) then
       call computeEventFilter(evfiltype,evfilspecs,gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,evfilter,errmsg)
       if (.level.errmsg == 2) goto 10
    endif
!
!  node radii
!
    rnod => getPointerDoubleRadiiExternalRadialNodes(gfk%exnod)
!
!  set frequency index limits and compute dt
!  desired dt = 0.25*tlen/(2*nf2)
!
    if (compute_spectra) then
       allocate(uspre(2*NCOMP,gfk%nf2),uspim(2*NCOMP,gfk%nf2))
    else
       dtp = 0.25d0/(2.d0*gfk%nf2*gfk%df)
       call newNsampDtFrequencyTime(gfk%df,dtp,nsamp,dt)
       allocate(urs(nsamp))
    endif
    tlen = 1.d0/gfk%df
!
    nproc = getNumProcBoundaryPoints(bp)       ! number of SPECFEM processes
!--------------------------------------------------------------------------------------
!  Loop over events
!--------------------------------------------------------------------------------------
    do while (nextEventSeismicEventList(event_list,event))
       call new(errmsg1,myname)
       if (single_event) then                                             ! selected event desired
          if (.not. ((.evid.event).equal.evsel)) cycle
       endif
   !
   !  open file for synthetics output
   !
       if (compute_spectra) then
          synfile = trim(synseispath)//trim(.evid.event)//'_spec.hdf'
       else
          synfile = trim(synseispath)//trim(.evid.event)//'_seis.hdf'
          timing_file = trim(synseispath)//'timing_'//trim(.evid.event)
       endif
       call createFileParallelAccessHDFWrapper(trim(synfile),fidseis,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       call setXferprpIndependentHDFWrapper(xferprp,errmsg)
       if (.level.errmsg == 2) goto 10
       if (compute_spectra) then
          call writeStringAttributeHDFWrapper(fidseis,'SyntheticsType','spectra',errmsg1)
       else
          call writeStringAttributeHDFWrapper(fidseis,'SyntheticsType','seismograms',errmsg1)
       endif
       if (.level.errmsg1 == 2) goto 11
   !
   !  write nproc as attribute to output file
   !
       id = [nproc]; call aria%assoc1d(id)
       call writeArrayAttributeHDFWrapper(fidseis,"nproc",aria,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       call aria%deassoc()
   !
       if (myrank == 0) call printSeismicEvent(event)
       tanf = .tanfdp.(.otime.event)                          ! origin time of event in seconds after midnight
       thetas = 0.5d0*mc_pid-.slatrad.event                   ! geographical source coordinates
       phis = .slonrad.event
       rs = gfk%rearth-(.sdepth.event)                        ! true source radius
       call getForceMomentSeismicEvent(event,fomt)            ! force or moment components depending on istyp
   !
   !  read source time function for this event
   !
       if (convolve_stf) then
          stffile = trim(stfpath)//'/'//trim(.evid.event)//'/stf.txt'
          call readAsciiSynseisIO(2,trim(stffile),nstf,tanfstf,dtstf,dummy1,dummy2,stf,ierr,op = .true.,cl = .true.)
          if (ierr /= 0) then; errstr = 'cannot read source wavelet'; goto 1; endif
          if (abs(dtstf-dt) > 1.d-6) then
             call add(errmsg1,2,'DT of source wavelet differs from that of synthetice',myname)
             goto 11
          endif
          if (tanfstf > 0.d0) then
             call add(errmsg1,2,'expect zero or negative start time for source wavelet',myname)
             goto 11
          endif
       !
       !  do no longer integrate source wavelet here, use source wavelet itself (source wavelet is moment rate)
       !  this implies that SPECFEM calculates particle velocity instead of particle displacement
       !  hence semd-synthetics contain particle velocity
       !
       !  read true phase window length from file if convolution with the STF is desired
       !  also read scattering window length which is later written to HDF
       !
          winlenfile = trim(stfpath)//'/'//trim(.evid.event)//'/phasewinlen.txt'
          open(3,file = trim(winlenfile),status = 'old',iostat = ierr)
          if (ierr /= 0) then; errstr = 'cannot open winlenfile: '//trim(winlenfile); goto 1; endif
          read(3,*) twinlen
          close(3)
          scatlenfile = trim(stfpath)//'/'//trim(.evid.event)//'/scatwinlen.txt'
          open(4,file = trim(scatlenfile),status = 'old',iostat = ierr)
          if (ierr /= 0) then; errstr = 'cannot open scatlenfile: '//trim(scatlenfile); goto 1; endif
          read(4,*) tscatlen
          close(4)
       endif
   !
   !  get closest source node to true source radius
   !      
       js = minloc(abs(gfk%rsnod-rs),1)
       if (myrank == 0) then
          write(6,'(a,f12.3)') 'True source radius: ',rs
          write(6,'(a,i5,f12.3)') 'Closest source node: ',js,gfk%rsnod(js)
       endif
       rs = gfk%rsnod(js)
   !
   !  get source info for travel time calculation
   !
       call setSourceInfoTravelTimes(raytable,rs,errmsg1)
       if (.level.errmsg1 == 2) goto 11
   !
   !  open associated GFK data file
   !
       write(cjs,'(i3.3)') js
       call openFileParallelAccessHDFWrapper(gfkbasename+'.'+cjs,fidgfk,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       if (.not. checkIdentityGreenFKSpectra(gfk,fidgfk,errmsg1)) then
          call add(errmsg1,2,'GFK data file inconsistent with meta file',myname)
          goto 11
       endif
       call readMetaGreenFKSpectra(gfk,fidgfk,'dataSpecificIntegers',errmsg1)
       if (.level.errmsg1 == 2) goto 11
   !
   !  create a group for event in HDF output and write event info
   !
       call createEventSynseisHDF(event,fidseis,evgrid,errmsg1)
       if (.level.errmsg1 == 2) goto 11
   !
   !  add df, sigma and length of spectrum (nf1,nf2) as attribute to file
   !
       if (compute_spectra) then
          d = [real(gfk%df),real(gfk%sigma)]; call arra%assoc1d(d)
          call writeArrayAttributeHDFWrapper(fidseis,"df_sigma",arra,errmsg1)
          if (.level.errmsg1 == 2) goto 11
          call arra%deassoc()
          id = [gfk%nf1,gfk%nf2]; call aria%assoc1d(id)
          call writeArrayAttributeHDFWrapper(fidseis,"nf1_nf2",aria,errmsg1)
          if (.level.errmsg1 == 2) goto 11
          call aria%deassoc()
       endif
   !
   !  add source node radius closest to true source radius as attribute to file
   !
       d = [real(rs)]; call arra%assoc1d(d)
       call writeArrayAttributeHDFWrapper(fidseis,"closestSourceNodeUsed",arra,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       call arra%deassoc(); deallocate(d)
   !
   !  precompute Legendre functions at unique epicentral distances on vertical boundaries
   !
       call findUniqueDistancesBoundaryPoints(bp,'v',thetas,phis,myrank,uniquedis,phi,point2dis,delmin,delmax)
       call precomputeLegendreFunctions(ierr)
       if(ierr == 2) goto 1
   !
       if (autoTimeshift) then
           !
           !  minimum epicentral distance from box to event
           !  get travel time and compute tred
           !  write as attribute to HDF if compute_spectra
           !
           call getReceiverTravelTimes(raytable,rnod(1),delmin,ttmin,errmsg1)
           if (.level.errmsg1 == 2) goto 11
           tred = ttmin-tsbuf
       else
           tred = .tshift.event
           ttmin = tred + tsbuf
       endif
   !
   !  maximum epicentral distance from box to event
   !  get travel time and compute tend = ttmax+twinlen
   !  to avoid scattering out of time window
   !  write as attribute to HDF
   !
       call getReceiverTravelTimes(raytable,rnod(size(rnod)),delmax,ttmax,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       tend = ttmax+twinlen
       if (tend-tred > tlen) then
          tend = tred+tlen
          twinlen = tend-ttmax
       end if
       d = [real(twinlen),real(tsbuf),real(ttmin),real(tred),real(tend),real(tscatlen)]; call arra%assoc1d(d)
       call writeArrayAttributeHDFWrapper(fidseis,"timing",arra,errmsg1)
       if (.level.errmsg1 == 2) goto 11
   !
   !  write timing to separate text file
   !
       if (myrank == 0) then
          open(unit=1,file=trim(timing_file),form='FORMATTED',status='UNKNOWN',action='WRITE')
          write(1,'(6f12.3)') twinlen,tsbuf,ttmin,tred,tend,tscatlen
          close(1)
       end if
   !
       if (myrank == 0) then
          write(6,'(a,2f15.3)') 'Minimum distance and travel time: ',delmin/mc_deg2radd,ttmin
          write(6,'(a,2f15.3)') 'Maximum distance and travel time: ',delmax/mc_deg2radd,ttmax
       endif
   !
   !  time shift filter, only for seismograms
   !
       if (.not. compute_spectra) then
          allocate(tsfil(gfk%nf2))
          tsfil = 0.d0
          do jf = gfk%nf1,gfk%nf2
             tsfil(jf) = zexp(mc_cid*mc_two_pid*(jf-1)*gfk%df*tred)*exp(+gfk%sigma*tred)
          enddo
       endif
   !
   !  set size of seismogram datasets per bp_point
   !
       if (compute_spectra) then
          countvt = [2*NCOMP,1,gfk%nf2]                 ! counts to be written per boundary point
       else
          nstrue = int((tend-tred)/dt)+1
          if (nstrue > nsamp) nstrue = nsamp
          if (convolve_stf) then                        ! check that nstrue > nstf
             if (nstrue < nstf) then
                call add(errmsg1,2,'Synthetics should be longer than source wavelet',myname)
                goto 11
             endif
          endif
          countvt = [2*NCOMP,1,nstrue]                  ! account for left shift and cutoff at end
          allocate(seis(2*NCOMP,nstrue))                ! allocate space for seismograms with true length
       endif
   !
   !  add dt and true length of seismogram to event group
   !
       if (.not. compute_spectra) then
          d = [real(dt),(nstrue-1)*real(dt)]; call arra%assoc1d(d)
          call writeArrayAttributeHDFWrapper(fidseis,"samplingIntervalLength",arra,errmsg1)
          if (.level.errmsg1 == 2) goto 11
          call arra%deassoc()
       endif
   !
   !  full proc-specific datasets: veltrac(6,nbpproc,nsamp), all boundaries
   !
       allocate(dimsall(3),procid(0:nproc-1),dset(0:nproc-1),dsp(0:nproc-1))
       allocate(dsetre(0:nproc-1),dsetim(0:nproc-1))
       allocate(dsetcoor(0:nproc-1),dspcoor(0:nproc-1))
       do ip = 0,nproc-1
          write(cip,'(i6.6)') ip
          call h5gcreate_f(fidseis,'proc_'+cip,procid(ip),ierr)
      !
      !  data set for seismograms
      !
          if (.not. compute_spectra) then
             dimsall = [2*NCOMP,getNumPerProcBoundaryPoints(bp,ip),nstrue]
             call h5screate_simple_f(3,dimsall,dsp(ip),ierr)
             if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
             call h5dcreate_f(procid(ip),"velocityTraction",H5T_NATIVE_REAL,dsp(ip),dset(ip),ierr)
             if (ierr < 0) then; print *,'h5dcreate_f veltrac'; goto 1; endif
      !
      !  data sets for spectra
      !
          else
             dimsall = [2*NCOMP,getNumPerProcBoundaryPoints(bp,ip),gfk%nf2]
             call h5screate_simple_f(3,dimsall,dsp(ip),ierr)
             if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
             call h5dcreate_f(procid(ip),"velocityTractionReal",H5T_NATIVE_REAL,dsp(ip),dsetre(ip),ierr)
             if (ierr < 0) then; print *,'h5dcreate_f veltracre'; goto 1; endif
             call h5dcreate_f(procid(ip),"velocityTractionImag",H5T_NATIVE_REAL,dsp(ip),dsetim(ip),ierr)
             if (ierr < 0) then; print *,'h5dcreate_f veltracim'; goto 1; endif
          endif
      !
      !  data set for coordinates
      !
          dimsall = [8,getNumPerProcBoundaryPoints(bp,ip),1]
          call h5screate_simple_f(2,dimsall(1:2),dspcoor(ip),ierr)
          if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
          call h5dcreate_f(procid(ip),"coordinates",H5T_NATIVE_REAL,dspcoor(ip),dsetcoor(ip),ierr)
          if (ierr < 0) then; print *,'h5dcreate_f coor'; goto 1; endif
      !
      !  nbpvert(ip) and nbpbot(ip) as attributes to proc-group
      !
          id = [ getSidesNumPerProcBoundaryPoints(bp,ip), getBottomNumPerProcBoundaryPoints(bp,ip)]
          call aria%assoc1d(id)
          call writeArrayAttributeHDFWrapper(procid(ip),"numBoundaryPointsPerProc",aria,errmsg1)
          if (.level.errmsg1 == 2) goto 11
          call aria%deassoc(); deallocate(id)
       enddo
   !-------------------------------------------------------------------------------------------------------
   !  Loop over points on vertical boundaries
   !
       jeprev = 0
       do while (iterateVerticalBoundaryPoints(bp,cbp,myrank,numtasks))
          ibp = getCounterSingleBoundaryPoint(cbp)
          idis = point2dis(ibp)                        ! distance index of bp
          call new(errmsg2,myname)
      !
      !  read entire Green spectra for given receiver node if new
      !
          re = cbp%r
          je = getNodeIndexFromRadiusGreenFKSpectra(gfk,re)
          if (je /= jeprev) then
             if (associated(gfkspr)) deallocate(gfkspr)
             if (associated(gfkspi)) deallocate(gfkspi)
             if (dabs(rnod(je)-re) > bp%rtol) then
                call add(errmsg2,2,'gfk nodes and BP radii seem not to fit',myname)
                print *,'True radius: ',re,' Node radius: ',rnod(je)
                goto 12
             endif
             re = rnod(je)
             call readDataNodeGreenFKSpectra(gfk,fidgfk,je,gfkspr,gfkspi,errmsg2)
             if (.level.errmsg2 == 2) goto 12
             if (myrank == 0) write(6,'(a,i5,f12.3)') 'Radial level: ',je,re
          endif
          jeprev = je
      !
      !  compute synthetic seismograms, allocates zsp, strainsp and tensionsp
      !
          call computeSynthetics(re,je,alf(idis),uniquedis(idis),phi(ibp),errmsg2)
          if (.level.errmsg2 == 2) goto 12
      !
      !  write hyperslab to process group (all components and time samples for given BP)
      !
          call getProcOffsetSidesBoundaryPoints(bp,ibp,off,iproc)
          offset = [0,off,0]
          call writeHyperslab(offset,iproc,uniquedis(idis),phi(ibp),errmsg2)
          if (.level.errmsg2 == 2) goto 12
      !
      !  clean up
      !
          deallocate(zsp,tensionsp,strainsp)
          call dealloc(errmsg2)
       enddo                                         ! end vertical boundary points loop
       do idis = 1,nd
          call dealloc(alf(idis))
       enddo
       deallocate(alf,uniquedis,phi,point2dis)
   !--------------------------------------------------------------------------------------------
   !  End of vertical boundaries
   !--------------------------------------------------------------------------------------------
       call barrier(mpisup)
   !--------------------------------------------------------------------------------------------
   !  Deal with bottom boundary
   !--------------------------------------------------------------------------------------------
   !  read GFK spectra for bottom radius
   !
       if (associated(gfkspr)) deallocate(gfkspr)
       if (associated(gfkspi)) deallocate(gfkspi)
       call readDataNodeGreenFKSpectra(gfk,fidgfk,1,gfkspr,gfkspi,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       if (myrank == 0) write(6,'(a,f12.3)') 'Radial level: ',rnod(1)
   !
   !  Loop over points on bottom boundary
   !
       do while (iterateBottomBoundaryPoints(bp,cbp,myrank,numtasks))
          call new(errmsg2,myname)
          ibp = getCounterSingleBoundaryPoint(cbp)
          call getEpicentralCoordinatesSingleBoundaryPoint(cbp,thetas,phis,delta,xi)
          call computeLegendreFunctions(alfb,gfk%nwnmax-1,2,delta,errmsg2)
      !
      !  compute synthetic seismograms, allocates zsp, strainsp and tensionsp
      !
          call computeSynthetics(rnod(1),1,alfb,delta,xi,errmsg2)
          if (.level.errmsg2 == 2) goto 12
      !
      !  write hyperslab (all components and time samples for given BP)
      !
          call getProcOffsetBottomBoundaryPoints(bp,ibp,off,iproc)
          offset = [0,off,0]
          call writeHyperslab(offset,iproc,delta,xi,errmsg2)
          if (.level.errmsg2 == 2) goto 12
          deallocate(zsp,tensionsp,strainsp)
          call dealloc(errmsg2)          
          call dealloc(alfb)
          if (mod(ibp,1000) == 0) print *,'Point: ',ibp,myrank,cbp%r,cbp%lat,cbp%lon
       enddo
   !--------------------------------------------------------------------------------------------
   !  End of bottom boundary
   !--------------------------------------------------------------------------------------------
       call barrier(mpisup)
   !
   !  clean up stuff allocated in event loop
   !      
       do ip = 0,nproc-1
          call h5sclose_f(dsp(ip),ierr)
          if (.not. compute_spectra) then
             call h5dclose_f(dset(ip),ierr)
          else
             call h5dclose_f(dsetre(ip),ierr)
             call h5dclose_f(dsetim(ip),ierr)
          endif
          call h5sclose_f(dspcoor(ip),ierr)
          call h5dclose_f(dsetcoor(ip),ierr)
          call h5gclose_f(procid(ip),ierr)
       enddo
       deallocate(dimsall,procid,dset,dsetre,dsetim,dsp,dsetcoor,dspcoor)
       call h5gclose_f(evgrid,ierr)
       call h5fclose_f(fidseis,ierr)
   !
       deallocate(gfkspr,gfkspi)
       deallocate(fomt)
       if (allocated(tsfil)) deallocate(tsfil)
       if (allocated(seis)) deallocate(seis)
       if (allocated(stf)) deallocate(stf)

       call h5fclose_f(fidgfk,ierr)
       call h5pclose_f(xferprp,ierr)                                                        ! close xfer property list
       call dealloc(errmsg1)
    enddo                                       ! end event loop
!-----------------------------------------------------------------------------------
!  deallocation
!
    call h5close_f(ierr)
    call dealloc(event_list)
    call dealloc(errmsg)
    call dealloc(mpisup)
    call dealloc(gfk)
    call dealloc(raytable)
    if (allocated(urs)) deallocate(urs)
    if (allocated(uspre)) deallocate(uspre,uspim)
    if (allocated(evfilter)) deallocate(evfilter)
!
!  treat error conditions    
!
1   if (ierr /= 0) then
       print *,'ERROR: ',trim(errstr)
       call abort(mpisup)
    endif
10  if (.level.errmsg == 2) then
       call print(errmsg)
       call abort(mpisup)
    endif
11  if (.level.errmsg1 == 2) then
       call print(errmsg1)
       call abort(mpisup)
    endif
12  if (.level.errmsg2 == 2) then
       call print(errmsg2)
       call abort(mpisup)
    endif
!----------------------------------------------------------------------
    contains
!----------------------------------------------------------------------
!  compute normal stress in Specfem box centered cartesian coordinates
!  get elastic constants at BP radius 
!  calculate stress tensor from strain tensor (spherical)          
!  rotate from epicentral ZRT into SPECFEM box centered cartesian components
!  mulitply with normal vector to get normal stresses
!  gfk:      Green FK object
!  jf:       frequency index
!  je:       radial node index
!  epidis:   epicentral distance in rad
!  phi:      epicentral azimuth (counted from south over east) in rad
!  thetas:   geographical colatitude of source (rad)
!  phis:     geographical longitude of source (rad)
!  thetac:   geographical colatitude of box center (rad)
!  phic:     geographical longitude of box center (rad)
!  nx,ny,nz: normal vector (box components)
!  zsp:      displacements (epicentral RLT) on input, box centered cartesian on output (XYZ)
!  strainsp: spherical strain tensor components in epicentral system
!  tension:  (out) traction vector in box centered cartesian
!
    subroutine computeTension(gfk,jf,je,epidis,phi,thetas,phis,thetac,phic,nx,ny,nz,zsp,strainsp,tension)
    type (green_fk_spectra) :: gfk
    integer :: jf,je
    double precision :: epidis,phi,thetas,phis,thetac,phic,nx,ny,nz
    double complex, dimension(:) :: zsp,tension
    double complex, dimension(:,:) :: strainsp
    double complex, dimension(NELCON) :: elcon
    double complex, dimension(NCOMP) :: uc,ug,ub
    double complex, dimension(NCOMP,NCOMP) :: tc,tg,tb,stress
!
    call getComplexConstantsSelectedExternalRadialNodes(gfk%exnod,gfk%fref,&
         (jf-1)*gfk%df,gfk%attnmode,je,elcon)
    call sphericalStressSingleFrequencyDisplacement(strainsp,elcon,stress)
    
    call vectorLCFromLSAxesRotation(epidis,phi,zsp,uc)          ! epicentral cartesian comps
    call vectorGCFromLCAxesRotation(thetas,phis,uc,ug)          ! global cartesian comps
    call vectorLCfromGCAxesRotation(thetac,phic,ug,ub)          ! box centered cartesian comps
    call vectorRCfromLCAxesRotation(0.5d0*mc_pid,ub,zsp)        ! box centered cartesian with x-axis pointing east and y north
    
    call tensorLCFromLSAxesRotation(epidis,phi,stress,tc)       ! epicentral cartesian comps
    call tensorGCFromLCAxesRotation(thetas,phis,tc,tg)          ! global cartesian comps
    call tensorLCfromGCAxesRotation(thetac,phic,tg,tb)          ! box centered cartesian comps
    call tensorRCfromLCAxesRotation(0.5d0*mc_pid,tb,stress)     ! box centered cartesian with x-axis pointing east and y north
    tension = matmul(stress,[nx,ny,nz])
    end subroutine computeTension
!-------------------------------------------------------------------------------------------------
!  Precompute Legendre functions at unique epicentral distances
!
    subroutine precomputeLegendreFunctions(ierr)
    integer :: ierr
    nd = size(uniquedis)
    allocate(alf(nd))
    do idis = 1,nd
       call new(errmsg2,myname)
       call computeLegendreFunctions(alf(idis),gfk%nwnmax-1,2,uniquedis(idis),errmsg2)
       if (.level.errmsg2 == 2) then
          ierr = 2; return
       endif
       call dealloc(errmsg2)
    enddo
    end subroutine precomputeLegendreFunctions
!-----------------------------------------------------------------------------------------------------
!  compute synthetics
!
    subroutine computeSynthetics(re,je,alf,delta,xi,errmsg)
    type (legendre_functions) :: alf
    double precision :: re,delta,xi
    type (error_message) :: errmsg
    integer :: je
    double complex :: zom
!
!  compute RLT displacement spectra, allocates zsp and strainsp
!
    call sphericalHarmonicsDisplacementSpectra(gfk,gfkspr,gfkspi,fomt,&
         alf,re,delta,xi,zsp,errmsg,strainsp)
    if (.level.errmsg == 2) return
!
!  get normal stresses, transforms everything to XYZ specfem coordinates
!  differentiate displacements to get velocity (needed for injection)
!
    allocate(tensionsp(gfk%nf2,NCOMP))
    tensionsp = 0.d0
    do jf = gfk%nf1,gfk%nf2
       zom = dcmplx(mc_two_pid*(jf-1)*gfk%df,-gfk%sigma)
       call computeTension(gfk,jf,je,delta,xi,thetas,phis,&
            thetac,phic,cbp%nx,cbp%ny,cbp%nz,zsp(jf,:),strainsp(jf,:,:),tensionsp(jf,:))
       if (.not. compute_spectra) then
          zsp(jf,:) = zsp(jf,:)*evfilter(jf)*mc_cid*zom*tsfil(jf)      ! multiply with event filter and i(om-i*sig)
          tensionsp(jf,:) = tensionsp(jf,:)*evfilter(jf)*tsfil(jf)     ! multiply with event filter and time shift
       endif
    enddo
!  
!  transform to time domain
!  convolve with STF (moment rate function)
!  explanation:
!   Gemini computes displacement spectra for a delta source wavelet (Green displacement).
!   Since the source wavelet represents moment rate (synthetic displacement fitted to observed velocity)
!   a convolution with the moment rate yields ground velocity.
!   
    if (.not. compute_spectra) then
       do ic = 1,NCOMP
          call transformFrequencyTime(gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,zsp(:,ic),nsamp,dt,urs,errmsg)
          if (.level.errmsg == 2) return
          if (convolve_stf) then
             call convolveWithSourceWavelet(nstrue,nstf,tanfstf,dt,stf,urs)
          endif
          if (.level.errmsg == 2) return
          seis(ic,:) = urs(1:nstrue)
       !
          call transformFrequencyTime(gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,tensionsp(:,ic),nsamp,dt,urs,errmsg)
          if (.level.errmsg == 2) return
          if (convolve_stf) then
             call convolveWithSourceWavelet(nstrue,nstf,tanfstf,dt,stf,urs)
          endif
          if (.level.errmsg == 2) return
          seis(ic+NCOMP,:) = urs(1:nstrue)
       enddo
    else
       do ic = 1,NCOMP
          uspre(ic,:) = real(zsp(:,ic))
          uspim(ic,:) = aimag(zsp(:,ic))
          uspre(ic+NCOMP,:) = real(tensionsp(:,ic))
          uspim(ic+NCOMP,:) = aimag(tensionsp(:,ic))
       enddo
    endif
    end subroutine computeSynthetics
!-----------------------------------------------------------------------------------------------------
!  write synthetics to hyperslab
!  returns pseudo mesh coordinates of boundary point
!
    subroutine writeHyperslab(offset,iproc,delta,xi,errmsg)
    integer(hsize_t), dimension(:) :: offset
    integer :: iproc
    double precision :: delta,xi
    type (error_message) :: errmsg
    double precision :: xg,yg,zg,xr,yr,zr,x,y,z
    real, dimension(8) :: d
!
!  write seismograms
!
    if (.not. compute_spectra) then
       call arra%assoc2d(seis)
       call writeArrayHDFWrapper(procid(iproc),'velocityTraction',arra,errmsg,&
            xferprp = xferprp,ds = dset(iproc),offset = offset,count = countvt)
       if (.level.errmsg == 2) return
       call arra%deassoc()
!
!  write real and imaginary part of spectra
!
    else
       call arra%assoc2d(uspre)
       call writeArrayHDFWrapper(procid(iproc),'velocityTractionReal',arra,errmsg,&
            xferprp = xferprp,ds = dsetre(iproc),offset = offset,count = countvt)
       if (.level.errmsg == 2) return
       call arra%deassoc()
       call arra%assoc2d(uspim)
       call writeArrayHDFWrapper(procid(iproc),'velocityTractionImag',arra,errmsg,&
            xferprp = xferprp,ds = dsetim(iproc),offset = offset,count = countvt)
       if (.level.errmsg == 2) return
       call arra%deassoc()
    endif

!  compute local cartesian coordinates and write spherical, cartesian and epicentral
!  coordinates to HDF, using the same ordering as for the seismograms
!
    call coordinatesLCfromLSAxesRotation(cbp%r,0.5d0*mc_pid-cbp%lat*mc_deg2radd,cbp%lon*mc_deg2radd,xg,yg,zg)  ! from geographical to global cartesian
    call coordinatesLCfromGCAxesRotation(thetac,phic,xg,yg,zg,xr,yr,zr)      ! from global cartesian to box cartesian
    call coordinatesRCfromLCAxesRotation(+0.5*mc_pid,xr,yr,zr,x,y,z)         ! from box cartesian to rotated box cartesian (x-axis east)
    d = [x,y,z-gfk%rearth,cbp%r,cbp%lat,cbp%lon,delta/mc_deg2radd,xi/mc_deg2radd]
    call arra%assoc1d(d)
    call writeArrayHDFWrapper(procid(iproc),'coordinates',arra,errmsg,&
         xferprp = xferprp,ds = dsetcoor(iproc),offset = offset(1:2),count = [8_hsize_t,1_hsize_t])
    if (.level.errmsg == 2) return
    call arra%deassoc()
    end subroutine writeHyperslab
!--------------------------------------------------------------------------
!  convolve synthetics with source wavelet (see noncausfalt.txt)
!
    subroutine convolveWithSourceWavelet(nd,nstf,tanf,dt,stf,urs)
    integer :: nstf,nd
    double precision :: tanf,dt
    real, dimension(:) :: stf,urs
    real, dimension(:), allocatable :: f
    integer :: nhm,k,n,j
!
    nhm = nint(abs(tanf)/dt)                                     ! Number of samples of h_j with j < 1
    allocate(f(nd))
    do k = 1,nd                                                  ! Loop over samples of synthetics
       f(k) = 0.d0
       do n = max(1,k-nd+nhm+1),min(k+nhm,nstf)                  ! Loop over samples of stf
          j = n-nhm
          f(k) = f(k)+dt*urs(k-j+1)*stf(n)
       enddo
    enddo
    urs(1:nd) = f
    deallocate(f)
    end subroutine convolveWithSourceWavelet
!--------------------------------------------------------------------------
!  apply a cos**2 end taper to phase window and zero beyond
!  n:      total number of samples
!  nta:    first sample of taper
!  ntl:    length of taper in samples
!  u:      1D-array for seismogram
!
    subroutine endTaperPhaseWindow(n,nta,ntl,u,errmsg)
    integer :: n,nta,ntl
    real, dimension(:) :: u
    type (error_message) :: errmsg
    integer :: i
    real :: taparg
!
    if (nta+ntl > n) then
       call add(errmsg,2,'taper exceeds length of time series','taperPhaseWindow')
       return
    end if
    do i = nta+1,nta+ntl
       taparg = real(i-nta)/real(ntl)*mc_pid/2.d0
       u(i) = u(i)*cos(taparg)**2
    end do
    u(nta+ntl+1:n) = 0.0
    end subroutine endTaperPhaseWindow
!----------------------------------------------------------------------------
!   apply a front taper inside the buffer time range
!   use 2/3 of buffer length as taper length
!   start taper at nzero
!   zero before nzero
!
    subroutine frontTaperPhaseWindow(nzero,tsbuf,dt,u)
    integer :: nzero
    double precision :: tsbuf,dt
    real, dimension(:) :: u
    real :: taparg
    integer :: i,ntl
!
    ntl = 0.6666*nint(tsbuf/dt)
    do i = nzero+1,nzero+ntl
       taparg = real(ntl+nzero-i)/real(ntl)*mc_pid/2.d0
       u(i) = u(i)*cos(taparg)**2
    end do
    u(1:nzero) = 0.0
    end subroutine frontTaperPhaseWindow
!-----------------------------------------------------------------------------------------------------
end program
