! ===============================================================================
!  Main program for the computation of synthetic seismograms
! ===============================================================================
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
!  Compute synthetic seismograms from Green FK spectra
!  User provides a station file, event file, filter
!  files for each event and filter files for each station.
!  Seismograms are written in ASCII into separate files for
!  each station and component.
!  Assume case of ONE receiver depth and varying source depths.
!  Typically, receiver will be at the surface.
!  Expects presence of a filter file for each event with file name filter_evid.
!  Expects presence of a filter file for each station and component 
!  with file name filter_staname_comp.
!-----------------------------------------------------------------------------
program computeSyntheticSeismograms
       use hdf5
    use mathConstants
    use argumentParser
    use greenFKSpectra
    use legendreFunctions
    use besselFunctions
    use displacementSpectra          ! NCOMP defined here
    use axesRotation
    use readEventStationFile
    use eventStationGeometry
    use eventFilter
    use frequencyTime
    use errorMessage
    use inputParameter
    use string
    use synseisHDF
    use travelTimes
    use asciiSynseisIO
    implicit none
    type (argument_parser) :: ap
    type (green_fk_spectra) :: gfk
    type (error_message) :: errmsg,errmsg1,errmsg2
    type (input_parameter) :: inpar
    type (seismic_event) :: event
    type (seismic_station) :: station
    type (seismic_event_list) :: event_list
    type (seismic_network) :: station_list
    type (legendre_functions) :: alf
    type (any_rank_real_array) :: ara
    type (raytable_travel_times) :: raytable
    real, dimension(:,:), pointer :: gfkspr,gfkspi
    real, dimension(:,:), allocatable :: urs
    real, dimension(:), allocatable :: stf
    double precision, dimension(:), allocatable :: fomt
    real, dimension(:,:,:), pointer :: p3
    double precision, dimension(:), pointer :: evfilspecs
    double precision :: dt,dtp,zs,rs,epidis,propdir,phi,tanf,tanfstf,tred,tsbuf,twinlen,tt,tend,tlen
    real :: dtstf
    double precision, dimension(:,:), allocatable :: besselj
    double complex, dimension(:,:,:), allocatable :: stafilter
    double complex, dimension(:), allocatable :: evfilter, tsfil
    double complex, dimension(:,:), allocatable :: zsp
    double complex :: uz,un,ue,zom
    logical :: autoTimeshift,split,convolve_stf,single_event,write_ascii
    integer(hid_t) :: fid,filid,fidseis,evgrid
    integer :: nsamp,nstat,is,ic,jr,jf,ierr,ios,nstrue,nstf
    character(len=max_length_string) :: eventfile,stationfile,evfiltype,&
          gfkfile,parfile,file_station_filter,synseisfile,seismoType,stffile,&
          source_type,raytablefile,dsvbasename,stfpath,evsel,dummy1,dummy2,asciifile
    character(len=max_length_string), dimension(:), pointer :: word
    character (len=3) :: comp = 'ZNE'
    character(len=max_length_string) :: errstr
    character (len=80), dimension(2) :: stype_keywords
    character (len=27) :: myname = 'computeSyntheticSeismograms'
    character (len=80), dimension(13) :: par_keys
!
!  keywords for input parameters
!
    data par_keys/'DSVBASENAME','FILE_SYNSEIS','FILE_EVENT_LIST','FILE_STATION_LIST',&
         'FILE_STATION_FILTER','EVENT_FILTER_TYPE','EVENT_FILTER_SPECS','SEISMO_TYPE',&
         'AUTOMATIC_TIMESHIFT','TIME_SHIFT_BUFFER','PHASE_WINDOW_LENGTH','FILE_RAY_TABLE',&
         'SOURCE_TYPE'/
    data stype_keywords/'velocity','displacement'/
!---------------------------------------------------------------------------------------------
    call new(errmsg,myname)
!
    call init(ap,myname,'Compute synthetic seismograms using GreenFKSpectraForSynthetics. '//new_line('a')//&
              'Work through event list and write seismograms into one large HDF file')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call addOption(ap,'-split',.false.,'Write seismograms into separate event subfolders. '//new_line('a')//&
                   'Take the dirname part set in the Gemini parfile as the parent of the event subfolders. '//new_line('a')//&
                   'Create event subfolders if not present. '//&
                   'Use the non-dirname part as filename','sval')
    call addOption(ap,'-stf',.true.,'Path to parent of event subfolders. '//new_line('a')//&
                   'Convolve with source time function in subfolders. '//new_line('a')//&
                   'Assume name of file is stf.txt','sval','None')
    call addOption(ap,'-single',.true.,'Run only for specified event in event list','sval','None')
    call addOption(ap,'-ascii',.true.,"Base name of ASCII seismograms. "//new_line('a')//&
                   "If -split, write to output event subfolder using this basename. "//new_line('a')//&
                   "If -single write to run folder using this basename",'sval','gemini_synthetics.asc')
    call parse(ap)
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    parfile = ap.sval.'parfile'
    split = ap.optset.'-split'
    stfpath = ap.sval.'-stf'
    convolve_stf = ap.optset.'-stf'
    evsel = ap.sval.'-single'
    single_event = ap.optset.'-single'
    write_ascii = ap.optset.'-ascii'
    asciifile = ap.sval.'-ascii'
    if (split .and. single_event) then
       call add(errmsg,2,'Either use -split or -single but not both',myname)
       goto 10
    endif
    if ((write_ascii .and. .not. split) .and. (write_ascii .and. .not. single_event)) then
       call add(errmsg,2,'Use of -ascii option only together with -split or -single',myname)
       goto 10
    endif
    call document(ap); call dealloc(ap)
!-------------------------------------------------------------
!  read input parameters from parfile
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 10
    source_type = inpar.sval.'SOURCE_TYPE'
    dsvbasename = inpar.sval.'DSVBASENAME'
    gfkfile = dsvbasename+"."+source_type
    synseisfile = inpar.sval.'FILE_SYNSEIS'
    eventfile = inpar.sval.'FILE_EVENT_LIST'
    stationfile = inpar.sval.'FILE_STATION_LIST'
    file_station_filter = inpar.sval.'FILE_STATION_FILTER'
    evfiltype = inpar.sval.'EVENT_FILTER_TYPE'
    evfilspecs => dvecp(inpar,'EVENT_FILTER_SPECS',4,ios)
    seismoType = inpar.sval.'SEISMO_TYPE'
    if (ios /= 0) then
       call add(errmsg,2,'Input for event filter specs cannot be read',myname)
       goto 10
    endif
    autoTimeshift = inpar.lval.'AUTOMATIC_TIMESHIFT'
    if (autoTimeshift) then
       raytablefile = inpar.sval.'FILE_RAY_TABLE'
       tsbuf = inpar.dval.'TIME_SHIFT_BUFFER'
       twinlen = inpar.dval.'PHASE_WINDOW_LENGTH'
    endif
    call printInputParameter(inpar)
    call dealloc(inpar)
!-----------------------------------------------------------------------------
!  check validity of seismo types
!
    if (.not. (seismoType.equal.stype_keywords(1)) .and. .not.(seismoType.equal.stype_keywords(2))) then
       call add(errmsg,2,'Invalid seismo type',myname)
       call print(errmsg); stop
    endif
!----------------------------------------------------------------
!  read event and station file
!  use ASKI convention for event and station file
!
    call createEventListFromEventFile(eventfile,1,'ASKI_events',event_list,errmsg)
    if (.level.errmsg == 2) goto 10
    call createStationListFromStationFile(stationfile,1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) goto 10
    nstat = .nstat.station_list
!
!  check whether csys in station and event file are consistent
!
    if (.csys.event_list /= .csys.station_list) then
       call add(errmsg,2,'CSYS inconsistent in event and station file',myname)
       goto 10
    endif
!------------------------------------------------------------------
!  open frequency wavenumber spectra file
!  having all necessary source depths    
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
    call openFileRoHDFWrapper(gfkfile,fid,errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'identity',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'reals',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'integers',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'numberOfWavenumbersForFrequency',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'externalNodes',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'dataSpecificIntegers',errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read ray table file
!
    if (autoTimeshift) then
       call readTableTravelTimes(raytable,raytablefile,errmsg)
       if (.level.errmsg == 2) goto 10
    endif
!
!  allocate space
!
    allocate(stafilter(gfk%nf2,NCOMP,nstat))
!----------------------------------------------------------------
!  Loop over stations to read station filters
!
    call openFileRoHDFWrapper(file_station_filter,filid,errmsg)
    if (.level.errmsg == 2) goto 10
    do while (nextStationSeismicNetwork(station_list,station,is))
       call readArrayHDFWrapper(filid,trim((.netcode.station)+'.'+(.staname.station)),ara,errmsg)
       if (.level.errmsg == 2) goto 10
       p3 => ara%get3d()
       stafilter(gfk%nf1:gfk%nf2,1:NCOMP,is) = dcmplx(p3(1,:,1:NCOMP),p3(2,:,1:NCOMP))
       deallocate(p3)
    enddo                    ! next station
    call h5fclose_f(filid,ierr)
    if (ierr < 0) then; errstr = 'h5fclose'; goto 1; endif
!
!  Calculate event filter (common for all events, comprises bandpass filtering)
!
    call computeEventFilter(evfiltype,evfilspecs,gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,evfilter,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  properties of time series and allocate
!
    tlen = 1.d0/gfk%df
    dtp = 0.25d0/(2.d0*gfk%nf2*gfk%df)
    call newNsampDtFrequencyTime(gfk%df,dtp,nsamp,dt)
    allocate(urs(nsamp,NCOMP))
!
    if (.not. split) call createFileHDFWrapper(synseisfile,fidseis,errmsg)    ! open one synseis file for all events
    if (.level.errmsg == 2) goto 10
!---------------------------------------------------------------------------
!  Loop over events
!---------------------------------------------------------------------------
    do while (nextEventSeismicEventList(event_list,event))
       if (single_event) then                                             ! selected event desired
          if (.not. ((.evid.event).equal.evsel)) cycle
       endif
       call new(errmsg1,myname)
       call printSeismicEvent(event)
       if ((.styp.event) .ne. gfk%istyp) then
          call add(errmsg1,2,'GFk source type inconsistent with event source type',myname)
          goto 11
       endif
       if (split) then                                                    ! open synseis file event wise
          word => getWordsString(synseisfile,'/',errmsg1)                 ! write into subfolder named according to evid
          if (.level.errmsg1 == 2) goto 11                                ! create folder if not there
          call system('mkdir -p '//trim(word(1))//'/'//trim(.evid.event),ierr)
          if (ierr /= 0) then; errstr = 'cannot create event subfolder'; goto 1; endif
          call createFileHDFWrapper(trim(word(1))//'/'//trim(.evid.event)//'/'//trim(word(2)),fidseis,errmsg1)
          if (.level.errmsg1 == 2) goto 11
          if (write_ascii) then
             dummy1 = trim(word(1))//'/'//trim(.evid.event)//'/'//trim(asciifile)
             open(7, file = trim(dummy1)//'.'//comp(1:1), status = 'unknown',iostat = ierr)
             open(8, file = trim(dummy1)//'.'//comp(2:2), status = 'unknown',iostat = ierr)
             open(9, file = trim(dummy1)//'.'//comp(3:3), status = 'unknown',iostat = ierr)
             if (ierr /= 0) then; errstr = 'cannot open ASCII seismogram file'; goto 1; endif
          endif
          deallocate(word)
       endif
       if (single_event .and. write_ascii) then
             open(7, file = trim(asciifile)//'.'//comp(1:1), status = 'unknown',iostat = ierr)
             open(8, file = trim(asciifile)//'.'//comp(2:2), status = 'unknown',iostat = ierr)
             open(9, file = trim(asciifile)//'.'//comp(3:3), status = 'unknown',iostat = ierr)
             if (ierr /= 0) then; errstr = 'cannot open ASCII seismogram file'; goto 1; endif
       endif
       call createEventSynseisHDF(event,fidseis,evgrid,errmsg1)
       if (.level.errmsg1 == 2) goto 11
       zs = .sdepth.event                                           ! in m for spherical system and m for cartesian one
       tanf = .tanfdp.(.otime.event)                                ! origin time of event in seconds after midnight
       call getForceMomentSeismicEvent(event,fomt)                  ! force or moment components depending on istyp
       jr = getNodeIndexFromRadiusGreenFKSpectra(gfk,gfk%rearth-zs) ! node index of source depth, assume meter for radii
       rs = getDoubleRadiusSelectedExternalRadialNodes(gfk%exnod,jr)
       write(6,'(a,f12.3,a,f12.3)') 'Desired source radius: ',gfk%rearth-zs,' True radius: ',rs
    !
    !  read entire Green spectra for given source node
    !
       call readDataNodeGreenFKSpectra(gfk,fid,jr,gfkspr,gfkspi,errmsg1)
       if (.level.errmsg1 == 2) goto 11
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
             call add(errmsg,2,'expect zero or negative start time for source wavelet',myname)
             goto 11
          endif
       endif
    !
    !  Loop over stations
    !
       do while (nextStationSeismicNetwork(station_list,station,is))
          call new(errmsg2,myname)
          if (.level.errmsg2 == 2) goto 12
       !
       !  get epicentral distances and and receiver angles (from south over east or from x to y)
       !  as seen at the event location
       !  get propagation directions as seen at receiver (from S to E or from x to y)
       !
          select case (.csys.event)
          case ('S')
             call deltaMeterPhiRadSphericalEventStationGeometry(event,station,gfk%rearth,epidis,phi)
             call propdirRadSphericalEventStationGeometry(event,station,propdir)
             print *,trim(.staname.station),epidis,phi*(180.d0/mc_pid)
          case ('C')
             call deltaMeterPhiRadHalfspaceEventStationGeometry(event,station,epidis,phi)
             propdir = phi
          case default
             call add(errmsg2,2,'Unknown coordinate system',myname)
             if (.level.errmsg2 == 2) goto 12
          end select
       !
       !  compute Bessel functions or Legendre functions and RLT displacement spectra
       !
          if (gfk%global == 0) then
             call computeBesselFunctions(gfk%nwnmax,gfk%dwn,epidis,besselj)
             call besselDisplacementSpectra(gfk,gfkspr,gfkspi,fomt,besselj,gfk%rse,epidis/gfk%rearth,phi,zsp,errmsg2)
             if (.level.errmsg2 == 2) goto 12
          else
             call computeLegendreFunctions(alf,gfk%nwnmax-1,2,epidis/gfk%rearth,errmsg2)
             call sphericalHarmonicsDisplacementSpectra(gfk,gfkspr,gfkspi,fomt,alf,gfk%rse,epidis/gfk%rearth,phi,zsp,errmsg2)
             if (.level.errmsg2 == 2) goto 12
          endif
       !
          if (autoTimeshift .and. gfk%global == 1) then
             call getTravelTimes(raytable,rs,gfk%rse,epidis/gfk%rearth,tt,errmsg2)
             if (.level.errmsg2 == 2) goto 12
             tred = tt-tsbuf
             tend = tt+twinlen
             if (tend-tred > tlen) tend = tred+tlen
             nstrue = int((tend-tred)/dt)+1
             if (nstrue > nsamp) nstrue = nsamp
             print *,'Travel time: ',tt,' Reduction time: ',tred,' End time: ',tend,' True number of samples: ',nstrue
          else
             tred = getTimeshiftSeismicEvent(event)
             nstrue = nsamp
             print *,' Reduction time: ',tred
          endif
       !
       !  time shift filter
       !
          allocate(tsfil(gfk%nf2))
          tsfil = 0.d0
          do jf = gfk%nf1,gfk%nf2
             tsfil(jf) = zexp(mc_cid*mc_two_pid*(jf-1)*gfk%df*tred)*exp(+gfk%sigma*tred)
          enddo
       !
       !  multiply with filters
       !
          do ic = 1,NCOMP
             zsp(:,ic) = zsp(:,ic)*evfilter*stafilter(:,ic,is)*tsfil
          enddo
       !  
       !  rotate into ZNE components and transform to time domain
       !   
          do jf = gfk%nf1,gfk%nf2
             if (seismoType == "velocity") then
                zom = dcmplx(mc_two_pid*(jf-1)*gfk%df,-gfk%sigma)
                zsp(jf,:) = zsp(jf,:)*mc_cid*zom
             else if (seismoType == "displacement") then
                continue
             endif
             call vectorZNEfromRLTAxesRotation(propdir,zsp(jf,1),zsp(jf,2),zsp(jf,3),uz,un,ue)
             zsp(jf,1) = uz; zsp(jf,2) = un; zsp(jf,3) = ue
          enddo
          do ic = 1,NCOMP
             call transformFrequencyTime(gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,zsp(:,ic),nsamp,dt,urs(:,ic),errmsg2)
             if (.level.errmsg == 2) goto 12
          enddo
       !
       !  convolve with source time function if desired
       !
          if (convolve_stf) then
             if (nstrue < nstf) then
                call add(errmsg2,2,'Synthetics should be longer than source wavelet',myname)
             endif
             call convolveWithSourceWavelet(nstrue,nstf,tanfstf,dt,stf,urs)
          endif
       !
       !  write to HDF (and ASCII if desired)
       !
          call writeStationSynseisHDF(station,fidseis,nstrue,tanf+tred,dt,'ZNE',urs(1:nstrue,:),errmsg2)
          if (.level.errmsg2 == 2) goto 12
          if (write_ascii) then
             do ic = 1,NCOMP
                call writeAsciiSynseisIO(6+ic,asciifile,nstrue,tanf+tred,real(dt),.staname.station,comp(ic:ic),&
                                         urs(1:nstrue,ic),ierr,op = .false.,cl = .false.)
                if (ierr < 0) then; errstr = "cannot write ASCII synseis file"; goto 1; endif
             enddo
          endif
          if (gfk%global == 1) call dealloc(alf)
          if (gfk%global == 0) deallocate(besselj)
          deallocate(zsp,tsfil)
          call dealloc(errmsg2)
       enddo                                                                         ! end station loop
       deallocate(gfkspr,gfkspi,fomt)
       if (allocated(stf)) deallocate(stf)
       call dealloc(errmsg1)
       call h5gclose_f(evgrid,ierr)
       if (ierr < 0) then; errstr = "error closing event group"; goto 1; endif
       if (split) then                                                                ! close synseis file if event specific
          call h5fclose_f(fidseis,ierr)
          if (ierr < 0) then; errstr = "error closing seis file"; goto 1; endif
          if (write_ascii) then; close(7); close(8); close(9); endif
       endif
    enddo                                                                             ! end event loop
!
!  deallocation
!
    call h5fclose_f(fid,ierr)
    if (ierr < 0) then; errstr = "error closing gfk file"; goto 1; endif
    if (.not. split) then
       call h5fclose_f(fidseis,ierr)
       if (ierr < 0) then; errstr = "error closing seis file"; goto 1; endif
    endif
    call h5close_f(ierr)
!
    deallocate(urs)
    deallocate(stafilter,evfilter)
    call dealloc(event_list); call dealloc(station_list)
    call dealloc(gfk); call dealloc(errmsg)
!
!  treat errors
!
 1  if (ierr < 0) then
       call add(errmsg,2,trim(errstr),myname)
       call print(errmsg)
       stop
    endif
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif
11  if (.level.errmsg1 == 2) then
       call print(errmsg1)
    endif
12  if (.level.errmsg2 == 2) then
       call print(errmsg2)
    endif
!--------------------------------------------------------------------------
contains
!--------------------------------------------------------------------------
!  convolve synthetics with source wavelet (see noncausfalt.txt)
!
    subroutine convolveWithSourceWavelet(nd,nstf,tanf,dt,stf,urs)
    integer :: nstf,nd
    double precision :: tanf,dt
    real, dimension(:) :: stf
    real, dimension(:,:) :: urs
    real, dimension(:), allocatable :: f
    integer :: nhm,k,n,j,ic
!
    nhm = nint(abs(tanf)/dt)                                     ! Number of samples of h_j with j < 1
    do ic = 1,NCOMP
       allocate(f(nd))
       do k = 1,nd                                               ! Loop over samples of synthetics
          f(k) = 0.d0
          do n = max(1,k-nd+nhm+1),min(k+nhm,nstf)               ! Loop over samples of stf
             j = n-nhm
             f(k) = f(k)+dt*urs(k-j+1,ic)*stf(n)
          enddo
       enddo
       urs(1:nd,ic) = f
       deallocate(f)
    enddo
    end subroutine convolveWithSourceWavelet
!
end program
