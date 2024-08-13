! ===============================================================================
!  Compute impulse response for event filter
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
!  Compute impules response of event filter
!-----------------------------------------------------------------------------
program computeSyntheticsForASKI
       use hdf5
    use mathConstants
    use argumentParser
    use greenFKSpectra
    use readEventStationFile
    use frequencyTime
    use eventFilter
    use errorMessage
    use inputParameter
    use fileUnitHandler
    use string
    use synseisHDF
    implicit none
    type (argument_parser) :: ap
    type (green_fk_spectra) :: gfk
    type (seismic_event) :: event
    type (seismic_station) :: station
    type (seismic_event_list) :: event_list
    type (seismic_network) :: station_list
    type (error_message) :: errmsg,errmsg1,errmsg2
    type (input_parameter) :: inpar
    real, dimension(:,:), allocatable :: urs
    double complex, dimension(:), allocatable :: evfilter,zsp
    double complex :: zom
    double precision, dimension(:), pointer :: evfilspecs
    double precision :: dt,dtp
    integer(hid_t) :: fidseis,evgrid,fidmeta
    integer :: nsamp,is,jf
    integer :: ierr
    character(len=max_length_string) :: eventfile,parfile,synseisfile,gfkmetafile,evfiltype,stationfile
    character(len=max_length_string) :: errstr
    character (len=24) :: myname = 'testButterworthImpulsRsponse'
    character (len=80), dimension(6) :: par_keys
!
!  keywords for input parameters
!
    data par_keys/'FILE_SYNTHETICS','FILE_EVENT_LIST','EVENT_FILTER_TYPE','EVENT_FILTER_SPECS',&
         'FILE_GFK_META','FILE_STATION_LIST'/
!--------------------------------------------------------------------------------------------------
    call init(ap,myname,'Compute Butterworth filter impulse reponse')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call parse(ap)
    parfile = ap.sval.'parfile'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
!-----------------------------------------------------------------------------
    call document(ap)
    call dealloc(ap)
!
    call new(errmsg,myname)
!-------------------------------------------------------------
!  read input parameters from parfile
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 1
    gfkmetafile = inpar.sval.'FILE_GFK_META'
    synseisfile = inpar.sval.'FILE_SYNTHETICS'
    eventfile = inpar.sval.'FILE_EVENT_LIST'
    stationfile = inpar.sval.'FILE_STATION_LIST'
    evfiltype = inpar.sval.'EVENT_FILTER_TYPE'
    evfilspecs => dvecp(inpar,'EVENT_FILTER_SPECS',4,ierr)
    if (ierr /= 0) then
       call add(errmsg,2,'Input for event filter specs cannot be read',myname)
       goto 10
    endif
    call printInputParameter(inpar)
    call dealloc(inpar)
!----------------------------------------------------------------
!  open HDF environment    
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read event and station file
!  use ASKI convention for both files
!
    call createEventListFromEventFile(eventfile,1,'ASKI_events',event_list,errmsg)
    if (.level.errmsg == 2) goto 10
    call createStationListFromStationFile(stationfile,1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  check whether csys in station and event file are consistent
!
    if (.csys.event_list /= .csys.station_list) then
       call add(errmsg,2,'CSYS inconsistent in event and station file',myname)
       goto 10
    endif
!
!  read GFK meta data, equal for all events and source depths
!
    call openFileRoHDFWrapper(gfkmetafile,fidmeta,errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'reals',errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fidmeta,'integers',errmsg)
    if (.level.errmsg == 2) goto 10
    call h5fclose_f(fidmeta,ierr)
    if (ierr < 0) then; errstr = 'failed to close GFK meta file'; goto 1; endif
!
!  Calculate event filter
!
    call computeEventFilter(evfiltype,evfilspecs,gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,evfilter,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  properties of time series and allocate
!
    dtp = 0.25d0/(2.d0*gfk%nf2*gfk%df)
    call newNsampDtFrequencyTime(gfk%df,dtp,nsamp,dt)
    allocate(urs(nsamp,1),zsp(gfk%nf2))
    zsp = 0.d0
!
!  open file for synthetics output
!    
    call createFileHDFWrapper(synseisfile,fidseis,errmsg)
    if (.level.errmsg == 2) goto 10
!--------------------------------------------------------------------------------------
!  Loop over events
!--------------------------------------------------------------------------------------
    do while (nextEventSeismicEventList(event_list,event))
       call new(errmsg1,myname)
       call printSeismicEvent(event)
   !   
   !  create a group for event in HDF output and write event info
   !
       call createEventSynseisHDF(event,fidseis,evgrid,errmsg1)
       if (.level.errmsg1 == 2) goto 11
    !
    !  Loop over stations
    !
       do while (nextStationSeismicNetwork(station_list,station,is))
          call new(errmsg2,myname)
          if (.level.errmsg2 == 2) goto 12
          if (is == 1) then
             zsp(:) = evfilter
          else if (is == 2) then
             do jf = gfk%nf1,gfk%nf2
                zom = dcmplx(mc_two_pid*(jf-1)*gfk%df,-gfk%sigma)
                zsp(jf) = mc_cid*zom*evfilter(jf)
             enddo
          endif
          call transformFrequencyTime(gfk%nf1,gfk%nf2,gfk%df,gfk%sigma,zsp,nsamp,dt,urs(:,1),errmsg2)
          if (.level.errmsg == 2) goto 12
       !
       !  write to HDF
       !
          call writeStationSynseisHDF(station,evgrid,nsamp,0.d0,dt,'Z',urs,errmsg2)
          if (.level.errmsg2 == 2) goto 12
          call dealloc(errmsg2)
       enddo                                                                         ! end station loop
       call dealloc(errmsg1)
       call h5gclose_f(evgrid,ierr)
       if (ierr < 0) then; errstr = "error closing event group"; goto 1; endif
    enddo                                                                             ! end event loop
!
!  deallocation
!
    call h5fclose_f(fidseis,ierr)
    if (ierr < 0) then; errstr = "error closing seis file"; goto 1; endif
    call h5close_f(ierr)
!
    deallocate(urs)
    deallocate(evfilter)
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
!
end program
