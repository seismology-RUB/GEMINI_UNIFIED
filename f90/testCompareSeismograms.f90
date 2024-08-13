! ============================================================================
!  Test: Compare synthetic seismogram gathers
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
!  Used for testing code after changes. Compares seismograms obtained
!  from a reference run and a repetition of it after changing the code.
!--------------------------------------------------------------------------
program testCompareSeismogramGather
    use readEventStationFile
    use asciiSynseisIO
    use argumentParser
    use errorMessage
    use inputParameter
    use string
    implicit none
    type (argument_parser) :: ap
    type (error_message) :: errmsg
    type (input_parameter) :: inpar
    type (seismic_event) :: event
    type (seismic_station) :: station
    type (seismic_event_list) :: event_list
    type (seismic_network) :: station_list
    real, dimension(:), allocatable :: urs,ursref
    character(len=22) :: myname = 'testCompareSeismograms'
    integer :: ic,nsamp,ios,cnt,cntsf
    real :: dt,rmssyn,rmsref,rmsdiff
    double precision :: tanfdp
    character(len=max_length_string) :: setup,eventfile,stationfile,gemini_parfile
    character(len=max_length_string) :: comps,path_synthetic_seis,staname,synpath,refpath
    character (len=1) :: component
    character (len=80), dimension(4) :: par_keys
    data par_keys/'PATH_SYNTHETIC_SEIS','FILE_EVENT_LIST','FILE_STATION_LIST','COMPONENTS'/
  !-------------------------------------------------------------------------------
    call init(ap,myname,'Compare seismogram gathers for testing purposes')
    call addPosarg(ap,'setup_dir','sval','Setup Folder (slash terminated)')
    call parse(ap)
    setup = ap.sval.'setup_dir'
    if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); stop; endif
    call document(ap); call dealloc(ap)
       
    call new(errmsg,myname)
  !-------------------------------------------------------------
  !  read input parameters from gemini_parfile
  !
    gemini_parfile = setup+'parfile_syn'
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,gemini_parfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    path_synthetic_seis = setup+(inpar.sval.'PATH_SYNTHETIC_SEIS')
    eventfile = setup+(inpar.sval.'FILE_EVENT_LIST')
    stationfile = setup+(inpar.sval.'FILE_STATION_LIST')
    comps = inpar.sval.'COMPONENTS'
    call printInputParameter(inpar)
    call dealloc(inpar)
  !----------------------------------------------------------------
  !  read event and station file
  !  use ASKI convention for event and station file
  !
    call createEventListFromEventFile(eventfile,1,'ASKI_events',event_list,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    call createStationListFromStationFile(stationfile,1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
  !
  !  check whether csys in station and event file are consistent
  !
    if (.csys.event_list /= .csys.station_list) then
       call add(errmsg,2,'CSYS inconsistent in event and station file',myname)
       call print(errmsg); stop
    endif
  !
  !  Loop over events
  !
    cnt = 0; cntsf = 0
    do while (nextEventSeismicEventList(event_list,event))
   !
   !  Loop over components and stations
   !
       do ic = 1,len_trim(comps)
          do while (nextStationSeismicNetwork(station_list,station))
         !
         !  read new synthetic data
         !
             synpath = path_synthetic_seis+(.evid.event)+'_'+(.staname.station)+'_'+comps(ic:ic)
             call readAsciiSynseisIO(1,synpath,nsamp,tanfdp,dt,staname,component,urs,ios)
             if (ios /= 0) then
                call add(errmsg,2,'Synthetic file could not be opened ->'+trim(synpath),myname)
                call print(errmsg); stop
             endif
             if (.not.(component.equal.comps(ic:ic))) then
                call add(errmsg,2,'Expected component does not match component in data file',myname)
                call print(errmsg); stop
             endif
         !
         !  read reference data
         !
             refpath = setup+'refseis/'+(.evid.event)+'_'+(.staname.station)+'_'+comps(ic:ic)
             call readAsciiSynseisIO(1,refpath,nsamp,tanfdp,dt,staname,component,ursref,ios)
             if (ios /= 0) then
                call add(errmsg,2,'Synthetic reference file could not be opened ->'+trim(refpath),myname)
                call print(errmsg); stop
             endif
             if (.not.(component.equal.comps(ic:ic))) then
                call add(errmsg,2,'Expected component does not match component in ref data file',myname)
                call print(errmsg); stop
             endif
          !
          !  compute rms of syn, ref, and difference
          !
             rmssyn = sqrt(sum(urs**2/nsamp))
             rmsref = sqrt(sum(ursref**2/nsamp))
             rmsdiff = sqrt(sum((urs-ursref)**2/nsamp))
             write(6,'(a,2e15.2)') trim((.evid.event)+'_'+(.staname.station)+'_'+comps(ic:ic)),rmsdiff/rmsref,rmsdiff/rmssyn
          !
          !  count cases and successful ones
          !
             if (rmsdiff/rmsref < 1.e-4) cntsf = cntsf+1
             cnt = cnt+1
             deallocate(urs,ursref)
          enddo
       enddo
    enddo
    write(6,'(i5,a,i5,a)') cntsf,' of ',cnt,' assertions passed'
    call dealloc(event_list)
    call dealloc(station_list)
end program
