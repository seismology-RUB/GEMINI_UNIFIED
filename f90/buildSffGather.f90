! ==============================================================================
!  Build an seismic gather in SFF format from indiidual seismograms
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
!---------------------------------------------------------------------
!  Collect individual synthetic seismograms to build a seismic gather
!  according to specific selection criteria
!-----------------------------------------------------------------------------
 program buildSffGather
    use argumentParser
    use readEventStationFile
    use asciiSynseisIO
    use errorMessage
    use sffHeader
    use sffDatablock
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
    type (sff_header) :: head
    type (sff_datablock) :: dbl
    type (date_time) :: orgtime
    real, dimension(:), allocatable :: urs
    integer :: ic,nsamp,ios,j
    real :: dt
    double precision :: tanfdp
    logical :: last
    character(len=max_length_string) :: text,eventfile,stationfile,gemini_parfile, &
         path_synthetic_seis,gather_type,staname,gatherfile
    character(len=max_length_string) :: comps
    character (len=1) :: component
    character (len=21) :: myname = 'buildSeismogramGather'
    character (len=80), dimension(4) :: par_keys
!
!  keywords for input parameters
!
    data par_keys/'PATH_SYNTHETIC_SEIS','FILE_EVENT_LIST','FILE_STATION_LIST','COMPONENTS'/
!-------------------------------------------------------------------------------
    call init(ap,myname,'Build up a seismic gather from individual synthetic traces')
    call addPosarg(ap,'gemini_parfile','sval','GFDSVSEIS related parameter file')
    call addOption(ap,'-gt',.true.,'gather type','sval','SOURCE_GATHER')
    call parse(ap)
    gemini_parfile = ap.sval.'gemini_parfile'
    gather_type = ap.sval.'-gt'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call document(ap); call dealloc(ap)
!
    call new(errmsg,myname)
!-------------------------------------------------------------
!  read input parameters from gemini_parfile
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,gemini_parfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    path_synthetic_seis = inpar.sval.'PATH_SYNTHETIC_SEIS'
    eventfile = inpar.sval.'FILE_EVENT_LIST'
    stationfile = inpar.sval.'FILE_STATION_LIST'
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
!  check whether csys im station and event file are consistent
!
    if (.csys.event_list /= .csys.station_list) then
       call add(errmsg,2,'CSYS inconsistent in event and station file',myname)
       call print(errmsg); stop
    endif
!---------------------------------------------------------------------
!  Bifurcate according to gather type
!  For SOURCE_GATHER, collect all traces of one event
!
    select case (trim(gather_type))
    case ('SOURCE_GATHER')
   !
   !  Loop over events
   !
       do while (nextEventSeismicEventList(event_list,event))
          gatherfile = path_synthetic_seis+(.evid.event)+'_gather.sff'
          call basicSFFHeader(head)
          call writeSFFHeader(head,2,gatherfile)
          orgtime = .otime.event
          last = .false.
   !
   !  Loop over components and stations
   !
          do ic = 1,len_trim(comps)
             j = 0
             do while (nextStationSeismicNetwork(station_list,station))
                j = j+1
                if (j == (.nstat.station_list) .and. ic == len_trim(comps)) last = .true.
                text = path_synthetic_seis+(.evid.event)+'_'+(.staname.station)+'_'+comps(ic:ic)
                call readAsciiSynseisIO(1,text,nsamp,tanfdp,dt,staname,component,urs,ios)
                if (ios /= 0) then
                   call add(errmsg,2,'Synthetic file could not be opened ->'+trim(text),myname)
                   call print(errmsg); stop
                endif
                if (.not.(component.equal.comps(ic:ic))) then
                   call add(errmsg,2,'Expected component does not match component in data file',myname)
                   call print(errmsg); stop
                endif
       !
       !  add seismogram to output gather file
       !
                call basicSFFDatablock(dbl,nsamp,dt,urs,'NSP',.staname.station,'LH'+component,.true.)
                call modifyDateSFFDatablock(dbl,.year.orgtime,.month.orgtime,.day.orgtime)
                call modifyTanfSFFDatablock(dbl,real(.tanfdp.orgtime))
                call writeSFFDataBlock(dbl,2,last)
                deallocate(urs)
             enddo
          enddo
       enddo
    case default
       call add(errmsg,2,'Sorry: Gather type is not implemented',myname)
       call print(errmsg); stop
    end select
!
!  clean up
!
    call dealloc(event_list); call dealloc(station_list)
    call print(errmsg); call dealloc(errmsg)
 end program buildSffGather
