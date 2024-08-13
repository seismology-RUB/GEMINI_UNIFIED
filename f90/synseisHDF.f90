! =================================================================================
!  Routines for reading and writing synthetic seismograms in HDF
! =================================================================================
!----------------------------------------------------------------------------
!   Copyright 2019 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!--------------------------------------------------------------------------------------
!   Reading and writing synthetic seismograms from or to a HDF file.
!   Events are created as a subgroup of the file.
!   Event info is written as attributes of the subgroup.
!   Synthetics seismograms for stations are written as datasets in the event subgroup.
!   Station info is written as attibutes of the dataset.
!   Reading routines recover event type, station type, components and seismograms.
!-------------------------------------------------------------------------------------
module synseisHDF
       use hdf5
    use constants
    use string
    use errorMessage
    use seismicEvent
    use seismicStation
    use dateTime
    use hdfWrapper
    use anyRankRealArray
    use anyRankIntegerArray
    implicit none
contains
!--------------------------------------------------------------
!  Create event sub group to parent
!  and write event info as attributes
!
    subroutine createEventSynseisHDF(event,fid,evgrid,errmsg)
    type(seismic_event) :: event
    integer(hid_t) :: fid
    integer(hid_t), intent(out) :: evgrid
    type(error_message) :: errmsg
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    real, dimension(:), allocatable, target :: d
    integer, dimension(:), allocatable, target :: id
    integer :: ierr
    character (len=date_time_string_length) :: dtstr
    character (len=21) :: myname = "createEventSynseisHDF"
!
    call h5gcreate_f(fid,'Event',evgrid,ierr)
    if (ierr < 0) goto 1
!
    call writeStringAttributeHDFWrapper(evgrid,'Event_id',trim(.evid.event),errmsg)
    if (.level.errmsg == 2) return
!
    dtstr = convertToFullTimestringDateTime(.otime.event)
    call writeStringAttributeHDFWrapper(evgrid,'Date_and_time',dtstr,errmsg)
    if (.level.errmsg == 2) return
!
    call writeStringAttributeHDFWrapper(evgrid,'Coordinate_system',.csys.event,errmsg) 
    if (.level.errmsg == 2) return
!
    id = (/.styp.event/); call aria%assoc1d(id)
    call writeArrayAttributeHDFWrapper(evgrid,'Source_type',aria,errmsg)
    if (.level.errmsg == 2) return; call aria%deassoc()
!
    d = (/.mag.event/); call arra%assoc1d(d)
    call writeArrayAttributeHDFWrapper(evgrid,'Magnitude',arra,errmsg)
    if (.level.errmsg == 2) return; call arra%deassoc()
!
    d = (/.slat.event,.slon.event,.sdepth.event/); call arra%assoc1d(d)
    call writeArrayAttributeHDFWrapper(evgrid,'Coordinates',arra,errmsg)
    if (.level.errmsg == 2) return; call arra%deassoc()
!
    if (.styp.event == 0) then
       d = .force.event; call arra%assoc1d(d)
       call writeArrayAttributeHDFWrapper(evgrid,'Force',arra,errmsg)
    else
       d = .momten.event; call arra%assoc1d(d)
       call writeArrayAttributeHDFWrapper(evgrid,'Moment',arra,errmsg)
    endif
    if (.level.errmsg == 2) return
    call arra%deassoc()

    d = (/.tshift.event/)
    call arra%assoc1d(d)
    call writeArrayAttributeHDFWrapper(evgrid,'Timeshift',arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%deassoc()
 !
    return
 1  call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine createEventSynseisHDF
!--------------------------------------------------------------
!  Open event sub group to file
!  and read event info from attributes
!
    subroutine openEventSynseisHDF(fid,evid,evgrid,event,errmsg)
    integer(hid_t) :: fid
    character (len=char_len_evid), intent(out) :: evid
    integer(hid_t), intent(out) :: evgrid
    type(seismic_event), intent(out) :: event
    type(error_message) :: errmsg
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    integer :: ierr,istyp,strlen
    character (len=date_time_string_length) :: dtstr
    character (len=1) :: csys
    integer, dimension(:), pointer :: id
    real, dimension(:), pointer :: d
    real :: mag,slon,slat,sdepth,timeshift
    type(date_time) :: otime
    character (len=19) :: myname = "openEventSynseisHDF"
    character(len=max_length_string) :: errstr,cval
!
    call h5gopen_f(fid,'Event',evgrid,ierr)                      ! open subgroup
    if (ierr < 0) then; errstr = "h5gopen"; goto 1; endif
!
    call readStringAttributeHDFWrapper(evgrid,'Event_id',cval,strlen,errmsg)
    if (.level.errmsg == 2) return
    evid = cval(1:strlen)
!
    call readStringAttributeHDFWrapper(evgrid,'Date_and_time',cval,strlen,errmsg)
    if (.level.errmsg == 2) return
    dtstr = cval(1:strlen)
    call createFromFullTimestringDateTime(otime,dtstr)
!
    call readStringAttributeHDFWrapper(evgrid,'Coordinate_system',cval,strlen,errmsg)
    if (.level.errmsg == 2) return
    csys = cval(1:strlen)
!
    call readArrayAttributeHDFWrapper(evgrid,'Magnitude',arra,errmsg)
    if (.level.errmsg == 2) return
    d => arra%get1d(); mag = d(1); deallocate(d)
!
    call readArrayAttributeHDFWrapper(evgrid,'Coordinates',arra,errmsg)
    if (.level.errmsg == 2) return
    d => arra%get1d(); slat = d(1); slon = d(2); sdepth = d(3); deallocate(d)
!
    call readArrayAttributeHDFWrapper(evgrid,'Source_type',aria,errmsg)
    if (.level.errmsg == 2) return
    id => aria%get1d(); istyp = id(1); deallocate(id)
!
    call readArrayAttributeHDFWrapper(evgrid,'Timeshift',arra,errmsg)
    if (.level.errmsg == 2) return
    d => arra%get1d(); timeshift = d(1); deallocate(d)
!
    if (istyp == 0) then                                                  ! moment or force
       call readArrayAttributeHDFWrapper(evgrid,'Force',arra,errmsg)
    else
       call readArrayAttributeHDFWrapper(evgrid,'Moment',arra,errmsg)
    endif
    if (.level.errmsg == 2) return
    d => arra%get1d()
    call createForceMomentSeismicEvent(event,evid,csys,slat,slon,sdepth,mag,otime,timeshift,istyp,d)
    deallocate(d)
 !
    return
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine openEventSynseisHDF
!--------------------------------------------------------------
!  Create station dataset in file
!  and write station info as attributes
!
    subroutine writeStationSynseisHDF(station,fid,nsamp,tanf,dt,comps,urs,errmsg,xferprp)
    type(seismic_station) :: station
    integer(hid_t) :: fid
    integer :: nsamp
    double precision :: tanf,dt
    character (len=*) :: comps
    real, dimension(:,:), target :: urs
    real, dimension(:), allocatable, target :: d
    type (any_rank_real_array) :: ara
    type(error_message) :: errmsg
    integer(hid_t), optional :: xferprp
    integer :: ierr
    integer(hid_t) :: dsetid,xfer
    character (len=22) :: myname = "writeStationSynseisHDF"
    character(len=max_length_string) :: errstr
!
    if (present(xferprp)) then; xfer = xferprp; else; xfer = H5P_DEFAULT_F; endif
!
    call ara%assoc2d(urs)
    call writeArrayHDFWrapper(fid,trim(.netcode.station)//'.'//trim(.staname.station),ara,errmsg,ds = dsetid,xferprp=xfer)
    if (.level.errmsg == 2) return; call ara%deassoc()
!
    call writeStringAttributeHDFWrapper(dsetid,'Staname',trim(.staname.station),errmsg)
    if (.level.errmsg == 2) return
!
    call writeStringAttributeHDFWrapper(dsetid,'Network',trim(.netcode.station),errmsg)
    if (.level.errmsg == 2) return
!
    call writeStringAttributeHDFWrapper(dsetid,'Components',trim(comps),errmsg)
    if (.level.errmsg == 2) return
!
    call writeStringAttributeHDFWrapper(dsetid,'Coordinate_system',.csys.station,errmsg)
    if (.level.errmsg == 2) return
!
    d = (/.lat.station,.lon.station,.alt.station/); call ara%assoc1d(d)
    call writeArrayAttributeHDFWrapper(dsetid,'Coordinates',ara,errmsg)
    if (.level.errmsg == 2) return; call ara%deassoc()
!
    d = (/real(tanf),real(dt)/); call ara%assoc1d(d)
    call writeArrayAttributeHDFWrapper(dsetid,'Time',ara,errmsg)
    if (.level.errmsg == 2) return; call ara%deassoc()
!
    deallocate(d)
    call h5dclose_f(dsetid,ierr)
    if (ierr < 0) then; errstr = 'h5dclose'; goto 1; endif
!
    return
 1  call add(errmsg,2,errstr,myname)
    end subroutine writeStationSynseisHDF
!---------------------------------------------------------------------
!  Read station dataset in evgrid
!  and read station info from attributes
!
    subroutine readStationSynseisHDF(fid,stdsetname,station,comps,tanf,dt,urs,errmsg)
    integer(hid_t) :: fid
    character(len=*) :: stdsetname
    type(seismic_station) :: station
    character (len=:), allocatable :: comps
    double precision :: tanf,dt
    real, dimension(:,:), pointer :: urs
    type(error_message) :: errmsg
    type (any_rank_real_array) :: ara
    integer :: strlen1,strlen2,strlen3,strlen4,ierr
    integer(hid_t) :: dsetid
    real, dimension(:), pointer :: d
    real :: lat,lon,alt
    character (len=char_len_sta) :: staname
    character (len=char_len_netcode) :: netcode
    character (len=1) :: csys
    character (len=21) :: myname = "readStationSynseisHDF"
    character (len=max_length_string) :: errstr,cval
!    
    call readArrayHDFWrapper(fid,trim(stdsetname),ara,errmsg,ds = dsetid)
    if (.level.errmsg == 2) return
    urs => ara%get2d(); call ara%deassoc()
!
    call readStringAttributeHDFWrapper(dsetid,'Staname',cval,strlen1,errmsg)
    if (.level.errmsg == 2) return
    staname = cval(1:strlen1)
!
    call readStringAttributeHDFWrapper(dsetid,'Network',cval,strlen2,errmsg)
    if (.level.errmsg == 2) return
    netcode = cval(1:strlen2)
!
    call readStringAttributeHDFWrapper(dsetid,'Components',cval,strlen3,errmsg)
    comps = cval(1:strlen3)
    if (.level.errmsg == 2) return
!
    call readStringAttributeHDFWrapper(dsetid,'Coordinate_system',cval,strlen4,errmsg)
    if (.level.errmsg == 2) return
    csys = cval(1:strlen4)
!
    call readArrayAttributeHDFWrapper(dsetid,'Coordinates',ara,errmsg)
    if (.level.errmsg == 2) return
    d => ara%get1d(); lat = d(1); lon = d(2); alt = d(3); deallocate(d)
!
    call readArrayAttributeHDFWrapper(dsetid,'Time',ara,errmsg)
    if (.level.errmsg == 2) return
    d => ara%get1d(); tanf = dble(d(1)); dt = dble(d(2)); deallocate(d)
!
    call h5dclose_f(dsetid,ierr)
    if (ierr < 0) then; errstr = "h5dclose"; goto 1; endif
!
    call createSeismicStation(station,csys(1:strlen4),staname(1:strlen1),lat,lon,alt = alt,netcode = netcode(1:strlen2))
!
    return
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine readStationSynseisHDF
!-------------------------------------------------------------
!  Get number of stations of event
!
    subroutine getNumberStationsSynseisHDF(fid,nsta,errmsg)
    integer(hid_t) :: fid
    type(error_message) :: errmsg
    integer :: nsta
    integer :: storage_type, nd, max_corder,ierr
    character (len=27) :: myname = "getNumberStationsSynseisHDF"
!
    call h5gget_info_f(fid,storage_type,nd,max_corder,ierr)
    if (ierr /= 0) goto 1
    nsta = nd-1
1   if (ierr < 0) then
       call add(errmsg,2,'Failed in '+myname,myname)
    endif
    end subroutine getNumberStationsSynseisHDF
!-------------------------------------------------------------
!  Get station names = dataset names in file
!
    subroutine getDatasetNamesSynseisHDF(fid,stdsetname,nsta,errmsg)
    integer(hid_t) :: fid
    character(len=char_len_sta+char_len_netcode+1), dimension(:), allocatable :: stdsetname
    integer :: nsta
    type(error_message) :: errmsg
    integer :: ierr
    integer(hsize_t) :: n
    integer :: storage_type,nd,max_corder,js
    character (len=char_len_sta) :: linkname
    character (len=29) :: myname = "getDatasetNamesSynseisHDF"
!
    call h5gget_info_f(fid,storage_type,nd,max_corder,ierr)
    if (ierr /= 0) goto 1
    nsta = nd-1
    print *,'Number of data sets in file: ',nsta
    allocate(stdsetname(nsta))
!
!  Loop through all links in file, get dataset names but exclude the event group
!  Link indexing starts at zero.
!
    js = 0
    do n = 0,nd-1
       call h5lget_name_by_idx_f(fid,'.',H5_INDEX_NAME_F,H5_ITER_INC_F,n,linkname,ierr)
       if (ierr /= 0) goto 1
       if (equalString(linkname,'Event')) then
          print *,'Event group link found'
       else
          js = js+1
          stdsetname(js) = linkname
       end if
    end do
1   if (ierr /= 0) then
       call add(errmsg,2,'Failed in '+myname,myname)
    endif
    end subroutine getDatasetNamesSynseisHDF
!-------------------------------------------------------------------
!  Get entire station object of given station in file
!
    subroutine getStationInfoSynseisHDF(fid,stdsetname,station,errmsg)
    integer(hid_t) :: fid
    character(len=*) :: stdsetname
    type (seismic_station) :: station
    type(error_message) :: errmsg
    character (len=char_len_netcode) :: netcode
    character (len=char_len_sta) :: staname
    character (len=max_length_string) :: errstr,cval
    integer(hid_t) :: dset
    integer :: ierr,strlen1,strlen2,strlen3
    type (any_rank_real_array) :: ara
    real, dimension(:), pointer :: d
    real :: lat,lon,alt
    character (len=1) :: csys
    character (len=24) :: myname = "getStationInfoSynseisHDF"
!
    call h5dopen_f(fid,trim(stdsetname),dset,ierr)
    if (ierr < 0) then; errstr = 'h5dopen'; goto 1; endif
!
    call readStringAttributeHDFWrapper(dset,'Staname',cval,strlen3,errmsg)
    if (.level.errmsg == 2) return
    staname = cval(1:strlen3)
!
    call readStringAttributeHDFWrapper(dset,'Network',cval,strlen1,errmsg)
    if (.level.errmsg == 2) return
    netcode = cval(1:strlen1)
!
    call readStringAttributeHDFWrapper(dset,'Coordinate_system',csys,strlen2,errmsg)
    if (.level.errmsg == 2) return
    csys = cval(1:strlen2)
!
    call readArrayAttributeHDFWrapper(dset,'Coordinates',ara,errmsg)
    if (.level.errmsg == 2) return
    d => ara%get1d(); lat = d(1); lon = d(2); alt = d(3); deallocate(d)
!
    call h5dclose_f(dset,ierr)
    if (ierr < 0) then; errstr = "h5dclose"; goto 1; endif
!
!  need to pass one more character than string length to provide space for end-of-string character
!
    call createSeismicStation(station,csys(1:strlen2),staname(1:strlen3),lat,lon,alt = alt,netcode = netcode(1:strlen1))
    return
!
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine getStationInfoSynseisHDF
!
end module synseisHDF
