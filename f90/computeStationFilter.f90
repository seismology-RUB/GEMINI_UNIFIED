! ===============================================================================
!  Main program for the computation of station filters
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
!--------------------------------------------------------------------------
!   Compute station filters for all stations in station_file.
!   This must be done before computing synthetic seismograms.
!   Parameters for this program are taken from the PARFILE_SYN parameter file.
!   Code currently offers trivial unit filter. Instrument responses
!   should be added later.
!   The code opens the GFDSV fk-spectra file to look up the required frequency
!   range for which filter values are computed. Note all filter values
!   are computed at frequencies with constant negative imaginary part
!   also taken from the GFDSV fk-spectra file.
!---------------------------------------------------------------------------
program computeStationFilter
       use hdf5
    use mathConstants
    use argumentParser
    use string
    use errorMessage
    use inputParameter
    use readEventStationFile
    use seismicNetwork
    use hdfWrapper
    use anyRankRealArray
    use greenFKSpectra
!
    implicit none
    type (argument_parser) :: ap
    type (error_message) :: errmsg
    type (seismic_station) :: station
    type (seismic_network) :: station_list
    type (input_parameter) :: inpar
    type (green_fk_spectra) :: gfk
    type (any_rank_real_array) :: ara
!    
    integer, parameter :: NC = 3
    integer :: nf,j,ierr
    integer(hid_t) :: fid
    logical :: use_meta
    real, dimension(:,:,:), allocatable :: hfil
    character(len=20) :: myname = 'computeStationFilter'
    character(len=max_length_string) :: stationfile,parfile,filfile,dsvbasename
    character(len=max_length_string) :: gfkfile,filtertype,errstr,sourcetype,dsetname
    character(len=80), dimension(5) :: par_keys
!
!  keywords for main parameters
!
    data par_keys/'FILE_STATION_LIST','FILE_STATION_FILTER',&
         'DSVBASENAME','STATION_FILTER_TYPE','SOURCE_TYPE'/
!------------------------------------------------------------------
    call init(ap,myname,'Compute station filters at frequencies specified in Green FK spectra for all stations')
    call addPosarg(ap,'parfile_synseis','sval','Parameter file needed for synthetic seismograms')
    call addOption(ap,'-meta',.false.,'read Green FK meta file')
    call parse(ap)
    parfile = ap.sval.'parfile_synseis'
    use_meta = ap.optset.'-meta'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call document(ap); call dealloc(ap)
!-------------------------------------------------
    call new(errmsg,myname)
!------------------------------------------------------------
!  read input parameters from parameter file
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    stationfile = inpar.sval.'FILE_STATION_LIST'
    dsvbasename = inpar.sval.'DSVBASENAME'
    sourcetype = inpar.sval.'SOURCE_TYPE'
    if (use_meta) then
       gfkfile = dsvbasename+'.meta'
    else
       gfkfile = dsvbasename+'.'+sourcetype
    endif
    filfile = inpar.sval.'FILE_STATION_FILTER'
    filtertype = inpar.sval.'STATION_FILTER_TYPE'
    call print(inpar); call dealloc(inpar)
!------------------------------------------------------------------------------
!  read station file
!
    call createStationListFromStationFile(stationfile,1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  open HDF environment
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
!-------------------------------------------------------------------------------
!  extract frequency information from GFDSV fk-spectra file
!  accessible through members of gfk
!
    call openFileRoHDFWrapper(gfkfile,fid,errmsg)
    if (.level.errmsg == 2) goto 10
    call readMetaGreenFKSpectra(gfk,fid,'integers',errmsg)
    if (.level.errmsg == 2) goto 10
    nf = gfk%nf2-gfk%nf1+1
    call h5fclose_f(fid,ierr)
!
!  output file for filter data
!
    call createFileHDFWrapper(filfile,fid,errmsg)
    if (.level.errmsg == 2) goto 10
!--------------------------------------------------------------------------------
!  loop over stations
!
    allocate(hfil(2,nf,NC))   
    call ara%assoc3d(hfil)
    do while (nextStationSeismicNetwork(station_list,station))
       do j = 1,NC
          select case (filtertype)
          case ('UNIT')
             hfil(1,:,j) = 1.0; hfil(2,:,j) = 0.0
          case default
             call add(errmsg,2,'Unknown station filter type specified',myname)
          end select
       enddo
       dsetname = (.netcode.station)+'.'+(.staname.station)                            ! next component
       call writeArrayHDFWrapper(fid,trim(dsetname),ara,errmsg)
       if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    enddo                                                                          ! next station
    call ara%deassoc()   
    call h5fclose_f(fid,ierr)
    if (ierr < 0) then; errstr = 'h5fclose'; goto 1; endif
    call h5close_f(ierr)
    if (ierr < 0) then; errstr = 'h5close'; goto 1; endif
    deallocate(hfil)
!
!  treat HFD error
!
 1  if (ierr < 0) then
       call add(errmsg,2,trim(errstr),myname)
       call print(errmsg)
       stop
    endif
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif
!
 end program computeStationFilter
