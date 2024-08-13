! ===============================================================================
!  Main program for the computation of event filters
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
!   Compute event filters for all events in event_file.
!   This must be done before computing synthetic seismograms.
!   Parameters for this program are taken from the GFDSVSEIS parameter file.
!   Code currently offers a Butterworth bandpass filter. Other options
!   may be added later. Selection will be done by keyword in the parameter 
!   file.
!   The code opens the GFDSV fk-spectra file to look up the required frequency
!   range for which filter values are computed. Note all filter values
!   are computed at frequencies with constant negative imaginary part
!   also taken from the GFDSV fk-spectra file.
!---------------------------------------------------------------------------
program computeEventFilter
      use hdf5
    use mathConstants
    use argumentParser
    use string
    use errorMessage
    use inputParameter
    use readEventStationFile
    use seismicEventList
    use hdfWrapper
    use anyRankRealArray
    use greenFKSpectra
    use butterworthFilter
!
    implicit none
    type (argument_parser) :: ap
    type (error_message) :: errmsg
    type (seismic_event) :: event
    type (seismic_event_list) :: event_list
    type (input_parameter) :: inpar
    type (any_rank_real_array) :: ara
!    
    integer :: nf,nf1,nf2,nord,ios,ierr
    integer(hid_t) :: fid
    double precision :: fc,sigma,df
    real, dimension(:,:), allocatable :: hfil
    real, dimension(:), pointer :: evfilspecs
    double complex, dimension(:), allocatable :: hlow,hhigh
    character(len=18) :: myname = 'computeEventFilter'
    character(len=max_length_string) :: eventfile,parfile,filfile
    character(len=max_length_string) :: gfkfile,evfiltype,errstr
    character(len=80), dimension(5) :: par_keys
!
!  keywords for main parameters
!
    data par_keys/'FILE_EVENT_LIST','FILE_EVENT_FILTER',&
         'GFKFILE','EVENT_FILTER_TYPE','EVENT_FILTER_SPECS'/
!------------------------------------------------------------------
    call init(ap,myname,'Compute event filters at frequencies specified in Green FK spectra for all events')
    call addPosarg(ap,'parfile_synseis','sval','Parameter file needed for synthetic seismograms')
    call parse(ap)
    parfile = ap.sval.'parfile_synseis'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call document(ap); call dealloc(ap)
!-------------------------------------------------
    call new(errmsg,myname)
    nullify(evfilspecs)
!------------------------------------------------------------
!  read input parameters from parameter file
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    eventfile = inpar.sval.'FILE_EVENT_LIST'
    filfile = inpar.sval.'FILE_EVENT_FILTER'
    gfkfile = inpar.sval.'GFKFILE'
    evfiltype = inpar.sval.'EVENT_FILTER_TYPE'
    evfilspecs => rvecp(inpar,'EVENT_FILTER_SPECS',4,ios)
    if (ios /= 0) then
       call add(errmsg,2,'Input for event filter specs cannot be read',myname)
       call print(errmsg); stop
    endif
    call print(inpar); call dealloc(inpar)
!------------------------------------------------------------------------------
!  read event file
!
    call createEventListFromEventFile(eventfile,1,'ASKI_events',event_list,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!-------------------------------------------------------------------------------
!  extract frequency information from GFDSV fk-spectra file
!  accessible through members of gfk
!
    call extractFrequencyInformationGreenFKSpectra(gfkfile,df,sigma,nf1,nf2,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    nf = nf2-nf1+1
!
!  open HDF environment and open filter HDF file
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    call createFileHDFWrapper(filfile,fid,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!--------------------------------------------------------------------------------
!  loop over events
!
    allocate(hfil(2,nf))
    call ara%assoc2d(hfil)   
    do while (nextEventSeismicEventList(event_list,event))
   !
   !  bifurcation according to event filter type
   !  BUTTERWORTH_BANDPASS
   !
       select case (evfiltype)
       case ('BUTTERWORTH_BANDPASS')
   !
   !  calculate frequency response of HIGH-pass Butterworth filter
   !  at selected complex frequencies with fixed imaginary part
   !
          nord = nint(evfilspecs(1))
          fc   = evfilspecs(2)
          call highPassComplexOmegaButterworthFilter(nord,nf,fc,(nf1-1)*df,df,sigma,hhigh)
   !
   !  calculate frequency response of LOW-pass Butterworth filter
   !  at selected complex frequencies with fixed imaginary part
   !
          nord = nint(evfilspecs(3))
          fc   = evfilspecs(4)
          call lowPassComplexOmegaButterworthFilter(nord,nf,fc,(nf1-1)*df,df,sigma,hlow)
   !
   !  multiply high and low pass and write to file
   !
          hfil(1,:) = real(hhigh*hlow)
          hfil(2,:) = aimag(hhigh*hlow)
          deallocate(hlow,hhigh)
   !
   !  default case if filter type is not implemented or invalid
   !          
       case default
          call add(errmsg,2,'Unknown event filter type specified',myname)
       end select
       call writeArrayHDFWrapper(fid,trim(.evid.event),ara,errmsg)
       if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    enddo                            ! next event
    call ara%deassoc()   
    call h5fclose_f(fid,ierr)
    if (ierr < 0) then; errstr = 'h5fclose'; goto 1; endif
    call h5close_f(ierr)
    if (ierr < 0) then; errstr = 'h5close'; goto 1; endif
!
!  treat HFD error
!
 1  if (ierr < 0) then
       call add(errmsg,2,trim(errstr),myname)
       call print(errmsg)
       stop
    endif
!
    deallocate(hfil)
    end program computeEventFilter
