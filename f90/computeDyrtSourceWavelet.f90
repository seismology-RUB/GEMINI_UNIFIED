! ==========================================================================================
!  Main program for the computation of source wavelets using ray travel times and amplitudes
! ==========================================================================================
!----------------------------------------------------------------------------
!   Copyright 2021 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!----------------------------------------------------------------------------------
program computeDyrtSourceWavelet
    use mathConstants
    use constants
    use argumentParser
    use inputParameter
    use string
    use miniSeed
    use errorMessage
    use hdfWrapper
    use timeSeries
    use fourierSpectrum
    use frequencyTime
    use readEventStationFile
    use asciiSynseisIO
    implicit none
!
    type (input_parameter) :: inpar
    type (error_message) :: errmsg,errmsg1
    type (mini_seed) :: ms
    type (argument_parser) :: ap
    type (seismic_event_list) :: event_list
    type (seismic_event) :: event
    type (seismic_network) :: station_list
    type (seismic_station) :: station
    type (seismic_station), dimension(:), allocatable :: useful_station
    type (time_series) :: stf,htstf,draw,fts0,fts1,fts2,fts3,fts4,fts5,fts6
    type (time_series), dimension(:), pointer :: dat,conv,untruncated_dat,avcb_dat,avcb_conv
    type (fourier_spectrum) :: ftstf
    type (fourier_spectrum), dimension(:), pointer :: ftdat
    type (any_rank_real_array) :: arra
    integer(kind=8) :: fid
    integer :: nsta,j,ierr,ntraces,n,nsyn,nrp,ndec,nrun,count_event,nfft
    integer, dimension(:),allocatable :: idx
    double precision, dimension(:), pointer :: evfilspecs
    double precision, dimension(:), allocatable :: chitr,snr,rmsdat,rmssyn,rmsnoise,rmsdaconv
    double precision, dimension(:), allocatable :: azre,azim,tanf
    double precision :: chi2,ta,dtdata,tsbuf,twinlen,esig,enoise
    double precision :: snrbound,misfitbound,rmsbound,visual_check_bound,taplen
    double precision :: rmsrefdat,frms,eps,azr,azi,anr,ani,aer,aei
    real, dimension(:), pointer :: d
    real :: tt
    logical :: ex,single_event,scaleflag,lastrun
    character (len=char_len_netcode+char_len_sta+1) :: nsname
    character (len=char_len_netcode+char_len_sta+1), dimension(:), allocatable :: stnetnam,avcbnam
    character (len=char_len_netcode+char_len_sta+1), dimension(:), pointer :: stexclude
    character(len=max_length_string) :: parfile,eventfile,evfiltype,stationfile,ttampfile,&
                                        path_to_events,path_to_data,datafile,excludefile,winlenfile,&
                                        tracefile,path_measured_data,errstr,logfile,evsel
    character (len=2) :: crun
    character (len=24) :: myname = 'computeDyrtSourceWavelet'
    character (len=80), dimension(16) :: par_keys
    integer, parameter :: lulog = 7
!
!  keywords for input parameters
!
    data par_keys/'FILE_EVENT_LIST','EVENT_FILTER_TYPE','EVENT_FILTER_SPECS','FILE_STATION_LIST',&
                  'DT_DATA','NFFT','TIME_SHIFT_BUFFER','FILE_SYN_TTIME_AMPLITUDE','STF_SNR_BOUND',&
                  'PATH_MEASURED_SEIS','DATA_PATH_TO_EVENTS','DATA_SUBEVENT_PATH','STF_MISFIT_BOUND',&
                  'STF_VISUAL_CHECK_BOUND','STF_RMS_BOUND','END_TAPER_LENGTH'/
!------------------------------------------------------------------------------------------------
    call init(ap,myname,'Compute source wavelet from ray synthetic traces and observed traces')
    call addPosarg(ap,'parfile','sval','DYRT parameter file')
    call addPosarg(ap,'pwl','dval','Phase window length of data and synthetics')
    call addOption(ap,'-single',.true.,'Run for single event','sval','None')
    call parse(ap)
    parfile = ap.sval.'parfile'
    twinlen = ap.dval.'pwl'
    evsel = ap.sval.'-single'
    single_event = ap.optset.'-single'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
!-------------------------------------------------------------
!  read input parameters from parfile
!
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 1
    eventfile = inpar.sval.'FILE_EVENT_LIST'
    evfiltype = inpar.sval.'EVENT_FILTER_TYPE'
    evfilspecs => dvecp(inpar,'EVENT_FILTER_SPECS',4,ierr)
    stationfile = inpar.sval.'FILE_STATION_LIST'
    dtdata = inpar.dval.'DT_DATA'
    nfft = inpar.ival.'NFFT'
    ttampfile = inpar.sval.'FILE_SYN_TTIME_AMPLITUDE'
    tsbuf = inpar.dval.'TIME_SHIFT_BUFFER'
    taplen = inpar.dval.'END_TAPER_LENGTH'
    snrbound = inpar.dval.'STF_SNR_BOUND'
    misfitbound = inpar.dval.'STF_MISFIT_BOUND'
    rmsbound = inpar.dval.'STF_RMS_BOUND'
    visual_check_bound = inpar.dval.'STF_VISUAL_CHECK_BOUND'
    path_measured_data = inpar.sval.'PATH_MEASURED_SEIS'
    path_to_events = inpar.sval.'DATA_PATH_TO_EVENTS'
    path_to_data = inpar.sval.'DATA_SUBEVENT_PATH'
    if (ierr /= 0) then
       call add(errmsg,2,'Input for event filter specs cannot be read',myname)
       goto 1
    endif
!----------------------------------------------------------------
!  open HDF environment
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read event file
!
    call createEventListFromEventFile(trim(eventfile),1,'dummy',event_list,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read station file (applicable for all events)
!
    call createStationListFromStationFile(trim(stationfile),1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    nsta = .nstat.station_list
!
!  open ttime and amplitude file
!
    call openFileRoHDFWrapper(ttampfile,fid,errmsg)
    if (.level.errmsg == 2) goto 10
!--------------------------------------------------------------------------------------
!  Loop over events
!--------------------------------------------------------------------------------------
    count_event = 0
    do while (nextEventSeismicEventList(event_list,event))
       if (single_event) then                               ! selected event desired
          if ((.evid.event).equal.evsel) then
             count_event = count_event+1
          else
             cycle
          endif
       else
          count_event = count_event+1
       endif
       call new(errmsg1,myname)
       call printSeismicEvent(event)
       logfile = trim(path_measured_data)//trim(.evid.event)//'/computeRaySourceWavelet.log'
       open(lulog,file = trim(logfile), status = 'unknown',iostat = ierr)
       if (ierr /= 0) then; errstr = 'cannot open logfile: '//trim(logfile); goto 1; endif
       call document(ap,lulog)
       call printInputParameter(inpar,lulog)
    !
    !  write true phase window length to file for later use with injection seismograms
    !
       winlenfile = trim(path_measured_data)//trim(.evid.event)//'/phasewinlen.txt'
       open(1,file = trim(winlenfile), status = 'unknown',iostat = ierr)
       if (ierr /= 0) then; errstr = 'cannot open winlenfile: '//trim(winlenfile); goto 1; endif
       write(1,*) twinlen
       close(1)
    !
    !  read exclude list if present, else set stexclude to null pointer
    !
       excludefile = trim(path_measured_data)//trim(.evid.event)//'/excluded-stations'
       inquire(file = trim(excludefile),EXIST = ex)
       if (ex) then
          call readExcludeFile(trim(excludefile),stexclude,ierr)
          if (ierr /= 0) then; errstr = 'reading of exclude file not successful'; goto 1; endif
       else
          nullify(stexclude)
          write(lulog,*) 'There is no excluded-stations file yet'
       endif
    !
    !  allocate station specific variables
    !
       allocate(dat(nsta),ftdat(nsta),untruncated_dat(nsta),stnetnam(nsta),snr(nsta),useful_station(nsta),&
                rmsdat(nsta),rmssyn(nsta),azre(nsta),azim(nsta),tanf(nsta),&
                rmsnoise(nsta),idx(nsta))
    !
    !  Loop over stations
    !
       ntraces = 0
       do while (nextStationSeismicNetwork(station_list,station,j))
          nsname = trim(.netcode.station)//'.'//(.staname.station)
          if (checkExcludeStations(trim(nsname),stexclude)) then
             write(lulog,'(a,a)') trim(nsname),' excluded'
             cycle
          endif
       !
       !  read data file if it exists
       !
          datafile = trim(path_to_events)//trim(.evid.event)//'/'//trim(path_to_data)
          datafile = trim(datafile)//trim(nsname)
          tracefile = trim(datafile)//'..HHZ'
          inquire(file = trim(tracefile),EXIST = ex)
          if (.not. ex) then; tracefile = trim(datafile)//'.00.HHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (.not. ex) then; tracefile = trim(datafile)//'.01.HHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (.not. ex) then; tracefile = trim(datafile)//'.02.HHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (.not. ex) then; tracefile = trim(datafile)//'..CHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (.not. ex) then; tracefile = trim(datafile)//'.00.BHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (.not. ex) then; tracefile = trim(datafile)//'.00.CHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (.not. ex) then; tracefile = trim(datafile)//'..BHZ'; inquire(file = trim(tracefile),EXIST = ex); endif
          if (ex) then
             write(lulog,*) trim(tracefile)
             call extractTimeWindowMiniSeed(ms,trim(tracefile),errmsg1)
             if (.level.errmsg1 == 2) goto 20
             call printMiniSeed(ms,lulog)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !  read travel time and complex ray amplitude for desired station and component here
          !
             call readRealScalarHDFWrapper(fid,trim(.evid.event)//'/'//trim(nsname)//'/'//'travel_time',tt,errmsg1)
             if (.level.errmsg1 == 2) goto 20
             call readArrayHDFWrapper(fid,trim(.evid.event)//'/'//trim(nsname)//'/'//'vectorial_amplitude',arra,errmsg1)
             if (.level.errmsg1 == 2) goto 20
             d => arra%get1d()
             azr = d(1)*1.e-15; azi = d(4)*1.e-15; anr = d(2); ani = d(5); aer = d(3); aei = d(6); deallocate(d)
             ta = .tanfdp.(.otime.event)+tt-tsbuf
             write(lulog,'(a,f15.3,a,e15.3)') 'Travel time: ',tt,' R-Amplitude(3): ',azr
             write(lulog,'(a,f15.3,a,f15.3)') 'Start time syn: ',ta,' Start time data: ',.tanfdp.ms
             write(lulog,'(a,f15.3,a,f15.3)') 'End   time syn: ',ta+tsbuf+twinlen,' End time data: ',.tanfdp.ms+(.nsamp.ms-1)*.dt.ms
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !
          !  check if data are available for entire time window from tt-tsbuf to tt+twinlen
          !  exclude station of not
          !
             if (.tanfdp.ms > ta) then                          ! Data start too late
                write(lulog,*) 'Data start too late: exclude station: ',trim(nsname)
                cycle
             endif
             if (.tanfdp.ms+(.nsamp.ms-1)*.dt.ms < ta+tsbuf+twinlen) then   ! Data end too early
                write(lulog,*) 'Data end too early: exclude station: ',trim(nsname)
                cycle
             endif
             ntraces = ntraces+1
             stnetnam(ntraces) = trim(nsname)
             useful_station(ntraces) = station
             azre(ntraces) = azr
             azim(ntraces) = azi
             tanf(ntraces) = ta
          else
             write(lulog,*) 'Data file for station '//trim(nsname)//' not found'
             cycle
          endif
       !
       !  create a time series object for observed data
       !
          call createDPLinkTimeSeries(draw,.nsamp.ms,.tanfdp.ms,.dt.ms,.trace.ms)
       !
       !  check if data are OK
       !
          if (isnan(energyUnnormalizedTimeSeries(draw))) then
             write(lulog,*) 'Raw data of '//trim(nsname)//' seem corrupt, exclude station'
             ntraces = ntraces-1
             goto 3
          endif
       !
       !  filter data, take filter specs from parfile_dyrt
       !
          fts0 = detrendTimeSeries(draw)
          fts1 = hanningTaperTimeSeries(fts0,0.05*lengthTimeSeries(draw))
          fts2 = lowPassButterworthRecursiveTimeSeries(fts1,evfilspecs(4),int(evfilspecs(3)))
          fts3 = highPassButterworthRecursiveTimeSeries(fts2,evfilspecs(2),int(evfilspecs(1)))
       !
       !  interpolate filtered data to new dt in time window ta, ta+tsbuf+twinlen
       !  skip time series on error
       !  full availability of data in time window was checked before
       !
          nsyn = nint((tsbuf+twinlen)/dtdata)
          fts4 = interpolateTimeSeries(fts3,nsyn,ta,dtdata,errmsg1)
          if (.level.errmsg1 == 2) then
             call reduceLevelErrorMessage(errmsg1)
             call print(errmsg1)
             write(lulog,*) 'Cannot interpolate trace: exclude '//trim(nsname)
             ntraces = ntraces-1
             goto 3
          endif
       !
       !  apply a taper of taplen to the end of the interpolated data
       !
          fts6 = hanningTaperTailTimeSeries(fts4,real(taplen))
       !
       !  compute energy of signal after onset and before onset and compare
       !  ta and .tanfdp.fts6 should be equal
       !
           n = 0.9*int((ta+tsbuf-.tanfdp.fts6)/dtdata)
           enoise = energyUnnormalizedTimeSeries(fts6,n2 = n)/n
           esig = energyUnnormalizedTimeSeries(fts6,n1 = n+1)/(.nsamp.fts6-n)
           snr(ntraces) = sqrt(esig/enoise)
           rmsnoise(ntraces) = sqrt(enoise)
           if (snr(ntraces) < snrbound) then
              write(lulog,*) 'SNR too low, exclude station: ',trim(nsname),' SNR = ',snr(ntraces)
              ntraces = ntraces-1
              goto 3
           endif
       !
       !  deep copy interpolated and tapered data to data traces array
       !
          call copyTimeSeries(fts6,dat(ntraces))
          rmsdat(ntraces) = rmsTimeSeries(fts6)
       !
       !  Fourier transform of interpolated and windowed data
       !
          call createFromTimeSeriesFourierSpectrum(ftdat(ntraces),dat(ntraces),ns=nfft)
       !
       !  downsample untruncated data and write to file, not used for other purposes
       !
          ndec = nint(dtdata/.dt.ms)
          fts5 = downsampleTimeseries(fts3,ndec)
          call copyTimeSeries(fts5,untruncated_dat(ntraces))
          rmssyn(ntraces) = hypot(azr,azi)
       !
 3        call dealloc(fts0); call dealloc(fts1); call dealloc(fts2); call dealloc(fts3)
          call dealloc(fts4); call dealloc(fts5); call dealloc(fts6)
          call dealloc(draw); call dealloc(ms)
       enddo                                      ! End station loop
       write(lulog,'(a,i6)') "Number of data traces: ",ntraces
       write(lulog,'(a,f15.5)') "DF = ",.df.ftdat(1)
    !
    !  write rms to ascii file
    !
       open(1,file = trim(path_measured_data)//trim(.evid.event)//'/rms.txt',iostat = ierr)
       if (ierr /= 0) then; errstr = 'cannot open rms.txt file'; goto 1; endif
       do j = 1,ntraces
          write(1,'(a,e15.3,e15.3)') stnetnam(j),rmsdat(j),rmssyn(j)
       enddo
       close(1)
    !
    !  take average of data energy of CH stations except TORNY as a reference
    !  this is only used for checking response problems
    !
       write(lulog,'(a)') "Swiss stations for checking response problems"
       rmsrefdat = 0.d0
       n = 0                ! count available Swiss stations
       do j = 1,ntraces
          if (stnetnam(j)(1:2) == 'CH' .and. trim(stnetnam(j)) /= 'CH.TORNY') then
             rmsrefdat = rmsrefdat+rmsdat(j)
             write(lulog,'(a,e15.3)') stnetnam(j),rmsdat(j)
             n = n+1
          endif
       enddo
       if (n == 0) then
          call add(errmsg1,2,'no CH station available above SNR bound',myname)
          call print(errmsg1)
          goto 2
       endif
       rmsrefdat = rmsrefdat/n
       write(lulog,'(a,2e15.3)') 'RMS-average of CH stations: ',rmsrefdat
    !
    !  look for response problems by forming ratio of rms of stations and reference rmsrefdat
    !  keep station if ratio is close to a power of ten from -12 to +12 in steps of 3
    !  correct response by scaling station data by power of ten
    !  if not write station name to specific file for later exclusion
    !
       open(1,file = trim(path_measured_data)//trim(.evid.event)//'/response_problem_stations.txt',iostat = ierr)
       if (ierr /= 0) then; errstr = 'cannot open response_problem_stations.txt file'; goto 1; endif
       eps = 0.3                ! accept deviations by a factor of 2 (= 10**0.3)
       nrp = 0                  ! count stations with response problems
       write(lulog,'(a)') "Log(rmsdat/rmsrefdat) based on reference CH data:"
       do j = 1,ntraces
          scaleflag = .false.
          frms = rmsdat(j)/rmsrefdat
          do n = -12,12,3
             if (abs(log10(frms)+n*1.d0) < eps) then
                call scalarMultiplyTimeSeries(dat(j),10.d0**n)
                call scalarMultiplyTimeSeries(untruncated_dat(j),10.d0**n)
                call scalarMultiplyFourierSpectrum(ftdat(j),10.d0**n)
                rmsdat(j) = rmsdat(j)*10.d0**n
                rmsnoise(j) = rmsnoise(j)*10.d0**n
                if (n /= 0) write(lulog,'(a,f12.3,a,i4)') trim(stnetnam(j)),log10(frms),' corrected with log factor: ',n
                scaleflag = .true.
                exit
             endif
          enddo
          if (.not. scaleflag) then
             write(1,'(a)') stnetnam(j)
             write(lulog,'(a,f12.3)') trim(stnetnam(j)),log10(frms)
             nrp = nrp+1
             idx(nrp) = j               ! store index of station with response problem
          endif
       enddo
       close(1)
    !
    !  sort out stations with problems
    !
       nrun = 0
       lastrun = .false.
 4     nrun = nrun+1
       write(crun,'(i2.2)') nrun
       write(lulog,'(a)') '---------------------------------------------------------------------'
       write(lulog,'(a,i6)') '                   R U N     ',nrun
       write(6,'(a,i6)') '                   R U N     ',nrun
       write(lulog,'(a)') '---------------------------------------------------------------------'
       do j = 1,nrp
          call dealloc(dat(idx(j)))                ! deallocate problem trace
          call dealloc(ftdat(idx(j)))                ! deallocate problem trace
          call dealloc(untruncated_dat(idx(j)))
          do n = idx(j),ntraces-1
             dat(n) = dat(n+1)
             ftdat(n) = ftdat(n+1)
             untruncated_dat(n) = untruncated_dat(n+1)
             stnetnam(n) = stnetnam(n+1)
             rmsdat(n) = rmsdat(n+1)
             rmssyn(n) = rmssyn(n+1)
             snr(n) = snr(n+1)
             rmsnoise(n) = rmsnoise(n+1)
             useful_station(n) = useful_station(n+1)
             azre(n) = azre(n+1)
             azim(n) = azim(n+1)
             tanf(n) = tanf(n+1)
          enddo
          ntraces = ntraces-1
          idx = idx-1                             ! correct indices because a bad one was taken out
       enddo
       write(lulog,'(a,i6)') "Number of available data traces: ",ntraces
    !----------------------------------------------------------------------------------------
    !  Compute Fourier transform of source wavelet
    !
       call sourceWavelet(ntraces,nsyn,dtdata,tanf,tsbuf,ftdat,azre,azim,&
                          ftstf,stf,htstf,chitr,chi2,conv,rmsdaconv,errmsg)
       if (.level.errmsg1 == 2) goto 20
    ! 
       write(lulog,'(a,f12.1)') 'Total misfit/Total energy in percent: ',chi2*100.d0
       write(6,'(a,f12.1)') 'Total misfit/Total energy in percent: ',chi2*100.d0
    !
    !  write trace fits and snr and rms ratio to file
    !
       if (.not. lastrun) then
          open(1,file = trim(path_measured_data)//trim(.evid.event)//'/run_'//crun//'_trace_misfits_snr.txt',iostat = ierr)
       else
          open(1,file = trim(path_measured_data)//trim(.evid.event)//'/trace_misfits_snr.txt',iostat = ierr)
       endif
       if (ierr /= 0) then; errstr = 'cannot open trace_misfit_snr.txt file'; goto 1; endif
       do j = 1,ntraces
          write(1,'(a,f12.1,2f12.3)') stnetnam(j),chitr(j)*100.d0,snr(j),rmsdaconv(j)
       enddo
       close(1)
    !
    !  write stations with too high misfit to log file and eliminate them
    !
       nrp = 0                ! count stations with too high misfit
       idx = 0                ! reset index array to zeros
       write(lulog,'(a)') 'Stations with too high misfit'
       do j = 1,ntraces
          if (chitr(j) > misfitbound) then
             write(lulog,'(a,e15.6)') stnetnam(j),chitr(j)*100.d0
             nrp = nrp+1
             idx(nrp) = j
          endif
       enddo
       close(1)
       if (nrp > 0 .and. .not. lastrun) then
          do j = 1,ntraces
             call dealloc(conv(j))
          enddo
          deallocate(chitr,conv,rmsdaconv)
          lastrun = .false.
          goto 4
       endif
    !
    !  write stations with too high or low rmsdat/rmsconv to log file and eliminate
    !
       nrp = 0                ! count stations with too high amplitude deviations
       idx = 0                ! reset index array to zeros
       write(lulog,'(a)') 'Stations with too high amplitude deviation'
       do j = 1,ntraces
          if (abs(log10(rmsdaconv(j))) > rmsbound) then
             write(lulog,'(a,f12.3)') stnetnam(j),rmsdaconv(j)
             nrp = nrp+1
             idx(nrp) = j
          endif
       enddo
       close(1)
       if (nrp > 0 .and. .not. lastrun) then
          do j = 1,ntraces
             call dealloc(conv(j))
          enddo
          deallocate(chitr,conv,rmsdaconv)
          lastrun = .false.
          goto 4
       endif
    !
    !  collect untruncated traces with misfit above visual check limit
    !
       nrp = 0
       idx = 0
       do j = 1,ntraces
          if (chitr(j) > visual_check_bound) then
             nrp = nrp+1
             idx(nrp) = j
          end if
       end do
       if (nrp > 0) then
          allocate(avcb_dat(nrp),avcb_conv(nrp),avcbnam(nrp))
          do j = 1,nrp
             call copyTimeSeries(untruncated_dat(idx(j)),avcb_dat(j))
             call copyTimeSeries(conv(idx(j)),avcb_conv(j))
             avcbnam(j) = stnetnam(idx(j))
          end do
       endif
       if (nrp > 0) then
          call writeAsciiGather(1,trim(path_measured_data)//trim(.evid.event)//'/above_visual_check_bound_data_UP',&
                                avcb_dat,avcbnam,nrp,errmsg1)
          if (.level.errmsg1 == 2) goto 20
          call writeAsciiGather(1,trim(path_measured_data)//trim(.evid.event)//'/above_visual_check_bound_conv_UP',&
                                avcb_conv,avcbnam,nrp,errmsg1)
          if (.level.errmsg1 == 2) goto 20
       endif
       do j = 1,nrp
          call dealloc(avcb_dat(j))
          call dealloc(avcb_conv(j))
       end do
       if (nrp > 0) deallocate(avcb_dat,avcb_conv,avcbnam)
    !
    !  write Fourier spectrum of source wavelet to file
    !
       open(1,file = trim(path_measured_data)//trim(.evid.event)//'/ftstf.txt',iostat = ierr)
       write(1,'(2i8,f16.8)') .nfa.ftstf,.nfb.ftstf,.df.ftstf
       do j = 1,.nf.ftstf
          write(1,'(2e20.8)') real(ftstf%y(j)),aimag(ftstf%y(j))
       enddo
       close(1)
    !
    !  write stf to file
    !
       tracefile = trim(path_measured_data)//trim(.evid.event)//'/stf.txt'
       call writeAsciiSynseisIO(1,tracefile,.nsamp.stf,.tanfdp.stf,real(.dt.stf),&
                                'Source_Wavelet','X',.trace.stf,ierr,.true.,.true.)
       if (ierr /= 0) then; errstr = 'cannot write stf file'; goto 1; endif
    !
    !  write Hilbert transform of stf to file
    !
       tracefile = trim(path_measured_data)//trim(.evid.event)//'/htstf.txt'
       call writeAsciiSynseisIO(1,tracefile,.nsamp.htstf,.tanfdp.htstf,real(.dt.htstf),&
                                'HT-Source_Wavelet','X',.trace.htstf,ierr,.true.,.true.)
       if (ierr /= 0) then; errstr = 'cannot write htstf file'; goto 1; endif
    !
    !  write final cleaned station file for FWI (ASKI station list)
    !
       open(1,file = trim(path_measured_data)//trim(.evid.event)//'/ASKI_clean_station_file',iostat = ierr)
       if (ierr /= 0) then; errstr = 'cannot open ASKI_clean_station_file'; goto 1; endif
       write(1,'(a)') .csys.useful_station(1)
       do j = 1,ntraces
          write(1,'(a,a5,2f15.5,f12.3)') .staname.useful_station(j),.netcode.useful_station(j),&
                                         .lat.useful_station(j),.lon.useful_station(j),.alt.useful_station(j)
       enddo
       close(1)
    !
    !  write a block for the data model space info file to be built later
    !
       call writeDataModelSpaceInfoBlock(1,trim(path_measured_data)//trim(.evid.event)//'/ASKI_dmspace_block',&
                                         ntraces,.evid.event,stnetnam,rmsnoise,errmsg1)
       if (.level.errmsg1 == 2) goto 20
    !
    !  write windowed and filtered data and convolved synthetics to event subfolder (one file per trace)
    !
       call system('mkdir -p '//trim(path_measured_data)//trim(.evid.event)//'/per_trace',ierr)
       call writeAsciiPerTrace(1,trim(path_measured_data)//trim(.evid.event)//'/per_trace/data_',dat,stnetnam,ntraces,errmsg1)
       if (.level.errmsg1 == 2) goto 20
       call writeAsciiPerTrace(1,trim(path_measured_data)//trim(.evid.event)//'/per_trace/conv_',conv,stnetnam,ntraces,errmsg1)
       if (.level.errmsg1 == 2) goto 20
    !
    !  write ensemble of data, untruncated_data and convolved synthetics to one gather file
    !
       call writeAsciiGather(1,trim(path_measured_data)//trim(.evid.event)//'/gather_data_UP',dat,stnetnam,ntraces,errmsg1)
       if (.level.errmsg1 == 2) goto 20
       call writeAsciiGather(1,trim(path_measured_data)//trim(.evid.event)//'/gather_untrunc_data_UP',&
                             untruncated_dat,stnetnam,ntraces,errmsg1)
       if (.level.errmsg1 == 2) goto 20
       call writeAsciiGather(1,trim(path_measured_data)//trim(.evid.event)//'/gather_conv_UP',conv,stnetnam,ntraces,errmsg1)
       if (.level.errmsg1 == 2) goto 20
    !
    !  clean up
    !
 2     do j = 1,ntraces
          call dealloc(dat(j))
          call dealloc(untruncated_dat(j))
          call dealloc(ftdat(j))
       enddo
       call dealloc(stf); call dealloc(htstf); call dealloc(ftstf)
       if (associated(stexclude)) deallocate(stexclude)
       deallocate(dat,untruncated_dat,ftdat,stnetnam,snr,useful_station,rmsdat,rmssyn,rmsnoise,rmsdaconv,tanf,azre,azim,idx)
       if (associated(conv)) then
          do j = 1,ntraces
             call dealloc(conv(j))
          enddo
          deallocate(conv)
       endif
       if (allocated(chitr)) deallocate(chitr)
       call dealloc(errmsg1)
    !
    !  close logfile
    !
       close(lulog)
       if (single_event .and. count_event > 0) exit
    enddo
    if (single_event .and. count_event == 0) then
       call add(errmsg,2,'Selected event is not in event list',myname)
       goto 10
    else
       print *,'Worked on ',count_event,' events'
    endif
    call closeFileHDFWrapper(fid,errmsg)
    call dealloc(errmsg)
    call dealloc(inpar)
!-----------------------------------------------------------------------------------
!  treat errors
!
 1  if (ierr < 0) then
       call add(errmsg,2,trim(errstr),myname)
       call print(errmsg)
       stop
    endif
 10 if (.level.errmsg == 2) then
       call print(errmsg)
    endif
 20 if (.level.errmsg1 == 2) then
       call print(errmsg1)
    endif
!
 contains
!----------------------------------------------------------------------------------------------
!  compute Fourier transform of source wavelet, source wavelet itself and its Hilbert transform
!  Processed data and their FT-transforms have all same length
!  Returns fourier_spectrum and time_series
!
    subroutine sourceWavelet(ntraces,nsyn,dtdata,tanf,tsbuf,ftdat,azre,azim,&
                             ftstf,stf,htstf,chitr,chi2,conv,rmsdaconv,errmsg)
    integer :: ntraces,nsyn
    double precision :: dtdata
    type (fourier_spectrum), dimension(:) :: ftdat
    double precision, dimension(:) :: azre,azim,tanf
    type (fourier_spectrum) :: ftstf
    type (time_series) :: stf,htstf
    double precision :: chi2,tsbuf
    double precision, dimension(:), allocatable :: chitr,rmsdaconv
    type (time_series), dimension(:), pointer :: conv
    type (error_message) :: errmsg
    double complex, dimension(:), allocatable :: hf,sp
    real, dimension(:), allocatable :: ts
    integer :: nf,nf1,nf2,jf,k,nsamp
    double precision :: f,df,sumaz2,om,endata,sum_endata,ensyn
    double complex :: csum,syn,dat
!
    nf1 = .nfa.ftdat(1)
    nf2 = .nfb.ftdat(1)
    df = .df.ftdat(1)
    nf = .nf.ftdat(1)
    sumaz2 = sum(azre(1:ntraces)**2)+sum(azim(1:ntraces)**2)
    allocate(hf(nf),sp(nf2))
    do jf = nf1,nf2
       csum = (0.d0,0.d0)
       om = 2.d0*mc_pid*(jf-1)*df
       do k = 1,ntraces
          csum = csum + dcmplx(azre(k),-azim(k))*ftdat(k)%y(jf-nf1+1)*exp(mc_cid*om*tsbuf)
       end do
       hf(jf-nf1+1) = csum/sumaz2
       sp(jf) = csum/sumaz2
    end do
 !
 !  Misfits
 !
    allocate(chitr(ntraces),rmsdaconv(ntraces))
    chi2 = 0.d0
    sum_endata = 0.d0
    do k=1,ntraces
       chitr(k) = 0.d0
       endata = 0.d0
       ensyn = 0.d0
       do jf = nf1,nf2
          f = (jf-1)*df
          om = 2.d0*mc_pid*f
          syn = dcmplx(azre(k),azim(k))*hf(jf-nf1+1)*exp(-mc_cid*om*tsbuf)
          dat = ftdat(k)%y(jf-nf1+1)
          chitr(k) = chitr(k) + abs(dat-syn)**2
          endata = endata + abs(dat)**2
          ensyn = ensyn + abs(syn)**2
       end do
       sum_endata = sum_endata+endata
       chitr(k) = chitr(k)/endata
       rmsdaconv(k) = sqrt(endata/ensyn)
       chi2 = chi2+chitr(k)*endata
    end do
    chi2 = chi2/sum_endata
 !
 !  create Fourier spectrum object
 !
    call createFromDataFourierSpectrum(ftstf,nf1,nf2,df,hf)
 !
 !  shift time zero of hf to -tsbuf to see if STF is nonzero at negative times
 !   
    do jf = nf1,nf2
       om = 2.d0*mc_pid*(jf-1)*df
       sp(jf) = sp(jf)*exp(-mc_cid*om*tsbuf)
    enddo
 !
 !  backtransform of hf to stf
 !
    nsamp = 2*nf2
    allocate(ts(nsamp))
    call transformFrequencyTime(nf1,nf2,df,0.d0,sp,nsamp,dtdata,ts,errmsg)
    if (.level.errmsg == 2) return
    call createDPFromDataTimeSeries(stf,nsamp,-tsbuf,dtdata,ts)
 !
 !  backtransform to Hilbert transform of stf = htstf
 !
    call transformFrequencyTime(nf1,nf2,df,0.d0,-mc_cid*sp,nsamp,dtdata,ts,errmsg)
    if (.level.errmsg == 2) return
    call createDPFromDataTimeSeries(htstf,nsamp,-tsbuf,dtdata,ts)
 !
 !  prediction for data (reuse sp-array)
 !
    allocate(conv(ntraces))
    do k = 1,ntraces
       do jf = nf1,nf2
          om = 2.d0*mc_pid*(jf-1)*df
          sp(jf) = dcmplx(azre(k),azim(k))*hf(jf-nf1+1)*exp(-mc_cid*om*tsbuf)
       enddo
       call transformFrequencyTime(nf1,nf2,df,0.d0,sp,nsamp,dtdata,ts,errmsg)
       if (.level.errmsg == 2) return
       call createDPFromDataTimeSeries(conv(k),nsyn,tanf(k),dtdata,ts)
    enddo
 !
    deallocate(hf,sp,ts) 
 !
    end subroutine sourceWavelet
!----------------------------------------------------------------------------------
!  routine to read file with stations to be excluded
!
    subroutine readExcludeFile(filename,p,ierr)
    character (len=*) :: filename
    character (len=char_len_netcode+char_len_sta+1), pointer, dimension(:) :: p
    character (len=max_length_string) :: line
    integer :: ierr,n
    open(1,file = trim(filename),status = 'old',action = 'read',iostat = ierr)
    if (ierr /= 0 ) return
    n = 0
    allocate(p(20))
    do while(.true.)
  1     read(1,'(a)',end = 11) line
       n = n+1
       p(n) = trim(line)
       if (n == size(p)) then
          p => reallocate(p,size(p)+20)
       endif
    enddo
 11 p => reallocate(p,n)
    close(1)
    end subroutine readExcludeFile
!-------------------------------------------------------------------------------------
!  check whether stnetname is in exclude list
!
    function checkExcludeStations(name,p) result(found)
    character (len=*) :: name
    character (len=char_len_netcode+char_len_sta+1), dimension(:), pointer :: p
    logical :: found
    integer :: j
    j = 0
    found = .false.
    if (.not. associated(p)) return
    do while (j < size(p))
       j = j+1
       if (name.equal.p(j)) then
          found = .true.
          exit
       endif
    enddo
    end function checkExcludeStations
!-------------------------------------------------------------------------------------
!  write a collection of time series by AsciiSynseisIO to individual files
!
    subroutine writeAsciiPerTrace(lu,basename,ts,stnetnam,ntraces,errmsg)
    integer :: lu
    character (len=*) :: basename
    type (time_series), dimension(:), pointer :: ts
    character (len=*), dimension(:) :: stnetnam
    integer :: ntraces
    type (error_message) :: errmsg
    integer :: j,ierr
    character (len=max_length_string) :: tracefile
!
    do j = 1,ntraces
       tracefile = trim(basename)//stnetnam(j)+'_UP'
       call writeAsciiSynseisIO(lu,trim(tracefile),.nsamp.ts(j),.tanfdp.ts(j),real(.dt.ts(j)),&
                                stnetnam(j),'UP',.trace.ts(j),ierr,.true.,.true.)
       if (ierr /= 0) then
          call add(errmsg,2,'cannot write trace','writeAsciiPerTrace')
          return
       endif
    enddo
    end subroutine writeAsciiPerTrace
!-------------------------------------------------------------------------------------
!  write a collection of time series by AsciiSynseisIO to one gather file
!
    subroutine writeAsciiGather(lu,filename,ts,stnetnam,ntraces,errmsg)
    integer :: lu
    character (len=*) :: filename
    type (time_series), dimension(:), pointer :: ts
    character (len=*), dimension(:) :: stnetnam
    integer :: ntraces
    type (error_message) :: errmsg
    integer :: j,ierr
    logical :: do_open
!
    do_open = .true.
    do j = 1,ntraces
       call writeAsciiSynseisIO(lu,filename,.nsamp.ts(j),.tanfdp.ts(j),real(.dt.ts(j)),&
                                stnetnam(j),'UP',.trace.ts(j),ierr,do_open,.false.)
       if (ierr /= 0) then
          call add(errmsg,2,'cannot write gather','writeAsciiGather')
          return
       endif
       do_open=.false.
    enddo
    close(lu)
    end subroutine writeAsciiGather
!-----------------------------------------------------------------------------
!  write a data model space info file for specific event (only data section)
!  used for transformation into freqeuncy domain
!
    subroutine writeDataModelSpaceInfoBlock(lu,filename,ntraces,evid,stnetnam,rmsnoise,errmsg)
    integer :: lu
    character (len=*) :: filename
    integer :: ntraces
    character (len=*) :: evid
    character (len=*), dimension(:) :: stnetnam
    double precision, dimension(:) :: rmsnoise
    type (error_message) :: errmsg
!
    open(lu,file = trim(filename),iostat = ierr)
    if (ierr /= 0) then
       call add(errmsg1,2,'cannot open data_model_space_info block file',myname)
       return
    endif
    write(lu,'(i6)') ntraces
    do j = 1,ntraces
        write(lu,'(a,1x,a,d15.3)') evid,stnetnam(j),1.d0/rmsnoise(j)
    enddo
    close(lu)
    end subroutine writeDataModelSpaceInfoBlock
!
end program computeDyrtSourceWavelet
