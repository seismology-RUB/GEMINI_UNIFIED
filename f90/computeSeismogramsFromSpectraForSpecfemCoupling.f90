! ===============================================================================
!  Main program for the computation of synthetic seismograms from precomputed
!  frequency spectra for SPECFEM coupling
! ===============================================================================
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
!  Compute synthetic seismograms from precomputed displacement and traction spectra
!  at boundary points of SPECFEM's computational box for a selected event.
!  User provides event ID, file with displacement and traction spectra and filter
!  specifications for the event. Seismograms (displacements and normal stresses)
!  are written to HDF for each boundary point.
!-----------------------------------------------------------------------------
program computeSeismogramsFromSpectraForSpecfemCoupling
    use hdf5
    use mathConstants
    use argumentParser
    use frequencyTime
    use eventFilter
    use synseisHDF
    use errorMessage
    use inputParameter
    use asciiSynseisIO
    use string
    use mpiSupport
    implicit none
    type (mpi_support) :: mpisup
    type (argument_parser) :: ap
    type (error_message) :: errmsg,errmsg1
    type (input_parameter) :: inpar
    type (seismic_event) :: event
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    real, dimension(:), pointer :: d,coor
    real, dimension(:,:), pointer :: uspre,uspim
    real, dimension(:,:), allocatable :: seis
    real, dimension(:), allocatable :: urs,stf
    real :: dtstf
    double complex :: zom
    double complex, dimension(:), allocatable :: evfilter,tsfil,zsp
    double precision, dimension(:), pointer :: evfilspecs
    double precision :: tanfstf,dt,dtp,tred,tend,df,sigma,tzero,ttmin,tt,twinlen,tsbuf
    integer, dimension(:), pointer :: id
    integer(hid_t) :: fidseis,fidspec,evgrid_spec,evgrid_seis,xferprp,dset,dsp,dsetcoor,dspcoor,procid,procid_spec
    integer(hsize_t), dimension(3) :: countvt,offset,dimsall
    integer(hsize_t), dimension(2) :: dimsslab
    integer :: slen,nsamp,nf1,nf2,nproc,ibp,ic,ierr,jf,myrank,nbpside,nbpbot,ip,numtasks,nstrue,nstf
    logical :: convolve_stf
    character(len=max_length_string) :: parfile,synseisfile,evfiltype,specfile,syntype,specpath
    character(len=max_length_string) :: errstr,stfpath,dummy1,dummy2,stffile
    character(len=char_len_evid) :: evid,eventid
    character (len=6) :: cip
    character (len=47) :: myname = 'computeSeismogramsFromSpectraForSpecfemCoupling'
    character (len=80), dimension(3) :: par_keys
!
!  keywords for input parameters
!
    data par_keys/'PATH_SYNTHETICS','EVENT_FILTER_TYPE','EVENT_FILTER_SPECS'/
!--------------------------------------------------------------------------------------------------
    call init(ap,myname,'Compute synthetic seismograms from precomputed displacement spectra for SPECFEM coupling')
    call addPosarg(ap,'eventid','sval',' Event ID')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call addOption(ap,'-stf',.true.,'Convolve with source time function taken from event subfolder of given path. '//&
                   'Assume name of file is stf.txt','sval','None')
    call parse(ap)
    parfile = ap.sval.'parfile'
    eventid = ap.sval.'eventid'
    stfpath = ap.sval.'-stf'
    convolve_stf = ap.optset.'-stf'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
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
    specpath = inpar.sval.'PATH_SYNTHETICS'
    evfiltype = inpar.sval.'EVENT_FILTER_TYPE'
    evfilspecs => dvecp(inpar,'EVENT_FILTER_SPECS',4,ierr)
    if (ierr /= 0) then
       call add(errmsg,2,'Input for event filter specs cannot be read',myname)
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
!  open precomputed displacement spectra HDF file
!
    specfile = trim(specpath)//trim(eventid)//'_spec.hdf'
    call openFileParallelAccessHDFWrapper(trim(specfile),fidspec,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  check type
!
    call readStringAttributeHDFWrapper(fidspec,'SyntheticsType',syntype,slen,errmsg)
    if (.level.errmsg == 2) goto 10
    if (.not. (syntype(1:slen).equal.'spectra')) then
       errstr = 'Displacement spectra file has wrong type!'; goto 10
    endif
!
!  read SPECFEM nproc as attribute from spectra file
!
    call readArrayAttributeHDFWrapper(fidspec,"nproc",aria,errmsg)
    if (.level.errmsg == 2) goto 10
    id => aria%get1d(); nproc = id(1); deallocate(id)
!
!  read event info from specfile, opens event group, return event
!
    call openEventSynseisHDF(fidspec,evid,evgrid_spec,event,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read df,sigma,nf1,nf2 from specfile
!
    call readArrayAttributeHDFWrapper(fidspec,"df_sigma",arra,errmsg)
    if (.level.errmsg == 2) goto 10
    d => arra%get1d()
    df = d(1); sigma = d(2); deallocate(d)
    call readArrayAttributeHDFWrapper(fidspec,"nf1_nf2",aria,errmsg)
    if (.level.errmsg == 2) goto 10
    id => aria%get1d()
    nf1 = id(1); nf2 = id(2); deallocate(id)
!
!  time shift filter
!
    call readArrayAttributeHDFWrapper(fidspec,"timing",arra,errmsg)
    if (.level.errmsg == 2) goto 10
    d => arra%get1d()
    twinlen = d(1); tsbuf = d(2); ttmin = d(3); tred = d(4); tend = d(5); deallocate(d)
    allocate(tsfil(nf2))
    tsfil = 0.d0
    do jf = nf1,nf2
       tsfil(jf) = zexp(mc_cid*mc_two_pid*(jf-1)*df*tred)*exp(+sigma*tred)
    enddo
!
!  allocate space for seismograms, set dt, nsamp
!  and reduced number of samples nstrue
!
    dtp = 0.25d0/(2.d0*nf2*df)
    call newNsampDtFrequencyTime(df,dtp,nsamp,dt)
    allocate(urs(nsamp))
    allocate(zsp(nf2))
    nstrue = int((tend-tred)/dt)+1
    if (convolve_stf) then                        ! check that nstrue > nstf
       if (nstrue < nstf) then
          call add(errmsg,2,'Synthetics should be longer than source wavelet',myname)
          goto 10
       endif
    endif
    allocate(seis(6,nstrue))
!
!  Calculate event filter, only if seismograms are computed
!
    call computeEventFilter(evfiltype,evfilspecs,nf1,nf2,df,sigma,evfilter,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  open file for seismograms output, write type
!
    synseisfile = trim(specpath)//trim(eventid)//'_seis.hdf'
    call createFileParallelAccessHDFWrapper(trim(synseisfile),fidseis,errmsg)
    if (.level.errmsg == 2) goto 10
    call setXferprpIndependentHDFWrapper(xferprp,errmsg)
    if (.level.errmsg == 2) goto 10
    call writeStringAttributeHDFWrapper(fidseis,'SyntheticsType','seismograms',errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read source time function for this event
!
    if (convolve_stf) then
       stffile = trim(stfpath)//'/'//trim(eventid)//'/stf.txt'
       call readAsciiSynseisIO(2,trim(stffile),nstf,tanfstf,dtstf,dummy1,dummy2,stf,ierr,op = .true.,cl = .true.)
       if (ierr /= 0) then; errstr = 'cannot read source wavelet'; goto 1; endif
       if (abs(dtstf-dt) > 1.d-6) then
          call add(errmsg,2,'DT of source wavelet differs from that of synthetice',myname)
          goto 10
       endif
       if (tanfstf > 0.d0) then
          call add(errmsg,2,'expect zero or negative start time for source wavelet',myname)
          goto 10
       endif
    endif
!
!  create a group for event in HDF output and write event info
!
    call createEventSynseisHDF(event,fidseis,evgrid_seis,errmsg)
    if (.level.errmsg1 == 2) goto 10
!
!  write dt, seismogram length
!
    allocate(d(2))
    d = [real(dt),(nstrue-1)*real(dt)]; call arra%assoc1d(d)
    call writeArrayAttributeHDFWrapper(fidseis,"samplingIntervalLength",arra,errmsg)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
!
!  write timing (twinlen,tsbuf,ttmin,tred,tend)
!
    allocate(d(4))
    d = [real(twinlen),real(tsbuf),real(ttmin),real(tred),real(tend)]; call arra%assoc1d(d)
    call writeArrayAttributeHDFWrapper(fidseis,"timing",arra,errmsg)
    if (.level.errmsg == 2) goto 10
    call arra%deassoc(); deallocate(d)
!
!  Loop over boundary points per SPECFEM process
!
    do ip = 0,nproc-1
       write(cip,'(i6.6)') ip
       call h5gopen_f(fidspec,'proc_'+cip,procid_spec,ierr)
       if (ierr < 0) then; errstr = "h5gopen"; goto 1; endif
       call readArrayAttributeHDFWrapper(procid_spec,"numBoundaryPointsPerProc",aria,errmsg)
       if (.level.errmsg == 2) goto 10
       id => aria%get1d(); call aria%deassoc()
       nbpside = id(1); nbpbot = id(2); deallocate(id)
       if (nbpside+nbpbot == 0) then
          call h5gclose_f(procid_spec,ierr)
          if (ierr < 0) then; errstr = "h5gopen"; goto 1; endif
          cycle
       endif
   !
   !  configure output data set, all components, all boundary points, reduced time samples
   !  create proc-group, data space and dataset
   !
       call h5gcreate_f(fidseis,'proc_'+cip,procid,ierr)
       dimsall = [6,nbpside+nbpbot,nstrue]
       call h5screate_simple_f(3,dimsall,dsp,ierr)
       if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
       call h5dcreate_f(procid,"velocityTraction",H5T_NATIVE_REAL,dsp,dset,ierr)
       if (ierr < 0) then; print *,'h5dcreate_f veltrac'; goto 1; endif
    !
    !  output data set for coordinates
    !
       dimsall = [9,nbpside+nbpbot,1]
       call h5screate_simple_f(2,dimsall(1:2),dspcoor,ierr)
       if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
       call h5dcreate_f(procid,"coordinates",H5T_NATIVE_REAL,dspcoor,dsetcoor,ierr)
       if (ierr < 0) then; print *,'h5dcreate_f coor'; goto 1; endif
   !
   !  write number of boundary points to output file
   !  allocate id because it is a pointer
   !
       allocate(id(2)); id = [nbpside,nbpbot]; call aria%assoc1d(id)
       call writeArrayAttributeHDFWrapper(procid,"numBoundaryPointsPerProc",aria,errmsg)
       call aria%deassoc(); deallocate(id)
   !
   !  Loop over boundary points
   !
       do ibp = myrank+1,nbpside+nbpbot,numtasks
          call new(errmsg1,myname)
      !
      !  read out displacement and traction spectra for boundary point from specfile
      !
          offset = [0,ibp-1,0]
          dimsslab = [6,nf2]
          countvt = [6,1,nf2]
          call readArrayHDFWrapper(procid_spec,"velocityTractionReal",arra,&
               errmsg1,xferprp = xferprp,dimsslab = dimsslab,offset = offset,count = countvt)
          if (.level.errmsg1 == 2) goto 11
          uspre => arra%get2d(); call arra%deassoc()
          call readArrayHDFWrapper(procid_spec,"velocityTractionImag",arra,&
               errmsg1,xferprp = xferprp,dimsslab = dimsslab,offset = offset,count = countvt)
          if (.level.errmsg1 == 2) goto 11
          uspim => arra%get2d(); call arra%deassoc()
      !
      !  read coordinates from spectra file
      !
          call readArrayHDFWrapper(procid_spec,"coordinates",arra,errmsg1,&
               xferprp = xferprp,dimsslab = [9_hsize_t],offset = offset(1:2),count = [9_hsize_t,1_hsize_t])
          if (.level.errmsg1 == 2) goto 11
          d => arra%get1d(); call arra%deassoc()
          tt = d(9); tzero = tt-ttmin
      !
      !  calculate seismograms for all components
      !
          do ic = 1,3
             zsp = 0.d0
             do jf = nf1,nf2
                zom = dcmplx(mc_two_pid*(jf-1)*df,-sigma)
                zsp(jf) = dcmplx(uspre(ic,jf),uspim(ic,jf))*evfilter(jf)*tsfil(jf)*mc_cid*zom
             enddo
             call transformFrequencyTime(nf1,nf2,df,sigma,zsp,nsamp,dt,urs,errmsg1)
             if (.level.errmsg1 == 2) goto 11
             if (convolve_stf) then
                call convolveWithSourceWavelet(nstrue,nstf,tanfstf,dt,stf,urs)
             endif
             seis(ic,:) = urs(1:nstrue)
             do jf = nf1,nf2
                zsp(jf) = dcmplx(uspre(ic+3,jf),uspim(ic+3,jf))*evfilter(jf)*tsfil(jf)
             enddo
             call transformFrequencyTime(nf1,nf2,df,sigma,zsp,nsamp,dt,urs,errmsg1)
             if (.level.errmsg1 == 2) goto 11
             if (convolve_stf) then
                call convolveWithSourceWavelet(nstrue,nstf,tanfstf,dt,stf,urs)
             endif
             seis(ic+3,:) = urs(1:nstrue)
          enddo     ! ic
          seis(:,1:int(tzero/dt)+1) = 0.d0          ! zero samples before arrival
      !
      !  write seismograms to output file
      !
          offset = [0,ibp-1,0]
          countvt = [6,1,nstrue]
          call arra%assoc2d(seis)
          call writeArrayHDFWrapper(procid,'velocityTraction',arra,errmsg1,&
               ds = dset,offset = offset,count = countvt)
          if (.level.errmsg1 == 2) goto 11
          call arra%deassoc()
      !
      !  write coordinates to output file
      !
          allocate(coor(9))
          coor(1:8) = d(1:8); coor(9) = tt
          call arra%assoc1d(coor)
          call writeArrayHDFWrapper(procid,'coordinates',arra,errmsg1,&
               ds = dsetcoor,offset = offset(1:2),count = [9_hsize_t,1_hsize_t])
          if (.level.errmsg1 == 2) goto 11
          call arra%deassoc()
          call dealloc(errmsg1)
          deallocate(uspre,uspim,coor,d)
       enddo     ! ibp
       call h5sclose_f(dsp,ierr)
       call h5dclose_f(dset,ierr)
       call h5sclose_f(dspcoor,ierr)
       call h5dclose_f(dsetcoor,ierr)
       call h5gclose_f(procid,ierr)
       call h5gclose_f(procid_spec,ierr)
    enddo    ! ip
!
!  clean up
!
    call h5gclose_f(evgrid_spec,ierr)
    call h5gclose_f(evgrid_seis,ierr)
    call h5fclose_f(fidseis,ierr)
    call h5fclose_f(fidspec,ierr)
    call h5close_f(ierr)
    deallocate(urs,zsp,seis,tsfil)
    if (allocated(stf)) deallocate(stf)
    if (allocated(evfilter)) deallocate(evfilter)
    call dealloc(mpisup)
!
!  treat error conditions
!
1   if (ierr < 0) then
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
!--------------------------------------------------------------------------
contains
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
!
end program
