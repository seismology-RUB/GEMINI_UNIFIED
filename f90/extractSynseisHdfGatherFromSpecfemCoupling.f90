! ===============================================================================
!  Main program for extracting synseis hdf seismogram gathers from a
!  Specfem coupling seismogram hdf file
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
!  User provides GEMINI parfile and event ID. Seismograms (displacements and normal stresses)
!  are written to HDF for selected boundary points.
!-----------------------------------------------------------------------------
program extractSynseisHdfGatherFromSpecfemCoupling
   use hdf5
   use mathConstants
   use argumentParser
   use synseisHDF
   use errorMessage
   use inputParameter
   use string
   use heapSort
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (seismic_event) :: event
   type (seismic_station) :: station
   type (any_rank_real_array) :: arra
   type (any_rank_integer_array) :: aria
   type (input_parameter) :: inpar
   real, dimension(:), pointer :: d
   real, dimension(:,:), pointer :: coor
   real, dimension(:,:), allocatable :: xyzw,xyze,xyzs,xyzn
   real, dimension(:,:,:), allocatable :: ursw,urse,urss,ursn
   double precision :: dt,tlen,xpm,ypm,zpm,r,beta,phi,rearth
   double precision :: twinlen,tsbuf,ttmin,tred,tend
   double precision, dimension(:), pointer :: pm_spacing,boxc
   integer, dimension(:), allocatable :: idxw,idxe,idxs,idxn
   integer, dimension(:), pointer :: id,boxdim
   integer(hid_t) :: fidseis,evgrid,procid,fidin,evoutid
   integer :: slen,nsamp,ierr,ibp,nbpvert,nbpbot,ip,nproc,nvseis,NGLLZ,cntw,cnte,cnts,cntn,step
   character(len=max_length_string) :: injection_file,injection_path,seisfile,syntype,parfile,coorfile,statfile,errstr
   character (len=char_len_evid) :: eventid,evid
   character (len=6) :: cibp,cip
   character (len=42) :: myname = 'extractSynseisHdfGatherFromSpecfemCoupling'
   character (len=80), dimension(5) :: par_keys
!
!  keywords for input parameters
!
   data par_keys/'MESH_SPACING','EXTERNAL_NODES_REARTH','BOX_CENTER','BOX_DIMENSIONS','PATH_SYNTHETICS'/
!--------------------------------------------------------------------------------------------------
   call init(ap,myname,'Extract synseis hdf gather from specfem coupling seismograms for selected event')
   call addPosarg(ap,'parfile','sval',' Gemini parameter file')
   call addPosarg(ap,'eventid','sval',' Event ID')
   call addOption(ap,'-step',.true.,'stepping for seismogram selection','ival','1')
   call parse(ap)
   parfile = ap.sval.'parfile'
   evid = ap.sval.'eventid'
   step = ap.ival.'-step' 
   if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
   call dealloc(ap)
!-----------------------------------------------------------------------------
   call new(errmsg,myname)
!-------------------------------------------------------------
!  read input parameters from parfile
!
   call createKeywordsInputParameter(inpar,par_keys)
   call readSubroutineInputParameter(inpar,1,parfile,errmsg)
   if (.level.errmsg == 2) goto 10
   pm_spacing => dvecp(inpar,'MESH_SPACING',3,ierr)
   rearth = inpar.dval.'EXTERNAL_NODES_REARTH'
   boxc => dvecp(inpar,'BOX_CENTER',2,ierr)
   boxdim => ivecp(inpar,'BOX_DIMENSIONS',3,ierr)
   injection_path = inpar.sval.'PATH_SYNTHETICS'
   injection_file = trim(injection_path)//trim(evid)//'_seis.hdf'
   NGLLZ = 5
!----------------------------------------------------------------
!  open HDF environment
!
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 10
!
!  open specfem coupling seismogram HDF file
!
   call openFileRoHDFWrapper(injection_file,fidin,errmsg)
   if (.level.errmsg == 2) goto 10
!
!  check type
!
   call readStringAttributeHDFWrapper(fidin,'SyntheticsType',syntype,slen,errmsg)
   if (.level.errmsg == 2) goto 10
   if (.not. (syntype(1:slen).equal.'seismograms')) then
      errstr = 'specfem coupling file has wrong type!'; goto 10
   endif
!
!  read event info from specfem coupling file, opens event group, return event
!
   call openEventSynseisHDF(fidin,eventid,evgrid,event,errmsg)
   if (.level.errmsg == 2) goto 10
   print *,'Extracting injection seimsograms for event: ',trim(eventid)
   if (.not. equalString(eventid,evid)) then
      call add(errmsg,2,'Injection file name and event ID are inconsistent',myname)
      goto 10
   endif
!
!  read dt and seismogram length
!
   call readArrayAttributeHDFWrapper(fidin,"samplingIntervalLength",arra,errmsg)
   if (.level.errmsg == 2) goto 10
   d => arra%get1d()
   dt = d(1); tlen = d(2); deallocate(d)
   nsamp = nint(tlen/dt)+1
   print *,'DT = ',dt,' NSAMP = ',nsamp
!
!  read number of SPECFEM proceses
!
   call readArrayAttributeHDFWrapper(fidin,"nproc",aria,errmsg)
   if (.level.errmsg == 2) goto 10
   id => aria%get1d()
   nproc = id(1); deallocate(id)
!
!  read timing
!
   call readArrayAttributeHDFWrapper(fidin,"timing",arra,errmsg)
   if (.level.errmsg == 2) goto 10
   d => arra%get1d()
   twinlen = d(1); tsbuf = d(2); ttmin = d(3); tred = d(4); tend = d(5)
   deallocate(d)
!
!  First collect all points and seismograms on vertical profiles separately
!
   nvseis = boxdim(3)*NGLLZ
   allocate(ursw(nsamp,6,nvseis),urse(nsamp,6,nvseis),urss(nsamp,6,nvseis),ursn(nsamp,6,nvseis))
   allocate(xyzw(nvseis,8),xyze(nvseis,8),xyzs(nvseis,8),xyzn(nvseis,8))
   allocate(idxw(nvseis),idxe(nvseis),idxs(nvseis),idxn(nvseis))
   cntw = 0; cnte = 0; cnts = 0; cntn = 0
!
!  loop over SPECFEM processes subgroups
!
   do ip = 0,nproc-1
      write(cip,'(i6.6)') ip
      call h5gopen_f(fidin,'proc_'//cip,procid,ierr)
      if (ierr < 0) then; errstr = "h5gopen "; goto 1; endif
   !
   !  read number of existing boundary points from specfem coupling file
   !
      call readArrayAttributeHDFWrapper(procid,"numBoundaryPointsPerProc",aria,errmsg)
      if (.level.errmsg == 2) goto 10
      id => aria%get1d(); call aria%deassoc()
      nbpvert = id(1); nbpbot = id(2)
      deallocate(id)
      print *,'Process: ',ip,' NVERT = ',nbpvert,' NBOT = ',nbpbot
   !
   !  read coordinates from specfem coupling file
   !
      call readArrayHDFWrapper(procid,"coordinates",arra,errmsg)
      if (.level.errmsg == 2) goto 10
      coor => arra%get2d(); call arra%deassoc()
   !
   !  Loop over points on vertical boundaries of this process
   !
      do ibp = 1,nbpvert
      !
      !  pseudo mesh coordinates
      !
         r = coor(4,ibp)
         beta = asin(coor(2,ibp)/r)
         phi = asin(coor(1,ibp)/(r*cos(beta)))
         xpm = rearth*phi
         ypm = rearth*beta
         zpm = r-rearth
      !
      !  check if on WE boundary and in center element north of y=0
      !
         if (abs(1.d0-0.5*pm_spacing(2)/ypm) < 1.d-3) then
            if (coor(6,ibp) < boxc(2)) then                           ! W side boundary
               cntw = cntw+1
               xyzw(cntw,:) = [xpm,ypm,zpm,dble(coor(4:8,ibp))]
               call readInjectionSeismogram(ibp,nsamp,ursw(:,:,cntw),errmsg)
               print *,'West boundary: ',ip,ibp,xpm,ypm,zpm
            else if (coor(6,ibp) > boxc(2)) then                      ! E side boundary
               cnte = cnte+1
               xyze(cnte,:) = [xpm,ypm,zpm,dble(coor(4:8,ibp))]
               call readInjectionSeismogram(ibp,nsamp,urse(:,:,cnte),errmsg)
               print *,'East boundary: ',ip,ibp,xpm,ypm,zpm
            end if
         end if
      !
      !  check if on NS boundary and in center element east of x=0
      !
         if (abs(1.d0-0.5*pm_spacing(1)/xpm) < 1.d-3) then
            if (coor(5,ibp) < boxc(1)) then                            ! S side boundary
               cnts = cnts+1
               xyzs(cnts,:) = [xpm,ypm,zpm,dble(coor(4:8,ibp))]
               call readInjectionSeismogram(ibp,nsamp,urss(:,:,cnts),errmsg)
               print *,'South boundary: ',ip,ibp,xpm,ypm,zpm
            else if (coor(5,ibp) > boxc(1)) then                       ! N side boundary
               cntn = cntn+1
               xyzn(cntn,:) = [xpm,ypm,zpm,dble(coor(4:8,ibp))]
               call readInjectionSeismogram(ibp,nsamp,ursn(:,:,cntn),errmsg)
               print *,'North boundary: ',ip,ibp,xpm,ypm,zpm
            end if
         end if
      end do                           ! boundary points
      deallocate(coor)
      call h5gclose_f(procid,ierr)
   end do                              ! process group
!
!  Sort xyz according to z
!
   call heapSort2D(xyzw,idxw,1,3)
   call heapSort2D(xyze,idxe,1,3)
   call heapSort2D(xyzs,idxs,1,3)
   call heapSort2D(xyzn,idxn,1,3)
!
!  open synseis hdf gather file and create event
!
   seisfile = 'vbseis_'//trim(eventid)//'.hdf'
   call createFileHDFWrapper(seisfile,fidseis,errmsg)
   if (.level.errmsg == 2) goto 10
   call createEventSynseisHDF(event,fidseis,evoutid,errmsg)
   if (.level.errmsg == 2) goto 10
!
!  open output file for coordinates and travel time
!
   coorfile = 'vbcoor_'//trim(eventid)//'.txt'
   open(1,file=trim(coorfile),iostat = ierr)
   if (ierr /= 0) then; errstr = 'can not open coorfile '//trim(coorfile); goto 1; endif
!
!  open output station file
!
   statfile = 'vbstat_'//trim(eventid)//'.txt'
   open(2,file=trim(statfile),iostat = ierr)
   if (ierr /= 0) then; errstr = 'can not open statfile '//trim(statfile); goto 1; endif
   write(2,'(a)') 'S'
!
!  Loop through sorted points and write seismograms to file
!
   do ibp = step,cntw,step
      write(cibp,'(i3.3)') ibp
      call createSeismicStation(station,'C','W'//cibp,real(xyzw(ibp,1)),real(xyzw(ibp,2)),alt=real(xyzw(ibp,3)),netcode='VB')
      call writeStationSynseisHDF(station,fidseis,nsamp,tred,dt,'XYZPQR',ursw(:,:,idxw(ibp)),errmsg)
      if (.level.errmsg == 2) goto 10
      write(1,'(a,4f16.2,4f15.6)') 'W'//cibp,xyzw(ibp,:)
      write(2,'(a10,a6,2f15.5,f15.1)') 'W'//cibp,'VB',xyzw(ibp,5),xyzw(ibp,6),xyzw(ibp,3)
   end do
   do ibp = step,cnte,step
      write(cibp,'(i3.3)') ibp
      call createSeismicStation(station,'C','E'//cibp,real(xyze(ibp,1)),real(xyze(ibp,2)),alt=real(xyze(ibp,3)),netcode='VB')
      call writeStationSynseisHDF(station,fidseis,nsamp,tred,dt,'XYZPQR',urse(:,:,idxe(ibp)),errmsg)
      if (.level.errmsg == 2) goto 10
      write(1,'(a,4f16.2,4f15.6)') 'E'//cibp,xyze(ibp,:)
      write(2,'(a10,a6,2f15.5,f15.1)') 'E'//cibp,'VB',xyze(ibp,5),xyze(ibp,6),xyze(ibp,3)
   end do
   do ibp = step,cnts,step
      write(cibp,'(i3.3)') ibp
      call createSeismicStation(station,'C','S'//cibp,real(xyzs(ibp,1)),real(xyzs(ibp,2)),alt=real(xyzs(ibp,3)),netcode='VB')
      call writeStationSynseisHDF(station,fidseis,nsamp,tred,dt,'XYZPQR',urss(:,:,idxs(ibp)),errmsg)
      if (.level.errmsg == 2) goto 10
      write(1,'(a,4f16.2,4f15.6)') 'S'//cibp,xyzs(ibp,:)
      write(2,'(a10,a6,2f15.5,f15.1)') 'S'//cibp,'VB',xyzs(ibp,5),xyzs(ibp,6),xyzs(ibp,3)
   end do
   do ibp = step,cntn,step
      write(cibp,'(i3.3)') ibp
      call createSeismicStation(station,'C','N'//cibp,real(xyzn(ibp,1)),real(xyzn(ibp,2)),alt=real(xyzn(ibp,3)),netcode='VB')
      call writeStationSynseisHDF(station,fidseis,nsamp,tred,dt,'XYZPQR',ursn(:,:,idxn(ibp)),errmsg)
      if (.level.errmsg == 2) goto 10
      write(1,'(a,4f16.2,4f15.6)') 'N'//cibp,xyzn(ibp,:)
      write(2,'(a10,a6,2f15.5,f15.1)') 'N'//cibp,'VB',xyzn(ibp,5),xyzn(ibp,6),xyzn(ibp,3)
   end do

!  clean up
!
   deallocate(pm_spacing,boxc,boxdim,xyzw,xyze,xyzs,xyzn)
   deallocate(ursw,urse,urss,ursn)
   deallocate(idxw,idxe,idxs,idxn)
   call dealloc(errmsg)
   close(1)
   close(2)
   call h5gclose_f(evgrid,ierr)
   call h5gclose_f(evoutid,ierr)
   call h5fclose_f(fidseis,ierr)
   call h5fclose_f(fidin,ierr)
   call h5close_f(ierr)
!
!  treat error conditions
!
1  if (ierr /= 0) then
      print *,'ERROR: ',trim(errstr)
   endif
10 if (.level.errmsg == 2) then
      call print(errmsg)
   endif
!---------------------------------------------------------------------
   contains
! --------------------------------------------------------------------
!  read injection seismogram
!
   subroutine readInjectionSeismogram(ibp,nsamp,urs,errmsg)
      integer :: ibp,nsamp
      real, dimension(:,:) :: urs
      type (error_message) :: errmsg
      integer(hsize_t), dimension(3) :: countvt,offset
      integer(hsize_t), dimension(2) :: dimsslab
      type (any_rank_real_array) :: arra
      real, dimension(:,:), pointer :: seis
      integer :: ic
   !
      offset = [0,ibp-1,0]
      dimsslab = [6,nsamp]
      countvt = [6,1,nsamp]
      call readArrayHDFWrapper(procid,"velocityTraction",arra,errmsg,&
                               dimsslab = dimsslab,offset = offset,count = countvt)
      if (.level.errmsg == 2) return
      seis => arra%get2d(); call arra%deassoc()
      do ic = 1,6
         urs(:,ic) = seis(ic,:)
      enddo
      deallocate(seis)
   end subroutine readInjectionSeismogram
!
end program extractSynseisHdfGatherFromSpecfemCoupling
