! ===============================================================================
!  Main program for computing the kernel Green tensor on a grid
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
!   Compute Green tensor from receiver positions at prescribed grid points
!   for later calculation of sensitivity kernels
!-----------------------------------------------------------------------------
 program computeKernelGreenTensor
    use mathConstants
    use argumentParser
    use greenFKSpectra
    use singleFrequencyDisplacement
    use wavenumberIntegrals
    use readEventStationFile
    use chunkCubedSphere
    use geminiEarthModel
    use coordinateTransform
    use errorMessage
    use fileUnitHandler
    use asciiDataIO
    use inputParameter
    use mpiSupport
    use string
    implicit none
    integer, parameter :: nc = 3
    type (argument_parser) :: ap
    type (green_fk_spectra) :: gfk
    type (error_message) :: errmsg,errmsg2
    type (file_unit_handler) :: fuh
    type (input_parameter) :: inpar
    type (seismic_station) :: station
    type (seismic_network) :: station_list
    type (chunk_cubed_sphere) :: ccs
    type (gemini_earth_model) :: gem
    type (mpi_support) :: mpisup
    complex, dimension(nc) :: zsp,zspdr,zspdt,zspdf
    complex, dimension(nc,nc) :: gradus,gepicar
    complex, dimension(nc,nc,nc) :: gradgepicar,gradgcar
    complex, dimension(:), allocatable :: wnint
    complex, dimension(:,:), allocatable :: green_wn_spectra
    complex, dimension(:,:,:), allocatable :: gout,gstrain
    real, dimension(:), allocatable :: epidis,propdir_rad,phi_rad,beta
    real, dimension(:,:,:), allocatable :: besselj,tm
    real, dimension(:), pointer :: xr,yr
    real, dimension(:,:), pointer :: tmlsrc,tmrcgc,tmlsrcfo
    real, dimension(3,3) :: tmfo
    real, dimension(3) :: force
    real :: df,r,dum,latrad,lonrad
    real :: wp_lon,wp_lat,wp_clon,wp_clat,wp_rot
    integer, dimension(:), allocatable :: kint
    integer :: nf1,nf2,nf,jf,ic,jc,ig,ng,jrs,ios,luout,n,is
    integer :: isp,jsp,j,offset,nwint,jr,wp_nlon,wp_nlat
    integer, dimension(3,3) :: stridx = reshape((/ 1,6,5,6,2,4,5,4,3 /),shape = (/3,3/))  ! means stridx(2,3) = 4
    character (len=max_length_string) :: text,stationfile,dsvbasename,parfile,gemfile,outfilebase
    character (len=character_length_staname) :: staname
    character (len=nc) :: comps = 'ZRT'
    character (len=27) :: myname = 'computeKernelGreenTensor'
    character (len=80), dimension(11) :: gemini_par_keys
    character (len=3) :: cjrs,cjf
!
!  keywords for input parameters
!
    data gemini_par_keys/'GFDSV_BASENAME','FILE_STATION_LIST', &
          & 'GEMINI_EARTH_MODEL','OUTFILE_BASE', &
          & 'WAVEFIELD_POINTS_CLON','WAVEFIELD_POINTS_CLAT','WAVEFIELD_POINTS_WLON', &
          & 'WAVEFIELD_POINTS_WLAT','WAVEFIELD_POINTS_NLON','WAVEFIELD_POINTS_NLAT', &
          & 'WAVEFIELD_POINTS_ROT'/
!----------------------------------------------------------------------------------
    call init(ap,myname,'Compute kernel Green tensor using Green FK-spectra for ASKI')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call parse(ap)
    parfile = ap.sval.'parfile'
    if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); stop; endif
    call document(ap); call dealloc(ap)
!
    call new(fuh,20)
    call new(errmsg,myname)
!
!  create MPI support object
!
    call new(mpisup)
!-------------------------------------------------------------
!  read input parameters from parfile
!
    call createKeywordsInputParameter(inpar,gemini_par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
    dsvbasename = inpar.sval.'GFDSV_BASENAME'
    stationfile = inpar.sval.'FILE_STATION_LIST'
    gemfile = inpar.sval.'GEMINI_EARTH_MODEL'
    outfilebase = inpar.sval.'OUTFILE_BASE'
    wp_clon = inpar.rval.'WAVEFIELD_POINTS_CLON'
    wp_clat = inpar.rval.'WAVEFIELD_POINTS_CLAT'
    wp_lon = inpar.rval.'WAVEFIELD_POINTS_WLON'
    wp_lat = inpar.rval.'WAVEFIELD_POINTS_WLAT'
    wp_nlon = inpar.ival.'WAVEFIELD_POINTS_NLON'
    wp_nlat = inpar.ival.'WAVEFIELD_POINTS_NLAT'
    wp_rot  = inpar.rval.'WAVEFIELD_POINTS_ROT'
    call printInputParameter(inpar)
    call dealloc(inpar)
!----------------------------------------------------------------
!  read gemini earth model
!
    call readGeminiEarthModel(gem,1,gemfile,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); call abort(mpisup); endif
!----------------------------------------------------------------
!  read station file
!  use ASKI convention for event and station file
!
    call createStationListFromStationFile(stationfile,1,'ASKI_stations',station_list,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); call abort(mpisup); endif
!---------------------------------------------------------------------------
!  generate wavefield points on unit sphere
!
    call createChunkCubedSphere(ccs,wp_clon,wp_clat,wp_rot,wp_lon,wp_lat,wp_nlon,wp_nlat,errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); call abort(mpisup); endif
    ng = .ng.ccs
    xr => .th.ccs; yr => .ph.ccs
!---------------------------------------------------------------------------------
!  overwrite earth model file and append ng to it (needed by kernelReferenceModel)
!
    call writeGeminiEarthModel(gem,1,gemfile,ng,errmsg)
!------------------------------------------------------------------------
!  allocate reusable space
!
    allocate(epidis(ng),propdir_rad(ng),phi_rad(ng),beta(ng))
!
!  allocate displacement and strain arrays
!
    allocate(gout(ng,nc,3),gstrain(ng,6,3))           ! Output arrays
    allocate(tm(ng,nc,nc))                            ! Transformation matrix 
!-----------------------------------------------------------------------------
!  Loop over stations
!
    do while (nextStationSeismicNetwork(station_list,station,is))
       if (.csys.station == 'C') then
          print *,'This program is designed for spherical applications'
          print *,'Please provide event files with spherical coordinates!'
          call abort(mpisup)
       endif
       latrad = .latrad.station
       lonrad = .lonrad.station
       staname = .staname.station
    !
    !  transformation matrix between station cartesian basis and global cartesian basis
    !  tmlsrcfo needed for force basis vector ordering
    !
       tmrcgc => localToGlobalCartesianOnSphereCoordinateTransform(latrad,lonrad)
       tmlsrcfo => localSphericalToGlobalCartesianOnSphereCoordinateTransform(0.5*mc_pi,0.0)
       tmfo = matmul(tmrcgc,tmlsrcfo)
       deallocate(tmlsrcfo)
    !
    !  open frequency wavenumber spectra file for receiver (single force at surface)    
    !
       jrs = getNodeIndexFromDepthGeminiEarthModel(gem,0.0)
       write(cjrs,'(i3.3)') jrs
       text = dsvbasename+'.FORCE.'+cjrs
       call new(errmsg2,myname)
       call readMetaGreenFKSpectra(gfk,fuh,text,errmsg2)
       if (.level.errmsg2 == 2) then; call print(errmsg2); stop; endif
       if (gfk%istyp > 0) then
          call add(errmsg2,2,'You need FK-spectra for a force',myname)
          call print(errmsg2); call abort(mpisup)
       endif
       call dealloc(errmsg2)
    !
    !  set frequency index limits
    !
       nf1 = gfk%nf1; nf2 = gfk%nf2
       nf = nf2-nf1+1; df = gfk%df
    !
    !  write basic information on frequencies and grid to outfilebase+'_'+staname+'.meta'
    !
       luout = get(fuh)
       open(luout,file = outfilebase+'_'+staname+'.meta', iostat = ios, form = 'formatted')
       if (ios /= 0) then
          call new(errmsg2,2,'Could not open output file: '+outfilebase+'_'+staname+'.meta',myname)
          call print(errmsg2); call abort(mpisup)
       endif
       write(luout,'(a)') trim(staname)
       write(luout,*) df,nf1,nf,ng,gfk%nnod,gfk%sigma
       close(luout); call add(fuh,luout)
    !
    !  get array of wavenumver integral indices as required by dsvmask and istyp
    !   
       call getIndexArrayWavenumberIntegrals(gfk,nwint,kint)
       allocate(wnint(maxIndexWavenumberIntegrals(gfk%istyp)))
    !
    !  space for Bessel functions
    !
       allocate(besselj(gfk%nwnmax,0:2,ng))
    !
    !  Loop over grid points
    !
       do ig = 1,ng
       !
       !  get distances and angles to wavefield points (from south over east)
       !  as seen at the station location
       !  get propagation directions as seen at wavefield points (from S to E)
       !
          call geo2epi(xr(ig),yr(ig),0.5*mc_pi-latrad,lonrad,beta(ig),phi_rad(ig))
          call geo2epi(0.5*mc_pi-latrad,lonrad,xr(ig),yr(ig),dum,propdir_rad(ig))
          epidis(ig) = beta(ig)*gfk%rearth
       !
       !  precompute Bessel functions for epicentral distance in km
       !   
          call besselWavenumberIntegrals(gfk%nwnmax,gfk%dwn,beta(ig),besselj(:,0:2,ig))
       !
       !  Transformation matrix from local ZRT to global XYZ system
       !  tmlsrc: 
       !    matrix transforming basis vectors of epicentral ZRT system at wavefield point 
       !    to cartesian basis vectors of SEZ system at station
       !  tmrcgc:
       !    matrix transforming basis vectors of SEZ system at station 
       !    to basis vectors of global XYZ system 
       !  tmrcgc*tmlsrc: 
       !    total matrix transforming from spherical basis vectors of epicentral ZRT
       !    to cartesian basis vectors of global XYZ system
       !
          tmlsrc => localSphericalToGlobalCartesianOnSphereCoordinateTransform(0.5*mc_pi-beta(ig),phi_rad(ig))
          tm(ig,:,:) = matmul(tmrcgc,tmlsrc)
          deallocate(tmlsrc)
       enddo                       ! end grid point loop
    !
    !  start frequency loop to compute displacement vs frequency
    !  for all grid points and depth points
    !
       do jf = nf1,nf2
       !
       !  open outut file for kernel Green tensor and Green strains
       !  a separate file for each frequency
       !
          luout = get(fuh)
          write(cjf,'(i3.3)') jf-1     ! ASKI uses f = jf*df
          open(luout,file = outfilebase+'_'+staname+'.'+cjf, iostat = ios, form = 'unformatted')
          if (ios /= 0) then
             call new(errmsg2,2,'Could not open output file: '+outfilebase+'_'+staname+'.'+cjf,myname)
             call print(errmsg2); call abort(mpisup)
          endif
       !
       !  radial node loop
       !
          do jr = 1,size(gfk%rnod)
             call new(errmsg2,myname)
             call readDataNodeGreenFKSpectra(gfk,jf,jr,green_wn_spectra,errmsg2)
             if (.level.errmsg2 == 2) then; call print(errmsg2); call abort(mpisup); endif
             call dealloc(errmsg2)
          !
          !  Loop over grid points
          !
             do ig = 1,ng
             !
             !  first compute wavenumber integrals
             !
                do j = 1,nwint
                   if (gfk%istyp == 1) then
                      isp = moment_component_jump_from_wavenumber_integral(kint(j),1)
                      jsp = moment_component_jump_from_wavenumber_integral(kint(j),2)
                      offset = positionDataGreenFKSpectra(gfk,isp,jsp)
                      call momentWavenumberIntegrals(gfk,green_wn_spectra(:,offset),kint(j),epidis(ig),&
                                                   & besselj(:,0:,ig),0.2,wnint(kint(j)))
                   else
                      isp = force_component_jump_from_wavenumber_integral(kint(j),1)
                      jsp = force_component_jump_from_wavenumber_integral(kint(j),2)
                      offset = positionDataGreenFKSpectra(gfk,isp,jsp)
                      call forceWavenumberIntegrals(gfk,green_wn_spectra(:,offset),kint(j),epidis(ig),&
                                                  & besselj(:,0:,ig),0.2,wnint(kint(j)))
                   endif
                enddo
             !
             !  then compute components and derivatives of Green tensor
             !
                r = gfk%rnod(jr)
                do n = 1,3                        ! force loop
                   force = 0.; force(n) = 1.0
                   do ic = 1,nc
                      call forceSingleFrequencyDisplacement(comps(ic:ic),propdir_rad(ig), &
                           & phi_rad(ig),force,wnint,zsp(ic))
                      call forceRadDerivSingleFrequencyDisplacement(comps(ic:ic),phi_rad(ig), &
                           & force,wnint,zspdr(ic))
                      call forceThetaDerivsingleFrequencyDisplacement(comps(ic:ic),r,phi_rad(ig), &
                           & force,wnint,zspdt(ic))
                      call forcePhiDerivSingleFrequencyDisplacement(comps(ic:ic),r,beta(ig),phi_rad(ig), &
                           & force,wnint,zspdf(ic))
                   enddo
                !
                !  calculate spherical displacement gradient matrix
                !
                   call sphericalGradientSingleFrequencyDisplacement(zsp,zspdr,zspdt,zspdf,r,beta(ig),gradus)
                !
                !  conversion to global cartesian coordinates (displacement index and derivative)
                !
                   gepicar(:,n) = matmul(tm(ig,:,:),zsp(:))
                   gradgepicar(:,:,n) = matmul(tm(ig,:,:),matmul(gradus,transpose(tm(ig,:,:))))
                enddo                      ! end force loop
             !
             !  conversion of force index n to global cartesian:
             !  G_in = G_is R_sn and d_j G_in = d_j G_is R_sn
             !
                gout(ig,:,:) = matmul(gepicar,transpose(tmfo))    ! transpose correct here ??
                do jc = 1,3
                   gradgcar(jc,:,:) = matmul(gradgepicar(jc,:,:),transpose(tmfo))
                enddo
             !
             !  calculation of Green strains
             !
                do ic = 1,nc          ! component loop
                   do jc = ic,nc      ! derivative loop
                      gstrain(ig,stridx(ic,jc),:) = 0.5*(gradgcar(jc,ic,:)+gradgcar(ic,jc,:))
                   enddo
                enddo
             enddo                         ! end grid point loop
             deallocate(green_wn_spectra)
          !
          !  write Green tensor for this frequency and radial node to file
          !
             do n = 1,3
                do ic = 1,nc
                   write(luout) (gout(ig,ic,n),ig = 1,ng)
                enddo
                do ic = 1,6
                   write(luout) (gstrain(ig,ic,n),ig = 1,ng)
                enddo
             enddo
          enddo                              ! end radial node loop
          close(luout); call add(fuh,luout)
       enddo                                 ! end frequency loop
       deallocate(tmrcgc,besselj,wnint,kint)
       call dealloc(gfk,fuh)
    enddo                                    ! end event loop
!
!  deallocation
!
    call print(errmsg)
    deallocate(gout,gstrain)
    deallocate(epidis,propdir_rad,phi_rad)
    call dealloc(station_list)
    call dealloc(errmsg); call dealloc(fuh)
    call dealloc(gem); call dealloc(ccs)
!
 end program
