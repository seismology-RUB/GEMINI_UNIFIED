! ====================================================================================
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
!----------------------------------------------------------------------
!  This module provides routines for coupling Gemini to Specfem.
!  There are 4 boundaries. On each one, the GLL-points lie on a regular grid
!  spanned by distinct radii and distinct points at the surface. This implies
!  that there is only one set of lateral positions for all radii, and only
!  one set of radii for all lateral positions.
!  If the boundaries were exactly plane, each one could be characterized
!  by one single normal vector. This is approximately true but there are
!  sytematic deviations of about 4 percent on the northern and southern boudaries (Why?).
!  However, normal vectors are identical for fixed lateral position and varying radius
!  to high accuracy (1.e-5). Therefore, we assume here that there is only one set of normals
!  for all radii. This set has the same size and order as the set of lateral points.
!----------------------------------------------------------------------------------------
module boundaryPoints
       use hdf5
    use axesRotation
    use mathConstants
    use locatePoint
    use asciiDataIO
    use errorMessage
    use hdfWrapper
    use heapSort
    use string
    implicit none
    integer, parameter :: NGLL = 5                          ! number of GLL points in one dimension
    integer, parameter :: NGNOD = 4                         ! number of vertices in SPECFEM element
    integer, parameter :: NBOUND = 5                        ! number of boundaries (WESNB)
    integer, parameter :: NBPIDX = 7                        ! number of boundary point indices for use with iterator
    real, parameter :: REARTH = 6371.0                      ! earth radius in km
    interface dealloc
        module procedure deallocBoundaryPoints
    end interface dealloc
    type boundary_points
        integer :: dimx,dimy,dimz                           ! number of elements in computational box along coordinates
        real :: latcenter,loncenter                         ! latitude and longitude of center of computational box (degrees)
        double precision :: rtol,dtol                       ! radial and epidistance tolerances for defining unique points
        integer :: nbptot                                   ! total number of boundary points
        integer, dimension(:), allocatable :: nbpvert       ! number of vertical sides boundary points per process
        integer, dimension(:), allocatable :: nbpbot        ! number of bottom boundary points per process
        integer, dimension(:), allocatable :: sortidx       ! index array telling original position assembled array of vertical b-points 
        real, dimension(:,:), allocatable :: rlatlonnvside  ! radially sorted array of spherical coordinates of vertical b-points and cartesian normals        
        real, dimension(:,:), allocatable :: rlatlonnvbot   ! array of spherical coordinates of bottom b-points and cartesian normals
        real, dimension(:), pointer :: runique              ! radii of GLL grid (km) (excluding those occuring twice)
    end type boundary_points
    type single_boundary_point
        integer :: current                                  ! current point count
        double precision :: r,lat,lon                       ! coordinates of point (r,lat,lon or z,x,y)
        double precision :: nx,ny,nz                        ! normal vector at point
    end type single_boundary_point
contains
!-----------------------------------------------------------------------
!  Deallocate    
!
    subroutine deallocBoundaryPoints(this)
    type (boundary_points) :: this
    if (allocated(this%nbpvert)) deallocate(this%nbpvert)
    if (allocated(this%nbpbot)) deallocate(this%nbpbot)
    if (allocated(this%sortidx)) deallocate(this%sortidx)
    if (allocated(this%rlatlonnvside)) deallocate(this%rlatlonnvside)
    if (allocated(this%rlatlonnvbot)) deallocate(this%rlatlonnvbot)
    if (associated(this%runique)) deallocate(this%runique)
    end subroutine deallocBoundaryPoints
!---------------------------------------------------------------------------------
!  Read boundary data from HDF file
!
    subroutine readBoundaryPoints(this,hdffile,errmsg,parallel)
    type (boundary_points) :: this
    character (len=*) :: hdffile
    type (error_message) :: errmsg
    logical, optional :: parallel
    integer(hid_t) :: fid
    integer :: ierr
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
    real, dimension(:), pointer :: d
    integer, dimension(:), pointer :: id
    integer(hsize_t), dimension(:), allocatable :: dims
!
    if (present(parallel) .and. parallel) then
       call openFileParallelAccessHDFWrapper(trim(hdffile),fid,errmsg)
    else
       call openFileRoHDFWrapper(trim(hdffile),fid,errmsg)
    endif
    if (.level.errmsg == 2) return
!
    call readArrayAttributeHDFWrapper(fid,"centerLatLon",arra,errmsg)
    if (.level.errmsg == 2) return
    d => arra%get1d()
    this%latcenter = d(1); this%loncenter = d(2)
    call arra%dealloc()
!
    call readArrayAttributeHDFWrapper(fid,"tolerances",arra,errmsg)
    if (.level.errmsg == 2) return
    d => arra%get1d()
    !this%rtol = d(1); this%dtol = 0
    this%rtol = d(1); this%dtol = d(2)*mc_deg2radd
    !this%rtol = d(1); this%dtol = 0.0001*mc_deg2radd
    call arra%dealloc()
!
    call readArrayAttributeHDFWrapper(fid,"boxElementDimensions",aria,errmsg)
    if (.level.errmsg == 2) return
    id => aria%get1d()
    this%dimx = id(1); this%dimy = id(2); this%dimz = id(3)
    call aria%dealloc()
    !
    call readArrayHDFWrapper(fid,"nbpvert_per_proc",aria,errmsg)
    if (.level.errmsg == 2) return
    call aria%getDims(dims)
    allocate(this%nbpvert(0:dims(1)-1))
    this%nbpvert = aria%get1d(); call aria%dealloc(); deallocate(dims)
!
    call readArrayHDFWrapper(fid,"nbpbot_per_proc",aria,errmsg)
    if (.level.errmsg == 2) return
    call aria%getDims(dims)
    allocate(this%nbpbot(0:dims(1)-1))
    this%nbpbot = aria%get1d(); call aria%dealloc(); deallocate(dims)
    this%nbptot = sum(this%nbpvert)+sum(this%nbpbot)
!
    call readArrayHDFWrapper(fid,"gllUniqueRadii",arra,errmsg)
    if (.level.errmsg == 2) return
    call arra%getDims(dims)
    allocate(this%runique(dims(1)))                   ! allocate required because it is a pointer
    this%runique = arra%get1d(); call arra%dealloc(); deallocate(dims)
!
    call readArrayHDFWrapper(fid,"gllRadLatLonNormalsSides",arra,errmsg)
    if (.level.errmsg == 2) return
    this%rlatlonnvside = arra%get2d(); call arra%dealloc()
!
    call readArrayHDFWrapper(fid,"gllSortingIndicesSides",aria,errmsg)
    if (.level.errmsg == 2) return
    this%sortidx = aria%get1d(); call aria%dealloc()
!
    call readArrayHDFWrapper(fid,"gllRadLatLonNormalsBottom",arra,errmsg)
    if (.level.errmsg == 2) return
    this%rlatlonnvbot = arra%get2d(); call arra%dealloc()
!
    call h5fclose_f(fid,ierr)
    end subroutine readBoundaryPoints
!--------------------------------------------------------------------------------
!  Write boundary points to HDF    
!
    subroutine writeBoundaryPoints(this,hdffile,errmsg)
    type (boundary_points) :: this
    character (len=*) :: hdffile
    type (error_message) :: errmsg
    integer(hid_t) :: fid
    integer :: ierr
    type (any_rank_real_array) :: arra
    type (any_rank_integer_array) :: aria
!
    call createFileHDFWrapper(trim(hdffile),fid,errmsg)
    if (.level.errmsg == 2) return
!
!  write geographical coordinates of box center
!
    call arra%assoc1d((/ this%latcenter,this%loncenter /))
    call writeArrayAttributeHDFWrapper(fid,"centerLatLon",arra,errmsg)
    if (.level.errmsg == 2) return
!
!  write radial and epicentral distance tolerances
!
    call arra%assoc1d((/ real(this%rtol), real(this%dtol) /))
    call writeArrayAttributeHDFWrapper(fid,"tolerances",arra,errmsg)
    if (.level.errmsg == 2) return
!
!  write number of elements per direction of computational box
!
    call aria%assoc1d((/ this%dimx,this%dimy,this%dimz /))
    call writeArrayAttributeHDFWrapper(fid,"boxElementDimensions",aria,errmsg)
    if (.level.errmsg == 2) return
!
!  write number of points on vertical side boundaries per Specfem process
!
    call aria%assoc1d(this%nbpvert)
    call writeArrayHDFWrapper(fid,"nbpvert_per_proc",aria,errmsg)
    if (.level.errmsg == 2) return
!
!  write number of points on bottom side boundaries per Specfem process
!
    call aria%assoc1d(this%nbpbot)
    call writeArrayHDFWrapper(fid,"nbpbot_per_proc",aria,errmsg)
    if (.level.errmsg == 2) return
!
!  write unique radii to HDF    
!
    call arra%assoc1d(this%runique)
    call writeArrayHDFWrapper(fid,"gllUniqueRadii",arra,errmsg)
    if (.level.errmsg == 2) return
!
!  write set of points and normals on vertical boundaries to HDF    
!
    call arra%assoc2d(this%rlatlonnvside)
    call writeArrayHDFWrapper(fid,"gllRadLatLonNormalsSides",arra,errmsg)
    if (.level.errmsg == 2) return
!
!  write sorting indices telling original positions of
!  points on vertical boundaries in concatenated SPECFEM
!  files to HDF    
!
    call aria%assoc1d(this%sortidx)
    call writeArrayHDFWrapper(fid,"gllSortingIndicesSides",aria,errmsg)
    if (.level.errmsg == 2) return
!
!  write set of bottom points and normals to HDF    
!
    call arra%assoc2d(this%rlatlonnvbot)
    call writeArrayHDFWrapper(fid,"gllRadLatLonNormalsBottom",arra,errmsg)
    if (.level.errmsg == 2) return
!
    call aria%deassoc()
    call arra%deassoc()
!
    call h5fclose_f(fid,ierr)
    end subroutine writeBoundaryPoints
!---------------------------------------------------------------------------------
!  Compute boundary point data from a given list containing cartesian coordinates
!  of GLL-points relative to center of box and outward normal vectors
!  Uses output files procXXXXXX_absorb_dsm and procXXXXXX_normal.txt in DATABASES_MPI    
!
    subroutine createFromFileBoundaryPoints(this,database,nproc,dimx,dimy,dimz,latc,lonc,rtol,dtol,errmsg)
    type (boundary_points) :: this
    character (len=*) :: database
    integer :: nproc
    integer :: dimx,dimy,dimz
    double precision :: latc,lonc
    double precision :: rtol,dtol
    type (error_message) :: errmsg
    character (len=6) :: cip
    character (len = max_length_string) :: filename
    character (len=37) :: myname = "createFromFileCartesianBoundaryPoints"
    real, dimension(:,:), pointer :: p
    double precision, dimension(:,:), allocatable :: pcar
    double precision :: thetac,phic,xg,yg,zg,xr,yr,zr,r,theta,phi
    integer :: ip,nbptop,off,nbp
!
    call addTrace(errmsg,myname)
    this%dimx = dimx; this%dimy = dimy; this%dimz = dimz
    this%latcenter = latc; this%loncenter = lonc
    this%dtol = dtol*mc_deg2radd    !  /  Es ist wahrscheinlich sinnvoll die konvertierung von deg nach rad schon im vorbereitungsschritt zu machen
                                    !dann muss diese codezeile entsprechend geändert werden und die konvertierung muss dann aus der auslesefunktion herausgenommen werden.
    this%rtol = rtol
!
!  get number of boundary points to read
!  remember number of points per process
!
    allocate(this%nbpvert(0:nproc-1))
    allocate(this%nbpbot(0:nproc-1))
    this%nbptot = 0
    do ip = 0,nproc-1
       write(cip,'(i6.6)') ip
       open(1,file = database+'/proc'+cip+'_absorbing_boundaries_for_gemini.txt')
       read(1,*) this%nbpvert(ip)
       read(1,*) this%nbpbot(ip)
       read(1,*) nbptop
       close(1)
       this%nbpvert(ip) = this%nbpvert(ip)*NGLL*NGLL
       this%nbpbot(ip) = this%nbpbot(ip)*NGLL*NGLL
       write(6,'(a,3i10)') 'process: ',ip,this%nbpvert(ip),this%nbpbot(ip)
    enddo
    this%nbptot = sum(this%nbpvert)+sum(this%nbpbot)
   !
   !  read cartesian coordinates of GLL-points (x,y,z) and normals for vertical side boundaries
   !
    nbp = sum(this%nbpvert)
    allocate(pcar(3,nbp))
    allocate(this%rlatlonnvside(6,nbp))         ! radius, lat, lon and normal in one array
    allocate(this%sortidx(nbp))
    off = 0
    do ip = 0,nproc-1
       if (this%nbpvert(ip) /= 0) then
          write(cip,'(i6.6)') ip
          filename = database+'/proc'+cip+'_absorbing_boundaries_for_gemini.txt'
          call readRealMatrixAsciiDataIO(trim(filename),1,p,errmsg,ncol = 6,nrow = this%nbpvert(ip),nskip = 3)
          if (.level.errmsg == 2) then
             return
          endif
          pcar(:,off+1:off+this%nbpvert(ip)) = transpose(dble(p(:,1:3)))
          pcar(3,off+1:off+this%nbpvert(ip)) = pcar(3,off+1:off+this%nbpvert(ip))+REARTH*1000.
           ! z = 0 at r = 0
          this%rlatlonnvside(4:6,off+1:off+this%nbpvert(ip)) = transpose(dble(p(:,4:6)))    ! normals only
          off = off+this%nbpvert(ip)
          deallocate(p)
          write(6,'(a,i6)') "Read vertical boundary data for process: ",ip
       else
          write(6,'(a,i6)') "No vertical boundary data for process: ",ip
       endif
    enddo
!
!  compute spherical coordinates radlatlon in km,degrees
!  first rotate box cartesian coordinates by -90 degrees around z to get x-> south, y->east
!  then convert to global cartesian with x through Greenwich meridian and z through north pole
!  then convert to geographical coordinates r, lat, lon in degrees
!
    thetac = (90.-latc)*mc_deg2radd
    phic = lonc*mc_deg2radd
    do ip = 1,sum(this%nbpvert)
       call coordinatesRCfromLCAxesRotation(-0.5*mc_pid,pcar(1,ip),pcar(2,ip),pcar(3,ip),xr,yr,zr)
       call coordinatesGCfromLCAxesRotation(thetac,phic,xr,yr,zr,xg,yg,zg)
       call coordinatesLSfromLCAxesRotation(xg,yg,zg,r,theta,phi)
       this%rlatlonnvside(1,ip) = r
       this%rlatlonnvside(2,ip) = (0.5*mc_pid-theta)/mc_deg2radd
       this%rlatlonnvside(3,ip) = phi/mc_deg2radd
    enddo
    deallocate(pcar)
    write(6,'(a)') "Minimal and maximal (r,lat,lon) values for vertical boundaries"
    write(6,'(3f12.3)') minval(this%rlatlonnvside(1:3,:),2)
    write(6,'(3f12.3)') maxval(this%rlatlonnvside(1:3,:),2)
!
!  sort points according to values of radius, overwrites rlatlonnvside
!
    call heapSort2D(this%rlatlonnvside,this%sortidx,dimsort = 2,is = 1)
!
!  find unique radii of GLL-points, for which we compute gfk-spectra
!
    call findUniqueRadiiBoundaryPoints(this)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                                  BOTTOM BOUNDARY POINTS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  read cartesian coordinates of GLL-points (x,y,z) and normals for bottom boundary
!
    allocate(pcar(3,sum(this%nbpbot)))
    allocate(this%rlatlonnvbot(6,sum(this%nbpbot)))   ! r,lat,lon and normal in one array
    off = 0
    do ip = 0,nproc-1
       if (this%nbpbot(ip) /= 0) then
          write(cip,'(i6.6)') ip
          filename = database+'/proc'+cip+'_absorbing_boundaries_for_gemini.txt'
          call readRealMatrixAsciiDataIO(trim(filename),1,p,errmsg,ncol = 6,nrow = this%nbpbot(ip),nskip = 3+this%nbpvert(ip))
          if (.level.errmsg == 2) then
             return
          endif
          pcar(:,off+1:off+this%nbpbot(ip)) = transpose(dble(p(:,1:3)))
          pcar(3,off+1:off+this%nbpbot(ip)) = pcar(3,off+1:off+this%nbpbot(ip))+REARTH*1000.
          this%rlatlonnvbot(4:6,off+1:off+this%nbpbot(ip)) = transpose(dble(p(:,4:6)))       ! normals only
          off = off+this%nbpbot(ip)
          deallocate(p)
          write(6,'(a,i6)') "Read bottom boundary data for process: ",ip 
       else
          write(6,'(a,i6)') "No bottom boundary data for process: ",ip
       endif
    enddo
!
!  compute spherical coordinates lat,lon in degrees
!  first rotate box cartesian coordinates by -90 degrees around z to get x-> south, y->east
!  then convert to global cartesian with x through Greenwich meridian and z through north pole
!  then convert to geographical coordinates r, lat, lon in degrees
!
    do ip = 1,sum(this%nbpbot)
       call coordinatesRCfromLCAxesRotation(-0.5*mc_pid,pcar(1,ip),pcar(2,ip),pcar(3,ip),xr,yr,zr)
       call coordinatesGCfromLCAxesRotation(thetac,phic,xr,yr,zr,xg,yg,zg)
       call coordinatesLSfromLCAxesRotation(xg,yg,zg,r,theta,phi)
       this%rlatlonnvbot(1,ip) = r
       this%rlatlonnvbot(2,ip) = (0.5*mc_pid-theta)/mc_deg2radd
       this%rlatlonnvbot(3,ip) = phi/mc_deg2radd
    enddo
    write(6,'(a)') "Minimal and maximal (r,lat,lon)values for bottom boundary"
    write(6,'(3f12.3)') minval(this%rlatlonnvbot(1:3,:),2)
    write(6,'(3f12.3)') maxval(this%rlatlonnvbot(1:3,:),2)
    deallocate(pcar)
!    
    end subroutine createFromFileBoundaryPoints
!------------------------------------------------------------------------
!  Find unique radii of boundary points
!
    subroutine findUniqueRadiiBoundaryPoints(this)
    type (boundary_points) :: this
    real :: r
    integer :: i,count,nr
!
    nr = this%dimz*NGLL
    allocate(this%runique(nr))
    r = this%rlatlonnvside(1,1)
    this%runique(1) = r
    count = 1
    do i = 2,sum(this%nbpvert)
       if (abs(r-this%rlatlonnvside(1,i)) > this%rtol) then
          count = count+1
          if (count > size(this%runique)) then
             this%runique => reallocate(this%runique,size(this%runique)+nr)
          endif
          r = this%rlatlonnvside(1,i)
          this%runique(count) = r
       endif
    enddo
    this%runique => reallocate(this%runique,count)
    end subroutine findUniqueRadiiBoundaryPoints
!-----------------------------------------------------------------------------------------
!  Find set of unique epicentral distances for all points on vertical or bottom boundaries
!  for a given source location
!  Returns an array of unique distances and
!  an array of indices mapping the pont index to the correspoding distance index
!  so we know which distance belongs to som point on the boundary.
!  bkey:         single character specifying either vertical or bottom boundaries (v,b)
!  thetas:       co-latitude of source in rad
!  phis:         longitude of source in rad
!  dtol:         distance tolerance in rad
!  uniquedis:    1D-array of unique epicentral distances
!  xi:           receiver azimuth at source counted from south over east (all points)
!  point2dis:    1D-array of indices: distance_index = point2dis(point_index) 
!  delmin:       minimal epicentral distance of boundary points (optional)
!  delmax:       maximal epicentral distance of boundary points (optional)
!
    subroutine findUniqueDistancesBoundaryPoints(this,bkey,thetas,phis,myrank,uniquedis,xi,point2dis,delmin,delmax)
    type (boundary_points) :: this
    character (len=1) :: bkey
    double precision :: thetas,phis
    double precision, dimension(:), pointer :: uniquedis
    integer, dimension(:), allocatable :: point2dis
    double precision, dimension(:,:), allocatable :: delta
    double precision, dimension(:), allocatable :: xi
    integer, dimension(:), allocatable :: idx
    double precision :: xg,yg,zg,rr,del,x,y,z,r,thetab,phib
    double precision, optional :: delmin,delmax
    integer :: nbp,nd,ip,count,myrank
!
!  calculate epicentral distance in rad
!
    if (bkey == 'v') then
       nbp = sum(this%nbpvert)
    else
       nbp = sum(this%nbpbot)
    endif
    allocate(delta(1,nbp),xi(nbp))
    allocate(idx(nbp))
    allocate(point2dis(nbp))
    do ip = 1,nbp
       if (bkey == 'v') then
          r = this%rlatlonnvside(1,ip)
          thetab = 0.5d0*mc_pid-this%rlatlonnvside(2,ip)*mc_deg2radd
          phib = this%rlatlonnvside(3,ip)*mc_deg2radd
       else
          r = this%rlatlonnvbot(1,ip)
          thetab = 0.5d0*mc_pid-this%rlatlonnvbot(2,ip)*mc_deg2radd
          phib = this%rlatlonnvbot(3,ip)*mc_deg2radd
       endif
       call coordinatesLCfromLSAxesRotation(r,thetab,phib,xg,yg,zg)         ! global cartesian coordinates of BP
       call coordinatesLCfromGCAxesRotation(thetas,phis,xg,yg,zg,x,y,z)     ! cartesian coordinates of BP with origin at source
       call coordinatesLSfromLCAxesRotation(x,y,z,rr,delta(1,ip),xi(ip))    ! spherical coordinates of BP with origin at source
    enddo
!
!   minimum and maximum epicentral distance
!
    if (present(delmin)) delmin = minval(delta(1,:))
    if (present(delmax)) delmax = maxval(delta(1,:))
!
!  sort distances
!
    call heapSort2D(delta,idx,dimsort = 2,is = 1)
!
!  find unique distances
!
    nd = this%dimz*NGLL
    allocate(uniquedis(nd))
    del = delta(1,1)
    uniquedis(1) = del
    count = 1
    do ip = 1,nbp
       if (abs(delta(1,ip)-del) > this%dtol) then
          count = count+1
          if (count > size(uniquedis)) then
             uniquedis => reallocate(uniquedis,size(uniquedis)+nd)
          endif
          del = delta(1,ip)
          uniquedis(count) = del
       endif
       point2dis(idx(ip)) = count
    enddo
    uniquedis => reallocate(uniquedis,count)
    deallocate(idx,delta)
    end subroutine findUniqueDistancesBoundaryPoints
!---------------------------------------------------------------------------
!  Iterator over GLL-points on vertical boundaries sorted according to radius
!
    function iterateVerticalBoundaryPoints(this,cbp,offst,stp) result(next)
    type (boundary_points) :: this
    type (single_boundary_point) :: cbp
    integer, optional :: stp,offst
    logical :: next
    integer, save :: ip = 1, count = 0
    integer :: step,offset
!
    if (present(stp)) then; step = stp; else; step = 1; endif
    if (present(offst)) then; offset = offst; else; offset = 0; endif
!
!  initialy add offset, later add step
!
    if (count == 0) then
       ip = ip+offset
    else
       ip = ip+step
    endif
    if (ip > sum(this%nbpvert)) then
       next = .false.
       ip = 1
       count = 0
       return
    endif
    next = .true.        
    count = count+1
    cbp%current = ip
    cbp%r = this%rlatlonnvside(1,ip)
    cbp%lat = this%rlatlonnvside(2,ip)
    cbp%lon = this%rlatlonnvside(3,ip)
    cbp%nx = this%rlatlonnvside(4,ip)
    cbp%ny = this%rlatlonnvside(5,ip)
    cbp%nz = this%rlatlonnvside(6,ip)
    end function iterateVerticalBoundaryPoints
!-----------------------------------------------------------------------
!  Iterator over GLL-points on bottom boundary ordered according to
!  elements
!
    function iterateBottomBoundaryPoints(this,cbp,offst,stp) result(next)
    type (boundary_points) :: this
    type (single_boundary_point) :: cbp
    integer, optional :: stp,offst
    logical :: next
    integer, save :: ibp = 1,count = 0
    integer :: step,offset
!
    if (present(stp)) then; step = stp; else; step = 1; endif
    if (present(offst)) then; offset = offst; else; offset = 0; endif
!
!  initialy add offset, later add step
!
    if (count == 0) then
       ibp = ibp+offset
    else
       ibp = ibp+step
    endif
!
    if (ibp > sum(this%nbpbot)) then
       next = .false.
       ibp = 1
       count = 0
       return
    endif
    next = .true.
!
!  calculate global element count    
!
    count = count+1
    cbp%current = ibp
    cbp%r = this%rlatlonnvbot(1,ibp)
    cbp%lat = this%rlatlonnvbot(2,ibp)
    cbp%lon = this%rlatlonnvbot(3,ibp)
    cbp%nx = this%rlatlonnvbot(4,ibp)
    cbp%ny = this%rlatlonnvbot(5,ibp)
    cbp%nz = this%rlatlonnvbot(6,ibp)
    end function iterateBottomBoundaryPoints
!----------------------------------------------------------------
!  Get total number of boundary points
!
    function getNumTotalBoundaryPoints(this) result(res)
    type (boundary_points) :: this
    integer :: res
    res = this%nbptot
    end function getNumTotalBoundaryPoints
!----------------------------------------------------------------
!  Get total number of points on vertical boundaries
!
    function getSidesNumTotalBoundaryPoints(this) result(res)
    type (boundary_points) :: this
    integer :: res
    res = sum(this%nbpvert)
    end function getSidesNumTotalBoundaryPoints
!----------------------------------------------------------------
!  Get number of bottom grid points on bottom boundary
!
    function getBottomNumTotalBoundaryPoints(this) result(nlat)
    type (boundary_points) :: this
    integer :: nlat
    nlat = sum(this%nbpbot)
    end function getBottomNumTotalBoundaryPoints
!----------------------------------------------------------------
!  Get number of total boundary points per process
!
    function getNumPerProcBoundaryPoints(this,ip) result(res)
    type (boundary_points) :: this
    integer :: res,ip
    res = this%nbpvert(ip)+this%nbpbot(ip)
    end function getNumPerProcBoundaryPoints
!----------------------------------------------------------------
!  Get number of sides boundary points per process
!
    function getSidesNumPerProcBoundaryPoints(this,ip) result(res)
    type (boundary_points) :: this
    integer :: res,ip
    res = this%nbpvert(ip)
    end function getSidesNumPerProcBoundaryPoints
!----------------------------------------------------------------
!  Get number of sides boundary points per process
!
    function getBottomNumPerProcBoundaryPoints(this,ip) result(res)
    type (boundary_points) :: this
    integer :: res,ip
    res = this%nbpbot(ip)
    end function getBottomNumPerProcBoundaryPoints
!----------------------------------------------------------------
!  Get offset in proc specific HDF for given bp-index
!  for vertical sides boundary points    
!  ip:     index of radially sorted bp-point 
!  offset: true offset in original SPECFEM process specific bp-file
!  iproc:  process rank owning bp-point
    
    subroutine getProcOffsetSidesBoundaryPoints(this,ip,offset,iproc)
    type (boundary_points) :: this
    integer :: ip,iproc,offset,iloc
!
!  use unsorted index (original one in SPECFEM file)
!
    iloc = this%sortidx(ip)
    iproc = 0
    do while(iloc > this%nbpvert(iproc))
       iloc = iloc-this%nbpvert(iproc)
       iproc = iproc+1
    enddo
    offset = iloc-1
    end subroutine getProcOffsetSidesBoundaryPoints
!----------------------------------------------------------------
!  Get offset in proc specific HDF for given bp-index
!  for bottom boundary points    
!  ip:     index of bottom bp-point (starts at 1) 
!  offset: true offset in original SPECFEM process specific bp-file
!  iproc:  process rank owning bp-point in SPECFEM
    
    subroutine getProcOffsetBottomBoundaryPoints(this,ip,offset,iproc)
    type (boundary_points) :: this
    integer :: ip,iproc,offset,iloc
!
    iloc = ip
    iproc = 0
    do while(iloc > this%nbpbot(iproc))
       iloc = iloc-this%nbpbot(iproc)
       iproc = iproc+1
    enddo
    offset = iloc-1+this%nbpvert(iproc)  ! append to vertical side BPs
    end subroutine getProcOffsetBottomBoundaryPoints
!----------------------------------------------------------------
!  Get center of box in radians
!
    subroutine getBoxCenterBoundaryPoints(this,thetac,phic)
    type (boundary_points) :: this
    double precision :: thetac,phic 
    thetac = (90.d0-this%latcenter)*mc_deg2radd
    phic = this%loncenter*mc_deg2radd
    end subroutine getBoxCenterBoundaryPoints
!----------------------------------------------------------------
!  Get number of processes 
!
    function getNumProcBoundaryPoints(this) result(res)
    type (boundary_points) :: this
    integer :: res
    res = size(this%nbpvert)
    end function getNumProcBoundaryPoints
!----------------------------------------------------------------
!  Get current boundary point counter
!
    function getCounterSingleBoundaryPoint(this) result(res)
    type (single_boundary_point) :: this
    integer :: res
    res = this%current
    end function getCounterSingleBoundaryPoint
!---------------------------------------------------------------------------
!  get epicentral coordinates of single boundary point in rad
!  thetas:       co-latitude of source in rad
!  phis:         longitude of source in rad
!
    subroutine getEpicentralCoordinatesSingleBoundaryPoint(this,thetas,phis,delta,xi)
    type (single_boundary_point) :: this
    double precision :: delta,xi,thetas,phis
    double precision :: r,thetab,phib,xg,yg,zg,x,y,z,rr
    r = this%r
    thetab = 0.5d0*mc_pid-this%lat*mc_deg2radd
    phib = this%lon*mc_deg2radd
    call coordinatesLCfromLSAxesRotation(r,thetab,phib,xg,yg,zg)           ! global cartesian coordinates of BP
    call coordinatesLCfromGCAxesRotation(thetas,phis,xg,yg,zg,x,y,z)       ! cartesian coordinates of BP with origin at source
    call coordinatesLSfromLCAxesRotation(x,y,z,rr,delta,xi)                ! spherical coordinates of BP with origin at source
    end subroutine getEpicentralCoordinatesSingleBoundaryPoint
!---------------------------------------------------------------------------
!  get geographical coordinates of center of box
!  from geographical and local cartesian coordinates
!
    subroutine locationBoxCenter(r,lat,phi,x,y,z,latc,lonc)
    double precision :: r,lat,phi,x,y,z,latc,lonc
    double precision :: theta,xs,ys,zs,delta,xi,thetac,phis
    double precision :: cd,ct,sd,st,cp,cx,sx
!
    theta = (90.d0-lat)*mc_pid/180.d0
    phi = phi*mc_pid/180.d0
    x = x*1.d-3
    y = y*1.d-3
    z = z*1.d-3
!
!  switch to x south and y east and z counted from center of sphere
!
    xs = -y
    ys = x
    zs = z+6371.d0
    call coordinatesLSfromLCAxesRotation(xs,ys,zs,r,delta,xi)
!
!  get longitude difference between C and P
!
    cd = cos(delta); ct = cos(theta); sd = sin(delta); st = sin(theta)
    cx = cos(xi); sx = sin(xi)
    phis = asin(-sd*sx/st)
    cp = cos(phis)
!
!  get theta of origin of local box
!
    thetac = asin((cd**2-ct**2)/(st*cd*cp+sd*ct*cx))
    latc = 0.5d0*mc_pid-thetac
    lonc = phi+phis
    end subroutine locationBoxCenter
!
end module boundaryPoints
