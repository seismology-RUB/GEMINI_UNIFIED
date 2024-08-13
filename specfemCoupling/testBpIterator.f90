program testBpIterator
       use hdf5
    use errorMessage
    use hdfWrapper
    use boundaryPoints
    use mathConstants
    implicit none
    integer :: ierr
    type (error_message) :: errmsg
    type (boundary_points) :: bp
    type (single_boundary_point) :: cbp
    double precision :: thetac,phic,theta,phi,r
    double precision :: xg,yg,zg,xr,yr,zr,x,y,z
    integer :: ibp,iproc,off
!
!  open HDF-file for output
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) then; call print(errmsg); stop; endif
!
!  read boundary points from HDF
!
    call readBoundaryPoints(bp,"specfem_boundary_data.hdf",errmsg)
    if (.level.errmsg == 2) goto 10
!
!  test iterators
!
    thetac = 0.5d0*mc_pid-bp%latcenter*mc_deg2radd
    phic = bp%loncenter*mc_deg2radd
    do while (iterateVerticalBoundaryPoints(bp,cbp,0,1))
       ibp = cbp%current
       r = cbp%r
       theta = 0.5d0*mc_pid-cbp%lat*mc_deg2radd
       phi = cbp%lon*mc_deg2radd
       call getProcOffsetSidesBoundaryPoints(bp,ibp,off,iproc)
       call coordinatesLCfromLSAxesRotation(r,theta,phi,xg,yg,zg)
       call coordinatesLCfromGCAxesRotation(thetac,phic,xg,yg,zg,xr,yr,zr)
       call coordinatesRCfromLCAxesRotation(+0.5*mc_pid,xr,yr,zr,x,y,z)
       if (iproc == 4) then
          write(6,'(3i10,6f12.3)') ibp,iproc,off+1,r,cbp%lat,cbp%lon,x,y,z-6371.d0
       endif
    enddo
    do while (iterateBottomBoundaryPoints(bp,cbp,0,1))
       ibp = cbp%current
       r = cbp%r
       theta = 0.5d0*mc_pid-cbp%lat*mc_deg2radd
       phi = cbp%lon*mc_deg2radd
       call getProcOffsetBottomBoundaryPoints(bp,ibp,off,iproc)
       call coordinatesLCfromLSAxesRotation(r,theta,phi,xg,yg,zg)
       call coordinatesLCfromGCAxesRotation(thetac,phic,xg,yg,zg,xr,yr,zr)
       call coordinatesRCfromLCAxesRotation(+0.5*mc_pid,xr,yr,zr,x,y,z)
       if (iproc == 4) then
          write(6,'(3i10,6f12.3)') ibp,iproc,off+1,r,cbp%lat,cbp%lon,x,y,z-6371.d0
       endif
    enddo
!        
!  close HDF stuff and dealloc
!
    call h5close_f(ierr)
    call dealloc(bp)
!
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif    
!
end program testBpIterator
