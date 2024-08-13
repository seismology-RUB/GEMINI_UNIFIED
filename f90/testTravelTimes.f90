program testTravelTimes
   use hdfWrapper
   use argumentParser
   use travelTimes
   use mathConstants
   use errorMessage
   implicit none
   type (argument_parser) :: ap
   type (raytable_travel_times) :: raytable
   type (error_message) :: errmsg
   type (any_rank_integer_array) :: aria
   type (any_rank_real_array) :: arra
   real, dimension(:), allocatable :: d
   integer, dimension(:), allocatable :: id
   integer :: ndel,nr,j,k
   integer (kind=8) :: hdf_xferprp,fid
   real, dimension(:,:), allocatable :: ttg,pg
   double precision :: tt,deltain,dstep,dr,rs,re,delmin,delmax,rmin,rmax,slow
   character(len=max_length_string) :: tablefile,ttgridfile
   character(len=15) :: myname = 'testTravelTimes'
!----------------------------------------------------------------------------
   call init(ap,myname,'Travel times in a vertical plane')
   call addPosarg(ap,'ttable','sval','travel time table file')
   call addOption(ap,'-gridfile',.true.,'Name of HDF tt grid file','sval','ttgrid.hdf')
   call addOption(ap,'-rs',.true.,'Source radius','dval','6371000.0')
   call addOption(ap,'-delmin',.true.,'Min distance (deg)','dval','0.0')
   call addOption(ap,'-delmax',.true.,'Max distance (deg)','dval','120.0')
   call addOption(ap,'-ndel',.true.,'Number of distances','ival','100')
   call addOption(ap,'-rmin',.true.,'Min distance (deg)','dval','5800000.0')
   call addOption(ap,'-rmax',.true.,'Max distance (deg)','dval','6371000.0')
   call addOption(ap,'-nr',.true.,'Number of radii','ival','100')
   call parse(ap)
!
   tablefile = ap.sval.'ttable'
   ttgridfile = ap.sval.'-gridfile'
   rs = ap.dval.'-rs'
   delmin = mc_deg2radd*(ap.dval.'-delmin')
   delmax = mc_deg2radd*(ap.dval.'-delmax')
   ndel = ap.ival.'-ndel'
   rmin = ap.dval.'-rmin'
   rmax = ap.dval.'-rmax'
   nr = ap.ival.'-nr'
   if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); stop; endif
   call document(ap)
   call dealloc(ap)
!-----------------------------------------------------------------------------
   call new(errmsg,myname)
!
!  open HDF environment
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 10
   call setXferprpIndependentHDFWrapper(hdf_xferprp,errmsg)
   if (.level.errmsg == 2) goto 10
!-----------------------------------------------------------------------------
   call readTableTravelTimes(raytable,tablefile,errmsg)
   if (.level.errmsg == 2) goto 10
!
   allocate(ttg(ndel,nr),pg(ndel,nr))
   dstep = (delmax-delmin)/(ndel-1)
   dr = (rmax-rmin)/(nr-1)
   call setSourceInfoTravelTimes(raytable,rs,errmsg)
   if (.level.errmsg == 2) goto 10
!
!  walk through grid
!
   do k = 1,nr
      re = rmin+(k-1)*dr
      do j = 1,ndel
         deltain = delmin+(j-1)*dstep
         call getReceiverTravelTimes(raytable,re,deltain,tt,errmsg,slow)
         ttg(j,k) = tt
         pg(j,k) = slow
         if (.level.errmsg == 2) goto 10
      end do
   enddo
!
!  write regular grid to HDF file
!
   call createFileHDFWrapper(trim(ttgridfile),fid,errmsg)
   if (.level.errmsg == 2) return
!
   id = [ndel,nr]
   call aria%assoc1d(id)
   call writeArrayAttributeHDFWrapper(fid,"dimensions",aria,errmsg)
   call aria%deassoc(); deallocate(id)
   if (.level.errmsg == 2) goto 10
!
   d = [real(dstep),real(dr)]
   call arra%assoc1d(d)
   call writeArrayAttributeHDFWrapper(fid,"grid_spacings_rad_m",arra,errmsg)
   call arra%deassoc(); deallocate(d)
   if (.level.errmsg == 2) goto 10
!
   d = [real(delmin),real(rmin)]
   call arra%assoc1d(d)
   call writeArrayAttributeHDFWrapper(fid,"minimum_values_rad_m",arra,errmsg)
   call arra%deassoc(); deallocate(d)
   if (.level.errmsg == 2) goto 10
!
   call arra%assoc2d(ttg)
   call writeArrayHDFWrapper(fid,'ttime',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 10
!
   call arra%assoc2d(pg)
   call writeArrayHDFWrapper(fid,'slowness',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 10
!
   call closeFileHDFWrapper(fid,errmsg)
   if (.level.errmsg == 2) return
!
   call closeEnvironmentHDFWrapper(errmsg)
   deallocate(ttg,pg)
   call dealloc(errmsg)
!
!  error exit
!
10 if (.level.errmsg == 2) then
      call print(errmsg)
   endif
end program testTravelTimes
