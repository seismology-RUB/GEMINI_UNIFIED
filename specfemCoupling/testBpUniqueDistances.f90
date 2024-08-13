program testBpUniqueDistances
       use hdf5
    use argumentParser
    use inputParameter
    use errorMessage
    use hdfWrapper
    use boundaryPoints
    implicit none
    type (argument_parser) :: ap
    type (input_parameter) :: inpar
    type (error_message) :: errmsg
    type (boundary_points) :: bp
    double precision, dimension(:), pointer :: uniquedis
    double precision, dimension(:), allocatable :: phi
    integer, dimension(:), allocatable :: point2dis
    integer :: ierr,j
    double precision :: dtol
    character (len=max_length_string) :: database,bpfile,parfile
    character (len=80), dimension(3) :: par_keys
    character (len=21) :: myname = "testBpUniqueDistances"
!
!  keywords for input parameters
!
    data par_keys/'SPECFEM_DATABASES_MPI','BOUNDARY_DATA_OUTPUT','DISTANCE_TOLERANCE'/
!-----------------------------------------------------------------------------------------------------
    call init(ap,myname,'test unique distances calculation')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call parse(ap)
    parfile = ap.sval.'parfile'
    if (.level.(.errmsg.ap) == 2) then
       call usage(ap); call print(.errmsg.ap); stop;
    endif
    call document(ap); call dealloc(ap)
!
    call new(errmsg,myname)
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 10
    database = inpar.sval.'SPECFEM_DATABASES_MPI'
    bpfile = inpar.sval.'BOUNDARY_DATA_OUTPUT'
    dtol = (inpar.rval.'DISTANCE_TOLERANCE')*mc_deg2radd
!
!  open HDF-file for output
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
!
!  read boundary points from HDF
!
    call readBoundaryPoints(bp,bpfile,errmsg)
    if (.level.errmsg == 2) goto 10
!
    call findUniqueDistancesBoundaryPoints(bp,'v',0.5d0*mc_pid,-0.5d0*mc_pid,dtol,uniquedis,phi,point2dis)
    do j = 1,size(uniquedis)
       write(6,'(i6,f12.3)') j,uniquedis(j)/mc_deg2radd
    enddo
    deallocate(uniquedis,phi,point2dis)
!
    call findUniqueDistancesBoundaryPoints(bp,'b',0.5d0*mc_pid,-0.5d0*mc_pid,dtol,uniquedis,phi,point2dis)
    do j = 1,size(uniquedis)
       write(6,'(i6,f12.3)') j,uniquedis(j)/mc_deg2radd
    enddo
    deallocate(uniquedis,phi,point2dis)
!        
!  close HDF stuff and dealloc
!
    call h5close_f(ierr)
    call dealloc(bp)
!
!  treat errors
!
10  if (.level.errmsg == 2) then
       call print(errmsg)
    endif
end program testBpUniqueDistances
