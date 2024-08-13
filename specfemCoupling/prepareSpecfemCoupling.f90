program prepareSpecfemCoupling
       use hdf5
    use asciiDataIO
    use argumentParser
    use inputParameter
    use errorMessage
    use hdfWrapper
    use boundaryPoints
    use string
    implicit none
    integer, dimension(3) :: dims
    double precision, dimension(2) :: boxc
    integer :: ierr,nproc
    double precision :: rtol,dtol
    type (argument_parser) :: ap
    type (error_message) :: errmsg
    type (boundary_points) :: bp
    type (input_parameter) :: inpar
    character (len=max_length_string) :: database,outfile,parfile
    character (len=80), dimension(7) :: par_keys
    character (len=22) :: myname = "prepareSpecfemCoupling"
!
!  keywords for input parameters
!
    data par_keys/'SPECFEM_DATABASES_MPI','BOUNDARY_DATA','BOX_DIMENSIONS','BOX_CENTER',&
         'NPROC','RADIAL_TOLERANCE','DISTANCE_TOLERANCE'/
!-----------------------------------------------------------------------------------------------------
    call init(ap,myname,'Prepare boundary points for specfem coupling')
    call addPosarg(ap,'parfile','sval',' Parameter file')
    call parse(ap)
    parfile = ap.sval.'parfile'
    if (.level.(.errmsg.ap) == 2) then
       call usage(ap); call print(.errmsg.ap); stop
    endif
    call document(ap); call dealloc(ap)
!
    call new(errmsg,myname)
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,1,parfile,errmsg)
    if (.level.errmsg == 2) goto 10
    database = inpar.sval.'SPECFEM_DATABASES_MPI'
    outfile = inpar.sval.'BOUNDARY_DATA'
    dims = ivec(inpar,'BOX_DIMENSIONS',3,ierr)
    boxc = dvec(inpar,'BOX_CENTER',2,ierr)
    nproc = inpar.ival.'NPROC'
    rtol = inpar.dval.'RADIAL_TOLERANCE'
    dtol = inpar.dval.'DISTANCE_TOLERANCE'
!
    call createFromFileBoundaryPoints(bp,database,nproc,dims(1),dims(2),dims(3),boxc(1),boxc(2),rtol,dtol,errmsg)
    if (.level.errmsg == 2) goto 10
!
!  open HDF-file for output
!
    call openEnvironmentHDFWrapper(errmsg)
    if (.level.errmsg == 2) goto 10
!
!  write boundary points to HDF
!
    call writeBoundaryPoints(bp,outfile,errmsg)
    if (.level.errmsg == 2) goto 10
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
!
end program prepareSpecfemCoupling
