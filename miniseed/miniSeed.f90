!----------------------------------------------------------------------------
!> \brief  Module to extract time window from mseed file
!----------------------------------------------------------------------------
 module miniSeed
 use iso_c_binding
	use dateTime
	use errorMessage
	implicit none
	interface
		function unpackTraceMSeed(mst,n,dt,sdate,tfs,tns,network,station,channel,location,y,mstfree) &
		                               &  bind(c,name="unpackTraceMSeed")
		use, intrinsic :: iso_c_binding
		integer (kind = c_int) :: unpackTraceMSeed
		type (c_ptr), value :: mst                   ! value attribute is important !!!
		character(kind = c_char), dimension(*) :: network,station,channel,location,sdate
		integer(kind = c_int ) :: n,tfs,tns
		real(kind = c_double) :: dt
		real(kind = c_float), dimension(*) :: y
		integer (kind = c_int) :: mstfree
		end function unpackTraceMSeed
	end interface
	interface
		function readTraceMSeed(filename,sdt,edt,nsamp) bind(c,name="readTraceMSeed")
		use, intrinsic :: iso_c_binding
		type (c_ptr) :: readTraceMSeed
		character (kind = c_char), dimension(*) :: filename,sdt,edt
		integer (kind = c_int) :: nsamp
		end function readTraceMSeed
	end interface
	interface
		function readFullTraceMSeed(filename,nsamp) bind(c,name="readFullTraceMSeed")
		use, intrinsic :: iso_c_binding
		type (c_ptr) :: readFullTraceMSeed
		character (kind = c_char), dimension(*) :: filename
		integer (kind = c_int) :: nsamp
		end function readFullTraceMSeed
	end interface
	interface
		function writeTraceMSeed(filename,idata,nsamp,samprate,timestr, &
		                        & network,station,channel,location) bind(c,name="writeTraceMSeed")
		use, intrinsic :: iso_c_binding
		integer (kind = c_int) :: writeTraceMSeed
		character (kind = c_char), dimension(*) :: filename,network,station,channel,location,timestr
		integer (kind = c_int) :: nsamp
		integer (kind = c_int), dimension(*) :: idata
		real (kind = c_double) :: samprate
		end function writeTraceMSeed
	end interface
	interface
		function readTraceGroupMSeed(filename, numtraces) bind(c,name="readTraceGroupMSeed")
		use, intrinsic :: iso_c_binding
		type (c_ptr) :: readTraceGroupMSeed
		character (kind = c_char), dimension(*) :: filename
		integer (kind = c_int) :: numtraces 
		end function readTraceGroupMSeed
	end interface
	interface
		function getFirstTraceMSeed(mstg,nsamp) bind(c,name="getFirstTraceMSeed")
		use, intrinsic :: iso_c_binding
		type (c_ptr) :: getFirstTraceMSeed
		type (c_ptr), value :: mstg                   ! value attribute is important !!!
		integer (kind = c_int) :: nsamp
		end function getFirstTraceMSeed
	end interface
	interface
		function getNextTraceMSeed(mstcur,nsamp) bind(c,name="getNextTraceMSeed") 
		use, intrinsic :: iso_c_binding
		type (c_ptr) :: getNextTraceMSeed
		type (c_ptr), value :: mstcur                  ! value attribute is important !!!
		integer (kind = c_int) :: nsamp
		end function getNextTraceMSeed
	end interface
	interface operator (.tanf.); module procedure tanfMiniSeed; end interface
	interface operator (.tanfdp.); module procedure tanfDPMiniSeed; end interface
	interface operator (.nsamp.); module procedure nsampMiniSeed; end interface
	interface operator (.tfs.); module procedure tfsMiniSeed; end interface
	interface operator (.tns.); module procedure tnsMiniSeed; end interface
	interface operator (.dt.); module procedure dtMiniSeed; end interface
	interface operator (.network.); module procedure networkMiniSeed; end interface
	interface operator (.station.); module procedure stationMiniSeed; end interface
	interface operator (.location.); module procedure locationMiniSeed; end interface
	interface operator (.sdate.); module procedure sdateMiniSeed; end interface
	interface operator (.channel.); module procedure channelMiniSeed; end interface
	interface operator (.trace.); module procedure traceMiniSeed; end interface
	interface operator (.tstart.); module procedure tstartMiniSeed; end interface
	interface operator (.tend.); module procedure tendMiniSeed; end interface
	interface dealloc; module procedure deallocMiniSeed; end interface
	type mini_seed
		private
		integer :: n
		double precision :: dt
		type (date_time) :: tstart
		character (len=2) :: network
		character (len=5) :: station
		character (len=3) :: channel
		character (len=2) :: location
		real, dimension(:), pointer :: y => null()
		integer, dimension(:), pointer :: iy => null()
		logical :: link
	end type
!
 contains
!-------------------------------------------------------------------------------
!> \brief Create a mini_seed object from data already in memory
!
	subroutine createLinkMiniSeed(this,n,dt,tstart,network,station,channel,location,iy)
	type (mini_seed) :: this
	integer :: n
	double precision :: dt
	type (date_time) :: tstart
	character (len=*) :: network,station,channel,location
	integer, dimension(:), target :: iy
	this%n = n; this%dt = dt; this%tstart = copyDateTime(tstart)
	this%network = network; this%station = station
	this%channel = channel; this%location = location
	this%iy => iy
	this%link = .true.
	end subroutine createLinkMiniSeed
!-------------------------------------------------------------------------------------------
!> \brief Create mini_seed object from unpacked data (internal use)
!
	subroutine createFromUnpackedMiniSeed(this,c_n,c_dt,c_sdate,c_tfs,c_tns,c_network,c_station,c_channel,c_location,c_y)
	type (mini_seed) :: this
	real (kind = c_float), dimension(:) :: c_y
	character (kind = c_char), dimension(:) :: c_network,c_station,c_channel,c_location
	character (kind = c_char), dimension(:) :: c_sdate
	integer (kind = c_int):: c_n,c_tfs,c_tns
	real (kind = c_double) :: c_dt
	integer :: tfs,tns,hh,mm,ss,year,month,day,doy
	character (len=8) :: sdate
!
	allocate(this%y(c_n))
	this%y = c_y
	this%dt = c_dt; this%n = c_n
	tns = c_tns; tfs = c_tfs
	this%link = .false.
	call convertCharStringMiniSeed(c_sdate,sdate,.false.)	
	call convertCharStringMiniSeed(c_network,this%network,.false.)
	call convertCharStringMiniSeed(c_station,this%station,.false.)
	call convertCharStringMiniSeed(c_channel,this%channel,.false.)
	call convertCharStringMiniSeed(c_location,this%location,.false.)
	call hmsFromTsecTimeUtils(tfs,hh,mm,ss)
	read(sdate,'(i4,i2,i2)') year,month,day
	doy = getDayOfYearTimeUtils(year,month,day)
	call new(this%tstart,year,doy,hh,mm,ss,tns)
	end subroutine createFromUnpackedMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Append mini_seed object to existing one
!> \param this Mini seed object to be extended, modified after the call
!> \param ms Mini seed to be appended
!
	subroutine appendMiniSeed(this,ms,errmsg)
	type (mini_seed) :: this,ms
	type (error_message) :: errmsg
	type (date_time) :: tend,msend,tenddt
	integer :: nval,i
	double precision :: err
!
	tend = tendMiniSeed(this)
	tenddt = addDateTime(tend,this%dt)     ! one dt after end of this
	msend = tendMiniSeed(ms)
!
!  condition to be satisfied: msend > tend+dt AND ms%tstart <= tend+dt
!
	if (msend > tenddt .and. (ms%tstart < tenddt .or. ms%tstart == tenddt)) then
		continue
	else
		print *,'appendMiniSeed: time windows do not fit together'
		call printDateTime(ms%tstart); call printDateTime(msend)
		call printDateTime(this%tstart); call printDateTime(tend)
		call add(errmsg,2,'time windows do not fit together','appendMiniSeed')
		return
	endif
!
!  sampling intervals between end of ms and end of this
!
	nval = intervalsBetweenDateTime(msend,tend,this%dt,err)
	if (err > 1.d-6 .and. err <= 1.d-3) then
		print *,'appendMiniSeed: WARNING: Timing inconsistency of ',err,' seconds'
	else if (err > 1.d-3) then
		print *,'appendMiniSeed: ERROR: Timing inconsistency of ',err,' seconds'
		call add(errmsg,2,'Timing inconsistency > 1.e-3 s','appendMiniSeed')
		return
	endif
!
!  add last nval samples of ms to this
!
	this%y => reallocate(this%y,this%n+nval)
	do i = 1,nval
		this%y(this%n+i) = ms%y(ms%n-nval+i)
	enddo
	this%n = this%n+nval
!	print *,'appendMiniSeed: appended ',nval,' samples'
!
	end subroutine appendMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get channel
!
	character (len=3) function channelMiniSeed(this)
	type (mini_seed), intent(in) :: this
	channelMiniSeed = this%channel
	end function channelMiniSeed
!--------------------------------------------------------------------------------------
!> \brief  Convert a string to a c_char array
!!  cut string if array is too small
!!  append c_null_char if array is too large
!
	subroutine convertStringCharMiniSeed(string,array)
	character (len=*) :: string
	character (kind = c_char), dimension(:) :: array
	integer :: n,i
!
	n = min(size(array)-1,len_trim(string))    ! leave space for \0 at end of array
	if (size(array) < len_trim(string)) then
		print *,'Warning'
		print *,'convertStringCharMiniSeed: character array too short to take ',trim(string)
		print *,'Some characters will be lost.'
	endif
	array = c_null_char                        ! fill with \0
	do i=1,n
		array(i) = string(i:i)
	enddo
	end subroutine convertStringCharMiniSeed
!-----------------------------------------------------------------------
!> \brief Convert a c_char array to a string
!!  fill string up to its length
!
	subroutine convertCharStringMiniSeed(array,string,warn)
	character (kind = c_char), dimension(:) :: array
	character (len=*) :: string
	integer :: n,i,ic
	logical :: warn
!
!  fill string with blanks
!
	do i=1,len(string); string(i:i) = ' '; enddo
!
	n = min(size(array)-1,len(string))          ! assume \0 at end of string
	if (warn .and. size(array) > len(string)) then
		print *,'Warning'
		print *,'convertCharStringMiniSeed: character array too large for string'
		print *,'Some characters will be lost.'
	endif
	do i = 1,n
		ic = ichar(array(i))
		if ((ic >= 65 .and. ic <=90) .or. (ic >=48 .and. ic <= 57)) then    ! characters between A-Z or 0-9  
			string(i:i) = array(i)
		else
			exit
		endif
	enddo
	end subroutine convertCharStringMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Deallocate mini_seed
!
	subroutine deallocMiniSeed(this)
	type (mini_seed) :: this
	if (associated(this%y)) then
		if (this%link) then; nullify(this%y); else; deallocate(this%y); endif
	endif
	if (associated(this%iy)) then
		if (this%link) then; nullify(this%iy); else; deallocate(this%iy); endif
	endif
	end subroutine deallocMiniSeed
!--------------------------------------------------------------------------------------
!> \brief  Deep copy mini seed object
!> \param this Mini seed object
!> \return copy of mini seed object
!
	function deepCopyMiniSeed(this) result(copy)
	type (mini_seed) :: this,copy
!
	copy%n = this%n; copy%dt = this%dt
	copy%tstart = copyDateTime(this%tstart)
	copy%network = this%network
	copy%station = this%station
	copy%channel = this%channel
	copy%location = this%location
	allocate(copy%y(this%n))
	copy%y = this%y
	copy%link = .false.
	end function deepCopyMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get sampling interval dt
!
	double precision function dtMiniSeed(this)
	type (mini_seed), intent(in) :: this
	dtMiniSeed = this%dt
	end function dtMiniSeed
!---------------------------------------------------------------------------
!> \brief Extend a pointer array of mini_seed-objects
!
	function extendArrayMiniSeed(array,n) result(newarray)
	type (mini_seed), dimension(:), pointer :: array
	type (mini_seed), dimension(:), pointer :: newarray
	integer :: n,nold,i
!
	allocate(newarray(n))
	if (.not. associated(array)) return
	nold = min(size(array),n)
	newarray(1:nold) = array(1:nold)
	do i = 1,nold
		call unlinkDataMiniSeed(array(i))
	enddo
	do i = nold+1,size(array)
		call dealloc(array(i))
	enddo
	deallocate(array)
	end function extendArrayMiniSeed
!-------------------------------------------------------------------------------
!> \brief  Extract time window from mseed file and fill mini_seed object
!> \param this mini_seed object
!> \param filename Name of mini seed file
!> \param sdt Start time string of form yyyymmdd_hhmmss_nnnnnnnnn
!> \param edt End time string of form yyyymmdd_hhmmss_nnnnnnnnn
!> \param verbose print debugging information if true
!
	subroutine extractTimeWindowMiniSeed(this,filename,errmsg,sdt,edt,verbose)
	type (mini_seed) :: this
	character (len=*) :: filename
	character (len=*), optional :: sdt,edt
	logical, optional :: verbose
	type (error_message) :: errmsg
	logical :: debug
	real (kind = c_float), dimension(:), allocatable :: c_y
	character (kind = c_char), dimension(11) :: c_network,c_station,c_channel,c_location
	character (kind = c_char), dimension(26) :: c_sdt,c_edt
	character (kind = c_char), dimension(:), allocatable :: c_filename
	character (kind = c_char), dimension(9) :: c_sdate
	integer (kind = c_int):: c_n,c_tfs,c_tns,c_nsamp
	real (kind = c_double) :: c_dt
	type (c_ptr) :: mstup
	integer (kind = c_int) :: ierr,mstfree = 1          ! free mst after unpacking
!
	if (present(verbose)) then; debug = verbose; else; debug = .false.; endif
!
!  convert strings to c_char array
!
	allocate(c_filename(len_trim(filename)+1))
	call convertStringCharMiniSeed(filename,c_filename)
	if (present(sdt)) then; call convertStringCharMiniSeed(sdt,c_sdt); endif
	if (present(edt)) then; call convertStringCharMiniSeed(edt,c_edt); endif
!
!  read time window
!
	if (present(sdt) .and. present(edt)) then
		mstup = readTraceMSeed(c_filename,c_sdt,c_edt,c_nsamp)
		if (debug) print *,'extractTimeWindowMiniSeed: desired start time: ',trim(sdt)
		if (debug) print *,'extractTimeWindowMiniSeed: desired end   time: ',trim(edt)
	else
		mstup = readFullTraceMseed(c_filename,c_nsamp)
	endif
	if (c_nsamp > 0) then
		allocate(c_y(c_nsamp))
		ierr = unpackTraceMSeed(mstup,c_n,c_dt,c_sdate,c_tfs,c_tns,c_network,c_station,c_channel,c_location,c_y,mstfree)
		if (ierr < 0) then
			call add(errmsg,2,'Non-usable sample type','extractTimeWindowMiniSeed')
			deallocate(c_y)
			return
		endif
	else
		call add(errmsg,2,'No record in time window for file: '//trim(filename),'extractTimeWindowMiniSeed')
		deallocate(c_filename)
		return
	endif
	call createFromUnpackedMiniSeed(this,c_n,c_dt,c_sdate,c_tfs,c_tns,c_network,c_station,c_channel,c_location,c_y)
	if (debug) then
		print *,'extractTimeWindowMiniSeed: Station = ',trim(this%station),' Network = ',trim(this%network), &
	                 & ' Channel = ',this%channel,' Location = ',this%location
		print *,'extractTimeWindowMiniSeed: True start time: ',trim(convertToFullTimestringDateTime(this%tstart))
		print *,'extractTimeWindowMiniSeed: True end   time: ',trim(convertToFullTimestringDateTime(tendMiniSeed(this)))
	endif
!
	deallocate(c_filename)
	if (allocated(c_y)) deallocate(c_y)
	end subroutine extractTimeWindowMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Extract trace group from mini seed file
!> \param filename Name of mini seed input file
!> \param numtraces Number of traces read from file (output)
!> \param mstg C-Pointer to trace group
!
	subroutine extractTraceGroupMiniSeed(filename,numtraces,mstg,errmsg)
	character (len=*) :: filename
	integer :: numtraces
	type (c_ptr) :: mstg
	type (error_message) :: errmsg
	character (kind = c_char), dimension(:), allocatable :: c_filename
	integer (kind = c_int) :: ntr
!
	allocate(c_filename(len_trim(filename)+1))
	call convertStringCharMiniSeed(filename,c_filename)
	mstg = readTraceGroupMSeed(c_filename,ntr)
	numtraces = ntr
	deallocate(c_filename)
	if (numtraces == 0) then
		call add(errmsg,2,'No traces found in file '//trim(filename),'extractTraceGroupMiniSeed')
		return
	endif
	end subroutine extractTraceGroupMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Iterator through traces in TraceGroup
!
	function nextTraceInGroupMiniSeed(this,mstg,debug,errmsg) result(next)
	type (mini_seed) :: this
	type (c_ptr), value :: mstg
	logical :: debug
	type (error_message) :: errmsg
	logical :: next
	real (kind = c_float), dimension(:), allocatable :: c_y
	character (kind = c_char), dimension(11) :: c_network,c_station,c_channel,c_location
	character (kind = c_char), dimension(9) :: c_sdate
	integer (kind = c_int):: c_n,c_tfs,c_tns,ierr,c_nsamp
	real (kind = c_double) :: c_dt
	integer (kind = c_int) :: mstfree = 0              ! do not free mst after unpacking
	integer :: call_count = 0
	type (c_ptr) :: mst,mstnext
	save :: call_count,mst
!
	call_count = call_count+1
	call new(errmsg,'nextTraceInGroupMiniSeed')
	if (call_count == 1) then
		mst = getFirstTraceMSeed(mstg,c_nsamp)
	else
		mstnext = getNextTraceMSeed(mst,c_nsamp)
		mst = mstnext
	endif
	if (.not. c_associated(mst)) then
		next = .false.
		return
	endif
	next = .true.
	if (c_nsamp > 0) then
		allocate(c_y(c_nsamp))
		ierr = unpackTraceMSeed(mst,c_n,c_dt,c_sdate,c_tfs,c_tns,c_network,c_station,c_channel,c_location,c_y,mstfree)
		if (ierr < 0) then
			call new(errmsg,2,'Non-usable sample type','nextTraceInGroupMiniSeed')
			deallocate(c_y)
			return
		endif
	else
		call new(errmsg,2,'Trace record is empty','nextTraceInGroupMiniSeed')
		return
	endif
	call createFromUnpackedMiniSeed(this,c_n,c_dt,c_sdate,c_tfs,c_tns,c_network,c_station,c_channel,c_location,c_y)
	if (debug) then
		print *,'nextTraceInGroupMiniSeed: Station = ',trim(this%station),' Network = ',trim(this%network), &
	                 & ' Channel = ',this%channel,' Location = ',this%location
		print *,'nextTraceInGroupMiniSeed: True start time: ',trim(convertToFullTimestringDateTime(this%tstart))
		print *,'nextTraceInGroupMiniSeed: True end   time: ',trim(convertToFullTimestringDateTime(tendMiniSeed(this)))
	endif
!
	if (allocated(c_y)) deallocate(c_y)
	end function nextTraceInGroupMiniSeed
!--------------------------------------------------------------------------------------
!> \brief print information about miniSeed object
!
	subroutine printMiniSeed(this,lu)
	type (mini_seed), intent(in) :: this
	integer, optional :: lu
	integer :: llu
	if (present(lu)) then; llu = lu; else; llu = 6; endif
	write(lu,*) 'Net: ',this%network,' Station: ',this%station,' Channel: ',this%channel,' Location: ',this%location
	write(lu,*) 'n = ',this%n,' dt = ',this%dt
	write(lu,'(a$)') 'Start time: '; call printDateTime(this%tstart,lu)
	end subroutine
!--------------------------------------------------------------------------------------
!> \brief  Get location code
!
	character (len=2) function locationMiniSeed(this)
	type (mini_seed), intent(in) :: this
	locationMiniSeed = this%location
	end function locationMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get network
!
	character (len=2) function networkMiniSeed(this)
	type (mini_seed), intent(in) :: this
	networkMiniSeed = this%network
	end function networkMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get number of samples
!
	integer function nsampMiniSeed(this)
	type (mini_seed), intent(in) :: this
	nsampMiniSeed = this%n
	end function nsampMiniSeed
!--------------------------------------------------------------------------------------
!> \brief  Get start date in form yyyymmdd
!
	function sdateMiniSeed(this) result(sdate)
	type (mini_seed), intent(in) :: this
	character (len=8) :: sdate
	sdate = convertToSdateDateTime(this%tstart)
	end function sdateMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get station
!
	character (len=5) function stationMiniSeed(this)
	type (mini_seed), intent(in) :: this
	stationMiniSeed = this%station
	end function stationMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get start time in seconds after midnight (single)
!
	real function tanfMiniSeed(this)
	type (mini_seed), intent(in) :: this
	tanfMiniSeed = getTfsDateTime(this%tstart)+getTnsDateTime(this%tstart)*1.e-9
	end function tanfMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get start time in seconds after midnight (double)
!
	double precision function tanfDPMiniSeed(this)
	type (mini_seed), intent(in) :: this
	tanfDPMiniSeed = getTfsDateTime(this%tstart)*1.d0+getTnsDateTime(this%tstart)*1.d-9
	end function tanfDPMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get start time as date_time object
!
	function tstartMiniSeed(this) result(tstart)
	type (mini_seed), intent(in) :: this
	type (date_time) :: tstart
	tstart = this%tstart
	end function tstartMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get end time as date_time object (time of last sample)
!
	function tendMiniSeed(this) result(tend)
	type (mini_seed), intent(in) :: this
	type (date_time) :: tend
	double precision :: tadd
	integer :: dt,dtns
	tadd = (this%n-1)*this%dt
	dt = int(tadd); dtns = int((tadd-dt)*1000000000)
	tend = addDateTime(this%tstart,dt,dtns)
	end function tendMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get start time after midnight in full seconds
!
	integer function tfsMiniSeed(this)
	type (mini_seed), intent(in) :: this
	tfsMiniSeed = getTfsDateTime(this%tstart)
	end function tfsMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Get fractional part of start time in nanseconds
!
	integer function tnsMiniSeed(this)
	type (mini_seed), intent(in) :: this
	tnsMiniSeed = getTnsDateTime(this%tstart)
	end function tnsMiniSeed
!-------------------------------------------------------------------------------------
!> \brief Get pointer to data
!
	function traceMiniSeed(this)
	real, dimension(:), pointer :: traceMiniSeed
	type (mini_seed), intent(in) :: this
	traceMiniSeed => this%y
	end function traceMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Nullify pointers to  mini_seed data
!
	subroutine unlinkDataMiniSeed(this)
	type (mini_seed) :: this
	if (associated(this%y)) nullify(this%y)
	if (associated(this%iy)) nullify(this%iy)
	end subroutine unlinkDataMiniSeed
!--------------------------------------------------------------------------------------
!> \brief Write mini_seed object to file
!
	subroutine writeMiniSeed(this,filename)
	type  (mini_seed) :: this
	character (len=*) :: filename
	character (len=24) :: timestr
	character (kind = c_char), dimension(11) :: c_network,c_station,c_channel,c_location
	character (kind = c_char), dimension(133) :: c_filename
	character (kind = c_char), dimension(25) :: c_timestr
	real (kind = c_double) :: c_rate
	integer (kind = c_int) :: c_n,c_ierr
	integer (kind = c_int), dimension(:), pointer :: c_y 
!
	call convertStringCharMiniSeed(this%network,c_network)	
	call convertStringCharMiniSeed(this%station,c_station)	
	call convertStringCharMiniSeed(this%channel,c_channel)	
	call convertStringCharMiniSeed(this%location,c_location)
	call convertStringCharMiniSeed(filename,c_filename)
	timestr = convertToSeedTimestringDateTime(this%tstart)
	call convertStringCharMiniSeed(timestr,c_timestr)
	c_rate = 1./this%dt
	c_n = this%n
	c_y => this%iy
	c_ierr = writeTraceMSeed(c_filename,c_y,c_n,c_rate,c_timestr,c_network,c_station,c_channel,c_location)
	end subroutine writeMiniSeed
!
 end module
