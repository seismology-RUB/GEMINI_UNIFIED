!--------------------------------------------------
! \brief Module for source time function inversion
!-------------------------------------------------
 module sourceTimeFunction
	use timeSeries
	implicit none
!
 contains
!---------------------------------------------
!  Invert for source time function by fitting
!  data to Green seismograms
!  Assumes that both data and Green seismograms have same dt,
!  same length and same time points and have been filtered
!  in an identical way
!
!  dat:	array of data time series objects
!  green:	array of Green function time series objects with identical dt
!  nstf:	number of samples of source time function
!  ntap: number of samples of end of data taper
!  sigsm:	smoothing sigma
!  mag:	magnification factor for upweighting small signals
!  chi2:	a posteriori misfit
!  stf:  time series with source wavelet
!  conv: array of time series from convolution of synthetics with STF
!
	function createSourceTimeFunction(lulog,dat,syn,ntraces,nstfpos,nstfneg,ntap,sigsm,mag,chi2,conv,chitr,errmsg) result(stf)
	type (time_series) :: stf
	type (time_series), dimension(:) :: dat,syn
	integer :: ntraces,lulog
	double precision :: sigsm,mag,chi2
	type (error_message) :: errmsg
   type (time_series), dimension(:), pointer :: conv
	double precision, dimension(:,:), pointer :: bm,bmori
	double precision, dimension(:), allocatable :: d,h,work,zz,s,dori,endata,chitr
	real, dimension(:), pointer :: sigma
	integer, dimension(:), allocatable :: indx	
	double precision :: dt,dts,dtd,rnorm
	integer :: nstf,nstfpos,nstfneg,ntap,k,ka,ntot,nsamp,nsampmax,mode
	character (len=24) :: myname = 'createSourceTimeFunction'
!
!  check that data and synthetics have same amount of traces
!
	if (size(dat) /= size(syn)) then
		print *,'different number of traces in data and synthetics'
		stop
	endif
!
!  check that data traces have common dt
!
   dtd = .dt.dat(1)
   do k = 2,ntraces
      if (abs(.dt.dat(k)-dtd) > 1.e-6) then
         call add(errmsg,2,'Data traces do not have identical sampling interval',myname)
         print *,'Trace: ',k,', dt = ',.dt.dat(k),' Trace 1: dt = ',dtd
         return
      endif
   enddo
!
!  check that synthetic traces have common dt
!
   dts = .dt.syn(1)
   do k = 2,ntraces
      if (abs(.dt.syn(k)-dts) > 1.e-6) then
         call add(errmsg,2,'Synthetic traces do not have identical sampling interval',myname)
         return
      endif
   enddo
!
   nstf = nstfpos+nstfneg
!
!  check that data and syn dt, nsamp, tanf is equal
!
   ntot = 0
   nsampmax = 0
   do k = 1,ntraces
      if (abs(.dt.dat(k)-.dt.syn(k)) > 1.e-6) then
         print *,'Trace: ',k
         call add(errmsg,2,'Data and Synthetic traces do not have identical sampling interval',myname)
         return
      endif
      if (abs(.tanfdp.dat(k)-.tanfdp.syn(k)) > 1.e-6) then
         print *,'Trace: ',k
         call add(errmsg,2,'Data and Synthetic traces do not have identical start time',myname)
         return
      endif
      if (.nsamp.dat(k) /= .nsamp.syn(k)) then
         print *,'Trace: ',k
         call add(errmsg,2,'Data and Synthetic traces do not have identical number of samples',myname)
         return
      endif
      if (.nsamp.dat(k) < nstf) then
         print *,'Trace: ',k,.nsamp.dat(k),nstf
         call add(errmsg,2,'Data traces shorter than source time function',myname)
         return
      endif
      ntot = ntot+.nsamp.dat(k)
      nsampmax = max(nsampmax,.nsamp.dat(k))
   enddo
!
   dt = dtd
!
!  calculate sigmas, normalize data and calculate system matrix
!
	allocate(bm(ntot+nstf,nstf))              ! space for system matrix
	allocate(bmori(ntot,nstf))                ! space for system matrix without smoothing constraints
	allocate(d(ntot+nstf))                    ! space for sytem data vector
	allocate(dori(ntot))                      ! space for original data vector without appended zeros
   allocate(s(nsampmax))                     ! space for convolved synthetics
   allocate(conv(ntraces))                   ! space for convolved time series
   allocate(endata(ntraces))                 ! space for energy of data
	endata = 0.
	bm = 0.d0
	d = 0.d0
	ka = 0
	do k = 1,ntraces
      nsamp = .nsamp.dat(k)
		sigma => sigmaEnvelopeTimeSeries(dat(k),real(mag),5)            ! nexp = 5
		call matrixDataSourceTimeFunction(dat(k),syn(k),sigma,nstf,nstfneg,ntap,d(ka+1:ka+nsamp),bm(ka+1:ka+nsamp,1:nstf))
		endata(k) = energyTimeSeries(dat(k),sigma)
		ka = ka+nsamp
      deallocate(sigma)
	enddo
	write(lulog,*) 'Total number of data samples: ',ntot
	write(lulog,*) 'Samples of source time function, positive side, negative side: ',nstf,nstfpos,nstfneg
	write(lulog,*) 'Number of samples of end of data taper: ',ntap
   write(lulog,*) 'Maximum of BM-matrix before smoothing added: ',maxval(bm)
!
!  create a copy of data vector and system matrix before adding zeros and smoothing constraints, respectively
!  for reconstruction of data by convolution of synthetics and Stf and for misfit calculation
!
   bmori = bm(1:ntot,1:nstf)
   dori = d(1:ntot)
!
!  append smoothing matrix and zeros to data vector
!
	call smoothingMatrixSourceTimeFunction(nstf,ntot,sigsm,bm(ntot+1:ntot+nstf,1:nstf))
   write(lulog,*) 'Maximum of BM-matrix after smoothing added: ',maxval(bm)
!
!  call non-negative least squares routine
!  bm and d are modified after return from nnls to QA and Qd (Q orthogonal matrix)
!
	allocate(h(nstf),work(nstf),indx(nstf),zz(ntot+nstf))
	call nnls(bm,ntot+nstf,ntot+nstf,nstf,d,h,rnorm,work,zz,indx,mode)
	if (mode > 1) then
		print *,'Non-negative least square routine not successful'
		print *,'mode = ',mode
	endif
!
!  calculate new misfit over energy (endata contains a factor dt and also 1/sigma)
!  fill result into time series objects
!
	call createFromDataTimeSeries(stf,nstf,-nstfneg*dt,dt,real(h))
	chi2 = 0.d0
	allocate(chitr(ntraces))
	ka = 0
	do k = 1,ntraces
      nsamp = .nsamp.dat(k)
      sigma => sigmaEnvelopeTimeSeries(dat(k),real(mag),5)            ! nexp = 5
	   s = matmul(bmori(ka+1:ka+nsamp,:),h)
	   chitr(k) = dt*sum((dori(ka+1:ka+nsamp)-s(1:nsamp))**2)/endata(k)
	   chi2 = chi2+chitr(k)*endata(k)
      call createFromDataTimeSeries(conv(k),nsamp,.tanfdp.dat(k),dt,real(s*sigma))
      ka = ka+nsamp
      deallocate(sigma)
   enddo
   chi2 = chi2/sum(endata)
!
	deallocate(d,dori,s,bm,bmori,h,work,indx,zz,endata)
	end function createSourceTimeFunction
!------------------------------------------------------------------------------
!  set up matrix connecting data and source time function: d = B h made out of
!  samples of a Green seismogram just for one trace
!  see noncausfalt.txt for explanation
!
	subroutine matrixDataSourceTimeFunction(dat,green,sigma,nstf,nhm,ntap,d,bm)
	type (time_series) :: dat,green
	real, dimension(:) :: sigma
	double precision, dimension(:,:) :: bm
	double precision, dimension(:) :: d
	real, dimension(:), pointer :: g
	real :: dt,tanfd,tanfg,tapval,taparg
	integer :: nstf,nhm,ntap,k,n,ng,nd,j
!
!  data vector
!
	nd = .nsamp.dat
	tanfd = .tanf.dat
	d(1:nd) = .trace.dat/sigma
!
!  part of matrix with green samples
!  accounts for the fact that data may start later than synthetics
!
	g => .trace.green
	ng = .nsamp.green
	dt = .dt.green
	tanfg = .tanf.green
!
	do k = 1,nd                                               ! Loop over samples of data
		tapval = 1.0
		if (k > nd-ntap) then                                  ! apply taper at end of data prediction
			taparg = real(k-nd+ntap)/real(ntap)*0.5*mc_pi
			tapval = cos(taparg)**2
		end if
		do n = max(1,k-nd+nhm+1),min(k+nhm,nstf)               ! Loop over samples of stf
		   j = n-nhm
			bm(k,n) = dt*g(k-j+1)/sigma(k)*tapval
		enddo
	enddo
	end subroutine matrixDataSourceTimeFunction
!--------------------------------------------------------------------------------
!  smoothing part of system matrix
!
	subroutine smoothingMatrixSourceTimeFunction(nstf,ntot,sigsm,sm)
	integer :: nstf,l,ntot
	double precision :: sigsm,sigeff
	double precision, dimension(:,:) :: sm
	sigeff = sigsm*sqrt(real(ntot))
	sm(1,1) = (2.d0+20.0)/sigeff
	sm(1,2) = -1.d0/sigeff
	do l = 2,nstf-1
		sm(l,l)=(2.d0+0.02)/sigeff
		sm(l,l-1)=-1.d0/sigeff
		sm(l,l+1)=-1.d0/sigeff
	enddo
	sm(nstf,nstf) = (2.d0+20.0)/sigeff
	sm(nstf,nstf-1) = -1.d0/sigeff
	end subroutine smoothingMatrixSourceTimeFunction
!----------------------------------------------------------------------------
!  write source time function to ascii file
!
	subroutine writeAsciiSourceTimeFunction(lu,filename,stf)
	type (time_series) :: stf
	integer :: lu,i,nstf
	character (len=*) :: filename
	real, dimension(:), pointer :: s
	real :: dt
!
	dt = .dt.stf; nstf = .nsamp.stf
	s => traceTimeSeries(stf)
	open(lu,file=filename)
	write(lu,*) nstf,dt
	do i=1,nstf
		write(lu,'(f15.6,e20.10)') (i-1)*dt,s(i)
	enddo
	close(lu)
	end subroutine writeAsciiSourceTimeFunction
!----------------------------------------------------------------------------
!  calculate moment magnitude time series
!
	function momentMagnitudeSourceTimeFunction(sumstf,rmsc) result(mw)
	type (time_series) :: sumstf
	real :: rmsc
	real, dimension(:), pointer :: mw
!
	allocate(mw(.nsamp.sumstf))
	mw = 0.0
	where (.trace.sumstf > 0.) mw = 0.6666666*alog10((.trace.sumstf)*rmsc*1.e7)-10.7
	end function momentMagnitudeSourceTimeFunction
!
 end module sourceTimeFunction
