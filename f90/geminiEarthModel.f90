! ==============================================================================
!  Earth model module for use with kernel computation
! ==============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 is free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 is distributed
!   in the hope that they will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------
!   Deals with earth model information at external radial nodes.
!   Needed for kernel computation.
!-----------------------------------------------------------------------------
module geminiEarthModel
    use errorMessage
    use locatePoint
    use string
    implicit none
    interface dealloc; module procedure deallocGeminiEarthModel; end interface
    type gemini_earth_model
       character (len = max_length_string) :: name              !< name of earth model
       character (len = max_length_string) :: attmode           !< attenuation mode string
       integer :: aniflag                                       !< isotropic (0), transversely isotropic (1)
       double precision :: fref                                 !< reference frequency
       double precision :: rearth                               !< earth radius
       integer :: nnod                                          !< number of nodes
       real, dimension(:), allocatable :: rnod                  !< radial nodes in kilometers
       complex, dimension(:,:), allocatable :: zelcon           !< complex A,C,F,L,N,Kappa,Mue
       real, dimension(:), allocatable :: rho                   !< density
       real, dimension(:), allocatable :: qkinv                 !< inverse Qkappa, qk = 1/Qk
       real, dimension(:), allocatable :: qminv                 !< inverse Qmu, qm = 1/Qm
    end type gemini_earth_model
!
 contains
!----------------------------------------------------------------------------
!  Read gemini earth model from ASCII file
!  Used by computeKernelWavefield and computeKernelGreenTensor
!
    subroutine readGeminiEarthModel(this,lu,filename,errmsg)
    type (gemini_earth_model) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len=20) :: myname = 'readGeminiEarthModel'
    integer :: i,j,ierr,nnod
    real, dimension(14) :: zelcon_tmp
!
    call addTrace(errmsg,myname)
    open(lu,file = filename, status = 'old', iostat = ierr)
    if (ierr /= 0) then
        call add(errmsg,2,'cannot open '+filename,myname)
        return
    endif
    read(lu,'(a)') this%name
    read(lu,'(a)') this%attmode
    read(lu,*) this%aniflag,this%fref,this%rearth
    read(lu,*) nnod
    this%nnod = nnod
!
!  allocate space for parameters and nodes and read
!    
    allocate(this%rnod(nnod),this%rho(nnod),this%qkinv(nnod),this%qminv(nnod))
    allocate(this%zelcon(nnod,7))
    do i = 1,nnod
       read(lu,*) this%rnod(i),this%rho(i),this%qkinv(i),this%qminv(i)
       read(lu,*) zelcon_tmp(:) !(this%zelcon(i,j),j = 1,7)
       do j = 1,7
          this%zelcon(i,j) = cmplx(zelcon_tmp(2*j-1),zelcon_tmp(2*j))
       end do
    enddo
    close(lu)
    end subroutine readGeminiEarthModel
!--------------------------------------------------------------------------------------
!  Write Gemini earth model adding final line with number of lateral grid nodes
!  Used by computeKernelWavefield and computeKernelGreenTensor
!
    subroutine writeGeminiEarthModel(this,lu,filename,ng,errmsg)
    type (gemini_earth_model) :: this
    integer :: lu
    character (len=*) :: filename
    integer :: ng
    type (error_message) :: errmsg
    integer :: i,j,ierr
    character (len=21) :: myname = 'writeGeminiEarthModel'
!
    call addTrace(errmsg,myname)
    open(lu,file = filename,iostat = ierr)
    if (ierr /= 0) then
        call add(errmsg,2,'cannot open '+filename,myname)
        return
    endif
    write(lu,'(a)') trim(this%name)
    write(lu,'(a)') trim(this%attmode)
    write(lu,'(i4,2e15.6)') this%aniflag,this%fref,this%rearth
    write(lu,'(i5)') this%nnod
    do i = 1,this%nnod
       write(lu,'(4e15.6)') this%rnod(i),this%rho(i),this%qkinv(i),this%qminv(i)
       write(lu,'(14e15.6)') (this%zelcon(i,j),j = 1,7)
    enddo
    write(lu,'(i7)') ng
    close(lu)
    end subroutine writeGeminiEarthModel
!--------------------------------------------------------------------------------------------
! Deallocate gemini_earth_model
!
    subroutine deallocGeminiEarthModel(this)
    type (gemini_earth_model) :: this
    if (allocated(this%rnod)) deallocate(this%rnod)
    if (allocated(this%rho)) deallocate(this%rho)
    if (allocated(this%qkinv)) deallocate(this%qkinv)
    if (allocated(this%qminv)) deallocate(this%qminv)
    if (allocated(this%zelcon)) deallocate(this%zelcon)
    end subroutine deallocGeminiEarthModel
!----------------------------------------------------------------------------
!> \brief Get index of node closest to given depth in km
!
    function getNodeIndexFromDepthGeminiEarthModel(this,zs) result(jr)
    type (gemini_earth_model) :: this
    real :: zs,rs
    integer :: jr
!
    rs = this%rearth-zs
    jr = locate(rs,this%nnod,this%rnod)              ! index of node below rs
    if (jr == 0) jr = 1                              ! use deepest node if rs is below it
    if (jr < this%nnod) then                         ! take node closest to rs
       if (abs(rs-this%rnod(jr)) > abs(rs-this%rnod(jr+1))) jr = jr+1
    endif
    end function getNodeIndexFromDepthGeminiEarthModel
 end module
