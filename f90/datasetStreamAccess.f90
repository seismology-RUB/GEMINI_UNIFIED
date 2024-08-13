!--------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini Unified
!
!   Gemini Unified is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini Unified are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini Unified.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------
!> \brief Handle structured stream input and output (dataset part)

!> \author Wolfgang Friederich

!> \par Description
!>  Support for stream access file reading and writing.
!!  File tree starts with a root group that may contain further groups
!!  and datasets. Each subgroup can contain subgroups and datasets
!!  itself. Idea is to create groups and datasets and to build up a
!!  tree under the root group which contains links to all subgroups
!!  and datasets. Datasets are written to file, file positions are stored
!!  in the group and dataset objects. Finally, the complete tree
!!  with information about the data is also written to file in a way
!!  that by reading the file the complete information and data can
!!  be retrieved.
!<---------------------------------------------------------------
module datasetStreamAccess
    use fileStreamAccess
    use flexibleType
    use kindDefinitions
    implicit none
    interface new
        module procedure createDatasetStreamAccess
    end interface
    interface dealloc
        module procedure fullyDeallocDatasetStreamAccess
    end interface
    interface writeDatasetVectorStreamAccess
        module procedure writeDatasetIntegerVectorStreamAccess
        module procedure writeDatasetRealVectorStreamAccess
        module procedure writeDatasetDoubleVectorStreamAccess
        module procedure writeDatasetComplexVectorStreamAccess
        module procedure writeDatasetDoubleComplexVectorStreamAccess
        module procedure writeDatasetCharVectorStreamAccess
        module procedure writeDatasetFlexibleVectorStreamAccess
    end interface
    interface writeDataset2DArrayStreamAccess
        module procedure writeDatasetReal2DArrayStreamAccess
        module procedure writeDatasetDouble2DArrayStreamAccess
        module procedure writeDatasetComplex2DArrayStreamAccess
        module procedure writeDatasetDoubleComplex2DArrayStreamAccess
    end interface
    interface readDatasetVectorStreamAccess
        module procedure readDatasetIntegerVectorStreamAccess
        module procedure readDatasetRealVectorStreamAccess
        module procedure readDatasetDoubleVectorStreamAccess
        module procedure readDatasetComplexVectorStreamAccess
        module procedure readDatasetDoubleComplexVectorStreamAccess
        module procedure readDatasetCharVectorStreamAccess
        module procedure readDatasetFlexibleVectorStreamAccess
    end interface
    interface readDataset2DArrayStreamAccess
        module procedure readDatasetReal2DArrayStreamAccess
        module procedure readDatasetDouble2DArrayStreamAccess
        module procedure readDatasetComplex2DArrayStreamAccess
        module procedure readDatasetDoubleComplex2DArrayStreamAccess
    end interface
    type data_stream_access
        integer (longint) :: filepos_start                 !< position where data start
        integer :: rank                                    !< dimensionality of data array
        integer, dimension(:), pointer :: dims =>  null()  !< number of data
        integer :: datatype                                !< integer coded primitive type of data
    end type
!
 contains
!----------------------------------------------------------------------
!> \brief Create a dataset
!> \param this a data_stream_access object
!> \param rank rank of data array
!> \param dims array of integers specifiying size in each rank
!> \param datatype type of data in array using conventions in \link primitiveTypeEncoding.f90 primitiveTypeEncoding \endlink
!> \par
!> Specify rank, dims and datatype but nothing else
!! which is done later on writing to file
!<
    subroutine createDatasetStreamAccess(this,rank,dims,datatype)
    type (data_stream_access) :: this
    integer :: rank
    integer, dimension(:) :: dims
    integer :: datatype
!
    this%rank = rank
    this%datatype = datatype
    allocate(this%dims(rank))
    this%dims = dims
    end subroutine createDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Fully deallocate a dataset
!> \param this data_stream_access object
!
    subroutine fullyDeallocDatasetStreamAccess(this)
    type (data_stream_access) :: this
    if (associated(this%dims)) deallocate(this%dims)
    end subroutine fullyDeallocDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Nullify the pointers contained in a dataset object
!> \param this data_stream_access object
!
    subroutine nullifyDatasetStreamAccess(this)
    type (data_stream_access) :: this
    if (associated(this%dims)) nullify(this%dims)
    end subroutine nullifyDatasetStreamAccess
!----------------------------------------------------------------------
!> \brief Reallocate a dataset pointer array
!> \param p pointer to an array of datasets
!> \param n number of dataset objects to be allocated
!> \return newp pointer to the new array of datasets
!> \par
!> The contents of the old array \a p is copied into the new one.
!! Afterwards, the old array \a p is deallocated.
!<
!--------------------------------------------------------------------
    function reallocateDatasetStreamAccess(p, n) result(newp)
    type (data_stream_access), pointer, dimension(:) :: p, newp
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(newp(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    newp(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocateDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Read info part of dataset (without data)
!
    subroutine readInfoDataStreamAccess(this,fda)
    type (data_stream_access) :: this
    type (file_stream_access) :: fda
    integer :: lu
    integer (longint) :: filepos
!
    lu = getFileUnitStreamAccess(fda)
!    filepos = getCurrentFilePositionStreamAccess(fda)
!    this%filepos_start = filepos
    filepos = this%filepos_start
    read(lu,pos = filepos) this%rank                  ! read dataset info
    allocate(this%dims(this%rank))
    read(lu,pos = filepos+kind(1)) this%dims,this%datatype
    call setCurrentFilePositionStreamAccess(fda,filepos+kind(1)*(1+this%rank+1))
!    filepos = IncrementCurrentFilePositionStreamAccess(fda,kind(1)*(1+this%rank+1))
    end subroutine readInfoDataStreamAccess
!----------------------------------------------------------------------
!> \brief Return size of data vector in data set
!
    function getVectorSizeDataStreamAccess(this) result(n)
    type (data_stream_access) :: this
    integer :: n
    n = this%dims(1)
    end function getVectorSizeDataStreamAccess
!----------------------------------------------------------------------
!> \brief  Write an integer vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d integer data array
!
    subroutine writeDatasetIntegerVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    integer, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1+size(d))*kind(1))
    end subroutine writeDatasetIntegerVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a real vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d real data array
!
    subroutine writeDatasetRealVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.0))
    end subroutine writeDatasetRealVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data array
!
    subroutine writeDatasetDoubleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.d0))
    end subroutine writeDatasetDoubleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a complex vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d complex data array
!
    subroutine writeDatasetComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.0))
    end subroutine writeDatasetComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double complex vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double complex data array
!
    subroutine writeDatasetDoubleComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.d0))
    end subroutine writeDatasetDoubleComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a character vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d character data array
!
    subroutine writeDatasetCharVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    character(len=*), dimension(:) :: d
    integer :: clen
    integer :: lu
    integer (longint) :: filepos
!
    clen = len(d(1))
    if (clen > 80) then
        print *,'writeDatasetCharVectorStreamAccess: string length greater 80 !'
        stop
    endif
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,clen,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+2)*kind(1)+size(d)*clen)
    end subroutine writeDatasetCharVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a flexible type vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d flexible data array
!
    subroutine writeDatasetFlexibleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    type (flexible), dimension(:) :: d
    integer :: lu,nbytes,j
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,T_FLEXIBLE
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    do j=1,size(d)
        call writeSAFlexibleType(d(j),lu,filepos,nbytes)
        filepos = IncrementCurrentFilePositionStreamAccess(fda,nbytes)
    enddo
    end subroutine writeDatasetFlexibleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a real 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d real data 2D-array
!
    subroutine writeDatasetReal2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.0))
    end subroutine writeDatasetReal2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data 2D-array
!
    subroutine writeDatasetDouble2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.d0))
    end subroutine writeDatasetDouble2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Write a complex 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d complex data 2D-array
!
    subroutine writeDatasetComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.0))
    end subroutine writeDatasetComplex2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double complex 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data 2D-array
!
    subroutine writeDatasetDoubleComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.d0))
    end subroutine writeDatasetDoubleComplex2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read an integer vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d integer data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetIntegerVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    integer, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetIntegerVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_INTEGER) &
        print *, "WARNING in readDatasetIntegerVectorStreamAccess: this dataset is NOT ",&
        "of type integer, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1))
    end subroutine readDatasetIntegerVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a real vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d real data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetRealVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetRealVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_REAL) &
        print *, "WARNING in readDatasetRealVectorStreamAccess: this dataset is NOT ",&
        "of type real, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.0))
    end subroutine readDatasetRealVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDoubleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDoubleVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE) &
        print *, "WARNING in readDatasetDoubleVectorStreamAccess: this dataset is NOT ",&
        "of type double, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.d0))
    end subroutine readDatasetDoubleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a complex vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d complex data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetComplexVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_COMPLEX) &
        print *, "WARNING in readDatasetComplexVectorStreamAccess: this dataset is NOT ",&
        "of type complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.0))
    end subroutine readDatasetComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double complex vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double complex data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDoubleComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDoubleComplexVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE_COMPLEX) &
        print *, "WARNING in readDatasetDoubleComplexVectorStreamAccess: this dataset is NOT ",&
        "of type double complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.d0))
    end subroutine readDatasetDoubleComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a character vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d character data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetCharVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    character (len=80), dimension(:), pointer :: d
    integer :: clen,i,lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype,clen
    else
        read(lu,pos=filepos) rank,dims,datatype,clen
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetCharVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_CHAR) &
        print *, "WARNING in readDatasetCharVectorStreamAccess: this dataset is NOT ",&
        "of type character, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+2)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) (d(i)(1:clen),i=1,size(d))
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*clen)
    end subroutine readDatasetCharVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a flexible type vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d flexible data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetFlexibleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    type (flexible), dimension(:), pointer :: d
    integer :: lu,nbytes,j
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetFlexibleVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_FLEXIBLE) &
        print *, "WARNING in readDatasetFlexibleVectorStreamAccess: this dataset is NOT ",&
        "of flexible type, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    do j=1,size(d)
        call readSAFlexibleType(d(j),lu,filepos,nbytes)
        filepos = IncrementCurrentFilePositionStreamAccess(fda,nbytes)
    enddo
    end subroutine readDatasetFlexibleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a real 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetReal2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetReal2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_REAL) &
        print *, "WARNING in readDatasetReal2DArrayStreamAccess: this dataset is NOT ",&
        "of type real, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.0))
    end subroutine readDatasetReal2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDouble2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDouble2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE) &
        print *, "WARNING in readDatasetDouble2DArrayStreamAccess: this dataset is NOT ",&
        "of type double, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.d0))
    end subroutine readDatasetDouble2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read a complex 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d complex data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetComplex2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_COMPLEX) &
        print *, "WARNING in readDatasetComplex2DArrayStreamAccess: this dataset is NOT ",&
        "of type complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.0))
    end subroutine readDatasetComplex2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double complex 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDoubleComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDoubleComplex2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE_COMPLEX) &
        print *, "WARNING in readDatasetDoubleComplex2DArrayStreamAccess: this dataset is NOT ",&
        "of type double complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.d0))
    end subroutine readDatasetDoubleComplex2DArrayStreamAccess
!--------------------------------------------------------------------------------
!> \brief Print Info about a dataset
!
    subroutine printDsetInfoStreamAccess(this,idx)
    type (data_stream_access) :: this
    integer :: idx
    character(len=29) :: fmt
    write(fmt,'(a25,i1,a3)') '(a2,i7,a2,i5,a2,i5,a2,1x,',this%rank,'i5)'
    write(*,fmt) '| ',idx,' |',this%rank,' |',this%datatype,' |',this%dims
    end subroutine printDsetInfoStreamAccess
!
 end module
