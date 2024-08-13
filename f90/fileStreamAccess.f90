!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------
!> \brief Handle structured stream input and output (file part)

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
module fileStreamAccess
    use kindDefinitions
    implicit none
    interface dealloc
        module procedure deallocFileStreamAccess
    end interface
!
    type file_stream_access
        private
        integer :: lu                            !< Fortran identifier of file
        integer (longint) :: current_file_pos    !< current position of file pointer
        character (len=132) :: filename          !< name of file
    end type
 contains
!---------------------------------------------------------------------------
!> \brief Create a new stream access file object. Open file for writing
!> \param this file_stream_access object
!> \param lu Fortran file identifier
!> \param filename Name of stream access file
!
    integer function createFileStreamAccess(this,lu,filename) result(ios)
    type (file_stream_access) :: this
    integer :: lu
    character (len=*) :: filename
    logical :: exflag
!
    inquire(file=filename, exist = exflag)   ! delete file if it exists
    if (exflag) then
        open(lu,file=filename,access = 'stream', form='unformatted')
        close(lu,status = 'delete')
    endif
    open(lu,file=filename,form='unformatted',access='stream',status='new',iostat=ios)
    if (ios /= 0) return
    this%lu = lu
    this%current_file_pos = 1+bit_size(this%current_file_pos)/8      !  leave space for a long integer
    this%filename = filename
    end function createFileStreamAccess
!-----------------------------------------------------------------
!> \brief Open an existing stream access file for reading
!> \param this file_stream_access object
!> \param lu Fortran file identifier
!> \param filename Name of stream access file
!
    integer function openFileStreamAccess(this,lu,filename) result(ios)
    type (file_stream_access) :: this
    integer :: lu
    character (len=*) :: filename
    open(lu,file=filename,form='unformatted',access='stream',status='old',iostat=ios)
    if (ios /= 0) return
    this%lu = lu
    read(lu,pos = 1) this%current_file_pos
    end function openFileStreamAccess
!-----------------------------------------------------------------
!> \brief Deallocate a direct access file object and close file
!> \param this file_stream_access object
!
    subroutine deallocFileStreamAccess(this)
    type (file_stream_access) :: this
    close(this%lu)
    end subroutine deallocFileStreamAccess
!----------------------------------------------------------------
!> \brief Get the logical unit of the file
!> \param this file_stream_access object
!
    integer function getFileUnitStreamAccess(this)
    type (file_stream_access) :: this
    getFileUnitStreamAccess = this%lu
    end function getFileUnitStreamAccess
!----------------------------------------------------------------
!> \brief Get the current file position
!> \param this file_stream_access object
!
    function getCurrentFilePositionStreamAccess(this) result(filepos)
    type (file_stream_access) :: this
    integer (longint) :: filepos
    filepos = this%current_file_pos
    end function getCurrentFilePositionStreamAccess
!----------------------------------------------------------------
!> \brief Get the file name
!> \param this file_stream_access object
!
    function getFileNameStreamAccess(this) result(fname)
    type (file_stream_access) :: this
    character(len=132) :: fname
    fname = this%filename
    end function getFileNameStreamAccess
!----------------------------------------------------------------
!> \brief Increment current file position
!> \param this file_stream_access object
!> \param n number of bytes to move
!
    function IncrementCurrentFilePositionStreamAccess(this,n) result(filepos)
    type (file_stream_access) :: this
    integer (longint) :: filepos
    integer :: n
    this%current_file_pos = this%current_file_pos + n
    filepos = this%current_file_pos
    end function IncrementCurrentFilePositionStreamAccess
!----------------------------------------------------------------
!> \brief Set current file position
!> \param this file_stream_access object
!> \param filepos desired position of file pointer
!
    subroutine setCurrentFilePositionStreamAccess(this,filepos)
    type (file_stream_access) :: this
    integer (longint) :: filepos
    this%current_file_pos = filepos
    end subroutine setCurrentFilePositionStreamAccess
!
 end module fileStreamAccess
