! ====================================================================
!  Control parameters for Bulirsch-Stoer integration engine
! ====================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   along with GEMIN_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------
!  Defines controlling parameters for the  Bulirsch-Stoer integration engine.
!  Also sets up scheme for step-size subdivision and Richardson extrapolation.
!-------------------------------------------------------------------------------
module parametersBulirschStep
    implicit none
!-------------------------------------------------------------------------------------------------------------
    double precision, parameter :: grow_bulirsch_step = 1.25   ! growth factor for step size (shrink*growth=1)
    double precision, parameter :: shrink_bulirsch_step = 0.8  ! shrink factor for step size (shrink*growth=1)
    integer, parameter :: secrecy_bulirsch_step = 0            ! print screen output if secrecy <= this value
    integer, parameter :: maxdiv_bulirsch_step = 24            ! smallest considered step is h/maxsub
    integer, parameter :: maxtry_bulirsch_step = 7             ! maximum number of subdivision attempts
    integer, parameter :: mintry_bulirsch_step = 3             ! minimum required number of subdivision attempts
    integer, parameter :: opttry_bulirsch_step = 6             ! optimal number of subdivision attempts
    integer, parameter :: kgv_bulirsch_step = 48               ! smallest common multiple of the subdivisions
    integer, dimension(7), parameter :: &
         seqdiv_bulirsch_step = (/2,4,6,8,12,16,24/)           ! subdivision sequence (only valid for maxtry=7 and maxdiv=24)
    integer, dimension(8,7), parameter :: &
         ixdiv_bulirsch_step = reshape( &                      ! only valid for maxtry=7 and maxdiv=24
         (/ 0,24,48, 0, 0, 0, 0, 0, &                          ! = m*48/nstep for m = 0,nstep
           12,36, 0, 0, 0, 0, 0, 0, &                          ! fractions of interval h expressed in h/48
            8,16,32,40, 0, 0, 0, 0, &                          ! at which new evaluation of system matrix
            6,18,30,42, 0, 0, 0, 0, &                          ! is needed. Zeros are just a fill in for
            4,20,28,44, 0, 0, 0, 0, &                          ! not required values
            3, 9,15,21,27,33,39,45, &
            2,10,14,22,26,34,38,46 /),(/8,7/))
    integer, dimension(7), parameter :: &                      ! number of new system matrix evaluations required
         nxdiv_bulirsch_step = (/3,2,4,4,4,8,8/)               ! per subdivision attempt
end module
