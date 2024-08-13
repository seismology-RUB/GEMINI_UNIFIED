/*----------------------------------------------------------------------------
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
 ---------------------------------------------------*/
/* these are mimicking functions to implement sunfortran routines */

#include <stddef.h>
#include <time.h>

void idate_(int *iarray);
void itime_(int *iarray);
/*
 *--------------------------------------
 *  function mimicking fortran idate
 *--------------------------------------
*/
void idate_(iarray)
    int *iarray;
{
    struct tm *datime;
    time_t zeit;
    zeit = time((time_t *) NULL);
    datime = localtime(&zeit);
    iarray[0] = datime->tm_mday;
    iarray[1] = datime->tm_mon + 1;
    iarray[2] = datime->tm_year;
}
/*
 *--------------------------------------
 *  function mimicking fortran itime
 *--------------------------------------
*/
void itime_(iarray)
    int *iarray;
{
    struct tm *datime;
    time_t zeit;
    zeit = time((time_t *) NULL);
    datime = localtime(&zeit);
    iarray[0] = datime->tm_hour;
    iarray[1] = datime->tm_min;
    iarray[2] = datime->tm_sec;
}
