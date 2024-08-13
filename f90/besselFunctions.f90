! ==============================================================================
!  Compute Bessel and Neumann functions of zeroth and first order
! ==============================================================================
!----------------------------------------------------------------------------
!   The code below was taken from Netlib (www.netlib.org)
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
!   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------
module besselFunctions
   implicit none
contains
!--------------------------------------------------------------------------
!  Compute values of Bessel functions at wavenumber points up to order 2
!
!  nwn:         number of wabenumbers
!  dwn:         wavenumber spacing (in 1/m)
!  x:           epicentral distance in m
!  besselj:     array with Bessel functions
!
    subroutine computeBesselFunctions(nwn,dwn,x,besselj)
    integer :: nwn
    double precision :: dwn,x
    double precision, dimension(:,:), allocatable :: besselj
    double precision :: wn
    integer :: j
!
    allocate(besselj(nwn,0:2))
    besselj(1,0) = 1.d0
    besselj(1,1) = 0.d0
    besselj(1,2) = 0.d0
    wn = dwn
    do j = 2,nwn
        besselj(j,0) = dbessj0(wn*x)
        besselj(j,1) = dbessj1(wn*x)
        besselj(j,2) = -besselj(j,0)+2./(wn*x)*besselj(j,1)
        wn = wn+dwn
    enddo
    end subroutine computeBesselFunctions
!--------------------------------------------------------------------
!  Bessel function of zeroth order for |x| <= XMAX (see CALJY0)
!
    REAL FUNCTION BESSJ0(X)
    INTEGER JINT
    REAL  X, RESULT
!
    JINT=0
    CALL CALJY0(X,RESULT,JINT)
    BESSJ0 = RESULT
    END FUNCTION BESSJ0
!--------------------------------------------------------------------
!  Neumann function of zeroth order for 0 < |x| <= XMAX (see CALJY0)
!    
    REAL FUNCTION BESSY0(X)
    INTEGER JINT
    REAL  X, RESULT
!
    JINT=1
    CALL CALJY0(X,RESULT,JINT)
    BESSY0 = RESULT
    END FUNCTION BESSY0
!--------------------------------------------------------------------
!  Bessel function of first order for |x| <= XMAX (see CALJY1)
!
    FUNCTION BESSJ1(X)
    INTEGER JINT
    REAL BESSJ1,RESULT,X
!
    JINT=0
    CALL CALJY1(X,RESULT,JINT)
    BESSJ1 = RESULT
    END FUNCTION BESSJ1
!--------------------------------------------------------------------
!  Neumann function of first order for 0 < |x| <= XMAX (see CALJY0)
!    
    FUNCTION BESSY1(X)
    INTEGER JINT
    REAL BESSY1,RESULT,X
!
    JINT=1
    CALL CALJY1(X,RESULT,JINT)
    BESSY1 = RESULT
    END FUNCTION BESSY1
!--------------------------------------------------------------------
!   Subroutine CALJY0 really doing work
!
!   This routine computes zero-order Bessel functions of the first and
!   second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!   for Y0, and |X| <= XMAX for J0. This routine is intended for internal
!   use only, all computations being concentrated in
!   this one routine.
!
!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the 
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESSJ0(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1. 
!   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).
!
!     Approximate values for some important machines are
!
!                          eps      XMAX     XSMALL      XINF  
!
!  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
!
!*******************************************************************
!*******************************************************************
!
! Error Returns
!
!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY0 is called with a negative or zero argument.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, COS, LOG, SIN, SQRT
!
!
!  Latest modification: June 2, 1989
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division 
!          Argonne National Laboratory
!          Argonne, IL 60439
!--------------------------------------------------------------------
    SUBROUTINE CALJY0(ARG,RESULT,JINT)
    INTEGER :: I,JINT
    REAL :: ARG,AX,CONS,DOWN,EIGHT,FIVE5,FOUR,ONE,ONEOV8,PI2,PJ0,   &
            PJ1,PLG,PROD,PY0,PY1,PY2,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1, &
            QY2,Q0,Q1,RESJ,RESULT,R0,R1,SIXTY4,THREE,TWOPI,TWOPI1,  &
            TWOPI2,TWO56,UP,W,WSQ,XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,   &
            XJ1,XJ01,XJ02,XJ11,XJ12,XY,XY0,XY01,XY02,XY1,XY11,XY12, &
            XY2,XY21,XY22,Z,ZERO,ZSQ
    DIMENSION :: PJ0(7),PJ1(8),PLG(4),PY0(6),PY1(7),PY2(8),P0(6),P1(6), &
            QJ0(5),QJ1(7),QLG(4),QY0(5),QY1(6),QY2(7),Q0(5),Q1(5)
!-------------------------------------------------------------------
!  Mathematical constants
!    CONS = ln(.5) + Euler's gamma
!-------------------------------------------------------------------
    DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0E0,1.0E0,3.0E0,4.0E0,8.0E0/,  &
         FIVE5,SIXTY4,ONEOV8,P17/5.5E0,64.0E0,0.125E0,1.716E-1/,  &
         TWO56,CONS/256.0E0,-1.1593151565841244881E-1/,           &
         PI2,TWOPI/6.3661977236758134308E-1,6.2831853071795864769E0/, &
         TWOPI1,TWOPI2/6.28125E0,1.9353071795864769253E-3/
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
      DATA XMAX/1.34E+08/,XSMALL/1.22E-09/,XINF/1.7E+38/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
      DATA XJ0/2.4048255576957727686E+0/,XJ1/5.5200781102863106496E+0/, &
           XY0/8.9357696627916752158E-1/,XY1/3.9576784193148578684E+0/, &
           XY2/7.0860510603017726976E+0/, &
           XJ01/ 616.0E+0/, XJ02/-1.4244423042272313784E-03/, &
           XJ11/1413.0E+0/, XJ12/ 5.4686028631064959660E-04/, &
           XY01/ 228.0E+0/, XY02/ 2.9519662791675215849E-03/, &
           XY11/1013.0E+0/, XY12/ 6.4716931485786837568E-04/, &
           XY21/1814.0E+0/, XY22/ 1.1356030177269762362E-04/
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
      DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02, &
               -5.4989956895857911039E+02,3.5687548468071500413E+02/ 
      DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02, &
               -3.3442903192607538956E+02,1.7843774234035750207E+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
      DATA PJ0/6.6302997904833794242E+06,-6.2140700423540120665E+08, &
               2.7282507878605942706E+10,-4.1298668500990866786E+11, &
              -1.2117036164593528341E-01, 1.0344222815443188943E+02, &
              -3.6629814655107086448E+04/
      DATA QJ0/4.5612696224219938200E+05, 1.3985097372263433271E+08, &
               2.6328198300859648632E+10, 2.3883787996332290397E+12, &
               9.3614022392337710626E+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
      DATA PJ1/4.4176707025325087628E+03, 1.1725046279757103576E+04, &
               1.0341910641583726701E+04,-7.2879702464464618998E+03, &
              -1.2254078161378989535E+04,-1.8319397969392084011E+03, &
               4.8591703355916499363E+01, 7.4321196680624245801E+02/
      DATA QJ1/3.3307310774649071172E+02,-2.9458766545509337327E+03, &
               1.8680990008359188352E+04,-8.4055062591169562211E+04, &
               2.4599102262586308984E+05,-3.5783478026152301072E+05, &
              -2.5258076240801555057E+01/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  3.0
!--------------------------------------------------------------------
      DATA PY0/1.0102532948020907590E+04,-2.1287548474401797963E+06, &
               2.0422274357376619816E+08,-8.3716255451260504098E+09, &
               1.0723538782003176831E+11,-1.8402381979244993524E+01/
      DATA QY0/6.6475986689240190091E+02, 2.3889393209447253406E+05, &
               5.5662956624278251596E+07, 8.1617187777290363573E+09, &
               5.8873865738997033405E+11/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
!        3.0  <  |X|  <=  5.5
!--------------------------------------------------------------------
      DATA PY1/-1.4566865832663635920E+04, 4.6905288611678631510E+06, &
               -6.9590439394619619534E+08, 4.3600098638603061642E+10, &
               -5.5107435206722644429E+11,-2.2213976967566192242E+13, &
                1.7427031242901594547E+01/
      DATA QY1/ 8.3030857612070288823E+02, 4.0669982352539552018E+05, &
                1.3960202770986831075E+08, 3.4015103849971240096E+10, &
                5.4266824419412347550E+12, 4.3386146580707264428E+14/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
!        5.5  <  |X|  <=  8.0
!--------------------------------------------------------------------
      DATA PY2/ 2.1363534169313901632E+04,-1.0085539923498211426E+07, &
                2.1958827170518100757E+09,-1.9363051266772083678E+11, &
               -1.2829912364088687306E+11, 6.7016641869173237784E+14, &
               -8.0728726905150210443E+15,-1.7439661319197499338E+01/
      DATA QY2/ 8.7903362168128450017E+02, 5.3924739209768057030E+05, &
                2.4727219475672302327E+08, 8.6926121104209825246E+10, &
                2.2598377924042897629E+13, 3.9272425569640309819E+15, &
                3.4563724628846457519E+17/
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
      DATA P0/3.4806486443249270347E+03, 2.1170523380864944322E+04, &
              4.1345386639580765797E+04, 2.2779090197304684302E+04, &
              8.8961548424210455236E-01, 1.5376201909008354296E+02/
      DATA Q0/3.5028735138235608207E+03, 2.1215350561880115730E+04, &
              4.1370412495510416640E+04, 2.2779090197304684318E+04, &
              1.5711159858080893649E+02/
      DATA P1/-2.2300261666214198472E+01,-1.1183429920482737611E+02, &
              -1.8591953644342993800E+02,-8.9226600200800094098E+01, &
              -8.8033303048680751817E-03,-1.2441026745835638459E+00/
      DATA Q1/1.4887231232283756582E+03, 7.2642780169211018836E+03, &
              1.1951131543434613647E+04, 5.7105024128512061905E+03, &
              9.0593769594993125859E+01/
!-------------------------------------------------------------------
!  Check for error conditions
!-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. (ARG .LE. ZERO)) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) GO TO 800
      IF (AX .LE. XSMALL) THEN
         IF (JINT .EQ. 0) THEN
               RESULT = ONE
            ELSE
               RESULT = PI2 * (LOG(AX) + CONS)
         END IF
         GO TO 2000
      END IF
!-------------------------------------------------------------------
!  Calculate J0 for appropriate interval, preserving
!     accuracy near the zero of J0
!-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(5) * ZSQ + PJ0(6)) * ZSQ + PJ0(7)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            WSQ = ONE - ZSQ / SIXTY4
            XNUM = PJ1(7) * WSQ + PJ1(8)
            XDEN = WSQ + QJ1(7)
            DO 220 I = 1, 6
               XNUM = XNUM * WSQ + PJ1(I)
               XDEN = XDEN * WSQ + QJ1(I)
  220       CONTINUE
            PROD = (AX + XJ1) * ((AX - XJ11/TWO56) - XJ12)
      END IF
      RESULT = PROD * XNUM / XDEN
      IF (JINT .EQ. 0) GO TO 2000
!-------------------------------------------------------------------
!  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!    where xn is a zero of Y0
!-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE IF (AX .LE. FIVE5) THEN
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
         ELSE
            UP = (AX-XY21/TWO56)-XY22
            XY = XY2
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
!-------------------------------------------------------------------
!  Now calculate Y0 for appropriate interval, preserving
!     accuracy near the zero of Y0
!-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            XNUM = PY0(6) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 5
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE IF (AX .LE. FIVE5) THEN
            XNUM = PY1(7) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 6
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
         ELSE
            XNUM = PY2(8) * ZSQ + PY2(1)
            XDEN = ZSQ + QY2(1)
            DO 380 I = 2, 7
               XNUM = XNUM * ZSQ + PY2(I)
               XDEN = XDEN * ZSQ + QY2(I)
  380       CONTINUE
      END IF
      RESULT = RESJ + UP * DOWN * XNUM / XDEN
      GO TO 2000
!-------------------------------------------------------------------
!  Calculate J0 or Y0 for |ARG|  >  8.0
!-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AX / TWOPI
      W = AINT(W) + ONEOV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(5) * ZSQ + P0(6)
      XDEN = ZSQ + Q0(5)
      UP = P1(5) * ZSQ + P1(6)
      DOWN = ZSQ + Q1(5)
      DO 850 I = 1, 4
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = SQRT(PI2/AX) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = SQRT(PI2/AX) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
 2000 RETURN
      END
!--------------------------------------------------------------------
!   Subroutine CALJY1 really doing work
!
!   This routine computes first-order Bessel functions of the first and
!   second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!   for Y1, and |X| <= XMAX for J1.
!
!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the 
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESSJ1(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1. 
!   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).
!
!     Approximate values for some important machines are
!
!                          eps      XMAX     XSMALL      XINF  
!
!  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
!
!*******************************************************************
!
! Error Returns
!
!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY1 is called with a negative or zero argument.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, COS, LOG, SIN, SQRT
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division 
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: November 10, 1987
!
!-------------------------------------------------------------------
    SUBROUTINE CALJY1(ARG,RESULT,JINT)
    INTEGER I,JINT
    REAL ::  ARG,AX,DOWN,EIGHT,FOUR,HALF,PI2,PJ0,PJ1,PLG,PROD,PY0, &
         PY1,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,Q0,Q1,RESJ,RESULT, &
         RTPI2,R0,R1,THROV8,TWOPI,TWOPI1,TWOPI2,TWO56,UP,W,WSQ, &
         XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,XJ1,XJ01,XJ02,XJ11,XJ12, &
         XY,XY0,XY01,XY02,XY1,XY11,XY12,Z,ZERO,ZSQ 
    DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(7),PY1(9),P0(6),P1(6), &
         QJ0(5),QJ1(7),QLG(4),QY0(6),QY1(8),Q0(6),Q1(6)
!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
    DATA EIGHT/8.0E0/, &
         FOUR/4.0E0/,HALF/0.5E0/,THROV8/0.375E0/, &
         PI2/6.3661977236758134308E-1/,P17/1.716E-1/ &
         TWOPI/6.2831853071795864769E+0/,ZERO/0.0E0/, &
         TWOPI1/6.28125E0/,TWOPI2/1.9353071795864769253E-03/ &
         TWO56/256.0E+0/,RTPI2/7.9788456080286535588E-1/
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
    DATA XMAX/8.19E+03/,XSMALL/1.22E-09/,XINF/1.7E+38/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
    DATA XJ0/3.8317059702075123156E+0/,XJ1/7.0155866698156187535E+0/, &
          XY0/2.1971413260310170351E+0/,XY1/5.4296810407941351328E+0/, &
          XJ01/ 981.0E+0/, XJ02/-3.2527979248768438556E-04/, &
          XJ11/1796.0E+0/, XJ12/-3.8330184381246462950E-05/, &
          XY01/ 562.0E+0/, XY02/ 1.8288260310170351490E-03/, &
          XY11/1390.0E+0/, XY12/-6.4592058648672279948E-06/
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02, &
         -5.4989956895857911039E+02,3.5687548468071500413E+02/
    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02, &
         -3.3442903192607538956E+02,1.7843774234035750207E+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
    DATA PJ0/9.8062904098958257677E+05,-1.1548696764841276794E+08, &
         6.6781041261492395835E+09,-1.4258509801366645672E+11, &
         -4.4615792982775076130E+03, 1.0650724020080236441E+01, &
         -1.0767857011487300348E-02/
    DATA QJ0/5.9117614494174794095E+05, 2.0228375140097033958E+08, &
         4.2091902282580133541E+10, 4.1868604460820175290E+12, &
         1.0742272239517380498E+03/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
    DATA PJ1/4.6179191852758252280E+00,-7.1329006872560947377E+03, &
         4.5039658105749078904E+06,-1.4437717718363239107E+09, &
         2.3569285397217157313E+11,-1.6324168293282543629E+13, &
         1.1357022719979468624E+14, 1.0051899717115285432E+15/
    DATA QJ1/1.1267125065029138050E+06, 6.4872502899596389593E+08, &
         2.7622777286244082666E+11, 8.4899346165481429307E+13, &
         1.7128800897135812012E+16, 1.7253905888447681194E+18, &
         1.3886978985861357615E+03/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
    DATA PY0/2.2157953222280260820E+05,-5.9157479997408395984E+07, &
         7.2144548214502560419E+09,-3.7595974497819597599E+11, &
         5.4708611716525426053E+12, 4.0535726612579544093E+13, &
         -3.1714424660046133456E+02/
    DATA QY0/8.2079908168393867438E+02, 3.8136470753052572164E+05, &
         1.2250435122182963220E+08, 2.7800352738690585613E+10, &
         4.1272286200406461981E+12, 3.0737873921079286084E+14/
!--------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
!        4.0  <  |X|  <=  8.0
!--------------------------------------------------------------------
    DATA PY1/ 1.9153806858264202986E+06,-1.1957961912070617006E+09, &
         3.7453673962438488783E+11,-5.9530713129741981618E+13, &
         4.0686275289804744814E+15,-2.3638408497043134724E+16, &
         -5.6808094574724204577E+18, 1.1514276357909013326E+19, &
         -1.2337180442012953128E+03/
    DATA QY1/ 1.2855164849321609336E+03, 1.0453748201934079734E+06, &
         6.3550318087088919566E+08, 3.0221766852960403645E+11, &
         1.1187010065856971027E+14, 3.0837179548112881950E+16, &
         5.6968198822857178911E+18, 5.3321844313316185697E+20/
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
    DATA P0/-1.0982405543459346727E+05,-1.5235293511811373833E+06, &
         -6.6033732483649391093E+06,-9.9422465050776411957E+06, &
         -4.4357578167941278571E+06,-1.6116166443246101165E+03/
    DATA Q0/-1.0726385991103820119E+05,-1.5118095066341608816E+06, &
         -6.5853394797230870728E+06,-9.9341243899345856590E+06, &
         -4.4357578167941278568E+06,-1.4550094401904961825E+03/
    DATA P1/ 1.7063754290207680021E+03, 1.8494262873223866797E+04, &
         6.6178836581270835179E+04, 8.5145160675335701966E+04, &
         3.3220913409857223519E+04, 3.5265133846636032186E+01/
    DATA Q1/ 3.7890229745772202641E+04, 4.0029443582266975117E+05, &
         1.4194606696037208929E+06, 1.8194580422439972989E+06, &
         7.0871281941028743574E+05, 8.6383677696049909675E+02/
!-------------------------------------------------------------------
!  Check for error conditions
!-------------------------------------------------------------------
    AX = ABS(ARG)
    IF ((JINT .EQ. 1) .AND. ((ARG .LE. ZERO) .OR. &
         ((ARG .LT. HALF) .AND. (AX*XINF .LT. PI2)))) THEN
       RESULT = -XINF
       GO TO 2000
    ELSE IF (AX .GT. XMAX) THEN
       RESULT = ZERO
       GO TO 2000
    END IF
    IF (AX .GT. EIGHT) THEN
       GO TO 800
    ELSE IF (AX .LE. XSMALL) THEN
       IF (JINT .EQ. 0) THEN
          RESULT = ARG * HALF
       ELSE
          RESULT = -PI2 / AX
       END IF
       GO TO 2000
    END IF
!-------------------------------------------------------------------
!  Calculate J1 for appropriate interval, preserving
!     accuracy near the zero of J1
!-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(7) * ZSQ + PJ0(6)) * ZSQ + PJ0(5)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ARG * ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            XNUM = PJ1(1)
            XDEN = (ZSQ + QJ1(7)) * ZSQ + QJ1(1)
            DO 220 I = 2, 6
               XNUM = XNUM * ZSQ + PJ1(I)
               XDEN = XDEN * ZSQ + QJ1(I)
  220       CONTINUE
            XNUM = XNUM * (AX - EIGHT) * (AX + EIGHT) + PJ1(7)
            XNUM = XNUM * (AX - FOUR) * (AX + FOUR) + PJ1(8)
            PROD = ARG * ((AX - XJ11/TWO56) - XJ12) * (AX + XJ1)
      END IF
      RESULT = PROD * (XNUM / XDEN)
      IF (JINT .EQ. 0) GO TO 2000
!-------------------------------------------------------------------
!  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
!    where xn is a zero of Y1
!-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
!-------------------------------------------------------------------
!  Now calculate Y1 for appropriate interval, preserving
!     accuracy near the zero of Y1
!-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            XNUM = PY0(7) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 6
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE
            XNUM = PY1(9) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 8
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
      END IF
      RESULT = RESJ + (UP*DOWN/AX) * XNUM / XDEN
      GO TO 2000
!-------------------------------------------------------------------
!  Calculate J1 or Y1 for |ARG|  >  8.0
!-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AINT(AX/TWOPI) + THROV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(6)
      XDEN = ZSQ + Q0(6)
      UP = P1(6)
      DOWN = ZSQ + Q1(6)
      DO 850 I = 1, 5
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = (RTPI2/SQRT(AX)) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = (RTPI2/SQRT(AX)) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
      IF ((JINT .EQ. 0) .AND. (ARG .LT. ZERO)) RESULT = -RESULT
 2000 RETURN
      END
!=======================================================================================
!                         DOUBLE PRECISION VERSION
!=======================================================================================
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for Bessel functions
!   of the first kind of order zero for arguments  |X| <= XMAX
!   (see comments heading CALJY0).
!
!--------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DBESSJ0(X)
      INTEGER JINT
      DOUBLE PRECISION  X, RESULT
!--------------------------------------------------------------------
      JINT=0
      CALL DCALJY0(X,RESULT,JINT)
      DBESSJ0 = RESULT
      RETURN
      END
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for Bessel functions
!   of the second kind of order zero for arguments 0 < X <= XMAX
!   (see comments heading CALJY0).
!
!--------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DBESSY0(X)
      INTEGER JINT
      DOUBLE PRECISION  X, RESULT
      JINT=1
      CALL DCALJY0(X,RESULT,JINT)
      DBESSY0 = RESULT
      RETURN
      END
!--------------------------------------------------------------------
!
!   This subprogram computes approximate values for Bessel functions
!   of the first kind of order zero for arguments  |X| <= XMAX
!   (see comments heading CALJY1).
!
!--------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DBESSJ1(X)
      INTEGER JINT
      DOUBLE PRECISION RESULT,X
      JINT=0
      CALL DCALJY1(X,RESULT,JINT)
      DBESSJ1 = RESULT
      RETURN
      END
!--------------------------------------------------------------------
!
!   This subprogram computes approximate values for Bessel functions
!   of the second kind of order zero for arguments 0 < X <= XMAX
!   (see comments heading CALJY1).
!
!--------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DBESSY1(X)
      INTEGER JINT
      DOUBLE PRECISION RESULT,X
      JINT=1
      CALL DCALJY1(X,RESULT,JINT)
      DBESSY1 = RESULT
      RETURN
      END
!-----------------------------------------------------------------------
!---------------------------------------------------------------------
!
!   This packet computes zero-order Bessel functions of the first and
!   second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!   for Y0, and |X| <= XMAX for J0.  It contains two function-type
!   subprograms,  DBESSJ0  and  DBESSY0,  and one subroutine-type
!   subprogram,  DCALJY0.  The calling statements for the primary
!   entries are:
!
!           Y = DBESSJ0(X)
!   and
!           Y = DBESSY0(X),
!
!   where the entry points correspond to the functions J0(X) and Y0(X),
!   respectively.  The routine  CALJY0  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY0  with
!   the statement
!           CALL DCALJY0(ARG,RESULT,JINT),
!   where the parameter usage is as follows:
!
!      Function                  Parameters for CALJY0
!       call              ARG             RESULT          JINT
!
!     DBESSJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
!     DBESSY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1
!
!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the 
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ0(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!*******************************************************************
!
!   Explanation of machine-dependent constants
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1. 
!   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).
!
!     Approximate values for some important machines are
!
!                          eps      XMAX     XSMALL      XINF  
!
!  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
!
!*******************************************************************
!*******************************************************************
!
!  Error Returns
!
!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY0 is called with a negative or zero argument.
!
!
!  Intrinsic functions required are:
!
!     ABS, AINT, COS, LOG, SIN, SQRT
!
!
!  Latest modification: June 2, 1989
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division 
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!--------------------------------------------------------------------
      SUBROUTINE DCALJY0(ARG,RESULT,JINT)
      INTEGER I,JINT
      DOUBLE PRECISION &
            ARG,AX,CONS,DOWN,EIGHT,FIVE5,FOUR,ONE,ONEOV8,PI2,PJ0, &
            PJ1,PLG,PROD,PY0,PY1,PY2,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1, &
            QY2,Q0,Q1,RESJ,RESULT,R0,R1,SIXTY4,THREE,TWOPI,TWOPI1, &
            TWOPI2,TWO56,UP,W,WSQ,XDEN,XINF,XMAX,XNUM,XSMALL,XJ0, &
            XJ1,XJ01,XJ02,XJ11,XJ12,XY,XY0,XY01,XY02,XY1,XY11,XY12, &
            XY2,XY21,XY22,Z,ZERO,ZSQ
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(6),PY1(7),PY2(8),P0(6),P1(6), &
                QJ0(5),QJ1(7),QLG(4),QY0(5),QY1(6),QY2(7),Q0(5),Q1(5)
!-------------------------------------------------------------------
!  Mathematical constants
!    CONS = ln(.5) + Euler's gamma
!-------------------------------------------------------------------
      DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0D0,1.0D0,3.0D0,4.0D0,8.0D0/, &
           FIVE5,SIXTY4,ONEOV8,P17/5.5D0,64.0D0,0.125D0,1.716D-1/, &
           TWO56,CONS/256.0D0,-1.1593151565841244881D-1/, &
           PI2,TWOPI/6.3661977236758134308D-1,6.2831853071795864769D0/, &
           TWOPI1,TWOPI2/6.28125D0,1.9353071795864769253D-3/
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
      DATA XJ0/2.4048255576957727686D+0/,XJ1/5.5200781102863106496D+0/, &
          XY0/8.9357696627916752158D-1/,XY1/3.9576784193148578684D+0/, &
          XY2/7.0860510603017726976D+0/, &
          XJ01/ 616.0D+0/, XJ02/-1.4244423042272313784D-03/, &
          XJ11/1413.0D+0/, XJ12/ 5.4686028631064959660D-04/, &
          XY01/ 228.0D+0/, XY02/ 2.9519662791675215849D-03/, &
          XY11/1013.0D+0/, XY12/ 6.4716931485786837568D-04/, &
          XY21/1814.0D+0/, XY22/ 1.1356030177269762362D-04/
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02, &
              -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02, &
              -3.3442903192607538956D+02,1.7843774234035750207D+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
      DATA PJ0/6.6302997904833794242D+06,-6.2140700423540120665D+08, &
              2.7282507878605942706D+10,-4.1298668500990866786D+11, &
             -1.2117036164593528341D-01, 1.0344222815443188943D+02, &
             -3.6629814655107086448D+04/
      DATA QJ0/4.5612696224219938200D+05, 1.3985097372263433271D+08, &
              2.6328198300859648632D+10, 2.3883787996332290397D+12, &
              9.3614022392337710626D+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
      DATA PJ1/4.4176707025325087628D+03, 1.1725046279757103576D+04, &
              1.0341910641583726701D+04,-7.2879702464464618998D+03, &
             -1.2254078161378989535D+04,-1.8319397969392084011D+03, &
              4.8591703355916499363D+01, 7.4321196680624245801D+02/
      DATA QJ1/3.3307310774649071172D+02,-2.9458766545509337327D+03, &
              1.8680990008359188352D+04,-8.4055062591169562211D+04, &
              2.4599102262586308984D+05,-3.5783478026152301072D+05, &
             -2.5258076240801555057D+01/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  3.0
!--------------------------------------------------------------------
      DATA PY0/1.0102532948020907590D+04,-2.1287548474401797963D+06, &
              2.0422274357376619816D+08,-8.3716255451260504098D+09, &
              1.0723538782003176831D+11,-1.8402381979244993524D+01/
      DATA QY0/6.6475986689240190091D+02, 2.3889393209447253406D+05, &
              5.5662956624278251596D+07, 8.1617187777290363573D+09, &
              5.8873865738997033405D+11/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
!        3.0  <  |X|  <=  5.5
!--------------------------------------------------------------------
      DATA PY1/-1.4566865832663635920D+04, 4.6905288611678631510D+06, &
              -6.9590439394619619534D+08, 4.3600098638603061642D+10, &
              -5.5107435206722644429D+11,-2.2213976967566192242D+13, &
               1.7427031242901594547D+01/
      DATA QY1/ 8.3030857612070288823D+02, 4.0669982352539552018D+05, &
               1.3960202770986831075D+08, 3.4015103849971240096D+10, &
               5.4266824419412347550D+12, 4.3386146580707264428D+14/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
!        5.5  <  |X|  <=  8.0
!--------------------------------------------------------------------
      DATA PY2/ 2.1363534169313901632D+04,-1.0085539923498211426D+07, &
               2.1958827170518100757D+09,-1.9363051266772083678D+11, &
              -1.2829912364088687306D+11, 6.7016641869173237784D+14, &
              -8.0728726905150210443D+15,-1.7439661319197499338D+01/
      DATA QY2/ 8.7903362168128450017D+02, 5.3924739209768057030D+05, &
               2.4727219475672302327D+08, 8.6926121104209825246D+10, &
               2.2598377924042897629D+13, 3.9272425569640309819D+15, &
               3.4563724628846457519D+17/
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
      DATA P0/3.4806486443249270347D+03, 2.1170523380864944322D+04, &
             4.1345386639580765797D+04, 2.2779090197304684302D+04, &
             8.8961548424210455236D-01, 1.5376201909008354296D+02/
      DATA Q0/3.5028735138235608207D+03, 2.1215350561880115730D+04, &
             4.1370412495510416640D+04, 2.2779090197304684318D+04, &
             1.5711159858080893649D+02/
      DATA P1/-2.2300261666214198472D+01,-1.1183429920482737611D+02, &
             -1.8591953644342993800D+02,-8.9226600200800094098D+01, &
             -8.8033303048680751817D-03,-1.2441026745835638459D+00/
      DATA Q1/1.4887231232283756582D+03, 7.2642780169211018836D+03, &
             1.1951131543434613647D+04, 5.7105024128512061905D+03, &
             9.0593769594993125859D+01/
!-------------------------------------------------------------------
!  Check for error conditions
!-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. (ARG .LE. ZERO)) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) GO TO 800
      IF (AX .LE. XSMALL) THEN
         IF (JINT .EQ. 0) THEN
               RESULT = ONE
            ELSE
               RESULT = PI2 * (LOG(AX) + CONS)
         END IF
         GO TO 2000
      END IF
!-------------------------------------------------------------------
!  Calculate J0 for appropriate interval, preserving
!     accuracy near the zero of J0
!-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(5) * ZSQ + PJ0(6)) * ZSQ + PJ0(7)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            WSQ = ONE - ZSQ / SIXTY4
            XNUM = PJ1(7) * WSQ + PJ1(8)
            XDEN = WSQ + QJ1(7)
            DO 220 I = 1, 6
               XNUM = XNUM * WSQ + PJ1(I)
               XDEN = XDEN * WSQ + QJ1(I)
  220       CONTINUE
            PROD = (AX + XJ1) * ((AX - XJ11/TWO56) - XJ12)
      END IF
      RESULT = PROD * XNUM / XDEN
      IF (JINT .EQ. 0) GO TO 2000
!-------------------------------------------------------------------
!  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!    where xn is a zero of Y0
!-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE IF (AX .LE. FIVE5) THEN
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
         ELSE
            UP = (AX-XY21/TWO56)-XY22
            XY = XY2
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
!-------------------------------------------------------------------
!  Now calculate Y0 for appropriate interval, preserving
!     accuracy near the zero of Y0
!-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            XNUM = PY0(6) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 5
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE IF (AX .LE. FIVE5) THEN
            XNUM = PY1(7) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 6
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
         ELSE
            XNUM = PY2(8) * ZSQ + PY2(1)
            XDEN = ZSQ + QY2(1)
            DO 380 I = 2, 7
               XNUM = XNUM * ZSQ + PY2(I)
               XDEN = XDEN * ZSQ + QY2(I)
  380       CONTINUE
      END IF
      RESULT = RESJ + UP * DOWN * XNUM / XDEN
      GO TO 2000
!-------------------------------------------------------------------
!  Calculate J0 or Y0 for |ARG|  >  8.0
!-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AX / TWOPI
      W = AINT(W) + ONEOV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(5) * ZSQ + P0(6)
      XDEN = ZSQ + Q0(5)
      UP = P1(5) * ZSQ + P1(6)
      DOWN = ZSQ + Q1(5)
      DO 850 I = 1, 4
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = SQRT(PI2/AX) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = SQRT(PI2/AX) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
 2000 RETURN
!---------- Last line of DCALJY0 ----------
      END
      SUBROUTINE DCALJY1(ARG,RESULT,JINT)
!---------------------------------------------------------------------
!
!   This packet computes first-order Bessel functions of the first and
!   second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!   for Y1, and |X| <= XMAX for J1.  It contains two function-type
!   subprograms,  DBESSJ1  and  DBESSY1,  and one subroutine-type
!   subprogram,  DCALJY1.  The calling statements for the primary
!   entries are:
!
!           Y = DBESSJ1(X)
!   and
!           Y = DBESSY1(X),
!
!   where the entry points correspond to the functions J1(X) and Y1(X),
!   respectively.  The routine  CALJY1  is intended for internal packet
!   use only, all computations within the packet being concentrated in
!   this one routine.  The function subprograms invoke  CALJY1  with
!   the statement
!           CALL DCALJY1(ARG,RESULT,JINT),
!   where the parameter usage is as follows:
!
!      Function                  Parameters for CALJY1
!       call              ARG             RESULT          JINT
!
!     DBESSJ1(ARG)     |ARG| .LE. XMAX       J1(ARG)          0
!     DBESSY1(ARG)   0 .LT. ARG .LE. XMAX    Y1(ARG)          1
!
!   The main computation uses unpublished minimax rational
!   approximations for X .LE. 8.0, and an approximation from the 
!   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
!   New York, 1968, for arguments larger than 8.0   Part of this
!   transportable packet is patterned after the machine-dependent
!   FUNPACK program BESJ1(X), but cannot match that version for
!   efficiency or accuracy.  This version uses rational functions
!   that are theoretically accurate to at least 18 significant decimal
!   digits for X <= 8, and at least 18 decimal places for X > 8.  The
!   accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-
!   dependent constants.
!
!*******************************************************************
!
!   Explanation of machine-dependent constants
!
!   XINF   = largest positive machine number
!   XMAX   = largest acceptable argument.  The functions AINT, SIN
!            and COS must perform properly for  ABS(X) .LE. XMAX.
!            We recommend that XMAX be a small integer multiple of
!            sqrt(1/eps), where eps is the smallest positive number
!            such that  1+eps > 1. 
!   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
!            to machine precision for all  ABS(X) .LE. XSMALL.
!            We recommend that  XSMALL < sqrt(eps)/beta, where beta
!            is the floating-point radix (usually 2 or 16).
!
!     Approximate values for some important machines are
!
!                          eps      XMAX     XSMALL      XINF  
!
!    C 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
!  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
!  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
!  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
!  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
!  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
!  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
!
!*******************************************************************
!*******************************************************************
!
!  Error Returns
!
!  The program returns the value zero for  X .GT. XMAX, and returns
!    -XINF when BESLY1 is called with a negative or zero argument.
!
!
!  Intrinsic functions required are:
!
!     ABS, AINT, COS, LOG, SIN, SQRT
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division 
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: November 10, 1987
!
!--------------------------------------------------------------------
      INTEGER I,JINT
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(7),PY1(9),P0(6),P1(6), &
               QJ0(5),QJ1(7),QLG(4),QY0(6),QY1(8),Q0(6),Q1(6)
      DOUBLE PRECISION &
        ARG,AX,DOWN,EIGHT,FOUR,HALF,PI2,PJ0,PJ1,PLG,PROD,PY0, &
        PY1,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,Q0,Q1,RESJ,RESULT, &
        RTPI2,R0,R1,THROV8,TWOPI,TWOPI1,TWOPI2,TWO56,UP,W,WSQ, &
        XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,XJ1,XJ01,XJ02,XJ11,XJ12, &
        XY,XY0,XY01,XY02,XY1,XY11,XY12,Z,ZERO,ZSQ
!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
      DATA EIGHT/8.0D0/, &
          FOUR/4.0D0/,HALF/0.5D0/,THROV8/0.375D0/, &
          PI2/6.3661977236758134308D-1/,P17/1.716D-1/ &
          TWOPI/6.2831853071795864769D+0/,ZERO/0.0D0/, &
          TWOPI1/6.28125D0/,TWOPI2/1.9353071795864769253D-03/ &
          TWO56/256.0D+0/,RTPI2/7.9788456080286535588D-1/
!-------------------------------------------------------------------
!  Machine-dependent constants
!-------------------------------------------------------------------
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
!-------------------------------------------------------------------
!  Zeroes of Bessel functions
!-------------------------------------------------------------------
      DATA XJ0/3.8317059702075123156D+0/,XJ1/7.0155866698156187535D+0/, &
          XY0/2.1971413260310170351D+0/,XY1/5.4296810407941351328D+0/, &
          XJ01/ 981.0D+0/, XJ02/-3.2527979248768438556D-04/, &
          XJ11/1796.0D+0/, XJ12/-3.8330184381246462950D-05/, &
          XY01/ 562.0D+0/, XY02/ 1.8288260310170351490D-03/, &
          XY11/1390.0D+0/, XY12/-6.4592058648672279948D-06/
!-------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a)
!--------------------------------------------------------------------
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02, &
              -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02, &
              -3.3442903192607538956D+02,1.7843774234035750207D+02/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
      DATA PJ0/9.8062904098958257677D+05,-1.1548696764841276794D+08, &
            6.6781041261492395835D+09,-1.4258509801366645672D+11, &
           -4.4615792982775076130D+03, 1.0650724020080236441D+01, &
           -1.0767857011487300348D-02/
      DATA QJ0/5.9117614494174794095D+05, 2.0228375140097033958D+08, &
            4.2091902282580133541D+10, 4.1868604460820175290D+12, &
            1.0742272239517380498D+03/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
!-------------------------------------------------------------------
      DATA PJ1/4.6179191852758252280D+00,-7.1329006872560947377D+03, &
            4.5039658105749078904D+06,-1.4437717718363239107D+09, &
            2.3569285397217157313D+11,-1.6324168293282543629D+13, &
            1.1357022719979468624D+14, 1.0051899717115285432D+15/
      DATA QJ1/1.1267125065029138050D+06, 6.4872502899596389593D+08, &
            2.7622777286244082666D+11, 8.4899346165481429307D+13, &
            1.7128800897135812012D+16, 1.7253905888447681194D+18, &
            1.3886978985861357615D+03/
!-------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
!        XSMALL  <  |X|  <=  4.0
!--------------------------------------------------------------------
      DATA PY0/2.2157953222280260820D+05,-5.9157479997408395984D+07, &
              7.2144548214502560419D+09,-3.7595974497819597599D+11, &
              5.4708611716525426053D+12, 4.0535726612579544093D+13, &
             -3.1714424660046133456D+02/
      DATA QY0/8.2079908168393867438D+02, 3.8136470753052572164D+05, &
              1.2250435122182963220D+08, 2.7800352738690585613D+10, &
              4.1272286200406461981D+12, 3.0737873921079286084D+14/
!--------------------------------------------------------------------
!  Coefficients for rational approximation of
!    (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
!        4.0  <  |X|  <=  8.0
!--------------------------------------------------------------------
      DATA PY1/ 1.9153806858264202986D+06,-1.1957961912070617006D+09, &
               3.7453673962438488783D+11,-5.9530713129741981618D+13, &
               4.0686275289804744814D+15,-2.3638408497043134724D+16, &
              -5.6808094574724204577D+18, 1.1514276357909013326D+19, &
              -1.2337180442012953128D+03/
      DATA QY1/ 1.2855164849321609336D+03, 1.0453748201934079734D+06, &
               6.3550318087088919566D+08, 3.0221766852960403645D+11, &
               1.1187010065856971027D+14, 3.0837179548112881950D+16, &
               5.6968198822857178911D+18, 5.3321844313316185697D+20/
!-------------------------------------------------------------------
!  Coefficients for Hart,s approximation,  |X| > 8.0
!-------------------------------------------------------------------
      DATA P0/-1.0982405543459346727D+05,-1.5235293511811373833D+06, &
              -6.6033732483649391093D+06,-9.9422465050776411957D+06, &
              -4.4357578167941278571D+06,-1.6116166443246101165D+03/
      DATA Q0/-1.0726385991103820119D+05,-1.5118095066341608816D+06, &
              -6.5853394797230870728D+06,-9.9341243899345856590D+06, &
              -4.4357578167941278568D+06,-1.4550094401904961825D+03/
      DATA P1/ 1.7063754290207680021D+03, 1.8494262873223866797D+04, &
               6.6178836581270835179D+04, 8.5145160675335701966D+04, &
               3.3220913409857223519D+04, 3.5265133846636032186D+01/
      DATA Q1/ 3.7890229745772202641D+04, 4.0029443582266975117D+05, &
               1.4194606696037208929D+06, 1.8194580422439972989D+06, &
               7.0871281941028743574D+05, 8.6383677696049909675D+02/
!-------------------------------------------------------------------
!  Check for error conditions
!-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. ((ARG .LE. ZERO) .OR. &
        ((ARG .LT. HALF) .AND. (AX*XINF .LT. PI2)))) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) THEN
            GO TO 800
         ELSE IF (AX .LE. XSMALL) THEN
            IF (JINT .EQ. 0) THEN
                  RESULT = ARG * HALF
               ELSE
                  RESULT = -PI2 / AX
            END IF
            GO TO 2000
      END IF
!-------------------------------------------------------------------
!  Calculate J1 for appropriate interval, preserving
!     accuracy near the zero of J1
!-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(7) * ZSQ + PJ0(6)) * ZSQ + PJ0(5)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ARG * ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            XNUM = PJ1(1)
            XDEN = (ZSQ + QJ1(7)) * ZSQ + QJ1(1)
            DO 220 I = 2, 6
               XNUM = XNUM * ZSQ + PJ1(I)
               XDEN = XDEN * ZSQ + QJ1(I)
  220       CONTINUE
            XNUM = XNUM * (AX - EIGHT) * (AX + EIGHT) + PJ1(7)
            XNUM = XNUM * (AX - FOUR) * (AX + FOUR) + PJ1(8)
            PROD = ARG * ((AX - XJ11/TWO56) - XJ12) * (AX + XJ1)
      END IF
      RESULT = PROD * (XNUM / XDEN)
      IF (JINT .EQ. 0) GO TO 2000
!-------------------------------------------------------------------
!  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
!    where xn is a zero of Y1
!-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
!-------------------------------------------------------------------
!  Now calculate Y1 for appropriate interval, preserving
!     accuracy near the zero of Y1
!-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            XNUM = PY0(7) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 6
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE
            XNUM = PY1(9) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 8
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
      END IF
      RESULT = RESJ + (UP*DOWN/AX) * XNUM / XDEN
      GO TO 2000
!-------------------------------------------------------------------
!  Calculate J1 or Y1 for |ARG|  >  8.0
!-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AINT(AX/TWOPI) + THROV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(6)
      XDEN = ZSQ + Q0(6)
      UP = P1(6)
      DOWN = ZSQ + Q1(6)
      DO 850 I = 1, 5
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = (RTPI2/SQRT(AX)) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = (RTPI2/SQRT(AX)) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
      IF ((JINT .EQ. 0) .AND. (ARG .LT. ZERO)) RESULT = -RESULT
 2000 RETURN
!---------- Last card of DCALJY1 ----------
      END      
end module besselFunctions
