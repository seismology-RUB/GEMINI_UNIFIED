! ==============================================================================
!  Compute complex elastic constants from Q and frequency
! ==============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016,2020 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!--------------------------------------------------------------------
!  Module for computing complex frequency dependent elastic constants 
!  from tabulated density, velocities and Qs in earth model files
!  which are valid for some reference frequency.
!
!  alfv: vertical P-wave velocity
!  alfh: horizontal P-wave velocity
!  betv: vertical S-wave velocity
!  beth: horizontal S-wave velocity
!  eta: anisotropy parameter eta
!  qkinv: inverse qkappa (1/qk)
!  qminv: inverse qmue (1/qm)
!  Output: complex constants A,C,F,L,N, equivalent isotropic kappa and mue
!--------------------------------------------------------------------
 module complexElasticConstants
    use mathConstants
    use string
    implicit none
    integer, parameter :: NELCON = 7
 contains
!---------------------------------------------------------------
!  convert tabulated transversely isotropic constants to
!  complex frequency dependent ones
!  zelcon = elcon*(fratio)**(2*phi/pi)*exp(i*phi)
!  phi = atan(qinv)
!
    subroutine complexElasticVTIConstants(fratio,ro,alfv,alfh,betv,beth,eta,qkinv,qminv,&
                                        & zpa,zpc,zpf,zpl,zpn,zkap,zmue)
    double precision fratio,ro,alfv,alfh,betv,beth,eta,qkinv,qminv,phik,phim
    double precision epa,epc,epn,epl,epf,kappa,mue
    double precision da,dc,df,dl,dn,lamp2mue,lambda
    double complex zpa,zpc,zpn,zpl,zpf,zkap,zmue,zlamp2mue,zlambda,zqk,zqm
!
    phik = datan(qkinv)
    phim = datan(qminv)
    zqk = dcmplx(dcos(phik),dsin(phik))*fratio**(2*phik/mc_pid)
    zqm = dcmplx(dcos(phim),dsin(phim))*fratio**(2*phim/mc_pid)
!
!  solid case
!
    if (betv .gt. 1.e-6) then
        epa=ro*alfh*alfh
        epc=ro*alfv*alfv
        epn=ro*beth*beth
        epl=ro*betv*betv
        epf=eta*(epa-2.*epl)
!
!  isotropic equivalents
!
        kappa=(epc+4.*epa-4.*epn+4.*epf)/9.d0
        mue=(epc+epa+6.*epl+5.*epn-2.*epf)/15.d0
!
!  purely anisotropic parts
!
        lamp2mue=kappa+4.*mue/3.
        lambda=kappa-2.*mue/3.
        da=epa-lamp2mue
        dc=epc-lamp2mue
        df=epf-lambda
        dn=epn-mue
        dl=epl-mue
!
!  complex isotropic constants
!
        zkap=kappa*zqk
        zmue=  mue*zqm
        zlamp2mue=zkap+4.*zmue/3.
        zlambda=zkap-2.*zmue/3.
!
        zpa=zlamp2mue+da
        zpc=zlamp2mue+dc
        zpn=zmue+dn
        zpl=zmue+dl
        zpf=zlambda+df
    else
        kappa=ro*alfv*alfv
        zkap=kappa*zqk
        zmue=0.d0
        zpa=zkap
        zpc=zkap
        zpf=zkap
        zpl=0.d0
        zpn=0.d0
    endif
    end subroutine complexElasticVTIConstants
!---------------------------------------------------------------
!  Calculate complex elastic constants from elastic ones
!
    subroutine complexFromElasticConstants(fratio,attmode,elcon,qinv,zelcon)
    double precision :: fratio
    character (len=*) :: attmode
    double complex, dimension(NELCON) :: elcon,zelcon
    double precision, dimension(:) :: qinv
    double precision :: myqkinv,myqminv,myfratio,phik,phim
    double complex :: epa,epc,epf,epl,epn,kappa,mue
    double complex :: lamp2mue,lambda,da,dc,df,dl,dn,zlamp2mue,zlambda,zqk,zqm,zkap,zmue
!
    epa = elcon(1); epc = elcon(2); epf = elcon(3); epl = elcon(4); epn = elcon(5)
    kappa = elcon(6); mue = elcon(7)
!
    if (attmode.equal.'ELASTIC') then
       zelcon = elcon
       return
    else if (attmode.equal.'ATTENUATION_ONLY') then              ! fratio = 1, qinv /= 0
       myfratio = 1.d0; myqkinv = qinv(1); myqminv = qinv(2)
    else if (attmode.equal.'DISPERSION_ONLY') then               ! qinv = 0
       myfratio = fratio; myqkinv = 0.d0; myqminv = 0.d0
    else if (attmode.equal.'ATTENUATION_AND_DISPERSION') then
       myfratio = fratio; myqkinv = qinv(1); myqminv = qinv(2)
    else
       print *,'complexFromElasticConstants: unknown attenuation mode'
       stop
    endif
!
    phik = datan(myqkinv)
    phim = datan(myqminv)
    zqk = dcmplx(dcos(phik),dsin(phik))*myfratio**(2*phik/mc_pid)
    zqm = dcmplx(dcos(phim),dsin(phim))*myfratio**(2*phim/mc_pid)
!
!  purely anisotropic parts
!
    lamp2mue = kappa+4.*mue/3.
    lambda = kappa-2.*mue/3.
    da = epa-lamp2mue
    dc = epc-lamp2mue
    df = epf-lambda
    dn = epn-mue
    dl = epl-mue
!
!  complex isotropic constants
!
    zkap = kappa*zqk
    zmue =   mue*zqm
    zlamp2mue = zkap+4.*zmue/3.
    zlambda = zkap-2.*zmue/3.
!
    zelcon(1) = zlamp2mue+da
    zelcon(2) = zlamp2mue+dc
    zelcon(3) = zlambda+df
    zelcon(4) = zmue+dn
    zelcon(5) = zmue+dl
    zelcon(6) = zkap
    zelcon(7) = zmue
    end subroutine complexFromElasticConstants
 end module
