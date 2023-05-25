!**********************************************************************
! Copyright 2016,2017                                                 *
! Christine Groot Zwaaftink                                           *
!                                                                     *
! This file is part of FLEXDUST.                                      *
!                                                                     *
! FLEXDUST is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXDUST is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXDUST.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************
subroutine getSizeDistribution(inisize, inifrac, typeSizeDistr, numberSpecies)
implicit none

!***********************************************************************
integer :: numberSpecies, typeSizeDistr
real(8), dimension(0:numberSpecies-1):: inisize
real(8), dimension(0:numberSpecies-1, 0:3):: inifrac
!***********************************************************************

if(typeSizeDistr.eq.1) then
!Define particle sizes of 15 species (micrometers, not used in model but given for completeness and might be used if a function to write species files is added
!***********************************************************************
inisize(0) = 2.993795e-01
inisize(1) = 1.051838e+00
inisize(2) = 2.029619e+00
inisize(3) = 3.123928e+00
inisize(4) = 4.240266e+00
inisize(5) = 5.190362e+00
inisize(6) = 6.577171e+00
inisize(7) = 8.237095e+00
inisize(8) = 1.046172e+01
inisize(9) = 1.322873e+01
inisize(10)= 1.996804e+01
inisize(11)= 3.038633e+01
inisize(12)= 3.795915e+01
inisize(13)= 4.542401e+01
inisize(14)= 5.637308e+01
!************************************************************************

!Define mass fraction per particle size
!***********************************************************************
inifrac(0,0) = 8.796671e-06
inifrac(1,0) = 9.398633e-05
inifrac(2,0) = 2.109303e-04
inifrac(3,0) = 5.798745e-04
inifrac(4,0) = 6.266683e-04
inifrac(5,0) = 9.421178e-04
inifrac(6,0) = 2.427683e-03
inifrac(7,0) = 2.823664e-03
inifrac(8,0) = 8.120885e-03
inifrac(9,0) = 9.579366e-03
inifrac(10,0) = 8.400262e-02
inifrac(11,0) = 1.206319e-01
inifrac(12,0) = 1.116229e-01
inifrac(13,0) = 2.194719e-01
inifrac(14,0) = 4.388567e-01


inifrac(0,1) = 8.818364e-05
inifrac(1,1) = 5.542972e-04
inifrac(2,1) = 9.297076e-04
inifrac(3,1) = 2.063497e-03
inifrac(4,1) = 1.927442e-03
inifrac(5,1) = 2.620314e-03
inifrac(6,1) = 5.986410e-03
inifrac(7,1) = 6.235843e-03
inifrac(8,1) = 1.587306e-02
inifrac(9,1) = 1.669190e-02
inifrac(10,1) = 1.173472e-01
inifrac(11,1) = 1.383223e-01
inifrac(12,1) = 1.148907e-01
inifrac(13,1) = 2.063497e-01
inifrac(14,1) = 3.701194e-01


inifrac(0,2) = 7.949126e-03
inifrac(1,2) = 1.589825e-02
inifrac(2,2) = 1.430843e-02
inifrac(3,2) = 2.066773e-02
inifrac(4,2) = 1.430843e-02
inifrac(5,2) = 1.589825e-02
inifrac(6,2) = 2.861685e-02
inifrac(7,2) = 2.384738e-02
inifrac(8,2) = 4.769475e-02
inifrac(9,2) = 3.974563e-02
inifrac(10,2) = 1.828299e-01
inifrac(11,2) = 1.430843e-01
inifrac(12,2) = 9.538951e-02
inifrac(13,2) = 1.430843e-01
inifrac(14,2) = 2.066773e-01

!************************************************************************

elseif(typeSizeDistr.eq.2)then
    !Particle size distribution with equally spaced size bins
    inisize(0) = 0.2
    inisize(1) = 2.2
    inisize(2) = 4.2
    inisize(3) = 6.2
    inisize(4) = 8.2
    inisize(5) = 10.2
    inisize(6) = 12.2
    inisize(7) = 14.2
    inisize(8) = 16.2
    inisize(9) = 18.2

    inifrac(0,0) = 0.44!0.02
    inifrac(1,0) = 0.13!0.05
    inifrac(2,0) = 0.11!0.10
    inifrac(3,0) = 0.10!0.13
    inifrac(4,0) = 0.08!0.16
    inifrac(5,0) = 0.06!0.16
    inifrac(6,0) = 0.04!0.15
    inifrac(7,0) = 0.02!0.11
    inifrac(8,0) = 0.01!0.08
    inifrac(9,0) = 0.01!0.04
    
    inifrac(:,1)=inifrac(:,0)
    inifrac(:,2)=inifrac(:,0)

elseif(typeSizeDistr.eq.3)then
    !Particle size distribution with more smaller size bins and less larger bins
    inisize(0) = 0.04
    inisize(1) = 0.22
    inisize(2) = 0.71
    inisize(3) = 1.30
    inisize(4) = 2.06
    inisize(5) = 3.53
    inisize(6) = 6.10
    inisize(7) = 8.63
    inisize(8) = 12.25
    inisize(9) = 17.32
    
    inifrac(0,0) = 0.0030
    inifrac(1,0) = 0.0163
    inifrac(2,0) = 0.0327
    inifrac(3,0) = 0.0556
    inifrac(4,0) = 0.0820
    inifrac(5,0) = 0.1605
    inifrac(6,0) = 0.2209
    inifrac(7,0) = 0.2425
    inifrac(8,0) = 0.1516
    inifrac(9,0) = 1-sum(inifrac(0:8,0))!0.0342
      
    inifrac(:,1)=inifrac(:,0)
    inifrac(:,2)=inifrac(:,0)
    
else 
    print*, 'Unknown size distribution'
    stop
endif

end subroutine getSizeDistribution