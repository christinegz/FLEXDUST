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
subroutine dustEmmission_MB95_tau(em_mass, soilF, lat,dxdy_degr, inNestNr, ix_wind, iy_wind, ix_wind_n, iy_wind_n,&
                                time_int, scalingFactor, mobilisationThreshold,shearStressThres, inClay, f_clay_tmp, em_flux, &
                                ustarVersion, gridarea)

use par_mod
use com_mod

implicit none

real     :: mfac, em_mass, lat,  dxdy_degr, em_flux, soilF
real*8   :: scalingFactor
real     :: u_star_local, rho_air, eta, f_clay_tmp, f_clay, stress_local
real     :: shearStressThres ! Threshold shear stress
real     :: mobilisationThreshold ! Threshold friction velocity
integer  :: ix_wind, iy_wind, ix_wind_n, iy_wind_n, time_int, inNestNr
real    :: ylatp, ylatm, gridarea, cosfactp, hzone, cosfactm
real    :: mass_tmp_ustar, mass_tmp_stress
logical :: inClay, test, ustarVersion


        if(inNestNr.eq.0)then
            u_star_local=ustar(ix_wind,iy_wind,1,1)
            rho_air=rho(ix_wind,iy_wind,1,1)
            stress_local=surfstr(ix_wind,iy_wind,1,1) 
        else
            u_star_local=ustarn(ix_wind_n,iy_wind_n,1,1,inNestNr)
            rho_air=rhon(ix_wind_n,iy_wind_n,1,1,inNestNr)
            stress_local=surfstrn(ix_wind_n,iy_wind_n,1,1,inNestNr) 
        endif

        !******************************************************************
        !Calculate emitted mass (em) per time step if the threshold for saltation is exceeded
        !******************************************************************	
        test=.false.
        if(ustarVersion .and.u_star_local.ge.mobilisationThreshold)then
            test=.true.
        endif
        
        if(.not.ustarVersion .and. stress_local.ge.shearStressThres)then
             test=.true.
        endif
        
        if(test)then	         
                !get clay fraction
                !*************************************
                if(inClay)then
                    f_clay=real(f_clay_tmp/100)
                    if(f_clay.lt.0)then
                        !No infomation on clay content but there is bare soil, assume 0.05
                        f_clay=0.05
                    endif
                    if(f_clay.gt.0.2)then ! not defined for larger values
                        f_clay=0.15
                    endif
                else
                    f_clay=0.05 !Soil grid point outside map of clay and sand (probably Antarctic)
                    print*, 'Bare soil outside clay/sand map at lat: ', lat , 'assumed f_clay=0.05.'
                endif
                        
                !multiplication factor accounting for time step, area and conversion to kg
                !**********************************************************
                mfac=3600.0*time_int*gridarea	
                                        
                !coefficient eta
                !****************************************
                eta=100.*10**(0.134*f_clay*100-6)
                        
                !Emitted mass in kg for shear stress
                mass_tmp_stress=scalingFactor*eta*rho_air/ga*(sqrt(stress_local/rho_air))**3.0*&
                (1-shearStressThres/stress_local)&
                *(1+sqrt(shearStressThres)/sqrt(stress_local))
                
                !Emitted mass in kg for ustar
                mass_tmp_ustar=scalingFactor*eta*rho_air/ga*u_star_local**3*&
                     (1-mobilisationThreshold**2/u_star_local**2)*(1+mobilisationThreshold/u_star_local)

                if(ustarVersion)then
                    em_mass=em_mass+mass_tmp_ustar*soilF*mfac   !kg
                    em_flux=em_flux+mass_tmp_ustar*soilF*3600.*time_int !kg/m2
                else
                    em_mass=em_mass+mass_tmp_stress*soilF*mfac
                    em_flux=em_flux+mass_tmp_stress*soilF*3600.*time_int 
                endif
              
    	endif
        !******************************************************************	

end subroutine dustEmmission_MB95_tau
