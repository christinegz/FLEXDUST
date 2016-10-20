subroutine dustEmmission_MB95_tau(em_mass, soilF, lat,dxdy_degr, inNestNr, ix_wind, iy_wind, ix_wind_n, iy_wind_n,&
                                time_int, scalingFactor, mobilisationThreshold,shearStressThres, inClay, f_clay_tmp, em_flux, &
                                ustarVersion)

use par_mod
use com_mod

implicit none

real     :: mfac, em_mass, soilF, lat,  dxdy_degr, em_flux
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
                        f_clay=0.0 !No infomation on clay content
                    endif
                    if(f_clay.gt.0.2)then ! not defined, use f_clay=0.1
                        f_clay=0.1
                    endif
                else
                    f_clay=0.1 !Should never occur since clay is available at each grid point with soil
                print*, 'Not in Clay'
		endif
        
                !Calculate area copied from FLEXPART
                !****************************************
                ylatp=lat+dxdy_degr
                ylatm=lat
                if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
                    hzone=dxdy_degr*r_earth*pi180
                else
                    cosfactp=cos(ylatp*pi180)
                    cosfactm=cos(ylatm*pi180)
                    if (cosfactp.lt.cosfactm) then
                        hzone=sqrt(1-cosfactp**2)- &
                        sqrt(1-cosfactm**2)
                        !print*, hzone
                        hzone=hzone*r_earth
                    else
                        hzone=sqrt(1-cosfactm**2)- &
                        sqrt(1-cosfactp**2)
                        hzone=hzone*r_earth
                    endif
                endif
                
                gridarea=2.*pi*r_earth*hzone*dxdy_degr/360.
                !****************************************
                
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
                
                     
!                 !Print for debugging
!                if(mass_tmp_ustar/mass_tmp_stress.ne.1.0000)then
!                    print*, 'rho= ',rho_air,', tau= ', stress_local,', tau_t=', shearStressThres,', u*=', u_star_local,&
!                      ', u*t=', mobilisationThreshold, ', mass_tau=', mass_tmp_stress, ', mass_ustar=', mass_tmp_ustar, ', ratio=',&
!                      mass_tmp_ustar/mass_tmp_stress, ustarVersion
!                endif
                
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
