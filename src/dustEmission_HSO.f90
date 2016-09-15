subroutine dustEmmission_HSO(em, soilF, lat,dxdy_degr_landuse, inNestNr, ix_wind, iy_wind,  ix_wind_n, iy_wind_n,&
                                time_step, scalingFactor, mobilisationThreshold)

use par_mod
use com_mod

implicit none

real     :: mfac, em, soilF, lat, scalingFactor, mass, mflux, mobilisationThreshold, vel, dxdy_degr_landuse
integer  :: ix_wind, iy_wind,  ix_wind_n, iy_wind_n, time_step, inNestNr

        !Calculate velocity (vel) from 10 m u- and v-component of the wind.
        !******************************************************************
        if(inNestNr.eq.0)then
                !normal field
                vel=sqrt(u10(ix_wind,iy_wind,1,1)*u10(ix_wind,iy_wind,1,1)&
                        +v10(ix_wind,iy_wind,1,1)*v10(ix_wind,iy_wind,1,1))
        else
                !nested field
                vel=sqrt(u10n(ix_wind_n,iy_wind_n,1,1,inNestNr)*u10n(ix_wind_n,iy_wind_n,1,1,inNestNr)&
                        +v10n(ix_wind_n,iy_wind_n,1,1,inNestNr)*v10n(ix_wind_n,iy_wind_n,1,1,inNestNr))
        endif
        !******************************************************************
        !Calculate emitted mass (em) per time step if the wind velocity exceeds threshold velocity
        !******************************************************************	
	if(vel.ge.mobilisationThreshold)then	
		
		!multiplication factor accounting for time step, area and conversion to kg
		!**********************************************************
		mfac=3600.0*time_step*dxdy_degr_landuse*111.0e3*&
				dxdy_degr_landuse*cos(pi/180.0*lat)*111.0e3/1e9	
				
		!Emitted mass in kg when mfac is used (otherwise mass flux in mug m-2 s-2)
		mass=scalingFactor*(vel-mobilisationThreshold)*vel*vel*soilF*mfac
		mflux=scalingFactor*(vel-mobilisationThreshold)*vel*vel*soilF
		em=em+mass

	endif
        !******************************************************************	
end subroutine dustEmmission_HSO