subroutine writeGrid(grid_filename, grid, nx_grid, ny_grid)
    !writeGrid: a subroutine to write a grid of size (nx_landuse, ny_landuse) to grid_filename (simple txt file, mostly for debugging)
    implicit none
    !***************************************************************
    character(*) :: grid_filename
    integer :: nx_grid, ny_grid, ix, iy
    real :: grid(0:nx_grid - 1, 0:ny_grid - 1)
    integer, parameter :: fileunit = 123
    !***************************************************************

    !Open the output file
    !***************************************************************
    open(unit = fileunit, file = grid_filename, action = 'write')

    !loop through grid, write per column
    !***************************************************************
    do ix = 0, nx_grid - 1
        do iy = 0, ny_grid - 1
            write(fileunit, '(E10.2)', advance = 'no') grid(ix, iy)
        end do
        write(fileunit, *)
    end do

    !close output file
    !***************************************************************
    close(fileunit)

end subroutine writeGrid

subroutine writeGridBin(grid_filename, grid, nx_grid, ny_grid)
    !writeGrid: a subroutine to write a grid of size (nx_landuse, ny_landuse) to grid_filename (simple txt file, mostly for debugging)
    use dust_mod
    implicit none

    !***************************************************************
    character(*) :: grid_filename
    integer :: nx_grid, ny_grid
    real :: grid(0:nx_grid - 1, 0:ny_grid - 1)
    integer, parameter :: fileunit = 123
    !***************************************************************

    !Open the output file
    !***************************************************************
    open(unit = fileunit, file = grid_filename, form = 'unformatted', access = 'stream')
    write(fileunit) lat_bottom
    write(fileunit) lon_left
    write(fileunit) dx_dy_out
    write(fileunit) ny_lat_out
    write(fileunit) nx_lon_out
    write(fileunit) time_step
    !write(*,*) 'Time step: ', time_step
    write(fileunit) grid(:,:)

    !close output file
    !***************************************************************
    close(fileunit)

end subroutine writeGridBin

subroutine writeGridInteger(grid_filename, grid, nx_grid, ny_grid)
    !writeGrid: a subroutine to write a grid of size (nx_landuse, ny_landuse) to grid_filename (simple txt file, mostly for debugging)
    implicit none
    !***************************************************************
    character(*) :: grid_filename
    integer :: nx_grid, ny_grid, ix, iy
    integer :: grid(0:nx_grid - 1, 0:ny_grid - 1)
    integer, parameter :: fileunit = 123
    !***************************************************************

    !Open the output file
    !***************************************************************
    open(unit = fileunit, file = grid_filename, action = 'write')

    !loop through grid, write per column
    !***************************************************************
    do ix = 0, nx_grid - 1
        do iy = 0, ny_grid - 1
            write(fileunit, '(I8)', advance = 'no') grid(ix, iy)
        end do
        !print*, grid(ix,:)
        write(fileunit, *)
    end do

    !close output file
    !***************************************************************
    close(fileunit)

end subroutine writeGridInteger


subroutine getGridPointWind(lat, lon, inNestNr, ix_wind, iy_wind, ix_wind_n, iy_wind_n)
    !getGridPoint: subroutine to obtain the gridpoint (ix_wind,iy_wind) in a windfield(inNestNr)
    !that corresponds to the point (lon,lat)
    use par_mod
    use com_mod
    implicit none

    !***************************************************************
    integer :: ix_wind, iy_wind, ix_wind_n, iy_wind_n, i_nest
    real :: lat, lon, nest_max_lat, nest_max_lon
    integer :: inNestNr
    !***************************************************************

    ix_wind_n = 0
    iy_wind_n = 0
    ix_wind = 0
    iy_wind = 0


    !First test if point is inside one of the nested grids
    !***************************************************************
    inNestNr = 0


    do i_nest = 1, numbnests
        nest_max_lat = ylat0n(i_nest) + nyn(i_nest) * dyn(i_nest)
        nest_max_lon = xlon0n(i_nest) + nxn(i_nest) * dxn(i_nest)

        if (lat .ge. ylat0n(i_nest) .and. lat .lt. nest_max_lat) then
            if (lon .ge. xlon0n(i_nest) .and. lon .lt. nest_max_lon) then
                !point is inside nested grid number 'inNestNr'
                !*****************************************************
                inNestNr = i_nest
                ix_wind_n = (lon - xlon0n(i_nest))/dxn(i_nest)
                iy_wind_n = (lat - ylat0n(i_nest))/dyn(i_nest)

                !If array boundaries are exceeded point was outside nest anyway
                !*****************************************************
                if (ix_wind_n .lt. 0) inNestNr = 0
                if (ix_wind_n .gt. nxn(i_nest)) inNestNr = 0
                if (iy_wind_n .lt. 0) inNestNr = 0
                if (iy_wind_n .gt. nyn(i_nest)) inNestNr = 0
                !  if(inNestNr.eq.1)then
                !      write(*,*) lat,lon,ix_wind_n, iy_wind_n, ustarn(ix_wind_n,iy_wind_n,1,1,1), uun(ix_wind_n,iy_wind_n,1,1,1)
                !  endif
            endif
        endif
    end do
    !***************************************************************

    !Get ix and iy in normal wind field
    !***************************************************************
    ix_wind = (lon - xlon0)/dx
    iy_wind = (lat - ylat0)/dy
    !check that array boundaries are not exceeded
    !***********************************************************
    if (ix_wind .lt. 0) ix_wind = ix_wind + (nx - 1)
    if (ix_wind .gt. nx - 1) ix_wind = ix_wind - (nx - 1)
    if (iy_wind .lt. 0) iy_wind = iy_wind + (ny - 1)
    if (iy_wind .gt. ny - 1) iy_wind = iy_wind - (ny - 1)

end subroutine getGridPointWind

subroutine getGridPointClay(lat, lon, inClay, ix_clay, iy_clay)
    use dust_mod
    implicit none
    integer :: ix_clay, iy_clay, ix_clay_tmp, iy_clay_tmp
    real :: lat, lon
    logical :: inClay

    inClay = .true.
    ix_clay_tmp = (lon - xlon0_c)/dx_c
    iy_clay_tmp = (lat - ylat0_c)/dy_c
    !write(*,*) iy_clay_tmp, lat, ylat0_c, dy_c
    !if array boundaries are exceeded point is outside grid with claycontent
    !***********************************************************************
    if (ix_clay_tmp .lt. 0 .or. ix_clay_tmp .gt. nx_c - 1 .or. iy_clay_tmp .lt. 0 .or. iy_clay_tmp .gt. ny_c - 1) then
        inClay = .false.
        !don't change ix&iy
    else
        !save location x and y
        ix_clay = ix_clay_tmp
        ix_clay = min(ix_clay, nx_c - 1)
        iy_clay = iy_clay_tmp
    endif
    !write(*,*) 'in clay: ',lat, lon, ix_clay, iy_clay, inClay
end subroutine getGridPointClay


subroutine getGridPointLandUse(lat, lon, ix_lu, iy_lu, dxdy_degr_landuse)
    implicit none
    integer :: ix_lu, iy_lu
    real :: lat, lon, dxdy_degr_landuse

    ix_lu = (lon - (-180))/dxdy_degr_landuse
    iy_lu = (lat - (-90))/dxdy_degr_landuse

    !write(*,*) lat, lon, ix_lu,iy_lu
end subroutine getGridPointLandUse

subroutine getGridPointLandUse_n(lat, lon, ix_lu_n, iy_lu_n, inLU_n)
    use dust_mod
    use par_mod
    implicit none
    integer :: ix_lu_n, iy_lu_n
    real :: lat, lon
    logical :: inLU_n

    inLU_n = .false.

    if (lon .ge. xlon0_n(1).and. lon .le. (xlon0_n(1) + nx_landuse_n(1) * dxdy_landuse_n(1)) &
        .and. lat .ge. ylat0_n(1) .and. lat .le. (ylat0_n(1) + ny_landuse_n(1) * dxdy_landuse_n(1)))then

        inLU_n = .true.

        ix_lu_n = (lon - xlon0_n(1))/dxdy_landuse_n(1)
        iy_lu_n = (lat - ylat0_n(1))/dxdy_landuse_n(1)
        !print*, lat, lon, ix_lu_n, iy_lu_n

        if (iy_lu_n .gt. ny_landuse_n(1))then
            inLU_n = .false.
        endif

        if (ix_lu_n .gt. nx_landuse_n(1))then
            inLU_n = .false.
        endif

    else
        inLU_n = .false.

    endif
end subroutine getGridPointLandUse_n

subroutine getGridPointErosionClass(lat, lon, ix, iy, inGrid)
    use dust_mod
    use par_mod
    implicit none
    integer :: ix, iy
    real :: lat, lon
    logical :: inGrid

    inGrid = .false.

    if (lon .ge. xlon0_erC.and. lon .le. (xlon0_erC + nx_erC * dxdy_erC) &
        .and. lat .ge. ylat0_erC .and. lat .le. (ylat0_erC + ny_erC * dxdy_erC))then

        inGrid = .true.

        ix = (lon - xlon0_erC)/dxdy_erC
        iy = (lat - ylat0_erC)/dxdy_erC


        if (iy .gt. ny_erC)then
            inGrid = .false.
        endif

        if (ix .gt. nx_erC)then
            inGrid = .false.
        endif

    else
        inGrid = .false.

    endif
end subroutine getGridPointErosionClass

subroutine getGridPoints(ix_lu, iy_lu, inLU_n, ix_lu_n, iy_lu_n, ix_wind, iy_wind, inClayGrid, ix_clay, iy_clay, &
    inNestNr, ix_wind_n, iy_wind_n, inErosionClGrid, ix_erClass, iy_erClass)
    use dust_mod
    use com_mod
    implicit none
    integer :: ix, iy !ix and iy in output grid
    real    :: lat_out, lon_out
    logical :: inLU_n(0:nx_lon_out - 1, 0:ny_lat_out - 1), inClayGrid(0:nx_lon_out - 1, 0:ny_lat_out - 1) !to check if a point in the output grid is available in the landuse and clay grid (not global)
    logical :: inErosionClGrid(0:nx_lon_out - 1, 0:ny_lat_out - 1)  !to check if a point in the output grid is available in the erosion class grid(not global)
    integer :: ix_lu_n(0:nx_lon_out - 1, 0:numbnests_landuse) !x-coordinate in land use nest, started to prepare for multiple nested fields but never tested
    integer :: iy_lu_n(0:ny_lat_out - 1, 0:numbnests_landuse)
    integer :: ix_lu(0:nx_lon_out - 1), ix_wind(0:nx_lon_out - 1), ix_clay(0:nx_lon_out - 1), ix_erClass(0:nx_lon_out - 1) !Gives x-coordinates in landuse/wind/clay/erosion grid for coordinate ix in output grid
    integer :: iy_lu(0:ny_lat_out - 1), iy_wind(0:ny_lat_out - 1), iy_clay(0:ny_lat_out - 1), iy_erClass(0:ny_lat_out - 1)
    integer :: inNestNr(0:nx_lon_out - 1, 0:ny_lat_out - 1) !check if outpunt point is in nested wind field
    integer :: tmp_inNestNr, tmp_ix_wind_n, tmp_iy_wind_n
    integer :: ix_wind_n(0:nx_lon_out - 1, 0:numbnests)
    integer :: iy_wind_n(0:ny_lat_out - 1, 0:numbnests)

    !initialize
    inLu_n(:,:)=.false.
    inNestNr(:,:)=0
    ix_lu_n(:,:)=0
    iy_lu_n(:,:)=0
    
    !Get information on location all grids with respect to output grid
    !*************************************************************
    do iy = 0, ny_lat_out - 1
        do ix = 0, nx_lon_out - 1

            lat_out = lat_bottom + iy * dx_dy_out
            lon_out = lon_left + ix * dx_dy_out

            !Determine position of output point in the several available grids
            call getGridPointLandUse(lat_out, lon_out, ix_lu(ix), iy_lu(iy), dxdy_degr_landuse)
            if (numbnests_landuse .gt. 0)then
                if( numbnests_landuse .eq. 1)then
                    call getGridPointLandUse_n(lat_out, lon_out, iy_lu_n(iy,numbnests_landuse), iy_lu_n(iy, numbnests_landuse), &
                    inLU_n(ix, iy))
                else
                    print*, 'Need to adjust source code grid point land use for more than 1 nests!'
                    stop
                endif
            endif
            call getGridPointWind(lat_out, lon_out, tmp_inNestNr, ix_wind(ix), iy_wind(iy), tmp_ix_wind_n, tmp_iy_wind_n)
            call getGridPointClay(lat_out, lon_out, inClayGrid(ix, iy), ix_clay(ix), iy_clay(iy))
            if (applyClassErosion)then
                call getGridPointErosionClass(lat_out, lon_out, ix_erClass(ix), iy_erClass(iy), inErosionClGrid(ix, iy))
            endif

            inNestNr(ix, iy) = tmp_inNestNr
            if (tmp_inNestNr .gt. 0)then
                ix_wind_n(ix, tmp_inNestNr) = tmp_ix_wind_n
                iy_wind_n(iy, tmp_inNestNr) = tmp_iy_wind_n
            else
                ix_wind_n(ix, tmp_inNestNr) = 0
                iy_wind_n(iy, tmp_inNestNr) = 0
            endif
        end do
    end do
    !*************************************************************
end subroutine getGridPoints

subroutine writeSummary(totalEmission)
    use dust_mod
    implicit none

    integer, parameter :: fileunit = 123
    real :: totalEmission

    !Open the output file
    !***************************************************************
    open(unit = fileunit, file = summary_file, action = 'write')

    !****************************************************************
    write(fileunit, *) 'INPUT'
    write(fileunit, *) 'Wind fields:', ECMWF_input
    write(fileunit, *) 'Nested wind fields:', ECMWF_input_nest
    write(fileunit, *) 'Time step wind fields:', time_step_wind
    write(fileunit, *) 'Land use file:', landuse_file
    write(fileunit,  *) 'Nested landuse: ',landuse_file_n

    write(fileunit, *) 'OUTPUT'
    write(fileunit, *) 'Modelled period:', start_date_day, start_date_hour, ' +', releaseDays, ' days'
    write(fileunit, *) 'Time resolution:', time_step
    write(fileunit, *) 'Grid location:', lat_bottom, lon_left
    write(fileunit, *) 'Grid resolution:', dx_dy_out
    write(fileunit, *) 'Grid dimension:', ny_lat_out, nx_lon_out

    write(fileunit, *) 'MODEL SETTINGS'
    write(fileunit, *) 'Mobility threshold:', mobThreshold
    write(fileunit, *) 'Size distribution:', typeSizeDistr
    write(fileunit, *) 'Junge index:', Junge_index
    write(fileunit, *) 'Scaling factor:', scalingFactor
    write(fileunit, *) 'Emission scheme:', emissionModel
    
    write(fileunit, *) 'Topographic control erodibility: ', EROSION_TOPO
    write(fileunit, *) 'Precipitation blocking mobil.:',PRECIP_BLOCK
    write(fileunit, *) 'Soil moisture:', SOILMOSITURE_DEP
    write(fileunit, *) 'Corr. land/sea:', correctLSM_SNOW

    write(fileunit, *) 'Total emitted mass:', totalEmission
    close(fileunit)
end subroutine writeSummary
