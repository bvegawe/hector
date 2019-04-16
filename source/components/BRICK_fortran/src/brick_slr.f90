!===============================================================================
!  Subroutines to run BRICK
!===============================================================================
!   Description to come...
!===============================================================================
! Copyright 2016 Tony Wong, Alexander Bakker
! This file is part of BRICK (Building blocks for Relevant Ice and Climate
! Knowledge). BRICK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! BRICK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
!===============================================================================

module brick

    USE global
    USE gsic_magicc
    USE simple
    USE brick_te
    USE brick_tee
    USE dais

    implicit none
    private

! variables
    real(DP) :: tstep
    logical  :: use_aisfastdyn
! public subroutines
    public :: brick_step_forward, init_brick

contains

!===============================================================================
subroutine init_brick(nstep, tstep_in, &
                      beta0_gsic_magicc_in, V0_gsic_magicc_in, n_gsic_magicc_in, &
                      Teq_gsic_magicc_in, Gs0_gsic_magicc_in, sl_gsic_init_out, &
                      a_te_in, b_te_in, invtau_te_in, V0_te_in, sl_te_init_out, &
                      c_tee_in, a_tee_in, rho_tee_in, sa_tee_in, luse_tee, &
                      a_simple_in, b_simple_in, alpha_simple_in, beta_simple_in, &
                      V0_simple_in, sl_gis_init_out, vol_gis_init_out, &
                      luse_aisfastdyn, parameters_dais_in, sl_ais_init_out, & 
                      rad_ais_init_out, vol_ais_init_out, disint_ais_init_out, &
                      sl_init_out)
!  =========================================================================
!   Initialize the BRICK parameters and initial variables
!  =========================================================================

    integer(i4b), intent(IN)  :: nstep
    real(DP), intent(IN) :: tstep_in
    !!real(DP), intent(IN) :: temp_forcing_current
    real(DP), intent(IN) :: beta0_gsic_magicc_in
    real(DP), intent(IN) :: V0_gsic_magicc_in
    real(DP), intent(IN) :: n_gsic_magicc_in
    real(DP), intent(IN) :: Teq_gsic_magicc_in
    real(DP), intent(IN) :: Gs0_gsic_magicc_in
    real(DP), intent(IN) :: a_te_in
    real(DP), intent(IN) :: b_te_in
    real(DP), intent(IN) :: invtau_te_in
    real(DP), intent(IN) :: V0_te_in
    real(DP), intent(IN) :: c_tee_in
    real(DP), intent(IN) :: a_tee_in
    real(DP), intent(IN) :: rho_tee_in
    real(DP), intent(IN) :: sa_tee_in
    integer(i4b), intent(IN) :: luse_tee
    real(DP), intent(IN) :: a_simple_in
    real(DP), intent(IN) :: b_simple_in
    real(DP), intent(IN) :: alpha_simple_in
    real(DP), intent(IN) :: beta_simple_in
    real(DP), intent(IN) :: V0_simple_in
    integer(i4b), intent(IN) :: luse_aisfastdyn
    real(DP), dimension(23), intent(IN) :: parameters_dais_in

    real(DP), intent(OUT) :: sl_gsic_init_out
    real(DP), intent(OUT) :: sl_te_init_out
    real(DP), intent(OUT) :: sl_gis_init_out
    real(DP), intent(OUT) :: vol_gis_init_out
    real(DP), intent(OUT) :: sl_ais_init_out
    real(DP), intent(OUT) :: rad_ais_init_out
    real(DP), intent(OUT) :: vol_ais_init_out
    real(DP), intent(OUT) :: disint_ais_init_out
    real(DP), intent(OUT) :: sl_init_out

    real(DP) :: sea_level_noAIS

! Assign values to model parameters, and initialize values
    tstep = tstep_in

! GSIC-MAGICC
    call init_gsic_magicc(tstep, beta0_gsic_magicc_in, V0_gsic_magicc_in, n_gsic_magicc_in, &
                          Teq_gsic_magicc_in, Gs0_gsic_magicc_in, sl_gsic_init_out)

! TE
    IF (luse_tee .EQ. 1) THEN 
        call init_brick_tee( c_tee_in, a_tee_in, rho_tee_in, sa_tee_in, &
                        V0_te_in, sl_te_init_out)
    ELSE
        call init_brick_te( tstep, a_te_in, b_te_in, invtau_te_in, V0_te_in, &
                        sl_te_init_out)
    END IF

! GIS-SIMPLE
    call init_simple(tstep, a_simple_in, b_simple_in, alpha_simple_in, &
                     beta_simple_in, V0_simple_in, vol_gis_init_out)
    sl_gis_init_out = V0_simple_in - vol_gis_init_out

! AIS-DAIS
    ! Sea level relative to 1961-1990 observational average (see ../obs
    ! directory for the calculation).
    sea_level_noAIS = sl_gsic_init_out + sl_te_init_out + sl_gis_init_out - 0.1225

    IF (luse_aisfastdyn .EQ. 1) THEN
        use_aisfastdyn = .TRUE.
    ELSE
        use_aisfastdyn = .FALSE.
    END IF

    call init_dais(tstep, parameters_dais_in, sea_level_noAIS, use_aisfastdyn, &
                   rad_ais_init_out, vol_ais_init_out, disint_ais_init_out )
    sl_ais_init_out = 0.0d0

! GMSL
    sl_init_out = sl_gsic_init_out + sl_te_init_out + sl_gis_init_out + sl_ais_init_out

end subroutine init_brick
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_step_forward(nstep, temp_forcing_previous, delta_ocheat_previous, luse_tee, &
                              sl_gsic_previous, sl_gsic_current, &
                              sl_te_previous, sl_te_current, &
                              sl_gis_previous, vol_gis_previous, sl_gis_current, vol_gis_current, &
                              a_anto, b_anto, slope_Ta2Tg, intercept_Ta2Tg, &
                              sl_ais_previous, rad_ais_previous, vol_ais_previous, & 
                              sl_ais_current, rad_ais_current, vol_ais_current, &
                              disint_ais_current, &
                              sl_previous, sl_current)
!------------------------------------------------------------------------------
! Calculate current state from previous state
! 
! Input:
!  nstep                current time step
!  temp_forcing_previous previous time step temperature forcing (rel to 1850) [deg C]
!  delta_ocheat_previous previous time step ocean heat change [J]
!  luse_tee             whether to use tee (1) or te (0)
!  sl_gsic_previous     previous time step cumulative GSIC contribution to SL [m]
!  sl_te_previous       previous time step cumulative TE contribution to SL [m]
!  vol_gis_previous     previous time step Greenland ice sheet volume [m SLE]
!  sl_gis_previous      previous time step cumulative GIS contribution to SL [m]
!  a_anto               sensitivity of Toc to Tg
!  b_anto               Toc for Tg=0 deg C
!  Tfrz                 freezing temperature of ocean water (deg C)
!  slope_Ta2Tg          slope of the regression of Tg as linear function of Ta
!  intercept_Ta2Tg      intercept of the regression of Tg as linear function of Ta
!  vol_ais_previous     previous time step Antarctice ice sheet volume [m^3]
!  sl_ais_previous      previous time step cumulative AIS contribution to SL [m]
!
! Output:
!  sl_gsic_current      current time step cumulative GSIC contribution to SL [m]
!  sl_te_current        current time step cumulative TE contribution to SL [m]
!  vol_gis_current      current time step Greenland ice sheet volume [m SLE]
!  sl_gis_current       current time step cumulative GIS contribution to SL [m]
!  rad_ais_current      current time step Antarctic ice sheet radius [m]
!  vol_ais_current      current time step Antarctic ice sheet volume [m^3]
!  disint_ais_current   current time step Antarctic ice sheet disintegration [m SLE]
!  sl_ais_current       current time step cumulative AIS contribution to SL [m]
!  sl_current           current time step total sea level [m]
!------------------------------------------------------------------------------

    implicit none

    integer(i4b), intent(IN)  :: nstep
    real(DP), intent(IN)  :: temp_forcing_previous
    real(DP), intent(IN)  :: delta_ocheat_previous
    integer(i4b), intent(IN)  :: luse_tee
    real(DP), intent(IN)  :: sl_gsic_previous
    real(DP), intent(IN)  :: sl_te_previous
    real(DP), intent(IN)  :: sl_gis_previous
    real(DP), intent(IN)  :: vol_gis_previous
    real(DP), intent(IN)  :: a_anto
    real(DP), intent(IN)  :: b_anto
    real(DP), intent(IN)  :: slope_Ta2Tg
    real(DP), intent(IN)  :: intercept_Ta2Tg
    real(DP), intent(IN)  :: sl_ais_previous
    real(DP), intent(IN)  :: rad_ais_previous
    real(DP), intent(IN)  :: vol_ais_previous
    real(DP), intent(IN)  :: sl_previous

    real(DP), intent(OUT) :: sl_gsic_current
    real(DP), intent(OUT) :: sl_te_current
    real(DP), intent(OUT) :: sl_gis_current
    real(DP), intent(OUT) :: vol_gis_current
    real(DP), intent(OUT) :: sl_ais_current
    real(DP), intent(OUT) :: rad_ais_current
    real(DP), intent(OUT) :: vol_ais_current
    real(DP), intent(OUT) :: disint_ais_current
    real(DP), intent(OUT) :: sl_current

    real(DP) :: change_sea_level_noAIS, sea_level_noAIS_previous
    real(DP) :: ctmp, Toc_previous, Ta_previous

! Start the show.

! GSIC-MAGICC
    call gsic_magicc_step_forward(temp_forcing_previous, sl_gsic_previous, sl_gsic_current)

! TE
    IF(luse_tee .EQ. 1) THEN
        call brick_tee_step_forward(delta_ocheat_previous, sl_te_previous, sl_te_current)
    ELSE
        call brick_te_step_forward(temp_forcing_previous, sl_te_previous, sl_te_current)
    END IF

! GIS-SIMPLE
! Temperatures are taken relative to the 1961-1990 observational average (see
! ../obs directory for the calculation).
    call simple_step_forward(temp_forcing_previous - 0.3224, vol_gis_previous, vol_gis_current) !for hadcrut
!    call simple_step_forward(temp_forcing_previous - 0.3209, vol_gis_previous, vol_gis_current) !for gistemp
    sl_gis_current = sl_gis_previous + (vol_gis_previous - vol_gis_current)

! AIS-DAIS
    ! sea-level changes, fingerprinted to Antarctic local
    change_sea_level_noAIS = 1.0d0*(sl_gsic_current - sl_gsic_previous) + &
                             1.0d0*(sl_te_current   - sl_te_previous)   + &
                             1.0d0*(sl_gis_current  - sl_gis_previous)
    ! Sea level relative to 1961-1990 observational average (see ../obs
    ! directory for the calculation).
    sea_level_noAIS_previous = 1.0d0*sl_gsic_previous + &
                               1.0d0*sl_te_previous   + &
                               1.0d0*sl_gis_previous  - &
                               0.1225

    ! scale temperatures, accounting for relative to 1850
    ctmp = (Tfrz-b_anto)/a_anto
    Toc_previous = Tfrz + ((a_anto*temp_forcing_previous + b_anto - Tfrz)/(1.0d0 + DEXP(-temp_forcing_previous + ctmp)))
    Ta_previous = (temp_forcing_previous - intercept_Ta2Tg)/slope_Ta2Tg

    call dais_step( Ta_previous, sea_level_noAIS_previous, Toc_previous, &
                    change_sea_level_noAIS, rad_ais_current, vol_ais_current, &
                    disint_ais_current )
    sl_ais_current = sl_ais_previous + &
                     (57.0d0 - sl_ais_previous)*(1.0d0 - (vol_ais_current/vol_ais_previous))
    
! GMSL
    sl_current = sl_gsic_current + sl_te_current + sl_gis_current + sl_ais_current

end subroutine brick_step_forward
!------------------------------------------------------------------------------

END MODULE brick