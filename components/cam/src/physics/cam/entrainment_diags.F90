module entrainment_diags

  ! integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
  use shr_kind_mod,  only: r8 => shr_kind_r8

  implicit none
  
  private
  
  public  :: entrainment_diags_eam
  
contains

subroutine entrainment_diags_eam (pver, &
     dp, z, &
     p, t, qv, ql, &
     dtdt, dqvdt, dqldt, &
     ! dtdt_adv, dqvdt_adv, dqldt_adv, &
     dtdt_radheat, &
     rain_prod, &
     hflux_srf, qvflx_srf, &
     pblh, &
     E_theta, E_q )

  use physconst, only: rair, cpair, rga, latvap

  integer,  intent(in)	:: pver							! number of model levels
  real(r8), intent(in)  :: dp(pver), z(pver)					! pressure thickness (Pa), height (m)
  real(r8), intent(in)  :: p(pver), t(pver), qv(pver), ql(pver)			! pressure (Pa), temperature (K), water vapor, liquid water mixing ratios (kg kg^-1)
  real(r8), intent(in)  :: dtdt(pver), dqvdt(pver), dqldt(pver)	! total tendencies thereof (multiply previous line by s^-1)
  ! real(r8), intent(in)  :: dtdt_adv(pver), dqvdt_adv(pver), dqldt_adv(pver)	! advective tendencies from dycore (same units as previous line)
  real(r8), intent(in)  :: dtdt_radheat(pver)					! total radiative heating (K s^-1)
  real(r8), intent(in)  :: rain_prod(pver)					! warm rain cloud water sink (kg kg^-1 s^-1)
  real(r8), intent(in)  :: hflux_srf, qvflx_srf					! surface sensible heat (W m^-2) and water vapor (kg kg^-1 m^-2) flux
  real(r8), intent(in)  :: pblh							! PBL top height (m)
  real(r8), intent(out) :: E_theta, E_q						! entrainment mass flux (kg m^-2 s^-1) diagnosed from theta_l and q_t budget equations

  !! working variables
  real(r8) :: theta(pver), theta_l(pver), q_t(pver)
  real(r8) :: dtheta_dt(pver), dtheta_l_dt(pver), dq_t_dt(pver)
  integer :: ipbl
  
  real(r8), parameter :: p0 = 1.0e5_r8 ! reference pressure for potential temperature
  
  !! calculate the (conserved) budget variables 
  theta(:) = t(:) * ((p0 / p(:)) ** (rair/cpair))
  theta_l(:) = theta(:) - latvap / cpair * ql(:)
  q_t(:) = qv(:) + ql(:)
  !! and their tendencies
  dtheta_dt(:) = dtdt(:) * ((p0 / p(:)) ** (rair/cpair))
  dtheta_l_dt(:) = dtheta_dt(:) - latvap / cpair * dqldt(:)
  dq_t_dt(:) = dqvdt(:) + dqldt(:)
  
  !! find inversion level
  ipbl = minloc(abs(z - pblh), dim = 1)

  !! call the actual entrainment code
  call entrain_diag_actual(rga, cpair, ipbl, pver, pver, &
       p, theta_l, q_t, &
       dtheta_l_dt, dq_t_dt, &
       dp, &
       dtdt_radheat, rain_prod, &
       hflux_srf, qvflx_srf, &
       E_theta, E_q )
  
end subroutine entrainment_diags_eam

subroutine entrain_diag_actual (rga, cpair, ipbl, isfc, pver, &
     p, theta_l, q_t, &
     dtheta_l_dt, dq_t_dt, &
     dp, &
     radheat, rain_prod, &
     shflx, qflx, &
     E_theta, E_q )

  real(r8), intent(in)  :: rga, cpair						! 1/g, c_p of dry air
  integer,  intent(in)  :: ipbl, isfc						! index of PBL top, surface
  integer,  intent(in)	:: pver							! number of model levels
  real(r8), intent(in)  :: p(pver), theta_l(pver), q_t(pver)			! pressure (Pa), temperature (K), water vapor, liquid water mixing ratios (kg kg^-1)
  real(r8), intent(in)  :: dtheta_l_dt(pver), dq_t_dt(pver)			! total tendencies thereof (multiply previous line by s^-1)
  real(r8), intent(in)  :: dp(pver)						! pressure thickness (Pa)
  ! real(r8), intent(in)  :: dtdt_adv(pver), dqvdt_adv(pver), dqldt_adv(pver)	! advective tendencies from dycore (same units as previous line)
  real(r8), intent(in)  :: radheat(pver)					! total radiative heating (K s^-1)
  real(r8), intent(in)  :: rain_prod(pver)					! warm rain cloud water sink (kg kg^-1 s^-1)
  real(r8), intent(in)  :: shflx, qflx
  real(r8), intent(out) :: E_theta, E_q						! entrainment mass flux (kg m^-2 s^-1) diagnosed from theta_l and q_t budget equations

  !! working variables
  real(r8) :: p_pbl, p_sfc		! PBL top and lowest model level pressure
  real(r8) :: theta_l_plus, q_t_plus, theta_l_0, q_t_0, delta_theta, delta_q
  real(r8) :: rhoH, rhoH_dtheta_l_0_dt, rhoH_dq_t_0_dt
  real(r8) :: delta_F_cp, sh_cp
  real(r8) :: delta_R
  real(r8) :: E_theta_delta_theta, E_q_delta_q
  
  ! if (ipbl == isfc .or. ipbl == 1) then
  ! return error
  p_pbl = p(ipbl)
  p_sfc = maxval(p)
  ! omega_pbl = omega(ipbl)
  ! calculate thermodynamics at the base of the FT
  theta_l_plus = theta_l(ipbl - 1)
  q_t_plus = q_t(ipbl - 1)
  ! calculate PBL averages
  rhoH = sum(dp(ipbl:isfc) * rga) ! FIXME: why does this expression give different results? (p_sfc - p_pbl) / 9.8
  rhoH_dtheta_l_0_dt = sum(dp(ipbl:isfc) * rga * dtheta_l_dt(ipbl:isfc))
  rhoH_dq_t_0_dt = sum(dp(ipbl:isfc) * rga * dq_t_dt(ipbl:isfc))
  theta_l_0 = sum(dp(ipbl:isfc) * rga * theta_l(ipbl:isfc)) / rhoH
  q_t_0 = sum(dp(ipbl:isfc) * rga * q_t(ipbl:isfc)) / rhoH
  delta_F_cp = -sum(dp(ipbl:isfc) * rga * radheat(ipbl:isfc))  ! $\Delta F$ is radiative *cooling* in our weird sign convention
  sh_cp = shflx / cpair
  delta_R = sum(dp(ipbl:isfc) * rga * rain_prod(ipbl:isfc)) ! nonnegative-definite, i.e., must be tendencies of precip, not cloud water
  ! solve budget equations for the entrainment terms
  delta_theta = theta_l_plus - theta_l_0
  E_theta_delta_theta = rhoH_dtheta_l_0_dt + delta_F_cp - sh_cp
  E_theta = E_theta_delta_theta / delta_theta
  delta_q = q_t_plus - q_t_0
  E_q_delta_q = rhoH_dq_t_0_dt - qflx + delta_R
  E_q = E_q_delta_q / delta_q
end subroutine entrain_diag_actual

end module entrainment_diags
