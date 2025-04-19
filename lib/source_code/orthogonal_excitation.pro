PRO orthogonal_excitation
  ;define the kernel
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  particle_size = 24.4924;PBS, [25.2579], 1% Glycerol
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
  ;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  Hx0=7.0;mT/mm
  Hz0=7.0;mT/mm

  N_tp = 1000;
  frequency = 20.; kHz
  delta_t = 0.0001 ; ms
  t_array = FINDGEN(N_tp)*delta_t
  Hx_array = Hx0 * SIN(t_array*frequency*2*!PI)
  iplot, t_array*10., Hx_array, xtitle='time [ms]', ytitle='Hx field'

  H_total_array = SQRT(Hx_array^2 + Hz0^2)
  iplot, t_array*10., H_total_array, xtitle='time [ms]', ytitle='H total field'

  M_array = Msat_kAm*(1./TANH(beta_particle*H_total_array) - 1./(beta_particle*H_total_array))
  iplot, t_array*10., M_array, xtitle='time [ms]', ytitle='M total'
  ;iplot, H_total_array, M_array, xtitle='H [mT]', ytitle='M total'

  Mz_array = M_array * Hz0 / H_total_array
  iplot, t_array*10., Mz_array, xtitle='time [ms]', ytitle='M z'

  signal_z_array = Mz_array[1:N_tp-1] - Mz_array[0:N_tp-2]
  iplot, t_array[0:N_tp-2]*10., signal_z_array, xtitle='time [ms]', ytitle='signal z'

END