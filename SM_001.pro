PRO SM_001
  
  Time = findgen(500)/100.     ;(1  ms)
  
  G_x = 2.0 ; mT/mm
  FOV = 20 ;mm
  DC_array = (findgen(6)-3.) * (FOV/6.) * G_x
  
  signal_FFT_3rd_array = DBLARR(6)*0.
  
  FOR i=0,5 DO BEGIN
    
    A = 10.; mT
    f = 1. ; kHz
    H = -A*cos(2.*!PI*f*Time) + DC_array[i] ;(mT)

    Msat_T = 0.551 ; Time/u0
    Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
    particle_size = 30. ;nm
    beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
    M = Msat_kAm*(1/tanh(beta*H) - 1/(beta*H)) ;(kA/m)

    signal = -(M[1:*] - M[0:498])/(Time[1]-Time[0])

    u = ABS(FFT(signal))
    t_ele = N_Elements(signal)
    delta_t = Time[1]-Time[0] ; us; samplying rate: 1MHz
    X = (FINDGEN((t_ele - 1)/2) + 1)
    is_N_even = (t_ele MOD 2) EQ 0
    if (is_N_even) then $
      frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
    else $
      frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)

    signal_FFT_3rd_array[i] = u[15]
    
  ENDFOR
  
  iplot, signal_FFT_3rd_array

 
 

END

;@GJ, 2023/9/7
PRO SM_002
  Time = findgen(500)/100.     ;(1  ms)
  G_x = 2.0 ; mT/mm
  FOV = 20 ;mm
  G_y = 2.0
  DC_array_x = DBLARR(6,6)*0.
  FOR i=0,5 DO DC_array_x[i,*] = (i-3.) * (FOV/6.) * G_x
  DC_array_y = DBLARR(6,6)*0.
  FOR j=0,5 DO DC_array_y[*,j] = (j-3.) * (FOV/6.) * G_y
  
  SM_FFT_3rd_array = DBLARR(6,6)*0.
  FOR i=0,5 DO BEGIN
    FOR j=0,5 DO BEGIN
      A = 10.; mT
      f = 1. ; kHz
      H_dyn = -A*cos(2.*!PI*f*Time)
      H = H_dyn + SQRT(DC_array_x[i,j]^2 + DC_array_y[i,j]^2)

      Msat_T = 0.551 ; Time/u0
      Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
      particle_size = 30. ;nm
      beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
      M = Msat_kAm*(1/tanh(beta*H) - 1/(beta*H)) ;(kA/m)

      signal = -(M[1:*] - M[0:498])/(Time[1]-Time[0])
;      signal_x = signal * (H_dyn + DC_array_x[i,j])/SQRT((H_dyn + DC_array_x[i,j])^2 + DC_array_y[i,j]^2)
;      signal_y = signal * DC_array_y[i,j]/SQRT((H_dyn + DC_array_x[i,j])^2 + DC_array_y[i,j]^2)

      u = ABS(FFT(signal))
      t_ele = N_Elements(signal)
      delta_t = Time[1]-Time[0] ; us; samplying rate: 1MHz
      X = (FINDGEN((t_ele - 1)/2) + 1)
      is_N_even = (t_ele MOD 2) EQ 0
      if (is_N_even) then $
        frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
      else $
        frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)

      SM_FFT_3rd_array[i,j] = u[15]
    ENDFOR
  ENDFOR
  
  iimage, SM_FFT_3rd_array

END