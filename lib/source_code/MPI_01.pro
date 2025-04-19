FUNCTION LANGEVIN_L, alpha
  ; 自定义郎之万函数计算，处理小alpha的数值稳定性
  threshold = 1e-3  ; 定义泰勒展开的阈值
  result = alpha * 0  ; 初始化输出数组

  ; 分离小alpha和大alpha的情况
  small = WHERE(ABS(alpha) LT threshold, count)
  IF count GT 0 THEN BEGIN
    x = alpha[small]
    ; 泰勒展开：L(α) ≈ α/3 - α^3/45 + 2α^5/945
    result[small] = x/3 - x^3/45 + 2*x^5/945
  ENDIF

  large = WHERE(ABS(alpha) GE threshold AND alpha NE 0, count_large)
  IF count_large GT 0 THEN BEGIN
    x = alpha[large]
    ; 直接计算coth(α) - 1/α
    coth = COSH(x)/SINH(x)  ; IDL内置的双曲函数
    result[large] = coth - 1/x
  ENDIF

  RETURN, result
END

PRO MPI_LANGEVIN_SIMULATION

  ; 常数定义
  mu0 = 4 * !PI * 1e-7       ; 真空磁导率
  kB = 1.38e-23              ; 玻尔兹曼常数
  T = 300                    ; 温度（K）
  Msat = 4.8e5               ; 饱和磁化强度（A/m，磁铁矿）

  ; 输入三种粒子直径（单位：米）
  diameters = [10e-9, 20e-9, 30e-9]

  ; 生成对称的磁场范围（-1e5 到 1e5 A/m）
  H = (FINDGEN(200) - 99.5) * 1e3  ; 200点对称分布

  ; 预存磁化曲线数据
  M = FLTARR(N_ELEMENTS(H), N_ELEMENTS(diameters))

  ; 计算各粒子的磁化强度
  FOR i = 0, N_ELEMENTS(diameters)-1 DO BEGIN
    d = diameters[i]
    ms = (!PI/6) * d^3 * Msat        ; 磁矩计算
    alpha = (mu0 * ms * H) / (kB * T)
    M[*,i] = LANGEVIN_L(alpha)      ; 调用自定义函数
  ENDFOR

  ; 绘制磁化曲线
  WINDOW, TITLE='Magnetization Curves (Corrected)'
  PLOT, H, M[*,0], XTITLE='H (A/m)', YTITLE='M/M_s', COLOR='ffffff'x;, /YNZERO
  OPLOT, H, M[*,1], COLOR='00ffff'x;, NAME='20nm'
  OPLOT, H, M[*,2], COLOR='ff00ff'x;, NAME='30nm'
;  LEGEND, POSITION=[0.1,0.9]

  ; 导出数据到文本文件（可选）
  ; WRITE_CSV, 'magnetization.csv', H, M, HEADER='H, M_10nm, M_20nm, M_30nm'

END



PRO doubao_MPI_Langevin

;  以下是一个IDL程序示例，用于展示朗之万函数在MPI（磁粒子成像）中对图像分辨率影响的模拟。这个示例假设你已经了解MPI的基本原理以及朗之万函数在其中的作用。

;  idl

  ; 定义图像尺寸
  nx = 256
  ny = 256

  ; 生成网格
  x = findgen(nx) - floor(nx / 2)
  y = findgen(ny) - floor(ny / 2)
;  xx, yy = meshgrid(x, y)
  
  xx = DBLARR(nx, ny) * 0.
  yy = DBLARR(nx, ny) * 0.
  FOR i=0, nx-1 DO BEGIN
    FOR j=0, ny-1 DO BEGIN
      xx[i, j] = x[i]
      yy[i, j] = y[j]
    ENDFOR
  ENDFOR
  
  ; 假设的粒子浓度分布（简单高斯分布）
  sigma = 20
  particle_density = exp(-(xx^2 + yy^2) / (2 * sigma^2))

  ; 定义磁场强度范围
  B_max = 1.0
  B = findgen(100) * B_max / 99.0

  ; 玻尔兹曼常数和温度（假设值）
  k_B = 1.38e - 23
  T = 300.0

  ; 磁矩（假设值）
  mu = 1.0e - 23

  ; 计算朗之万函数
  L = 1./TANH(mu * B / (k_B * T)) - (k_B * T) / (mu * B)
  ;deal with finite
  index_nan = where(~FINITE(L), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(L, B, t[index_nan], /NAN)
    L[index_nan] = RESULT
  ENDIF

  ; 定义不同的磁场梯度（模拟对分辨率的影响）
  gradients = [0.1, 0.5, 1.0]

  ; 循环计算不同磁场梯度下的“成像结果”
  for i = 0, n_elements(gradients) - 1 do begin
    gradient = gradients[i]
    ; 模拟磁场分布
    B_field = gradient * sqrt(xx^2 + yy^2)
    ; 根据朗之万函数计算响应
    response = particle_density * INTERPOL(B_field, sqrt(xx^2 + yy^2), L)
    ; 显示结果
    device, DECOMPOSED = 0
    window, 1, TITLE = 'MPI Image with Gradient=' + string(gradient)
    tv, response
  endfor
;
;  代码说明：
;
;  图像设置：
;
;  - 定义图像的尺寸  nx  和  ny ，并生成对应的二维网格  xx  和  yy 。
;
;  粒子浓度分布：
;
;  - 假设粒子浓度分布为高斯分布，使用  sigma  控制分布的宽度。
;
;  朗之万函数计算：
;
;  - 定义磁场强度范围  B 。
;
;  - 设置玻尔兹曼常数  k_B 、温度  T  和磁矩  mu  等物理参数。
;
;  - 根据朗之万函数的定义  L = coth(x) - 1/x （其中  x = mu * B / (k_B * T) ）计算朗之万函数值。
;
;  磁场梯度与成像模拟：
;
;  - 定义不同的磁场梯度  gradients 。
;
;  - 对于每个磁场梯度，模拟磁场分布  B_field 。
;
;  - 根据粒子浓度分布和朗之万函数的响应，计算最终的“成像结果”  response 。
;
;  - 使用  tv  函数显示不同磁场梯度下的成像结果，以观察朗之万函数通过磁场梯度对图像分辨率的影响。
;
;  请注意，这是一个简化的模拟，实际的MPI系统涉及更多复杂的物理和信号处理过程。你可能需要根据具体的研究需求进一步调整和扩展代码。
;  
END


PRO MPI_01

  t = 4*!PI/100 * FINDGEN(100)
  ch1 = sin(t)
  ch2 = sin(t+!PI/2.)
  
  iplot, t, ch1+ch2
  
  plot, t, ch1+ch2
  oplot, t, ch1
  oplot, t, ch2
END

PRO MPI_02
  t = 2*!PI/1000 * FINDGEN(1000)*5.
  ;y = 2*!PI/100 * FINDGEN(100)
  
  ch1=sin(t)
  ch2=sin(18./13.*t)
  
  plot, ch1, ch2
END

PRO MPI_03

  fy = 98000. ;Hz
  fz = 99000. ;Hz
  
  t = FINDGEN(10000)/10000./fy * 30.
  
  y = sin(2*!PI*fy*t)
  z = sin(2*!PI*fz*t)
  
  plot, y, z
END

;@GJ, 2022/11/30, system matrix
PRO MPI04_system_matrix
  
  ;CT projection
  A=[[1.,0.,1.,0.],[0.,1.,0.,1.],[1.,1.,0.,0.],[0.,0.,1.,1.]]
  x = [[3.],[2.],[5.],[10.]]
  b=[[8.],[12.],[5.],[15.]]
  A_invert = INVERT(A, /DOUBLE)
  rec_image = A_invert ## b
  print, 'invert: ', rec_image
  rec_image_1d = art_func(A,b)
  print, 'ART: ', rec_image_1d
  
  ;MRI encoding
  A=[[1.,1.,1.,1.],[-1.,0.,-1.,0.],[-1.,-1.,0.,0.],[1.,1.,-1.,-1.]]
  b=[[20.],[-8.],[-5.],[-10.]]
  A_invert = INVERT(A, /DOUBLE)
  rec_image = A_invert ## b
  print, 'invert: ', rec_image ;good result
  rec_image_1d = art_func(A,b)
  print, 'ART: ', rec_image_1d
  
  ;MPI FFL encoding
  Angle=[[0.,0.,0.,0.],[10.,20.,50.,100.],[20.,40.,100.,200.],[30.,60.,150.,300.]]
  A = COS(Angle/180.*!PI)
  b = A ## x
  A_invert = INVERT(A, /DOUBLE)
  rec_image = A_invert ## b
  print, 'invert: ', rec_image ;good result
  rec_image_1d = art_func(A,b)
  print, 'ART: ', rec_image_1d
  
  ;@GJ, 2022/12/1, extrapolation of system matrices
  A=[[1.,0.,1.,0.,0.],[0.,1.,0.,1.,0.],[1.,1.,0.,0.,1.],[0.,0.,1.,1.,0.],[0.,0.,0.,0.,0.]]
  x_est = [[3.],[2.],[5.],[10.],[20.]]
  img = A ## x_est
  b=[[8.],[12.],[25.],[15.],[0.]]
  rec_image_1d = art_func(A,b)
  print, rec_image_1d
  print, INVERT(A) ## b
  
  
  ;calculate the invert matrix of the signal Kernal
  ;  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  Gradient_invert = INVERT(signalKernal_2d_array, /DOUBLE)
  result = Gradient_invert # signalKernal_2d_array
  iimage, result, title='Unit invert matrix *'

  ;GJ, 2021/9/5, do calculation based ART method
  ;rec_image_1d_art = ART_cal(signal_3rdHar_1d, signalKernal_2d_array)
  ;GJ, 2021/9/25, find solution using Liang Xiaofeng's code
  ;A=[[1,0,1,0],[0,1,0,1],[1,1,0,0],[0,0,1,1]]
  ;b=[[8],[12],[5],[15]]
  ;rec_image_1d = art_func(A,b)
  ;  image = DBLARR(256) * 0.
  ;  image[100:102, 127:129] = 256
  ;  image[130:132, 127:129] = 256
  ;  signal_1d = signalKernal_2d_array # REFORM(image[*, x_size/2])
  ;  rec_image_1d = art_func(signalKernal_2d_array,signal_1d)
  ;end of test
  
END

PRO Langevin_Function_simple
  H_array=(FINDGEN(201)/10-10.)*2.-0.1
  
  ;plot, H_array, 1./(TANH(H_array))
  iplot, H_array, -1./H_array, PSYM=2, yrange=[-15,15]
  iplot, H_array, 1./(TANH(H_array))-1./H_array, PSYM=3, XTITLE='H [mT]', YTITLE='M', TITLE='Langevin Function';, yrange=[-15,15]
END

;@GJ, 2022/10/27, real langevin function with particle size
PRO Langevin_function_complicated, particle_size
  IF N_ELEMENTS(particle_size) EQ 0 THEN particle_size = 20. ;nm
  particle_size *= 1.0

  H_array=(FINDGEN(201)/10-10.)*2.-0.1
  
  Msat_T = 0.551 ; T/u0
  u0 = 4.0 * !PI * (0.1^4) ; T*m/kA
  Msat_kAm = Msat_T/u0; kA/m
  m_moment = 1./6. * !PI * (particle_size^3) * Msat_kAm ; 10^-27 kA*m2
  k_B = 1.380469 ; 10^-23 J/K
  T_p = 36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = m_moment / (K_B * T_p_kelvin) / (10.^4) ; in kA*m2/J or 1/mT
  M = Msat_kAm * (1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
  
  iplot, H_array, M, PSYM=3, XTITLE='H [mT]', YTITLE='M', XRANGE=[-20,20], TITLE='Langevin Function '+STRTRIM(STRING(particle_size),1)+' nm'
END

;@GJ, 2023/1/16, calculating Bc
PRO Bc_estimation, Bc, hydro_particle_size, w0

  IF N_ELEMENTS(particle_size) EQ 0 THEN core_particle_size = 10.;20. ;nm
  core_particle_size *= 1.0
  IF core_particle_size GT hydro_particle_size THEN core_particle_size = hydro_particle_size
  
  ;assume the core and hydro are the same
  ;@GJ, 2023/1/16, check the assumption
  ;core_particle_size = hydro_particle_size
  
  Msat_T = 0.551 ; T/u0
  u0 = 4.0 * !PI * (0.1^4) ; T*m/kA
  Msat_kAm = Msat_T/u0; kA/m
  m_moment = 1./6. * !PI * (core_particle_size^3) * Msat_kAm ; 10^-27 kA*m2
  
  k_B = 1.380469 ; 10^-23 J/K
  T_p = 36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = m_moment / (K_B * T_p_kelvin) / (10.^4) ; in kA*m2/J or 1/mT
  
  IF N_ELEMENTS(tao_B) EQ 0 THEN tao_B = 20. ;us
  basic_hydro_particle_size = 20. ;nm
  tao_B = tao_B / (basic_hydro_particle_size^3) * (hydro_particle_size^3)
  
  IF N_ELEMENTS(w0) EQ 0 THEN w0 = 60. ;kHz
  w0 *= 1.0
  Bc = 2.0 * (1./beta_particle)*(tao_B/1000.)*w0
  print, 'Bc: ', Bc, ' mT'
  
END

PRO BC_size_calc
  
  N_ELE = 101.
  particle_size_array = FINDGEN(N_ELE)+10.0
  BC_array = DBLARR(N_ELE)*0.
  w0=60. ; kHz
  FOR i=0, N_ELE-1 DO BEGIN
    Bc_estimation, Bc, particle_size_array[i], w0
    BC_array[i] = BC
  ENDFOR
  iplot, particle_size_array, BC_array, PSYM=3, /YLOG, XTITLE='particle size [nm]', YTITLE='Bc [mT]', XRANGE=[0,110], TITLE='Bc: Critical Field at '+STRTRIM(STRING(FLOOR(w0)),1)+' kHz'

END

;@GJ, 2023/1/16, rotational frequency drift with normal particle size distribution
PRO RFD_frequency_range_based_on_size_range
  
  w0=30.;60. ; kHz
  
  N_ele = 1000.
  B_Bc = (FINDGEN(N_ele*1.2+1)/N_ele)
  wD = B_Bc*0.
  FOR i=0, N_ele*1.2 DO BEGIN
    IF B_Bc[i] LT 1.0 THEN wD[i] = w0 * (1. - SQRT(1. - B_Bc[i]^2)) ELSE wD[i]=w0
  ENDFOR
  
  
  ave_particle_size = 20.
  stdev_particle_size = 0.5
  Bc_estimation, Bc_max, ave_particle_size+stdev_particle_size, w0
  Bc_estimation, Bc_min, ave_particle_size-stdev_particle_size, w0
  Bc_estimation, Bc_mean, ave_particle_size, w0
  
  window, 1
  plot, B_Bc*Bc_max, wD, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[0., w0*1.1], xtitle='B [mT]', ytitle='wD [kHz]', title='Rotational Frequency Drift at '+STRTRIM(STRING(FLOOR(w0)),1)+' kHz'
  oplot, B_Bc*Bc_min, wD, LINESTYLE=5, COLOR = '0000FF'x;, COLOR='red'
  oplot, B_Bc*Bc_mean, wD, LINESTYLE=3, COLOR = '00FF00'x;, COLOR='red'
  
  iplot, B_Bc*Bc_min, wD, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[0., w0*1.1], xtitle='B [mT]', ytitle='wD [kHz]', title='Rotational Frequency Drift at '+STRTRIM(STRING(FLOOR(w0)),1)+' kHz'
  iplot, B_Bc*Bc_max, wD, BACKGROUND = 'FFFFFF'x, COLOR = 0, XRANGE=[0,25], YRANGE=[0., w0*1.1], xtitle='B [mT]', ytitle='wD [kHz]', title='Rotational Frequency Drift at '+STRTRIM(STRING(FLOOR(w0)),1)+' kHz'
  
  window, 2
  FOR i=0, N_ele*1.1-1 DO IF B_Bc[i]*Bc_max GT 20. THEN break
  plot, B_Bc*Bc_max, wD, BACKGROUND = 'FFFFFF'x, COLOR = 0, XRANGE=[0., B_Bc[i]*Bc_max], YRANGE=[0, wD[i]], xtitle='B [mT]', ytitle='wD [kHz]', title='Rotational Frequency Drift at '+STRTRIM(STRING(FLOOR(w0)),1)+' kHz'
  oplot, B_Bc*Bc_min, wD, LINESTYLE=5, COLOR = '0000FF'x;, COLOR='red'
  oplot, B_Bc*Bc_mean, wD, LINESTYLE=3, COLOR = '00FF00'x;, COLOR='red'
END

;@GJ, 2022/11/01, rotational frequency drift equation
;input: w0: excitation frequency in kHz
PRO RFD_equation, w0, wD, B_Bc, N_ele
  N_ele = 1000.
  B_Bc = (FINDGEN(N_ele+1)/N_ele)
  
  IF N_ELEMENTS(w0) EQ 0 THEN w0 = 5.; kHz
  wD = w0 * (1. - SQRT(1. - B_Bc^2))
  iplot, B_Bc, wD, title='RFD', xtitle='B/Bc', ytitle='wD [kHz]'
;  iplot, B_Bc, wD, /YLOG, title='RFD', xtitle='B/Bc', ytitle='ALOG wD/w0'
  
  wD_gradient = wD[1:*] - wD[0:N_ele-1]
;  iplot, B_Bc[1:*], wD_gradient, xtitle='B/Bc', ytitle='wD/w0 gradient', title='RFD gradient'
;  iplot, B_Bc[1:*], wD_gradient, /YLOG, xtitle='B/Bc', ytitle='ALOG wD/w0 gradient', title='RFD gradient'

  wD_accele = wD_gradient[1:*] - wD_gradient[0:N_ele-2]
;  iplot, B_Bc[2:*], wD_accele, /YLOG, xtitle='B/Bc', ytitle='wD/w0 acceleration', title='RFD acceleration'


END

PRO RFD_gauss_FID
  
  w0 = 50.; kHz
  B_Bc = 0.33;0.85;0.33
  wD = w0 * (1. - SQRT(1. - B_Bc^2))
  
  N = 10000.
  wD_array = DBLARR(N)*0. + wD
  seed = !NULL
  Result = (RANDOMU(seed, N) - 0.5) * 2. 
  wD_array = wD_array + Result * WD * 0.05
;  iplot, wD_array

  t_ele = 1./w0 * 1000. * 100.
  delta_t = 5.0 ;us
  t_array = FINDGEN(t_ele) * delta_t   ;us
  signal_array = DBLARR(t_ele)*0.
  noise_level = 100.
  
  FOR t_ind=0, t_ele-1 DO BEGIN
    FOR j=0, N-1 DO BEGIN
      ;diffusion
      seed = !NULL
      noise_angle = (RANDOMU(seed)-0.5)/90.*!PI*noise_level
      signal_t_j = 1.0 * cos(2. * !PI * wD_array[j] * t_array[t_ind] / 1000. + noise_angle)
      signal_array[t_ind] += signal_t_j
    ENDFOR
  ENDFOR
  
;  iplot, t_array/1000., signal_array, xtitle='time [ms]', ytitle='M'
    
END

;@GJ, 2023/3/10, save the sm for xiaofeng liang with different particle size
PRO RFD_orth_encoding_1d_SM_save
  
  N_particles = 5
  N_ele = 256;400.
  Particle_size_array = FINDGEN(N_particles)*10.+10.
  FOR k=0, N_particles-1 DO BEGIN
    particle_size = Particle_size_array[k]
    RFD_orth_encoding_1d_ellipsoidal, particle_size, N_ele, A
    filename_System = 'C:\D_drive\MPI_Tianjie\RotationalDriftMPI\xiaofeng_thesis\systemMatrix'+STRTRIM(STRING(FLOOR(particle_size)),1)+'.dat'
    ;write file
    openw,unit,filename_System,/get_lun
    FOR i=0, N_ele-1 DO BEGIN
      FOR j=0, N_ele-1 DO BEGIN
        ;; write the data points
        WRITEU, unit, A[i, j]
      ENDFOR
    ENDFOR
    ;; close file and END execution
    FREE_LUN, unit
    print, 'k=',k
  ENDFOR

END

;@GJ, 2022/12/30, do the weak gradient field based encoding
;@GJ, 2023/3/10, generate for Xiaofeng Liang
;A is the system matrix
PRO RFD_orth_encoding_1d_ellipsoidal, particle_size, N_ele, A
  B_0 = 10.0; mT
  G_z = 0.04; mT/mm
  IF N_ELEMENTS(N_ele) EQ 0 THEN N_ele = 400.
  FOV = 200. ;mm; 20 cm
  z = FINDGEN(N_ele)/N_ele*FOV - FOV/2.
  H_z = G_z * z + B_0
  alpha_array = H_z/B_0
  H_array = SQRT(H_z^2 + B_0^2)

  ;@GJ, 2022/12/30
  ;calculate M0, the net magnetization
  Msat_T = 0.551 ; T/u0
  u0 = 4.0 * !PI * (0.1^4) ; T*m/kA
  Msat_kAm = Msat_T/u0; kA/m
  IF N_ELEMENTS(particle_size) EQ 0 THEN particle_size = 20. ;nm
  m_moment = 1./6. * !PI * (particle_size^3) * Msat_kAm ; 10^-27 kA*m2
  k_B = 1.380469 ; 10^-23 J/K
  T_p = 36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = m_moment / (K_B * T_p_kelvin) / (10.^4) ; in kA*m2/J or 1/mT
  M_0_array = Msat_kAm * (1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
;  iplot, z, M_0_array, xtitle='FOV (mm)', ytitle='M0'
  
;  ;shepp_logan_phantom_2d
;  image = BYTSCL(shepp_logan_phantom_2d(N_ele))
;;  iimage, image, title='original'
;  image_1d = DBLARR(N_ele)*0.
;  FOR i=0, N_ele-1 DO BEGIN
;    image_1d[i] = MEAN(image[*,i])
;  ENDFOR
  image_1d=image_phantom_res(N_ele)
;  image_1d[0:188]=0
;  image_1d[218:*]=0
  
  M_0_1d_array = image_1d * M_0_array
   
  ;by assuming B/Bc = 0.5;0.85
  B_c = B_0/0.5
  w_0 = 2./1000.; kHz/1000 = MHz
  delta_t = 1.0; us
  N_times = 100. / w_0 / delta_t
;  N_times = 250. / w_0 / delta_t
  ;N_times = 500.
  time_array = FINDGEN(N_times) * delta_t
  t_ele = N_times
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
  else $
    freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)
  
  Mx_t_array = time_array*0.
  Mz_t_array = time_array*0.
  signal_array = time_array*0.
  ;@GJ, 2022/12/30, initiate the M_t_array
  Mxp_t_array = DBLARR(N_times, N_ele)*0.
  Mzp_t_array = DBLARR(N_times, N_ele)*0.
  SM_Mxp_t_array = DBLARR(N_times, N_ele)*0.
  SM_Mzp_t_array = DBLARR(N_times, N_ele)*0.
  ;t=0, initial net magnetization
  Mxp_t_array[0,*] = 0
  Mzp_t_array[0,*] = M_0_1d_array
  SM_Mxp_t_array[0,*] = 0
  SM_Mzp_t_array[0,*] = M_0_1d_array
  Mx_t_array[0] = TOTAL(Mxp_t_array[0,*])
  Mz_t_array[0] = TOTAL(Mzp_t_array[0,*])
  signal_t_array = DBLARR(N_times, N_ele)*0.
  theta_t_array = DBLARR(N_times, N_ele)*0.
  w_D_t_array = DBLARR(N_times, N_ele)*0.
  FOR i=0, N_times-2 DO BEGIN
    FOR j=0, N_ele-1 DO BEGIN
      theta_temp = theta_t_array[i, j]
      ;@GJ,2022/12/30, calculate the total B^2
      B_square = H_z[j]^2*cos(theta_temp)^2 + B_0^2*sin(theta_temp)^2
      IF B_square LT B_c^2 THEN BEGIN
        w_D = w_0 * (1. - SQRT(1. - B_square/(B_c^2)))
      ENDIF ELSE BEGIN
        w_D = w_0
      ENDELSE
      theta_temp += 2. * !PI * w_D * delta_t
      SM_Mzp_t_array[i+1,j] = M_0_array[j] * cos(theta_temp)
      SM_Mxp_t_array[i+1,j] = M_0_array[j] * sin(theta_temp)
      Mzp_t_array[i+1,j] = image_1d[j] * SM_Mzp_t_array[i+1,j]
      Mxp_t_array[i+1,j] = image_1d[j] * SM_Mxp_t_array[i+1,j]
      theta_t_array[i+1, j] = theta_temp
      w_D_t_array[i+1,j] = w_D
    ENDFOR
    Mx_t_array[i+1] = TOTAL(Mxp_t_array[i+1,*])
    Mz_t_array[i+1] = TOTAL(Mzp_t_array[i+1,*])
  ENDFOR
  
  ;based on ART fot reconstruction
  ;A_invert = INVERT(A, /DOUBLE)
  ;rec_image = A_invert ## b
  ;print, 'invert: ', rec_image
  A = SM_Mzp_t_array[0:(N_ele-1)*50.:50., 0:N_ele-1]
  b = Mz_t_array[0:(N_ele-1)*50.:50.]
  rec_image_1d = art_func(A,b)
  iplot, rec_image_1d
  print, 'ART: ', rec_image_1d
  
  ;calculate the signal
  signalz_t_array = Mz_t_array[1:N_times-1]-Mz_t_array[0:N_times-2]
  signalx_t_array = Mx_t_array[1:N_times-1]-Mx_t_array[0:N_times-2]
  iplot, time_array[0:N_times-2], signalz_t_array, xtitle='time [us]', ytitle='signal z direction'
  iplot, time_array[0:N_times-2], signalx_t_array, xtitle='time [us]', ytitle='signal x direction'
  
  ;fast FFT
  aaa = FFT(COMPLEX(TOTAL(signalz_t_array, /CUMULATIVE), TOTAL(signalx_t_array, /CUMULATIVE)))
  powerSpectrum = ABS(aaa)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  iplot, freq[0:N_ele-1], scaledPowerSpect[0:N_ele-1], title='reconstructed'
  iplot, z, M_0_1d_array, xtitle='FOV (mm)', ytitle='M0'

  ;@GJ, 2022/11/2, do another interpolation
;  Power_inter = INTERPOL(scaledPowerSpect[0:t_ele/2-1], B_Bc_int, B_Bc_array)
END

PRO RFD_gradient_encoding_resolution

  w0 = 2.; kHz
  N_period_acq = 250.;500;250.;150.
  delta_w = w0/(N_period_acq*2.)
  Bc =  10. ;mT
  gradient = 0.040 ;mT/mm
  
  ;B/Bc range: 0.3~0.7
  FFL_ele = 800
  B_Bc_range = [0.1, 0.9]
  B_Bc_array = FINDGEN(FFL_ele)/(FFL_ele-1)*(B_Bc_range[1]-B_Bc_range[0])+B_Bc_range[0]
  wD_array = w0 * (1. - SQRT(1-B_Bc_array^2))
  
  delta_x = delta_w / w0 * SQRT(1-B_Bc_array^2) / B_Bc_array * Bc / gradient
  
  iplot, B_Bc_array, delta_x, xtitle='B/Bc', ytitle='Resolution [mm]', title='Resolution vs B/Bc'

END

;@GJ, 2022/11/2, use gradient encoding and rotation drift to scan and reconstruct the image
;@GJ, 2022/11/2, the drift frequency is added
;@GJ, 2023/1/2, 1Dimage was reconstructed
;@GJ, 2023/1/2, 2D image was reconstructed
PRO RFD_gradient_encoding_1d
  ;excitation field amplitude 5 kHz
  w0 = 2.; kHz
  RFD_equation, w0, wD, B_Bc, N_ele
  
  ;B/Bc range: 0.3~0.7
  FFL_ele = 800
  B_Bc_range = [0.5, 0.9]
  B_Bc_array = FINDGEN(FFL_ele)/(FFL_ele-1)*(B_Bc_range[1]-B_Bc_range[0])+B_Bc_range[0]
  wD_array = w0 * (1. - SQRT(1-B_Bc_array^2))
  print, 'B range: ', MIN(B_Bc_array), ' ~ ', MAX(B_Bc_array)
  ;plot the array
  iplot, B_Bc_array, wD_array, title='RFD', xtitle='B/Bc', ytitle='wD [kHz]'
  
  ;get the sample image
  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), FFL_ele, FFL_ele))
  
  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image = BYTSCL(CONGRID(REFORM(phantom_image[0,*,*]), FFL_ele, FFL_ele))
  
  ;shepp_logan_phantom_2d
  image = BYTSCL(shepp_logan_phantom_2d(FFL_ele))
  
  FOV = 200. ; mm
  delta_x = FOV/FFL_ele ;mm
  x_array = FINDGEN(FFL_ele)/(FFL_ele-1)*FOV
  image_1d=image_phantom_res(FFL_ele)
  image_1d[0:188]=0
  image_1d[200:*]=0;image_1d[218:*]=0
  image_1d[230:235]=MAX(image_1d[0:200])
  image_1d=SHIFT(image_1d, 500)
  image_1d[*] = 0
  FOR i=0,FFL_ele-1 DO BEGIN
    IF i mod 50 EQ 0 THEN image_1d[i]=255
  ENDFOR
  
  image_1d[*] = 0
;  image_1d[764:765]=255
;  image_1d[768:769]=255
  image_1d[464:465]=255
  image_1d[468:469]=255
  image[*, 0] = REFORM(image_1d)
  ;iimage, image, title='original'
  IF N_ELEMENTS(image_1d) GT 10 THEN iplot, x_array, image_1d, xtitle='mm', title='original 1d image' ELSE iimage, image, title='original'
  
;  t_ele = 30.*FFL_ele
  delta_t = 5.0 ; us; samplying rate: 200kHz;1MHz
  N_tp_period = 1000./ w0 / delta_t
  N_period_acq = 250.;500;250.;150.
  t_ele = N_tp_period * N_period_acq
  t_array = FINDGEN(t_ele)*delta_t
  rec_array = image*0.
  
  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
  else $
    freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)
  
  ;@GJ, 2022/11/2, converting freq to B/Bc
  B_Bc_int = INTERPOL(B_Bc, wD, freq[0:t_ele/2-1])
  
  FOR i=0, FFL_ele-1 DO BEGIN
;    i=FFL_ele/2
;    i=234
    c_array = REFORM(image[*, i])
    ;IF i EQ 260 THEN BEGIN c_array = REFORM(image_1d)
    signal_array = DBLARR(t_ele)*0.
    IF max(c_array) EQ 0 THEN CONTINUE
    ;debug
;    c_array = c_array*0.
;    c_array[200] = 1.
;    IF i EQ 260 THEN iplot, c_array
    FOR t_ind = 0, t_ele-1 DO BEGIN
      FOR j=0, FFL_ele-1 DO BEGIN
        signal_t_j = c_array[j] * cos(2. * !PI * wD_array[j] * t_array[t_ind] / 1000.) * delta_x
;        print, 'angle = ', 2. * !PI * wD_array[j] * t_array[t_ind] / 1000. * 180.
        signal_array[t_ind] += signal_t_j
;        iplot, t_array, cos(2. * !PI * wD_array[122] * t_array / 1000.), xtitle='time [us]'
      ENDFOR
    ENDFOR
    
;    iplot, t_array, signal_array, xtitle='time [us]', ytitle='signal', title='1D signal in K-space'
    ;fast FFT
    aaa = FFT(signal_array, /inverse, /double)
    powerSpectrum = ABS(aaa)^2
    scaledPowerSpect = ALOG10(powerSpectrum)
    ;@GJ, 2022/11/2, do another interpolation
    Power_inter = INTERPOL(scaledPowerSpect[0:t_ele/2-1], B_Bc_int, B_Bc_array)
    IF i EQ 0 THEN BEGIN
      iplot, t_array/1000., signal_array, xtitle='time [ms]', ytitle='signal', title='1D signal in K-space ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
      iplot, freq[0:t_ele/2-1], scaledPowerSpect[0:t_ele/2-1], xtitle='frequency [kHz]', title='FFT ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
;      iplot, B_Bc_int, scaledPowerSpect[0:t_ele/2-1], xtitle='B/Bc', title='FFT'
      iplot, B_Bc_array, Power_inter, xtitle='B/Bc', title='Reconstructed 1D image ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
      IF N_ELEMENTS(image_1d) GT 10 THEN return
    ENDIF
    ;save the reconstructed impaage
    rec_array[*,i] = Power_inter
    print, 'i = ', i
  ENDFOR
  
  iimage, rec_array, title='reconstructed image'
  
  correlation_coe = CORRELATE(rec_array, image)
  print, 'correlation coefficient = ', correlation_coe
;  MSE_coe[i] = calMSE(REFORM(test_image[i,*,*]), REFORM(rec_image[i,*,*]))
;  RMSE_coe[i] = Sqrt(Total((REFORM(test_image[i,*,*]) - REFORM(rec_image[i,*,*]))^2)/N_Elements(REFORM(rec_image[i,*,*])))

END

;@GJ, 2022/11/2, use gradient encoding and rotation drift to scan and reconstruct the image
;@GJ, 2022/11/2, the drift frequency is added
;@GJ, 2023/1/2, 1Dimage was reconstructed
;@GJ, 2023/1/2, 2D image was reconstructed
PRO RFD_gradient_encoding_2d_FFL
  ;excitation field amplitude 5 kHz
  w0 = 2.; kHz
  RFD_equation, w0, wD, B_Bc, N_ele

  ;B/Bc range: 0.3~0.7
  FFL_ele = 800
  B_Bc_range = [0.5, 0.9]
  B_Bc_array = FINDGEN(FFL_ele)/(FFL_ele-1)*(B_Bc_range[1]-B_Bc_range[0])+B_Bc_range[0]
  wD_array = w0 * (1. - SQRT(1-B_Bc_array^2))
  print, 'B range: ', MIN(B_Bc_array), ' ~ ', MAX(B_Bc_array)
  ;plot the array
  iplot, B_Bc_array, wD_array, title='RFD', xtitle='B/Bc', ytitle='wD [kHz]'

  ;get the sample image
  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), FFL_ele, FFL_ele))

  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image = BYTSCL(CONGRID(REFORM(phantom_image[0,*,*]), FFL_ele, FFL_ele))

  ;shepp_logan_phantom_2d
  image = BYTSCL(shepp_logan_phantom_2d(FFL_ele))

  FOV = 200. ; mm
  delta_x = FOV/FFL_ele ;mm
  x_array = FINDGEN(FFL_ele)/(FFL_ele-1)*FOV
  iimage, image, title='original'

  ;  t_ele = 30.*FFL_ele
  delta_t = 5.0 ; us; samplying rate: 200kHz;1MHz
  N_tp_period = 1000./ w0 / delta_t
  N_period_acq = 250.;500;250.;150.
  t_ele = N_tp_period * N_period_acq
  t_array = FINDGEN(t_ele)*delta_t
  rec_array = image*0.

  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
  else $
    freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)

  ;@GJ, 2022/11/2, converting freq to B/Bc
  B_Bc_int = INTERPOL(B_Bc, wD, freq[0:t_ele/2-1])

  FOR i=0, FFL_ele-1 DO BEGIN
    ;    i=FFL_ele/2
    ;    i=234
    c_array = REFORM(image[*, i])
    ;IF i EQ 260 THEN BEGIN c_array = REFORM(image_1d)
    signal_array = DBLARR(t_ele)*0.
    IF max(c_array) EQ 0 THEN CONTINUE
    ;debug
    ;    c_array = c_array*0.
    ;    c_array[200] = 1.
    ;    IF i EQ 260 THEN iplot, c_array
    FOR t_ind = 0, t_ele-1 DO BEGIN
      FOR j=0, FFL_ele-1 DO BEGIN
        signal_t_j = c_array[j] * cos(2. * !PI * wD_array[j] * t_array[t_ind] / 1000.) * delta_x
        ;        print, 'angle = ', 2. * !PI * wD_array[j] * t_array[t_ind] / 1000. * 180.
        signal_array[t_ind] += signal_t_j
        ;        iplot, t_array, cos(2. * !PI * wD_array[122] * t_array / 1000.), xtitle='time [us]'
      ENDFOR
    ENDFOR

    ;    iplot, t_array, signal_array, xtitle='time [us]', ytitle='signal', title='1D signal in K-space'
    ;fast FFT
    aaa = FFT(signal_array, /inverse, /double)
    powerSpectrum = ABS(aaa)^2
    scaledPowerSpect = ALOG10(powerSpectrum)
    ;@GJ, 2022/11/2, do another interpolation
    Power_inter = INTERPOL(scaledPowerSpect[0:t_ele/2-1], B_Bc_int, B_Bc_array)
    IF i EQ 33 THEN BEGIN
      iplot, t_array/1000., signal_array, xtitle='time [ms]', ytitle='signal', title='1D signal in K-space ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
      iplot, freq[0:t_ele/2-1], scaledPowerSpect[0:t_ele/2-1], xtitle='frequency [kHz]', title='FFT ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
      ;      iplot, B_Bc_int, scaledPowerSpect[0:t_ele/2-1], xtitle='B/Bc', title='FFT'
      iplot, B_Bc_array, Power_inter, xtitle='B/Bc', title='Reconstructed 1D image ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
    ENDIF
    ;save the reconstructed impaage
    rec_array[*,i] = Power_inter
    print, 'i = ', i
  ENDFOR

  iimage, rec_array, title='reconstructed image'

  correlation_coe = CORRELATE(rec_array, image)
  print, 'correlation coefficient = ', correlation_coe
  ;  MSE_coe[i] = calMSE(REFORM(test_image[i,*,*]), REFORM(rec_image[i,*,*]))
  ;  RMSE_coe[i] = Sqrt(Total((REFORM(test_image[i,*,*]) - REFORM(rec_image[i,*,*]))^2)/N_Elements(REFORM(rec_image[i,*,*])))

END

;@GJ, 2022/11/2, use gradient encoding and rotation drift to scan and reconstruct the image
;@GJ, 2022/11/2, the drift frequency is added
;@GJ, 2023/1/2, 1Dimage was reconstructed
;@GJ, 2023/1/2, 2D image was reconstructed
;@GJ, 2023/1/2, 2D image using radon transform was reconstructed
PRO RFD_gradient_encoding_2d_Radon
  ;excitation field amplitude 5 kHz
  w0 = 2.; kHz
  RFD_equation, w0, wD, B_Bc, N_ele

  ;B/Bc range: 0.3~0.7
  FFL_ele = 400
  
  ;get the sample image
  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), FFL_ele, FFL_ele))
  
  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image = BYTSCL(CONGRID(REFORM(phantom_image[0,*,*]), FFL_ele, FFL_ele))

  ;shepp_logan_phantom_2d
  image = BYTSCL(shepp_logan_phantom_2d(FFL_ele))
  
  ;vessel image
  p_fn = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20230102c\MIP_Unet\002.jpeg'
  READ_JPEG, p_fn, phantom_image
  image = BYTSCL(CONGRID(phantom_image, FFL_ele, FFL_ele))
  
  FOV = 200. ; mm
  delta_x = FOV/FFL_ele ;mm
  x_array = FINDGEN(FFL_ele)/(FFL_ele-1)*FOV
  iimage, image, title='original'
  ;@GJ, 2023/1/2, do the Radon transform
  Radon_result_ori = RADON(image, RHO=rho, THETA=theta)
  Radon_result_temp = Radon_result_ori * 0.
  iimage, Radon_result_ori, title='original sinogram'
  FOR ij = 0, N_ELEMENTS(theta)-1 DO Radon_result_temp[ij, *] = Filter_FFT(REFORM(Radon_result_ori[ij, *]))
  backproject_rec_dire = RADON(Radon_result_temp, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject_rec_dire, title='directly rec'
  
  ;@GJ, 2022/1/2, set up the B and Bc based on the Radon sinogram
  B_Bc_range = [0.5, 0.9]
  B_Bc_array = FINDGEN(N_ELEMENTS(rho))/(N_ELEMENTS(rho)-1)*(B_Bc_range[1]-B_Bc_range[0])+B_Bc_range[0]
  wD_array = w0 * (1. - SQRT(1-B_Bc_array^2))
  print, 'B range: ', MIN(B_Bc_array), ' ~ ', MAX(B_Bc_array)
  ;plot the array
  iplot, B_Bc_array, wD_array, title='RFD', xtitle='B/Bc', ytitle='wD [kHz]'
  
  ;  t_ele = 30.*FFL_ele
  delta_t = 5.0 ; us; samplying rate: 200kHz;1MHz
  N_tp_period = 1000./ w0 / delta_t
  N_period_acq = 250.;500;250.;150.
  t_ele = N_tp_period * N_period_acq
  t_array = FINDGEN(t_ele)*delta_t
  rec_array = Radon_result_ori*0.
  rec_array_filter = Radon_result_ori*0.

  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
  else $
    freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)

  ;@GJ, 2022/11/2, converting freq to B/Bc
  B_Bc_int = INTERPOL(B_Bc, wD, freq[0:t_ele/2-1])

  FOR i=0, N_ELEMENTS(theta)-1 DO BEGIN
    ;    i=FFL_ele/2
    ;    i=234
    c_array = REFORM(Radon_result_ori[i, *])
    ;IF i EQ 260 THEN BEGIN c_array = REFORM(image_1d)
    signal_array = DBLARR(t_ele)*0.
    IF max(c_array) EQ 0 THEN CONTINUE
    ;debug
    ;    c_array = c_array*0.
    ;    c_array[200] = 1.
    ;    IF i EQ 260 THEN iplot, c_array
    FOR t_ind = 0, t_ele-1 DO BEGIN
      FOR j=0, N_ELEMENTS(rho)-1 DO BEGIN
        signal_t_j = c_array[j] * cos(2. * !PI * wD_array[j] * t_array[t_ind] / 1000.) * delta_x
        ;        print, 'angle = ', 2. * !PI * wD_array[j] * t_array[t_ind] / 1000. * 180.
        signal_array[t_ind] += signal_t_j
        ;        iplot, t_array, cos(2. * !PI * wD_array[122] * t_array / 1000.), xtitle='time [us]'
      ENDFOR
    ENDFOR

    ;    iplot, t_array, signal_array, xtitle='time [us]', ytitle='signal', title='1D signal in K-space'
    ;fast FFT
    aaa = FFT(signal_array, /inverse, /double)
    powerSpectrum = ABS(aaa)^2
    scaledPowerSpect = ALOG10(powerSpectrum)
    ;@GJ, 2022/11/2, do another interpolation
    Power_inter = INTERPOL(scaledPowerSpect[0:t_ele/2-1], B_Bc_int, B_Bc_array)
    IF i EQ 33 THEN BEGIN
      iplot, c_array, xtitle='rho', title='1D sinogram'
      iplot, t_array/1000., signal_array, xtitle='time [ms]', ytitle='signal', title='1D signal in K-space ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
      iplot, freq[0:t_ele/2-1], scaledPowerSpect[0:t_ele/2-1], xtitle='frequency [kHz]', title='FFT ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
      ;      iplot, B_Bc_int, scaledPowerSpect[0:t_ele/2-1], xtitle='B/Bc', title='FFT'
      iplot, B_Bc_array, Power_inter, xtitle='B/Bc', title='Reconstructed 1D sinogram ('+STRTRIM(FLOOR(N_period_acq), 2)+' periods)'
    ENDIF
    ;save the reconstructed impaage
    rec_array[i, *] = REFORM(Power_inter)
    rec_array_filter[i, *] = Filter_FFT(REFORM(rec_array[i, *]))
    print, 'i = ', i, '/', N_ELEMENTS(theta)-1
  ENDFOR
  
  backproject_rec = RADON(rec_array_filter, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, rec_array, title='reconstructed sinogram'
  iimage, backproject_rec, title='reconstructed image'

  correlation_coe = CORRELATE(backproject_rec, image)
  print, 'correlation coefficient = ', correlation_coe
  ;  MSE_coe[i] = calMSE(REFORM(test_image[i,*,*]), REFORM(rec_image[i,*,*]))
  ;  RMSE_coe[i] = Sqrt(Total((REFORM(test_image[i,*,*]) - REFORM(rec_image[i,*,*]))^2)/N_Elements(REFORM(rec_image[i,*,*])))

END


PRO RFD_gradient_encoding_echoes_delta_t_array
  
  n_delta_t = 6.
  delta_t_array = FINDGEN(n_delta_t)+1. ; us
  correlation_coe_array = delta_t_array * 0.
  
  N_echo = 1
  FOR i=0, n_delta_t-1 DO BEGIN
    RFD_gradient_encoding_echoes, N_echo, delta_t_array[i], correlation_coe_temp
    correlation_coe_array[i] = correlation_coe_temp
    print, 'Sampling rate: ', 1./delta_t_array[i], ' MHz'
  ENDFOR
  
  iplot, REVERSE(1./delta_t_array), REVERSE(correlation_coe_array), xtitle='Sampling Rate [MHz]', ytitle='Coefficient', title='Correlation Coeff'
END

;@GJ, 2022/11/3, generate the echo signals and K-space image
PRO RFD_gradient_encoding_echoes, N_echo, delta_t, correlation_coe, noise_level, TE_ms
  ;excitation field amplitude 5 kHz
  w0 = 30.; kHz
  RFD_equation, w0, wD, B_Bc, N_ele

  ;B/Bc range: 0.3~0.7
  FFL_ele = 400;400;100;400
  wD_array = wD[N_ele/2-FFL_ele/2+1:N_ele/2+FFL_ele/2]
  B_Bc_array = B_Bc[N_ele/2-FFL_ele/2+1:N_ele/2+FFL_ele/2]
  ;plot the array
  iplot, B_Bc_array, wD_array, title='RFD', xtitle='B/Bc', ytitle='wD [kHz]'

  ;get the sample image
  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), FFL_ele, FFL_ele))

  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image = BYTSCL(CONGRID(REFORM(phantom_image[0,*,*]), FFL_ele, FFL_ele))
  
  ;shepp_logan_phantom_2d
  ;image = BYTSCL(shepp_logan_phantom_2d(FFL_ele))
  iimage, image, title='original'
  
  IF N_ELEMENTS(delta_t) EQ 0 THEN delta_t = 5.0 ; us; samplying rate: 1MHz
  ;@GJ, set the TE
  IF N_ELEMENTS(TE_ms) EQ 0 THEN BEGIN
    t_ele = N_ele ; number of t_ele in FID
    TE_ms = 2 * t_ele * delta_t / 1000. ; in ms
  ENDIF ELSE BEGIN
    IF TE_ms LT 5*delta_t THEN BEGIN
      t_ele = N_ele ; number of t_ele in FID
      TE_ms = 2 * t_ele * delta_t / 1000. ; in ms
    ENDIF ELSE BEGIN
      t_ele = TE_ms * 1000. / 2. / delta_t ;number of t_ele in FID
    ENDELSE
  ENDELSE
  
  print, 'TE = ', TE_ms, ' ms'
  
  ;@set the N_echo
  IF N_ELEMENTS(N_echo) EQ 0 THEN BEGIN
    N_echo = 1 ; at least one echo
  ENDIF ELSE BEGIN
    IF N_echo LT 0 THEN BEGIN
      N_echo = 1 ; at least one echo
    ENDIF
  ENDELSE
  
  IF N_ELEMENTS(noise_level) EQ 0 THEN BEGIN
    noise_level = 10.
  ENDIF ELSE BEGIN
    IF noise_level LE 0 THEN BEGIN
      noise_level = 10 ; at least 10
    ENDIF
  ENDELSE
  
;  N_echo=0
  signal_array = DBLARR(t_ele+2*t_ele*N_echo)*0.
  FOV = 200. ; mm
  delta_x = FOV/FFL_ele ;mm
  t_array = FINDGEN(t_ele+2*t_ele*N_echo)*delta_t
  rec_array = image*0.
  k_space_array = image*0.
  phase_array = image*0.

  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
  else $
    freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)

  ;@GJ, 2022/11/2, converting freq to B/Bc
  B_Bc_int = INTERPOL(B_Bc, wD, freq[0:t_ele/2-1])

  FOR i=0, FFL_ele-1 DO BEGIN
;  FOR i=FFL_ele/4, FFL_ele/4 DO BEGIN
    ;    i=FFL_ele/2
    ;    i=234
    c_array = REFORM(image[*, i])
    signal_array *= 0.
    IF max(c_array) EQ 0 THEN CONTINUE
    ;debug
    ;    c_array = c_array*0.
    ;    c_array[200] = 1.
    ;    IF i EQ 260 THEN iplot, c_array
    FOR t_ind = 0, t_ele-1 DO BEGIN
      FOR j=0, FFL_ele-1 DO BEGIN
        angle = 2. * !PI * wD_array[j] * t_array[t_ind] / 1000.
        seed = !NULL
        noise_angle = (RANDOMU(seed)-0.5)/90.*!PI*noise_level
        angle += noise_angle
        signal_t_j = c_array[j] * cos(angle) * delta_x
        ;        print, 'angle = ', 2. * !PI * wD_array[j] * t_array[t_ind] / 1000. * 180.
        signal_array[t_ind] += signal_t_j
        ;        iplot, t_array, cos(2. * !PI * wD_array[122] * t_array / 1000.), xtitle='time [us]'
        ;@GJ, record the last phase angle
        phase_array[j] = angle
      ENDFOR
    ENDFOR
    
    ;@GJ, generate the echo signals
    IF N_echo GE 1 THEN BEGIN
      FOR k=0, N_echo-1 DO BEGIN
        FOR t_ind = 0, 2*t_ele-1 DO BEGIN
          ;set the signal calculation
          FOR j=0, FFL_ele-1 DO BEGIN
            ;if k is even, use minus
            IF (k mod 2) EQ 0 THEN BEGIN
              angle = phase_array[j] - 2. * !PI * wD_array[j] * t_array[t_ind] / 1000.
            ENDIF ELSE BEGIN
              ;if k is odd, use plus
              angle = phase_array[j] + 2. * !PI * wD_array[j] * t_array[t_ind] / 1000.
            ENDELSE
            seed = !NULL
            noise_angle = (RANDOMU(seed)-0.5)/90.*!PI*noise_level
            angle += noise_angle
            signal_t_j = c_array[j] * cos(angle) * delta_x
            signal_array[(2*k+1)*t_ele + t_ind] += signal_t_j
            ;save the phase angle at the end of the TE cycle
            IF t_ind EQ 2*t_ele-2 THEN phase_array[j] = angle
          ENDFOR
        ENDFOR
        IF k EQ 0 THEN k_space_array[*,i] = CONGRID(signal_array[t_ele:3*t_ele-1], FFL_ele)
      ENDFOR
    ENDIF
    
    ;    iplot, t_array, signal_array, xtitle='time [us]', ytitle='signal', title='1D signal in K-space'
    ;fast FFT
    aaa = FFT(signal_array[0:t_ele-1])
    powerSpectrum = ABS(aaa)^2
    scaledPowerSpect = ALOG10(powerSpectrum)
    ;@GJ, 2022/11/2, do another interpolation
    Power_inter = INTERPOL(scaledPowerSpect[0:t_ele/2-1], B_Bc_int, B_Bc_array)
    IF i EQ FFL_ele/2 THEN BEGIN
      iplot, t_array/1000., signal_array, xtitle='time [ms]', xrange=[-1, MAX(t_array/1000.)+1], ytitle='signal', title='1D signal in K-space'
      iplot, t_array[0:t_ele-1]/1000., signal_array[0:t_ele-1], xtitle='time [ms]', xrange=[-1, MAX(t_array[0:t_ele-1]/1000.)+1], ytitle='signal', title='1D dephasing signal'
      iplot, freq[0:t_ele/2-1], scaledPowerSpect[0:t_ele/2-1], xtitle='frequency [kHz]', title='FFT'
      ;      iplot, B_Bc_int, scaledPowerSpect[0:t_ele/2-1], xtitle='B/Bc', title='FFT'
      iplot, B_Bc_array, Power_inter, xtitle='B/Bc', title='FFT in B/Bc'
    ENDIF
    ;save the reconstructed impaage
    rec_array[*,i] = Power_inter
    IF i mod 50 EQ 0 THEN print, 'i = ', i
  ENDFOR
  
  IF N_echo GE 1 THEN iimage, k_space_array, title='K space'
  iimage, rec_array, title='reconstructed image'
  correlation_coe = CORRELATE(rec_array, image)
  print, 'sampling rate: ', delta_t, 'us'
  print, 'correlation coefficient = ', correlation_coe

END

;Author:Guang Jia
;2020-05-06

PRO FFT_SENSE

  ;读取图像
  file =  DIALOG_PICKFILE(/READ, FILTER = '*.bmp')
  binary_img = CONGRID(READ_BMP(file), 256, 256)
;  img01 = image(binary_img,LAYOUT=[2,2,1],title='Original')
  iimage, binary_img, title='Original'
  

  ;傅里叶变换
  N_col=N_elements(binary_img[*,0])
  M_row=N_elements(binary_img[0,*])
  M=M_row
  N=N_col
  
  ;fast FFT
  aaa = FFT(binary_img, /center)
  powerSpectrum = ABS(aaa)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  iimage, scaledPowerSpect, title='FFT'
  bbb = FFT(aaa, /inverse)
  powerSpectrum_inv = ABS(bbb)^2
;  scaledPowerSpect_inv = ALOG10(powerSpectrum_inv)
  iimage, powerSpectrum_inv, title='inverse FFT'
  
  aaa_SENSE = COMPLEXARR(N_col/2, M_row)
  FOR i=0, N_col/2-1 DO BEGIN
    aaa_SENSE[i, *] = aaa[i*2, *]
  ENDFOR
  powerSpectrum_SENSE = ABS(aaa_SENSE)^2
  scaledPowerSpect_SENSE = ALOG10(powerSpectrum_SENSE)
  iimage, scaledPowerSpect_SENSE, title='FFT SENSE'
  bbb_SENSE = FFT(aaa_SENSE, /inverse)
  powerSpectrum_inv_SENSE = ABS(bbb_SENSE)^2
  iimage, REAL_PART(powerSpectrum_inv_SENSE), title='inverse FFT SENSE'
  
  n_dft=complexarr(N,M)
  n_dft_tp=complexarr(N,M)
  dft_img=complexarr(N,M)
  ;列（一维DFT）
  FOR v=0,N-1,1 DO BEGIN
    FOR x=0,M-1,1 DO BEGIN
      fxv=complex(0.0d,0.0d)
      FOR y=0,N-1,1 DO BEGIN
        fxv=fxv+binary_img[y,x]*complex(cos(2*!pi*v*y/N),-1*sin(2*!pi*v*y/N))
      ENDFOR
      n_dft_tp[v,x]=fxv
    ENDFOR
  ENDFOR
  ;行（一维DFT）
  FOR u=0,M-1,1 DO BEGIN
    FOR y=0,N-1,1 DO BEGIN
      fxv=complex(0.0d,0.0d)
      FOR x=0,M-1,1 DO BEGIN
        fxv=fxv+n_dft_tp[y,x]*complex(cos(2*!pi*u*x/M),-1*sin(2*!pi*u*x/M))
      ENDFOR
      n_dft[y,u]=fxv
    ENDFOR
  ENDFOR
  n_dft=n_dft/(M*N)

  powerSpectrum = ABS(n_dft)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  img02 = IMAGE(scaledPowerSpect,/CURRENT,LAYOUT=[2,2,2],title='MyFFT')

  ;自带FFT
  ffTransform = FFT(binary_img, /center)
  powerSpectrum = ABS(ffTransform)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  img03 = IMAGE(scaledPowerSpect,/CURRENT,LAYOUT=[2,2,3],title='FFT')



END

PRO DFTSELF

  ;读取图像
  file =  DIALOG_PICKFILE(/READ, FILTER = '*.bmp')
  binary_img = CONGRID(READ_BMP(file), 256, 256)
  ;  img01 = image(binary_img,LAYOUT=[2,2,1],title='Original')
  iimage, binary_img, title='Original'


  ;傅里叶变换
  N_col=N_elements(binary_img[*,0])
  M_row=N_elements(binary_img[0,*])
  M=M_row
  N=N_col

  ;fast FFT
  aaa = FFT(binary_img, /center)
  powerSpectrum = ABS(aaa)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  iimage, scaledPowerSpect, title='FFT'
  bbb = FFT(aaa, /inverse)
  powerSpectrum_inv = ABS(bbb)^2
  ;  scaledPowerSpect_inv = ALOG10(powerSpectrum_inv)
  iimage, powerSpectrum_inv, title='inverse FFT'

  aaa_SENSE = COMPLEXARR(N_col/2, M_row)
  FOR i=0, N_col/2-1 DO BEGIN
    aaa_SENSE[i, *] = aaa[i*2, *]
  ENDFOR
  powerSpectrum_SENSE = ABS(aaa_SENSE)^2
  scaledPowerSpect_SENSE = ALOG10(powerSpectrum_SENSE)
  iimage, scaledPowerSpect_SENSE, title='FFT SENSE'
  bbb_SENSE = FFT(aaa_SENSE, /inverse)
  powerSpectrum_inv_SENSE = ABS(bbb_SENSE)^2
  iimage, REAL_PART(powerSpectrum_inv_SENSE), title='inverse FFT SENSE'

  n_dft=complexarr(N,M)
  n_dft_tp=complexarr(N,M)
  dft_img=complexarr(N,M)
  ;列（一维DFT）
  FOR v=0,N-1,1 DO BEGIN
    FOR x=0,M-1,1 DO BEGIN
      fxv=complex(0.0d,0.0d)
      FOR y=0,N-1,1 DO BEGIN
        fxv=fxv+binary_img[y,x]*complex(cos(2*!pi*v*y/N),-1*sin(2*!pi*v*y/N))
      ENDFOR
      n_dft_tp[v,x]=fxv
    ENDFOR
  ENDFOR
  ;行（一维DFT）
  FOR u=0,M-1,1 DO BEGIN
    FOR y=0,N-1,1 DO BEGIN
      fxv=complex(0.0d,0.0d)
      FOR x=0,M-1,1 DO BEGIN
        fxv=fxv+n_dft_tp[y,x]*complex(cos(2*!pi*u*x/M),-1*sin(2*!pi*u*x/M))
      ENDFOR
      n_dft[y,u]=fxv
    ENDFOR
  ENDFOR
  n_dft=n_dft/(M*N)

  powerSpectrum = ABS(n_dft)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  img02 = IMAGE(scaledPowerSpect,/CURRENT,LAYOUT=[2,2,2],title='MyFFT')

  ;自带FFT
  ffTransform = FFT(binary_img, /center)
  powerSpectrum = ABS(ffTransform)^2
  scaledPowerSpect = ALOG10(powerSpectrum)
  img03 = IMAGE(scaledPowerSpect,/CURRENT,LAYOUT=[2,2,3],title='FFT')

  ;反傅立叶变换
  ;列（一维DFT）
  FOR v=0,N-1,1 DO BEGIN
    FOR x=0,M-1,1 DO BEGIN
      fxv=complex(0.0d,0.0d)
      FOR y=0,N-1,1 DO BEGIN
        fxv=fxv+n_dft[y,x]*complex(cos(2*!pi*v*y/N),sin(2*!pi*v*y/N))
      ENDFOR
      n_dft_tp[v,x]=fxv
    ENDFOR
  ENDFOR
  ;行（一维DFT）
  FOR u=0,M-1,1 DO BEGIN
    FOR y=0,N-1,1 DO BEGIN
      fxv=complex(0.0d,0.0d)
      FOR x=0,M-1,1 DO BEGIN
        fxv=fxv+n_dft_tp[y,x]*complex(cos(2*!pi*u*x/M),sin(2*!pi*u*x/M))
      ENDFOR
      dft_img[y,u]=fxv
    ENDFOR
  ENDFOR
  ;dft_img=dft_img/(M*N)

  img04 = IMAGE(dft_img,/CURRENT,LAYOUT=[2,2,4],title='DFT to IMAGE')

END

;GJ, 2022/5/1, the calculation of AUC_flat
PRO AUC_flat
  
  H = (FINDGEN(401))/20.
  particle_size = 30.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  M = (1./TANH(beta_particle*H) - 1./(beta_particle*H))
;  M[200] = 0.
  iplot, H, M
  
  M_prime = H[0:399]*(M[1:400]-M[0:399])
  iplot, H[0:399], M_prime
  
  
  window, 1
  plot, H[0:399], M[0:399]
  oplot, H[0:399], 30*M_prime
  
  iplot, H[0:399], M[0:399] - M_prime
  
END

;FFT and magnetization curve
PRO SingalResponse

  X = 2*!PI/100 * FINDGEN(600)
  H = -COS(X)*20.
  M = 1./TANH(H) - 1./H
  iplot, H
  iplot, M
  
  plot, FFT(M)
  iplot, ABS(FFT(M[1:*])), xrange=[0,100]
  

END


PRO SingalResponse_7point
  ;!P.COLOR = '000000'x
  
  X = 2*!PI/100 * FINDGEN(600)
  H = -COS(X)*30.
  M = 1./TANH(H) - 1./H
;  iplot, H
;  iplot, M

;  plot, X, H, LINESTYLE=1
  plot, X, M, BACKGROUND = 'FFFFFF'x, COLOR = 0
  
  H2=H+10.
  M2=1./TANH(H2) - 1./H2
  oplot, X, M2, LINESTYLE=5, COLOR = '0000FF'x;, COLOR='red'
;  plot, FFT(M)
;  iplot, ABS(FFT(M[1:*])), xrange=[0,100]

  H3=H+20.
  M3=1./TANH(H3) - 1./H3
  oplot, X, M3, LINESTYLE=2, COLOR='FF0000'x
;  plot, FFT(M)
  
  ;plot signal
  window, 1
  Nelem = N_ELEMENTS(M)
  plot, X[1:*], M[1:*]-M[0:Nelem-1], BACKGROUND = 'FFFFFF'x, COLOR = 0
  oplot, X[1:*], M2[1:*]-M2[0:Nelem-1], LINESTYLE=5, COLOR = '0000FF'x
  oplot, X[1:*], M3[1:*]-M3[0:Nelem-1], LINESTYLE=2, COLOR='FF0000'x
  
  
  ;plot FFT
  window, 2
  plot, ABS(FFT(M[1:*])), xrange=[0,100], BACKGROUND = 'FFFFFF'x, COLOR = 0
  oplot, ABS(FFT(M2[1:*])), LINESTYLE=5, COLOR = '0000FF'x 
  oplot, ABS(FFT(M3[1:*])),  LINESTYLE=2, COLOR='FF0000'x
  
  
  window, 3
  plot, X, M, BACKGROUND = 'FFFFFF'x, COLOR = 0
  
  H_m2=H-10.
  M_m2=1./TANH(H_m2) - 1./H_m2
  oplot, X, M_m2, LINESTYLE=5, COLOR = '0000FF'x;, COLOR='red'


  H_m3=H-20.
  M_m3=1./TANH(H_m3) - 1./H_m3
  oplot, X, M_m3, LINESTYLE=2, COLOR='FF0000'x
  ;plot, FFT(M)
  
  ;plot signal
  window, 4
  Nelem = N_ELEMENTS(M)
  plot, X[1:*], M[1:*]-M[0:Nelem-1], BACKGROUND = 'FFFFFF'x, COLOR = 0
  oplot, X[1:*], M_m2[1:*]-M_m2[0:Nelem-1], LINESTYLE=5, COLOR = '0000FF'x
  oplot, X[1:*], M_m3[1:*]-M_m3[0:Nelem-1], LINESTYLE=2, COLOR='FF0000'x


  ;plot FFT
  window, 5
  plot, ABS(FFT(M[1:*])), xrange=[0,100], BACKGROUND = 'FFFFFF'x, COLOR = 0
  oplot, ABS(FFT(M_m2[1:*])), LINESTYLE=5, COLOR = '0000FF'x
  oplot, ABS(FFT(M_m3[1:*])),  LINESTYLE=2, COLOR='FF0000'x
END

;; write a dat file
;; major modifications in writing process
;; compared to Ulf
;;GJ, 2021/07/28, write signal as 1d dat file
pro write_MPI_signal_dat_file
  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20211109b\MPI_signal_1d.dat'
  Ndata=512
  X = 2.*!PI/Ndata * FINDGEN(Ndata+1)
  T = 1./16.5; us 10-6 sec
  H = -COS(X)*30.
  M = 1.0*(1./TANH(H) - 1./H)
  signal = BYTSCL(M[1:*]-M[0:Ndata-1])
  plot, signal
  
  ;write file
  openw,unit,filename,/get_lun
    ;; write the data points
  WRITEU, unit, signal
  
  ;; close file and END execution
  FREE_LUN, unit
END

pro read_MPI_signal_dat_file

  filename = dialog_pickfile(TITLE='Select MPI signal dat File', FILTER='*.dat', /MUST_EXIST)
  OPENR, unit, filename, /GET_LUN
  newdata = BYTARR(512)
  READU, unit, newdata
  FREE_LUN, unit
END

;GJ, 2021/7/2
;GJ, 2021/7/23, H is changed to -cos function
;signal vs Magnetic field strength
PRO SingalResponse_7point_different_driveAmp
  ;!P.COLOR = '000000'x
  
  ;pre-defined frequency: 33 KHz
  ;sampling frequency 16.5 MHz
  ;sampling time interval 0.06 us (10-6 sec)
  X = 2*!PI/100. * FINDGEN(600);/500
  T = 1./16.5; us 10-6 sec
  ;FFT frequency preparation; GJ, 2021/7/2
  N = 500.
  XN = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  ;freq unit is 0.033 MHz = 33 KHz
  if (is_N_even) then $
    freq = [0.0, XN, N/2, -N/2 + XN]/(N*T) $
  else $
    freq = [0.0, XN, -(N/2 + 1) + XN]/(N*T)
  
  H_array = [10., 30., 90., 3.33, 1.]
  threeH_array = H_array*0.
  
  H = -COS(X)*10.
  M = 1.0*(1./TANH(H) - 1./H)
;  plot, H, M, title=''
  ;  iplot, H
  ;  iplot, M
  
  window, 0
  ;  plot, X, H, LINESTYLE=1
  plot, X, M, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[-1.2, 1.2], title='Magnetization'


  H2=H*3
  M2=1./TANH(H2) - 1./H2
  oplot, X, M2, LINESTYLE=0, COLOR = '0000FF'x;, COLOR='red'
  ;  plot, FFT(M)
  ;  iplot, ABS(FFT(M[1:*])), xrange=[0,100]

  H3=H*9
  M3=1./TANH(H3) - 1./H3
  oplot, X, M3, LINESTYLE=0, COLOR='FF0000'x
  ;  plot, FFT(M)

  ;plot signal
  window, 1
  Nelem = N_ELEMENTS(M)
  plot, X[1:*], M[1:*]-M[0:Nelem-1], COLOR=0,BACKGROUND = 'FFFFFF'x, YRANGE=[-0.4, 0.4], title='Signal'
  oplot, X[1:*], M2[1:*]-M2[0:Nelem-1], LINESTYLE=0, COLOR = '0000FF'x
  oplot, X[1:*], M3[1:*]-M3[0:Nelem-1], LINESTYLE=0, COLOR = 'FF0000'x
  

  ;plot FFT
  window, 2
  ;freq is 33 KHz
  plot, (ABS(FFT(M[1:*])))[0:20], BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Spectrum' ;freq[0:8], 
  oplot, (ABS(FFT(M2[1:*])))[0:20], LINESTYLE=0, COLOR = '0000FF'x
  oplot, (ABS(FFT(M3[1:*])))[0:20],  LINESTYLE=0, COLOR='FF0000'x
  threeH_array[0] = (ABS(FFT(M[1:*])))[18]
  threeH_array[1] = (ABS(FFT(M2[1:*])))[18]
  threeH_array[2] = (ABS(FFT(M3[1:*])))[18]

  window, 3
  plot, X, M, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Magnetization'

  H_m2=H/3.
  M_m2=1./TANH(H_m2) - 1./H_m2
  oplot, X, M_m2, LINESTYLE=0, COLOR = '0000FF'x;, COLOR='red'


  H_m3=H/10.
  M_m3=1./TANH(H_m3) - 1./H_m3
  oplot, X, M_m3, LINESTYLE=0, COLOR='FF0000'x
  ;plot, FFT(M)

  ;plot signal
  window, 4
  Nelem = N_ELEMENTS(M)
  plot, X[1:*], M[1:*]-M[0:Nelem-1], BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Signal'
  oplot, X[1:*], M_m2[1:*]-M_m2[0:Nelem-1], LINESTYLE=0, COLOR = '0000FF'x
  oplot, X[1:*], M_m3[1:*]-M_m3[0:Nelem-1], LINESTYLE=0, COLOR='FF0000'x


  ;plot FFT
  window, 5
  plot, (ABS(FFT(M[1:*])))[0:20], BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Spectrum'
  oplot, (ABS(FFT(M_m2[1:*])))[0:20], LINESTYLE=0, COLOR = '0000FF'x
  oplot, (ABS(FFT(M_m3[1:*])))[0:20],  LINESTYLE=0, COLOR='FF0000'x
  
  threeH_array[3] = (ABS(FFT(M_m2[1:*])))[18]
  threeH_array[4] = (ABS(FFT(M_m3[1:*])))[18]
  
  window, 6
  plot, H_array, threeH_array, PSYM=4, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='3 times Harmonic Component'
  
END

;GJ, 2021/7/2
;GJ, 2021/7/23, H is changed to -cos function
;GJ, 2021/8/19, Harmonics versus excitation field
;signal vs Magnetic field strength
PRO Harmonics_vs_H
  ;!P.COLOR = '000000'x

  ;pre-defined frequency: 3.3 KHz
  ;sampling frequency 16.5 MHz
  ;sampling time interval 0.06 us (10-6 sec)
  X = 2*!PI/1000. * FINDGEN(6001)
  T = 1./16.5; us 10-6 sec
  time = FINDGEN(N_ELEMENTS(X)-1)*T
  ;FFT frequency preparation; GJ, 2021/7/2
  N = N_ELEMENTS(X)-1;500.
  XN = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  ;freq unit is 0.033 MHz = 33 KHz
  if (is_N_even) then $
    freq = [0.0, XN, N/2, -N/2 + XN]/(N*T) $
  else $
    freq = [0.0, XN, -(N/2 + 1) + XN]/(N*T)

  HN = 100
  H_array = (FINDGEN(HN)+1.)/HN*30.;from 0 to 99 [10., 30., 90., 3.33, 1.]
  firstHar_array = H_array*0.
  thirdHar_array = H_array*0.
  fifthHar_array = H_array*0.
  seventhHar_array = H_array*0.
  PeakAmp_array = H_array*0.
  fortynineHar_array = H_array*0.

  beta_Resovist = 2.08 ;for beta*H (mT)
  FOR i=0, HN-1 DO BEGIN
    H = -COS(X)*H_array[i]
    M = 1.0*(1./TANH(beta_Resovist* H) - 1./(beta_Resovist* H))
    
    Nelem = N_ELEMENTS(M)
    Signal = M[1:*]-M[0:Nelem-2]
;    IF i LT 20 THEN BEGIN
;      M = SMOOTH(M, 3, /NAN) ;to avoid the unsmoothed points
;      Signal = SMOOTH(signal, 3, /NAN) ;to avoid the unsmoothed points
;    ENDIF
    Signal_ext=Signal;[Signal, Signal, Signal, Signal, Signal, Signal];, Signal, Signal, Signal, Signal]
    fft_comp = ABS(FFT(Signal_ext))
    firstHar_array[i] = fft_comp[6]
    thirdHar_array[i] = fft_comp[18];/fft_comp[6]
    fifthHar_array[i] = fft_comp[30]
    seventhHar_array[i] = fft_comp[42]
    fortynineHar_array[i] = fft_comp[294]
    PeakAmp_array[i] = MAX(Signal)

    IF i EQ 0 THEN BEGIN
      window, 0
;      plot, time, Signal, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[-0.4, 0.4], title='Signal'
      plot, (freq*N*T)[0:200], fft_comp[0:200], BACKGROUND = 'FFFFFF'x, COLOR = i*500, YRANGE=[-0.0005, 0.004], title='Fourier Coefficients';XRANGE=[0,70], 
    ENDIF ELSE BEGIN
      oplot, (freq*N*T)[0:200], fft_comp[0:200], LINESTYLE=0, COLOR = i*500
;      oplot, time, Signal, LINESTYLE=0, COLOR = i*500
    ENDELSE
    wait, 0.1
  ENDFOR
  
  window, 1
  plot, H_array, firstHar_array, PSYM=4, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Harmonics'
  oplot, H_array, thirdHar_array, PSYM=2, COLOR = 500
  oplot, H_array, fifthHar_array, PSYM=1, COLOR = 1000
  oplot, H_array, seventhHar_array, PSYM=5, COLOR = 1500
  oplot, H_array, fortynineHar_array, PSYM=7, COLOR = 2000
   ; Compute the second degree polynomial fit to the data:
   result = POLY_FIT(H_array, thirdHar_array, 5, MEASURE_ERRORS=measure_errors, SIGMA=sigma, YFIT=thirdHar_y)
    ; Print the coefficients:
    PRINT, 'Coefficients: ', result
    PRINT, 'Standard errors: ', sigma
;    window, 2
;    plot, H_array, thirdHar_array
;    oplot, H_array, thirdHar_array, PSYM=4, COLOR = 500
;    oplot, H_array, thirdHar_array, PSYM=4, COLOR = 500
     
  iplot, H_array, thirdHar_array, PSYM=4, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='3rd Harmonic'
  iplot, H_array, fortynineHar_array, PSYM=4, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='49th Harmonic'
  iplot, H_array, PeakAmp_array, PSYM=4, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Peak Amplitude'


  
END

;GJ, 2021/8/21, 3rd harmonic equation, not fitting well
PRO Harmonics_equation

  HN = 100
  H_array = FINDGEN(HN);from 0 to 99 [10., 30., 90., 3.33, 1.]
  firstHar_array = (20.*H_array - H_array^3)/30.; + H_array^5*5./1908.
  thirdHar_array = H_array^3/30. - H_array^5*5./1272.
  
  window, 0
  plot, H_array, firstHar_array, PSYM=4, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Harmonics'
  oplot, H_array, thirdHar_array, PSYM=2, COLOR = 500

END

;GJ, 2021/8/21, gradient changes
;
PRO Gradient_changes_half_cycle
  ;pre-defined frequency: 33 KHz
  ;sampling frequency 16.5 MHz
  ;sampling time interval 0.06 us (10-6 sec)
  X = 2*!PI/100. * FINDGEN(50)
  T = 1./16.5; us 10-6 sec
  time = FINDGEN(N_ELEMENTS(X)-1)*T
  z_gradient = FINDGEN(256)+1
  x_gradient = CEIL(z_gradient/15)*10.
  y_gradient = CEIL(z_gradient/30)*10.
  a0_shape = FINDGEN(256*50)*0. + 20.
  z_shape = FINDGEN(256*50)+1
  z_gradient_shape = FINDGEN(256*50)+1
  x_shape = FINDGEN(256*50)+1
  y_shape = FINDGEN(256*50)+1
  FOR i=0, 255 DO BEGIN
    IF i MOD 2 THEN BEGIN
      a0_shape[i*50:(i*50+49)] = 20.*cos(X)
      z_shape[i*50:(i*50+49)] = z_gradient[i]*cos(X)
      z_gradient_shape[i*50:(i*50+49)] = z_gradient[i]
      x_shape[i*50:(i*50+49)] = x_gradient[i]*cos(X)
      y_shape[i*50:(i*50+49)] = y_gradient[i]*cos(X)
    ENDIF ELSE BEGIN
      a0_shape[i*50:(i*50+49)] = 20.*cos(X+!PI)
      z_shape[i*50:(i*50+49)] = z_gradient[i]*cos(X+!PI)
      z_gradient_shape[i*50:(i*50+49)] = z_gradient[i]
      x_shape[i*50:(i*50+49)] = x_gradient[i]*cos(X+!PI)
      y_shape[i*50:(i*50+49)] = y_gradient[i]*cos(X+!PI)
    ENDELSE
  ENDFOR
  
  iplot, z_shape[0:2000], YRANGE=[-80, 80], title='z-gradient Field'
  iplot, z_gradient_shape[0:2000], YRANGE=[-80, 80], title='z-gradient'
  iplot, x_shape[0:2000], YRANGE=[-80, 80], title='x-gradient Field'
  iplot, y_shape[0:2000], YRANGE=[-80, 80], title='y-gradient Field'
  iplot, a0_shape[0:2000], YRANGE=[-80, 80], title='baseline A0'
END

;GJ, 2021/8/21, gradient changes
;GJ, 2021/8/22, full cycle for another gradient change
;GJ, 2021/8/23, using 50mT/m as the 1.5T MRI scanner
PRO Gradient_changes_full_cycle
  ;pre-defined frequency: 3.3 KHz
  ;sampling frequency 1.65 MHz;;;;;;;;16.5 MHz
  ;sampling time interval 0.6 us (10-6 sec)
  Nsampling = 500;samples per cycle
  X = 2.*!PI/Nsampling * FINDGEN(Nsampling)
  
  Nx=32
  delta_T = 1./1.65; us 10-6 sec
  time = FINDGEN(Nx*Nsampling)*delta_T/1000. ; in ms
  max_gradient = 50.; mT/m
  min_gradient = -50.; mT/m
  FOV = 0.2;m
  z_gradient = (FINDGEN(Nx))/(Nx-1)*max_gradient*2. - max_gradient
;  x_gradient = CEIL(z_gradient/5)*10.
;  y_gradient = CEIL(z_gradient/10)*10.
  A0 = 15.; mT
  a0_shape = FINDGEN(Nx*Nsampling)*0.
  z_shape = FINDGEN(Nx*Nsampling)*0.
  z_gradient_shape_field_FOV = FINDGEN(Nx*Nsampling)*0.
  total_gradient_field_FOV = FINDGEN(Nx*Nsampling)*0.;x=FOV, 0.2m
;  x_shape = FINDGEN(256*100)+1
;  y_shape = FINDGEN(256*100)+1
  FOR i=0, Nx-1 DO BEGIN
    a0_shape[i*Nsampling:(i*Nsampling+Nsampling-1)] = -A0*cos(X)
    z_shape[i*Nsampling:(i*Nsampling+Nsampling-1)] = -z_gradient[i]*cos(X)
    z_gradient_shape_field_FOV[i*Nsampling:(i*Nsampling+Nsampling-1)] = -z_gradient[i]*cos(X)*FOV
    total_gradient_field_FOV[i*Nsampling:(i*Nsampling+Nsampling-1)] = -A0*cos(X) -z_gradient[i]*cos(X)*FOV/2.
;    x_shape[i*100:(i*100+99)] = x_gradient[i]*cos(X)
;    y_shape[i*100:(i*100+99)] = y_gradient[i]*cos(X)
  ENDFOR
  
  iplot, time, a0_shape, YRANGE=[-40, 40], XRANGE=[0, max(time)], ytitle='field [mT]', xtitle='time [ms]', title='baseline Field A(t)'
  iplot, time, z_shape, YRANGE=[-80, 80], XRANGE=[0, max(time)], ytitle='gradient [mT/m]', xtitle='time [ms]', title='Gradient G(t)'
  iplot, time, z_gradient_shape_field_FOV, YRANGE=[-40, 40], XRANGE=[0, max(time)], ytitle='field [mT]', xtitle='time [ms]', title='Gradient Field H(t) at Edge'
  iplot, time, total_gradient_field_FOV, YRANGE=[-40, 40], XRANGE=[0, max(time)], ytitle='field [mT]', xtitle='time [ms]', title='Combined Field H(t) at Edge'
;  iplot, x_shape[0:2000], YRANGE=[-80, 80], title='x-gradient Field'
;  iplot, y_shape[0:2000], YRANGE=[-80, 80], title='y-gradient Field'
  
END


;evaluate the relationship between magnetic field strength and signal FWHM
;GJ, Xidian University, 2021/6/15
PRO Magnetic_field_Signal_relationship

  N=512
  H_array = (FINDGEN(N)+1.);*2.
  FWHM_array = DBLARR(N)*0.
  Amplitude_array = DBLARR(N)*0.
  AUC_array = DBLARR(N)*0.
  Mrange_array = DBLARR(N)*0.

  X = 2*!PI/1000. * FINDGEN(1200)
  FOR i=0, N_ELEMENTS(H_array)-1 DO BEGIN
    H = -COS(X)*H_array[i]*0.1
    M = 1.0*(1./TANH(H) - 1./H)
    Mrange_array[i]=MAX(M, /NAN)-MIN(M, /NAN)
    Nelem = N_ELEMENTS(M)
    signal = M[1:*]-M[0:Nelem-1]
    
    IF i LT 20 THEN BEGIN
      M = SMOOTH(M, 40, /NAN) ;to avoid the unsmoothed points
      signal = SMOOTH(signal, 40, /NAN) ;to avoid the unsmoothed points
    ENDIF
    
;    IF i EQ 0 THEN BEGIN
;      plot, X, M, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[-1.1, 1.1]
;    ENDIF ELSE BEGIN
;      oplot, X, M, COLOR = 0;, YRANGE=[-1.1, 1.1]
;    ENDELSE
;    wait, 0.4
;    oplot, X[1:*], 10.*(M[1:*]-M[0:Nelem-1]), COLOR = 0;, YRANGE=[-0.05, 0.05]
;    wait, 0.4
    
    
    time_range = where(X GT 3.14-1. AND X LT 3.14+1., NCounts)
    signal_range = signal[time_range+1]
    AUC_array[i] = TOTAL(signal, /NAN)
    initial = 0; 0; 100
    IF i GE initial AND i LE 455 THEN BEGIN
      IF i EQ initial THEN BEGIN
        plot, time_range, signal_range, title='Signal', COLOR = 0, YRANGE=[0, 0.1], BACKGROUND = 'FFFFFF'x;, BACKGROUND = 'FFFFFF'x
      ENDIF ELSE BEGIN
        oplot, time_range, signal_range, COLOR = i*500
      ENDELSE
      wait, 0.01
    ENDIF
  
 ;   result = POLY_FIT(time_range, signal_range, 15, MEASURE_ERRORS=measure_errors, SIGMA=sigma, YFIT=yfit)
    max_signal = MAX(signal_range, maxind)
    Amplitude_array[i] = max_signal
    lower_signal_sub = ABS(signal_range[0:maxind]-max_signal/2)
    higher_signal_sub = ABS(signal_range[maxind+1:*]-max_signal/2)
  ;  plot, time_range, [lower_signal_sub, higher_signal_sub], title='signal minus half_max', BACKGROUND = 'FFFFFF'x, COLOR = 0
    min_lower_signal_sub = MIN(lower_signal_sub, lower_index)
    min_higher_signal_sub = MIN(higher_signal_sub, higher_index)
    FWHM = higher_index + maxind - lower_index
    FWHM_array[i] = FWHM*(X[1]-X[0])
;    print, signal
  ENDFOR

;;plot FWHM
;;window, 1
;;iplot, H_array[100:200], ALOG(FWHM_array[100:200]), title='H vs ALOG(FWHM)' 
;iplot, H_array[100:200], FWHM_array[100:200], title='H vs FWHM'
;iplot, H_array, FWHM_array, title='H vs FWHM' 
;;print, transpose(H_array)
;;print, transpose(FWHM_array)
;;check their relationship

;plot AUC
;window, 2
iplot,  H_array[200:455], AUC_array[200:455], YRANGE=[0., 1.], title='H vs Area of Peak'
iplot, H_array, AUC_array, title='H vs AUC'

;window, 3
;iplot,  H_array[100:200], AUC_array[100:200]/Mrange_array[100:200], YRANGE=[0.,0.5], title='H vs Area of Peak/M Range'
;print, transpose(AUC_array)

;plot amplitude
iplot, H_array[200:455], Amplitude_array[200:455], title='H vs Amplitude'
iplot, H_array, Amplitude_array, title='H vs Amplitude'

END

;GJ, 2021/8/26; check the amplitude;
;GJ, 2021/9/12; check the 3rd harmonic with a flat pulse portion
;GJ, 2021/10/1; fix some bugs with a flat pulse portion
PRO signal_cal_vs_H_flat_offset_particle_info
  
;  ;plan 1
;  N = 61.;H array elements
;  H_array=FINDGEN(N)/2.-15.;-6 to 6 mT/u0
  
  ;plan 2
  N = 25.;H array elements
  H_array=FINDGEN(N)/2.-(N-1.)/4.;-6 to 6 mT/u0

  gradient = 10. ; mT/cm
  FOV = 20.
  N_offset = FOV*2.+1.; H offset
  H_offset_array=(FINDGEN(N_offset)-FLOOR(N_offset/2))*gradient*0.5
    
  N_flat = 10;
  flat_array = (0.5 + FINDGEN(N_flat)/N_flat/2.)*100.;[0.5, 0.6, 0.7, 0.8, 0.9, 0.96]
  N_particle_size = 10
  particle_size_array = 20. + FINDGEN(N_particle_size)*2.;[21.0, 25.0, 27.0, 32.0, 35.0]

  Amp_array=FINDGEN(N, N_flat, N_particle_size, N_offset)*0.
  thirdHar_array=FINDGEN(N, N_flat, N_particle_size, N_offset)*0.
  AUC_array=FINDGEN(N, N_flat, N_particle_size, N_offset)*0.
  filtered_Amp_array=FINDGEN(N, N_flat, N_particle_size, N_offset)*0.
  filtered_thirdHar_array=FINDGEN(N, N_flat, N_particle_size, N_offset)*0.
  filtered_AUC_array=FINDGEN(N, N_flat, N_particle_size, N_offset)*0.
  
  FOR m=0, N_offset-1 DO BEGIN
    FOR k=0, N_particle_size-1 DO BEGIN
      FOR j=0, N_flat-1 DO BEGIN
        FOR i=0, N-1 DO BEGIN
          flat_portion = flat_array[j]/100.;0.4;5;8;96;0.5;0.96;5
          H_value_shape_freq = [H_array[i], flat_portion, 2.5, H_offset_array[m]]; 0.96 [H field in mT, flat_portion, frequency in kHz]
          ;H_value_shape_freq = [1.0, 0.98, 2.5]
          plotYes_temp = 0;
          tracer_info = [0, particle_size_array[k], plotYes_temp] ;30 nm is the particle size [type, particle size in nm]
          a = signal_vs_H(H_value_shape_freq,tracer_info,0);,/plotYes)
          ;    H_flat_portion = 0.6;9
          ;    a = signal_vs_H(H_array[i],H_flat_portion,0)
          ;Amp_array[i] = a[0, 0]
          thirdHar_array[i, j, k, m] = a[0, 3]
          AUC_array[i, j, k, m] = a[0, 1]
          ;filtered_Amp_array[i] = a[1, 0]
          filtered_thirdHar_array[i, j, k, m] = a[1, 3]
          filtered_AUC_array[i, j, k, m] = a[1, 1]
          print, 'i=',i,', j=',j,', k=',k, 'm=',m
          ;    wait, 0.5
        ENDFOR
      ENDFOR
    ENDFOR
  ENDFOR
    
  peak_filtered_AUC_array=FINDGEN(N_flat, N_particle_size, N_offset)*0.
  peakH_filtered_AUC_array=FINDGEN(N_flat, N_particle_size, N_offset)*0.
  
  FOR m=0, N_offset-1 DO BEGIN
    FOR k=0, N_particle_size-1 DO BEGIN
      FOR j=0, N_flat-1 DO BEGIN
        maxind = 0;initiate
        peak_filtered_AUC_array[j, k, m] = MAX(filtered_AUC_array[*, j, k, m], maxind)
        peakH_filtered_AUC_array[j, k, m] =H_array[maxind]
      ENDFOR
    ENDFOR
  ENDFOR

  
  cgsurface, peak_filtered_AUC_array[*,*,FOV], flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=6, /Brewer, TITLE='Peak Signal AUC', $
    XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'Peak Signal AUC [A.U.]', ZRANGE=[0,1000];,/CONSTRAIN_ASPECT

  cgsurface, peakH_filtered_AUC_array[*,*,FOV], flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=10, /Brewer, TITLE='Field Amplitude with Peak Signal AUC', $
    XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'Field Amplitude [mT/u0]', ZRANGE=[0,10];,/CONSTRAIN_ASPECT


;  ;iplot, H_array, Amp_array, title='Amplitude', xtitle='H [mT]'
;  iplot, H_array, thirdHar_array, title='3rd Harmonic', xtitle='H [mT]'
;  ;iplot, H_array, filtered_Amp_array, title='Filtered Amplitude', xtitle='H [mT]'
;  iplot, H_array, filtered_thirdHar_array, title='Filtered 3rd Harmonic', xtitle='H [mT]'
  
  ;compare the 3rd harmonic with and without relaxation
;  window, 1
;  plot, H_array, thirdHar_array, title='3rd Harmonic', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x, COLOR = 0
;  oplot, H_array, filtered_thirdHar_array, PSYM=2, color=3500
  
  ;compare the 3rd harmonic with and without relaxation
  j=3
  cgsurface, REFORM(filtered_AUC_array[*, j, *, FOV]), H_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Flat portion: '+STRING(flat_array[j], format='(I2)')+'%]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Size [nm]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT
  
  j=6
  cgsurface, REFORM(filtered_AUC_array[*, j, *, FOV]), H_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Flat portion: '+STRING(flat_array[j], format='(I2)')+'%]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Size [nm]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT
  
  j=9
  cgsurface, REFORM(filtered_AUC_array[*, j, *, FOV]), H_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Flat portion: '+STRING(flat_array[j], format='(I2)')+'%]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Size [nm]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

    
  ;compare the 3rd harmonic with and without relaxation
  k=3
  cgsurface, REFORM(filtered_AUC_array[*, *, k, FOV]), H_array, flat_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Particle size: '+STRING(particle_size_array[k], format='(I2)')+' nm]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Flat ratio [%]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  k=6
  cgsurface, REFORM(filtered_AUC_array[*, *, k, FOV]), H_array, flat_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Particle size: '+STRING(particle_size_array[k], format='(I2)')+' nm]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Flat ratio [%]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  k=9
  cgsurface, REFORM(filtered_AUC_array[*, *, k, FOV]), H_array, flat_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Particle size: '+STRING(particle_size_array[k], format='(I2)')+' nm]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Flat ratio [%]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  
;  mySurface = SURFACE(REFORM(filtered_AUC_array[*, 0, *]), H_array, particle_size_array, $
;    TITLE='Relaxation Decay Signal AUC')
;  ; Overlay the contour data.
;  myContour = CONTOUR(REFORM(filtered_AUC_array[*, 0, *]), H_array, particle_size_array, $
;    PLANAR=0, /OVERPLOT);N_LEVELS=10, 
;  ax = mySurface.AXES
;  ax[0].TITLE = 'H [mT/u0]'
;  ax[1].TITLE = 'Size [nm]'
;  ax[2].TITLE = 'Signal AUC [A.U.]'
;
;
;  window, 2
;  plot, H_array, filtered_AUC_array[*, 0, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[0])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 0, k], thick = 1, color=350000*(k+1), PSYM=0
;
;  window, 3
;  plot, H_array, filtered_AUC_array[*, 2, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[2])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 1, k], thick = 1, color=350000*(k+1), PSYM=0
;  
;  window, 4
;  plot, H_array, filtered_AUC_array[*, 4, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[4])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 2, k], thick = 1, color=350000*(k+1), PSYM=0
;
;  window, 5
;  plot, H_array, filtered_AUC_array[*, 6, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[6])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 3, k], thick = 1, color=350000*(k+1), PSYM=0
;
;  window, 6
;  plot, H_array, filtered_AUC_array[*, 8, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[8])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 4, k], thick = 1, color=350000*(k+1), PSYM=0
;
;  window, 7
;  plot, H_array, filtered_AUC_array[*, 9, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[9])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 5, k], thick = 1, color=350000*(k+1), PSYM=0
;  
;  window, 8
;  plot, H_array, filtered_AUC_array[*, N_flat-1, 0], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[0])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 0], thick = 1, color=350000*(j+1), PSYM=0
;
;  window, 9
;  plot, H_array, filtered_AUC_array[*, N_flat-1, 3], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[3])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 1], thick = 1, color=350000*(j+1), PSYM=0
;   
;  window, 10
;  plot, H_array, filtered_AUC_array[*, N_flat-1, 6], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[6])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 2], thick = 1, color=350000*(j+1), PSYM=0
;  
;  window, 11
;  plot, H_array, filtered_AUC_array[*, N_flat-1, 9], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[9])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 3], thick = 1, color=350000*(j+1), PSYM=0


END

;GJ, 2021/8/26; check the amplitude;
;GJ, 2021/9/12; check the 3rd harmonic with a flat pulse portion
;GJ, 2021/10/1; fix some bugs with a flat pulse portion
;GJ, 2021/11/18; if offset is set, then calculate the signal with different offsets
PRO signal_features_vs_H_test, Hoffset=Hoffset
  
  ;calculate the offset if Hoffest is set and return
  IF KEYWORD_SET(Hoffset) THEN BEGIN
    signal_cal_vs_H_flat_offset_particle_info
    RETURN
  ENDIF
  
  N = 61.;H array elements
  H_array=FINDGEN(N)/2.-(N-1.)/4.;-15 to 15 mT/u0
  
  N_flat = 10;
  flat_array = (0.5 + FINDGEN(N_flat)/N_flat/2.)*100.;[0.5, 0.6, 0.7, 0.8, 0.9, 0.96]
  N_particle_size = 10
  particle_size_array = 20. + FINDGEN(N_particle_size)*2.;[21.0, 25.0, 27.0, 32.0, 35.0]

  Amp_array=FINDGEN(N, N_flat, N_particle_size)*0.
  thirdHar_array=FINDGEN(N, N_flat, N_particle_size)*0.
  AUC_array=FINDGEN(N, N_flat, N_particle_size)*0.
  filtered_Amp_array=FINDGEN(N, N_flat, N_particle_size)*0.
  filtered_thirdHar_array=FINDGEN(N, N_flat, N_particle_size)*0.
  filtered_AUC_array=FINDGEN(N, N_flat, N_particle_size)*0.

  FOR k=0, N_particle_size-1 DO BEGIN
    FOR j=0, N_flat-1 DO BEGIN
      FOR i=0, N-1 DO BEGIN
        flat_portion = flat_array[j]/100.;0.4;5;8;96;0.5;0.96;5
        H_value_shape_freq = [H_array[i], flat_portion, 2.5]; 0.96 [H field in mT, flat_portion, frequency in kHz]
        ;H_value_shape_freq = [1.0, 0.98, 2.5]
        plotYes = 0;
        tracer_info = [0, particle_size_array[k], plotYes] ;30 nm is the particle size [type, particle size in nm]
        a = signal_vs_H(H_value_shape_freq,tracer_info,0)
        ;    H_flat_portion = 0.6;9
        ;    a = signal_vs_H(H_array[i],H_flat_portion,0)
        ;Amp_array[i] = a[0, 0]
        thirdHar_array[i, j, k] = a[0, 3]
        AUC_array[i, j, k] = a[0, 1]
        ;filtered_Amp_array[i] = a[1, 0]
        filtered_thirdHar_array[i, j, k] = a[1, 3]
        filtered_AUC_array[i, j, k] = a[1, 1]
        print, 'i=',i,', j=',j,', k=',k
        ;    wait, 0.5
      ENDFOR
    ENDFOR
  ENDFOR

  peak_filtered_AUC_array=FINDGEN(N_flat, N_particle_size)*0.
  peakH_filtered_AUC_array=FINDGEN(N_flat, N_particle_size)*0.
  FOR k=0, N_particle_size-1 DO BEGIN
    FOR j=0, N_flat-1 DO BEGIN
      maxind = 0;initiate
      peak_filtered_AUC_array[j, k] = MAX(filtered_AUC_array[*, j, k], maxind)
      peakH_filtered_AUC_array[j, k] =H_array[maxind]
    ENDFOR
  ENDFOR

  cgsurface, peak_filtered_AUC_array, flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=6, /Brewer, TITLE='Peak Signal AUC', $
    XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'Peak Signal AUC [A.U.]', ZRANGE=[0,1000];,/CONSTRAIN_ASPECT

  cgsurface, peakH_filtered_AUC_array, flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=10, /Brewer, TITLE='Field Amplitude with Peak Signal AUC', $
    XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'Field Amplitude [mT/u0]', ZRANGE=[0,10];,/CONSTRAIN_ASPECT


  ;  ;iplot, H_array, Amp_array, title='Amplitude', xtitle='H [mT]'
  ;  iplot, H_array, thirdHar_array, title='3rd Harmonic', xtitle='H [mT]'
  ;  ;iplot, H_array, filtered_Amp_array, title='Filtered Amplitude', xtitle='H [mT]'
  ;  iplot, H_array, filtered_thirdHar_array, title='Filtered 3rd Harmonic', xtitle='H [mT]'

  ;compare the 3rd harmonic with and without relaxation
  ;  window, 1
  ;  plot, H_array, thirdHar_array, title='3rd Harmonic', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  ;  oplot, H_array, filtered_thirdHar_array, PSYM=2, color=3500

  ;compare the 3rd harmonic with and without relaxation
  j=3
  cgsurface, REFORM(filtered_AUC_array[*, j, *]), H_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Flat portion: '+STRING(flat_array[j], format='(I2)')+'%]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Size [nm]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  j=6
  cgsurface, REFORM(filtered_AUC_array[*, j, *]), H_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Flat portion: '+STRING(flat_array[j], format='(I2)')+'%]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Size [nm]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  j=9
  cgsurface, REFORM(filtered_AUC_array[*, j, *]), H_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Flat portion: '+STRING(flat_array[j], format='(I2)')+'%]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Size [nm]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT


  ;compare the 3rd harmonic with and without relaxation
  k=3
  cgsurface, REFORM(filtered_AUC_array[*, *, k]), H_array, flat_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Particle size: '+STRING(particle_size_array[k], format='(I2)')+' nm]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Flat ratio [%]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  k=6
  cgsurface, REFORM(filtered_AUC_array[*, *, k]), H_array, flat_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Particle size: '+STRING(particle_size_array[k], format='(I2)')+' nm]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Flat ratio [%]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT

  k=9
  cgsurface, REFORM(filtered_AUC_array[*, *, k]), H_array, flat_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='Signal AUC [Particle size: '+STRING(particle_size_array[k], format='(I2)')+' nm]', $
    XTITLE = 'H [mT/u0]', YTITLE = 'Flat ratio [%]', ZTITLE = 'Signal AUC [A.U.]', ZRANGE=[-1000,1000];,/CONSTRAIN_ASPECT


  ;  mySurface = SURFACE(REFORM(filtered_AUC_array[*, 0, *]), H_array, particle_size_array, $
  ;    TITLE='Relaxation Decay Signal AUC')
  ;  ; Overlay the contour data.
  ;  myContour = CONTOUR(REFORM(filtered_AUC_array[*, 0, *]), H_array, particle_size_array, $
  ;    PLANAR=0, /OVERPLOT);N_LEVELS=10,
  ;  ax = mySurface.AXES
  ;  ax[0].TITLE = 'H [mT/u0]'
  ;  ax[1].TITLE = 'Size [nm]'
  ;  ax[2].TITLE = 'Signal AUC [A.U.]'
  ;
  ;
  ;  window, 2
  ;  plot, H_array, filtered_AUC_array[*, 0, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[0])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 0, k], thick = 1, color=350000*(k+1), PSYM=0
  ;
  ;  window, 3
  ;  plot, H_array, filtered_AUC_array[*, 2, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[2])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 1, k], thick = 1, color=350000*(k+1), PSYM=0
  ;
  ;  window, 4
  ;  plot, H_array, filtered_AUC_array[*, 4, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[4])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 2, k], thick = 1, color=350000*(k+1), PSYM=0
  ;
  ;  window, 5
  ;  plot, H_array, filtered_AUC_array[*, 6, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[6])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 3, k], thick = 1, color=350000*(k+1), PSYM=0
  ;
  ;  window, 6
  ;  plot, H_array, filtered_AUC_array[*, 8, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[8])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 4, k], thick = 1, color=350000*(k+1), PSYM=0
  ;
  ;  window, 7
  ;  plot, H_array, filtered_AUC_array[*, 9, N_particle_size-1], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [flat_part = '+STRING(flat_array[9])+']', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR k=0, N_particle_size-2 DO oplot, H_array, filtered_AUC_array[*, 5, k], thick = 1, color=350000*(k+1), PSYM=0
  ;
  ;  window, 8
  ;  plot, H_array, filtered_AUC_array[*, N_flat-1, 0], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[0])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 0], thick = 1, color=350000*(j+1), PSYM=0
  ;
  ;  window, 9
  ;  plot, H_array, filtered_AUC_array[*, N_flat-1, 3], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[3])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 1], thick = 1, color=350000*(j+1), PSYM=0
  ;
  ;  window, 10
  ;  plot, H_array, filtered_AUC_array[*, N_flat-1, 6], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[6])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 2], thick = 1, color=350000*(j+1), PSYM=0
  ;
  ;  window, 11
  ;  plot, H_array, filtered_AUC_array[*, N_flat-1, 9], YRANGE=[-1000,1000], PSYM=0, thick=1, color=0, title='Relaxation AUC [particle size = '+STRING(particle_size_array[9])+'nm]', xtitle='H [mT/u0]', BACKGROUND = 'FFFFFF'x
  ;  FOR j=0, N_flat-2 DO oplot, H_array, filtered_AUC_array[*, j, 3], thick = 1, color=350000*(j+1), PSYM=0


END

;GJ, 2021/11/25, plot the pulsed excitations with different flat portions
PRO plot_pulsed_H_shapes

  H_value = 5.0;mT
  H_offset = 0.0
  ;Sampling frequency
  frequency = 2.5;10;33.;3.3; kHz
  H_flat_portion_array = [0.5, 0.6, 0.7, 0.8, 0.9]; 
  
  ;define the period based on frequency
  period = 1.0/frequency; ms
  Nsample = 1000.

  delta_t = period/Nsample * 1000; us

  ;define the cycle index
  X = 2*!PI/Nsample * FINDGEN(1.2*NSample)
  Nelem = N_ELEMENTS(X)
  ;GJ, 2021/11/25, define the time series
  T = FINDGEN(Nelem-1) * delta_t ; us
  H_array = DBLARR(5, Nelem)
  
  FOR i=0,4 DO BEGIN
    H_curve_portion = 1.- H_flat_portion_array[i]
    flat_Nsample = FLOOR(NSample/4.*H_flat_portion_array[i])
    H1 = (DBLARR(flat_Nsample)*0.+H_value) * (-1.)
    H2 = -REVERSE(FINDGEN(NSample/4-flat_Nsample)/(NSample/4-flat_Nsample-1)*H_value)
    H3 = (FINDGEN(NSample/4-flat_Nsample+1)/(NSample/4-flat_Nsample))*H_value
    H4 = -H1
    H1st = [H1, H2, H3[1:*], H4]
    H_cycle = [H1st, REVERSE(H1st)]
    H = [H_cycle, H_cycle[0:199]]
    H_array[i, *] = H
  ENDFOR
  
  ; Create the plot.
  cgPlot, T[0:999], H_array[0, 0:999], Thick=2, Linestyle = 0, Color='sea green', YRANGE=[-MAX(H)*1.1, MAX(H)*1.5], YStyle=1, $
    XTitle='Time / us', YTitle='Applied Field / mT'
  cgPlot, T[0:999], H_array[1, 0:999],  Thick=2, Linestyle = 1, Color='orange red', /Overplot
  cgPlot, T[0:999], H_array[2, 0:999],  Thick=2, Linestyle = 2, Color='goldenrod', /Overplot
  cgPlot, T[0:999], H_array[3, 0:999],  Thick=2, Linestyle = 3, Color='dark orchid', /Overplot
  cgPlot, T[0:999], H_array[4, 0:999],  Thick=2, Linestyle = 4, Color='blue', /Overplot

  ; Create the legend with NASA Astronomy routine AL_LEGEND.
  items = ['50%', '60%','70%','80%','90%']
  Linestyles = [0, 1, 2, 3, 4]
  ;psyms = [-15, -16,-17,-18]
  colors = ['sea green', 'orange red','goldenrod','dark orchid', 'blue']

  ; Location of legend in data coordinates.
  yloc = (!Y.CRange[1] - !Y.CRange[0]) * 0.95 + !Y.CRange[0]
  xloc = (!X.CRange[1] - !X.CRange[0]) * 0.05 + !X.CRange[0]

  ; Add the legend.
  cgLegend, Title=items, Thick=2, CHARTHICK=1, Linestyles = Linestyles, Color=colors,$
    Location=[xloc,yloc], /Data
   
   
   
END


;GJ, 2021/8/25, core function
;GJ, 2021/9/2, adding signal output and options
;GJ, 2021/9/12, add pulsed excitation
;GJ, 2021/10/17, add more parameters
;calcualte signal and feature out of H
;
;input: H_value_shape_freq = [H_value, flat_portion, frequency]
;       tracer_info = [tracer_type, particle_size in nm, plotYes]
; if s_and_f = 0, output feature; s_and_f = 1, output signal and filtered signal
; if signal_array_for_feature is given, and s_and_f = 0, output feature
;       
;output: signal_feature [2, 5] or signal[4, 1199]
;   signal_feature[0, 0:4] = [Amplitude, AUC, FWHM, thirdHar, H_value]
;   signal_feature[1, 0:4] = [filteredAmplitude, filteredAUC, filteredFWHM, filteredthirdHar, H_value]
;   signal[0, 0:1199] is the time in us
;   signal[1, 0:1199] is the non-adiabatic signal without relaxation
;   signal[2, 0:1199] is the adiabatic signal with relaxation
;   signal[3, 0:1199] is zero

;example 1:
;test_signal = signal_vs_H(H_info, tracer_info, 1, signal_array_output) ; test_signal is the signal calculated from H_info and tracer_info
;bb = signal_vs_H(H_info, tracer_info, 0, test_signal); test_signal as the input, bb is the features from aa, not from H_info or tracer_info
;cc = signal_vs_H(H_info, tracer_info, 0); calculate the feature based on H and tracer info

;example 2:
;H_value_shape_freq = [H_value, flat_portion, frequency]; 0.96 [H field in mT, flat_portion, frequency in kHz]
;plotYes = 0;
;tracer_info = [tracer_type, particle_size, plotYes] ;trace type: 0 is for resovist, 30 nm is the particle size [type, particle size in nm], plot or not
;para = signal_vs_H(H_value_shape_freq,tracer_info,0)
;;make plotYes as a keyword
;example 3:
;signal_1d[i] = (signal_vs_H(H_value_shape_freq,tracer_info,0, signal_total))[1, 1];get the filtered signal AUC
FUNCTION signal_vs_H, H_value_shape_freq, tracer_info, s_and_f, signal_array_for_feature, plotYes=plotYes
  
  IF N_ELEMENTS(H_value_shape_freq) EQ 1 THEN BEGIN
    H_value = H_value_shape_freq
    H_offset = 0.0
    ;Sampling frequency
    frequency = 3.3;10;33.;3.3; kHz
    IF N_ELEMENTS(tracer_info) EQ 1 THEN BEGIN
      H_flat_portion = tracer_info;the old version
    ENDIF ELSE BEGIN
      H_flat_portion = 0.9; by defaulty
    ENDELSE
  ENDIF ELSE BEGIN
    H_value = H_value_shape_freq[0]
    H_flat_portion = H_value_shape_freq[1]
    frequency = H_value_shape_freq[2]
    ;add the H offset component
    IF N_ELEMENTS(H_value_shape_freq) GT 3 THEN H_offset = H_value_shape_freq[3] ELSE H_offset = 0.0
  ENDELSE

  
  IF N_ELEMENTS(tracer_info) EQ 1 THEN BEGIN
    H_flat_portion = tracer_info;the old version
    tracer_type = 0; by default is Resovist
    particle_size = 30. ; in nm
;    plotYes=0
  ENDIF ELSE BEGIN
    tracer_type = tracer_info[0]
    particle_size = tracer_info[1]
;    IF N_ELEMENTS(tracer_info) EQ 3 THEN plotYes = tracer_info[2] ELSE plotYes=0
  ENDELSE
  
;  print, 'H value: ', H_value, ' mT'
;  print, 'Flat ratio: ', H_flat_portion
;  print, 'Frequency: ', frequency, 'kHz'
;  print, 'Particle size: ', particle_size, 'nm'
  
  ;define the period based on frequency
  period = 1.0/frequency; ms
  Nsample = 1000.
  
  delta_t = period/Nsample * 1000; us
  
  ;define the cycle index
  X = 2*!PI/Nsample * FINDGEN(1.2*NSample)
  Nelem = N_ELEMENTS(X)
  ;GJ, 2021/11/25, define the time series
  T = FINDGEN(Nelem-1) * delta_t ; us
 
  ;in case the H_value is zero
  IF ABS(H_value) LT 0.001 THEN BEGIN
    ;T = FINDGEN(1.2*Nsample-1) * delta_t ; us
    M = FINDGEN(1.2*Nsample)*0.
    signal = M[1:*]-M[0:1.2*Nsample-2]
    filteredsignal = signal
    ;only output t, signal, filtered signal
    signal_array_for_feature = DBLARR(4, N_ELEMENTS(T))*0.
    IF s_and_f EQ 1 THEN BEGIN
      ;s_and_f = 0, output feature; s_and_f = 1, output signal and filtered signal
      signal_array_for_feature = ([TRANSPOSE(T), TRANSPOSE(signal), TRANSPOSE(filteredsignal),TRANSPOSE(T*0.)]);
    ENDIF
    return, DBLARR(2,5)*0.; return the parameters with 0    
  ENDIF
 
  ;the portion of curve is flat, the rest is linear
  ;H_flat_portion = 0.9
  IF H_flat_portion GT 0.01 AND H_flat_portion LT 0.9999 THEN BEGIN
    H_curve_portion = 1.- H_flat_portion
    flat_Nsample = FLOOR(NSample/4.*H_flat_portion)
    H1 = (DBLARR(flat_Nsample)*0.+H_value) * (-1.)
    H2 = -REVERSE(FINDGEN(NSample/4-flat_Nsample)/(NSample/4-flat_Nsample-1)*H_value)
    H3 = (FINDGEN(NSample/4-flat_Nsample+1)/(NSample/4-flat_Nsample))*H_value
    H4 = -H1
    H1st = [H1, H2, H3[1:*], H4]
    H_cycle = [H1st, REVERSE(H1st)]
    H = [H_cycle, H_cycle[0:199]]
  ENDIF ELSE BEGIN
    ;this is a cosine curve
    H = -COS(X)*H_value
  ENDELSE
  ;GJ, 2021/11/17
  H = H + H_offset
  
  ;GJ, 2021/11/25, plot the flat portion
  ;  plot, T, -H, title='H field', xtitle='Time / us', ytitle='Applied FIeld / mT/u0', BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[-MAX(H)*1.1, MAX(H)*1.1]
  
  ;Perimag, Synomg, and Perimag
  ;collection of different beta values
;  size_array = [15., 20., 30., 40.]
;  beta_array = [4.16/25.0, 4.16/12.0, 4.16/3.0, 4.16/2.0]
;  beta_interpo = INTERPOL(beta_array, size_array, particle_size)
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  ;Msat = 446.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  ;beta_particle = !PI/6.*Msat_kAm*0.0001*(particle_size^3)/309.65/1.380469; in 1/mT
  FWHM_particle = 4.16/beta_particle/3.5; mT/mm
  ;m_particle = (Msat*10000./4.) * (particle_size^3)/6. /1000000. ;KA*mm2*10-27
;  M_particle = 446.*!PI/6.*(particle_size^3)/1000000. ;KA*mm2*10-27

  IF KEYWORD_SET(plotYes) THEN BEGIN
    print, 'particle size =', particle_size, ' [nm]'
    ;  print, 'beta1 = ', beta_particle1, ' [1/mT]'
    print, 'beta = ', beta_particle, ' [1/mT]'
    print, 'FWHM = ', FWHM_particle, ' [mm]'
    ;print, 'magnetic moment m = ', m_particle, ' [KA*mm2*10-27]'
  ENDIF

  
;  beta_Resovist = 4.16/2.0; 2.08 ;for beta*H (mT); 40nm
;  beta_UW33 = 4.16/4.16/25.0; 2.08 ;for beta*H (mT); 15nm
;  beta_40nm = 4.16/2.0;
;  beta_30nm = 4.16/3.0;
;  beta_20nm = 4.16/12.0;
;  beta_15nm = 4.16/25.0;
  
  CASE tracer_type OF
    0: M = Msat_kAm*(1./TANH(beta_particle*H) - 1./(beta_particle*H))
    1: M = Msat_kAm*(1./TANH(beta_UW33*H) - 1./(beta_UW33*H))
;    2: M = Msat_kAm*(1./TANH(beta_Perimag*H) - 1./(beta_Perimag*H))
;    3: M = Msat_kAm*(1./TANH(Synomg*H) - 1./(Synomg*H))
;    4: M = Msat_kAm*(1./TANH(beta_VivoTrax*H) - 1./(beta_VivoTrax*H))
    ELSE: M = Msat_kAm*(1./TANH(beta_particle*H) - 1./(beta_particle*H))
  ENDCASE

  ;deal with finite
  index_nan = where(~FINITE(M), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(M, T, T[index_nan], /NAN)
    M[index_nan] = RESULT
  ENDIF

  Mrange = MAX(M, /NAN)-MIN(M, /NAN)
  signal = (M[1:*]-M[0:Nelem-1])/delta_t; kA/m/us
  signal_sign = 1
  IF TOTAL(Signal[0:Nsample/2-1]) LT 0 THEN BEGIN
    signal = -signal
    signal_sign = -1
  ENDIF
;  signal = SMOOTH(M[1:*]-M[0:Nelem-1],3,/NAN)
  

  ;GJ, 20219/2, added filtered signal
  CASE tracer_type OF
    0: relaxation_time = relaxation_pulsed_excitation(2.*ABS(H_value), particle_size);, /plotYes); 30. is the particle size
    1: relaxation_time = relaxation_time_UW33(2.*ABS(H_value), frequency)
;    2: relaxation_time = relaxation_time_Perimag(ABS(H_value), frequency)
;    3: relaxation_time = relaxation_time_Synomg(ABS(H_value), frequency)
;    4: relaxation_time = relaxation_time_VivoTrax(ABS(H_value), frequency)
    ELSE: relaxation_time = relaxation_pulsed_excitation(2.*ABS(H_value), particle_size);, /plotYes); 30. is the particle size
  ENDCASE
  
  IF N_ELEMENTS(tracer_info) EQ 3 THEN relaxation_time = tracer_info[2];us
  
  ;relaxation_time = [30.0]; us ;DBLARR(5, 100)
  ;  signal_nat = non_adiabatic(T, signal, relaxtion_time)
  filteredsignal = non_adiabatic(T, signal, relaxation_time)
  IF KEYWORD_SET(plotYes) THEN BEGIN
      print, 'H_value =', H_value, '[mT]', 'frequency =',  frequency, '[kHz]'
      print, 'Relaxation_time = ', relaxation_time, '[us]'
      print, 'particle size = ', particle_size, '[nm]'

    ;GJ, 2021/9/15, plot the filtered signal
      window, 0
      plot, T, H/MAX(H)*MAX(signal), thick = 1, title='signal　at H='+STRING(H_value)+ '[mT]', xtitle='time [us]', ytitle='A.U.', BACKGROUND = 'FFFFFF'x, COLOR = 0
      oplot, T, signal*signal_sign, LINESTYLE=3, thick = 2, color=35000
      oplot, T, filteredsignal*signal_sign, thick = 2, color=3500;/MAX(filteredsignal)*MAX(signal)
    ;    wait, 0.1
  ENDIF


  IF s_and_f EQ 1 THEN BEGIN
    ;s_and_f = 0, output feature; s_and_f = 1, output signal and filtered signal
    signal_array_for_feature = TRANSPOSE([[T], [signal*signal_sign], [filteredsignal*signal_sign], [H[0:Nelem-2]]]);
  ENDIF ELSE BEGIN
    ;signal_array_for_feature is given, and s_and_f = 0, output feature
    IF N_ELEMENTS(signal_array_for_feature) GT 10 THEN BEGIN
      Signal = REFORM(signal_array_for_feature[1, *])
      filteredsignal = REFORM(signal_array_for_feature[2, *])
      IF TOTAL(Signal[0:Nsample/2-1]) LT 0 THEN BEGIN
        signal = -signal
        filteredsignal = -filteredsignal
        signal_sign = -1
      ENDIF
    ENDIF
  ENDELSE

  ;using 6 repeating for FFT analysis
  Signal_ext = [Signal[0:Nsample-1], Signal[0:Nsample-1], Signal[0:Nsample-1], Signal[0:Nsample-1], Signal[0:Nsample-1], Signal[0:Nsample-1]]
  N_ext = N_ELEMENTS(Signal_ext)
  T_ext = FINDGEN(N_ext) * delta_t ; us
  fft_comp = ABS(FFT(Signal_ext))
  thirdHar = fft_comp[18];get the 3rd harmonic

  ;  IF H_value LT 20 THEN BEGIN
  ;    M = SMOOTH(M, 40, /NAN) ;to avoid the unsmoothed points
  ;    signal = SMOOTH(signal, 40, /NAN) ;to avoid the unsmoothed points
  ;  ENDIF

  X_range = where(X GT 0. AND X LT 3.14, NCounts)
  signal_range = signal[X_range]
  AUC = TOTAL(signal_range, /NAN)*delta_t ;Area under the curve, total signal * time interval;; kA/m
  
  ;2021/11/8, GJ, fix the bug of irregular signal FWHM calculation
  max_signal = MAX(signal_range, maxind)
  Amplitude = max_signal
  IF maxind GT 2 AND maxind LT N_ELEMENTS(signal_range)-3 THEN BEGIN
    lower_signal_sub = ABS(signal_range[0:maxind]-max_signal/2)
    higher_signal_sub = ABS(signal_range[maxind+1:*]-max_signal/2)
    ;  plot, X_range, [lower_signal_sub, higher_signal_sub], title='signal minus half_max', BACKGROUND = 'FFFFFF'x, COLOR = 0
    min_lower_signal_sub = MIN(lower_signal_sub, lower_index)
    min_higher_signal_sub = MIN(higher_signal_sub, higher_index)
    FWHM_temp = higher_index + maxind - lower_index
    FWHM = FWHM_temp*delta_t ; in us
  ENDIF ELSE BEGIN
    FWHM = 0.
  ENDELSE
  
  
  IF signal_sign EQ -1 THEN BEGIN
    signal_feature = [-Amplitude, -AUC, FWHM, -thirdHar, H_value]
  ENDIF ELSE BEGIN
    signal_feature = [Amplitude, AUC, FWHM, thirdHar, H_value]
  ENDELSE
    
;    print, 'H = ', signal_feature[4]
;    print, 'Amp, AUC(us), FWHM(us), 3rd Harmonic = ', signal_feature[0:3]


  ;GJ, 2021/8/26
  ;add relaxation effects
  IF H_flat_portion GT 0.01 AND H_flat_portion LT 0.96 THEN BEGIN
    ;only calculate the relaxation part
    filteredsignal_range = filteredsignal[(NSample/2-flat_Nsample-2):(NSample/2+flat_Nsample-1)]
  ENDIF ELSE BEGIN
    filteredsignal_range = filteredsignal[X_range]
  ENDELSE
  filteredAUC = TOTAL(filteredsignal_range, /NAN)*delta_t ;Area under the curve, total signal * time interval;; kA/m

  ;2021/11/8, GJ, fix the bug of irregular filtered signal FWHM calculation
  max_filteredsignal = MAX(filteredsignal_range, filteredmaxind)
  filteredAmplitude = max_filteredsignal
  IF filteredmaxind GT 2 AND filteredmaxind LT N_ELEMENTS(filteredsignal_range)-3 THEN BEGIN
    lower_filteredsignal_sub = ABS(filteredsignal_range[0:filteredmaxind]-max_filteredsignal/2)
    higher_filteredsignal_sub = ABS(filteredsignal_range[filteredmaxind+1:*]-max_filteredsignal/2)
    ;  plot, X_range, [lower_signal_sub, higher_signal_sub], title='signal minus half_max', BACKGROUND = 'FFFFFF'x, COLOR = 0
    min_lower_filteredsignal_sub = MIN(lower_filteredsignal_sub, lower_filteredindex)
    min_higher_filteredsignal_sub = MIN(higher_filteredsignal_sub, higher_filteredindex)
    filteredFWHM_temp = higher_filteredindex + filteredmaxind - lower_filteredindex
    filteredFWHM = filteredFWHM_temp*delta_t ; in us
  ENDIF ELSE BEGIN
    filteredFWHM = 0.
  ENDELSE


  ;3rd harmonic
  ;using 6 repeating for FFT analysis
  filteredSignal_ext = [filteredsignal[0:Nsample-1], filteredsignal[0:Nsample-1], filteredsignal[0:Nsample-1], filteredsignal[0:Nsample-1], filteredsignal[0:Nsample-1], filteredsignal[0:Nsample-1]]
  N_ext = N_ELEMENTS(filteredSignal_ext)
  T_ext = FINDGEN(N_ext) * delta_t ; us
  filteredfft_comp = ABS(FFT(filteredSignal_ext))
  filteredthirdHar = filteredfft_comp[18];get the 3rd harmonic
  
  IF signal_sign EQ -1 THEN BEGIN
    filteredsignal_feature = [-filteredAmplitude, -filteredAUC, filteredFWHM, -filteredthirdHar, H_value]
  ENDIF ELSE BEGIN
    filteredsignal_feature = [filteredAmplitude, filteredAUC, filteredFWHM, filteredthirdHar, H_value]
  ENDELSE
  ;filteredsignal_feature = [filteredAmplitude, filteredAUC, filteredFWHM, filteredthirdHar, H_value]
;    print, 'H = ', filteredsignal_feature[4]
;    print, 'Filtered Amp, AUC(us), FWHM(us), 3rd Harmonic = ', filteredsignal_feature[0:3]

  return, TRANSPOSE([[signal_feature], [filteredsignal_feature]]);
  ;  return, filteredsignal_feature
  ;Amplitude
  ;AUC, signal * time interval us
  ;FWHM, time interval in us
  ;thirdHar
END

;@GJ, 2022/7/6, Medical Physics paper to be finished
FUNCTION visocity_glycero_temp, Temperature, Glycerol
  Temperature_array = [25., 30., 35., 40., 45.]
  Glycerol_array = [0.3, 4.1, 7.6, 11.0, 14.2, 21.1, 24.6, 28.0, 31.2, 34.5, 37.7, 41.0, 44.2, 47.4, 50.5]
  Viscosity_array = DBLARR(5, 15)*0.
  Viscosity_array[0, *] = [0.90, 1.01, 1.13, 1.26, 1.40, 1.80, 2.06, 2.36, 2.71, 3.13, 3.63, 4.25, 5.00, 5.93, 7.06]
  RETURN, 0
END

;particle_size = 24.4924;PBS, [25.2579], 1% Glycerol
;T_p = 20.;36.5; degree
;@GJ, 2023/2/12, calculating Brownian relaxation based on pulsed excitation
;@GJ, 2023/3/10, relaxation time with field changes
PRO SynamagD_relaxation_time_Brownian_Neel, eta, B_array, Neel_time_B_array, Brownian_time_eta_B_array
  eta = FINDGEN(401)/100.*1.e-3 ; Pa.s
  IF N_ELEMENTS(D_h) EQ 0 THEN D_h = 54.;
  print, 'D_h: ', D_h, ' nm'
  D_h *= 1.e-9 ; m
  V_h = 1./6. * !PI * D_h^3
  k_B = 1.380469e-23; J/K
  ;T_p = 25.; degree
  T_p_kelvin = 300.;273.15 + T_p ; in Kelvin

  tao_B = 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
  print, 'tao_B: ', tao_B, ' us'
  ;  iplot, eta, tao_B, xtitle='Viscosity [mPa.s]', ytitle='Brownian Relaxation Time [us]'

  M_s = 374000.; J/m3 T (Ms = 374000 J/m3 T)
  IF N_ELEMENTS(D_c) EQ 0 THEN D_c = 16.3; nm
  D_c *= 1.e-9 ;m
  ;D_c = 5.e-9 ;m
  V_c =  1./6. * !PI * D_c^3

  B_array = (FINDGEN(101.)/10.+0.1)*1.e-3 ;T
  alpha_array = M_s * V_c * B_array / (k_B * T_p_kelvin)

  coth_alpha = 1./(TANH(alpha_array))
  slope_array = (1. + alpha_array^2 - alpha_array^2* coth_alpha^2) / (alpha_array*coth_alpha-1.)
  print, 'slope _array: ', slope_array[0]
  ;1.01 mPa.s
  print, 'eta = ', eta[101]
  ;  iplot, B_array*1000., tao_B[101]*slope_array, title='1.0 mPa.s', xtitle='Field Amplitude [mT]', ytitle='Brownian Relaxation Time [us]'
  print, 'Brownian(1.01 mPa.s, 4.5 mT) us: ', tao_B[101]*slope_array[44]
  ;  Brownian_time = tao_B[101]*slope_array[0]
  Brownian_time_eta_B_array = DBLARR(N_ELEMENTS(eta), N_ELEMENTS(B_array))
  FOR i=0, N_ELEMENTS(eta)-1 DO BEGIN
    FOR j=0, N_ELEMENTS(B_array)-1 DO BEGIN
      Brownian_time_eta_B_array[i,j] = tao_B[i]*slope_array[j]
    ENDFOR
  ENDFOR


  beta_c = V_c / (k_B * T_p_kelvin); nm3/J
  ;  K_c_nm = 20.e-6; J/nm
  K_c_m = 2350. * k_B / V_c; J/m
  sigma_c = 2350./T_p_kelvin
  ;h_array = alpha_array/(2.*sigma_c)
  h_array = M_s * B_array / (2.*K_c_m)
  alpha_prim = 0.1
  gamma_c = 1.75e11; rad/sT
  tao_N0 = beta_c * (1.+alpha_prim^2) * (M_s) / (2.*gamma_c*alpha_prim) * 1.e6 ;us
  print, 'tao Neel 0 = ', tao_N0
  ele_temp = (1.-h_array)*exp(-sigma_c*(1.-h_array)^2) + (1.+h_array)*exp(-sigma_c*(1.+h_array)^2)
  Neel_slope_array = SQRT(!PI) * sigma_c^(-1.5) / (1.-h_array^2) / ele_temp

  ;1.01 mPa.s
  print, 'eta = ', eta[101]
  ;  iplot, B_array*1000., tao_N0*Neel_slope_array, title='1.0 mPa.s', xtitle='Field Amplitude [mT]', ytitle='Neel Relaxation Time [us]'
  print, 'Neel(4.5 mT) us: ', tao_N0*Neel_slope_array[44]
  ;  Neel_time = tao_N0*Neel_slope_array[0]
  Neel_time_B_array = tao_N0*Neel_slope_array

  ;  window, 2
  ;  plot, B_array*1000., tao_B[101]*slope_array, BACKGROUND = 'FFFFFF'x, COLOR = 0, yrange=[0, 100.], xtitle='Field Amplitude [mT]', ytitle='Relaxation time [us]', title='Neel & Brownian vs Field'
  ;  oplot, B_array*1000., tao_N0*Neel_slope_array, LINESTYLE=2, color=3500
END

;particle_size = 24.4924;PBS, [25.2579], 1% Glycerol
;T_p = 20.;36.5; degree
;@GJ, 2023/2/12, calculating Brownian relaxation based on pulsed excitation
;@GJ, 2023/3/10, relaxation time with field changes
PRO relaxation_time_Brownian_Neel, D_c, D_h, eta, B_array, Neel_time_B_array, Brownian_time_eta_B_array
  eta = FINDGEN(401)/100.*1.e-3 ; Pa.s
  IF N_ELEMENTS(D_h) EQ 0 THEN D_h = 60.;
  print, 'D_h: ', D_h, ' nm'
  D_h *= 1.e-9 ; m  
  V_h = 1./6. * !PI * D_h^3
  k_B = 1.380469e-23; J/K
  T_p = 25.; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  
  tao_B = 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
  print, 'tao_B: ', tao_B, ' us'
;  iplot, eta, tao_B, xtitle='Viscosity [mPa.s]', ytitle='Brownian Relaxation Time [us]'
  
  M_s = 474000.; J/m3 T (Ms = 474000 J/m3 T)
  IF N_ELEMENTS(D_c) EQ 0 THEN D_c = 16; nm
  D_c *= 1.e-9 ;m
  ;D_c = 5.e-9 ;m
  V_c =  1./6. * !PI * D_c^3
  
  B_array = (FINDGEN(101.)/10.+0.1)*1.e-3 ;T
  alpha_array = M_s * V_c * B_array / (k_B * T_p_kelvin)
   
  coth_alpha = 1./(TANH(alpha_array))
  slope_array = (1. + alpha_array^2 - alpha_array^2* coth_alpha^2) / (alpha_array*coth_alpha-1.)
  print, 'slope _array: ', slope_array[0]   
  ;1.01 mPa.s
  print, 'eta = ', eta[101]
;  iplot, B_array*1000., tao_B[101]*slope_array, title='1.0 mPa.s', xtitle='Field Amplitude [mT]', ytitle='Brownian Relaxation Time [us]'
  print, 'Brownian(1.01 mPa.s, 4.5 mT) us: ', tao_B[101]*slope_array[44]
;  Brownian_time = tao_B[101]*slope_array[0]
  Brownian_time_eta_B_array = DBLARR(N_ELEMENTS(eta), N_ELEMENTS(B_array))
  FOR i=0, N_ELEMENTS(eta)-1 DO BEGIN
    FOR j=0, N_ELEMENTS(B_array)-1 DO BEGIN
      Brownian_time_eta_B_array[i,j] = tao_B[i]*slope_array[j]
    ENDFOR
  ENDFOR
  
  
  beta_c = V_c / (k_B * T_p_kelvin); nm3/J
;  K_c_nm = 20.e-6; J/nm
  K_c_m = 2.e4; J/m
  sigma_c = K_c_m * beta_c
  ;h_array = alpha_array/(2.*sigma_c)
  h_array = M_s * B_array / (2.*K_c_m)
  alpha_prim = 0.1
  gamma_c = 1.75e11; rad/sT
  tao_N0 = beta_c * (1.+alpha_prim^2) * (M_s) / (2.*gamma_c*alpha_prim) * 1.e6 ;us
  print, 'tao Neel 0 = ', tao_N0
  ele_temp = (1.-h_array)*exp(-sigma_c*(1.-h_array)^2) + (1.+h_array)*exp(-sigma_c*(1.+h_array)^2)
  Neel_slope_array = SQRT(!PI) * sigma_c^(-1.5) / (1.-h_array^2) / ele_temp
  
  ;1.01 mPa.s
  print, 'eta = ', eta[101]
;  iplot, B_array*1000., tao_N0*Neel_slope_array, title='1.0 mPa.s', xtitle='Field Amplitude [mT]', ytitle='Neel Relaxation Time [us]'
  print, 'Neel(4.5 mT) us: ', tao_N0*Neel_slope_array[44]
;  Neel_time = tao_N0*Neel_slope_array[0]
  Neel_time_B_array = tao_N0*Neel_slope_array
  
;  window, 2
;  plot, B_array*1000., tao_B[101]*slope_array, BACKGROUND = 'FFFFFF'x, COLOR = 0, yrange=[0, 100.], xtitle='Field Amplitude [mT]', ytitle='Relaxation time [us]', title='Neel & Brownian vs Field'
;  oplot, B_array*1000., tao_N0*Neel_slope_array, LINESTYLE=2, color=3500
END

;@GJ, 2023/4/10, calculating Neel and Brownian with different Dc and Dh
;@GJ, 2023/4/27, update the code
PRO relaxation_time_Dc_full_field_approximation

D_c_array = FINDGEN(251)/10.+5.
Neel_array = D_c_array*0.
Brownian_array = D_c_array*0.

FOR i=0, N_ELEMENTS(D_c_array)-1 DO BEGIN
  D_h = 50.
  relaxation_time_Brownian_Neel, D_c_array[i], D_h, eta, B_array, Neel_time_B_array, Brownian_time_eta_B_array
  ;1.01 mPas, 4.5 mT
  Neel_array[i] = Neel_time_B_array[44];4.5 mT
  Brownian_array[i] = Brownian_time_eta_B_array[101, 44];4.5 mT
ENDFOR

iplot, D_c_array, Neel_array, /ylog, xtitle='D_c [nm]', ytitle='Relaxation Time [us]', title='Neel Relaxation Time vs D_c (4.5 mT)'
iplot, D_c_array, Brownian_array, xtitle='D_c [nm]', ytitle='Relaxation Time [us]', title='Brownian Relaxation Time vs D_c (1.01 mPa.s, 4.5 mT)'

window, 1
plot, D_c_array, Neel_array, /ylog, BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='D_c [nm]', ytitle='Relaxation time [us]', title='Neel & Brownian vs D_c (1.01 mPa.s, 4.5 mT)'
oplot, D_c_array, Brownian_array, LINESTYLE=2, color=3500

;plot the 16.5 nm, 50 nm at different field
D_h = 50.
relaxation_time_Brownian_Neel, 16.5, D_h, eta, B_array, Neel_time_B_array, Brownian_time_eta_B_array
window, 2
plot, B_array*1000., Brownian_time_eta_B_array[101,*], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Field Amplitude [mT]', ytitle='Relaxation time [us]', title='Neel & Brownian vs Field Amplitude (1.01 mPa.s, 16.5nm+50nm)'
oplot, B_array*1000., Neel_time_B_array, LINESTYLE=2, color=3500

;plot the  16.5 nm, 50 nm at different viscosity at 4.5 mT
window, 3
plot, eta*1000., Brownian_time_eta_B_array[*,9], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Viscosity [mPa.s]', ytitle='Relaxation time [us]', title='Brownian vs Field Amplitude (1, 2, 4.5 mT, 16.5nm+50nm)'
oplot, eta*1000., Brownian_time_eta_B_array[*,19], LINESTYLE=2, color=3500
oplot, eta*1000., Brownian_time_eta_B_array[*,44], LINESTYLE=2, color=3500

slope_Browian_eta = DBLARR(N_ELEMENTS(B_array))
sensitivity_Browian_eta = DBLARR(N_ELEMENTS(B_array))
FOR i=0, N_ELEMENTS(B_array)-1 DO BEGIN
  fit_eta_Brownian = LINFIT(eta*1000., REFORM(Brownian_time_eta_B_array[*,i])) ;us/mPa.s
  slope_Browian_eta[i] = fit_eta_Brownian[1];us/mPa.s
  sensitivity_Browian_eta[i] = fit_eta_Brownian[1]/Brownian_time_eta_B_array[91,i]*100.;1/mPa.s %
ENDFOR
;window, 4
;plot, B_array, slope_Browian_eta
iplot, B_array*1000., slope_Browian_eta, xtitle='Field Amplitude [mT]', ytitle='Slope [us/mPa.s]', title='Brownian Relaxation Time Slope'
iplot, B_array*1000., sensitivity_Browian_eta, xtitle='Field Amplitude [mT]', ytitle='Sensitivity [1/mPa.s %]', title='Brownian Relaxation Time Sensitivity [ref 0.91 mPa.s]'
END


;based on Medical Physics paper to interpolate relaxation time
;GJ, 2021/9/1
FUNCTION relaxation_time_Resovist, H_field, Frequency
  
  freq_array = [4.5, 9.3, 12.2, 25.0]
  H_array = [5.0, 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30.]
  relax_array = DBLARR(4, 11) * 0.
  relax_array[0, *] = [7.9, 6.8, 6.2, 5.8, 5.4, 5.1, 4.9, 4.6, 4.3, 4.2, 3.9]
  relax_array[1, *] = [3.9, 3.4, 3.2, 3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3]
  relax_array[2, *] = [2.5, 2.3, 2.3, 2.2, 2.0, 2.0, 1.9, 1.8, 1.7, 1.6, 1.6]
  relax_array[3, *] = [1.2, 1.3, 1.2, 1.2, 1.2, 1.1, 1.0, 1.0, 1.0, 1.0, 0.9]

  ;H_field = 14.4
  ;Frequency = 33.
  relax_H_array = DBLARR(4) * 0.
  relax_H_array[0] = INTERPOL(relax_array[0, *], H_array, H_field)
  relax_H_array[1] = INTERPOL(relax_array[1, *], H_array, H_field)
  relax_H_array[2] = INTERPOL(relax_array[2, *], H_array, H_field)
  relax_H_array[3] = INTERPOL(relax_array[3, *], H_array, H_field)
  
  relaxation_time = INTERPOL(relax_H_array, freq_array, Frequency)
  IF relaxation_time LT 0. THEN relaxation_time = 0.
;  window, 1
;  plot, freq_array, relax_H_array, XRANGE=[0, 30], YRANGE=[0, 10], XTITLE='Frequency [kHz]', YTITLE='Relaxation Time [us]', title = 'Resovist: '+STRING(relaxation_time,format='(f10.4)')+'uT ('+STRING(H_field,format='(f10.4)') + 'mT, '+STRING(Frequency,format='(f10.4)')+'kHz)'
;  FOR i = 0, 10 DO BEGIN
;    oplot, freq_array, REFORM(relax_array[*, i])
;  ENDFOR
;  oplot, freq_array, relax_H_array, PSYM=1
;  oplot, [Frequency], [relaxation_time], PSYM=2
;  
;  window, 2
;  plot, H_array, relax_array[0, *], XRANGE=[0, 40], YRANGE=[0, 10], XTITLE='Drive Field [mT]', YTITLE='Relaxation Time [us]', title = 'Resovist: '+STRING(relaxation_time,format='(f10.4)')+'uT ('+STRING(H_field,format='(f10.4)') + 'mT, '+STRING(Frequency,format='(f10.4)')+'kHz)
;  oplot, H_array, relax_array[1, *]
;  oplot, H_array, relax_array[2, *]
;  oplot, H_array, relax_array[3, *]
;  oplot, [H_field], [relax_H_array[0]], PSYM=1
;  oplot, [H_field], [relax_H_array[1]], PSYM=1
;  oplot, [H_field], [relax_H_array[2]], PSYM=1
;  oplot, [H_field], [relax_H_array[3]], PSYM=1
;  oplot, [H_field], [relaxation_time], PSYM=2
  
  RETURN, relaxation_time
END


;GJ & Zhe Wang, 2021/10/18, calculate relaxation time based on the Figure 38 in the patent
FUNCTION relaxation_pulsed_excitation, H_field, particle_size, plotYes=plotYes
  H_array = [-30,-20,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30]
  particle_size_array = [21.0, 25.0, 27.0, 32.0]
  relax_array = DBLARR(4, 30) * 0.
  relax_array[0, *] = [2.5,3.05,3.635,4.7,4.95,5.1,5.3,5.34,5.4,5.1,4.8,4.2,4.05,4,3.9,3.8,3.7,3.6,3.6,3.55,3.5,3.4,3.3,3.2,3.1,3.0,2.9,2.8,2.75,2.5]
  relax_array[1, *] = [2.5,3.05,3.2,4.1,4.6,5.1,5.9,7.2,8.4,10.1,12.5,14.5,15.6,14.6,13.4,11.6,10.5,9.3,8.230,7.4,6.8,6.2,5.9,5.6,5.3,5.1,4.9,4.7,4.39,3]
  relax_array[2, *] = [2.5,3.05,3.635,4.9,5.6,6.5,7.5,9,11.8,15.4,18.4,21.8,23.6,22.9,21.3,19.0,16.9,15.1,13.15,11.8,10.5,9.5,8.5,7.7,7.1,6.5,6.1,5.8,5.43,3]
  relax_array[3, *] = [2.5,3.05,3.2,4.1,4.3,4.6,5.4,6.2,7.8,9.8,18.5,24.8,28.4,31,28.9,26.2,23.5,21.2,18.8,17.0,15.3,14.1,13.0,11.9,11.0,10.1,9.3,8.5,7.88,3.76]
 
  relax_H_array = DBLARR(4) * 0.
  relax_H_array[0] = INTERPOL(relax_array[0, *], H_array, H_field)
  relax_H_array[1] = INTERPOL(relax_array[1, *], H_array, H_field)
  relax_H_array[2] = INTERPOL(relax_array[2, *], H_array, H_field)
  relax_H_array[3] = INTERPOL(relax_array[3, *], H_array, H_field)

  relaxation_time = INTERPOL(relax_H_array, particle_size_array, particle_size)
  IF relaxation_time LT 0. THEN relaxation_time = 0.
  IF KEYWORD_SET(plotYes) THEN BEGIN
      window, 1
      plot, H_array, relax_array[0, *],XRANGE=[-30, 30], YRANGE=[0, 35], XTITLE='Particle size [nm]', YTITLE='Relaxation Time [us]', title = STRING(relaxation_time,format='(f10.4)')+'us ('+STRING(H_field,format='(f10.4)') + 'mT, '+STRING(particle_size,format='(f6.1)')+'nm)'
      oplot, H_array, relax_array[1, *]
      oplot, H_array, relax_array[2, *]
      oplot, H_array, relax_array[3, *]
      FOR i = 0, 3 DO oplot, [H_field], REFORM(relax_H_array[i]), PSYM=1
      oplot, [H_field], [relaxation_time], PSYM=2
  ENDIF

  RETURN, relaxation_time

END


;based on Medical Physics paper to interpolate relaxation time
;GJ, 2021/9/1
FUNCTION relaxation_time_UW33, H_field, Frequency

  freq_array = [4.5, 9.3, 12.2, 25.0]
  H_array = [5.0, 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30.]
  relax_array = DBLARR(4, 11) * 0.
  relax_array[0, *] = [3.5, 3.0, 2.8, 2.5, 2.3, 2.1, 1.9, 1.8, 1.7, 1.6, 1.5]
  relax_array[1, *] = [2.4, 1.7, 1.5, 1.4, 1.2, 1.1, 1.1, 1.1, 1.0, 1.0, 1.0]
  relax_array[2, *] = [1.4, 1.2, 1.0, 0.9, 0.8, 0.7, 0.7, 0.7, 0.6, 0.6, 0.6]
  relax_array[3, *] = [0.7, 0.6, 0.5, 0.5, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5]

  ;  H_field = 14.4
  ;  Frequency = 3.3
  relax_H_array = DBLARR(4) * 0.
  relax_H_array[0] = INTERPOL(relax_array[0, *], H_array, H_field)
  relax_H_array[1] = INTERPOL(relax_array[1, *], H_array, H_field)
  relax_H_array[2] = INTERPOL(relax_array[2, *], H_array, H_field)
  relax_H_array[3] = INTERPOL(relax_array[3, *], H_array, H_field)

  relaxation_time = INTERPOL(relax_H_array, freq_array, Frequency)
  IF relaxation_time LT 0. THEN relaxation_time = 0.
  window, 1
  plot, freq_array, relax_H_array, XRANGE=[0, 30], YRANGE=[0, 5], XTITLE='Frequency [kHz]', YTITLE='Relaxation Time [us]', title = 'UW33: '+STRING(relaxation_time,format='(f10.4)')+'uT ('+STRING(H_field,format='(f10.4)') + 'mT, '+STRING(Frequency,format='(f10.4)')+'kHz)'
  FOR i = 0, 10 DO BEGIN
    oplot, freq_array, REFORM(relax_array[*, i])
  ENDFOR
  oplot, freq_array, relax_H_array, PSYM=1
  oplot, [Frequency], [relaxation_time], PSYM=2
  
  window, 2
  plot, H_array, relax_array[0, *], XRANGE=[0, 40], YRANGE=[0, 5], XTITLE='Drive Field [mT]', YTITLE='Relaxation Time [us]', title = 'UW33: '+STRING(relaxation_time,format='(f10.4)')+'uT ('+STRING(H_field,format='(f10.4)') + 'mT, '+STRING(Frequency,format='(f10.4)')+'kHz)
  oplot, H_array, relax_array[1, *]
  oplot, H_array, relax_array[2, *]
  oplot, H_array, relax_array[3, *]
  oplot, [H_field], [relax_H_array[0]], PSYM=1
  oplot, [H_field], [relax_H_array[1]], PSYM=1
  oplot, [H_field], [relax_H_array[2]], PSYM=1
  oplot, [H_field], [relax_H_array[3]], PSYM=1
  oplot, [H_field], [relaxation_time], PSYM=2

  RETURN, relaxation_time
END

;GJ, 2021/8/26, signal with non-adiabatic theory
;signal_nat = non_adiabatic(T, signal, relaxtion_time)
FUNCTION non_adiabatic, T, signal, relaxation_time
  IF relaxation_time GT 0. THEN BEGIN
  
    Nsample = N_ELEMENTS(T)
    distfun = DINDGEN(Nsample)   ;  1D Euclidean dist function.
    distfun = distfun < (Nsample - distfun)
    distfun[0] = 1d-4
    tao = relaxation_time[0]
    r = (1./tao) * exp(-distfun * (T[1]-T[0])/tao)
    filter = distfun * 0.
    filter[Nsample/2:*] = r[0:N_ELEMENTS(filter[Nsample/2:*])-1]
    filter=REVERSE(filter, /OVERWRITE)
    ;GJ, 2021/11/8, fix the bug of signal[*] is not zero because signal is more than a cycle. signal[0:999] is a cycle and mean is zero
;    filteredsignal_temp = CONVOL([DBLARR(Nsample)*0+MEAN(signal[0:4]), signal, DBLARR(Nsample)*0+MEAN(signal[0:Nsample-1])], filter, /NORMALIZE)
;    filteredsignal = filteredsignal_temp[Nsample:(2*Nsample-1)]
    
    filteredsignal_temp = CONVOL([DBLARR(Nsample)*0+MEAN(signal[0:4]), signal], filter, /NORMALIZE, /CENTER, /EDGE_ZERO)
    filteredsignal = filteredsignal_temp[Nsample:*]
    ;  oplot, T, filteredsignal, LINESTYLE=3
  ENDIF ELSE BEGIN
    filteredsignal = signal
  ENDELSE

  RETURN, filteredsignal;signal_nat

END

;GJ, 2021/8/26, signal with non-adiabatic theory
;signal_nat = non_adiabatic(T, signal, relaxtion_time)
FUNCTION dec_non_adiabatic, T, signal_ref, signal_tao
  maxind_signal_ref = 249; MAX(signal_ref, maxind_signal_ref)
  signal_ref_max = MAX(signal_ref, maxind_signal_ref)
  signal_tao_max = MAX(signal_tao, maxind_signal_tao)
  signal_max_diff = signal_ref_max - signal_tao_max
  signal_max_diff_ratio = signal_max_diff/signal_ref_max

  delay_time = maxind_signal_tao - maxind_signal_ref
  
  ;move the peak back to the center
  signal_temp = CONGRID(signal_tao[0:maxind_signal_tao], 250)
  
  ;interpolate to get the peak signal
  ;signal_interpo = INTERPOL(signal_temp[220:239], T[220:239], T[240:249])
  ;signal_temp[240:249] = signal_interpo
  
  ;based on AUC is constant, rescale the signal
  AUC_signal = TOTAL(signal_tao[0:499], /NAN)/2.
  AUC_signal_temp = TOTAL(signal_temp, /NAN)
  signal_scaled = signal_temp * AUC_signal / AUC_signal_temp
  
  ;based on the relationship between delay time and signal peak difference to do calibration
;  coeff = [0.65716375, -0.025619823]
;  recov_signal_max = signal_tao_max / (1.0- (coeff[0] + coeff[1]*delay_time))
;  scaled_ratio = recov_signal_max - signal_scaled[249]
;  FOR i = 0, 9 DO signal_scaled[240+i] = signal_scaled[240+i] + (i+1)/10. * scaled_ratio
;
;  ;rescale back to the AUC
;  REPEAT BEGIN      
;    AUC_temp = TOTAL(signal_scaled, /NAN)
;    ratio = (AUC_temp - AUC_signal)/AUC_temp
;  ;  signal_scaled[0:249] = signal_scaled[0:249]*(1.0 - ratio)
;    N_pixels = FLOOR(250 * ratio)
;    IF N_pixels GT 2 THEN BEGIN
;      signal_temp = CONGRID(signal_scaled, 250-N_pixels)
;      signal_scaled = DBLARR(250)*0.
;      signal_scaled[N_pixels:249] = signal_temp
;    ENDIF
;    AUC_temp = TOTAL(signal_scaled, /NAN)
;    ratio = (AUC_temp - AUC_signal)/AUC_temp
;  ENDREP UNTIL ratio LT 0.01
;
;  
;  
  ;make a complete signal
  dec_signal = [signal_scaled, REVERSE(signal_scaled), -signal_scaled, -REVERSE(signal_scaled), signal_scaled[0:198]]

  RETURN, dec_signal;signal_nat

END

;GJ, 2021/6/30, simulation is much faster
;GJ, 2021/7/11, 1d V-shape field encoding and deconding + 1d FBP
;GJ, 2021/7/28, file is moved to default folder
;GJ, 2021/8/25, using MRI graidents to do signal encoding
;GJ, 2021/8/31, adding gradient values based on angles (0 to 180 degrees)
PRO aimis3d_MPT_2d_Gradient_field

;    filename = '.\MPI_images\mp.bmp'
;    image_temp2 = read_bmp(filename)
;    image = MAX(image_temp2) - image_temp2

;    ;get the sample image
;    filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
;    image_temp2 = read_image(filename2)
;    image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256))

  ;  filename2 = 'C:\D_drive\MPI_Tianjie\MPT_english.bmp'
  ;;  filename2 = 'C:\D_drive\MPI_Tianjie\MagneticPT_english.bmp'
  ;  image_temp2 = read_bmp(filename2)
  ;  image = MAX(image_temp2) - image_temp2

  filename2='F:\MPI_Tianjie\MIP\000.tif'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256), MAX=130)
;      image = image * 0.
;      image[100:102, 127:129] = 256
;      image[130:132, 127:129] = 256

  ;setup the 1d image encoding and decoding kernal
  iimage,  image, title='Original'
  Matrix_size = N_ELEMENTS(image[0,*]); causing invert calcuation error
  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  x_size = N_ELEMENTS(image[*,0])
  ;set up the encoding magnetic fields
  H0 = 0.000;5;15 ;T
  G = 0.075; T/m
  delta_G = G*2./Matrix_size
  FOV = 0.2; m
  x_res = FOV/x_size
  G_array = (FINDGEN(Matrix_size)-Matrix_size/2)*delta_G
  ;G_array = DBLARR(Matrix_size)*0. - G
  Hx = ((FINDGEN(x_size*2)-x_size)*G*x_res + H0) * 1000. ;mT
  x_axis = (FINDGEN(x_size*2)-x_size)*x_res * 1000. ;mm
  iplot, x_axis, Hx
  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = (((FINDGEN(x_size)-x_size)+i)*G/3.*x_res+H0) * 1000. ; mT
 
   
  ;define the signal kernal
  signalKernal_2d_array = H_2d_array*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      ;H_flat_portion = 0.0
      H_flat_portion = 0.9
      ;a = signal_vs_H(H_array[i],H_flat_portion,0)
      ;      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j], 0))[0, 3]; for 3rd harmonic; [0] for Amplitude
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[1, 1] ;3rd harmonic without relaxation
    ENDFOR
  ENDFOR
  
  ;calculate the invert matrix of the signal Kernal
  ;  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  Gradient_invert = INVERT(signalKernal_2d_array, /DOUBLE)
  result = Gradient_invert # signalKernal_2d_array
  iimage, result, title='Unit invert matrix *'
  
  ;GJ, 2021/9/5, do calculation based ART method
  ;rec_image_1d_art = ART_cal(signal_3rdHar_1d, signalKernal_2d_array)
  ;GJ, 2021/9/25, find solution using Liang Xiaofeng's code
  ;A=[[1,1,1,1],[1,2,-1,4],[2,-3,-1,-5],[3,1,2,11]]
  ;b=[[5],[-2],[-2],[0]]
  ;rec_image_1d = art_func(A,b)
;  image = DBLARR(256) * 0.
;  image[100:102, 127:129] = 256
;  image[130:132, 127:129] = 256
;  signal_1d = signalKernal_2d_array # REFORM(image[*, x_size/2])
;  rec_image_1d = art_func(signalKernal_2d_array,signal_1d)
  ;end of test
  
  

  ;;Calculate and display the Radon transform:
  Radon_result_ori = RADON(image, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)
  Radon_result_temp = NOISE_PICK(Radon_result_ori, 0.1, ITER=6)
  ; Display the images side by side:
;  IIMAGE, Radon_result_ori, VIEW_GRID=[2,1], VIEW_TITLE='Original', $
;    DIMENSIONS=[850, 550], WINDOW_TITLE='NOISE_PICK Example', $
;    /NO_SAVEPROMPT
;  IIMAGE, Radon_result_temp, /VIEW_NEXT, VIEW_TITLE='Noisy'
  
  new_Radon_result = CONGRID(Radon_result_temp, 180, 256)
  new_rho = CONGRID(rho, 256)
  new_theta = CONGRID(theta, 180)
  backproject0 = RADON(new_Radon_result, /BACKPROJECT, RHO=new_rho, THETA=new_theta)
;  iimage, backproject0, title='BP without filtered'
  ;calculate gradient
  theta_deg = new_theta/(!PI)*180.
  N_theta = N_ELEMENTS(new_theta)
  phi = 90.
  
  ;define and calculate gradients
  G_spherical_array = DBLARR(N_theta, Matrix_size, 3) * 0. 
  Gxyz_array = DBLARR(N_theta, Matrix_size, 3) * 0.
  FOR i = 0, N_theta-1 DO BEGIN
    FOR j = 0, Matrix_size-1 DO BEGIN
      G_spherical_array[i, j, 0] = G_array[j]
      G_spherical_array[i, j, 1] = theta_deg[i]
      G_spherical_array[i, j, 2] = phi
      G_xyz = Gradient_spherical_to_rectangular(G_array[j], theta_deg[i], phi)
      Gxyz_array[i, j, 0] = G_xyz[0]
      Gxyz_array[i, j, 1] = G_xyz[1]
      Gxyz_array[i, j, 2] = G_xyz[2]
    ENDFOR
  ENDFOR
  iimage, G_spherical_array[*,*,0], title='G'
  iimage, G_spherical_array[*,*,1], title='theta'
  iimage, G_spherical_array[*,*,2], title='phi'
  iimage, Gxyz_array[*,*,0], title='Gx'
  iimage, Gxyz_array[*,*,1], title='Gy'
  iimage, Gxyz_array[*,*,2], title='Gz'
  
;  ;calculate H field
;  Hxyz_array = DBLARR(N_theta, Matrix_size, Matrix_size, Matrix_size, Matrix_size) * 0.
;  FOR i = 0, N_theta-1 DO BEGIN
;    FOR j = 0, Matrix_size-1 DO BEGIN
;      H_volume = H_volume_from_G(G_array[j], theta_deg[i], phi, [Matrix_size, Matrix_size, Matrix_size], FOV)
;      Hxyz_array[i, j, *, *, *] = H_volume
;    ENDFOR
;  ENDFOR
  
  rec_result_filter = new_Radon_result*0.
  rec_tao_result_filter = new_Radon_result*0.
  new_Radon_result_filter = new_Radon_result*0.
  signal_tao_Amplitude_1d = DBLARR(Matrix_size)*0.
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  signal_AUC_1d = DBLARR(Matrix_size)*0.
  signal_3rdHar_1d = DBLARR(Matrix_size)*0.
  signal_tao_3rdHar_1d = DBLARR(Matrix_size)*0.
  signal_tao_AUC_1d = DBLARR(Matrix_size)*0.
  feature_diff_signal_tao = DBLARR(Matrix_size*Matrix_size,3)*0.
  FOR i=0, N_ELEMENTS(new_theta)-1 DO BEGIN
    image_1d = CONGRID(REFORM(new_Radon_result[i, *]), Matrix_size)
    
    FOR j = 0, Matrix_size-1 DO BEGIN
      FOR k = 0, Matrix_size-1 DO BEGIN
        H_flat_portion = 0.9
        system_temp = signal_vs_H(H_2d_array[j, k],H_flat_portion,1)
        ;print, j, k, '=', where(~finite(system_temp[2,*]))
        signal_temp = system_temp*image_1d[k]/Matrix_size
        IF k EQ 0 THEN BEGIN
          ;first time, create one
          signal_average = signal_temp
        ENDIF ELSE BEGIN
          ;the rest, just add together
          signal_average = signal_average + signal_temp
        ENDELSE
      ENDFOR
      ;the T was corrected
      signal_average[0, *] = system_temp[0, *]
      

;      ;GJ, 2021/9/3, do relaxation time correction
;      feature_diff_signal_tao[(i*Matrix_size+j), *] = feature_diff_signal_tao(signal_average)
;      IF TOTAL(signal_average[2,*]) GT 3 THEN BEGIN
;        deconv_signal_tao = signal_deconv_correction(signal_average)
;        signal_average[2, *] = deconv_signal_tao
;      ENDIF

      
;      ;plot the signal and signal with relaxation
;      window, 3
;      plot, signal_average[0,200:300], signal_average[1,200:300], title='signal measured', XTITLE='time [us]', YTITLE='signal'
;      oplot, signal_average[0,200:300], signal_average[2,200:300], PSYM=2

      
      ;evaluate the features after correction
      H_flat_portion = 0.9
      feature_temp = (signal_vs_H(10, H_flat_portion, 0, signal_average))
      signal_Amplitude_1d[j] = feature_temp[0, 0]
      signal_AUC_1d[j] = feature_temp[0, 1]
      signal_3rdHar_1d[j] = feature_temp[0, 3]
      signal_tao_Amplitude_1d[j] = feature_temp[1, 0]
      signal_tao_AUC_1d[j] = feature_temp[1, 1]
      ;do correction based on the calibration curve
      signal_tao_3rdHar_1d[j] = feature_temp[1, 3]
      
      ;zero the signal average for next cycle
      signal_average = signal_average * 0.
    ENDFOR
;    signal_Amplitude_1d = signalKernal_2d_array # image_1d
    rec_image_1d = Gradient_invert # signal_3rdHar_1d
    ;GJ, 2021/9/5, do calculation based ART method
    ;rec_image_1d_art = ART_cal(signal_3rdHar_1d, signalKernal_2d_array)
    ;GJ, 2021/9/25, find solution using Liang Xiaofeng's code
    ;A=[[1,1,1,1],[1,2,-1,4],[2,-3,-1,-5],[3,1,2,11]]
    ;b=[[5],[-2],[-2],[0]]
    ;rec_image_1d = art_func(A,b)
    ;rec_image_1d = art_func(signalKernal_2d_array,signal_3rdHar_1d)

    rec_remp = CONGRID(REFORM(rec_image_1d), N_ELEMENTS(new_rho))
    ;do filter
    rec_result_filter[i, *] = Filter_FFT(rec_remp)
    
    ;with relaxation filter
    ;signal_tao_3rdHar_1d = -0.75014898 + 1.0546962*signal_tao_3rdHar_1d
    rec_tao_image_1d = Gradient_invert # signal_tao_3rdHar_1d
    rec_tao_remp = CONGRID(REFORM(rec_tao_image_1d), N_ELEMENTS(new_rho))
    ;do filter
    rec_tao_result_filter[i, *] = Filter_FFT(rec_tao_remp)

    
    new_Radon_result_filter[i, *] = Filter_FFT(REFORM(new_Radon_result[i, *]))
    print, 'i = ', i


    window, 4
    image_1d_ori = CONGRID(REFORM(new_Radon_result[i, *]), 256)
    plot, image_1d_ori, BACKGROUND = 'FFFFFF'x, COLOR = 0
    oplot, rec_image_1d/MEAN(rec_image_1d)*MEAN(image_1d_ori), PSYM=2, color=3500
;    oplot, rec_tao_image_1d/MEAN(ABS(rec_tao_image_1d))*MEAN(image_1d_ori), PSYM=4, color=3500
;;    wait, 0.1
;
;    window, 5
;    plot, signal_AUC_1d, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='AUC'
;    oplot, signal_tao_AUC_1d, PSYM=2, color=3500
;    
;    window, 6
;    plot, signal_Amplitude_1d, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='Peak Amplitude'
;    oplot, signal_tao_Amplitude_1d, PSYM=2, color=3500
;    
;    window, 7
;    plot, signal_3rdHar_1d, BACKGROUND = 'FFFFFF'x, COLOR = 0, title='3rd Har'
;    oplot, signal_tao_3rdHar_1d, color=3500;, PSYM=2
    ;coeff = poly_fit(signal_tao_3rdHar_1d, signal_3rdHar_1d, 1, yfit=yfit)
    ;print, coeff
    ;-0.75014898
    ;1.0546962
  ENDFOR
 
;  window, 8
;  plot, feature_diff_signal_tao[*, 1], feature_diff_signal_tao[*, 0], xtitle='delay time', ytitle='signal peak diff'
  
  
  backproject1 = RADON(new_Radon_result_filter, /BACKPROJECT, RHO=new_rho, THETA=new_theta)
  iimage, backproject1, title='FBP'
  


  backproject2 = RADON(rec_result_filter, /BACKPROJECT, RHO=new_rho, THETA=new_theta)
  iimage, backproject2, title='MPT using MRI gradient encoding'
  
  backproject3 = RADON(rec_tao_result_filter, /BACKPROJECT, RHO=new_rho, THETA=new_theta)
  iimage, backproject3, title='MPT using MRI gradient encoding'
  
  ;GJ, 2021/8/26
  radon_doc, image, backproject2
  
  ;GJ, 2021/9/2
;  radon_doc, image, backproject3
END

FUNCTION MPT_1d_Transform, image_1d

  Matrix_size = N_ELEMENTS(image_1d); causing invert calcuation error
  
  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  H_gradient = FINDGEN(Matrix_size) ;0 to 255
  ratios = FINDGEN(Matrix_size)/(Matrix_size-1)
  
  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = FINDGEN(Matrix_size)*ratios[i] + 200.
  
  signalKernal_2d_array = H_2d_array*0.
  
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  signal_AUC_1d = DBLARR(Matrix_size)*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signal = (signal_vs_H(H_2d_array[i, j], H_flat_portion, 0))*image_1d[j]
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j], H_flat_portion, 0))[0, 0] ;Amplitude
      signal_Amplitude_1d[i] = signal_Amplitude_1d[i] + signal[0, 0]
      signal_AUC_1d[i] = signal_AUC_1d[i] + signal[0, 1]
    ENDFOR
  ENDFOR
  
  ;;reconstruction
  ;signalKernal_2d_array[Matrix_size/2-1, *] = signalKernal_2d_array[0, *]
  ;signalKernal_2d_array[Matrix_size/2, *] = signalKernal_2d_array[0, *]
  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
;  meas = signalKernal_2d_array # image_1d
  test = Gradient_invert # signal_Amplitude_1d
;  result = Gradient_invert # (signalKernal_2d_array)
  
  plot, test, color=500
  oplot, image_1d
  ;iimage, test
  RETURN, test
END

FUNCTION Filter_FFT_temp, image_1d
  filename2 = 'C:\D_drive\MPI_Tianjie\MagneticPT_english.bmp'
  image_temp2 = read_bmp(filename2)
  image = MAX(image_temp2) - image_temp2
  ;;Calculate and display the Radon transform:
  Radon_result = RADON(image, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  N = N_ELEMENTS(rho)
  T = 1.
  X = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
  else $
    freq = [0.0, X, -(N/2 + 1) + X]/(N*T)
    
  FOR i=0, N_ELEMENTS(theta)-1 DO BEGIN
    image_1d = REFORM(Radon_result[i, *])
    P_FFT = FFT(image_1d)
    filtered_P_FFT = P_FFT * (ABS(freq)) / MAX((ABS(freq)))
    ;filtered_P_FFT = P_FFT
    ;filtered_P_FFT[1:*] = P_FFT[1:*] * (ABS(freq))[1:*]
    plot, freq, ABS(FFT(image_1d))
    oplot, freq, ABS(filtered_P_FFT), color=33500
    test = FFT(filtered_P_FFT, /INVERSE)
  ENDFOR
  RETURN, ABS(test)
END

;GJ, 2021/7/25
;Debye relaxation kernel
FUNCTION Debye_Kernel, tao, kernel_time
  ;tao1 = 1.0 ;us
  ;tao2 = 8.0 ;us
  ;tao = tao1
  ;33KHz,
  ;N = 33.* 1000. / 16.5; time points
  N = N_ELEMENTS(kernel_time)
  kernel = DBLARR(N)*0.
  ;
  ;
  ;采样频率sampling frequency 16 MHz
  d_t = kernel_time[1]-kernel_time[0] ;us
  time = FINDGEN(N) * d_t  ; us
  
  ;r1 = (1./tao1) * exp(-time/tao1)
  ;r2 = (1./tao2) * exp(-time/tao2)
  r = (1./tao) * exp(-time/tao)
  kernel[0:N/2] = REVERSE(r[0:N/2])
  
  ;
  ;; Create a sampled signal with random noise
  ;signal = SIN((FINDGEN(1000)/35.0))
  ;signal = (signal GT 0) * signal
  ;; Convolve the filter and signal using block convolution
  ;kernel=r[0:100]
  ;;result = BLK_CON(filter, signal)
  ;result = CONVOL(signal, kernel)
  
  RETURN, kernel
END

;GJ, 2021/7/25
;calculate the effect of relaxation to signal
;2021/7/26, deconvolution using FFT
PRO relaxation_Signal

  ;pre-defined frequency: 33 KHz
  ;sampling frequency 16.5 MHz
  ;sampling time interval 0.06 us (10-6 sec)
  ;X = 2*!PI/500. * FINDGEN(600)
  
  T = 1./16.5; us 10-6 sec
  X = FINDGEN(600) * T
  f = 0.033; MHz

  H = -COS(2*!PI* X * f)*10.
  H = H * 3.
  M = 1.0*(1./TANH(H) - 1./H)
  Nelem = N_ELEMENTS(M)
  signal = M[1:*]-M[0:Nelem-1]
  time = X[1:*]
  
  tao=0.5
  kernel = Debye_Kernel(tao, time[0:200])
  window, 1
  plot, kernel, BACKGROUND = 'FFFFFF'x, COLOR = 0
  
  signal_na_temp = CONVOL([DBLARR(300)*0, signal], kernel, /NORMALIZE);, CENTER=0, /NORMALIZE)
  signal_na = signal_na_temp[300:*]
  window, 2
  plot, time, signal, BACKGROUND = 'FFFFFF'x, COLOR = 0
  oplot, time, signal_na, LINESTYLE=0, COLOR = '0000FF'x
  
  ;GJ, 2021/7/26
  ;FFT-based convolution
  signal_FFT = FFT(signal, /CENTER)
  kernel_FFT = FFT(Debye_Kernel(tao, time), /CENTER)
;  signal_kernel_FFT = ABS(signal_FFT) * ABS(kernel_FFT)
  signal_kernel_FFT = signal_FFT * kernel_FFT
  signal_na_FFT = REAL_PART(FFT(signal_kernel_FFT, /INVERSE, /CENTER))
  oplot, time, signal_na_FFT, LINESTYLE=1, COLOR = '0000FF'x
  
  dec_kernel = Debye_Kernel(tao, time[0:200])
  dec_signal_na_temp = CONVOL([DBLARR(300)*0, signal_na], REVERSE(dec_kernel), /NORMALIZE);, CENTER=0, /NORMALIZE)
  dec_signal_na = dec_signal_na_temp[300:*]
  window, 3
  plot, time, signal, BACKGROUND = 'FFFFFF'x, COLOR = 0
  oplot, time, dec_signal_na, LINESTYLE=0, COLOR = '0000FF'x
  
;  ld_dec_signal_na = Filter_FFT(dec_signal_na)
;  oplot, time, ld_dec_signal_na, LINESTYLE=1, COLOR = '0000FF'x
  
END

;GJ, 2021/9/4 analyze signal_tao feature
;signal_tao_diff = feature_diff_signal_tao(signal_average)
FUNCTION feature_diff_signal_tao, signal_average
  T = REFORM(signal_average[0, *])
  signal_ref = REFORM(signal_average[1, *])
  signal_tao = REFORM(signal_average[2, *])
  signal_ref_max = MAX(signal_ref, maxind_signal_ref)
  signal_tao_max = MAX(signal_tao, maxind_signal_tao)
  
  signal_max_diff = signal_ref_max - signal_tao_max
  signal_max_diff_ratio = signal_max_diff/signal_ref_max
  
  delay_time = maxind_signal_tao - maxind_signal_ref
  
  signal_range = signal_tao[0:499]
  lower_signal_tao_sub = ABS(signal_range[0:maxind_signal_tao]-signal_tao_max/2)
  higher_signal_tao_sub = ABS(signal_range[maxind_signal_tao+1:*]-signal_tao_max/2)
  min_lower_signal_tao_sub = MIN(lower_signal_tao_sub, lower_index)
  min_higher_signal_tao_sub = MIN(higher_signal_tao_sub, higher_index)
  FWHM_temp = higher_index + (maxind_signal_tao - lower_index)
  FWHM = FWHM_temp * (T[1]-T[0]) ; in us
  
  FWHM_asym = (higher_index*1.0) / (maxind_signal_tao*1.0 - lower_index*1.0)
  
  feature_signal_tao = [signal_max_diff_ratio, delay_time, FWHM_asym]
  RETURN, feature_signal_tao
END


;GJ, 2021/9/3 do signal deconvolution for correction
;deconv_signal_tao = signal_deconv_correction(signal_average)
FUNCTION signal_deconv_correction, signal_average
  T = REFORM(signal_average[0, *])
  signal_ref = REFORM(signal_average[1, *])
  signal_tao = REFORM(signal_average[2, *])
  
;  window, 1
;  plot, T, signal_ref
;  oplot, T, signal_tao
;  relaxation_time = 1.0
;  relaxation_time = 
  deconv_signal_tao = dec_non_adiabatic(T, signal_ref, signal_tao)
  RETURN, deconv_signal_tao
END
;GJ, 2021/7/25
;Debye relaxation kernel
PRO Debye_Kernel_demo, tao, time
  tao1 = 1.0 ;us
  tao2 = 8.0 ;us
  tao = tao1
  ;33KHz,
  N = 33.* 1000. / 16.5; time points

  ;
  ;
  ;采样频率sampling frequency 16 MHz
  d_t = 1./16.5 ;us
  time = FINDGEN(N) * d_t  ; us

  r1 = (1./tao1) * exp(-time/tao1)
  r2 = (1./tao2) * exp(-time/tao2)
  r = (1./tao) * exp(-time/tao)


  ; Create a sampled signal with random noise
  signal = SIN((FINDGEN(1000)/35.0))
  signal = (signal GT 0) * signal
  ; Convolve the filter and signal using block convolution
  kernel=r[0:100]
  ;result = BLK_CON(filter, signal)
  result = CONVOL(signal, kernel)


;  RETURN, r
END

FUNCTION Filter_FFT, image_1d, filter_type

  N = N_ELEMENTS(image_1d)
  T = 0.1
  X = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  if (is_N_even) then $
    freqs = [0.0, X, N/2, -N/2 + X]/(N*T) $
  else $
    freqs = [0.0, X, -(N/2 + 1) + X]/(N*T)

  P_FFT = FFT(image_1d)
  filtered_P_FFT = P_FFT
  
  print, 'filter type: ', filter_type
  
  case filter_type of
    'Ramp': begin
      filter = abs(freqs)
    end
    'Shepp-Logan':begin
      filter = abs(freqs)/(1.+abs(freqs)^2)
    end
    'Hamming': begin
      filter = 0.54 + 0.46 * cos(2 *!PI * freqs / (N - 1))
      filter = filter * abs(freqs)
    end
;    'hann': begin
;      filter = 0.5 * (1 + cos(2 *!PI * freqs / (N - 1)))
;      filter = filter * abs(freqs)
;    end
    'Cosine': begin
      filter = cos(!PI * freqs / (N - 1))
      filter = filter * abs(freqs)
    end
    else: begin
      filter = abs(freqs)
    end
  endcase
  filtered_P_FFT[1:*] = P_FFT[1:*] * filter[1:*]
    
  ;filtered_P_FFT = P_FFT * (ABS(freq)) / MAX((ABS(freq)))
  abs_filtered_P_FFT = ABS(FFT(filtered_P_FFT, /INVERSE))
  return, abs_filtered_P_FFT
  
;  ;adjust the mean value, to the calibration, GJ, 2021/7/3
;  mean_1d = MEAN(image_1d)
;  stdev_1d = STDDEV(image_1d)
;  mean_filtered = MEAN(abs_filtered_P_FFT)
;  stdev_filtered = STDDEV(abs_filtered_P_FFT)
;  result = (abs_filtered_P_FFT-mean_filtered)/stdev_filtered*stdev_1d + mean_1d
;  
;  RETURN, result
END

;GJ, 2021/6/30, simulation is much faster
;GJ, 2021/7/4, simulation with single wire
PRO MPT_2d_single_wire

  ;  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20210627c\mp.bmp'
  ;  image_temp2 = read_bmp(filename)
  ;  image = MAX(image_temp2) - image_temp2

 ;   filename2 = 'C:\D_drive\MPI_Tianjie\MPT_english.bmp'
    filename2 = 'C:\D_drive\MPI_Tianjie\MagneticPT_english.bmp'
    image_temp2 = read_bmp(filename2)
    image = MAX(image_temp2) - image_temp2

;  filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
;  image_temp2 = read_image(filename2)
;  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256), MAX=130)

  ;setup the 1d image reconstruction kernal
  Matrix_size = N_ELEMENTS(image[0,*]); causing invert calcuation error
  Matrix_size_double = 512; causing invert calcuation error
  field_amp = DBLARR(Matrix_size_double,Matrix_size_double)*0. ; 100 rows, 201 columns
  x_field = DBLARR(Matrix_size_double,Matrix_size_double)*0. ; 100 rows, 201 columns
  y_field = DBLARR(Matrix_size_double,Matrix_size_double)*0. ; 100 rows, 201 columns
  ;;calculate the magnetic field
  coeff=1000.
  FOR i=0., Matrix_size_double-1. DO BEGIN
    FOR j=0., Matrix_size_double-1. DO BEGIN
      dist = SQRT((Matrix_size_double-i)^2 + (j+1)^2)
      field_amp[i,j] = coeff/dist
      x_field[i,j] = field_amp[i,j] * ((j+1) / dist)
      y_field[i,j] = field_amp[i,j] * ((Matrix_size_double-i) / dist)
    ENDFOR
  ENDFOR
  
  ;take the corner part of the field
  H_2d_array = x_field[0:Matrix_size-1, 0:Matrix_size-1]
  
  
;  H_2d_array = DBLARR(Matrix_size, Matrix_size)
;  ratios = FINDGEN(Matrix_size)/(Matrix_size-1)
;  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = FINDGEN(Matrix_size)*ratios[i] + 200.

  ;define the signal kernal
  signalKernal_2d_array = H_2d_array*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[0, 0]
    ENDFOR
  ENDFOR
  ;calculate the invert matrix of the signal Kernal
  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  result = Gradient_invert # (signalKernal_2d_array)
  iimage, result, title='inverse matrix multiplication'
  
  iimage,  image, title='Original'
  ;;Calculate and display the Radon transform:
  Radon_result = RADON(image, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)

  rec_result = Radon_result*0.
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  FOR i=0, N_ELEMENTS(theta)-1 DO BEGIN
    image_1d = CONGRID(REFORM(Radon_result[i, *]), 256)
    signal_Amplitude_1d = signalKernal_2d_array # image_1d
    rec_image_1d = Gradient_invert # signal_Amplitude_1d
    rec_remp = CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
    ;do filter
    rec_result[i, *] = Filter_FFT(rec_remp);rec_remp; Filter_FFT(rec_remp)
    Radon_result[i, *] = Filter_FFT(REFORM(Radon_result[i, *]))
    print, 'i = ', i
    plot, REFORM(Radon_result[i, *]), color=500
    oplot, CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
  ENDFOR

  backproject1 = RADON(Radon_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject1, title='FBP'

  backproject2 = RADON(rec_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject2, title='MPT+FBP'
END


;GJ, 2021/6/30, simulation is much faster
;GJ, 2021/7/4, simulation with single wire
;GJ, 2021/7/11, adding slice selection
PRO MPT_2d_single_wire_slice_selection

  ;  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20210627c\mp.bmp'
  ;  image_temp2 = read_bmp(filename)
  ;  image = MAX(image_temp2) - image_temp2

;  ;   filename2 = 'C:\D_drive\MPI_Tianjie\MPT_english.bmp'
;  filename2 = 'C:\D_drive\MPI_Tianjie\MagneticPT_english.bmp'
;  image_temp2 = read_bmp(filename2)
;  image = MAX(image_temp2) - image_temp2

  ;get the sample image
  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256))
  image = image * 0.
  image[127:129, 127:129] = 255
  image[100:103, 127:129] = 255
  
  iimage, image, title='original'
  ;  filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
  ;  image_temp2 = read_image(filename2)
  ;  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256), MAX=130)

  ;setup the 1d image reconstruction kernal
  Matrix_size = N_ELEMENTS(image[0,*]); 
  Matrix_size_double = 2*Matrix_size; + 20; 512; 
  field_amp = DBLARR(Matrix_size_double,Matrix_size_double)*0. ; 100 rows, 201 columns
  x_field = DBLARR(Matrix_size_double,Matrix_size_double)*0. ; 100 rows, 201 columns
  y_field = DBLARR(Matrix_size_double,Matrix_size_double)*0. ; 100 rows, 201 columns
  ;;calculate the magnetic field
  coeff=10000.
  FOR i=0., Matrix_size_double-1. DO BEGIN
    FOR j=0., Matrix_size_double-1. DO BEGIN
      dist = SQRT((Matrix_size_double/2.-i)^2 + (j+100)^2)
      ;
      ;H_flat_portion = 0.0
      ;;field_amp[i,j] = (signal_vs_H(coeff/dist),H_flat_portion,0)[0, 0] ;get the signal of single particle
      field_amp[i,j] = coeff/dist ;this is the H field
      x_field[i,j] = field_amp[i,j] * ((j+100) / dist) ;get the x direction signal of single particle
      y_field[i,j] = field_amp[i,j] * ((Matrix_size_double/2.-i) / dist) ;get the y direction signal of signal particle
    ENDFOR
  ENDFOR

  x_signal_i = DBLARR(Matrix_size)
  x_baseline_signal = DBLARR(Matrix_size)
  y_signal_i = DBLARR(Matrix_size)
  y_baseline_signal = DBLARR(Matrix_size)
  FOR i=0, Matrix_size-1 DO BEGIN
    x_field_matrix = x_field[(Matrix_size-i):(Matrix_size*2-i-1), 0:(Matrix_size-1)]
    x_signal_i[i] = TOTAL(x_field_matrix * image, /NAN)
    x_baseline_signal[i] = TOTAL(x_field_matrix, /NAN)*MEAN(image, /NAN)
    
    y_field_matrix = y_field[(Matrix_size-i):(Matrix_size*2-i-1), 0:(Matrix_size-1)]
    y_signal_i[i] = TOTAL(y_field_matrix * image, /NAN)
    y_baseline_signal[i] = TOTAL(y_field_matrix, /NAN)*MEAN(image, /NAN)
  ENDFOR
  
  window, 0
  plot, x_signal_i, color=500
  oplot, x_baseline_signal
  
  window, 1
  plot, y_signal_i, color=500
  oplot, y_baseline_signal
  
  y_single_signal = DBLARR(Matrix_size-2)
  x_single_signal = DBLARR(Matrix_size-2)
  y_baseline_single_signal = DBLARR(Matrix_size-2)
  x_baseline_single_signal = DBLARR(Matrix_size-2)
  FOR i=1, Matrix_size-2 DO BEGIN
    y_single_signal[i-1] = 0.5 * (y_signal_i[i-1] + y_signal_i[i+1]) - y_signal_i[i]
    x_single_signal[i-1] = x_signal_i[i] - 0.5 * (x_signal_i[i-1] + x_signal_i[i+1])
    
    y_baseline_single_signal[i-1] = 0.5 * (y_baseline_signal[i-1] + y_baseline_signal[i+1]) - y_baseline_signal[i]
    x_baseline_single_signal[i-1] = x_baseline_signal[i] - 0.5 * (x_baseline_signal[i-1] + x_baseline_signal[i+1])
  ENDFOR
  
  window, 2
  plot, x_single_signal, color=500
  oplot, x_baseline_single_signal
  
  window, 3
  plot, y_single_signal, color=500
  oplot, y_baseline_single_signal
  

  
  
  ;take the corner part of the field
  H_2d_array = x_field[0:Matrix_size-1, 0:Matrix_size-1]


  ;  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  ;  ratios = FINDGEN(Matrix_size)/(Matrix_size-1)
  ;  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = FINDGEN(Matrix_size)*ratios[i] + 200.

  ;define the signal kernal
  signalKernal_2d_array = H_2d_array*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[0, 0]
    ENDFOR
  ENDFOR
  ;calculate the invert matrix of the signal Kernal
  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  result = Gradient_invert # (signalKernal_2d_array)
  iimage, result, title='inverse matrix multiplication'

  iimage,  image, title='Original'
  ;;Calculate and display the Radon transform:
  Radon_result = RADON(image, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)

  rec_result = Radon_result*0.
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  FOR i=0, N_ELEMENTS(theta)-1 DO BEGIN
    image_1d = CONGRID(REFORM(Radon_result[i, *]), 256)
    signal_Amplitude_1d = signalKernal_2d_array # image_1d
    rec_image_1d = Gradient_invert # signal_Amplitude_1d
    rec_remp = CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
    ;do filter
    rec_result[i, *] = Filter_FFT(rec_remp);rec_remp; Filter_FFT(rec_remp)
    Radon_result[i, *] = Filter_FFT(REFORM(Radon_result[i, *]))
    print, 'i = ', i
    plot, REFORM(Radon_result[i, *]), color=500
    oplot, CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
  ENDFOR

  backproject1 = RADON(Radon_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject1, title='FBP'

  backproject2 = RADON(rec_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject2, title='MPT+FBP'
END


;GJ, 2021/7/10, simulation on 2d (x-z)
PRO MPT_xz_XidianUniv
  
  ;get the sample image
  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256))
;  image = image * 0. + 128.
  iimage,  image, title='Original'
  ;setup the 1d image reconstruction kernal
  x_size = N_ELEMENTS(image[*,0])
  z_size = N_ELEMENTS(image[0,*])
  
  ;set up the encoding magnetic fields
  H0 = 0.02 ;T
  Gx = 0.1; T/m
  Gz  =0.1; T/m
  FOV = 0.2; m
  x_res = FOV/x_size
  z_res = FOV/z_size
  Hx = (ABS(FINDGEN(x_size*2)-x_size)*Gx*x_res + H0) * 1000. ;mT
  Hz = (ABS(FINDGEN(z_size*2)-z_size)*Gz*z_res + H0) * 1000. ;mT
  x_axis = FINDGEN(x_size*2)-x_size*x_res
  
  ;window, 0
  iplot, x_axis, Hx, title='x-direction H'
  H_flat_portion = 0.0
  FOR l=0, x_size*2-1 DO Hx[l] = (signal_vs_H(Hx[l],H_flat_portion,0))[0, 0]
  FOR l=0, z_size*2-1 DO Hz[l] = (signal_vs_H(Hz[l],H_flat_portion,0))[0, 0]
  
  ;calculat the signal
  signal = DBLARR(x_size, z_size) * 0.
  H_matrix = DBLARR(x_size, z_size) * 0.
  FOR i=0, x_size-1 DO BEGIN
    FOR k=0, z_size-1 DO BEGIN
      ;encoding magnetic field
      sel_Hx = Hx[i:(i+255)]
      sel_Hz = Hz[k:(k+255)]
      
      ;put the H into a matrix
      H_matrix = H_matrix * 0.
      FOR m=0, x_size-1 DO H_matrix[m, *] = H_matrix[m, *] + sel_Hz
      FOR n=0, z_size-1 DO H_matrix[*, n] = H_matrix[*, n] + sel_Hx
      
;      ;calculat the signal
;      
      signal[i, k] = TOTAL(H_matrix*image, /NAN)
      
    ENDFOR
    print, 'i, k = ', i, k
  ENDFOR
  
  iimage, signal, title='Signal'
  
END

;GJ, 2021/6/30, simulation is much faster
;GJ, 2021/7/11, 1d V-shape field encoding and deconding + 1d FBP
;GJ, 2021/7/28, file is moved to default folder
PRO aimis3d_MPT_2d_V_field

  filename = '.\MPI_images\mp.bmp'
  image_temp2 = read_bmp(filename)
  image = MAX(image_temp2) - image_temp2

  ;  ;get the sample image
  ;  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
  ;  image_temp2 = read_image(filename2)
  ;  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256))

  ;  filename2 = 'C:\D_drive\MPI_Tianjie\MPT_english.bmp'
  ;;  filename2 = 'C:\D_drive\MPI_Tianjie\MagneticPT_english.bmp'
  ;  image_temp2 = read_bmp(filename2)
  ;  image = MAX(image_temp2) - image_temp2

  ;  filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
  ;  image_temp2 = read_image(filename2)
  ;  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256), MAX=130)
  ;  image = image * 0.
  ;  image[100:102, 127:129] = 256
  ;  image[130:132, 127:129] = 256

  ;setup the 1d image encoding and decoding kernal
  iimage,  image, title='Original'
  Matrix_size = N_ELEMENTS(image[0,*]); causing invert calcuation error
  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  x_size = N_ELEMENTS(image[*,0])
  ;set up the encoding magnetic fields
  H0 = 0.02 ;T
  Gx = 0.1; T/m
  FOV = 0.2; m
  x_res = FOV/x_size
  Hx = (ABS(FINDGEN(x_size*2)-x_size)*Gx*x_res + H0) * 1000. ;mT
  x_axis = (FINDGEN(x_size*2)-x_size)*x_res * 1000. ;mm
  iplot, x_axis, Hx
  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = Hx[i:i+Matrix_size-1]

  ;define the signal kernal
  signalKernal_2d_array = H_2d_array*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[0, 0]
    ENDFOR
  ENDFOR
  ;calculate the invert matrix of the signal Kernal
  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  result = Gradient_invert # signalKernal_2d_array
  iimage, result, title='invert matrix *'

  ;;Calculate and display the Radon transform:
  Radon_result_ori = RADON(image, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)
  Radon_result = NOISE_PICK(Radon_result_ori, 0.1, ITER=6)
  ; Display the images side by side:
  IIMAGE, Radon_result_ori, VIEW_GRID=[2,1], VIEW_TITLE='Original', $
    DIMENSIONS=[850, 550], WINDOW_TITLE='NOISE_PICK Example', $
    /NO_SAVEPROMPT
  IIMAGE, Radon_result, /VIEW_NEXT, VIEW_TITLE='Noisy'

  window, 2
  TVSCL, radon_result,  .08, .08, /NORMAL
  PLOT, theta, rho, POSITION=[0.08, 0.08, 1, 0.98], /NODATA, /NOERASE, XSTYLE=9,YSTYLE=9,XTITLE='Theta', YTITLE='R'
  rec_result = Radon_result*0.
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  FOR i=0, N_ELEMENTS(theta)-1 DO BEGIN
    image_1d = CONGRID(REFORM(Radon_result[i, *]), 256)
    signal_Amplitude_1d = signalKernal_2d_array # image_1d
    rec_image_1d = Gradient_invert # signal_Amplitude_1d
    rec_remp = CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
    ;do filter
    rec_result[i, *] = Filter_FFT(rec_remp)
    Radon_result[i, *] = Filter_FFT(REFORM(Radon_result[i, *]))
    print, 'i = ', i
    ;window, 3
    ;plot, REFORM(Radon_result[i, *]), color=500
    ;oplot, CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))

    window, 4
    image_1d_ori = CONGRID(REFORM(Radon_result_ori[i, *]), 256)
    plot, image_1d_ori, color=500
    oplot, rec_image_1d, linestyle=1
  ENDFOR

  backproject1 = RADON(Radon_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject1, title='FBP'

  backproject2 = RADON(rec_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject2, title='1d V-shape MPT+FBP'
END



;GJ, 2021/6/30, simulation is much faster
;GJ, 2021/7/11, 1d V-shape field encoding and deconding + 1d FBP
PRO MPT_2d_V_field
  
  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20210627c\mp.bmp'
  image_temp2 = read_bmp(filename)
  image = MAX(image_temp2) - image_temp2
  
;  ;get the sample image
;  filename2='C:\D_drive\MPI_Tianjie\XidianUniversity_ori.png'
;  image_temp2 = read_image(filename2)
;  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256))
  
;  filename2 = 'C:\D_drive\MPI_Tianjie\MPT_english.bmp'
;;  filename2 = 'C:\D_drive\MPI_Tianjie\MagneticPT_english.bmp'
;  image_temp2 = read_bmp(filename2)
;  image = MAX(image_temp2) - image_temp2
  
;  filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
;  image_temp2 = read_image(filename2)
;  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256), MAX=130)
;  image = image * 0.
;  image[100:102, 127:129] = 256
;  image[130:132, 127:129] = 256
  
  ;setup the 1d image encoding and decoding kernal
  iimage,  image, title='Original'
  Matrix_size = N_ELEMENTS(image[0,*]); causing invert calcuation error
  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  x_size = N_ELEMENTS(image[*,0])
  ;set up the encoding magnetic fields
  H0 = 0.02 ;T
  Gx = 0.1; T/m
  FOV = 0.2; m
  x_res = FOV/x_size
  Hx = (ABS(FINDGEN(x_size*2)-x_size)*Gx*x_res + H0) * 1000. ;mT
  x_axis = (FINDGEN(x_size*2)-x_size)*x_res * 1000. ;mm
  iplot, x_axis, Hx
  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = Hx[i:i+Matrix_size-1]

  ;define the signal kernal
  signalKernal_2d_array = H_2d_array*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[0, 0]
    ENDFOR
  ENDFOR
  ;calculate the invert matrix of the signal Kernal
  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  result = Gradient_invert # signalKernal_2d_array
  iimage, result, title='invert matrix *'
  
  ;;Calculate and display the Radon transform:
  Radon_result_ori = RADON(image, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)
  Radon_result = NOISE_PICK(Radon_result_ori, 0.1, ITER=6)
  ; Display the images side by side:
  IIMAGE, Radon_result_ori, VIEW_GRID=[2,1], VIEW_TITLE='Original', $
    DIMENSIONS=[850, 550], WINDOW_TITLE='NOISE_PICK Example', $
    /NO_SAVEPROMPT
  IIMAGE, Radon_result, /VIEW_NEXT, VIEW_TITLE='Noisy'
  
  window, 2
  TVSCL, radon_result,  .08, .08, /NORMAL 
  PLOT, theta, rho, POSITION=[0.08, 0.08, 1, 0.98], /NODATA, /NOERASE, XSTYLE=9,YSTYLE=9,XTITLE='Theta', YTITLE='R'
  rec_result = Radon_result*0.
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  FOR i=0, N_ELEMENTS(theta)-1 DO BEGIN
    image_1d = CONGRID(REFORM(Radon_result[i, *]), 256)
    signal_Amplitude_1d = signalKernal_2d_array # image_1d
    rec_image_1d = Gradient_invert # signal_Amplitude_1d
    rec_remp = CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
    ;do filter
    rec_result[i, *] = Filter_FFT(rec_remp)
    Radon_result[i, *] = Filter_FFT(REFORM(Radon_result[i, *]))
    print, 'i = ', i
    ;window, 3
    ;plot, REFORM(Radon_result[i, *]), color=500
    ;oplot, CONGRID(REFORM(rec_image_1d), N_ELEMENTS(rho))
    
    window, 4
    image_1d_ori = CONGRID(REFORM(Radon_result_ori[i, *]), 256)
    plot, image_1d_ori, color=500
    oplot, rec_image_1d, linestyle=1
  ENDFOR
  
  backproject1 = RADON(Radon_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject1, title='FBP'

  backproject2 = RADON(rec_result, /BACKPROJECT, RHO=rho, THETA=theta)
  iimage, backproject2, title='1d V-shape MPT+FBP'
END

PRO MPT_oneD_Simulation_SingleWire
  Matrix_size = 512; causing invert calcuation error
  field_amp = DBLARR(Matrix_size,Matrix_size)*0. ; 100 rows, 201 columns
  x_field = DBLARR(Matrix_size,Matrix_size)*0. ; 100 rows, 201 columns
  y_field = DBLARR(Matrix_size,Matrix_size)*0. ; 100 rows, 201 columns
  ;;calculate the magnetic field
  coeff=1000.
  FOR i=0., Matrix_size-1. DO BEGIN
    FOR j=0., Matrix_size-1. DO BEGIN
      dist = SQRT((Matrix_size/2-i)^2 + (j+1)^2)
      field_amp[i,j] = coeff/dist
      x_field[i,j] = field_amp[i,j] * ((j+1) / dist)
      y_field[i,j] = field_amp[i,j] * ((Matrix_size/2-i) / dist)
    ENDFOR
  ENDFOR

  Matrix_size = 256
  H_2d_array = x_field[0:Matrix_size-1, 0:Matrix_size-1]

  iimage, H_2d_array

  signalKernal_2d_array = H_2d_array*0.

  image_1d = DBLARR(Matrix_size)*0.
  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20210627c\mp.bmp'
  image_temp = read_bmp(filename)
  image = MAX(image_temp) - image_temp
  FOR i=0, Matrix_size-1 DO image_1d[i] = TOTAL(image[i,*])
  iplot, image_1d, YRANGE=[0, MAX(image_1d)*1.1], title='original image'
 
 
  signalKernal_2d_array = H_2d_array*0.

  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  signal_AUC_1d = DBLARR(Matrix_size)*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signal = signal_vs_H(H_2d_array[i, j],H_flat_portion,0)*image_1d[j]
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[0, 0]
      signal_Amplitude_1d[i] = signal_Amplitude_1d[i] + signal[0, 0]
      signal_AUC_1d[i] = signal_AUC_1d[i] + signal[0, 1]
    ENDFOR
  ENDFOR

  ;;reconstruction
  test_N = Matrix_size-1
  Gradient_invert = LA_INVERT(signalKernal_2d_array[0:test_N, 0:test_N], /DOUBLE, status =status)
  result = Gradient_invert # (signalKernal_2d_array[0:test_N, 0:test_N])
  iimage, result
  
  ;  meas = signalKernal_2d_array # image_1d
  test = Gradient_invert # signal_Amplitude_1d
  plot, test, color=500
  oplot, image_1d


END

PRO MPT_oneD_Simulation
  Matrix_size = 256; causing invert calcuation error
  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  H_gradient = FINDGEN(Matrix_size) ;0 to 255
  ratios = FINDGEN(Matrix_size)/(Matrix_size-1)
  ;FOR i=0, Matrix_size/2.-1 DO BEGIN
  ;  H_2d_array[i, *] = REVERSE(FINDGEN(Matrix_size))*(REVERSE(ratios))[i] + 200.
  ;ENDFOR
  
  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = FINDGEN(Matrix_size)*ratios[i] + 200.
  
  iimage, H_2d_array
  
  signalKernal_2d_array = H_2d_array*0.
  
  image_1d = DBLARR(Matrix_size)*0.
  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20210627c\mp.bmp'
  image_temp = read_bmp(filename)
  image = MAX(image_temp) - image_temp
  FOR i=0, Matrix_size-1 DO image_1d[i] = TOTAL(image[i,*])
  
  ;image_1d[50:55] = 1.
  ;image_1d[56:60] = 2.
  ;image_1d[61:65] = 3.
  ;image_1d[66:70] = 2.
  ;image_1d[71:75] = 1.
  ;image_1d[120:135] = 1.
  ;image_1d[*] = 1.
  iplot, image_1d, YRANGE=[0, MAX(image_1d)*1.1], title='original image'
  
  signal_Amplitude_1d = DBLARR(Matrix_size)*0.
  signal_AUC_1d = DBLARR(Matrix_size)*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signal = signal_vs_H(H_2d_array[i, j],H_flat_portion,0)*image_1d[j]
      signalKernal_2d_array[i, j] = (signal_vs_H(H_2d_array[i, j],H_flat_portion,0))[0, 0]
      signal_Amplitude_1d[i] = signal_Amplitude_1d[i] + signal[0, 0]
      signal_AUC_1d[i] = signal_AUC_1d[i] + signal[0, 1]
    ENDFOR
  ;  print, 'signal_Amplitude_1d = ', signal_Amplitude_1d
      IF i EQ 0 THEN BEGIN
        plot, signal_Amplitude_1d, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[-1.1, MAX(signal_Amplitude_1d)*1.1], title='Signal vs Gradient'
      ENDIF ELSE BEGIN
        oplot, signal_Amplitude_1d, COLOR = 0;, YRANGE=[-1.1, 1.1]
      ENDELSE
  ENDFOR
  
  iplot, signal_Amplitude_1d, title='Amplitude vs Gradient', YRANGE=[0, MAX(signal_Amplitude_1d)*1.1]
  iplot, signal_AUC_1d, title='AUC vs Gradient', YRANGE=[0, MAX(signal_AUC_1d)*1.1]
  
  ;;reconstruction
  ;signalKernal_2d_array[Matrix_size/2-1, *] = signalKernal_2d_array[0, *]
  ;signalKernal_2d_array[Matrix_size/2, *] = signalKernal_2d_array[0, *]
  Gradient_invert = LA_INVERT(signalKernal_2d_array, /DOUBLE, status =status)
  meas = signalKernal_2d_array # image_1d
  test = Gradient_invert # signal_Amplitude_1d
  result = Gradient_invert # (signalKernal_2d_array)
  iplot, test, title='reconstructed image', YRANGE=[0, MAX(test)*1.1]
  ;iimage, result
  
  ;ssignalKernal_2d_array = signalKernal_2d_array[0:110, 0:110]
  ;sGradient_invert = LA_INVERT(ssignalKernal_2d_array, /DOUBLE, status =status)
  ;result = sGradient_invert # (ssignalKernal_2d_array)
  ;stest = sGradient_invert # signal_Amplitude_1d[0:100]
  ;iimage, result
  
  ;G_invert = MoorePmatrix(H_2d_array-MEAN(H_2d_array))
  ;;G_invert = INVERT(H_2d_array, /DOUBLE)
  ;result = G_invert # (H_2d_array-MEAN(H_2d_array))
  ;test2 = G_invert # signal_Amplitude_1d
  ;plot, test2
END

PRO MPT_twoD_Simulation
  Matrix_size = 256; causing invert calcuation error
  H_2d_array_x = DBLARR(Matrix_size, Matrix_size)
  H_2d_array_y = DBLARR(Matrix_size, Matrix_size)
  H_gradient = FINDGEN(Matrix_size) ;0 to 255
  ratios = FINDGEN(Matrix_size)/(Matrix_size-1)
  
  ;H_2d_array_x and H_2d_array_y they are the same
  FOR i=0, Matrix_size-1 DO H_2d_array_x[i, *] = FINDGEN(Matrix_size)*ratios[i] + 200.
  FOR j=0, Matrix_size-1 DO H_2d_array_y[*, j] = FINDGEN(Matrix_size)*ratios[j] + 200.

;  iimage, H_2d_array_y

  signalKernal_2d_array_x = H_2d_array_x*0.
  signalKernal_2d_array_y = H_2d_array_y*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_flat_portion = 0.0
      signalKernal_2d_array_x[i, j] = (signal_vs_H(H_2d_array_x[i, j],H_flat_portion,0))[0, 0]
      signalKernal_2d_array_y[i, j] = (signal_vs_H(H_2d_array_y[i, j],H_flat_portion,0))[0, 0]
    ENDFOR
  ENDFOR

;  image_1d = DBLARR(Matrix_size)*0.
  filename = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20210627c\mp.bmp'
  image_temp = read_bmp(filename)
  MPI_image = MAX(image_temp) - image_temp
  ;MPI_image = MPI_image*0.
  ;MPI_image[Matrix_size/2-5:Matrix_size/2+5, Matrix_size/2-5:Matrix_size/2+5] = 10

  
  signal_Amplitude_temp = signalKernal_2d_array_y # MPI_image
  signal_Amplitude = signalKernal_2d_array_x # signal_Amplitude_temp
  ;signal_Amplitude = signalKernal_2d_array_x # MPI_image
  
  iimage, signal_Amplitude, title='Amplitude vs Gradient'
  
  ;;reconstruction
  Gradient_invert_x = LA_INVERT(signalKernal_2d_array_x, /DOUBLE, status =status)
  Gradient_invert_y = LA_INVERT(signalKernal_2d_array_y, /DOUBLE, status =status)
  
  Gradient_invert_x = Gradient_invert_x/MAX(gradient_invert_x)
  Gradient_invert_y = Gradient_invert_y/MAX(gradient_invert_y)
  ;meas = signalKernal_2d_array # image_1d
  ;test = Gradient_invert_y # Gradient_invert_x # signal_Amplitude
  recon_temp = Gradient_invert_x # signal_Amplitude
  recon = Gradient_invert_y # recon_temp
  recon_test = Gradient_invert_y # signal_Amplitude_temp
  result_x = Gradient_invert_x # (signalKernal_2d_array_x)
  result_y = Gradient_invert_y # (signalKernal_2d_array_y)
  iimage, recon, title='reconstructed image'
  iimage, recon_test, title='reconstructed image ori'
  ;iimage, result_x
  ;iimage, result_y
END

;process the total signal and separate a single layer
;GJ, Xidian University, 2021/6/26
PRO Signal_Single_layer

  N=40   ;number of magnetic field
  Ntp=1200   ;number of time points
  H_array = (FINDGEN(N)+1.)*10.
  FWHM_array = DBLARR(N)*0.
  AUC_array = DBLARR(N)*0.
  Mrange_array = DBLARR(N)*0.
  T_Signal_array = DBLARR(N, Ntp-1)*0.

  X = 2*!PI/1000. * FINDGEN(Ntp)
  FOR i=0, N-1 DO BEGIN
    H = -COS(X)*H_array[i]*0.1
    M = 1.0*(1./TANH(H) - 1./H)
    Mrange_array[i]=MAX(M, /NAN)-MIN(M, /NAN)
    Nelem = N_ELEMENTS(M)
    signal = M[1:*]-M[0:Nelem-1]
    Nelem = N_ELEMENTS(M)

    ;calcaluate the total signal
    IF i EQ 0 THEN BEGIN
      T_Signal_array[i,*] = signal
    ENDIF ELSE BEGIN
      T_Signal_array[i,*] = T_Signal_array[i-1,*]+signal
    ENDELSE

        
;        IF i EQ 0 THEN BEGIN
;          window, 0
;          plot, X, M, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[-1.1, 1.1]
;        ENDIF ELSE BEGIN
;          wset, 0
;          oplot, X, M, COLOR = 0;, YRANGE=[-1.1, 1.1]
;        ENDELSE
;        wait, 0.4
    ;    oplot, X[1:*], 10.*(M[1:*]-M[0:Nelem-1]), COLOR = 0;, YRANGE=[-0.05, 0.05]
    ;    wait, 0.4


    time_range = where(X GT 5.3-2. AND X LT 7.3+2., NCounts)
    signal_range = signal[time_range+1]
    AUC_array[i] = TOTAL(T_signal_array[i,*], /NAN)
    T_signal_range = T_Signal_array[i, time_range+1]
   
    IF i EQ 0 THEN BEGIN
       window, 1
       plot, time_range, T_signal_range, title='a signal peak', BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[0, 0.03]
    ENDIF ELSE BEGIN
      wset, 1
      oplot, time_range, T_signal_range, COLOR = i*100
    ENDELSE
    wait, 0.2

    ;   result = POLY_FIT(time_range, signal_range, 15, MEASURE_ERRORS=measure_errors, SIGMA=sigma, YFIT=yfit)
    max_signal = MAX(T_signal_range, maxind)
    lower_signal_sub = ABS(T_signal_range[0:maxind]-max_signal/2)
    higher_signal_sub = ABS(T_signal_range[maxind+1:*]-max_signal/2)
    ;  plot, time_range, [lower_signal_sub, higher_signal_sub], title='signal minus half_max', BACKGROUND = 'FFFFFF'x, COLOR = 0
    min_lower_signal_sub = MIN(lower_signal_sub, lower_index)
    min_higher_signal_sub = MIN(higher_signal_sub, higher_index)
    FWHM = higher_index + maxind - lower_index
    FWHM_array[i] = FWHM*(X[1]-X[0])
    ;    print, signal
  ENDFOR

  window, 2
  FOR i=0, N-2 DO BEGIN
    IF i EQ 0 THEN BEGIN
      plot, X[1:*], T_Signal_array[N-1-i,*]-T_Signal_array[N-2-i,*], title='a signal total', BACKGROUND = 'FFFFFF'x, COLOR = 0, XRANGE=[0, 1]
    ENDIF ELSE BEGIN
      oplot, X[1:*], T_Signal_array[N-1-i,*]-T_Signal_array[N-2-i,*], COLOR = i*100
    ENDELSE
    wait, 0.8
  ENDFOR
 
 window, 3
 M=14
 plot, X[1:*], T_Signal_array[M,*], title='a signal total', BACKGROUND = 'FFFFFF'x, COLOR = 0;, XRANGE=[1,2]
 oplot, X[1:*], T_Signal_array[M,*]-T_Signal_array[4,*], color=1222
 oplot, X[1:*], T_Signal_array[4,*], color=1222
 
  ;oplot, X[1:*], T_Signal_array[N-1,*]-T_Signal_array[20,*], color=1222
  print, 'test'
 ;  window, 1
;  plot, H_array, FWHM_array, title='H vs FWHM'
;  print, transpose(H_array)
;  print, transpose(FWHM_array)
;  ;check their relationship
;
;  window, 4
  iplot,  H_array, AUC_array, title='H vs Area of Peak'
;
;
;  window, 3
;  plot,  H_array, AUC_array/Mrange_array, title='H vs Area of Peak/M Range'
;  print, transpose(AUC_array)
END

;GJ, 2021/11/28, field free slice
PRO FFS_parallel_coil_plans

  field_amp = DBLARR(401,201)*0. ; 100 rows, 201 columns
  x_field = DBLARR(401,201)*0. ; 100 rows, 201 columns
  y_field = DBLARR(401,201)*0. ; 100 rows, 201 columns

  coeff=100.
  FOR i=0, 400 DO BEGIN
    FOR j=0, 200 DO BEGIN
      dist = SQRT((i - 200.)^2 + (j - 100.)^2)
      field_amp[i,j] = 1./dist;*coeff
      x_field[i,j] = field_amp[i,j] * (j - 100.)
      y_field[i,j] = -field_amp[i,j] * (i - 200.)
    ENDFOR
  ENDFOR
;  iimage, field_amp, title = 'Field Amplitude of single line'
;  iimage, x_field, title = 'x-direction Field of single line'
;  iimage, y_field, title = 'y-direction Field of single line'
  
  sum_x_field = DBLARR(99,99)*0.
  sum_y_field = DBLARR(99,99)*0.
  sum_field_amp = DBLARR(99,99)*0.
  FOR k=100, 300, 50 DO BEGIN
    sum_x_field = sum_x_field + x_field[k-49:k+49, 1:99] + x_field[k-49:k+49, 101:199]
    sum_y_field = sum_y_field + y_field[k-49:k+49, 1:99] + y_field[k-49:k+49, 101:199]
  ENDFOR
  
  sum_field_amp = SQRT(sum_x_field^2 + sum_y_field^2)
  
  iimage, sum_field_amp, title='Field Amplitude of Line Coils'
  isurface, sum_field_amp, title='Field Amplitude of Line Coils'
  contour, sum_field_amp, levels=FINDGEN(100)/99.*max(sum_field_amp), title='Field Amplitude'
  iimage, sum_x_field, title = 'x-direction Field of Line Coils'
  iimage, sum_y_field, title = 'y-direction Field of Line Coils'
  
  iplot, sum_field_amp[*, 50], title='Field Amplitude on x-axis'
  iplot, sum_field_amp[50, *], title='Field Amplitude on y-axis'
  
  
END


;GJ, 2021/11/25, 4 parallel lines with the same field direction
PRO Anderson_coils_magnetic_field

  field_amp = DBLARR(201,201)*0. ; 100 rows, 201 columns
  x_field = DBLARR(201,201)*0. ; 100 rows, 201 columns
  y_field = DBLARR(201,201)*0. ; 100 rows, 201 columns

  coeff=100.
  FOR i=0, 200 DO BEGIN
    FOR j=0, 200 DO BEGIN
      dist = SQRT((i - 100.)^2 + (j - 100.)^2)
      field_amp[i,j] = 1./dist;*coeff
      x_field[i,j] = field_amp[i,j] * (j - 100.)
      y_field[i,j] = -field_amp[i,j] * (i - 100.)
    ENDFOR
  ENDFOR
  iimage, field_amp, title = 'Field Amplitude of single line'
  iimage, x_field, title = 'x-direction Field of single line'
  iimage, y_field, title = 'y-direction Field of single line'
  
  sum_x_field = DBLARR(100,100)*0.
  sum_y_field = DBLARR(100,100)*0.
  sum_field_amp = DBLARR(100,100)*0.
  sum_x_field = x_field[101:200, 0:99] - x_field[0:99, 0:99] + x_field[0:99, 101:200] - x_field[101:200, 101:200] 
  sum_y_field = y_field[101:200, 0:99] - y_field[0:99, 0:99] + y_field[0:99, 101:200] - y_field[101:200, 101:200] 
  
  sum_field_amp = SQRT(sum_x_field^2 + sum_y_field^2)
  
  iimage, sum_field_amp, title='Field Amplitude of Anderson coils'
  isurface, sum_field_amp, title='Field Amplitude of Anderson coils'
  contour, sum_field_amp, levels=FINDGEN(100)/99.*max(sum_field_amp), title='Field Amplitude'
  iimage, sum_x_field, title = 'x-direction Field of Anderson coils'
  iimage, sum_y_field, title = 'y-direction Field of Anderson coils'
  
  iplot, sum_field_amp[*, 50], title='Field Amplitude on x-axis'
  iplot, sum_field_amp[50, *], title='Field Amplitude on y-axis'
END

PRO coils_magnetic_field

field_amp = DBLARR(201,100)*0. ; 100 rows, 201 columns
x_field = DBLARR(201,100)*0. ; 100 rows, 201 columns
y_field = DBLARR(201,100)*0. ; 100 rows, 201 columns

;;check the coordiantes
;x_field[0,0]=255
;x_field[100,0] = 255
;iimage, x_field

;;calculate the magnetic field
coeff=1000.
FOR i=0, 200 DO BEGIN
  FOR j=0, 99 DO BEGIN
    dist = SQRT((100.-i)^2 + j^2)
    field_amp[i,j] = 1./dist*coeff
    x_field[i,j] = field_amp[i,j] * (j / dist)
    y_field[i,j] = field_amp[i,j] * ((100.-i) / dist)    
  ENDFOR
ENDFOR

;iimage, BYTSCL(field_amp, /NAN)
;iimage, BYTSCL(x_field)
;iimage, BYTSCL(y_field)
;iimage, field_amp

N=100
sum_field_amp = DBLARR(N,N)*0.
sum_x = DBLARR(N,N)*0.
sum_y = DBLARR(N,N)*0.

coil_ind = FINDGEN(N)
ratio_ind = FINDGEN(N)*0.

X = 2.*!PI/N * FINDGEN(N)
ratio_ind = SIN(X)+2.

FOR k=0, N-1 DO BEGIN
  ;;ratio = ABS(50. - SQRT(k*(N-1-k)))/50.+2.
  
  ;ratio = ABS(50. - k)/100. + 1.
  ;ratio = k/(N-1.)
  ratio = 1.
  ;ratio = ratio_ind[k] 
  sum_field_amp = sum_field_amp + ratio*field_amp[k:k+99, *]
  sum_x = sum_x + ratio*x_field[k:k+99, *]
  sum_y = sum_y + ratio*y_field[k:k+99, *]
  
ENDFOR


plot, coil_ind, ratio_ind, BACKGROUND = 'FFFFFF'x, COLOR = 0, YRANGE=[0, MAX(ratio_ind)+0.2], title='coil electric field profile'
iimage, BYTSCL(sum_field_amp, /NAN), title='field amplitude'
iimage, BYTSCL(sum_x, /NAN), title='x direction'
iimage, BYTSCL(sum_y, /NAN), title='y direction'
central_y = sum_y[33:67, 1:34]
iimage, BYTSCL(central_y, /NAN), title='y direction'

END


Function MoorePmatrix, matrix

  nc = N_Elements(matrix[*, 0])    ;求矩阵列数

  nr = N_Elements(matrix[0, *])    ;求矩阵行数

  rank = 0               ;矩阵的秩初始化

  row = IntArr(nc)   ;满秩列的标记

  tempspectral = FLtArr(nc, nr)  ;定义变换矩阵

  For i = 0, nc-1 Do tempspectral[i, *] = matrix [i, *]

  ;对矩阵高斯约化，获得行变换最简式

  For i=0, nc-1 Do Begin

    ;寻找列主元即其位置

    For j=rank, nr-1 Do Begin

      me=0.0

      me=Max(Abs(tempspectral[i,j:nr-1]))

      If me LE 0.0001 Then Break

      For k=j, nr-1 Do Begin

        If Abs(tempspectral[i,k]) EQ me Then Break

      EndFor

      ;交换行主元

      temp=tempspectral[*,j]

      tempspectral[*,j]=tempspectral[*,k]

      tempspectral[*,k]=temp

      ;消元

      m=tempspectral[i, j]

      tempspectral[*, j] = tempspectral[*,j]/m

      For k=0, nr-1 Do Begin

        If abs(tempspectral[i,k]) LE 0.0001 Then Break

        If k EQ j Then Continue

        n=tempspectral[i,k]

        tempspectral[i:nc-1,k]=tempspectral[i:nc-1,k]-n*tempspectral[i:nc-1,j]

      EndFor

      row[rank] = i

      rank = rank + 1

    EndFor

  EndFor

  ;确定分解矩阵

  rmatrix = FltArr(rank, nr) ;定义列满秩矩阵

  For i = 0, rank-1 Do rmatrix[i, *] = matrix[row[i], *]

  cmatrix = FltArr(nc, rank) ;定义行满秩矩阵

  cmatrix = tempspectral[*, 0:rank-1]

  ;利用满秩分解求广义逆

  mpmatrix = Transpose(cmatrix) ## Invert(cmatrix ## Transpose(cmatrix))## Invert(Transpose(rmatrix) ## rmatrix) ## Transpose(rmatrix)

  return, mpmatrix

End

;GJ, 2021/8/30, get 3D magnetic field based on gradient
;input: G in mT/m
;       theta in degree
;       phi in degree
FUNCTION H_volume_from_G, G, theta, phi, H_volume_size, FOV
  ;G = 5.0
  ;theta = 0.
  ;phi = 90.
  ;H_volume_size = [256, 256, 256]
  ;FOV = 0.20 ; m
  pixel_size = FOV/H_volume_size
  Gradient = Gradient_spherical_to_rectangular(G, theta, phi) ; mT/m
  H_volume = DBLARR(H_volume_size[0], H_volume_size[1], H_volume_size[2])*0.

  ;calculate from x direction
  FOR i = 0, H_volume_size[0]-1 DO BEGIN
    H_value = (i - H_volume_size[0]/2.) * pixel_size[0] * Gradient[0]
    H_volume[i, *, *] = H_volume[i, *, *] + H_value
  ENDFOR

  ;calculate from y direction
  FOR j = 0, H_volume_size[1]-1 DO BEGIN
    H_value = (j - H_volume_size[1]/2.) * pixel_size[1] * Gradient[1]
    H_volume[*, j, *] = H_volume[*, j, *] + H_value
  ENDFOR

  ;calculate from z direction
  FOR k = 0, H_volume_size[2]-1 DO BEGIN
    H_value = (k - H_volume_size[2]/2.) * pixel_size[2] * Gradient[2]
    H_volume[*, *, k] = H_volume[*, *, k] + H_value
  ENDFOR

;  iimage, REFORM(H_volume[H_volume_size[0]/2, *, *])
;  iimage, REFORM(H_volume[*, H_volume_size[1]/2, *])
;  iimage, REFORM(H_volume[*, *, H_volume_size[2]/2])

  RETURN, H_volume ; in mT
END

;GJ, 2021/8/30, get magnetic field based on gradient field encoding
;input: G in mT/m
;       theta in degree
;       phi in degree
;If H_volume_size (input) is 1d, eg. [256], H_volume (output) is [256, 256] with 256 encoding steps, each setp has 1d magentic field
;If H_volume_size (input) is 2d, eg. [256, 256], H_volume (output) is [256, 256, 256] with 256*256 encoding steps, each setp has 1d magentic field
;If H_volume_size (input) is 3d, eg. [256, 256, 256], H_volume (output) is [256, 256, 256, 256] with 256*256*256 encoding steps, each setp has 1d magentic field
FUNCTION H_from_G_encoding1, Gmax, H_volume_size, FOV
  ;Gmax = 5.0
  ;theta = 0.
  ;phi = 90.
  ;H_volume_size = [256]; [256, 256, 256]
  ;FOV = 0.20 ; m
  pixel_size = FOV/H_volume_size
;  Gradient = Gradient_spherical_to_rectangular(Gmax, theta, phi) ; mT/m
;  H_volume = DBLARR(H_volume_size[0], H_volume_size[1], H_volume_size[2])*0.
  
  IF N_ELEMENTS(H_volume_size) EQ 1 THEN BEGIN
    
    ;calculate the gradient
    theta = 0.
    phi = 90.
    Gradient_max = Gradient_spherical_to_rectangular(Gmax, theta, phi) ; mT/m
    
    H_volume = DBLARR(H_volume_size[0], H_volume_size[0])*0.
    delta_Gradient = Gradient_max[0]*2./(H_volume_size[0]-1)
    
    ;calculate H volume from x direction
    FOR i = 0, H_volume_size[0]-1 DO BEGIN
      Gradient = (i - H_volume_size[0]/2.) * delta_Gradient
      FOR j = 0, H_volume_size[0]-1 DO BEGIN
        H_value = j * pixel_size[0] * Gradient
        H_volume[i, j] = H_volume[i, j] + H_value
      ENDFOR
    ENDFOR
    
    iimage, H_volume
  ENDIF
  
;  ;calculate from x direction
;  FOR i = 0, H_volume_size[0]-1 DO BEGIN
;    H_value = (i - H_volume_size[0]/2.) * pixel_size[0] * Gradient[0]
;    H_volume[i, *, *] = H_volume[i, *, *] + H_value
;  ENDFOR
;  
;  ;calculate from y direction
;  FOR j = 0, H_volume_size[1]-1 DO BEGIN
;    H_value = (j - H_volume_size[1]/2.) * pixel_size[1] * Gradient[1]
;    H_volume[*, j, *] = H_volume[*, j, *] + H_value
;  ENDFOR
;  
;  ;calculate from z direction
;  FOR k = 0, H_volume_size[2]-1 DO BEGIN
;    H_value = (k - H_volume_size[2]/2.) * pixel_size[2] * Gradient[2]
;    H_volume[*, *, k] = H_volume[*, *, k] + H_value
;  ENDFOR
;  
;  iimage, REFORM(H_volume[H_volume_size[0]/2, *, *])
;  iimage, REFORM(H_volume[*, H_volume_size[1]/2, *])
;  iimage, REFORM(H_volume[*, *, H_volume_size[2]/2])
  
  RETURN, H_volume ; in mT
END

;GJ, 2021/8/29, convert spherical to rectangluar coordinates
;input: G is the gradient amplitude; theta in degree, phi in degree
FUNCTION Gradient_spherical_to_rectangular, G, theta, phi
  Gx = G * COS(!PI/180. * theta) * SIN(!PI/180. * phi)
  IF ABS(Gx) LT 0.0001 THEN Gx = 0.
  Gy = G * SIN(!PI/180. * theta) * SIN(!PI/180. * phi)
  IF ABS(Gy) LT 0.0001 THEN Gy = 0.
  Gz = G * COS(!PI/180. * phi)
  IF ABS(Gz) LT 0.0001 THEN Gz = 0.
  RETURN, [Gx, Gy, Gz]
END

;GJ, 2021/8/29, convert rectangular to spherical coordiantes
;input: Gx, Gy, Gz are MRI gradients along x, y, z directions
FUNCTION Gradient_rectangular_to_spherical, Gx, Gy, Gz
  G = SQRT(Gx^2 + Gy^2 + Gz^2)
  theta = 180/!PI * ATAN(Gy, Gx)
  phi = 180/!PI * ATAN(SQRT(Gx^2 + Gy^2), Gz)
  RETURN, [G, theta, phi];G is the gradient amplitude; theta in degree, phi in degree
END


;GJ, 2021/8/26, combination of gradient fields
PRO Gradient_excitation_3d
  N_ele = 256
  x_image = DBLARR(N_ele, N_ele, N_ele)
  
  FOR i=0, N_ele-1 DO x_image[i, *, *] = 0.5*i;ABS(N_ele/2-i)
  
  y_image = DBLARR(N_ele, N_ele, N_ele)
  ;FOR j=0, N_ele-1 DO y_image[*, j] = ABS(N_ele/2-j)
  FOR j=0, N_ele-1 DO y_image[*, j, *] = 0;j*3.
  
  z_image = DBLARR(N_ele, N_ele, N_ele)
  FOR k=0, N_ele-1 DO z_image[*, *, k] = 0.7*(N_ele-k)
  
  new_image = x_image + y_image + z_image
  iimage, REFORM(new_image[*,*, 127])
  iimage, REFORM(new_image[127,*, *])
  iimage, REFORM(new_image[*, 127, *])


END

PRO radon_doc, input_array, backproject1

  array = CONGRID(input_array, 128, 128)

  DEVICE, DECOMPOSED=0

  ;Create an image with a ring plus random noise:

;  x = (LINDGEN(128,128) MOD 128) - 63.5
;
;  y = (LINDGEN(128,128)/128) - 63.5
;
;  radius = SQRT(x^2 + y^2)
;
;  array = (radius GT 40) AND (radius LT 50)
;
;  array = array + RANDOMU(seed,128,128)

  ;Create display window, set graphics properties:

  WINDOW, XSIZE=440,YSIZE=700, TITLE='Magnetic Particle Tomography'

  !P.BACKGROUND = 255 ; white

  !P.COLOR = 0 ; black

  !P.FONT=2
  
  LOADCT, 27
  
  ERASE

  XYOUTS, .05, .94, 'Original Image', /NORMAL

  ;Display the image. 255b changes black values to white:

  TVSCL, array, .05, .75, /NORMAL

  ;Calculate and display the Radon transform:

  XYOUTS, .05, .70, 'Multiple-angle Projections', /NORMAL

  result = RADON(array, RHO=rho, THETA=theta)

  TVSCL, 255b - result, .08, .32, /NORMAL

  PLOT, theta, rho, /NODATA, /NOERASE, $

    POSITION=[0.08,0.32, 1, 0.68], $

    XSTYLE=9,YSTYLE=9,XTITLE='Theta', YTITLE='R'

  ;Find the Radon backprojection and display the output:

  XYOUTS, .05, .21, 'MPT Image', /NORMAL

  ;backproject = RADON(result, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject = CONGRID(backproject1, 128, 128)

  TVSCL, backproject, .05, .02, /NORMAL
  
  DEVICE, DECOMPOSED=1

END

PRO radon_test_angle

  DEVICE, DECOMPOSED=0

  ;Create an image with a ring plus random noise:

    x = (LINDGEN(256,256) MOD 256) - 127.5
  
    y = (LINDGEN(256,256)/256) - 127.5
  
    radius = SQRT(x^2 + y^2)
  
    plate = (radius LT 30); AND (radius LT 100)
    array = BYTSCL(SHIFT(plate, 60, 0) + SHIFT(plate, -60, 0))
;    plate1 = (radius LT 50)
    
  
;    array = array + RANDOMU(seed,256,256)
;    iimage, array, title='original'
  ;Create display window, set graphics properties:
;  result = RADON(array, RHO=rho, THETA=theta);, NRHO=256, NTHETA=256)
  
  NRHO = 137;;137
  NTHETA = 181
  array = CONGRID(array, NRHO/2, NRHO/2)
  FOV = 190.;
  n_x = NRHO/2
  n_y = NRHO/2
  d_x = FOV/n_x
  d_y = FOV/n_y
  iimage, array, title='original'
;  new_result = CONGRID(result, NTHETA, NRHO)
;  new_rho = CONGRID(rho, NRHO)
;  new_theta = CONGRID(theta, NTHETA)
  new_result = RADON(array, RHO=rho, THETA=theta, NRHO=NRHO, NTHETA=NTHETA)
  iimage, new_result, title='Sinogram'
  
  ;@GJ, 2025/2/1, do the FFL-based signal simulation
  sinogram_R = new_result * 0.
  sinogram_I = new_result * 0.
  sinogram_STED = new_result * 0.
  G3RI_equation_calculation, x_array, H_x_array, G_3R_psf, G_3I_psf, G_3STED_psf
  FOR i=0, NTHETA-1 DO BEGIN
    sinogram_R[i, *] = CONVOL(REFORM(new_result[i, *]), G_3R_psf[100:150])
    sinogram_I[i, *] = CONVOL(REFORM(new_result[i, *]), G_3I_psf[100:150])
    sinogram_STED[i, *] = CONVOL(REFORM(new_result[i, *]), G_3STED_psf[100:150])
  ENDFOR
  
  backproject = RADON(new_result, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_R = RADON(sinogram_R, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_I = RADON(sinogram_I, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_STED = RADON(sinogram_STED, /BACKPROJECT, RHO=rho, THETA=theta)
;  iimage, backproject, title='BP'
  
  filtered_result = new_result * 0.
  filtered_resultR = new_result * 0.
  filtered_resultI = new_result * 0.
  filtered_resultSTED = new_result * 0.
  FOR i=0, NTHETA-1 DO BEGIN
    filtered_result[i, *] = Filter_FFT(new_result[i, *])
    filtered_resultR[i, *] = Filter_FFT(sinogram_R[i, *])
    filtered_resultI[i, *] = Filter_FFT(sinogram_I[i, *])
    filtered_resultSTED[i, *] = Filter_FFT(sinogram_STED[i, *])
  ENDFOR
  filtered_backproject = RADON(filtered_result, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectR = RADON(filtered_resultR, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectI = RADON(filtered_resultI, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectSTED = RADON(filtered_resultSTED, /BACKPROJECT, RHO=rho, THETA=theta)
;  iimage, filtered_backproject, title='filtered BP'
  
  DEVICE, DECOMPOSED = 0
  LOADCT, 0
  WINDOW, 0, XSIZE = 4*filtered_result[0], YSIZE = 2*filtered_result[1], TITLE = 'Sinogram'
  TV, BYTSCL(new_result), 0
  TV, BYTSCL(sinogram_R), 1
  TV, BYTSCL(sinogram_I), 2
  TV, BYTSCL(sinogram_STED), 3
  TV, BYTSCL(filtered_result), 4
  TV, BYTSCL(filtered_resultR), 5
  TV, BYTSCL(filtered_resultI), 6
  TV, BYTSCL(filtered_resultSTED), 7
  
  DEVICE, DECOMPOSED = 0
  LOADCT, 0
  WINDOW, 1, XSIZE = 4*filtered_backproject[0], YSIZE = 2*filtered_backproject[1], TITLE = 'FBP'
  TV, BYTSCL(backproject), 0
  TV, BYTSCL(backproject_R), 1
  TV, BYTSCL(backproject_I), 2
  TV, BYTSCL(backproject_STED), 3
  TV, BYTSCL(filtered_backproject), 4
  TV, BYTSCL(filtered_backprojectR), 5
  TV, BYTSCL(filtered_backprojectI), 6
  TV, BYTSCL(filtered_backprojectSTED), 7
END

;;2021/10/31, https://idlastro.gsfc.nasa.gov/ftp/pro/image/correl_images.pro
;;Read an image:
;image = BYTSCL(DIST(400))
;; Add noise to the image:
;image_noisy = NOISE_HURL(image, 0.2)
;a=correl_images(image, image_noisy, XSHIFT=0, YSHIFT=0)
function correl_images, image_A, image_B, XSHIFT = x_shift, $
  YSHIFT = y_shift,   $
  XOFFSET_B = x_offset, $
  YOFFSET_B = y_offset, $
  REDUCTION = reducf, $
  MAGNIFICATION = Magf, $
  NUMPIX=numpix, MONITOR=monitor
  ;+
  ; NAME:
  ; CORREL_IMAGES
  ; PURPOSE:
  ;       Compute the 2-D cross-correlation function of two images
  ; EXPLANATION:
  ;       Computes the 2-D cross-correlation function of two images for
  ;       a range of (x,y) shifting by pixels of one image relative to the other.
  ;
  ; CALLING SEQUENCE:
  ;       Result = CORREL_IMAGES( image_A, image_B,
  ;                        [XSHIFT=, YSHIFT=, XOFFSET_B=, YOFFSET_B=, REDUCTION=,
  ;                        MAGNIFICATION=, /NUMPIX, /MONITOR  )
  ;
  ; INPUTS:
  ;       image_A, image_B = the two images of interest.
  ;
  ; OPTIONAL INPUT KEYWORDS:
  ;       XSHIFT = the + & - shift to be applied in X direction, default=7.
  ;       YSHIFT = the Y direction + & - shifting, default=7.
  ;
  ;       XOFFSET_B = initial X pixel offset of image_B relative to image_A.
  ;       YOFFSET_B = Y pixel offset, defaults are (0,0).
  ;
  ;       REDUCTION = optional reduction factor causes computation of
  ;                       Low resolution correlation of bin averaged images,
  ;                       thus faster. Can be used to get approximate optimal
  ;                       (x,y) offset of images, and then called for successive
  ;                       lower reductions in conjunction with CorrMat_Analyze
  ;                       until REDUCTION=1, getting offset up to single pixel.
  ;
  ;       MAGNIFICATION = option causes computation of high resolution correlation
  ;                       of magnified images, thus much slower.
  ;                       Shifting distance is automatically = 2 + Magnification,
  ;                       and optimal pixel offset should be known and specified.
  ;                       Optimal offset can then be found to fractional pixels
  ;                       using CorrMat_Analyze( correl_images( ) ).
  ;
  ;       /NUMPIX - if set, causes the number of pixels for each correlation
  ;                       to be saved in a second image, concatenated to the
  ;                       correlation image, so Result is fltarr( Nx, Ny, 2 ).
  ;       /MONITOR causes the progress of computation to be briefly printed.
  ;
  ; OUTPUTS:
  ;       Result is the cross-correlation function, given as a matrix.
  ;
  ; PROCEDURE:
  ;       Loop over all possible (x,y) shifts, compute overlap and correlation
  ;       for each shift. Correlation set to zero when there is no overlap.
  ;
  ; MODIFICATION HISTORY:
  ;       Written, July,1991, Frank Varosi, STX @ NASA/GSFC
  ;       Use ROUND instead of NINT, June 1995, Wayne Landsman HSTX
  ;       Avoid divide by zero errors, W. Landsman HSTX April 1996
  ; Remove use of !DEBUG    W. Landsman   June 1997
  ;       Subtract mean of entire image before computing correlation, not just
  ;          mean of overlap region   H. Ebeling/W. Landsman   June 1998
  ;       Always REBIN() using floating pt arithmetic W. Landsman  Nov 2007
  ;
  ;-
  compile_opt idl2
  if N_params() LT 2 then begin
    print,'Syntax  -  Result = CORREL_IMAGES( image_A, image_B,'
    print,'[         XSHIFT=, YSHIFT=, XOFFSET_B=, YOFFSET_B=, REDUCTION=, '
    print,'          MAGNIFICATION=, /NUMPIX, /MONITOR  )'
    return,-1
  endif

  simA = size( image_A )
  simB = size( image_B )
  do_int = (simA[3] LE 3) or (simA[3] GE 12) or $
    (simB[3] LE 3) or (simB[3] GE 12)

  if (simA[0] LT 2) OR (simB[0] LT 2) then begin
    message,"first two arguments must be images",/INFO,/CONTIN
    return,[-1]
  endif

  if N_elements( x_offset ) NE 1 then x_offset=0
  if N_elements( y_offset ) NE 1 then y_offset=0

  if N_elements( x_shift ) NE 1 then x_shift = 7
  if N_elements( y_shift ) NE 1 then y_shift = 7
  x_shift = abs( x_shift )
  y_shift = abs( y_shift )

  if keyword_set( reducf ) then begin

    reducf = fix( reducf ) > 1
    if keyword_set( monitor ) then $
      print,"Reduction = ",strtrim( reducf, 2 )
    simA = simA/reducf
    LA = simA * reducf -1 ;may have to drop edges of images.
    simB = simB/reducf
    LB = simB * reducf -1

    if do_int then begin

      imtmp_A = Rebin( float( image_A[ 0:LA[1], 0:LA[2] ]),  $
        simA[1], simA[2] )
      imtmp_B = Rebin( float( image_B[ 0:LB[1], 0:LB[2] ]),  $
        simB[1], simB[2] )
    endif else begin
      imtmp_A =Rebin( image_A[ 0:LA[1], 0:LA[2] ], simA[1], simA[2] )
      imtmp_B =Rebin( image_B[ 0:LB[1], 0:LB[2] ], simB[1], simB[2] )
    endelse

    xoff = round ( x_offset/reducf )
    yoff = round ( y_offset/reducf )
    xs = x_shift/reducf
    ys = y_shift/reducf

    return, correl_images( imtmp_A, imtmp_B, XS=xs,YS=ys,$
      XOFF=xoff, YOFF=yoff, $
      MONITOR=monitor, NUMPIX=numpix )

  endif else if keyword_set( Magf ) then begin

    Magf = fix( Magf ) > 1
    if keyword_set( monitor ) then $
      print,"Magnification = ",strtrim( Magf, 2 )
    simA = simA*Magf
    simB = simB*Magf

    imtmp_A = rebin( image_A, simA[1], simA[2], /SAMPLE )
    imtmp_B = rebin( image_B, simB[1], simB[2], /SAMPLE )

    xoff = round( x_offset*Magf )
    yoff = round( y_offset*Magf )

    return, correl_images( imtmp_A, imtmp_B, XS=Magf+2, YS=Magf+2,$
      XOFF=xoff, YOFF=yoff, $
      MONITOR=monitor, NUMPIX=numpix )
  endif

  Nx = 2 * x_shift + 1
  Ny = 2 * y_shift + 1
  if keyword_set( numpix ) then Nim=2 else Nim=1

  correl_mat = fltarr( Nx, Ny, Nim )

  xs = round( x_offset ) - x_shift
  ys = round( y_offset ) - y_shift

  sAx = simA[1]-1
  sAy = simA[2]-1
  sBx = simB[1]-1
  sBy = simB[2]-1
  meanA = total( image_A )/(simA[1]*simA[2])
  meanB = total( image_B )/(simB[1]*simB[2])

  for y = 0, Ny-1 do begin  ;compute correlation for each y,x shift.

    yoff = ys + y
    yAmin = yoff > 0
    yAmax = sAy < (sBy + yoff)
    yBmin = (-yoff) > 0
    yBmax = sBy < (sAy - yoff)    ;Y overlap

    if (yAmax GT yAmin) then begin

      for x = 0, Nx-1 do begin

        xoff = xs + x
        xAmin = xoff > 0
        xAmax = sAx < (sBx + xoff)
        xBmin = (-xoff) > 0
        xBmax = sBx < (sAx - xoff)   ;X overlap

        if (xAmax GT xAmin) then begin

          im_ov_A = image_A[ xAmin:xAmax, yAmin:yAmax ]
          im_ov_B = image_B[ xBmin:xBmax, yBmin:yBmax ]
          Npix = N_elements( im_ov_A )

          if N_elements( im_ov_B ) NE Npix then begin
            message,"overlap error: # pixels NE",/INFO,/CONT
            print, Npix, N_elements( im_ov_B )
          endif

          im_ov_A = im_ov_A - meanA
          im_ov_B = im_ov_B - meanB
          totAA = total( im_ov_A * im_ov_A )
          totBB = total( im_ov_B * im_ov_B )

          if (totAA EQ 0) or (totBB EQ 0) then $
            correl_mat[x,y] = 0.0 else $
            correl_mat[x,y] = total( im_ov_A * im_ov_B ) / $
            sqrt( totAA * totBB )

          if keyword_set( numpix ) then correl_mat[x,y,1] = Npix
        endif

      endfor
    endif

    if keyword_set( monitor ) then print, Ny-y, FORM="($,i3)"
  endfor

  if keyword_set( monitor ) then print," "

  return, correl_mat
end

;@GJ, 2024/1/18, read zhongwei'data 
PRO read_ZW_FM_data
  
  filename = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\8mT2-9.txt'
  HAC = 8.;
;  filename = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\14.2mT2-9.txt'
; HAC = 14.2;
  ;read the first line and get the size title
  OPENR, lun, filename, /GET_LUN
  data_string=''
  READF, lun, data_string, FORMAT='(%"%s")'
  ;close reading
  FREE_LUN, lun
  temp_data=DOUBLE(STRSPLIT(data_string, STRING(9B), /extract))
  z_gradient = 1.717;3.434
  x_gradient = 1.717
  z_x_ratio = z_gradient / x_gradient
  n_ele = SQRT(N_ELEMENTS(temp_data))
  x_n_pixels = 256
  z_n_pixels = x_n_pixels * z_x_ratio
  harmonic_maps = DBLARR(8, x_n_pixels, x_n_pixels * z_x_ratio)
  OPENR, lun, filename, /GET_LUN
  FOR i=0, 7 DO BEGIN
    IF EOF(lun) EQ 0 THEN BEGIN
      data_string=''
      READF, lun, data_string, FORMAT='(%"%s")'
      temp_data=DOUBLE(STRSPLIT(data_string, STRING(9B), /extract))
      harmonic_maps[i, *, *] = CONGRID(REFORM(temp_data, n_ele, n_ele), x_n_pixels, z_n_pixels, /INTERP)
    ENDIF
  ENDFOR
  ;close reading
  FREE_LUN, lun
  loadct, 8, RGB_TABLE=RBG_value_8
  iimage, REFORM(harmonic_maps[0, *, *]), title='A2 (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_8
  iimage, ALOG(REFORM(harmonic_maps[0, *, *])), title='ALOG(A2) (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_8
  loadct, 3, RGB_TABLE=RBG_value_3
  iimage, REFORM(harmonic_maps[1, *, *]), title='A3 (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_3
  iimage, ALOG(REFORM(harmonic_maps[1, *, *])), title='ALOG(A3) (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_3

  ;find the max location
  A3_max = MAX(REFORM(harmonic_maps[1, *, *]), max_loc)
  max_loc_2d = ARRAY_INDICES(REFORM(harmonic_maps[1, *, *]), max_loc)

  FOV = 20. ;mm
  z_field_range = z_gradient * FOV
  x_field_range = x_gradient * FOV
  z_field_per_pixel = z_field_range / z_n_pixels
  x_field_per_pixel = x_field_range / x_n_pixels
  x_field = (FINDGEN(x_n_pixels) - max_loc_2d[0]) * x_field_per_pixel
  z_field = (FINDGEN(z_n_pixels) - max_loc_2d[1]) * z_field_per_pixel
  z_pixel_mm = z_field / z_gradient
  x_pixel_mm = x_field / x_gradient
  
  Donut_diameter_array = DBLARR(8, x_n_pixels)
  harmonic_map_w_Donut = harmonic_maps
  FOR i=0, 7 DO BEGIN
    FOR j=0, x_n_pixels-1 DO BEGIN
      temp = REFORM(harmonic_maps[i, j, *])
      max_value = MAX(temp, max_ind)
      Donut_diameter_array[i, j] = 2.*ABS(z_pixel_mm[max_ind] - z_pixel_mm[max_loc_2d[1]])
      harmonic_map_w_Donut[i, j, max_ind] = MAX(harmonic_maps)
      diff_ind = max_ind - max_loc_2d[1]
      IF diff_ind GT 0 AND ABS(max_ind-2*diff_ind) LT z_n_pixels-1 THEN harmonic_map_w_Donut[i, j, max_ind-2*diff_ind] = MAX(harmonic_maps)
      IF diff_ind LT 0 AND ABS(max_ind-2*diff_ind) LT z_n_pixels-1 THEN harmonic_map_w_Donut[i, j, max_ind-2*diff_ind] = MAX(harmonic_maps)
    ENDFOR
  ENDFOR
  
  iplot, x_field/HAC*100., REFORM(Donut_diameter_array[0, *]), color='green', thick=2, xtitle='HDC/HAC [%]', ytitle='Donut Diameter [mm]', title='A2 Donut (HAC='+STRING(HAC, format='(f4.1)')+' mT)'
  iplot, x_field/HAC*100., REFORM(Donut_diameter_array[1, *]), color='red', thick=2, xtitle='HDC/HAC [%]', ytitle='Donut Diameter [mm]', title='A3 Donut (HAC='+STRING(HAC, format='(f4.1)')+' mT)'
  iimage, ALOG(REFORM(harmonic_map_w_Donut[0, *, *])), title='ALOG(A2) with Donut (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_8
  iimage, ALOG(REFORM(harmonic_map_w_Donut[1, *, *])), title='ALOG(A3) with Donut (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_3
 
  ;50 loc
;  x_50_loc = max_loc_2d[0] + HAC*0.5/x_field_per_pixel
;  iplot, z_pixel_mm, REFORM(harmonic_maps[1, x_50_loc, *])
;  
;  ;100 loc
;  x_110_loc = max_loc_2d[0] + HAC*1.1/x_field_per_pixel
;  iplot, z_pixel_mm, REFORM(harmonic_maps[1, x_110_loc, *]), /overplot
;  
;  iplot, z_pixel_mm, REFORM(harmonic_maps[1, 200, *]), color='red', /overplot
  
  IF HAC GT 12 THEN x_ind = 236 ELSE x_ind = 200
  temp = BYTSCL(REFORM(harmonic_maps[0, *, *]))
  temp[x_ind, *] = 255.
  iimage, temp, title='A2 (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_8
  iplot, z_pixel_mm, REFORM(harmonic_maps[0, x_ind, *]), color='green', thick=2, xtitle='z direction [mm]', ytitle='Intensity [a.u.]', title='A2 (HAC='+STRING(HAC, format='(f4.1)')+' mT)'
  
  IF HAC GT 12 THEN x_ind = 248 ELSE x_ind = 220
  temp = BYTSCL(REFORM(harmonic_maps[1, *, *]))
  temp[x_ind, *] = 255.
  iimage, temp, title='A3 (HAC='+STRING(HAC, format='(f4.1)')+' mT)', RGB_TABLE=RBG_value_3
  iplot, z_pixel_mm, REFORM(harmonic_maps[1, x_ind, *]), color='red', thick=2, xtitle='z direction [mm]', ytitle='Intensity [a.u.]', title='A3 (HAC='+STRING(HAC, format='(f4.1)')+' mT)'

END


;GJ, 2021/11/25
;read the results & plot the data
PRO read_data_resolution_vs_H_test
  cd, current=old_dir
  filename = DIALOG_PICKFILE(FILTER = '*.xls', title='Please select your saved data', PATH=old_dir)
  
  ;define the matrix
  N_flat = 10;
  N_particle_size = 10;
  flat_array = DBLARR(N_flat)*0.
  N_particle_size = 10
  particle_size_array = DBLARR(N_particle_size)*0.

  ;define the quality eavluation array
  FWHM_array = DBLARR(N_flat, N_particle_size)*0.
  correlation_array = DBLARR(N_flat, N_particle_size)*0.
  rmse_array = DBLARR(N_flat, N_particle_size)*0.
  ;svd_index_array = DBLARR(N_flat, N_particle_size)*0.
  correlation_image_array = DBLARR(N_flat, N_particle_size)*0.
  rmse_image_array = DBLARR(N_flat, N_particle_size)*0.
  
  ;read the first line and get the size title
  OPENR, lun, filename, /GET_LUN
  title=''
  READF, lun, title, FORMAT='(%"%s")'
  temp_size = STRMID(title, 5, 5, /REVERSE_OFFSET)
  size_title = STRTRIM(temp_size, 1)
  
  ;read the data
  data_string = ''
  READF, lun, data_string, FORMAT='(%"%s")'
  temp_data=STRSPLIT(data_string, ' ', /extract)
  particle_size_array = DOUBLE(temp_data[1:*])
  FOR i=0L, N_flat-1 DO BEGIN
    READF, lun, data_string, FORMAT='(%"%s")'
    temp_data=STRSPLIT(data_string, ' ', /extract)
    flat_array[i] = DOUBLE(temp_data[0])
    FWHM_array[i, *] = DOUBLE(temp_data[1:*])
  ENDFOR
  
  ;read the data
  READF, lun, title, FORMAT='(%"%s")'
  data_string = ''
  READF, lun, data_string, FORMAT='(%"%s")'
  temp_data=STRSPLIT(data_string, ' ', /extract)
  particle_size_array = DOUBLE(temp_data[1:*])
  FOR i=0L, N_flat-1 DO BEGIN
    READF, lun, data_string, FORMAT='(%"%s")'
    temp_data=STRSPLIT(data_string, ' ', /extract)
    correlation_array[i, *] = DOUBLE(temp_data[1:*])
  ENDFOR
  
  ;read the data
  READF, lun, title, FORMAT='(%"%s")'
  data_string = ''
  READF, lun, data_string, FORMAT='(%"%s")'
  temp_data=STRSPLIT(data_string, ' ', /extract)
  particle_size_array = DOUBLE(temp_data[1:*])
  FOR i=0L, N_flat-1 DO BEGIN
    READF, lun, data_string, FORMAT='(%"%s")'
    temp_data=STRSPLIT(data_string, ' ', /extract)
    rmse_array[i, *] = DOUBLE(temp_data[1:*])
  ENDFOR
  
  ;if there is not the end, continue reading the data
  IF EOF(lun) EQ 0 THEN BEGIN
    ;there is image reconstruction
    image_plot = 1
    
    ;read the data
    READF, lun, title, FORMAT='(%"%s")'
    data_string = ''
    READF, lun, data_string, FORMAT='(%"%s")'
    temp_data=STRSPLIT(data_string, ' ', /extract)
    particle_size_array = DOUBLE(temp_data[1:*])
    FOR i=0L, N_flat-1 DO BEGIN
      READF, lun, data_string, FORMAT='(%"%s")'
      temp_data=STRSPLIT(data_string, ' ', /extract)
      correlation_image_array[i, *] = DOUBLE(temp_data[1:*])
    ENDFOR
    
    ;read the data
    READF, lun, title, FORMAT='(%"%s")'
    data_string = ''
    READF, lun, data_string, FORMAT='(%"%s")'
    temp_data=STRSPLIT(data_string, ' ', /extract)
    particle_size_array = DOUBLE(temp_data[1:*])
    FOR i=0L, N_flat-1 DO BEGIN
      READF, lun, data_string, FORMAT='(%"%s")'
      temp_data=STRSPLIT(data_string, ' ', /extract)
      rmse_image_array[i, *] = DOUBLE(temp_data[1:*])
    ENDFOR
  ENDIF ELSE BEGIN
    ;there is no reconstructed image
    image_plot = 0
  ENDELSE
  
  ;close reading
  FREE_LUN, lun
  
  ;plot the data
  cgsurface, FWHM_array, flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=6, /Brewer, TITLE='FWHM (Bar Phantom '+size_title+')', $
    XTITLE = 'Flat portion / %', YTITLE = 'Size / nm', ZTITLE = 'FWHM / mm', ZRANGE=[MIN(FWHM_array)-0.2,MAX(FWHM_array)+0.2];,/CONSTRAIN_ASPECT

  cgsurface, correlation_array, flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=2, /Brewer, TITLE='Correlation (Bar Phantom '+size_title+')', $
    XTITLE = 'Flat portion / %', YTITLE = 'Size / nm', ZTITLE = 'Correlation', ZRANGE=[MIN(correlation_array)-0.04,MAX(correlation_array)+0.04];,/CONSTRAIN_ASPECT

  cgsurface, ALOG(rmse_array), flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='ALOG RMSE (Bar Phantom '+size_title+')', $
    XTITLE = 'Flat portion / %', YTITLE = 'Size / nm', ZTITLE = 'ALOG RMSE', ZRANGE=[MIN(ALOG(rmse_array))-0.04,MAX(ALOG(rmse_array))+0.04];,/CONSTRAIN_ASPECT

  IF image_plot EQ 1 THEN BEGIN
    cgsurface, correlation_image_array, flat_array, particle_size_array, $
      /Shaded, /Elevation_Shading, CTable=3, /Brewer, TITLE='Correlation (Image '+size_title+')', $
      XTITLE = 'Flat portion / %', YTITLE = 'Size / nm', ZTITLE = 'Correlation', ZRANGE=[MIN(correlation_image_array)-0.04,MAX(correlation_image_array)+0.04];,/CONSTRAIN_ASPECT

    cgsurface, ALOG(rmse_image_array), flat_array, particle_size_array, $
      /Shaded, /Elevation_Shading, CTable=5, /Brewer, TITLE='ALOG RMSE (Image '+size_title+')', $
      XTITLE = 'Flat portion / %', YTITLE = 'Size / nm', ZTITLE = 'ALOG RMSE', ZRANGE=[MIN(ALOG(rmse_image_array))-0.04,MAX(ALOG(rmse_image_array))+0.04];,/CONSTRAIN_ASPECT
  ENDIF

  window, 4
  plot, flat_array, FWHM_array[*, N_particle_size-1], YRANGE=[MIN(FWHM_array)-0.2,MAX(FWHM_array)+0.2], PSYM=0, thick=1, color=0, xtitle='flat portion / %', ytitle='FWHM / mm', title='FWHM (Bar Phantom '+size_title+')', BACKGROUND = 'FFFFFF'x
  FOR k=0, N_particle_size-2 DO oplot, flat_array, FWHM_array[*, k], thick = 1, color=350000*(k+1), PSYM=0
    
END


;GJ, 2021/8/26; check the amplitude;
;GJ, 2021/9/12; check the 3rd harmonic with a flat pulse portion
;GJ, 2021/10/1; fix some bugs with a flat pulse portion
;GJ, 2021/11/8; resolution_vs_H_test, 128, /saveDat is the way to evaluate the phantom
;               resolution_vs_H_test, 128, /image, /saveDat is the way to evaluate the 1d image reconstruction
;GJ, 2021/11/10; resolution_vs_H_test, 128, /image, /radonTrans, /saveDat is the way to evaluate 2d image reconstruction
;                 /radonTrans is for 2D images
;GJ, 2022/4/26, evaluate SVD the singular values for 1D
;GJ, 2022/4/28, calculate the SVD decay for 1D
;example:
;         resolution_vs_H_test, 64, /image, /saveDat, /saveImage
;         resolution_vs_H_test, 64, 180, /image, /radonTrans, /saveDat, /saveImage
;         resolution_vs_H_test, 64, 180, /image, /radonTrans, /rasterScan, /saveDat, /saveImage; both radonTrans and rasterScan are for line by line scan
;         resolution_vs_H_test, 64, 180, /saveDat, /saveImage ;this is for bar phantom
;         resolution_vs_H_test, 64, 180, /radonTrans, /saveDat, /saveImage, /sl2d
;         resolution_vs_H_test, 128, 180, /radonTrans, /saveDat, /saveImage, /sl2d
;         resolution_vs_H_test, 128, 180, /saveDat, /saveImage, /sl2d; for projected 1d 
;         resolution_vs_H_test, 128, 180, /saveDat, /saveImage; for projected 1d
;         resolution_vs_H_test, 64, 180, /radonTrans, /rasterScan, /saveDat, /saveImage, /sl2d; for parallel line of sl2d phantom
PRO resolution_vs_H_test, matrix_size, N_angles, image=image, radonTrans=radonTrans, rasterScan = rasterScan, saveDat=saveDat, plotYes=plotYes, saveImage=saveImage, sl2d=sl2d
  
  ;pre-set the matrix_size and N_angles
  IF N_ELEMENTS(matrix_size) EQ 0 THEN matrix_size = 64;128;;192;128;64;128;192;256; move this out of the sequence
  IF N_ELEMENTS(N_angles) EQ 0 THEN N_angles = 180;;128;;192;128;64;128;192;256; move this out of the sequence
  size_title = STRCOMPRESS(STRING(matrix_size, format='(I6)'), /REMOVE_ALL)
  
  ;define number of flat portion and particle sizes
  N_flat = 10;
  flat_array = (0.5 + FINDGEN(N_flat)/N_flat/2.)*100.;[0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
  N_particle_size = 10
  particle_size_array = 20. + FINDGEN(N_particle_size)*2.;[21.0, 25.0, 27.0, 32.0, 35.0]
  
  IF KEYWORD_SET(image) THEN BEGIN
    cd, current=old_dir
    i=0;
    REPEAT BEGIN
      filename2 = DIALOG_PICKFILE(FILTER = ['*.tif', '*.jpeg', '*.jpg', '*.png', '*.dcm'], title='Please select your image', PATH=old_dir)
      i++
      IF (i GT 4) THEN BEGIN
        test=DIALOG_MESSAGE('Program Cancelled!', /ERROR)
        RETURN;
      ENDIF
    ENDREP UNTIL (STRLEN(filename2) GT 0)
    image_temp2 = read_image(filename2)
    IF STRPOS(filename2, '.dcm') GT 0 THEN BEGIN
      size_image = SIZE(image_temp2, /dim)
      mask_image = IMAGE_THRESHOLD(image_temp2, THRESHOLD=tf, /MAXENTROPY)
      image_temp2 = image_temp2*mask_image
    ENDIF
    IF size(image_temp2, /N_DIMENSIONS) EQ 2 THEN BEGIN
      image_2d = BYTSCL(CONGRID(image_temp2, Matrix_size, Matrix_size));, MAX=130)
    ENDIF ELSE BEGIN
      image_2d = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), Matrix_size, Matrix_size));, MAX=130)
    ENDELSE
    
    ;      image = image * 0.
    ;      image[127, 127] = 256
    ;      image[130:132, 127:129] = 256
    iimage,  image_2d, title='Original'
    image_1d_temp = DBLARR(Matrix_size)*0.
    FOR i = 0, Matrix_size-1 DO image_1d_temp[i] = TOTAL(image_2d[i, *])
    image_1d = BYTSCL(image_1d_temp);*10.
    IF KEYWORD_SET(plotYes) THEN plot, image_1d, title='Projected 1d image'
  ENDIF
  
  ;define the quality eavluation array
  FWHM_array = DBLARR(N_flat, N_particle_size)*0.
  correlation_array = DBLARR(N_flat, N_particle_size)*0.
  rmse_array = DBLARR(N_flat, N_particle_size)*0.
  svd_index_array = DBLARR(N_flat, N_particle_size)*0.
  correlation_image_array = DBLARR(N_flat, N_particle_size)*0.
  rmse_image_array = DBLARR(N_flat, N_particle_size)*0.
  svd_slope_array = DBLARR(N_flat, N_particle_size)*0.
  
  IF KEYWORD_SET(image) THEN BEGIN
    mySurface = SURFACE(correlation_image_array, flat_array, particle_size_array, $
      TITLE='Correlation (Image'+size_title+')')
  ENDIF ELSE BEGIN
    mySurface = SURFACE(FWHM_array, flat_array, particle_size_array, $
      TITLE='FWHM (Bar Phantom'+size_title+')')
  ENDELSE
  ; Overlay the contour data.
  ax = mySurface.AXES
  ax[0].TITLE = 'Flat ratio [%]'
  ax[1].TITLE = 'Size [nm]'
  ax[2].TITLE = 'FWHM [mm]'
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='yellow', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start

  FOR k=N_particle_size-1, 0, -1 DO BEGIN
    ;update the progress bar
    count=(N_particle_size-k)/(TOTAL(N_particle_size))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    
    FOR j=N_flat-1, 0, -1 DO BEGIN
      flat_portion = flat_array[j]/100.;0.4;5;8;96;0.5;0.96;5
      H_value_shape_freq = [5.0, flat_portion, 2.5]; [H field in mT, flat_portion, frequency in kHz]
      ;H_value_shape_freq = [1.0, 0.98, 2.5]
      plotYes = 0;
      tracer_info = [0, particle_size_array[k], plotYes] ;30 nm is the particle size [type, particle size in nm]
      sting_particle_size = STRING(particle_size_array[k], format='(I2)')
      string_flat_portion = STRING(flat_array[j], format='(f4.1)')
      string_time = STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(f10.2)')
      IF KEYWORD_SET(image) THEN BEGIN
        IF KEYWORD_SET(radonTrans) THEN BEGIN
          radonTrans=1
          IF KEYWORD_SET(rasterScan) THEN BEGIN
            rasterScan=1
            FWHM_correl_info = gradientMPI_2d_simulation(H_value_shape_freq, tracer_info, matrix_size, N_angles, image_2d_rec, image_2d, radonTrans, rasterScan);, /plotYes)
          ENDIF ELSE BEGIN
            rasterScan=0
            FWHM_correl_info = gradientMPI_2d_simulation(H_value_shape_freq, tracer_info, matrix_size, N_angles, image_2d_rec, image_2d, radonTrans);, /plotYes)
          ENDELSE
        ENDIF ELSE BEGIN
          FWHM_correl_info = gradientMPI_1d_simulation(H_value_shape_freq, tracer_info, matrix_size, image_1d_tvrd, image_1d, image_1d_rec);, /plotYes)
        ENDELSE
      ENDIF ELSE BEGIN
        IF KEYWORD_SET(sl2d) THEN BEGIN
          IF KEYWORD_SET(radonTrans) THEN BEGIN
            radonTrans=1
            IF KEYWORD_SET(rasterScan) THEN BEGIN
              rasterScan=1
              FWHM_correl_info = gradientMPI_2d_simulation(H_value_shape_freq, tracer_info, matrix_size, N_angles, image_2d_rec, radonTrans, rasterScan)
            ENDIF ELSE BEGIN
              rasterScan=0
              FWHM_correl_info = gradientMPI_2d_simulation(H_value_shape_freq, tracer_info, matrix_size, N_angles, image_2d_rec, radonTrans)
            ENDELSE
          ENDIF ELSE BEGIN
            radonTrans=0
            FWHM_correl_info = gradientMPI_1d_simulation(H_value_shape_freq, tracer_info, matrix_size, image_1d_tvrd)
          ENDELSE
        ENDIF ELSE BEGIN
          bar_phantom_only=1
          plotYes=1
          FWHM_correl_info = gradientMPI_1d_simulation(H_value_shape_freq, tracer_info, matrix_size, image_1d_tvrd, bar_phantom_only, plotYes)
        ENDELSE
      ENDELSE
      
      ;save the images
      IF KEYWORD_SET(saveImage) THEN BEGIN
        IF N_ELEMENTS(filename2) GT 0 THEN BEGIN
          saveImageDir = filename2+'_'+size_title
          FILE_MKDIR, filename2+'_'+size_title
          filename = FILE_BASENAME(filename2)
        ENDIF ELSE BEGIN
          IF KEYWORD_SET(sl2d) THEN BEGIN
            cd, current=old_dir
            saveImageDir = old_dir+'\sl_phantom_'+size_title
            FILE_MKDIR, old_dir+'\sl_phantom_'+size_title
            filename = 'sl_phantom'
          ENDIF ELSE BEGIN
            cd, current=old_dir
            saveImageDir = old_dir+'\bar_phantom_'+size_title
            FILE_MKDIR, old_dir+'\bar_phantom_'+size_title
            filename = 'bar_phantom'
          ENDELSE
          
        ENDELSE
        aSaveImageFile = saveImageDir+'\'+filename+'_'+size_title+'_Size'+sting_particle_size+'_flat'+string_flat_portion+'_'+string_time+'.jpg'
        IF N_ELEMENTS(image_2d_rec) GT 1 THEN WRITE_JPEG, aSaveImageFile, BYTSCL(image_2d_rec)
        IF N_ELEMENTS(image_1d_tvrd) GT 1 THEN WRITE_JPEG, aSaveImageFile, BYTSCL(image_1d_tvrd), TRUE=1
      ENDIF

      ;the first 4 are bar_image, the last 2 are image
      ;[FWHM, max_correl, svd_index, rms_error, correl_image, rms_image_error]
      FWHM_array[j, k] = FWHM_correl_info[0]
      correlation_array[j, k] = FWHM_correl_info[1]
      svd_index_array[j, k] = FWHM_correl_info[2]
      rmse_array[j, k] = FWHM_correl_info[3]
      correlation_image_array[j, k] = FWHM_correl_info[4]
      rmse_image_array[j, k] = FWHM_correl_info[5]
      svd_slope_array[j, k] = FWHM_correl_info[6]
      print, 'k =', k, ', j =', j, ', correlation =', FWHM_correl_info[1],  ', correlation (image) =', FWHM_correl_info[4]
      print, 'particle size = ', particle_size_array[k], 'nm, flat portion = ', flat_portion
      IF KEYWORD_SET(image) THEN BEGIN
        mySurface.SetData, correlation_image_array, flat_array, particle_size_array
      ENDIF ELSE BEGIN
        mySurface.SetData, FWHM_array, flat_array, particle_size_array
      ENDELSE
    ENDFOR
  ENDFOR
  
  
  
  ;destroy the progress bar
  progressbar -> Destroy
  
  IF KEYWORD_SET(saveDat) THEN BEGIN
    ;Select a file to save exported data, the default opened directory is the one you select just now
    ;cd, current=old_dir
    ;aSaveFile=DIALOG_PICKFILE(FILTER = ['*.txt', '*.xls', '*.xlsx'], title='Save Your Output as EXCEL file', PATH=old_dir)
    ;Result = FILE_DIRNAME(filename2)
    IF N_ELEMENTS(filename2) GT 0 THEN BEGIN
      saveImageDir = filename2+'_'+size_title
      FILE_MKDIR, filename2+'_'+size_title
      filename = FILE_BASENAME(filename2)
    ENDIF ELSE BEGIN
      cd, current=old_dir
      saveImageDir = old_dir+'\bar_phantom_'+size_title
      FILE_MKDIR, old_dir+'\bar_phantom_'+size_title
      filename = 'bar_phantom'
    ENDELSE
    ;save the surface plot file and close
    aSaveFile_surface = saveImageDir+'\'+filename+'_'+size_title+'_surface'+STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(f10.2)')+'.jpg'
    mySurface.Save, aSaveFile_surface
    mySurface.close
    
    ;save the results as dat file
    aSaveFile = saveImageDir+'\'+filename+'_'+size_title+'_'+STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(f10.2)')+'.xls'
    ;// open an ASCII file for output
    ;IF STRCMP(STRMID(aSaveFile, 2, 3, /REVERSE_OFFSET), 'txt',/FOLD_CASE)   THEN BEGIN
    IF (STRLEN(aSaveFile) NE 0)  THEN BEGIN
      GET_LUN, U
      openw,U,aSaveFile
      printf, U, 'FWHM bar phantom (Matrix size:'+STRING(matrix_size, format='(I5)')+')'
      printf, U, FORMAT = '(11(F8.3, %"\t"))', [0, particle_size_array]
      FOR i=0L, N_flat-1 DO BEGIN
        printf, U, FORMAT = '(11(F8.3, %"\t"))', [flat_array[i], REFORM(FWHM_array[i, *])]
      ENDFOR

      printf, U, 'Correlation bar phantom (Matrix size:'+STRING(matrix_size, format='(I5)')+')'
      printf, U, FORMAT = '(11(F8.3, %"\t"))', [0, particle_size_array]
      FOR i=0L, N_flat-1 DO BEGIN
        printf, U, FORMAT = '(11(F8.3, %"\t"))', [flat_array[i], REFORM(correlation_array[i, *])]
      ENDFOR

      printf, U, 'RMSE bar phantom (Matrix size:'+STRING(matrix_size, format='(I5)')+')'
      printf, U, FORMAT = '(11(F8.3, %"\t"))', [0, particle_size_array]
      FOR i=0L, N_flat-1 DO BEGIN
        printf, U, FORMAT = '(11(F11.3, %"\t"))', [flat_array[i], REFORM(rmse_array[i, *])]
      ENDFOR

      IF KEYWORD_SET(image) THEN BEGIN
        printf, U, 'Correlation image (Matrix size:'+STRING(matrix_size, format='(I5)')+')'
        printf, U, FORMAT = '(11(F8.3, %"\t"))', [0, particle_size_array]
        FOR i=0L, N_flat-1 DO BEGIN
          printf, U, FORMAT = '(11(F8.3, %"\t"))', [flat_array[i], REFORM(correlation_image_array[i, *])]
        ENDFOR

        printf, U, 'RMSE image (Matrix size:'+STRING(matrix_size, format='(I5)')+')'
        printf, U, FORMAT = '(11(F8.3, %"\t"))', [0, particle_size_array]
        FOR i=0L, N_flat-1 DO BEGIN
          printf, U, FORMAT = '(11(F11.3, %"\t"))', [flat_array[i], REFORM(rmse_image_array[i, *])]
        ENDFOR
      ENDIF

      FREE_LUN, U
    ENDIF
  ENDIF
  
  ;GJ, 2022/4/28, calculate SVD slope
  cgsurface, svd_slope_array, flat_array, particle_size_array, $
    /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='SVD slope (Bar Phantom'+size_title+')', $
    XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'SVD slope', ZRANGE=[MIN(svd_slope_array)-0.2,MAX(svd_slope_array)+0.2];;,/CONSTRAIN_ASPECT

  
  IF KEYWORD_SET(plotYes) THEN BEGIN
    ;plot the data
    cgsurface, FWHM_array, flat_array, particle_size_array, $
      /Shaded, /Elevation_Shading, CTable=6, /Brewer, TITLE='FWHM (Bar Phantom'+size_title+')', $
      XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'FWHM', ZRANGE=[MIN(FWHM_array)-0.2,MAX(FWHM_array)+0.2];,/CONSTRAIN_ASPECT

    cgsurface, correlation_array, flat_array, particle_size_array, $
      /Shaded, /Elevation_Shading, CTable=2, /Brewer, TITLE='Correlation (Bar Phantom'+size_title+')', $
      XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'Correlation', ZRANGE=[MIN(correlation_array)-0.04,MAX(correlation_array)+0.04];,/CONSTRAIN_ASPECT

    cgsurface, ALOG(rmse_array), flat_array, particle_size_array, $
      /Shaded, /Elevation_Shading, CTable=4, /Brewer, TITLE='ALOG RMSE (Bar Phantom'+size_title+')', $
      XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'ALOG RMSE', ZRANGE=[MIN(ALOG(rmse_array))-0.04,MAX(ALOG(rmse_array))+0.04];,/CONSTRAIN_ASPECT
  
  
    IF KEYWORD_SET(image) THEN BEGIN
      cgsurface, correlation_image_array, flat_array, particle_size_array, $
        /Shaded, /Elevation_Shading, CTable=3, /Brewer, TITLE='Correlation (Image'+size_title+')', $
        XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'Correlation', ZRANGE=[MIN(correlation_image_array)-0.04,MAX(correlation_image_array)+0.04];,/CONSTRAIN_ASPECT

      cgsurface, ALOG(rmse_image_array), flat_array, particle_size_array, $
        /Shaded, /Elevation_Shading, CTable=5, /Brewer, TITLE='ALOG RMSE (Image'+size_title+')', $
        XTITLE = 'Flat ratio [%]', YTITLE = 'Size [nm]', ZTITLE = 'ALOG RMSE', ZRANGE=[MIN(ALOG(rmse_image_array))-0.04,MAX(ALOG(rmse_image_array))+0.04];,/CONSTRAIN_ASPECT
    ENDIF

    window, 4
    plot, flat_array, FWHM_array[*, N_particle_size-1], YRANGE=[MIN(FWHM_array)-0.2,MAX(FWHM_array)+0.2], PSYM=0, thick=1, color=0, xtitle='flat portion [%]', ytitle='FWHM', title='FWHM (Bar Phantom'+size_title+')', BACKGROUND = 'FFFFFF'x
    FOR k=0, N_particle_size-2 DO oplot, flat_array, FWHM_array[*, k], thick = 1, color=350000*(k+1), PSYM=0

  ENDIF
  
    ;GJ, 2021/11/6, close the plot
    ;mySurface.close

    
 
END

;GJ, 2021/10/07, do the 1-d simulation
;GJ, 2021/10/29, if phantom is set, the linepairs is generated
;GJ, 2021/11/03, noise level is fixed as 1/10000 of the signal based on literature, IJNM_010_3097_15$MPIFuture%P
;GJ, 2021/11/05, increase the speed by saving the field and signal AUC relationship table
;GJ, 2021/11/08, include the cases for bar_phantom_only and image_only or both
;GJ, 2022/4/26, evaluate SVD the singular values
;GJ, 2022/4/28, calculate the SVD decay
;if outputDat is set, the data is output as dat file
;example: FWHM_info = gradientMPI_1d_simulation(/bar_phantom_only)
;example: FWHM_info = gradientMPI_1d_simulation();calculate the Shepp-Logan phantom image
FUNCTION gradientMPI_1d_simulation, H_value_shape_freq, tracer_info, matrix_size, image_1d_tvrd, image_1d, image_1d_rec, sv_inv, V, U, image_only=image_only, bar_phantom_only=bar_phantom_only, saveDat=saveDat, plotYes=plotYes
  
  ;GJ, 2021/10/29
  IF N_ELEMENTS(matrix_size) EQ 0 THEN Matrix_size = 256;512;128;256;192;256;128;64;512;256;128;64;256;64;128;256;128;64;512;256;512;256;128;64;256;512;256;256;64;64;128;180

  IF keyword_set(bar_phantom_only) THEN BEGIN
    ;image_1d=image_phantom_res(Matrix_size)
    bar_image_1d=image_phantom_res(Matrix_size, freq_array);, /plotYes)
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(image_1d) LT 10 THEN BEGIN 
      ;filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
      ;image_temp = read_image(filename2)
      ;image = BYTSCL(CONGRID((REFORM(image_temp[0,*,*])), Matrix_size, Matrix_size));, MAX=130)
      image = BYTSCL(shepp_logan_phantom_2d(Matrix_size))
      ;      image = image * 0.
      ;      image[127, 127] = 256
      ;      image[130:132, 127:129] = 256
      IF KEYWORD_SET(plotYes) THEN iimage,  image, title='Original'
      image_1d_temp = DBLARR(Matrix_size)*0.
      FOR i = 0, Matrix_size-1 DO image_1d_temp[i] = TOTAL(image[i, *])
      image_1d = BYTSCL(image_1d_temp);*10.
    ENDIF ELSE BEGIN
      Matrix_size = N_ELEMENTS(image_1d)
    ENDELSE
    bar_image_1d=image_phantom_res(Matrix_size, freq_array);, /plotYes)
  ENDELSE
  
  ;set up the encoding magnetic fields
  G = 0.025; T/m
  FOV = 0.2; m
  x_res = FOV/(Matrix_size-1.)
  x_loc_array = (FINDGEN(Matrix_size*2)-Matrix_size+0.5)*x_res ;calculate the total number of locations
  H_array = x_loc_array * G * 1000. ;get the corresponding H based on gradient and location
  IF KEYWORD_SET(plotYes) THEN iplot, x_loc_array*1000., H_array, xtitle='Location [mm]', ytitle='Field Amplitude [mT/u0]', title='H array'
  
  ;calculate the signal from single nanoparticle and signal features
  ;set the scan parameters for calculation
  IF N_ELEMENTS(H_value_shape_freq) EQ 0 THEN BEGIN
    ;H_flat_portion = 0.0
    H_flat_portion = 0.9;0.8;0.9
    H_frequency = 2.5; kHz
    ;H_value_shape_freq = [1.0, 0.98, 2.5]
  ENDIF ELSE BEGIN
    H_flat_portion = H_value_shape_freq[1]
    H_frequency = H_value_shape_freq[2]
  ENDELSE
  
  IF N_ELEMENTS(tracer_info) EQ 0 THEN BEGIN
    plotYes = 0
    tracer_type = 0.
    tracer_size = 28. ; nm
    tracer_info = [tracer_type, tracer_size, plotYes]
  ENDIF

  FOR i=0, Matrix_size*2-1 DO BEGIN
    H_value_shape_freq = [H_array[i], H_flat_portion, H_frequency]; 0.96
    feature = signal_vs_H(H_value_shape_freq, tracer_info, 1, test_signal)
    IF i EQ 0 THEN BEGIN
      size_feature = size(feature, /dim)
      feature_array = DBLARR(Matrix_size*2, size_feature[0], size_feature[1])
      size_signal = size(test_signal, /dim)
      signal_array = DBLARR(Matrix_size*2, size_signal[0], size_signal[1])
    ENDIF
    
    feature_array[i, *, *] = feature
    signal_array[i, *, *] = test_signal
  ENDFOR
  
  ;determine noise level
  noise_level = ABS(MAX(bar_image_1d))/100000.
  ;calculate bar phantom
  ;define the signal kernal
  signalKernal_2d_array = DBLARR(Matrix_size, Matrix_size)*0.
  
  ;if the sv_ind has been determined
  IF N_ELEMENTS(sv_inv) GT 2 THEN GOTO, pos_directRec

  bar_signal_1d = DBLARR(Matrix_size)*0.
  bar_signal_total = DBLARR(size_signal[0], size_signal[1])*0.
  ;  filteredsignal_total = DBLARR(3, 1199)*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      signalKernal_2d_array[i, j] = feature_array[i+j, 1, 1] ;Decay signal AUC for pulsed excitation
      bar_signal_total[0,*] = signal_array[i+j, 0, *] ;this is the time
      ;save the H curve
      ;noise is a random number between -1 and 1
      seed = !NULL
      noise_array = noise_level*RANDOMU(seed, N_ELEMENTS(bar_signal_total[0,*]))
      ;noise_array = RANDOMU(seed, N_ELEMENTS(bar_signal_total[0,*]))*0.0
      test_signal[1, *] = bar_image_1d[j] * test_signal[1, *] + noise_array ;this is the non-adiabatic signal
      test_signal[2, *] = bar_image_1d[j] * test_signal[2, *] + noise_array ;adiabatic signal with relaxation
      bar_signal_total[1, *] = bar_signal_total[1, *] + bar_image_1d[j] * signal_array[i+j, 1, *] + noise_array ;non-adiabatic signal
      bar_signal_total[2, *] = bar_signal_total[2, *] + bar_image_1d[j] * signal_array[i+j, 2, *] + noise_array ;adiabatic signal with relaxation
      bar_signal_total[3, *] = bar_signal_total[3, *] + signal_array[i+j, 3, *] ;magentic field
    ENDFOR
    IF MAX(bar_signal_total[3,*]) GT 0.001 THEN BEGIN
      IF KEYWORD_SET(plotYes) THEN BEGIN
        window, 2
;        plot, bar_signal_total[0,*], bar_signal_total[3,*]/MAX(bar_signal_total[3,*])*MAX(bar_signal_total[1,*]), BACKGROUND = 'FFFFFF'x, COLOR = 0, XTITLE='time [us]', YTITLE='Signal [A.U.]', title=STRING(i) + ' encoding step'
;        oplot, bar_signal_total[0,*], bar_signal_total[1,*], LINESTYLE=3, thick = 2, color=35000
;        oplot, bar_signal_total[0,*], bar_signal_total[2,*], thick = 2, color=3500;/MAX(filteredsignal)*MAX(signal)
        H_array = bar_signal_total[3,*]
        signal_array2 = bar_signal_total[1,*]
        signal_array1 = bar_signal_total[2,*]
        cgplot, bar_signal_total[0,*], H_array/MAX(H_array)*MAX(signal_array2), LINESTYLE = 2, thick = 2, title='signal', xtitle='time [us]', ytitle='Signal \ A.U.', BACKGROUND = 'FFFFFF'x, COLOR = 0
        cgplot, bar_signal_total[0,*], signal_array2, thick = 2, LINESTYLE = 0, color='blue',/overplot
        cgplot, bar_signal_total[0,*], signal_array1, thick = 2, color='red',/overplot
      ENDIF
    ENDIF
    ;wait, 0.5
    ;calculate the AUC
    IF ABS(H_value_shape_freq[0]) LT 0.01 THEN H_value_shape_freq[0] = 5.
    bar_signal_1d[i] = (signal_vs_H(H_value_shape_freq,tracer_info,0, bar_signal_total))[1, 1];get the AUC
    bar_signal_total = bar_signal_total*0.
  ENDFOR
  
  ; Compute the Singular Value Decomposition:
  SVDC, signalKernal_2d_array, W, U, V
  
  IF KEYWORD_SET(plotYes) THEN iplot, w, thick=2, title='SVD of System Matrix', xtitle='Component Index', ytitle='Singular Value', /ylog, xrange=[1,N_ELEMENTS(w)]
;  iplot, alog(w), title='SVD Analysis', xtitle='Index', ytitle='Alog W'
  
  poly_fit_result = POLY_FIT(FINDGEN(N_ELEMENTS(w))+1, ALOG(w), 10,sigma=sigma, yfit=yw)
  svd_slope = ABS(poly_fit_result[1])

  IF KEYWORD_SET(plotYes) THEN BEGIN
    plot, alog(w), title='10th order ploy fit'
    oplot, yw
  END
;  plot, alog(w), title='10th order ploy fit'
;  oplot, yw
  
  ;plot, ABS(alog(w)-yw), title='5th order ploy fit Error'
  max_error = MAX(ABS(alog(w)-yw), maxind_polyfit)
  svd_cutoff = maxind_polyfit
  print, 'svd_cutoff: ', svd_cutoff
;  
;  x_test = FINDGEN(210)+1
;  y_test = ALOG(w[0:209])
;  weights = y_test
;  A = [MAX(y_test),-0.1,MEAN(y_test[200:209])]
;  ;Compute the parameters.
;  ywfit = CURVEFIT(x_test, y_test, weights, A, SIGMA, FUNCTION_NAME='gfunct')
;  plot, x_test, y_test
;  oplot, x_test, ywfit
;  
  ;do loop iterative reconstruction for convergence reconstruction
  sv_inv = FLTARR(Matrix_size, Matrix_size)*0.
  bar_y_rec = FLTARR(Matrix_size, Matrix_size)*0.
  bar_y0 = DBLARR(Matrix_size)*0.
  ;error_dist = DBLARR(Matrix_size)*0.
  correl_array = DBLARR(Matrix_size)*0.
  FOR K = 0, Matrix_size-1 DO BEGIN
    sv_inv[K,K] = 1./W[K]
    bar_y = V ## sv_inv ## TRANSPOSE(U) ## bar_signal_1d
    ;bar_y[WHERE(bar_y LT 0, /NULL)] = 0
    ;bar_y[WHERE(bar_y GT 255, /NULL)] = 255
    bar_y_rec[K, *] = REFORM(bar_y)
    correl_array[K] = CORRELATE(bar_y, bar_image_1d)
    ;error_dist[k-1] = TOTAL(ABS(y-y0))/Matrix_size
  ENDFOR
  IF KEYWORD_SET(plotYes) THEN BEGIN
    iimage, bar_y_rec, title='Reconstrcted Bar Phantom'
    window, 3
    plot, correl_array, title='bar phantom correlation'
  ENDIF

  max_correl = MAX(correl_array[FLOOR(Matrix_size*2/4)-1:Matrix_size-1], svd_index)
;  ;GJ, 2021/11/8
  svd_index = FLOOR(Matrix_size*2/4)-1 + svd_index
;  svd_index = Matrix_size-1
  ;print, 'SVD index (bar phantom) = ', svd_index
  ;print, 'correlation (bar phantom) = ', max_correl
  ;remove the higher components
  IF svd_index LT Matrix_size-1 THEN sv_inv[svd_index+1:Matrix_size-1] = 0.
  bar_image_rec = REFORM(bar_y_rec[svd_index, *])
  ;calculate rms_error
  rms_error = Sqrt(Total((bar_image_rec - bar_image_1d)^2)/N_Elements(bar_image_1d))
  IF KEYWORD_SET(plotYes) THEN BEGIN
    window, 5
    plot, bar_image_rec, thick=1, title='SVD reconstruction'+STRING(Matrix_size), color=0, BACKGROUND = 'FFFFFF'x
    oplot, bar_image_1d, color=3500000, THICK=2, LINESTYLE=2
    ;do MTF calculation
  ENDIF
;  window, 5
;  plot, bar_image_rec, thick=1, title='SVD reconstruction'+STRING(Matrix_size), color=0, BACKGROUND = 'FFFFFF'x
;  oplot, bar_image_1d, color=3500000, THICK=2, LINESTYLE=2
  image_1d_tvrd = TVRD(TRUE=1)
  
  IF N_ELEMENTS(freq_array) GE 1 THEN FWHM = MTF_cal(bar_image_1d, bar_image_rec, freq_array, FOV)
  
  pos_directRec:
    ;if bar_phantom_only is set, calculate bar phantom's MTF
  ;if bar_phantom_only is not set, reconstruct the image
  IF keyword_set(bar_phantom_only) EQ 0 THEN BEGIN
    ;define the signal kernal
    signal_1d = DBLARR(Matrix_size)*0.
    signal_total = DBLARR(size_signal[0], size_signal[1])*0.
    FOR i=0, Matrix_size-1 DO BEGIN
      FOR j=0, Matrix_size-1 DO BEGIN
        signal_total[0,*] = signal_array[i+j, 0, *] ;this is the time
        ;save the H curve
        ;noise is a random number between -1 and 1
        seed = !NULL
        noise_array = noise_level*RANDOMU(seed, size_signal[1])
        test_signal[1, *] = image_1d[j] * test_signal[1, *] + noise_array ;this is the non-adiabatic signal
        test_signal[2, *] = image_1d[j] * test_signal[2, *] + noise_array ;adiabatic signal with relaxation
        signal_total[1, *] = signal_total[1, *] + image_1d[j] * signal_array[i+j, 1, *] + noise_array ;non-adiabatic signal
        signal_total[2, *] = signal_total[2, *] + image_1d[j] * signal_array[i+j, 2, *] + noise_array ;adiabatic signal with relaxation
        signal_total[3, *] = signal_total[3, *] + signal_array[i+j, 3, *] ;magentic field
      ENDFOR
      IF MAX(signal_total[3,*]) GT 0.001 THEN BEGIN
        IF KEYWORD_SET(plotYes) THEN BEGIN
          window, 5
;          plot, signal_total[0,*], signal_total[3,*]/MAX(signal_total[3,*])*MAX(signal_total[1,*]), BACKGROUND = 'FFFFFF'x, COLOR = 0, XTITLE='time [us]', YTITLE='Signal [A.U.]', title=STRING(i) + ' encoding step'
;          oplot, signal_total[0,*], signal_total[1,*], LINESTYLE=3, thick = 2, color=35000
;          oplot, signal_total[0,*], signal_total[2,*], thick = 2, color=3500;/MAX(filteredsignal)*MAX(signal)
          H_array = signal_total[3,*]
          signal_array2 = signal_total[1,*]
          signal_array1 = signal_total[2,*]
          cgplot, signal_total[0,*], H_array/MAX(H_array)*MAX(signal_array2), LINESTYLE = 2, thick = 2, title='signal', xtitle='time [us]', ytitle='Signal \ A.U.', BACKGROUND = 'FFFFFF'x, COLOR = 0
          cgplot, signal_total[0,*], signal_array2, thick = 2, LINESTYLE = 0, color='blue',/overplot
          cgplot, signal_total[0,*], signal_array1, thick = 2, color='red',/overplot
        ENDIF
      ENDIF
      ;wait, 0.5
      ;calculate the AUC
      IF ABS(H_value_shape_freq[0]) LT 0.01 THEN H_value_shape_freq[0] = 5.
      signal_1d[i] = (signal_vs_H(H_value_shape_freq,tracer_info,0, signal_total))[1, 1];get the AUC
      signal_total = signal_total*0.
    ENDFOR

    ;reconstruct image based on SVD
    image_1d_y = V ## sv_inv ## TRANSPOSE(U) ## signal_1d
    image_1d_rec = REFORM(image_1d_y)
    correl_image = CORRELATE(image_1d_rec, image_1d)
    
    IF KEYWORD_SET(plotYes) THEN BEGIN
      window, 6
      plot, image_1d_rec, thick=1, title='SVD reconstruction (size='+STRING(Matrix_size)+')', color=0, BACKGROUND = 'FFFFFF'x
      oplot, image_1d, color=3500000, THICK=2;, PSYM=2

    ENDIF
    window, 6
    plot, image_1d_rec, thick=1, title='SVD reconstruction (size='+STRING(Matrix_size)+')', color=0, BACKGROUND = 'FFFFFF'x
    oplot, image_1d, color=3500000, THICK=2;, PSYM=2
    print, 'correlation of 1d reconstructed image = ', correl_image
    image_1d_tvrd = TVRD(TRUE=1)
          
    ;calculate rms_error
    rms_image_error = Sqrt(Total((image_1d_rec - image_1d)^2)/N_Elements(image_1d))
    ;do MTF calculation
    ;IF N_ELEMENTS(freq_array) GE 1 THEN MTF_cal, image_1d, rec_image, freq_array, FOV
  ENDIF

  IF keyword_set(saveDat) THEN BEGIN
    cd, current=old_dir
    ;write the system matrix into a dat file
      filename_System = old_dir+'\LXF_signalKernal_2d_array.dat'
      ;write file
      openw,unit,filename_System,/get_lun
      FOR i=0, Matrix_size-1 DO BEGIN
        FOR j=0, Matrix_size-1 DO BEGIN
          ;; write the data points
          WRITEU, unit, signalKernal_2d_array[i, j]
        ENDFOR
      ENDFOR
      ;; close file and END execution
      FREE_LUN, unit
    
      ;write the signal into a dat file
      filename_1D = old_dir+'\LXF_signal_1d.dat'
      ;write file
      openw,unit,filename_1D,/get_lun
      FOR i=0, Matrix_size-1 DO BEGIN
        ;; write the data points
        WRITEU, unit, signal_1d[i]
      ENDFOR
      ;; close file and END execution
      FREE_LUN, unit
    
      ;write the signal into a dat file
      filename_1D = old_dir+'\LXF_image_1d.dat'
      ;write file
      openw,unit,filename_1D,/get_lun
      FOR i=0, Matrix_size-1 DO BEGIN
        ;; write the data points
        WRITEU, unit, image_1d[i]
      ENDFOR
      ;; close file and END execution
      FREE_LUN, unit
  ENDIF
  
  ;output the parameters
  ;the first 4 are bar_image, the last 2 are image
  ;[FWHM, max_correl, svd_index, rms_error, correl_image, rms_image_error]
  IF keyword_set(bar_phantom_only) THEN BEGIN
    RETURN, [FWHM, max_correl, svd_index, rms_error, 0, 0, svd_slope]
  ENDIF ELSE BEGIN
    IF keyword_set(image_only) THEN BEGIN
      RETURN, [0, 0, 0, 0, correl_image, rms_image_error, svd_slope]
    ENDIF ELSE BEGIN
      IF N_ELEMENTS(FWHM) GT 0 THEN BEGIN
        RETURN, [FWHM, max_correl, svd_index, rms_error, correl_image, rms_image_error, svd_slope]
      ENDIF ELSE BEGIN
        RETURN, [0, 0, 0, 0, correl_image, rms_image_error, svd_slope]
      ENDELSE
    ENDELSE
  ENDELSE
end

;GJ, 2021/11/10, do 2d image simulation based on 1d simulation
;GJ, 
;example: FWHM_info = gradientMPI_2d_simulation();calculate the Shepp-Logan phantom image
FUNCTION gradientMPI_2d_simulation, H_value_shape_freq, tracer_info, matrix_size, N_angles, image_2d_rec, image_2d, saveDat=saveDat, plotYes=plotYes, radonTrans=radonTrans, rasterScan=rasterScan

  ;GJ, 2021/10/29
  ;pre-determine matrix size and number of projection angles
  IF N_ELEMENTS(matrix_size) EQ 0 THEN Matrix_size = 256;192;256;128;64;512;256;128;64;256;64;128;256;128;64;512;256;512;256;128;64;256;512;256;256;64;64;128;180
  IF N_ELEMENTS(N_angles) EQ 0 THEN N_angles = 180

  IF N_ELEMENTS(image_2d) LT 10 THEN BEGIN
    ;filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
    ;image_temp2 = read_image(filename2)
    ;image_2d = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), Matrix_size, Matrix_size));, MAX=130)
    image_2d = BYTSCL(shepp_logan_phantom_2d(Matrix_size))
  ENDIF ELSE BEGIN
    Matrix_size = N_ELEMENTS(image_2d[0,*])
  ENDELSE
  
  ;if radonTrans is set, use radon transform, else use line by line raster scan
  IF KEYWORD_SET(radonTrans) THEN BEGIN
    IF KEYWORD_SET(rasterScan) THEN BEGIN
      ;raster scan
      N_angles = matrix_size
      new_Radon_result = image_2d
    ENDIF ELSE BEGIN
      ;;Calculate and display the Radon transform:
      Radon_result_temp = RADON(image_2d, RHO=rho, THETA=theta);, NTHETA=256, NRHO=256)
      IF KEYWORD_SET(plotYes) THEN iimage, Radon_result_temp, title='Radon Transform'

      ;congrid
      new_Radon_result = CONGRID(Radon_result_temp, N_angles, Matrix_size)
      new_rho = CONGRID(rho, Matrix_size)
      new_theta = CONGRID(theta, N_angles)
      ;  backproject0 = RADON(new_Radon_result, /BACKPROJECT, RHO=new_rho, THETA=new_theta)
      window, 8, XSIZE=N_angles, YSIZE=Matrix_size, title='Radon Transform'
      tvscl, new_Radon_result
    ENDELSE
  ENDIF
  
  rec_result = new_Radon_result*0.
  rec_result_filter = new_Radon_result*0.
  FOR i=0, N_angles-1 DO BEGIN
    image_1d = REFORM(new_Radon_result[i, *])
    FWHM_correl_info = gradientMPI_1d_simulation(H_value_shape_freq, tracer_info, matrix_size, image_1d_tvrd, image_1d, rec_image_1d, sv_inv, V, U);, /plotYes)
    rec_temp = CONGRID(REFORM(rec_image_1d), Matrix_size)
    rec_result[i, *] = rec_temp
    ;do filter
    rec_result_filter[i, *] = Filter_FFT(rec_temp)
    print, 'i = ', i, ', total = ', N_angles 
  ENDFOR
  
  ;if radonTrans is set, use radon transform, else use line by line raster scan
  IF KEYWORD_SET(radonTrans) THEN BEGIN
    IF KEYWORD_SET(rasterScan) THEN BEGIN
      image_2d_rec = rec_result
    ENDIF ELSE BEGIN
      image_2d_rec = RADON(rec_result_filter, /BACKPROJECT, RHO=new_rho, THETA=new_theta)
    ENDELSE
  ENDIF
   
  IF KEYWORD_SET(plotYes) THEN iimage, image_2d_rec, title='MPI using MRI gradient encoding'
  window, 9, XSIZE=2*Matrix_size, YSIZE=Matrix_size, title='MPI using MRI gradient encoding'
  tvscl, image_2d, 0
  tvscl, image_2d_rec, 1
  
  ;use correlate_images to calculate correlation
  correl_image = correl_images(image_2d, image_2d_rec, XSHIFT=0, YSHIFT=0)
  print, 'correlation of 2d image = ', correl_image
  ;calculate rms_error
  rms_image_error = Sqrt(Total((image_2d_rec - image_2d)^2)/N_Elements(image_2d))
  
  RETURN, [0,0,0,0, correl_image, rms_image_error]
END


;GJ, 2021/10/25, generate the phantom with different lps/mm
;freq_arr[10,3] include the 1,2,3,... pixels, the starting index and ending index
FUNCTION image_phantom_res, size, freq_array, plotYes=plotYes

  ;size = 512
  phantom_image = DBLARR(size)*0.
  freq_array_temp=DBLARR(100,3)*0.
  
  n_lps = 4
  n_gap = 4
  n_freq = 0
  start_ind = 4+n_gap
  FOR i = 1, 100 DO BEGIN
    length = i*2*n_lps; five lps
    temp_part = DBLARR(length)*0.
    FOR j = 0, length/i-1, 2 DO temp_part[j*i:(j+1)*i-1] = 255
    IF KEYWORD_SET(plotYes) EQ 1 THEN plot, temp_part
    ending_index = start_ind+length-1
    IF start_ind LT size-1 THEN BEGIN
      IF ending_index LT size-1 THEN BEGIN
        phantom_image[start_ind:ending_index] = temp_part
        freq_array_temp[i-1, *] = [i, start_ind, ending_index]
        n_freq++
      ENDIF
;      ENDIF ELSE BEGIN
;        length_temp = size - start_ind
;        phantom_image[start_ind:*] = temp_part[0:length_temp-1]
;        freq_array[i-1, *] = [i, start_ind, size-1]
;      ENDELSE
      start_ind = ending_index+n_gap
      IF KEYWORD_SET(plotYes) EQ 1 THEN print, 'start index = ', start_ind
    ENDIF
  ENDFOR
  
  IF KEYWORD_SET(plotYes) EQ 1 THEN plot, phantom_image
  
  phantom_image_2d = DBLARR(size, size)*0.
  FOR i=0, size-1 DO BEGIN
    phantom_image_2d[*, i] = phantom_image
  ENDFOR
  IF KEYWORD_SET(plotYes) EQ 1 THEN iimage, phantom_image_2d
  
  freq_array = freq_array_temp[0:n_freq-1, *]
  IF KEYWORD_SET(plotYes) EQ 1 THEN print, freq_array
  return, phantom_image
END

;GJ, 2021/10/27, modulation transfer function
FUNCTION MTF_cal, phantom_image_1d, rec_image_1d, freq_array, FOV, plotYes=plotYes
  
  IF N_ELEMENTS(phantom_image_1d) LT 10 THEN BEGIN
    ;get the MTF phantom image
    Matrix_size = 256;512;512;1024;256;4102;2056;
    phantom_image_1d = image_phantom_res(Matrix_size, freq_array=freq_array, plotYes=0)
  ENDIF ELSE BEGIN
    Matrix_size = N_ELEMENTS(phantom_image_1d)
  ENDELSE
  
  phantom_image_2d = DBLARR(Matrix_size, Matrix_size)*0.
  FOR i=0, Matrix_size-1 DO BEGIN
    phantom_image_2d[*, i] = phantom_image_1d
  ENDFOR

  IF N_ELEMENTS(rec_image_1d) LT 10 THEN BEGIN
    ;generate the noise image
    ;rec_image_2d = SMOOTH(phantom_image_2d, 10)
    rec_image_2d = NOISE_PICK(phantom_image_2d, 0.2, ITER=5)
    rec_image_1d = DBLARR(Matrix_size)*0.
    FOR i=0, Matrix_size-1 DO BEGIN
      rec_image_1d[i] = MEAN(rec_image_2d[i,*])
    ENDFOR
  ENDIF ELSE BEGIN
    rec_image_2d = DBLARR(Matrix_size, Matrix_size)*0.
    FOR i=0, Matrix_size-1 DO BEGIN
      rec_image_2d[*, i] = rec_image_1d
    ENDFOR
  ENDELSE

  ; Display the images side by side:
  IF KEYWORD_SET(plotYes) THEN BEGIN
    IIMAGE, phantom_image_2d, VIEW_GRID=[2,1], VIEW_TITLE='Original', $
      DIMENSIONS=[Matrix_size*2, Matrix_size], WINDOW_TITLE='NOISE_SCATTER Example', $
      /NO_SAVEPROMPT
    IIMAGE, rec_image_2d, /VIEW_NEXT, VIEW_TITLE='Noisy'
    window, 1
    plot, rec_image_1d, color='000000'x, background='FFFFFF'x
    oplot, phantom_image_1d, color='FF0000'x, LINESTYLE=2
  ENDIF

   
  pixel = FOV*1000./Matrix_size ;FOV from m to mm
  n_freq = N_ELEMENTS(freq_array[*,0])
  freq = DBLARR(n_freq)
  modulation = DBLARR(n_freq)
  FOR i=0, n_freq-1 DO BEGIN
    freq[i] = 1./(freq_array[i,0]*2.*pixel)
    temp_phantom = phantom_image_1d[freq_array[i,1]:freq_array[i,2]]
    bright_ind = WHERE(temp_phantom GT 128, count, complement=black_ind, ncomplement=c_count)
    diff_phantom = MEAN(temp_phantom[bright_ind])-MEAN(temp_phantom[black_ind])
    temp_rec = rec_image_1d[freq_array[i,1]:freq_array[i,2]]
    temp_rec[WHERE(temp_rec LT 0, /NULL)] = 0
    ;bar_y[WHERE(bar_y GT 255, /NULL)] = 255
;    IF i EQ 2 THEN BEGIN
;      window, 2
;      plot, temp_phantom, color='000000'x, background='FFFFFF'x, title='freq='+STRING(freq[i])
;      oplot, temp_rec, color='FF0000'x
;    ENDIF
    diff_rec = MEAN(temp_rec[bright_ind])-MEAN(temp_rec[black_ind])
    total_rec = MEAN(temp_rec[bright_ind])+MEAN(temp_rec[black_ind])
    modulation[i] = 1.0*diff_rec/total_rec*100.
  ENDFOR
  ;XRANGE=[0.,10.], 
  IF KEYWORD_SET(plotYes) THEN BEGIN
    iplot, freq, modulation, title='MTF', YRANGE=[0.,100.], XTITLE='Frequency [lp/mm]', YTITLE='Modulation [%]'
    ;show the modulation results and print
    window, 3
    plot, freq, modulation, title='MTF', YRANGE=[0.,100.], XTITLE='Frequency [lp/mm]', YTITLE='Modulation [%]', color='000000'x, background='FFFFFF'x
    print, 'modulation = ', modulation
  ENDIF
  
  ;calculate PSF
  N=256
  freq_grid = FINDGEN(N)/(N-1.)*MAX(freq)
  T_lpmm = freq_grid[1] - freq_grid[0]
  mod_grid = INTERPOL(modulation, freq, freq_grid)
  ;plot, freq_grid, mod_grid
  psf = ABS(FFT(mod_grid))
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  if (is_N_even) then dist_grid = [0.0, X, N/2, -N/2 + X]/(N*T_lpmm) else dist_grid = [0.0, X, -(N/2 + 1) + X]/(N*T_lpmm)
  T_mm = dist_grid[1] - dist_grid[0]
  dist_grid_int = FINDGEN(N)/(N-1.)*dist_grid[10]
  psf_grid_int = INTERPOL(psf[0:10], dist_grid[0:10], dist_grid_int)
  IF KEYWORD_SET(plotYes) THEN BEGIN
    window, 4
    plot, dist_grid_int, psf_grid_int, title='PSF', XTITLE='distance [mm]', color='000000'x, background='FFFFFF'x
  ENDIF
  min_psf_FWHM = MIN(ABS(psf_grid_int-MAX(psf)/2.), psf_FWHM_ind)
  FWHM = 2. * dist_grid_int[psf_FWHM_ind]
  RETURN, FWHM
;  RETURN, [TRANSPOSE(freq), TRANSPOSE(modulation)]
END

;@GJ, 2023/1/18, calculate RFD pSNR
;@GJ, 2023/1/27, using images from SUN MENG
PRO cal_RFD_pSNR_SSIM_MRA
  fn_ori='C:\D_drive\MPI_Tianjie\RotationalDriftMPI\MRA_ori.png'
  fn_rec='C:\D_drive\MPI_Tianjie\RotationalDriftMPI\MRA_rec.png'

;  @GJ, 2023/1/27, images from SUN MENG
;  fn_rec='C:\D_drive\MPI_Tianjie\RotationalDriftMPI\MRA_sunmeng01.png'
  
  im_ori=read_png(fn_ori)
  iimage, im_ori, title='original image'
  im_rec=read_png(fn_rec)
  IF SIZE(im_rec, /N_DIMENSIONS) GT 2 THEN BEGIN
    im_ori_size = SIZE(im_ori, /DIMENSIONS)
    im_rec = CONGRID(REFORM(im_rec[0,*,*]),im_ori_size[0],im_ori_size[1])
  ENDIF
  iimage, im_rec, title='reconstructed image'
  
  MSE=calMSE(im_ori, im_rec);72.1795
  print, 'MSE: ', MSE
  pSNR=calPSNR(im_ori, MSE);68.0337
  print, 'pSNR: ', pSNR, ' dB'
  
  mssim = SSIM(im_ori, im_rec)
  print, 'SSIM: ', mssim
END

;function calSSIM, image1, image2
;  image1 = DOUBLE(BYTSCL(image1))
;  image2 = DOUBLE(BYTSCL(image2))
;;  image1 = DOUBLE(image1)
;;  image2 = DOUBLE(image2)
;  ; 计算结构相似指数 (SSIM)
;  ; 首先定义一些参数
;  K1 = 0.01
;  K2 = 0.03
;  L = 255.0
;  C1 = (K1 * L)^2
;  C2 = (K2 * L)^2
;
;  ; 计算均值
;  mu1 = CONVOL(image1, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE)
;  mu2 = CONVOL(image2, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE)
;
;  ; 计算方差
;  sigma1_sq = CONVOL(image1^2, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE) - mu1^2
;  sigma2_sq = CONVOL(image2^2, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE) - mu2^2
;  sigma12 = CONVOL(image1 * image2, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE) - mu1 * mu2
;
;  ; 计算 SSIM
;  num = (2 * mu1 * mu2 + C1) * (2 * sigma12 + C2)
;  den = (mu1^2 + mu2^2 + C1) * (sigma1_sq + sigma2_sq + C2)
;  ssim = num / den
;  ssim_mean = MEAN(ssim)
;
;  return, ssim_mean
;end

;用于计算原始影像与经过处理后影像的平均绝对误差 MAE
function calMAE, image_or, image_deal
  nDim = size(image_or,  /n_dimension)
  nDim2 = size(image_deal, /n_dimension)
  if nDim le 2 and nDim2 le 2 then begin
    valid_data_or = image_or[where(~finite(image_or, /nan))]
    valid_data_deal = image_deal[where(~finite(image_deal, /nan))]
    if n_elements(valid_data_or) eq n_elements(valid_data_deal) then begin
      nums = n_elements(valid_data_or)
      MAE = total(abs(valid_data_or - valid_data_deal)) / nums
    endif else begin
      print, 'vaild data in 2 image is not equal'
      return, 0
    endelse
  endif else begin
    print, 'image is not le 2 dims.'
    return, 0
  endelse
  return, MAE
end

;;用于计算原始影像与经过处理后影像的SSIM
;function calSSIM, image_or, image_deal
;  nDim = size(image_or,  /n_dimension)
;  nDim2 = size(image_deal, /n_dimension)
;  if nDim le 2 and nDim2 le 2 then begin
;    valid_data_or = image_or[where(~finite(image_or, /nan))]
;    valid_data_deal = image_deal[where(~finite(image_deal, /nan))]
;    if n_elements(valid_data_or) eq n_elements(valid_data_deal) then begin
;      nums = n_elements(valid_data_or)
;      MAE = total(abs(valid_data_or - valid_data_deal)) / nums
;      ; 计算结构相似指数 (SSIM)
;      ; 首先定义一些参数
;      K1 = 0.01
;      K2 = 0.03
;      L = 255.0
;      C1 = (K1 * L)^2
;      C2 = (K2 * L)^2
;
;      ; 计算均值
;      mu1 = CONVOL(image_or, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE)
;      mu2 = CONVOL(image_deal, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE)
;
;      ; 计算方差
;      sigma1_sq = CONVOL(image_or^2, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE) - mu1^2
;      sigma2_sq = CONVOL(image_deal^2, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE) - mu2^2
;      sigma12 = CONVOL(image_or * image_deal, REPLICATE(1.0, 3, 3) / 9.0, /EDGE_TRUNCATE) - mu1 * mu2
;
;      ; 计算 SSIM
;      num = (2 * mu1 * mu2 + C1) * (2 * sigma12 + C2)
;      den = (mu1^2 + mu2^2 + C1) * (sigma1_sq + sigma2_sq + C2)
;      ssim = num / den
;      ssim_mean = MEAN(ssim)
;    endif else begin
;      print, 'vaild data in 2 image is not equal'
;      return, 0
;    endelse
;  endif else begin
;    print, 'image is not le 2 dims.'
;    return, 0
;  endelse
;  return, ssim_mean
;end

;用于计算原始影像与经过处理后影像的均方误差MSE
function calMSE, image_or, image_deal
  image_or = DOUBLE(BYTSCL(image_or))
  image_deal = DOUBLE(BYTSCL(image_deal))
  nDim = size(image_or,  /n_dimension)
  nDim2 = size(image_deal, /n_dimension)
  if nDim le 2 and nDim2 le 2 then begin
    valid_data_or = image_or[where(~finite(image_or, /nan))]
    valid_data_deal = image_deal[where(~finite(image_deal, /nan))]
    if n_elements(valid_data_or) eq n_elements(valid_data_deal) then begin
      nums = n_elements(valid_data_or)
      MSE = total((valid_data_or - valid_data_deal)*(valid_data_or - valid_data_deal)) / nums
    endif else begin
      print, 'vaild data in 2 image is not equal'
      return, 0
    endelse
  endif else begin
    print, 'image is not le 2 dims.'
    return, 0
  endelse
  return, MSE
end

;用于计算单波段影像的峰值信噪比PSNR
function calPSNR, image, MSE
  image = DOUBLE(BYTSCL(image))
  nDim = size(image, /n_dimension)
  if nDim le 2 then begin
    image_max = max(image, /nan)
    PSNR = 20 * alog10(image_max / sqrt(MSE))
    return, PSNR
  endif else begin
    print, 'image dims is not le 2.'
    return, 0
  endelse
end

;----------------------------------------------------------------------
; This is an implementation of the algorithm for calculating the
; Structural SIMilarity (SSIM) index between two images
;
; This is the IDL implementation by Dr. Christiaan Boersma
; (Christiaan.Boersma@nasa.gov) that has been ported from the Matlab
; version available at http://www.cns.nyu.edu/~lcv/ssim/ssim_index.m.
;
; Please refer to the following paper
;
; Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
; quality assessment: From error visibility to structural similarity,"
; IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612,
; Apr. 2004.
;
;----------------------------------------------------------------------
; Permission to use, copy, or modify this software and its documentation
; for educational and research purposes only and without fee is hereby
; granted, provided that this copyright notice and the original authors'
; names appear on all copies and supporting documentation. This program
; shall not be used, rewritten, or adapted as the basis of a commercial
; software or hardware product without first obtaining permission of the
; authors. The authors make no representations about the suitability of
; this software for any purpose. It is provided "as is" without express
; or implied warranty.
;----------------------------------------------------------------------
;
;Input : (1) img1: the first image being compared
;        (2) img2: the second image being compared
;        (3) K: constants in the SSIM index formula (see the above
;            reference). defualt value: K = [0.01 0.03]
;        (4) window: local window for statistics (see the above
;            reference). default widnow is Gaussian given by
;            window = fspecial('gaussian', 11, 1.5);
;        (5) L: dynamic range of the images. default: L = 255
;
;Output: (1) mssim: the mean SSIM index value between 2 images.
;            If one of the images being compared is regarded as
;            perfect quality, then mssim can be considered as the
;            quality measure of the other image.
;            If img1 = img2, then mssim = 1.
;        (2) SSIM_MAP: the SSIM index map of the test image.
;
;Basic Usage:
;   Given 2 test images img1 and img2, whose dynamic range is 0-255
;
;   mssim = SSIM(img1, img2)
;
;Advanced Usage:
;   User defined parameters. For example
;
;   K = [0.05 0.05]
;   window = INTARR(8,8) + 1
;   L = 100
;   mssim = SSIM(img1, img2, K, window, L, SSIM_MAP=ssim_map)
;========================================================================

FUNCTION SSIM,img1,img2,K,window,L,SSIM_MAP=SSIM_MAP

  COMPILE_OPT IDL2

  ON_ERROR, 2

  IF N_PARAMS() LT 2 THEN $
    MESSAGE,"NEED TWO IMAGES"

  img1 = DOUBLE(BYTSCL(img1))
  img2 = DOUBLE(BYTSCL(img2))
  dim1 = SIZE(img1, /DIMENSIONS)

  IF NOT ARRAY_EQUAL(dim1, SIZE(img2, /DIMENSIONS)) THEN $
    MESSAGE,"DIMENSIONS OF BOTH IMAGES DO NOT AGREE"

  IF N_ELEMENTS(dim1) NE 2 THEN $
    MESSAGE,"IMAGE HAS WRONG NUMBER OF DIMENSIONS"

  IF N_PARAMS() LT 3 THEN $
    K = [0.01, 0.03]

  IF K[0] LT 0 OR K[1] LT 0 THEN $
    MESSAGE,"K'S CANNOT BE NEGATIVE"

  K = DOUBLE(K)

  IF N_PARAMS() LT 4 THEN $
    window = GAUSSIAN_FUNCTION([1.5D, 1.5], WIDTH=11)

  dimw = SIZE(window, /DIMENSIONS)

  IF dimw[0] * dimw[1] LT 4 OR $
    dimw[0] GT dim1[0] OR $
    dimw[1] GT dim1[1] OR $
    dimw[0] MOD 2 EQ 0 OR $
    dimw[1] MOD 2 EQ 0 THEN $
    MESSAGE,"WINDOW EITHER TOO SMALL, TOO BIG OR EVEN"

  IF dim1[0] LT dimw[0] OR dim1[1] LT dimw[1] THEN $
    MESSAGE,"IMAGE IS TOO SMALL FOR WINDOW"

  window /= TOTAL(window)

  IF N_PARAMS() LT 5 THEN $
    L = 255

  L = DOUBLE(L)

  C = (K * L)^2

  img1 = DOUBLE(img1)

  img2 = DOUBLE(img2)

  mu1 = CONVOL(img1, window, /EDGE_ZERO)

  mu2 = CONVOL(img2, window, /EDGE_ZERO)

  mu1_sq = mu1^2

  mu2_sq = mu2^2

  mu1_mu2 = mu1 * mu2

  sigma1_sq = CONVOL(img1^2, window, /EDGE_ZERO) - mu1_sq

  sigma2_sq = CONVOL(img2^2, window, /EDGE_ZERO) - mu2_sq

  sigma12 = CONVOL(img1 * img2, window, /EDGE_ZERO) - mu1_mu2

  IF C[0] GT 0 AND C[1] GT 0 THEN $
    ssim_map = ((2D * mu1_mu2 + C[0]) * (2D * sigma12 + C[1])) / $
    ((mu1_sq + mu2_sq + C[0]) * (sigma1_sq + sigma2_sq + C[1])) $
  ELSE BEGIN
    numerator1 = 2D * mu1_mu2 + C[0]
    numerator2 = 2D * sigma12 + C[1]
    denominator1 = mu1_sq + mu2_sq + C[0]
    denominator2 = sigma1_sq + sigma2_sq + C[1]
    ssim_map = MAKE_ARRAY(dim1[0], dim1[1], VALUE=1D)
    index = WHERE(denominator1 * denominator2 GT 0)
    ssim_map[index] = (numerator1[index] * numerator2[index]) / $
      (denominator1[index] * denominator2[index])
    index = WHERE(denominator1 NE 0 AND denominator2 EQ 0)
    ssim_map[index] = numerator1[index] / $
      denominator1[index]
  ENDELSE

  RETURN,MEAN(ssim_map[dimw[0]/2:-dimw[0]/2-1,dimw[1]/2:-dimw[1]/2-1])
END

; 用于计算单波段影像的信息熵H
function calEntroy, image
  nDim = size(image, /n_dimension)
  if nDim le 2 then begin
    pixNums = histogram(image, binsize=1, /nan)
    nums = total(pixNums, /nan)
    pby = pixNums*1d / nums
    imageEntroy = -total(pby*alog(pby) / alog(2), /nan)
    return, imageEntroy
  endif else begin
    print, 'image dims is not le 2.'
    return, 0
  endelse
end

;计算影像辐射质量改善因子IF
function calIF, image_or, image_deal
  mean_or = mean(image_or, dimension=1, /nan)
  row_num = n_elements(mean_or)
  sum_or = 0
  for i=1, row_num-1 do begin
    number = (mean_or[i] - mean_or[i-1]) * (mean_or[i] - mean_or[i-1])
    sum_or = sum_or + number
  endfor

  mean_deal = mean(image_deal, dimension=1, /nan)
  sum_deal = 0
  for i=1, row_num-1 do begin
    number = (mean_deal[i] - mean_deal[i-1]) * (mean_deal[i] - mean_deal[i-1])
    sum_deal = sum_deal + number
  endfor
  IF_cal = 10 * alog(total(sum_or / sum_deal))
  return, IF_cal
end


;@LXF, Liangxiaofeng,2021/11/16,2d shepp-logan phantom function
;input: P=shepp_logan_2d(N)
;        iimage,P
;output:   N*N,2d shepp-logan picture
;examples:
;;P=shepp_logan_phantom_2d(256)
;iimage,P
;P=shepp_logan_phantom_2d(512, /plotYes)
;P=shepp_logan_phantom_2d(256, /plotYes)
;P=shepp_logan_phantom_2d(64, /plotYes)
function shepp_logan_phantom_2d, matrix_size, plotYes=plotYes
  
  IF N_ELEMENTS(matrix_size) GT 0 THEN N=matrix_size ELSE N=256
 
  I=make_array(N,N,/float,value=0)

  e=[[1, 0.69, 0.92, 0, 0, 0],$
    [-0.8, 0.6624, 0.8740, 0, -0.0184, 0],$
    [-0.2, 0.1100, 0.3100, 0.22, 0, -18],$
    [-0.2, 0.1600, 0.4100, -0.22, 0, 18],$
    [0.1, 0.2100, 0.2500, 0, 0.35, 0],$
    [0.1, 0.0460, 0.0460, 0, 0.1, 0],$
    [0.1, 0.0460, 0.0460, 0, -0.1, 0],$
    [0.1, 0.0460, 0.0230, -0.08, -0.605, 0],$
    [0.1, 0.0230, 0.0230, 0, -0.606, 0],$
    [0.1, 0.0230, 0.0460, 0.06, -0.605, 0]]
  grey=transpose(e[0,*])     ;密度
  x_a=transpose(e[1,*])      ;x半轴长
  y_a=transpose(e[2,*])      ;y半轴长
  x_b=transpose(e[3,*])      ;x坐标
  y_b=transpose(e[4,*])      ;y坐标
  a=transpose(e[5,*])        ;旋转角度
  phi = a*!pi/180     ;转换成弧度
  ;旋转矩阵
  T11=cos(phi)
  T12=sin(phi)
  T21=-sin(phi)
  T22=cos(phi)

  for k1=0,N-1 do begin
    for k2=0,N-1 do begin
      x0 = (k2-N/2.0)/(N/2.0)
      y0 = (k1-N/2.0)/(N/2.0)
      XX=T11##x0+T12##y0;
      YY=T21##x0+T22##y0
      ;x = cos(phi)*x0+sin(phi)*y0
      ;y = -sin(phi)*x0+cos(phi)*y0
      x=XX-x_b
      y=YY-y_b
      ellipsoid = dblarr(10)*0.
      ind = dblarr(10)*0.
      ellipsoid=x^2/(x_a^2)+y^2/(y_a^2);
      ind = ellipsoid le 1.0
      grayval=total(grey*ind)
      I[k1,k2]=I[k1,k2]+grayval
      ;print,grayval
    endfor
  endfor
  if KEYWORD_SET(plotYes) THEN iimage,transpose(I)
  return,transpose(I)
end

;@LXF,2021/11/19,3d shepp-logan phantom function
;input: P=shepp_logan_3d(N)
;
;output:   N*N*N,3d shepp-logan
;examples:
;;P=shepp_logan_phantom_3d(256)
;b=P[*,N/2,*]
;iimage,reform(b)        ;The output is a slice in the x direction
;b=P[N/2,*,*]
;iimage,reform(b)        ;The output is a slice in the y direction
;b=P[*,*,N/2]
;iimage,reform(b)        ;The output is a slice in the z direction
;P=shepp_logan_phantom_3d(512)
;P=shepp_logan_phantom_3d(256)
;P=shepp_logan_phantom_3d(64)


function shepp_logan_phantom_3d, matrix_size;, plotYes=plotYes

  IF N_ELEMENTS(matrix_size) GT 0 THEN N=matrix_size ELSE N=256
  ;N=128
  e=[[0,       0,      0,   0.69,   0.92,   0.81, 0,   0,  0,     1],$
    [0,      -0.0184, 0,   0.6624, 0.8740, 0.78, 0,   0,  0,  -0.8],$
    [0.22,    0,      0,   0.11,   0.31,   0.22, -18, 0,  10, -0.2],$
    [-0.22,   0,      0,   0.16,   0.41,   0.28, 18,  0,  10, -0.2],$
    [0,      0.35,  -0.15, 0.21,   0.25,   0.41, 0,   0,   0,  0.1],$
    [0,      0.1,    0.25, 0.046,  0.046,  0.05, 0,   0,   0,  0.1],$
    [0,     -0.1,    0.25, 0.046,  0.046,  0.05, 0,   0,   0,  0.1],$
    [-0.08, -0.605,  0,    0.046,  0.023,  0.05, 0,   0,   0,  0.1],$
    [0,     -0.606,  0,    0.023,  0.023,  0.02, 0,   0,   0,  0.1],$
    [0.06,  -0.605,  0,    0.023,  0.046,  0.02, 0,   0,   0,  0.1]]
  I=make_array(N,N,N,/float,value=0)
  x_a=transpose(e[0,*])         ;x的坐标
  y_a=transpose(e[1,*])         ;y的坐标
  z_a=transpose(e[2,*])         ;z的坐标
  x_b=transpose(e[3,*])               ;x方向的半长轴
  y_b=transpose(e[4,*])
  z_b=transpose(e[5,*])
  phi=transpose(e[6,*])
  gamma=transpose(e[7,*])
  theta=transpose(e[8,*])
  rho=transpose(e[9,*])          ;密度
  phi1=phi*!pi/180
  gamma1=gamma*!pi/180
  theta1=theta*!pi/180

  ;计算旋转矩阵

  T11=cos(theta1)*cos(phi1)-cos(gamma1)*sin(phi1)*sin(theta1)
  T12=-(cos(theta1)*sin(phi1)+cos(gamma1)*cos(phi1)*sin(theta1))
  T13=sin(theta1)*sin(gamma1)
  T21=-sin(theta1)*cos(phi1)-cos(gamma1)*sin(phi1)*cos(theta1)
  T22=-(-sin(theta1)*sin(phi1)+cos(gamma1)*cos(phi1)*cos(theta1))
  T23=cos(theta1)*sin(gamma1)
  T31=sin(gamma1)*sin(phi1)
  T32=-sin(gamma1)*cos(phi1)
  T33=cos(gamma1)

  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start

  for k1=0,N-1 do begin
    ;update the progress bar
    count=(k1*1.0/N)*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    
    for k2=0,N-1 do begin
      for k3=0,N-1 do begin         ;k1,k2,k3等价于x，y，z
        x0=(k2-N/2.0)/(N/2.0)       ;归一化坐标，将其转化为x，y，z坐标
        y0=(k1-N/2.0)/(N/2.0)
        z0=(k3-N/2.0)/(N/2.0)
        XX=T11*x0+T12*y0+T13*z0
        YY=T21*x0+T22*y0+T23*z0
        ZZ=T31*x0+T32*y0+T33*z0
        x=XX-x_a
        y=YY-y_a
        z=ZZ-z_a
        ellipsoid=x^2/x_b^2+y^2/y_b^2+z^2/z_b^2
        ind=ellipsoid le 1
        grayval=total(rho*ind)
        I[k1,k2,k3]=I[k1,k2,k3]+grayval
        ;print,I
      endfor
    endfor
  endfor
 
 
  ;destroy the progress bar
  progressbar -> Destroy
  
  return,I




end


PRO simulation, base, current, length

  length = length * 0.01

  field_amp = DBLARR(201, 201) * 0.0 ; 磁场强度

  x_field = DBLARR(201, 201) * 0.0 ; x方向场强

  y_field = DBLARR(201, 201) * 0.0 ; y方向场强

  FOR i = 0, 200 DO BEGIN

    FOR j = 0, 200 DO BEGIN

      dist = SQRT((i - 100.0) ^ 2 + (j - 100.0) ^ 2) * length; 距离

      IF dist EQ 0 THEN field_amp[i, j] = 0 ELSE field_amp[i, j] = current / dist;

      x_field[i, j] = field_amp[i, j] * (j - 100.0) / length

      y_field[i, j] = field_amp[i, j] * (i - 100.0) / length

    ENDFOR

  ENDFOR

  ;iimage, field_amp, title = 'Field Amplitude of single line'

  ;iimage, x_field, title = 'x-direction Field of single line'

  ;iimage, y_field, title = 'y-direction Field of single line'
  sum_x_field = DBLARR(100, 100) * 0.0

  sum_y_field = DBLARR(100, 100) * 0.0

  sum_field_amp = DBLARR(100, 100) * 0.0



  sum_x_field = x_field[101:200, 0:99] - x_field[0:99, 0:99] + x_field[0:99, 101:200] - x_field[101:200, 101:200]

  sum_y_field = y_field[101:200, 0:99] - y_field[0:99, 0:99] + y_field[0:99, 101:200] - y_field[101:200, 101:200]

  sum_field_amp = SQRT(sum_x_field ^ 2 + sum_y_field ^ 2)



  draw1 = widget_info(base, find_by_uname='draw1') ; 通过widget_info获取组件返回组件id

  widget_control, draw1, get_value=window1 ; 获取组件draw1的值并赋值给window1

  wset, window1

  surface, field_amp[50:150, 50:150], zrange=[0, 200], title='Field Amplitude of single line', background='ffffff'x, color='0000ff'x



  draw2 = widget_info(base, find_by_uname='draw2') ; 通过widget_info获取组件返回组件id

  widget_control, draw2, get_value=window2 ; 获取组件draw2的值并赋值给window3

  wset, window2

  surface, sum_field_amp, zrange=[0, 500], title='Field Amplitude of Anderson coils', background='ffffff'x, color='00ff00'x



  ;contour, sum_field_amp, levels = FINDGEN(100) / 99.0 * MAX(sum_field_amp), title = 'Field Amplitude'

  ;iimage, sum_x_field, title = 'x-direction Field Amplitude of Anderson coils'

  ;iimage, sum_y_field, title = 'y-direction Field Amplitude of Anderson coils'



END
PRO slider_event_current, event

  widget_control, event.id, get_value=current ; 获取组件slider的值并赋值



  slider_length = widget_info(event.top, find_by_uname='slider_length')

  widget_control, slider_length, get_value=length



  print, length

  print, current

  simulation, event.top, current, length ; 执行上述函数

END
PRO slider_event_length, event

  widget_control, event.id, get_value=length ; 获取组件slider的值并赋值
  slider_current = widget_info(event.top, find_by_uname='slider_current')

  widget_control, slider_current, get_value=current



  print, length

  print, current

  simulation, event.top, current, length ; 执行上述函数

END
PRO homework5
  current = 50.0 ; 电流

  length = 50.0; 中心点距离正方形区域顶点的长度



  base = widget_base(/column) ; 主组件, 其子组件按列分布(If this keyword is included, the base lays out its children in columns.)

  menubase = widget_base(base, /column) ; 菜单组件, 其子组件按列分布

  base_parameter = widget_base(menubase, /row) ; 参数组件, 其子组件按行分布



  label_current = widget_label(base_parameter, value='Current : ')

  slider_current = widget_slider(base_parameter, xsize=512, event_pro='slider_event_current', uname='slider_current')

  label_length = widget_label(base_parameter, value='length : ')

  slider_length = widget_slider(base_parameter, xsize=512, event_pro='slider_event_length', uname='slider_length')

  ;slider_length = widget_slider(base_parameter, minimum=0.1, maximum=1, xsize=512, event_pro='slider_event_length', uname='slider_length')



  drawbase = widget_base(base, /column) ; 图像主组件, 其子组件按列分布

  drawbase1 = widget_base(drawbase, /row) ; 图像组件1, 其子组件按行分布

  draw1 = widget_draw(drawbase1, xsize=512, ysize=400, retain=2, uname='draw1')

  draw2 = widget_draw(drawbase1, xsize=512, ysize=400, retain=2, uname='draw2')





  widget_control, base, /realize

  widget_control, base, set_uvalue=current

  widget_control, base, set_uvalue=length

  widget_control, slider_current, set_value=current

  widget_control, slider_length, set_value=length

  simulation, base, current, length

  xmanager, 'homework5', base



END

;GJ, 2022/01/02
;compare the signal of sonymag and perimag
;GJ, 2022/1/3, return the stucture and compare
pro signal_compare
 filename1 = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20211228\5mT-10%-Perimag-方波激励\peri-gly10.txt'
 filename2 = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20211228\5mT-10%-Synomag-方波激励\syno_gly10.txt'
 fit_1 = read_pulsed_signal_dat(filename1); s1 [500,3]
 fit_2 = read_pulsed_signal_dat(filename2); s2 [500,3]
 H_array = fit_1.H; REFORM(s2[*,0])
 signal_array1 = fit_1.signal; REFORM(s1[*,1])
 signal_array2 = fit_2.signal; REFORM(s2[*,1])
 
 window, 6
 cgplot, -H_array/MAX(H_array)*MAX(signal_array2), LINESTYLE = 2, thick = 2, title='signal', xtitle='time [us]', ytitle='Signal \ A.U.', BACKGROUND = 'FFFFFF'x, COLOR = 0
 cgplot, signal_array2, thick = 2, LINESTYLE = 0, color='blue',/overplot
 cgplot, signal_array1, thick = 2, color='red',/overplot
 
END

PRO color_MPI_simulation

;  filename1 = 'C:\D_drive\MPI_Tianjie\Xinfeng\SPimage\p.png'
;  image_temp1 = read_image(filename1)
;  
;  p_image = DBLARR(4,4)
;  FOR i=0,3 DO BEGIN
;    FOR j=0,3 DO BEGIN
;      p_image[i,j] = TOTAL(image_temp1[*,i,j])
;    ENDFOR
;  ENDFOR
;  
;  filename2 = 'C:\D_drive\MPI_Tianjie\Xinfeng\SPimage\s.png'
;  image_temp2 = read_image(filename2)
;
;  s_image = DBLARR(4,4)
;  FOR i=0,3 DO BEGIN
;    FOR j=0,3 DO BEGIN
;      s_image[i,j] = TOTAL(image_temp2[*,i,j])
;    ENDFOR
;  ENDFOR
   
   p_image = DBLARR(4,4)*0.
   p_image[3, 0:2] = 1.
   p_image[2, 0] = 1. & p_image[2, 2:3] = 1.
   p_image[1, 0:2] = 1.
   p_image[0,0] = 1.
   p_image = TRANSPOSE(p_image)
;   iimage, p_image
   par_name = 'peri_gly1_'
   rec_p_image = DBLARR(4,4)*0.
   signal_p_ptr = PTRARR(4, /ALLOCATE_HEAP)
   FOR i=0, 3 DO BEGIN
     image_1d = REFORM(p_image[i, *])
     p_rec_result = system_matrix_calibration(par_name, image_1d)
     rec_p_image[i, *] = p_rec_result.y_rec
     *(signal_p_ptr[i]) = p_rec_result.signal_vs_time
   ENDFOR
   iimage, rec_p_image, title='reconstructed peri image'
   
   s_image = DBLARR(4,4)*0.
   s_image[3, 1:3] = 1.
   s_image[2, 1] = 1.
   s_image[1, 2] = 1.
   s_image[0,0:2] = 1.
   s_image = TRANSPOSE(s_image)
;   iimage, TRANSPOSE(s_image)
   par_name = 'syno_gly1_'
   rec_s_image = DBLARR(4,4)*0.
   ;define the signal array
   signal_s_ptr = PTRARR(4, /ALLOCATE_HEAP)
   FOR i=0, 3 DO BEGIN
    image_1d = REFORM(s_image[i, *])
    s_rec_result = system_matrix_calibration(par_name, image_1d)
    rec_s_image[i, *] = s_rec_result.y_rec
    *(signal_s_ptr[i]) = s_rec_result.signal_vs_time
   ENDFOR
   
   iimage, rec_s_image, title='reconstructed syno image'
   
;   window, 2
;   plot, signal_vs_time[0, *], YRANGE=[MIN(signal_vs_time), MAX(signal_vs_time)]
;   oplot, signal_vs_time[1, *]
;   oplot, signal_vs_time[2, *]
;   oplot, signal_vs_time[3, *]

   ;mix the p and s original image
   ps_image = s_image + p_image
   iimage, ps_image, title='original P+S image'
   signal_ps_ptr = PTRARR(4, /ALLOCATE_HEAP)
   FOR i=0, 3 DO BEGIN
    *(signal_ps_ptr[i]) = *(signal_p_ptr[i]) + *(signal_s_ptr[i])
    temp_signal_arr = *(signal_p_ptr[i]) + *(signal_s_ptr[i])
    FOR j=0, 3 DO BEGIN
      temp_signal = REFORM(temp_signal_arr[j, *])
      AUC_ps = TOTAL(temp_signal)
      AUC_p = TOTAL((*(signal_p_ptr[i]))[j, *])
      AUC_s = TOTAL((*(signal_s_ptr[i]))[j, *])
      
      Y = TOTAL(temp_signal, /CUMULATIVE)
      X = FINDGEN(N_ELEMENTS(Y))
      ;tri-exponential decay
      C = [-5.0, -1./8.5, -5.0, -1./40., -5.0, -1./100., MAX(Y)]
      fixc = [1, 0, 1, 0, 1, 0, 1]
      ;Compute the parameters.
      yfit_triex = CURVEFIT(X, Y, weights, C, SIGMA, FITA=fixc, FUNCTION_NAME='gfunct_triex')
      ;Print the parameters returned in A.
      PRINT, 'Function parameters (Triex): ', C
      tri_tao = 1./ABS([C[1], C[3], c[5]])
      tri_tao_sort = tri_tao[SORT(tri_tao)]
      tri_tao_1 = tri_tao_sort[0]
      tri_tao_2 = tri_tao_sort[1]
      tri_tao_3 = tri_tao_sort[2]
      print, 'Relaxation time_1 = ', tri_tao_1, ' [us]'
      print, 'Relaxation time_2 = ', tri_tao_2, ' [us]'
      print, 'Relaxation time_3 = ', tri_tao_3, ' [us]'
    ENDFOR
   ENDFOR
   
   par_name = 'syno_gly1_'
   rec_ps_image = DBLARR(4,4)*0.
   FOR i=0, 3 DO BEGIN
     ps_rec_result = system_matrix_calibration(par_name, image_1d, *(signal_ps_ptr[i]))
     rec_ps_image[i, *] = ps_rec_result.y_rec
   ENDFOR
   
   iimage, rec_ps_image, title='reconstructed P+S image'
   
   PTR_FREE, signal_s_ptr
   PTR_FREE, signal_p_ptr
   PTR_FREE, signal_ps_ptr
   
END

;GJ, 2022/01/09
;system_matrix_calibration
FUNCTION system_matrix_calibration, par_name, image_1d, synth_signal_vs_time

 ;par_name = 'syno_gly1_'
 ;par_name = 'peri_gly1_'
 dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\'
 dir_name = 'D:\MPI_MPT\Xinfeng\AWR_meas\'
 
 filename1 = dir_name + par_name + '2_5mT.txt'
 fit_2_5 = read_pulsed_signal_dat(filename1)
 
 filename2 = dir_name + par_name + '5mT.txt'
 fit_5 = read_pulsed_signal_dat(filename2)
 
 filename3 = dir_name + par_name + '7_5mT.txt'
 fit_7_5 = read_pulsed_signal_dat(filename3)
 
 N_tp = N_ELEMENTS(fit_7_5.signal_n)
 ;get AUC array and set the calibrtaion to system_matrix
 AUC_array = [fit_7_5.AUC_n, fit_5.AUC_n, fit_2_5.AUC_n, 0, fit_2_5.AUC_p, fit_5.AUC_p, fit_7_5.AUC_p]
 signal_array = TRANSPOSE([[fit_7_5.signal_n], [fit_5.signal_n], [fit_2_5.signal_n], [fit_2_5.signal_n*0.], [fit_2_5.signal_p], [fit_5.signal_p], [fit_7_5.signal_p]])
 system_matrix = DBLARR(4, 4)
 ;image_1d = [0.5, 0.6, 1.3, 0.5]
 signal_1d = DBLARR(4)*0.
 FOR i=0, 3 DO system_matrix[i, *] = AUC_array[i:i+3]
 ;print, 'system matrix: ', system_matrix
 
 ;calculate the signal based on the each voxel
 signal_vs_time = DBLARR(4, N_tp)
 FOR i=0, 3 DO BEGIN
  FOR j=0, 3 DO BEGIN
    signal_vs_time[i, *] += image_1d[j]*signal_array[i+j, *]
  ENDFOR
 ENDFOR
 
 ;if there is input of synthized signal, then use this one
 IF N_ELEMENTS(synth_signal_vs_time) GT 10 THEN BEGIN
  signal_vs_time = synth_signal_vs_time
 ENDIF
 
 ; window, 2
; plot, signal_vs_time[0, *], YRANGE=[MIN(signal_vs_time), MAX(signal_vs_time)]
; oplot, signal_vs_time[1, *]
; oplot, signal_vs_time[2, *]
; oplot, signal_vs_time[3, *]
 
 FOR i=0, 3 DO BEGIN
  signal_1d[i] = TOTAL(signal_vs_time[i, *])*fit_7_5.delta_t
 ENDFOR
 signal_1d_simple = system_matrix ## image_1d
 diff_signal = signal_1d_simple - signal_1d
 ;print, 'signal diff: ', diff_signal

 ;check the revert matrix
 sm_invert = INVERT(system_matrix, /DOUBLE)
 y_rec = sm_invert ## signal_1d

 ;calcualte the correlation coeff
 correl = CORRELATE(y_rec, image_1d)
 ;print, 'image diff: ', y_rec - image_1d
 print, 'correlation coeff: ', correl
 
 rec_result = { $
 y_rec: y_rec, $
 signal_vs_time: signal_vs_time $
 }
 
 RETURN, rec_result
END

;GJ, 2022/1/19
;output all fitted results
PRO batch_pulsed_fitting
  
  ;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\syno_gly40_7_5mT.txt'
  ;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\peri_gly40_7_5mT.txt'
  ;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\peri_gly3_7_5mT.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\peri_gly1_7_5mT.txt'
;    signalOnly = 1
;    test_fit_para = read_pulsed_signal_dat(test_fn)
;    test_fit_para = read_pulsed_signal_dat(test_fn, signalOnly=signalOnly)
  
  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\'
  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\'
  dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)
  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  output_filename = dir_name+'fit_para20220517.xls'
  
  ;write file
  openw,unit,output_filename,/get_lun
  first_line =['filename', 'AUC_hold', 'M_diff_hold', 'size', 'tao_Debye_mono', 'tao_mono', 'tao_biex_1', 'tao_biex_2', 'tao_biex_1_comp', 'tao_biex_2_comp', 'tao_biex_12combined']
  printf, unit, FORMAT = '(15(A, %"\t"))', first_line

  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  H_peak_array_p = DBLARR(nrfile)*0.
  para_p_ptr = PTRARR(nrfile, /ALLOCATE_HEAP)
  FOR i=0L, nrfile-1L DO BEGIN
    filename = FILE_BASENAME(a_file[i], '.txt')
    fit_para = read_pulsed_signal_dat(a_file[i])
    *(para_p_ptr[i]) = fit_para
    H_peak_array_p[i] = fit_para.H_peak_p
    printf, unit, FORMAT = '(15(A, %"\t"))', filename, fit_para.AUC_n, fit_para.M_diff_n, $
      fit_para.particle_size_mono, fit_para.tao_debye_mono, fit_para.tao, fit_para.bi_tao, fit_para.bi_tao_comp, fit_para.bi_tao_combine
  ENDFOR

  FREE_LUN, unit

END



;GJ, 2022/1/27, do batch fitting of sine excitation
PRO sine_debye_fitting


  ;; write the data points
  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\500Hz-10mT-Synomag\202112211.txt'
  ;test_fit_para = read_sine_signal_dat(test_fn)
;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\1000Hz-5mT-Synomag\202112219.txt'
  
;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\8000Hz-5mT-Synomag\202112211.txt'
  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\8000Hz-10mT-Synomag\202112215.txt'
  ;test_fit_para = read_sine_signal_dat(test_fn)
  
;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\5000Hz-10mT-Synomag\202112211.txt'
  
  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\2000Hz-10mT-Synomag\202112219.txt'
;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\2000Hz-5mT-Synomag\202112211.txt'
  test_fit_para = read_sine_signal_dat(test_fn)
  H = test_fit_para.H
  Nelem = N_ELEMENTS(H)
  particle_size = 30.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  FWHM_particle = 4.16/beta_particle/3.5; mT/mm
  M = Msat_kAm*(1./TANH(beta_particle*H) - 1./(beta_particle*H))
  Max_H = MAX(H, max_ind_H)
  frequency = test_fit_para.frequency
  N_period = 1000000/(frequency)
  shift_step = N_period/2-max_ind_H
  signal_langevin_temp = -(M[1:*]-M[0:Nelem-1])/test_fit_para.delta_t; kA/m/us
  signal_meas = test_fit_para.signal[0:Nelem-2]
  time = test_fit_para.time[0:Nelem-2]
  signal_langevin = signal_langevin_temp/MEAN(ABS(signal_langevin_temp))*MEAN(ABS(signal_meas))
;  signal_langevin = SHIFT(signal_langevin, -shift_step)
  window, 1
  plot, time, signal_langevin, title='Langevin & Signal, '+STRTRIM(STRING(frequency),1)+' Hz', XTITLE='time [us]', YTITLE='signal [a.u.]' 
  oplot, time, signal_meas, line=2
  oplot, time, ABS(H)
  
  signal_meas[0:Nelem/2-1] = MEAN(signal_meas)
  signal_langevin[0:Nelem/2-1] = MEAN(signal_langevin)
;  time = time[0:N_ELEMENTS(signal_meas)-1]
  
  relaxation_time = FINDGEN(100)+1.
  correlation_array = DBLARR(100)*0.
  MSE_array = DBLARR(100)*0.+10000.
  FOR i = 0, 99 DO BEGIN
   signal_nat = non_adiabatic(time, signal_langevin, relaxation_time[i])
   correlation_array[i] = CORRELATE(signal_nat, signal_meas)
   MSE_array[i] = calMSE(signal_nat, signal_meas)
  ENDFOR
  
  ;max_corr = MAX(correlation_array, max_ind)
  min_MSE = MIN(MSE_array, max_ind)
  print, 'relaxation time = ', relaxation_time[max_ind]
  signal_nat = non_adiabatic(time, signal_langevin, relaxation_time[max_ind])
  ;  oplot, time, signal_nat, line=3
  corre_Mono = CORRELATE(signal_meas, signal_nat)
  PRINT, 'correlation coeff (Mono): ', corre_Mono
  
  window, 2
  cgplot, time, signal_langevin, LINESTYLE = 2, thick = 2, title='Debye Mono-ex Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='signal [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  cgplot, time, signal_meas, thick = 2, LINESTYLE = 0, color='blue',/overplot
  cgplot, time, signal_nat, thick = 2, color='red',/overplot
  
  ;calculate the difference
  max_signal = max(signal_meas, max_ind_s)
  max_H = max(H[Nelem/2:*], max_ind_H)
  time_diff = 5;ABS(max_ind_s - max_ind_H - Nelem/4.)
  print, 'time diff = ', time_diff, ' us'
  N_Neel = ROUND(time_diff*10.)*2
  Neel_array = (-FINDGEN(N_Neel)/10.+N_Neel/20.)*test_fit_para.delta_t
  N_Brownian = ROUND(60)
  Brownian_array = (FINDGEN(N_Brownian)+time_diff)*test_fit_para.delta_t*2.
  Neel_ratio = FINDGEN(101)/100.
  correlation_array = DBLARR(N_Neel, N_Brownian, 101)*0.
  MSE_array = DBLARR(N_Neel, N_Brownian, 101)*0.+10000.
  FOR i = 0, N_Neel-1 DO BEGIN
    signal_Neel = non_adiabatic(time, signal_langevin, Neel_array[i])
    FOR j = 0, N_Brownian-1 DO BEGIN
      signal_Brownian = non_adiabatic(time, signal_langevin, Brownian_array[j])
      FOR k = 0, 100 DO BEGIN
        signal_nat = Neel_ratio[k]*signal_Neel + (1.-Neel_ratio[k]) * signal_Brownian
        correlation_array[i, j, k] = CORRELATE(signal_nat, signal_meas)
        MSE_array[i, j, k] = calMSE(signal_nat, signal_meas)
      ENDFOR
    ENDFOR
  ENDFOR
  
;  max_corr = MAX(correlation_array, max_ind_corr)
  min_MSE = MIN(MSE_array, max_ind_corr)
  max_index = ARRAY_INDICES(correlation_array, max_ind_corr)
  print, 'Neel time = ', Neel_array[max_index[0]]
  print, 'Brownian time = ', Brownian_array[max_index[1]]
  print, 'Neel ratio = ', Neel_ratio[max_index[2]]
  
  signal_Neel = non_adiabatic(time, signal_langevin, Neel_array[max_index[0]])
  signal_Brownian = non_adiabatic(time, signal_langevin, Brownian_array[max_index[1]])
  signal_nat = Neel_ratio[max_index[2]]*signal_Neel + (1.-Neel_ratio[max_index[2]]) * signal_Brownian
  ;  oplot, time, signal_nat, line=3
  corre_Biex = CORRELATE(signal_meas, signal_nat)
  PRINT, 'correlation coeff (Biex): ', corre_Biex
  
  window, 3
  plot, time, signal_meas, title='Biex Relaxation Fitting, '+STRTRIM(STRING(frequency),1)+' Hz'
  oplot, time, Neel_ratio[max_index[2]]*signal_Neel, LINESTYLE = 1
  oplot, time, (1.-Neel_ratio[max_index[2]]) * signal_Brownian, LINESTYLE = 2
  oplot, time, signal_nat, LINESTYLE = 3
  
  window, 4
  cgplot, time, signal_langevin, LINESTYLE = 2, thick = 2, title='Debye Biex Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='signal [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  cgplot, time, signal_meas, thick = 2, LINESTYLE = 0, color='blue',/overplot
  cgplot, time, signal_nat, thick = 2, color='red',/overplot
  
  ;calculate the dual-Debye model
  ;2022/1/30
  N_Neel = 20
  Neel_array = FINDGEN(N_Neel) + 1.
  N_Brownian = 20
  Brownian_array = FINDGEN(N_Brownian)*5. + 20.
  N_Neel_Comp = 11
  Neel_Comp_array = FINDGEN(N_Neel_Comp)/10.
  N_Brownian_Comp = 11
  Brownian_Comp_array = FINDGEN(N_Brownian_Comp)/10.
  N_coeff = 11
  coeff_array = FINDGEN(N_coeff)/5.+0.2
  correlation_array = DBLARR(N_Neel, N_Brownian, N_Neel_Comp, N_Brownian_Comp, N_coeff)*0.
  MSE_array = DBLARR(N_Neel, N_Brownian, N_Neel_Comp, N_Brownian_Comp, N_coeff)*0.+10000.
  FOR i=0, N_Neel-1 DO BEGIN
    signal_Neel = non_adiabatic(time, signal_langevin, Neel_array[i])
    FOR j=0, N_Brownian-1 DO BEGIN
      signal_Brownian = non_adiabatic(time, signal_langevin, Brownian_array[j])
      FOR k=0, N_Neel_Comp-1 DO BEGIN
        a = Neel_Comp_array[k]
        FOR l=0, N_Brownian_Comp-1 DO BEGIN
          b = Brownian_Comp_array[l]
          c = 1. - a - b
          IF c LT 0 THEN CONTINUE; jump to next iteration
          FOR m=0, N_coeff-1 DO BEGIN
            signal_nat = coeff_array[m] * (a*signal_Neel + b * signal_Brownian + c * signal_langevin)
            correlation_array[i, j, k, l, m] = CORRELATE(signal_nat, signal_meas)
            MSE_array[i, j, k, l, m] = calMSE(signal_nat, signal_meas)
          ENDFOR
        ENDFOR
      ENDFOR
    ENDFOR
  ENDFOR
  
;  max_corr = MAX(correlation_array, max_ind_corr)
  min_MSE = MIN(MSE_array, max_ind_corr)
  max_index = ARRAY_INDICES(correlation_array, max_ind_corr)
  a = Neel_Comp_array[max_index[2]]
  b = Brownian_Comp_array[max_index[3]]
  coeff = coeff_array[max_index[4]]
  print, '(Dual Debye) Neel time = ', Neel_array[max_index[0]]
  print, '(Dual Debye) Brownian time = ', Brownian_array[max_index[1]]
  print, '(Dual Debye) Neel component = ', a
  print, '(Dual Debye) Brownian component = ', b
  print, '(Dual Debye) Adiabatic component = ', (1.-a-b)
  print, '(Dual Debye) Neel ratio = ', a/(a + b)
  
  signal_Neel = non_adiabatic(time, signal_langevin, Neel_array[max_index[0]])
  signal_Brownian = non_adiabatic(time, signal_langevin, Brownian_array[max_index[1]])
  signal_nat = coeff * (a*signal_Neel + b * signal_Brownian + (1.-a-b) * signal_langevin)
  ;  oplot, time, signal_nat, line=3
  corre_DualDebye = CORRELATE(signal_meas, signal_nat)
  PRINT, '(Dual Debye) Correlation coeff: ', corre_DualDebye

  window, 5
  cgplot, time, signal_langevin, LINESTYLE = 2, thick = 2, title='Dual Debye Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='signal [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  cgplot, time, signal_meas, thick = 2, LINESTYLE = 0, color='blue',/overplot
  cgplot, time, signal_nat, thick = 2, color='red',/overplot
END

;GJ, 2022/2/16, generate .001 file for Upen
PRO batch_generate_UPEN_001_files

  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_fields\'
  dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)

  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  
  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  FOR i=0L, nrfile-1L DO BEGIN
    filename_base = FILE_BASENAME(a_file[i], '.txt')
    fit_para = read_pulsed_signal_dat(a_file[i])
    print, 'i = ', i, ', filename = ', filename_base
  ENDFOR


END

;GJ, 2022/2/15, calculate the UPEN eps file
;GJ, 2022/7/7, add sensitivity calculation
PRO batch_read_pulsed_sensitivity, sState

  dir_name = 'C:\UpenWin\bin\MPI01\'
  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly1\'
  
  ;do the folder selection
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    dir_name = DIALOG_PICKFILE(PATH=sState.directory, /MUST_EXIST, TITLE="txt files", /DIRECTORY)
  ENDIF ELSE BEGIN
    dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)
  ENDELSE

  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  IF STRLEN(dir_name) EQ 0 THEN RETURN
  ;@GJ, 2022/7/2, update the directory
  IF N_ELEMENTS(sState) GT 0 THEN sState.directory = dir_name

  ;GJ, 2022/2/6, redo the fitting and write the results
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  output_filename_temp = dir_name+current_time+'_sensitivity.xls'

  ;Select a file to save exported data, the default opened directory is the one you select just now
  output_filename=DIALOG_PICKFILE(FILE=output_filename_temp, FILTER = ['*.xls'], title='Save Your Sensitivity Output as EXCEL file', PATH=dir_name)
  ;// open an ASCII file for output
  ;IF STRCMP(STRMID(aSaveFile, 2, 3, /REVERSE_OFFSET), 'txt',/FOLD_CASE)   THEN BEGIN
  IF (STRLEN(output_filename) EQ 0)  THEN BEGIN
    infowarning = DIALOG_MESSAGE('No sensitivity output excel file selected!', /ERROR)
    RETURN
  ENDIF

  ;write file
  openw,unit,output_filename,/get_lun
  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  IF nrfile EQ 0 THEN BEGIN
    infowarning = DIALOG_MESSAGE('No txt file found!', /ERROR)
    RETURN
  ENDIF

  ;@GJ, 2022/7/7, define the title
  IF nrfile EQ 9 THEN BEGIN
    GlyCon = [1.,3.,5.,10.,15.,20.,30.,40.,50.]
    VisArray = [0.912, 0.959, 1.010, 1.153, 1.331, 1.542, 2.157, 3.181, 5.041]
    parameter0_title = 'Gly Con [%]'
    parameter1_title = 'Viscosity [mPa.s]'
    first_line =['filename', parameter0_title, parameter1_title, 'Debye time [us]', 'Neel time [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 10 THEN BEGIN
    GlyCon = [1.,3.,5.,10.,20.,30.,40.,50.,60.,70.]
    VisArray = [0.912, 0.959, 1.010, 1.153, 1.542, 2.157, 3.181, 5.041, 8.823, 17.96]
    parameter0_title = 'Glycerol Con [%]'
    parameter1_title = 'Viscosity [mPa.s]'
    first_line =['filename', parameter0_title, parameter1_title, 'Debye time [us]', 'Neel time [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 20 THEN BEGIN
    FieldAmp = FINDGEN(nrfile)*0.5+0.5
    parameter0_title = 'Field Amplitude [mT]'
    x_array = FieldAmp
    first_line =['filename', parameter0_title, 'Debye time [us]', 'Neel time [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 4 THEN BEGIN
    IF STRPOS(a_file[0], '20220720') NE -1 THEN BEGIN
      TemArray = [4., 25., 40., 50.]
      parameter0_title = 'Temperature [degree]'
      x_array = TemArray
    ENDIF ELSE BEGIN
      GelArray = [5., 10., 20., 30.]
      parameter0_title = 'Gelatin Con [%]'
      x_array = GelArray
    ENDELSE
    first_line =['filename', parameter0_title, 'Debye time [us]', 'Neel time [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian Percentage [%]']
  ENDIF
  ;for temperature
  IF nrfile EQ 3 THEN BEGIN
    IF STRPOS(a_file[0], '20220720') NE -1 THEN BEGIN
      TemArray = [4., 25., 50.]
      parameter0_title = 'Temperature [degree]'
      x_array = TemArray
    ENDIF ELSE BEGIN
      TemArray = [25.,37.,52.]
      parameter0_title = 'Temperature [degree]'
      x_array = TemArray
    ENDELSE
    first_line =['filename', parameter0_title, 'Debye time [us]', 'Neel time [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian Percentage [%]']
  ENDIF

  ;write the 1st line
  printf, unit, FORMAT = '(10(A, %"\t"))', first_line

  N_tao = 1001
  data_array = DBLARR(nrfile+2, N_tao)*0.
  time_Debye = DBLARR(nrfile)*0.
  time_Neel = DBLARR(nrfile)*0.
  time_Brownian = DBLARR(nrfile)*0.
  ratio_Neel = DBLARR(nrfile)*0.
  ratio_Brownian = DBLARR(nrfile)*0.

  FOR i=0L, nrfile-1L DO BEGIN
    filename = FILE_BASENAME(a_file[i], '.txt')
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      sState.filename = filename
      sState.directory = FILE_DIRNAME(a_file[i])
      WIDGET_CONTROL, sState.wFilenameField, SET_VALUE=sState.filename+', '+STRTRIM(STRING(i+1),1)+'/'+STRTRIM(STRING(nrfile),1)+' analyzing...'
    ENDIF
    fit_para = read_pulsed_signal_dat(a_file[i], sState)
    time_Debye[i] = fit_para.tao_debye_mono
    time_Neel[i] = fit_para.bi_tao[1]
    time_Brownian[i] = fit_para.bi_tao[0]
    ratio_Neel[i] = fit_para.bi_tao_comp[1]
    ratio_Brownian[i] = fit_para.bi_tao_comp[0]
    IF nrfile EQ 3 OR nrfile EQ 20 OR nrfile EQ 4 THEN BEGIN
      printf, unit, FORMAT = '(10(A, %"\t"))', filename, x_array[i], time_Debye[i], time_Neel[i], ratio_Neel[i]*100., time_Brownian[i], ratio_Brownian[i]*100.
    ENDIF
    IF nrfile EQ 9 OR nrfile EQ 10 THEN BEGIN
      printf, unit, FORMAT = '(10(A, %"\t"))', filename, GlyCon[i], VisArray[i], time_Debye[i], time_Neel[i], ratio_Neel[i]*100., time_Brownian[i], ratio_Brownian[i]*100.
    ENDIF
  ENDFOR

  ;  ;GJ, 2022/2/16, write the results
  ;  output_filename = dir_name+'UpenEps_data.dat'
  ;  openw,unit,output_filename,/get_lun
  ;  FOR l=0, nrfile+1 DO BEGIN
  ;    printf, unit, REFORM(data_array[l, *])
  ;  ENDFOR
  ;  FREE_LUN, unit
  ;


  ; Scale the height for drawing zonal winds and load program colors.
  cgLoadCT, 0
  
  IF nrfile EQ 20 OR nrfile EQ 3 OR nrfile EQ 4 THEN BEGIN
    window, 0, title='Relaxation times'
    plot, x_array, time_Brownian, xtitle=parameter0_title, LINESTYLE=1, THICK=2
    oplot, x_array, time_Neel
    oplot, x_array, time_Debye, LINESTYLE=3, THICK=2
    window, 1, title='Ratio of relaxation time'
    plot, x_array, ratio_Neel, xtitle=parameter0_title, yrange=[0, 1]
    oplot, x_array, ratio_Brownian, LINESTYLE=1, THICK=2
    
    ;calculate the sensitivity
    fit_time_Debye = LINFIT(x_array, time_Debye)
    sensitivity_time_Debye = ABS(fit_time_Debye[1])/x_array[0]*100.
    print, 'Debye time sensitivity = ', sensitivity_time_Debye
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Debye time sensitivity = ', sensitivity_time_Debye

    
    fit_time_Brownian = LINFIT(x_array, time_Brownian)
    sensitivity_time_Brownian = ABS(fit_time_Brownian[1])/x_array[0]*100.
    print, 'Brownian time sensitivity = ', sensitivity_time_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian time sensitivity = ', sensitivity_time_Brownian

    fit_time_Neel = LINFIT(x_array, time_Neel)
    sensitivity_time_Neel = ABS(fit_time_Neel[1])/x_array[0]*100.
    print, 'Neel time sensitivity = ', sensitivity_time_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel time sensitivity = ', sensitivity_time_Neel

    fit_ratio_Brownian = LINFIT(x_array, ratio_Brownian)
    sensitivity_ratio_Brownian = ABS(fit_ratio_Brownian[1])/x_array[0]*100.
    print, 'Brownian ratio sensitivity = ', sensitivity_ratio_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian ratio sensitivity = ', sensitivity_ratio_Brownian

    fit_ratio_Neel = LINFIT(x_array, ratio_Neel)
    sensitivity_ratio_Neel = ABS(fit_ratio_Neel[1])/x_array[0]*100.
    print, 'Neel ratio sensitivity (FieldAmp) = ', sensitivity_ratio_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel ratio sensitivity = ', sensitivity_ratio_Neel
  ENDIF
  IF nrfile EQ 9 OR nrfile EQ 10 THEN BEGIN
    window, 0, title='Relaxation times'
    plot, VisArray, time_Brownian, xtitle='Viscosity [mPa.s]', ytitle='Relaxation time [us]', LINESTYLE=1, THICK=2
    oplot, VisArray, time_Neel
    oplot, VisArray, time_Debye, LINESTYLE=3, THICK=2
    window, 1, title='Ratio of relaxation time'
    plot, VisArray, ratio_Neel, xtitle='Viscosity [mPa.s]', ytitle='Ratio [%]', yrange=[0, 1]
    oplot, VisArray, ratio_Brownian, LINESTYLE=1, THICK=2
    GlyConColors = BytScl(GlyCon, Min=0, Max=MAX(GlyCon)+10, Top=MAX(GlyCon)+10)

    ;calculate the sensitivity
    fit_time_Debye = LINFIT(VisArray, time_Debye)
    sensitivity_time_Debye = ABS(fit_time_Debye[1])/VisArray[0]*100.
    print, 'Debye time sensitivity (Viscosity) = ', sensitivity_time_Debye
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Debye time sensitivity (Viscosity) = ', sensitivity_time_Debye

    ;calculate the sensitivity
    fit_time_Brownian = LINFIT(VisArray, time_Brownian)
    sensitivity_time_Brownian = ABS(fit_time_Brownian[1])/VisArray[0]*100.
    print, 'Brownian time sensitivity (Viscosity) = ', sensitivity_time_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian time sensitivity (Viscosity) = ', sensitivity_time_Brownian

    fit_time_Neel = LINFIT(VisArray, time_Neel)
    sensitivity_time_Neel = ABS(fit_time_Neel[1])/VisArray[0]*100.
    print, 'Neel time sensitivity (Viscosity) = ', sensitivity_time_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel time sensitivity (Viscosity) = ', sensitivity_time_Neel

    fit_ratio_Brownian = LINFIT(VisArray, ratio_Brownian)
    sensitivity_ratio_Brownian = ABS(fit_ratio_Brownian[1])/VisArray[0]*100.
    print, 'Brownian ratio sensitivity (Viscosity) = ', sensitivity_ratio_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian ratio sensitivity (Viscosity) = ', sensitivity_ratio_Brownian

    fit_ratio_Neel = LINFIT(VisArray, ratio_Neel)
    sensitivity_ratio_Neel = ABS(fit_ratio_Neel[1])/VisArray[0]*100.
    print, 'Neel ratio sensitivity (Viscosity) = ', sensitivity_ratio_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel ratio sensitivity (Viscosity) = ', sensitivity_ratio_Neel
  ENDIF

  ;close the unit
  FREE_LUN, unit

END

;GJ, 2022/2/15, calculate the UPEN eps file
;GJ, 2022/7/7, add sensitivity calculation
PRO batch_read_UPEN_eps, sState
  
  dir_name = 'C:\UpenWin\bin\MPI01\'
  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly1\eps\'
  IF N_ELEMENTS(sState) GT 0 THEN dir_name = DIALOG_PICKFILE(PATH=sState.directory, /MUST_EXIST, TITLE="eps files", /DIRECTORY)
  
  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  IF STRLEN(dir_name) EQ 0 THEN RETURN
  ;@GJ, 2022/7/2, update the directory
  IF N_ELEMENTS(sState) GT 0 THEN sState.directory = dir_name
  
  ;GJ, 2022/2/6, redo the fitting and write the results
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  output_filename_temp = dir_name+current_time+'_sensitivity.xls'

  ;Select a file to save exported data, the default opened directory is the one you select just now
  output_filename=DIALOG_PICKFILE(FILE=output_filename_temp, FILTER = ['*.xls'], title='Save Your Output as EXCEL file', PATH=dir_name)
  ;// open an ASCII file for output
  ;IF STRCMP(STRMID(aSaveFile, 2, 3, /REVERSE_OFFSET), 'txt',/FOLD_CASE)   THEN BEGIN
  IF (STRLEN(output_filename) EQ 0)  THEN BEGIN
    infowarning = DIALOG_MESSAGE('No output excel file selected!', /ERROR)
    RETURN
  ENDIF
  
  ;write file
  openw,unit,output_filename,/get_lun
  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.eps', count=nrfile)
  IF nrfile EQ 0 THEN BEGIN
    infowarning = DIALOG_MESSAGE('No eps file found!', /ERROR)
    RETURN
  ENDIF
  
  ;@GJ, 2022/7/7, define the title 
  IF nrfile EQ 9 THEN BEGIN
    GlyCon = [1.,3.,5.,10.,15.,20.,30.,40.,50.]
    VisArray = [0.912, 0.959, 1.010, 1.153, 1.331, 1.542, 2.157, 3.181, 5.041]
    parameter0_title = 'Gly Con [%]'
    parameter1_title = 'Viscosity [mPa.s]'
    first_line =['filename', parameter0_title, parameter1_title, 'Neel time [us]', 'Neel FWHM [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian FWHM [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 10 THEN BEGIN
    GlyCon = [1.,3.,5.,10.,20.,30.,40.,50.,60.,70.]
    VisArray = [0.912, 0.959, 1.010, 1.153, 1.542, 2.157, 3.181, 5.041, 8.823, 17.96]
    parameter0_title = 'Gly Con [%]'
    parameter1_title = 'Viscosity [mPa.s]'
    first_line =['filename', parameter0_title, parameter1_title, 'Neel time [us]', 'Neel FWHM [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian FWHM [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 20 THEN BEGIN
    FieldAmp = FINDGEN(nrfile)*0.5+0.5
    parameter0_title = 'Field Amplitude [mT]'
    x_array = FieldAmp
    first_line =['filename', parameter0_title, 'Neel time [us]', 'Neel FWHM [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian FWHM [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 3 THEN BEGIN
    IF STRPOS(a_file[0], '20220720') NE -1 THEN BEGIN
      TemArray = [4., 25., 50.]
      parameter0_title = 'Temperature [degree]'
      x_array = TemArray
    ENDIF ELSE BEGIN
      TemArray = [25.,37.,52.]
      inverse_TemKArray = 1./(273.15 + TemArray); 1/K
      parameter0_title = 'Temperature [degree]'
      x_array = TemArray
    ENDELSE
    first_line =['filename', parameter0_title, 'Neel time [us]', 'Neel FWHM [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian FWHM [us]', 'Brownian Percentage [%]']
  ENDIF
  IF nrfile EQ 4 THEN BEGIN
    IF STRPOS(a_file[0], '20220720') NE -1 THEN BEGIN
      TemArray = [4., 25., 40., 50.]
      parameter0_title = 'Temperature [degree]'
      x_array = TemArray
    ENDIF ELSE BEGIN
      GelArray = [5., 10., 20., 30.]
      parameter0_title = 'Gelatin Con [%]'
      x_array = GelArray
    ENDELSE
    first_line =['filename', parameter0_title, 'Neel time [us]', 'Neel FWHM [us]', 'Neel Percentage [%]', 'Brownian time [us]', 'Brownian FWHM [us]', 'Brownian Percentage [%]']
  ENDIF
  
  ;write the 1st line
  printf, unit, FORMAT = '(10(A, %"\t"))', first_line
  
  data_array_para = read_UPEN_eps(a_file[0], sState)
  N_tao = data_array_para.N_times
  data_array = DBLARR(nrfile+2, N_tao)*0.
  time_Neel = DBLARR(nrfile)*0.
  time_Brownian = DBLARR(nrfile)*0.
  FWHM_Neel = DBLARR(nrfile)*0.
  FWHM_Brownian = DBLARR(nrfile)*0.
  ratio_Neel = DBLARR(nrfile)*0.
  ratio_Brownian = DBLARR(nrfile)*0.

  FOR i=0L, nrfile-1L DO BEGIN
    filename = FILE_BASENAME(a_file[i], '.eps')
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      sState.filename = filename
      sState.directory = FILE_DIRNAME(a_file[i])
      WIDGET_CONTROL, sState.wFilenameField, SET_VALUE=sState.filename+', '+STRTRIM(STRING(i+1),1)+'/'+STRTRIM(STRING(nrfile),1)+' analyzing...'
    ENDIF
    ;@GJ, 2022/4/13, replace the old program with new function
    data_array_para = read_UPEN_eps(a_file[i], sState)
    
    IF i EQ 0 THEN data_array[0,*] = REFORM(data_array_para.times[0:N_tao-1])
    data_array[i+2,*] = REFORM(data_array_para.amp_times[0:N_tao-1])
    time_Neel[i] = data_array_para.time_Neel
    time_Brownian[i] = data_array_para.time_Brownian
    FWHM_Neel[i] = data_array_para.FWHM_Neel
    FWHM_Brownian[i] = data_array_para.FWHM_Brownian
    ratio_Neel[i] = data_array_para.ratio_Neel
    ratio_Brownian[i] = data_array_para.ratio_Brownian
    IF nrfile EQ 20 OR nrfile EQ 3 OR nrfile EQ 4 THEN BEGIN
      printf, unit, FORMAT = '(10(A, %"\t"))', filename, x_array[i], time_Neel[i], FWHM_Neel[i], ratio_Neel[i]*100., time_Brownian[i], FWHM_Brownian[i], ratio_Brownian[i]*100.
    ENDIF
    IF nrfile EQ 9 OR nrfile EQ 10 THEN BEGIN
      printf, unit, FORMAT = '(10(A, %"\t"))', filename, GlyCon[i], VisArray[i], time_Neel[i], FWHM_Neel[i], ratio_Neel[i]*100., time_Brownian[i], FWHM_Brownian[i], ratio_Brownian[i]*100.
    ENDIF
    
    ;@GJ, 2022/10/01, read Jiaming' PDCO file
    filename_hist = dir_name+filename+'_hist.txt'
    file_ex = FILE_TEST(filename_hist)
    IF file_ex EQ 1 THEN BEGIN
      i_temp=0
      OPENR, unit1, filename_hist, /GET_LUN
      temp=''
      READF, unit1, temp, FORMAT='(%"%s")'
      WHILE ~ EOF(unit1) DO BEGIN
        READF, unit1, temp, FORMAT='(%"%s")'
        temp_spl = STRSPLIT(temp, ' ', count=n_spl, /EXTRACT)
        ;save the result
        IF i_temp EQ 0 THEN BEGIN
          data_array_temp = TRANSPOSE(DOUBLE(temp_spl))
        ENDIF ELSE BEGIN
          data_array_temp = [data_array_temp, TRANSPOSE(DOUBLE(temp_spl))]
        ENDELSE
        i_temp++
      ENDWHILE
      FREE_LUN, unit1
      ;@GJ, 2022/10/1, rescale
      data_array_temp[*,1] = data_array_temp[*,1]/MAX(data_array_temp[*,1])*100. 
      ;@GJ, 2022/7/6, calculate the relaxation times spectrum by interpolation
      amp_times = INTERPOL(REFORM(data_array_temp[*,1]), REFORM(data_array_temp[*,0]), data_array_para.times)
      data_array[i+2,*] = REFORM(amp_times[0:N_tao-1])
      iplot, data_array_para.times, amp_times, xtitle='relaxation time [us]', ytitle='amp'
    ENDIF
ENDFOR
      
;  ;GJ, 2022/2/16, write the results
;  output_filename = dir_name+'UpenEps_data.dat'
;  openw,unit,output_filename,/get_lun
;  FOR l=0, nrfile+1 DO BEGIN
;    printf, unit, REFORM(data_array[l, *])
;  ENDFOR
;  FREE_LUN, unit
;  

  
  ; Get dimension of wind data set.
  tao_distr = data_array[2:*,*]
  dims = Size(tao_distr, /Dimensions)

  ; Scale the height for drawing zonal winds and load program colors.
  cgLoadCT, 0
  IF nrfile EQ 20 OR nrfile EQ 3 OR nrfile EQ 4 THEN BEGIN
    window, 1, title='Relaxation times'
    plot, x_array, time_Brownian, xtitle='Field Amplitude [mT]', LINESTYLE=1, THICK=2
    oplot, x_array, time_Neel
    window, 2, title='FWHM of relaxation time'
    plot, x_array, FWHM_Neel, xtitle='Field Amplitude [mT]', yrange=[0, MAX([FWHM_Neel, FWHM_Brownian])]
    oplot, x_array, FWHM_Brownian, LINESTYLE=1, THICK=2
    window, 3, title='Ratio of relaxation time'
    plot, x_array, ratio_Neel, xtitle='Field Amplitude [mT]', yrange=[0, 1]
    oplot, x_array, ratio_Brownian, LINESTYLE=1, THICK=2
    FieldAmpColors = BytScl(x_array, Min=0, Max=MAX(x_array)*2+1, Top=MAX(x_array)*2.+1)
    cgLoadCT, 33, NColors=MAX(x_array)*2+1
    ;calculate the sensitivity
    fit_time_Neel = LINFIT(x_array, time_Neel)
    sensitivity_time_Neel = ABS(fit_time_Neel[1])/x_array[0]*100.
    print, 'Neel time sensitivity = ', sensitivity_time_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel time sensitivity = ', sensitivity_time_Neel
    
    fit_time_Brownian = LINFIT(x_array, time_Brownian)
    sensitivity_time_Brownian = ABS(fit_time_Brownian[1])/x_array[0]*100.
    print, 'Brownian time sensitivity = ', sensitivity_time_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian time sensitivity = ', sensitivity_time_Brownian

    fit_ratio_Neel = LINFIT(x_array, ratio_Neel)
    sensitivity_ratio_Neel = ABS(fit_ratio_Neel[1])/x_array[0]*100.
    print, 'Neel ratio sensitivity = ', sensitivity_ratio_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel ratio sensitivity = ', sensitivity_ratio_Neel
    
    fit_ratio_Brownian = LINFIT(x_array, ratio_Brownian)
    sensitivity_ratio_Brownian = ABS(fit_ratio_Brownian[1])/x_array[0]*100.
    print, 'Brownian ratio sensitivity = ', sensitivity_ratio_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian ratio sensitivity = ', sensitivity_ratio_Brownian
    
    fit_FWHM_Neel = LINFIT(x_array, FWHM_Neel)
    sensitivity_FWHM_Neel = ABS(fit_FWHM_Neel[1])/x_array[0]*100.
    print, 'Neel FWHM sensitivity = ', sensitivity_FWHM_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel FWHM sensitivity = ', sensitivity_FWHM_Neel
    
    fit_FWHM_Brownian = LINFIT(x_array, FWHM_Brownian)
    sensitivity_FWHM_Brownian = ABS(fit_FWHM_Brownian[1])/x_array[0]*100.
    print, 'Brownian FWHM sensitivity = ', sensitivity_FWHM_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian FWHM sensitivity = ', sensitivity_FWHM_Brownian
  ENDIF
  IF nrfile EQ 9 OR nrfile EQ 10 THEN BEGIN
    window, 1, title='Relaxation times'
    plot, VisArray, time_Brownian, xtitle='Viscosity [mPa.s]', ytitle='Relaxation time [us]', LINESTYLE=1, THICK=2
    oplot, VisArray, time_Neel
    window, 2, title='FWHM of relaxation time'
    plot, VisArray, FWHM_Neel, xtitle='Viscosity [mPa.s]', ytitle='FWHM [us]', yrange=[0, MAX([FWHM_Neel, FWHM_Brownian])]
    oplot, VisArray, FWHM_Brownian, LINESTYLE=1, THICK=2
    window, 3, title='Ratio of relaxation time'
    plot, VisArray, ratio_Neel, xtitle='Viscosity [mPa.s]', ytitle='Ratio [%]', yrange=[0, 1]
    oplot, VisArray, ratio_Brownian, LINESTYLE=1, THICK=2
    GlyConColors = BytScl(GlyCon, Min=0, Max=MAX(GlyCon)+10, Top=MAX(GlyCon)+10)
    cgLoadCT, 33, NColors=MAX(GlyCon)+10
    
    ;calculate the sensitivity
    fit_time_Brownian = LINFIT(VisArray, time_Brownian)
    sensitivity_time_Brownian = ABS(fit_time_Brownian[1])/VisArray[0]*100.
    print, 'Brownian time sensitivity (Viscosity) = ', sensitivity_time_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian time sensitivity (Viscosity) = ', sensitivity_time_Brownian
    
    fit_time_Neel = LINFIT(VisArray, time_Neel)
    sensitivity_time_Neel = ABS(fit_time_Neel[1])/VisArray[0]*100.
    print, 'Neel time sensitivity (Viscosity) = ', sensitivity_time_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel time sensitivity (Viscosity) = ', sensitivity_time_Neel
    
    fit_FWHM_Brownian = LINFIT(VisArray, FWHM_Brownian)
    sensitivity_FWHM_Brownian = ABS(fit_FWHM_Brownian[1])/VisArray[0]*100.
    print, 'Brownian FWHM sensitivity (Viscosity) = ', sensitivity_FWHM_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian FWHM sensitivity (Viscosity) = ', sensitivity_FWHM_Brownian

    fit_FWHM_Neel = LINFIT(VisArray, FWHM_Neel)
    sensitivity_FWHM_Neel = ABS(fit_FWHM_Neel[1])/VisArray[0]*100.
    print, 'Neel FWHM sensitivity (Viscosity) = ', sensitivity_FWHM_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel FWHM sensitivity (Viscosity) = ', sensitivity_FWHM_Neel

    fit_ratio_Brownian = LINFIT(VisArray, ratio_Brownian)
    sensitivity_ratio_Brownian = ABS(fit_ratio_Brownian[1])/VisArray[0]*100.
    print, 'Brownian ratio sensitivity (Viscosity) = ', sensitivity_ratio_Brownian
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Brownian ratio sensitivity (Viscosity) = ', sensitivity_ratio_Brownian

    fit_ratio_Neel = LINFIT(VisArray, ratio_Neel)
    sensitivity_ratio_Neel = ABS(fit_ratio_Neel[1])/VisArray[0]*100.
    print, 'Neel ratio sensitivity (Viscosity) = ', sensitivity_ratio_Neel
    printf, unit, FORMAT = '(10(A, %"\t"))', 'Neel ratio sensitivity (Viscosity) = ', sensitivity_ratio_Neel
  ENDIF
  
  ;close the unit
  FREE_LUN, unit
  
  ; Open a display window.
  cgDisplay

  ; Set up plotting parameters.
  thick = (!D.Name EQ 'PS') ? 6 : 2
  plotPosition = [0.15, 0.125, 0.95, 0.75]

  ; Draw the initial plot, no data.
  ;plot the result
  ; Reformat pressure and height as row vectors.
  tao_array = Reform(data_array[0,*])
  IF nrfile EQ 20 THEN BEGIN
    ; Draw a gray background so colors stand out more.
    cgColorFill, Position=plotPosition, Color='gray'
    cgPlot, tao_array[0:1000], tao_distr[0, 0:1000], XRange=[0,100], YRange=[0.01,150],YStyle=1, YLOG=1, $
      YTitle='Amp [A.U.]', XTitle='Relaxation time [us]', XStyle=1, $
      XTicks=10, /NoData, /NoErase, $
      Position=plotPosition, Label='Relaxation Spectral Analysis'
    ; Overplot the zonal winds in colors scaled by height.
    FOR j=0,N_ELEMENTS(x_array)-1 DO cgPlotS, tao_array[0:1000], tao_distr[j,0:1000]+(j+2.)*10.^(-2), tao_array[0:1000], Color=ROUND(x_array[j])*2, Thick=thick
    ; Add a color bar.
    cgColorbar, NColors=MAX(x_array)*2+1, Range=[0,MAX(x_array)], Title=parameter0_title, $
      TLocation='top', Format='(i0)', Position=[0.125, 0.86, 0.9, 0.90]
  ENDIF
  IF nrfile EQ 3 OR nrfile EQ 4 THEN BEGIN
    ; Draw a gray background so colors stand out more.
    cgColorFill, Position=plotPosition, Color='gray'
    cgPlot, tao_array[0:1000], tao_distr[0, 0:1000], XRange=[0,100], YRange=[0.01,150],YStyle=1, YLOG=1, $
      YTitle='Amp [A.U.]', XTitle='Relaxation time [us]', XStyle=1, $
      XTicks=10, /NoData, /NoErase, $
      Position=plotPosition, Label='Relaxation Spectral Analysis'
    ; Overplot the zonal winds in colors scaled by height.
    FOR j=0,N_ELEMENTS(x_array)-1 DO cgPlotS, tao_array[0:1000], tao_distr[j,0:1000]+(j+2.)*10.^(-2), tao_array[0:1000], Color=ROUND(x_array[j])*2, Thick=thick
    ; Add a color bar.
    cgColorbar, NColors=MAX(x_array)*2+1, Range=[0,MAX(x_array)], Title=parameter0_title, $
      TLocation='top', Format='(i0)', Position=[0.125, 0.86, 0.9, 0.90]
  ENDIF
  IF nrfile EQ 9 OR nrfile EQ 10 THEN BEGIN
    ; Draw a gray background so colors stand out more.
    cgColorFill, Position=plotPosition, Color='gray'
    cgPlot, tao_array[0:1000], tao_distr[0, 0:1000], XRange=[0,100], YRange=[0.01,150],YStyle=1, YLOG=1, $
      YTitle='Amp [A.U.]', XTitle='Relaxation time [us]', XStyle=1, $
      XTicks=10, /NoData, /NoErase, $
      Position=plotPosition, Label='Relaxation Spectral Analysis'; (' + filename + ')'
    ; Overplot the zonal winds in colors scaled by height.
    FOR j=0,N_ELEMENTS(GlyConColors)-1 DO cgPlotS, tao_array[0:1000], tao_distr[j,0:1000]+(j+2.)*10.^(-2), tao_array[0:1000], Color=GlyConColors[j], Thick=thick
    ; Add a color bar.
    ;GJ, 2022/3/29, correct the spelling of glycerol
    cgColorbar, NColors=MAX(GlyCon)+10, Range=[0,MAX(GlyCon)+10], Title='Glycerol Concentration [%]', $
      TLocation='top', Format='(i0)', Position=[0.125, 0.86, 0.9, 0.90]
  ENDIF

  
  ;  ;GJ, 2022/2/16, write the results
  ;  output_filename = dir_name+'UpenEps_data.dat'
END

;@GJ, 2022/4/13, modify the versions
;data_array = read_UPEN_eps()
FUNCTION read_UPEN_eps, filename, sState
  ;GJ, 2022/6/19, clean the plot in sState
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    new_pic = BYTARR(3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF

  ;GJ, 2022/3/8, nonzero gilename
  IF (N_ELEMENTS(filename) EQ 0) THEN BEGIN
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      filename = DIALOG_PICKFILE(PATH=sState.directory, /MUST_EXIST, FILTER='*.eps', TITLE="eps file")
    ENDIF ELSE BEGIN
      filename = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\Xinfeng\20220616\', /MUST_EXIST, FILTER='*.eps', TITLE="eps file")
    ENDELSE
    IF (STRLEN(filename) EQ 0) THEN BEGIN
      info=dialog_message('No file name!', /information)
      RETURN, 0
    ENDIF
  ENDIF ELSE BEGIN
     IF (STRLEN(filename) EQ 0) THEN BEGIN
       IF N_ELEMENTS(sState) GT 0 THEN BEGIN
         filename = DIALOG_PICKFILE(PATH=sState.directory, /MUST_EXIST, FILTER='*.eps', TITLE="eps file")
       ENDIF ELSE BEGIN
         filename = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\Xinfeng\20220616\', /MUST_EXIST, FILTER='*.eps', TITLE="eps file")
       ENDELSE
       IF (STRLEN(filename) EQ 0) THEN BEGIN
         info=dialog_message('No file name!', /information)
         RETURN, 0
       ENDIF
     ENDIF
  ENDELSE
  
  IF STRLEN(filename) GT 1 THEN BEGIN
    filename_info = FILE_INFO(filename)
    IF (filename_info.exists EQ 0) THEN BEGIN
      info=dialog_message('Wrong file name!', /information)
      RETURN, 0
    ENDIF
    IF (STRLEN(filename) EQ 0) THEN BEGIN
      info=dialog_message('No file name!', /information)
      RETURN, 0
    ENDIF
  ENDIF
  
  ;replace the sstate.directory
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    sState.directory = FILE_DIRNAME(filename)
    sState.filename = FILE_BASENAME(filename)
  ENDIF
  
  OPENR, lun, filename, /get_lun
  
  ;read the first 1042 lines, extract the position of log coordinates
  temp_str=''
  str_arr = ['(.01)', '(.1)', '(1)', '(1O)', '(1OO)']
  log_N = N_ELEMENTS(str_arr)
  coor_arr = DBLARR(log_N)*0.
  for k=0,1042 do begin
    readf, lun, temp_str
    for i=0, log_N-1 do begin
      result = STRPOS(temp_str, str_arr[i]+' stringwidth pop')
      IF result NE -1 THEN BEGIN
        print, 'found', i, ': ', temp_str
        coor_arr[i] = (STRSPLIT(temp_str, ' ', /EXTRACT))[0]
      ENDIF
    endfor
  endfor
  print, coor_arr
  
  N_tao = 100
  data_x = DBLARR(N_tao)*0.
  data_y = DBLARR(N_tao)*0.
  ;starting from 1043, data in black
  for j=0,N_tao-1 do begin
    READF, lun, temp1, temp2, FORMAT='(%"%f %f")'
    data_x[j] = temp1
    data_y[j] = temp2
  endfor
  FREE_LUN, lun
  
  ;GJ, 2022/4/13, to rescale and get the correct time
  tao_x = 10.^((data_x-coor_arr[2])/(coor_arr[4]-coor_arr[2])*2.)
  ;GJ, 2022/7/14, reverse the orders
  tao_x = REVERSE(tao_x, /OVERWRITE)
  data_y = REVERSE(data_y, /OVERWRITE)
  
  ;GJ, 2022/7/14, extract the time featuers
  ;extract Neel from 0-20us
  NB_sepind = WHERE(tao_x LT 20, NB_separa)
  max_N = MAX(data_y[0:NB_separa], maxind_N)
  time_Neel = tao_x[maxind_N]
  half_find = WHERE(data_y[0:NB_separa] GE max_N/2., count_N)
  FWHM_Neel = count_N
  NB_separa = (maxind_N + FWHM_Neel*3.)
  AUC_Neel = 0.
  FOR i=0, NB_separa DO AUC_Neel += data_y[i]*(tao_x[i+1]-tao_x[i])
  
  ;extract Brownian from NB_separa to 150 ms
  max_B = MAX(data_y[NB_separa:*], maxind_B)
  time_Brownian = tao_x[maxind_B+NB_separa]
  half_find = WHERE(data_y[NB_separa:*] GE max_B/2., count_B)
  IF count_B LT 5 THEN BEGIN
    FWHM_Brownian = 0
  ENDIF ELSE BEGIN
    FWHM_Brownian = count_B*(tao_x[maxind_B+NB_separa+1]-tao_x[maxind_B+NB_separa-1])/2.
  ENDELSE
  AUC_Brownian = 0
  FOR i=0, N_ELEMENTS(data_y)-NB_separa-2 DO AUC_Brownian += data_y[NB_separa+i]*(tao_x[NB_separa+i+1]-tao_x[NB_separa+i])

  ratio_Neel = AUC_Neel/(AUC_Neel + AUC_Brownian)
  ratio_Brownian = AUC_Brownian/(AUC_Neel + AUC_Brownian)
  IF ratio_Brownian LT 0.001 THEN time_Brownian = 0.
  print, 'Neel: ', time_Neel, ' [us] with FWHM ', ROUND(FWHM_Neel), ' us with ', ratio_Neel*100., '%'
  print, 'Brownian: ', time_Brownian, ' [us] with FWHM ', ROUND(FWHM_Brownian), ' us with ', ratio_Brownian*100., '%'
  ;
  ;plot the eps file
;  max_y = MAX(data_y, maxind)
  filename_base = FILE_BASENAME(filename, '.eps')
;  iplot, tao_x, data_y, XTITLE='tau [us]', YTITLE='Amp [A.U.]', XRANGE=[0,100], YRANGE=[0,100], TITLE='Multi-Exponential Analysis (' + filename_base + ')'
  IF OBJ_VALID(p) THEN p.close
  p = PLOT(tao_x, data_y, 'b-2', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[0,150], XLOG=1, YRANGE=[0,100], $; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='Multi-Exponential Analysis (' + filename_base + ')', XTITLE='tau [us]', YTITLE='Amp [A.U.]')
  p.name = ' Spectrum'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.7])
  result00 = p.CopyWindow(BORDER=10, RESOLUTION=300)

  ;display in the main window
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    ;save as png file
    WRITE_PNG, sState.directory+'\'+filename_base+'_eps.png', result00
    new_pic = CONGRID(result00, 3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF ELSE BEGIN
    ;save as png file
    filename_dir = FILE_DIRNAME(filename)
    WRITE_PNG, filename_dir+'\'+filename_base+'_eps.png', result00
  ENDELSE
  ;close the plot
  IF OBJ_VALID(p) THEN p.close
  
  ;@GJ, 2022/7/14, validate the spectrum
  ;load original decay signal
  fn_len = STRLEN(filename_base)
  filename_dir = FILE_DIRNAME(filename)
  file_001 = filename_dir+'\'+STRMID(filename_base, 0, fn_len-1)+'.001'
  data_array = read_UPEN_001_file(file_001)
  time = REFORM(data_array[*,0])
  original_signal = REFORM(data_array[*,1])
  window, 4
  noise = MEAN(original_signal[N_ELEMENTS(time)/2-1-N_ELEMENTS(time)/50:N_ELEMENTS(time)/2-1])
  plot, time, original_signal

  ;N_tao syn
  syn_tao_signal = original_signal*0.
  FOR i=1, N_tao-1 DO syn_tao_signal += data_y[i]*EXP(-time/tao_x[i])
  syn_tao_signal = syn_tao_signal/MAX(syn_tao_signal)*MAX(original_signal)
  oplot, time, syn_tao_signal, LINESTYLE=3, thick=2
  
  ;@GJ, 2022/7/6, calculate the relaxation times spectrum by interpolation
  N_times = 1501.
  times = FINDGEN(N_times)/10.
  amp_times = INTERPOL(data_y, tao_x, times)
  amp_times[0] = 0
  
  ;only plot the time from 0 to 100 us
  data_array_para = { $
    tao_x: tao_x, $
    data_y: data_y, $
    N_times: N_times, $
    times: times, $
    amp_times: amp_times, $
    time_Neel: time_Neel, $
    time_Brownian: time_Brownian, $
    FWHM_Neel: FWHM_Neel, $
    FWHM_Brownian: FWHM_Brownian, $
    ratio_Neel: ratio_Neel, $
    ratio_Brownian: ratio_Brownian $
    }
  
  ;write file
  N_int = 10000.
  time_array = FINDGEN(N_int)/100.
  signal_array = FINDGEN(N_int)*0.
  FOR i=0, N_ELEMENTS(tao_x)-1 DO BEGIN
    IF data_y[i] GT 0.01 THEN BEGIN
      temp_signal = data_y[i]*exp(-time_array/tao_x[i])
      signal_array += temp_signal
    ENDIF
  ENDFOR
  ;add noise to signal_array
  signal_array_noise = signal_array*0.
  FOR k=0, N_ELEMENTS(signal_array)-1 DO BEGIN
    seed = !NULL
    randomValue = RANDOMU(seed)/10.
    signal_array_noise[k] = signal_array[k]*(1.+randomValue)
  ENDFOR
;  signal_array_noise = NOISE_HURL(signal_array, 0.001)
  
;  iplot, time_array, signal_array, title='decay signal plot w/o noise'
;  iplot, time_array, signal_array_noise, /overplot
;  file_dir = FILE_DIRNAME(filename)
;  write_UPEN_001_file, file_dir+'\'+filename_base+'highresN.001', time_array[10:*], signal_array_noise[10:*]
  RETURN, data_array_para
END

;@GJ, 2022/7/27, analyze signal system matrix based on the number of field strengths
PRO read_signal_sav_lr
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="sav file", FILTER='*.sav')
  RESTORE, sav_name
  ;GJ, 2022/7/27, save the image
  time_rise_flat = REFORM(signal_matrix[0, *])
  N_fields = N_ELEMENTS(signal_matrix[*,0])-1
  N_tp = N_ELEMENTS(signal_matrix[0,*])
  window, 0, title='signal with 400 averages'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix[i, *]), YRANGE=[MIN(signal_matrix[1:*, *]), MAX(signal_matrix[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix[i, *])
  ENDFOR
 
  window, 1, title='signal with single acquisition'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix_1ave[i, *]), YRANGE=[MIN(signal_matrix_1ave[1:*, *]), MAX(signal_matrix_1ave[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix_1ave[i, *])
  ENDFOR
  
  window, 2, title='signal with 10 averages'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix_10ave[i, *]), YRANGE=[MIN(signal_matrix_10ave[1:*, *]), MAX(signal_matrix_10ave[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix_10ave[i, *])
  ENDFOR
  
  window, 3, title='signal with 50 averages'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix_50ave[i, *]), YRANGE=[MIN(signal_matrix_50ave[1:*, *]), MAX(signal_matrix_50ave[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix_50ave[i, *])
  ENDFOR
 
  window, 4, title='signal with 100 averages'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix_100ave[i, *]), YRANGE=[MIN(signal_matrix_100ave[1:*, *]), MAX(signal_matrix_100ave[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix_100ave[i, *])
  ENDFOR
 
  window, 5, title='signal with 200 averages'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix_200ave[i, *]), YRANGE=[MIN(signal_matrix_200ave[1:*, *]), MAX(signal_matrix_200ave[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix_200ave[i, *])
  ENDFOR
 
  window, 6, title='signal with 300 averages'
  FOR i=1, N_fields DO BEGIN
    IF i EQ 1 THEN plot, time_rise_flat, REFORM(signal_matrix_300ave[i, *]), YRANGE=[MIN(signal_matrix_300ave[1:*, *]), MAX(signal_matrix_300ave[1:*, *])] ELSE oplot, time_rise_flat, REFORM(signal_matrix_300ave[i, *])
  ENDFOR


  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image_lr = CONGRID(REFORM(phantom_image[0,*,*]), N_fields, N_fields)
  
;  image_lr = image_lr*0.
;  image_lr[N_fields/2, N_fields/2] = 255
  iimage, image_lr, title='original image'
  
  signal_2d = DBLARR(N_fields, N_tp)*0.
  signal_2d_1ave = DBLARR(N_fields, N_tp)*0.
  signal_2d_10ave = DBLARR(N_fields, N_tp)*0.
  signal_2d_50ave = DBLARR(N_fields, N_tp)*0.
  signal_2d_100ave = DBLARR(N_fields, N_tp)*0.
  signal_2d_200ave = DBLARR(N_fields, N_tp)*0.
  signal_2d_300ave = DBLARR(N_fields, N_tp)*0.
  FOR i=0, N_fields-1 DO BEGIN
    signal_temp = DBLARR(N_tp)*0.
    image_1d = REFORM(image_lr[i, *])
    FOR j=0, N_fields-1 DO BEGIN
      signal_2d[i, *] += image_1d[j] * REFORM(signal_matrix[j+1, *])
      signal_2d_1ave[i, *] += image_1d[j] * REFORM(signal_matrix_1ave[j+1, *])
      signal_2d_10ave[i, *] += image_1d[j] * REFORM(signal_matrix_10ave[j+1, *])
      signal_2d_50ave[i, *] += image_1d[j] * REFORM(signal_matrix_50ave[j+1, *])
      signal_2d_100ave[i, *] += image_1d[j] * REFORM(signal_matrix_100ave[j+1, *])
      signal_2d_200ave[i, *] += image_1d[j] * REFORM(signal_matrix_200ave[j+1, *])
      signal_2d_300ave[i, *] += image_1d[j] * REFORM(signal_matrix_300ave[j+1, *])
    ENDFOR
  ENDFOR
  
  cur_dir = FILE_DIRNAME(sav_name)
  filename_System = cur_dir+'\systemMatrix_2d_array.dat'
  ;write file
  offset=10
  openw,unit,filename_System,/get_lun
  FOR i=0, N_fields-1 DO BEGIN
    FOR j=0, N_fields-1 DO BEGIN
      ;; write the data points
      WRITEU, unit, (signal_matrix[1:*, 0+offset:offset+N_fields-1])[i, j]
    ENDFOR
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit
  
  filename_System = cur_dir+'\signal_1d_array.dat'
  ;write file
  openw,unit,filename_System,/get_lun
  FOR i=0, N_fields-1 DO BEGIN
    j=6
 ;   FOR j=0, N_fields-1 DO BEGIN
      ;; write the data points
      WRITEU, unit, (REFORM(signal_2d[*, 0+offset:offset+N_fields-1]))[i, j]
 ;   ENDFOR
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit
  
  iimage, signal_2d, title='signal'
  offset=10
  Gradient_invert = LA_INVERT(signal_matrix[1:*, 0+offset:offset+N_fields-1], /DOUBLE, status =status)
  result = Gradient_invert ## (signal_matrix[1:*, 0+offset:offset+N_fields-1])
  iimage, result, title='inverse * ori'
  rec_image_lr = image_lr*0.
  rec_image_lr_1ave = image_lr*0.
  rec_image_lr_10ave = image_lr*0.
  rec_image_lr_50ave = image_lr*0.
  rec_image_lr_100ave = image_lr*0.
  rec_image_lr_200ave = image_lr*0.
  rec_image_lr_300ave = image_lr*0.
  FOR i=0, N_fields-1 DO BEGIN
    rec_image_lr[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d[i, 0+offset:offset+N_fields-1]))); REFORM(Gradient_invert ## REFORM(signal_2d[i, 0+offset:offset+N_fields-1]))
    rec_image_lr_1ave[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
    rec_image_lr_10ave[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d_10ave[i, 0+offset:offset+N_fields-1]))); REFORM(Gradient_invert ## REFORM(signal_2d_10ave[i, 0+offset:offset+N_fields-1]))
    rec_image_lr_50ave[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d_50ave[i, 0+offset:offset+N_fields-1])));REFORM(Gradient_invert ## REFORM(signal_2d_50ave[i, 0+offset:offset+N_fields-1]))
    rec_image_lr_100ave[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d_100ave[i, 0+offset:offset+N_fields-1])));REFORM(Gradient_invert ## REFORM(signal_2d_100ave[i, 0+offset:offset+N_fields-1]))
    rec_image_lr_200ave[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d_200ave[i, 0+offset:offset+N_fields-1])));REFORM(Gradient_invert ## REFORM(signal_2d_200ave[i, 0+offset:offset+N_fields-1]))
    rec_image_lr_300ave[i, *] = REFORM(art_func(signal_matrix[1:*, 0+offset:offset+N_fields-1], REFORM(signal_2d_300ave[i, 0+offset:offset+N_fields-1])));REFORM(Gradient_invert ## REFORM(signal_2d_300ave[i, 0+offset:offset+N_fields-1]))
  ENDFOR
  
;  iimage, rec_image_lr, title='from ART'
  all_rec_image_lr = DBLARR(N_fields*8+10*7, N_fields)
  all_rec_image_lr[0:N_fields-1, *] = rec_image_lr_1ave
  all_rec_image_lr[N_fields+10:2*N_fields-1+10, *] = rec_image_lr_10ave
  all_rec_image_lr[2*N_fields+20:3*N_fields-1+20, *] = rec_image_lr_50ave
  all_rec_image_lr[3*N_fields+30:4*N_fields-1+30, *] = rec_image_lr_100ave
  all_rec_image_lr[4*N_fields+40:5*N_fields-1+40, *] = rec_image_lr_200ave
  all_rec_image_lr[5*N_fields+50:6*N_fields-1+50, *] = rec_image_lr_300ave
  all_rec_image_lr[6*N_fields+60:7*N_fields-1+60, *] = rec_image_lr
  all_rec_image_lr[7*N_fields+70:8*N_fields-1+70, *] = image_lr
  iimage, all_rec_image_lr, title='reconstructed image ART'
  
  
  ;@GJ, 2022/7/29, interpolate to 128 pixels
  N_steps = 128.
  image_hr = CONGRID(REFORM(phantom_image[0,*,*]), N_steps, N_steps)
  iimage, image_hr, title='ori (high resolution)'
  field_array = (FINDGEN(N_steps)+1.)/N_steps*10.
  signal_2d_hr = DBLARR(N_steps, N_tp)*0.
  signal_matrix_hr = CONGRID(signal_matrix[1:N_fields, *], N_steps, N_tp)
  signal_matrix_hr_1ave = CONGRID(signal_matrix_1ave[1:N_fields, *], N_steps, N_tp)
  FOR i=0, N_steps-1 DO BEGIN
    signal_temp = DBLARR(N_tp)*0.
    image_1d_hr = REFORM(image_hr[i, *])
    FOR j=0, N_steps-1 DO BEGIN
      signal_2d_hr[i, *] += image_1d_hr[j] * REFORM(signal_matrix_hr_1ave[j, *])
    ENDFOR
  ENDFOR
  system_matrix = signal_matrix_hr[*, 0+offset:offset+N_steps-1]
  Gradient_invert_hr = LA_INVERT(system_matrix, /DOUBLE, status =status)
;  Gradient_invert_hr = INVERT(system_matrix, /DOUBLE)
  IF status GT 0 THEN result=dialog_message('unsuccessful inverse cal', /ERROR)
  result_hr = Gradient_invert_hr ## system_matrix
  iimage, result_hr, title='inverse * ori (high resolution)'
  
  ;@GJ, 2022/7/29, reconstruct
  rec_image_hr = image_hr*0.
  FOR i=0, N_steps-1 DO rec_image_hr[i, *] = REFORM(art_func(system_matrix, REFORM(signal_2d_hr[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert_hr ## REFORM(signal_2d_hr[i, 0+offset:offset+N_steps-1]))
  iimage, rec_image_hr, title='reconstructed image ART (high resolution)'
  
  ; Compute the Singular Value Decomposition:
  SVDC, system_matrix, W, U, V
  iplot, alog(w), title='w log plot'
  sv_inv = system_matrix*0.
  FOR k=0, 20 DO sv_inv[K,K] = 1./W[K]
  rec_image_hr_svdc = image_hr*0.
  FOR i=0, N_steps-1 DO rec_image_hr_svdc[i, *] = REFORM(V ## sv_inv ## TRANSPOSE(U) ## REFORM(signal_2d_hr[i, 0+offset:offset+N_steps-1]))
  iimage, rec_image_hr_svdc, title='reconstructed image svdc (high resolution)'

END

;@GJ, 2022/11/18, orthogonal excitation
;@GJ, 2022/11/19, adjusting the field strength
PRO orthogonal_excitation_diff_H
  ;define the kernel
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  particle_size = 24.4924;PBS, [25.2579], 1% Glycerol
  relaxation_time = 25.;0.;0.5;5.;0.5;2.;5.;0.;5.;2. ;us
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
  ;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  Hx0 = 10.;2.0;10.0;  mT; AC
  Hz0 = 2.0;10.0;  mT; DC
  Hz_offset = 0.;10.; mT
  Hz0_array = FINDGEN(21)/20.*Hz0 + Hz_offset
  Hz0_array[0] += 0.1
  delta_t = 0.1 ; us
  ;  t_range = 1000.; us
  ;  N_tp = t_range/delta_t;
  frequency = 5.;20.; kHz
  N_tp = 1./frequency*1000./delta_t*2.5; 2.5 periods
  t_array = FINDGEN(N_tp)*delta_t
  signal_z_array = DBLARR(N_tp-1, N_ELEMENTS(Hz0_array))*0.
  Hx_array = Hx0 * SIN(t_array*frequency/1000.*2*!PI)
  iplot, t_array/1000., Hx_array, xtitle='time [ms]', ytitle='Hx field'

  FOR i=0, N_ELEMENTS(Hz0_array)-1 DO BEGIN
    H_total_array = SQRT(Hx_array^2 + Hz0_array[i]^2)
    ;iplot, t_array/1000., H_total_array, xtitle='time [ms]', ytitle='H total field'

    M_array = Msat_kAm*(1./TANH(beta_particle*H_total_array) - 1./(beta_particle*H_total_array))
    ;iplot, t_array/1000., M_array, xtitle='time [ms]', ytitle='M total'
    ;iplot, H_total_array, M_array, xtitle='H [mT]', ytitle='M total'

    M_nat_array = non_adiabatic(t_array, M_array, relaxation_time)
    IF i EQ N_ELEMENTS(Hz0_array)-1 THEN BEGIN
      window, 1
      plot, H_total_array, M_array, YRANGE=[MIN(M_array), MAX(M_array)], xtitle='H [mT]', ytitle='M total', title='M-H curve', BACKGROUND='FFFFFF'x, COLOR='000000'x
      oplot, H_total_array, M_nat_array, LINESTYLE=3, COLOR='0000FF'x
    ENDIF
    M_array = M_nat_array

    Mz_array = M_array * Hz0_array[i] / H_total_array
    IF i EQ N_ELEMENTS(Hz0_array)-1 THEN BEGIN
      window, 2
      ; Create the plot
      plot = PLOT(Hx_array, Mz_array, "r4D-", TITLE="Mz-Hx curve")
      ; Set some properties
      plot.SYM_INCREMENT = 10
      plot.SYM_COLOR = "blue"
      plot.SYM_FILLED = 1
      plot.SYM_FILL_COLOR = 0
      plot.xtitle='H_AC_x [mT]'
      plot.ytitle='Transverse Mz'
      ;      plot, Hx_array, Mz_array, YRANGE=[MIN(Mz_array), MAX(Mz_array)], xtitle='H_AC_x [mT]', ytitle='Transverse Mz', title='M-H curve', BACKGROUND='FFFFFF'x, COLOR='000000'x
      ;      oplot, H_total_array, M_nat_array, LINESTYLE=3, COLOR='0000FF'x
    ENDIF
    ;iplot, t_array/1000., Mz_array, xtitle='time [ms]', ytitle='M z'

    signal_z_array[*,i] = Mz_array[1:N_tp-1] - Mz_array[0:N_tp-2]
    ;iplot, t_array[0:N_tp-2]/1000., signal_z_array, xtitle='time [ms]', ytitle='signal z'

    ;plot, t_array[0:N_tp-2]/1000., signal_z_array, xtitle='time [ms]', ytitle='Signal z', BACKGROUND='FFFFFF'x, COLOR='000000'x
    ;oplot, t_array[0:N_tp-2]/1000., signal_z_array*0., LINESTYLE=3, COLOR='FF0000'x
    ;oplot, t_array[0:N_tp-2]/1000., Hx_array[0:N_tp-2]/MAX(Hx_array)*MAX(signal_z_array), COLOR='0000FF'x

  ENDFOR

  FOR i=0, N_ELEMENTS(Hz0_array)-1 DO BEGIN
    IF i EQ 0 THEN BEGIN
      window, 3
      plot, t_array[0:N_tp-2]/1000., REFORM(signal_z_array[*,i]), YRANGE=[MIN(signal_z_array), MAX(signal_z_array)], xtitle='time [ms]', ytitle='signal z', BACKGROUND='FFFFFF'x, COLOR='000000'x
      oplot, t_array[0:N_tp-2]/1000., REFORM(signal_z_array[*,i])*0., LINESTYLE=3, COLOR='FF0000'x
      oplot, t_array[0:N_tp-2]/1000., Hx_array[0:N_tp-2]/MAX(Hx_array)*MAX(signal_z_array), COLOR='FF0000'x
    ENDIF ELSE BEGIN
      oplot, t_array[0:N_tp-2]/1000., REFORM(signal_z_array[*,i]), COLOR='000033'x*i
    ENDELSE
    wait, 0.1
  ENDFOR

END

;@GJ, 2022/11/18, orthogonal excitation
;@GJ, 2022/11/19, adjusting the field strength
PRO orthogonal_excitation_diff_relax, frequency, particle_size, signal_z_array, Mz_array, relaxation_time_array, Hx0, Hz0, signal_x_array, Mx_array
  ;define the kernel
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  IF N_ELEMENTS(particle_size) EQ 0 THEN particle_size = 10.;24.4924;PBS, [25.2579], 1% Glycerol
  IF N_ELEMENTS(relaxation_time_array) EQ 0 THEN BEGIN
    relaxation_time_array = FINDGEN(101)/100. * 50.
    print, 'relaxation time array: ', relaxation_time_array
  ENDIF
  
  ;relaxation_time = 5.;0.5;2.;5.;0.;5.;2. ;us
  T_p = 25.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
  ;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  IF N_ELEMENTS(Hx0) EQ 0 THEN Hx0 = 10.;2.0;10.0;  mT; AC
  IF N_ELEMENTS(Hz0) EQ 0 THEN Hz0 = 2.0;10.0;  mT; DC
  Hz_offset = 0.;10.; mT
  ;Hz0_array = FINDGEN(21)/20.*Hz0 + Hz_offset
  ;Hz0_array[0] += 0.1
  Hz0_array = [Hz0]
  delta_t = 1.0;0.1 ; us
;  t_range = 1000.; us
;  N_tp = t_range/delta_t;
  IF N_ELEMENTS(frequency) EQ 0 THEN frequency = 2.5;5.;20.; kHz
  N_tp = 1./frequency*1000./delta_t*1.0;2.5; 2.5 periods
  t_array = FINDGEN(N_tp)*delta_t
  signal_z_array = DBLARR(N_tp-1, N_ELEMENTS(relaxation_time_array))*0.
  Mz_array = DBLARR(N_tp, N_ELEMENTS(relaxation_time_array))*0.
  signal_x_array = DBLARR(N_tp-1, N_ELEMENTS(relaxation_time_array))*0.
  Mx_array = DBLARR(N_tp, N_ELEMENTS(relaxation_time_array))*0.
;  signal_z_array = DBLARR(N_tp-1, N_ELEMENTS(Hz0_array))*0.
  Hx_array = Hx0 * COS(t_array*frequency/1000.*2*!PI)
  
  ;iplot, t_array/1000., Hx_array, xtitle='time [ms]', ytitle='Hx field'
  
  i=0
  FOR j=0, N_ELEMENTS(relaxation_time_array)-1 DO BEGIN
;  FOR i=0, N_ELEMENTS(Hz0_array)-1 DO BEGIN
    H_total_array = SQRT(Hx_array^2 + Hz0_array[i]^2)
    ;iplot, t_array/1000., H_total_array, xtitle='time [ms]', ytitle='H total field'

    M_array = Msat_kAm*(1./TANH(beta_particle*H_total_array) - 1./(beta_particle*H_total_array))
    ;iplot, t_array/1000., M_array, xtitle='time [ms]', ytitle='M total'
    ;iplot, H_total_array, M_array, xtitle='H [mT]', ytitle='M total'
    
    IF relaxation_time_array[j] LT 0.001 THEN BEGIN
      M_nat_array = non_adiabatic(t_array, M_array, relaxation_time_array[j])
    ENDIF ELSE BEGIN
      M_nat_array = M_array
    ENDELSE
    
;    IF j EQ N_ELEMENTS(relaxation_time_array)-1 THEN BEGIN
;;    IF i EQ N_ELEMENTS(Hz0_array)-1 THEN BEGIN
;     ; window, 1
;     ; plot, H_total_array, M_array, YRANGE=[MIN(M_array), MAX(M_array)], xtitle='H [mT]', ytitle='M total', title='M-H curve', BACKGROUND='FFFFFF'x, COLOR='000000'x
;    ;  oplot, H_total_array, M_nat_array, LINESTYLE=2, color='000000'x, PSYM=5
;    ENDIF

    Mz_array[*,j] = M_nat_array * Hz0_array[i] / H_total_array
    Mz_array[*,j] -= MIN(Mz_array[*,j])
    Mx_array[*,j] = M_nat_array * Hx_array / H_total_array
;    IF j EQ N_ELEMENTS(relaxation_time_array)-1 THEN BEGIN
;;    IF i EQ N_ELEMENTS(Hz0_array)-1 THEN BEGIN
;     
;      ; Create the plot
;      plot = PLOT(Hx_array, Mz_array, "r4D-", TITLE="Mz-Hx curve")
;      ; Set some properties
;      plot.SYM_INCREMENT = 10
;      plot.SYM_COLOR = "blue"
;      plot.SYM_FILLED = 1
;      plot.SYM_FILL_COLOR = 0
;      plot.xtitle='H_AC_x [mT]'
;      plot.ytitle='Transverse Mz'
;      window, 2
;      plot, Hx_array, Mz_array, YRANGE=[MIN(Mz_array), MAX(Mz_array)*1.5], xtitle='H_AC_x [mT]', ytitle='Transverse Mz', title='Mz-Hx curve', BACKGROUND='FFFFFF'x, COLOR='000000'x
;      oplot, H_total_array, M_nat_array, LINESTYLE=3, COLOR='0000FF'x
;    ENDIF ELSE BEGIN
;     ; oplot, Hx_array, Mz_array, COLOR='000033'x*j
;    ENDELSE
;    ;iplot, t_array/1000., Mz_array, xtitle='time [ms]', ytitle='M z'

    ;signal_z_array[*,i] = Mz_array[1:N_tp-1] - Mz_array[0:N_tp-2]
    signal_z_array[*,j] = (Mz_array[1:N_tp-1,j] - Mz_array[0:N_tp-2,j])/delta_t
    signal_x_array[*,j] = (Mx_array[1:N_tp-1,j] - Mx_array[0:N_tp-2,j])/delta_t
    ;iplot, t_array[0:N_tp-2]/1000., signal_z_array, xtitle='time [ms]', ytitle='signal z'

    ;plot, t_array[0:N_tp-2]/1000., signal_z_array, xtitle='time [ms]', ytitle='Signal z', BACKGROUND='FFFFFF'x, COLOR='000000'x
    ;oplot, t_array[0:N_tp-2]/1000., signal_z_array*0., LINESTYLE=3, COLOR='FF0000'x
    ;oplot, t_array[0:N_tp-2]/1000., Hx_array[0:N_tp-2]/MAX(Hx_array)*MAX(signal_z_array), COLOR='0000FF'x

  ENDFOR
  
;  FOR j=0, N_ELEMENTS(relaxation_time_array)-1 DO BEGIN
;;  FOR i=0, N_ELEMENTS(Hz0_array)-1 DO BEGIN
;    IF j EQ 0 THEN BEGIN
;      window, 3
;      plot, t_array[0:N_tp-2]/1000., REFORM(signal_z_array[*,j]), YRANGE=[MIN(signal_z_array), MAX(signal_z_array)], xtitle='time [ms]', ytitle='signal z', BACKGROUND='FFFFFF'x, COLOR='000000'x
;      oplot, t_array[0:N_tp-2]/1000., REFORM(signal_z_array[*,j])*0., LINESTYLE=3, COLOR='FF0000'x
;      oplot, t_array[0:N_tp-2]/1000., Hx_array[0:N_tp-2]/MAX(Hx_array)*MAX(signal_z_array), COLOR='FF0000'x
;    ENDIF ELSE BEGIN
;      oplot, t_array[0:N_tp-2]/1000., REFORM(signal_z_array[*,j]), COLOR='000033'x*j
;    ENDELSE
;    wait, 0.1
;  ENDFOR
  
END




;@GJ, 2022/8/20, spectrum imaging
PRO relaxation_spectrum_imaging_3mT
  ;open the baseline 25 deg of Synomag-D
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Baseline sav file", FILTER='*.sav')
  RESTORE, sav_name
  signal_matrix_0 = signal_matrix
  signal_matrix_sp_0 = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 50 degree of Synomag-D
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T50deg\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Hyperthermia sav file", FILTER='*.sav')
  RESTORE, sav_name
  signal_matrix_50d = signal_matrix
  signal_matrix_sp_50d = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 70% glycerol of Synomag-D
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly70\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Plaque sav file", FILTER='*.sav')
  RESTORE, sav_name
  signal_matrix_g70 = signal_matrix
  signal_matrix_sp_g70 = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 30% Gelatin of Synomag-D
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Catheter sav file", FILTER='*.sav')
  RESTORE, sav_name
  signal_matrix_ge30 = signal_matrix
  signal_matrix_sp_ge30 = signal_matrix_sp

  ;@GJ, 2022/7/27, save the image
  time_rise_flat = REFORM(signal_matrix_0[0, *])
  N_fields = N_ELEMENTS(signal_matrix_0[*,0])-1
  N_tp = N_ELEMENTS(signal_matrix_0[0,*])
  N_steps = 100

  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image_hr = CONGRID(REFORM(phantom_image[0,*,*]), N_steps, N_steps)
  ;  iimage, image_hr, title='original image (high resolution)'
  cur_dir = FILE_DIRNAME(sav_name)
  filename_ori_image = cur_dir+'\ori_image.png'
  WRITE_PNG, filename_ori_image, image_hr
  filename_System = cur_dir+'\original_image.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit1, image_hr[i, j]
  FREE_LUN, unit1

  filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
  image_bmp = READ_BMP(filename_bmp)
  d_phantom = REFORM(image_bmp[0,*,*])
  ;  iimage, d_phantom
  size_dp = SIZE(d_phantom, /dim)
  d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], N_steps, N_steps, /center, /minus_one)
  d_phantom_sq[*,0] = d_phantom_sq[*,1]
  d_phantom_sq[*,N_steps-1] = d_phantom_sq[*,N_steps-2]
  d_phantom_sq[N_steps-1,*] = d_phantom_sq[N_steps-2,*]
  ;extract catheter part
  catheter_index = WHERE(d_phantom_sq GT 240, catheter_count)
  catheter_index_2d = ARRAY_INDICES(d_phantom_sq, catheter_index)
  d_phantom_sq[0:40,*] = 0
  catheter_image = d_phantom_sq*0.
  catheter_image[catheter_index] = 256
  ;  iimage, image_hr+catheter_image, title='baseline + catheter'
  ;extract plaque part
  plaque_index = WHERE(d_phantom_sq GT 200, plauqe_count)
  plaque_index_2d = ARRAY_INDICES(d_phantom_sq, plaque_index)
  d_phantom_sq[*,26:*] = 0
  plaque_image = d_phantom_sq*0.
  plaque_image[plaque_index] = 256
  ; iimage, image_hr+plaque_image, title='baseline + plaque'
  ;extract thermia part
  thermia_index = WHERE(d_phantom_sq GT 0, thermia_count)
  thermia_index_2d = ARRAY_INDICES(d_phantom_sq, thermia_index)
  thermia_image = d_phantom_sq*0.
  thermia_image[thermia_index] = 256
  ; iimage, image_hr+thermia_image, title='baseline + thermia'
  total_image = image_hr+catheter_image+plaque_image+thermia_image
  
  field_array = (FINDGEN(N_steps)+1.)/N_steps*10.
  signal_matrix_200ave = REFORM(MEAN(signal_matrix_sP_0[*, 0:199, *], DIMENSION=2));i

  ;50deg
  signal_matrix_200ave_50d = REFORM(MEAN(signal_matrix_sP_50d[*, 0:199, *], DIMENSION=2));i

  ;70% glycerol
  signal_matrix_200ave_g70 = REFORM(MEAN(signal_matrix_sP_g70[*, 0:199, *], DIMENSION=2));i

  ;30% gelatin
  signal_matrix_200ave_ge30 = REFORM(MEAN(signal_matrix_sP_ge30[*, 0:199, *], DIMENSION=2));i

  ;field index to scan the MNPs
  field_excitation = 3.0;6.0;mT
  field_index = field_excitation/0.5-1 ; 6 mT
  
  ;  signals from each pixel
  signal_2d_200ave = DBLARR(N_steps, N_steps, N_tp)*0.
  signal_2d_200ave_baseline = DBLARR(N_steps, N_steps, N_tp)*0.
  FOR i=0, N_steps-1 DO BEGIN
    FOR j=0, N_steps-1 DO BEGIN
      flag = 0
      FOR k=0, N_ELEMENTS(catheter_index)-1 DO BEGIN
        IF i EQ catheter_index_2d[0, k] AND j EQ catheter_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 1;0;1; for catheter
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(plaque_index)-1 DO BEGIN
        IF i EQ plaque_index_2d[0, k] AND j EQ plaque_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 2; for plaque
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(thermia_index)-1 DO BEGIN
        IF i EQ thermia_index_2d[0, k] AND j EQ thermia_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 3; for thermia
          break
        ENDIF
      ENDFOR
      IF flag EQ 1 THEN BEGIN
        ;catheter
        signal_2d_200ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_200ave_ge30[field_index, *])
      ENDIF
      IF flag EQ 2 THEN BEGIN
        ;plaque
        signal_2d_200ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_200ave_g70[field_index, *])
      ENDIF
      IF flag EQ 3 THEN BEGIN
        ;thermia
        signal_2d_200ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_200ave_50d[field_index, *])
      ENDIF
      IF flag EQ 0 THEN BEGIN
        signal_2d_200ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_200ave[field_index, *])
      ENDIF
      ;;the baseline signal
      signal_2d_200ave_baseline[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_200ave[field_index, *])
    ENDFOR
  ENDFOR
  
  ;save the result
  ;save the decay signal
  flat_ratio = 0.9
  flat_index = N_tp*(1-flat_ratio)
  cur_dir = FILE_DIRNAME(sav_name)
  signal_preConv_decay_2d_200ave = signal_2d_200ave[*,*,flat_index+1:*]
  filename_System = cur_dir+'\image_preConv_signal3mT.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(signal_preConv_decay_2d_200ave[0,0,*])-1 DO WRITEU, unit1, signal_preConv_decay_2d_200ave[i, j, k]
  FREE_LUN, unit1

  ;save the baseline decay signal
  signal_preConv_decay_2d_200ave_baseline = signal_2d_200ave_baseline[*,*,flat_index+1:*]
  filename_System = cur_dir+'\image_preConv_signal_baseline3mT.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(signal_preConv_decay_2d_200ave_baseline[0,0,*])-1 DO WRITEU, unit1, signal_preConv_decay_2d_200ave_baseline[i, j, k]
  FREE_LUN, unit1
   
  ;define the kernel
  kernel_FFP = DBLARR(N_steps-1, N_steps-1)*0.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  particle_size = 24.4924;PBS, [25.2579], 1% Glycerol
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  Gx=7.0;mT/mm
  Gy=7.0;mT/mm
  FOV=50.;mm
  delta_xy = FOV/N_steps
  x_array = DBLARR(N_steps-1)
  FOR i=0, N_steps-2 DO BEGIN
    center_ind=N_steps/2-1
    x_cood = (i-center_ind)*delta_xy
    x_array[i] = x_cood
    FOR j=0, N_steps-2 DO BEGIN
      y_cood = (j-center_ind)*delta_xy
      ;define field ratio
      field_total = SQRT(Gx^2*x_cood^2 + Gy^2*y_cood^2 + field_excitation^2)
      field_ratio = field_excitation/field_total
;      print, 'field ratio = ', field_ratio
      
      ;define M ratio
      M_excitation = Msat_kAm*(1./TANH(beta_particle*field_excitation) - 1./(beta_particle*field_excitation))
      M_total = Msat_kAm*(1./TANH(beta_particle*field_total) - 1./(beta_particle*field_total))
      M_ratio = M_excitation/M_total
;      print, 'M ratio = ', M_ratio
      
      ;calculate the kernel
      kernel_FFP[i, j] = M_ratio * field_ratio
    ENDFOR
  ENDFOR
  iimage, kernel_FFP
  iplot, x_array, REFORM(kernel_FFP[*,N_steps/2-1]), xtitle='distance from FFP center [mm]', ytitle='ws', title='Weighting Factor'
  
  FFP_signal_2d_200ave = DBLARR(N_steps, N_steps, N_tp)*0.
  FFP_signal_2d_200ave_baseline = DBLARR(N_steps, N_steps, N_tp)*0.
  FOR k=0, N_tp-1 DO BEGIN
    print, 'k =', k, ' out of', N_tp
    
    image_temp = CONVOL(REFORM(signal_2d_200ave[*,*,k]), kernel_FFP, /NORMALIZE, /CENTER, /EDGE_ZERO)
    FFP_signal_2d_200ave[*,*,k] = REFORM(image_temp)
    
    image_temp = CONVOL(REFORM(signal_2d_200ave_baseline[*,*,k]), kernel_FFP, /NORMALIZE, /CENTER, /EDGE_ZERO)
    FFP_signal_2d_200ave_baseline[*,*,k] = REFORM(image_temp)
  ENDFOR
  
  ;save the decay signal
  FFP_signal_decay_2d_200ave = FFP_signal_2d_200ave[*,*,flat_index+1:*]
  filename_System = cur_dir+'\image_FFP_signal3mT.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(FFP_signal_decay_2d_200ave[0,0,*])-1 DO WRITEU, unit1, FFP_signal_decay_2d_200ave[i, j, k]
  FREE_LUN, unit1
  
  ;save the baseline decay signal
  FFP_signal_decay_2d_200ave_baseline = FFP_signal_2d_200ave_baseline[*,*,flat_index+1:*]
  filename_System = cur_dir+'\image_FFP_signal_baseline3mT.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(FFP_signal_decay_2d_200ave_baseline[0,0,*])-1 DO WRITEU, unit1, FFP_signal_decay_2d_200ave_baseline[i, j, k]
  FREE_LUN, unit1
  
  filename_sav = cur_dir+'\image_FFP_decay_signals3mT.sav'
  save, FFP_signal_decay_2d_200ave, FFP_signal_decay_2d_200ave_baseline, image_hr, total_image, filename=filename_sav
  ;fn='C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\image_FFP_decay_signals.sav'
  ;restore, fn
  
  FFP_signal_debye_2d_200ave = FFP_signal_2d_200ave[*,*,0:flat_index]
  FFP_signal_debye_2d_200ave_baseline = FFP_signal_2d_200ave_baseline[*,*,0:flat_index]
  filename_sav = cur_dir+'\image_FFP_debye_signals3mT.sav'
  save, FFP_signal_debye_2d_200ave, FFP_signal_debye_2d_200ave_baseline, image_hr, total_image, filename=filename_sav

END


;@GJ, 2022/8/20, spectrum imaging
PRO relaxation_spectrum_imaging, field_excitation, N_periods, title_temp
  ;open the baseline 25 deg of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Baseline sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\20220817121631signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_0 = signal_matrix
  signal_matrix_sp_0 = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 50 degree of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T50deg\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Hyperthermia sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T50deg\20220813100000signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_50d = signal_matrix
  signal_matrix_sp_50d = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 70% glycerol of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly70\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Plaque sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly70\20220813121126signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_g70 = signal_matrix
  signal_matrix_sp_g70 = signal_matrix_sp
  
  ;@GJ, 2022/9/16, adding the baseline as g1
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly1\20220813120203signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_g1 = signal_matrix
  signal_matrix_sp_g1 = signal_matrix_sp
  
  ;@GJ, 2022/9/16, adding the baseline as g1
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly40\20220916162759signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_g40 = signal_matrix
  signal_matrix_sp_g40 = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 30% Gelatin of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Catheter sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\20220817100728signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_ge30 = signal_matrix
  signal_matrix_sp_ge30 = signal_matrix_sp

  ;@GJ, 2022/7/27, save the image
  time_rise_flat = REFORM(signal_matrix_0[0, *])
  N_fields = N_ELEMENTS(signal_matrix_0[*,0])-1
  N_tp = N_ELEMENTS(signal_matrix_0[0,*])
  N_steps = 100
  
  ;@GJ, 2023/5/9, concentration map, could be the same
;  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw2.jpg'
  READ_JPEG, p_fn, phantom_image
  ;@GJ, 2022/9/15, removing the tumor region
;  phantom_image[*, 320:*, 0:370] = 0.
  image_hr = CONGRID(REFORM(phantom_image[0,*,*]), N_steps, N_steps)
  
;  ;@GJ, 2023/5/9, adding
;  p_fn = 'C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20230509b\vesselMIP_ori.png'
;  READ_PNG, p_fn, phantom_image
;  image_hr = CONGRID(phantom_image, N_steps, N_steps)
  
  ;  iimage, image_hr, title='original image (high resolution)'
  cur_dir = FILE_DIRNAME(p_fn)
  filename_ori_image = cur_dir+'\ori_image.png'
  WRITE_PNG, filename_ori_image, image_hr
  filename_System = cur_dir+'\original_image.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit1, image_hr[i, j]
  FREE_LUN, unit1

;  filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
  filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw2.bmp'
  image_bmp = READ_BMP(filename_bmp)
  d_phantom = REFORM(image_bmp[0,*,*])
  ;  iimage, d_phantom
  size_dp = SIZE(d_phantom, /dim)
  d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], N_steps, N_steps, /center, /minus_one)
  d_phantom_sq[*,0] = d_phantom_sq[*,1]
  d_phantom_sq[*,N_steps-1] = d_phantom_sq[*,N_steps-2]
  d_phantom_sq[N_steps-1,*] = d_phantom_sq[N_steps-2,*]
  ;extract catheter part
  catheter_index = WHERE(d_phantom_sq GT 240, catheter_count)
  catheter_index_2d = ARRAY_INDICES(d_phantom_sq, catheter_index)
  d_phantom_sq[0:40,*] = 0
  catheter_image = d_phantom_sq*0.
  catheter_image[catheter_index] = 256
  ;  iimage, image_hr+catheter_image, title='baseline + catheter'
  ;extract plaque part
  plaque_index = WHERE(d_phantom_sq GT 200, plauqe_count)
  plaque_index_2d = ARRAY_INDICES(d_phantom_sq, plaque_index)
  d_phantom_sq[*,26:*] = 0
  plaque_image = d_phantom_sq*0.
  plaque_image[plaque_index] = 256
  ; iimage, image_hr+plaque_image, title='baseline + plaque'
  ;extract thermia part
  thermia_index = WHERE(d_phantom_sq GT 0, thermia_count)
  thermia_index_2d = ARRAY_INDICES(d_phantom_sq, thermia_index)
  thermia_image = d_phantom_sq*0.
  thermia_image[thermia_index] = 256
  ; iimage, image_hr+thermia_image, title='baseline + thermia'
  total_image = image_hr+catheter_image+plaque_image+thermia_image

  IF N_ELEMENTS(N_periods) EQ 0 THEN BEGIN
    N_periods = 200
  ENDIF ELSE BEGIN
    IF N_periods LT 1 OR N_periods GT 400 THEN BEGIN
      N_periods = 200
    ENDIF
  ENDELSE
  print, '# of periods = ', N_periods
  
  signal_matrix_ave = REFORM(MEAN(signal_matrix_sP_g1[*, 0:N_periods-1, *], DIMENSION=2));i

  ;50deg
  signal_matrix_ave_50d = REFORM(MEAN(signal_matrix_sP_50d[*, 0:N_periods-1, *], DIMENSION=2));i

  ;70% glycerol
  signal_matrix_ave_g70 = REFORM(MEAN(signal_matrix_sP_g70[*, 0:N_periods-1, *], DIMENSION=2));i
  
  ;70% glycerol
  signal_matrix_ave_g40 = REFORM(MEAN(signal_matrix_sP_g40[*, 0:N_periods-1, *], DIMENSION=2));i
  
  ;30% gelatin
  signal_matrix_ave_ge30 = REFORM(MEAN(signal_matrix_sP_ge30[*, 0:N_periods-1, *], DIMENSION=2));i

  ;field index to scan the MNPs
  IF N_ELEMENTS(field_excitation) EQ 0 THEN BEGIN
    ;the default is 6 mT
    field_excitation = 6.0;mT
  ENDIF ELSE BEGIN
    IF field_excitation LT 0 OR field_excitation GT 10 THEN field_excitation = 6.0;mT
  ENDELSE
  field_index = field_excitation/0.5-1 ; 6 mT

  ;  signals from each pixel
  signal_2d_ave = DBLARR(N_steps, N_steps, N_tp)*0.
  signal_2d_ave_baseline = DBLARR(N_steps, N_steps, N_tp)*0.
  FOR i=0, N_steps-1 DO BEGIN
    FOR j=0, N_steps-1 DO BEGIN
      flag = 0
      FOR k=0, N_ELEMENTS(catheter_index)-1 DO BEGIN
        IF i EQ catheter_index_2d[0, k] AND j EQ catheter_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 1;0;1; for catheter
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(plaque_index)-1 DO BEGIN
        IF i EQ plaque_index_2d[0, k] AND j EQ plaque_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 2; for plaque
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(thermia_index)-1 DO BEGIN
        IF i EQ thermia_index_2d[0, k] AND j EQ thermia_index_2d[1, k] THEN BEGIN
          ;found and break
          ;removing the thermia region
          flag = 0;3;0;3; for thermia
          break
        ENDIF
      ENDFOR
      IF flag EQ 1 THEN BEGIN
        ;catheter
        signal_2d_ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_ave_ge30[field_index, *])
      ENDIF
      IF flag EQ 2 THEN BEGIN
        ;plaque
        signal_2d_ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_ave_g40[field_index, *])
      ENDIF
      IF flag EQ 3 THEN BEGIN
        ;thermia
        signal_2d_ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_ave_50d[field_index, *])
      ENDIF
      IF flag EQ 0 THEN BEGIN
        signal_2d_ave[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_ave[field_index, *])
      ENDIF
      ;;the baseline signal
      signal_2d_ave_baseline[i, j, *] = image_hr[i, j] * REFORM(signal_matrix_ave[field_index, *])
    ENDFOR
  ENDFOR

  ;save the result
  ;save the decay signal
  flat_ratio = 0.9
  flat_index = N_tp*(1-flat_ratio)
;  cur_dir = FILE_DIRNAME(sav_name)
;  signal_preConv_decay_2d_ave = signal_2d_ave[*,*,flat_index+1:*]
;  filename_System = cur_dir+'\image_preConv_signal.dat'
;  openw,unit1,filename_System,/get_lun
;  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(signal_preConv_decay_2d_ave[0,0,*])-1 DO WRITEU, unit1, signal_preConv_decay_2d_ave[i, j, k]
;  FREE_LUN, unit1

  ;save the baseline decay signal
;  signal_preConv_decay_2d_ave_baseline = signal_2d_ave_baseline[*,*,flat_index+1:*]
;  filename_System = cur_dir+'\image_preConv_signal_baseline.dat'
;  openw,unit1,filename_System,/get_lun
;  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(signal_preConv_decay_2d_ave_baseline[0,0,*])-1 DO WRITEU, unit1, signal_preConv_decay_2d_ave_baseline[i, j, k]
;  FREE_LUN, unit1

  ;define the kernel
  kernel_FFP = DBLARR(N_steps-1, N_steps-1)*0.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  particle_size = 24.4924;PBS, [25.2579], 1% Glycerol
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  Gx=7.0;mT/mm
  Gy=7.0;mT/mm
  FOV=50.;mm
  delta_xy = FOV/N_steps
  x_array = DBLARR(N_steps-1)
  FOR i=0, N_steps-2 DO BEGIN
    center_ind=N_steps/2-1
    x_cood = (i-center_ind)*delta_xy
    x_array[i] = x_cood
    FOR j=0, N_steps-2 DO BEGIN
      y_cood = (j-center_ind)*delta_xy
      ;define field ratio
      field_total = SQRT(Gx^2*x_cood^2 + Gy^2*y_cood^2 + field_excitation^2)
      field_ratio = field_excitation/field_total
      ;      print, 'field ratio = ', field_ratio

      ;define M ratio
      M_excitation = Msat_kAm*(1./TANH(beta_particle*field_excitation) - 1./(beta_particle*field_excitation))
      M_total = Msat_kAm*(1./TANH(beta_particle*field_total) - 1./(beta_particle*field_total))
      M_ratio = M_excitation/M_total
      ;      print, 'M ratio = ', M_ratio

      ;calculate the kernel
      kernel_FFP[i, j] = M_ratio * field_ratio
    ENDFOR
  ENDFOR
;  iimage, kernel_FFP
;  iplot, x_array, REFORM(kernel_FFP[*,N_steps/2-1]), xtitle='distance from FFP center [mm]', ytitle='ws', title='Weighting Factor'
  x=x_array
  y=x
  data = kernel_FFP;/MAX(kernel_FFP)
  ; Load program colors.
  cgLoadCT, 25;, /Brewer, /Reverse
  ; Create a position for the surface in the display window.
  pos = [0.1, 0.1, 0.9, 0.96]
  ; With surfaces, the charsize often needs to be adjusted in PostScript.
  thisCharSize = (!D.Name EQ 'PS') ? cgDefCharsize()*1.75 : cgDefCharsize()
  ; Draw the fluence surface plot as a shaded surface.
  cgSurf, data, x, y, RotX=50, RotZ=25, Position=pos, /Elevation, ZStyle=4, $
    /Shaded, ZTickformat='(A1)', XTitle='y Location (mm)', YTitle='z Location (mm)', $
    Charsize=thisCharSize
  ; Overlay the surface lines and add a skirt.
;  cgSurf, data, x, y, RotX=50, RotZ=25, Position=pos, Skirt=Min(data), $
;    XTickformat='(A1)', YTickformat='(A1)', ZStyle=4, /NoErase, Color='gray', $
;    Charsize=thisCharSize
;  ; Add a color bar to the plot.
  cgColorbar, Range=[0,1], Title='Signal Weighting Factor', TLocation='top', $
    Position=[0.86, 0.15, 0.91, 0.85], Charsize=thisCharSize*0.8 
 
  ;GJ, 2022/9/15, calculate the average
  one_fifth_index = WHERE(data GT 0.1, Ncount)
  ;GJ, 2022/9/18, to avoid the kernel is too small
  IF Ncount LT 16 THEN Ncount = 16.
  ;GJ, 2022/9/18, fix the bug to put 1 in the center of the kernel
  kernel_FFP = data[N_steps/2-FIX(sqrt(Ncount))/2-1:N_steps/2+FIX(sqrt(Ncount))/2-1, N_steps/2-FIX(sqrt(Ncount))/2-1:N_steps/2+FIX(sqrt(Ncount))/2-1]
;  kernel_FFP = kernel_FFP*0.
;  kernel_FFP[FIX(sqrt(Ncount))/2, FIX(sqrt(Ncount))/2] = 1.
;  iimage, kernel_FFP, title='current kernel shape'
  
  FFP_signal_2d_ave = DBLARR(N_steps, N_steps, N_tp)*0.
  FFP_signal_2d_ave_baseline = DBLARR(N_steps, N_steps, N_tp)*0.
  FOR k=0, N_tp-1 DO BEGIN
;    print, 'FFP convolution timepoint #: k =', k, ' out of', N_tp
;    image_temp = CONVOL(REFORM(signal_2d_ave[*,*,k]), kernel_FFP, /NORMALIZE, /CENTER, /EDGE_ZERO)
    image_temp = CONVOL_FFT(REFORM(signal_2d_ave[*,*,k]), kernel_FFP)
    FFP_signal_2d_ave[*,*,k] = REFORM(image_temp)

;    image_temp = CONVOL(REFORM(signal_2d_ave_baseline[*,*,k]), kernel_FFP, /NORMALIZE, /CENTER, /EDGE_ZERO)
;    FFP_signal_2d_ave_baseline[*,*,k] = REFORM(image_temp)
  ENDFOR

  ;save the decay signal
  FFP_signal_decay_2d_ave = FFP_signal_2d_ave[*,*,flat_index+1:*]
  filename_System = cur_dir+'\'+title_temp+'_FFP_signal_decay.dat'
  filename_dir = FILE_DIRNAME(filename_System)
  IF FILE_TEST(filename_dir) NE 1 THEN FILE_MKDIR, filename_dir
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(FFP_signal_decay_2d_ave[0,0,*])-1 DO WRITEU, unit1, FFP_signal_decay_2d_ave[i, j, k]
  FREE_LUN, unit1

  ;save the baseline decay signal
  
;  filename_System = cur_dir+'\image_FFP_signal_baseline.dat'
;  openw,unit1,filename_System,/get_lun
;  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO FOR k=0, N_ELEMENTS(FFP_signal_decay_2d_ave_baseline[0,0,*])-1 DO WRITEU, unit1, FFP_signal_decay_2d_ave_baseline[i, j, k]
;  FREE_LUN, unit1
  
  FFP_signal_decay_2d_ave_baseline = FFP_signal_2d_ave_baseline[*,*,flat_index+1:*]
;  filename_sav = cur_dir+'\image_FFP_decay_signals.sav'
;  save, field_excitation, FFP_signal_2d_ave, FFP_signal_decay_2d_ave, FFP_signal_decay_2d_ave_baseline, image_hr, total_image, filename=filename_sav
  ;fn='C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\image_FFP_decay_signals.sav'
  ;restore, fn

  FFP_signal_debye_2d_ave = FFP_signal_2d_ave[*,*,0:flat_index]
  FFP_signal_debye_2d_ave_baseline = FFP_signal_2d_ave_baseline[*,*,0:flat_index]
  filename_sav = cur_dir+'\'+title_temp+'_FFP_signals.sav'
  save, field_excitation, N_periods, flat_index, d_phantom_sq, N_tp, FFP_signal_2d_ave, FFP_signal_2d_ave_baseline, image_hr, total_image, catheter_index, plaque_index, thermia_index, filename=filename_sav
  
;  RETURN, result_array
END

;@GJ, 2022/8/31, save all FFP signals
PRO batch_FFP_spectrum_save_N_averages
  
  num_k = 41;3;41
  k_array = FINDGEN(num_k)*10./(num_k-1)*40.
  k_array[0]=1
  
  FOR k=0, N_ELEMENTS(k_array)-1 DO BEGIN
    field_excitation = 4.5;6.0;
    N_periods = k_array[k]
    IF field_excitation LT 10 THEN BEGIN
      title_temp = STRING(field_excitation, format='(f3.1)')+'mT\'+STRTRIM(FIX(N_periods), 1)+'period'
    ENDIF ELSE BEGIN
      title_temp = STRING(field_excitation, format='(f4.1)')+'mT\'+STRTRIM(FIX(N_periods), 1)+'period'
    ENDELSE
    relaxation_spectrum_imaging, field_excitation, N_periods, title_temp
  ENDFOR

END

;@GJ, 2022/8/31, save all FFP signals
PRO batch_FFP_spectrum_save_field_amplitudes

  num_k = 20.
  k_array = FINDGEN(num_k)/2. + 0.5

  FOR k=0, N_ELEMENTS(k_array)-1 DO BEGIN
    field_excitation = k_array[k];
    N_periods = 200
    IF field_excitation LT 10 THEN BEGIN
      title_temp = STRTRIM(FIX(N_periods), 1)+'period\'+STRING(field_excitation, format='(f3.1)')+'mT'
    ENDIF ELSE BEGIN
      title_temp = STRTRIM(FIX(N_periods), 1)+'period\'+STRING(field_excitation, format='(f4.1)')+'mT'
    ENDELSE
    relaxation_spectrum_imaging, field_excitation, N_periods, title_temp
  ENDFOR

END

;@GJ, 2022/8/30
PRO SNR_analysis_spectrum_imaging
  
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\6mT\'
  sav_file_dir = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, /DIRECTORY, TITLE='FFP signal save file directory')
  
  save_file_list = FILE_SEARCH(sav_file_dir, '*_FFP_signals.sav', count=nrfile)
  num_k = nrfile;11;41
  N_periods_array = DBLARR(num_k)*0.
  SNR_vessel_array = DBLARR(num_k)*0.
  CNR_catheter_array = DBLARR(num_k)*0.
  CNR_plaque_array = DBLARR(num_k)*0.
  CNR_thermia_array = DBLARR(num_k)*0.
  FOR k=0, nrfile-1 DO BEGIN
    RESTORE, save_file_list[k]
    N_periods_array[k] = N_periods
    ;SNR calculation
    ;get the magnetization or signal averages
    bg_image = image_hr * 0.
    bg_image_baseline = image_hr * 0.
    FOR i=0, N_ELEMENTS(image_hr[*,0])-1 DO BEGIN
      FOR j=0, N_ELEMENTS(image_hr[0,*])-1 DO BEGIN
        bg_image[i, j] = MEAN(FFP_signal_2d_ave[i, j, *])
        bg_image_baseline[i, j] = MEAN(FFP_signal_2d_ave_baseline[i, j, *])
      ENDFOR
    ENDFOR

    IF N_periods EQ 20 OR N_periods EQ 201 OR N_periods EQ 50 THEN iimage, bg_image, title='reconstructed image ('+STRTRIM(N_periods, 1)+') times average'

    signal_catheter = MEAN(bg_image[catheter_index])
    signal_plaque = MEAN(bg_image[plaque_index])
    signal_thermia = MEAN(bg_image[thermia_index])
    noise = STDDEV(bg_image[0:20, 65:*])
    temp_image_hr = image_hr * 0.
    temp_image_hr[20:45, 65:*] = image_hr[20:45, 65:*]
    vessel_index = WHERE(temp_image_hr GT 170, vessel_count)
    signal_vessel = MEAN(bg_image[vessel_index])
    signal_vessel_array = DBLARR(N_tp) * 0.
    FOR k_tp=0, N_tp-1 DO BEGIN
      temp_signal = REFORM(FFP_signal_2d_ave[*, *, k_tp])
      signal_vessel_array[k_tp] = MEAN(temp_signal[vessel_index])
    ENDFOR
;    iplot, signal_vessel_array, title='vessel signal'
    
    filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
    filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw2.bmp'
    image_bmp = READ_BMP(filename_bmp)
    d_phantom = REFORM(image_bmp[0,*,*])
    ;  iimage, d_phantom
    size_dp = SIZE(d_phantom, /dim)
    d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], 100, 100, /center, /minus_one)
    d_phantom_sq[thermia_index]=0
    d_phantom_sq[vessel_index]=0
    d_phantom_sq[plaque_index]=0
    d_phantom_sq[catheter_index]=0
    tumor_unheated_index = WHERE(d_phantom_sq EQ 38 OR d_phantom_sq EQ 152, tumor_unheated_count)
    signal_tumor_unheated = MEAN(bg_image[tumor_unheated_index])
    print, 'N_periods: ', N_periods
    SNR_vessel = signal_vessel/noise
    print, 'vessel SNR: ', SNR_vessel
    CNR_catheter = ABS(signal_vessel-signal_catheter)/noise
    print, 'catheter CNR: ', CNR_catheter
    CNR_plaque = ABS(signal_vessel-signal_plaque)/noise
    print, 'plaque CNR: ', CNR_plaque
    CNR_thermia = ABS(signal_tumor_unheated-signal_thermia)/noise
    print, 'hyperthermia CNR: ', CNR_thermia
    
    result_array = [N_periods, SNR_vessel, CNR_catheter, CNR_plaque, CNR_thermia]
    SNR_vessel_array[k] = result_array[1]
    CNR_catheter_array[k] = result_array[2]
    CNR_plaque_array[k] = result_array[3]
    CNR_thermia_array[k] = result_array[4]
    
  ENDFOR
  
  sort_ind = SORT(N_periods_array)
  iplot, N_periods_array[sort_ind], SNR_vessel_array[sort_ind], title='SNR vessel', xtitle='# of periods', ytitle='SNR'
  iplot, N_periods_array[sort_ind], CNR_catheter_array[sort_ind], title='CNR catheter', xtitle='# of periods', ytitle='CNR'
  iplot, N_periods_array[sort_ind], CNR_plaque_array[sort_ind], title='CNR plaque', xtitle='# of periods', ytitle='CNR'
  iplot, N_periods_array[sort_ind], CNR_thermia_array[sort_ind], title='CNR thermia', xtitle='# of periods', ytitle='CNR'

  plot, N_periods_array[sort_ind], SNR_vessel_array[sort_ind], YRANGE = [0, CEIL(MAX(SNR_vessel_array))], title='SNR vessel', xtitle='# of periods', ytitle='SNR'
  oplot, N_periods_array[sort_ind], CNR_catheter_array[sort_ind]
  oplot, N_periods_array[sort_ind], CNR_plaque_array[sort_ind]
  oplot, N_periods_array[sort_ind], CNR_thermia_array[sort_ind]
END

;@GJ, 2022/8/31, calculate SNR on relaxation time maps
;@GJ, 2022/9/20, modify the codes for signal evaluation
PRO SNR_analysis_debye_relaxation_maps

  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\6mT\'
  sav_file_dir = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, /DIRECTORY, TITLE='FFP signal save file directory')

  save_file_list = FILE_SEARCH(sav_file_dir, '*_FFP_signals.sav', count=nrfile)
  num_k = nrfile;11;41
  N_periods_array = DBLARR(num_k)*0.
  N_fieldAmp_array = DBLARR(num_k)*0.
  S_vessel_array = DBLARR(num_k)*0.
  S_catheter_array = DBLARR(num_k)*0.
  S_plaque_array = DBLARR(num_k)*0.
  S_thermia_array = DBLARR(num_k)*0.
  w_flag=1
  FOR k=0, nrfile-1 DO BEGIN
    RESTORE, save_file_list[k]
    N_periods_array[k] = N_periods
    N_fieldAmp_array[k] = field_excitation
    
    FFP_signal_debye_2d_ave = FFP_signal_2d_ave[*,*,0:flat_index]
    FFP_signal_debye_2d_ave_baseline = FFP_signal_2d_ave_baseline[*,*,0:flat_index]
    
    N_steps=N_ELEMENTS(FFP_signal_debye_2d_ave[*,0,0])
    N_tp=N_ELEMENTS(FFP_signal_debye_2d_ave[0,0,*])
    Debye_image=DBLARR(N_steps,N_steps)*0.
    bg_image=DBLARR(N_steps,N_steps)*0.

    ;GJ, 2022/2/11, define debby tao
    N_debye_tao = 600; 200
    ;GJ, 2022/1/30, calculate Debye model
    relaxation_time = FINDGEN(N_debye_tao)/10.+1.

    N_particle_size = 1
    particle_size_array = [24.4924]

    delta_t = 1.0
    time_rise = FINDGEN(N_tp+50.)*delta_t
    H_max = field_excitation
    H_array = [DBLARR(50)*0.-H_max, FINDGEN(N_tp)/(N_tp-1)*2*H_max-H_max]
    FOR i=0, N_steps-1 DO BEGIN
      FOR j=0, N_steps-1 DO BEGIN
        IF MEAN(FFP_signal_debye_2d_ave[i,j,*]) LT 0.001 THEN continue
        ;do the biexponential curve fit
        ;Provide an initial guess of the function’s parameters.
        M_meas_temp=TOTAL(REFORM(FFP_signal_debye_2d_ave[i,j,*]), /CUMULATIVE) * delta_t
        M_meas = [DBLARR(50)*0., M_meas_temp] -(MAX(M_meas_temp)-MIN(M_meas_temp))/2.
        bg_image[i,j] = TOTAL(REFORM(FFP_signal_debye_2d_ave[i,j,*])) * delta_t

        MSE_array = DBLARR(N_particle_size, N_debye_tao)*0.+10000.
        correlation_array = DBLARR(N_particle_size, N_debye_tao)*0.
        FOR kp=0, N_particle_size-1 DO BEGIN
          particle_size = particle_size_array[kp]
          Msat_T = 0.551 ; T/u0
          Msat_kAm = 0.551/4./!PI*10000.; kA/m
          T_p = 20.;36.5; degree
          T_p_kelvin = 273.15 + T_p ; in Kelvin
          beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;          beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
          FWHM_particle = 4.16/beta_particle/3.5; mT/mm
          M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))

          FOR l = 0, N_debye_tao-1 DO BEGIN
            M_nat = non_adiabatic(time_rise, M, relaxation_time[l])
;            M_nat -= MIN(M_nat)
            M_nat = M_nat/MEAN(M_nat)*MEAN(M_meas)
            correlation_array[kp, l] = CORRELATE(M_nat, M_meas)
            MSE_array[kp, l] = calMSE(M_nat, M_meas)
          ENDFOR
        ENDFOR

        min_MSE = MIN(MSE_array, min_ind)
        min_ind2d = ARRAY_INDICES(MSE_array, min_ind)

        tao_debye = relaxation_time[min_ind2d[1]]
;        print, '[',i,', ',j,']: (Debye) relaxation time = ', tao_debye

        Debye_image[i,j] = tao_debye
        
        ;
;        IF i EQ 62 AND j EQ 62 THEN BEGIN
        IF i EQ 30 AND j EQ 55 THEN BEGIN
          print, 'plaque, min ind', min_ind2d[1]
          iplot, relaxation_time, MSE_array, xtitle='relaxation time [us]', ytitle='MSE'
          M_nat_tao = non_adiabatic(time_rise, M, tao_debye)
          M_nat_tao = M_nat_tao/MEAN(M_nat_tao)*MEAN(M_meas)
          window, 0
          plot, time_rise, M_meas, xtitle='time [us]', ytitle='M'
          oplot, time_rise, M_nat
          oplot, time_rise, M_nat_tao, LINESTYLE=2
        ENDIF
      ENDFOR
    ENDFOR
    
    IF N_periods EQ 20 OR N_periods EQ 50 OR N_periods EQ 200 THEN BEGIN
      N_steps_new = 512
      Debye_image_new = CONGRID(Debye_image, N_steps_new, N_steps_new, /INTERP, /CENTER)
      bg_image_new = CONGRID(bg_image, N_steps_new, N_steps_new, /INTERP, /CENTER)
      image_hr_new = CONGRID(image_hr, N_steps_new, N_steps_new, /INTERP, /CENTER)
      iimage, Debye_image_new, title='Debye time map '+STRTRIM(FIX(N_periods), 1)
      window, w_flag
      cgDisplay, N_steps_new, N_steps_new, title='Debye time map '+STRTRIM(FIX(N_periods), 1), free=1
      cgImage, BYTSCL(bg_image_new), CTIndex=0
      mask_image = IMAGE_THRESHOLD(bg_image_new, THRESHOLD=tf, /MAXENTROPY)
      ;@GJ, 2022/9/18, mask out the non-vessel background
      nonZero_ind = WHERE(bg_image_new GT tf, N_nz,  COMPLEMENT=zero_ind, NCOMPLEMENT=N_z)
      min_Debye_new = MIN(Debye_image_new[nonZero_ind])
      max_Debye_new = MAX(Debye_image_new[nonZero_ind])
      print, 'Debye range: ', min_Debye_new, ', ', max_Debye_new
;      cgImage, HIST_EQUAL(Debye_image_new, MINV=min_Debye_new, MAXV=max_Debye_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      cgImage, BYTSCL(Debye_image_new, MIN=min_Debye_new, MAX=max_Debye_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
      IF N_z GT 1 THEN Debye_image_new[Zero_ind] = 0
;      cgImage, BYTSCL(Debye_image_new, MIN=min_Debye_new, MAX=max_Debye_new, TOP=255), CTIndex=33, Transparent=50;, Missing_Value=0;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
      cgImage, BYTSCL(Debye_image_new, MIN=0, MAX=40, TOP=255), CTIndex=33, Transparent=50;, Missing_Value=0;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
      w_flag++
      
      ;save the data Debye map
      filename_System = sav_file_dir+'\Debye_'+STRTRIM(FIX(N_periods), 1)+'.dat'
      filename_dir = FILE_DIRNAME(filename_System)
      IF FILE_TEST(filename_dir) NE 1 THEN FILE_MKDIR, filename_dir
      openw,unit1,filename_System,/get_lun
      FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit1, Debye_image[i,j]
      FREE_LUN, unit1
    ENDIF
    
    signal_catheter = MEAN(Debye_image[catheter_index])
    signal_plaque = MEAN(Debye_image[plaque_index])
    signal_thermia = MEAN(Debye_image[thermia_index])
    noise = STDDEV(Debye_image[0:20, 65:*])
    temp_image_hr = image_hr * 0.
    temp_image_hr[20:45, 65:*] = image_hr[20:45, 65:*]
    vessel_index = WHERE(temp_image_hr GT 170, vessel_count)
    signal_vessel = MEAN(Debye_image[vessel_index])
    
    filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
    filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw2.bmp'
    image_bmp = READ_BMP(filename_bmp)
    d_phantom = REFORM(image_bmp[0,*,*])
    ;  iimage, d_phantom
    size_dp = SIZE(d_phantom, /dim)
    d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], 100, 100, /center, /minus_one)
    d_phantom_sq[thermia_index]=0
    d_phantom_sq[vessel_index]=0
    d_phantom_sq[plaque_index]=0
    d_phantom_sq[catheter_index]=0
    tumor_unheated_index = WHERE(d_phantom_sq EQ 38 OR d_phantom_sq EQ 152, tumor_unheated_count)
    signal_tumor_unheated = MEAN(Debye_image[tumor_unheated_index])
    print, 'N_periods: ', N_periods, '; field amplitude: ', field_excitation
    print, 'vessel: ', signal_vessel
    print, 'catheter: ', signal_catheter
    print, 'plaque: ', signal_plaque
;    print, 'hyperthermia CNR: ', CNR_thermia

    S_vessel_array[k] = signal_vessel
    S_catheter_array[k] = signal_catheter
    S_plaque_array[k] = signal_plaque
;    S_thermia_array[k] = signal_thermia
  ENDFOR
  
  IF ABS(MAX(N_periods_array)-MIN(N_periods_array)) GT 1. THEN BEGIN
    sort_ind = SORT(N_periods_array)
    window, w_flag
    plot, N_periods_array[sort_ind], S_vessel_array[sort_ind], YRANGE = [0, CEIL(MAX(S_vessel_array))], title='S Debye', xtitle='# of periods', ytitle='relaxation time'
    oplot, N_periods_array[sort_ind], S_catheter_array[sort_ind], LINESTYLE=1
    oplot, N_periods_array[sort_ind], S_plaque_array[sort_ind], LINESTYLE=2
;    oplot, N_periods_array[sort_ind], S_thermia_array[sort_ind]
  ENDIF ELSE BEGIN
    sort_ind = SORT(N_fieldAmp_array)
    window, w_flag
    plot, N_fieldAmp_array[sort_ind], S_vessel_array[sort_ind], YRANGE = [0, CEIL(MAX(S_vessel_array))], title='S Debye', xtitle='Field amplitude', ytitle='relaxation time'
    oplot, N_fieldAmp_array[sort_ind], S_catheter_array[sort_ind], LINESTYLE=1
    oplot, N_fieldAmp_array[sort_ind], S_plaque_array[sort_ind], LINESTYLE=2
;    oplot, N_fieldAmp_array[sort_ind], S_thermia_array[sort_ind]
  ENDELSE
  
  current_time = string(SYSTIME(/JULIAN, /UTC), FORMAT='(f8.0)')
  output_filename = sav_file_dir+current_time+'xls'
  ;write file
  openw,unit,output_filename,/get_lun
  first_line =['Field Amplitude', 'N_periods', 'vessel', 'catheter', 'plaque', 'hyperthermia']
  printf, unit, FORMAT = '(10(A, %"\t"))', first_line
  FOR i=0L, nrfile-1L DO BEGIN
    printf, unit, FORMAT = '(10(A, %"\t"))', N_fieldAmp_array[sort_ind[i]], N_periods_array[sort_ind[i]], S_vessel_array[sort_ind[i]], S_catheter_array[sort_ind[i]], S_plaque_array[sort_ind[i]], S_thermia_array[sort_ind[i]]
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit

END


;@GJ, 2022/9/1, calculate SNR on NB relaxation time maps
PRO SNR_analysis_NB_relaxation_maps

  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\4.5mT\'
  sav_file_dir = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, /DIRECTORY, TITLE='FFP signal save file directory')

  save_file_list = FILE_SEARCH(sav_file_dir, '*_FFP_signals.sav', count=nrfile)
  num_k = nrfile;11;41
  N_periods_array = DBLARR(num_k)*0.
  N_fieldAmp_array = DBLARR(num_k)*0.
  S_vessel_array_Neel = DBLARR(num_k)*0.
  S_catheter_array_Neel = DBLARR(num_k)*0.
  S_plaque_array_Neel = DBLARR(num_k)*0.
  S_thermia_array_Neel = DBLARR(num_k)*0.
  
  S_vessel_array_Brownian = DBLARR(num_k)*0.
  S_catheter_array_Brownian = DBLARR(num_k)*0.
  S_plaque_array_Brownian = DBLARR(num_k)*0.
  S_thermia_array_Brownian = DBLARR(num_k)*0.
  w_flag=0
  FOR k=0, nrfile-1 DO BEGIN
    RESTORE, save_file_list[k]
    N_periods_array[k] = N_periods
    N_fieldAmp_array[k] = field_excitation
    
    FFP_signal_decay_2d_ave = FFP_signal_2d_ave[*,*,flat_index+1:*]
    FFP_signal_decay_2d_ave_baseline = FFP_signal_2d_ave_baseline[*,*,flat_index+1:*]

    N_steps=N_ELEMENTS(FFP_signal_decay_2d_ave[*,0,0])
    N_tp=N_ELEMENTS(FFP_signal_decay_2d_ave[0,0,*])
    delta_t = 1.0
    Neel_image=DBLARR(N_steps,N_steps)*0.
    Brownian_image=DBLARR(N_steps,N_steps)*0.
    bg_image=DBLARR(N_steps,N_steps)*0.
   
    FOR i=0, N_steps-1 DO BEGIN
      FOR j=0, N_steps-1 DO BEGIN
        IF MEAN(FFP_signal_decay_2d_ave[i,j,*]) LT 0.001 THEN continue
        ;do the biexponential curve fit
        ;Provide an initial guess of the function’s parameters.
        Y=TOTAL(REFORM(FFP_signal_decay_2d_ave[i,j,*]), /CUMULATIVE) * delta_t
        IF i EQ 62 AND j EQ 62 THEN BEGIN
          ;save the 001
          signal_array = REFORM(FFP_signal_decay_2d_ave[i,j,*])
          X = FINDGEN(N_tp) * delta_t
          N_int = 8000.
          noise = MEAN(signal_array[N_tp-1-24:N_tp-1])
          y_data_ori = ABS(signal_array - noise)*1000.
          X_data_int = FINDGEN(N_int+1)/N_int*MAX(X)
          y_data_int = INTERPOL(y_data_ori, X, X_data_int)
          ;GJ, 2022/5/10, should not include the t=0 signal
          write_UPEN_001_file, save_file_list[k], X_data_int[1:*], y_data_int[1:*], file_basename(save_file_list[k], '.sav')+'_plaque'
        ENDIF
        IF i EQ 21 AND j EQ 41 THEN BEGIN
          ;save the 001
          signal_array = REFORM(FFP_signal_decay_2d_ave[i,j,*])
          X = FINDGEN(N_tp) * delta_t
          N_int = 8000.
          noise = MEAN(signal_array[N_tp-1-24:N_tp-1])
          y_data_ori = ABS(signal_array - noise)*1000.
          X_data_int = FINDGEN(N_int+1)/N_int*MAX(X)
          y_data_int = INTERPOL(y_data_ori, X, X_data_int)
          ;GJ, 2022/5/10, should not include the t=0 signal
          write_UPEN_001_file, save_file_list[k], X_data_int[1:*], y_data_int[1:*], file_basename(save_file_list[k], '.sav')+'_catheter'
        ENDIF
        bg_image[i,j]=TOTAL(REFORM(FFP_signal_decay_2d_ave[i,j,*])) * delta_t
        weights=1./Y
        X=FINDGEN(N_tp)*delta_t
        B = [-MAX(Y)/2.0, -0.01, -MAX(Y)/2.0, -0.1, MIN(Y), 0.]; ;@GJ, 2022/7/10, B = [-5.0, -0.01, -5.0, -0.1, MAX(Y), 0.]
        FITA = [1,1,1,1,1,0]
        ;Compute the parameters.
        yfit_biex = CURVEFIT(X, Y, weights, B, SIGMA, FITA=FITA, FUNCTION_NAME='gfunct_biex')
        bi_tao_1 = 1./ABS(B[1])
        bi_tao_2 = 1./ABS(B[3])
        print, 'pixel = [',i,',',j,',]'
        print, 'Relaxation time_1 = ', bi_tao_1, ' [us]'
        print, 'Relaxation time_2 = ', bi_tao_2, ' [us]'
        Neel_image[i,j] = MIN([bi_tao_1, bi_tao_2])
        Brownian_image[i,j] = MAX([bi_tao_1, bi_tao_2])
      ENDFOR
    ENDFOR
    
    IF N_periods EQ 20 OR N_periods EQ 201 OR N_periods EQ 50 THEN BEGIN
      N_steps_new = 512
      Neel_image_new = CONGRID(Neel_image, N_steps_new, N_steps_new, /INTERP, /CENTER)
      Brownian_image_new = CONGRID(Brownian_image, N_steps_new, N_steps_new, /INTERP, /CENTER)
      bg_image_new = CONGRID(bg_image, N_steps_new, N_steps_new, /INTERP, /CENTER)
      iimage, Neel_image_new, title='Neel time map '+STRTRIM(FIX(N_periods), 1)
      iimage, Brownian_image_new, title='Brownian time map '+STRTRIM(FIX(N_periods), 1)
      nonZero_ind = WHERE(bg_image_new GT 2400, N_nz)
      window, 1+w_flag*2
      cgDisplay, N_steps_new, N_steps_new, title='Neel time map '+STRTRIM(FIX(N_periods), 1), free=1
      cgImage, BYTSCL(bg_image_new), CTIndex=0
      min_Neel_new = MIN(Neel_image_new[nonZero_ind])
      max_Neel_new = MAX(Neel_image_new[nonZero_ind])
      print, 'Neel range: ', min_Neel_new, ', ', max_Neel_new
;      cgImage, HIST_EQUAL(Neel_image_new, MINV=min_Neel_new, MAXV=max_Neel_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      cgImage, BYTSCL(Neel_image_new, MIN=min_Neel_new, MAX=max_Neel_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      cgImage, BYTSCL(Neel_image_new, MIN=min_Neel_new, MAX=max_Neel_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
       cgImage, BYTSCL(Neel_image_new, MIN=0, MAX=40, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'

      window, 2+w_flag*2
      cgDisplay, N_steps_new, N_steps_new, title='Brownian time map '+STRTRIM(FIX(N_periods), 1), free=1
      cgImage, BYTSCL(bg_image_new), CTIndex=0
      min_Brownian_new = MIN(Brownian_image_new[nonZero_ind])
      max_Brownian_new = MAX(Brownian_image_new[nonZero_ind])
      print, 'Brownian range: ', min_Brownian_new, ', ', max_Brownian_new
;      cgImage, BYTSCL(Brownian_image_new, MIN=min_Brownian_new, MAX=max_Brownian_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
      cgImage, BYTSCL(Brownian_image_new, MIN=0, MAX=40, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      cgImage, BYTSCL(Brownian_image_new, MIN=min_Brownian_new, MAX=max_Brownian_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      cgImage, HIST_EQUAL(Brownian_image_new, MINV=min_Brownian_new, MAXV=max_Brownian_new, TOP=255), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'   
      w_flag++
      
      ;save the data Neel map
      filename_System = sav_file_dir+'\Neel_'+STRTRIM(FIX(N_periods), 1)+'.dat'
      filename_dir = FILE_DIRNAME(filename_System)
      IF FILE_TEST(filename_dir) NE 1 THEN FILE_MKDIR, filename_dir
      openw,unit1,filename_System,/get_lun
      FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit1, Neel_image[i,j]
      FREE_LUN, unit1
      
      ;save the data Brownian map
      filename_System = sav_file_dir+'\Brownian_'+STRTRIM(FIX(N_periods), 1)+'.dat'
      filename_dir = FILE_DIRNAME(filename_System)
      IF FILE_TEST(filename_dir) NE 1 THEN FILE_MKDIR, filename_dir
      openw,unit2,filename_System,/get_lun
      FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit2, Brownian_image[i,j]
      FREE_LUN, unit2
    ENDIF
    
    ;Neel SNR and CNRs
    signal_catheter_Neel = MEAN(Neel_image[catheter_index])
    signal_plaque_Neel = MEAN(Neel_image[plaque_index])
    signal_thermia_Neel = MEAN(Neel_image[thermia_index])
    noise_Neel = STDDEV(Neel_image[0:20, 65:*])
    temp_image_hr = image_hr * 0.
    temp_image_hr[20:45, 65:*] = image_hr[20:45, 65:*]
    vessel_index = WHERE(temp_image_hr GT 170, vessel_count)
    signal_vessel_Neel = MEAN(Neel_image[vessel_index])
    
    filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
    filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw2.bmp'
    image_bmp = READ_BMP(filename_bmp)
    d_phantom = REFORM(image_bmp[0,*,*])
    ;  iimage, d_phantom
    size_dp = SIZE(d_phantom, /dim)
    d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], 100, 100, /center, /minus_one)
    d_phantom_sq[thermia_index]=0
    d_phantom_sq[vessel_index]=0
    d_phantom_sq[plaque_index]=0
    d_phantom_sq[catheter_index]=0
    tumor_unheated_index = WHERE(d_phantom_sq EQ 38 OR d_phantom_sq EQ 152, tumor_unheated_count)
    signal_tumor_unheated_Neel = MEAN(Neel_image[tumor_unheated_index])
    print, 'N_periods: ', N_periods, '; field amplitude: ', field_excitation
    print, 'vessel Neel: ', signal_vessel_Neel
    print, 'catheter Neel: ', signal_catheter_Neel
    print, 'plaque Neel: ', signal_plaque_Neel
    print, 'hyperthermia Neel: ', signal_thermia_Neel

;    result_array = [N_periods, SNR_vessel_Neel, CNR_catheter_Neel, CNR_plaque_Neel, CNR_thermia_Neel]
    S_vessel_array_Neel[k] = signal_vessel_Neel
    S_catheter_array_Neel[k] = signal_catheter_Neel
    S_plaque_array_Neel[k] = signal_plaque_Neel
    S_thermia_array_Neel[k] = signal_thermia_Neel
    
    ;Brownian SNR and CNRs
    signal_catheter_Brownian = MEAN(Brownian_image[catheter_index])
    signal_plaque_Brownian = MEAN(Brownian_image[plaque_index])
    signal_thermia_Brownian = MEAN(Brownian_image[thermia_index])
    noise_Brownian = STDDEV(Brownian_image[0:20, 65:*])
    signal_vessel_Brownian = MEAN(Brownian_image[vessel_index])
    signal_tumor_unheated_Brownian = MEAN(Brownian_image[tumor_unheated_index])
    print, 'N_periods: ', N_periods, '; field amplitude: ', field_excitation
    print, 'vessel Brownian: ', signal_vessel_Brownian
    print, 'catheter Brownian: ', signal_catheter_Brownian
    print, 'plaque Brownian: ', signal_plaque_Brownian
    print, 'hyperthermia Brownian: ', signal_thermia_Brownian

;    result_array = [N_periods, SNR_vessel_Brownian, CNR_catheter_Brownian, CNR_plaque_Brownian, CNR_thermia_Brownian]
    S_vessel_array_Brownian[k] = signal_vessel_Brownian
    S_catheter_array_Brownian[k] = signal_catheter_Brownian
    S_plaque_array_Brownian[k] = signal_plaque_Brownian
    S_thermia_array_Brownian[k] = signal_thermia_Brownian
  ENDFOR


  IF ABS(MAX(N_periods_array)-MIN(N_periods_array)) GT 1. THEN BEGIN
    sort_ind = SORT(N_periods_array)
    window, 3
    plot, N_periods_array[sort_ind], S_vessel_array_Neel[sort_ind], YRANGE = [0, CEIL(MAX(S_vessel_array_Brownian))], title='S Neel', xtitle='# of periods', ytitle='relaxation time'
    oplot, N_periods_array[sort_ind], S_catheter_array_Neel[sort_ind], LINESTYLE=1;, YRANGE = [0, CEIL(MAX(S_plaque_array_Neel))], title='SNR Neel', xtitle='# of periods', ytitle='SNR'
    oplot, N_periods_array[sort_ind], S_plaque_array_Neel[sort_ind], LINESTYLE=2
;    oplot, N_periods_array[sort_ind], S_thermia_array_Neel[sort_ind]
    window, 4
    plot, N_periods_array[sort_ind], S_vessel_array_Brownian[sort_ind], YRANGE = [0, CEIL(MAX(S_plaque_array_Brownian))], title='S Brownian', xtitle='# of periods', ytitle='relaxation time'
    oplot, N_periods_array[sort_ind], S_catheter_array_Brownian[sort_ind], LINESTYLE=1;, YRANGE = [0, CEIL(MAX(S_plaque_array_Brownian))], title='SNR Brownian', xtitle='# of periods', ytitle='SNR'
    oplot, N_periods_array[sort_ind], S_plaque_array_Brownian[sort_ind], LINESTYLE=2
;    oplot, N_periods_array[sort_ind], S_thermia_array_Brownian[sort_ind]
  ENDIF ELSE BEGIN
    sort_ind = SORT(N_fieldAmp_array)
    window, 3
    plot, N_fieldAmp_array[sort_ind], S_vessel_array_Neel[sort_ind], YRANGE = [0, CEIL(MAX(S_vessel_array_Brownian))], title='S Neel', xtitle='field amplitude', ytitle='relaxation time'
    oplot, N_fieldAmp_array[sort_ind], S_catheter_array_Neel[sort_ind], LINESTYLE=1;, YRANGE = [0, CEIL(MAX(S_plaque_array_Neel))], title='SNR Neel', xtitle='# of periods', ytitle='SNR'
    oplot, N_fieldAmp_array[sort_ind], S_plaque_array_Neel[sort_ind], LINESTYLE=2
;    oplot, N_fieldAmp_array[sort_ind], S_thermia_array_Neel[sort_ind]
    window, 4
    plot, N_fieldAmp_array[sort_ind], S_vessel_array_Brownian[sort_ind], YRANGE = [0, CEIL(MAX(S_plaque_array_Brownian))], title='S Brownian', xtitle='field amplitude', ytitle='relaxation time'
    oplot, N_fieldAmp_array[sort_ind], S_catheter_array_Brownian[sort_ind], LINESTYLE=1;, YRANGE = [0, CEIL(MAX(S_plaque_array_Brownian))], title='SNR Brownian', xtitle='# of periods', ytitle='SNR'
    oplot, N_fieldAmp_array[sort_ind], S_plaque_array_Brownian[sort_ind], LINESTYLE=2
  ENDELSE


  current_time = string(SYSTIME(/JULIAN, /UTC), FORMAT='(f8.0)')
  output_filename = sav_file_dir+current_time+'xls'
  ;write file
  openw,unit,output_filename,/get_lun
  first_line =['Field Amplitude', 'N_periods', 'Neel vessel', 'Neel catheter', 'Neel plaque', 'Neel hyperthermia', 'Brownian vessel', 'Brownian catheter', 'Brownian plaque', 'Brownian hyperthermia']
  printf, unit, FORMAT = '(12(A, %"\t"))', first_line
  FOR i=0L, nrfile-1L DO printf, unit, FORMAT = '(12(A, %"\t"))', N_fieldAmp_array[sort_ind[i]], N_periods_array[sort_ind[i]], $
      S_vessel_array_Neel[sort_ind[i]], S_catheter_array_Neel[sort_ind[i]], S_plaque_array_Neel[sort_ind[i]], S_thermia_array_Neel[sort_ind[i]], $
      S_vessel_array_Brownian[sort_ind[i]], S_catheter_array_Brownian[sort_ind[i]], S_plaque_array_Brownian[sort_ind[i]], S_thermia_array_Brownian[sort_ind[i]]
  ;; close file and END execution
  FREE_LUN, unit

END


FUNCTION debye_signal_analysis, fn
  ;6mT
  IF N_ELEMENTS(fn) EQ 0 THEN fn='C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\image_FFP_debye_signals.sav'
  restore, fn
  H_max = field_excitation
;  H_max = 6;mT
  
;  ;3mT
;  fn='C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\image_FFP_debye_signals3mT.sav'
;  H_max = 3;mT


  N_steps=N_ELEMENTS(FFP_signal_debye_2d_200ave[*,0,0])
  N_tp=N_ELEMENTS(FFP_signal_debye_2d_200ave[0,0,*])
  Debye_image=DBLARR(N_steps,N_steps)*0.
  
  ;GJ, 2022/2/11, define debby tao
  N_debye_tao = 600; 200
  ;GJ, 2022/1/30, calculate Debye model
  relaxation_time = FINDGEN(N_debye_tao)/10.+1.

  N_particle_size = 1
  particle_size_array = [24.4924]
  
  delta_t = 1.0
  time_rise = FINDGEN(N_tp)*delta_t
  H_array = FINDGEN(N_tp)/(N_tp-1)*2*H_max-H_max
  FOR i=0, N_steps-1 DO BEGIN
    FOR j=0, N_steps-1 DO BEGIN
      IF MEAN(FFP_signal_debye_2d_200ave[i,j,*]) LT 0.001 THEN continue
      ;do the biexponential curve fit
      ;Provide an initial guess of the function’s parameters.
      M_meas=TOTAL(REFORM(FFP_signal_debye_2d_200ave[i,j,*]), /CUMULATIVE) * delta_t
      
      MSE_array = DBLARR(N_particle_size, N_debye_tao)*0.+10000.
      correlation_array = DBLARR(N_particle_size, N_debye_tao)*0.
      FOR k=0, N_particle_size-1 DO BEGIN
        particle_size = particle_size_array[k]
        Msat_T = 0.551 ; T/u0
        Msat_kAm = 0.551/4./!PI*10000.; kA/m
        T_p = 20.;36.5; degree
        T_p_kelvin = 273.15 + T_p ; in Kelvin
        beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;        beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
        FWHM_particle = 4.16/beta_particle/3.5; mT/mm
        M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))

        FOR l = 0, N_debye_tao-1 DO BEGIN
          M_nat = non_adiabatic(time_rise, M, relaxation_time[l])
          M_nat -= MIN(M_nat)
          M_nat = M_nat/MEAN(M_nat)*MEAN(M_meas)
          correlation_array[k, l] = CORRELATE(M_nat, M_meas)
          MSE_array[k, l] = calMSE(M_nat, M_meas)
        ENDFOR
      ENDFOR
      
      min_MSE = MIN(MSE_array, min_ind)
      min_ind2d = ARRAY_INDICES(MSE_array, min_ind)

      tao_debye = relaxation_time[min_ind2d[1]]
      print, '[',i,', ',j,']: (Debye) relaxation time = ', tao_debye
      
      Debye_image[i,j] = tao_debye
    ENDFOR
  ENDFOR
  
  iimage, Debye_image, title='Debye time map'
  
  p_b_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\rec_image200_perpendicular_baseline.png'
  h_b_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\rec_image200_horizontal_baseline.png'
  image_hr = read_png(p_b_fn)/2.+read_png(h_b_fn)/2.
 
  cur_dir = FILE_DIRNAME(fn)
  cgDisplay, N_steps, N_steps, title='Debye', free=1
  cgImage, BYTSCL(image_hr), CTIndex=0
  cgImage, BYTSCL(Debye_image), CTIndex=33, Transparent=50;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;  cgImage, output='PNG', outfilename=cur_dir+'\Debye_map1.png'

  RETURN, Debye_image
END

;@GJ, 2022/8/22, do a biexponential analysis
FUNCTION decay_signal_analysis, fn
  ;6mT
  IF N_ELEMENTS(fn) EQ 0 THEN fn='C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\image_FFP_decay_signals.sav'
  ;;3mT
  ;fn='C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\image_FFP_decay_signals3mT.sav'
  
  restore, fn
  
  N_steps=N_ELEMENTS(FFP_signal_decay_2d_200ave[*,0,0])
  N_tp=N_ELEMENTS(FFP_signal_decay_2d_200ave[0,0,*])
  delta_t = 1.0
  Neel_image=DBLARR(N_steps,N_steps)*0.
  Brownian_image=DBLARR(N_steps,N_steps)*0.
  Neel_image_baseline=DBLARR(N_steps,N_steps)*0.
  Brownian_image_baseline=DBLARR(N_steps,N_steps)*0.
  
  FOR i=0, N_steps-1 DO BEGIN
    FOR j=0, N_steps-1 DO BEGIN
      IF MEAN(FFP_signal_decay_2d_200ave[i,j,*]) LT 0.001 THEN continue
      ;do the biexponential curve fit
      ;Provide an initial guess of the function’s parameters.
      Y=TOTAL(REFORM(FFP_signal_decay_2d_200ave[i,j,*]), /CUMULATIVE) * delta_t
      weights=1./Y
      X=FINDGEN(N_tp)*delta_t
      B = [-MAX(Y)/2.0, -0.01, -MAX(Y)/2.0, -0.1, MIN(Y), 0.]; ;@GJ, 2022/7/10, B = [-5.0, -0.01, -5.0, -0.1, MAX(Y), 0.]
      FITA = [1,1,1,1,1,0]
      ;Compute the parameters.
      yfit_biex = CURVEFIT(X, Y, weights, B, SIGMA, FITA=FITA, FUNCTION_NAME='gfunct_biex')
      bi_tao_1 = 1./ABS(B[1])
      bi_tao_2 = 1./ABS(B[3])
      print, 'pixel = [',i,',',j,',]'
      print, 'Relaxation time_1 = ', bi_tao_1, ' [us]'
      print, 'Relaxation time_2 = ', bi_tao_2, ' [us]'
      Neel_image[i,j] = MIN([bi_tao_1, bi_tao_2])
      Brownian_image[i,j] = MAX([bi_tao_1, bi_tao_2])
      
      ;do the biexponential curve fit
      ;Provide an initial guess of the function’s parameters.
      Y=TOTAL(REFORM(FFP_signal_decay_2d_200ave_baseline[i,j,*]), /CUMULATIVE) * delta_t
      weights=1./Y
      X=FINDGEN(N_tp)*delta_t
      B = [-MAX(Y)/2.0, -0.01, -MAX(Y)/2.0, -0.1, MIN(Y), 0.]; ;@GJ, 2022/7/10, B = [-5.0, -0.01, -5.0, -0.1, MAX(Y), 0.]
      FITA = [1,1,1,1,1,0]
      ;Compute the parameters.
      yfit_biex = CURVEFIT(X, Y, weights, B, SIGMA, FITA=FITA, FUNCTION_NAME='gfunct_biex')
      bi_tao_1 = 1./ABS(B[1])
      bi_tao_2 = 1./ABS(B[3])
      print, 'pixel = [',i,',',j,',]'
      print, 'Relaxation time_1 = ', bi_tao_1, ' [us]'
      print, 'Relaxation time_2 = ', bi_tao_2, ' [us]'
      Neel_image_baseline[i,j] = MIN([bi_tao_1, bi_tao_2])
      Brownian_image_baseline[i,j] = MAX([bi_tao_1, bi_tao_2])
    ENDFOR
  ENDFOR
  iimage, Neel_image, title='Neel time map'
  iimage, Brownian_image, title='Brownian time map'
  ;iimage, ALOG(Brownian_image), title='Alog Brownian time map'
  iimage, Neel_image_baseline, title='Neel time map baseline'
  iimage, Brownian_image_baseline, title='Brownian time map baseline'
  
  p_b_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\rec_image200_perpendicular_baseline.png'
  h_b_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\rec_image200_horizontal_baseline.png'
  image_hr = read_png(p_b_fn)/2.+read_png(h_b_fn)/2.

  cgDisplay, N_steps, N_steps, title='Brownian', free=1
  cgImage, BYTSCL(image_hr), CTIndex=0
  cgImage, BYTSCL(Brownian_image), CTIndex=33, Transparent=50
  
  cgDisplay, N_steps, N_steps, title='Neel', free=1
  cgImage, BYTSCL(image_hr), CTIndex=0
  cgImage, BYTSCL(Neel_image), CTIndex=33, Transparent=50
  
  RETURN, [Neel_image, Brownian_image]
  
;  neel_bmp_fn='C:\D_drive\MPI_Tianjie\Xinfeng\jiamingana\001.bmp'
;  Neel_image_jiaming=read_bmp(neel_bmp_fn, /rgb)
;  cgDisplay, N_steps, N_steps, title='Neel Jiaming', free=1
;  cgImage, BYTSCL(image_hr), CTIndex=0
;  cgImage, BYTSCL(Neel_image_jiaming), CTIndex=33, Transparent=50
  ;  total_image = image_hr+catheter_image+thermia_image+plaque_image
  ;  Image_Blend, BYTSCL(total_image), BYTSCL(Brownian_image), Colortable=25;, blendTitle='Neel'

  ;  ;  Image_Blend, BYTSCL(REFORM(rec_image[0,*,*])), BYTSCL(REFORM(rec_image[2,*,*]), MAX=30), Colortable=25
  ;
  ;  ;GJ, 2022/5/10, should not include the t=0 signal
  ;  filename = cur_dir+'\pixel6515_FFP_signal'
  ;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
  ;
  ;  y_data_int_baseline = INTERPOL(REFORM(FFP_signal_2d_200ave_baseline[65, 15, 26:*]), FINDGEN(N_tp-26), X_data_int)
  ;  ;GJ, 2022/5/10, should not include the t=0 signal
  ;  filename = cur_dir+'\pixel6515_FFP_signal_baseline'
  ;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int_baseline[1:*]
;  cur_dir = FILE_DIRNAME(fn)
;  N_int = 8000.
;  X_data_int = FINDGEN(N_int+1)/N_int*MAX(X)
;
;  ;hyperthermai
;  i=62
;  j=18
;  signal_array = REFORM(FFP_signal_decay_2d_200ave[i,j,*])
;  y_data_ori = signal_array*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
;  ;  ;GJ, 2022/5/10, should not include the t=0 signal
;  filename = cur_dir+'\pixel6218thermia_FFP_signal'
;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
;
;  ;catheter
;  i=20
;  j=38
;  signal_array = REFORM(FFP_signal_decay_2d_200ave[i,j,*])
;  y_data_ori = signal_array*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
;  ;  ;GJ, 2022/5/10, should not include the t=0 signal
;  filename = cur_dir+'\pixel2038cath_FFP_signal'
;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
;  
;  ;vessel
;  i=30
;  j=62
;  signal_array = REFORM(FFP_signal_decay_2d_200ave[i,j,*])
;  y_data_ori = signal_array*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
;  ;  ;GJ, 2022/5/10, should not include the t=0 signal
;  filename = cur_dir+'\pixel3062vessel_FFP_signal'
;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
;  
;  ;plaque
;  i=60
;  j=61
;  signal_array = REFORM(FFP_signal_decay_2d_200ave[i,j,*])
;  y_data_ori = signal_array*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
;  ;  ;GJ, 2022/5/10, should not include the t=0 signal
;  filename = cur_dir+'\pixel6061plaque_FFP_signal'
;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
;  
;  ;tumor
;  i=62
;  j=36
;  signal_array = REFORM(FFP_signal_decay_2d_200ave[i,j,*])
;  y_data_ori = signal_array*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
;  ;  ;GJ, 2022/5/10, should not include the t=0 signal
;  filename = cur_dir+'\pixel6236tumor_FFP_signal'
;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
;  
;  ;thermia
;  i=62
;  j=18
;  signal_array = REFORM(FFP_signal_decay_2d_200ave[i,j,*])
;  y_data_ori = signal_array*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
;  ;  ;GJ, 2022/5/10, should not include the t=0 signal
;  filename = cur_dir+'\pixel6218thermia_FFP_signal'
;  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*]
END

;@GJ, 2022/7/27, analyze signal system matrix based on the number of field strengths
;200 averages as system matrix, 200 averages as signal matrix
;@GJ, 2022/8/13, add the catheter with 50 degree
;@GJ, 2022/8/13, open the 70% glycerol of Synomag-D
;@GJ, 2022/8/18, add the catheter with 30% Gelatin, by changing the flag, catheter or other can be evaluated
;@GJ, 2022/9/22, change the vessel to 1% glycerol, the plaque to 40% glycerol. Only calculate catheter

PRO read_signal_sav_hrNew_SP_catheter_plaque_thermia, N_periods, cur_dir, correlation_array, MSE_array, stripe_stdev_h, stripe_stdev_p
  ;open the baseline 25 deg of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Baseline sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\20220817121631signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_0 = signal_matrix
  signal_matrix_sp_0 = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 50 degree of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T50deg\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Hyperthermia sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T50deg\20220813100000signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_50d = signal_matrix
  signal_matrix_sp_50d = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 70% glycerol of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly70\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Plaque sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly70\20220813121126signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_g70 = signal_matrix
  signal_matrix_sp_g70 = signal_matrix_sp

  ;@GJ, 2022/9/16, adding the baseline as g1
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly1\20220813120203signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_g1 = signal_matrix
  signal_matrix_sp_g1 = signal_matrix_sp

  ;@GJ, 2022/9/16, adding the baseline as g1
;  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_Gly40\20220916162759signalNewSP.sav'
;  RESTORE, sav_name
;  signal_matrix_g40 = signal_matrix
;  signal_matrix_sp_g40 = signal_matrix_sp

  ;@GJ, 2022/8/13, open the 30% Gelatin of Synomag-D
  ;directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\'
  ;sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Catheter sav file", FILTER='*.sav')
  sav_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220815_Ge\Synomag_Ge30\20220817100728signalNewSP.sav'
  RESTORE, sav_name
  signal_matrix_ge30 = signal_matrix
  signal_matrix_sp_ge30 = signal_matrix_sp

  ;@GJ, 2022/7/27, save the image
  time_rise_flat = REFORM(signal_matrix_0[0, *])
  N_fields = N_ELEMENTS(signal_matrix_0[*,0])-1
  N_tp = N_ELEMENTS(signal_matrix_0[0,*])
  N_steps = 100

  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  ;@GJ, 2022/9/15, removing the tumor region
  phantom_image[*, 320:*, 0:370] = 0.
  image_hr = CONGRID(REFORM(phantom_image[0,*,*]), N_steps, N_steps)
  ;  iimage, image_hr, title='original image (high resolution)'
  
  base_fn = FILE_BASENAME(p_fn, 'jpg')
  FILE_MKDIR, 'C:\D_drive\MPI_Tianjie\Xinfeng\reconst\'+base_fn+'\'
  cur_dir = 'C:\D_drive\MPI_Tianjie\Xinfeng\reconst\'+base_fn+'\'
  filename_ori_image = cur_dir+'\ori_image.png'
  WRITE_PNG, filename_ori_image, image_hr
  filename_System = cur_dir+'\original_image.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit1, image_hr[i, j]
  FREE_LUN, unit1
  
  filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
  image_bmp = READ_BMP(filename_bmp)
  d_phantom = REFORM(image_bmp[0,*,*])
  ;  iimage, d_phantom
  size_dp = SIZE(d_phantom, /dim)
  d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], N_steps, N_steps, /center, /minus_one)
  d_phantom_sq[*,0] = d_phantom_sq[*,1]
  d_phantom_sq[*,N_steps-1] = d_phantom_sq[*,N_steps-2]
  d_phantom_sq[N_steps-1,*] = d_phantom_sq[N_steps-2,*]
  ;extract catheter part
  catheter_index = WHERE(d_phantom_sq GT 240, catheter_count)
  catheter_index_2d = ARRAY_INDICES(d_phantom_sq, catheter_index)
  d_phantom_sq[0:40,*] = 0
  catheter_image = d_phantom_sq*0.
  catheter_image[catheter_index] = 256
  ;  iimage, image_hr+catheter_image, title='baseline + catheter'
  ;extract plaque part
  plaque_index = WHERE(d_phantom_sq GT 200, plauqe_count)
  plaque_index_2d = ARRAY_INDICES(d_phantom_sq, plaque_index)
  d_phantom_sq[*,26:*] = 0
  plaque_image = d_phantom_sq*0.
  plaque_image[plaque_index] = 256
  ; iimage, image_hr+plaque_image, title='baseline + plaque'
  ;extract thermia part
  thermia_index = WHERE(d_phantom_sq GT 0, thermia_count)
  thermia_index_2d = ARRAY_INDICES(d_phantom_sq, thermia_index)
  thermia_image = d_phantom_sq*0.
  thermia_image[thermia_index] = 256
  ; iimage, image_hr+thermia_image, title='baseline + thermia'
  total_image = image_hr+catheter_image;+plaque_image+thermia_image
  
  field_array = (FINDGEN(N_steps)+1.)/N_steps*10.
  
  IF N_ELEMENTS(N_periods) EQ 0 THEN N_periods = 200
  
  ;1% glycerol as the baseline
  signal_matrix_200ave = REFORM(MEAN(signal_matrix_sP_g1[*, 0:N_periods-1, *], DIMENSION=2));i
  signal_matrix_hr_200ave = CONGRID(signal_matrix_200ave, N_steps, N_tp, /CENTER, /INTERP)
  ;  signal_matrix_hr_400ave = CONGRID(signal_matrix_0[1:N_fields, *], N_steps, N_tp, /CENTER, /INTERP)

  ;50deg
  signal_matrix_200ave_50d = REFORM(MEAN(signal_matrix_sP_50d[*, 0:N_periods-1, *], DIMENSION=2));i
  signal_matrix_hr_200ave_50d = CONGRID(signal_matrix_200ave_50d, N_steps, N_tp, /CENTER, /INTERP)
  
  ;70% glycerol
  signal_matrix_200ave_g70 = REFORM(MEAN(signal_matrix_sP_g70[*, 0:N_periods-1, *], DIMENSION=2));i
  signal_matrix_hr_200ave_g70 = CONGRID(signal_matrix_200ave_g70, N_steps, N_tp, /CENTER, /INTERP)

;  ;40% glycerol
;  signal_matrix_200ave_g40 = REFORM(MEAN(signal_matrix_sP_g40[*, 0:N_periods-1, *], DIMENSION=2));i
;  signal_matrix_hr_200ave_g40 = CONGRID(signal_matrix_200ave_g40, N_steps, N_tp, /CENTER, /INTERP)

  ;30% gelatin
  signal_matrix_200ave_ge30 = REFORM(MEAN(signal_matrix_sP_ge30[*, 0:N_periods-1, *], DIMENSION=2));i
  signal_matrix_hr_200ave_ge30 = CONGRID(signal_matrix_200ave_ge30, N_steps, N_tp, /CENTER, /INTERP)


  ;  signal_2d_hr_h_400ave = DBLARR(N_steps, N_tp)*0.
  signal_2d_hr_h_200ave = DBLARR(N_steps, N_tp)*0.
  signal_2d_hr_h_200ave_baseline = DBLARR(N_steps, N_tp)*0.
  FOR i=0, N_steps-1 DO BEGIN
    image_1d_hr_h = REFORM(image_hr[*, i])
    FOR j=0, N_steps-1 DO BEGIN
      ;      signal_2d_hr_h_400ave[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_400ave[j, *])
      flag = 0
      FOR k=0, N_ELEMENTS(catheter_index)-1 DO BEGIN
        IF j EQ catheter_index_2d[0, k] AND i EQ catheter_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 1;0;1; for catheter
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(plaque_index)-1 DO BEGIN
        IF j EQ plaque_index_2d[0, k] AND i EQ plaque_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 0;2; for plaque
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(thermia_index)-1 DO BEGIN
        IF j EQ thermia_index_2d[0, k] AND i EQ thermia_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 0;3;0;3; for thermia
          break
        ENDIF
      ENDFOR
      IF flag EQ 1 THEN BEGIN
        ;catheter
        signal_2d_hr_h_200ave[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_200ave_ge30[j, *])
      ENDIF
      IF flag EQ 2 THEN BEGIN
        ;plaque
        signal_2d_hr_h_200ave[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_200ave_g40[j, *])
      ENDIF
      IF flag EQ 3 THEN BEGIN
        ;thermia
        signal_2d_hr_h_200ave[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_200ave_50d[j, *])
      ENDIF
      IF flag EQ 0 THEN BEGIN
        signal_2d_hr_h_200ave[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_200ave[j, *])
      ENDIF
      ;;the baseline signal
      signal_2d_hr_h_200ave_baseline[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_200ave[j, *])
    ENDFOR
  ENDFOR

  ;  signal_2d_hr_p_400ave = DBLARR(N_steps_h, N_tp)*0.
  signal_2d_hr_p_200ave = DBLARR(N_steps, N_tp)*0.
  signal_2d_hr_p_200ave_baseline = DBLARR(N_steps, N_tp)*0.
  FOR i=0, N_steps-1 DO BEGIN
    image_1d_hr_p = REFORM(image_hr[i, *]); perpendicular direction
    FOR j=0, N_steps-1 DO BEGIN
      flag = 0
      FOR k=0, N_ELEMENTS(catheter_index)-1 DO BEGIN
        IF i EQ catheter_index_2d[0, k] AND j EQ catheter_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 1;0;1; for catheter
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(plaque_index)-1 DO BEGIN
        IF i EQ plaque_index_2d[0, k] AND j EQ plaque_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 0;2; for plaque
          break
        ENDIF
      ENDFOR
      FOR k=0, N_ELEMENTS(thermia_index)-1 DO BEGIN
        IF i EQ thermia_index_2d[0, k] AND j EQ thermia_index_2d[1, k] THEN BEGIN
          ;found and break
          flag = 0;3; for thermia
          break
        ENDIF
      ENDFOR
      IF flag EQ 1 THEN BEGIN
        ;catheter
        signal_2d_hr_p_200ave[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr_200ave_ge30[j, *])
      ENDIF
      IF flag EQ 2 THEN BEGIN
        ;plaque
        signal_2d_hr_p_200ave[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr_200ave_g40[j, *])
      ENDIF
      IF flag EQ 3 THEN BEGIN
        ;thermia
        signal_2d_hr_p_200ave[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr_200ave_50d[j, *])
      ENDIF
      IF flag EQ 0 THEN BEGIN
        signal_2d_hr_p_200ave[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr_200ave[j, *])
      ENDIF
      ;the baseline signal
      signal_2d_hr_p_200ave_baseline[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr_200ave[j, *])
    ENDFOR
  ENDFOR


  ;@GJ, 2022/7/31, generate the system matrix
  offset=10
  ;  system_matrix_400ave = signal_matrix_hr_400ave[*, 0+offset:offset+N_steps-1]
  ;@GJ, 2022/9/22, use the half 200 periods as system matrix
  ;@GJ, 2022/9/24, use the 1st 200 periods as system matrix
  signal_matrix_200ave_h = REFORM(MEAN(signal_matrix_sP_g1[*, 0:199, *], DIMENSION=2));i
  signal_matrix_hr_200ave_h = CONGRID(signal_matrix_200ave_h, N_steps, N_tp, /CENTER, /INTERP)
  system_matrix_200ave = signal_matrix_hr_200ave_h[*, 0+offset:offset+N_steps-1]

  ;save the results
  ;@GJ, 2022/8/3, find the current directory
  ;  ;write file
  string_ave = STRTRIM(FIX(N_periods), 1)
  filename_System = cur_dir+'\systemMatrix_2d_array_200ave.dat'
  filename_Signal_h = cur_dir+'\signal_2d_h_array_'+ string_ave+'ave.dat'
  filename_Signal_p = cur_dir+'\signal_2d_p_array_'+ string_ave+'ave.dat'
  filename_Signal_h_baseline = cur_dir+'\signal_2d_h_array_'+ string_ave+'ave_baseline.dat'
  filename_Signal_p_baseline = cur_dir+'\signal_2d_p_array_'+ string_ave+'ave_baseline.dat'
  openw,unit1,filename_System,/get_lun
  openw,unit2,filename_Signal_h,/get_lun
  openw,unit3,filename_Signal_p,/get_lun
  openw,unit4,filename_Signal_h_baseline,/get_lun
  openw,unit5,filename_Signal_p_baseline,/get_lun
  FOR i=0, N_steps-1 DO BEGIN
    FOR j=0, N_steps-1 DO BEGIN
      ;; write the data points
      WRITEU, unit1, (system_matrix_200ave)[i, j]
      WRITEU, unit2, (REFORM(signal_2d_hr_h_200ave[*, 0+offset:offset+N_steps-1]))[i, j]
      WRITEU, unit3, (REFORM(signal_2d_hr_p_200ave[*, 0+offset:offset+N_steps-1]))[j, i]
      WRITEU, unit4, (REFORM(signal_2d_hr_h_200ave_baseline[*, 0+offset:offset+N_steps-1]))[i, j]
      WRITEU, unit5, (REFORM(signal_2d_hr_p_200ave_baseline[*, 0+offset:offset+N_steps-1]))[j, i]
    ENDFOR
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit1, unit2, unit3, unit4, unit5

  ;@GJ, 2022/8/3, start reconstruction
  ;@GJ, 2022/8/3, horizontal direction
  ;  rec_image_hr_h_400ave = image_hr*0.
  rec_image_hr_h_200ave = image_hr*0.
  rec_image_hr_h_200ave_baseline = image_hr*0.
  FOR i=0, N_steps-1 DO BEGIN
    ;@GJ, 2022/8/3, system_matrix_400ave as the standard system matrix
    ;   rec_image_hr_h_400ave[*, i] = REFORM(art_func(system_matrix_400ave_h, REFORM(signal_2d_hr_h_400ave[i, 0+offset:offset+N_steps-1]))); REFORM(Gradient_invert ## REFORM(signal_2d[i, 0+offset:offset+N_fields-1]))
    rec_image_hr_h_200ave[*, i] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_h_200ave[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
    rec_image_hr_h_200ave_baseline[*, i] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_h_200ave_baseline[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
  ENDFOR
;  iimage, rec_image_hr_h_200ave, title='200 ave horizontal'
  filename_rec_image = cur_dir+'\rec_image'+ string_ave+'_horizontal.png'
  WRITE_PNG, filename_rec_image, BYTSCL(rec_image_hr_h_200ave, MAX=255, MIN=0)
  filename_rec_image_baseline = cur_dir+'\rec_image'+ string_ave+'_horizontal_baseline.png'
  WRITE_PNG, filename_rec_image_baseline, BYTSCL(rec_image_hr_h_200ave_baseline, MAX=255, MIN=0)

  ;perpendicular direction
  ;  rec_image_hr_p_400ave = image_hr*0.
  rec_image_hr_p_200ave = image_hr*0.
  rec_image_hr_p_200ave_baseline = image_hr*0.
  FOR i=0, N_steps-1 DO BEGIN
    ;@GJ, 2022/8/3, system_matrix_400ave as the standard system matrix
    ;    rec_image_hr_p_400ave[i, *] = REFORM(art_func(system_matrix_400ave_p, REFORM(signal_2d_hr_p_400ave[i, 0+offset:offset+N_steps-1]))); REFORM(Gradient_invert ## REFORM(signal_2d[i, 0+offset:offset+N_fields-1]))
    rec_image_hr_p_200ave[i, *] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_p_200ave[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
    rec_image_hr_p_200ave_baseline[i, *] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_p_200ave_baseline[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
  ENDFOR
;  iimage, rec_image_hr_p_200ave, title='200 ave perpendicular'
  filename_rec_image = cur_dir+'\rec_image'+ string_ave+'_perpendicular.png'
  WRITE_PNG, filename_rec_image, BYTSCL(rec_image_hr_p_200ave, MAX=255, MIN=0)
  filename_rec_image_baseline = cur_dir+'\rec_image'+ string_ave+'_perpendicular_baseline.png'
  WRITE_PNG, filename_rec_image_baseline, BYTSCL(rec_image_hr_p_200ave_baseline, MAX=255, MIN=0)

  ;  iimage, rec_image_lr, title='from ART'
  all_rec_image_hr = DBLARR(N_steps*4+10*3, N_steps*3+10*2)
  ;save the reconstructed
  all_rec_image_hr[0:N_steps-1, 2*N_steps+20:*] = rec_image_hr_h_200ave_baseline
  all_rec_image_hr[1*N_steps+10:2*N_steps-1+10, 2*N_steps+20:*] = rec_image_hr_p_200ave_baseline
  all_rec_image_hr[2*N_steps+20:3*N_steps-1+20, 2*N_steps+20:*] = (rec_image_hr_h_200ave_baseline + rec_image_hr_p_200ave_baseline)/2.
  all_rec_image_hr[3*N_steps+30:4*N_steps-1+30, 2*N_steps+20:*] = image_hr
  ;save the reconstructed baseline
  all_rec_image_hr[0:N_steps-1, N_steps+10:2*N_steps-1+10] = rec_image_hr_h_200ave
  all_rec_image_hr[1*N_steps+10:2*N_steps-1+10, N_steps+10:2*N_steps-1+10] = rec_image_hr_p_200ave
  all_rec_image_hr[2*N_steps+20:3*N_steps-1+20, N_steps+10:2*N_steps-1+10] = (rec_image_hr_h_200ave + rec_image_hr_p_200ave)/2.
  all_rec_image_hr[3*N_steps+30:4*N_steps-1+30, N_steps+10:2*N_steps-1+10] = image_hr+catheter_image;+plaque_image;+thermia_image;+catheter_image;
  ;save the difference
  all_rec_image_hr[*, 0:N_steps-1] = BYTSCL(all_rec_image_hr[*, 2*N_steps+20:*]-all_rec_image_hr[*, N_steps+10:2*N_steps-1+10], MAX=50, MIN=-50)
;  iimage, all_rec_image_hr, title='reconstructed images ART'+ string_ave
  filename_rec_image1 = cur_dir+'\rec_image'+ string_ave+'_ART.png'
  WRITE_PNG, filename_rec_image1, BYTSCL(all_rec_image_hr, MAX=255, MIN=0)
  
  correlation_array = DBLARR(6)*0.
  MSE_array = DBLARR(6)*0.
  correlation_array[0] = CORRELATE(image_hr, rec_image_hr_h_200ave_baseline)
  MSE_array[0] = calMSE(image_hr, rec_image_hr_h_200ave_baseline)
  correlation_array[1] = CORRELATE(image_hr, rec_image_hr_p_200ave_baseline)
  MSE_array[1] = calMSE(image_hr, rec_image_hr_p_200ave_baseline)
  correlation_array[2] = CORRELATE(image_hr, (rec_image_hr_h_200ave_baseline + rec_image_hr_p_200ave_baseline)/2.)
  MSE_array[2] = calMSE(image_hr, (rec_image_hr_h_200ave_baseline + rec_image_hr_p_200ave_baseline)/2.)
  correlation_array[3] = CORRELATE(image_hr, rec_image_hr_h_200ave)
  MSE_array[3] = calMSE(image_hr, rec_image_hr_h_200ave)
  correlation_array[4] = CORRELATE(image_hr, rec_image_hr_p_200ave)
  MSE_array[4] = calMSE(image_hr, rec_image_hr_p_200ave)
  correlation_array[5] = CORRELATE(image_hr, (rec_image_hr_h_200ave + rec_image_hr_p_200ave)/2.)
  MSE_array[5] = calMSE(image_hr, (rec_image_hr_h_200ave + rec_image_hr_p_200ave)/2.)
;  print, 'correlation: ', correlation_array
;  print, 'MSE: ', MSE_array
  ;  average_arr = [1, 400]
  correlation_HandP = correlation_array[5]
  MSE_HandP = MSE_array[5]
  ;
;  iplot, correlation_array, title = 'Correlation Coefficient', xtitle='# of signal averages', ytitle='Coefficient'
;  iplot, MSE_array, title = 'MSE', xtitle='# of signal averages', ytitle='MSE'
  ;
  
  rec_image_sub_h = rec_image_hr_h_200ave_baseline-rec_image_hr_h_200ave
  stripe_stdev_h = STDDEV(rec_image_sub_h[*, MIN(catheter_index_2d[1,*]):MAX(catheter_index_2d[1,*])])
  rec_image_sub_p = rec_image_hr_p_200ave_baseline-rec_image_hr_p_200ave
  stripe_stdev_p = STDDEV(rec_image_sub_h[MIN(catheter_index_2d[0,*]):MAX(catheter_index_2d[0,*]), *])
  print, 'N periods for average: ', N_periods
  print, 'stripe stddev horizontal: ', stripe_stdev_h
  print, 'stripe stddev perpendicular: ', stripe_stdev_p
END

;@GJ, 2022/9/24, batch analyze the stripe aritfacts
;@GJ, 2023/05/09, reconstruct images using different phantom
PRO batch_read_signal_sav_hrNew_SP

  num_k = 21;3;41
  k_array = FINDGEN(num_k)*10./(num_k-1)*20.
  k_array[0]=1
  
  stripe_stdev_h_array = DBLARR(num_k)
  stripe_stdev_p_array = DBLARR(num_k)
  correlation_h_array = DBLARR(num_k)
  correlation_p_array = DBLARR(num_k)
  FOR k=0, N_ELEMENTS(k_array)-1 DO BEGIN
    read_signal_sav_hrNew_SP_catheter_plaque_thermia, k_array[k], cur_dir, correlation_array, MSE_array, stripe_stdev_h, stripe_stdev_p
    stripe_stdev_h_array[k] = stripe_stdev_h
    stripe_stdev_p_array[k] = stripe_stdev_p
    correlation_h_array[k] = correlation_array[0]
    correlation_p_array[k] = correlation_array[1]
  ENDFOR

  iplot, k_array, stripe_stdev_h_array, title = 'Horizontal Stripe Stddev', xtitle='# of periods for signal averages', ytitle='Stripe Stddev'
  iplot, k_array, stripe_stdev_p_array, title = 'Vertical Stripe Stddev', xtitle='# of periods for signal averages', ytitle='Stripe Stddev'
  
  current_time = string(SYSTIME(/JULIAN, /UTC), FORMAT='(f8.0)')
  output_filename = cur_dir+'\StripeAnalysis_'+current_time+'xls'
  ;write file
  openw,unit,output_filename,/get_lun
  first_line =['N_periods', 'stripe_stdev_horizontal', 'stripe_stdev_vertical', 'correlation_horizontal', 'correlation_vertical']
  printf, unit, FORMAT = '(10(A, %"\t"))', first_line
  FOR k=0L, N_ELEMENTS(k_array)-1 DO BEGIN
    printf, unit, FORMAT = '(10(A, %"\t"))', k_array[k], stripe_stdev_h_array[k], stripe_stdev_p_array[k], correlation_h_array[k], correlation_p_array[k]
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit
  
END


PRO read_signal_baseline_reconstruction
  ;open the baseline 25 deg of Synomag-D
  directory = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220720_T2\T25deg\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="Baseline sav file", FILTER='*SP.sav')
  RESTORE, sav_name
  signal_matrix_0 = signal_matrix
  signal_matrix_sp_0 = signal_matrix_sp

  ;@GJ, 2022/7/27, save the image
  time_rise_flat = REFORM(signal_matrix_0[0, *])
  N_fields = N_ELEMENTS(signal_matrix_0[*,0])-1
  N_tp = N_ELEMENTS(signal_matrix_0[0,*])
  N_steps = 100

  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
  image_hr = CONGRID(REFORM(phantom_image[0,*,*]), N_steps, N_steps)
  ;  iimage, image_hr, title='original image (high resolution)'
  cur_dir = FILE_DIRNAME(p_fn)
  base_p_fn = FILE_BASENAME(p_fn, 'jpg')
  filename_ori_image = cur_dir+'\ori_image.png'
  WRITE_PNG, filename_ori_image, image_hr
  filename_System = cur_dir+'\original_image.dat'
  openw,unit1,filename_System,/get_lun
  FOR i=0, N_steps-1 DO FOR j=0, N_steps-1 DO WRITEU, unit1, image_hr[i, j]
  FREE_LUN, unit1

;  signal_matrix_200ave = REFORM(MEAN(signal_matrix_sP_0[*, 0:199, *], DIMENSION=2));i
;  signal_matrix_hr_200ave = CONGRID(signal_matrix_200ave, N_steps, N_tp, /CENTER, /INTERP)
  ;@GJ, 2022/7/31, generate the system matrix
  offset=10
  signal_matrix_200ave_h = REFORM(MEAN(signal_matrix_sP_0[*, 0:199, *], DIMENSION=2));i
  signal_matrix_hr_200ave_h = CONGRID(signal_matrix_200ave_h, N_steps, N_tp, /CENTER, /INTERP)
  system_matrix_200ave = signal_matrix_hr_200ave_h[*, 0+offset:offset+N_steps-1]
  
  k_array=INDGEN(41)*10
  k_array[0]=1
  MSE_array = DBLARR(41)*0.
  SNR_array = DBLARR(41)*0.
  correlation_array = DBLARR(41)*0.
  
;  FOR k=0, N_ELEMENTS(k_array)-1 DO BEGIN
   FOR k=20, 20 DO BEGIN
    
    signal_matrix = REFORM(MEAN(signal_matrix_sP_0[*, 0:(k_array[k]-1), *], DIMENSION=2));i
    signal_matrix_hr = CONGRID(signal_matrix, N_steps, N_tp, /CENTER, /INTERP)
    
    signal_2d_hr_h = DBLARR(N_steps, N_tp)*0.
;    signal_2d_hr_h_200ave_baseline = DBLARR(N_steps, N_tp)*0.
    FOR i=0, N_steps-1 DO BEGIN
      image_1d_hr_h = REFORM(image_hr[*, i])
      FOR j=0, N_steps-1 DO BEGIN
        signal_2d_hr_h[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr[j, *]) 
        ;;the baseline signal
;        signal_2d_hr_h_200ave_baseline[i, *] += image_1d_hr_h[j] * REFORM(signal_matrix_hr_200ave[j, *])
      ENDFOR
    ENDFOR
  
    signal_2d_hr_p = DBLARR(N_steps, N_tp)*0.
;    signal_2d_hr_p_200ave_baseline = DBLARR(N_steps, N_tp)*0.
    FOR i=0, N_steps-1 DO BEGIN
      image_1d_hr_p = REFORM(image_hr[i, *]); perpendicular direction
      FOR j=0, N_steps-1 DO BEGIN
        signal_2d_hr_p[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr[j, *])
        ;the baseline signal
;        signal_2d_hr_p_200ave_baseline[i, *] += image_1d_hr_p[j] * REFORM(signal_matrix_hr_200ave[j, *])
      ENDFOR
    ENDFOR
  
    ;@GJ, 2022/8/3, start reconstruction
    ;@GJ, 2022/8/3, horizontal direction
    rec_image_hr_h = image_hr*0.
;    rec_image_hr_h_200ave_baseline = image_hr*0.
    FOR i=0, N_steps-1 DO BEGIN
      ;@GJ, 2022/8/3, system_matrix_200ave as the standard system matrix
      rec_image_hr_h[*, i] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_h[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
;      rec_image_hr_h_200ave_baseline[*, i] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_h_200ave_baseline[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
    ENDFOR
    filename_rec_image = cur_dir+'\'+STRTRIM(k_array[k], 1)+'_rec_image'+STRTRIM(k_array[k], 1)+'_horizontal.png'
    WRITE_PNG, filename_rec_image, BYTSCL(rec_image_hr_h, MAX=255, MIN=0)
;    filename_rec_image_baseline = cur_dir+'\rec_image200_horizontal_baseline.png'
;    WRITE_PNG, filename_rec_image_baseline, BYTSCL(rec_image_hr_h_200ave_baseline, MAX=255, MIN=0)
  
    ;perpendicular direction
    rec_image_hr_p = image_hr*0.
;    rec_image_hr_p_200ave_baseline = image_hr*0.
    FOR i=0, N_steps-1 DO BEGIN
      ;@GJ, 2022/8/3, system_matrix_200ave as the standard system matrix
      rec_image_hr_p[i, *] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_p[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
 ;     rec_image_hr_p_200ave_baseline[i, *] = REFORM(art_func(system_matrix_200ave, REFORM(signal_2d_hr_p_200ave_baseline[i, 0+offset:offset+N_steps-1])));REFORM(Gradient_invert ## REFORM(signal_2d_1ave[i, 0+offset:offset+N_fields-1]))
    ENDFOR
    filename_rec_image = cur_dir+'\'+STRTRIM(k_array[k], 1)+'_rec_image'+'_perpendicular.png'
    WRITE_PNG, filename_rec_image, BYTSCL(rec_image_hr_p, MAX=255, MIN=0)
    filename_rec_image = cur_dir+'\'+STRTRIM(k_array[k], 1)+'_rec_image'+STRTRIM(k_array[k], 1)+'_HandP.png'
    WRITE_PNG, filename_rec_image, BYTSCL(rec_image_hr_h/2. + rec_image_hr_p/2., MAX=255, MIN=0)
;    filename_rec_image_baseline = cur_dir+'\rec_image200_perpendicular_baseline.png'
;    WRITE_PNG, filename_rec_image_baseline, BYTSCL(rec_image_hr_p_200ave_baseline, MAX=255, MIN=0)
    
    correlation_array[k] = CORRELATE(image_hr, (rec_image_hr_h/2. + rec_image_hr_p/2.))
    MSE_array[k] = calMSE(image_hr, (rec_image_hr_h/2. + rec_image_hr_p/2.))
  
  ENDFOR  
  print, 'correlation: ', correlation_array
  print, 'MSE: ', MSE_array
  ;
  iplot, k_array, correlation_array, title = 'Correlation Coefficient', xtitle='# of signal averages', ytitle='Coefficient'
  iplot, k_array, MSE_array, title = 'MSE', xtitle='# of signal averages', ytitle='MSE'
  ;

END

;@GJ, 2022/12/5, do image reconstruction based on orthogonal excitation
PRO orthDC_excitation_image_reconstruction

  directory = 'C:\D_drive\MPI_Tianjie\ACDC_Orth\synomag70_20221119\Low_vis\'
  sav_name = DIALOG_PICKFILE(PATH=directory, /MUST_EXIST, TITLE="sav file", FILTER='*.sav')
  RESTORE, sav_name
  
;  save, N_tp, signal_array_selected, FOV, gradient, H_DC_offset, H_DC_range

  N_steps = 200
  SM_hr = CONGRID(signal_array_selected, N_steps, N_tp*N_steps/100, /CENTER, /INTERP)
  iimage, SM_hr, title='system matrix based on orthogonal excitation'
  max_SM_hr = MAX(SM_hr, maxind)
  ind = ARRAY_INDICES(SM_hr, maxind)
  ;define the matrix
  SM = SM_hr[*, ind[1]-N_steps/2+1:ind[1]+N_steps/2]
  iimage, SM, title='system matrix'
  p_fn = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\concentration_bw.jpg'
  READ_JPEG, p_fn, phantom_image
;  ;@GJ, 2022/9/15, removing the tumor region
;  phantom_image[*, 320:*, 0:370] = 0.
   image_hr = CONGRID(REFORM(phantom_image[0,*,*]), N_steps, N_steps)
;  p_fn = 'C:\D_drive\MPI_Tianjie\ACDC_Orth\MIP_Unet\002.jpeg'
;  READ_JPEG, p_fn, phantom_image
;  image_hr = CONGRID(phantom_image, N_steps, N_steps)
  iimage, image_hr, title='original'
  
  T_space = image_hr * 0.
  rec_image = image_hr * 0.
  FOR i=0, N_steps-1 DO BEGIN
    FFL = REFORM(image_hr[*,i])
    signal = SM ## FFL
    T_space[i, *] = signal
    rec_image_1d = art_func(SM,signal)
    rec_image[*,i] = rec_image_1d
    print, 'i=', i
  ENDFOR
  iimage, rec_image, title='reconstructed'
  iimage, T_space, title='T space'
  
END

;GJ, 2022/1/27, do batch fitting of sine excitation
;GJ, 2022/1/30, change the name for MH curves
PRO batch_pulsed_MH_curves_fitting, sState


 ;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\syno_gly40_7_5mT.txt'
  ;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\peri_gly40_7_5mT.txt'
  ;  test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\peri_gly3_7_5mT.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\peri_gly1_7_5mT.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20220408\10mT\NE50_gly005.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20220408\10mT\NE50_gly010.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20220408\10mT\NE50_gly100.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20220408\10mT\NE25_gly100.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20220408\10mT\Synomag_gly100.txt'
;    test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\测量20220408\7_5mT\NE50_100.txt'
;    test_fit_para = read_pulsed_signal_dat(test_fn)
;    signalOnly = 1
;    test_fit_para = read_pulsed_signal_dat(test_fn)
;    test_fit_para = read_pulsed_signal_dat(test_fn, signalOnly=signalOnly)
  

;  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_fields\'
  ;dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\'
;  IF N_ELEMENTS(dir_name) EQ 0 THEN BEGIN
;    dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\'
;  ENDIF
  dir_name = DIALOG_PICKFILE(PATH=sState.directory, /MUST_EXIST, TITLE="txt files", /DIRECTORY)
  
  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  IF STRLEN(dir_name) EQ 0 THEN RETURN
  ;@GJ, 2022/7/2, update the directory
  sState.directory = dir_name
  
  ;GJ, 2022/2/6, redo the fitting and write the results
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  output_filename_temp = dir_name+current_time+'.xls'
  
  ;Select a file to save exported data, the default opened directory is the one you select just now
  output_filename=DIALOG_PICKFILE(FILE=output_filename_temp, FILTER = ['*.xls'], title='Save Your Output as EXCEL file', PATH=dir_name)
  ;// open an ASCII file for output
  ;IF STRCMP(STRMID(aSaveFile, 2, 3, /REVERSE_OFFSET), 'txt',/FOLD_CASE)   THEN BEGIN
  IF (STRLEN(output_filename) EQ 0)  THEN BEGIN
    infowarning = DIALOG_MESSAGE('No output excel file selected!', /ERROR)
    RETURN  
  ENDIF
  
  ;write file
  openw,unit,output_filename,/get_lun
  
  IF STRPOS(dir_name, '20221119') NE -1 THEN BEGIN
    print, 'this is an orthogonal excitation!'
    first_line =['filename', 'number', 'Hx AC [mT]', 'Hz DC [mT]', 'relaxation time [us]', 'delay phase angle [deg]', 'signal ratio', 'correlation coefficient']
  ENDIF ELSE BEGIN
    first_line =['filename', 'peak H', 'AUC_slew', 'AUC_hold', 'AUC_total', 'AUC_slew_ratio', 'size', 'tao_Debye_mono', 'tao_mono', 'tao_biex_1', 'tao_biex_2', 'tao_biex_1_comp', 'tao_biex_2_comp', 'tao_biex_12combined', 'correlation_debye', 'correlation_biex']
  ENDELSE
  printf, unit, FORMAT = '(16(A, %"\t"))', first_line
  
  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  ;@GJ, 2022/11/30, define the pointer about the parameters
  para_p_ptr = PTRARR(nrfile, /ALLOCATE_HEAP)
  signal_max_array = DBLARR(nrfile)*0.
  DC_field_array = DBLARR(nrfile)*0.
  IF STRPOS(dir_name, '20221119') NE -1 THEN BEGIN
    FOR i=0L, nrfile-1L DO BEGIN
      filename = FILE_BASENAME(a_file[i], '.txt')
      sState.filename = filename
      sState.directory = FILE_DIRNAME(a_file[i])
      WIDGET_CONTROL, sState.wFilenameField, SET_VALUE=sState.filename+', '+STRTRIM(STRING(i+1),1)+'/'+STRTRIM(STRING(nrfile),1)+' analyzing...'
      fit_para = read_pulsed_signal_dat(a_file[i], sState)
      *(para_p_ptr[i]) = fit_para
      signal_max_array[i] = MAX(ABS(fit_para.signal))
      DC_field_array[i] = fit_para.Hz0
      printf, unit, FORMAT = '(16(A, %"\t"))', filename, FIX(STRMID(filename, 1, 2, /REVERSE_OFFSET)), fit_para.Hx0, fit_para.Hz0, fit_para.relaxation_orth, fit_para.delay_phase_angle, fit_para.signal_ratio, fit_para.correlation_coeff
    ENDFOR
    FREE_LUN, unit
    
    items = STRARR(nrfile);
    psyms = [-15, -16, -17, -18, -19, 20, -21, -22, -23]
    colors = ['pink', 'green', 'red3', 'red4', 'blue', 'red5', 'red6', 'red7', 'red8']
    FOR i=0L, 4, 1L DO BEGIN
      fit_para = *(para_p_ptr[(SORT(DC_field_array))[i]])
      items[i] = STRTRIM(STRING(fit_para.Hz0),1)+' mT'
      IF i EQ 0L THEN BEGIN
        xtitle = 'time [ms]'
        ytitle = 'signal [A.U.]'
        title = 'signal curves'
        window, 6
        cgPlot, fit_para.time/1000., fit_para.filtered_signal, Color=colors[i MOD N_ELEMENTS(colors)], $; BACKGROUND='black', $; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], $
          Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, YRANGE=[-MAX(signal_max_array)*1.1, MAX(signal_max_array)*1.1]
         cgPlots, fit_para.time/1000., fit_para.fitted_signal, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2, Linestyle=2
        ENDIF ELSE BEGIN
        cgPlots, fit_para.time/1000., fit_para.filtered_signal, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2 ; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], Thick=1
        cgPlots, fit_para.time/1000., fit_para.fitted_signal, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2, Linestyle=2
      ENDELSE
    ENDFOR
    ; Location of legend in data coordinates.
    yloc = (!Y.CRange[1] - !Y.CRange[0]) * 0.95 + !Y.CRange[0]
    xloc = (!X.CRange[1] - !X.CRange[0]) * 0.05 + !X.CRange[0]
    ; Add the legend.
    IF i LT N_ELEMENTS(colors)-1 THEN BEGIN
      cgLegend, Title=items[i], Lines=lines, Color=colors[i mod N_ELEMENTS(colors)], Thick=1, $
        Location=[xloc,yloc], /Box, /Data
    ENDIF
    
    ;@GJ 2022/12/5, do the gradient encoding
    N_tp = N_ELEMENTS(fit_para.filtered_signal)
    gradient = 45.; mT/m
    FOV = 0.2; m
    H_DC_range = gradient * FOV ; mT
    H_DC_offset = 1.0 ; mT
    signal_array_selected = DBLARR(FLOOR(H_DC_range), N_tp)
    FOR i=0, FLOOR(H_DC_range)-1 DO BEGIN
      fit_para = *(para_p_ptr[(SORT(DC_field_array))[i]])
      IF fit_para.Hz0 EQ H_DC_offset+i THEN BEGIN
        signal_array_selected[i, *] = REFORM(fit_para.filtered_signal)
      ENDIF
    ENDFOR
    save, FILENAME = dir_name+current_time+'signal_SM_hr.sav', N_tp, signal_array_selected, FOV, gradient, H_DC_offset, H_DC_range
    
    RETURN
  ENDIF
  ;@GJ, 2022/11/30, the end of signal plot of orthogonal excitation
  
  
  H_max_array_p = DBLARR(nrfile)*0.
  H_max_array = DBLARR(nrfile)*0.
  M_max_array = DBLARR(nrfile)*0.
  FOR i=0L, nrfile-1L DO BEGIN
    IF i EQ nrfile THEN break
    filename = FILE_BASENAME(a_file[i], '.txt')
    sState.filename = filename
    sState.directory = FILE_DIRNAME(a_file[i])
    WIDGET_CONTROL, sState.wFilenameField, SET_VALUE=sState.filename+', '+STRTRIM(STRING(i+1),1)+'/'+STRTRIM(STRING(nrfile),1)+' analyzing...'
    fit_para = read_pulsed_signal_dat(a_file[i], sState)
    IF N_TAGS(fit_para) EQ 0 THEN BEGIN
      printf, unit, FORMAT = '(2(A, %"\t"))', filename, 'error!!!'
      
      ;@GJ, 2022/6/27, in case the dat file is wrongly read
      IF i EQ 0 THEN a_file = a_file[1:*]
      IF i NE 0 AND i LT nrfile-2L THEN a_file = [a_file[0:i-1], a_file[i+1:*]]
      i--
      nrfile--
      IF i EQ nrfile-1L THEN BEGIN
        a_file = a_file[0:nrfile-2L]
        break
      ENDIF
      continue
    ENDIF
    *(para_p_ptr[i]) = fit_para
    H_max_array_p[i] = MAX(fit_para.H)
    M_max_array[i] = MAX(fit_para.M)
    H_max = ROUND(fit_para.H_peak_p*10.)/10.
    H_max_array[i] = H_max
    IF STRPOS(dir_name, 'WangQianData') NE -1 THEN BEGIN
      printf, unit, FORMAT = '(16(A, %"\t"))', filename, H_max, fit_para.AUC_slew, fit_para.AUC_hold, fit_para.AUC_total, fit_para.AUC_slew_ratio, $
        fit_para.particle_size, fit_para.tao_debye_mono
    ENDIF ELSE BEGIN
      printf, unit, FORMAT = '(16(A, %"\t"))', filename, H_max, fit_para.AUC_slew, fit_para.AUC_hold, fit_para.AUC_total, fit_para.AUC_slew_ratio, $
        fit_para.particle_size_mono, fit_para.tao_debye_mono, fit_para.tao, fit_para.bi_tao, fit_para.bi_tao_comp, fit_para.bi_tao_combine, $
        fit_para.corre_biex, fit_para.corre_Debye_Mono
    ENDELSE   
  ENDFOR
  
  IF STRPOS(dir_name, 'WangQianData') NE -1 THEN BEGIN
    FREE_LUN, unit
    RETURN
  ENDIF
   
  ; Create the legend with NASA Astronomy routine AL_LEGEND.
;  items = ['0.01', '0.03', '0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5']
;  items = ['2.5 mT', '5 mT', '7.5 mT']
  items = STRARR(nrfile);
  psyms = [-15, -16, -17, -18, -19, 20, -21, -22, -23]
  colors = ['pink', 'green', 'red3', 'red4', 'blue', 'red5', 'red6', 'red7', 'red8']
  ;GJ, 2022/2/6, for langevin function and size estimation
  X = DBLARR(nrfile*2.)
  Y = DBLARR(nrfile*2.)
  
  ;GJ, 2022/2/18, plot Neel and Brownian vs H
  Neel_array = DBLARR(nrfile)*0.
  Brownian_array = DBLARR(nrfile)*0.
  signal_p_max_array = DBLARR(nrfile)*0.
  signal_p_min_array = DBLARR(nrfile)*0.
  
;  ;GJ, 2022/2/11, if there are multiple H values, then do particle size estimation. Else return
;  IF nrfile GT N_ELEMENTS(colors) THEN BEGIN
;    FREE_LUN, unit
;    RETURN
;  ENDIF
  
 
  FOR i=nrfile-1L, 0L, -1 DO BEGIN
    fit_para = *(para_p_ptr[i])
    ;GJ, 2022/2/7, get the avearge peak value H
    X[i*2] = fit_para.H_peak_n
    X[i*2+1] = fit_para.H_peak_p
    Y[i*2] = MIN(fit_para.M)
    Y[i*2+1] = MAX(fit_para.M)
    items[i] = STRTRIM(STRING(H_max_array[i], format=''), 1) + ' mT'
    
    ;GJ, 2022/2/18, plot Neel and Brownian vs H
    Neel_array[i] = fit_para.bi_tao[1]
    Brownian_array[i] = fit_para.bi_tao[0]
    signal_p_max_array[i] = MAX(fit_para.signal_p)
    signal_p_min_array[i] = MIN(fit_para.signal_p)
    
    IF i EQ nrfile-1L THEN BEGIN
      xtitle = 'H [mT/u0]'
      ytitle = 'M [A.U.]'
      filename_arr = STRSPLIT(filename, '_', /EXTRACT)
      IF N_ELEMENTS(filename_arr) GT 1 THEN BEGIN
        title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+' ' + filename_arr[1] +')';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
      ENDIF ELSE BEGIN
        title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+')';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
      ENDELSE
      window, 6
      cgPlot, fit_para.H, fit_para.M, Color=colors[i MOD N_ELEMENTS(colors)], $; BACKGROUND='black', $; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], $
        Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, XRANGE=[-MAX(H_max_array_p)*1.1, MAX(H_max_array_p)*1.1], YRANGE=[-MAX(M_max_array)*1.1, MAX(M_max_array)*1.1]
      ;GJ, 2022/2/9, plot the fitted M curve
;      cgPlots, fit_para.debye_H, fit_para.debye_M, Color=colors[i], Thick=2
;      cgPlots, fit_para.relax_H, fit_para.relax_mono_Mfit, Color='black', Thick=2, Linestyle=1
      cgPlots, fit_para.relax_H, fit_para.relax_biex_Mfit, Color='black', Thick=2, Linestyle=1
      cgPlots, fit_para.debye_H, fit_para.debye_mono_Mfit, Color='black', Thick=2, Linestyle=2
    ENDIF ELSE BEGIN
      cgPlots, fit_para.H, fit_para.M, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2 ; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], Thick=1
      ;GJ, 2022/2/9, plot the fitted M curve
;      cgPlots, fit_para.debye_H, fit_para.debye_M, Color=colors[i], Thick=2
;      cgPlots, fit_para.relax_H, fit_para.relax_mono_Mfit, Color='black', Thick=2, Linestyle=1
      cgPlots, fit_para.relax_H, fit_para.relax_biex_Mfit, Color='black', Thick=2, Linestyle=1
      cgPlots, fit_para.debye_H, fit_para.debye_mono_Mfit, Color='black', Thick=2, Linestyle=2
    ENDELSE

;    ; Location of legend in data coordinates.
;    yloc = (!Y.CRange[1] - !Y.CRange[0]) * 0.95 + !Y.CRange[0]
;    xloc = (!X.CRange[1] - !X.CRange[0]) * 0.05 + !X.CRange[0]
;
;;    ; Add the legend.
;    IF i LT N_ELEMENTS(colors)-1 THEN BEGIN
;      cgLegend, Title=items[i], Lines=lines, Color=colors[i mod N_ELEMENTS(colors)], Thick=1, $
;        Location=[xloc,yloc], /Box, /Data
;    ENDIF
  ENDFOR
  
;  ;@GJ, 2022/6/26, plot the decay signals
;  sort_ind = SORT(H_max_array)
;  FOR i=nrfile-1L, 0L, -1 DO BEGIN
;    fit_para = *(para_p_ptr[sort_ind[i]])
;
;    IF i EQ nrfile-1L THEN BEGIN
;      xtitle = 't [us]'
;      ytitle = 'signal [A.U.]'
;      filename_arr = STRSPLIT(filename, '_', /EXTRACT)
;      IF N_ELEMENTS(filename_arr) GT 1 THEN BEGIN
;        title = 'Field-flat Signal Curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+' ' + filename_arr[1]+')'; + STRING(particle_size,format='(f10.4)')+' nm)';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
;      ENDIF ELSE BEGIN
;        title = 'Field-flat Signal Curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+')';+ STRING(particle_size,format='(f10.4)')+' nm)';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
;      ENDELSE
;      window, 7
;      cgPlot, fit_para.time_flat, fit_para.signal_p, Color=colors[i MOD N_ELEMENTS(colors)], $; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], $
;        Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, XRANGE=[0., MAX(fit_para.time_flat)*1.1], YRANGE=[MIN(signal_p_min_array)*1.1, MAX(signal_p_max_array)*1.1]
;    ENDIF ELSE BEGIN
;      cgPlots, fit_para.time_flat, fit_para.signal_p, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2 ; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], Thick=1
;    ENDELSE
;    
;  ENDFOR
;  ; Location of legend in data coordinates.
;  yloc = (!Y.CRange[1] - !Y.CRange[0]) * 0.95 + !Y.CRange[0]
;  xloc = (!X.CRange[1] - !X.CRange[0]) * 0.65 + !X.CRange[0]
;  ;    ; Add the legend.
;  IF nrfile LT N_ELEMENTS(colors) THEN n_item = nrfile ELSE n_item = N_ELEMENTS(colors)
;  cgLegend, Title=items[0:n_item-1], Color=colors[0:n_item-1], Thick=1, $
;     Location=[xloc,yloc], /Box, /Data
  
  ;@GJ, 2022/7/2, the H max are the same then return
  IF MAX(H_max_array)-MIN(H_max_array) LT 0.3 THEN RETURN
  
  ;@GJ, 2022/6/29, a color plot of signals
  ; Scale the height for drawing zonal winds and load program colors.
  sort_ind = SORT(H_max_array)
  HColors = BytScl(H_max_array[sort_ind]*2., Min=MIN(H_max_array), Max=MAX(H_max_array)*2.+1, Top=MAX(H_max_array)*2.+1)
  ;@GJ, 2022/6/30, correct the color bar
  cgLoadCT, 33, NColors=MAX(H_max_array)*2.+1
  ; Open a display window.
  cgDisplay
  ; Set up plotting parameters.
  thick = (!D.Name EQ 'PS') ? 6 : 2
  plotPosition = [0.125, 0.125, 0.9, 0.75]
  ; Draw a gray background so colors stand out more.
  cgColorFill, Position=plotPosition, Color='gray'
  ; Draw the initial plot, no data.
  FOR i=nrfile-1L, 0L, -1 DO BEGIN
    fit_para = *(para_p_ptr[sort_ind[i]])
    IF i EQ nrfile-1L THEN BEGIN
      cgPlot, fit_para.time_flat, fit_para.signal_p, XRange=[0,200], YRange=[0,0.2],YStyle=1, $
        YTitle='Signal [A.U.]', XTitle='time [us]', XStyle=1, $
        XTicks=10, /NoData, /NoErase, $
        Position=plotPosition, Label='Signals in field flat phase'
    ENDIF
    IF N_ELEMENTS((fit_para.time_flat)) LT 201 THEN nele = N_ELEMENTS((fit_para.time_flat)) ELSE nele=201
    IF (H_max_array[sort_ind[i]]-MIN(H_max_array))*2. LT N_ELEMENTS(HColors) THEN BEGIN
      cgPlotS, (fit_para.time_flat)[0:nele-1], (fit_para.signal_p)[0:nele-1]+(H_max_array[sort_ind[i]]+0.5)*0.01, fit_para.time_flat, Color=HColors[(H_max_array[sort_ind[i]]-MIN(H_max_array))*2.], Thick=thick
      ;    print, 'H max = ', H_max_array[sort_ind[i]], ', signal min = ', (H_max_array[sort_ind[i]]+0.5)*0.01
    ENDIF ELSE BEGIN
      ;@GJ, 2023/4/13, handling this error
      FREE_LUN, unit
      RETURN
    ENDELSE
  ENDFOR
  ; Add a color bar.
  ;GJ, 2022/3/29, correct the spelling of glycerol
  ;GJ, 2022/6/30, correct the range for colorbar label
  cgColorbar, NColors=MAX(H_max_array)*2.+1, Range=[0,MAX(H_max_array)], Title='Trapezoid Wave Field Amplitude [mT]', $
    TLocation='top', Format='(i0)', Position=[0.125, 0.86, 0.9, 0.90]
  
  ;@GJ, 2022/7/26, plot the signal during rise and flat phases
  ; Open a display window.
  cgDisplay
  ; Set up plotting parameters.
  thick = (!D.Name EQ 'PS') ? 6 : 2
  plotPosition = [0.125, 0.125, 0.9, 0.75]
  ; Draw a gray background so colors stand out more.
  cgColorFill, Position=plotPosition, Color='gray'
  ; Draw the initial plot, no data.
  FOR i=nrfile-1L, 0L, -1 DO BEGIN
    fit_para = *(para_p_ptr[sort_ind[i]])
    IF i EQ nrfile-1L THEN BEGIN
      cgPlot, fit_para.time_rise_flat, fit_para.signal_rf_p, XRange=[0,200], YRange=[0,0.2],YStyle=1, $
        YTitle='signal [A.U.]', XTitle='time [us]', XStyle=1, $
        XTicks=10, /NoData, /NoErase, $
        Position=plotPosition, Label='Signals in field rise and flat phases'
    ENDIF
    IF N_ELEMENTS((fit_para.time_rise_flat)) LT 201 THEN nele = N_ELEMENTS((fit_para.time_rise_flat)) ELSE nele=201
    cgPlotS, (fit_para.time_rise_flat)[0:nele-1], (fit_para.signal_rf_p)[0:nele-1]+(H_max_array[sort_ind[i]]+0.5)*0.01, fit_para.time_rise_flat, Color=HColors[(H_max_array[sort_ind[i]]-MIN(H_max_array))*2.], Thick=thick
    ;    print, 'H max = ', H_max_array[sort_ind[i]], ', signal min = ', (H_max_array[sort_ind[i]]+0.5)*0.01
  ENDFOR
  ; Add a color bar.
  ;GJ, 2022/3/29, correct the spelling of glycerol
  ;GJ, 2022/6/30, correct the range for colorbar label
  cgColorbar, NColors=MAX(H_max_array)*2.+1, Range=[0,MAX(H_max_array)], Title='Trapezoid Wave Field Amplitude [mT]', $
    TLocation='top', Format='(i0)', Position=[0.125, 0.86, 0.9, 0.90]

  cgLoadCT, 0
  
  ;GJ, 2022/7/27, save the image
  signal_matrix = DBLARR(nrfile+1, N_ELEMENTS((fit_para.time_rise_flat)))
  signal_matrix[0, *] = fit_para.time_rise_flat
  signal_matrix_sP = DBLARR(nrfile, fit_para.times, N_ELEMENTS((fit_para.time_rise_flat)))
  FOR i=nrfile-1L, 0L, -1 DO BEGIN
    fit_para = *(para_p_ptr[sort_ind[i]])
    IF N_ELEMENTS((fit_para.signal_rf_p)) EQ N_ELEMENTS(signal_matrix[sort_ind[i]+1, *]) THEN BEGIN
      signal_matrix[sort_ind[i]+1, *] = fit_para.signal_rf_p
      signal_matrix_sP[sort_ind[i], *, *] = fit_para.signal_rf_p_sP
    ENDIF
  ENDFOR
  iimage, signal_matrix
  SAVE, FILENAME = dir_name+current_time+'signalNewSP.sav', signal_matrix, signal_matrix_sP
  
  ;GJ, 2022/2/11, if there are multiple H values, then do particle size estimation. Else return
  IF MAX(ABS(H_max_array_p))-MIN(ABS(H_max_array_p)) LT 0.5 THEN BEGIN
    FREE_LUN, unit
    RETURN
  ENDIF

  ;GJ, 2022/6/22, if there is only one particle
  info=dialog_message('MNP size estimation, continue?', /question)
  IF STRCMP(info, 'No') THEN BEGIN
    FREE_LUN, unit
    RETURN
  ENDIF
  
  ;GJ, 2022/2/18, plot Neel and Brownian vs H
  H_array_int = FINDGEN(21);in mT
  ; Interpolate.
  Neel_array_int = INTERPOL(Neel_array[sort_ind], H_max_array[sort_ind], H_array_int, /QUADRATIC)
  Brownian_array_int = INTERPOL(Brownian_array[sort_ind], H_max_array[sort_ind], H_array_int, /QUADRATIC)
  time_diff_array = ABS(Neel_array_int-Brownian_array_int)
  min_td = MIN(time_diff_array[5:*], min_td_ind)
  H_equal = H_array_int[min_td_ind+5]
  print, 'Brownian = Neel: H = ', H_equal, ' mT'
  window, 9
  xtitle1 = 'H [mT/u0]'
  ytitle1 = 'tao [us]'
  title1 = 'Relaxation times ('+file_basename(filename, 'txt')+', taoB = taoN at '+STRTRIM(STRING(H_equal, format=''), 1) +' mT)';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
  items1 = ['Brownian', 'Neel']
  psyms1 = [-15, -16]
  colors1 = ['red', 'green']
  Lines1 = [0, 2]
  brownian_inf = WHERE( ~FINITE(brownian_array_int), count_b_inf)
  Neel_inf = WHERE( ~FINITE(Neel_array_int), count_n_inf)
  IF count_b_inf EQ 0 AND count_n_inf EQ 0 THEN BEGIN
    cgPlot, H_array_int, Brownian_array_int, Title=title1, XTitle=xtitle1, YTitle=ytitle1, Thick=2, Linestyle=Lines1[1], Color=colors1[0], PSym=psyms1[0]
    cgPlots, H_array_int, Neel_array_int, Thick=2, Linestyle=Lines1[1], Color=colors1[1], PSym=psyms1[1]
    ; Location of legend in data coordinates.
    yloc = (!Y.CRange[1] - !Y.CRange[0]) * 0.95 + !Y.CRange[0]
    xloc = (!X.CRange[1] - !X.CRange[0]) * 0.65 + !X.CRange[0]
    ; Add the legend.
    cgLegend, Title=items1, Lines=lines1, Color=colors1, Thick=2, $;, PSym=psyms1, $
      Location=[xloc,yloc], /Box, /Data
    ;
  ENDIF

  

  
  ;GJ, 2022/2/6
  ;initial guess of the particle size
  ;@GJ, 2022/7/21, fit the MH from both large and small size
  particle_size_Large = 25.;19.;40.;50.;19.;40.;50;70.;19.;30.
  particle_size_Small = 4.;19.;40.;50.;19.;40.;50;70.;19.;30.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle_Large = Msat_T*(particle_size_Large^3) /24./1.380469/T_p_kelvin ; in 1/mT
  beta_particle_Small = Msat_T*(particle_size_Small^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle_Large = Msat_T*(particle_size_Large^3) /24./1.380469/309.65 ; in 1/mT
;  beta_particle_Small = Msat_T*(particle_size_Small^3) /24./1.380469/309.65 ; in 1/mT
  FWHM_particle_Large = 4.16/beta_particle_Large/3.5; mT/mm
  FWHM_particle_Small = 4.16/beta_particle_Small/3.5; mT/mm
  ;M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
  
  ;sort
  sort_ind = SORT(X)
  X = X[sort_ind]
  Y = Y[sort_ind]
  A_Large = [MAX(Y), beta_particle_Large]
  A_Small = [MAX(Y), beta_particle_Small]
  weights = Y * 0. + 1. 
  ;Compute the parameters.
  yfit_Large = CURVEFIT(X, Y, weights, A_Large, SIGMA, CHISQ=csq_Large, FUNCTION_NAME='langevin_funct')
  yfit_Small = CURVEFIT(X, Y, weights, A_Small, SIGMA, CHISQ=csq_Small, FUNCTION_NAME='langevin_funct')
  
  ;Print the parameters returned in A.
;  PRINT, 'Function parameters: ', A
  IF csq_Large LT csq_Small THEN BEGIN
    A = A_Large
    yfit = yfit_Large
  ENDIF ELSE BEGIN
    A = A_Small
    yfit = yfit_Small
  ENDELSE
  particle_size = (A[1]*24.*1.380469*309.65/Msat_T)^(1./3.)
  print, 'particle size (From Langevin Function): ', particle_size, ' nm'
  
  printf, unit, FORMAT = '(10(A, %"\t"))', 'particle size (nm)'
  printf, unit, FORMAT = '(10(A, %"\t"))', particle_size
  ;; close file and END execution
  FREE_LUN, unit
  
  window, 10
  n_ele_x = 100
  X_array = (FINDGEN(n_ele_x+1)/n_ele_x-0.5)*MAX(X)*2.
  X_array = [X_array[0]-1., X_array, X_array[n_ele_x]+1]
  yfit_array = A[0]*(1./TANH(A[1]*X_array) - 1./(A[1]*X_array))
  yfit_array[n_ele_x/2+1]=0.
  plot, X_array, yfit_array, XTITLE='H [mT/u0]', YTITLE='M [A.U.]', TITLE='M-H curve (' + filename_arr[0] + STRING(particle_size,format='(f10.4)')+' nm)'
  oplot, X, y, PSYM=1
  
  ;@GJ, 2022/6/23, plot the MH curve and Langevin
  FOR i=nrfile-1L, 0L, -1 DO BEGIN
    fit_para = *(para_p_ptr[i])
    ;GJ, 2022/2/7, get the avearge peak value H
    X[i*2] = fit_para.H_peak_n
    X[i*2+1] = fit_para.H_peak_p
    Y[i*2] = MIN(fit_para.M)
    Y[i*2+1] = MAX(fit_para.M)
    H_max = ROUND(fit_para.H_peak_p*10.)/10.
    items[i] = STRTRIM(STRING(H_max, format=''), 1) + ' mT'

    ;GJ, 2022/2/18, plot Neel and Brownian vs H
    Neel_array[i] = fit_para.bi_tao[1]
    Brownian_array[i] = fit_para.bi_tao[0]

    IF i EQ nrfile-1L THEN BEGIN
      xtitle = 'H [mT/u0]'
      ytitle = 'M [A.U.]'
      filename_arr = STRSPLIT(filename, '_', /EXTRACT)
      IF N_ELEMENTS(filename_arr) GT 1 THEN BEGIN
        title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+' ' + filename_arr[1] + STRING(particle_size,format='(f10.4)')+' nm)';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
      ENDIF ELSE BEGIN
        title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+ STRING(particle_size,format='(f10.4)')+' nm)';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
      ENDELSE
      window, 11
      cgPlot, fit_para.H, fit_para.M, Color=colors[i MOD N_ELEMENTS(colors)], $; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], $
        Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, XRANGE=[-MAX(H_max_array_p)*1.1, MAX(H_max_array_p)*1.1], YRANGE=[-MAX(M_max_array)*1.1, MAX(M_max_array)*1.1]
     ENDIF ELSE BEGIN
      cgPlots, fit_para.H, fit_para.M, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2 ; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], Thick=1
    ENDELSE
  ENDFOR
  cgPlots, X_array, yfit_array, Color='black', Thick=2 
  
END

;GJ, 2022/1/27, do batch fitting of sine excitation
;GJ, 2022/1/30, change the name for MH curves
PRO batch_sine_MH_curves_fitting


  ;   test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\500Hz-10mT-Synomag\202112216.txt'
  ;   test_fit_para = read_sine_signal_dat(test_fn)
  ;
  ;   test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\8000Hz-10mT-Synomag\202112215.txt'
  ;   test_fit_para = read_sine_signal_dat(test_fn)

  dir_name = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas_sinewave_Synomag\'
  dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)
  
  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  
  current_time = string(SYSTIME(/JULIAN, /UTC), FORMAT='(f8.0)')
  output_filename = dir_name+current_time+'xls'
;  output_filename = dir_name+'sine_para20220131.xls'
  
  
  ;write file
  openw,unit,output_filename,/get_lun
  first_line =['filename', 'tao', 'tao_1', 'tao_2', 'tao_1_comp', 'tao_2_comp']
  printf, unit, FORMAT = '(10(A, %"\t"))', first_line
  
  a = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  para_p_ptr = PTRARR(nrfile, /ALLOCATE_HEAP)
  FOR i=0L, nrfile-1L DO BEGIN
    filename = FILE_BASENAME(a[i], '.txt')
    fit_para = read_sine_signal_dat(a[i])
    printf, unit, FORMAT = '(10(A, %"\t"))', filename, fit_para.tao, fit_para.bi_tao, fit_para.bi_tao_comp
    *(para_p_ptr[i]) = fit_para
  ENDFOR
  
  ;; close file and END execution
  FREE_LUN, unit
  
  ; Create the legend with NASA Astronomy routine AL_LEGEND.
  items = ['0.01', '0.03', '0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5']
  psyms = [-15, -16, -17, -18, -19, 20, -21, -22, -23]
  colors = ['blue', 'green', 'pink', 'red3', 'red4', 'red5', 'red6', 'red7', 'red8']
  FOR i=0L, nrfile-1L DO BEGIN
    fit_para = *(para_p_ptr[i])
    IF i EQ 0 THEN BEGIN
      xtitle = 'H [mT/u0]'
      ytitle = 'M [A.U.]'
      H_max = ROUND(MAX(fit_para.M_H_array[0, *]))
      title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, ' + STRTRIM(STRING(H_max), 1) + ' mT)'
      cgPlot, fit_para.M_H_array[0,*], fit_para.M_H_array[1, *]/MAX(fit_para.M_H_array[1, *]), Color=colors[i], $; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], $
        Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, YRANGE=[-1.2, 1.2]
    ENDIF ELSE BEGIN
      cgPlots, fit_para.M_H_array[0,*], fit_para.M_H_array[1, *]/MAX(fit_para.M_H_array[1, *]), Color=colors[i], Thick=2 ; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], Thick=1
    ENDELSE
    
    ; Location of legend in data coordinates.
    yloc = (!Y.CRange[1] - !Y.CRange[0]) * 0.95 + !Y.CRange[0]
    xloc = (!X.CRange[1] - !X.CRange[0]) * 0.05 + !X.CRange[0]

    ; Add the legend.
    cgLegend, Title=items, Lines=lines, Color=colors, Thick=1, $
      Location=[xloc,yloc], /Box, /Data
  ENDFOR
  
  
  
  PTR_FREE, para_p_ptr
END

;GJ, 2021/12/24
;read data from Xin Feng
;sine excitation, relaxation induced signal decay
;Perimag, Synomag_D with different viscosity
;GJ, 2022/1/26, read sine signal
;GJ, 2022/1/26, automatically determine the frequency
;GJ, 2023/5/8, including Wang Qian's data about intracellular MNPs
;GJ, 2023/8/18, analyzing Li Lei's data from Nanobubble with DC offsets (30 mT + 24 mT), F:\MPI_Tianjie\Nanorubble_MPI\DC_offset
; fit_para = read_sine_signal_dat()
FUNCTION read_sine_signal_dat, filename

  IF N_ELEMENTS(filename) EQ 0 THEN filename = dialog_pickfile(TITLE='Select Signal dat File', FILTER='*.txt', /MUST_EXIST, PATH='C:\D_drive\MPI_Tianjie\NanoRod_MPI\20030602Relax\01_SineExc\')
  
  ;read the txt file 1st time to get the frequency and period
  OPENR, unit, filename, /GET_LUN
  i_temp = 0LL;
;  H_auc = 0;
  WHILE i_temp LT 10000 DO BEGIN
    READF, unit, a_temp, FORMAT='(%"%f")'
;    H_auc += a_temp;
;    IF ABS(H_auc) LT ABS(a_temp/100.) THEN break;
    IF i_temp EQ 0 THEN BEGIN
      H_array_temp = a_temp
    ENDIF ELSE BEGIN
      H_array_temp = [H_array_temp, [a_temp]]
    ENDELSE
    i_temp ++;
  ENDWHILE
  FREE_LUN, unit
  
  FOR N_period = 50, i_temp/5-1, 25 DO BEGIN
    temp = [H_array_temp[0], H_array_temp[N_period], H_array_temp[2*N_period], H_array_temp[3*N_period], H_array_temp[4*N_period]]
    IF STDDEV(temp) LT 0.05 THEN break;
  ENDFOR
  ;N_period = frequency
  frequency = 1000000/N_period/100*100
  print, 'frequency = ', frequency, ' Hz'
  

  ;read the txt file 2nd time for calculation
  OPENR, unit, filename, /GET_LUN
  delta_t = 1.0
  i_temp = 0LL;
  ;get the period
;  frequency = 8000;500;
  N_period = 1000000/frequency
  i = FINDGEN(N_period)
  H_array = DBLARR(N_period)
  signal_array = DBLARR(N_period)
  signal_array_high_error = DBLARR(N_period)
  signal_array_low_error = DBLARR(N_period)
  spec_array = DBLARR(N_period)
  times = 1LL;400LL;400LL;400LL;20LL;400LL
  WHILE ~ EOF(unit) DO BEGIN
    IF STRPOS(filename, '20221119') NE -1 THEN BEGIN
      READF, unit, a_temp, noise_temp, b_temp, H_temp, M_temp, FORMAT='(%"%f %f %f %f %f")'
      c_temp = b_temp
    ENDIF ELSE BEGIN
      IF STRPOS(filename, 'WangQianData') NE -1 OR STRPOS(filename, '202306') NE -1 OR STRPOS(filename, '20230817') NE -1 THEN BEGIN
        IF STRPOS(filename, '202306') NE -1 OR STRPOS(filename, '20230817') NE -1 THEN BEGIN
          READF, unit, a_temp, sig_temp, b_temp, c_temp, H_temp, M_temp, FORMAT='(%"%f %f %f %f %f %f")'
          ;excite_temp = b_temp
        ENDIF ELSE BEGIN
          READF, unit, a_temp, b_temp, excite_temp, c_temp, H_temp, M_temp, FORMAT='(%"%f %f %f %f %f %f")'
          ;excite_temp = b_temp
        ENDELSE
        
      ENDIF ELSE BEGIN
        READF, unit, a_temp, b_temp, c_temp, H_temp, M_temp, FORMAT='(%"%f %f %f %f %f")'
      ENDELSE
    ENDELSE
    ;IF ABS(c_temp) LT 0.001 THEN BREAK
    
    IF i_temp LT N_period THEN BEGIN
      signal_array_high_error[i_temp MOD N_period] = b_temp
      signal_array_low_error[i_temp MOD N_period] = b_temp
    ENDIF ELSE BEGIN
      signal_array_high_error[i_temp MOD N_period] = MAX([signal_array_high_error[i_temp MOD N_period], b_temp])
      signal_array_low_error[i_temp MOD N_period] = MIN([signal_array_low_error[i_temp MOD N_period], b_temp])
    ENDELSE
    H_array[i_temp MOD N_period] += a_temp
    signal_array[i_temp MOD N_period] += b_temp
    spec_array[i_temp MOD N_period] += c_temp
    
    ;save the M-H curve
    IF M_temp NE 0 THEN BEGIN
      IF i_temp EQ 0 THEN BEGIN
        M_H_array = TRANSPOSE([H_temp, M_temp])
      ENDIF ELSE BEGIN
        M_H_array = [M_H_array, TRANSPOSE([H_temp, M_temp])]
      ENDELSE
    ENDIF
    
    i_temp++
    IF i_temp GT N_period*times THEN BREAK
    ;save the signal as previous one
    ;      previous_signal = b_temp
    ;      print, 'index = ', i_temp
  ENDWHILE
  FREE_LUN, unit
  
  ;transpose M_H_array as [2, ***] matrix
  M_H_array = TRANSPOSE(M_H_array)
  ;plot, M_H curve
;  window, 0
  xtitle = 'H [mT/u0]'
  ytitle = 'M [A/s]'
  title = 'M-H curve ('+STRTRIM(STRING(frequency),1)+' Hz)'
  IF STRPOS(filename, '202306') NE -1 THEN BEGIN
    cgPlot, M_H_array[0,*], -M_H_array[1, *], Color='red', PSym=-16,  SymSize=0.5, SymColor='red', $
      Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, /window

  ENDIF ELSE BEGIN
    cgPlot, M_H_array[0,*], M_H_array[1, *], Color='red', PSym=-16,  SymSize=0.5, SymColor='red', $
      Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, /window

  ENDELSE
  
;  iplot, M_H_array[0,*], M_H_array[1, *], Color='red', xrange=[-40.,40.], yrange=[-1.e-5, 1.e-5], Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2
  
   
  ;make sure there is not enough data for times
  times = i_temp/N_period
  H_array /= (times*1.0)
  signal_array /= (times*1.0)
  spec_array /= (times*1.0)

  ;find the flat portion index
  H_array_diff = ABS(H_array[1:*]-H_array[0:N_period-1])
  mm=where(H_array_diff LT 0.2, count)
  ;15, 238
  ;265,488
  signal_array_n = signal_array[0:N_period/2-1]
  signal_array_p = signal_array[N_period/2:N_period-1]

;  window, 1
  ; Set up variables for the plot. Normally, these values would be
  ; passed into the program as positional and keyword parameters.
  xtitle = 'time [us]'
  ytitle = 'A.U.'
  title = 'H and Signal ('+STRTRIM(STRING(frequency),1)+' Hz)'
  position = [0.125, 0.125, 0.9, 0.925]
  time = i * delta_t
  
  ;@GJ, 2023/8/18, to calculate M based on the signal's integral
  max_H_array = MAX(H_array, max_H_ind)
  time_s = SHIFT(time, -N_ELEMENTS(time)/4-max_H_ind)
  signal_array_s = SHIFT(signal_array, N_ELEMENTS(time)/4-max_H_ind)
  H_array_s = SHIFT(H_array, N_ELEMENTS(time)/4-max_H_ind)
  M_array_temp = TOTAL(-signal_array_s, /cumulative)*delta_t
  M_array_s = M_array_temp - (MAX(M_array_temp)-MIN(M_array_temp))/2.
  iplot, time, signal_array_s, yrange=[-0.2, 0.2], xtitle=xtitle, color='blue', thick=2, ytitle=ytitle, title=title
  iplot, time, H_array_s, xtitle=xtitle, color='green', thick=2, ytitle=ytitle, title=title
  iplot, -H_array_s, -M_array_s, xtitle='H [mT]', ytitle='M', title='M-H curve', color='red', thick=3, xrange=[-60.,60.], yrange=[-10.,10.]
  
  ; Set up a "window" for the plot. The PostScript output will have
  ; the same aspect ratio as the graphics window on the display.
  cgDisplay, 600, 500, Title='H and Signal '+STRTRIM(STRING(frequency),1)+' Hz'
  ; Draw the line plot with no data
  cgPlot, time, signal_array, Title=title, XTitle=xtitle, YTitle=ytitle, $
    Position=position, /NoData, YRange=[MIN(signal_array_low_error), MAX(signal_array_high_error)], YStyle=1
  ; Fill in the error estimates.
  cgColorFill, [time, Reverse(time), time[0]], $
    [signal_array_high_error, Reverse(signal_array_low_error), signal_array_high_error[0]], $
    Color='sky blue'
  ; Draw the line plot with no data
  cgPlotS, time, -H_array/MEAN(ABS(H_array))*MEAN(ABS(signal_array)), Color='blue', $
    Thick=2, LINESTYLE = 2

   ;@GJ, 2023/6/5

  IF STRPOS(filename, '202306') NE -1 OR STRPOS(filename, '20221020') NE -1 OR STRPOS(filename, '20230817') NE -1 THEN RETURN, M_H_array

   
  ;  cgplot, -H_array/MAX(H_array)*MAX(signal_array), LINESTYLE = 2, thick = 2, title='signal', xtitle='time [us]', ytitle='A.U.', BACKGROUND = 'FFFFFF'x, COLOR = 0
  ;  cgplot, signal_array, thick = 2, color='red',/overplot
  
  ;calculate the adiabatic signal
  Max_H = MAX(H_array, max_ind)
  shift_step = N_period/2-max_ind-1
  H_array_temp = SHIFT(H_array, shift_step)
  M_adiabatic_a = INTERPOL(M_H_array[1,0:N_period/2-1], M_H_array[0,0:N_period/2-1], H_array_temp[0:N_period/2-1])
  M_adiabatic_b = INTERPOL(M_H_array[1,N_period/2:*], M_H_array[0,N_period/2:*], H_array_temp[N_period/2:*])
  M_adiabatic = SHIFT([M_adiabatic_b, M_adiabatic_a], -shift_step)
  
  S_adiabatic = M_adiabatic[1:*] - M_adiabatic[0:N_period-1]
  cgPlotS, time[0:N_period-2], S_adiabatic/MEAN(ABS(S_adiabatic))*MEAN(ABS(signal_array)), Color='green', PSym=-16,  SymSize=0.5, SymColor='green', $
    Thick=2;, LINESTYLE = 3
  cgPlotS, time, signal_array, Color='red', PSym=-16,  SymSize=0.5, SymColor='red', $
    Thick=1, LINESTYLE = 3
  
  ;iplot, H_array, title='Pulsed Excitation'
  ;iplot, signal_array, title='Signal'

  ; Use Count to get the number of nonzero elements:
  ;  index = WHERE(H_array GT MAX(H_array)/2., count)
  ; Only subscript the array if it is safe:
  ;  IF count NE 0 THEN end_ind = MAX(index)
  
 
  ;GJ, 2022/1/30, calculate Debye and dual_Debye model
  Nelem = N_ELEMENTS(H_array)
  particle_size = 27.;19.;40.;50.;19.;40.;50;70.;19.;30.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  FWHM_particle = 4.16/beta_particle/3.5; mT/mm
  M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
  Max_H = MAX(H_array, max_ind_H)
  shift_step = N_period/2-max_ind_H
  signal_langevin_temp = -(M[1:*]-M[0:Nelem-1])/delta_t; kA/m/us
  signal_meas = signal_array[0:Nelem-2]
  time_half = time[0:Nelem-2]
  signal_langevin = signal_langevin_temp/MEAN(ABS(signal_langevin_temp))*MEAN(ABS(signal_meas))
  ;  signal_langevin = SHIFT(signal_langevin, -shift_step)
  window, 0
  plot, time_half, signal_langevin, title='Langevin & Signal, '+STRTRIM(STRING(frequency),1)+' Hz', XTITLE='time [us]', YTITLE='signal [a.u.]'
  oplot, time_half, signal_meas, line=2
  oplot, time_half, ABS(H_array)
  signal_meas[0:Nelem/2-1] = MEAN(signal_meas)
  signal_langevin[0:Nelem/2-1] = MEAN(signal_langevin)
  
  ;GJ, 2022/1/30, calculate Debye model
  relaxation_time = FINDGEN(100)+1.
  correlation_array = DBLARR(100)*0.
  MSE_array = DBLARR(100)*0.+10000.
  FOR i = 0, 99 DO BEGIN
    signal_nat = non_adiabatic(time_half, signal_langevin, relaxation_time[i])
    signal_nat = signal_nat/MEAN(signal_nat)*MEAN(signal_meas)
    correlation_array[i] = CORRELATE(signal_nat, signal_meas)
    MSE_array[i] = calMSE(signal_nat, signal_meas)
  ENDFOR

  ;max_corr = MAX(correlation_array, max_ind)
  min_MSE = MIN(MSE_array, max_ind)
  print, 'relaxation time = ', relaxation_time[max_ind]
  signal_nat = non_adiabatic(time_half, signal_langevin, relaxation_time[max_ind])
  signal_nat = signal_nat/MEAN(signal_nat)*MEAN(signal_meas)
  ;  oplot, time, signal_nat, line=3
  corre_Mono = CORRELATE(signal_meas, signal_nat)
  PRINT, 'correlation coeff (Mono): ', corre_Mono

  window, 1
  cgplot, time_half, signal_langevin, LINESTYLE = 2, thick = 2, title='Debye Mono-ex Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='signal [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  cgplot, time_half, signal_meas, thick = 2, LINESTYLE = 0, color='blue',/overplot
  cgplot, time_half, signal_nat, thick = 2, color='red',/overplot
  
  ;GJ, 2022/1/30, calculate dual_Debye model
  ;calculate the difference
  max_signal = max(signal_meas, max_ind_s)
  max_H = max(H_array[Nelem/2:*], max_ind_H)
  time_diff = 5;ABS(max_ind_s - max_ind_H - Nelem/4.)
  print, 'time diff = ', time_diff, ' us'
  N_Neel = ROUND(time_diff*10.)*2
  Neel_array = (-FINDGEN(N_Neel)/10.+N_Neel/20.)*delta_t*2.
  N_Brownian = ROUND(100)
  Brownian_array = (FINDGEN(N_Brownian))*delta_t*2.
  Neel_ratio = FINDGEN(101)/100.
  correlation_array = DBLARR(N_Neel, N_Brownian, 101)*0.
  MSE_array = DBLARR(N_Neel, N_Brownian, 101)*0.+10000.
  FOR i = 0, N_Neel-1 DO BEGIN
    signal_Neel = non_adiabatic(time_half, signal_langevin, Neel_array[i])
    FOR j = 0, N_Brownian-1 DO BEGIN
      signal_Brownian = non_adiabatic(time_half, signal_langevin, Brownian_array[j])
      FOR k = 0, 100 DO BEGIN
        signal_nat = Neel_ratio[k]*signal_Neel + (1.-Neel_ratio[k]) * signal_Brownian
        signal_nat = signal_nat/MEAN(signal_nat)*MEAN(signal_meas)
        correlation_array[i, j, k] = CORRELATE(signal_nat, signal_meas)
        MSE_array[i, j, k] = calMSE(signal_nat, signal_meas)
      ENDFOR
    ENDFOR
  ENDFOR

  ;  max_corr = MAX(correlation_array, max_ind_corr)
  min_MSE = MIN(MSE_array, max_ind_corr)
  max_index = ARRAY_INDICES(correlation_array, max_ind_corr)
  print, '(Dual Debye) Neel time = ', Neel_array[max_index[0]]
  print, '(Dual Debye) Brownian time = ', Brownian_array[max_index[1]]
  print, '(Dual Debye) Neel ratio = ', Neel_ratio[max_index[2]]

  signal_Neel = non_adiabatic(time_half, signal_langevin, Neel_array[max_index[0]])
  signal_Brownian = non_adiabatic(time_half, signal_langevin, Brownian_array[max_index[1]])
  signal_nat = Neel_ratio[max_index[2]]*signal_Neel + (1.-Neel_ratio[max_index[2]]) * signal_Brownian
  ;  oplot, time_half, signal_nat, line=3
  corre_Biex = CORRELATE(signal_meas, signal_nat)
  PRINT, 'correlation coeff (Biex): ', corre_Biex

  window, 2
  plot, time_half, signal_meas, title='Dual Debye Relaxation Fitting, '+STRTRIM(STRING(frequency),1)+' Hz'
  oplot, time_half, Neel_ratio[max_index[2]]*signal_Neel, LINESTYLE = 1
  oplot, time_half, (1.-Neel_ratio[max_index[2]]) * signal_Brownian, LINESTYLE = 2
  oplot, time_half, signal_nat, LINESTYLE = 3

  window, 3
  cgplot, time_half, signal_langevin, LINESTYLE = 2, thick = 2, title='Dual Debye Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='signal [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  cgplot, time_half, signal_meas, thick = 2, LINESTYLE = 0, color='blue',/overplot
  cgplot, time_half, signal_nat, thick = 2, color='red',/overplot
  
  fit_result = {$
    frequency:frequency, $
    time:time, $
    H:H_array, $
    signal:signal_array, $
    signal_p:signal_array_p, $ ;positive signal
    signal_n:signal_array_n, $ ;negative signal
    delta_t:delta_t, $
    M_adiabatic: M_adiabatic, $
    S_adiabatic: S_adiabatic, $
    M_H_array: M_H_array, $
    tao: relaxation_time[max_ind], $
    bi_tao: [Neel_array[max_index[0]], Brownian_array[max_index[1]]], $
    bi_tao_comp: [Neel_ratio[max_index[2]], 1.-Neel_ratio[max_index[2]]] $
  }
  RETURN, fit_result;[[H_array], [signal_array], [spec_array]]
END


;GJ, 2021/12/24
;read data from Xin Feng
;pulsed excitation, relaxation induced signal decay
;Perimag, Synomag_D with different viscosity
;calculate AUC, Neel relaxation and Brownian relaxation constants based on mono-exponential and bi-exponential model
;GJ, 2022/1/3, return the structure including all
;GJ, 2022/1/9, calculate AUC_p and AUC_n for positive and negative signals
;GJ, 2022/1/25, pack the result as a structure and fix the bug of bi-expo fitting
;GJ, 2022/1/25, add high_error and low_error to the signal plot
;GJ, 2022/1/26, automatically determine the frequency
;GJ, 2022/2/4, add a keyword 'signalOnly' for signal decay curve fitting, if there is no keyword, M is used for curve fitting
;GJ, 2022/6/16, add plot of magnetization
;GJ, 2022/11/21, add the orthogonal data file
;GJ, 2022/11/21, 
;test_fn = 'C:\D_drive\MPI_Tianjie\ACDC_Orth\synomag70_20221119\2022081016.txt'
;test_fit_para = read_pulsed_signal_dat(test_fn)
;example:
;test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\syno_gly40_7_5mT.txt'
; test_fit_para = read_pulsed_signal_dat(test_fn)
; 
; signalOnly=1
;test_fit_para = read_pulsed_signal_dat(test_fn, signalOnly=signalOnly)
;test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\AWR_meas\syno_gly1\syno_gly1_7_5mT.txt'
;test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220523\2022052321.txt'
;test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220524_T_5mT_90flat\20220524\2022052418.txt'
;test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220523_Gelatin_90flat\20220523\2022052325.txt'
;test_fn = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\Synomag_7.5mT\20220616121.txt'
; test_fit_para = read_pulsed_signal_dat(test_fn)

FUNCTION read_pulsed_signal_dat, filename, sState, particle_size, signalOnly=signalOnly
  ;GJ, 2022/6/19, clean the plot in sState
  ;GJ, 2022/11/30, remove the cleaning
;  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
;    new_pic = BYTARR(3, sState.xdim, sState.ydim)
;    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
;    WSET, drawpicid
;    TVSCL, new_pic, TRUE=1
;  ENDIF
  
  ;GJ, 2022/3/8, nonzero gilename
  IF (N_ELEMENTS(filename) EQ 0) THEN BEGIN
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      filename = DIALOG_PICKFILE(PATH=sState.directory, /MUST_EXIST, FILTER='*.txt', TITLE="txt file")
    ENDIF ELSE BEGIN
      filename = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\Xinfeng\20220616', /MUST_EXIST, FILTER='*.txt', TITLE="txt file")
    ENDELSE
    IF (STRLEN(filename) EQ 0) THEN BEGIN
      info=dialog_message('No file name!', /information)
      RETURN, 0
    ENDIF
  ENDIF
  IF STRLEN(filename) GT 1 THEN BEGIN
    filename_info = FILE_INFO(filename)
    IF (filename_info.exists EQ 0) THEN BEGIN
      info=dialog_message('Wrong file name!', /information)
      RETURN, 0
    ENDIF
    IF (STRLEN(filename) EQ 0) THEN BEGIN
      info=dialog_message('No file name!', /information)
      RETURN, 0
    ENDIF
  ENDIF
  
  ;replace the sstate.directory
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    sState.directory = FILE_DIRNAME(filename)
    sState.filename = FILE_BASENAME(filename)
  ENDIF
  
  ;@GJ, 2023/4/13, checking file size
  finfo_Result = FILE_INFO(filename)
  IF finfo_Result.size LT 500000. THEN RETURN, 0
  
  ;check whether the file is readable
  OPENR, unit, filename, /GET_LUN
  test_temp=''
  READF, unit, test_temp
  IF STRPOS(test_temp, '/') NE -1 THEN BEGIN
    FREE_LUN, unit
    info=dialog_message(file_basename(filename)+' is not readable!', /information)
    RETURN, 0
  ENDIF ELSE BEGIN
    FREE_LUN, unit
  ENDELSE


  ;read the txt file 1st time to get the frequency and period
  OPENR, unit, filename, /GET_LUN
  i_temp = 0LL;
  ;  H_auc = 0;
  WHILE (i_temp LT 10000) AND ~EOF(unit) DO BEGIN
    READF, unit, a_temp, FORMAT='(%"%f")'
    ;    H_auc += a_temp;
    ;    IF ABS(H_auc) LT ABS(a_temp/100.) THEN break;
    IF i_temp EQ 0 THEN BEGIN
      H_array_temp = a_temp
    ENDIF ELSE BEGIN
      H_array_temp = [H_array_temp, [a_temp]]
    ENDELSE
    i_temp ++;
  ENDWHILE
  FREE_LUN, unit

  ;GJ, 2022/6/19, fix the condition of searching frequency
  MAX_H=MAX(H_array_temp[0:500], max_ind_H)
  H_array_accu = TOTAL(H_array_temp, /cumulative)
  FOR N_period = 50, i_temp/5-1, 25 DO BEGIN
    temp = [H_array_temp[max_ind_H], H_array_temp[N_period+max_ind_H], H_array_temp[2*N_period+max_ind_H], H_array_temp[3*N_period+max_ind_H]]
    temp_accu = [H_array_accu[max_ind_H], H_array_accu[N_period+max_ind_H], H_array_accu[2*N_period+max_ind_H], H_array_accu[3*N_period+max_ind_H]]
    IF (STDDEV(temp) LT ABS(MAX(H_array_temp))*0.1) AND (STDDEV(temp_accu) LT ABS(MAX(temp_accu))*0.3) THEN break;
  ENDFOR
  ;@GJ, 2023/4/10, increasing 
  N_period = CEIL(N_period/100.)*100.
  IF STRPOS(filename, '20230602') NE -1 THEN N_period = 500.
  frequency = 1000000/N_period/100*100
  ; frequency = 2000
  print, 'frequency = ', frequency, ' Hz'

  ;read the txt file 2nd time for calculation
  OPENR, unit, filename, /GET_LUN
  delta_t = 1.0 ;us
  i_temp = 0LL;
  N_period = 1000000/frequency
  i = FINDGEN(N_period)
  H_array = DBLARR(N_period)
  signal_array = DBLARR(N_period)
  signal_array_high_error = DBLARR(N_period)
  signal_array_low_error = DBLARR(N_period)
  ;GJ, 2022/5/9, include high and low error
  spec_array = DBLARR(N_period)
  spec_array_high_error = DBLARR(N_period)
  spec_array_low_error = DBLARR(N_period)
  M_H_array = DBLARR(2, N_period)
  times = 400LL;400LL;400LL;20LL;400LL
  ;define the signal for every single period
  signal_array_sP = DBLARR(times, N_period)
  WHILE ~ EOF(unit) DO BEGIN
    ;      READF, unit, a_temp, b_temp, c_temp, FORMAT='(%"%f %f %f")'
    temp=''
    READF, unit, temp, FORMAT='(%"%s")'
    temp_spl = STRSPLIT(temp, STRING(9B), count=n_spl, /EXTRACT)
    IF n_spl GE 2 THEN BEGIN
      a_temp = temp_spl[0] ;field strength
      b_temp = temp_spl[1] ;signal
      IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, '20220720') NE -1 OR STRPOS(filename, '20220721') NE -1 OR STRPOS(filename, '20220815') NE -1 OR STRPOS(filename, '202304') NE -1  OR STRPOS(filename, '20230602') NE -1 THEN b_temp = temp_spl[2]
      IF n_spl EQ 3 THEN BEGIN
        c_temp = temp_spl[2] ;spec after FT
;        IF STRPOS(filename_base, '20220720') NE -1 OR STRPOS(filename_base, '20220721') NE -1 THEN b_temp = c_temp
        ;spec_array[i_temp MOD N_period] += c_temp
      ENDIF
      IF n_spl EQ 4 THEN BEGIN
        M_H_array[0, i_temp MOD N_period] = temp_spl[2]
        M_H_array[1, i_temp MOD N_period] = temp_spl[3]
      ENDIF
      ;@GJ, 2023/5/8, to process Wang QIan's Intracellualr data MPN
      IF n_spl EQ 6 THEN BEGIN
        c_temp = temp_spl[3] ;spec after FT
        IF ABS(temp_spl[4]) GT 0 OR ABS(temp_spl[5]) GT 0 THEN BEGIN
          M_H_array[0, i_temp MOD N_period] = temp_spl[4]
          M_H_array[1, i_temp MOD N_period] = temp_spl[5]
        ENDIF
      ENDIF
    ENDIF
    IF i_temp GT N_period*times THEN BREAK
    IF i_temp LT N_period THEN BEGIN
      signal_array_high_error[i_temp MOD N_period] = b_temp
      signal_array_low_error[i_temp MOD N_period] = b_temp
      ;include the data reading on 2022/5/8
      IF n_spl EQ 3 THEN BEGIN
        spec_array_high_error[i_temp MOD N_period] = c_temp
        spec_array_low_error[i_temp MOD N_period] = c_temp
      ENDIF
    ENDIF ELSE BEGIN
      signal_array_high_error[i_temp MOD N_period] = MAX([signal_array_high_error[i_temp MOD N_period], b_temp])
      signal_array_low_error[i_temp MOD N_period] = MIN([signal_array_low_error[i_temp MOD N_period], b_temp])
      ;include the data reading on 2022/5/8
      IF n_spl EQ 3 THEN BEGIN
        spec_array_high_error[i_temp MOD N_period] = MAX([spec_array_high_error[i_temp MOD N_period], c_temp])
        spec_array_low_error[i_temp MOD N_period] = MIN([spec_array_low_error[i_temp MOD N_period], c_temp])
      ENDIF
    ENDELSE

    ;check the readout
    ;      print, 'H '+STRING(i_temp)+' value = ', a_temp

    H_array[i_temp MOD N_period] += a_temp
    signal_array[i_temp MOD N_period] += b_temp
    ;@GJ, 2022/8/4, save all signal for every peroid
    IF i_temp/N_period LT times THEN signal_array_sP[i_temp/N_period, i_temp MOD N_period] = b_temp
    ;include the data reading on 2022/5/8
    IF n_spl EQ 3 THEN BEGIN
      spec_array[i_temp MOD N_period] += c_temp
    ENDIF

    i_temp++
    ;save the signal as previous one
    ;      previous_signal = b_temp
    ;       print, 'index = ', i_temp
    ;@GJ, 2022/11/30, check the signal is right
;    IF i_temp mod N_period EQ 0 THEN BEGIN
;      IF i_temp/N_period EQ 0 THEN BEGIN
;        window, 1
;        plot, signal_array/MAX(ABS(signal_array))
;      ENDIF ELSE BEGIN
;        oplot, signal_array/MAX(ABS(signal_array))
;      ENDELSE
;    ENDIF
  ENDWHILE
  FREE_LUN, unit

  ;make sure there is not enough data for times
  times = i_temp/N_period

  H_array /= (times*1.0)
  signal_array /= (times*1.0)
  ;  signal_array_high_error /= (times*1.0)
  ;  signal_array_low_error /= (times*1.0)
  spec_array /= (times*1.0)
  ;  spec_array_high_error /= (times*1.0)
  ;  spec_array_low_error /= (times*1.0)
  time = i * delta_t

  ;extract filename and ratio for calculation
  filename_base = file_basename(filename, '.txt')
  IF STRPOS(filename_base, '20220616') NE -1 THEN BEGIN
    IF STRPOS(filename_base, '_') NE -1 THEN BEGIN
      ratio_array=[200., 200., 100., 100., 100., 50., 50., 50., 50., 50., 50., 50., 50., 20., 20., 20., 20., 20., 20., 20.]
      filename_num = STRMID(filename_base, 0, 3)
      ratio = ratio_array[(FIX(filename_num)-1)]
    ENDIF ELSE BEGIN
      ratio_array=[100., 100., 100., 50., 50., 50., 50., 50., 50., 50., 50., 20., 20., 20., 20., 20., 20., 20., 200., 200.]
      filename_num = STRMID(filename_base, 8)
      ratio = ratio_array[(FIX(filename_num)-1)/10]
    ENDELSE
  ENDIF
  
  ;@GJ, 2022/7/21, repeated measurement of temperature  
  IF STRPOS(filename_base, '20220720') NE -1 OR STRPOS(filename_base, '20220721') NE -1 THEN BEGIN
    ratio_array=[500., 200., 200., 200., 100., 100., 100., 100., 100., 50., 50., 50., 50., 50., 50., 50., 20., 20., 20., 20.]
    IF STRPOS(filename_base, '_') NE -1 THEN BEGIN
      filename_num = STRMID(filename_base, 0, 3)
      ratio = ratio_array[(FIX(filename_num)-1)]
    ENDIF ELSE BEGIN
      filename_num = STRMID(filename_base, 8)
      ratio = ratio_array[(FIX(filename_num)-1)/10]
    ENDELSE
  ENDIF
  
  ;@GJ, 2022/8/15, repeated measurement of Gelatin
  IF STRPOS(filename_base, '20220815') NE -1 THEN BEGIN
    ratio_array=[500., 200., 200., 200., 100., 100., 100., 100., 100., 50., 50., 50., 50., 50., 50., 50., 20., 20., 20., 20.]
    IF STRPOS(filename_base, '_') NE -1 THEN BEGIN
      filename_num = STRMID(filename_base, 0, 3)
      ratio = ratio_array[(FIX(filename_num)-1)]
    ENDIF ELSE BEGIN
      filename_num = STRMID(filename_base, 8)
      ratio = ratio_array[(FIX(filename_num)-1) MOD 20]
    ENDELSE
  ENDIF
  
  IF STRPOS(filename_base, '20220523') NE -1 THEN BEGIN
    ratio_array=[100., 100., 100., 50., 50., 50., 50., 50., 50., 50., 20., 20., 20.]
    IF STRPOS(filename_base, '_') NE -1 THEN BEGIN
      filename_num = STRMID(filename_base, 0, 2)
      ratio = ratio_array[(FIX(filename_num)-1)]
    ENDIF ELSE BEGIN
      filename_num = STRMID(filename_base, 8)
      ratio = ratio_array[(FIX(filename_num)-1)/4]
    ENDELSE
  ENDIF
  
  ;@GJ, 2022/7/28, do the calculation
  IF N_ELEMENTS(ratio_array) GT 2 THEN BEGIN
    print, 'filename = ', filename
    print, 'ratio = ', ratio
    signal_array /= (ratio/MIN(ratio_array))
    signal_array_sP /= (ratio/MIN(ratio_array))
    signal_array_high_error /= (ratio/MIN(ratio_array))
    signal_array_low_error /= (ratio/MIN(ratio_array))
    spec_array /= (ratio/MIN(ratio_array))
    spec_array_high_error /= (ratio/MIN(ratio_array))
    spec_array_low_error /= (ratio/MIN(ratio_array))
  ENDIF
  ;
  ;@GJ, 2022/6/28
  ;only read the 2nd line as signal
  ;  ;repalce the raw data by smooth filtered data
  ;  IF CORRELATE(signal_array, spec_array) GT 0.8 THEN BEGIN
  ;    signal_array = spec_array
  ;    signal_array_high_error = spec_array_high_error
  ;    signal_array_low_error = spec_array_low_error
  ;  ENDIF
  

  ;GJ, 2022/5/9, automatically calculate flat ratio
  zero_H = MIN(ABS(H_array[15:N_period*0.8]), zero_ind)
  zero_ind += 15
  temp_H_array = H_array[zero_ind-10:zero_ind+10]
  ind_temp = time[0:N_ELEMENTS(temp_H_array)-1]
  linfit_result = LINFIT(ind_temp, temp_H_array)
  slew_range = 4.0 * ROUND(MAX(H_array) / ABS(linfit_result[1]))
  slew_range = slew_range - (slew_range MOD 5)
  flat_ratio_temp = 1. - slew_range*1.0/N_period
  ;@GJ, 2022/6/30, fix the flat ratio to 90% or 80% for accurate analysis
  flat_ratio = FLOOR(flat_ratio_temp*10.)/10.
  ;GJ, 2022/8/15, make sure the flat ratio is appropriate
  ;  IF flat_ratio LT 0 THEN flat_ratio = 0.9
;  IF STRPOS(filename, '20230602') NE -1 THEN flat_ratio = 0.95
  print, 'Flat portion = ', flat_ratio
  IF flat_ratio LT 0.0001 OR flat_ratio GT 1.0 THEN BEGIN
    IF OBJ_VALID(p) THEN p.close
    p = PLOT(time, H_array/MAX(H_array)*MAX(signal_array), 'g--2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
      XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='Signal and M, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='M [A.U.]')
    p.name = ' H'
    p1 = PLOT(time, Signal_array, '.r-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
    p1.name = ' Signal'
    l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
    ;p.save, filename+'01.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
    result_signal_fig = p.CopyWindow(BORDER=10, RESOLUTION=300)
    ;save the signal plot
    fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
    IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
    WRITE_PNG, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_fitplot.png', result_signal_fig
    IF OBJ_VALID(p) THEN p.close
    return, 0
  ENDIF
  ;flat_ratio=0.90;0.80;0.90

  ;GJ, 2022/5/12, shift the data
  ;calculate shift elements
  shift_ele = 0.
  ;shift the maximum to the beginning
  max_H_array = MAX(H_array, max_ind)
  shift_ele += max_ind
    ;@GJ, 2022/11/21, orthogonal
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN max_ind = 0
  temp_H_array = SHIFT(H_array, -max_ind)
  ;find the zero H index
  zero_H_array_first = MIN(ABS(temp_H_array[0:N_period/2-1]), zero_ind_first)

  ;check the slope
  IF temp_H_array[zero_ind_first+1]-temp_H_array[zero_ind_first] GT 0 THEN BEGIN
    ;@GJ, 2022/11/21, orthogonal
    IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN zero_ind_first = 0
    temp_H_array = SHIFT(temp_H_array, -zero_ind_first)
    shift_ele += zero_ind_first
  ENDIF ELSE BEGIN
    zero_ind_second = zero_ind_first + N_period/2
    ;@GJ, 2022/11/21, orthogonal
    IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN zero_ind_second = 0
    temp_H_array = SHIFT(temp_H_array, -zero_ind_second)
    shift_ele += zero_ind_second
  ENDELSE

  ;1/4 shift
  shift_1_4 = ROUND(N_period * (1. - flat_ratio)/4.);+1
    ;@GJ, 2022/11/21, orthogonal
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN shift_1_4 = 0
  temp_H_array = SHIFT(temp_H_array, shift_1_4)
  shift_ele -= shift_1_4

  ;do correction by shifting and changing the sign of H_array
  H_array = -temp_H_array
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN H_array = temp_H_array ELSE H_array = -temp_H_array
  ;@GJ, 2022/11/21, orthogonal
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN shift_ele=0
  signal_array = SHIFT(signal_array, -shift_ele)
  
  ;@GJ, 2023/4/10, adjusting the signal array
  ;assuming the max appears in the middle
  max_signal_array = MAX(signal_array, max_signal_ind)
  IF max_signal_ind LT N_period/2-1 THEN BEGIN
    signal_array = SHIFT(signal_array, N_period/2)
    signal_array_low_error = SHIFT(signal_array_low_error, N_period/2)
    signal_array_high_error = SHIFT(signal_array_high_error, N_period/2)
  ENDIF
  
  ;@GJ, 2022/12/3, do FFT and filter the 15, 30, 70 kHz errors
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
    t_ele = N_ELEMENTS(signal_array)
    X = (FINDGEN((t_ele - 1)/2) + 1)
    is_N_even = (t_ele MOD 2) EQ 0
    if (is_N_even) then $
      freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
    else $
      freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)
    signal_FFT = FFT(signal_array)
    logpower = ABS(signal_FFT)^2;ALOG10(ABS(signal_FFT)^2) ; log of Fourier power spectrum
;    iplot, freq, logpower
    ;15 kHz 6*freq
    ind = WHERE(freq EQ 6.*frequency/1000., count)
    ;force it to be 0
    IF count EQ 1 THEN signal_FFT[ind] = 0;(signal_FFT[ind-1]+signal_FFT[ind+1])/2.
    ;30 kHz
    ind = WHERE(freq EQ 12.*frequency/1000., count)
    IF count EQ 1 THEN signal_FFT[ind] = (signal_FFT[ind-1]+signal_FFT[ind+1])/2.
    ;70 kHz
    ind = WHERE(freq EQ 28.*frequency/1000., count)
    IF count EQ 1 THEN signal_FFT[ind] = (signal_FFT[ind-1]+signal_FFT[ind+1])/2.
    ;-15 kHz
    ind = WHERE(freq EQ -6.*frequency/1000., count)
    IF count EQ 1 THEN signal_FFT[ind] = 0; = (signal_FFT[ind-1]+signal_FFT[ind+1])/2.
    ;-30 kHz
    ind = WHERE(freq EQ -12.*frequency/1000., count)
    IF count EQ 1 THEN signal_FFT[ind] = (signal_FFT[ind-1]+signal_FFT[ind+1])/2.
    ;-70 kHz
    ind = WHERE(freq EQ -28.*frequency/1000., count)
    IF count EQ 1 THEN signal_FFT[ind] = (signal_FFT[ind-1]+signal_FFT[ind+1])/2.
    filtered_signal_array = FFT(signal_FFT, /inverse)
;    plot, signal_array
;    oplot, filtered_signal_array, LINESTYLE=1
    unfiltered_signal_array = signal_array
    ;update the signal
    signal_array = filtered_signal_array
  ENDIF
  signal_array_sP = SHIFT(signal_array_sP, 0, -shift_ele) ;only shift the 2nd dimension in signal_array_sP
  ;spec_array = SHIFT(spec_array, shift_ele)
  signal_array_high_error = SHIFT(signal_array_high_error, -shift_ele)
  signal_array_low_error = SHIFT(signal_array_low_error, -shift_ele)

  ;  IF OBJ_VALID(p) THEN p.close
  ;  p = PLOT(time, H_array/MAX(H_array)*MAX(signal_array_high_error), 'g--2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
  ;    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='Signal and M, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='M [A.U.]')
  ;  p.name = ' H'
  ;  p1 = PLOT(time, Signal_array, '.r-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  ;  p1.name = ' Signal'
  ;  p2 = PLOT(time, signal_array_low_error, '.b-2', LINESTYLE = 3, /OVERPLOT, /SYM_FILLED)
  ;  p3 = PLOT(time, signal_array_high_error, '.b-2', LINESTYLE = 3, /OVERPLOT, /SYM_FILLED)
  ;  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  ;  ;p.save, filename+'01.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  ;  result_signal_fig = p.CopyWindow(BORDER=10, RESOLUTION=300)
  ;;  ;@GJ, 2022/4/9, switch negative and postive parts if necessary
  ;  IF (MEAN(H_array[0:N_period/2-1]) GT MEAN(H_array[N_period/2:*])) THEN BEGIN
  ;    H_array = SHIFT(H_array, N_period/2)
  ;    signal_array = SHIFT(signal_array, N_period/2)
  ;    ;spec_array = SHIFT(spec_array, N_period/2)
  ;    signal_array_high_error = SHIFT(signal_array_high_error, N_period/2)
  ;    signal_array_low_error = SHIFT(signal_array_low_error, N_period/2)
  ;  ENDIF

  ;find the flat portion index
  half_rise_time = ROUND(N_period * (1. - flat_ratio)/2.)
  start_ind = half_rise_time;*2.
  print, 'start index: ', start_ind
  signal_array_slew_n = signal_array[0:start_ind-1]
  signal_array_slew_p = signal_array[N_period/2:start_ind+N_period/2-1]
  IF start_ind GT (N_period/2-1-N_period/50) THEN BEGIN
    ;save the signal plot
    fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
    IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
    WRITE_PNG, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_fitplot.png', result_signal_fig
    IF OBJ_VALID(p) THEN p.close
    RETURN, 0
  ENDIF
  signal_array_n = signal_array[start_ind:N_period/2-1-N_period/50]
  signal_array_p = signal_array[start_ind+N_period/2:N_period-1-N_period/50]
  ;  AUC = ABS(TOTAL(signal_array_p)) * delta_t
  AUC_slew_p = ABS(TOTAL(signal_array_slew_p)) * delta_t
  AUC_slew_n = ABS(TOTAL(signal_array_slew_n)) * delta_t
  AUC_hold_p = ABS(TOTAL(signal_array_p)) * delta_t
  AUC_hold_n = ABS(TOTAL(signal_array_n)) * delta_t
  AUC_total_p = AUC_slew_p + AUC_hold_p
  AUC_total_n = AUC_slew_n + AUC_hold_n
  AUC_hold = (AUC_hold_p+AUC_hold_n)/2.
  AUC_slew = (AUC_slew_p+AUC_slew_n)/2.
  AUC_total = (AUC_total_p+AUC_total_n)/2.
  AUC_slew_ratio = AUC_slew/AUC_total*100.
  PRINT, 'AUC slew_p: ', AUC_slew_p, ' [signal*us]'
  PRINT, 'AUC hold_p: ', AUC_hold_p, ' [signal*us]'
  PRINT, 'AUC slew_n: ', AUC_slew_n, ' [signal*us]'
  PRINT, 'AUC hold_n: ', AUC_hold_n, ' [signal*us]'

  ;find the peak H
  H_peak_n = MEAN(H_array[start_ind:N_period/2-1])
  H_peak_p = MEAN(H_array[start_ind+N_period/2:N_period-1])

  ;  window, 0
  ;  ; Set up variables for the plot. Normally, these values would be
  ;  ; passed into the program as positional and keyword parameters.
  ;  xtitle = 'time [us]'
  ;  ytitle = 'Signal [A.U.]'
  ;  title = 'H and Signal'
  ;  position = [0.125, 0.125, 0.9, 0.925]
  ;
  ;  ; Set up a "window" for the plot. The PostScript output will have
  ;  ; the same aspect ratio as the graphics window on the display.
  ;  cgDisplay, 600, 500, Title='H and Signal'
  ;  ; Draw the line plot with no data
  ;  cgPlot, time, signal_array, Title=title, XTitle=xtitle, YTitle=ytitle, $
  ;    Position=position, /NoData, YRange=[MIN(signal_array_low_error), MAX(signal_array_high_error)], YStyle=1;, $
  ; ;   OUTPUT='JPEG', OUTFILENAME=filename
  ;  ; Fill in the error estimates.
  ;  cgColorFill, [time, Reverse(time), time[0]], $
  ;    [signal_array_high_error, Reverse(signal_array_low_error), signal_array_high_error[0]], $
  ;    Color='sky blue'
  ;  ; Draw the line plot with no data
  ;  cgPlotS, time, H_array/MAX(H_array)*MAX(signal_array), Color='blue', $
  ;    Thick=2, LINESTYLE = 2
  ;  cgPlotS, time, signal_array, Color='red', PSym=-16,  SymSize=0.5, SymColor='red', $
  ;    Thick=2


  ;  ;GJ, 2022/6/18, save the figure
  ;  fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
  ;  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  ;  WRITE_JPEG, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_curve.jpg', TVRD(), QUALITY=25

  ;calculate the M by intergrating signal by time
  ;temp_M_array = TOTAL(signal_array-MEAN(signal_array), /CUMULATIVE) * delta_t
  ;shift the M_array to middle
  ;temp_M_array = temp_M_array; - MEAN(temp_M_array)

  percentile_90_n = cgPercentiles(signal_array[0:N_period/2-1], Percentiles=0.90)
  percentile_10_p = cgPercentiles(signal_array[N_period/2:*], Percentiles=0.10)
  tail_noise_n = STDDEV(signal_array[N_period/2-50:N_period/2-1])
  tail_noise_p = STDDEV(signal_array[N_period-50:N_period-1])
  tail_mean_n = MEAN(signal_array[N_period/2-50:N_period/2-1])
  tail_mean_p = MEAN(signal_array[N_period-50:N_period-1])
  ;  M_array_n = TOTAL(signal_array[0:N_period/2-1]-MIN([percentile_90_n, tail_mean_n]), /CUMULATIVE) * delta_t
  ;  M_array_p = TOTAL(signal_array[N_period/2:*]-MAX([percentile_10_p, tail_mean_p]), /CUMULATIVE) * delta_t
  M_array_n = TOTAL(signal_array[0:N_period/2-1]-percentile_90_n[0], /CUMULATIVE) * delta_t
  M_array_p = TOTAL(signal_array[N_period/2:*]-percentile_10_p[0], /CUMULATIVE) * delta_t
  M_array = [M_array_n, M_array_n[N_period/2-1]+M_array_p]
  M_array = M_array - MEAN(M_array)
  ;@GJ, 2022/11/21, add the orthogonal 
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
    M_array = TOTAL(signal_array, /CUMULATIVE) * delta_t
    M_array = M_array - MIN(M_array)
  ENDIF
  IF STRPOS(filename, 'WangQianData') NE -1 THEN M_array = M_array - MEAN(M_array)
  
  ;M_array = temp_M_array
  ;M_array[start_ind:N_period/2-1] = M_array_n
  ;M_array[start_ind+N_period/2:N_period-1] = M_array_p

  M_array_n =M_array[start_ind:N_period/2-1-N_period/50]
  M_array_p = M_array[start_ind+N_period/2:N_period-1-N_period/50]

  M_array_n_diff = MAX(M_array_n) - MIN(M_array_n)
  M_array_p_diff = MAX(M_array_p) - MIN(M_array_p)
  PRINT, 'M_diff hold: ', M_array_p_diff, ' [signal*us]'

  ;GJ, 2022/6/16, plot the original signal and M
  ;
  ;if p is available
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
    max_H = MAX(H_array[0:N_period/2-1], start_ind)
    H_array = SHIFT(H_array, -start_ind)
    Signal_array = SHIFT(Signal_array, -start_ind)
    ;@GJ, 2022/12/3, shift the unfiltered signal
    unfiltered_signal_array = SHIFT(unfiltered_signal_array, -start_ind)
    signal_array_low_error = SHIFT(signal_array_low_error, -start_ind)
    signal_array_high_error = SHIFT(signal_array_high_error, -start_ind)
    M_array = SHIFT(M_array, -start_ind)
    ;start_ind = N_period/4.
  ENDIF
  
  
  IF OBJ_VALID(p) THEN p.close
  ori_title = 'Signal and M, '+STRTRIM(STRING(frequency),1)+' Hz'
  
  IF STRPOS(filename, 'WangQianData') NE -1  OR STRPOS(filename, '20230602') NE -1 THEN BEGIN
    p = PLOT(time, H_array, 'g-.2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title=ori_title, xtitle='time [us]', ytitle='Signal, H and M [A.U.]')
  ENDIF ELSE BEGIN
    p = PLOT(time, -H_array, 'g-.2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
      XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title=ori_title, xtitle='time [us]', ytitle='Signal, H and M [A.U.]')
  ENDELSE
  
    p.name = ' H'
  IF N_ELEMENTS(unfiltered_signal_array) EQ 0 THEN unfiltered_signal_array = Signal_array - MEAN(Signal_array)
  p1 = PLOT(time, unfiltered_signal_array/MAX(ABS(unfiltered_signal_array))*MAX(H_array), 'c.-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p1.name = ' Signal'
  IF ABS(M_array[start_ind]) LT 0.2 THEN BEGIN
    ;max_signal_plot = MAX(Signal_array/ABS(Signal_array[start_ind])*ABS(H_array[start_ind]))
    p4 = PLOT(time, M_array/MAX(ABS(M_array))*MAX(H_array), '.r-2', LINESTYLE = 3, /OVERPLOT, /SYM_FILLED)
  ENDIF ELSE BEGIN
    p4 = PLOT(time, M_array/MAX(ABS(M_array))*MAX(H_array), '.r-2', LINESTYLE = 3, /OVERPLOT, /SYM_FILLED)
  ENDELSE
  p4.name = ' M'
  p5 = PLOT(time, M_array*0., '.-1', LINESTYLE = 0, /OVERPLOT)
  IF STRPOS(filename, '20221119') NE -1 THEN BEGIN
    p5.name = ' zero'
    l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  ENDIF
  ;add the error range
  p2 = PLOT(time, (signal_array_low_error-MEAN(Signal_array_low_error))/MAX(ABS(unfiltered_signal_array))*MAX(H_array), '.b-2', LINESTYLE = 3, /OVERPLOT, /SYM_FILLED)
  p3 = PLOT(time, (signal_array_high_error-MEAN(Signal_array_high_error))/MAX(ABS(unfiltered_signal_array))*MAX(H_array), '.b-2', LINESTYLE = 3, /OVERPLOT, /SYM_FILLED)
  fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  p.save, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_signalM.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  result00 = p.CopyWindow(BORDER=10, RESOLUTION=300)
  size_single = SIZE(result00, /DIMENSIONS)
  ;put these figures together as one png file
  size_single = SIZE(result00, /DIMENSIONS)
  result_fig = BYTARR(3, 3*size_single[1], 2*size_single[2])
  result_fig[*, 0:size_single[1]-1, size_single[2]:2*size_single[2]-1] = CONGRID(result00, 3, size_single[1], size_single[2])
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF
  
  ;@GJ, 2022/11/21, do curve fitting for relaxation time estimation
  ;@GJ, 2023/5/8, add the process of sine excitation
  IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
    print, 'this is an orthogonal or sine excitation!'
    ;@GJ, 2022/11/28, calculate the phase delay angle
    zeroHac = MIN(ABS(H_array[N_period/4.:N_period*3./4.]), zeroHac_ind)
    maxM = MAX(M_array[N_period/4.:N_period*3./4.], maxM_ind)
    delay_phase_angle = (maxM_ind - zeroHac_ind) * 360. / N_period; in degree
    print, 'delay phase angle: ', delay_phase_angle, ' degree'
    
    ;@GJ, 2022/11/28, calculate the signal with relaxation time
    Hx0 = MAX(H_array)
    ;based on Yanjun's description
    Hz0 = FIX(STRMID(filename, 5, 2, /REVERSE_OFFSET)) mod 12
    IF Hz0 EQ 0 THEN Hz0 = 12
    ;calculate the signal adiabatic
    particle_size_array = FINDGEN(161)*0.5+5.
    relaxation_time_array = FINDGEN(501)/500. * 50.
    ;print, 'relaxation time array: ', relaxation_time_array
    dist_array = DBLARR(N_ELEMENTS(relaxation_time_array), N_ELEMENTS(particle_size_array))
    corr_array = dist_array*0.
    maxM_nat_array = dist_array*0.
    M_ratio_array = dist_array*0.
    temp_signal_nat_array = DBLARR(N_period-1, N_ELEMENTS(relaxation_time_array), N_ELEMENTS(particle_size_array))
    temp_M_nat_array = DBLARR(N_period-1, N_ELEMENTS(relaxation_time_array), N_ELEMENTS(particle_size_array))
    FOR k = 0, N_ELEMENTS(particle_size_array)-1 DO BEGIN
      particle_size = particle_size_array[k]
      orthogonal_excitation_diff_relax, frequency/1000., particle_size, signal_z_array, M_z_array, relaxation_time_array, Hx0, DOUBLE(Hz0), signal_x_array, M_x_array
      IF N_ELEMENTS(relaxation_time_array) EQ 1 THEN BEGIN
        i=0
        maxM_z = MAX(REFORM(M_z_array[N_period/4.:N_period*3./4., i]), maxM_z_ind)
        M_ratio = 1./ maxM_z * maxM
        temp_signal = REFORM(signal_z_array[*,i]) * M_ratio
      ENDIF

      IF STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
        Hz0=0.
        orthogonal_excitation_diff_relax, frequency/1000., particle_size, signal_z_array, M_z_array, relaxation_time_array, Hx0, DOUBLE(Hz0), signal_x_array, M_x_array
        temp_signal = signal_x_array
      ENDIF
      
      FOR j=0, N_ELEMENTS(relaxation_time_array)-1 DO BEGIN
        temp_signal_nat_array[*,j, k] = non_adiabatic(time[0:N_period-2], temp_signal[0:N_period-2], relaxation_time_array[j])
        maxM = MAX(M_array[N_period/4.:N_period*3./4.], maxM_ind)
        temp_M_nat = TOTAL(REFORM(temp_signal_nat_array[*,j, k]), /CUMULATIVE)*delta_t
        maxM_nat_array[j, k] = MAX(temp_M_nat[N_period/4.:N_period*3./4.], maxM_nat_ind)
        dist_array[j, k] = maxM_ind-maxM_nat_ind

        ;find the negative indices of negative H
        neg_H_ind = WHERE(H_array LE 0, count_H)
        ;do the corrections
        IF count_H GT 0 THEN maxM_nat_array[j, k] -= MIN(temp_M_nat[neg_H_ind])
        M_ratio_array[j, k] = 1./ maxM_nat_array[j, k] * maxM
        temp_signal_nat_array[*,j, k] *= M_ratio_array[j, k]
        fitted_M = TOTAL(REFORM(temp_signal_nat_array[*,j, k]), /CUMULATIVE) * delta_t
        IF count_H GT 0 THEN BEGIN
          temp_M_nat_array[*,j, k] = (fitted_M - MIN(fitted_M[neg_H_ind]))
        ENDIF ELSE BEGIN
          temp_M_nat_array[*,j, k] = (fitted_M - MIN(fitted_M))
        ENDELSE
        corr_array[j, k] = CORRELATE(temp_signal, REFORM(temp_signal_nat_array[*,j, k]))
        IF STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
          corr_array[j, k] = CORRELATE(signal_array[0:N_period-2], REFORM(temp_signal_nat_array[*,j, k]))
        ENDIF
      ENDFOR
    ENDFOR
        
;    iplot, relaxation_time_array, corr_array, title='correlation'
    ;iplot, ABS(dist_array), title='zero point dist'
    min_dist = MIN(ABS(dist_array), min_dist_ind)
    IF (min_dist_ind GE 5) AND (min_dist_ind LE N_ELEMENTS(relaxation_time_array)-6) THEN BEGIN
      max_cor = MAX(corr_array[min_dist_ind-5:min_dist_ind+5], max_cor_ind)
      IF max_cor_ind NE 0 AND max_cor_ind NE 10 THEN BEGIN
        max_cor_ind += min_dist_ind-5
      ENDIF ELSE BEGIN
        max_cor_ind = min_dist_ind
      ENDELSE
    ENDIF ELSE BEGIN
      max_cor_ind = min_dist_ind
    ENDELSE
    
    IF STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
      max_cor = MAX(corr_array, max_cor_ind)
      max_cor_ind_2d = ARRAY_INDICES(corr_array, max_cor_ind)
      IF max_cor_ind_2d[0] EQ 0 AND max_cor_ind_2d[1] EQ N_ELEMENTS(particle_size_array)-1 THEN BEGIN
        max_cor_ind_2d = [0,0]
      ENDIF
    ENDIF
    
    print, 'relaxation time: ', relaxation_time_array[max_cor_ind_2d[0]]
    print, 'particle size: ', particle_size_array[max_cor_ind_2d[1]]
    fitted_signal = REFORM(temp_signal_nat_array[*, max_cor_ind_2d[0], max_cor_ind_2d[1]])
    fitted_signal *= (STDDEV(Signal_array)/STDDEV(fitted_signal))
    fitted_M = REFORM(temp_M_nat_array[*,max_cor_ind_2d[0], max_cor_ind_2d[1]])
    fitted_M -= MEAN(fitted_M)
    fitted_M *= (STDDEV(M_array)/STDDEV(fitted_M))
    
    ;@GJ, 2022/11/27, save the fitted signal
    IF OBJ_VALID(p) THEN p.close
    fit_title = 'Fitted Signal(tao:'+STRTRIM(STRING(relaxation_time_array[max_cor_ind_2d[0]]),1)+'us,Size:'+STRTRIM(STRING(particle_size_array[max_cor_ind_2d[1]]),1)+'nm';Hac:'+STRTRIM(STRING(Hx0),1)+'mT,Hdc:'+STRTRIM(STRING(Hz0),1)+'mT)'
    p = PLOT(time[0:N_period-2], Signal_array, 'r-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
      XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title=fit_title, xtitle='time [us]', ytitle='Signal')
    p.name = ' Signal'
    p1 = PLOT(time[0:N_period-2], fitted_signal, 'b.-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
    p1.name = ' Fitted'
    IF STRPOS(filename, '20221119') NE -1 OR STRPOS(filename, 'WangQianData') NE -1 THEN BEGIN
      ;don't plot p2
    ENDIF ELSE BEGIN
      p2 = PLOT(time[0:N_period-2], unfiltered_signal_array, 'g.-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
      p2.name = ' Original
    ENDELSE
    p.save, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_FitSignal.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
    result00 = p.CopyWindow(BORDER=10, RESOLUTION=300)
    ;put these figures together as one png file
    result_fig[*, size_single[1]:2*size_single[1]-1, size_single[2]:2*size_single[2]-1] = CONGRID(result00, 3, size_single[1], size_single[2])
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
      WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
      WSET, drawpicid
      TVSCL, new_pic, TRUE=1
    ENDIF
    
    ;@GJ, 2022/11/27, save the fitted M-H
    IF OBJ_VALID(p) THEN p.close
    p = PLOT(H_array, M_array, 'r-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
      XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='M-H (tao:'+STRTRIM(STRING(relaxation_time_array[max_cor_ind_2d[0]]),1)+'us,Size:'+STRTRIM(STRING(particle_size_array[max_cor_ind_2d[1]]),1)+'nm', xtitle='Hac [mT]', ytitle='M [A.U.]')
    p.name = ' M'
    ;@GJ, 2022/11/25, minus the zero part
    p1 = PLOT(H_array, fitted_M, 'b.-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
    p1.name = ' Fitted'
    p.save, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_MH.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
    result00 = p.CopyWindow(BORDER=10, RESOLUTION=300)
    ;put these figures together as one png file
    result_fig[*, 2*size_single[1]:3*size_single[1]-1, size_single[2]:2*size_single[2]-1] = CONGRID(result00, 3, size_single[1], size_single[2])
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
      WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
      WSET, drawpicid
      TVSCL, new_pic, TRUE=1
    ENDIF
    
    ;@GJ, 2022/11/29, plot the error signals
    IF OBJ_VALID(p) THEN p.close
    error_title = 'Fit Error';(tao:'+STRTRIM(STRING(relaxation_time_array[max_cor_ind]),1)+'us';,Freq:'+STRTRIM(STRING(frequency),1)+'Hz,Hac:'+STRTRIM(STRING(Hx0),1)+'mT,Hdc:'+STRTRIM(STRING(Hz0),1)+'mT)'
    error_signal = Signal_array-fitted_signal
    p = PLOT(time[0:N_period-2], error_signal, 'r-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
      XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title=error_title, xtitle='time [us]', ytitle='Signal')
    p.name = ' Error'
    p.save, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_FitError.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
    result00 = p.CopyWindow(BORDER=10, RESOLUTION=300)
    ;put these figures together as one png file
    result_fig[*, 0:size_single[1]-1, 0:size_single[2]-1] = CONGRID(result00, 3, size_single[1], size_single[2])
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
      WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
      WSET, drawpicid
      TVSCL, new_pic, TRUE=1
    ENDIF
    
    ; Frequency for FFT in kHz, unit
    ; N is an integer giving the number of elements in a particular dimension
    ; T is a floating-point number giving the sampling interval
    t_ele = N_ELEMENTS(error_signal)
    X = (FINDGEN((t_ele - 1)/2) + 1)
    is_N_even = (t_ele MOD 2) EQ 0
    if (is_N_even) then $
      freq = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t*0.001) $
    else $
      freq = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t*0.001)
    error_signal_FFT = FFT(error_signal)
    error_powerSpectrum = ABS(error_signal_FFT)^2
    ;@GJ, 2022/11/29, plot the error signals
    IF OBJ_VALID(p) THEN p.close
    error_title = 'Fit Error FFT'
    p = PLOT(freq[0:29], error_powerSpectrum[0:29], 'r-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
      XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title=error_title, xtitle='Freq [kHz]', ytitle='Power')
    p.name = ' Error FFT'
    p.save, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_FitErrorFFT.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
    result00 = p.CopyWindow(BORDER=10, RESOLUTION=300)
    ;put these figures together as one png file
    result_fig[*, size_single[1]:2*size_single[1]-1, 0:size_single[2]-1] = CONGRID(result00, 3, size_single[1], size_single[2])
    IF N_ELEMENTS(sState) GT 0 THEN BEGIN
      new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
      WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
      WSET, drawpicid
      TVSCL, new_pic, TRUE=1
    ENDIF
    
    fit_result = {$
      frequency:frequency, $
      Hx0:Hx0, $
      Hz0:Hz0, $
      time:time[0:N_period-2], $
      H:H_array, $
      M:M_array, $
      H_peak_n:H_peak_n, $
      H_peak_p:H_peak_p, $
      delay_phase_angle:delay_phase_angle, $
      fitted_M:fitted_M, $
      relaxation_orth:relaxation_time_array[max_cor_ind_2d[0]], $
      particle_size:particle_size_array[max_cor_ind_2d[1]],$
      particle_size_mono:particle_size_array[max_cor_ind_2d[1]],$
      tao_debye_mono:relaxation_time_array[max_cor_ind_2d[0]], $
      corre_biex:0, $
      signal_ratio:M_ratio_array[max_cor_ind], $
      correlation_coeff:corr_array[max_cor_ind], $
      signal:signal_array, $
      AUC_total_p: AUC_total_p,$
      AUC_total_n: AUC_total_n,$
      AUC_slew_p: AUC_slew_p, $
      AUC_slew_n: AUC_slew_n, $
      AUC_hold_p: AUC_hold_p,$
      AUC_hold_n: AUC_hold_n,$
      AUC_slew: AUC_slew, $
      AUC_hold: AUC_hold, $
      AUC_total: AUC_total, $
      AUC_slew_ratio: AUC_slew_ratio, $
      unfiltered_signal: unfiltered_signal_array, $
      filtered_signal: signal_array, $
      fitted_signal:fitted_signal}

    ;close the window
    IF OBJ_VALID(p) THEN p.close
    
    ;@GJ, 2022/11/27, return the fitted result
    RETURN, fit_result
  ENDIF
  
  ;define X and Y for curve fitting
  IF KEYWORD_SET(signalOnly) THEN BEGIN
    Y = ABS(signal_array_n)-MEAN(signal_array)
  ENDIF ELSE BEGIN
    IF MIN(M_array_p) LT 0 THEN Y = M_array_p-2.*MIN(M_array_p) ELSE Y = M_array_p;[start_ind:N_period/2-1]
  ENDELSE
  X = time[0:N_ELEMENTS(Y)-1]

  ;GJ, 2022/2/9, set for the curve fit
  relax_H = -H_array[start_ind:N_period/2-1]
  relax_M = Y

  ;Define a vector of weights.
  IF KEYWORD_SET(signalOnly) THEN weights = Y ELSE weights = 1.0/Y;

  ;Provide an initial guess of the function’s parameters.
  IF KEYWORD_SET(signalOnly) THEN A = [MAX(Y),-0.1,MIN(Y)] ELSE A = [-10.0,-0.1,MAX(Y)]

  ;Compute the parameters.
  yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct')

  relax_mono_Mfit = yfit
  ;Print the parameters returned in A.
  PRINT, 'Function parameters: ', A

  tao = 1./ABS(A[1])
  print, 'Relaxation time = ', tao, ' [us]'

  ;do the biexponential curve fit
  ;Provide an initial guess of the function’s parameters.
  IF KEYWORD_SET(signalOnly) THEN B = [MAX(Y)/2.0, -0.01, MAX(Y)/2.0, -0.1, MIN(Y), 0.] ELSE B = [-MAX(Y)/2.0, -0.01, -MAX(Y)/2.0, -0.1, MIN(Y), 0.]; ;@GJ, 2022/7/10, B = [-5.0, -0.01, -5.0, -0.1, MAX(Y), 0.]
  FITA = [1,1,1,1,1,0]
  ;Compute the parameters.
  yfit_biex = CURVEFIT(X, Y, weights, B, SIGMA, FITA=FITA, FUNCTION_NAME='gfunct_biex')
  ;calcculate the dy/dt for signal
  bx_1 = EXP(B[1] * (X-B[5]))
  bx_2 = EXP(B[3] * (X-B[5]))
  ;  yfit_biex_dev = B[0] * B[1] * bx_1 + B[2] * B[3] * bx_2
  yfit_biex_dev = yfit_biex[1:*]-yfit_biex[0:N_ELEMENTS(X)-2]

  x_ext =[-REVERSE(X),X[1:*]]
  bx_1 = EXP(B[1] * (X_ext-B[5]))
  bx_2 = EXP(B[3] * (X_ext-B[5]))
  y_ext = B[0] * bx_1 + B[2] * bx_2 + B[4]
  y_ext_min = MIN(ABS(y_ext+B[4]), y_ext_min_ind)
  t0 = x_ext[y_ext_min_ind]

  ;change Y back by adding up
  IF MIN(M_array_p) LT 0 THEN BEGIN
    Y = Y + 2.*MIN(M_array_p)
    yfit = yfit + 2.*MIN(M_array_p)
    yfit_biex = yfit_biex + 2.*MIN(M_array_p)
  ENDIF

  relax_biex_Mfit = yfit_biex
  ;Print the parameters returned in A.
  PRINT, 'Function parameters (Bi-ex): ', B

  corre_biex = CORRELATE(y, yfit_biex)
  PRINT, 'correlation coeff (Bi-ex): ', corre_biex

  ;GJ, 2022/2/11, plot the fitting
  ;if p is available
  IF OBJ_VALID(p) THEN p.close
  p = PLOT(X, Y, 'o-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='Exponential Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='M [A.U.]')
  p.name = ' Original data'
  p1 = PLOT(X, yfit, 'c-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p1.name = ' Mono-ex fitting'
  p2 = PLOT(X, yfit_biex, 'r-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p2.name = ' Bi-ex fitting'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.5])
  ;p.save, filename+'01.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  result01 = p.CopyWindow(BORDER=10, RESOLUTION=300)
  result_fig[*, 2*size_single[1]:3*size_single[1]-1, 0:size_single[2]-1] = CONGRID(result01, 3, size_single[1], size_single[2])
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF

  ;GJ, 2022/2/10, do correctin to the coefficients for correction Neel composition
  B[0] = B[0]*exp(ABS(t0)*B[1])/2.
  B[2] = B[2]*exp(ABS(t0)*B[3])/2.

  bi_tao_1 = 1./ABS(B[1])
  bi_tao_1_comp = ABS(B[0]/(B[0]+B[2]))
  bi_tao_2 = 1./ABS(B[3])
  bi_tao_2_comp = ABS(B[2]/(B[0]+B[2]))
  print, 'Relaxation time_1 = ', bi_tao_1, ' [us] with',  bi_tao_1_comp*100., '%'
  print, 'Relaxation time_2 = ', bi_tao_2, ' [us] with',  bi_tao_2_comp*100., '%'

  ;GJ, 2022/2/10, Calculate the combined tao
  bi_tao_combine = 1./(1./bi_tao_1*bi_tao_1_comp + 1./bi_tao_2*bi_tao_2_comp)

  ;GJ, 2022/2/11, define the particle size
  ;N_particle_size = 151
  ;particle_size_array = FINDGEN(N_particle_size)/(N_particle_size-1)*15. + 20.

  ;find the particle size as 24 nm for Synomag
  N_particle_size = 1
  ;Synomag-D particle size as measured
  ;@GJ, 2022/8/15, using Synamag-D in PBS as the basic particle size
  particle_size_array = [24.4924];PBS, [25.2579], 1% Glycerol
  IF STRPOS(filename_base, '20220616') NE -1 THEN particle_size_array = [24.4924];PBS;[25.2579], 1% Glycerol;[23.9678];[24.0]; synomag-D
  IF STRPOS(filename_base, '20220523') NE -1 THEN particle_size_array = [9.3797]; Gelatin at 2022/5/23
;  IF STRPOS(filename_base, '20220815') NE -1 THEN particle_size_array = [7.8833]; for Gelatin at 2022/8/15
  
  ;GJ, 2022/2/11, define debye tao
  N_debye_tao = 600; 200
  ;GJ, 2022/1/30, calculate Debye model
  relaxation_time = FINDGEN(N_debye_tao)/10.+1.

  MSE_array = DBLARR(N_particle_size, N_debye_tao)*0.+10000.
  correlation_array = DBLARR(N_particle_size, N_debye_tao)*0.

  FOR j=0, N_particle_size-1 DO BEGIN
    ;GJ, 2022/2/1, Chinese New Year, Debye relaxation model
    Nelem = N_ELEMENTS(H_array)
    IF ARG_PRESENT(particle_size) EQ 0 THEN BEGIN
      ;      particle_size = 27.;19.;40.;50.;19.;40.;50;70.;19.;30.
      particle_size = particle_size_array[j]
    ENDIF ELSE BEGIN
      ;stop the circulation
      j= N_particle_size-1
    ENDELSE
    Msat_T = 0.551 ; T/u0
    Msat_kAm = 0.551/4./!PI*10000.; kA/m
    T_p = 20.;36.5; degree
    T_p_kelvin = 273.15 + T_p ; in Kelvin
    beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
;    beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
    FWHM_particle = 4.16/beta_particle/3.5; mT/mm
    M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))

    start_s = N_period/2. - half_rise_time*2.
    IF start_s LT 0 THEN BEGIN
      ;save the signal plot
      fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
      IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
      IF N_ELEMENTS(result_signal_fig) GT 10 THEN WRITE_PNG, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_fitplot.png', result_signal_fig
      IF OBJ_VALID(p) THEN p.close
      RETURN, 0
    ENDIF
    ;GJ, 2022/2/9, add min_ind to end_s
    end_s = N_period/2. + half_rise_time
    ;define X and Y for curve fitting
    IF KEYWORD_SET(signalOnly) THEN BEGIN
      signal_langevin_temp = (M[1:*]-M[0:Nelem-1])/delta_t; kA/m/us
      signal_meas = signal_array[start_s:end_s];[0:Nelem-2]
    ENDIF ELSE BEGIN
      signal_langevin_temp = M
      signal_meas =  M_array[start_s:end_s]
    ENDELSE
    signal_langevin = signal_langevin_temp[start_s:end_s]/MEAN(ABS(signal_langevin_temp))*MEAN(ABS(signal_meas))
    time_half = time[0:end_s-start_s];[0:Nelem-2]

    ;    window, 4
    ;    plot, time_half, signal_langevin, YRANGE=[MIN([signal_langevin, signal_meas])-0.2, MAX([signal_langevin, signal_meas])+0.2], title='Langevin & Signal, '+STRTRIM(STRING(frequency),1)+' Hz', XTITLE='time [us]', YTITLE='signal [a.u.]'
    ;    oplot, time_half, signal_meas, line=2
    ;    oplot, time_half, ABS(H_array[start_s:end_s])

    FOR i = 0, N_debye_tao-1 DO BEGIN
      signal_nat = non_adiabatic(time_half, signal_langevin, relaxation_time[i])
      signal_nat = signal_nat/MEAN(signal_nat)*MEAN(signal_meas)
      correlation_array[j, i] = CORRELATE(signal_nat, signal_meas)
      MSE_array[j, i] = calMSE(signal_nat, signal_meas)
    ENDFOR

  ENDFOR

  ;max_corr = MAX(correlation_array, max_ind)
  min_MSE = MIN(MSE_array, min_ind)
;  iplot, REFORM(MSE_array[0,*]), title='MSE array debye'
  min_ind2d = ARRAY_INDICES(MSE_array, min_ind)

  particle_size_mono = particle_size_array[min_ind2d[0]]
  print, 'Particle size = ', particle_size_mono, ' nm'
  tao_debye_mono = relaxation_time[min_ind2d[1]]
  print, '(Debye) relaxation time = ', tao_debye_mono
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size_mono^3) /24./1.380469/T_p_kelvin ; in 1/mT
;  beta_particle = Msat_T*(particle_size_mono^3) /24./1.380469/309.65 ; in 1/mT
  FWHM_particle = 4.16/beta_particle/3.5; mT/mm
  M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
  IF KEYWORD_SET(signalOnly) THEN BEGIN
    signal_langevin_temp = (M[1:*]-M[0:Nelem-1])/delta_t; kA/m/us
    signal_meas = signal_array[start_s:end_s];[0:Nelem-2]
  ENDIF ELSE BEGIN
    signal_langevin_temp = M
    signal_meas =  M_array[start_s:end_s]
  ENDELSE
  signal_langevin = signal_langevin_temp[start_s:end_s]/MEAN(ABS(signal_langevin_temp))*MEAN(ABS(signal_meas))
  signal_nat = non_adiabatic(time_half, signal_langevin, tao_debye_mono)
  signal_nat_mono = signal_nat/MEAN(signal_nat)*MEAN(signal_meas)
  ;  oplot, time, signal_nat, line=3
  corre_Debye_Mono = CORRELATE(signal_meas, signal_nat_mono)
  PRINT, 'Debye correlation coeff (Mono): ', corre_Debye_Mono

  ;  window, 4
  ;  cgplot, time_half, signal_meas, LINESTYLE = 0, thick = 2, YRANGE=[MIN([signal_langevin, signal_meas, signal_nat])-0.2, MAX([signal_langevin, signal_meas, signal_nat])+0.2], title='Debye Mono-ex Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='M [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 'blue';
  ;  cgplot, time_half, signal_nat_mono, thick = 2, LINESTYLE = 1, color='red',/overplot
  ;  cgplot, time_half, signal_langevin, thick = 2, color=0 ,/overplot
  IF OBJ_VALID(p) THEN p.close
  p = PLOT(time_half, signal_meas, 'o-2', AXIS_STYLE=1, FONT_SIZE=9, $;LAYOUT=[2,2,3], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='Debye Mono-ex Relaxation, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='M [A.U.]')
  p.name = ' Original data'
  p1 = PLOT(time_half, signal_nat_mono, '+c-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p1.name = ' Debye model fitting'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.7, 0.8])
  ;p.save, filename+'02.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  result02 = p.CopyWindow(BORDER=10, RESOLUTION=300)
  result_fig[*, size_single[1]:2*size_single[1]-1, 0:size_single[2]-1] = CONGRID(result02, 3, size_single[1], size_single[2])
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF

  ;GJ, 2022/2/5, define M-H curve
  ;  iplot, M_H_array[0,*], M_H_array[1,*]
  IF OBJ_VALID(p) THEN p.close
  p = PLOT(H_array, M_array, 'o-2', AXIS_STYLE=1, FONT_SIZE=9, $;LAYOUT=[2,2,1], $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], XRANGE=[MIN(H_array)*1.1, MAX(H_array)*1.1], YRANGE=[MIN(M_array)*1.1, MAX(M_array)*1.1], title='M-H Curve', xtitle='H [mT/u0]', ytitle='M [A.U.]')
  ;GJ, 2022/2/11, plot the fitting
  p.name = ' Original'
  p1 = PLOT(relax_H, relax_biex_Mfit, '*r-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p1.name = ' Bi-ex'
  debye_H = H_array[start_s:end_s]
  debye_M = M_array[start_s:end_s]
  debye_mono_Mfit = signal_nat_mono
  p2 = PLOT(debye_H, debye_mono_Mfit, '+c-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p2.name = ' Debye'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.4, 0.9])
  ;p.save, filename+'03.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  result03 = p.CopyWindow(BORDER=10, RESOLUTION=300)
  result_fig[*, 0:size_single[1]-1, 0:size_single[2]-1] = CONGRID(result03, 3, size_single[1], size_single[2])
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF

  ;plot the H and signal
  IF OBJ_VALID(p) THEN p.close
  p = PLOT(time, H_array/MAX(H_array)*MAX(signal_array), 'hg-2', AXIS_STYLE=1, FONT_SIZE=9, YRANGE=[-ABS(MAX(signal_array))*1.2,ABS(MAX(signal_array))*1.2], $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='H and Signal, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='Signal & H [A.U.]')
  p.name = ' H'
  p1 = PLOT(time, signal_array-MEAN(signal_array), 'o-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p1.name = ' Signal'
  debye_mono_Mfit_dev = (signal_nat_mono[1:end_s-start_s-1]-signal_nat_mono[0:end_s-start_s-2])/delta_t
  p2 = PLOT(time[start_s:end_s-2], debye_mono_Mfit_dev, '+c-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p2.name = ' Debye'
  p3 = PLOT(time[N_period/2+start_ind:*], yfit_biex_dev, '*r-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p3.name = ' Bi-ex'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.46])
  ;p.save, filename+'04.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  result04 = p.CopyWindow(BORDER=10, RESOLUTION=300)


  ;GJ, 2022/6/16, plot the signal during the inversion recovery process
  IF OBJ_VALID(p) THEN p.close
  p = PLOT(time, H_array/ABS(H_array[start_ind])*ABS(M_array[start_ind]), 'hg-2', AXIS_STYLE=1, FONT_SIZE=9, YRANGE=[-ABS(MAX(M_array))*1.2,ABS(MAX(M_array))*1.2], $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='H and M, '+STRTRIM(STRING(frequency),1)+' Hz', xtitle='time [us]', ytitle='M & H [A.U.]')
  p.name = ' H'
  p1 = PLOT(time, M_array, 'o-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p1.name = ' M'
  debye_H = H_array[start_s:end_s]
  debye_M = M_array[start_s:end_s]
  debye_mono_Mfit = signal_nat_mono
  p2 = PLOT(time[start_s:end_s-2], debye_mono_Mfit, '+c-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p2.name = ' Debye'
  p3 = PLOT(time[N_period/2+start_ind:*], relax_biex_Mfit, '*r-2', LINESTYLE = 2, /OVERPLOT, /SYM_FILLED)
  p3.name = ' Bi-ex'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.4, 0.9])
  ;p.save, filename+'03.png', BORDER=10, RESOLUTION=300, /TRANSPARENT
  result05 = p.CopyWindow(BORDER=10, RESOLUTION=300)


  result_fig[*, size_single[1]:2*size_single[1]-1, size_single[2]:2*size_single[2]-1] = CONGRID(result04, 3, size_single[1], size_single[2])
  result_fig[*, 2*size_single[1]:3*size_single[1]-1, size_single[2]:2*size_single[2]-1] = CONGRID(result05, 3, size_single[1], size_single[2])
  IF N_ELEMENTS(sState) GT 0 THEN BEGIN
    new_pic = CONGRID(result_fig, 3, sState.xdim, sState.ydim)
    WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
    WSET, drawpicid
    TVSCL, new_pic, TRUE=1
  ENDIF

  ;write the whole plots
  fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(filename, '.txt')+'_fitplot.png', result_fig

  ;write the HDF5 file
  fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'hdf5\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  file = fitplot_dir+FILE_BASENAME(filename, '.txt')+'_plot_out.h5'
  fid = H5F_CREATE(file)
  ;; get data type and space, needed to create the dataset
  datatype_id = H5T_IDL_CREATE(result_fig)
  dataspace_id = H5S_CREATE_SIMPLE(size(result_fig,/DIMENSIONS))
  ;; create dataset in the output file
  dataset_id = H5D_CREATE(fid, 'dataset',datatype_id,dataspace_id)
  ;; write data to dataset
  H5D_WRITE,dataset_id,result_fig
  ;; close all open identifiers
  H5D_CLOSE,dataset_id
  H5S_CLOSE,dataspace_id
  H5T_CLOSE,datatype_id
  H5F_CLOSE,fid


  ;@GJ, 2022/6/9, save the result
  ;  picture = result_fig
  ;write file
  ;write file
  N_int = 8000.

  ;  ;@GJ, 2022/6/30 for 20220616131
  base_filename = FILE_BASENAME(filename, '.txt')
  filename_3peak_list = ['2022061602', '2022061606', '2022061607', '2022061610', '2022061619', '2022061631', '20220616098', '20220616099', '2022061637', '_2022072062', '_2022072104', '_2022072101']
  flag=0
  FOR i=0, N_ELEMENTS(filename_3peak_list)-1  DO BEGIN
    IF STRPOS(base_filename, filename_3peak_list[i]) NE -1 THEN BEGIN
      flag=1
      break
    ENDIF
  ENDFOR
  ;@GJ, 2022/9/1, always do the noise subtraction
  ;always do this process by substracting the noise
;  IF flag EQ 1 THEN BEGIN
    noise_n = MEAN(signal_array[N_period/2-1-N_period/50:N_period/2-1])
    y_data_ori = ABS(signal_array_n - noise_n)*1000.
;  ENDIF ELSE BEGIN
;    y_data_ori = -signal_array_n*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
;  ENDELSE
  X_data_int = FINDGEN(N_int+1)/N_int*MAX(X)
  y_data_int = signal_array_n
  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
  ;GJ, 2022/5/10, should not include the t=0 signal
  write_UPEN_001_file, filename, X_data_int[1:*], y_data_int[1:*], base_filename

  ;GJ, 2022/5/15, write the whole half-cycle signal
  write_UPEN_001_file, filename, time[0:N_period/2-1-N_period/50], -signal_array[0:N_period/2-1-N_period/50], base_filename+'half_cycle'

  file = fitplot_dir+FILE_BASENAME(filename, '.txt')+'_001_out.h5'
  fid = H5F_CREATE(file)
  ;; get data type and space, needed to create the dataset
  data = [X_data_int[1:*], y_data_int[1:*]]
  datatype_id = H5T_IDL_CREATE(data)
  dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))
  ;; create dataset in the output file
  dataset_id = H5D_CREATE(fid, 'dataset',datatype_id,dataspace_id)
  ;; write data to dataset
  H5D_WRITE,dataset_id,data
  ;; close all open identifiers
  H5D_CLOSE,dataset_id
  H5S_CLOSE,dataspace_id
  H5T_CLOSE,datatype_id
  H5F_CLOSE,fid

  ;GJ, 2022/6/16, close the plot window
  IF OBJ_VALID(p) THEN p.close

  ;  ;GJ, 2022/2/13, write the data for testing
  ;  filename = 'C:\D_drive\AIMIS_3D\LAPLACE\multiT2-master\multiT2-master\y_test2.dat'
  ;  openw,unit,filename,/get_lun
  ;  FOR i=0, N_ELEMENTS(Y)-1 DO BEGIN
  ;    ;; write the data points
  ;    WRITEU, unit, (ABS(signal_array_n)-MEAN(signal_array))[i]
  ;  ENDFOR
  ;  ;; close file and END execution
  ;  FREE_LUN, unit
  
  signal_rf_p_sP = DBLARR(times, N_period/2)
  FOR i=0, times-1 DO BEGIN
    signal_rf_p_sP[i, *] = REFORM(signal_array_sP[i, N_period/2:N_period-1]-MEAN(signal_array_sP[i, N_period-1-N_period/50:N_period-1]))
  ENDFOR
  signal_rf_p_1ave = REFORM(signal_rf_p_sP[0, *])
  signal_rf_p_1aveNew = REFORM(signal_rf_p_sP[1, *])
  signal_rf_p_1aveNew2 = REFORM(signal_rf_p_sP[100, *])
  signal_rf_p_1aveNew3 = REFORM(signal_rf_p_sP[200, *])
  signal_rf_p_1aveNew4 = REFORM(signal_rf_p_sP[300, *])
  signal_rf_p_1aveNew5 = REFORM(signal_rf_p_sP[399, *])
  signal_rf_p_10ave = REFORM(MEAN(signal_rf_p_sP[0:9, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_10aveNew = REFORM(MEAN(signal_rf_p_sP[10:19, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_50ave = REFORM(MEAN(signal_rf_p_sP[0:49, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_50aveNew = REFORM(MEAN(signal_rf_p_sP[50:99, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_100ave = REFORM(MEAN(signal_rf_p_sP[0:99, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_100aveNew = REFORM(MEAN(signal_rf_p_sP[100:199, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_200ave = REFORM(MEAN(signal_rf_p_sP[0:199, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_200aveNew = REFORM(MEAN(signal_rf_p_sP[200:399, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_300ave = REFORM(MEAN(signal_rf_p_sP[0:299, *], DIMENSION=1, /DOUBLE))
  signal_rf_p_400ave = REFORM(MEAN(signal_rf_p_sP[0:399, *], DIMENSION=1, /DOUBLE))
  
  fit_result = {$
    frequency:frequency, $
    time:time, $
    H:H_array, $
    H_peak_n:H_peak_n, $
    H_peak_p:H_peak_p, $
    M:M_array, $
    relax_H:relax_H, $
    relax_M:relax_M, $
    relax_mono_Mfit:relax_mono_Mfit, $
    relax_biex_Mfit:relax_biex_Mfit, $
    debye_H:debye_H, $
    debye_M:debye_M, $
    debye_mono_Mfit:debye_mono_Mfit, $
    tao_debye_mono:tao_debye_mono, $
    particle_size_mono:particle_size_mono, $ 
    signal:signal_array, $
    times: times, $
    signal_rf_p_sP: signal_rf_p_sP, $
    signal_rf_p_1ave: signal_rf_p_1ave, $
    signal_rf_p_10ave: signal_rf_p_10ave, $
    signal_rf_p_50ave: signal_rf_p_50ave, $
    signal_rf_p_100ave: signal_rf_p_100ave, $
    signal_rf_p_200ave: signal_rf_p_200ave, $
    signal_rf_p_300ave: signal_rf_p_300ave, $
    signal_rf_p_400ave: signal_rf_p_400ave, $
    signal_rf_p_1aveNew: signal_rf_p_1aveNew, $
    signal_rf_p_1aveNew2: signal_rf_p_1aveNew2, $
    signal_rf_p_1aveNew3: signal_rf_p_1aveNew3, $
    signal_rf_p_1aveNew4: signal_rf_p_1aveNew4, $
    signal_rf_p_1aveNew5: signal_rf_p_1aveNew5, $
    signal_rf_p_10aveNew: signal_rf_p_10aveNew, $
    signal_rf_p_50aveNew: signal_rf_p_50aveNew, $
    signal_rf_p_100aveNew: signal_rf_p_100aveNew, $
    signal_rf_p_200aveNew: signal_rf_p_200aveNew, $
    signal_p:signal_array_p, $ ;positive signal
    signal_n:signal_array_n, $ ;negative signal
    signal_rf_n:signal_array[0:N_period/2-1]-MEAN(signal_array[N_period/2-1-N_period/50:N_period/2-1]), $ ;signal from rise to flat
    signal_rf_p:signal_array[N_period/2:N_period-1]-MEAN(signal_array[N_period-1-N_period/50:N_period-1]), $ ;signal from rise to flat
    time_rise_flat:time[0:N_period/2-1], $ ;time during rise and flat phases
    time_flat: X, $ ;time for field flat phase
    time_rise: time_half, $ ; time for field rise phase
    M_diff_p:M_array_p_diff, $
    M_diff_n:M_array_n_diff, $
    delta_t:delta_t, $
    AUC_total_p: AUC_total_p,$
    AUC_total_n: AUC_total_n,$
    AUC_slew_p: AUC_slew_p, $
    AUC_slew_n: AUC_slew_n, $
    AUC_hold_p: AUC_hold_p,$
    AUC_hold_n: AUC_hold_n,$
    AUC_slew: AUC_slew, $
    AUC_hold: AUC_hold, $
    AUC_total: AUC_total, $
    AUC_slew_ratio: AUC_slew_ratio, $
    spec:spec_array, $
    corre_biex: corre_biex, $
    corre_Debye_Mono: corre_Debye_Mono, $
    tao:tao, $
    bi_tao:[bi_tao_1, bi_tao_2], $
    bi_tao_comp: [bi_tao_1_comp, bi_tao_2_comp], $
    bi_tao_combine: bi_tao_combine, $
    mono_para:A, bi_para:B $
    }
  
  ;close the window  
  IF OBJ_VALID(p) THEN p.close
  
  RETURN, fit_result;[[H_array], [signal_array], [spec_array]]
END

PRO read_data_H5
  ; Open the HDF5 file.
  file = 'C:\D_drive\MPI_Tianjie\Xinfeng\20220524_T_5mT_90flat\20220524\20220524\hdf5\2022052434_plot_out.h5'

  file_id = H5F_OPEN(file)

  ; Open the image dataset within the file.
  ; This is located within the /images group.
  ; We could also have used H5G_OPEN to open up the group first.
  dataset_id = H5D_OPEN(file_id, 'dataset')
  ; Read in the actual image data.
  image = H5D_READ(dataset_id1)
  H5S_CLOSE, dataset_id
  H5F_CLOSE, file_id

  ; Display the data.
  xdim=1200
  ydim=600
  DEVICE, DECOMPOSED=0
  WINDOW, 1, xsize=xdim, ysize=ydim
  new_pic = CONGRID(image, 3, xdim, ydim)
  TVSCL, new_pic, TRUE=1

END

;GJ, 2022/4/14, write and read 001 files for multi-color MPI
PRO write_UPEN_001_file, filename, X_data_int, y_data_int, base_filename
  
  ;GJ, 2022/2/15, for UPEN program
  IF N_ELEMENTS(base_filename) NE 0 THEN BEGIN
    ;filename = 'C:\UpenWin\bin\MPI01\test\'+base_filename+'.001'
    UPEN001_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'UPEN001\'
    IF FILE_TEST(UPEN001_dir, /directory) LT 1 THEN FILE_MKDIR, UPEN001_dir
    filename001 = UPEN001_dir+base_filename+'.001'
  ENDIF ELSE BEGIN
    filename001 = filename+'.001'
  ENDELSE
  openw,unit1,filename001,/get_lun
  first_line =['time', 'Signal']
  printf, unit1, FORMAT = '(2(A, %"  "))', first_line
  FOR i=0, N_elements(X_data_int)-1 DO BEGIN
    ;; write the data points
    printf, unit1, FORMAT = '(2(A, %"  "))', X_data_int[i], y_data_int[i]
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit1
  
  print, 'finish writing a file:'+filename001
END

;@GJ, 2022/4/14, read UPEN 001 file for multi-color MPI
;data_array = read_UPEN_001_file()
;data_array is N*2 array
FUNCTION read_UPEN_001_file, filename
  
  IF N_ELEMENTS(filename) EQ 0 THEN filename = dialog_pickfile(TITLE='Select UPEN 001 File', FILTER='*.001', /MUST_EXIST, PATH='C:\UpenWin\bin\MPI01\')
  
  OPENR, lun, filename, /get_lun

  ;read the first 1042 lines, extract the position of log coordinates
  temp_str=''
  readf, lun, temp_str
  
  ;define the array
  N_int = 10000.
  data_x = DBLARR(N_int)*0.
  data_y = DBLARR(N_int)*0.
  j=0
  WHILE ~ EOF(lun) DO BEGIN
    READF, lun, temp1, temp2, FORMAT='(%"%f %f")'
    data_x[j] = temp1
    data_y[j] = temp2
    j++
  ENDWHILE
  FREE_LUN, lun
  
  data_array = [[data_x[0:j-1]], [data_y[0:j-1]]]

  RETURN, data_array
END

;GJ, 2022/4/14, to mix different particles, viscosity, temperature
PRO multi_color_MPI
  
  ;load color1
  data_array1 = read_UPEN_001_file()
  plot, REFORM(data_array1[*,0]), REFORM(data_array1[*,1])
  
  ;load color2
  data_array2 = read_UPEN_001_file()
  data_array2[*,1] = data_array2[*,1]/MAX(data_array2[*,1])*MAX(data_array1[*,1])
  oplot, REFORM(data_array2[*,0]), REFORM(data_array2[*,1])
  
  ;mix them
  mix_data_array = (data_array1+data_array2)/2.
  oplot, REFORM(mix_data_array[*,0]), REFORM(mix_data_array[*,1])

  output_filename = 'C:\UpenWin\bin\MPI01\mixed.001'
  write_UPEN_001_file, output_filename, REFORM(mix_data_array[*,0]), REFORM(mix_data_array[*,1])
  
   
END

PRO gfunct, X, A, F, pder

  bx = EXP(A[1] * X)
  F = A[0] * bx + A[2]

  ;If the procedure is called with four parameters, calculate the
  ;partial derivatives.

  IF N_PARAMS() GE 4 THEN $
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]

END

;@GJ, 2022/8/20, A[0]*EXP(-A[2] * X)+A[1]*EXP(-A[3] * X)+A[4]
PRO gfunct_biexdecay, X, A, F, pder

  bx1 = EXP(-A[2] * X)
  F1 = A[0] * bx1
  
  bx2 = EXP(-A[3] * X)
  F2 = A[1] * bx2
  
  F = F1 + F2 + A[4]
  ;If the procedure is called with four parameters, calculate the
  ;partial derivatives.

  IF N_PARAMS() GE 4 THEN $
    pder = [[bx1], [bx2], [-A[0] * X * bx1], [-A[1] * X * bx2], [replicate(1.0, N_ELEMENTS(X))]]

END

;A[0] = alpha
;A[1] = beta
;define the Langevin function
PRO langevin_funct, X, A, F, pder
  
  
  L = 1./TANH(A[1]*X) - 1./(A[1]*X)
  L_prime = 1./(A[1]*X)^2 - 1./(SINH(A[1]*X)^2)
  
  F = A[0] * L

  ;If the procedure is called with four parameters, calculate the
  ;partial derivatives.

  IF N_PARAMS() GE 4 THEN $
    pder = [[L], [A[0] * X * L_prime]]

END

PRO gfunct_biex, X, A, F, pder

  bx_1 = EXP(A[1] * (X-A[5]))
  bx_2 = EXP(A[3] * (X-A[5]))
  F = A[0] * bx_1 + A[2] * bx_2 + A[4]

  ;If the procedure is called with four parameters, calculate the
  ;partial derivatives.

  IF N_PARAMS() GE 4 THEN $
    pder = [[bx_1], [A[0] * X * bx_1], [bx_2], [A[2] * X * bx_2], [replicate(1.0, N_ELEMENTS(X))], [A[0] * A[1] * bx_1 + A[2] * A[3] * bx_2]]

END


PRO gfunct_triex, X, A, F, pder

  bx_1 = EXP(A[1] * X)
  bx_2 = EXP(A[3] * X)
  bx_3 = EXP(A[5] * X)
  F = A[0] * bx_1 + A[2] * bx_2 + A[4] * bx_3 + A[6]

  ;If the procedure is called with four parameters, calculate the
  ;partial derivatives.

  IF N_PARAMS() GE 4 THEN $
    pder = [[bx_1], [A[0] * X * bx_1], [bx_2], [A[2] * X * bx_2], [bx_3], [A[4] * X * bx_3], [replicate(1.0, N_ELEMENTS(X))]]

END


PRO iLaplace
  
  N_period = 1000
  delta_t  = 1.0
  t = FINDGEN(N_period)/10.*delta_t
  tao = 6.4
  y = 5. * EXP(t / tao)
;  y = t
  plot, t, y
  
  F = DBLARR(N_period)
  s = (FINDGEN(N_period)+1.)/N_period/delta_t
  FOR i=0, N_period-1 DO IF s[i] GT 1./tao THEN F[i] = TOTAL(y * EXP(-s[i] * t) * delta_t)
  plot, s, F, xrange=[0,1]

END

;GJ, 2022/4/7, generate .001 file for Upen based signal_b
PRO batch_generate_UPEN_001_files_ADC

  dir_name = 'C:\D_drive\MPI_Tianjie\T1rho_MPI\b-valuevssignal\'
;  dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)

  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN

  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.*', count=nrfile)
  FOR i=0L, nrfile-1L DO BEGIN
    filename_base = FILE_BASENAME(a_file[i], '.001')
    fit_para = read_b_value_generate_UPEN_001_file(a_file[i])
    print, 'i = ', i+1, ', filename = ', filename_base
  ENDFOR

END

;GJ, 2022/4/7, read signal vs b value
FUNCTION read_b_value_generate_UPEN_001_file, filename
  
;  filename = 'C:\D_drive\MPI_Tianjie\T1rho_MPI\b-valuevssignal\V20_scan1.001'

  
  ;GJ, 2022/3/8, nonzero gilename
  IF STRLEN(filename) GT 1 THEN BEGIN
    filename_info = FILE_INFO(filename)
    IF (filename_info.exists EQ 0) THEN BEGIN
      info=dialog_message('Wrong file name!', /information)
      RETURN, 0
    ENDIF
    IF (STRLEN(filename) EQ 0) THEN BEGIN
      info=dialog_message('No file name!', /information)
      RETURN, 0
    ENDIF
  ENDIF
  
  
  ;read the txt file 1st time to get the frequency and period
  OPENR, unit, filename, /GET_LUN
  header=''
  READF, unit, header, FORMAT='(%"%s")'
  i=0;
  WHILE ~ EOF(unit) DO BEGIN
    READF, unit, a_temp, b_temp, FORMAT='(%"%f  %f")'
    IF i EQ 0 THEN BEGIN
      b_array = [a_temp]
      signal_array = [b_temp]
    ENDIF ELSE BEGIN
      b_array = [b_array, a_temp]
      signal_array = [signal_array, b_temp]
    ENDELSE
    i++
  ENDWHILE
  ;; close file and END execution
  FREE_LUN, unit
  
  base_filename = FILE_BASENAME(filename, '.txt')
  print, 'b_array: ', b_array
  iplot, b_array, signal_array, title='Signal vs b value ('+base_filename+')'
  
  ;write file
  N_int = 8000.
  y_data_ori = signal_array;*1000. ;(ABS(signal_array_n)-MEAN(signal_array))*100000.
  X = b_array
  X_data_int = FINDGEN(N_int+1)/N_int*MAX(X)
  y_data_int = INTERPOL(y_data_ori, X, X_data_int)
  ;GJ, 2022/2/15, for UPEN program
  output_filename = 'C:\D_drive\UpenWin\bin\ADC_b\'+base_filename+'.001'
  openw,unit1,output_filename,/get_lun
  first_line =['b_value', 'Signal']
  printf, unit1, FORMAT = '(2(A, %"  "))', first_line
  FOR i=5, N_int DO BEGIN
    ;; write the data points
    printf, unit1, FORMAT = '(2(A, %"  "))', X_data_int[i], y_data_int[i]
  ENDFOR
  ;; close file and END execution
  FREE_LUN, unit1
  
  RETURN, 1
END

;@GJ, 2022/4/16, RGB represents AUC, Neel, and Brownian relaxation times respectively
PRO multi_color_phantom, mc_phantom_image, square_phantom_image
  ;read the black-white image
  filename_bmp = 'C:\D_drive\MPI_Tianjie\FID_MPI_iLT\dp_bw.bmp'
  image_bmp = READ_BMP(filename_bmp)
  d_phantom = REFORM(image_bmp[0,*,*])
;  iimage, d_phantom
  size_dp = SIZE(d_phantom, /dim)
  d_phantom_sq = CONGRID(d_phantom[size_dp[0]-size_dp[1]:*, *], 128, 128)
  ;extract catheter part
  catheter_index = WHERE(d_phantom_sq GT 240, catheter_count)
  catheter_index_2d = ARRAY_INDICES(d_phantom_sq, catheter_index)
  catheter_image = d_phantom_sq*0.
  catheter_image[catheter_index] = 1
  
  ;define the phantom image
  phantom_image = image_bmp[0:2, *, *] * 0.
  
  ;extract catheter part
  catheter_index = WHERE(d_phantom GT 240, catheter_count)
  d_phantom[catheter_index] = 0.
  catheter_part = d_phantom*0.
  catheter_part[catheter_index] = 1
;  iimage, catheter_part, title='Catheter'
  para = [10., 10., 25.]; [AUC, tau_Neel, tau_Brownian]
  phantom_image[0, *, *] += catheter_part*para[0] ;AUC
  phantom_image[1, *, *] += catheter_part*para[1] ;tau_Neel
  phantom_image[2, *, *] += catheter_part*para[2] ;tau_Brownian

  ;extract plaque part
  plaque_index = WHERE(d_phantom GT 200, plaque_count)
  IF plaque_count NE 0 THEN BEGIN
    d_phantom[plaque_index] = 0.
    plaque_temp = d_phantom*0.
    plaque_temp[plaque_index] = 1
    labelImg = LABEL_REGION(plaque_temp, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    plaque_part = FLOAT(labelImg EQ sort_re[1])
;    iimage, plaque_part, title='Plaque'
    para = [10., 7., 60.]; [AUC, tau_Neel, tau_Brownian]
    phantom_image[0, *, *] += plaque_part*para[0] ;AUC
    phantom_image[1, *, *] += plaque_part*para[1] ;tau_Neel
    phantom_image[2, *, *] += plaque_part*para[2] ;tau_Brownian
  ENDIF
 
  ;extract vessel part
  vessel_index = WHERE(d_phantom GT 160, vessel_count)
  IF vessel_count NE 0 THEN BEGIN
    d_phantom[vessel_index] = 0.
    vessel_temp = d_phantom*0.
    vessel_temp[vessel_index] = 1
    labelImg = LABEL_REGION(vessel_temp, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    vessel_part = FLOAT(labelImg EQ sort_re[1])
;    iimage, vessel_part, title='Vessel'
    para = [10., 7., 15.]; [AUC, tau_Neel, tau_Brownian]
    phantom_image[0, *, *] += vessel_part*para[0] ;AUC
    phantom_image[1, *, *] += vessel_part*para[1] ;tau_Neel
    phantom_image[2, *, *] += vessel_part*para[2] ;tau_Brownian
  ENDIF
  
  ;extract tumorRimUpper part
  tumorRimUpper_index = WHERE(d_phantom GT 125, tumorRimUpper_count)
  IF tumorRimUpper_count NE 0 THEN BEGIN
    d_phantom[tumorRimUpper_index] = 0.
    tumorRimUpper_temp = d_phantom*0.
    tumorRimUpper_temp[tumorRimUpper_index] = 1
    labelImg = LABEL_REGION(tumorRimUpper_temp, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    tumorRimUpper_part = FLOAT(labelImg EQ sort_re[1])
;    iimage, tumorRimUpper_part, title='tumorRimUpper'
    para = [15., 7., 15.]; [AUC, tau_Neel, tau_Brownian]
    phantom_image[0, *, *] += tumorRimUpper_part*para[0] ;AUC
    phantom_image[1, *, *] += tumorRimUpper_part*para[1] ;tau_Neel
    phantom_image[2, *, *] += tumorRimUpper_part*para[2] ;tau_Brownian
  ENDIF
  
  ;extract tumorRimLower part
  tumorRimLower_index = WHERE(d_phantom GT 85, tumorRimLower_count)
  IF tumorRimLower_count NE 0 THEN BEGIN
    d_phantom[tumorRimLower_index] = 0.
    tumorRimLower_temp = d_phantom*0.
    tumorRimLower_temp[tumorRimLower_index] = 1
    labelImg = LABEL_REGION(tumorRimLower_temp, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    tumorRimLower_part = FLOAT(labelImg EQ sort_re[1])
;    iimage, tumorRimLower_part, title='tumorRimLower'
    para = [15., 3., 8.]; [AUC, tau_Neel, tau_Brownian]
    phantom_image[0, *, *] += tumorRimLower_part*para[0] ;AUC
    phantom_image[1, *, *] += tumorRimLower_part*para[1] ;tau_Neel
    phantom_image[2, *, *] += tumorRimLower_part*para[2] ;tau_Brownian
  ENDIF
  
  ;extract tumorCenterLower part
  tumorCenterLower_index = WHERE(d_phantom GT 45, tumorCenterLower_count)
  IF tumorCenterLower_count NE 0 THEN BEGIN
    d_phantom[tumorCenterLower_index] = 0.
    tumorCenterLower_temp = d_phantom*0.
    tumorCenterLower_temp[tumorCenterLower_index] = 1
    labelImg = LABEL_REGION(tumorCenterLower_temp, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    tumorCenterLower_part = FLOAT(labelImg EQ sort_re[1])
;    iimage, tumorCenterLower_part, title='tumorCenterLower'
    para = [5., 3., 8.]; [AUC, tau_Neel, tau_Brownian]
    phantom_image[0, *, *] += tumorCenterLower_part*para[0] ;AUC
    phantom_image[1, *, *] += tumorCenterLower_part*para[1] ;tau_Neel
    phantom_image[2, *, *] += tumorCenterLower_part*para[2] ;tau_Brownian
  ENDIF

  ;extract tumorCenterUpper part
  tumorCenterUpper_index = WHERE(d_phantom GT 0, tumorCenterUpper_count)
  IF tumorCenterUpper_count NE 0 THEN BEGIN
    d_phantom[tumorCenterUpper_index] = 0.
    tumorCenterUpper_temp = d_phantom*0.
    tumorCenterUpper_temp[tumorCenterUpper_index] = 1
    labelImg = LABEL_REGION(tumorCenterUpper_temp, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    tumorCenterUpper_part = FLOAT(labelImg EQ sort_re[1])
;    iimage, tumorCenterUpper_part, title='tumorCenterUpper'
    para = [5., 7., 15.]; [AUC, tau_Neel, tau_Brownian]
    phantom_image[0, *, *] += tumorCenterUpper_part*para[0] ;AUC
    phantom_image[1, *, *] += tumorCenterUpper_part*para[1] ;tau_Neel
    phantom_image[2, *, *] += tumorCenterUpper_part*para[2] ;tau_Brownian
  ENDIF
  
  ;remove the gap pixels between two parts
  size_ima = SIZE(phantom_image, /dim)
  mc_phantom_image = phantom_image[0:2, 2:(size_ima[1]-5), 5:(size_ima[2]-3)]
  size_mc = SIZE(mc_phantom_image, /dim)
  ;skip the edge pixels
  FOR j=2, size_mc[1]-3 DO BEGIN
    FOR k=2, size_mc[2]-3 DO BEGIN
      IF mc_phantom_image[0, j, k] EQ 0 THEN BEGIN
        mc_phantom_image[0, j, k] = MEDIAN(mc_phantom_image[0, j-1:j+1, k-1:k+1])
        mc_phantom_image[1, j, k] = MEDIAN(mc_phantom_image[1, j-1:j+1, k-1:k+1])
        mc_phantom_image[2, j, k] = MEDIAN(mc_phantom_image[2, j-1:j+1, k-1:k+1])
      ENDIF
    ENDFOR
  ENDFOR
;  iimage, mc_phantom_image, title='Multi-color Phantom Image'
  
  ;crop the image as a squre matrix
  square_phantom_image = mc_phantom_image[0:2, size_mc[1]-size_mc[2]:size_mc[1]-1, 0:size_mc[2]-1]
;  iimage, square_phantom_image, title='Multi-color Phantom Image'
  
;  iimage, REFORM(square_phantom_image[0,*,*]), title='Concentration Phantom Image'
END

;@GJ, 2022/4/17, simulate the signal
PRO multi_color_signal_generation, test_image, rec_image, correlation_coe, MSE_coe, RMSE_coe, noise_level
  
  ;generate the image
  multi_color_phantom, mc_phantom_image, square_phantom_image
  size_sqmc = SIZE(square_phantom_image, /dim)
  
  matrix_dim = 128
  test_image = CONGRID(mc_phantom_image, size_sqmc[0], matrix_dim, matrix_dim)
  iimage, test_image, title='Test Phantom Image ('+STRING(matrix_dim)+'*'+STRING(matrix_dim)+')'
  
  ;reconstructed image
  rec_image = test_image*0.
  
  ;generate the relaxation-induced decay signal pixel by pixel
  N_int = 2500. ; in 0.1 us
  data_x = (FINDGEN(N_int)+5.)/10.
  delta_x = data_x[5]-data_x[4]
  data_y = DBLARR(N_int)*0.
  FOR j=0, matrix_dim-1 DO BEGIN
    FOR k=0, matrix_dim-1 DO BEGIN
      IF test_image[0, j, k] NE 0 THEN BEGIN
        ;[AUC, tau_Neel, tau_Brownian]
        AUC = test_image[0, j, k]
        tau_N = test_image[1, j, k]
        tau_B = test_image[2, j, k]
        IF tau_B GT 0 THEN BEGIN
          temp_data_y = (1./tau_N*EXP(-data_x/tau_N) + 1./tau_B*EXP(-data_x/tau_B))
          data_y = temp_data_y/(TOTAL(temp_data_y)*delta_x)*AUC
        ENDIF ELSE BEGIN
          temp_data_y = (1./tau_N*EXP(-data_x/tau_N))
          data_y = temp_data_y/(TOTAL(temp_data_y)*delta_x)*AUC
        ENDELSE
        
        IF N_ELEMENTS(noise_level) EQ 0 THEN noise_level=0.1
        
        ;add noise to data_y
        data_y_noise = data_y*0.
        FOR i=0, N_ELEMENTS(data_y)-1 DO BEGIN
          seed = !NULL
          randomValue = RANDOMU(seed)*noise_level
          data_y[i] = data_y[i]*(1.+randomValue)
        ENDFOR
        
;        plot, data_x, data_y
        print, 'j =', j, ', k=', k
;        wait, 1
;        print, 'tao_Neel: ', tau_N, '[us]'
;        print, 'tao_Brownian: ', tau_B, '[us]'
        
        ;start the reconstruction
        X = data_x
        Y = data_y
        ;do the biexponential curve fit
        ;Provide an initial guess of the function’s parameters.
        B = [MAX(Y)/2.0, -0.01, MAX(Y)/2.0, -0.1, MIN(Y), 0.]
        FITA = [1,1,1,1,1,0]
        ;Compute the parameters.
        yfit_biex = CURVEFIT(X, Y, weights, B, SIGMA, FITA=FITA, FUNCTION_NAME='gfunct_biex')
        corre_biex = CORRELATE(y, yfit_biex)
;        PRINT, 'correlation coeff (Bi-ex): ', corre_biex

;        ;GJ, 2022/2/11, plot the fitting
;        p = PLOT(X, Y, 'o-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
;          XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], title='Exponential Relaxation, '+STRTRIM(STRING(2000),1)+' Hz', xtitle='time [us]', ytitle='M [A.U.]')
;        p.name = ' Original data'
;        p2 = PLOT(X, yfit_biex, '*r-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
;        p2.name = ' Bi-ex fitting'
;        l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.5])

        ;GJ, 2022/2/10, do correctin to the coefficients for correction Neel composition
        B[0] = B[0]/ABS(B[1])
        B[2] = B[2]/ABS(B[3])
        
        bi_tao_1 = 1./ABS(B[1])
        bi_tao_1_comp = ABS(B[0]/(B[0]+B[2]))
        bi_tao_2 = 1./ABS(B[3])
        bi_tao_2_comp = ABS(B[2]/(B[0]+B[2]))
;        print, 'Relaxation time_1 = ', bi_tao_1, ' [us] with',  bi_tao_1_comp*100., '%'
;        print, 'Relaxation time_2 = ', bi_tao_2, ' [us] with',  bi_tao_2_comp*100., '%'
;        wait, 1
        
        ;set the reconstructed result
        rec_image[0, j, k] = TOTAL(data_y)*delta_x
        rec_image[1, j, k] = bi_tao_2
        rec_image[2, j, k] = bi_tao_1
      ENDIF
    ENDFOR
  ENDFOR
  
  iimage, rec_image, title='reconstrcted image'
  
  correlation_coe = DBLARR(3)*0.
  MSE_coe = DBLARR(3)*0.
  RMSE_coe = DBLARR(3)*0.
  FOR i=0, 2 DO BEGIN
    correlation_coe[i] = CORRELATE(REFORM(test_image[i,*,*]), REFORM(rec_image[i,*,*]))
    MSE_coe[i] = calMSE(REFORM(test_image[i,*,*]), REFORM(rec_image[i,*,*]))
    RMSE_coe[i] = Sqrt(Total((REFORM(test_image[i,*,*]) - REFORM(rec_image[i,*,*]))^2)/N_Elements(REFORM(rec_image[i,*,*])))
  ENDFOR
  
  print, correlation_coe
  print, MSE_coe
  print, RMSE_coe
;  Image_Blend, BYTSCL(REFORM(rec_image[0,*,*])), BYTSCL(REFORM(rec_image[1,*,*]), MAX=12), Colortable=25;, blendTitle='Neel'
  
;  Image_Blend, BYTSCL(REFORM(rec_image[0,*,*])), BYTSCL(REFORM(rec_image[2,*,*]), MAX=30), Colortable=25
END

;@LXF, 2022/7/30, code for ART algebraic reconstruction technique (ART)
function art_func,A,b

  n=N_Elements(b)
  X0=FltArr(1,n)
  ;X0=transpose(X0)
  ;print,X0
  X=X0
  k=0
  k1=1000;1000
  e=2
  e0=0.001
  ;while(SQRT(e) GE e0) do begin
  while(k le k1) do begin
    for i=0,n-1 do begin
      unitLen = norm(A[*,i],/double)
      ;unitLen = Total(A[*,i]##transpose(A[*,i]))
      ;unitLen = SQRT(A[*,i]##transpose(A[*,i]))
      ;print,unitLen
      d = (b[i]-A[*,i]##X)/unitLen
      ;print,d
      Xf = X
      X = X + (d##(A[*,i]/unitLen))
      ;print,A[*,i]/unitLen
      ;print,d*transpose(A[*,i]/unitLen)
      ;e = norm(Xf-X)
    endfor
    k=k+1
  endwhile
  print,k
  return,X
end


;朗之万方程求导
PRO Langevin_derivative_Function
  x1=FINDGEN(101)/10*2.+0.1

  ;PLOT, x1, 1./(TANH(x1))-1./x1

  x2=(FINDGEN(101)/10-10.)*2.-0.1
  ;PLOT, x2, 1./(TANH(x2))-1./x2

  x12=[x1,x2]
  x=x12[sort(x12)]
  iplot, x, 1./(TANH(x))-1./x, PSYM=3, title='Langevin';, yrange=[-15,15]
  ;  plot, x, 1./(TANH(x))
  ;  iplot, x, -1./x, PSYM=2, yrange=[-15,15]
  iplot, x, 1./(x*x)-1./(SINH(x)*SINH(x)), PSYM=3, title='Langevin Derivative';, yrange=[-15,15]


END

;@GJ, 2023/5/4, gradient-sweeping MPI for high-resolution MPI
PRO Langevin_derivative_PSF_estimation, particle_size, H_max, H_array, PSF, H_array_100mT, PSF_100mT

  IF N_ELEMENTS(particle_size) EQ 0 THEN particle_size = 20. ;nm
  particle_size *= 1.0

  n_H = 600.
  H_array = FINDGEN(n_H+1)/(n_H/2.) - 1.
  IF N_ELEMENTS(H_max) EQ 0 THEN H_max = 20.
  H_array *= H_max

  n_H_100mT = 10000.
  H_array_100mT = FINDGEN(n_H_100mT+1)/(n_H_100mT/2.) - 1.
  H_array_100mT *= (H_max * 10.) ; 100 mT

  Msat_T = 0.551 ; T/u0
  u0 = 4.0 * !PI * (0.1^4) ; T*m/kA
  Msat_kAm = Msat_T/u0; kA/m
  m_moment = 1./6. * !PI * (particle_size^3) * Msat_kAm ; 10^-27 kA*m2
  k_B = 1.380469 ; 10^-23 J/K
  T_p = 36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = m_moment / (K_B * T_p_kelvin) / (10.^4) ; in kA*m2/J or 1/mT
  M = Msat_kAm * (1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
  M_100mT = Msat_kAm * (1./TANH(beta_particle*H_array_100mT) - 1./(beta_particle*H_array_100mT))
  M[n_H/2] = (M[n_H/2-1]+M[n_H/2+1])/2.
  M_100mT[n_H_100mT/2] = (M_100mT[n_H_100mT/2-1]+M_100mT[n_H_100mT/2+1])/2.
  PSF = M[1:*] - M[0:n_H]
  PSF_100mT = M_100mT[1:*] - M_100mT[0:n_H_100mT]
;    iplot, H_array, M, PSYM=3, XTITLE='H [mT]', YTITLE='M', XRANGE=[-H_max,H_max], TITLE='Langevin Function '+STRTRIM(STRING(particle_size),1)+' nm'
;    iplot, H_array[0:n_H-1], PSF, PSYM=3, XTITLE='H [mT]', YTITLE='PSF', XRANGE=[-H_max,H_max], TITLE='PSF '+STRTRIM(STRING(particle_size),1)+' nm'
;    iplot, H_array_100mT[0:n_H_100mT-1], PSF_100mT
END

;@GJ, 2023/8/9, 2014 Nobel Prize on fluorescence imaging, resolution enhancement based on sequential imaging (RESI)
PRO RESI_donut_pattern
  
  particle_size = 30. ;nm
  H_max = [10., 12.] ; mT
  Langevin_derivative_PSF_estimation, particle_size, H_max[0], H_array1, PSF1, H_array_100mT1, PSF_100mT1
  Langevin_derivative_PSF_estimation, particle_size, H_max[1], H_array2, PSF2, H_array_100mT2, PSF_100mT2
  n_H_100mT1 = N_ELEMENTS(H_array_100mT1)
  n_H_100mT2 = N_ELEMENTS(H_array_100mT2)
  gradient = 2.; mT/mm
  gradient2 = gradient / MAX(H_array_100mT1) * MAX(H_array_100mT2)
  print, 'gradient: ', gradient, gradient2, 'T/m'

  iplot, H_array_100mT1[0:n_H_100mT1-2]/gradient, PSF_100mT1/MAX(PSF_100mT1), color='blue', xtitle='x [mm]', ytitle='PSF', xrange=[-5.,5.], title='PSF'
  iplot, H_array_100mT2[0:n_H_100mT2-2]/gradient2, PSF_100mT2/MAX(PSF_100mT2), color='green', /overplot
  
  result1 = INTERPOL(H_array_100mT1[0:n_H_100mT1/2]/gradient, PSF_100mT1[0:n_H_100mT1/2]/MAX(PSF_100mT1), 0.5)
  FWHM1 = ABS(result1) * 2.
  result2 = INTERPOL(H_array_100mT2[0:n_H_100mT2/2]/gradient2, PSF_100mT2[0:n_H_100mT2/2]/MAX(PSF_100mT2), 0.5)
  FWHM2 = ABS(result2) * 2.
  
  PSF_diff = PSF_100mT1/MAX(PSF_100mT1) - PSF_100mT2/MAX(PSF_100mT2)
  iplot, H_array_100mT1[0:n_H_100mT2-2]/gradient, PSF_diff, color='cyan', thick=2, /overplot
  
;  Langevin_derivative_PSF_estimation, particle_size, H_max[0], H_array0, PSF0, H_array_100mT0, PSF_100mT0
;  iplot, H_array_100mT1[0:n_H_100mT2-2], PSF_100mT0/MAX(PSF_100mT0), color='pink', /overplot
  
;  iplot, PSF_diff, PSF_100mT0
  
  maxPD = MAX(PSF_diff[0:n_H_100mT1/2], maxind)
  maxPD2 = MAX(PSF_diff[n_H_100mT1/2:*], maxind2)
  maxPD_PSF_100mT2 = PSF_100mT2[maxind]/MAX(PSF_100mT2)
  
  tao = -maxPD / ALOG(0.001)
;  tao = -maxPD / ALOG(maxPD_PSF_100mT2)
  new_PSF_diff = exp(-PSF_diff/tao)
;  iplot, psf_diff, new_psf_diff, title='mapping relationship'
  new_PSF_diff[0:maxind] = 0;PSF_diff[0:maxind]
  new_PSF_diff[n_H_100mT1/2+maxind2:*] = 0;PSF_diff[n_H_100mT1/2+maxind2:*]
  iplot, H_array_100mT1[0:n_H_100mT2-2]/gradient, new_PSF_diff, color='red', /overplot
  
  new_result = INTERPOL(H_array_100mT1[0:n_H_100mT2/2]/gradient, new_PSF_diff[0:n_H_100mT2/2], 0.5)
  new_FWHM = ABS(new_result) * 2.
  print, 'FWHM: ', FWHM1, FWHM2, new_FWHM, 'mm'
  
  
  
;  iplot, H_array_100mT1[0:n_H_100mT2-2], new_PSF_diff-PSF_100mT2/MAX(PSF_100mT2), color='black', /overplot
END


;@GJ, 2023/5/4, Gradient sweeping MPI
;@GJ, 2023/8/10, Resolution ehancement based on two gradients
PRO RESI_gradient_sweeping_MPI_one_dimension

  N_H_max = 10.
  H_max_array = FINDGEN(N_H_max)+10.
  H_max = H_max_array[0]
  Langevin_derivative_PSF_estimation, particle_size, H_max, H_array, PSF0, H_array_100mT, PSF_100mT
  n_H = N_ELEMENTS(H_array)
  PSF_array = FINDGEN(N_H_max, n_H-1)
  n_H_100mT = N_ELEMENTS(H_array_100mT)
  PSF_array_100mT = DBLARR(N_H_max, n_H_100mT-1)*0.
  H_array_H_max = DBLARR(N_H_max, n_H)*0.
  H_array_H_max_100mT = DBLARR(N_H_max, n_H_100mT)*0.
  FOR i=0, N_H_max-1 DO BEGIN
    Langevin_derivative_PSF_estimation, particle_size, H_max_array[i], H_array, PSF, H_array_100mT, PSF_100mT
    PSF_array[i,*] = PSF;/MAX(PSF)
    PSF_array_100mT[i,*] = PSF_100mT;/MAX(PSF_100mT)
    H_array_H_max[i,*] = H_array
    H_array_H_max_100mT[i,*] = H_array_100mT
  ENDFOR

  ;we assume FOV 4 mm
  FOV = 40. ;mm
  gradient_field_array = H_max_array / (FOV/2.) ;3.0e3 ;mT/mm
  print, 'gradient field array: ', gradient_field_array
  print, 'H max array: ', H_max_array

  window, 0
  plot, H_array_H_max[N_H_max-1,0:n_H-1]/gradient_field_array[N_H_max-1], REFORM(PSF_array[N_H_max-1,*]), XRANGE=[-FOV/2., FOV/2.], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Location [mm]', ytitle='PSF', TITLE='PSF '+STRTRIM(STRING(particle_size),1)+' nm'
  FOR i=1, N_H_max-1 DO oplot, H_array_H_max[N_H_max-1-i, 0:n_H-1]/gradient_field_array[N_H_max-1-i], REFORM(PSF_array[N_H_max-1-i,*]), LINESTYLE=2, color=3500
  ;  oplot, H_array[0:n_H-1], PSF3, LINESTYLE=3, color=6000

  window, 1
  plot, H_array_H_max_100mT[N_H_max-1,0:n_H_100mT-1]/gradient_field_array[N_H_max-1], REFORM(PSF_array_100mT[N_H_max-1,*]), XRANGE=[-FOV/2., FOV/2.], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Location [mm]', ytitle='PSF 100 mT', TITLE='PSF 100 mT (Particle Size: '+STRTRIM(STRING(particle_size),1)+' nm)'
  FOR i=1, N_H_max-1 DO oplot, H_array_H_max_100mT[N_H_max-1-i, 0:n_H_100mT-1]/gradient_field_array[N_H_max-1-i], REFORM(PSF_array_100mT[N_H_max-1-i,*]), LINESTYLE=2, color=3500
  
  ;@GJ, 2023/8/10, analyzing the diff
  PSF_diff = REFORM(PSF_array_100mT[N_H_max-1-9,*])/MAX(PSF_array_100mT[N_H_max-1-9,*])-REFORM(PSF_array_100mT[N_H_max-1,*])/MAX(PSF_array_100mT[N_H_max-1,*])
  oplot, H_array_H_max_100mT[N_H_max-1, 0:n_H_100mT-1]/gradient_field_array[N_H_max-1], PSF_diff, color='0000x'
  
  maxPD = MAX(PSF_diff[0:n_H_100mT/2], maxind)
  maxPD2 = MAX(PSF_diff[n_H_100mT/2:*], maxind2)
  
  tao = -maxPD / ALOG(0.001)
  new_PSF_diff = exp(-PSF_diff/tao)
  new_PSF_diff[0:maxind] = 0;PSF_diff[0:maxind]
  new_PSF_diff[n_H_100mT/2+maxind2:*] = 0;PSF_diff[n_H_100mT1/2+maxind2:*]
  oplot, H_array_H_max_100mT[N_H_max-1, 0:n_H_100mT-1]/gradient_field_array[N_H_max-1], new_PSF_diff, color='0000x'

  ;multiple pixels
  image_dim = 256.
  image_1d = DBLARR(image_dim) * 0.
  image_1d[100:103] = 1.
  image_1d[106:110] = 1.
  image_1d[115:118] = 1.
  image_MPI = DBLARR(N_H_max+3, image_dim) * 0.
  FOR j=0, N_H_max-1 DO BEGIN
    temp_image_1d = DBLARR(n_H_100mT-1) * 0.
    FOR i=0, image_dim-1 DO BEGIN
      shift_dist = (i - image_dim/2.) / image_dim * FOV
      unit_dist = (H_array_H_max_100mT[j, 1]-H_array_H_max_100mT[j, 0])/gradient_field_array[j]
      shift_dist_ele = shift_dist / unit_dist
      IF i EQ 0 THEN shift_dist_ele_edge = shift_dist_ele
      temp_image_1d += image_1d[i] * SHIFT(REFORM(PSF_array_100mT[j,*]), shift_dist_ele)
    ENDFOR
    image_MPI[j,*] = CONGRID(temp_image_1d[(n_H_100mT-1)/2-shift_dist_ele:(n_H_100mT-1)/2+shift_dist_ele], image_dim)
  ENDFOR
  
  imageMPI_1 = image_MPI[N_H_max-1,*]/MAX(image_MPI[N_H_max-1,*])-MIN(image_MPI[N_H_max-1,*]/MAX(image_MPI[N_H_max-1,*]))
  imageMPI_0 = image_MPI[0,*]/MAX(image_MPI[0,*])-MIN(image_MPI[0,*]/MAX(image_MPI[0,*]))
  image_MPI[N_H_max, *] = ABS(imageMPI_1 - imageMPI_0)
  iplot, imageMPI_1
  iplot, imageMPI_0, /overplot
  iplot, REFORM(image_MPI[N_H_max, *]), color='blue', /overplot
  interv = 5
  FOR i=interv, image_dim-1-interv DO BEGIN
    neg_tre = image_MPI[N_H_max, i]-image_MPI[N_H_max, i-interv]
    pos_tre = image_MPI[N_H_max, i]-image_MPI[N_H_max, i+interv]
    IF neg_tre LT 0 AND pos_tre LT 0 THEN image_MPI[N_H_max+1, i] = exp(-REFORM(image_MPI[N_H_max, i])/tao)
  ENDFOR
  

  iplot, REFORM(image_MPI[N_H_max+1, *]), color='red', /overplot
  image_MPI[N_H_max+1, *] *= MAX(image_MPI)
  image_MPI[N_H_max+2, *] = image_1d*MAX(image_MPI)
  iimage, image_MPI
  window, 2
  plot, image_1d/MAX(image_1d)*MAX(image_MPI), YRANGE=[MIN(image_MPI), MAX(image_MPI)], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Index', ytitle='Intensity', TITLE='Pixel Intensity '+STRTRIM(STRING(particle_size),1)+' nm'
  FOR j=0, N_H_max-1 DO BEGIN
    oplot, REFORM(image_MPI[j,*]), LINESTYLE=2, color=3500
  ENDFOR

  ;calculate relative change and slope and fitting
  image_MPI_norm = image_MPI * 0.
  slope_array = image_1d * 0.
  chisq_array = image_1d * 0.
  FOR i=0, image_dim-1 DO BEGIN
    image_MPI_norm[*,i] = image_MPI[*,i]/image_MPI[0,i]
    result = LINFIT(H_max_array, REFORM(image_MPI_norm[0:N_H_max-1,i]), chisq=chisq_temp)
    slope_array[i] = result[1]
    chisq_array[i] = chisq_temp
  ENDFOR
  ;plot the slope and linearity
  iplot, slope_array, title = 'slope array'
  iplot, chisq_array, title = 'chisq array'

  window, 3
  plot, H_max_array, REFORM(image_MPI_norm[0:N_H_max-1,0]), YRANGE=[MIN(image_MPI_norm), MAX(image_MPI_norm)], LINESTYLE=2, BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='H max [mT]', ytitle='Signal Ratio', TITLE='Relative Signal '+STRTRIM(STRING(particle_size),1)+' nm'
  FOR i=1, image_dim-1 DO BEGIN
    oplot, H_max_array, REFORM(image_MPI_norm[0:N_H_max-1,i]), LINESTYLE=2, color=3500
    IF image_1d[i] NE 0 THEN oplot, H_max_array, REFORM(image_MPI_norm[0:N_H_max-1,i]), LINESTYLE=0, color=6500, THICK=2
  ENDFOR

END

;@GJ, 2023/5/4, Gradient sweeping MPI by analyzing two points
PRO gradient_sweeping_MPI_two_point

  N_H_max = 20.
  H_max_array = FINDGEN(N_H_max)+20.
  H_max = H_max_array[0]
  Langevin_derivative_PSF_estimation, particle_size, H_max, H_array, PSF0, H_array_100mT, PSF_100mT
  n_H = N_ELEMENTS(H_array)
  PSF_array = FINDGEN(N_H_max, n_H-1)
  window, 0
  FOR i=0, N_H_max-1 DO BEGIN
    Langevin_derivative_PSF_estimation, particle_size, H_max_array[i], H_array, PSF, H_array_100mT, PSF_100mT
    PSF_array[i,*] = PSF
  ENDFOR

  ;@GJ, 2023/5/4, averaging two points
  PSF_array_temp = PSF_array
  PSF_array = (SHIFT(PSF_array_temp, 0, n_H/20) + SHIFT(PSF_array_temp, 0, -n_H/20))/2.
  window, 0
  plot, H_array[0:n_H-1], REFORM(PSF_array[N_H_max-1,*]), XRANGE=[-H_max/3,H_max/3], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Location [mm]', ytitle='PSF', TITLE='PSF '+STRTRIM(STRING(particle_size),1)+' nm'
  FOR i=1, N_H_max-1 DO oplot, H_array[0:n_H-1], REFORM(PSF_array[N_H_max-1-i,*]), LINESTYLE=2, color=3500
  ;  oplot, H_array[0:n_H-1], PSF3, LINESTYLE=3, color=6000

  PSF_array_norm = PSF_array*0.
  FOR j=0, N_H-2 DO PSF_array_norm[*,j] = PSF_array[*,j]/PSF_array[0,j]

  window, 1
  plot, H_max_array, REFORM(PSF_array_norm[*,0]), LINESTYLE=3, YRANGE=[MIN(PSF_array_norm), MAX(PSF_array_norm)], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='H max [mT]', ytitle='PSF', TITLE='Relative Signal '+STRTRIM(STRING(particle_size),1)+' nm'
  FOR j=0, n_H/2.-1,10 DO BEGIN
    oplot, H_max_array, REFORM(PSF_array_norm[*,j]), LINESTYLE=3, color=3500
    ;    wait, 1
  ENDFOR
  oplot, H_max_array, REFORM(PSF_array_norm[*,n_H/2.-1-n_H/20]), color=0.

END


PRO read_MRI_Reconstruction_txt, filename, data_array

  ;read filename
  IF N_ELEMENTS(filename) EQ 0 THEN BEGIN
    filename =  DIALOG_PICKFILE(/READ, PATH='C:\D_drive\MPI_Tianjie\Shen_MPI_simulator\', FILTER = '*.txt')
  ENDIF ELSE BEGIN
    WHILE STRPOS(filename, '.txt') EQ -1 DO BEGIN
      filename =  DIALOG_PICKFILE(/READ, PATH='C:\D_drive\MPI_Tianjie\Shen_MPI_simulator\', FILTER = '*.txt')
    ENDWHILE
  ENDELSE
  
  data_dim = 121-2
  OPENR, lun, filename, /get_lun
  data_array = DBLARR(data_dim, data_dim) * 0.
  
  FOR i=0L, data_dim-1 DO BEGIN
    data_string=''
    READF, lun, data_string, FORMAT='(%"%s")'
    temp_data=STRSPLIT(data_string, ' ', /extract)
    data_array[*, data_dim-1-i] = DOUBLE(temp_data)
  ENDFOR
  
  ; Close the file and free the file unit
  FREE_LUN, lun
;  iimage, data_array, title='reconstucted'

END

PRO analyze_MRI_Rec_data

  dir_name = 'C:\D_drive\MPI_Tianjie\Shen_MPI_simulator\'
;  dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)
;  IF N_ELEMENTS(dir_name) EQ 0 THEN RETURN
  ;count the number of files
  a_file = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  read_MRI_Reconstruction_txt, a_file[0], data_array0
  mask_image = IMAGE_THRESHOLD(data_array0, THRESHOLD=tf, /MEAN)
  data_array_gs = DBLARR(nrfile, (SIZE(data_array0))[1], (SIZE(data_array0))[2]) * 0.
  gradient_array = DBLARR(nrfile) * 0.
  FOR i=0L, nrfile-1L DO BEGIN
    read_MRI_Reconstruction_txt, a_file[i], data_array_temp
    data_array_gs[i, *, *] = data_array_temp
    gradient_array[i] = DOUBLE(FILE_BASENAME(a_file[i], '.txt'))/10.
  ENDFOR
  
  ivolume, data_array_gs

  image_ori_fn = 'C:\D_drive\MPI_Tianjie\Shen_MPI_simulator\MPIRF\dist\example.png'
  READ_PNG, image_ori_fn, image_ori 
  image_2d = CONGRID(REFORM(image_ori[0,*,*]),(SIZE(data_array0))[1], (SIZE(data_array0))[2])
  iimage, image_2d, title='ori'
  iimage, data_array_temp, title='MPI '+STRTRIM(STRING(gradient_array[nrfile-1]),1)+' T/m'
  
  ind=0
  slope_array = image_2d * 0.-2.
  chisq_array = image_2d * 0.
  FOR i=0, (SIZE(data_array0))[1]-1 DO BEGIN
    FOR j=0, (SIZE(data_array0))[2]-1 DO BEGIN
      IF mask_image[i,j] GT 0 THEN BEGIN
        signal_array = REFORM(data_array_gs[*, i, j])/data_array_gs[0, i, j]
        result = LINFIT(gradient_array, signal_array, chisq=chisq_temp)
        slope_array[i,j] = result[1]
        chisq_array[i,j] = chisq_temp
        IF ind EQ 0 THEN BEGIN
          window, 2
          plot, gradient_array, signal_array, YRANGE=[0.5,1.5], BACKGROUND = 'FFFFFF'x, COLOR = 0, xtitle='Gradient [T/m]', ytitle='Signal Ratio', TITLE='Pixel Intensity'
        ENDIF ELSE BEGIN
          IF image_2d[i,j] GT 128 THEN BEGIN
            color_temp=0&linstyle_temp=0
          ENDIF ELSE BEGIN
            color_temp=3500&linstyle_temp=2
          ENDELSE
          oplot, gradient_array, signal_array, LINESTYLE=linstyle_temp, color=color_temp
        ENDELSE
        ind++
;        wait, 0.1
      ENDIF
    ENDFOR
  ENDFOR
  
  iimage, slope_array, title='slope'
  iimage, chisq_array, title='chisq'
  
END

;@GJ, 2023/5/13, fusion of MPI and CT of liver
PRO MPI_CT_fusion

  fn_MRI = 'C:\Program Files\Exelis\IDL83\examples\demo\demodata\mri.sav'
  restore, fn_MRI
  fn_PET = 'C:\Program Files\Exelis\IDL83\examples\demo\demodata\pet.sav'
  restore, fn_PET
  
  
END



;GJ, 2023/5/26, Nanorod magnetic field calculation
PRO Nanorod_magnetic_field

  field_amp = DBLARR(1001,1001)*0. ; 100 rows, 201 columns
  x_field = DBLARR(1001,1001)*0. ; 100 rows, 201 columns
  y_field = DBLARR(1001,1001)*0. ; 100 rows, 201 columns
  d = 10.; nm
  coeff=100.
  FOR i=0, 1000 DO BEGIN
    FOR j=0, 1000 DO BEGIN
      dist = SQRT((i - 500.)^2 + (j - 500.)^2)*0.01 ;nm
      field_amp[i,j] = (d/dist)^3;*coeff
      sin_theta = (i - 500.)/SQRT((i - 500.)^2 + (j - 500.)^2)
      cos_theta = (j - 500.)/SQRT((i - 500.)^2 + (j - 500.)^2)
      x_field[i,j] = field_amp[i,j] * 3 * sin_theta * cos_theta
      y_field[i,j] = field_amp[i,j] * (3 * cos_theta^2 - 1)
    ENDFOR
  ENDFOR

  x_field = ROT(x_field, 90, 1.0, /INTERP)
  y_field = ROT(y_field, 90, 1.0, /INTERP)
  ;  iimage, field_amp, title = 'Field Amplitude of sphere'
  ;  iimage, x_field, title = 'x-direction Field of sphere'
  ;  iimage, y_field, title = 'y-direction Field of sphere'

  sum_x_field = x_field
  sum_y_field = y_field
  sum_field_amp = DBLARR(1001,1001)*0.

  FOR i=1, 100 DO BEGIN
    sum_x_field += SHIFT(x_field, i, 0) + SHIFT(x_field, -i, 0)
    sum_y_field += SHIFT(y_field, i, 0) + SHIFT(y_field, -i, 0)
  ENDFOR

  sum_field_amp = SQRT(sum_x_field^2 + sum_y_field^2)

  iimage, sum_field_amp, title='Field Amplitude of Nanorod'
  ;  isurface, sum_field_amp, title='Field Amplitude of Anderson coils'
  contour, sum_field_amp, levels=FINDGEN(100)/99.*max(sum_field_amp), title='Field Amplitude'
  iimage, sum_x_field, title = 'x-direction Field of Nanorod'
  iimage, sum_y_field, title = 'y-direction Field of Nanorod'

  ;  iplot, sum_field_amp[*, 50], title='Field Amplitude on x-axis'
  ;  iplot, sum_field_amp[50, *], title='Field Amplitude on y-axis'
END

;@GJ, 2023/5/29, analyze the twin images based on gradient matching
PRO twin_image_analyze
  ;loading DICOM image and get informations
;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\dandu teshutanzhen.dcm'
;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'

  image_dcm = read_dicom(fn_dcm)
  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = 256
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = 256

  OBJ_DESTROY, obj
  
  ;calculate the gradient map
  gradient_map = DBLARR(4, cols, rows) * 0.
  window_width = 10
  x_row = FINDGEN(window_width+1)*pixelSp[0]
  y_col = FINDGEN(window_width+1)*pixelSp[1]
  ;select points
  mask_points = image_dcm * 0.
  FOR i= window_width, cols-window_width-1 DO BEGIN
    FOR j= window_width, rows-window_width-1 DO BEGIN
      ;flag the pixel
      flag = 0
      ;calcualte gradient
      result = LINFIT(x_row, REVERSE(REFORM(image_dcm[i-window_width:i,j])), chisq=chisq_temp)
      gradient_map[0, i, j] = result[1]
      IF gradient_map[0, i, j] GT 0 THEN flag = 1
      result = LINFIT(x_row, REFORM(image_dcm[i:i+window_width,j]), chisq=chisq_temp)
      gradient_map[1, i, j] = result[1]
      IF gradient_map[1, i, j] GT 0 THEN flag = 1
      result = LINFIT(y_col, REVERSE(REFORM(image_dcm[i, j-window_width:j])), chisq=chisq_temp)
      gradient_map[2, i, j] = result[1]
      IF gradient_map[2, i, j] GT 0 THEN flag = 1
      result = LINFIT(y_col, REFORM(image_dcm[i, j:j+window_width]), chisq=chisq_temp)
      gradient_map[3, i, j] = result[1]
      IF gradient_map[3, i, j] GT 0 THEN flag = 1
      
      IF flag EQ 0 THEN mask_points[i, j] = 1
    ENDFOR
  ENDFOR

;  iimage, REVERSE(image_dcm,2), title='original image'  
;  iimage, REVERSE(REFORM(gradient_map[0,*,*]),2), title='x- gradient map'
;  iimage, REVERSE(REFORM(gradient_map[1,*,*]),2), title='x+ gradient map'
;  iimage, REVERSE(REFORM(gradient_map[2,*,*]),2), title='y- gradient map'
;  iimage, REVERSE(REFORM(gradient_map[3,*,*]),2), title='y+ gradient map'
  
  iimage, mask_points, title='points for pairs'
  
  ;Estimate twin's distance based on particle's coercivity
  coercivity = 32.2; mT
  gradient = 5.7; mT/mm
  dist_twin = 2.*coercivity/gradient; in mm
  dist_twin_pixels = dist_twin/pixelSp[0]
  
  ;Looking for pairs
  true_pixel_map = image_dcm * 0.
  ind = WHERE(mask_points GT 0, count)
    
  IF count GT 2 THEN BEGIN
    ind_2d = ARRAY_INDICES(mask_points, ind)
    ; Compute the Euclidean distance between each point.
    DISTANCE = DISTANCE_MEASURE(ind_2d, /MATRIX)
    dis_ind = WHERE(DISTANCE GT 45 AND DISTANCE LT 46, dis_count)
    dis_ind_2d = ARRAY_INDICES(DISTANCE, dis_ind)
    FOR i = 0, dis_count-1 DO BEGIN
      true_pixel_map[ind[dis_ind_2d[0,i]]] = 10
    ENDFOR
    iimage, DISTANCE, title='Distance between points'

    iimage, mask_points+true_pixel_map, title='twin points'
    
    ;Find the paired pixel
    pair_pixel_map = image_dcm * 0.
    index_pixel = 1
    FOR i = 0, count-1 DO BEGIN
      flag = 0
      FOR j = 0, count-1 DO BEGIN
        IF DISTANCE[i , j] GT 45 AND DISTANCE[i,j] LT 46 AND ind_2d[0,j] GT ind_2d[0,i] AND pair_pixel_map[ind[j]] EQ 0 AND pair_pixel_map[ind[i]] EQ 0 THEN BEGIN
          IF flag EQ 0 THEN BEGIN
            pair_a_x = ind_2d[0, i]
            pair_a_y = ind_2d[1, i]
            pair_b_x = ind_2d[0, j]
            pair_b_y = ind_2d[1, j]
            pair_a = ind_2d[*, i]
          ENDIF ELSE BEGIN
            pair_b_x = [pair_b_x, ind_2d[0, j]]
            pair_b_y = [pair_b_y, ind_2d[1, j]]
          ENDELSE
          ;increasing flag by 1
          flag++
        ENDIF
      ENDFOR
      IF flag GT 1 THEN BEGIN
        pair_b = [MEDIAN(pair_b_x), MEDIAN(pair_b_y)]
        pair_pixel_map[pair_a_x, pair_a_y] = index_pixel
        pair_pixel_map[pair_b_x, pair_b_y] = index_pixel
        pair_pixel_map[(pair_b_x+pair_a_x)/2, (pair_b_y+pair_a_y)/2] = index_pixel
        index_pixel++
      ENDIF
    ENDFOR
    iimage, pair_pixel_map, title='pair points'
    
  ENDIF
 
END

;@GJ, 2023/6/15, do the calculation simulation
PRO TWIN_scan_simulation

  filename2='C:\D_drive\MPI_Tianjie\MIP\000.tif'
  image_temp2 = read_image(filename2)
  image = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 64, 64), MAX=130)
  ;      image = image * 0.
  ;      image[100:102, 127:129] = 256
  ;      image[130:132, 127:129] = 256

  ;setup the 1d image encoding and decoding kernal
;  iimage,  image, title='Original'
  Matrix_size = N_ELEMENTS(image[0,*]); causing invert calcuation error
  H_2d_array = DBLARR(Matrix_size, Matrix_size)
  x_size = N_ELEMENTS(image[*,0])
  ;set up the encoding magnetic fields
  H0 = 0.000;5;15 ;T
  G = 0.075; T/m
  delta_G = G*2./Matrix_size
  FOV = 0.2; m
  x_res = FOV/x_size
  G_array = (FINDGEN(Matrix_size)-Matrix_size/2)*delta_G
  ;G_array = DBLARR(Matrix_size)*0. - G
  Hx = ((FINDGEN(x_size*2)-x_size)*G*x_res + H0) * 1000. ;mT
  x_axis = (FINDGEN(x_size*2)-x_size)*x_res * 1000. ;mm
;  iplot, x_axis, Hx
  FOR i=0, Matrix_size-1 DO H_2d_array[i, *] = (((FINDGEN(x_size)-x_size)+i)*G/3.*x_res+H0) * 1000. ; mT


  H_flat_portion = 0.0
  H_frequency = 2.0; kHz
  H_info = [13.0, H_flat_portion, H_frequency]
  tracer_info = [0, 30., 50.] ;type, size, relaxation time
  s_f = 1
  signalKernal_2d_array = (signal_vs_H(H_info,tracer_info, s_f, s_array))[1, 1] ;3rd harmonic without relaxation
  
  print, signalKernal_2d_array
  iplot, REFORM(s_array[0,*]), REFORM(s_array[1,*]), title='original signal'
  iplot, REFORM(s_array[0,*]), REFORM(s_array[2,*]), title='filtered signal'
  iplot, REFORM(s_array[0,*]), REFORM(s_array[3,*]), title='H curve'


END

;@GJ, 2023/6/16, simulating the 1D image of twin imaging of magnetic nanoparticles
PRO TWIN_1D_simulation_MNP
  ;define the image
  Matrix_size = 256
  x_size = Matrix_size
  image_1d = DBLARR(Matrix_size) * 0.
  image_1d[Matrix_size/2-1]  = 256
  
  ;define the gradient field
  G = -5.7; T/m
  FOV = 0.06; m
  x_res = FOV/Matrix_size
  x_array = (FINDGEN(Matrix_size)-Matrix_size/2)*x_res
  H_G = G * x_array * 1000. ; mT
  
  ;define the AC excitation field
  H_offset = 0.0;5;15 ;mT
  H_max = 20.0; mT
  frequency = 20. ;kHz
  ;define the period based on frequency
  period = 1.0/frequency; ms
  Nsample = 1000.
  delta_t = period/Nsample * 1000; us
  ;GJ, 2021/11/25, define the time series
  t = FINDGEN(Nsample) * delta_t ; us
  ;this is a cosine curve
  H_t = COS(2.*!PI*frequency*t/1000.)*H_max
  
  ;define the basic Langevin function parameters
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  particle_size = 30.;10.;50.;30.;20.;30. ; nm
  ;Msat = 446.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
  
  ;calculate the magnetization
  M_array = DBLARR(Nsample)*0.
  FOR i=0, Nsample-1 DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_temp = H_t[i] + H_G[j]
      M_temp = Msat_kAm*(1./TANH(beta_particle*H_temp) - 1./(beta_particle*H_temp))
      M_array[i] = M_array[i] + M_temp * image_1d[j]
    ENDFOR
  ENDFOR
  
  ;deal with finite
  index_nan = where(~FINITE(M_array), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(M_array, t, t[index_nan], /NAN)
    MM_array[index_nan] = RESULT
  ENDIF
  
  ;calculate the signal
  signal_array = (M_array[1:*]-M_array[0:Nsample-1])/delta_t; kA/m/us
  relaxation_time = 3. ;us
  IF relaxation_time GT 0 THEN signal_array = non_adiabatic(t, signal_array, relaxation_time)
  iplot, t[0:Nsample-2], signal_array, xtitle='time [us]', ytitle='signal'
  
  ;X-space reconstruction
  ;calculate the position of FFP
  xs_t_array = H_t/ABS(G)
  xs_t_velocity =  (xs_t_array[1:*]-xs_t_array[0:Nsample-1])/delta_t
;  iplot, t[0:Nsample-1], xs_t_array, xtitle='time [us]', ytitle='FFR location'
;  iplot, t[0:Nsample-2], xs_t_velocity, xtitle='time [us]', ytitle='FFR velocity'
  image_native = signal_array/ABS(G)/xs_t_velocity
  index_zero = WHERE(ABS(xs_t_velocity) LT 0.1, count)
  IF count GT 0 THEN image_native[index_zero] = 0.
  image_rec_a = CONGRID(REFORM(image_native[Matrix_size/10.:Nsample/2-1]), Matrix_size)
  image_rec_b = REVERSE(CONGRID(REFORM(image_native[Nsample/2:Nsample-1-Matrix_size/10.]), Matrix_size))
  image_rec_ab = (image_rec_a+image_rec_b)/2.
;  iplot, image_rec_a, title='native image 1d a'
;  iplot, image_rec_b, title='native image 1d b'
  
  image_combine = DBLARR(Matrix_size, Matrix_size)
  FOR i=0, Matrix_size/4.-1 DO image_combine[i, *] = image_1d
  FOR i=Matrix_size/4., Matrix_size/4.*2.-1 DO image_combine[i, *] = BYTSCL(image_rec_a)
  FOR i=Matrix_size/4.*2., Matrix_size/4.*3.-1 DO image_combine[i, *] = BYTSCL(image_rec_b)
  FOR i=Matrix_size/4.*3., Matrix_size-1 DO image_combine[i, *] = BYTSCL(image_rec_ab)
  iimage, TRANSPOSE(image_combine), title='combined image'
END

;calculate langevin function
PRO Langevin_pulse_nanorubble, beta

  dir_name = 'F:\MPI_Tianjie\Nanorubble_MPI\20230602Relax\02_PulsedExc\'
;  dir_name = DIALOG_PICKFILE(PATH=dir_name, /MUST_EXIST, TITLE="txt files", /DIRECTORY)

  a = FILE_SEARCH(dir_name, '*.txt', count=nrfile)
  M_array = DBLARR(nrfile*2) * 0.
  H_array = DBLARR(nrfile*2) * 0.
  items = STRARR(nrfile);
  psyms = [-15, -16, -17, -18, -19, 20, -21, -22, -23]
  colors = ['pink', 'green', 'red3', 'red4', 'blue', 'red5', 'red6', 'red7', 'red8']

  FOR i=nrfile-1L, 0L, -1 DO BEGIN
    filename = FILE_BASENAME(a[i], '.txt')
    fit_para = read_pulsed_signal_dat(a[i])
    ;    IF H_array[i*2+1] LT 6.9 THEN coe = 1./50. ELSE coe = 1./50.;/20.
    coe = 1.
  
    nele = N_ELEMENTS(fit_para.H)
    H_array[i*2] = fit_para.H[nele/2 -1 ]
    H_array[i*2+1] = fit_para.H[nele - 1]
  
    M_array[i*2] = fit_para.M[nele/2 -1 ] * coe
    M_array[i*2+1] = fit_para.M[nele - 1] * coe
    
    IF i EQ nrfile-1L THEN BEGIN
      xtitle = 'H [mT/u0]'
      ytitle = 'M [A.U.]'
      filename_arr = STRSPLIT(filename, '_', /EXTRACT)
      IF N_ELEMENTS(filename_arr) GT 1 THEN BEGIN
        title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+' ' + filename_arr[1] +')';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
      ENDIF ELSE BEGIN
        title = 'M-H curves ('+STRTRIM(STRING(fit_para.frequency), 1)+' Hz, '+filename_arr[0]+')';, ' + STRTRIM(STRING(H_max, format=''), 1) + ' mT)'
      ENDELSE
      window, 6
      cgPlot, fit_para.H, fit_para.M * coe, Color=colors[i MOD N_ELEMENTS(colors)], $; BACKGROUND='black', $; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], $
        Title=title, XTitle=xtitle, YTitle=ytitle, Thick=2, XRANGE=[-MAX(H_array)*1.1, MAX(H_array)*1.1], YRANGE=[-MAX(M_array)*1.1, MAX(M_array)*1.1]
      ;GJ, 2022/2/9, plot the fitted M curve
      ;      cgPlots, fit_para.debye_H, fit_para.debye_M, Color=colors[i], Thick=2
      ;      cgPlots, fit_para.relax_H, fit_para.relax_mono_Mfit, Color='black', Thick=2, Linestyle=1
;      cgPlots, fit_para.relax_H, fit_para.relax_biex_Mfit * coe, Color='black', Thick=2, Linestyle=1
;      cgPlots, fit_para.debye_H, fit_para.debye_mono_Mfit * coe, Color='black', Thick=2, Linestyle=2
    ENDIF ELSE BEGIN
      cgPlots, fit_para.H, fit_para.M * coe, Color=colors[i MOD N_ELEMENTS(colors)], Thick=2 ; PSym=psyms[i],  SymSize=0.5, SymColor=colors[i], Thick=1
      ;GJ, 2022/2/9, plot the fitted M curve
      ;      cgPlots, fit_para.debye_H, fit_para.debye_M, Color=colors[i], Thick=2
      ;      cgPlots, fit_para.relax_H, fit_para.relax_mono_Mfit, Color='black', Thick=2, Linestyle=1
;      cgPlots, fit_para.relax_H, fit_para.relax_biex_Mfit * coe, Color='black', Thick=2, Linestyle=1
;      cgPlots, fit_para.debye_H, fit_para.debye_mono_Mfit * coe, Color='black', Thick=2, Linestyle=2
    ENDELSE

  ENDFOR
  
  sort_ind = SORT(H_array)
  H_array = H_array[sort_ind]
  M_array = M_array[sort_ind]
  
  
  cgplots, H_array, M_array, Color='black', Thick=2, Linestyle=2
END

;@GJ, 2023/6/28, doing a simple simulation based on 1d twin imaging
PRO TWIN_2D_simulation_nanorubble_single_point

  TWIN_1D_simulation_nanorubble_single_point, G, FOV, x_res, Matrix_size, image_rec_size, image_native, image_rec_ab, image_rec_ab_final
  
  help, image_native
  iplot, image_native
  
  x_size = Matrix_size
  x_array = (FINDGEN(x_size)-x_size/2)*x_res
  H_G = G * x_array * 1000. ; mT

  ;define ghe y gradient field
  y_size = Matrix_size;/8.;4;Matrix_size;4;Matrix_size/8;
  y_res = x_res
  FOV_y = y_res * y_size
  y_array = (FINDGEN(y_size)-y_size/2)*y_res
  H_y = G * y_array * 1000. ; mT
  
  temp_native = CONGRID(image_rec_ab_final, Matrix_size)
  iplot, temp_native, title='1d'
  image_rec_2d = DBLARR(Matrix_size, Matrix_size)
  Hx_comp_array = DBLARR(Matrix_size, Matrix_size)
;  Hx_comp_array = DBLARR(Matrix_size, Matrix_size)
  FOR k=0, y_size-1 DO BEGIN
    FOR j=0, x_size-1 DO BEGIN
      Hx = H_G[j]
      Hy = H_y[k]
      H_temp = SQRT(Hx^2 + Hy^2)
;      Hx_comp_array[j, k] = Hx/H_temp
      Hx_comp_array[j, k] = 1.-ABS(Hy/H_temp); EXP(-ABS(Hy/H_temp))
      image_rec_2d[j, k] = temp_native[j] * Hx_comp_array[j, k]
    ENDFOR
    IF Hy EQ 0 THEN iplot, image_rec_2d[*, k], title='1d in 2d rec'
  ENDFOR
  
  iimage, Hx_comp_array, title='x component'
  iimage, image_rec_2d, title='reconstructed 2d'
;    ;  FOR k=y_size/2, y_size/2 DO BEGIN
;    ;    k = y_size/2
;    ;calculate the magnetization
;    M_array = DBLARR(Nsample+1)*0.
;    Mx_array = DBLARR(Nsample+1)*0.
;    My_array = DBLARR(Nsample+1)*0.
;    H_temp_array = DBLARR(Nsample+1, x_size)*0.
;    Hx_comp = DBLARR(Nsample+1, x_size)*0.
;    Hy_comp = DBLARR(Nsample+1, x_size)*0.
;    FOR i=0, Nsample DO BEGIN
;      FOR j=0, x_size-1 DO BEGIN
;        ;calculate the field amplitude
;        H_x = H_t[i] + H_G[j]
;        H_temp = SIGNUM(H_x) * SQRT(H_x^2 + H_y[k]^2)
;        H_temp_array[i, j] = H_temp
;        Hx_comp[i, j] = H_x/H_temp
;        Hy_comp[i, j] = H_y[k]/H_temp
;
;        ;only calculate the non-zero pixels
;        IF image_1d[j] EQ 0 THEN CONTINUE ;j=127 is nonzero
;
;        ;calculate the magnetization
;        IF H_t[i+1]-H_t[i] GT 0 THEN BEGIN
;          M_temp = INTERPOL(M_curve_inc, H_curve_inc, H_temp)
;          ;          M_temp = INTERPOL(M_curve_inc, H_curve_inc, H_x)
;        ENDIF ELSE BEGIN
;          M_temp = INTERPOL(M_curve_dec, H_curve_dec, H_temp)
;          ;          M_temp = INTERPOL(M_curve_dec, H_curve_dec, H_x)
;        ENDELSE
;        M_array[i] = M_array[i] + M_temp * image_1d[j]
;      ENDFOR
;    ENDFOR
  

END


;@GJ, 2023/6/16, simulating the 1D image of twin imaging of magnetic nanorubbles
PRO TWIN_1D_simulation_nanorubble_single_point, G, FOV, x_res, Matrix_size, image_rec_size, image_native, image_rec_ab, image_rec_ab_final
  ;define the image
  Matrix_size = 256
  x_size = Matrix_size
  image_1d = DBLARR(Matrix_size) * 0.
  image_1d[Matrix_size/2-1]  = 256

  ;define the gradient field
  G = -5.7; T/m
  FOV = 0.06; m
  x_res = FOV/Matrix_size
  x_array = (FINDGEN(Matrix_size)-Matrix_size/2)*x_res
  H_G = G * x_array * 1000. ; mT

  ;define the AC excitation field
  H_offset = 0.0;5;15 ;mT
  H_max = 130.; mT
  frequency = 20. ;kHz
  ;define the period based on frequency
  period = 1.0/frequency; ms
  Nsample = 1000.
  delta_t = period/Nsample * 1000; us
  ;GJ, 2021/11/25, define the time series
  t = FINDGEN(Nsample+2) * delta_t ; us
  ;this is a cosine curve
  H_t = COS(2.*!PI*frequency*t/1000.)*H_max
  
;  fn = 'F:\MPI_Tianjie\NanoRod_MPI\20030602Relax\01_SineExc\2023060503.txt'
;  M_H_array = read_sine_signal_dat(fn)
;  H_curve = REFORM(M_H_array[0,*])
;  Nele_MH = N_ELEMENTS(H_curve)
;  M_curve = -REFORM(M_H_array[1,*])
;  MIN_H = MIN(H_curve, min_ind)
;  H_curve = SHIFT(H_curve, -min_ind)
;  M_curve = SHIFT(M_curve, -min_ind)
;;  M_H_array = SHIFT(M_H_array, 0, -min_ind)
;  H_curve_inc = H_curve[0:Nele_MH/2-1]
;  M_curve_inc = M_curve[0:Nele_MH/2-1]
;  H_curve_dec = H_curve[Nele_MH/2:*]
;  M_curve_dec = M_curve[Nele_MH/2:*]
  
  M_H_array=read_VSM_file(filename)
  H_curve = REFORM(M_H_array[0,*])
  Nele_MH = N_ELEMENTS(H_curve)
  M_curve = REFORM(M_H_array[1,*])
  H_curve_dec = H_curve[0:Nele_MH/2-1]
  M_curve_dec = M_curve[0:Nele_MH/2-1]
  H_curve_inc = H_curve[Nele_MH/2:*]
  M_curve_inc = M_curve[Nele_MH/2:*]
  
  ;calculate the magnetization
  M_array = DBLARR(Nsample+1)*0.
  FOR i=0, Nsample DO BEGIN
    FOR j=0, Matrix_size-1 DO BEGIN
      H_temp = H_t[i] + H_G[j]
      ;M_temp = Msat_kAm*(1./TANH(beta_particle*H_temp) - 1./(beta_particle*H_temp))
;      IF i LT Nsample AND ABS(H_temp) LT ABS(MIN_H) THEN BEGIN
        IF H_t[i+1]-H_t[i] GT 0 THEN BEGIN
          M_temp = INTERPOL(M_curve_inc, H_curve_inc, H_temp)
        ENDIF ELSE BEGIN
          M_temp = INTERPOL(M_curve_dec, H_curve_dec, H_temp)
        ENDELSE
 ;     ENDIF ELSE BEGIN
 ;       IF H_temp GT ABS(MIN_H) THEN M_temp = MAX(M_curve) ELSE M_temp = MIN(M_curve)
 ;     ENDELSE

      M_array[i] = M_array[i] + M_temp * image_1d[j]
    ENDFOR
  ENDFOR

  ;deal with finite
  index_nan = where(~FINITE(M_array), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(M_array, t, t[index_nan], /NAN)
    MM_array[index_nan] = RESULT
  ENDIF

  ;calculate the signal
  signal_array = (M_array[1:*]-M_array[0:Nsample])/delta_t; kA/m/us
  relaxation_time = 3. ;us
  IF relaxation_time GT 0 THEN signal_array = non_adiabatic(t, signal_array, relaxation_time)
  iplot, t[0:Nsample-1], signal_array, xtitle='time [us]', ytitle='signal'

  ;X-space reconstruction
  ;calculate the position of FFP
  xs_t_array = H_t/ABS(G)
  xs_t_velocity =  (xs_t_array[1:*]-xs_t_array[0:Nsample-1])/delta_t
  ;  iplot, t[0:Nsample], xs_t_array, xtitle='time [us]', ytitle='FFR location'
  ;  iplot, t[0:Nsample-1], xs_t_velocity, xtitle='time [us]', ytitle='FFR velocity'
  image_native = signal_array/ABS(G)/xs_t_velocity
  index_zero = WHERE(ABS(xs_t_velocity) LT 0.1, count)
  IF count GT 0 THEN image_native[index_zero] = 0.
  image_rec_a = CONGRID(image_native, Matrix_size*2)
  image_rec_b = REVERSE(image_rec_a)
  image_rec_ab = (image_rec_a+image_rec_b)/2.
  ;  iplot, image_rec_a, title='native image 1d a'
  ;  iplot, image_rec_b, title='native image 1d b'
  
  image_rec_size = ROUND(Matrix_size*0.8)
  image_combine = DBLARR(Matrix_size, image_rec_size)
  FOR i=0, Matrix_size/4.-1 DO image_combine[i, *] = image_1d[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1]/MAX(image_1d[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])*MAX(image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])
  FOR i=Matrix_size/4., Matrix_size/4.*2.-1 DO image_combine[i, *] = image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1]/MAX(image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])*MAX(image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])
  FOR i=Matrix_size/4.*2., Matrix_size/4.*3.-1 DO image_combine[i, *] = image_rec_b[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1]/MAX(image_rec_b[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])*MAX(image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])
  FOR i=Matrix_size/4.*3., Matrix_size-1 DO image_combine[i, *] = image_rec_ab[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1]/MAX(image_rec_ab[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])*MAX(image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])
  iimage, BYTSCL(TRANSPOSE(image_combine)), title='combined image'
  
  image_rec_ab_final = image_rec_ab[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1]/MAX(image_rec_ab[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])*MAX(image_rec_a[Matrix_size*0.1:Matrix_size*0.1+image_rec_size-1])
END


;@GJ, 2023/6/16, simulating the 2D image of twin imaging of magnetic nanorubbles
PRO TWIN_2D_simulation_nanorubble_single_point_old
  ;define the image
  Matrix_size = 256
  image_1d = DBLARR(Matrix_size) * 0.
  image_1d[Matrix_size/2-1]  = 256

  ;define the x gradient field
  G = -5.7; T/m
  FOV = 0.06; m
  x_size = Matrix_size
  x_res = FOV/Matrix_size
  x_array = (FINDGEN(x_size)-x_size/2)*x_res
  H_G = G * x_array * 1000. ; mT
  
  ;define ghe y gradient field
  y_size = Matrix_size;/8.;4;Matrix_size;4;Matrix_size/8;
  y_res = x_res
  FOV_y = y_res * y_size
  y_array = (FINDGEN(y_size)-y_size/2)*y_res
  H_y = G * y_array * 1000. ; mT
  
  ;define the AC excitation field
  H_coercivity = 13.;
  H_offset = 0.0;5;15 ;mT
  H_max = 130.; mT
  frequency = 20. ;kHz
  ;define the period based on frequency
  period = 1.0/frequency; ms
  Nsample = 1000.
  delta_t = period/Nsample * 1000; us
  ;GJ, 2021/11/25, define the time series
  t = FINDGEN(Nsample+2) * delta_t ; us
  ;this is a cosine curve
  H_t = COS(2.*!PI*frequency*t/1000.)*H_max

;  fn = 'F:\MPI_Tianjie\NanoRod_MPI\20030602Relax\01_SineExc\2023060503.txt'
;  M_H_array = read_sine_signal_dat(fn)
;  H_curve = REFORM(M_H_array[0,*])
;  Nele_MH = N_ELEMENTS(H_curve)
;  M_curve = -REFORM(M_H_array[1,*])
;  MIN_H = MIN(H_curve, min_ind)
;  H_curve = SHIFT(H_curve, -min_ind)
;  M_curve = SHIFT(M_curve, -min_ind)
;  ;  M_H_array = SHIFT(M_H_array, 0, -min_ind)
;  H_curve_inc = H_curve[0:Nele_MH/2-1]
;  M_curve_inc = M_curve[0:Nele_MH/2-1]
;  H_curve_dec = H_curve[Nele_MH/2:*]
;  M_curve_dec = M_curve[Nele_MH/2:*]

    M_H_array=read_VSM_file(filename)
    H_curve = REFORM(M_H_array[0,*])
    Nele_MH = N_ELEMENTS(H_curve)
    M_curve = REFORM(M_H_array[1,*])
    H_curve_dec = H_curve[0:Nele_MH/2-1]
    M_curve_dec = M_curve[0:Nele_MH/2-1]
    H_curve_inc = H_curve[Nele_MH/2:*]
    M_curve_inc = M_curve[Nele_MH/2:*]
  
  imagex_2d = DBLARR(Matrix_size, Matrix_size)*0.
  imagey_2d = DBLARR(Matrix_size, Matrix_size)*0.
  image_2d = DBLARR(Matrix_size, Matrix_size)*0.
  FOR k=0, y_size-1 DO BEGIN
;  FOR k=y_size/2, y_size/2 DO BEGIN
;    k = y_size/2
    ;calculate the magnetization
    M_array = DBLARR(Nsample+1)*0.
    Mx_array = DBLARR(Nsample+1)*0.
    My_array = DBLARR(Nsample+1)*0.
    H_temp_array = DBLARR(Nsample+1, x_size)*0.
    Hx_comp = DBLARR(Nsample+1, x_size)*0.
    Hy_comp = DBLARR(Nsample+1, x_size)*0.
    FOR i=0, Nsample DO BEGIN
      FOR j=0, x_size-1 DO BEGIN
        ;calculate the field amplitude
        H_x = H_t[i] + H_G[j]
        H_temp = SIGNUM(H_x) * SQRT(H_x^2 + H_y[k]^2)
        H_temp_array[i, j] = H_temp
        Hx_comp[i, j] = H_x/H_temp
        Hy_comp[i, j] = H_y[k]/H_temp

        ;only calculate the non-zero pixels
        IF image_1d[j] EQ 0 THEN CONTINUE ;j=127 is nonzero
        
        ;calculate the magnetization
        IF H_t[i+1]-H_t[i] GT 0 THEN BEGIN
          M_temp = INTERPOL(M_curve_inc, H_curve_inc, H_temp)
;          M_temp = INTERPOL(M_curve_inc, H_curve_inc, H_x)
        ENDIF ELSE BEGIN
          M_temp = INTERPOL(M_curve_dec, H_curve_dec, H_temp)
;          M_temp = INTERPOL(M_curve_dec, H_curve_dec, H_x)
        ENDELSE
        M_array[i] = M_array[i] + M_temp * image_1d[j]
      ENDFOR
    ENDFOR

    relaxation_time = 3. ;us
    IF relaxation_time GT 0 THEN BEGIN
      M_array = non_adiabatic(t, M_array, relaxation_time)
    ENDIF

    ;calculate the signal
    signal_array = (M_array[1:*]-M_array[0:Nsample])/delta_t; kA/m/us
    
    signalx_array = signal_array * REFORM(Hx_comp[0:Nsample-1, 0])
    signaly_array = signal_array * REFORM(Hy_comp[0:Nsample-1, 0])
    
    ;X-space reconstruction
    ;calculate the position of FFP
    xs_t_array = H_t/ABS(G)
    xs_t_velocity =  ABS((xs_t_array[1:*]-xs_t_array[0:Nsample-1]))/delta_t
    ;  iplot, t[0:Nsample], xs_t_array, xtitle='time [us]', ytitle='FFR location'
    ;  iplot, t[0:Nsample-1], xs_t_velocity, xtitle='time [us]', ytitle='FFR velocity'
    imagex_native = signalx_array/ABS(G)/xs_t_velocity
    index_zero = WHERE(ABS(xs_t_velocity) LT 0.1, count)
    IF count GT 0 THEN imagex_native[index_zero] = 0.
    imagex_rec_a = CONGRID(imagex_native, Matrix_size*2)
    imagex_rec_b = REVERSE(imagex_rec_a)
    imagex_rec_ab = (imagex_rec_a+imagex_rec_b)/2.
    ;  iplot, image_rec_a, title='native image 1d a'
    ;  iplot, image_rec_b, title='native image 1d b'
    imagex_2d[k, *] = imagex_rec_ab[0:Matrix_size-1]
    
    imagey_native = signaly_array/ABS(G)/xs_t_velocity
    index_zero = WHERE(ABS(xs_t_velocity) LT 0.1, count)
    IF count GT 0 THEN imagey_native[index_zero] = 0.
    imagey_rec_a = CONGRID(imagey_native, Matrix_size*2)
    imagey_rec_b = REVERSE(imagey_rec_a)
    imagey_rec_ab = (imagey_rec_a+imagey_rec_b)/2.
    imagey_2d[k, *] = imagey_rec_ab[0:Matrix_size-1]
    
    image_rec_a = SQRT(imagex_rec_a^2 + imagey_rec_a^2)
    image_rec_b = SQRT(imagex_rec_b^2 + imagey_rec_b^2)
    image_rec_ab = (image_rec_a+image_rec_b)/2.
    image_2d[k, *] = image_rec_ab[0:Matrix_size-1]
    print, k+1, ' out of ', y_size
    IF k EQ y_size/2-1 THEN iplot, imagex_rec_ab
  ENDFOR
  
    iimage, TRANSPOSE(imagex_2D), title='x-channel 2D image'
    iimage, TRANSPOSE(imagey_2D), title='y-channel 2D image'
    iimage, TRANSPOSE(image_2D), title='multi-channel 2D image'
;    iplot, image_2d

END

;@GJ, 2023/6/25, DICOM two contrast agent plots
PRO line_profile_plot_diff_conc


  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
  fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\NC_COMPASS\')
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF
  
  
  image_dcm = read_dicom(fn_dcm)
;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm, MIN=30000)
  it_fig = TVRD(TRUE=1)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'LSP\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig

  
  
  CURSOR, X, Y, /DEVICE
  width_in_pixel = 60./pixelSp[0]
    
  pixel_profile = image_dcm[X-width_in_pixel/2.:X+width_in_pixel/2., Y]
  line_cor = FINDGEN(width_in_pixel) *  pixelSp[0]
  
  ; Create the plot
  p = PLOT(line_cor, pixel_profile, 'k-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[600, 300], title=' ', xtitle='Sample Distance [mm]', ytitle='MPI Signal [A.U.]')
  result_signal_fig = p.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_lineSurfactPlot.png', result_signal_fig
  IF OBJ_VALID(p) THEN p.close
  
  

END

;@GJ, 2023/6/25, DICOM two contrast agent plots
PRO line_profile_plot_3_types


  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
  fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\NC_COMPASS\WangQianData\20230627 4geyizu de 102550130n\')
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF


  image_dcm = read_dicom(fn_dcm)
  ;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm);, MIN=30000)
  it_fig = TVRD(TRUE=1)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'LSP\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig



  CURSOR, X, Y, /DEVICE
  width_in_pixel = 54./pixelSp[0]

  pixel_profile = image_dcm[X-width_in_pixel/2.:X+width_in_pixel/2., Y]
  line_cor = FINDGEN(width_in_pixel) *  pixelSp[0]

  ; Create the plot
  p = PLOT(line_cor, pixel_profile, 'k-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[600, 300], title=' ', xtitle='Sample Distance [mm]', ytitle='MPI Signal [A.U.]')
  result_signal_fig = p.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_lineSurfactPlot.png', result_signal_fig
  IF OBJ_VALID(p) THEN p.close

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 1, XSIZE = width_in_pixel, YSIZE = rows
  TV, BYTSCL(image_dcm[X-width_in_pixel/2.:X+width_in_pixel/2., *]);, MIN=30000)
  it_fig = TVRD(TRUE=1)
  ;save the signal plot
  ;fitplot_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'LSP\'
  ;IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageCrop.png', it_fig

END

;@GJ, 2023/6/25, DICOM two contrast agent plots
PRO line_profile_plot_3_types_vivoTrax


  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
  fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\NC_COMPASS\WangQianData\20230627 4geyizu de 102550130n\')
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF


  image_dcm = read_dicom(fn_dcm)
  ;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm);, MIN=30000)
  it_fig = TVRD(TRUE=1)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'LSP\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig



  CURSOR, X, Y, /DEVICE
  width_in_pixel = 54./pixelSp[0]

  pixel_profile = image_dcm[X-width_in_pixel/6.:X+width_in_pixel*5./6., Y]
  line_cor = FINDGEN(width_in_pixel) *  pixelSp[0]

  ; Create the plot
  p = PLOT(line_cor, pixel_profile, 'k-2', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[600, 300], title=' ', xtitle='Sample Distance [mm]', ytitle='MPI Signal [A.U.]')
  result_signal_fig = p.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_lineSurfactPlot.png', result_signal_fig
  IF OBJ_VALID(p) THEN p.close

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 1, XSIZE = width_in_pixel, YSIZE = rows
  TV, BYTSCL(image_dcm[X-width_in_pixel/6.:X+width_in_pixel*5./6., *]);, MIN=30000)
  it_fig = TVRD(TRUE=1)
  ;save the signal plot
  ;fitplot_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'LSP\'
  ;IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageCrop.png', it_fig

END

;@GJ, 2023/6/26, reading the VSM file from Haohao Yan
;data=read_VSM_file()
;help, data
FUNCTION read_VSM_file, filename
  
  filename = 'F:\MPI_Tianjie\Nanorubble_MPI\VSM\N240_21-28mg.txt'
;  filename = 'F:\MPI_Tianjie\Nanorubble_MPI\VSM\N200_23-61mg.txt'


  ;GJ, 2023/6/26, nonzero gilename
  IF STRLEN(filename) GT 1 THEN BEGIN
    filename_info = FILE_INFO(filename)
    IF (filename_info.exists EQ 0) THEN BEGIN
      info=dialog_message('Wrong file name!', /information)
      RETURN, 0
    ENDIF
    IF (STRLEN(filename) EQ 0) THEN BEGIN
      info=dialog_message('No file name!', /information)
      RETURN, 0
    ENDIF
  ENDIF


  ;read the txt file 1st time to get the frequency and period
  OPENR, unit, filename, /GET_LUN
  header=''
  WHILE ~ EOF(unit) DO BEGIN
    READF, unit, header, FORMAT='(%"%s")'
    print, header
    IF STRPOS(header, '-') NE -1 AND STRPOS(header, 'mg') NE -1 AND N_ELEMENTS(weight) EQ 0 THEN BEGIN
      startpos = STRPOS(header, '-')
      endpos = STRPOS(header, 'mg')
      weight = DOUBLE(STRMID(header, startpos+1, endpos-startpos-1))/1000. ; weight in g
      print, 'weight: ', weight, ' g'
    ENDIF
    IF STRCMP(header, 'Field(G)', 8, /FOLD_CASE) THEN break;
  ENDWHILE
  
  ;start reading the data
  i=0;
  WHILE ~ EOF(unit) DO BEGIN
;    print, 'i = ', i
    data_string = ''
    READF, unit, data_string, FORMAT='(%"%s")'
;    print, 'data_string', data_string, ' end'
    temp_data=STRSPLIT(STRCOMPRESS(data_string), ' ', /extract)
    IF STRLEN(temp_data[0]) EQ 0 THEN break;
    IF i EQ 0 THEN BEGIN
      Field_array = [double(temp_data[0])]
      M_array = [double(temp_data[1])]
    ENDIF ELSE BEGIN
      Field_array = [Field_array, double(temp_data[0])]
      M_array = [M_array, double(temp_data[1])]
    ENDELSE
    i++
  ENDWHILE
  ;; close file and END execution
  FREE_LUN, unit
  
  IF N_ELEMENTS(weight) EQ 0 THEN weight = 1.0; g
  print, 'weight: ', weight, ' g'
  
;  iplot, Field_array * 0.1, M_array/weight, title='M-H curve from VSM', xtitle='Field [mT]', ytitle='Moment [emu/g]'
  
;  iplot, Field_array * 0.1, M_array/weight, xrange=[-150., 150.], title='M-H curve from VSM', xtitle='Field [mT]', ytitle='Moment [emu/g]'
  
;  iplot, Field_array * 0.1, M_array/weight, xrange=[-15., 15.], title='M-H curve from VSM', xtitle='Field [mT]', ytitle='Moment [emu/g]'

  data_array = [[Field_array * 0.1], [M_array/weight]]
  
  
  
  RETURN, TRANSPOSE(data_array)

END

;@GJ, 2023/6/30, plot the M-H curve of Wang Qian's data
PRO wq_sine_MH_plot_cell
  
  fn1 = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\20221020test1mps 10 25 50 130nm的铁粒子频谱信号\20221020\2022102012.txt'
  M_H_array1 = read_sine_signal_dat(fn1)
  H_curve1 = REFORM(M_H_array1[0,*])
  M_curve1 = REFORM(M_H_array1[1,*])
  ; Create the plot
  p1 = PLOT(H_curve1, M_curve1, 'r-3', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='M (A.U.)'); Create the plot
  p1.name = ' VivoTrax'
    
  fn2 = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\20221020test1mps 10 25 50 130nm的铁粒子频谱信号\20221020\2022102013.txt'
  M_H_array2 = read_sine_signal_dat(fn2)
  H_curve2 = REFORM(M_H_array2[0,*])
  M_curve2 = REFORM(M_H_array2[1,*])
  p2 = PLOT(H_curve2, M_curve2, '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p2.name = ' Perimag'
  
  fn3 = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\20221020test1mps 10 25 50 130nm的铁粒子频谱信号\20221020\2022102014.txt'
  M_H_array3 = read_sine_signal_dat(fn3)
  H_curve3 = REFORM(M_H_array3[0,*])
  M_curve3 = REFORM(M_H_array3[1,*])
  p3 = PLOT(H_curve3, M_curve3, '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p3.name = ' Synomag-D'

  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.55, 0.8])
  result00 = p1.CopyWindow(BORDER=10, RESOLUTION=300)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn1, /MARK_DIRECTORY)+'fitplot\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_curveCell.png', result00
;  IF OBJ_VALID(p1) THEN p1.close

  
  MIN_H1 = MIN(H_curve1, min_ind)
  H_curve1 = SHIFT(H_curve1, -min_ind)
  M_curve1 = SHIFT(M_curve1, -min_ind)
  Nele_MH = N_ELEMENTS(H_curve1)
  H_curve1_inc = H_curve1[0:Nele_MH/2-1]
  M_curve1_inc = M_curve1[0:Nele_MH/2-1]
  M_curve1_inc_der = (M_curve1[1:Nele_MH/2]-M_curve1[0:Nele_MH/2-1])/(H_curve1[1:Nele_MH/2]-H_curve1[0:Nele_MH/2-1])
  H_curve1_dec = H_curve1[Nele_MH/2:*]
  M_curve1_dec = M_curve1[Nele_MH/2:*]
  M_curve1_dec_der = (M_curve1[Nele_MH/2:*]-M_curve1[Nele_MH/2-1:Nele_MH-2])/(H_curve1[Nele_MH/2:*]-H_curve1[Nele_MH/2-1:Nele_MH-2])
  ;  iplot, H_curve1_dec, M_curve1_dec_der
  zero_H_curve1_inc = MIN(ABS(H_curve1_inc), zero_ind_inc)
  zero_H_curve1_dec = MIN(ABS(H_curve1_dec), zero_ind_dec)
  comb_H_curve1 = [H_curve1_inc[0:zero_ind_inc], H_curve1_dec[0:zero_ind_dec]]
  comb_M_curve1_der = [M_curve1_inc_der[0:zero_ind_inc], M_curve1_dec_der[0:zero_ind_dec]]
  comb_H_sort_ind = SORT(comb_H_curve1)
  N_ele_comb = N_ELEMENTS(comb_M_curve1_der)
  p4 = PLOT((comb_H_curve1(comb_H_sort_ind))[5:N_ele_comb-6], (comb_M_curve1_der(comb_H_sort_ind))[5:N_ele_comb-6], 'r-3', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[-9.,9.],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='dM/dH (A.U.)'); Create the plot
  p4.name = ' VivoTrax'
  
  MIN_H2 = MIN(H_curve2, min_ind)
  H_curve2 = SHIFT(H_curve2, -min_ind)
  M_curve2 = SHIFT(M_curve2, -min_ind)
  Nele_MH = N_ELEMENTS(H_curve2)
  H_curve2_inc = H_curve2[0:Nele_MH/2-1]
  M_curve2_inc = M_curve2[0:Nele_MH/2-1]
  M_curve2_inc_der = (M_curve2[1:Nele_MH/2]-M_curve2[0:Nele_MH/2-1])/(H_curve2[1:Nele_MH/2]-H_curve2[0:Nele_MH/2-1])
  H_curve2_dec = H_curve2[Nele_MH/2:*]
  M_curve2_dec = M_curve2[Nele_MH/2:*]
  M_curve2_dec_der = (M_curve2[Nele_MH/2:*]-M_curve2[Nele_MH/2-1:Nele_MH-2])/(H_curve2[Nele_MH/2:*]-H_curve2[Nele_MH/2-1:Nele_MH-2])
  ;  iplot, H_curve2_inc, M_curve2_inc_der
  zero_H_curve2_inc = MIN(ABS(H_curve2_inc), zero_ind_inc)
  zero_H_curve2_dec = MIN(ABS(H_curve2_dec), zero_ind_dec)
  comb_H_curve2 = [H_curve2_inc[0:zero_ind_inc], H_curve2_dec[0:zero_ind_dec]]
  comb_M_curve2_der = [M_curve2_inc_der[0:zero_ind_inc], M_curve2_dec_der[0:zero_ind_dec]]
  comb_H_sort_ind = SORT(comb_H_curve2)
  N_ele_comb = N_ELEMENTS(comb_M_curve2_der)
  p5 = PLOT((comb_H_curve2(comb_H_sort_ind))[5:N_ele_comb-6], (comb_M_curve2_der(comb_H_sort_ind))[5:N_ele_comb-6], '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p5.name = ' Perimag'
  
  MIN_H3 = MIN(H_curve3, min_ind)
  H_curve3 = SHIFT(H_curve3, -min_ind)
  M_curve3 = SHIFT(M_curve3, -min_ind)
  Nele_MH = N_ELEMENTS(H_curve3)
  H_curve3_inc = H_curve3[0:Nele_MH/2-1]
  M_curve3_inc = M_curve3[0:Nele_MH/2-1]
  M_curve3_inc_der = (M_curve3[1:Nele_MH/2]-M_curve3[0:Nele_MH/2-1])/(H_curve3[1:Nele_MH/2]-H_curve3[0:Nele_MH/2-1])
  H_curve3_dec = H_curve3[Nele_MH/2:*]
  M_curve3_dec = M_curve3[Nele_MH/2:*]
  M_curve3_dec_der = (M_curve3[Nele_MH/2:*]-M_curve3[Nele_MH/2-1:Nele_MH-2])/(H_curve3[Nele_MH/2:*]-H_curve3[Nele_MH/2-1:Nele_MH-2])
  ;  iplot, H_curve3_inc, M_curve3_inc_der
  zero_H_curve3_inc = MIN(ABS(H_curve3_inc), zero_ind_inc)
  zero_H_curve3_dec = MIN(ABS(H_curve3_dec), zero_ind_dec)
  comb_H_curve3 = [H_curve3_inc[0:zero_ind_inc], H_curve3_dec[0:zero_ind_dec]]
  comb_M_curve3_der = [M_curve3_inc_der[0:zero_ind_inc], M_curve3_dec_der[0:zero_ind_dec]]
  comb_H_sort_ind = SORT(comb_H_curve3)
  N_ele_comb = N_ELEMENTS(comb_M_curve3_der)
  p6 = PLOT((comb_H_curve3(comb_H_sort_ind))[5:N_ele_comb-6], (comb_M_curve3_der(comb_H_sort_ind))[5:N_ele_comb-6], '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p6.name = ' Synomag-D'
  
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  result00 = p4.CopyWindow(BORDER=10, RESOLUTION=300)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn1, /MARK_DIRECTORY)+'fitplot\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_der_curveCell.png', result00
  ;  IF OBJ_VALID(p4) THEN p1.close
  
  H_curve1_ind = FINDGEN(Nele_MH)/10. - 10.
  M_curve1_ind = INTERPOL(M_curve1_inc, H_curve1_inc, H_curve1_ind)
  M_curve1_ind_der = (M_curve1_ind[1:Nele_MH-2]-M_curve1_ind[0:Nele_MH-1])/(H_curve1_ind[1:Nele_MH-2]-H_curve1_ind[0:Nele_MH-1])
  p7 = PLOT(H_curve1_ind[5:Nele_MH-7], SMOOTH(M_curve1_ind_der[4:Nele_MH-7], 3), 'r-3', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[-15.,15.],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='M (A.U.)'); Create the plot
  p7.name = ' VivoTrax'
  H_curve2_ind = FINDGEN(Nele_MH)/10. - 10.
  M_curve2_ind = INTERPOL(M_curve2_inc, H_curve2_inc, H_curve2_ind)
  M_curve2_ind_der = (M_curve2_ind[1:Nele_MH-2]-M_curve2_ind[0:Nele_MH-1])/(H_curve2_ind[1:Nele_MH-2]-H_curve2_ind[0:Nele_MH-1])
  p8 = PLOT(H_curve2_ind[5:Nele_MH-7], SMOOTH(M_curve2_ind_der[4:Nele_MH-7], 3), '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p8.name = ' Perimag'
  H_curve3_ind = FINDGEN(Nele_MH)/10. - 10.
  M_curve3_ind = INTERPOL(M_curve3_inc, H_curve3_inc, H_curve3_ind)
  M_curve3_ind_der = (M_curve3_ind[1:Nele_MH-2]-M_curve3_ind[0:Nele_MH-1])/(H_curve3_ind[1:Nele_MH-2]-H_curve3_ind[0:Nele_MH-1])
  p9 = PLOT(H_curve3_ind[5:Nele_MH-7], SMOOTH(M_curve3_ind_der[4:Nele_MH-7], 3), '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p9.name = ' Synomag-D'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  result00 = p7.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_der_curveCell_ss.png', result00
  
  p10 = PLOT(H_curve1_ind[5:Nele_MH-7], BYTSCL(SMOOTH(M_curve1_ind_der[4:Nele_MH-7], 3)), 'r-3', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[-15.,15.],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='M (A.U.)'); Create the plot
  p10.name = ' VivoTrax'
  p11 = PLOT(H_curve2_ind[5:Nele_MH-7], BYTSCL(SMOOTH(M_curve2_ind_der[4:Nele_MH-7], 3)), '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p11.name = ' Perimag'
  p12 = PLOT(H_curve3_ind[5:Nele_MH-7], BYTSCL(SMOOTH(M_curve3_ind_der[4:Nele_MH-7], 3)), '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p12.name = ' Synomag-D'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  result00 = p10.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_der_curveCell_ssScale.png', result00
   

END

;@GJ, 2023/6/30, plot the M-H curve of Wang Qian's data
PRO wq_sine_MH_plot_SPIO

  fn1 = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\20221020test1mps 10 25 50 130nm的铁粒子频谱信号\20221020\202210202.txt'
  M_H_array1 = read_sine_signal_dat(fn1)
  H_curve1 = REFORM(M_H_array1[0,*])
  M_curve1 = REFORM(M_H_array1[1,*])
  ; Create the plot
  p1 = PLOT(H_curve1, M_curve1, 'r-3', AXIS_STYLE=1, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='M (A.U.)'); Create the plot
  p1.name = ' VivoTrax'
  
  fn2 = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\20221020test1mps 10 25 50 130nm的铁粒子频谱信号\20221020\202210203.txt'
  M_H_array2 = read_sine_signal_dat(fn2)
  H_curve2 = REFORM(M_H_array2[0,*])
  M_curve2 = REFORM(M_H_array2[1,*])
  p2 = PLOT(H_curve2, M_curve2, '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p2.name = ' Perimag'

  fn3 = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\20221020test1mps 10 25 50 130nm的铁粒子频谱信号\20221020\202210204.txt'
  M_H_array3 = read_sine_signal_dat(fn3)
  H_curve3 = REFORM(M_H_array3[0,*])
  M_curve3 = REFORM(M_H_array3[1,*])
  p3 = PLOT(H_curve3, M_curve3, '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p3.name = ' Synomag-D'
  
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.55, 0.8])
  result00 = p1.CopyWindow(BORDER=10, RESOLUTION=300)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn1, /MARK_DIRECTORY)+'fitplot\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_curveSPIO.png', result00
  ;  IF OBJ_VALID(p1) THEN p1.close
  
  
  MIN_H1 = MIN(H_curve1, min_ind)
  H_curve1 = SHIFT(H_curve1, -min_ind)
  M_curve1 = SHIFT(M_curve1, -min_ind)
  Nele_MH = N_ELEMENTS(H_curve1)
  H_curve1_inc = H_curve1[0:Nele_MH/2-1]
  M_curve1_inc = M_curve1[0:Nele_MH/2-1]
  M_curve1_inc_der = (M_curve1[1:Nele_MH/2]-M_curve1[0:Nele_MH/2-1])/(H_curve1[1:Nele_MH/2]-H_curve1[0:Nele_MH/2-1])
  H_curve1_dec = H_curve1[Nele_MH/2:*]
  M_curve1_dec = M_curve1[Nele_MH/2:*]
  M_curve1_dec_der = (M_curve1[Nele_MH/2:*]-M_curve1[Nele_MH/2-1:Nele_MH-2])/(H_curve1[Nele_MH/2:*]-H_curve1[Nele_MH/2-1:Nele_MH-2])
  ;  iplot, H_curve1_dec, M_curve1_dec_der
  zero_H_curve1_inc = MIN(ABS(H_curve1_inc), zero_ind_inc)
  zero_H_curve1_dec = MIN(ABS(H_curve1_dec), zero_ind_dec)
  comb_H_curve1 = [H_curve1_inc[0:zero_ind_inc], H_curve1_dec[0:zero_ind_dec]]
  comb_M_curve1_der = [M_curve1_inc_der[0:zero_ind_inc], M_curve1_dec_der[0:zero_ind_dec]]
  comb_H_sort_ind = SORT(comb_H_curve1)
  N_ele_comb = N_ELEMENTS(comb_M_curve1_der)
  p4 = PLOT((comb_H_curve1(comb_H_sort_ind))[5:N_ele_comb-6], (comb_M_curve1_der(comb_H_sort_ind))[5:N_ele_comb-6], 'r-3', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[-9.,9.],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='dM/dH (A.U.)'); Create the plot
  p4.name = ' VivoTrax'
  
  MIN_H2 = MIN(H_curve2, min_ind)
  H_curve2 = SHIFT(H_curve2, -min_ind)
  M_curve2 = SHIFT(M_curve2, -min_ind)
  Nele_MH = N_ELEMENTS(H_curve2)
  H_curve2_inc = H_curve2[0:Nele_MH/2-1]
  M_curve2_inc = M_curve2[0:Nele_MH/2-1]
  M_curve2_inc_der = (M_curve2[1:Nele_MH/2]-M_curve2[0:Nele_MH/2-1])/(H_curve2[1:Nele_MH/2]-H_curve2[0:Nele_MH/2-1])
  H_curve2_dec = H_curve2[Nele_MH/2:*]
  M_curve2_dec = M_curve2[Nele_MH/2:*]
  M_curve2_dec_der = (M_curve2[Nele_MH/2:*]-M_curve2[Nele_MH/2-1:Nele_MH-2])/(H_curve2[Nele_MH/2:*]-H_curve2[Nele_MH/2-1:Nele_MH-2])
  ;  iplot, H_curve2_inc, M_curve2_inc_der
  zero_H_curve2_inc = MIN(ABS(H_curve2_inc), zero_ind_inc)
  zero_H_curve2_dec = MIN(ABS(H_curve2_dec), zero_ind_dec)
  comb_H_curve2 = [H_curve2_inc[0:zero_ind_inc], H_curve2_dec[0:zero_ind_dec]]
  comb_M_curve2_der = [M_curve2_inc_der[0:zero_ind_inc], M_curve2_dec_der[0:zero_ind_dec]]
  comb_H_sort_ind = SORT(comb_H_curve2)
  N_ele_comb = N_ELEMENTS(comb_M_curve2_der)
  p5 = PLOT((comb_H_curve2(comb_H_sort_ind))[5:N_ele_comb-6], (comb_M_curve2_der(comb_H_sort_ind))[5:N_ele_comb-6], '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p5.name = ' Perimag'
  
  MIN_H3 = MIN(H_curve3, min_ind)
  H_curve3 = SHIFT(H_curve3, -min_ind)
  M_curve3 = SHIFT(M_curve3, -min_ind)
  Nele_MH = N_ELEMENTS(H_curve3)
  H_curve3_inc = H_curve3[0:Nele_MH/2-1]
  M_curve3_inc = M_curve3[0:Nele_MH/2-1]
  M_curve3_inc_der = (M_curve3[1:Nele_MH/2]-M_curve3[0:Nele_MH/2-1])/(H_curve3[1:Nele_MH/2]-H_curve3[0:Nele_MH/2-1])
  H_curve3_dec = H_curve3[Nele_MH/2:*]
  M_curve3_dec = M_curve3[Nele_MH/2:*]
  M_curve3_dec_der = (M_curve3[Nele_MH/2:*]-M_curve3[Nele_MH/2-1:Nele_MH-2])/(H_curve3[Nele_MH/2:*]-H_curve3[Nele_MH/2-1:Nele_MH-2])
  ;  iplot, H_curve3_inc, M_curve3_inc_der
  zero_H_curve3_inc = MIN(ABS(H_curve3_inc), zero_ind_inc)
  zero_H_curve3_dec = MIN(ABS(H_curve3_dec), zero_ind_dec)
  comb_H_curve3 = [H_curve3_inc[0:zero_ind_inc], H_curve3_dec[0:zero_ind_dec]]
  comb_M_curve3_der = [M_curve3_inc_der[0:zero_ind_inc], M_curve3_dec_der[0:zero_ind_dec]]
  comb_H_sort_ind = SORT(comb_H_curve3)
  N_ele_comb = N_ELEMENTS(comb_M_curve3_der)
  p6 = PLOT((comb_H_curve3(comb_H_sort_ind))[5:N_ele_comb-6], (comb_M_curve3_der(comb_H_sort_ind))[5:N_ele_comb-6], '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p6.name = ' Synomag-D'
  
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  result00 = p4.CopyWindow(BORDER=10, RESOLUTION=300)
  ;save the signal plot
  fitplot_dir = FILE_DIRNAME(fn1, /MARK_DIRECTORY)+'fitplot\'
  IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_der_curveSPIO.png', result00
  ;  IF OBJ_VALID(p4) THEN p1.close
  
  H_curve1_ind = FINDGEN(Nele_MH)/10. - 10.
  M_curve1_ind = INTERPOL(M_curve1_inc, H_curve1_inc, H_curve1_ind)
  M_curve1_ind_der = (M_curve1_ind[1:Nele_MH-2]-M_curve1_ind[0:Nele_MH-1])/(H_curve1_ind[1:Nele_MH-2]-H_curve1_ind[0:Nele_MH-1])
  p7 = PLOT(H_curve1_ind[5:Nele_MH-7], SMOOTH(M_curve1_ind_der[4:Nele_MH-7], 3), 'r-3', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[-15.,15.],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='M (A.U.)'); Create the plot
  p7.name = ' VivoTrax'
  H_curve2_ind = FINDGEN(Nele_MH)/10. - 10.
  M_curve2_ind = INTERPOL(M_curve2_inc, H_curve2_inc, H_curve2_ind)
  M_curve2_ind_der = (M_curve2_ind[1:Nele_MH-2]-M_curve2_ind[0:Nele_MH-1])/(H_curve2_ind[1:Nele_MH-2]-H_curve2_ind[0:Nele_MH-1])
  p8 = PLOT(H_curve2_ind[5:Nele_MH-7], SMOOTH(M_curve2_ind_der[4:Nele_MH-7], 3), '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p8.name = ' Perimag'
  H_curve3_ind = FINDGEN(Nele_MH)/10. - 10.
  M_curve3_ind = INTERPOL(M_curve3_inc, H_curve3_inc, H_curve3_ind)
  M_curve3_ind_der = (M_curve3_ind[1:Nele_MH-2]-M_curve3_ind[0:Nele_MH-1])/(H_curve3_ind[1:Nele_MH-2]-H_curve3_ind[0:Nele_MH-1])
  p9 = PLOT(H_curve3_ind[5:Nele_MH-7], SMOOTH(M_curve3_ind_der[4:Nele_MH-7], 3), '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p9.name = ' Synomag-D'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  result00 = p7.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_der_curveSPIO_ss.png', result00
  
  p10 = PLOT(H_curve1_ind[5:Nele_MH-7], BYTSCL(SMOOTH(M_curve1_ind_der[4:Nele_MH-7], 3)), 'r-3', AXIS_STYLE=1, FONT_SIZE=9, XRANGE=[-15.,15.],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[700, 700], title=' ', xtitle='H (mT)', ytitle='M (A.U.)'); Create the plot
  p10.name = ' VivoTrax'
  p11 = PLOT(H_curve2_ind[5:Nele_MH-7], BYTSCL(SMOOTH(M_curve2_ind_der[4:Nele_MH-7], 3)), '.g-3', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p11.name = ' Perimag'
  p12 = PLOT(H_curve3_ind[5:Nele_MH-7], BYTSCL(SMOOTH(M_curve3_ind_der[4:Nele_MH-7], 3)), '.-3', COLOR='#7F007F', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p12.name = ' Synomag-D'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  result00 = p10.CopyWindow(BORDER=10, RESOLUTION=300)
  WRITE_PNG, fitplot_dir+FILE_BASENAME(fn1, '.txt')+'_MH_der_curveSPIO_ssScale.png', result00
   
   
END


;@GJ, 2023/7/8, resolution versus K subset by averaging
PRO resolution_vs_K_average


  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
  ;fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\NC_COMPASS\')
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF


  image_dcm = read_dicom(fn_dcm)
  ;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm, MIN=3000)
  it_fig = TVRD(TRUE=1)
  ResoK_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'ResoK\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN FILE_MKDIR, ResoK_dir
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig

;  CURSOR, X, Y, /DEVICE
  X = 127;
  print, 'X = ', X, ' (cols:', cols, ')'
  image_left = DOUBLE(image_dcm[0:X, *])/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *])/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(FLOOR(image_left))
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageLeft.png', it_fig
  
  flag_k=0;
  FOR i=0, X-1 DO BEGIN
    FOR j=0, rows-1 DO BEGIN
      IF image_left[i, j] GT 1 THEN BEGIN
        IF flag_k EQ 0 THEN BEGIN
          pixel_array = DBLARR(FLOOR(image_left[i, j]), 2)
          pixel_array[*,0] = i
          pixel_array[*,1] = j
        ENDIF ELSE BEGIN
          temp_pixel_array = DBLARR(FLOOR(image_left[i, j]), 2)
          temp_pixel_array[*,0] = i
          temp_pixel_array[*,1] = j
          pixel_array =[pixel_array, temp_pixel_array]       
        ENDELSE
        flag_k++
        print, 'flag_k = ', flag_k
      ENDIF
    ENDFOR
  ENDFOR
;  print, 'flag_k = ', k
  Pixels_N = N_ELEMENTS(pixel_array[*,0])
  
  K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 1000]
  FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
    K = K_array[k_ind]
    new_image = image_left * 0.
    FOR i=0, 500 DO BEGIN
      seed = !NULL
      d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
      x_ave = MEAN(pixel_array[d1, 0])
      y_ave = MEAN(pixel_array[d1, 1])
      new_image[x_ave, y_ave] = 1.
    ENDFOR
    ;  iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
    ; R, G, B color table - make everything red
    r = REPLICATE(255b,256)
    r[0:112] = 0
    g = REPLICATE(0b,256)
    b = g
    WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_ResoK'+STRTRIM(FIX(K), 1)+'.png', BYTSCL(new_image);, r, g, b
  ENDFOR

END

;@GJ, 2023/8/28, analyzing the Lion Eye image
PRO exp_transform_Lion_eye_image
  
  fn_dir = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\'
  fn_name = 'Mammo4im_rec_Lion.png'
  fn = fn_dir + fn_name
  image_sim = DOUBLE(read_png(fn))
  image_sim[*, 0:326] = 0.
  iimage, BYTSCL(image_sim), title='Original'
  
  rows_S = N_ELEMENTS(image_sim[0,*])
  cols_S = N_ELEMENTS(image_sim[*,0])
  pixelSp = [1.0, 1.0]
  
  tao = -MAX(image_sim) / ALOG(0.001)
  sepa = 10.
  thres = 18.
  exp_image = image_sim * 0.
  FOR i=sepa, cols_S-sepa-1 DO BEGIN
    FOR j=sepa, rows_S-sepa-1 DO BEGIN
      IF image_sim[i, j] GT thres/6. THEN BEGIN
        gradient_x_left = image_sim[i, j]-image_sim[i-sepa, j]
        gradient_x_right = image_sim[i, j]-image_sim[i+sepa, j]
        IF gradient_x_left LT 0 AND gradient_x_right LT 0 THEN BEGIN
          exp_image[i,j] = EXP(-image_sim[i, j]/tao)
        ENDIF
      ENDIF
    ENDFOR
  ENDFOR
 
 iimage, BYTSCL(exp_image), title='Exponential'
  
END

;@GJ, 2023/8/23, phantom shape simulation
PRO phantom_shape_simulation

  fn_dir = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\'
  ;fn_name = 'S_shape.png'
;  fn_name = 'letter32.png'
;  fn_name = 'SUA_letter24.png'
;  fn_name = 'Mammo2.png'
  fn_name = 'Mammo4.png'
;  fn_name = 'Mammo5.jpg'
  fn_name = 'MRA000.png'
  fn = fn_dir + fn_name
  
  IF QUERY_PNG(fn) EQ 1 THEN BEGIN
    image_temp = read_png(fn)
    image_sim_temp = DOUBLE(REFORM(image_temp[0,*,*]))
    image_sim = MAX(image_sim_temp) - image_sim_temp
    IF STRCMP(fn_name, 'MRA', 3) THEN image_sim = DOUBLE(image_sim LT 192) * 255.
    rows_S = N_ELEMENTS(image_sim[0,*])
    cols_S = N_ELEMENTS(image_sim[*,0])
    pixelSp = [1.0, 1.0]
  ENDIF ELSE BEGIN
    return
  ENDELSE
  ;@GJ, modify the phantom data
;  image_sim *= 0.
;  image_sim[rows_S/2-3:rows_S/2+3, cols_S/2-30:cols_S/2+30] = 255
;  image_sim[rows_S/2-30:rows_S/2+30, cols_S/2-3:cols_S/2+3] = 255
  
  ;display the image phantom
  WINDOW, 0, XSIZE = cols_S, YSIZE = rows_S
  TV, BYTSCL(image_sim)  
  
  ;@GJ, 2023/8/24, eye of Sauron
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\005_0626YHHprobe_2023-06-26_15-22-01\combined_0000_combined_image_0000.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF
  IF QUERY_DICOM(fn_dcm) EQ 1 THEN BEGIN
    image_dcm = read_dicom(fn_dcm)
    ;  iimage, image_dcm, title='original'
    obj = OBJ_NEW('IDLffDICOM', fn_dcm)
    pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
    pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
    ; Get the row & column size of the image(s):
    temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
    temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

    OBJ_DESTROY, obj
  ENDIF ELSE BEGIN
    return
  ENDELSE

  X = 188.;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)

  image_sim_Sauron = CONVOL(image_sim, image_left, /EDGE_ZERO, /NORMALIZE)
;  image_sim_Sauron = CONVOL_FFT(image_sim, image_left)
  iimage, image_sim_Sauron, title='Eye of Sauron'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Sauron.png', BYTSCL(image_sim_Sauron)

  im_ori = image_sim
  im_rec = image_sim_Sauron
  MSE=calMSE(im_ori, im_rec);72.1795
  print, 'Sauron Eye MSE: ', MSE
  pSNR=calPSNR(im_ori, MSE);68.0337
  print, 'Sauron Eye pSNR: ', pSNR, ' dB'
  mssim = SSIM(im_ori, im_rec)
  print, 'Sauron Eye SSIM: ', mssim
  
;  ;get the image background noise
;  ;do pixel searching
;  hit_map = DBLARR(cols_S, rows_S) * 0.
;  min_sim = MAX(image_sim_Sauron)
;  image_sim_threshold = 0.5 * MAX(image_sim_Sauron)
;  twin_distance = 10
;  FOR k=1, twin_distance DO BEGIN
;    FOR i=twin_distance, cols_S-twin_distance-1 DO BEGIN
;      FOR j=twin_distance/2., rows_S-twin_distance/2.-1 DO BEGIN
;        IF image_sim_Sauron[i, j] GT image_sim_threshold*0.8 THEN BEGIN
;          gradient_x_left = image_sim_Sauron[i, j]-image_sim_Sauron[i-twin_distance, j]
;          gradient_x_right = image_sim_Sauron[i, j]-image_sim_Sauron[i+twin_distance, j]
;          IF gradient_x_left LT 0 AND gradient_x_right LT 0 THEN BEGIN
;            hit_map[i, j] = 1.
;            IF min_sim GT image_sim_Sauron[i, j] THEN min_sim = image_sim_Sauron[i, j]
;          ENDIF
;        ENDIF
;      ENDFOR
;    ENDFOR
;  ENDFOR
;  image_sim_threshold = 0.5 * MAX(image_sim_Sauron)
;  min_Sauron = MIN(image_sim_Sauron[where(image_sim_Sauron GT image_sim_threshold)])
;  max_Sauron = MAX(image_sim_Sauron[where(image_sim_Sauron GT image_sim_threshold)])
;  tao = -(max_Sauron-min_Sauron) / ALOG(0.001)
;  exp_Sauron = image_sim_Sauron * 0.
;  FOR i=0, cols_S-1 DO BEGIN
;    FOR j=0, rows_S-1 DO BEGIN
;      IF image_sim_Sauron[i,j] GT min_Sauron THEN exp_Sauron[i,j] = EXP(-(image_sim_Sauron[i,j]-min_Sauron)/tao)
;    ENDFOR
;  ENDFOR
;  iimage, image_filter, title='filter image'
  
;  max_Sauron = MAX(image_sim_Sauron, max_Sauron_ind)
;  tao = -(max_Sauron-min_sim) / ALOG(0.001)
;  exp_Sauron = image_sim_Sauron * 0.
;  FOR i=twin_distance, cols_S-twin_distance-1 DO BEGIN
;    FOR j=twin_distance/2., rows_S-twin_distance/2.-1 DO BEGIN
;      IF image_sim_Sauron[i,j]-min_sim THEN exp_Sauron[i,j] = EXP(-(image_sim_Sauron[i,j]-min_sim)/tao)
;    ENDFOR
;  ENDFOR
;  iimage, exp_Sauron, title='Exp Eye of Sauron'
;  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_SauronExp.png', BYTSCL(exp_Sauron)
  
;  image_sim_Sauron_vert = CONVOL(image_sim, ROTATE(image_left, 1), /EDGE_ZERO, /NORMALIZE)
;  iimage, image_sim_Sauron_vert, title='Eye of Sauron (Rot 90)'
;  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Sauron_vert.png', BYTSCL(image_sim_Sauron_vert)
;
;  im_ori = image_sim
;  im_rec = image_sim_Sauron_vert
;  MSE=calMSE(im_ori, im_rec);72.1795
;  print, 'Sauron_vert Eye MSE: ', MSE
;  pSNR=calPSNR(im_ori, MSE);68.0337
;  print, 'Sauron_vert Eye pSNR: ', pSNR, ' dB'
;  mssim = SSIM(im_ori, im_rec)
;  print, 'Sauron_vert Eye SSIM: ', mssim
;
;  max_Sauron_vert = MAX(image_sim_Sauron_vert, max_Sauron_vert_ind)
;  tao_vert = -max_Sauron_vert / ALOG(0.001)
;  exp_Sauron_vert = EXP(-image_sim_Sauron_vert/tao_vert)
;  exp_Sauron_vert[where(image_sim_Sauron_vert LT 8.8, /null)] = 0.
;  iimage, exp_Sauron_vert, title='Exp Eye of Sauron'
;  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_ExpSauron_vert.png', BYTSCL(exp_Sauron_vert)

;  ;half Sauron eye
;  X = 99.;
;  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
;  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
;  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
;  WINDOW, 1, XSIZE = X, YSIZE = rows
;  TV, BYTSCL(image_left)
;
;  image_sim_Half_Sauron = CONVOL(image_sim, image_left, /EDGE_ZERO, /NORMALIZE)
;  iimage, image_sim_Half_Sauron, title='Half Eye of Sauron'
;  write_png, fn_dir + 'im_rec_Half_Sauron.png', BYTSCL(image_sim_Half_Sauron)
;  
;  im_ori = image_sim
;  im_rec = image_sim_Half_Sauron
;  MSE=calMSE(im_ori, im_rec);72.1795
;  print, 'Half Sauron Eye MSE: ', MSE
;  pSNR=calPSNR(im_ori, MSE);68.0337
;  print, 'Half Sauron Eye pSNR: ', pSNR, ' dB'
;  mssim = SSIM(im_ori, im_rec)
;  print, 'Half Sauron Eye SSIM: ', mssim
  
  ;@GJ, lion eye
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF
  IF QUERY_DICOM(fn_dcm) EQ 1 THEN BEGIN
    image_dcm = read_dicom(fn_dcm)
    ;  iimage, image_dcm, title='original'
    obj = OBJ_NEW('IDLffDICOM', fn_dcm)
    pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
    pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
    ; Get the row & column size of the image(s):
    temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
    temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

    OBJ_DESTROY, obj
  ENDIF ELSE BEGIN
    return
  ENDELSE
  
;  X = 127;244; 0;127;
;  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
;  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
;  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
;  WINDOW, 1, XSIZE = X, YSIZE = rows
;  TV, BYTSCL(image_left)
;
;  image_sim_Half_Lion = CONVOL(image_sim, image_left, /EDGE_ZERO, /NORMALIZE)
;  iimage, image_sim_Half_Lion, title='Half Lion Eye'
;  write_png, fn_dir + 'im_rec_Half_Lion.png', BYTSCL(image_sim_Half_Lion)
;  
;  im_ori = image_sim
;  im_rec = image_sim_Half_Lion
;  MSE=calMSE(im_ori, im_rec);72.1795
;  print, 'Half Lion Eye MSE: ', MSE
;  pSNR=calPSNR(im_ori, MSE);68.0337
;  print, 'Half Lion Eye pSNR: ', pSNR, ' dB'
;  mssim = SSIM(im_ori, im_rec)
;  print, 'Half Lion Eye SSIM: ', mssim

  X = 244; 0;127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)

  image_sim_Lion = CONVOL_FFT(image_sim, image_left);, /EDGE_ZERO, /NORMALIZE)
  iimage, image_sim_Lion, title='Lion Eye'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion.png', BYTSCL(image_sim_Lion)

  im_ori = image_sim
  im_rec = image_sim_Lion
  MSE=calMSE(im_ori, im_rec);72.1795
  print, 'Lion Eye MSE: ', MSE
  pSNR=calPSNR(im_ori, MSE);68.0337
  print, 'Lion Eye pSNR: ', pSNR, ' dB'
  mssim = SSIM(im_ori, im_rec)
  print, 'Lion Eye SSIM: ', mssim
  
  image_sim_Lion_vert = CONVOL_FFT(image_sim, ROTATE(image_left, 1));, /EDGE_ZERO, /NORMALIZE)
  iimage, image_sim_Lion_vert, title='Eye of Lion (Rot 90)'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion_vert.png', BYTSCL(image_sim_Lion_vert)
  iimage, image_sim_Lion_vert+image_sim_Lion, title='Eye of Lion (Added)'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion_Added.png', BYTSCL(image_sim_Lion_vert+image_sim_Lion)
  iimage, image_sim_Lion_vert-image_sim_Lion, title='Eye of Lion (Subtracted)'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion_Subtracted.png', BYTSCL(image_sim_Lion_vert-image_sim_Lion)
  
  image_sim_Lion_comb = image_sim_Lion * 0.
  FOR i=0, cols_S-1 DO FOR j=0, rows_S-1 DO image_sim_Lion_comb[i, j] = MIN([image_sim_Lion[i,j],image_sim_Lion_vert[i,j]])
  iimage, image_sim_Lion_comb, title='Eye of Lion (Minimized)'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion_Minimized.png', BYTSCL(image_sim_Lion_comb)
  
  ;standard MNP
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\004_duizhaotu\duizhaotu.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF
  IF QUERY_DICOM(fn_dcm) EQ 1 THEN BEGIN
    image_dcm = read_dicom(fn_dcm)
    ;  iimage, image_dcm, title='original'
    obj = OBJ_NEW('IDLffDICOM', fn_dcm)
    pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
    pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
    ; Get the row & column size of the image(s):
    temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
    temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

    OBJ_DESTROY, obj
  ENDIF ELSE BEGIN
    return
  ENDELSE

  X = 200;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 2, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)

  image_sim_MNP = CONVOL_FFT(image_sim, image_left);, /EDGE_ZERO, /NORMALIZE)
  iimage, image_sim_MNP, title='Standard MNP'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_MNP.png', BYTSCL(image_sim_MNP)
  
  im_ori = image_sim
  im_rec = image_sim_MNP
  MSE=calMSE(im_ori, im_rec);72.1795
  print, 'MNP MSE: ', MSE
  pSNR=calPSNR(im_ori, MSE);68.0337
  print, 'MNP pSNR: ', pSNR, ' dB'
  mssim = SSIM(im_ori, im_rec)
  print, 'MNP SSIM: ', mssim
END

;@GJ, 2023/8/23, phantom shape simulation with Lion Eyes
PRO phantom_shape_simulation_lion_left_right_eye

  fn_dir = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\'
  ;fn_name = 'S_shape.png'
  fn_name = 'letter32.png'
  fn = fn_dir + fn_name

  IF QUERY_PNG(fn) EQ 1 THEN BEGIN
    image_temp = read_png(fn)
    image_sim = DOUBLE(REFORM(image_temp[3,*,*]))
    rows_S = N_ELEMENTS(image_sim[0,*])
    cols_S = N_ELEMENTS(image_sim[*,0])
    pixelSp = [1.0, 1.0]
  ENDIF ELSE BEGIN
    return
  ENDELSE
  ;@GJ, modify the phantom data
;  image_sim *= 0.
  ;  image_sim[rows_S/2-3:rows_S/2+3, cols_S/2-30:cols_S/2+30] = 255
  ;image_sim[rows_S/2-30:rows_S/2+30, cols_S/2-3:cols_S/2+3] = 255
;  image_sim[rows_S/2-3:rows_S/2+3, cols_S/2-3:cols_S/2+3] = 255
  
  ;display the image phantom
  WINDOW, 0, XSIZE = cols_S, YSIZE = rows_S
  TV, BYTSCL(image_sim)

  ;@GJ, lion eye
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF
  IF QUERY_DICOM(fn_dcm) EQ 1 THEN BEGIN
    image_dcm = read_dicom(fn_dcm)
    ;  iimage, image_dcm, title='original'
    obj = OBJ_NEW('IDLffDICOM', fn_dcm)
    pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
    pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
    ; Get the row & column size of the image(s):
    temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
    temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

    OBJ_DESTROY, obj
  ENDIF ELSE BEGIN
    return
  ENDELSE

  X = 127;244; 0;127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = image_dcm;DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_left[X:*,*] = 0.
  image_right = image_dcm;DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  image_right[0:X,*] = 0.
  WINDOW, 1, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_left)

  image_sim_left = CONVOL_FFT(image_sim, image_left);, /EDGE_ZERO, /NORMALIZE)
  iimage, image_sim_left, title='Left Lion Eye'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_left_Lion.png', BYTSCL(image_sim_left)

  im_ori = image_sim
  im_rec = image_sim_left
  MSE=calMSE(im_ori, im_rec);72.1795
  print, 'Left Lion Eye MSE: ', MSE
  pSNR=calPSNR(im_ori, MSE);68.0337
  print, 'Left Lion Eye pSNR: ', pSNR, ' dB'
  mssim = SSIM(im_ori, im_rec)
  print, 'Left Lion Eye SSIM: ', mssim

  image_sim_right = CONVOL_FFT(image_sim, image_right);, /EDGE_ZERO, /NORMALIZE)
  iimage, image_sim_right, title='Right Lion Eye'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_right_Lion.png', BYTSCL(image_sim_right)

  im_ori = image_sim
  im_rec = image_sim_right
  MSE=calMSE(im_ori, im_rec);72.1795
  print, 'Right Lion Eye MSE: ', MSE
  pSNR=calPSNR(im_ori, MSE);68.0337
  print, 'Right Lion Eye pSNR: ', pSNR, ' dB'
  mssim = SSIM(im_ori, im_rec)
  print, 'Right Lion Eye SSIM: ', mssim
 
  iimage, image_sim_left - image_sim_right, title='Combined Line Eye 1'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion_cobmined1.png', BYTSCL(image_sim_left - image_sim_right)
  
  iimage, image_sim_left + image_sim_right, title='Combined Line Eye 2'
  write_png, fn_dir + FILE_BASENAME(fn_name, '.png') + 'im_rec_Lion_cobmined2.png', BYTSCL(image_sim_left + image_sim_right)

END


;@GJ, 2023/7/8, resolution versus K subset by maximal likehood
;@GJ, 2023/7/11, standard deviation, error bar plotting and sigma fitting curve is added
;@GJ, 2023/8/23, add the images for simulation
PRO resolution_vs_K_max_sigma_plot_single_center


  ;loading DICOM image and get informations
;    fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
  fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
;  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF 
  
;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\smallMethod_data\002.png'
;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\smallMethod_data\001.png'
;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\SUT_data\C20230726\005.png'
;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\SUT_data\line20230726\005.png'
;   fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\im_rec_MNP.png'
;   fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\im_rec_Lion.png'
;   fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\im_rec_Sauron.png'
;   fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\im_rec_SauronVertical.png'
   fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\im_rec_Half_Lion.png'
  IF QUERY_DICOM(fn_dcm) EQ 1 THEN BEGIN
    image_dcm = read_dicom(fn_dcm)
    ;  iimage, image_dcm, title='original'
    obj = OBJ_NEW('IDLffDICOM', fn_dcm)
    pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
    pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
    ; Get the row & column size of the image(s):
    temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
    temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

    OBJ_DESTROY, obj
  ENDIF
  
  IF QUERY_PNG(fn_dcm) EQ 1 THEN BEGIN
    image_temp = read_png(fn_dcm)
    IF (SIZE(image_temp))[0] EQ 2 THEN image_dcm = DOUBLE(image_temp) ELSE image_dcm = DOUBLE(REFORM(image_temp[0,*,*]))
    rows = N_ELEMENTS(image_dcm[0,*])
    cols = N_ELEMENTS(image_dcm[*,0])
    pixelSp = [1.0, 1.0]
  ENDIF

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm)
  it_fig = TVRD(TRUE=1)
  
  separa = 0;2;5;2;0;1;2.;8.;10.;5.;0.;10.
  ResoK_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'ResoK_max'+STRTRIM(FIX(separa), 1)+'\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN FILE_MKDIR, ResoK_dir
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig
  
  ;@GJ, adding two together
;  image_dcm += SHIFT(image_dcm, separa, 0.)
  
  ;  CURSOR, X, Y, /DEVICE
  X = cols-1;0;127;
;  X = 127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageLeft.png', it_fig


  Pixels_N = N_ELEMENTS(image_left)
  print, '# pixels: ', pixels_N

  K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000]
  ;  K_array = [100]
  ;@GJ, 2023/7/11, the rounds of calculagint the standard deviation
  N_rounds = 10.;00.
  sigma_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  x_stddev_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  y_stddev_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  x_variance_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  y_variance_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  x_skewness_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  y_skewness_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  x_kurtosis_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  y_kurtosis_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  localization_image_setting = DBLARR(N_ELEMENTS(K_array), N_rounds, N_ELEMENTS(image_left[*,0]), N_ELEMENTS(image_left[0,*]))
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
      K = K_array[k_ind]
      new_image = image_left * 0.
      rep_N = 500000
      loc_x = DBLARR(rep_N) * 0.
      loc_y = DBLARR(rep_N) * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(image_left[d1], maxind)
        new_image[d1[maxind]] += 1.
        ind = ARRAY_INDICES(new_image, d1[maxind])
        loc_x[i] = ind[0] * pixelSp[0]
        loc_y[i] = ind[1] * pixelSp[1]
      ENDFOR
      ;calculate the covariance of the shrinked images
      sigma_array[k_ind, i_round] = SQRT(0.5 * ABS(C_CORRELATE(loc_x, loc_y, 0, /COVARIANCE, /DOUBLE)))
      x_stddev_array[k_ind, i_round] = STDDEV(loc_x)
      y_stddev_array[k_ind, i_round] = STDDEV(loc_y)
      x_result = MOMENT(loc_x)
      y_result = MOMENT(loc_y)
      x_variance_array[k_ind, i_round] = x_result[1]
      y_variance_array[k_ind, i_round] = y_result[1]
      x_skewness_array[k_ind, i_round] = x_result[2]
      y_skewness_array[k_ind, i_round] = y_result[2]
      x_kurtosis_array[k_ind, i_round] = x_result[3]
      y_kurtosis_array[k_ind, i_round] = y_result[3]
      localization_image_setting[k_ind, i_round, *, *] = new_image
      ;    window, 3
      ;    plot, loc_x, loc_y
;        iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
      ; R, G, B color table - make everything red
      IF i_round EQ 0 THEN BEGIN
        r = REPLICATE(255b,256)
        r[0:112] = 0
        g = REPLICATE(0b,256)
        b = g
        WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_ResoK'+STRTRIM(FIX(K), 1)+'.png', BYTSCL(new_image);, r, g, b
      ENDIF
    ENDFOR
    IF i_round EQ 0 THEN BEGIN
;      print, 'sigma_array: ', sigma_array
;      iplot, K_array[1:*], sigma_array[1:*, i_round], /xlog, ytitle='sigma (mm)', xtitle='K', title='Localization precision vs K (1 round)'
      iplot, K_array[1:*], sigma_array[1:*, i_round], yrange=[-1,1], color='red', ytitle='sigma (mm)', xtitle='K', title='Localization precision vs K ('+STRTRIM(FIX(separa), 1)+')'
;      iplot, K_array[1:*], x_stddev_array[1:*, i_round], color='green', /overplot
;      iplot, K_array[1:*], y_stddev_array[1:*, i_round], color='blue', /overplot
      iplot, K_array[1:*], x_stddev_array[1:*, i_round]-y_stddev_array[1:*, i_round], color='green', /overplot
;      iplot, K_array[1:*], x_variance_array[1:*, i_round], color='green', psym=1, linestyle=1, /overplot
;      iplot, K_array[1:*], y_variance_array[1:*, i_round], color='blue', psym=1, linestyle=1, /overplot
      iplot, K_array[1:*], x_variance_array[1:*, i_round]-y_variance_array[1:*, i_round], color='green', psym=1, linestyle=1, /overplot
;      iplot, K_array[1:*], x_skewness_array[1:*, i_round], color='green', psym=2, linestyle=2, /overplot
;      iplot, K_array[1:*], y_skewness_array[1:*, i_round], color='blue', psym=2, linestyle=2, /overplot
      iplot, K_array[1:*], x_skewness_array[1:*, i_round]-y_skewness_array[1:*, i_round], color='green', psym=2, linestyle=2, /overplot
;      iplot, K_array[1:*], x_kurtosis_array[1:*, i_round], color='green', psym=4, linestyle=3, /overplot
;      iplot, K_array[1:*], y_kurtosis_array[1:*, i_round], color='blue', psym=4, linestyle=3, /overplot
      iplot, K_array[1:*], x_kurtosis_array[1:*, i_round]-y_kurtosis_array[1:*, i_round], color='green', psym=4, linestyle=3, /overplot
    ENDIF
  ENDFOR
  
  ;analyz the density versus K
  w_size = 5;20
  flag = 0LL
  rec_image = image_left * 0.
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR i = 0, N_ELEMENTS(image_left[*,0])-w_size-1, w_size DO BEGIN
      FOR j = 0, N_ELEMENTS(image_left[0,*])-w_size-1, w_size DO BEGIN
        density_array = K_array * 0.
        FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
          density_array[k_ind] = TOTAL(localization_image_setting[k_ind, i_round, i:i+w_size-1, j:j+w_size-1])
        ENDFOR
;        IF flag EQ 0 THEN BEGIN
;          iplot, K_array, density_array/density_array[0], /xlog, /ylog, title='pixel density vs K'
;        ENDIF ELSE BEGIN
;          iplot, K_array, density_array/density_array[0], /overplot
;        ENDELSE
        flag++
;        print, 'flag: ', flag
        max_density = MAX(density_array, maxInd)
        rec_image[i:i+w_size-1, j:j+w_size-1] = maxInd
      ENDFOR
    ENDFOR
  ENDFOR
  
  iimage, rec_image, title='reconstructed image'
  ;@GJ, 2023/7/11, calculate the mean and standardev
  average_sigma = DBLARR(N_ELEMENTS(K_array))
  stddev_sigma = DBLARR(N_ELEMENTS(K_array))
  FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
    average_sigma[k_ind] = MEAN(sigma_array[k_ind, *])
    stddev_sigma[k_ind] = STDDEV(sigma_array[k_ind, *])    
  ENDFOR
  iplot, ALOG(K_array), ALOG(average_sigma), ytitle='ALOG sigma (mm)', xtitle='ALOG K', title='Localization precision vs K ('+STRTRIM(FIX(separa), 1)+')'
;  iplot, K_array, stddev_sigma, ytitle='stddev sigma (mm)', xtitle='K', xrange=[0, 200], title='Localization precision vs K'
  
  result = LINFIT(ALOG(K_array), ALOG(average_sigma), MEASURE_ERRORS=ALOG(stddev_sigma))
  print, 'fitting result: ', result
  sigma_fit1 = EXP(result[0]) * (K_array^result[1])

  ;do squr fitting
  initial_sigma = EXP(MEAN(ALOG(average_sigma) + 0.5 * ALOG(K_array)))
  sigma_fit2 = initial_sigma / SQRT(K_array)
  print, 'initial sigma: ', initial_sigma
   
  
  ; Set up variables for the plot. Normally, these values would be
  ; passed into the program as positional and keyword parameters.
  xtitle = 'K'
  ytitle = 'sigma (mm)'
  title = 'Localization precision vs K ('+STRTRIM(FIX(separa), 1)+')'
  position = [0.125, 0.125, 0.9, 0.925]
  thick = (!D.Name EQ 'PS') ? 3 : 1
  ; Set up a "window" for the plot. The PostScript output will have
  ; the same aspect ratio as the graphics window on the display.
  cgDisplay, 600, 1000, Title='Localization Plot'
  ; Draw the line plot.
  cgPlot, K_array[1:*], average_sigma[1:*], Color='black', PSym=16, SymColor='black', $
    SymSize=1.0, Thick=thick, Title=title, XTitle=xtitle, YTitle=ytitle, $
    Position=position, YRange=[0, 1.5], XRange=[-10,1010], YStyle=1, $
    ERR_YLow=stddev_sigma[1:*], ERR_YHigh=stddev_sigma[1:*], ERR_Color='black'
  cgPlots, K_array[1:N_ELEMENTS(K_array)-2], sigma_fit2[1:N_ELEMENTS(K_array)-2], thick=3, Color='blu5'
;  cgPlots, K_array[1:N_ELEMENTS(K_array)-2], sigma_fit1[1:N_ELEMENTS(K_array)-2], Color='red5'


  
END

;@GJ, 2023/7/8, resolution versus K subset by maximal likehood
;@GJ, 2023/7/11, standard deviation, error bar plotting and sigma fitting curve is added
;@GJ, 2023/7/15, k-means clustering, separating the twin pixels clusters
;@GJ, 2023/7/16, estimating the real location of the magnetic nanorubbles based on mirror symmetry
PRO resolution_vs_K_max_twin_clustering_mirror


  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
    fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
  ;fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF


  image_dcm = read_dicom(fn_dcm)
  ;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj
  
  ;@GJ, adding two together
  ;image_dcm += SHIFT(image_dcm, 50, 0)
  
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm, MIN=3000)
  it_fig = TVRD(TRUE=1)
  ResoK_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'ResoK_max\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN FILE_MKDIR, ResoK_dir
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig

  ;  CURSOR, X, Y, /DEVICE
    X = cols-1;0;127;
  ;X = 127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageLeft.png', it_fig


  Pixels_N = N_ELEMENTS(image_left)
  print, '# pixels: ', pixels_N

  K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000]
  ;  K_array = [100]
  ;@GJ, 2023/7/11, the rounds of calculagint the standard deviation
  N_rounds = 3;1;10;3;10;4;100.;1.;100.
  sigma_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  rep_N = 500.
  loc_x = DBLARR(N_rounds, N_ELEMENTS(K_array), rep_N) * 0.
  loc_y = DBLARR(N_rounds, N_ELEMENTS(K_array), rep_N) * 0.
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = N_ELEMENTS(K_array)-1, 0, -1 DO BEGIN
      K = K_array[k_ind]
      new_image = image_left * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(image_left[d1], maxind)
        new_image[d1[maxind]] = 1.
        ind = ARRAY_INDICES(new_image, d1[maxind])
        loc_x[i_round, k_ind, i] = ind[0]
        loc_y[i_round, k_ind, i] = ind[1]
      ENDFOR
      
      array = REFORM([loc_x[i_round,k_ind,*], loc_y[i_round,k_ind,*]])
      ; Compute cluster weights, using three clusters:
      weights = CLUST_WTS(array, N_CLUSTERS = 2)
      ; Compute the classification of each sample:
      result = CLUSTER(array, weights, N_CLUSTERS = 2)
      IF i_round EQ 0 AND k_ind EQ N_ELEMENTS(K_array)-1 THEN BEGIN
      ; Plot each cluster using a different symbol:
        IPLOT, array[*, WHERE(result eq 0)], $
          LINESTYLE = 6, SYM_INDEX = 2
        IPLOT, array[*, WHERE(result eq 1)], /OVERPLOT, $
          LINESTYLE = 6, SYM_INDEX = 4
      ENDIF
      
      ;caculate the pixels
      IF i_round EQ 0 THEN new_image[array[0, WHERE(result eq 1)], array[1, WHERE(result eq 1)]]=2
      twin1_pt = [MEAN(array[0, WHERE(result eq 0)]), MEAN(array[1, WHERE(result eq 0)])]
      twin2_pt = [MEAN(array[0, WHERE(result eq 1)]), MEAN(array[1, WHERE(result eq 1)])]
      center_pt = (twin1_pt + twin2_pt)/2.
      IF i_round EQ 0 THEN new_image[twin1_pt[0], twin1_pt[1]]=3
      IF i_round EQ 0 THEN new_image[twin2_pt[0], twin2_pt[1]]=4
      center_dist = SQRT((twin1_pt[0]-twin2_pt[0])^2 + (twin1_pt[1]-twin2_pt[1])^2)
      print, 'K: ', k_ind
      print, 'cluster center dist=', center_dist*pixelSp[0], ' mm'

      ;@GJ, 2023/7/15, move the two cluster to the center by checking the existance
      test_image = new_image * 0.
      flag = 0.
      FOR i_rep=0, rep_N-1 DO BEGIN
        IF result[i_rep] EQ 0 THEN BEGIN
          x_vec_0 = loc_x[i_round, k_ind, i_rep] - twin1_pt[0]
          y_vec_0 = loc_y[i_round, k_ind, i_rep] - twin1_pt[1]
          
          ;the twin location based on x-direction mirroring
          twin_loc_x = twin2_pt[0] + (-x_vec_0)
          twin_loc_y = twin2_pt[1] + y_vec_0
          IF twin_loc_x GE 0 AND twin_loc_x LE cols-1 AND twin_loc_y GE 0 AND twin_loc_y LT rows-1 THEN BEGIN
            IF new_image[twin_loc_x, twin_loc_y] GT 0 THEN BEGIN
              real_loc_x = center_pt[0] + x_vec_0
              real_loc_y = center_pt[1] + y_vec_0
              test_image[real_loc_x, real_loc_y] = 1.
              IF flag EQ 0 THEN real_loc = [real_loc_x, real_loc_y] ELSE real_loc = [[real_loc], [real_loc_x, real_loc_y]]
              flag++
              IF i_round EQ 0 THEN new_image[real_loc_x, real_loc_y] = 5
            ENDIF
          ENDIF  
        ENDIF
      ENDFOR
;      IF K EQ 1000 THEN iimage, new_image, title='image with merging'
      
      IF N_ELEMENTS(real_loc[0,*]) GT 1 THEN sigma_array[k_ind, i_round] = SQRT(0.5 * ABS(C_CORRELATE(REFORM(real_loc[0,*])*pixelSp[0], REFORM(real_loc[1, *])*pixelSp[1], 0, /COVARIANCE, /DOUBLE)))
      ; R, G, B color table - make everything red
      IF i_round EQ 0 THEN BEGIN
        r = REPLICATE(255b,256)
        r[0:112] = 0
        g = REPLICATE(0b,256)
        b = g
        WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_ResoK'+STRTRIM(FIX(K), 1)+'.png', BYTSCL(new_image);, r, g, b
        IF k_ind EQ N_ELEMENTS(K_array)-1 THEN iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
        ;iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
      ENDIF
    ENDFOR
    IF i_round EQ 0 THEN BEGIN
      ;      print, 'sigma_array: ', sigma_array
      iplot, K_array[1:*], sigma_array[1:*, i_round], /xlog, ytitle='sigma (mm)', xtitle='K', title='Localization precision vs K (1 round)'
      iplot, K_array[1:*], sigma_array[1:*, i_round], ytitle='sigma (mm)', xtitle='K', title='Localization precision vs K (1 round)'
    ENDIF
  ENDFOR
  

  
  ;@GJ, 2023/7/11, calculate the mean and standardev
  average_sigma = DBLARR(N_ELEMENTS(K_array))
  stddev_sigma = DBLARR(N_ELEMENTS(K_array))
  FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
    average_sigma[k_ind] = MEAN(sigma_array[k_ind, *])
    stddev_sigma[k_ind] = STDDEV(sigma_array[k_ind, *])
  ENDFOR
  iplot, ALOG(K_array), ALOG(average_sigma), ytitle='ALOG sigma (mm)', xtitle='ALOG K', title='Localization precision vs K'
  ;  iplot, K_array, stddev_sigma, ytitle='stddev sigma (mm)', xtitle='K', xrange=[0, 200], title='Localization precision vs K'

  result = LINFIT(ALOG(K_array), ALOG(average_sigma), MEASURE_ERRORS=ALOG(stddev_sigma))
  print, 'fitting result: ', result
  sigma_fit1 = EXP(result[0]) * (K_array^result[1])

  ;do sqrt fitting
  initial_sigma = EXP(MEAN(ALOG(average_sigma) + 0.5 * ALOG(K_array)))
  sigma_fit2 = initial_sigma / SQRT(K_array)
  print, 'initial sigma: ', initial_sigma


  ; Set up variables for the plot. Normally, these values would be
  ; passed into the program as positional and keyword parameters.
  xtitle = 'K'
  ytitle = 'sigma (mm)'
  title = 'Localization precision vs K'
  position = [0.125, 0.125, 0.9, 0.925]
  thick = (!D.Name EQ 'PS') ? 3 : 1
  ; Set up a "window" for the plot. The PostScript output will have
  ; the same aspect ratio as the graphics window on the display.
  cgDisplay, 600, 1000, Title='Localization Plot'
  ; Draw the line plot.
  cgPlot, K_array[1:*], average_sigma[1:*], Color='black', PSym=16, SymColor='black', $
    SymSize=1.0, Thick=thick, Title=title, XTitle=xtitle, YTitle=ytitle, $
    Position=position, YRange=[0, 1.5], XRange=[-10,1010], YStyle=1, $
    ERR_YLow=stddev_sigma[1:*], ERR_YHigh=stddev_sigma[1:*], ERR_Color='black'
  cgPlots, K_array[1:N_ELEMENTS(K_array)-2], sigma_fit2[1:N_ELEMENTS(K_array)-2], thick=3, Color='blu5'
  ;  cgPlots, K_array[1:N_ELEMENTS(K_array)-2], sigma_fit1[1:N_ELEMENTS(K_array)-2], Color='red5'
END


;@GJ, 2023/7/8, resolution versus K subset by maximal likehood
;@GJ, 2023/7/11, standard deviation, error bar plotting and sigma fitting curve is added
;@GJ, 2023/7/24, analyzing images with multiple centers
PRO resolution_vs_K_max_sigma_plot_multiple_centers


  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\duizhaotu.dcm'
    fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
;  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\dandu teshutanzhen.dcm'
;  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\4_7_2023-05-05_10-47-57.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF


  image_dcm = read_dicom(fn_dcm)
  ;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4

  ;  CURSOR, X, Y, /DEVICE
    X = cols-1;0;127;
  ;X = 127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  ResoK_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'ResoK_max_mc\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN  FILE_MKDIR, ResoK_dir
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageLeft.png', it_fig

  Pixels_N = N_ELEMENTS(image_left)
  print, '# pixels: ', pixels_N

  K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000]
  ;  K_array = [100]
  ;@GJ, 2023/7/11, the rounds of calculagint the standard deviation
  N_rounds = 5;00.
  sigma_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
      K = K_array[k_ind]
      new_image = image_left * 0.
      rep_N = 500
      loc_x = DBLARR(rep_N) * 0.
      loc_y = DBLARR(rep_N) * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(image_left[d1], maxind)
        new_image[d1[maxind]] = 1.
        ind = ARRAY_INDICES(new_image, d1[maxind])
        IF k_ind EQ N_ELEMENTS(K_array)-1 THEN BEGIN
          print, 'K=', K
          image_left[ind[0], ind[1]] = 0.
        ENDIF
        loc_x[i] = ind[0] * pixelSp[0]
        loc_y[i] = ind[1] * pixelSp[1]
      ENDFOR
      IF i_round EQ 0 THEN BEGIN
        r = REPLICATE(255b,256)
        r[0:112] = 0
        g = REPLICATE(0b,256)
        b = g
        WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_ResoK'+STRTRIM(FIX(K), 1)+'.png', BYTSCL(new_image);, r, g, b
;        iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
      ENDIF
    ENDFOR
    WINDOW, 2, XSIZE = X, YSIZE = rows
    TV, BYTSCL(image_left)
    it_fig = TVRD(TRUE=1)
    WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_'+STRTRIM(FIX(i_round), 1)+'_imageLeft01.png', it_fig
  ENDFOR
  
  
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
      K = K_array[k_ind]
      new_image = image_left * 0.
      rep_N = 500
      loc_x = DBLARR(rep_N) * 0.
      loc_y = DBLARR(rep_N) * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(image_left[d1], maxind)
        new_image[d1[maxind]] = 1.
        ind = ARRAY_INDICES(new_image, d1[maxind])
        IF k_ind EQ N_ELEMENTS(K_array)-1 THEN BEGIN
          print, 'K=', K
          image_left[ind[0], ind[1]] = 0.
        ENDIF
        loc_x[i] = ind[0] * pixelSp[0]
        loc_y[i] = ind[1] * pixelSp[1]
      ENDFOR
      IF i_round EQ 0 THEN BEGIN
        r = REPLICATE(255b,256)
        r[0:112] = 0
        g = REPLICATE(0b,256)
        b = g
        WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_ResoK'+STRTRIM(FIX(K), 1)+'.png', BYTSCL(new_image);, r, g, b
        ;        iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
      ENDIF
    ENDFOR
    WINDOW, 3, XSIZE = X, YSIZE = rows
    TV, BYTSCL(image_left)
    it_fig = TVRD(TRUE=1)
    WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_2nd_'+STRTRIM(FIX(i_round), 1)+'_imageLeft01.png', it_fig
  ENDFOR
END

;@GJ, 2023/7/24, 3D segmentation and sampling based on maximum likelihood modeling
PRO ThreeD_DICOM_vol_MLM

;  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\01_20230308_130nmiv3d_ori\01DICOM_sorted\MPI_combined_ori\'
;  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\02_20230311_130nmiv3d_ori\Dicom_sorted\MPI_combined_ori\'
;  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\03_20230404_25nmiv2h\Dicom_sorted\MPI_combined_ori\'
;  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\04_20230404_25nmSv2h\Dicom_sorted\MPI_combined_ori\'
;  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\05_AFLD\DICOM_sorted\MPI_668e_2023_20230517_031\'
;  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\09_20230603NO1\DICOM_sorted\MPI_7c62_2023_20230517_004\'
  dcm_dir = 'F:\MPI_Tianjie\NC_COMPASS\WangQianData\10_20230603NO2\DICOM_sorted\MPI_3049_2023_20230517_001\'
  image_read, dcm_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
  print, 'direction1 =', direction
  print, 'origin_loc1 =', origin_loc
  ResoK_dir = dcm_dir+'ResoK_max_mc\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN FILE_MKDIR, REsoK_dir
  
  Pixels_N = N_ELEMENTS(vol_HU_cube)
  K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000LL, 10000LL, 15000LL, 20000LL, 25000LL, 30000LL, 35000LL, 40000LL]
  ;  K_array = [100]
  ;@GJ, 2023/7/11, the rounds of calculagint the standard deviation
  
  N_rounds = 1;5;00.
  sigma_array = DBLARR(N_ELEMENTS(K_array), N_rounds) * 0.
  Nvolume = DBLARR(N_ELEMENTS(K_array), N_rounds) * 0.

  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = 0, N_ELEMENTS(K_array)-1 DO BEGIN
      K = K_array[k_ind]
      new_image_vol = vol_HU_cube * 0.
      rep_N = 20000
      loc_x = DBLARR(rep_N) * 0.
      loc_y = DBLARR(rep_N) * 0.
      loc_z = DBLARR(rep_N) * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(vol_HU_cube[d1], maxind)
        new_image_vol[d1[maxind]] = 1.
        ind = ARRAY_INDICES(new_image_vol, d1[maxind])
        loc_x[i] = ind[0] * vol_HU_cube_resolution
        loc_y[i] = ind[1] * vol_HU_cube_resolution
        loc_z[i] = ind[2] * vol_HU_cube_resolution
      ENDFOR
      IF i_round EQ 0 THEN BEGIN
        r = REPLICATE(255b,256)
        r[0:112] = 0
        g = REPLICATE(0b,256)
        b = g
        new_image_z = REFORM(vol_HU_cube[*,*,0]) * 0.
        new_image_x = REFORM(vol_HU_cube[0,*,*]) * 0.
        new_image_y = REFORM(vol_HU_cube[*,0,*]) * 0.
        FOR l = 0, N_ELEMENTS(new_image_z[*,0])-1 DO BEGIN
          FOR m = 0, N_ELEMENTS(new_image_z[0,*])-1 DO BEGIN
            new_image_z[l, m] = MAX(new_image_vol[l, m, *])
            new_image_x[l, m] = MAX(new_image_vol[*, l, m])
            new_image_y[l, m] = MAX(new_image_vol[l, *, m])
          ENDFOR
        ENDFOR
        ind = where(new_image_vol GT 0, Counts)
        Nvolume[k_ind, i_round] = Counts
        WRITE_PNG, ResoK_dir+STRTRIM(FIX(rep_N), 1)+'_MLM_x_'+STRTRIM(FIX(K, type=3), 1)+'.png', BYTSCL(new_image_x);, r, g, b
        WRITE_PNG, ResoK_dir+STRTRIM(FIX(rep_N), 1)+'_MLM_y_'+STRTRIM(FIX(K, type=3), 1)+'.png', BYTSCL(new_image_y);, r, g, b
        WRITE_PNG, ResoK_dir+STRTRIM(FIX(rep_N), 1)+'_MLM_z_'+STRTRIM(FIX(K, type=3), 1)+'.png', BYTSCL(new_image_z);, r, g, b
        ;        iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K, type=3), 1)+')'
      ENDIF
    ENDFOR
  ENDFOR
  
  iplot, K_array, Nvolume, title='# pixels vs K'
;  iplot, K_array, Nvolume, /xlog, title='# pixels vs K'
  ;ivolume, new_image_vol
END



;@GJ, 2023/7/31, synthesis magnetization based on twin images
PRO magnetization_twin_synthesis
  
  ;find the image and center
  twin_image_center_line, image_left, pixelSp, twin1_pt, twin2_pt, center_pt
  
  ;get the line profile
  linefit_result = LINFIT([twin1_pt[0], twin2_pt[0]], [twin1_pt[1], twin2_pt[1]])
  
  x = FINDGEN(N_ELEMENTS(image_left[*,0]))
  y = linefit_result[0] + x * linefit_result[1]
  
  temp_image = image_left
  temp_image[x, y] = MAX(image_left)
  iimage, image_left + temp_image
  
  line_ind = SQRT((x * pixelSp[0])^2 + (y * pixelSp[1])^2)
  line_intensity = image_left[x, y]
 
  line_signal_pos = line_intensity
  line_signal_pos[0:center_pt[0]] = 0.
  ;@GJ, 2023/8/1, normalization based on FFP velocity
  line_signal_pos *= SIN(!PI/MAX(x) * x)
  line_signal_integral_pos = TOTAL(line_signal_pos, /cumulative)
  iplot, line_ind, line_signal_integral_pos, title='M-H curve'
  
  line_signal_neg = line_intensity
  line_signal_neg[center_pt[0]:*] = 0.
  ;@GJ, 2023/8/1, normalization based on FFP velocity
  line_signal_neg *= SIN(!PI/MAX(x) * x+!PI)
  line_signal_integral_neg = REVERSE(TOTAL(REVERSE(line_signal_neg), /cumulative) + TOTAL(line_signal_pos))
  iplot, line_ind, line_signal_integral_neg, /overplot
  
  M_range = MAX(line_signal_integral_pos)-MIN(line_signal_integral_pos)
  iplot, line_ind, line_signal_integral_pos-0.5*m_range, title='M-H curve'
  iplot, line_ind, line_signal_integral_neg-0.5*m_range, /overplot
  
  
  ;GJ, 2022/1/30, calculate Debye model
  Nelem = N_ELEMENTS(x)
  H_0 = 50.
  H_array = x / MAX(x) * H_0 - H_0/2.; mT
  particle_size = 27.;19.;40.;50.;19.;40.;50;70.;19.;30.
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000.; kA/m
  T_p = 20.;36.5; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  beta_particle = Msat_T*(particle_size^3) /24./1.380469/T_p_kelvin ; in 1/mT
  ;  beta_particle = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
  FWHM_particle = 4.16/beta_particle/3.5; mT/mm
  M = Msat_kAm*(1./TANH(beta_particle*H_array) - 1./(beta_particle*H_array))
  M[WHERE(~FINITE(M))] = 0.
  signal_langevin_temp = M[1:*]-M[0:Nelem-1];/delta_t; kA/m/us
  signal_meas = line_signal_pos[0:Nelem-2]
  signal_langevin = signal_langevin_temp/MAX(ABS(signal_langevin_temp))*MAX(ABS(signal_meas))
  
  delta_t = 1.0; us
  time = x * delta_t
  time_half = time[0:Nelem-2]
  window, 0
  plot, time_half, signal_langevin, title='Langevin & Signal', XTITLE='time [us]', YTITLE='signal [a.u.]'
  oplot, time_half, signal_meas, line=2
  oplot, time_half, ABS(H_array)

  ;GJ, 2022/1/30, calculate Debye model
  relaxation_time = FINDGEN(100)+1.
  correlation_array = DBLARR(100)*0.
  MSE_array = DBLARR(100)*0.+10000.
  FOR i = 0, 99 DO BEGIN
    signal_nat = non_adiabatic(time_half, signal_langevin, relaxation_time[i])
    signal_nat = signal_nat/MAX(signal_nat)*MAX(signal_meas)
    correlation_array[i] = CORRELATE(signal_nat, signal_meas)
    MSE_array[i] = calMSE(signal_nat, signal_meas)
  ENDFOR

  ;max_corr = MAX(correlation_array, max_ind)
  min_MSE = MIN(MSE_array, max_ind)
  print, 'relaxation time = ', relaxation_time[max_ind]
  signal_nat = non_adiabatic(time_half, signal_langevin, relaxation_time[max_ind])
  signal_nat = signal_nat/MAX(signal_nat)*MAX(signal_meas)
  ;  oplot, time, signal_nat, line=3
  corre_Mono = CORRELATE(signal_meas, signal_nat)
  PRINT, 'correlation coeff (Mono): ', corre_Mono

  iplot, time_half, signal_langevin, LINESTYLE = 2, thick = 2, title='Debye Relaxation', xtitle='time [us]', ytitle='signal [a.u.]', BACKGROUND = 'FFFFFF'x, COLOR = 0
  iplot, time_half, signal_meas, thick = 2, LINESTYLE = 0, color='blue',/overplot
  iplot, time_half, signal_nat, thick = 2, color='red',/overplot

END


PRO random_sampling_single_center, image_left, shrink_image, loc_x, loc_y

  Pixels_N = N_ELEMENTS(image_left)
  ;  print, '# pixels: ', pixels_N

  ;K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000]
  K_array = [2000]
  ;@GJ, 2023/7/11, the rounds of calculagint the standard deviation
  N_rounds = 1;3;1;10;3;10;4;100.;1.;100.
  sigma_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  rep_N = 1000.
  loc_x = DBLARR(N_rounds, N_ELEMENTS(K_array), rep_N) * 0.
  loc_y = DBLARR(N_rounds, N_ELEMENTS(K_array), rep_N) * 0.
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = N_ELEMENTS(K_array)-1, 0, -1 DO BEGIN
      K = K_array[k_ind]
      shrink_image = image_left * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(image_left[d1], maxind)
        shrink_image[d1[maxind]] = 1.
        ind = ARRAY_INDICES(shrink_image, d1[maxind])
        loc_x[i_round, k_ind, i] = ind[0]
        loc_y[i_round, k_ind, i] = ind[1]
      ENDFOR
    ENDFOR
  ENDFOR

  sigma = SQRT(0.5 * ABS(C_CORRELATE(loc_x, loc_y, 0, /COVARIANCE, /DOUBLE)))
  ;  print, 'sigma', sigma
END


;@GJ, 2023/7/31, calculating center and alinement
PRO twin_image_center_line, image_left, pixelSp, twin1_pt, twin2_pt, center_pt, mask

  ;loading DICOM image and get informations
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\4_7_2023-05-05_10-47-57.dcm'
  ;  fn_dcm = 'F:\MPI_Tianjie\NC_COMPASS\LineProfile\019_20230603yuanyeNO13\projections\combined_0000_tangential_image_0000.dcm'
  ;  fn_dcm = 'C:\D_drive\MPI_Tianjie\NanoRod_MPI\DICOM_images\duizhaotu.dcm'
  ;fn_dcm = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF


  image_dcm = read_dicom(fn_dcm)
  ;  iimage, image_dcm, title='original'
  obj = OBJ_NEW('IDLffDICOM', fn_dcm)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  ; Get the row & column size of the image(s):
  temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm[0,*])
  temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
  IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm[*,0])

  OBJ_DESTROY, obj

  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = cols, YSIZE = rows
  TV, BYTSCL(image_dcm, MIN=3000)
  it_fig = TVRD(TRUE=1)
  ResoK_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'ResoK_max\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN FILE_MKDIR, ResoK_dir
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_image.png', it_fig

  ;  CURSOR, X, Y, /DEVICE
  X = cols-1;0;127;
  ;X = 127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  ;@GJ, 2023/8/3, adding mask
  IF N_ELEMENTS(mask) THEN image_left *= mask
  image_right = DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  WINDOW, 1, XSIZE = X, YSIZE = rows
  TV, BYTSCL(image_left)
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_imageLeft.png', it_fig


  Pixels_N = N_ELEMENTS(image_left)
  print, '# pixels: ', pixels_N

  ;K_array = [1, 5, 10, 30, 50, 80, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 2000]
  K_array = [1000]
  ;@GJ, 2023/7/11, the rounds of calculagint the standard deviation
  N_rounds = 1;3;1;10;3;10;4;100.;1.;100.
  sigma_array = DBLARR(N_ELEMENTS(K_array), N_rounds)
  rep_N = 500.
  loc_x = DBLARR(N_rounds, N_ELEMENTS(K_array), rep_N) * 0.
  loc_y = DBLARR(N_rounds, N_ELEMENTS(K_array), rep_N) * 0.
  FOR i_round = 0, N_rounds-1 DO BEGIN
    FOR k_ind = N_ELEMENTS(K_array)-1, 0, -1 DO BEGIN
      K = K_array[k_ind]
      new_image = image_left * 0.
      FOR i=0, rep_N-1 DO BEGIN
        seed = !NULL
        d1 = ULONG64(DOUBLE(Pixels_N) * RANDOMU(Seed, K, /DOUBLE))
        ;      x_ave = MEAN(pixel_array[d1, 0])
        ;      y_ave = MEAN(pixel_array[d1, 1])
        Max_value = MAX(image_left[d1], maxind)
        new_image[d1[maxind]] = 1.
        ind = ARRAY_INDICES(new_image, d1[maxind])
        loc_x[i_round, k_ind, i] = ind[0]
        loc_y[i_round, k_ind, i] = ind[1]
      ENDFOR

      array = REFORM([loc_x[i_round,k_ind,*], loc_y[i_round,k_ind,*]])
      ; Compute cluster weights, using three clusters:
      weights = CLUST_WTS(array, N_CLUSTERS = 2)
      ; Compute the classification of each sample:
      result = CLUSTER(array, weights, N_CLUSTERS = 2)
      IF i_round EQ 0 AND k_ind EQ N_ELEMENTS(K_array)-1 THEN BEGIN
        ; Plot each cluster using a different symbol:
        IPLOT, array[*, WHERE(result eq 0)], $
          LINESTYLE = 6, SYM_INDEX = 2
        IPLOT, array[*, WHERE(result eq 1)], /OVERPLOT, $
          LINESTYLE = 6, SYM_INDEX = 4
      ENDIF

      ;caculate the pixels
      IF i_round EQ 0 THEN new_image[array[0, WHERE(result eq 1)], array[1, WHERE(result eq 1)]]=2
      twin1_pt = [MEAN(array[0, WHERE(result eq 0)]), MEAN(array[1, WHERE(result eq 0)])]
      twin2_pt = [MEAN(array[0, WHERE(result eq 1)]), MEAN(array[1, WHERE(result eq 1)])]
      center_pt = (twin1_pt + twin2_pt)/2.
      IF i_round EQ 0 THEN new_image[twin1_pt[0], twin1_pt[1]]=3
      IF i_round EQ 0 THEN new_image[twin2_pt[0], twin2_pt[1]]=4
      center_dist = SQRT((twin1_pt[0]-twin2_pt[0])^2 + (twin1_pt[1]-twin2_pt[1])^2)
      print, 'K: ', k_ind
      print, 'cluster center dist=', center_dist*pixelSp[0], ' mm'
      print, 'center point: ', center_pt

      ;@GJ, 2023/7/15, move the two cluster to the center by checking the existance
      test_image = new_image * 0.
      flag = 0.
      FOR i_rep=0, rep_N-1 DO BEGIN
        IF result[i_rep] EQ 0 THEN BEGIN
          x_vec_0 = loc_x[i_round, k_ind, i_rep] - twin1_pt[0]
          y_vec_0 = loc_y[i_round, k_ind, i_rep] - twin1_pt[1]

          ;the twin location based on x-direction mirroring
          twin_loc_x = twin2_pt[0] + (-x_vec_0)
          twin_loc_y = twin2_pt[1] + y_vec_0
          IF twin_loc_x GE 0 AND twin_loc_x LE cols-1 AND twin_loc_y GE 0 AND twin_loc_y LT rows-1 THEN BEGIN
            IF new_image[twin_loc_x, twin_loc_y] GT 0 THEN BEGIN
              real_loc_x = center_pt[0] + x_vec_0
              real_loc_y = center_pt[1] + y_vec_0
              test_image[real_loc_x, real_loc_y] = 1.
              IF flag EQ 0 THEN real_loc = [real_loc_x, real_loc_y] ELSE real_loc = [[real_loc], [real_loc_x, real_loc_y]]
              flag++
              IF i_round EQ 0 THEN new_image[real_loc_x, real_loc_y] = 5
            ENDIF
          ENDIF
        ENDIF
      ENDFOR
      ;      IF K EQ 1000 THEN iimage, new_image, title='image with merging'

      ; R, G, B color table - make everything red
      IF i_round EQ 0 THEN BEGIN
        r = REPLICATE(255b,256)
        r[0:112] = 0
        g = REPLICATE(0b,256)
        b = g
        WRITE_PNG, ResoK_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_ResoK'+STRTRIM(FIX(K), 1)+'.png', BYTSCL(new_image);, r, g, b
        ;IF k_ind EQ N_ELEMENTS(K_array)-1 THEN iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
        iimage, new_image, title='resolution image (K = '+STRTRIM(FIX(K), 1)+')'
      ENDIF
    ENDFOR
  ENDFOR

END

;@GJ, 2023/8/22, do pattern analysis based on donut shape
PRO twin_donut_pattern_analysis
  ;  ;find the image and center
  ;  twin_image_center_line, image_left, pixelSp_0, twin1_pt, twin2_pt, center_pt
  ;
  ;
  ;  ;get the twin distance
  ;  twin_distance_0 = ABS(ROUND(twin1_pt[0] - twin2_pt[0]))
  ;  IF twin_distance_0 LT 30 THEN twin_distance_0 = 46

  image_dcm_fn = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
  IF N_ELEMENTS(image_dcm_fn) EQ 0 THEN image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  IF STRLEN(image_dcm_fn) EQ 0 THEN image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  ;   image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\002_dandu tanzhen\dandu tanzhen.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\004_duizhaotu\duizhaotu.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\005_0626YHHprobe_2023-06-26_15-22-01\combined_0000_combined_image_0000.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\006_0626YHHprobe_2023-06-26_15-32-27\combined_0000_combined_image_0000.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\007_0626YHHprobe_2023-06-26_15-42-33\combined_0000_combined_image_0000.dcm'
  image_left = DOUBLE(read_dicom(image_dcm_fn))
  IF STRPOS(image_dcm_fn, '2023-06-26_15-42-33') NE -1 THEN image_left = ROTATE(image_left, 1)
  IF STRPOS(image_dcm_fn, '001_3teshu') NE -1 THEN image_left = ROTATE(image_left, 2)
  obj = OBJ_NEW('IDLffDICOM', image_dcm_fn)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp_1=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  OBJ_DESTROY, obj
  
  ;get the image size
  image_size = SIZE(image_left, /dim)
  XSIZE = image_size[0]
  YSIZE = image_size[1]
  
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 23
  WINDOW, 0, XSIZE = XSIZE, YSIZE = YSIZE
  TV, BYTSCL(image_left)
  it_fig = TVRD(TRUE=1)
  Donut_dir = FILE_DIRNAME(image_dcm_fn, /MARK_DIRECTORY)+'Donut_Shape\'
  IF FILE_TEST(Donut_dir, /directory) LT 1 THEN FILE_MKDIR, Donut_dir
  WRITE_PNG, Donut_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_image_left.png', it_fig
  
  ;@GJ, 2023/8/11, surface plot
  ; Load color tables.
  window, 1
  cgLoadCT, 33, RGB_Table=rainbowPalette
  ;  cgSurf, image_left[0:195,*], ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
  ;  window, 4
  ;  cgSurf, image_left[196:*,*], ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
  ;  window, 5
  cgSurf, image_left, ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, Donut_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'PSF_surf_original.png', it_fig
  
  ;image modification
  image_shift = SHIFT(image_left, 20., 0)
  new_image = image_left + image_shift
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 23
  WINDOW, 2, XSIZE = XSIZE, YSIZE = YSIZE
  TV, BYTSCL(new_image)
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, Donut_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_image_2pixel_shift.png', it_fig
  
  iimage, image_left, title='original'
  iimage, new_image, title='2 pixel shift'
  
  iimage, new_image - image_left, title='difference'
  

END





;@GJ, 2023/8/16, compare the lion eye and evil eye (Eye of Sauron) difference
PRO twin_ISTED_comparison

  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
  twin_ISTED_calculation, image_dcm_fn, line_cor1, line_int_x_direc1, line_int_y_direc1, ISTED_map_PSF1

  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\005_0626YHHprobe_2023-06-26_15-22-01\combined_0000_combined_image_0000.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  twin_ISTED_calculation, image_dcm_fn, line_cor2, line_int_x_direc2, line_int_y_direc2, ISTED_map_PSF2
  
  iplot, line_cor1, line_int_x_direc1-MIN(line_int_x_direc1), thick=2, color='blue', xtitle='distance [mm]', ytitle='Intensity', title='PSF (z direction)'
  iplot, line_cor2, line_int_x_direc2-MIN(line_int_x_direc2), thick=2, color='red', /overplot
  
  iplot, line_cor1, line_int_y_direc1-MIN(line_int_y_direc1), thick=2, color='blue', xtitle='distance [mm]', ytitle='Intensity', title='PSF (x direction)'
  iplot, line_cor2, line_int_y_direc2-MIN(line_int_y_direc2), thick=2, color='red', /overplot
  
END

;@GJ, 2023/8/16, intensity calculation based on PNAS2000, fluorescence microscopy
PRO twin_ISTED_calculation, image_dcm_fn, line_cor, line_int_x_direc, line_int_y_direc, ISTED_map_PSF
  
  IF N_ELEMENTS(image_dcm_fn) EQ 0 THEN BEGIN
    image_dcm_fn = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
    IF N_ELEMENTS(image_dcm_fn) EQ 0 THEN image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
    IF STRLEN(image_dcm_fn) EQ 0 THEN image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  ENDIF
  
  ;   image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\002_dandu tanzhen\dandu tanzhen.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\004_duizhaotu\duizhaotu.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\005_0626YHHprobe_2023-06-26_15-22-01\combined_0000_combined_image_0000.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\006_0626YHHprobe_2023-06-26_15-32-27\combined_0000_combined_image_0000.dcm'
  ;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\007_0626YHHprobe_2023-06-26_15-42-33\combined_0000_combined_image_0000.dcm'
  image_left = DOUBLE(read_dicom(image_dcm_fn))
  IF STRPOS(image_dcm_fn, '2023-06-26_15-42-33') NE -1 THEN image_left = ROTATE(image_left, 1)
  IF STRPOS(image_dcm_fn, '001_3teshu') NE -1 THEN image_left = ROTATE(image_left, 2)
  iimage, image_left, title='original image'
  obj = OBJ_NEW('IDLffDICOM', image_dcm_fn)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp_1=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  OBJ_DESTROY, obj

  ;get the image size
  image_size = SIZE(image_left, /dim)
  XSIZE = image_size[0]
  YSIZE = image_size[1]
  
  twin_distance = 46.;twin_distance_0 * pixelSp_0[0] / pixelSp_1[0]
  ;win_distance = 4.48 / pixelSp_1[0]
  print, 'twin_distance: ', twin_distance

  ;get the image background noise
  image_left_threshold = 0.3 * MAX(image_left)
  
  ;do pixel searching
  ISTED_map = DBLARR(XSIZE, YSIZE) * 0.
  FOR i=twin_distance, XSIZE-twin_distance-1 DO BEGIN
    FOR j=twin_distance/2., YSIZE-twin_distance/2.-1 DO BEGIN
      IF image_left[i, j] GT image_left_threshold*0.1 THEN BEGIN
        I_integral = 0.
        FOR k=0, twin_distance/2. DO BEGIN
          r = k * pixelSp_1[0]
          dr = pixelSp_1[0]
          dS = (image_left[i-k, j] + image_left[i+k, j]) * r * dr / 2.
          I_integral += ds
        ENDFOR
        ISTED_map[i, j] = image_left[i, j] / I_integral
      ENDIF
    ENDFOR
  ENDFOR
  
  iimage, ISTED_map, title='I_STED map'
  
  maxISTED = MAX(ISTED_map, maxind)

  tao = -maxISTED / ALOG(0.001)
  new_ISTED_map = exp(-ISTED_map/tao)
  
  ;@GJ, 2023/8/12, close gaps in the mask
  ;strucElem = REPLICATE(1, 3, 3)
  ;temp_ISTED_map = DILATE(ERODE((ISTED_map EQ 0), strucElem), strucElem)
  new_ISTED_map[WHERE(ISTED_map EQ 0, /NULL)] = 0      
  iimage, new_ISTED_map, title='N map'
  
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 4
  WINDOW, 0, XSIZE = XSIZE, YSIZE = YSIZE
  TV, BYTSCL(new_ISTED_map)
  it_fig = TVRD(TRUE=1)
  ISTED_dir = FILE_DIRNAME(image_dcm_fn, /MARK_DIRECTORY)+'ISTED\'
  IF FILE_TEST(ISTED_dir, /directory) LT 1 THEN FILE_MKDIR, ISTED_dir
  WRITE_PNG, ISTED_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_ISTED_map.png', it_fig
  CURSOR, X, Y, /DEVICE
  line_cor = FINDGEN(twin_distance*3.)*pixelSp_1[0]
  n_line_cor = N_ELEMENTS(line_cor)
  ext_ISTED_map = DBLARR(XSIZE+2*n_line_cor, YSIZE+2*n_line_cor)*0.
  ext_ISTED_map[n_line_cor:XSIZE+n_line_cor-1, n_line_cor:YSIZE+n_line_cor-1] = new_ISTED_map
  line_int_x_direc = ext_ISTED_map[X+n_line_cor/2.:X+n_line_cor*3./2.-1, Y+n_line_cor]
  line_int_y_direc = ext_ISTED_map[X+n_line_cor, Y+n_line_cor/2.:Y+n_line_cor*3./2.-1]
  
  print, 'image x size', XSIZE*pixelSp_1[0], ' mm'
  ;calculate the FWHM
  maxX = MAX(line_int_x_direc, maxXind)
  minX = MIN(line_int_x_direc, minXind)
  half_x_1 = INTERPOL(line_cor[0:maxXind], line_int_x_direc[0:maxXind], minX+0.5*(maxX-minX))
  half_x_2 = INTERPOL(line_cor[maxXind:*], line_int_x_direc[maxXind:*], minX+0.5*(maxX-minX))
  FWHM_x = half_x_2 - half_x_1
  print, 'FWHM (x): ', FWHM_x, ' mm'
  
  ;calculate the FWHM
  maxY = MAX(line_int_y_direc, maxYind)
  minY = MIN(line_int_y_direc, minYind)
  half_y_1 = INTERPOL(line_cor[0:maxYind], line_int_y_direc[0:maxYind], minY+0.5*(maxY-minY))
  half_y_2 = INTERPOL(line_cor[maxYind:*], line_int_y_direc[maxyind:*], minY+0.5*(maxY-minY))
  FWHM_y = half_y_2 - half_y_1
  print, 'FWHM (y): ', FWHM_y, ' mm'
 
  ISTED_map_PSF = ext_ISTED_map[X+n_line_cor/2.:X+n_line_cor*3./2.-1, Y+n_line_cor/2.:Y+n_line_cor*3./2.-1]
  window, 1, XSIZE = n_line_cor, YSIZE = n_line_cor
  LOADCT, 33
  TV, BYTSCL(ISTED_map_PSF)
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ISTED_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_ISTED_PSF.png', it_fig
  window, 2
  cgLoadCT, 33, RGB_Table=rainbowPalette
  cgSurf, ISTED_map_PSF, ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
  it_fig = TVRD(TRUE=1)
  WRITE_PNG, ISTED_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_PSF_surf.png', it_fig
END

;@GJ, 2023/8/3, do pixel searching based on the twin imaging to find the position
;@GJ, 2023/8/12, fill the gaps and holes in the mask when searching the minimum
PRO twin_pattern_hit_MLM_localization
;  ;find the image and center
;  twin_image_center_line, image_left, pixelSp_0, twin1_pt, twin2_pt, center_pt
;  
;  
;  ;get the twin distance
;  twin_distance_0 = ABS(ROUND(twin1_pt[0] - twin2_pt[0]))
;  IF twin_distance_0 LT 30 THEN twin_distance_0 = 46
  
  image_dcm_fn = DIALOG_PICKFILE(FILTER = ['*.dcm'], title='Please select your image', PATH='F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\')
  IF N_ELEMENTS(image_dcm_fn) EQ 0 THEN image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
  IF STRLEN(image_dcm_fn) EQ 0 THEN image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
;   image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\002_dandu tanzhen\dandu tanzhen.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\004_duizhaotu\duizhaotu.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\001_3teshu\4_7_2023-05-05_10-47-57.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\005_0626YHHprobe_2023-06-26_15-22-01\combined_0000_combined_image_0000.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\006_0626YHHprobe_2023-06-26_15-32-27\combined_0000_combined_image_0000.dcm'
;  image_dcm_fn = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\007_0626YHHprobe_2023-06-26_15-42-33\combined_0000_combined_image_0000.dcm'
  image_left = DOUBLE(read_dicom(image_dcm_fn))
  IF STRPOS(image_dcm_fn, '2023-06-26_15-42-33') NE -1 THEN image_left = ROTATE(image_left, 1)
  IF STRPOS(image_dcm_fn, '001_3teshu') NE -1 THEN image_left = ROTATE(image_left, 2)
  iimage, image_left, title='original image'
  obj = OBJ_NEW('IDLffDICOM', image_dcm_fn)
  pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
  pixelSp_1=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
  OBJ_DESTROY, obj
  
  ;@GJ, 2023/8/11, surface plot
  ; Load color tables.
  window, 3
  cgLoadCT, 33, RGB_Table=rainbowPalette
;  cgSurf, image_left[0:195,*], ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
;  window, 4
;  cgSurf, image_left[196:*,*], ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
;  window, 5
  cgSurf, image_left, ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
  it_fig = TVRD(TRUE=1)
  ResoK_dir = FILE_DIRNAME(image_dcm_fn, /MARK_DIRECTORY)+'TwinPattern_ResoK\'
  IF FILE_TEST(ResoK_dir, /directory) LT 1 THEN FILE_MKDIR, ResoK_dir
  WRITE_PNG, ResoK_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'PSF_surf_original.png', it_fig
  

  twin_distance = 46.;twin_distance_0 * pixelSp_0[0] / pixelSp_1[0]
  win_distance = 4.48 / pixelSp_1[0]
  print, 'twin_distance: ', twin_distance

  ;get the image background noise
  image_left_threshold = 0.3 * MAX(image_left)
  
  
  ;get the image size
  image_size = SIZE(image_left, /dim)
  XSIZE = image_size[0]
  YSIZE = image_size[1]

  ;do pixel searching
  hit_map = DBLARR(XSIZE, YSIZE) * 0.
  FOR i=twin_distance, XSIZE-twin_distance-1 DO BEGIN
    FOR j=twin_distance/2., YSIZE-twin_distance/2.-1 DO BEGIN
      IF image_left[i, j] GT image_left_threshold*0.1 THEN BEGIN
        twin1 = image_left[i-twin_distance/2., j]
        twin2 = image_left[i+twin_distance/2., j]
        gradient_minus_m = twin1-image_left[i-twin_distance, j]
        gradient_minus = twin1 - image_left[i, j]
        gradient_plus = twin2 - image_left[i, j]
        gradient_plus_p = twin2 - image_left[i+twin_distance, j]
        gradient_y_m = 1.;image_left[i, j] - image_left[i, j-twin_distance/2.]
        gradient_y_p = 1.;image_left[i, j] - image_left[i, j+twin_distance/2.]
        IF gradient_minus_m GT 0 AND gradient_minus GT 0 AND gradient_plus GT 0 AND gradient_plus_p GT 0 AND gradient_y_m*gradient_y_p GT 0 THEN hit_map[i, j] = ABS(twin1 - twin2) / MAX([twin1, twin2])
      ENDIF
    ENDFOR
  ENDFOR
  
  iimage, hit_map, title='hit map'
  
  IF MAX(hit_map) EQ 0 THEN BEGIN
    info=DIALOG_MESSAGE('The image did not have nanorubbles!', /INFORMATION)
    RETURN
  ENDIF
  
  ;@GJ, 2023/8/12, plot the image based on the STED method
;  tao = -MAX(image_left) / ALOG(0.001)
;  image_STED = exp(-image_w_mask/tao)
;  iimage, image_STED, title='image STED'
  ;@GJ, 2023/8/12, dilate and process the image
  mask_hit = hit_map GT 0.00001
  ;@GJ, 2023/8/12, close gaps in the mask
  strucElem = REPLICATE(1, 3, 3)
  mask_hit = ERODE(DILATE(TEMPORARY(mask_hit), strucElem), strucElem)
  
  ;separate the images
  labelImg = LABEL_REGION(mask_hit, ALL_NEIGHBORS=allNeighbors, ULONG=ulong)
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
  sort_re=REVERSE(SORT(h))
  final_image = image_left * 0.
  final_mask = image_left * 0.
  twin_image = image_left * 0.
  sigma_array = DBLARR(N_ELEMENTS(h)-1) * 0.
  Npixels_array = DBLARR(N_ELEMENTS(h)-1) * 0.
  FOR i=1, N_ELEMENTS(h)-1 DO BEGIN
    center_part = FLOAT(labelImg EQ sort_re[i])
    
    ;@GJ, 2023/8/12, jump if there is no central minimum
    center_ind = WHERE(labelImg EQ sort_re[i], center_count)
    center_ind2d = ARRAY_INDICES(labelImg, center_ind)
    center_middle = [MEAN(center_ind2d[0,*]),MEAN(center_ind2d[1,*])]
    line_profile = image_left[MIN(center_ind2d[0,*]):MAX(center_ind2d[0,*]), center_middle[1]]
    line_min = MIN(line_profile, minind)
    ;if there is no minimum in the middle, continue to next iteration
    IF minind EQ 0 OR minind EQ N_ELEMENTS(line_profile)-1 THEN continue
    
    ;@GJ, 2023/8/12, removing the middle hole, fill the hole
    anti_center_part = FLOAT(labelImg NE sort_re[i])
    anti_labelImg = LABEL_REGION(anti_center_part, ALL_NEIGHBORS=allNeighbors, ULONG=ulong)
    anti_h = HISTOGRAM(anti_labelImg, REVERSE_INDICES=anti_r)
    IF N_ELEMENTS(anti_h) GE 1 THEN BEGIN
      anti_sort_re=REVERSE(SORT(anti_h))
      FOR anti_i=1, N_ELEMENTS(anti_h)-1 DO BEGIN
        center_part_add = FLOAT(anti_labelImg EQ anti_sort_re[anti_i])
        center_part += center_part_add
        center_part = FLOAT(center_part GE 1)
      ENDFOR
    ENDIF
    
    left_part_mask = SHIFT(center_part, -twin_distance/2.)
    random_sampling_single_center, image_left*left_part_mask, shrink_image1, loc_x1, loc_y1
    right_part_mask = SHIFT(center_part, twin_distance/2.)
    random_sampling_single_center, image_left*right_part_mask, shrink_image2, loc_x2, loc_y2
    
    ;calculate mirror points
    twin1_pt = [MEAN(loc_x1), MEAN(loc_y1)]
    twin2_pt = [MEAN(loc_x2), MEAN(loc_y2)]
    center_pt = (twin1_pt + twin2_pt)/2.
    center_dist = SQRT((twin1_pt[0]-twin2_pt[0])^2 + (twin1_pt[1]-twin2_pt[1])^2)
;    print, 'i: ', i
;    print, 'cluster center dist=', center_dist*pixelSp_1[0], ' mm'
;    print, 'center point: ', center_pt

    ;@GJ, 2023/7/15, move the two cluster to the center by checking the existance
    flag = 0.
    rep_N = N_ELEMENTS(loc_x1)
    FOR i_rep=0, rep_N-1 DO BEGIN
        x_vec_1 = loc_x1[0,0,i_rep] - twin1_pt[0]
        y_vec_1 = loc_y1[0,0,i_rep] - twin1_pt[1]
        intensity_1 = image_left[loc_x1[0,0,i_rep], loc_y1[0,0,i_rep]]
        
        ;the twin location based on x-direction mirroring
        twin_loc_x2 = twin2_pt[0] + (-x_vec_1)
        twin_loc_y2 = twin2_pt[1] + y_vec_1
        IF twin_loc_x2 GE 0 AND twin_loc_x2 LE XSIZE-1 AND twin_loc_y2 GE 0 AND twin_loc_y2 LT YSIZE-1 THEN BEGIN
          intensity_2 = image_left[twin_loc_x2, twin_loc_y2]
          IF shrink_image2[twin_loc_x2, twin_loc_y2] GT 0 THEN BEGIN
            real_loc_x = center_pt[0] + x_vec_1
            real_loc_y = center_pt[1] + y_vec_1
            intensity = (intensity_1 + intensity_2)/2.
            IF flag EQ 0 THEN real_loc_int = [real_loc_x, real_loc_y, intensity] ELSE real_loc_int = [[real_loc_int], [real_loc_x, real_loc_y, intensity]]
            flag++
;            print, 'flag =', flag
          ENDIF
        ENDIF
      ENDFOR
      Npixels_array[i-1] = N_ELEMENTS(real_loc_int[0,*])
      IF  Npixels_array[i-1] GT rep_N*0.5 AND MEAN(real_loc_int[2,*]) GT image_left_threshold THEN BEGIN
        sigma_array[i-1] = SQRT(0.5 * ABS(C_CORRELATE(real_loc_int[0,*]*pixelSp_1[0], real_loc_int[1,*]*pixelSp_1[1], 0, /COVARIANCE, /DOUBLE)))
        twin_image += (shrink_image1 + shrink_image2) * i
        twin_image[real_loc_int[0,*], real_loc_int[1,*]] = i+5.
        final_image[real_loc_int[0,*], real_loc_int[1,*]] = real_loc_int[2,*]
        final_mask[real_loc_int[0,*], real_loc_int[1,*]] = 1
        x_center = STRTRIM(STRING(ROUND(MEAN(real_loc_int[0,*]))),1)
        y_center = STRTRIM(STRING(ROUND(MEAN(real_loc_int[1,*]))),1)
        print, '[', x_center,',',y_center,'] sigma: ', sigma_array[i-1], ' mm'
      ENDIF
    ENDFOR
    
    iplot, Npixels_array, title = 'N pixels'
    iplot, sigma_array, title = 'sigma'

    
    ; R, G, B color table - make everything red
    r = REPLICATE(255b,256)
    r[0:112] = 0
    g = REPLICATE(0b,256)
    b = g
    
    WRITE_PNG, ResoK_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_01ori.png', BYTSCL(image_left);, r, g, b
    WRITE_PNG, ResoK_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_02hitmap.png', BYTSCL(hit_map);, r, g, b
    WRITE_PNG, ResoK_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_03twin.png', BYTSCL(twin_image);, r, g, b
    iimage, twin_image, title='twin image'
    WRITE_PNG, ResoK_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_04enhanced.png', BYTSCL(final_image);, r, g, b
    iimage, final_image, title='enhanced image'
    
    IF MAX(final_image) EQ 0 THEN BEGIN
      info=DIALOG_MESSAGE('The image did not have nanorubbles!', /INFORMATION)
      RETURN
    ENDIF
    
    window, 6
    cgLoadCT, 33, RGB_Table=rainbowPalette
    cgSurf, final_image+(1.-final_mask)* MIN(final_image), ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
    it_fig = TVRD(TRUE=1)
    WRITE_PNG, ResoK_dir+FILE_BASENAME(image_dcm_fn, '.dcm')+'_PSF_surf.png', it_fig
END

;@GJ, 2023/8/13, twin imaging based on 2nd harmonics
PRO twin_2nd_harmonics
  A = 20.; mT
  A_DC = 6.;0.;-6.; mT
  particle_size = 30.; nm
  f = 1. ;kHz
  tau = 1.;20.; ms
  plotYes = 0;1
  curveColor = 'blue'
  sine_excitation_harmonics, A, A_DC, particle_size,f,tau,H_2, M_2,frequency, u, plotYes, curveColor
;  iplot, frequency/f, u, XRANGE=[0, 9],background = 'ffffff'x, color = '000000'x, xtitle='n of f'
;  tau = 29.;ms
;  curveColor = 'red'
;  sine_excitation_harmonics, A, A_DC, particle_size,f,tau,H_2, M_2,frequency, u, plotYes, curveColor
  
  FOV = 200;mm
  gradient = 2. ; mT/mm
  DC_n = 1001
  
  A_DC_array = (FINDGEN(DC_n)/(DC_n-1)-0.5) * (FOV * gradient) / 2.
  distance_array = A_DC_array / gradient
  u_2nd = DBLARR(DC_n) * 0.
  u_3rd = DBLARR(DC_n) * 0.
  u_5th = DBLARR(DC_n) * 0.
  FOR i=0, DC_n-1 DO BEGIN
    plotYes = 0
    sine_excitation_harmonics, A, A_DC_array[i], particle_size,f,tau,H_2, M_2,frequency, u, plotYes
    u_2nd[i] = u[10]
    u_3rd[i] = u[15]
    u_5th[i] = u[25]
  ENDFOR
 
;  iplot, distance_array, u_2nd, color='red', xrange=[-20., 20.], xtitle='distance [mm]', title='2nd and 3rd harmonics'
;  iplot, distance_array, u_3rd, color='blue', /overplot
  
  maxu2nd = MAX(ABS(u_2nd[0:DC_n/2]), maxind)
  maxu2nd2 = MAX(ABS(u_2nd[DC_n/2:*]), maxind2)

  tao = -maxu2nd / ALOG(0.001)
  new_u_2nd = exp(-ABS(u_2nd)/tao)
  new_u_2nd[0:maxind] = 0;PSF_diff[0:maxind]
  new_u_2nd[DC_n/2+maxind2:*] = 0;PSF_diff[n_H_100mT1/2+maxind2:*]
  iplot, distance_array, new_u_2nd*maxu2nd, color='red', xrange=[-20., 20.], xtitle='distance [mm]', title='harmonics difference'
  
  iplot, distance_array, ABS(u_2nd), color='green', /overplot
  iplot, distance_array, ABS(u_3rd), color='blue', /overplot
  iplot, distance_array, ABS(u_5th), color='yellow', /overplot
  
  ;add another point
  offset_dist = 0.1;3.;2.;1.;3.;;0.1;2.0; mm
  offset_A_DC = offset_dist * gradient
  u_2nd_offset = DBLARR(DC_n) * 0.
  u_3rd_offset = DBLARR(DC_n) * 0.
  u_5th_offset = DBLARR(DC_n) * 0.
  FOR i=0, DC_n-1 DO BEGIN
    plotYes = 0
    sine_excitation_harmonics, A, A_DC_array[i]+offset_A_DC, particle_size,f,tau,H_2, M_2,frequency, u, plotYes
    u_2nd_offset[i] = u[10]
    u_3rd_offset[i] = u[15]
    u_5th_offset[i] = u[25]
  ENDFOR
  
  u_2nd_2p = u_2nd + u_2nd_offset
  u_3rd_2p = u_3rd + u_3rd_offset
  u_5th_2p = u_3rd + u_5th_offset
  new_u_2nd_2p = exp(-ABS(u_2nd_2p)/tao)
  new_u_2nd_2p[0:maxind] = 0;PSF_diff[0:maxind]
  new_u_2nd_2p[DC_n/2+maxind2:*] = 0;PSF_diff[n_H_100mT1/2+maxind2:*]
  iplot, distance_array, new_u_2nd_2p*maxu2nd, color='red', xrange=[-20., 20.], xtitle='distance [mm]', title='exp 2nd harmonics'
;  iplot, distance_array, ABS(u_2nd), color='green', xrange=[-20., 20.], xtitle='distance [mm]', title='2nd and 3rd harmonics'
  iplot, distance_array, ABS(u_2nd_2p), color='green', /overplot
  iplot, distance_array, ABS(u_3rd_2p), color='blue', /overplot
  iplot, distance_array, ABS(u_5th_2p), color='yellow', /overplot
END

PRO sine_excitation_harmonics, A, A_DC, particle_size,f,tau,H_2, M_2,frequency, u, plotYes, curvecolor
  tau_ms = tau/1000.

  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  Time = findgen(500)/100.     ;(1  ms)
  ;@GJ, 2022/12/14, changing to offset
  H = -A*cos(2.*!PI*f*Time)+A_DC ;(mT)
;  H = A*sin(2.*!PI*f*Time)+A_DC ;(mT)
  IF plotYes EQ 1 THEN iplot, Time, H, YRANGE=[MIN(H)*1.2, MAX(H)*1.2], thick=3, background = 'ffffff'x, color = 'green', xtitle='time (ms)', ytitle='H (mT)', title='Excitation'
  
  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
  ;  LangzwX = (FINDGEN(401)-200.)/10.
  ;  Langzw = Msat_kAm*(1./TANH(beta*LangzwX) - 1./(beta*LangzwX))
  ;  Langzw[200] = 0

  M = Msat_kAm*(1./tanh(beta*H) - 1./(beta*H)) ;(kA/m)
  R = (1./tau_ms)*exp(-Time/tau_ms)
  n = N_elements(M)
  signal = M[1:*] - M[0:n-1]
  ;@GJ, 2022/12/14, removing the NaN points
  index_nan = where(~FINITE(signal), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(signal, Time, Time[index_nan], /NAN)
    signal[index_nan] = RESULT
  ENDIF
  kernel = Debye_Kernel_sin(tau_ms, Time[0:100])


  M = convol([dblarr(n_elements(kernel)/2)*0, M, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
  M = M[n_elements(kernel)/2 : n_elements(M) - n_elements(kernel)/2 - 1]
  IF plotYes EQ 1 THEN iplot, Time, M, thick=3, background = 'ffffff'x, color = curvecolor, xtitle='time (ms)', ytitle='Magnetization', title='Magnetization Curve'

  signal_new = convol([dblarr(n_elements(kernel)/2)*0, signal, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
  signal_new = signal_new[n_elements(kernel)/2 : n_elements(signal) + n_elements(kernel)/2]
  IF plotYes EQ 1 THEN iplot, Time, signal_new, background = 'ffffff'x, color = curvecolor, thick=3, xtitle='time (ms)', ytitle='Signal', title='Signal Curve'
  
  H_2 = H[N_Elements(H)/2 : N_Elements(H)-1]  ; ��ȡһ�� ��ֹ��һ�����ھ����ֶ�����
  M_2 = M[N_Elements(M)/2 : N_Elements(M)-1]
  IF plotYes EQ 1 THEN iplot, H_2, M_2, background = 'ffffff'x, color = curvecolor, xrange=[MIN(H_2)*1.2, MAX(H_2)*1.2], thick=3, xtitle='magnetic field intensity (mT)', ytitle='Magnetization', title='M-H Curve'
  
  u = FFT(signal_new)
  ;frequency = findgen(N_Elements(u))
  ;@GJ, 2022/12/14, redefining the frequency
  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  t_ele = N_Elements(signal_new)
  delta_t = Time[1]-Time[0] ; us; samplying rate: 1MHz
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
  else $
    frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
  IF plotYes EQ 1 THEN BEGIN
    iplot, frequency/f, ABS(u), XRANGE=[0, 9], thick=3, background = 'ffffff'x, color = curvecolor, xtitle='n of f', ytitle='Amplitude', title='FFT'
    iplot, frequency/f, REAL_PART(u), color='red', PSYM=3, /overplot
    iplot, frequency/f, IMAGINARY(u), color='green', PSYM=4, /overplot
  ENDIF
  
END

;@GJ, 2023/8/15, do twin imaging based on sine excitation of nanorubbles
PRO sine_excitation_twin_imaging, A_twin, A_DC_twin, time_twin, signal, plotYes
  A = 70.; mT
  A_DC = 0.;-6.; mT
  particle_size = 30.; nm
  f = 1. ;kHz
  tau = 77.;20.; ms
  plotYesSine = 0
  curveColor = 'blue'
  sine_excitation_harmonics, A, A_DC, particle_size,f,tau,H_2, M_2, frequency, u, plotYesSine, curveColor
  
  IF N_ELEMENTS(A_twin) EQ 0 THEN A_twin = 10.; mT
  IF N_ELEMENTS(A_DC_twin) EQ 0 THEN A_DC_twin = 5.;10.;5.; mT
  min_A = -A_twin + A_DC_twin
  max_A = A_twin + A_DC_twin

  
  IF ABS(min_A) LE A AND ABS(max_A) LE A THEN BEGIN
    rise_H = H_2[25:75]
    rise_M = M_2[25:75]
    decr_H = H_2[75:125]
    decr_M = M_2[75:125]
    N_time = 500
    Time = findgen(N_time+2)/100.     ;(1  ms)
    H_twin = A_twin*sin(2.*!PI*f*Time)+A_DC_twin ;(mT)
    IF plotYes EQ 1 THEN iplot, Time, H_twin, yrange=[-A, A], xtitle='time [ms]', ytitle='H', color='green', thick = 3
    M_twin = H_twin * 0.
    M_twin_rise_temp = 0
    M_twin_decr_temp = 0
    
    FOR i=0, N_time DO BEGIN
      IF H_twin[i+1]-H_twin[i] GT 0 THEN BEGIN
        M_twin[i] = INTERPOL(rise_M, rise_H, H_twin[i])
        M_twin_rise_temp =[M_twin_rise_temp, M_twin[i]]
      ENDIF ELSE BEGIN
        M_twin[i] = INTERPOL(decr_M, decr_H, H_twin[i])
        M_twin_decr_temp =[M_twin_decr_temp, M_twin[i]]
      ENDELSE
    ENDFOR
    IF plotYes EQ 1 THEN BEGIN
      iplot, H_2, M_2, xtitle='H [mT]', ytitle='M', color='black', thick = 1, title='M-H curve'
      iplot, H_twin[0:N_time], M_twin[0:N_time], color='blue', thick = 1, /overplot
    ENDIF
    
    M_twin_rise = M_twin_rise_temp[1:*]
    M_twin_rise_range = MAX(M_twin_rise) - MIN(M_twin_rise)
    M_twin_rise_mean = (MAX(M_twin_rise) + MIN(M_twin_rise))/2.
    M_twin_decr = M_twin_decr_temp[1:*]
    M_twin_decr_range = MAX(M_twin_decr) - MIN(M_twin_decr)
    M_twin_decr_mean = (MAX(M_twin_decr) + MIN(M_twin_decr))/2.
    
    IF M_twin_rise_range LT M_twin_decr_range THEN BEGIN
      M_twin_decr_ratio = 1. / M_twin_decr_range * M_twin_rise_range
      M_twin_rise_ratio = 1.
    ENDIF ELSE BEGIN
      M_twin_decr_ratio = 1.
      M_twin_rise_ratio = 1. / M_twin_rise_range * M_twin_decr_range
    ENDELSE
    
    FOR i=0, N_time DO BEGIN
      IF H_twin[i+1]-H_twin[i] GT 0 THEN BEGIN
        M_twin[i] = (M_twin[i] - M_twin_rise_mean) * M_twin_rise_ratio
      ENDIF ELSE BEGIN
        M_twin[i] = (M_twin[i] - M_twin_decr_mean) * M_twin_decr_ratio
      ENDELSE
    ENDFOR
    
    IF plotYes EQ 1 THEN iplot, H_twin[0:N_time], M_twin[0:N_time], color='blue', thick = 3, /overplot
    
    signal = (M_twin[1:N_time] - M_twin[0:N_time-1]) / (Time[0]-Time[1])
    time_twin = time[0:N_time-1]
    IF plotYes EQ 1 THEN iplot, time_twin, signal, yrange=[-3000.,3000.], color='red', thick = 3, xtitle='time [ms]', ytitle='Signal', title ='Signal curve'
    
  ENDIF
END

;@GJ, 2023/8/15, estimation with different DC offset
;@GJ, 2023/8/18, 30 mT + 24 mT to mimick Li Lei's data
PRO twin_sine_excitation
  
  
  A_twin = 30.
  A_DC_twin_array = FINDGEN(101)/10.*2.4
  signal_max_array = A_DC_twin_array
  FOR i=0, 100 DO BEGIN
    IF i EQ 0 OR i EQ 100 THEN plotYes = 1
    sine_excitation_twin_imaging, A_twin, A_DC_twin_array[i], time_twin, signal, plotYes
    signal_max_array[i] = MAX(signal)
    plotYes = 0
  ENDFOR
  
  iplot, A_DC_twin_array, signal_max_array, color='cyan', thick = 3, xtitle='H_DC [mT]', ytitle='Maximal Signal', title ='Signal Max vs H_DC'
  
END

;@GJ, 2023/8/29, evaluate signal ratio and RESI ratio of different fields and particle sizes
;
PRO RESI_different_H_and_particle
  A_DC = 0.
  f = 1
  tau = 20  ;(ms)
  gradient_field = 40. ; (*0.1mT/mm)
  
  N_A = 400.
  N_P = 500.
  A_array = FINDGEN(N_A)/10.+1.
  Particle_size_array = FINDGEN(N_P)/10.+5.
  
  signal_ratio_array = DBLARR(N_A, N_P) * 0.
  RESI_ratio_array = DBLARR(N_A, N_P) * 0.
  RESI_reso_array = DBLARR(N_A, N_P) * 0.
  Dpp_reso_array = DBLARR(N_A, N_P) * 0.
  STD_reso_array = DBLARR(N_A, N_P) * 0.
  FOR i=0, N_A-1 DO BEGIN
    FOR j=0, N_P-1 DO BEGIN
      A = A_array[i]
      particle_size = Particle_size_array[j]
      sin_simulation, -1, A, A_DC, particle_size, f, tau, gradient_field, signal_ratio, RESI_ratio, RESI_reso, STD_reso, Dpp_reso
      signal_ratio_array[i,j]= signal_ratio
      RESI_ratio_array[i,j] = RESI_ratio
      RESI_reso_array[i,j] = RESI_reso
      Dpp_reso_array[i,j] = Dpp_reso
      STD_reso_array[i,j] = STD_reso
    ENDFOR
  ENDFOR
  
    isurface, signal_ratio_array, A_array, Particle_size_array, xtitle='Field Amplitude [mT]', ytitle='Particle Size [nm]', ztitle='Signal Ratio'
  
    isurface, RESI_ratio_array, A_array, Particle_size_array, zrange=[0.,15.], xtitle='Field Amplitude [mT]', ytitle='Particle Size [nm]', ztitle='RESI ratio'
  
    isurface, RESI_reso_array, A_array, Particle_size_array, xtitle='Field Amplitude [mT]', ytitle='Particle Size [nm]', ztitle='RESI Reso [mm]'

    isurface, Dpp_reso_array, A_array, Particle_size_array, xtitle='Field Amplitude [mT]', ytitle='Particle Size [nm]', ztitle='Dpp Reso [mm]'

    isurface, STD_reso_array, A_array, Particle_size_array, xtitle='Field Amplitude [mT]', ytitle='Particle Size [nm]', ztitle='STD Reso [mm]'
  
    print, 'RESI ratio: ', MAX(RESI_ratio_array)
  
  FOR i=0, 4 DO BEGIN
    IF i EQ 0 THEN BEGIN
      image = RESI_reso_array
      cbTitle = 'RESI Reso [mm]'
      title = 'RESI Resolution'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      maxValue = 0.5;1.;Ceil(Max(image[WHERE(FINITE(image))]))
    ENDIF
    IF i EQ 1 THEN BEGIN
      image = signal_ratio_array
      cbTitle = 'Signal Ratio'
      title = 'Signal Ratio'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 1.;
    ENDIF
    IF i EQ 2 THEN BEGIN
      image = RESI_ratio_array
      cbTitle = 'RESI Ratio'
      title = 'RESI Ratio'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
;      maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 30;15
    ENDIF
    IF i EQ 3 THEN BEGIN
      image = STD_reso_array
      cbTitle = 'STD Reso [mm]'
      title = 'STD Resolution'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 10;15
    ENDIF
    
    IF i EQ 4 THEN BEGIN
      image = Dpp_reso_array
      cbTitle = 'Dpp Reso [mm]'
      title = 'Dpp Resolution'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 2.;10;15
    ENDIF
    ; Set up variables for the contour plot. Normally, these values
    ; would be passed into the program as positional and keyword parameters.
;    image = RESI_reso_array
    
    nLevels = 10
    xtitle = 'Field Amplitude [mT]'
    ytitle = 'Particle Size [nm]'
    position =   [0.125, 0.125, 0.9, 0.800]
    cbposition = [0.125, 0.865, 0.9, 0.895]


    ; Set up a "window" for the plot. The PostScript output will have
    ; the same aspect ratio as the graphics window on the display.
;    window, i
    cgDisplay, 600, 650, Title=title, /free

    ; Set up colors for contour plot.
    cgLoadCT, 33, CLIP=[30,255];, /reverse

    ; Display the image on the display. Keep its aspect ratio.
    cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
      /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
      XRange=[MIN(A_array),MAX(A_array)], YRange=[MIN(Particle_size_array), MAX(Particle_size_array)], /Keep_Aspect

    ; Overplot the contours.
    contourLevels = cgConLevels(image, NLevels=10, MinValue=minValue, MaxValue=maxValue)
    cgContour, image, Levels=contourLevels, /OnImage, Color='charcoal'

    ; Draw the color bar. Fit it to the location of the image.
    cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
      Title=cbTitle, TLocation='Top', /Fit

  ENDFOR
  


END

;@GJ, 2023/9/1, evaluate signal ratio and RESI ratio of different relaxation time and particle sizes
;
PRO RESI_different_tau_and_particle_size
  A = 20.
  A_DC = 0.
  f = 1.
  gradient_field = 40. ; (*0.1mT/mm)

  N_tau = 501.
  N_P = 500.
  Tau_array = FINDGEN(N_tau)/10.+1. ;us 1~50us
  Particle_size_array = FINDGEN(N_P)/10.+5.

  signal_ratio_array = DBLARR(N_tau, N_P) * 0.
  RESI_ratio_array = DBLARR(N_tau, N_P) * 0.
  RESI_reso_array = DBLARR(N_tau, N_P) * 0.
  STD_reso_array = DBLARR(N_tau, N_P) * 0.
  Dpp_reso_array = DBLARR(N_tau, N_P) * 0.
  FOR i=0, N_tau-1 DO BEGIN
    FOR j=0, N_P-1 DO BEGIN
      tau = tau_array[i]
      particle_size = Particle_size_array[j]
      sin_simulation, -1, A, A_DC, particle_size, f, tau, gradient_field, signal_ratio, RESI_ratio, RESI_reso, STD_reso, Dpp_reso
      signal_ratio_array[i,j]= signal_ratio
      RESI_ratio_array[i,j] = RESI_ratio
      RESI_reso_array[i,j] = RESI_reso
      STD_reso_array[i,j] = STD_reso
      Dpp_reso_array[i,j] = Dpp_reso
    ENDFOR
  ENDFOR

  isurface, signal_ratio_array, tau_array, Particle_size_array, xtitle='Relaxation Time [us]', ytitle='Particle Size [nm]', ztitle='Signal Ratio'

  isurface, RESI_ratio_array, tau_array, Particle_size_array, xtitle='Relaxation Time [us]', ytitle='Particle Size [nm]', ztitle='RESI ratio'

  isurface, RESI_reso_array, tau_array, Particle_size_array, xtitle='Relaxation Time [us]', ytitle='Particle Size [nm]', ztitle='RESI Reso [mm]'

  isurface, Dpp_reso_array, tau_array, Particle_size_array, xtitle='Relaxation Time [us]', ytitle='Particle Size [nm]', ztitle='Dpp Reso [mm]'

  isurface, STD_reso_array, tau_array, Particle_size_array, xtitle='Relaxation Time [us]', ytitle='Particle Size [nm]', ztitle='STD Reso [mm]'

  print, 'RESI ratio: ', MAX(RESI_ratio_array)

  FOR i=0, 4 DO BEGIN
    IF i EQ 0 THEN BEGIN
      image = RESI_reso_array
      cbTitle = 'RESI Reso [mm]'
      title = 'RESI Resolution'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      maxValue = 0.5;1.;Ceil(Max(image[WHERE(FINITE(image))]))
    ENDIF
    IF i EQ 1 THEN BEGIN
      image = signal_ratio_array
      cbTitle = 'Signal Ratio'
      title = 'Signal Ratio'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 1.;
    ENDIF
    IF i EQ 2 THEN BEGIN
      image = RESI_ratio_array
      cbTitle = 'RESI Ratio'
      title = 'RESI Ratio'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;      maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 50;15
    ENDIF
    IF i EQ 3 THEN BEGIN
      image = STD_reso_array
      cbTitle = 'STD Reso [mm]'
      title = 'STD Resolution'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 10;15
    ENDIF
    
    IF i EQ 4 THEN BEGIN
      image = Dpp_reso_array
      cbTitle = 'Dpp Reso [mm]'
      title = 'Dpp Resolution'
      minValue = MAX([0., Floor(Min(image[WHERE(FINITE(image))]))])
      ;maxValue = Ceil(Max(image[WHERE(FINITE(image))]))
      maxValue = 2.;10;15
    ENDIF

    ; Set up variables for the contour plot. Normally, these values
    ; would be passed into the program as positional and keyword parameters.
    ;    image = RESI_reso_array

    nLevels = 10
    xtitle = 'Relaxation Time [us]'
    ytitle = 'Particle Size [nm]'
    position =   [0.125, 0.125, 0.9, 0.800]
    cbposition = [0.125, 0.865, 0.9, 0.895]


    ; Set up a "window" for the plot. The PostScript output will have
    ; the same aspect ratio as the graphics window on the display.
    ;    window, i
    cgDisplay, 600, 650, Title=title, /free

    ; Set up colors for contour plot.
    cgLoadCT, 33, CLIP=[30,255];, /reverse

    ; Display the image on the display. Keep its aspect ratio.
    cgImage, image, Stretch=1, MinValue=minValue, MaxValue=maxValue, $
      /Axes, XTitle=xtitle, YTitle=ytitle, Position=position, $
      XRange=[MIN(tau_array),MAX(tau_array)], YRange=[MIN(Particle_size_array), MAX(Particle_size_array)], /Keep_Aspect

    ; Overplot the contours.
    contourLevels = cgConLevels(image, NLevels=10, MinValue=minValue, MaxValue=maxValue)
    cgContour, image, Levels=contourLevels, /OnImage, Color='charcoal'

    ; Draw the color bar. Fit it to the location of the image.
    cgColorbar, Position=cbposition, Range=[MinValue, MaxValue], $
      Title=cbTitle, TLocation='Top', /Fit

  ENDFOR



END

;@GJ, 2023/8/30, Evaluating the merging of Lion Eyes
PRO RESI_lion_left_right_eye_merging

  ;@GJ, lion eye
  fn_dcm = 'F:\MPI_Tianjie\Nanorubble_MPI\DICOM_images\003_dandu teshutanzhen\dandu teshutanzhen.dcm'
  afn_des = FILE_DIRNAME(fn_dcm)
  IF QUERY_DICOM(fn_dcm) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = fn_dcm
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(fn_dcm)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
    fn_dcm=fn_MPI_new
  ENDIF
  IF QUERY_DICOM(fn_dcm) EQ 1 THEN BEGIN
    image_dcm_ori = read_dicom(fn_dcm)
    ;  iimage, image_dcm, title='original'
    obj = OBJ_NEW('IDLffDICOM', fn_dcm)
    pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
    pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
    ; Get the row & column size of the image(s):
    temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = N_ELEMENTS(image_dcm_ori[0,*])
    temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = N_ELEMENTS(image_dcm_ori[*,0])

    OBJ_DESTROY, obj
  ENDIF ELSE BEGIN
    return
  ENDELSE
  
;  iimage, image_dcm_ori, title='lion eye image'
  angle = 180./!PI*ATAN((114.-112.)/(119.-104.))
  center = [(114.+112.)/2., (119.+104.)/2.]
  
  image_dcm = ROT(image_dcm_ori, angle/2., 1.0, center[0], center[1], /INTERP)
  iimage, image_dcm, title='corrected Lion eye'
  
  X = (113.+158.)/2.;127;244; 0;127;
  print, 'X =', X, ' (cols:', cols,'), pixel space ', pixelSp, 'mm'
  image_left = DOUBLE(image_dcm);DOUBLE(image_dcm[0:X, *]);/DOUBLE(MAX(image_dcm[0:X, *]))*5.
  image_left[X:*,*] = 0.
  image_right = DOUBLE(image_dcm);DOUBLE(image_dcm[X:*, *]);/DOUBLE(MAX(image_dcm[X:*, *]))*5.
  image_right[0:X,*] = 0.
;  WINDOW, 1, XSIZE = cols, YSIZE = rows
;  TV, BYTSCL(image_left)
;  WINDOW, 2, XSIZE = cols, YSIZE = rows
;  TV, BYTSCL(image_right)
  
  ;read phantom image
  fn_dir = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\'
  ;fn_name = 'S_shape.png'
  ;  fn_name = 'letter32.png'
  ;  fn_name = 'SUA_letter24.png'
  ;  fn_name = 'Mammo2.png'
  fn_name = 'Mammo4.png'
  ;  fn_name = 'Mammo5.jpg'
  fn_name = 'MRA000.png'
  fn = fn_dir + fn_name
  
  IF QUERY_PNG(fn) EQ 1 THEN BEGIN
    image_temp = read_png(fn)
    image_sim_temp = DOUBLE(REFORM(image_temp[0,*,*]))
    image_sim = MAX(image_sim_temp) - image_sim_temp
    IF STRCMP(fn_name, 'MRA', 3) THEN image_sim = DOUBLE(image_sim LT 192) * 255.
    rows_S = N_ELEMENTS(image_sim[0,*])
    cols_S = N_ELEMENTS(image_sim[*,0])
    pixelSp = [1.0, 1.0]
  ENDIF ELSE BEGIN
    return
  ENDELSE
  ;@GJ, modify the phantom data
  ;  image_sim *= 0.
  ;  image_sim[rows_S/2-3:rows_S/2+3, cols_S/2-30:cols_S/2+30] = 255
  ;  image_sim[rows_S/2-30:rows_S/2+30, cols_S/2-3:cols_S/2+3] = 255
  
  ;display the image phantom
  WINDOW, 0, XSIZE = cols_S, YSIZE = rows_S
  TV, BYTSCL(image_sim)
  
  n_shift_array = INDGEN((158.-113.)/2.)+10
  n_shift_array = FINDGEN(10.)+3.
  FOR k=0, N_ELEMENTS(n_shift_array)-1 DO BEGIN
    n_shift = n_shift_array[k]
    ;define the directory
    fitplot_dir = FILE_DIRNAME(fn_dcm, /MARK_DIRECTORY)+'RESI_shift'+STRTRIM(STRING(FLOOR(n_shift)),1)+'\'
    IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
    
    image_merge = SHIFT(image_left, n_shift, 0) + SHIFT(image_right, -n_shift, 0)
    image_merge[X-n_shift+1:X+n_shift-1, *] /= 2.
;    iimage, image_merge[X-n_shift+1:X+n_shift-1, *], title='Added'

    ;do convolution
    left_kernel = ROT(SHIFT(image_left, n_shift, 0), 0, 0.1)
    ;phantom_left = CONVOL(image_sim, left_kernel, /EDGE_ZERO, /NORMALIZE)
    phantom_left = CONVOL_FFT(image_sim, left_kernel)
    right_kernel = ROT(SHIFT(image_right, -n_shift, 0), 0, 0.1)
    iimage, left_kernel+right_kernel, title='Kernel'
    write_png, fitplot_dir + FILE_BASENAME(fn_name, '.png') + 'kernel_merge.png', BYTSCL(left_kernel+right_kernel)
    ;phantom_right = CONVOL(image_sim, right_kernel, /EDGE_ZERO, /NORMALIZE)
    phantom_right = CONVOL_FFT(image_sim, right_kernel)
    phantom_merge = phantom_left + phantom_right
    iimage, phantom_merge, title='Phantom Merge'
    write_png, fitplot_dir + FILE_BASENAME(fn_name, '.png') + 'phantom_merge.png', BYTSCL(phantom_merge)
    phantom_subtr = ABS(phantom_left - phantom_right)
    iimage, phantom_subtr, title='Phantom Subtr'
    write_png, fitplot_dir + FILE_BASENAME(fn_name, '.png') + 'phantom_subtr.png', BYTSCL(phantom_subtr)
    tao = -MAX(phantom_subtr) / ALOG(0.001)
    phantom_subtr_exp = exp(-phantom_subtr/tao)
    iimage, phantom_subtr_exp, title='Phantom Subtr Exp'
    write_png, fitplot_dir + FILE_BASENAME(fn_name, '.png') + 'phantom_subtr_exp.png', BYTSCL(phantom_subtr_exp)
    ;rescale the image based on the expo transform
    merged_exp_phantom = phantom_merge
    FOR xi=0, cols_S-1 DO FOR yi=0, rows_S-1 DO merged_exp_phantom[xi,yi] *= phantom_subtr_exp[xi,yi]/MAX(phantom_subtr_exp)
    iimage, merged_exp_phantom, title='Merged Phantom Subtr Exp'
    write_png, fitplot_dir + FILE_BASENAME(fn_name, '.png') + '_merged_phantom_subtr_exp.png', BYTSCL(merged_exp_phantom)
    
    image_RESI = ABS(left_kernel - right_kernel)
;    iimage, image_RESI[X-n_shift+1:X+n_shift-1, *], title='RESI Subtracted'
    min_RESI = MIN(image_RESI[X-n_shift+1:X+n_shift-1, *])
    max_RESI = MAX(image_RESI[X-n_shift+1:X+n_shift-1, *])
    ;@GJ, 2023/8/30, calculating the exponential
    tao = -max_RESI / ALOG(0.001)
    exp_RESI_map = exp(-image_RESI[X-n_shift+1:X+n_shift-1, *]/tao)
    min_exp_RESI = MIN(exp_RESI_map)
    max_exp_RESI = MAX(exp_RESI_map)
   
    ;do surface plot of the FOV
    merged_exp_RESI_image = image_merge[X-n_shift+1:X+n_shift-1, *]
    FOR xi=0, N_ELEMENTS(merged_exp_RESI_image[*,0])-1 DO BEGIN
      FOR yi=0, N_ELEMENTS(merged_exp_RESI_image[0,*])-1 DO BEGIN
        merged_exp_RESI_image[xi,yi] *= exp_RESI_map[xi,yi]/max_exp_RESI
      ENDFOR
    ENDFOR
    window, 1
    cgLoadCT, 33, RGB_Table=rainbowPalette
    cgSurf, image_merge[X-n_shift+1:X+n_shift-1, *], ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
    it_fig = TVRD(TRUE=1)
    WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'ori_PSF_surface.png', it_fig
    window, 2
    cgSurf, merged_exp_RESI_image, ROTX=10, ROTZ=20, ZCharsize=cgDefCharsize()*1.25, /Elevation, Palette=rainbowPalette
    it_fig = TVRD(TRUE=1)
    WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'exp_PSF_surface.png', it_fig
    
    ;@GJ, 2023/9/5, do convolution and check results
    
    
;    FOR i=100, 0, -10 DO BEGIN
;      transp = i
;      ;  Image_Blend, BYTSCL(image_merge[X-n_shift+1:X+n_shift-1, *]), BYTSCL(image_RESI[X-n_shift+1:X+n_shift-1, *]), Colortable=25, blendTitle='RESI'
;      WINDOW, 3, XSIZE = 2*n_shift-1, YSIZE = rows, title='RESI Image'
;      ;  cgDisplay, 2*n_shift-1, rows, title='RESI Image', free=1
;      cgImage, BYTSCL(image_merge[X-n_shift+1:X+n_shift-1, *]), CTIndex=0
;      
;      ;    transp = 50
;      cgImage, BYTSCL(image_RESI[X-n_shift+1:X+n_shift-1, *], MIN=min_RESI, MAX=max_RESI, TOP=255), CTIndex=33, Transparent=transp;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      it_fig = TVRD(TRUE=1)
;      ;save the signal plot
;      WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_RESI_'+STRTRIM(STRING(i),1)+'.png', it_fig
;      WDELETE, 3
;
;      WINDOW, 4, XSIZE = 2*n_shift-1, YSIZE = rows, title='expRESI Image'
;      ;  cgDisplay, 2*n_shift-1, rows, title='RESI Image', free=1
;      cgImage, BYTSCL(image_merge[X-n_shift+1:X+n_shift-1, *]), CTIndex=0
;      transp = i
;      cgImage, BYTSCL(exp_RESI_map, MIN=min_exp_RESI, MAX=max_exp_RESI, TOP=255), CTIndex=33, Transparent=transp;, output='PNG', outfilename=cur_dir+'\Debye_map1.png'
;      it_fig = TVRD(TRUE=1)
;      WRITE_PNG, fitplot_dir+FILE_BASENAME(fn_dcm, '.dcm')+'_expRESI_'+STRTRIM(STRING(i),1)+'.png', it_fig
;      WDELETE, 4
;    ENDFOR
    
;    IF k EQ N_ELEMENTS(n_shift_array)-1 THEN BEGIN
     IF n_shift EQ 19 THEN BEGIN
      ;  iplot, FINDGEN(2*n_shift-1)*pixelSp[0], image_merge[X-n_shift+1:X+n_shift-1, rows/2]-MIN(image_merge[X-n_shift+1:X+n_shift-1, rows/2]), color='blue', xtitle='Distance [mm]', ytitle='PSF', title='PSF'
      ;  iplot, FINDGEN(2*n_shift-1)*pixelSp[0], image_RESI[X-n_shift+1:X+n_shift-1, rows/2],color='green', PSYM=3, /overplot
      ;  iplot, FINDGEN(2*n_shift-1)*pixelSp[0], exp_RESI_map[*, rows/2]/max_exp_RESI*(max(image_merge[X-n_shift+1:X+n_shift-1, *])-MIN(image_merge[X-n_shift+1:X+n_shift-1, rows/2])),color='red', thick=2, /overplot
      ;
      ;
      ; Display the first plot.
      plot1 = PLOT(FINDGEN(2*n_shift-1)*pixelSp[0], image_merge[X-n_shift+1:X+n_shift-1, rows/2]-MIN(image_merge[X-n_shift+1:X+n_shift-1, rows/2]), 'b2', xtitle='Distance [mm]', ytitle='PSF', title='PSF', NAME='Standard')
      ; Display the second plot.
      plot2 = PLOT(FINDGEN(2*n_shift-1)*pixelSp[0], image_RESI[X-n_shift+1:X+n_shift-1, rows/2], /OVERPLOT, 'g2', NAME='RESI')
      plot3 = PLOT(FINDGEN(2*n_shift-1)*pixelSp[0], exp_RESI_map[*, rows/2]/max_exp_RESI*(max(image_merge[X-n_shift+1:X+n_shift-1, *])-MIN(image_merge[X-n_shift+1:X+n_shift-1, rows/2])), /OVERPLOT, 'r2', NAME='expRESI')
      ; Add the legend.
      leg = LEGEND(TARGET=[plot1,plot2,plot3], POSITION=[185,0.9], /DATA, /AUTO_TEXT_COLOR)
    ENDIF
  ENDFOR
  
  
  
END


PRO RESI_4_angles_analysis
  testfile = DIALOG_PICKFILE(FILTER='*.csv', path='C:\D_drive\MPI_Tianjie\RESI_twin_imaging\')
  aaa= read_csv(testfile)
  n_tags=N_TAGS(aaa)
  tag_array = TAG_NAMES(aaa)
  
  data = DBLARR(n_tags, N_ELEMENTS(aaa.(tag_array[0])))
  FOR i=0, 29 DO BEGIN
    data[i,*] = aaa.(tag_array[i])
  ENDFOR
  iimage, data
END

PRO RESI_offset
  A=20.;20.;16.;4.;20.;20.;4.;20.;
;  A_DC = 0.
;  RESI_donut_Orth, A, A_DC, Harmonic_3rd_Amp_array_n, Harmonic_3rd_Phase_array_n, Harmonic_2nd_Amp_array_n, Harmonic_2nd_Phase_array_n, Donut_n, Donut2_n, phantom_2D_array, z_axis, filename
  ;    iimage, Harmonic_3rd_Amp_array_n, title='No Offset 3rd Amp image'
  
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  
;  A_DC_array = FINDGEN(30)
;  FOR i=16, 21-1 DO BEGIN
;    A_DC = A_DC_array[i]
;    RESI_donut_Orth, A, A_DC, Harmonic_3rd_Amp_array_n, Harmonic_3rd_Phase_array_n, Harmonic_2nd_Amp_array_n, Harmonic_2nd_Phase_array_n, Donut_n, Donut2_n, phantom_2D_array, z_axis, filename
;    iimage, Harmonic_3rd_Amp_array_n, title='No Offset 3rd Amp image'
;  ENDFOR
;  
;  print, 'finished circling'
    
  A_DC_s = 19.;20.;20.;21.;9.;A*0.9;19.;14.;A-0.5;1.;0.5;2.;1.;0.5;
  RESI_donut_Orth, A, A_DC_s, Harmonic_3rd_Amp_array_s, PSF_Harmonic_3rd_Amp_array_s, rec_A3_Amp_s, Harmonic_3rd_Phase_array_s, Harmonic_2nd_Amp_array_s, Harmonic_2nd_Phase_array_s, Donut_s, Donut2_s, phantom_2D_array, PSF_phantom, z_axis, filename
  ;@GJ, 2023/10/26, based on the searching results
  A_DC_l = 21.;-20.;22.;16.;20.;A+0.5;1.;0.5;1.;0.5;
  RESI_donut_Orth, A, A_DC_l, Harmonic_3rd_Amp_array_l, PSF_Harmonic_3rd_Amp_array_l, rec_A3_Amp_l, Harmonic_3rd_Phase_array_l, Harmonic_2nd_Amp_array_l, Harmonic_2nd_Phase_array_l, Donut_l, Donut2_l, phantom_2D_array, PSF_phantom, z_axis, filename

  IF OBJ_VALID(p) THEN p.close
  m = N_ELEMENTS(phantom_2D_array[*,0])/2.
  delta_Donut = Donut_s - Donut_l
  p = PLOT(z_axis, REVERSE(delta_Donut[m, *], 2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='AmpDonut3 [AU]', xtitle='z axis [mm]', title='A3Donut (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC_s, format='(f5.1)')+'-'+STRING(A_DC_l, format='(f5.1)')+'mT')
  p.name = ' Amp*sinPhi'
  p1 = PLOT(z_axis, DOUBLE(phantom_2D_array[m, *])/MAX(phantom_2D_array[m, *])*2.*STDDEV(delta_Donut[m, *])+MIN(delta_Donut[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p1.name = ' Phantom'
  l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  ;      fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
  ;      IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
  dir_name = '.\RESI_donut_pics\AC'+STRING(A, format='(f5.1)')+'mT_'
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  output_filename_temp = dir_name+current_time+'_delta_Harmonic_3rd_AmpDonut_curve.png'
  p.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  ;close the plot
  p.close
  
  IF OBJ_VALID(p0) THEN p0.close
  delta_phase = Harmonic_3rd_Phase_array_s - Harmonic_3rd_Phase_array_l
  p0 = PLOT(z_axis, REVERSE(delta_phase[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Phi3 [deg]', xtitle='z axis [mm]', title='Phi3 (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC_s, format='(f5.1)')+'-'+STRING(A_DC_l, format='(f5.1)')+'mT')
  p0.name = ' Phase'
  p01 = PLOT(z_axis, DOUBLE(phantom_2D_array[m, *])/MAX(phantom_2D_array[m, *])*2.*STDDEV(delta_phase[m, *])+MIN(delta_phase[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p01.name = ' Phantom'
  l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  output_filename_temp = dir_name+current_time+'_delta_Harmonic_3rd_Phase_curve.png'
  p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  p0.close
  
  IF OBJ_VALID(p1) THEN p1.close
  delta_amp = Harmonic_3rd_Amp_array_s - Harmonic_3rd_Amp_array_l
  p1 = PLOT(z_axis, REVERSE(delta_amp[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='A3  [AU]', xtitle='z axis [mm]', title='A3 (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC_s, format='(f5.1)')+'-'+STRING(A_DC_l, format='(f5.1)')+'mT')
  p1.name = ' Amp'
  p11 = PLOT(z_axis, DOUBLE(phantom_2D_array[m, *])/MAX(phantom_2D_array[m, *])*2.*STDDEV(delta_amp[m, *])+MIN(delta_amp[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p11.name = ' Phantom'
  l1 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  output_filename_temp = dir_name+current_time+'_delta_Harmonic_3rd_Amp_curve.png'
  p1.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  
  IF OBJ_VALID(p2) THEN p2.close
  delta_PSF_amp = PSF_Harmonic_3rd_Amp_array_s - PSF_Harmonic_3rd_Amp_array_l
  p2 = PLOT(z_axis, REVERSE(delta_PSF_amp[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='A3  [AU]', xtitle='z axis [mm]', title='A3 PSF (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC_s, format='(f5.1)')+'-'+STRING(A_DC_l, format='(f5.1)')+'mT')
  p2.name = ' PSF Amp'
  p21 = PLOT(z_axis, DOUBLE(PSF_phantom[m, *])/MAX(PSF_phantom[m, *])*2.*STDDEV(delta_PSF_amp[m, *])+MIN(delta_PSF_amp[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p21.name = ' PSF'
  l2 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  output_filename_temp = dir_name+current_time+'_delta_PSF_Harmonic_3rd_Amp_curve.png'
  p2.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT

  IF OBJ_VALID(p0) THEN p0.close
  delta_phase2 = Harmonic_2nd_Phase_array_s - Harmonic_2nd_Phase_array_l
  p0 = PLOT(z_axis, REVERSE(delta_phase2[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Phi2 (deg)', xtitle='z axis [mm]', title='Phi2 (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC_s, format='(f5.1)')+'-'+STRING(A_DC_l, format='(f5.1)')+'mT')
  p0.name = ' Phase'
  p01 = PLOT(z_axis, DOUBLE(phantom_2D_array[m, *])/MAX(phantom_2D_array[m, *])*2.*STDDEV(delta_phase2[m, *])+MIN(delta_phase2[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p01.name = ' Phantom'
  l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
  output_filename_temp = dir_name+current_time+'_delta_Harmonic_2nd_Phase_curve.png'
  p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  p0.close
  
  ;@GJ, 2023/11/1, save the original image
  iimage, BYTSCL(REVERSE(phantom_2D_array,2)), title='Original image'
  WRITE_PNG, dir_name+current_time+'_Original_image.png', BYTSCL(REVERSE(phantom_2D_array,2))
  
  subtract_Donut = Donut_s - Donut_l
  scaled_subtract_Donut = BYTSCL(subtract_Donut, max=MAX(subtract_Donut[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Donut[m/2:m+m/2-1, m/2:m+m/2-1]))
 ; iimage, scaled_subtract_Donut, title='substracted 3rd Donut image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_AmpDonut_Subtracted.png', scaled_subtract_Donut
  
;  subtract_Donut2 = Donut2_s - Donut2_l
;  scaled_subtract_Donut2 = BYTSCL(subtract_Donut2, max=MAX(subtract_Donut2[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Donut2[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Donut2, title='substracted 2nd Donut image'
;  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_AmpDonut_Subtracted.png', scaled_subtract_Donut2

  subtract_Amp3 = Harmonic_3rd_Amp_array_s - Harmonic_3rd_Amp_array_l
  scaled_subtract_Amp3 = BYTSCL(subtract_Amp3, max=MAX(subtract_Amp3[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Amp3[m/2:m+m/2-1, m/2:m+m/2-1]))
  iimage, scaled_subtract_Amp3, title='substracted 3rd Amp image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_Amp_Subtracted.png', scaled_subtract_Amp3
 
  subtract_PSF_Amp3 = PSF_Harmonic_3rd_Amp_array_s - PSF_Harmonic_3rd_Amp_array_l
  scaled_subtract_PSF_Amp3 = BYTSCL(subtract_PSF_Amp3, max=MAX(subtract_PSF_Amp3[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_PSF_Amp3[m/2:m+m/2-1, m/2:m+m/2-1]))
  iimage, scaled_subtract_PSF_Amp3, title='substracted PSF 3rd Amp image'
  WRITE_PNG, dir_name+current_time+'_A3_PSF_Amp_Subtracted.png', scaled_subtract_PSF_Amp3
  
  ;
  z_ele = N_Elements(z_axis)
  delta_z = z_axis[1] - z_axis[0] ; mm; samplying rate: 1MHz
  Z = (FINDGEN((z_ele - 1)/2) + 1)
  is_N_even = (z_ele MOD 2) EQ 0
  if (is_N_even) then $
    frequency_z = [0.0, Z, z_ele/2, -z_ele/2 + Z]/(z_ele*delta_z) $
  else $
    frequency_z = [0.0, Z, -(z_ele/2 + 1) + Z]/(z_ele*delta_z)
  PSF_A3_Amp_FFT = ABS(FFT(subtract_PSF_Amp3))
  A3_Amp_FFT = FFT(subtract_Amp3)
  A3_Amp_FFT_div = A3_Amp_FFT/(PSF_A3_Amp_FFT);+MAX(A3_Amp_FFT)*0.01)
  rec_A3_Amp = FFT(A3_Amp_FFT_div, /INVERSE)
  p0 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(PSF_A3_Amp_FFT[m, *]))[SORT(frequency_z)], 'r-2', xrange=[0, MAX(frequency_z)], YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Contrast', xtitle='Spatial Frequancy [lp/mm]', title='FFT (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC_s, format='(f5.1)')+'-'+STRING(A_DC_l, format='(f5.1)')+'mT')
  p0.name = ' PSF MTF'
  p01 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A3_Amp_FFT[m, *])))[SORT(frequency_z)], 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p01.name = ' Phantom FFT'
  p02 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A3_Amp_FFT_div[m, *])))[SORT(frequency_z)], 'g-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p02.name = ' Phantom Rec'
  l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  output_filename_temp = dir_name+current_time+'_delta_Harmonic_3rd_Amp_FFT_curves.png'
  p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  
  
  
  iimage, BYTSCL(A3_Amp_FFT, max=MAX(A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1])), title='Delta Phantom 3rd Amp (FFT)'
  iimage, BYTSCL(PSF_A3_Amp_FFT, max=MAX(PSF_A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(PSF_A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1])), title='Delta PSF 3rd Amp (FFT)'
  iimage, BYTSCL(rec_A3_Amp, max=MAX(rec_A3_Amp[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(rec_A3_Amp[m/2:m+m/2-1, m/2:m+m/2-1])), title='Delta Rec Phantom 3rd Amp (FFT)'
  WRITE_PNG, dir_name+current_time+'_delta_A3_Amp_phantom_FFT.png', BYTSCL(A3_Amp_FFT, max=MAX(A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1]))
  WRITE_PNG, dir_name+current_time+'_delta_A3_Amp_PSF_FFT.png', BYTSCL(PSF_A3_Amp_FFT, max=MAX(PSF_A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(PSF_A3_Amp_FFT[m/2:m+m/2-1, m/2:m+m/2-1]))
  WRITE_PNG, dir_name+current_time+'_delta_A3_Amp_phantom_rec_FFT.png', BYTSCL(rec_A3_Amp, max=MAX(rec_A3_Amp[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(rec_A3_Amp[m/2:m+m/2-1, m/2:m+m/2-1]))
  
  rec_A3_Amp_substraction = rec_A3_Amp_s - rec_A3_Amp_l
  iimage, BYTSCL(rec_A3_Amp_substraction, max=MAX(rec_A3_Amp_substraction[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(rec_A3_Amp_substraction[m/2:m+m/2-1, m/2:m+m/2-1])), title='Subtraction Rec Phantom 3rd Amp (FFT)'
  WRITE_PNG, dir_name+current_time+'_subtraction_A3_Amp_phantom_FFT.png', BYTSCL(rec_A3_Amp_substraction, max=MAX(rec_A3_Amp_substraction[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(rec_A3_Amp_substraction[m/2:m+m/2-1, m/2:m+m/2-1]))

;  subtract_Amp2 = Harmonic_2nd_Amp_array_s - Harmonic_2nd_Amp_array_l
;  scaled_subtract_Amp2 = BYTSCL(subtract_Amp2, max=MAX(subtract_Amp2[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Amp2[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Amp2, title='substracted 2nd Amp image'
;  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_Amp_Subtracted.png', scaled_subtract_Amp2
  
  subtract_Phase3 = Harmonic_3rd_Phase_array_s - Harmonic_3rd_Phase_array_l
  scaled_subtract_Phase3 = BYTSCL(subtract_Phase3, max=MAX(subtract_Phase3[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Phase3[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Phase3, title='substracted 3rd Phase image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_Phase_Subtracted.png', scaled_subtract_Phase3
  
  subtract_Phase3_A3 = subtract_Phase3/TOTAL(subtract_Phase3)*(subtract_Amp3-MIN(subtract_Amp3))
  scaled_subtract_Phase3_A3 = BYTSCL(subtract_Phase3_A3, max=MAX(subtract_Phase3_A3[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Phase3_A3[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Phase3_A3, title='substracted 3rd Phase*Amp image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_PhaseAmp_Subtracted.png', scaled_subtract_Phase3_A3
 
  subtract_Amp2 = Harmonic_2nd_Amp_array_s - Harmonic_2nd_Amp_array_l
  scaled_subtract_Amp2 = BYTSCL(subtract_Amp2, max=MAX(subtract_Amp2[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Amp2[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Amp2, title='substracted 2nd Amp image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_Amp_Subtracted.png', scaled_subtract_Amp2

  subtract_Phase2 = Harmonic_2nd_Phase_array_s - Harmonic_2nd_Phase_array_l
  scaled_subtract_Phase2 = BYTSCL(subtract_Phase2, max=MAX(subtract_Phase2[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Phase2[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Phase2, title='substracted 2nd Phase image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_Phase_Subtracted.png', scaled_subtract_Phase2
  
  subtract_Phase2_A2 = subtract_Phase2/TOTAL(subtract_Phase2)*(subtract_Amp2-MIN(subtract_Amp2))
  scaled_subtract_Phase2_A2 = BYTSCL(subtract_Phase2_A2, max=MAX(subtract_Phase2_A2[m/2:m+m/2-1, m/2:m+m/2-1]), MIN=MIN(subtract_Phase2_A2[m/2:m+m/2-1, m/2:m+m/2-1]))
;  iimage, scaled_subtract_Phase2_A2, title='substracted 2nd Phase*Amp image'
  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_PhaseAmp_Subtracted.png', scaled_subtract_Phase2_A2
  
  

END

PRO RESI_offset_CP
  A = 20.
  A_DC = A/2.;
  RESI_donut_Orth, A, A_DC, Harmonic_3rd_Amp_array_s, PSF_Harmonic_3rd_Amp_array_s, Harmonic_3rd_Phase_array_s, Harmonic_2nd_Amp_array_s, Harmonic_2nd_Phase_array_s, Donut_s, Donut2_s, phantom_2D_array, z_axis, filename
  A_DC = A/2.+0.5;
  RESI_donut_Orth, A, A_DC, Harmonic_3rd_Amp_array_l, PSF_Harmonic_3rd_Amp_array_l, Harmonic_3rd_Phase_array_l, Harmonic_2nd_Amp_array_l, Harmonic_2nd_Phase_array_l, Donut_l, Donut2_l, phantom_2D_array, z_axis, filename

  iimage, Harmonic_3rd_Amp_array_s - Harmonic_3rd_Amp_array_l, title='substracted 3rd Amp image'
;  iimage, Harmonic_3rd_Phase_array_s - Harmonic_3rd_Phase_array_l, title='substracted 3rd Phase image'
  
;  iimage, Harmonic_2nd_Amp_array_s - Harmonic_2nd_Amp_array_l, title='substracted 2nd Amp image'
;  iimage, Harmonic_3rd_Phase_array_s - Harmonic_3rd_Phase_array_l, title='substracted 3rd Phase image'

END

PRO RESI_0offset
  A = 20.
  A_DC = 0.;A/2.;
  RESI_donut_Orth, A, A_DC, Harmonic_3rd_Amp_array_s, PSF_Harmonic_3rd_Amp_array_s, Harmonic_3rd_Phase_array_s, Harmonic_2nd_Amp_array_s, Harmonic_2nd_Phase_array_s, Donut_s, Donut2_s, phantom_2D_array, z_axis, filename

END

;@GJ, 2023/10/25, search for donut
PRO Donut_search_offset
  f = 1
  gradient_field = 20. ; (*0.1mT/mm)
  gradient_field_y = gradient_field
  gradient_field_z = gradient_field
  signal_ratio = 0.;
  RESI_ratio = 1.;
  STD_reso = 1.;
  RESI_reso = 1.;
  Dpp_reso = 1.;
  noise_level = 0.01;
  N_points = 1.;2.;1.;2.;2.;1.;2.;
  D_points = 5.;20.;100.;20.;20.;5.;30.;10.; *0.1mm
  fov_phantom = 250.;
  n_angles = 4.;
  angle = 90.;
  shift_perc = 43.;100.;

  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m

  n_Time = 1000;5000.
  Time = findgen(n_Time)/1000.     ;(1  ms)
  

  ;@GJ, 2023/9/8, plot the 2D PSF figure
  FOV = 2. * 20. / (gradient_field / 10.) ; in mm
  n_x = FLOOR(n_Time/(5.*2.*f)) ;this is corrected number
  print, 'FOV = ', STRING(FOV,format='(f6.2)')+' mm'
  ;@GJ, 2022/12/14, changing to offset
  n_y = n_x
  n_z = n_x
  resolution = [FOV/n_x, FOV/n_y, FOV/n_z]
  print, 'resolution: ', resolution
  y_axis = FINDGEN(n_y) * resolution[1]
  z_axis = FINDGEN(n_z) * resolution[2]

  ;@GJ, 2023/10/11, design the phantom
  phantom_1D_array = DBLARR(n_y*4.) * 0.
  N_D_points = FLOOR(D_points*0.1/resolution[1])
  Npoints_ind = 0.
  FOR i=0, n_y*4.-1 DO BEGIN
    IF (i MOD N_D_points) EQ 0 AND Npoints_ind LT N_points THEN  BEGIN
      phantom_1D_array[i] = 1;i+1;1.
      Npoints_ind++
    ENDIF
  ENDFOR
  phantom_1D_array = SHIFT(phantom_1D_array, n_y/2.-(N_D_points*(N_points-1)/2.))
  ;  iplot, z_axis, phantom_1D_array[0:n_y-1], xtitle='z axis [mm]', title='1D phantom'
  phantom_2D_array = DBLARR(n_y, n_z) * 0.
  phantom_2D_array[n_y/2, *] = phantom_1D_array[0:n_y-1]
  ;  iimage, phantom_2D_array

  

  particle_size_array = [30.];DBLARR(50)+1.
  tau_array = [20.];FINDGEN(50)+1.
  tau_ms_array = tau_array/1000.
  A_array = [20.];FINDGEN(30)+1.;[16]
  A_DC_array = FINDGEN(30);/2.;-10.
  
  donut_radius = DBLARR(N_ELEMENTS(particle_size_array), N_ELEMENTS(tau_ms_array), N_ELEMENTS(A_array), N_ELEMENTS(A_DC_array))
  donut2_radius = donut_radius * 0.
  donut_thickness = donut_radius * 0.
  donut2_thickness = donut2_radius * 0.
  donut_outer_FWHM = donut_radius * 0.
  non_donut_FWHM = donut_radius * 0.
  
  FOR id=0, N_ELEMENTS(particle_size_array)-1 DO BEGIN
    FOR jd=0, N_ELEMENTS(tau_ms_array)-1 DO BEGIN
      FOR kd=0, N_ELEMENTS(A_array)-1 DO BEGIN
        FOR ld=0, N_ELEMENTS(A_DC_array)-1 DO BEGIN
          print, 'A_ind, A_DC_ind=', kd, ',', ld
          particle_size = particle_size_array[id]
          tau_ms = tau_ms_array[jd]
          A = A_array[kd]
          A_DC = A_DC_array[ld]
          print, 'A, A_DC=', A, ',', A_DC, '(mT)'
          beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
          R = (1/tau_ms)*exp(-Time/tau_ms)
          kernel = Debye_Kernel_sin(tau_ms, Time[0:100])

          H_y_array_L = (FINDGEN(n_y*2.)-FLOOR(n_y)) * resolution[1] * (gradient_field_y/10.)
          H_z_array_L = (FINDGEN(n_z*2.)-FLOOR(n_z)) * resolution[2] * (gradient_field_z/10.)
          signal_x_array = DBLARR(n_y, n_z, n_Time-1) * 0.
          Harmonic_3rd_Amp_array = DBLARR(n_y, n_z) * 0.
          Harmonic_3rd_Phase_array = DBLARR(n_y, n_z) * 0.
          Harmonic_3rd_Real_array = DBLARR(n_y, n_z) * 0.
          Harmonic_3rd_Imaginary_array = DBLARR(n_y, n_z) * 0.
          Harmonic_2nd_Amp_array = DBLARR(n_y, n_z) * 0.
          Harmonic_2nd_Phase_array = DBLARR(n_y, n_z) * 0.
          m = n_y/2.
          FOR n=0, n_z-1 DO BEGIN
            H_y_array = H_y_array_L[m:m+n_y-1]
            H_z_array = H_z_array_L[n:n+n_y-1]
            M_x = DBLARR(n_Time, n_y, n_z) * 0.
            FOR i=0, n_y-1 DO BEGIN
              FOR j=0, n_z-1 DO BEGIN
                IF phantom_2D_array[i, j] GT 0 THEN BEGIN
                  ;print, 'm, n=', m, ',', n
                  ;@GJ, 2023/10/20, adding the DC field
                  H_x = -A*cos(2.*!PI*f*Time) + A_DC
                  H_y = H_y_array[i]
                  H_z = H_z_array[j]
                  H_yz = SQRT(H_y^2 + H_z^2)
                  H_total = SIGNUM(H_x) * SQRT(H_x^2 + H_yz^2)

                  ;calculate M
                  M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
                  index_nan = where(~FINITE(M_total), count)
                  IF count GT 0 THEN BEGIN
                    RESULT = INTERPOL(M_total, Time, Time[index_nan], /NAN)
                    M_total[index_nan] = RESULT
                  ENDIF

                  ;do vector splitting
                  ;@GJ, 2023/9/23, simple calculation with the signal intensity
                  M_x[*,i,j] = ABS(M_total) * H_x/ABS(H_total) * phantom_2D_array[i, j]

                  signal_x = -(M_x[1:*, i, j] - M_x[0:n_Time-2, i, j])/(Time[1]-Time[0])
                  ;@GJ, 2022/12/14, removing the NaN points
                  index_nan_x = where(~FINITE(signal_x), count_x)
                  IF count_x GT 0 THEN BEGIN
                    RESULT_x = INTERPOL(signal_x, Time, Time[index_nan_x], /NAN)
                    signal_x[index_nan_x] = RESULT_x
                  ENDIF
                  signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
                  signal_new_x = signal_new_x[n_elements(kernel)/2 : n_elements(signal_x) + n_elements(kernel)/2]
               
                  signal_x_array[m, n, *] += signal_new_x
                ENDIF
              ENDFOR
            ENDFOR
            
            t_ele = N_Elements(signal_x_array[m, n, *])
            delta_t = (Time[1]-Time[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
            X = (FINDGEN((t_ele - 1)/2) + 1)
            is_N_even = (t_ele MOD 2) EQ 0
            if (is_N_even) then $
              frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
            else $
              frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
            base_freq = MIN(ABS(frequency-1.0), base_ind)

            signal_FFT_x = FFT(signal_x_array[m, n, *])
            ;@GJ, 2023/10/20, get the harmonic amp and phase
            Harmonic_3rd_Amp_array[m, n] = ABS(signal_FFT_x[base_ind*3.])
            Harmonic_3rd_Real_array[m, n] = Real_part(signal_FFT_x[base_ind*3.])
            Harmonic_3rd_Imaginary_array[m, n] = Imaginary(signal_FFT_x[base_ind*3.])
            Harmonic_3rd_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*3.]), Real_part(signal_FFT_x[base_ind*3.]), /phase)

          ENDFOR
          
          ;calculate the donut diameter and thickness
          donut = Harmonic_3rd_Phase_array[m, *];Harmonic_3rd_Amp_array[m, *] * SIN(Harmonic_3rd_Phase_array[m, *])
          
          IF MEAN(donut) LT 0 THEN donut *= -1
          maxD = MAX(donut[0:n_z/2-1], maxD_ind)
          IF maxD_ind LT n_z/2-1 AND donut[n_z/2-1] LT donut[n_z/2-3]  THEN BEGIN
            ;there is a donut
            donut_radius[id, jd, kd, ld] = (n_z/2-1 - maxD_ind) * resolution[2] ; in mm
            donut_thickness[id, jd, kd, ld] = maxD - donut[n_z/2-1] ; in random unit
            FWHM_ind = INTERPOL(z_axis[0:n_z/2-1], donut[0:n_z/2-1], 0.5*maxD)
            donut_outer_FWHM[id, jd, kd, ld] = 2 * (n_z/2-1 - FWHM_ind)
          ENDIF ELSE BEGIN
            FWHM_ind = INTERPOL(z_axis[0:n_z/2-1], donut[0:n_z/2-1], 0.5*maxD)
            non_donut_FWHM[id, jd, kd, ld] = 2 * (n_z/2-1 - FWHM_ind)
            donut_outer_FWHM[id, jd, kd, ld] = 2 * (n_z/2-1 - FWHM_ind)
          ENDELSE

         
        ENDFOR
 
        ;@GJ, 2023/11/1, plot the result
        IF OBJ_VALID(p) THEN p.close
        p = PLOT(A_DC_array, REFORM(donut_radius[id, jd, kd, *]), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
          XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Size [mm]', xtitle='A_DC [mT]', title='Donut A3 (AC'+STRING(A, format='(f5.1)')+'mT, Particle'+STRING(Particle_Size, format='(i2)')+'nm)')
        p.name = ' Inner Radius'
        p1 = PLOT(A_DC_array, REFORM(donut_thickness[id, jd, kd, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
        p1.name = ' Inner Thickness'
        l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.5, 0.85])
        ;      fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
        ;      IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
        dir_name = '.\RESI_donut_pics\Par'+STRING(Particle_Size, format='(i3)')+'nm_AC'+STRING(A, format='(f5.1)')+'mT_'
        temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
        ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
        current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
        output_filename_temp = dir_name+current_time+'_Harmonic_3rd_AmpDonut_curve.png'
        p.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
        ;close the plot
        
        IF OBJ_VALID(p0) THEN p0.close
        p0 = PLOT(A_DC_array, REFORM(donut2_radius[id, jd, kd, *]), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
          XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Size [mm]', xtitle='A_DC [mT]', title='Donut A2 (AC'+STRING(A, format='(f5.1)')+'mT, Particle'+STRING(Particle_Size, format='(i2)')+'nm)')
        p0.name = ' Inner Radius'
        p01 = PLOT(A_DC_array, REFORM(donut2_thickness[id, jd, kd, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
        p01.name = ' Inner Thickness'
        l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.5, 0.85])
        ;      fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
        ;      IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
        dir_name = '.\RESI_donut_pics\Par'+STRING(Particle_Size, format='(i3)')+'nm_AC'+STRING(A, format='(f5.1)')+'mT_'
        temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
        ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
        current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
        output_filename_temp = dir_name+current_time+'_Harmonic_2nd_AmpDonut_curve.png'
        p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
        ;close the plot

        
;        IF id+jd+kd EQ 0 THEN BEGIN
;          iplot, A_DC_array, REFORM(donut_radius[0, 0, 0, *]), color='red', xtitle='A_DC [mT]', ytitle='Donut Radius [mm]', title='Donut Radius and non-Donut FWHM'
;          iplot, A_DC_array, donut_thickness[0, 0, 0, *], color='blue', /overplot
;;          iplot, A_DC_array, donut2_radius[0, 0, 0, *], color='blue', /overplot
;          ;iplot, A_DC_array, non_donut_FWHM[0, 0, 0, *], color='blue', /overplot
; ;         iplot, A_DC_array, donut_thickness[0, 0, 0, *], color='blue', xtitle='A_DC [mT]', ytitle='Donut Thickness [AU]', title='Donut Thickness'
;        ENDIF
        
      ENDFOR
      
      IF id+jd EQ 0 THEN BEGIN
        isurface, REFORM(donut_radius[0, 0, *, *]), A_array, A_DC_array, xtitle='A [mT]', ytitle='A_DC [mT]', ztitle='Donut Radius [mm]', title='DonutA3 Radius'
 ;       isurface, REFORM(non_donut_FWHM[0, 0, *, *]), A_array, A_DC_array, xtitle='A [mT]', ytitle='A_DC [mT]', ztitle='Donut Thickness [AU]', title = 'Donut Thickness'
        isurface, REFORM(donut_thickness[0, 0, *, *]), A_array, A_DC_array, xtitle='A [mT]', ytitle='A_DC [mT]', ztitle='Donut Thickness [AU]', title = 'DonutA3 Thickness'
        isurface, REFORM(donut2_radius[0, 0, *, *]), A_array, A_DC_array, xtitle='A [mT]', ytitle='A_DC [mT]', ztitle='Donut Radius [mm]', title='DonutA2 Radius'
        isurface, REFORM(donut2_thickness[0, 0, *, *]), A_array, A_DC_array, xtitle='A [mT]', ytitle='A_DC [mT]', ztitle='Donut Thickness [AU]', title='DonutA2 Thickness'
      ENDIF
    ENDFOR
  ENDFOR

;  isurface, donut_diameter[0, 0, *, *], A_array, A_DC_array, xtitle='A', ytitle='A_DC', ztitle='Donut Diameter'
;  iplot, A_DC_array, donut_diameter[0, 0, *, *], color='red', xtitle='A_DC [mT]', ytitle='Donut Diameter [mm]', title='Donut Diameter, A=20 mT'
;  iplot, A_DC_array, donut_thickness[0, 0, *, *], color='blue', xtitle='A_DC [mT]', ytitle='Donut Thickness', title='Donut Thickness, A=20 mT'
END

PRO Phase_contrast_MPI
  ;initial setting
  A = 19.;4.;20.
  ;  A_DC = 8.;10.;21.;10.;(*0.1mT)
  particle_size = 30.;5.;30
  f = 1.;kHz
  tau = 20.;0.;0.01;30.;1.;30.;3.;20.;3.;20  ;(us)
  flat_portion = 50 ;
  gradient_field = 20. ; (*0.1mT/mm)
  gradient_field_y = gradient_field
  gradient_field_z = gradient_field
  signal_ratio = 0.;
  RESI_ratio = 1.;
  STD_reso = 1.;
  RESI_reso = 1.;
  Dpp_reso = 1.;
  noise_level = 0.01;
  N_points = 2.;2.;3.;2.;2.;2.;2.;3.;3.;3.;2.;2.;3.;3.;1.;2.;3.;2.;2.;2.;0.;0.;3.;2.;1.;2.;2.;1.;2.;
  D_points = 10.;5;20.;10.;20.;20.;100.;20.;20.;5.;30.;10.; *0.1mm
  fov_phantom = 250.;
  n_angles = 4.;
  angle = 90.;
  shift_perc = 43.;100.;

  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m

  n_Time = 1000;5000.
  Time = findgen(n_Time)/1000.     ;(1  ms)
  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT

  tau_ms = tau/1000.
  R = (1/tau_ms)*exp(-Time/tau_ms)
  kernel = Debye_Kernel_sin(tau_ms, Time[0:100])

  ;@GJ, 2023/9/8, plot the 2D PSF figure
  FOV = 2. * 20. / (gradient_field / 10.) ; in mm
  n_x = FLOOR(n_Time/(5.*2.*f)) ;this is corrected number
  print, 'FOV = ', STRING(FOV,format='(f6.2)')+' mm'
  ;@GJ, 2022/12/14, changing to offset
  n_y = n_x
  n_z = n_x
  resolution = [FOV/n_x, FOV/n_y, FOV/n_z]
  print, 'resolution: ', resolution
  x_axis = FINDGEN(n_x) * resolution[0]
  y_axis = FINDGEN(n_y) * resolution[1]
  z_axis = FINDGEN(n_z) * resolution[2]
  H_DC_array = (FINDGEN(n_x*2.+1)-n_x) * resolution[0] * (gradient_field / 10.); + 10.

  ;@GJ, 2023/10/11, design the phantom
  phantom_1D_array = DBLARR(n_y*4.) * 0.
  N_D_points = FLOOR(D_points*0.1/resolution[1])
  Npoints_ind = 0.
  FOR i=0, n_y*4.-1 DO BEGIN
    IF (i MOD N_D_points) EQ 0 AND Npoints_ind LT N_points THEN  BEGIN
      phantom_1D_array[i] = 1.;i+1;1.
      Npoints_ind++
    ENDIF
  ENDFOR
  phantom_1D_array = (SHIFT(phantom_1D_array, n_y/2.-(N_D_points*(N_points-1)/2.)))[0:n_x-1]
  iplot, x_axis, phantom_1D_array, xtitle='x [mm]', ytitle='Signal'
  
  ;two phases for subtraction
  n_H_G = 2.
  H_G_array = [0., 1.] ; mT
  A3_Amp_array = DBLARR(n_H_G, n_x)
  A3_Phase_array = DBLARR(n_H_G, n_x)
  FOR i=0, n_H_G-1 DO BEGIN
    H_G = H_G_array[i]
    ;FFP movement
    signal_x_array = DBLARR(n_x, n_Time-1) * 0.
    FOR j=0, n_x-1 DO BEGIN
      A_DC_array = H_DC_array[j:n_x-1+j]
      ;signal acquisition
      M_x = DBLARR(n_Time, n_x) * 0.
      FOR k=0, n_x-1 DO BEGIN
        A_DC = A_DC_array[k]
        H_x = -A*cos(2.*!PI*f*Time) + A_DC
        H_y = 0.
        H_z = H_G
        H_yz = SQRT(H_y^2 + H_z^2) + 14.;12.;14.;4.;7.;9.; offset
        H_total = SIGNUM(H_x) * SQRT(H_x^2 + H_yz^2)
;        print, 'Hy, Hz', H_y, H_z, ' mT'
        
        ;calculate M
        M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
        index_nan = where(~FINITE(M_total), count)
        IF count GT 0 THEN BEGIN
          RESULT = INTERPOL(M_total, Time, Time[index_nan], /NAN)
          M_total[index_nan] = RESULT
        ENDIF
        
        ;do vector splitting
        ;@GJ, 2023/9/23, simple calculation with the signal intensity
        M_x[*,k] = ABS(M_total) * H_x/ABS(H_total) * phantom_1D_array[k]
        
        signal_x = -(M_x[1:*, k] - M_x[0:n_Time-2, k])/(Time[1]-Time[0])
        ;@GJ, 2022/12/14, removing the NaN points
        index_nan_x = where(~FINITE(signal_x), count_x)
        IF count_x GT 0 THEN BEGIN
          RESULT_x = INTERPOL(signal_x, Time, Time[index_nan_x], /NAN)
          signal_x[index_nan_x] = RESULT_x
        ENDIF
        signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
        signal_new_x = signal_new_x[n_elements(kernel)/2 : n_elements(signal_x) + n_elements(kernel)/2]
        ;            signal_new_x[n_Time/2:*] = 0.
        signal_x_array[j, *] += signal_new_x
      ENDFOR
      
      t_ele = N_Elements(signal_x_array[j, *])
      delta_t = (Time[1]-Time[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
      X = (FINDGEN((t_ele - 1)/2) + 1)
      is_N_even = (t_ele MOD 2) EQ 0
      if (is_N_even) then $
        frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
      else $
        frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
      base_freq = MIN(ABS(frequency-1.0), base_ind)

      signal_FFT_x = FFT(signal_x_array[j, *])
      ;@GJ, 2023/10/20, get the harmonic amp and phase
      H_ind = 3.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      A3_Amp_array[i,j] = ABS(signal_FFT_x[base_ind*H_ind])
      A3_Phase_array[i,j] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*H_ind]), Real_part(signal_FFT_x[base_ind*H_ind]), /phase)
    ENDFOR
;    isurface, signal_x_array
;    iimage, BYTSCL(signal_x_array)
  ENDFOR
  
  iplot, x_axis, A3_Phase_array[0,*], color='blue', title='A3 Phase'
  iplot, x_axis, A3_Phase_array[1,*], color='red', /overplot
  iplot, x_axis, A3_Phase_array[1,*] - A3_Phase_array[0,*], title='A3 Phase Diff'

  iplot, x_axis, A3_Amp_array[0,*], color='blue', title='A3 Amp'
  iplot, x_axis, A3_Amp_array[1,*], color='red', /overplot
  iplot, x_axis, A3_Amp_array[1,*] - A3_Amp_array[0,*], title='A3 Amp Diff'

END

;@GJ, 2023/10/25, search for donut
;@GJ, 2023/11/15, do raster scan
;@GJ, 2023/11/17, calculate donut size
;@GJ, 2023/11/23, using the parameters is Zhong Jing's poster
;@GJ, 2023/11/25, save all 21 harmonics
;@GJ, 2023/12/16, output all A3 Amp and Phase, calculate donut radius and FWHM
PRO Mag_Phase_DC_G_scan
  f = 10.;kHz
  gradient_field = 11.36 ; (*0.1mT/mm)
  gradient_field_y = gradient_field
  gradient_field_z = gradient_field

  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m

  n_Time = 1000;5000.
  Time = findgen(n_Time)/1000.     ;(1  ms)

  particle_size_array = [30.];[100.];[30.];DBLARR(50)+1.
  tau_array = [20.];[0.1];[20.];FINDGEN(50)+1.
  tau_ms_array = tau_array/1000.
  A_array = FINDGEN(20)/1.+1.;[10.];FINDGEN(11.)/2.+8.;FINDGEN(200)/10.+1.;[20.];[20.];FINDGEN(200)/10.+1.;/10.+17.;[16]
  A_DC_array = FINDGEN(161)/4.-20.;(FINDGEN(61)-30.)/30.*A_array[0]*4./3.;FINDGEN(30);/10.;/2.;-10.
  A_G_array = FINDGEN(161)/4.-20.;30.*A_array[0]*4./3.;/10.;/2.;-10.
  
  n_Harmonic_Amp_array = 21
  Harmonic_Amp_array = DBLARR(n_Harmonic_Amp_array, N_ELEMENTS(particle_size_array), N_ELEMENTS(tau_ms_array), N_ELEMENTS(A_array), N_ELEMENTS(A_DC_array), N_ELEMENTS(A_G_array))
  Harmonic_Phase_array = Harmonic_Amp_array * 0.
  Harmonic_Imaginary_array = Harmonic_Amp_array * 0.
  Harmonic_Real_array = Harmonic_Amp_array * 0.
  donut_Harmonic_radius = DBLARR(n_Harmonic_Amp_array, N_ELEMENTS(particle_size_array), N_ELEMENTS(tau_ms_array), N_ELEMENTS(A_array), N_ELEMENTS(A_DC_array))
  donut_Harmonic_thick = donut_Harmonic_radius * 0.
  non_donut_Harmonic_FWHM = donut_Harmonic_radius * 0.
  non_donut_Harmonic_Height = donut_Harmonic_radius * 0.

  A3_Amp_array = DBLARR(N_ELEMENTS(particle_size_array), N_ELEMENTS(tau_ms_array), N_ELEMENTS(A_array), N_ELEMENTS(A_DC_array), N_ELEMENTS(A_G_array))
  A3_Phase_array = A3_Amp_array * 0.
  A3_Real_array = A3_Amp_array * 0.
  A3_Imaginary_array = A3_Amp_array * 0.
  donut_A3_radius = DBLARR(N_ELEMENTS(particle_size_array), N_ELEMENTS(tau_ms_array), N_ELEMENTS(A_array), N_ELEMENTS(A_DC_array))
  donut_A3_thick = DBLARR(N_ELEMENTS(particle_size_array), N_ELEMENTS(tau_ms_array), N_ELEMENTS(A_array), N_ELEMENTS(A_DC_array))

  FOR id=0, N_ELEMENTS(particle_size_array)-1 DO BEGIN
    particle_size = particle_size_array[id]
    FOR jd=0, N_ELEMENTS(tau_ms_array)-1 DO BEGIN
      tau_ms = tau_ms_array[jd]
      ;@GJ, 2023/11/6, calculate the tau and relaxation kernel
      beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
      R = (1/tau_ms)*exp(-Time/tau_ms)
      kernel = Debye_Kernel_sin(tau_ms, Time[0:100])
      FOR kd=0, N_ELEMENTS(A_array)-1 DO BEGIN
        A = A_array[kd]
        FOR ld=0, N_ELEMENTS(A_DC_array)-1 DO BEGIN
          A_DC = A_DC_array[ld]
          print, 'A, A_DC=', A, ',', A_DC, '(mT)'
          FOR md=0, N_ELEMENTS(A_G_array)-1 DO BEGIN
;            print, 'A_ind, A_DC_ind, A_G_ind=', kd, ',', ld, ',', md
            H_x = -A*cos(2.*!PI*f*Time) + A_DC
;            H_x = A*sin(2.*!PI*f*Time) + A_DC
            H_G = A_G_array[md]
            H_total = SIGNUM(H_x) * SQRT(H_x^2 + H_G^2)

            ;calculate M
            M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
            index_nan = where(~FINITE(M_total), count)
            IF count GT 0 THEN BEGIN
              RESULT = INTERPOL(M_total, Time, Time[index_nan], /NAN)
              M_total[index_nan] = RESULT
            ENDIF

            ;do vector splitting
            ;@GJ, 2023/9/23, simple calculation with the signal intensity
            M_x = ABS(M_total) * H_x/ABS(H_total)

            signal_x = -(M_x[1:*] - M_x[0:n_Time-2])/(Time[1]-Time[0])
            ;@GJ, 2022/12/14, removing the NaN points
            index_nan_x = where(~FINITE(signal_x), count_x)
            IF count_x GT 0 THEN BEGIN
              RESULT_x = INTERPOL(signal_x, Time, Time[index_nan_x], /NAN)
              signal_x[index_nan_x] = RESULT_x
            ENDIF
            signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
            signal_new_x = signal_new_x[n_elements(kernel)/2 : n_elements(signal_x) + n_elements(kernel)/2]

            t_ele = N_Elements(signal_new_x)
            delta_t = (Time[1]-Time[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
            X = (FINDGEN((t_ele - 1)/2) + 1)
            is_N_even = (t_ele MOD 2) EQ 0
            if (is_N_even) then $
              frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
            else $
              frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
            base_freq = MIN(ABS(frequency-1.0), base_ind)

            signal_FFT_x = FFT(signal_new_x)
            ;@GJ, 2023/10/20, get the harmonic amp and phase
            FOR hd = 0, n_Harmonic_Amp_array-1 DO BEGIN
              ;@GJ, 2023/11/9, calculate the A2 components
              Harmonic_Amp_array[hd, id, jd, kd, ld, md] = ABS(signal_FFT_x[base_ind*(hd+1)])
              Harmonic_Real_array[hd, id, jd, kd, ld, md] = Real_part(signal_FFT_x[base_ind*(hd+1)])
              Harmonic_Imaginary_array[hd, id, jd, kd, ld, md] = Imaginary(signal_FFT_x[base_ind*(hd+1)])
              Harmonic_Phase_array[hd, id, jd, kd, ld, md] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*(hd+1)]), Real_part(signal_FFT_x[base_ind*(hd+1)]), /phase)
            ENDFOR
           ENDFOR
          
          FOR hd = 0, n_Harmonic_Amp_array-1 DO BEGIN
            ;calculate the donut diameter and thickness
            zero_ind=WHERE(ABS(A_G_Array) EQ 0)
            donut_Harmonic = REFORM(Harmonic_Amp_array[hd, id, jd, kd, ld, zero_ind:*])
            temp_A_G_array = A_G_array[zero_ind:*]
            maxD_Harmonic = MAX(donut_Harmonic, maxD_ind_Harmonic)
            IF maxD_ind_Harmonic GT 2. THEN BEGIN
              donut_Harmonic_radius[hd, id, jd, kd, ld] = temp_A_G_array[maxD_ind_Harmonic]
              donut_Harmonic_thick[hd, id, jd, kd, ld] = donut_Harmonic[maxD_ind_Harmonic] - donut_Harmonic[0]
            ENDIF ELSE BEGIN
              ;non-donut, calculate FWHM
              FWHM_ind = INTERPOL(temp_A_G_array, donut_Harmonic, 0.5*maxD_Harmonic)
              non_donut_Harmonic_FWHM[hd, id, jd, kd, ld] = 2. * FWHM_ind
              non_donut_Harmonic_Height[hd, id, jd, kd, ld] = maxD_Harmonic
            ENDELSE
          ENDFOR
        ENDFOR
      ENDFOR
    ENDFOR
  ENDFOR
  
  dir_name = '.\RESI_donut_pics\'
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  LOADCT, 8
  w_size = 512
  x = CONGRID(A_DC_array, w_size)/A_array[N_ELEMENTS(A_array)/2-1]
  y = CONGRID(A_G_array, N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size)/A_array[N_ELEMENTS(A_array)/2-1]
  
  FOR hd = 0, n_Harmonic_Amp_array-1 DO BEGIN
    temp = BYTSCL(REFORM(Harmonic_Amp_array[hd, 0, 0, N_ELEMENTS(A_array)/2-1, *, *]))
    temp_int = CONGRID(temp, w_size, N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size, /INTERP)
    WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
    PX = !X.WINDOW * !D.X_VSIZE
    PY = !Y.WINDOW * !D.Y_VSIZE
    SX=PX[1]-PX[0]+1
    SY=PY[1]-PY[0]+1
    TVSCL,CONGRID(temp_int,SX,SY),PX[0],PY[0]
    CONTOUR,temp_int,X,Y,NLEVELS=10,$
      XSTYLE=1,YSTYLE=1,$
      C_LINESTYLE=[1,0],$
      C_THICK=[1,1,1,1,1,3],$
      TITLE=STRING(hd+1, format='(i2)')+' Harmonic Amp',$
      XTITLE ='H_DC/H_AC', YTITLE='H_G/H_AC', $
      /NoErase
    it_fig = TVRD(TRUE=1)
    ;save the image plot
    output_filename_temp = dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_Harmonic_Amp_'+STRING(hd+1, format='(i02)')+'.png'
    WRITE_PNG, output_filename_temp, it_fig
    wdelete, 0
  ENDFOR
  
  ;@GJ, 2023/12/10, plot the surface NC Fig. 2
;  x_plot=FINDGEN(10)+1; # harnomics
;  y_plot=A_G_array[N_ELEMENTS(A_G_array)/2:*]; Gradient arrays
;  z_plot=REFORM(Harmonic_Amp_array[0:9,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,N_ELEMENTS(A_G_array)/2:*])
;  my=surface(z_plot, x_plot, y_plot, title='Harmonics Amp')
;  cgSURFACE, z_plot, X_plot, Y_plot, xtitle='Harmonics', ytitle='H_G [mT]', ztitle='A3 Amp', title='Harmonics Amp'
;;  myContour = CONTOUR(z_plot, x_plot, y_plot, N_LEVELS=15, /ZVALUE, PLANAR=0, OVERPLOT=my)
;  iplot, y_plot, REFORM(Harmonic_Amp_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,N_ELEMENTS(A_G_array)/2:*]), PSYM=3
;  ;@GJ, 2023/12/10, plot the surface NC Fig. 2
;  iplot, A_G_array, REFORM(Harmonic_Amp_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*]), LINESTYLE=0, color='red', xtitle='H_G', ytitle='A3 Phi3(A_DC:'+STRING(A_DC_array[N_ELEMENTS(A_DC_array)*3./4.+4], format='(i02)')+' mT)'
;  iplot, A_G_array, REFORM(Harmonic_Phase_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*])/MAX(REFORM(Harmonic_Phase_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*]))*MAX(REFORM(Harmonic_Amp_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*])), LINESTYLE=2, color='blue',/overplot
;  iplot, A_G_array, REFORM(Harmonic_Real_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*])/MAX(REFORM(Harmonic_Real_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*]))*MAX(REFORM(Harmonic_Amp_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*])), LINESTYLE=3, /overplot, color='green'
;  iplot, A_G_array, REFORM(Harmonic_Imaginary_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*])/MAX(REFORM(Harmonic_Imaginary_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*]))*MAX(REFORM(Harmonic_Amp_array[2,0,0,0,N_ELEMENTS(A_DC_array)*3./4.+4,*])), LINESTYLE=4, /overplot, color='black'
;
;  iplot, A_DC_array, Harmonic_Phase_array[0, 0, 0, N_ELEMENTS(A_array)/2-1, *, N_ELEMENTS(A_G_array)/2], PSYM = 4, xrange=[0, MAX(A_DC_array)], xtitle='H_DC [mT]', ytitle='A1 Phase [A.U.]', title='Phase Plot'
;  iplot, A_DC_array, Harmonic_Phase_array[1, 0, 0, N_ELEMENTS(A_array)/2-1, *, N_ELEMENTS(A_G_array)/2], PSYM = 4, xrange=[0, MAX(A_DC_array)], xtitle='H_DC [mT]', ytitle='A2 Phase [A.U.]', title='Phase Plot'
;  iplot, A_DC_array, Harmonic_Phase_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, *, N_ELEMENTS(A_G_array)/2], PSYM = 4, xrange=[0, MAX(A_DC_array)], xtitle='H_DC [mT]', ytitle='A3 Phase [A.U.]', title='Phase Plot'
;  iplot, A_DC_array, Harmonic_Phase_array[3, 0, 0, N_ELEMENTS(A_array)/2-1, *, N_ELEMENTS(A_G_array)/2], PSYM = 4, xrange=[0, MAX(A_DC_array)], xtitle='H_DC [mT]', ytitle='A4 Phase [A.U.]', title='Phase Plot'
  
  LOADCT, 2
  WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
  FOR hd = 0, n_Harmonic_Amp_array-1 DO BEGIN
    temp = BYTSCL(REFORM(Harmonic_Phase_array[hd, 0, 0, N_ELEMENTS(A_array)/2-1, *, *]))
    temp_int = CONGRID(temp, w_size, N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size)
    WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
    PX = !X.WINDOW * !D.X_VSIZE
    PY = !Y.WINDOW * !D.Y_VSIZE
    SX=PX[1]-PX[0]+1
    SY=PY[1]-PY[0]+1
    TVSCL,CONGRID(temp_int,SX,SY),PX[0],PY[0]
    CONTOUR,temp_int,X,Y,NLEVELS=10,$
      XSTYLE=1,YSTYLE=1,$
      C_LINESTYLE=[1,0],$
      C_THICK=[1,1,1,1,1,3],$
      TITLE=STRING(hd+1, format='(i2)')+' Harmonic Phase',$
      XTITLE ='H_DC/H_AC', YTITLE='H_G/H_AC', $
      /NoErase
    it_fig = TVRD(TRUE=1)
    ;save the image plot
    output_filename_temp = dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_Harmonic_Phase_'+STRING(hd+1, format='(i02)')+'.png'
    WRITE_PNG, output_filename_temp, it_fig
    wdelete, 0
  ENDFOR

  LOADCT, 3
  WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
  FOR hd = 0, n_Harmonic_Amp_array-1 DO BEGIN
    temp = BYTSCL(REFORM(Harmonic_Real_array[hd, 0, 0, N_ELEMENTS(A_array)/2-1, *, *]))
    temp_int = CONGRID(temp, w_size, N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size)
    WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
    PX = !X.WINDOW * !D.X_VSIZE
    PY = !Y.WINDOW * !D.Y_VSIZE
    SX=PX[1]-PX[0]+1
    SY=PY[1]-PY[0]+1
    TVSCL,CONGRID(temp_int,SX,SY),PX[0],PY[0]
    CONTOUR,temp_int,X,Y,NLEVELS=10,$
      XSTYLE=1,YSTYLE=1,$
      C_LINESTYLE=[1,0],$
      C_THICK=[1,1,1,1,1,3],$
      TITLE=STRING(hd+1, format='(i2)')+' Harmonic Real',$
      XTITLE ='H_DC/H_AC', YTITLE='H_G/H_AC', $
      /NoErase
    it_fig = TVRD(TRUE=1)
    ;save the image plot
    output_filename_temp = dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_Harmonic_Real_'+STRING(hd+1, format='(i02)')+'.png'
    WRITE_PNG, output_filename_temp, it_fig
    wdelete, 0
  ENDFOR

  LOADCT, 1
  WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
  FOR hd = 0, n_Harmonic_Amp_array-1 DO BEGIN
    temp = BYTSCL(REFORM(Harmonic_Imaginary_array[hd, 0, 0, N_ELEMENTS(A_array)/2-1, *, *]))
    temp_int = CONGRID(temp, w_size, N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size)
    WINDOW, 0, XSIZE = w_size, YSIZE = N_ELEMENTS(A_G_array)/ N_ELEMENTS(A_DC_array)*w_size
    PX = !X.WINDOW * !D.X_VSIZE
    PY = !Y.WINDOW * !D.Y_VSIZE
    SX=PX[1]-PX[0]+1
    SY=PY[1]-PY[0]+1
    TVSCL,CONGRID(temp_int,SX,SY),PX[0],PY[0]
    CONTOUR,temp_int,X,Y,NLEVELS=10,$
      XSTYLE=1,YSTYLE=1,$
      C_LINESTYLE=[1,0],$
      C_THICK=[1,1,1,1,1,3],$
      TITLE=STRING(hd+1, format='(i2)')+' Harmonic Real',$
      XTITLE ='H_DC/H_AC', YTITLE='H_G/H_AC', $
      /NoErase
    it_fig = TVRD(TRUE=1)
    ;save the image plot
    output_filename_temp = dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_Imaginary_Phase_'+STRING(hd+1, format='(i02)')+'.png'
    WRITE_PNG, output_filename_temp, it_fig
    wdelete, 0
  ENDFOR
  
 ; DEVICE, DECOMPOSED=1
  
;  iimage, BYTSCL(REFORM(Harmonic_Amp_array[1, 0, 0, N_ELEMENTS(A_array)/2-1, *, *])), title='A2 amp (AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f4.1)')+'mT)'
;  iimage, BYTSCL(REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, *, *])), title='A3 amp (AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f4.1)')+'mT)'
;  iimage, BYTSCL(REFORM(Harmonic_Amp_array[3, 0, 0, N_ELEMENTS(A_array)/2-1, *, *])), title='A4 amp (AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f4.1)')+'mT)'
  
;  isurface, REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, *, *]), A_DC_array, A_G_array, xtitle='H_DC [mT]', ytitle='H_G [mT]', ztitle='A3 Amp [a.u.]', title='A3 Amp'
  
;  DEVICE, DECOMPOSED=0
  LOADCT, 1
;  WINDOW, 0, XSIZE = w_size, YSIZE = w_size
  H_Gy_array = A_G_array
  H_Gz_array = A_G_array
  n_y = N_ELEMENTS(H_Gy_array)
  n_z = N_ELEMENTS(H_Gz_array)
  FOR k=0, N_ELEMENTS(A_DC_array)-1 DO BEGIN
    A_DC = A_DC_array[k]
    focal_image_0 = DBLARR(N_ELEMENTS(H_Gy_array), N_ELEMENTS(H_Gz_array)) * 0.
    focal_image_0_Phase = DBLARR(N_ELEMENTS(H_Gy_array), N_ELEMENTS(H_Gz_array)) * 0.
    A3_amp_array_0 = REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, k, *])
    A3_Phase_array_0 = REFORM(Harmonic_Phase_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, k, *])
    FOR i=0, N_ELEMENTS(H_Gy_array)-1 DO BEGIN
      FOR j=0, N_ELEMENTS(H_Gz_array)-1 DO BEGIN
        H_G = SQRT(H_Gy_array[i]^2 + H_Gz_array[j]^2)
        focal_image_0[i,j] = INTERPOL(A3_amp_array_0, A_G_array, H_G)
        focal_image_0_Phase[i,j] = INTERPOL(A3_Phase_array_0, A_G_array, H_G)
      ENDFOR
    ENDFOR
    IF A_DC LT 0 THEN BEGIN
      focal_image_fn = dir_name+current_time+'_A3Spot_A_DCm'+STRING(ABS(A_DC_array[k]), format='(f05.1)')+'mT.png'
      focal_image_Phase_fn = dir_name+current_time+'_A3PhaseSpot_A_DCm'+STRING(ABS(A_DC_array[k]), format='(f05.1)')+'mT.png'
    ENDIF ELSE BEGIN
      focal_image_fn = dir_name+current_time+'_A3Spot_A_DC'+STRING(ABS(A_DC_array[k]), format='(f05.1)')+'mT_A3Spot.png'
      focal_image_Phase_fn = dir_name+current_time+'_A3PhaseSpot_A_DC'+STRING(ABS(A_DC_array[k]), format='(f05.1)')+'mT.png'
    ENDELSE
    WRITE_PNG, focal_image_fn, BYTSCL(focal_image_0, MAX=MAX(focal_image_0[n_y*0.1:n_y*0.8, n_z*0.1:n_z*0.8]), MIN=MIN(focal_image_0[n_y*0.1:n_y*0.8, n_z*0.1:n_z*0.8]))
    WRITE_PNG, focal_image_Phase_fn, BYTSCL(focal_image_0_Phase, MAX=MAX(focal_image_0_Phase[n_y*0.1:n_y*0.8, n_z*0.1:n_z*0.8]), MIN=MIN(focal_image_0_Phase[n_y*0.1:n_y*0.8, n_z*0.1:n_z*0.8]))
  ENDFOR
  DEVICE, DECOMPOSED=1
  
  ;@GJ, 2023/11/26, calculate the image
  A_DC = A_array[N_ELEMENTS(A_array)/2-1]
  result = MIN(ABS(A_DC_array-A_DC), A_DC_ind)
  A3_amp_array_0 = REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, N_ELEMENTS(A_DC_array)/2-1, *])
  A3_amp_array_T = REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, A_DC_ind, *])
  A3_amp_array_M = REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, A_DC_ind-4., *])
  A3_amp_array_P = REFORM(Harmonic_Amp_array[2, 0, 0, N_ELEMENTS(A_array)/2-1, A_DC_ind+4., *])
  focal_image_0 = DBLARR(N_ELEMENTS(H_Gy_array), N_ELEMENTS(H_Gz_array))
  focal_image_T = focal_image_0 * 0.
  focal_image_M = focal_image_0 * 0.
  focal_image_P = focal_image_0 * 0.
  FOR i=0, N_ELEMENTS(H_Gy_array)-1 DO BEGIN
    FOR j=0, N_ELEMENTS(H_Gz_array)-1 DO BEGIN
      H_G = SQRT(H_Gy_array[i]^2 + H_Gz_array[j]^2)
      focal_image_0[i,j] = INTERPOL(A3_amp_array_0, A_G_array, H_G)
      focal_image_T[i,j] = INTERPOL(A3_amp_array_T, A_G_array, H_G)
      focal_image_M[i,j] = INTERPOL(A3_amp_array_M, A_G_array, H_G)
      focal_image_P[i,j] = INTERPOL(A3_amp_array_P, A_G_array, H_G)
    ENDFOR
  ENDFOR
;  iimage, BYTSCL(focal_image_0), title='A3 Focal Spot Image (H_DC: 0 mT)'
;  iimage, BYTSCL(focal_image_T), title='A3 Focal Spot Image (H_DC: 10 mT)'
;  iimage, BYTSCL(focal_image_M), title='A3 Focal Spot Image (H_DC: 9 mT)'
;  iimage, BYTSCL(focal_image_P, MIN=MIN(focal_image_P[n_y/2-5:n_y/2+5, n_z/2-5:n_z/2+5])), title='A3 Focal Spot Image (H_DC: 11 mT)'
;  iimage, BYTSCL(focal_image_M-focal_image_P), title='A3 Focal Spot Image Subtraction (H_DC: 9-11 mT)'
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A3Spot_DC_0.png', BYTSCL(focal_image_0)
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A3Spot_DC_T.png', BYTSCL(focal_image_T)
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A3Spot_DC_M.png', BYTSCL(focal_image_M)
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A3Spot_DC_P.png', BYTSCL(focal_image_P, MIN=MIN(focal_image_P[n_y/2-5:n_y/2+5, n_z/2-5:n_z/2+5]))
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A3Spot_Subt.png', BYTSCL(focal_image_M-focal_image_P)

  ;@GJ, 2023/11/26, calculate the image A2
  A2_amp_array_0 = REFORM(Harmonic_Amp_array[1, 0, 0, N_ELEMENTS(A_array)/2-1, N_ELEMENTS(A_DC_array)/2-1, *])
  A2_amp_array_T = REFORM(Harmonic_Amp_array[1, 0, 0, N_ELEMENTS(A_array)/2-1, A_DC_ind, *])
  A2_amp_array_M = REFORM(Harmonic_Amp_array[1, 0, 0, N_ELEMENTS(A_array)/2-1, A_DC_ind-4., *])
  A2_amp_array_P = REFORM(Harmonic_Amp_array[1, 0, 0, N_ELEMENTS(A_array)/2-1, A_DC_ind+4., *])
  A2_focal_image_0 = DBLARR(N_ELEMENTS(H_Gy_array), N_ELEMENTS(H_Gz_array))
  A2_focal_image_T = A2_focal_image_0 * 0.
  A2_focal_image_M = A2_focal_image_0 * 0.
  A2_focal_image_P = A2_focal_image_0 * 0.
  FOR i=0, N_ELEMENTS(H_Gy_array)-1 DO BEGIN
    FOR j=0, N_ELEMENTS(H_Gz_array)-1 DO BEGIN
      H_G = SQRT(H_Gy_array[i]^2 + H_Gz_array[j]^2)
      A2_focal_image_0[i,j] = INTERPOL(A2_amp_array_0, A_G_array, H_G)
      A2_focal_image_T[i,j] = INTERPOL(A2_amp_array_T, A_G_array, H_G)
      A2_focal_image_M[i,j] = INTERPOL(A2_amp_array_M, A_G_array, H_G)
      A2_focal_image_P[i,j] = INTERPOL(A2_amp_array_P, A_G_array, H_G)
    ENDFOR
  ENDFOR
;  iimage, BYTSCL(A2_focal_image_0), title='A2 Focal Spot Image (H_DC: 0 mT)'
;  iimage, BYTSCL(A2_focal_image_T), title='A2 Focal Spot Image (H_DC: 10 mT)'
;  iimage, BYTSCL(A2_focal_image_M), title='A2 Focal Spot Image (H_DC: 9 mT)'
;  iimage, BYTSCL(A2_focal_image_P, MIN=MIN(A2_focal_image_P[n_y/2-5:n_y/2+5, n_z/2-5:n_z/2+5])), title='A2 Focal Spot Image (H_DC: 11 mT)'
;  iimage, BYTSCL(A2_focal_image_M-A2_focal_image_P), title='A2 Focal Spot Image Subtraction (H_DC: 9-11 mT)'
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A2Spot_DC_0.png', BYTSCL(A2_focal_image_0)
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A2Spot_DC_T.png', BYTSCL(A2_focal_image_T)
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A2Spot_DC_M.png', BYTSCL(A2_focal_image_M)
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A2Spot_DC_P.png', BYTSCL(A2_focal_image_P, MIN=MIN(A2_focal_image_P[n_y/2-5:n_y/2+5, n_z/2-5:n_z/2+5]))
  WRITE_PNG, dir_name+current_time+'_AC'+STRING(A_array[N_ELEMENTS(A_array)/2-1], format='(f04.1)')+'mT_A2Spot_Subt.png', BYTSCL(A2_focal_image_M-A2_focal_image_P)
  
  ;@GJ, 2023/12/16, delete the variable
  Harmonic_Amp_array = 0.
  Harmonic_Real_array = 0.
  Harmonic_Imaginary_array = 0.
  Harmonic_Phase_array = 0.
  A3_Amp_array = 0.
  A3_Phase_array = 0.
  A3_Real_array = 0.
  A3_Imaginary_array = 0.
  
  DC_zero_ind = WHERE(ABS(A_DC_Array) EQ 0)
  cgloadCT, 2, RGB_Table=blugrn_palette
  ;isurface, REFORM(donut_Harmonic_radius[2, 0, 0, *, DC_zero_ind:*]), A_array, A_DC_array[DC_zero_ind:*], xtitle='H_AC [mT]', ytitle='H_DC [mT]', ztitle='Donut Radius [a.u.]', title='Donut Radius', COLOR='yellow'
  ;isurface, REFORM(non_donut_Harmonic_FWHM[2, 0, 0, *, DC_zero_ind:*]), A_array, A_DC_array[DC_zero_ind:*], xtitle='H_AC [mT]', ytitle='H_DC [mT]', ztitle='Non-Donut FWHM [a.u.]', title='Non-Donut FWHM', COLOR='light blue'
  cgsurface, REFORM(donut_Harmonic_radius[2, 0, 0, *, DC_zero_ind:*]), A_array, A_DC_array[DC_zero_ind:*], CTable=2, xtitle='H_AC [mT]', ytitle='H_DC [mT]', ztitle='Donut Radius [mT]', title='Donut Radius', /Shaded, /Elevation_Shading
  cgcolorbar, Palette=blugrn_palette, window=0, Range=[Min(donut_Harmonic_radius[2, 0, 0, *, DC_zero_ind:*]),Max(donut_Harmonic_radius[2, 0, 0, *, DC_zero_ind:*])], title='Donut Radius [mT]'
  cgsurface, REFORM(donut_Harmonic_thick[2, 0, 0, *, DC_zero_ind:*]), A_array, A_DC_array[DC_zero_ind:*], CTable=2, xtitle='H_AC [mT]', ytitle='H_DC [mT]', ztitle='Donut Thickness [a.u.]', title='Donut Thickness', /Shaded, /Elevation_Shading
  cgcolorbar, Palette=blugrn_palette, window=1, Range=[Min(donut_Harmonic_thick[2, 0, 0, *, DC_zero_ind:*]),Max(donut_Harmonic_thick[2, 0, 0, *, DC_zero_ind:*])], title='Donut Thickness [a.u.]'
  cgsurface, REFORM(non_donut_Harmonic_FWHM[2, 0, 0, *, DC_zero_ind:*]), A_array, A_DC_array[DC_zero_ind:*], CTable=2, xtitle='H_AC [mT]', ytitle='H_DC [mT]', ztitle='Bright Focal FWHM [mT]', title='Bright Focal FWHM', /Shaded, /Elevation_Shading
  cgcolorbar, Palette=blugrn_palette, window=2, Range=[Min(non_donut_Harmonic_FWHM[2, 0, 0, *, DC_zero_ind:*]),Max(non_donut_Harmonic_FWHM[2, 0, 0, *, DC_zero_ind:*])], title='Bright Focal FWHM [mT]'
  cgsurface, REFORM(non_donut_Harmonic_Height[2, 0, 0, *, DC_zero_ind:*]), A_array, A_DC_array[DC_zero_ind:*], CTable=2, xtitle='H_AC [mT]', ytitle='H_DC [mT]', ztitle='Bright Focal Max Signal [a.u.]', title='Bright Focal Max Signal', /Shaded, /Elevation_Shading
  cgcolorbar, Palette=blugrn_palette, window=3, Range=[Min(non_donut_Harmonic_Height[2, 0, 0, *, DC_zero_ind:*]),Max(non_donut_Harmonic_Height[2, 0, 0, *, DC_zero_ind:*])], title='Bright Focal Max Signal [a.u.]'

END



;@GJ, 2023/10/22, RESI based on donut and critical point (NC)
PRO RESI_donut_Orth, A, A_DC, Harmonic_3rd_Amp_array, PSF_Harmonic_3rd_Amp_array, rec_A3_Amp, Harmonic_3rd_Phase_array, Harmonic_2nd_Amp_array, Harmonic_2nd_Phase_array, Donut, Donut2, phantom_2D_array, PSF_phantom, z_axis, filename
  
  ;initial setting
  ;A = 4.;20.
;  A_DC = 8.;10.;21.;10.;(*0.1mT)
  particle_size = 30.;5.;30
  f = 1
  tau = 20.;0.;0.01;30.;1.;30.;3.;20.;3.;20  ;(us)
  flat_portion = 50 ;
  gradient_field = 20. ; (*0.1mT/mm)
  gradient_field_y = gradient_field
  gradient_field_z = gradient_field
  signal_ratio = 0.;
  RESI_ratio = 1.;
  STD_reso = 1.;
  RESI_reso = 1.;
  Dpp_reso = 1.;
  noise_level = 0.01;
  N_points = 5.;2.;1.;4.;2.;3.;2.;2.;2.;2.;3.;3.;3.;2.;2.;3.;3.;1.;2.;3.;2.;2.;2.;0.;0.;3.;2.;1.;2.;2.;1.;2.;
  D_points = 5.;10.;10.;5;20.;10.;20.;20.;100.;20.;20.;5.;30.;10.; *0.1mm
  fov_phantom = 250.;
  n_angles = 4.;
  angle = 90.;
  shift_perc = 43.;100.;
  
  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  
  n_Time = 1000;5000.
  Time = findgen(n_Time)/1000.     ;(1  ms)
  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
  
  tau_ms = tau/1000.
  R = (1/tau_ms)*exp(-Time/tau_ms)
  kernel = Debye_Kernel_sin(tau_ms, Time[0:100])
  
  ;@GJ, 2023/9/8, plot the 2D PSF figure
  FOV = 2. * 20. / (gradient_field / 10.) ; in mm
  n_x = FLOOR(n_Time/(5.*2.*f)) ;this is corrected number
  print, 'FOV = ', STRING(FOV,format='(f6.2)')+' mm'
  ;@GJ, 2022/12/14, changing to offset
  n_y = n_x
  n_z = n_x
  resolution = [FOV/n_x, FOV/n_y, FOV/n_z]
  print, 'resolution: ', resolution
  y_axis = FINDGEN(n_y) * resolution[1]
  z_axis = FINDGEN(n_z) * resolution[2]

  ;@GJ, 2023/10/11, design the phantom
  phantom_1D_array = DBLARR(n_y*4.) * 0.
  N_D_points = FLOOR(D_points*0.1/resolution[1])
  Npoints_ind = 0.
  FOR i=0, n_y*4.-1 DO BEGIN
    IF (i MOD N_D_points) EQ 0 AND Npoints_ind LT N_points THEN  BEGIN
      phantom_1D_array[i] = 1.;i+1;1.
      Npoints_ind++
    ENDIF
  ENDFOR
  phantom_1D_array = SHIFT(phantom_1D_array, n_y/2.-(N_D_points*(N_points-1)/2.))
;  iplot, z_axis, phantom_1D_array[0:n_y-1], xtitle='z axis [mm]', title='1D phantom'
  phantom_2D_array_temp = DBLARR(n_y, n_z) * 0.
;  FOR kl=0, n_y-1 DO phantom_2D_array_temp[kl, *] = phantom_1D_array[0:n_y-1]
  phantom_2D_array_temp[n_y/2, *] = phantom_1D_array[0:n_y-1]
  
  ;@GJ, 2023/11/18, PSF phantom definition
  PSF_phantom = DBLARR(n_y, n_z) * 0.
  PSF_phantom[n_y/2-1, n_z/2-1] = 1.
  
  ;@GJ, 2023/11/1, if filename is not ready, read a new phantom
  IF N_ELEMENTS(filename) NE 0 THEN BEGIN
    IF STRLEN(filename) GT 5 AND (FILE_INFO(filename)).exists EQ 1 AND N_ELEMENTS(phantom_2D_array) GT 5 THEN BEGIN
      phantom_1D_array[0:n_y-1] = phantom_2D_array[n_y/2, *]
    ENDIF ELSE BEGIN
      phantom_1D_array[0:n_y-1] = phantom_2D_array_temp[n_y/2, *]
      phantom_2D_array = phantom_2D_array_temp
    ENDELSE
  ENDIF ELSE BEGIN
    ;  iimage, phantom_2D_array, title='Phantom Image'
    phantom_dir='.\lib\icon_pics\phantom\'
    filename = dialog_pickfile(TITLE='Select Phantom Image File', FILTER=['*.jpg', '*.png'], /MUST_EXIST, PATH=phantom_dir)
    ;  IF STRLEN(filename) LT 5 THEN return
    ;  IF (FILE_INFO(filename)).exists EQ 0 THEN return
    ;filename = 'C:\D_drive\MPI_Tianjie\RotationalDriftMPI\MRA_ori.png'
    IF QUERY_PNG(filename) EQ 1 THEN BEGIN
      image_temp = read_png(filename)
      IF (SIZE(image_temp))[0] EQ 2 THEN BEGIN
        image_sim_temp = DOUBLE(REFORM(image_temp[*,*]))
        ;      IF STRPOS(filename, 'MRA') NE -1 THEN image_sim_ori = DOUBLE(image_sim_temp GT 78) * 255.
        image_sim_ori = image_sim_temp
        IF STRPOS(filename, 'MRA') NE -1 THEN image_sim_ori[WHERE(image_sim_temp LT 78)] = 0.
      ENDIF ELSE BEGIN
        image_sim_temp = DOUBLE(REFORM(image_temp[0,*,*]))
        image_sim_ori = MAX(image_sim_temp) - image_sim_temp
        IF STRPOS(filename, 'MRA') NE -1 THEN image_sim_ori = DOUBLE(image_sim_ori LT 192) * 255.
      ENDELSE
      rows_S = N_ELEMENTS(image_sim_ori[0,*])
      cols_S = N_ELEMENTS(image_sim_ori[*,0])
      FOV_sim = fov_phantom; mm
      pixelSp_sim = [FOV_sim/DOUBLE(cols_S), FOV_sim/DOUBLE(rows_S)]
      print, 'pixelSp_sim: ', pixelSp_sim
      phantom_2D_array = BYTSCL(image_sim_ori[cols_S/2-n_y/2:cols_S/2+n_y/2-1, rows_S/2-n_z/2:rows_S/2+n_z/2-1])
      min_p = 4.5*MEAN(phantom_2D_array)
      max_p = MAX(phantom_2D_array)
      phantom_2D_array = DOUBLE(BYTSCL(phantom_2D_array, MIN=min_p, MAX=max_p))
      phantom_1D_array[0:n_y-1] = phantom_2D_array[n_y/2, *]
    ENDIF ELSE BEGIN
      phantom_1D_array[0:n_y-1] = phantom_2D_array_temp[n_y/2, *]
      phantom_2D_array = phantom_2D_array_temp
    ENDELSE
  ENDELSE
  dir_name = '.\RESI_donut_pics\AC'+STRING(A, format='(f5.1)')+'mT_'
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  ;@GJ, 2023/11/1, save the original image
  iimage, REVERSE(BYTSCL(phantom_2D_array),2), title='Original image'
  WRITE_PNG, dir_name+current_time+'_Original_image.png', REVERSE(BYTSCL(phantom_2D_array),2)
  
  H_y_array_L = (FINDGEN(n_y*2.)-FLOOR(n_y)) * resolution[1] * (gradient_field_y/10.)
  H_z_array_L = (FINDGEN(n_z*2.)-FLOOR(n_z)) * resolution[2] * (gradient_field_z/10.)
  signal_x_array = DBLARR(n_y, n_z, n_Time-1) * 0.
  Harmonic_3rd_Amp_array = DBLARR(n_y, n_z) * 0.
  Harmonic_3rd_Phase_array = DBLARR(n_y, n_z) * 0.
  Harmonic_5th_Amp_array = DBLARR(n_y, n_z) * 0.
  Harmonic_5th_Phase_array = DBLARR(n_y, n_z) * 0.
  Harmonic_7th_Amp_array = DBLARR(n_y, n_z) * 0.
  Harmonic_7th_Phase_array = DBLARR(n_y, n_z) * 0.
  Harmonic_2nd_Amp_array = DBLARR(n_y, n_z) * 0.
  Harmonic_2nd_Phase_array = DBLARR(n_y, n_z) * 0.
  Donut = DBLARR(n_y, n_z) * 0.
  Donut2 = DBLARR(n_y, n_z) * 0.
  PSF_signal_x_array = DBLARR(n_y, n_z, n_Time-1) * 0.
  PSF_Harmonic_3rd_Amp_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_3rd_Phase_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_5th_Amp_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_5th_Phase_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_7th_Amp_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_7th_Phase_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_2nd_Amp_array = DBLARR(n_y, n_z) * 0.
  PSF_Harmonic_2nd_Phase_array = DBLARR(n_y, n_z) * 0.

;  FOR m=0, n_y-1 DO BEGIN
  FOR m=n_y/2, n_y/2 DO BEGIN
    FOR n=0, n_z-1 DO BEGIN
      H_y_array = H_y_array_L[m:m+n_y-1]
      H_z_array = H_z_array_L[n:n+n_y-1]
      M_x = DBLARR(n_Time, n_y, n_z) * 0.
      PSF_M_x = DBLARR(n_Time, n_y, n_z) * 0.
      FOR i=0, n_y-1 DO BEGIN
        FOR j=0, n_z-1 DO BEGIN
          IF phantom_2D_array[i, j] GT 0 OR PSF_phantom[i, j] GT 0 THEN BEGIN
            print, 'm, n=', m, ',', n
            ;@GJ, 2023/10/20, adding the DC field
            H_x = -A*cos(2.*!PI*f*Time) + A_DC
            H_y = H_y_array[i]
            H_z = H_z_array[j]; + 14.
            H_yz = SQRT(H_y^2 + H_z^2); + 14.;12.;14.;4.;7.;9.; offset
            H_total = SIGNUM(H_x) * SQRT(H_x^2 + H_yz^2)
            print, 'Hy, Hz', H_y, H_z, ' mT'

            ;calculate M
            M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
            index_nan = where(~FINITE(M_total), count)
            IF count GT 0 THEN BEGIN
              RESULT = INTERPOL(M_total, Time, Time[index_nan], /NAN)
              M_total[index_nan] = RESULT
            ENDIF

            ;do vector splitting
            ;@GJ, 2023/9/23, simple calculation with the signal intensity
            M_x[*,i,j] = ABS(M_total) * H_x/ABS(H_total) * phantom_2D_array[i, j]
            PSF_M_x[*,i,j] = ABS(M_total) * H_x/ABS(H_total) * PSF_phantom[i, j]

            signal_x = -(M_x[1:*, i, j] - M_x[0:n_Time-2, i, j])/(Time[1]-Time[0])
            ;@GJ, 2022/12/14, removing the NaN points
            index_nan_x = where(~FINITE(signal_x), count_x)
            IF count_x GT 0 THEN BEGIN
              RESULT_x = INTERPOL(signal_x, Time, Time[index_nan_x], /NAN)
              signal_x[index_nan_x] = RESULT_x
            ENDIF
            signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
            signal_new_x = signal_new_x[n_elements(kernel)/2 : n_elements(signal_x) + n_elements(kernel)/2]
            
            ;@GJ, 2023/11/18, add noise to the signal
            noise_level = ABS(MAX(signal_new_x))/1000.
            seed = !NULL
            noise_array = noise_level*RANDOMU(seed, N_ELEMENTS(signal_new_x))
;            signal_new_x[n_Time/2:*] = 0.
            signal_x_array[m, n, *] += signal_new_x + noise_array
            
            PSF_signal_x = -(PSF_M_x[1:*, i, j] - PSF_M_x[0:n_Time-2, i, j])/(Time[1]-Time[0])
            ;@GJ, 2022/12/14, removing the NaN points
            PSF_index_nan_x = where(~FINITE(PSF_signal_x), PSF_count_x)
            IF PSF_count_x GT 0 THEN BEGIN
              PSF_RESULT_x = INTERPOL(PSF_signal_x, Time, Time[PSF_index_nan_x], /NAN)
              PSF_signal_x[index_nan_x] = PSF_RESULT_x
            ENDIF
            PSF_signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, PSF_signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
            PSF_signal_new_x = PSF_signal_new_x[n_elements(kernel)/2 : n_elements(PSF_signal_x) + n_elements(kernel)/2]
            ;            signal_new_x[n_Time/2:*] = 0.
            PSF_signal_x_array[m, n, *] += PSF_signal_new_x

          ENDIF
        ENDFOR
      ENDFOR
      
      t_ele = N_Elements(signal_x_array[m, n, *])
      delta_t = (Time[1]-Time[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
      X = (FINDGEN((t_ele - 1)/2) + 1)
      is_N_even = (t_ele MOD 2) EQ 0
      if (is_N_even) then $
        frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
      else $
        frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
      base_freq = MIN(ABS(frequency-1.0), base_ind)

      signal_FFT_x = FFT(signal_x_array[m, n, *])
      ;@GJ, 2023/10/20, get the harmonic amp and phase
      H_ind = 3.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      Harmonic_3rd_Amp_array[m, n] = ABS(signal_FFT_x[base_ind*H_ind])
      Harmonic_3rd_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*H_ind]), Real_part(signal_FFT_x[base_ind*H_ind]), /phase)
;      IF Harmonic_3rd_Phase_array[m, n] LT 0 THEN BEGIN
;        Harmonic_3rd_Amp_array[m, n] *= -1.
;;        Harmonic_3rd_Phase_array[m, n] *= -1.
;      ENDIF
      H_ind = 5.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      Harmonic_5th_Amp_array[m, n] = ABS(signal_FFT_x[base_ind*H_ind])
      Harmonic_5th_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*H_ind]), Real_part(signal_FFT_x[base_ind*H_ind]), /phase)
;      IF Harmonic_5th_Phase_array[m, n] LT 0 THEN BEGIN
;        Harmonic_5th_Amp_array[m, n] *= -1.
;        ;        Harmonic_3rd_Phase_array[m, n] *= -1.
;      ENDIF
      H_ind = 7.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      Harmonic_7th_Amp_array[m, n] = ABS(signal_FFT_x[base_ind*H_ind])
      Harmonic_7th_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*H_ind]), Real_part(signal_FFT_x[base_ind*H_ind]), /phase)
;      IF Harmonic_7th_Phase_array[m, n] LT 0 THEN BEGIN
;        Harmonic_7th_Amp_array[m, n] *= -1.
;        ;        Harmonic_3rd_Phase_array[m, n] *= -1.
;      ENDIF
      Harmonic_2nd_Amp_array[m, n] = ABS(signal_FFT_x[base_ind*2.])
      Harmonic_2nd_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*2.]), Real_part(signal_FFT_x[base_ind*2.]), /phase)
      IF Harmonic_2nd_Phase_array[m, n] LT 0 THEN BEGIN
        Harmonic_2nd_Amp_array[m, n] *= -1.
        Harmonic_2nd_Phase_array[m, n] *= -1.
      ENDIF
      
      PSF_signal_FFT_x = FFT(PSF_signal_x_array[m, n, *])
      ;@GJ, 2023/10/20, get the harmonic amp and phase
      H_ind = 3.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      PSF_Harmonic_3rd_Amp_array[m, n] = ABS(PSF_signal_FFT_x[base_ind*H_ind])
      PSF_Harmonic_3rd_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(PSF_signal_FFT_x[base_ind*H_ind]), Real_part(PSF_signal_FFT_x[base_ind*H_ind]), /phase)
;      IF PSF_Harmonic_3rd_Phase_array[m, n] LT 0 THEN BEGIN
;        PSF_Harmonic_3rd_Amp_array[m, n] *= -1.
;        PSF_Harmonic_3rd_Phase_array[m, n] *= -1.
;      ENDIF
      H_ind = 5.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      PSF_Harmonic_5th_Amp_array[m, n] = ABS(PSF_signal_FFT_x[base_ind*H_ind])
      PSF_Harmonic_5th_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(PSF_signal_FFT_x[base_ind*H_ind]), Real_part(PSF_signal_FFT_x[base_ind*H_ind]), /phase)
;      IF PSF_Harmonic_5th_Phase_array[m, n] LT 0 THEN BEGIN
;        PSF_Harmonic_5th_Amp_array[m, n] *= -1.
;        PSF_Harmonic_5th_Phase_array[m, n] *= -1.
;      ENDIF
      H_ind = 7.;7.;3.;7.;9.;7.;3.;7.;5.;3.;5.;7.
      PSF_Harmonic_7th_Amp_array[m, n] = ABS(PSF_signal_FFT_x[base_ind*H_ind])
      PSF_Harmonic_7th_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(PSF_signal_FFT_x[base_ind*H_ind]), Real_part(PSF_signal_FFT_x[base_ind*H_ind]), /phase)
;      IF PSF_Harmonic_7th_Phase_array[m, n] LT 0 THEN BEGIN
;        PSF_Harmonic_7th_Amp_array[m, n] *= -1.
;        PSF_Harmonic_7th_Phase_array[m, n] *= -1.
;      ENDIF
      PSF_Harmonic_2nd_Amp_array[m, n] = ABS(PSF_signal_FFT_x[base_ind*2.])
      PSF_Harmonic_2nd_Phase_array[m, n] = 180./!PI*ATAN(Imaginary(PSF_signal_FFT_x[base_ind*2.]), Real_part(PSF_signal_FFT_x[base_ind*2.]), /phase)
      IF PSF_Harmonic_2nd_Phase_array[m, n] LT 0 THEN BEGIN
        PSF_Harmonic_2nd_Amp_array[m, n] *= -1.
        PSF_Harmonic_2nd_Phase_array[m, n] *= -1.
      ENDIF
;      print, 'm, n=', m, ',', n
    ENDFOR
    
    
    ;calculate the donut diameter and thickness
    donut[m, *] = Harmonic_3rd_Amp_array[m, *] * SIN(Harmonic_3rd_Phase_array[m, *])
    IF MEAN(donut[m, *]) LT 0 THEN donut[m, *] *= -1
    
    
    donut2[m, *] = Harmonic_2nd_Amp_array[m, *] * SIN(Harmonic_2nd_Phase_array[m, *])
    IF MEAN(donut2[m, *]) LT 0 THEN donut2[m, *] *= -1


    IF m EQ n_y/2 AND MAX(phantom_1D_array[0:n_y-1]) GT 0 THEN BEGIN
      IF OBJ_VALID(p) THEN p.close
      p = PLOT(z_axis, REVERSE(Donut[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
        XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Amp3Donut [A.U.]', xtitle='z axis [mm]', title='A3Donut (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
      p.name = ' Amp*sinPhi'
      p1 = PLOT(z_axis, phantom_1D_array[0:n_y-1]/MAX(phantom_1D_array[0:n_y-1])*0.5*MAX(Donut[m, *])-MIN(Donut[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
      p1.name = ' Phantom'
      l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
;      fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
;      IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
      dir_name = '.\RESI_donut_pics\AC'+STRING(A, format='(f5.1)')+'mT_DC'+STRING(A_DC, format='(f5.1)')+'mT_'
      temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
      ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
      current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
      output_filename_temp = dir_name+current_time+'_Harmonic_3rd_AmpDonut_curve.png'
      p.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
      ;close the plot
      p.close
      
;      IF OBJ_VALID(p) THEN p.close
;      p = PLOT(z_axis, REVERSE(Donut2[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
;        XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Amp2Donut', xtitle='z axis [mm]', title='A2Donut (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT')
;      p.name = ' Amp'
;      p1 = PLOT(z_axis, phantom_1D_array[0:n_y-1]/MAX(phantom_1D_array[0:n_y-1])*0.5*MAX(Donut2[m, *])-MIN(Donut2[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
;      p1.name = ' Phantom'
;      l = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
;      ;      fitplot_dir = FILE_DIRNAME(filename, /MARK_DIRECTORY)+'fitplot\'
;      ;      IF FILE_TEST(fitplot_dir, /directory) LT 1 THEN FILE_MKDIR, fitplot_dir
;      dir_name = '.\RESI_donut_pics\AC'+STRING(A, format='(f5.1)')+'mT_DC'+STRING(A_DC, format='(f5.1)')+'mT_'
;      temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
;      ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
;      current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
;;      output_filename_temp = dir_name+current_time+'_Harmonic_2nd_AmpDonut_curve.png'
;      p.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
      ;close the plot
      
      IF OBJ_VALID(p0) THEN p0.close
      p0 = PLOT(z_axis, REVERSE(Harmonic_2nd_Phase_array[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
        XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Phi2 (deg)', xtitle='z axis [mm]', title='Phi2 (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
      p0.name = ' Phase'
      p01 = PLOT(z_axis, phantom_1D_array[0:n_y-1]/MAX(phantom_1D_array[0:n_y-1])*2.*STDDEV(Harmonic_2nd_Phase_array[m, *])+MIN(Harmonic_2nd_Phase_array[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
      p01.name = ' Phantom'
      l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
      output_filename_temp = dir_name+current_time+'_Harmonic_2nd_Phase_curve.png'
      p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
      p0.close
      
      IF OBJ_VALID(p1) THEN p1.close
      p1 = PLOT(z_axis, REVERSE(Harmonic_3rd_Amp_array[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
        XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='A3 [A.U.]', xtitle='z axis [mm]', title='A3 (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
      p1.name = ' Amp'
      p11 = PLOT(z_axis, phantom_1D_array[0:n_y-1]/MAX(phantom_1D_array[0:n_y-1])*2.*STDDEV(Harmonic_3rd_Amp_array[m, *])+MIN(Harmonic_3rd_Amp_array[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
      p11.name = ' Phantom'
      l1 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
      output_filename_temp = dir_name+current_time+'_Harmonic_3rd_Amp_curve.png'
      p1.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
      
      IF OBJ_VALID(p0) THEN p0.close
      p0 = PLOT(z_axis, REVERSE(Harmonic_3rd_Phase_array[m, *],2), 'r-2', YSTYLE=3, FONT_SIZE=9,$; LAYOUT=[2,2,2], /CURRENT, $
        XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Phi3 [deg]', xtitle='z axis [mm]', title='Phi3 (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
      p0.name = ' Phase'
      p01 = PLOT(z_axis, phantom_1D_array[0:n_y-1]/MAX(phantom_1D_array[0:n_y-1])*2.*STDDEV(Harmonic_3rd_Phase_array[m, *])+MIN(Harmonic_3rd_Phase_array[m, *]), 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
      p01.name = ' Phantom'
      l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.4])
      output_filename_temp = dir_name+current_time+'_Harmonic_3rd_Phase_curve.png'
      p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
      p0.close
;      iplot, z_axis, Harmonic_3rd_Amp_array[m, *]/MAX(Harmonic_3rd_Amp_array[m, *]), color='red', xtitle='z axis [mm]', title='3rd Amp'
;      iplot, z_axis, Harmonic_3rd_Phase_array[m, *]/MAX(Harmonic_3rd_Phase_array[m, *]), color='blue', /overplot
;      iplot, z_axis, Harmonic_2nd_Amp_array[m, *], xtitle='z axis [mm]', title='2nd Amp'
    ENDIF

  ENDFOR
  
  iimage, BYTSCL(Harmonic_3rd_Amp_array, max=MAX(Harmonic_3rd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(Harmonic_3rd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='3rd Amp'
;  iimage, BYTSCL(Donut, max=MAX(Donut[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(Donut[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='3rd Amp Donut'
;  iimage, BYTSCL(Harmonic_3rd_Phase_array, max=MAX(Harmonic_3rd_Phase_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), min=MIN(Harmonic_3rd_Phase_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='3rd Phase'
  ;save as png file
  dir_name = '.\RESI_donut_pics\AC'+STRING(A, format='(f5.1)')+'mT_DC'+STRING(A_DC, format='(f5.1)')+'mT_'
  temp_current_time = STRING(systime(/julian), format='(c(CYI,CMOI,CDI,CHI,CMI,CSI))')
  ;  current_time = STRCOMPRESS(temp_current_time, /REMOVE_ALL)
  current_time = STRJOIN(STRSPLIT(temp_current_time, /EXTRACT), '0')
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_AmpDonut.png', BYTSCL(Donut, max=MAX(Donut[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(Donut[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_Amp.png', BYTSCL(Harmonic_3rd_Amp_array, max=MAX(Harmonic_3rd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(Harmonic_3rd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_PSF_Harmonic_3rd_Amp.png', BYTSCL(PSF_Harmonic_3rd_Amp_array, max=MAX(PSF_Harmonic_3rd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_Harmonic_3rd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_Harmonic_3rd_Phase.png', BYTSCL(Harmonic_3rd_Phase_array, max=MAX(Harmonic_3rd_Phase_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), min=MIN(Harmonic_3rd_Phase_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_AmpDonut.png', BYTSCL(Donut2, max=MAX(Donut2[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(Donut2[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_Amp.png', BYTSCL(Harmonic_2nd_Amp_array, max=MAX(Harmonic_2nd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), min=MIN(Harmonic_2nd_Amp_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_Harmonic_2nd_Phase.png', BYTSCL(Harmonic_2nd_Phase_array, max=MAX(Harmonic_2nd_Phase_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), min=MIN(Harmonic_2nd_Phase_array[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
;  iimage, Harmonic_2nd_Phase_array, title='2nd Phase'
  
  z_ele = n_z
  delta_z = resolution[2] ; us; samplying rate: 1MHz
  Z = (FINDGEN((z_ele - 1)/2) + 1)
  is_N_even = (z_ele MOD 2) EQ 0
  if (is_N_even) then $
    frequency_z = [0.0, Z, z_ele/2, -z_ele/2 + Z]/(z_ele*delta_z) $
  else $
    frequency_z = [0.0, Z, -(z_ele/2 + 1) + Z]/(z_ele*delta_z)
  Phantom_1D_FFT = ABS(FFT(phantom_1D_array[0:n_y-1]))
  iplot, frequency_z[SORT(frequency_z)], Phantom_1D_FFT, xrange=[0, MAX(frequency_z)]
;  iplot, Phantom_1D_FFT,
  PSF_A3_Amp_FFT = ABS(FFT(PSF_Harmonic_3rd_Amp_array))
  A3_Amp_FFT = FFT(Harmonic_3rd_Amp_array)
  A3_Amp_FFT_div = A3_Amp_FFT/(PSF_A3_Amp_FFT);+MAX(A3_Amp_FFT)*0.01)
  rec_A3_Amp = FFT(A3_Amp_FFT_div, /INVERSE)
  p0 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(PSF_A3_Amp_FFT[n_y/2, *]))[SORT(frequency_z)], 'r-2', xrange=[0, MAX(frequency_z)], YSTYLE=3, FONT_SIZE=9, yrange=[0, 30],$; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Contrast', xtitle='Spatial Frequancy [lp/mm]', title='3rd FFT (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
  p0.name = ' PSF MTF'
  p01 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A3_Amp_FFT[n_y/2, *])))[SORT(frequency_z)], 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p01.name = ' Phantom FFT'
  p02 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A3_Amp_FFT_div[n_y/2, *])))[SORT(frequency_z)], 'g-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p02.name = ' Phantom Rec'
  l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  output_filename_temp = dir_name+current_time+'_Harmonic_3rd_Amp_FFT_curves.png'
  p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  iimage, BYTSCL(A3_Amp_FFT, max=MAX(A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='Phantom 3rd Amp (FFT)'
  iimage, BYTSCL(PSF_A3_Amp_FFT, max=MAX(PSF_A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='PSF 3rd Amp (FFT)'
  iimage, BYTSCL(rec_A3_Amp, max=MAX(rec_A3_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(rec_A3_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='Rec Phantom 3rd Amp'
  WRITE_PNG, dir_name+current_time+'_A3_Amp_phantom_FFT.png', BYTSCL(A3_Amp_FFT, max=MAX(A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_A3_Amp_PSF_FFT.png', BYTSCL(PSF_A3_Amp_FFT, max=MAX(PSF_A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_A3_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_A3_Amp_phantom_rec.png', BYTSCL(rec_A3_Amp, max=MAX(rec_A3_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(rec_A3_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))


  PSF_A5_Amp_FFT = ABS(FFT(PSF_Harmonic_5th_Amp_array))
  A5_Amp_FFT = FFT(Harmonic_5th_Amp_array)
  A5_Amp_FFT_div = A5_Amp_FFT/(PSF_A5_Amp_FFT);+MAX(A5_Amp_FFT)*0.01)
  rec_A5_Amp = FFT(A5_Amp_FFT_div, /INVERSE)
  p0 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(PSF_A5_Amp_FFT[n_y/2, *]))[SORT(frequency_z)], 'r-2', xrange=[0, MAX(frequency_z)], YSTYLE=3, FONT_SIZE=9, yrange=[0, 30], $; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Contrast', xtitle='Spatial Frequancy [lp/mm]', title='5th FFT (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
  p0.name = ' PSF MTF'
  p01 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A5_Amp_FFT[n_y/2, *])))[SORT(frequency_z)], 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p01.name = ' Phantom FFT'
  p02 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A5_Amp_FFT_div[n_y/2, *])))[SORT(frequency_z)], 'g-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p02.name = ' Phantom Rec'
  l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  output_filename_temp = dir_name+current_time+'_Harmonic_5th_Amp_FFT_curves.png'
  p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  iimage, BYTSCL(A5_Amp_FFT, max=MAX(A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='Phantom 5th Amp (FFT)'
  iimage, BYTSCL(PSF_A5_Amp_FFT, max=MAX(PSF_A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='PSF 5th Amp (FFT)'
  iimage, BYTSCL(rec_A5_Amp, max=MAX(rec_A5_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(rec_A5_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='Rec Phantom 5th Amp'
  WRITE_PNG, dir_name+current_time+'_A5_Amp_phantom_FFT.png', BYTSCL(A5_Amp_FFT, max=MAX(A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_A5_Amp_PSF_FFT.png', BYTSCL(PSF_A5_Amp_FFT, max=MAX(PSF_A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_A5_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_A5_Amp_phantom_rec.png', BYTSCL(rec_A5_Amp, max=MAX(rec_A5_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(rec_A5_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))


  PSF_A7_Amp_FFT = ABS(FFT(PSF_Harmonic_7th_Amp_array))
  A7_Amp_FFT = FFT(Harmonic_7th_Amp_array)
  A7_Amp_FFT_div = A7_Amp_FFT/(PSF_A7_Amp_FFT);+MAX(A7_Amp_FFT)*0.01)
  rec_A7_Amp = FFT(A7_Amp_FFT_div, /INVERSE)
  p0 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(PSF_A7_Amp_FFT[n_y/2, *]))[SORT(frequency_z)], 'r-2', xrange=[0, MAX(frequency_z)], YSTYLE=3, FONT_SIZE=9, yrange=[0, 30], $; LAYOUT=[2,2,2], /CURRENT, $
    XMINOR=2, YMINOR=2, /SYM_FILLED, DIM=[400, 300], ytitle='Contrast', xtitle='Spatial Frequancy [lp/mm]', title='7th FFT (AC'+STRING(A, format='(f5.1)')+'mT, DC'+STRING(A_DC, format='(f5.1)')+'mT)')
  p0.name = ' PSF MTF'
  p01 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A7_Amp_FFT[n_y/2, *])))[SORT(frequency_z)], 'b-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p01.name = ' Phantom FFT'
  p02 = PLOT(frequency_z[SORT(frequency_z)], (REFORM(ABS(A7_Amp_FFT_div[n_y/2, *])))[SORT(frequency_z)], 'g-2', /OVERPLOT, /SYM_FILLED, SYM_SIZE=0.75)
  p02.name = ' Phantom Rec'
  l0 = LEGEND(SHADOW=0,LINESTYLE='none', POSITION=[0.85, 0.8])
  output_filename_temp = dir_name+current_time+'_Harmonic_7th_Amp_FFT_curves.png'
  p0.save, output_filename_temp, BORDER=10, RESOLUTION=300, /TRANSPARENT
  iimage, BYTSCL(A7_Amp_FFT, max=MAX(A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='Phantom 7th Amp (FFT)'
  iimage, BYTSCL(PSF_A7_Amp_FFT, max=MAX(PSF_A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='PSF 7th Amp (FFT)'
  iimage, BYTSCL(rec_A7_Amp, max=MAX(rec_A7_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(rec_A7_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1])), title='Rec Phantom 7th Amp'
  WRITE_PNG, dir_name+current_time+'_A7_Amp_phantom_FFT.png', BYTSCL(A7_Amp_FFT, max=MAX(A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_A7_Amp_PSF_FFT.png', BYTSCL(PSF_A7_Amp_FFT, max=MAX(PSF_A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(PSF_A7_Amp_FFT[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  WRITE_PNG, dir_name+current_time+'_A7_Amp_phantom_rec.png', BYTSCL(rec_A7_Amp, max=MAX(rec_A7_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]), MIN=MIN(rec_A7_Amp[n_y/4:n_y*3/4-1, n_z/4:n_z*3/4-1]))
  
  iplot, REFORM(PSF_Harmonic_3rd_Amp_array[n_y/2, *]), title='3rd PSF'
  iplot, REFORM(PSF_Harmonic_5th_Amp_array[n_y/2, *]), title='5th PSF'
  iplot, REFORM(PSF_Harmonic_7th_Amp_array[n_y/2, *]), title='7th PSF'

;  print, 'test'
    
END

;@GJ, 2023/12/13, calculate the theoretical
;usage:
;   ;@GJ, 2024/7/8, n=2, plot the odd
;n=2
;chebyshev_MPI, n, z_array, orth_curve
PRO chebyshev_MPI, n, z_array, orth_curve

  n_x = 512
  n_z = 512
  Gz_A = FINDGEN(n_z)/128-2.
  kernel_Gz_A = Gz_A * 0.
;  n=3.
  FOR k=0, n_z-1 DO IF ABS(Gz_A[k]) LT 1.0 THEN kernel_Gz_A[k] = SIN(n*ACOS(Gz_A[k]))/SIN(ACOS(Gz_A[k]))*SQRT(1-Gz_A[k]^2)
;  iplot, kernel_Gz_A
  
  Sz_array = DBLARR(n_x, n_z) * 0.
  Sz_array_kernel = DBLARR(n_x, n_z) * 0.
  Sz_array_xzNonCOnv = DBLARR(n_x, n_z) * 0.
  x_array = FINDGEN(n_x)-255.
  z_array = FINDGEN(n_z)-255.
  FOR i=0, n_x-1 DO BEGIN
    base_array = (x_array[i])^2/((x_array[i])^2 + z_array^2)^1.5
    Sz_array[i, *] = CONVOL(base_array, kernel_Gz_A, /center, /EDGE_ZERO, /NORMALIZE)
    Sz_array_kernel[i, *] = kernel_Gz_A
    Sz_array_xzNonCOnv[i, *] = base_array
  ENDFOR
  
  FOR k=0, n_z-1 DO Sz_array[255, k] = (Sz_array[254, k] + Sz_array[256, k])/2.
  
;  print, 'test'
 ; iimage, TRANSPOSE(BYTSCL(Sz_array)), title='Sz'
 ; iimage, TRANSPOSE(BYTSCL(Sz_array_xzNonCOnv, max=8)), title='Sz xz non-convol'
 ; iimage, TRANSPOSE(BYTSCL(Sz_array_kernel)), title='Sz kernel'
  
  dir_name = '.\RESI_donut_pics\'
  ;@GJ, 2023/11/1, save the original image
  WRITE_PNG, dir_name+STRING(n, format='(i02)')+'_Sz.png', TRANSPOSE(BYTSCL(Sz_array))
  WRITE_PNG, dir_name+STRING(n, format='(i02)')+'_base.png', BYTSCL(TRANSPOSE((Sz_array_xzNonCOnv)), max=0.02)
  WRITE_PNG, dir_name+STRING(n, format='(i02)')+'_chebyshev.png', TRANSPOSE(BYTSCL(Sz_array_kernel))
  
  orth_curve = REFORM(ABS(Sz_array[*, 255]))
;  iplot, orth_curve, title='orth'

END

;@GJ, 2023/12/15, modeling the excitation with offset
PRO orth_curve_set_compare
  flag=0
  FOR n=1,7 DO BEGIN
    ;@GJ, 2024/7/8, n=2, plot the odd
    ;n=2
    chebyshev_MPI, n, z_array, orth_curve
    nterms = 6.;4.
    orth_curve_fit = GAUSSFIT(z_array, orth_curve, coeff, NTERMS=nterms)
    print, 'n=', n
    IF flag EQ 0 THEN BEGIN
      iplot, z_array, orth_curve/MAX(orth_curve), yrange=[0., 1.], LINESTYLE=flag+1, thick=2, title='Harmonics and GaussFit'
      iplot, z_array, orth_curve_fit/MAX(orth_curve_fit), color='red', thick=1, /overplot
    ENDIF ELSE BEGIN
      iplot, z_array, orth_curve/MAX(orth_curve), LINESTYLE=flag+1, thick=2, /overplot
      iplot, z_array, orth_curve_fit/MAX(orth_curve_fit), color='red', thick=1, /overplot
    ENDELSE
    flag++
  ENDFOR

END



Function All_Tags, structure, rootname

  ; This is a function that recursively searches through
  ; a structure tree, finding ALL of the structure's field names.
  ; It returns a pointer to an array of pointers, each pointing
  ; to the names of structure fields.

  IF N_Elements(rootname) EQ 0 THEN rootname = '.' ELSE $
    rootname = StrUpCase(rootname) + '.'
  names = Tag_Names(structure)
  retValue = Ptr_New(rootname + names)

  ; If any of the fields are structures, report them too.

  FOR j=0,N_Elements(names)-1 DO BEGIN
    ok = Execute('s = Size(structure.' + names[j] + ')')
    IF s[s[0]+1] EQ 8 THEN BEGIN
      newrootname = rootname + names[j]
      theseNames = Call_Function('All_Tags', $
        structure.(j), newrootname)
      retValue = [[retValue],[theseNames]]
    ENDIF
  ENDFOR

  RETURN, retValue
END
;-------------------------------------------------------------------



FUNCTION Get_Tags, structure, rootname

  ; This function returns the names of all structure fields
  ; in the structure as a string array. The names are given
  ; as valid structure names from the root structure name,
  ; which can be passed in along with the structure itself.

  On_Error, 1

  ; Check parameters.

  CASE N_Params() OF

    0: BEGIN
      Message, 'Structure argument is required.'
    ENDCASE

    1: BEGIN
      rootname = ''
      s = Size(structure)
      IF s[s[0]+1] NE 8 THEN $
        Message, 'Structure argument is required.'
    ENDCASE

    2: BEGIN
      s = Size(structure)
      IF s[s[0]+1] NE 8 THEN $
        Message, 'Structure argument is required.'
      s = Size(rootname)
      IF s[s[0]+1] NE 7 THEN $
        Message, 'Root Name parameter must be a STRING'
    ENDCASE

  ENDCASE

  tags = All_Tags(structure, rootname)

  ; Extract and free the first pointer.

  retval = [*tags[0,0]]
  Ptr_Free, tags[0,0]

  ; Extract and free the the rest of the pointers.

  s = Size(tags)
  FOR j=1,s[2]-1 DO BEGIN
    retval = [retval, *tags[0,j]]
    Ptr_Free, tags[0,j]
  ENDFOR
  Ptr_Free, tags

  ; Return the structure names.

  RETURN, retval
END

;@GJ, 2024/10/23, analyze csv file converted by Yongchen
PRO LEO_simulation_data_analysis_Zhao_Zhiyuan_data
  
  B_offset_file = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\Dp20nm_Bac10mT\B_bias.csv'
  temp_B = READ_CSV(B_offset_file, HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
  offset_field_array = temp_B.(0)
  min=MIN(ABS(offset_field_array-12.), minind)
;  N_B_offset = minind
  N_B_offset = N_ELEMENTS(offset_field_array)
  
  ques = DIALOG_MESSAGE("Select a folder?", /question)
  IF STRCMP(ques, 'Yes') THEN BEGIN
    real_dir = dialog_pickfile(TITLE='Select a directory', /DIRECTORY, /MUST_EXIST, PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\')
    ;count the number of files
    real_files = FILE_SEARCH(real_dir, '*kHz.csv', count=nrfile)
    title_temp = STRSPLIT(real_dir, '\', /EXTRACT, COUNT=n_dir)
    title = title_temp[n_dir-1]
  ENDIF ELSE BEGIN
    real_file = dialog_pickfile(TITLE='Select a csv file', /MUST_EXIST, PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\', FILTER='*.csv')
    ;count the number of files
    real_files = [real_file];['C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\Dp25nm_Bac10mT\f20kHz.csv']
    nrfile = N_ELEMENTS(real_files)
    title = FILE_BASENAME(real_file, '.csv')
  ENDELSE
  
  FOR i_file = 0, nrfile-1 DO BEGIN
    ;  real_file = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\', /MUST_EXIST, TITLE="Real csv file", FILTER='*.csv')
    
    real_data = READ_CSV(real_files[i_file], HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
    n_tag = N_TAGS(real_data)
    n_datapoints = N_ELEMENTS(real_data.(0))

    khz_pos = STRPOS(FILE_BASENAME(real_files[i_file], '.csv'), 'kHz')
    f_pos = STRPOS(FILE_BASENAME(real_files[i_file], '.csv'), 'f')
    freq_base = DOUBLE(STRMID(FILE_BASENAME(real_files[i_file], '.csv'), f_pos+1, khz_pos-f_pos-1));20. ;kHz
    print, real_files[i_file]
    print, 'base freq: ', freq_base
    
    start_ind = 1.;1. / freq_base / 2. / ((real_data.(1))[1000] - (real_data.(1))[0]) + 1.;0.;n_datapoints - 5000. * 8. - 2. - 2500. ;starting from the middle
    t_array = DBLARR(n_datapoints - start_ind)
    t_array = ((real_data.(1))[start_ind:*] - (real_data.(1))[start_ind]) * 1000. ;ms

    T = ((real_data.(1))[start_ind+1000] - (real_data.(1))[start_ind]) ;ms
    N = N_ELEMENTS(t_array)
    X = (FINDGEN((N - 1)/2) + 1)
    is_N_even = (N MOD 2) EQ 0
    if (is_N_even) then $
      freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
    else $
      freq = [0.0, X, -(N/2 + 1) + X]/(N*T)
    freq_delta = 1./(N*T)

    line_abs = DBLARR(N_B_offset, N/2.) * 0.
    line_real = DBLARR(N_B_offset, N/2.) * 0.
    line_imag = DBLARR(N_B_offset, N/2.) * 0.
    line_phase = DBLARR(N_B_offset, N/2.) * 0.

    start_col = 2
    FOR i=start_col, N_B_offset+start_col-1 DO BEGIN
      temp_signal = -((real_data.(i))[start_ind:*] - (real_data.(i))[(start_ind-1):(n_datapoints-2)])/T
      temp_signal_fft = FFT(temp_signal)
      line_abs[i-start_col, *] = ABS(temp_signal_fft[0:N/2.-1])
      line_real[i-start_col, *] = Real_part(temp_signal_fft[0:N/2.-1])
      line_imag[i-start_col, *] = Imaginary(temp_signal_fft[0:N/2.-1])
      line_phase[i-start_col, *] = 180./!PI*ATAN(temp_signal_fft[0:N/2.-1], /phase)
      IF nrfile EQ 1 THEN BEGIN
        IF i EQ start_col THEN iplot, freq[0:N/2.-1]/freq_base, line_abs[i-start_col, *], color='red', xrange=[0, 10], xtitle='freq', ytitle='Amp', title='FFT', /NO_SAVEPROMPT ELSE iplot, freq[0:N/2.-1]/freq_base, line_abs[i-start_col, *], color=35000*i, /overplot
      ENDIF
    ENDFOR
    
    IF nrfile EQ 1 THEN BEGIN
      ;plot the 1st LEO
      ind_1st = CEIL(freq_base*1./freq_delta)
      max_r = MAX(line_abs[*, ind_1st])
      iplot, REFORM(line_real[*, ind_1st])/max_r, REFORM(line_imag[*, ind_1st])/max_r, xrange=[-1.1, 1.1], yrange=[-1.1, 1.1], thick=2, color='blue', xtitle='real', ytitle='imag', title='LEOs ('+ title+')', /NO_SAVEPROMPT
      
      ;plot the 3rd LEO
      ind_3rd = CEIL(freq_base*3./freq_delta)
      max_r = MAX(line_abs[*, ind_3rd]);MAX(ABS([line_real[*, ind_3rd], line_imag[*, ind_3rd]]))
      iplot, REFORM(line_real[*, ind_3rd])/max_r, REFORM(line_imag[*, ind_3rd])/max_r, thick=2, color='red', /OVERPLOT
      
      ;plot the 5th LEO
      ind_5th = CEIL(freq_base*5./freq_delta)
      max_r = MAX(line_abs[*, ind_5th])
      iplot, REFORM(line_real[*, ind_5th])/max_r, REFORM(line_imag[*, ind_5th])/max_r, thick=2, color='green', /OVERPLOT
    ENDIF ELSE BEGIN
      
      ;only plot the 3rd LEO
      ind_3rd = CEIL(freq_base*3./freq_delta)
      max_r = MAX(line_abs[*, ind_3rd]);MAX(ABS([line_real[*, ind_3rd], line_imag[*, ind_3rd]]))
      IF i_file EQ 0 THEN BEGIN
        iplot, REFORM(line_real[*, ind_3rd])/max_r, REFORM(line_imag[*, ind_3rd])/max_r, xrange=[-0.2, 1.1], yrange=[-0.2, 1.1], thick=2, color='red', xtitle='real', ytitle='imag', title='G3 CEOs ('+ title+')', /NO_SAVEPROMPT
      ENDIF ELSE BEGIN
        iplot, REFORM(line_real[*, ind_3rd])/max_r, REFORM(line_imag[*, ind_3rd])/max_r, xrange=[-0.2, 1.1], yrange=[-0.2, 1.1], thick=2, color=35000*i_file, /OVERPLOT
      ENDELSE
    ENDELSE
  ENDFOR
  
END

PRO GnRI_2D_PSF, GnR_array, GnI_array, nrfile, offset_field_array, Regular_PSF_array, Donut_PSF_array, STED_PSF_array, title, signal_peak_g_array, fwhm_g_array, signal_peak_STED_array, fwhm_STED_array
  
  n_GnRI = N_ELEMENTS(GnR_array[0,*])
  fwhm_g_array = DBLARR(n_GnRI)
  signal_peak_g_array = DBLARR(n_GnRI)
  fwhm_STED_array = DBLARR(n_GnRI)
  signal_peak_STED_array = DBLARR(n_GnRI)
  n_y = 256.;128;
  Regular_PSF_array = DBLARR(n_y, n_GnRI)
  Donut_PSF_array = DBLARR(n_y, n_GnRI)
  STED_PSF_array = DBLARR(n_y, n_GnRI)
  h_y_array = (FINDGEN(n_y)-n_y/2) / (n_y/2) * offset_field_array[N_ELEMENTS(offset_field_array)-1]
  H_y_2D_array = DBLARR(n_y, n_y)
  H_z_2D_array = DBLARR(n_y, n_y)
  FOR m=0, n_y-1 DO BEGIN
    H_y_2D_array[m, *] = H_y_array
    H_z_2D_array[*, m] = H_y_array
  ENDFOR
  
  FOR i=0, nrfile-1 DO BEGIN
    GnR = REFORM(GnR_array[*, i])
    GnI = REFORM(GnI_array[*, i])
    Harmonic_Real_yz_2D_array = DBLARR(n_y, n_y) * 0.
    Harmonic_Imag_yz_2D_array = DBLARR(n_y, n_y) * 0.
    FOR m=0, n_y-1 DO BEGIN
      FOR n=0, n_y-1 DO BEGIN
        H_temp = SQRT(H_y_2D_array[m, n]^2 + H_z_2D_array[m, n]^2)
        ;      IF ABS(H_temp-A_G*0.1) GT 0.02 THEN BEGIN
        Harmonic_Real_yz_2D_array[m, n] = INTERPOL(GnR, offset_field_array, H_temp, /SPLINE)
        Harmonic_Imag_yz_2D_array[m, n] = INTERPOL(GnI, offset_field_array, H_temp, /SPLINE)
        ;      ENDIF
      ENDFOR
    ENDFOR
    
    ;@GJ, 2025/4/10, calculate the dog factor and STED
    Regular_PSF = REFORM(Harmonic_Real_yz_2D_array[n_y/2, *])
    Regular_PSF_array[*, i] = Regular_PSF
    Donut_PSF = REFORM(Harmonic_Imag_yz_2D_array[n_y/2, *])
    Donut_PSF_array[*, i] = Donut_PSF
    max_Donut_PSF = MAX(ABS(Donut_PSF), max_Donut_PSF_id)
    ;calculate dog factor
    dog_factor = ABS(Regular_PSF[max_Donut_PSF_id]) / max_Donut_PSF
    sted_PSF = Regular_PSF - Donut_PSF * dog_factor
    STED_PSF_array[*, i] = STED_PSF
    offset_field_fov_imag_1d = REFORM(H_y_2D_array[n_y/2, *])
    yfit_psf_g = GAUSSFIT(offset_field_fov_imag_1d, Regular_PSF, coeff_g, NTERMS=4)
    fwhm_g_array[i] = 2*SQRT(2*ALOG(2))*coeff_g[2]
    signal_peak_g_array[i] = MAX(ABS(Regular_PSF))
    yfit_psf_STED = GAUSSFIT(offset_field_fov_imag_1d, sted_PSF, coeff_STED, NTERMS=4)
    fwhm_STED_array[i] = 2*SQRT(2*ALOG(2))*coeff_STED[2]
    signal_peak_STED_array[i] = MAX(ABS(STED_PSF))
;    print, 'Gaussian Peak: ', signal_peak_g_array[i], ' A.U., STED Peak: ', signal_peak_STED_array[i], ' A.U.'
;    print, 'Gaussian FWHM: ', fwhm_g_array[i], ' mT, STED FWHM: ', fwhm_STED_array[i], ' mT'
;    iplot, offset_field_fov_imag_1d, Regular_PSF, color='red', SYM_INDEX=0, THICK=2, LINESTYLE=0, xtitle='x-y [mT]', ytitle='G3R_I [A.U.]', title='PSF', /NO_SAVEPROMPT
;    iplot, offset_field_fov_imag_1d, Donut_PSF, color='green', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
;    iplot, offset_field_fov_imag_1d, sted_PSF, color='blue', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
    
    IF i EQ 0 THEN BEGIN
      iimage, BYTSCL(Harmonic_Real_yz_2D_array, max=MAX(GnR_array)), view_title=title[i]+' Real', RGB_TABLE=3, VIEW_GRID=[2,nrfile], DIMENSIONS=[n_y*2.2, n_y*nrfile*1.2], WINDOW_TITLE='GnRI', /NO_SAVEPROMPT
      iimage, BYTSCL(Harmonic_Imag_yz_2D_array, max=MAX(GnI_array)), view_title=title[i]+' Imag', RGB_TABLE=8, /VIEW_NEXT
    ENDIF ELSE BEGIN
      iimage, BYTSCL(Harmonic_Real_yz_2D_array, max=MAX(GnR_array)), view_title=title[i]+' Real', RGB_TABLE=3, /VIEW_NEXT
      iimage, BYTSCL(Harmonic_Imag_yz_2D_array, max=MAX(GnI_array)), view_title=title[i]+' Imag', RGB_TABLE=8, /VIEW_NEXT
    ENDELSE

  ENDFOR
  
END


;@GJ, 2025/4/7, analyze csv file converted by Yongchen and convert to CEOs (curves of excitation offset)
PRO CEO_simulation_data_analysis_Zhao_Zhiyuan_data

  B_offset_file = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\Dp20nm_Bac10mT\B_bias.csv'
  temp_B = READ_CSV(B_offset_file, HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
  offset_field_array = temp_B.(0)
  min=MIN(ABS(offset_field_array-12.), minind)
  ;N_B_offset = minind
  N_B_offset = N_ELEMENTS(offset_field_array)
  
  ques = DIALOG_MESSAGE("Select a folder?", /question)
  IF STRCMP(ques, 'Yes') THEN BEGIN
    real_dir = dialog_pickfile(TITLE='Select a directory', /DIRECTORY, /MUST_EXIST, PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\')
    ;count the number of files
    real_files = FILE_SEARCH(real_dir, '*kHz.csv', count=nrfile)
    title_temp = STRSPLIT(real_dir, '\', /EXTRACT, COUNT=n_dir)
    title = title_temp[n_dir-1]
    ;@GJ, 2025/4/9, define the rotated RI
    GnR = DBLARR(N_B_offset, nrfile)
    GnI = DBLARR(N_B_offset, nrfile)
    title_array=STRARR(nrfile)
  ENDIF ELSE BEGIN
    real_file = dialog_pickfile(TITLE='Select a csv file', /MUST_EXIST, PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\', FILTER='*.csv')
    ;count the number of files
    real_files = [real_file];['C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\data_file\Dp25nm_Bac10mT\f20kHz.csv']
    nrfile = N_ELEMENTS(real_files)
    title = FILE_BASENAME(real_file, '.csv')
    ;@GJ, 2025/4/9, define the rotated RI from 3rd to 9th
    GnR = DBLARR(N_B_offset, 4)
    GnI = DBLARR(N_B_offset, 4)
    title_array=STRARR(4)
  ENDELSE
  
  GnRI_ind = 0.
  FOR i_file = 0, nrfile-1 DO BEGIN
    ;  real_file = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\', /MUST_EXIST, TITLE="Real csv file", FILTER='*.csv')

    real_data = READ_CSV(real_files[i_file], HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
    n_tag = N_TAGS(real_data)
    n_datapoints = N_ELEMENTS(real_data.(0))

    khz_pos = STRPOS(FILE_BASENAME(real_files[i_file], '.csv'), 'kHz')
    f_pos = STRPOS(FILE_BASENAME(real_files[i_file], '.csv'), 'f')
    freq_base = DOUBLE(STRMID(FILE_BASENAME(real_files[i_file], '.csv'), f_pos+1, khz_pos-f_pos-1));20. ;kHz
    print, real_files[i_file]
    print, 'base freq: ', freq_base

    start_ind = 1.;1. / freq_base / 2. / ((real_data.(1))[1000] - (real_data.(1))[0]) + 1.;0.;n_datapoints - 5000. * 8. - 2. - 2500. ;starting from the middle
    t_array = DBLARR(n_datapoints - start_ind)
    t_array = ((real_data.(1))[start_ind:*] - (real_data.(1))[start_ind]) * 1000. ;ms

    T = ((real_data.(1))[start_ind+1000] - (real_data.(1))[start_ind]) ;ms
    N = N_ELEMENTS(t_array)
    X = (FINDGEN((N - 1)/2) + 1)
    is_N_even = (N MOD 2) EQ 0
    if (is_N_even) then $
      freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
    else $
      freq = [0.0, X, -(N/2 + 1) + X]/(N*T)
    freq_delta = 1./(N*T)

    line_abs = DBLARR(N_B_offset, N/2.) * 0.
    line_real = DBLARR(N_B_offset, N/2.) * 0.
    line_imag = DBLARR(N_B_offset, N/2.) * 0.
    line_phase = DBLARR(N_B_offset, N/2.) * 0.

    start_col = 2
    FOR i=start_col, N_B_offset+start_col-1 DO BEGIN
      temp_signal = -((real_data.(i))[start_ind:*] - (real_data.(i))[(start_ind-1):(n_datapoints-2)])/T
      temp_signal_fft = FFT(temp_signal)
      line_abs[i-start_col, *] = ABS(temp_signal_fft[0:N/2.-1])
      line_real[i-start_col, *] = Real_part(temp_signal_fft[0:N/2.-1])
      line_imag[i-start_col, *] = Imaginary(temp_signal_fft[0:N/2.-1])
      line_phase[i-start_col, *] = 180./!PI*ATAN(temp_signal_fft[0:N/2.-1], /phase)
      IF nrfile EQ 1 THEN BEGIN
        IF i EQ start_col THEN iplot, freq[0:N/2.-1]/freq_base, line_abs[i-start_col, *], color='red', xrange=[0, 10], xtitle='freq', ytitle='Amp', title='FFT', /NO_SAVEPROMPT ELSE iplot, freq[0:N/2.-1]/freq_base, line_abs[i-start_col, *], color=35000*i, /overplot
      ENDIF
    ENDFOR
    
    ;@GJ, 2025/4/7, do the theoretical calculation
    ;@GJ, 2025/1/12, calculate the tau_ms
    viscosity = 1.0;26.9; mPa.s
    T_p = 25.; degree
    eta = viscosity / 1000. ; Pa.s
    particle_size = 25.;130.;70.;30. ; in nm
    D_h = particle_size ;  print, 'D_h: ', D_h, ' nm'
    D_h *= 1.e-9 ; m
    V_h = 1./6. * !PI * D_h^3
    k_B = 1.380469e-23; J/K
    T_p_kelvin = 273.15 + T_p ; in Kelvin
    Msat_T = 0.551 ; Time/u0
    beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
    B_AC = 10.; mT
    f = 25.0; kHz
    w = 2. * !PI * (f / 10^3) ; in us-1
;    tau_B_0 = 1. / SQRT(1. + 0.126 * (beta * B_AC)^1.72) * 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
    tau_B_0 = 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
    tau_B_array = DBLARR(N_B_offset)
    B_array = SQRT(offset_field_array^2 + B_AC^2)
    for i_tau=0, N_B_offset-1 DO tau_B_array[i_tau] = fields_brownian_time_calc(tau_B_0, beta, B_AC, 0., B_array[i_tau])
    print, 'tau_us_0: ', tau_B_array
    signal_ratio_array = B_AC/B_array
    
    IF nrfile EQ 1 THEN BEGIN
      ;plot the 1st LEO
      ind_1st = CEIL(freq_base*1./freq_delta)
      angle = ATAN(line_imag[0, ind_1st], line_real[0, ind_1st], /phase)
      temp_imag = -REFORM(line_real[*, ind_1st]) * sin(angle) + REFORM(line_imag[*, ind_1st]) * cos(angle)
      temp_real = REFORM(line_imag[*, ind_1st]) * sin(angle) + REFORM(line_real[*, ind_1st]) * cos(angle)
      iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), /nodata, xrange=[-0.1, 1.1], yrange=[-0.1, 0.4], thick=2, sym_index = 0, linestyle = 0, color='dark red', xtitle='Real', ytitle='Imag', title='CEOs ('+ title+')', /NO_SAVEPROMPT
      iplot, [-0.1, 1.1], [0, 0], linestyle = 2, /overplot
      iplot, [0, 0], [-0.1, 1.1], linestyle = 2, /overplot
      n = 1.
      tan_angle = n * w * (tau_B_array[0] - tau_B_array) / (1. + n^2 * w^2 * tau_B_array * tau_B_array[0])
      fit_imag = signal_ratio_array * line_abs[0, ind_1st] * sin(atan(tan_angle))
      fit_real = signal_ratio_array * line_abs[0, ind_1st] * cos(atan(tan_angle))
;      iplot, fit_real/MAX(fit_real), fit_imag/MAX(fit_imag)*MAX(temp_imag)/MAX(temp_real), color='dark red', /overplot
     
      ;plot the 3rd LEO
      ind_3rd = CEIL(freq_base*3./freq_delta)
      angle = ATAN(line_imag[0, ind_3rd], line_real[0, ind_3rd], /phase)
      temp_imag = -REFORM(line_real[*, ind_3rd]) * sin(angle) + REFORM(line_imag[*, ind_3rd]) * cos(angle)
      temp_real = REFORM(line_imag[*, ind_3rd]) * sin(angle) + REFORM(line_real[*, ind_3rd]) * cos(angle)
      iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), thick=2, sym_index = 2, linestyle = 0, color='blue', /OVERPLOT
      title_array[0] = '3rd Harmonic'
      GnR[*, 0] = temp_real
      GnI[*, 0] = temp_imag
;      GnRI_2D_PSF, temp_real, temp_imag, offset_field_array, Regular_psf, Donut_PSF, Harmonic_Real_yz_2D_array, Harmonic_Imag_yz_2D_array, title
      n = 3.
      tan_angle = n * w * (tau_B_array[0] - tau_B_array) / (1. + n^2 * w^2 * tau_B_array * tau_B_array[0])
      fit_imag = signal_ratio_array * line_abs[0, ind_3rd] * sin(atan(tan_angle))
      fit_real = signal_ratio_array * line_abs[0, ind_3rd] * cos(atan(tan_angle))
      iplot, fit_real/MAX(ABS(fit_real)), fit_imag/MAX(fit_imag)*MAX(temp_imag)/MAX(temp_real), color='blue', /overplot

      ;plot the 5th LEO
      ind_5th = CEIL(freq_base*5./freq_delta)
      angle = ATAN(line_imag[0, ind_5th], line_real[0, ind_5th], /phase)
      temp_imag = -REFORM(line_real[*, ind_5th]) * sin(angle) + REFORM(line_imag[*, ind_5th]) * cos(angle)
      temp_real = REFORM(line_imag[*, ind_5th]) * sin(angle) + REFORM(line_real[*, ind_5th]) * cos(angle)
      iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), thick=2, sym_index = 2, linestyle = 0, color='green', /OVERPLOT
      title_array[1] = '5th Harmonic'
      GnR[*, 1] = temp_real
      GnI[*, 1] = temp_imag
;      title='5th Harmonic'
;      GnRI_2D_PSF, temp_real, temp_imag, offset_field_array, Regular_psf, Donut_PSF, Harmonic_Real_yz_2D_array, Harmonic_Imag_yz_2D_array, title
      n = 5.
      tan_angle = n * w * (tau_B_array[0] - tau_B_array) / (1. + n^2 * w^2 * tau_B_array * tau_B_array[0])
      fit_imag = signal_ratio_array * line_abs[0, ind_5th] * sin(atan(tan_angle))
      fit_real = signal_ratio_array * line_abs[0, ind_5th] * cos(atan(tan_angle))
      iplot, fit_real/MAX(fit_real), fit_imag/MAX(fit_imag)*MAX(temp_imag)/MAX(temp_real), color='green', /overplot

      ;plot the 7th LEO
      ind_7th = CEIL(freq_base*7./freq_delta)
      angle = ATAN(line_imag[0, ind_7th], line_real[0, ind_7th], /phase)
      temp_imag = -REFORM(line_real[*, ind_7th]) * sin(angle) + REFORM(line_imag[*, ind_7th]) * cos(angle)
      temp_real = REFORM(line_imag[*, ind_7th]) * sin(angle) + REFORM(line_real[*, ind_7th]) * cos(angle)
      iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), thick=2, sym_index = 2, linestyle = 0, color='red', /OVERPLOT
      title_array[2] = '7th Harmonic'
      GnR[*, 2] = temp_real
      GnI[*, 2] = temp_imag
;      title='7th Harmonic'
;      GnRI_2D_PSF, temp_real, temp_imag, offset_field_array, Regular_psf, Donut_PSF, Harmonic_Real_yz_2D_array, Harmonic_Imag_yz_2D_array, title
      n = 7.
      tan_angle = n * w * (tau_B_array[0] - tau_B_array) / (1. + n^2 * w^2 * tau_B_array * tau_B_array[0])
      fit_imag = signal_ratio_array * line_abs[0, ind_7th] * sin(atan(tan_angle))
      fit_real = signal_ratio_array * line_abs[0, ind_7th] * cos(atan(tan_angle))
      iplot, fit_real/MAX(fit_real), fit_imag/MAX(fit_imag)*MAX(temp_imag)/MAX(temp_real), color='red', /overplot
    
    
      ;plot the 9th LEO
      ind_9th = CEIL(freq_base*9./freq_delta)
      angle = ATAN(line_imag[0, ind_9th], line_real[0, ind_9th], /phase)
      temp_imag = -REFORM(line_real[*, ind_9th]) * sin(angle) + REFORM(line_imag[*, ind_9th]) * cos(angle)
      temp_real = REFORM(line_imag[*, ind_9th]) * sin(angle) + REFORM(line_real[*, ind_9th]) * cos(angle)
      iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), thick=2, sym_index = 2, linestyle = 0, color='dark orange', /OVERPLOT
      title_array[3] = '9th Harmonic'
      GnR[*, 3] = temp_real
      GnI[*, 3] = temp_imag
;      title='9th Harmonic'
;      GnRI_2D_PSF, temp_real, temp_imag, offset_field_array, Regular_psf, Donut_PSF, Harmonic_Real_yz_2D_array, Harmonic_Imag_yz_2D_array, title
      n = 9.
      tan_angle = n * w * (tau_B_array[0] - tau_B_array) / (1. + n^2 * w^2 * tau_B_array * tau_B_array[0])
      fit_imag = signal_ratio_array * line_abs[0, ind_9th] * sin(atan(tan_angle))
      fit_real = signal_ratio_array * line_abs[0, ind_9th] * cos(atan(tan_angle))
      iplot, fit_real/MAX(fit_real), fit_imag/MAX(fit_imag)*MAX(temp_imag)/MAX(temp_real), color='dark orange', /overplot
      
      GnRI_2D_PSF, GnR, GnI, 4, offset_field_array, Regular_psf_array, Donut_PSF_array, STED_PSF_array, title_array, signal_peak_g_array, fwhm_g_array, signal_peak_STED_array, fwhm_STED_array
      print, 'Gaussian Peak: ', signal_peak_g_array, ' A.U.'
      print, 'STED Peak: ', signal_peak_STED_array, ' A.U.'
      print, 'Gaussian FWHM: ', fwhm_g_array, ' mT'
      print, 'STED FWHM: ', fwhm_STED_array, ' mT'
;    iplot, offset_field_fov_imag_1d, Regular_PSF, color='red', SYM_INDEX=0, THICK=2, LINESTYLE=0, xtitle='x-y [mT]', ytitle='G3R_I [A.U.]', title='PSF', /NO_SAVEPROMPT
;    iplot, offset_field_fov_imag_1d, Donut_PSF, color='green', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
;    iplot, offset_field_fov_imag_1d, sted_PSF, color='blue', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT


    ENDIF ELSE BEGIN

      ;only plot the 3rd LEO
      ind_3rd = CEIL(freq_base*3./freq_delta)
      angle = ATAN(line_imag[0, ind_3rd], line_real[0, ind_3rd], /phase)
      temp_imag = -REFORM(line_real[*, ind_3rd]) * sin(angle) + REFORM(line_imag[*, ind_3rd]) * cos(angle)
      temp_real = REFORM(line_imag[*, ind_3rd]) * sin(angle) + REFORM(line_real[*, ind_3rd]) * cos(angle)
     IF i_file EQ 0 THEN BEGIN
        iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), xrange=[-0.1, 1.1], yrange=[-0.1, 1.1], thick=2, color='red', xtitle='G3R', ytitle='G3I', title='CEOs ('+ title+')', /NO_SAVEPROMPT
        iplot, [-0.1, 1.1], [0, 0], linestyle = 2, /overplot
        iplot, [0, 0], [-0.1, 1.1], linestyle = 2, /overplot
      ENDIF ELSE BEGIN
        iplot, temp_real/MAX(temp_real), temp_imag/MAX(temp_real), thick=2, color=35000*i_file, /OVERPLOT
      ENDELSE
      
      
;      IF (freq_base-8.) mod 4 EQ 0 THEN BEGIN
        ;@GJ, 2025/4/9, save the result
        title_array[GnRI_ind] = freq_base
        GnR[*, GnRI_ind] = temp_real
        GnI[*, GnRI_ind] = temp_imag
        GnRI_ind++
;      ENDIF  
    ENDELSE
  ENDFOR
  
  IF nrfile NE 1 THEN BEGIN
;    GnRI_2D_PSF, GnR, GnI, MIN([4, nrfile]), offset_field_array, Regular_psf_array, Donut_PSF_array, STED_PSF_array, title_array, signal_peak_g_array, fwhm_g_array, signal_peak_STED_array, fwhm_STED_array
    GnRI_2D_PSF, GnR, GnI, nrfile, offset_field_array, Regular_psf_array, Donut_PSF_array, STED_PSF_array, title_array, signal_peak_g_array, fwhm_g_array, signal_peak_STED_array, fwhm_STED_array
    print, 'Gaussian Peak: ', signal_peak_g_array, ' A.U.'
    print, 'STED Peak: ', signal_peak_STED_array, ' A.U.'
    print, 'Gaussian FWHM: ', fwhm_g_array, ' mT'
    print, 'STED FWHM: ', fwhm_STED_array, ' mT'
  ENDIF


END

;@GJ, 2024/9/14, analyze Dr. Zhou's csv data
PRO MPI_zhao_data_analyze
  real_file = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhaoSimu\', /MUST_EXIST, TITLE="Real csv file", FILTER='*.csv')
  real_data = READ_CSV(real_file, HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
  n_tag = N_TAGS(real_data)
  n_datapoints = N_ELEMENTS(real_data.(0))
  offset_field_array = DOUBLE(SedHeader[2:12])
  
  start_ind = n_datapoints - 5000. * 8. - 2. - 2500. ;starting from the middle
  t_array = ((real_data.(1))[start_ind:*] - (real_data.(1))[start_ind]) * 1000. ;ms
  T = ((real_data.(1))[start_ind+1000] - (real_data.(1))[start_ind]) ;ms
  N = N_ELEMENTS(t_array)
  X = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
  else $
    freq = [0.0, X, -(N/2 + 1) + X]/(N*T)
  
  freq_base = 20. ;kHz
  line_abs = DBLARR(11, N) * 0.
  line_real = DBLARR(11, N) * 0.
  line_imag = DBLARR(11, N) * 0.
  line_phase = DBLARR(11, N) * 0.
  FOR i=2, 12 DO BEGIN
    temp_signal = -((real_data.(i))[start_ind:*] - (real_data.(i))[(start_ind-1):(n_datapoints-2)])/T
    temp_signal_fft = FFT(temp_signal)
    line_abs[i-2, *] = ABS(temp_signal_fft)
    line_real[i-2, *] = Real_part(temp_signal_fft)
    line_imag[i-2, *] = Imaginary(temp_signal_fft)
    line_phase[i-2, *] = 180./!PI*ATAN(temp_signal_fft, /phase)
    IF i EQ 2 THEN iplot, freq/freq_base, line_abs[i-2, *], color='red', xrange=[0, 10], xtitle='freq', ytitle='Amp', title='FFT', /NO_SAVEPROMPT ELSE iplot, freq/freq_base, line_abs[i-2, *], color=35000*i, /overplot
  ENDFOR
  
  ind_1st = freq_base*1./freq[1]
  max_r = MAX(ABS([line_real[*, ind_1st], line_imag[*, ind_1st]]))
  iplot, REFORM(line_real[*, ind_1st])/max_r, REFORM(line_imag[*, ind_1st])/max_r, thick=2, color='blue', xtitle='real', ytitle='imag', title='LEOs', /NO_SAVEPROMPT
  
  ind_3rd = freq_base*3./freq[1]
  max_r = MAX(ABS([line_real[*, ind_3rd], line_imag[*, ind_3rd]]))
  iplot, REFORM(line_real[*, ind_3rd])/max_r, REFORM(line_imag[*, ind_3rd])/max_r, thick=2, color='red', /overplot; xtitle='real', ytitle='imag', title='3rd LEO', /NO_SAVEPROMPT

  ind_5th = freq_base*5./freq[1]
  max_r = MAX(ABS([line_real[*, ind_5th], line_imag[*, ind_5th]]))
  iplot, REFORM(line_real[*, ind_5th])/max_r, REFORM(line_imag[*, ind_5th])/max_r, thick=2, color='green', /overplot; xrange=[-max_r, max_r], yrange=[-max_r, max_r], xtitle='real', ytitle='imag', title='5th LEO', /NO_SAVEPROMPT

;  iplot, t_array, (real_data.(5))[start_ind:*]
  
END


;@GJ, 2024/7/31, read csv data from Tianshu Li
PRO MPI_tianshu_data_analyze
  real_file = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\baishi_tianshuli\', /MUST_EXIST, TITLE="Real csv file", FILTER='*.csv')
  real_data = READ_CSV(real_file, HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
  n_tag = N_TAGS(real_data)
  real_image = DBLARR(n_tag, N_ELEMENTS(real_data.(0)))
  for i=0, n_tag-1 do real_image[i,*]=REFORM(real_data.(i))
;  real_image = shift(real_image, -5, -11)
;  iimage, shift(real_image, -5, -11), title='Real'
  
  file_path = FILE_DIRNAME(real_file)
  imag_file = DIALOG_PICKFILE(PATH=file_path, /MUST_EXIST, TITLE="Imaginary csv file", FILTER='*.csv')
  imag_data = READ_CSV(imag_file, HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
  n_tag = N_TAGS(imag_data)
  imag_image = DBLARR(n_tag, N_ELEMENTS(imag_data.(0)))
  for i=0, n_tag-1 do imag_image[i,*]=REFORM(imag_data.(i))
;  imag_image = shift(imag_image, -5, -11)
;  iimage, shift(imag_image, -5, -11), title='Imag'
  
  amp_image =  DBLARR(n_tag, N_ELEMENTS(imag_data.(0)))
  phase_image = DBLARR(n_tag, N_ELEMENTS(imag_data.(0)))
  for i=0, n_tag-1 do begin
    for j=0, N_ELEMENTS(imag_data.(0))-1 do begin
      temp_real = real_image[i, j]
      temp_imag = imag_image[i, j]
      amp_image[i, j] = SQRT(temp_real^2 + temp_imag^2)
      phase_image[i, j] = 180./!PI*ATAN(temp_imag, temp_real)
      IF phase_image[i, j] LT 0 THEN phase_image[i, j] += 360. 
    endfor
  endfor
  
  iimage, amp_image, view_title='3rd Amp', VIEW_GRID=[5,7], DIMENSIONS=[n_tag*6., N_ELEMENTS(imag_data.(0))*8.], WINDOW_TITLE='3rd Harmonic', /NO_SAVEPROMPT
  iimage, real_image, view_title='3rd Real', /VIEW_NEXT
  iimage, imag_image, view_title='3rd Imag', /VIEW_NEXT
  iimage, BYTSCL(phase_image), view_title='3rd Phase', /VIEW_NEXT
  iimage, 255.-BYTSCL(phase_image), view_title='3rd -Phase', /VIEW_NEXT
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (i*10.)/180.
    new_image = -real_image * sin(phi) + imag_image * cos(phi)
    iimage, BYTSCL(new_image), view_title='3rd New', /VIEW_NEXT
  ENDFOR
  
  i=5
  phi = -!PI * (i*10.)/180.
  new_image_imag = -real_image * sin(phi) + imag_image * cos(phi)
  new_image_real = imag_image * sin(phi) + real_image * cos(phi)
  n=15
  mean_imag = DBLARR(n)
  mean_real = DBLARR(n)
  linear_ratio = DBLARR(n)
  FOR line = 21, 21+n-1 DO BEGIN; No DC offset
    mean_imag[line-21] = MEAN(new_image_imag[line,*])
    mean_real[line-21] = MEAN(new_image_real[line,*])
    fit_result = LINFIT(new_image_real, new_image_imag, CHISQR=fit_chisqr)
    linear_ratio[line-21] = fit_chisqr
    max_imag = MAX(ABS(new_image_imag))
    max_real = MAX(ABS(new_image_real))
    wait, 2
    IF line EQ 21 THEN BEGIN
      iplot, new_image_real[line,*], new_image_imag[line,*], xtitle='Real', ytitle='Imag', title='Scatter Plot', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='red'
    ENDIF ELSE BEGIN
      iplot, new_image_real[line,*], new_image_imag[line,*], xrange=[-max_real, max_real], yrange=[-max_imag, max_imag], SYM_INDEX=2, LINESTYLE=0, color='blue', /overplot
    ENDELSE
  ENDFOR
  
  ;@GJ, 2024/8/31, plot the phase under different DC offset
  iplot, mean_real, mean_imag, xtitle='Real', ytitle='Imag', yerror=linear_ratio/MAX(linear_ratio)*30., title='Mean Plot vs DC Offset'
  
  line=25; with DC offset
  i=23
  phi = -!PI * (i*10.)/180.
  new_image_imag = -real_image * sin(phi) + imag_image * cos(phi)
  max_imag = MAX(ABS(new_image_imag))
  iplot, new_image_imag[line,*], xtitle='Pixel No.', ytitle='Imag', title='Line Profile (with DC Offset)', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='blue'
  
  line = 21; DC offset
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (i*10.)/180.
    new_image_imag = -real_image * sin(phi) + imag_image * cos(phi)
    new_image_real = imag_image * sin(phi) + real_image * cos(phi)
    IF i EQ 1 THEN BEGIN
      max_imag = MAX(ABS(new_image_imag))
      max_real = MAX(ABS(new_image_real))
      iplot, new_image_real[line,*], new_image_imag[line,*], xtitle='Real', ytitle='Imag', title='3rd Harmonic (No DC Offset)', xrange=[-max_real, max_real], yrange=[-max_imag, max_imag], SYM_INDEX=2, LINESTYLE=0, color='red'
;      iplot, new_image_imag[line,*], xtitle='Pixel No.', ytitle='Imag', title='Line Profile (No DC Offset)', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='red'
    ENDIF ELSE BEGIN
      max_imag = MAX([max_imag, MAX(ABS(new_image_imag))])
      max_real = MAX([max_real, MAX(ABS(new_image_real))])
;      iplot, new_image_real[line,*], new_image_imag[line,*], xtitle='Real', ytitle='Imag', xrange=[-max_real, max_real], yrange=[-max_imag, max_imag], SYM_INDEX=2, LINESTYLE=0, color='blue', /overplot
;      iplot, new_image_imag[line,*], xtitle='Pixel No.', ytitle='Imag', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='blue', /overplot
    ENDELSE
  ENDFOR

  line = 25 ; DC offset
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (i*10.)/180.
    new_image_imag = -real_image * sin(phi) + imag_image * cos(phi)
    new_image_real = imag_image * sin(phi) + real_image * cos(phi)
    IF i EQ 1 THEN BEGIN
      iplot, new_image_real[line,*], new_image_imag[line,*], xtitle='Real', ytitle='Imag', title='3rd Harmonic (with DC Offset)', xrange=[-max_real, max_real], yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='red'
;      iplot, new_image_imag[line,*], xtitle='Pixel No.', ytitle='Imag', title='Line Profile (with DC Offset)', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='red'
    ENDIF ELSE BEGIN
;      iplot, new_image_real[line,*], new_image_imag[line,*], xtitle='Real', ytitle='Imag', SYM_INDEX=3, LINESTYLE=0, color='blue', /overplot
;      iplot, new_image_imag[line,*], xtitle='Pixel No.', ytitle='Imag', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='blue', /overplot
    ENDELSE
  ENDFOR
  
  iimage, amp_image, view_title='3rd Amp', VIEW_GRID=[5,7], DIMENSIONS=[n_tag*6., N_ELEMENTS(imag_data.(0))*8.], WINDOW_TITLE='3rd Harmonic', /NO_SAVEPROMPT
  iimage, real_image, view_title='3rd Real', /VIEW_NEXT
  iimage, imag_image, view_title='3rd Imag', /VIEW_NEXT
  iimage, BYTSCL(phase_image), view_title='3rd Phase', /VIEW_NEXT
  iimage, 255.-BYTSCL(phase_image), view_title='3rd -Phase', /VIEW_NEXT
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (50.+i)/180.
    new_image = -real_image * sin(phi) + imag_image * cos(phi)
    iimage, BYTSCL(new_image), view_title='3rd New', /VIEW_NEXT
  ENDFOR
  
  gap = 3.
  combine_image = DBLARR(4*(n_tag+gap), N_ELEMENTS(imag_data.(0)))
  combine_image[0:n_tag-1, *] = BYTSCL(amp_image)
  combine_image[n_tag+gap:2*n_tag+gap-1, *] = BYTSCL(real_image)
  combine_image[2*n_tag+2*gap:3*n_tag+2*gap-1, *] = BYTSCL(imag_image)
  combine_image[3*n_tag+3*gap:4*n_tag+3*gap-1, *] = BYTSCL(phase_image)
  
  iimage, combine_image, title='combined', /NO_SAVEPROMPT
  ;@GJ, 2024/7/24, save the results
  combine_image_file = file_path + '\combined_image.jpg'
  WRITE_JPEG, combine_image_file, combine_image
  
  ;@GJ, 2024/8/10, check the line profile of central
  iimage, phase_image, title='Phase Image', /NO_SAVEPROMPT
  iplot, phase_image[23,*], color='red', /NO_SAVEPROMPT
  iplot, phase_image[21,*], color='green', /overplot
  iplot, phase_image[23,*]-phase_image[21,*], color='blue', /overplot
    
  ;@GJ, 2024/8/6, calculate the pfps
  map_phase_flip_line = phase_image - shift(phase_image, 1, 0)
  mask_pfp = (ABS(map_phase_flip_line) GT 100.)
  map_phase_flip_line *= mask_pfp
  
  ;@GJ, 2024/8/6, shift the pfps
  temp_image = amp_image * 0.
  size_image = size(amp_image, /dim)
  A_G_half = 2.
  ind_mask = WHERE(mask_pfp EQ 1, n_ind)
  IF n_ind GT 1 THEN BEGIN
    FOR i=0, n_ind-1 DO BEGIN
      ind_2d = ARRAY_INDICES(amp_image, ind_mask[i])
      IF map_phase_flip_line[ind_mask[i]] GT 0 THEN BEGIN ;LT for frequency-mixing, GT for single freq
        shift_range = ind_2d[0] + A_G_half
        IF shift_range LT  size_image[0] THEN temp_image[shift_range, ind_2d[1]] = 1.
      ENDIF ELSE BEGIN
        shift_range = ind_2d[0] - A_G_half
        IF shift_range GE  0 THEN temp_image[shift_range, ind_2d[1]] = 1.
      ENDELSE
    ENDFOR
  ENDIF
  
  ;@GJ, display the temp_image
  iimage, BYTSCL(temp_image*amp_image), /NO_SAVEPROMPT
END

;@GJ, 20254/14, do the data analysis
PRO MPI_TianshuLi_JingZhong_data_analysis_202504145

  ;@GJ, 2025/4/14
  real_file = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\baishi_tianshuli\20250414_PSF\', /MUST_EXIST, TITLE="Real csv file", FILTER='*.csv')
;  real_file = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\baishi_tianshuli\20250414_PSF\tidu6 jili10.csv'
;  real_file = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\baishi_tianshuli\20250414_PSF\5mT.csv'
  filename = FILE_BASENAME(real_file, '.csv')
  sp = STRPOS(filename, 'mT')
  excitation = float(STRMID(filename, 0, sp))
  real_data = READ_CSV(real_file, HEADER=SedHeader, N_TABLE_HEADER=1, TABLE_HEADER=SedTableHeader)
  n_tag = N_TAGS(real_data)
  real_image = DBLARR(n_tag, N_ELEMENTS(real_data.(0)))
  for i=0, n_tag-1 do real_image[i,*]=REFORM(real_data.(i))
  max_int = MAX(ABS(real_image), max_loc)
  max_loc_2d = ARRAY_INDICES(ABS(real_image), max_loc)
  temp_real_imag = SHIFT(real_image, -max_loc_2d[0]+n_tag/2, -max_loc_2d[1]+N_ELEMENTS(real_data.(0))/2)
  real_image = temp_real_imag
  iimage, BYTSCL(ABS(real_image), MAX=max_int/3), RGB_TABLE=8, title='PSF'
  max_int = MAX(ABS(real_image), max_loc)
  max_loc_2d = ARRAY_INDICES(ABS(real_image), max_loc)
;  iplot, ABS(real_image[max_loc_2d[0],*])/MAX(ABS(real_image[max_loc_2d[0],*])), xtitle='Pixel Index', ytitle='Signal (A.U.)', color='red', thick=2, title='Line Profiles', /nosaveplot
;  iplot, ABS(real_image[max_loc_2d[0]+5,*])/MAX(ABS(real_image[max_loc_2d[0]+5,*]))*0.35, color='green', thick=2, /overplot
  first_plot = 0
;  FOR i=max_loc_2d[0]+1, n_tag-3 DO BEGIN
;;    iplot, ABS(real_image[36,*])/MAX(ABS(real_image[36,*]))*0.35, color='green', thick=2, /overplot
;    IF first_plot NE 0 THEN BEGIN
;      iplot, ABS(real_image[i-1,*])/MAX(ABS(real_image[i-1,*]))*0.35, color='white', thick=2, /overplot
;      iplot, ABS(real_image[max_loc_2d[0],*])/MAX(ABS(real_image[max_loc_2d[0],*])), color='red', thick=2, /overplot
;    ENDIF
;    iplot, ABS(real_image[i,*])/MAX(ABS(real_image[i,*]))*0.35, color='green', thick=2, /overplot
;    first_plot = 1
;    wait, 1
;  ENDFOR
  
  n_z = N_ELEMENTS(real_data.(0))
  donut_radius = DBLARR(n_tag) * 0.
  donut_thickness = DBLARR(n_tag) * 0.
  resolution = 40. / n_z;
  FOR i=6, n_tag-7 DO BEGIN
    donut = REFORM(real_image[i,0:n_z/2])
    maxD = MAX(donut, maxD_ind)
    IF maxD_ind LT n_z/2 AND donut[n_z/2] LT donut[n_z/2-2]  THEN BEGIN
      ;there is a donut
      donut_radius[i] = (n_z/2 - maxD_ind) * resolution; in mm
      donut_thickness[i] = maxD - donut[n_z/2] ; in random unit
    ENDIF
  ENDFOR
  gradient_field = 0.85; mT/mm
  iplot, (findgen(n_tag)-n_z/2.)*resolution*gradient_field/excitation*100., donut_radius, color='red', thick=2, xrange=[-500,500], xmajor=11, xtitle='H_DC/H_AC (%)', ytitle='Donut Radius (mm)', title='Donut Radius (' + STRMID(filename, 0, sp+2)+')'
  iplot, (findgen(n_tag)-n_z/2.)*resolution*gradient_field/excitation*100., donut_thickness, color='blue', thick=2, xrange=[-500,500], xmajor=11, xtitle='H_DC/H_AC (%)', ytitle='Donut Thickness (A.U.)', title='Donut Thickness (' + STRMID(filename, 0, sp+2)+')'
  
END