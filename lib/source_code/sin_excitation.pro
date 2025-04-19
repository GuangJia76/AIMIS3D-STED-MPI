function rl_deconv1d, signal, psf, iterations=iterations, epsilon=epsilon
  ; 输入参数检查
;  if (size(signal, /dimensions)[0] ne 1) then $
;    message, 'Signal must be 1D array'
;  if (size(psf, /dimensions)[0] ne 1) then $
;    message, 'PSF must be 1D array'
  ; 转换为浮点类型以支持高精度计算
  signal = float(signal)
  psf = float(psf)
  ; PSF归一化（确保能量守恒）
  psf = psf / total(psf)
  ; 初始化估计值为输入信号（加速收敛）
  u = signal
  ; 设置默认参数
  if keyword_set(iterations) EQ 0 then iterations = 50
  if keyword_set(epsilon) EQ 0 then epsilon = 1e-10
  ; 主迭代循环
  for i = 1L, iterations do begin
    ; 前向卷积：u ⊗ psf
    conv_step = convol(u, psf, /center, /edge_zero)
    ; 计算比值：data / (u ⊗ psf + ε)
    ratio = signal / (conv_step + epsilon)
    ; 反向卷积：ratio ⊗ psf_flipped（对称PSF无需翻转）
    backward_conv = convol(ratio, reverse(psf), /center, /edge_zero)
    ; 更新估计：u *= backward_conv
    u = u * backward_conv
  endfor
  return, u
end

PRO RL_Decon_1d

  ; 生成测试信号（理想脉冲）
  n = 128
  signal = fltarr(n)
  signal[n/2.] = 1.0
  signal[n*4./5.] = 1.0

  ; 创建对称高斯PSF
  x = findgen(n) - n/2
  psf = exp(-(x/5.0)^2)
  psf = psf / total(psf)  ; 归一化

  ; 模拟模糊信号
  blurred = convol(signal, psf, /center, /normalize, /edge_zero)

  ; 执行RL解卷积
  restored = rl_deconv1d(blurred, psf, iterations=100)

  ; 可视化结果
;  legend=['Blurred', 'Restored', 'Original']
  iplot, blurred, color='blue', linestyle=1, title='Blurred Signal', insert_legend='blurred'
  iplot, restored, color='red', linestyle=2, /overplot, insert_legend='restored'
  iplot, signal, color='green', linestyle=3, /overplot, insert_legend='signal'


END

function Debye_Kernel_sin, tau, kernel_time
  N = N_ELEMENTS(kernel_time)
  kernel = DBLARR(N)*0.
  d_t = kernel_time[1]-kernel_time[0] ;us
  time = FINDGEN(N) * d_t  ; us
  r = (1./tau) * exp(-time/tau)
  kernel[0:N/2] = REVERSE(r[0:N/2])

  return, kernel
end

;;@GJ, 2025/2/6, calculate the pSNR, SSIM, MSE of two images
;pro image_pSNR_SSIM_RMSE, image1, image2, rmse, psnr, ssim_mean
;
;;  在IDL（Interactive Data Language）中，可以通过以下步骤来计算两幅黑白图像的峰值信噪比（PSNR）、结构相似指数（SSIM）和均方根误差（RMSE）。以下是示例代码：
;;
;;  idl
;;
;  ; 读取两幅图像
;;  image1 = READ_IMAGE('image1.fits')
;;  image2 = READ_IMAGE('image2.fits')
;
;  ; 计算均方根误差 (RMSE)
;  rmse = SQRT(MEAN((image1 - image2)^2))
;
;  ; 计算峰值信噪比 (PSNR)
;  max_pixel = MAX(image1)
;  psnr = 10.0 * alog10(max_pixel^2 / rmse^2)
;
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
;  PRINT, '均方根误差 (RMSE): ', rmse
;  PRINT, '峰值信噪比 (PSNR): ', psnr
;  PRINT, '结构相似指数 (SSIM): ', ssim_mean
;;  代码说明：
;;
;;  图像读取：使用 READ_IMAGE 函数读取两幅黑白图像。这里假设图像格式为FITS，如果是其他格式，需要相应修改读取函数。
;;
;;  均方根误差 (RMSE)：通过计算两幅图像对应像素差值的平方的均值，再取平方根得到RMSE。
;;
;;  峰值信噪比 (PSNR)：利用RMSE和图像的最大像素值计算PSNR。
;;
;;  结构相似指数 (SSIM)：
;;
;;  - 定义一些参数，如 K1 、 K2 和 L 。
;;
;;  - 使用卷积计算图像的均值、方差和协方差。
;;
;;  - 根据公式计算SSIM，并取其均值作为最终结果。
;;
;;  请根据实际图像格式和路径修改 READ_IMAGE 函数中的文件名。
;
;end


;@GJ, 2025/1/24, STED-MPI theory equation
pro G3RI_equation_calculation, x_array, H_x_array, G_3R_psf, G_3I_psf, G_3STED_psf, parameters
  
  ;@GJ, 2025/1/12, calculate the tau_ms
  viscosity = 1.; mPa.s
  T_p = 25.; degree
  eta = viscosity / 1000. ; Pa.s
  particle_size = parameters.particle_size;particle_size = 60. ; in nm
  D_h = particle_size ;  print, 'D_h: ', D_h, ' nm'
  D_h *= 1.e-9 ; m
  V_h = 1./6. * !PI * D_h^3
  k_B = 1.380469e-23; J/K
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  Msat_T = 0.551 ; Time/u0
  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT

  tau_B = 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
  print, 'tau_us_0: ', tau_B
  
  ;@GJ, 2025/1/24, stimulation parameters
  A = parameters.A;5.8 ; mT
  f = parameters.f;2.0; kHz
  w = 2. * !PI * (f / 10^3) ; in us-1

  ;@GJ, 2025/1/24, scan parameters
  n_x = 256.
  fov = 190. ; mm
  resolution = fov / n_x ; mm
  gradient_field = 1.13 ; mT/mm
  H_x_array = (FINDGEN(n_x)-FLOOR(n_x/2.)) * resolution * gradient_field
  x_array = (FINDGEN(n_x)-FLOOR(n_x/2.)) * resolution
  G_3R_array = DBLARR(n_x) * 0.
  G_3I_array = DBLARR(n_x) * 0.
  tau_us_0 = fields_brownian_time_calc(tau_B, beta, A, 0., 0.)
  
  n = 3.
  for i=0, n_x-1 DO BEGIN
    tau_us = fields_brownian_time_calc(tau_B, beta, A, 0., H_x_array[i])
    factor = (1. + (n * w * tau_us)^2) * SQRT(1. + (n * w * tau_us_0)^2)
    G_3R_array[i] = f * (1. + (n * w)^2 * tau_us * tau_us_0) / factor
    G_3I_array[i] = f * n * w * (tau_us_0 -  tau_us)
  endfor
  
  iplot, x_array, G_3R_array/MAX(G_3R_array), yrange=[-0.1, 1.1], color='red', thick=2, xtitle='x-y [mm]', ytitle='R3R & R3I [A.U.]', title='Debye R3', /NO_SAVEPROMPT
  iplot, x_array, G_3I_array/MAX(G_3I_array), color='green', thick=2, /overplot
  
  ; n_x is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((n_x - 1)/2) + 1)
  is_N_even = (n_x MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, n_x/2, -n_x/2 + X]/(n_x*resolution) $
  else $
    freq = [0.0, X, -(n_x/2 + 1) + X]/(n_x*resolution)
;  iplot, freq, ABS(FFT(G_3I_array/MAX(G_3I_array))), xrange=[0, 20./fov], color='blue', xtitle='Freq [1/mm]', ytitle='Amp', title='FFT of G3', /NO_SAVEPROMPT
  
  alpha_3 = 3.2;10.;2.4;10.;2.4;5.0
  beta_3 = 0.32;0.57;1.0
  G_3_ad = beta_3 * alpha_3 / (x_array^2 + alpha_3^2)
  iplot, x_array, G_3_ad/MAX(G_3_ad), color='black', thick=2, linestyle=1, /overplot
  
  iplot, x_array, G_3_ad, color='black', thick=2, /nodata, xtitle='x-y [mm]', ytitle='G3R & G3I [A.U.]', title='G3R and G3I', /NO_SAVEPROMPT
  iplot, x_array, G_3_ad*G_3R_array/MAX(G_3R_array), color='red', thick=2, /overplot
  iplot, x_array, G_3_ad*G_3I_array/MAX(G_3I_array), color='green', thick=2, /overplot
;  iplot, freq, ABS(FFT(G_3_ad)), xrange=[0, 20./fov], color='blue', xtitle='Freq [1/mm]', ytitle='Amp', title='FFT of G3ad', /NO_SAVEPROMPT
;  iplot, freq, ABS(FFT(G_3_ad*G_3I_array)), xrange=[0, 20./fov], color='blue', xtitle='Freq [1/mm]', ytitle='Amp', title='FFT of G3I', /NO_SAVEPROMPT
;  print, 'G3RI plot done'
  
  ;@GJ, 2025/1/31, calculate the theoretical STED curve
  G_3I_array_new = DBLARR(n_x) * 0.
  n = 3.
  for i=0, n_x-1 DO BEGIN
    tau_us = fields_brownian_time_calc(tau_B, beta, A, 0., H_x_array[i])
    ratio_IR = n * w * (tau_us_0 -  tau_us) / (1. + (n * w)^2 * tau_us * tau_us_0)
    G_3I_array_new[i] = G_3_ad[i]*G_3R_array[i] * ratio_IR
  endfor
;  iplot, G_3_ad*G_3R_array, G_3_ad*G_3I_array, xrange=[0, MAX(G_3_ad*G_3R_array)*1.1], yrange=[0, MAX(G_3_ad*G_3R_array)*1.1], color='red', thick=1, xtitle='G3R [A.U.]', ytitle='G3I [A.U.]', title='STED curve with Fitting', /NO_SAVEPROMPT
;  iplot, G_3_ad*G_3R_array, G_3I_array_new, color='green', linestyle=1, thick=2, /overplot
  
  ;@GJ, 2025/1/26, calculate the filtration kernel
  ;SGJ, 2025/1/27, calculte FWHM and negative
  n_STED = 100.
  STED_factor_array = DBLARR(n_STED) * 0.
  FWHM_array = DBLARR(n_STED) * 0.
  FWHM_array_Gfit = DBLARR(n_STED) * 0.
  neg_array = DBLARR(n_STED) * 0.
  thre_ind = 0.
  FOR i=0, n_STED-1 DO BEGIN
    alpha = i / 80.
    STED_factor_array[i] = alpha / MAX(G_3_ad * G_3I_array) * MAX(G_3_ad * G_3R_array)
    G_3STED = G_3_ad * (G_3R_array - STED_factor_array[i] * G_3I_array)
    IF MIN(G_3STED)-G_3STED[0] GE 0 THEN thre_ind = i
    FWHM_array[i] = 2. * ABS(INTERPOL(x_array[0:n_x/2-1], G_3STED[0:n_x/2-1], MAX(G_3STED)/2.))
    neg_array[i] = ABS(MIN(G_3STED)/MAX(G_3STED))
    yfit_psf_STED = GAUSSFIT(x_array, G_3STED, coeff_STED, NTERMS=4)
    FWHM_array_Gfit[i] = 2*SQRT(2*ALOG(2))*coeff_STED[2]; 
    ;  print, 'Gaussian PSF Fit Result: ', coeff_g[0:4-1]
;    print, 'Gaussian FWHM: ', 2*SQRT(2*ALOG(2))*coeff_g[2], ' mT, ', 2*SQRT(2*ALOG(2))*coeff_g[2] / gradient_x, ' mm'
;    iplot, offset_field_fov_imag_1d, yfit_psf_g, LINSTYLE=0, thick=4, color='red', /OVERPLOT
;    yfit_psf_STED = GAUSSFIT(offset_field_fov_imag_1d, new_3rd_fov_imag_sted_1d, coeff_STED, NTERMS=4)
    ;  print, 'STED PSF Fit Result: ', coeff_STED[0:4-1]
;    print, 'STED FWHM: ', 2*SQRT(2*ALOG(2))*coeff_STED[2], ' mT, ', 2*SQRT(2*ALOG(2))*coeff_STED[2] / gradient_x, ' mm'
;    iplot, offset_field_fov_imag_1d, yfit_psf_sted, LINSTYLE=0, thick=4, color='blue', /OVERPLOT
;    yfit_psf_d = yfit_psf_g - yfit_psf_STED
;    iplot, offset_field_fov_imag_1d, yfit_psf_d, LINSTYLE=0, thick=4, color='green', /OVERPLOT
    IF i EQ 0 THEN BEGIN
;      iplot, x_array, G_3STED, color='red', thick=1, xtitle='x-direction [mm]', ytitle='G3STED [A.U.]', title='G3 STED', /NO_SAVEPROMPT
;      iplot, x_array, yfit_psf_STED, color='red', /overplot
      ;iplot, freq, ABS(FFT(G_3STED)), color='red', thick=1, xrange=[0, 20./fov], xtitle='Freq [1/mm]', ytitle='Amp', title='FFT of G3STED', /NO_SAVEPROMPT
    ENDIF ELSE BEGIN
;      IF i mod 10 EQ 0 THEN BEGIN
;        iplot, x_array, G_3STED, color=350000*(i+1), /overplot
;        iplot, x_array, yfit_psf_STED, color=350000*(i+1), /overplot
;      ENDIF
      ;iplot, freq, ABS(FFT(G_3STED)), color=350000*(i+1), /overplot
    ENDELSE
;    wait, 0.5
;    print, 'i = ', i
  ENDFOR  
  
;  iplot, STED_factor_array * MAX(G_3_ad * G_3I_array) / MAX(G_3_ad * G_3R_array), FWHM_array, color='red', thick=2, xtitle='STED factor', ytitle='FWHM (mm) & Neg', title='G3STED FWHM', /NO_SAVEPROMPT
;  iplot, STED_factor_array * MAX(G_3_ad * G_3I_array) / MAX(G_3_ad * G_3R_array), FWHM_array_Gfit, color='pink', thick=2, /overplot
;  iplot, STED_factor_array * MAX(G_3_ad * G_3I_array) / MAX(G_3_ad * G_3R_array), neg_array/MAX(neg_array)*MAX(FWHM_array), color='green', thick=2, /overplot
  
;  min = MIN(ABS(STED_factor_array * MAX(G_3_ad * G_3I_array) / MAX(G_3_ad * G_3R_array) - 0.618), min_ind)
;  print, '0.618 ind: ', min_ind
;  iplot, G_3_ad*G_3R_array, G_3_ad*G_3I_array, color='green', thick=2, linestyle=0, xrange=[0, MAX(G_3_ad*G_3R_array)*1.1], yrange=[0, MAX(G_3_ad*G_3R_array)*1.1], xtitle='G3 Real', ytitle='G3 Imag', title='STED curves', /NO_SAVEPROMPT
;  iplot, G_3_ad*G_3R_array, G_3I_array_new, color='red', linestyle=2, thick=3, /overplot
;  iplot, G_3_ad*G_3R_array, STED_factor_array[min_ind] * G_3_ad*G_3I_array, color='blue', thick=2, /overplot
  
  ;@GJ, 2025/2/1, output the PSF
  G_3R_psf = G_3_ad * G_3R_array
  G_3I_psf = G_3_ad * STED_factor_array[thre_ind] * G_3I_array
  G_3STED_psf = G_3_ad * (G_3R_array - STED_factor_array[thre_ind] * G_3I_array)
  print, 'Positive ind: ', thre_ind

;  G_3I_psf = G_3_ad * STED_factor_array[min_ind] * G_3I_array
;  G_3STED_psf = G_3_ad * (G_3R_array - STED_factor_array[min_ind] * G_3I_array)
;  FOR i=0, 100 DO BEGIN
;    alpha = i / 30.
;    G_3F = G_3_ad*(G_3R_array / MAX(G_3R_array) - alpha * G_3I_array / MAX(G_3I_array))
;    IF i EQ 0 THEN BEGIN
;;      iplot, x_array, G_3F, color='red', thick=1, xtitle='x-direction [mm]', ytitle='G3STED [A.U.]', title='G3STED', /NO_SAVEPROMPT
;      iplot, freq, ABS(FFT(G_3F)), color='red', thick=1, xrange=[0, 100./fov], xtitle='Freq [1/mm]', ytitle='Amp', title='FFT of G3STED', /NO_SAVEPROMPT
;    ENDIF ELSE BEGIN
;;      iplot, x_array, G_3F, color=350000*(i+1), /overplot
;      iplot, freq, ABS(FFT(G_3F)), color=350000*(i+1), /overplot
;    ENDELSE
;    wait, 0.2
;    print, 'i = ', i
;  ENDFOR
END 

;function image_pSNR_SSIM_RMSE, im_ori, im_rec
;  im_ori_d = DOUBLE(im_ori)
;  im_rec_d = DOUBLE(im_rec)
;  MSE=calMSE(im_ori_d, im_rec_d);72.1795
;  print, 'MSE: ', MSE
;  pSNR=calPSNR(im_ori_d, MSE);68.0337
;  print, 'pSNR: ', pSNR, ' dB'
;  mssim = SSIM(im_ori_d, im_rec_d)
;  print, 'SSIM: ', mssim
;  return, [SQRT(MSE), pSNR, mSSIM]
;end

;@GJ, 2025/2/2, do a radon evaluation based on G3 R, G3 I, and G3 STED
PRO radon_G3STED, wid, phantomListValue, filter_type, parameters
  
  DEVICE, DECOMPOSED=0

  ;Create an image with a ring plus random noise:

  x = (LINDGEN(256,256) MOD 256) - 127.5

  y = (LINDGEN(256,256)/256) - 127.5

  radius = SQRT(x^2 + y^2)

  plate = (radius LT 30 AND radius GT 10); AND (radius LT 100)
  IF phantomListValue EQ 3 THEN array = BYTSCL(SHIFT(plate, 70, 0) + SHIFT(plate, -70, 0))
  plate = (radius LT 30); AND radius GT 10); AND (radius LT 100)
  IF phantomListValue EQ 2 THEN array = BYTSCL(SHIFT(plate, 70, 0) + SHIFT(plate, -70, 0))

;  array = BYTSCL(SHIFT(plate, 70, 0) + SHIFT(plate, -70, 0))

  ;    array = array + RANDOMU(seed,256,256)
  ;    iimage, array, title='original'
  ;Create display window, set graphics properties:
  ;  result = RADON(array, RHO=rho, THETA=theta);, NRHO=256, NTHETA=256)
  IF phantomListValue EQ 6 THEN BEGIN
    fn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\G3_STED\PMB_G.png'
    im = READ_PNG(fn)
    array = BYTSCL(REFORM(im[0,*,*]))
  ENDIF
  IF phantomListValue EQ 8 THEN BEGIN
    fn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\G3_STED\PMB_bar2.png'
    im = READ_PNG(fn)
    array = BYTSCL(REFORM(im[0,*,*]))
  ENDIF
  IF phantomListValue EQ 7 THEN BEGIN
    fn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\G3_STED\PMB_bar3.png'
    im = READ_PNG(fn)
    array = BYTSCL(REFORM(im[0,*,*]))
  ENDIF
  
  ;@GJ, 2025/2/6, do different plot  
  ;shepp_logan_phantom_2d
  IF phantomListValue EQ 1 THEN array = BYTSCL(shepp_logan_phantom_2d(256))
  
  ;get the sample image
  IF phantomListValue EQ 0 THEN BEGIN
    filename2='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\G3_STED\XidianUniversity_ori.png'
    image_temp2 = read_image(filename2)
    array = BYTSCL(CONGRID((REFORM(image_temp2[0,*,*])), 256, 256))
  ENDIF

  IF phantomListValue EQ 4 THEN BEGIN
    p_fn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\G3_STED\concentration_bw.jpg'
    READ_JPEG, p_fn, phantom_image
    array = BYTSCL(CONGRID(REFORM(phantom_image[0,*,*]), 256, 256))
  ENDIF

  IF phantomListValue EQ 5 THEN BEGIN
    p_fn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\G3_STED\catheter_vasculature_thermia02.png'
    phantom_image_thermia = read_image(p_fn)
    array = BYTSCL(CONGRID(phantom_image_thermia, 256, 256))
  ENDIF

  NRHO = parameters.N_points; 237.;137;;137
  NTHETA = parameters.N_angles;27;1. * CEIL(!PI * NRHO/2.);127.;181.;27.;181
;  slider_N_points = widget_info(wid, find_by_uname='N_points')
;  widget_control, slider_N_points, set_value = NRHO
;  parameters.N_points = NRHO
;  slider_N_angles = widget_info(wid, find_by_uname='N_angles')
;  widget_control, slider_N_angles, set_value = NTHETA
;  parameters.N_angles = NTHETA
  FOV = parameters.FOV_phantom;190.;
;  slider_FOV_phantom = widget_info(wid, find_by_uname='FOV_phantom')
;  widget_control, slider_FOV_phantom, set_value = FOV
;  parameters.FOV_phantom = FOV
  array = CONGRID(array, NRHO/2, NRHO/2)
  resolution = FOV/NRHO
  n_x = NRHO/2
  n_y = NRHO/2
  d_x = FOV/n_x
  d_y = FOV/n_y
  ;@GJ, 2025/2/4, define the frequency
  rho = (FINDGEN((NRHO - 1)/2) + 1)
  is_N_even = (NRHO MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, rho, NRHO/2, -NRHO/2 + rho]/(NRHO*resolution) $
  else $
    freq = [0.0, rho, -(NRHO/2 + 1) + rho]/(NRHO*resolution)
;  iimage, array, title='original'
  ;  new_result = CONGRID(result, NTHETA, NRHO)
  ;  new_rho = CONGRID(rho, NRHO)
  ;  new_theta = CONGRID(theta, NTHETA)
  new_result = RADON(array, RHO=rho, THETA=theta, NRHO=NRHO, NTHETA=NTHETA)
  new_result += RANDOMU(seed, NTHETA, NRHO)
;  iimage, new_result, title='Sinogram'

  ;@GJ, 2025/2/1, do the FFL-based signal simulation
  sinogram_R = new_result * 0.
  sinogram_I = new_result * 0.
  sinogram_STED = new_result * 0.
  sinogram_STED_decon = new_result * 0.
  G3RI_equation_calculation, x_array, H_x_array, G_3R_psf, G_3I_psf, G_3STED_psf, parameters
  G_3R_psf_int = CONGRID(G_3R_psf, NRHO, /CENTER, /INTERP)
  G_3I_psf_int = CONGRID(G_3I_psf, NRHO, /CENTER, /INTERP)
  G_3STED_psf_int = CONGRID(G_3STED_psf, NRHO, /CENTER, /INTERP)
  FOR i=0, NTHETA-1 DO BEGIN
    sinogram_R[i, *] = SHIFT(FFT(FFT(REFORM(new_result[i, *])) * FFT(G_3R_psf_int), /INVERSE), NRHO/2)
    sinogram_R[i, *] += RANDOMU(seed, 1, NRHO)
;    CONVOL(REFORM(new_result[i, *]), G_3R_psf[100:150])
;    sinogram_I[i, *] = CONVOL(REFORM(new_result[i, *]), G_3I_psf[100:150])
;    sinogram_STED[i, *] = CONVOL(REFORM(new_result[i, *]), G_3STED_psf[100:150])
    sinogram_I[i, *] = SHIFT(FFT(FFT(REFORM(new_result[i, *])) * FFT(G_3I_psf_int), /INVERSE), NRHO/2)
    sinogram_I[i, *] += RANDOMU(seed, 1, NRHO)
    sinogram_STED[i, *] = SHIFT(FFT(FFT(REFORM(new_result[i, *])) * FFT(G_3STED_psf_int), /INVERSE), NRHO/2)
    sinogram_STED[i, *] += RANDOMU(seed, 1, NRHO)
    
    ;@GJ, 2025/2/3, do the deconvolution of G3STED
    temp_psf_FFT = FFT(G_3STED_psf_int)
    k = MAX(ABS(temp_psf_FFT))/2.
    wiener_filter = CONJ(temp_psf_FFT) / (ABS(temp_psf_FFT)^2 + k)
;    thre_decon = MAX(ABS(temp_psf_FFT))/5.
;    psf_FFT = temp_psf_FFT * (ABS(temp_psf_FFT) GT thre_decon) + (ABS(temp_psf_FFT) LT thre_decon)
;    temp_sinogram_STED_decon = SHIFT(FFT(FFT(REFORM(sinogram_STED[i, *])) * wiener_filter, /INVERSE), NRHO/2)
;    sinogram_STED_decon[i, *] = temp_sinogram_STED_decon; / MAX(temp_sinogram_STED_decon) * MAX(sinogram_STED[i, *])
    ;@GJ, 2025/3/23, do the RL 1d deconvolution
    blurred = reform(sinogram_STED[i, *])
    restored = rl_deconv1d(blurred, G_3STED_psf_int, iterations=100)
    sinogram_STED_decon[i, *] = restored
    
    ;@GJ, 2025/2/2, plot the sinograms
    IF i EQ 0 THEN BEGIN
      iplot, sinogram_R[i, *], color='red', thick=2, yrange=[-MAX(sinogram_R[i, *])/5., MAX(sinogram_R[i, *])*6./5.], xrange=[-NRHO/5, NRHO*6/5], /NO_SAVEPROMPT
      iplot, sinogram_I[i, *], color='green', thick=2, /overplot
      iplot, sinogram_STED[i, *]/MAX(sinogram_STED[i, *])*MAX(sinogram_R[i, *]), color='blue', thick=2, /overplot
      iplot, sinogram_STED_decon[i, *]/MAX(sinogram_STED_decon[i, *])*MAX(sinogram_STED[i, *]), color='cyan', thick=2, /overplot
    ENDIF
  ENDFOR
  
  backproject = RADON(new_result, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_R = RADON(sinogram_R, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_I = RADON(sinogram_I, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_STED = RADON(sinogram_STED, /BACKPROJECT, RHO=rho, THETA=theta)
  backproject_STED_decon = RADON(sinogram_STED_decon, /BACKPROJECT, RHO=rho, THETA=theta)
  ;  iimage, backproject, title='BP'

  filtered_result = new_result * 0.
  filtered_resultR = new_result * 0.
  filtered_resultI = new_result * 0.
  filtered_resultSTED = new_result * 0.
  filtered_resultSTED_decon = new_result * 0.
;  filter_type = 'ram-lak'
;  filter_type = 'cosine'
;  filter_type = 'hamming'
;  filter_type = 'hann'
  FOR i=0, NTHETA-1 DO BEGIN
    filtered_result[i, *] = Filter_FFT(new_result[i, *], filter_type)
    filtered_resultR[i, *] = Filter_FFT(sinogram_R[i, *], filter_type)
    filtered_resultI[i, *] = Filter_FFT(sinogram_I[i, *], filter_type)
    filtered_resultSTED[i, *] = Filter_FFT(sinogram_STED[i, *], filter_type)
    filtered_resultSTED_decon[i, *] = Filter_FFT(sinogram_STED_decon[i, *], filter_type)
  ENDFOR
  filtered_backproject = RADON(filtered_result, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectR = RADON(filtered_resultR, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectI = RADON(filtered_resultI, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectSTED = RADON(filtered_resultSTED, /BACKPROJECT, RHO=rho, THETA=theta)
  filtered_backprojectSTED_decon = RADON(filtered_resultSTED_decon, /BACKPROJECT, RHO=rho, THETA=theta)
  ;  iimage, filtered_backproject, title='filtered BP'
  im_ori = CONGRID(array, NRHO, NRHO)
  im_rec = CONGRID(backproject, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(DOUBLE(BYTSCL(im_ori)), DOUBLE(BYTSCL(im_rec))))
  pSNR = calPSNR(DOUBLE(BYTSCL(im_ori)), RMSE^2)
  quality_backproject_gd = [RMSE, pSNR, mSSIM]
;  quality_backproject_gd = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(backproject_R, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(DOUBLE(BYTSCL(im_ori)), DOUBLE(BYTSCL(im_rec))))
  pSNR = calPSNR(DOUBLE(BYTSCL(im_ori)), RMSE^2)
  quality_backproject_R = [RMSE, pSNR, mSSIM]
;  quality_backproject_R = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(backproject_I, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(DOUBLE(BYTSCL(im_ori)), DOUBLE(BYTSCL(im_rec))))
  pSNR = calPSNR(DOUBLE(BYTSCL(im_ori)), RMSE^2)
  quality_backproject_I = [RMSE, pSNR, mSSIM]
;  quality_backproject_I = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(backproject_STED, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(DOUBLE(BYTSCL(im_ori)), DOUBLE(BYTSCL(im_rec))))
  pSNR = calPSNR(DOUBLE(BYTSCL(im_ori)), RMSE^2)
  quality_backproject_STED = [RMSE, pSNR, mSSIM]
;  quality_backproject_STED = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(backproject_STED_decon, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(DOUBLE(BYTSCL(im_ori)), DOUBLE(BYTSCL(im_rec))))
  pSNR = calPSNR(DOUBLE(BYTSCL(im_ori)), RMSE^2)
  quality_backproject_STED_decon = [RMSE, pSNR, mSSIM]
;  quality_backproject_STED = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(filtered_backproject, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(im_ori, im_rec))
  pSNR = calPSNR(im_ori, RMSE^2)
  quality_result_gd = [RMSE, pSNR, mSSIM]
;  quality_result_gd = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(filtered_backprojectR, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(im_ori, im_rec))
  pSNR = calPSNR(im_ori, RMSE^2)
  quality_result_R = [RMSE, pSNR, mSSIM]
;  quality_result_R = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(filtered_backprojectI, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(im_ori, im_rec))
  pSNR = calPSNR(im_ori, RMSE^2)
  quality_result_I = [RMSE, pSNR, mSSIM]
;  quality_result_I = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(filtered_backprojectSTED, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(im_ori, im_rec))
  pSNR = calPSNR(im_ori, RMSE^2)
  quality_result_STED = [RMSE, pSNR, mSSIM]
; quality_result_STED = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  im_rec = CONGRID(filtered_backprojectSTED_decon, NRHO, NRHO)
  mSSIM = SSIM(im_ori, im_rec)
  RMSE = SQRT(calMSE(im_ori, im_rec))
  pSNR = calPSNR(im_ori, RMSE^2)
  quality_result_STED_decon = [RMSE, pSNR, mSSIM]
; quality_result_STED_decon = image_pSNR_SSIM_RMSE(im_ori, im_rec)
  
;  iplot, [quality_result_gd[0], quality_result_R[0], quality_result_I[0], quality_result_STED[0]], SYM_INDEX=4, SYM_SIZE=2, linestyle=1, title='pSNR'
;  iplot, [quality_result_gd[1], quality_result_R[1], quality_result_I[1], quality_result_STED[1]], SYM_INDEX=4, SYM_SIZE=2, linestyle=1, title='SSIM'
;  iplot, [quality_result_gd[2], quality_result_R[2], quality_result_I[2], quality_result_STED[2]], SYM_INDEX=4, SYM_SIZE=2, linestyle=1, title='MSE'

  NRHO = MIN([120, NRHO])
  DEVICE, DECOMPOSED = 0
  LOADCT, 4
  WINDOW, 0, XSIZE = 6*NRHO, YSIZE = 6*NRHO, TITLE = 'Sinogram & FBP'
  TV, CONGRID(BYTSCL(new_result), NRHO, NRHO), 1
  TV, CONGRID(BYTSCL(sinogram_R), NRHO, NRHO), 2
  TV, CONGRID(BYTSCL(sinogram_I), NRHO, NRHO), 3
  TV, CONGRID(BYTSCL(sinogram_STED), NRHO, NRHO), 4
  TV, CONGRID(BYTSCL(sinogram_STED_decon), NRHO, NRHO), 5
  TV, CONGRID(BYTSCL(filtered_result), NRHO, NRHO), 19
  TV, CONGRID(BYTSCL(filtered_resultR), NRHO, NRHO), 20
  TV, CONGRID(BYTSCL(filtered_resultI), NRHO, NRHO), 21
  TV, CONGRID(BYTSCL(filtered_resultSTED), NRHO, NRHO), 22
  TV, CONGRID(BYTSCL(filtered_resultSTED_decon), NRHO, NRHO), 23
  ;@GJ, 2025/2/2, show the reconstructed images
  LOADCT, 3
  TV, CONGRID(BYTSCL(array), NRHO, NRHO), 6
  TV, CONGRID(BYTSCL(array), NRHO, NRHO), 24
;  TV, CONGRID(BYTSCL(array), NRHO, NRHO), 25
  TV, CONGRID(BYTSCL(backproject), NRHO, NRHO), 7
  TV, CONGRID(BYTSCL(backproject_R), NRHO, NRHO), 8
  TV, CONGRID(BYTSCL(backproject_I), NRHO, NRHO), 9
  TV, CONGRID(BYTSCL(backproject_STED), NRHO, NRHO), 10
  TV, CONGRID(BYTSCL(backproject_STED_decon), NRHO, NRHO), 11
  TV, CONGRID(BYTSCL(filtered_backproject), NRHO, NRHO), 25
  TV, CONGRID(BYTSCL(filtered_backprojectR), NRHO, NRHO), 26
  TV, CONGRID(BYTSCL(filtered_backprojectI), NRHO, NRHO), 27
  TV, CONGRID(BYTSCL(filtered_backprojectSTED), NRHO, NRHO), 28
  TV, CONGRID(BYTSCL(filtered_backprojectSTED_decon), NRHO, NRHO), 29

  ;display the results
  xyouts,0.4*NRHO,3.7*NRHO,'RMSE', color = 'ffffff'x, /DEVICE
  xyouts,0.4*NRHO,3.4*NRHO,'PSNR', color = 'ffffff'x, /DEVICE
  xyouts,0.4*NRHO,3.1*NRHO,'SSIM', color = 'ffffff'x, /DEVICE
  
  xyouts,1.2*NRHO,3.7*NRHO,STRING(quality_backproject_gd[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,1.2*NRHO,3.4*NRHO,STRING(quality_backproject_gd[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,1.2*NRHO,3.1*NRHO,STRING(quality_backproject_gd[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,2.2*NRHO,3.7*NRHO,STRING(quality_backproject_R[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,2.2*NRHO,3.4*NRHO,STRING(quality_backproject_R[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,2.2*NRHO,3.1*NRHO,STRING(quality_backproject_R[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,3.2*NRHO,3.7*NRHO,STRING(quality_backproject_I[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,3.2*NRHO,3.4*NRHO,STRING(quality_backproject_I[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,3.2*NRHO,3.1*NRHO,STRING(quality_backproject_I[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,4.2*NRHO,3.7*NRHO,STRING(quality_backproject_STED[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,4.2*NRHO,3.4*NRHO,STRING(quality_backproject_STED[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,4.2*NRHO,3.1*NRHO,STRING(quality_backproject_STED[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,5.2*NRHO,3.7*NRHO,STRING(quality_backproject_STED_decon[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,5.2*NRHO,3.4*NRHO,STRING(quality_backproject_STED_decon[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,5.2*NRHO,3.1*NRHO,STRING(quality_backproject_STED_decon[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  ;display the results
  xyouts,0.4*NRHO,0.7*NRHO,'RMSE', color = 'ffffff'x, /DEVICE
  xyouts,0.4*NRHO,0.4*NRHO,'PSNR', color = 'ffffff'x, /DEVICE
  xyouts,0.4*NRHO,0.1*NRHO,'SSIM', color = 'ffffff'x, /DEVICE

  xyouts,1.2*NRHO,0.7*NRHO,STRING(quality_result_gd[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,1.2*NRHO,0.4*NRHO,STRING(quality_result_gd[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,1.2*NRHO,0.1*NRHO,STRING(quality_result_gd[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,2.2*NRHO,0.7*NRHO,STRING(quality_result_R[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,2.2*NRHO,0.4*NRHO,STRING(quality_result_R[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,2.2*NRHO,0.1*NRHO,STRING(quality_result_R[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,3.2*NRHO,0.7*NRHO,STRING(quality_result_I[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,3.2*NRHO,0.4*NRHO,STRING(quality_result_I[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,3.2*NRHO,0.1*NRHO,STRING(quality_result_I[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,4.2*NRHO,0.7*NRHO,STRING(quality_result_STED[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,4.2*NRHO,0.4*NRHO,STRING(quality_result_STED[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,4.2*NRHO,0.1*NRHO,STRING(quality_result_STED[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  xyouts,5.2*NRHO,0.7*NRHO,STRING(quality_result_STED_decon[0],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,5.2*NRHO,0.4*NRHO,STRING(quality_result_STED_decon[1],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  xyouts,5.2*NRHO,0.1*NRHO,STRING(quality_result_STED_decon[2],format='(f8.4)'), color = 'ffffff'x, /DEVICE
  
  DEVICE, DECOMPOSED = 0
  WINDOW, 1, XSIZE = 5*NRHO, YSIZE = 2*NRHO, TITLE = 'Sinogram & FBP'
  LOADCT, 4
;  TV, CONGRID(BYTSCL(new_result), NRHO, NRHO), 4
  TV, CONGRID(BYTSCL(sinogram_R), NRHO, NRHO), 1
  TV, CONGRID(BYTSCL(sinogram_I), NRHO, NRHO), 2
  TV, CONGRID(BYTSCL(sinogram_STED), NRHO, NRHO), 3
  TV, CONGRID(BYTSCL(sinogram_STED_decon), NRHO, NRHO), 4
  ;@GJ, 2025/2/2, show the reconstructed images
  LOADCT, 3
  TV, CONGRID(BYTSCL(array), NRHO, NRHO), 5
;  TV, CONGRID(BYTSCL(filtered_backproject), NRHO, NRHO), 8
  TV, CONGRID(BYTSCL(filtered_backprojectR), NRHO, NRHO), 6
  TV, CONGRID(BYTSCL(filtered_backprojectI), NRHO, NRHO), 7
  TV, CONGRID(BYTSCL(filtered_backprojectSTED), NRHO, NRHO), 8
  TV, CONGRID(BYTSCL(filtered_backprojectSTED_decon), NRHO, NRHO), 9

  DEVICE, DECOMPOSED = 1
END


;@GJ, 2024/9/10, calculate the magnetic field dependent relaxation time
;@GJ, 2025/4/7, correct the perp H as DC^2 + AC^2
FUNCTION fields_brownian_time_calc, tau_B, beta, H_AC, H_DC_par, H_DC_perp
  xi_H_AC = beta * H_AC
  coeff_H_AC = 1.;/SQRT(1. + 0.126 * xi_H_AC^1.72)

  xi_H_DC_par = beta * H_DC_par
  IF xi_H_DC_par NE 0 THEN BEGIN
    Lprime_H_DC_par = 1./(xi_H_DC_par^2) - 1./(sinh(xi_H_DC_par)^2)
    L_H_DC_par = 1./tanh(xi_H_DC_par) - 1./xi_H_DC_par
    coeff_H_DC_par = xi_H_DC_par * Lprime_H_DC_par / L_H_DC_par
  ENDIF ELSE BEGIN
    coeff_H_DC_par = 1.
  ENDELSE

  xi_H_DC_perp = beta * H_DC_perp;SQRT(H_DC_perp^2 + H_AC^2)
  IF xi_H_DC_perp NE 0 THEN BEGIN
    L_H_DC_perp = 1./tanh(xi_H_DC_perp) - 1./xi_H_DC_perp
    coeff_H_DC_perp = 2. * L_H_DC_perp / (xi_H_DC_perp - L_H_DC_perp)
  ENDIF ELSE BEGIN
    coeff_H_DC_perp = 1.
  ENDELSE

  tau_B_H = tau_B * coeff_H_AC * coeff_H_DC_perp; * coeff_H_DC_par
  RETURN, tau_B_H; brwonian_relaxation_time
END

;@GJ, 2023/12/10, including the orghogonal Gradient field
pro sin_simulation,wid,A,A_DC,A_G,particle_size,f,tau,gradient_field,signal_ratio,RESI_ratio,resolution_exp_signal,resolution_signal_new_sum, Dpp_reso, drawXSize, drawYSize, noise_level, filename, N_points, D_points, fov_phantom, n_angles, angle, shift_perc
  
;  print, 'wid:', wid
  tau_ms = tau/1000.
  
  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  N = 5000.;
  Time = findgen(5000)/1000.     ;(1  ms)
  N_period = 1./ f / (Time[1]-Time[0])
  Time_period = 1. / f
  N_cycles = N / N_period
  ;@GJ, 2022/12/14, changing to offset
  H = -A*cos(2.*!PI*f*Time)+A_DC*0.1 ;(mT)
  H_0 = -A*cos(2.*!PI*f*Time)
;  H = A*sin(2.*!PI*f*Time)+A_DC*0.1 ;(mT)
;  H_0 = A*sin(2.*!PI*f*Time)

  ;@GJ, 2023/12/10, Orthogonal field is not zero
  IF A_G GT 0 THEN BEGIN
    H_total = SIGNUM(H) * SQRT(H^2 + (A_G*0.1)^2)
  ENDIF

  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
;  LangzwX = (FINDGEN(401)-200.)/10.
;  Langzw = Msat_kAm*(1./TANH(beta*LangzwX) - 1./(beta*LangzwX))
;  Langzw[200] = 0

  M = Msat_kAm*(1/tanh(beta*H) - 1/(beta*H)) ;(kA/m)
  IF ~FINITE(M[0]) THEN M[0] = M[N_period]
  M_0 = Msat_kAm*(1/tanh(beta*H_0) - 1/(beta*H_0)) ;(kA/m)
  IF ~FINITE(M_0[0]) THEN M_0[0] = M_0[N_period]
  ;@GJ, 2023/12/10, Orthogonal field is not zero
  IF A_G GT 0 THEN BEGIN
    M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
    ;do vector splitting
    ;@GJ, 2023/9/23, simple calculation with the signal intensity
    M = ABS(M_total) * H/ABS(H_total)
  ENDIF
  
  R = (1/tau_ms)*exp(-Time/tau_ms)
  n = N_elements(M)
  signal = M[1:*] - M[0:n-1]
  signal_0 = M_0[1:*] - M_0[0:n-1]
  ;@GJ, 2022/12/14, removing the NaN points
  index_nan = where(~FINITE(signal), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(signal, Time, Time[index_nan], /NAN)
    signal[index_nan] = RESULT
  ENDIF
  index_nan = where(~FINITE(signal_0), count)
  IF count GT 0 THEN BEGIN
    RESULT = INTERPOL(signal_0, Time, Time[index_nan], /NAN)
    signal_0[index_nan] = RESULT
  ENDIF
  ;@GJ, 2025/1/12, calculate the tau_ms
  viscosity = 1.; mPa.s
  T_p = 25.; degree
  eta = viscosity / 1000. ; Pa.s
  D_h = particle_size ; in nm
;  print, 'D_h: ', D_h, ' nm'
  D_h *= 1.e-9 ; m
  V_h = 1./6. * !PI * D_h^3
  k_B = 1.380469e-23; J/K
  T_p_kelvin = 273.15 + T_p ; in Kelvin
  tau_B = 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
  print, 'tau_us_0: ', tau_B
  tau_us = fields_brownian_time_calc(tau_B, beta, A, A_DC*0.1, A_G*0.1)
  tau_us_0 = fields_brownian_time_calc(tau_B, beta, A, 0., 0.)
  print, 'tau_us with offset: ', tau_us
  slider_tau_sin = widget_info(wid, find_by_uname='slider_tao_sin')
  widget_control, slider_tau_sin, set_value = tau_us
  tau_ms = tau_us / 1000.
  tau_ms_0 = tau_us_0 / 1000.
  widget_control, wid, get_uvalue = parameters
  parameters.tau = tau_us;
  widget_control, wid, set_uvalue = parameters
  kernel = Debye_Kernel_sin(tau_ms, Time[0:100])
  kernel_0 = Debye_Kernel_sin(tau_ms_0, Time[0:100])
  
  
  M = convol([dblarr(n_elements(kernel)/2)*0, M, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
  M = M[n_elements(kernel)/2 : n_elements(M) - n_elements(kernel)/2 - 1]
  M_0 = convol([dblarr(n_elements(kernel_0)/2)*0, M_0, dblarr(n_elements(kernel_0)/2)*0], kernel_0, /NORMALIZE)
  M_0 = M_0[n_elements(kernel_0)/2 : n_elements(M_0) - n_elements(kernel_0)/2 - 1]

  signal_new = convol([dblarr(n_elements(kernel)/2)*0, signal, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
  signal_new = signal_new[n_elements(kernel)/2 : n_elements(signal) + n_elements(kernel)/2]
  signal_new_0 = convol([dblarr(n_elements(kernel_0)/2)*0, signal_0, dblarr(n_elements(kernel_0)/2)*0], kernel_0, /NORMALIZE)
  signal_new_0 = signal_new_0[n_elements(kernel_0)/2 : n_elements(signal_0) + n_elements(kernel_0)/2]

  IF wid GT 0 THEN BEGIN
    !P.FONT = 2
    
    cgLoadCT, 34
    ;ͨ��widget_info�����������unameΪ��ʶ���в��ң��������id
    draw1 = widget_info(wid, find_by_uname='draw1')
    widget_control, draw1, get_value = drawWindow1
    drawSize = widget_info(draw1, /GEOMETRY)
    ;wset����draw1���Ϊ�滭����
    wset, drawWindow1
    ;cgScaleVector���ڽ������������Ԫ������Ϊ����ݷ�Χ��
    colors = cgScaleVector(Findgen(N_Elements(H)), Min(H), Max(H))
    ;Value_Locate���ڽ�  data ��  colors(��������)��ʾ�� ��ͬ��dataֵ�в�ͬ��colors
    elevColors = Value_Locate(colors, H)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    plot, Time, H, /NoData, XRANGE=[0, 2.*Time_period], background = 'ffffff'x, color = '000000'x, xtitle='time (ms)'
    xyouts,5,drawSize.ysize*0.95,'magnetic field (mT)', color = '000000'x, /DEVICE
    FOR j=0,2.*N_Elements(H)/N_cycles $
      DO cgPlotS, [Time[j], Time[j+1]], [H[j], H[j+1]], Color=elevColors[j+1]
    
    draw2 = widget_info(wid, find_by_uname='draw2')
    widget_control, draw2, get_value = drawWindow2
    wset, drawWindow2
  ENDIF

  H_2 = H[N_Elements(H)/2 : N_Elements(H)-1]  ; ��ȡһ�� ��ֹ��һ�����ھ����ֶ�����
  M_2 = M[N_Elements(M)/2 : N_Elements(M)-1]
  
  H_2_0 = H_0[N_Elements(H_0)/2 : N_Elements(H_0)-1]  ; ��ȡһ�� ��ֹ��һ�����ھ����ֶ�����
  M_2_0 = M_0[N_Elements(M_0)/2 : N_Elements(M_0)-1]
  
  IF wid GT 0 THEN BEGIN
    colors = cgScaleVector(Findgen(N_Elements(M_2)), Min(M_2), Max(M_2))
    elevColors = Value_Locate(colors, M_2)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    plot, H_2, M_2, /NoData, background = 'ffffff'x, color = '000000'x, xtitle='magnetic field (mT)'
    xyouts,5,drawSize.ysize*0.95,'magnetization (kA/m)', color = '000000'x, /DEVICE
    FOR j=0,N_Elements(H_2)-2 $
      DO cgPlotS, [H_2[j], H_2[j+1]], [M_2[j], M_2[j+1]], Color=elevColors[j+1]

;    draw3 = widget_info(wid, find_by_uname='draw3')
;    widget_control, draw3, get_value = drawWindow3
;    wset, drawWindow3
;    colors = cgScaleVector(Findgen(N_Elements(M)), Min(M), Max(M))
;    elevColors = Value_Locate(colors, M)
;    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
;    plot, Time, M, /NoData, background = 'ffffff'x, color = '000000'x, xtitle='time (ms)'
;    xyouts,5,drawSize.ysize*0.95,'magnetization (kA/m)', color = '000000'x, /DEVICE
;    FOR j=0,N_Elements(M)-2 $
;      DO cgPlotS, [Time[j], Time[j+1]], [M[j], M[j+1]], Color=elevColors[j+1]

    draw3 = widget_info(wid, find_by_uname='draw3')
    widget_control, draw3, get_value = drawWindow3
    wset, drawWindow3
    colors = cgScaleVector(Findgen(N_Elements(signal_new)), Min(signal_new), Max(signal_new))
    elevColors = Value_Locate(colors, signal_new)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    plot, Time, signal_new, /NoData, XRANGE=[0, 2.*Time_period], YRANGE=[-1.1*MAX(signal_new_0), 1.1*MAX(signal_new_0)], background = 'ffffff'x, color = '000000'x, xtitle='time (ms)'
    xyouts,5,drawSize.ysize*0.95,'signal (kA/m/ms)', color = '000000'x, /DEVICE
    FOR j=0,2.*N_Elements(signal_new)/N_cycles $
      DO cgPlotS, [Time[j], Time[j+1]], [signal_new[j], signal_new[j+1]], Color=elevColors[j+1]

  ENDIF
  
  ;frequency = findgen(N_Elements(u))
  ;@GJ, 2022/12/14, redefining the frequency
  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  t_ele = N_Elements(signal_new)
  delta_t = (Time[1]-Time[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
  else $
    frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
    
  u = ABS(FFT(signal_new))
  u_real = REAL_PART(FFT(signal_new))
  u_imag = IMAGINARY(FFT(signal_new))
  u_0 = ABS(FFT(signal_new_0))
  u_real_0 = REAL_PART(FFT(signal_new_0))
  u_imag_0 = IMAGINARY(FFT(signal_new_0))
  base_freq = MIN(ABS(frequency-1.0), base_ind)
  signal_FFT = FFT(signal_new)
  signal_FFT_real = REAL_PART(signal_FFT)
  signal_FFT_imag = IMAGINARY(signal_FFT)
  signal_FFT_0 = FFT(signal_new_0)
  angle_0 = ATAN(IMAGINARY(signal_FFT_0[3.*base_ind]), REAL_PART(signal_FFT_0[3.*base_ind]), /phase) - !PI/2.
  u_real_d = -u_real * sin(angle_0) + u_imag * cos(angle_0)
  u_imag_d = -u_imag * sin(angle_0) - u_real * cos(angle_0)
  signal_FFT_real_d = -signal_FFT_real * sin(angle_0) + signal_FFT_imag * cos(angle_0)
  signal_FFT_imag_d = -signal_FFT_imag * sin(angle_0) - signal_FFT_real * cos(angle_0)
  signal_FFT = COMPLEX(signal_FFT_real_d, signal_FFT_imag_d)
  
  freqRangeMax = 10.
  min_freq = MIN(ABS(frequency-freqRangeMax), freqRangeMax_Ind)
  IF wid GT 0 THEN BEGIN
    colors = cgScaleVector(Findgen(N_Elements(u)), Min(u), Max(u))
    elevColors = Value_Locate(colors, u)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))

    colors_real_d = cgScaleVector(Findgen(N_Elements(u_real_d)), Min(u_real_d), Max(u_real_d))
    elevColors_real_d = Value_Locate(colors_real_d, u_real_d)
    elevColors_real_d = Byte(Round(cgScaleVector(elevColors_real_d, 0, 255)))

    colors_imag_d = cgScaleVector(Findgen(N_Elements(u_imag_d)), Min(u_imag_d), Max(u_imag_d))
    elevColors_imag_d = Value_Locate(colors_imag_d, u_imag_d)
    elevColors_imag_d = Byte(Round(cgScaleVector(elevColors_imag_d, 0, 255)))

;    plot, frequency/f, u, /NoData, XRANGE=[0, freqRangeMax],background = 'ffffff'x, color = '000000'x, xtitle='n of f'
;    xyouts,20,drawSize.ysize*0.95,'u', color = '000000'x, /DEVICE
;    FOR j=0,freqRangeMax_Ind do $
;      cgPlotS, [frequency[j], frequency[j+1]], [u[j], u[j+1]], Color=elevColors[j+1]
;    cgPlotS, [frequency[0], frequency[freqRangeMax_Ind]], [0, 0], Color=black
    
    ;draw the real part
    draw4a = widget_info(wid, find_by_uname='draw4a')
    widget_control, draw4a, get_value = drawWindow4a
    wset, drawWindow4a
    plot, frequency, u_real_d, XRANGE=[0, freqRangeMax], YRANGE=[-3, 2], background = 'ffffff'x, color = '000000'x, xtitle='n of f'
    xyouts,20,drawSize.ysize*0.4,'Harmonic Real', color = '000000'x, /DEVICE
    for j=0,freqRangeMax_Ind do $
      cgPlotS, [frequency[j], frequency[j+1]], [u_real_d[j], u_real_d[j+1]], linestyle=0, thick=2, Color=elevColors_real_d[j+1]
    cgPlotS, [frequency[0], frequency[freqRangeMax_Ind]], [0, 0], Color=black
    
    ;draw the imaginary part
    draw4b = widget_info(wid, find_by_uname='draw4b')
    widget_control, draw4b, get_value = drawWindow4b
    wset, drawWindow4b
    plot, frequency, u_imag_d, /NoData, XRANGE=[0, freqRangeMax], background = 'ffffff'x, color = '000000'x, xtitle='n of f'
    xyouts,20,drawSize.ysize*0.4,'Harmonic Imaginary', color = '000000'x, /DEVICE
    for j=0,freqRangeMax_Ind-2 do $
      cgPlotS, [frequency[j], frequency[j+1]], [u_imag_d[j], u_imag_d[j+1]], linestyle=0, thick=2, Color=elevColors_imag_d[j+1]
    cgPlotS, [frequency[0], frequency[freqRangeMax_Ind]], [0, 0], Color=black
  ENDIF
  
  ;@GJ, 2025/1/14, update the a3_array
  widget_control, wid, get_uvalue = parameters
  a3_array = [*(parameters.a3_array), signal_FFT[3.*base_ind]]
  parameters.a3_array = PTR_NEW(a3_array)
  widget_control, wid, set_uvalue = parameters

  
  signal_FFT[base_ind] *= 0.
  signal_new_no_base = REAL_PART(FFT(signal_FFT, /inverse))
  
  IF wid GT 0 THEN BEGIN
    draw6 = widget_info(wid, find_by_uname='draw6')
    widget_control, draw6, get_value = drawWindow6
    wset, drawWindow6
    colors = cgScaleVector(Findgen(N_Elements(signal_new_no_base)), Min(signal_new_no_base), Max(signal_new_no_base))
    elevColors = Value_Locate(colors, signal_new_no_base)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    plot, H, signal_new_no_base,/NoData, background = 'ffffff'x, color = '000000'x, xtitle='magnetic field (mT)'
    xyouts,5,drawSize.ysize*0.95,'signal (kA/m/ms)', color = '000000'x, /DEVICE
    FOR j=0,N_Elements(signal_new_no_base)-2 $
      DO cgPlotS, [H[j], H[j+1]], [signal_new_no_base[j], signal_new_no_base[j+1]], Color=elevColors[j+1]
  ENDIF ELSE BEGIN
    ;plot the signal with and without base
    plot1 = PLOT(Time, signal_new, 'r2', xtitle='time (ms)', ytitle='signal (kA/m/ms)', title='Signal w/o Base', NAME='Signal')
    ; Display the second plot.
    plot2 = PLOT(Time, signal_new_no_base, /OVERPLOT, 'b2', NAME='Signal no Base')
    leg = LEGEND(TARGET=[plot1,plot2], POSITION=[185,0.9], /DATA, /AUTO_TEXT_COLOR)
  ENDELSE
  
  
  half_cycle = 1. / (2. * f) / (Time[1]-Time[0])
  ;@GJ, 2023/9/15, signal original
  signal_new_rev = signal_new[0:half_cycle-1] * 0.
  FOR k=0, half_cycle-1 DO signal_new_rev[k] = INTERPOL(signal_new[half_cycle:2.*half_cycle-1], H[half_cycle:2.*half_cycle-1], H[k])
  signal_new_diff = signal_new[0:half_cycle-1] + signal_new_rev
  signal_new_sum = signal_new[0:half_cycle-1] - signal_new_rev
  ;@GJ, 2023/9/15, signal new
  signal_new_rev_no_base = signal_new_no_base[0:half_cycle-1] * 0.
  FOR k=0, half_cycle-1 DO signal_new_rev_no_base[k] = INTERPOL(signal_new_no_base[half_cycle:2.*half_cycle-1], H[half_cycle:2.*half_cycle-1], H[k])
  signal_new_diff_no_base_no_shift = signal_new_no_base[0:half_cycle-1] + signal_new_rev_no_base
;  signal_new_sum_no_base = signal_new_no_base[0:half_cycle-1] - signal_new_rev_no_base
  
  ;@GJ, 2023/9/12, calculate Dpp
  min_signal_new_diff = MIN(signal_new_diff_no_base_no_shift, min_ind)
  max_signal_new_diff = MAX(signal_new_diff_no_base_no_shift, max_ind)
  Dpp_H = ABS(H[max_ind] - H[min_ind]) ; mT
  Dpp_reso = Dpp_H/(gradient_field*0.1) ; mm
;  print, 'Dpp (RESI):'+STRING(Dpp_H,format='(f6.2)')+' mT, '+STRING(Dpp_reso,format='(f6.2)')+' mm'
  
  Dpp_ind = ABS(max_ind - min_ind)
;  print, 'Dpp (index):'+STRING(Dpp_ind)
  shift_ele = FLOAT(shift_perc) / 200. * Dpp_ind

  shift_ind_array = FINDGEN(Dpp_ind*2);-Dpp_ind/2
  diff_abs_max_array = DBLARR(Dpp_ind*2) * 0.
  sum_max_array = DBLARR(Dpp_ind*2) * 0.
  Dpp_array = DBLARR(Dpp_ind*2) * 0.
  diff_abs_array = DBLARR(half_cycle, Dpp_ind*2) * 0.
  sum_array = DBLARR(half_cycle, Dpp_ind*2) * 0.
  product_array = DBLARR(half_cycle, Dpp_ind*2) * 0.
  FOR i=0, 2*Dpp_ind-1 DO BEGIN
    ;@GJ, 2023/9/27, fixing the positive part, only move the negative part
    signal_positive = SHIFT(signal_new[0:half_cycle-1], -shift_ind_array[i]/2)
    signal_negative = SHIFT(signal_new_rev, shift_ind_array[i]/2)
    signal_diff_abs_temp = ABS(signal_positive + signal_negative)
    product_signal = signal_positive*(-signal_negative)
    product_array[*,i] = product_signal
    max_diff = MAX(signal_positive + signal_negative, max_diff_ind)
    min_diff = MIN(signal_positive + signal_negative, min_diff_ind)
    Dpp_array[i] = ABS(H[max_diff_ind] - H[min_diff_ind])/(gradient_field*0.1) ; mm
    signal_sum_temp = signal_positive - signal_negative
    diff_abs_array[*,i] = signal_diff_abs_temp
    sum_array[*,i] = signal_sum_temp
    diff_abs_max_array[i] = MAX(signal_diff_abs_temp)
    sum_max_array[i] = MAX(signal_sum_temp)
  ENDFOR
  ;@GJ, 2023/9/16, plot the shift results
;  iplot, shift_ind_array, diff_abs_max_array, title='max', xtitle='shift index [# pixels]', ytitle='diff abs max value'
;  iplot, shift_ind_array, sum_max_array, /overplot
;  iplot, shift_ind_array, sum_max_array, title='max', xtitle='shift index [# pixels]', ytitle='sum max value'
;  iplot, shift_ind_array, Dpp_array, title='Dpp', xtitle='shift index [# pixels]', ytitle='Dpp [mm]'
;  isurface, diff_abs_array, Time[0:half_cycle-1], shift_ind_array, title='diff abs curves', xtitle='Time [ms]', ytitle='shift index [# pixels]', ztitle='diff abs'
;  isurface, sum_array, Time[0:half_cycle-1], shift_ind_array, title='sum curves', xtitle='Time [ms]', ytitle='shift index [# pixels]', ztitle='sum'
;  
  ;@GJ, 2023/9/13, TAURUS analysis
  t_ele_TAURUS = N_ELEMENTS(signal_new_no_base[0:half_cycle-1])
  X = (FINDGEN((t_ele_TAURUS - 1)/2) + 1)
  is_N_even = (t_ele_TAURUS MOD 2) EQ 0
  if (is_N_even) then $
    frequency_TAURUS = [0.0, X, t_ele_TAURUS/2, -t_ele_TAURUS/2 + X]/(t_ele_TAURUS*delta_t) $
  else $
    frequency_TAURUS = [0.0, X, -(t_ele_TAURUS/2 + 1) + X]/(t_ele_TAURUS*delta_t)
  S_pos = FFT(signal_new_no_base[0:half_cycle-1])
  S_neg = FFT(-signal_new_rev_no_base)
  tao_f = -IMAGINARY((CONJ(S_pos)+S_neg)/(2.*!PI*frequency_TAURUS*(CONJ(S_pos)-S_neg)))
  tao_TAURUS = TOTAL(ABS(S_pos[1:*])*tao_f[1:*])/TOTAL(ABS(S_pos[1:*]))
;  print, 'taurus: ' +STRING(tao_TAURUS,format='(f6.2)')+' us'
  
  ;@GJ, 2023/10/1, calculate signal with shift
  signal_new_diff_no_base = SHIFT(signal_new[0:half_cycle-1], -shift_ele) + SHIFT(signal_new_rev, shift_ele)
  signal_new_sum_no_base = SHIFT(signal_new[0:half_cycle-1], -shift_ele) - SHIFT(signal_new_rev, shift_ele)
  signal_ratio = Max(signal_new_diff_no_base)/Max(signal_new_no_base) ; max signal ratio
;  print, 'Signal ratio (RESI/STD):'+STRING(signal_ratio,format='(f6.2)')
  
  ;@GJ, 2023/10/10, calculate the signal slope
  ;1st, a bar
;  signal_1st_temp = signal_new_sum_no_base * 0.
;  FOR i=0,3 DO signal_1st_temp += SHIFT(signal_new_sum_no_base, -i) + SHIFT(signal_new_sum_no_base, i)
;  signal_1st_temp /= i
;  signal_2nd_temp = SHIFT(signal_new_sum_no_base, -0) + SHIFT(signal_new_sum_no_base, 0)
;  signal_3rd_temp = SHIFT(signal_new_sum_no_base, -2) + SHIFT(signal_new_sum_no_base, 2)
;  iplot, H[0:half_cycle-1], signal_1st_temp, PSYM=4, title='TWIN signal curve'
;;  iplot, H[0:half_cycle-1], signal_2nd_temp, color='blue', /overplot
;  iplot, H[0:half_cycle-1], signal_3rd_temp, color='red', /overplot
;  
;  iplot, H[0:half_cycle-2], signal_1st_temp[1:half_cycle-1]-signal_1st_temp[0:half_cycle-2], PSYM=4, title='TWIN signal gradient'
;;  iplot, H[0:half_cycle-2], signal_2nd_temp[1:half_cycle-1]-signal_2nd_temp[0:half_cycle-2], color='blue', /overplot
;  iplot, H[0:half_cycle-2], signal_3rd_temp[1:half_cycle-1]-signal_3rd_temp[0:half_cycle-2], color='red', /overplot
;  ;
;  ;@GJ, 2023/10/10, calculate the signal slope
;  ;1st, a bar
;  signaldiff_1st_temp = signal_new_diff_no_base * 0.
;  FOR i=0,3 DO signaldiff_1st_temp += SHIFT(signal_new_diff_no_base, -i) + SHIFT(signal_new_diff_no_base, i)
;  signaldiff_1st_temp /= i
;  signaldiff_2nd_temp = SHIFT(signal_new_diff_no_base, -0) + SHIFT(signal_new_diff_no_base, 0)
;  signaldiff_3rd_temp = SHIFT(signal_new_diff_no_base, -2) + SHIFT(signal_new_diff_no_base, 2)
;  iplot, H[0:half_cycle-1], signaldiff_1st_temp, PSYM=4, title='TWIN signal curve'
;;  iplot, H[0:half_cycle-1], signaldiff_2nd_temp, color='blue', /overplot
;  iplot, H[0:half_cycle-1], signaldiff_3rd_temp, color='red', /overplot
;
;  iplot, H[0:half_cycle-2], signaldiff_1st_temp[1:half_cycle-1]-signaldiff_1st_temp[0:half_cycle-2], PSYM=4, title='TWIN signal gradient'
;;  iplot, H[0:half_cycle-2], signaldiff_2nd_temp[1:half_cycle-1]-signaldiff_2nd_temp[0:half_cycle-2], color='blue', /overplot
;  iplot, H[0:half_cycle-2], signaldiff_3rd_temp[1:half_cycle-1]-signaldiff_3rd_temp[0:half_cycle-2], color='red', /overplot
  
  IF wid GT 0 THEN BEGIN
    colors = cgScaleVector(Findgen(N_Elements(signal_new_diff_no_base)), Min(signal_new_diff_no_base), Max(signal_new_diff_no_base))
    elevColors = Value_Locate(colors, signal_new_diff_no_base)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    FOR j=0, half_cycle-2 DO BEGIN
      cgPlotS, [H[j], H[j+1]], [signal_new_diff_no_base[j], signal_new_diff_no_base[j+1]], Color=elevColors[j+1]
    ENDFOR
;    xyouts,110,drawSize.ysize*0.38,'Signal ratio (RESI/STD):'+STRING(signal_ratio,format='(f6.2)'), color = '000000'x, /DEVICE
  ENDIF  
    
  ;@GJ, 2023/9/12, calculate signal without shift
  signal_new_diff_no_shift = signal_new[0:half_cycle-1] + signal_new_rev
  signal_new_sum_no_shift = signal_new[0:half_cycle-1] - signal_new_rev
  abs_signal_new_diff_no_shift = ABS(signal_new_diff_no_base_no_shift)
  relax_signal_no_shift = signal_new_sum_no_shift - abs_signal_new_diff_no_shift
  relax_signal_no_shift[WHERE(relax_signal_no_shift LT 0)] = 0.
  
  ;@GJ, 2023/10/1, calculate the enhanced resolution
  ;@GJ, 2023/9/12, do simple subtraction based signal_new_sum_no_base
  abs_signal_new_diff = ABS(signal_new_diff_no_base)
  ;@GJ, 2023/10/1, calculate the nearest max
  FOR k=2, half_cycle/2-2 DO BEGIN
    IF abs_signal_new_diff[half_cycle/2+k]-abs_signal_new_diff[half_cycle/2+k-1] LT 0 THEN BEGIN
      first_max = abs_signal_new_diff[half_cycle/2+k-1]
      BREAK
    ENDIF
  ENDFOR
  
  ;@GJ, 2023/10/1
;  relax_signal = relax_signal_no_shift - abs_signal_new_diff / ABS(first_max) * MAX(relax_signal_no_shift)
  ;@GJ, 2023/10/11, do the scaling and subtraction
  relax_signal = signal_new_sum_no_base - abs_signal_new_diff / ABS(first_max) * MAX(signal_new_sum_no_base)
  relax_signal[WHERE(relax_signal LT 0)] = 0.
  IF wid GT 0 THEN BEGIN
    cgPlotS, H[0:half_cycle-1], relax_signal/MAX(relax_signal)*Max(signal_new_no_base)*1.2, color='red'
    cgPlotS, H[0:half_cycle-1], signal_new_sum_no_base/MAX(signal_new_sum_no_base)*Max(signal_new_no_base)*1.2, color='blue'
 ;   cgPlotS, H[0:half_cycle-1], signal_new_sum_no_base/MAX(signal_new_sum_no_base)*Max(signal_new_no_base), color='black'
  ENDIF
  
  ;@GJ, 2023/8/25, calculate FWHM of the signal
  max_relax_signal = MAX(relax_signal, max_relax_signal_ind)
  relax_signal_lower_H = INTERPOL(H[0:max_relax_signal_ind], relax_signal[0:max_relax_signal_ind], 0.5*max_relax_signal)
  relax_signal_higher_H = INTERPOL(H[max_relax_signal_ind:half_cycle-1], relax_signal[max_relax_signal_ind:half_cycle-1], 0.5*max_relax_signal)
  FWHM_relax_signal_H = relax_signal_higher_H - relax_signal_lower_H
  resolution_relax_signal = FWHM_relax_signal_H/(gradient_field*0.1)
  
  ;@GJ, 2023/8/25, calculate FWHM of the signal without delay
  max_signal = MAX(signal_new_sum_no_base, max_signal_sum_ind)
  signal_lower_H = INTERPOL(H[0:max_signal_sum_ind], signal_new_sum_no_base[0:max_signal_sum_ind], 0.5*MAX(signal_new_sum_no_base))
  signal_higher_H = INTERPOL(H[max_signal_sum_ind:half_cycle-1], signal_new_sum_no_base[max_signal_sum_ind:half_cycle-1], 0.5*MAX(signal_new_sum_no_base))
  FWHM_signal_sum_H = signal_higher_H - signal_lower_H
  resolution_signal_sum = FWHM_signal_sum_H/(gradient_field*0.1)
  
  ;@GJ, 2023/8/25, calculate FWHM of the signal
  max_signal_new_sum = MAX(signal_new_sum_no_base, max_signal_new_sum_ind)
  signal_new_lower_H = INTERPOL(H[0:max_signal_new_sum_ind], signal_new_sum_no_base[0:max_signal_new_sum_ind], 0.5*MAX(signal_new_sum_no_base))
  signal_new_higher_H = INTERPOL(H[max_signal_new_sum_ind:half_cycle-1], signal_new_sum_no_base[max_signal_new_sum_ind:half_cycle-1], 0.5*MAX(signal_new_sum))
  FWHM_signal_new_sum_H = signal_new_higher_H - signal_new_lower_H
  resolution_signal_new_sum = FWHM_signal_new_sum_H/(gradient_field*0.1)
  RESI_ratio = FWHM_signal_new_sum_H/FWHM_relax_signal_H ;RESI ratio about resolution enhancement ratio
;  print, 'FWHM (STD):'+STRING(FWHM_signal_new_sum_H,format='(f6.2)')+' mT, '+STRING(resolution_signal_new_sum,format='(f6.2)')+' mm'
;  print, 'FWHM (M-H):'+STRING(FWHM_signal_sum_H,format='(f6.2)')+' mT, '+STRING(resolution_signal_sum,format='(f6.2)')+' mm'
;  print, 'FWHM (RESI):'+STRING(FWHM_relax_signal_H,format='(f6.2)')+' mT, '+STRING(resolution_relax_signal,format='(f6.2)')+' mm'
;  print, 'RESI Ratio:'+STRING(FWHM_signal_new_sum_H/FWHM_relax_signal_H,format='(f6.2)')
  
;  IF wid GT 0 THEN BEGIN
;    xyouts,110,drawSize.ysize*0.31,'FWHM (STD):'+STRING(FWHM_signal_new_sum_H,format='(f6.2)')+' mT, '+STRING(resolution_signal_new_sum,format='(f6.2)')+' mm', color = '000000'x, /DEVICE
;    xyouts,110,drawSize.ysize*0.24,'FWHM (RESI):'+STRING(FWHM_relax_signal_H,format='(f6.2)')+' mT, '+STRING(resolution_relax_signal,format='(f6.2)')+' mm, RESI:'+STRING(RESI_ratio,format='(f6.2)'), color = '000000'x, /DEVICE
;  ENDIF 
  
  
  IF wid GT 0 THEN BEGIN
;    xyouts,75,drawSize.ysize*0.41,'FWHM (STD):'+STRING(resolution_signal_new_sum,format='(f6.2)')+' mm', color = 'FF0000'x, /DEVICE
;    xyouts,75,drawSize.ysize*0.34,'FWHM (M-H):'+STRING(resolution_signal_sum,format='(f6.2)')+' mm', color = '000000'x, /DEVICE
    xyouts,75,drawSize.ysize*0.34,'FWHM (STD):'+STRING(resolution_signal_new_sum,format='(f6.2)')+' mm', color = 'FF0000'x, /DEVICE
    xyouts,75,drawSize.ysize*0.27,'FWHM (RESI):'+STRING(resolution_relax_signal,format='(f6.2)')+' mm', color = '0000FF'x, /DEVICE
    xyouts,75,drawSize.ysize*0.20,'RESI (RESI/STD):'+STRING(RESI_ratio,format='(f6.2)'), color = '0000FF'x, /DEVICE
    
    draw6 = widget_info(wid, find_by_uname='draw6')
    widget_control, draw6, get_value = drawWindow6
    wset, drawWindow6
    colors = cgScaleVector(Findgen(N_Elements(M)), Min(M), Max(M))
    elevColors = Value_Locate(colors, M)
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    plot, Time, M, /NoData, XRANGE=[0, 2.*Time_period], background = 'ffffff'x, color = '000000'x, xtitle='time (ms)'
    xyouts,5,drawSize.ysize*0.95,'magnetization (kA/m)', color = '000000'x, /DEVICE
    FOR j=0,2.*N_Elements(H)/N_cycles $
      DO cgPlotS, [Time[j], Time[j+1]], [M[j], M[j+1]], Color=elevColors[j+1]
    
  ENDIF 
  
  ;@GJ, 2023/9/8, plot the 2D PSF figure
  FOV = 2. * A / (gradient_field / 10.) ; in mm
  n_x = FLOOR(N_ELEMENTS(Time)/(5.*2.*f)) ;this is corrected number
;  print, 'FOV = ', STRING(FOV,format='(f6.2)')+' mm'
  ;@GJ, 2022/12/14, changing to offset
  n_y = n_x
  resolution = [FOV/n_x, FOV/n_y]
;  print, 'resolution: ', resolution
  y_axis = FINDGEN(n_y) * resolution[1]

  
  ;@GJ, 2023/10/11, design the phantom
  phantom_1D_array = DBLARR(n_x*4.) * 0.
  N_D_points = FLOOR(D_points*0.1/resolution[0])
  Npoints_ind = 0.
  FOR i=0, n_x*4.-1 DO BEGIN
    IF (i MOD N_D_points) EQ 0 AND Npoints_ind LT N_points THEN  BEGIN
      phantom_1D_array[i] = 1;i+1;1.
      Npoints_ind++
    ENDIF
  ENDFOR
  phantom_1D_array = SHIFT(phantom_1D_array, n_x/2.-(N_D_points*(N_points-1)/2.))
;  is_Npoints_even = (N_points MOD 2) EQ 0
;  if (is_Npoints_even) then begin
;    phantom_1D_array = SHIFT(phantom_1D_array, n_x/2.-(N_D_points*(N_points-1)/2.)-N_D_points/2.)
;  endif else begin
;    phantom_1D_array = SHIFT(phantom_1D_array, n_x/2.-(N_D_points*(N_points-1)/2.))
;  endelse
;  iplot, phantom_1D_array[0:n_x-1], title='1D phantom'
  
  H_y_array = (FINDGEN(n_y)-FLOOR(n_y/2.)) * resolution[1] * (gradient_field/10.)
  signal_array = DBLARR(N_ELEMENTS(Time), n_y) * 0.
  signal_x_array = DBLARR(N_ELEMENTS(Time), n_y) * 0.
  n_harmonic = 10.
  Harmonic_Amp_x_array = DBLARR(n_harmonic, n_y) * 0.
  Harmonic_Phase_x_array = DBLARR(n_harmonic, n_y) * 0.
  Harmonic_Real_x_array = DBLARR(n_harmonic, n_y) * 0.
  Harmonic_Imag_x_array = DBLARR(n_harmonic, n_y) * 0.
  signal_y_array = DBLARR(N_ELEMENTS(Time), n_y) * 0.
  ;@GJ, 2025/1/18, for calculating the x-y curve in 3rd harmonic complex plane  
  tau_ms_array = DBLARR(n_y) * 0.
  theta_0 = 0.
  tau_ms_0 = 0.
  M_x = DBLARR(N_ELEMENTS(Time), n_y) * 0.
  M_y = DBLARR(N_ELEMENTS(Time), n_y) * 0.
  FOR i=0, n_y-1 DO BEGIN
    ;@GJ, 2023/10/20, adding the DC field
    H_x = -A*cos(2.*!PI*f*Time) + A_DC*0.1
    x_velocity = A*2.*!PI*f*sin(2.*!PI*f*Time)
    H_y = H_y_array[i]
    H_total = SIGNUM(H_x) * SQRT(H_x^2 + H_y^2)

    ;calculate M
    M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
    index_nan = where(~FINITE(M_total), count)
    IF count GT 0 THEN BEGIN
      RESULT = INTERPOL(M_total, Time, Time[index_nan], /NAN)
      M_total[index_nan] = RESULT
    ENDIF

    ;do vector splitting
    ;@GJ, 2023/9/23, simple calculation
    M_x[*,i] = ABS(M_total) * H_x/ABS(H_total)
    M_y[*,i] = ABS(M_total) * H_y/ABS(H_total)

    n = N_elements(Time)
    signal_x = (M_x[1:*, i] - M_x[0:n-2, i]);/(Time[1]-Time[0])
    signal_y = (M_y[1:*, i] - M_y[0:n-2, i]);/(Time[1]-Time[0])
    ;@GJ, 2022/12/14, removing the NaN points
    index_nan_x = where(~FINITE(signal_x), count_x)
    IF count_x GT 0 THEN BEGIN
      RESULT_x = INTERPOL(signal_x, Time, Time[index_nan_x], /NAN)
      signal_x[index_nan_x] = RESULT_x
    ENDIF
    index_nan_y = where(~FINITE(signal_y), count_y)
    IF count_y GT 0 THEN BEGIN
      RESULT_y = INTERPOL(signal_y, Time, Time[index_nan_y], /NAN)
      signal_y[index_nan_y] = RESULT_y
    ENDIF
    
    ;@GJ, 2023/10/11, calculate the phantom 1D
    signal_x_phantom = signal_x * 0.
 ;   IF i EQ FLOOR(n_y/2.) THEN iplot, signal_x_phantom, title='signal'
    FOR m=0, n_x-1 DO BEGIN
      IF phantom_1D_array[m] NE 0 THEN BEGIN
        signal_x_phantom += SHIFT(signal_x, m-n_x/2.)*phantom_1D_array[m]
 ;       IF i EQ FLOOR(n_y/2.) THEN iplot, signal_x_phantom, /overplot
      ENDIF
    ENDFOR
    signal_x = signal_x_phantom/TOTAL(phantom_1D_array[0:n_x-1])
    
    tau_us = fields_brownian_time_calc(tau_B, beta, A, A_DC*0.1, H_y)
    tau_ms = tau_us / 1000.
    ;calculate the STED curve in 3rd harmonic complex plane
    tau_ms_array[i] = tau_ms
    IF H_y EQ 0 THEN BEGIN
      tau_ms_0 = tau_ms
      theta_0 = atan(3. * (2. * !PI * f) * tau_ms)
    ENDIF
    kernel = Debye_Kernel_sin(tau_ms, Time[0:100])
    signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
    signal_new_x = signal_new_x[n_elements(kernel)/2 : n_elements(signal_x) + n_elements(kernel)/2]
    signal_new_y = convol([dblarr(n_elements(kernel)/2)*0, signal_y, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
    signal_new_y = signal_new_y[n_elements(kernel)/2 : n_elements(signal_y) + n_elements(kernel)/2]
    t_ele = N_Elements(signal_new_x)
    ;delta_t = (Time[1]-Time[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
    X = (FINDGEN((t_ele - 1)/2) + 1)
    is_N_even = (t_ele MOD 2) EQ 0
    if (is_N_even) then $
      frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
    else $
      frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
    base_freq = MIN(ABS(frequency-1.0), base_ind)
    
    signal_FFT_x = FFT(signal_new_x)
    ;@GJ, 2023/10/20, get the harmonic amp and phase
    FOR m=0, n_harmonic-1 DO BEGIN
      Harmonic_Amp_x_array[m, i] = ABS(signal_FFT_x[base_ind*(m+1)])
      Harmonic_Phase_x_array[m, i] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*(m+1)]), Real_part(signal_FFT_x[base_ind*(m+1)]), /phase);IMAGINARY(signal_FFT_x[base_ind*(m+1)])
      Harmonic_Real_x_array[m, i] = REAL_PART(signal_FFT_x[base_ind*(m+1)])
      Harmonic_Imag_x_array[m, i] = IMAGINARY(signal_FFT_x[base_ind*(m+1)])
    ENDFOR
    
    signal_FFT_x[base_ind] *= 0.
    signal_new_no_base_x = REAL_PART(FFT(signal_FFT_x, /inverse))
    FOR j=0, N_ELEMENTS(Time)-1 DO BEGIN
      IF x_velocity[j] NE 0 THEN BEGIN
        signal_x_array[j,i] = signal_new_no_base_x[j]/(gradient_field/10.*x_velocity[j])
        signal_y_array[j,i] = signal_new_y[j]/(gradient_field/10.*x_velocity[j])
      ENDIF
    ENDFOR
  ENDFOR
  
  ;@GJ, 2025/1/18, calculate the STED curve in 3rd harmonic complex plane
  ;@GJ, 2025/1/18, for calculating the x-y curve in 3rd harmonic complex plane
;  nwtau = 3. * (2.*!PI * f) * tau_ms_array
;  theta_array = atan(nwtau/(1.+nwtau^2), 1./(1.+nwtau^2), /phase)
;  harmonic3_real_array = A / SQRT(A^2 + H_y_array^2) * f * cos(theta_array[125] - theta_array)
;  harmonic3_imag_array = A / SQRT(A^2 + H_y_array^2) * f * sin(theta_array[125] - theta_array)
;;  tan_factor = n_harmonic * f * (tau_ms_0 - tau_ms_array) / (1. + (n_harmonic * f)^2 * tau_ms_0 * tau_ms_array)
;  harmonic3_imag_equa = harmonic3_real_array * tan(theta_array[125] - theta_array)
;;  iplot, harmonic3_real_array, harmonic3_imag_array, thick=1, title='3rd Harmonic theoretical', /NO_SAVEPROMPT
;  iplot, harmonic3_real_array, harmonic3_imag_equa, /overplot, linestyle=2, thick=3, color='red'
;  iplot, H_y_array, harmonic3_real_array/MAX(harmonic3_real_array), title='3rd Harmonic theoretical', /NO_SAVEPROMPT
;  iplot, H_y_array, harmonic3_imag_array/MAX(harmonic3_imag_array), /overplot, color='blue'
    
  angle_0 = ATAN(Harmonic_Imag_x_array[2, n_y/2], Harmonic_Real_x_array[2, n_y/2], /phase) + 1./2.*!PI
;  angle_0 = 0.
  FOR i=0, n_y-1 DO BEGIN
    temp_real = -Harmonic_Real_x_array[*, i] * sin(angle_0) + Harmonic_Imag_x_array[*, i] * cos(angle_0)
    temp_imag = Harmonic_Imag_x_array[*, i] * sin(angle_0) + Harmonic_Real_x_array[*, i] * cos(angle_0)
    Harmonic_Real_x_array[*, i] = -temp_real
    Harmonic_Imag_x_array[*, i] = temp_imag
  ENDFOR
  
;  ;@GJ, 2025/1/18, check the subtraction for filtered BP
;  iplot, H_y_array, Harmonic_Real_x_array[2, *]/MAX(Harmonic_Real_x_array[2, *]), title='Harmonic 3', /NO_SAVEPROMPT
;  iplot, H_y_array, Harmonic_Imag_x_array[2, *]/MAX(Harmonic_Imag_x_array[2, *]), /overplot, color='blue'
  
  ;@GJ, 2025/1/31, calculate the equations of STED and G3ad
  Harmonic_Imag_x_array_new = REFORM(Harmonic_Imag_x_array[2, *]) * 0.
  n = 3.
  w = 2.*!PI * f
  FOR i=0, n_y-1 DO BEGIN
    ratio_IR = n * w * (tau_ms_array[n_y/2] -  tau_ms_array[i]) / (1. + (n * w)^2 * tau_ms_array[i] * tau_ms_array[n_y/2])
    Harmonic_Imag_x_array_new[i] = Harmonic_Real_x_array[2, i] * ratio_IR
  ENDFOR
  ;@GJ, 2025/1/31, calculate the G3ad form
  Harmonic3rd_Amp_x_array_new = SQRT(REFORM(Harmonic_Real_x_array[2, *])^2 + REFORM(Harmonic_Imag_x_array[2, *])^2)
  G3_ad_inv = f / Harmonic3rd_Amp_x_array_new / SQRT(1. + (n*w*tau_ms_array)^2)
  ; Define an n-element vector of Poisson measurement errors:
  measure_errors = SQRT(ABS(G3_ad_inv))
  result = LINFIT(H_y_array^2, G3_ad_inv, MEASURE_ERRORS=measure_errors)
  alpha_3 = SQRT(result[0]/result[1])
  beta_3 = SQRT(result[0]*result[1])
;  iplot, H_y_array^2, G3_ad_inv, xtitle='H_y^2', ytitle='1/G3_ad', title='G3_ad'
  print, 'alpha 3 :', alpha_3
  print, 'beta 3:', beta_3
;  iplot, Harmonic_Real_x_array[2, *], Harmonic_Imag_x_array[2, *], xrange=[0, MAX(Harmonic_Real_x_array[2, *])*1.1], yrange=[0, MAX(Harmonic_Real_x_array[2, *])*1.1], color='red', thick=1, xtitle='G3R [A.U.]', ytitle='G3I [A.U.]', title='STED with Fitting', /NO_SAVEPROMPT
;  iplot, Harmonic_Real_x_array[2, *], Harmonic_Imag_x_array_new, color='green', linestyle=1, thick=2, /overplot
    
  ;@GJ, 2025/1/12, draw the A3 complex
  IF wid GT 0 THEN BEGIN
    draw5 = widget_info(wid, find_by_uname='draw5')
    widget_control, draw5, get_value = drawWindow5
    wset, drawWindow5
    colors = cgScaleVector(Findgen(N_Elements(a3_array)), Min(IMAGINARY(a3_array)), Max(IMAGINARY(a3_array)))
    elevColors = Value_Locate(colors, IMAGINARY(a3_array))
    elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
    IF N_points EQ 1 THEN BEGIN
      max_a3 = MAX(IMAGINARY(a3_array), max_ind)
      temp_imag = INTERPOL(Harmonic_Imag_x_array[2, *],  Harmonic_Real_x_array[2, *], REAL_PART(a3_array[max_ind]))
      scale_factor = temp_imag / IMAGINARY(a3_array[max_ind])
      IF scale_factor LT 0. OR ~FINITE(scale_factor) THEN scale_factor = 1.
      plot, Harmonic_Real_x_array[2, *], Harmonic_Imag_x_array[2, *]/scale_factor, xrange=[0, 1.1*MAX(REAL_PART(a3_array))], linestyle=1, thick=2, background = 'ffffff'x, color = '000000'x, xtitle='R3', ytitle='I3'
      ;@GJ, 2025/1/31, add the fitted by the ratio of I to R
      oplot, Harmonic_Real_x_array[2, *], Harmonic_Imag_x_array_new/scale_factor, linestyle=0, thick=1, color = '111111'x
    ENDIF ELSE BEGIN
      plot, REAL_PART(a3_array), IMAGINARY(a3_array), /nodata, psym=4, background = 'ffffff'x, color = '000000'x, xtitle='R3', ytitle='I3'
    ENDELSE
    
    xyouts,20,drawSize.ysize*0.95,'3rd Harmonic', color = '000000'x, /DEVICE
    FOR j=0,N_ELEMENTS(a3_array)-2 do $
      cgPlotS, [REAL_PART(a3_array[j]),REAL_PART(a3_array[j])], [IMAGINARY(a3_array[j]),IMAGINARY(a3_array[j])], psym=4, Color=elevColors[j]
    cgPlotS, [REAL_PART(signal_FFT[3.*base_ind]),REAL_PART(signal_FFT[3.*base_ind])], [IMAGINARY(signal_FFT[3.*base_ind]),IMAGINARY(signal_FFT[3.*base_ind])], psym=2, Color=elevColors[N_ELEMENTS(a3_array)-1]
  ENDIF
  
  max_signal_y = MAX(signal_y_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8], maxind_y)
  signal_y_array(WHERE(signal_y_array LT -ABS(max_signal_y))) = 0.
  signal_y_array(WHERE(signal_y_array GT ABS(max_signal_y))) = 0.
  signal_y_array = signal_y_array + noise_level*max_signal_y*randomn(1, N_ELEMENTS(Time), n_y)
;  iimage, signal_y_array[0:2*n_x-1, 0:n_y-1], title='y image'
  
  max_signal_x = MAX(signal_x_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8], maxind_x)
  signal_x_array(WHERE(signal_x_array LT -ABS(max_signal_x))) = 0.
  signal_x_array(WHERE(signal_x_array GT ABS(max_signal_x))) = 0.
  signal_x_array = signal_x_array + noise_level*max_signal_x*randomn(1, N_ELEMENTS(Time), n_x)
;  iimage, signal_x_array[0:2*n_x-1, 0:n_y-1], title='x image'
  
  merged_image_x = SHIFT(signal_x_array[0:n_x-1, *], -shift_ele/half_cycle*n_x, 0) + SHIFT(REVERSE(signal_x_array[n_x:2*n_x-1, *], 1), shift_ele/half_cycle*n_x, 0)
  merged_image_y = SHIFT(signal_y_array[0:n_x-1, *], -shift_ele/half_cycle*n_x, 0) + SHIFT(REVERSE(signal_y_array[n_x:2*n_x-1, *], 1), shift_ele/half_cycle*n_x, 0)
;  iimage, merged_image, title='merged image'
  
  subtracted_image_x = SHIFT(signal_x_array[0:n_x-1, *], -shift_ele/half_cycle*n_x, 0) - SHIFT(REVERSE(signal_x_array[n_x:2*n_x-1, *], 1), shift_ele/half_cycle*n_x, 0)
  subtracted_image_y = SHIFT(signal_y_array[0:n_x-1, *], -shift_ele/half_cycle*n_x, 0) - SHIFT(REVERSE(signal_y_array[n_x:2*n_x-1, *], 1), shift_ele/half_cycle*n_x, 0)
  
  RESI_image_x = ABS(SHIFT(signal_x_array[0:n_x-1, *], -shift_ele/half_cycle*n_x, 0) - SHIFT(REVERSE(signal_x_array[n_x:2*n_x-1, *], 1), shift_ele/half_cycle*n_x, 0))
  RESI_image_y = ABS(SHIFT(signal_y_array[0:n_x-1, *], -shift_ele/half_cycle*n_x, 0) - SHIFT(REVERSE(signal_y_array[n_x:2*n_x-1, *], 1), shift_ele/half_cycle*n_x, 0))
  RESI_image = RESI_image_x/MAX(RESI_image_x) + RESI_image_y/MAX(RESI_image_y)
  
  ;@GJ, 2023/10/1, do the calculations without shift
  merged_image_x_no_shift = signal_x_array[0:n_x-1, *] + REVERSE(signal_x_array[n_x:2*n_x-1, *], 1)
  merged_image_y_no_shift = signal_y_array[0:n_x-1, *] + REVERSE(signal_y_array[n_x:2*n_x-1, *], 1)
  ;subtracted
  subtracted_image_x_no_shift = signal_x_array[0:n_x-1, *] - REVERSE(signal_x_array[n_x:2*n_x-1, *], 1)
  subtracted_image_y_no_shift = signal_y_array[0:n_x-1, *] - REVERSE(signal_y_array[n_x:2*n_x-1, *], 1)
  ;RESI image
  RESI_image_x_no_shift = ABS(signal_x_array[0:n_x-1, *] - REVERSE(signal_x_array[n_x:2*n_x-1, *], 1))
  RESI_image_y_no_shift = ABS(signal_y_array[0:n_x-1, *] - REVERSE(signal_y_array[n_x:2*n_x-1, *], 1))
  RESI_image_no_shift = RESI_image_x_no_shift/MAX(RESI_image_x_no_shift) + RESI_image_y_no_shift/MAX(RESI_image_y_no_shift)
  relax_image_x_no_shift = merged_image_x_no_shift-RESI_image_x_no_shift
  relax_image_x_no_shift[WHERE(relax_image_x_no_shift LT 0)] = 0.
  
  ;do the calculation
  ;@GJ, 2023/10/1, calculate the nearest max
  FOR k=5, n_x/2-2, 3 DO BEGIN
    IF RESI_image_x[n_x/2+k, n_y/2]-RESI_image_x[n_x/2+k-3, n_y/2] LT 0 THEN BEGIN
      first_max = RESI_image_x[half_cycle/2+k-3, n_y/2]
      BREAK
    ENDIF
  ENDFOR

  ;@GJ, 2023/10/1
  relax_image_x = relax_image_x_no_shift - RESI_image_x / ABS(first_max) * MAX(relax_image_x_no_shift)
  ;relax_image_x = merged_image_x-RESI_image_x
  relax_image_x[WHERE(relax_image_x LT 0)] = 0.
;  iimage, RESI_image, title='RESI image'
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 1
    
    draw3_r = widget_info(wid, find_by_uname='draw3_r')
    widget_control, draw3_r, get_value = drawWindow3_r
    wset, drawWindow3_r
    TVSCL, BYTSCL(CONGRID(ROTATE(merged_image_x,1), drawXSize, drawYSize), MAX=MAX(merged_image_x[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
    xyouts,drawSize.ysize*0.25,drawSize.ysize*0.9,'Standard', color = 'FFFFFF'x, /DEVICE
    ;@GJ, 2025/1/7, plot the multi-channel x-space
    ;iimage, CONGRID(ROTATE(merged_image_x,1), drawXSize, drawXSize)
    ;iimage, CONGRID(ROTATE(merged_image_x,1), drawXSize, drawXSize) + CONGRID(ROTATE(merged_image_x,0), drawXSize, drawXSize), title='Multi-Channel X-space'
    
    
    LOADCT, 3    
    draw3_rr = widget_info(wid, find_by_uname='draw3_rr')
    widget_control, draw3_rr, get_value = drawWindow3_rr
    wset, drawWindow3_rr
    TVSCL, BYTSCL(CONGRID(ROTATE(subtracted_image_x,1), drawXSize, drawYSize), MAX=MAX(subtracted_image_x[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
    xyouts,drawSize.ysize*0.15,drawSize.ysize*0.9,'Twin Subtracted', color = 'FFFFFF'x, /DEVICE
    
    draw6_r = widget_info(wid, find_by_uname='draw6_r')
    widget_control, draw6_r, get_value = drawWindow6_r
    wset, drawWindow6_r
    TVSCL, BYTSCL(CONGRID(ROTATE(RESI_image_x,1), drawXSize, drawYSize), MAX=MAX(RESI_image_x[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
    xyouts,drawSize.ysize*0.1,drawSize.ysize*0.9,'ABS (Twin Subtracted)', color = 'FFFFFF'x, /DEVICE
    
    LOADCT, 1
    draw6_rr = widget_info(wid, find_by_uname='draw6_rr')
    widget_control, draw6_rr, get_value = drawWindow6_rr
    wset, drawWindow6_rr
    TVSCL, BYTSCL(CONGRID(ROTATE(relax_image_x,1), drawXSize, drawYSize), MAX=MAX((merged_image_x-RESI_image_x)[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
    xyouts,drawSize.ysize*0.3,drawSize.ysize*0.9,'RESI', color = 'FFFFFF'x, /DEVICE
    
    ;@GJ, 2023/10/11, plot the RESI and phantom    
;    IF STRLEN(filename) LT 5 OR (FILE_INFO(filename)).exists EQ 0 THEN BEGIN
;      iimage, ROTATE(Harmonic_Amp_x_array,1), MAX=MAX(Harmonic_Amp_x_array[n_harmonic*0.1:n_harmonic*0.8, n_y*0.1:n_y*0.8]), title='Real Part'
;      iimage, ROTATE(Harmonic_Phase_x_array,1), MAX=MAX(Harmonic_Phase_x_array[n_harmonic*0.1:n_harmonic*0.8, n_y*0.1:n_y*0.8]), title='Imaginary Part'
;      iplot, y_axis, REFORM(Harmonic_Amp_x_array[2,*]), color='blue', xtitle='z axis [mm]', title='3rd Harmonic Amp'
;      iplot, y_axis, REFORM(Harmonic_Phase_x_array[2,*]), color='red', xtitle='z axis [mm]', title='3rd Harmonic Phase';/overplot
      H_y_2D_array = DBLARR(n_y, n_y)
      H_z_2D_array = DBLARR(n_y, n_y)
      FOR m=0, n_y-1 DO BEGIN
        H_y_2D_array[m, *] = H_y_array
        H_z_2D_array[*, m] = H_y_array
      ENDFOR
      Harmonic_Amp_yz_2D_array = DBLARR(n_y, n_y) * 0.
      Harmonic_Phase_yz_2D_array = DBLARR(n_y, n_y) * 0.
      Harmonic_Real_yz_2D_array = DBLARR(n_y, n_y) * 0.
      Harmonic_Imag_yz_2D_array = DBLARR(n_y, n_y) * 0.
      FOR m=0, n_y-1 DO BEGIN
        FOR n=0, n_y-1 DO BEGIN
          H_temp = SQRT(H_y_2D_array[m, n]^2 + H_z_2D_array[m, n]^2)
          IF ABS(H_temp-A_G*0.1) GT 0.02 THEN BEGIN
            Harmonic_Amp_yz_2D_array[m, n] = INTERPOL(Harmonic_Amp_x_array[2,*], H_y_array, H_temp)
            Harmonic_Phase_yz_2D_array[m, n] = INTERPOL(Harmonic_Phase_x_array[2,*], H_y_array, H_temp)
            Harmonic_Real_yz_2D_array[m, n] = INTERPOL(Harmonic_Real_x_array[2,*], H_y_array, H_temp)
            Harmonic_Imag_yz_2D_array[m, n] = INTERPOL(Harmonic_Imag_x_array[2,*], H_y_array, H_temp)
          ENDIF
        ENDFOR
      ENDFOR
;      iimage, Harmonic_Amp_yz_2D_array, title='Amp'
;      iimage, Harmonic_Phase_yz_2D_array, title='Phase'
      IF wid GT 0 THEN BEGIN
        DEVICE, DECOMPOSED = 0, RETAIN = 2
        LOADCT, 1
        draw3_0 = widget_info(wid, find_by_uname='draw3_0')
        widget_control, draw3_0, get_value = drawWindow3_0
        wset, drawWindow3_0
;        TVSCL, BYTSCL(CONGRID(Harmonic_Amp_yz_2D_array, drawXSize, drawYSize), MAX=MAX(Harmonic_Amp_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]), MIN=MIN(Harmonic_Amp_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
;        xyouts,drawSize.ysize*0.25,drawSize.ysize*0.9,'A3 Amp '+STRING(angle,format='(i4)'), color = 'FFFFFF'x, /DEVICE
        TVSCL, BYTSCL(CONGRID(Harmonic_Real_yz_2D_array, drawXSize, drawYSize), MAX=MAX(Harmonic_Real_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]), MIN=MIN(Harmonic_Real_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
        xyouts,drawSize.ysize*0.25,drawSize.ysize*0.9,'A3 Real '+STRING(angle,format='(i4)'), color = 'FFFFFF'x, /DEVICE

        LOADCT, 3
        draw6_0 = widget_info(wid, find_by_uname='draw6_0')
        widget_control, draw6_0, get_value = drawWindow6_0
        wset, drawWindow6_0
;        TVSCL, BYTSCL(CONGRID(Harmonic_Phase_yz_2D_array, drawXSize, drawYSize), MAX=MAX(Harmonic_Phase_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]), MIN=MIN(Harmonic_Phase_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
;        xyouts,drawSize.ysize*0.25,drawSize.ysize*0.9,'A3 Phase '+STRING(angle,format='(i4)'), color = 'FFFFFF'x, /DEVICE
        TVSCL, BYTSCL(CONGRID(Harmonic_Imag_yz_2D_array, drawXSize, drawYSize), MAX=MAX(Harmonic_Imag_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]), MIN=MIN(Harmonic_Imag_yz_2D_array[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]))
        xyouts,drawSize.ysize*0.25,drawSize.ysize*0.9,'A3 Imag '+STRING(angle,format='(i4)'), color = 'FFFFFF'x, /DEVICE

        DEVICE, DECOMPOSED = 1
      ENDIF
;      iplot, REFORM(Harmonic_Phase_x_array[*,n_y/2]), title='Phase'
;      iplot, phantom_1D_array[0:n_x-1], title='1D phantom'
;      relax_image_x_label = relax_image_x
;      relax_image_x_label[*, n_y/2] += phantom_1D_array*MAX(relax_image_x[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8])
;      iimage, ROTATE(relax_image_x_label,1), MAX=MAX(relax_image_x_label[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]), title='RESI and Phantom'
;      RESI_image_x_label = RESI_image_x
;      RESI_image_x_label[*, n_y/2] += phantom_1D_array*MAX(RESI_image_x[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8])
;      iimage, ROTATE(RESI_image_x_label,1), MAX=MAX(RESI_image_x_label[n_x*0.1:n_x*0.8, n_y*0.1:n_y*0.8]), title='Abs Sum and Phantom'
    ENDIF
    
    DEVICE, DECOMPOSED=1
;  ENDIF
  
  IF STRLEN(filename) LT 5 THEN return
  IF (FILE_INFO(filename)).exists EQ 0 THEN return
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
  ENDIF ELSE BEGIN
    image_sim_ori = DBLARR(512, 512)
    image_sim_ori[245:246, 200:280] = 1.
    image_sim_ori[258:259, 200:280] = 1.
    image_sim_ori[200:280, 245:246] = 1.
    image_sim_ori[265:270, 265:270] = 1.
    image_sim_ori = ROT(image_sim_ori, 119, 1, /INTERP)
    shift_dis = 100.;200.;100.;100.
    image_sim_ori = SHIFT(image_sim_ori, shift_dis, shift_dis)
;    image_sim_ori = SHIFT(image_sim_ori, shift_dis, 0)
    rows_S = N_ELEMENTS(image_sim_ori[0,*])
    cols_S = N_ELEMENTS(image_sim_ori[*,0])
    FOV_sim = fov_phantom; mm
    pixelSp_sim = [FOV_sim/DOUBLE(cols_S), FOV_sim/DOUBLE(rows_S)]
    print, 'pixelSp_sim: ', pixelSp_sim
  ENDELSE
  
  IF wid GT 0 THEN BEGIN
;    DEVICE, DECOMPOSED = 0, RETAIN = 2
;    LOADCT, 1
    draw3_r = widget_info(wid, find_by_uname='draw3_r')
    widget_control, draw3_r, get_value = drawWindow3_r
    wset, drawWindow3_r
    TVSCL, BYTSCL(CONGRID(image_sim_ori, drawXSize, drawYSize), MAX=MAX(image_sim_ori[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawSize.ysize*0.35,drawSize.ysize*0.9,'Phantom', color = 'FFFFFF'x, /DEVICE
;    DEVICE, DECOMPOSED = 1
  ENDIF
  
  ;@GJ, 2023/9/13, repeated scans
  scale_ratio = pixelSp_sim[0]/resolution[0]
  print, 'scale ratio: ', scale_ratio
  left_kernel_x_no_shift = CONGRID(REVERSE(signal_x_array[n_y:2*n_y-1, *], 1),n_x/scale_ratio,n_y/scale_ratio)
  ;@GJ, 2023/9/29, shift the kernel by elements
  left_kernel_x = SHIFT(left_kernel_x_no_shift, shift_ele/half_cycle*n_x/scale_ratio)
  right_kernel_x_no_shift = CONGRID(signal_x_array[0:n_y-1, *],n_x/scale_ratio,n_y/scale_ratio)
  right_kernel_x = SHIFT(right_kernel_x_no_shift, -shift_ele/half_cycle*n_x/scale_ratio)
  
  ;@GJ, 2023/9/23, make the size smaller than the sim image
  IF n_x/scale_ratio GT cols_S THEN BEGIN
    left_kernel_x_no_shift = left_kernel_x_no_shift[n_x/scale_ratio/2-cols_S/3:n_x/scale_ratio/2+cols_S/3, n_y/scale_ratio/2-rows_S/3:n_y/scale_ratio/2+rows_S/3]
    right_kernel_x_no_shift = right_kernel_x_no_shift[n_x/scale_ratio/2-cols_S/3:n_x/scale_ratio/2+cols_S/3, n_y/scale_ratio/2-rows_S/3:n_y/scale_ratio/2+rows_S/3]
    left_kernel_x = left_kernel_x[n_x/scale_ratio/2-cols_S/3:n_x/scale_ratio/2+cols_S/3, n_y/scale_ratio/2-rows_S/3:n_y/scale_ratio/2+rows_S/3]
    right_kernel_x = right_kernel_x[n_x/scale_ratio/2-cols_S/3:n_x/scale_ratio/2+cols_S/3, n_y/scale_ratio/2-rows_S/3:n_y/scale_ratio/2+rows_S/3]
  ENDIF
  resolution *= scale_ratio
  print, '(kernel) Cols*Rows:', n_x/scale_ratio, ' *', n_y/scale_ratio
  print, 'Resolution:', STRING(resolution[1],format='(f6.2)'), ' mm'
  print, 'FOV:', STRING(rows_S * resolution[1],format='(f6.2)'), ' mm'
  print, 'Cols*Rows:', STRING(cols_S,format='(f6.2)'), ' *', STRING(rows_S,format='(f6.2)')
  Dpp_H = ABS(H[max_ind] - H[min_ind]) ; mT
  ;@GJ, 2023/9/15, calculate number of pixels
  Dpp_N_pixels = Dpp_reso / resolution[1];
  print, 'Dpp (# pixels):'+STRING(Dpp_N_pixels,format='(f6.2)')
  
;  left_kernel_x = SHIFT(left_kernel_x_ori, j, 0)
;  right_kernel_x = SHIFT(right_kernel_x_ori, -j, 0)
    
; iimage, ABS(left_kernel_x-right_kernel_x)
    
;  n_angles = 4.;20.;4.;20.;4.;12.;4.
;  angle_array = [0,1,3,2] * 180. / n_angles
  phantom_merge = image_sim_ori * 0.
  phantom_subtr = image_sim_ori * 0.
  phantom_relax = DBLARR(cols_S, rows_S, n_angles) * 0.
  phantom_relax_no_shift = DBLARR(cols_S, rows_S, n_angles) * 0.
  phantom_relax_rec = image_sim_ori * 0.
  FOR i=0, n_angles-1 DO BEGIN
    ;do rotation
    IF n_angles GT 1 THEN angle = (n_angles-1 - i) * 180. / n_angles + 90.
    print, 'angle = ', angle, ' degree'
    image_sim = ROT(image_sim_ori, angle, 1.0, /INTERP)

    ;do convolution
    ;resolution become 5 mm
    phantom_left_x = CONVOL_FFT(image_sim, left_kernel_x)
    ;@GJ, 2023/9/15, adding noise to the image
    phantom_left_x = phantom_left_x + noise_level*MEAN(phantom_left_x)*randomn(1, cols_S, rows_S);NOISE_PICK(phantom_left_x, 0.1, ITER=10);NOISE_SCATTER(phantom_left_x, LEVELS=[0.1,0.1])

    phantom_right_x = CONVOL_FFT(image_sim, right_kernel_x)
    ;@GJ, 2023/9/15, adding noise to the image
    phantom_right_x = phantom_right_x + noise_level*MEAN(phantom_right_x)*randomn(1, cols_S, rows_S);NOISE_PICK(phantom_right_x, 0.1, ITER=10);NOISE_SCATTER(phantom_right_x, LEVELS=[0.1,0.1])

    phantom_merge_x = phantom_left_x + phantom_right_x
    phantom_merge_x = ROT(phantom_merge_x, -angle, 1.0, /INTERP)
    phantom_subtr_x = ABS(phantom_left_x - phantom_right_x)
    phantom_subtr_x = ROT(phantom_subtr_x, -angle, 1.0, /INTERP)
    phantom_relax_x = phantom_merge_x - phantom_subtr_x
    phantom_relax_x[WHERE(phantom_relax_x LT 0)] = 0.
   
    ;@GJ, 2023/10/1, no shift
    ;;do convolution
    ;resolution become 5 mm
    phantom_left_x_no_shift = CONVOL_FFT(image_sim, left_kernel_x_no_shift)
    ;@GJ, 2023/9/15, adding noise to the image
    phantom_left_x_no_shift = phantom_left_x_no_shift + noise_level*MEAN(phantom_left_x_no_shift)*randomn(1, cols_S, rows_S);NOISE_PICK(phantom_left_x, 0.1, ITER=10);NOISE_SCATTER(phantom_left_x, LEVELS=[0.1,0.1])

    phantom_right_x_no_shift = CONVOL_FFT(image_sim, right_kernel_x_no_shift)
    ;@GJ, 2023/9/15, adding noise to the image
    phantom_right_x_no_shift = phantom_right_x_no_shift + noise_level*MEAN(phantom_right_x_no_shift)*randomn(1, cols_S, rows_S);NOISE_PICK(phantom_right_x, 0.1, ITER=10);NOISE_SCATTER(phantom_right_x, LEVELS=[0.1,0.1])

    phantom_merge_x_no_shift = phantom_left_x_no_shift + phantom_right_x_no_shift
    phantom_merge_x_no_shift = ROT(phantom_merge_x_no_shift, -angle, 1.0, /INTERP)
    phantom_subtr_x_no_shift = ABS(phantom_left_x_no_shift - phantom_right_x_no_shift)
    phantom_subtr_x_no_shift = ROT(phantom_subtr_x_no_shift, -angle, 1.0, /INTERP)
    phantom_relax_x_no_shift = phantom_merge_x_no_shift - phantom_subtr_x_no_shift
    phantom_relax_x_no_shift[WHERE(phantom_relax_x_no_shift LT 0)] = 0. 
    phantom_relax_no_shift[*,*,i] = phantom_relax_x_no_shift
     
    phantom_relax[*,*,i] = phantom_relax_x
    ;@GJ, 2023/9/15, by averaging the values giving the best result
    phantom_relax_rec += REFORM(phantom_relax[*,*,i])

    IF wid GT 0 THEN BEGIN
      DEVICE, DECOMPOSED = 0, RETAIN = 2
      LOADCT, 1
      draw3_rr = widget_info(wid, find_by_uname='draw3_rr')
      widget_control, draw3_rr, get_value = drawWindow3_rr
      wset, drawWindow3_rr
      phantom_merge_x[WHERE(phantom_merge_x LT 0)] = 0.
      TVSCL, BYTSCL(CONGRID(phantom_merge_x, drawXSize, drawYSize), MAX=MAX(phantom_merge_x[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
      xyouts,drawSize.ysize*0.25,drawSize.ysize*0.9,'Standard '+STRING(angle,format='(i4)'), color = 'FFFFFF'x, /DEVICE

      LOADCT, 3
      draw6_r = widget_info(wid, find_by_uname='draw6_r')
      widget_control, draw6_r, get_value = drawWindow6_r
      wset, drawWindow6_r
      TVSCL, BYTSCL(CONGRID(phantom_subtr_x, drawXSize, drawYSize), MAX=MAX(phantom_subtr_x[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
      xyouts,drawSize.ysize*0.1,drawSize.ysize*0.9,'ABS (Twin Subtracted) '+STRING(angle,format='(i4)'), color = 'FFFFFF'x, /DEVICE
      
;      LOADCT, 1
      draw6_rr = widget_info(wid, find_by_uname='draw6_rr')
      widget_control, draw6_rr, get_value = drawWindow6_rr
      wset, drawWindow6_rr
      TVSCL, BYTSCL(CONGRID(REFORM(phantom_relax_rec), drawXSize, drawYSize), MAX=MAX(phantom_relax_rec[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
      xyouts,drawSize.ysize*0.3,drawSize.ysize*0.9,'RESI', color = 'FFFFFF'x, /DEVICE
      DEVICE, DECOMPOSED = 1
    ENDIF
    wait, 0.5
  ENDFOR


  
  
;  iimage, phantom_relax_rec, title='reconstrcted image'

  ; Clean up.
  !P.Multi = 0  
end

pro slider_A_DC_event, ev
  widget_control, ev.id, get_value = A_DC
  widget_control, ev.top, get_uvalue = parameters
  parameters.A_DC = A_DC;
  IF A_DC EQ 0 THEN parameters.a3_array = PTR_NEW(COMPLEX(0., 0., /double));
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
end

pro slider_A_G_event, ev
  widget_control, ev.id, get_value = A_G
  widget_control, ev.top, get_uvalue = parameters
  parameters.A_G = A_G;
  IF A_G EQ 0 THEN parameters.a3_array = PTR_NEW(COMPLEX(0., 0., /double));
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
end

pro slider_A_sin_event, ev
  widget_control, ev.id, get_value = A
  widget_control, ev.top, get_uvalue = parameters
  parameters.A = A;
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
end

pro slider_f_sin_event, ev
  widget_control, ev.id, get_value = f
  widget_control, ev.top, get_uvalue = parameters
  parameters.f = f;
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
end

pro slider_particle_size_sin_event, ev
  widget_control, ev.id, get_value = particle_size
  widget_control, ev.top, get_uvalue = parameters
  parameters.particle_size = particle_size;
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
end

pro slider_tao_sin_event, ev
  widget_control, ev.id, get_value = tau
  widget_control, ev.top, get_uvalue = parameters
  parameters.tau = tau;
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end
pro slider_gradient_field_sin_event, ev
  widget_control, ev.id, get_value = gradient_field
  widget_control, ev.top, get_uvalue = parameters
  parameters.gradient_field = gradient_field;
  widget_control, ev.top, set_uvalue = parameters
;  print, 'gradient field: ', parameters.gradient_field, '*0.1 mT/mm'
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
end
pro PhantomSelection_event, ev
  widget_control, ev.top, get_uvalue = parameters
  filename = dialog_pickfile(TITLE='Select Phantom Image File', FILTER=['*.jpg', '*.png'], /MUST_EXIST, PATH=parameters.directory)
  ;read phantom image
;  fn_dir = 'F:\MPI_Tianjie\Nanorubble_MPI\phantomSimulation\'
;  ;fn_name = 'S_shape.png'
;  ;  fn_name = 'letter32.png'
;  ;  fn_name = 'SUA_letter24.png'
;  ;  fn_name = 'Mammo2.png'
;  fn_name = 'Mammo4.png'
;  ;  fn_name = 'Mammo5.jpg'
;  fn_name = 'MRA000.png'
;;  fn_name = 'S_shape.png'
;  fn_dir = 'C:\D_drive\MPI_Tianjie\RotationalDriftMPI\'
;  fn_name = 'MRA_ori.png'
;  filename = fn_dir + fn_name

  IF STRLEN(filename) GT 1 THEN BEGIN
    parameters.filename = filename
    parameters.directory = FILE_DIRNAME(filename)
    sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
  ENDIF ELSE BEGIN
    parameters.filename = filename
    sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
  ENDELSE
  widget_control, ev.top, set_uvalue = parameters
end

pro Droplist_Phantom_Name_event, ev
  widget_control, ev.top, get_uvalue = parameters
  listValue = WIDGET_INFO( ev.id, /DROPLIST_SELECT)
  parameters.phantom_name = listValue
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  widget_control, ev.top, set_uvalue = parameters
end

pro Droplist_Filter_Type_event, ev
  widget_control, ev.top, get_uvalue = parameters
  listValue = WIDGET_INFO( ev.id, /DROPLIST_SELECT)
  parameters.filter_type = (parameters.FilterTypeList)[listValue]
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  widget_control, ev.top, set_uvalue = parameters
end

pro slider_Npoints_event, ev
  widget_control, ev.id, get_value = N_points
  widget_control, ev.top, get_uvalue = parameters
  parameters.N_points = N_points;
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro slider_Dpoints_event, ev
  widget_control, ev.id, get_value = D_points
  widget_control, ev.top, get_uvalue = parameters
  parameters.D_points = D_points;
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro slider_fov_event, ev
  widget_control, ev.id, get_value = fov_phantom
  widget_control, ev.top, get_uvalue = parameters
  parameters.fov_phantom = fov_phantom;
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro slider_n_angles_event, ev
  widget_control, ev.id, get_value = n_angles
  widget_control, ev.top, get_uvalue = parameters
  parameters.n_angles = n_angles;
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro button_n_angles_event, ev
  ;widget_control, ev.id, get_value = n_angles
  widget_control, ev.top, get_uvalue = parameters
  Msat = 0.551/(4.0 * !PI * 1.e-7); / (4.e-7 * !PI)
;  k_B = 1.380469 ; 10^-23 J/K
;  u0 = 4.0 * !PI * (0.1^4) ; T*m/kA
  FWHM = 34. * 1.38e-23 * 300. /(parameters.gradient_field/10. * !PI * Msat * (parameters.particle_size * 1.e-9)^3)
  Kmax = 1. / 2. / (FWHM * 1000.)
  NTHETA = !PI * Kmax * parameters.fov_phantom ;1. * CEIL(!PI * NRHO/2.);127.;181.;27.;181
  IF NTHETA GT 180 THEN NTHETA = 180. 
  slider_N_angles = widget_info(ev.top, find_by_uname='N_angles')
  widget_control, slider_N_angles, set_value = CEIL(NTHETA)
  parameters.N_angles = CEIL(NTHETA)
  widget_control, ev.top, set_uvalue = parameters
  radon_G3STED, ev.top, parameters.phantom_name, parameters.filter_type, parameters
  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro slider_angle_event, ev
  widget_control, ev.id, get_value = angle
  widget_control, ev.top, get_uvalue = parameters
  parameters.angle = angle;
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro slider_shift_event, ev
  widget_control, ev.id, get_value = shift_perc
  widget_control, ev.top, get_uvalue = parameters
  parameters.shift_perc = shift_perc;
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation,ev.top, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc

end

pro quit_event, ev
  widget_control, ev.top, set_uvalue = parameters
  ;  Destroy widget hierarchy.
  ;
  WIDGET_CONTROL, ev.top, /DESTROY

  RETURN
end
;@GJ, 2023/8/13, Zhang Yifei helped developing this part
;@GJ, 2023/12/10, Fully modified to include the offset DC and Gradient fields for focal modulation MPI
pro Sin_Excitation

  ;��ʼֵ
  A = 4.
  A_DC = 0.;(*0.1mT)
  A_G = 0.;[*0.1mT]
  particle_size = 30.
  f = 2.;5.;25.;
  tau = 20  ;(ms)
  flat_portion = 50 ;
  gradient_field = 11.3;20. ; (*0.1mT/mm)
  signal_ratio = 0.;
  RESI_ratio = 1.;
  STD_reso = 1.;
  RESI_reso = 1.;
  Dpp_reso = 1.;
  noise_level = 0.01;
  N_points = 132.;1.;2.;
  D_points = 20.; *0.1mm
  fov_phantom = 190.;250.;
  n_angles = 27;4.;
  angle = 90.;
  shift_perc = 43.;100.;
  a3_array = PTR_NEW([COMPLEX(0, 0, /double)])
  device, decomposed=1
  device, get_screen_size = screenSize
  baseXSize = screenSize[0] * 0.9
  baseYSize = screenSize[1] * 0.9

  base = widget_base(title='Sin Excitation & Focal Modulation of MNPs', xsize = baseXSize, ysize = baseYSize, /column, MBAR=barBase)
  
  ;  Create the quit button
  ;
  wFileButton = WIDGET_BUTTON(barBase, VALUE= 'Phantom', /MENU)

  wFileSelectButton = WIDGET_BUTTON(wFileButton, $
    VALUE='Select a phantom', event_pro = 'PhantomSelection_event')
  wQuitButton = WIDGET_BUTTON(wFileButton, $
    VALUE='Quit', event_pro = 'quit_event')

  menubase = widget_base(base, xsize = baseXSize/3, /row)

  sliderbase = widget_base(menubase, /column)
  ;base_A
  base_A = widget_base(sliderbase, /row)
  label_A = widget_label(base_A, value = "Excitation wave amplitude(mT) :")
  slider_A_sin = widget_slider(base_A, xsize=250, MINIMUM=1, MAXIMUM=50, event_pro = 'slider_A_sin_event')
  ;base_A_DC
  base_A_DC = widget_base(sliderbase, /row)
  label_A_DC = widget_label(base_A_DC, value = "DC offset amplitude(*0.1mT) :")
  slider_A_DC = widget_slider(base_A_DC, xsize=250, MINIMUM=0, MAXIMUM=500, event_pro = 'slider_A_DC_event')
  ;base_A_G, 2023/12/10, adding the slide bar about gradient field
  base_A_G = widget_base(sliderbase, /row)
  label_A_G = widget_label(base_A_G, value = "Orth Gradient amplitude(*0.1mT) :")
  slider_A_G = widget_slider(base_A_G, xsize=250, MINIMUM=0, MAXIMUM=500, event_pro = 'slider_A_G_event')
  ;base_f
  base_f = widget_base(sliderbase, /row)
  label_f = widget_label(base_f, value =  "Excitation wave frequency(kHz):" )
  slider_f_sin = widget_slider(base_f, xsize=250, MINIMUM=1, MAXIMUM=10, event_pro = 'slider_f_sin_event')
  ;base_particle_size
  base_particle_size = widget_base(sliderbase, /row)
  label_particle_size = widget_label(base_particle_size, value = "Particle size(nm)             :")
  slider_particle_size_sin = widget_slider(base_particle_size, xsize=250, MINIMUM=10, MAXIMUM=100, event_pro = 'slider_particle_size_sin_event')
  ;base_tao
  base_tao = widget_base(sliderbase, /row)
  label_tao = widget_label(base_tao, value = "relaxation time(us)           :")
  slider_tao_sin = widget_slider(base_tao, xsize=250, MINIMUM=1, MAXIMUM=100, event_pro = 'slider_tao_sin_event', uname='slider_tao_sin')
  ;base_gradient field
  base_gradient_field = widget_base(sliderbase, /row)
  label_gradient_field = widget_label(base_gradient_field, value = "gradient field(*0.1mT/mm)           :")
  slider_gradient_field_sin = widget_slider(base_gradient_field, xsize=250, MINIMUM=1, MAXIMUM=100, event_pro = 'slider_gradient_field_sin_event')


  picbase = widget_base(menubase, /column)
  MPIjpg_path = '.\lib\icon_pics\cosExcitation.jpg'
  READ_JPEG, MPIjpg_path, MPIjpg, true=1
  picwindow = WIDGET_draw(picbase, uvalue='pic_win', xsize = (size(MPIjpg))[2], ysize = (size(MPIjpg))[3])
  
  phantombase = widget_base(menubase, /column)
  base_Phantom = widget_base(phantombase, /row)
;  base_PhantomSel = widget_base(base_Phantom, /row)
;  wFileSelectButton = WIDGET_BUTTON(base_PhantomSel, value='Select a phantom image file', event_pro='PhantomSelection_event')
  ;phantom name
  base_PhantomName = widget_base(base_Phantom, /row)
  wPhantomName = WIDGET_LABEL(base_PhantomName, VALUE='Phantom name')
  Phantom_Name = 'Xidian University'
  PhantomNameList = [Phantom_Name, 'Shepp-Logan', 'Two Plates', 'Hollow Plates', 'Vessel', 'Vessel Thermia', 'G', 'Two Bars', 'Two Close Bars']
  wPhantomDroplist = WIDGET_DROPLIST(base_PhantomName, $
    VALUE=PhantomNameList, $
    UVALUE='PhantomNameList', event_pro = 'droplist_phantom_name_event')
  base_FilterType = widget_base(base_Phantom, /row)
  wFilterType = WIDGET_LABEL(base_FilterType, VALUE='FBP kernel')
  Filter_type = 'Ramp'
  FilterTypeList = [Filter_type, 'Shepp-Logan', 'Hamming', 'Cosine']
  wFilterTypeDroplist = WIDGET_DROPLIST(base_FilterType, $
    VALUE=FilterTypeList, $
    UVALUE='FilterTypeList', event_pro = 'droplist_filter_type_event')

  ;N points selection
  base_Npoints = widget_base(phantombase, /row)
  label_Npoints = widget_label(base_Npoints, value = "N points:")
  slider_Npoints = widget_slider(base_Npoints, xsize=250, MINIMUM=1, MAXIMUM=500, uname='N_points', event_pro = 'slider_Npoints_event')
  ;Distance between points selection
  base_Dpoints = widget_base(phantombase, /row)
  label_Dpoints = widget_label(base_Dpoints, value = "Distance between points (*0.1mm):")
  slider_Dpoints = widget_slider(base_Dpoints, xsize=250, MINIMUM=1, MAXIMUM=100, uname='D_points', event_pro = 'slider_Dpoints_event')
  ;FOV selection
  base_FOV = widget_base(phantombase, /row)
  label_FOV = widget_label(base_FOV, value = "FOV (mm):")
  slider_FOV = widget_slider(base_FOV, xsize=250, MINIMUM=1, MAXIMUM=500, uname='FOV_phantom', event_pro = 'slider_fov_event')
  ;n_angles selection
  base_n_angles = widget_base(phantombase, /row)
  label_n_angles = widget_label(base_n_angles, value = "# of angles:")
  slider_n_angles = widget_slider(base_n_angles, xsize=250, MINIMUM=1, MAXIMUM=200, uname='N_angles', event_pro = 'slider_n_angles_event')
  button_n_angles = widget_button(base_n_angles, value='Auto', event_pro='button_n_angles_event')
  ;angle selection
  base_angle = widget_base(phantombase, /row)
  label_angle = widget_label(base_angle, value = "Angle (degree):")
  slider_angle = widget_slider(base_angle, xsize=250, MINIMUM=0, MAXIMUM=360, uname='Angle_deg', event_pro = 'slider_angle_event')
  ;shift selection
  base_shift = widget_base(phantombase, /row)
  label_shift = widget_label(base_shift, value = "Shift (%):")
  slider_shift = widget_slider(base_shift, xsize=250, MINIMUM=-200, MAXIMUM=200, uname='Shift_perc', event_pro = 'slider_shift_event')
  
  ;��ͼ��
  drawXSize = baseXSize / 6
  drawYSize = baseYSize / 3.1
  drawbase = widget_base(base, /column)
  drawbase1 = widget_base(drawbase, /row)
  draw1 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw1')
  draw2 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw2')
  draw3 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw3')
  draw3_0 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw3_0')
  draw3_r = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw3_r')
  draw3_rr = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw3_rr')
  drawbase2 = widget_base(drawbase, /row)
  draw4 = widget_base(drawbase2, /column)
  draw4a = widget_draw(draw4, xsize = drawXSize, ysize = drawYSize/2.1, retain=2, uname='draw4a')
  draw4b = widget_draw(draw4, xsize = drawXSize, ysize = drawYSize/2.1, retain=2, uname='draw4b')
  draw5 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw5')
  draw6 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw6')
  draw6_0 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw6_0')
  draw6_r = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw6_r')
  draw6_rr = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw6_rr')

  widget_control, base, /realize

  widget_control, picwindow, get_value = drawPic
  wset, drawPic
  TVSCL, MPIjpg, true=1
  
  parameters = {A:A, A_DC:A_DC, A_G:A_G, portion:flat_portion, particle_size:particle_size, f:f, tau:tau, gradient_field:gradient_field, a3_array:a3_array, signal_ratio:signal_ratio, RESI_ratio:RESI_ratio, STD_reso:STD_reso, RESI_reso:RESI_reso, Dpp_reso:Dpp_reso, drawXSize:drawXSize, drawYSize:drawYSize, directory:'.\lib\icon_pics\phantom\', noise_level:noise_level, filename:'', phantom_name: 0, filter_type:filter_type, FilterTypeList:FilterTypeList, fov_phantom:fov_phantom, n_angles:n_angles, angle:angle, shift_perc:shift_perc, N_points:N_points, D_points:D_points}
  widget_control, base, set_uvalue = parameters

  ;����slider�����value��ʼֵ
  widget_control, slider_A_sin, set_value = A
  widget_control, slider_A_DC, set_value = A_DC
  widget_control, slider_A_G, set_value = A_G
  widget_control, slider_f_sin, set_value = f
  widget_control, slider_particle_size_sin, set_value = particle_size
  widget_control, slider_tao_sin, set_value = tau
  widget_control, slider_gradient_field_sin, set_value = gradient_field
  widget_control, slider_Npoints, set_value = N_points
  widget_control, slider_Dpoints, set_value = D_points
  widget_control, slider_FOV, set_value = fov_phantom
  widget_control, slider_n_angles, set_value = n_angles
  widget_control, slider_angle, set_value = angle
  widget_control, slider_shift, set_value = shift_perc
  sin_simulation, base, parameters.A, parameters.A_DC, parameters.A_G, parameters.particle_size, parameters.f, parameters.tau, parameters.gradient_field, parameters.signal_ratio, parameters.RESI_ratio, parameters.RESI_reso, parameters.STD_reso, parameters.Dpp_reso, parameters.drawXSize, parameters.drawYSize, parameters.noise_level, parameters.filename, parameters.N_points, parameters.D_points, parameters.fov_phantom, parameters.n_angles, parameters.angle, parameters.shift_perc
  radon_G3STED, base, parameters.phantom_name, parameters.filter_type, parameters
  xmanager, 'Sin_Excitation', base
end