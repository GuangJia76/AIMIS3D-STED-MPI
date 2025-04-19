;@GJ, 2024/6/8, keyhole method
PRO keyhole_method, wid, parameters, index_A_DC, n_x, n_z, gradient_ratio, image_sim_ori

  ;@GJ, 2024/6/4, find the signal rato for noise DB
  signal_ratio = DBLARR(6)
  FOR k=0, 5 DO signal_ratio[k] = MEAN((*parameters.Sx_array)[k * 2. + 2, 0, *, *])
  signal_ratio /= signal_ratio[0]
  noise_level_percentage = 1. / (10.^(parameters.noise_level/10.)*signal_ratio) * 100.

  ;@GJ, 2024/6/7, calculate the image with noise
  focal_coor = FINDGEN(n_z) - n_z/2.
  focal_image_0 = DBLARR(n_x, n_z)
  total_comb_image = DBLARR(4*n_x, 7*n_z)
  A_G_Array = *parameters.A_G_Array
  zero_ind=WHERE(ABS(A_G_Array) EQ 0)
  temp_A_G_array = A_G_array[zero_ind:*]
  h_FWHM_0 = DBLARR(6) * 0.
  h_MTF_0 = DBLARR(6) * 0.
  MTF_peak_1st = DBLARR(6) * 0.
  noisy_image_FFT_div_0 = DCOMPLEXARR(6, n_x, n_z)
  line_profile_FFT = DBLARR(7, n_x/2.+1) * 0.
  line_profile_nCount = DBLARR(7, n_x/2.+1) * 0.
  FOR k=0, 5 DO BEGIN
    harmonics_n = k * 2. + 3
    curr_Sx_array = REFORM((*parameters.Sx_array)[harmonics_n-1, parameters.harmonics_comp, *, *])

    ;calculate the middle without offset
    focal_line = REFORM(curr_Sx_array[*, index_A_DC])
    FOR i=0, n_x-1 DO BEGIN
      FOR j=0, n_z-1 DO BEGIN
        focal_image_0[i, j] = INTERPOL([focal_line, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
      ENDFOR
    ENDFOR

    ;@GJ, 2024/6/7, calculate the central frequency FWHM
    donut_Harmonic_0 = REFORM(focal_line[zero_ind:*])
    maxD_Harmonic_0 = MAX(donut_Harmonic_0, maxD_ind_Harmonic_0)
    FWHM_ind_0 = INTERPOL(temp_A_G_array, donut_Harmonic_0, 0.5*maxD_Harmonic_0)
    IF finite(FWHM_ind_0) THEN h_FWHM_0[k] = 2. * FWHM_ind_0
    h_MTF_0[k] = 1./h_FWHM_0[k]
    A_temp = image_sim_ori * 0.
    B_temp = image_sim_ori * 0.
    MTF_peak_1st[k] = parameters.FOV / (h_FWHM_0[k]/(parameters.gradient_field_yz_current * 0.1)) / n_x
    A_temp[n_x/2.*(1.+MTF_peak_1st[k])+1, n_z/2.+1] = 255.
    ;    A_temp[drawXSize/2.*(1.-MTF_peak)+1, drawYSize/2.+1] = 255.
    FOR i_r = 0, 360 DO B_temp += ROT(A_temp, i_r, 1., /INTERP)
    MTF_indices_1st = WHERE(B_temp GT 128)

    image_sim_conv = CONVOL_FFT(image_sim_ori, focal_image_0)
    noisy_image_sim_conv = NOISE_SCATTER(BYTSCL(image_sim_conv), LEVELS=[noise_level_percentage[k]/100.,noise_level_percentage[k]/100.])/255.*MAX(image_sim_conv); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)

    focal_image_0_FFT = ABS(FFT(focal_image_0,/center,/double))
    noisy_image_sim_conv_FFT = FFT(noisy_image_sim_conv,/center,/double)
;    noisy_image_FFT_div_0[k, *, *] = noisy_image_sim_conv_FFT/focal_image_0_FFT
    noisy_image_FFT_div = noisy_image_sim_conv_FFT/focal_image_0_FFT
    noisy_image_FFT_div_0[k, *, *] = noisy_image_FFT_div
    
    ;@GJ, 2024/6/13, calculate the line profile
    FOR i=0, n_x-1 DO BEGIN
      FOR j=0, n_z-1 DO BEGIN
        dist_pixel = FLOOR(SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2))
        IF dist_pixel LE N_ELEMENTS(line_profile_FFT[k, *])-1 THEN BEGIN
          line_profile_nCount[k, dist_pixel]++
          line_profile_FFT[k, dist_pixel] += ALOG(ABS(noisy_image_FFT_div_0[k, i, j]))
        ENDIF
      ENDFOR
    ENDFOR
    line_profile_FFT[k, *] = line_profile_FFT[k, *]/FLOAT(line_profile_nCount[k, *])
    IF k EQ 0 THEN iplot, line_profile_FFT[k, *], /NO_SAVEPROMPT ELSE iplot, line_profile_FFT[k, *], /overplot
    
    noisy_image_FFT_div[where(focal_image_0_FFT LT 0.001*MAX(focal_image_0_FFT))] = COMPLEX(0, 0)
    ;  iimage, BYTSCL(ALOG(noisy_image_FFT_div)), title='1st Noisy FFT Div'
    rec_noisy_image = ABS(FFT(noisy_image_FFT_div, /INVERSE,/center,/double))

    ;@GJ, 2024/6/7, display the central frequency
    temp_display =  REFORM(BYTSCL(ALOG(noisy_image_FFT_div_0[k, *, *]),/NAN))
    ;    temp_display[MTF_indices_1st] = 255
    combined_image =  [temp_display, BYTSCL(rec_noisy_image,/NAN)]
    ;    iimage, combined_image, title=STRING(harmonics_n, format='(I0)')+' harmonics'
    total_comb_image[0:2*n_x-1, (7-k-1)*n_z:(7-k)*n_z-1] = combined_image
  ENDFOR

  ;@GJ< 2024/6/7, reconbine the k-space using keyhole method
  noisy_image_FFT_div_com = REFORM(noisy_image_FFT_div_0[0., *, *]) * 0.
  min_image_FFT_div_com = REFORM(noisy_image_FFT_div_0[0., *, *]) * 0.
  mask_image = DBLARR(n_x, n_z) * 0. + 100.
  min_mask_image = DBLARR(n_x, n_z) * 0.
  FOR i=0, n_x-1 DO BEGIN
    FOR j=0, n_z-1 DO BEGIN
      ;@GJ, calculate the minimal part
      FFT_abs_array_temp = ABS(REFORM(noisy_image_FFT_div_0[*, i, j]))
      min_array = MIN(FFT_abs_array_temp, min_array_ind)
      min_mask_image[i,j] = min_array_ind + 1.
      min_image_FFT_div_com[i,j] = noisy_image_FFT_div_0[min_array_ind, i, j]
      
      dist_pixel = SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)
      harmonic_ind = INTERPOL(FINDGEN(6), MTF_peak_1st*n_x/2.*1.5, dist_pixel)
      IF harmonic_ind LE 0 THEN BEGIN
        noisy_image_FFT_div_com[i,j] = noisy_image_FFT_div_0[0, i, j]
        mask_image[i,j] = 0.
      ENDIF ELSE BEGIN
        IF harmonic_ind LT 5 THEN BEGIN
          noisy_image_FFT_div_com[i,j] = noisy_image_FFT_div_0[CEIL(harmonic_ind), i, j]
          mask_image[i,j] = CEIL(harmonic_ind)
        ENDIF
      ENDELSE
    ENDFOR
  ENDFOR

  ;@GJ, 2024/6/13, calculate the line profile
  FOR i=0, n_x-1 DO BEGIN
    FOR j=0, n_z-1 DO BEGIN
      dist_pixel = FLOOR(SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2))
      IF dist_pixel LE N_ELEMENTS(line_profile_FFT[6, *])-1 THEN BEGIN
        line_profile_nCount[6, dist_pixel]++
        line_profile_FFT[6, dist_pixel] += ALOG(ABS(min_image_FFT_div_com[i, j]))
      ENDIF
    ENDFOR
  ENDFOR
  line_profile_FFT[6, *] = line_profile_FFT[6, *]/FLOAT(line_profile_nCount[6, *])
  iplot, line_profile_FFT[6, *], color='red', /overplot
  ;@GJ, 2024/6/14, smooth the profiles with searching for the minimum
  min_value = MIN(SMOOTH(REFORM(line_profile_FFT[6, *]), 10), min_ind)
  print, 'minind', min_ind
  FOR i=0, n_x-1 DO BEGIN
    FOR j=0, n_z-1 DO BEGIN
      dist_pixel = SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)
      IF dist_pixel GT min_ind*1. THEN BEGIN
        min_image_FFT_div_com[i,j] = COMPLEX(0., 0.)
        min_mask_image[i,j] = 0.
      ENDIF
    ENDFOR
  ENDFOR
  
  ;@GJ, 2024/6/14, sorting the plots
  iimage, min_mask_image, title='min mask'
  stat_harmonics_ind = [3,5,7,9,11,13]
  stat_harmonics = stat_harmonics_ind * 0.
  FOR k=0, 5 DO BEGIN
    ind = WHERE(min_mask_image EQ k+1, count)
    stat_harmonics[k] = count
  ENDFOR
  iplot, stat_harmonics_ind, stat_harmonics, color='red', xtitle='Harmonic #', ytitle='Histogram', title='Fourier Selection', /NO_SAVEPROMPT

  keyhole_layer = mask_image * 0.
  rec_layer = mask_image * 0.
  FOR k=0, 5 DO BEGIN
    
    ;@GJ, 2024/6/13, do the line profile
    mask_image = DBLARR(n_x, n_z) * 0. + 100.
    FOR i=0, n_x-1 DO BEGIN
      FOR j=0, n_z-1 DO BEGIN
        dist_pixel = SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)
        min_value = MIN(REFORM(line_profile_FFT[k, *]), min_ind)
        IF dist_pixel LE min_ind*1. THEN mask_image[i,j] = k
      ENDFOR
    ENDFOR
    ind = WHERE(mask_image LE k, count, COMPLEMENT=ind_com)
    IF count GT 0 THEN BEGIN
      temp_focal = REFORM(noisy_image_FFT_div_0[k, *, *])
      temp_focal[ind_com] = COMPLEX(0., 0.)
      keyhole_layer +=  temp_focal
      
      ;@GJ, 2024/6/12, do the image reconstruction of a certain layer
      rec_noisy_image = ABS(FFT(temp_focal, /INVERSE,/center,/double))
      
      ;@GJ, 2024/6/7, display the central frequency with mask
      temp_display =  REFORM(BYTSCL(ALOG(temp_focal),/NAN))
      ;    temp_display[MTF_indices_1st] = 255
      combined_image =  [temp_display, BYTSCL(rec_noisy_image,/NAN)]
      ;    iimage, combined_image, title=STRING(harmonics_n, format='(I0)')+' harmonics'
      total_comb_image[2*n_x:4*n_x-1, (7-k-1)*n_z:(7-k)*n_z-1] = combined_image
      rec_layer += rec_noisy_image
    ENDIF
  ENDFOR
  
  ;@GJ, 2024/6/12, save the keyhole layer
  total_comb_image[2*n_x:3*n_x-1, 0:n_z-1] = REFORM(BYTSCL(ALOG(keyhole_layer),/NAN))
  total_comb_image[3*n_x:4*n_x-1, 0:n_z-1] = BYTSCL(rec_layer,/NAN)

  ;@GJ, 2024/6/7, reconstructed the images
  temp_display =  BYTSCL(ALOG(min_image_FFT_div_com),/NAN)
  rec_noisy_image = ABS(FFT(min_image_FFT_div_com, /INVERSE,/center,/double))
  combined_image =  [temp_display, BYTSCL(rec_noisy_image,/NAN)]
  total_comb_image[0:2*n_x-1, 0:n_z-1] = combined_image

  iimage, total_comb_image[0:2*n_x-1, *], title='Harmonics & Reconstructed'
  psf_filename=parameters.rec_directory+file_basename(parameters.filename)+'_keyhole_Rec_total.jpg'
  IF (FILE_INFO(parameters.rec_directory)).exists EQ 0 THEN FILE_MKDIR, parameters.rec_directory
  WRITE_JPEG, psf_filename, BYTSCL(total_comb_image[0:2*n_x-1, *],/NAN)
  psf_filename_rec=parameters.rec_directory+file_basename(parameters.filename)+'_keyhole_Rec.jpg'
  WRITE_JPEG, psf_filename_rec, BYTSCL(total_comb_image[0:2*n_x-1, 0:n_z-1],/NAN)

END

; docformat = 'rst'
;+
; This is an example program to demonstrate how to create additional Y axes for a plot
; with IDL function graphics routines.
;
; :Categories:
;    Graphics
;
; :Examples:
;    Save the program as "additional_axes_plot_fg.pro" and run it like this::
;       IDL> .RUN additional_axes_plot_fg
;
; :Author:
;    FANNING SOFTWARE CONSULTING::
;       David W. Fanning
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: david@idlcoyote.com
;       Coyote's Guide to IDL Programming: http://www.idlcoyote.com
;
; :History:
;     Change History::
;        Written, 14 January 2014 by David W. Fanning.
;        Written, 1 June 2024 by Guang Jia.

;-
PRO Additional_Axes_Plot_GJ, x, data_1, data_2, data_3, x_title, data_1_title, data_2_title, data_3_title, WINDOW=aWindow

  ; Do this in decomposed color.
  cgSetColorState, 1, CURRENT=currentColorState

  ; Create some data.
;  data_1 = cgScaleVector(cgDemodata(17), 0.0, 1.0)
;  data_2 = cgScaleVector(cgDemodata(17), 0.0, 1000.0)
;  data_3 = (Findgen(101)+1) / 5

  ; Specify colors.
  red = cgColor('red7', /Row, /Triple)
  grn = cgColor('grn7', /Row, /Triple)
  blu = cgColor('blu7', /Row, /Triple)

  ; Leave room in the window for three axes.
  position = [0.15, 0.15, 0.85, 0.85]
  thick = 2

  ; Open a window and return its reference to the user.
  aWindow = Window(WINDOW_TITLE="Resolution vs SNR")

  ; Add the first set of data and the left axes.
  Plot1 = Plot(x, data_1, Color=red, Thick=thick, Position=position, symbol='+', $
    XTITLE=x_title, Ytitle=data_1_title, /Current)
  ay1 = Plot1['Axis 1']
  ay1.Color = red
  ay3 = Plot1['Axis 3']
  ay3.hide = 1

  ; Add the second set of data. Turn off the axes.
  Plot2 = Plot(x, data_2, /Current, AXIS_STYLE=0, Position=position, symbol='tu', $
    Color=blu, LineStyle=1, Thick=thick)

  ; Add the right axis
  Axis2 = Axis('Y', Location=[(Plot2.xrange)[1], 0, 0], Target=Plot2, TextPos=1, $
    Color=blu, Title=data_2_title, TickDir=1)

;  ; Add the third set of data. Turn off the axes
;  Plot3 = Plot(x, data_3, /Current, AXIS_STYLE=0, Position=position, $
;    Color=grn, LineStyle=1, Thick=thick)
;
;  ; Add a third axis and make sure this data is scaled accordiningly.
;  xr = Plot3.xrange
;  Axis3 = Axis('Y', Location=[xr[1]+0.25*(xr[1]-xr[0]),0,0], Target=Plot3, TextPos=1, $
;    Color=grn, Title=data_3_title, /Log, TickDir=1)
;
  ; Incoming color state.
  cgSetColorState, currentColorState

END ;*****************************************************************

;@GJ, 2024/3/22, plotting the figures in manuscript
PRO STED_MPI_simulation
  
  r=FINDGEN(300)/30.
  iplot, r, r^2/(r^2 + 2.^2)^1.5, /NO_SAVEPROMPT

  H_offset_array = FINDGEN(300)/30.
  s_z_array = DBLARR(300,300)
  FOR i=0,299 DO BEGIN
    s_z_array[i,*]  = r^2/(r^2 + H_offset_array[i]^2)^1.5
  ENDFOR
  
;  iimage, s_z_array
  particle_size = 30.;nm
  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  Time = findgen(1000)/1000.     ;(1  ms)
  n_Time = N_ELEMENTS(Time)
  beta = Msat_T*(particle_size^3) /24./1.380469/309.65 ; in 1/mT
;  tau_ms = parameters.tau/1000.
  ;ideal particles
  s_z_image = DBLARR(512, 512) * 0.
  s_z_L_image = DBLARR(512, 512) * 0.
  M_deri = DBLARR(512, 512) * 0.
  H_offset_G = 100.;100.
  FOR i=0, 511 DO BEGIN
    FOR j=0, 511 DO BEGIN
      x = (i - 255)/1.
      y = (j - 255)/1.
      r = x^2 + y^2
      s_z_image[i,j] = r^2/(r^2 + H_offset_G^2)^1.5
      ita = beta * SQRT(r^2 + H_offset_G^2)
      M_deri[i,j] = (1./ita^2 - 1./SINH(ita)^2)
      s_z_L_image[i,j] = M_deri[i,j] * s_z_image[i,j]
    ENDFOR
  ENDFOR
  iimage, s_z_image, title='sz image', /NO_SAVEPROMPT
  iplot, s_z_image[255,*], color='blue', thick=3, title='Ideal', xrange=[0,512], /NO_SAVEPROMPT
  iimage, s_z_L_image, title='sz Langevin image', /NO_SAVEPROMPT
  iplot, s_z_L_image[255,*], color='green', thick=3, title='Langevin', xrange=[0,512], /NO_SAVEPROMPT
  iimage, M_deri, title='M deri', /NO_SAVEPROMPT
  
  chebyshev_array = DBLARR(512)*0.
  n=3.
  H_offset_G = 100.;100.
  H_offset_A_array = (FINDGEN(512)-256)/256.*1.5
  FOR i=0, 511 DO BEGIN
    IF ABS(H_offset_A_array[i]) LT 1.0 THEN chebyshev_array[i] = SIN(n*ACOS(H_offset_A_array[i]));/SIN(ACOS(H_offset_A_array[i]))*SQRT(1-H_offset_A_array[i]^2));Chebyshev
  ENDFOR
  iplot, H_offset_A_array, chebyshev_array, /NO_SAVEPROMPT
;  print, 'test'
  
  s_array = DBLARR(512, 512, 512)
  base_array = DBLARR(512, 512, 512)
  H_offset_G_array = (FINDGEN(512)-256)/256.*H_offset_G*1.5
  M_deri_array = DBLARR(512, 512, 512) * 0.
  H_z_comp_array = DBLARR(512, 512, 512) * 0.
  FOR i=0, 511 DO BEGIN
    print, 'i,j', i,j
    FOR j=0, 511 DO BEGIN
      FOR k=0, 511 DO BEGIN
        x = (i - 255)/1.
        y = (j - 255)/1.
        r = x^2 + y^2
        H_z_comp_array[i,j,k] = r^2/(r^2 + H_offset_G_array[k]^2)^1.5
        ita = SQRT(r^2 + H_offset_G_array[k]^2)
        M_deri_array[i,j,k] = (1./ita^2 - 1./SINH(ita)^2)
        base_array[i,j,k] = M_deri_array[i,j,k] * H_z_comp_array[i,j,k]
      ENDFOR
      temp_conv = CONVOL(REFORM(base_array[i,j,*]), REFORM(chebyshev_array), /center, /EDGE_ZERO, /NORMALIZE)
      s_array[i,j,*] = temp_conv
    ENDFOR
  ENDFOR
  
  print, H_offset_A_array[256*1.8], 'HDC/HAC'
  iimage, REFORM(s_array[*,*,256*1.8]),title='conv'
  iimage, ALOG(REFORM(s_array[256,*,*])),title='ALOG(s_3_z)'
  print, 'test'
  
END


function Debye_Kernel_sin_FMMPI, tau, kernel_time
  N = N_ELEMENTS(kernel_time)
  kernel = DBLARR(N)*0.
  d_t = kernel_time[1]-kernel_time[0] ;us
  time = FINDGEN(N) * d_t  ; us
  r = (1./tau) * exp(-time/tau)
  kernel[0:N/2] = REVERSE(r[0:N/2])

  return, kernel
end

FUNCTION ReverseIndices_FMMPI, ri, index, COUNT=count

  Compile_Opt idl2

  ; Error handling.
  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /CANCEL
;    void = Error_Message()
    count = 0
    RETURN, -1
  ENDIF

  ; Need two parameters.
  IF N_Params() NE 2 THEN Message, 'Two positional parameters, REVERSE_INDICES and INDEX, are required.'

  ; Return the indices, if there are any.
  IF ri[index] NE ri[index+1] THEN BEGIN
    indices = ri[ri[index]:ri[index+1]-1]
    count = N_Elements(indices)
  ENDIF ELSE BEGIN
    indices = -1
    count = 0
  ENDELSE

  RETURN, indices

END

FUNCTION AddPoissonNoise_FMMPI, image, PERCENTNOISE=percentNoise, SEED=seed
  SetDefaultValue, image, cgDemoData(20)
  SetDefaultValue, percentNoise, 5.0
  h = cgHistogram(image, LOCATIONS=xbin, REVERSE_INDICES=ri, NBINS=100)
  noiseImage = image
  FOR j=0,N_Elements(h)-1 DO BEGIN
    currentMean = xbin[j]
;    print, 'currentMean = ', currentMean
    indices = ReverseIndices_FMMPI(ri, j, COUNT=count)
    IF count GT 0 THEN BEGIN
;      print, 'POISSON: ', currentMean
      IF currentMean EQ 0 THEN currentMean++
      noiseImage[indices] = RandomU(seed, count, POISSON=currentMean, /DOUBLE)
      newMean = SqRt(currentMean) / (percentNoise/100.0)

      ; Adjust the standard deviation to the required level.
      noiseImage[indices] += newMean - currentMean

      ; Now adjust it back to starting the mean, but with noise added.
      noiseImage[indices] *= currentMean / newMean
    ENDIF
  ENDFOR
  RETURN, noiseImage
END

FUNCTION brownian_time_calc, parameters
  
  eta = (parameters.viscosity) / 1000. ; Pa.s
  D_h = parameters.particle_size ; in nm
  print, 'D_h: ', D_h, ' nm'
  D_h *= 1.e-9 ; m
  V_h = 1./6. * !PI * D_h^3
  k_B = 1.380469e-23; J/K
  T_p = parameters.T_p; degree
  T_p_kelvin = 273.15 + T_p ; in Kelvin

  tau_B = 3. * eta * V_h / (k_B * T_p_kelvin) * 1.e6 ; in us
  print, 'tao=u_B: ', tau_B, ' us'
;  parameters.tau = tau_B
  RETURN, tau_B; brwonian_relaxation_time
END

;@GJ, 2024/9/10, calculate the magnetic field dependent relaxation time
FUNCTION fields_brownian_time_calc, tau_B, beta, H_AC, H_DC_par, H_DC_perp
  xi_H_AC = beta * H_AC
  coeff_H_AC = 1./SQRT(1. + 0.126 * xi_H_AC^1.72)
  
  xi_H_DC_par = beta * H_DC_par
  IF xi_H_DC_par NE 0 THEN BEGIN
    Lprime_H_DC_par = 1./(xi_H_DC_par^2) - 1./(sinh(xi_H_DC_par)^2)
    L_H_DC_par = 1./tanh(xi_H_DC_par) - 1./xi_H_DC_par
    coeff_H_DC_par = xi_H_DC_par * Lprime_H_DC_par / L_H_DC_par
  ENDIF ELSE BEGIN
    coeff_H_DC_par = 1.
  ENDELSE
  
  xi_H_DC_perp = beta * H_DC_perp
  IF xi_H_DC_perp NE 0 THEN BEGIN
    L_H_DC_perp = 1./tanh(xi_H_DC_perp) - 1./xi_H_DC_perp
    coeff_H_DC_perp = 2. * L_H_DC_perp / (xi_H_DC_perp - L_H_DC_perp)
  ENDIF ELSE BEGIN
    coeff_H_DC_perp = 1.
  ENDELSE
  
  tau_B_H = tau_B * coeff_H_AC * coeff_H_DC_perp; * coeff_H_DC_par
  RETURN, tau_B_H; brwonian_relaxation_time
END


;@GJ, 2024/4/18, optimization of the dog factor
PRO Dog_factor_optimization, wid, parameters
  ;@GJ, 2024/4/18, find the optimal dog_factor
  dog_factor_array = (FINDGEN(100)+1.)/10.
;  STED_chisq_array = DOUBLE(dog_factor_array * 0.)
  error_fit_array = DOUBLE(dog_factor_array * 0.)
  ele_fit_array = DOUBLE(dog_factor_array * 0.)
  
  min_0 = MIN(ABS((*parameters.focal_line)[0:parameters.n_z/2-1]-0.5*MAX(*parameters.focal_line)), half_ind)
  FOR i_dog=0, N_ELEMENTS(dog_factor_array)-1 DO BEGIN
    temp_STED_PSF = *parameters.focal_line-*parameters.focal_line_1*dog_factor_array[i_dog]
    error_fit_array[i_dog] = MEAN(ABS(temp_STED_PSF[0:half_ind]))
    IF MIN(temp_STED_PSF) LT 0 THEN ele_fit_array[i_dog] = MAX(*parameters.focal_line) ELSE ele_fit_array[i_dog] = temp_STED_PSF[0]
  ENDFOR
  min_error = MIN(error_fit_array, min_ind_1)
  min_ele = MIN(ele_fit_array, min_ind_2)
  min_ind = (min_ind_1 + min_ind_2) / 2. 
  ;    iplot, STED_chisq_array/MAX(STED_chisq_array), /NO_SAVEPROMPT
;  iplot, dog_factor_array, STED_chisq_array/MAX(STED_chisq_array)*MAX(FWHM_array), xtitle='dog factor', ytitle='Chisq', title='Chisq'
;  iplot, dog_factor_array, FWHM_array, /overplot, color='red';xtitle='dog factor', ytitle='FWHM [mT]', title='Resolution vs Dog Factor'
  ;    min_FWHM = MIN(FWHM_array, min_ind)
;  min_STED_chisq_array = MIN(STED_chisq_array, min_ind)
  
;  yfit = GAUSSFIT(*parameters.A_G_array, *parameters.focal_line-*parameters.focal_line_1*dog_factor_array[min_ind], coeff, CHISQ=temp_chisq, NTERMS=4)
;  print, 'optimal dog factor and resolution: ',  dog_factor_array[min_ind], FWHM_array[min_ind]
;  iplot, *parameters.A_G_array, *parameters.focal_line-*parameters.focal_line_1*dog_factor_array[min_ind], title='Minimal STED PSF', xtitle='A_G [mT]'
;  iplot, *parameters.A_G_array, yfit, color='red', /overplot
;  iplot, *parameters.A_G_array, *parameters.focal_line/MAX(*parameters.focal_line)*coeff[0], color='blue', /overplot
  *parameters.focal_STED = *parameters.focal_line-*parameters.focal_line_1*dog_factor_array[min_ind]
  parameters.optimal_dog_factor = dog_factor_array[min_ind] * 10
  
  ;@GJ, 2024/4/18
  widget_dog_factor = widget_info(wid, find_by_uname='dog_factor')
  widget_control, widget_dog_factor, set_value = parameters.optimal_dog_factor
  
  parameters.dog_factor = dog_factor_array[min_ind] * 10
END

;@GJ, 2024/4/21, optimization of the reconstructed image
PRO Rec_sub_threshold_optimization, wid, parameters
  
  ;create the progress bar
  progressbar = Obj_New('progressbar', Color='cyan', Text='0%', Title='Optimizing...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  
  SSIM_array = DBLARR(256) * 0.
  ; Place the progress bar on the display.
  progressbar -> Start
  
  FOR i=0, N_ELEMENTS(SSIM_array)-1 DO BEGIN
    ;update the progress bar
    count=(i+1.)/N_ELEMENTS(SSIM_array)*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    
    temp_image_sub = BYTSCL(*parameters.rec_noisy_image_sub)
    temp_image_sub[WHERE(temp_image_sub LT i, /null)] = 0.
    SSIM_array[i] = SSIM(BYTSCL(*parameters.image_sim_ori), BYTSCL(temp_image_sub))
  ENDFOR
  
  ;@GJ, 2024/4/21, find the largest ssim
  max_ssim = MAX(SSIM_array, max_ind)
  IF max_ind GT 250. THEN max_ind = 0.;
;  iplot, SSIM_array, title='SSIM array', /NO_SAVEPROMPT
  
  parameters.phantom_rec_sub_threshold_opt = max_ind
  parameters.phantom_rec_sub_threshold =  max_ind
  
  ;@GJ, 2024/4/18
  widget_threshold = widget_info(wid, find_by_uname='threshold')
  widget_control, widget_threshold, set_value = parameters.phantom_rec_sub_threshold_opt
  
  ;destroy the progress bar
  progressbar -> Destroy
  
END

;@GJ, 2024/5/3, write log and parameters
PRO write_results_FM_MPI, parameters
;  psf_filename=parameters.rec_directory+'log.xls'
  current_time = string(SYSTIME(/JULIAN, /UTC), FORMAT='(f8.0)')
  output_filename = parameters.rec_directory+'log_STED.xls'
  ;write file
  IF (file_info(output_filename)).exists EQ 0 THEN BEGIN
    openw,unit,output_filename,/get_lun
    first_line =['Particle Type', 'Particle Size [nm]', 'Harmonic #', 'Harmonic Type 0 Amp', 'AC [mT]', 'DC Donut Threshold [%]', 'Initial Donut Radius [mT]', 'Donut Slope [mT/mT]']
    printf, unit, FORMAT = '(10(A, %"\t"))', first_line
  ENDIF ELSE BEGIN
    openw,unit,output_filename,/get_lun, /append
  ENDELSE
  printf, unit, FORMAT = '(10(A, %"\t"))', parameters.particle_type, parameters.particle_size, parameters.harmonics_n, parameters.harmonics_comp, parameters.A, parameters.A_DC_threshold_donut, parameters.donut_radius0, parameters.donut_slope
  ;; close file and END execution
  FREE_LUN, unit
END


PRO Signal_calculation_FM_MPI, parameters
  
;  print, 'parameters.rec_directory: ', parameters.rec_directory
  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  Time = findgen(1000)/1000.     ;(1  ms)
  n_Time = N_ELEMENTS(Time)
  beta = Msat_T*((parameters.particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
  ;  print, 'beta: ', beta
  ;@GJ, 2023/12/26, calculate the period in ms
  T_period = 1./parameters.f; in ms
  n_T_period = T_period / 0.001

  ;@GJ, 2023/12/21, generate PSF
  delta_x = parameters.FOV/parameters.n_x
  delta_z = parameters.FOV/parameters.n_z
;  gradient_field_x = parameters.gradient_field_x * 0.1 ;mT/mm
  gradient_field_x = parameters.A * 2. * (parameters.A_DC_max/100.) / parameters.FOV ;calculating the gradient field
  gradient_field_yz = parameters.gradient_field_yz * 0.1 ;mT/mm
  print, 'gradient_x: ', gradient_field_x
  print, 'gradient_yz: ', gradient_field_yz

  kernel_Gx_A = FINDGEN(parameters.harmonics_n_max, parameters.n_x) * 0.
  FOR n=0, parameters.harmonics_n_max-1 DO BEGIN
    FOR k=0, parameters.n_z-1 DO BEGIN
      Gx_A = (k - parameters.n_z/2.) * delta_x * gradient_field_x / parameters.A
      IF ABS(Gx_A) LT 1.0 THEN BEGIN
        kernel_Gx_A[n, k] = SIN((n+1.)*ACOS(-Gx_A));ABS(SIN((n+1.)*ACOS(Gx_A))/SIN(ACOS(Gx_A))*SQRT(1-Gx_A^2));Chebyshev
      ENDIF
    ENDFOR
  ENDFOR
  ;  iimage, kernel_Gz_A

  Sx_array = DBLARR(parameters.harmonics_n_max+2, parameters.harmonics_comp_N, parameters.n_x, parameters.n_z) * 0.
  Sx_array_kernel = DBLARR(parameters.harmonics_n_max+2, parameters.harmonics_comp_N, parameters.n_x, parameters.n_z) * 0.
  Sx_array_xzNonCOnv = DBLARR(parameters.harmonics_n_max+2, parameters.harmonics_comp_N, parameters.n_x, parameters.n_z) * 0.
  IF parameters.AC_DC_mode EQ 0. THEN BEGIN
    ;AC excitation mode
    x_array = (FINDGEN(parameters.n_x))/(parameters.n_x) * (parameters.AC_ex_max - parameters.AC_ex_min) * parameters.A + parameters.AC_ex_min * parameters.A
  ENDIF ELSE BEGIN
    ;DC offset mode
    x_array = (FINDGEN(parameters.n_x)-parameters.n_x/2.) * delta_x * gradient_field_x
  ENDELSE
  
  z_array = (FINDGEN(parameters.n_z)-parameters.n_z/2.) * delta_z * gradient_field_yz
  
  IF parameters.particle_type LE 2 THEN BEGIN
    ;      ; Create the progress bar.
    IF parameters.particle_type EQ 0 THEN progressbar = Obj_New('progressbar', Color='pink', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    IF parameters.particle_type EQ 1 THEN progressbar = Obj_New('progressbar', Color='red', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    IF parameters.particle_type EQ 2 THEN progressbar = Obj_New('progressbar', Color='yellow', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
      ; Place the progress bar on the display.
    progressbar -> Start
    FOR i=0, parameters.n_z-1 DO BEGIN
      count=(i+1.)/parameters.n_z*100.0
      progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
      ;0: Ideal chebyshev
      IF parameters.particle_type EQ 0 THEN base_array = Msat_kAm/T_period*(z_array[i])^2/((z_array[i])^2 + x_array^2)^1.5
      ;1: SFMIO chebyshev
      IF parameters.particle_type EQ 1 THEN BEGIN
        ;@GJ, 2023/12/24, calculating SFMIO particles from small method
        H_coercivity = parameters.H_coercivity; mT
        part_1 = Msat_kAm/T_period*(z_array[i])^2/ABS((SQRT((z_array[i])^2 + x_array^2)-H_coercivity))^3
        part_2 = Msat_kAm/T_period*(z_array[i])^2/ABS((SQRT((z_array[i])^2 + x_array^2)+H_coercivity))^3
        base_array =  (part_1 + part_2) / 2.
      ENDIF
      ;2: Langevin chebyshev
      IF parameters.particle_type EQ 2 THEN BEGIN
        ;@GJ, 2023/12/24, calculate the absolute of H
        base_array = x_array * 0.
        H_abs_array = x_array * 0.
        zero_k = -1
        zero_i = -1
        FOR k=0, parameters.n_x-1 DO BEGIN
          H_abs = SQRT(z_array[i]^2 + x_array[k]^2)
          H_abs_array[k] = H_abs
          IF H_abs NE 0 THEN BEGIN
            part_1 = -2.*beta*Msat_kAm/T_period*(1./(beta*H_abs)^2 - 1./(sinh(beta*H_abs)^2))*x_array[k]^2/H_abs^2
            part_2 = -2.*Msat_kAm/T_period*(1/tanh(beta*H_abs) - 1/(beta*H_abs)) * (z_array[i])^2/H_abs^3
            IF ABS(part_1/part_2) GT 0.5 THEN print, 'beta part/2nd part = ', part_1/part_2
;            part_3 = -2.*Msat_kAm/T_period*(1./(beta*z_array[i])^2 - 1./(sinh(beta*z_array[i])^2))*z_array[k]^2/H_abs^3
          ENDIF ELSE BEGIN
            zero_k = k
            zero_i = 1
            ;@GJ, 2023/12/24, if the H is zero, do special analysis
            part_1 = 0.;-2.*beta*Msat_kAm/T_period*(1./3.)*x_array[k]^2/H_abs^2
            part_2 = 0.; * (z_array[i])^2/((z_array[i])^2 + x_array[k]^2)^1.5
;            part_3 = -2.*Msat_kAm/T_period*(1./3.);*z_array[k]^2/H_abs^3
          ENDELSE
          base_array[k] = part_1 + part_2
        ENDFOR
        
        ;@GJ, 2024/3/27, do interpolate
        IF zero_i GT 0 THEN BEGIN
          base_array[zero_k] = (base_array[zero_k-1]+base_array[zero_k+1])/2.
          base_array[zero_k] = INTERPOL(base_array, FINDGEN(parameters.n_x), zero_k, /LSQUADRATIC)
        ENDIF
        
      ENDIF
      FOR n=0, parameters.harmonics_n_max-1 DO BEGIN
        ;signal
        temp_signal = CONVOL(base_array, REFORM(kernel_Gx_A[n, *]), /center, /EDGE_ZERO, /NORMALIZE)
        (*parameters.Sx_array)[n, 0, i, *] = ABS(temp_signal) ;Amp
        temp_phase = temp_signal/ABS(temp_signal)*90.*180./!PI
        IF temp_phase LT 0 THEN temp_phase += 360.
        (*parameters.Sx_array)[n, 1, i, *] = temp_phase ; Phase
        (*parameters.Sx_array)[n, 2, i, *] = 0. ; Real
        (*parameters.Sx_array)[n, 3, i, *] = temp_signal ;Imaginary
        
        ;kernel
        (*parameters.Sx_array_kernel)[n, 0, i, *] = REFORM(kernel_Gx_A[n, *])
        (*parameters.Sx_array_kernel)[n, 1, i, *] = 0.
        (*parameters.Sx_array_kernel)[n, 2, i, *] = 0.
        (*parameters.Sx_array_kernel)[n, 3, i, *] = REFORM(kernel_Gx_A[n, *])
        ;nonconv
        (*parameters.Sx_array_xzNonCOnv)[n, 0, i, *] = ABS(base_array)
        (*parameters.Sx_array_xzNonCOnv)[n, 1, i, *] = 0.
        (*parameters.Sx_array_xzNonCOnv)[n, 2, i, *] = 0.
        (*parameters.Sx_array_xzNonCOnv)[n, 3, i, *] = base_array
      ENDFOR

      Amp_dH_dt = parameters.A * 2. * !PI * parameters.f
      (*parameters.Sx_array)[parameters.harmonics_n_max+1, 0, i, *] = 0;-(base_array * Amp_dH_dt); Unfiltered signal
      
;      ;GJ, 2024/4/1, calculating the filtered signal
;      signal_new_x = -(base_array * Amp_dH_dt)
;      t_ele = N_Elements(signal_new_x)
;      delta_t = (Time[1]-Time[0]) * DOUBLE(parameters.f) ; us; samplying rate: 1MHz
;      X = (FINDGEN((t_ele - 1)/2) + 1)
;      is_N_even = (t_ele MOD 2) EQ 0
;      if (is_N_even) then $
;        frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
;      else $
;        frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
;      base_freq = MIN(ABS(frequency-1.0), base_ind)
;      signal_FFT_x = FFT(signal_new_x)
;      signal_FFT_x[base_ind] = COMPLEX(0., 0.)
;      signal_FFT_x_filtered = FFT(signal_FFT_x, /INVERSE)
      (*parameters.Sx_array)[parameters.harmonics_n_max, 0, i, *] = 0; ABS(signal_FFT_x_filtered);filtered signal
    ENDFOR
    ;destroy the progress bar
    progressbar -> Destroy
  ENDIF
  
  ;@GJ, 2023/12/25, for simulating particles
  IF parameters.particle_type GT 2 THEN BEGIN
    ;      ; Create the progress bar.
    IF parameters.particle_type EQ 3 THEN progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Simulating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    IF parameters.particle_type EQ 4 THEN progressbar = Obj_New('progressbar', Color='cyan', Text='0%', Title='Simulating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    IF parameters.particle_type EQ 5 THEN progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Simulating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    ;      ; Place the progress bar on the display.
    progressbar -> Start
    flag=0
    FOR ld=0, parameters.n_x-1 DO BEGIN
      count=(ld+1.)/parameters.n_x*100.0
      progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
      
      ;print, 'A, A_DC=', parameters.A, ',', A_DC_ld, '(mT)'
      ;IF ABS(A_DC_ld) GT 1.5*A THEN continue
      FOR md=0, parameters.n_z-1 DO BEGIN
        H_G = z_array[md]
        ;IF ABS(H_G) GT 1.5*A THEN continue
        IF parameters.AC_DC_mode EQ 0. THEN BEGIN
          ;AC ex mode
          AC_ex_1d = x_array[ld]
          H_x = AC_ex_1d*cos(2.*!PI*parameters.f*Time)
        ENDIF ELSE BEGIN
          ;DC offset mode
          A_DC_ld = x_array[ld]
          H_x = parameters.A*cos(2.*!PI*parameters.f*Time) + A_DC_ld
        ENDELSE
        ;@GJ, 2024/9/17, modify the A_DC to A.
        
        ;            H_x = A*sin(2.*!PI*f*Time) + A_DC
        H_total = SQRT(H_x^2 + H_G^2) ;SIGNUM(H_x) * 
        
        ;Langevin particlee
        IF parameters.particle_type EQ 3 OR parameters.particle_type EQ 4 THEN BEGIN
          ;calculate M
          M_total = Msat_kAm*(1/tanh(beta*H_total) - 1/(beta*H_total)) ;(kA/m)
          index_nan = where(~FINITE(M_total), count)
          IF count GT 0 THEN BEGIN
            RESULT = INTERPOL(M_total, Time, Time[index_nan], /NAN)
            M_total[index_nan] = RESULT
          ENDIF

          ;do vector splitting
          ;@GJ, 2023/9/23, simple calculation with the signal intensity
          ;@GJ, 2024/8/11, do the calculation based on the Jie He's MatLab code
          M_x = M_total * H_x / H_total
          ;@GJ, 2024/9/9, calculate the adiabatic signal
          signal_x_adiabatic = -(M_x[1:*] - M_x[0:n_Time-2])/(Time[1]-Time[0])
          
          ;@GJ, 2024/4/10, Moving Debye convolution to the magnetization
          IF parameters.particle_type EQ 4 THEN BEGIN
;            tau_ms = relaxation_time_Resovist(ABS(H_G), parameters.f)
;            tau_us = relaxation_pulsed_excitation(ABS(H_G), parameters.particle_size)
;            tau_ms = tau_us / 1000.
;            tau_ms = parameters.tau/1000.
            IF parameters.AC_DC_mode EQ 0. THEN BEGIN
              tau_ms = fields_brownian_time_calc(parameters.tau/1000., beta, AC_ex_1d, 0., ABS(H_G))
            ENDIF ELSE BEGIN
              tau_ms = fields_brownian_time_calc(parameters.tau/1000., beta, parameters.A, ABS(A_DC_ld), ABS(H_G))
            ENDELSE
;            ;@GJ, 2024/9/17, remove the A_DC_ld
            (*parameters.tau_ms_array)[md, ld] = tau_ms 
;            IF md EQ parameters.n_x/2 THEN print, 'tau_ms: ', tau_ms, " ms"
            IF tau_ms GT 0.00001 THEN BEGIN
              R = (1/tau_ms)*exp(-Time/tau_ms)
              kernel = Debye_Kernel_sin_FMMPI(tau_ms, Time[0:100])
              M_x_D = convol(M_x, kernel, /EDGE_MIRROR, /NORMALIZE) ;@GJ, 2024/8/25, add normalize to make the convol better
              M_x = M_x_D
            ENDIF
          ENDIF
          
          ;@GJ, 2024/4/27, calculate the signal along x direction
          signal_x = -(M_x[1:*] - M_x[0:n_Time-2])/(Time[1]-Time[0])
          signal_new_x = signal_x
        ENDIF
        
        ;SFMIOs simulation, @GJ, 2024/04/27
        IF parameters.particle_type EQ 5 THEN BEGIN
          ;initialize the signal
          M_x = H_x * 0.
          signal_new_x = H_x[1:*] * 0.
          signal_x = signal_new_x

          ;@GJ, both surpassing the coercivity
          IF (-ABS(parameters.A)+A_DC_ld LT -parameters.H_coercivity) AND (ABS(parameters.A)+A_DC_ld GT parameters.H_coercivity) THEN BEGIN
            FOR i=0, FLOOR(n_Time/n_T_period-1) DO BEGIN
              ;calculate the 1st signal location
              min_1st = MIN(ABS(H_x[i*n_T_period:i*n_T_period+n_T_period/2-1]+parameters.H_coercivity), min_ind_1st)
              signal_new_x[i*n_T_period+min_ind_1st] = Msat_kAm * 2./T_period * ABS(H_x[i*n_T_period+min_ind_1st])/ABS(H_total[i*n_T_period+min_ind_1st])
              ;calculate the 2nd signal location
              min_2nd = MIN(ABS(H_x[i*n_T_period+n_T_period/2:i*n_T_period+n_T_period-1]-parameters.H_coercivity), min_ind_2nd)
              signal_new_x[i*n_T_period+n_T_period/2+min_ind_2nd] = -Msat_kAm * 2./T_period * ABS(H_x[i*n_T_period+n_T_period/2+min_ind_2nd])/ABS(H_total[i*n_T_period+n_T_period/2+min_ind_2nd])
            ENDFOR
            ;
            IF A_DC_ld EQ parameters.H_coercivity THEN BEGIN
              IF flag EQ 0 THEN iplot, signal_new_x, /NO_SAVEPROMPT ELSE iplot, signal_new_x, /overplot
              flag++
            ENDIF
          ENDIF
        ENDIF
        
        ;@GJ, 2022/12/14, removing the NaN points
        index_nan_x = where(~FINITE(signal_x), count_x)
        IF count_x GT 0 AND count_x LT n_Time/10. THEN BEGIN
          RESULT_x = INTERPOL(signal_x, Time[0:n_Time-2], Time[index_nan_x], /NAN)
          index_nan_Result_x = where(~FINITE(Result_x), count_Result_x)
          IF count_Result_x EQ 0 THEN BEGIN
            signal_new_x[index_nan_x] = RESULT_x
            ;@GJ, 2024/9/9, calculate the adiabatic signal
            signal_x_adiabatic[index_nan_x] = RESULT_x
          ENDIF
        ENDIF
;        IF parameters.particle_type EQ 3 THEN BEGIN
;          signal_new_x = signal_x
;        ENDIF
;        ;@GJ, 2024/4/5, removing a bug about too small tau
;        IF parameters.particle_type EQ 4 AND tau_ms GT 0.001 THEN BEGIN
;          R = (1/tau_ms)*exp(-Time/tau_ms)
;          kernel = Debye_Kernel_sin_FMMPI(tau_ms, Time[0:100])
;          signal_new_x = convol([dblarr(n_elements(kernel)/2)*0, signal_x, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
;          signal_new_x = signal_new_x[n_elements(kernel)/2 : n_elements(signal_x) + n_elements(kernel)/2]
;        ENDIF ELSE BEGIN
;          signal_new_x = signal_x
;        ENDELSE

        t_ele = N_Elements(signal_new_x)
        delta_t = (Time[1]-Time[0]) * DOUBLE(parameters.f) ; us; samplying rate: 1MHz
        X = (FINDGEN((t_ele - 1)/2) + 1)
        is_N_even = (t_ele MOD 2) EQ 0
        if (is_N_even) then $
          frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
        else $
          frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
        base_freq = MIN(ABS(frequency-1.0), base_ind)
        
        ;@GJ, 2024/8/16, add noise to the signal
        noise_level_percentage = 1. / (10.^(parameters.noise_level/10.)) * 100.
;        print, 'noise level =  ', noise_level_percentage, '%'
        IF noise_level_percentage GT 1 THEN BEGIN
          min_signal = MIN(signal_new_x)
          noisy_signal_new_x = min_signal + AddPoissonNoise_FMMPI(signal_new_x-min_signal, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(bar_signal_1d_0), LEVELS=[noise_level_percentage, noise_level_percentage]);/255.*MAX(bar_signal_1d_0); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
          ;@GJ, 2024/9/9, calculate the noisy adiabatic signal
          noisy_signal_x_adiabatic = min_signal + AddPoissonNoise_FMMPI(signal_x_adiabatic-min_signal, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(bar_signal_1d_0), LEVELS=[noise_level_percentage, noise_level_percentage]);/255.*MAX(bar_signal_1d_0); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
        ENDIF ELSE BEGIN
          noisy_signal_new_x = signal_new_x
          ;@GJ, 2024/9/9, calculate the noisy adiabatic signal
          noisy_signal_x_adiabatic = signal_x_adiabatic
        ENDELSE
        signal_FFT_x = FFT(noisy_signal_new_x)
        signal_FFT_x_adiabatic = FFT(noisy_signal_x_adiabatic)

        ;signal_FFT_x = FFT(signal_new_x)
        
;        ;@GJ, 2024/8/11, plot the result
;        IF ABS(H_G) LT 0.001 THEN BEGIN
;          IF ABS(ABS(A_DC_ld) - parameters.A*2./3.) LT 0.1 OR ABS(A_DC_ld) LT 0.01 THEN BEGIN
;            iplot, Time, H_x, title='H', xtitle='Time [ms]', ytitle='H [mT]', /NO_SAVEPROMPT
;            iplot, Time[1:*], signal_new_x, title='Signal', xtitle='Time [ms]', ytitle='Signal [A.U.]', /NO_SAVEPROMPT
;;            iplot, frequency/base_ind, ABS(signal_FFT_x), xrange=[0, 10], title='FFT', xtitle='Freq', ytitle='FFT', /NO_SAVEPROMPT
;            iplot, frequency/frequency[base_ind], 180./!PI*ATAN(signal_FFT_x, /phase), xrange=[0, 10], title='FFT', xtitle='Freq', ytitle='Angle', /NO_SAVEPROMPT
;            print, 'md ', md, '; ld ', ld
;            print, 'H_G ', z_array[md], '; H_DC ', x_array[ld], '; 3rd Phase ', 180./!PI*ATAN(signal_FFT_x[base_ind*2.7], /phase)
;          ENDIF
;        ENDIF
          
        ;@GJ, 2023/10/20, get the harmonic amp and phase
        FOR n=0, parameters.harmonics_n_max-1 DO BEGIN
          IF base_ind*(n+1) LT N_ELEMENTS(signal_FFT_x) THEN BEGIN
            (*parameters.Sx_array)[n, 0, md, ld] = ABS(signal_FFT_x[base_ind*(n+1)])
            ;@GJ, 2024/1/5, correcting the phase calculation
            ;(*parameters.Sx_array)[n, 1, md, ld] = 180./!PI*ATAN(Imaginary(signal_FFT_x[base_ind*(n+1)]), Real_part(signal_FFT_x[base_ind*(n+1)]))
            ;IF (*parameters.Sx_array)[n, 1, md, ld] GT 180. THEN (*parameters.Sx_array)[n, 1, md, ld] = 360. - (*parameters.Sx_array)[n, 1, md, ld]
            ;IF (*parameters.Sx_array)[n, 1, md, ld] LT -180. THEN (*parameters.Sx_array)[n, 1, md, ld] = 360. + (*parameters.Sx_array)[n, 1, md, ld]

            ;@GJ, 2024/8/11, calculate the phase in a standard way
            temp_phase = 180./!PI*ATAN(signal_FFT_x[base_ind*(n+1)], /phase)
            IF temp_phase LT 0 THEN temp_phase += 360.
            (*parameters.Sx_array)[n, 1, md, ld] = temp_phase
            (*parameters.Sx_array)[n, 2, md, ld] = Real_part(signal_FFT_x[base_ind*(n+1)])
            (*parameters.Sx_array)[n, 3, md, ld] = Imaginary(signal_FFT_x[base_ind*(n+1)])

            ;@GJ, save the adiabatic part
            (*parameters.Sx_adiabatic_array)[n, 0, md, ld] = ABS(signal_FFT_x_adiabatic[base_ind*(n+1)])
            temp_phase = 180./!PI*ATAN(signal_FFT_x_adiabatic[base_ind*(n+1)], /phase)
            IF temp_phase LT 0 THEN temp_phase += 360.
            (*parameters.Sx_adiabatic_array)[n, 1, md, ld] = temp_phase
            (*parameters.Sx_adiabatic_array)[n, 2, md, ld] = Real_part(signal_FFT_x_adiabatic[base_ind*(n+1)])
            (*parameters.Sx_adiabatic_array)[n, 3, md, ld] = Imaginary(signal_FFT_x_adiabatic[base_ind*(n+1)])
          ENDIF
        ENDFOR
        
        ;@GJ, 2024/8/12, save any harmonic
        IF parameters.harmonics_n-FLOOR(parameters.harmonics_n) NE 0 THEN BEGIN
          harmonics_n = FLOOR(parameters.harmonics_n)
          (*parameters.Sx_array)[harmonics_n-1, 0, md, ld] = ABS(signal_FFT_x[base_ind*(parameters.harmonics_n)])
           temp_phase = 180./!PI*ATAN(signal_FFT_x[base_ind*(parameters.harmonics_n)], /phase)
          IF temp_phase LT 0 THEN temp_phase += 360.
          (*parameters.Sx_array)[harmonics_n-1, 1, md, ld] = temp_phase
          (*parameters.Sx_array)[harmonics_n-1, 2, md, ld] = Real_part(signal_FFT_x[base_ind*(parameters.harmonics_n)])
          (*parameters.Sx_array)[harmonics_n-1, 3, md, ld] = Imaginary(signal_FFT_x[base_ind*(parameters.harmonics_n)])
          
          ;@GJ, save the adiabatic part
          (*parameters.Sx_adiabatic_array)[harmonics_n-1, 0, md, ld] = ABS(signal_FFT_x_adiabatic[base_ind*(parameters.harmonics_n)])
          temp_phase = 180./!PI*ATAN(signal_FFT_x_adiabatic[base_ind*(parameters.harmonics_n)], /phase)
          IF temp_phase LT 0 THEN temp_phase += 360.
          (*parameters.Sx_adiabatic_array)[harmonics_n-1, 1, md, ld] = temp_phase
          (*parameters.Sx_adiabatic_array)[harmonics_n-1, 2, md, ld] = Real_part(signal_FFT_x_adiabatic[base_ind*(parameters.harmonics_n)])
          (*parameters.Sx_adiabatic_array)[harmonics_n-1, 3, md, ld] = Imaginary(signal_FFT_x_adiabatic[base_ind*(parameters.harmonics_n)])
        ENDIF   
        signal_FFT_x[base_ind] = COMPLEX(0., 0.)
        signal_FFT_x_filtered = FFT(signal_FFT_x, /INVERSE)
        (*parameters.Sx_array)[parameters.harmonics_n_max, 0, md, ld] = MAX(ABS(signal_FFT_x_filtered))
        (*parameters.Sx_array)[parameters.harmonics_n_max+1, 0, md, ld] = MAX(ABS(noisy_signal_new_x)); Unfiltered signal
        
        ;@GJ, save the adiabatic part
        signal_FFT_x_adiabatic[base_ind] = COMPLEX(0., 0.)
        signal_FFT_x_adiabatic_filtered = FFT(signal_FFT_x_adiabatic, /INVERSE)
        (*parameters.Sx_adiabatic_array)[parameters.harmonics_n_max, 0, md, ld] = MAX(ABS(signal_FFT_x_adiabatic_filtered))
        (*parameters.Sx_adiabatic_array)[parameters.harmonics_n_max+1, 0, md, ld] = MAX(ABS(noisy_signal_x_adiabatic)); Unfiltered signal

      ENDFOR
    ENDFOR
    ;@GJ, 2024/8/12, display the 3rd phase image
    ;iimage, REFORM((*parameters.Sx_array)[2, 1, *, *]), title='3rd Phase Image'
    ;destroy the progress bar
    progressbar -> Destroy
  ENDIF
  
  ;@GJ, 2024/2/11, interpolate the NaN elements
  ;  iimage,Sx_array_xzNonCOnv, title='base'
  FOR n=0, parameters.harmonics_n_max+1 DO BEGIN
    FOR m=0, parameters.harmonics_comp_N-1 DO BEGIN
      temp_Sx_array = REFORM((*parameters.Sx_array)[n, m, *, *])
      index_nan = where(~FINITE(temp_Sx_array), count)
      IF count GT 0 THEN BEGIN
        temp_Sx_array[index_nan] = 0.
;        FOR j=0, count-1 DO BEGIN
;          ind_2d = ARRAY_INDICES(temp_Sx_array, index_nan[j])
;          temp_Sx_array[index_nan[j]] = (temp_Sx_array[ind_2d[0]-1, ind_2d[1]] + temp_Sx_array[ind_2d[0]+1, ind_2d[1]])/2.
;        ENDFOR
        (*parameters.Sx_array)[n, m, *, *] = temp_Sx_array
      ENDIF
    ENDFOR
  ENDFOR
  
  ;@GJ, 2024/3/27, plotting the curves
  harmonics_n = FLOOR(parameters.harmonics_n)
  Sx_array = REFORM((*parameters.Sx_array)[harmonics_n-1, *, *, *])
  ind_z = WHERE(ABS(z_array) EQ 0)
  ind_x = WHERE(ABS(x_array) EQ 0)
  IF parameters.particle_type LE 2 THEN colorRBG = [255,235,177]
  IF parameters.particle_type EQ 3 THEN colorRBG = [255, 0, 0]
  IF parameters.particle_type EQ 4 THEN colorRBG = [0, 176, 240]
;  iplot, x_array[ind_x:*]/parameters.A, REFORM(Sx_array[0,ind_z,ind_x:*]), LINESTYLE = 0, /FILL_BACKGROUND, FILL_TRANSPARENC=50, FILL_COLOR=colorRBG;[255,235,177];
;  iplot, x_array[ind_x:*]/parameters.A, REFORM(Sx_array[1,ind_z,ind_x:*])/MAX(ABS(Sx_array[1,ind_z,ind_x:*]))*MAX(ABS(Sx_array[0,ind_z,ind_x:*])), LINESTYLE = 1, /overplot
;  iplot, x_array[ind_x:*]/parameters.A, REFORM(Sx_array[2,ind_z,ind_x:*]), LINESTYLE = 2, /overplot
;  iplot, x_array[ind_x:*]/parameters.A, REFORM(Sx_array[3,ind_z,ind_x:*]), LINESTYLE = 3, /overplot
;  
;@GJ, 2024,/4/18, plot the data
;  iplot, x_array/parameters.A, REFORM(Sx_array[0,ind_z,*]), LINESTYLE = 0, XRANGE=[MIN(x_array/parameters.A), MAX(x_array/parameters.A)], YRANGE=[-MAX(Sx_array[0,*,*])*1.1, MAX(Sx_array[0,*,*])*1.1], XTITLE='HDC/HAC', YTITLE='A3', /FILL_BACKGROUND, FILL_TRANSPARENC=50, FILL_COLOR=colorRBG;[255,235,177];
;  iplot, x_array/parameters.A, REFORM(Sx_array[1,ind_z,*])/MAX(ABS(Sx_array[1,ind_z,*]))*MAX(ABS(Sx_array[0,ind_z,*])), LINESTYLE = 1, /overplot
;  iplot, x_array/parameters.A, REFORM(Sx_array[2,ind_z,*]), LINESTYLE = 2, /overplot
;  iplot, x_array/parameters.A, REFORM(Sx_array[3,ind_z,*]), LINESTYLE = 3, /overplot
;  
;  max_z = MAX(Sx_array[0,*,214.], max_ind_z)
;  iplot, x_array/parameters.A, REFORM(Sx_array[0,max_ind_z,*]), LINESTYLE = 0, XRANGE=[MIN(x_array/parameters.A), MAX(x_array/parameters.A)], YRANGE=[-MAX(Sx_array[0,*,*])*1.1, MAX(Sx_array[0,*,*])*1.1], XTITLE='HDC/HAC', YTITLE='A3', /FILL_BACKGROUND, FILL_TRANSPARENC=50, FILL_COLOR=colorRBG+[0,0,100];[255,235,177];
;  iplot, x_array/parameters.A, REFORM(Sx_array[1,max_ind_z,*])/MAX(ABS(Sx_array[1,ind_z,*]))*MAX(ABS(Sx_array[0,ind_z,*])), LINESTYLE = 1, /overplot
;  iplot, x_array/parameters.A, REFORM(Sx_array[2,max_ind_z,*]), LINESTYLE = 2, /overplot
;  iplot, x_array/parameters.A, REFORM(Sx_array[3,max_ind_z,*]), LINESTYLE = 3, /overplot
;  
;  iplot, z_array, REFORM(Sx_array[0,*,224]), title='PSF', xrange=[-25, 25], XTITLE='HG [mT]', YTITLE='A3', /FILL_BACKGROUND, FILL_TRANSPARENC=50, FILL_COLOR=[255,0,0];
;  iplot, z_array, REFORM(Sx_array[1,*,224])/MAX(ABS(Sx_array[1,*,224]))*MAX(ABS(Sx_array[0,*,224])), LINESTYLE = 1, /overplot
;  iplot, z_array, REFORM(Sx_array[2,*,224]), LINESTYLE = 2, /overplot
;  iplot, z_array, REFORM(Sx_array[3,*,224]), LINESTYLE = 3, /overplot

;@GJ, 2024,/4/18, end plot the data  
;
;  
;  iplot, z_array/parameters.A, REFORM(Sx_array[0,*,ind_x]), LINESTYLE = 0, /FILL_BACKGROUND, FILL_TRANSPARENC=50, FILL_COLOR=colorRBG;[255,235,177];
;  iplot, z_array/parameters.A, REFORM(Sx_array[1,*,ind_x])/MAX(ABS(Sx_array[1,*,ind_x]))*MAX(ABS(Sx_array[0,*,ind_x])), LINESTYLE = 1, /overplot
;  iplot, z_array/parameters.A, REFORM(Sx_array[2,*,ind_x]), LINESTYLE = 2, /overplot
;  iplot, z_array/parameters.A, REFORM(Sx_array[3,*,ind_x]), LINESTYLE = 3, /overplot
  
END

;@GJ, 2023/12/10, including the orghogonal Gradient field
pro sin_simulation_Chebyshev,wid,parameters
  
  print, 'wid:', wid
  tau = brownian_time_calc(parameters)
  parameters.tau = tau
  widget_tau = widget_info(wid, find_by_uname='tau')
  widget_control, widget_tau, set_value = tau
  
  IF parameters.AC_DC_mode EQ 0. THEN BEGIN
    ;AC excitation mode
    widget_label_A_DC = widget_info(wid, find_by_uname='label_A_DC')
    widget_control, widget_label_A_DC, set_value = "AC (Red, %AC)"
        
    widget_label_A_DC_1 = widget_info(wid, find_by_uname='label_A_DC_1')
    widget_control, widget_label_A_DC_1, set_value = "AC (Green, %AC)"
    
    widget_A_DC = widget_info(wid, find_by_uname='A_DC')
    widget_control, widget_A_DC, set_value = parameters.AC_ex*100., set_slider_min=parameters.AC_ex_min*100., set_slider_max=parameters.AC_ex_max*100. 

    widget_A_DC_1 = widget_info(wid, find_by_uname='A_DC_1')
    widget_control, widget_A_DC_1, set_value = parameters.AC_ex_1*100., set_slider_min=parameters.AC_ex_min*100., set_slider_max=parameters.AC_ex_max*100.
  ENDIF ELSE BEGIN
    ;DC offset mode
    widget_label_A_DC = widget_info(wid, find_by_uname='label_A_DC')
    widget_control, widget_label_A_DC, set_value = "DC Offset (Red, %AC)"

    widget_label_A_DC_1 = widget_info(wid, find_by_uname='label_A_DC_1')
    widget_control, widget_label_A_DC_1, set_value = "DC Offset(Green, %AC)"
   
    ;@GJ, 2024/2/12, update the widget A_DC and A_DC_1
    widget_A_DC = widget_info(wid, find_by_uname='A_DC')
    widget_control, widget_A_DC, set_value = parameters.A_DC, set_slider_min=-(parameters.A_DC_max-1), set_slider_max=parameters.A_DC_max-1
    
    widget_A_DC_1 = widget_info(wid, find_by_uname='A_DC_1')
    widget_control, widget_A_DC_1, set_value = parameters.A_DC_1, set_slider_min=-(parameters.A_DC_max-1), set_slider_max=parameters.A_DC_max-1
  ENDELSE
  
  tau_ms = parameters.tau/1000.;@GJ, 2023/12/24, tau in us is converted to ms
  A = parameters.A
  A_DC = parameters.A_DC/100.*A; DC selection
  A_DC_1 = parameters.A_DC_1/100.*A; DC selection
  gradient_field_x = parameters.A * 2. * (parameters.A_DC_max/100.) / parameters.FOV ;calculating the gradient field
  ;gradient_field_x = parameters.gradient_field_x * 0.1 ;mT/mm
  gradient_field_yz = parameters.gradient_field_yz * 0.1 ;mT/mm
  gradient_field_yz_current = parameters.gradient_field_yz_current * 0.1 ; mT
  
  ;@GJ, 2023/12/26, packing to a function for signal calculation
  n_x = parameters.n_x
  n_z = parameters.n_z
  FOV = parameters.FOV
  delta_x = FOV/n_x
  delta_z = FOV/n_z
  IF parameters.recalc_flag EQ 1. THEN BEGIN
    Signal_calculation_FM_MPI, parameters
    ;@GJ, 2024/9/9, recalculating the original angle
    parameters.shift_angle_renew = 0.
  ENDIF
;  parameters.rec_directory = parameters.rec_directory+'rec'+STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(f0)')+'\'
  
  rows_S = n_x
  cols_S = n_z
  IF parameters.AC_DC_mode EQ 0. THEN BEGIN
    ;AC excitation mode
    A_DC_array = FINDGEN(n_x)/n_x * (parameters.AC_ex_max - parameters.AC_ex_min) + parameters.AC_ex_min
    min_A_DC = MIN(ABS(A_DC_array-parameters.AC_ex), index_A_DC)
    index_A_DC += 1.
    min_A_DC_1 = MIN(ABS(A_DC_array-parameters.AC_ex_1), index_A_DC_1)
    index_A_DC_1 += 1.
  ENDIF ELSE BEGIN
    ;DC offset mode
    A_DC_array = (FINDGEN(n_x) - n_x/2.) * delta_x * gradient_field_x
    index_A_DC = FLOOR(A_DC / (gradient_field_x * delta_x)) + n_x / 2.
    index_A_DC_1 = FLOOR(A_DC_1 / (gradient_field_x * delta_x)) + n_x / 2.
  ENDELSE

  harmonics_n = FLOOR(parameters.harmonics_n)
  harmonics_comp = parameters.harmonics_comp
  harmonics_n_max = parameters.harmonics_n_max
  ;@GJ, 2024/2/11, selecting the filtered, unfiltered or harmonic signal
  IF parameters.signal_type EQ 0 THEN BEGIN
    curr_Sx_array = REFORM((*parameters.Sx_array)[harmonics_n_max, 0, *, *])
    curr_Sx_array_real = REFORM((*parameters.Sx_array)[harmonics_n_max, 2, *, *])
    curr_Sx_array_imag = REFORM((*parameters.Sx_array)[harmonics_n_max, 3, *, *])
    
    curr_Sx_adiabatic_array = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n_max, 0, *, *])
    curr_Sx_adiabatic_array_real = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n_max, 2, *, *])
    curr_Sx_adiabatic_array_imag = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n_max, 3, *, *])
  ENDIF
  IF parameters.signal_type EQ 1 THEN BEGIN
    curr_Sx_array = REFORM((*parameters.Sx_array)[harmonics_n_max+1, 0, *, *])
    curr_Sx_array_real = REFORM((*parameters.Sx_array)[harmonics_n_max+1, 2, *, *])
    curr_Sx_array_imag = REFORM((*parameters.Sx_array)[harmonics_n_max+1, 3, *, *])
    
    curr_Sx_adiabatic_array = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n_max+1, 0, *, *])
    curr_Sx_adiabatic_array_real = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n_max+1, 2, *, *])
    curr_Sx_adiabatic_array_imag = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n_max+1, 3, *, *])
  ENDIF
  IF parameters.signal_type EQ 2 THEN BEGIN
    curr_Sx_array = REFORM((*parameters.Sx_array)[harmonics_n-1, harmonics_comp, *, *])
    curr_Sx_array_real = REFORM((*parameters.Sx_array)[harmonics_n-1, 2, *, *])
    curr_Sx_array_imag = REFORM((*parameters.Sx_array)[harmonics_n-1, 3, *, *])
    
    curr_Sx_adiabatic_array = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n-1, harmonics_comp, *, *])
    curr_Sx_adiabatic_array_real = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n-1, 2, *, *])
    curr_Sx_adiabatic_array_imag = REFORM((*parameters.Sx_adiabatic_array)[harmonics_n-1, 3, *, *])
  ENDIF
  
  focal_coor = FINDGEN(n_z) - n_z/2.
  ;@GJ, 2023/12/31, calculate the radius and thick
  ;calculate the donut diameter and thickness
  A_G_Array = focal_coor * delta_z * gradient_field_yz
  *parameters.A_G_Array = A_G_Array
  
  ;@GJ, 2024/9/8, plot the LEO (line of excitation offset)
  i_DC=0
  DC_line_real = REFORM(curr_Sx_array_real[*, n_x/2 - i_DC])
  DC_line_imag = REFORM(curr_Sx_array_imag[*, n_x/2 - i_DC])
  DC_line_adiabatic_real = REFORM(curr_Sx_adiabatic_array_real[*, n_x/2 - i_DC])
  DC_line_adiabatic_imag = REFORM(curr_Sx_adiabatic_array_imag[*, n_x/2 - i_DC])
  
  nwtao = harmonics_n * (parameters.f * 2. * !PI) * REFORM((*parameters.tau_ms_array)[*, n_x/2 - i_DC])
  debye_f3_real = 0.5 * parameters.f / (1. + nwtao^2)
  debye_f3_imag = 0.5 * parameters.f * nwtao / (1. + nwtao^2)
  DC_line_real_equ = -DC_line_adiabatic_real * debye_f3_real + DC_line_adiabatic_imag * debye_f3_imag
  DC_line_imag_equ = DC_line_adiabatic_real * debye_f3_imag + DC_line_adiabatic_imag * debye_f3_real
  ;@GJ, 2024/9/9, do the plotting
  IF parameters.phase_scatter_plot EQ 1 THEN BEGIN
    max_real = MAX(ABS([DC_line_adiabatic_real, DC_line_real]))
    max_imag = MAX(ABS([DC_line_adiabatic_imag, DC_line_real]))
    max_abs = SQRT(MAX_real^2 + max_imag^2)
    min_real = MIN([DC_line_adiabatic_real, DC_line_real])
    min_imag = MIN([DC_line_adiabatic_imag, DC_line_real])
    iplot, DC_line_real/max_abs, DC_line_imag/max_abs, xtitle='Real', ytitle='Imag', title='LEO w/o Adiabatic', xrange=[-0.2, 1.1], yrange=[-0.2, 1.1], LINESTYLE=0, THICK=2, color='red', /NO_SAVEPROMPT
    iplot, DC_line_adiabatic_real/max_abs, DC_line_adiabatic_imag/max_abs, LINESTYLE=0, THICK=2, color='blue', /overplot
    ;@GJ, 2024/9/13, fix the bug of display negative real/negative imag curves
    iplot, DC_line_real_equ/MAX(ABS(DC_line_real_equ))*MAX(ABS(DC_line_real))/max_abs, DC_line_imag_equ/MAX(ABS(DC_line_imag_equ))*MAX(ABS(DC_line_imag))/max_abs, LINESTYLE=3, THICK=2, color='green', /overplot

;    DC_real_imag = DC_line_real / DC_line_imag
;    result1 = LINFIT(nwtao, DC_real_imag, yfit=yfit1)
;    iplot, nwtao, DC_real_imag, SYM_INDEX=2, LINESTYLE=0, xtitle='nwtao', ytitle='Real/Imaginary', title='FFT of Debye Relaxation Real', color='red', /NO_SAVEPROMPT
;    iplot, nwtao, yfit1, color='blue', /overplot
;    print, 'slope Real: ', Result1[1]
;    DC_line_imag_adia = DC_line_imag / DC_line_adiabatic_imag
;    result2 = LINFIT(wnwtao, DC_line_imag_adia, yfit=yfit2)
;    iplot, wnwtao, DC_line_imag_adia, SYM_INDEX=2, LINESTYLE=0, xtitle='wnwtao', ytitle='Imaginary/Adia', title='FFT of Debye Relaxation Imag', color='red', /NO_SAVEPROMPT
;    iplot, wnwtao, yfit2, color='blue', /overplot
;    print, 'slope Imag: ', Result2[1]
    
    ;@GJ, 2024/9/11, check the relaxation time map
    iimage, ROTATE((*parameters.tau_ms_array)*1000., 1), title='Debye Relaxation Map [us]'
    iplot, A_G_array, REFORM((*parameters.tau_ms_array)[*, n_x/2 - i_DC]), xtitle='HG [mT]', ytitle='Relaxation Time [ms]', title='Relaxation Time', /NO_SAVEPROMPT
  ENDIF
  
  ;@GJ, 2024/9/8, calculate the shift angle
  IF parameters.shift_angle_renew EQ 0 THEN BEGIN
    ;@GJ, 2024/9/12, if the type is for best SNR
    IF parameters.shift_angle_type EQ 0 THEN BEGIN
      current_angle = 180/!PI*atan(curr_Sx_array_imag[n_x/2, n_x/2], curr_Sx_array_real[n_x/2, n_x/2])
      parameters.shift_angle = 90. - current_angle
    ENDIF
    ;@GJ, 2024/9/12, if the type is for best resolution
    IF parameters.shift_angle_type EQ 1 THEN BEGIN
      current_angle = 180/!PI*atan(MEAN(DC_line_imag), MEAN(DC_line_real))
      parameters.shift_angle = 180 - current_angle
    ENDIF
    widget_text_shift_angle = widget_info(wid, find_by_uname='shift_angle')
    widget_control, widget_text_shift_angle, set_value = STRING(parameters.shift_angle, format='(i4)')
  ENDIF

  IF parameters.shift_angle NE 0 THEN BEGIN
    new_curr_Sx_array_real = -curr_Sx_array_imag * sin(parameters.shift_angle/180.*!PI) + curr_Sx_array_real * cos(parameters.shift_angle/180.*!PI)
    new_curr_Sx_array_imag = curr_Sx_array_imag * cos(parameters.shift_angle/180.*!PI) + curr_Sx_array_real * sin(parameters.shift_angle/180.*!PI)
    curr_Sx_array_real = new_curr_Sx_array_real
    curr_Sx_array_imag = new_curr_Sx_array_imag
    IF harmonics_comp EQ 1 THEN curr_Sx_array += parameters.shift_angle
    IF harmonics_comp EQ 2 THEN curr_Sx_array = new_curr_Sx_array_real
    IF harmonics_comp EQ 3 THEN curr_Sx_array = new_curr_Sx_array_imag
    
    new_curr_Sx_adiabatic_array_real = -curr_Sx_adiabatic_array_imag * sin(parameters.shift_angle/180.*!PI) + curr_Sx_adiabatic_array_real * cos(parameters.shift_angle/180.*!PI)
    new_curr_Sx_adiabatic_array_imag = curr_Sx_adiabatic_array_imag * cos(parameters.shift_angle/180.*!PI) + curr_Sx_adiabatic_array_real * sin(parameters.shift_angle/180.*!PI)
    curr_Sx_adiabatic_array_real = new_curr_Sx_adiabatic_array_real
    curr_Sx_adiabatic_array_imag = new_curr_Sx_adiabatic_array_imag
    IF harmonics_comp EQ 1 THEN curr_Sx_adiabatic_array += parameters.shift_angle
    IF harmonics_comp EQ 2 THEN curr_Sx_adiabatic_array = new_curr_Sx_adiabatic_array_real
    IF harmonics_comp EQ 3 THEN curr_Sx_adiabatic_array = new_curr_Sx_adiabatic_array_imag
  ENDIF
  
  ;@GJ, 2024/6/4, find the signal rato for noise DB
  signal_ratio = REFORM((*parameters.Sx_array)[*, 0, 0, 0])
  FOR i=0, N_ELEMENTS(signal_ratio)-1 DO signal_ratio[i] = MEAN((*parameters.Sx_array)[i, 0, *, *])
  signal_ratio /= signal_ratio[0]
  current_signal_ratio = signal_ratio[harmonics_n-1]
;  iplot, signal_ratio[0:harmonics_n_max-1], title='Signal ratio', xtitle='Harmonic #'
  
  ;@GJ, 2024/5/5, removing the middle dark line
  IF parameters.signal_type LE 2 AND MEAN(curr_Sx_array[n_z/2, *]) LT MEAN(curr_Sx_array[n_z/2-1, *]) THEN BEGIN
    curr_Sx_array[n_z/2, *] =(curr_Sx_array[n_z/2-1, *] + curr_Sx_array[n_z/2+1, *]) / 2.
    FOR i_ind=0, n_x-1 DO curr_Sx_array[n_z/2, i_ind] = INTERPOL(REFORM(curr_Sx_array[*, i_ind]), FINDGEN(n_z), n_z/2, /LSQUADRATIC)
  ENDIF
  ;@GJ, 2024/9/9, adiabatic part
  IF parameters.signal_type LE 2 AND MEAN(curr_Sx_adiabatic_array[n_z/2, *]) LT MEAN(curr_Sx_adiabatic_array[n_z/2-1, *]) THEN BEGIN
    curr_Sx_adiabatic_array[n_z/2, *] =(curr_Sx_adiabatic_array[n_z/2-1, *] + curr_Sx_adiabatic_array[n_z/2+1, *]) / 2.
    FOR i_ind=0, n_x-1 DO curr_Sx_adiabatic_array[n_z/2, i_ind] = INTERPOL(REFORM(curr_Sx_adiabatic_array[*, i_ind]), FINDGEN(n_z), n_z/2, /LSQUADRATIC)
  ENDIF
  ;@GJ, 2024/3/26, if there is no signal, jump
  IF STDEV(curr_Sx_array) LT 0.0001 AND MEAN(curr_Sx_array) LT 0.0001 THEN BEGIN
    warning = DIALOG_MESSAGE('No values for display!', /ERROR)
    GOTO, plot_End_spot
  ENDIF
  
;  iimage, curr_Sx_array, title='harmonic map'
  
  focal_line = REFORM(curr_Sx_array[*, n_x - index_A_DC]) ;@GJ, 2024/8/12, fix the bug of selecting focal line
  *parameters.focal_line = focal_line
  focal_line_1 = REFORM(curr_Sx_array[*, n_x - index_A_DC_1])  ;@GJ, 2024/8/12, fix the bug of selecting focal line
  *parameters.focal_line_1 = focal_line_1
  
  focal_line_real = REFORM(curr_Sx_array_real[*, n_x - index_A_DC]) ;@GJ, 2024/8/12, fix the bug of selecting focal line
  *parameters.focal_line_real = focal_line_real
  focal_line_1_real = REFORM(curr_Sx_array_real[*, n_x - index_A_DC_1])  ;@GJ, 2024/8/12, fix the bug of selecting focal line
  *parameters.focal_line_1_real = focal_line_1_real
  
  focal_line_imag = REFORM(curr_Sx_array_imag[*, n_x - index_A_DC]) ;@GJ, 2024/8/12, fix the bug of selecting focal line
  *parameters.focal_line_imag = focal_line_imag
  focal_line_1_imag = REFORM(curr_Sx_array_imag[*, n_x - index_A_DC_1])  ;@GJ, 2024/8/12, fix the bug of selecting focal line
  *parameters.focal_line_1_imag = focal_line_1_imag
  
  IF parameters.phase_scatter_plot EQ 1 THEN BEGIN
    mean_imag = DBLARR(n_x/2-1)
    mean_real = DBLARR(n_x/2-1)
    linear_ratio = DBLARR(n_x/2-1)

    ;@GJ, 2024/11/8, plot the MRF fingerprinting
    iimage, REFORM(curr_Sx_array[*, *]), title='Signal (MRF)'
    FOR m=0, 98 DO BEGIN
      ind_m = 2.-m/100.
      IF m EQ 0 THEN BEGIN
        iplot, REFORM(curr_Sx_array_real[n_z/ind_m, *]), REFORM(curr_Sx_array_imag[n_z/ind_m, *]), color='red', xtitle='R3', ytitle='I3', title='Magnetic Particle Fingerprinting', /NO_SAVEPROMPT
      ENDIF ELSE BEGIN
        iplot, REFORM(curr_Sx_array_real[n_z/ind_m, *]), REFORM(curr_Sx_array_imag[n_z/ind_m, *]), color=350000*(m+1), /overplot
      ENDELSE
    ENDFOR
    
    ;@GJ, 2024/12/10, plot the PSF and calculate the PSF
    FWHM_g = DBLARR(n_x) * 0.
    FWHM_STED =  DBLARR(n_x) * 0.
    filter_peak = DBLARR(n_x) * 0.
    FOR i_DC=0, n_x-1 DO BEGIN
      DC_line = REFORM(curr_Sx_array[*, i_DC])
      DC_line_real = REFORM(curr_Sx_array_real[*, i_DC])
      DC_line_imag = REFORM(curr_Sx_array_imag[*, i_DC])
      
      ;@GJ, 2024/12/10, do STED imaging
      real_1d = DC_line_real
      imag_1d = DC_line_imag
      result = POLY_FIT(real_1d, imag_1d, 6, yfit=imag_1d_fit)
      ;rotate to get the donut-shaped focal spot
      max_real = MAX(ABS(real_1d), maxid)
      angle = 180./!PI*ATAN(imag_1d_fit[maxid], real_1d[maxid], /phase)
      ;@GJ, 2024/8/23, plot the rotated
      phi_g = angle / 180. * !PI - 0.5 * !PI
      real_1d_g = imag_1d_fit * sin(phi_g) + real_1d * cos(phi_g)
      imag_1d_g = -real_1d * sin(phi_g) + imag_1d_fit * cos(phi_g)
      max_real_g = MAX(ABS(real_1d_g), max_real_id_g)
      ;  max_imag_g = MAX(ABS(imag_1d_g), max_imag_id_g)
      phi_d = angle / 180. * !PI
      real_1d_d = imag_1d_fit * sin(phi_d) + real_1d * cos(phi_d)
      imag_1d_d = -real_1d * sin(phi_d) + imag_1d_fit * cos(phi_d)
      ;  max_real_d = MAX(ABS(real_1d_d), max_real_id_d)
      max_imag_d = MAX(ABS(imag_1d_d), max_imag_id_d)
      ;calculate dog factor
      dog_factor = ABS(imag_1d_g[max_real_id_g]) / max_imag_d

      ;calculate the FWHM
      yfit_psf_g = GAUSSFIT(A_G_array, imag_1d_g, coeff_g, NTERMS=5,  ESTIMATES= coeff_g)
      print, 'Gaussian FWHM: ', 2*SQRT(2*ALOG(2))*coeff_g[2], ' mT, ', 2*SQRT(2*ALOG(2))*coeff_g[2] / gradient_field_yz, ' mm'
      yfit_psf_STED = GAUSSFIT(A_G_array, imag_1d_g - imag_1d_d*dog_factor, coeff_STED, NTERMS=4)
      ;  print, 'STED PSF Fit Result: ', coeff_STED[0:4-1]
      print, 'STED FWHM: ', 2*SQRT(2*ALOG(2))*coeff_STED[2], ' mT, ', 2*SQRT(2*ALOG(2))*coeff_STED[2] / gradient_field_yz, ' mm'
      yfit_psf_d = yfit_psf_g - yfit_psf_STED
      
      IF i_DC EQ 0 OR i_DC EQ n_x-2 THEN BEGIN
        iplot, real_1d, imag_1d, xtitle='Real', ytitle='Imag', color='purple', SYM_INDEX=1, LINESTYLE=6, title='PSF Scatter Plot', /NO_SAVEPROMPT
        iplot, real_1d, imag_1d_fit, COLOR='purple', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
        iplot, real_1d_g, imag_1d_g, COLOR='green', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
        iplot, real_1d_d, imag_1d_d, COLOR='red', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
        iplot, real_1d_d, imag_1d_d*dog_factor, COLOR='red', SYM_INDEX=0, THICK=2, LINESTYLE=2, /OVERPLOT
        iplot, real_1d_d, imag_1d_g - imag_1d_d*dog_factor, COLOR='blue', SYM_INDEX=0, THICK=2, LINESTYLE=1, /OVERPLOT
        iplot, A_G_array, imag_1d_g - imag_1d_d*dog_factor, SYM_INDEX=1, LINESTYLE=6, xtitle='Offset Field [mT]', ytitle='Imag', color='blue', title='STED-MPI', /NO_SAVEPROMPT
        iplot, A_G_array, yfit_psf_sted, LINSTYLE=0, thick=2, color='blue', /OVERPLOT
        iplot, A_G_array, imag_1d_g, SYM_INDEX=1, LINESTYLE=6, color='green', /OVERPLOT
        iplot, A_G_array, yfit_psf_g, LINSTYLE=0, thick=2, color='green', /OVERPLOT
        iplot, A_G_array, imag_1d_d*dog_factor, SYM_INDEX=1, LINESTYLE=6, color='red', /OVERPLOT
        iplot, A_G_array, yfit_psf_d, LINSTYLE=0, thick=2, color='red', /OVERPLOT
        
        ;@GJ, 2024/12/18, define the filtration in FBP for STED-MPI
        ; Frequency for FFT in 1/mT, unit
        ; N is an integer giving the number of elements in a particular dimension
        ; T is a floating-point number giving the sampling interval
        g_ele = N_Elements(A_G_array)
        delta_g = (A_G_array[1]-A_G_array[0]) ; us; samplying rate: 1MHz
        X = (FINDGEN((g_ele - 1)/2) + 1)
        is_N_even = (g_ele MOD 2) EQ 0
        if (is_N_even) then $
          frequency_g = [0.0, X, g_ele/2, -g_ele/2 + X]/(g_ele*delta_g) $
        else $
          frequency_g = [0.0, X, -(g_ele/2 + 1) + X]/(g_ele*delta_t)
        FOR dog_index = 0, 20 DO BEGIN
          temp_dog_factor = dog_index / 10. * dog_factor
          IF dog_index EQ 0 THEN BEGIN
            iplot, frequency_g, ABS(FFT(imag_1d_g - imag_1d_d*temp_dog_factor)), SYM_INDEX=0, LINESTYLE=0, xrange=[-2., 2], xtitle='Freq [1/mT]', ytitle='Amp', color='blue', title='Filtration Spec', /NO_SAVEPROMPT
          ENDIF ELSE BEGIN
            iplot, frequency_g, ABS(FFT(imag_1d_g - imag_1d_d*temp_dog_factor)), SYM_INDEX=0, LINESTYLE=0, color=350000*(dog_index+1), /OVERPLOT
          ENDELSE
        ENDFOR
        print, 'filtration in FBP'
        

       
      ENDIF
      
      max_filter = MAX(ABS(FFT(imag_1d_g - imag_1d_d*dog_factor)), max_ind)
      filter_peak[i_DC] = ABS(frequency_g[max_ind])
      FWHM_g[i_DC] = 2*SQRT(2*ALOG(2))*coeff_g[2]
      FWHM_STED[i_DC] =  2*SQRT(2*ALOG(2))*coeff_STED[2]
    ENDFOR
    iplot, A*A_DC_array, FWHM_g, xtitle='Magnetic Field [mT]', ytitle='FWHM [mT]', color='green', SYM_INDEX=1, LINESTYLE=6, title='FWHM vs Field Amplitude', /NO_SAVEPROMPT
    iplot, A*A_DC_array, FWHM_STED, color='blue', SYM_INDEX=1, LINESTYLE=6, /overplot
    iplot, A*A_DC_array, filter_peak, xtitle='Magnetic Field [mT]', ytitle='Peak Freq [1/mT]', color='red', SYM_INDEX=1, LINESTYLE=6, title='Peak Freq vs Field Amplitude', /NO_SAVEPROMPT

    ;@GJ, 2024/8/31, plot the excitation string
    FOR i_DC=0, n_x/2-2 DO BEGIN
      DC_line = REFORM(curr_Sx_array[*, n_x/2 - i_DC])
      DC_line_real = REFORM(curr_Sx_array_real[*, n_x/2 - i_DC])
      DC_line_imag = REFORM(curr_Sx_array_imag[*, n_x/2 - i_DC])
      mean_real[i_DC] = MAX(DC_line_real)
      mean_imag[i_DC] = MEAN(DC_line_imag)
      fit_result = LINFIT(DC_line_real, DC_line_imag, CHISQR=fit_chisqr)
      linear_ratio[i_DC] = fit_chisqr
   
;      wait, 1
      IF i_DC EQ 0 THEN BEGIN
        max_real = MAX(ABS(DC_line_real))
        max_imag = MAX(ABS(DC_line_imag))
        max_r = MAX([max_real, max_imag])
        iplot, DC_line_real, DC_line_imag, xtitle='Real', ytitle='Imag', title='Scatter Plot', color='red', /NO_SAVEPROMPT ;xrange=[-max_r, max_r], yrange=[-max_r, max_r], SYM_INDEX=3, LINESTYLE=0, 
      ENDIF ELSE BEGIN
        max_real = MAX([max_real, MAX(ABS(DC_line_real))])
        max_imag = MAX([max_imag, MAX(ABS(DC_line_imag))])
        max_r = MAX([max_real, max_imag])
        iplot, DC_line_real, DC_line_imag, color=350000*(i_DC+1), /overplot;xrange=[-max_r, max_r], yrange=[-max_r, max_r], SYM_INDEX=2, LINESTYLE=0, 
      ENDELSE
    ENDFOR
    ;@GJ, 2024/8/31, plot the phase under different DC offset
;    iplot, mean_real, mean_imag, xtitle='Max Real', ytitle='Imag', yerror=linear_ratio/MAX(linear_ratio)*30., title='Mean Plot vs DC Offset', /NO_SAVEPROMPT
    
;    ;@GJ, 2024/9/8, plot the LEO (line of excitation offset)
;    i_DC=0
;    DC_line_real = REFORM(curr_Sx_array_real[*, n_x/2 - i_DC])
;    DC_line_imag = REFORM(curr_Sx_array_imag[*, n_x/2 - i_DC])
;    nwtao = harmonics_n * (parameters.f * 2. * !PI) * REFORM((*parameters.tau_ms_array)[*, n_x/2 - i_DC])
;    DC_real_imag = DC_line_real / DC_line_imag
;    result = LINFIT(nwtao, DC_real_imag, yfit=yfit)
;    iplot, nwtao, DC_line_real / DC_line_imag, SYM_INDEX=2, LINESTYLE=0, xtitle='nwtao', ytitle='Real/Imaginary', title='FFT of Debye Relaxation After Rotation', color='red', /NO_SAVEPROMPT
;    iplot, nwtao, yfit, color='blue', /overplot
;    print, 'slope after rotation: ', Result[1]
;    iplot, A_G_array,(DC_line_real / DC_line_imag)/(harmonics_n * parameters.f * REFORM((*parameters.tau_ms_array)[*, n_x/2 - i_DC])), xrange=[-MAX(A_G_array)*1.1, MAX(A_G_array)*1.1], yrange=[0, MAX([harmonics_n+2, 5])], xtitle='HG [mT]', ytitle='Delay Ratio', title='LEO Slope', /NO_SAVEPROMPT
  ENDIF
 
;  IF parameters.shift_angle NE 0 THEN BEGIN
;    focal_line_real = -focal_line_imag * sin(parameters.shift_angle/180.*!PI) + focal_line_real * cos(parameters.shift_angle/180.*!PI)
;    focal_line_imag = focal_line_imag * cos(parameters.shift_angle/180.*!PI) + focal_line_real * sin(parameters.shift_angle/180.*!PI)
;    focal_line_1_real = -focal_line_1_imag * sin(parameters.shift_angle/180.*!PI) + focal_line_1_real * cos(parameters.shift_angle/180.*!PI)
;    focal_line_1_imag = focal_line_1_imag * cos(parameters.shift_angle/180.*!PI) + focal_line_1_real * sin(parameters.shift_angle/180.*!PI)
;  ENDIF
  

  
  ;@GJ, 2024/5/16, plot the FFT of focal line, calculate
  focal_line_FFT = ABS(FFT(focal_line))
  focal_line_1_FFT = ABS(FFT(focal_line_1))
  focal_line_sub_FFT = ABS(FFT(focal_line - focal_line_1 * parameters.dog_factor/10.))
  ;@GJ, 2022/12/14, redefining the frequency
  ; Frequency for FFT in kHz, unit
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  t_ele = N_Elements(focal_line)
  f=1.
  delta_t = (focal_coor[1]-focal_coor[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
  X = (FINDGEN((t_ele - 1)/2) + 1)
  is_N_even = (t_ele MOD 2) EQ 0
  if (is_N_even) then $
    frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
  else $
    frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
;  iplot, frequency, focal_line_FFT, xrange=[0, 20./(t_ele*delta_t)], color='red', title='FFT of focal spots'
;  iplot, frequency, focal_line_1_FFT, color='blue', /overplot
;  iplot, frequency, focal_line_sub_FFT, color='green', /overplot

  
;  ;@GJ, 2024/4/18, optimizing dog factor
;  Dog_factor_optimization, wid, parameters

  zero_ind=WHERE(ABS(A_G_Array) EQ 0)
  donut_Harmonic_0 = REFORM(focal_line[zero_ind:*])
  temp_A_G_array = A_G_array[zero_ind:*]
  maxD_Harmonic_0 = MAX(donut_Harmonic_0, maxD_ind_Harmonic_0)
  IF maxD_ind_Harmonic_0 GT 3. THEN BEGIN
    donut_flag_0 = 1
    donut_radius_0 = temp_A_G_array[maxD_ind_Harmonic_0]
    donut_thick_0 = donut_Harmonic_0[maxD_ind_Harmonic_0] - donut_Harmonic_0[0]
  ENDIF ELSE BEGIN
    ;non-donut, calculate FWHM
    FWHM_ind_0 = INTERPOL(temp_A_G_array, donut_Harmonic_0, 0.5*maxD_Harmonic_0)
    donut_flag_0 = 0
    donut_radius_0 = 2. * FWHM_ind_0
    donut_thick_0 = maxD_Harmonic_0
  ENDELSE
  
  h_FWHM_0 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_MTF_0 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_radius_0 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_energy_0 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_FWHM_1 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_MTF_1 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_radius_1 = DBLARR(parameters.harmonics_n_max+2)*0.
  h_energy_1 = DBLARR(parameters.harmonics_n_max+2)*0.
  ;@GJ, 2024/5/21
  flag_i = 0;
  FOR h_i=0, parameters.harmonics_n_max+1 DO BEGIN
    IF h_i EQ 0 AND N_ELEMENTS(h_odd_index) EQ 0 THEN h_odd_index = [h_i+1] ELSE IF (h_i MOD 2) EQ 0 THEN h_odd_index = [h_odd_index, h_i+1]
    IF h_i EQ 1 AND N_ELEMENTS(h_even_index) EQ 0 THEN h_even_index = [h_i+1] ELSE IF (h_i MOD 2) EQ 1 THEN h_even_index = [h_even_index, h_i+1]    
    
    ;@GJ, 2024/5/21, plot the middle one
    signal_array_0_noDC = REFORM((*parameters.Sx_array)[h_i, harmonics_comp, zero_ind:*, 0])
    t_ele = N_Elements(signal_array_0_noDC)
    f=1.
    delta_t = (temp_A_G_array[1]-temp_A_G_array[0]) * DOUBLE(f) ; us; samplying rate: 1MHz
    X = (FINDGEN((t_ele - 1)/2) + 1)
    is_N_even = (t_ele MOD 2) EQ 0
    if (is_N_even) then $
      frequency = [0.0, X, t_ele/2, -t_ele/2 + X]/(t_ele*delta_t) $
    else $
      frequency = [0.0, X, -(t_ele/2 + 1) + X]/(t_ele*delta_t)
    
;    IF flag_i EQ 0 THEN BEGIN
;      iplot, frequency, ABS(FFT(signal_array_0_noDC)), xrange=[0, 20./(t_ele*delta_t)], color='red', title='FFT of focal spots'
;      flag_i++
;    ENDIF ELSE BEGIN
;      iplot, frequency, ABS(FFT(signal_array_0_noDC)), color='green', /overplot
;    ENDELSE
    
    signal_array_0 = REFORM((*parameters.Sx_array)[h_i, harmonics_comp, zero_ind:*, index_A_DC-1])
    ;@GJ, 2024/5/31, calculate the energy
    h_energy_0[h_i] = MEAN(ABS(signal_array_0))
    maxD_Harmonic_0 = MAX(signal_array_0, maxD_ind_Harmonic_0)
    IF maxD_ind_Harmonic_0 GT 3. THEN BEGIN
      h_radius_0[h_i] = temp_A_G_array[maxD_ind_Harmonic_0]
    ENDIF; ELSE BEGIN
      ;non-donut, calculate FWHM
      FWHM_ind_0 = INTERPOL(temp_A_G_array, signal_array_0, 0.5*maxD_Harmonic_0)
      IF finite(FWHM_ind_0) THEN h_FWHM_0[h_i] = 2. * FWHM_ind_0
      h_MTF_0[h_i] = 1./h_FWHM_0[h_i]
;    ENDELSE

    signal_array_1 = REFORM((*parameters.Sx_array)[h_i, harmonics_comp, zero_ind:*, index_A_DC_1-1])
    ;@GJ, 2024/5/31, calculate the energy
    h_energy_1[h_i] = MEAN(ABS(signal_array_1))
    maxD_Harmonic_1 = MAX(signal_array_1, maxD_ind_Harmonic_1)
    IF maxD_ind_Harmonic_1 GT 3. THEN BEGIN
      h_radius_1[h_i] = temp_A_G_array[maxD_ind_Harmonic_1]
    ENDIF; ELSE BEGIN
      ;non-donut, calculate FWHM
      FWHM_ind_1 = INTERPOL(temp_A_G_array, signal_array_1, 0.5*maxD_Harmonic_1)
      IF finite(FWHM_ind_1) THEN h_FWHM_1[h_i] = 2. * FWHM_ind_1
      h_MTF_1[h_i] = 1./h_FWHM_1[h_i]
;    ENDELSE
  ENDFOR
;  iplot, h_odd_index, h_radius_0[h_odd_index-1], psym=0, /NO_SAVEPROMPT
;  iplot, h_odd_index, h_FWHM_0[h_odd_index-1], psym=1, /overplot
;  iplot, h_even_index, h_radius_0[h_even_index-1], psym=2, color='ff0000'x, /overplot
;  iplot, h_even_index, h_FWHM_0[h_even_index-1], psym=3, color='00ff00'x, /overplot
  
  ;@GJ, 2024/2/11, plotting the harmonics based on even or odd indices
  drawXsize = parameters.drawXSize
  drawYsize = parameters.drawYSize
  IF wid GT 0 AND parameters.signal_type NE 2 THEN BEGIN
    draw3 = widget_info(wid, find_by_uname='draw3')
    widget_control, draw3, get_value = drawWindow3
    wset, drawWindow3
;    wdelete, drawWindow3
    tvscl, BYTARR(drawXsize, drawYsize)*0.
    
    LOADCT, 1
    draw13 = widget_info(wid, find_by_uname='draw13')
    widget_control, draw13, get_value = drawWindow13
    wset, drawWindow13
;    wdelete, drawWindow13
    tvscl, BYTARR(drawXsize, drawYsize)*0.
  ENDIF
  IF wid GT 0 THEN BEGIN ;AND parameters.signal_type EQ 2
    ;    DEVICE, DECOMPOSED = 0, RETAIN = 2
    ;    LOADCT, 3
    draw3 = widget_info(wid, find_by_uname='draw3')
    widget_control, draw3, get_value = drawWindow3
    wset, drawWindow3
    ;display odd harmonics
    IF harmonics_n MOD 2 EQ 1 THEN BEGIN
      maxy = MAX([h_radius_0[h_odd_index-1], h_FWHM_0[h_odd_index-1]])
      miny = MIN([h_radius_0[h_odd_index-1], h_FWHM_0[h_odd_index-1]])
      plot, h_odd_index, h_FWHM_0[h_odd_index-1], /nodata, thick=1, xrange=[1, parameters.harmonics_n_max], xstyle=1, yrange=[0., maxy], xtitle='Harmonic Index', title='FWHM (Odd)', ytitle='A_G [mT]';, title='PSF';, ytitle='Signal [A.U.]';background = 'ffffff'x, color = '000000'x,
      oplot, h_odd_index, h_FWHM_0[h_odd_index-1], thick=1, color = '0000ff'x
      oplot, h_odd_index, h_radius_0[h_odd_index-1], thick=1, linestyle=1, color = '0000ff'x
      oplot, [harmonics_n], [h_FWHM_0[harmonics_n-1]], thick=1, psym=6, color = '0000ff'x
      oplot, [harmonics_n], [h_radius_0[harmonics_n-1]], thick=1, psym=6, color = '0000ff'x
      ;@GJ, 2024/6/1, plotting the MTF and SNR
;      Additional_Axes_Plot_GJ, h_odd_index[0:6], (h_MTF_0[h_odd_index-1])[0:6], (h_energy_0[h_odd_index-1])[0:6], 1./(h_radius_0[h_odd_index-1])[0:6], 'Harmonic Index', 'Resolution [1/mm]', 'SNR', 'Donut Size [1/mm]'
    ENDIF ELSE BEGIN
      ;display even harmonics
      maxy = MAX([h_radius_0[h_even_index-1], h_FWHM_0[h_even_index-1]])
      miny = MIN([h_radius_0[h_even_index-1], h_FWHM_0[h_even_index-1]])
      plot, h_even_index, h_radius_0[h_even_index-1], /nodata, thick=1, xrange=[1, parameters.harmonics_n_max], xstyle=1, yrange=[0., maxy], xtitle='Harmonic Index', title='FWHM (Even)', ytitle='A_G [mT]';, title='PSF';, ytitle='Signal [A.U.]';background = 'ffffff'x, color = '000000'x,
      oplot, h_even_index, h_radius_0[h_even_index-1], thick=1, linestyle=1, color = '0000ff'x
      oplot, h_even_index, h_FWHM_0[h_even_index-1], thick=1, color = '0000ff'x
      oplot, [harmonics_n], [h_FWHM_0[harmonics_n-1]], thick=1, psym=6, color = '0000ff'x
      oplot, [harmonics_n], [h_radius_0[harmonics_n-1]], thick=1, psym=6, color = '0000ff'x
      ;@GJ, 2024/6/1, plotting the MTF and SNR
;      Additional_Axes_Plot_GJ, h_even_index[0:6], (h_MTF_0[h_even_index-1])[0:6], (h_energy_0[h_even_index-1])[0:6], 1./(h_radius_0[h_even_index-1])[0:6], 'Harmonic Index', 'Resolution [1/mm]', 'SNR', 'Donut Size [1/mm]'
    ENDELSE

    LOADCT, 1
    draw13 = widget_info(wid, find_by_uname='draw13')
    widget_control, draw13, get_value = drawWindow13
    wset, drawWindow13
    ;display odd harmonics
    IF harmonics_n MOD 2 EQ 1 THEN BEGIN
      maxy = MAX([h_radius_1[h_odd_index-1], h_FWHM_1[h_odd_index-1]])
      miny = MIN([h_radius_1[h_odd_index-1], h_FWHM_1[h_odd_index-1]])
      plot, h_odd_index, h_FWHM_1[h_odd_index-1], /nodata, thick=1, xrange=[1, parameters.harmonics_n_max], xstyle=1, yrange=[0., maxy], xtitle='Harmonic Index', title='Radius (Odd)', ytitle='A_G [mT]';, title='PSF';, ytitle='Signal [A.U.]';background = 'ffffff'x, color = '000000'x,
      oplot, h_odd_index, h_FWHM_1[h_odd_index-1], thick=1, linestyle=1, color = '00ff00'x
      oplot, h_odd_index, h_radius_1[h_odd_index-1], thick=1, color = '00ff00'x
      oplot, [harmonics_n], [h_FWHM_1[harmonics_n-1]], thick=1, psym=6, color = '00ff00'x
      oplot, [harmonics_n], [h_radius_1[harmonics_n-1]], thick=1, psym=6, color = '00ff00'x
    ENDIF ELSE BEGIN
      ;display even harmonics
      maxy = MAX([h_radius_1[h_even_index-1], h_FWHM_1[h_even_index-1]])
      miny = MIN([h_radius_1[h_even_index-1], h_FWHM_1[h_even_index-1]])
      plot, h_even_index, h_radius_1[h_even_index-1], /nodata, thick=1, xrange=[1, parameters.harmonics_n_max], xstyle=1, yrange=[0., maxy], xtitle='Harmonic Index', title='Radius (Even)', ytitle='A_G [mT]';, title='PSF';, ytitle='Signal [A.U.]';background = 'ffffff'x, color = '000000'x,
      oplot, h_even_index, h_radius_1[h_even_index-1], thick=1, color = '00ff00'x
      oplot, h_even_index, h_FWHM_1[h_even_index-1], thick=1, linestyle=1, color = '00ff00'x
      oplot, [harmonics_n], [h_FWHM_1[harmonics_n-1]], thick=1, psym=6, color = '00ff00'x
      oplot, [harmonics_n], [h_radius_1[harmonics_n-1]], thick=1, psym=6, color = '00ff00'x
    ENDELSE
    DEVICE, DECOMPOSED = 1
  ENDIF
   
  ;@GJ, 2024/1/2, calculating the donut array
  curr_Sx_array_donut = curr_Sx_array;REFORM((*parameters.Sx_array)[harmonics_n-1, harmonics_comp, *, *])
  FOR id=0, n_z-1 DO BEGIN
    focal_line_temp = REFORM(curr_Sx_array[*, id])
    donut_Harmonic_temp = REFORM(focal_line_temp[zero_ind:*])
    maxD_Harmonic_temp = MAX(donut_Harmonic_temp, maxD_ind_Harmonic_temp)
    ;@GJ, 2024/1/2, if non-donut, removing it
    IF maxD_ind_Harmonic_temp LE 3. THEN curr_Sx_array_donut[*, id] = MIN(curr_Sx_array)
  ENDFOR
  
  ;@GJ, 2023/12/31, calculate the radius and thick
  donut_Harmonic_1 = REFORM(focal_line_1[zero_ind:*])
  maxD_Harmonic_1 = MAX(donut_Harmonic_1, maxD_ind_Harmonic_1)
  donut_flag_1 = 0
  donut_radius_1 = 0.
  donut_thick_1 = 0.
  IF maxD_ind_Harmonic_1 GT 3. THEN BEGIN
    donut_flag_1 = 1
    donut_radius_1 = temp_A_G_array[maxD_ind_Harmonic_1]
    donut_thick_1 = donut_Harmonic_1[maxD_ind_Harmonic_1] - donut_Harmonic_1[0]
  ENDIF ELSE BEGIN
    ;non-donut, calculate FWHM
    FWHM_ind_1 = INTERPOL(temp_A_G_array, donut_Harmonic_1, 0.5*maxD_Harmonic_1)
    IF finite(FWHM_ind_1) THEN BEGIN
      donut_flag_1 = 0
      donut_radius_1 = 2. * FWHM_ind_1
      donut_thick_1 = maxD_Harmonic_1
    ENDIF  
  ENDELSE
  
  ;@GJ, 2023/12/31, defining a focal ratio
;  focal_ratio = 1./maxD_Harmonic_1*donut_Harmonic_0[maxD_ind_Harmonic_1]
;  focal_ratio = total(donut_Harmonic_0)/total(donut_Harmonic_1)
;  focal_ratio = 1.
;  donut_Harmonic_1 = donut_Harmonic_1 * focal_ratio
  ;@GJ, 2023/1/4, calculating the subtracted dont harmonic
  ori_x = temp_A_G_array[maxD_ind_Harmonic_1:*]
  ori_y = donut_Harmonic_0[maxD_ind_Harmonic_1:*]
  don_y = donut_Harmonic_1[maxD_ind_Harmonic_1:*]
;  y_1_new = INTERPOL(ori_y, don_y, donut_Harmonic_1)
;  donut_Harmonic_s = donut_Harmonic_0-y_1_new;donut_Harmonic_0 - donut_Harmonic_1 * focal_ratio
  y_1_new = donut_Harmonic_1 * parameters.dog_factor/10.
  donut_Harmonic_s = donut_Harmonic_0 - donut_Harmonic_1 * parameters.dog_factor/10.
  ;  iplot, donut_Harmonic_s, title='focal spot', /NO_SAVEPROMPT
  ;@GJ, 2023/12/31, calculate the radius and thick
  ;@GJ, 2023/1/4, calculating the conversion factor
  donut_flag_s = 0
  donut_radius_s = 0.
  donut_thick_s = 0.
  maxD_Harmonic_s = MAX(donut_Harmonic_s, maxD_ind_Harmonic_s)
  IF maxD_ind_Harmonic_s GT 3. THEN BEGIN
    donut_flag_s = 1
    donut_radius_s = temp_A_G_array[maxD_ind_Harmonic_s]
    donut_thick_s = donut_Harmonic_s[maxD_ind_Harmonic_s] - donut_Harmonic_s[0]
  ENDIF ELSE BEGIN
    ;non-donut, calculate FWHM
    FWHM_ind_s = INTERPOL(temp_A_G_array, donut_Harmonic_s, 0.5*maxD_Harmonic_s)
    IF finite(FWHM_ind_s) THEN BEGIN
      donut_flag_s = 0
      donut_radius_s = 2. * FWHM_ind_s
      donut_thick_s = maxD_Harmonic_s
    ENDIF
  ENDELSE
  
  ;@GJ, 2024/1/10, plot the radius and FWHM
  A_DC_FWHM = DBLARR(n_x)*0.
  A_DC_radius = DBLARR(n_x)*0.
  A_DC_thickness = DBLARR(n_x)*0.
  FOR A_DC_i=0, n_x-1 DO BEGIN
    signal_array = curr_Sx_array[zero_ind:*, A_DC_i];REFORM((*parameters.Sx_array)[harmonics_n-1, harmonics_comp, zero_ind:*, A_DC_i])
    maxD_Harmonic = MAX(signal_array, maxD_ind_Harmonic)
    IF maxD_ind_Harmonic GE 1. THEN BEGIN
      A_DC_radius[A_DC_i] = temp_A_G_array[maxD_ind_Harmonic]
      A_DC_thickness[A_DC_i] = maxD_Harmonic - signal_array[n_z/2-1]
    ENDIF ELSE BEGIN
      ;non-donut, calculate FWHM
      FWHM_ind = INTERPOL(temp_A_G_array, signal_array, 0.5*maxD_Harmonic)
      IF finite(FWHM_ind) THEN A_DC_FWHM[A_DC_i] = 2. * FWHM_ind
    ENDELSE
  ENDFOR
  
  ;@GJ, 2024/5/1, calculating the threshold of donut
  min_90 = MIN(ABS(A_DC_array/A*100. - 90.), min_90_ind)
  FOR i_90 = min_90_ind, n_x-1 DO BEGIN
    IF A_DC_FWHM[i_90] EQ 0 THEN BEGIN
      parameters.A_DC_threshold_donut = A_DC_array[i_90]/A*100.
      print, 'threshold: ', parameters.A_DC_threshold_donut
      ;@GJ, 2024/5/5, calculate the slope
      linfit_result = LINFIT(A_DC_array[i_90:*]-A_DC_array[i_90], A_DC_radius[i_90:*])
      parameters.donut_radius0 = A_DC_radius[i_90];linfit_result[0];
      parameters.donut_slope = linfit_result[1];MEAN(A_DC_radius[i_90:*]-A_DC_radius[i_90])/MEAN(A_DC_array[i_90:*]-A_DC_array[i_90]); mT/mT;/linfit_result[1]; slope
;      iplot, A_DC_array[i_90:*]-A_DC_array[i_90], A_DC_radius[i_90:*], color='black', /NO_SAVEPROMPT
;      iplot, A_DC_array[i_90:*]-A_DC_array[i_90], parameters.donut_slope*(A_DC_array[i_90:*]-A_DC_array[i_90])+A_DC_radius[i_90], color='red', /overplot
      break;
    ENDIF
  ENDFOR
  
;  nterms = 3.;6.;4.
;  gauss_fit_0 = GAUSSFIT(temp_A_G_array, donut_Harmonic_0, coeff, NTERMS=nterms)
;  x_1 = temp_A_G_array[maxD_ind_Harmonic_1]
;  y_1 = donut_Harmonic_1[maxD_ind_Harmonic_1]
;  y_0 = donut_Harmonic_0[maxD_ind_Harmonic_1]
;  gauss_fit_1 = GAUSSFIT(temp_A_G_array[maxD_ind_Harmonic_1:*], donut_Harmonic_1[maxD_ind_Harmonic_1:*], coeff_1, NTERMS=nterms)
;  x_1_new = SQRT(-ALOG(donut_Harmonic_1/coeff_1[0])*2.*coeff_1[2])-coeff_1[1]
;  z_1_new = (x_1_new-coeff[1])/coeff[2]
;  y_1_new = coeff[0] * EXP(-z_1_new^2/2.)
    
;  iplot, REFORM(parameters.Sx_array[255, *]), title='Cheybshev ('+STRTRIM(STRING(particle_size),1)+' nm)', /NO_SAVEPROMPT
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 0
    draw1 = widget_info(wid, find_by_uname='draw1')
    widget_control, draw1, get_value = drawWindow1
    wset, drawWindow1
;    Sx_array_focal = curr_Sx_array
;    Sx_array_focal[*, ] = MAX(curr_Sx_array)
    width_ratio = (n_x*gradient_field_x)/(n_z*gradient_field_yz)
    temp_image = ROTATE(CONGRID(curr_Sx_array, drawXSize, drawYSize), 1)
;    IF gradient_field_x LT gradient_field_yz THEN BEGIN
;      temp_image = ROTATE(CONGRID(curr_Sx_array[(n_x/2.-0.5*width_ratio*n_x):n_x/2.+0.5*width_ratio*n_x, *], drawXSize, drawYSize),1)
;    ENDIF ELSE BEGIN
;      temp_map = DBLARR(n_z*width_ratio, n_z)
;      temp_map[(n_z*width_ratio/2.-0.5*n_z):n_z*width_ratio/2.+0.5*n_z-1, *] = curr_Sx_array
;      temp_image = ROTATE(CONGRID(temp_map, drawXSize, drawYSize),1)
;      ;temp_image = ROTATE(CONGRID(curr_Sx_array[*, (n_z/2.-0.5*width_ratio*n_z):n_z/2.+0.5*width_ratio*n_z], drawXSize, drawYSize),1)
;    ENDELSE
    ;temp_image = ROTATE(CONGRID(curr_Sx_array, drawXSize, drawYSize),1)
    ;@GJ， 2025/1/7， plot the multi-channel X-space
    ;iimage, (CONGRID(curr_Sx_array, drawXSize, drawYSize) + temp_image)/2.
    ;red color for A_DC
    topClr_1 = !D.TABLE_SIZE - 1
    TVLCT, 255, 0, 0, topClr_1
    ;green color for A_DC_1
    topClr_2 = !D.TABLE_SIZE - 2
    TVLCT, 0, 255, 0, topClr_2
    topClr_3 = !D.TABLE_SIZE - 3
    TVLCT, 255, 0, 255, topClr_3
    temp_image = BYTSCL(temp_image, MAX=MAX(curr_Sx_array[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(curr_Sx_array[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), TOP=(topClr_3-1))
    ;@GJ, 2024/7/11, plot the sum curves
;    iplot, REFORM(temp_image[92,*]), color='blue', /NO_SAVEPROMPT
    max_temp = MAX(REFORM(temp_image[92,*]))
    sum_curve = DBLARR(183)
    temp_image_t = temp_image
    temp_image_t[0:162,*]=0
    for i=0,182 do sum_curve[i]= TOTAL(temp_image_t[*,i])
;    iplot, sum_curve/MAX(sum_curve)*max_temp, color='red',/overplot
    
    gradient_yz_ratio = parameters.gradient_field_yz_current / parameters.gradient_field_yz
    IF gradient_yz_ratio LT 1 THEN BEGIN
      temp_image[*, drawXSize/2. * (1 + gradient_yz_ratio)] = topClr_3
      temp_image[*, drawXSize/2. * (1 - gradient_yz_ratio)] = topClr_3
    ENDIF
    temp_image[drawYSize / 2. + (index_A_DC-n_x/2.)/n_x*drawYSize, *] = topClr_1
    temp_image[drawYSize / 2. + (index_A_DC_1-n_x/2.)/n_x*drawYSize, *] = topClr_2
    TV, temp_image
    IF harmonics_comp EQ 0 THEN map_title = 'Harmonic Amp'
    IF harmonics_comp EQ 1 THEN map_title = 'Harmonic Phase'
    IF harmonics_comp EQ 2 THEN map_title = 'Harmonic Real'
    IF harmonics_comp EQ 3 THEN map_title = 'Harmonic Imaginary'
    xyouts,drawysize*0.23,drawysize*0.9, map_title, color = (topClr_3-1), /DEVICE
    FOV_harmonic_map = parameters.A * 2. * (parameters.A_DC_max/100.) / gradient_field_yz_current
    xyouts,drawysize*0.03,drawysize*0.8,'Displayed FOV: '+STRING(FOV_harmonic_map, format='(f7.2)')+' mm', color = (topClr_3-1), /DEVICE
    
    DEVICE, DECOMPOSED = 1
    draw11 = widget_info(wid, find_by_uname='draw11')
    widget_control, draw11, get_value = drawWindow11
    wset, drawWindow11
    maxy = MAX([A_DC_radius, A_DC_FWHM]);50.;
    miny = MIN([A_DC_radius, A_DC_FWHM])

    plot, A_DC_array/A*100., A_DC_FWHM, /nodata, thick=2, xrange=[0, parameters.A_DC_max], xstyle=1, yrange=[0., maxy], xtitle='A_DC/A [%]', title='FWHM & Radius', ytitle='A_G [mT]';, yrange=[0., 50.];, title='PSF';, ytitle='Signal [A.U.]';background = 'ffffff'x, color = '000000'x,
    oplot, A_DC_array/A*100., A_DC_FWHM, thick=1, color = '0000ff'x
    oplot, A_DC_array/A*100., A_DC_radius, thick=1, color = '00ff00'x
;    oplot, A_DC_array/A*100., A_DC_thickness/MAX(A_DC_thickness)*maxy, thick=1, linestyle=2, color = 'ffff00'x
    oplot, [A_DC_array[index_A_DC-1]/A*100.], [A_DC_FWHM[index_A_DC-1]], thick=1, psym=4, color = '0000ff'x
;    oplot, [A_DC_array[index_A_DC]/A*100.], [A_DC_radius[index_A_DC]], thick=2, psym=6, color = '0000ff'x
;    oplot, [A_DC_array[index_A_DC_1]/A*100.], [A_DC_FWHM[index_A_DC_1]], thick=2, psym=6, color = 'ffb037'x
    oplot, [A_DC_array[index_A_DC_1-1]/A*100.], [A_DC_radius[index_A_DC_1-1]], thick=1, psym=6, color = '00ff00'x
;    oplot, [A_DC_array[index_A_DC_1]/A*100.], [A_DC_thickness[index_A_DC_1]]/MAX(A_DC_thickness)*maxy, thick=1, psym=6, color = 'ffff00'x

;    ;    Sx_array_focal = curr_Sx_array
;    ;    Sx_array_focal[*, ] = MAX(curr_Sx_array)
;    IF gradient_field_x LT gradient_field_yz THEN BEGIN
;      width_ratio = (n_x*gradient_field_x)/(n_z*gradient_field_yz)
;      temp_image = ROTATE(CONGRID(curr_Sx_array_donut[(n_x/2.-0.5*width_ratio*n_x):n_x/2.+0.5*width_ratio*n_x, *], drawXSize, drawYSize),1)
;    ENDIF ELSE BEGIN
;      width_ratio = (n_z*gradient_field_yz)/(n_x*gradient_field_x)
;      temp_image = ROTATE(CONGRID(curr_Sx_array_donut[*, (n_z/2.-0.5*width_ratio*n_z):n_z/2.+0.5*width_ratio*n_z], drawXSize, drawYSize),1)
;    ENDELSE
;;    temp_image = ROTATE(CONGRID(curr_Sx_array, drawXSize, drawYSize),1)
;    temp_image[drawYSize / 2. + (index_A_DC_1-n_x/2.)/n_x*drawYSize, *] = MAX(curr_Sx_array_donut)
;    TVSCL, BYTSCL(temp_image, MAX=MAX(curr_Sx_array_donut[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(curr_Sx_array_donut[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;    xyouts,drawysize*0.30,drawysize*0.9,'Donut Map', color = 'FFFFFF'x, /DEVICE
;    xyouts,drawysize*0.03,drawysize*0.8,'Displayed FOV: '+STRING(width_ratio*FOV, format='(f6.2)')+' mm', color = 'FFFFFF'x, /DEVICE
    ;    DEVICE, DECOMPOSED = 1
    
    ;@GJ, 2023/12/31, plotting the PSF
    draw21 = widget_info(wid, find_by_uname='draw21')
    widget_control, draw21, get_value = drawWindow21
    wset, drawWindow21
    maxy = MAX([donut_Harmonic_0, donut_Harmonic_1, donut_Harmonic_s])
    miny = MIN([donut_Harmonic_0, donut_Harmonic_1, donut_Harmonic_s])
    min_abs = MIN(ABS(donut_Harmonic_0-0.1*MAX(donut_Harmonic_0)), min_abs_ind)
    IF min_abs_ind LT n_z/10. THEN min_abs_ind = n_z/10.*2.
    plot, A_G_array, focal_line, /nodata, thick=1, xrange=[-temp_A_G_array[min_abs_ind], temp_A_G_array[min_abs_ind]], xstyle=1, yrange=[miny, maxy], xtitle='A_G [mT]', title='PSFs';, title='PSF';, ytitle='Signal [A.U.]';background = 'ffffff'x, color = '000000'x, 
    oplot, A_G_array, focal_line, thick=1, color = '0000ff'x
    oplot, A_G_array, focal_line_1, thick=1, color = '00ff00'x
    oplot, A_G_array, focal_line-focal_line_1*parameters.dog_factor/10., thick=1, color = 'ffb037'x
;    oplot, A_G_array[maxD_ind_Harmonic_1:*], gauss_fit_1, thick=2, color = 'ff00ff'x
    oplot, A_G_array, focal_line_1*parameters.dog_factor/10., thick=1, linestyle=1, color = '00ff00'x
;    oplot, A_G_array, , thick=2, color = '00ff00'x
  ENDIF
 
  ;@GJ, 2024/5/10, adjust the PSF based on the gradient_field_yz_current
  gradient_ratio = parameters.gradient_field_yz_current / parameters.gradient_field_yz
  focal_image_0 = curr_Sx_array * 0.
  focal_image_0_real = curr_Sx_array * 0.
  focal_image_0_imag = curr_Sx_array * 0.
  focal_image_1 = curr_Sx_array * 0.
  focal_image_1_real = curr_Sx_array * 0.
  focal_image_1_imag = curr_Sx_array * 0.
  FOR i=0, n_x-1 DO BEGIN
    FOR j=0, n_z-1 DO BEGIN
      focal_image_0[i, j] = INTERPOL([focal_line, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
      focal_image_0_real[i, j] = INTERPOL([focal_line_real, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
      focal_image_0_imag[i, j] = INTERPOL([focal_line_imag, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
      focal_image_1[i, j] = INTERPOL([focal_line_1, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
      focal_image_1_real[i, j] = INTERPOL([focal_line_1_real, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
      focal_image_1_imag[i, j] = INTERPOL([focal_line_1_imag, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio, /LSQUADRATIC)
    ENDFOR
  ENDFOR
  focal_image_sub = focal_image_0 - focal_image_1*parameters.dog_factor/10.;INTERPOL(ori_y, don_y, focal_image_1)
  
  gradient_ratio_display = gradient_ratio * parameters.FOV_current / parameters.FOV
  focal_image_0_display = focal_image_0 * 0.
  focal_image_1_display = focal_image_1 * 0.
  FOR i=0, n_x-1 DO BEGIN
    FOR j=0, n_z-1 DO BEGIN
      focal_image_0_display[i, j] = INTERPOL([focal_line, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio_display, /LSQUADRATIC)
      focal_image_1_display[i, j] = INTERPOL([focal_line_1, 0, 0, 0], [focal_coor, n_x*gradient_ratio, 1.5*n_x*gradient_ratio, 2*n_x*gradient_ratio], SQRT((i-n_x/2.)^2 + (j-n_z/2.)^2)*gradient_ratio_display, /LSQUADRATIC)
    ENDFOR
  ENDFOR
  focal_image_sub_display = focal_image_0_display - focal_image_1_display*parameters.dog_factor/10.;INTERPOL(ori_y, don_y, focal_image_1)
     
  ;@GJ, 2024/1/22, plotting the images
  test_ima = focal_image_0[128-50:128+50, 128-20:128+20]
  little = SHIFT(test_ima, 3, 0) + SHIFT(test_ima, -3, 0)
  small = SHIFT(test_ima, 5, 0) + SHIFT(test_ima, -5, 0)
  large = SHIFT(test_ima, 8, 0) + SHIFT(test_ima, -8, 0)
;  iimage, little, 'close'
;  iimage, small, 'near'
;  iimage, large, 'far'
;  iplot, little(*, 20), /NO_SAVEPROMPT
;  iplot, small(*, 20), /NO_SAVEPROMPT
;  iplot, large(*, 20), /NO_SAVEPROMPT
  
  ;calculate the system matrix
  SM_0 = DBLARR(n_z, n_z)
  SM_1 = DBLARR(n_z, n_z)
  FOR j=0, n_z-1 DO BEGIN
    FOR k=0, n_z-1 DO BEGIN
      pos = -k + j
      SM_0[j,k] = INTERPOL(focal_line, focal_coor, pos)
      SM_1[j,k] = INTERPOL(focal_line_1, focal_coor, pos)
    ENDFOR
  ENDFOR
  
;  SM_1 *= MAX(SM_0)/MAX(SM_1)
  
  ;calculate the SVD
  ; Compute the Singular Value Decomposition:
  LA_SVD, SM_0, W0, U0, V0
;  iplot, w0, thick=2, title='SVD of System Matrix 0', xtitle='Component Index', ytitle='Singular Value', /ylog, xrange=[1,N_ELEMENTS(w0)], /NO_SAVEPROMPT
  poly_fit_result0 = POLY_FIT(FINDGEN(N_ELEMENTS(w0))+1, ALOG(w0), 10,sigma=sigma, yfit=yw0)
  svd_slope0 = ABS(poly_fit_result0[1])
  print, 'svd_slope 0: ', svd_slope0
  
  lA_SVD, SM_1, W1, U1, V1
;  iplot, w1, thick=2, title='SVD of System Matrix 1', xtitle='Component Index', ytitle='Singular Value', /ylog, xrange=[1,N_ELEMENTS(w1)], /NO_SAVEPROMPT
  poly_fit_result1 = POLY_FIT(FINDGEN(N_ELEMENTS(w1))+1, ALOG(w1), 10,sigma=sigma, yfit=yw1)
  svd_slope1 = ABS(poly_fit_result1[1])
  print, 'svd_slope 1: ', svd_slope1
   
  ;SM_2 = [[SM_0], [SM_1]]
  SM_2 = SM_0 - SM_1*parameters.dog_factor/10.
  LA_SVD, SM_2, W2, U2, V2
  x_index = FINDGEN(N_ELEMENTS(w2))+1.
  poly_fit_result2 = POLY_FIT(FINDGEN(N_ELEMENTS(w2))+1, ALOG(w2), 10,sigma=sigma, yfit=yw1)
  svd_slope2 = ABS(poly_fit_result2[1])
  print, 'svd_slope 1: ', svd_slope2
;  iplot, x_index, w2, thick=2, color='green', title='SVD of System Matrix ('+STRING(parameters.A_DC, format='(i0)')+'%+'+STRING(parameters.A_DC_1, format='(i0)')+'%)', xtitle='Component Index', ytitle='Singular Value', /ylog, /NO_SAVEPROMPT
;  iplot, x_index, w0, /ylog, color='red', xrange=[1,N_ELEMENTS(w2)], /overplot
;  iplot, x_index, w1, /ylog, color='blue', xrange=[1,N_ELEMENTS(w2)], /overplot
;  iplot, x_index, (w2-w0)/MAX(w0)*100., thick=2, color='green', title='SVD Diff of System Matrix ('+STRING(parameters.A_DC, format='(i0)')+'%+'+STRING(parameters.A_DC_1, format='(i0)')+'%)', xtitle='Component Index [a.u.]', ytitle='SV Diff [%]';, /xlog, /NO_SAVEPROMPT
;  IF wid GT 0 THEN BEGIN
;    ;@GJ, 2023/12/31, plotting the PSF
;    draw25 = widget_info(wid, find_by_uname='draw25')
;    widget_control, draw25, get_value = drawWindow25
;    wset, drawWindow25
;;    plot, x_index, w2, /nodata, title='SVD of System Matrix ('+STRING(parameters.A_DC, format='(i0)')+'%+'+STRING(parameters.A_DC_1, format='(i0)')+'%)', xtitle='Component Index', ytitle='Singular Value', /xlog, /ylog
;    plot, x_index, w2, /nodata, title='SVD of SM', xtitle='Component Index', ytitle='Singular Value', /ylog
;    oplot, x_index, w0, color='0000ff'x
;    oplot, x_index, w1, color='00ff00'x
;    oplot, x_index, w2, color='ffb037'x;, thick=2
;  ENDIF
  
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 3
    draw4 = widget_info(wid, find_by_uname='draw4')
    widget_control, draw4, get_value = drawWindow4
    wset, drawWindow4
;    TVSCL, BYTSCL(CONGRID(SM_0, drawXSize, drawYSize), MAX=MAX(SM_0[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(SM_0[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;    xyouts,drawysize*0.35,drawysize*0.9,'SM (Regular)', color = 'FFFFFF'x, /DEVICE
    focal_image_0_FFT_alog = ALOG(ABS(FFT(focal_image_0,/center,/double)))
    MTF_peak_1st = parameters.FOV / (donut_radius_0/gradient_field_yz_current) / drawXSize
    temp_display = BYTSCL(CONGRID(focal_image_0_FFT_alog, drawXSize, drawYSize), MAX=MAX(focal_image_0_FFT_alog[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
   ;@GJ, 2024/6/3, labe the central frequency
    A_temp = temp_display * 0.
    B_temp = temp_display * 0.
    IF drawXSize/2.*(1.+MTF_peak_1st)+1 LE drawXSize-1 THEN A_temp[drawXSize/2.*(1.+MTF_peak_1st)+1, drawYSize/2.+1] = 255.
;    A_temp[drawXSize/2.*(1.-MTF_peak)+1, drawYSize/2.+1] = 255.
    FOR i_r = 0, 360 DO B_temp += ROT(A_temp, i_r, 1., /INTERP)
    MTF_indices_1st = WHERE(B_temp GT 128)
    temp_display[MTF_indices_1st] = 255
    TVSCL, temp_display
    xyouts,drawysize*0.30,drawysize*0.9,'1st PSF logFFT', color = 'FFFFFF'x, /DEVICE

    LOADCT, 8
    draw14 = widget_info(wid, find_by_uname='draw14')
    widget_control, draw14, get_value = drawWindow14
    wset, drawWindow14
;    TVSCL, BYTSCL(CONGRID(SM_1, drawXSize, drawYSize), MAX=MAX(SM_1[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(SM_1[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;    xyouts,drawysize*0.35,drawysize*0.9,'SM (Donut)', color = 'FFFFFF'x, /DEVICE
    focal_image_1_FFT_alog = ALOG(ABS(FFT(focal_image_1,/center,/double)))
    TVSCL, BYTSCL(CONGRID(focal_image_1_FFT_alog, drawXSize, drawYSize), MAX=MAX(focal_image_1_FFT_alog[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.30,drawysize*0.9,'2nd PSF logFFT', color = 'FFFFFF'x, /DEVICE

    LOADCT, 64
    draw24 = widget_info(wid, find_by_uname='draw24')
    widget_control, draw24, get_value = drawWindow24
    wset, drawWindow24
    focal_image_sub_FFT_alog = ALOG(ABS(FFT(focal_image_sub,/center,/double)))
    STED_FWHM = donut_radius_s/gradient_field_yz_current
;    ;    STED_MTF = 1./STED_FWHM
    MTF_peak_sub = parameters.FOV / STED_FWHM / drawXSize
    temp_display = BYTSCL(CONGRID(focal_image_sub_FFT_alog, drawXSize, drawYSize), MAX=MAX(focal_image_sub_FFT_alog[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    ;@GJ, 2024/6/3, labe the central frequency
    A_temp = temp_display * 0.
    B_temp = temp_display * 0.
    IF drawXSize/2.*(1.+MTF_peak_sub)+1 GT 0 AND drawXSize/2.*(1.+MTF_peak_sub)+1 LT FLOOR(drawXSize) THEN A_temp[drawXSize/2.*(1.+MTF_peak_sub)+1, drawYSize/2.+1] = 255.
    ;    A_temp[drawXSize/2.*(1.-MTF_peak)+1, drawYSize/2.+1] = 255.
    FOR i_r = 0, 360 DO B_temp += ROT(A_temp, i_r, 1., /INTERP)
    MTF_indices_sub = WHERE(B_temp GT 128)
    temp_display[MTF_indices_sub] = 0
    TVSCL, temp_display
;    TVSCL, 255-BYTSCL(CONGRID(SM_2, drawXSize, drawYSize), MAX=MAX(SM_2[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(SM_2[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;    xyouts,drawysize*0.35,drawysize*0.9,'SM (Subtr)', color = '000000'x, /DEVICE
;    TVSCL, BYTSCL(CONGRID(focal_image_sub_FFT_alog, drawXSize, drawYSize), MAX=MAX(focal_image_sub_FFT_alog[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.30,drawysize*0.9,'Sub PSF logFFT', color = 'FFFFFF'x, /DEVICE


    DEVICE, DECOMPOSED = 1
  ENDIF
  
  ;do loop iterative reconstruction for convergence reconstruction
  bar_image_1d=image_phantom_res(n_z, freq_array)
  bar_image_2d=DBLARR(n_z, n_z)
  FOR K = 0, n_z-1 DO bar_image_2d[K,*]=bar_image_1d
  bar_signal_1d_0 = SM_0 ## bar_image_1d
  bar_signal_1d_1 = SM_1 ## bar_image_1d
  bar_signal_1d_2 = bar_signal_1d_0 - bar_signal_1d_1*parameters.dog_factor/10.
;  bar_signal_1d_2 = SM_2 ## bar_image_1d
  ;@GJ, 2024/6/4, adding the current signal ratio to reflect the harmonics
;  noise_level_percentage = 1. / (10.^(parameters.noise_level/20.)) * 100.
;  noise_level_percentage = 1. / (10.^(parameters.noise_level/10.)) * 100.
  noise_level_percentage = 1. / (10.^(parameters.noise_level/10.)*current_signal_ratio) * 100.
  print, 'noise level =  ', noise_level_percentage, '%'
  noisy_bar_signal_1d_0 = bar_signal_1d_0;AddPoissonNoise_FMMPI(bar_signal_1d_0, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(bar_signal_1d_0), LEVELS=[noise_level_percentage, noise_level_percentage]);/255.*MAX(bar_signal_1d_0); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
  noisy_bar_signal_1d_1 = bar_signal_1d_1;AddPoissonNoise_FMMPI(bar_signal_1d_1, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(bar_signal_1d_1), LEVELS=[noise_level_percentage, noise_level_percentage]);/255.*MAX(bar_signal_1d_1);AddPoissonNoise_FMMPI(image_sim_conv_1, PercentNoise=noise_level_percentage, Seed=seed)
  noisy_bar_signal_1d_2 = [[noisy_bar_signal_1d_0], [noisy_bar_signal_1d_1]]
;  noisy_bar_signal_1d_2 = NOISE_SCATTER(BYTSCL(bar_signal_1d_0), LEVELS=[0., noise_level_percentage])
  noisy_bar_signal_1d_2 = noisy_bar_signal_1d_0 - noisy_bar_signal_1d_1*parameters.dog_factor/10.
  
  sv_inv_0 = FLTARR(n_z, n_z)*0.
  bar_y_rec_0 = FLTARR(n_z, n_z)*0.
  correl_array_0 = DBLARR(n_z)*0.
  sv_inv_1 = FLTARR(n_z, n_z)*0.
  bar_y_rec_1 = FLTARR(n_z, n_z)*0.
  correl_array_1 = DBLARR(n_z)*0.
  sv_inv_2 = FLTARR(n_z, n_z)*0.
  bar_y_rec_2 = FLTARR(n_z, n_z)*0.
  correl_array_2 = DBLARR(n_z)*0.
  FOR K = 0, n_z-1 DO BEGIN
    sv_inv_0[K,K] = 1./w0[K]
    bar_y_0 = V0 ## sv_inv_0 ## TRANSPOSE(U0) ## noisy_bar_signal_1d_0
    bar_y_0 = bar_y_0/MAX(bar_y_0)*MAX(bar_image_1d)
    bar_y_rec_0[K, *] = REFORM(bar_y_0)
    correl_array_0[K] = CORRELATE(bar_y_0, bar_image_1d)
    
    sv_inv_1[K,K] = 1./w1[K]
    bar_y_1 = V1 ## sv_inv_1 ## TRANSPOSE(U1) ## noisy_bar_signal_1d_1
    bar_y_1 = bar_y_1/MAX(bar_y_1)*MAX(bar_image_1d)
    bar_y_rec_1[K, *] = REFORM(bar_y_1)
    correl_array_1[K] = CORRELATE(bar_y_1, bar_image_1d)
    
    sv_inv_2[K,K] = 1./w2[K]
    bar_y_2 = V2 ## sv_inv_2 ## TRANSPOSE(U2) ## noisy_bar_signal_1d_2
    bar_y_2 = bar_y_2/MAX(bar_y_2)*MAX(bar_image_1d)
    bar_y_rec_2[K, *] = REFORM(bar_y_2)
    correl_array_2[K] = CORRELATE(bar_y_2, bar_image_1d)
  ENDFOR
  
  ;@GJ, 2023/12/22, loading the point source phantom by default
  image_sim_ori = DBLARR(n_x, n_z) * 0.
  ;add the bar phantom
  ;  bar_image_1d=image_phantom_res(n_z, freq_array)
  ;  image_sim_ori[n_x/2-30:n_x/2+30, *] = bar_image_1d
  image_sim_ori = bar_image_2d
  ;  image_sim_ori[n_x/2-30:n_x/2+30, n_z/2-10:n_z/2-5] = 1.
  ;  image_sim_ori[n_x/2-30:n_x/2+30, n_z/2+5:n_z/2+10] = 1.
  rows_S = N_ELEMENTS(image_sim_ori[0,*])
  cols_S = N_ELEMENTS(image_sim_ori[*,0])

  ;@GJ, 2023/12/22, loading the phantom image
  IF STRLEN(parameters.filename) GT 5 THEN BEGIN
    IF (FILE_INFO(parameters.filename)).exists EQ 1 THEN BEGIN
      IF QUERY_IMAGE(parameters.filename) EQ 1 THEN BEGIN
        ;        image_temp = read_png(parameters.filename)
        image_temp = read_image(parameters.filename)
        ;        iimage, image_temp
        IF STRPOS(parameters.filename, 'gj') NE -1 THEN BEGIN
          IF (SIZE(image_temp))[0] EQ 2 THEN BEGIN
            image_sim_temp = DOUBLE(REFORM(image_temp[*,*]))
            ;      IF STRPOS(filename, 'MRA') NE -1 THEN image_sim_ori = DOUBLE(image_sim_temp GT 78) * 255.
            image_sim_ori = CONGRID(image_sim_temp, n_x, n_z, /INTERP)
            IF STRPOS(parameters.filename, 'MRA') NE -1 THEN image_sim_ori[WHERE(image_sim_ori LT 78)] = 0.
          ENDIF ELSE BEGIN
            image_sim_temp = DOUBLE(REFORM(image_temp[0,*,*]))
            image_sim_ori = CONGRID(MAX(image_sim_temp) - image_sim_temp, n_x, n_z)
            IF STRPOS(parameters.filename, 'MRA') NE -1 THEN image_sim_ori = DOUBLE(image_sim_ori LT 192) * 255.
          ENDELSE
        ENDIF ELSE BEGIN
          IF (SIZE(image_temp))[0] EQ 2 THEN BEGIN
            image_sim_temp = DOUBLE(REFORM(image_temp[*,*]))
          ENDIF ELSE BEGIN
            image_sim_temp = DOUBLE(REFORM(image_temp[0,*,*]))
          ENDELSE
          image_sim_ori = CONGRID(image_sim_temp, n_x, n_z, /INTERP)
        ENDELSE
      ENDIF
    ENDIF
  ENDIF
  
  IF parameters.FOV_current LT parameters.FOV THEN BEGIN
    factor_ima = (parameters.FOV_current) / (parameters.FOV)
    temp_n_x = FLOOR(n_x*factor_ima/2.)
    temp_n_z = FLOOR(n_z*factor_ima/2.)
    temp_image_sim = CONGRID(image_sim_ori, temp_n_x*2., temp_n_z*2., /INtERP)
    image_sim_ori *= 0.
    image_sim_ori[n_x/2-temp_n_x+1:n_x/2+temp_n_x, n_z/2-temp_n_z+1:n_z/2+temp_n_z] = temp_image_sim
  ENDIF  
  
  ;@GJ, 2023/1/13, calculate the max coorelation index
   correl_max_0 = MAX(correl_array_0, correl_max_ind_0)
   correl_max_1 = MAX(correl_array_1, correl_max_ind_1)
   correl_max_2 = MAX(correl_array_2, correl_max_ind_2)
;  iimage, bar_y_rec_0, title='Bright Reconstrcted Bar Phantom'
;  iimage, bar_y_rec_1, title='Donut Reconstrcted Bar Phantom'
;  iimage, bar_y_rec_2, title='Combined Reconstrcted Bar Phantom'
;  iplot, correl_array_0, title='bar phantom correlation'
;  
;@GJ, 2024/4/19, plot the correlatation array
;  IF wid GT 0 THEN BEGIN
;    ;@GJ, 2023/12/31, plotting the PSF
;    draw26 = widget_info(wid, find_by_uname='draw26')
;    widget_control, draw26, get_value = drawWindow26
;    wset, drawWindow26
;    ;    plot, x_index, w2, /nodata, title='SVD of System Matrix ('+STRING(parameters.A_DC, format='(i0)')+'%+'+STRING(parameters.A_DC_1, format='(i0)')+'%)', xtitle='Component Index', ytitle='Singular Value', /xlog, /ylog
;    plot, x_index, correl_array_2, /nodata, title='Correlation', yrange=[0., 1.], xtitle='# of SV', ytitle='Correlation Coeff';, /ylog
;    oplot, x_index, correl_array_2, thick=1, color='ffb037'x
;    oplot, x_index, correl_array_0, thick=1, color='0000ff'x
;    oplot, x_index, correl_array_1, thick=1, color='00ff00'x
;  ENDIF
;@GJ, 2024/4/19, plot the correlatation array, end
  
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 3
    draw2 = widget_info(wid, find_by_uname='draw2')
    widget_control, draw2, get_value = drawWindow2
    wset, drawWindow2
    TVSCL, BYTSCL(CONGRID(focal_image_0_display, drawXSize, drawYSize), MAX=MAX(focal_image_0_display[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.10,drawysize*0.9,'1st PSF (FOV='+STRING(parameters.FOV_current, format='(f5.1)')+' mm)', color = 'FFFFFF'x, /DEVICE
    IF donut_flag_s EQ 1 THEN label_title = 'Radius:' ELSE label_title = 'FWHM:'
    IF ABS(donut_radius_0/gradient_field_yz_current) LT 999. THEN BEGIN
      xyouts, drawysize*0.25,drawysize*0.8,label_title+STRING(donut_radius_0/gradient_field_yz_current, format='(f6.2)')+' mm', color = 'FFFFFF'x, /DEVICE
    ENDIF ELSE BEGIN
      xyouts, drawysize*0.25,drawysize*0.8,label_title+STRING(donut_radius_0/gradient_field_yz_current, format='(e10.2)')+' mm', color = 'FFFFFF'x, /DEVICE
    ENDELSE
    psf_filename=parameters.rec_directory+file_basename(parameters.filename)+'_gaussian_PSF.tif'
    IF (FILE_INFO(parameters.rec_directory)).exists EQ 0 THEN FILE_MKDIR, parameters.rec_directory
    write_tiff, psf_filename, BYTSCL(focal_image_0), /long

    LOADCT, 8
    draw12 = widget_info(wid, find_by_uname='draw12')
    widget_control, draw12, get_value = drawWindow12
    wset, drawWindow12
    TVSCL, BYTSCL(CONGRID(focal_image_1_display, drawXSize, drawYSize), MAX=MAX(focal_image_1_display[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.10,drawysize*0.9,'2nd PSF (FOV='+STRING(parameters.FOV_current, format='(f5.1)')+' mm)', color = 'FFFFFF'x, /DEVICE
    IF donut_flag_1 EQ 1 THEN label_title = 'Radius:' ELSE label_title = 'FWHM:'
    IF ABS(donut_radius_1/gradient_field_yz_current) LT 999. THEN BEGIN
      xyouts, drawysize*0.25,drawysize*0.8,label_title+STRING(donut_radius_1/gradient_field_yz_current, format='(f6.2)')+' mm', color = 'FFFFFF'x, /DEVICE
    ENDIF ELSE BEGIN
      xyouts, drawysize*0.25,drawysize*0.8,label_title+STRING(donut_radius_1/gradient_field_yz_current, format='(e10.2)')+' mm', color = 'FFFFFF'x, /DEVICE
    ENDELSE
    psf_filename=parameters.rec_directory+file_basename(parameters.filename)+'_donut_PSF.tif'
    write_tiff, psf_filename, BYTSCL(focal_image_1), /long
    
    ;subtraction
    LOADCT, 64
    draw22 = widget_info(wid, find_by_uname='draw22')
    widget_control, draw22, get_value = drawWindow22
    wset, drawWindow22
    TVSCL, 255-BYTSCL(CONGRID(focal_image_sub_display, drawXSize, drawYSize), MAX=MAX(focal_image_sub_display[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.25,drawysize*0.9,'PSF (subtracted)', color = '000000'x, /DEVICE
    IF donut_flag_s EQ 1 THEN label_title = 'Radius:' ELSE label_title = 'FWHM:'
    IF ABS(donut_radius_s/gradient_field_yz_current) LT 999. THEN BEGIN
      xyouts, drawysize*0.25,drawysize*0.8,label_title+STRING(donut_radius_s/gradient_field_yz_current, format='(f6.2)')+' mm', color = '000000'x, /DEVICE
    ENDIF ELSE BEGIN
      xyouts, drawysize*0.25,drawysize*0.8,label_title+STRING(donut_radius_s/gradient_field_yz_current, format='(e10.2)')+' mm', color = '000000'x, /DEVICE
    ENDELSE
    psf_filename=parameters.rec_directory+file_basename(parameters.filename)+'_STED_PSF.tif'
    write_tiff, psf_filename, BYTSCL(focal_image_sub), /long
    
    DEVICE, DECOMPOSED = 1
  ENDIF
  
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
;    LOADCT, 3
;    draw3 = widget_info(wid, find_by_uname='draw3')
;    widget_control, draw3, get_value = drawWindow3
;    wset, drawWindow3
;;    TVSCL, BYTSCL(CONGRID(image_sim_ori, drawXSize, drawYSize), MAX=MAX(image_sim_ori[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;;    xyouts,drawysize*0.35,drawysize*0.9,'Phantom', color = 'FFFFFF'x, /DEVICE
;    
;    LOADCT, 1
;    draw13 = widget_info(wid, find_by_uname='draw13')
;    widget_control, draw13, get_value = drawWindow13
;    wset, drawWindow13
;;    TVSCL, BYTSCL(CONGRID(image_sim_ori, drawXSize, drawYSize), MAX=MAX(image_sim_ori[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;;    xyouts,drawysize*0.35,drawysize*0.9,'Phantom', color = 'FFFFFF'x, /DEVICE
    
    LOADCT, 0;8
    topClr_1 = !D.TABLE_SIZE - 1
    TVLCT, 255, 255, 0, topClr_1
    draw23 = widget_info(wid, find_by_uname='draw23')
    widget_control, draw23, get_value = drawWindow23
    wset, drawWindow23
    temp_image_sim_ori = BYTSCL(CONGRID(image_sim_ori, drawXSize, drawYSize), MAX=MAX(image_sim_ori[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), TOP=(topClr_1-1))
    IF parameters.FOV_current LT parameters.FOV THEN BEGIN
      factor_ima = (parameters.FOV_current) / (parameters.FOV)
      temp_n_x = FLOOR(drawXSize*factor_ima/2.)
      temp_n_z = FLOOR(drawYSize*factor_ima/2.)
      temp_image_sim_ori[drawXSize/2-temp_n_x+1, drawYSize/2-temp_n_z+1:drawYSize/2+temp_n_z] = topClr_1
      temp_image_sim_ori[drawXSize/2+temp_n_x, drawYSize/2-temp_n_z+1:drawYSize/2+temp_n_z] = topClr_1
      temp_image_sim_ori[drawXSize/2-temp_n_x+1:drawXSize/2+temp_n_x, drawYSize/2-temp_n_z+1] = topClr_1
      temp_image_sim_ori[drawXSize/2-temp_n_x+1:drawXSize/2+temp_n_x, drawYSize/2+temp_n_z] = topClr_1
    ENDIF
    
    TVSCL, temp_image_sim_ori
    xyouts,drawysize*0.25,drawysize*0.9,'Phantom ('+STRING(parameters.FOV, format='(I0)')+'cm)', color = topClr_1-1, /DEVICE
    
    DEVICE, DECOMPOSED = 1
  ENDIF
  
  ;@GJ, 2023/12/22, do the image convolution
;  image_sim_conv = CONVOL(image_sim_ori, focal_image_0, /NORMALIZE, /CENTER, /EDGE_ZERO)
  IF STRPOS(parameters.filename, '.tif') NE -1 THEN BEGIN
    phantom_filename=parameters.rec_directory+file_basename(parameters.filename)
  ENDIF ELSE BEGIN
    phantom_filename=parameters.rec_directory+file_basename(parameters.filename)+'.tif'
  ENDELSE
  write_tiff, phantom_filename, BYTSCL(image_sim_ori), /long
  
  ;@GJ, 2024/8/25, do the redon rec simulation
  IF parameters.radon_rec EQ 1 THEN BEGIN
    ;@GJ, do FBP and Radon of the phantom
    ;;Calculate and display the Radon transform:
    Radon_result = RADON(image_sim_ori, RHO=rho, THETA=theta, NTHETA=256);, NRHO=256)
    delta_rho = 0.5 * SQRT(delta_x^2 + delta_z^2)

    Radon_result_conv_0_real = Radon_result*0.
    Radon_result_conv_0_imag = Radon_result*0.
    Radon_result_conv_0_amp = Radon_result*0.
    Radon_result_conv_0_phase = Radon_result*0.
    Radon_result_conv_1_real = Radon_result*0.
    Radon_result_conv_1_imag = Radon_result*0.
    Radon_result_conv_1_amp = Radon_result*0.
    Radon_result_conv_1_phase = Radon_result*0.
    n_theta = N_ELEMENTS(theta)
    n_rho = N_ELEMENTS(rho)
    n_slice = N_ELEMENTS(rho)
    FOR it = 0, n_theta-1 DO BEGIN
      ;@GJ, 2024/8/21, do the Radon transform based on 2D image
      Radon_result_2d = DBLARR(n_rho, n_slice) * 0.
      FOR is = 0, 11 DO Radon_result_2d[*, n_slice/2-5+is] = REFORM(Radon_result[it, *])
      ;@GJ, 2024/8/16, do the real and imag conv and then calculate the phase
      Radon_result_conv_0_real[it, *] = (CONVOL_FFT(Radon_result_2d, CONGRID(focal_image_0_real, n_x*SQRT(2.), n_z*SQRT(2.))))[*, n_slice/2]
      Radon_result_conv_0_imag[it, *] = (CONVOL_FFT(Radon_result_2d, CONGRID(focal_image_0_imag, n_x*SQRT(2.), n_z*SQRT(2.))))[*, n_slice/2]
      Radon_result_conv_1_real[it, *] = (CONVOL_FFT(Radon_result_2d, CONGRID(focal_image_1_real, n_x*SQRT(2.), n_z*SQRT(2.))))[*, n_slice/2]
      Radon_result_conv_1_imag[it, *] = (CONVOL_FFT(Radon_result_2d, CONGRID(focal_image_1_imag, n_x*SQRT(2.), n_z*SQRT(2.))))[*, n_slice/2]
      print, '(convol) i_theta:', it, ' of n_theta', n_theta
    ENDFOR
    Radon_result_conv_0_amp = SQRT(Radon_result_conv_0_real^2 + Radon_result_conv_0_imag^2)
    Radon_result_conv_0_phase = 180./!PI*ATAN(Radon_result_conv_0_imag, Radon_result_conv_0_real)
    Radon_result_conv_1_amp = SQRT(Radon_result_conv_1_real^2 + Radon_result_conv_1_imag^2)
    Radon_result_conv_1_phase = 180./!PI*ATAN(Radon_result_conv_1_imag, Radon_result_conv_1_real)
    ;subtraction
    Radon_result_conv_real_sub = Radon_result_conv_0_real -  Radon_result_conv_1_real * parameters.dog_factor/10.
    Radon_result_conv_imag_sub = Radon_result_conv_0_imag -  Radon_result_conv_1_imag * parameters.dog_factor/10.
    Radon_result_conv_amp_sub = Radon_result_conv_0_amp -  Radon_result_conv_1_amp * parameters.dog_factor/10.
    Radon_result_conv_phase_sub = Radon_result_conv_0_phase -  Radon_result_conv_1_phase * parameters.dog_factor/10.

    ;@GJ, 2024/8/21, display the signals
    IIMAGE, BYTSCL(CONGRID(Radon_result, 850/3, 550)), VIEW_GRID=[3,1], VIEW_TITLE='Radon', $
      DIMENSIONS=[850, 550], WINDOW_TITLE='Radon', $
      /NO_SAVEPROMPT
    iimage, BYTSCL(CONGRID(Radon_result_conv_amp_sub, 850/3, 550)), /VIEW_NEXT, VIEW_TITLE='Radon amp'
    iimage, BYTSCL(CONGRID(Radon_result_conv_phase_sub, 850/3, 550)), /VIEW_NEXT, VIEW_TITLE='Radon phase'

    ;@GJ, 2024/8/21, do the filtered BP
    Radon_result_filtered = Radon_result*0.
    FOR it = 0, n_theta-1 DO BEGIN
      ;get the radon line
      Radon_result_1d = REFORM(Radon_result[it, *])
      Radon_result_filtered[it, *] = Filter_FFT(Radon_result_1d)

      ;get the filtered for DC_0
      Radon_result_conv_real_sub[it, *] = Filter_FFT(REFORM(Radon_result_conv_real_sub[it, *]))
      Radon_result_conv_imag_sub[it, *] = Filter_FFT(REFORM(Radon_result_conv_imag_sub[it, *]))
      Radon_result_conv_amp_sub[it, *] = Filter_FFT(REFORM(Radon_result_conv_amp_sub[it, *]))
      Radon_result_conv_phase_sub[it, *] = Filter_FFT(REFORM(Radon_result_conv_phase_sub[it, *]))
      ;    print, '(filter) i_theta:', it, ' of n_theta', n_theta
    ENDFOR

    ;@GJ, 2024/8/20, plot the FBP results
    backproject_result = RADON(Radon_result_filtered, /BACKPROJECT, RHO=rho, THETA=theta)
    backproject_sub_amp = RADON(Radon_result_conv_amp_sub, /BACKPROJECT, RHO=rho, THETA=theta)
    backproject_sub_phase = RADON(Radon_result_conv_phase_sub, /BACKPROJECT, RHO=rho, THETA=theta)

    IIMAGE, BYTSCL(CONGRID(backproject_result, 850/3, 550)), VIEW_GRID=[3,1], VIEW_TITLE='FBP Radon', $
      DIMENSIONS=[850, 550], WINDOW_TITLE='FBP', $
      /NO_SAVEPROMPT
    iimage, BYTSCL(CONGRID(backproject_sub_amp, 850/3, 550)), /VIEW_NEXT, VIEW_TITLE='FBP amp'
    iimage, BYTSCL(CONGRID(backproject_sub_phase, 850/3, 550)), /VIEW_NEXT, VIEW_TITLE='FBP phase'

    ;@GJ, 2024/8/20, delete variables to release memory
    DELVAR, Radon_result_conv_0_real, Radon_result_conv_0_imag, Radon_result_conv_0_amp, Radon_result_conv_0_phase
    DELVAR, Radon_result_conv_1_real, Radon_result_conv_1_imag, Radon_result_conv_1_amp, Radon_result_conv_1_phase
    DELVAR, Radon_result, Radon_result_conv_real_sub, Radon_result_conv_imag_sub, Radon_result_conv_amp_sub, Radon_result_conv_phase_sub
    DELVAR, backproject_sub_amp, backproject_sub_phase, backproject_result, Radon_result_filtered
  ENDIF
  
  ;@GJ, 2024/8/16, do the real and imag conv and then calculate the phase
  image_sim_conv_real = CONVOL_FFT(image_sim_ori, focal_image_0_real)
;  image_sim_conv_real = CONVOL(image_sim_ori, focal_image_0_real, /EDGE_ZERO, /NORMALIZE)
  image_sim_conv_imag = CONVOL_FFT(image_sim_ori, focal_image_0_imag);
;  image_sim_conv_imag = CONVOL(image_sim_ori, focal_image_0_imag, /EDGE_ZERO, /NORMALIZE)
  IF parameters.phase_scatter_plot EQ 1 THEN BEGIN
    ;@GJ, 2024/8/17, plot the real and imag as scatter
    ;@GJ, 2024/8/17, plot the real and imag
    iplot, focal_line_real, focal_line_imag, xtitle='Real', ytitle='Imag', SYM_INDEX=1, LINESTYLE=6, COLOR='red', title='Scatter Plot DC1 & 2', /NO_SAVEPROMPT
;    iplot, focal_image_0_real, focal_image_0_imag, SYM_INDEX=4, LINESTYLE=6, COLOR='pink', /OVERPLOT
;    parameters.shift_angle = 30.
    shift_0_real = -focal_line_imag * sin(-parameters.shift_angle/180.*!PI) + focal_line_real * cos(-parameters.shift_angle/180.*!PI)
    shift_0_imag = focal_line_imag * cos(-parameters.shift_angle/180.*!PI) + focal_line_real * sin(-parameters.shift_angle/180.*!PI)
;    iplot, shift_0_real, shift_0_imag, SYM_INDEX=5, LINESTYLE=6, COLOR='cyan', /OVERPLOT
;    iplot, image_sim_conv_real, image_sim_conv_imag, SYM_INDEX=5, LINESTYLE=6, COLOR='purple', /OVERPLOT
  ENDIF
  ;@GJ, 2024/8/16, real or imag
  image_sim_conv = CONVOL_FFT(image_sim_ori, focal_image_0)
  IF harmonics_comp EQ 0 THEN BEGIN
    image_sim_conv = SQRT(image_sim_conv_real^2 + image_sim_conv_imag^2)
  ENDIF
  IF harmonics_comp EQ 1 THEN BEGIN
    image_sim_conv = 180./!PI*ATAN(image_sim_conv_imag, image_sim_conv_real)
  ENDIF
  
  image_sim_conv_1_real = CONVOL_FFT(image_sim_ori, focal_image_1_real)
  image_sim_conv_1_imag = CONVOL_FFT(image_sim_ori, focal_image_1_imag)
  ;@GJ, 2024/8/17, plot the real and imag as scatter
  IF parameters.phase_scatter_plot EQ 1 THEN BEGIN
    iplot, focal_line_1_real, focal_line_1_imag, SYM_INDEX=1, LINESTYLE=6, COLOR='blue', /OVERPLOT
    shift_1_real = -focal_line_1_imag * sin(-parameters.shift_angle/180.*!PI) + focal_line_1_real * cos(-parameters.shift_angle/180.*!PI)
    shift_1_imag = focal_line_1_imag * cos(-parameters.shift_angle/180.*!PI) + focal_line_1_real * sin(-parameters.shift_angle/180.*!PI)
;    iplot, shift_1_real, shift_1_imag, SYM_INDEX=5, LINESTYLE=6, COLOR='cyan', /OVERPLOT
;   iplot, focal_image_1_real, focal_image_1_imag, SYM_INDEX=1, LINESTYLE=6, COLOR='green', /OVERPLOT
 ;   iplot, image_sim_conv_1_real, image_sim_conv_1_imag, SYM_INDEX=3, LINESTYLE=6, COLOR='green', /OVERPLOT
  ENDIF
  ;@GJ, 2024/8/16, real or imag
  image_sim_conv_1 = CONVOL_FFT(image_sim_ori, focal_image_1)
  IF harmonics_comp EQ 0 THEN BEGIN
    image_sim_conv_1 = SQRT(image_sim_conv_1_real^2 + image_sim_conv_1_imag^2)
  ENDIF
  IF harmonics_comp EQ 1 THEN BEGIN
    image_sim_conv_1 = 180./!PI*ATAN(image_sim_conv_1_imag, image_sim_conv_1_real)
  ENDIF

  ;@GJ, 2024/8/20, delete the variables to release memory
  DELVAR, image_sim_conv_real, image_sim_conv_imag
  DELVAR, image_sim_conv_1_real, image_sim_conv_1_imag

  image_sim_conv_sub = image_sim_conv -  image_sim_conv_1 * parameters.dog_factor/10.
  IF wid GT 0  AND parameters.particle_type GT 2 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 3
    draw5 = widget_info(wid, find_by_uname='draw5')
    widget_control, draw5, get_value = drawWindow5
    wset, drawWindow5
    TVSCL, BYTSCL(CONGRID(image_sim_conv, drawXSize, drawYSize), MAX=MAX(image_sim_conv[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.15,drawysize*0.9,'1st Convoluted Image', color = 'FFFFFF'x, /DEVICE
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(image_sim_conv));72.1795
    ;    print, 'Rec 1 MSE: ', MSE
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
    ;    print, 'Rec 1 pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(image_sim_conv))
    ;    print, 'Rec 1 SSIM: ', mssim
;    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
        
    LOADCT, 8
    draw15 = widget_info(wid, find_by_uname='draw15')
    widget_control, draw15, get_value = drawWindow15
    wset, drawWindow15
    TVSCL, BYTSCL(CONGRID(image_sim_conv_1, drawXSize, drawYSize), MAX=MAX(image_sim_conv_1[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.15,drawysize*0.9,'2nd Convoluted Image', color = 'FFFFFF'x, /DEVICE
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(image_sim_conv_1));72.1795
    ;    print, 'Rec 1 MSE: ', MSE
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
    ;    print, 'Rec 1 pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(image_sim_conv_1))
    ;    print, 'Rec 1 SSIM: ', mssim
;    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
    
    ;subtraction
    LOADCT, 64
    draw25 = widget_info(wid, find_by_uname='draw25')
    widget_control, draw25, get_value = drawWindow25
    wset, drawWindow25
    TVSCL, 255-BYTSCL(CONGRID(image_sim_conv_sub, drawXSize, drawYSize), MAX=MAX(image_sim_conv_sub[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.10,drawysize*0.9,'Convoluted (subtracted)', color = '000000'x, /DEVICE
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(image_sim_conv_sub));72.1795
;    print, 'MSE: ', MSE
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
;    print, 'pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(image_sim_conv_sub))
;    print, 'SSIM: ', mssim
;    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = '000000'x, /DEVICE

    DEVICE, DECOMPOSED = 1
  ENDIF ELSE BEGIN
  ;@GJ, 2023/1/14, display the kernel and chebyshev
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 0
    draw5 = widget_info(wid, find_by_uname='draw5')
    widget_control, draw5, get_value = drawWindow5
    wset, drawWindow5
    temp_kernel_image = REFORM((*parameters.Sx_array_xzNonCOnv)[harmonics_n-1, harmonics_comp, *, *])
    temp_image = ROTATE(CONGRID(temp_kernel_image, drawXSize, drawYSize),1)
;    IF gradient_field_x LT gradient_field_yz THEN BEGIN
;      width_ratio = (n_x*gradient_field_x)/(n_z*gradient_field_yz)
;      temp_image = ROTATE(CONGRID(temp_kernel_image[(n_x/2.-0.5*width_ratio*n_x):n_x/2.+0.5*width_ratio*n_x, *], drawXSize, drawYSize),1)
;    ENDIF ELSE BEGIN
;      width_ratio = (n_z*gradient_field_yz)/(n_x*gradient_field_x)
;      temp_image = ROTATE(CONGRID(temp_kernel_image[*, (n_z/2.-0.5*width_ratio*n_z):n_z/2.+0.5*width_ratio*n_z], drawXSize, drawYSize),1)
;    ENDELSE
    TVSCL, temp_image;BYTSCL(CONGRID(temp_kernel_image, drawXSize, drawYSize), MAX=MAX(temp_kernel_image[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(temp_kernel_image[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.40,drawysize*0.9,'Base', color = 'FFFFFF'x, /DEVICE

    LOADCT, 0
    draw15 = widget_info(wid, find_by_uname='draw15')
    widget_control, draw15, get_value = drawWindow15
    wset, drawWindow15
    temp_Chebyshev_image = REFORM((*parameters.Sx_array_kernel)[harmonics_n-1, harmonics_comp, *, *])
    temp_image = ROTATE(CONGRID(temp_Chebyshev_image, drawXSize, drawYSize), 1)
;    IF gradient_field_x LT gradient_field_yz THEN BEGIN
;      width_ratio = (n_x*gradient_field_x)/(n_z*gradient_field_yz)
;      temp_image = ROTATE(CONGRID(temp_Chebyshev_image[(n_x/2.-0.5*width_ratio*n_x):n_x/2.+0.5*width_ratio*n_x, *], drawXSize, drawYSize),1)
;    ENDIF ELSE BEGIN
;      width_ratio = (n_z*gradient_field_yz)/(n_x*gradient_field_x)
;      temp_image = ROTATE(CONGRID(temp_Chebyshev_image[*, (n_z/2.-0.5*width_ratio*n_z):n_z/2.+0.5*width_ratio*n_z], drawXSize, drawYSize),1)
;    ENDELSE
    TVSCL, temp_image;BYTSCL(CONGRID(temp_Chebyshev_image, drawXSize, drawYSize), MAX=MAX(temp_Chebyshev_image[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(temp_Chebyshev_image[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.36,drawysize*0.9,'Kernel', color = 'FFFFFF'x, /DEVICE

    DEVICE, DECOMPOSED = 1
  ENDELSE
  
  ;@GJ, 2023/12/23, adding Poission noise
  seed = -3L ; Only so that the example here matches your result!
;  Print, "Mean and standard deviation of original image is: " $
;    + StrTrim(Mean(image_sim_conv),2) + ", " + StrTrim(StdDev(image_sim_conv),2)
  ;@GJ, 2023/12/23, converting noise level dB to percent Noise;@GJ, 2024/8/16, add the noisy
;  noisy_image_sim_conv_real = AddPoissonNoise_FMMPI(image_sim_conv_real, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(image_sim_conv_real), LEVELS=[noise_level_percentage/100.,noise_level_percentage/100.])*MEAN(image_sim_conv_real); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
;  noisy_image_sim_conv_imag = AddPoissonNoise_FMMPI(image_sim_conv_imag, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(image_sim_conv_imag), LEVELS=[noise_level_percentage/100.,noise_level_percentage/100.])*MEAN(image_sim_conv_imag); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
  noisy_image_sim_conv = image_sim_conv;AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(image_sim_conv), LEVELS=[noise_level_percentage/100.,noise_level_percentage/100.])*MEAN(image_sim_conv); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
;  IF harmonics_comp EQ 0 THEN BEGIN
;    noisy_image_sim_conv = AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed);SQRT(noisy_image_sim_conv_real^2 + noisy_image_sim_conv_imag^2)
;  ENDIF
;  IF harmonics_comp EQ 1 THEN BEGIN
;    noisy_image_sim_conv = AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed);180./!PI*ATAN(noisy_image_sim_conv_imag, noisy_image_sim_conv_real)
;  ENDIF
  ;@GJ, 2024/8/16, add the noisy
;  noisy_image_sim_conv_1_real = AddPoissonNoise_FMMPI(image_sim_conv_1_real, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(image_sim_conv_1_real), LEVELS=[noise_level_percentage/100.,noise_level_percentage/100.])*MEAN(image_sim_conv_1_real); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
;  noisy_image_sim_conv_1_imag = AddPoissonNoise_FMMPI(image_sim_conv_1_imag, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(image_sim_conv_1_imag), LEVELS=[noise_level_percentage/100.,noise_level_percentage/100.])*MEAN(image_sim_conv_1_imag); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
  noisy_image_sim_conv_1 = image_sim_conv_1;AddPoissonNoise_FMMPI(image_sim_conv_1, PercentNoise=noise_level_percentage, Seed=seed);NOISE_SCATTER(BYTSCL(image_sim_conv_1), LEVELS=[noise_level_percentage/100.,noise_level_percentage/100.])*MEAN(image_sim_conv_1); AddPoissonNoise_FMMPI(image_sim_conv, PercentNoise=noise_level_percentage, Seed=seed)
;  IF harmonics_comp EQ 0 THEN BEGIN
;    noisy_image_sim_conv_1 = AddPoissonNoise_FMMPI(image_sim_conv_1, PercentNoise=noise_level_percentage, Seed=seed);SQRT(noisy_image_sim_conv_1_real^2 + noisy_image_sim_conv_1_imag^2)
;  ENDIF
;  IF harmonics_comp EQ 1 THEN BEGIN
;    noisy_image_sim_conv_1 = AddPoissonNoise_FMMPI(image_sim_conv_1, PercentNoise=noise_level_percentage, Seed=seed);180./!PI*ATAN(noisy_image_sim_conv_1_imag, noisy_image_sim_conv_1_real)
;  ENDIF
  ;@GJ, 2024/8/16, do DOG subtraction
  noisy_image_sim_conv_sub = noisy_image_sim_conv - noisy_image_sim_conv_1 * parameters.dog_factor/10.
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 3
    draw5 = widget_info(wid, find_by_uname='draw5')
    widget_control, draw5, get_value = drawWindow5
    wset, drawWindow5
    TVSCL, BYTSCL(CONGRID(noisy_image_sim_conv, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.30,drawysize*0.9,'1st Noisy Image', color = 'FFFFFF'x, /DEVICE
    
    noisy_image_sim_conv_FFT = ALOG(ABS(FFT(noisy_image_sim_conv, /CENTER)))
    draw6 = widget_info(wid, find_by_uname='draw6')
    widget_control, draw6, get_value = drawWindow6
    wset, drawWindow6
    temp_display = BYTSCL(CONGRID(noisy_image_sim_conv_FFT, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv_FFT[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    temp_display[MTF_indices_1st] = 255
    TVSCL, temp_display
    xyouts,drawysize*0.30,drawysize*0.9,'1st Noisy logFFT', color = 'FFFFFF'x, /DEVICE

    
    noisy_filename=parameters.rec_directory+file_basename(parameters.filename)+'_gaussian_noisy.tif'
    write_tiff, noisy_filename, BYTSCL(noisy_image_sim_conv), /long
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv));72.1795
;    print, 'Rec 1 MSE: ', MSE
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
;    print, 'Rec 1 pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv))
;    print, 'Rec 1 SSIM: ', mssim
;    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
    
    LOADCT, 8
    draw15 = widget_info(wid, find_by_uname='draw15')
    widget_control, draw15, get_value = drawWindow15
    wset, drawWindow15
    TVSCL, BYTSCL(CONGRID(noisy_image_sim_conv_1, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv_1[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.30,drawysize*0.9,'2nd Noisy Image', color = 'FFFFFF'x, /DEVICE

    noisy_image_sim_conv_1_FFT = ALOG(ABS(FFT(noisy_image_sim_conv_1, /CENTER)))
    draw16 = widget_info(wid, find_by_uname='draw16')
    widget_control, draw16, get_value = drawWindow16
    wset, drawWindow16
    TVSCL, BYTSCL(CONGRID(noisy_image_sim_conv_1_FFT, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv_1_FFT[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.30,drawysize*0.9,'2nd Noisy logFFT', color = 'FFFFFF'x, /DEVICE

    noisy_filename=parameters.rec_directory+file_basename(parameters.filename)+'_donut_noisy.tif'
    write_tiff, noisy_filename, BYTSCL(noisy_image_sim_conv_1), /long
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv_1));72.1795
;    print, 'Rec 1 MSE: ', MSE
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
;    print, 'Rec 1 pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv_1))
;    print, 'Rec 1 SSIM: ', mssim
;    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
    
    ;subtraction
    LOADCT, 64
    draw25 = widget_info(wid, find_by_uname='draw25')
    widget_control, draw25, get_value = drawWindow25
    wset, drawWindow25
    TVSCL, 255-BYTSCL(CONGRID(noisy_image_sim_conv_sub, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv_sub[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.10,drawysize*0.9,'Noisy Image (subtracted)', color = '000000'x, /DEVICE

    noisy_image_sim_conv_sub_FFT = ALOG(ABS(FFT(noisy_image_sim_conv_sub, /CENTER)))
    draw26 = widget_info(wid, find_by_uname='draw26')
    widget_control, draw26, get_value = drawWindow26
    wset, drawWindow26
    temp_display = BYTSCL(CONGRID(noisy_image_sim_conv_sub_FFT, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv_sub_FFT[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    temp_display[MTF_indices_sub] = 0
    TVSCL, temp_display
    xyouts,drawysize*0.20,drawysize*0.9,'Subtracted Noisy logFFT', color = 'FFFFFF'x, /DEVICE
 
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv_sub));72.1795
    noisy_filename=parameters.rec_directory+file_basename(parameters.filename)+'_STED_noisy.tif'
    write_tiff, noisy_filename, BYTSCL(noisy_image_sim_conv_sub), /long
    print, 'noise MSE: ', MSE
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
    print, 'noise pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv_sub))
    print, 'noise SSIM: ', mssim
;   xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = '000000'x, /DEVICE

;    LOADCT, 8
;    draw25 = widget_info(wid, find_by_uname='draw25')
;    widget_control, draw25, get_value = drawWindow25
;    wset, drawWindow25
;    TVSCL, BYTSCL(CONGRID(noisy_image_sim_conv_sub, drawXSize, drawYSize), MAX=MAX(noisy_image_sim_conv_sub[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
;    xyouts,drawysize*0.10,drawysize*0.9,'Noisy Image (subtracted)', color = 'FFFFFF'x, /DEVICE
;    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv_sub));72.1795
;    print, 'noise MSE: ', MSE
;    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
;    print, 'noise pSNR: ', pSNR, ' dB'
;    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(noisy_image_sim_conv_sub))
;    print, 'noise SSIM: ', mssim
;    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE    
    
    DEVICE, DECOMPOSED = 1
  ENDIF
  
  ;@GJ, 2023/12/23, deconvolution
  ;@GJ, 2023/12/31, use un-noised image
;  noisy_image_sim_conv = image_sim_conv
  focal_image_0_FFT = ABS(FFT(focal_image_0,/center,/double))
  noisy_image_sim_conv_FFT = FFT(noisy_image_sim_conv,/center,/double)
  noisy_image_FFT_div = noisy_image_sim_conv_FFT/focal_image_0_FFT
  noisy_image_FFT_div[where(focal_image_0_FFT LT 0.001*MAX(focal_image_0_FFT))] = COMPLEX(0, 0)
;  iimage, BYTSCL(ALOG(noisy_image_FFT_div)), title='1st Noisy FFT Div'
  rec_noisy_image = ABS(FFT(noisy_image_FFT_div, /INVERSE,/center,/double))
  
  ;GJ, 2024/6/7, keyhole method for precise image reconstruction 
;  keyhole_method, wid, parameters, index_A_DC, n_x, n_z, gradient_ratio, image_sim_ori
 
  ;@GJ, 2023/12/31, use un-noised image
;  noisy_image_sim_conv_1 = image_sim_conv_1
  focal_image_1_FFT = ABS(FFT(focal_image_1,/center,/double))
  noisy_image_1_sim_conv_FFT = FFT(noisy_image_sim_conv_1,/center,/double)
  noisy_image_1_FFT_div = noisy_image_1_sim_conv_FFT/focal_image_1_FFT
  noisy_image_1_FFT_div[where(focal_image_1_FFT LT 0.001*MAX(focal_image_1_FFT))] = COMPLEX(0, 0)
  rec_noisy_image_1 = ABS(FFT(noisy_image_1_FFT_div, /INVERSE,/center,/double))

  ;@GJ, 2023/12/31, use un-noised image
;  noisy_image_sim_conv_sub = image_sim_conv_sub
  focal_image_sub_FFT = ABS(FFT(focal_image_sub,/center,/double))
  noisy_image_sub_sim_conv_FFT = FFT(noisy_image_sim_conv_sub,/center,/double)
  noisy_image_sub_FFT_div = noisy_image_sub_sim_conv_FFT/focal_image_sub_FFT
  noisy_image_sub_FFT_div[where(focal_image_sub_FFT LT 0.001*MAX(focal_image_sub_FFT))] = COMPLEX(0, 0)
  rec_noisy_image_sub = ABS(FFT(noisy_image_sub_FFT_div, /INVERSE,/center,/double))
  *parameters.rec_noisy_image_sub = rec_noisy_image_sub
  *parameters.image_sim_ori = image_sim_ori
  ;@GJ, 2024/4/20, adding a threshold to remove the background
;  temp_result = IMAGE_THRESHOLD(BYTSCL(rec_noisy_image), THRESHOLD=tf, /MAXENTROPY)
;  rec_noisy_image=BYTSCL(rec_noisy_image, MIN=tf/2.);[WHERE(BYTSCL(rec_noisy_image) LT tf/2.)] = MIN(rec_noisy_image)
;  temp_result_1 = IMAGE_THRESHOLD(BYTSCL(rec_noisy_image_1), THRESHOLD=tf_1, /MAXENTROPY)
;  rec_noisy_image_1=BYTSCL(rec_noisy_image_1, MIN=tf_1/2.);rec_noisy_image_1[WHERE(BYTSCL(rec_noisy_image_1) LT tf_1/2.)] = MIN(rec_noisy_image_1)
  ;Rec_sub_threshold_optimization, wid, parameters
  temp_image_sub = BYTSCL(rec_noisy_image_sub)
  temp_image_sub[WHERE(temp_image_sub LT parameters.phantom_rec_sub_threshold, /null)] = 0.
;  temp_result_sub = IMAGE_THRESHOLD(temp_image, THRESHOLD=tf_sub, /MAXENTROPY, /INVERT)
;  temp_image[WHERE(temp_image LT tf_sub[0]/16., /null)] = 0.
  print, 'threshold:', parameters.phantom_rec_sub_threshold
  rec_noisy_image_sub = DOUBLE(temp_image_sub);temp_image(WHERE( temp_image LT tf_sub/3.)];rec_noisy_image_sub[WHERE(BYTSCL(rec_noisy_image_sub) LT tf_sub/2.)] = MIN(rec_noisy_image_sub)
  
  
;  rec_noisy_image = (IMAGE_THRESHOLD(BYTSCL(rec_noisy_image), THRESHOLD=tf, /MAXENTROPY)*rec_noisy_image)
;  rec_noisy_image_1 = (IMAGE_THRESHOLD(BYTSCL(rec_noisy_image_1), THRESHOLD=tf, /MAXENTROPY)*rec_noisy_image_1)
;  rec_noisy_image_sub = (IMAGE_THRESHOLD(BYTSCL(rec_noisy_image_sub), THRESHOLD=tf, /MAXENTROPY)*rec_noisy_image_sub)
;  rec_noisy_image_sub = rec_noisy_image - rec_noisy_image_1 * focal_ratio

;  rec_noisy_image = bar_y_rec_0
;  IF correl_max_ind_0-1 GE 0 AND correl_max_ind_0+1 LT n_x-1 THEN rec_noisy_image[correl_max_ind_0-1:correl_max_ind_0+1, *] *= 3.;MAX(rec_noisy_image)
;  rec_noisy_image_1 = bar_y_rec_1
;  IF correl_max_ind_1-1 GE 0 AND correl_max_ind_1+1 LT n_x-1 THEN rec_noisy_image_1[correl_max_ind_1-1:correl_max_ind_1+1, *] *= 3.;MAX(rec_noisy_image_1)
;  rec_noisy_image_sub =  bar_y_rec_2
;  IF correl_max_ind_2-1 GE 0 AND correl_max_ind_2+1 LT n_x-1 THEN rec_noisy_image_sub[correl_max_ind_2-1:correl_max_ind_2+1, *] *= 3.;MAX(rec_noisy_image_sub)
  IF wid GT 0 THEN BEGIN
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 3
    draw7 = widget_info(wid, find_by_uname='draw7')
    widget_control, draw7, get_value = drawWindow7
    wset, drawWindow7
    TVSCL, BYTSCL(CONGRID(rec_noisy_image, drawXSize, drawYSize), MAX=MAX(rec_noisy_image[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(rec_noisy_image[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.15,drawysize*0.9,'Reconstructed (Regular)', color = 'FFFFFF'x, /DEVICE
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(rec_noisy_image));72.1795
    print, 'Rec 0 MSE: ', MSE
    print, 'Rec 0 RMSE: ', STRING(SQRT(MSE), format='(f6.2)')
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
    print, 'Rec 0 pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(rec_noisy_image))
    print, 'Rec 0 SSIM: ', mssim
    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
;    xyouts, drawysize*0.15,drawysize*0.8,'Max Corr Coef:'+STRING(correl_max_0, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE

    
    LOADCT, 8
    draw17 = widget_info(wid, find_by_uname='draw17')
    widget_control, draw17, get_value = drawWindow17
    wset, drawWindow17
    TVSCL, BYTSCL(CONGRID(rec_noisy_image_1, drawXSize, drawYSize), MAX=MAX(rec_noisy_image_1[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(rec_noisy_image_1[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.15,drawysize*0.9,'Reconstructed (Donut)', color = 'FFFFFF'x, /DEVICE
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(rec_noisy_image_1));72.1795
    print, 'Rec 1 MSE: ', MSE
    print, 'Rec 1 RMSE: ', STRING(SQRT(MSE), format='(f6.2)')
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
    print, 'Rec 1 pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(rec_noisy_image_1))
    print, 'Rec 1 SSIM: ', mssim
    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
;    xyouts, drawysize*0.15,drawysize*0.8,'Max Corr Coef:'+STRING(correl_max_1, format='(f6.3)'), color = 'FFFFFF'x, /DEVICE
    
    LOADCT, 64
    draw27 = widget_info(wid, find_by_uname='draw27')
    widget_control, draw27, get_value = drawWindow27
    wset, drawWindow27
    TVSCL, 255-BYTSCL(CONGRID(rec_noisy_image_sub, drawXSize, drawYSize), MAX=MAX(rec_noisy_image_sub[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]), MIN=MIN(rec_noisy_image_sub[cols_S*0.1:cols_S*0.8, rows_S*0.1:rows_S*0.8]))
    xyouts,drawysize*0.15,drawysize*0.9,'Reconstructed (Subtr)', color = '000000'x, /DEVICE
    MSE=calMSE(BYTSCL(image_sim_ori), BYTSCL(rec_noisy_image_sub));72.1795
    print, 'Rec sub MSE: ', MSE
    print, 'Rec sub RMSE: ', STRING(SQRT(MSE), format='(f6.2)')
    pSNR=calPSNR(BYTSCL(image_sim_ori), MSE);68.0337
    print, 'Rec sub pSNR: ', pSNR, ' dB'
    mssim = SSIM(BYTSCL(image_sim_ori), BYTSCL(rec_noisy_image_sub))
    print, 'Rec sub SSIM: ', mssim
    xyouts, drawysize*0.02,drawysize*0.8,'PSNR: '+STRING(pSNR, format='(f6.2)')+', SSIM:'+STRING(mSSIM, format='(f6.3)'), color = '000000'x, /DEVICE    
;    xyouts, drawysize*0.15,drawysize*0.8,'Max Corr Coef:'+STRING(correl_max_2, format='(f6.3)'), color = '000000'x, /DEVICE
    
    DEVICE, DECOMPOSED = 1
  ENDIF
  ; Clean up.
  !P.Multi = 0
  
  ;@GJ, 2024/5/3, write the results into log
  write_results_FM_MPI, parameters
  
  plot_End_spot: print, 'no plot'
  
end

pro FMMPI_slider_A_DC_event, ev
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.AC_DC_mode EQ 0. THEN BEGIN
    ;AC excitation mode
    widget_control, ev.id, get_value = AC_ex
    parameters.AC_ex = AC_ex/100.;
  ENDIF ELSE BEGIN
    ;DC offset mode
    widget_control, ev.id, get_value = A_DC
    parameters.A_DC = A_DC;
  ENDELSE
  
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_A_DC_1_event, ev
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.AC_DC_mode EQ 0. THEN BEGIN
    ;AC excitation mode
    widget_control, ev.id, get_value = AC_ex_1
    parameters.AC_ex_1 = AC_ex_1/100.;
  ENDIF ELSE BEGIN
    ;DC offset mode
    widget_control, ev.id, get_value = A_DC_1
    parameters.A_DC_1 = A_DC_1;
  ENDELSE
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

;@GJ, 2024/4/19, calculate the automated dog factor
pro FMMPI_button_Dog_factor_auto_event, ev
  widget_control, ev.top, get_uvalue = parameters
  
  Dog_factor_optimization, ev.top, parameters
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_Dog_factor_event, ev
  widget_control, ev.id, get_value = Dog_factor
  widget_control, ev.top, get_uvalue = parameters
  parameters.Dog_factor = Dog_factor;
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_A_sin_event, ev
  widget_control, ev.id, get_value = A
  widget_control, ev.top, get_uvalue = parameters
  parameters.A = A;
  parameters.gradient_field_x = parameters.A * 2. * (parameters.A_DC_max/100.) / parameters.FOV; (*0.1mT/mm)
  parameters.gradient_field_yz = parameters.gradient_field_x * 10.;2.5;20. ; (*0.1mT/mm)
  parameters.recalc_flag = 1.
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_f_sin_event, ev
  widget_control, ev.id, get_value = f
  widget_control, ev.top, get_uvalue = parameters
  parameters.f = f;
  parameters.recalc_flag = 1.
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_harmonics_n_event, ev
  widget_control, ev.id, get_value = harmonics_n
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.harmonics_n NE FLOOR(parameters.harmonics_n) THEN BEGIN
    parameters.harmonics_n = harmonics_n;
    widget_text_harmonics_n = widget_info(ev.top, find_by_uname='text_harmonics_n')
    widget_control, widget_text_harmonics_n, set_value = STRING(harmonics_n, format='(f3.1)')
    parameters.recalc_flag = 1.
  ENDIF ELSE BEGIN
    parameters.harmonics_n = harmonics_n;
    widget_text_harmonics_n = widget_info(ev.top, find_by_uname='text_harmonics_n')
    widget_control, widget_text_harmonics_n, set_value = STRING(harmonics_n, format='(f3.1)')
    parameters.recalc_flag = 0.
  ENDELSE
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_text_harmonics_n_event, ev
  widget_control, ev.id, get_value = harmonics_n
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.harmonics_n NE FLOAT(harmonics_n[0]) THEN BEGIN
    parameters.harmonics_n = FLOAT(harmonics_n[0])
    ;  print, 'harmonic number', FLOAT(harmonics_n[0]);
    parameters.recalc_flag = 1.
  ENDIF ELSE BEGIN
    parameters.recalc_flag = 0.
  ENDELSE
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

function FMMPI_droplist_particle_type_event, ev
  particle_type = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.particle_type NE particle_type THEN BEGIN
    parameters.particle_type = particle_type;
    parameters.recalc_flag = 1.
  ENDIF
  
  IF particle_type EQ 1 THEN BEGIN
    widget_coercivity = widget_info(ev.top, find_by_uname='coercivity')
    widget_control, widget_coercivity, SENSITIVE=1
  ENDIF ELSE BEGIN
    widget_coercivity = widget_info(ev.top, find_by_uname='coercivity')
    widget_control, widget_coercivity, SENSITIVE=0
  ENDELSE
  IF particle_type EQ 4 THEN BEGIN
    widget_viscosity = widget_info(ev.top, find_by_uname='viscosity')
    widget_control, widget_viscosity, SENSITIVE=1
    widget_temperature = widget_info(ev.top, find_by_uname='temperature')
    widget_control, widget_temperature, SENSITIVE=1
  ENDIF ELSE BEGIN
    widget_viscosity = widget_info(ev.top, find_by_uname='viscosity')
    widget_control, widget_viscosity, SENSITIVE=0
    widget_temperature = widget_info(ev.top, find_by_uname='temperature')
    widget_control, widget_temperature, SENSITIVE=0
  ENDELSE
  
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
  return, 0
end

function FMMPI_droplist_Radon_rec_event, ev
  radon_rec = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.radon_rec NE radon_rec THEN BEGIN
    parameters.radon_rec = radon_rec;
    parameters.recalc_flag = 0.
    sin_simulation_Chebyshev,ev.top, parameters
  ENDIF

  widget_control, ev.top, set_uvalue = parameters
  return, 0
end

function FMMPI_droplist_phase_scatter_plot_event, ev
  phase_scatter_plot = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.phase_scatter_plot NE phase_scatter_plot THEN BEGIN
    parameters.phase_scatter_plot = phase_scatter_plot;
    parameters.recalc_flag = 0.
    IF phase_scatter_plot EQ 1 THEN sin_simulation_Chebyshev,ev.top, parameters
  ENDIF

  widget_control, ev.top, set_uvalue = parameters
  return, 0
end

function FMMPI_droplist_AC_DC_mode_event, ev
  AC_DC_mode = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.AC_DC_mode NE AC_DC_mode THEN BEGIN
    parameters.AC_DC_mode = AC_DC_mode;
    parameters.recalc_flag = 1.
    sin_simulation_Chebyshev,ev.top, parameters
  ENDIF

  widget_control, ev.top, set_uvalue = parameters
  return, 0
end

function FMMPI_droplist_shift_angle_type_event, ev
  shift_angle_type = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  IF parameters.shift_angle_type NE shift_angle_type THEN BEGIN
    parameters.shift_angle_type = shift_angle_type;
    parameters.recalc_flag = 0.
    parameters.shift_angle_renew = 0.
    sin_simulation_Chebyshev,ev.top, parameters
  ENDIF

  widget_control, ev.top, set_uvalue = parameters
  return, 0
end

function FMMPI_droplist_harmonics_comp_event, ev
  harmonics_comp = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  parameters.harmonics_comp = harmonics_comp;
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
  return, 0
end

function FMMPI_droplist_signal_type_event, ev
  signal_type = WIDGET_INFO(ev.id, /DROPLIST_SELECT)
  widget_control, ev.top, get_uvalue = parameters
  parameters.signal_type = signal_type;
  IF signal_type NE 2 THEN BEGIN
    widget_control, parameters.wHarmonicsCompList, SENSITIVE=0
    widget_control, parameters.slider_harmonics_n, SENSITIVE=0
  ENDIF ELSE BEGIN
    widget_control, parameters.wHarmonicsCompList, SENSITIVE=1
    widget_control, parameters.slider_harmonics_n, SENSITIVE=1
  ENDELSE
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
  return, 0
end

pro FMMPI_slider_particle_size_sin_event, ev
  widget_control, ev.id, get_value = particle_size
  widget_control, ev.top, get_uvalue = parameters
  parameters.particle_size = particle_size;
  parameters.recalc_flag = 1.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_shift_angle_event, ev
  widget_control, ev.id, get_value = shift_angle
  widget_control, ev.top, get_uvalue = parameters
  parameters.shift_angle = shift_angle;
  parameters.recalc_flag = 0.
  parameters.shift_angle_renew = 1.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_viscosity_sin_event, ev
  widget_control, ev.id, get_value = viscosity
  widget_control, ev.top, get_uvalue = parameters
  parameters.viscosity = viscosity;
  IF parameters.particle_type EQ 4 THEN parameters.recalc_flag = 1. ELSE parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_temperature_sin_event, ev
  widget_control, ev.id, get_value = T_p
  widget_control, ev.top, get_uvalue = parameters
  parameters.T_p = T_p;
  IF parameters.particle_type EQ 4 THEN parameters.recalc_flag = 1. ELSE parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_tau_sin_event, ev
  widget_control, ev.id, get_value = tau
  widget_control, ev.top, get_uvalue = parameters
  parameters.tau = tau;
  IF parameters.particle_type EQ 4 THEN parameters.recalc_flag = 1. ELSE parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_coercivity_sin_event, ev
  widget_control, ev.id, get_value = H_coercivity
  widget_control, ev.top, get_uvalue = parameters
  parameters.H_coercivity = H_coercivity;
  IF parameters.particle_type EQ 1 THEN parameters.recalc_flag = 1. ELSE parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

;pro FMMPI_slider_gradient_field_x_event, ev
;  widget_control, ev.id, get_value = gradient_field_x
;  widget_control, ev.top, get_uvalue = parameters
;  parameters.gradient_field_x = gradient_field_x;
;  parameters.recalc_flag = 1.
;  widget_control, ev.top, set_uvalue = parameters
;;  print, 'gradient field: ', parameters.gradient_field, '*0.1 mT/mm'
;  sin_simulation_Chebyshev,ev.top, parameters  
;end

pro FMMPI_slider_gradient_field_yz_event, ev
  widget_control, ev.id, get_value = gradient_field_yz_current
  widget_control, ev.top, get_uvalue = parameters
  parameters.gradient_field_yz_current = gradient_field_yz_current;
  parameters.recalc_flag = 0.
  widget_control, ev.top, set_uvalue = parameters
  ;  print, 'gradient field: ', parameters.gradient_field, '*0.1 mT/mm'
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_PhantomSelection_event, ev
  widget_control, ev.top, get_uvalue = parameters
  filename = dialog_pickfile(TITLE='Select Phantom Image File', FILTER=['*.jpg', '*.png', '*.tif'], /MUST_EXIST, PATH=parameters.directory)
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
  
  ;@GJ, 2024/4/21, reset the optimization factor
  parameters.phantom_rec_sub_threshold = 0.
  parameters.phantom_rec_sub_threshold_opt = 0.
  ;@GJ, 2024/4/18
  widget_threshold = widget_info(ev.top, find_by_uname='threshold')
  widget_control, widget_threshold, set_value = parameters.phantom_rec_sub_threshold_opt

  IF STRLEN(filename) GT 1 THEN BEGIN
    parameters.filename = filename
    parameters.directory = FILE_DIRNAME(filename)
    parameters.recalc_flag = 0.
    parameters.phase_scatter_plot=0
    sin_simulation_Chebyshev,ev.top, parameters
  ENDIF ELSE BEGIN
    parameters.filename = 'gj_barPattern.tif'
    parameters.recalc_flag = 0.
    parameters.phase_scatter_plot=0
    sin_simulation_Chebyshev,ev.top, parameters
  ENDELSE
  
  
  widget_control, ev.top, set_uvalue = parameters
end



pro FMMPI_OutputSelection_event, ev
  widget_control, ev.top, get_uvalue = parameters
  rec_directory = dialog_pickfile(TITLE='Select Output Directory', /DIRECTORY, /MUST_EXIST, PATH=parameters.rec_directory)
  IF STRLEN(rec_directory) GT 1 THEN BEGIN
    parameters.rec_directory = rec_directory
    parameters.recalc_flag = 0.
    parameters.phase_scatter_plot=0
    sin_simulation_Chebyshev,ev.top, parameters
  ENDIF
  widget_control, ev.top, set_uvalue = parameters
end

pro FMMPI_slider_fov_event, ev
  widget_control, ev.id, get_value = FOV_current
  widget_control, ev.top, get_uvalue = parameters
  parameters.FOV_current = FOV_current;
  parameters.recalc_flag = 0.
  parameters.phase_scatter_plot=0
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_noise_level_event, ev
  widget_control, ev.id, get_value = noise_level
  widget_control, ev.top, get_uvalue = parameters
  parameters.noise_level = noise_level;
  parameters.recalc_flag = 1.
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_slider_threshold_event, ev
  widget_control, ev.id, get_value = phantom_rec_sub_threshold
  widget_control, ev.top, get_uvalue = parameters
  parameters.phantom_rec_sub_threshold = phantom_rec_sub_threshold;
  parameters.recalc_flag = 0.
  parameters.phase_scatter_plot=0
  widget_control, ev.top, set_uvalue = parameters
  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_button_threshold_auto_event, ev
  widget_control, ev.top, get_uvalue = parameters

  Rec_sub_threshold_optimization, ev.top, parameters
  parameters.recalc_flag = 0.
  parameters.phase_scatter_plot=0
  widget_control, ev.top, set_uvalue = parameters

  sin_simulation_Chebyshev,ev.top, parameters
end

pro FMMPI_button_A_batch_event, ev
  widget_control, ev.top, get_uvalue = parameters
  
  FOR index=1, 50 DO BEGIN
    parameters.A = index
    widget_A_sin = widget_info(ev.top, find_by_uname='AC_field')
    widget_control, widget_A_sin, set_value = parameters.A
    parameters.recalc_flag = 1.
    widget_control, ev.top, set_uvalue = parameters
    sin_simulation_Chebyshev,ev.top, parameters
  ENDFOR
end

pro FMMPI_button_particle_size_batch_event, ev
  widget_control, ev.top, get_uvalue = parameters

  FOR index=10, 100 DO BEGIN
    parameters.particle_size = index
    widget_particle_size = widget_info(ev.top, find_by_uname='particle_size')
    widget_control, widget_particle_size, set_value = parameters.particle_size
    parameters.recalc_flag = 1.
    widget_control, ev.top, set_uvalue = parameters
    sin_simulation_Chebyshev,ev.top, parameters
  ENDFOR
end

;@GJ, 2024/2/12, selecting the percentage by mouse click
pro FMMPI_draw1_event, ev
;  widget_control, ev.id, get_value = noise_level
  widget_control, ev.top, get_uvalue = parameters
  ;press the button
  IF ev.Type EQ 0 THEN BEGIN
    ;left mouse button
    IF ev.Press EQ 1 THEN BEGIN
      IF parameters.AC_DC_mode EQ 0. THEN BEGIN
        ;AC excitation mode
        parameters.AC_ex = (ev.x - (parameters.DRAWXSIZE/2.))/(parameters.DRAWXSIZE/2.) * (parameters.AC_ex_max-parameters.AC_ex_min) + parameters.AC_ex_min
        print, ev.x, ev.y, 'left', ', parameters.AC_ex ', parameters.AC_ex
      ENDIF ELSE BEGIN
        ;DC offset mode
        parameters.A_DC = (ev.x - (parameters.DRAWXSIZE/2.))/(parameters.DRAWXSIZE/2.) * parameters.A_DC_max
        print, ev.x, ev.y, 'left', ', parameters.A_DC ', parameters.A_DC
      ENDELSE
    ENDIF
    
    ;right mouse button
    IF ev.Press EQ 4 THEN BEGIN
      IF parameters.AC_DC_mode EQ 0. THEN BEGIN
        ;AC excitation mode
        parameters.AC_ex_1 = (ev.x - (parameters.DRAWXSIZE/2.))/(parameters.DRAWXSIZE/2.) * (parameters.AC_ex_max-parameters.AC_ex_min) + parameters.AC_ex_min
        print, ev.x, ev.y, 'left', ', parameters.AC_ex_1 ', parameters.AC_ex_1
      ENDIF ELSE BEGIN
        ;DC offset mode
        parameters.A_DC_1 = (ev.x - (parameters.DRAWXSIZE/2.))/(parameters.DRAWXSIZE/2.) * parameters.A_DC_max
        print, ev.x, ev.y, 'right', ', parameters.A_DC ', parameters.A_DC_1
      ENDELSE
    ENDIF
    
    parameters.recalc_flag = 0.
    sin_simulation_Chebyshev,ev.top, parameters
  ENDIF
 
  widget_control, ev.top, set_uvalue = parameters

end


pro FMMPI_quit_event, ev
  widget_control, ev.top, set_uvalue = parameters
  ;  Destroy widget hierarchy.
  ;
  WIDGET_CONTROL, ev.top, /DESTROY

  RETURN
end
;@GJ, 2023/8/13, Zhang Yifei helped developing this part
;@GJ, 2023/12/10, Fully modified to include the offset DC and Gradient fields for focal modulation MPI
;@GJ, 2023/12/26, Included Chebyshev and simulated MNPs
;@GJ, 2024/2/13, Mouse button click for selecting A_DC and A_DC_1
;@GJ, 2024/9/12, Fix the bug of mouse selection of A_DC and A_DC_1; Add a text widget to select any harmonics, e.g. 3.1 or 2.7
pro Sin_Excitation_Chebyshev_MPI

  ;@GJ, 2024/6/6, particles
  particle_type = 4.;5.;4.;4.;type
  particle_size = 30.;37.;23.;1.;30.;nm
  tau = 20  ;(us)
  viscosity = 1.; mPa.s
  T_p = 25.; degree
  ;  tau = 20  ;(us)
  H_coercivity = 13.; coercivity in mT

  ;@GJ, 2024/6/6, excitation oddֵ
  A = 8.;25.;7.;25.;0.;5.;20.;mT
  f = 2.;10;kHz
  A_DC = 50.;60.;59.;0.;42.;95.;%,@GJ, 2023/12/23, changing to percentage20.;mT
  A_DC_1 = 120.;105.;;105.;50.;;%additional DC field
  Dog_factor = 6.;9.;37.;23.;19.; /10.Donut * Dog factor 
  A_DC_max = 201.;%@GJ, 2023/12/30,
  harmonics_n = 2.;7.;3.
  harmonics_comp = 0.; Amplitude
  noise_level = 40.;60.;@GJ, 2023/12/22, in dB5.; noise level
  FOV_current = 140.;50.;15.;100.; mm
  
  ;GJ, 2024/6/6, keyhole
  A = 8.;25.;7.;25.;0.;5.;20.;mT
  f = 2.;10;kHz
  A_DC = 0.;60.;59.;0.;42.;95.;%,@GJ, 2023/12/23, changing to percentage20.;mT
  A_DC_1 = 53.;105.;;105.;50.;;%additional DC field
  Dog_factor = 16.;9.;37.;23.;19.; /10.Donut * Dog factor
  A_DC_max = 201.;%@GJ, 2023/12/30,
  harmonics_n = 3.;2.;7.;3.
  harmonics_comp = 0.; Amplitude
  noise_level = 40.;60.;@GJ, 2023/12/22, in dB5.; noise level
  FOV_current = 140.;50.;15.;100.; mm
  
  ;GJ, 2024/8/16, phase imaging
  A = 10.;6.;10.;25.;7.;25.;0.;5.;20.;mT
  f = 8.;5.;2.;10;kHz
  A_DC = 0.;60.;59.;0.;42.;95.;%,@GJ, 2023/12/23, changing to percentage20.;mT
  A_DC_1 = 75.;20.;60.;105.;;105.;50.;;%additional DC field
  Dog_factor = 10.;16.;9.;37.;23.;19.; /10.Donut * Dog factor
  A_DC_max = 201.;%@GJ, 2023/12/30,
  harmonics_n = 3.;2.;7.;3.
  harmonics_comp = 0.;3.;1.; Phase
  noise_level = 40.;60.;@GJ, 2023/12/22, in dB5.; noise level
    
  ;@GJ, 2024/6/6, FOV setting
  ;@GJ, 2024/10/8, large bore MPI
  FOV = 200.;150.;400.; mm
  FOV_current = 20.;40.;60.;140.;50.;15.;100.; mm
  n_x = 256.;512.;
  n_z = 256.;512.;
  ;@GJ, 2023/12/30, adjusting gradient field x
  gradient_field_x = A * 2. * (A_DC_max/100.) / FOV; (*0.1mT/mm)
  gradient_field_yz = gradient_field_x * 10.;2.5;20. ; (*0.1mT/mm)
  gradient_field_yz_current = 20.;3.;1.;20.; * 0.1 mT
  recalc_flag = 1.;@GJ, 2023/12/26, whether to recalculate
  harmonics_comp_N = 4.
  harmonics_n_max = 13.;31.
  signal_type = 2.;1.;2.
  Sx_array = DBLARR(harmonics_n_max+2, harmonics_comp_N, n_x, n_z) * 0.
  Sx_adiabatic_array = DBLARR(harmonics_n_max+2, harmonics_comp_N, n_x, n_z) * 0.
  Sx_array_kernel = DBLARR(harmonics_n_max+2, harmonics_comp_N, n_x, n_z) * 0.
  Sx_array_xzNonConv = DBLARR(harmonics_n_max+2, harmonics_comp_N, n_x, n_z) * 0.
  focal_line = DBLARR(n_z)
  focal_line_1 = DBLARR(n_z)
  focal_STED = DBLARR(n_z)
  A_G_Array = DBLARR(n_z)
  optimal_dog_factor = Dog_factor
  phantom_rec_sub_min=0.
  phantom_rec_sub_max=255.
  phantom_rec_sub_threshold=0.
  phantom_rec_sub_threshold_opt=0.
  radon_rec = 0.
  phase_scatter_plot = 0.
  AC_DC_mode = 0.;1.;
  AC_ex_min = .1;
  AC_ex_max = 5.;
  AC_ex = 1.; excitation AC amplitude
  AC_ex_1 = 1.1; excitation AC amplitude
  shift_angle = 280. ;70.
  shift_angle_type = 1. ; best resolution, 0 best SNR
  shift_angle_renew = 0. ;calculate original data
  
  device, decomposed=1
  device, get_screen_size = screenSize
  baseXSize = screenSize[0] / 4.5 * 4.3
  baseYSize = screenSize[1] * 0.85

  base = widget_base(title='STED-MPI Simulation', xsize = baseXSize, ysize = baseYSize, /column, MBAR=barBase)
  
  ;  Create the quit button
  ;
  wFileButton = WIDGET_BUTTON(barBase, VALUE= 'Phantom', /MENU)

  wFileSelectButton = WIDGET_BUTTON(wFileButton, $
    VALUE='Select a phantom', event_pro = 'FMMPI_PhantomSelection_event')
  wQuitButton = WIDGET_BUTTON(wFileButton, $
    VALUE='Quit', event_pro = 'FMMPI_quit_event')

  menubase = widget_base(base, xsize = baseXSize/5, /row)
  
  sliderbase0 = widget_base(menubase, /column,/frame)
  ;@GJ, 2024/8/25, base_particle_type
  base_particle_type = widget_base(sliderbase0, /column)
  wParticleLabel = WIDGET_LABEL(base_particle_type, $
    value='Particle Type')
  ParticleTypeList = ['Ideal (Chebyshev)', 'SFMIOs (Chebyshev)', 'Langevin (Chebyshev)', 'Langevin MNPs (Simulation)', 'Langevin+Debye (Simulation)', 'SFMIOs (Simulation)']
  wParticleDroplist = WIDGET_DROPLIST(base_particle_type, $
    value=ParticleTypeList, $
    event_func = 'FMMPI_droplist_particle_type_event')
  widget_control, wParticleDroplist, SET_DROPLIST_SELECT=particle_type
  base_Radon_rec = widget_base(sliderbase0, /column)
  wRadonRecLabel = WIDGET_LABEL(base_Radon_rec, $
    value='Radon Recon')
  wRadonRecDroplist = WIDGET_DROPLIST(base_Radon_rec, $
    value=['No', 'Yes'], $
    event_func = 'FMMPI_droplist_Radon_rec_event')
  widget_control, wRadonRecDroplist, SET_DROPLIST_SELECT=radon_rec
  base_phase_scatter_plot = widget_base(sliderbase0, /column)
  wPhaseScatterPlotLabel = WIDGET_LABEL(base_phase_scatter_plot, $
    value='Phase Scatter Plot')
  wPhaseScatterPlotlist = WIDGET_DROPLIST(base_phase_scatter_plot, $
    value=['No', 'Yes'], $
    event_func = 'FMMPI_droplist_Phase_scatter_plot_event')
  widget_control, wPhaseScatterPlotlist, SET_DROPLIST_SELECT=phase_scatter_plot
  
  base_AC_DC_mode = widget_base(sliderbase0, /column)
  wACDEModeLabel = WIDGET_LABEL(base_AC_DC_mode, $
    value='Excitation Mode')
  wPhaseScatterPlotlist = WIDGET_DROPLIST(base_AC_DC_mode, $
    value=['AC Amp', 'DC Offset'], $
    event_func = 'FMMPI_droplist_AC_DC_mode_event')
  widget_control, wPhaseScatterPlotlist, SET_DROPLIST_SELECT=AC_DC_mode
    
  sliderbase0P = widget_base(menubase, /column,/frame)
  ;base_particle_size
  base_particle_size = widget_base(sliderbase0P, /row)
  label_particle_size = widget_label(base_particle_size, value = "Particle size (nm)")
  slider_particle_size_sin = widget_slider(base_particle_size, xsize=100, MINIMUM=10, MAXIMUM=100, uname='particle_size', event_pro = 'FMMPI_slider_particle_size_sin_event')
  button_particle_size_batch = widget_button(base_particle_size, value='Batch', event_pro='FMMPI_button_particle_size_batch_event')
  ;base_viscosity
  base_viscosity = widget_base(sliderbase0P, /row)
  label_viscosity = widget_label(base_viscosity, value = "Viscosity (mPa.s)")
  slider_viscosity_sin = widget_slider(base_viscosity, xsize=150, MINIMUM=1, MAXIMUM=50, uname='viscosity', event_pro = 'FMMPI_slider_viscosity_sin_event')
  ;base temperature
  base_temperature = widget_base(sliderbase0P, /row)
  label_temperature = widget_label(base_temperature, value = "Temperature (degree)")
  slider_temperature_sin = widget_slider(base_temperature, xsize=150, MINIMUM=0, MAXIMUM=100, uname='temperature', event_pro = 'FMMPI_slider_temperature_sin_event')
  ;base tau
  base_tau = widget_base(sliderbase0P, /row)
  label_tau = widget_label(base_tau, value = "Relaxation time (us)")
  slider_tau_sin = widget_slider(base_tau, xsize=150, MINIMUM=0., MAXIMUM=100, uname='tau', event_pro = 'FMMPI_slider_tau_sin_event')
  ;base coercivity
  base_coercivity = widget_base(sliderbase0P, /row)
  label_coercivity = widget_label(base_coercivity, value = "SFMIO Coercivity (mT)")
  slider_coercivity_sin = widget_slider(base_coercivity, xsize=150, MINIMUM=0., MAXIMUM=100, uname='coercivity', event_pro = 'FMMPI_slider_coercivity_sin_event')
  
  sliderbase = widget_base(menubase, /column, /BASE_ALIGN_RIGHT, /frame)
  ;base_A
  base_A = widget_base(sliderbase, /row)
  label_A = widget_label(base_A, value = "AC field (mT)")
  slider_A_sin = widget_slider(base_A, xsize=100, MINIMUM=1, MAXIMUM=50, uname='AC_field', event_pro = 'FMMPI_slider_A_sin_event')
  button_A_batch = widget_button(base_A, value='Batch', event_pro='FMMPI_button_A_batch_event')
  ;base_f
  base_f = widget_base(sliderbase, /row)
  label_f = widget_label(base_f, value =  "AC frequency (kHz)" )
  slider_f_sin = widget_slider(base_f, xsize=150, MINIMUM=1, MAXIMUM=100, event_pro = 'FMMPI_slider_f_sin_event')
  ;base A_DC
  base_A_DC = widget_base(sliderbase, /row, uname='base_A_DC')
  label_A_DC = widget_label(base_A_DC, value = "DC offset (Red, %AC)", uname='label_A_DC')
  slider_A_DC_sin = widget_slider(base_A_DC, xsize=150, MINIMUM=-(A_DC_max-1), MAXIMUM=A_DC_max-1, uname='A_DC', event_pro = 'FMMPI_slider_A_DC_event')
  ;base A_DC_1
  base_A_DC_1 = widget_base(sliderbase, /row, uname='base_A_DC_1')
  label_A_DC_1 = widget_label(base_A_DC_1, value = "DC offset (Green, %AC)", uname='label_A_DC_1')
  slider_A_DC_1_sin = widget_slider(base_A_DC_1, xsize=150, MINIMUM=-(A_DC_max-1), MAXIMUM=A_DC_max-1, uname='A_DC_1', event_pro = 'FMMPI_slider_A_DC_1_event')
;  ;base AC_ex
;  base_AC_ex = widget_base(sliderbase, /row, uname='base_AC_ex')
;  label_AC_ex = widget_label(base_AC_ex, value = "AC excitation (Red, %AC)")
;  slider_AC_ex_sin = widget_slider(base_AC_ex, xsize=150, MINIMUM=AC_ex_min*100., MAXIMUM=AC_ex_max*100., uname='AC_ex', event_pro = 'FMMPI_slider_AC_ex_event')
;  widget_control, base_AC_ex, show=0
;  ;base AC_ex_1
;  base_AC_ex_1 = widget_base(sliderbase, /row, uname='base_AC_ex_1')
;  label_AC_ex_1 = widget_label(base_AC_ex_1, value = "AC excitation (Green, %AC)")
;  slider_AC_ex_1_sin = widget_slider(base_AC_ex_1, xsize=150, MINIMUM=AC_ex_min*100., MAXIMUM=AC_ex_max*100., uname='AC_ex_1', event_pro = 'FMMPI_slider_AC_ex_1_event')
;  widget_control, base_AC_ex_1, show=0
  
  ;base Dog factor
  base_Dog_factor = widget_base(sliderbase, /row)
  label_Dog_factor = widget_label(base_Dog_factor, value = "Dog factor (/10)")
  slider_Dog_factor = widget_slider(base_Dog_factor, xsize=120, MINIMUM=1, MAXIMUM=100, uname='dog_factor', event_pro = 'FMMPI_slider_Dog_factor_event')
  button_Dog_opt = widget_button(base_Dog_factor, value='Auto', event_pro='FMMPI_button_Dog_factor_auto_event')
;  ;base_gradient field x
;  base_gradient_field_x = widget_base(sliderbase, /row)
;  label_gradient_field_x = widget_label(base_gradient_field_x, value = "gradient excitation direc (*0.1mT/mm):")
;  slider_gradient_field_x = widget_slider(base_gradient_field_x, xsize=150, MINIMUM=1, MAXIMUM=100, event_pro = 'FMMPI_slider_gradient_field_x_event')

  sliderbase2 = widget_base(menubase, /column, /frame)
  ;base_signal_type
  base_comp_signal = widget_base(sliderbase2, /col)
  base_signal_type_comp = widget_base(base_comp_signal, /row)
  wSignalTypeLabel = WIDGET_LABEL(base_signal_type_comp, $
    value='Signal Type')
  SignalTypeList = ['Filtered', 'Unfiltered', 'Harmonic']
  wSignalTypelist = WIDGET_DROPLIST(base_signal_type_comp, $
    value=SignalTypeList, $
    event_func = 'FMMPI_droplist_signal_type_event')
  widget_control, wSignalTypelist, SET_DROPLIST_SELECT=signal_type
   ;base_harmonics_comp
   base_harmonics_comp = widget_base(base_comp_signal, /row)
   wHarmonicsCompLabel = WIDGET_LABEL(base_harmonics_comp, $
     value='Harmonic Type')
   HarmonicsCompList = ['Amplitude', 'Phase', 'Real', 'Imaginary']
   wHarmonicsComplist = WIDGET_DROPLIST(base_harmonics_comp, $
     value=HarmonicsCompList, $
     event_func = 'FMMPI_droplist_harmonics_comp_event')
   widget_control, wHarmonicsComplist, SET_DROPLIST_SELECT=harmonics_comp
   ;shift_angle_type
   base_shift_angle_type = widget_base(base_comp_signal, /row)
   wShiftAngleTypeLabel = WIDGET_LABEL(base_shift_angle_type, $
     value='Shift Angle Type')
   wShiftAngleTypelist = WIDGET_DROPLIST(base_shift_angle_type, $
     value=['0 SNR', '1 Resolution'], $
     event_func = 'FMMPI_droplist_shift_angle_type_event')
   widget_control, wShiftAngleTypelist, SET_DROPLIST_SELECT=shift_angle_type
   
   ;Shift angle
   base_shift_angle = widget_base(base_comp_signal, /row)
   label_shift_angle = widget_label(base_shift_angle, value = "Shift Angle (degree)")
   slider_shift_angle = widget_slider(base_shift_angle, xsize=150, MINIMUM=0, MAXIMUM=360, uname='shift_angle', event_pro = 'FMMPI_slider_shift_angle_event')

  ;base_harmonics_n
  base_harmonics_n = widget_base(sliderbase2, /row)
  label_harmonics_n = widget_label(base_harmonics_n, value = "Harmonics No")
  slider_harmonics_n = widget_slider(base_harmonics_n, xsize=120, MINIMUM=1, MAXIMUM=harmonics_n_max, event_pro = 'FMMPI_slider_harmonics_n_event')
  text_harmonics_n = widget_text(base_harmonics_n, /EDITABLE, xsize=5, uname='text_harmonics_n', event_pro = 'FMMPI_text_harmonics_n_event')
    
  phantombase = widget_base(menubase, /column, /BASE_ALIGN_RIGHT, /frame)
  base_selection = widget_base(phantombase, /row)
  wFileSelectButton = WIDGET_BUTTON(base_selection, value='Select phantom image', event_pro='FMMPI_PhantomSelection_event')
  wFileSelectButton_output = WIDGET_BUTTON(base_selection, value='Select output path', event_pro='FMMPI_OutputSelection_event')
  ;base_noise_level
  base_noise_level = widget_base(phantombase, /row)
  label_noise_level = widget_label(base_noise_level, value = "Noise (dB)")
  slider_noise_level = widget_slider(base_noise_level, xsize=150, MINIMUM=1, MAXIMUM=100, event_pro = 'FMMPI_slider_noise_level_event')
  ;base_FOV
  base_FOV = widget_base(phantombase, /row)
  label_FOV = widget_label(base_FOV, value = "FOV (mm)")
  slider_FOV = widget_slider(base_FOV, xsize=150, MINIMUM=1, MAXIMUM=FOV, event_pro = 'FMMPI_slider_FOV_event')
  ;base_gradient field yz
  base_gradient_field_yz = widget_base(phantombase, /row)
  label_gradient_field_yz = widget_label(base_gradient_field_yz, value = "Gradient (*0.1T/m):")
  slider_gradient_field_yz = widget_slider(base_gradient_field_yz, xsize=150, MINIMUM=1, MAXIMUM=100, event_pro = 'FMMPI_slider_gradient_field_yz_event')
  base_threshold = widget_base(phantombase, /row)
  label_threshold = widget_label(base_threshold, value = "Threshold")
  slider_threshold = widget_slider(base_threshold, xsize=120, MINIMUM=phantom_rec_sub_min, MAXIMUM=phantom_rec_sub_max, uname='threshold', event_pro = 'FMMPI_slider_threshold_event')
  button_threshold_opt = widget_button(base_threshold, value='Auto', event_pro='FMMPI_button_threshold_auto_event')


  drawXSize = baseXSize / 7.5
  drawYSize = drawXSize;baseYSize / 3.5
  drawbase = widget_base(base, /column)
  drawbase1 = widget_base(drawbase, /row)
  draw1 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw1', /BUTTON_EVENTS, event_pro = 'FMMPI_draw1_event')
  draw2 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw2')
  draw3 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw3')
  draw4 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw4')
  draw5 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw5')
  draw6 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw6')
  draw7 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw7')
  drawbase2 = widget_base(drawbase, /row)
  draw11 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw11')
  draw12 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw12')
  draw13 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw13')
  draw14 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw14')
  draw15 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw15')
  draw16 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw16')
  draw17 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw17')
  drawbase3 = widget_base(drawbase, /row)
  draw21 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw21')
  draw22 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw22')
  draw23 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw23')
  draw24 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw24')
  draw25 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw25')
  draw26 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw26')
  draw27 = widget_draw(drawbase3, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw27')
  
  widget_control, base, /realize
  
  parameters = { $
    A:A, $
    f:f, $
    A_DC:A_DC, $
    A_DC_1:A_DC_1, $
    A_DC_max:A_DC_max, $
    A_DC_threshold_donut:A_DC_1, $ ;@GJ, 2024/5/1, threshold for donut appearance
    donut_radius0:0., $ ;@GJ, 2024/5/4, the initial radius;
    donut_slope:0., $ ;@GJ, 2024/5/4, the slope
    Dog_factor:Dog_factor, $
    focal_line:PTR_NEW(focal_line), $
    focal_line_1:PTR_NEW(focal_line_1), $
    focal_line_real:PTR_NEW(focal_line), $
    focal_line_1_real:PTR_NEW(focal_line_1), $
    focal_line_imag:PTR_NEW(focal_line), $
    focal_line_1_imag:PTR_NEW(focal_line_1), $
    focal_STED:PTR_NEW(focal_STED), $
    A_G_Array:PTR_NEW(A_G_Array), $
    optimal_dog_factor:optimal_dog_factor, $
    phantom_rec_sub_min:phantom_rec_sub_min, $
    phantom_rec_sub_max:phantom_rec_sub_max, $
    phantom_rec_sub_threshold:phantom_rec_sub_threshold, $
    phantom_rec_sub_threshold_opt:phantom_rec_sub_threshold_opt, $
    rec_noisy_image_sub:PTR_NEW(rec_noisy_image_sub), $
    radon_rec:radon_rec, $
    phase_scatter_plot:phase_scatter_plot, $
    AC_DC_mode:AC_DC_mode, $; for different excitation mode
    AC_ex_min:AC_ex_min, $; @GJ, 2024/10/7, 
    AC_ex_max:AC_ex_max, $; @GJ, 2024/10/7, 
    AC_ex:AC_ex, $; @GJ, 2024/10/7, 
    AC_ex_1:AC_ex_1, $; @GJ, 2024/10/7, 
    shift_angle: shift_angle, $
    shift_angle_type: shift_angle_type, $
    shift_angle_renew:shift_angle_renew, $
    image_sim_ori:PTR_NEW(image_sim_ori), $
    particle_type:particle_type, $
    particle_size:particle_size, $
    viscosity:viscosity, $
    T_p:T_p, $
    tau:tau, $
    H_coercivity:H_coercivity, $
    noise_level:noise_level, $
    harmonics_n:harmonics_n, $
    harmonics_comp_N: harmonics_comp_N, $
    harmonics_n_max:harmonics_n_max, $
    harmonics_comp:harmonics_comp, $
    signal_type: signal_type, $
    gradient_field_x:gradient_field_x, $
    gradient_field_yz:gradient_field_yz, $
    gradient_field_yz_current:gradient_field_yz_current, $
    drawXSize:drawXSize, $
    drawYSize:drawYSize, $
    directory:'.\lib\icon_pics\phantom\', $
    rec_directory:'.\lib\icon_pics\phantom\rec'+STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(i0)')+'\', $
    filename:'gj_barPattern.tif', $
    FOV_current:FOV_current, $
    FOV:FOV, $
    n_x:n_x, $ = 512
    n_z:n_z, $ = 512
    Sx_array:PTR_NEW(Sx_array), $
    Sx_adiabatic_array:PTR_NEW(Sx_adiabatic_array), $ ; @GJ, 2024/9/9, storage of adiabatic signal
    Sx_array_kernel:PTR_NEW(Sx_array_kernel), $ = DBLARR(harmonics_n_max, harmonics_comp_N, n_x, n_z) * 0.
    Sx_array_xzNonCOnv: PTR_NEW(Sx_array_xzNonCOnv), $= DBLARR(harmonics_n_max, harmonics_comp_N, n_x, n_z) * 0.
    tau_ms_array:PTR_NEW(DBLARR(n_x, n_z)), $ ; the array of relaxation time
    recalc_flag:recalc_flag, $
    draw1:draw1, $
    wSignalTypelist:wSignalTypelist, $
    wHarmonicsCompList:wHarmonicsCompList, $
    slider_shift_angle:slider_shift_angle, $
    slider_harmonics_n:slider_harmonics_n, $
    slider_coercivity_sin:slider_coercivity_sin $
    }
    
  widget_control, base, set_uvalue = parameters

  ;����slider�����value��ʼֵ
  widget_control, slider_A_sin, set_value = A
  widget_control, slider_f_sin, set_value = f
  widget_control, slider_A_DC_sin, set_value = A_DC
  widget_control, slider_A_DC_1_sin, set_value = A_DC_1
  widget_control, slider_Dog_factor, set_value = Dog_factor
  widget_control, slider_particle_size_sin, set_value = particle_size
  widget_control, slider_viscosity_sin, set_value = viscosity
  widget_control, slider_temperature_sin, set_value = T_p
  ;@GJ, 2024/1/1, calculate the tau
  tau = brownian_time_calc(parameters)
  parameters.tau = tau
  widget_control, slider_tau_sin, set_value = tau, SENSITIVE=0
  widget_control, slider_coercivity_sin, set_value = H_coercivity
  IF particle_type NE 1 THEN widget_control, slider_coercivity_sin, SENSITIVE=0
  widget_control, slider_noise_level, set_value = noise_level
  widget_control, slider_harmonics_n, set_value = harmonics_n
  widget_control, slider_shift_angle, set_value = shift_angle
  widget_control, text_harmonics_n, set_value = STRING(harmonics_n, format='(f3.1)')
  widget_control, slider_FOV, set_value = FOV_current
  ;@GJ, 2023/12/30, removing the gradient field x for automatic calculation
;  widget_control, slider_gradient_field_x, set_value = gradient_field_x
  widget_control, slider_gradient_field_yz, set_value = gradient_field_yz_current
  widget_control, slider_threshold, set_value = phantom_rec_sub_threshold_opt
  sin_simulation_Chebyshev, base, parameters
  xmanager, 'Sin_Excitation Cheybshev', base
end