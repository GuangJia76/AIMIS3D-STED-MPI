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


;�������,�������¹�ʽ��
;ev Ϊ�¼�ID
;������       H = -Acos(2*pi*f*t) 
;��֮�򷽳�   M = 1/tanh(beta*H) - 1/(beta*H)
;�źų�ԥ       R = (1/tau)*exp(-T/tau)
pro simulation,wid,A,particle_size,f,tau
  Msat_T = 0.551 ; T/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  T = findgen(200)/40     ;(1  ms)
  H = -A*cos(2*!PI*f*T) ;(mT)
  
  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT
  LangzwX = (FINDGEN(401)-200.)/10.
  Langzw = Msat_kAm*(1./TANH(beta*LangzwX) - 1./(beta*LangzwX))
  Langzw[200] = 0
  
  M = 1/tanh(beta*H) - 1/(beta*H) ;(kA/m)
  R = (1/tau)*exp(-T/tau)
  n = N_elements(M)
  signal = M[1:*] - M[0:n-1]    ; 
  r = (1/tau)*exp(-T/tau)
  kernel = dblarr(n-1)*0
  kernel[0:*] = r[0:n-2]
  signal_new = convol(signal, kernel, /EDGE_ZERO)
  
  !P.FONT = 0
  DEVICE,SET_FONT = "����*16"
  cgLoadCT, 34
  ;ͨ��widget_info�����������unameΪ��ʶ���в��ң��������id
  draw1 = widget_info(wid, find_by_uname='draw1')
  widget_control, draw1, get_value = drawWindow1
  ;wset����draw1���Ϊ�滭����
  wset, drawWindow1
  ;cgScaleVector���ڽ������������Ԫ������Ϊ����ݷ�Χ��
  colors = cgScaleVector(Findgen(N_Elements(H)), Min(H), Max(H))
  ;Value_Locate���ڽ�  data ��  colors(��������)��ʾ�� ��ͬ��dataֵ�в�ͬ��colors
  elevColors = Value_Locate(colors, H)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  cgplot, T, H, /NoData, YRANGE=[-100, 100], background = 'ffffff'x, color = '000000'x, xtitle='time(ms)'
  xyouts,5,360,' magnetic!c field!cintensity!c  (mT)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(H)-2 $
    DO cgPlotS, [T[j], T[j+1]], [H[j], H[j+1]], Color=elevColors[j]
  
  
  draw2 = widget_info(wid, find_by_uname='draw2')
  widget_control, draw2, get_value = drawWindow2
  wset, drawWindow2
  colors = cgScaleVector(Findgen(N_Elements(Langzw)), Min(Langzw), Max(Langzw))
  elevColors = Value_Locate(colors, Langzw)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  cgplot, LangzwX, Langzw, /NoData, background = 'ffffff'x, color = '000000'x, xtitle='magnetization(kA/m)'
  xyouts,5,360,'magnetic!c field!cintensity!c (mT)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(Langzw)-2 $
    DO cgPlotS, [LangzwX[j], LangzwX[j+1]], [Langzw[j], Langzw[j+1]], Color=(elevColors[j]+elevColors[j+1])/2
  
  
  draw3 = widget_info(wid, find_by_uname='draw3')
  widget_control, draw3, get_value = drawWindow3
  wset, drawWindow3
  colors = cgScaleVector(Findgen(N_Elements(M)), Min(M), Max(M))
  elevColors = Value_Locate(colors, M)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  cgplot, T, M, /NoData, background = 'ffffff'x, color = '000000'x, xtitle='time(ms)'
  xyouts,5,360,'magnetization!c   (kA/m)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(M)-2 $
    DO cgPlotS, [T[j], T[j+1]], [M[j], M[j+1]], Color=(elevColors[j]+elevColors[j+1])/2
  
  draw4 = widget_info(wid, find_by_uname='draw4')
  widget_control, draw4, get_value = drawWindow4
  wset, drawWindow4
  colors = cgScaleVector(Findgen(N_Elements(signal)), Min(signal), Max(signal))
  elevColors = Value_Locate(colors, signal)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  cgplot, T, signal, /NoData, background = 'ffffff'x, color = '000000'x, xtitle='time(ms)'
  xyouts,5,360,'magnetization!c   (kA/m)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(signal)-2 $
    DO cgPlotS, [T[j], T[j+1]], [signal[j], signal[j+1]], Color=(elevColors[j]+elevColors[j+1])/2
  
  draw5 = widget_info(wid, find_by_uname='draw5')
  widget_control, draw5, get_value = drawWindow5
  wset, drawWindow5
  u = ABS(FFT(signal))
  frequency = findgen(N_Elements(u))
  colors = cgScaleVector(Findgen(N_Elements(u)), Min(u), Max(u))
  elevColors = Value_Locate(colors, u)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  cgplot, frequency, u,  XRANGE=[0, 100],background = 'ffffff'x, color = '000000'x, xtitle='frequency'
  xyouts,20,350,'u', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(u)/2 $
    DO cgPlotS, [frequency[j], frequency[j+1]], [u[j], u[j+1]], Color=(elevColors[j]+elevColors[j+1])/2
  
  ; Clean up.
  !P.Multi = 0

end

pro slider_A_event, ev
  widget_control, ev.id, get_value = A
  widget_control, ev.top, get_uvalue = parameters
  parameters.A = A;
  widget_control, ev.top, set_uvalue = parameters

  simulation,ev.top, parameters.A, parameters.particle_size, parameters.f, parameters.tau
end

pro slider_particle_size_event, ev
  widget_control, ev.id, get_value = particle_size
  widget_control, ev.top, get_uvalue = parameters
  parameters.particle_size = particle_size;
  widget_control, ev.top, set_uvalue = parameters
  
  simulation,ev.top, parameters.A, parameters.particle_size, parameters.f, parameters.tau
end

pro slider_f_event, ev
  widget_control, ev.id, get_value = f
  widget_control, ev.top, get_uvalue = parameters
  parameters.f = f;
  widget_control, ev.top, set_uvalue = parameters

  simulation,ev.top, parameters.A, parameters.particle_size, parameters.f, parameters.tau
end

pro MPI_simulation
  ;��ʼֵ
  A = 20
  particle_size = 30
  f = 1
  tau = 2
  parameters = {A:A, particle_size:particle_size, f:f, tau:tau}
  

  
  base = widget_base(title='MPI', /column)
  
  menubase = widget_base(base, /row)
  
  sliderbase = widget_base(menubase, /column)
  ;base_A
  base_A = widget_base(sliderbase, /row)
  label_A = widget_label(base_A, value = "Excitation wave amplitude(mT) : ")
  slider_A = widget_slider(base_A, xsize=400, MINIMUM=1, event_pro = 'slider_A_event')
  ;base_particle_size
  base_particle_size = widget_base(sliderbase, /row)
  label_particle_size = widget_label(base_particle_size, value = "Particle size(nm)             :")
  slider_particle_size = widget_slider(base_particle_size, xsize=400, MINIMUM=10, MAXIMUM=100, event_pro = 'slider_particle_size_event')
  ;base_f
  base_f = widget_base(sliderbase, /row)
  label_f = widget_label(base_f, value =  "Excitation wave frequency(kHz):" )
  slider_f = widget_slider(base_f, xsize=400, MINIMUM=1, MAXIMUM=10, event_pro = 'slider_f_event')
  
  picbase = widget_base(menubase, /column)
  MPIjpg_path = ".\lib\icon_pics\MPI_sim.jpg"
  READ_JPEG, MPIjpg_path, MPIjpg, true=1
  picwindow = WIDGET_draw(picbase, uvalue='pic_win', xsize = (size(MPIjpg))[2], ysize = (size(MPIjpg))[3])

  ;��ͼ��
  drawXSize = 550
  drawYSize = 390
  drawbase = widget_base(base, /column)
  drawbase1 = widget_base(drawbase, /row)
  draw1 = widget_draw(drawbase1, xsize=drawXSize, ysize=drawYSize, retain=2, uname='draw1')
  draw2 = widget_draw(drawbase1, xsize=drawXSize, ysize=drawYSize, retain=2, uname='draw2')
  draw3 = widget_draw(drawbase1, xsize=drawXSize, ysize=drawYSize, retain=2, uname='draw3')
  drawbase2 = widget_base(drawbase, /row)
  draw4 = widget_draw(drawbase2, xsize=drawXSize, ysize=drawYSize, retain=2, uname='draw4')
  draw5 = widget_draw(drawbase2, xsize=drawXSize, ysize=drawYSize, retain=2, uname='draw5')
  

  widget_control, base, /realize
  
  widget_control, picwindow, get_value = drawPic
  wset, drawPic
  TVSCL, MPIjpg, true=1
  
  widget_control, base, set_uvalue = parameters
  
  ;����slider�����value��ʼֵ
  widget_control, slider_A, set_value = A
  widget_control, slider_particle_size, set_value = particle_size
  widget_control, slider_f, set_value = f
  simulation, base, parameters.A, parameters.particle_size, parameters.f, parameters.tau
  xmanager, 'MPI_simulation', base
end