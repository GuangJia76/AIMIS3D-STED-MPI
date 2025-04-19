function Debye_Kernel, tau, kernel_time
  N = N_ELEMENTS(kernel_time)
  kernel = DBLARR(N)*0.
  d_t = kernel_time[1]-kernel_time[0] ;us
  time = FINDGEN(N) * d_t  ; us
  r = (1./tau) * exp(-time/tau)
  kernel[0:N/2] = REVERSE(r[0:N/2])

  return, kernel
end

pro pulse_simulation,wid,A,flat_portion,particle_size,f,tau

  tau = tau/1000.
  flat_portion = flat_portion/100.
  Msat_T = 0.551 ; Time/u0
  Msat_kAm = 0.551/4./!PI*10000. ; kA/ m
  Time = findgen(500)/100.     ;(1  ms)
  
  num_T = 100/f
  num_T1 = num_T/2.*flat_portion
  num_T2 = num_T/2. - num_T1
  num_T3 = (num_T-num_T/2.)*flat_portion
  num_T4 = (num_T-num_T/2.) - num_T3
  H1 = -A * (dblarr(num_T1)+1.)
  H2 = findgen(num_T2) / (num_T2-1) * 2 * A - A
  H3 = A * (dblarr(num_T3)+1.)
  H4 = -(findgen(num_T4) / (num_T4-1) * 2 * A - A)
  tmp_H = [H1, H2, H3, H4]
  H = []  ;(mT)
  for i=1,5*f do $
    H = [H, tmp_H]
  num_H = n_elements(H)
;  repeat begin
;    print, n_elements(H)
;    H = [H, tmp_H]
;    num_H = n_elements(H)
;  endrep until (num_H GT 500)
;  H = H[0 : 499]

  
  beta = Msat_T*((particle_size*1.)^3) /24./1.380469/309.65 ; in 1/mT

  M = Msat_kAm*(1/tanh(beta*H) - 1/(beta*H)) ;(kA/m)
  for i=0,n_elements(M)-1 do $
    if FINITE(M[i], /NAN) then M[i]=0
    
  R = (1/tau)*exp(-Time/tau)
  n = N_elements(M)
  signal = M[1:*] - M[0:n-1]
  
  kernel = Debye_Kernel(tau, Time[0:100])
  
  
  M = convol([dblarr(n_elements(kernel)/2)*0, M, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
  M = M[n_elements(kernel)/2 : n_elements(M) - n_elements(kernel)/2 - 1]

  signal_new = convol([dblarr(n_elements(kernel)/2)*0, signal, dblarr(n_elements(kernel)/2)*0], kernel, /NORMALIZE)
  signal_new = signal_new[n_elements(kernel)/2 : n_elements(signal) + n_elements(kernel)/2]
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
  plot, Time, H, /NoData, YRANGE=[-100, 100], background = 'ffffff'x, color = '000000'x, xtitle='time(ms)'
  xyouts,5,drawSize.ysize*0.95,'magnetic field intensity (mT)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(H)-2 $
    DO cgPlotS, [Time[j], Time[j+1]], [H[j], H[j+1]], Color=elevColors[j+1]


  draw2 = widget_info(wid, find_by_uname='draw2')
  widget_control, draw2, get_value = drawWindow2
  wset, drawWindow2
  H_2 = H[N_Elements(H)/2 : N_Elements(H)-1]  ; ��ȡһ�� ��ֹ��һ�����ھ����ֶ�����
  M_2 = M[N_Elements(M)/2 : N_Elements(M)-1]
  colors = cgScaleVector(Findgen(N_Elements(M_2)), Min(M_2), Max(M_2))
  elevColors = Value_Locate(colors, M_2)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  plot, H_2, M_2, /NoData,  background = 'ffffff'x, color = '000000'x, xtitle='magnetic field intensity (mT)'
  xyouts,5,drawSize.ysize*0.95,'magnetization(kA/m)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(H_2)-2 $
    DO cgPlotS, [H_2[j], H_2[j+1]], [M_2[j], M_2[j+1]], Color=elevColors[j+1]

  draw3 = widget_info(wid, find_by_uname='draw3')
  widget_control, draw3, get_value = drawWindow3
  wset, drawWindow3
  colors = cgScaleVector(Findgen(N_Elements(M)), Min(M), Max(M))
  elevColors = Value_Locate(colors, M)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  plot, Time, M, /NoData, background = 'ffffff'x, color = '000000'x, xtitle='time(ms)'
  xyouts,5,drawSize.ysize*0.95,'magnetization (kA/m)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(M)-2 $
    DO cgPlotS, [Time[j], Time[j+1]], [M[j], M[j+1]], Color=elevColors[j+1]

  draw4 = widget_info(wid, find_by_uname='draw4')
  widget_control, draw4, get_value = drawWindow4
  wset, drawWindow4
  colors = cgScaleVector(Findgen(N_Elements(signal_new)), Min(signal_new), Max(signal_new))
  elevColors = Value_Locate(colors, signal_new)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  plot, Time, signal,/NoData, background = 'ffffff'x, color = '000000'x, xtitle='time(ms)'
  xyouts,5,drawSize.ysize*0.95,'signal (kA/m/ms)', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(signal_new)-2 $
    DO cgPlotS, [Time[j], Time[j+1]], [signal_new[j], signal_new[j+1]], Color=elevColors[j+1]

  draw5 = widget_info(wid, find_by_uname='draw5')
  widget_control, draw5, get_value = drawWindow5
  wset, drawWindow5
  u = ABS(FFT(signal_new))
  frequency = findgen(N_Elements(u))
  colors = cgScaleVector(Findgen(N_Elements(u)), Min(u), Max(u))
  elevColors = Value_Locate(colors, u)
  elevColors = Byte(Round(cgScaleVector(elevColors, 0, 255)))
  plot, frequency, u, /NoData, XRANGE=[0, N_Elements(u)/2], background = 'ffffff'x, color = '000000'x, xtitle='frequency'
  xyouts,20,drawSize.ysize*0.95,'u', color = '000000'x, /DEVICE
  FOR j=0,N_Elements(u)/2-1 do $
    cgPlotS, [frequency[j], frequency[j+1]], [u[j], u[j+1]], Color=elevColors[j+1]
  cgPlotS, [frequency[0], frequency[N_Elements(u)/2]], [0, 0], Color=black

  ; Clean up.
  !P.Multi = 0

end

pro slider_A_event, ev
  widget_control, ev.id, get_value = A
  widget_control, ev.top, get_uvalue = parameters
  parameters.A = A;
  widget_control, ev.top, set_uvalue = parameters

  pulse_simulation,ev.top, parameters.A,parameters.portion,parameters.particle_size, parameters.f, parameters.tau
end

pro slider_portion_event, ev
  widget_control, ev.id, get_value = portion
  widget_control, ev.top, get_uvalue = parameters
  parameters.portion = portion;
  widget_control, ev.top, set_uvalue = parameters

  pulse_simulation,ev.top, parameters.A,parameters.portion,parameters.particle_size, parameters.f, parameters.tau
end

pro slider_f_event, ev
  widget_control, ev.id, get_value = f
  widget_control, ev.top, get_uvalue = parameters
  parameters.f = f;
  widget_control, ev.top, set_uvalue = parameters

  pulse_simulation,ev.top, parameters.A,parameters.portion,parameters.particle_size, parameters.f, parameters.tau
end

pro slider_particle_size_event, ev
  widget_control, ev.id, get_value = particle_size
  widget_control, ev.top, get_uvalue = parameters
  parameters.particle_size = particle_size;
  widget_control, ev.top, set_uvalue = parameters

  pulse_simulation,ev.top, parameters.A,parameters.portion,parameters.particle_size, parameters.f, parameters.tau
end

pro slider_tao_event, ev
  widget_control, ev.id, get_value = tau
  widget_control, ev.top, get_uvalue = parameters
  parameters.tau = tau;
  widget_control, ev.top, set_uvalue = parameters

  pulse_simulation,ev.top, parameters.A,parameters.portion,parameters.particle_size, parameters.f, parameters.tau
end

pro pulse_Excitation

  ;��ʼֵ
  A = 20
  particle_size = 30
  f = 1
  tau = 20  ;(ms)
  flat_portion = 50
  parameters = {A:A, portion:flat_portion, particle_size:particle_size, f:f, tau:tau}
  device, decomposed=1
  device, get_screen_size = screenSize
  baseXSize = screenSize[0] * 0.8
  baseYSize = screenSize[1] * 0.8

  base = widget_base(title='Pluse Excitation', xsize = baseXSize, ysize = baseYSize, /column)

  menubase = widget_base(base, xsize = baseXSize/3, /row)

  sliderbase = widget_base(menubase, /column)
  ;base_A
  base_A = widget_base(sliderbase, /row)
  label_A = widget_label(base_A, value = "Excitation wave amplitude(mT) :")
  slider_A = widget_slider(base_A, xsize=400, MINIMUM=1, event_pro = 'slider_A_event')
  ;base_portion  flat_portion
  base_portion = widget_base(sliderbase, /row)
  label_portion = widget_label(base_portion, value = "flat_portion(%)               :")
  slider_portion = widget_slider(base_portion, xsize=400, MINIMUM=50, MAXIMUM=95, event_pro = 'slider_portion_event')
  ;base_f
  base_f = widget_base(sliderbase, /row)
  label_f = widget_label(base_f, value =  "Excitation wave frequency(kHz):" )
  slider_f = widget_slider(base_f, xsize=400, MINIMUM=1, MAXIMUM=5, event_pro = 'slider_f_event')
  ;base_particle_size
  base_particle_size = widget_base(sliderbase, /row)
  label_particle_size = widget_label(base_particle_size, value = "Particle size(nm)             :")
  slider_particle_size = widget_slider(base_particle_size, xsize=400, MINIMUM=10, MAXIMUM=100, event_pro = 'slider_particle_size_event')
  ;base_tao
  base_tao = widget_base(sliderbase, /row)
  label_tao = widget_label(base_tao, value = "relaxation time(us)           :")
  slider_tao = widget_slider(base_tao, xsize=400, MINIMUM=1, MAXIMUM=100, event_pro = 'slider_tao_event')


  picbase = widget_base(menubase, /column)
  MPIjpg_path = './lib/icon_pics/pulsedExcitation.jpg'
  READ_JPEG, MPIjpg_path, MPIjpg, true=1
  picwindow = WIDGET_draw(picbase, uvalue='pic_win', xsize = (size(MPIjpg))[2], ysize = (size(MPIjpg))[3])

  ;��ͼ��
  drawXSize = baseXSize / 3.1
  drawYSize = baseYSize / 3
  drawbase = widget_base(base, /column)
  drawbase1 = widget_base(drawbase, /row)
  draw1 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw1')
  draw2 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw2')
  draw3 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw3')
  drawbase2 = widget_base(drawbase, /row)
  draw4 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw4')
  draw5 = widget_draw(drawbase2, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw5')


  widget_control, base, /realize

  widget_control, picwindow, get_value = drawPic
  wset, drawPic
  TVSCL, MPIjpg, true=1

  widget_control, base, set_uvalue = parameters

  ;����slider�����value��ʼֵ
  widget_control, slider_A, set_value = A
  widget_control, slider_f, set_value = f
  widget_control, slider_particle_size, set_value = particle_size
  widget_control, slider_tao, set_value = tau
  widget_control, slider_portion, set_value = flat_portion
  pulse_simulation, base, parameters.A,parameters.portion,parameters.particle_size, parameters.f, parameters.tau
  xmanager, 'pulse_Excitation', base
end