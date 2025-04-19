
pro hk_simulation,wid,A

  field_amp = DBLARR(201,201)*0. ; 100 rows, 201 columns
  x_field = DBLARR(201,201)*0. ; 100 rows, 201 columns
  y_field = DBLARR(201,201)*0. ; 100 rows, 201 columns

  coeff=100.
  FOR i=0, 200 DO BEGIN
    FOR j=0, 200 DO BEGIN
      dist = SQRT((i - 100.)^2 + (j - 100.)^2)
      field_amp[i,j] = 1./dist;*coeff
      x_field[i,j] = field_amp[i,j] * (j - 100.)
      y_field[i,j] = field_amp[i,j] * (i - 100.)
    ENDFOR
  ENDFOR

  sum_x_field = DBLARR(100,100)*0.
  sum_y_field = DBLARR(100,100)*0.
  sum_field_amp = DBLARR(100,100)*0.
  sum_x_field = x_field[101:200, 0:99] - A * x_field[0:99, 0:99] + x_field[0:99, 101:200] - x_field[101:200, 101:200];按钮转动变换符号，电流方向
  sum_y_field = y_field[101:200, 0:99] - A * y_field[0:99, 0:99] + y_field[0:99, 101:200] - y_field[101:200, 101:200]
  sum_field_amp = SQRT(sum_x_field^2 + sum_y_field^2)

  
  draw1 = widget_info(wid, find_by_uname='draw1')
  widget_control, draw1, get_value = drawWindow1
  drawSize = widget_info(draw1, /GEOMETRY)

  wset, drawWindow1
 contour, sum_field_amp,  levels=FINDGEN(100)/99.*max(sum_field_amp),title='Field Amplitude';整体场强分布图
 
 ;iimage, sum_x_field,title ='x-direction Field of Anderson coil '
 ;iimage, sum_y_field,title ='y-direction Field of Anderson coil '
  plot ,sum_field_amp[50,*], title ='field Amplitude on y-axis
 plot ,sum_field_amp[*,50], title ='field Amplitude on x-axis'
 
 !P.Multi = [0,3,1]

end

pro slider_A_hk_event, ev
  widget_control, ev.id, get_value = A
  widget_control, ev.top, get_uvalue = parameters
  parameters.A = A;
  widget_control, ev.top, set_uvalue = parameters

  hk_simulation,ev.top, parameters.A
end

pro hk_Anderson_coil

  ;��ʼֵ
  A = 1
  parameters = {A:A}
  
  device, decomposed=1
  device, get_screen_size = screenSize
  baseXSize = screenSize[0] * 2.5
  baseYSize = screenSize[1] * 1

  base = widget_base(title='Anderson Coil', xsize = baseXSize/2.5, ysize = baseYSize, /column)

  menubase = widget_base(base, xsize = baseXSize/3, /row)

  sliderbase = widget_base(menubase, /column)
  ;base_A
  base_A = widget_base(sliderbase, /column)
  label_A = widget_label(base_A, value = "The magnitude of the current, negative represent the direction of the current:")
  slider_A_hk = widget_slider(base_A, xsize=400, MINIMUM=-10,MAXIMUM=10 ,event_pro = 'slider_A_hk_event')
  
  picbase = widget_base(menubase, /column)
  MPIjpg_path = './lib/icon_pics/MPI_coil.jpg'
  READ_JPEG, MPIjpg_path, MPIjpg, true=1
  picwindow = WIDGET_draw(picbase, uvalue='pic_win', xsize = (size(MPIjpg))[2], ysize = (size(MPIjpg))[3])

  
  ;��ͼ��
  drawXSize = baseXSize / 3.1
  drawYSize = baseYSize / 3
  drawbase = widget_base(base, /column)
  drawbase1 = widget_base(drawbase, /row)
  draw1 = widget_draw(drawbase1, xsize = drawXSize, ysize = drawYSize, retain=2, uname='draw1')



  widget_control, base, /realize
  
  widget_control, picwindow, get_value = drawPic
  wset, drawPic
  TVSCL, MPIjpg, true=1

  widget_control, base, set_uvalue = parameters

  ;����slider�����value��ʼֵ
  widget_control, slider_A_hk, set_value = A

  hk_simulation, base, parameters.A
  xmanager, 'Excitation', base
end
