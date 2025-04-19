pro final_act

  read,x_interpol,prompt='输入磁场强度(khz):'
  read,y_interpol,prompt='输入粒子大小(mt):'

  a=''
  read,a,prompt='粒子类型'
  R="resovist"
  U="uw33"
  a=strlowcase(a)
  c=strcmp(a,r)
  d=strcmp(a,u)


  x = [4.5,9.3,12.2,25.0]
  y = [5.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30]

  if c=1
  then begin
    z = DBLARR(4,11) * 0
    z[0,*] = [7.80,6.80,6.10,5.75,5.20,5.05,4.70,4.40,4.25,4.05,3.90]
    z[1,*] = [3.85,3.35,3.15,3.00,2.90,2.65,2.55,2.35,2.25,2.20,2.15]
    z[2,*] = [2.45,2.35,2.30,2.10,2.05,1.95,1.90,1.80,1.65,1.55,1.65]
    z[3,*] = [1.20,1.30,1.20,1.25,1.20,1.15,1.10,1.05,1.00,0.90,0.80]
  end if
  if d=1
  then begin
    z = DBLARR(4,11) * 0
    z[0,*] = [3.45,3.05,2.80,2.35,2.25,2.10,1.90,1.80,1.70,1.60,1.55]
    z[1,*] = [2.30,1.65,1.50,1.35,1.25,1.20,1.15,1.15,1.10,1.10,1.10]
    z[2,*] = [1.35,1.15,0.95,0.85,0.80,0.75,0.70,0.70,0.65,0.70,0.65]
    z[3,*] = [0.60,0.55,0.50,0.45,0.40,0.45,0.50,0.50,0.45,0.45,0.45]
  endif

  g = DBLARR(4) * 0
  g[0] = interpol(z[0,*],y,y_interpol)
  g[1] = interpol(z[1,*],y,y_interpol)
  g[2] = interpol(z[2,*],y,y_interpol)
  g[3] = interpol(z[3,*],y,y_interpol)

  a = interpol(g,x,x_interpol)

  window, 1
  plot, x, g, XRANGE=[0, 30], YRANGE=[0, 10], background='ffffff'x,color=0,XTITLE='Frequency [kHz]', YTITLE='Relaxation Time [us]', title = 'uw33: '+STRING(a,format='(f10.4)')+'uT ('+STRING(y_interpol,format='(f10.4)') + 'mT, '+STRING(x_interpol,format='(f10.4)')+'kHz)'
  FOR i = 0, 10 DO BEGIN
    oplot, x, REFORM(z[*, i]),color=0
  ENDFOR
  oplot, x, g, color=0,PSYM=1
  oplot, [x_interpol], [a], color=0,PSYM=2

  window, 2
  plot, y, z[0, *], XRANGE=[0, 40], YRANGE=[0, 10],background='ffffff'x,color=0, XTITLE='Drive Field [mT]', YTITLE='Relaxation Time [us]', title = 'uw33: '+STRING(a,format='(f10.4)')+'uT ('+STRING(y_interpol,format='(f10.4)') + 'mT, '+STRING(x_interpol,format='(f10.4)')+'kHz)
  oplot, y, z[1, *],color=0
  oplot, y, z[2, *],color=0
  oplot, y, z[3, *],color=0
  oplot, [y_interpol], [g[0]],color=0, PSYM=1
  oplot, [y_interpol], [g[1]],color=0, PSYM=1
  oplot, [y_interpol], [g[2]],color=0, PSYM=1
  oplot, [y_interpol], [g[3]],color=0, PSYM=1
  oplot, [y_interpol], [a],color=0, PSYM=2

  print , a

  result = SFIT(z, 3)
  WINDOW,3, XSIZE = 800, YSIZE = 400
  !P.MULTI = [0, 2, 1]
  DEVICE, DECOMPOSED=0
  !P.BACKGROUND = 255
  !P.COLOR = 0
  SURFACE, z, X, Y
  SURFACE, result, X, Y


end