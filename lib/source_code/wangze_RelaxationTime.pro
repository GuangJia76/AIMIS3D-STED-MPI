pro interpol

read,x_interpol,prompt='输入磁场强度(khz):'
read,y_interpol,prompt='输入粒子大小(mt):'


x = [4.5,9.3,12.2,25.0]
y = [5.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30]
z = DBLARR(4,11) * 0
z[0,*] = [7.80,6.80,6.10,5.75,5.20,5.05,4.70,4.40,4.25,4.05,3.90]
z[1,*] = [3.85,3.35,3.15,3.00,2.90,2.65,2.55,2.35,2.25,2.20,2.15]
z[2,*] = [2.45,2.35,2.30,2.10,2.05,1.95,1.90,1.80,1.65,1.55,1.65]
z[3,*] = [1.20,1.30,1.20,1.25,1.20,1.15,1.10,1.05,1.00,0.90,0.80]

g = DBLARR(4) * 0
g[0] = interpol(z[0,*],y,y_interpol)
g[1] = interpol(z[1,*],y,y_interpol)
g[2] = interpol(z[2,*],y,y_interpol)
g[3] = interpol(z[3,*],y,y_interpol)

a = interpol(g,x,x_interpol)

window, 1
plot, x, g, XRANGE=[0, 30], YRANGE=[0, 10], background='ffffff'x,color=0,XTITLE='Frequency [kHz]', YTITLE='Relaxation Time [us]', title = 'Resovist: '+STRING(a,format='(f10.4)')+'uT ('+STRING(y_interpol,format='(f10.4)') + 'mT, '+STRING(x_interpol,format='(f10.4)')+'kHz)'
FOR i = 0, 10 DO BEGIN
oplot, x, REFORM(z[*, i]),color=0
ENDFOR
oplot, x, g, color=0,PSYM=1
oplot, [x_interpol], [a], color=0,PSYM=2
;
window, 2
plot, y, z[0, *], XRANGE=[0, 40], YRANGE=[0, 10],background='ffffff'x,color=0, XTITLE='Drive Field [mT]', YTITLE='Relaxation Time [us]', title = 'Resovist: '+STRING(a,format='(f10.4)')+'uT ('+STRING(y_interpol,format='(f10.4)') + 'mT, '+STRING(x_interpol,format='(f10.4)')+'kHz)
oplot, y, z[1, *],color=0
oplot, y, z[2, *],color=0
oplot, y, z[3, *],color=0
oplot, [y_interpol], [g[0]],color=0, PSYM=1
oplot, [y_interpol], [g[1]],color=0, PSYM=1
oplot, [y_interpol], [g[2]],color=0, PSYM=1
oplot, [y_interpol], [g[3]],color=0, PSYM=1
oplot, [y_interpol], [a],color=0, PSYM=2

print , a
end




PRO idl_spec_sym

  A = [1,2,3,4,5,6]
  B = [1,2,3,8,5,6]
  C = DBLARR(6)*0.
  Resovist = 1
  FOR　i=0, 5 DO BEGIN
    C[i] = wangze(A[i], B[i], /Resovist)
  ENDFOR

  window, 1
  plot,A,C,TITLE='IDL',XTITLE='Amplitude',YTITLE='Relaxation',psym=4
  window, 2
  plot,B,C,TITLE='IDL',XTITLE='Frequency',YTITLE='Relaxation',psym=4
END


pro dd

  x = [1.3,2.6,5.1]
  y0 = [16,7.9,4.0]
  y1 = [6.6,3.3,1.7]
  y2 = [3.4,1.7,0.86]
  y3 = [2.0,1.0,0.5]
  y4 = [0.83,0.41,0.21]
  y5 = [0.42,0.21,0.11]

  spline_p, x, y0, xr, yr0
  spline_p, x, y1, xr1, yr1
  spline_p, x, y2, xr2, yr2
  spline_p, x, y3, xr3, yr3
  spline_p, x, y4, xr4, yr4
  spline_p, x, y5, xr5, yr5
  window , 1 , xsize = 400 ,ysize = 400
  plot,xr,yr0, XRANGE=[1, 5], YRANGE=[0, 10.0],background='ffffff'x,color=0, XTITLE='Gradient', YTITLE='Resolution', title = ''
  oplot,xr1,yr1,color=0
  oplot,xr2,yr2,color=0
  oplot,xr3,yr3,color=0
  oplot,xr4,yr4,color=0
  oplot,xr5,yr5,color=0

end


pro abc

  read,x_interpol,prompt='输入磁场强度(khz):'
  read,y_interpol,prompt='输入粒子大小(mt):'


  x = [4.5,9.3,12.2,25.0]
  y = [5.0,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30]
  z = DBLARR(4,11) * 0
  z[0,*] = [3.45,3.05,2.80,2.35,2.25,2.10,1.90,1.80,1.70,1.60,1.55]
  z[1,*] = [2.30,1.65,1.50,1.35,1.25,1.20,1.15,1.15,1.10,1.10,1.10]
  z[2,*] = [1.35,1.15,0.95,0.85,0.80,0.75,0.70,0.70,0.65,0.70,0.65]
  z[3,*] = [0.60,0.55,0.50,0.45,0.40,0.45,0.50,0.50,0.45,0.45,0.45]

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
end
