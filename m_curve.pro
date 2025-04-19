PRO M_curve

t = FINDGEN(1001) ;us
f = 1.; kHz
H = 10.*sin(t/1000. * f * 2.* !PI)

M = 1./tanh(H) - 1./H
iplot, t, H
iplot, t, M, color='red', /overplot

iplot, H, M, title='M-H curve', xtitle='H', ytitle='M'

signal = -(M[1:*] - M[0:999])/(t[1]-t[0])

iplot, t, signal
END