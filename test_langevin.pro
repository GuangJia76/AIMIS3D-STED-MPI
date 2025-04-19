PRO abc

print, 'hello new world!'

x = (FINDGEN(401)-200.)/10.
y = 1./tanh(x) - 1./x
y[200] = 0.

iplot, x, y

END

PRO test_Langevin

;print, 'hello old world!'
print, 'calculate Langevin'

x = (FINDGEN(401)-200.)/10.

y1 = 1./tanh(x)

iplot, x, y1

y2 = -1./x

iplot, x, y2

y = y1 + y2

y[200] = 0.

iplot, x, y

END