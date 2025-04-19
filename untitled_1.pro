PRO test_aa

  a = 2.
  b = 6.
  c = 8./25.
  
  x = FINDGEN(1000.)-500.
  
  y = a * x^2 + b * x + c
  
  iplot, x, y, color='red'


END