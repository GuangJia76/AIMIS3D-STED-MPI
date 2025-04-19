pro slider_event, ev
  widget_control, ev.id, get_value = a
  widget_control, ev.id, get_uvalue = info
  print, a 
  text = widget_info(ev.top, find_by_uname = 'text')
  widget_control, text, set_value = string(a + info.b + info.c)
  
  draw = widget_info(ev.top, find_by_uname = 'draw')
  widget_control, draw, get_value = window
  wset, window
  X = findgen(100)/10
  Y = sin(X)*a
  plot, X, Y
  
  draw2 = widget_info(ev.top, find_by_uname = 'draw2')
  widget_control, draw2, get_value = window2
  wset, window2
  H = 1/tanh(Y) - 1/Y
  plot, X, H
end

pro test
  b = 10
  c = 20
  info = {b:b, c:c}
  base = widget_base(xsize = 1000, ysize = 1000, title = 'abc', /column)
  
    slider = widget_slider(base, event_pro = 'slider_event', uvalue = info)
    text = widget_text(base, uname = 'text', /editable)
    draw = widget_draw(base, xsize = 400, ysize = 300, uname = 'draw')
    draw2 = widget_draw(base, xsize = 400, ysize = 300, uname = 'draw2')
  widget_control, base, /realize
  xmanager, 'test', base
end