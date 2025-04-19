; docformat = 'rst'
;+
; This is an example program to demontrate how to create a filled contour plot
; with Coyote Graphics routines.
;
; :Categories:
;    Graphics
;
; :Examples:
;    Save the program as "filled_contour_plot_1.pro" and run it like this::
;       IDL> .RUN filled_contour_plot_1
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
;        Written, 23 January 2013 by David W. Fanning.
;
; :Copyright:
;     Copyright (c) 2013, Fanning Software Consulting, Inc.
;-
PRO Filled_Contour_Plot_1

  ; Example Gaussian data.
  data = cgDemoData(26)

  ; Set up variables for the contour plot. Normally, these values
  ; would be passed into the program as positional and keyword parameters.
  minValue = Floor(Min(data))
  maxValue = Ceil(Max(data))
  nLevels = 10
  xtitle = 'X Axis'
  ytitle = 'Y Axis'
  position =   [0.125, 0.125, 0.9, 0.800]
  cbposition = [0.125, 0.865, 0.9, 0.895]
  cbTitle = 'Data Value'

  ; Set up a "window" for the plot. The PostScript output will have
  ; the same aspect ratio as the graphics window on the display.
  cgDisplay, 600, 500, Title='Filled Contour Plot 1'

  ; Set up colors for contour plot.
  cgLoadCT, 33, NColors=nlevels, Bottom=1, CLIP=[30,255]

  ; Draw the filled contour plot.
  contourLevels = cgConLevels(data, NLevels=10, MinValue=minValue)
  cgContour, data, /Fill, Levels=contourLevels, C_Colors=Bindgen(nLevels)+1B, $
    /OutLine, Position=position, XTitle=xtitle, YTitle=ytitle

  ; Draw the color bar.
  cgColorbar, NColors=nlevels, Bottom=1, Position=cbposition, $
    Range=[MinValue, MaxValue], Divisions=nlevels, /Discrete, $
    Title=cbTitle, TLocation='Top'

END ;*****************************************************************

; This main program shows how to call the program and produce
; various types of output.

; Display the contour plot in a graphics window.
Filled_Contour_Plot_1

; Display the contour plot in a resizeable graphics window.
cgWindow, 'Filled_Contour_Plot_1', WXSize=600, WYSize=500, $
  WTitle='Filled Contour Plot in Resizeable Graphics Window'

; Create a PostScript file.
cgPS_Open, 'filled_contour_plot_1.ps'
Filled_Contour_Plot_1
cgPS_Close

; Create a PNG file with a width of 600 pixels.
cgPS2Raster, 'filled_contour_plot_1.ps', /PNG, Width=600

END