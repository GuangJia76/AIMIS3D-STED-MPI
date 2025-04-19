; docformat = 'rst'
;+
; This is an example program to demonstrate how to create a zonal wind plot
; with Coyote Graphics routines.
;
; :Categories:
;    Graphics
;    
; :Examples:
;    Save the program as "zonal_wind_plot.pro" and run it like this::
;       IDL> .RUN zonal_wind_plot
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
;        Written, 10 February 2013 by David W. Fanning.
;
; :Copyright:
;     Copyright (c) 2013, Fanning Software Consulting, Inc.
;-
PRO Zonal_Wind_Plot
 
   ; Set up variables for the plot. Normally, these values would be 
   ; passed into the program as positional and keyword parameters.
   ; Data file obtained here: http://www.ncl.ucar.edu/Document/Manuals/Getting_Started/Examples/gsun11n.shtml.
   file = 'C:\D_drive\AIMIS_3D\coyote\zonal_wind_plot\cocos_island.dat'
   
   ; Read the ASCII data file. First column is pressure. Second column is height.
   ; Next 12 columns are zonal winds for each month.
   OpenR, lun, file, /Get_Lun
   lines = File_Lines(file)
   data = Fltarr(14, lines)
   ReadF, lun, data
   Free_Lun, lun
   
   ; Reformat pressure and height as row vectors.
   pressure = Reform(data[0,*])
   height = Reform(data[1,*])
   
   ; Get dimension of wind data set.
   monthlyWinds = data[2:*,*]
   dims = Size(monthlyWinds, /Dimensions)
   
   ; Create a month vector for plotting purposes
   time = Indgen(12)
   tick_names = cgMonths(time+1, /FirstLetter, /Abbreviation)

   ; Scale the height for drawing zonal winds and load program colors.
   heightColors = BytScl(height, Min=0, Max=23, Top=23)
   cgLoadCT, 33, NColors=24

   ; Open a display window.
   cgDisplay

   ; Set up plotting parameters.
   thick = (!D.Name EQ 'PS') ? 6 : 2
   plotPosition = [0.125, 0.125, 0.9, 0.75]
   
   ; Draw a gray background so colors stand out more.
   cgColorFill, Position=plotPosition, Color='gray'
   
   ; Draw the initial plot, no data.
   cgPlot, time, monthlyWinds[*, 0], YRange=[-20,10], YStyle=1, $
       YTitle='Amplitude (m/s)', XTitle='Month', XStyle=1, $
       XTicks=11, XTickName=tick_names, /NoData, /NoErase, $
       Position=plotPosition, Label='Zonal Winds Cocos Island'
       
   ; Overplot the zonal winds in colors scaled by height.
   FOR j=0,dims[1]-1 DO cgPlotS, time, monthlyWinds[*,j], time, $
       Color=heightColors[j], Thick=thick
       
   ; Add a color bar.
   cgColorbar, Divisions=24, NColors=24, Range=[0,24], Title='Height (km)', $
       TLocation='top', /Discrete, Format='(i0)', Position=[0.125, 0.86, 0.9, 0.90]
       
END ;*****************************************************************

; This main program shows how to call the program and produce
; various types of output.

  ; Display the plot in a graphics window.
  Zonal_Wind_Plot
  
  ; Display the plot in a resizeable graphics window.
  cgWindow, 'Zonal_Wind_Plot', WBackground='White', $
     WTitle='Zonal Wind Plot in a Resizeable Graphics Window'
  
  ; Create a PostScript file.
  cgPS_Open, 'zonal_wind_plot.ps'
  Zonal_Wind_Plot
  cgPS_Close
  
  ; Create a PNG file with a width of 600 pixels.
  cgPS2Raster, 'zonal_wind_plot.ps', /PNG, Width=600

END