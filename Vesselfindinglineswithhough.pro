;  $Id: //depot/Release/ENVI53_IDL85/idl/idldir/examples/doc/image/findinglineswithhough.pro#1 $

;  Copyright (c) 2005-2015, Exelis Visual Information Solutions, Inc. All
;       rights reserved.
; 
PRO VesselFindingLinesWithHough

; Import the image from file.
;file = FILEPATH('rockland.png', $
;   SUBDIRECTORY = ['examples', 'data'])
;file = 'C:\D_drive\人民医院\pic1.png'
file = 'C:\D_drive\AIMIS_3D\images\test3.png'
image = READ_PNG(file)

; Determine size of image.
imageSize = SIZE(image, /DIMENSIONS)

; Initialize the TrueColor display.
DEVICE, DECOMPOSED = 1

;; Create a window and display the original image.
;WINDOW, 0, XSIZE = imageSize[1], YSIZE = imageSize[2], $
;   TITLE = 'Rockland, Maine'
;TV, image, TRUE = 1
;iimage, image

; Use the image from green channel to provide outlines
; of shapes.
intensity = REFORM(image[1, *, *])

; Determine size of intensity image.
intensitySize = SIZE(intensity, /DIMENSIONS)

; Mask intensity image to highlight power lines.
mask = intensity GT 150

; Initialize the remaining displays.
DEVICE, DECOMPOSED = 0
LOADCT, 0

;; Create another window and display the masked image.
;WINDOW, 1, XSIZE = intensitySize[0], $
;   YSIZE = intensitySize[1], $
;   TITLE = 'Mask to Locate Power Lines'
;TVSCL, mask
iimage, BYTSCL(mask)

; Transform mask.
transform = HOUGH(mask, RHO = rho, THETA = theta)

; Define the size and offset parameters for the
; transform displays.
displaySize = [256, 256]
offset = displaySize/3

; Reverse color table to clarify lines.
TVLCT, red, green, blue, /GET
TVLCT, 255 - red, 255 - green, 255 - blue

;; Create another window and display the Hough transform
;; with axes.
;WINDOW, 2, XSIZE = displaySize[0] + 1.5*offset[0], $
;   YSIZE = displaySize[1] + 1.5*offset[1], $
;   TITLE = 'Hough Transform'
;TVSCL, CONGRID(transform, displaySize[0], $
;   displaySize[1]), offset[0], offset[1]
;PLOT, theta, rho, /XSTYLE, /YSTYLE, $
;   TITLE = 'Hough Transform', XTITLE = 'Theta', $
;   YTITLE = 'Rho', /NODATA, /NOERASE, /DEVICE, $
;   POSITION = [offset[0], offset[1], $
;   displaySize[0] + offset[0], $
;   displaySize[1] + offset[1]], CHARSIZE = 1.5, $
;   COLOR = !P.BACKGROUND

;iimage, transform
; Scale transform to obtain just the power lines.
thres = MAX(transform)-40
transform = (TEMPORARY(transform) - thres) > 0
;iimage, transform
labelImg = LABEL_REGION(transform, ALL_NEIGHBORS=allNeighbors)
;find the maximum number of connected reions
labelmax=max(labelImg)
h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
sort_re=REVERSE(SORT(h))
print,'max', max(h)
hmax=n_elements(sort)
maskend=transform*0
back_maskend = mask*0
; Determine size of label image.
labelImgSize = SIZE(labelImg, /DIMENSIONS)
N_ele=labelmax;15
if labelmax ge N_ele then begin
  n_grows=N_ele
  for j=1,N_ele do begin
    maskendmid = labelImg EQ sort_re[j]
    maskindex = WHERE(labelImg EQ sort_re[j], ncount)
    maskindex2d = ARRAY_INDICES(labelImg, maskindex)
    x_index = MEDIAN(maskindex2d[0,*])
    y_index = MEDIAN(maskindex2d[1,*])
    temp = maskendmid*0.
    temp[x_index, y_index] = 1
    backprojectionmid = HOUGH(maskendmid*temp, /BACKPROJECT, $
      RHO = rho, THETA = theta, $
      NX = intensitySize[0], NY = intensitySize[1])
    if j EQ 6 THEN BEGIN
;      iimage, maskendmid
      ;iimage, mask*backprojectionmid
      backprojectionmid[intensitySize[0]/2, intensitySize[1]/2] = 255
      iimage, backprojectionmid + maskend
      ;determin the line's rho and theta
      theta_x_index = MEAN(maskindex2d[0,*])/labelImgSize[0] * 180; * !PI
      rho_y_index = (0.5 - MEAN(maskindex2d[1,*])/labelImgSize[1]) * 2 * MAX(rho)
      print, 'theta =', theta_x_index, ' degree, rho =', rho_y_index

      maskend+=maskendmid*j
    endif
    back_maskend = back_maskend + mask*backprojectionmid
;    maskend+=maskendmid*j
  endfor
endif

iimage, back_maskend
;iimage, maskend

; Create another window and display the scaled transform.
WINDOW, 3, XSIZE = displaySize[0] + 1.5*offset[0], $
   YSIZE = displaySize[1] + 1.5*offset[1], $
   TITLE = 'Scaled Hough Transform'
TVSCL, REVERSE(CONGRID(maskend, displaySize[0], $
   displaySize[1]),2), offset[0], offset[1]
PLOT, theta*180./!PI, rho, /XSTYLE, /YSTYLE, $
   TITLE = 'Scaled Hough Transform', XTITLE = 'Theta', $
   YTITLE = 'Rho', /NODATA, /NOERASE, /DEVICE, $
   POSITION = [offset[0], offset[1], $
   displaySize[0] + offset[0], $
   displaySize[1] + offset[1]], CHARSIZE = 1.5, $
   COLOR = !P.BACKGROUND

; Backproject to compare with original image.
;backprojection = HOUGH(transform, /BACKPROJECT, $
backprojection = HOUGH(maskend, /BACKPROJECT, $
   RHO = rho, THETA = theta, $
   NX = intensitySize[0], NY = intensitySize[1])

;; Create another window and display the results.
;WINDOW, 4, XSIZE = intensitySize[0], $
;   YSIZE = intensitySize[1], $
;   TITLE = 'Resulting Power Lines'
;TVSCL, backprojection

;iimage, backprojection

END