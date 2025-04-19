;+
; NAME: struc_vlineobject__define
;
;
;
; PURPOSE: define the "struc_vlineobject_info" structure. Define vessel line object for automatic brain vessel
; segmentation based on anatomical structure, A1, A2, M1, P1, etc.
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
;; GJ, 2020/6/16
;
PRO  struc_vlineobject__define
  
  s ={STRUC_VLINEOBJECT,      $
    ID:0L, $;the ID of vessel line segment
    sort_ID:0L, $;the ID after sorting
    rho:0L, $;the angle of ρ is the perpendicular distance from the origin of a line at an angle θ
    theta:0L, $;the angle θ will be limited to 0 <= θ < π which could result in negative ρ values.
    length:0L, $;the length of the line segment
    midpoint_loc:[0,0], $;the location of the middle point
    endpoint1_loc:[0,0], $;the location of the 1st endpoint
    endpoint1_Ncross:0L, $;the number of cross lines at 1st endpoint
    endpoint1_Mcross:DBLARR(10, 4), $;max 10 cross lines, including ID, angle, location x, location y 
    endpoint2_loc:[0,0], $;the location of the 2nd endpoint
    endpoint2_Ncross:0L, $;the number of cross lines at 1st endpoint
    endpoint2_Mcross:DBLARR(10, 4) $;max 10 cross lines, including ID, angle, location x, location y
  }
END

pro struc_vlineobject::info
  self_tag = tag_names({struc_vlineobject})
  for i= 0, n_elements(self_tag) - 1 do print,self_tag(i),' = ',self.(i)
end ; of method struc_vlineobject.info

;**********************************************************************************

pro struc_vlineobject::set,_extra = extra

  self_tag = tag_names({struc_vlineobject})
  extra_tag = tag_names(extra)

  for i= 0, n_elements(self_tag) - 1 do begin
    j= 0
    repeat begin
      treffer = self_tag(i) eq extra_tag(j)
      if treffer then self.(i) = extra.(j) else j = j+1
    endrep until treffer or (j ge n_elements(extra_tag))
  endfor
end ; of method struc_vlineobject::set

;**********************************************************************************
;finding vessel line using hough transformation
;GJ, 2020/6/16
;
PRO vlines_hough_transformj, mask, vline_info, back_maskend;_s

  ; Import the image from file.
  ;file = FILEPATH('rockland.png', $
  ;   SUBDIRECTORY = ['examples', 'data'])
  ;file = 'C:\D_drive\人民医院\pic1.png'

;  file = 'C:\D_drive\AIMIS_3D\images\test5.png';4.png';3.png'
;  image = READ_PNG(file)

  ; Determine size of image.
  imageSize = SIZE(mask, /DIMENSIONS)
  
  ;2020/7/14, debugging
  ;iimage, mask
;  ; Initialize the TrueColor display.
;  DEVICE, DECOMPOSED = 1

;  ;; Create a window and display the original image.
;  WINDOW, 0, XSIZE = imageSize[0], YSIZE = imageSize[1], $
;     TITLE = 'MIP of brain vessel'
;  TV, REVERSE(BYTSCL(mask),2);, TRUE = 1
;  ;iimage, image

;  ; Use the image from green channel to provide outlines
;  ; of shapes.
;  intensity = REFORM(image[1, *, *])

  ; Determine size of intensity image.
  intensitySize = SIZE(mask, /DIMENSIONS)

;  iimage, BYTSCL(mask)

  ; Transform mask.
  transform = HOUGH(mask, RHO = rho, THETA = theta)

  ; Define the size and offset parameters for the
  ; transform displays.
  displaySize = [256, 256]
  offset = displaySize/3

;  ; Reverse color table to clarify lines.
;  TVLCT, red, green, blue, /GET
;  TVLCT, 255 - red, 255 - green, 255 - blue

;  iimage, transform;BYTSCL(transform)
  ; Scale transform to obtain just the power lines.
  thres = MAX(transform)/2;-23;23;60;50;40
  print, 'threshold = ', thres
  print, 'max = ', MAX(transform)
  transform = (TEMPORARY(transform) - thres) > 0
  
;  ;2020/7/14, check whether the threshold affect
;  backprojectionmid_trans = HOUGH(transform, /BACKPROJECT, $
;    RHO = rho, THETA = theta, $
;    NX = intensitySize[0], NY = intensitySize[1])
;  iimage, BYTSCL(backprojectionmid_trans)
  
  ; Create another window and display the Hough transform
  ; with axes.
  WINDOW, 1, XSIZE = displaySize[0] + 1.5*offset[0], $
    YSIZE = displaySize[1] + 1.5*offset[1], $
    TITLE = 'Hough Transform'
  TVSCL, CONGRID(transform, displaySize[0], $
    displaySize[1]), offset[0], offset[1]
  PLOT, theta, rho, /XSTYLE, /YSTYLE, $
    TITLE = 'Hough Transform', XTITLE = 'Theta', $
    YTITLE = 'Rho', /NODATA, /NOERASE, /DEVICE, $
    POSITION = [offset[0], offset[1], $
    displaySize[0] + offset[0], $
    displaySize[1] + offset[1]], CHARSIZE = 1.5, $
    COLOR = 255;!P.BACKGROUND
  
  labelImg = LABEL_REGION(transform, ALL_NEIGHBORS=allNeighbors)
  ;find the maximum number of connected reions
  labelmax=max(labelImg)
  ;define the vessel line information structure
  vline_info = replicate({struc_vlineobject},labelmax)
  
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
  sort_re=REVERSE(SORT(h))
  print,'max', max(h)
  hmax=n_elements(sort_re)
  maskend=transform*0
  back_maskend = mask*0
  back_maskend_s = mask*0
  ; Determine size of label image.
  labelImgSize = SIZE(labelImg, /DIMENSIONS)
  N_ele=labelmax;15
  if labelmax ge N_ele then begin
    n_grows=N_ele
    for j=1,N_ele do begin
      maskendmid = labelImg EQ sort_re[j]
      maskindex = WHERE(labelImg EQ sort_re[j], ncount)
      maskindex2d = ARRAY_INDICES(labelImg, maskindex)
      
      ;to find the rho and theta based on the central point
      x_index = MEDIAN(maskindex2d[0,*]) ;define the mid central point location
      y_index = MEDIAN(maskindex2d[1,*])
      temp = maskendmid*0.
      temp[x_index, y_index] = 1 ;define the central point
      ;a central point's hough transform
      backprojectionmid_s = HOUGH(maskendmid*temp, /BACKPROJECT, $
        RHO = rho, THETA = theta, $
        NX = intensitySize[0], NY = intensitySize[1])
;      if j EQ 6 THEN BEGIN
        ;      iimage, maskendmid
        ;iimage, mask*backprojectionmid
        backprojectionmid_s[intensitySize[0]/2, intensitySize[1]/2] = 255
        ;iimage, backprojectionmid_s + maskend
        ;determin the line's rho and theta
        theta_x_index = MEAN(maskindex2d[0,*])/labelImgSize[0] * 180; * !PI
        rho_y_index = (0.5 - MEAN(maskindex2d[1,*])/labelImgSize[1]) * 2 * MAX(rho)
        vline_info[j-1].ID = j
        vline_info[j-1].sort_ID = 0 ;has not been sorted
        vline_info[j-1].rho = rho_y_index
        vline_info[j-1].theta = theta_x_index
        print, 'theta =', theta_x_index, ' degree, rho =', rho_y_index
        
        ;find the line segment
        backprojectionmid = HOUGH(maskendmid, /BACKPROJECT, $
          RHO = rho, THETA = theta, $
          NX = intensitySize[0], NY = intensitySize[1])
        vline_mask = mask*backprojectionmid
;        iimage, vline_mask
        label_vlines = LABEL_REGION(vline_mask, ALL_NEIGHBORS=allNeighbors)
        ;find the maximum number of connected reions
        label_vlines_max=max(label_vlines)
        vline_h = HISTOGRAM(label_vlines, REVERSE_INDICES=r)
        vline_sort_re=REVERSE(SORT(vline_h))
        print,'max', max(vline_h)
        vline_hmax=n_elements(vline_sort_re)
        vline_maskendmid = label_vlines EQ vline_sort_re[1]
         
        vline_segment_ind = WHERE(vline_maskendmid, ncount_vline)
        vline_segment2d = ARRAY_INDICES(mask, vline_segment_ind)
       
        ;determine the line's length and endpoint locations
        x_length = MAX(vline_segment2d[0,*], max_x_ind)-MIN(vline_segment2d[0,*], min_x_ind)
        y_length = MAX(vline_segment2d[1,*], max_y_ind)-MIN(vline_segment2d[1,*], min_y_ind)
        vline_info[j-1].length = SQRT(x_length^2 + y_length^2)
        vline_info[j-1].midpoint_loc = [MEAN(vline_segment2d[0,*]), MEAN(vline_segment2d[1,*])]
        IF min_x_ind LT max_x_ind THEN BEGIN
          vline_info[j-1].endpoint1_loc = [vline_segment2d[*, min_x_ind]]
          vline_info[j-1].endpoint2_loc = [vline_segment2d[*, max_x_ind]]
        ENDIF ELSE BEGIN
          vline_info[j-1].endpoint1_loc = [vline_segment2d[*, max_x_ind]]
          vline_info[j-1].endpoint2_loc = [vline_segment2d[*, min_x_ind]]
        ENDELSE
;        help, vline_info[j-1]
;        maskend+=maskendmid*j
;      endif
      back_maskend_s = back_maskend_s + (mask GT 0)*(backprojectionmid_s GT 0)*j
      back_maskend = back_maskend + (mask GT 0)*(backprojectionmid GT 0)
      ;    maskend+=maskendmid*j
    endfor
  endif
  
;  iimage, back_maskend
;  iimage, back_maskend_s
  print, 'rho = ', vline_info.rho
  print, 'theta = ', vline_info.theta
  print, 'length = ', vline_info.length
  ;GJ 2020/7/5, plot the midpoint
  ; Create another window and display the vessel segment.
  WINDOW, 2, TITLE = 'vessel line midpoints'
  plot, vline_info.midpoint_loc[0,*], vline_info.midpoint_loc[1,*], PSYM = 1
;  flag=0
  FOR i=0, N_ele-1 DO BEGIN
    IF vline_info[i].length GT 1 THEN BEGIN
;      IF flag EQ 0 THEN BEGIN
        oplot, vline_info[i].midpoint_loc[0,*], vline_info[i].midpoint_loc[1,*], PSYM = 4
;        flag=1
;      ENDIF ELSE BEGIN
;        oplot, vline_info[i].midpoint_loc[0,*], vline_info[i].midpoint_loc[1,*], PSYM = 4
;      ENDELSE       
    ENDIF
  ENDFOR
  
  ;GJ 2020/7/9, straight line commectivity
  sort_ID = 1; set the first sort ID as 1
  ;start sorting based on staight line connectivity
  FOR i=0, N_ele-1 DO BEGIN
    IF vline_info[i].sort_ID EQ 0 THEN BEGIN
      vline_info[i].sort_ID = sort_ID
      ;change the ID to sort_ID on mask
      imask_index = WHERE(back_maskend EQ vline_info[i].ID, ncount)
      IF ncount GT 0 THEN back_maskend[imask_index] = vline_info[i].sort_ID
      ;self increment by 1
      sort_ID++;
    ENDIF
    average_length = MEAN(vline_info.length)
    IF vline_info[i].length GT average_length THEN BEGIN
      line_connectivity = BYTARR(N_ele)*0.
      FOR j=0, N_ele-1 DO BEGIN
        ;一键三联GJ 2020/7/9
        mid_point_mid = (vline_info[i].midpoint_loc + vline_info[j].midpoint_loc)/2
        mid_point_mid_l = (vline_info[i].midpoint_loc + mid_point_mid)/2
        mid_point_mid_2 = (vline_info[j].midpoint_loc + mid_point_mid)/2
        IF (mask[mid_point_mid[0], mid_point_mid[1]] GT 0) AND (mask[mid_point_mid_l[0], mid_point_mid_l[1]] GT 0) AND (mask[mid_point_mid_2[0], mid_point_mid_2[1]] GT 0) THEN BEGIN
          line_connectivity[j] = 1 ;straight line is connected
        ENDIF
      ENDFOR
      
      ;find the longest length
      line_length = line_connectivity*vline_info.length
      max_length = MAX(line_length, max_length_ind)
      key_ID = vline_info[max_length_ind].ID
      FOR j=0, N_ele-1 DO BEGIN
        IF line_connectivity[j] EQ 1 THEN BEGIN
          vline_info[j].sort_ID = key_ID
          mask_index = WHERE(back_maskend EQ vline_info[j].ID, ncount)
          IF ncount GT 0 THEN back_maskend[mask_index] = key_ID
        ENDIF
      ENDFOR
    ENDIF
  ENDFOR
;  iimage, mask
  
;  iimage, back_maskend
  ;GJ, 2020/7/10, check single pixel's value, which should equal to the value of non-zero 8-nearest pixels
;  FOR i=0,10 DO BEGIN
    maskindex = WHERE(back_maskend GT 0, ncount)
    IF ncount GT 0 THEN BEGIN
       temp_back_maskend = back_maskend*0.
       maskindex2d = ARRAY_INDICES(back_maskend, maskindex)
       FOR i=0, ncount-1 DO BEGIN
         nearest = back_maskend[(maskindex2d[0, i]-1):(maskindex2d[0, i]+1), (maskindex2d[1, i]-1):(maskindex2d[1, i]+1)]
         value = nearest[1,1]
         nearest[1,1] = 0
         ;it is not a single pixel
         IF MAX(nearest) GT 0 THEN BEGIN
           pdf = histogram(nearest, MIN=1, LOCATIONS=xbin)
           maxpdf = MAX(pdf, maxpdf_ind)
           temp_back_maskend[maskindex2d[0, i], maskindex2d[1, i]] = xbin[maxpdf_ind]
         ENDIF ELSE BEGIN
           ;;it is not a single pixel
           temp_back_maskend[maskindex2d[0, i], maskindex2d[1, i]] = 0
         ENDELSE
       ENDFOR
       back_maskend = temp_back_maskend
    ENDIF
    
    maskindex = WHERE(back_maskend GT 0, ncount)
    IF ncount GT 0 THEN BEGIN
      temp_back_maskend = back_maskend*0.
      maskindex2d = ARRAY_INDICES(back_maskend, maskindex)
      FOR i=0, ncount-1 DO BEGIN
        nearest = back_maskend[(maskindex2d[0, i]-1):(maskindex2d[0, i]+1), (maskindex2d[1, i]-1):(maskindex2d[1, i]+1)]
        value = nearest[1,1]
        pdf = histogram(nearest, MIN=1, LOCATIONS=xbin)
        IF MAX(nearest) GT 0 AND pdf[value-1] EQ 1 THEN BEGIN
          maxpdf = MAX(pdf, maxpdf_ind)
          temp_back_maskend[maskindex2d[0, i], maskindex2d[1, i]] = xbin[maxpdf_ind]
        ENDIF
     ENDFOR
     
     ;count the fixed pixels
     temp_maskindex = WHERE(temp_back_maskend GT 0, temp_ncount)
     IF temp_ncount GT 0 THEN BEGIN
       FOR j=0, temp_ncount-1 DO BEGIN
        back_maskend[temp_maskindex[j]] = temp_back_maskend[temp_maskindex[j]]
       ENDFOR
     ENDIF
    ENDIF
;  ENDFOR

  
  ;separate the un-connected regions
  maskindex = WHERE(back_maskend GT 0, ncount)
  max_ind = MAX(back_maskend)+1
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      temp_back_maskend = back_maskend EQ i
      labelImg = LABEL_REGION(temp_back_maskend, ALL_NEIGHBORS=allNeighbors)
      IF MAX(labelImg) GT 1 THEN BEGIN
        FOR j=2,MAX(labelImg) DO BEGIN
          labelind = WHERE(labelImg EQ j, jcount)
          IF jcount GT 0 THEN BEGIN
            back_maskend[labelind] = max_ind
            max_ind++
          ENDIF
        ENDFOR 
      ENDIF
    ENDFOR
  ENDIF
  
;  iimage, back_maskend
  
  ;2020/7/12; A2 segmentation
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF incount LT 6 THEN BEGIN
          ;less than 45 degree
          IF (ymid-127 GT 0) AND (ABS(xmid-127) LT ABS(ymid-127)*0.5) THEN back_maskend[imaskindex]=100
        ENDIF ELSE BEGIN
          xrange = MAX(imaskindex2d[0,*]) - MIN(imaskindex2d[0,*])
          yrange = MAX(imaskindex2d[1,*]) - MIN(imaskindex2d[1,*])
          ;less than 11.3degree, y>5*x
          IF (yrange GT 2.*xrange) AND (ABS(xmid-127) LT ABS(ymid-127)) AND (ymid-127 GT 0) THEN back_maskend[imaskindex]=100
        ENDELSE
      ENDIF
    ENDFOR
  ENDIF
  ;2020/7/13, A2 segmentation
  ;100 is for A2 segmentation
  a2_index = WHERE(back_maskend EQ 100, a2_ncount)
  IF a2_ncount GT 0 THEN BEGIN
    a2_index2d = ARRAY_INDICES(back_maskend, a2_index)
  ENDIF ELSE BEGIN
    ;a2 could not be found
    ;error handling
    a2_index2d = [127,127]
  ENDELSE
  ;calculate the lower-right corner coordinates
  x0 = MAX(a2_index2d[0, *])
  y0 = MIN(a2_index2d[1, *])
  ;start resorting A2
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF (ymid GT y0) AND (ABS(ymid-y0) GT ABS(xmid-x0)) THEN back_maskend[imaskindex]=100
      ENDIF
    ENDFOR
  ENDIF
  
  
  ;2020/7/12; Right A1 segmentation
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF incount LT 4 THEN BEGIN
          IF (ABS(ymid-127) LT 20) AND (xmid-127 GT 0) AND (ABS(xmid-127) LT 10) THEN back_maskend[imaskindex]=110
        ENDIF ELSE BEGIN
          xrange = MAX(imaskindex2d[0,*], xmax_ind) - MIN(imaskindex2d[0,*], xmin_ind)
          y_xrange = imaskindex2d[1, xmax_ind] - imaskindex2d[1, xmin_ind]
          yrange = MAX(imaskindex2d[1,*], ymax_ind) - MIN(imaskindex2d[1,*], ymin_ind)
          ;less than 11.3degree, y>5*x
          IF (y_xrange LT 0) AND (ABS(y_xrange) LT 1.3*xrange) AND (ABS(ymid-127) LT 20) AND (xmid-127 GT 0) THEN back_maskend[imaskindex]=110
        ENDELSE
      ENDIF
    ENDFOR
  ENDIF
  
  ;2020/7/13, right M1 and ICA segmentation
  ;110 is for right A1 segmentation
  right_a1_index = WHERE(back_maskend EQ 110, right_a1_ncount)
  IF right_a1_ncount GT 0 THEN BEGIN
    right_a1_index2d = ARRAY_INDICES(back_maskend, right_a1_index)
  ENDIF ELSE BEGIN
    ;right a1 could not be found
    ;error handling
    right_a1_index2d = [127,127]
  ENDELSE
  ;calculate the lower-right corner coordinates
  x0 = MAX(right_a1_index2d[0, *])
  y0 = MIN(right_a1_index2d[1, *])
  ;start sorting M1
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF (xmid GT x0) AND (ymid GT y0) AND (ABS(ymid-y0) LT 1.2*ABS(xmid-x0)) THEN back_maskend[imaskindex]=120
      ENDIF
    ENDFOR
  ENDIF
  ;start sorting ICA
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF (xmid GT x0) AND (ymid LT y0) THEN back_maskend[imaskindex]=130
      ENDIF
    ENDFOR
  ENDIF
  
  
  ;2020/7/12; Left A1 segmentation
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF incount LT 4 THEN BEGIN
          IF (ABS(ymid-127) LT 20) AND (xmid-127 LT 0) AND (ABS(xmid-127) LT 10) THEN back_maskend[imaskindex]=120
        ENDIF ELSE BEGIN
          xrange = MAX(imaskindex2d[0,*], xmax_ind) - MIN(imaskindex2d[0,*], xmin_ind)
          y_xrange = imaskindex2d[1, xmax_ind] - imaskindex2d[1, xmin_ind]
          yrange = MAX(imaskindex2d[1,*], ymax_ind) - MIN(imaskindex2d[1,*], ymin_ind)
          ;less than 11.3degree, y>5*x
          IF (y_xrange GT 0) AND (ABS(y_xrange) LT 1.3*xrange) AND (ABS(ymid-127) LT 20) AND (xmid-127 LT 0) THEN back_maskend[imaskindex]=140
        ENDELSE
      ENDIF
    ENDFOR
  ENDIF
  
  ;2020/7/13, left M1 and ICA segmentation
  ;140 is for left A1 segmentation
  left_a1_index = WHERE(back_maskend EQ 140, left_a1_ncount)
  IF left_a1_ncount GT 0 THEN BEGIN
    left_a1_index2d = ARRAY_INDICES(back_maskend, left_a1_index)
  ENDIF ELSE BEGIN
    ;left a1 could not be found
    ;error handling
    left_a1_index2d = [127,127]
  ENDELSE
  ;calculate the lower-left corner coordinates
  x0 = MIN(left_a1_index2d[0, *])
  y0 = MIN(left_a1_index2d[1, *])
  ;start sorting M1
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF (xmid LT x0) AND (ymid GT y0) AND (ABS(ymid-y0) LT 1.2*ABS(xmid-x0)) THEN back_maskend[imaskindex]=150
      ENDIF
    ENDFOR
  ENDIF
  ;start sorting ICA
  maskindex = WHERE(back_maskend GT 0, ncount)
  IF ncount GT 0 THEN BEGIN
    pdf = histogram(back_maskend, MIN=1, LOCATIONS=xbin)
    FOR i=1, MAX(xbin) DO BEGIN
      imaskindex = WHERE(back_maskend EQ i, incount)
      IF incount GT 0 THEN BEGIN
        imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
        xmid = MEAN(imaskindex2d[0,*])
        ymid = MEAN(imaskindex2d[1,*])
        IF (xmid LT x0) AND (ymid LT y0) THEN back_maskend[imaskindex]=160
      ENDIF
    ENDFOR
  ENDIF
  
  ;GJ, 2020/7/14 add the missing pixels
  vessel_num = 7
  vessel_center = DBLARR(2,vessel_num)*0.
  distance = DBLARR(vessel_num)*0.
  FOR i=0, vessel_num-1 DO BEGIN
    imaskindex = WHERE(back_maskend EQ 100+i*10, ncount)
    IF ncount GT 0 THEN BEGIN
      imaskindex2d = ARRAY_INDICES(back_maskend, imaskindex)
      xmid = MEAN(imaskindex2d[0,*])
      ymid = MEAN(imaskindex2d[1,*])
      vessel_center[0, i] = xmid
      vessel_center[1, i] = xmid
    ENDIF
  ENDFOR
  ;add together
  temp_mask = mask + back_maskend
  ;find the missing pixels
  missing_maskindex = WHERE(temp_mask EQ 1, missing_count)
  IF missing_count GT 0 THEN BEGIN
    missing_maskindex2d = ARRAY_INDICES(back_maskend, missing_maskindex)
    FOR j=0, missing_count-1 DO BEGIN
      ;calculate distance to vessel center
      FOR k=0, vessel_num-1 DO BEGIN
        distance[k] = SQRT((missing_maskindex2d[0,j]-vessel_center[0, k])^2 + (missing_maskindex2d[1,j]-vessel_center[1, k])^2)
      ENDFOR
      ;pick the closest vessel
      min_dis = MIN(distance, min_ind)
      ;assign the value
      back_maskend[missing_maskindex[j]] = 100 + min_ind*10
    ENDFOR
  ENDIF
  
;  iimage, back_maskend
;  print, 'id   sorted_id'
;  FOR i=0, N_ele-1 DO BEGIN
;    print, vline_info[i].ID, vline_info[i].sort_ID
;  ENDFOR
  
;  ;calculate the distance
;  dist_matrix = DBLARR(N_ele)
;  FOR i=0, N_ele-1 DO BEGIN
;    FOR j=0, N_ele-1 DO BEGIN
;      diff = vline_info[i].midpoint_loc - vline_info[j].midpoint_loc
;      dist_matrix[j] = SQRT(diff[0]^2 + diff[1]^2)
;    ENDFOR
;    sort_dist = SORT(dist_matrix)
;    j_ind = sort_dist[1]
;;    print, 'vline_info[i].sort_ID NE vline_info[j_ind].sort_ID', vline_info[i].sort_ID NE vline_info[j_ind].sort_ID
;    IF vline_info[i].sort_ID NE vline_info[j_ind].sort_ID THEN BEGIN
;      max_length = MAX([vline_info[i].length, vline_info[j_ind].length], max_ind)
;      smallID = ([vline_info[i].sort_ID, vline_info[j_ind].sort_ID])[max_ind]
;      mask_index = WHERE((back_maskend EQ vline_info[i].sort_ID) or (back_maskend EQ vline_info[j_ind].sort_ID), ncount)
;;      print, 'dist_matrix = ', dist_matrix
;;      print, 'sort_dist = ', sort_dist
;;      print, 'j_ind = ', j_ind
;      print, 'ncount 1= ', ncount
;;      print, 'small ID 1= ', smallID
;      IF ncount GT 0 THEN BEGIN
;        back_maskend[mask_index] = smallID
;        vline_info[i].sort_ID = smallID
;        vline_info[j_ind].sort_ID = smallID
;      ENDIF
;    ENDIF
;;    j_ind = sort_dist[2]
;;    ;    print, 'vline_info[i].sort_ID NE vline_info[j_ind].sort_ID', vline_info[i].sort_ID NE vline_info[j_ind].sort_ID
;;    IF vline_info[i].sort_ID NE vline_info[j_ind].sort_ID THEN BEGIN
;;      max_length = MAX([vline_info[i].length, vline_info[j_ind].length], max_ind)
;;      smallID = ([vline_info[i].sort_ID, vline_info[j_ind].sort_ID])[max_ind]
;;      mask_index = WHERE((back_maskend EQ vline_info[i].sort_ID) or (back_maskend EQ vline_info[j_ind].sort_ID), ncount)
;;;      print, 'dist_matrix = ', dist_matrix
;;;      print, 'sort_dist = ', sort_dist
;;;      print, 'j_ind = ', j_ind
;;      print, 'ncount 2= ', ncount
;;;      print, 'small ID = ', smallID
;;      IF ncount GT 0 THEN BEGIN
;;        back_maskend[mask_index] = smallID
;;        vline_info[i].sort_ID = smallID
;;        vline_info[j_ind].sort_ID = smallID
;;      ENDIF
;;    ENDIF
;  ENDFOR
  
;  iimage, back_maskend_s
;  iimage, BYTSCL(back_maskend)
;  iimage, maskend



;  ; Create another window and display the scaled transform.
;  WINDOW, 3, XSIZE = displaySize[0] + 1.5*offset[0], $
;    YSIZE = displaySize[1] + 1.5*offset[1], $
;    TITLE = 'Scaled Hough Transform'
;  TVSCL, REVERSE(CONGRID(maskend, displaySize[0], $
;    displaySize[1]),2), offset[0], offset[1]
;  PLOT, theta*180./!PI, rho, /XSTYLE, /YSTYLE, $
;    TITLE = 'Scaled Hough Transform', XTITLE = 'Theta', $
;    YTITLE = 'Rho', /NODATA, /NOERASE, /DEVICE, $
;    POSITION = [offset[0], offset[1], $
;    displaySize[0] + offset[0], $
;    displaySize[1] + offset[1]], CHARSIZE = 1.5, $
;    COLOR = !P.BACKGROUND

;  ; Backproject to compare with original image.
;  ;backprojection = HOUGH(transform, /BACKPROJECT, $
;  backprojection = HOUGH(maskend, /BACKPROJECT, $
;    RHO = rho, THETA = theta, $
;    NX = intensitySize[0], NY = intensitySize[1])
;
;  ;; Create another window and display the results.
;  ;WINDOW, 4, XSIZE = intensitySize[0], $
;  ;   YSIZE = intensitySize[1], $
;  ;   TITLE = 'Resulting Power Lines'
;  ;TVSCL, backprojection
;
;  ;iimage, backprojection

END


pro test_struc_vlineobject

  IF N_ELEMENTS(struc_vlineobject) EQ 0 THEN BEGIN
    ;; create the struc_vlineobject structure
    cur_vline=0
    ;; create a struc_vlineobject structure
    struc_vlineobject = {struc_vlineobject}
    struc_vlineobject.hide_Status =+1
  ENDIF ELSE BEGIN
    p=replicate({struc_vlineobject},n+1)
    p[0:n-1]=struc_vlineobject[0:n-1]
    cur_vline=n
    struc_vlineobject=p
  ENDELSE


end