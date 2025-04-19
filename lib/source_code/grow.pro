
;  img=vol_HU_iso[*,*,20]
;  
;  imgDims = SIZE(img, /DIMENSIONS)
;  x0 = 234
;
;  y0 = 169
;
;  x = INDGEN(16*16) MOD 16 + x0
;
;  y = INDGEN(16*16) / 16 + y0
;
;  roiPixels = x + y * imgDims[0]


pro growmax,a,thmin,thmax,b

vol_HU_iso=a


threshArray = (vol_HU_iso ge thmin and vol_HU_iso le thmax)

; Establish error handler. When errors occur, the index of the
; error is returned in the variable Error_status:
CATCH, Error_status
;This statement begins the error handler:
IF Error_status NE 0 THEN BEGIN
  ulong = 1
  labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
  CATCH, /CANCEL
ENDIF
  labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
  labelmax=max(labelImg)
  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
  h[0]=0
  
 print,'max', max(h)
 hmax=where(h eq max(h))
 if hmax ge 0 then begin
  maskend = labelImg EQ hmax[0]
  pos=where(maskend  eq 1)

  b=maskend
  endif else begin
    
    b=-1
  endelse
  end
  
  
pro growNOTmax15,a,thmin,thmax,maskend,n_grows

  rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]
  ;  vol_HU_iso=a
  maskend=a*0
  n_grows = 0

  threshArray = (a ge thmin and a le thmax)


  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    ulong = 1
    labelImg = LABEL_REGION(threshArray, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    CATCH, /CANCEL
  ENDIF
  labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
  labelmax=max(labelImg)
  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)


  sort_re=REVERSE(SORT(h))

  print,'max', max(h)
  hmax=n_elements(sort)
  if labelmax ge 15 then begin
    n_grows=15
    ;mymodel = OBJ_NEW('IDLgrModel')
    for j=1,15 do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
      ;      SHADE_VOLUME, maskend,0.9, v, p
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      rgbi=j mod 7
      ;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
      ;
      ;      mymodel->Add, oPoly

    endfor
    ;    xobjview_XDU, mymodel, BACKGROUND = [0, 0, 0]

    ;    maskend = labelImg EQ hmax[0]
    ;    pos=where(maskend  eq 1)

    ;    b=maskend
  endif

  if   labelmax lt 15 then begin

    if   labelmax le 0 then  begin

      ;b=-1
      n_grows = 0
    endif
    ;mymodel = OBJ_NEW('IDLgrModel')
    n_grows = labelmax
    for j=1,labelmax do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
      ;      SHADE_VOLUME, maskend,0.9, v, p
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      rgbi=j mod 7
      ;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
      ;
      ;      mymodel->Add, oPoly

    endfor
    ;    xobjview_XDU, mymodexobjview_XDUROUND = [0, 0, 0]

    ;    b=maskend

  endif
  
  maskend = threshArray - maskend
  print, 'n_grows = ', n_grows
    
end

  
  
 pro growpixel,a,x,y,z,thmin,b

  vol_HU_iso=a
  HUsize = SIZE(vol_HU_iso)
  roiPixels = FLOOR(x) + FLOOR(y) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  threshArray = (vol_HU_iso ge thmin); and vol_HU_iso le thmax)

  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    ulong = 1
    labelImg = LABEL_REGION(threshArray,  ALL_NEIGHBORS=allNeighbors, ULONG=ulong)
    CATCH, /CANCEL
  ENDIF
  labelImg = LABEL_REGION(threshArray,  ALL_NEIGHBORS=allNeighbors, ULONG=ulong)
  
  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
  
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
  
 

hsize=size(h)
print,'hszie',hsize[1]
flag=0
  for i=1,hsize[1]-1 do begin
;    print,'i',i
   maskend = labelImg EQ i
    
     pos=where(maskend  eq 1)
    
    posresult=where(pos eq roiPixels )
     if (posresult ge 0) then begin
      
      b=maskend
      flag=1
      break
;      SHADE_VOLUME, b,0.9, v, p
;
;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;      oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedVertices , POLYGON=p)
;      mymodel1 = OBJ_NEW('IDLgrModel')
;      mymodel1->Add, oPoly1
;      XOBJVIEW_XDU, mymodel1, BACKGROUND = [0, 0, 0]

     endif
     
  endfor
if  (flag eq 0) then begin
  print,'not fond'
  b=-1
endif

;  print,'max', max(h)
;  hmax=where(h eq max(h))
;  maskend = labelImg EQ hmax[0]
;  pos=where(maskend  eq 1)
;
;  b=maskend

end 
  
  
  
 pro growpixel8,a,x,y,z,thmin,b,notb
 print,'x',FLOOR(x)
  print,'y',FLOOR(y)
   print,'z',FLOOR(z)
roiPixels=lonarr(27)
  vol_HU_iso=a
  HUsize = SIZE(vol_HU_iso)
  roiPixels[0] = FLOOR(x-1) + FLOOR(y+1) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[1] = FLOOR(x-1) + FLOOR(y) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[2] = FLOOR(x-1) + FLOOR(y-1) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[3] = FLOOR(x) + FLOOR(y+1) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[4] = FLOOR(x) + FLOOR(y) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[5] = FLOOR(x) + FLOOR(y-1) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[6] = FLOOR(x+1) + FLOOR(y+1) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[7] = FLOOR(x+1) + FLOOR(y) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  roiPixels[8] = FLOOR(x+1) + FLOOR(y-1) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  
  roiPixels[9] = FLOOR(x-1) + FLOOR(y+1) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[10] = FLOOR(x-1) + FLOOR(y) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[11] = FLOOR(x-1) + FLOOR(y-1) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[12] = FLOOR(x) + FLOOR(y+1) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[13] = FLOOR(x) + FLOOR(y) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[14] = FLOOR(x) + FLOOR(y-1) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[15] = FLOOR(x+1) + FLOOR(y+1) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[16] = FLOOR(x+1) + FLOOR(y) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  roiPixels[17] = FLOOR(x+1) + FLOOR(y-1) * HUsize[2]+FLOOR(z-1)*HUsize[1]*HUsize[2]
  
  roiPixels[18] = FLOOR(x-1) + FLOOR(y+1) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[19] = FLOOR(x-1) + FLOOR(y) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[20] = FLOOR(x-1) + FLOOR(y-1) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[21] = FLOOR(x) + FLOOR(y+1) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[22] = FLOOR(x) + FLOOR(y) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[23] = FLOOR(x) + FLOOR(y-1) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[24] = FLOOR(x+1) + FLOOR(y+1) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[25] = FLOOR(x+1) + FLOOR(y) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  roiPixels[26] = FLOOR(x+1) + FLOOR(y-1) * HUsize[2]+FLOOR(z+1)*HUsize[1]*HUsize[2]
  
  
 ; roiPixels = FLOOR(x) + FLOOR(y) * HUsize[2]+FLOOR(z)*HUsize[1]*HUsize[2]
  threshArray = (vol_HU_iso ge thmin); and vol_HU_iso le thmax)

  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    ULONG = 1
    labelImg = LABEL_REGION(threshArray,  ALL_NEIGHBORS=allNeighbors, /ULONG)
    CATCH, /CANCEL
  ENDIF
  labelImg = LABEL_REGION(threshArray,  ALL_NEIGHBORS=allNeighbors, ULONG=ulong)

  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)

  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)



  hsize=size(h)
  print,'hszie',hsize[1]
  flag=0
  for i=1,hsize[1]-1 do begin
;    print,'i',i
    maskend = labelImg EQ i

    pos=where(maskend  eq 1)

    ;posresult=where(pos eq roiPixels )
    posresult= SETINTERSECTION(pos , roiPixels)
    if (posresult[0] ge 0) then begin

      b=maskend
      notb = threshArray AND (1-b)
      flag=1
      break
      ;      SHADE_VOLUME, b,0.9, v, p
      ;
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedVertices , POLYGON=p)
      ;      mymodel1 = OBJ_NEW('IDLgrModel')
      ;      mymodel1->Add, oPoly1
      ;      xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]

    endif

  endfor
  if  (flag eq 0) then begin
    print,'not fond'
    b=-1
    notb = threshArray
  endif

  ;  print,'max', max(h)
  ;  hmax=where(h eq max(h))
  ;  maskend = labelImg EQ hmax[0]
  ;  pos=where(maskend  eq 1)
  ;
  ;  b=maskend

end 
  
;  
;  SHADE_VOLUME, maskend,0.9, v, p
;;FOR i=1, labelmax do  begin
;smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedVertices , POLYGON=p)
;mymodel1 = OBJ_NEW('IDLgrModel')
;mymodel1->Add, oPoly1
;xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]
;
;;     flag=0
;;     if(h[i] le 5) then continue
;
;     maskImg = labelImg EQ 1
;     pos=where(maskImg  eq 1)
;;endfor
;  ;img=vol_HU_iso
;
;  imgDims = SIZE(img, /DIMENSIONS)
;
;
;
;  ; Define original region pixels.
;
; 
;
;
;  ; Grow the region.
;
;  newROIPixels = REGION_GROW(img, roiPixels, $
;
;    THRESHOLD=[400,2000])
;
;
;
;  ; Draw the original image
;
;  im1 = IMAGE(img, LAYOUT=[1,2,1], DIM=[700, 900])
;
;  p = POLYGON(x0+[0,15,15,0], y0+[0,0,15,15], $
;
;    COLOR='Lime', FILL_COLOR='Lime', /DATA)
;
;
;
;  ; Color the grown region, and draw the new image.
;
;  img1 = img
;
;  img1[newROIPixels] = 255b
;
;  imgTrue = REBIN(img, imgDims[0], imgDims[1], 3)
;
;  imgTrue[*,*,2] = img1
;
;  im2 = IMAGE(imgTrue, LAYOUT=[1,2,2], /CURRENT)


pro grow

  afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  print,'pro11'
  IF STRLEN(afn) NE 0 THEN BEGIn
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    IF nrfile NE 0 THEN BEGIN
      obj = OBJ_NEW('IDLffDICOM', a[0]);, /NO_PIXEL_DATA)
      ; Get the row & column size of the image(s):
      rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
      cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
      OBJ_DESTROY, obj
      vol_HU = DBLARR(rows, cols, nrfile)
      voxelsize = DBLARR(nrfile, 6) * 0.
      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOM', a[i])
        result=QUERY_DICOM(a[i], info)

        count=(i+1.)/nrfile*100.0


        IF result NE 0 THEN BEGIN
          ; Get the image data
          array=obj->GetValue('7fe0'x, '0010'x)
          vPixels=*array[0]
          PTR_FREE, array

          ; Get the row & column size of the image(s):
          rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
          cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
          
          spacing = FLOAT(*(obj->GetValue('0018'x,'0088'x,/NO_COPY))[0])
          sliceTh = FLOAT(*(obj->GetValue('0018'x,'0050'x,/NO_COPY))[0])
          pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
          pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
          ptPosition = STRING(*(obj->GetValue('0018'x,'5100'x,/NO_COPY))[0])
          sliceLo = STRING(*(obj->GetValue('0020'x,'1041'x,/NO_COPY))[0])
          
          ;确认图像分辨率读取正确
          ;          print, "voxel size: ", i, ':', voxelsize[i, *]
          IF sliceTh GT 0. THEN BEGIN
            voxelsize[i, 0] = pixelSp[0]
            voxelsize[i, 1] = pixelSp[1]
            IF spacing LE 0. THEN voxelsize[i, 2] = sliceTh + spacing ELSE voxelsize[i, 2] = spacing ;sliceTh +
            voxelsize[i, 3] = rows
            voxelsize[i, 4] = cols
          ENDIF

          voxelsize[i, 5] = sliceLo
          ;确认图像分辨率读取正确
          ;  print, "voxel size: ", i, ':', voxelsize[i, *]

          ;convert grey pixel to HU
          ;rescale intercept
          rescale_intercept = FLOAT(*(obj->GetValue('0028'x,'1052'x,/NO_COPY))[0])
          rescale_slope = FLOAT(*(obj->GetValue('0028'x,'1053'x,/NO_COPY))[0])
          ;          IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
          ;          ;rescale slope
          ;          IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 1.
          OBJ_DESTROY, obj
          vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)

          ;smooth the image and store to 3D
          vol_HU[*,*,i] = vPixels_HU
        ENDIF
      ENDFOR


      ;得到正确的层厚
      sliceDistance = ABS(voxelsize[0, 5] - voxelsize[1, 5])

      ;;减少图像尺寸
      IF FLOOR(voxelsize[0, 3]*voxelsize[0, 0]/sliceDistance) LT 512 THEN BEGIN
        vol_HU_iso = CONGRID(vol_HU,  FLOOR(voxelsize[0, 3]*voxelsize[0, 0]/sliceDistance), FLOOR(voxelsize[0, 3]*voxelsize[0, 0]/sliceDistance), nrfile)
      ENDIF ELSE BEGIN
        FOV = voxelsize[0, 3]*voxelsize[0, 0]
        length = nrfile * sliceDistance
        vol_HU_iso = CONGRID(vol_HU,  512, 512, FLOOR(length/(FOV/512.)))
      ENDELSE


      vol_HU_iso = CONGRID(vol_HU,  FLOOR(voxelsize[0, 3]*voxelsize[0, 0]/voxelsize[0, 2]), FLOOR(voxelsize[0, 4]*voxelsize[0, 1]/voxelsize[0, 2]), nrfile)

    ENDIF
  ENDIF




  growmaxfive,vol_HU_iso,250,2000,b
  ;growpixel,vol_HU_iso,238,292,20,400,2000,b
  sizeb=size(b)
  if(sizeb[0] eq 0) then begin
    print,'not found'
  endif  else begin

;    SHADE_VOLUME, b,0.9, v, p
;
;    smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;    oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedVertices , POLYGON=p)
;    mymodel1 = OBJ_NEW('IDLgrModel')
;    mymodel1->Add, oPoly1
;    xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]
  endelse

end

pro growmax15,a,thmin,thmax,maskend,n_grows
  rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]
  ;  vol_HU_iso=a
  maskend=a*0
  n_grows = 0

  threshArray = (a ge thmin and a le thmax)


  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    ulong = 1
    labelImg = LABEL_REGION(threshArray, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    CATCH, /CANCEL
  ENDIF
  labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
  labelmax=max(labelImg)
  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)


  sort_re=REVERSE(SORT(h))

  print,'max', max(h)
  hmax=n_elements(sort)
  if labelmax ge 15 then begin
    n_grows=15
    ;mymodel = OBJ_NEW('IDLgrModel')
    for j=1,15 do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
      ;      SHADE_VOLUME, maskend,0.9, v, p
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      rgbi=j mod 7
      ;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
      ;
      ;      mymodel->Add, oPoly

    endfor
    ;    xobjview_XDU, mymodel, BACKGROUND = [0, 0, 0]

    ;    maskend = labelImg EQ hmax[0]
    ;    pos=where(maskend  eq 1)

    ;    b=maskend
  endif

  if   labelmax lt 15 then begin

    if   labelmax le 0 then  begin

      ;b=-1
      n_grows = 0
    endif
    ;mymodel = OBJ_NEW('IDLgrModel')
    n_grows = labelmax
    for j=1,labelmax do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
      ;      SHADE_VOLUME, maskend,0.9, v, p
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      rgbi=j mod 7
      ;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
      ;
      ;      mymodel->Add, oPoly

    endfor
    ;    xobjview_XDU, mymodexobjview_XDUROUND = [0, 0, 0]

    ;    b=maskend

  endif

  print, 'n_grows = ', n_grows

end

pro growmaxten,a,thmin,thmax,maskend,n_grows
  rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]
  ;  vol_HU_iso=a
  maskend=a*0
  n_grows = 0

  threshArray = (a ge thmin and a le thmax)


  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    ulong = 1
    labelImg = LABEL_REGION(threshArray, $
      ALL_NEIGHBORS=allNeighbors, $
      ULONG=ulong)
    CATCH, /CANCEL
  ENDIF
  labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
  labelmax=max(labelImg)
  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)


  sort_re=REVERSE(SORT(h))

  print,'max', max(h)
  hmax=n_elements(sort)
  if labelmax ge 10 then begin
    n_grows=10
    ;mymodel = OBJ_NEW('IDLgrModel')
    for j=1,10 do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
      ;      SHADE_VOLUME, maskend,0.9, v, p
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      rgbi=j mod 7
      ;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
      ;
      ;      mymodel->Add, oPoly

    endfor
    ;    xobjview_XDU, mymodel, BACKGROUND = [0, 0, 0]

    ;    maskend = labelImg EQ hmax[0]
    ;    pos=where(maskend  eq 1)

    ;    b=maskend
  endif

  if   labelmax lt 10 then begin

    if   labelmax le 0 then  begin

      ;b=-1
      n_grows = 0
    endif
    ;mymodel = OBJ_NEW('IDLgrModel')
    n_grows = labelmax
    for j=1,labelmax do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
      ;      SHADE_VOLUME, maskend,0.9, v, p
      ;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
      ;      rgbi=j mod 7
      ;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
      ;
      ;      mymodel->Add, oPoly

    endfor
    ;    xobjview_XDU, mymodexobjview_XDUROUND = [0, 0, 0]

    ;    b=maskend

  endif

  print, 'n_grows = ', n_grows

end

pro growmaxfive,a,thmin,thmax,maskend,n_grows
  rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]
;  vol_HU_iso=a
  maskend=a*0
  n_grows = 0

  threshArray = (a ge thmin and a le thmax)


  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    ulong = 1
    labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
    CATCH, /CANCEL
  ENDIF
  labelImg = LABEL_REGION(threshArray, $
    ALL_NEIGHBORS=allNeighbors, $
    ULONG=ulong)
  labelmax=max(labelImg)
  ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
  h = HISTOGRAM(labelImg, REVERSE_INDICES=r)


  sort_re=REVERSE(SORT(h))

  print,'max', max(h)
  hmax=n_elements(sort)
  if labelmax ge 5 then begin
    n_grows=5
    ;mymodel = OBJ_NEW('IDLgrModel')
    for j=1,5 do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
;      SHADE_VOLUME, maskend,0.9, v, p
;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;      rgbi=j mod 7
;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
;
;      mymodel->Add, oPoly

    endfor
;    xobjview_XDU, mymodel, BACKGROUND = [0, 0, 0]

    ;    maskend = labelImg EQ hmax[0]
    ;    pos=where(maskend  eq 1)

;    b=maskend
  endif

  if   labelmax lt 5 then begin

    if   labelmax le 0 then  begin

      ;b=-1
      n_grows = 0
    endif
    ;mymodel = OBJ_NEW('IDLgrModel')
    n_grows = labelmax
    for j=1,labelmax do begin
      maskendmid = labelImg EQ sort_re[j]
      maskend+=maskendmid*j
;      SHADE_VOLUME, maskend,0.9, v, p
;      smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;      rgbi=j mod 7
;      oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
;
;      mymodel->Add, oPoly

    endfor
;    xobjview_XDU, mymodexobjview_XDUROUND = [0, 0, 0]

;    b=maskend

  endif

 print, 'n_grows = ', n_grows

end






function setintersection, a, b
  ;
  compile_opt StrictArr

  minab = min(a, Max=maxa) > min(b, Max=maxb) ;Only need intersection of ranges
  maxab = maxa < maxb
  ;
  ; If either set is empty, or their ranges don't intersect: result = NULL.
  if maxab lt minab or maxab lt 0 then return, -1
  r = where((histogram(a, Min=minab, Max=maxab) ne 0) and $
    (histogram(b, Min=minab, Max=maxab) ne 0), count)

  if count eq 0 then return, -1 else return, r + minab
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1222
pro change_color,verts,maxdistance_vol_hu_cube

  print,max(maxdistance_vol_hu_cube)
  iimage,maxdistance_vol_hu_cube[*,*,67]
  help,verts,maxdistance_vol_hu_cube
  temp=size(verts)
  print, 'x', min(verts[0,*]), 'to ', max(verts[0,*])
  print, 'y', min(verts[1,*]), 'to ', max(verts[1,*])
  print, 'z', min(verts[2,*]), 'to ', max(verts[2,*])
  ;    for i=0,temp[2]-1 do begin
  for i=1500,1540 do begin
    index=verts[*,i]
    ;      a=floor(index[0]*3)
    ;      b=floor(index[1]*3)
    ;      c=floor(index[2]*3)
    a=CEIL(index[0])
    b=CEIL(index[1])
    c=CEIL(index[2])
    number=maxdistance_vol_hu_cube[a,b,c]
    print,'a:',a,'b:',b,'c:',c,',  number',number

  endfor

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1222
