Pro separate

  ;;去掉小链接
 
  afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  print,'pro11'
  IF STRLEN(afn) NE 0 THEN BEGIn
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    IF nrfile NE 0 THEN BEGIN
      obj = OBJ_NEW('IDLffDICOM', a[0])
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


;SHADE_VOLUME, vol_HU_iso,160, v2, p2
INTERVAL_VOLUME, vol_HU_iso, 160, 2000, v2, p2
p2 = TETRA_SURFACE(v2, p2)
v2 = MESH_SMOOTH(v2, p2, ITERATIONS=200)
;p2 = TETRA_SURFACE(v2, p2)


oPoly2 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], v2, POLYGON=p2)
mymodel2 = OBJ_NEW('IDLgrModel')
mymodel2->Add, oPoly2
;
XOBJVIEW, mymodel2, BACKGROUND = [0, 0, 0]




print,'pro1'
findmax, vol_HU_iso,nrfile,vol_HU_iso_max


vol_HU_iso_max[0, *, *] = (vol_HU_iso_max[0, *, *] GT 160) * 160
vol_HU_iso_max[(size(vol_HU_iso_max))[1]-1, *, *] = (vol_HU_iso_max[(size(vol_HU_iso_max))[1]-1, *, *] GT 160) * 160
vol_HU_iso_max[*, 0, *] = (vol_HU_iso_max[*, 0, *] GT 160) * 160
vol_HU_iso_max[*, (size(vol_HU_iso_max))[2]-1, *] = (vol_HU_iso_max[*, (size(vol_HU_iso_max))[2]-1, *] GT 160) * 160
vol_HU_iso_max[*, *, 0] = (vol_HU_iso_max[*, *, 0] GT 160) * 160
vol_HU_iso_max[*, *, (size(vol_HU_iso_max))[3]-1] = (vol_HU_iso_max[*, *, (size(vol_HU_iso_max))[3]-1] GT 160) * 160
 

 SHADE_VOLUME, vol_HU_iso_max,160, v1, p1
 oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], v1, POLYGON=p1)
 mymodel1 = OBJ_NEW('IDLgrModel')
 mymodel1->Add, oPoly1
  XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]

  i=0

;division1, vol_HU_iso,nrfile
i=0
END

;  SHADE_VOLUME, arr,0.9, v1, p1


;  oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], v1, POLYGON=p1)
;  mymodel1 = OBJ_NEW('IDLgrModel')
;  mymodel1->Add, oPoly1
;
;  XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]

Pro findmax,a,b,c, threshold
  vol_HU_iso=a
vol_HU_iso_max=a
  nrfile=b
binaryVol = vol_HU_iso GE threshold
sz=size(binaryVol[*,*,0])
xsize=sz[1]
ysize=sz[2]


wholeone=bytarr(xsize,ysize)+1
 for nf=0,nrfile-1 do begin

   labelImg = LABEL_REGION(binaryVol[*,*,nf],/ALL_NEIGHBORS)
   labelmax=max(labelImg)
   ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
   h = HISTOGRAM(labelImg, REVERSE_INDICES=r)

  print,h
  
  reversesot=REVERSE(sort(h))
  
  hmaxi1=reversesot[1]
  if labelmax gt 1 then begin
  hmaxi2=reversesot[2]
  
  if h[ hmaxi1] gt (h[hmaxi2]*2) then begin
 maskImg = labelImg EQ hmaxi1
  pos=where(maskImg  eq 1)
  endif else begin
maskImg1 = labelImg EQ hmaxi1
pos1=where(maskImg1  eq 1)
maskImg2 = labelImg EQ hmaxi2
pos2=where(maskImg2  eq 1)
maskImg = maskImg1 +maskImg2 
  pos=pos1+pos2
  endelse
  
  endif else begin
    
    
    maskImg = labelImg EQ hmaxi1
  endelse
;  hmax=max(h[1:labelmax])
;  print,hmax
; i=where(h eq  hmax)
;  print,h[i]
;
;    maskImg = labelImg EQ i[0]
;    pos=where(maskImg  eq 1)


 ;maskImg = labelImg EQ hmaxi1

vol_HU_iso_max[*,*,nf]=maskImg and vol_HU_iso[*,*,nf]

revermaskImg=(wholeone  xor maskImg)*(-1000)
vol_HU_iso_max[*,*,nf]=vol_HU_iso_max[*,*,nf]+revermaskImg
temp= pos
endfor
h=0

c=vol_HU_iso_max


end


Pro findmax2,a,b,c

  vol_HU_iso=a
  vol_HU_iso_max=a
  nrfile=b
  binaryVol = vol_HU_iso GE -200.
  sz=size(binaryVol[*,*,0])
  xsize=sz[1]
  ysize=sz[2]



  for nf=0,nrfile-1 do begin

labelImg = LABEL_REGION(binaryVol[*,*,nf],/ALL_NEIGHBORS)
 labelmax=max(labelImg)
 
h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)

 FOR i=1, labelmax do  begin

   flag=0
   if(h[i] le 5) then continue

   maskImg = labelImg EQ i
   pos=where(maskImg  eq 1)

for j=nf+1, nrfile-1 do  begin

 if(flag eq 0)then begin
  comp10=pos
  comp11=pos+1
  comp12=pos-1
  comp13=pos+xsize
  comp14=pos-xsize
endif else begin
  pos=poslist.toarray()
  comp10=pos
  comp11=pos+1
  comp12=pos-1
  comp13=pos+xsize
  comp14=pos-xsize
  ;obj_destroy,poslist
endelse
flag=1
labelImg1 = LABEL_REGION(binaryVol[*,*,j],/ALL_NEIGHBORS)
labelmax1=max(labelImg1)
poslist=list()
FOR ii=1, labelmax1 do  begin
  maskImg1 = labelImg1 EQ ii
  h1 = HISTOGRAM(labelImg1, REVERSE_INDICES=r)
  if(h1[ii] le 5) then continue
pos1=where(maskImg1  eq 1)

RESULT0=SETINTERSECTION(POS1,COMP10)
result1=setintersection(pos1,comp11)
result2=setintersection(pos1,comp12)
result3=setintersection(pos1,comp13)
result4=setintersection(pos1,comp14)
result=(result0[0] ge 0)or(result1[0] ge 0)or(result2[0] ge 0)or(result3[0] ge 0)or(result4[0] ge 0)
if(result0[0] ge 0)or(result1[0] ge 0)or(result2[0] ge 0)or(result3[0] ge 0)or(result4[0] ge 0) then begin
  nflist.add,j
  hilist.add,ii

endif
endfor
endfor
endfor
endfor
end

PRO division1,a,b
vol_HU_iso=a
nrfile=b
 rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]

 binaryVol = vol_HU_iso GE -200.
 ;maskimgwhole=arr*0
 sz=size(binaryVol[*,*,0])
 xsize=sz[1]
 ysize=sz[2]


 SHADE_VOLUME, binaryVol,0.9, v1, p1

 ;
 oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], v1, POLYGON=p1)
 mymodel1 = OBJ_NEW('IDLgrModel')
 mymodel1->Add, oPoly1
 ;
 XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]







 maxh=0
 maxnf=0



 nflist=list()
 hilist=list()
 mymodel = OBJ_NEW('IDLgrModel')
 ;nrfile=60
 ;nf=60
 for nf=0,nrfile-1 do begin


   print,'nf',nf
   print,'circl1'
   ;for nf=0,nrfile-1 do begin
   ;nf=0
   labelImg = LABEL_REGION(binaryVol[*,*,nf],/ALL_NEIGHBORS)
   labelmax=max(labelImg)
   ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
   h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
   ; print,'h',h
   nflist.add,-1
   hilist.add,-1
   FOR i=1, labelmax do  begin

     print,'i',i
     ; print,'circl2'

     flag=0
     if(h[i] le 5) then continue
     rem=0
     remc=0
     d=0
     ; poslist=list()
     ; if(h[i] le 40)  then continue
     for p=0,N_ELEMENTS(nflist)-1  do  begin
       remc=(nflist[p] eq nf  ) and(hilist[p] eq i )and(nf ne 0)

       if(remc eq 1) then begin
         rem=1
       endif
       if (rem eq 1)then break
     endfor


     if (rem eq 1)then continue



     nflist.add,nf
     hilist.add,i
     ;      mainRegion = WHERE(HISTOGRAM(labelImg) EQ  h[i])
     maskImg = labelImg EQ i
     pos=where(maskImg  eq 1)

     for j=nf+1, nrfile-1 do  begin
       print,'j',j
       ; print,'circl3'
       c=0
       ;d=0
       cc=0
       ; print,'dd'
       ;print,'j',j
       if(flag eq 0)then begin
         comp10=pos
         comp11=pos+1
         comp12=pos-1
         comp13=pos+xsize
         comp14=pos-xsize
       endif else begin
         pos=poslist.toarray()
         comp10=pos
         comp11=pos+1
         comp12=pos-1
         comp13=pos+xsize
         comp14=pos-xsize
         ;obj_destroy,poslist
       endelse
       flag=1
       labelImg1 = LABEL_REGION(binaryVol[*,*,j],/ALL_NEIGHBORS)
       labelmax1=max(labelImg1)
       poslist=list()
       FOR ii=1, labelmax1 do  begin

         maskImg1 = labelImg1 EQ ii
         h1 = HISTOGRAM(labelImg1, REVERSE_INDICES=r)
         if(h1[ii] le 5) then continue

         pos1=where(maskImg1  eq 1)
         ; print,'pos',pos
         ; print,'pos1',pos1
         RESULT0=SETINTERSECTION(POS1,COMP10)
         result1=setintersection(pos1,comp11)
         result2=setintersection(pos1,comp12)
         result3=setintersection(pos1,comp13)
         result4=setintersection(pos1,comp14)
         result=(result0[0] ge 0)or(result1[0] ge 0)or(result2[0] ge 0)or(result3[0] ge 0)or(result4[0] ge 0)
         ;print,'result',result
         if(result0[0] ge 0)or(result1[0] ge 0)or(result2[0] ge 0)or(result3[0] ge 0)or(result4[0] ge 0) then begin
           nflist.add,j
           hilist.add,ii
           c=1
           d=1
           for posi=0,N_elements(pos1)-1 do begin

             poslist.add,pos1[posi]

           endfor
         endif
       endfor
       if(c eq 0) then break
     endfor
     nflist.add,-1
     hilist.add,-1
   endfor
   nflist.add,-2
   hilist.add,-2
 endfor
 print,'nflist',nflist
 print,'hilist',hilist



 nflistmid=list()
 hilistmid=list()
 remeber=list()


 nfarr=nflist.ToArray()
 hiarr=hilist.ToArray()
 listresult=where(nfarr eq -1,countlist)

 for lisi=0,countlist-2 do begin
   listmidnf=list()
   listmidhi=list()
   nflistmid.add,-1
   hilistmid.add,-1

   lisia=0
   if(n_elements(remeber) gt 0)  then begin
     remeberarr=remeber.ToArray()
     RESULTlisi=where(remeberarr eq lisi)
     if(RESULTlisi[0] ge 0)then begin
       lisia=1
     endif
   endif
   print,'lisia',lisia
   if (lisia eq 1)then continue


   compnf1=nfarr(listresult[lisi]+1:listresult[lisi+1]-1)
   comphi1=hiarr(listresult[lisi]+1:listresult[lisi+1]-1)
   ;nflistmid.add,nflist(listresult[lisi]+1:listresult[lisi+1]-1)
   ;hilistmid.add,hilist(listresult[lisi]+1:listresult[lisi+1]-1)
   for ia=0,n_elements(compnf1)-1 do begin
     listmidnf.add,compnf1[ia]
     listmidhi.add,comphi1[ia]
   endfor
   ;  print,'compnf1',compnf1
   ;  print,'comphi1',comphi1
   ;print,'compnf1[0]',compnf1[0]
   if(compnf1[0] lt 0)then  continue
   for m=lisi+1 , countlist-2 do begin
     print,'m',m


     a=0
     compnf1=listmidnf.ToArray()
     comphi1=listmidhi.ToArray()
     ; print,'compnf1',compnf1
     ;  print,'comphi1',comphi1
     compnf2=nfarr(listresult[m]+1:listresult[m+1]-1)
     comphi2=hiarr(listresult[m]+1:listresult[m+1]-1)
     if(compnf2[0] lt 0)then  continue
     ; print,'compnf2',compnf2
     ; print,'comphi2',comphi2


     for i=0,n_elements(compnf1)-1 do begin

       p1=compnf1[i]
       q1=comphi1[i]

       ; print,'p1',p1
       ; print,'q1',q1
       for j=0,n_elements(compnf2)-1 do begin

         p2=compnf2[j]
         q2=comphi2[j]
         ; print,'p2',p2
         ; print,'q2',q2
         if( p1 eq p2 )  and(q1 eq  q2  )then  begin

           a=1
           print,a
         endif

         if(a eq 1)then break
       endfor
       if(a eq 1)then break
     endfor

     ; print,'a',a









     ;  resultnf=setintersection(compnf1,compnf2)
     ;  resulthi=setintersection(comphi1,comphi2)
     ;  a=((resultnf[0] ge 0 )and (resulthi[0] ge 0))? 1:0

     ;print,'a',a

     if(a eq 1)then begin

       for ia=long(0),n_elements(compnf2)-1 do begin
         ;print,ia
         listmidnf.add,compnf2[ia]
         listmidhi.add,comphi2[ia]
       endfor
       remeber.add,m
       ;  print,'listmidnf',listmidnf
       ; print,'listmidhi',listmidhi
     endif
     ; print,'remeber',remeber
   endfor
   nflistmid=nflistmid+listmidnf
   hilistmid=hilistmid+listmidhi
   ;print,'next'
 endfor
 ;print,'nflistmid',nflistmid
 ; print,'hilistmid',hilistmid

 ; listresult2=where(nflist eq -2,countlist2)
 ;
 listresult=where(nflistmid eq -1,countlist)
 ; print,countlist

 num=0
 for lisi=0,countlist-2 do begin
   print,'lisi',lisi
   if (listresult[lisi+1]-listresult[lisi] le 1)then continue
   maskimgwhole=binaryVol*0.
   cutnf=nflistmid(listresult[lisi]+1:listresult[lisi+1]-1)
   cuthi=hilistmid(listresult[lisi]+1:listresult[lisi+1]-1)
   cutnfarr=cutnf.toarray()
   cuthiarr=cuthi.toarray()
   sortarr=sort(cutnfarr)
   cutnfarr=cutnfarr[sortarr]

   cuthiarr=cuthiarr[sortarr]
   ; print,'cutnf',cutnfarr
   ;print,'cuthi',cuthiarr
   if(cutnf[0] lt 0) then continue
   for i=0,n_elements(cutnf)-1 do begin
     ;
     labelImgend = LABEL_REGION(binaryVol[*,*,cutnf[i]],/ALL_NEIGHBORS)
     ; hend= HISTOGRAM(labelImgend, MIN=1, REVERSE_INDICES=r)

     ;
     ;      mainRegionend = WHERE(HISTOGRAM(labelImgend) EQ hend[cuthi[i]])
     maskend= labelImgend EQ cuthi[i]
     ;
     maskimgwhole[*,*,cutnf[i]]= maskimgwhole[*,*,cutnf[i]]+ maskend
   endfor
   SHADE_VOLUME, maskimgwhole,0.9, v, p
   ;
   ;     ; Smooth entire mesh to reduce the effect of the irregular vertex.
   smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
   rgbi=lisi mod 7
   oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
   ;    ;mymodel = OBJ_NEW('IDLgrModel')
   mymodel->Add, oPoly
   ; XOBJVIEW, mymodel, BACKGROUND = [0, 0, 0]
   num=num+1
 endfor
 XOBJVIEW, mymodel, BACKGROUND = [0, 0, 0]
 ;
 ;
 ;  ;print,'countlist',countlist
 print,'num',num
 ii=0
 ;
 ;
 
 
 
end

PRO division2,a,b,c
  vol_HU_iso=a
  nrfile=b
  rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]

  binaryVol = vol_HU_iso GE -200.
  ;maskimgwhole=arr*0
  sz=size(binaryVol[*,*,0])
  xsize=sz[1]
  ysize=sz[2]


  SHADE_VOLUME, binaryVol,0.9, v1, p1

  ;
  oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], v1, POLYGON=p1)
  mymodel1 = OBJ_NEW('IDLgrModel')
  mymodel1->Add, oPoly1
  ;
  XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]







  maxh=0
  maxnf=0



  nflist=list()
  hilist=list()
  mymodel = OBJ_NEW('IDLgrModel')
  ;nrfile=60
  ;nf=60
  for nf=0,nrfile-1 do begin


    print,'nf',nf
    print,'circl1'
    ;for nf=0,nrfile-1 do begin
    ;nf=0
    labelImg = LABEL_REGION(binaryVol[*,*,nf],/ALL_NEIGHBORS)
    labelmax=max(labelImg)
    ; h = HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
    h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
    ; print,'h',h
    nflist.add,-1
    hilist.add,-1
    FOR i=1, labelmax do  begin

      print,'i',i
      ; print,'circl2'

      flag=0
      if(h[i] le 5) then continue
      rem=0
      remc=0
      d=0
      ; poslist=list()
      ; if(h[i] le 40)  then continue
      for p=0,N_ELEMENTS(nflist)-1  do  begin
        remc=(nflist[p] eq nf  ) and(hilist[p] eq i )and(nf ne 0)

        if(remc eq 1) then begin
          rem=1
        endif
        if (rem eq 1)then break
      endfor


      if (rem eq 1)then continue



      nflist.add,nf
      hilist.add,i
      ;      mainRegion = WHERE(HISTOGRAM(labelImg) EQ  h[i])
      maskImg = labelImg EQ i
      pos=where(maskImg  eq 1)

      for j=nf+1, nrfile-1 do  begin
        print,'j',j
        ; print,'circl3'
        c=0
        ;d=0
        cc=0
        ; print,'dd'
        ;print,'j',j
        if(flag eq 0)then begin
          comp10=pos
          comp11=pos+1
          comp12=pos-1
          comp13=pos+xsize
          comp14=pos-xsize
        endif else begin
          pos=poslist.toarray()
          comp10=pos
          comp11=pos+1
          comp12=pos-1
          comp13=pos+xsize
          comp14=pos-xsize
          ;obj_destroy,poslist
        endelse
        flag=1
        labelImg1 = LABEL_REGION(binaryVol[*,*,j],/ALL_NEIGHBORS)
        labelmax1=max(labelImg1)
        poslist=list()
        FOR ii=1, labelmax1 do  begin

          maskImg1 = labelImg1 EQ ii
          h1 = HISTOGRAM(labelImg1, REVERSE_INDICES=r)
          if(h1[ii] le 5) then continue

          pos1=where(maskImg1  eq 1)
          ; print,'pos',pos
          ; print,'pos1',pos1
          RESULT0=SETINTERSECTION(POS1,COMP10)
          result1=setintersection(pos1,comp11)
          result2=setintersection(pos1,comp12)
          result3=setintersection(pos1,comp13)
          result4=setintersection(pos1,comp14)
          result=(result0[0] ge 0)or(result1[0] ge 0)or(result2[0] ge 0)or(result3[0] ge 0)or(result4[0] ge 0)
          ;print,'result',result
          if(result0[0] ge 0)or(result1[0] ge 0)or(result2[0] ge 0)or(result3[0] ge 0)or(result4[0] ge 0) then begin
            nflist.add,j
            hilist.add,ii
            c=1
            d=1
            for posi=0,N_elements(pos1)-1 do begin

              poslist.add,pos1[posi]

            endfor
          endif
        endfor
        if(c eq 0) then break
      endfor
      nflist.add,-1
      hilist.add,-1
    endfor
    nflist.add,-2
    hilist.add,-2
  endfor
  print,'nflist',nflist
  print,'hilist',hilist



  nflistmid=list()
  hilistmid=list()
  remeber=list()


  nfarr=nflist.ToArray()
  hiarr=hilist.ToArray()
  listresult=where(nfarr eq -1,countlist)

  for lisi=0,countlist-2 do begin
    listmidnf=list()
    listmidhi=list()
    nflistmid.add,-1
    hilistmid.add,-1

    lisia=0
    if(n_elements(remeber) gt 0)  then begin
      remeberarr=remeber.ToArray()
      RESULTlisi=where(remeberarr eq lisi)
      if(RESULTlisi[0] ge 0)then begin
        lisia=1
      endif
    endif
    print,'lisia',lisia
    if (lisia eq 1)then continue


    compnf1=nfarr(listresult[lisi]+1:listresult[lisi+1]-1)
    comphi1=hiarr(listresult[lisi]+1:listresult[lisi+1]-1)
    ;nflistmid.add,nflist(listresult[lisi]+1:listresult[lisi+1]-1)
    ;hilistmid.add,hilist(listresult[lisi]+1:listresult[lisi+1]-1)
    for ia=0,n_elements(compnf1)-1 do begin
      listmidnf.add,compnf1[ia]
      listmidhi.add,comphi1[ia]
    endfor
    ;  print,'compnf1',compnf1
    ;  print,'comphi1',comphi1
    ;print,'compnf1[0]',compnf1[0]
    if(compnf1[0] lt 0)then  continue
    for m=lisi+1 , countlist-2 do begin
      print,'m',m


      a=0
      compnf1=listmidnf.ToArray()
      comphi1=listmidhi.ToArray()
      ; print,'compnf1',compnf1
      ;  print,'comphi1',comphi1
      compnf2=nfarr(listresult[m]+1:listresult[m+1]-1)
      comphi2=hiarr(listresult[m]+1:listresult[m+1]-1)
      if(compnf2[0] lt 0)then  continue
      ; print,'compnf2',compnf2
      ; print,'comphi2',comphi2


      for i=0,n_elements(compnf1)-1 do begin

        p1=compnf1[i]
        q1=comphi1[i]

        ; print,'p1',p1
        ; print,'q1',q1
        for j=0,n_elements(compnf2)-1 do begin

          p2=compnf2[j]
          q2=comphi2[j]
          ; print,'p2',p2
          ; print,'q2',q2
          if( p1 eq p2 )  and(q1 eq  q2  )then  begin

            a=1
            print,a
          endif

          if(a eq 1)then break
        endfor
        if(a eq 1)then break
      endfor

      ; print,'a',a









      ;  resultnf=setintersection(compnf1,compnf2)
      ;  resulthi=setintersection(comphi1,comphi2)
      ;  a=((resultnf[0] ge 0 )and (resulthi[0] ge 0))? 1:0

      ;print,'a',a

      if(a eq 1)then begin

        for ia=long(0),n_elements(compnf2)-1 do begin
          ;print,ia
          listmidnf.add,compnf2[ia]
          listmidhi.add,comphi2[ia]
        endfor
        remeber.add,m
        ;  print,'listmidnf',listmidnf
        ; print,'listmidhi',listmidhi
      endif
      ; print,'remeber',remeber
    endfor
    nflistmid=nflistmid+listmidnf
    hilistmid=hilistmid+listmidhi
    ;print,'next'
  endfor
  ;print,'nflistmid',nflistmid
  ; print,'hilistmid',hilistmid

  ; listresult2=where(nflist eq -2,countlist2)
  ;
  listresult=where(nflistmid eq -1,countlist)
  ; print,countlist

  num=0
  hmax=0
  imax=0
  for lisi=0,countlist-2 do begin
  
    if (listresult[lisi+1]-listresult[lisi] le 1)then continue
    maskimgwhole=binaryVol*0.
    if (listresult[lisi+1]-listresult[lisi] ge hmax) then begin
      
     hmax=listresult[lisi+1]-listresult[lisi]
     imax=lisi
    endif
    endfor
    lisi=imax
    cutnf=nflistmid(listresult[lisi]+1:listresult[lisi+1]-1)
    cuthi=hilistmid(listresult[lisi]+1:listresult[lisi+1]-1)
    cutnfarr=cutnf.toarray()
    cuthiarr=cuthi.toarray()
    sortarr=sort(cutnfarr)
    cutnfarr=cutnfarr[sortarr]

    cuthiarr=cuthiarr[sortarr]
    ; print,'cutnf',cutnfarr
    ;print,'cuthi',cuthiarr
   ; if(cutnf[0] lt 0) then continue
   
   
   
   nfmax=max(cutnfarr,min=nfmin)
    finalmask=binaryVol*0.
 ; order=ptr_new()
   for i=nfmin,nfmax do begin
  
   
     labelImg = LABEL_REGION(binaryVol[*,*,i],/ALL_NEIGHBORS)
       

      hend= HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
   

    nfresult=where(cutnfarr eq i )

ll=0
    
    for j=0,n_elements(nfresult)-1 do begin


      k=nfresult[j]
  print,'cuthiarr[k]',cuthiarr[k]   
      maskImg = labelImg EQ cuthiarr[k]
      

      finalmask[*,*,i]= finalmask[*,*,i]+maskImg 
  
    endfor
      
      maskimgwhole[*,*,i]=finalmask[*,*,i] and  vol_HU_iso[*,*,i]
  

    ll=0
    endfor
   
   
   
   
   
   
;    for i=0,n_elements(cutnf)-1 do begin
;      ;
;       labelImg = LABEL_REGION(binaryVol[*,*,cutnf[i]],/ALL_NEIGHBORS)
;      ;labelImgend = LABEL_REGION(binaryVol[*,*,cutnf[i]],/ALL_NEIGHBORS)
;       hend= HISTOGRAM(labelImg, MIN=1, REVERSE_INDICES=r)
;       reversesot=REVERSE(sort( hend))
;     hmaxhi=reversesot[1]
;      ;
;      ;      mainRegionend = WHERE(HISTOGRAM(labelImgend) EQ hend[cuthi[i]])
;     ; maskend= labelImgend EQ cuthi[i]
;      maskImg = labelImg EQ cuthi[i]
;      ;
;      
;      maskimgwhole[*,*,cutnf[i]]=maskImg and vol_HU_iso[*,*,cutnf[i]]
;      ;maskimgwhole[*,*,cutnf[i]]= maskimgwhole[*,*,cutnf[i]]+ maskend
;    endfor
;    SHADE_VOLUME, maskimgwhole,0.9, v, p
;    ;
;    ;     ; Smooth entire mesh to reduce the effect of the irregular vertex.
;    smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;    rgbi=lisi mod 7
;    oPoly = OBJ_NEW('IDLgrPolygon', COLOR = rgb[rgbi:rgbi+2], smoothedVertices , POLYGON=p)
;    ;    ;mymodel = OBJ_NEW('IDLgrModel')
;    mymodel->Add, oPoly
;    ; XOBJVIEW, mymodel, BACKGROUND = [0, 0, 0]
;    num=num+1
; ; endfor
;  XOBJVIEW, mymodel, BACKGROUND = [0, 0, 0]
  ;
  ;
  ;  ;print,'countlist',countlist

  ii=0
  ;
  ;

c=maskimgwhole

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

