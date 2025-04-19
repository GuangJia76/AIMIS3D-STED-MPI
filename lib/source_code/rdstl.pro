pro rdstl
  file='D:\astro2.stl'
  readstl,file,smoothedOutverts,Outconn1
  
    
 
  
 IF N_ELEMENTS(smoothedOutverts) GT 10 THEN BEGIN
                  v2=smoothedOutverts
                  ma=max(v2)
                  mi=min(v2)
                  bias=(ma-mi)/2
                  rl=mi+bias
                  rr=ma+bias
                  cc=[(-rl)/(rr-rl),1/(rr-rl)]
                  xma = max((v2)[0,*])
                  xmi = min((v2)[0,*])
                  xmid = 0.5*(xma+xmi)*cc[1]
                  yma = max((v2)[1,*])
                  ymi = min((v2)[1,*])
                  ymid = 0.5*(yma+ymi)*cc[1]
                  zma = max((v2)[2,*])
                  zmi = min((v2)[2,*])
                  zmid = 0.5*(zma+zmi)*cc[1]

          color = [[254,206,180],[216,50,46],[230,197,183],[215,220,221],[255,255,255]]
          result = round(randomu(s, 1) * 5) 
          oPoly = obj_new('IDLgrPolygon',smoothedOutverts,poly=Outconn1,color=color[*, result],$
            shading=0,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
            
            
    mymodel1 = OBJ_NEW('IDLgrModel')
    mymodel1->Add, oPoly
    ;  ;
    XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]
    
    
    
    

  ENDIF
  
;  arr=bytarr(5,5,5)
;  arr[1,1,1]=1
;  arr[1,3,1]=1
;  arr[3,1,1]=1
;  arr[3,3,1]=1
;  arr[2,2,2]=1
;  SHADE_VOLUME, arr,0.9, v1, p1
;
;  help,v1
;  help,p1
;  
;  numberVertices = MESH_DECIMATE(v1, p1, p2, VERTICES = v2)
;  help,v2
;  help,p2
;  
;  
;  writestl,file,v2,p2
end

pro readstl,file,point,conn
  ;file='D:\flower.stl'

  openr,lun,file,/get_lun
  hh=''
  readf,lun,hh 
  readf,lun,hh
  str=strtrim(hh,1)
  isascii=strcmp(str,'facet',5)
  if isascii eq 1 then begin
    print,'ASCII format'
    readascii,file,point,conn
  endif else begin
    print,'binary format'
    readbinary,file,point,conn
  endelse

end




pro readstl_old,inputfile,point,conn
file=inputfile
  openr,file_lun,file,/get_lun

  filename=indgen(80,/byte)
  readu,file_lun,filename
  ;point_lun,file_lun,80
  ;facennum=indgen(1,/long)
  facenum=lonarr(1)
  readu,file_lun,facenum
  ;eachtriangle=fltarr(3,facenum*4)
  normal=fltarr(3,1)
  readpoint=fltarr(3,3)
  k=0

  readu,file_lun,normal
  readu,file_lun,readpoint
  ;print,readpoint
  point=readpoint
  Attributebytecountend=uintarr(1)
  readu,file_lun,Attributebytecountend

  ;print,Attributebytecountend
  print,'111'
  for i=1l,facenum[0]-1 do begin
    readu,file_lun,normal
    readu,file_lun,readpoint
    ;print,readpoint
    point=[[point],[readpoint]]
    readu,file_lun,Attributebytecountend
  endfor
  help,point
  sizeponint=size(point)
  print,sizeponint[2]
  p=lonarr(sizeponint[2])
  for m=0l,sizeponint[2]-1 do begin
    p[m]=m
  endfor
  ;print,p
sizeconn=long(n_elements(p)/3+n_elements(p))
conn=lonarr(sizeconn)
j=0l
i=0l
 for q=0l,long(facenum[0])-1  do begin
    conn[j]=3l
    conn[j+1]=p[i]
    conn[j+2]=p[i+1]
    conn[j+3]=p[i+2]

    j=j+4
    i=i+3

  endfor
;  oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], point, POLYGON=conn)
;  mymodel1 = OBJ_NEW('IDLgrModel')
;  mymodel1->Add, oPoly1
;  ;
;  XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]
   free_lun, file_lun

end

Pro readascii,inputfile,point,conn
  file=inputfile
  openr,lun,file,/get_lun
  hh=''

  list=LIST()
  list2=LIST()
  while(~eof(lun))do begin
    readf,lun,hh
    if (strpos(hh, 'vertex')) ge 0 then begin
      list2.add,hh
    endif
    list.add,hh
  endwhile
  Result = list.Count( list)
  Result2 = N_ELEMENTS(list2)
  point=fltarr(3,Result2)
  print,Result2
  facenum=long(Result2/3)
  p=lonarr(Result2)
  for m=0l, Result2-1 do begin
    p[m]=m
    split=STRSPLIT(list2[m],/EXTRACT)
    point[0,m]=split[1]
    point[1,m]=split[2]
    point[2,m]=split[3]
  endfor

  j=0l
  i=0l
  sizeconn=long(Result2/3+Result2)
  conn=lonarr(sizeconn)
  for q=0l,long(facenum)-1  do begin
    conn[j]=3l
    conn[j+1]=p[i]
    conn[j+2]=p[i+1]
    conn[j+3]=p[i+2]

    j=j+4
    i=i+3

  endfor


  ;oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], point, POLYGON=conn)
  ;mymodel1 = OBJ_NEW('IDLgrModel')
  ;mymodel1->Add, oPoly1
  ;
  ;XOBJVIEW, mymodel1, BACKGROUND = [0, 0, 0]




  free_lun,lun


end




pro readbinary,inputfile,newpoint,conn
  file=inputfile
  openr,file_lun,file,/get_lun

;  header_info = STRARR(1)
;  readu,file_lun,header_info
;  print, 'header: ', header_info
  
  point_lun,file_lun,80
  ;facennum=indgen(1,/long)
  facenum=lonarr(1)
  readu,file_lun,facenum
  print,'long ', facenum
  ;eachtriangle=fltarr(3,facenum*4)
  normal=fltarr(3,1)
  readpoint=fltarr(3,3)
  
;  print, 'ulong ', ulong(facenum)
;  print, 'long64 ', long64(facenum)
  
  ;  k=0
  ;
  ;  readu,file_lun,normal
  ;  readu,file_lun,readpoint
  ;print,readpoint
  ;point=fltarr(3,3)
  Attributebytecountend=uintarr(1)
  ;  readu,file_lun,Attributebytecountend

  ;print,Attributebytecountend
  print,'start reading stl file, patches # = ', facenum
  j=0l
  sizeconn=(facenum)*4
  conn=lonarr(sizeconn)

  ;thisFormat='(12F0)'
  array1 = Fltarr(25, facenum[0]/2-1)
  Readu, file_lun, array1;, thisFormat
  print, array1[*,0:5]
  print, 'finish'
  
  point_lun,file_lun,96+38*1
  array2 = Fltarr(25, facenum[0]/2-1)
  Readu, file_lun, array2;, thisFormat
  print, array2[*,0:5]
  
  help, array1
  help, array2
  point = FLTARR(9, facenum[0])
  FOR i=0l,facenum[0]-6 do begin
;    PRINT, 'i = ', i
    IF (i mod 2) EQ 0 THEN point[*, i] = array1[3:11, i/2] ELSE point[*, i] = array2[3:11, i/2]
    k=3*i
    conn[j]=3l
    conn[j+1]=k
    conn[j+2]=k+1
    conn[j+3]=k+2
    j=j+4


  ENDFOR
  
  newpoint=REFORM(point, 3, facenum[0]*3)
  
;  print, newpoint[*, 0:5]
;  
  newpoint[0, *] = newpoint[0, *] - MEDIAN(newpoint[0, *])
  newpoint[1, *] = newpoint[1, *] - MEDIAN(newpoint[1, *])
  newpoint[2, *] = newpoint[2, *] - MEDIAN(newpoint[2, *])
;  
;  print, 'new', newpoint[*, 0:5]
;  print, MEAN(newpoint[0, *])
;  for i=1l,facenum[0]/2-2, 2 do array1[*, i] = array2[*, i-1]
;  print, 'finish2'
;
;    
;    IF (i mod 2) EQ 0 THEN BEGIN
;      IF i EQ 0 THEN BEGIN
;        point=REFORM(array1[3:11, i], 3, 3)
;      ENDIF ELSE BEGIN
;        point=[[point],[REFORM(array2[3:11, i], 3, 3)]]
;      ENDELSE
;    ENDIF ELSE BEGIN
;      point=[[point],[REFORM(array2[3:11, i-1], 3, 3)]]
;    ENDELSE
;    k=3*i
;    conn[j]=3l
;    conn[j+1]=k
;    conn[j+2]=k+1
;    conn[j+3]=k+2
;    j=j+4
;
;    count=(i+1.)/facenum[0]*100.0
;    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
;;    print, 'i = ', i    
;
;  
;;  for i=0l,5-1 do begin
;;    ;    readu,file_lun,normal
;;    point_lun,file_lun,96+50*i
;;    readu,file_lun,readpoint
;;    ;print,readpoint
;;    if (i eq 0l ) then begin
;;      point=readpoint
;;    endif else begin
;;      point=[[point],[readpoint]]
;;    endelse
;;    
;;    print, 'readpoint = ', readpoint
;;    count=(i+1.)/facenum[0]*100.0
;;    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
;;
;;    ; print,point
;;    ;    readu,file_lun,Attributebytecountend
;;    k=3*i
;;    conn[j]=3l
;;    conn[j+1]=k
;;    conn[j+2]=k+1
;;    conn[j+3]=k+2
;;    j=j+4
;;
;;;    print,'i',i
;;
;;  endfor

  free_lun, file_lun

  ; help,point
  ;  sizeponint=size(point)
  ;  print,sizeponint[2]
  ;  p=lonarr(sizeponint[2])
  ;  for m=0l,sizeponint[2]-1 do begin
  ;    p[m]=m
  ;  endfor
  ;  print,'p',p

  ;  conn=lonarr(sizeconn)
  ;  j=0l
  ;  i=0l
  ;  for q=0l,long(facenum[0])-1  do begin
  ;    conn[j]=3l
  ;    conn[j+1]=p[i]
  ;    conn[j+2]=p[i+1]
  ;    conn[j+3]=p[i+2]
  ;
  ;    j=j+4
  ;    i=i+3
  ;
  ;  endfor
  ;print,'conn',conn



end
pro  trai,input
  a=1
  for i=0,N_elements(input)-2 do begin
    a=(input[i] eq 3)?1:0
    if (a eq 0) then begin

      print,'i',i

    endif
    i=i+3
  endfor
  i=0
end
