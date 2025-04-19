
pro browser_convert_bmp, afn,  afn_des, tagarr, tagsNumber, valuearr
COMMON counter,counter
;afn is the .dcm filename
; Include IDL system procedures and functions when profiling:

;PROFILER, /SYSTEM


;;;DICOM Files Section
fnstar=STRTRIM(afn)+'*'
a=findfile(fnstar)
;print, 'a=(L12)', a


;Temp_value_arr=
;@m.bat to run Klaus' codes, we will call a procedure: fdecomp
FOR i=0L,N_ELEMENTS(a)-1 DO BEGIN

    ;Call procedure in Klaus' code
    fdecomp, a[i], disk, dir, name, qual, version

    IF STRLEN(name) EQ 0 THEN BEGIN

       ;efg is the last two digits of a a[i]
      efg= STRMID(a[i],STRLEN(a[i])-2,2)
      IF (efg NE ".\") THEN BEGIN
         ; the a[i] is not the current directory .\ or ..\

         ;a[i] is a directory, do recursive calling
         browser_convert_bmp,a[i], afn_des, tagarr, tagsNumber, valuearr
       ENDIF
    ENDIF ELSE BEGIN
    ;a[i] is a filename, not a directory
       d=STRMID(a[i],STRLEN(a[i])-3,3)
       ;IF STRLEN(qual) NE 0 THEN BEGIN
         counter=counter+1
       ;ENDIF ELSE BEGIN
       ; FILE_COPY, a[i], afn_des;"E:\images\AEE_Novartis\AEE_Bruno_Morgan\SUTTON_04_JOHN_4525446\dyn_output\"
       IF QUERY_DICOM(a[i]) THEN BEGIN
         obj = OBJ_NEW('IDLffDICOM', a[i])
         result=QUERY_DICOM(a[i], info)
         IF result NE 0 THEN BEGIN
;          obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
;          order=1
;          vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
          ; Get the image data
          array=obj->GetValue('7fe0'x, '0010'x)
          vPixels=*array[0]
          PTR_FREE, array
          ; select monochrome image
          ;select the last sample data (need to be added here).
;          IF samples gt 1 THEN BEGIN
;              pixeldata=vPixels
;;              a=dialog_message('multiple sample per pixel, add codes here, line 2129, jia.11@osu.edu. file_io.pro')
;          ENDIF ELSE BEGIN
;              pixeldata=vPixels
;          ENDELSE
          image_temp=vPixels
  
          ;; Series Number
          series_number = STRING(*(obj->GetValue('0020'x,'0011'x,/NO_COPY))[0])
;         IF (obj->QueryValue('0020,0011')) EQ 2 THEN series_number = obj->GetValue('0020,0011') ELSE series_number = ''
  
          ;create the folder structure to copy the bmp images
          dirstring=STRSPLIT(dir, '\',/EXTRACT)
          afn_des_sub=afn_des+STRCOMPRESS(dirstring[N_ELEMENTS(dirstring)-1],/REMOVE_ALL)+'\'
          FILE_MKDIR, afn_des_sub
          new_filename=afn_des_sub+name+".jpg"
          s=size(image_temp)
          IF S[0] EQ 2 THEN BEGIN
              window, 0, XSIZE=s[1], YSIZE=s[2]
              tvscl,image_temp, /ORDER
              new_image = TVRD(ORDER=0,TRUE=1)
              WRITE_JPEG, new_filename, new_image, TRUE=1;, /RGB;"E:\images\AEE_Novartis\AEE_Bruno_Morgan\SUTTON_04_JOHN_4525446\dyn_output\"
              WDELETE, 0
             ;WRITE_BMP, new_filename, image_temp
          ENDIF ELSE BEGIN
              IF s[0] EQ 3 THEN BEGIN
                  window, 0, XSIZE=s[2], YSIZE=s[3]
                  tv, image_temp,  TRUE=1, /ORDER
                  new_image = TVRD(ORDER=0,TRUE=1)
                  WRITE_JPEG, new_filename, new_image, TRUE=1;, /RGB;"E:\images\AEE_Novartis\AEE_Bruno_Morgan\SUTTON_04_JOHN_4525446\dyn_output\"
                  WDELETE, 0
                 ;WRITE_BMP, new_filename, image_temp
             ENDIF
          ENDELSE
  
         ENDIF
        OBJ_DESTROY, obj
      ENDIF
    ENDELSE
ENDFOR

; Retrieve the profiling results:
;PROFILER, /REPORT

end

;
; ��GJ $Id: convert_to_bmp_images.pro, 2004/03/21 11:19: Guang Jia $
;
;
;+
; NAME:
;   convert_to_bmp_images
;
; PURPOSE:
;   This procedure converts the dicom images to bitmap images with their orginal size and image display order
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;   convert_to_bmp_images
;
; OPTIONAL INPUTS:
;
;
;
; MODIFICATION HISTORY:
;   Written by:  GJ, March, 2004
;-


;
PRO convert_to_bmp_images
COMMON counter,counter

counter=0
afn= DIALOG_PICKFILE(/MUST_EXIST,/DIRECTORY, TITLE="Dicom Image Directory:")
IF STRLEN(afn) NE 0 THEN BEGIN
    afn_des= DIALOG_PICKFILE(PATH=afn,/MUST_EXIST,/DIRECTORY, TITLE="Destination Directory:")
    IF STRLEN(afn_des) NE 0 THEN BEGIN
        browser_convert_bmp,afn, afn_des, tagarr, TagsNumber, valuearr
    ENDIF
ENDIF
END