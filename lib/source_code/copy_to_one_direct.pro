
pro browse_copy, afn,  afn_des, tagarr, tagsNumber, valuearr
COMMON counter,counter
;afn is the .dcm filename
; Include IDL system procedures and functions when profiling:

;PROFILER, /SYSTEM


;;;DICOM Files Section
fnstar=STRTRIM(afn)+'*'
a=findfile(fnstar)
;print, 'a=(L12)', a

; Create the progress bar.
progressbar = Obj_New('progressbar', Color='Chocolate', Text='0%', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
; Place the progress bar on the display.
progressbar -> Start

;Temp_value_arr=
;@m.bat to run Klaus' codes, we will call a procedure: fdecomp
FOR i=0L,N_ELEMENTS(a)-1 DO BEGIN

    count=(i+1.)/N_ELEMENTS(a)*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;Call procedure in Klaus' code
    fdecomp, a[i], disk, dir, name, qual, version

    IF STRLEN(name) EQ 0 THEN BEGIN

       ;efg is the last two digits of a a[i]
      efg= STRMID(a[i],STRLEN(a[i])-2,2)
      IF (efg NE ".\") THEN BEGIN
         ; the a[i] is not the current directory .\ or ..\

         ;a[i] is a directory, do recursive calling
         browse_copy,a[i], afn_des, tagarr, tagsNumber, valuearr
       ENDIF
    ENDIF ELSE BEGIN
    ;a[i] is a filename, not a directory
       dirarr=STRSPLIT(dir, '\', /EXTRACT)
       fn=dirarr[N_ELEMENTS(dirarr)-1]
       fnarr=STRSPLIT(fn, '_', /EXTRACT)
       new_fn=STRJOIN([fnarr[0:(N_ELEMENTS(fnarr)-1)], name],'_')
       fdecomp, afn_des, disk1, dir1, name1, qual1, version
       fcomp, fn_des, disk1, dir1, new_fn, qual, version
       d=STRMID(a[i],STRLEN(a[i])-3,3)
       FILE_COPY, a[i], fn_des
    ENDELSE
ENDFOR

progressbar -> Destroy
; Retrieve the profiling results:
;PROFILER, /REPORT

end

;
; ��GJ $Id: ADD_DCM_TYPENAME_v01.pro, 2003/08/05 11:00:26 Guang Jia $
;
;
;+
; NAME:
;   ADD_DCM_TYPENAME_v01
;
; PURPOSE:
;   This procedure add ".dcm" to the files under a directory. e.g. the data from
;   Houston (2003/08/05).
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;   ADD_DCM_TYPENAME_v01
;
; OPTIONAL INPUTS:
;
;
;
; MODIFICATION HISTORY:
;   Written by:  GJ, August 2003
;-


;version 01 for Houston data
PRO copy_to_one_direct
COMMON counter,counter

counter=0
afn= DIALOG_PICKFILE(/MUST_EXIST,/DIRECTORY, TITLE='Source Directory:')
IF (STRLEN(afn) NE 0) THEN BEGIN
select_des:    afn_des= DIALOG_PICKFILE(PATH=afn, /MUST_EXIST,/DIRECTORY, TITLE="Destination Directory:")
    IF (STRLEN(afn_des) NE 0) THEN BEGIN
        afn_len=STRLEN(afn)
        IF STRCMP(afn, afn_des, afn_len, /FOLD_CASE) THEN BEGIN
            afnError=DIALOG_MESSAGE('Dest folder cannot be the same as source folder or subfolder! Please reselect dest folder.', /ERROR)
            GOTO, select_des
        ENDIF

        browse_copy,afn, afn_des, tagarr, TagsNumber, valuearr
    ENDIF
ENDIF
END