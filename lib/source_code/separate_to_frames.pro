
FUNCTION STRCOMPRESSJ, stringa

words=STRSPLIT(stringa, ' ',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, ',',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '-',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '\',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '/',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, ':',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '*',/EXTRACT)
stringa=STRJOIN(words, 'star')

words=STRSPLIT(stringa, '?',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '"',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '<',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '>',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '^',/EXTRACT)
stringa=STRJOIN(words, '')

words=STRSPLIT(stringa, '|',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '.',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '[',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, ']',/EXTRACT)
stringa=STRJOIN(words, '_')

words=STRSPLIT(stringa, '_',/EXTRACT)
stringa=STRJOIN(words, '_')

return, stringa
END

pro browse_separate, afn,  afn_des, afn_des_output
COMMON counter,counter
;afn is the .dcm filename
; Include IDL system procedures and functions when profiling:

;PROFILER, /SYSTEM
;Establish error handler
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    Error_status = 0
    IF i LT N_ELEMENTS(a)-1 THEN BEGIN
        i=i+1
        GOTO, Next
    ENDIF ELSE BEGIN
      ;2019/1/29 GJ, remove the bug
;        OBJ_DESTROY, obj
        progressbar -> Destroy
        return
    ENDELSE
ENDIF

;;;DICOM Files Section
fnstar=STRTRIM(afn)+'*'
a=findfile(fnstar)
;print, 'a=(L12)', a

; Create the progress bar.
progressbar = Obj_New('progressbar', Color='Tomato', Text='0%', Title='Sorting...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
; Place the progress bar on the display.
progressbar -> Start

FOR i=0L,N_ELEMENTS(a)-1 DO BEGIN

Next:    count=(i+1.)/N_ELEMENTS(a)*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

    ;Call procedure in Klaus' code
    fdecomp, a[i], disk, dir, name, qual, version
    
    ;@GJ, 2023/5/18, Initialize the modified file name
    fn_MPI_new = ''
    
    IF STRLEN(name) EQ 0 THEN BEGIN
        ;efg is the last two digits of a a[i]
        efg= STRMID(a[i],STRLEN(a[i])-2,2)
         ; the a[i] is not the current directory .\ or ..\
        ;a[i] is a directory, do recursive calling
        result=FILE_INFO(a[i])
        IF (efg NE '.\') AND (result.directory EQ 1) THEN browse_separate, a[i], afn_des, afn_des_output
    ENDIF ELSE BEGIN
        ;@GJ, 2023/5/17, Try adding meta type file group length 
        IF QUERY_DICOM(a[i]) EQ 0 THEN BEGIN
          data_1 = bytarr(158)
          data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]
          
          fn_MPI = a[i]
          data_2 = read_binary(fn_MPI, DATA_START=132)

          ; Use Count to get the number of nonzero elements:
          index = WHERE(data_2 EQ 8, count)
          ; Only subscript the array if it is safe:
          ;define the File Meta Group Element Length 0002,0000
          IF count NE 0 THEN BEGIN
            FOR j=0, count-1 DO BEGIN
              IF data_2[index[j]+1] EQ 0 THEN BEGIN
                data_1[140]=BYTE(index[j])+14B
                break
              ENDIF
            ENDFOR
          ENDIF
          
          ;@GJ, save the modified dicom file         
          pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
          cur_dir = STRMID(afn_des, 0, pos_dir+1)
          fn_MPI_new = cur_dir+FILE_BASENAME(a[i])
          OPENW, U, fn_MPI_new,/GET_LUN
          WRITEU, U, data_1
          WRITEU, U, data_2
          FREE_LUN, U
;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
          ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
          a[i]=fn_MPI_new
        ENDIF
        
        
        ;a[i] is a filename, not a directory
        IF QUERY_DICOM(a[i]) THEN BEGIN
                       
            counter=counter+1

            obj = OBJ_NEW('IDLffDICOM', a[i])
            
               ;;Patient Name
               IF PTR_VALID(obj->GetValue('0010'x,'0010'x,/NO_COPY)) THEN patient_name = STRING(*(obj->GetValue('0010'x,'0010'x,/NO_COPY))[0]) ELSE patient_name='test'
               IF STRLEN(patient_name) GT 4 THEN patient_name=STRMID(patient_name, 0, 4)
;               IF (obj->QueryValue('0010,0010')) EQ 2 THEN patient_name = obj->GetValue('0010,0010') ELSE patient_name= ''

               ;;Patient ID
               IF PTR_VALID(obj->GetValue('0010'x,'0020'x,/NO_COPY)) THEN patient_ID = STRING(*(obj->GetValue('0010'x,'0020'x,/NO_COPY))[0]) ELSE patient_ID='test'
               IF STRLEN(patient_ID) GT 4 THEN patient_ID=STRMID(patient_ID, 0, 4)
               ;patient_ID = STRING(*(obj->GetValue('0010'x,'0020'x,/NO_COPY))[0])
;               IF (obj->QueryValue('0010,0020')) EQ 2 THEN patient_ID = obj->GetValue('0010,0020') ELSE patient_ID = ''

               ;;Study description
               IF PTR_VALID(obj->GetValue('0008'x,'1030'x,/NO_COPY)) THEN study_descr = STRING(*(obj->GetValue('0008'x,'1030'x,/NO_COPY))[0]) ELSE study_descr='test'
               IF STRLEN(study_descr) GT 4 THEN study_descr=STRMID(study_descr, 0, 4)
               ;study_descr = STRING(*(obj->GetValue('0008'x,'1030'x,/NO_COPY))[0])
;              IF (obj->QueryValue('0008,1030')) EQ 2 THEN study_descr = obj->GetValue('0008,1030') ELSE study_descr = ''

               ;;Series Date
               IF PTR_VALID(obj->GetValue('0008'x,'0020'x,/NO_COPY)) THEN study_date = STRING(*(obj->GetValue('0008'x,'0020'x,/NO_COPY))[0]) ELSE study_date='20230517'
               ;study_date = STRING(*(obj->GetValue('0008'x,'0020'x,/NO_COPY))[0])
;              IF (obj->QueryValue('0008,0020')) EQ 2 THEN study_date = obj->GetValue('0008,0020') ELSE study_date = ''

                ;; Series Number
                IF PTR_VALID(obj->GetValue('0020'x,'0011'x,/NO_COPY)) THEN series_number = STRING(*(obj->GetValue('0020'x,'0011'x,/NO_COPY))[0]) ELSE series_number='S001'
                ;series_number = STRING(*(obj->GetValue('0020'x,'0011'x,/NO_COPY))[0])
;                IF (obj->QueryValue('0020,0011')) EQ 2 THEN series_number = obj->GetValue('0020,0011') ELSE series_number = ''

               ;; Series Description
               IF PTR_VALID(obj->GetValue('0008'x,'103e'x,/NO_COPY)) THEN series_descr = STRING(*(obj->GetValue('0008'x,'103e'x,/NO_COPY))[0]) ELSE series_descr='test'
               IF STRLEN(series_descr) GT 4 THEN series_descr=STRMID(series_descr, 0, 4)
               ;series_descr = STRING(*(obj->GetValue('0008'x,'103e'x,/NO_COPY))[0])
;                IF (obj->QueryValue('0008,103e')) EQ 2 THEN series_descr = obj->GetValue('0008,103e') ELSE series_descr = ''

               ;; Series Time
               IF PTR_VALID(obj->GetValue('0008'x,'0031'x,/NO_COPY)) THEN series_time = STRING(*(obj->GetValue('0008'x,'0031'x,/NO_COPY))[0]) ELSE series_time='000000'
               ;series_time = STRING(*(obj->GetValue('0008'x,'0031'x,/NO_COPY))[0])
;                IF (obj->QueryValue('0008,0031')) EQ 2 THEN series_time = obj->GetValue('0008,0031') ELSE series_time = ''
              ;
              ; P modal
               p_modal_temp = STRING(*(obj->GetValue('0008'x,'0060'x,/NO_COPY))[0])
               p_modal = STRCOMPRESS(p_modal_temp, /REMOVE_ALL)
               
               ;@GJ,2023/5/24, adding the modality
                new_des_sub1=p_modal+'_'+STRCOMPRESS(STRING(BYTE(STRING(patient_ID)))+'_'+STRING(BYTE(patient_name))+'_'+STRING(study_date)+'_'+STRING(study_descr))
                new_des_sub2=STRCOMPRESS('S'+STRING(series_number, FORMAT='(I5.4)')+'_'+STRMID(STRING(series_time), 0, 6)+'_'+STRING(series_descr)+'\',/REMOVE_ALL) ;
                new_des_sub1=STRCOMPRESSJ(new_des_sub1)
                new_des_sub2=STRCOMPRESSJ(new_des_sub2)
                afn_des_sub=afn_des+new_des_sub1+STRMID(afn_des,0,1,/REVERSE_OFFSET)+new_des_sub2

                IF STRLEN(afn_des_sub) GT 240 THEN BEGIN
                    new_des_sub1=STRCOMPRESS(STRING(BYTE(patient_name))+'_'+STRING(study_date))
                    afn_des_sub=afn_des+new_des_sub1+STRMID(afn_des,0,1,/REVERSE_OFFSET)+new_des_sub2
                ENDIF

                ;check the existance of directory
                dircheck=FILE_INFO(afn_des_sub)
                IF dircheck.exists NE 1 THEN BEGIN
                    FILE_MKDIR, afn_des_sub
                    afn_des_output=[afn_des_output, afn_des_sub]
                ENDIF
                FILE_COPY, a[i], afn_des_sub

                fdecomp, afn_des_sub+'\'+name, disk1, dir1, name1, qual1, version
                fcomp, newA, disk1, dir1, name, qual, version
                IF STRCMP(STRMID(newA, 0, 1, /REVERSE_OFFSET),'.') THEN newA=STRMID(newA, 0, STRLEN(newA)-1)
                IF (STRCMP(qual, 'dcm', /FOLD_CASE) EQ 0) THEN BEGIN
                    IF (file_info(newA+'.dcm')).exists EQ 0 THEN BEGIN
                      FILE_MOVE, newA, newA+'.dcm', /OVERWRITE
                    ENDIF ELSE BEGIN
                      systm=STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(f12.2)')
                      FILE_MOVE, newA, newA+systm+'.dcm', /OVERWRITE
                    ENDELSE
                      newA=newA+'.dcm'
                ENDIF

                fdecomp, newA, disk2, dir2, name2, qual2, version
                ; Image Number
                IF PTR_VALID(obj->GetValue('0020'x,'0013'x,/NO_COPY)) THEN value = STRING(*(obj->GetValue('0020'x,'0013'x,/NO_COPY))[0]) ELSE value = 0
                ;value = DOUBLE(*(obj->GetValue('0020'x,'0013'x,/NO_COPY))[0])
;                IF (obj->QueryValue('0020,0013')) EQ 2 THEN value = obj->GetValue('0020,0013') ELSE value = 0
                image_number= MAX(value)
               ; Rename the file by adding image number if there is no image number
                imageNStr=STRCOMPRESS(STRING(image_number, FORMAT='(I5.5)'),/REMOVE_ALL)
                result=STRCMP(name2, imageNStr, 5, /FOLD_CASE)
                IF result EQ 0 THEN BEGIN
                    Newname2=imageNStr+'_'+name2
                    fcomp, Newai, disk2, dir2, Newname2, qual2, version
                    result=FILE_INFO(Newai)
                    IF result.exists EQ 0 THEN BEGIN
                      FILE_MOVE, newA, Newai, /OVERWRITE
                    ENDIF ELSE BEGIN
                      systm=STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(I12)')
                      FILE_MOVE, newA, Newai+systm, /OVERWRITE
                    ENDELSE
                    
                ENDIF

        ENDIF
        OBJ_DESTROY, obj
        IF (FILE_INFO(fn_MPI_new)).directory EQ 0 AND (FILE_INFO(fn_MPI_new)).exists EQ 1 THEN FILE_DELETE, fn_MPI_new, /QUIET
    ENDELSE
ENDFOR


progressbar -> Destroy

; Retrieve the profiling results:
;PROFILER, /REPORT

end



pro browse_separate_by_volume, afn
COMMON counter,counter
;afn is the .dcm filename
; Include IDL system procedures and functions when profiling:

;PROFILER, /SYSTEM
;PROFILER, /SYSTEM
;Establish error handler
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    Error_status = 0
    IF i LT N_ELEMENTS(a)-1 THEN BEGIN
        i=i+1
        GOTO, Next
    ENDIF ELSE BEGIN
        progressbar -> Destroy
        return
    ENDELSE
ENDIF

;;;DICOM Files Section
fnstar=STRTRIM(afn)+'*'
a=findfile(fnstar)
;print, 'a=(L12)', a

; Create the progress bar.
progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Separating by volume...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
; Place the progress bar on the display.
progressbar -> Start

FOR i=0L,N_ELEMENTS(a)-1 DO BEGIN

    ;Call procedure in Klaus' code
Next:    fdecomp, a[i], disk, dir, name, qual, version

    IF STRLEN(name) EQ 0 THEN BEGIN
        ;efg is the last two digits of a a[i]
        efg= STRMID(a[i],STRLEN(a[i])-2,2)
         ; the a[i] is not the current directory .\ or ..\
        ;a[i] is a directory, do recursive calling
        result=FILE_INFO(a[i])
        IF (efg NE '.\') AND (result.directory EQ 1) THEN browse_separate_by_volume, a[i]
    ENDIF ELSE BEGIN
        obj = OBJ_NEW('IDLffDICOM', a[i]);, /NO_PIXEL_DATA)
        OBJ_DESTROY, obj
        ;a[i] is a filename, not a directory
        IF QUERY_DICOM(a[i]) THEN BEGIN
            dummy=browse_directory(a[i],fnarr_in,tError)
            result=sort_DCM_by_ser_and_ima_number(fnarr_in, fnarr_out, numSlices, tError)
            IF (result EQ 0) THEN BEGIN
                numTimePT=N_ELEMENTS(fnarr_out)/numSlices
                IF (numTimePT GT 1) THEN BEGIN
                    ATagARR=STRARR(numTimePT)
                    FOR j=0, numTimePT-1 DO BEGIN

                        count=(j+1.)/numTimePT*100.0
                        progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

                        temp=STRLEN(STRCOMPRESS(STRING(numTimePT),/REMOVE_ALL))
                        msx=STRCOMPRESS(STRING(temp), /REMOVE_ALL)
                        int1='('+'I'+msx+'.'+msx+')'
                        ;afn_sub=afn+'V'+STRING(j+1, FORMAT=int1)+'\'
                        ;print out trigger time
                        FOR kl=0, numSlices-1 DO BEGIN
                            obj = OBJ_NEW('IDLffDICOM', fnarr_out[j*numSlices+kl])
                            temp = STRING(*(obj->GetValue('0018'x,'0024'x,/NO_COPY))[0])
                            tempTag = STRCOMPRESSJ(temp)
                            IF kl EQ 0 THEN ATag=STRING(tempTag) ELSE ATag=ATag+','+STRING(tempTag)
                            OBJ_DESTROY, obj
                        ENDFOR
                        ATagARR[j]='volume_'+ STRING(j+1, FORMAT=int1) +": " + ATag
                        print, ATagARR[j]
                        afn_sub=afn+'V'+STRING(j+1, FORMAT=int1)+'_'+STRCOMPRESS(tempTag, /REMOVE_ALL)+'\'
                        FILE_MKDIR, afn_sub
                        FILE_MOVE, fnarr_out[(j*numSlices):((j+1)*numSlices-1)], afn_sub
                    ENDFOR
                    ;2005 Jun 27, for Xiangyu to output trigger time
                    ;SaveResultfn=DIALOG_PICKFILE(TITLE='save Trigger Time in every volume',PATH=afn_sub,FILTER='*.txt')
                    ;IF STRLEN(SaveResultfn) NE 0 THEN write_triggertime, SaveResultfn, TTARR
                ENDIF
            ENDIF
        ENDIF
        BREAK
    ENDELSE
ENDFOR

progressbar -> Destroy

; Retrieve the profiling results:
;PROFILER, /REPORT

end

;TTarray is trigger time array
PRO write_triggertime, filename, TTarray

N=N_ELEMENTS(TTarray)
openw,unit,filename,/get_lun

printf,unit,STRCOMPRESS(STRING(N),/REMOVE_ALL)+' volumes: trigger times',FORMAT='(%"%s")'
FOR i=0, N-1 DO printf,unit,TTarray[i]

;; close file and END execution
FREE_LUN, unit
END


;
; ��GJ $Id: Separate images by series, 2003/08/05 11:00:26 Guang Jia $
;
;
;+
; NAME:
;   Separate_to_series
;
; PURPOSE:
;   This procedure separates images based on their series number
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;   Separate_to_series
;
; OPTIONAL INPUTS:
;
;
;
; MODIFICATION HISTORY:
;   Written by:  GJ, August 2003
;   Feb, 13, 2004, GJ (jia.11@osu.edu)
;   Jan, 12, 2005, GJ (jia.11@osu.edu) separate the images from series description
;   April, 20, 2005, GJ (jia.11@osu.edu) Jun suggested to further separate the images by volume within a series
;-


;version 01 for Houston data
PRO Separate_to_series, afn, afn_des, afn_des_output = afn_des_output
COMMON counter,counter

counter=0
;afn= DIALOG_PICKFILE(/MUST_EXIST,/DIRECTORY, TITLE="Source Directory:")
IF STRLEN(afn) NE 0 THEN BEGIN
select_des: IF N_ELEMENTS(afn_des) EQ 0 THEN afn_des= DIALOG_PICKFILE(PATH=afn, /MUST_EXIST,/DIRECTORY, TITLE="Destination DCM file Directory:")
    IF STRLEN(afn_des) NE 0 THEN BEGIN
        afn_len=STRLEN(afn)
        IF STRCMP(afn, afn_des, afn_len, /FOLD_CASE) THEN BEGIN
            afnError=DIALOG_MESSAGE('Dest folder cannot be the same as source folder or subfolder! Please reselect dest folder.', /ERROR)
            GOTO, select_des
        ENDIF

       afn_des_output=''
        browse_separate,afn, afn_des, afn_des_output

        ;text="Separation by Series Done. Continue: separate by volumes within one series?"
        ;result=DIALOG_MESSAGE(text, /QUESTION, TITLE="continue: YES or NO")
        ;IF STRCMP(result, 'YES', /FOLD_CASE) THEN BEGIN
        ;    FOR i=1, N_ELEMENTS(afn_des_output)-1 DO browse_separate_by_volume, afn_des_output[i]
        ;ENDIF
    ENDIF
ENDIF

END


;
; ��GJ $Id: Separate images by series, 2003/08/05 11:00:26 Guang Jia $
;
;
;+
; NAME:
;   Separate_to_frames
;
; PURPOSE:
;   This procedure separates images based on their frame
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;   Separate_to_frames
;
; OPTIONAL INPUTS:
;
;
;
; MODIFICATION HISTORY:
;   Written by:  GJ, August 2003
;   Feb, 13, 2004, GJ (jia.11@osu.edu)
;   Jan, 12, 2005, GJ (jia.11@osu.edu) separate the images from series description
;   April, 20, 2005, GJ (jia.11@osu.edu) Jun suggested to further separate the images by volume within a series
;-


;version 01 for Houston data
PRO Separate_to_frames
COMMON counter,counter

counter=0
afn= DIALOG_PICKFILE(/MUST_EXIST,/DIRECTORY, TITLE="Source Directory:")
IF STRLEN(afn) NE 0 THEN BEGIN
    browse_separate_by_volume, afn
ENDIF

END
