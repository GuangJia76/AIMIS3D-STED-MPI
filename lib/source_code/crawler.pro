PRO Tags_Define, tagarr_MR, descrarr_MR, Notes_MR, tagarr_PETCT, descrarr_PETCT, Notes_PETCT

;Define the tags you want to search
TagsNumber = 24L
tagarr = lonarr(TagsNumber,2)
descrarr = strarr(TagsNumber + 8)
tagarr[0,0] = '0010'x
tagarr[0,1] = '0010'x ;Patient name
descrarr[0] = "PatientName"
tagarr[1,0] = '0010'x
tagarr[1,1] = '0020'x ;Patient ID
descrarr[1] = "PatientID"
tagarr[2,0] = '0010'x
tagarr[2,1] = '0030'x ;Patient date of birth
descrarr[2] = "PatientDOB"
tagarr[3,0] = '0010'x
tagarr[3,1] = '0040'x ;Patient sex
descrarr[3] = "PatientSex"
tagarr[4,0] = '0010'x
tagarr[4,1] = '1010'x ;Patient age
descrarr[4] = "PatientAge"
tagarr[5,0] = '0010'x
tagarr[5,1] = '1030'x ;Patient weight
descrarr[5] = "PatientWeight"
tagarr[6,0] = '0008'x
tagarr[6,1] = '0023'x ;Study date
descrarr[6] = "StudyDate"
tagarr[7,0] = '0008'x
tagarr[7,1] = '1030'x ;Study Description
descrarr[7] = "StudyDesc"
tagarr[8,0] = '0020'x
tagarr[8,1] = '0011'x ;Series number
descrarr[8] = "SE#"
tagarr[9,0] = '0008'x
tagarr[9,1] = '103e'x ;Series Description
descrarr[9] = "SeDesc"
tagarr[10,0] = '0008'x
tagarr[10,1] = '0031'x ;Series time
descrarr[10] = "Time"
tagarr[11,0] = '0020'x
tagarr[11,1] = '0010'x ;Study ID
descrarr[11] = "StudyID"
tagarr[12,0] = '0018'x
tagarr[12,1] = '0080'x ;Reptition Time
descrarr[12] = "TR"
tagarr[13,0] = '0018'x
tagarr[13,1] = '0081'x ;Echo Time
descrarr[13] = "TE"
tagarr[14,0] = '0018'x
tagarr[14,1] = '1314'x ;Flip Angle
descrarr[14] = "Flip"
tagarr[15,0] = '0018'x
tagarr[15,1] = '0088'x ;Spacing between slices
descrarr[15] = "Spacing"
tagarr[16,0] = '0018'x
tagarr[16,1] = '0050'x ;Slice Thickness
descrarr[16] = "SliceTh"
tagarr[17,0] = '0028'x
tagarr[17,1] = '0030'x ;Pixel Spacing
descrarr[17] = "PixelSp"
tagarr[18,0] = '0028'x
tagarr[18,1] = '0010'x ;Rows
descrarr[18] = "Rows"
tagarr[19,0] = '0028'x
tagarr[19,1] = '0011'x ;Columns
descrarr[19] = "Cols"
tagarr[20,0] = '0018'x
tagarr[20,1] = '0095'x ;Pixel Bandwidth
descrarr[20] = "PixelBW [Hz/pixel]"
tagarr[21,0] = '0018'x
tagarr[21,1] = '0083'x ;Number of averages
descrarr[21] = "NEX"
tagarr[22,0] = '2001'x
tagarr[22,1] = '1003'x ;Diffusion B-Factor
descrarr[22] = "b-values [s/mm2]"
tagarr[23,0] = '2001'x
tagarr[23,1] = '1004'x ;Diffusion direction
descrarr[23] = "Diffusion Direction"

descrarr[TagsNumber]   = "BW"
descrarr[TagsNumber+1] = "FirstImage"
descrarr[TagsNumber+2] = "BriefDir"
descrarr[TagsNumber+3] = "#Images"
descrarr[TagsNumber+4] = "#Par"
descrarr[TagsNumber+5] = "#Rep"
descrarr[TagsNumber+6] = "Temporal Res[s]"
descrarr[TagsNumber+7] = "IsChecked?"
Notes='Notes: SE#=Series Number,Flip=Flip Angle, Spacing=Spacing Between Slices, SliceTh=Slice Thickness, PixelSp=Pixel Spacing, BW=Band Width (kHz), #Part=Number of Partitions, #Rep=Number of Repetitions'
tagarr_MR=tagarr
descrarr_MR=descrarr
Notes_MR=Notes

;define PET CT
TagsNumber = 24L
tagarr = lonarr(TagsNumber,2)
descrarr = strarr(TagsNumber + 8)
tagarr[0,0] = '0010'x
tagarr[0,1] = '0010'x ;Patient name
descrarr[0] = "PatientName"
tagarr[1,0] = '0010'x
tagarr[1,1] = '0020'x ;Patient ID
descrarr[1] = "PatientID"
tagarr[2,0] = '0010'x
tagarr[2,1] = '0030'x ;Patient date of birth
descrarr[2] = "PatientDOB"
tagarr[3,0] = '0010'x
tagarr[3,1] = '0040'x ;Patient sex
descrarr[3] = "PatientSex"
tagarr[4,0] = '0010'x
tagarr[4,1] = '1010'x ;Patient age
descrarr[4] = "PatientAge"
tagarr[5,0] = '0010'x
tagarr[5,1] = '1030'x ;Patient weight
descrarr[5] = "PatientWeight"
tagarr[6,0] = '0008'x
tagarr[6,1] = '0023'x ;Study date
descrarr[6] = "StudyDate"
tagarr[7,0] = '0008'x
tagarr[7,1] = '1030'x ;Study Description
descrarr[7] = "StudyDesc"
tagarr[8,0] = '0020'x
tagarr[8,1] = '0011'x ;Series number
descrarr[8] = "SE#"
tagarr[9,0] = '0008'x
tagarr[9,1] = '103e'x ;Series Description
descrarr[9] = "SeDesc"
tagarr[10,0] = '0008'x
tagarr[10,1] = '0031'x ;Series time
descrarr[10] = "Time"
tagarr[11,0] = '0020'x
tagarr[11,1] = '0010'x ;Study ID
descrarr[11] = "StudyID"
tagarr[12,0] = '0010'x
tagarr[12,1] = '1020'x ;Patient Height
descrarr[12] = "PatientHeight"
tagarr[13,0] = '0018'x
tagarr[13,1] = '1072'x ;Dose Start Time
descrarr[13] = "DoseStart"
tagarr[14,0] = '0018'x
tagarr[14,1] = '1074'x ;Dose
descrarr[14] = "Dose"
tagarr[15,0] = '0018'x
tagarr[15,1] = '0050'x ;Slice Thickness
descrarr[15] = "SliceTh"
tagarr[16,0] = '0028'x
tagarr[16,1] = '0030'x ;Pixel Spacing
descrarr[16] = "PixelSp"
tagarr[17,0] = '0028'x
tagarr[17,1] = '0010'x ;Rows
descrarr[17] = "Rows"
tagarr[18,0] = '0028'x
tagarr[18,1] = '0011'x ;Columns
descrarr[18] = "Cols"
tagarr[19,0] = '0018'x
tagarr[19,1] = '1075'x ;Half Life
descrarr[19] = "HalfLife"
tagarr[20,0] = '0054'x
tagarr[20,1] = '1001'x ;Units
descrarr[20] = "Units"
tagarr[21,0] = '0054'x
tagarr[21,1] = '1102'x ;Decay Correction
descrarr[21] = "DecayCorrection"
tagarr[22,0] = '0018'x
tagarr[22,1] = '0060'x ;kvP
descrarr[22] = "KVP[PeakKV]"
tagarr[23,0] = '0018'x
tagarr[21,1] = '1151'x ;x-Ray Tube Current
descrarr[23] = "x-rayTubeCurrent[mA]"

descrarr[TagsNumber]   = "Spacing"
descrarr[TagsNumber+1] = "FirstImage"
descrarr[TagsNumber+2] = "BriefDir"
descrarr[TagsNumber+3] = "#Images"
descrarr[TagsNumber+4] = "#Par"
descrarr[TagsNumber+5] = "#Rep"
descrarr[TagsNumber+6] = "Temporal Res[s]"
descrarr[TagsNumber+7] = "IsChecked?"
Notes='Notes: SE#=Series Number, DoseStart=Dose Start Time, Spacing=Spacing Between Slices, SliceTh=Slice Thickness, PixelSp=Pixel Spacing, #Part=Number of Partitions, #Rep=Number of Repetitions'
tagarr_PETCT=tagarr
descrarr_PETCT=descrarr
Notes_PETCT=Notes

END

pro browse_diffTE, afn, tagarr_MR, tagarr_PETCT, valuearr, problemReport, MRorPETCT
COMMON counter,counter

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

;afn is the .dcm filename
; Include IDL system procedures and functions when profiling:

;PROFILER, /SYSTEM


;;;DICOM Files Section
fnstar=STRTRIM(afn)+'*'
a=findfile(fnstar)
;print, 'a=(L12)', a

; Create the progress bar.
progressbar = Obj_New('progressbar', Color='Orange', Text='0%', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
; Place the progress bar on the display.
progressbar -> Start

;Go through each files under a directory
FOR i=0L,N_ELEMENTS(a)-1 DO BEGIN

        ;Call procedure in Klaus' code
Next:    fdecomp, a[i], disk, dir, name, qual, version
    finfo=file_info(a[i])

    IF STRLEN(name) EQ 0 THEN BEGIN

       ;efg is the last two digits of a a[i]
      efg= STRMID(a[i],STRLEN(a[i])-2,2)
      IF (efg NE ".\") AND (finfo.directory EQ 1) THEN browse_diffTE,a[i],tagarr_MR, tagarr_PETCT, valuearr, problemReport, MRorPETCT
         ; the a[i] is not the current directory .\ or ..\

         ;a[i] is a directory, do recursive calling

    ENDIF ELSE BEGIN

    count=(i+1.)/N_ELEMENTS(a)*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

     ;//Obtain the value of the tags in the DICOM image files
     ;;check the validity of dicom image
     ;result=QUERY_DICOM(a[i])
     result=1; here does not check the dicom validity
     first=STRMID(qual,0,1)
     IF (result NE 1) AND (BYTE(first) GE 48) AND (BYTE(first) LE 57) THEN BEGIN
        ;assume file always has 3-letter legal extension, if extension are digits, assume it is dcm file
        problemReport = [[problemReport], [a[i] + " could not be read as dcm (" + STRCOMPRESS(STRING(finfo.SIZE)) + " Bytes)"]]
     ENDIF ELSE BEGIN
        counter=counter+1

        obj = OBJ_NEW('IDLffDICOM', a[i]);, /NO_PIXEL_DATA)
        ;Assume dicom images must have these 6 parameters
        Temp_PatientName = STRING(*(obj->GetValue('0010'x,'0010'x,/NO_COPY))[0])
        Temp_PatientID = STRING(*(obj->GetValue('0010'x,'0020'x,/NO_COPY))[0])
        Temp_StudyDate = STRING(*(obj->GetValue('0008'x,'0023'x,/NO_COPY))[0])
        Temp_StudyDesc = STRING(*(obj->GetValue('0008'x,'1030'x,/NO_COPY))[0])
        Temp_SeriesNumber = STRING(*(obj->GetValue('0020'x,'0011'x,/NO_COPY))[0])
        Temp_SeriesDesc = STRING(*(obj->GetValue('0008'x,'103e'x,/NO_COPY))[0])
        Temp_SeriesTime = '';STRING(*(obj->GetValue('0008'x,'0031'x,/NO_COPY))[0])
;        IF (obj->QueryValue('0010,0010')) EQ 2 THEN Temp_PatientName = obj->GetValue('0010,0010') ELSE Temp_PatientName = '---'
;        IF (obj->QueryValue('0010,0020')) EQ 2 THEN Temp_PatientID = obj->GetValue('0010,0020') ELSE Temp_PatientID = '---'
;        IF (obj->QueryValue('0008,0020')) EQ 2 THEN Temp_StudyDate = obj->GetValue('0008,0020') ELSE Temp_StudyDate = '---'
;        IF (obj->QueryValue('0008,1030')) EQ 2 THEN Temp_StudyDesc = obj->GetValue('0008,1030') ELSE Temp_StudyDesc = '---'
;        IF (obj->QueryValue('0020,0011')) EQ 2 THEN Temp_SeriesNumber = obj->GetValue('0020,0011') ELSE Temp_SeriesNumber = '---'
;        IF (obj->QueryValue('0008,103e')) EQ 2 THEN Temp_SeriesDesc = obj->GetValue('0008,103e') ELSE Temp_SeriesDesc = '---'
;        IF (obj->QueryValue('0008,0031')) EQ 2 THEN Temp_SeriesTime = obj->GetValue('0008,0031') ELSE Temp_SeriesTime = '---'

            number_series=(SIZE(valuearr))[1]
            FOR k=0,number_series-1 DO BEGIN
                IF  (Temp_PatientName EQ valuearr[k,0]) AND $ ; Patient Name
                    (Temp_PatientID EQ valuearr[k,1]) AND $ ; Patient ID
                    (Temp_StudyDate EQ valuearr[k,6]) AND $ ; Study date
                    (Temp_StudyDesc EQ valuearr[k,7]) AND $ ; Study description
                    (Temp_SeriesNumber EQ valuearr[k,8]) AND $  ; Series number
                    (Temp_SeriesDesc EQ valuearr[k,9]) AND $  ; Series description
                    (Temp_SeriesTime EQ valuearr[k,10]) $  ; Series time
                 THEN GOTO, jump_dicom
            ENDFOR

jump_dicom: IF k GE number_series THEN BEGIN
                TR = obj->GetValue('0018'x,'0080'x,/NO_COPY);(SIZE(TR,/N_DIMENSIONS) NE 0)(obj->QueryValue('0018,0080')) NE 2
                Modality=*(obj->GetValue('0008'x,'0060'x,/NO_COPY))[0]
                IF (STRCMP(Modality, 'MR', 2, /FOLD_CASE) EQ 0) THEN BEGIN ;check modality, GJ 2020/04/13
                    MRorPETCT=0
                    tagarr=tagarr_PETCT
                    TagsNumber=(SIZE(tagarr))[1]
                    ;Build an array of values of tags, the last element of which is the file name
                    Temp_value = strarr(1,TagsNumber + 9)
                    for j = 0, TagsNumber-1 do begin
                      temp = obj->GetValue(tagarr[0,0],tagarr[0,1],/NO_COPY)
                      IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
                        value = STRING(*(obj->GetValue(tagarr[0,0],tagarr[0,1],/NO_COPY))[0]);obj->GetValue(tagarr[0],/NO_COPY)
                      ENDIF ELSE BEGIN
                        value= '---'
                      ENDELSE
                      temp = obj->GetValue(tagarr[j,0],tagarr[j,1],/NO_COPY)
                      IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
                        Temp_value[0,j] = STRJOIN(STRING(*(obj->GetValue(tagarr[j,0],tagarr[j,1],/NO_COPY))[0]),'\')
                      ENDIF ELSE BEGIN
                        Temp_value[0,j]= '---'
                      ENDELSE
;                      ;IF (obj->QueryValue(tagarr[0])) EQ 2 THEN value = obj->GetValue(tagarr[0]) ELSE value= '---'
;                        value = STRING(*(obj->GetValue(tagarr[0,0],tagarr[0,1],/NO_COPY))[0]);obj->GetValue(tagarr[0],/NO_COPY)
;                        Temp_value[0,j] = STRJOIN(STRING(*(obj->GetValue(tagarr[j,0],tagarr[j,1],/NO_COPY))[0]),'\')
;;                        IF (obj->QueryValue(tagarr[j])) EQ 2 THEN Temp_value[0,j] = STRJOIN(STRING(obj->GetValue(tagarr[j])),'\') ELSE Temp_value[0,j]= '---'
                    endfor

                    ;Do shorten the directories, only select the last 5 digits in the lower 2 levels of directory
                    afn_split=strsplit(afn, '\', COUNT=afn_split_num, /EXTRACT)
                    afn_brief=strmid(afn_split(afn_split_num-2), 4, 5, /reverse_offset) + "\" + strmid(afn_split(afn_split_num-1), 4, 5, /reverse_offset)
                    Temp_value[0,TagsNumber+1] = a[i]
                    Temp_value[0,TagsNumber+2] = afn_brief
                    Temp_value[0,TagsNumber+3] = 1
                    Temp_value[0,TagsNumber+7] = '(  )'

                    ;Add this seires
                    ;IF k GT 1 THEN BEGIN
                    Valuearr=[Valuearr, Temp_value]
                    ;ENDIF ELSE BEGIN
                    ;  Valuearr=Temp_value
                    ;ENDELSE

                    dummy=browse_directory(a[i],fn_files,status)
                    status=sort_DCM_by_ser_and_ima_number(fn_files, fn_files, sl_slab, tError)
                    IF status NE 1 THEN BEGIN
                        slice_TimePT=N_ELEMENTS(fn_files)/sl_slab
                        Valuearr[k, TagsNumber+4]=sl_slab; number of slices
                        ima1_obj=OBJ_NEW('IDLffDICOM', fn_files[0]);, /NO_PIXEL_DATA)
                        ima2_obj=OBJ_NEW('IDLffDICOM', fn_files[1]);, /NO_PIXEL_DATA)
                        sliceLocation1 = STRING(*(ima1_obj->GetValue('0020'x,'1041'x,/NO_COPY))[0])
                        sliceLocation2 = STRING(*(ima2_obj->GetValue('0020'x,'1041'x,/NO_COPY))[0])
;                        IF (ima1_obj->QueryValue('0020,1041')) EQ 2 THEN sliceLocation1 = ima1_obj->GetValue('0020,1041') ELSE sliceLocation1 = 0.
;                        IF (ima2_obj->QueryValue('0020,1041')) EQ 2 THEN sliceLocation2 = ima2_obj->GetValue('0020,1041') ELSE sliceLocation2 = 0.
                        Spacing=ABS(DOUBLE(sliceLocation1)-DOUBLE(sliceLocation2))
                        Valuearr[k, TagsNumber]=Spacing
                        OBJ_DESTROY, ima1_obj
                        OBJ_DESTROY, ima2_obj
                    ENDIF ELSE BEGIN
                        Valuearr[k, TagsNumber+4]='---'; number of slices
                        Valuearr[k, TagsNumber]='---'
                    ENDELSE
                ENDIF ELSE BEGIN
                    MRorPETCT=1
                    tagarr=tagarr_MR
                    TagsNumber=(SIZE(tagarr))[1]
                    ;Build an array of values of tags, the last element of which is the file name
                    Temp_value = strarr(1,TagsNumber + 9)
                    for j = 0, TagsNumber-1 do begin
                      temp = obj->GetValue(tagarr[0,0],tagarr[0,1],/NO_COPY)
                      IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
                        value = STRING(*(obj->GetValue(tagarr[0,0],tagarr[0,1],/NO_COPY))[0]);obj->GetValue(tagarr[0],/NO_COPY)
                      ENDIF ELSE BEGIN
                        value= '---'
                      ENDELSE
                      temp = obj->GetValue(tagarr[j,0],tagarr[j,1],/NO_COPY)
                      IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
                        Temp_value[0,j] = STRJOIN(STRING(*(obj->GetValue(tagarr[j,0],tagarr[j,1],/NO_COPY))[0]),'\')
                      ENDIF ELSE BEGIN
                        Temp_value[0,j]= '---'
                      ENDELSE
                      
;                      print, j, '   ', Temp_value[0,j]
;                      print, 'test'
                      ;IF (obj->QueryValue(tagarr[0])) EQ 2 THEN value = obj->GetValue(tagarr[0]) ELSE value= '---'
                      ;  IF (obj->QueryValue(tagarr[j])) EQ 2 THEN Temp_value[0,j] = STRJOIN(STRING(obj->GetValue(tagarr[j])),'\') ELSE Temp_value[0,j]= '---'
                    endfor

                    ;Bandwidth
                    Image_Bandwidth='unknown'
;                    IF (obj->QueryValue('0018,0095')) EQ 2 THEN PixelBandwidth = obj->GetValue('0018,0095') ELSE PixelBandwidth = 0.
                    temp = obj->GetValue('0018'x,'0095'x,/NO_COPY)
                    IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
                      PixelBandwidth = DOUBLE(*(obj->GetValue('0018'x,'0095'x,/NO_COPY))[0])
                    ENDIF ELSE BEGIN
                      PixelBandwidth = 0.
                    ENDELSE
                    IF (PixelBandwidth NE 0.) THEN BEGIN
                        ;IF (obj->QueryValue('0018,1312')) EQ 2 THEN PhaseEncoDirection = obj->GetValue('0018,1312') ELSE PhaseEncoDirection = '---'
                        temp = obj->GetValue('0018'x,'1312'x,/NO_COPY)
                        IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
                            PhaseEncoDirection = STRING(*(obj->GetValue('0018'x,'1312'x,/NO_COPY))[0])
                            IF STRCMP(PhaseEncoDirection, 'ROW', 3, /FOLD_CASE) THEN BEGIN
                                Temp_value[0,TagsNumber]= (PixelBandwidth) * FLOAT(Temp_value[0,18]) / 1000.0
                            ENDIF ELSE BEGIN
                                IF STRCMP(PhaseEncoDirection, 'COL', 3, /FOLD_CASE) THEN Temp_value[0,TagsNumber]= (PixelBandwidth) *FLOAT(Temp_value[0,17])/1000.0
                            ENDELSE
                        ENDIF
                    ENDIF

                    ;Do shorten the directories, only select the last 5 digits in the lower 2 levels of directory
                    afn_split=strsplit(afn, '\', COUNT=afn_split_num, /EXTRACT)
                    afn_brief=strmid(afn_split(afn_split_num-2), 4, 5, /reverse_offset) + "\" + strmid(afn_split(afn_split_num-1), 4, 5, /reverse_offset)
                    Temp_value[0,TagsNumber+1] = a[i]
                    Temp_value[0,TagsNumber+2] = afn_brief
                    Temp_value[0,TagsNumber+3] = 1
                    Temp_value[0,TagsNumber+7] = '(  )'

                    ;Add this seires
                    ;IF k GT 1 THEN BEGIN
                    Valuearr=[Valuearr, Temp_value]
                    ;ENDIF ELSE BEGIN
                    ;  Valuearr=Temp_value
                    ;ENDELSE

                    dummy=browse_directory(a[i],fn_files,status)
                    status=sort_DCM_by_ser_and_ima_number(fn_files, fn_files, sl_slab, tError)
                    Valuearr[k, TagsNumber+4]=sl_slab; number of slices
                    numTimept=N_ELEMENTS(fn_files)/sl_slab
                    Valuearr[k, TagsNumber+5]=numTimept
                    IF numTimept GT 1 THEN BEGIN
                        ima1_obj=OBJ_NEW('IDLffDICOM', fn_files[0]);, /NO_PIXEL_DATA)
                        ima2_obj=OBJ_NEW('IDLffDICOM', fn_files[sl_slab]);, /NO_PIXEL_DATA)
                        ;IF (ima1_obj->QueryValue('0008,0070')) EQ 2 THEN manufacturer = ima1_obj->GetValue('0008,0070') ELSE manufacturer = '---'
                        manufacturer = STRING(*(ima1_obj->GetValue('0008'x,'0070'x,/NO_COPY))[0])
                        IF STRCMP(manufacturer, "GE", 2, /FOLD_CASE) THEN BEGIN
                            ;IF (ima1_obj->QueryValue('0018,1060')) EQ 2 THEN trigger_time1 = ima1_obj->GetValue('0018,1060') ELSE trigger_time1 = 0.
                            ;IF (ima2_obj->QueryValue('0018,1060')) EQ 2 THEN trigger_time2 = ima2_obj->GetValue('0018,1060') ELSE trigger_time2 = 0.
                            trigger_time1 = STRING(*(ima1_obj->GetValue('0018'x,'1060'x,/NO_COPY))[0])
                            trigger_time2 = STRING(*(ima2_obj->GetValue('0018'x,'1060'x,/NO_COPY))[0])
                            temporalRes=(DOUBLE(trigger_time2)-DOUBLE(trigger_time1))/1000. ; ms->s
                            Valuearr[k, TagsNumber+6]=temporalRes
                        ENDIF ELSE BEGIN
                            IF STRCMP(manufacturer, "SIEMENS", 2, /FOLD_CASE) THEN BEGIN
                                ;IF (ima1_obj->QueryValue('0008,0033')) EQ 2 THEN image_time1 = ima1_obj->GetValue('0008,0033') ELSE image_time1 = 0.
                                ;IF (ima2_obj->QueryValue('0008,0033')) EQ 2 THEN image_time2 = ima2_obj->GetValue('0008,0033') ELSE image_time2 = 0.
                                image_time1 = STRING(*(ima1_obj->GetValue('0008'x,'0033'x,/NO_COPY))[0])
                                image_time2 = STRING(*(ima2_obj->GetValue('0008'x,'0033'x,/NO_COPY))[0])
                                temporalRes=(DOUBLE(image_time2)-DOUBLE(image_time1))
                                Valuearr[k, TagsNumber+6]=temporalRes
                            ENDIF
                        ENDELSE
                        OBJ_DESTROY, ima1_obj
                        OBJ_DESTROY, ima2_obj

                        ;find out different TE
                        ;IF STRCMP(manufacturer, "GE", 2, /FOLD_CASE) THEN BEGIN
                        TEarr1=STRARR(numTimept)
                        TEarray='' ;to save different TEs
                        FOR j=0L, numTimept-1L DO BEGIN
                            ima1_obj=OBJ_NEW('IDLffDICOM', fn_files[j*sl_slab]);, /NO_PIXEL_DATA)
                            echo_time1 =  STRING(*(ima1_obj->GetValue('0018'x,'0081'x,/NO_COPY))[0])
;                           IF (ima1_obj->QueryValue('0018,0081')) EQ 2 THEN echo_time1 = ima1_obj->GetValue('0018,0081') ELSE echo_time1 = 0.
                            TEarr1[j]=STRING(echo_time1, FORMAT='(F6.2)')
                            IF j EQ 0 THEN BEGIN
                                TEarray=TEarray+STRING(echo_time1, FORMAT='(F6.2)')
                            ENDIF ELSE BEGIN
                                IF STRCMP(TEarr1[j], TEarr1[j-1],  /FOLD_CASE) NE 1 THEN TEarray=TEarray+','+STRING(echo_time1, FORMAT='(F6.2)')
                            ENDELSE
                            OBJ_DESTROY, ima1_obj
                        ENDFOR
                        Valuearr[k,13]=STRCOMPRESS(TEarray,/REMOVE_ALL) ; put multiple TEs to TE array

                     ;find out different b-values
                        ;IF STRCMP(manufacturer, "GE", 2, /FOLD_CASE) THEN BEGIN
                        b_valuearr1=STRARR(numTimept)
                        b_valuearray='' ;to save different TEs
                        FOR j=0L, numTimept-1L DO BEGIN
                            ima1_obj=OBJ_NEW('IDLffDICOM', fn_files[j*sl_slab]);, /NO_PIXEL_DATA)
                            b_value1 = STRING(*(ima1_obj->GetValue('2001'x,'1003'x,/NO_COPY))[0])
    ;                        IF (ima1_obj->QueryValue('2001,1003')) EQ 2 THEN b_value1 = ima1_obj->GetValue('2001,1003') ELSE b_value1 = 0.
                            b_valuearr1[j]=STRING(b_value1, FORMAT='(F6.1)')
                            IF j EQ 0 THEN BEGIN
                                b_valuearray=b_valuearray+STRING(b_value1, FORMAT='(F6.1)')
                            ENDIF ELSE BEGIN
                                IF STRCMP(b_valuearr1[j], b_valuearr1[j-1],  /FOLD_CASE) NE 1 THEN b_valuearray=b_valuearray+','+STRING(b_value1, FORMAT='(F6.1)')
                            ENDELSE
                            OBJ_DESTROY, ima1_obj
                        ENDFOR
                        Valuearr[k,22]=STRCOMPRESS(b_valuearray,/REMOVE_ALL) ; put multiple b_values to b_value array

                        ;find out different flip angle
                        ;IF STRCMP(manufacturer, "GE", 2, /FOLD_CASE) THEN BEGIN
                        FAarr1=STRARR(numTimept)
                        FAarray='' ;to save different FAs
                        FOR j=0L, numTimept-1L DO BEGIN
                            ima1_obj=OBJ_NEW('IDLffDICOM', fn_files[j*sl_slab]);, /NO_PIXEL_DATA)
                            flip_angle1 = STRING(*(ima1_obj->GetValue('0018'x,'1314'x,/NO_COPY))[0])
                            ;IF (ima1_obj->QueryValue('0018,1314')) EQ 2 THEN flip_angle1 = ima1_obj->GetValue('0018,1314') ELSE flip_angle1 = 0.
                            FAarr1[j]=STRING(flip_angle1, FORMAT='(I6.0)')
                            IF j EQ 0 THEN BEGIN
                                FAarray=FAarray+STRING(flip_angle1, FORMAT='(I6.0)')
                            ENDIF ELSE BEGIN
                                IF STRCMP(FAarr1[j], FAarr1[j-1],  /FOLD_CASE) NE 1 THEN FAarray=FAarray+','+STRING(flip_angle1, FORMAT='(I6.0)')
                            ENDELSE
                            OBJ_DESTROY, ima1_obj
                        ENDFOR
                        Valuearr[k,14]=STRCOMPRESS(FAarray,/REMOVE_ALL) ; put multiple FAs to FA array

                    ENDIF
                ENDELSE
             ENDIF ELSE BEGIN
                Valuearr[k, TagsNumber+3]=Valuearr[k, TagsNumber+3]+1
                dir_arr=STRSPLIT(Valuearr[k, TagsNumber+1], '*', COUNT=numDir, /EXTRACT)
                flag=0
                FOR countNumDir=0, numDir-1 DO BEGIN
                    fdecomp, dir_arr[countNumDir], disk_k, dir_k, name_k, qual_k, version_k
                    IF ((STRCMP(disk_k, disk, /FOLD_CASE) EQ 1) AND (STRCMP(dir_k, dir, /FOLD_CASE) EQ 1)) THEN flag=1
                ENDFOR
                IF (flag EQ 0) THEN BEGIN
                    Valuearr[k, TagsNumber+1]=Valuearr[k, TagsNumber+1]+"*" + a[i]
                    a_split=strsplit(disk+dir, '\', COUNT=a_split_num, /EXTRACT)
                    a_brief=strmid(a_split(a_split_num-2), 4, 5, /reverse_offset) + "\" + strmid(a_split(a_split_num-1), 4, 5, /reverse_offset)
                    Valuearr[k, TagsNumber+2]=Valuearr[k, TagsNumber+2]+"*" + a_brief
                ENDIF
             ENDELSE

             ;Destruction of the object
             OBJ_DESTROY, obj
       ENDELSE
     ENDELSE
ENDFOR

progressbar -> Destroy
; Retrieve the profiling results:
;PROFILER, /REPORT

END


;+
; NAME:
;   DICOM_Crawler
; PURPOSE:
;   Browse dcm file directories and report important scanning parameters
; EXPLANATION:
;     IDL V5.5 or later!
;
; CALLING SEQUENCE:
;   DICOM_Crawler
;
; INPUT:
;
;
; INPUT-OUTPUT:
;
;
;
; OUTPUT:
;
;
; EXAMPLE AND USAGE:
;   1. run this program.
;   2. a popup window is for you to select the directory where your dcm images are saved. Click "OK" after selection.
;   3. the program will go through every dcm file (with .dcm or without suffix) and read out the scanning parameters.
;   4. finally, a popup window is for you to save the browser results into a .txt or .xls file
;
;   ---------------------------------------------
;   Use Microsoft Excel to organize your browse result
;   5. open MS Excel
;   6. in Excel, "Data" --> "Import External Data" --> "Import Data ...", a popup window can allow you select the txt file
;      you created in step 7.
;   7. a "Text Import Wizard -Step 1 of 3" window jumps out. You can directly select "Finish".
;   8. now you can do anything in Excel to the browse results of dcm images and save them as .xls file.
;
;
; HISTORY
;   9th version, ��GJ(jia.11@osu.edu), Sep, 09, 2003
;               include the file names without dcm suffix.
;               expand the situation that several directories have the same sereis .
;   10th version,��GJ(jia.11@osu.edu), Sep, 10, 2003
;               add the error handler of reading dcm files.
;   11th version,��GJ(jia.11@osu.edu), Dec, 22, 2003
;               add patient DOB, age, sex, weight in case of mixture of patients.
;   12th version,��GJ(jia.11@osu.edu), Jan, 08, 2004
;               1. add the problem report at the end of the file for some images have zero bytes.
;               2. optimize the browse program, remove the limit that filename must be **.dcm.
;   13th version,��GJ(jia.11@osu.edu), Feb, 05, 2004
;               find out temporal sampling data, i.e. temporal resolution from DICOM tag field 0018,1060 ("Trigger Time")
;                   from GE scanner if possible
;   14th version,��GJ(jia.11@osu.edu), March, 24, 2004
;               let the program check the different TE and save them
;   15th version,��GJ(jia.11@osu.edu), March, 26, 2004
;               a minor changes about "browse_diffTE", different TE can be read in both GE and Siemens images
;   16th version,��GJ(jia.11@osu.edu), Jan, 11, 2005
;               Include codes to take out specific tags from PET and CT images
PRO Crawler
COMMON counter,counter

Tags_Define, tagarr_MR, descrarr_MR, Notes_MR, tagarr_PETCT, descrarr_PETCT, Notes_PETCT
;descrarr[TagsNumber+4] = "Bandwidth [kHz]"
problemReport = "!!!Problem Report: "

;Print the tags in the scree
print, 'Tags will be searched:'
;print, format = '(%"%4.4Z\t%4.4Z")', tagarr

counter=0
;Initialize value array

valuearr=STRARR(1, N_ELEMENTS(descrarr_MR)+1)
valuearr[0,8]=0 ; Initialize SeriesNumber
valuearr[0,25]=0 ; Initialize Number of Images
afn= DIALOG_PICKFILE(/MUST_EXIST,/DIRECTORY, title="Select Image Folder")
IF STRLEN(afn) NE 0 THEN BEGIN

    ;IF MRorPETCT=1, MR; if MRorPETCT=0, PETCT
    browse_diffTE,afn, tagarr_MR, tagarr_PETCT, valuearr, problemReport, MRorPETCT

    IF MRorPETCT EQ 1 THEN BEGIN
        descrarr=descrarr_MR
        Notes=Notes_MR
    ENDIF ELSE BEGIN
        descrarr=descrarr_PETCT
        Notes=Notes_PETCT
    ENDELSE
    ;Select a file to save exported data, the default opened directory is the one you select just now
    aSaveFile=DIALOG_PICKFILE(FILTER = ['*.txt', '*.xls', '*.xlsx'], title='Save Your Output as EXCEL file', PATH=afn)
    ;// open an ASCII file for output
    ;IF STRCMP(STRMID(aSaveFile, 2, 3, /REVERSE_OFFSET), 'txt',/FOLD_CASE)   THEN BEGIN
    IF (STRLEN(aSaveFile) NE 0)  THEN BEGIN
       GET_LUN, U
       openw,U,aSaveFile
       ;printf, 5, 'Clinical Protocal:'
       ;Save the descriptions of tags as the first row in the output ASCII file
       ;Mix "A" format code and "\" character escapes
       ;(see Building IDL applications, pp. 190)
       printf, U, FORMAT = '(50(A, %"\t"))', descrarr
       ;;Save this value array to output .txt file
       Num_series=(SIZE(valuearr))[1]
       valuearr=TRANSPOSE(valuearr)
       FOR i=1L, Num_series-1 DO BEGIN
         printf, U, FORMAT = '(50(A, %"\t"))', valuearr[*, i]
       ENDFOR
       ;print, 'valuearr=',valuearr

;;Print the directories for volume analysis, by filtering the series description
       ;PRINTF, U, '+++++++++++++++++++++++++++++'
       ;PRINTF, U, 'T1w and T2w for volume analysis'
       ;FOR i=1L, Num_series-1 DO BEGIN
       ; IF (STRPOS(valuearr[4, i], 'Ax') NE -1) AND (STRPOS(valuearr[4, i], '2D') NE -1) AND (STRPOS(valuearr[4, i], 'SE') NE -1) THEN BEGIN
       ;     printf, U, FORMAT = '(50(A, %"\t"))', valuearr[*, i]
       ; ENDIF
       ;ENDFOR

;;Print the directories for diameter analysis, by filtering the sereis description
       ;PRINTF, U, '+++++++++++++++++++++++++++++'
       ;PRINTF, U, 'Axial and Coronal for diameter analysis'
       ;FOR i=1L, Num_series-1 DO BEGIN
    ;   IF (STRPOS(valuearr[4, i], '2D') NE -1) AND (STRPOS(valuearr[4, i], 'SE') NE -1) AND (STRPOS(valuearr[4, i], 'Fast') NE -1) THEN BEGIN
;          printf, U, FORMAT = '(50(A, %"\t"))', valuearr[*, i]
;      ENDIF
;     ENDFOR

       printf, U, Notes
       printf, U, ""
       printf, U, problemReport
       FREE_LUN, U
    ENDIF
ENDIF
END