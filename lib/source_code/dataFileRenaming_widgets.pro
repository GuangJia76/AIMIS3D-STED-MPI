;+
;  FILE:
;       dataFileRenaming_widgets.pro
;
;  CALLING SEQUENCE: dataFileRenaming_widgets
;
;  PURPOSE:
;       This datefile renaming.
;
;  MAJOR TOPICS: Widgets
;
;  CATEGORY:
;       IDL datefile renaming System
;
;  INTERNAL FUNCTIONS and PROCEDURES:
;       pro dataFileRenaming_widgetsEvent            -  Event handler
;       pro dataFileRenaming_widgetsCleanup          -  Cleanup
;       pro dataFileRenaming_widgets                 -  Main procedure
;
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       widgets.tip
;       pro datefile renaming_gettips            - Read the tip file and create widgets
;
;  REFERENCE: IDL Reference Guide, IDL User's Guide
;
;  NAMED STRUCTURES:
;       none.
;
;  COMMON BLOCS:
;       none.
;
;  MODIFICATION HISTORY:
;   2022/7/6, @GJ
;-

; -----------------------------------------------------------------------------
;
;  Purpose:  Function returns the 3 angles of a space three 1-2-3
;            given a 3 x 3 cosine direction matrix
;            else -1 on failure.
;
;  Definition :  Given 2 sets of dextral orthogonal unit vectors
;                (a1, a2, a3) and (b1, b2, b3), the cosine direction matrix
;                C (3 x 3) is defined as the dot product of:
;
;                C(i,j) = ai . bi  where i = 1,2,3
;
;                A column vector X (3 x 1) becomes X' (3 x 1)
;                after the rotation as defined as :
;
;                X' = C X
;
;                The space three 1-2-3 means that the x rotation is first,
;                followed by the y rotation, then the z.
;
function angle3123, $
    cosMat           ; IN: cosine direction matrix (3 x 3)

    ;  Verify the input parameters
    ;
    if (N_PARAMS() ne 1) then begin
        PRINT,'Error in angle3123: 1 parameters must be passed.'
        RETURN, -1
    endif
    sizec = size(cosMat)
    if (sizec[0] ne 2) then begin
        PRINT,'Error, the input matrix must be of dimension 2'
        RETURN, -1
    endif
    if ((sizec[1] ne 3) or (sizec[2] ne 3)) then begin
        PRINT,'Error, the input matrix must be 3 by 3'
        RETURN, -1
    endif

    ;  Compute the 3 angles (in degrees)
    ;
    cosMat = TRANSPOSE(cosMat)
    angle = FLTARR(3)
    angle[1] = -cosMat[2,0]
    angle[1] = ASIN(angle[1])
    c2 = COS(angle[1])
    if (ABS(c2) lt 1.0e-6) then begin
        angle[0] = ATAN(-cosMat[1,2], cosMat[1,1])
        angle[2] = 0.0
    endif else begin
        angle[0] = ATAN( cosMat[2,1], cosMat[2,2])
        angle[2] = ATAN( cosMat[1,0], cosMat[0,0])
    endelse
    angle = angle * (180.0/!DPI)

    RETURN, angle

end    ;   of angle3123

; -----------------------------------------------------------------------------
;
;  Purpose:  Function returns the cosine direction matrix (3 x 3)
;            given the space three 1-2-3 rotation angles(i.e. rotation around
;            x axis, followed by Y axis, then z axis),
;            else -1 on failure.
;
;  Definition :  Given 2 sets of dextral orthogonal unit vectors
;                (a1, a2, a3) and (b1, b2, b3), the cosine direction matrix
;                C (3 x 3) is defined as the dot product of:
;
;                C(i,j) = ai . bi  where i = 1,2,3
;
;                A column vector X (3 x 1) becomes X' (3 x 1)
;                after the rotation as defined as :
;
;                X' = C X
;
function space3123, $
    theta, $        ; IN: angle of rotation around the x axis(in degrees)
    phi, $          ; IN: angle of rotation around the y axis(in degrees)
    gamma           ; IN: angle of rotation around the z axis(in degrees)

    ;  Verify the input parameters.
    ;
    if (N_PARAMS() ne 3) then begin
        PRINT,'Error in space3123: 3 parameters must be passed.'
        RETURN, -1
    endif

    cosMat = FLTARR(3, 3)

    ;  Transform the angle in radians.
    ;
    rTheta = theta * !DPI / 180.0
    rPhi = Phi * !DPI / 180.0
    rGamma = Gamma * !DPI / 180.0

    cos1 = COS(rTheta)
    cos2 = COS(rPhi)
    cos3 = COS(rGamma)
    sin1 = SIN(rTheta)
    sin2 = SIN(rPhi)
    sin3 = SIN(rGamma)

    ;  Compute the cosine direction matrix.
    ;
    cosMat[0,0] = cos2*cos3
    cosMat[1,0] = cos2*sin3
    cosMat[2,0] = -sin2
    cosMat[0,1] = (sin1*sin2*cos3) - (cos1*sin3)
    cosMat[1,1] = (sin1*sin2*sin3) + (cos1*cos3)
    cosMat[2,1] = sin1*cos2
    cosMat[0,2] = (cos1*sin2*cos3) + (sin1*sin3)
    cosMat[1,2] = (cos1*sin2*sin3) - (sin1*cos3)
    cosMat[2,2] = cos1*cos2

    RETURN, cosMat

end    ;   of space3123

; -----------------------------------------------------------------------------
;
;  Purpose:  Draw view.  Some platforms throw math errors that are
;            beyond our control.  Supress the printing of those errors.
;
pro dataFileRenaming_widgetsDraw, sState

;  Flush and print any accumulated math errors
;
void = check_math(/print)

;  Silently accumulate any subsequent math errors, unless we are debuggung.
;
orig_except = !except
!except = ([0, 2])[keyword_set(sState.debug)]

;  Draw.
;
;sState.drawWindowID->Draw, sState.oView

;  Silently (unless we are debuggung) flush any accumulated math errors.
;
void = check_math(PRINT=keyword_set(sState.debug))

;  Restore original math error behavior.
;
!except = orig_except
end

;@GJ, 2022/6/5, define the nomanclature
PRO data_naming_nomanclature, sState
  
  ;check month and day for single digit
  IF sState.mm LT 10 THEN BEGIN
    IF sState.dd LT 10 THEN BEGIN
      data_name_temp = STRING(sState.yy)+'0'+STRING(sState.mm)+'0'+STRING(sState.dd)
    ENDIF ELSE BEGIN
      data_name_temp = STRING(sState.yy)+'0'+STRING(sState.mm)+STRING(sState.dd)
    ENDELSE
  ENDIF ELSE BEGIN
    IF sState.dd LT 10 THEN BEGIN
      data_name_temp = STRING(sState.yy)+STRING(sState.mm)+'0'+STRING(sState.dd)
    ENDIF ELSE BEGIN
      data_name_temp = STRING(sState.yy)+STRING(sState.mm)+STRING(sState.dd)
    ENDELSE
  ENDELSE
  
  IF sState.ParticleConShow EQ 0 THEN BEGIN
    data_name_temp = data_name_temp+'_'+sState.ParticleName
  ENDIF ELSE BEGIN
    data_name_temp = data_name_temp+'_'+sState.ParticleName+'c'+STRING(sState.ParticleConStr)
  ENDELSE
  
  IF STRCMP(sState.BufferName, 'None', 4, /FOLD_CASE) NE 1 THEN BEGIN
    data_name_temp = data_name_temp+'_'+sState.BufferName+'c'+STRING(sState.BufferConStr)
  ENDIF
  
  IF STRCMP(sState.Shape, 'None', 4, /FOLD_CASE) NE 1 THEN BEGIN
    IF sState.FlatPortionShow GT 0 THEN BEGIN
      data_name_temp = data_name_temp+'_'+sState.Shape+'f'+STRING(sState.FrequencyStr)+'ka'+STRING(sState.AmpStr)+'fp'+STRING(sState.FlatPortionStr)+'%'
    ENDIF ELSE BEGIN
      data_name_temp = data_name_temp+'_'+sState.Shape+'f'+STRING(sState.FrequencyStr)+'ka'+STRING(sState.AmpStr)
    ENDELSE
  ENDIF
  
  IF sState.TempShow EQ 1 THEN BEGIN
    data_name_temp = data_name_temp+'_t'+STRING(sState.TempStr)+'dg'
  ENDIF
  
  IF sState.Rep GT 0 THEN BEGIN
    data_name_temp = data_name_temp+'_r'+STRING(sState.Rep)
  ENDIF
  
  data_name = STRCOMPRESS(data_name_temp, /REMOVE_ALL)
  sState.data_name = data_name
  WIDGET_CONTROL, sState.wTextField, SET_VALUE=STRING(sState.data_name)
end

;
; -----------------------------------------------------------------------------
;
;  Purpose:  Event handler
;
pro dataFileRenaming_widgetsEvent, $
    sEvent     ; IN: event structure

    ;  Quit the application using the close box.
    ;
    if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ $
        'WIDGET_KILL_REQUEST') then begin
        WIDGET_CONTROL, sEvent.top, /DESTROY
        RETURN
    endif

    ;  Get the info structure from top-level base.
    ;
    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

    ;  Determine which event.
    ;
    WIDGET_CONTROL, sEvent.id, GET_UVALUE=eventval

    ;  Take the following action based on the corresponding event.
    ;
    case eventval of


        'TABLE' : begin

            ;  Insertion of character.
            ;
            if (sEvent.type eq 0) then begin

                ;  Take action only when the carriage return character
                ;  is typed.
                ;
                if (sEvent.ch EQ 13) then begin
                    WIDGET_CONTROL, sState.wDataTable, GET_VALUE=data

                    sState.oTextureImage->SetProperty, DATA=bytscl(data)

                    ;  Draw the surface.
                    ;
                    if(sState.textureFlag eq 1) then begin
                        sState.oSimpleSurface->SetProperty, $
                            TEXTURE_MAP=sState.oTextureImage
                    endif
                    sState.oSimpleSurface->SetProperty, DATAZ=data
                    dataFileRenaming_widgetsDraw, sState

                endif
            endif
        end        ;   of TABLE

       'POPINFO' : begin
           widget_control, sState.wTextField,  get_value=textVal
           void = DIALOG_MESSAGE(['This is a message displayed by the ', $
                                 'DIALOG_MESSAGE routine.  It can be used ', $
                                 'to display messages that a user must ',$
                                 'acknowledge before the program continues.', $
                                 '', $
                                 'For example, the text from the Editable', $
                                 'Text field is:', $
                                 textVal])

        end           ;  of POPINFO
       
        'YearSLIDER': begin
          WIDGET_CONTROL, sState.wDateyearSlider, GET_VALUE=yy_temp
          ;update the year
          sState.yy = yy_temp
          print, 'yy=',sState.yy
          data_naming_nomanclature, sState
        end
        
        'MonthSLIDER': begin
          WIDGET_CONTROL, sState.wDatemonthSlider, GET_VALUE=mm_temp
          ;update the month
          sState.mm = mm_temp
          print, 'mm=',sState.mm
          data_naming_nomanclature, sState
        end
        
        'DaySLIDER': begin
          WIDGET_CONTROL, sState.wDatedaySlider, GET_VALUE=dd_temp
          ;update the day
          sState.dd = dd_temp
          print, 'dd=',sState.dd
          data_naming_nomanclature, sState
        end
        
        'ParticleList': begin
          WIDGET_CONTROL, /HOURGLASS
          listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)
          case listValue of
            ;  Synomag particle.
            ;
            0 : begin
              sState.ParticleName = sState.ParticleNameList[0]
              data_naming_nomanclature, sState
            end    ;    of 0
            
            ;  Perimag particle.
            ;
            1 : begin
              sState.ParticleName = sState.ParticleNameList[1]
              data_naming_nomanclature, sState
            end    ;    of 1
            
            ;  NE10 particle.
            ;
            2 : begin
              sState.ParticleName = sState.ParticleNameList[2]
              data_naming_nomanclature, sState
            end    ;    of 2
            
            ;  NE25 particle.
            ;
            3 : begin
              sState.ParticleName = sState.ParticleNameList[3]
              data_naming_nomanclature, sState
            end    ;    of 3

            ;  NE50 particle.
            ;
            4 : begin
              sState.ParticleName = sState.ParticleNameList[4]
              data_naming_nomanclature, sState
            end    ;    of 4
            
            ;  Resovist particle.
            ;
            5 : begin
              sState.ParticleName = sState.ParticleNameList[5]
              data_naming_nomanclature, sState
            end    ;    of 5
            
            ;  Vivotrax particle.
            ;
            6 : begin
              sState.ParticleName = sState.ParticleNameList[6]
              data_naming_nomanclature, sState
            end    ;    of 6
            
            ;  Dox-CACN particle from Chem, Song
            ;
            7 : begin
              sState.ParticleName = sState.ParticleNameList[7]
              data_naming_nomanclature, sState
            end    ;    of 7  

            ;  FeCo@C-PEG particle from NBE, Song
            ;
            8 : begin
              sState.ParticleName = sState.ParticleNameList[8]
              data_naming_nomanclature, sState
            end    ;    of 8     
            
            ;  IO particle from Nano Letters
            ;
            9 : begin
              sState.ParticleName = sState.ParticleNameList[9]
              data_naming_nomanclature, sState
            end    ;    of 9
            
            ;  MIO from Nanoletters
            ;
            10 : begin
              sState.ParticleName = sState.ParticleNameList[10]
              data_naming_nomanclature, sState
            end    ;    of 10  
            
            ;  Mag3300 from Nanoletters
            ;
            11 : begin
              sState.ParticleName = sState.ParticleNameList[11]
              data_naming_nomanclature, sState
            end    ;    of 11               

            ;  Ferumoxytol from Nanoletters
            ;
            12 : begin
              sState.ParticleName = sState.ParticleNameList[12]
              data_naming_nomanclature, sState
            end    ;    of 12

          endcase
          WIDGET_CONTROL, sState.wParticleTextField, SET_VALUE=sState.ParticleName
        end
        
        'ParticleNameText': begin
          WIDGET_CONTROL, sState.wParticleTextField, GET_VALUE=temp
          sState.ParticleName = temp
          data_naming_nomanclature, sState
;          print, temp
        end
        
        'ParticleConSLIDER': begin
          WIDGET_CONTROL, sState.wParticleConSlider, GET_VALUE=ParticleCon_temp
          WIDGET_CONTROL, sState.wParticleConTextField, SET_VALUE=STRING(ParticleCon_temp)
          ;update the particle concentration
          sState.ParticleCon = ParticleCon_temp
          sState.ParticleConStr = STRING(ParticleCon_temp)
          data_naming_nomanclature, sState
        end
        
        'ParticleConText': begin
          WIDGET_CONTROL, sState.wParticleConTextField, GET_VALUE=ParticleCon_temp
          sState.ParticleCon = FLOAT(ParticleCon_temp)
          sState.ParticleConStr = STRING(ParticleCon_temp, format='(f5.1)')
          WIDGET_CONTROL, sState.wParticleConSlider, SET_VALUE=LONG(ParticleCon_temp)
          data_naming_nomanclature, sState
        end
        
        'ConList': begin
          WIDGET_CONTROL, /HOURGLASS

          listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)

          case listValue of
            
            ;  hide style.
            ;
            0 : begin
              sState.ParticleConShow = 0
              data_naming_nomanclature, sState
            end    ;    of 0
            
            ;  hide style.
            ;
            1 : begin
              sState.ParticleConShow = 1
              data_naming_nomanclature, sState
            end    ;    of 1

          endcase     ;  of listValue
        end
        
        'BufferList': begin
          WIDGET_CONTROL, /HOURGLASS
          listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)
          case listValue of
            ;  None buffer.
            ;
            0 : begin
              sState.BufferName = sState.BufferNameList[0]
           
              data_naming_nomanclature, sState
            end    ;    of 0

            ;  Glycero buffer.
            ;
            1 : begin
              sState.BufferName = sState.BufferNameList[1]
              data_naming_nomanclature, sState
            end    ;    of 1

            2 : begin
              sState.BufferName = sState.BufferNameList[2]
              data_naming_nomanclature, sState
            end    ;    of 2
          
            3 : begin
              sState.BufferName = sState.BufferNameList[3]
              data_naming_nomanclature, sState
            end    ;    of 3
            
            4 : begin
              sState.BufferName = sState.BufferNameList[4]
              data_naming_nomanclature, sState
            end    ;    of 4

            5 : begin
              sState.BufferName = sState.BufferNameList[5]
              data_naming_nomanclature, sState
            end    ;    of 5            

            6 : begin
              sState.BufferName = sState.BufferNameList[6]
              data_naming_nomanclature, sState
            end    ;    of 6

            7 : begin
              sState.BufferName = sState.BufferNameList[7]
              data_naming_nomanclature, sState
            end    ;    of 7

            8 : begin
              sState.BufferName = sState.BufferNameList[8]
              data_naming_nomanclature, sState
            end    ;    of 8
          endcase
          WIDGET_CONTROL, sState.wBufferTextField, SET_VALUE=sState.BufferName
        end
        
        'BufferNameText': begin
          WIDGET_CONTROL, sState.wBufferTextField, GET_VALUE=temp
          sState.BufferName = temp
          data_naming_nomanclature, sState
          ;          print, temp
        end
        
        'BufferConSLIDER': begin
          WIDGET_CONTROL, sState.wBufferConSlider, GET_VALUE=BufferCon_temp
          WIDGET_CONTROL, sState.wBufferConTextField, SET_VALUE=STRING(BufferCon_temp)
          ;update the particle concentration
          sState.BufferCon = BufferCon_temp
          sState.BufferConStr = STRING(BufferCon_temp)
          data_naming_nomanclature, sState
        end

        'BufferConText': begin
          WIDGET_CONTROL, sState.wBufferConTextField, GET_VALUE=BufferCon_temp
          sState.BufferCon = FLOAT(BufferCon_temp)
          sState.BufferConStr = STRING(BufferCon_temp, format='(f5.1)')
          WIDGET_CONTROL, sState.wBufferConSlider, SET_VALUE=LONG(BufferCon_temp)
          data_naming_nomanclature, sState
        end
        
        'ShapeList': begin
          WIDGET_CONTROL, /HOURGLASS
          listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)
          case listValue of
            ;  Pulsed particle.
            ;
            0 : begin
              sState.Shape = sState.ShapeList[0]
              sState.FlatPortionShow = 1
              WIDGET_CONTROL, sState.wFlatPortionDroplist, SET_VALUE='Yes'
              data_naming_nomanclature, sState
            end    ;    of 0

            ;  Sine particle.
            ;
            1 : begin
              sState.Shape = sState.ShapeList[1]
              sState.FlatPortionShow = 0
              WIDGET_CONTROL, sState.wFlatPortionDroplist, SET_VALUE='No'
              data_naming_nomanclature, sState
            end    ;    of 1

            ;  triangular particle.
            ;
            2 : begin
              sState.Shape = sState.ShapeList[2]
              sState.FlatPortionShow = 0
              WIDGET_CONTROL, sState.wFlatPortionDroplist, SET_VALUE='No'
              data_naming_nomanclature, sState
            end    ;    of 2

            ;  circular particle.
            ;
            3 : begin
              sState.Shape = sState.ShapeList[3]
              sState.FlatPortionShow = 0
              WIDGET_CONTROL, sState.wFlatPortionDroplist, SET_VALUE='No'
              data_naming_nomanclature, sState
            end    ;    of 3

            ;  None particle.
            ;
            4 : begin
              sState.Shape = sState.ShapeList[4]
              sState.FlatPortionShow = 0
              WIDGET_CONTROL, sState.wFlatPortionDroplist, SET_VALUE='No'
              data_naming_nomanclature, sState
            end    ;    of 4

          endcase
          WIDGET_CONTROL, sState.wShapeTextField, SET_VALUE=sState.Shape
        end

        'ShapeText': begin
          WIDGET_CONTROL, sState.wShapeTextField, GET_VALUE=temp
          sState.Shape = temp
          data_naming_nomanclature, sState
          ;          print, temp
        end
        
        'FrequencySLIDER': begin
          WIDGET_CONTROL, sState.wFrequencySlider, GET_VALUE=Frequency_temp
          WIDGET_CONTROL, sState.wFrequencyTextField, SET_VALUE=STRING(Frequency_temp)
          ;update the temperature
          sState.Frequency = Frequency_temp
          sState.FrequencyStr = STRING(Frequency_temp)
          data_naming_nomanclature, sState
        end

        'FrequencyText': begin
          WIDGET_CONTROL, sState.wFrequencyTextField, GET_VALUE=Frequency_temp
          WIDGET_CONTROL, sState.wFrequencySlider, SET_VALUE=STRING(Frequency_temp)
          ;update the temperature
          sState.Frequency = Frequency_temp
          sState.FrequencyStr = STRING(Frequency_temp, format='(f5.1)')
          data_naming_nomanclature, sState
        end
        
        'AmpSLIDER': begin
          WIDGET_CONTROL, sState.wAmpSlider, GET_VALUE=Amp_temp
          WIDGET_CONTROL, sState.wAmpTextField, SET_VALUE=STRING(Amp_temp)
          ;update the temperature
          sState.Amp = Amp_temp
          sState.AmpStr = STRING(Amp_temp)
          data_naming_nomanclature, sState
        end

        'AmpText': begin
          WIDGET_CONTROL, sState.wAmpTextField, GET_VALUE=Amp_temp
          WIDGET_CONTROL, sState.wAmpSlider, SET_VALUE=STRING(Amp_temp)
          ;update the temperature
          sState.Amp = Amp_temp
          sState.AmpStr = STRING(Amp_temp, format='(f5.1)')
          data_naming_nomanclature, sState
        end

        'FlatPortionSLIDER': begin
          WIDGET_CONTROL, sState.wFlatPortionSlider, GET_VALUE=FlatPortion_temp
          WIDGET_CONTROL, sState.wFlatPortionTextField, SET_VALUE=STRING(FlatPortion_temp)
          ;update the temperature
          sState.FlatPortion = FlatPortion_temp
          sState.FlatPortionStr = STRING(FlatPortion_temp)
          data_naming_nomanclature, sState
        end

        'FlatPortionText': begin
          WIDGET_CONTROL, sState.wFlatPortionTextField, GET_VALUE=FlatPortion_temp
          WIDGET_CONTROL, sState.wFlatPortionSlider, SET_VALUE=STRING(FlatPortion_temp)
          ;update the temperature
          sState.FlatPortion = FlatPortion_temp
          sState.FlatPortionStr = STRING(FlatPortion_temp, format='(f5.1)')
          data_naming_nomanclature, sState
        end        
        
        'FlatPortionList': begin
          WIDGET_CONTROL, /HOURGLASS

          listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)

          case listValue of

            ;  show style.
            ;
            0 : begin
              sState.FlatPortionShow = 1
              data_naming_nomanclature, sState
            end    ;    of 0

            ;  hide style.
            ;
            1 : begin
              sState.FlatPortionShow = 0
              data_naming_nomanclature, sState
            end    ;    of 1

          endcase     ;  of listValue
          
        end
        
        'TempSLIDER': begin
          WIDGET_CONTROL, sState.wTempSlider, GET_VALUE=Temp_temp
          WIDGET_CONTROL, sState.wTempTextField, SET_VALUE=STRING(Temp_temp)
          ;update the temperature
          sState.Temp = Temp_temp
          sState.TempStr = STRING(Temp_temp)
          data_naming_nomanclature, sState
        end
        
        'TempText': begin
          WIDGET_CONTROL, sState.wTempTextField, GET_VALUE=Temp_temp
          WIDGET_CONTROL, sState.wTempSlider, SET_VALUE=STRING(Temp_temp)
          ;update the temperature
          sState.Temp = Temp_temp
          sState.TempStr = STRING(Temp_temp, format='(f5.1)')
          data_naming_nomanclature, sState
        end
        
        'TempShowList': begin
          WIDGET_CONTROL, /HOURGLASS

          listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)

          case listValue of

            ;  hide style.
            ;
            0 : begin
              sState.TempShow = 0
              data_naming_nomanclature, sState
            end    ;    of 0

            ;  hide style.
            ;
            1 : begin
              sState.TempShow = 1
              data_naming_nomanclature, sState
            end    ;    of 1

          endcase     ;  of listValue
        end

        'RepSLIDER': begin
          WIDGET_CONTROL, sState.wRepSlider, GET_VALUE=Rep_temp
          sState.Rep = LONG(Rep_temp)
          data_naming_nomanclature, sState
         ;print, 'rep = ', sState.Rep
        end
        
        'FileSelection': begin
          filename = dialog_pickfile(TITLE='Select MPI signal dat File', FILTER='*.txt', /MUST_EXIST, PATH=sState.directory)          
          IF STRLEN(filename) GT 1 THEN BEGIN
            sState.filename = filename
            sState.directory = FILE_DIRNAME(filename)
            WIDGET_CONTROL, sState.wFilenameField, SET_VALUE=FILE_BASENAME(filename, '.txt')

            ;GJ, 2022/6/18, parse file name
            datalist = STRSPLIT(FILE_BASENAME(filename, '.txt'), '_', /EXTRACT)

            ;find the date
            IF STRLEN(datalist[0]) EQ 6 THEN BEGIN
              yy_temp = FIX(STRMID(datalist[0], 0, 2))
              IF yy_temp GT 0 THEN BEGIN
                WIDGET_CONTROL, sState.wDateyearSlider, SET_VALUE=yy_temp
                ;update the year
                sState.yy = yy_temp
              ENDIF
              mm_temp = FIX(STRMID(datalist[0], 2, 2))
              IF mm_temp GT 0 THEN BEGIN
                WIDGET_CONTROL, sState.wDatemonthSlider, SET_VALUE=mm_temp
                ;update the month
                sState.mm = mm_temp
              ENDIF
              dd_temp = FIX(STRMID(datalist[0], 4, 2))
              IF dd_temp GT 0 THEN BEGIN
                WIDGET_CONTROL, sState.wDatedaySlider, SET_VALUE=dd_temp
                ;update the month
                sState.dd = dd_temp
              ENDIF
            ENDIF

            ;find the particle name
            IF N_ELEMENTS(datalist) GT 1 THEN BEGIN
              FOR i=0, N_ELEMENTS(sState.ParticleNameList)-1 DO BEGIN
                result_pN = STRPOS(datalist[1], sState.ParticleNameList[i])
                IF result_pN NE -1 THEN BEGIN
                  sState.ParticleName = sState.ParticleNameList[i]
                  WIDGET_CONTROL, sState.wParticleTextField, SET_VALUE=sState.ParticleName
                  WIDGET_CONTROL, sState.wParticleDroplist, SET_DROPLIST_SELECT=i
                ENDIF
              ENDFOR
            ENDIF

            ;find the buffer name
            FOR j=0, N_ELEMENTS(sState.BufferNameList)-1 DO BEGIN
              result_bN = STRPOS(sState.filename, sState.BufferNameList[j])
              IF result_bN NE -1 THEN BEGIN
                sState.BufferName = sState.BufferNameList[j]
                WIDGET_CONTROL, sState.wBufferTextField, SET_VALUE=sState.BufferName
                WIDGET_CONTROL, sState.wBufferDroplist, SET_DROPLIST_SELECT=j
              ENDIF
            ENDFOR

            ;find the pulse shape
            FOR k=0, N_ELEMENTS(sState.ShapeList)-1 DO BEGIN
              result_s = STRPOS(sState.filename, sState.ShapeList[k])
              IF result_s NE -1 THEN BEGIN
                sState.Shape = sState.ShapeList[k]
                WIDGET_CONTROL, sState.wShapeTextField, SET_VALUE=sState.Shape
                WIDGET_CONTROL, sState.wShapeDroplist, SET_DROPLIST_SELECT=k
              ENDIF
            ENDFOR

            ;check the temperature, shown
            TempShow=0
            FOR i=0, N_ELEMENTS(datalist)-1 DO BEGIN
              result_t = STRMATCH(datalist[i], 't*dg')
              IF result_t EQ 1 THEN BEGIN
                WIDGET_CONTROL, sState.wTempShowDroplist, SET_DROPLIST_SELECT=1
                sState.TempShow = 1
                TempShow = 1
              ENDIF
            ENDFOR
            IF TempShow EQ 0 THEN BEGIN
              WIDGET_CONTROL, sState.wTempShowDroplist, SET_DROPLIST_SELECT=0
              sState.TempShow = 0
            ENDIF

            ;update the filename
            data_naming_nomanclature, sState

            widget_control, sState.wModifyButton, sensitive=1
            widget_control, sState.wAnalyzeButton, sensitive=1
            widget_control, sState.wFileModifyButton, sensitive=1
            widget_control, sState.wFileAnalyzeButton, sensitive=1
          ENDIF

        end
        
        'Modify': begin
          WIDGET_CONTROL, sState.wTextField, GET_VALUE=data_name
          ;GJ, 2022/6/16, add name modification history log to the end of the file
          OPENR, lun, sState.filename, /get_lun
          ;read the first line
          temp_str=''
          readf, lun, temp_str
          FREE_LUN, lun
          add_str = SYSTIME()+': '+FILE_BASENAME(sState.filename, '.txt') + ' was changed to ' + data_name
          OPENW, lun, sState.filename, /APPEND, /get_lun
          WRITEU, lun, add_str
          FREE_LUN, lun
          FILE_MOVE, sState.filename, sState.directory+'\'+data_name+'.txt'
          sState.filename = sState.directory+data_name+'.txt'
          WIDGET_CONTROL, sState.wFilenameField, SET_VALUE='Renaming Done! Please select a new file'
          widget_control, sState.wModifyButton, sensitive=0
          widget_control, sState.wFileModifyButton, sensitive=0
        end
        
        'Analyze': begin
          WIDGET_CONTROL, sState.wTextField, GET_VALUE=data_name
          test_fit_para = read_pulsed_signal_dat(sState.filename, sState)
;          iimage, picture
;          new_pic = CONGRID(picture, 3, sState.xdim, sState.ydim)
;          WIDGET_CONTROL, sState.wDrawPic, GET_VALUE=drawpicid
;          WSET, drawpicid
;          TVSCL, new_pic, TRUE=1
          WIDGET_CONTROL, sState.wFilenameField, SET_VALUE=FILE_BASENAME(sState.filename)+' analyse Done!'
          widget_control, sState.wAnalyzeButton, sensitive=0
          widget_control, sState.wFileAnalyzeButton, sensitive=0
        end
        
        'BatchAnalyze': begin
          batch_pulsed_MH_curves_fitting, sState
          WIDGET_CONTROL, sState.wFilenameField, SET_VALUE='Batch Analyse Done or Cancelled! Please select a new folder'
        end

        'BatchEps': begin
          batch_read_UPEN_eps, sState
          WIDGET_CONTROL, sState.wFilenameField, SET_VALUE='Batch Ploting Done or Cancelled! Please select a new folder'
        end
        
        'ReadEps': begin
          data_array = read_UPEN_eps('', sState)
          WIDGET_CONTROL, sState.wFilenameField, SET_VALUE='Eps ploting Done or Cancelled! Please select a new eps file'
        end
        
        'BatchSens': begin
          ;@GJ, 2022/7/9, calculate sensitivity based on reading txt files
          batch_read_pulsed_sensitivity, sState
          WIDGET_CONTROL, sState.wFilenameField, SET_VALUE='Sensitivity Done or Cancelled! Please select a new folder'
        end
        
        
        'STYLELIST' : begin

            WIDGET_CONTROL, /HOURGLASS

            listValue = WIDGET_INFO( sEvent.id, /DROPLIST_SELECT)

            case listValue of

                ;  Shaded style.
                ;
                0 : begin
                    sState.oSimpleSurface->SetProperty, STYLE=2
                    dataFileRenaming_widgetsDraw, sState
                end    ;    of 0

                ;  Wire style.
                ;
                1 : begin
                    sState.oSimpleSurface->SetProperty, STYLE=1
                    dataFileRenaming_widgetsDraw, sState
                end    ;    of 1

                ;  Lego solid style.
                ;
                2 : begin
                    sState.oSimpleSurface->SetProperty, STYLE=6
                    dataFileRenaming_widgetsDraw, sState
                end    ;    of 2

            endcase     ;  of listValue

        end       ;   of DROPLIST

        ;  Handle the list event, Choose between 4 color scenarios.
        ;
        'COLORLIST' : begin

            WIDGET_CONTROL, /HOURGLASS

            listValue = WIDGET_INFO( sEvent.id, /LIST_SELECT)

            case listValue of

                ;  Texture Map
                ;
                0 : begin
                    sState.oSimpleSurface->SetProperty, $
                        TEXTURE_MAP=sState.oTextureImage
                    sState.oSimpleSurface->SetProperty, COLOR=[230,230,230], $
                       BOTTOM=[64, 192, 128]
                    dataFileRenaming_widgetsDraw, sState
                    sState.textureFlag=1
                end    ;  of 0

                ;  White.
                ;
                1 : begin
                    sState.oSimpleSurface->SetProperty, $
                        TEXTURE_MAP=OBJ_NEW()
                    sState.oSimpleSurface->SetProperty, COLOR=[200,200,200]
                    dataFileRenaming_widgetsDraw, sState
                    sState.textureFlag=0
                end    ;  of 1

                ;  Yellow.
                ;
                2 : begin
                    sState.oSimpleSurface->SetProperty, $
                        TEXTURE_MAP=OBJ_NEW()
                    sState.oSimpleSurface->SetProperty, COLOR=[200,200,0]
                    dataFileRenaming_widgetsDraw, sState
                    sState.textureFlag=0
                end    ;  of 2

                ;  Red.
                ;
                3 : begin
                    sState.oSimpleSurface->SetProperty, $
                        TEXTURE_MAP=OBJ_NEW()
                    sState.oSimpleSurface->SetProperty, COLOR=[200,0,0]
                    dataFileRenaming_widgetsDraw, sState
                    sState.textureFlag=0
                end    ;  of 3

            endcase

        end     ;   of LIST

        ;  Handle the rotation.
        ;
        'SLIDER' : begin

            WIDGET_CONTROL, sState.wXSlider, GET_VALUE=xDegree
            WIDGET_CONTROL, sState.wYSlider, GET_VALUE=yDegree
            WIDGET_CONTROL, sState.wZSlider, GET_VALUE=zDegree
            matFinal = FLTARR(3,3)
            matFinal = space3123(xDegree, yDegree, zDegree)

            sState.oRotationModel->GetProperty, TRANSFORM=t
            tempMat = FLTARR(3,3)
            tempMat[0:2, 0:2] = TRANSPOSE(t[0:2, 0:2])
            tempMat = TRANSPOSE(tempMat)
            rotMat = matFinal # tempMat

            ;  Find the Euler parameters 'e4' of rotMat
            ;  which is the rotation it takes to go from
            ;  the original (t) to the final (matFinal).
            ;
            e4 = 0.5 * SQRT(1.0 + rotMat[0,0] + $
                rotMat[1,1] + rotMat[2,2])

            ;  Find the unit vector of the single rotation axis
            ;  and the angle of rotation.
            ;
            if (e4 eq 0) then begin
                if (rotMat[0,0] eq 1) then begin
                    axisRot = [1, 0, 0]
                endif else if(rotMat[1,1] eq 1) then begin
                    axisRot = [0, 1, 0]
                endif else begin
                    axisRot = [0, 0, 1]
                endelse
                angleRot = 180.0
            endif else begin
                e1 = (rotMat[2,1] - rotMat[1,2])/(4.0*e4)
                e2 = (rotMat[0,2] - rotMat[2,0])/(4.0*e4)
                e3 = (rotMat[1,0] - rotMat[0,1])/(4.0*e4)
                modulusE = SQRT(e1*e1 + e2*e2 +e3*e3)
                if(modulusE eq 0.0) then begin
                    WIDGET_CONTROL, sEvent.top, $
                        SET_UVALUE=sState, /NO_COPY
                    RETURN
                endif
                axisRot = FLTARR(3)
                axisRot[0] = e1/modulusE
                axisRot[1] = e2/modulusE
                axisRot[2] = e3/modulusE
                angleRot = (2.0 * ACOS(e4)) * 180 / !DPI
            endelse

            for i = 0, 2 do begin
                if(ABS(axisRot[i]) lt 1.0e-6) then axisRot[i]=1.0e-6
            endfor
            sState.oRotationModel->Rotate, axisRot, angleRot
            dataFileRenaming_widgetsDraw, sState

        end    ;  of SLIDER

        'DRAW': begin

            ;  Expose.
            ;
            if (sEvent.type eq 4) then begin
                dataFileRenaming_widgetsDraw, sState
                WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
                RETURN
            endif


            ;  Handle trackball update
            ;
            bHaveTransform = sState.oTrack->Update(sEvent, TRANSFORM=qmat )
            if (bHaveTransform NE 0) then begin
                sState.oRotationModel->GetProperty, TRANSFORM=t
                mt = t # qmat
                sState.oRotationModel->SetProperty,TRANSFORM=mt
            endif

            ;  Button press.
            ;
            if (sEvent.type eq 0) then begin
                sState.btndown = 1B
                sState.drawWindowID->SetProperty, QUALITY=2
                WIDGET_CONTROL, sState.wDraw, /DRAW_MOTION
            endif    ;     of Button press

            ;  Button motion.
            ;
            if ((sEvent.type eq 2) and (sState.btndown eq 1B)) then begin

                if (bHaveTransform) then begin
                    ;  Reset the x, y, z axis angle values.
                    ;
                    sState.oRotationModel->GetProperty, TRANSFORM=transform
                    tempMat = FLTARR(3,3)
                    xyzAngles = FLTARR(3)
                    tempMat[0:2, 0:2] = transform[0:2, 0:2]
                    xyzAngles = Angle3123(tempMat)
                    WIDGET_CONTROL, sState.wXSlider, SET_VALUE=xyzAngles[0]
                    WIDGET_CONTROL, sState.wYSlider, SET_VALUE=xyzAngles[1]
                    WIDGET_CONTROL, sState.wZSlider, SET_VALUE=xyzAngles[2]
                    dataFileRenaming_widgetsDraw, sState
                endif

            endif     ;   of Button motion

            ;  Button release.
            ;
            if (sEvent.type eq 1) then begin
                if (sState.btndown EQ 1b) then begin
                    sState.drawWindowID->SetProperty, QUALITY=2
                    dataFileRenaming_widgetsDraw, sState
                endif
                sState.btndown = 0B
                WIDGET_CONTROL, sState.wDraw, DRAW_MOTION=0
            endif
        end                 ;     of DRAW
        'TEXT': begin
        end


        "QUIT": begin

            ;  Restore the info structure before destroying event.top
            ;
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY

            ;  Destroy widget hierarchy.
            ;
            WIDGET_CONTROL, sEvent.top, /DESTROY

            RETURN
        end

        ELSE :  begin
            PRINT, 'Case Statement found no matches'
        end

    endcase

    ; Restore the info structure
    ;
    WIDGET_CONTROL, sEvent.top, Set_UValue=sState, /No_Copy
end               ; of dataFileRenaming_widgetsEvent

; -----------------------------------------------------------------------------
;
;  Purpose:  Cleanup procedure
;
pro dataFileRenaming_widgetsCleanup, $
    wTopBase      ; IN: top level base associated with the cleanup

    ;  Get the color table saved in the window's user value.
    ;
    WIDGET_CONTROL, wTopBase, GET_UVALUE=sState,/No_Copy

    ;  Destroy the top objects
    ;
    OBJ_DESTROY, sState.oStaticModel
;    OBJ_DESTROY, sState.oContainer
;    OBJ_DESTROY, sState.oTextureImage
;    OBJ_DESTROY, sState.oPalette

    ;  Restore the previous color table.
    ;
;    TVLCT, sState.colorTable

    ;  Map the group leader base if it exists.
    ;
    if (WIDGET_INFO(sState.groupBase, /VALID_ID)) then $
        WIDGET_CONTROL, sState.groupBase, /MAP

end   ; of dataFileRenaming_widgetsCleanup


; -----------------------------------------------------------------------------
;
;  Purpose:  Main procedure of the widgets datefile renaming
;
pro dataFileRenaming_widgets, $
    GROUP=group, $     ; IN: (opt) group identifier
    RECORD_TO_FILENAME=record_to_filename, $
    DEBUG=debug, $     ; IN: (opt)
    APPTLB = appTLB    ; OUT: (opt) TLB of this application

    ; Check the validity of the group identifier
    ;
    ngroup = N_ELEMENTS(group)
    if (ngroup NE 0) then begin
        check = WIDGET_INFO(group, /VALID_ID)
        if (check NE 1) then begin
            print,'Error, the group identifier is not valid'
            print, 'Return to the main application'
            RETURN
        endif
        groupBase = group
    endif else groupBase = 0L


    ;  Get the screen size.
    ;
    Device, GET_SCREEN_SIZE = screenSize

    ;  Set up dimensions of the drawing (viewing) area.
    ;
    ydim = screenSize[0]*0.25
    xdim = ydim * 2.17

    ;  Make the system have a maximum of 256 colors
    ;
    numcolors = !d.N_COLORS
    if( (( !D.NAME EQ 'X') or (!D.NAME EQ 'MAC')) $
        and (!d.N_COLORS GE 256L)) then $
        DEVICE, PSEUDO_COLOR=8

    DEVICE, DECOMPOSED=0, BYPASS_TRANSLATION=0

    ;  Get the current color table
    ;
    TVLCT, savedR, savedG, savedB, /GET

    ;  Build color table from color vectors
    ;
    colorTable = [[savedR],[savedG],[savedB]]

    ;  Load an initial file .
    ;
    file = 'abnorm.dat'
    demo_getdata, NewImage, FILENAME=file, /TWO_DIM
    data=SMOOTH(NewImage[*,*,7], 3, /EDGE_TRUNCATE)

    ; Define a main widget base.
    ;
    if (N_ELEMENTS(group) EQ 0) then begin
        wTopBase = WIDGET_BASE(TITLE="Data File Nomanclature (MPS)", /COLUMN, $
            /TLB_KILL_REQUEST_EVENTS, $
            MBAR=barBase)
    endif else begin
        wTopBase = WIDGET_BASE(TITLE="Data File Nomanclature (MPS)", /COLUMN, $
            /TLB_KILL_REQUEST_EVENTS, $
            GROUP_LEADER=group, $
            MBAR=barBase)
    endelse

        ;  Create the quit button
        ;
        wFileButton = WIDGET_BUTTON(barBase, VALUE= 'File', /MENU)
            
            wFileSelectButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Select a file', UVALUE='FileSelection')
            wFileModifyButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Modify filename', UVALUE='Modify')
            widget_control, wFileModifyButton, sensitive=0 
            wFileAnalyzeButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Analyze a file', UVALUE='Analyze')
            widget_control, wFileAnalyzeButton, sensitive=0 
            wBatchAnalyzeButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Batch Analyze files', UVALUE='BatchAnalyze')
            wQuitButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Quit', UVALUE='QUIT')


        ;  Create the first child of the graghic base.
        ;
        wSubBase = WIDGET_BASE(wTopBase, COLUMN=5)


        w1stBase = WIDGET_BASE(wSubBase, /BASE_ALIGN_CENTER, $
          XPAD=5, YPAD=5, $
          /FRAME, /COLUMN)
        wDateSliderBase = WIDGET_BASE(w1stBase, /COLUMN, $
          /FRAME, YPAD=8, XPAD=8)
        ; Initial rotation values, also used below when
        ; initializing model
        current_yy = STRING(systime(/julian), format='(c(CYI))')
        yy = FIX(STRMID(current_yy, 2, 2))
        current_month = STRING(systime(/julian), format='(c(CMOI))')
        mm = FIX(current_month)
        current_day = STRING(systime(/julian), format='(c(CDI))')
        dd = FIX(current_day)
        wDateSliderLabel = WIDGET_LABEL(wDateSliderBase, $
          VALUE='Date', /ALIGN_CENTER)
        wDateyearSlider = WIDGET_SLIDER(wDateSliderBase, $
          UVALUE='YearSLIDER', $
          VALUE=yy, MINIMUM=10, MAXIMUM=99)
        wDateSlideryearLabel = WIDGET_LABEL(wDateSliderBase, $
          VALUE='year')
        wDatemonthSlider = WIDGET_SLIDER(wDateSliderBase, $
          UVALUE='MonthSLIDER', $
          VALUE=mm, MINIMUM=1, MAXIMUM=12)
        wDateSlidermonthLabel = WIDGET_LABEL(wDateSliderBase, $
          VALUE='month')
        wDateDaySlider = WIDGET_SLIDER(wDateSliderBase, $
          VALUE=dd, MINIMUM=1, MAXIMUM=31, $
          UVALUE='DaySLIDER')
        wDateSliderDayLabel = WIDGET_LABEL(wDateSliderBase, $
          VALUE='Day')
          
        w2ndBase = WIDGET_BASE(wSubBase, /BASE_ALIGN_CENTER, $
          XPAD=5, YPAD=5, $
          /FRAME, /COLUMN)
        wParticleBase = WIDGET_BASE(w2ndBase, /COLUMN, $
          /FRAME, YPAD=8, XPAD=8)
        ;@GJ, 2022/06/01, particle name
        wParticleLabel = WIDGET_LABEL(wParticleBase, $
          VALUE='Particle name')
        ParticleName = 'Synomag'
        ParticleNameList = [ParticleName, 'Perimag', 'NE10', 'NE25', 'NE50', 'Resovist', 'VivoTrak', 'Dox-CACN', 'FeCo@C-PEG', 'IO', 'MIO', 'Mag3300', 'Ferumxytol']
        wParticleDroplist = WIDGET_DROPLIST(wParticleBase, $
          VALUE=ParticleNameList, $
          UVALUE='ParticleList')
        wParticleTextField = widget_text(wParticleBase, $
          /EDITABLE, $
          UVALUE='ParticleNameText', $
          VALUE=ParticleName)
        wParticleConLabel = WIDGET_LABEL(wParticleBase, $
          VALUE='Concentration (mg/ml Fe)')
        ParticleCon = 1
        wParticleConSlider = WIDGET_SLIDER(wParticleBase, $
          VALUE=ParticleCon, MINIMUM=0, MAXIMUM=100, $
          UVALUE='ParticleConSLIDER')
        wParticleConTextField = widget_text(wParticleBase, $
          /EDITABLE, $
          UVALUE='ParticleConText', $
          VALUE=STRING(ParticleCon))
        ParticleConShow=0
        wParticleConShowLabel = WIDGET_LABEL(wParticleBase, $
          VALUE='Concenration display:')
        wParticleConDroplist = WIDGET_DROPLIST(wParticleBase, $
          VALUE=['No', 'Yes'], $
          UVALUE='ConList')
          
        w3rdBase = WIDGET_BASE(wSubBase, /BASE_ALIGN_CENTER, $
          XPAD=5, YPAD=5, $
          /FRAME, /COLUMN)
        wBufferBase = WIDGET_BASE(w3rdBase, /COLUMN, $
          /FRAME, YPAD=8, XPAD=8)
        ;@GJ, 2022/06/01, buffer name
        wBufferLabel = WIDGET_LABEL(wBufferBase, $
          VALUE='Buffer name')
        BufferName='None'
        BufferNameList = [BufferName, 'Glycero', 'Gelatin', 'Sugar', 'PBS', 'ATP', 'Water', 'NaCl', 'DMEM']
        wBufferDroplist = WIDGET_DROPLIST(wBufferBase, $
          VALUE=BufferNameList, $
          UVALUE='BufferList')
        wBufferTextField = widget_text(wBufferBase, $
          /EDITABLE, $
          UVALUE='BufferNameText', $
          VALUE=BufferName)
        BufferCon = 1
        wBufferConLabel = WIDGET_LABEL(wBufferBase, $
          VALUE='Buffer (%)')
        wBufferConSlider = WIDGET_SLIDER(wBufferBase, $
          VALUE=BufferCon, MINIMUM=0, MAXIMUM=100, $
          UVALUE='BufferConSLIDER')
        wBufferConTextField = widget_text(wBufferBase, $
          /EDITABLE, $
          UVALUE='BufferConText', $
          VALUE=STRING(BufferCon))
          
        w4thBase = WIDGET_BASE(wSubBase, /BASE_ALIGN_CENTER, $
          XPAD=5, YPAD=5, $
          /FRAME, /COLUMN)
        wFieldBase = WIDGET_BASE(w4thBase, /COLUMN, $
          /FRAME, YPAD=8, XPAD=8)
        ;@GJ, 2022/06/01, buffer name
        ;@GJ, 2022/06/01, excitation field shape
        Shape = 'pulsed'
        ShapeList = [Shape, 'sine', 'triangular', 'circular', 'None']
        wShapeLabel = WIDGET_LABEL(wFieldBase, $
          VALUE='Excitation Field')
        wShapeDroplist = WIDGET_DROPLIST(wFieldBase, $
          VALUE=ShapeList, $
          UVALUE='ShapeList')
        wShapeTextField = widget_text(wFieldBase, $
          /EDITABLE, $
          UVALUE='ShapeText', $
          VALUE=Shape)
        Frequency=2
        wFrequencyConLabel = WIDGET_LABEL(wFieldBase, $
          VALUE='Frequency (kHz)')
        wFrequencySlider = WIDGET_SLIDER(wFieldBase, $
          VALUE=Frequency, MINIMUM=0, MAXIMUM=100, $
          UVALUE='FrequencySLIDER')
        wFrequencyTextField = widget_text(wFieldBase, $
          /EDITABLE, $
          UVALUE='FrequencyText', $
          VALUE=STRING(Frequency))
        Amp=1
        wAmpLabel = WIDGET_LABEL(wFieldBase, $
          VALUE='Amplitude (mT)')
        wAmpSlider = WIDGET_SLIDER(wFieldBase, $
          VALUE=Amp, MINIMUM=0, MAXIMUM=200, $
          UVALUE='AmpSLIDER')
        wAmpTextField = widget_text(wFieldBase, $
          /EDITABLE, $
          UVALUE='AmpText', $
          VALUE=STRING(Amp))
        FlatPortion=90
        wFlatPortionLabel = WIDGET_LABEL(wFieldBase, $
          VALUE='Flat Portion (%)')
        wFlatPortionSlider = WIDGET_SLIDER(wFieldBase, $
          VALUE=FlatPortion, MINIMUM=0, MAXIMUM=100, $
          UVALUE='FlatPortionSLIDER')
        wFlatPortionTextField = widget_text(wFieldBase, $
          /EDITABLE, $
          UVALUE='FlatPortionText', $
          VALUE=STRING(FlatPortion))
        FlatPortionShow=1
        wFlatPortionShowLabel = WIDGET_LABEL(wFieldBase, $
          VALUE='Flat Portion display:')
        wFlatPortionDroplist = WIDGET_DROPLIST(wFieldBase, $
          VALUE=['Yes', 'No'], $
          UVALUE='FlatPortionList')
        
        
        w5thBase = WIDGET_BASE(wSubBase, /BASE_ALIGN_CENTER, $
          XPAD=5, YPAD=5, $
          /FRAME, /COLUMN)
        wTempBase = WIDGET_BASE(w5thBase, /COLUMN, $
          /FRAME, YPAD=8, XPAD=8)
        ;@GJ, 2022/06/01, temperature
        Temp=25
        wTempLabel = WIDGET_LABEL(wTempBase, $
          VALUE='Temperature (degree)')
        wTempSlider = WIDGET_SLIDER(wTempBase, $
          UVALUE='TempSLIDER', $
          VALUE=Temp, MINIMUM=-100, MAXIMUM=100)
        wTempTextField = widget_text(wTempBase, $
          /EDITABLE, $
          UVALUE='TempText', $
          VALUE=STRING(Temp))
        TempShow=0
        wTempShowLabel = WIDGET_LABEL(wTempBase, $
          VALUE='Temperature display:')
        wTempShowDroplist = WIDGET_DROPLIST(wTempBase, $
          VALUE=['No', 'Yes'], $
          UVALUE='TempShowList')
        ;@GJ, 2022/06/01, repeation
        Rep=0
        wRepLabel = WIDGET_LABEL(wTempBase, $
          VALUE='Repetition (times)')
        wRepSlider = WIDGET_SLIDER(wTempBase, $
          UVALUE='RepSLIDER', $
          VALUE=Rep, MINIMUM=0, MAXIMUM=100)
          
         wLowBase = WIDGET_BASE(wTopBase, /BASE_ALIGN_CENTER, $
          XPAD=5, YPAD=5, $
          /FRAME, /ROW) 

        
            ;  Create a base for the left column.
            ;
            wLeftBase = WIDGET_BASE(wLowBase, /BASE_ALIGN_CENTER, $
               XPAD=5, YPAD=5, $
               /FRAME, /COLUMN)

             IF mm LT 10 THEN BEGIN
               IF dd LT 10 THEN BEGIN
                 data_name_temp = STRING(yy)+'0'+STRING(mm)+'0'+STRING(dd)+'_'+ParticleName+'_'+Shape+'f'+STRING(Frequency)+'ka'+STRING(Amp)+'fp'+STRING(FlatPortion)+'%'
               ENDIF ELSE BEGIN
                 data_name_temp = STRING(yy)+'0'+STRING(mm)+STRING(dd)+'_'+ParticleName+'_'+Shape+'f'+STRING(Frequency)+'ka'+STRING(Amp)+'fp'+STRING(FlatPortion)+'%'
               ENDELSE
             ENDIF ELSE BEGIN
               IF dd LT 10 THEN BEGIN
                 data_name_temp = STRING(yy)+STRING(mm)+'0'+STRING(dd)+'_'+ParticleName+'_'+Shape+'f'+STRING(Frequency)+'ka'+STRING(Amp)+'fp'+STRING(FlatPortion)+'%'
               ENDIF ELSE BEGIN
                 data_name_temp = STRING(yy)+STRING(mm)+STRING(dd)+'_'+ParticleName+'_'+Shape+'f'+STRING(Frequency)+'ka'+STRING(Amp)+'fp'+STRING(FlatPortion)+'%'
               ENDELSE
              ENDELSE
            data_name = STRCOMPRESS(data_name_temp, /REMOVE_ALL)
            wStatusLabel = WIDGET_LABEL(wLeftBase, $
                VALUE='Editable Name')
            
            wTextField = widget_text(wLeftBase, $
                /EDITABLE, XSIZE=50, $
                UVALUE='TEXT', $
                VALUE=data_name)
       
              wFilenameField = widget_text(wLeftBase, $
                /EDITABLE, XSIZE=50, $
                UVALUE='Filename', $
                VALUE='Please select a file!')
            
            wButton1Base = WIDGET_BASE(wLeftBase, /BASE_ALIGN_CENTER, $
                XPAD=5, YPAD=5, $
                /FRAME, /ROW)
            wFileSelectionButton = WIDGET_BUTTON(wButton1Base, $
                VALUE='File Selection', UVALUE='FileSelection')

            wModifyButton = WIDGET_BUTTON(wButton1Base, $
                VALUE='Modify FileName', UVALUE='Modify')
            widget_control, wModifyButton, sensitive=0     
           
            wButton2Base = WIDGET_BASE(wLeftBase, /BASE_ALIGN_CENTER, $
              XPAD=5, YPAD=5, $
              /FRAME, /ROW)
            wAnalyzeButton = WIDGET_BUTTON(wButton2Base, $
              VALUE='Analyze File', UVALUE='Analyze')
            widget_control, wAnalyzeButton, sensitive=0
            
            wBatchAnalyzeButton = WIDGET_BUTTON(wButton2Base, $
              VALUE='Batch Analyze', UVALUE='BatchAnalyze')
            widget_control, wBatchAnalyzeButton, sensitive=1
            
            
            wButton4Base = WIDGET_BASE(wLeftBase, /BASE_ALIGN_CENTER, $
              XPAD=5, YPAD=5, $
              /FRAME, /ROW)
            wBatchSensButton = WIDGET_BUTTON(wButton4Base, $
              VALUE='Batch analyze for senxitivity', UVALUE='BatchSens')
            widget_control, wBatchSensButton, sensitive=1

            wButton3Base = WIDGET_BASE(wLeftBase, /BASE_ALIGN_CENTER, $
              XPAD=5, YPAD=5, $
              /FRAME, /ROW)
            wReadEpsButton = WIDGET_BUTTON(wButton3Base, $
              VALUE='Read eps file', UVALUE='ReadEps')
            widget_control, wReadEpsButton, sensitive=1
            wBatchEpsButton = WIDGET_BUTTON(wButton3Base, $
              VALUE='Batch Plot eps files', UVALUE='BatchEps')
            widget_control, wBatchEpsButton, sensitive=1            
            

            
            
            wRightBase = WIDGET_BASE(wLowBase, /BASE_ALIGN_CENTER, $
              XPAD=5, YPAD=5, $
              /FRAME, /ROW)
            
            wDrawPic = WIDGET_DRAW(wRightBase, $
                XSIZE=xdim, YSIZE=ydim, UVALUE='DRAWPIC')
;            ;  Create a base for the right column.
;            ;
;            wRightBase = WIDGET_BASE(wSubBase, /COLUMN)
;
;                wTableBase = WIDGET_BASE(wRightBase, /COLUMN, $
;                    /FRAME, YPAD=0, XPAD=0)
;                    sz = SIZE(data)
;                    wLabel = WIDGET_LABEL(wTableBase, $
;                        VALUE='View or modify surface data with the table widget')
;                    wDataTable = WIDGET_TABLE(wTableBase, $
;                        UVALUE='TABLE', $
;                        VALUE=data, /EDITABLE, /ALL_EVENTS, $
;                        XSIZE=sz[1], YSIZE=sz[2], $
;                        X_SCROLL_SIZE=5, Y_SCROLL_SIZE=3 )
;
;
;                wDrawBase = WIDGET_BASE(wRightBase, /COLUMN, $
;                    /FRAME, UVALUE=-1)
;
;                    wDraw = WIDGET_DRAW(wDrawBase, $
;                        XSIZE=xdim, YSIZE=ydim, /BUTTON_EVENTS, $
;                        /EXPOSE_EVENTS, UVALUE='DRAW', $
;                        RETAIN=0, $
;                        GRAPHICS_LEVEL=2)
;
;        ;  Create tips texts.
        ;
        wStatusBase = WIDGET_BASE(wTopBase, MAP=0, /ROW)

    ;  Realize the widget hierarchy.
    ;
    WIDGET_CONTROL, wTopBase, /REALIZE

    ;  Returns the top level base in the appTLB keyword.
    ;
    appTLB = wTopBase

    ;  Get the tips
    ;
    sText = demo_getTips(demo_filepath('widgets.tip', $
                             SUBDIR=['examples','demo', 'demotext']), $
                         wTopBase, $
                         wStatusBase)

    WIDGET_CONTROL, wTopBase, SENSITIVE=0

    ; Determine the window value of plot window, wDraw1.
    ;
;    WIDGET_CONTROL, wDraw, GET_VALUE=drawWindowID


    ;  Set the view such that the surface is
    ;  contained within a box defined by
    ;  by radx and rady ( view normal coordinates).
    ;
    sqr3 = SQRT(3)/1.5   ; length of a diagonal in a cube
    radx = SQRT(2.0)       ; viewport in normal coordinates
    rady = SQRT(2.0)

    ;  Select the normal corrdinates of the
    ;  orthogonal axes location within the volume.
    ;  example:
    ;  middle :  axesX, axesY, axesZ = -0.5
    ;  lowest data values : axesX, axesY, axesZ = 0.0
    ;  highest data values : axesX, axesY, axesZ = -1.0
    ;
    axesX = -0.5
    axesY = -0.5
    axesZ = -0.5

    xMargin = (1.0-sqr3)/2.0
    yMargin = (1.0-sqr3)/2.0
    xv = ((xMargin)*radx + axesX)
    yv = ((yMargin)*rady + axesY)
    width = 1.0 - 2.0 * xMargin* radx
    height = 1.0 - 2.0 * yMargin * radY

    myview = [xv, yv, width, height]

    ;  Create view.
    ;
    oView = OBJ_NEW('idlgrview', PROJECTION=2, EYE=3, $
        ZCLIP=[1.5, -1.5], VIEWPLANE_RECT=myview, COLOR=[0,0,0])

    ;  Create model.
    ;
    oStaticModel = OBJ_NEW('idlgrmodel')
    oMovableModel = OBJ_NEW('idlgrmodel')
    oRotationModel = OBJ_NEW('idlgrmodel')
    oScalingModel = OBJ_NEW('idlgrmodel')
    oTranslationModel = OBJ_NEW('idlgrmodel')

    oStaticModel->Add, oMovableModel
    oMovableModel->Add, oRotationModel
    oRotationModel->Add, oScalingModel
    oScalingModel->Add, oTranslationModel

    sc = 0.7
    oStaticModel->Scale, sc, sc, sc

    ;  Create light.
    ;
    oLight1 = OBJ_NEW('idlgrLight', TYPE=1, INTENSITY=0.5, $
        LOCATION=[1.5, 0, 1])
    oLight2 = OBJ_NEW('idlgrLight', TYPE=0, $
        INTENSITY=0.75)
    oStaticModel->Add, oLight1
    oStaticModel->Add, oLight2

    ;  Compute coordinate conversion to normalize.
    ;
    z = data
    sz = SIZE(z)
    maxx = sz[1] - 1
    maxy = sz[2] - 1
    maxz = MAX(z,min=minz)
    xs = [axesX,1.0/maxx]
    ys = [axesY,1.0/maxy]
    minz2 = minz - 1
    maxz2 = maxz + 1
    zs = [(-minz2/(maxz2-minz2))+axesZ, 1.0/(maxz2-minz2)]

    oPalette = OBJ_NEW('idlgrPalette')
    oPalette->LOADCT, 25
    oTextureImage = OBJ_NEW('idlgrImage', BYTSCL(z), PALETTE=oPalette)

    ;  Create the surface.
    ;
    oSimpleSurface = OBJ_NEW('IDLgrSurface', data, $
        TEXTURE_MAP=oTextureImage, $
        STYLE=2, SHADING=1, $
        /USE_TRIANGLES, $
        ;COLOR=[60,60,255], $
        ;BOTTOM=[64,192,128], $
        COLOR=[230, 230, 230], BOTTOM=[64, 192, 128], $
        XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)

    oTranslationModel->Add, oSimpleSurface

    ;  Rotate the original display.
    ;
    ;  Place the model in the view.
    ;
    oView->Add, oStaticModel

    ;  Add the trackball object for interactive change
    ;  of the scene orientation
    ;
    oTrack = OBJ_NEW('Trackball', [xdim/2.0, ydim/2.0], xdim/2.0)

    oContainer = OBJ_NEW('IDLgrContainer')
    oContainer->Add, oView
    oContainer->Add, oTrack

    ;  Create the info structure
    ;
    sState={ $
        yy: yy, $
        mm: mm, $
        dd: dd, $
        ParticleName: ParticleName, $
        ParticleNameList: ParticleNameList, $
        ParticleCon: ParticleCon, $
        ParticleConStr: STRING(ParticleCon), $
        ParticleConShow: ParticleConShow, $
        BufferName: BufferName, $
        BufferNameList: BufferNameList, $
        BufferCon: BufferCon, $
        BufferConStr: STRING(BufferCon), $
        Shape: Shape, $
        ShapeList: ShapeList, $
        Frequency: Frequency, $
        FrequencyStr: STRING(Frequency), $
        Amp: Amp, $
        AmpStr: STRING(Amp), $
        FlatPortion: FlatPortion, $
        FlatPortionStr: STRING(FlatPortion), $
        FlatPortionShow: FlatPortionShow, $
        Temp: Temp, $
        TempStr: STRING(Temp), $
        TempShow: TempShow, $
        Rep: Rep, $
        data_name: data_name, $
        filename: '', $
        directory: 'C:\D_drive\MPI_Tianjie\Xinfeng\20220616\', $
        n_modify: 0, $
        xdim: xdim, $
        ydim: ydim, $
        wDateyearSlider: wDateyearSlider, $
        wDatemonthSlider: wDatemonthSlider, $
        wDateDaySlider:  wDateDaySlider, $
        wParticleDroplist: wParticleDroplist, $
        wParticleTextField: wParticleTextField, $
        wParticleConSlider: wParticleConSlider, $
        wParticleConTextField: wParticleConTextField, $
        wParticleConDroplist: wParticleConDroplist, $
        wBufferDroplist: wBufferDroplist, $
        wBufferTextField: wBufferTextField, $
        wBufferConSlider: wBufferConSlider, $
        wBufferConTextField: wBufferConTextField, $
        wShapeDroplist: wShapeDroplist, $
        wShapeTextField:wShapeTextField, $
        wFrequencySlider: wFrequencySlider, $
        wFrequencyTextField: wFrequencyTextField, $
        wAmpSlider: wAmpSlider, $
        wAmpTextField: wAmpTextField, $
        wFlatPortionSlider: wFlatPortionSlider, $
        wFlatPortionTextField: wFlatPortionTextField, $
        wFlatPortionDroplist: wFlatPortionDroplist, $
        wTempSlider: wTempSlider, $
        wTempTextField: wTempTextField, $
        wTempShowDroplist: wTempShowDroplist, $
        wRepSlider: wRepSlider, $
        wFilenameField: wFilenameField, $
        wModifyButton: wModifyButton, $
        wAnalyzeButton: wAnalyzeButton, $
        wFileSelectButton: wFileSelectButton, $
        wFileModifyButton: wFileModifyButton, $
        wFileAnalyzeButton: wFileAnalyzeButton, $
        wBatchAnalyzeButton: wBatchAnalyzeButton, $
        wReadEpsButton: wReadEpsButton, $ ;@GJ, 2022/7/4, read a single eps file
        wBatchEpsButton: wBatchEpsButton, $  ;add a button to plot all eps files in filename order
        wBatchSensButton: wBatchSensButton, $ ;GJ, 2022/7/9, read txt files and calculate sensitivity
        wDrawPic: wDrawPic, $
        BtnDown: 0, $                           ; mouse button down flag
        OTrack: oTrack, $                       ; Trackball object
        OContainer: oContainer, $               ; Container object
        OView: oView, $                         ; View object
        OStaticModel: oStaticModel, $           ; Models objects
        ORotationModel: oRotationModel, $
        OSimpleSurface: oSimpleSurface, $       ; Surface object
 ;       WDataTable: wDataTable, $               ; Widget table ID
 ;       WXSlider: wXSlider, $                   ; Widget sliders ID
;        WYSlider: wYSlider, $
;        WZSlider: wZSlider, $
;        DrawWindowID: drawWindowID, $           ; Window ID
;        WDraw: wDraw, $                         ; Widget draw ID
        wTextField: wTextField, $               ; C-Widget text ID
;        OTextureImage: oTextureImage, $         ; Texture image object
;        OPalette: oPalette, $                   ; Palette for image
;        textureFlag: 1, $                       ; Texture mapping flag
;        colorTable: colorTable, $               ; color table to restore
        debug: keyword_set(debug), $            ; debug flag
        groupBase: groupBase $                  ; Base of Group Leader
    }

    dataFileRenaming_widgetsDraw, sState

    ;  Register the info structure in the user value of the top-level base
    ;
    WIDGET_CONTROL, wTopBase, SET_UVALUE=sState, /NO_COPY

    WIDGET_CONTROL, wTopBase, SENSITIVE=1

    ;  Map the top level base.
    ;
    WIDGET_CONTROL, wTopBase, MAP=1

    XMANAGER, "dataFileRenaming_widgets", wTopBase, $
        /NO_BLOCK, $
        EVENT_HANDLER="dataFileRenaming_widgetsEvent", CLEANUP="dataFileRenaming_widgetsCleanup"

end   ;  main procedure


