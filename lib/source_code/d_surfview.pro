; $Id: //depot/Release/ENVI51_IDL83/idl/idldir/examples/demo/demosrc/d_surfview.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  FILE:
;       d_surfview.pro
;
;  CALLING SEQUENCE: d_surfview
;
;  PURPOSE:
;       This application demonstrates the capabilities of a surface.
;
;  MAJOR TOPICS: Visualization
;
;  CATEGORY:
;       IDL Demo System
;
;  INTERNAL FUNCTIONS and PROCEDURES:
;       fun d_surfviewToggleOffOn          - Toggle between 'off' and 'on'
;       pro d_surfviewAnimateSurface50     - Animate the surface
;       pro d_surfviewAddTracePoint        - Add a tracing point
;       pro d_surfviewEvent                - Event handler
;       pro d_surfviewCleanup              - Cleanup
;       pro d_surfview                     - Main procedure
;
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       pro idlexrotator__define    - Create a trackball-like object
;       pro demo_gettips            - Read the tip file and create widgets
;       pro demo_puttips            - Change tips text
;       surfview.tip                - "Tips" Text
;
;  REFERENCE: IDL Reference Guide, IDL User's Guide
;
;  MODIFICATION HISTORY:
;       9/96,   DD   - Written.
;       10/96,  DAT  - New GUI.
;       1998    PCS  - Rotation constrictions.  Changes to GUI.
;-
;----------------------------------------------------------------------------
;
;  Purpose :  Toggle the string 'activate' and 'deactivate' of a widget
;             value.  Returns -1 on failure.
;
function d_surfviewToggleOffOn, $
    widgetID      ; IN: widget identifier for which its value is the name

    WIDGET_CONTROL, widgetID, GET_VALUE=name

    ;  Toggle the 'activate' and 'deactivate' strings.
    ;  Otherwise, the function retuns -1.
    ;
    offPosition = STRPOS(name, 'Deactivate ')
    if (offPosition ne -1) then begin
        name = 'Activate ' + STRMID(name, 11)
        returnFlag = 0
    endif else begin
        onPosition = STRPOS(name, 'Activate ')
        if (onPosition ne -1) then begin
            name = 'Deactivate ' + STRMID(name, 9)
            returnFlag = 1
        endif else begin
            returnFlag = -1
        endelse
    endelse

    WIDGET_CONTROL, widgetID, SET_VALUE=name

    RETURN, returnFlag

end

;----------------------------------------------------------------------------
;
;  Purpose:  Animate the surface object.
;
pro d_surfviewAnimateSurface, $
    sState, $            ; IN: sState structure
    scenario             ; IN: index indicating the specified animation

    ;  Get the x and y dimension of the surface object.
    ;
    xdim = sState.xSizeData
    ydim = sState.ySizeData

    case scenario of

        ;  This animation deforms the middle of the surface back an forth
        ;  by translating the x and y coordinates.
        ;
        0 : begin
            x = FLTARR(xdim, ydim)
            for i = 0, 10 do begin
                floatIndex = FLOAT(i)
                for j = 0, xdim-1 do x[j, *] = j
                for j = 0, ydim-1 do x[*, j] = x[*,j] $
                    + floatIndex*SIN(j*(!PI/20.0))
                sState.oSurface->SetProperty, DATAX = x
                sState.oWindow->Draw, sState.oView
            endfor
            for i = 9, 0, -1 do begin
                floatIndex = FLOAT(i)
                for j = 0, xdim-1 do x[j, *] = j
                for j = 0, ydim-1 do x[*, j] = x[*,j] $
                    + floatIndex*SIN(j*(!pi/20.0))
                sState.oSurface->SetProperty, DATAX=x
                sState.oWindow->Draw, sState.oView
            endfor
        end

    endcase

end   ; of Animate

;----------------------------------------------------------------------------
;
;  Purpose:  Add a data point to the trace path.
;
pro d_surfviewAddTracePoint, $
    oTrace, $        ;  IN: polyline object of the trace
    xyzVector, $     ;  IN: x, y, z, vector of the new data point
    status           ;  OUT:  -1: failure, 1: success

    ;  Verify the validity of oTrace.
    ;
    if (OBJ_VALID(oTrace) eq 0 ) then begin
        PRINT,'Error in d_surfviewAddTracePoint: ' + $
            'invalid trace polyline object.'
        status = -1
        RETURN
    endif

    ;  Verify that xyzVector contains 3 elements.
    ;
    sizexyz = SIZE(xyzVector)
    if ((sizexyz[1] ne 3)) then begin
        PRINT,'Error in d_surfviewAddTracePoint:'
        PRINT,'The dimension of xyzVector must be 3.'
        status = -1
        RETURN
    endif

    ;  Get the size of data and the current number of points.
    ;
    oTrace->GetProperty, DATA=data, POLYLINE=connectivityList
    sizeData = size(data)
    nCurrentPoints = connectivityList[0]

    ;  Return if the number of points exceeds the data size.
    ;
    if (nCurrentPoints GT sizeData[2]-2) then begin
        PRINT,'Error in d_surfviewAddTracePoint:'
        PRINT,'The number of data points exceeds the dimension of data'
        status = -1
        RETURN
    endif

    ;  Add the data point.
    ;
    data[0,nCurrentPoints] = xyzVector[0]
    data[1,nCurrentPoints] = xyzVector[1]
    data[2,nCurrentPoints] = xyzVector[2]
    connectivityList[nCurrentPoints+1] = nCurrentPoints
    connectivityList[nCurrentPoints+2] = -1
    connectivityList[0] = nCurrentPoints + 1
    oTrace->SetProperty, DATA=data, POLYLINE=connectivityList

    status = 1

end        ;  of d_surfviewAddDataPoint

;----------------------------------------------------------------------------
;
;  Purpose:  Handle the event.
;
pro d_surfviewEvent, $
    sEvent         ; IN: event structure

    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
    demo_record, sEvent, filename=sState.record_to_filename
    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY

    ;  Quit the application using the close box.
    ;
    if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ $
        'WIDGET_KILL_REQUEST') then begin
        WIDGET_CONTROL, sEvent.top, /DESTROY
        RETURN
    endif

    WIDGET_CONTROL, sEvent.id, GET_UVAL=uval

    case uval of

        ;  Animate the surface accordingly to scenario number 0.
        ;
        'ANIMATE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            d_surfviewAnimateSurface, sState, 0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ; of ANIMATE

        ;  Change the font.
        ;
        'FONT' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            case sEvent.index of
                0: sState.oFont->SetProperty, NAME='Helvetica'
                1: sState.oFont->SetProperty, NAME='Times'
                2: begin
                     if (!VERSION.OS_FAMILY eq "Windows") then $
                         sState.oFont->SetProperty, NAME='Courier New' $
                     else sState.oFont->SetProperty, NAME='Courier'
                end    ; of case 2
            endcase
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ; of FONT

        ;  Set the shading to flat.
        ;
        'FLAT' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oSurface->SetProperty, SHADING=0
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wFlatButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wGouraudButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of FLAT

        ;  Set the shading to gouraud.
        ;
        'GOURAUD' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oSurface->SetProperty, SHADING=1
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wFlatButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wGouraudButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of GOURAUD

        ;  Set the style to point.
        ;
        'POINT' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, SENSITIVE=0
            sState.oSurface->SetProperty, STYLE=0
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of POINT

        ;  Set the style to wire.
        ;
        'WIRE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, SENSITIVE=0
            sState.oSurface->SetProperty, STYLE=1
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of POINT

        ;  Set the style to solid.
        ;
        'SOLID' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, /SENSITIVE

            ;  Make linestyle solid, because it will show when dragging
            ;  at low quality.
            ;
            j = d_surfviewToggleOffOn(sState.wLineStyleButton)
            case j of
                 0: sState.oSurface->SetProperty, LINESTYLE=0 ; solid
                 1: if d_surfviewToggleOffOn(sState.wLineStyleButton) ne 0 $
                    then $
                        PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
                -1: PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
            endcase

            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oSurface->SetProperty, STYLE=2
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SOLID

        ;  Set the style to ruled xz.
        ;
        'RULEDXZ' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, SENSITIVE=0
            sState.oSurface->SetProperty, STYLE=3
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of RULEDXZ

        ;  Set the style to ruled yz.
        ;
        'RULEDYZ' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, SENSITIVE=0
            sState.oSurface->SetProperty, STYLE=4
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of RULEDYZ

        ;  Set the style to lego wire.
        ;
        'LEGOWIRE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, SENSITIVE=0
            sState.oSurface->SetProperty, STYLE=5
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of LEGOWIRE

        ;  Set the style to lego solid.
        ;
        'LEGOSOLID' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

            ;  Make linestyle solid, because it will show when dragging
            ;  at low quality.
            ;
            j = d_surfviewToggleOffOn(sState.wLineStyleButton)
            case j of
                 0: sState.oSurface->SetProperty, LINESTYLE=0 ; solid
                 1: if d_surfviewToggleOffOn(sState.wLineStyleButton) ne 0 $
                     then $
                        PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
                -1: PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
            endcase
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wTracingButton, SENSITIVE=0
            sState.oSurface->SetProperty, STYLE=6
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wPointButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSolidButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledXZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wRuledYZButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoWireButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wLegoSolidButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wHiddenButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wLineStyleButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wShadingButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of LEGOSOLID

        ;  Set maximum or minimum value.
        ;
        'SHOWMAXMIN': begin   ; this handles event from either slider
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            WIDGET_CONTROL, sState.wMinSlider, GET_VALUE=min
            WIDGET_CONTROL, sState.wMaxSlider, GET_VALUE=max

            sState.oSurface->SetProperty, MIN_VALUE=min
            sState.oSurface->SetProperty, MAX_VALUE=max
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SHOWMIN

        ;  Set scaling
        ;
        'SCALING': begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            WIDGET_CONTROL, sState.wScalingSlider, GET_VALUE=scale

            scale = 0.75 + FLOAT(scale) / 100.0
            scalep = scale*100.0
            scalingString = 'Scaling : ' + STRING(scalep, $
                FORMAT='(f5.1)') + ' %'
            WIDGET_CONTROL, sState.wScalingLabel, $
                SET_VALUE=scalingString

            transform = [[scale, 0, 0, 0.0], [0, scale, 0, 0.0], $
                [0, 0, scale, 0.0], [0, 0, 0, 1]]
            sState.oScalingModel->SetProperty, TRANSFORM=transform
            sState.oTraceScalingModel->SetProperty, TRANSFORM=transform
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SCALING

        'CONSTRAINT': begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            case sEvent.index of
                0: begin ; unconstrained rotation
                    sState.oRotationModel->SetProperty, $
                        CONSTRAIN=0
                    sState.oTraceRotationModel->SetProperty, $
                        CONSTRAIN=0
                    end
                1: begin
                    sState.oRotationModel->SetProperty, $
                        AXIS=0, /CONSTRAIN
                    sState.oTraceRotationModel->SetProperty, $
                        AXIS=0, /CONSTRAIN
                    end
                2: begin
                    sState.oRotationModel->SetProperty, $
                        AXIS=1, /CONSTRAIN
                    sState.oTraceRotationModel->SetProperty, $
                        AXIS=1, /CONSTRAIN
                    end
                3: begin
                    sState.oRotationModel->SetProperty, $
                        AXIS=2, /CONSTRAIN
                    sState.oTraceRotationModel->SetProperty, $
                        AXIS=2, /CONSTRAIN
                    end
                4: begin
                    sState.oRotationModel->SetProperty, $
                        AXIS=0, CONSTRAIN=2
                    sState.oTraceRotationModel->SetProperty, $
                        AXIS=0, CONSTRAIN=2
                    end
                5: begin
                    sState.oRotationModel->SetProperty, $
                        AXIS=1, CONSTRAIN=2
                    sState.oTraceRotationModel->SetProperty, $
                        AXIS=1, CONSTRAIN=2
                    end
                6: begin
                    sState.oRotationModel->SetProperty, $
                        AXIS=2, CONSTRAIN=2
                    sState.oTraceRotationModel->SetProperty, $
                        AXIS=2, CONSTRAIN=2
                    end
            endcase

            sState.oRotationModel->GetProperty, CONSTRAIN=constrain
            sState.oRedAxis->GetProperty, HIDE=red_axis_hidden

            if (constrain eq 2) and red_axis_hidden then begin
                WIDGET_CONTROL, sEvent.top, /HOURGLASS
                WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
                sState.oRedAxis->SetProperty, HIDE=0
                sState.oGreenAxis->SetProperty, HIDE=0
                sState.oBlueAxis->SetProperty, HIDE=0
                sState.oWindow->Draw, sState.oView
                WIDGET_CONTROL, sState.wHideAxes, $
                    SET_VALUE='Hide Axes'
                WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            end
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        'HOTKEY' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            case STRUPCASE(sEvent.ch) of
                'U': begin ; unconstrained rotation
                    sState.oRotationModel->SetProperty, CONSTRAIN=0
                    sState.oTraceRotationModel->SetProperty, CONSTRAIN=0
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=0
                    end
                'X': begin
                    sState.oRotationModel->SetProperty, AXIS=0, $
                        /CONSTRAIN
                    sState.oTraceRotationModel->SetProperty, AXIS=0, $
                        /CONSTRAIN
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=1
                    end
                'Y': begin
                    sState.oRotationModel->SetProperty, AXIS=1, $
                        /CONSTRAIN
                    sState.oTraceRotationModel->SetProperty, AXIS=1, $
                        /CONSTRAIN
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=2
                    end
                'Z': begin
                    sState.oRotationModel->SetProperty, AXIS=2, $
                        /CONSTRAIN
                    sState.oTraceRotationModel->SetProperty, AXIS=2, $
                        /CONSTRAIN
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=3
                    end
                'R': begin
                    sState.oRotationModel->SetProperty, AXIS=0, $
                        CONSTRAIN=2
                    sState.oTraceRotationModel->SetProperty, AXIS=0, $
                        CONSTRAIN=2
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=4
                    end
                'G': begin
                    sState.oRotationModel->SetProperty, AXIS=1, $
                        CONSTRAIN=2
                    sState.oTraceRotationModel->SetProperty, AXIS=1, $
                        CONSTRAIN=2
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=5
                    end
                'B': begin
                    sState.oRotationModel->SetProperty, AXIS=2, $
                        CONSTRAIN=2
                    sState.oTraceRotationModel->SetProperty, AXIS=2, $
                        CONSTRAIN=2
                    WIDGET_CONTROL, sState.wConstraintsDroplist, $
                        SET_DROPLIST_SELECT=6
                    end
                'H': begin ; Toggle hide-show rotation axes
                    id = sState.wHideAxes
                    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
                    d_surfviewEvent, {id: id, top: sEvent.top, handler:0l}
                    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
                    end
                else:
            endcase

            sState.oRotationModel->GetProperty, CONSTRAIN=constrain
            sState.oRedAxis->GetProperty, HIDE=red_axis_hidden

            if (constrain eq 2) $
            and red_axis_hidden $
            and STRUPCASE(sEvent.ch) ne 'H' $
            then begin
                WIDGET_CONTROL, sEvent.top, /HOURGLASS
                WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
                sState.oRedAxis->SetProperty, HIDE=0
                sState.oGreenAxis->SetProperty, HIDE=0
                sState.oBlueAxis->SetProperty, HIDE=0
                sState.oWindow->Draw, sState.oView
                WIDGET_CONTROL, sState.wHideAxes, $
                    SET_VALUE='Hide Axes'
                WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            end

            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        ;  Hide the Rotation Axes
        ;
        'HIDEAXES': begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oRedAxis->GetProperty, HIDE=red_axis_hidden
            sState.oRedAxis->SetProperty, HIDE=1-red_axis_hidden
            sState.oGreenAxis->SetProperty, HIDE=1-red_axis_hidden
            sState.oBlueAxis->SetProperty, HIDE=1-red_axis_hidden
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wHideAxes, $
                SET_VALUE=(['Show','Hide'])[red_axis_hidden]+' Axes'
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        ;  Reset the initial orientation of the surface
        ;
        'RESETTRANSFORM': begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oRotationModel->SetProperty, $
                TRANSFORM=sState.initTransformRotation
            sState.oTraceRotationModel->SetProperty, $
                TRANSFORM=sState.initTransformRotation
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of RESETTRANSFORM

        ;  Toggle off and on the texture mapping.
        ;
        'TEXTURE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            j = d_surfviewToggleOffOn(sEvent.id)
            case j of

                0: begin
                    sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
                end    ;   of  0

                1: begin
                    sState.oSurface->SetProperty, $
                        TEXTURE_MAP=sState.oTextureImage

                    ;  Turn off vertex coloring.
                    ;
                    if d_surfviewToggleOffOn(sState.wVertexButton) eq 1 $
                    then $
                        void = d_surfviewToggleOffOn(sState.wVertexButton)
                    sState.oSurface->SetProperty, VERT_COLORS=0

                    if d_surfviewToggleOffOn(sState.wTracingMaskButton) eq 1 $
                    then $
                        void = d_surfviewToggleOffOn( $
                            sState.wTracingMaskButton $
                            )
                    ;  Turn off tracing mode.
                    ;
                    if d_surfviewToggleOffOn(sState.wTracingModeButton) eq 1 $
                    then $
                        void = d_surfviewToggleOffOn( $
                            sState.wTracingModeButton $
                            )
                    sState.tracingMode = 0
                    sState.oTracePolyline->SetProperty, /HIDE

                    if not sState.tracingMode then $
                        WIDGET_CONTROL, sState.wStyleButton, SENSITIVE=1

                    sState.oWindow->Draw, sState.oView
                    demo_putTips, sState, ['infor','instr'], [11,12], /LABEL
                end    ;   of  1

                -1: PRINT, 'Error in d_surfviewToggleOffOn/texture.'

            endcase
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of TEXTURE

        ;  Toggle off and on the vertex colors.
        ;
        'VERTEXCOLOR' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            j = d_surfviewToggleOffOn(sEvent.id)
            case j of
                 0: sState.oSurface->SetProperty, VERT_COLORS=0
                 1: begin

                    ;  Turn off texture mapping.
                    ;
                    if d_surfviewToggleOffOn(sState.wTextureButton) eq 1 $
                    then $
                        void = d_surfviewToggleOffOn(sState.wTextureButton) $
                    else begin
                        sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
                        sState.oSurface->GetProperty, STYLE=style
                        WIDGET_CONTROL, sState.wTracingButton, $
                            SENSITIVE=([0,1])[style eq 2];Solid can be traced
                    end

                    ;  Turn on vertex colors.
                    ;
                    sState.oSurface->SetProperty, $
                        VERT_COLORS=sState.vertexColors
                end
                -1: PRINT, 'Error in d_surfviewToggleOffOn/vertexcolors.'
            endcase
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of VERTEXCOLOR

        ;  Toggle off and on the hidden points and lines.
        ;
        'HIDE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            j = d_surfviewToggleOffOn(sEvent.id)
            if (j eq -1) then begin
                 PRINT, 'Error in d_surfviewToggleOffOn/hide.'
                 RETURN
            endif
            sState.oSurface->SetProperty, HIDDEN_LINES=j
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of HIDE

        ;  Toggle between solid and dash linestyles.
        ;
        'LINESTYLE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            j = d_surfviewToggleOffOn(sEvent.id)
            case j of
                 0: sState.oSurface->SetProperty, LINESTYLE=0 ; solid
                 1: sState.oSurface->SetProperty, LINESTYLE=4
                -1: PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
            endcase
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of LINESTYLE

        ;  Show no skirt.
        ;
        'SKIRTNONE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oSurface->SetProperty, SHOW_SKIRT = 0
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wSkirtNoneButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wSkirt10Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt20Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt30Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SKIRTNONE

        ;  Set skirt to -0.5.
        ;
        'SKIRT10' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oSurface->GetProperty, ZRANGE=zrange
            sState.oSurface->SetProperty, SKIRT=zrange[0], /SHOW_SKIRT
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wSkirtNoneButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt10Button, SENSITIVE=0
            WIDGET_CONTROL, sState.wSkirt20Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt30Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SKIR10

        ;  Set skirt to 0.0.
        ;
        'SKIRT20' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oSurface->GetProperty, ZRANGE=zrange
            sState.oSurface->SetProperty, $
                SKIRT=(zrange[1] - zrange[0]) / 2. + zrange[0], $
                /SHOW_SKIRT
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wSkirtNoneButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt10Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt20Button, SENSITIVE=0
            WIDGET_CONTROL, sState.wSkirt30Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SKIR20

        ;  Set skirt to 0.5.
        ;
        'SKIRT30' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            sState.oSurface->GetProperty, ZRANGE=zrange
            sState.oSurface->SetProperty, SKIRT=zrange[1], /SHOW_SKIRT
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sState.wSkirtNoneButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt10Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt20Button, SENSITIVE=1
            WIDGET_CONTROL, sState.wSkirt30Button, SENSITIVE=0
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of SKIR30

        ;  Set drag quality to low.
        ;
        'LOW' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            sState.dragq = 0
            WIDGET_CONTROL, sState.wLowButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wMediumButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHighButton, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of LOW

        ;  Set drag quality to medium.
        ;
        'MEDIUM' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            sState.dragq = 1
            WIDGET_CONTROL, sState.wLowButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wMediumButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wHighButton, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of MEDIUM

        ;  Set drag quality to high.
        ;
        'HIGH' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            sState.dragq = 2
            WIDGET_CONTROL, sState.wLowButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wMediumButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wHighButton, SENSITIVE=0
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end        ;  of HIGH

        ;  Draw only the surface within the tracing
        ;  contour (trace mask on),
        ;  or draw the whole surface (trace mask off).
        ;  Toggle between these 2 options.
        ;
        'TRACING_MASK' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /HOURGLASS

            case d_surfviewToggleOffOn(sEvent.id) of
                0: begin
                    sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
                    if sState.tracingMode eq 0 then begin
                        WIDGET_CONTROL, sState.wStyleButton, SENSITIVE=1
                        WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=1
                        demo_putTips, $
                            sState, $
                            ['infor','instr'], $
                            [11,12], $
                            /LABEL
                    end
                end
                1: begin
                    sState.oTracePolyline->GetProperty, $
                        DATA=data, POLYLINE=connectivityList
                    sState.oTracingMask->GetProperty, DATA=idata
                    siz = SIZE(idata)
                    mask = BYTARR(siz[1], siz[2])
                    idata = BYTARR(siz[1], siz[2], 4) + 255b
                    if (connectivityList[0] ge 3) then begin
                        ;  Get the number of current points minus 1
                        ;
                        nCurrentPointsM1 = connectivityList[0] - 1
                        x = data[0, 0:nCurrentPointsM1] * 10.0
                        y = data[1, 0:nCurrentPointsM1] * 10.0
                        fill = POLYFILLV(x, y, siz[1], siz[2])
                        mask[*, *] = 0
                        if fill[0] gt -1 then $
                            mask[fill] = 255
                    endif else begin
                        mask[*, *] = 255
                    endelse
                    idata[*, *, 3] = mask
                    sState.oTracingMask->SetProperty, DATA=idata
                    sState.oSurface->SetProperty, $
                        TEXTURE_MAP=sState.oTracingMask

                    if d_surfviewToggleOffOn(sState.wTextureButton) eq 1 $
                    then $
                        void = d_surfviewToggleOffOn(sState.wTextureButton)

                    WIDGET_CONTROL, sEvent.top, /HOURGLASS

                    sState.oWindow->Draw, sState.oView
                    WIDGET_CONTROL, sState.wStyleButton, SENSITIVE=0
                    if sState.tracingMode eq 1 then begin
                        demo_putTips, $
                            sState, $
                            ['reset','erase'], $
                            [11,12], $
                            /LABEL
                    end else begin
                        demo_putTips, $
                            sState, $
                            ['mask2','disp2'], $
                            [11,12], $
                            /LABEL
                    end
                end          ; of 1
                -1: PRINT, 'Error in d_surfviewToggleOffOn/tracingMask.'
            endcase

            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end         ; of TRACING_MASK

        ;  Enable (on) or disable (off) the tracing mode.
        ;  Toggle between these 2 options.
        ;
        'TRACING_MODE' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            j = d_surfviewToggleOffOn(sEvent.id)
            case j of
                -1: begin
                     PRINT, 'Error in d_surfviewToggleOffOn/tracingMode.'
                     RETURN
                end
                1: begin
                    sState.tracingMode = 1
                    sState.oTracePolyline->GetProperty, $
                        POLYLINE=connectivityList
                    connectivityList[0] = 0
                    connectivityList[1] = -1
                    sState.oTracePolyline->SetProperty, $
                        POLYLINE=connectivityList, $
                        HIDE=0
                    if d_surfviewToggleOffOn(sState.wTextureButton) eq 1 $
                    then $
                        void = d_surfviewToggleOffOn(sState.wTextureButton)
                    sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
                    sState.oWindow->Draw, sState.oView

                    WIDGET_CONTROL, sState.wTracingMaskButton, SENSITIVE=0

                    WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=0
                    WIDGET_CONTROL, SsTATE.wAnimateButton, SENSITIVE=0
                    WIDGET_CONTROL, sState.wStyleButton, SENSITIVE=0
                    demo_putTips, sState, ['regio','right'], [11,12], /LABEL
                end
                0: begin
                    sState.tracingMode = 0
                    sState.oTracePolyline->SetProperty, /HIDE
                    sState.oWindow->Draw, sState.oView
                    sState.oSurface->GetProperty, TEXTURE_MAP=texture_map
                    if OBJ_VALID(texture_map) then begin
                        demo_putTips, $
                            sState, $
                            ['mask2','disp2'], $
                            [11,12], $
                            /LABEL
                    endif else begin
                        WIDGET_CONTROL, sState.wStyleButton, SENSITIVE=1
                        demo_putTips, $
                            sState, $
                            ['infor','instr'], $
                            [11,12], $
                            /LABEL
                    end
                    WIDGET_CONTROL, sState.wTextureButton, SENSITIVE=1
                    WIDGET_CONTROL, sState.wAnimateButton, SENSITIVE=1

                end
            endcase
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end         ; of TRACING_MODE

        ;  Handle the event that occurs in the drawing (viewing) area.
        ;  These are : expose and mouse button(press, motion, release).
        ;
        'DRAW': begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

            ;  Expose.
            ;
            if (sEvent.type eq 4) then begin
                sState.oWindow->draw, sState.oView
                WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
                RETURN
            endif

            ;  User rotation.
            ;
            void = sState.oTraceRotationModel->Update(sEvent)
            if sState.oRotationModel->Update(sEvent) then $
                sState.oWindow->Draw, sState.oView

            ;  Button press.
            ;
            if (sEvent.type EQ 0) then begin

                ;  Right button press.
                ;
                if (sEvent.press EQ 4) then begin
                    pick = sState.oWindow->PickData(sState.oView, $
                    sState.oSurface, [sEvent.x, sEvent.y], dataxyz)
                    if (pick ne 0) then begin
                         statusString = STRING(dataxyz[0], $
                             dataxyz[1],dataxyz[2], $
                             FORMAT='("X=", F6.2,' + $
                             ' ", Y=",F6.2,", Z=",F6.2)')
                        demo_putTips, sState, 'locat', 11, /LABEL
                        demo_putTips, sState, statusString, 12
                        if (sState.tracingMode eq 1) then begin
                            d_surfviewAddTracePoint, sState.oTracePolyline, $
                                dataxyz, check
                            sState.firstPoint = dataxyz
                            if (check eq -1) then begin
                                PRINT, $
                                    'Error in d_surfviewAddTracePoint/Draw.'
                            endif
                            sState.oWindow->Draw, sState.oView
                        endif
                        sState.btndown = 4b
                        WIDGET_CONTROL, sState.wDraw, /DRAW_MOTION
                    endif else begin
                        demo_putTips, sState, "Data point:In background", 12
                    endelse
                endif else begin
                    sState.btndown = 1b
                    sState.oWindow->SetProperty, QUALITY=sState.dragq
                    WIDGET_CONTROL, sState.wDraw, /DRAW_MOTION
                endelse
            endif    ; ev.typ EQ 0

            ;  Button motion and tracing mode
            ;
            if ((sEvent.type eq 2) and (sState.btndown eq 4b)) then begin
                pick = sState.oWindow->PickData(sState.oView, $
                sState.oSurface, [sEvent.x, sEvent.y], dataxyz)
                if (pick ne 0) then begin
                     statusString = STRING(dataxyz[0], $
                         dataxyz[1],dataxyz[2], $
                         FORMAT='("X=", F6.2,' + $
                         ' ", Y=",F6.2,", Z=",F6.2)')
                        demo_putTips, sState, 'locat', 11, /LABEL
                        demo_putTips, sState, statusString, 12
                    if (sState.tracingMode ne 0) then begin
                        d_surfviewAddTracePoint, sState.oTracePolyline, $
                            dataxyz, check
                        if (check eq -1) then begin
                            PRINT,'Error in d_surfviewAddTracePoint/Draw.'
                        endif
                        sState.oWindow->Draw, sState.oView
                    endif
                endif else begin
                    demo_putTips, sState, "Data point:In background", 12
                endelse
            endif

            ;  Button release.
            ;
            if (sEvent.type eq 1) then begin
                if (sState.btndown EQ 1b) then begin
                        sState.oWindow->SetProperty, QUALITY=2 ; High
                        if sState.dragq ne 2 then begin
                            sState.oSurface->GetProperty, STYLE=style
                            if style eq 5 then $ ; Lego
                                widget_control, /hourglass
                            if style eq 6 then $ ; Filled Lego
                                widget_control, /hourglass
                            sState.oWindow->Draw, sState.oView
                        endif
                endif else if ((sState.btndown EQ 4b) $
                    AND (sState.tracingmode EQ 1) ) then begin
                    d_surfviewAddTracePoint, sState.oTracePolyline, $
                        sState.firstPoint, check
                    if (check eq -1) then begin
                        PRINT,'Error in d_surfviewAddTracePoint/Draw.'
                    endif
                    sState.oWindow->Draw, sState.oView
                    sState.tracingMode = 1

                    sState.oTracePolyline->GetProperty, POLYLINE=polyline
                    if polyline[0] ge 3 then begin
                        WIDGET_CONTROL, sState.wTracingMaskButton, SENSITIVE=1
                        demo_putTips, sState, $
                            ['mask1','displ'], [11,12], /LABEL
                        if d_surfviewToggleOffOn(sState.wTracingMaskButton) $
                        eq 1 $
                        then $
                            void = d_surfviewToggleOffOn( $
                                sState.wTracingMaskButton $
                                )
                    end

                endif
                sState.btndown = 0b
                WIDGET_CONTROL, sState.wDraw, DRAW_MOTION = 0
            endif
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end   ; of DRAW

        ;  Quit the application.
        ;
        'QUIT' : begin
           WIDGET_CONTROL, sEvent.top, /DESTROY
           RETURN
        end   ; of QUIT

    endcase   ; of case uval of

    if XREGISTERED('demo_tour') eq 0 then begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
        WIDGET_CONTROL, sState.wHotKeyReceptor, /INPUT_FOCUS
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
    end
end     ; of event handler

;----------------------------------------------------------------------------
;
;  Purpose:  Restore the previous color table and
;            destroy the top objects.
;
pro d_surfviewCleanup, $
    wTopBase          ;  IN: top level base ID.
    WIDGET_CONTROL, wTopBase, GET_UVALUE=sState, /NO_COPY

    ;Clean up heap variables.
    ;
    for i=0,n_tags(sState)-1 do begin
        case size((sState).(i), /TNAME) of
            'POINTER': $
                ptr_free, (sState).(i)
            'OBJREF': $
                obj_destroy, (sState).(i)
            else:
        endcase
    end

    ;  Silently flush any accumulated math errors.
    ;
    void = check_math()

    ;  Restore math error behavior.
    ;
    !except = sState.orig_except

    ;  Restore the color table
    ;
    TVLCT, sState.colorTable

    ;  Map the group leader base if it exists.
    ;
    if (WIDGET_INFO(sState.groupBase, /VALID_ID)) then $
        WIDGET_CONTROL, sState.groupBase, /MAP

end   ;  of d_surfviewCleanup



;----------------------------------------------------------------------------
;
;  Purpose:  Display a surface.
;
pro d_surfview, $
;    ALT_FUNC=ALT_FUNC, $       ; IN: (opt) Alternative function : sine dist
    TRANSPARENT=transparent, $ ; IN: (opt) Transparent across a plane
    GROUP=group, $             ; IN: (opt) group identifier
    DEBUG=debug, $             ; IN: (opt) debug mode
    RECORD_TO_FILENAME=record_to_filename, $
    APPTLB = appTLB            ; OUT: (opt) TLB of this application

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

    ;  Set up dimensions of the drawing (viewing) area.
    ;
    device, GET_SCREEN_SIZE=scr
    xdim = scr[0]*0.6
    ydim = xdim*0.8

    ;  Get the current color vectors to restore
    ;  when this application is exited.
    TVLCT, savedR, savedG, savedB, /GET

    ; Build color table from color vectors
    ;
    colorTable = [[savedR],[savedG],[savedB]]

    ;  Load an initial file .
    ;
    file = 'abnorm.dat'
    demo_getdata, NewImage, FILENAME=file, /TWO_DIM
    z=NewImage[*,*,7]
    z=smooth(z, 3, /EDGE_TRUNCATE)
    siz = SIZE(z)
    z=REBIN(z, siz[1]/2, siz[2]/2)

    ;  Create widgets.
    ;
    if (N_ELEMENTS(group) EQ 0) then begin
        wTopBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            TLB_FRAME_ATTR=1, MBAR=barBase, $
            UNAME='d_surfview:tlb', $
            TITLE="Surface Objects")
    endif else begin
        wTopBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            TLB_FRAME_ATTR=1, MBAR=barBase, $
            GROUP_LEADER=group, $
            UNAME='d_surfview:tlb', $
            TITLE="Surface Objects")
    endelse

        ;  Create the menu bar. It contains the file,
        ;  edit, and help buttons.
        ;
        wFileButton = WIDGET_BUTTON(barBase, VALUE = 'File', /MENU)
            wQuitButton = WIDGET_BUTTON(wFileButton, VALUE='Quit', $
                UNAME='d_surfview:quit', $
                UVALUE='QUIT')

        ;  Create the menu bar item Edit
        ;  that has the shade and style options.
        ;
        wOptionButton = WIDGET_BUTTON(barBase, VALUE = 'Options', /MENU)

            ;  Select the plot shading.
            ;
            wShadingButton = WIDGET_BUTTON(wOptionButton, $
                VALUE='Shading', UVALUE='SHADING', MENU=1)

                wFlatButton = WIDGET_BUTTON(wShadingButton, $
                   VALUE='Flat', UVALUE='FLAT', UNAME='d_surfview:flat')

                wGouraudButton = WIDGET_BUTTON(wShadingButton, $
                   VALUE='Gouraud', UVALUE='GOURAUD', $
                   UNAME='d_surfview:gouraud')

            ;  Select the plot style.
            ;
            wStyleButton = WIDGET_BUTTON(wOptionButton, $
                VALUE='Style', UVALUE='STYLE', /MENU)

                wPointButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Point', UVALUE='POINT', UNAME='d_surfview:point')

                wWireButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Wire', UVALUE='WIRE', UNAME='d_surfview:wire')

                wSolidButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Solid', UVALUE='SOLID', UNAME='d_surfview:solid')

                wRuledXZButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Ruled XZ', UVALUE='RULEDXZ', $
                        UNAME='d_surfview:ruledxz')

                wRuledYZButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Ruled YZ', UVALUE='RULEDYZ', $
                        UNAME='d_surfview:ruledyz')

                wLegoWireButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Lego Wire', UVALUE='LEGOWIRE', $
                    UNAME='d_surfview:legowire')

                wLegoSolidButton = WIDGET_BUTTON(wStyleButton, $
                    VALUE='Lego Solid', UVALUE='LEGOSOLID', $
                    UNAME='d_surfview:legosolid')

            ;  Select the skirt value.
            ;
            wSkirtButton = WIDGET_BUTTON(wOptionButton, $
                VALUE='Skirt', UVALUE='SKIRT', /MENU)

                ; Remove the skirt
                ;
                wSkirtNoneButton = WIDGET_BUTTON(wSkirtButton, $
                    VALUE='None', UVALUE='SKIRTNONE', $
                    UNAME='d_surfview:skirtnone')

                ; Skirt10 is  -0.5.
                ;
                wSkirt10Button = WIDGET_BUTTON(wSkirtButton, $
                    VALUE='z=6', UVALUE='SKIRT10', $
                    UNAME='d_surfview:skirt6')

                ; Skirt20 is  0.0
                ;
                wSkirt20Button = WIDGET_BUTTON(wSkirtButton, $
                    VALUE='z=102', UVALUE='SKIRT20', $
                    UNAME='d_surfview:skirt102')

                ; Skirt30 is  0.5
                ;
                wSkirt30Button = WIDGET_BUTTON(wSkirtButton, $
                    VALUE='z=198', UVALUE='SKIRT30', $
                    UNAME='d_surfview:skirt198')

            ;  Set up the drag quality. Low is wire, medium is
            ;  polygons, high is smoothed polygons.
            ;
            wDragButton = Widget_Button(wOptionButton, $
                VALUE="Drag Quality", UVALUE='DRAGQ', /MENU)

                wLowButton = WIDGET_BUTTON(wDragButton, $
                    VALUE='Low', UVALUE='LOW', UNAME='d_surfview:lowdrag')

                wMediumButton = WIDGET_BUTTON(wDragButton, $
                    VALUE='Medium', UVALUE='MEDIUM', $
                    UNAME='d_surfview:meddrag')

                wHighButton = WIDGET_BUTTON(wDragButton, $
                    VALUE='High', UVALUE='HIGH', $
                    UNAME='d_surfview:highdrag')

            ;  Allows to trace a contour on the surface, and
            ;  then to display only the surface contained within
            ;  that contour.
            ;
            wTracingButton = WIDGET_BUTTON(wOptionButton, $
                VALUE='Tracing', UVALUE='TRACING', /MENU)

                wTracingModeButton = WIDGET_BUTTON(wTracingButton, $
                    VALUE='Activate Outlining', $
                    UVALUE='TRACING_MODE', $
                    UNAME='d_surfview:tracing_mode')

                wTracingMaskButton = WIDGET_BUTTON(wTracingButton, $
                    VALUE='Activate Trace Mask', UVALUE='TRACING_MASK', $
                    UNAME='d_surfview:tracing_mask')

            ;  Toggle between showing or not showing
            ;  the hidden points and lines.
            ;
            wHiddenButton = WIDGET_BUTTON(wOptionButton, $
                VALUE="Activate Hiding", UVALUE='HIDE', $
                UNAME='d_surfview:hidden')

            ;  Toggle between showing or not showing
            ;  the vertices in colors.
            ;
            wVertexButton = WIDGET_BUTTON(wOptionButton, $
                VALUE="Activate Vertex Colors", UVALUE='VERTEXCOLOR', $
                UNAME='d_surfview:vertcolor')

            ;  Toggle between showing or not showing
            ;  the texture mapping.
            ;
            wTextureButton = WIDGET_BUTTON(wOptionButton, $
                VALUE="Deactivate Texture Map", UVALUE='TEXTURE', $
                UNAME='d_surfview:texture')

            ;  Toggle between a solid or dash line style.
            ;
            wLineStyleButton = WIDGET_BUTTON(wOptionButton, $
                VALUE="Activate Line Style", UVALUE='LINESTYLE', $
                UNAME='d_surfview:linestyle')

        ;  Create the menu bar item Edit
        ;  that has the shade and style options.
        ;
        wViewButton = WIDGET_BUTTON(barBase, VALUE = 'View', /MENU)

            wAnimateButton = WIDGET_BUTTON(wViewButton, $
                VALUE="Animate", UVALUE='ANIMATE', $
                UNAME='d_surfview:animate')

            wResetButton = WIDGET_BUTTON(wViewButton, $
                VALUE="Reset Orientation", UVALUE='RESETTRANSFORM', $
                UNAME='d_surfview:reset_transfrm')


        ;  Create a sub base of the top base (wBase).
        ;
        wSubBase = WIDGET_BASE(wTopBase, COLUMN=2)

            ;  Create the left Base that contains the functionality buttons.
            ;  Here the only button is to animate the object.
            ;
            wLeftbase = WIDGET_BASE(wSubBase, $ ;/BASE_ALIGN_CENTER, $
                COLUMN=1)

                    minValue = MIN(z, MAX=maxValue)
                    wMaxLabel = WIDGET_LABEL(wLeftBase, $
                        VALUE='Data Maximum:')

                    wMaxSlider = WIDGET_SLIDER(wLeftBase, $
                        MINIMUM=minValue, $
                        MAXIMUM=maxValue, VALUE=maxValue, $
                        UVALUE='SHOWMAXMIN', $
                        UNAME='d_surfview:max_slider')

                    wMinLabel = WIDGET_LABEL(wLeftBase, $
                        VALUE='Data Minimum:')

                    wMinSlider = WIDGET_SLIDER(wLeftBase, $
                        MINIMUM=minValue, $
                        MAXIMUM=maxValue, VALUE=minValue, $
                        UVALUE='SHOWMAXMIN', $
                        UNAME='d_surfview:min_slider')

                wScalingBase = WIDGET_BASE(wLeftBase, $
                    /COLUMN, YPAD=10)

                    percent = 100
                    scalingString = 'Scaling : ' + STRING(percent, $
                        FORMAT='(f5.1)') + ' %'
                    wScalingLabel = WIDGET_LABEL(wScalingBase, $
                        VALUE=scalingString)

                    wScalingSlider = WIDGET_SLIDER(wScalingBase, $
                        MINIMUM=0, $
                        MAXIMUM=50, VALUE=25, $
                        /SUPPRESS_VALUE, $
                        UVALUE='SCALING', $
                        UNAME='d_surfview:zoom')

                void = WIDGET_LABEL(wLeftBase, value='Constrain Rotations:')
                wConstraintsDroplist = WIDGET_DROPLIST(wLeftBase, $
                    VALUE=['Unconstrained', $
                        "about Screen X", $
                        "about Screen Y", $
                        "about Screen Z", $
                        "about Data X (Red)", $
                        "about Data Y (Green)", $
                        "about Data Z (Blue)" $
                        ], $
                    UVALUE='CONSTRAINT', $
                    UNAME='d_surfview:constraint')

                wHideAxes = WIDGET_BUTTON(wLeftBase, $
                    value='Show Axes', $
                    UVALUE='HIDEAXES', $
                    UNAME='d_surfview:hide_axes')

            ;  Create the right Base that has the drawing area.
            ;
            wRightbase = WIDGET_BASE(wSubBase)

                wDraw = WIDGET_DRAW(wRightBase, $
                    GRAPHICS_LEVEL=2, $
                    XSIZE=xdim, YSIZE=ydim, /BUTTON_EVENTS, $
                    UVALUE='DRAW', RETAIN=0, /EXPOSE_EVENT, $
                    UNAME='d_surfview:draw')
                wHotKeyReceptor = WIDGET_TEXT(wRightBase, $
                    /ALL_EVENTS, $
                    UVALUE='HOTKEY', $
                    UNAME='d_surfview:hotkey')

        ;  Create tips texts.
        ;
        wStatusBase = WIDGET_BASE(wTopBase, MAP=0, /ROW)

    ;  Now the widget have been created, realize it.
    ;
    WIDGET_CONTROL, wTopBase, /REALIZE

    ; Returns the top level base to the APPTLB keyword.
    ;
    appTLB = wtopBase

    ;  Get the tips
    ;
    sText = demo_getTips(demo_filepath('surfview.tip', $
                             SUBDIR=['examples','demo', 'demotext']), $
                         wTopBase, $
                         wStatusBase)

    ;  Grab the window id of the drawable.
    ;
    WIDGET_CONTROL, wDraw, GET_VALUE=oWindow

    WIDGET_CONTROL, wTopBase, SENSITIVE=0


    bias = -0.5
    ;  Compute viewplane rectangle based on aspect ratio.
    ;
    aspect = float(xdim)/float(ydim)
    if (aspect > 1) then $
        myview = [(1.0-aspect)/2.0+bias, 0.0+bias, aspect, 1.0] $
    else $
        myview = [0.0+bias, (1.0-(1.0/aspect))/2.0+bias, 1.0, (1.0/aspect)]

    ;  Create view object.
    ;
    oView = OBJ_NEW('IDLgrView', PROJECTION=1, EYE=3, $
        ZCLIP=[1.5,-1.5], VIEW=myview, COLOR=[0, 0, 0])

    ;  Make the text location to be centered.
    ;
    textLocation = [myview[0]+0.5*myview[2], myview[1]+0.5*myview[3]]

    ;  Create and display the PLEASE WAIT text.
    ;
    oFont = OBJ_NEW('IDLgrFont', 'Helvetica', SIZE=18)
    oText = OBJ_NEW('IDLgrText', $
        'Starting up  Please wait...', $
        ALIGN=0.5, $
        LOCATION=textLocation, $
        COLOR=[255,255,0], FONT=oFont)

    ;  Create model.
    ;
    oTopModel = OBJ_NEW('IDLgrModel')
    oScalingModel = OBJ_NEW('IDLgrModel')
    oRotationModel = OBJ_NEW('IDLexRotator', $
        [xdim/2.0, ydim/2.0], ydim/2.0)

    oTopModel->Add, oScalingModel
    oScalingModel->Add, oRotationModel

    ;  To avoid conflict with the surface, trace outlines will be
    ;  offset slightly.
    ;
    oTraceScalingModel = OBJ_NEW('IDLgrModel')
    oTraceRotationModel = OBJ_NEW('IDLexRotator', $
        [xdim/2.0, ydim/2.0], ydim/2.0)
    oTraceOffset = OBJ_NEW('IDLgrModel')

    oTopModel->Add, oTraceOffset
    oTraceOffset->Add, oTraceScalingModel
    oTraceScalingModel->Add, oTraceRotationModel

    oTraceOffset->Translate, 0, 0, .01 ; Value determined empirically.

    ;  Place the model in the view.
    ;
    oView->Add, oTopModel

    ;  Scale the top model to fit the viewing area.
    ;
    sct = 0.6
    oTopModel->Scale, sct, sct, sct

    oTopModel->Add, oText

    ;  Draw the starting up screen.
    ;
    oWindow->Draw, oView

    ;  Compute coordinate conversion to normalize.
    ;
    sz = SIZE(z)
    maxx = sz[1] - 1
    maxy = sz[2] - 1
    maxz = MAX(z, MIN=minz)
    xs = [0+bias, 1.0/maxx]
    ys = [0+bias, 1.0/maxy]
    minz2 = minz - 1
    maxz2 = maxz + 1
    zs = [-minz2/(maxz2-minz2)+bias, 1.0/(maxz2-minz2)]

    ;  For height-fields, use the following vertex colors.
    ;
    vertexColors = BYTARR(3, sz[1]*sz[2], /NOZERO)
    cbins= $
        [[255,   0, 0],$
         [255,  85, 0],$
         [255, 170, 0],$
         [255, 255, 0],$
         [170, 255, 0],$
         [85,  255, 0],$
         [0,   255, 0]]

    zi = ROUND(z/float(maxz) * 6.0)
    vertexColors[*, *] = cbins[*, zi]

    oPalette = OBJ_NEW('IDLgrPalette')
    oPalette->LOADCT, 13

    oTextureImage = OBJ_NEW('IDLgrImage', $
        BYTSCL(REBIN(z,256,256)), PALETTE=oPalette)

    ;  Create the tracing objects.
    ;
    workData = BYTARR((size(z))[1]*10, (size(z))[2]*10, 4) + 255b
    oTracingMask = OBJ_NEW('IDLgrImage', workData, INTERLEAVE=2)
    tracingData = FLTARR(3, 1024)
    tracingConnectivityList = LONARR(1024)
    tracingConnectivityList[0] = 0
    tracingConnectivityList[1] = -1

    ;  Create polyline object in the same space as the surface.
    ;
    oTracePolyline = OBJ_NEW('IDLgrPolyline', $
        tracingData, POLYLINES=tracingConnectivityList, $
        COLOR=[255, 0, 0], $
        XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs, $
        THICK=3)
    oTraceRotationModel->Add, oTracePolyline

    ;  Create surface object.
    ;
    oSurface = OBJ_NEW('IDLgrSurface', $
        z, STYLE=2, $
        SHADING=1, $
        /USE_TRIANGLES, $
        TEXTURE_MAP=oTextureImage, $
        COLOR=[230, 230, 230], BOTTOM=[64, 192, 128], $
        XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)

    oRotationModel->Add, oSurface

    ;  Create axes objects
    ;
    stick = [-.5, .5] * 1.5
    oRedAxis = OBJ_NEW('IDLgrPolyline', $
        stick, [0,0], [0,0], COLOR=[255,0,0], /HIDE)
    oGreenAxis = OBJ_NEW('IDLgrPolyline', $
        [0,0], stick, [0,0], COLOR=[0,255,0], /HIDE)
    oBlueAxis = OBJ_NEW('IDLgrPolyline', $
        [0,0], [0,0], stick, COLOR=[0,0,255], /HIDE)
    oRotationModel->Add, oRedAxis
    oRotationModel->Add, oGreenAxis
    oRotationModel->Add, oBlueAxis

    ;  Create a light.
    ;
    oSunLight = OBJ_NEW('IDLgrLight', LOCATION=[1.5, 0, 1], $
        TYPE=1, INTENSITY=0.5)
    otopModel->Add, oSunLight
    oSunLight = OBJ_NEW('IDLgrLight', TYPE=0, $
        INTENSITY=0.75)
    otopModel->Add, oSunLight

    ; Rotate to standard view for first draw.
    ;
    oRotationModel->Rotate,      [1, 0, 0], -90
    oTraceRotationModel->Rotate, [1, 0, 0], -90
    oRotationModel->Rotate, [0, 1, 0], 30
    oRotationModel->Rotate, [1, 0, 0], 30
    oTraceRotationModel->Rotate, [0, 1, 0], 30
    oTraceRotationModel->Rotate, [1, 0, 0], 30

    oRotationModel->GetProperty, TRANSFORM=initTransformRotation
    oScalingModel->GetProperty, TRANSFORM=initTransformScaling

    if not keyword_set(record_to_filename) then $
        record_to_filename = ''

    ;  Create the sState
    ;
    sState = { $
        colorTable: colorTable, $
        center: xdim/2., $                           ; Center of view
        radius: ydim/2, $                            ; Radius
        xSizeData: sz[1], $                          ; x data dimension
        ySizeData: sz[2], $                          ; x data dimension
        btndown: 0b, $                               ; Botton down flag
        pt0: FLTARR(3), $                            ; Initial point
        pt1: FLTARR(3), $                            ; Final point
        dragq: 2, $                                  ; Drag quality
        firstPoint: FLTARR(3), $
        wDraw: wDraw, $                              ; Widget draw
        oTopModel: oTopModel, $                      ; Top model
        oScalingModel: oScalingModel, $
        oTraceScalingModel: oTraceScalingModel, $
        oRotationModel: oRotationModel, $
        oTraceRotationModel: oTraceRotationModel, $
        oSurface: oSurface, $                        ; Surface object
        oRedAxis: oRedAxis, $
        oGreenAxis: oGreenAxis, $
        oBlueAxis: oBlueAxis, $
        initTransformScaling: initTransformScaling, $
        initTransformRotation: initTransformRotation, $
        oView: oView, $                              ; Main view object
        oFont: oFont, $                              ; Font object
        oText: oText, $                              ; Text object
        vertexColors: vertexColors, $                ; Vertex colors(RGB)
        oTracePolyline: oTracePolyline, $            ; Trace object
        oTracingMask: oTracingMask, $                ; Tracing mask object
        oTextureImage: oTextureImage, $              ; Texture image object
        oPalette: oPalette, $                        ; Palette for image
        tracingMode: 0, $                            ; 0=off,1=on
        oWindow: oWindow, $                          ; Window object
        wTopBase : wTopBase, $                       ; Top level base
        wShadingButton : wShadingButton, $
        wFlatButton : wFlatButton, $                 ; Shading options
        wGouraudButton : wGouraudButton, $
        wStyleButton: wStyleButton, $
        wPointButton : wPointButton, $             ; Styles options
        wWireButton : wWireButton, $
        wSolidButton : wSolidButton, $
        wRuledXZButton : wRuledXZButton, $
        wRuledYZButton : wRuledYZButton, $
        wLegoWireButton : wLegoWireButton, $
        wLegoSolidButton : wLegoSolidButton, $
        wSkirtNoneButton : wSkirtNoneButton, $       ; Skirt options
        wSkirt10Button : wSkirt10Button, $
        wSkirt20Button : wSkirt20Button, $
        wSkirt30Button : wSkirt30Button, $
        wLowButton : wLowButton, $                   ; Drag quality options
        wMediumButton : wMediumButton, $
        wHighButton : wHighButton, $
        wMinSlider: wMinSlider, $                    ; Sliders IDs
        wMaxSlider: wMaxSlider, $
        wMinLabel: wMinLabel, $                      ; Sliders label IDs
        wMaxLabel: wMaxLabel, $
        wScalingSlider: wScalingSlider, $
        wScalingLabel: wScalingLabel, $
        wTracingButton : wTracingButton, $           ; Tracing options
        wTracingModeButton : wTracingModeButton, $
        wTracingMaskButton : wTracingMaskButton, $
        wHiddenButton : wHiddenButton, $             ; Hidden option
        wHideAxes: wHideAxes, $
        wAnimateButton: wAnimateButton, $
        wHotKeyReceptor: wHotKeyReceptor, $
        wConstraintsDroplist: wConstraintsDroplist, $
        wTextureButton : wTextureButton, $           ; Texture mapping option
        wVertexButton: wVertexButton, $              ; Vertex coloring option
        sText: sText, $                              ; Text tips structure
        wLineStyleButton : wLineStyleButton, $       ; Linestyle option
        orig_except : !except, $
        record_to_filename: record_to_filename, $
        debug: keyword_set(debug), $
        groupBase: groupBase $                       ; Base of Group Leader
    }

    ;  Unless we are debugging, accumulate any math warnings silently.
    ;
    !except = ([0, 2])[keyword_set(debug)]

    ;  Desensitize the defaults buttons.
    ;
    WIDGET_CONTROL, wGouraudButton, SENSITIVE=0
    WIDGET_CONTROL, wSolidButton, SENSITIVE=0
    WIDGET_CONTROL, wSkirtNoneButton, SENSITIVE=0
    WIDGET_CONTROL, wHighButton, SENSITIVE=0
    WIDGET_CONTROL, wHiddenButton, SENSITIVE=0
    WIDGET_CONTROL, wLineStyleButton, SENSITIVE=0
    WIDGET_CONTROL, wTracingMaskButton, SENSITIVE=0

    WIDGET_CONTROL, wTopBase, SET_UVALUE=sState, /NO_COPY

    WIDGET_CONTROL, wTopBase, /HOURGLASS

    otopModel->Remove, oText
    oWindow->Draw, oView

    WIDGET_CONTROL, wTopBase, SENSITIVE=1

    XMANAGER, 'd_surfview', wTopBase, EVENT_HANDLER='d_surfviewEvent', $
        CLEANUP='d_surfviewCleanup', /NO_BLOCK

    WIDGET_CONTROL, wHotKeyReceptor, /INPUT_FOCUS

end
