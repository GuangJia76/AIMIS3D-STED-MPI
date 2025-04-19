;$Id: //depot/Release/ENVI53_IDL85/idl/idldir/examples/demo/demosrc/d_flythru.pro#1 $
;
;  Copyright (c) 1997-2015, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  FILE:
;       d_flythru.pro
;
;  CALLING SEQUENCE: d_flythru
;
;  PURPOSE:
;       Shows texture mapping and an interactive flight simulation.
;
;
;  MAJOR TOPICS: Visualization
;
;  CATEGORY:
;       IDL Demo System
;
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       pro trackball__define    -  Create the trackball object
;       pro demo_gettips         - Read the tip file and create widgets
;       flythru.tip
;       elevbin.dat
;       elev_t.jpg
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
;       1/97,   DAT   - Written.
;       4/98,   ACY   - Added interactive flying controls and functionality.
;       2/2000  PCS   - New algorithm for interspersing timer events with
;                           other events.
;-
;----------------------------------------------------------------------------
;
;


FUNCTION JSTICKDISPLAY::UPDATE, delta
   compile_opt hidden

   if delta[0] LT 0 then delta[0]=delta[0] > (- 0.5) else $
                         delta[0]=delta[0] <  0.5
   if delta[1] LT 0 then delta[1]=delta[1] > (- 0.5) else $
                         delta[1]=delta[1] <  0.5
   self.arrow->setproperty, data=[[0,0],[delta[0],delta[1]]]
   return, 0
END

FUNCTION JSTICKDISPLAY::INIT, size
    compile_opt hidden

    IF (N_ELEMENTS(size) NE 2) THEN BEGIN
        PRINT, 'Jstickdisplay: size must be a two-dimensional array.'
        RETURN, 0
    ENDIF

    status = self->idlgrmodel::init(_extra=_e)
    self.size = size
    n=20
    x1 = .9*(cos(findgen(n)/(n-1)*2*!pi))
    y1 = .9*(sin(findgen(n)/(n-1)*2*!pi))
    oRing = obj_new('IDLgrPolyline', x1, y1, thick=2,color=[0,255,255])
    ringmodel = OBJ_NEW('IDLgrModel')
    ringmodel->add, oRing
    self.ringmodel = ringmodel
    self->add, ringmodel

    x2=[0.,0.]
    y2=[0.,0.]
    oArrow = obj_new('IDLgrPolyline', x2, y2, thick=3, color=[255,0,0])
    arrowmodel = OBJ_NEW('IDLgrModel')
    arrowmodel->add, oArrow
    self.arrowmodel = arrowmodel
    self->add, arrowmodel
    self.arrow = oArrow

    RETURN, 1
END

PRO jstickdisplay__define
    compile_opt hidden

    struct_hide, {jstickdisplay, $
            inherits idlgrmodel, $
            size: LONARR(2), $
            ringmodel: OBJ_NEW(), $
            arrowmodel: OBJ_NEW(), $
            arrow: OBJ_NEW(), $
            delta: FLTARR(2) $
           }
END



;  Purpose:  flight crashed into surface

pro d_flythruCrash, $
    sState       ; IN
    compile_opt hidden

    sState.omTop->SetProperty,   HIDE=1
    sState.omFixed->SetProperty,   HIDE=1
    sState.omCrash->SetProperty, HIDE=0
    demo_draw, sState.oWindow, sState.oView, debug=sState.debug

    WIDGET_CONTROL, sState.wResetButton, SENSITIVE=1
    WIDGET_CONTROL, sState.wStopButton, SENSITIVE=0
    WIDGET_CONTROL, sState.wPlayButton, SENSITIVE=1
end


;  Purpose:  Stop the flight

pro d_flythruStop, $
    sEvent       ; IN: event sturcture
    compile_opt hidden

    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

    sState.stopFlag = 1

    WIDGET_CONTROL, sState.wStartButton, SENSITIVE=1
    WIDGET_CONTROL, sState.wFileButton, SENSITIVE=1
    WIDGET_CONTROL, sState.wStopButton, SENSITIVE=0
    WIDGET_CONTROL, sState.wResetButton, SENSITIVE=1
    WIDGET_CONTROL, sState.wPlayButton, SENSITIVE=1
    WIDGET_CONTROL, sState.wClearButton, SENSITIVE=1
    WIDGET_CONTROL, sEvent.top, /CLEAR_EVENTS

    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
end

;  Purpose:  reset the surface to its initial position and orientation.

pro d_flythruReset, $
    sEvent       ; IN: event sturcture
    compile_opt hidden

    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

    sState.omTop->SetProperty,   HIDE=0
    sState.omFixed->SetProperty,   HIDE=0
    sState.omCrash->SetProperty, HIDE=1

    sState.omScale->SetProperty, $
        TRANSFORM=sState.initScaleTM
    sState.omTop->SetProperty, $
        TRANSFORM=sState.initTopTM
    demo_draw, sState.oWindow, sState.oView, debug=sState.debug

    WIDGET_CONTROL, sState.wStartButton, SENSITIVE=1
    WIDGET_CONTROL, sState.wResetButton, SENSITIVE=0
    WIDGET_CONTROL, sEvent.top, /CLEAR_EVENTS

    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
end


;----------------------------------------------------------------------------
;
;  Purpose:  Event handler for surface controls popup
;
;  Note: id of main base is stored in uvalue to allow retrieving state
;
pro d_flythruSurfBaseEvent, $
    sEvent
    compile_opt hidden

    if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ  $
        'WIDGET_KILL_REQUEST') then begin
        ; don't really destroy this popup.  Just unmap it
        ; so that it can be mapped again later if desired.
        WIDGET_CONTROL, sEvent.top, MAP=0
        RETURN
    endif

    WIDGET_CONTROL, sEvent.top, GET_UVALUE=wMainTopBase

    WIDGET_CONTROL, wMainTopBase, GET_UVALUE=sState, /NO_COPY

    WIDGET_CONTROL, sEvent.id, GET_UVALUE= uvalue

    case uvalue of

        'SURFCNTRL_UNMAP': widget_control, sState.wSurfCntrlsTlb, MAP=0

        'TEXTURE': begin
           if (sEvent.select) then begin
              sState.oSurface->SetProperty, TEXTURE_MAP=sState.oImage
              sSTate.oSurfacePaths->SetProperty, /HIDE
              WIDGET_CONTROL, sState.wInterpButton, /SENS
           endif else begin
              sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
              sSTate.oSurfacePaths->SetProperty, HIDE=0
              WIDGET_CONTROL, sState.wInterpButton, SENS=0
           endelse
           demo_draw, sState.oWindow, sState.oView, debug=sState.debug
        end

        'INTERP_BILINEAR': begin
           if (sEvent.select GT 0 AND NOT sState.bilinWarned) then begin
              tmp=dialog_message(['Use of this operation ' + $
                                   'may be unacceptably slow ' + $
                                   'due to data size or machine speed.', $
                                   'Do you want to continue ?'], $
                                   /QUESTION)
              sState.bilinWarned=1
           endif else tmp='Yes'
           if (sEvent.select EQ 0) OR $
                 (sEvent.select EQ 1 AND tmp EQ 'Yes') then begin
              sState.oSurface->SetProperty, $
                 TEXTURE_INTERP=sEvent.select
              demo_draw, sState.oWindow, sState.oView, debug=sState.debug
           endif

            if tmp eq 'No' then $
                WIDGET_CONTROL, sEvent.id, SET_BUTTON=0

        end
        'VERT_EXAG': begin
           sState.vertExag = sEvent.value
           maxz = MAX(sState.z1, MIN=minz)
           zs = [-minz/(maxz-minz), 1.0/(maxz-minz)]/(11-sState.vertExag)
           sState.oSurface->SetProperty, ZCOORD_CONV=zs
           sState.oSurfacePaths->SetProperty, ZCOORD_CONV=zs
           sState.zs = zs
           demo_draw, sState.oWindow, sState.oView, debug=sState.debug
        end
        'ZGRID': begin
           sState.oGridPaths->SetProperty, $
              HIDE=1-sEvent.select
           demo_draw, sState.oWindow, sState.oView, debug=sState.debug
        end
        'BOX': begin
           sState.oBoxPaths->SetProperty, $
              HIDE=1-sEvent.select
           demo_draw, sState.oWindow, sState.oView, debug=sState.debug
        end
    endcase

    if XREGISTERED('demo_tour') eq 0 then $
        WIDGET_CONTROL, sState.wHotkeyReceptor, /INPUT_FOCUS

    ; put the state back in the main top level base
    WIDGET_CONTROL, wMainTopBase, SET_UVALUE=sState, /NO_COPY

end
;----------------------------------------------------------------------------
pro d_flythruAdvance, sState
compile_opt hidden

    case 1 of
        sState.curPathPlaying and sState.btndown eq 0: begin
            if (sState.curPathIndx GE N_ELEMENTS(*(sState.curPath))) $
            then $
                sState.curPathIndx = 0
            if (sState.stopFlag NE 1) then begin
                sState.omTop->SetProperty, TRANSFORM= $
                    *((*(sState.curPath))[sState.curPathIndx])
                sState.curPathIndx = sState.curPathIndx + 1
                demo_draw, $
                    sState.oWindow, $
                    sState.oView, $
                    DEBUG=sState.debug
            endif
        end

        sState.btndown eq 0: begin
            zMoveDelta = 0.001 * sState.airspeed
            rollFactor = 1
            pitchFactor = 1
            sState.rollDelta= $
                sState.xdelta * sState.rollSensitivity * rollFactor
            sState.roll = sState.rollDelta
            sState.pitchDelta= $
                sState.ydelta * sState.pitchSensitivity * pitchFactor
            sState.pitch = sState.pitchDelta

            ; Roll
            ; Move to center of desired rotation
            sState.omTop->Translate, -.5, 0, 0
            ; Do the rotations
            sState.omTop->Rotate, [0,0,1], sState.roll
            ; throw in a little spin (yaw) also to help move around
            ; in the view
            yawFactor =3
            sState.omTop->Rotate, [0,1,0], (sState.roll)/yawFactor
            ; Move back
            sState.omTop->Translate, .5, 0, 0

            ; Pitch
            ; Move to center of desired rotation
            sState.omTop->Translate, 0, 0, sState.eyePos
            ; Do the rotations
            sState.omTop->Rotate, [1,0,0], -(sState.pitch)
            ; Move back
            sState.omTop->Translate, 0, 0, -(sState.eyePos)

            ; Move in Z to simulate forward movement
            sState.omTop->Translate, 0, 0, zMoveDelta

            ; save the latest movement
            sState.omTop->GetProperty, TRANSFORM=tm
            *(sState.curPath) = [ $
                *(sState.curPath), PTR_NEW(tm, /NO_COPY) $
                ]

            if (sState.ignoreCrashes EQ 0) then begin
                ; check for impact
                ctm = sState.oSurface->GetCTM()
                ctm_inv = INVERT(ctm, status)
                curSurfIntercept = [0.5,0.5,sState.eyePos, 1] # ctm_inv
                surfCoord = ROUND(curSurfIntercept[0:1])
                sizeZ = size(sState.z1)

                if (surfCoord[0] GE 0 AND surfCoord[0] LT sizeZ[1]) AND $
                   (surfCoord[1] GE 0 AND surfCoord[1] LT sizeZ[2]) $
                THEN BEGIN
                    unscaledSurfZ = sState.z1[surfCoord[0], surfCoord[1]]
                    scaledSurfZ = sState.zs[0] + sState.zs[1] * $
                          (sState.z1[surfCoord[0], surfCoord[1]])

                    sState.oWindow->GetProperty, $
                        SCREEN_DIMENSIONS=screen_dims

                    curSurfIntercept[2] = sState.zs[0] + $
                        sState.zs[1] * curSurfIntercept[2]

                    deltaZ = curSurfIntercept[2] - scaledSurfZ
                    if sState.debug then begin
                        print, 'i: ', surfCoord[0], $
                            ', j: ', surfCoord[1], $
                            ', z[i,j]: ', sState.z1[surfCoord[0], $
                                surfCoord[1]], $
                            ', scaled z: ', scaledSurfZ, $
                            ', mapped eyePos: ', curSurfIntercept[2], $
                            ', deltaZ: ', deltaz
                    endif

                    if (sState.deltaZValid) AND $
                    ((sState.deltaZ GT 0 AND deltaZ LT 0) $
                        OR $
                     (sState.deltaZ LT 0 AND deltaZ GT 0)) THEN BEGIN
                        ; impact
                        sState.deltaZValid=0
                        sState.stopFlag = 1
;                       WIDGET_CONTROL, sEvent.top, /CLEAR_EVENTS
                        d_flythruCrash, sState
                        RETURN
                    endif else $
                        sState.deltaZValid=1
                    sState.deltaZ = deltaZ
                endif
            endif

            demo_draw, sState.oWindow, sState.oView, debug=sState.debug
        end
        else:
    endcase
end
;----------------------------------------------------------------------------
;
;  Purpose:  Event handler.
;
pro d_flythruEvent, $
    sEvent
    compile_opt hidden

    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
    demo_record, $
        sEvent, $
        'd_flythruEvent', $
        filename=sState.record_to_filename, $
        CW=sState.wCrashButton
    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY

    if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ  $
        'WIDGET_KILL_REQUEST') then begin
        WIDGET_CONTROL, sEvent.top, /DESTROY
        RETURN
    endif

    WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
    IF (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ 'WIDGET_TIMER') then begin
        if sState.stopFlag EQ 1 and sState.btndown eq 0 then begin
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
            d_flythruStop, sEvent
            RETURN
        endif
        d_flythruAdvance, sState

        ; OK, now we have handled the timer event, so there is not
        ; one pending.  If other events are pending, handle one of
        ; those with the following call to WIDGET_EVENT.
        ; (BAD_ID handles the case where this application has been
        ; closed.)

        wDraw = sState.wDraw
        wSurfCntrlsTlb = sState.wSurfCntrlsTlb
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY

        if XREGISTERED('demo_tour') gt 0 then $
            RETURN

        void = WIDGET_EVENT( $
            [sEvent.top, wSurfCntrlsTlb], $
            BAD_ID=bad_id, $
            /NOWAIT $
            )
        if bad_id ne 0 then $
            RETURN

        ; Launch a new timer event.
        WIDGET_CONTROL, wDraw, TIMER=0
        RETURN
    endif

    WIDGET_CONTROL, sEvent.id, GET_UVALUE=uvalue

    case uvalue of
        'HOTKEY': begin
            case STRUPCASE(sEvent.ch) of
                ' ': begin
                    case sState.stopFlag of
                        0: begin
                            sState.omScale->SetProperty, $
                                TRANSFORM=sState.initScaleTM
                            sState.omTop->SetProperty, $
                                TRANSFORM=sState.initTopTM

                            demo_draw, sState.oWindow, sState.oView, $
                                debug=sState.debug
                        end
                        else: begin
                            wStartButton = sState.wStartButton
                            wResetButton = sState.wResetButton
                            tlb = sEvent.top
                            sState.omCrash->GetProperty, HIDE=hide_crashed

                            WIDGET_CONTROL, $
                                sEvent.top, $
                                SET_UVALUE=sState, $
                                /NO_COPY
                            if hide_crashed eq 0 then $
                                d_flythruEvent, { $
                                    id: wResetButton, $
                                    top: tlb, $
                                    handler: 0L $
                                    }
                            d_flythruEvent, { $
                                id: wStartButton, $
                                top: tlb, $
                                handler: 0L $
                                }
                            WIDGET_CONTROL, $
                                sEvent.top, $
                                GET_UVALUE=sState, $
                                /NO_COPY
                        end
                    endcase
                end
                else:
            endcase
        end

        'RESET' : begin
            sState.stopFlag = 1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
            d_flythruReset, sEvent
            RETURN
        end

        'PLAY' : begin
            sState.omCrash->GetProperty, HIDE=tmpHide
            ; if crashed, reset to view surface object
            if (tmpHide EQ 0) then begin
               WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
               d_flythruReset, sEvent
               WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            endif
            WIDGET_CONTROL, sState.wStartButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wStopButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wResetButton, SENSITIVE=0
            widget_control, sState.wPlayButton, SENS=0
            widget_control, sState.wClearButton, SENS=0
            sState.stopFlag = 0
            sState.curPathPlaying = 1
            WIDGET_CONTROL, sState.wDraw, TIMER=0
        end

        'CLEAR' : begin
            widget_control, sState.wClearButton, SENS=0
            widget_control, sState.wPlayButton, SENS=0
            for i = 0, N_ELEMENTS(*(sState.curPath)) - 1 do begin
              PTR_FREE, (*(sState.curPath))[i]
            endfor
            PTR_FREE, sState.curPath
            sState.curPathIndx = 0
        end

        ;  Handle the event generated within the drawing area
        ;
        'GO': begin
            sState.stopFlag = 0
            widget_control, sState.wStartButton, SENS=0
            widget_control, sState.wStopButton, SENS=1
            widget_control, sState.wResetButton, SENS=0
            widget_control, sState.wPlayButton, SENS=0
            widget_control, sState.wClearButton, SENS=0
            sState.omTop->GetProperty, TRANSFORM = initTopTM
            if (NOT PTR_VALID(sState.curPath))  then begin
               sState.curpath=PTR_NEW(PTRARR(1))
               *(sState.curPath)[0] = PTR_NEW(initTopTM)
            endif else begin
               *(sState.curPath)=[*(sState.curPath), PTR_NEW(initTopTM)]
            endelse
            WIDGET_CONTROL, sState.wDraw, TIMER=0
        end
        'STOP' : begin
            sState.stopFlag = 1
            sState.curPathPlaying = 0
        end
        'AIRSPEED': begin
            sState.airspeed = sEvent.value
        end
        'PITCH_SENS': begin
            sState.pitchSensitivity = sEvent.value
         end
        'ROLL_SENS': begin
            sState.rollSensitivity = sEvent.value
         end
        'DRAW': begin
            if XREGISTERED('demo_tour') eq 0 then $
                WIDGET_CONTROL, sState.wHotkeyReceptor, /INPUT_FOCUS

            ; convert the event coords to normalized coords with 0,0 at
            ; center of screen.
            delta=[sEvent.x/sState.xdim, sEvent.y/sState.ydim] - 0.5

            ; update the jstick
            temp = sState.oJstickDisplay->Update( delta )
            demo_draw, sState.oWindow2, sState.oView2, debug=sState.debug
            sState.xDelta = delta[0]
            sState.yDelta = delta[1]

            ;  Expose.
            ;
            if (sEvent.type eq 4) then begin
                demo_draw, sState.oWindow, sState.oView, debug=sState.debug
                WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
                RETURN
            endif

            ;  Handle trackball update
            ;

            ;  Button press.
            ;
            if (sEvent.type EQ 0) then begin
                if sEvent.press EQ 1 OR sEvent.press EQ 4 then begin
                   sState.btndown = 1B
                   sState.oWindow->setproperty, QUALITY=0
                   ; for trackball update call, save which button pressed
                   sState.whichButton = sEvent.press
                   ;WIDGET_CONTROL, sState.wDraw, /DRAW_MOTION
                endif else begin
                   ; middle mouse button is used for "pause"
                   if sState.stopFlag EQ 0 then begin
                      ; currently flying, stop flying
                      sState.stopFlag = 1
                      sState.curPathPlaying = 0
                   endif else begin
                      ; currently stopped, start flying
                      sState.stopFlag = 0
                      widget_control, sState.wStartButton, SENS=0
                      widget_control, sState.wStopButton, SENS=1
                      widget_control, sState.wResetButton, SENS=0
                      widget_control, sState.wPlayButton, SENS=0
                      widget_control, sState.wClearButton, SENS=0
                      sState.omTop->GetProperty, TRANSFORM = initTopTM
                      if (NOT PTR_VALID(sState.curPath))  then begin
                         sState.curpath=PTR_NEW(PTRARR(1))
                         *(sState.curPath)[0] = PTR_NEW(initTopTM)
                      endif else begin
                         *(sState.curPath)= $
                            [*(sState.curPath), PTR_NEW(initTopTM)]
                      endelse
                      WIDGET_CONTROL, sState.wDraw, TIMER=0
                   endelse
                endelse
            endif

            ; test for movement
            bHaveTransform = sState.oTrack->Update( $
                sEvent, $
                TRANSFORM=qmat, $
                MOUSE=sState.whichButton, $
                TRANSLATE=(sState.whichButton EQ 4) $
                )
            if (bHaveTransform NE 0) then begin
                sState.omTop->Translate, -.5, 0, 0
                sState.omTop->GetProperty, TRANSFORM=t
                mt = t # qmat
                sState.omTop->SetProperty,TRANSFORM=mt
                sState.omTop->Translate,  .5, 0, 0
            endif

            ;  Button motion.
            ;
            if sEvent.type eq 2 then begin
                case 1 of
                    sState.btndown eq 1b: begin
                        if (bHaveTransform) then begin
                            demo_draw, sState.oWindow, sState.oView, $
                                debug=sState.debug
                        endif
                    end
                    sState.stopFlag eq 0: begin
                        ; Some machines (i.e. Silicone Graphics')
                        ; seem to favor mouse events over timer
                        ; events.  On these machines, timer
                        ; events are postponed as long as the
                        ; user is rolling the mouse.  In order
                        ; for the plane to continue flying (rather
                        ; than freeze) while the user is rolling
                        ; the mouse on those machines, we manually
                        ; advance the plane here.
                        d_flythruAdvance, sState
                    end
                    else:
                endcase
            endif

            ;  Button release.
            ;
            if (sEvent.type eq 1) then begin
                if (sState.btndown EQ 1b) then begin
                    sState.oWindow->SetProperty, QUALITY=2
                    demo_draw, sState.oWindow, sState.oView, $
                        debug=sState.debug
                endif
                sState.btndown = 0b

                ; Invalidate the deltaZ, since we have moved drastically
                ; with trackball
                sState.deltaZ=0
                sState.deltaZValid=0
            endif

        end   ;   of DRAW


        'IGNORE_CRASHES': sState.ignoreCrashes=sEvent.select

        'SURFCNTRL_MAP'  : widget_control, sState.wSurfCntrlsTlb, /MAP
        'DRAW2':

        'QUIT' : BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sEvent.top, /DESTROY
            RETURN
        end   ; of QUIT

    endcase
    WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY

end

;-----------------------------------------------------------------
;
;    PURPOSE : cleanup procedure. restore colortable, destroy objects.
;
pro d_flythruCleanup, wTopBase
compile_opt hidden

    WIDGET_CONTROL, wTopBase, GET_UVALUE=sState, /NO_COPY

    ;  Destroy the top objects.
    ;
    OBJ_DESTROY, sState.oView
    OBJ_DESTROY, sState.omTop
    OBJ_DESTROY, sState.omSurf
    OBJ_DESTROY, sState.omScale
    OBJ_DESTROY, sState.omRotZ
    OBJ_DESTROY, sState.omRotX
    if (PTR_VALID(sState.curPath))  then begin
       for i = 0, N_ELEMENTS(*(sState.curPath)) - 1 do begin
          PTR_FREE, (*(sState.curPath))[i]
       endfor
    endif
    PTR_FREE, sState.curPath
    OBJ_DESTROY, sState.oSurface
    OBJ_DESTROY, sState.oSurfacePaths
    OBJ_DESTROY, sState.oImage
    OBJ_DESTROY, sState.oTrack
    OBJ_DESTROY, sState.oText
    OBJ_DESTROY, sState.oFont

    OBJ_DESTROY, sState.oView2
    OBJ_DESTROY, sState.oJstickDisplay
    OBJ_DESTROY, sState.oContainer

    ;  Restore the color table.
    ;
    TVLCT, sState.colorTable

    if WIDGET_INFO(sState.groupBase, /VALID_ID) then $
        WIDGET_CONTROL, sState.groupBase, /MAP

end   ;  of d_flythruCleanup

;-----------------------------------------------------------------
;
;    PURPOSE : show the texture mapping capability
;
PRO d_flythru_gj, $
    SURFACE=z, $
    RECORD_TO_FILENAME=record_to_filename, $
    IMAGE = image, $             ; If set, display texture map,
                                 ; otherwise surface paths
    INTERP_BILINEAR = interp_bilinear, $
    NO_GRID = no_grid, $         ; IN: (opt) hide the additional grid at z=0
    NO_BOX = no_box, $           ; IN: (opt) hide the center box
    TEXCOORDS = texcoords, $     ; IN: (opt) texture coordinates
    DEBUG = debug, $             ; IN: (opt)
    GROUP=group, $               ; IN: (opt) group identifier
    APPTLB = appTLB              ; OUT: (opt) TLB of this application

    ; Check the validity of the group identifier.
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
    xdim = screenSize[0]*0.6
    ydim = xdim*0.8

    ;  Get the current color vectors to restore
    ;  when this application is exited.
    ;
    TVLCT, savedR, savedG, savedB, /GET

    ;  Build color table from color vectors
    ;
    colorTable = [[savedR],[savedG],[savedB]]

    ;  Create widgets.
    ;
    if (N_ELEMENTS(group) EQ 0) then begin
        wTopBase = WIDGET_BASE(/COLUMN, $
            TITLE="Satellite and Topographic Data Flythrough", $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            UNAME='d_flythru:tlb', $
            TLB_FRAME_ATTR=1, MBAR=barBase)
    endif else begin
        wTopBase = WIDGET_BASE(/COLUMN, $
            TITLE="Satellite and Topographic Data Flythrough", $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            GROUP_LEADER=group, $
            UNAME='d_flythru:tlb', $
            TLB_FRAME_ATTR=1, MBAR=barBase)
    endelse

        ;  Create the menu bar. It contains the file/quit,
        ;  edit/ shade-style, help/about.
        ;
        wFileButton = WIDGET_BUTTON(barBase, VALUE='File', /MENU)

            wQuitButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Quit', UVALUE='QUIT')

         ;  Create a sub base of the top base (wTopBase)
         ;
         subBase = WIDGET_BASE(wTopBase, /ROW)

             ;  Create the left Base that contains GUI controls.
             ;
             wLeftbase = WIDGET_BASE(subBase, /ALIGN_CENTER, /COLUMN)

                wFlightbase = WIDGET_BASE(wLeftBase, /ALIGN_CENTER, $
                    /COLUMN, /FRAME)

                    label = widget_label(wFlightBase, $
                        VALUE='Flight Controls')

                    wStartButton = WIDGET_BUTTON(wFlightBase, $
                        VALUE='Start', UVALUE='GO', $
                        UNAME='d_flythru:Start')

                    wStopButton = WIDGET_BUTTON(wFlightBase, $
                        VALUE='Stop', UVALUE='STOP', $
                        UNAME='d_flythru:Stop')
                    widget_control, wStopButton, SENS=0

                    wResetButton = WIDGET_BUTTON(wFlightBase, $
                        VALUE='Reset', UVALUE='RESET')
                    widget_control, wResetButton, SENS=0

                    wPlayButton = WIDGET_BUTTON(wFlightBase, $
                        VALUE='Replay Path', UVALUE='PLAY')
                    widget_control, wPlayButton, SENS=0
                    wClearButton = WIDGET_BUTTON(wFlightBase, $
                        VALUE='Clear Path', UVALUE='CLEAR')
                    widget_control, wClearButton, SENS=0

                    airspeed = 7
                    airspeedSlider = widget_slider(wFlightBase, MIN=-20,$
                        MAX=20,VALUE=airspeed, $
                        uname='d_flythru:airspeed', $
                        title="Airspeed", $
                        uval='AIRSPEED')
                    pitchSensitivity=1
                    pitchSlider = widget_slider(wFlightBase, MIN=1,MAX=10,$
                        VALUE=pitchSensitivity, $
                        uname='d_flythru:pitch', $
                        title="Pitch Sensitivity", $
                        uval='PITCH_SENS')
                    rollSensitivity=1
                    rollSlider = widget_slider(wFlightBase, MIN=1,MAX=10, $
                        VALUE=rollSensitivity, $
                        title="Roll Sensitivity", $
                        uname='d_flythru:roll', $
                        uval='ROLL_SENS')

                    wcrashbase = WIDGET_BASE(wFlightBase, /ALIGN_CENTER, $
                        /COLUMN, /NONEXCLUS)

                        wCrashButton = WIDGET_BUTTON(wcrashBase, $
                            VALUE='Ignore Crashes', $
                            UVALUE='IGNORE_CRASHES', $
                            UNAME='d_flythru:Ignore Crashes')

                    xdim2 = (ydim2 = 100)
                    wDraw2 = WIDGET_DRAW(wFlightBase, $
                        XSIZE=xdim2, YSIZE=ydim2, $
                        GRAPHICS_LEVEL=2, $
                        UVALUE='DRAW2', $
                        UNAME='d_flythru:joydraw', $
                        /EXPOSE_EVENTS, RETAIN=0 )

                wSurfaceButton = WIDGET_BUTTON(wLeftBase, $
                    VALUE='Surface Controls', UVALUE='SURFCNTRL_MAP')

            ;  Create the right Base that has the drawing area
            ;
            wRightBase = WIDGET_BASE(subBase)

                ;  Draw area.
                ;
                wDraw = WIDGET_DRAW(wRightBase, $
                    GRAPHICS_LEVEL=2, $
                    XSIZE=xdim, YSIZE=ydim, /BUTTON_EVENTS, $
                    /MOTION_EVENTS, $
                    UVALUE='DRAW', $
                    UNAME='d_flythru:draw', $
                    RETAIN=0, /EXPOSE_EVENT)
                wHotKeyReceptor = WIDGET_TEXT(wRightBase, $
                    /ALL_EVENTS, $
                    UVALUE='HOTKEY', $
                    UNAME='d_flythru:hotkey')

        ;  Create tips texts.
        ;
        wStatusBase = WIDGET_BASE(wTopBase, MAP=0, /ROW)

    ;  Create popup Surface Controls top-level base.
    ;
    wSurfCntrlsTlb = WIDGET_BASE(MAP=0, $
        XOFFSET=100, YOFFSET=100, $
        GROUP=wTopBase, /FLOATING, $
        /TLB_KILL_REQUEST_EVENTS, $
        TITLE='Surface Controls')

            ;  To avoid color flashing that can occur when the
            ;  user clicks anywhere in wSurfCntrlsTlb on
            ;  256-color Windows displays, hide an IDLgrWindow
            ;  in wSurfCntrlsTlb.
            ;
            void = WIDGET_DRAW( $
                WIDGET_BASE(wSurfCntrlsTlb, MAP=0), $
                GRAPHICS_LEVEL=2)

            ;  Populate Surface Controls panel with widgets.
            ;
            wCntrlsBase = WIDGET_BASE(wSurfCntrlsTlb, /ALIGN_CENTER, $
                /COLUMN)
                wButtonBase = WIDGET_BASE(wCntrlsBase, /ALIGN_CENTER, $
                    /COLUMN, /NONEXCLUS)

                    wTextureButton = WIDGET_BUTTON(wButtonBase, $
                        VALUE='Texture Map', UVALUE='TEXTURE')

                    wInterpButton = WIDGET_BUTTON(wButtonBase, $
                        VALUE='Bilinear Interp.', $
                        UVALUE='INTERP_BILINEAR')

                    wGridButton = WIDGET_BUTTON(wButtonBase, $
                        VALUE='Grid at Z=0', UVALUE='ZGRID')

                    wBoxButton = WIDGET_BUTTON(wButtonBase, $
                        VALUE='Center Box', UVALUE='BOX')

                vertExag=7
                vertExagSlider = widget_slider(wCntrlsBase, $
                    MIN=1,MAX=10, $
                    VALUE=vertExag, $
                    title="Vertical Exag.", $
                    uval='VERT_EXAG')

                wSurfDoneButton = WIDGET_BUTTON(wCntrlsBase, $
                    VALUE='Done', UVALUE='SURFCNTRL_UNMAP')

    ;  Realize the base widgets.
    ;
    WIDGET_CONTROL, wTopBase, /REALIZE
    WIDGET_CONTROL, wSurfCntrlsTlb, /REALIZE

    ;  Returns the top level base in the appTLB keyword
    ;
    appTLB = wTopBase

    ;  Get the tips
    ;
    sText = demo_getTips(demo_filepath('flythru.tip', $
                             SUBDIR=['examples','demo', 'demotext']), $
                         wTopBase, $
                         wStatusBase)


    WIDGET_CONTROL, wTopBase, SENSITIVE=0

    ;  Grab the window id of the drawable.
    ;
    WIDGET_CONTROL, wDraw, GET_VALUE=oWindow
    WIDGET_CONTROL, wDraw2, GET_VALUE=oWindow2

    ;  Create view.
    ;
    myview = [0, 0, 1, 1]
    eyePos = 2.0
    zclip = [eyePos-.01,-2.0]
    oView = OBJ_NEW('idlgrview', $
       ; PROJECTION=2, EYE=3, ZCLIP=[2.9,-10.0], $
        PROJECTION=2, EYE=eyePos, ZCLIP=zclip, $
        VIEWPLANE_RECT=myview, COLOR=[0, 0, 0])

    ;  center the startup text
    ;
    ;textLocation = [myview[0]+0.5*myview[2], myview[1]+0.5*myview[3]]

    ;  Create and display the PLEASE WAIT text.
    ;
    ;oFont = OBJ_NEW('IDLgrFont', 'Helvetica', SIZE=10)
    ;oText = OBJ_NEW('IDLgrText', $
    ;    'Starting up  Please wait...', $
    ;    ALIGN=0.5, $
    ;    LOCATION=textLocation, $
    ;    COLOR=[255,255,0], FONT=oFont)
    xcoords = [52b,72b,69b,82b,69b,0b,66b,69b,0b,68b,82b]
    ycoords  = [65b,71b,79b,78b,83b,0b,72b,69b,82b,69b,1b]
    locations = [xcoords, ycoords]
    oFont = OBJ_NEW('IDLgrFont', 'Helvetica', SIZE=10)
    oText = OBJ_NEW('IDLgrText', $
        string(locations+32b), $
        ALIGN=0.5, $
        LOCATION=[.5,.5, -6], $
        BASELINE=[1,0,0], $   ;default
        COLOR=[255,0,0], FONT=oFont)


    ;  Create model.
    ;
    omTop = OBJ_NEW('idlgrmodel')
    omSurf = OBJ_NEW('idlgrmodel')
    omScale = OBJ_NEW('idlgrmodel')
    omTrans = OBJ_NEW('idlgrmodel')
    omRotZ = OBJ_NEW('idlgrmodel')
    omRotX = OBJ_NEW('idlgrmodel')
    omTransPaths = OBJ_NEW('idlgrmodel')

    omFixed = OBJ_NEW('idlgrmodel')
    omCrash = OBJ_NEW('idlgrmodel', HIDE=1)

    omTop->Add, omSurf
    omSurf->Add, omScale
    omScale->Add, omTrans
    omTrans->Add, omRotZ
    omRotZ->Add, omRotX
    omRotX->Add, omTransPaths

    ;  Place the model in the view.
    ;
    oView->Add, omTop
    oView->Add, omFixed   ; for objects which don't move
    oView->Add, omCrash

    oContainer = OBJ_NEW('IDLgrContainer')
    oContainer->Add, oView

    ;  Add the trackball object for interactive change
    ;  of the scene orientation
    ;
    oTrack = OBJ_NEW('Trackball', [xdim/2.0, ydim/2.0], xdim/2.0)
    oContainer->Add, oTrack

    ;  Add the text.
    ;
    omTop->Add, oText

    ;  Draw the starting up screen.
    ;
    ;oWindow->Draw, oView

    ;  Surface data is read from elevation data file.
    ;
    IF N_ELEMENTS(z) EQ 0 then begin
       z = BYTARR(64,64, /NOZERO)
       OPENR, lun, demo_filepath('elevbin.dat', $
           SUBDIR=['examples','data']), $
           /GET_LUN
       READU, lun, z
       FREE_LUN, lun
       z = REVERSE(TEMPORARY(z), 2)

       ;  Create texture map.
       ;
       READ_JPEG, demo_filepath('elev_t.jpg', $
           SUBDIR=['examples','data']), $
           idata, TRUE=3
       idata = REVERSE(TEMPORARY(idata), 2)
       
       ;GJ 2019/8/12, test texture mapping
       imagefile='D:\AIMIS_3D\lucy2.jpg';kelly.jpg';lucy2.jpg';kelly.jpg'; gj.jpg
       read_jpeg, imagefile, idata1
       
       idata1_0 = CONGRID(idata1, 3, 512, 512)
       FOR i=0,511 DO BEGIN
        FOR j=0, 511 DO BEGIN
          idata[i, j, *] = idata1_0[*, i, j]
        ENDFOR
       ENDFOR

       oImage = OBJ_NEW('IDLgrImage', idata, INTERLEAVE=2)
    endif else begin
       if keyword_set(image) then $
          oImage = OBJ_NEW('IDLgrImage', image)
    endelse

    sz = SIZE(z)
    minx = (miny = 0)
    maxx = sz[1] - 1
    maxy = sz[2] - 1
    maxz = MAX(z, MIN=minz)

    ;  Compute coordinate conversion to normalize.
    ;
    xs = [-minx/(maxx-minx), 1.0/(maxx-minx)]
    ys = [-miny/(maxy-miny), 1.0/(maxy-miny)]
    zs = [-minz/(maxz-minz), 1.0/(maxz-minz)]/(11-vertExag)

    ;  Create the surface object.
    ;
    oSurface = OBJ_NEW('IDLgrSurface', z, $
          STYLE=2, $
          SHADING=1, $
          COLOR=[255,255,255], $
          TEXTURE_MAP=oImage, $
          XCOORD_CONV=xs, $
          YCOORD_CONV=ys, $
          ZCOORD_CONV=zs )

    oSurfacePaths = OBJ_NEW('IDLgrSurface', z, $
       STYLE=1, $
       XCOORD_CONV=xs, $
       YCOORD_CONV=ys, $
       ZCOORD_CONV=zs )

    ; If user supplied surface there might not be a corresponding image
    IF (n_elements(oImage) GT 0) then begin
       oSurface->SetProperty, TEXTURE_MAP=oImage
       oSurfacePaths->SetProperty, /HIDE
       ; this is cool but slow on most machines
       oSurface->SetProperty, TEXTURE_INTERP=keyword_set(interp_bilinear)
       WIDGET_CONTROL, wTextureButton, /SET_BUTTON
    endif else begin
       oImage=OBJ_NEW()
       oSurface->SetProperty, TEXTURE_MAP=oImage
       oSurfacePaths->SetProperty, HIDE=0
       WIDGET_CONTROL, wTextureButton, SENS=0
       WIDGET_CONTROL, wInterpButton, SENS=0
    endelse

    gridSize=5
    ysGrid = (xsGrid = ([0, 1.0/gridSize]*gridSize)-.5)
    oGridPaths = OBJ_NEW('IDLgrSurface', intarr(gridSize,gridSize), $
       STYLE=1, $
       COLOR=[255,255,0] , $
       XCOORD_CONV=xsGrid, $
       YCOORD_CONV=ysGrid )
    if keyword_set(NO_GRID) then begin
       oGridPaths->SetProperty, /HIDE
    endif else widget_control, wGridButton, /SET_BUTTON

    lo = .499
    hi = .501
    ;zs = (ys = (xs = [0, 1]))         ; polyline coords already normalized
    boxCoordsX= [lo,hi,hi,lo,lo]
    boxCoordsY= [lo,lo,hi,hi,lo]
    boxCoordsZ= fltarr(5)+zclip[0]-.05
    oBoxPaths = OBJ_NEW('IDLgrPolyline', $
       boxCoordsX, boxCoordsY, boxCoordsZ, $
       LINESTYLE=0, $
       COLOR=[0,255,255])
    if keyword_set(NO_BOX) then begin
       oBoxPaths->SetProperty, /HIDE
    endif else widget_control, wBoxButton, /SET_BUTTON
    omFixed->Add, oBoxPaths


    omRotX->Add, oSurface
    omTransPaths->Add, oSurfacePaths
    omRotX->Add, oGridPaths

    z1 = z
    oSurface->SetProperty, DATAZ=z1
    oSurfacePaths->SetProperty, DATAZ=z1
    ; translate paths so they show on top of surface
    omTransPaths->Translate, 0, 0, .002

    ;  Create a light.
    ;
    oLight1 = OBJ_NEW('IDLgrLight', $
        TYPE=0, INTENSITY=0.25, COLOR=[255,255,255])
    omTop->Add, oLight1
    ;oLight2 = OBJ_NEW('IDLgrLight', LOCATION=[-2,2,-2], TYPE=1)
    oLight2 = OBJ_NEW('IDLgrLight', LOCATION=[0,20,10], TYPE=1)
    omTop->Add, oLight2

    lineCoordsX = [.5,0.450,0.453,0.370,0.354,0.335,0.255,0.189,0.207,0.115]
    lineCoordsY = [.5,0.523,0.569,0.558,0.662,0.722,0.697,0.782,0.914,0.928]
    olinePaths1 = OBJ_NEW('IDLgrPolyline', $
       lineCoordsX, lineCoordsY, $
       LINESTYLE=0, $
       COLOR=[255,255,255])
    omCrash->Add, olinePaths1

    lineCoordsX = [.5,0.523,0.524,0.481,0.493,0.531,0.535,0.488,0.569,0.644]
    lineCoordsY = [.5,0.549,0.623,0.646,0.720,0.771,0.863,0.917,0.944,0.944]
    olinePaths2 = OBJ_NEW('IDLgrPolyline', $
       lineCoordsX, lineCoordsY, $
       LINESTYLE=0, $
       COLOR=[255,255,255])
    omCrash->Add, olinePaths2

    lineCoordsX = [.5,0.531,0.552,0.563,0.595,0.649,0.701,0.778,0.865,0.962]
    lineCoordsY = [.5,0.491,0.521,0.456,0.507,0.500,0.447,0.507,0.447,0.403]
    olinePaths3 = OBJ_NEW('IDLgrPolyline', $
       lineCoordsX, lineCoordsY, $
       LINESTYLE=0, $
       COLOR=[255,255,255])
    omCrash->Add, olinePaths3

    lineCoordsX = [.5,0.503,0.547,0.531,0.566,0.622,0.632,0.783,0.851,0.925]
    lineCoordsY = [.5,0.456,0.421,0.363,0.322,0.303,0.243,0.229,0.162,0.118]
    olinePaths4 = OBJ_NEW('IDLgrPolyline', $
       lineCoordsX, lineCoordsY, $
       LINESTYLE=0, $
       COLOR=[255,255,255])
    omCrash->Add, olinePaths4

    lineCoordsX = [.5,0.465,0.441,0.455,0.438,0.385,0.392,0.267,0.170,0.122]
    lineCoordsY = [.5,0.435,0.366,0.315,0.282,0.269,0.201,0.169,0.130,0.194]
    olinePaths5 = OBJ_NEW('IDLgrPolyline', $
       lineCoordsX, lineCoordsY, $
       LINESTYLE=0, $
       COLOR=[255,255,255])
    omCrash->Add, olinePaths5

    ;  Set the initial view.
    ;
    zScale = 1.
    omScale->Scale, 1.0, 1.0, zScale
    ;omRotZ->Rotate, [1,0,0],-35
    omRotZ->Rotate, [1,0,0],-90
    ;omTrans->Translate, 0,0.25,0
    omTop->Translate, 0,0.25,0.5

    ;  Get the initial transformation.
    ;
    omTop->GetProperty, TRANSFORM = initTopTM
    omScale->GetProperty, TRANSFORM = initScaleTM
    omRotZ->GetProperty, TRANSFORM = initZRotationTM
    omRotX->GetProperty, TRANSFORM = initXRotationTM


    oView2=OBJ_NEW('IDLgrView')
    oModel2 = OBJ_NEW('IDLgrModel')
    oView2->Add, oModel2
    oJstickDisplay = OBJ_NEW('jstickdisplay', [xdim2,ydim2])
    oModel2->Add, oJstickDisplay
    demo_draw, oWindow2, oView2, debug=keyword_set(debug)

    if n_elements(record_to_filename) eq 0 then $
        record_to_filename = ''

    ;  Create the state structure.
    ;
    sState = { $
        xs: xs, $
        ys: ys, $
        zs: zs, $
        deltaZ: 0.0, $
        deltaZValid: 0B, $
        btndown: 0b, $                      ; 0 = not pressed, 1= pressed
        wDraw: wDraw, $                     ; Widget draw ID
        omTop: omTop, $      ; Models
        omCrash: omCrash, $
        omFixed: omFixed, $
        omSurf: omSurf, $
        omScale: omScale, $
        omTrans: omTrans, $
        omRotZ: omRotZ, $
        omRotX: omRotX, $
        curPath: PTR_NEW(), $
        curPathIndx: 0L, $
        curPathPlaying: 0B, $
        eyePos: eyePos, $
        OSurface: oSurface, $               ; Surface object
        OSurfacePaths: oSurfacePaths, $     ; Surface object
        OGridPaths: oGridPaths, $           ; Surface object
        OBoxPaths: oBoxPaths, $
        OView: oView, $                     ; View object
        oTrack: oTrack, $
        WTopBase : wTopbase, $              ; Top level base
        WFileButton : wFileButton, $        ; Buttons ID
        WStartButton : wStartButton, $
        WStopButton : wStopButton, $
        WResetButton : wResetButton, $
        WPlayButton : wPlayButton, $
        WClearButton : wClearButton, $
        WQuitButton : wQuitButton, $
        WCrashButton: wCrashButton, $
        wSurfCntrlsTlb : wSurfCntrlsTlb, $  ; Popup dialog top-level base
        wInterpButton : wInterpButton, $
        wGridButton : wGridButton, $
        wBoxButton : wBoxButton, $
        wHotkeyReceptor: wHotkeyReceptor, $
        OImage: oImage, $                   ; Image object
        ColorTable: colorTable, $           ; Color table to restore at exit
        Z1: z1, $                           ; Surface height data
        InitScaleTM: initScaleTM, $         ; Initial scale matrix
        InitTopTM: initTopTM, $             ; Initial transformation matrix
        InitZRotationTM: initZRotationTM, $ ; Initial rotation of surface
        InitXRotationTM: initXRotationTM, $ ; Initial rotation of surface
        OWindow: oWindow, $                 ; Window object
        OText: oText, $                     ; Starting up text object
        OFont: oFont, $                     ; Starting up text font object
        StopFlag: 1, $
        groupBase: groupBase, $             ; Base of Group Leader
        oWindow2: oWindow2, $
        oView2: oView2, $
        oContainer: oContainer, $
        oJstickDisplay: oJstickDisplay, $
        airspeed: airspeed, $
        xDelta: 0.0, $
        yDelta: 0.0,  $
        rollDelta: 0.0, $
        rollSensitivity: rollSensitivity, $
        roll: 0.0, $
        pitchDelta: 0.0, $
        pitchSensitivity: pitchSensitivity, $
        pitch: 0.0, $
        vFwd: [0.0,0.0,1.0], $
        vUp: [0.0,1.0,0.0], $
        whichButton: 1B, $
        vertExag: vertExag, $
        ignoreCrashes: 0b, $
        bilinWarned: 0B, $
        xdim: xdim, $
        ydim: ydim, $
        record_to_filename: record_to_filename, $
        debug: keyword_set(debug) $
     }

    WIDGET_CONTROL, wTopBase, SET_UVALUE=sState, /NO_COPY

    ; save the id of the top level base to allow surfbase's event handler
    ; to access the state
    WIDGET_CONTROL, wSurfCntrlsTlb, SET_UVALUE=wTopBase

    WIDGET_CONTROL, wTopBase, SENSITIVE=1

    ;  Remove the starting up text.
    ;
    ;omTop->Remove, oText

    ;  Draw the screen.
    ;
    demo_draw, oWindow, oView, debug=keyword_set(debug)

    XMANAGER, 'd_flythruPopUp', wSurfCntrlsTlb, $
        EVENT_HANDLER='d_flythruSurfBaseEvent', $
        /JUST_REG, /NO_BLOCK

    XMANAGER, 'd_flythru', wTopBase, $
        EVENT_HANDLER='d_flythruEvent', $
        /NO_BLOCK, $
        CLEANUP='d_flythruCleanup'

end   ;   of  d_flythru



PRO test
z=bytarr(128,128)*0+1
z[31,31] = 100
imagefile='D:\AIMIS_3D\kelly.jpg';lucy2.jpg';kelly.jpg';lucy2.jpg';kelly.jpg'; gj.jpg
read_jpeg, imagefile, idata1
image=idata1

d_flythru_gj, SURFACE=z, image=image

END
