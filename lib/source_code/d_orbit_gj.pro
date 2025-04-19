; $Id: //depot/Release/ENVI53_IDL85/idl/idldir/examples/demo/demosrc/d_orbit.pro#1 $
;
;  Copyright (c) 1997-2015, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  FILE:
;       d_orbit.pro
;
;  CALLING SEQUENCE: d_orbit
;
;  PURPOSE:
;       Shows an orbiting satellite. The forcings are the Earth's
;       central body and J2.
;
;  MAJOR TOPICS: Visualization.
;
;  CATEGORY:
;       IDL Demo System
;
;  INTERNAL FUNCTIONS and PROCEDURES:
;       fun d_orbitGetRotation         - Rotation values of the satellite
;       fun d_orbitCreateSatellite     - create the satellite object
;       fun d_orbitFindE               - Compute the eccentric anomaly
;       fun d_orbitFindPosition        - Compute the satelitte's ephemeris.
;       pro d_orbitEvent               - Event handler
;       pro d_orbitCleanup             - Cleanup
;       pro d_orbit                    - Main procedure
;
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       pro orb__define.pro     - Create an orb object
;       pro trackball__define   - Create the trackball object
;       pro demo_gettips        - Read the tip file and create widgets
;       orbit.tip
;
;  REFERENCE: 1) Fundamentals of celestial mechanics (2nd ed.)
;                J.M.A. Danby
;                Willmann-Bell editor
;                ISBN 0-943396-20-4
;             2) Methods of orbit determination
;                P. R. Escobal
;                R. E. Kreiger publishing company
;                ISBN 0-88275-319-3
;
;  NAMED STRUCTURES:
;       none.
;
;  COMMON BLOCS:
;       none.
;
;  MODIFICATION HISTORY: Written by DAT,RSI,  July 1996 
;
;-
;---------------------------------------------------------------------------
;
;  Purpose:  Returns  3 rotation angle increment in degrees
;
;
function d_orbitGetRotation, $
    time         ; IN: time in sec passed the epoch

    angles = FLTARR(3)
    xAngle = 1.0d-3 * time
    yAngle = 4.0d-2 * time
    zAngle = 3.4d-3 * time

    angles = [ xAngle, yAngle, zAngle]

    RETURN, angles

end

;---------------------------------------------------------------------------
;
;  Purpose:  Create a satellite object and return it into a model
;
;
function d_orbitCreateSatellite

    mainModel = OBJ_NEW('IDLgrModel')

    ;  Create bloc1 : main satellite body in yellow.
    ;
    xp=[-0.07, 0.07, 0.07,-0.07, $
        -0.07, 0.07, 0.07,-0.07]
    yp=[-0.40,-0.40, 0.40, 0.40, $
        -0.40,-0.40, 0.40, 0.40]
    zp=[ 0.07, 0.07, 0.07, 0.07, $
        -0.07,-0.07,-0.07,-0.07]
    bloc1Vertices = fltarr(3,8)
    bloc1Vertices = [ [xp], [yp], [zp] ]
    bloc1Vertices = TRANSPOSE(bloc1Vertices)
    bloc1Mesh = [ [4,0,1,2,3], $
                  [4,1,5,6,2], $
                  [4,4,7,6,5], $
                  [4,0,3,7,4], $
                  [4,3,2,6,7], $
                  [4,0,4,5,1] ]
    bloc1List = FLTARR(3,24)
    bloc1NewMesh = bloc1Mesh
    j = 0
    for i = 0, 5 do begin
        bloc1List[0:2, i*4+0] = bloc1Vertices[0:2, bloc1Mesh[i*5+1]]
        bloc1List[0:2, i*4+1] = bloc1Vertices[0:2, bloc1Mesh[i*5+2]]
        bloc1List[0:2, i*4+2] = bloc1Vertices[0:2, bloc1Mesh[i*5+3]]
        bloc1List[0:2, i*4+3] = bloc1Vertices[0:2, bloc1Mesh[i*5+4]]
        bloc1NewMesh[*,i] = [4, j+0, j+1, j+2, j+3]
        j = j + 4
    endfor
    oBloc1 = OBJ_NEW('IDLgrPolygon', bloc1List, $
        POLYGONS=bloc1NewMesh, COLOR=[255,255,0] )
    mainModel->Add, oBloc1

    ;  Create bloc2 : solar array panel.
    ;
    xp=[0.20, 0.60, 0.60, 0.20, $
        0.20, 0.60, 0.60, 0.20]
    yp=[-0.20,-0.20, 0.20, 0.20, $
        -0.20,-0.20, 0.20, 0.20]
    zp=[ 0.02, 0.02, 0.02, 0.02, $
        -0.02,-0.02,-0.02,-0.02]
    bloc2Vertices = fltarr(3,8)
    bloc2Vertices = [ [xp], [yp], [zp] ]
    bloc2Vertices = TRANSPOSE(bloc2Vertices)
    bloc2Mesh = [ [4,0,1,2,3], $
                  [4,1,5,6,2], $
                  [4,4,7,6,5], $
                  [4,0,3,7,4], $
                  [4,3,2,6,7], $
                  [4,0,4,5,1] ]
    bloc2List = FLTARR(3,24)
    bloc2NewMesh = bloc1Mesh
    j = 0
    for i = 0, 5 do begin
        bloc2List[0:2, i*4+0] = bloc2Vertices[0:2, bloc2Mesh[i*5+1]]
        bloc2List[0:2, i*4+1] = bloc2Vertices[0:2, bloc2Mesh[i*5+2]]
        bloc2List[0:2, i*4+2] = bloc2Vertices[0:2, bloc2Mesh[i*5+3]]
        bloc2List[0:2, i*4+3] = bloc2Vertices[0:2, bloc2Mesh[i*5+4]]
        bloc2NewMesh[*,i] = [4, j+0, j+1, j+2, j+3]
        j = j + 4
    endfor
    oBloc2 = OBJ_NEW('IDLgrPolygon', bloc2List, $
        POLYGONS=bloc2NewMesh, COLOR=[100,100,255] )
    mainModel->Add, oBloc2

    ;  Create bloc3 : the other solar array panel.
    ;
    xp=[-0.20, -0.60, -0.60, -0.20, $
        -0.20, -0.60, -0.60, -0.20]
    yp=[-0.20,-0.20, 0.20, 0.20, $
        -0.20,-0.20, 0.20, 0.20]
    zp=[ 0.02, 0.02, 0.02, 0.02, $
        -0.02,-0.02,-0.02,-0.02]
    bloc3Vertices = fltarr(3,8)
    bloc3Vertices = [ [xp], [yp], [zp] ]
    bloc3Vertices = TRANSPOSE(bloc3Vertices)
    bloc3Mesh = [ [4,0,1,2,3], $
                  [4,1,5,6,2], $
                  [4,4,7,6,5], $
                  [4,0,3,7,4], $
                  [4,3,2,6,7], $
                  [4,0,4,5,1] ]
    bloc3List = FLTARR(3,24)
    bloc3NewMesh = bloc3Mesh
    j = 0
    for i = 0, 5 do begin
        bloc3List[0:2, i*4+0] = bloc3Vertices[0:2, bloc3Mesh[i*5+1]]
        bloc3List[0:2, i*4+1] = bloc3Vertices[0:2, bloc3Mesh[i*5+2]]
        bloc3List[0:2, i*4+2] = bloc3Vertices[0:2, bloc3Mesh[i*5+3]]
        bloc3List[0:2, i*4+3] = bloc3Vertices[0:2, bloc3Mesh[i*5+4]]
        bloc3NewMesh[*,i] = [4, j+0, j+1, j+2, j+3]
        j = j + 4
    endfor
    oBloc3 = OBJ_NEW('IDLgrPolygon', bloc3List, $
        POLYGONS=bloc3NewMesh, COLOR=[100,100,255] )
    mainModel->Add, oBloc3

    ;  Create the 4 white line connecting the main satellite
    ;  body to the solar array panels.
    ;
    x = [0.07, 0.20]
    y = [0.00, 0.18]
    z = [0.0, 0.0]
    oLine1 = OBJ_NEW('IDLgrPolyline', x, y, z, $
        COLOR=[255, 255, 255])
    mainModel->Add, oLine1

    x = [0.07, 0.20]
    y = [-0.00, -0.18]
    z = [0.0, 0.0]
    oLine2 = OBJ_NEW('IDLgrPolyline', x, y, z, $
        COLOR=[255, 255, 255])
    mainModel->Add, oLine2

    x = [-0.07, -0.20]
    y = [0.00, 0.18]
    z = [0.0, 0.0]
    oLine3 = OBJ_NEW('IDLgrPolyline', x, y, z, $
        COLOR=[255, 255, 255])
    mainModel->Add, oLine3

    x = [-0.07, -0.20]
    y = [-0.00, -0.18]
    z = [0.0, 0.0]
    oLine4 = OBJ_NEW('IDLgrPolyline', x, y, z, $
        COLOR=[255, 255, 255])
    mainModel->Add, oLine4

    mainModel->Scale, 0.5, 0.5, 0.5

    RETURN, mainModel

end

;--------------------------------------------------------------------------
;
;  Purpose:  Find the eccentric anomaly
;
function d_orbitFindE, $
    e, $        ; IN:  eccentricity (dimensionless)
    M, $        ; IN:  Mean anomaly
    tolerance   ; tolerance of iterative process

    M = M MOD (2.0*!DPI)
    diffx = 2.0 * tolerance
    xold = M
    while(diffx gt tolerance) do begin
        fx = xold - e*sin(xold) - M
        dfx = 1.0 - e*cos(xold)
        xnew = xold - fx/dfx
        diffx = xnew - xold
        xold = xnew
    endwhile
    RETURN, xnew
end

;--------------------------------------------------------------------------
;
;  Purpose:  Function returns the satellite position in
;            the geocentric coordinates system given the
;            time (in seconds) passed the epoch.
;
function d_orbitFindPosition, $
    sOrbit, $  ; IN:  structure containing the orbital elements
    time       ; IN:  time in seconds passed the epoch

    rad = !DPI / 180.0
    sini = sin(sOrbit.i * rad)
    cosi = cos(sOrbit.i * rad)

    ;  Compute the 3 orbital elements affected by J2 and
    ;  compute their new values.
    ;
    M = sOrbit.M0*rad + sOrbit.dMdt * time
    omega = sOrbit.omega0*rad + sOrbit.dOmegadt * time
    w = sOrbit.w0*rad + sOrbit.dwdt * time

    ;  Initialize the p and q vector for the
    ;  transoformation matrix form the orbital plane to
    ;  earth-fixed coordinate system of coordinates.
    ;
    p = FLTARR(3,3)
    q = FLTARR(3)
    cosOmega = cos(omega)
    sinOmega = sin(omega)
    cosw = cos(w)
    sinw = sin(w)
    p[0,0] = cosw*cosOmega - sinw*sinOmega*cosi
    p[0,1] = cosw*sinOmega + sinw*cosOmega*cosi
    p[0,2] = sinw*sini
    p[1,0] = - sinw*cosOmega - cosw*sinOmega*cosi
    p[1,1] = - sinw*sinOmega + cosw*cosOmega*cosi
    p[1,2] = cosw*sini
    p[2,0] = sinOmega*sini
    p[2,1] = -cosOmega*sini
    p[2,2] = cosi

    ;  Find the eccentric anomaly using the 'Newton' method.
    ;
    tolerance = 1.0d-5
    EAnomaly = d_orbitFindE(sOrbit.e, M, tolerance)

    ;  Find the cartesian coordinates in the orbital plane.
    ;
    pp = FLTARR(3)
    pp[0] = sOrbit.a*(cos(eAnomaly) - sOrbit.e)
    pp[1] = sOrbit.a * SQRT(1.- sOrbit.e*sOrbit.e) * sin(eAnomaly)
    pp[2] = 0.0

    ;  Find the satellite position in Earth fixed coordinate system
    ;  Note : x axis points toward the  ernal equinox
    ;         z axis points toward the Earth body-fixed z axis
    ;         y axis is orthonormal to x and z
    ;
    xf = FLTARR(3)
    xf = pp # p

    RETURN, xf
end

;--------------------------------------------------------------------------
;
;  Purpose:  Handle the events
;
pro d_orbitEvent, $
     sEvent       ; IN:  event structure

    ;  Quit the application using the close box.
    ;
    if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ $
        'WIDGET_KILL_REQUEST') then begin
        WIDGET_CONTROL, sEvent.top, /DESTROY
        RETURN
    endif

    WIDGET_CONTROL, sEvent.id, GET_UVALUE=uValue

    case uValue of
        ;  Start the animation.
        ;
        'START' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wRightBase, TIMER=sState.timer
            sState.stopFlag = 0
            WIDGET_CONTROL, sState.wStopButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wResetButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wStartButton, SENSITIVE=0
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        ;  Handle the widget timer event. This event create another
        ;  widget timer event, hence creating a loop. The loop is
        ;  broken by pushing the 'STOP' button.
        ;
        'TIMER' : begin

            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

            ;  Reset if time greater than 1 d7 seconds (115 days).
            ;
            if (sState.time GT 1.0d7) then begin
                WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
                WIDGET_CONTROL, sState.wStopButton, SENSITIVE=0
                WIDGET_CONTROL, sState.wResetButton, SENSITIVE=0
                WIDGET_CONTROL, sState.wStartButton, SENSITIVE=1
                sState.oEarthRotationModel->SetProperty, $
                    TRANSFORM=sState.initEarthTM
                sState.oSatRotationModel->SetProperty, $
                    TRANSFORM=sState.initSatTM
                sState.oTopModel->SetProperty, $
                    TRANSFORM=sState.initTopTM
                sState.oSky->Rotate,[0,0,1], -sState.skyAngle
                delta = sState.initSatPosition-sState.Position
                sState.oSatModel->Translate, delta[1], delta[2], delta[0]
                sState.oWindow->draw, sState.oView
                sState.time = 0.0
                sState.position = sState.initSatPosition
                WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
                RETURN
            endif

            ;  Process the animation only if the 'STOP' button has
            ;  not been pushed.
            ;
            if (sState.stopFlag EQ 0) then begin

                ;  Increment the time.
                ;
                sState.time = sState.time + sState.timeIncrement

                ;  Given the time (in sec.) passed the epoch, compute
                ;  The satellite ephemeris (x, y, z, position).
                ;
                newPosition = d_orbitFindPosition(sState.sOrbit, sState.time)
                newPosition = newPosition/6378137.0
                delta = newPosition - sState.position
                sState.Position = newPosition


                ;  Translate the satellite.
                ;  Notice the axis correspondance :
                ;  orbit coord. system          graphic cood. system
                ;
                ;          x                              z
                ;          y                              x
                ;          z                              y
                ;  In order to make this transofrmation, the translation
                ;  is done as follows :
                ;
                sState.oSatModel->Translate, delta[1], delta[2], delta[0]

                ;  Rotate the satellite.
                ;
                angles = d_orbitGetRotation(sState.timeIncrement)
                sState.oSatRotationModel->Rotate, [1,0,0], angles[1]
                sState.oSatRotationModel->Rotate, [0,1,0], angles[1]
                sState.oSatRotationModel->Rotate, [0,0,1], angles[2]

                ;  Rotate the earth and keep track of the sky rotation
                ;
                earthRotationRate = 0.004178074   ; degrees per sec
                earthAngle = earthRotationRate*sState.timeIncrement
                sState.oEarthRotationModel->Rotate, [0,1,0], earthAngle
                sState.oSky->Rotate,[0,0,1], -earthAngle*1.5
                sState.skyAngle = sState.skyAngle - (earthAngle*1.5)
                sState.skyAngle = (sState.skyAngle MOD 360)

                ;  Draw the new frame.
                ;
                sState.oWindow->Draw, sState.oView

                ;  Create a timer event to form the animation loop
                ;
                WIDGET_CONTROL, sState.wRightBase, TIMER=sState.timer

            endif

            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
         end

        ;  Handle events of the drawing area.
        ;
        'DRAWING' : begin

            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY

            ;  Expose.
            ;
            if (sEvent.type EQ 4) then begin
                sState.oWindow->draw, sState.oView
                WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
                RETURN
            endif

            ;  Handle trackball update.
            ;
            bHaveTransform = sState.oTrack->Update(sEvent, TRANSFORM=qmat )
            if (bHaveTransform NE 0) then begin
                sState.oTopModel->GetProperty, TRANSFORM=t
                mt = t # qmat
                sState.oTopModel->SetProperty,TRANSFORM=mt
            endif

            ;  Handle button press.
            ;
            if (sEvent.type EQ 0) then begin
                sState.btndown = 1b
                sState.oWindow->SetProperty, QUALITY=0
                WIDGET_CONTROL,sState.wDraw, /DRAW_MOTION
            endif

            ;  Handle button motion.
            ;
            if (sEvent.type EQ 2) then begin
                if (bHaveTransform) then $
                    sState.oWindow->Draw, sState.oView
            endif

            ;  Handle button release
            ;
            if (sEvent.type EQ 1) then begin
                if (sState.btndown EQ 1b) then begin
                    sState.oWindow->SetProperty, QUALITY=2
                    sState.oWindow->Draw, sState.oview
                endif
                sState.btndown = 0b
                WIDGET_CONTROL, sState.wDraw, DRAW_MOTION=0
            endif

            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        ;  Stop the animation.
        ;
        'STOP' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            sState.stopFlag = 1
            WIDGET_CONTROL, sState.wStopButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wResetButton, SENSITIVE=1
            WIDGET_CONTROL, sState.wStartButton, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
            RETURN
        end

        ;  Reset the initial state (time and orientation of the objects.)
        ;
        'RESET' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=0
            WIDGET_CONTROL, sState.wStopButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wResetButton, SENSITIVE=0
            WIDGET_CONTROL, sState.wStartButton, SENSITIVE=1
            sState.oEarthRotationModel->SetProperty, $
                TRANSFORM=sState.initEarthTM
            sState.oSatRotationModel->SetProperty, $
                TRANSFORM=sState.initSatTM
            sState.oTopModel->SetProperty, $
                TRANSFORM=sState.initTopTM
            sState.oSky->Rotate,[0,0,1], -sState.skyAngle
            delta = sState.initSatPosition-sState.Position
            sState.oSatModel->Translate, delta[1], delta[2], delta[0]
            sState.oWindow->draw, sState.oView
            sState.time = 0.0
            sState.position = sState.initSatPosition
            WIDGET_CONTROL, sState.wTopBase, SENSITIVE=1
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
            RETURN
        end

        ;  Scale the satellite.
        ;
        'SCALING' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wScalingSlider, GET_VALUE=scale

            scale = 0.75 + FLOAT(scale) / 100.0
            scalep = scale*100.0
            scalingString = STRING(scalep, FORMAT='(f5.1)') + ' %'
            WIDGET_CONTROL, sState.wScalingLabel2, $
                SET_VALUE=scalingString

            transform = [[scale, 0, 0, 0.0], [0, scale, 0, 0.0], $
            [0, 0, scale, 0.0], [0, 0, 0, 1]]
            sState.oSatScalingModel->SetProperty, TRANSFORM=transform
            WIDGET_CONTROL, sState.wTopBase, /HOURGLASS
            sState.oWindow->Draw, sState.oView
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        ;  Determine the time step between each redraw.
        ;  A larger time step will make the satellite to
        ;  translate at a faster rate but the  animation
        ;  will be less 'smoother'.
        ;
        'STEP' : begin
            WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
            WIDGET_CONTROL, sState.wStepSlider, GET_VALUE=timeIncrement
            sState.timeIncrement = timeIncrement
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        end

        ;  Quit this application.
        ;
        'QUIT' : begin
            WIDGET_CONTROL, sEvent.top, /DESTROY
        end

        ELSE :   ;  do nothing

    endcase
end

;--------------------------------------------------------------------------
;
;  Purpose:  Destroy the top objects and restore the previous
;            color table
;
pro d_orbitCleanup, $
    wTopBase        ;  IN: top level base identifier

    WIDGET_CONTROL, wTopBase, GET_UVALUE=sState, /NO_COPY

    ;  Destroy the top objects
    ;
    OBJ_DESTROY, sState.oView
    OBJ_DESTROY, sState.oTrack
    OBJ_DESTROY, sState.oText
    OBJ_DESTROY, sState.oFont
    OBJ_DESTROY, sState.oContainer

    ;  Restore the color table
    ;
    TVLCT, sState.colorTable

    if WIDGET_INFO(sState.groupBase, /VALID_ID) then $
        WIDGET_CONTROL, sState.groupBase, /MAP

end   ;  of d_orbitCleanup

;--------------------------------------------------------------------------
;
;  Purpose:  Show a satellite orbiting the Earth
;
pro d_orbit_gj, $
    GROUP=group, $         ; IN: (opt) semi-major axis (meters)
    AXIS=a, $              ; IN: (opt) semi-major axis (meters)
    INCLINATION=i, $       ; IN: (opt) inclination (deg.)
    ECCENTRICITY=e, $      ; IN: (opt) eccentricity
    MEANANOMALY=M0, $      ; IN: (opt) mean anomaly (deg.)
    NODE=omega0, $         ; IN: (opt) longitude of ascending node (deg.)
    PERIGEE=w0, $          ; IN: (opt) argument of perigee (deg.)
    SELECTION=select, $    ; IN: (opt) 1 = circular orbit
                           ;           2 = elliptical orbit (default)
                           ;           3 = polar orbit
    RECORD_TO_FILENAME=record_to_filename, $
    APPTLB = appTLB        ; OUT: (opt) TLB of this application


    ;  Check the validity of the group identifier.
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

    ;  The default values of the orbital elements are those of
    ;  the elliptical orbit
    ;
    ;  NOTE :  In this code, the value range of the orbital elements are :
    ;
    ;  a     :  9200 to 11000 km
    ;  e     :  0  to 0.3
    ;  i     :  0  to 180
    ;  M     :  0 to 360
    ;  omega :  0 to 360
    ;  w     :  0 to 360
    ;

    ;  Initialize the Keplerian orbit elements
    ;  at the epoch time here
    ;

    ;  Semi-major axis (in meters).
    ;
    if (N_ELEMENTS(a) NE 0) then begin
        a = FLOAT(a)
        if ( (a LT 9200000.0) OR (a GT 11000000.0) ) then $
            a = 9200000.0d0
    endif else begin
        a = 9200000.0d0
    endelse

    ;  Eccentricity (dimensionless).
    ;
    if (N_ELEMENTS(e) NE 0) then begin
        e = FLOAT(e)
        if ( (e LT 0.0) OR (e GT 0.3) ) then $
            e = 0.3
    endif else begin
        e = 0.3
    endelse

    ;  Inclination (in degrees).
    ;
    if (N_ELEMENTS(i) NE 0) then begin
        i = FLOAT(i)
        if ( (i LT 0.0) OR (i GT 180.0) ) then $
            i = 75.0
    endif else begin
        i = 75.0
    endelse

    ;  Mean anomaly (in degrees).
    ;
    if (N_ELEMENTS(M0) NE 0) then begin
        M0 = FLOAT(M0)
        if ( (M0 LT 0.0) OR (M0 GT 360.0) ) then $
            M0 = 38.0
    endif else begin
        M0 = 38.0
    endelse

    ;  Longitude of ascending node (in degrees).
    ;
    if (N_ELEMENTS(omega0) NE 0) then begin
        omega0 = FLOAT(omega0)
        if ( (omega0 LT 0.0) OR (omega0 GT 360.0) ) then $
            omega0 = 64.0
    endif else begin
        omega0 = 64.0
    endelse

    ;  Argument of perigee (in degrees).
    ;
    if (N_ELEMENTS(w0) NE 0) then begin
        w0 = FLOAT(w0)
        if ( (w0 LT 0.0) OR (w0 GT 360.0) ) then $
            w0 = 28.0
    endif else begin
        w0 = 28.0
    endelse

    ;  Orbit selection
    ;
    ;  1 = circular orbit
    ;  2 = elliptical orbit (default)
    ;  3 = polar orbit
    ;
    if (N_ELEMENTS(select) NE 0) then begin
        if (select EQ 1) then begin
            a = 9200000.0d0
            e = 0.0
            i = 85.0
            M0 = 38.0
            omega0 = 90.0
            w0 = 28.0
        endif else if (select EQ 3) then begin
            a = 9200000.0d0
            e = 0.3
            i = 90.0
            M0 = 0.0
            omega0 = 180.0d0
            w0 = 90.0
        endif else begin
            a = 9200000.0d0
            e = 0.3
            i = 75.0
            M0 = 38.0
            omega0 = 64.0
            w0 = 28.0
        endelse
    endif

    ;  Get the current color vectors to restore
    ;  when this application is exited.
    TVLCT, savedR, savedG, savedB, /GET

    ; Build color table from color vectors
    ;
    colorTable = [[savedR],[savedG],[savedB]]

    ;  Set up dimensions of the drawing (viewing) area.
    ;
    DEVICE, GET_SCREEN_SIZE=screenSize
    xdim = screenSize[0]*0.6
    ydim = xdim*0.8

    ;  Create widgets.
    ;
    if (N_ELEMENTS(group) NE 0) then begin
        wTopBase = WIDGET_BASE(/COLUMN, $
            TLB_FRAME_ATTR=1, MBAR=barBase, $
            GROUP_LEADER=group, $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            TITLE="Orbiting Satellite")
    endif else begin
        wTopBase = WIDGET_BASE(/COLUMN, $
            TLB_FRAME_ATTR=1, MBAR=barBase, $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            TITLE="Orbiting Satellite")
    endelse

        ;  Create the menu bar. It contains the file,
        ;  edit, and help buttons.
        ;
        wFileButton = WIDGET_BUTTON(barBase, VALUE = 'File', /MENU)

            wQuitButton = WIDGET_BUTTON(wFileButton, VALUE='Quit', $
                UVALUE='QUIT')

        ;  Create a sub base of the top base (wBase).
        ;
        wSubBase = WIDGET_BASE(wTopBase, COLUMN=2)

            ;  Create the left Base that contains the functionality buttons.
            ;  Here the only button is to animate the object.
            ;
            wLeftbase = WIDGET_BASE(wSubBase,/BASE_ALIGN_CENTER, $
                COLUMN=1)

                wAnimationBase = WIDGET_BASE(wLeftBase, $
                    /COLUMN, YPAD=10)

                    wAnimationLabel = WIDGET_LABEL(wAnimationBase, $
                        VALUE='Animation')

                    wStartButton = WIDGET_BUTTON(wAnimationBase, $
                        VALUE="Start", UVALUE='START')

                    wStopButton = WIDGET_BUTTON(wAnimationBase, $
                        VALUE="Stop", UVALUE='STOP')

                    wResetButton = WIDGET_BUTTON(wAnimationBase, $
                        VALUE="Reset", UVALUE='RESET')

                wScalingBase = WIDGET_BASE(wLeftBase, $
                    /COLUMN, YPAD=10)

                    wScalingLabel1 = WIDGET_LABEL(wScalingBase, $
                        VALUE='Satellite Scaling')

                    percent = 100
                    scalingString = STRING(percent, FORMAT='(f5.1)') + ' %'

                    wScalingLabel2 = WIDGET_LABEL(wScalingBase, $
                        VALUE=scalingString)

                    wScalingSlider = WIDGET_SLIDER(wScalingBase, $
                        MINIMUM=0, $
                        MAXIMUM=50, VALUE=25, $
                        /SUPPRESS_VALUE, $
                        UVALUE='SCALING')


                wStepBase = WIDGET_BASE(wLeftBase, $
                    /COLUMN, YPAD=10)

                    wStepLabel1 = WIDGET_LABEL(wStepBase, $
                        VALUE='Animation Step')

                    wStepLabel2 = WIDGET_LABEL(wStepBase, $
                        VALUE='in Seconds')

                    wStepSlider = WIDGET_SLIDER(wStepBase, $
                        MINIMUM=25, $
                        MAXIMUM=100, VALUE=100, $
                        UVALUE='STEP')

            ;  Create the right Base that has the drawing area.
            ;
            wRightbase = WIDGET_BASE(wSubBase, COLUMN=1, $
                UVALUE='TIMER')

                wDraw = WIDGET_DRAW(wRightBase, $
                    GRAPHICS_LEVEL=2, $
                    XSIZE=xdim, YSIZE=ydim, /BUTTON_EVENTS, $
                    UVALUE='DRAWING', RETAIN=0, /EXPOSE_EVENT)

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
    sText = demo_getTips(demo_filepath('orbit.tip', $
                             SUBDIR=['examples','demo', 'demotext']), $
                         wTopBase, $
                         wStatusBase)

    WIDGET_CONTROL, wTopBase, SENSITIVE=0

    ;  Grab the window id of the drawable.
    ;
    WIDGET_CONTROL, wDraw, GET_VALUE=oWindow

    ;  Compute viewplane rectangle based on aspect ratio.
    ;
    aspect = float(xdim)/float(ydim)
    myview = [-1.0,-1.0,2.0,2.0]
    myview = [-1.5,-1.5,3.0,3.0]
    if (aspect > 1) then begin
        myview[0] = myview[0] - ((aspect-1.0)*myview[2])/2.0
        myview[2] = myview[2] * aspect
    endif else begin
        myview[1] = myview[1] - (((1.0/aspect)-1.0)*myview[3])/2.0
        myview[3] = myview[3] * aspect
    endelse

    ; Create view object.
    ;
    oView = OBJ_NEW('idlgrview', PROJECTION=2, EYE=4, $
        ZCLIP=[2.0,-2.0], VIEW=myview, COLOR=[0, 0, 0])

    ;  Create a centered starting up text.
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


    ; Create models and its tre structure.
    ;
    oTopModel = obj_new('idlgrmodel')
        oEarthModel = OBJ_NEW('idlgrmodel')
        oEarthRotationModel = OBJ_NEW('idlgrmodel')
        oEarthScalingModel = OBJ_NEW('idlgrmodel')
        oSatModel = OBJ_NEW('idlgrmodel')
        oSatRotationModel = OBJ_NEW('idlgrmodel')
        oSatScalingModel = OBJ_NEW('idlgrmodel')

    ;  Place the model in the view.
    ;
    oView->Add, oTopModel

    oTopModel->Add, oText
    oWindow->Draw, oView

    oTopModel->Add, oEarthModel
    oTopModel->Add, oSatModel
    oEarthModel->Add, oEarthRotationModel
    oEarthRotationModel->Add, oEarthScalingModel
    oSatModel->Add, oSatRotationModel
    oSatRotationModel->Add, oSatScalingModel
    scale = 2.0
    oEarthScalingModel->Scale, scale, scale, scale

    ;  Scale the top model to fit the viewing area.
    ;
    sct = 0.7
    oTopModel->Scale, sct, sct, sct

    ;  Create lights
    ;
    oLIght1 = OBJ_NEW('IDLgrLight', DIRECTION=[2, 2, 5], $
        TYPE=2, INTENSITY=0.25)

    oLIght2 = OBJ_NEW('IDLgrLight', $
        TYPE=0, INTENSITY=0.7)

    oTopModel->Add, oLight1
    oTopModel->Add, oLight2

    ;  Create the earth
    ;
    oImageArray = OBJARR(2)
    READ_JPEG, demo_filepath("earth.jpg", $
        SUBDIR=['examples','demo','demodata']), $
        image, TRUE=1
      ;GJ 2019/8/12, test texture mapping
      imagefile='D:\AIMIS_3D\gj.jpg'
      read_jpeg, imagefile, image3
    oImageArray[0] = OBJ_NEW('IDLgrImage', image3, HIDE=1)
    READ_JPEG, demo_filepath("cloud.jpg", $
        SUBDIR=['examples','demo','demodata']), $
        image2, TRUE=1
    sizeImage = SIZE(image2)
    rgba = BYTARR(4, sizeImage[2], sizeImage[3])
    rgba[0, *, *] = image2[0, *, *]
    rgba[1, *, *] = image2[1, *, *]
    rgba[2, *, *] = image2[2, *, *]
    rgba[3, *, *] = (FLOAT(image2[0, *, *]) $
        + FLOAT(image2[1, *, *]) + FLOAT(image2[2, *, *]) )/3.0
    oImageArray[1] = OBJ_NEW('IDLgrImage', rgba, HIDE=1)

    oEarthScalingModel->Add, oImageArray[0]
    oEarthScalingModel->Add, oImageArray[1]

    oSky = OBJ_NEW('orb', COLOR=[255, 255, 255], RADIUS=0.26, $
        DENSITY=0.8, /TEX_COORDS, TEXTURE_MAP=oImageArray[0])

    oPlanet = OBJ_NEW('orb', COLOR=[255, 255, 255], RADIUS=0.25, $
        DENSITY=0.8, /TEX_COORDS, TEXTURE_MAP=oImageArray[1])

    oEarthScalingModel->Add, oPlanet
    oEarthScalingModel->Add, oSky
    oEarthRotationModel->Rotate, [1,0,0], 90
    oEarthRotationModel->Rotate, [0,1,0], 90
    oEarthRotationModel->GetProperty, TRANSFORM=initEarthTM
    oSatRotationModel->GetProperty, TRANSFORM=initSatTM
    oTopModel->GetProperty, TRANSFORM=initTopTM

    ;  Invert texture coordinates.
    ;
    oPlanet->GetProperty, POBJ=p
    p->GetProperty, TEXTURE_COORD=t
    
    oPlanet->getproperty, POLYGONS=conn
    oPlanet->getproperty, DATA=verts

    t[0,*] = 1.0 - t[0,*]
    p->SetProperty, TEXTURE_COORD=t

    oWindow->SetProperty, QUALITY=2

    ;  Create a satellite object placed into a model,
    ;  return that model.
    ;
    oSatelliteModel = d_orbitCreateSatellite()
    oSatScalingModel->Add, oSatelliteModel

    ;  Initialize the Keplerian orbit elements
    ;  at the epoch time here
    ;
;    a =  9000000.0     ; semi-major axis in meters
;    e = 0.0            ; eccentricity (dimensionless)
;    i =  80.0          ; inclination in degrees
;    M0 =  38.0         ; mean anomaly at the epoch time (degree)
;    omega0 = 90.0     ; Longitude of ascending node at the epoch (degree)
;    w0 = 28.0          ; argument of perigee (degree)

    rad = !DPI / 180.0
    sini = sin(i * rad)
    cosi = cos(i * rad)

    ;  Compute the secular orbit perturbations due to J2.
    ;
    J2 = 1.0822d-3
    u = 3.9860044d14
    n = SQRT(u / a^3)
    kterm = 3.0*n * J2 / $
        (2.0* a * a * (1.0 - e*e)^2)
    dwdt = -kterm * ( 2.5*sini^2 -2.0 )
    dOmegadt = -kterm * cosi
    dMdt = n - kterm * SQRT(1.0 - e*e) * $
        ( 1.5*sini^2 - 1.0 )

    ;  Place the orbital elements and the secular rate due to J2
    ;  into a structure.
    ;
    sOrbit = { $
        A: a, $                        ; orbital elements , semi-major axis
        E: e, $                        ; eccentricity
        i: i, $                        ; inclination
        M0:m0, $                       ; mean anomaly
        Omega0: omega0, $              ; longitude of ascending node
        w0: w0, $                      ; argument of perigee
        DMdt: dMdt, $                  ; secular rate due to J2: mean anomaly
        dOmegadt: dOmegadt, $          ; longitude of ascending node
        dwdt: dwdt $                   ; argument of perigee
    }

    ;  Translate the satellite and place it to its
    ;  location at the epoch time
    ;
    position = d_orbitFindPosition(sOrbit, 0.0d0)
    newPosition = position/6378137.0

    ;  The orbit coordinate system does not correspond to
    ;  the graphic coordinates system. The linear transformation
    ;  is done by translating the satellite (oSatModel) as
    ;  follows :
    ;    orbit coord. system    graphic coord. system
    ;         x                      z
    ;         y                      x
    ;         z                      y
    ;
    oSatModel->Translate, newPosition[1], newPosition[2], newPosition[0]

    ;  Add the trackball object for interactive change
    ;  of the scene orientation
    ;
    oTrack = OBJ_NEW('Trackball', [xdim/2.0, ydim/2.0], xdim/2.0)

    ocontainer = OBJ_NEW('IDLgrContainer')
    oContainer->Add, oView
    oContainer->Add, oTrack

    sState = { $
        ColorTable: colorTable, $            ; Color table to restore
        Btndown:0b, $                        ; 0= not presses, pressed otherwise
        SOrbit: SOrbit, $                    ; Orbital elements structure
        InitEarthTM: initEarthTM, $          ; Initial Orientation of Earth
        InitSatTM: initSatTM, $              ; Initial Satellite transformation
        InitTopTM: initTopTM, $              ; Initial transformationo of top model
        InitSatPosition: newPosition, $      ; Initial Position of the satellite
        StopFlag: 0, $                       ; 0= animation not stopped, 1=stopped
        Position: newPosition, $             ; Updated position of the satellite
        Timer: 0.1, $                        ; Time delay between frame
        Time: 0.00d0, $                      ; Accumulated time (sec.)
        TimeIncrement: 50.0, $               ; Time between frames (sec.)
        SkyAngle: 0.0, $                     ; rotation angle of the clouds (sky)
        WTopBase: wTopBase, $                ; Top level base IDs
        WRightBase: wRightBase, $            ; Right base ID
        WStartButton: wStartButton, $        ; Animation buttons IDs
        WResetButton: wResetButton, $
        WStopButton: wStopButton, $
        WScalingSlider: wScalingSlider, $    ; Sliders IDs
        WScalingLabel2: wScalingLabel2, $    ; Sliders IDs
        WStepSlider: wStepSlider, $
        WDraw: wDraw, $                      ; Wdiget draw ID
        OView: oView, $                      ; View object
        OWindow: oWindow, $                  ; Window object
        OSatModel: oSatModel, $              ; Satellite models
        OSatRotationModel: oSatRotationModel, $
        OSatScalingModel: oSatScalingModel, $
        OEarthRotationModel: oEarthRotationModel, $ ; Earth models
        OSky: oSky, $                        ; Clouds image object
        OTrack: oTrack, $                    ; Trackball object
        OFont: oFont, $                      ; Font object
        OText: oText, $                      ; Text object
        OContainer: oContainer, $            ; Container object
        OTopModel:oTopModel, $               ; Top model
        groupBase: groupBase $               ; Base of Group Leader
    }     

    WIDGET_CONTROL, wTopBase, SET_UVALUE=sState, /NO_COPY

    WIDGET_CONTROL, wStartButton, SENSITIVE=0
    WIDGET_CONTROL, wResetButton, SENSITIVE=0

    ;  Draw the first view.
    ;
    WIDGET_CONTROL, wTopBase, /HOURGLASS
    oTopModel->Remove, oText
    oWindow->Draw, oView
    
    WIDGET_CONTROL, wTopBase, SENSITIVE=1

    ;  Animate at start up.
    ;
    pseudoEvent = { $
       ID: wRightBase, $
       TOP: wTopBase, $
       HANDLER: wTopBase $
    }

    d_orbitEvent, pseudoEvent

    XMANAGER, 'd_orbit', wTopBase, EVENT_HANDLER='d_orbitEvent', $
        /NO_BLOCK, $
        CLEANUP='d_orbitCleanup'

end
