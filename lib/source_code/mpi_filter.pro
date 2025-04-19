; $Id: //depot/Release/ENVI51_IDL83/idl/idldir/examples/demo/demosrc/mpi_filter.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  FILE:
;       mpi_filter.pro
;
;  CALLING SEQUENCE: mpi_filter
;
;  PURPOSE:
;       Shows the MPI filtering technique.
;
;  MAJOR TOPICS: Data analysis and plotting.
;
;  CATEGORY:
;       IDL Demo System
;
;  INTERNAL FUNCTIONS and PROCEDURES:
;       pro mpi_filter                - Main procedure
;       pro mpi_filterSetSliders      - Create the sliders
;       pro mpi_filterLoadData        - Read and load new data set
;       pro mpi_filterResetFilter     - Compute the new filter function
;       pro mpi_filterColorMapCallback - Redraw after colormap change
;       pro mpi_filterResetView       - Show the updated plots
;       pro mpi_filterEvent           - Event handler
;       pro mpi_filterCleanup         - Event handler
;       pro mpi_filterGetImage        - Load a new image
;       pro mpi_filterGetImageEvent   - Event handler for mpi_filterGetImage
;       pro mpi_filterGetImageCleanup - Event handler for mpi_filterGetImage
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       pro demo_gettips        - Read the tip file and create widgets
;       filter.tip
;       MPI_data.txt
;       abnorm.dat
;       alie.dat
;       cereb.dat
;       chirp.dat
;       damp_sn.dat
;       damp_sn2.dat
;       dirty_sine.dat
;       galaxy.dat
;       jet.dat
;
;  REFERENCE: IDL Reference Guide, IDL User's Guide
;
;  NAMED STRUCTURES:
;       none.
;
;  COMMON BLOCS:
;       filterCommon
;
;  MODIFICATION HISTORY:  Written by:  DS, RSI,
;                         Modified by : DAT  12/96, Change GUI
;                         Modified by: GJ, 2021/7/27, add MPI relaxation filter
;-
; -----------------------------------------------------------------------------
;
;  PURPOSE:  This routine sets the sliders that are data dependent.
;
pro mpi_filterSetSliders, $
    sliderBase     ; IN: slider base identifier

    COMMON filterCommon, filterbase, image, imagetrans, view, $
      viewsize, order, cutoff, bandwidth, filter, filteredImage, $
      filteredSpectrum, nviews, bandtype, wcutoffslider, $
      wbandwidthslider, filtertype, $
      wFunctionButton,wOrderButton,wTypeButton, $
      wImageButton, wIntensityButton, wCrossButton, wShadedButton, $
      currentRed, distfun

    ;  Determine the dimensions of the selected data.
    ;
    dimensions = SIZE(image)

    if (N_PARAMS() NE 0) then begin
   
        ;  Define the width of the sliders as a function of hardware display
        ;  size (VIEWSIZE) and number of plot windows (NVIEWS). VIEWSIZE is
        ;  set at 192 for horizontal hardware display sizes less than 1155
        ;  and 256 otherwise.
        ;
        sliderWidth = (NVIEWS * VIEWSIZE) + ((NVIEWS + 1) * 20) - 670

        ;  Define a floating-point slider to adjust the frequency cutoff value.
        ;
        cutoff = 20
;        wCutOffSlider = WIDGET_SLIDER(sliderBase, VALUE=cutoff, $
;            MINIMUM=1, MAXIMUM=dimensions[1]/2., $
;            UVALUE="FILTERCUTOFF", TITLE="Frequency Cutoff")
                       
        ;  Define a floating-point slider to adjust the bandwidth value.
        ;
        bandWidth = 10
;        wBandWidthSlider = WIDGET_SLIDER(sliderBase, VALUE=bandWidth, $
;            MINIMUM=1, MAXIMUM=dimensions[1]/4., $
;            UVALUE="FILTERBANDWIDTH", TITLE="Bandwidth")

        ;  Desensitize the frequency bandwidth slider.
        ;
;        WIDGET_CONTROL, wBandWidthSlider, SENSITIVE=0

    endif else begin

        ;  Initialize the frequency bandwidth value.
        ;
 ;       WIDGET_CONTROL, wBandWidthSlider, SET_VALUE=10

    endelse

end        ;   of  mpi_filterSetSliders 

;+
; NAME:
;	mpi_filterGetImage
;
; PURPOSE:
;	This procedure retrieves an image file from the examples/data directory.
;
; CATEGORY:
;	Examples General.
;
; CALLING SEQUENCE:
;	mpi_filterGetImage, newdata
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	newdata:The image data is returned in this variable.
;
; COMMON BLOCKS:
;	GF:
;
; RESTRICTIONS:
;	This procedure can only be used to reference the example data.
;
; PROCEDURE:
;
; EXAMPLE:
;	mpi_filterGetImage, frames, DESCRIPTION=description, $
;		DIMENSIONS=dimensions, SUBDIR = ["examples","data"], $
;		/THREE_DIM, TITLE="Select Animation", /ASSOC_IT
;
; MODIFICATION HISTORY:  Written by:  WSO, RSI, January 1995
;-


;------------------------------------------------------
;   procedure mpi_filterGetImageEvent
;------------------------------------------------------

pro mpi_filterGetImageEvent, event

COMMON GF, selection

   WIDGET_CONTROL, event.id, GET_UVALUE = selected

   CASE selected OF

     "FILELST": BEGIN
         selection = event.index
           ; If double clicking on a list item - emulate OK button
         IF event.clicks EQ 2 THEN BEGIN
            selection = selection + 1
            WIDGET_CONTROL, event.top, SENSITIVE = 0
            WIDGET_CONTROL, event.top, /DESTROY
         ENDIF
         END

     "OK": BEGIN
         selection = selection + 1
         WIDGET_CONTROL, event.top, SENSITIVE = 0
         WIDGET_CONTROL, event.top, /DESTROY
         END

     "CANCEL": BEGIN
         WIDGET_CONTROL, event.top, /DESTROY
         selection = 0
         END

   ENDCASE

END
;---------- end of procedure mpi_filterGetImageEvent ------------

;-----------------------------------------------------------
;
;  PUPOSE Cleanup procedure.
;
pro mpi_filterGetImageCleanup, $
    tlb

    WIDGET_CONTROL, tlb, GET_UVALUE=sState, /NO_COPY

    ;  Restore the previous color table.
    ;
    TVLCT, sState.colorTable

end   ;   of mpi_filterGetImageCleanup


;------------------------------------------------------
;   procedure mpi_filterGetImage
;------------------------------------------------------

PRO mpi_filterGetImage, newdata, DESCRIPTION = description, $
      GROUP=group, $
      DIMENSIONS = dimensions, ONE_DIM = one_dim, TWO_DIM = two_dim, $
      THREE_DIM = three_dim, TITLE = title, FILENAME = filename, $
      OFILENAME = ofilename, ASSOC_IT = assoc_it, $
      SUBDIRECTORY = subdirectory

   COMMON GF, selection

   if n_elements(group) eq 0 then group = 0L

   ;  Get the current color vectors to restore when this application is exited.
   ;
   TVLCT, savedR, savedG, savedB, /GET

   ;  Build color table from color vectors.
   ;
   colorTable = [[savedR],[savedG],[savedB]]

   name = ''
   dim = LONARR(3)
   des = ''
   del = ''
   numfiles = 0L
   one_mask = 1
   two_mask = 2
   three_mask = 4

   newdata = 0
   description = 0
   dimensions = 0

   IF (KEYWORD_SET(filename)) THEN BEGIN

      OPENR, unit, DEMO_FILEPATH("MPI_data.txt", SUBDIRECTORY=subdirectory), $
         /GET_LUN

      READF, unit, numfiles

      IF (numfiles NE 0) THEN BEGIN
         goodindex = 0
         nameindex = 0

         WHILE((nameindex LT (numfiles)) AND (goodindex EQ 0)) DO BEGIN

            READF, unit, name, dim, des, del

            IF (name EQ filename) THEN $
               goodindex = 1
            IF (del NE '*') THEN $
               MESSAGE, "* delimiter not found in MPI_data.txt"

            ENDWHILE

         FREE_LUN, unit

         IF (goodindex NE 0) THEN BEGIN
            ofilename = filename

            OPENR, unit, DEMO_FILEPATH(filename, SUBDIRECTORY=subdirectory), $
               /GET_LUN, /BLOCK

            IF KEYWORD_SET(assoc_it) THEN $
               newdata = ASSOC(unit, BYTARR(dim[0],dim[1])) $
            ELSE BEGIN
               newdata = BYTARR(dim[0], dim[1], dim[2])
               READU, unit, newdata
               FREE_LUN, unit
               ENDELSE

            description = des
            dimensions = dim

            ENDIF
         ENDIF
      ENDIF $
   ELSE $
      IF ((XRegistered("mpi_filterGetImage") EQ 0)) THEN BEGIN

         FILTER = 0

         IF KEYWORD_SET(one_dim) THEN $
            filter = filter OR one_mask

         IF KEYWORD_SET(two_dim) THEN $
            filter = filter OR two_mask OR three_mask

         IF KEYWORD_SET(three_dim) THEN $
            filter = filter OR three_mask

         IF (FILTER EQ 0) THEN $
            filter = one_mask + two_mask + three_mask

         OPENR, unit, DEMO_FILEPATH("MPI_data.txt", SUBDIRECTORY=subdirectory), $
            /GET_LUN

         READF, unit, numfiles

         IF (numfiles NE 0) THEN BEGIN

            names = STRARR(numfiles)
            descriptions = STRARR(numfiles)
            dimensions = LONARR(3, numfiles)
            goodindex = 0

            FOR nameindex = 0, numfiles - 1 DO BEGIN

               READF, unit, name, dim, des, del

               tempfilt = 0

               IF dim[0] GT 1 THEN $
                  tempfilt = one_mask

               IF dim[1] GT 1 THEN $
                  tempfilt = two_mask

                ;  Restrict animations to more than 2 images
               IF dim[2] GT 2 THEN $
                  tempfilt = three_mask

               IF ((tempfilt AND filter) NE 0) THEN BEGIN
                  names[goodindex] = name
                  dimensions[*, goodindex] = dim
                  descriptions[goodindex] = des
                  goodindex = goodindex + 1
                  ENDIF

               IF (del NE '*') THEN $
                  button = $
                      DIALOG_MESSAGE("* delimiter not found in MPI_data.txt", $
                      /ERROR)

               ENDFOR
            FREE_LUN, unit

            neworder = SORT(names[0:goodindex-1])
            names = names[neworder]
            descriptions = descriptions[neworder]
            dimensions = dimensions[*, neworder]

            IF (NOT(KEYWORD_SET(title))) THEN $
               title =  "Please Select an Image"

            
            loadbase = WIDGET_BASE(TITLE=title, $
                                   /COLUMN, TLB_FRAME_ATTR=1, $
                                   GROUP_LEADER=group, $
                                   MODAL = group ne 0L, $
                                   YOFFSET=40, XOFFSET=20)

            loadlist = WIDGET_LIST(loadbase, VALUE=descriptions, $
                         UVALUE="FILELST", YSIZE=10)

            buttonbase = WIDGET_BASE(loadbase, /ROW)

            loadok = WIDGET_BUTTON(buttonbase, VALUE="  OK  ", UVALUE="OK")

            loadcancel = WIDGET_BUTTON(buttonbase, VALUE="Cancel", $
                           UVALUE="CANCEL")

            WIDGET_CONTROL, loadbase, /REALIZE

            sState = {ColorTable:colorTable}

            WIDGET_CONTROL, loadbase, SET_UVALUE=sState, /NO_COPY
  

            selection = 0

            WIDGET_CONTROL, loadlist, SET_LIST_SELECT=selection

            Xmanager, "mpi_filterGetImage", loadbase, GROUP_LEADER=GROUP, $
              CLEANUP='mpi_filterGetImageCleanup', $
              EVENT_HANDLER="mpi_filterGetImageEvent"

            IF (selection NE 0) THEN BEGIN

               IF ((NOT(KEYWORD_SET(three_dim))) AND $
                     (dimensions[2, selection-1] NE 1)) THEN $
                  dimensions[2, selection-1] = 1

               OPENR, unit, DEMO_FILEPATH(names[selection - 1], $
                                          SUBDIRECTORY=subdirectory), $
                  /GET_LUN, /BLOCK

               ofilename = names[selection - 1]

               IF KEYWORD_SET(assoc_it) THEN $
                  newdata = ASSOC(unit, BYTARR(dimensions[0, selection-1], $
                              dimensions[1, selection-1])) $
               ELSE BEGIN
                  newdata = BYTARR(dimensions[0, selection - 1], $
                              dimensions[1, selection - 1], $
                              dimensions[2, selection - 1])

                  READU, unit, newdata

                  FREE_LUN, unit

                  ENDELSE

               description = descriptions[selection - 1]
               dimensions = dimensions[*,selection - 1]

               ENDIF
            ENDIF
         ENDIF

END
;------------- end of procedure mpi_filterGetImage --------------

; -----------------------------------------------------------------------------
;
;  PURPOSE:  This routine gets the name of the new data file to be
;            loaded using the procedure getfile.  It then loads the
;            data, resets the controls, and computes the new FFT of
;            the new data.  It also erases the old filtered image
;            and the old filtered power spectrum.
;   
;  KEYWORDS:
;      NAME  - if name is specified, the user is not
;            prompted with a choice and name is used in its
;            place.
;
pro mpi_filterLoadData, $
    NEW=newFile, $       ; IN: (opt) indicate a new file to load
    DATADESC=datadesc, $ ; OUT (opt) Data set description
    STATUS=status, $     ; OUT: 1 = successful, 0 = failure
    GROUP = group        ; IN: (opt) Group leader identifier

    COMMON filterCommon

    ;  Load a new signal or image from the "data" subdirectory.
    ;
    if N_ELEMENTS(group) eq 0 then group = 0L

    mpi_filterGetImage, newdata, DESCRIPTION=datadesc, $
      GROUP=group, $
      DIMENSIONS=datadim, /ONE_DIM, $
      /TWO_DIM, TITLE="Please Select Sample For Filtering", $
      SUBDIRECTORY=["examples", "data"]

    ;  Was the mpi_filterGetImage successful.
    ;
    fileOK = (SIZE(newdata))[0] GT 0
    status = fileOK

    ;  If mpi_filterGetImage was successful
    ;
    if (fileOK) then begin
        if (WIDGET_INFO(filterbase, /VALID_ID)) then begin
            TLBTitle = 'MPI filtering : ' + datadesc
            WIDGET_CONTROL, filterBase, TLB_SET_TITLE= TLBTitle
        endif

        if (KEYWORD_SET(newdata)) then begin
            if (datadim[1] NE 1) then begin ;2D?
                if (datadim[0] LT (VIEWSIZE / 2)) then begin
                    datadim[0] = datadim[0] * (VIEWSIZE / datadim[0])
                    datadim[1] = datadim[1] * (VIEWSIZE / datadim[1])
                    image = REBIN(newdata, datadim[0], datadim[1])
                endif else $
                    if (datadim[0] GT VIEWSIZE) then begin
                        datadim = datadim < [ viewsize, viewsize]
                        image = CONGRID(newdata, datadim[0], datadim[1])
                endif else begin
                    image = newdata
                endelse
                image = BYTSCL(image, TOP=!D.TABLE_SIZE-1)
            endif else  begin   ;1D case
                image = newdata
            endelse
        endif
   
        if (NOT(KEYWORD_SET(newFile))) then $
            mpi_filterSetSliders
      
        ;  Compute the forward fast MPI transform.
        ;
        imagetrans = (FFT(image, -1))      
        filteredImage = 0
        filteredSpectrum = 0
    endif

end        ;   of   mpi_filterLoadData


pro mpi_filterColorMapCallback

COMMON filterCommon

; For visuals with static colormaps, update the graphics
; after a change by XLOADCT.  This routine will be specified
; as the UPDATECALLBACK routine and will be called by XLOADCT.
; See the call of XLOADCT, below.
if ((COLORMAP_APPLICABLE(redrawRequired) GT 0) AND $
    (redrawRequired GT 0)) then begin
   for i = 0, NVIEWS - 1 do mpi_filterResetView, i 
endif

end

; -----------------------------------------------------------------------------
;
;  PURPOSE:  This procedure computes the correct filter depending
;            on the bandpass type.  Here a high pass filter is
;            one minus the low pass filter.
;
pro mpi_filterResetFilter

COMMON filterCommon

;;;; help, order, cutoff, bandwidth, filtertype, bandtype

WIDGET_CONTROL, filterbase, /HOURGLASS ;  Display busy while computing filter.

    ;  Reset the filtered image and the filtered spectrum since
    ;  changing the filter affects them.
    ;
filteredImage = 0
filteredSpectrum = 0

imagesize = SIZE(image)         ;Obtain dimensions
nx = imagesize[1]
ny = imagesize[2]
is2D = imagesize[0] gt 1        ;TRUE if an image
   
if (N_ELEMENTS(distfun) NE N_ELEMENTS(image)) then begin ;New distfun?
    if is2D then begin          ; If 2 dimensional...
        distfun = double( DIST(nx))
    endif else begin            ; Else 2D case...
        distfun = DINDGEN(nx)   ;  1D Euclidean dist function.
        distfun = distfun < (nx - distfun)
    endelse
    distfun[0] = 1d-4           ; Avoid division by 0
endif                           ; New dist fcn

case filtertype of
        ;  Define Butterworth filter types.
        ;
    0: BEGIN  
        if (bandtype EQ 0) then begin ; Compute lowpass filter.
          tao = order*3.
          r = (1./tao) * exp(-distfun[*, 0]/tao)
          filter = distfun * 0.
          if is2D EQ 0 then begin
            ;filter = 1.0d / (1.0d + (distfun / cutoff)^(2 * order)) 
            filter[nx/2:*] = r[0:nx/2-1]
            filter=REVERSE(filter, /OVERWRITE)
          endif else begin
            for mi = nx/2-10, nx/2+10 do filter[nx/2:*, mi] = r[0:nx/2-1]
            filter=REVERSE(filter, 1, /OVERWRITE)
          endelse
        endif else if (bandtype EQ 1) then begin ; Compute highpass filter.
            filter = 1.0d / (1.0d + (cutoff / distfun)^(2 * order))
        endif else begin        ;Bandpass or bandreject
            filter = distfun^2 - cutoff^2 ;Dist squared
            zeroes = WHERE(abs(filter) EQ 0.0, count)
            if (count NE 0) then filter[zeroes] = 1d-4 ;Avoid divide by 0
            filter = 1.0d / (1.0d + $ ; Compute bandreject filter.
                             ((distfun * bandwidth) / filter) ^ (2 * order))
            
            if (bandtype EQ 2) then $ ; Compute bandpass filter.
              filter = 1.0 - filter
        endelse
end                             ;     of  0
      

1: BEGIN                        ;  Define Exponential filter types.
    if (bandtype EQ 0) then begin ; Compute lowpass filter.
        filter = EXP(-(distfun / cutoff)^order)
    endif else if (bandtype EQ 1) then begin ;Highpass
        filter = 1.0d - EXP(-(distfun / cutoff)^order)
    endif else begin ;  Bandpass / reject,  avoid underflow in EXP function
        filter = (distfun^2 - cutoff^2) / (distfun * bandwidth)
        filter = EXP(-(filter ^ (2 * order) < 25))               
        if (bandtype EQ 3) then $ ; Compute bandreject filter.
          filter = 1.0d0 - filter
    endelse
endcase                         ;   exponential
      
2: BEGIN                        ;  Define Ideal Mathematical filter types.
    case bandtype of
        0: filter = (distfun LE cutoff) ;Low pass
        1: filter = (distfun GE cutoff) ;High pass
        2: filter = abs(distfun - cutoff) le (bandwidth / 2.0) ;Bandpass
        3: filter = abs(distfun - cutoff) gt (bandwidth / 2.0) ;Bandreject
    endcase
endcase                         ;   of  2
   
else: MESSAGE, "mpi_filterResetFilter: Bad Filter Type"
      
endcase                         ;   of filtertype
   
for i = 0, NVIEWS - 1 do begin
    if (view[i].viewtype GT 1) then begin ; Redraw  windows that change
        mpi_filterResetView, i 
    endif
endfor
end           ;   of  mpi_filterResetFilter



; -----------------------------------------------------------------------------
;
;  PURPOSE:  This procedure resets the view number passed in.  It 
;            is assumed that the structure associated with the 
;            passed in view has been altered in view type so it
;            needs to be redrawn.  In the case of one dimensional
;            data, there is no filter intensity image and cross -
;            sections of the filter are the same as filter plots
;            so they just call those routines.  Once the view has
;            been updated, the button that controls that view is
;            renamed to reflect the new view type.
;
pro mpi_filterResetView, $
    windowIndex  ; IN: window indes identifier

    COMMON filterCommon

    ;  Set the active draw window.
    ;
    WSET, view[windowIndex].windownum

    ;  If the view is not "saved" then erase the contents of the draw window.
    ;
    if (view[windowIndex].viewtype NE 8) then begin 
        ERASE, COLOR=100
    endif
   
    ;  Set the character scaling factor.
    ;
    DEVICE, GET_SCREEN_SIZE=screenSize
    if (screenSize[0] LT 850) then charscale = 5.5/!D.X_CH_SIZE $
    else charscale = 7.5/!D.X_CH_SIZE

    imagesize = SIZE(image)
    nx = imagesize[1]
    ny = imagesize[2]
    is2D = imagesize[0] gt 1    ;TRUE if an image
  
    ;  Set the frequency in the x axis.
    ;
    x = FINDGEN(nx) - nx/2

    ;  Make sure to redraw the original image
    ;  with the correct color table if necessary.
    ;
    if (!D.N_COLORS GT 256) then begin

        TVLCT, red, green, blue, /GET
        if ( TOTAL(ABS(red - currentRed)) GT 0 ) then begin
            TVLCT, red, green, blue
            currentRed = red
            WSET, view[0].windownum
            if (imagesize[0] GT 1) then begin
                TVSCL, image, (VIEWSIZE - nx) / 2, (VIEWSIZE - ny) / 2 
            endif else begin
                PLOT, image, /XSTYLE, /YSTYLE, $
                CHARSIZE=charscale
            endelse
        endif
        WSET, view[windowIndex].windownum
    endif

   

    case view[windowIndex].viewtype of 
         
        ;  Redraw original signal or image data.
        ;
        0 : begin 
            if is2D then begin
                TVSCL, image, (VIEWSIZE - nx) / 2, (VIEWSIZE - ny) / 2 
            endif else begin
                PLOT, image, /XSTYLE, /YSTYLE, $
                XTITLE='Time', $
                CHARSIZE=charscale
            endelse
        end    ;   of   0
   
        ;  Redraw logarithmic power spectrum.
        ;
        1 : begin 
            if (is2D) then begin
                TVSCL, SHIFT(ALOG(ABS(imagetrans) > 1.0e-10), nx / 2, ny / 2),$
                    (VIEWSIZE - nx) / 2, (VIEWSIZE - ny) / 2 
             endif else begin
                 PLOT, x, SHIFT(ABS(imagetrans), nx/2), $
                   XTITLE='Frequency', YTITLE='Power', $
                   /YLOG, $
                   /XSTYLE, charsize = charscale
             endelse
        end    ;   of   1

        ;  Redraw filter.
        ;
        2 : begin 
            if (is2D) then begin
                grid_size = 32
                thinfilter = CONGRID(filter, grid_size, grid_size)
                SURFACE, thinfilter,$
                    XSTYLE=4, YSTYLE=4, ZSTYLE=4
                 ; SURFACE, SHIFT(thinfilter, grid_size/2, grid_size/2),$
                  ;  XSTYLE=4, YSTYLE=4, ZSTYLE=4
            endif else begin
;                PLOT, x, SHIFT(FLOAT(FFT(filter)), nx / 2), $
;                    XTITLE='Frequency', $
;                    XSTYLE=2, YSTYLE=2, $
;                    CHARSIZE=charscale
                  PLOT, x, filter, $
                    XTITLE='Time', $
                    XSTYLE=2, YSTYLE=2, $
                    CHARSIZE=charscale
            endelse
            EMPTY
        end    ;   of   2
   
        ;  Redraw filter intensity.
        ;
        3 : begin 
            if (is2D) then begin
                 TVSCL, filter, $
                   (VIEWSIZE - nx) / 2, (VIEWSIZE - ny) / 2 
            endif 
        end    ;   of   3
   
        ;  Redraw filter cross section.
        ;
        4 : begin 
            if (is2D) then begin
                PLOT, x, filter[*, ny/2], /XSTYLE, $
                    XTITLE='Time', $
                    YSTYLE=2, $
                    CHARSIZE=charscale
                EMPTY
            endif
        end    ;   of   4
            
        ;  Redraw filter shaded surface.
        ;
        5 : begin 
            grid_size = 48
            if (is2D) then begin
                thinfilter = CONGRID(filter, grid_size, grid_size)
                SHADE_SURF, thinfilter, $
                     XSTYLE=4, YSTYLE=4, ZSTYLE=4
            endif

        end    ;   of   5
          
        ;  Redraw filtered logarithmic power spectrum.
        ;
        6 : begin 
            if (N_ELEMENTS(filteredSpectrum) LE 1) then begin
                filteredSpectrum = imagetrans * FLOAT(FFT(filter))
            endif

            if (is2D) then begin
                TVSCL, SHIFT(ALOG(ABS(filteredSpectrum)> 1e-10), $
                             nx / 2, ny / 2), $
                  (VIEWSIZE - nx) / 2, (VIEWSIZE - ny) / 2 
            endif else begin
                PLOT, x, SHIFT(ABS(FLOAT(filteredSpectrum)) > 1.0E-12, $
                               nx/2), XSTYLE=1, YSTYLE=1, $
                  XTITLE='Frequency', YTITLE='Power', $
                  /YLOG, $
                  CHARSIZE=charscale
            endelse
        end    ;   of   6

        ; Redraw filtered signal or image.
        ;
        7 : begin 

            if (N_ELEMENTS(filteredImage) LE 1) then begin
                filteredImage = FLOAT(FFT(imagetrans * FLOAT(FFT(filter)), 1))
            endif

            if (is2D) then begin
             ;   filteredSpectrum = imagetrans * FLOAT(FFT(filter))
             ;   filteredImage = FFT(filteredSpectrum, /INVERSE)
                filteredImage = CONVOL(DOUBLE(Image), filter[nx/2-20:nx/2+20, ny/2-20:ny/2+20], /NORMALIZE)
                TVSCL, filteredImage, (VIEWSIZE - nx) / 2, (VIEWSIZE - ny) / 2 
            endif else begin
                filteredImage_temp = CONVOL([DBLARR(600)*0+MEAN(Image), Image, DBLARR(600)*0+MEAN(Image)], filter, /NORMALIZE)
                filteredImage = filteredImage_temp[600:(600+nx-1)]
                PLOT, filteredImage, $
                    XTITLE='Time', $
                    CHARSIZE=charscale
            endelse
        end    ;   of   7
            
        ;  Saved view. (Do nothing)
        ;
        8 : 
   
        ELSE:
      
    endcase

    ;  If needed, empty the view.
    ;
    EMPTY

end          ;   of mpi_filterResetView


; -----------------------------------------------------------------------------
;
;  PURPOSE:  This is the main event handler for the Filter demo.
;
pro mpi_filterEvent, $
    event              ; IN: event structure

    COMMON filterCommon

    if (TAG_NAMES(event, /STRUCTURE_NAME) EQ $
        'WIDGET_KILL_REQUEST') then begin
        WIDGET_CONTROL, event.top, /DESTROY
        RETURN
    endif

    WIDGET_CONTROL, event.id, GET_UVALUE=retval

    if (N_ELEMENTS(retval) EQ 0) then retval = event.value

    ;  Take the following action based on the corresponding event.
    ;
    case retval of

        ;  Load new signal or image data.
        ;
        "Open..." : begin
            mpi_filterLoadData, GROUP = event.top
            mpi_filterResetFilter
            for i = 0, NVIEWS - 1 do begin
                if (view[i].viewtype LT 2) then begin
                    mpi_filterResetView, i
                endif
            endfor

            ;  Viewindex is always 1 now ( 2nd window..).
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 7) then begin
                view[viewindex].viewtype = 7
                WIDGET_CONTROL, wImageButton, SET_VALUE='Filtered signal'
                mpi_filterResetView, viewindex
            endif

            ; Desensitize  plots options if the image is 1D.
            ; 
            imagesize = size(image)

            WIDGET_CONTROL, wIntensityButton, SENSITIVE=imagesize[0] NE 1
            WIDGET_CONTROL, wCrossButton, SENSITIVE=imagesize[0] NE 1
            WIDGET_CONTROL, wShadedButton, SENSITIVE=imagesize[0] NE 1

        end     ;   of  Open...

        ; Load pre-defined color table tool.
        ;
        "Color Palette" : XLOADCT, /SILENT, GROUP=event.top, $
                                   UPDATECALLBACK='mpi_filterColorMapCallback'

        ; Reset the filter if the order value changed.
        ;
        "ORDERSELECTION" : begin
            if (event.index+1 NE order) then begin
                order = event.index+1
                mpi_filterResetFilter
            endif
        end     ;   of  ORDERSELECTION

        ; Reset the filter if frequency cutoff has changed.
        ;
        "FILTERCUTOFF" : begin  
            cutoff = event.value
            mpi_filterResetFilter
        end     ;   of FILTERCUTOFF

        ; Reset the filter if frequency bandwidth has changed and also if 
        ; bandwidth matters with current filter type.
        ;
        "FILTERBANDWIDTH" : begin
            bandwidth = event.value
            mpi_filterResetFilter
        end     ;   of FILTERBANDWIDTH

        ;  Select the filter selection.
        ;
        "FILTERSELECTION" : begin
            if (filtertype NE event.index) then begin
                filtertype = event.index
                mpi_filterResetFilter
            endif
        end     ;   of  FILTERSELECTION

        ;  Select the  band.
        ;
        "BANDSELECTION" : begin
            if (bandtype NE event.index) THEN BEGIN
                bandtype = event.index
                mpi_filterResetFilter

                ;  If not bandpass and not bandreject,
                ;  sensitize or desensitize the bandwidth slider.
                ;
                WIDGET_CONTROL, wBandWidthSlider, $
                  SENSITIVE= (bandtype eq 2) or (bandtype eq 3)
            endif
        end     ;   of  BANDSELECTION

        ;  Quit this application.
        ;
        "Quit": begin
	    ;  Deallocate variables in common.
            ;  Destroy widget hierarchy.
            WIDGET_CONTROL, event.top, /DESTROY
        end     ;   of  QUIT

        ;  View the original signal.
        ;
        "Original" : begin

            ;  Viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 0) then begin
                view[viewindex].viewtype = 0
                WIDGET_CONTROL, wImageButton, SET_VALUE='Original Signal'
                mpi_filterResetView, viewindex
            endif
        end     ;   of Original

        ;  View the log power spectrum
        ;
        "Logpower" : begin
            ;  Viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 1) then begin
                view[viewindex].viewtype = 1
                WIDGET_CONTROL, wImageButton, SET_VALUE='Log Power Spectrum'
                mpi_filterResetView, viewindex
            endif
        end     ;   of Logpower

        ;  View the filter plot
        ;
        "Filterplot" : begin
            ;  Viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 2) then begin
                view[viewindex].viewtype = 2
                WIDGET_CONTROL, wImageButton, SET_VALUE='Filter Plot'
                mpi_filterResetView, viewindex
            endif
        end     ;   of Filterplot

        ;  View the filter intensity as an image.
        ;
        "Intensity" : begin
            ; viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if( view[viewindex].viewtype NE 3) then begin
                view[viewindex].viewtype = 3
                WIDGET_CONTROL, wImageButton, $
                    SET_VALUE='Filter Intensity Image'
                mpi_filterResetView, viewindex
            endif
        end     ;   of  Intensity

        ;  View the filter cross section plot.
        ;
        "Cross" : begin
            ; viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if( view[viewindex].viewtype NE 4) then begin
                view[viewindex].viewtype = 4
                WIDGET_CONTROL, wImageButton, SET_VALUE='Filter Cross Section'
                mpi_filterResetView, viewindex
            endif
        end     ;   of  Intensity

        ;  View the filter as a shaed plot.
        ;
        "Shaded" : begin
            ; viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 5) then begin
               view[viewindex].viewtype = 5
               WIDGET_CONTROL, wImageButton, $
                   SET_VALUE='Shaded Surface of Filter'
               mpi_filterResetView, viewindex
            endif
        end     ;   of  Shaded

        ;  View the log filtered power spectrum
        ;
        "Logfilter" : begin
            ; viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 6) then begin
               view[viewindex].viewtype = 6
               WIDGET_CONTROL, wImageButton, $
                   SET_VALUE='Log Filtered Power Spectrum'
               mpi_filterResetView, viewindex
            endif
        end     ;   of  Logfilter

        ;  View the filtered (processed) signal.
        ;
        "Filtersignal" : begin
            ; viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 7) then begin
                view[viewindex].viewtype = 7
                WIDGET_CONTROL, wImageButton, SET_VALUE='Filtered Signal'
                mpi_filterResetView, viewindex
            endif
        end     ;   of Filtersignal

        ;  View the saved image (or plot).
        ;
        "Save" : begin
            ; viewindex is always 1 now ( 2nd window..)
            ;
            viewindex = 1
            if (view[viewindex].viewtype NE 8) then begin
               view[viewindex].viewtype = 8
               WIDGET_CONTROL, wImageButton, SET_VALUE='Save View'
               mpi_filterResetView, viewindex
            endif
        end     ;   of  Save

        ELSE:

    endcase  ;  of  retval

end       ;   of  mpi_filterEvent

; -----------------------------------------------------------------------------
;
;  PURPOSE:  Cleanup procedure.
;
pro mpi_filterCleanup, wFilterWindow
    COMMON filterCommon

    ;  Get the color table saved in the window's user value.
    ;
    WIDGET_CONTROL, wFilterWindow, GET_UVALUE=previousState
   
    ;  Restore the previous color table.
    ;
    TVLCT, previousState.colorTable

    ;  Restore previous margins.
    ;
    !X.MARGIN = previousState.xMargins
    !Y.MARGIN = previousState.yMargins

    ;  Clean up images in common.
    ;
    image = 0 & imagetrans=0 & filter = 0 & filteredImage = 0
    filteredSpectrum = 0
    distfun = 0
   
    if widget_info(previousState.groupBase, /valid) then $
        widget_control, previousState.groupBase, /map

end    ;   of mpi_filterCleanup

; -----------------------------------------------------------------------------
;
;  PURPOSE:  Main procedure of filter
;
pro mpi_filter, $
    Group= group, $      ; IN: (opt) group leader Identifer
    FILENAME=filename, $ ; IN: (opt) Data filename to view at atart up.
    RECORD_TO_FILENAME=record_to_filename, $
    APPTLB= appTLB       ; OUT: (opt) top level base of this application

    COMMON filterCommon

    if n_elements(group) eq 0 then group = 0L

    ;  Only one instance at a time.
    ;
    if (XRegistered("mpi_filter")) then begin
        GOTO, FINISH
    endif
   
    ;  Get the current color vectors to restore when this application is exited.
    ;
    TVLCT, savedR, savedG, savedB, /GET

    ;  Build color table from color vectors.
    ;
    colorTable = [[savedR],[savedG],[savedB]]
   
    ;  Save the current margins in order to restore it when exiting.
    ;
    xMargins = !X.MARGIN
    yMargins = !Y.MARGIN
      
    ;  Save items to be restored on exit in a structure.
    ;
    previousState = { $
        colorTable:colorTable, $
        xMargins:xMargins, $
        yMargins:yMargins, $
        groupBase: group $
    }
   
    !X.MARGIN = [10, 4]
    !Y.MARGIN = [4, 2]

    ;  Initialize filter states.   
    ;
    filteredImage = 0  
    filteredSpectrum = 0
   
    ;  Initialize the filter parameters.
    ;  
    order = 2        ; order   
    filtertype = 0   ; Butterworth filter 
    cutoff = 20.     ; Cutoff frequency
    bandwidth = 10.  ; Filter bandwidth
    bandtype = 0     ; Low pass filter (0=low pass , 1=high pass, 
                     ;                  2=band pass, 2=band reject).
   
    ;  Determine the hardware device size in pixels.
    ;
    DEVICE, GET_SCREEN=screendims
   
    if screendims[0] ge 1100 then i = 2 $	;General screen sizing
    else if screendims[0] ge 800 then i = 1 $
    else i = 0
   
    ;  Number of views and view size in pixels (device).
    ;  Try to have as many factors of 2 as possible for FFT efficiency.
    viewsize = ([ 224, 256, 320])[i]
    nviews = ([ 2, 2, 2 ])[i]

    fileSet = N_ELEMENTS(filename)

    ;  Was the mpi_filterGetImage successful.
    ;
    fileOK = (SIZE(newdata))[0] GT 0
    status = fileOK

    filterbase = -1L

    if (fileSet NE 0) then begin

        mpi_filterGetImage, newdata, DESCRIPTION=datadesc, $
            FILENAME=filename, $
            DIMENSIONS=datadim, /ONE_DIM, $
            /TWO_DIM, TITLE="Please Select Sample For Filtering", $
            SUBDIRECTORY = ["examples", "data"]

        ;  Was the Getdata2 successful.
        ;
        fileOK = (SIZE(newdata))[0] GT 0
        status = fileOK

        ;  If mpi_filterGetImage was successful
        ;
        if (fileOK) then begin
              if (KEYWORD_SET(newdata)) then begin
                if (datadim[1] NE 1) then begin
                    if (datadim[0] LT (VIEWSIZE / 2)) then begin
                        datadim[0] = datadim[0] * (VIEWSIZE / datadim[0])
                        datadim[1] = datadim[1] * (VIEWSIZE / datadim[1])
                        image = REBIN(newdata, datadim[0], datadim[1])
                    endif else $
                        if (datadim[0] GT VIEWSIZE) then begin
                            datadim = datadim < [ viewsize, viewsize]
                            image = CONGRID(newdata, datadim[0], datadim[1])
                    endif else begin
                        image = newdata
                    endelse
                    image = BYTSCL(image, TOP=!D.TABLE_SIZE-1)
                endif else  begin
                    image = newdata
                endelse
              endif


             ;  Compute the forward fast MPI transform.
             ;
            imagetrans = (FFT(image, -1))
      
            filteredImage = 0
            filteredSpectrum = 0
        endif

        status = 1

    endif else begin

        ;  Load the initial data.
        ;
        mpi_filterLoadData, NEW=1, DATADESC=datadesc, STATUS=status 

    endelse

    if (status EQ 0) then begin

        ;  Restore the previous color table.
        ;
        TVLCT, previousState.colorTable

        ;  Restore previous margins.
        ;
        !X.MARGIN = previousState.xMargins
        !Y.MARGIN = previousState.yMargins

        RETURN
    endif

    ;  Set the Top level base title string.
    ;
    TLBTitle = 'MPI Filtering : ' + datadesc

    ;  Get the screen size.
    ;
    DEVICE, GET_SCREEN_SIZE = screenSize

    ;  Create the starting up message.
    ;
    drawbase = demo_startmes(GROUP=group)

    ;  Load a particular default color table.
    ;
    LOADCT, 0, /SILENT

    TVLCT, red, green, blue, /GET
    currentRed = red

    ;  If no data was selected, drop out of routine.
    ;
    if (NOT(KEYWORD_SET(image))) then $
        GOTO, FINISH

    ;  Define an array of view structures.
    ;
    view = REPLICATE({vinst, viewtype:0b, $
        windowid:0L, windownum:0}, NVIEWS)
   
    ;  The default view settings that are used
    ;  for differing numbers of views.
    ;
    view.viewtype = ([0,7,2,1,6])[0:nviews-1]

    ;  The types of  views we can show.
    ;
    viewoptions = [ "Original Signal", "Log Power Spectrum", "Filter Plot", $
	"Filter Intensity Image", "Filter Cross Section", $
	"Shaded Surface of Filter", "Log Filtered Power Spectrum", $
	"Filtered Signal", "Save View" ]
   
    ;  Define a main filter base.
    ;
    filterBase = WIDGET_BASE(TITLE=TLBTitle, $
                             /COLUMN, $
                             MAP=0, $
                             /TLB_KILL_REQUEST_EVENTS, $
                             GROUP_LEADER=group, $
                             MBAR=mbar, UVALUE=previousState, TLB_FRAME_ATTR=1)

    menubar = CW_PDMENU(mbar, /RETURN_NAME, /MBAR, /HELP, $
        ['1\File', '0\Open...', '2\Quit', '1\Edit', '2\Color Palette'])

        ;  Define a subbase.
        ;
        subbase1 = WIDGET_BASE(filterBase, COLUMN = 2)
            ;  Define the base for the left column
            ;  It contains all the buttons and sliders...
            ;
            leftBase = WIDGET_BASE(subbase1, /COLUMN)

                ; Define the filter base.
                ; It has all the 3 filter options ( function, order, type)
                ;
                filterOptBase = WIDGET_BASE(leftBase,/FRAME,/COLUMN)
         
                    ; Define the filter options buttons.
                    ;
                    wFilterLbl = WIDGET_LABEL(filterOptBase, $
                        VALUE='Filter', /ALIGN_CENTER)

                    wFunctionLbl = WIDGET_LABEL(filterOptBase, $
                        VALUE='Function', /ALIGN_LEFT)

                    wFunctionButton = WIDGET_DROPLIST(filterOptBase, $
                        UVALUE='FILTERSELECTION',$
                        /ALIGN_LEFT, $
                        VALUE=['Relaxation']);, 'Exponential','Ideal'])

                    wOrderLbl = WIDGET_LABEL(filterOptBase, VALUE='Constant [us]', $
                       /ALIGN_LEFT)

                    wOrderButton = WIDGET_DROPLIST(filterOptBase, $
                        UVALUE='ORDERSELECTION',$
                        /ALIGN_LEFT, Value=['1', '2','3','4','5','6'])

                        ;  Set the order button's initial value.
                        ;
                        WIDGET_CONTROL, wOrderButton, $
                        SET_DROPLIST_SELECT=order-1

;                    wTypeLbl = WIDGET_LABEL(filterOptBase, VALUE='Type', $
;                        /ALIGN_LEFT)
;
;                    wTypeButton = WIDGET_DROPLIST(filterOptBase, $
;                        UVALUE='BANDSELECTION',$
;                        /ALIGN_LEFT,  $
;                        VALUE =["Low Pass", "High Pass", $
;                        "Band Pass", "Band Reject"])

                sliderBase = WIDGET_BASE(Leftbase, /COLUMN)
    
                    ;  Define the floating-point sliders
                    ;  for cutoff frequency and bandwidth. 
                    ;
;                    mpi_filterSetSliders, sliderBase


            ;  Define a sub-base belonging to filterBase
            ;  which contains the viewing areas.
            ;
            viewbase = WIDGET_BASE(subBase1, /ROW)
   
                ;  Define the view windows.
                ;
                for viewindex = 0, NVIEWS - 1 do begin

                    tempbase = WIDGET_BASE(viewbase, $
                        /FRAME, /COLUMN, /BASE_ALIGN_CENTER)
   
                        view[viewindex].windowid = WIDGET_DRAW(tempbase, $
                            XSIZE=VIEWSIZE, YSIZE=VIEWSIZE, RETAIN=2)

                        ;  Create the filter method drop list.
                        ;
                        if (viewIndex EQ 0) then begin
                            wOriginalLabel = WIDGET_LABEL(tempbase, $
                                VALUE='Original Image', $
                                /ALIGN_CENTER)
                        endif else begin

                            wImageButton = WIDGET_BUTTON(tempBase, $
                                VALUE='Log Filtered Power Spectrum', $
                                UVALUE= viewindex, Menu = 1)

                                wOriginalButton = WIDGET_BUTTON(wImageButton, $
                                    VALUE='Original Signal', $
                                    UVALUE='Original')

                                wLogpowerButton = WIDGET_BUTTON(wImageButton, $
                                    VALUE='Log Power Spectrum', $
                                    UVALUE='Logpower')

                                wFilterplotButton = $
                                    WIDGET_BUTTON(wImageButton, $
                                    VALUE='Filter Plot', $
                                    UVALUE='Filterplot')

                                wIntensityButton = WIDGET_BUTTON(wImageButton, $
                                    VALUE='Filter Intensity Image', $
                                    UVALUE='Intensity')

                                wCrossButton = WIDGET_BUTTON(wImageButton, $
                                    VALUE='Filter Cross Section', $
                                    UVALUE='Cross')

                                wShadedButton = WIDGET_BUTTON(wImageButton, $
                                    VALUE='Shaded Surface of Filter', $
                                    UVALUE='Shaded')

                                wLogfilterButton = WIDGET_BUTTON(wImageButton, $
                                    VALUE='Log filtered Power Spectrum', $
                                    UVALUE='Logfilter')

                                wFiltersignalButton = $
                                    WIDGET_BUTTON(wImageButton, $
                                    VALUE='Filtered Signal', $
                                    UVALUE='Filtersignal')

;                                wSaveButton = WIDGET_BUTTON(wImageButton, $
;                                    VALUE='Save View', $
;                                    UVALUE='Save')

                        endelse
                endfor

        ;  Create tips texts.
        ;
        wStatusBase = WIDGET_BASE(filterBase, MAP=0, /ROW)

    ;  Display HOURGLASS during setup.
    ;
    WIDGET_CONTROL, filterBase, /REALIZE, /HOURGLASS

    ;  Place this application top level base ID into appTLB.
    ;
    appTLB = filterBase

    ;  Get the tips
    ;
    sText = demo_getTips('./lib/source_code/mpifilter.tip', $
                         filterBase, $
                         wStatusBase)

    ;  Select Butterworth as default filter function.
    ;
    WIDGET_CONTROL, wImageButton, SET_VALUE='Filtered Signal'
    WIDGET_CONTROL, wFunctionButton, SET_DROPLIST_SELECT = 0
   
    ;  Select lowpass as default filter type.
    ;
;    WIDGET_CONTROL, wTypeButton, SET_DROPLIST_SELECT = 0
   
    ;  Get the drawing area numbers and set the view structure array so 
    ;  these values can be referenced when drawing.
    ;

    for i = 0, NVIEWS - 1 do begin
        WIDGET_CONTROL, view[i].windowid, GET_VALUE=temp
        view[i].windownum = temp
    endfor

    ;  Draw the default views that don't rely on the filter and then reset 
    ;  the filter which automatically draws the filter's dependent views.
    ;
    mpi_filterResetFilter
    mpi_filterResetView, 0
   
    ; Desensitize  plots options if the image is 1D.
    ; 
    imagesize = size(image)
    if (imagesize[0] EQ 1) then begin
        WIDGET_CONTROL, wIntensityButton, SENSITIVE=0
        WIDGET_CONTROL, wCrossButton, SENSITIVE=0
        WIDGET_CONTROL, wShadedButton, SENSITIVE=0
    endif

    ;  Destroy the starting up window.
    ;
    WIDGET_CONTROL, drawbase, /DESTROY

    ;  Map the top level base.
    ;
    WIDGET_CONTROL, filterBase, MAP=1

    XMANAGER, "mpi_filter", filterBase, $
        EVENT_HANDLER="mpi_filterEvent", $
        CLEANUP="mpi_filterCleanup", $
        /NO_BLOCK

    FINISH:

end     ;    of  mpi_filter
