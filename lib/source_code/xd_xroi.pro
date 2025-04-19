;
; $Id: //depot/Release/ENVI51_IDL83/idl/idldir/lib/utilities/xd_xroi.pro#1 $
;
; Copyright (c) 1999-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;   xd_xroi
;
; PURPOSE:
;   This procedure serves as utility for generating, displaying, and
;   manipulating regions of interest.
;
; CATEGORY:
;   Utilites. Object graphics.
;
; CALLING SEQUENCE:
;   xd_xroi, [ImageData] [,r] [,g] [,b]
;
; OPTIONAL INPUT/OUTPUT:
;   ImageData: An array representing an image to be displayed.  Can be:
;       [m,n]   -- 8-bit image.
;       [3,m,n] -- 24-bit image
;       [m,3,n] -- 24-bit image
;       [m,n,3] -- 24-bit image
;
;       If ImageData is not supplied, the user will be prompted for
;       a file via dialog_pickfile.  On output, ImageData will be
;       set to the current image data.  (The current image data
;       can be different than the input image data if the
;       user Imported an image via xd_xroi's "Import Image..."
;       menu button.)
;   R: An array of bytes representing red color table values.  On
;       input, these values are applied to the image if the image
;       is 8-bit.  To get the red color table values for the image
;       on output from xd_xroi, set this argument to a named variable.
;       (If the image is 24-bit, this argument will output a 256-
;       element byte array containing the values given at input,
;       or bindgen(256) if R was undefined on input.)
;   G: An array of bytes representing green color table values.  On
;       input these values are applied to the image if the image
;       is 8-bit.  To get the green color table values for the image
;       on output from xd_xroi, set this argument to a named variable.
;       (If the image is 24-bit, this argument will output a 256-
;       element byte array containing the values given at input,
;       or bindgen(256) if G was undefined on input.)
;   B: An array of bytes representing blue color table values.  On
;       input these values are applied to the image if the image
;       is 8-bit.  To get the blue color table values for the image
;       on output from xd_xroi, set this argument to a named variable.
;       (If the image is 24-bit, this argument will output a 256-
;       element byte array containing the values given at input,
;       or bindgen(256) if B was undefined on input.)
;
; INPUT KEYWORD PARAMETERS:
;   BLOCK: If set, block IDL command line.
;   MODAL: If set, block other IDL widgets from receiving events.
;   GROUP: Group leader widget.
;   FLOATING: Set this keyword -- along with the GROUP_LEADER
;       keyword -- to create a "floating" top-level base widget.
;       If the windowing system provides Z-order control, floating
;       base widgets appear above the base specified as their group
;       leader. If the windowing system does not provide Z-order
;       control, the FLOATING keyword has no effect.
;   RENDERER:  Set this keyword to an integer value to indicate which
;                  graphics renderer to use when drawing objects within
;                  the window.  Valid values are:
;                      0 = Platform native OpenGL
;                      1 = IDL's software implementation (the default)
;   REGIONS_IN: Array of IDLgrROI references.  This keyword also accepts
;       -1, or OBJ_NEW() (Null object) to indicate that there are no
;       regions_in.
;   TOOLS:  Set this keyword to a string (or vector of strings) from the
;       list below to indicate which types of ROI manipulation tools
;       should be supported when the ROI example is run. If more than
;       one string is specified, a series of bitmap buttons will appear
;       at the top of the ROI example widget (to the right of the
;       fixed set of bitmap buttons used for saving regions, displaying
;       region information, and copying to clipboard). If only one
;       string is specified, no additional bitmap buttons will appear,
;       and the manipulation mode is implied by the given string. The
;       default is to activate all types of manipulation tools.
;
;       Valid strings:
;
;       'Translate-Scale' Translation and scaling.  Mouse down selects
;       a region, mouse motion either scales or translates (depending
;       upon where the mouse down occurred).
;
;       'Rectangle' Rectangle ROI drawing.  Mouse down begins a
;       rectangular region, mouse motion repositions a corner of the
;       rectangle, mouse up finishes the region.
;
;       'Ellipse' Ellipse ROI drawing.  Mouse down begins an
;       elliptical region, mouse motion sizes the ellipse, mouse up
;       finishes the region.
;
;       'Freehand Draw' Freehand ROI drawing. Mouse down begins a
;       region, mouse motion adds vertices to the region (following
;       the path of the mouse), mouse up finishes the region.
;
;       'CircularMask Draw' Circular Mask ROI drawing. Mouse down begins a
;       region, mouse motion adds vertices to the region (following
;       the path of the mouse), mouse up finishes the region.
;
;       'Polygon Draw'  Polygon ROI drawing. Mouse down begins a
;       region, mouse motion temporarily locates the next vertex,
;       mouse up fixes the location of that vertex, double click finishes
;       the region.
;
;       'Selection' ROI selection. Mouse down/up selects the nearest
;       region. The nearest vertex in that region is identified with a
;       yellow crosshair symbol.
;
;   X_SCROLL_SIZE: Set this keyword to the width of the drawing area.
;       The default (and also the maximum) is the width of the image.
;
;   Y_SCROLL_SIZE: Set this keyword to the height of the drawing area.
;       The default (and also the maximum) is the height of the image.
;
; OUTPUT KEYWORD PARAMETERS:
;   REGIONS_OUT: Array of IDLgrROI references.  This keyword is assigned
;       the null object refernce if there are no regions_out.
;   REJECTED: Those regions_in that are not in regions_out.  This keyword
;       is assigned the null object reference if no regions_in are rejected
;       by the user.
;   STATISTICS: Set this argument to a named variable to receive an array
;       of anonymous structures, one for each ROI that is valid when this
;       routine returns.  The structures will contain: {
;           count: 0UL, $       ; Number of pixels in region.
;           minimum: 0.0, $     ; Pixel value.
;           maximum: 0.0, $     ; Pixel value.
;           mean: 0.0, $        ; Mean pixel value.
;           stddev: 0.0 $       ; Standard deviation of pixel values.
;           }
;       When this routine exits, if ImageData is 24-bit at that time,
;       or if there are no valid regions of interest at that time,
;       parameter STATISTICS will be undefined.
;   ROI_GEOMETRY: Set this argument to a named variable to receive an array
;       of anonymous structures, one for each ROI that is valid when this
;       routine returns.  The structures will contain: {
;           area: 0.0, $
;           centroid: FLTARR(3), $
;           perimeter: 0.0 $
;           }
;       If there are no valid regions of interest when this routine returns,
;       ROI_GEOMETRY will be undefined.
;
; INPUT/OUTPUT KEYWORD PARAMETERS:
;   ROI_SELECT_COLOR: A 3-element byte array indicating the color of
;       ROI outlines when they are selected.  [r,g,b]
;   ROI_COLOR: A 3-element byte array indicating the color of
;       ROI outlines when they are not selected.  [r,g,b]
;
; EXAMPLE:
;   To sweep through a series of images, maintaining a single list of
;   regions as we go:
;
;       images = randomu(seed, 5, 5, 2)
;       images = bytscl(congrid(images, 400, 400, 5))
;       for i=0,4 do begin
;           xd_xroi, $
;               images[*, *, i], $
;               r, g, b, $
;               regions_in=regions, $
;               regions_out=regions, $
;               roi_select_color=roi_select_color, $
;               roi_color=roi_color, $
;               rejected=rejected, $
;               /bloc
;           obj_destroy, rejected
;           endfor
;       obj_destroy, regions
;
; MODIFICATION HISTORY:
;   Written by: DD, Aug 1999
;   6/2000, PCS - Added features and changed from example to xd_xroi command.
;   7/2001, DLD - Added rectangle, ellipse, and translate-scale tools.
;                 Added support for region growing.
;                 Added support for histogram plots for RGB images.
;                 Improve support for RGB image interleaving.
;   DLD, Sept 2002: Cleanup LIVE tools environment if necessary.
;   CT, Sept 2002: Add tooltips. Change to exclusive buttons.
;                  Disable context menu for segmented polygon ROI.
;                  Add X_SCROLL_SIZE and Y_SCROLL_SIZE.
;                  Use APP_SCROLL so we can handle huge images.
;   SM, June 2003: For histograms change LIVE_PLOT to iPlot
;   CT, Sept 2003: Destroy iPlots when widget exits
;-
pro xd_xroiLoadct__Cleanup, wID
  COMPILE_OPT hidden

  WIDGET_CONTROL, wID, GET_UVALUE=pState
  if PTR_VALID((*pState).pLeaderState) then begin
    (*(*pState).pLeaderState).wLoadCT = -1
  end
  PTR_FREE, pState
end
;----------------------------------------------------------------------------
pro xd_xroiLoadCT_event, event
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, event.top, GET_UVALUE=pState

  case 1 of
    TAG_NAMES(event, /STRUCTURE_NAME) eq 'CW_PALETTE_EDITOR_PM': begin
      (*(*pState).pLeaderState).oImage->SetProperty, $
        PALETTE=(*(*pState).pLeaderState).oPalette
      (*(*pState).pLeaderState).oWindow->Draw, $
        (*(*pState).pLeaderState).oView
    endcase
    event.id eq (*pState).wOK: begin
      WIDGET_CONTROL, event.top, /DESTROY
    end
    event.id eq (*pState).wCancel: begin
      (*(*pState).pLeaderState).oPalette->SetProperty, $
        RED_VALUES=(*pState).red_values, $
        GREEN_VALUES=(*pState).green_values, $
        BLUE_VALUES=(*pState).blue_values
      (*(*pState).pLeaderState).oImage->SetProperty, $
        PALETTE=(*(*pState).pLeaderState).oPalette
      (*(*pState).pLeaderState).oWindow->Draw, $
        (*(*pState).pLeaderState).oView
      WIDGET_CONTROL, event.top, /DESTROY
    endcase
    else:
  endcase
end
;----------------------------------------------------------------------------
pro xd_xroiLoadCT, group_leader
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, group_leader, GET_UVALUE=pLeaderState

  if WIDGET_INFO((*pLeaderState).wLoadCT, /VALID_ID) then begin
    WIDGET_CONTROL, (*pLeaderState).wLoadCT, /SHOW
    RETURN
  end

  prefix = 'xd_xroiLoadCT:'
  tlb = WIDGET_BASE( $  ; Top-Level Base.
    /COLUMN, $
    TITLE='Image Color Table', $
    GROUP_LEADER=group_leader, $
    UNAME=prefix + 'tlb', $
    MODAL=(*pLeaderState).modal $
    )

  (*pLeaderState).wPaletteEdit = CW_PALETTE_EDITOR( $
    tlb, $
    DATA=(*pLeaderState).oPalette, $
    UNAME=prefix + 'pal_edit', $
    /FRAME $
    )
  wRowBase = WIDGET_BASE(tlb, /ROW, /GRID)
  wOK = WIDGET_BUTTON(wRowBase, VALUE='OK', UNAME=prefix + 'OK')
  wCancel = WIDGET_BUTTON( $
    wRowBase, $
    VALUE='Cancel', $
    UNAME=prefix + 'cancel' $
    )

  (*pLeaderState).oPalette->GetProperty, $
    RED_VALUES=red_values, $
    GREEN_VALUES=green_values, $
    BLUE_VALUES=blue_values

  if not (*pLeaderState).modal then $
    WIDGET_CONTROL, tlb, MAP=0

  WIDGET_CONTROL, tlb, /REALIZE

  if not (*pLeaderState).modal then begin
    tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
    leader_geom = WIDGET_INFO(group_leader, /GEOMETRY)
    DEVICE, GET_SCREEN_SIZE=screen_size

    x = leader_geom.scr_xsize[0] + leader_geom.xoffset
    x = x < (screen_size[0] - tlb_geom.scr_xsize)
    x = x > 0

    WIDGET_CONTROL, tlb, TLB_SET_XOFFSET=x
    WIDGET_CONTROL, tlb, MAP=1
  endif

  WIDGET_CONTROL, tlb, SET_UVALUE=PTR_NEW({ $
    red_values: red_values, $
    green_values: green_values, $
    blue_values: blue_values, $
    wOK: wOK, $
    wCancel:wCancel, $
    pLeaderState: pLeaderState $
  })

  (*pLeaderState).wLoadCT = tlb
  XMANAGER, 'xd_xroiLoadCT', tlb, CLEANUP='xd_xroiLoadct__Cleanup', /NO_BLOCK
end


;------------------------------------------------------------------------------
pro xd_xroi__BaseResize, sEvent

  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  ; Compute the change in width and height of the base widget.
  deltaX = sEvent.x - (*pState).scr_xsize
  deltaY = sEvent.y - (*pState).scr_ysize
  if (deltaX eq 0) and (deltaY eq 0) then $
    return

  ; Retrieve the current draw size.
  draw_geom = WIDGET_INFO((*pState).wDraw, /GEOMETRY)

  ; On Motif, the geom.xsize may include room for the scrollbar even if
  ; no scrollbar is present. So restrict size to <= to the virtual canvas.
  draw_xsize = draw_geom.xsize < draw_geom.draw_xsize
  draw_ysize = draw_geom.ysize < draw_geom.draw_ysize

  ; Size must be > toolbar size and less than virtual canvas.
  newXsize = (*pState).toolbar_xsize > (draw_xsize + deltaX) $
    < draw_geom.draw_xsize

  ; Size must be > some minimum and less than virtual canvas.
  newYsize = 1 > (draw_ysize + deltaY) < draw_geom.draw_ysize


  ; Because of a Motif bug, if the base widget is made wider, but the
  ; draw widget is already at its maximum, then the resize is ignored,
  ; and base doesn't shrink back to fit the draw widget. To avoid this,
  ; we do two resizes, one slightly larger than the other.
  ; This forces the base to be resized.
  if (!version.os_family eq 'unix') then $
    WIDGET_CONTROL, (*pState).wDraw, $
    XSIZE=newXsize+1, YSIZE=newYsize

  ; Set viewport size.
  WIDGET_CONTROL, (*pState).wDraw, $
    XSIZE=newXsize, YSIZE=newYsize

  ; Get viewport offset (scrollbar positions).
  WIDGET_CONTROL, (*pState).wDraw, GET_DRAW_VIEW=viewport

  xoffset = viewport[0]
  yoffset = viewport[1]

  ; Decrease the viewplane offset if we can now display more of the
  ; image (or the entire image).
  xend = viewport[0] + newXsize
  if (xend gt draw_geom.draw_xsize) then $
    xoffset = (xoffset - (xend - draw_geom.draw_xsize)) > 0
  yend = viewport[1] + newYsize
  if (yend gt draw_geom.draw_ysize) then $
    yoffset = (yoffset - (yend - draw_geom.draw_ysize)) > 0

  (*pState).oView->SetProperty, $
    VIEWPLANE_RECT=[xoffset, yoffset, newXsize, newYsize]

  ; Retrieve the new base size and cache it.
  base_geom = WIDGET_INFO(sEvent.top, /GEOMETRY)

  ; On Motif, the widget resize event.y includes the menubar height,
  ; which corresponds to geom.scr_ysize.
  ; On Windows, the widget resize event.y does not include the menubar,
  ; so we need to cache geom.ysize instead.
  if (!version.os_family eq 'unix') then begin
    (*pState).scr_xsize = base_geom.scr_xsize
    (*pState).scr_ysize = base_geom.scr_ysize
  endif else begin
    (*pState).scr_xsize = base_geom.xsize
    (*pState).scr_ysize = base_geom.ysize
  endelse

end


;------------------------------------------------------------------------------
pro xd_xroi_event, sEvent

  COMPILE_OPT idl2, hidden

  if (TAG_NAMES(sEvent, /STRUC) eq 'WIDGET_BASE') then begin
    xd_xroi__BaseResize, sEvent
    ; We're done with this event.
    return
  endif

  WIDGET_CONTROL, sEvent.id, GET_UVALUE=uval
  print, uval

  case uval of

    'DRAW': begin
      ; Handle all events in the draw area.

      case sEvent.type of
        ; Button Press
        0: xd_xroi__ButtonPress, sEvent

        ; Button Release
        1: xd_xroi__ButtonRelease, sEvent

        ; Motion
        2: xd_xroi__Motion, sEvent

        ; Viewport moved (scrollbars)
        3: xd_xroi__Viewport, sEvent

        ; Expose
        4: xd_xroi__Expose, sEvent

        else: begin
        end
      endcase
    end

    'ADDLAYER': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).Modify_ROI = 'ADD'
      ; Get current list of ROI names.
      oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
      if (nROIs gt 0) then begin
        for i=0, nROIs-1 do begin
          (*pState).oROIModel->Remove, oROIs[i]
          (*pState).oROIGroup->Remove, oROIs[i]
          if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
          OBJ_DESTROY, oROIs[i]
        endfor

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        (*pState).oWindow->Draw, (*pState).oView
        ; Pick a new selected region.
        ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
        ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
        xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
      endif
    end

    'REMOVELAYER': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).Modify_ROI = 'REMOVE'
      ; Get current list of ROI names.
      oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
      if (nROIs gt 0) then begin
        for i=0, nROIs-1 do begin
          (*pState).oROIModel->Remove, oROIs[i]
          (*pState).oROIGroup->Remove, oROIs[i]
          if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
          OBJ_DESTROY, oROIs[i]
        endfor

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        (*pState).oWindow->Draw, (*pState).oView
        ; Pick a new selected region.
        ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
        ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
        xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
      endif
    end

    'BORDERSELECT': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).Modify_ROI = 'BORDER'
      ; Get current list of ROI names.
      oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
      if (nROIs gt 0) then begin
        for i=0, nROIs-1 do begin
          (*pState).oROIModel->Remove, oROIs[i]
          (*pState).oROIGroup->Remove, oROIs[i]
          if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
          OBJ_DESTROY, oROIs[i]
        endfor

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        (*pState).oWindow->Draw, (*pState).oView
        ; Pick a new selected region.
        ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
        ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
        xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
      endif
    end

    'ADDRED': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).Modify_ROI = 'ADDNEW1'
      ; Get current list of ROI names.
      oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
      if (nROIs gt 0) then begin
        for i=0, nROIs-1 do begin
          (*pState).oROIModel->Remove, oROIs[i]
          (*pState).oROIGroup->Remove, oROIs[i]
          if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
          OBJ_DESTROY, oROIs[i]
        endfor

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        (*pState).oWindow->Draw, (*pState).oView
        ; Pick a new selected region.
        ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
        ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
        xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
      endif
    end

    'ADDGREEN': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).Modify_ROI = 'ADDNEW2'
      ; Get current list of ROI names.
      oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
      if (nROIs gt 0) then begin
        for i=0, nROIs-1 do begin
          (*pState).oROIModel->Remove, oROIs[i]
          (*pState).oROIGroup->Remove, oROIs[i]
          if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
          OBJ_DESTROY, oROIs[i]
        endfor

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        (*pState).oWindow->Draw, (*pState).oView
        ; Pick a new selected region.
        ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
        ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
        xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
      endif
    end

    'CONTINUOUS': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).CONTORINDE='CONTINUOUS'
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
    end

    'DISCRETE': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).CONTORINDE='DISCRETE'
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
    end

    'THRESHOLD_min': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      WIDGET_CONTROL, (*pState).wThresholdSlider_min, GET_VALUE=threshold_min
      WIDGET_CONTROL, (*pState).wThresholdSlider_max, GET_VALUE=threshold_max
      (*pState).oImage->SetProperty, DATA=BYTSCL(*((*pState).pImg), MAX=threshold_max, MIN=threshold_min)
      (*pState).oWindow->Draw, (*pState).oView
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
    end

    'THRESHOLD_max': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      WIDGET_CONTROL, (*pState).wThresholdSlider_min, GET_VALUE=threshold_min
      WIDGET_CONTROL, (*pState).wThresholdSlider_max, GET_VALUE=threshold_max
      (*pState).oImage->SetProperty, DATA=BYTSCL(*((*pState).pImg), MAX=threshold_max, MIN=threshold_min)
      (*pState).oWindow->Draw, (*pState).oView
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
    end

    'THRESHOLD_radii': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      WIDGET_CONTROL, (*pState).wThresholdSlider_radii, GET_VALUE=Radii

      (*pState).Radii = Radii
      oROI = (*pState).oCurrROI
      style_closed = 2
      nPts = Radii * 4
      a = FINDGEN(nPts) * ( (2 * !PI) / (nPts-1) )
      oROI->GetProperty, DATA = data
      xImage = (MAX(data[0,*])+MIN(data[0,*]))/2.
      yImage = (MAX(data[1,*])+MIN(data[1,*]))/2.
      newX = COS(a) * Radii + xImage
      newY = SIN(a) * Radii + yImage
      newZ = REPLICATE(0.0, nPts)
      style = style_closed

      oROI->GetProperty, N_VERTS=nVerts
      oROI->ReplaceData, newX, newY, newZ, START=0, FINISH=nVerts-1
      oROI->SetProperty, STYLE=style
      ;      xd_xroi__DeleteSelectedROI, pState
      ;      (*pState).oWindow->Draw, (*pState).oView
      (*pState).oCurrROI = oROI
      ;      (*pState).oModel->Add, oROI
      (*pState).oWindow->Draw, (*pState).oView
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
    end



    ;    'LAYERMODE': begin
    ;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
    ;      case sEvent.value of
    ;        0: (*pState).n_layers = 1
    ;        1: (*pState).n_layers = 5
    ;        2: (*pState).n_layers = 11
    ;        3: (*pState).n_layers = 21
    ;        4: (*pState).n_layers = 51
    ;        5: (*pState).n_layers = 101
    ;        6: (*pState).n_layers = 201
    ;        7: (*pState).n_layers = (*pState).xd_sState.imageSize
    ;      endcase
    ;;      print, (*pState).n_layers
    ;    end

    'TRANSLATE-SCALE': begin
      ; Translate/Scale tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'TRANSLATE-SCALE' then begin
        (*pState).mode = 'TRANSLATE-SCALE'

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif

        ; Set the translate/scale selection visual as current.
        (*pState).oSelVisual = (*pState).oTransScaleVisual
        xd_xroi__ReshapeSelectionVisual, pState, (*pState).oSelROI

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'RECTANGLE': begin
      ; Rectangle ROI tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'RECTANGLE' then begin
        (*pState).mode = 'RECTANGLE'

        ; Get current list of ROI names.
        oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
        if (nROIs gt 0) then begin
          for i=0, nROIs-1 do begin
            (*pState).oROIModel->Remove, oROIs[i]
            (*pState).oROIGroup->Remove, oROIs[i]
            if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
            OBJ_DESTROY, oROIs[i]
          endfor
        endif

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'ELLIPSE': begin
      ; Ellipse ROI tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'ELLIPSE' then begin
        (*pState).mode = 'ELLIPSE'

        ; Get current list of ROI names.
        oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
        if (nROIs gt 0) then begin
          for i=0, nROIs-1 do begin
            (*pState).oROIModel->Remove, oROIs[i]
            (*pState).oROIGroup->Remove, oROIs[i]
            if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
            OBJ_DESTROY, oROIs[i]
          endfor
        endif

        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'FREEPOLY': begin
      ; Freehand ROI tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'FREEHAND DRAW' then begin
        (*pState).mode = 'FREEHAND DRAW'

        ; Get current list of ROI names.
        oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
        if (nROIs gt 0) then begin
          for i=0, nROIs-1 do begin
            (*pState).oROIModel->Remove, oROIs[i]
            (*pState).oROIGroup->Remove, oROIs[i]
            if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
            OBJ_DESTROY, oROIs[i]
          endfor

          (*pState).oWindow->Draw, (*pState).oView
          ; Pick a new selected region.
          ;              WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
          ;              WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
          xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
        endif


        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'CIRCULARMASKPOLY': begin
      ; Freehand ROI tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'CIRCULARMASK DRAW' then begin
        (*pState).mode = 'CIRCULARMASK DRAW'

        ; Get current list of ROI names.
        oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
        if (nROIs gt 0) then begin
          for i=0, nROIs-1 do begin
            (*pState).oROIModel->Remove, oROIs[i]
            (*pState).oROIGroup->Remove, oROIs[i]
            if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
            OBJ_DESTROY, oROIs[i]
          endfor
        endif


        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'SEGPOLY': begin
      ; Segmented ROI tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'POLYGON DRAW' then begin
        (*pState).mode = 'POLYGON DRAW'

        ; Get current list of ROI names.
        oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
        if (nROIs gt 0) then begin
          for i=0, nROIs-1 do begin
            (*pState).oROIModel->Remove, oROIs[i]
            (*pState).oROIGroup->Remove, oROIs[i]
            if OBJ_VALID((*pState).oRegionsOut) then (*pState).oRegionsOut->Remove, oROIs[i]
            OBJ_DESTROY, oROIs[i]
          endfor
        endif

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif


        (*pState).oImage0->GetProperty, DATA=Mask_ori
        (*pState).oImage1->SetProperty, DATA=Mask_ori
        (*pState).oImage2->SetProperty, DATA=Mask_ori*0.
        (*pState).new_mask_vol_cube = (*pState).original_mask_vol_cube

        (*pState).bTempSegment = 1b
        oROI = (*pState).oCurrROI
        oROI->setProperty, ALPHA_CHANNEL=1.0
        oROI->setProperty, STYLE=1
        ;            print, 'style = ', style
        oROI->GetProperty, N_VERTS = nv
        oROI->RemoveData, Count= nv
        (*pState).oWindow->Draw, (*pState).oView
      endif
    end



    'HOTKEY' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      case STRUPCASE(sEvent.ch) of
        '-': begin ; unconstrained rotation
          WIDGET_CONTROL, (*pState).wThresholdSlider_radii, GET_VALUE=Radii
          IF Radii GT 2 THEN BEGIN
            (*pState).Radii = Radii-1
            ;print, 'Radii = ', Radii-1
            WIDGET_CONTROL, (*pState).wThresholdSlider_radii, SET_VALUE=(*pState).Radii
          ENDIF
        end
        '=': begin
          WIDGET_CONTROL, (*pState).wThresholdSlider_radii, GET_VALUE=Radii
          IF Radii LT 50 THEN BEGIN
            (*pState).Radii = Radii+1
            ;print, 'Radii = ', Radii+1
            WIDGET_CONTROL, (*pState).wThresholdSlider_radii, SET_VALUE=(*pState).Radii
          ENDIF
          ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
        end
        '+': begin
          WIDGET_CONTROL, (*pState).wThresholdSlider_radii, GET_VALUE=Radii
          IF Radii LT 50 THEN BEGIN
            (*pState).Radii = Radii+1
            ;print, 'Radii = ', Radii+1
            WIDGET_CONTROL, (*pState).wThresholdSlider_radii, SET_VALUE=(*pState).Radii
          ENDIF
          ;WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
        end
        else:
      endcase
      oROI = (*pState).oCurrROI
      oROI->GetProperty, DATA = data
      style_closed = 2
      nPts = (*pState).Radii * 4
      a = FINDGEN(nPts) * ( (2 * !PI) / (nPts-1) )
      xImage = (MAX(data[0,*])+MIN(data[0,*]))/2.
      yImage = (MAX(data[1,*])+MIN(data[1,*]))/2.
      ;            print, 'xMax = ', MAX(data[0,*]), ', yMax = ', MAX(data[1,*])
      ;            print, 'xMin = ', MIN(data[0,*]), ', yMin = ', MIN(data[1,*])
      ;            print, 'xImage = ', xImage, ', yImage = ', yImage
      newX = COS(a) * (*pState).Radii + xImage
      newY = SIN(a) * (*pState).Radii + yImage
      newZ = REPLICATE(0.0, nPts)
      style = style_closed

      oROI->GetProperty, N_VERTS=nVerts
      oROI->ReplaceData, newX, newY, newZ, START=0, FINISH=nVerts-1
      oROI->SetProperty, STYLE=style
      (*pState).oCurrROI = oROI
      (*pState).oWindow->Draw, (*pState).oView

      ; If button down, append a vertex.
      if ((*pState).bButtonDown NE 0) then begin
        ; Set the region as current.
        xd_xroi__SetROI, pState, oROI, /UPDATE_LIST, $
          /SET_LIST_SELECT
        oROI->SetProperty, STYLE=2
        ;        modify_mask_by_ROI, pState
        (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac
      endif else begin
        (*pState).oWindow->Draw, (*pState).oView
      endelse

    end



    'PICK': begin
      ; Pick tool selected.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      if (*pState).mode ne 'SELECTION' then begin
        (*pState).mode = 'SELECTION'

        ; Disable old selection visual, if any.
        oSelVisual = (*pState).oSelVisual
        if (OBJ_VALID(oSelVisual) ne 0) then begin
          oSelVisual->SetProperty, /HIDE
          (*pState).oSelVisual = OBJ_NEW()
        endif

        ; Set the vertex picking visual as the current selection visual.
        (*pState).oSelVisual = (*pState).oPickVisual
        xd_xroi__ReshapeSelectionVisual, pState, (*pState).oSelROI

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'FLIP': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).oImage->GetProperty, ORDER=order
      (*pState).oImage->SetProperty, ORDER=1-KEYWORD_SET(order)
      (*pState).oWindow->Draw, (*pState).oView
    end

    'IMPORT': begin
      nDims = 0
      while nDims eq 0 do begin
        if not DIALOG_READ_IMAGE( $
          RED=red, $
          GREEN=green, $
          BLUE=blue, $
          DIALOG_PARENT=sEvent.top, $
          IMAGE=image, $
          QUERY=query, $
          TITLE='Import Image File' $
          ) $
          then $
          RETURN
        if image[0] eq -1 then $
          RETURN ; User "opened" non-image file or directory.

        ; Verify that the selected image has valid dimensions.
        nDims = SIZE(image, /N_DIMENSIONS)
        case nDims of
          2: begin ; 8-bit image
            dimensions = SIZE(image, /DIMENSIONS)
            image_is_8bit = 1b
          end

          3: begin ; RGB image
            allDims = SIZE(image, /DIMENSIONS)
            iInterleave = (WHERE(allDims eq 3))[0]
            case iInterleave of
              0: dimensions = allDims[1:2]
              1: dimensions = [allDims[0],allDims[2]]
              2: dimensions = allDims[0:1]
              else: begin
                void = DIALOG_MESSAGE(/INFORM, $
                  ['Image must have dimensions', $
                  '[3,m,n], [m,3,n], or [m,n,3].'], $
                  DIALOG_PARENT=sEvent.top)
                nDims = 0
                image = 0
              end
            endcase
            image_is_8bit = 0b
          end

          else: begin
            void = DIALOG_MESSAGE(/INFORM, $
              'Image data must have 2 or 3 dimensions.', $
              DIALOG_PARENT=sEvent.top)
            nDims = 0
            image = 0
          end
        endcase
      endwhile

      ; The imported image has been accepted.  Store it and update.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      (*pState).oModel->Remove, (*pState).oImage
      OBJ_DESTROY, (*pState).oImage
      (*pState).image_is_8bit = image_is_8bit
      if image_is_8bit then $
        *(*pState).pImg = image
      (*pState).oImage = OBJ_NEW('IDLgrImage', TEMPORARY(image), $
        INTERLEAVE=iInterleave)

      if (image_is_8bit and (query.has_palette eq 0)) then begin
        query.has_palette = 1
        red = BINDGEN(256)
        green = red
        blue = red
      endif

      if query.has_palette then begin
        (*pState).oPalette->SetProperty, $
          RED_VALUES=red, $
          GREEN_VALUES=green, $
          BLUE_VALUES=blue
        if (image_is_8bit and $
          WIDGET_INFO((*pState).wLoadCT, /VALID_ID)) then begin
          WIDGET_CONTROL, (*pState).wPaletteEdit, $
            SET_VALUE=(*pState).oPalette
          WIDGET_CONTROL, (*pState).wLoadCT, GET_UVALUE=p
          (*p).red_values = red
          (*p).green_values = green
          (*p).blue_values = blue
        endif
      endif
      (*pState).oImage->SetProperty, PALETTE=(*pState).oPalette

      if image_is_8bit then begin
        WIDGET_CONTROL, (*pState).wLoadCTButton, SENSITIVE=1
      endif else begin
        if WIDGET_INFO((*pState).wLoadCT, /VALID) then $
          WIDGET_CONTROL, (*pState).wLoadCT, /DESTROY
        WIDGET_CONTROL, (*pState).wLoadCTButton, SENSITIVE=0
      endelse


      ; Change the virtual canvas.
      WIDGET_CONTROL, (*pState).wDraw, $
        DRAW_XSIZE=query.dimensions[0], $
        DRAW_YSIZE=query.dimensions[1]


      ; Shrink the viewport if the image is smaller.
      draw_geom = WIDGET_INFO((*pState).wDraw, /GEOM)
      newXsize = (query.dimensions[0] < draw_geom.xsize)
      ; If new image is bigger than toolbar, make sure draw window is also that big.
      if (query.dimensions[0] gt (*pState).toolbar_xsize) then $
        newXsize = newXsize > (*pState).toolbar_xsize
      newYsize = query.dimensions[1] < draw_geom.ysize

      ; Change the viewport.
      WIDGET_CONTROL, (*pState).wDraw, XSIZE=newXsize, YSIZE=newYsize


      ; Retrieve the new base size and cache it.
      base_geom = WIDGET_INFO(sEvent.top, /GEOMETRY)
      (*pState).scr_xsize = base_geom.scr_xsize
      (*pState).scr_ysize = base_geom.scr_ysize

      ; Get the new geometry.
      draw_geom = WIDGET_INFO((*pState).wDraw, /GEOM)

      ; Set the new viewplane rect, depending upon the new viewport.
      WIDGET_CONTROL, (*pState).wDraw, GET_DRAW_VIEW=viewport

      (*pState).oView->SetProperty, $
        VIEWPLANE_RECT=[viewport[0], viewport[1], $
        draw_geom.xsize < draw_geom.draw_xsize, draw_geom.ysize]

      (*pState).oModel->Add, (*pState).oImage, POSITION=0
      (*pState).oWindow->Draw, (*pState).oView
      if OBJ_VALID((*pState).oSelROI) then begin
        xd_xroi__SetROI, pState, (*pState).oSelROI
      endif
    end

    'LOADCT': begin
      xd_xroiLoadCT, sEvent.top
    end

    'PICKCOLOR': begin
      xd_xroiPickColor, sEvent.top
    end

    'ROI_INFO': begin
      ; ROI Info menu item selected.  Open the ROI Info dialog.
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
    end

    'ROI_DELETE': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      xd_xroi__DeleteSelectedROI, pState

      (*pState).oWindow->Draw, (*pState).oView
    end

    'ROI_HISTOGRAM': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      xd_xroi__HistogramSelectedROI, pState, GROUP=sEvent.top
    end

    'ROI_GROW_BY_THRESHOLD': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      xd_xroi__GrowSelectedROI, pState

      (*pState).oWindow->Draw, (*pState).oView
    end

    'ROI_GROW_BY_STDDEV': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      xd_xroi__GrowSelectedROI, pState, /STDDEV

      (*pState).oWindow->Draw, (*pState).oView
    end

    'ROI_GROW_PROPERTIES': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      xd_xroiGrowProps, pState, GROUP_LEADER=sEvent.top
    end

    'SLIDER': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, (*pState).wSlider, GET_VALUE=location_value
      WIDGET_CONTROL, (*pState).wThresholdSlider_min, GET_VALUE=threshold_min
      WIDGET_CONTROL, (*pState).wThresholdSlider_max, GET_VALUE=threshold_max

      imageSize = (*pState).imageSize * 3
      ;vol_HU_cube = (*pState).vol_HU_cube
      IF STRCMP((*pState).direction, 'xy') THEN BEGIN
        xyImage = DBLARR(imageSize, imageSize)
        (*pState).location3D[2] = location_value - 1
        location3D = (*pState).location3D
        temp_xyImage = CONGRID((*pState).svol_HU_cube[*,*, location3D[2]], imageSize, imageSize, 1)
        xyImage[*, *] = temp_xyImage[*, *, 0]
        Img = REVERSE(xyImage, 2)
      ENDIF

      IF STRCMP((*pState).direction, 'yz') THEN BEGIN
        yzImage = DBLARR(imageSize, imageSize)
        (*pState).location3D[0] = location_value - 1
        location3D = (*pState).location3D
        temp_yzImage = CONGRID((*pState).svol_HU_cube[location3D[0], *,*], 1, imageSize, imageSize)
        yzImage[*, *] = temp_yzImage[0, *, *]
        Img = yzImage
      ENDIF

      IF STRCMP((*pState).direction, 'xz') THEN BEGIN
        xzImage = DBLARR(imageSize, imageSize)
        (*pState).location3D[1] = location_value - 1
        location3D = (*pState).location3D
        temp_xzImage = CONGRID((*pState).svol_HU_cube[*, location3D[1], *], imageSize, 1, imageSize)
        xzImage[*, *] = temp_xzImage[*, 0, *]
        Img = xzImage
      ENDIF

      ;;;;;;;;;;;;;;;;;;;;;;hjx1129
      IF STRCMP((*pState).direction, 'label') THEN BEGIN
        xxImage = DBLARR(imageSize, imageSize)
        (*pState).location3D[1] = location_value - 1
        location3D = (*pState).location3D
        temp_xxImage = CONGRID((*pState).svol_HU_cube[*, location3D[1], *], imageSize, 1, imageSize)
        xxImage[*, *] = temp_xxImage[*, 0, *]
        Img = xxImage
      ENDIF
      ;;;;;;;;;;;;;;;;;;;;;;hjx1129
      ;
      ;;;;;;;;;;;;;;;;;;;;;;hjx1212
      IF STRCMP((*pState).direction, 'distance') THEN BEGIN
        wwImage = DBLARR(imageSize, imageSize)
        (*pState).location3D[1] = location_value - 1
        location3D = (*pState).location3D
        Img = wwImage
      ENDIF
      ;;;;;;;;;;;;;;;;;;;;;;hjx1212

      ;
      *((*pState).pImg) = Img
      (*pState).oImage->SetProperty, DATA=BYTSCL(*((*pState).pImg), MAX=threshold_max, MIN=threshold_min)

      currdim = SIZE((*pState).svol_HU_cube, /DIMENSIONS)
      maskdim = SIZE((*pState).new_mask_vol_cube, /DIMENSIONS)
      Mask = DBLARR(imageSize, imageSize)


      IF STRCMP((*pState).direction, 'xy') THEN BEGIN
        Mask3D = REVERSE((CONGRID((*pState).new_mask_vol_cube[*,*,location3D[2]*maskdim[2]/currdim[2]], $
          imageSize, imageSize,1) GT 0)*250, 2)
        Mask[*, *] = Mask3d[*, *, 0]
        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
        (*pState).oImage0->SetProperty, DATA=Mask
        (*pState).oImage1->SetProperty, DATA=Mask
        (*pState).oImage2->SetProperty, DATA=Mask*0.
        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129

      ENDIF

      IF STRCMP((*pState).direction, 'yz') THEN BEGIN
        Mask3D = (CONGRID((*pState).new_mask_vol_cube[location3D[0]*maskdim[0]/currdim[0],*,*], 1, $
          imageSize, imageSize) GT 0)*140
        Mask[*, *] = Mask3d[0, *, *]
        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
        (*pState).oImage0->SetProperty, DATA=Mask
        (*pState).oImage1->SetProperty, DATA=Mask
        (*pState).oImage2->SetProperty, DATA=Mask*0.
        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
      ENDIF

      IF STRCMP((*pState).direction, 'xz') THEN BEGIN
        Mask3D = (CONGRID((*pState).new_mask_vol_cube[*,location3D[1]*maskdim[1]/currdim[1],*], $
          imageSize, 1, imageSize) GT 0)*210
        Mask[*, *] = Mask3d[*, 0, *]
        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
        (*pState).oImage0->SetProperty, DATA=Mask
        (*pState).oImage1->SetProperty, DATA=Mask
        (*pState).oImage2->SetProperty, DATA=Mask*0.
        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
      ENDIF
      ;        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
      IF STRCMP((*pState).direction, 'label') THEN BEGIN
        Mask3D = (CONGRID((*pState).new_mask_vol_cube[*,location3D[1]*maskdim[1]/currdim[1],*], $
          imageSize, 1, imageSize) GT 0)*210
        Mask[*, *] = Mask3d[*, 0, *]
        label=label_region(Mask[*,*])
        Mask_size=size(Mask)
        new_mask=intarr(Mask_size[1],Mask_size[2])
        image_col=Mask_size[1]
        loadct,27
        if (max(label) NE 0)then begin

          for j=1,max(label) do begin
            location=where(label eq j)
            col=location mod image_col
            row=location/image_col
            a=255
            mask[min(col):max(col),max(row)]=a
            mask[min(col):max(col),min(row)]=a
            mask[min(col),min(row):max(row)]=a
            mask[max(col),min(row):max(row)]=a

          endfor
        endif
        (*pState).oImage0->SetProperty, DATA=Mask
        (*pState).oImage1->SetProperty, DATA=mask
        (*pState).oImage2->SetProperty, DATA=Mask*0.

      ENDIF
      ;        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
      ;
      ;        ;;;;;;;;;;;;;;;;;;;;;;;;hjx1212
      IF STRCMP((*pState).direction, 'distance') THEN BEGIN

        Mask3D = (CONGRID((*pState).new_mask_vol_cube[*,location3D[1]*maskdim[1]/currdim[1],*], $
          imageSize, 1, imageSize) GT 0)*210
        Mask[*, *] = Mask3d[*, 0, *]

        mask_size=size(mask)
        image_col=Mask_size[1]
        image_row=Mask_size[2]
        w=[[1,1,1],[1,1,1],[1,1,1]]
        minnum=min(mask)
        new_mask=intarr(image_col,image_row)

        while(total(mask) ne image_col*image_row*minnum) do begin
          for i=0,image_col-1 do begin
            for j=0,image_row-1 do begin
              if (mask[i,j] ne minnum) then begin
                new_mask[i,j]=new_mask[i,j]+20
              endif
            endfor
          endfor
          mask=ERODE(mask,w)


        endwhile

        (*pState).oImage0->SetProperty, DATA=new_Mask
        (*pState).oImage1->SetProperty, DATA=new_mask
        (*pState).oImage2->SetProperty, DATA=new_mask*0.

      ENDIF

      ;;;;;;;;;;;;;;;;;;;;;;;;hjx1212

      (*pState).oWindow->Draw, (*pState).oView
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=pState, /NO_COPY
    end

    else: begin
    end

  endcase


  if xregistered('demo_tour') eq 0 then begin
    widget_control, sEvent.top, get_uvalue=pState, /no_copy
    widget_control, (*pState).wHotKeyReceptor, /input_focus
    widget_control, sEvent.top, set_uvalue=pState, /no_copy
  end


end

;------------------------------------------------------------------------------
pro xd_xroi__ButtonPress, sEvent
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  ; Convert from viewport coordinates to image coordinates.
  WIDGET_CONTROL, (*pState).wDraw, GET_DRAW_VIEW=viewport
  xImage = sEvent.x + viewport[0]
  yImage = sEvent.y + viewport[1]

  case (*pState).mode of
    'TRANSLATE-SCALE': begin
      if (sEvent.press eq 1) then begin  ; Left mouse button.
        ; Temporarily hide the image for selection.
        (*pState).oImage->SetProperty, HIDE=1

        ; Temporarily hide all regions for selection.
        (*pState).oROIModel->SetProperty, /HIDE

        ; Check if a selection visual was hit.
        oSel = (*pState).oWindow->Select((*pState).oView, $
          [sEvent.x, sEvent.y])
        selType = SIZE(oSel, /TYPE)
        if (selType eq 11) then begin  ; Object reference?
          ; A selection visual handle was hit.
          (*pState).oSelHandle = oSel[0]

          oSel[0]->GetProperty, NAME=handleName
          case handleName of
            'SCALE_LL': bSaveROIData = 1b
            'SCALE_LR': bSaveROIData = 1b
            'SCALE_UL': bSaveROIData = 1b
            'SCALE_UR': bSaveROIData = 1b
            else: bSaveROIData = 0b
          endcase
          if (bSaveROIData ne 0) then begin
            oROI = (*pState).oSelROI
            if (OBJ_VALID(oROI) ne 0) then begin
              oROI->GetProperty, DATA=roiData, $
                ROI_XRANGE=xrange, ROI_YRANGE=yrange
              ; Translate to origin.
              roiData[0,*] = roiData[0,*] - xrange[0]
              roiData[1,*] = roiData[1,*] - yrange[0]
              (*pState).pSavedROIData = PTR_NEW(roiData)
              (*pState).savedROIXRange = xrange
              (*pState).savedROIYRange = yrange
            endif else $
              (*pState).pSavedROIData = PTR_NEW()
          endif
        endif else $
          (*pState).oSelHandle = OBJ_NEW()

        ; Restore the regions.
        (*pState).oROIModel->SetProperty, HIDE=0

        ; If a selection visual was not hit...
        if (OBJ_VALID((*pState).oSelHandle) eq 0) then begin

          ; Temporarily hide the selection visual.
          (*pState).oTransScaleVisual->SetProperty, /HIDE

          ; Check if a region was hit.
          oSel = (*pState).oWindow->Select((*pState).oView, $
            DIMENSIONS=[16,16], $
            [sEvent.x, sEvent.y])

          selType = SIZE(oSel, /TYPE)
          if (selType eq 11) then begin  ; Object reference?
            oROI = oSel[0]

            ; Mark the region as selected.
            xd_xroi__SetROI, pState, oROI, /SET_LIST_SELECT
          endif else $
            ; Mark no regions as being currently selected.
            xd_xroi__SetROI, pState, OBJ_NEW()
        endif

        ; Restore the image.
        (*pState).oImage->SetProperty, HIDE=0

        (*pState).bButtonDown = 1b
        (*pState).buttonXY = [sEvent.x, sEvent.y]

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'RECTANGLE': begin
      if (sEvent.press eq 1) then begin  ; Left mouse button.
        ; Change the color of the previously selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

        ; Create a new rectangle region.
        oROI = OBJ_NEW('IDLgrROI', $
          COLOR=(*pState).sel_rgb, $
          STYLE=0 $
          )

        (*pState).oCurrROI = oROI
        (*pState).oModel->Add, oROI

        ; Set initial corner for the rectangle.
        oROI->AppendData, [xImage, yImage, 0]

        (*pState).oWindow->Draw, (*pState).oView

        (*pState).bButtonDown = 1b
        (*pState).buttonXY = [xImage, yImage]
      endif
    end

    'ELLIPSE': begin
      if (sEvent.press eq 1) then begin  ; Left mouse button.
        ; Change the color of the previously selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

        ; Create a new ellipse region.
        oROI = OBJ_NEW('IDLgrROI', $
          COLOR=(*pState).sel_rgb, $
          STYLE=0 $
          )

        (*pState).oCurrROI = oROI
        (*pState).oModel->Add, oROI

        ; Initialize the ellipse as a single point.
        oROI->AppendData, [xImage, yImage, 0]

        (*pState).oWindow->Draw, (*pState).oView

        (*pState).bButtonDown = 1b
        (*pState).buttonXY = [xImage, yImage]
      endif
    end

    'FREEHAND DRAW': begin
      if (sEvent.press eq 1) then begin
        oROI = (*pState).oCurrROI
        if (OBJ_VALID(oROI) eq 0) then begin
          oOldSelROI = (*pState).oSelROI
          if (OBJ_VALID(oOldSelROI) ne 0) then $
            oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

          oROI = OBJ_NEW('IDLgrROI', $
            COLOR=(*pState).sel_rgb, $
            STYLE=1 $
            )

          (*pState).oCurrROI = oROI
          (*pState).oModel->Add, oROI
        endif

        oROI->setProperty, ALPHA_CHANNEL=1.0
        oROI->GetProperty, N_VERTS = nv
        oROI->RemoveData, Count= nv
        oROI->AppendData, [xImage, yImage, 0]

        (*pState).oWindow->Draw, (*pState).oView
        (*pState).bButtonDown = sEvent.clicks
      endif
    end

    'CIRCULARMASK DRAW': begin
      if (sEvent.press eq 1) then begin
        oROI = (*pState).oCurrROI
        if (OBJ_VALID(oROI) eq 0) then begin
          oOldSelROI = (*pState).oSelROI
          if (OBJ_VALID(oOldSelROI) ne 0) then $
            oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

          oROI = OBJ_NEW('IDLgrROI', $
            COLOR=(*pState).sel_rgb, $
            STYLE=1 $
            )

          (*pState).oCurrROI = oROI
          (*pState).oModel->Add, oROI
        endif

        ;            oROI->AppendData, [xImage, yImage, 0]
        oROI = (*pState).oCurrROI
        if (OBJ_VALID(oROI) EQ 0) then return
        style_closed = 2
        Radii = (*pState).Radii
        nPts = Radii * 4
        a = FINDGEN(nPts) * ( (2 * !PI) / (nPts-1) )
        newX = COS(a) * Radii + xImage
        newY = SIN(a) * Radii + yImage
        newZ = REPLICATE(0.0, nPts)
        style = style_closed
        oROI->GetProperty, N_VERTS=nVerts
        oROI->ReplaceData, newX, newY, newZ, START=0, FINISH=nVerts-1
        oROI->SetProperty, STYLE=style

        (*pState).bButtonDown = sEvent.clicks
        print, '(*pState).bButtonDown = ', (*pState).bButtonDown
        ; If button down, append a vertex.
        if ((*pState).bButtonDown NE 0) then begin
          ; Set the region as current.
          xd_xroi__SetROI, pState, oROI, /UPDATE_LIST, $
            /SET_LIST_SELECT
          oROI->SetProperty, STYLE=2

          ;get the size of oImage
          (*pState).oImage->getProperty, DATA=Img
          roi_mask_dims = SIZE(Img, /DIMENSIONS)
          mask_ROI = (*pState).oSelROI -> ComputeMask(DIMENSIONS = roi_mask_dims)
          ;update the mask
          (*pState).oImage1->GetProperty, DATA=mask_old
          IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN mask = mask_old OR mask_ROI
          IF STRCMP((*pState).Modify_ROI, 'ADD') THEN mask = mask_old OR mask_ROI
          IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN mask = ((mask_old GT 100) AND (1-(mask_ROI GT 100)))*255
          
          (*pState).oImage1->SetProperty, DATA=Mask
          (*pState).oImage2->getProperty, DATA=mask_ROI_old
          (*pState).oImage2->setProperty, DATA=mask_ROI_old+mask_ROI

          ;define alpha channel for storing location info
          ac = 1.0
          IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN ac = 0.000012
          IF STRCMP((*pState).Modify_ROI, 'ADD') THEN ac = 0.000004
          IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN ac = 0.000008

          IF STRCMP((*pState).direction, 'xy') THEN Location = (*pState).location3D[2]
          IF STRCMP((*pState).direction, 'yz') THEN Location = (*pState).location3D[0]
          IF STRCMP((*pState).direction, 'xz') THEN Location = (*pState).location3D[1]

          aac=Location/10000.+ac
          (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
          (*pState).oWindow->Draw, (*pState).oView
          (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac
        endif else begin
          (*pState).oWindow->Draw, (*pState).oView
        endelse


      endif
    end

    ; Segmented ROI
    'POLYGON DRAW': begin
      if (sEvent.press eq 1) then begin
        oROI = (*pState).oCurrROI
        if (OBJ_VALID(oROI) eq 0) then begin
          oOldSelROI = (*pState).oSelROI
          if (OBJ_VALID(oOldSelROI) ne 0) then $
            oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

          oROI = OBJ_NEW('IDLgrROI', $
            COLOR=(*pState).sel_rgb, $
            STYLE=1 $
            )
          (*pState).oCurrROI = oROI
          (*pState).oModel->Add, oROI
          (*pState).oWindow->Draw, (*pState).oView
        endif

        if (sEvent.clicks eq 1) then begin
          ; If dragging a temporary segment, start a new one.
          ; Otherwise, append a new vertex.
          if ((*pState).bTempSegment eq 1) then $
            (*pState).bTempSegment = 0 $
          else $
            oROI->AppendData, [xImage, yImage, 0]
        endif

        (*pState).bButtonDown = sEvent.clicks
      endif
    end

    ; Pick ROI.
    'SELECTION': begin
      (*pState).oImage->SetProperty, HIDE=1
      oSel = (*pState).oWindow->Select((*pState).oView, $
        DIMENSIONS=[32,32], $
        [sEvent.x, sEvent.y])
      (*pState).oImage->SetProperty, HIDE=0
      selType = SIZE(oSel, /TYPE)
      if (selType eq 11) then begin  ; Object reference?
        oROI = oSel[0]
        xd_xroi__SetROI, pState, oROI, /SET_LIST_SELECT
      endif
      (*pState).bButtonDown = 1b

      if (OBJ_VALID((*pState).oSelROI) ne 0) then begin
        ; Pick nearest vertex in ROI.
        oROI = (*pState).oSelROI

        vertIndx = oROI->PickVertex((*pState).oWindow, $
          (*pState).oView, $
          [sEvent.x, sEvent.y])
        if (vertIndx ge 0) then begin

          ; Show the pick vertex.
          oPickPolyline = (*pState).oPickVisual->Get()
          oPickPolyline->SetProperty, HIDE=0

          (*pState).oPickVisual->Reset
          oROI->GetProperty, DATA=vertData
          selVert = vertData[*,vertIndx]
          (*pState).oPickVisual->Translate, $
            selVert[0], selVert[1], selVert[2]

          (*pState).oWindow->Draw, (*pState).oView
        endif
      endif
      if (selType eq 11) then begin ; Object reference?
        if not (*pState).modal then $
          xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
      endif
    end
  endcase
end

;------------------------------------------------------------------------------
function xd_xroi__GenerateName, state
  COMPILE_OPT hidden
  ;
  ;Purpose: generate a unique string for naming an ROI.
  ;
  oROIs = state.oROIModel->Get(/ALL, COUNT=count)
  if count gt 0 then begin
    current_names = STRARR(count)
    for i=0,count-1 do begin
      oROIs[i]->GetProperty, NAME=name
      current_names[i] = name
    endfor
  endif

  done = 0b
  id = 0
  while not done do begin
    id = id + 1
    name = 'Region '+ STRTRIM(id, 2)
    name_is_duplicate = 0b
    for i=0,N_ELEMENTS(current_names)-1 do begin
      if STRUPCASE(STRTRIM(current_names[i], 2)) $
        eq STRUPCASE(name) then begin
        name_is_duplicate = 1b
      endif
    end
    if not name_is_duplicate then begin
      done = 1b
    end
  endwhile

  return, name
end

;------------------------------------------------------------------------------
pro xd_xroi__ButtonRelease, sEvent
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  ; Convert from viewport coordinates to image coordinates.
  WIDGET_CONTROL, (*pState).wDraw, GET_DRAW_VIEW=viewport
  xImage = sEvent.x + viewport[0]
  yImage = sEvent.y + viewport[1]


  case (*pState).mode of
    'TRANSLATE-SCALE': begin
      if (sEvent.release ne 1) then break  ; Left mouse button.

      ; Restore region statistics.
      if ((*pState).bButtonDown eq 2) then begin
        ; Determine if ROI Info dialog is currently realized.
        haveInfoBase = WIDGET_INFO((*pState).wROIInfo, /VALID_ID)
        if (haveInfoBase NE 0) then begin
          WIDGET_CONTROL, (*pState).wROIInfo, GET_UVALUE=sState

          xd_xroiInfo__SetStats, sState, (*pState).oSelROI

          modify_mask_by_ROI, pState
          (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
          (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
          (*pState).oWindow->Draw, (*pState).oView
          (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac
        endif
      endif

      ; Free any pointers to saved ROI data.
      if (PTR_VALID((*pState).pSavedROIData) ne 0) then $
        PTR_FREE, (*pState).pSavedROIData
      (*pState).pSavedROIData = PTR_NEW()

      ; Restore button state.
      (*pState).bButtonDown = 0b
    end

    'RECTANGLE': begin
      if (sEvent.release ne 1) then break  ; Left mouse button.
      if ((*pState).bButtonDown ne 1) then break  ; button was down

      ; Reset button down state.
      (*pState).bButtonDown = 0b

      oROI = (*pState).oCurrROI
      if (not OBJ_VALID(oROI)) then break

      ; Ensure that the rectangle has at 4 vertices.
      oROI->GetProperty, DATA=roiData
      if ((N_ELEMENTS(roiData)/3) eq 4) then begin
        ; The rectangle region is valid.  Give it a name
        ; and add it to the appropriate containers.
        oROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oModel->Remove, oROI
        (*pState).oROIModel->Add, oROI
        (*pState).oROIGroup->Add, oROI
        if OBJ_VALID((*pState).oRegionsOut) then $
          (*pState).oRegionsOut->Add, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ; Activate appropriate tool buttons.
        ;                if lmgr(/demo) ne 1 then begin
        ;                    WIDGET_CONTROL, (*pState).wSaveButton, $
        ;                        SENSITIVE=1
        ;                    WIDGET_CONTROL, (*pState).wSaveToolButton, $
        ;                        SENSITIVE=1
        ;                endif

        ; Set the region as current.
        xd_xroi__SetROI, pState, oROI, /UPDATE_LIST, $
          /SET_LIST_SELECT

        oROI->SetProperty, STYLE=2

        modify_mask_by_ROI, pState
        (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac

        ; If this is the first region, bring up the
        ; region information dialog.
        if ((*pState).bFirstROI eq 1b) then begin
          WIDGET_CONTROL, /HOURGLASS
          (*pState).bFirstROI = 0b
          xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
        endif
      endif else begin
        ; Fewer than 4 vertices; delete.
        (*pState).oModel->Remove, oROI
        OBJ_DESTROY, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ; Reset color of formerly selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, $
          COLOR=(*pState).sel_rgb

        (*pState).oWindow->Draw, (*pState).oView
      endelse

    end

    'ELLIPSE': begin
      if (sEvent.release ne 1) then break  ; Left mouse button only.
      if ((*pState).bButtonDown ne 1) then break  ; button was down

      ; Reset button down state.
      (*pState).bButtonDown = 0b

      oROI = (*pState).oCurrROI
      if (not OBJ_VALID(oROI)) then break

      ; Ensure that the ellipse has at least 4 vertices.
      oROI->GetProperty, DATA=roiData
      if ((N_ELEMENTS(roiData)/3) ge 4) then begin
        ; The ellipse region is valid.  Give it a name
        ; and add it to the appropriate containers.
        oROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oModel->Remove, oROI
        (*pState).oROIModel->Add, oROI
        (*pState).oROIGroup->Add, oROI
        if OBJ_VALID((*pState).oRegionsOut) then $
          (*pState).oRegionsOut->Add, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ;                ; Activate appropriate tool buttons.
        ;                if lmgr(/demo) ne 1 then begin
        ;                    WIDGET_CONTROL, (*pState).wSaveButton, $
        ;                        SENSITIVE=1
        ;                    WIDGET_CONTROL, (*pState).wSaveToolButton, $
        ;                        SENSITIVE=1
        ;                endif

        ; Set the region as current.
        xd_xroi__SetROI, pState, oROI, /UPDATE_LIST, $
          /SET_LIST_SELECT

        oROI->SetProperty, STYLE=2

        modify_mask_by_ROI, pState
        (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac

        ; If this is the first region, bring up the
        ; region information dialog.
        if ((*pState).bFirstROI eq 1b) then begin
          WIDGET_CONTROL, /HOURGLASS
          (*pState).bFirstROI = 0b
          xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
        endif
      endif else begin
        ; Fewer than 4 vertices; delete.
        (*pState).oModel->Remove, oROI
        OBJ_DESTROY, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ; Reset color of formerly selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, $
          COLOR=(*pState).sel_rgb

        (*pState).oWindow->Draw, (*pState).oView
      endelse
    end

    ; Freehand ROI
    'FREEHAND DRAW': begin
      if (sEvent.release ne 1) then break
      if ((*pState).bButtonDown ne 1) then break  ; Left mouse button.

      ; Reset button down state.
      (*pState).bButtonDown = 0b
      (*pState).bTempSegment = 0b

      ; End ROI
      oROI = (*pState).oCurrROI
      if (not OBJ_VALID(oROI)) then break

      ; Ensure that the region has at 3 vertices.
      oROI->GetProperty, DATA=roiData
      if ((N_ELEMENTS(roiData)/3) ge 3) then begin

        ; The region is valid.  Give it a name and add
        ; it to the appropriate containers.
        oROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oModel->Remove, oROI
        (*pState).oROIModel->Add, oROI
        (*pState).oROIGroup->Add, oROI
        if OBJ_VALID((*pState).oRegionsOut) then begin
          (*pState).oRegionsOut->Add, oROI
        end
        (*pState).oCurrROI = OBJ_NEW()

        ;                ; Activate appropriate tool buttons.
        ;                if lmgr(/demo) ne 1 then begin
        ;                    WIDGET_CONTROL, (*pState).wSaveButton, $
        ;                        SENSITIVE=1
        ;                    WIDGET_CONTROL, $
        ;                        (*pState).wSaveToolButton, $
        ;                        SENSITIVE=1
        ;                endif

        ; Set the region as current.
        xd_xroi__SetROI, pState, oROI, /UPDATE_LIST, $
          /SET_LIST_SELECT

        oROI->SetProperty, STYLE=2

        modify_mask_by_ROI, pState
        (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac

        ; If this is the first region, bring up the
        ; region information dialog.
        if ((*pState).bFirstROI eq 1b) then begin
          WIDGET_CONTROL, /HOURGLASS
          (*pState).bFirstROI = 0b
          xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
        endif

      endif else begin
        ; Fewer than 3 vertices; delete.
        (*pState).oModel->Remove, oROI
        OBJ_DESTROY, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ; Reset color of formerly selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, $
          COLOR=(*pState).sel_rgb

        (*pState).oWindow->Draw, (*pState).oView

      endelse
    end  ; FREEHAND DRAW


    ; Circular Mask ROI
    'CIRCULARMASK DRAW': begin
      if (sEvent.release ne 1) then break
      if ((*pState).bButtonDown ne 1) then break  ; Left mouse button.

      ; Reset button down state.
      (*pState).bButtonDown = 0b
      (*pState).bTempSegment = 0b

      ; End ROI
      oROI = (*pState).oCurrROI
      if (not OBJ_VALID(oROI)) then break

      ; Ensure that the region has at 3 vertices.
      oROI->GetProperty, DATA=roiData
      if ((N_ELEMENTS(roiData)/3) ge 3) then begin

        ; The region is valid.  Give it a name and add
        ; it to the appropriate containers.
        oROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oModel->Remove, oROI
        (*pState).oROIModel->Add, oROI
        (*pState).oROIGroup->Add, oROI
        if OBJ_VALID((*pState).oRegionsOut) then begin
          (*pState).oRegionsOut->Add, oROI
        end
        (*pState).oCurrROI = OBJ_NEW()

        ;by comparing the image difference to extract ROI
        ;            IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN extract_mask_diff_to_ROI, pState
        extract_mask_diff_to_ROI, pState
        (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac
        ;            print, 'alpha channel = ', aac

        IF STRCMP((*pState).direction, 'xy') THEN Location = (*pState).location3D[2]
        IF STRCMP((*pState).direction, 'yz') THEN Location = (*pState).location3D[0]
        IF STRCMP((*pState).direction, 'xz') THEN Location = (*pState).location3D[1]


        dims = SIZE((*pState).new_mask_vol_cube, /DIMENSIONS)
        (*pState).oImage2->getProperty, DATA=update_mask
        update_mask = CONGRID(update_mask, dims[0], dims[1])
        i=Location
        IF STRCMP((*pState).direction, 'xy') THEN BEGIN
          (*pState).new_mask_vol_cube[*,*,i*dims[2]/(*pState).imageSize] = REVERSE(update_mask, 2)
        ENDIF

        IF STRCMP((*pState).direction, 'yz') THEN BEGIN
          (*pState).new_mask_vol_cube[i*dims[0]/(*pState).imageSize,*,*] = update_mask
        ENDIF

        IF STRCMP((*pState).direction, 'xz') THEN BEGIN
          (*pState).new_mask_vol_cube[*,i*dims[1]/(*pState).imageSize,*] = update_mask
        ENDIF


        ; If this is the first region, bring up the
        ; region information dialog.
        if ((*pState).bFirstROI eq 1b) then begin
          WIDGET_CONTROL, /HOURGLASS
          (*pState).bFirstROI = 0b
          xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
        endif

      endif else begin
        ; Fewer than 3 vertices; delete.
        (*pState).oModel->Remove, oROI
        OBJ_DESTROY, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ; Reset color of formerly selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, $
          COLOR=(*pState).sel_rgb

        (*pState).oWindow->Draw, (*pState).oView

      endelse

      ;          (*pState).oImage1->getProperty, DATA=Mask_new
      ;          iimage, Mask_new

      ;after processing, replot circular ROI
      oROI = OBJ_NEW('IDLgrROI', $
        COLOR=(*pState).sel_rgb, $
        STYLE=1 $
        )

      style_closed = 2
      Radii = (*pState).Radii
      nPts = Radii * 4
      a = FINDGEN(nPts) * ( (2 * !PI) / (nPts-1) )
      newX = COS(a) * Radii + xImage
      newY = SIN(a) * Radii + yImage
      newZ = REPLICATE(0.0, nPts)
      style = style_closed

      oROI->GetProperty, N_VERTS=nVerts
      oROI->ReplaceData, newX, newY, newZ, START=0, FINISH=nVerts-1
      oROI->SetProperty, STYLE=style
      (*pState).oCurrROI = oROI
      (*pState).oModel->Add, oROI
    end  ; CIRCULARMASK DRAW


    ; Segmented ROI
    'POLYGON DRAW': begin
      ; Double-click or right mouse-button up.
      if not ((sEvent.release eq 1 and (*pState).bButtonDown eq 2) or $
        (sEvent.release eq 4)) then break

      ; Reset button down state.
      (*pState).bButtonDown = 0b

      oROI = (*pState).oCurrROI
      if (not OBJ_VALID(oROI)) then break

      value = $
        'x:' + STRING(xImage, FORMAT='(i5)') + '  ' + $
        'y:' + STRING(yImage, FORMAT='(i5)')
      if (*pState).image_is_8bit then begin
        value = value + '  z:' + STRING( $
          (*(*pState).pImg)[xImage, yImage], $
          FORMAT='(i4)' $
          )
      endif

      WIDGET_CONTROL, (*pState).wStatus, SET_VALUE=value

      ; Ensure that the region has at 3 vertices.
      oROI->GetProperty, DATA=roiData
      if ((N_ELEMENTS(roiData)/3) ge 3) then begin
        ; The region is valid.  Give it a name and add it
        ; to the appropriate containers.
        oROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oModel->Remove, oROI
        (*pState).oROIModel->Add, oROI
        (*pState).oROIGroup->Add, oROI
        if OBJ_VALID((*pState).oRegionsOut) then begin
          (*pState).oRegionsOut->Add, oROI
        end
        (*pState).oCurrROI = OBJ_NEW()

        ;                ; Activate appropriate tool buttons.
        ;                if not lmgr(/demo) then begin
        ;                    WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=1
        ;                    WIDGET_CONTROL, (*pState).wSaveToolButton, $
        ;                        SENSITIVE=1
        ;                endif

        ; Set the region as current.
        xd_xroi__SetROI, pState, oROI, /UPDATE_LIST, $
          /SET_LIST_SELECT

        oROI->SetProperty, STYLE=2

        modify_mask_by_ROI, pState
        (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac

        ; If this is the first region, bring up the
        ; region information dialog.
        if ((*pState).bFirstROI eq 1b) then begin
          WIDGET_CONTROL, /HOURGLASS
          (*pState).bFirstROI = 0b
          xd_xroiInfo, pState, GROUP_LEADER=sEvent.top
        endif

        ; Reset state.
        (*pState).bTempSegment = 0b

      endif else begin
        ; Fewer than 3 vertices; delete.
        (*pState).oModel->Remove, oROI
        OBJ_DESTROY, oROI
        (*pState).oCurrROI = OBJ_NEW()

        ; Reset color of formerly selected ROI.
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, COLOR=(*pState).sel_rgb
        (*pState).oWindow->Draw, (*pState).oView

        ; Reset state.
        (*pState).bTempSegment = 0b
      endelse

      ; This is a special case. For a segmented polygon, we don't want
      ; to check for the right mouse button below because this will
      ; pop up the context menu. So just bail out here.
      return

    end   ; POLYGON DRAW

    ; Pick ROI.
    'SELECTION': begin
      (*pState).bButtonDown = 0b
    end
  endcase


  ; On a right mouse up, perform a selection, then display the
  ; appropriate context menu.
  if (sEvent.release eq 4) then begin
    ; Temporarily hide the image for selection.
    (*pState).oImage->SetProperty, HIDE=1

    ; Temporarily hide the selection visual.
    if (OBJ_VALID((*pState).oSelVisual) ne 0) then begin
      (*pState).oSelVisual->GetProperty, HIDE=oldHide
      (*pState).oSelVisual->SetProperty, /HIDE
    endif

    ; Check if a region was hit.
    oSel = (*pState).oWindow->Select((*pState).oView, $
      DIMENSIONS=[16,16], $
      [sEvent.x, sEvent.y])

    ; Restore the selection visual.
    if (OBJ_VALID((*pState).oSelVisual) ne 0) then $
      (*pState).oSelVisual->SetProperty, HIDE=oldHide

    selType = SIZE(oSel, /TYPE)
    ;        if (selType eq 11) then begin  ; Object reference?
    ;            oROI = oSel[0]
    ;
    ;            ; Mark the region as selected.
    ;            xd_xroi__SetROI, pState, oROI, /SET_LIST_SELECT
    ;
    ;
    ;            oROI->SetProperty, STYLE=2
    ;
    ;            modify_mask_by_ROI, pState
    ;                   (*pState).oSelROI->getProperty, ALPHA_CHANNEL=aac
    ;                   (*pState).oSelROI->setProperty, ALPHA_CHANNEL=1.0
    ;                   (*pState).oWindow->Draw, (*pState).oView
    ;                   (*pState).oSelROI->setProperty, ALPHA_CHANNEL=aac
    ;
    ;        endif else $
    ;            xd_xroi__SetROI, pState, OBJ_NEW()


    ; Restore the image.
    (*pState).oImage->SetProperty, HIDE=0

    (*pState).oWindow->Draw, (*pState).oView

    ;        if (OBJ_VALID((*pState).oSelROI) ne 0) then $
    ;            WIDGET_DISPLAYCONTEXTMENU, sEvent.id, sEvent.x, sEvent.y, $
    ;                (*pState).wROIContextMenu

  endif

end


;------------------------------------------------------------------------------
pro xd_xroi__Motion, sEvent
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  ON_ERROR, (*pState).debug ? 0 : 2

  ; Convert from viewport coordinates to image coordinates.
  WIDGET_CONTROL, (*pState).wDraw, GET_DRAW_VIEW=viewport
  xImage = sEvent.x + viewport[0]
  yImage = sEvent.y + viewport[1]

  (*pState).oImage->GetProperty, DIMENSIONS=dimensions
  event_is_in_image = $
    xImage ge 0 and $
    yImage ge 0 and $
    xImage lt dimensions[0] and $
    yImage lt dimensions[1]

  value = $
    'x:' + STRING(xImage, FORMAT='(i5)') + '  ' + $
    'y:' + STRING(yImage, FORMAT='(i5)')
  if (*pState).image_is_8bit then begin
    if event_is_in_image then begin
      value = value + '  z:' + STRING( $
        (*(*pState).pImg)[xImage, yImage], $
        FORMAT='(i6)' $
        )
    endif else begin
      ;          print, 'event is not in image!'
    endelse
  endif

  if event_is_in_image then begin
    append = ' '
    if OBJ_VALID((*pState).oCurrROI) then begin
      if (*pState).mode eq 'POLYGON DRAW' then begin
        append = '  Double-click to finish'
      endif
    endif else begin
      if (*pState).oROIGroup->Count() gt 0 then begin
        c = (*pState).oROIGroup->ContainsPoints(xImage, yImage)
        case c[0] of
          0: append = '  (Outside)'
          1: append = '  (Inside)'
          2: append = '  (On Edge)'
          3: append = '  (On Vertex)'
        endcase
      endif
    endelse
    value = value + append
  endif

  WIDGET_CONTROL, (*pState).wStatus, SET_VALUE=value

  case (*pState).mode of
    'TRANSLATE-SCALE': begin
      if ((*pState).bButtonDown ne 0) then begin
        oROI = (*pState).oSelROI

        ; If this is the first motion event since the button press,
        ; then temporarily disable the region statistics (until
        ; the button release.
        if ((*pState).bButtonDown eq 1) then begin
          ; Determine if ROI Info dialog is currently realized.
          haveInfoBase = WIDGET_INFO((*pState).wROIInfo, /VALID_ID)
          if (haveInfoBase NE 0) then begin
            WIDGET_CONTROL, (*pState).wROIInfo, GET_UVALUE=sState

            xd_xroiInfo__SetStats, sState, OBJ_NEW()
          endif

          (*pState).bButtonDown = 2
        endif

        ; First check if the mouse down occurred within a
        ; selection visual handle.
        oSelHandle = (*pState).oSelHandle
        if (OBJ_VALID(oSelHandle) ne 0) then begin
          oSelHandle->GetProperty, NAME=handleName

          if (handleName eq 'TRANSLATE') then begin
            ; Mouse down occurred within translation box.
            xd_xroi__TranslateROI, pState, oROI, sEvent.x, sEvent.y
          endif else begin
            ; Mouse down must have occurred in a scale handle.
            xd_xroi__ScaleROI, pState, oROI, sEvent.x, sEvent.y
          endelse

        endif else begin
          ; Translate currently selected region, if any.
          xd_xroi__TranslateROI, pState, oROI, sEvent.x, sEvent.y
        endelse

        (*pState).oWindow->Draw, (*pState).oView

        ; Store new button location.
        (*pState).buttonXY = [sEvent.x, sEvent.y]
      endif
    end

    'RECTANGLE': begin
      oROI = (*pState).oCurrROI
      if (OBJ_VALID(oROI) EQ 0) then return

      ; If button down, reposition rectangle corner.
      if ((*pState).bButtonDown NE 0) then begin

        style_point = 0
        style_line = 1
        style_closed = 2

        x0 = (*pState).buttonXY[0]
        y0 = (*pState).buttonXY[1]
        x1 = xImage
        y1 = yImage
        if (x0 eq x1) then begin
          if (y0 eq y1) then begin
            newBox = [[x0,y0,0.0]]
            style = style_point
          endif else begin
            newBox = [[x0,y0,0.0], [x0,y1,0.0]]
            style = style_line
          endelse
        endif else if (y0 eq y1) then begin
          newBox = [[x0,y0,0.0], [x1,y0,0.0]]
          style = style_line
        endif else begin
          newBox = [[x0,y0,0.0],[x1,y0,0.0],[x1,y1,0.0],[x0,y1,0.0]]
          style = style_closed
        endelse

        oROI->GetProperty, N_VERTS=nVerts
        oROI->ReplaceData, newBox, START=0, FINISH=nVerts-1
        oROI->SetProperty, STYLE=style

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    'ELLIPSE': begin
      oROI = (*pState).oCurrROI
      if (OBJ_VALID(oROI) EQ 0) then return

      ; If button down, reposition radii.
      if ((*pState).bButtonDown NE 0) then begin

        style_point = 0
        style_line = 1
        style_closed = 2

        x0 = (*pState).buttonXY[0]
        y0 = (*pState).buttonXY[1]
        x1 = xImage
        y1 = yImage

        if (x0 eq x1) then begin
          if (y0 eq y1) then begin
            newX = [x0]
            newY = [y0]
            newZ = [0.0]
            style = style_point
          endif else begin
            vertRad = (y1 gt y0) ? (y1-y0) : (y0-y1)
            newX = [x0,x0]
            newY = [y0-vertRad,y0+vertRad]
            newZ = [0.0,0.0]
            style = style_line
          endelse
        endif else if (y0 eq y1) then begin
          horizRad = (x1 gt x0) ? (x1-x0) : (x0-x1)
          newX = [x0-horizRad,x0+horizRad]
          newY = [y0,y0]
          newZ = [0.0,0.0]
          style = style_line
        endif else begin
          horizRad = (x1 gt x0) ? (x1-x0) : (x0-x1)
          vertRad = (y1 gt y0) ? (y1-y0) : (y0-y1)

          ; Number of vertices is dependent upon the greater
          ; of the two radii.
          nPts = (horizRad > vertRad) * 4
          a = FINDGEN(nPts) * ( (2 * !PI) / (nPts-1) )
          newX = COS(a) * horizRad + x0
          newY = SIN(a) * vertRad + y0
          newZ = REPLICATE(0.0, nPts)
          style = style_closed
        endelse

        oROI->GetProperty, N_VERTS=nVerts
        oROI->ReplaceData, newX, newY, newZ, START=0, FINISH=nVerts-1
        oROI->SetProperty, STYLE=style

        (*pState).oWindow->Draw, (*pState).oView
      endif
    end

    ; Freehand ROI
    'FREEHAND DRAW': begin
      oROI = (*pState).oCurrROI
      if (OBJ_VALID(oROI) EQ 0) then return

      ; If button down, append a vertex.
      if ((*pState).bButtonDown NE 0) then begin
        oROI->AppendData, [xImage, yImage]
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).bTempSegment = 1b
      endif
    end

    ; Circular Mask ROI
    'CIRCULARMASK DRAW': begin

      ;check circular ROI size
      oROI = (*pState).oCurrROI
      if (OBJ_VALID(oROI) eq 0) then begin
        oOldSelROI = (*pState).oSelROI
        if (OBJ_VALID(oOldSelROI) ne 0) then $
          oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

        oROI = OBJ_NEW('IDLgrROI', $
          COLOR=(*pState).sel_rgb, $
          STYLE=1 $
          )

        (*pState).oCurrROI = oROI
        (*pState).oModel->Add, oROI
      endif


      oROI = (*pState).oCurrROI
      if (OBJ_VALID(oROI) EQ 0) then return
      style_closed = 2
      Radii = (*pState).Radii
      nPts = Radii * 4
      a = FINDGEN(nPts) * ( (2 * !PI) / (nPts-1) )
      newX = COS(a) * Radii + xImage
      newY = SIN(a) * Radii + yImage
      newZ = REPLICATE(0.0, nPts)
      style = style_closed

      oROI->GetProperty, N_VERTS=nVerts
      oROI->ReplaceData, newX, newY, newZ, START=0, FINISH=nVerts-1
      oROI->SetProperty, STYLE=style


      ; If button down, append a vertex.
      if ((*pState).bButtonDown NE 0) then begin

        ;get the size of oImage
        (*pState).oImage->getProperty, DATA=Img
        roi_mask_dims = SIZE(Img, /DIMENSIONS)
        mask_ROI = oROI -> ComputeMask(DIMENSIONS = roi_mask_dims)
        ;update the mask
        (*pState).oImage1->GetProperty, DATA=mask_old
        IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN mask = mask_old OR mask_ROI
        IF STRCMP((*pState).Modify_ROI, 'ADD') THEN mask = mask_old OR mask_ROI
        IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN mask = ((mask_old GT 100) AND (1-(mask_ROI GT 100)))*255

        (*pState).oImage1->SetProperty, DATA=Mask
        (*pState).oImage2->getProperty, DATA=mask_ROI_old
        (*pState).oImage2->setProperty, DATA=mask_ROI_old+mask_ROI
        (*pState).oWindow->Draw, (*pState).oView
      endif else begin
        (*pState).oWindow->Draw, (*pState).oView
      endelse
    end

    ; Segmented ROI
    'POLYGON DRAW': begin
      oROI = (*pState).oCurrROI
      if (OBJ_VALID(oROI) EQ 0) then return

      ;            print, '(*pState).bTempSegment = ', (*pState).bTempSegment
      ; Replace the final vertex with current mouse location.
      if ((*pState).bTempSegment eq 0) then begin
        oROI->AppendData, [xImage, yImage]
        (*pState).oWindow->Draw, (*pState).oView
        (*pState).bTempSegment = 1b
      endif else begin
        oROI->ReplaceData, [xImage, yImage]
        (*pState).oWindow->Draw, (*pState).oView
      endelse
    end

    ; Pick ROI/Vertex.
    'SELECTION': begin
      if ((*pState).bButtonDown ne 0) then begin
        ; Pick nearest vertex in ROI.
        oROI = (*pState).oSelROI

        if not OBJ_VALID(oROI) then $
          RETURN

        vertIndx = oROI->PickVertex((*pState).oWindow, $
          (*pState).oView, $
          [sEvent.x, sEvent.y])
        if (vertIndx ge 0) then begin
          ; Show the pick vertex.
          oPickPolyline = (*pState).oPickVisual->Get()
          oPickPolyline->SetProperty, HIDE=0

          (*pState).oPickVisual->Reset
          oROI->GetProperty, DATA=vertData
          selVert = vertData[*,vertIndx]
          (*pState).oPickVisual->Translate, $
            selVert[0], selVert[1], selVert[2]

          (*pState).oWindow->Draw, (*pState).oView
        endif
      endif
    end
  endcase

end

;------------------------------------------------------------------------------
pro xd_xroi__Expose, sEvent
  COMPILE_OPT idl2, hidden

  ; Handle expose event in the draw area.

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  ; If we know our draw time is long, then set the hourglass.
  if ((*pState).draw_time gt 0.1) then $
    WIDGET_CONTROL, /HOURGLASS

  ; With APP_SCROLL, the default is RETAIN=1 (system backing), so we only
  ; get expose events when the window is resized. We use these expose
  ; events to time the drawing speed, so we know whether to use a
  ; hourglass or not.
  t = systime(1)
  (*pState).oWindow->Draw, (*pState).oView
  (*pState).draw_time = systime(1) - t

end


;------------------------------------------------------------------------------
pro xd_xroi__Viewport, sEvent

  compile_opt idl2, hidden

  ; Handle viewport move (scroll) event in the draw area.

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  draw_geom = WIDGET_INFO((*pState).wDraw, /GEOM)

  ; On Motif, the geom.xsize may include room for the scrollbar even if
  ; no scrollbar is present. So restrict size to <= to the virtual canvas.
  draw_xsize = draw_geom.xsize < draw_geom.draw_xsize
  draw_ysize = draw_geom.ysize < draw_geom.draw_ysize

  (*pState).oView->SetProperty, $
    VIEWPLANE_RECT=[sEvent.x, sEvent.y, draw_xsize, draw_ysize]

  ; If we know our draw time is long, then set the hourglass.
  if ((*pState).draw_time gt 0.1) then $
    WIDGET_CONTROL, /HOURGLASS

  (*pState).oWindow->Draw, (*pState).oView
end


;------------------------------------------------------------------------------
; Copy the entire image to the system clipboard.
;
function xd_xroi__Copy, sEvent
  COMPILE_OPT idl2, hidden


  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState


  (*pState).oWindow->GetProperty, $
    DIMENSIONS=dims, RESOLUTION=res

  ; Retrieve the image dimensions.
  (*pState).oImage->GetProperty, DIMENSIONS=imagedims

  ; Maximum clipboard size according to MesaGL.
  clipmax = 2048

  ; If too big, clip the image and throw up a warning.
  if ((imagedims[0] gt clipmax) or (imagedims[1] gt clipmax)) then begin
    imagedims = imagedims < clipmax
    dummy = DIALOG_MESSAGE([ $
      'Image dimensions exceeds clipboard maximum of ' + $
      STRTRIM(clipmax, 2) + '.', $
      'Clipboard copy will be truncated.'], $
      DIALOG_PARENT=sEvent.top)
  endif

  WIDGET_CONTROL, /HOURGLASS

  ; Cache the viewplane_rect so we can restore it.
  (*pState).oView->GetProperty, VIEWPLANE_RECT=viewplane_rect

  ; Set the viewplane to include the entire image.
  (*pState).oView->SetProperty, $
    VIEWPLANE_RECT=[0,0,imagedims[0],imagedims[1]]

  oClipboard = OBJ_NEW('IDLgrClipboard', $
    DIMENSIONS=imagedims, $
    RESOLUTION=res)
  oClipboard->Draw, (*pState).oView
  OBJ_DESTROY, oClipboard

  ; Restore the viewplane_rect to its original (possibly cropped) value.
  (*pState).oView->SetProperty, $
    VIEWPLANE_RECT=viewplane_rect

  RETURN, 0 ; "Swallow" event.
end

;------------------------------------------------------------------------------
pro xd_xroiInfo__SetStats, sState, oROI
  COMPILE_OPT idl2, hidden

  ON_ERROR, KEYWORD_SET(sState.debug) ? 0 : 2

  ; Update the statistics in the ROI Info dialog for the given ROI.

  pParentState = sState.pParentState

  if (OBJ_VALID(oROI) eq 0) then begin
    WIDGET_CONTROL, sState.wLabelArea, SET_VALUE="N/A"
    WIDGET_CONTROL, sState.wLabelPerim, SET_VALUE="N/A"
    WIDGET_CONTROL, sState.wLabelNPixel, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelMin, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelMax, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelMean, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelStdDev, SET_VALUE='N/A'
    return
  endif

  ; Compute geometrical statistics.
  result = oROI->ComputeGeometry(AREA=area, PERIMETER=perim)
  if (result ne 0) then begin
    WIDGET_CONTROL, sState.wLabelArea, SET_VALUE=STRTRIM(area, 2)
    WIDGET_CONTROL, sState.wLabelPerim, SET_VALUE=STRTRIM(perim, 2)
  endif else begin
    WIDGET_CONTROL, sState.wLabelArea, SET_VALUE="N/A"
    WIDGET_CONTROL, sState.wLabelPerim, SET_VALUE="N/A"
  endelse

  ; Compute pixel statistics.
  (*pParentState).oImage->GetProperty, DATA=img, DIMENSIONS=dims
  mask = oROI->ComputeMask(DIMENSIONS=dims, MASK_RULE=2)
  if N_ELEMENTS(mask) GT 0 then begin
    if SIZE(img, /N_DIMENSIONS) gt 2 then begin ; 24-bit image.
      IMAGE_STATISTICS, BYTARR(dims[0],dims[1]), MASK=mask, $
        COUNT=pxlCount
      WIDGET_CONTROL, $
        sState.wLabelNPixel, $
        SET_VALUE=STRTRIM(pxlCount,2)
      WIDGET_CONTROL, sState.wLabelMin, SET_VALUE='N/A'
      WIDGET_CONTROL, sState.wLabelMax, SET_VALUE='N/A'
      WIDGET_CONTROL, sState.wLabelMean, SET_VALUE='N/A'
      WIDGET_CONTROL, sState.wLabelStdDev, SET_VALUE='N/A'
    end else begin
      IMAGE_STATISTICS, img, MASK=mask, COUNT=pxlCount, $
        MEAN=pxlMean, STDDEV=pxlStdDev, $
        MINIMUM=pxlMin, MAXIMUM=pxlMax
      WIDGET_CONTROL, sState.wLabelNPixel, SET_VALUE=STRTRIM(pxlCount,2)
      WIDGET_CONTROL, sState.wLabelMin, SET_VALUE=STRTRIM(pxlMin,2)
      WIDGET_CONTROL, sState.wLabelMax, SET_VALUE=STRTRIM(pxlMax,2)
      WIDGET_CONTROL, sState.wLabelMean, SET_VALUE=STRTRIM(pxlMean,2)
      WIDGET_CONTROL, sState.wLabelStdDev, SET_VALUE=STRTRIM(pxlStdDev,2)
    endelse
  endif else begin
    WIDGET_CONTROL, sState.wLabelNPixel, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelMin, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelMax, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelMean, SET_VALUE='N/A'
    WIDGET_CONTROL, sState.wLabelStdDev, SET_VALUE='N/A'
  endelse
end

;------------------------------------------------------------------------------
pro xd_xroi__SetROI, pState, oROI, $
  SET_LIST_SELECT=set_list_select, $ ; IN: Select in list?
  UPDATE_LIST=update_list            ; IN: Add to list?

  COMPILE_OPT idl2, hidden

  ; Set the given ROI as the currently selected ROI.

  ; Determine if ROI Info dialog is currently realized.
  haveInfoBase = WIDGET_INFO((*pState).wROIInfo, /VALID_ID)
  if (haveInfoBase NE 0) then $
    WIDGET_CONTROL, (*pState).wROIInfo, GET_UVALUE=sState

  ; Reset the list of names in the region list.
  if (KEYWORD_SET(update_list) ne 0) then begin
    if (haveInfoBase ne 0) then begin
      oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
      if (nROIs gt 0) then begin
        roiNames = STRARR(nROIs)
        for i=0,nROIs-1 do begin
          oROIs[i]->GetProperty, NAME=name
          roiNames[i] = name
        endfor
      endif else $
        roiNames = ['']

      WIDGET_CONTROL, sState.wList, SET_VALUE=roiNames
    endif
  endif

  oOldSelROI = (*pState).oSelROI
  ; Reset formerly selected ROI's color.
  if (OBJ_VALID(oOldSelROI) ne 0) then $
    oOldSelROI->SetProperty, COLOR=(*pState).roi_rgb

  if (OBJ_VALID(oROI) ne 0) then begin
    ; Set newly selected ROI's color.
    oROI->SetProperty, COLOR=(*pState).sel_rgb

    ; Reshape the current selection visual, if any, to match the
    ; selected region.
    xd_xroi__ReshapeSelectionVisual, pState, oROI

    if (haveInfoBase ne 0) then begin
      ; Update name text widget.
      oROI->GetProperty, NAME=name
      WIDGET_CONTROL, sState.wName, SET_VALUE=name, SENSITIVE=1

      ; Update selected list item.
      if (KEYWORD_SET(set_list_select) NE 0) then begin
        result = (*pState).oROIModel->IsContained(oROI, $
          POSITION=pos)
        if (result ne 0) then $
          WIDGET_CONTROL, sState.wList, SET_LIST_SELECT=pos
      endif

      ; Update delete and histogram button sensitivity.
      WIDGET_CONTROL, sState.wDeleteButton, SENSITIVE=1
      WIDGET_CONTROL, sState.wHistButton, SENSITIVE=1
    endif
  endif else begin
    ; New ROI is not valid; desensitize appropriate widgets.
    if (haveInfoBase ne 0) then begin
      WIDGET_CONTROL, sState.wName, SET_VALUE=' ', SENSITIVE=0
      WIDGET_CONTROL, sState.wDeleteButton, SENSITIVE=0
      WIDGET_CONTROL, sState.wHistButton, SENSITIVE=0
      WIDGET_CONTROL, sState.wList, SET_LIST_SELECT=-1
    endif

    ; Hide any selection visuals.
    if (OBJ_VALID((*pState).oSelVisual) ne 0) then $
      (*pState).oSelVisual->SetProperty, /HIDE
  endelse
  (*pState).oSelROI = oROI

  ; Update the statistics in the ROI Info dialog if appropriate.
  if (haveInfoBase ne 0) then $
    xd_xroiInfo__SetStats, sState, oROI
end
;------------------------------------------------------------------------------
function xd_xroi__CreatePickVisual
  oModel = OBJ_NEW('IDLgrModel', HIDE=1)

  ; Crosshair polyline.
  oPolyline = OBJ_NEW('IDLgrPolyline', [[-8,0],[8,0],[0,-8],[0,8]], $
    POLYLINES=[2,0,1,2,2,3], COLOR=[255,255,0])

  oModel->Add, oPolyline

  return, oModel
end
;------------------------------------------------------------------------------
function xd_xroi__CreateTransScaleVisual, COLOR=color
  oModel = OBJ_NEW('IDLgrModel', /HIDE)

  transpWhiteImage = BYTARR(4,2,2)
  transpWhiteImage[0:2,*,*] = 255
  transpWhiteImage[3,*,*] = 0
  oTranspImage = OBJ_NEW('IDLgrImage', transpWhiteImage, INTERLEAVE=0, $
    /HIDE)
  oModel->Add, oTranspImage ; Add for cleanup

  ; Create scale handles.
  oScaleBox = OBJ_NEW('IDLgrPolygon', [-5,5,5,-5],[-5,-5,5,5], $
    TEXTURE_MAP=oTranspImage, TEXTURE_COORD=[[0,0],[1,0],[1,1],[0,1]], $
    NAME='SCALE_BOX', COLOR=[255,255,255])
  oScaleBoxOutline = OBJ_NEW('IDLgrPolygon', [-5,5,5,-5],[-5,-5,5,5], $
    STYLE=1, NAME='SCALE_BOX_OUTLINE', COLOR=color)

  oScaleLLModel = OBJ_NEW('IDLgrModel', /SELECT_TARGET, NAME='SCALE_LL')
  oScaleLLModel->Add, oScaleBox
  oScaleLLModel->Add, oScaleBoxOutline
  oModel->Add, oScaleLLModel

  oScaleLRModel = OBJ_NEW('IDLgrModel', /SELECT_TARGET, NAME='SCALE_LR')
  oScaleLRModel->Add, oScaleBox, /ALIAS
  oScaleLRModel->Add, oScaleBoxOutline, /ALIAS
  oModel->Add, oScaleLRModel

  oScaleULModel = OBJ_NEW('IDLgrModel', /SELECT_TARGET, NAME='SCALE_UL')
  oScaleULModel->Add, oScaleBox, /ALIAS
  oScaleULModel->Add, oScaleBoxOutline, /ALIAS
  oModel->Add, oScaleULModel

  oScaleURModel = OBJ_NEW('IDLgrModel', /SELECT_TARGET, NAME='SCALE_UR')
  oScaleURModel->Add, oScaleBox, /ALIAS
  oScaleURModel->Add, oScaleBoxOutline, /ALIAS
  oModel->Add, oScaleURModel

  ; Create translation bounding box.
  oTransModel = OBJ_NEW('IDLgrModel', /SELECT_TARGET, NAME='TRANSLATE')
  oTransBoxOutline = OBJ_NEW('IDLgrPolygon', STYLE=1, LINESTYLE=1, $
    NAME='TRANSLATE_BOX_OUTLINE',  COLOR=color)
  oTransModel->Add, oTransBoxOutline
  oModel->Add, oTransModel

  ; Save object references for easy access.
  sState = {$
    oScaleLLModel: oScaleLLModel, $
    oScaleLRModel: oScaleLRModel, $
    oScaleULModel: oScaleULModel, $
    oScaleURModel: oScaleURModel, $
    oScaleBoxOutline: oScaleBoxOutline, $
    oTransModel: oTransModel, $
    oTransBoxOutline: oTransBoxOutline $
  }
  oModel->SetProperty, UVALUE=PTR_NEW(sState,/NO_COPY)

  return, oModel
end
;------------------------------------------------------------------------------
pro xd_xroi__ReshapeSelectionVisual, pState, oROI

  COMPILE_OPT idl2, hidden

  ; If the region is not valid, hide any selection visuals and return.
  if (OBJ_VALID(oROI) eq 0) then begin
    if (OBJ_VALID((*pState).oSelVisual) ne 0) then $
      (*pState).oSelVisual->SetProperty, HIDE=1
    RETURN
  endif

  case (*pState).mode of
    'TRANSLATE-SCALE': begin
      ; Update translate/scale selection visual to match
      ; selected region's bounding box.
      oROI->GetProperty, ROI_XRANGE=xrange, ROI_YRANGE=yrange
      x0 = MIN(xrange, MAX=x1)
      y0 = MIN(yrange, MAX=y1)
      newBox = [[x0,y0,0.0],[x1,y0,0.0],[x1,y1,0.0],[x0,y1,0.0]]

      (*pState).oTransScaleVisual->GetProperty, UVALUE=pTSState

      (*pTSState).oTransModel->Reset
      (*pTSState).oTransBoxOutline->SetProperty, DATA=newBox

      (*pTSState).oScaleLLModel->Reset
      (*pTSState).oScaleLLModel->Translate, x0, y0, 0
      (*pTSState).oScaleLRModel->Reset
      (*pTSState).oScaleLRModel->Translate, x1, y0, 0
      (*pTSState).oScaleULModel->Reset
      (*pTSState).oScaleULModel->Translate, x0, y1, 0
      (*pTSState).oScaleURModel->Reset
      (*pTSState).oScaleURModel->Translate, x1, y1, 0
      (*pTSState).oScaleLLModel->SetProperty, UVALUE=[x0,y0]
      (*pTSState).oScaleURModel->SetProperty, UVALUE=[x1,y1]

      ; Make sure the selection visual is visible.
      (*pState).oTransScaleVisual->SetProperty, HIDE=0
    end

    'SELECTION': begin
      ; Make sure the pick visual is visible.
      (*pState).oPickVisual->SetProperty, HIDE=0

      ; For now, hide vertex within the visual since we do not
      ; know where to position it.
      oPickPolyline = (*pState).oPickVisual->Get()
      oPickPolyline->SetProperty, HIDE=1
    end

    else: begin
    end
  endcase
end
;------------------------------------------------------------------------------
pro xd_xroi__SetSelectionVisualColor, pState, rgb
  ; Set the translation-scale selection visual color.
  (*pState).oTransScaleVisual->GetProperty, UVALUE=pTSState
  (*pTSState).oScaleBoxOutline->SetProperty, COLOR=rgb
  (*pTSState).oTransBoxOutline->SetProperty, COLOR=rgb

  ; Note: currently, the pick selection visual color remains fixed.
end
;------------------------------------------------------------------------------
pro xd_xroi__GrowSelectedROI, pState, STDDEV=stddev

  COMPILE_OPT idl2, hidden



  ; If no currently selected regions, return.
  oROI = (*pState).oSelROI
  if (OBJ_VALID(oROI) eq 0) then $
    RETURN

  ; If image is RGB, either convert to a luminosity image or
  ; use the red, green, or blue channel.
  if ((*pState).image_is_8bit eq 0) then begin
    (*pState).oImage->GetProperty, DATA=rgbData, INTERLEAVE=iInterleave

    allDims = SIZE(rgbData, /DIMENSIONS)
    case iInterleave of
      0: begin
        imgDimensions = allDims[1:2]
        if ((*pState).growRGBUseLuminosity ne 0) then begin
          r = REFORM(rgbData[0, *, *],imgDimensions)
          g = REFORM(rgbData[1, *, *],imgDimensions)
          b = REFORM(rgbData[2, *, *],imgDimensions)
          imgData = BYTE((0.3 * r) + (0.59 * g) + (0.11 * b))
        endif else $
          imgData = REFORM(rgbData[(*pState).growRGBChannel, *, *], $
          imgDimensions)
      end
      1: begin
        imgDimensions = [allDims[0], allDims[2]]
        if ((*pState).growRGBUseLuminosity ne 0) then begin
          r = REFORM(rgbData[*, 0, *],imgDimensions)
          g = REFORM(rgbData[*, 1, *],imgDimensions)
          b = REFORM(rgbData[*, 2, *],imgDimensions)
          imgData = BYTE((0.3 * r) + (0.59 * g) + (0.11 * b))
        endif else $
          imgData = REFORM(rgbData[*, (*pState).growRGBChannel, *], $
          imgDimensions)
      end
      2: begin
        imgDimensions = [allDims[0:1]]
        if ((*pState).growRGBUseLuminosity ne 0) then begin
          r = REFORM(rgbData[*, *, 0],imgDimensions)
          g = REFORM(rgbData[*, *, 1],imgDimensions)
          b = REFORM(rgbData[*, *, 2],imgDimensions)
          imgData = BYTE((0.3 * r) + (0.59 * g) + (0.11 * b))
        endif else $
          imgData = REFORM(rgbData[*, *, (*pState).growRGBChannel], $
          imgDimensions)
      end
    endcase
  endif else begin
    (*pState).oImage->GetProperty, DATA=imgData
    imgDimensions = SIZE(imgData, /DIMENSIONS)
  endelse

  ; Compute a mask for the region.
  roiMask = oROI->ComputeMask(DIMENSIONS=imgDimensions)
  roiPixels = WHERE(roiMask ne 0, count)
  if (count eq 0) then return

  ; Grow the region.
  if (KEYWORD_SET(stddev) ne 0) then begin
    stddev_mult = (*pState).growStdDevMult
    growPixels = REGION_GROW(imgData, roiPixels, $
      STDDEV_MULTIPLIER=stddev_mult, $
      ALL_NEIGHBORS=(*pState).growAllNeighbors)
  endif else begin
    if ((*pState).growThreshROI eq 0) then $
      thresh = [(*pState).growThreshMin, (*pState).growThreshMax]
    growPixels = REGION_GROW(imgData, roiPixels, $
      THRESHOLD=thresh, $
      ALL_NEIGHBORS=(*pState).growAllNeighbors)
  endelse

  if ((N_ELEMENTS(growPixels) eq 1) and (growPixels[0] eq -1)) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      'Grown region contains no pixels.', $
      DIALOG_PARENT=(*pState).wBase)
    RETURN
  endif

  ; Create a mask for the grown region.
  growMask = BYTARR(imgDimensions[0], imgDimensions[1])
  growMask[growPixels] = 255

  ; Contour the grown region's mask.
  CONTOUR, growMask, LEVELS=255, PATH_INFO=pathInfo, PATH_XY=pathXY, $
    /PATH_DATA_COORD
  if (N_ELEMENTS(pathInfo) eq 0) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      ['Grown region could not', $
      'be successfully contoured.'], $
      DIALOG_PARENT=(*pState).wBase)
    RETURN
  endif

  ; Select all exterior contours.
  iids = WHERE(pathInfo[*].high_low eq 1, nContours)
  if (nContours eq  0) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      ['Grown region results in', $
      'no exterior contours.'], $
      DIALOG_PARENT=(*pState).wBase)
    RETURN
  endif

  ; Create a temporary region group to store grown regions.
  oGrowROIGroup = OBJ_NEW('IDLanROIGroup')

  ; Add a new region for each contour.
  iNAccept = 0
  for i=0,nContours-1 do begin
    contId = iids[i]
    iStart = pathInfo[contId].offset
    iFinish = iStart + pathInfo[contId].n - 1

    if (iNAccept gt 0) then begin
      ; Determine if (the first vertex of) this region is contained
      ; within any of the previously accepted regions.
      ; If so, bypass.
      ; Otherwise, accept it as a new valid region.
      containStatus = oGrowROIGroup->ContainsPoints( $
        pathXY[0,iStart], pathXY[0,iStart])
      bAccept = (containStatus[0] eq 0) ? 1b : 0b
    endif else $
      bAccept = 1b

    ; If requested, do not accept if beyond maximum count.
    if (((*pState).growAcceptAll eq 0) and $
      (iNAccept ge (*pState).growMaxCount)) then $
      bAccept = 0b

    if (bAccept ne 0) then begin
      ; Grab the contour path.
      growROIData = FLTARR(3,pathInfo[contId].n)
      growROIData[0,*] = pathXY[0, iStart:iFinish]
      growROIData[1,*] = pathXY[1, iStart:iFinish]

      ; Create a new region.
      oGrowROI = OBJ_NEW('IDLgrROI', growROIData, STYLE=2)

      ; If requested, do not accept if less than minimum area.
      if ((*pState).growAcceptAll eq 0) then begin
        isOk = oGrowROI->ComputeGeometry(AREA=a)
        if (a lt (*pState).growMinArea) then $
          bAccept = 0b
      endif

      ; The contour has passed all acceptance criteria.
      if (bAccept ne 0) then begin
        oGrowROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oROIModel->Add, oGrowROI
        (*pState).oROIGroup->Add, oGrowROI
        if OBJ_VALID((*pState).oRegionsOut) then $
          (*pState).oRegionsOut->Add, oGrowROI

        oGrowROIGroup->Add, oGrowROI

        xd_xroi__SetROI, pState, oGrowROI, /UPDATE_LIST, /SET_LIST_SELECT
        iNAccept = iNAccept + 1
      endif else $
        OBJ_DESTROY, oGrowROI
    endif
  endfor

  if (iNAccept eq 0) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      ['No contours of the grown region', $
      'passed the acceptance criteria.'], $
      DIALOG_PARENT=(*pState).wBase)
    OBJ_DESTROY, oGrowROIGroup
    RETURN
  endif

  ; Select the first grown ROI.
  oGrowROI =  oGrowROIGroup->Get(POSITION=0)
  xd_xroi__SetROI, pState, oGrowROI, /SET_LIST_SELECT

  ; Free temporary region group.
  oGrowROIGroup->Remove, /ALL
  OBJ_DESTROY, oGrowROIGroup

  ;    ; Activate appropriate tool buttons.
  ;    if lmgr(/demo) ne 1 then begin
  ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=1
  ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=1
  ;    endif

end
;------------------------------------------------------------------------------
pro xd_xroi__HistogramSelectedROI, pState, GROUP=group

  ; Nothing to do if no currently selected regions.
  oROI = (*pState).oSelROI
  if (OBJ_VALID(oROI) eq 0) then $
    RETURN

  oROI->GetProperty, NAME=roiName

  ; If image is RGB, prepare to plot all three channels.
  if ((*pState).image_is_8bit eq 0) then begin
    (*pState).oImage->GetProperty, DATA=rgbData, INTERLEAVE=iInterleave
    allDims = SIZE(rgbData, /DIMENSIONS)
    case iInterleave of
      0: imgDimensions = allDims[1:2]
      1: imgDimensions = [allDims[0], allDims[2]]
      2: imgDimensions = [allDims[0:1]]
    endcase
    maxChannel = 2
    rgbNames = ['Red', 'Green', 'Blue']
  endif else begin
    ; If 8-bit image, single plot.
    (*pState).oImage->GetProperty, DATA=data
    imgDimensions = SIZE(data, /DIMENSIONS)
    maxChannel = 0
    rgbNames = ['Black']
  endelse

  ; Determine where the region mask is non-zero.
  w = WHERE(oROI->ComputeMask(DIMENSIONS=imgDimensions))

  havePlot = 0
  for i=0,maxChannel do begin
    if ((*pState).image_is_8bit eq 0) then begin
      case iInterleave of
        0: data = REFORM(rgbData[i, *, *], imgDimensions)
        1: data = REFORM(rgbData[*, i, *], imgDimensions)
        2: data = REFORM(rgbData[*, *, i], imgDimensions)
      endcase
      name = roiName + ' ' + rgbNames[i]
    endif else $
      name = roiName

    if w[0] ne -1 then begin
      case rgbNames[i] of
        'Black': color=[0,0,0]
        'Red': color=[255,0,0]
        'Green': color=[0,255,0]
        'Blue': color=[0,0,255]
        else:
      endcase
      ; Plot the histogram.
      iPlot, $
        HISTOGRAM(data[w]), $
        XTITLE="Bins", $
        YTITLE="Occurrences", $
        NAME=name, $
        COLOR=color, $
        OVERPLOT=i gt 0, $
        /NO_SAVEPROMPT
    endif
  endfor

  oSys = _IDLitSys_GetSystem()
  oTool = oSys->_GetCurrentTool()
  oPlot = (oTool->GetSelectedItems())[0]
  oPlot->GetProperty, _PARENT=oParent
  oParent->Select
  void = oTool->DoAction('OPERATIONS/INSERT/LEGEND')

  ; Add this tool to the list for auto cleanup.
  if (PTR_VALID((*pState).pTools)) then begin
    *(*pState).pTools = [*(*pState).pTools, oTool]
  endif else begin
    (*pState).pTools = PTR_NEW(oTool)
  endelse

end


;------------------------------------------------------------------------------
pro xd_xroi__DeleteSelectedROI, pState

  COMPILE_OPT idl2, hidden

  ; Nothing to do if no regions available.
  nROIs = (*pState).oROIModel->Count()
  if (nROIs eq 0) then return

  ; Nothing to do if no currently selected regions.
  oROI = (*pState).oSelROI
  if (OBJ_VALID(oROI) eq 0) then return

  ; Get position of selected region.
  result =(*pState).oROIModel->IsContained(oROI, POSITION=pos)
  if (result eq 0) then return

  ; Remove and, if appropriate, destroy selected region.
  (*pState).oROIModel->Remove, oROI
  (*pState).oROIGroup->Remove, oROI
  if OBJ_VALID((*pState).oRegionsOut) then $
    (*pState).oRegionsOut->Remove, oROI

  if OBJ_VALID((*pState).oRegionsIn) then begin
    if (*pState).oRegionsIn->IsContained(oROI) then begin
      if OBJ_VALID((*pState).oRejected) then $
        (*pState).oRejected->Add, oROI
    endif else $
      OBJ_DESTROY, oROI
  endif else $
    OBJ_DESTROY, oROI

  ; Pick a new selected region.
  if (nROIs eq 1) then begin
    ;        WIDGET_CONTROL, (*pState).wSaveButton, SENSITIVE=0
    ;        WIDGET_CONTROL, (*pState).wSaveToolButton, SENSITIVE=0
    xd_xroi__SetROI, pState, OBJ_NEW(), /UPDATE_LIST
  endif else begin
    oROI = (*pState).oROIModel->Get(POSITION=((pos-1) > 0))
    xd_xroi__SetROI, pState, oROI, /SET_LIST_SELECT, /UPDATE_LIST
  endelse
end
;------------------------------------------------------------------------------
pro xd_xroi__ScaleROI, pState, oROI, sEventX, sEventY

  COMPILE_OPT idl2, hidden

  ; If the region is not valid, nothing to do.
  if (OBJ_VALID(oROI) eq 0) then $
    RETURN

  ; Retrieve original bounding box corners.
  (*pState).oTransScaleVisual->GetProperty, UVALUE=pTSState
  (*pTSState).oScaleLLModel->GetProperty, UVALUE=xyLL
  (*pTSState).oScaleURModel->GetProperty, UVALUE=xyUR
  x0 = xyLL[0]
  x1 = xyUR[0]
  y0 = xyLL[1]
  y1 = xyUR[1]

  ; Compute deltas relative to previous mouse location.
  dx = sEventX - (*pState).buttonXY[0]
  dy = sEventY - (*pState).buttonXY[1]

  ; Compute new bounding box corners.
  (*pState).oSelHandle->GetProperty, NAME=handleName
  case handleName of
    'SCALE_LL': begin
      newx0 = x0 + dx
      newy0 = y0 + dy

      bSwapX = (newx0 gt x1) ? 1b : 0b
      bSwapY = (newy0 gt y1) ? 1b : 0b

      if (bSwapX) then begin
        newx0 = x1
        newx1 = newx0
        if (bSwapY) then $
          (*pState).oSelHandle = (*pTSState).oScaleURModel $
        else $
          (*pState).oSelHandle = (*pTSState).oScaleLRModel
      endif else $
        newx1 = x1

      if (bSwapY) then begin
        newy0 = y1
        newy1 = newy0
        if (not bSwapX) then $
          (*pState).oSelHandle = (*pTSState).oScaleULModel
      endif else $
        newy1 = y1
    end

    'SCALE_LR': begin
      newx1 = x1 + dx
      newy0 = y0 + dy

      bSwapX = (x0 gt newx1) ? 1b : 0b
      bSwapY = (newy0 gt y1) ? 1b : 0b

      if (bSwapX) then begin
        newx0 = newx1
        newx1 = x0
        if (bSwapY) then $
          (*pState).oSelHandle = (*pTSState).oScaleULModel $
        else $
          (*pState).oSelHandle = (*pTSState).oScaleLLModel
      endif else $
        newx0 = x0

      if (bSwapY) then begin
        newy1 = newy0
        newy0 = y1
        if (not bSwapX) then $
          (*pState).oSelHandle = (*pTSState).oScaleURModel
      endif else $
        newy1 = y1
    end

    'SCALE_UL': begin
      newx0 = x0 + dx
      newy1 = y1 + dy

      bSwapX = (newx0 gt x1) ? 1b : 0b
      bSwapY = (y0 gt newy1) ? 1b : 0b

      if (bSwapX) then begin
        newx1 = newx0
        newx0 = x1
        if (bSwapY) then $
          (*pState).oSelHandle = (*pTSState).oScaleLRModel $
        else $
          (*pState).oSelHandle = (*pTSState).oScaleURModel
      endif else $
        newx1 = x1

      if (bSwapY) then begin
        newy0 = newy1
        newy1 = y0
        if (not bSwapX) then $
          (*pState).oSelHandle = (*pTSState).oScaleLLModel
      endif else $
        newy0 = y0
    end

    'SCALE_UR': begin
      newx1 = x1 + dx
      newy1 = y1 + dy

      bSwapX = (x0 gt newx1) ? 1b : 0b
      bSwapY = (y0 gt newy1) ? 1b : 0b

      if (bSwapX) then begin
        newx0 = newx1
        newx1 = x0
        if (bSwapY) then $
          (*pState).oSelHandle = (*pTSState).oScaleLLModel $
        else $
          (*pState).oSelHandle = (*pTSState).oScaleULModel
      endif else $
        newx0 = x0

      if (bSwapY) then begin
        newy0 = newy1
        newy1 = y0
        if (not bSwapX) then $
          (*pState).oSelHandle = (*pTSState).oScaleLRModel
      endif else $
        newy0 = y0
    end

    else: begin
      newx0 = x0
      newx1 = x1
      newy0 = y0
      newy1 = y1
      bSwapX = 0b
      bSwapY = 0b
    end
  endcase

  ; Compute scale factors relative to original ranges.
  origXRange = (*pState).savedROIXRange
  origYRange = (*pState).savedROIYRange
  sx = float(newx1-newx0+1) / float(origXRange[1]-origXRange[0]+1)
  sy = float(newy1-newy0+1) / float(origYRange[1]-origYRange[0]+1)

  ; Swap the original data as needed.
  newROIData = *(*pState).pSavedROIData
  if (bSwapX) then begin
    newROIData[0,*] = (origXRange[1]-origXRange[0]) - newROIData[0,*]
  endif
  if (bSwapY) then $
    newROIData[1,*] = (origYRange[1]-origYRange[0]) - newROIData[1,*]
  if (bSwapX or bSwapY) then $
    *(*pState).pSavedROIData = newROIData

  ; Scale the original data, and store result.
  newROIData[0,*] = newROIData[0,*] * sx + newx0
  newROIData[1,*] = newROIData[1,*] * sy + newy0
  oROI->ReplaceData, newROIData

  ; Translate the scale selection visual handles.
  (*pTSState).oScaleLLModel->Reset
  (*pTSState).oScaleLLModel->Translate, newx0, newy0, 0
  (*pTSState).oScaleLRModel->Reset
  (*pTSState).oScaleLRModel->Translate, newx1, newy0, 0
  (*pTSState).oScaleULModel->Reset
  (*pTSState).oScaleULModel->Translate, newx0, newy1, 0
  (*pTSState).oScaleURModel->Reset
  (*pTSState).oScaleURModel->Translate, newx1, newy1, 0

  ; Keep track of the new bounding box corners.
  (*pTSState).oScaleLLModel->SetProperty, UVALUE=[newx0,newy0]
  (*pTSState).oScaleURModel->SetProperty, UVALUE=[newx1,newy1]

  ; Reset the translation box (selection visual).
  newBox = [[newx0,newy0],[newx1,newy0],[newx1,newy1],[newx0,newy1]]
  (*pTSState).oTransModel->Reset
  (*pTSState).oTransBoxOutline->SetProperty, DATA=newBox
end
;------------------------------------------------------------------------------
pro xd_xroi__TranslateROI, pState, oROI, sEventX, sEventY

  COMPILE_OPT idl2, hidden

  ; If the region is not valid, nothing to do.
  if (OBJ_VALID(oROI) eq 0) then $
    RETURN

  ; Compute deltas relative to previous mouse location.
  dx = sEventX - (*pState).buttonXY[0]
  dy = sEventY - (*pState).buttonXY[1]

  ; Translate the region.
  oROI->Translate, dx, dy

  ; Move translate/scale selection visual along with the region.
  (*pState).oTransScaleVisual->GetProperty, UVALUE=pTSState
  (*pTSState).oTransModel->Translate, dx, dy, 0

  (*pTSState).oScaleLLModel->Translate, dx, dy, 0
  (*pTSState).oScaleLRModel->Translate, dx, dy, 0
  (*pTSState).oScaleULModel->Translate, dx, dy, 0
  (*pTSState).oScaleURModel->Translate, dx, dy, 0

  ; Update bounding box corners.
  (*pTSState).oScaleLLModel->GetProperty, UVALUE=oldXY
  (*pTSState).oScaleLLModel->SetProperty, UVALUE=oldXY+[dx,dy]
  (*pTSState).oScaleURModel->GetProperty, UVALUE=oldXY
  (*pTSState).oScaleURModel->SetProperty, UVALUE=oldXY+[dx,dy]
end
;------------------------------------------------------------------------------
pro xd_xroiHistPlot_event, event
end
;------------------------------------------------------------------------------
pro xd_xroiInfo_event, sEvent
  COMPILE_OPT idl2, hidden

  ; Handle events in the ROI Info dialogs.

  ; Handle any kill events.
  if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST') $
    then begin
    WIDGET_CONTROL, sEvent.top, /DESTROY
    RETURN
  endif

  WIDGET_CONTROL, sEvent.id, GET_UVALUE=uval
  case uval of
    'ROI_LIST': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState

      nROIs = (*pParentState).oROIModel->Count()
      if (nROIs eq 0) then begin
        xd_xroi__SetROI, pParentState, OBJ_NEW()
        return
      endif

      oROI = (*pParentState).oROIModel->Get(POSITION=sEvent.index)
      xd_xroi__SetROI, pParentState, oROI

      (*pParentState).oWindow->Draw, (*pParentState).oView

    endcase
    'NAME': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      oSelROI = (*pParentState).oSelROI
      if (OBJ_VALID(oSelROI) ne 0) then begin
        WIDGET_CONTROL, sEvent.id, GET_VALUE=name
        oSelROI->SetProperty, NAME=name[0]

        ; Update list widget names.
        oROIs = (*pParentState).oROIModel->Get(/ALL, COUNT=nROIs)
        if (nROIs gt 0) then begin
          roiNames = STRARR(nROIs)
          for i=0,nROIs-1 do begin
            oROIs[i]->GetProperty, NAME=name
            roiNames[i] = name
          endfor
        endif else $
          roiNames = ['']

        WIDGET_CONTROL, sState.wList, SET_VALUE=roiNames

        result = (*pParentState).oROIModel->IsContained(oSelROI,  $
          POSITION=pos)
        if (result ne 0) then $
          WIDGET_CONTROL, sState.wList, SET_LIST_SELECT=pos
      endif
    endcase
    'DELETE': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState

      xd_xroi__DeleteSelectedROI, pParentState

      (*pParentState).oWindow->Draw, (*pParentState).oView
    endcase
    'HISTOGRAM': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState

      xd_xroi__HistogramSelectedROI, pParentState, GROUP=sEvent.top


    endcase
    'CLOSE': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      WIDGET_CONTROL, sEvent.top, /DESTROY
    end
    else:
  endcase
end

;------------------------------------------------------------------------------
pro xd_xroiInfo, pParentState, GROUP_LEADER=group

  COMPILE_OPT idl2, hidden

  ; Open the ROI Info dialog.

  ; Check if already created.  If so, return.
  if (WIDGET_INFO((*pParentState).wROIInfo, /VALID_ID) NE 0) then begin
    WIDGET_CONTROL,(*pParentState).wROIInfo, /SHOW
    return
  endif

  ; Get current list of ROI names.
  oROIs = (*pParentState).oROIModel->Get(/ALL, COUNT=nROIs)
  if (nROIs gt 0) then begin
    roiNames = STRARR(nROIs)
    for i=0,nROIs-1 do begin
      oROIs[i]->GetProperty, NAME=name
      roiNames[i] = name
    endfor
  endif else $
    roiNames = ['']

  ; Determine which ROI, if any, is currently selected.
  oSelROI = (*pParentState).oSelROI
  pos = -1L
  if (OBJ_VALID(oSelROI) ne 0) then $
    result = (*pParentState).oROIModel->IsContained(oSelROI, POSITION=pos)

  title = N_ELEMENTS(*(*pParentState).pTitle) eq 0 ? $
    'ROI' : *(*pParentState).pTitle
  title = title + ' Information'

  prefix = "xd_xroiInfo:"
  wBase = WIDGET_BASE( $
    TITLE=title, $
    GROUP_LEADER=group, $
    /COLUMN, $
    /TLB_KILL_REQUEST_EVENTS, $
    FLOATING=(*pParentState).floating, $
    UNAME=prefix + 'tlb', $
    MODAL=(*pParentState).modal $
    )

  wLabel = WIDGET_LABEL(wBase, VALUE='Regions of Interest:', /ALIGN_LEFT)
  wRowBase = WIDGET_BASE(wBase, /ROW)
  wList = WIDGET_LIST(wRowBase, XSIZE=10, YSIZE=5, VALUE=roiNames, $
    UVALUE='ROI_LIST', UNAME=prefix + 'list')
  wLabelBase = WIDGET_BASE(wRowBase, /COL)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='Area:', /ALIGN_LEFT)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='Perimeter:', /ALIGN_LEFT)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='# Pixels:', /ALIGN_LEFT)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='Minimum:', /ALIGN_LEFT)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='Maximum:', /ALIGN_LEFT)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='Mean:', /ALIGN_LEFT)
  wLabel = WIDGET_LABEL(wLabelBase, VALUE='Std. Dev.:', /ALIGN_LEFT)
  wStatBase = WIDGET_BASE(wRowBase, /COL)
  wLabelArea = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wLabelPerim = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wLabelNPixel = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wLabelMin = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wLabelMax = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wLabelMean = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wLabelStdDev = WIDGET_LABEL(wStatBase, VALUE='N/A       ', /ALIGN_LEFT, $
    /DYNAMIC_RESIZE)
  wRowBase = WIDGET_BASE(wBase, /ROW)
  wLabel = WIDGET_LABEL(wRowBase, VALUE='Name:')
  if (pos ge 0) then $
    currName = roiNames[pos] $
  else $
    currName = ' '
  wName = WIDGET_TEXT(wRowBase, VALUE=currName, XSIZE=20, YSIZE=1, $
    /EDITABLE, SENSITIVE=(pos ge 0), UVALUE='NAME', $
    UNAME=prefix + 'name_text')
  wDeleteButton = WIDGET_BUTTON(wRowBase, VALUE='Delete ROI', $
    UNAME=prefix + 'delete', $
    UVALUE='DELETE', SENSITIVE=(pos ge 0))
  wRowBase = WIDGET_BASE(wBase, /ROW)
  wButton = WIDGET_BUTTON(wRowBase, VALUE='Close', UVALUE='CLOSE', $
    UNAME=prefix + 'close')
  wHistButton = WIDGET_BUTTON( $
    wRowBase, $
    VALUE='Histogram', $
    UVALUE='HISTOGRAM', $
    UNAME=prefix + 'histogram', $
    SENSITIVE=(pos ge 0) $
    )

  (*pParentState).wROIInfo = wBase

  sState = {pParentState: pParentState, $
    wList: wList, $
    wName: wName, $
    wDeleteButton: wDeleteButton, $
    wHistButton: wHistButton, $
    wLabelArea: wLabelArea, $
    wLabelPerim: wLabelPerim, $
    wLabelNPixel: wLabelNPixel, $
    wLabelMin: wLabelMin, $
    wLabelMax: wLabelMax, $
    wLabelMean: wLabelMean, $
    wLabelStdDev: wLabelStdDev, $
    debug: (*pParentState).debug $
  }
  WIDGET_CONTROL, wBase, SET_UVALUE=sState

  if not (*pParentState).modal then $
    WIDGET_CONTROL, wBase, MAP=0

  WIDGET_CONTROL, wBase, /REALIZE

  if not (*pParentState).modal then begin
    wbase_geom = WIDGET_INFO(wBase, /GEOMETRY)
    leader_geom = WIDGET_INFO(group, /GEOMETRY)
    DEVICE, GET_SCREEN_SIZE=screen_size

    x = leader_geom.scr_xsize[0] + leader_geom.xoffset
    x = x < (screen_size[0] - wBase_geom.scr_xsize)
    x = x > 0

    y = leader_geom.yoffset
    y = y < (screen_size[1] - wBase_geom.scr_xsize)
    y = y > 0

    WIDGET_CONTROL, wBase, TLB_SET_XOFFSET=x, TLB_SET_YOFFSET=y
    WIDGET_CONTROL, wBase, MAP=1
  endif

  if (pos ge 0) then begin
    WIDGET_CONTROL, wList, SET_LIST_SELECT=pos
    xd_xroiInfo__SetStats, sState, oSelROI
  endif

  XMANAGER, "xd_xroiInfo", wBase, /NO_BLOCK
end

;------------------------------------------------------------------------------
pro xd_xroiGrowProps_event, sEvent
  COMPILE_OPT idl2, hidden

  ; Handle events in the Region Grow property dialog.

  ; Handle any kill events.
  if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST') $
    then begin
    WIDGET_CONTROL, sEvent.top, /DESTROY
    RETURN
  endif

  WIDGET_CONTROL, sEvent.id, GET_UVALUE=uval
  case uval of
    'SEARCH_ALL_NEIGHBORS': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growAllNeighbors = sEvent.select
    end
    'SEARCH_IMMEDIATE_NEIGHBORS': begin
      ; Toggle handled by 'SEARCH_ALL_NEIGHBORS'.
    end
    'THRESH_MIN': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growThreshMin = sEvent.value
    end
    'THRESH_MAX': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growThreshMax = sEvent.value
    end
    'THRESH_ROI': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growThreshROI = sEvent.select

      if (sEvent.select ne 0) then begin
        ; Make the appropriate fields insensitive.
        WIDGET_CONTROL, sState.wThreshMinField, SENSITIVE=0
        WIDGET_CONTROL, sState.wThreshMaxField, SENSITIVE=0
      endif else begin
        ; Make the appropriate fields sensitive.
        WIDGET_CONTROL, sState.wThreshMinField, SENSITIVE=1
        WIDGET_CONTROL, sState.wThreshMaxField, SENSITIVE=1
      endelse
    end
    'STDDEV_MULT': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growStdDevMult = sEvent.value
    end
    'USE_LUMINOSITY': begin
      if (sEvent.select ne 0) then begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
        pParentState = sState.pParentState
        (*pParentState).growRGBUseLuminosity = 1
      endif
    end
    'USE_RED': begin
      if (sEvent.select ne 0) then begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
        pParentState = sState.pParentState
        (*pParentState).growRGBUseLuminosity = 0
        (*pParentState).growRGBChannel = 0
      endif
    end
    'USE_GREEN': begin
      if (sEvent.select ne 0) then begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
        pParentState = sState.pParentState
        (*pParentState).growRGBUseLuminosity = 0
        (*pParentState).growRGBChannel = 1
      endif
    end
    'USE_BLUE': begin
      if (sEvent.select ne 0) then begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
        pParentState = sState.pParentState
        (*pParentState).growRGBUseLuminosity = 0
        (*pParentState).growRGBChannel = 2
      endif
    end
    'ACCEPT_ALL': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growAcceptAll = sEvent.select

      if (sEvent.select ne 0) then begin
        ; Make the appropriate fields insensitive.
        WIDGET_CONTROL, sState.wMaxCountField, SENSITIVE=0
        WIDGET_CONTROL, sState.wMinAreaField, SENSITIVE=0
        WIDGET_CONTROL, sState.wAreaUnitLabel, SENSITIVE=0
      endif else begin
        ; Make the appropriate fields sensitive.
        WIDGET_CONTROL, sState.wMaxCountField, SENSITIVE=1
        WIDGET_CONTROL, sState.wMinAreaField, SENSITIVE=1
        WIDGET_CONTROL, sState.wAreaUnitLabel, SENSITIVE=1
      endelse
    end
    'MAX_COUNT': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growMaxCount = sEvent.value
    end
    'MIN_AREA': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      pParentState = sState.pParentState
      (*pParentState).growMinArea = sEvent.value
    end
    'CLOSE': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState
      WIDGET_CONTROL, sEvent.top, /DESTROY
    end
  endcase
end
;------------------------------------------------------------------------------
pro xd_xroiGrowProps, pParentState, GROUP_LEADER=group

  COMPILE_OPT idl2, hidden

  ; Open the Region Grow property dialog.

  ; Check if already created.  If so, return.
  if (WIDGET_INFO((*pParentState).wROIGrowProps, /VALID_ID) NE 0) then begin
    WIDGET_CONTROL,(*pParentState).wROIGrowProps, /SHOW
    return
  endif

  prefix = "xd_xroiGrow:"
  title = 'Region Grow Properties'
  wBase = WIDGET_BASE( $
    TITLE=title, $
    GROUP_LEADER=group, $
    /COLUMN, $
    /TLB_KILL_REQUEST_EVENTS, $
    FLOATING=(*pParentState).floating, $
    UNAME=prefix + 'tlb', $
    MODAL=(*pParentState).modal $
    )

  wRowBase = WIDGET_BASE(wBase, /ROW)
  wLabel = WIDGET_LABEL(wRowBase, VALUE='Pixel search method:')
  wExclusiveBase = WIDGET_BASE(wRowBase, /ROW, /EXCLUSIVE)
  wImmediateNeighborButton = WIDGET_BUTTON(wExclusiveBase, $
    VALUE='4-neighbor', $
    UVALUE='SEARCH_IMMEDIATE_NEIGHBORS', $
    UNAME=prefix+'immediate_neighbors')
  wAllNeighborButton = WIDGET_BUTTON(wExclusiveBase, $
    VALUE='8-neighbor', $
    UVALUE='SEARCH_ALL_NEIGHBORS', $
    UNAME=prefix+'all_neighbors')

  wFrameBase = WIDGET_BASE(wBase, /COLUMN, FRAME=1)
  wRowBase = WIDGET_BASE(wFrameBase, /ROW)
  wLabel = WIDGET_LABEL(wRowBase, VALUE='Threshold Range:  ')
  wThreshMinField = CW_FIELD(wRowBase, $
    TITLE='Min: ', $
    VALUE=(*pParentState).growThreshMin, $
    /INTEGER, /ALL_EVENTS, XSIZE=3, $
    UVALUE='THRESH_MIN', UNAME=prefix+'thresh_min')
  wThreshMaxField = CW_FIELD(wRowBase, $
    TITLE='Max: ', $
    VALUE=(*pParentState).growThreshMax, $
    /INTEGER, /ALL_EVENTS, XSIZE=3, $
    UVALUE='THRESH_MAX', UNAME=prefix+'thresh_max')

  wNonExclusiveBase = WIDGET_BASE(wFrameBase, /ROW, /NONEXCLUSIVE)
  wThreshROIButton = WIDGET_BUTTON(wNonExclusiveBase, $
    VALUE='Use source ROI threshold', UVALUE='THRESH_ROI', $
    UNAME=prefix+'thresh_roi')

  wStdDevField = CW_FIELD(wFrameBase, $
    TITLE='Standard deviation multiplier: ', $
    /FLOATING, /ALL_EVENTS, XSIZE=8, $
    VALUE=(*pParentState).growStdDevMult, $
    UVALUE='STDDEV_MULT', UNAME=prefix+'stddev_mult')

  wFrameBase = WIDGET_BASE(wBase, /COLUMN, FRAME=1)
  wRowBase = WIDGET_BASE(wFrameBase, /ROW)
  wLabel = WIDGET_LABEL(wRowBase, VALUE='For RGB images, use:')
  wExclusiveBase = WIDGET_BASE(wFrameBase, /COLUMN, /EXCLUSIVE)
  wUseLuminosityButton = WIDGET_BUTTON(wExclusiveBase, $
    VALUE='Luminosity', $
    UVALUE='USE_LUMINOSITY', $
    UNAME=prefix+'use_luminosity')
  wUseRedButton = WIDGET_BUTTON(wExclusiveBase, $
    VALUE='Red channel', $
    UVALUE='USE_RED', $
    UNAME=prefix+'use_red')
  wUseGreenButton = WIDGET_BUTTON(wExclusiveBase, $
    VALUE='Green channel', $
    UVALUE='USE_GREEN', $
    UNAME=prefix+'use_green')
  wUseBlueButton = WIDGET_BUTTON(wExclusiveBase, $
    VALUE='Blue channel', $
    UVALUE='USE_BLUE', $
    UNAME=prefix+'use_blue')

  wFrameBase = WIDGET_BASE(wBase, /COLUMN, FRAME=1)
  wRowBase = WIDGET_BASE(wFrameBase, /ROW)
  wLabel = WIDGET_LABEL(wRowBase, VALUE='Acceptance criteria:')

  wMaxCountField = CW_FIELD(wFrameBase, $
    TITLE='Maximum number of regions: ', $
    VALUE=(*pParentState).growMaxCount, $
    /INTEGER, /ALL_EVENTS, XSIZE=8, $
    UVALUE='MAX_COUNT', UNAME=prefix+'max_count')

  wRowBase = WIDGET_BASE(wFrameBase, /ROW)
  wMinAreaField = CW_FIELD(wRowBase, $
    TITLE='Minimum area per region:   ', $
    VALUE=(*pParentState).growMinArea, $
    /FLOATING, /ALL_EVENTS, XSIZE=8, $
    UVALUE='MIN_AREA', UNAME=prefix+'min_area')
  wAreaUnitLabel = WIDGET_LABEL(wRowBase, VALUE='(device units)')

  wNonExclusiveBase = WIDGET_BASE(wFrameBase, /ROW, /NONEXCLUSIVE)
  wAcceptAllButton = WIDGET_BUTTON(wNonExclusiveBase, $
    VALUE='Accept all regions', UVALUE='ACCEPT_ALL', $
    UNAME=prefix+'accept_all')

  wRowBase = WIDGET_BASE(wBase,/ROW)
  wButton = WIDGET_BUTTON(wRowBase, VALUE='Close', UVALUE='CLOSE', $
    UNAME=prefix+'close')

  (*pParentState).wROIGrowProps = wBase

  sState = {pParentState: pParentState, $
    wThreshMinField: wThreshMinField, $
    wThreshMaxField: wThreshMaxField, $
    wMaxCountField: wMaxCountField, $
    wMinAreaField: wMinAreaField, $
    wAreaUnitLabel: wAreaUnitLabel $
  }

  WIDGET_CONTROL, wBase, SET_UVALUE=sState
  if ((*pParentState).growAllNeighbors ne 0) then $
    WIDGET_CONTROL, wAllNeighborButton, SET_BUTTON=1 $
  else $
    WIDGET_CONTROL, wImmediateNeighborButton, SET_BUTTON=1

  if ((*pParentState).growThreshROI ne 0) then begin
    WIDGET_CONTROL, wThreshROIButton, SET_BUTTON=1
    WIDGET_CONTROL, wThreshMinField, SENSITIVE=0
    WIDGET_CONTROL, wThreshMaxField, SENSITIVE=0
  endif else begin
    WIDGET_CONTROL, wThreshROIButton, SET_BUTTON=0
    WIDGET_CONTROL, wThreshMinField, SENSITIVE=1
    WIDGET_CONTROL, wThreshMaxField, SENSITIVE=1
  endelse

  if ((*pParentState).growAcceptAll ne 0) then begin
    WIDGET_CONTROL, wAcceptAllButton, SET_BUTTON=1
    WIDGET_CONTROL, wMaxCountField, SENSITIVE=0
    WIDGET_CONTROL, wMinAreaField, SENSITIVE=0
    WIDGET_CONTROL, wAreaUnitLabel, SENSITIVE=0
  endif else begin
    WIDGET_CONTROL, wAcceptAllButton, SET_BUTTON=0
    WIDGET_CONTROL, wMaxCountField, SENSITIVE=1
    WIDGET_CONTROL, wMinAreaField, SENSITIVE=1
    WIDGET_CONTROL, wAreaUnitLabel, SENSITIVE=1
  endelse

  if ((*pParentState).growRGBUseLuminosity ne 0) then begin
    WIDGET_CONTROL, wUseLuminosityButton, SET_BUTTON=1
    WIDGET_CONTROL, wUseRedButton, SET_BUTTON=0
    WIDGET_CONTROL, wUseGreenButton, SET_BUTTON=0
    WIDGET_CONTROL, wUseBlueButton, SET_BUTTON=0
  endif else begin
    WIDGET_CONTROL, wUseLuminosityButton, SET_BUTTON=0
    case (*pParentState).growRGBChannel of
      0: begin
        WIDGET_CONTROL, wUseRedButton, SET_BUTTON=1
        WIDGET_CONTROL, wUseGreenButton, SET_BUTTON=0
        WIDGET_CONTROL, wUseBlueButton, SET_BUTTON=0
      end
      1: begin
        WIDGET_CONTROL, wUseRedButton, SET_BUTTON=0
        WIDGET_CONTROL, wUseGreenButton, SET_BUTTON=1
        WIDGET_CONTROL, wUseBlueButton, SET_BUTTON=0
      end
      2: begin
        WIDGET_CONTROL, wUseRedButton, SET_BUTTON=0
        WIDGET_CONTROL, wUseGreenButton, SET_BUTTON=0
        WIDGET_CONTROL, wUseBlueButton, SET_BUTTON=1
      end
    endcase
  endelse

  if not (*pParentState).modal then $
    WIDGET_CONTROL, wBase, MAP=0

  WIDGET_CONTROL, wBase, /REALIZE

  if not (*pParentState).modal then begin
    wbase_geom = WIDGET_INFO(wBase, /GEOMETRY)
    DEVICE, GET_SCREEN_SIZE=screen_size

    x = (screen_size[0] - wbase_geom.scr_xsize)/2
    x = x > 0

    y = (screen_size[1] - wbase_geom.scr_ysize)/2
    y = y > 0

    WIDGET_CONTROL, wBase, TLB_SET_XOFFSET=x, TLB_SET_YOFFSET=y
    WIDGET_CONTROL, wBase, MAP=1
  endif

  XMANAGER, "xd_xroiGrowProps", wBase, /NO_BLOCK
end

;------------------------------------------------------------------------------
function xd_xroi__Save, sEvent
  COMPILE_OPT idl2, hidden

  ; Save the current ROIs to a user-selected file.

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=pState

  oROIs = (*pState).oROIModel->Get(/ALL, COUNT=nROIs)
  if (nROIs gt 0) then begin
    fname = DIALOG_PICKFILE(GROUP=sEvent.top, FILE='regions.sav', $
      FILTER='*.sav', /WRITE)
    if (STRLEN(fname) gt 0) then begin
      xd_xroi__SetROI, pState, OBJ_NEW()
      (*pState).oROIModel->Remove, oROIs
      (*pState).oROIGroup->Remove, oROIs
      SAVE, oROIs, FILENAME=fname
      (*pState).oROIModel->Add, oROIs
      (*pState).oROIGroup->Add, oROIs
      if OBJ_VALID((*pState).oSelROI) then $
        xd_xroi__SetROI, pState, oSelROI, /SET_LIST_SELECT
    endif
  endif

  RETURN, 0 ; "Swallow" event.
end

;------------------------------------------------------------------------------
pro xd_xroi__Quit, sEvent
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, sEvent.top, /DESTROY
end

;------------------------------------------------------------------------------
pro xd_xroiPickColor_event, event
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, event.top, GET_UVALUE=pState
  pLeaderState = (*pState).pLeaderState

  case event.id of
    (*pState).wROISelColor: begin
      (*pLeaderState).sel_rgb = [event.r, event.g, event.b]
      xd_xroi__SetSelectionVisualColor, pLeaderState, $
        [event.r, event.g, event.b]
      if OBJ_VALID((*pLeaderState).oSelROI) then begin
        (*pLeaderState).oSelROI->SetProperty, $
          COLOR=(*pLeaderState).sel_rgb
        (*pLeaderState).oWindow->Draw, (*pLeaderState).oView
      endif
    endcase
    (*pState).wROIColor: begin
      (*pLeaderState).roi_rgb = [event.r, event.g, event.b]
      oROIs = (*pLeaderState).oROIModel->Get(/ALL, COUNT=count)
      for i=0,count-1 do begin
        oROIs[i]->SetProperty, COLOR=(*pLeaderState).roi_rgb
      endfor
      if OBJ_VALID((*pLeaderState).oSelROI) then begin
        (*pLeaderState).oSelROI->SetProperty, $
          COLOR=(*pLeaderState).sel_rgb
      endif
      if count gt 0 then begin
        (*pLeaderState).oWindow->Draw, (*pLeaderState).oView
      endif
    endcase
    (*pState).wOK: begin
      WIDGET_CONTROL, event.top, /DESTROY
    end
    (*pState).wCancel: begin
      oROIs = (*pLeaderState).oROIModel->Get(/ALL, COUNT=count)
      for i=0,count-1 do begin
        oROIs[i]->SetProperty, COLOR=(*pState).roi_rgb
      endfor
      if OBJ_VALID((*pLeaderState).oSelROI) then begin
        (*pLeaderState).oSelROI->SetProperty, COLOR=(*pState).sel_rgb
      endif
      xd_xroi__SetSelectionVisualColor, pLeaderState, (*pState).sel_rgb

      (*pLeaderState).sel_rgb = (*pState).sel_rgb
      (*pLeaderState).roi_rgb = (*pState).roi_rgb

      (*pLeaderState).oWindow->Draw, (*pLeaderState).oView
      WIDGET_CONTROL, event.top, /DESTROY
    endcase
    else:
  endcase
end

;------------------------------------------------------------------------------
pro xd_xroiPickColor__Cleanup, wID
  WIDGET_CONTROL, wID, GET_UVALUE=pState

  PTR_FREE, pState
end

;------------------------------------------------------------------------------
pro xd_xroiPickColor, group_leader
  COMPILE_OPT idl2, hidden

  WIDGET_CONTROL, group_leader, GET_UVALUE=pLeaderState

  if WIDGET_INFO((*pLeaderState).wPickColor, /VALID_ID) then begin
    WIDGET_CONTROL, (*pLeaderState).wPickColor, /SHOW
    RETURN
  end

  prefix = 'xd_xroiPickColor:'
  tlb = WIDGET_BASE( $  ; Top-Level Base.
    /COLUMN, $
    TITLE='ROI Outline Colors', $
    GROUP_LEADER=group_leader, $
    UNAME=prefix + 'tlb', $
    MODAL=(*pLeaderState).modal $
    )

  wRowBase = WIDGET_BASE(tlb, /ROW, XPAD=0, YPAD=0)
  wColumnBase = WIDGET_BASE(wRowBase, /COLUMN)
  void = WIDGET_LABEL(wColumnBase, VALUE='Selected Outline Color:')
  wROISelColor = CW_RGBSlider( $
    wColumnBase, $
    VALUE=(*pLeaderState).sel_rgb, $
    GRAPHICS_LEVEL=2, $
    /HLS, $
    UNAME=prefix + 'ROISelColor', $
    /FRAME $
    )
  wColumnBase = WIDGET_BASE(wRowBase, /COLUMN)
  void = WIDGET_LABEL(wColumnBase, VALUE='Unselected Outline Color:')
  wROIColor = CW_RGBSlider( $
    wColumnBase, $
    VALUE=(*pLeaderState).roi_rgb, $
    GRAPHICS_LEVEL=2, $
    /HLS, $
    UNAME=prefix + 'ROIColor', $
    /FRAME $
    )

  wRowBase = WIDGET_BASE(tlb, /ROW, /GRID)
  wOK = WIDGET_BUTTON(wRowBase, VALUE='OK', UNAME=prefix + 'OK')
  wCancel = WIDGET_BUTTON( $
    wRowBase, $
    VALUE='Cancel', $
    UNAME=prefix + 'cancel' $
    )

  if not (*pLeaderState).modal then $
    WIDGET_CONTROL, tlb, MAP=0

  WIDGET_CONTROL, tlb, /REALIZE

  if not (*pLeaderState).modal then begin
    tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
    leader_geom = WIDGET_INFO(group_leader, /GEOMETRY)
    DEVICE, GET_SCREEN_SIZE=screen_size

    x = leader_geom.xoffset
    x = x < (screen_size[0] - tlb_geom.scr_xsize)
    x = x > 0

    y = leader_geom.scr_ysize[0] + leader_geom.yoffset
    y = y < (screen_size[1] - tlb_geom.scr_ysize)
    y = y > 0

    WIDGET_CONTROL, tlb, TLB_SET_XOFFSET=x, TLB_SET_YOFFSET=y
    WIDGET_CONTROL, tlb, MAP=1
  endif

  WIDGET_CONTROL, tlb, SET_UVALUE=PTR_NEW({ $
    sel_rgb: (*pLeaderState).sel_rgb, $
    roi_rgb: (*pLeaderState).roi_rgb, $
    wROISelColor: wROISelColor, $
    wROIColor: wROIColor, $
    wOK: wOK, $
    wCancel:wCancel, $
    pLeaderState: pLeaderState $
  })

  (*pLeaderState).wPickColor = tlb
  XMANAGER, $
    'xd_xroiPickColor', $
    tlb, $
    CLEANUP='xd_xroiPickColor__Cleanup', $
    /NO_BLOCK
end

;------------------------------------------------------------------------------
pro xd_xroi__cleanup, wID
  COMPILE_OPT hidden

  WIDGET_CONTROL, wID, GET_UVALUE=pState

  if (*pState).demo_system then begin
    if WIDGET_INFO((*pState).group_leader, /VALID) then $
      WIDGET_CONTROL, (*pState).group_leader, /MAP
  end

  if OBJ_VALID((*pState).oSelROI) then $
    (*pState).oSelROI->SetProperty, COLOR=(*pState).roi_rgb

  if PTR_VALID((*pState).pR) then begin
    (*pState).oPalette->GetProperty, RED_VALUES=red_values
    *(*pState).pR = red_values
  end
  if PTR_VALID((*pState).pG) then begin
    (*pState).oPalette->GetProperty, GREEN_VALUES=green_values
    *(*pState).pG = green_values
  end
  if PTR_VALID((*pState).pB) then begin
    (*pState).oPalette->GetProperty, BLUE_VALUES=blue_values
    *(*pState).pB = blue_values
  end

  if PTR_VALID((*pState).pSelRGB) then begin
    *(*pState).pSelRGB = (*pState).sel_rgb
  end

  if PTR_VALID((*pState).pROIRGB) then begin
    *(*pState).pROIRGB = (*pState).roi_rgb
  end

  if (PTR_VALID((*pState).pTools)) then begin
    for i=0,N_ELEMENTS(*(*pState).pTools)-1 do begin
      oTool = (*(*pState).pTools)[i]
      if (OBJ_VALID(oTool)) then $
        void = oTool->DoAction('/SERVICES/SHUTDOWN')
    endfor
  endif

  ;
  ;   Clean up (*pState).pImg, if appropriate.
  ;
  ;   In order that *(*pState).pImg always be available when we need it, we
  ;   want to clean up (*pState).pImg at the last opportunity.  If our xd_xroi's
  ;   call to XMANAGER did not stop control flow, then by the time we get
  ;   here, this is that last opportunity.  To determine if our call to
  ;   XMANAGER did not stop control flow, we check whether (*pState).pR
  ;   has been cleaned up.
  ;
  if not PTR_VALID((*pState).pR) then begin
    PTR_FREE, (*pState).pImg
  end
  ;
  ;   Clean up other heap variables.
  ;
  if OBJ_VALID((*pState).oRegionsOut) then begin
    (*pState).oROIModel->Remove, /ALL
    (*pState).oROIGroup->Remove, /ALL
  endif
  if (OBJ_VALID((*pState).oTransScaleVisual) ne 0) then begin
    (*pState).oTransScaleVisual->GetProperty, UVALUE=pTSState
    PTR_FREE, pTSState
  endif
  for i=0,n_tags(*pState)-1 do begin
    case size((*pState).(i), /TNAME) of
      'POINTER': begin
        if (*pState).(i) ne (*pState).pR $
          and (*pState).(i) ne (*pState).pG $
          and (*pState).(i) ne (*pState).pB $
          and (*pState).(i) ne (*pState).pSelRGB $
          and (*pState).(i) ne (*pState).pROIRGB $
          and (*pState).(i) ne (*pState).pImg then begin
          ptr_free, (*pState).(i)
        endif
      endcase
      'OBJREF': begin
        if (*pState).(i) ne (*pState).oRegionsOut $
          and (*pState).(i) ne (*pState).oRejected $
          and (*pState).(i) ne (*pState).oSelROI $
          and (*pState).(i) ne (*pState).oRegionsIn then begin
          obj_destroy, (*pState).(i)
        endif
      endcase
      else:
    endcase
  endfor

  ; Cannot free up the state for modal, since we need to do more cleanup.
  if not (*pState).modal then $
    ptr_free, pState
end

;name: mask_display
;input: new_mask_vol_cube
;       location3D
;       imageSize
;
FUNCTION mask_display, new_mask_vol_cube, location3D, imageSize, direction
  IF STRCMP(direction, 'xy') THEN BEGIN
    ;print, xd_sState.location3D[2], 'location'
    ; Create a mask for the region.
    xyMask_small = new_mask_vol_cube[*,*,location3D[2]*(SIZE(new_mask_vol_cube))[3]/imageSize]
    xyMask3d = (REVERSE(CONGRID(xyMask_small, imageSize*3, imageSize*3, 1), 2)) GT 0
    xyMask2d = DBLARR(imageSize*3, imageSize*3)
    xyMask2d[*, *] = xyMask3d[*, *, 0]
    strucElem = REPLICATE(1, 3, 3)
    xyMask = ERODE(DILATE(TEMPORARY(xyMask2D), strucElem), strucElem)
    ; Contour the mask.
    oContour = OBJ_NEW('IDLgrContour', xyMask, COLOR =[255,127,127], /PLANAR, GEOMZ=0, C_VALUE=INDGEN(20), N_LEVEL = 1, /FILL, ALPHA_CHANNEL = 0.5)
  ENDIF

  IF STRCMP(direction, 'yz') THEN BEGIN
    ; Create a mask for the region.
    yzMask_small = new_mask_vol_cube[location3D[0]*(SIZE(new_mask_vol_cube))[1]/imageSize, *,*]
    yzMask3d = CONGRID(yzMask_small, 1, imageSize*3, imageSize*3) GT 0
    yzMask2d = DBLARR(imageSize*3, imageSize*3)
    yzMask2d[*, *] = yzMask3d[0, *, *]
    strucElem = REPLICATE(1, 3, 3)
    yzMask = ERODE(DILATE(TEMPORARY(yzMask2d), strucElem), strucElem)
    ; Contour the mask.
    oContour = OBJ_NEW('IDLgrContour', yzMask, COLOR =[127,255,127], /PLANAR, GEOMZ=0, C_VALUE=INDGEN(20), N_LEVEL = 1, /FILL, ALPHA_CHANNEL = 0.5)
  ENDIF

  IF STRCMP(direction, 'xz') THEN BEGIN
    ; Create a mask for the region.
    xzMask_small = new_mask_vol_cube[*, location3D[1]*(SIZE(new_mask_vol_cube))[2]/imageSize, *,*]
    xzMask3d = CONGRID(xzMask_small, imageSize*3, 1, imageSize*3) GT 0
    xzMask2d = DBLARR(imageSize*3, imageSize*3)
    xzMask2d[*, *] = xzMask3d[*,0, *]
    strucElem = REPLICATE(1, 3, 3)
    xzMask = ERODE(DILATE(TEMPORARY(xzMask2d), strucElem), strucElem)
    ; Contour the mask.
    oContour = OBJ_NEW('IDLgrContour', xzMask, COLOR =[255, 127, 12], /PLANAR, GEOMZ=0, C_VALUE=INDGEN(20), N_LEVEL = 1, /FILL, ALPHA_CHANNEL = 0.5)
  ENDIF
  RETURN, oContour
END


PRO extract_mask_diff_to_roi, pState
  (*pState).oImage2->getProperty, DATA=Mask_ROI
  (*pState).oImage1->getProperty, DATA=Mask_new
  (*pState).oImage0->getProperty, DATA=Mask_ori
  Mask_diff = Mask_new - Mask_ori
  IF MIN(Mask_diff) EQ -1. THEN Mask_diff = Mask_ori - Mask_new
  Mask_diff = (Mask_diff GT 100) OR (Mask_ROI GT 100)
  ;iimage, Mask_diff
  CONTOUR, Mask_diff, LEVELS=1, PATH_INFO=pathInfo, PATH_XY=pathXY, $
    /PATH_DATA_COORD

  if (N_ELEMENTS(pathInfo) eq 0) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      ['New mask region could not', $
      'be successfully contoured.'], $
      DIALOG_PARENT=(*pState).wBase)
    RETURN
  endif

  ; Select all exterior contours.
  iids = WHERE(pathInfo[*].high_low eq 1, nContours)
  if (nContours eq 0) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      ['Grown region results in', $
      'no exterior contours.'], $
      DIALOG_PARENT=(*pState).wBase)
    RETURN
  endif

  ; Create a temporary region group to store grown regions.
  oGrowROIGroup = OBJ_NEW('IDLanROIGroup')

  ; Add a new region for each contour.
  iNAccept = 0
  ;if the last one has a much small ROI, then remove
  IF pathInfo[nContours-1].N LT 0.1*pathInfo[0].N THEN nContours = nContours-1
  for i=0,nContours-1 do begin
    IF pathInfo[i].N LT 0.1*pathInfo[0].N THEN break
    contId = iids[i]
    iStart = pathInfo[contId].offset
    iFinish = iStart + pathInfo[contId].n - 1

    if (iNAccept gt 0) then begin
      ; Determine if (the first vertex of) this region is contained
      ; within any of the previously accepted regions.
      ; If so, bypass.
      ; Otherwise, accept it as a new valid region.
      containStatus = oGrowROIGroup->ContainsPoints( $
        pathXY[0,iStart], pathXY[0,iStart])
      bAccept = (containStatus[0] eq 0) ? 1b : 0b
    endif else $
      bAccept = 1b

    ; If requested, do not accept if beyond maximum count.
    if (((*pState).growAcceptAll eq 0) and $
      (iNAccept ge (*pState).growMaxCount)) then $
      bAccept = 0b

    if (bAccept ne 0) then begin
      ; Grab the contour path.
      growROIData = FLTARR(3,pathInfo[contId].n)
      growROIData[0,*] = pathXY[0, iStart:iFinish]
      growROIData[1,*] = pathXY[1, iStart:iFinish]

      oGrowROI = OBJ_NEW('IDLgrROI', growROIData, STYLE=2)

      ;update the mask
      IF STRCMP((*pState).Modify_ROI, 'BORDER') THEN ac = 0.000001
      IF STRCMP((*pState).Modify_ROI, 'ADD') THEN ac = 0.000004
      IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN ac = 0.000008
      IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN ac = 0.000012
      IF STRCMP((*pState).Modify_ROI, 'ADDNEW2') THEN ac = 0.000016

      IF STRCMP((*pState).CONTORINDE, 'CONTINUOUS') THEN linestyle = 1
      IF STRCMP((*pState).CONTORINDE, 'DISCRETE') THEN linestyle = 0
      oGrowROI->setProperty, LINESTYLE = linestyle

      location3D = (*pState).location3D

      IF STRCMP((*pState).direction, 'xy') THEN BEGIN
        Location = location3D[2]
        oGrowROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
      ENDIF

      IF STRCMP((*pState).direction, 'yz') THEN BEGIN
        Location = location3D[0]
        oGrowROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
      ENDIF

      IF STRCMP((*pState).direction, 'xz') THEN BEGIN
        Location = location3D[1]
        oGrowROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
      ENDIF

      ; If requested, do not accept if less than minimum area.
      if ((*pState).growAcceptAll eq 0) then begin
        isOk = oGrowROI->ComputeGeometry(AREA=a)
        if (a lt (*pState).growMinArea) then $
          bAccept = 0b
      endif

      ; The contour has passed all acceptance criteria.
      if (bAccept ne 0) then begin
        oGrowROI->SetProperty, NAME=xd_xroi__GenerateName(*pState)
        (*pState).oModel->Remove, oGrowROI
        (*pState).oROIModel->Add, oGrowROI
        (*pState).oROIGroup->Add, oGrowROI
        if OBJ_VALID((*pState).oRegionsOut) then $
          (*pState).oRegionsOut->Add, oGrowROI

        oGrowROIGroup->Add, oGrowROI

        xd_xroi__SetROI, pState, oGrowROI, /UPDATE_LIST, /SET_LIST_SELECT
        iNAccept = iNAccept + 1
      endif else $
        OBJ_DESTROY, oGrowROI
    endif
  endfor

  if (iNAccept eq 0) then begin
    void = DIALOG_MESSAGE(/INFORM, $
      ['No contours of the grown region', $
      'passed the acceptance criteria.'], $
      DIALOG_PARENT=(*pState).wBase)
    OBJ_DESTROY, oGrowROIGroup
    RETURN
  endif

  ; Select the first grown ROI.
  oGrowROI =  oGrowROIGroup->Get(POSITION=0)
  xd_xroi__SetROI, pState, oGrowROI, /SET_LIST_SELECT

  ; Free temporary region group.
  oGrowROIGroup->Remove, /ALL
  OBJ_DESTROY, oGrowROIGroup


end

;name: modify mask by roi
;input:
PRO modify_mask_by_ROI, pState


  ;update the mask
  IF STRCMP((*pState).Modify_ROI, 'BORDER') THEN ac = 0.000001
  IF STRCMP((*pState).Modify_ROI, 'ADD') THEN ac = 0.000004
  IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN ac = 0.000008
  IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN ac = 0.000012
  IF STRCMP((*pState).Modify_ROI, 'ADDNEW2') THEN ac = 0.000016
  
  IF STRCMP((*pState).CONTORINDE, 'CONTINUOUS') THEN linestyle = 1
  IF STRCMP((*pState).CONTORINDE, 'DISCRETE') THEN linestyle = 0
  (*pState).oSelROI->setProperty, LINESTYLE = linestyle

  (*pState).oImage->getProperty, DATA=Img
  roi_mask_dims = SIZE(Img, /DIMENSIONS)
  IF OBJ_VALID((*pState).oSelROI) EQ 0 THEN RETURN
  maskResult = (*pState).oSelROI -> ComputeMask(DIMENSIONS = roi_mask_dims)
  dims = SIZE((*pState).new_mask_vol_cube, /DIMENSIONS)
  location3D = (*pState).location3D

  IF STRCMP((*pState).direction, 'xy') THEN BEGIN
    Location = location3D[2]
    (*pState).oSelROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
  ENDIF

  IF STRCMP((*pState).direction, 'yz') THEN BEGIN
    Location = location3D[0]
    (*pState).oSelROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
  ENDIF

  IF STRCMP((*pState).direction, 'xz') THEN BEGIN
    Location = location3D[1]
    (*pState).oSelROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
  ENDIF

  i=Location
  pre_mask = BYTARR(dims[0], dims[1])
  IF STRCMP((*pState).direction, 'xy') THEN BEGIN
    pre_mask[*, *] = (*pState).new_mask_vol_cube[*,*,i*dims[2]/(*pState).imageSize]
    mask_ROI = ((REVERSE(CONGRID(maskResult, dims[0], dims[1]), 2)) GT 0)
    IF STRCMP((*pState).Modify_ROI, 'BORDER') THEN update_mask = mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADD') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADDNEW2') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN update_mask = pre_mask AND (1-mask_ROI)
    (*pState).new_mask_vol_cube[*,*,i*dims[2]/(*pState).imageSize] = update_mask
  ENDIF

  IF STRCMP((*pState).direction, 'yz') THEN BEGIN
    pre_mask[*, *] = (*pState).new_mask_vol_cube[i*dims[0]/(*pState).imageSize,*,*]
    mask_ROI = ((CONGRID(maskResult, dims[1], dims[2])) GT 0)
    IF STRCMP((*pState).Modify_ROI, 'BORDER') THEN update_mask = mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADD') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADDNEW2') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN update_mask = pre_mask AND (1-mask_ROI)
    (*pState).new_mask_vol_cube[i*dims[0]/(*pState).imageSize,*,*] = update_mask
  ENDIF

  IF STRCMP((*pState).direction, 'xz') THEN BEGIN
    pre_mask[*, *] = (*pState).new_mask_vol_cube[*,i*dims[1]/(*pState).imageSize,*]
    mask_ROI = ((CONGRID(maskResult, dims[0], dims[2])) GT 0)
    IF STRCMP((*pState).Modify_ROI, 'BORDER') THEN update_mask = mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADD') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADDNEW1') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'ADDNEW2') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP((*pState).Modify_ROI, 'REMOVE') THEN update_mask = pre_mask AND (1-mask_ROI)
    (*pState).new_mask_vol_cube[*,i*dims[1]/(*pState).imageSize,*] = update_mask
  ENDIF

  ;update the mask image
  currdim = SIZE((*pState).svol_HU_cube, /DIMENSIONS)
  maskdim = SIZE((*pState).new_mask_vol_cube, /DIMENSIONS)
  imageSize = (*pState).imageSize * 3
  Mask = DBLARR(imageSize, imageSize)
  IF STRCMP((*pState).direction, 'xy') THEN BEGIN
    Mask3D = REVERSE((CONGRID((*pState).new_mask_vol_cube[*,*,location3D[2]*maskdim[2]/currdim[2]], $
      imageSize, imageSize,1) GT 0)*250, 2)
    Mask[*, *] = Mask3d[*, *, 0]
  ENDIF

  IF STRCMP((*pState).direction, 'yz') THEN BEGIN
    Mask3D = (CONGRID((*pState).new_mask_vol_cube[location3D[0]*maskdim[0]/currdim[0],*,*], 1, $
      imageSize, imageSize) GT 0)*140
    Mask[*, *] = Mask3d[0, *, *]
  ENDIF

  IF STRCMP((*pState).direction, 'xz') THEN BEGIN
    Mask3D = (CONGRID((*pState).new_mask_vol_cube[*,location3D[1]*maskdim[1]/currdim[1],*], $
      imageSize, 1, imageSize) GT 0)*210
    Mask[*, *] = Mask3d[*, 0, *]
  ENDIF

  (*pState).oImage1->SetProperty, DATA=Mask
  ;GJ, this is used to store ROI
  (*pState).oImage2->SetProperty, DATA=maskResult

  ;    iimage, maskResult
  ;delete the variables
  Mask = !NULL
  Mask3D = !NULL

END

;name: modify mask by roi no pState
;input:
FUNCTION modify_mask_by_ROI_nopState, Modify_ROI, oImage, oSelROI, new_mask_vol_cube, mask_ROI, location3D, direction, imageSize, currdim

  ;print, 'currdim = ', currdim

  ;update the mask
  IF STRCMP(Modify_ROI, 'ADDNEW1') THEN ac = 0.000012
  IF STRCMP(Modify_ROI, 'ADD') THEN ac = 0.000004
  IF STRCMP(Modify_ROI, 'REMOVE') THEN ac = 0.000008

  IF STRCMP(Modify_ROI, 'BORDER') THEN ac = 0.000001
  IF STRCMP(Modify_ROI, 'ADDNEW2') THEN ac = 0.000016

  oImage->getProperty, DATA=Img
  roi_mask_dims = SIZE(Img, /DIMENSIONS)
  IF OBJ_VALID(oSelROI) EQ 0 THEN RETURN, 0
  maskResult = oSelROI -> ComputeMask(DIMENSIONS = roi_mask_dims)
  dims = SIZE(new_mask_vol_cube, /DIMENSIONS)
  location3D = location3D

  IF STRCMP(direction, 'xy') THEN BEGIN
    Location = location3D[2]
    oSelROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
  ENDIF

  IF STRCMP(direction, 'yz') THEN BEGIN
    Location = location3D[0]
    oSelROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
  ENDIF

  IF STRCMP(direction, 'xz') THEN BEGIN
    Location = location3D[1]
    oSelROI->setProperty, ALPHA_CHANNEL=Location/10000.+ac
  ENDIF

  i=Location
  pre_mask = BYTARR(dims[0], dims[1])
  IF STRCMP(direction, 'xy') THEN BEGIN
    pre_mask[*, *] = new_mask_vol_cube[*,*,i*dims[2]/imageSize]
    mask_ROI = ((REVERSE(CONGRID(maskResult, dims[0], dims[1]), 2)) GT 0)
    IF STRCMP(Modify_ROI, 'BORDER') THEN update_mask = mask_ROI
    IF STRCMP(Modify_ROI, 'ADD') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'ADDNEW1') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'ADDNEW2') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'REMOVE') THEN update_mask = pre_mask AND (1-mask_ROI)
    new_mask_vol_cube[*,*,i*dims[2]/imageSize] = update_mask
  ENDIF

  IF STRCMP(direction, 'yz') THEN BEGIN
    pre_mask[*, *] = new_mask_vol_cube[i*dims[0]/imageSize,*,*]
    mask_ROI = ((CONGRID(maskResult, dims[1], dims[2])) GT 0)
    IF STRCMP(Modify_ROI, 'BORDER') THEN update_mask = mask_ROI
    IF STRCMP(Modify_ROI, 'ADD') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'ADDNEW1') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'ADDNEW2') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'REMOVE') THEN update_mask = pre_mask AND (1-mask_ROI)
    new_mask_vol_cube[i*dims[0]/imageSize,*,*] = update_mask
  ENDIF

  IF STRCMP(direction, 'xz') THEN BEGIN
    pre_mask[*, *] = new_mask_vol_cube[*,i*dims[1]/imageSize,*]
    mask_ROI = ((CONGRID(maskResult, dims[0], dims[2])) GT 0)
    IF STRCMP(Modify_ROI, 'BORDER') THEN update_mask = mask_ROI
    IF STRCMP(Modify_ROI, 'ADD') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'ADDNEW1') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'ADDNEW2') THEN update_mask = pre_mask OR mask_ROI
    IF STRCMP(Modify_ROI, 'REMOVE') THEN update_mask = pre_mask AND (1-mask_ROI)
    new_mask_vol_cube[*,i*dims[1]/imageSize,*] = update_mask
  ENDIF

  ;update the mask image
  maskdim = SIZE(new_mask_vol_cube, /DIMENSIONS)
  imageSize = imageSize * 3
  Mask = DBLARR(imageSize, imageSize)
  IF STRCMP(direction, 'xy') THEN BEGIN
    Mask3D = REVERSE((CONGRID(new_mask_vol_cube[*,*,location3D[2]*maskdim[2]/currdim[2]], $
      imageSize, imageSize,1) GT 0)*250, 2)
    Mask[*, *] = Mask3d[*, *, 0]
  ENDIF

  IF STRCMP(direction, 'yz') THEN BEGIN
    Mask3D = (CONGRID(new_mask_vol_cube[location3D[0]*maskdim[0]/currdim[0],*,*], 1, $
      imageSize, imageSize) GT 0)*140
    Mask[*, *] = Mask3d[0, *, *]
  ENDIF

  IF STRCMP(direction, 'xz') THEN BEGIN
    Mask3D = (CONGRID(new_mask_vol_cube[*,location3D[1]*maskdim[1]/currdim[1],*], $
      imageSize, 1, imageSize) GT 0)*210
    Mask[*, *] = Mask3d[*, 0, *]
  ENDIF

  ;change mask_ROI to original size
  mask_ROI = maskResult
  ;  iimage, mask_ROI
  ;
  ;delete the variables
  Mask3D = !NULL

  RETURN, Mask

END

FUNCTION donut_shaped_roi, roi_vertices_inner, roi_vertices_outer

  ;; Build the new vertices file for the donut:
  size_inner_vertices = SIZE(roi_vertices_inner)
  n_inner_vertices = size_inner_vertices[2]
  size_outer_vertices = SIZE(roi_vertices_outer)
  n_outer_vertices = size_outer_vertices[2]
  n_donut_vertices = n_inner_vertices + n_outer_vertices + 2
  ;            roi_vertices_inner[0, n_inner_vertices - 1] = roi_vertices_inner[0, n_inner_vertices - 1] - 1
  y0 = roi_vertices_inner[1, n_inner_vertices - 1]
  x0 = roi_vertices_inner[0, n_inner_vertices - 1]
  FOR ii = 0, n_outer_vertices - 2 DO BEGIN
    x1 = roi_vertices_outer[0, ii]
    y1 = roi_vertices_outer[1, ii]
    x2 = roi_vertices_outer[0, ii + 1]
    y2 = roi_vertices_outer[1, ii + 1]
    IF ((y1 GT y0) OR (y2 GT y0)) AND (((x1 LT x0) AND (x2 GE x0)) OR ((x1 GE x0) AND (x2 LT x0))) THEN BEGIN
      point_1 = [x1, y1, ii]
      point_2 = [x2, y2, ii + 1]
      BREAK
    ENDIF
  ENDFOR
  x1 = DOUBLE(point_1[0])
  y1 = DOUBLE(point_1[1])
  x2 = DOUBLE(point_2[0])
  y2 = DOUBLE(point_2[1])
  ii = point_1[2]
  kk = (y1 - y2)/(x1 - x2)
  yo_x0 = ROUND(y1 + kk * (x0 - x1))
  ;            yo_x0_1 = ROUND(y1 + kk * (x0 - x1 - 1.D))
  inserted_outer_vertices = INTARR(2, n_outer_vertices + 1)
  inserted_outer_vertices[0:1, 0] = [x0, yo_x0]
  ;inserted_outer_vertices[1, 0] = yo_x0
  inserted_outer_vertices[0:1, n_outer_vertices] = [x0, yo_x0]
  ;inserted_outer_vertices[1, n_vertices_outer] = yo_x0
  inserted_outer_vertices[0:1, 1:n_outer_vertices - ii - 1] = roi_vertices_outer[0:1,ii + 1:n_outer_vertices - 1]
  inserted_outer_vertices[0:1, n_outer_vertices - ii:n_outer_vertices - 1] = roi_vertices_outer[0:1, 1:ii]
  inserted_outer_vertices = REVERSE(inserted_outer_vertices, 2)
  roi_vertices = LONARR(2, n_donut_vertices)
  roi_vertices[0:1, 0:n_inner_vertices - 1] = roi_vertices_inner
  roi_vertices[0:1, n_inner_vertices:n_donut_vertices - 2] = inserted_outer_vertices
  roi_vertices[0:1, n_donut_vertices - 1] = [x0, y0]

  RETURN, roi_vertices
END


;------------------------------------------------------------------------------
pro xd_xroi, $
  img, $                  ; IN/OUT: (opt) image data
  ;    xd_sState, $
  svol_HU_cube, threshold, new_mask_vol_cube, location3D, imageSize, icon_dir, $
  ObjectIndex_value=ObjectIndex_value, $;IN: (opt)
  direction = direction, $     ;'xy', 'yz', 'xz'
  r, $                    ; IN/OUT: (opt) Red.  256-element byte array.
  g, $                    ; IN/OUT: (opt) Green. 256-element byte array.
  b, $                    ; IN/OUT: (opt) Blue. 256-element byte array.
  ROI_COLOR=roi_color, $  ; IN/OUT: (opt) [r,g,b]
  ROI_SELECT_COLOR=roi_select_color, $  ; IN/OUT: (opt) [r,g,b]
  BLOCK=block, $          ; IN: (opt)
  FLOATING=floating, $    ; IN: (opt) property of self's top-level base.
  GROUP=group_leader, $   ; IN: (opt)
  MODAL=modal, $          ; IN: (opt)
  REGIONS_IN=regions_in, $        ; IN: (opt) array of objects
  REGIONS_OUT=regions_out, $      ; OUT: (opt) array of objects
  REJECTED=rejected, $            ; OUT: (opt) array of objects
  ROI_GEOMETRY=roi_geometry, $    ; OUT: (opt)
  STATISTICS=statistics, $        ; OUT: (opt)
  TITLE=title, $          ; IN: (opt)
  TOOLS=tools, $          ; IN: (opt)
  RENDERER=renderer, $    ; IN: (opt) 0=Native OpenGL, 1=Software
  UNAME=uname, $
  TEST=test, $
  APPTLB=wBase, $         ; IN: (opt) used by demo system
  RECORD_TO_FILENAME=record_to_filename, $ ; IN: (opt) used by demo system
  DEBUG=debug, $
  X_SCROLL_SIZE=xScrollSizeIn, $
  Y_SCROLL_SIZE=yScrollSizeIn

  ON_ERROR, KEYWORD_SET(debug) ? 0 : 2

  if N_ELEMENTS(group_leader) ne 0 then begin
    if not WIDGET_INFO(group_leader, /VALID_ID) then $
      MESSAGE, 'Specified Group Leader is not valid.'
  endif else begin
    if KEYWORD_SET(floating) then $
      MESSAGE, 'floating xd_xroi requires a group leader.'
    if KEYWORD_SET(modal) then $
      MESSAGE, 'modal xd_xroi requires a group leader.'
  endelse

  demo_system = ARG_PRESENT(wBase)

  if N_ELEMENTS(tools) eq 0 then begin
    _tools = STRUPCASE(['CircularMask Draw', 'Freehand Draw', 'Rectangle', 'Ellipse', $
      'Polygon Draw'])
  endif else begin
    _tools = STRUPCASE(tools)
  endelse

  if SIZE(_tools, /TNAME) ne 'STRING' then $
    MESSAGE, 'TOOLS must be a string or string array.'
  if N_ELEMENTS(_tools) gt 7 then $
    MESSAGE, 'TOOLS cannot have more than seven elements.'
  for i=0,N_ELEMENTS(_tools)-1 do begin
    case STRUPCASE(_tools[i]) of
      'CIRCULARMASK DRAW':
      'FREEHAND DRAW':
      'RECTANGLE':
      'ELLIPSE':
      'POLYGON DRAW':
      ;            'SELECTION':
      ;            'TRANSLATE-SCALE':
      else: MESSAGE, 'Unknown TOOLS value: ' + _tools[i]
    endcase
  endfor

  ; Default to software renderer, since many native OpenGL implementations
  ; are not optimized for image display.
  if (N_ELEMENTS(renderer) eq 0) then renderer = 1

  if KEYWORD_SET(test) then begin
    read_jpeg, filepath('muscle.jpg', SUBDIR=['examples','data']), $
      img, /grayscale
  end

  ; Verify that provided image (if any) has valid dimensions.
  if (N_ELEMENTS(img) gt 0) then begin
    nDims = SIZE(img, /N_DIMENSIONS)
    case nDims of
      2: begin ; 8-bit image
        dimensions = SIZE(img, /DIMENSIONS)
        image_is_8bit = 1b
      end

      3: begin ; RGB image
        allDims = SIZE(img, /DIMENSIONS)
        iInterleave = (WHERE(allDims eq 3))[0]
        case iInterleave of
          0: dimensions = allDims[1:2]
          1: dimensions = [allDims[0],allDims[2]]
          2: dimensions = allDims[0:1]
          else: MESSAGE, $
            'Image must have dimensions [3,m,n], [m,3,n], or [m,n,3].'
        endcase
        image_is_8bit = 0b
      end

      else: MESSAGE, 'Image data must have 2 or 3 dimensions.'
    endcase
  endif else $
    nDims = 0

  ; Prompt user for image file until success or cancel.
  while nDims eq 0 do begin
    if not DIALOG_READ_IMAGE( $
      RED=r, $
      GREEN=g, $
      BLUE=b, $
      IMAGE=img $
      ) $
      then $
      RETURN
    if img[0] eq -1 then $
      RETURN ; User "opened" a non-image file or directory.

    if r[0] eq -1 then begin ; e.g. IDL's examples/data/muscle.jpg
      PTR_FREE, PTR_NEW(TEMPORARY(r))
      PTR_FREE, PTR_NEW(TEMPORARY(g))
      PTR_FREE, PTR_NEW(TEMPORARY(b))
    end

    ; Verify that the selected image has valid dimensions.
    nDims = SIZE(img, /N_DIMENSIONS)
    case nDims of
      2: begin ; 8-bit image
        dimensions = SIZE(img, /DIMENSIONS)
        image_is_8bit = 1b
      end

      3: begin ; RGB image
        allDims = SIZE(img, /DIMENSIONS)
        iInterleave = (WHERE(allDims eq 3))[0]
        case iInterleave of
          0: dimensions = allDims[1:2]
          1: dimensions = [allDims[0],allDims[2]]
          2: dimensions = allDims[0:1]
          else: begin
            void = DIALOG_MESSAGE(/INFORM, $
              ['Image must have dimensions', $
              '[3,m,n], [m,3,n], or [m,n,3].'])

          end
        endcase
        image_is_8bit = 0b
      end

      else: begin
        void = DIALOG_MESSAGE(/INFORM, $
          'Image data must have 2 or 3 dimensions.')
        nDims = 0
        img = 0
      end
    endcase
  endwhile

  if N_ELEMENTS(r) eq 0 then $
    r = BINDGEN(256)
  if N_ELEMENTS(g) eq 0 then $
    g = BINDGEN(256)
  if N_ELEMENTS(b) eq 0 then $
    b = BINDGEN(256)
  oPalette = OBJ_NEW('IDLgrPalette', r, g, b)

  pR = PTR_NEW(r)
  pG = PTR_NEW(g)
  pB = PTR_NEW(b)

  if N_ELEMENTS(roi_select_color) eq 0 then $
    roi_select_color = BYTE([255, 0, 0])

  if N_ELEMENTS(roi_color) eq 0 then $
    roi_color = BYTE([0, 255, 255])

  pSelRGB = PTR_NEW(roi_select_color)
  pROIRGB = PTR_NEW(roi_color)

;
;  ;define threshold value, GJ
;  initThresholdValue_min = -500
;  initThresholdValue_max = 1000

  oImage = OBJ_NEW('IDLgrImage', BYTSCL(img, MIN=threshold[0], MAX=threshold[1]), PALETTE=oPalette, $
    INTERLEAVE=iInterleave)
  oPalette1 = OBJ_NEW('IDLgrPalette')
  oPalette1->Loadct, 13

  Mask = DBLARR(imageSize*3, imageSize*3)
  IF STRCMP(direction, 'xy') THEN BEGIN
    Mask3D = REVERSE((CONGRID(new_mask_vol_cube[*,*,location3D[2]*(SIZE(new_mask_vol_cube))[3]/imageSize], $
      dimensions[0], dimensions[1],1) GT 0)*250, 2)
    Mask[*, *] = Mask3d[*, *, 0]
    oImage0 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage1 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage2 = OBJ_NEW('IDLgrImage', Mask*0., $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
  ENDIF

  IF STRCMP(direction, 'yz') THEN BEGIN
    Mask3D = (CONGRID(new_mask_vol_cube[location3D[0]*(SIZE(new_mask_vol_cube))[1]/imageSize,*,*], 1, $
      dimensions[0], dimensions[1]) GT 0)*140
    Mask[*, *] = Mask3d[0, *, *]
    oImage0 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage1 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage2 = OBJ_NEW('IDLgrImage', Mask*0., $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
  ENDIF

  IF STRCMP(direction, 'xz') THEN BEGIN
    Mask3D = (CONGRID(new_mask_vol_cube[*,location3D[1]*(SIZE(new_mask_vol_cube))[2]/imageSize,*], $
      dimensions[0], 1, dimensions[1]) GT 0)*210
    Mask[*, *] = Mask3d[*, 0, *]
    oImage0 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage1 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage2 = OBJ_NEW('IDLgrImage', Mask*0., $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1129
  IF STRCMP(direction, 'label') THEN BEGIN
    Mask3D = (CONGRID(new_mask_vol_cube[*,location3D[1]*(SIZE(new_mask_vol_cube))[2]/imageSize,*], $
      dimensions[0], 1, dimensions[1]) GT 0)*210
    Mask[*, *] = Mask3d[*, 0, *]
    oImage0 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage1 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage2 = OBJ_NEW('IDLgrImage', Mask*0., $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1129

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1212
  IF STRCMP(direction, 'distance') THEN BEGIN
    Mask3D = (CONGRID(new_mask_vol_cube[*,location3D[1]*(SIZE(new_mask_vol_cube))[2]/imageSize,*], $
      dimensions[0], 1, dimensions[1]) GT 0)*210
    Mask[*, *] = Mask3d[*, 0, *]
    oImage0 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage1 = OBJ_NEW('IDLgrImage', Mask, $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
    oImage2 = OBJ_NEW('IDLgrImage', Mask*0., $
      BLEND_FUNCTION=[3,4], ALPHA_CHANNEL=0.30, PALETTE=oPalette1)
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1212

  ;GJ 2019/4/27; hiding oImag2
  ;oImage2 is used to store ROI mask
  ;      oImage2->SetProperty, Hide = 1

  xdim = dimensions[0]
  ydim = dimensions[1]

  ; Widget UNAMES assigned below are used internally to help
  ; automate testing of xd_xroi.  A prefix is given to the UNAMES
  ; to help ensure that they are unique from UNAMES in other programs.
  ;
  prefix = 'xd_xroi:'

  ; Retrieve screen dimensions.
  DEVICE, GET_SCREEN_SIZE=screen
  xScrollSize = (N_ELEMENTS(xScrollSizeIn) gt 0) ? xScrollSizeIn : $
    (screen[0] - 200)
  xScrollSize = xScrollSize < xdim
  if ((xScrollSize + 20) gt xdim) then $
    xScrollSize = xdim
  yScrollSize = (N_ELEMENTS(yScrollSizeIn) gt 0) ? yScrollSizeIn : $
    (screen[1] - 250)
  yScrollSize = yScrollSize < ydim
  if ((yScrollSize + 20) gt ydim) then $
    yScrollSize = ydim

  ; Create top level base.
  if KEYWORD_SET(modal) then begin
    wBase = WIDGET_BASE( $
      TITLE=N_ELEMENTS(title) eq 0 ? 'ROI' : title, $
      GROUP_LEADER=group_leader, $
      bitmap=icon_dir+'BASE.bmp', $
      /COLUMN, $
      UNAME=prefix + 'xd_xroi', $
      FLOATING=floating, $
      /MODAL)
    wMenuBase = WIDGET_BASE(wBase, /ROW, /FRAME)
  endif else begin
    wBase = WIDGET_BASE( $
      TITLE=N_ELEMENTS(title) eq 0 ? 'ROI' : title, $
      GROUP_LEADER=group_leader, $
      bitmap=icon_dir+'BASE.bmp', $
      /COLUMN, $
      UNAME=prefix + 'xd_xroi', $
      UVALUE='BASE', $
      FLOATING=floating, $
      MBAR=wMenuBase, $
      /TLB_SIZE_EVENTS, $
      XOFFSET=(screen[0]/2 - xScrollSize/2) > 0, $
      YOFFSET=(screen[1]/2 - yScrollSize/2 - 100) > 0)
  endelse

  ;    ; Populate the menus.
  ;    wFileMenu = WIDGET_BUTTON(wMenuBase, VALUE='File', /MENU)
  ;    void = WIDGET_BUTTON( $
  ;        wFileMenu, $
  ;        VALUE='Import Image...', $
  ;        UNAME=prefix + 'import_image', $
  ;        UVALUE='IMPORT' $
  ;        )
  ;    wSaveButton = WIDGET_BUTTON( $
  ;        wFileMenu, $
  ;        VALUE='Save ROIs...', $
  ;        SENSITIVE=0, $
  ;        UNAME='save_menu_bttn', $
  ;        EVENT_FUNC='xd_xroi__Save' $
  ;        )
  ;    wButton = WIDGET_BUTTON( $
  ;        wFileMenu, $
  ;        VALUE='Quit', $
  ;        /SEPARATOR, $
  ;        UNAME=prefix + 'quit', $
  ;        EVENT_PRO='xd_xroi__Quit' $
  ;        )
  ;
  ;    wEditMenu = WIDGET_BUTTON(wMenuBase, VALUE='Edit', /MENU)
  ;    wButton = WIDGET_BUTTON( $
  ;        wEditMenu, $
  ;        VALUE='Copy Image', $
  ;        UNAME=prefix + 'copy', $
  ;        EVENT_FUNC='xd_xroi__Copy' $
  ;        )
  ;    wLoadCTButton = WIDGET_BUTTON( $
  ;        wEditMenu, $
  ;        VALUE='Image Color Table...', $
  ;        UVALUE='LOADCT', $
  ;        UNAME=prefix + 'ct', $
  ;        SENSITIVE=image_is_8bit, $
  ;        /SEPARATOR $
  ;        )
  ;    void = WIDGET_BUTTON( $
  ;        wEditMenu, $
  ;        VALUE='ROI Outline Colors...', $
  ;        UNAME=prefix + 'roi_color', $
  ;        UVALUE='PICKCOLOR' $
  ;        )
  ;    wButton = WIDGET_BUTTON( $
  ;        wEditMenu, $
  ;        VALUE='ROI Information...', $
  ;        UVALUE='ROI_INFO', $
  ;        UNAME=prefix + 'info_menu_bttn', $
  ;        /SEPARATOR $
  ;        )
  ;    wButton = WIDGET_BUTTON( $
  ;        wEditMenu, $
  ;        VALUE='Region Grow Properties...', $
  ;        UVALUE='ROI_GROW_PROPERTIES', $
  ;        UNAME=prefix + 'grow_menu_bttn' $
  ;        )
  ;
  ;    ; Create the toolbar.
  wToolbars = WIDGET_BASE(wBase, /ROW, /FRAME)

  subdir = ['resource','bitmaps']

  wToolbarBase = WIDGET_BASE(wToolbars, /ROW, SPACE=0, /TOOLBAR)
  ;    wSaveToolButton = WIDGET_BUTTON( $
  ;        wToolbarBase, $
  ;        VALUE=FILEPATH('save.bmp', SUBDIR=subdir), $
  ;        /BITMAP, $
  ;        TOOLTIP='Save all ROIs', $
  ;        UVALUE='SAVE', $
  ;        SENSITIVE=0, $
  ;        UNAME=prefix + 'save_tool_bttn', $
  ;        EVENT_FUNC='xd_xroi__Save' $
  ;        )
  ;    void = WIDGET_BUTTON( $
  ;        wToolbarBase, $
  ;        VALUE=FILEPATH('prop.bmp', SUBDIR=subdir), $
  ;        /BITMAP, $
  ;        TOOLTIP='Open ROI Information window', $
  ;        UNAME=prefix + 'info_tool_bttn', $
  ;        UVALUE='ROI_INFO' $
  ;        )
  ;    void = WIDGET_LABEL(wToolbarBase, VALUE=' ')
  ;    void = WIDGET_BUTTON( $
  ;        wToolbarBase, $
  ;        VALUE=FILEPATH('copy.bmp', SUBDIR=subdir), $
  ;        /BITMAP, $
  ;        TOOLTIP='Copy Image', $
  ;        UVALUE='COPY', $
  ;        UNAME=prefix + 'copy_tool_bttn', $
  ;        EVENT_FUNC='xd_xroi__Copy' $
  ;        )
  ;    void = WIDGET_BUTTON( $
  ;        wToolbarBase, $
  ;        VALUE=FILEPATH('flipvert.bmp', SUBDIR=subdir), $
  ;        /BITMAP, $
  ;        TOOLTIP='Flip image', $
  ;        UNAME=prefix + 'flip', $
  ;        UVALUE='FLIP' $
  ;        )

  if N_ELEMENTS(_tools) gt 1 then begin

    ; Create the exclusive button toolbar.
    wExcToolbarBase = WIDGET_BASE(wToolbars, /ROW, $
      /EXCLUSIVE, SPACE=0, /TOOLBAR)

    for i=0,N_ELEMENTS(_tools)-1 do begin
      case STRUPCASE(_tools[i]) of
        'CIRCULARMASK DRAW': begin
          wCircularMaskPoly = WIDGET_BUTTON( $
            wExcToolbarBase, $
            VALUE=FILEPATH('circularmaskpoly.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Draw Circular Mask ROIs', $
            UNAME=prefix + 'circularmask_mode', $
            UVALUE='CIRCULARMASKPOLY' $
            )
        end

        'FREEHAND DRAW': begin
          wFreePoly = WIDGET_BUTTON( $
            wExcToolbarBase, $
            VALUE=FILEPATH('freehand.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Draw Freehand ROIs', $
            UNAME=prefix + 'freehand_mode', $
            UVALUE='FREEPOLY' $
            )
        end

        'RECTANGLE': begin
          wRectangle = WIDGET_BUTTON( $
            wExcToolbarBase, $
            VALUE=FILEPATH('rectangle.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Draw Rectangle ROIs', $
            UNAME=prefix + 'rectangle_mode', $
            UVALUE='RECTANGLE' $
            )
        end
        'ELLIPSE': begin
          wEllipse = WIDGET_BUTTON( $
            wExcToolbarBase, $
            VALUE=FILEPATH('ellipse.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Draw Ellipse ROIs', $
            UNAME=prefix + 'ellipse_mode', $
            UVALUE='ELLIPSE' $
            )
        end

        'POLYGON DRAW': begin
          wSegPoly = WIDGET_BUTTON($
            wExcToolbarBase, $
            VALUE=FILEPATH('polygon.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Draw Polygon ROIs', $
            UNAME=prefix + 'polygon_mode', $
            UVALUE='SEGPOLY' $
            )
        end
      endcase
    endfor


    WIDGET_CONTROL,  WIDGET_INFO(wExcToolbarBase, /CHILD), /SET_BUTTON

  endif

  wTopRowBase = widget_base(wBase,/row,/frame)

  wControlBase = widget_base(wTopRowBase, xpad=0, ypad=0, /column)

  wMaxMinBase = widget_base(wControlBase, xpad=0, ypad=0, /column, /frame)
  minValue = MIN(img)-1
  maxValue = MAX(img)+1
  wThresholdLabel_min = WIDGET_LABEL(wMaxMinBase, $
    VALUE='Threshold_min:')

;  print, 'minValue  = ', minValue
;  print, 'maxValue = ', maxValue
;  print, 'initThresholdValue_min = ', initThresholdValue_min
;  print, 'initThresholdValue_max = ', initThresholdValue_max
  print, 'min vol = ', MIN(svol_HU_cube)
  print, 'max vol = ', MAX(svol_HU_cube)
  print, 'threshold = ', threshold

  wThresholdSlider_min = WIDGET_SLIDER(wMaxMinBase, $
    MINIMUM=MIN(svol_HU_cube), /DRAG, $
    MAXIMUM=MAX(svol_HU_cube), VALUE=MAX([threshold[0],MIN(svol_HU_cube)]), $
    UVALUE='THRESHOLD_min', $
    UNAME='xd_xroi:threshold_min_slider')

  wThresholdLabel_max = WIDGET_LABEL(wMaxMinBase, $
    VALUE='Threshold_max:')

  wThresholdSlider_max = WIDGET_SLIDER(wMaxMinBase, $
    MINIMUM=MIN(svol_HU_cube), /DRAG, $
    MAXIMUM=MAX(svol_HU_cube), VALUE=MIN([threshold[1], MAX(svol_HU_cube)]), $
    UVALUE='THRESHOLD_max', $
    UNAME='xd_xroi:threshold_max_slider')

  wThresholdLabel_radii = WIDGET_LABEL(wMaxMinBase, $
    VALUE='Radii:')

  wThresholdSlider_radii = WIDGET_SLIDER(wMaxMinBase, $
    MINIMUM=1, /DRAG, $
    MAXIMUM=50, VALUE=20, $
    UVALUE='THRESHOLD_radii', $
    UNAME='xd_xroi:threshold_radii_slider')


  _tools1 = STRUPCASE(['ADDRED', 'ADDLAYER', 'REMOVELAYER']);, 'ADDGREEN'])
  if N_ELEMENTS(_tools1) gt 1 then begin

    ; Create the exclusive button toolbar.
    wModelModeRadio = WIDGET_BASE(wControlBase, /ALIGN_CENTER, /COLUMN, $
      /EXCLUSIVE, SPACE=0, /TOOLBAR)

    for i=0,N_ELEMENTS(_tools1)-1 do begin
      case STRUPCASE(_tools1[i]) of
        'ADDRED': begin
          wAddNewRadio = WIDGET_BUTTON( $
            wModelModeRadio, $
            VALUE=FILEPATH('addred.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Add new object', $
            UNAME='xd_xroi:addred', $
            UVALUE='ADDRED' $
            )
        end
        'ADDLAYER': begin
          wAddLayerRadio = WIDGET_BUTTON( $
            wModelModeRadio, $
            VALUE=FILEPATH('addlayer.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Add layer', $
            UNAME='xd_xroi:addlayer', $
            UVALUE='ADDLAYER' $
            )
        end
        'REMOVELAYER': begin
          wRemoveLayerRadio = WIDGET_BUTTON( $
            wModelModeRadio, $
            VALUE=FILEPATH('removelayer.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Remove layer', $
            UNAME='xd_xroi:removelayer', $
            UVALUE='REMOVELAYER' $
            )
        end
      endcase
    endfor

    WIDGET_CONTROL,  WIDGET_INFO(wModelModeRadio, /CHILD), /SET_BUTTON

    widget_control, wModelModeRadio, $
      set_uname='xd_xroi:modelmoderadio', sensitive = 1

    IF N_ELEMENTS(ObjectIndex_value) THEN BEGIN
      IF ObjectIndex_value EQ 0 THEN BEGIN
        widget_control, wAddLayerRadio, sensitive = 0
        widget_control, wRemoveLayerRadio, sensitive = 0
      ENDIF
    ENDIF
  endif


  _tools2 = STRUPCASE(['CONTINUOUS', 'DISCRETE']);, 'ADDGREEN'])
  if N_ELEMENTS(_tools2) gt 1 then begin

    ; Create the exclusive button toolbar.
    wContModeRadio = WIDGET_BASE(wControlBase, /ALIGN_CENTER, /COLUMN, $
      /EXCLUSIVE, SPACE=0, /TOOLBAR)

    for i=0,N_ELEMENTS(_tools2)-1 do begin
      case STRUPCASE(_tools2[i]) of
        'CONTINUOUS': begin
          wContRadio = WIDGET_BUTTON( $
            wContModeRadio, $
            VALUE=FILEPATH('continuous.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Continuous', $
            UNAME='xd_xroi:continuous', $
            UVALUE='CONTINUOUS' $
            )
        end
        'DISCRETE': begin
          wIndeRadio = WIDGET_BUTTON( $
            wContModeRadio, $
            VALUE=FILEPATH('discrete.bmp', ROOT_DIR=icon_dir), $
            /BITMAP, $
            TOOLTIP='Not continuous', $
            UNAME='xd_xroi:discrete', $
            UVALUE='DISCRETE' $
            )
        end
      endcase
    endfor

    WIDGET_CONTROL,  WIDGET_INFO(wContModeRadio, /CHILD), /SET_BUTTON
  endif

  wModifyControlBase = widget_base(wControlBase, $
    xpad=0, $
    ypad=0, $
    /column $
    )

  wDrawBaseL = widget_base(wTopRowBase,/COL, /frame)
  wDrawBase = widget_base(wDrawBaseL,/frame)
  ; Create draw area.
  wDraw = WIDGET_DRAW( $
    wDrawBase, $
    /APP_SCROLL, $
    X_SCROLL_SIZE=xScrollSize, $
    Y_SCROLL_SIZE=yScrollSize, $
    XSIZE=xdim, $
    YSIZE=ydim, $
    UVALUE='DRAW', $
    GRAPHICS_LEVEL=2, $
    RENDERER=renderer, $
    /BUTTON_EVENTS, $
    /EXPOSE_EVENTS, $
    /VIEWPORT_EVENTS, $
    UNAME=prefix + 'draw', $
    /MOTION_EVENTS $
    )

  ;GJ 2019/1/2, add + and - sign to adjust circular mask size
  wHotKeyReceptor = WIDGET_TEXT(wDrawBase, $
    /ALL_EVENTS, $
    UVALUE='HOTKEY', $
    UNAME=prefix + 'hotkey')

  ; Create a context menu.
  wROIContextMenu = WIDGET_BASE(wDraw, /CONTEXT_MENU)
  wButton = WIDGET_BUTTON(wROIContextMenu, VALUE='Delete', $
    UVALUE='ROI_DELETE')
  wButton = WIDGET_BUTTON(wROIContextMenu, VALUE='Plot Histogram', $
    UVALUE='ROI_HISTOGRAM', /SEPARATOR)
  wGrowROIMenu = WIDGET_BUTTON(wROIContextMenu, VALUE='Grow Region', $
    UVALUE='ROI_GROW', /MENU)
  wButton = WIDGET_BUTTON(wGrowROIMenu, VALUE='By threshold', $
    UVALUE='ROI_GROW_BY_THRESHOLD')
  wButton = WIDGET_BUTTON(wGrowROIMenu, VALUE='By std. dev. multiple', $
    UVALUE='ROI_GROW_BY_STDDEV')
  wButton = WIDGET_BUTTON(wGrowROIMenu, VALUE='Properties...', $
    UVALUE='ROI_GROW_PROPERTIES', /SEPARATOR)

  ; Create a status bar.
  case _tools[0] of
    'TRANSLATE-SCALE': begin
      ; Determine if the user provided regions.
      haveRegions = 0
      if (N_ELEMENTS(regions_in) gt 0) then begin
        if SIZE(regions_in, /TNAME) eq 'OBJREF' then begin
          if OBJ_VALID(regions_in[0]) then $
            haveRegions = 1
        endif
      endif

      if (haveRegions ne 0) then begin
        ; If user provided regions, explain how to translate/scale.
        value = $
          'Click left mouse to select an ROI; drag to translate/scale.'
      endif else begin
        ; Otherwise, advise selection of a drawing tool.
        value = $
          'Choose a tool to draw an ROI.'
      endelse
    end
    'RECTANGLE': value = $
      'Click and drag left mouse to draw a rectangle ROI.'
    'ELLIPSE': value = $
      'Click and drag left mouse to draw an ellipse ROI.'
    'FREEHAND DRAW': value = $
      'Click and drag left mouse to draw freehand ROI.'
    'CIRCULARMASK DRAW': value = $
      'Click and drag left mouse to draw circular mask ROI.'
    'POLYGON DRAW': value = $
      'Click mouse to draw segmented ROI; double click to finish.'
    'SELECTION': value = $
      'Click left mouse to pick region & mark nearest vertex.'
  endcase

  IF STRCMP(direction, 'xy') THEN BEGIN
    wSlider = WIDGET_SLIDER(wDrawBaseL, $
      MINIMUM=0, /DRAG, $
      MAXIMUM=imageSize, VALUE= location3D[2])
  ENDIF

  IF STRCMP(direction, 'yz') THEN BEGIN
    wSlider = WIDGET_SLIDER(wDrawBaseL, $
      MINIMUM=0, /DRAG, $
      MAXIMUM=imageSize, VALUE= location3D[0])
  ENDIF

  IF STRCMP(direction, 'xz') THEN BEGIN
    wSlider = WIDGET_SLIDER(wDrawBaseL, $
      MINIMUM=0, /DRAG, $
      MAXIMUM=imageSize, VALUE= location3D[1])
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;hjx1129
  IF STRCMP(direction, 'label') THEN BEGIN
    wSlider = WIDGET_SLIDER(wDrawBaseL, $
      MINIMUM=0, /DRAG, $
      MAXIMUM=imageSize, VALUE= location3D[2])
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;hjx1129

  ;;;;;;;;;;;;;;;;;;;;hjx1212
  IF STRCMP(direction, 'distance') THEN BEGIN
    wSlider = WIDGET_SLIDER(wDrawBaseL, $
      MINIMUM=0, /DRAG, $
      MAXIMUM=imageSize, VALUE= location3D[2])
  ENDIF
  ;;;;;;;;;;;;;;;;;;;;hjx1212
  WIDGET_CONTROL, wSlider, SET_UVALUE = 'SLIDER', SET_UNAME='xd_xroi:Slider'

  wStatus = WIDGET_LABEL( $
    wBase, $
    VALUE=value, $
    /DYNAMIC_RESIZE, $
    /ALIGN_LEFT $
    )

  WIDGET_CONTROL, wBase, /REALIZE
  WIDGET_CONTROL, wDraw, GET_VALUE=oWindow

  ; Create the graphics hierarchy.
  oView = OBJ_NEW('IDLgrView', VIEWPLANE_RECT=[0,0,xdim,ydim])
  oModel = OBJ_NEW('IDLgrModel')
  oView->Add, oModel
  oModel->Add, oImage
  oModel->Add, oImage1
  ;    oModel->Add, oImage2
  ; Add a container for ROIs.
  oROIModel = OBJ_NEW('IDLgrModel')
  oROIGroup = OBJ_NEW('IDLanROIGroup')
  oModel->Add, oROIModel


  ;    oContour = mask_display(new_mask_vol_cube, location3D, imageSize, direction)
  ;    oModel->Add, oContour

  ; Create a selection visual for vertex picking.
  oPickVisual = xd_xroi__CreatePickVisual()
  oModel->Add, oPickVisual

  ; Create a translate/scale model to be used as a selection visual.
  oTransScaleVisual = xd_xroi__CreateTransScaleVisual(COLOR=roi_select_color)
  oModel->Add, oTransScaleVisual

  oRejected = OBJ_NEW('IDL_Container')
  oRegionsOut = OBJ_NEW('IDL_Container')
  oRegionsIn = OBJ_NEW('IDLgrModel')

  ; The following "if" conditions used on REGIONS_IN are intended
  ; to casue xd_xroi to silently ignore REGIONS_IN when REGIONS_IN=-1
  ; or REGIONS_IN=OBJ_NEW().  This makes xd_xroi easier to use because
  ; the user can assign the result of a container ::Get directly
  ; to keyword REGIONS_IN without worrying about the case where
  ; ::Get returns a -1, and because the user can assign the result
  ; of a previous REGIONS_OUT to REGIONS_IN without worrying about
  ; the case where the previous REGIONS_OUT was OBJ_NEW().
  if N_ELEMENTS(regions_in) gt 0 then begin
    if SIZE(regions_in, /TNAME) eq 'OBJREF' then begin
      if OBJ_VALID(regions_in[0]) then begin
        oRegionsOut->Add, regions_in
        oROIModel->Add, regions_in
        oROIGroup->Add, regions_in
        oRegionsIn->Add, regions_in, /ALIAS
        ;                if not lmgr(/demo) then begin
        ;                    WIDGET_CONTROL, wSaveButton, SENSITIVE=1
        ;                    WIDGET_CONTROL, wSaveToolButton, SENSITIVE=1
        ;                endif
      endif
    endif
  endif

  toolbar_geom = WIDGET_INFO(wToolbarBase, /GEOMETRY)

  ; Start with width of main toolbar.
  toolbar_xsize = toolbar_geom.scr_xsize

  ; Add on width of extra tools toolbar.
  if (N_ELEMENTS(wExcToolbarBase) gt 0) then $
    toolbar_xsize = toolbar_xsize + $
    (WIDGET_INFO(wExcToolbarBase, /GEOM)).scr_xsize

  base_geom = widget_info(wBase, /GEOMETRY)

  pImg = PTR_NEW(img, /NO_COPY)

  ;    if lmgr(/demo) then begin
  ;        WIDGET_CONTROL, wSaveToolButton, SENSITIVE=0
  ;        WIDGET_CONTROL, wSaveButton, SENSITIVE=0
  ;    endif

  sState = {wBase: wBase, $
    toolbar_xsize: toolbar_xsize, $
    toolbar_ysize: toolbar_geom.scr_ysize, $
    scr_xsize:base_geom.scr_xsize, $
    scr_ysize:base_geom.scr_ysize, $
    wDraw: wDraw, $
    wHotKeyReceptor: wHotKeyReceptor, $;GJ, 2019/1/2
    wROIContextMenu: wROIContextMenu, $
    wROIInfo: -1L, $
    wROIGrowProps: -1L, $
    ;              wSaveButton: wSaveButton, $
    ;              wSaveToolButton: wSaveToolButton, $
    image_is_8bit: image_is_8bit, $
    wStatus: wStatus, $
    oWindow: oWindow, $
    oView: oView, $
    oImage: oImage, $
    oImage0: oImage0, $ ;GJ 2019/1/2 original mask
    oImage1: oImage1, $; GJ, new mask
    oImage2: oImage2, $; GJ, new mask 2019/4/27

    ;             xd_sState:xd_sState, $    ;GJ, 2018/5/17
    ;              vol_HU_cube: vol_HU_cube, $
    svol_HU_cube: svol_HU_cube, $
    new_mask_vol_cube: new_mask_vol_cube, $
    original_mask_vol_cube: new_mask_vol_cube, $
    location3D: location3D, $
    imageSize: imageSize, $
    wSlider:wSlider, $
    direction:direction, $
    ;              oContour:oContour, $
    n_layers:0, $;;xd_sState.imageSize, $
    Modify_ROI:'ADDNEW1', $;'ADD' by default, $ ;GJ, 2018/5/22  BORDER, ADD, REMOVE
    wModelModeRadio:wModelModeRadio, $;GJ 2019/1/2
    wContModeRadio:wContModeRadio, $;GJ 2019/4/21
    CONTORINDE:'CONTINUOUS', $
    ;              wModifyModeRadio:wModifyModeRadio, $
    ;              wLayerModeRadio:wLayerModeRadio, $
    wThresholdSlider_min: wThresholdSlider_min, $
    wThresholdSlider_max: wThresholdSlider_max, $
    wThresholdSlider_radii: wThresholdSlider_radii, $
    Radii:20, $
    oModel: oModel, $
    oROIModel: oROIModel, $
    oROIGroup: oROIGroup, $
    oRegionsOut: oRegionsOut, $
    oRegionsIn: oRegionsIn, $
    oRejected: oRejected, $
    oSelROI: OBJ_NEW(), $
    oPickVisual: oPickVisual, $
    oTransScaleVisual: oTransScaleVisual, $
    oSelVisual: OBJ_NEW(), $
    oSelHandle: OBJ_NEW(), $
    pSavedROIData: PTR_NEW(), $
    savedROIXRange: DBLARR(2), $
    savedROIYRange: DBLARR(2), $
    oCurrROI: OBJ_NEW(), $
    oPalette: oPalette, $
    pR: pR, $
    pG: pG, $
    pB: pB, $
    pSelRGB: pSelRGB, $
    pROIRGB: pROIRGB, $
    sel_rgb: roi_select_color, $
    roi_rgb: roi_color, $
    growAllNeighbors: 0, $
    growThreshMin: 0, $
    growThreshMax: 255, $
    growThreshROI: 1, $
    growStdDevMult: 1.0, $
    growMinArea: 0.0, $
    growMaxCount: 2, $
    growAcceptAll: 0, $
    growRGBUseLuminosity: 1, $
    growRGBChannel: 0, $
    mode: _tools[0], $
    bFirstROI: 1B, $
    bButtonDown: 0B, $
    bTempSegment: 0B, $
    buttonXY: LONARR(2), $
    floating: keyword_set(floating), $
    pTitle: PTR_NEW(title), $
    pTools: PTR_NEW(), $
    wLoadCT: -1L, $
    ;              wLoadCTButton: wLoadCTButton, $
    wPaletteEdit: -1L, $
    wPickColor: -1l, $
    pImg: pImg, $
    modal: KEYWORD_SET(modal), $
    demo_system: KEYWORD_SET(demo_system), $
    group_leader: $
    N_ELEMENTS(group_leader) gt 0 ? group_leader : -1, $
    draw_time: 0d, $
    debug: KEYWORD_SET(debug) $
  }

  pState = PTR_NEW(sState, /NO_COPY)
  WIDGET_CONTROL, wBase, SET_UVALUE=pState

  xd_xroi__Viewport, {ID: wDraw, TOP: wBase, HANDLER: wBase, TYPE: 3, $
    X: 0, Y: (ydim - yScrollSize) > 0}

  XMANAGER, $
    'xd_xroi', $
    wBase, $
    NO_BLOCK=KEYWORD_SET(block) EQ 0, $
    CLEANUP='xd_xroi__cleanup'

  ;GJ 2019/1/2, add hotkey events
  ;    WIDGET_CONTROL, wHotKeyReceptor, /INPUT_FOCUS


  if KEYWORD_SET(modal) then begin
    PTR_FREE, pState  ; now we can free the state.
  endif

  ;
  ;   Return R, G & B parameters.
  ;
  r = *pR
  g = *pG
  b = *pB
  ;
  ;   Return ROI outline color.
  ;
  roi_select_color = *pSelRGB
  roi_color = *pROIRGB
  ;
  ;   Return REGIONS_OUT parameter.
  ;
  regions_out = oRegionsOut->Get(/ALL, COUNT=count)
  if count eq 0 then begin
    regions_out = OBJ_NEW()
  end

  ;
  ;   Return REJECTED parameter.
  ;
  rejected = oRejected->Get(/ALL, COUNT=rej_count)
  if rej_count eq 0 then begin
    rejected = OBJ_NEW()
  end

  ;
  ;   Return ROI_GEOMETRY and STATISTICS parameters.
  ;
  if N_ELEMENTS(statistics) gt 0 then $
    PTR_FREE, PTR_NEW(TEMPORARY(statistics))
  if N_ELEMENTS(roi_geometry) gt 0 then $
    PTR_FREE, PTR_NEW(TEMPORARY(roi_geometry))
  if count gt 0 then begin
    roi_geometry = REPLICATE({ $
      area: 0.0, $
      centroid: FLTARR(3), $
      perimeter: 0.0 $
    }, count)
    if image_is_8bit then $
      statistics = REPLICATE({ $
      count: 0UL, $
      minimum: 0.0, $
      maximum: 0.0, $
      mean: 0.0, $
      stddev: 0.0 $
    }, count)

    for i=0,count-1 do begin
      if ARG_PRESENT(roi_geometry) then begin
        if not regions_out[i]->ComputeGeometry( $
          AREA=area, $
          CENTROID=centroid, $
          PERIMETER=perimeter $
          ) $
          then $
          MESSAGE, 'failed to compute ROI_GEOMETRY.'
        roi_geometry[i].area = area
        roi_geometry[i].centroid = centroid
        roi_geometry[i].perimeter = perimeter
      endif
      if arg_present(statistics) and image_is_8bit then begin
        mask = regions_out[i]->ComputeMask( $
          DIMENSIONS=dimensions, $
          MASK_RULE=2 $
          )
        if N_ELEMENTS(mask) le 0 then $
          MESSAGE, 'failed to compute STATISTICS.' $
        else begin
          IMAGE_STATISTICS, $
            *pImg, $
            MASK=mask, $
            COUNT=pxl_count, $
            MEAN=mean, $
            STDDEV=stddev, $
            MINIMUM=minimum, $
            MAXIMUM=maximum
          statistics[i].count = pxl_count
          statistics[i].minimum = minimum
          statistics[i].maximum = maximum
          statistics[i].mean = mean
          statistics[i].stddev = stddev
        endelse
      endif
    endfor
  endif

  if ARG_PRESENT(regions_out) then begin
    oRegionsOut->Remove, /ALL
  endif else begin
    if oRegionsIn->Count() gt 0 then begin
      oRegionsOut->Remove, oRegionsIn->Get(/ALL)
    endif
  endelse
  oRejected->Remove, /ALL

  OBJ_DESTROY, oRejected
  OBJ_DESTROY, oRegionsOut
  OBJ_DESTROY, oRegionsIn

  PTR_FREE, pR
  PTR_FREE, pG
  PTR_FREE, pB
  PTR_FREE, pSelRGB
  PTR_FREE, pROIRGB

  ;
  ;   Deterimine if our call to XMANAGER (above) stopped program control
  ;   flow.  To do this, we cannot simply test KEYWORD_SET(BLOCK) OR
  ;   KEYWORD_SET(MODAL) because  BLOCK does not stop control flow if
  ;   control flow is already stopped somewhere else.  Thus we test
  ;   whether a variable such as oImage has been cleaned up.  If it has,
  ;   then we know our call to XMANAGER stopped control flow.
  ;
  flow_was_stopped = OBJ_VALID(oImage) eq 0
  ;
  if ARG_PRESENT(img) then begin
    ;
    ;       Return IMG parameter.  One reason we do this here is that IMG
    ;       may have become initialized or changed if the user had opportunity
    ;       to utilize our "Import Image..." menu button before we got here.
    ;
    if (not image_is_8bit) or flow_was_stopped then begin
      img = TEMPORARY(*pImg) ; We wont be needing *pImg any longer.
    endif else begin
      img = *pImg
    end
  endif else begin
    if not image_is_8bit then begin
      PTR_FREE, PTR_NEW(TEMPORARY(*pImg))
    endif
  endelse
  ;
  ;   In order that *pImg be available when we need it, we want to clean up
  ;   pImg at the last opportunity.  If our call to XMANAGER (above) stopped
  ;   program control flow, then by the time we get here, this is that last
  ;   opportunity.
  ;
  if flow_was_stopped then begin
    PTR_FREE, pImg
  end

end
