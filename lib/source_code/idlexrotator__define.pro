;$Id: //depot/Release/ENVI51_IDL83/idl/idldir/examples/doc/utilities/idlexrotator__define.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved.
;
;+
; NAME:
;   IDLexRotator
;
; PURPOSE:
;   Provide a type of IDLgrModel that can be rotated via mouse.
;       Also, provide a way to optionally constrain rotations 
;       with respect to self's orientation (CONSTRAIN=2).
;
; CATEGORY:
;   IDL object examples.
;
; CREATION:
;   oRotator = OBJ_NEW('IDLexRotator', Center, Radius)
;
; METHODS:
;   The methods for IDLexRotator are the same as those of IDL's
;   Trackball and IDLgrModel classes, with the following exceptions:
;       INIT:
;           Requires Center and Radius arguments.
;           Keyword TRANSLATE not accepted.
;           Keyword CONSTRAIN:  If this keyword is set to 2, constrain
;               with respect to self's orientation (rather than screen
;               orientation).
;
;               Note: Rotations specified via the ROTATE method are
;               not subject to constraints.
;
;       RESET:
;           Keyword TRANSLATE not accepted.
;           Keyword PRESERVE_ROTATION added.  When this keyword
;               is set, RESET will not change self's current TRANSFORM
;               property.
;
;       SCALE:
;           Not accepted.  Error thrown at run-time.
;
;       TRANSLATE:
;           Not accepted.  Error thrown at run time.
;
;       UPDATE:
;           Keyword TRANSLATE not accepted.
;           Keyword TRANSFORM not accepted.
;
;       SETPROPERTY:
;           Keyword TRANSFORM: Matrix must be a valid 4x4 rotation matrix,
;           Keyword CENTER added.  Sets Trackball center.
;           Keyword RADIUS added.  Sets Trackball radius.
;           Keyword AXIS added.  Sets Trackball axis for constraints.
;               0=constrain to X, 1=y, 2=z
;           Keyword CONSTRAIN added.  Sets Trackball rotation constraint.
;               0=none. 2=constrain to AXIS in Model orientation.
;               else: constrain to AXIS in screen orientation.
;           Keyword MOUSE added.  Sets Trackball MOUSE property.
;
;       GETPROPERTY:
;           Keyword CENTER added.  Gets Trackball center.
;           Keyword RADIUS added.  Gets Trackball radius.
;           Keyword AXIS added.  Gets Trackball axis for constraints.
;               0=constrain to X, 1=y, 2=z
;           Keyword CONSTRAIN added.  Gets Trackball rotation constraint.
;               0=none. 2=constrain to AXIS in Model orientation.
;               else: constrain to AXIS in screen orientation.
;           Keyword MOUSE added.  Gets Trackball MOUSE property.
;
; EXAMPLE:
;   pro example_event, event
;   widget_control, event.top, get_uvalue=pState
;   if (*pState).oRotator->Update(event) then $
;       (*pState).oWindow->Draw, (*pState).oView
;   end
;   ;
;   pro example_cleanup, wid
;   widget_control, wid, get_uvalue=pState
;   obj_destroy, (*pState).oWindow
;   obj_destroy, (*pState).oView
;   ptr_free, pState
;   end
;   ;
;   pro example
;   
;   tlb = widget_base()
;   
;   xsize = 300
;   ysize = 300
;   wDraw = widget_draw( $
;       tlb, $
;       xsize=xsize, $
;       ysize=ysize, $
;       graphics_level=2, $
;       /button_events, $
;       /motion_events $
;       )
;   widget_control, tlb, /realize
;   widget_control, wDraw, get_value=oWindow
;   
;   oModel = obj_new('IDLgrModel')
;   oModel->Add, obj_new('IDLgrSurface', dist(50))
;   oModel->Scale, 1./50, 1./50, 1./50
;   oModel->Translate, -.5, -.5, -.25
;   
;   oRotator = obj_new('IDLexRotator', $
;       [xsize/2.,ysize/2.], $
;       xsize/2. $
;       )
;   oRotator->Add, oModel
;   
;   oView = obj_new('IDLgrView')
;   oView->Add, oRotator
;   
;   oWindow->Draw, oView
;   
;   widget_control, tlb, set_uvalue=ptr_new({ $
;       oRotator: oRotator, $
;       oWindow: oWindow, $
;       oView: oView $
;       })
;   
;   xmanager, 'example', tlb, cleanup='example_cleanup', /no_block
;   end
;
;  MODIFICATION HISTORY: Written by PCS, RSI, 7/1998
;-
function IDLexRotator__matrix, q
;
;Given a unit quaternion, q, return its 4x4 rotation matrix.
;
x = q[0]
y = q[1]
z = q[2]
w = q[3]
a = transpose([[ w^2+x^2-y^2-z^2, 2*(x*y+w*z), 2*(x*z-w*y), 0], $
               [ 2*(x*y-w*z), w^2-x^2+y^2-z^2, 2*(y*z+w*x), 0], $
               [ 2*(x*z+w*y), 2*(y*z-w*x), w^2-x^2-y^2+z^2, 0], $
               [           0,           0,               0 ,1]])
return, a
end
;----------------------------------------------------------------------------
function IDLexRotator__valid, matrix
;Purpose: Check that MATRIX is a valid 4x4 rotation matrix (no translation
;  or perspective).
;
;  The logic here is borrowed from code in IDL 5.1 cw_arcball.pro.  That
;  routine used 1.0e-6 as the literal against which imperfection is
;  compared.  Here we use 1.0e-3 as a more relaxed test case, to allow
;  for larger imperfections which can accumulate after many rotations.
;
siz = size(matrix)
if siz[0] ne 2 or siz[1] ne 4 or siz[2] ne 4 then $
    return, 0 ; not a 4x4 matrix.

if (max(abs(matrix - matrix # transpose(matrix) # matrix)) gt 1.0e-3) then $
    return, 0 ; not a rotation matrix.

return, 1 ; Valid.
end
;----------------------------------------------------------------------------
pro IDLexRotator::cleanup
self->Trackball::cleanup
self->IDLgrModel::cleanup
end
;----------------------------------------------------------------------------
pro IDLexRotator::SetProperty, $
    transform=transform, $
    center=center, $
    radius=radius, $
    axis=axis, $
    constrain=constrain, $
    mouse=mouse, $
    _extra=e

if n_elements(transform) gt 0 then begin
    if keyword_set(self.check) then begin
        if not IDLexRotator__valid(transform) then begin
            message, 'Transform is not a valid rotation matrix.'
            end
        end

    self.anchor_transform=transform
    self->IDLgrModel::SetProperty, transform=transform
    end

if n_elements(center) gt 0 then $
    self->Reset, $
        center, $
        self.radius, $
        constrain=self.constrain, $
        axis=self.axis, $
        /preserve_rotation
if n_elements(radius) gt 0 then $
    self->Reset, $
        self.center, $
        radius, $
        constrain=self.constrain, $
        axis=self.axis, $
        /preserve_rotation
if n_elements(axis) gt 0 then $
    self->Reset, $
        self.center, $
        self.radius, $
        constrain=self.constrain, $
        axis=axis, $
        /preserve_rotation
if n_elements(constrain) gt 0 then $
    self->Reset, $
        self.center, $
        self.radius, $
        constrain=constrain, $
        axis=self.axis, $
        /preserve_rotation
if n_elements(mouse) gt 0 then $
    self->Reset, $
        self.center, $
        self.radius, $
        constrain=self.constrain, $
        axis=self.axis, $
        /preserve_rotation, $
        mouse=mouse

self->IDLgrModel::SetProperty, _extra=e
end
;----------------------------------------------------------------------------
pro IDLexRotator::GetProperty, $
    center=center, $
    radius=radius, $
    axis=axis, $
    constrain=constrain, $
    mouse=mouse, $
    _ref_extra=e

center = self.center
radius = self.radius
axis = self.axis
constrain = self.constrain
mouse = self.mouse

self->IDLgrModel::GetProperty, _extra=e
end
;----------------------------------------------------------------------------
;IDLexRotator::UPDATE
;
;Purpose:
;   Rotate self's TRANSFORM porperty via WIDGET_DRAW events.
;
function IDLexRotator::Update, $
    event, $
    mouse=mouse

;print, 'event.type = ', event.type
;print, 'event.press = ', event.press
;print, 'self.btndown = ', self.btndown
;print, 'self.mouse = ', self.mouse
;print, '+++++++++++++++++++'

if not keyword_set(self.debug) then $
    on_error, 2 ; Return to caller on error.

forward_function trackball_constrain

changed_transform = 0 ;Initialize return value.

mouse=4
if n_elements(mouse) ne 0 then begin
    if (mouse ne 1) and (mouse ne 2) and (mouse ne 4) then begin
        message, 'Invalid value for mouse keyword.'
        end $
    else $
        self.mouse = mouse
    end
;
;Ignore non-Draw-Widget events.
;
if (tag_names(event, /structure_name) ne 'WIDGET_DRAW') then $
    return, changed_transform
;
case event.type of
    0: begin ; Button press.
        if (event.press eq self.mouse) then begin
;
;           Set self.PT0: location of this event on unit sphere.
;
            xy = ([event.x,event.y] - self.center) / self.radius
            r = total(xy^2)
            if (r gt 1.0) then $
                self.pt0 = [xy/sqrt(r) ,0.0] $
            else $
                self.pt0 = [xy,sqrt(1.0-r)]
            if self.constrain ne 0 then begin
                if self.constrain eq 2 then $
                    constraint_vec = transpose( $
                        self.anchor_transform[0:2, 0:2] $
                        ) $
                else $
                    constraint_vec = identity(3)
                self.pt0 = trackball_constrain( $
                    self.pt0, $
                    constraint_vec[*, self.axis] $
                    )
                end
;
            self.btndown = 1b
            end
        end

    2: begin ; Button motion.
        if self.btndown eq 1b then begin
;
;           Set self.PT1: location of this event on unit sphere.
;
            xy = ([event.x,event.y] - self.center) / self.radius
            r = total(xy^2)
            if (r gt 1.0) then $
                self.pt1 = [xy/sqrt(r) ,0.0] $
            else $
                self.pt1 = [xy,sqrt(1.0-r)]
;            print, 'r = ', r
;            print, 'pt0 = ', self.pt0
;            print, 'pt1 = ', self.pt1
            IF sqrt(TOTAL(self.pt0)) EQ 0 THEN self.pt0 = self.pt1
            if self.constrain ne 0 then begin
                if self.constrain eq 2 then $
                    constraint_vec = transpose($
                        self.anchor_transform[0:2, 0:2] $
                        ) $
                else $
                    constraint_vec = identity(3)
                self.pt1 = trackball_constrain( $
                    self.pt1, $
                    constraint_vec[*, self.axis] $
                    )
                end
;
;           Update self's transform only if the mouse button has actually
;           moved from PT0.
;
            if ((self.pt0[0] ne self.pt1[0]) or $
                (self.pt0[1] ne self.pt1[1]) or $
                (self.pt0[2] ne self.pt1[2])) then begin

                changed_transform = 1b
;
;               Imperfections accumulate each time anchor_transform
;               is updated.  Thus we don't test the validity of TRANSFORM
;               because it will fail after sufficiently many rotations.
;
                self.check = 0b ; Don't test the validity of TRANSFORM.
;
                self->IDLgrModel::SetProperty, $
                    Transform= $
                        self.anchor_transform # $
                            IDLexRotator__matrix([ $
                                crossp(self.pt0, self.pt1), $
                                total(self.pt0 * self.pt1) $
                                ])
                self.check = 1b
                end
            end
        end

    1: begin ; Button release.
        self->IDLgrModel::GetProperty, transform=transform
        self.anchor_transform = transform
        self.btndown = 0b
        end

    else: $
        return, changed_transform

    endcase

return, changed_transform
end
;----------------------------------------------------------------------------
function IDLexRotator::init, center, radius, $
    debug=debug, $
    transform=transform, $ ; Initial rotation matrix.
    _extra=e

catch, error_status
if error_status ne 0 then begin
    catch, /cancel
    print, !error_state.msg
    return, 0
    end
if keyword_set(debug) then $
    catch, /cancel

if n_params() ne 2 then $
    message, 'IDLexRotator::init requires Center and Radius arguments.'

if n_elements(transform) gt 0 then begin
    self->SetProperty, transform=transform
    end $
else $
    self.anchor_transform = identity(4)

self.debug = keyword_set(debug)

if not self->Trackball::init(center, radius, _extra=e) then $
    message, 'IDLexRotator: Failed to init Trackball part of self.'

if not self->IDLgrModel::init(_extra=e) then $
    message, 'IDLexRotator: Failed to init IDLgrModel part of self.'

return, 1 ; Success.

end
;----------------------------------------------------------------------------
pro IDLexRotator::Reset, $
    center, $
    radius, $
    _extra=e, $
    preserve_rotation=preserve_rotation

self->Trackball::Reset, center, radius, _extra=e
if not keyword_set(preserve_rotation) then begin
    self->IDLgrModel::Reset
    self.anchor_transform = identity(4)
    end
end
;----------------------------------------------------------------------------
pro IDLexRotator::Scale
message, 'IDLexRotator does not scale.'
end
;----------------------------------------------------------------------------
pro IDLexRotator::Translate
message, 'IDLexRotator does not translate.'
end
;----------------------------------------------------------------------------
pro IDLexRotator::Rotate, axis, angle, _extra=e
self->IDLgrModel::Rotate, axis, angle, _extra=e
self->GetProperty, transform=transform
self.anchor_transform = transform
end
;----------------------------------------------------------------------------
pro IDLexRotator__define

void = {IDLexRotator, $
    inherits Trackball, $
    inherits IDLgrModel, $
    anchor_transform: fltarr(4,4), $
    check: 0b, $
    debug: 0b $
    }

end



