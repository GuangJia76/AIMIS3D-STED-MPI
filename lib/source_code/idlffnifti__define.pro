; ############################################################################
;  IDLffNifti__define
;
;  Defines the IDLffNifti object which can be used to read/write
;   Nifti1 format
;
;  IDLffNifti::isNum 
;  IDLffNifti::Swap
;  IDLffNifti::SwapHeader
;  IDLffNifti::IDLtoNiftiType
;  IDLffNifti::NiftitoIDLType
;  IDLffNifti::IDLtoNiftiBitPix
;  IDLffNifti::SystemByteOrder
;  IDLffNifti::Destroy
;  IDLffNifti::Calcsrow
;
;  Set
;  IDLffNifti::SetDim_info
;  IDLffNifti::SetDim
;  IDLffNifti::SetIntent_para
;  IDLffNifti::SetDatatype
;  IDLffNifti::SetBitpix
;  IDLffNifti::SetSlice_start
;  IDLffNifti::SetPixdim
;  IDLffNifti::SetVoxOffset
;  IDLffNifti::SetScl_slope
;  IDLffNifti::SetScl_inter
;  IDLffNifti::SetScl_end
;  IDLffNifti::SetScl_code
;  IDLffNifti::SetXYZT_units
;  IDLffNifti::SetCalMin
;  IDLffNifti::SetCalMax
;  IDLffNifti::SetGlmin
;  IDLffNifti::SetGlmax
;  IDLffNifti::SetXform
;  IDLffNifti::SetQuatern_bcd
;  IDLffNifti::SetOrigin
;  IDLffNifti::SetEndian
;
;  Header
;  IDLffNifti::GetProperty
;  IDLffNifti::SetProperty, header     
;
;  Data
;  IDLffNifti::QueryData
;  IDLffNifti::Getdata
;  IDLffNifti::Setdata
;
;  Read/Write
;  IDLffNifti::OutputHeader
;  IDLffNifti::ReadHeaderFile
;  IDLffNifti::ReadData
;  IDLffNifti::ReadFile
;
;  IDLffNifti::WriteFile
;
;  IDLffNifti::Init
;  IDLffNifti::Cleanup
;
;   class = { IDLffNifti,               $ 
;             endian  : 0,              $ ;  1=big-endian, 0=little-endian (read/write)
;             hdr     : {header_nifti}, $ ;  header data
;             default : {header_nifti}, $ ;  default header values
;             data    : ptr_new()       $ ;  actual image data
;           }
; ############################################################################



; ----------------------------------------------------------------------------
; IDLffNifti::isNum 
;
; Assert that n is a number and of the proper type
; ----------------------------------------------------------------------------
function IDLffNifti::isNum, n, POSITIVE=pos, NEGATIVE=neg, ZERO=zero,   $
                            LIMITS=lim, TYPE=type, INTEGER=integer, FLOATING=floating
  compile_opt idl2, HIDDEN
  
  ;;  Assert that n is a number and of the proper type
  t = size(n, /TYPE)
  
  if where(t eq [0,7,8,10,11]) ne -1 then return, 0
 
  if n_elements(type) ne 0 then begin
    if t ne type[0] then  return, 0
  endif else begin
    z = keyword_set(integer) + 2*keyword_set(floating)
    case z of
      0: ; no other type keywords, ignore
      1: if (where(t eq [1,2,3,12,13,14,15]) eq -1) then  return, 0  ; INTEGER
      2: if (t ne 4) and (t ne 5) then  return, 0                    ; FLOATING
      3: message, 'Incompatible keywords: /INTEGER and /FLOATING.'
    endcase
  endelse
  
  ;
  ;  Use limits or other keywords
  ;
  if ~keyword_set(lim) then begin  
    ;
    ;  Assert that n is positive, negative, or zero, etc.
    ;
    z = keyword_set(pos) + 2*keyword_set(neg) + 4*keyword_set(zero)

    case z of
      0: ; no keywords set, ignored
      1: if (min(n) le 0.0) then  return, 0         ; positive only
      2: if (max(n) ge 0.0) then  return, 0         ; negative only
      3: if (where(n eq 0.0) ne -1) then return, 0  ; pos or neg, anything but zero
      4: if (where(n ne 0.0) ne -1) then return, 0  ; zero only
      5: if (min(n) lt 0.0) then  return, 0         ; pos or zero
      6: if (max(n) gt 0.0) then  return, 0         ; neg or zero
      7: ; any value, ignore these keywords
    endcase
  
  endif else begin
    ;
    ;  Assert that all elements are in [lim[0], lim[1]]
    ;
    if (n_elements(lim) ne 2) then message, 'LIMITS must be a two-element numeric vector.'
    if (where(size(lim,/TYPE) eq [0,7,8,10,11]) ne -1) then  $
               message, 'LIMITS must be a numeric vector.'
      
    lo = lim[0]
    hi = lim[1]
    if (lo ge hi) then message, 'LIMITS must be of the form [lo,hi] with lo > hi.'
      
    if (min(n) lt lo) or (max(n) gt hi) then return, 0
  endelse
  
  ;  All assertions pass
  return, 1
end


;-----------------------------------------------------------
;  IDLffNifti::Swap
;
;  @param n {in}{type=numeric}  The value to be byte swapped.
;-----------------------------------------------------------
pro IDLffNifti::Swap, n
  compile_opt idl2, HIDDEN
  
  case size(n,/TYPE) of
    2: byteorder, n, /SSWAP      ; short
    3: byteorder, n, /LSWAP      ; long
    4: byteorder, n, /LSWAP      ; float
    5: byteorder, n, /L64SWAP    ; double
    6: byteorder, n, /L64SWAP    ; complex
  ELSE: ; nothing, leave as is
  endcase
end


;-----------------------------------------------------------
;  IDLffNifti::SwapHeader
;+
;  Byte swap a HEADER structure.
;
;  @param h {in}{type=HEADER structure}
;    The structure to byte swap.
;-----------------------------------------------------------
pro IDLffNifti::SwapHeader, h
  compile_opt idl2, HIDDEN
  
  ;;  Must be a header structure
  if (size(h,/TYPE) ne 8) then  message, 'Argument must be a HEADER structure.'
  if (size(h,/SNAME) ne 'HEADER_NIFTI') then message, 'Argument must be a HEADER structure.'
  h =  swap_endian(h)
end


;-----------------------------------------------------------
;  IDLffNifti::IDLtoNiftiType
;+
;  @returns The Nifti type for a given IDL type, -1 if
;    Nifti does not support the IDL type.
;
;  @param itype {in}{type=integer}
;    The IDL type code.
;-----------------------------------------------------------
function IDLffNifti::IDLtoNiftiType, itype
  compile_opt idl2, HIDDEN

  nifti_type = [2,4,8,16,32,64,128,256,512,768,1024,1280,1536,1792,2048,2304]
  idl_type   = [1,2,3, 4, 9, 5, -1, -1, 12, 13,  14,  15,  -1,  -1,  -1,  -1]

  idx  = where(itype eq idl_type)
  type = nifti_type[idx]
  if (type lt 2) then stop

  return, type
end


;-----------------------------------------------------------
;  IDLffNifti::NiftitoIDLType
;+
;  @returns The IDL type for a given Nifti type.
;
;  @param atype {in}{type=integer}
;    The Nifti type code.
;-----------------------------------------------------------
function IDLffNifti::NiftitoIDLType, atype
  compile_opt idl2, HIDDEN

  nifti_type = [2,4,8,16,32,64,128,256,512,768,1024,1280,1536,1792,2048,2304]
;  idl_type   = [1,2,3, 4, 9, 5, -1, -1, 12, 13,  14,  15,  -1,  -1,  -1,  -1]
  idl_type   = [1,2,3, 4, 9, 5,  1, -1, 12, 13,  14,  15,  -1,  -1,  -1,  -1]

  idx  = (where(atype eq nifti_type))[0]
  type = idl_type[idx]
  if (type lt 1) then stop

  return, type
end


;-----------------------------------------------------------
;  IDLffNifti::IDLtoNiftiBitPix
;+
;  @returns The Nifti bits per pixel for a given IDL type, -1 if
;    Nifti does not support the IDL type.
;
;  @param itype {in}{type=integer}
;    The IDL type code.
;-----------------------------------------------------------
function IDLffNifti::IDLtoNiftiBitPix, itype
  compile_opt idl2, HIDDEN

  idl_type = [1, 2, 3, 4, 9, 5, -1, -1, 12, 13, 14, 15,-1,-1,-1,-1]
  bit_pix  = [8,16,32,32,64,64,  0,  0, 16, 32, 64, 64, 0, 0, 0, 0] 

  idx = (where(itype eq idl_type))[0]
  if (idx lt 0) then stop
  bitpix = bit_pix[idx]

  return, bitpix
end


;-----------------------------------------------------------
;  IDLffNifti::SystemByteOrder
;+
;  @returns True (1) if the system is big-endian, false (0)
;    if the system is little-endian.
;-----------------------------------------------------------
function IDLffNifti::SystemByteOrder
  compile_opt HIDDEN

  b = (a=123)
  byteorder, a, /SWAP_IF_LITTLE_ENDIAN
  return, (b ne a) ? 0 : 1
end


;-----------------------------------------------------------
;  IDLffNifti::Destroy
;+
;  Free memory used by this object.
;-----------------------------------------------------------
pro IDLffNifti::Destroy
  compile_opt idl2, HIDDEN
  
  ;  remove any image data
  if ptr_valid(self.data) then  ptr_free, self.data
end




;-----------------------------------------------------------
;  IDLffNifti::Calcsrow
;
;-----------------------------------------------------------
PRO IDLffNifti::CalcSrow, quatern=quatern_bcd, trans=quatern_xyz
  compile_opt idl2, HIDDEN

 if (n_elements(quatern_bcd) eq 3) then begin
     self.hdr.quatern_b = float(quatern_bcd[0])
     self.hdr.quatern_c = float(quatern_bcd[1])
     self.hdr.quatern_d = float(quatern_bcd[2])
 endif 
 if (n_elements(quatern_xyz) eq 3) then begin
     self.hdr.quatern_x = float(quatern_xyz[0])
     self.hdr.quatern_y = float(quatern_xyz[1])
     self.hdr.quatern_z = float(quatern_xyz[2])
 endif 

 b    = self.hdr.quatern_b
 c    = self.hdr.quatern_c
 d    = self.hdr.quatern_d
 qfac = self.hdr.pixdim[0]      ; should be -1 or 1
 qx   = self.hdr.quatern_x
 qy   = self.hdr.quatern_y
 qz   = self.hdr.quatern_z


 ;; compute a parameter from b,c,d */
  a = 1.0 - (b*b + c*c + d*d)   
  if (a lt 1.e-7) then begin        ; special case
     a = 1d / sqrt(b*b+c*c+d*d) 
     b *= a
     c *= a
     d *= a 
     ;; normalize (b,c,d) vector
     a = 0d              ; a = 0 ==> 180 degree rotation
  endif else a = sqrt(a) ; angle = 2*arccos(a)

  ;; load rotation matrix, including scaling factors for voxel sizes
  i = self.hdr.pixdim[1]
  j = self.hdr.pixdim[2]
  k = qfac * self.hdr.pixdim[3]
  R = [[a*a+b*b-c*c-d*d,     2*b*c-2*a*d,       2*b*d+2*a*c],  $
       [2*b*c+2*a*d,         a*a+c*c-b*b-d*d,   2*c*d-2*a*b],  $
       [2*b*d-2*a*c,         2*c*d+2*a*b,       a*a+d*d-c*c-b*b]] 
  R = R ## diag_matrix([i,j,k]) 
  T = [self.hdr.quatern_x,self.hdr.quatern_y,self.hdr.quatern_z]

  self.hdr.srow_x = [R[*,0],T[0]]
  self.hdr.srow_y = [R[*,1],T[1]]
  self.hdr.srow_z = [R[*,2],T[2]]

;; print
;;   print, '----- CalcSrow -----'
;;   print, self.hdr.srow_x
;;   print, self.hdr.srow_y
;;   print, self.hdr.srow_z
;; print

end





;-----------------------------------------------------------
;  IDLffNifti::CalcQuatern
;
; Convert from rotation matrix to quaternion form
; calculate quaternion b,c,d  from the sform matrix
;-----------------------------------------------------------
PRO IDLffNifti::CalcQuatern, smatrix=smatrix
  compile_opt idl2, HIDDEN

  ;; smatrix is of type [3,4] as stored in hdr.srow_x,hdr.srow_y,hdr.srow_z

  sz = size(smatrix,/dim)
  if (n_elements(sz) eq 2) and (n_elements(smatrix) eq 12) then begin
     R   = smatrix[0:2,0:2] 
     qx  = smatrix[3,0]
     qy  = smatrix[3,1]
     qz  = smatrix[3,2]
  endif else begin
     R = dblarr(3,3)
     R[*,0] = self.hdr.srow_x[0:2]
     R[*,1] = self.hdr.srow_y[0:2]
     R[*,2] = self.hdr.srow_z[0:2]
     qx     = self.hdr.srow_x[3]
     qy     = self.hdr.srow_y[3]
     qz     = self.hdr.srow_z[3]
  endelse

  ;; compute lengths of each column => these determine grid spacings
  xd = norm(R[0:2,0])
  yd = norm(R[0:2,1])
  zd = norm(R[0:2,2])

  ;; if a column length is zero, patch the trouble
  if (xd eq 0) then begin
     R[0:2,0] = [1.0, 0.0, 0.0]
     xd       = 1.0
  endif
  if (yd eq 0) then begin
     R[0:2,1] = [0.0, 1.0, 0.0]
     yd       = 1.0
  endif
  if (zd eq 0 ) then begin
     R[0:2,2] = [0.0, 0.0, 1.0]
     zd       = 1.0
  endif
  dx=xd   &   dy=yd   &   dz=zd

  ;;  normalize the columns
  R[0:2,0] /= xd
  R[0:2,1] /= yd
  R[0:2,2] /= zd

  ;; At this point, the matrix has normal columns, but we have to allow
  ;; for the fact that the hideous user may not have given us a matrix
  ;; with orthogonal columns.
  ;; So, now find the orthogonal matrix closest to the current matrix.
  ;; One reason for using the polar decomposition to get this
  ;; orthogonal matrix, rather than just directly orthogonalizing
  ;; the columns, is so that inputting the inverse matrix to R
  ;; will result in the inverse orthogonal matrix at this point.
  ;; If we just orthogonalized the columns, this wouldn't necessarily hold. */
  ;; P = nifti_mat33_polar(Q) ;  
;@GJ, 2023/5/19, error may happen when saving nifti files
;  R = m3x3_ortho_polar(R)        ; R is now converted to the closest orthog matrix 
  zd = determ(R)                 ; should be -1 or 1 


  if (zd gt 0) then qfac = 1.0  $   ; proper
  else begin                        ; improper ==> flip 3rd column
      qfac = -1.0
      R[2,*] *= -1
   endelse

   ;; load 3x3 matrix into local variables
   r11 = R[0,0]   &   r12 = R[1,0]   &   r13 = R[2,0]
   r21 = R[0,1]   &   r22 = R[1,1]   &   r23 = R[2,1]
   r31 = R[0,2]   &   r32 = R[1,2]   &   r33 = R[2,2]

   ;; now, compute quaternion parameters
   a = total(diag_matrix(r)) + 1.0   

   if( a gt 0.5 ) then begin      ; simplest case
     a = 0.50 * sqrt(a) 
     b = 0.25 * (r32-r23) / a 
     c = 0.25 * (r13-r31) / a 
     d = 0.25 * (r21-r12) / a 
  endif else begin                ; trickier case
     xd = 1.0 + r11 - (r22+r33)   ; 4*b*b
     yd = 1.0 + r22 - (r11+r33)   ; 4*c*c
     zd = 1.0 + r33 - (r11+r22)   ; 4*d*d 
     if (xd gt 1.0) then begin
        b = 0.50 * sqrt(xd)      
        c = 0.25 * (r12+r21) / b 
        d = 0.25 * (r13+r31) / b 
        a = 0.25 * (r32-r23) / b 
     endif else if (yd gt 1.0) then begin
        c = 0.5 * sqrt(yd) 
        b = 0.25* (r12+r21) / c 
        d = 0.25* (r23+r32) / c 
        a = 0.25* (r13-r31) / c 
     endif else begin
        d = 0.5 * sqrt(zd)      
        b = 0.25* (r13+r31) / d 
        c = 0.25* (r23+r32) / d 
        a = 0.25* (r21-r12) / d 
     endelse
     if (a lt 0.0) then begin
        b = -b
        c = -c
        d = -d
        a = -a
     endif
  endelse

  ;; copy results in structure
  self.hdr.pixdim[0]   = qfac
  self.hdr.quatern_b   = b            ; quaterion b,c,d
  self.hdr.quatern_c   = c
  self.hdr.quatern_d   = d
  self.hdr.srow_x[0:2] = R[*,0] * dx  ; this maybe changed do to othogonalization
  self.hdr.srow_y[0:2] = R[*,1] * dy
  self.hdr.srow_z[0:2] = R[*,2] * dz

 
  ;; set origin
  self -> SetOrigin,[qx,qy,qz]

  ;; set pixdim
  pixdim      = self.hdr.pixdim
  pixdim[1:3] = [dx,dy,dz]
  self -> SetPixdim, pixdim
 
;; print
;;   print, '----- CalcQuatern -----'
;;   print, self.hdr.srow_x
;;   print, self.hdr.srow_y
;;   print, self.hdr.srow_z
;;   print
;;   print, self.hdr.quatern_b, self.hdr.quatern_c, self.hdr.quatern_d 

;; print

END








; =========================================================
;  Set functions
; =========================================================


;-----------------------------------------------------------
;  IDLffNifti::SetDim_info
;+
;  @param dim_info {in}{type=byte} number (MRI slice ordering
;-----------------------------------------------------------
pro IDLffNifti::SetDim_info, dim_info
  compile_opt idl2, HIDDEN
  
  ; Must be integer
  if (~self->isNum(dim, /BYTE)) then message, 'DIM must be of type BYTE.'
  if (n_elements(dim) ne 1) then message, 'DIM must be an 1 element number.'
    
  self.hdr.dim_info = dim_info
end





;-----------------------------------------------------------
;  IDLffNifti::SetDim
;+
;  @param dim {in}{type=vector} Dimension vector.
;-----------------------------------------------------------
pro IDLffNifti::SetDim, dim
  compile_opt idl2, HIDDEN
  
  ; Must be of type byte
  if (~self->isNum(dim, /INTEGER)) then message, 'DIM must be integer.'
  if (n_elements(dim) ne 8) then message, 'DIM must be an 8 element vector.'
    
  self.hdr.dim = dim
end


;-----------------------------------------------------------
;  IDLffNifti::SetIntent_para
;  
;  @param dim {in}{type = float 3D-vector} 3 intent_p parameter
;-----------------------------------------------------------
pro IDLffNifti::SetIntent_para, intent_para
  compile_opt idl2, HIDDEN

  ; Must be of type float
  if (~self->isNum(dim, /FLOAT)) then message, '3 intent paramters must be of type float.'
  if (n_elements(dim) ne 8) then message, 'Intent parameter must be a 3 element vector.'
    
  ; intent is a 3 element float vector:
  ; [intent_p1,intent_p2,intent_p3]
  self.hdr.intent_p1   = intent_para[0]
  self.hdr.intent_p2   = intent_para[1]
  self.hdr.intent_p3   = intent_para[2]
end



;-----------------------------------------------------------
;  IDLffNifti::SetDatatype
;+
;  N.B. This function expects the Nifti datatype codes, not
;       IDL datatype codes.  See self->DT_XXX() functions.
;
;  @param datatype {in}{type=int16} Data type code for image pixels.
;-----------------------------------------------------------
pro IDLffNifti::SetDatatype, datatype
  compile_opt idl2, HIDDEN
  
  ; Must be numeric and scalar
  if (~self->isNum(datatype, /POSITIVE)) then message, 'DATATYPE must be numeric.'
  if (n_elements(datatype) ne 1) then message, 'DATATYPE must be a scalar.'
    
  ; Must also be a valid code
  if (where(datatype eq [0,1,2,4,8,16,32,64,128,255]) eq -1) then $
            message, 'DATATYPE must be a valid Nifti datatype code.'

  self.hdr.datatype = datatype
end


;-----------------------------------------------------------
;  IDLffNifti::SetBitpix
;+
;  @param bitpix {in}{type=int16} Bits per pixel.
;-----------------------------------------------------------
pro IDLffNifti::SetBitpix, bitpix
  compile_opt idl2, HIDDEN
  
  ; Must be numeric and scalar
  if (~self->isNum(bitpix, /POSITIVE)) then message, 'BITPIX must be numeric.'
  if (n_elements(bitpix) ne 1) then message, 'BITPIX must be a scalar.'
    
  ; Must also be 1, 8, 16, 32, or 64
  if (where(bitpix eq [1,8,16,32,64]) eq -1) then $
            message, 'BITPIX must be 1, 8, 16, 32, or 64.'
    
  self.hdr.bitpix = bitpix
end


;-----------------------------------------------------------
;  IDLffNifti::SetSlice_start
;+
;  @param Slice_start {in}{type=Integer} first slice index
;-----------------------------------------------------------
pro IDLffNifti::SetSlice_start, value
  compile_opt idl2, HIDDEN
  
  if (~self->isNum(value, /INTEGER)) or (n_elements(value) ne 1) then $
        message, 'Slice_start must be a single number of type Integer.'
    
  self.hdr.slice_start = value
end


;-----------------------------------------------------------
;  IDLffNifti::SetPixdim
;+
;  @param pixdim {in}{type=vector} Pixel dimensions vector.
;-----------------------------------------------------------
pro IDLffNifti::SetPixdim, pixdim
  compile_opt idl2, HIDDEN
  
  ;  Must be floating point array
  if (~self->isNum(pixdim, /FLOATING)) then message, 'PIXDIM must be a 32-bit float.'
  if (n_elements(pixdim) ne 8) then message, 'PIXDIM must be an 8 element vector.'

  self.hdr.pixdim = pixdim
  self.hdr.srow_x[0:2] *= pixdim[1] 
  self.hdr.srow_y[0:2] *= pixdim[2] 
  self.hdr.srow_z[0:2] *= pixdim[3] 
end


;-----------------------------------------------------------
;  IDLffNifti::SetVoxOffset
;+
;  @param vox_offset {in}{type=float32}  Offset, in bytes, 
;                    to the start of voxel data in the file.
;-----------------------------------------------------------
pro IDLffNifti::SetVoxOffset, vox_offset
  compile_opt idl2, HIDDEN
  
  ;  Must be a positive or zero integer scalar
  if (~self->isNum(vox_offset, /POS, /ZERO)) then message, 'VOX_OFFSET must be a number >= 0.'
  if (n_elements(vox_offset) ne 1) then message, 'VOX_OFFSET must be a scalar.'
    
  self.hdr.vox_offset = vox_offset
end


;-----------------------------------------------------------
;  IDLffNifti::SetScl_slope
;-----------------------------------------------------------
pro IDLffNifti::SetScl_slope, value
  compile_opt idl2, HIDDEN

 if (~self->isNum(value, /Floating,/positive)) or (n_elements(value) ne 1) then $
        message, 'Scl_slope must be a single positive number of type Float.'

  self.hdr.scl_slope = value
end

;-----------------------------------------------------------
;  IDLffNifti::SetScl_inter
;-----------------------------------------------------------
pro IDLffNifti::SetScl_inter, value
  compile_opt idl2, HIDDEN

 if (~self->isNum(value, /Floating)) or (n_elements(value) ne 1) then $
        message, 'Scl_inter must be a single number of type Float.'

  self.hdr.scl_inter = inter
end

;-----------------------------------------------------------
;  IDLffNifti::SetScl_end
;-----------------------------------------------------------
pro IDLffNifti::SetScl_end, value
  compile_opt idl2, HIDDEN
 if (~self->isNum(value, /Integer,/positive)) or (n_elements(value) ne 1) then $
        message, 'Scl_end must be a single positive number of type INTEGER.'

  self.hdr.scl_end = value
end

;-----------------------------------------------------------
;  IDLffNifti::SetScl_code
;-----------------------------------------------------------
pro IDLffNifti::SetScl_code, value
  compile_opt idl2, HIDDEN

if (~self->isNum(value, /BYTE,/positive)) or (n_elements(value) ne 1) then $
        message, 'Scl_end must be a single positive number of type BYTE.'

  self.hdr.scl_code = value
end




;-----------------------------------------------------------
;  IDLffNifti::SetXYZT_units
;-----------------------------------------------------------
pro IDLffNifti::SetXYZT_units, value
  compile_opt idl2, HIDDEN

if (~self->isNum(value, /BYTE,/positive)) or (n_elements(value) ne 1) then $
        message, 'xyzt_units must be a single positive number of type BYTE.'

  self.hdr.xyzt_units = value
end



;-----------------------------------------------------------
;  IDLffNifti::SetCalMin
;+
;  @param cal_min {in}{type=float32} Minimum calibration value.
;-----------------------------------------------------------
pro IDLffNifti::SetCalMin, value
  compile_opt idl2, HIDDEN

  if (~self->isNum(value, /Floating)) or (n_elements(value) ne 1) then $
        message, 'Cal_min must be a single number of type Float.'

  self.hdr.cal_min = value
end


;-----------------------------------------------------------
;  IDLffNifti::SetCalMax
;+
;  @param cal_max {in}{type=float32} Maximum calibration value.
;-----------------------------------------------------------
pro IDLffNifti::SetCalMax, value
  compile_opt idl2, HIDDEN

  if (~self->isNum(value, /Floating)) or (n_elements(value) ne 1) then $
        message, 'Cal_max must be a single number of type Float.'

  self.hdr.cal_max = value
end


;-----------------------------------------------------------
;  IDLffNifti::SetGlmin
;+
;  Minimum pixel value for entire file.
;
;-----------------------------------------------------------
pro IDLffNifti::SetGlmin, value
  compile_opt idl2, HIDDEN
 
  if (~self->isNum(value, /LONG)) or (n_elements(value) ne 1) then $
        message, 'Glmin must be a single number of type Float.'

  self.hdr.glmin = value
end

;-----------------------------------------------------------
;  IDLffNifti::SetGlmax
;+
;  Maximum pixel value for entire file.
;
;-----------------------------------------------------------
pro IDLffNifti::SetGlmax, value
  compile_opt idl2, HIDDEN

  if (~self->isNum(value, /LONG)) or (n_elements(value) ne 1) then $
        message, 'Glmax must be a single number of type Float.'

  self.hdr.glmax = value
end


;-----------------------------------------------------------
;  IDLffNifti::SetXform
;-----------------------------------------------------------
pro IDLffNifti::SetXform, Xform
  compile_opt idl2, HIDDEN
  
 if (~self->isNum(xform, /INTEGER)) or (n_elements(value) ne 2) then $
        message, '[qform, sform] must be a two element vector of type INTEGER.'
  qform_code      = xform[0]
  sform_code      = xform[1]
  self.hdr.qform_code = qform_code
  self.hdr.sform_code = sform_code
end


;-----------------------------------------------------------
;  IDLffNifti::SetQuatern_bcd
;-----------------------------------------------------------
pro IDLffNifti::SetQuatern_bcd, quatern
  compile_opt idl2, HIDDEN

 if (~self->isNum(quatern, /FLOATING)) or (n_elements(quatern) ne 3) then $
        message, '[quatern_b,quatern_c,quatern_d] must be a 3 element vector of type FLOAT.'

  qb = quatern[0]
  qc = quatern[1]
  qd = quatern[2]
  self.hdr.quatern_b = qb
  self.hdr.quatern_c = qc
  self.hdr.quatern_d = qd
  if (self.hdr.qform_code eq 0) then self.hdr.qform_code = 1 
  self->Calcsrow
  if (self.hdr.sform_code eq 0) then self.hdr.sform_code = 1 
end



;-----------------------------------------------------------
;  IDLffNifti::SetOrigin
;-----------------------------------------------------------
pro IDLffNifti::SetOrigin, origin
  compile_opt idl2, HIDDEN

 if (~self->isNum(origin)) or (n_elements(origin) ne 3) then $
        message, '[quatern_z,quatern_y,quatern_z] must be a 3 element vector.'

  self.hdr.quatern_x = origin[0]
  self.hdr.quatern_y = origin[1]
  self.hdr.quatern_z = origin[2]
  self->Calcsrow
end




;-----------------------------------------------------------
;  IDLffNifti::SetEndian
;+
;  @param endian {in}{type=integer}  Endian, 1 (big), 0 (little).
;-----------------------------------------------------------
pro IDLffNifti::SetEndian, endian
  compile_opt idl2, HIDDEN
  
  ;  Must be 0 or 1
  if (~self->isNum(endian, /INTEGER, LIMITS=[0,1])) then  $
    message, 'ENDIAN must be 0 (little) or 1 (big).'
    
  self.endian = endian
end



;-----------------------------------------------------------
;  IDLffNifti::GetProperty
;
;-----------------------------------------------------------
FUNCTION IDLffNifti::GetProperty
  compile_opt idl2

  return, self.hdr
END


;-----------------------------------------------------------
;  IDLffNifti::SetProperty
;+
;  This is a dangerous function.  *No* checking whatsoever
;  is performed for consistent values.
;
;  @param all {in}{type=header struct} Structure containing
;             the entire Nifti header.
;-----------------------------------------------------------
pro IDLffNifti::SetProperty, hdr
  compile_opt idl2, HIDDEN

  ; Must be a header structure
  if (size(hdr,/SNAME) ne 'HEADER_NIFTI') then message, 'Must be a Nifti header structure.'
    
  self.hdr = hdr
end



;-----------------------------------------------------------
;  IDLffNifti::QueryData
;+
;  @returns True (1) if there is valid image data, false (0)
;       otherwise.
;-----------------------------------------------------------
function IDLffNifti::QueryData
  return, ptr_valid(self.data)
end



;-----------------------------------------------------------
;  IDLffNifti::Getdata
;+
;  @returns The pixel data, or query image data status.
;
;  @keyword POINTER {in}{type=boolean}{optional}
;    If set, return a pointer to the image data.  This
;    is against "good" object-oriented design but these
;    files might be very large and direct access through a 
;    pointer will not result in copying megabytes of data.
;-
function IDLffNifti::Getdata, POINTER=pointer
  compile_opt idl2
  
  if ~ptr_valid(self.data) then message, 'No image data available.'
      
  if ~keyword_set(pointer) then begin
      return, *self.data ; this might require a lot of memory
  endif else begin
      return, self.data  ; memory friendly, but dangerous
  endelse
end


;-----------------------------------------------------------
;  IDLffNifti::Setdata
;+
;  Set the pixel data and associated header elements using
;  the input array.  Will set glmin, glmax, datatype, bitpix,
;  and dim but does not set orient.
;
;  @param d {in}{type=numeric array, 2, 3, or 4D}
;    The data to store.
;
;  @keyword NO_COPY {in}{type=boolean}{optional}
;    If set the original array input will be destroyed
;    once the data has been assigned to the heap.  Ie,
;    the NO_COPY keyword to ptr_new().
;-----------------------------------------------------------
pro IDLffNifti::Setdata, d, NO_COPY=no_copy0
  compile_opt idl2
  
  no_copy = keyword_set(no_copy0) ? no_copy0 : 0
  
  ;
  ;  Must be a numeric array that is supported by the
  ;  Nifti1 "standard".
  ;
  if (where(size(d,/TYPE) eq [1,2,3,4,5,6,12]) eq -1) then  $
    message, 'Input array not of a type supported by Nifti-1.'
    
  ;  Must be 2, 3, or 4 dimensional
  dims = size(d,/DIM)
  ndims = n_elements(dims)
  if (ndims lt 2) or (ndims gt 4) then  $
    message, 'Input array must be of 2, 3 or 4 dimensions.'
 
  ;  Remove any existing data
  if ptr_valid(self.data) then ptr_free, self.data
  
  
  ;  Set the dimensions
  self.hdr.dim[0] = size(d,/n_dimensions)       ; number of dimensions
  self.hdr.dim[1] = dims[0]                     ; X direction (columns)
  self.hdr.dim[2] = dims[1]                     ; Y direction (rows)
  self.hdr.dim[3] = (ndims ge 3) ? dims[2] : 1  ; Z direction (slices), 2D is one slice
  self.hdr.dim[4] = (ndims eq 4) ? dims[3] : 1  ; volumes (time points), 3D is one time point
  
  ;  Set the datatype
  self.hdr.datatype = self->IDLtoNiftiType(size(d,/TYPE))
  
  ;  Set bits per pixel
  self.hdr.bitpix = self->IDLtoNiftiBitPix(size(d,/TYPE))
  
  ;  Set the glmin and glmax
  self.hdr.glmin = min(d)
  self.hdr.glmax = max(d)
  
  ;  Set endian to the order this system uses
  self.endian = self->SystemByteOrder()

  ;  Set the data itself
  self.data = ptr_new(d, NO_COPY=no_copy)
end



;
;  File I/O methods:
;

;----------------------------------------------------------
;  IDLffNifti::OutputHeader
;----------------------------------------------------------
PRO IDLffNifti::OutputHeader, hdr, flun

  if (size(hdr,/type) ne 8) then hdr = self.hdr 
  if (n_elements(flun) ne 1) then flun=-1 
  tag   = tag_names(hdr)
  for i=0,n_tags(hdr)-1 do begin
      type = size(hdr.(i),/type)
      case type of
          1: begin
              n = n_elements(hdr.(i))
              if (n gt 1) then begin
                  info = strarr(n)
                  for ni=0, n-1 do info[ni] = string((hdr.(i))[ni])
                  info = strjoin(info)
              endif else info = string(format='(I0)',hdr.(i))
          end
          2: info = strjoin(string(format='(100I)',hdr.(i))) 
          3: info = strjoin(string(format='(100I)',hdr.(i))) 
          4: info = strjoin(string(format='(100F)',hdr.(i)))
      endcase
      printf, flun, format='(A30,": ",A)',tag[i],strtrim(strcompress(info),2)
  endfor
END




;----------------------------------------------------------
;  IDLffNifti::ReadHeaderFile
;+
;  @returns  An instance of {header} for the given .hdr file.
;
;  @param f {in}{type=string}
;    Pathname to the header file.
;
;  @keyword ENDIAN {out}{type=integer}
;    Set to the byte order of the header (1=big, 0=little).
;-----------------------------------------------------------
function IDLffNifti::ReadHeaderFile, f, ENDIAN=endian, compress=compress
  compile_opt idl2
  
  ;; Filename must really exist with read permission
  if (~file_test(f,/READ)) then message, 'Unable to read file ['+f+']'

  ;; nii or hdr
  ext = strmid(f,strlen(f)-3,3)
  if (ext eq 'hdr') then begin
      ;; assume hdr type of file 
      ;;  File must be 352 bytes long, no extended headers
      info = file_info(f)
      if (info.size ne 348) then  $
          message, 'Non-standard header files not supported, length must be 348 bytes.'
  endif
 
  ;; Read data into structure
  hdr = {header_nifti}
  openr, u, f, /GET_LUN , compress=compress
  readu, u, hdr
  free_lun, u
  
  ;  Determine file byte order and set local ENDIAN property.
  ;  Used by ReadData() to swap bytes properly.
  if (hdr.sizeof_hdr eq 348) then begin
      ;; system and header using same byte order
      endian = self->SystemByteOrder()
  endif else begin
    ; system and header using opposite byte order
    endian = ~self->SystemByteOrder()  ; store true order
    self->SwapHeader, hdr  ; swap to get to local system order
  endelse

  self.hdr = hdr
  
  ;  Return the structure
  return, hdr
end



;----------------------------------------------------------
;  IDLffNifti::ReadData
;+
;  @returns The data for the given file.
;
;  @param f {in}{type=string}
;    The pathname of the image file.
;
;  @param h {in}{type=HEADER structure}
;    The header file structure associated with this image data.
;
;
;  @keyword ENDIAN {in}{type=integer}
;    Byte order of the image data.
;-----------------------------------------------------------
;@private
;-
function IDLffNifti::ReadData, f, h, ENDIAN=endian, compress=compress
  compile_opt idl2, HIDDEN
  
  ;; check nii / img
  name = strmid(f,0,strlen(f)-3)
  ext  = strmid(f,strlen(f)-3,3)
  if (ext eq '.gz') then begin
     name = strmid(f,0,strlen(f)-6)
     ext = strmid(f,strlen(f)-6,6)
  endif

  case ext of 
      'hdr':    fname = name + 'img'
      'nii':    fname = name + 'nii'
      'nii.gz': fname = name + 'nii.gz'
      else : message, 'Unsupported file format ['+ f +']'
  endcase

  ;  Verify that there is a data file to read.
  if (~file_test(fname,/READ)) then message, 'Unable to read file ['+fname+']'
 
  ;; read data
  if (h.dim[4] eq 0) then h.dim[4] = 1

  if (h.dim[0] eq 128) then $    ; RGB type 24 bit  
       data = make_array(h.dim[1], h.dim[2], h.dim[0]-1, h.dim[3], h.dim[4],TYPE=self->NiftitoIDLType(h.datatype)) $
  else $  ; non-RGB
       data = make_array(h.dim[1], h.dim[2], h.dim[3], h.dim[4],TYPE=self->NiftitoIDLType(h.datatype))


  openr,lun,fname,/get_lun, compress=compress, swap_endian=endian
  finfo =  Fstat(lun)
  point_lun,lun,h.vox_offset    ;; for nii: move file pointer to data   
  readu,lun,data
  free_lun,lun

  if (h.dim[0] eq 128) then data = transpose(data,[2,0,1])   ; transpose 24 bit RGB image 


  self.hdr  = h
  self.data =  ptr_new(data)

  return, data
end


;-----------------------------------------------------------
;  IDLffNifti::ReadFile
;+
;  Reads an Nifti-1 file.  
;
;  @param f {in}{type=string vector}{optional}
;    Name or names of the files to load.  If not given
;    an open dialog box is used.
;
;  @keyword ORIENTATION {in}{type=string or integer}{optional}
;    Orientation to use when reading the file, overrides any
;    specified in the header.  If not given use the value in
;    the header.  Integer 0..5 or 'TRS','COR','SAG' with '-'
;    in front for flipped (values 3..5).
;
;  @keyword _EXTRA {in}{type=any}{optional}
;    Extra keywords for dialog_pickfile.
;-----------------------------------------------------------
pro IDLffNifti::ReadFile, f, _EXTRA=extra, compress=compress
  compile_opt idl2

  ; If f given?
  if (n_params() eq 0) then begin
    f = dialog_pickfile(/MUST_EXIST, FILTER=['*.hdr','*.nii'], _EXTRA=extra)
    if f[0] eq '' then return  ; user cancel
    f = f[0]
  endif else begin
    ;  Validate f
    if (size(f,/TYPE) ne 7) then message, 'Argument must be a string.'
    if (n_elements(f) ne 1) then message, 'Argument must be a single string.'
  endelse
  
  ;  Read the header file
  h = self->ReadHeaderFile(f, ENDIAN=endian, compress=compress)

  ;  Now read the data
  data = self->ReadData(f, h, ENDIAN=endian, compress=compress)
  
  ;  All read successfully, set necessary fields
  if ptr_valid(self.data) then  ptr_free, self.data
  self.data   = ptr_new(data,/NO_COPY)
  self.hdr    = h
  self.endian = endian
end



;-----------------------------------------------------------
;  IDLffNifti::WriteFile
;+
;  Output the image data to a file or set of files.
;
;  @param f {in}{type=string}{optional}
;    The output file name.  Base name only, no extensions.
;
;  @keyword _EXTRA {in}{type=any}{optional}
;    Extra keywords for dialog_pickfile.
;-----------------------------------------------------------
 PRO IDLffNifti::WriteFile, fname, just_header=just_header, _EXTRA=extra
   compile_opt idl2

   ;; Must have something to write
   if not keyword_set(just_header) and ~ptr_valid(self.data) then message, 'No image data available.'
    
   ;;  Get a filename
   if (n_params() eq 1) then begin
       if (size(fname,/TYPE) ne 7) then fname = (dialog_pickfile(_EXTRA=extra))[0]
       if (fname eq '') then begin
           message, 'Output filename must be a string.'
           return
       end
   endif

   ;; check nii / img
   name = strmid(fname,0,strlen(fname)-3)
   ext  = strmid(fname,strlen(fname)-3,3)
   ;; store magic code 
   if (ext eq 'hdr') then $
        self.hdr.magic = byte([110,105,49,0]) $  ; for hdr: "ni1\0"
   else self.hdr.magic = byte([110, 43,49,0])    ; for nii: "n+1\0"

   ; write header first 
   openw,lun,fname,/get_lun
   writeu,lun,self.hdr

   ;; write volume data
   if not keyword_set(just_header) then begin

       ;; write img data
       if (ext eq 'hdr') then begin ; separate file for img
           free_lun,lun
           openw,lun,name+'img',/get_lun
       endif else begin
           ;; write nii data
           point_lun,lun,self.hdr.vox_offset   
       endelse 
       writeu,lun,*self.data
   endif
   free_lun,lun

 END



;-----------------------------------------------------------
;  IDLffNifti::AppendFile
;+
;  Append image data to one file
;
; fname :   filename of the Nifti file 
; volume:   data to append  (only 1 volume allowed)
;
; NOTE:  there is no check if data fit to the rest of the data
;-----------------------------------------------------------
 PRO IDLffNifti::AppendFile, fname, vol, yreverse=yreverse
   compile_opt idl2

   ;; check input
   if (file_test(fname) ne 1) then begin
       msg  = ['ERROR: IDLffNifti::AppendFile','']
       msg  = [msg,'Cannot find file:',fname,'']
       tmp  = dialog_message(msg,/ERROR)
       return
   endif
   ndim = size(vol,/n_dim)
   if (ndim lt 3) or (ndim gt 4) then begin
       msg  = ['ERROR: IDLffNifti::AppendFile','']
       msg  = [msg,'Size of volume(s) MUST be 3 or 4 dimensional:','']
       tmp  = dialog_message(msg,/ERROR)
       return
   endif else begin
       if (ndim eq 3) then nvol=1 else nvol=(size(vol,/dim))[3]
   endelse

   ;; read header data
   hdr = self->ReadHeaderFile(fname)


   ;; check nii / img
   name = strmid(fname,0,strlen(fname)-3)
   ext  = strmid(fname,strlen(fname)-3,3)
   ;; store magic code 
   if (ext eq 'hdr') then $
        fname_img = name + '.img' $
   else fname_img = fname 

   
   ;; check img file with image data exists
   if (file_test(fname_img) ne 1) then begin
       msg  = ['ERROR: IDLffNifti::AppendFile','']
       msg  = [msg,'Cannot find file:',fname_img,'']
       tmp  = dialog_message(msg,/ERROR)
       return
   endif


   ;; check if we have to flip top-bottom
   data = reverse(vol,2)


   ;; append data 
   openw,lun,fname_img,/get_lun,/append
   writeu,lun, data


   ;; modify header data 
   hdr.dim[4] += nvol
   case hdr.dim[0] of 
       3:  hdr.dim[0] = 4
       4:  
       else: stop
   endcase
   point_lun,lun,0
   writeu,lun,hdr
   free_lun,lun   


 END









;
;  Lifecycle methods:
;

;-----------------------------------------------------------
;  IDLffNifti::Init
;+
;  Create an instance of this class.
;
;  @returns  True (1) if the class successfully created.
;
;  @param f {in}{type=string}{optional}
;    If present, must be a string giving the name of a 
;    header file.  This file (and its associated .img
;    file which must be in the same directory) are loaded.  
;
;  @keyword DATA {in}{type=array}{optional}
;    If present, the data from this array is used to
;    set up the object.  This overrides any given argument.
;-----------------------------------------------------------
function IDLffNifti::Init, f, DATA=data
  compile_opt idl2

  ;;  Catch any errors so memory isn't wasted
  CATCH, err
  if (err ne 0) then begin
      catch, /CANCEL
      self->Destroy             ; release any allocated memory
      return, 0                 ; object creation failed
  endif

  ;;  Set certain constant values, save as default
  self.hdr.sizeof_hdr  = long(348)
  self.hdr.extents     = 0L
  self.hdr.regular     = byte('r')
  self.hdr.dim_info    = 0b
  self.hdr.intent_p1   = 0.0
  self.hdr.intent_p2   = 0.0
  self.hdr.intent_p3   = 0.0
  self.hdr.intent_code = 0s
  self.hdr.slice_start = 0s
  self.hdr.pixdim[0]   = -1
  self.hdr.pixdim[1:3] =  1
  self.hdr.vox_offset  = 0.0
  self.hdr.scl_slope   = 1.0
  self.hdr.scl_inter   = 0.0
  self.hdr.scl_end     = 0s
  self.hdr.scl_code    = 0s
  self.hdr.xyzt_units  = 10b
  self.hdr.cal_max     = 0.0
  self.hdr.cal_min     = 0.0
  self.hdr.slice_duration  = 0.0
  self.hdr.toffset     = 0.0
   self.hdr.glmax       = 0L
  self.hdr.glmin       = 0L
  self.hdr.qform_code  = 1s  
  self.hdr.sform_code  = 1s
  self.hdr.quatern_b   =  0.5                  ; default rotation
  self.hdr.quatern_c   = -0.5                  ; default rotation
  self.hdr.quatern_d   = -0.5                  ; default rotation
  self.hdr.quatern_x   = -88                   ; default origin
  self.hdr.quatern_y   =  128                  ; default origin
  self.hdr.quatern_z   = -128                  ; default origin
  self.hdr.magic       = byte([110,105,49,0])
  
  ; calculate srow
  self->CalcSrow

  ;; store defaults into default header structure
  self.default = self.hdr
  
  ;; Use the system byte order by default
  self.endian = self->SystemByteOrder()
  
  ;; read data if filename is given
  if (size(f,/type) eq 7) then self->ReadFile, f

  ;; or store data (even overwrite after reading) into structure
  if keyword_set(data) then self->Setdata, data 
  
  ;;  have a valid object
  return, 1
end


;-----------------------------------------------------------
;  IDLffNifti::Cleanup
;
;  Destory the object.
;-----------------------------------------------------------
pro IDLffNifti::Cleanup
  compile_opt idl2
  
  self->Destroy  ; separate destroy method so it can be
                 ; called from Init if necessary
end


;-----------------------------------------------------------
;  IDLffNifti__define
;+
;  Define the class and other structures.
;-----------------------------------------------------------
pro idlffnifti__define
  compile_opt idl2  ;  N.B. matters here!

  ;
  ;  These structures define the contents of the .hdr file:
  ;
  
  ;
  ;  Nifti header structure - incorporates the
  ;    header_key, image_dimension, and data_history structures.
  ;
  void  = {header_nifti,                     $                         
             sizeof_hdr    : 348L,           $ ; int sizeof_hdr                    offset:   0
             data_type     : bytarr(10),     $ ; char data_type[10]                offset:   4
             db_name       : bytarr(18),     $ ; char db_name[18]                  offset:  14
             extents       : 0L,             $ ; extents                           offset:  32
             session_error : 0s,             $ ; short int session_error           offset:  36
             regular       : 0b,             $ ; char regular;                     offset:  38
             dim_info      : 0b,             $ ; MRI slice ordering                offset:  39

             ;;  Nifti image_dimension structure:
             dim         : intarr(8),        $ ; data array dimension vector       offset:  40
             intent_p1   : 0.0,              $ ; 1st intent parameter              offset:  56
             intent_p2   : 0.0,              $ ; 2nd intent parameter              offset:  60
             intent_p3   : 0.0,              $ ; 3rd intent parameter              offset:  64
             intent_code : 0s,               $ ; Nifti intent code                 offset:  68
             datatype    : 0s,               $ ; short int datatype                offset:  70
             bitpix      : 0s,               $ ; number of bits per voxel          offset:  72
             slice_start : 0s,               $ ; first slice index                 offset:  74
             pixdim      : fltarr(8),        $ ; float pixdim[8]                   offset:  76
             vox_offset  : 0.0,              $ ; offset into *.nii data file       offset: 108     
             scl_slope   : 1.0,              $ ; data scaling slope                offset: 112
             scl_inter   : 0.0,              $ ; data scaling offset               offset: 116
             scl_end     : 0s,               $ ; last slice index                  offset: 120
             scl_code    : 0b,               $ ; slice timing order                offset: 122
             xyzt_units  : 0b,               $ ; units of pixdim [1..4]            offset: 123
             cal_max     : 0.0,              $ ; max display intensity             offset: 124
             cal_min     : 0.0,              $ ; min display intensity             offset: 128
             slice_duration : 0.0,           $ ; time for 1 slice                  offset: 132
             toffset     : 0.0,              $ ; time axis shift                   offset: 136
             glmax       : 0L,               $ ; unused                            offset: 140
             glmin       : 0L,               $ ; unused                            offset: 144

             ;   struct data_history 
             descrip     : bytarr(80),       $ ; any text you like                 offset: 148  
             aux_file    : bytarr(24),       $ ; auxiliary filename                offset: 228
             qform_code  : 0s,               $ ; Nifti Xform code                  offset: 252
             sform_code  : 0s,               $ ; Nifti Xform code                  offset: 254
             quatern_b   : 0.0,              $ ; Quaternion b parameter            offset: 256
             quatern_c   : 0.0,              $ ; Quaternion c parameter            offset: 260
             quatern_d   : 0.0,              $ ; Quaternion d parameter            offset: 264
             quatern_x   : 0.0,              $ ; Quaternion x parameter            offset: 268
             quatern_y   : 0.0,              $ ; Quaternion y parameter            offset: 272
             quatern_z   : 0.0,              $ ; Quaternion z parameter            offset: 276
             srow_x      : fltarr(4),        $ ; 1st row affine transformation     offset: 280
             srow_y      : fltarr(4),        $ ; 2nd row affine transformation     offset: 296
             srow_z      : fltarr(4),        $ ; 3rd row affine transformation     offset: 312

             ;; 
             intent_name : bytarr(16),       $ ; name or meaning of data           offset: 328
             magic       : bytarr(4)         $ ;  !< MUST be "ni1\0" or "n+1\0".   offset: 344
                                               ; "ni1\0" = [110,105,49,0]
                                               ; "n+1\0" = [110, 43,49,0]

          }                                                                     ;; bytes total: 348
  

  ;
  ;  Class instance vars:
  ;
  class = { IDLffNifti,              $ 
            endian  : 0,              $ ;  1=big-endian, 0=little-endian (read/write)
            hdr     : {header_nifti}, $ ;  header data
            default : {header_nifti}, $ ;  default header values
            data    : ptr_new()       $ ;  actual image data
          }
end


;
;  End idlffanalyze__define.pro
;
