;+
; NAME: struc_modifyobject__define
;
;
;
; PURPOSE: define the "struc_modifyobject_info" structure. Define it as an
; object, because otherwise it has to be defined several times, and it
; is hard to maintain same tages everywhere. This keeps track of the
; several things needed
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
;; GJ, 2019/2/19
;
PRO  struc_modifyobject__define
  
  s ={STRUC_MODIFYOBJECT,      $
    modify_step:0L, $;record 100 steps for modifying 30 3D objects
    current_step_index:0L, $;record current step index
    exist_num: INTARR(100)*0, $;100 steps, 1 drawn; 0 not exist ; this is identical to N_ELEMENTS(xd_sState.oPoly)
    hide_status: INTARR(100, 30)*0, $;100 steps, hide=0 show; hide=1 hide
;    color_status: BYTARR(100, 30, 3), $; the color of each 3D object
    polygon_list: PTRARR(30), $; the polygon list
    p_ori_verts_list: PTRARR(30), $; the original verts of each 3D object
;    p_ori_conn_list: PTRARR(30), $; the original conn of each 3D object
;    verts_min_list: DBLARR(30, 3), $; min verts at x,y,z for whole centering and scaling
;    verts_max_list: DBLARR(30, 3), $; max verts at x,y,z for whole centering and scaling
    p_maskindex_list: PTRARR(30), $; record all non-zero voxels as a mask
    countmask_list: LONARR(30)*0 $; the number of voxels in the mask
  }
END

;GJ 2019/219
;s ={STRUC_MODIFYOBJECT,      $
;  modify_step:0L, $;record 100 steps for modifying 30 3D objects
;  exist_status: INTARR(100, 30)*0, $;100 steps, 1 drawn; 0 not exist
;  showhide_status: INTARR(100, 30)*0, $;100 steps, 1 show; 0 hide
;  smooth_status: FLTARR(100, 30), $; the smooth level of each 3D object
;  ;    smooth_status: FLTARR(100, 30), $; the smooth level history of each 3D object
;    transparent_status: FLTARR(100, 30), $; the smooth level history of each 3D object
;  color_status: BYTARR(100, 30, 3), $; the color of each 3D object
;  p_verts_list: PTRARR(30), $; the verts of each 3D object
;  p_conn_list: PTRARR(30), $; the conn of each 3D object
;  verts_min_list: DBLARR(30, 3), $; min verts at x,y,z
;  verts_max_list: DBLARR(30, 3), $; max verts at x,y,z
;  p_maskindex_list: PTRARR(30), $; record all non-zero voxels as a mask
;  countmask_list: LONARR(30)*0 $; the number of voxels in the mask
;}

pro struc_modifyobject::info
  self_tag = tag_names({struc_modifyobject})
  for i= 0, n_elements(self_tag) - 1 do print,self_tag(i),' = ',self.(i)
end ; of method struc_modifyobject.info

;**********************************************************************************

pro struc_modifyobject::set,_extra = extra

  self_tag = tag_names({struc_modifyobject})
  extra_tag = tag_names(extra)

  for i= 0, n_elements(self_tag) - 1 do begin
    j= 0
    repeat begin
      treffer = self_tag(i) eq extra_tag(j)
      if treffer then self.(i) = extra.(j) else j = j+1
    endrep until treffer or (j ge n_elements(extra_tag))
  endfor
end ; of method struc_modifyobject::set

;**********************************************************************************

pro test_struc_modifyobject

  IF N_ELEMENTS(struc_modifyobject) EQ 0 THEN BEGIN
    ;; create the struc_modifyobject structure
    cur_modify=0
    ;; create a struc_modifyobject structure
    struc_modifyobject = {struc_modifyobject}
    struc_modifyobject.hide_Status =+1
  ENDIF ELSE BEGIN
    n=N_ELEMENTS(struc_modifyobject)
    p=replicate({struc_modifyobject},n+1)
    p[0:n-1]=struc_modifyobject[0:n-1]
    cur_modify=n
    struc_modifyobject=p
  ENDELSE


end