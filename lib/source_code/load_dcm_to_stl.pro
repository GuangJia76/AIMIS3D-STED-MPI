
;;GJ, 2018/10/08
;moved from d_objworld21.pro
;;load images and output volume cube
PRO image_read, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal = p_modal

  ;  afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  ;  afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  IF STRLEN(afn) NE 0 THEN BEGIN
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    a_nii = FILE_SEARCH(afn, '*.nii', count=nrfile_nii)

    ;    b = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    ;    print, 'before nrfiles = ', nrfile
    sample_ind = 1.
    IF nrfile GT 5 THEN BEGIN
      ;GJ 2019/2/13 if there are two many files, only pick 1/3 of the files
      IF nrfile GT 400*3 THEN BEGIN
        sample_ind = 3
        nrfile = nrfile/sample_ind
        temp_a = STRARR(nrfile)
        FOR ai = 0, N_ELEMENTS(temp_a)-1 DO temp_a[ai] = a[ai*sample_ind] 
        a = temp_a
      ENDIF
      ;      IF nrfile LT 350 THEN BEGIN
      ;        a = b
      ;      ENDIF
      ;      IF nrfile GT 350 AND nrfile LE 600 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/2.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[2*k]
      ;      ENDIF
      ;      IF nrfile GT 600 AND nrfile LE 900 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/3.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[3*k]
      ;      ENDIF
      ;      IF nrfile GT 900 AND nrfile LE 1200 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/4.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[4*k]
      ;      ENDIF
      ;      IF nrfile GT 1200 AND nrfile LE 1500 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/5.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[5*k]
      ;      ENDIF
      ;      IF nrfile GT 1500 THEN BEGIN
      ;        infowarning = DIALOG_MESSAGE('Too many images, please double-check your images!', /ERROR)
      ;        RETURN
      ;      ENDIF
      ;      print, 'after nrfiles = ', nrfile
      
      IF QUERY_DICOM(a[0]) EQ 0 THEN BEGIN
        MPI_dicom_modify, a[0], fn_MPI_new
        obj = OBJ_NEW('IDLffDICOM', fn_MPI_new)
      ENDIF ELSE BEGIN
        obj = OBJ_NEW('IDLffDICOM', a[0])
      ENDELSE
     
      ; Get the row & column size of the image(s):
      temp = obj->GetValue('0028'x,'0010'x,/NO_COPY)
      IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0] ELSE rows = 256
      temp = obj->GetValue('0028'x,'0011'x,/NO_COPY)
      IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0] ELSE cols = 256

      OBJ_DESTROY, obj
      
      vol_HU = DBLARR(cols, rows, nrfile);GJ 2020/4/9, a big bug
      voxelsize = DBLARR(nrfile, 7) * 0.
      imagePos = DBLARR(nrfile, 3) * 0.
      imageOri = DBLARR(nrfile, 6) * 0.
      centerPos = DBLARR(nrfile, 3) * 0.
      origin_loc = DBLARR(1, 3) * 0.

      FOR i=0L, nrfile-1L DO BEGIN
        
        IF QUERY_DICOM(a[i]) EQ 0 THEN BEGIN
          MPI_dicom_modify, a[i], fn_MPI_new
          obj = OBJ_NEW('IDLffDICOM', fn_MPI_new)
          result=QUERY_DICOM(fn_MPI_new, info)
        ENDIF ELSE BEGIN
          obj = OBJ_NEW('IDLffDICOM', a[i])
          result=QUERY_DICOM(a[i], info)
        ENDELSE

        IF result NE 0 THEN BEGIN
          ; Get the row & column size of the image(s):
          rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
          cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
          order=1
          
          sliceTh = FLOAT(*(obj->GetValue('0018'x,'0050'x,/NO_COPY))[0])*sample_ind
          ;@GJ, 2023/5/15, calculating spacing
          temp = obj->GetValue('0018'x,'0088'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN spacing = FLOAT(*(obj->GetValue('0018'x,'0088'x,/NO_COPY))[0])*sample_ind ELSE spacing = sliceTh*sample_ind
          ;          spacing = FLOAT(*(obj->GetValue('0018'x,'0088'x,/NO_COPY))[0])
          pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
          pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
          temp = obj->GetValue('0018'x,'5100'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN ptPosition = STRING(*(obj->GetValue('0018'x,'5100'x,/NO_COPY))[0]) ELSE ptPosition = i
          temp = obj->GetValue('0018'x,'1041'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN sliceLo = STRING(*(obj->GetValue('0020'x,'1041'x,/NO_COPY))[0]) ELSE sliceLo = i*sliceTh
          temp = obj->GetValue('0010'x,'1010'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN p_age = STRING(*(obj->GetValue('0010'x,'1010'x,/NO_COPY))[0]) ELSE p_age='0'
          p_modal_temp = STRING(*(obj->GetValue('0008'x,'0060'x,/NO_COPY))[0])
          p_modal = STRCOMPRESS(p_modal_temp, /REMOVE_ALL)
          imagePos[i,*] = DOUBLE(STRSPLIT(*(obj->GetValue('0020'x,'0032'x,/NO_COPY))[0], '\',/extract));GJ, 2020/4/11, bug fixed
          imageOri[i,*] = DOUBLE(STRSPLIT(*(obj->GetValue('0020'x,'0037'x,/NO_COPY))[0], '\',/extract));GJ, 2020/4/11, bug fixed
          
          ;GJ, 2021/2/21, calculating center pos of each slice
          centerPos[i,*] = imagePos[i,*] + cols/2.*imageOri[i,0:2]*pixelSp[0] + rows/2.*imageOri[i,3:5]*pixelSp[1]
          origin_loc = origin_loc + centerPos[i,*]
          
          ;GJ, 2021/2/20
          ;image number
          value = 0
          value = DOUBLE(*(obj->GetValue('0020'x,'0013'x,/NO_COPY))[0])
          image_number= MAX(value)
          voxelsize[i, 6] = image_number
;          IF (obj->QueryValue('0018,0088')) EQ 2 THEN spacing = FLOAT(obj->GetValue('0018,0088')) ELSE spacing = 0.
;          IF (obj->QueryValue('0018,0050')) EQ 2 THEN sliceTh = FLOAT(obj->GetValue('0018,0050')) ELSE sliceTh = 0.
;          IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
;          IF (obj->QueryValue('0018,5100')) EQ 2 THEN ptPosition = STRING(obj->GetValue('0018,5100')) ELSE ptPosition = 'FFS'
;          IF (obj->QueryValue('0020,1041')) EQ 2 THEN sliceLo = STRING(obj->GetValue('0020,1041')) ELSE sliceLo = 0.
;          IF (obj->QueryValue('0010,1010')) EQ 2 THEN p_age = STRING(obj->GetValue('0010,1010')) ELSE p_age = '040Y'
;          IF (obj->QueryValue('0008,0060')) EQ 2 THEN p_modal = STRING(obj->GetValue('0008,0060')) ELSE p_modal = ''
;          IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
;          IF (obj->QueryValue('0020,0037')) EQ 2 THEN imageOri[i,*] = DOUBLE(obj->GetValue('0020,0037'))
          
;          print, 'Ori', TRANSPOSE(imageOri[i,*])
;          print, 'Pos', TRANSPOSE(imagePos[i,*])
          patient_age=FLOAT(STRMID(p_age, 0,3))
          OBJ_DESTROY, obj
          ;缁绢収鍠涢濠氬炊閹冨壖闁告帒妫滄ご鎼佹偝閸ヮ亶鍤㈤柛娆愮墬椤掓粎娑甸敓锟�          print, "voxel size: ", i, ':', voxelsize[i, *]
          IF sliceTh GT 0. THEN BEGIN
            voxelsize[i, 0] = pixelSp[0]
            voxelsize[i, 1] = pixelSp[1]
            IF spacing LE 0. THEN BEGIN
              voxelsize[i, 2] = sliceTh + spacing
            ENDIF ELSE BEGIN
              IF sliceTh GT 0. THEN voxelsize[i, 2] = sliceTh ELSE voxelsize[i, 2] = spacing ;sliceTh +
            ENDELSE
            voxelsize[i, 3] = rows
            voxelsize[i, 4] = cols
          ENDIF
          voxelsize[i, 5] = sliceLo
        ENDIF
      ENDFOR

      ;GJ 2021/2/21, get original loc by averaging the number of images
      origin_loc = origin_loc/FLOAT(nrfile)

      ;GJ, 2021/02/20, sort the array based on image_number
      ;GJ, 2021/03/05, fix the bug of DICOM images with same image number 0
      IF ABS(MAX(voxelsize[*, 6]) - MIN(voxelsize[*, 6])) GT 1 THEN BEGIN
        sort_imageNum_ind = SORT(voxelsize[*, 6])
      ENDIF ELSE BEGIN
        sort_imageNum_ind = SORT(voxelsize[*, 5])
      ENDELSE
      a=a[sort_imageNum_ind]
      voxelsize=voxelsize[sort_imageNum_ind, *]
      imagePos = imagePos[sort_imageNum_ind, *]
      imageOri = imageOri[sort_imageNum_ind, *]
      
      ;
      ;double-check the images
      sort_lo = SORT(voxelsize[*, 5])
      diff_lo = SHIFT(voxelsize[*, 5], -1) - voxelsize[*, 5]
      
      ;proceed however the error may happen
      IF ABS(MAX(diff_lo[0:nrfile-2L]) - MIN(diff_lo[0:nrfile-2L])) GT 0.001 THEN BEGIN
        infowarning = DIALOG_MESSAGE('Images do not have the same slice distance!', /ERROR)
;        RETURN
      ENDIF

;      ; Create the progress bar.
       progressbar = Obj_New('progressbar', Color='sky blue', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
;      ; Place the progress bar on the display.
       progressbar -> Start

      FOR i=0L, nrfile-1L DO BEGIN
        
        IF QUERY_DICOM(a[i]) EQ 0 THEN BEGIN
          MPI_dicom_modify, a[i], fn_MPI_new
          obj = OBJ_NEW('IDLffDICOM', fn_MPI_new)
          result=QUERY_DICOM(fn_MPI_new, info)
        ENDIF ELSE BEGIN
          obj = OBJ_NEW('IDLffDICOM', a[i])
          result=QUERY_DICOM(a[i], info)
        ENDELSE

        count=(i+1.)/nrfile*100.0
        progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

        IF result NE 0 THEN BEGIN
          ; Get the image data
          array=obj->GetValue('7fe0'x, '0010'x)
          ;GJ 2020/10/17, to add this to pick the real image, instead of small icon picture
          ;GJ 2021/2/13, pick data from the larger one
          vPixels=*array[N_ELEMENTS(array)-1]
          new_cols=SIZE(vPixels[*,0], /dim)
          vPixels=CONGRID(vPixels, cols, rows)
;          FOR ijk=N_ELEMENTS(array)-1,0,-1 DO BEGIN
;            vPixels=*array[ijk]
;            IF (SIZE(vPixels, /dim))[0] EQ cols THEN BREAK ;
;          ENDFOR
;          vPixels=*array[0]
          PTR_FREE, array

          ;convert grey pixel to HU
          ;rescale intercept
          temp = obj->GetValue('0028'x,'1052'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
            rescale_intercept = FLOAT(*(obj->GetValue('0028'x,'1052'x,/NO_COPY))[0])
          ENDIF ELSE BEGIN
            rescale_intercept = 0.
          ENDELSE
          
          temp = obj->GetValue('0028'x,'1053'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
            rescale_slope = FLOAT(*(obj->GetValue('0028'x,'1053'x,/NO_COPY))[0])
          ENDIF ELSE BEGIN
            rescale_slope = 1.
          ENDELSE
          
;          IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
;          ;rescale slope
;          IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 1.
          OBJ_DESTROY, obj
          vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)

          ;smooth the image and store to 3D
          vol_HU[*,*,i] = vPixels_HU
        ENDIF
      ENDFOR
      
      ;@GJ, 2023/3/28, for MPI change back to pixel intensity
      ;@GJ, 2023/5/16, for MPI mulitple 1000 times for normal display
      IF STRCMP(p_modal, 'MPI', 3, /FOLD_CASE) THEN vol_HU = (vol_HU-DOUBLE(rescale_intercept))/DOUBLE(rescale_slope)
      
      progressbar -> Destroy
      
      ;GJ 2021/2/13, fix the bug of calculating FOV
      voxelsize[*, 3] = voxelsize[*, 3]*new_cols/cols

      ;鐎电増顨滈崺灞筋灉閿濆洠锟介柣銊ュ閻即宕㈤敓锟�
      sliceDistance = ABS(voxelsize[nrfile/2-1, 5] - voxelsize[nrfile/2, 5])


      ImageOri_temp = imagePos[nrfile/2,*] - imagePos[nrfile/2-1,*]
      ImageOri_z = ImageOri_temp/SQRT(ImageOri_temp[0]^2 + ImageOri_temp[1]^2 + ImageOri_temp[2]^2)

      x_d = [TRANSPOSE(imageOri[nrfile/2,0:2])]
      y_d = [TRANSPOSE(imageOri[nrfile/2,3:5])]
      z_d_cal = CROSSP(x_d, y_d)
      z_d = [TRANSPOSE(ImageOri_z)]
      IF TOTAL(z_d*z_d_cal) LT 0. THEN BEGIN
        vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)
        direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(-ImageOri_z)]]
      ENDIF ELSE BEGIN
        direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(ImageOri_z)]]
      ENDELSE
;      origin_loc = imagePos[0,*]
      print, 'direction = ', direction
      

      ;濠碘�鍊归悘澶愬及椤栨繂澹栭柛蹇撶墣缁绘﹢鏁嶅畝瀣诞闁挎冻璁ｅЧ鏌ユ晬鐏炴儳惟闁搞儲鍎抽崕姘跺磹閹烘洘宕抽柡澶涙嫹
;      IF STRCMP(ptPosition, 'FFS') THEN vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)

;      IF (voxelsize[nrfile/2, 5] - voxelsize[nrfile/2-1, 5]) LT 0. THEN vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)
      
      FOV = voxelsize[nrfile/2-1, 3]*voxelsize[nrfile/2-1, 0]
      length = nrfile * sliceDistance
      IF FLOOR(length*256/FOV) LE 256 THEN BEGIN
        vol_HU_iso_resolution = FOV/256.
        vol_HU_iso = CONGRID(vol_HU,  256, 256, FLOOR(length*256/FOV))
      ENDIF ELSE BEGIN
        vol_HU_iso_resolution = length/256.
        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*256/length), FLOOR(FOV*256/length), 256)
      ENDELSE
      print, 'FOV: ', vol_HU_iso_resolution*(SIZE(vol_HU_iso,/DIMENSIONS))[0], 'mm'
      print, 'length: ', vol_HU_iso_resolution*(SIZE(vol_HU_iso,/DIMENSIONS))[2], 'mm'
;      IF FLOOR(length*512/FOV) LE 512 THEN BEGIN
;        vol_HU_iso_resolution = FOV/512.
;        vol_HU_iso = CONGRID(vol_HU,  512, 512, FLOOR(length*512/FOV))
;      ENDIF ELSE BEGIN
;        vol_HU_iso_resolution = length/512.
;        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*512/length), FLOOR(FOV*512/length), 512)
;      ENDELSE
;      
      ;delete vol_HU
      vol_HU = !NULL

      WHILE vol_HU_iso_resolution LT 0.2 DO BEGIN
        dims = size(vol_HU_iso)
        vol_HU_iso_resolution = vol_HU_iso_resolution*2.
        vol_HU_iso = CONGRID(vol_HU_iso, dims[1]/2, dims[2]/2, dims[3]/2)
      ENDWHILE
     
      ;閻忓繐鏁妎l_HU_iso濠靛鍋勯幊鍡涘炊鐎电硶鏁勯柣褝鎷烽柣鈺冾焾缂嶅绂嶆惔锝傛晞婵﹫鎷稨U=-1000, 闁活枌鍔嶅鐢告偨缁暜aw[0:2]
      max_dim = MAX((SIZE(vol_HU_iso))[1:3], max_ind)
      btm_value = MIN(vol_HU_iso)-5;, -1000])
      vol_HU_cube = DBLARR(max_dim, max_dim, max_dim)*0. + btm_value; -1000.
      IF max_ind EQ 2 THEN BEGIN
        vol_HU_cube[(max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, (max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, *] = vol_HU_iso
;        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[1])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDIF ELSE BEGIN
        vol_HU_cube[*, *, (max_dim-(SIZE(vol_HU_iso))[3])/2:(max_dim-(SIZE(vol_HU_iso))[3])/2+(SIZE(vol_HU_iso))[3]-1] = vol_HU_iso
 ;       origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[3])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDELSE

;     vol_HU_cube = vol_HU_cube - MIN(vol_HU_cube) + 1
;      print, 'min_value = ', MIN(vol_HU_cube)
;      print, 'max_value = ', MAX(vol_HU_cube)
      ;print, 'size(vol_HU_cube) = ', size(vol_HU_cube)
      size_vol = SIZE(vol_HU_cube, /dimensions)
      vol_HU_cube[*, *, 0] = btm_value
      vol_HU_cube[*, *, size_vol[2]-1] =  btm_value
      vol_HU_cube[*, 0, *] = btm_value
      vol_HU_cube[*, size_vol[1]-1, *] =  btm_value
      vol_HU_cube[0, *, *] = btm_value
      vol_HU_cube[size_vol[0]-1, *, *] =  btm_value
      ;GJ 2020/7/7 test the image cube
;      iimage, REFORM(vol_HU_cube[*,127,*]) 
;      direction_1=[[direction[0:2]], [direction[3:5]], [direction[6:8]]]
;      direction_0=[[1.,0,0], [0,1.,0], [0,0,1.]]
;      kernal = TRANSPOSE(direction_1)
;      R_inv = TRANSPOSE(direction_0) # INVERT(kernal)
;      angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;      angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;      angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;
;      IF ABS(cos(angle_beta*!PI/180.) - R_inv[2,2]) GT 0.0001 THEN angle_beta = 180. - angle_beta
;      IF ABS(sin(angle_beta*!PI/180.)*cos(angle_alpha*!PI/180.) - R_inv[1,2]) GT 0.0001 THEN angle_alpha = 180. - angle_alpha
;      IF ABS(-sin(angle_beta*!PI/180.)*cos(angle_gamma*!PI/180.) - R_inv[2,1]) GT 0.0001 THEN angle_gamma = 180. - angle_gamma
;
;      TRANSFORM_VOLUME, vol_HU_cube, rot_z, Rotation=[0, 0, angle_gamma]
;      TRANSFORM_VOLUME, rot_z, rot_x, Rotation=[angle_beta, 0, 0]
;      TRANSFORM_VOLUME, rot_x, vol_HU_cube, Rotation=[0, 0, angle_alpha]

;      R_inv = INVERT(curr_ori)
;      angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;      angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;      angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 0, angle_gamma])
;      rot_x = TRANSFORM_VOLUME(rot_z, MISSING = btm_value, Rotation=[angle_beta, 0, 0])
;      vol_HU_cube = TRANSFORM_VOLUME(rot_x, MISSING = btm_value, Rotation=[0, 0, angle_alpha])
      
;      print, 'Ori z', TRANSPOSE(ImageOri_z)
;      IF imageOri[nrfile/2,4] EQ 0 THEN imageOri[nrfile/2,4] = 0.00001
;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, Rotation=[0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
;      print, 'rot_z', [0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])]
;      IF imageOri[nrfile/2,1] EQ 0 THEN imageOri[nrfile/2,1] = 0.00001
;      rot_x = TRANSFORM_VOLUME(rot_z, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0])
;      print, 'rot_x', [180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0]
;      IF imageOri[nrfile/2,0] EQ 0 THEN imageOri[nrfile/2,0] = 0.00001
;      rot_y = TRANSFORM_VOLUME(rot_x, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0])
;      print, 'rot_y', [0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0]
;
;      vol_HU_cube = rot_y
;      TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
;      vol_HU_cube = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])


      
      ;delete vol_HU
      vol_HU_iso = !NULL
      ;for black blood
;     vol_HU_cube = MAX(vol_HU_cube) - vol_HU_cube


      vol_HU_cube_resolution = vol_HU_iso_resolution
      print, 'resolution', vol_HU_cube_resolution
    ENDIF ELSE BEGIN
      ;@GJ, 2023/4/4, loading nifti format data
      ;a_nii = FILE_SEARCH(afn, '*.nii', count=nrfile_nii)
      IF nrfile_nii EQ 0 THEN BEGIN
        print, 'no nii file'
        return
      ENDIF
      print, 'nii read'
      vol_HU = read_nifti(a_nii, header=hdr)
      FOV = hdr.dim[1]*hdr.pixdim[1]
      length = hdr.dim[3]*hdr.pixdim[3]
      IF FLOOR(length*256/FOV) LE 256 THEN BEGIN
        vol_HU_iso_resolution = FOV/256.
        vol_HU_iso = CONGRID(vol_HU,  256, 256, FLOOR(length*256/FOV))
      ENDIF ELSE BEGIN
        vol_HU_iso_resolution = length/256.
        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*256/length), FLOOR(FOV*256/length), 256)
      ENDELSE
      
      ;delete vol_HU
      vol_HU = !NULL

      WHILE vol_HU_iso_resolution LT 0.2 DO BEGIN
        dims = size(vol_HU_iso)
        vol_HU_iso_resolution = vol_HU_iso_resolution*2.
        vol_HU_iso = CONGRID(vol_HU_iso, dims[1]/2, dims[2]/2, dims[3]/2)
      ENDWHILE

      max_dim = MAX((SIZE(vol_HU_iso))[1:3], max_ind)
      btm_value = MIN(vol_HU_iso)-5;, -1000])
      vol_HU_cube = DBLARR(max_dim, max_dim, max_dim)*0. + btm_value; -1000.
      IF max_ind EQ 2 THEN BEGIN
        vol_HU_cube[(max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, (max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, *] = vol_HU_iso
      ENDIF ELSE BEGIN
        vol_HU_cube[*, *, (max_dim-(SIZE(vol_HU_iso))[3])/2:(max_dim-(SIZE(vol_HU_iso))[3])/2+(SIZE(vol_HU_iso))[3]-1] = vol_HU_iso
      ENDELSE

      size_vol = SIZE(vol_HU_cube, /dimensions)
      vol_HU_cube[*, *, 0] = btm_value
      vol_HU_cube[*, *, size_vol[2]-1] =  btm_value
      vol_HU_cube[*, 0, *] = btm_value
      vol_HU_cube[*, size_vol[1]-1, *] =  btm_value
      vol_HU_cube[0, *, *] = btm_value
      vol_HU_cube[size_vol[0]-1, *, *] =  btm_value

      ;delete vol_HU
      vol_HU_iso = !NULL
      vol_HU_cube_resolution = vol_HU_iso_resolution
      print, 'resolution', vol_HU_cube_resolution
      ;randomly set patient age
      patient_age = '28'
    ENDELSE
  ENDIF ELSE BEGIN
;    print, 'test'

    ;    return,[]
  ENDELSE
END

;@GJ, 2023/05/24, modify MPI non standard dicom file to a standard file
PRO MPI_dicom_modify, MPI_DCM_fn, fn_MPI_new

  IF QUERY_DICOM(MPI_DCM_fn) EQ 0 THEN BEGIN
    data_1 = bytarr(158)
    data_1[128:*]=[68B, 73B, 67B, 77B, 2B, 0B, 0B, 0B, 85B, 76B, 4B, 0B, 198B, 0B, 0B, 0B, 2B, 0B, 1B, 0B, 79B, 66B, 0B, 0B, 2B, 0B, 0B, 0B, 0B, 1B]

    fn_MPI = MPI_DCM_fn
    data_2 = read_binary(fn_MPI, DATA_START=132)

    ; Use Count to get the number of nonzero elements:
    index = WHERE(data_2 EQ 8, count)
    ; Only subscript the array if it is safe:
    ;define the File Meta Group Element Length 0002,0000
    IF count NE 0 THEN BEGIN
      FOR j=0, count-1 DO BEGIN
        IF data_2[index[j]+1] EQ 0 THEN BEGIN
          data_1[140]=BYTE(index[j])+14B
          break
        ENDIF
      ENDFOR
    ENDIF

    ;@GJ, save the modified dicom file
    afn_des = FILE_DIRNAME(fn_MPI)
    pos_dir = STRPOS(STRMID(afn_des, 0, STRLEN(afn_des)-1), '\', /REVERSE_SEARCH)
    cur_dir = STRMID(afn_des, 0, pos_dir+1)
    fn_MPI_new = cur_dir+FILE_BASENAME(MPI_DCM_fn)
    OPENW, U, fn_MPI_new,/GET_LUN
    WRITEU, U, data_1
    WRITEU, U, data_2
    FREE_LUN, U
    ;          IF QUERY_DICOM(fn_MPI_new) THEN FILE_COPY, fn_MPI_new, a[i],/ALLOW_SAME, /OVERWRITE
    ;@GJ, 2023/5/18, replace the sorted file as the new dicom file
;    a[i]=fn_MPI_new
  ENDIF
END


;;GJ, 2018/10/08
;moved from d_objworld21.pro
;;load images and output volume cube
PRO image_read_512, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal = p_modal

  ;  afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  ;  afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  IF STRLEN(afn) NE 0 THEN BEGIN
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)

    ;    b = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    ;    print, 'before nrfiles = ', nrfile
    IF nrfile NE 0 THEN BEGIN
      ;GJ 2019/2/13 if there are two many files, only pick 1/3 of the files
      IF nrfile GT 400 THEN BEGIN
        sample_ind = 3
        nrfile = nrfile/sample_ind
        temp_a = STRARR(nrfile)
        FOR ai = 0, N_ELEMENTS(temp_a)-1 DO temp_a[ai] = a[ai*sample_ind] 
        a = temp_a
      ENDIF
      ;      IF nrfile LT 350 THEN BEGIN
      ;        a = b
      ;      ENDIF
      ;      IF nrfile GT 350 AND nrfile LE 600 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/2.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[2*k]
      ;      ENDIF
      ;      IF nrfile GT 600 AND nrfile LE 900 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/3.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[3*k]
      ;      ENDIF
      ;      IF nrfile GT 900 AND nrfile LE 1200 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/4.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[4*k]
      ;      ENDIF
      ;      IF nrfile GT 1200 AND nrfile LE 1500 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/5.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[5*k]
      ;      ENDIF
      ;      IF nrfile GT 1500 THEN BEGIN
      ;        infowarning = DIALOG_MESSAGE('Too many images, please double-check your images!', /ERROR)
      ;        RETURN
      ;      ENDIF
      ;      print, 'after nrfiles = ', nrfile
      obj = OBJ_NEW('IDLffDICOM', a[0])
      
      ; Get the row & column size of the image(s):
      rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
      cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
      OBJ_DESTROY, obj
      
      vol_HU = DBLARR(rows, cols, nrfile)
      voxelsize = DBLARR(nrfile, 6) * 0.
      imagePos = DBLARR(nrfile, 3) * 0.
      imageOri = DBLARR(nrfile, 6) * 0.

      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOM', a[i])
        result=QUERY_DICOM(a[i], info)

        IF result NE 0 THEN BEGIN
          ; Get the row & column size of the image(s):
          rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
          cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
          order=1
          
          spacing = FLOAT(*(obj->GetValue('0018'x,'0088'x,/NO_COPY))[0])
          sliceTh = FLOAT(*(obj->GetValue('0018'x,'0050'x,/NO_COPY))[0])
          pixelSpa = STRING(*(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0])
          pixelSp=FLOAT(STRSPLIT(pixelSpa, '\', /extract))
          ptPosition = STRING(*(obj->GetValue('0018'x,'5100'x,/NO_COPY))[0])
          sliceLo = STRING(*(obj->GetValue('0020'x,'1041'x,/NO_COPY))[0])
          p_age = STRING(*(obj->GetValue('0010'x,'1010'x,/NO_COPY))[0])
          p_modal = STRING(*(obj->GetValue('0008'x,'0060'x,/NO_COPY))[0])
          imagePos[i,*] = DOUBLE(*(obj->GetValue('0020'x,'0032'x,/NO_COPY))[0])
          imageOri[i,*] = DOUBLE(*(obj->GetValue('0020'x,'0037'x,/NO_COPY))[0])
          
;          IF (obj->QueryValue('0018,0088')) EQ 2 THEN spacing = FLOAT(obj->GetValue('0018,0088')) ELSE spacing = 0.
;          IF (obj->QueryValue('0018,0050')) EQ 2 THEN sliceTh = FLOAT(obj->GetValue('0018,0050')) ELSE sliceTh = 0.
;          IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
;          IF (obj->QueryValue('0018,5100')) EQ 2 THEN ptPosition = STRING(obj->GetValue('0018,5100')) ELSE ptPosition = 'FFS'
;          IF (obj->QueryValue('0020,1041')) EQ 2 THEN sliceLo = STRING(obj->GetValue('0020,1041')) ELSE sliceLo = 0.
;          IF (obj->QueryValue('0010,1010')) EQ 2 THEN p_age = STRING(obj->GetValue('0010,1010')) ELSE p_age = '040Y'
;          IF (obj->QueryValue('0008,0060')) EQ 2 THEN p_modal = STRING(obj->GetValue('0008,0060')) ELSE p_modal = ''
;          IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
;          IF (obj->QueryValue('0020,0037')) EQ 2 THEN imageOri[i,*] = DOUBLE(obj->GetValue('0020,0037'))
          
;          print, 'Ori', TRANSPOSE(imageOri[i,*])
;          print, 'Pos', TRANSPOSE(imagePos[i,*])
          patient_age=FLOAT(STRMID(p_age, 0,3))
          OBJ_DESTROY, obj
          ;缁绢収鍠涢濠氬炊閹冨壖闁告帒妫滄ご鎼佹偝閸ヮ亶鍤㈤柛娆愮墬椤掓粎娑甸敓锟�          print, "voxel size: ", i, ':', voxelsize[i, *]
          IF sliceTh GT 0. THEN BEGIN
            voxelsize[i, 0] = pixelSp[0]
            voxelsize[i, 1] = pixelSp[1]
            IF spacing LE 0. THEN voxelsize[i, 2] = sliceTh + spacing ELSE voxelsize[i, 2] = spacing ;sliceTh +
            voxelsize[i, 3] = rows
            voxelsize[i, 4] = cols
          ENDIF

          voxelsize[i, 5] = sliceLo
        ENDIF
      ENDFOR

      ;double-check the images
      sort_lo = SORT(voxelsize[*, 5])
      diff_lo = SHIFT(voxelsize[*, 5], -1) - voxelsize[*, 5]
      
      ;proceed however the error may happen
      IF ABS(MAX(diff_lo[0:nrfile-2L]) - MIN(diff_lo[0:nrfile-2L])) GT 0.001 THEN BEGIN
        infowarning = DIALOG_MESSAGE('Images do not have the same slice distance!', /ERROR)
;        RETURN
      ENDIF

;      ; Create the progress bar.
       progressbar = Obj_New('progressbar', Color='sky blue', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
;      ; Place the progress bar on the display.
       progressbar -> Start

      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOM', a[i])
        result=QUERY_DICOM(a[i], info)

        count=(i+1.)/nrfile*100.0
        progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

        IF result NE 0 THEN BEGIN
          ; Get the image data
          array=obj->GetValue('7fe0'x, '0010'x)
          ;GJ 2020/10/17, to add this to pick the real image, instead of small icon picture
          FOR i=0, N_ELEMENTS(array)-1 DO BEGIN
            vPixels=*array[i]
            IF (SIZE(vPixels, /dim))[0] EQ cols THEN BREAK ;
          ENDFOR
;          vPixels=*array[0]
          PTR_FREE, array
          
          temp = obj->GetValue('0028'x,'1052'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
            rescale_intercept = FLOAT(*(obj->GetValue('0028'x,'1052'x,/NO_COPY))[0])
          ENDIF ELSE BEGIN
            rescale_intercept = 0.
          ENDELSE

          temp = obj->GetValue('0028'x,'1053'x,/NO_COPY)
          IF (SIZE(temp,/N_DIMENSIONS) NE 0) THEN BEGIN
            rescale_slope = FLOAT(*(obj->GetValue('0028'x,'1053'x,/NO_COPY))[0])
          ENDIF ELSE BEGIN
            rescale_slope = 1.
          ENDELSE

;          ;convert grey pixel to HU
;          ;rescale intercept
;          rescale_intercept = FLOAT(*(obj->GetValue('0028'x,'1052'x,/NO_COPY))[0])
;          rescale_slope = FLOAT(*(obj->GetValue('0028'x,'1053'x,/NO_COPY))[0])
;          IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
;          ;rescale slope
;          IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 1.
          OBJ_DESTROY, obj
          vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)

          ;smooth the image and store to 3D
          vol_HU[*,*,i] = vPixels_HU
        ENDIF
      ENDFOR

      progressbar -> Destroy
      
      ;鐎电増顨滈崺灞筋灉閿濆洠锟介柣銊ュ閻即宕㈤敓锟�
      sliceDistance = ABS(voxelsize[nrfile/2-1, 5] - voxelsize[nrfile/2, 5])


      ImageOri_temp = imagePos[nrfile/2,*] - imagePos[nrfile/2-1,*]
      ImageOri_z = ImageOri_temp/SQRT(ImageOri_temp[0]^2 + ImageOri_temp[1]^2 + ImageOri_temp[2]^2)

      x_d = [TRANSPOSE(imageOri[nrfile/2,0:2])]
      y_d = [TRANSPOSE(imageOri[nrfile/2,3:5])]
      z_d_cal = CROSSP(x_d, y_d)
      z_d = [TRANSPOSE(ImageOri_z)]
      IF TOTAL(z_d*z_d_cal) LT 0. THEN BEGIN
        vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)
        direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(-ImageOri_z)]]
      ENDIF ELSE BEGIN
        direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(ImageOri_z)]]
      ENDELSE
      origin_loc = imagePos[0,*]
      print, 'direction = ', direction
      

      ;濠碘�鍊归悘澶愬及椤栨繂澹栭柛蹇撶墣缁绘﹢鏁嶅畝瀣诞闁挎冻璁ｅЧ鏌ユ晬鐏炴儳惟闁搞儲鍎抽崕姘跺磹閹烘洘宕抽柡澶涙嫹
;      IF STRCMP(ptPosition, 'FFS') THEN vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)

;      IF (voxelsize[nrfile/2, 5] - voxelsize[nrfile/2-1, 5]) LT 0. THEN vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)
      
      FOV = voxelsize[nrfile/2-1, 3]*voxelsize[nrfile/2-1, 0]
      length = nrfile * sliceDistance
      IF FLOOR(length*512/FOV) LE 512 THEN BEGIN
        vol_HU_iso_resolution = FOV/512.
        vol_HU_iso = CONGRID(vol_HU,  512, 512, FLOOR(length*512/FOV))
      ENDIF ELSE BEGIN
        vol_HU_iso_resolution = length/512.
        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*512/length), FLOOR(FOV*512/length), 512)
      ENDELSE

;      IF FLOOR(length*512/FOV) LE 512 THEN BEGIN
;        vol_HU_iso_resolution = FOV/512.
;        vol_HU_iso = CONGRID(vol_HU,  512, 512, FLOOR(length*512/FOV))
;      ENDIF ELSE BEGIN
;        vol_HU_iso_resolution = length/512.
;        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*512/length), FLOOR(FOV*512/length), 512)
;      ENDELSE
;      
      ;delete vol_HU
      vol_HU = !NULL

      WHILE vol_HU_iso_resolution LT 0.4 DO BEGIN
        dims = size(vol_HU_iso)
        vol_HU_iso_resolution = vol_HU_iso_resolution*2.
        vol_HU_iso = CONGRID(vol_HU_iso, dims[1]/2, dims[2]/2, dims[3]/2)
      ENDWHILE
     
      ;閻忓繐鏁妎l_HU_iso濠靛鍋勯幊鍡涘炊鐎电硶鏁勯柣褝鎷烽柣鈺冾焾缂嶅绂嶆惔锝傛晞婵﹫鎷稨U=-1000, 闁活枌鍔嶅鐢告偨缁暜aw[0:2]
      max_dim = MAX((SIZE(vol_HU_iso))[1:3], max_ind)
      btm_value = MIN(vol_HU_iso)-5;, -1000])
      vol_HU_cube = DBLARR(max_dim, max_dim, max_dim)*0. + btm_value; -1000.
      IF max_ind EQ 2 THEN BEGIN
        vol_HU_cube[(max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, (max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, *] = vol_HU_iso
        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[1])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDIF ELSE BEGIN
        vol_HU_cube[*, *, (max_dim-(SIZE(vol_HU_iso))[3])/2:(max_dim-(SIZE(vol_HU_iso))[3])/2+(SIZE(vol_HU_iso))[3]-1] = vol_HU_iso
        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[3])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDELSE

;     vol_HU_cube = vol_HU_cube - MIN(vol_HU_cube) + 1
;      print, 'min_value = ', MIN(vol_HU_cube)
;      print, 'max_value = ', MAX(vol_HU_cube)
      ;print, 'size(vol_HU_cube) = ', size(vol_HU_cube)
      size_vol = SIZE(vol_HU_cube, /dimensions)
      vol_HU_cube[*, *, 0] = btm_value
      vol_HU_cube[*, *, size_vol[2]-1] =  btm_value
      vol_HU_cube[*, 0, *] = btm_value
      vol_HU_cube[*, size_vol[1]-1, *] =  btm_value
      vol_HU_cube[0, *, *] = btm_value
      vol_HU_cube[size_vol[0]-1, *, *] =  btm_value
;      direction_1=[[direction[0:2]], [direction[3:5]], [direction[6:8]]]
;      direction_0=[[1.,0,0], [0,1.,0], [0,0,1.]]
;      kernal = TRANSPOSE(direction_1)
;      R_inv = TRANSPOSE(direction_0) # INVERT(kernal)
;      angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;      angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;      angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;
;      IF ABS(cos(angle_beta*!PI/180.) - R_inv[2,2]) GT 0.0001 THEN angle_beta = 180. - angle_beta
;      IF ABS(sin(angle_beta*!PI/180.)*cos(angle_alpha*!PI/180.) - R_inv[1,2]) GT 0.0001 THEN angle_alpha = 180. - angle_alpha
;      IF ABS(-sin(angle_beta*!PI/180.)*cos(angle_gamma*!PI/180.) - R_inv[2,1]) GT 0.0001 THEN angle_gamma = 180. - angle_gamma
;
;      TRANSFORM_VOLUME, vol_HU_cube, rot_z, Rotation=[0, 0, angle_gamma]
;      TRANSFORM_VOLUME, rot_z, rot_x, Rotation=[angle_beta, 0, 0]
;      TRANSFORM_VOLUME, rot_x, vol_HU_cube, Rotation=[0, 0, angle_alpha]

;      R_inv = INVERT(curr_ori)
;      angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;      angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;      angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 0, angle_gamma])
;      rot_x = TRANSFORM_VOLUME(rot_z, MISSING = btm_value, Rotation=[angle_beta, 0, 0])
;      vol_HU_cube = TRANSFORM_VOLUME(rot_x, MISSING = btm_value, Rotation=[0, 0, angle_alpha])
      
;      print, 'Ori z', TRANSPOSE(ImageOri_z)
;      IF imageOri[nrfile/2,4] EQ 0 THEN imageOri[nrfile/2,4] = 0.00001
;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, Rotation=[0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
;      print, 'rot_z', [0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])]
;      IF imageOri[nrfile/2,1] EQ 0 THEN imageOri[nrfile/2,1] = 0.00001
;      rot_x = TRANSFORM_VOLUME(rot_z, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0])
;      print, 'rot_x', [180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0]
;      IF imageOri[nrfile/2,0] EQ 0 THEN imageOri[nrfile/2,0] = 0.00001
;      rot_y = TRANSFORM_VOLUME(rot_x, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0])
;      print, 'rot_y', [0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0]
;
;      vol_HU_cube = rot_y
;      TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
;      vol_HU_cube = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])


      
      ;delete vol_HU
      vol_HU_iso = !NULL
      ;for black blood
;     vol_HU_cube = MAX(vol_HU_cube) - vol_HU_cube


      vol_HU_cube_resolution = vol_HU_iso_resolution
      print, 'resolution', vol_HU_cube_resolution
    ENDIF
  ENDIF ELSE BEGIN
    ;    return,[]
  ENDELSE
END

;;GJ, 2018/10/08
;moved from d_objworld21.pro
;;load images and output volume cube
PRO image_read_512_DCE, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal = p_modal

  ;  afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  ;  afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  IF STRLEN(afn) NE 0 THEN BEGIN
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)

    ;    b = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    ;    print, 'before nrfiles = ', nrfile
    IF nrfile NE 0 THEN BEGIN
      ;GJ 2019/2/13 if there are two many files, only pick 1/3 of the files
      IF nrfile GT 400 THEN BEGIN
        sample_ind = 3
        nrfile = nrfile/sample_ind
        temp_a = STRARR(nrfile)
        FOR ai = 0, N_ELEMENTS(temp_a)-1 DO temp_a[ai] = a[ai*sample_ind]
        a = temp_a
      ENDIF
      obj = OBJ_NEW('IDLffDICOMex', a[0], /NO_PIXEL_DATA)
      obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
      OBJ_DESTROY, obj
      vol_HU = DBLARR(rows, cols, nrfile)
      voxelsize = DBLARR(nrfile, 6) * 0.
      imagePos = DBLARR(nrfile, 3) * 0.
      imageOri = DBLARR(nrfile, 6) * 0.

      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOMex', a[i])
        result=QUERY_DICOM(a[i], info)

        IF result NE 0 THEN BEGIN
          obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
          order=1

          IF (obj->QueryValue('0018,0088')) EQ 2 THEN spacing = FLOAT(obj->GetValue('0018,0088')) ELSE spacing = 0.
          IF (obj->QueryValue('0018,0050')) EQ 2 THEN sliceTh = FLOAT(obj->GetValue('0018,0050')) ELSE sliceTh = 0.
          IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
          IF (obj->QueryValue('0018,5100')) EQ 2 THEN ptPosition = STRING(obj->GetValue('0018,5100')) ELSE ptPosition = 'FFS'
          IF (obj->QueryValue('0020,1041')) EQ 2 THEN sliceLo = STRING(obj->GetValue('0020,1041')) ELSE sliceLo = 0.
          IF (obj->QueryValue('0010,1010')) EQ 2 THEN p_age = STRING(obj->GetValue('0010,1010')) ELSE p_age = '040Y'
          IF (obj->QueryValue('0008,0060')) EQ 2 THEN p_modal = STRING(obj->GetValue('0008,0060')) ELSE p_modal = ''
          IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
          IF (obj->QueryValue('0020,0037')) EQ 2 THEN imageOri[i,*] = DOUBLE(obj->GetValue('0020,0037'))

          patient_age=FLOAT(STRMID(p_age, 0,3))
          OBJ_DESTROY, obj
          ;缁绢収鍠涢濠氬炊閹冨壖闁告帒妫滄ご鎼佹偝閸ヮ亶鍤㈤柛娆愮墬椤掓粎娑甸敓锟�          print, "voxel size: ", i, ':', voxelsize[i, *]
          IF sliceTh GT 0. THEN BEGIN
            voxelsize[i, 0] = pixelSp[0]
            voxelsize[i, 1] = pixelSp[1]
            IF spacing LE 0. THEN voxelsize[i, 2] = sliceTh + spacing ELSE voxelsize[i, 2] = spacing ;sliceTh +
            voxelsize[i, 3] = rows
            voxelsize[i, 4] = cols
          ENDIF

          voxelsize[i, 5] = sliceLo
        ENDIF
      ENDFOR

      ;double-check the images
      sort_lo = SORT(voxelsize[*, 5])
      diff_lo = SHIFT(voxelsize[*, 5], -1) - voxelsize[*, 5]

      ;proceed however the error may happen
      IF ABS(MAX(diff_lo[0:nrfile-2L]) - MIN(diff_lo[0:nrfile-2L])) GT 0.001 THEN BEGIN
        infowarning = DIALOG_MESSAGE('Images do not have the same slice distance!', /ERROR)
        ;        RETURN
      ENDIF

      ;      ; Create the progress bar.
      progressbar = Obj_New('progressbar', Color='sky blue', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
      ;      ; Place the progress bar on the display.
      progressbar -> Start

      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOMex', a[i])
        result=QUERY_DICOM(a[i], info)

        count=(i+1.)/nrfile*100.0
        progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

        IF result NE 0 THEN BEGIN
          vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
          ;convert grey pixel to HU
          ;rescale intercept
          IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
          ;rescale slope
          IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 1.
          OBJ_DESTROY, obj
          vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)

          ;smooth the image and store to 3D
          vol_HU[*,*,i] = vPixels_HU
        ENDIF
      ENDFOR

      ;replacing with Amp parameter files
;      bfn='C:\D_drive\Bladder_UnetDCE\Bladderpatient\P04\B01_20090817\Dyn01\Dyn01_output\AUC\'
      bfn='C:\D_drive\Bladder_UnetDCE\Bladderpatient\P01\B01_20090708\Dyn01\dyn_output\AUC\'
      b = FILE_SEARCH(bfn, '*.dcm', count=nrfile)
      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOMex', b[i])
        result=QUERY_DICOM(b[i], info)

        IF result NE 0 THEN BEGIN
          vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
          ;convert grey pixel to HU
          ;rescale intercept
          IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
          ;rescale slope
          IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 1.
          OBJ_DESTROY, obj
          vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)

          ;smooth the image and store to 3D
          vol_HU[*,*,i] = vPixels_HU
        ENDIF
      ENDFOR


      progressbar -> Destroy

      ;鐎电増顨滈崺灞筋灉閿濆洠锟介柣銊ュ閻即宕㈤敓锟�
      sliceDistance = ABS(voxelsize[nrfile/2-1, 5] - voxelsize[nrfile/2, 5])


      ImageOri_temp = imagePos[nrfile/2,*] - imagePos[nrfile/2-1,*]
      ImageOri_z = ImageOri_temp/SQRT(ImageOri_temp[0]^2 + ImageOri_temp[1]^2 + ImageOri_temp[2]^2)

      x_d = [TRANSPOSE(imageOri[nrfile/2,0:2])]
      y_d = [TRANSPOSE(imageOri[nrfile/2,3:5])]
      z_d_cal = CROSSP(x_d, y_d)
      z_d = [TRANSPOSE(ImageOri_z)]
      IF TOTAL(z_d*z_d_cal) LT 0. THEN BEGIN
        vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)
        direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(-ImageOri_z)]]
      ENDIF ELSE BEGIN
        direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(ImageOri_z)]]
      ENDELSE
      origin_loc = imagePos[0,*]
      print, 'direction = ', direction


      FOV = voxelsize[nrfile/2-1, 3]*voxelsize[nrfile/2-1, 0]
      length = nrfile * sliceDistance
      IF FLOOR(length*512/FOV) LE 512 THEN BEGIN
        vol_HU_iso_resolution = FOV/512.
        vol_HU_iso = CONGRID(vol_HU,  512, 512, FLOOR(length*512/FOV))
      ENDIF ELSE BEGIN
        vol_HU_iso_resolution = length/512.
        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*512/length), FLOOR(FOV*512/length), 512)
      ENDELSE

      ;delete vol_HU
      vol_HU = !NULL

      max_dim = MAX((SIZE(vol_HU_iso))[1:3], max_ind)
      btm_value = MIN(vol_HU_iso)-5;, -1000])
      vol_HU_cube = DBLARR(max_dim, max_dim, max_dim)*0. + btm_value; -1000.
      IF max_ind EQ 2 THEN BEGIN
        vol_HU_cube[(max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, (max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, *] = vol_HU_iso
        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[1])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDIF ELSE BEGIN
        vol_HU_cube[*, *, (max_dim-(SIZE(vol_HU_iso))[3])/2:(max_dim-(SIZE(vol_HU_iso))[3])/2+(SIZE(vol_HU_iso))[3]-1] = vol_HU_iso
        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[3])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDELSE

      size_vol = SIZE(vol_HU_cube, /dimensions)
      vol_HU_cube[*, *, 0] = btm_value
      vol_HU_cube[*, *, size_vol[2]-1] =  btm_value
      vol_HU_cube[*, 0, *] = btm_value
      vol_HU_cube[*, size_vol[1]-1, *] =  btm_value
      vol_HU_cube[0, *, *] = btm_value
      vol_HU_cube[size_vol[0]-1, *, *] =  btm_value

      ;delete vol_HU
      vol_HU_iso = !NULL

      vol_HU_cube_resolution = vol_HU_iso_resolution
      print, 'resolution', vol_HU_cube_resolution
    ENDIF
  ENDIF ELSE BEGIN
    ;    return,[]
  ENDELSE
END



PRO image_read_neg, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal

  ;  afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  ;  afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  IF STRLEN(afn) NE 0 THEN BEGIn
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)

    ;    b = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    ;    print, 'before nrfiles = ', nrfile
    IF nrfile NE 0 THEN BEGIN
      ;      IF nrfile LT 350 THEN BEGIN
      ;        a = b
      ;      ENDIF
      ;      IF nrfile GT 350 AND nrfile LE 600 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/2.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[2*k]
      ;      ENDIF
      ;      IF nrfile GT 600 AND nrfile LE 900 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/3.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[3*k]
      ;      ENDIF
      ;      IF nrfile GT 900 AND nrfile LE 1200 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/4.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[4*k]
      ;      ENDIF
      ;      IF nrfile GT 1200 AND nrfile LE 1500 THEN BEGIN
      ;        nrfile = FLOOR(nrfile/5.)
      ;        a = STRARR(nrfile)
      ;        FOR k=0, nrfile-1 DO a[k] = b[5*k]
      ;      ENDIF
      ;      IF nrfile GT 1500 THEN BEGIN
      ;        infowarning = DIALOG_MESSAGE('Too many images, please double-check your images!', /ERROR)
      ;        RETURN
      ;      ENDIF
      ;      print, 'after nrfiles = ', nrfile
      obj = OBJ_NEW('IDLffDICOMex', a[0], /NO_PIXEL_DATA)
      obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
      OBJ_DESTROY, obj
      vol_HU = DBLARR(rows, cols, nrfile)
      voxelsize = DBLARR(nrfile, 6) * 0.
      imagePos = DBLARR(nrfile, 3) * 0.
      imageOri = DBLARR(nrfile, 6) * 0.

      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOMex', a[i])
        result=QUERY_DICOM(a[i], info)

        IF result NE 0 THEN BEGIN
          obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
          order=1

          IF (obj->QueryValue('0018,0088')) EQ 2 THEN spacing = FLOAT(obj->GetValue('0018,0088')) ELSE spacing = 0.
          IF (obj->QueryValue('0018,0050')) EQ 2 THEN sliceTh = FLOAT(obj->GetValue('0018,0050')) ELSE sliceTh = 0.
          IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
          IF (obj->QueryValue('0018,5100')) EQ 2 THEN ptPosition = STRING(obj->GetValue('0018,5100')) ELSE ptPosition = 'FFS'
          IF (obj->QueryValue('0020,1041')) EQ 2 THEN sliceLo = STRING(obj->GetValue('0020,1041')) ELSE sliceLo = 0.
          IF (obj->QueryValue('0010,1010')) EQ 2 THEN p_age = STRING(obj->GetValue('0010,1010')) ELSE p_age = '040Y'
          IF (obj->QueryValue('0008,0060')) EQ 2 THEN p_modal = STRING(obj->GetValue('0008,0060')) ELSE p_modal = ''
          IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
          IF (obj->QueryValue('0020,0037')) EQ 2 THEN imageOri[i,*] = DOUBLE(obj->GetValue('0020,0037'))

          ;          print, 'Ori', TRANSPOSE(imageOri[i,*])
          ;          print, 'Pos', TRANSPOSE(imagePos[i,*])
          patient_age=FLOAT(STRMID(p_age, 0,3))
          OBJ_DESTROY, obj
          ;缁绢収鍠涢濠氬炊閹冨壖闁告帒妫滄ご鎼佹偝閸ヮ亶鍤㈤柛娆愮墬椤掓粎娑甸敓锟�          print, "voxel size: ", i, ':', voxelsize[i, *]
          IF sliceTh GT 0. THEN BEGIN
            voxelsize[i, 0] = pixelSp[0]
            voxelsize[i, 1] = pixelSp[1]
            IF spacing LE 0. THEN voxelsize[i, 2] = sliceTh + spacing ELSE voxelsize[i, 2] = spacing ;sliceTh +
            voxelsize[i, 3] = rows
            voxelsize[i, 4] = cols
          ENDIF

          voxelsize[i, 5] = sliceLo
        ENDIF
      ENDFOR

      ;double-check the images
      sort_lo = SORT(voxelsize[*, 5])
      diff_lo = SHIFT(voxelsize[*, 5], -1) - voxelsize[*, 5]
      IF ABS(MAX(diff_lo[0:nrfile-2L]) - MIN(diff_lo[0:nrfile-2L])) GT 0.001 THEN BEGIN
        infowarning = DIALOG_MESSAGE('Please sort images first (Please use New -> Sort...)!', /ERROR)
        RETURN
      ENDIF

      ;      ; Create the progress bar.
      progressbar = Obj_New('progressbar', Color='sky blue', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
      ;      ; Place the progress bar on the display.
      progressbar -> Start

      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOMex', a[i])
        result=QUERY_DICOM(a[i], info)

        count=(i+1.)/nrfile*100.0
        progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

        IF result NE 0 THEN BEGIN
          vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
          ;convert grey pixel to HU
          ;rescale intercept
          IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
          ;rescale slope
          IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 1.
          OBJ_DESTROY, obj
          vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)

          ;smooth the image and store to 3D
          vol_HU[*,*,i] = vPixels_HU
        ENDIF
      ENDFOR

      progressbar -> Destroy

      ;鐎电増顨滈崺灞筋灉閿濆洠锟介柣銊ュ閻即宕㈤敓锟�
      sliceDistance = ABS(voxelsize[nrfile/2-1, 5] - voxelsize[nrfile/2, 5])

      ;濠碘�鍊归悘澶愬及椤栨繂澹栭柛蹇撶墣缁绘﹢鏁嶅畝瀣诞闁挎冻璁ｅЧ鏌ユ晬鐏炴儳惟闁搞儲鍎抽崕姘跺磹閹烘洘宕抽柡澶涙嫹
      ;      IF STRCMP(ptPosition, 'FFS') THEN vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)

      ;      IF (voxelsize[nrfile/2, 5] - voxelsize[nrfile/2-1, 5]) LT 0. THEN vol_HU = REVERSE(vol_HU, 3, /OVERWRITE)

      FOV = voxelsize[nrfile/2-1, 3]*voxelsize[nrfile/2-1, 0]
      length = nrfile * sliceDistance
      IF FLOOR(length*256/FOV) LE 256 THEN BEGIN
        vol_HU_iso_resolution = FOV/256.
        vol_HU_iso = CONGRID(vol_HU,  256, 256, FLOOR(length*256/FOV))
      ENDIF ELSE BEGIN
        vol_HU_iso_resolution = length/256.
        vol_HU_iso = CONGRID(vol_HU, FLOOR(FOV*256/length), FLOOR(FOV*256/length), 256)
      ENDELSE

      ;delete vol_HU
      vol_HU = !NULL

      WHILE vol_HU_iso_resolution LT 0.4 DO BEGIN
        dims = size(vol_HU_iso)
        vol_HU_iso_resolution = vol_HU_iso_resolution*2.
        vol_HU_iso = CONGRID(vol_HU_iso, dims[1]/2, dims[2]/2, dims[3]/2)
      ENDWHILE


      ImageOri_temp = imagePos[nrfile/2,*] - imagePos[nrfile/2-1,*]
      ImageOri_z = ImageOri_temp/SQRT(ImageOri_temp[0]^2 + ImageOri_temp[1]^2 + ImageOri_temp[2]^2)

      direction = [[TRANSPOSE(imageOri[nrfile/2,0:2])], [TRANSPOSE(imageOri[nrfile/2,3:5])], [TRANSPOSE(ImageOri_z)]]
      origin_loc = imagePos[0,*]

      ;閻忓繐鏁妎l_HU_iso濠靛鍋勯幊鍡涘炊鐎电硶鏁勯柣褝鎷烽柣鈺冾焾缂嶅绂嶆惔锝傛晞婵﹫鎷稨U=-1000, 闁活枌鍔嶅鐢告偨缁暜aw[0:2]
      max_dim = MAX((SIZE(vol_HU_iso))[1:3], max_ind)
      btm_value = MIN([MIN(vol_HU_iso), -1000])
      vol_HU_cube = DBLARR(max_dim, max_dim, max_dim)*0. + btm_value; -1000.
      IF max_ind EQ 2 THEN BEGIN
        vol_HU_cube[(max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, (max_dim-(SIZE(vol_HU_iso))[1])/2:(max_dim-(SIZE(vol_HU_iso))[1])/2+(SIZE(vol_HU_iso))[1]-1, *] = vol_HU_iso
        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[1])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDIF ELSE BEGIN
        vol_HU_cube[*, *, (max_dim-(SIZE(vol_HU_iso))[3])/2:(max_dim-(SIZE(vol_HU_iso))[3])/2+(SIZE(vol_HU_iso))[3]-1] = vol_HU_iso
        origin_loc = origin_loc - (max_dim-(SIZE(vol_HU_iso))[3])/2 * ImageOri_z * vol_HU_iso_resolution
      ENDELSE

      ;      print, 'min_value = ', MIN(vol_HU_cube)
      ;      print, 'max_value = ', MAX(vol_HU_cube)
      print, 'size(vol_HU_cube) = ', size(vol_HU_cube)


      ;      R_inv = INVERT(curr_ori)
      ;      angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
      ;      angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
      ;      angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
      ;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 0, angle_gamma])
      ;      rot_x = TRANSFORM_VOLUME(rot_z, MISSING = btm_value, Rotation=[angle_beta, 0, 0])
      ;      vol_HU_cube = TRANSFORM_VOLUME(rot_x, MISSING = btm_value, Rotation=[0, 0, angle_alpha])

      ;      print, 'Ori z', TRANSPOSE(ImageOri_z)
      ;      IF imageOri[nrfile/2,4] EQ 0 THEN imageOri[nrfile/2,4] = 0.00001
      ;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, Rotation=[0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
      ;      print, 'rot_z', [0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])]
      ;      IF imageOri[nrfile/2,1] EQ 0 THEN imageOri[nrfile/2,1] = 0.00001
      ;      rot_x = TRANSFORM_VOLUME(rot_z, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0])
      ;      print, 'rot_x', [180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0]
      ;      IF imageOri[nrfile/2,0] EQ 0 THEN imageOri[nrfile/2,0] = 0.00001
      ;      rot_y = TRANSFORM_VOLUME(rot_x, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0])
      ;      print, 'rot_y', [0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0]
      ;
      ;      vol_HU_cube = rot_y
      ;      TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
      ;      vol_HU_cube = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])



      ;delete vol_HU
      vol_HU_iso = !NULL
      ;for black blood
          vol_HU_cube = MAX(vol_HU_cube) - vol_HU_cube


      vol_HU_cube_resolution = vol_HU_iso_resolution
      print, 'resolution', vol_HU_cube_resolution
    ENDIF
  ENDIF ELSE BEGIN
    ;    return,[]
  ENDELSE
END

;GJ, gjia@xidian.edu.cn; transform vol_HU_cube_1 to a new volume based on vol_HU_cube_0
;GJ, 3 steps of fusing two volume cubes based on orientation, center point and resolution
;GJ, Euler angle: https://www.pianshen.com/article/48361722399/
;GJ, 2021/2/23
FUNCTION volume_fusion, vol_HU_cube_0, vol_HU_cube_resolution_0, direction_0, origin_loc_0, ori_vol_HU_cube_1, vol_HU_cube_resolution_1, direction_1, origin_loc_1
  
  ;GJ 2020/4/18, the probram changed ori_vol_HU_cube_1's value. Now fix the bug.
  vol_HU_cube_1 = ori_vol_HU_cube_1
  
  IF MAX(ABS(direction_0-direction_1)) LT 0.001 AND MAX(ABS(origin_loc_0-origin_loc_1)) LT 0.001 THEN RETURN, vol_HU_cube_1

  size_0 = SIZE(vol_HU_cube_0, /dimension)
  size_1 = SIZE(vol_HU_cube_1, /dimension)
  
  ;step 1
  ;GJ, 2021/2/23, adjust the volume size based on cube direction or orientation
  ;avoid the infimite or NaN
  kernal = TRANSPOSE(direction_0)
  R_inv = TRANSPOSE(direction_1) # INVERT(kernal)
  
  ;GJ, 2021/2/23
  ;https://www.pianshen.com/article/48361722399/
  ;z-y-z strictly follow
  angle_alpha = 180./!PI*ATAN(R_inv[1,2], R_inv[0,2])
  angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2), R_inv[2,2])
  angle_gamma = 180./!PI*ATAN(R_inv[2,1], -R_inv[2,0])
  
  IF ABS(angle_beta) LT 0.10 THEN BEGIN
    angle_alpha = 0.
    angle_gamma = 180./!PI*ATAN(-R_inv[0,1], R_inv[0,0])
  ENDIF

  IF ABS(angle_beta-180.) LT 0.10 THEN BEGIN
    angle_alpha = 0.
    angle_gamma = 180./!PI*ATAN(R_inv[0,1], -R_inv[0,0])
  ENDIF
  
  TRANSFORM_VOLUME, vol_HU_cube_1, rot_z, Rotation=[0, 0, angle_alpha]
  TRANSFORM_VOLUME, rot_z, rot_y, Rotation=[0, angle_beta, 0]
  TRANSFORM_VOLUME, rot_y, temp1_vol_HU_cube_1, Rotation=[0, 0, angle_gamma]
  
  ;step 2
  ;GJ, 2021/2/23, adjust the volume size based on center point
  center_loc_0 = origin_loc_0; + TRANSPOSE((direction_0[*,0]*(SIZE_0[0]-1.)/2.0 + direction_0[*,1]*(SIZE_0[1]-1.)/2. + direction_0[*,2]*(SIZE_0[2]-1.)/2.)) * vol_HU_cube_resolution_0
  center_loc_1 = origin_loc_1; + TRANSPOSE((direction_1[*,0]*(SIZE_1[0]-1.)/2.0 + direction_1[*,1]*(SIZE_1[1]-1.)/2. + direction_1[*,2]*(SIZE_1[2]-1.)/2.)) * vol_HU_cube_resolution_1
  center_diff = (center_loc_0 - center_loc_1)/vol_HU_cube_resolution_1
  ;  center_diff[2] = center_diff[2]+8.
  print, '1st center loc = ', center_loc_0
  print, '2nd center loc = ', center_loc_1

  TRANSFORM_VOLUME, temp1_vol_HU_cube_1, temp2_vol_HU_cube_1, Translate=center_diff
  
  ;step 3 
  ;GJ, 2021/2/22, adjust the volume size based on resolution
  print, 'resolution0:', vol_HU_cube_resolution_0
  print, 'resolution1:', vol_HU_cube_resolution_1
  ave_resolution = 0.5 * (vol_HU_cube_resolution_0 + vol_HU_cube_resolution_1)
  IF ABS(vol_HU_cube_resolution_0-vol_HU_cube_resolution_1) GT 0.01*ave_resolution THEN BEGIN
    FOV_0 = vol_HU_cube_resolution_0*size_0[0]
    FOV_1 = vol_HU_cube_resolution_1*size_1[0]
    new_size_1 = FLOOR(FOV_1/vol_HU_cube_resolution_0)
    temp3_vol_HU_cube_1 = CONGRID(temp2_vol_HU_cube_1, new_size_1, new_size_1, new_size_1)
    IF new_size_1 GT size_0[0] THEN BEGIN
      size_diff_low = FLOOR((new_size_1 - size_0[0])/2.)
      new_vol_HU_cube_2=temp3_vol_HU_cube_1[size_diff_low:(size_diff_low+size_0[0]-1), size_diff_low:(size_diff_low+size_0[0]-1), size_diff_low:(size_diff_low+size_0[0]-1)]
    ENDIF ELSE BEGIN
      size_diff_low = FLOOR(-1. * (new_size_1 - size_0[0])/2.)
      new_vol_HU_cube_2=vol_HU_cube_0*0. + MIN(temp3_vol_HU_cube_1)
      new_vol_HU_cube_2[size_diff_low:(size_diff_low+new_size_1-1), size_diff_low:(size_diff_low+new_size_1-1), size_diff_low:(size_diff_low+new_size_1-1)] = temp3_vol_HU_cube_1
    ENDELSE
  ENDIF ELSE BEGIN
    new_vol_HU_cube_2 = temp2_vol_HU_cube_1
  ENDELSE

  
  RETURN, new_vol_HU_cube_2
END


FUNCTION volume_fusion_old, vol_HU_cube_0, vol_HU_cube_resolution_0, direction_0, origin_loc_0, vol_HU_cube_1, vol_HU_cube_resolution_1, direction_1, origin_loc_1

new_vol_HU_cube_1 = vol_HU_cube_1*0. + MIN(vol_HU_cube_1)
size_0 = SIZE(vol_HU_cube_0, /dimension)
size_1 = SIZE(vol_HU_cube_1, /dimension)

FOR i=0, size_0[0]-1L DO BEGIN
  FOR j=0, size_0[1]-1L DO BEGIN
    FOR k=0, size_0[2]-1L DO BEGIN
      loc_0 = origin_loc_0 + TRANSPOSE((direction_0[*,0]*i + direction_0[*,1]*(size_0[1]-1-j) + direction_0[*,2]*k)) * vol_HU_cube_resolution_0
      
      kernal = TRANSPOSE(direction_1)
      loc_index = (loc_0 - origin_loc_1)/vol_HU_cube_resolution_1 # INVERT(kernal)
      loc_index[1] = size_1[1]-1-loc_index[1]
      
      loc_1 = origin_loc_1 + TRANSPOSE((direction_1[*,0]*loc_index[0] + direction_1[*,1]*(size_1[1]-1-loc_index[1]) + direction_1[*,2]*loc_index[2])) * vol_HU_cube_resolution_1
      
      IF loc_index[0] GT 0 AND loc_index[0] LT size_1[0]-1 THEN BEGIN
        IF loc_index[1] GT 0 AND loc_index[1] LT size_1[1]-1 THEN BEGIN
          IF loc_index[2] GT 0 AND loc_index[2] LT size_1[2]-1 THEN BEGIN
            new_vol_HU_cube_1[i, size_0[1]-1-j, k] = vol_HU_cube_1[loc_index[0], size_1[1]-1-loc_index[1], loc_index[2]]
          ENDIF
        ENDIF
      ENDIF
    ENDFOR
  ENDFOR
ENDFOR

RETURN, new_vol_HU_cube_1
END



;+
; NAME:
;       TRANSFORM_VOLUME
;
; PURPOSE:
;
;       The purpose of this program is to transform (e.g., rotate,
;       scale, and translate) a 3D array or volume.
;
; AUTHOR:
;
;       Martin Downing,
;       Clinical Research Physicist,
;       Grampian Orthopaedic RSA Research Centre,
;       Woodend Hospital, Aberdeen, AB15 6LS.
;       Pnone: 01224 556055 / 07903901612
;       Fa: 01224 556662
;       E-mail: m.downing@abdn.ac.uk
;
; CATEGORY:
;
;      Mathematics, graphics.
;
; CALLING SEQUENCE:
;
;      result = TRANSFORM_VOLUME( volume )
;
; INPUTS:
;
;       volume:    The 3D array or volume to be transformed.
;
; OPTIONAL KEYWORDS:
;
;      BUFFER_SIZE: To reduce memory overhead the routine processes the job in chunks, the number
;         of elements of which can be set using the BUFFER_SIZE keyword, set this keyword to
;         0 to force the whole array to be processed at one time. The default value is 128.
;
;      MISSING: The value to return for transformed values outside the bounds of
;         the volume. (Passed to the INTERPOLATE function.) Default is 0.
;
;      T3DMAT: The homogeneous transforamtion matrix. If this keyword is not present,
;         the following keywords can be used to create a homogeneous transformation matrix:
;
;         ROTATION - The rotation vector [rx,ry,rz]. The order of rotation is ZYX.
;         TRANSLATE - The translation vector [tx,ty,tz].
;         SCALE - The scale vector [sx,sy,sz].
;         CENTRE_ROTATION - The centre of rotation [cx,cy,cz].
;
; OUTPUTS:
;
;       result:    The transformed array or volume.
;
; COMMON BLOCKS:
;
;       None.
;
; DEPENDENCIES:
;
;       The program uses the library INTERPLOLATE routine, which currently (IDL 5.4)
;       uses linear interpolation. Note that the operation is performed in chunks,
;       each of which is independant of the result of the others, so the operation
;       could easiliy be parallelised.
;
; MODIFICATION HISTORY:
;
;       Written by: Martin Downing, 16 September 2001.
;       Added MISSING keyword. Removed INPLACE keyword. 25 Nov 2001. MD
;-
;******************************************************************************************;
;  Copyright (c) 2008, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
PRO Transform_Volume, volume, new_vol, Rotation=rotation, $
  Scale=scale, Translate=translate, Centre_Rotation=centre_rotation, $
  T3Dmat=t3dmat, Buffer_Size=buffer_size, Missing=missing

  ; Error handling.

  Catch, theError
  IF theError NE 0 THEN BEGIN
    Catch, /Cancel
    ok = Dialog_Message(!Error.State.Msg)
    RETURN
  ENDIF

  ; Find the dimensions of the volume.

  s = Size(volume)
  sx=s[1] & sy=s[2] & sz=s[3]
  st = sx*sy*sz

  vol_t = volume
  IF N_Elements(missing) THEN missing = 0

  ; Create a transform matrix, if one is not provided.

  IF N_Elements(t3dmat) EQ 0 THEN begin

    IF N_Elements(rotation) EQ 0 THEN rotation =[0,0,0]
    IF N_Elements(centre_rotation) EQ 0  THEN centre_rotation=[(sx-1)/2.0,(sy-1)/2.0,(sz-1)/2.0]
    IF N_Elements(translate) EQ 0 THEN translate =[0,0,0]
    IF N_Elements(scale) EQ 0 THEN scale =[1,1,1]

    T3D, /Reset, Translate = -centre_rotation
    T3D, Rotate=rotation
    T3D, Translate= centre_rotation + translate, Scale=scale
    t3dmat = !P.T

  ENDIF

  ; Check buffer size. The size 128 is optimim on my system, You may
  ; want to try other values.

  IF N_Elements(buffer_size) EQ 0 THEN buffer_size = 256
  IF buffer_size LE 0 THEN buffer_size = st

  ; Perform the transformations.

  FOR j=0L,(st-1),buffer_size DO BEGIN

    ; Account for possible odd last chunk.

    bufsize = buffer_size < (st-j)

    ; Generate volume coordinates by interpolating temporary array of volume indices.

    i = j + Lindgen(bufsize)
    coords = [ [(i MOD sx)],[((i / sx) MOD (sy))], [(i / (sx * sy))], [Replicate(1b, bufsize)]]
    coords = Temporary(coords) # t3dmat
    vol_t[j:j+bufsize-1] = Interpolate(volume, coords[*,0], coords[*,1], coords[*,2], Missing=missing)
  ENDFOR

  ; Return the transformed volume.

  new_vol = vol_t
END


;---------------------------------------------------------
PRO CINE

  ; Load the data.
  head = BytArr(80,100,57)
  file = Filepath(Subdirectory=['examples','data'], 'head.dat')
  OpenR, lun, file, /Get_Lun
  ReadU, lun, head
  Free_Lun, lun

  ; Set up the animation.
  XInterAnimate, Set=[240, 300, 37], /Showload

  ; Load the animation
  FOR j=0,36 DO BEGIN
    Transform_Volume, head, rotHead, Missing=0, Rotation=[0,0,(j*10) MOD 360]
    print, [0,0,(j*10) MOD 360]
    TVScl, Reform(Max(rotHead, DIMENSION=1))
    XInteranimate, Frame=j, Window=!D.Window
  ENDFOR
  XInterAnimate, 50

END
;---------------------------------------------------------



pro writestl,outputfile,vers,conn
  file=outputfile
  
  IF STRCMP(STRMID(file, STRLEN(file)-4, 4), '.stl') NE 1 THEN file = file + '.stl'
  numberVertices = MESH_DECIMATE(vers,conn, p2, VERTICES = v2)
  
  ;delete the variables
  vers = !NULL
  conn = !NULL
  ;  arr=bytarr(5,5,5)
  ;  arr[1,1,1]=1
  ;  arr[1,3,1]=1
  ;  arr[3,1,1]=1
  ;  arr[3,3,1]=1
  ;  arr[2,2,2]=1

  ;endfor
  ;
  ; rgb=[255,1,1,255,134,1,255,255,1,1,255,25,1,1,255,92,1,255,255,1,243]
  ;afn = 'C:\temp\MED_vis\S0004_130952_Recon3\'
  ;afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  ;IF STRLEN(afn) NE 0 THEN BEGIn
  ;  a = FILE_SEARCH(afn, '*.dcm', count=nrfile)
  ;  IF nrfile NE 0 THEN BEGIN
  ;    obj = OBJ_NEW('IDLffDICOMex', a[0], /NO_PIXEL_DATA)
  ;    obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
  ;    OBJ_DESTROY, obj
  ;    vol_HU = DBLARR(rows, cols, nrfile)
  ;    voxelsize = DBLARR(nrfile, 5) * 0.
  ;    FOR i=0L, nrfile-1L DO BEGIN       obj = OBJ_NEW('IDLffDICOMex', a[i])
  ;    result=QUERY_DICOM(a[i], info)
  ;    IF result NE 0 THEN BEGIN
  ;      obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
  ;      order=1
  ;      vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
  ;
  ;      IF (obj->QueryValue('0018,0088')) EQ 2 THEN spacing = obj->GetValue('0018,0088') ELSE spacing = 0.
  ;      IF (obj->QueryValue('0018,0050')) EQ 2 THEN sliceTh = obj->GetValue('0018,0050') ELSE sliceTh = 0.
  ;      IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = obj->GetValue('0028,0030') ELSE pixelSp = 0.
  ;      IF sliceTh GT 0. THEN BEGIN
  ;        voxelsize[i, 0] = pixelSp[0]
  ;        voxelsize[i, 1] = pixelSp[1]
  ;       ; voxelsize[i, 2] = sliceTh + spacing
  ;        voxelsize[i, 2] = sliceTh
  ;        voxelsize[i, 3] = rows
  ;        voxelsize[i, 4] = cols
  ;      ENDIF
  ;
  ;      ;convert grey pixel to HU
  ;      ;rescale intercept
  ;      IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
  ;      ;rescale slope
  ;      IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 0.
  ;      OBJ_DESTROY, obj
  ;      vPixels_HU = vPixels * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)
  ;
  ;      ;smooth the image and store to 3D
  ;      vol_HU[*,*,i] = vPixels_HU
  ;    ENDIF
  ;  ENDFOR
  ;
  ;  vol_HU_iso = CONGRID(vol_HU,  FLOOR(voxelsize[0, 3]*voxelsize[0, 0]/voxelsize[0, 2]), FLOOR(voxelsize[0, 4]*voxelsize[0, 1]/voxelsize[0, 2]), nrfile)
  ;
  ;ENDIF
  ;ENDIF

  ; SHADE_VOLUME, arr,0.9, v1, p1

  ;  help,v1
  ;  help,p1
  ;  numberVertices = MESH_DECIMATE(v1, p1, p2, VERTICES = v2)
  ;  help,v2
  ;  help,p2



  ;facenum=long(n_elements(p2)/4)
  facenum=long(n_elements(p2)/4)
  print,facenum
  normalandver=fltarr(3,facenum*4)

  normal=fltarr(3,1)
  point1=fltarr(3,1)
  point2=fltarr(3,1)
  point3=fltarr(3,1)
  ;print,'p2',p2
  ;print,'v2',v2

  p=0

  i=long(0)

  for j=0,facenum-1 do begin
    ;for j=long(0),1 do begin
    ;j=0


    order=long(p2[i+1:i+3] )
;    print,'order',order
    point1=v2[*,order[0]]
    point2=v2[*,order[1]]
    point3=v2[*,order[2]]
    ;    a= ( (point2[1]-point1[1])*(point3[2]-point1[2])-(point2[2]-point1[2])*(point3[1]-point1[1]) );
    ;    b = ( (point2[2]-point1[2])*(point3[0]-point1[0])-(point2[0]-point1[0])*(point3[2]-point1[2]) );
    ;    c = ( (point2[0]-point1[0])*(point3[1]-point1[1])-(point2[1]-point1[1])*(point3[0]-point1[0]) );

    ;;;;;;;;hxn20180530;;;;
    a= ( (point1[1]-point3[1])*(point2[2]-point3[2])-(point1[2]-point3[2])*(point2[1]-point3[1]) );
    b = ( (point1[2]-point3[2])*(point2[0]-point3[0])-(point2[2]-point3[2])*(point1[0]-point3[0]) );
    c = ( (point1[0]-point3[0])*(point2[1]-point3[1])-(point2[0]-point3[0])*(point1[1]-point3[1]) );



    ;;;;;;;;;;;;;;
    model=sqrt(a*a+b*b+c*c)
    normalandver[*,i]  =[a,b,c]/model
    normalandver[*,i+1]=point1
    normalandver[*,i+2]=point2
    normalandver[*,i+3]=point3
    i=i+4

    ; order=p2[i+1:i+3]
    ; print,'order',order
    ; point1=v2[*,order[0]]
    ; point2=v2[*,order[1]]
    ; point3=v2[*,order[2]]
    ; a= ( (point2[1]-point1[1])*(point3[2]-point1[2])-(point2[2]-point1[2])*(point3[1]-point1[1]) );
    ; b = ( (point2[2]-point1[2])*(point3[0]-point1[0])-(point2[0]-point1[0])*(point3[2]-point1[2]) );
    ; c = ( (point2[0]-point1[0])*(point3[1]-point1[1])-(point2[1]-point1[1])*(point3[0]-point1[0]) );
    ; normalandver[*,i]  =[a,b,c]
    ; normalandver[*,i+1]=point1
    ; normalandver[*,i+2]=point2
    ; normalandver[*,i+3]=point3



  endfor

  filename=strarr(80)
  filename[0]='l'
  filename[1]='i'
  filename[2]='u'
  filename[3:79]=' '
  sou="liubo"
  strput,filename,sou,1

  Attribute=0u


  ;
  oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], v2, POLYGON=p2)
  mymodel1 = OBJ_NEW('IDLgrModel')
  mymodel1->Add, oPoly1
  ;
  xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]
  
;  ;GJ 2019/8/12, test texture mapping
;  imagefile='D:\AIMIS_3D\gj.jpg'
;  read_jpeg, imagefile, image
;  oimage = OBJ_NEW('IDLgrImage', image, interleave=0, /interpolate)
;  n=size(v2,/dim)
;  oPoly2 = OBJ_NEW('IDLgrPolygon', v2, POLYGON=p2, COLOR = [255, 127, 127], texture_coord=FINDGEN(2, n[1])/n[1]*4., texture_map=oimage, /texture_interp)
;  mymodel2 = OBJ_NEW('IDLgrModel')
;  mymodel2->Add, oPoly2
;  xobjview_XDU, mymodel2, BACKGROUND = [0, 0, 0]
  
  
  openw,lun,file,/get_lun
  writeu,lun,filename

  writeu,lun,facenum
  i=long(0)
  for k=0,facenum-1 do begin

    writeu,lun,normalandver[*,i]
    writeu,lun,normalandver[*,i+1]
    writeu,lun,normalandver[*,i+2]
    writeu,lun,normalandver[*,i+3]
    writeu,lun,Attribute

    i=i+4

  endfor

  free_lun,lun


  q=0

end


;;GJ, 2018/10/08
;;Function name: load_to_stl
;Input: dcm_dir: dicom file directory
;       stl_fn: stl filename for stl file to be saved
;打开文件：
;传入变量： dcm文件夹名（字符串） 生成的stl存储文件夹名（字符串）
;函数功能：根据dcm文件，生成对应stl文件，存储在对应stl文件夹中。
;dcm_dir = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
;stl_fn = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\test.stl'
PRO load_dcm_to_stl, dcm_dir, stl_fn

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN dcm_dir = 'C:\Temp\Research\Med_visualization\4330003_LIU_YING_20180907_HR_MRA\S0201_165711_s3DI_MC_HR\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
  IF N_ELEMENTS(stl_fn) EQ 0 THEN stl_fn = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\test.stl'
  IF STRLEN(dcm_dir) NE 0 THEN BEGIN
    a = FILE_SEARCH(dcm_dir, '*.dcm', count=nrfile)
    IF nrfile EQ 0 THEN RETURN

    image_read, dcm_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
    
    print, 'direction =', direction
    print, 'origin_loc =', origin_loc
;    
;    kernal = TRANSPOSE(direction)
;    R_inv = INVERT(kernal)
;    angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;    angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;    angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;    rot_z = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 0, angle_gamma])
;    rot_x = TRANSFORM_VOLUME(rot_z, MISSING = btm_value, Rotation=[angle_beta, 0, 0])
;    vol_HU_cube = TRANSFORM_VOLUME(rot_x, MISSING = btm_value, Rotation=[0, 0, angle_alpha])
;    
;    size = SIZE(vol_HU_cube, /dimension)
;    loc = origin_loc + TRANSPOSE((direction[*,0]*(SIZE[0]/2-1) + direction[*,1]*(SIZE[1]/2-1) + direction[*,2]*(SIZE[2]/2-1))) * vol_HU_cube_resolution
;    new_loc = R_inv * loc

    
    ;      print, 'Ori z', TRANSPOSE(ImageOri_z)
    ;      IF imageOri[nrfile/2,4] EQ 0 THEN imageOri[nrfile/2,4] = 0.00001
    ;      rot_z = TRANSFORM_VOLUME(vol_HU_cube, Rotation=[0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
    ;      print, 'rot_z', [0, 0, 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])]
    ;      IF imageOri[nrfile/2,1] EQ 0 THEN imageOri[nrfile/2,1] = 0.00001
    ;      rot_x = TRANSFORM_VOLUME(rot_z, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0])
    ;      print, 'rot_x', [180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 0, 0]
    ;      IF imageOri[nrfile/2,0] EQ 0 THEN imageOri[nrfile/2,0] = 0.00001
    ;      rot_y = TRANSFORM_VOLUME(rot_x, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0])
    ;      print, 'rot_y', [0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 0]
    ;
    ;      vol_HU_cube = rot_y
    ;      TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[180./!PI*ATAN(ImageOri_z[1]/ImageOri_z[2]), 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])
    ;      vol_HU_cube = TRANSFORM_VOLUME(vol_HU_cube, MISSING = btm_value, Rotation=[0, 180./!PI*ATAN(imageOri[nrfile/2,2]/imageOri[nrfile/2,0]), 180./!PI*ATAN(imageOri[nrfile/2,3]/imageOri[nrfile/2,4])])


    
    dcm_dir1 = 'C:\Temp\Research\Med_visualization\4330003_LIU_YING_20180907_HR_MRA\S0401_170837_3D_T1WIVesselWall\'
    image_read, dcm_dir1, vol_HU_cube1, vol_HU_cube_resolution1, patient_age1, direction1, origin_loc1
    print, 'direction1 =', direction1
    print, 'origin_loc1 =', origin_loc1
;    
;    kernal = TRANSPOSE(direction1)
;    R_inv = INVERT(kernal)
;    angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;    angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;    angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;    rot_z = TRANSFORM_VOLUME(vol_HU_cube1, MISSING = btm_value, Rotation=[0, 0, angle_gamma])
;    rot_x = TRANSFORM_VOLUME(rot_z, MISSING = btm_value, Rotation=[angle_beta, 0, 0])
;    vol_HU_cube1 = TRANSFORM_VOLUME(rot_x, MISSING = btm_value, Rotation=[0, 0, angle_alpha])
;    
;    size1 = SIZE(vol_HU_cube1, /dimension)
;    loc1 = origin_loc1 + TRANSPOSE((direction1[*,0]*(SIZE1[0]/2-1) + direction1[*,1]*(SIZE1[1]/2-1) + direction1[*,2]*(SIZE1[2]/2-1))) * vol_HU_cube_resolution1
;    new_loc1 = R_inv * loc1

    
    vol_HU_cube_new1 = volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
    
    dcm_dir2 = 'C:\Temp\Research\Med_visualization\4330003_LIU_YING_20180907_HR_MRA\S1001_173516_3D_T1WIVesselWall\'
    image_read, dcm_dir2, vol_HU_cube2, vol_HU_cube_resolution2, patient_age2, direction2, origin_loc2
    print, 'direction2 =', direction2
    print, 'origin_loc2 =', origin_loc2

;    kernal = TRANSPOSE(direction2)
;    R_inv = INVERT(kernal)
;    angle_alpha = -180./!PI*ATAN(R_inv[0,2]/R_inv[1,2])
;    angle_beta = 180./!PI*ATAN(SQRT(R_inv[2,0]^2+R_inv[2,1]^2)/R_inv[2,2])
;    angle_gamma = 180./!PI*ATAN(R_inv[2,0]/R_inv[2,1])
;    rot_z = TRANSFORM_VOLUME(vol_HU_cube2, MISSING = btm_value, Rotation=[0, 0, angle_gamma])
;    rot_x = TRANSFORM_VOLUME(rot_z, MISSING = btm_value, Rotation=[angle_beta, 0, 0])
;    vol_HU_cube2 = TRANSFORM_VOLUME(rot_x, MISSING = btm_value, Rotation=[0, 0, angle_alpha])

    vol_HU_cube_new2 = volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, vol_HU_cube2, vol_HU_cube_resolution2, direction2, origin_loc2)



    sizeVHc = SIZE(vol_HU_cube)
    mean_image = IMAGE_THRESHOLD(BYTE(vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MEAN, /INVERT)
    threshold_min_value = tf[0];*0.5 ;GJ 2019/5/19
    print, 'threshold_min_value = ', threshold_min_value
    threshold_max_value = MAX(vol_HU_cube)
    print, 'threshold_max_value = ', threshold_max_value
    vol_HU_cube_min = MIN(vol_HU_cube)
    print, 'size(vol_HU_cube) = ', size(vol_HU_cube)

    vol_HU_cube_mask = (vol_HU_cube GE threshold_min_value) * (vol_HU_cube LE threshold_max_value)
    btm_value = MIN([MIN(vol_HU_cube), -1000])
    vol_HU_cube =  vol_HU_cube_mask * (threshold_min_value+1-btm_value) +btm_value
    
    ;delete the variable
    vol_HU_cube_mask = !NULL

    ;
    vol_HU_cube[0, *, *] = vol_HU_cube[0, *, *]*0. +btm_value
    vol_HU_cube[(size(vol_HU_cube))[1]-1, *, *] = vol_HU_cube[(size(vol_HU_cube))[1]-1, *, *]*0. +btm_value
    vol_HU_cube[*, 0, *] = vol_HU_cube[*, 0, *]*0. +btm_value
    vol_HU_cube[*, (size(vol_HU_cube))[2]-1, *] = vol_HU_cube[*, (size(vol_HU_cube))[2]-1, *]*0. +btm_value
    vol_HU_cube[*, *, 0] = vol_HU_cube[*, *, 0]*0. +btm_value
    vol_HU_cube[*, *, (size(vol_HU_cube))[3]-1] = vol_HU_cube[*, *, (size(vol_HU_cube))[3]-1]*0. +btm_value
    SHADE_VOLUME, vol_HU_cube,threshold_min_value, Outverts, Outconn1
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1);, ITERATIONS=200)   ;smooth the segmented volume; GJ 2017/12/22
    
    ;delete the vairables
    vol_HU_cube = !NULL
    Outverts = !NULL
;    oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], smoothedOutverts, POLYGON=Outconn1)
;    mymodel1 = OBJ_NEW('IDLgrModel')
;    mymodel1->Add, oPoly1
    ;
;    xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]
    

;    IF STRCMP(STRMID(stl_fn, STRLEN(stl_fn)-4, 4), '.stl') NE 1 THEN stl_fn = stl_fn + '.stl'
;    numberVertices = MESH_DECIMATE(smoothedOutverts, Outconn1, p2, VERTICES = v2)
    writestl,stl_fn,vol_HU_cube_resolution1, Outconn1
    
    ;delete the variables
    smoothedOutverts  = !NULL
    Outconn1 = !NULL
    
  ENDIF

END

;;GJ, 2019/12/03

FUNCTION volume_reslice, reslice2d, reslice2d_resolution_0, direction_0, origin_loc_0, vol_HU_cube_1, vol_HU_cube_resolution_1, direction_1, origin_loc_1

  new_reslice2d = reslice2d*0. + MIN(reslice2d)
  size_0 = SIZE(reslice2d, /dimension)
  size_1 = SIZE(vol_HU_cube_1, /dimension)

  FOR i=0, size_0[0]-1L DO BEGIN
    FOR j=0, size_0[1]-1L DO BEGIN
      k=0.
      loc_0 = origin_loc_0 + TRANSPOSE((direction_0[*,0]*i + direction_0[*,1]*(size_0[1]-1-j) + direction_0[*,2]*k)) * reslice2d_resolution_0

      kernal = TRANSPOSE(direction_1)
      loc_index = (loc_0 - origin_loc_1)/vol_HU_cube_resolution_1 # INVERT(kernal)
      loc_index[1] = size_1[1]-1-loc_index[1]

      loc_1 = origin_loc_1 + TRANSPOSE((direction_1[*,0]*loc_index[0] + direction_1[*,1]*(size_1[1]-1-loc_index[1]) + direction_1[*,2]*loc_index[2])) * vol_HU_cube_resolution_1

      IF loc_index[0] GT 0 AND loc_index[0] LT size_1[0]-1 THEN BEGIN
        IF loc_index[1] GT 0 AND loc_index[1] LT size_1[1]-1 THEN BEGIN
          IF loc_index[2] GT 0 AND loc_index[2] LT size_1[2]-1 THEN BEGIN
            new_reslice2d[i, size_0[1]-1-j] = vol_HU_cube_1[loc_index[0], size_1[1]-1-loc_index[1], loc_index[2]]
          ENDIF
        ENDIF
      ENDIF
    ENDFOR
  ENDFOR

  RETURN, new_reslice2d
END

;;GJ, 2019/12/03
;;Function name: reslice_dcm_volume
;Input: ori_dcm_dir: original dicom file directory
;       reslice_dcm_dir: resliced image directory based on these dicom images
;打开文件：
;传入变量： 原始dcm文件夹名（字符串）
;      重新切片的dcm文件夹名（字符串）
;函数功能： 生成的按照新dcm文件生产的图片。
;先选择亮血序列文件夹
;再选择黑血序列文件夹
;在黑血序列文件夹中自动保存dicom文件同名的jpeg文件
;dcm_dir = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
;stl_fn = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\test.stl'
PRO reslice_dcm_volume, dcm_dir, reslice_dcm_dir

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    ;    dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\2163399_ZHU_KUI_HE_20170824_MCA\S0301_173624_s3DI_MC_HR\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  IF N_ELEMENTS(reslice_dcm_dir) EQ 0 THEN BEGIN
    ;    reslice_dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\2163399_ZHU_KUI_HE_20170824_MCA\S0901_180949_3D_T1WIVesselWall\'
    reslice_dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select reslice2d dcm folder:')
    filearr=file_search(reslice_dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  IF STRLEN(reslice_dcm_dir) NE 0 THEN BEGIN
    a = FILE_SEARCH(reslice_dcm_dir, '*.dcm', count=nrfile)
    IF nrfile EQ 0 THEN RETURN

    image_read, dcm_dir, vol_HU_cube1, vol_HU_cube_resolution1, patient_age1, direction1, origin_loc1

    print, 'direction =', direction1
    print, 'origin_loc =', origin_loc1

    image_read, reslice_dcm_dir, vol_HU_cube0, vol_HU_cube_resolution, patient_age, direction, origin_loc0
    print, 'direction1 =', direction
    print, 'origin_loc1 =', origin_loc0

    ;      ; Create the progress bar.
    progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    ;      ; Place the progress bar on the display.
    progressbar -> Start

    imagePos = DBLARR(nrfile, 3) * 0.
    FOR i=0L, nrfile-1L DO BEGIN
      count=(i+1.)/nrfile*100.0
      progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

      obj = OBJ_NEW('IDLffDICOMex', a[i])
      result=QUERY_DICOM(a[i], info)
      IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
      IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
      OBJ_DESTROY, obj

      origin_loc = imagePos[i,*]
      reslice2d_resolution = pixelSp[0]
      reslice2d=read_dicom(a[i])
      new_reslice2d = volume_reslice(reslice2d, reslice2d_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
      WRITE_JPEG, a[i]+'.jpeg', BYTSCL(new_reslice2d), QUALITY=200
    ENDFOR

    progressbar -> Destroy
    ;    i=45
    ;    origin_loc = imagePos[i,*]
    ;    reslice2d_resolution = pixelSp[0]
    ;    reslice2d=read_dicom(a[i])
    ;    iimage, reslice2d
    ;
    ;    new_reslice2d = volume_reslice(reslice2d, reslice2d_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
    ;
    ;    iimage, new_reslice2d

  ENDIF

END

;;GJ, 2019/12/03
;;Function name: reslice_dcm_volume
;Input: ori_dcm_dir: original dicom file directory
;       reslice_dcm_dir: resliced image directory based on these dicom images
;打开文件：
;传入变量： 原始dcm文件夹名（字符串）
;      重新切片的dcm文件夹名（字符串）
;函数功能： 生成的按照新dcm文件生产的图片。
;先选择亮血序列文件夹
;再选择黑血序列文件夹
;在黑血序列文件夹中自动保存dicom文件同名的jpeg文件
;dcm_dir = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
;stl_fn = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\test.stl'
PRO reslice_dcm_volume_DCE, dcm_dir, reslice_dcm_dir

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    
    dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P01\B01_20090708\Dyn01\S2201_084710_DYNAMIC_FFE\V06_\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    ;dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P04\B01_20090817\Dyn01\S2301_083914_DYNAMIC_FFE\V09_\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    ;dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P05_neoadjuvant\B01_20090824\Dyn01\Dyn01_images\V06_\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    ;    dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  IF N_ELEMENTS(reslice_dcm_dir) EQ 0 THEN BEGIN
    reslice_dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P01\B01_20090708\raw_data\S1101_082126_T2_TSE_ax\'
    ;reslice_dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P04\B01_20090817\raw_data\Raw_data\S0801_075142_T2_TSE_ax\'
;    reslice_dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select reslice2d dcm folder:')
    filearr=file_search(reslice_dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  IF STRLEN(reslice_dcm_dir) NE 0 THEN BEGIN
    a = FILE_SEARCH(reslice_dcm_dir, '*.dcm', count=nrfile)
    IF nrfile EQ 0 THEN RETURN

    image_read_512_DCE, dcm_dir, vol_HU_cube1, vol_HU_cube_resolution1, patient_age1, direction1, origin_loc1

    print, 'direction =', direction1
    print, 'origin_loc =', origin_loc1

    image_read_512, reslice_dcm_dir, vol_HU_cube0, vol_HU_cube_resolution, patient_age, direction, origin_loc0
    print, 'direction1 =', direction
    print, 'origin_loc1 =', origin_loc0

    ;      ; Create the progress bar.
    progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    ;      ; Place the progress bar on the display.
    progressbar -> Start

    imagePos = DBLARR(nrfile, 3) * 0.
    FOR i=0L, nrfile-1L DO BEGIN
      count=(i+1.)/nrfile*100.0
      progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

      obj = OBJ_NEW('IDLffDICOMex', a[i])
      result=QUERY_DICOM(a[i], info)
      IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
      IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
      OBJ_DESTROY, obj

      origin_loc = imagePos[i,*]
      reslice2d_resolution = pixelSp[0]
      reslice2d=read_dicom(a[i])
      new_reslice2d = volume_reslice(reslice2d, reslice2d_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
      iF i LT 9 THEN BEGIN
        unet_file = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P01\B01_20090708\P01_unet\outline\outline'+STRING(i+41, format='(I2)')+'.png'
      ENDIF ELSE BEGIN
        unet_file = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P01\B01_20090708\P01_unet\outline\outline'+STRING(i+41, format='(I2)')+'.png'
      ENDELSE
;      iF i LT 9 THEN BEGIN
;        unet_file = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P04\B01_20090817\P04_unet\outline\outline'+STRING(i+1, format='(I1)')+'.png'
;      ENDIF ELSE BEGIN
;        unet_file = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P04\B01_20090817\P04_unet\outline\outline'+STRING(i+1, format='(I2)')+'.png'
;      ENDELSE
      READ_PNG, unet_file, image_png
      temp_unet_image = REVERSE(image_png, 2)
      
      WRITE_JPEG, a[i]+'.jpeg', BYTSCL(new_reslice2d), QUALITY=200
      new_reslice2d[where(temp_unet_image LT 20, /NULL)] = 0

      IF i NE -1 THEN BEGIN
;      IF i EQ 19 THEN BEGIN
;        iimage, new_reslice2d
;        iimage, reslice2d
        ;blend display the result
        s1=new_reslice2d/100.
        s1(where(S1 LT 30)) = 0;s1(where(S1 LT 300)) = 0
        
        s0=reslice2d
        backgroundImage=DOUBLE(S0)/DOUBLE(max(S0))*255.
        foregroundImage1=s1/200.*255.;foregroundImage1=s1/2000.*255.

        major=5
        max_ct=200;max_ct=2000
        min_ct=0
        filename1=a[i]+'_auc.jpeg';'unet.jpg'
        closeYes=0
        Image_Blend, REVERSE(backgroundImage,2), REVERSE(foregroundImage1,2), COLORTABLE=33, blendTitle='T2_AUC Image', major, max_ct, min_ct, filename1, closeYes, 'test'


      ENDIF
    ENDFOR

    progressbar -> Destroy
    ;    i=45
    ;    origin_loc = imagePos[i,*]
    ;    reslice2d_resolution = pixelSp[0]
    ;    reslice2d=read_dicom(a[i])
    ;    iimage, reslice2d
    ;
    ;    new_reslice2d = volume_reslice(reslice2d, reslice2d_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
    ;
    ;    iimage, new_reslice2d

  ENDIF

END


;;GJ, 2019/12/03
;;Function name: reslice_dcm_volume
;Input: ori_dcm_dir: original dicom file directory
;       reslice_dcm_dir: resliced image directory based on these dicom images
;打开文件：
;传入变量： 原始dcm文件夹名（字符串）
;      重新切片的dcm文件夹名（字符串）
;函数功能： 生成的按照新dcm文件生产的图片。
;先选择亮血序列文件夹
;再选择黑血序列文件夹
;在黑血序列文件夹中自动保存dicom文件同名的jpeg文件
;dcm_dir = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
;stl_fn = 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\test.stl'
PRO reslice_dcm_volume_DWI, dcm_dir, reslice_dcm_dir

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P06\B01_20090825\Raw_data\S1202_140542_dWIPDWI4bNSA3CLEAR\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
;    dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF
  
  IF N_ELEMENTS(reslice_dcm_dir) EQ 0 THEN BEGIN
    reslice_dcm_dir = 'C:\D_drive\Bladder_UnetDCE\Bladderpatient\P06\B01_20090825\Raw_data\S1001_135126_T2_TSE_ax\'
;    reslice_dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select reslice2d dcm folder:')
    filearr=file_search(reslice_dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF
  
  IF STRLEN(reslice_dcm_dir) NE 0 THEN BEGIN
    a = FILE_SEARCH(reslice_dcm_dir, '*.dcm', count=nrfile)
    IF nrfile EQ 0 THEN RETURN

    image_read_512, dcm_dir, vol_HU_cube1, vol_HU_cube_resolution1, patient_age1, direction1, origin_loc1

    print, 'direction =', direction1
    print, 'origin_loc =', origin_loc1

    image_read_512, reslice_dcm_dir, vol_HU_cube0, vol_HU_cube_resolution, patient_age, direction, origin_loc0
    print, 'direction1 =', direction
    print, 'origin_loc1 =', origin_loc0
    
    ;      ; Create the progress bar.
    progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
    ;      ; Place the progress bar on the display.
    progressbar -> Start
    
    imagePos = DBLARR(nrfile, 3) * 0.
    FOR i=0L, nrfile-1L DO BEGIN
      count=(i+1.)/nrfile*100.0
      progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
      
      obj = OBJ_NEW('IDLffDICOMex', a[i])
      result=QUERY_DICOM(a[i], info)
      IF (obj->QueryValue('0020,0032')) EQ 2 THEN imagePos[i,*] = DOUBLE(obj->GetValue('0020,0032'))
      IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = FLOAT(obj->GetValue('0028,0030')) ELSE pixelSp = [0., 0.]
      OBJ_DESTROY, obj
      
      origin_loc = imagePos[i,*]
      reslice2d_resolution = pixelSp[0]
      reslice2d=read_dicom(a[i])
      new_reslice2d = volume_reslice(reslice2d, reslice2d_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
      WRITE_JPEG, a[i]+'.jpeg', BYTSCL(new_reslice2d), QUALITY=200
      
;      IF i NE -1 THEN BEGIN
      IF i EQ 17 THEN BEGIN
        iimage, new_reslice2d
        iimage, reslice2d
        ;blend display the result
        s1=new_reslice2d
        s1(where(S1 LT 30)) = 0;s1(where(S1 LT 300)) = 0
        s0=reslice2d
        backgroundImage=DOUBLE(S0)/DOUBLE(max(S0))*255.
        foregroundImage1=s1/2000.*255.;foregroundImage1=s1/2000.*255.
        
        major=5
        max_ct=2000;max_ct=2000
        min_ct=0
        filename1=a[i]+'_auc.jpeg';'unet.jpg'
        closeYes=0
        Image_Blend, REVERSE(backgroundImage,2), REVERSE(foregroundImage1,2), COLORTABLE=33, blendTitle='T2_amp Image', major, max_ct, min_ct, filename1, closeYes, 'test'

        
      ENDIF
    ENDFOR
    
    progressbar -> Destroy
;    i=45
;    origin_loc = imagePos[i,*]
;    reslice2d_resolution = pixelSp[0] 
;    reslice2d=read_dicom(a[i])
;    iimage, reslice2d
;
;    new_reslice2d = volume_reslice(reslice2d, reslice2d_resolution, direction, origin_loc, vol_HU_cube1, vol_HU_cube_resolution1, direction1, origin_loc1)
;
;    iimage, new_reslice2d

  ENDIF

END

;MIP of breast CTA
;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21
;@Guang Jia, gjia@xidian.edu.cn, 2020/01/02, Generating 3 directions MIP without any rotations
PRO mip_prostate_MRA_3directions, dcm_dir, dcm_dir_c

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    ;dcm_dir = 'C:\D_drive\AIMIS_3D\images\ProstateS2001_111754_DYNAMIC\V01_\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    dcm_dir = 'C:\D_drive\AIMIS_3D\images\Bladderpatient\P06\B01_20090825\Dyn01\S3101_145707_DYNAMIC_FFE\V02_\'
    ;dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  image_read_512, dcm_dir, vol_HU_cube_b, vol_HU_cube_resolution, patient_age, direction, origin_loc

  IF N_ELEMENTS(dcm_dir_c) EQ 0 THEN BEGIN
    ;dcm_dir_c = 'C:\D_drive\AIMIS_3D\images\ProstateS2001_111754_DYNAMIC\V10_\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    dcm_dir_c = 'C:\D_drive\AIMIS_3D\images\Bladderpatient\P06\B01_20090825\Dyn01\S3101_145707_DYNAMIC_FFE\V30_\'
    ;dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir_c,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  image_read_512, dcm_dir_c, vol_HU_cube_c, vol_HU_cube_resolution, patient_age, direction, origin_loc

  vol_HU_cube = vol_HU_cube_c - vol_HU_cube_b
  vol_HU_cube[where(vol_HU_cube LT 0, /NULL)] = 0

  size_vol = SIZE(vol_HU_cube, /dimensions)
  mip_slice2d = DBLARR(size_vol[0], size_vol[1])



  ;xy plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR j=0, size_vol[1]-1 DO BEGIN
      mip_slice2d[i, j] = MAX(vol_HU_cube[i, j, *])
    ENDFOR
  ENDFOR
  FILE_MKDIR, dcm_dir+'MIP_Unet\'
  WRITE_JPEG,  dcm_dir+'MIP_Unet\000.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

  ;yz plane
  FOR j=0, size_vol[1]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[j, k] = MAX(vol_HU_cube[*, j, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\001.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

  ;yz plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[i, k] = MAX(vol_HU_cube[i, *, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\002.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

END


;MIP of breast CTA
;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21
;@Guang Jia, gjia@xidian.edu.cn, 2020/01/02, Generating 3 directions MIP without any rotations
PRO mip_breast_CTA_3directions, dcm_dir, dcm_dir_c

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    dcm_dir = 'C:\D_drive\Tandu_BreastCancer\32018_LI_FENG_LI_20181029_CHEST_32018\S0002_091545_C+20181029\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    ;dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  image_read_512, dcm_dir, vol_HU_cube_b, vol_HU_cube_resolution, patient_age, direction, origin_loc
  
  IF N_ELEMENTS(dcm_dir_c) EQ 0 THEN BEGIN
    dcm_dir_c = 'C:\D_drive\Tandu_BreastCancer\32018_LI_FENG_LI_20181029_CHEST_32018\S0006_092036_C+20181029\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    ;dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir_c,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  image_read_512, dcm_dir_c, vol_HU_cube_c, vol_HU_cube_resolution, patient_age, direction, origin_loc
  
  vol_HU_cube = vol_HU_cube_c - vol_HU_cube_b
  vol_HU_cube[where(vol_HU_cube LT 0, /NULL)] = 0
  
  size_vol = SIZE(vol_HU_cube, /dimensions)
  mip_slice2d = DBLARR(size_vol[0], size_vol[1])
  
  

  ;xy plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR j=0, size_vol[1]-1 DO BEGIN
      mip_slice2d[i, j] = MAX(vol_HU_cube[i, j, *])
    ENDFOR
  ENDFOR
  FILE_MKDIR, dcm_dir+'MIP_Unet\'
  WRITE_JPEG,  dcm_dir+'MIP_Unet\000.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

  ;yz plane
  FOR j=0, size_vol[1]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[j, k] = MAX(vol_HU_cube[*, j, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\001.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

  ;yz plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[i, k] = MAX(vol_HU_cube[i, *, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\002.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

END


;MIP of MRA
;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21
;@Guang Jia, gjia@xidian.edu.cn, 2020/01/02, Generating 3 directions MIP without any rotations
PRO mip_mra_3directions, dcm_dir

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    ;dcm_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0201_135520_s3DI_MC_HR\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

  image_read_512, dcm_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
  size_vol = SIZE(vol_HU_cube, /dimensions)
  mip_slice2d = DBLARR(size_vol[0], size_vol[1])
  
  ;xy plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR j=0, size_vol[1]-1 DO BEGIN
      mip_slice2d[i, j] = MAX(vol_HU_cube[i, j, *])
    ENDFOR
  ENDFOR
  FILE_MKDIR, dcm_dir+'MIP_Unet\'
  WRITE_JPEG,  dcm_dir+'MIP_Unet\000.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200
  
  ;yz plane
  FOR j=0, size_vol[1]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[j, k] = MAX(vol_HU_cube[*, j, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\001.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200
  
  ;yz plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[i, k] = MAX(vol_HU_cube[i, *, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\002.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200

END

;MIP of MRA
;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21
PRO mip_mra, dcm_dir

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    ;dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0201_135520_s3DI_MC_HR\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF
  
  image_read_512, dcm_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
  size_vol = SIZE(vol_HU_cube, /dimensions)
  mip_slice2d = DBLARR(size_vol[0], size_vol[1], 18)
  
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='red', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start


  FOR k=0, 17 DO BEGIN
    count=(k+1.)/18*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

    TRANSFORM_VOLUME, vol_HU_cube, temp_vol_HU_cube, Rotation=[10*k, 0, 0]
    FOR i=0, size_vol[0]-1 DO BEGIN
      FOR j=0, size_vol[1]-1 DO BEGIN
        mip_slice2d[i, j] = MAX(temp_vol_HU_cube[i, j, *])
        ENDFOR
      ENDFOR
      WRITE_JPEG,  dcm_dir+'x'+STRING(k, format='(I3.3)')+'.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200
  ENDFOR
  progressbar -> Destroy
  
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='yellow', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start
  FOR k=0, 17 DO BEGIN
    count=(k+1.)/18*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    TRANSFORM_VOLUME, vol_HU_cube, temp_vol_HU_cube, Rotation=[0, 10*k, 0]
    FOR i=0, size_vol[0]-1 DO BEGIN
      FOR j=0, size_vol[1]-1 DO BEGIN
        mip_slice2d[i, j] = MAX(temp_vol_HU_cube[i, j, *])
      ENDFOR
    ENDFOR
    WRITE_JPEG,  dcm_dir+'y'+STRING(k, format='(I3.3)')+'.jpeg', BYTSCL(mip_slice2d[*,*]), QUALITY=200
  ENDFOR
  
  progressbar -> Destroy
  
END


;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21 created.
;@Guang Jia, gjia@xidian.edu.cn, 2020/1/4, modified to be a real program
;@Guang Jia, gjia@xidian.edu.cn, 2020/3/7, modified to be a real program by adding a real reconstruction part
;@Guang Jia, gjia@xidian.edu.cn, 2020/3/15, modified to use CE instead of signal difference
;vessel wall reconstruction
PRO VW_reconstruction
  ;load the raw TOF image cube
  ;dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\3826893_WANG_QI_XI_20181031_HR_MRA\S0201_120213_s3DI_MC_HR\'
  ;dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\2163399_ZHU_KUI_HE_20170824_MCA\S0301_173624_s3DI_MC_HR\'
  dcm_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0201_135520_s3DI_MC_HR\'
  image_read, dcm_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc

  BB_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0401_140721_3D_T1WIVesselWall\'
  BBC_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S1401_144602_3D_T1WITRA+C\'
  ;BB_dir = 'D:\AIMIS_3D\images\MRI_BB\2163399_ZHU_KUI_HE_20170824_MCA\aS0901_180949_3D_T1WIVesselWall\'
  ;BBC_dir = 'D:\AIMIS_3D\images\MRI_BB\2163399_ZHU_KUI_HE_20170824_MCA\bS1301_183344_3D_T1WITRA+C\'
  ;load the raw image cube
  image_read, BB_dir, ori_BB_vol_HU_cube, BBvol_HU_cube_resolution, patient_age, BBdirection, BBorigin_loc
  image_read, BBC_dir, ori_BBC_vol_HU_cube, BBCvol_HU_cube_resolution, patient_age, BBCdirection, BBCorigin_loc
  BB_vol_HU_cube = volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, ori_BB_vol_HU_cube, BBvol_HU_cube_resolution, BBdirection, BBorigin_loc)
  BBC_vol_HU_cube = volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, ori_BBC_vol_HU_cube, BBCvol_HU_cube_resolution, BBCdirection, BBCorigin_loc)
  
  ;difference between BBC and BB
  ;remove negative signal intensities
  diff_vol_HU_cube = BBC_vol_HU_cube - BB_vol_HU_cube
  ;remove background
  diff_vol_HU_cube[where(BB_vol_HU_cube EQ -5, /NULL)] = 0
  diff_vol_HU_cube[where(BBC_vol_HU_cube EQ -5, /NULL)] = 0
  CE_vol_HU_cube = diff_vol_HU_cube/BB_vol_HU_cube
;  iimage, vol_HU_cube[*, *, 127]
;  iimage, BB_vol_HU_cube[*, *, 127]
;  iimage, BBC_vol_HU_cube[*, *, 127]
;  iimage, diff_vol_HU_cube[*, *, 127]
  
  unet_reconstruction_3directions, dcm_dir, vessel_mask_vol, smoothedOutverts, p
  oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [0,127,127], smoothedOutverts, POLYGON=p)
  
  ;dilate the unet slice2d
  size_vol = SIZE(vessel_mask_vol, /dimensions)
  radius = 1
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
  FOR i=0, size_vol[2]-1 DO BEGIN
   temp = DILATE(vessel_mask_vol[*,*,i], strucElem)
   CE_vol_HU_cube[*,*,i] = temp * CE_vol_HU_cube[*,*,i]
  ENDFOR
  
  ;use mesh_decimate to save stl file
  SHADE_VOLUME, CE_vol_HU_cube,0.01, v0, p0
  v02 = MESH_SMOOTH(v0, p0, ITERATIONS=50)
  
  oPoly2 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], v02, POLYGON=p0, ALPHA_CHANNEL=1)
  mymodel1 = OBJ_NEW('IDLgrModel')
  mymodel1->Add, oPoly1
  mymodel1->Add, oPoly2
  xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0], TITLE='Vessel Wall Enhancement'
  
END

;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21 created.
;@Guang Jia, gjia@xidian.edu.cn, 2020/1/4, modified to be a real program
;@Guang Jia, gjia@xidian.edu.cn, 2020/3/7, modified to be a real program by adding a real reconstruction part
;@Guang Jia, gjia@xidian.edu.cn, 2020/3/15, modified to use CE instead of signal difference
;@Guang Jia, gjia@xidian.edu.cn, 2020/3/19, if there is no U-net result, use threshold to do analysis
;vessel wall reconstruction
PRO VW_CE_reconstruction, TOF_dir, BBC_dir, BB_dir
  ;load the raw TOF image cube
  ;dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\3826893_WANG_QI_XI_20181031_HR_MRA\S0201_120213_s3DI_MC_HR\'
  ;dcm_dir = 'D:\AIMIS_3D\images\MRI_BB\2163399_ZHU_KUI_HE_20170824_MCA\S0301_173624_s3DI_MC_HR\'
  IF N_ELEMENTS(TOF_dir) EQ 0 THEN BEGIN
    temp_TOF_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0201_135520_s3DI_MC_HR\'
    TOF_dir = DIALOG_PICKFILE(PATH=temp_TOF_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select TOF dcm folder:')
    filearr=file_search(TOF_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF
  image_read, TOF_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
  
  ;Check whether there is any unet vessl segmentation result
  Unet_filearr=file_search(TOF_dir,'Unet*.png',count=Unet_num)
  IF Unet_num EQ 0 THEN BEGIN
    sizeVHc = SIZE(vol_HU_cube)
    mean_image = IMAGE_THRESHOLD(BYTSCL(vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
    image_statistics, vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
    vessel_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
    ;get rid of the maximum 15 connected regions
    growmaxfive,vol_HU_cube,vessel_value,img_max,vessel_mask_vol
    ;use mesh_decimate to save stl file
    SHADE_VOLUME, vessel_mask_vol,0.5, v, p
    smoothedOutverts = MESH_SMOOTH(v, p, ITERATIONS=50)
  ENDIF ELSE BEGIN
    unet_reconstruction_3directions, TOF_dir, vessel_mask_vol, smoothedOutverts, p
  ENDELSE
  oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [0,127,127], smoothedOutverts, POLYGON=p)
  
  IF N_ELEMENTS(BBC_dir) EQ 0 THEN BEGIN
    ;BBC_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S1401_144602_3D_T1WITRA+C\'
    BBC_dir = DIALOG_PICKFILE(PATH=TOF_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select BB+C dcm folder:')
    filearr=file_search(BBC_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF

;  ;read BB directory
;  IF N_ELEMENTS(BB_dir) EQ 0 THEN BEGIN
;    ;BB_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0401_140721_3D_T1WIVesselWall\'
;    BB_dir = DIALOG_PICKFILE(PATH=TOF_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Black Blood dcm folder:')
;    filearr=file_search(BB_dir,'*.dcm',count=num)
;    IF num EQ 0 THEN RETURN
;  ENDIF
;  ;load the raw image cube
;  image_read, BB_dir, ori_BB_vol_HU_cube, BBvol_HU_cube_resolution, patient_age, BBdirection, BBorigin_loc
;  BB_vol_HU_cube = volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, ori_BB_vol_HU_cube, BBvol_HU_cube_resolution, BBdirection, BBorigin_loc)
  
  ;load the BBC image cube
  image_read, BBC_dir, ori_BBC_vol_HU_cube, BBCvol_HU_cube_resolution, patient_age, BBCdirection, BBCorigin_loc
  BBC_vol_HU_cube = volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, ori_BBC_vol_HU_cube, BBCvol_HU_cube_resolution, BBCdirection, BBCorigin_loc)
  size_vol = SIZE(BBC_vol_HU_cube, /dimensions)
    
  ;check BBC
  sizeVHc = SIZE(BBC_vol_HU_cube)
  mean_image = IMAGE_THRESHOLD(BYTSCL(BBC_vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
  image_statistics, BBC_vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
  CE_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
  
  ;get rid of the maximum 15 connected regions
  growNOTmax15,BBC_vol_HU_cube,CE_value,img_max,BBC_vol_HU_cube_mask
  ;find number of connected regions left
  labelImg = LABEL_REGION(BBC_vol_HU_cube_mask, ALL_NEIGHBORS=allNeighbors, ULONG=ulong)
  ;find the maximum number of connected reions
  labelmax=max(labelImg)
  
  ;define the vessel wall mask of enhanced region
  temp_VW_BBC_vol_HU_cube_mask = BBC_vol_HU_cube *0.
  VW_BBC_vol_HU_cube_mask = BBC_vol_HU_cube *0.
  ;dilate the unet slice2d
  radius = 1
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
  FOR i=0, size_vol[2]-1 DO BEGIN
    temp = DILATE(vessel_mask_vol[*,*,i], strucElem)
    ;define the vessel wall mask of enhanced region
    temp_VW_BBC_vol_HU_cube_mask[*,*,i] = temp * BBC_vol_HU_cube_mask[*,*,i]
  ENDFOR
  
  ;defind the connected regions within vessl wall region
  VW_labelImg = labelImg * temp_VW_BBC_vol_HU_cube_mask
  FOR index = 1, labelmax DO BEGIN
    VW_indexNum = WHERE(VW_labelImg EQ index, VW_count)
    indexNum = WHERE(labelImg EQ index, count)
    IF (VW_count GT 1) AND (VW_count*1.0/count*1.0 GT 0.2) THEN BEGIN
      VW_BBC_vol_HU_cube_mask(indexNum) = 1
      print, 'count = ', count, 'count ratio = ', VW_count*1.0/count*1.0
    ENDIF
  ENDFOR
  
  ;generate suface
  SHADE_VOLUME, VW_BBC_vol_HU_cube_mask,0.5, v0, p0
  v0_smooth = MESH_SMOOTH(v0, p0, ITERATIONS=50)

  oPoly0 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], v0_smooth, POLYGON=p0, ALPHA_CHANNEL=1)
  mymodel1 = OBJ_NEW('IDLgrModel')
  mymodel1->Add, oPoly1
  mymodel1->Add, oPoly0
  xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0], TITLE='Vessel Wall Enhancement'
  
  ;use mesh_decimate to save stl file
;  SHADE_VOLUME, BBC_vol_HU_cube_mask,0.5, v00, p00
  SHADE_VOLUME, BBC_vol_HU_cube,CE_value, v00, p00
  v02 = MESH_SMOOTH(v00, p00, ITERATIONS=50)

  oPoly2 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], v02, POLYGON=p00, ALPHA_CHANNEL=1)
;  oPoly3 = OBJ_NEW('IDLgrPolygon',COLOR = [0,127,127], smoothedOutverts, POLYGON=p)
  mymodel2 = OBJ_NEW('IDLgrModel')
;  mymodel2->Add, oPoly3
  mymodel2->Add, oPoly2
  xobjview_XDU, mymodel2, BACKGROUND = [0, 0, 0], TITLE='Enhancement'

END


;@Guang Jia, gjia@xidian.edu.cn, 2019/12/21, creating this program
;@Guang Jia, gjia@xidian.edu.cn, 2020/01/01, using only 3 MIP projections for reconstruction
;                                            automatically deterniming threshold 
;@Guang Jia, gjia@xidian.edu.cn, 2020/01/02, using only 3 MIP projections for reconstruction without volume rotations
;@Guang Jia, gjia@xidian.edu.cn, 2020/03/07, Adpative reconstruction with different thresholds in a small cubes
;@Guang Jia, gjia@xidian.edu.cn, 2020/03/08, Rewrite the procedure, unet_reconstruction_3directions, and works well
;input: dcm_dir, mip_dir
;output: vessel_mask_vol, smoothedOutverts, p
PRO unet_reconstruction_3directions, dcm_dir, vessel_mask_vol, smoothedOutverts, p

  IF N_ELEMENTS(dcm_dir) EQ 0 THEN BEGIN
    dcm_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0201_135520_s3DI_MC_HR\'; 'E:\Research\MED_visualization\program_xd_xwt_3D\test_data\S0002_162200\'
    ;dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select Original dcm folder:')
    filearr=file_search(dcm_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF
  
  ;load the raw image cube
;  image_read_512, dcm_dir, origin_vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
  image_read, dcm_dir, origin_vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc
;  origin_vol_HU_cube =CONGRID(origin_vol_HU_cube, 512, 512, 512)
  vol_HU_cube = origin_vol_HU_cube
  size_vol = SIZE(vol_HU_cube, /dimensions)
  print, 'size: ', size_vol
  max_vol = ABS(MAX(vol_HU_cube))
  vessel_mask_vol_x = vol_HU_cube*0.
  vessel_mask_vol_y = vol_HU_cube*0.
  vessel_mask_vol_z = vol_HU_cube*0.
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='red', Text='0%', Title='Loading...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start
  
;  unet_slice2d=[size_vol[0], size_vol[1], 18]
;  IF N_ELEMENTS(mip_dir) EQ 0 THEN BEGIN
;    mip_dir='C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\MIP_Unet\'
;  ENDIF
  
  FOR m=0, 2 DO BEGIN
    count=(m+1.)/3.*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    READ_PNG,  dcm_dir+'Unet'+STRING(m, format='(I3.3)')+'.png', temp_slice2d
    temp_size = SIZE(temp_slice2d, /N_DIMENSIONS)
    IF temp_size EQ 3 THEN BEGIN
      temp = REFORM(temp_slice2d[0, *, *])
    ENDIF ELSE BEGIN
      temp = temp_slice2d
    ENDELSE
    unet_slice2d = CONGRID(temp, size_vol[0], size_vol[1])
    
    ;define the mip and shrink
    temp_mip = read_image(dcm_dir+STRING(m, format='(I3.3)')+'.jpeg')
    size_mip = SIZE(temp_mip, /dimensions)
    center_mip = temp_mip[FLOOR(size_mip[0]*2./5.):FLOOR(size_mip[0]*3./5.), FLOOR(size_mip[1]*2./5.):FLOOR(size_mip[1]*3./5.)]
    ;get the mask region of whole brain
    mask_mip = temp_mip GT MIN(center_mip)-5.
    ;fill the holes
    ; Create the shape operator:
    radius = 4
    strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
    mask_mip_dia = DILATE(mask_mip, strucElem)
    
    ;shrink the borders many times
    ;in order to remove bone
    radius = 2
    strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
    mask_mip_shrink = ERODE(mask_mip_dia, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    mask_mip_shrink = ERODE(mask_mip_shrink, strucElem)
    
    ;define the mask
    maskend=unet_slice2d*0
    temp_label=label_region(unet_slice2d*CONGRID(mask_mip_shrink, 256,256) GT 10, /all_neighbors)
    labelmax=max(temp_label)
    h = HISTOGRAM(temp_label, REVERSE_INDICES=r)
    sort_re=REVERSE(SORT(h))
    if labelmax ge 10 then begin
      for j=1,9 do begin
        maskendmid = temp_label EQ sort_re[j]
        image_statistics, unet_slice2d, mean=mean_unet, count=count_unet, mask=maskendmid
        IF (mean_unet GT 30) THEN maskend+=maskendmid*j
      endfor
    endif
    ;only take the largest ones
 ;   pick_label=(temp_label GT 0) AND (temp_label LT 8)
    unet_slice2d_mask = maskend*unet_slice2d GT 10
;    unet_slice2d_mask = unet_slice2d_dilated GT 10
    ;define the flag
    vessel_flag = 0
    ;define the vessel threshold
    vessel_thres = 0
    
  
    ;back-projection along x direction
    IF m EQ 0 THEN BEGIN
      
      ;for display
      temp_mip = DBLARR(size_vol[0], size_vol[1])
      
      ;i+, j+
      FOR i=0, size_vol[0]-1 DO BEGIN
       FOR j=0, size_vol[1]-1 DO BEGIN
        ;this is for display
        temp_mip[i,j] = MAX(vol_HU_cube[i, j, *])
        ;start processing each pixels
        IF unet_slice2d_mask[i,j] LT 1 THEN BEGIN
          ;non-vessel projections, background
          vessel_flag = 0
          vessel_mask_vol_z[i,j,*] = 0
        ENDIF ELSE BEGIN
          ;vessel projections
          ;get the first pixel's mip value as threshold as flag
          IF vessel_flag LT 1 THEN BEGIN
            vessel_flag = 1
            vessel_thres = MAX(vol_HU_cube[i, j, *])
          ENDIF
          ;check each pixel value and remove the background
          FOR k=0, size_vol[2]-1 DO BEGIN
            vessel_mask_vol_z[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
          ENDFOR
        ENDELSE        
       ENDFOR
     ENDFOR
    
     ;i+, j-
     FOR i=0, size_vol[0]-1 DO BEGIN
       FOR j=size_vol[1]-1, 0, -1 DO BEGIN
         ;this is for display
         temp_mip[i,j] = MAX(vol_HU_cube[i, j, *])
         ;start processing each pixels
         IF unet_slice2d_mask[i,j] LT 1 THEN BEGIN
           ;non-vessel projections, background
           vessel_flag = 0
           vessel_mask_vol_z[i,j,*] = 0
         ENDIF ELSE BEGIN
           ;vessel projections
           ;get the first pixel's mip value as threshold as flag
           IF vessel_flag LT 1 THEN BEGIN
             vessel_flag = 1
             vessel_thres = MAX(vol_HU_cube[i, j, *])
           ENDIF
           ;check each pixel value and remove the background
           FOR k=0, size_vol[2]-1 DO BEGIN
             vessel_mask_vol_z[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
           ENDFOR
         ENDELSE
       ENDFOR
     ENDFOR
    
    ;j+, i+
     FOR j=0, size_vol[1]-1 DO BEGIN
       FOR i=0, size_vol[0]-1 DO BEGIN
         ;this is for display
         temp_mip[i,j] = MAX(vol_HU_cube[i, j, *])
         ;start processing each pixels
         IF unet_slice2d_mask[i,j] LT 1 THEN BEGIN
           ;non-vessel projections, background
           vessel_flag = 0
           vessel_mask_vol_z[i,j,*] = 0
         ENDIF ELSE BEGIN
           ;vessel projections
           ;get the first pixel's mip value as threshold as flag
           IF vessel_flag LT 1 THEN BEGIN
             vessel_flag = 1
             vessel_thres = MAX(vol_HU_cube[i, j, *])
           ENDIF
           ;check each pixel value and remove the background
           FOR k=0, size_vol[2]-1 DO BEGIN
             vessel_mask_vol_z[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
           ENDFOR
         ENDELSE
       ENDFOR
     ENDFOR

     ;j+, i-
     FOR j=0, size_vol[1]-1 DO BEGIN
       FOR i=size_vol[0]-1, 0, -1 DO BEGIN
         ;this is for display
         temp_mip[i,j] = MAX(vol_HU_cube[i, j, *])
         ;start processing each pixels
         IF unet_slice2d_mask[i,j] LT 1 THEN BEGIN
           ;non-vessel projections, background
           vessel_flag = 0
           vessel_mask_vol_z[i,j,*] = 0
         ENDIF ELSE BEGIN
           ;vessel projections
           ;get the first pixel's mip value as threshold as flag
           IF vessel_flag LT 1 THEN BEGIN
             vessel_flag = 1
             vessel_thres = MAX(vol_HU_cube[i, j, *])
           ENDIF
           ;check each pixel value and remove the background
           FOR k=0, size_vol[2]-1 DO BEGIN
             vessel_mask_vol_z[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
           ENDFOR
         ENDELSE
       ENDFOR
     ENDFOR
     
     ;remove the pixels as straight lines
     FOR i=0, size_vol[0]-1 DO BEGIN
       FOR j=0, size_vol[1]-1 DO BEGIN
        IF TOTAL(vessel_mask_vol_z[i,j,*]) GT 0.1*size_vol[2] THEN vessel_mask_vol_z[i,j,*]=0
       ENDFOR
     ENDFOR
    ENDIF
    
    ;back-projection along x direction 
    IF m EQ 1 THEN BEGIN
      
      ;for display
      temp_mip = DBLARR(size_vol[1], size_vol[2])

      ;j+, k+
      FOR j=0, size_vol[1]-1 DO BEGIN
        FOR k=0, size_vol[2]-1 DO BEGIN
          ;this is for display
          temp_mip[j,k] = MAX(vol_HU_cube[*, j, k])
          ;start processing each pixels
          IF unet_slice2d_mask[j,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_x[*,j,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[*, j, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR i=0, size_vol[0]-1 DO BEGIN
              vessel_mask_vol_x[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
      
      ;j+, k-
      FOR j=0, size_vol[1]-1 DO BEGIN
        FOR k=size_vol[2]-1, 0, -1 DO BEGIN
          ;this is for display
          temp_mip[j,k] = MAX(vol_HU_cube[*, j, k])
          ;start processing each pixels
          IF unet_slice2d_mask[j,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_x[*,j,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[*, j, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR i=0, size_vol[0]-1 DO BEGIN
              vessel_mask_vol_x[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
     
      ;k+, j+
      FOR k=0, size_vol[2]-1 DO BEGIN
        FOR j=0, size_vol[1]-1 DO BEGIN
          ;this is for display
          temp_mip[j,k] = MAX(vol_HU_cube[*, j, k])
          ;start processing each pixels
          IF unet_slice2d_mask[j,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_x[*,j,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[*, j, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR i=0, size_vol[0]-1 DO BEGIN
              vessel_mask_vol_x[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
      
      ;k+, j-
      FOR k=0, size_vol[2]-1 DO BEGIN
        FOR j=size_vol[1]-1, 0, -1 DO BEGIN
          ;this is for display
          temp_mip[j,k] = MAX(vol_HU_cube[*, j, k])
          ;start processing each pixels
          IF unet_slice2d_mask[j,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_x[*,j,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[*, j, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR i=0, size_vol[0]-1 DO BEGIN
              vessel_mask_vol_x[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
     
     
      ;remove the pixels as straight lines
      FOR j=0, size_vol[1]-1 DO BEGIN
        FOR k=0, size_vol[2]-1 DO BEGIN
          IF TOTAL(vessel_mask_vol_x[*,j,k]) GT 0.1*size_vol[0] THEN vessel_mask_vol_x[*,j,k]=0
        ENDFOR
      ENDFOR

    ENDIF
 
 
    ;back-projection along y directgion
    IF m EQ 2 THEN BEGIN
      
      ;for display
      temp_mip = DBLARR(size_vol[0], size_vol[2])

      ;i+, k+
      FOR i=0, size_vol[0]-1 DO BEGIN
        FOR k=0, size_vol[2]-1 DO BEGIN
          ;this is for display
          temp_mip[i,k] = MAX(vol_HU_cube[i, *, k])
          ;start processing each pixels
          IF unet_slice2d_mask[i,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_y[i,*,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[i, *, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR j=0, size_vol[1]-1 DO BEGIN
              vessel_mask_vol_y[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR

      ;i+, k-
      FOR i=0, size_vol[0]-1 DO BEGIN
        FOR k=size_vol[2]-1, 0, -1 DO BEGIN
          ;this is for display
          temp_mip[i,k] = MAX(vol_HU_cube[i, *, k])
          ;start processing each pixels
          IF unet_slice2d_mask[i,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_y[i,*,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[i, *, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR j=0, size_vol[1]-1 DO BEGIN
              vessel_mask_vol_y[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
      
      ;k+, i+
      FOR k=0, size_vol[2]-1 DO BEGIN
        FOR i=0, size_vol[0]-1 DO BEGIN
          ;this is for display
          temp_mip[i,k] = MAX(vol_HU_cube[i, *, k])
          ;start processing each pixels
          IF unet_slice2d_mask[i,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_y[i,*,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[i, *, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR j=0, size_vol[1]-1 DO BEGIN
              vessel_mask_vol_y[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
      
      ;k+, i+
      FOR k=0, size_vol[2]-1 DO BEGIN
        FOR i=size_vol[0]-1, 0, -1 DO BEGIN
          ;this is for display
          temp_mip[i,k] = MAX(vol_HU_cube[i, *, k])
          ;start processing each pixels
          IF unet_slice2d_mask[i,k] LT 1 THEN BEGIN
            ;non-vessel projections, background
            vessel_flag = 0
            vessel_mask_vol_y[i,*,k] = 0
          ENDIF ELSE BEGIN
            ;vessel projections
            ;get the first pixel's mip value as threshold as flag
            IF vessel_flag LT 1 THEN BEGIN
              vessel_flag = 1
              vessel_thres = MAX(vol_HU_cube[i, *, k])
            ENDIF
            ;check each pixel value and remove the background
            FOR j=0, size_vol[1]-1 DO BEGIN
              vessel_mask_vol_y[i,j,k] += ((vol_HU_cube[i,j,k] LT vessel_thres) ? 0 : 1)
            ENDFOR
          ENDELSE
        ENDFOR
      ENDFOR
      
      
      ;remove the pixels as straight lines
      FOR i=0, size_vol[0]-1 DO BEGIN
        FOR k=0, size_vol[2]-1 DO BEGIN
          IF TOTAL(vessel_mask_vol_y[i,*,k]) GT 0.1*size_vol[1] THEN vessel_mask_vol_y[i,*,k]=0
        ENDFOR
      ENDFOR
      
    ENDIF
    
    ;blend display the result    
;    s1=unet_slice2d
    s1=unet_slice2d_mask*unet_slice2d
    s0=temp_mip
    backgroundImage=DOUBLE(S0)/DOUBLE(max(S0))*255.
    foregroundImage1=s1/40.*255.
    major=5
    max_ct=40
    min_ct=0
    filename1='unet.jpg'
    closeYes=0
    Image_Blend, REVERSE(backgroundImage,2), REVERSE(foregroundImage1,2), COLORTABLE=33, blendTitle='Unet Image', major, max_ct, min_ct, filename1, closeYes, 'test'
    
;    vol_HU_cube = vol_HU_cube + vessel_mask_vol
    
  ENDFOR
  
  progressbar -> Destroy
  
  sizeVHc = SIZE(origin_vol_HU_cube)
  mean_image = IMAGE_THRESHOLD(BYTSCL(origin_vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
  image_statistics, origin_vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
  origin_bone_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.

  SHADE_VOLUME, origin_vol_HU_cube,origin_bone_value, origin_v, origin_p
  origin_smoothedOutverts = MESH_SMOOTH(origin_v, origin_p, ITERATIONS=100)   ;smooth the segmented volume; GJ 2017/12/22
  oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], origin_smoothedOutverts, POLYGON=origin_p)
  mymodel1 = OBJ_NEW('IDLgrModel')
  mymodel1->Add, oPoly1
  xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0], TITLE='Threshold'


;  SHADE_VOLUME, vessel_mask_vol_x+vessel_mask_vol_y+vessel_mask_vol_z, 0.5, v, p
  ;keep only largest 15 parts
  growmax15,vessel_mask_vol_x+vessel_mask_vol_y+vessel_mask_vol_z,0.5,img_max,vessel_mask_vol
  SHADE_VOLUME, vessel_mask_vol,0.5, v, p

  smoothedOutverts = MESH_SMOOTH(v, p, ITERATIONS=100)   ;smooth the segmented volume; GJ 2017/12/22
  oPoly2 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedOutverts, POLYGON=p)
  mymodel2 = OBJ_NEW('IDLgrModel')
  mymodel2->Add, oPoly2
  xobjview_XDU, mymodel2, BACKGROUND = [0, 0, 0], TITLE='U-net'

END

PRO unet_blend_analysis
  mip_dir='C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S0201_135520_s3DI_MC_HR\S0201_135520_s3DI_MC_HR\'
  unet_dir='C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\mipGT_x\'
  FOR k=0, 1 DO BEGIN
    READ_JPEG, mip_dir+'x'+STRING(k, format='(I3.3)')+'.jpeg', mip_slice2d
    READ_PNG,  unet_dir+'Unetx'+STRING(k+1, format='(I3.3)')+'.png', unet_slice2d
    ;iimage, unet_slice2d
    unet_slice2d = CONGRID(unet_slice2d, 512, 512)
    FOR i=0,511 DO BEGIN
      FOR j=0, 511 DO BEGIN
        mip_slice2d[i,j] = (unet_slice2d[i,j] GT 6)?mip_slice2d[i,j]:0
        unet_slice2d[i,j] = (unet_slice2d[i,j] GT 6)?unet_slice2d[i,j]:0
      ENDFOR
    ENDFOR
    
    s0=mip_slice2d
    s1=unet_slice2d
    nonzeroMask = unet_slice2d NE 0
    
    ;image statistics
    IMAGE_STATISTICS, s0, COUNT = pixelNumber, $
      DATA_SUM = pixelTotal, MASK = nonzeroMask, $
      MAXIMUM = pixelMax, MEAN = pixelMean, $
      MINIMUM = pixelMin, STDDEV = pixelDeviation, $
      SUM_OF_SQUARES = pixelSquareSum, $
      VARIANCE = pixelVariance
    PRINT, ''
    PRINT, 'MASKED IMAGE STATISTICS:'
    PRINT, 'Total Number of Pixels = ', pixelNumber
    PRINT, 'Total of Pixel Values = ', pixelTotal
    PRINT, 'Maximum Pixel Value = ', pixelMax
    PRINT, 'Mean of Pixel Values = ', pixelMean
    PRINT, 'Minimum Pixel Value = ', pixelMin
    PRINT, 'Standard Deviation of Pixel Values = ', $
      pixelDeviation
    PRINT, 'Total of Squared Pixel Values = ', $
      pixelSquareSum
    PRINT, 'Variance of Pixel Values = ', pixelVariance
    
    ; Compute the image histogram, using the default bin size of 1.
    pdf = HISTOGRAM(s0, LOCATIONS=xbin, min=1)
    phisto = PLOT(xbin, pdf, /CURRENT, XRANGE=[0,255], $
      TITLE='Histogram', XTITLE='Pixel value', YTITLE='Frequency', $
      MAX_VALUE=5e5, AXIS_STYLE=1, COLOR='red')
    
    
    backgroundImage=DOUBLE(S0)/DOUBLE(max(S0))*255.
    foregroundImage1=s1/255.*255.

    major=5
    max_ct=255
    min_ct=0
    filename1='unet.jpg'
    closeYes=0

    Image_Blend, REVERSE(backgroundImage,2), REVERSE(foregroundImage1,2), COLORTABLE=33, blendTitle='Unet Image', major, max_ct, min_ct, filename1, closeYes, 'test'
  ENDFOR
END

