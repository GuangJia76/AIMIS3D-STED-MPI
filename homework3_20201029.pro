PRO homework2

  filename = 'C:\D_drive\AIMIS_3D\images\WeiShuFang_2020\t2_fse_tra_20200822_093702\00000012.dcm'

  ;filename = dialog_pickfile(FILTER='*.dcm')
  print, filename
  info = file_info(filename)
  fsize = info.size

  obj = OBJ_NEW('IDLffDICOM', filename)
  ; Get the row & column size of the image(s):
  rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
  cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
  bits = *(obj->GetValue('0028'x,'0100'x,/NO_COPY))[0]
  
  datasize = rows * cols * bits /8. ;in bytes
  
  patientname = *(obj->GetValue('0010'x,'0010'x,/NO_COPY))[0]
  patientage = *(obj->GetValue('0010'x,'1010'x,/NO_COPY))[0]
  
  hospital =  *(obj->GetValue('0008'x,'0080'x,/NO_COPY))[0]
  
  headersize = fsize - datasize
  print, 'headersize = ', headersize
  
  OBJ_DESTROY, obj
      
END

;
;   filename = 'C:\D_drive\AIMIS_3D\images\WeiShuFang_2020\t2_fse_tra_20200822_093702\00000012.dcm'
;  pinfo = homework2_function(filename)
;
FUNCTION homework2_function, filename

;  filename = 'C:\D_drive\AIMIS_3D\images\WeiShuFang_2020\t2_fse_tra_20200822_093702\00000012.dcm'

  ;filename = dialog_pickfile(FILTER='*.dcm')
  print, filename
  info = file_info(filename)
  fsize = info.size

  obj = OBJ_NEW('IDLffDICOM', filename)
  ; Get the row & column size of the image(s):
  rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
  cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
  bits = *(obj->GetValue('0028'x,'0100'x,/NO_COPY))[0]

  datasize = rows * cols * bits /8. ;in bytes

  patientname = *(obj->GetValue('0010'x,'0010'x,/NO_COPY))[0]
  patientage = *(obj->GetValue('0010'x,'1010'x,/NO_COPY))[0]

  hospital =  *(obj->GetValue('0008'x,'0080'x,/NO_COPY))[0]

  headersize = fsize - datasize
  print, 'headersize = ', headersize

  OBJ_DESTROY, obj
  
  return, [patientname, patientage, hospital, STRING(fsize), STRING(datasize), STRING(headersize)]
END

PRO homework3

;afn = 'C:\D_drive\AIMIS_3D\images\BinzhouMC\ssz_20200423\Sorted\2246658_YANG_HONG_QUAN_20190417_BYhead\S0005_205008_TOF_3D_MRA\'
afn = 'C:\D_drive\AIMIS_3D\images\CT_HHYY_spine_cases\test_case_04\'
a = FILE_SEARCH(afn, '*.dcm', count=nrfile)

firstimage = read_image(a[0])
size_image = SIZE(firstimage, /dim)

image_vol = DBLARR(size_image[0], size_image[1], nrfile)

FOR i=0, nrfile-1 DO BEGIN
  image_temp = read_image(a[i])
  image_vol[*, *, i] = image_temp
ENDFOR
 
;ivolume, image_vol 

;temp = image_vol[*,*,52]
;iimage, temp

;temp = image_vol[255,*,*]
;iimage, reform(temp)

;temp = image_vol[*,319,*]
;iimage, reform(temp)

;convert image intensity to CT number
filename = a[0]
obj = OBJ_NEW('IDLffDICOM', filename)
intercept = *(obj->GetValue('0028'x,'1052'x,/NO_COPY))[0]
slope = *(obj->GetValue('0028'x,'1053'x,/NO_COPY))[0]
HU_image_vol = image_vol * slope + intercept

;get image resolution
thickness = *(obj->GetValue('0018'x,'0050'x,/NO_COPY))[0]
pixelsp = *(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0]
resolu = STRSPLIT(pixelsp, '\', /EXTRACT)

obj_destroy, obj

;adjust z direction size
z_direc = nrfile * FLOAT(thickness)/FLOAT(resolu[0])
iso_image_vol = CONGRID(image_vol, size_image[0], size_image[1], FLOOR(z_direc))
;ntemp = iso_image_vol[255,*,*]
;iimage, reform(ntemp)

shade_volume, iso_image_vol, 400, vertex, connec
oPoly = OBJ_NEW('IDLgrPolygon', COLOR = [255, 0, 0], vertex, POLYGON=connec)
oModel = OBJ_NEW('IDLgrModel')
oModel->Add, oPoly
xobjview, oModel, BACKGROUND = [0, 0, 0]
  
END

PRO homework4

  ;afn = 'C:\D_drive\AIMIS_3D\images\BinzhouMC\ssz_20200423\Sorted\2246658_YANG_HONG_QUAN_20190417_BYhead\S0005_205008_TOF_3D_MRA\'
  afn = 'C:\D_drive\AIMIS_3D\images\CT_HHYY_spine_cases\test_case_04\'
  a = FILE_SEARCH(afn, '*.dcm', count=nrfile)

  firstimage = read_image(a[0])
  size_image = SIZE(firstimage, /dim)

  image_vol = DBLARR(size_image[0], size_image[1], nrfile)

  FOR i=0, nrfile-1 DO BEGIN
    image_temp = read_image(a[i])
    image_vol[*, *, i] = image_temp
  ENDFOR

  ;ivolume, image_vol

  ;temp = image_vol[*,*,52]
  ;iimage, temp

  ;temp = image_vol[255,*,*]
  ;iimage, reform(temp)

  ;temp = image_vol[*,319,*]
  ;iimage, reform(temp)

  ;convert image intensity to CT number
  filename = a[0]
  obj = OBJ_NEW('IDLffDICOM', filename)
  intercept = *(obj->GetValue('0028'x,'1052'x,/NO_COPY))[0]
  slope = *(obj->GetValue('0028'x,'1053'x,/NO_COPY))[0]
  HU_image_vol = image_vol * slope + intercept

  ;get image resolution
  thickness = *(obj->GetValue('0018'x,'0050'x,/NO_COPY))[0]
  pixelsp = *(obj->GetValue('0028'x,'0030'x,/NO_COPY))[0]
  resolu = STRSPLIT(pixelsp, '\', /EXTRACT)

  obj_destroy, obj

  ;adjust z direction size
  z_direc = nrfile * FLOAT(thickness)/FLOAT(resolu[0])
  iso_image_vol = CONGRID(image_vol, size_image[0], size_image[1], FLOOR(z_direc))
  ;ntemp = iso_image_vol[255,*,*]
  ;iimage, reform(ntemp)

  shade_volume, iso_image_vol, 400, vertex, connec
  oPoly = OBJ_NEW('IDLgrPolygon', COLOR = [255, 0, 0], vertex, POLYGON=connec)
  oModel = OBJ_NEW('IDLgrModel')
  oModel->Add, oPoly
  xobjview, oModel, BACKGROUND = [0, 0, 0]

END