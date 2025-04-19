FUNCTION findcenter,arr
  edge = size(arr,/DIMENSIONS) & edge = edge[0]
  e = edge-1
  xsum = 0UL & ysum = 0UL & zsum = 0UL
  count = 0
  
  for i = 0,e do begin
    for j = 0,e do begin
      for k = 0,e do begin
        if arr[i,j,k] GT 0 then begin
          xsum = xsum+i
          ysum = ysum+j
          zsum = zsum+k
          count++
        endif
      endfor
    endfor
  endfor
  x=xsum/count & y=ysum/count & z=zsum/count
  center=[x,y,z]
  
  RETURN,center
END

FUNCTION fixedcut,arr,center,loc
  x = center[0] & y = center[1] & z = center[2]
  dxs = loc[0] & dxe = loc[1]
  dys = loc[2] & dye = loc[3]
  dzs = loc[4] & dze = loc[5]
  xs = x+dxs & xe = x+dxe
  ys = y+dys & ye = y+dye
  zs = z+dzs & ze = z+dze
  edge = size(arr,/DIMENSIONS) & edge = edge[0]
  e = edge-1
  
  mask = BytArr(edge,edge,edge)
  mask[xs:xe,ys:ye,zs:ze]=1
  
  result = BytArr(edge,edge,edge)
  for i = 0,e do begin
    for j = 0,e do begin
      for k = 0,e do begin
        if (mask[i,j,k] GT 0) and (arr[i,j,k] GT 0) then result[i,j,k]=1
      endfor
    endfor
  endfor
  
  RETURN,result
END

pro fcfd
  arr = BytArr(100,100,100)
  arr[50:60,50:60,50:60] = 1
  arr[20:30,20:30,20:30] = 1

  INTERVAL_VOLUME, arr, 0.9,2, v, p
  p = TETRA_SURFACE(v, p)
  v = MESH_SMOOTH(v, p, ITERATIONS=200)
  oPoly = OBJ_NEW('IDLgrPolygon', COLOR = [255, 255, 255], v, POLYGON=p)
  model = OBJ_NEW('IDLgrModel')
  model->Add, oPoly
  XOBJVIEW, model, BACKGROUND = [0, 0, 0]
  
  ;找体心函数调用
  center = findcenter(arr)
  print,center
  
  ;固定切割函数调用
  loc = [-25,-15,-25,-15,-25,-15]
  result = fixedcut(arr,center,loc)
;  print,result[20,20,20]
;  print,result[50,50,50]
  
  INTERVAL_VOLUME, result, 0.9,2, v, p
  p = TETRA_SURFACE(v, p)
  v = MESH_SMOOTH(v, p)
  oPoly = OBJ_NEW('IDLgrPolygon', COLOR = [255, 255, 255], v, POLYGON=p)
  model = OBJ_NEW('IDLgrModel')
  model->Add, oPoly
  XOBJVIEW, model, BACKGROUND = [0, 0, 0]
END
