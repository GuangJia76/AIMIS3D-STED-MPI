;$Id: //depot/Release/ENVI51_IDL83/idl/idldir/examples/demo/demosrc/d_objworld2.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;;  Copyright (c) 2013-2021, Guang Jia, gjia@xidian.edu.cn. All
;       rights reserved. Unauthorized reproduction is prohibited.

;+
;  FILE:
;       d_objworld2.pro
;
;  CALLING SEQUENCE: d_objworld2
;
;  PURPOSE:
;       Visualization of medical image data.
;
;  MAJOR TOPICS: Visualization
;
;  CATEGORY:
;       
;
;  INTERNAL FUNCTIONS and PROCEDURES:
;       pro d_objworld2Event      -  Event handler
;       pro d_objworld2Cleanup    -  Cleanup
;       pro d_objworld2           -  Main procedure
;
;  EXTERNAL FUNCTIONS, PROCEDURES, and FILES:
;       pro trackball__define   -  Create the trackball object
;       pro IDLexModelManip__define - Define Model Manipulator
;       pro IDLexViewManip__define  - Define View Manipulator
;       pro demo_gettips        - Read the tip file and create widgets
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
;       1/97,   ACY   - adapted from objworld2, written by R.F.
;       7/98,   PCS   - changed GUI to require only the left
;                       mouse-button.  Changed code to use
;                       IDLexModelManip and IDLexViewManip.
;       2018/5/15, Guang Jia  -click a point at 3D model with a specific location, then 2D cross-sectional image will be displayed at the selected location
;       2018/5/17, Guang Jia  -xd_ROI was added to display cross-sectional images
;       2018/5/18a, Guang Jia  -xd_ROI can display mask ROIs
;-
;----------------------------------------------------------------------------
function d_objworld2ReadPalette, n

  filename = filepath('colors1.tbl', subdir=['resource', 'colors'])
  openr,lun,filename, /block, /get_lun
  ntables = 0b
  readu, lun, ntables
  tnum = n
  if (tnum LT 0) then tnum = 0
  if (tnum GE ntables) then tnum = ntables-1
  arr = bytarr(256)
  ctab = bytarr(3,256)
  point_lun, lun, tnum*768L + 1L
  readu, lun, arr
  ctab[0,*] = arr
  readu, lun, arr
  ctab[1,*] = arr
  readu, lun, arr
  ctab[2,*] = arr
  close,lun
  free_lun,lun

  return,ctab

end

;----------------------------------------------------------------------------
pro d_objworld2SceneAntiAlias, w, v, n
  widget_control, /hourglass
  case n of
    2 : begin
      jitter = [ $
        [ 0.246490,  0.249999], $
        [-0.246490, -0.249999] $
        ]
      njitter = 2
    end
    3 : begin
      jitter = [ $
        [-0.373411, -0.250550], $
        [ 0.256263,  0.368119], $
        [ 0.117148, -0.117570] $
        ]
      njitter = 3
    end
    8 : begin
      jitter = [ $
        [-0.334818,  0.435331], $
        [ 0.286438, -0.393495], $
        [ 0.459462,  0.141540], $
        [-0.414498, -0.192829], $
        [-0.183790,  0.082102], $
        [-0.079263, -0.317383], $
        [ 0.102254,  0.299133], $
        [ 0.164216, -0.054399] $
        ]
      njitter = 8
    end
    else : begin
      jitter = [ $
        [-0.208147,  0.353730], $
        [ 0.203849, -0.353780], $
        [-0.292626, -0.149945], $
        [ 0.296924,  0.149994] $
        ]
      njitter = 4
    end
  endcase

  w->GetProperty, dimension=d
  acc = fltarr(3, d[0], d[1])

  if obj_isa(v, 'IDLgrView') then begin
    nViews = 1
    oViews = objarr(1)
    oViews[0] = v
  end $
  else begin
    nViews = v->count()
    oViews = v->get(/all)
  end

  rects = fltarr(4, nViews)
  for j=0,nViews-1 do begin
    oViews[j]->GetProperty, viewplane_rect=viewplane_rect
    rects[*,j] = viewplane_rect
  end

  for i=0,njitter-1 do begin
    for j=0,nViews-1 do begin
      sc = rects[2,j] / float(d[0])
      oViews[j]->setproperty, view=[ $
        rects[0,j] + jitter[0,i] * sc, $
        rects[1,j] + jitter[1,i] * sc, $
        rects[2,j], $
        rects[3,j] $
        ]
    end

    demo_draw, w, v
    img = w->read()
    img->GetProperty ,data=data
    acc = acc + float(data)
    obj_destroy, img

    for j=0,nViews-1 do begin
      oViews[j]->setproperty, viewplane_rect=rects[*,j]
    end
  end

  acc = acc / float(njitter)

  o = obj_new('IDLgrImage',acc)

  v2 = obj_new('IDLgrView', view=[0,0,d[0],d[1]], proj=1)
  m = obj_new('IDLgrModel')
  m->add, o
  v2->add, m

  demo_draw, w, v2

  obj_destroy, v2

end

;----------------------------------------------------------------------------
function d_objworld2MakeView, xdim, ydim, uval, icon_dir, baseplate
  ;
  ;Compute viewplane rect based on aspect ratio.
  ;
  aspect = xdim / float(ydim)
  myview = [-1, -1, 2, 2] * sqrt(2)
  if (aspect > 1) then begin
    myview[0] = myview[0] - ((aspect-1.0)*myview[2])/2.0
    myview[2] = myview[2] * aspect
  end $
  else begin
    myview[1] = myview[1] - (((1.0/aspect)-1.0)*myview[3])/2.0
    myview[3] = myview[3] / aspect
  end
  ;
  ;Create view.
  ;
  v = obj_new('IDLgrView', $
    projection=2, $
    eye=3, $
    zclip=[1.5,-1.5], $
    dim=[xdim,ydim],$
    viewplane_rect=myview, $
    color=[0,0,0], $
    uvalue=uval $
    )
  ;

  ;
  ;Create model.
  ;
  gg = obj_new('IDLgrModel')
  g = obj_new('IDLgrModel')
  gg->add,g
  ;
  ;Create the base plate.
  ;
  b_verts = fltarr(3,5,5)
  b_conn = lonarr(5,16)
  vert_cols = bytarr(3,25)

  j = 0
  for i=0,15 do begin
    b_conn[0,i] = 4
    b_conn[1,i] = j
    b_conn[2,i] = j+1
    b_conn[3,i] = j+6
    b_conn[4,i] = j+5
    j = j + 1
    if (j MOD 5) EQ 4 then $
      j = j + 1
  end

  k = 0
  for i=0,4 do begin
    for j=0,4 do begin
      b_verts[0,i,j] = i
      b_verts[1,i,j] = j
      b_verts[2,i,j] = 0
      if (k EQ 1) then begin
        vert_cols[*, i+j*5] = [40,40,40]
      end $
      else begin
        vert_cols[*, i+j*5] = [255,255,255]-40
      end
      k = 1 - k
    end
  end

  b_verts[0,*,*] = (b_verts[0,*,*]-2)/2.0
  b_verts[1,*,*] = (b_verts[1,*,*]-2)/2.0
  ;baseplate = obj_new('IDLgrPolygon', $
  ;    b_verts, $
  ;    poly=b_conn, $
  ;    shading=0, $
  ;    vert_colors=vert_cols $
  ;    )
;  Dir_imgXD= FILE_WHICH('XidianUniv.jpg',/INCLUDE_CURRENT_DIR)
;  Dir_imgXWT=FILE_WHICH('xwt_xd.png',/INCLUDE_CURRENT_DIR)
  READ_JPEG, icon_dir+'AIMIS3D_icon.jpg', imageXD
;  imageXWT = READ_PNG(icon_dir+'xwt_bdiv.png')
  imageXWT = READ_PNG(icon_dir+'AIMIS3D5.png')
  dims_XD = SIZE(imageXD, /DIMENSIONS)
  dims_XWT = SIZE(imageXWT, /DIMENSIONS)
  ma=max(dims_XWT[1:2])
  mi=0.;min(dims_XD[1:2])
  bias=(ma-mi)/2.
  rl=mi+bias
  rr=ma+bias
  cc=[(-rl)/(rr-rl),1/(rr-rl)]
  baseplate=obj_new('IDLgrImage', imageXWT, ALPHA_CHANNEL=0.3, XCOORD_CONV=cc*1., YCOORD_CONV=cc*0.3, NAME = 'CoverImage');

  g->add, baseplate
  ;
;Define the object tree add point.
;
g->add, obj_new('IDLgrModel')
;
;Add some lights.
;
gg->add, obj_new('IDLgrLight', $
    loc=[2,2,5], $
    type=2, $ ; Directional (parallel rays).
    color=[255,255,255], $
    intensity=.5 $
    )
gg->add, obj_new('IDLgrLight', $
    type=0, $ ; Ambient.
    intensity=.5, $
    color=[255,255,255] $
    )
;
;Place the model in the view.
;
;gg->Rotate, [-1,0,0], 90
;gg->Rotate, [0,0,1], -145
v->add, gg
  ;
  return, v

end

;----------------------------------------------------------------------------
pro d_objworld2GetViewObjs, view, oWorldRotModel, oBasePlatePolygon, model_top
  ;
  ;Luckily, this portion of the hierarchy is fixed...
  ;
  gg = view->get()
  oWorldRotModel = gg->get(pos=0)
  oBasePlatePolygon = oWorldRotModel->get(pos=0)
  model_top = oWorldRotModel->get(pos=1)

end

;----------------------------------------------------------------------------
pro d_objworld2Cone,verts,conn,n

  verts = fltarr(3,n+1)
  verts[0,0] = 0.0
  verts[1,0] = 0.0
  verts[2,0] = 0.1
  t = 0.0
  tinc = (2.*!PI)/float(n)
  for i=1,n do begin
    verts[0,i] = 0.1*cos(t)
    verts[1,i] = 0.1*sin(t)
    verts[2,i] = -0.1
    t = t + tinc
  end
  conn = fltarr(4*n+(n+1))
  i = 0
  conn[0] = n
  for i=1,n do conn[i] = (n-i+1)
  j = n+1
  for i=1,n do begin
    conn[j] = 3
    conn[j+1] = i
    conn[j+2] = 0
    conn[j+3] = i + 1
    if (i EQ n) then conn[j+3] = 1
    j = j + 4
  end
end

;----------------------------------------------------------------------------
function d_objworld2MakeGroupObj,type,thefont,numObj

  oModel  = obj_new('IDLgrModel')
  case type of
     0: begin
      inputFile = './lib/stl_template/terracotta/terracotta1.stl'
      readstl,inputfile,smoothedOutverts,Outconn1
      print, 'finish .stl file loading'
      str = "terracotta warrior"
     end
     
     1: begin
       inputFile = './lib/stl_template/animal/hollow_pig.stl'
       readstl,inputfile,smoothedOutverts,Outconn1
       print, 'finish .stl file loading'
       str = "animal"
      
     end
     
     2: begin
       inputFile = './lib/stl_template/plant/open_shine_rose-Cube.stl'
       readstl,inputfile,smoothedOutverts,Outconn1
       print, 'finish .stl file loading'
       str = "plant"

     end
     
     3: begin
       inputFile = './lib/stl_template/mermaid/mermaid_t.stl'
       readstl,inputfile,smoothedOutverts,Outconn1
       print, 'finish .stl file loading'
       str = "plant"

     end

     4: begin
       inputFile = './lib/stl_template/humanbody/Explore_Organs.stl'
       readstl,inputfile,smoothedOutverts,Outconn1
       print, 'finish .stl file loading'
       str = "human body"
     end
     
     5: begin
       inputFile = './lib/stl_template/computerparts/Lucas_Computer.stl'
       readstl,inputfile,smoothedOutverts,Outconn1
       print, 'finish .stl file loading'
       str = "hardware"
     end
     
     6: begin
       inputFile = './lib/stl_template/key_ring/key_ring.stl'
       readstl,inputfile,smoothedOutverts,Outconn1
       print, 'finish .stl file loading'
       str = "hardware"
     end
     
     ;read stl rdstl from main program
     7 : begin
       inputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
         TITLE='Select stl File', FILTER='*.stl', /MUST_EXIST)
       ;      inputfile = 'E:\Research\MED_visualization\6_spineCases\Li  guizhen\LI GUIZHEN.stl'
       IF STRLEN(inputFile) GT 0 THEN BEGIN
         readstl,inputfile,smoothedOutverts,Outconn1
         N_ele = N_ELEMENTS(smoothedOutverts[0,*])
         smoothedOutverts[0,N_ele-30:N_ele-1] = smoothedOutverts[0,N_ele-31]
         smoothedOutverts[1,N_ele-30:N_ele-1] = smoothedOutverts[1,N_ele-31]
         smoothedOutverts[2,N_ele-30:N_ele-1] = smoothedOutverts[2,N_ele-31]
         print, 'finish .stl file loading'
         str = "stl file"
       ENDIF ELSE BEGIN
         return, 0
       ENDELSE
     end     
   endcase

     
        IF N_ELEMENTS(smoothedOutverts) GT 10 THEN BEGIN
          v2=smoothedOutverts
          ma=max(v2)
          mi=min(v2)
          bias=(ma-mi)/2
          rl=mi+bias
          rr=ma+bias
          cc=[(-rl)/(rr-rl),1/(rr-rl)]
          xma = max((v2)[0,*])
          xmi = min((v2)[0,*])
          xmid = 0.5*(xma+xmi)*cc[1]
          yma = max((v2)[1,*])
          ymi = min((v2)[1,*])
          ymid = 0.5*(yma+ymi)*cc[1]
          zma = max((v2)[2,*])
          zmi = min((v2)[2,*])
          zmid = 0.5*(zma+zmi)*cc[1]
          ;
          color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
          numObj = numObj+1
          s = obj_new('IDLgrPolygon',smoothedOutverts,poly=Outconn1,color=TRANSPOSE(color_rgb[numObj, *]),$
            shading=0,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
          ;this is used to label number of objects for color
          

        ENDIF
      
  IF OBJ_VALID(s) THEN BEGIN
    oModel->Add, s
    oModel->SetProperty, name=str
;    oModel->GetProperty, transform = old_transform
;    print, 'old transform = ', old_transform
;    oModel->Translate, 0.3,0.4,0
;    oModel->GetProperty, transform = new_transform
;    print, 'new transform = ', new_transform
    ; oModel->SetProperty, uvalue=str
    oModel->Rotate, [-1,0,0], 90
    oModel->Rotate, [0,0,1], 160
  ENDIF

  return, oModel
end


;----------------------------------------------------------------------------
function d_objworld2MakeObj,type,thefont

  oModel  = obj_new('IDLgrModel')
  case type of
      ;read stl rdstl from main program
    0 : begin
      inputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
        TITLE='Select stl File', FILTER='*.stl', /MUST_EXIST)
      ;      inputfile = 'E:\Research\MED_visualization\6_spineCases\Li  guizhen\LI GUIZHEN.stl'
      IF STRLEN(inputFile) GT 0 THEN BEGIN
        readstl,inputfile,smoothedOutverts,Outconn1
        print, 'finish .stl file loading'
        IF N_ELEMENTS(smoothedOutverts) GT 10 THEN BEGIN
                  v2=smoothedOutverts
                  ma=max(v2)
                  mi=min(v2)
                  bias=(ma-mi)/2
                  rl=mi+bias
                  rr=ma+bias
                  cc=[(-rl)/(rr-rl),1/(rr-rl)]
                  xma = max((v2)[0,*])
                  xmi = min((v2)[0,*])
                  xmid = 0.5*(xma+xmi)*cc[1]
                  yma = max((v2)[1,*])
                  ymi = min((v2)[1,*])
                  ymid = 0.5*(yma+ymi)*cc[1]
                  zma = max((v2)[2,*])
                  zmi = min((v2)[2,*])
                  zmid = 0.5*(zma+zmi)*cc[1]
          ;
          ;        xd_sState.oPoly1 =OBJ_NEW('IDLgrPolygon',smoothedOutverts,POLYGON=Outconn1);,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
          color = [[254,206,180],[216,50,46],[230,197,183],[215,220,221],[255,255,255]]
          result = round(randomu(s, 1) * 5)
          s = obj_new('IDLgrPolygon',smoothedOutverts,poly=Outconn1,color=color[*, result],$
            shading=0,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
          str = "stl"

        ENDIF
      ENDIF
    end
    1 : begin
      s = obj_new('orb',color=[255,0,0],radius=0.1,shading=1,$
        select_target=1)
      str = "Sphere"
    end
    2 : begin
      verts = [[-0.1,-0.1,-0.1],[0.1,-0.1,-0.1],[0.1,0.1,-0.1], $
        [-0.1,0.1,-0.1], $
        [-0.1,-0.1,0.1],[0.1,-0.1,0.1],[0.1,0.1,0.1],$
        [-0.1,0.1,0.1]]
      conn = [[4,3,2,1,0],[4,4,5,6,7],[4,0,1,5,4],[4,1,2,6,5], $
        [4,2,3,7,6],[4,3,0,4,7]]
      s = obj_new('IDLgrPolygon',verts,poly=conn,color=[0,255,0],$
        shading=0)
      str = "Cube"
    end
    3 : begin
      d_objworld2Cone,verts,conn,3
      s = obj_new('IDLgrPolygon',verts,poly=conn,$
        color=[0,255,255],shading=0)
      str = "Tetrahedron"
    end
    4 : begin
      d_objworld2Cone,verts,conn,20
      s = obj_new('IDLgrPolygon',verts,poly=conn,$
        color=[255,128,255],shading=1)
      str = "Cone"
    end
    5 : begin
      d_objworld2Cone,verts,conn,4
      l = obj_new('IDLgrPolygon',verts*0.5,poly=conn,$
        color=[100,255,100],shading=0)
      oModel->add,l
      l = obj_new('IDLgrPolyline',[[0,0,0],[0,0,-0.1]],$
        color=[100,255,100])
      oModel->add,l
      s = obj_new('IDLgrLight',loc=[0,0,0],dir=[0,0,-1],cone=40,$
        focus=0,type = 3,color=[100,255,100])
      str = "Green Light"
    end
    6 : begin
      ;  Surface data is read from elevation data file.
      e_height = BYTARR(64,64, /NOZERO)
      OPENR, lun, /GET_LUN, demo_filepath('elevbin.dat', $
        SUBDIR=['examples','data'])
      READU, lun, e_height
      FREE_LUN, lun
      zdata = e_height / (1.7 * max(e_height)) + .001
      xdata = (findgen(64)-32.0)/64.0
      ydata = (findgen(64)-32.0)/64.0
      s = obj_new('IDLgrSurface',zdata,shading=1,style=2,$
        datax=xdata,datay=ydata,color=[150,50,150])
      str = "Surface"
    end
    7 : begin
      ctab = d_objworld2ReadPalette(26)
      restore, demo_filepath('marbells.dat', $
        subdir=['examples','data'])
      image = bytscl(elev, min=2658, max=4241)
      image = image[8:*, *] ; Trim unsightly junk from left side.
      sz = size(image)
      img = bytarr(3,sz[1],sz[2])
      img[0,*,*] = ctab[0,image]
      img[1,*,*] = ctab[1,image]
      img[2,*,*] = ctab[2,image]
      oTextureImage = obj_new('IDLgrImage', $
        img, $
        loc=[0.0,0.0],$
        dim=[0.01,0.01], $
        hide=1 $
        )
      oModel->getProperty,uvalue=uval
      uval = n_elements(uval) GT 0 ? [uval,oTextureImage] : [oTextureImage]
      oModel->setProperty,uvalue=uval
      ;        oModel->add, oTextureImage
      xp=0.5
      yp=0.5*(72./92.)
      zp=0.1
      s=obj_new('IDLgrPolygon',$
        [[-xp,-yp,zp],[xp,-yp,zp],[xp,yp,zp],[-xp,yp,zp]],$
        texture_coord=[[0,0],[1,0],[1,1],[0,1]],$
        texture_map=oTextureImage, $
        color=[255,255,255] $
        )
      str = "Image"
    end
    8 : begin
      d_objworld2Cone, verts, conn, 4
      oModel->add, obj_new('IDLgrPolygon', $
        verts*0.5, $
        poly=conn, $
        color=[255,255,255], $
        shading=0 $
        )
      oModel->add, obj_new('IDLgrPolyline', $
        [[0,0,0], [0,0,-0.1]],$
        color=[255,255,255] $
        )
      s = obj_new('IDLgrLight', $
        loc=[0,0,0], $
        dir=[0,0,-1], $
        cone=20,$
        focus=0, $
        type=3, $
        color=[255,255,255] $
        )
      str = "White Light"
    end
    9 : begin
      s=obj_new('IDLgrText', $
        "XiDian University", $
        location=[0,0,0.001], $
        align=0.5, $
        color=[255,0,255], $
        font=thefont[0] $
        )
      str = "Text"
    end
    10 : begin
      ; Data for plot generated from Magnitude example
      ; in Chapter 13, "Signal Processing", of _Using IDL_.

      N = 1024       ; number of time samples in data set
      delt = 0.02    ; sampling interval in seconds

      U = -0.3 $     ;  DC component
        + 1.0 * Sin(2*!Pi* 2.8 *delt*FIndGen(N)) $ ;2.8 c/s comp
        + 1.0 * Sin(2*!Pi* 6.25*delt*FIndGen(N)) $ ;6.25 c/s comp
        + 1.0 * Sin(2*!Pi*11.0 *delt*FIndGen(N)); 11.0 c/s comp

      V = fft(U) ; compute spectrum v

      ; signal_x is [0.0, 1.0/(N*delt), ... , 1.0/(2.0*delt)]
      signal_x = FINDGEN(N/2+1) / (N*delt)

      mag = ABS(V[0:N/2]); magnitude of first half of v
      signal_y = 20*ALOG10(mag)

      ; phase not used here, included for completeness
      phi = ATAN(V[0:N/2]) ; phase of first half of v

      xc=[-0.5,1.0/25.0]
      yc=[0.5,1.0/80.0]
      s=obj_new('IDLgrPolygon', $
        [[-7,-90,-0.002],[30,-90,-0.002],$
        [30,10,-0.002],[-7,10,-0.002]],$
        color=[0,0,0], $
        xcoord_conv=xc, $
        ycoord_conv=yc $
        )
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        0, $
        range=[0.0,25.0],$
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        location=[0,-80.0], $
        color=[128,60,39], $
        ticklen=5, $
        /exact $
        )
      s->GetProperty,ticktext=tt
      tt->setproperty,font=thefont[3]
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        0, $
        range=[0.0,25.0], $
        /notext,$
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        location=[0.0,0.0], $
        color=[128,60,39], $
        ticklen=-5, $
        /exact $
        )
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        1, $
        range=[-80.0,0.0],$
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        color=[128,60,39], $
        ticklen=1.0, $
        /exact $
        )
      s->GetProperty,ticktext=tt
      tt->setproperty,font=thefont[3]
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        1, $
        range=[-80.0,0.0], $
        /notext,$
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        loc=[25.0,0.0], $
        color=[128,60,39], $
        ticklen=-1.0, $
        /exact $
        )
      oModel->add,s
      s=obj_new('idlgrplot', $
        signal_x, $
        signal_y, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        color=[0,255,255] $
        )
      str = "Plot"
    end
    11 : begin
      x=indgen(200)
      yexp = exp(-x*0.015)
      ysexp = exp(-x*0.015)*sin(x*0.1)
      dataz=fltarr(200,5)
      dataz[*,0] = yexp
      dataz[*,1] = yexp
      dataz[*,2] = REPLICATE(1.1,200)
      dataz[*,3] = ysexp-0.01
      dataz[*,4] = ysexp-0.01
      datay = fltarr(200,5)
      datay[*,0] = 0.0
      datay[*,1] = 1.0
      datay[*,2] = 0.0
      datay[*,3] = 0.0
      datay[*,4] = 1.0
      cbins = bytarr(3,60)
      for i=0,59 do begin
        color_convert, float(i)*4., 1., 1., r,g,b, /HSV_RGB
        cbins[*,59-i] = [r,g,b]
      end
      colors = bytarr(3,200*5)
      colors[0,0:599] = REPLICATE(80,3*200)
      colors[1,0:599] = REPLICATE(80,3*200)
      colors[2,0:599] = REPLICATE(200,3*200)
      colors[*,600:799] = cbins[*,(ysexp+1.0)*30.0]
      colors[*,800:999] = cbins[*,(ysexp+1.0)*30.0]
      xc = [-0.5,1.0/200.0]*0.8
      yc = [-0.5,1.0/1.0]*0.1
      zc = [-0.5,1.0/1.0]*0.4
      s=obj_new('IDLgrAxis', $
        0, $
        range=[0,200],$
        color=[255,255,255], $
        ticklen=0.2, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        2, $
        range=[-1.,1.],$
        color=[255,255,255], $
        ticklen=4, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      oModel->add,s
      s=obj_new('IDLgrSurface', $
        dataz, $
        style=2, $
        vert_colors=colors,$
        datay=datay, $
        max_value=1.05, $
        shading=1, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      oModel->add,s
      s=obj_new('IDLgrSurface', $
        dataz, $
        style=3, $
        color=[0,0,0],$
        datay=datay, $
        max_value=1.05, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      str = 'Ribbon Plot'
    end
    12 : begin
      dataz = dist(8)
      dataz[1,*] = -1
      dataz[3,*] = -1
      dataz[5,*] = -1
      dataz[*,1] = -1
      dataz[*,3] = -1
      dataz[*,5] = -1
      dataz = dataz + 1
      cbins=[ $
        [255,0,0],$
        [255,85,0],$
        [255,170,0],$
        [255,255,0],$
        [170,255,0],$
        [85,255,0],$
        [0,255,0] $
        ]
      colors = bytarr(3, 8*8)
      minz = min(dataz)
      maxz = max(dataz)
      zi = round((dataz - minz)/(maxz-minz) * 6.0)
      colors[*,*] = cbins[*,zi]
      xc = [-0.5,1.0/8.0]*0.4
      yc = [-0.5,1.0/8.0]*0.4
      zc = [0,1.0/8.0]*0.4
      s=obj_new('IDLgrAxis', $
        0, $
        range=[0,8], $
        major=5, $
        color=[255,255,255], $
        ticklen=0.2, $
        /exact, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        1, $
        range=[0,8], $
        major=5,$
        color=[255,255,255], $
        ticklen=0.2, $
        /exact, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      oModel->add,s
      s=obj_new('IDLgrAxis', $
        2, $
        range=[0,8], $
        major=5,$
        color=[255,255,255], $
        ticklen=0.2, $
        /exact, $
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      oModel->add,s
      s=obj_new('IDLgrSurface', $
        dataz, $
        STYLE=6, $
        VERT_COLORS=colors,$
        xcoord_conv=xc, $
        ycoord_conv=yc, $
        zcoord_conv=zc $
        )
      str = 'Bar Plot'
    end
    13 : begin
      ;  Surface data is read from elevation data file.
      e_height = bytarr(64,64, /nozero)
      openr, lun, /get_lun, demo_filepath('elevbin.dat', $
        subdir=['examples','data'])
      readu, lun, e_height
      free_lun, lun

      ;Get image data for texture map.
      read_jpeg, $
        demo_filepath( $
        'elev_t.jpg', $
        subdir=['examples','data'] $
        ), $
        e_image, $
        true=3

      zm=max(e_height)
      zdata = float(e_height)/(1.7*float(zm))
      xdata = (findgen(64)-32.0)/64.0
      ydata = (findgen(64)-32.0)/64.0
      oTextureImage = obj_new('IDLgrImage', $
        e_image, $ ;>255, $
        hide=1, $
        ;            dim=[0.01,0.01], $
        interleave=2 $
        )
      oModel->getProperty,uvalue=uval
      uval = n_elements(uval) GT 0 ? [uval,oTextureImage] : [oTextureImage]
      oModel->setProperty,uvalue=uval
      ;        oModel->add, oTextureImage
      s = obj_new('IDLgrSurface', $
        zdata, $
        shading=1, $
        style=2,$
        datax=xdata, $
        datay=ydata, $
        color=[255,255,255], $
        TEXTURE_MAP=oTextureImage $
        )
      str = "Textured Surface"
    end
  endcase
  IF OBJ_VALID(s) THEN BEGIN
    oModel->Add, s
    oModel->SetProperty, name=str
    ; oModel->SetProperty, uvalue=str
  ENDIF
  
  return, oModel
end
;----------------------------------------------------------------------------
pro d_objworld2NewMode, state, mode
  widget_control, /hourglass
  state.oModelMan->SetProperty, mode=mode
;  widget_control, state.wModelModeRadio, set_value=mode
end
;-----------------------------------------------
;;GJ, 2019/2/7, adding a parameter isaPolygon to return whether the object is a polygon or not-----------------------------
pro d_objworld2Add, state, oModel, isaPolygon, as_child=as_child
  ;
  ;Procedure d_objworldAdd: Add oModel to objworld's graphics tree.
  ;
  if keyword_set(as_child) then $
    state.selected->add, oModel $
  else $
    state.oCurrentTopModel->add, oModel

  state.oCurrentView->GetProperty, uvalue=view_uval
  *(state.model_lists[view_uval.num]) = $
    [oModel, *(state.model_lists[view_uval.num])]
  state.model_cycle_pos = 0
  ;
  ;Make the new object be the current selection...
  ;
  state.selected = oModel
  g = oModel->get(pos=0)
  IF OBJ_VALID(g) THEN BEGIN
    if (obj_isa(g,'IDLgrText')) then begin
    rect = state.win->gettextdimensions(g)
    endif
    
    if (obj_isa(g, 'IDLgrPolygon')) then begin
      isaPolygon = 1
    endif else begin
      isaPolygon = 0
    endelse
  ENDIF
  
  state.oModelMan->SetTarget, state.selected
  state.selected->GetProperty, name=s
  ; state.selected->GetProperty, uvalue=s
  str = "Current selection:" + s
  widget_control, state.text, set_value=str
  widget_control, state.wModelDeleteButton, sensitive=1
  ;widget_control, state.wAddChildButton, sensitive=1
  widget_control, state.wUnselectButton, sensitive=1
  widget_control, state.wModelModeRadio, sensitive=1
  widget_control, state.wSaveButton, sensitive=([1,0])[lmgr(/demo)]
  widget_control, state.wModelSelectButton, $
    sensitive= $
    n_elements(*(state.model_lists[view_uval.num])) gt 2

  demo_draw, state.win, state.scene, debug=state.debug
end

;----------------------------------------------------------------------------
Function d_objworld2ToggleState, wid

  widget_control, wid, get_value=name

  s = strpos(name,'(off)')
  if (s NE -1) then begin
    strput,name,'(on )',s
    ret = 1
  end $
  else begin
    s = strpos(name,'(on )')
    strput,name,'(off)',s
    ret = 0
  end

  widget_control, wid, set_value=name
  return,ret
end
;----------------------------------------------------------------------------
pro d_objworld2Cleanup, wTopBase

  widget_control, wtopbase, get_uvalue=state, /no_copy
  ;
  ;Remove manipulators.
  ;
  state.oModelMan->SetTarget,obj_new()
  state.oViewMan->SetTarget,obj_new(),state.win
  ;
  ;Clean up heap variables.
  ;
  for i=0,n_tags(state)-1 do begin
    case size(state.(i), /TNAME) of
      'POINTER': $
        ptr_free, state.(i)
      'OBJREF': $
        obj_destroy, state.(i)
      else:
    endcase
  end
  ;
  ;Restore the color table.
  ;
  tvlct, state.colortable

  if widget_info(state.groupbase, /valid_id) then $
    widget_control, state.groupbase, /map

end


;---------------------------------------------------------------------------
;'gesture control' code start from here
;author:  THL,    2020/3/31
;update current model's state which showed in current window according to the 'key' and 'action'
;action:  0-translate,1-rotate,2-scale
pro updateModelsByGestures,ev,key,action
  if key EQ -1 then return
  widget_control, ev.top, get_uvalue=state, /no_copy
  case key of
    0:begin
      state.oCurrentTopModel->scale,0.9,0.9,0.9   ;zoom out
    end

    1:begin
      state.oCurrentTopModel->scale,1.1,1.1,1.1    ;zoom in
    end

    2:begin
      ;counterclockwise rotation relative to X axis
      state.oCurrentTopModel->rotate,[1,0,0],-9
      state.oCurrentTopModel->rotate,[0,1,0],3    ; correct the deviation of Y axis
      ;counterclockwise rotation relative to X axis
    end

    3:begin
      ;clockwise rotation relative to X axis
      state.oCurrentTopModel->rotate,[1,0,0],9
      state.oCurrentTopModel->rotate,[0,1,0],-3   ;correct the deviation of Y axis
      ;clockwise rotation relative to X axis
    end

    4:begin
      ;counterclockwise rotation relative to Z axis
      state.oCurrentTopModel->rotate,[0,0,1],9
      ;counterclockwise rotation relative to Z axis
    end

    5:begin
      ;clockwise rotation relative to Z axis
      state.oCurrentTopModel->rotate,[0,0,1],-9
      ;clockwise rotation relative to Z axis
    end

    6:begin
      ;counterclockwise rotation relative to Y axis
      state.oCurrentTopModel->rotate,[0,1,0],-9
      state.oCurrentTopModel->rotate,[1,0,0],-3   ;correct the deviation of X axis
      ;counterclockwise rotation relative to Y axis
    end

    7:begin
      ;clockwise rotation relative to Y axis
      state.oCurrentTopModel->rotate,[0,1,0],9
      state.oCurrentTopModel->rotate,[1,0,0],3   ;correct the deviation of X axis
      ;clockwise rotation relative to Y axis
    end

    8:begin
      ;flip horizontal
      state.oCurrentTopModel->rotate,[0,1,0],180
      state.oCurrentTopModel->rotate,[0,0,1],-40    ;correct the deviation of Z axis
      ;flip horizontal
    end

    9:begin
      state.oCurrentTopModel->translate,0,0,0.1   ;translate towards front
    end

    10:begin
      state.oCurrentTopModel->translate,0,0,-0.1   ;translate towards backward
    end

    11:begin
      state.oCurrentTopModel->translate,0.1,-0.035,0     ;translate towards left  （+ correcting the deviation of Y axis）
    end

    12:begin
      state.oCurrentTopModel->translate,-0.1,0.035,0     ;translate towards right  （+ correcting the deviation of Y axis）
    end

    13:begin
      state.oCurrentTopModel->translate,0.035,0.1,0     ;translate towards down  （+ correcting the deviation of X axis）
    end

    14:begin
      state.oCurrentTopModel->translate,-0.035,-0.1,0     ;translate towards up  （+ correcting the deviation of X axis）
    end
    else: begin

    end
  endcase
  state.win->Draw, state.oCurrentView
  d_objworld2NewMode, state, action
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end

;author:  THL,    2020/10/7
;update current model's state which showed in current window according to the 'key' and 'action'
;the 'model' consists of four sub-models which located in 'west','north','east','south' respectively
;action:  0-translate,1-rotate,2-scale
pro updateFourModelsByGestures,ev,key,action
  if key EQ -1 then return
  widget_control, ev.top, get_uvalue=state, /no_copy
  state.oCurrentView->GetProperty, uvalue=view_uval
  for i =0,3  do begin
    case key of
      0:begin ;zoom out
        (*(state.model_lists[view_uval.num]))[i]->scale,0.9,0.9,0.9,/PREMULTIPLY
      end

      1:begin ;zoom in
        (*(state.model_lists[view_uval.num]))[i]->scale,1.1,1.1,1.1, /PREMULTIPLY
      end

      2:begin ;counterclockwise rotation relative to X axis
        (*(state.model_lists[view_uval.num]))[i]->rotate,[1,0,0],9, /PREMULTIPLY
      end

      3:begin ;clockwise rotation relative to X axis
        (*(state.model_lists[view_uval.num]))[i]->rotate,[1,0,0],-9, /PREMULTIPLY
      end

      4:begin ;counterclockwise rotation relative to Z axis
        (*(state.model_lists[view_uval.num]))[i]->rotate,[0,1,0],-9, /PREMULTIPLY
      end

      5:begin ;clockwise rotation relative to Z axis
        (*(state.model_lists[view_uval.num]))[i]->rotate,[0,1,0],9, /PREMULTIPLY
      end

      6:begin ;counterclockwise rotation relative to Y axis
        (*(state.model_lists[view_uval.num]))[i]->rotate,[0,0,1],9, /PREMULTIPLY
      end

      7:begin ;clockwise rotation relative to Y axis
        (*(state.model_lists[view_uval.num]))[i]->rotate,[0,0,1],-9, /PREMULTIPLY
      end

      8:begin ;flip horizontal
        (*(state.model_lists[view_uval.num]))[i]->rotate,[0,0,1],180, /PREMULTIPLY
      end

      9:begin ;translate towards front
        (*(state.model_lists[view_uval.num]))[i]->translate,0,-0.1,0, /PREMULTIPLY
      end

      10:begin  ;translate towards backward
        (*(state.model_lists[view_uval.num]))[i]->translate,0,0.1,0, /PREMULTIPLY
      end

      11:begin  ;translate towards left
        (*(state.model_lists[view_uval.num]))[i]->translate,-0.05,0,0, /PREMULTIPLY
      end

      12:begin  ;translate towards right
        (*(state.model_lists[view_uval.num]))[i]->translate,0.05,0,0, /PREMULTIPLY
      end

      13:begin  ;translate towards down
        (*(state.model_lists[view_uval.num]))[i]->translate,0,0,-0.0175,/PREMULTIPLY
      end

      14:begin  ;translate towards up
        (*(state.model_lists[view_uval.num]))[i]->translate,0,0,0.0175,/PREMULTIPLY
      end
      else: begin ;other

      end
    endcase
  endfor
  state.win->Draw, state.oCurrentView
  d_objworld2NewMode, state, action
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end

FUNCTION getKeyAndActionByGesturesClass, gc
  ;action: 0-translate,1-rotate,2-scale
  key=-1
  action=-1
  case gc of
    '9' OR '13':begin  ;9、13  zoom out
      key=0
      action=2
    end
    '8' OR '12':begin  ;8、12 zoom in
      key=1
      action=2
    end
    '21':begin  ;21 counterclockwise rotation relative to X axis
      key=2
      action=1
    end
    '20':begin  ;20 clockwise rotation relative to X axis
      key=3
      action=1
    end
    '11' OR '15':begin  ;11、15  counterclockwise rotation relative to Z axis
      key=4
      action=1
    end
    '10' OR '14':begin  ;10、14  clockwise rotation relative to Z axis
      key=5
      action=1
    end
    '31':begin  ;31 counterclockwise rotation relative to Y axis
      key=6
      action=1
    end
    '30':begin  ;30 clockwise rotation relative to Y axis
      key = 7
      action = 1
    end
    '56':begin  ;56  flip horizontal
      key=8
      action=1
    end
    '5' OR '62':begin   ;5、62   translate towards front
      key=9
      action=0
    end
    '6':begin   ;6   translate towards backward
      key=10
      action=0
    end
    '2':begin   ;2   translate towards left
      key=11
      action=0
    end
    '1':begin   ;1   translate towards right
      key=12
      action=0
    end
    '3':begin   ;3   translate towards down
      key=13
      action=0
    end
    '4':begin   ;4   translate towards up
      key=14
      action=0
    end
    else: begin

    end
  endcase
  return, [key,action]
END

;------------------------------------------------------------------------------
;THL,
;SocketListenerCallback
Pro SocketListenerCallback,ID,H
  Status = File_Poll_Input(H['listenerLUN'], Timeout = .1)
  If (Status) then Begin
    catch, error_status
    if error_status ne 0 then begin
      CATCH, /CANCEL
      PRINT, 'An error occured while reading socket message'
      close,ClientLUN
      free_lun,ClientLUN
    endif
    Socket, ClientLUN, Accept = H['listenerLUN'], /Get_LUN,Connect_Timeout = 30., Read_Timeout = 30.
    da=''
    readf,ClientLUN,da
    print,'recived socket message：'+da
    res=getKeyAndActionByGesturesClass(da)
    updateFourModelsByGestures,H['ev'],res[0],res[1]
    close,ClientLUN
    free_lun,ClientLUN
  Endif
  !null = Timer.Set(.1, 'SocketListenerCallback',H)
end
;'gesture control' code end until here
;---------------------------------------------------------------------------



;----------------------------------------------------------------------------
pro d_objworld2Event, ev

COMMON GJPOLY, smoothedOutverts, Outconn1, threshold_min_value, threshold_max_value

    widget_control, ev.top, get_uvalue=state, /no_copy
   if n_elements(state) then demo_record, ev, $
    filename=state.record_to_filename, $
    cw=state.wModelModeRadio
  widget_control, ev.top, set_uvalue=state, /no_copy
  ;
  ;Shutdown.
  ;
  if tag_names(ev, /structure_name) eq 'WIDGET_KILL_REQUEST' then begin
    ;; first destroy any objects held in the uvalues of oModels
    widget_control, ev.top, get_uvalue=state, /no_copy
    if n_elements(state) eq 0 then return
    FOR i=0,state.highest_view_count-1 DO BEGIN
      IF ptr_valid(state.model_lists[i]) THEN BEGIN
        modelarr = *state.model_lists[i]
        FOR j=0,n_elements(modelarr)-2 DO BEGIN
          IF obj_valid(modelarr[j]) THEN BEGIN
            modelarr[j]->getproperty,uvalue=uval
            IF obj_valid(uval) THEN obj_destroy,uval
          ENDIF
        ENDFOR
      ENDIF
    ENDFOR
    widget_control, ev.top, set_uvalue=state, /no_copy
    ;; destroy top level base
    widget_control,ev.top,/destroy
    return
  end
  ;
  ;If mouse buttons are down, only process draw widget events.
  ;
  widget_control, ev.top, get_uvalue=state, /no_copy
  widget_control, ev.id, get_uval=uval
  if n_elements(state) eq 0 then return
  if state.btndown eq 1 then begin
    if uval[0] eq 'DRAW' then begin
      if ev.type eq 0 then begin ; Button down event?...
        widget_control, ev.top, set_uvalue=state, /no_copy
        return ; ...ignore it.  A mouse button is already down.
      end
    end $
    else begin
      widget_control, ev.top, set_uvalue=state, /no_copy
      return
    end
  end
  widget_control, ev.top, set_uvalue=state, /no_copy
  ;
  ;Normal event handling.
  ;
  case uval[0] of
    'WEBSITE': BEGIN
      url = 'http://web.xidian.edu.cn/gjia/index.html'
      url2 = 'http://web.xidian.edu.cn/gjia/en/research.html'
      mg_open_url, url2
    END    
    'AboutAIMIS3D' : BEGIN
      ;
      ; Create a non-modal text widget with information about
      ; the code used in the demo.
      ;
      IF (NOT XREGISTERED('DemoHelp')) THEN BEGIN
        
        AIMIS3DText = [ $
          '', $
          '  MPI+ Software', $
          '  AI-based Medical Image Segmentation for 3D Pringing and Visualization', $
          '  Version: 1.0-beta (2020-5-19)', $
          '', $
          '  Author: Guang Jia, PhD, DABR, Professor', $
          '  Xidian University, School of Computer Science & Technology', $
          '  E-mail: gjia@xidian.edu.cn; 779574877@qq.com', $
          '  Wechat: guangjia1976, QQ: 779574877', $
          '  Skype: guangjia1976, Facebook: jia_11_osu@hotmail.com', $
          '  Phone: (+86)-13609118250', $
          "  Address: No. 2 South Taibai Road", $
          "  Xi'an, Shaanxi 710071, China", $ 
          '', $
          '  Supported by:', $
          '  Shaanxi Xinweitai Biological Technology Co., Ltd.', $
          "  Xi'an Key Lab of Big Data and Intelligient Vision", $
          "  Xi'an AIer Kelly Electronic Technology Co., Ltd.", $
          '  Tangdu Hospital Affiliated to Air Force Medical University', $
          "  Shaanxi Provincial People's Hospital", $
          "--------------------------------------------------------------------------", $
          "| In memory of our beloved friend and genius scientist Dr. Jian Z. Wang! |", $
          "--------------------------------------------------------------------------", $
          '']
        TextTLB = WIDGET_BASE(GROUP_LEADER = ev.Top, /COLUMN)
        WIDGET_CONTROL, TextTLB, TLB_SET_TITLE = 'About MPI+'
        T = WIDGET_TEXT(TextTLB, XSIZE = 90, YSIZE = 20, $
            VALUE = AIMIS3DText, /SCROLL)
          WIDGET_CONTROL, TextTLB, /REALIZE
      ENDIF
    END  
  
    'QUIT' : begin
      ;; first destroy any objects held in the uvalues of oModels
      widget_control, ev.top, get_uvalue=state, /no_copy
      FOR i=0,state.highest_view_count-1 DO BEGIN
        IF ptr_valid(state.model_lists[i]) THEN BEGIN
          modelarr = *state.model_lists[i]
          FOR j=0,n_elements(modelarr)-2 DO BEGIN
            IF obj_valid(modelarr[j]) THEN BEGIN
              modelarr[j]->getproperty,uvalue=uval
              IF obj_valid(uval) THEN obj_destroy,uval
            ENDIF
          ENDFOR
        ENDIF
      ENDFOR
      widget_control, ev.top, set_uvalue=state, /no_copy
      ;; destroy top level base
      widget_control, ev.top, /destroy
      return
    end
    
    ;GJ, 2019/2/7
    ;allow user to select different language
    'LANGUAGE': begin
      widget_control, ev.top, get_uvalue=state, /no_copy
;      languageind = where(state.language_subjects eq uval[1])
      state.icon_dir = './lib/icon_pics/'+uval[1]+'/'
      print, 'language = ', state.icon_dir
      ;store this selection to the language_file
      language_file = './lib/icon_pics/language_file.txt'
      IF FILE_TEST(language_file, /WRITE) THEN write_language_file,language_file,state.icon_dir
      
      WIDGET_CONTROL, state.wDCMsort, SET_VALUE=state.icon_dir+'dcmsort.bmp', /bitmap
      WIDGET_CONTROL, state.wDCMsum, SET_VALUE=state.icon_dir+'dcmsum.bmp', /bitmap
      WIDGET_CONTROL, state.wDCMone, SET_VALUE=state.icon_dir+'dcmone.bmp', /bitmap
      WIDGET_CONTROL, state.wDCMpict, SET_VALUE=state.icon_dir+'dcmpict.bmp', /bitmap
      WIDGET_CONTROL, state.wGestureControlButton, SET_VALUE=state.icon_dir+'gesture_control.bmp', /BITMAP   ;THL,2019/11/12
      WIDGET_CONTROL, state.wVoiceControlButton, SET_VALUE=state.icon_dir+'voice_control.bmp', /BITMAP   ;THL,2020/06/16
      WIDGET_CONTROL, state.wRotateRadio, SET_VALUE=state.icon_dir+'rotate.bmp', /bitmap
      WIDGET_CONTROL, state.wTranslateRadio, SET_VALUE=state.icon_dir+'translate.bmp', /bitmap
      WIDGET_CONTROL, state.wScaleRadio, SET_VALUE=state.icon_dir+'scale.bmp', /bitmap
      WIDGET_CONTROL, state.wOpenButton, SET_VALUE=state.icon_dir+'open.bmp', /BITMAP
      for i=0, n_elements(state.wOpen_listButton)-1 do WIDGET_CONTROL, state.wOpen_listButton[i], SET_VALUE=state.icon_dir+state.open_subjects[i]+'.bmp', /BITMAP
      WIDGET_CONTROL, state.wMyAddFolderButton, SET_VALUE=state.icon_dir+'addfolder.bmp', /BITMAP
      WIDGET_CONTROL, state.wSplitFive, SET_VALUE=state.icon_dir+'splitfive.bmp', /BITMAP
      WIDGET_CONTROL, state.wAddButton, SET_VALUE=state.icon_dir+'add.bmp', /bitmap
      for i=0, n_elements(state.wAdd_listButton)-1 do WIDGET_CONTROL, state.wAdd_listButton[i], SET_VALUE=state.icon_dir+state.addgroup_subjects[i]+'.bmp', /BITMAP
    ;  WIDGET_CONTROL, state.wAddChildButton, SET_VALUE=state.icon_dir+'addchild.bmp', /BITMAP
      WIDGET_CONTROL, state.wModelDeleteButton, SET_VALUE=state.icon_dir+'delete.bmp', /BITMAP
      WIDGET_CONTROL, state.wModelSelectButton, SET_VALUE=state.icon_dir+'select.bmp', /BITMAP
      WIDGET_CONTROL, state.wUnselectButton, SET_VALUE=state.icon_dir+'unselect.bmp', /BITMAP
      WIDGET_CONTROL, state.wModifyButton, SET_VALUE=state.icon_dir+'modifyselected.bmp', /BITMAP
      WIDGET_CONTROL, state.wModifyAllButton, SET_VALUE=state.icon_dir+'modifyall.bmp', /BITMAP
      WIDGET_CONTROL, state.wMyResetButton, SET_VALUE=state.icon_dir+'reset.bmp', /BITMAP

      widget_control, ev.top, set_uvalue=state, /no_copy
    end
    
    ;GJ and YFZ, 2022/5/16
    ;MPI simulation based on sine excitation with relaxation time
    'SinExcitation': begin
      Sin_Excitation
    end
    
    'Sin_Excitation_Chebyshev_MPI': begin
      Sin_Excitation_Chebyshev_MPI
    end
    
    'PulsedExcitation': begin
      Pulse_Excitation
    end
    
    ;GJ, 2021/7/27
    ;MPI simulation for relaxation kernel
    'RelaxationKernel': begin
      mpi_filter
      
    end
    
    'Anderson': begin
     hk_Anderson_coil
    end
    
    ;GJ, 2021/7/28
    ;MPI simulation for different H
    'Langevin': begin
      
      Langevin_derivative_Function
    end
    
    ;GJ, 2021/7/28
    ;MPI simulation for different H
    'Signal_H': begin

      SingalResponse_7point_different_driveAmp
    end
    
    ;GJ, 2021/7/28
    ;MPI simulation for MPT 2dV
    'MPI_FBP_Sim': begin

      aimis3d_MPT_2d_V_field
    end
    
    ;GJ, 2022/6/8
    ;MRT dat file renaming and analyzing
    'MPS_Data_Analyze': begin
      dataFileRenaming_widgets
    end
    
    ;GJ, 2023/8/5
    ;MRT dat file renaming and analyzing
    'Twin_DCM_Analyze': begin
      twin_pattern_hit_MLM_localization
    end
    
    'VRML' : begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      if (state.oCurrentView NE obj_new()) then begin
        file = dialog_pickfile( $
          /write, $
          file='untitled.wrl', $
          group=ev.top, $
          filter='*.wrl' $
          )
        if (file NE '') then begin
          widget_control, /hourglass
          state.win->GetProperty, $
            dimension=wdims, $
            resolution=res,$
            color_model=cm, $
            n_colors=icolors
          oVRML = obj_new('IDLgrVRML', $
            dimensions=wdims, $
            resolution=res, $
            color_model=cm, $
            n_colors=icolors $
            )
          oVRML->setproperty, filename=file
          demo_draw, oVRML, state.oCurrentView, debug=state.debug
          obj_destroy,oVRML
        end
      end
      widget_control, ev.top, set_uvalue=state, /no_copy
    end
    'CLIPBOARD' : begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      state.win->GetProperty, $
        dimension=wdims, $
        resolution=res, $
        color_model=cm, $
        n_colors=icolors
      oClipboard = obj_new('IDLgrClipboard', $
        dimensions=wdims, $
        resolution=res, $
        color_model=cm, $
        n_colors=icolors $
        )
      demo_draw, oClipboard, state.scene, debug=state.debug
      obj_destroy, oClipboard
      widget_control, ev.top, set_uvalue=state, /no_copy
    end
    'PRINT' : begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      oPrinter = obj_new('IDLgrPrinter')
      if (dialog_printersetup(oPrinter)) then begin
        if (dialog_printjob(oPrinter)) then begin
          oPrinter->GetProperty,resolution=res
          DPI = 2.54/float(res)
          state.win->GetProperty,resolution=res
          DPI = 2.54/float(res)

          ;Hack, swap from pixels to inches for views...
          state.win->GetProperty, dimension=wdims
          oViews = state.scene->get(/all)
          for i=0,n_elements(oViews)-1 do begin
            oViews[i]->IDLgrView::getproperty,$
              loc=loc,dim=vdim
            loc = loc/DPI
            vdim = vdim/DPI
            oViews[i]->IDLgrView::setproperty,$
              loc=loc, dim=vdim, units=1
          end

          ;...PRINT!...
          demo_draw, oPrinter, state.scene, debug=state.debug
          oPrinter->newdocument

          ;...and back to pixels
          for i=0,N_ELEMENTS(oViews)-1 do begin
            oViews[i]->IDLgrView::getproperty,$
              loc=loc,dim=vdim
            loc = loc*DPI
            vdim = vdim*DPI
            oViews[i]->IDLgrView::setproperty,$
              loc=loc,dim=vdim,units=0
          end

        end
      end
      obj_destroy,oPrinter
      widget_control, ev.top, set_uvalue=state, /no_copy
    end
    'AA' : begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      d_objworld2SceneAntiAlias, state.win, state.scene, 8
      widget_control, ev.top, set_uvalue=state, /no_copy
    end
    'LOAD' : begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      if ((state.selected ne obj_new()) and $
        (state.cur_tool ne 1)) then begin
      wDel=state.wMyResetButton
        ;wMyResetButton
        
        ;myWorldAdd,state
      end
      DB_ADD,state=state
      print, 'vol_HU_cube_resolution: ', state.vol_HU_cube_resolution
      print, 'patient age 5 = ', state.patient_age
       widget_control, ev.top, set_uvalue=state, /no_copy
      d_objworld2Event,{ $
        top: ev.top, $
        handler: 0L, $
        id: wDel $
      }   
    end
    'DBLOAD':begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      DB_LOAD,state
      ;GJ, 2018/12/11, if select nothing, do nothing
      IF ~PTR_VALID(state.vol_HU_cube) THEN BEGIN
        widget_control, ev.top, set_uvalue=state
        return
      ENDIF
      
      ;GJ, 2019/6/5, after loading the object, hiding baseplate with company icons
      state.obaseplate->SetProperty, hide = 1


      ;GJ 2020/5/19, updating folder
      temp_dir=*state.additionalDir
      widget_control, state.text1, set_value='folder = '+temp_dir[N_ELEMENTS(temp_dir)-1],/APPEND

      
;      print, 'patient age 6 = ', state.patient_age
      ;myWorldAdd,state 
       wDel=state.wMyResetButton
       widget_control, state.wMyAddFolderButton, sensitive=1
       widget_control, state.wModelModeRadio, sensitive = 1
       widget_control, state.wModifyButton,  sensitive = 1
       widget_control, state.wModifyAllButton, sensitive = 0
       
;       widget_control, state.wModelModeRadio, get_value=mode
       mode = 1
       widget_control, state.wRotateRadio, SET_BUTTON = 1
       state.oModelMan->SetProperty, mode=mode
       d_objworld2NewMode, state, mode
       widget_control, ev.top, set_uvalue=state
     IF(PTR_VALID(state.vol_HU_cube)) THEN BEGIN
      ;GJ, 2018/12/11, renewing wdraw list
      PTR_FREE, ((state.Out)[0])
      d_objworld2Event,{ $
          top: ev.top, $
          handler: 0L, $
          id: wDel $
          }
     ENDIF
   end
   
   ;GJ 2019/2/13 loading imges from CD or U disk
   'CDU':begin
     widget_control, ev.top, get_uvalue=state, /no_copy
  ;   DB_LOAD,state
     afn_disc = DIALOG_PICKFILE(PATH='E:\',/MUST_EXIST, /DIRECTORY, TITLE='Please select a disc:')
     afn = file_dirname(afn_disc)
     filearr=file_search(afn_disc,'00*.dcm',count=num)
     IF num EQ 0 THEN BEGIN
       afn_des = './data_temp/'
       filearr_des=file_search(afn_des,'*.dcm',count=num_des)
       IF num_des GT 0 THEN FILE_DELETE, afn_des, /RECURSIVE
       Separate_to_series, afn_disc, afn_des, afn_des_output = afn_des_output
       IF N_ELEMENTS(afn_des_output) LE 1 THEN BEGIN
         widget_control, ev.top, set_uvalue=state
         RETURN
       ENDIF
       afn = DIALOG_PICKFILE(PATH=afn_des_output[1],/MUST_EXIST, /DIRECTORY, TITLE='Select Series folder:')
       filearr=file_search(afn,'00*.dcm',count=num)
       IF num EQ 0 THEN BEGIN
         widget_control, ev.top, set_uvalue=state
         RETURN
       ENDIF
       image_read, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal
       ;确保vol_HU_iso不为空
       IF((size(vol_HU_cube))[0] LE 2) THEN BEGIN
         widget_control, ev.top, set_uvalue=state
         RETURN
       ENDIF
       file_delete, afn_des, /RECURSIVE, /QUIET
     ENDIF ELSE BEGIN
       image_read, afn_disc, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal
       ;确保vol_HU_iso不为空
       IF((size(vol_HU_cube))[0] LE 2) THEN BEGIN
         widget_control, ev.top, set_uvalue=state
         RETURN
       ENDIF
     ENDELSE
     ptr_free,state.vol_HU_cube
     state.vol_HU_cube = ptr_new(vol_HU_cube)
     state.vol_HU_cube_resolution = vol_HU_cube_resolution
     state.patient_age = patient_age
     state.additionalDir = ptr_new(afn)
     ;GJ 2020/5/19, updating folder
     temp_dir=*state.additionalDir
     widget_control, state.text1, set_value='folder = '+temp_dir[N_ELEMENTS(temp_dir)-1],/APPEND
     state.p_modal = p_modal
     print, 'patient age 2 = ', state.patient_age
     print, 'database resolution', state.vol_HU_cube_resolution
   
     ;GJ, 2018/12/11, if select nothing, do nothing
     IF ~PTR_VALID(state.vol_HU_cube) THEN BEGIN
       widget_control, ev.top, set_uvalue=state
       return
     ENDIF
     
     ;GJ, 2019/6/5, after loading the object, hiding baseplate with company icons
     state.obaseplate->SetProperty, hide = 1
     
     wDel=state.wMyResetButton
     widget_control, state.wMyAddFolderButton, sensitive=1
     widget_control, state.wModelModeRadio, sensitive = 1
     widget_control, state.wModifyButton,  sensitive = 1
     widget_control, state.wModifyAllButton, sensitive = 0
  
     mode = 1
     widget_control, state.wRotateRadio, SET_BUTTON = 1
     state.oModelMan->SetProperty, mode=mode
     d_objworld2NewMode, state, mode
     widget_control, ev.top, set_uvalue=state
     IF(PTR_VALID(state.vol_HU_cube)) THEN BEGIN
       ;GJ, 2018/12/11, renewing wdraw list
       PTR_FREE, ((state.Out)[0])
       d_objworld2Event,{ $
         top: ev.top, $
         handler: 0L, $
         id: wDel $
       }
     ENDIF
   end

    
     'SORT':begin
     widget_control, ev.top, get_uvalue=state, /no_copy
     afn = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/',/MUST_EXIST, TITLE="file choose", /DIRECTORY)
     IF STRLEN(afn) NE 0 THEN Separate_to_series, afn, afn_des
     wDel=state.wMyResetButton
     widget_control, ev.top, set_uvalue=state, /no_copy
     
;     d_objworld2Event,{ $
;       top: ev.top, $
;       handler: 0L, $
;       id: wDel $
;     }
   end
   
   'DCMSUM':begin
   widget_control, ev.top, get_uvalue=state, /no_copy
   Crawler
   wDel=state.wMyResetButton
   widget_control, ev.top, set_uvalue=state, /no_copy

   ;     d_objworld2Event,{ $
   ;       top: ev.top, $
   ;       handler: 0L, $
   ;       id: wDel $
   ;     }
   end
   
   'DCMONE':begin
   widget_control, ev.top, get_uvalue=state, /no_copy
   copy_to_one_direct
   wDel=state.wMyResetButton
   widget_control, ev.top, set_uvalue=state, /no_copy

   ;     d_objworld2Event,{ $
   ;       top: ev.top, $
   ;       handler: 0L, $
   ;       id: wDel $
   ;     }
   end
   
   'DCMPICT':begin
   widget_control, ev.top, get_uvalue=state, /no_copy
   convert_to_bmp_images
   wDel=state.wMyResetButton
   widget_control, ev.top, set_uvalue=state, /no_copy
  
   ;     d_objworld2Event,{ $
   ;       top: ev.top, $
   ;       handler: 0L, $
   ;       id: wDel $
   ;     }
   end
   
   ;THL,  2019/11/13  call the gesture_control.py
   'GESTURECONTROL':begin
     WORK_DIR=file_dirname(routine_filepath()) ;current working directory
     
     spawn,'python '+WORK_DIR+'/gesture_control.py',/nowait;,/hide
     widget_control, ev.top, get_uvalue=state, /no_copy
     widget_control, state.wGestureControlButton, sensitive = 0
     widget_control, ev.top, set_uvalue=state, /no_copy
     print,'the gesture control has been started successfully'

     port=6667
     SOCKET, listenerLUN, port, /GET_LUN, /LISTEN
     !null = Timer.Set (.1, 'SocketListenerCallback',Hash('ev', ev,'listenerLUN',listenerLUN))
   end

   ;THL,  2020/06/16  call the voice_control.py
   'VOICECONTROL':begin
     WORK_DIR=file_dirname(routine_filepath()) ;current working directory

     spawn,'python '+WORK_DIR+'\voice_control.py',/nowait,/hide
     widget_control, ev.top, get_uvalue=state, /no_copy
     widget_control, state.wVoiceControlButton, sensitive = 0
     widget_control, ev.top, set_uvalue=state, /no_copy
     print,'the voice control has been started successfully'

     port=6666
     SOCKET, listenerLUN, port, /GET_LUN, /LISTEN
     !null = Timer.Set (.1, 'SocketListenerCallback',Hash('ev', ev,'listenerLUN',listenerLUN))

   end
   
   
   'CRAWLER':begin
   widget_control, ev.top, get_uvalue=state, /no_copy
   Crawler
   wDel=state.wMyResetButton
   widget_control, ev.top, set_uvalue=state, /no_copy

   ;     d_objworld2Event,{ $
   ;       top: ev.top, $
   ;       handler: 0L, $
   ;       id: wDel $
   ;     }
 end

   'ADDFOLDER':begin
 widget_control, ev.top, get_uvalue=state, /no_copy
 
 IF PTR_VALID(state.additionalDir) THEN BEGIN
   IF STRLEN((*(state.additionalDir))[0]) GT 5 THEN BEGIN
     afn = DIALOG_PICKFILE(PATH=(*(state.additionalDir))[0],/MUST_EXIST, TITLE="additional folder choose", /DIRECTORY)
     IF STRLEN(afn) NE 0 THEN BEGIN
       n_folder = SIZE((*(state.additionalDir)), /dimension)
       samefolder = 0
       FOR i=0, n_folder[0]-1 DO BEGIN
         IF STRCMP((*(state.additionalDir))[i], afn) EQ 1 THEN samefolder = 1
       ENDFOR
       IF samefolder EQ 0 THEN BEGIN
         temp_dir = [*(state.additionalDir), afn]
         state.additionalDir = PTR_NEW(temp_dir)
         ;GJ 2020/5/19, updating folder
         widget_control, state.text1, set_value='folder = '+temp_dir[N_ELEMENTS(temp_dir)-1],/APPEND
       ENDIF
     ENDIF
     print, 'folders = ', *(state.additionalDir)
   ENDIF
 ENDIF
 
 wDel=state.wMyResetButton
 widget_control, ev.top, set_uvalue=state, /no_copy

 ;     d_objworld2Event,{ $
 ;       top: ev.top, $
 ;       handler: 0L, $
 ;       id: wDel $
 ;     }
 end
    
    'MODIFY':begin
    widget_control, ev.top, get_uvalue=state, /no_copy
    if state.selected ne obj_new() and PTR_VALID(state.info_tables) then begin
      info=*(state.info_tables)
      
      FOR m=0, N_ELEMENTS(info.str)-1 DO BEGIN
        IF STRCMP(state.selected.name,info.str[m]) EQ 1 THEN BEGIN
           HU_max_value=info.HU_max_value[m]
           HU_min_value=info.HU_min_value[m]
           select_mask_value = m
        ENDIF
      ENDFOR
      IF N_ELEMENTS(select_mask_Value) EQ 0 THEN BEGIN
        ;widget_control, state.text, set_value="No body part selection!"
        ;infowarning = DIALOG_MESSAGE('No body part picked! Please pick one!', /ERROR)
        ;after deleting all objects, then load new one
        GOTO, stlobjectmodify
      ENDIF ELSE BEGIN
        IF select_mask_value LE N_ELEMENTS(info.str)-1 THEN BEGIN
          oPolygon =state.selected->get(pos=0)
          oPolygon->getProperty,POLYGONS=POLYGONS
          oPolygon->getProperty,DATA=DATA
          oPolygon->getProperty,COLOR=COLOR
  
          ;define to the rotation center
          smoothedOutverts = DATA
          ma=max(smoothedOutverts)
          mi=min(smoothedOutverts)
          bias=(ma-mi)/2
          rl=mi+bias
          rr=ma+bias
          cc=[(-rl)/(rr-rl),1/(rr-rl)]
          xma = max((smoothedOutverts)[0,*])
          xmi = min((smoothedOutverts)[0,*])
          xmid = 0.5*(xma+xmi)*cc[1]
          yma = max((smoothedOutverts)[1,*])
          ymi = min((smoothedOutverts)[1,*])
          ymid = 0.5*(yma+ymi)*cc[1]
          zma = max((smoothedOutverts)[2,*])
          zmi = min((smoothedOutverts)[2,*])
          zmid = 0.5*(zma+zmi)*cc[1]
          sel_one={name:state.selected.name,$
            p_verts:PTR_NEW(p_verts),$
            p_conn:PTR_NEW(p_conn),$
            verts:DATA,$
            conn:POLYGONS, $
            color:COLOR,$
            xcoord_conv:[-xmid, cc[1]], $
            ycoord_conv:[-ymid, cc[1]], $
            zcoord_conv:[-zmid, cc[1]], $
            HU_min_value:HU_min_value, $
            HU_max_value:HU_max_value $
          }

          
          ;define the struc_modifyobject
          IF N_ELEMENTS(struc_modifyobject) EQ 0 THEN BEGIN
            ;; create the struc_modifyobject structure
            cur_patient=0
            ;; create a struc_modifyobject structure
            struc_modifyobject = {struc_modifyobject}
          ENDIF ELSE BEGIN
            n=N_ELEMENTS(struc_modifyobject)
            p=replicate({struc_modifyobject},n+1)
            p[0:n-1]=struc_modifyobject[0:n-1]
            cur_patient=n
            struc_modifyobject=p
          ENDELSE
          ;GJ, 2019/2/19
          ;; initialize the elements for THIS struc_modifyobject
          struc_modifyobject[cur_patient].modify_step = 0
          struc_modifyobject[cur_patient].current_step_index = 0
          struc_modifyobject[cur_patient].exist_num[0] = 1
          print, 'hide status = ', struc_modifyobject[cur_patient].hide_status
          struc_modifyobject[cur_patient].hide_status = 1
          struc_modifyobject[cur_patient].hide_status[0, 0] = 0
          color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
          struc_modifyobject[cur_patient].p_ori_verts_list[0] = PTR_NEW(DATA)
          ;define a new polygon
          vert_colors = DATA
          vert_colors[0,*] = color_rgb[1, 0]
          vert_colors[1,*] = color_rgb[1, 1]
          vert_colors[2,*] = color_rgb[1, 2]
          oPoly = OBJ_NEW('IDLgrPolygon', DATA, NAME='0', POLYGON=POLYGONS, STYLE=2, $
            SHADING=1, vert_colors=vert_colors, $
            XCOORD_CONV=sel_one.xcoord_conv, YCOORD_CONV=sel_one.ycoord_conv, ZCOORD_CONV=sel_one.zcoord_conv)
          struc_modifyobject[cur_patient].polygon_list[0] = PTR_NEW(oPoly)
          
          print, 'vol_HU_cube_resolution: ', state.vol_HU_cube_resolution
          IF PTR_VALID(state.growmaxfive_mask) THEN BEGIN
            ;; initialize the elements for THIS struc_modifyobject
            index = WHERE(*(state.growmaxfive_mask) EQ (select_mask_value+1), count)
            IF count GT 1 THEN BEGIN
              struc_modifyobject[cur_patient].p_maskindex_list[0] = PTR_NEW(index)
              struc_modifyobject[cur_patient].countmask_list[0]=count
            ENDIF
            
            d_surfview1,sel_one, *(state.vol_HU_cube), *(state.growmaxfive_mask) EQ (select_mask_value+1), *(state.growmaxfive_mask) EQ (select_mask_value+1), state.vol_HU_cube_resolution, struc_modifyobject, cur_patient, state.icon_dir, additionalDir = *(state.additionalDir)
          ENDIF ELSE BEGIN
            ;; initialize the elements for THIS struc_modifyobject
            index = WHERE(*(state.vol_HU_cube_mask) EQ 1, count)
            IF count GT 1 THEN BEGIN
              struc_modifyobject[cur_patient].p_maskindex_list[0] = PTR_NEW(index)
              struc_modifyobject[cur_patient].countmask_list[0]=count
            ENDIF
            
            d_surfview1,sel_one, *(state.vol_HU_cube), *(state.vol_HU_cube_mask), *(state.vol_HU_cube_mask_modify), state.vol_HU_cube_resolution, struc_modifyobject, cur_patient, state.icon_dir, additionalDir = *(state.additionalDir)
          ENDELSE
          
          ;
        ENDIF ELSE BEGIN
          infowarning = DIALOG_MESSAGE('No volume picked! Please pick one!', /ERROR)
        ENDELSE
      ENDELSE
      widget_control, ev.top, set_uvalue=state, /no_copy
    endif else begin
stlobjectmodify:      IF state.selected ne obj_new() THEN BEGIN
        oPolygon =state.selected->get(pos=0)
;        print, 'is a polygon 1 yes, 0 no: ', (obj_isa(oPolygon,'IDLgrPolygon'))
;        oPolygon->getProperty,ALL=allpro
;        alltagnames = TAG_NAMES(allpro)
;;        print, 'all ', alltagnames
;;        data_ind = STRMATCH(alltagnames, 'DATA')
;;        print, 'data ', data_ind
;        poly_ind = STRMATCH(alltagnames, 'POLYGONS')
;;        print, 'oply ', poly_ind
;        IF TOTAL(poly_ind) EQ 1 THEN BEGIN
         IF obj_isa(oPolygon,'IDLgrPolygon') THEN BEGIN 
          oPolygon->getProperty,POLYGONS=POLYGONS_0
          oPolygon->getProperty,DATA=DATA_0
          oPolygon->getProperty,COLOR=COLOR

;          smoothedOutverts = MESH_SMOOTH(DATA, POLYGONS,  ITERATIONS=50)
          result = MESH_DECIMATE(DATA_0, POLYGONS_0, POLYGONS, vertices=DATA)
          ;rescale all data to 256*256*256 matrix
          ;GJ 2019/1/30
          max_DATA = MAX(DATA)
          min_DATA = MIN(DATA)
          range_DATA = max_DATA - min_DATA
          vol_HU_cube_resolution = (range_DATA + 10.) / 256.
          DATA_scaled = (DATA - min_DATA) * 256. / (range_DATA + 10.)
          ;
          ;define to the rotation center
          smoothedOutverts = DATA_scaled

          vol_HU_cube = DBLARR(256, 256, 256)*0.-100.
          vol_HU_cube_mask = BYTARR(256, 256, 256)*0.
          vol_HU_cube_mask_modify = INTARR(256, 256, 256)*0.
          FOR ii=0, (SIZE(smoothedOutverts, /DIM))[1]-1 DO BEGIN
;            vol_HU_cube[smoothedOutverts[0, ii], smoothedOutverts[1, ii], smoothedOutverts[2, ii]] = 255.
            vol_HU_cube_mask[smoothedOutverts[0, ii], smoothedOutverts[1, ii], smoothedOutverts[2, ii]] = 1
;            vol_HU_cube_mask_modify[smoothedOutverts[0, ii], smoothedOutverts[1, ii], smoothedOutverts[2, ii]] = 255
          ENDFOR
          
          ;GJ, 2019/4/21
          ; Create the shape operator:
          S = REPLICATE(1, 3, 3, 3)
          ; "Opening" operator:
          vol_HU_cube_mask = DILATE(vol_HU_cube_mask, S);, S)
          vol_HU_cube_mask = ERODE(vol_HU_cube_mask, S) GT 0
          
          radius = 1
          strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
          
          FOR i=0, 255 DO BEGIN
            dilateImg = DILATE(REFORM(vol_HU_cube_mask[*, *, i]), strucElem)
            index1 = WHERE(dilateImg GT 0, oneCount)
            index = WHERE(dilateImg EQ 0, zeroCount)
            IF zeroCount GT 10 AND oneCount GT 10 THEN BEGIN
              index1_2d = ARRAY_INDICES(dilateImg, index1)
              x_range = [MIN(index1_2d[0,*]), MAX(index1_2d[0,*])]
              y_range = [MIN(index1_2d[1,*]), MAX(index1_2d[1,*])]
              index2d = ARRAY_INDICES(dilateImg, index)
              FOR j=0, zeroCount-1 DO BEGIN
                x = index2d[0, j]
                y = index2d[1, j]
                IF x GT x_range[0] AND x LT x_range[1] AND y GT y_range[0] AND y LT y_range[1] THEN BEGIN
                  x_m = TOTAL(dilateImg[0:x, y]);+TOTAL(vol_HU_cube_mask[0:x, y, i-1:i+1])
                  x_p = TOTAL(dilateImg[x:*, y]);+TOTAL(vol_HU_cube_mask[x:*, y, i-1:i+1])
                  y_m = TOTAL(dilateImg[x, 0:y]);+TOTAL(vol_HU_cube_mask[x, 0:y, i-1:i+1])
                  y_p = TOTAL(dilateImg[x, y:*]);+TOTAL(vol_HU_cube_mask[x, y:*, i-1:i+1])
                  IF x_m*x_p*y_m*y_p NE 0 THEN dilateImg[x,y]=1
                ENDIF                  
              ENDFOR
            ENDIF
            
            vol_HU_cube_mask[*, *, i] = ERODE(DILATE(dilateImg, strucElem), strucElem) GT 0
          ENDFOR
          
          ;GJ 2020/06/24, remove the small holes
          FOR i=0, 255 DO BEGIN
            temp_img = REFORM(vol_HU_cube_mask[i, *, *])
            vol_HU_cube_mask[i, *, *] = ERODE(DILATE(temp_img, strucElem), strucElem)
          ENDFOR
          
          ;GJ 2020/06/24, remove the small holes@
          FOR j=0, 255 DO BEGIN
            temp_img = REFORM(vol_HU_cube_mask[*, j, *])
            vol_HU_cube_mask[*, j, *] = ERODE(DILATE(temp_img, strucElem), strucElem)
          ENDFOR
          
;          ;GJ 2020/6/26, output the data as txt file
;          fn='C:\D_drive\Tangdu_NCS\lll_20200528\dm.txt'
;          openw, lun, fn, /get_lun, width=1000
;          FOR i=0, 255 DO BEGIN
;            printf, lun, vol_HU_cube_mask[*, *, i]
;          ENDFOR
;          free_lun, lun
          
          ;GJ, 2020/6/24, to convert stl file to IDL 3d objects
          SHADE_VOLUME, vol_HU_cube_mask,0.9, v, p          
          smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
          smoothedOutverts = smoothedVertices
          POLYGONS = p
          ma=max(smoothedOutverts)
          mi=min(smoothedOutverts)
          bias=(ma-mi)/2
          rl=mi+bias
          rr=ma+bias
          cc=[(-rl)/(rr-rl),1/(rr-rl)]
          xma = max((smoothedOutverts)[0,*])
          xmi = min((smoothedOutverts)[0,*])
          xmid = 0.5*(xma+xmi)*cc[1]
          yma = max((smoothedOutverts)[1,*])
          ymi = min((smoothedOutverts)[1,*])
          ymid = 0.5*(yma+ymi)*cc[1]
          zma = max((smoothedOutverts)[2,*])
          zmi = min((smoothedOutverts)[2,*])
          zmid = 0.5*(zma+zmi)*cc[1]
          sel_one={name:'stl',$
            p_verts:PTR_NEW(smoothedOutverts),$
            p_conn:PTR_NEW(POLYGONS),$
            verts:smoothedOutverts,$
            conn:POLYGONS, $
            color:COLOR,$
            xcoord_conv:[-xmid, cc[1]], $
            ycoord_conv:[-ymid, cc[1]], $
            zcoord_conv:[-zmid, cc[1]], $
            HU_min_value:128, $
            HU_max_value:255 $
          }
          
          index = WHERE(vol_HU_cube_mask GT 0, count)
          IF count GT 1 THEN BEGIN
            vol_HU_cube[index] = 255
            vol_HU_cube_mask_modify[index] = 255
          ENDIF
;          vol_HU_cube = ERODE(DILATE(vol_HU_cube, S, /GRAY), S, /GRAY)
;          vol_HU_cube_mask_modify = ERODE(DILATE(vol_HU_cube_mask_modify, S, /GRAY), S, /GRAY)
          
          addDIR = PTR_NEW('c:\')
          
          ;define the struc_modifyobject
          IF N_ELEMENTS(struc_modifyobject) EQ 0 THEN BEGIN
            ;; create the struc_modifyobject structure
            cur_patient=0
            ;; create a struc_modifyobject structure
            struc_modifyobject = {struc_modifyobject}
          ENDIF ELSE BEGIN
            n=N_ELEMENTS(struc_modifyobject)
            p=replicate({struc_modifyobject},n+1)
            p[0:n-1]=struc_modifyobject[0:n-1]
            cur_patient=n
            struc_modifyobject=p
          ENDELSE
          ;GJ, 2019/2/19
          ;; initialize the elements for THIS struc_modifyobject
          struc_modifyobject[cur_patient].modify_step = 0
          struc_modifyobject[cur_patient].current_step_index = 0
          struc_modifyobject[cur_patient].exist_num[0] = 1
          struc_modifyobject[cur_patient].hide_status = 1
          struc_modifyobject[cur_patient].hide_status[0, 0] = 0
          color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
          struc_modifyobject[cur_patient].p_ori_verts_list[0] = PTR_NEW(smoothedOutverts)
          ;define a new polygon
          vert_colors = smoothedOutverts
          vert_colors[0,*] = color_rgb[1, 0]
          vert_colors[1,*] = color_rgb[1, 1]
          vert_colors[2,*] = color_rgb[1, 2]
          oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME='0', POLYGON=POLYGONS, STYLE=2, $
            SHADING=1, vert_colors=vert_colors, $
            XCOORD_CONV=sel_one.xcoord_conv, YCOORD_CONV=sel_one.ycoord_conv, ZCOORD_CONV=sel_one.zcoord_conv)
          struc_modifyobject[cur_patient].polygon_list[0] = PTR_NEW(oPoly)
          
          index = WHERE(vol_HU_cube_mask EQ 1, count)
          IF count GT 1 THEN BEGIN
            struc_modifyobject[cur_patient].p_maskindex_list[0] = PTR_NEW(index)
            struc_modifyobject[cur_patient].countmask_list[0]=count
          ENDIF
          
          d_surfview1,sel_one, vol_HU_cube, vol_HU_cube_mask, vol_HU_cube_mask_modify, vol_HU_cube_resolution, struc_modifyobject, cur_patient, state.icon_dir, additionalDir = addDir
        ENDIF
        widget_control, ev.top, set_uvalue=state, /no_copy
      ENDIF ELSE BEGIN
        widget_control, state.text, set_value="No current selection!"
      ENDELSE
    endelse

  end

  'MODIFYALL': begin
    widget_control, ev.top, get_uvalue=state, /no_copy
    ;闁兼儳鍢茶ぐ鍣坥del_lists 闁告瑥锕弳杈ㄦ償閿燂拷
    state.oCurrentView->GetProperty, uvalue=view_uval
    count = n_elements( $
      *(state.model_lists[view_uval.num]) $
      )
    model_lists=*(state.model_lists[view_uval.num])
    state.oModelMan->SetTarget, obj_new()
    
    
    
;    print, 'n models = ', count
;    print, 'model_lists = ', model_lists
    IF ~PTR_VALID(state.info_tables) THEN BEGIN
      widget_control, ev.top, set_uvalue=state, /no_copy
      RETURN
    ENDIF  
    info=*(state.info_tables)
    IF count GE 9 THEN ncount = 5 ELSE ncount = count-3
    temp_maskind_list = PTRARR(ncount)
    temp_countmask_list = LONARR(ncount)
    FOR i=0,ncount-1 DO BEGIN
      model_lists[i]->GETproperty, all=allinfo
      
      FOR m=0, N_ELEMENTS(info.str)-1 DO BEGIN
        IF STRCMP(allinfo.name, info.str[m]) EQ 1 THEN BEGIN
          HU_max_value=info.HU_max_value[m]
          HU_min_value=info.HU_min_value[m]
          select_mask_value = m
          print, 'selected mask value = ', m
          IF PTR_VALID(state.growmaxfive_mask) AND N_ELEMENTS(temp_mask) EQ 0 THEN BEGIN
            temp_mask = *(state.growmaxfive_mask) EQ (select_mask_value+1)
            ;record the maskindex and countmask
            index = WHERE(temp_mask EQ 1, count)
            IF count GT 1 THEN BEGIN
              temp_maskind_list[i] = PTR_NEW(index)
              temp_countmask_list[i]=count
            ENDIF
          ENDIF ELSE BEGIN
            IF PTR_VALID(state.growmaxfive_mask) THEN BEGIN
              temp_mask = temp_mask OR (*(state.growmaxfive_mask) EQ (select_mask_value+1))
              ;record the maskindex and countmask
              index = WHERE(*(state.growmaxfive_mask) EQ (select_mask_value+1), count)
              IF count GT 1 THEN BEGIN
                temp_maskind_list[i] = PTR_NEW(index)
                temp_countmask_list[i]=count
              ENDIF
            ENDIF
          ENDELSE
        ENDIF
      ENDFOR
    ENDFOR
    
    oPolygon = model_lists[0]->get(pos=0)
    oPolygon->getProperty,POLYGONS=POLYGONS
    oPolygon->getProperty,DATA=DATA
    oPolygon->getProperty,COLOR=COLOR
    model_lists[0]->IDLgrModel::getproperty,parent=parent
    state.selected = parent
    str = "New view selected"
    help, parent, /struc
    ; point the oModelMan at the node...
    state.oModelMan->GetProperty,target=manip
    if (manip ne state.selected) then begin
      state.oModelMan->SetTarget,obj_new()
    end
    if ((state.selected ne state.oCurrentTopModel) and $
      (state.selected ne obj_new())) then begin
      state.oModelMan->SetTarget,state.selected
    end

    widget_control, state.text, set_value=str
    demo_draw, state.win, state.scene, debug=state.debug


  IF N_ELEMENTS(temp_mask) GT 1 THEN BEGIN
    SHADE_VOLUME, temp_mask,0.9, v, p
    smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;    oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedVertices , POLYGON=p)
;    mymodel1 = OBJ_NEW('IDLgrModel')
;    mymodel1->Add, oPoly1
;    xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]
    
    ;define to the rotation center
    smoothedOutverts = smoothedVertices
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]
    sel_one={name:allinfo.name,$
      p_verts:PTR_NEW(p_verts),$
      p_conn:PTR_NEW(p_conn),$
      verts:smoothedOutverts,$
      conn:p, $
      color:COLOR,$
      xcoord_conv:[-xmid, cc[1]], $
      ycoord_conv:[-ymid, cc[1]], $
      zcoord_conv:[-zmid, cc[1]], $
      HU_min_value:HU_min_value, $
      HU_max_value:HU_max_value $
    }
    
    ;define the struc_modifyobject
    IF N_ELEMENTS(struc_modifyobject) EQ 0 THEN BEGIN
      ;; create the struc_modifyobject structure
      cur_patient=0
      ;; create a struc_modifyobject structure
      struc_modifyobject = {struc_modifyobject}
    ENDIF ELSE BEGIN
      n=N_ELEMENTS(struc_modifyobject)
      p=replicate({struc_modifyobject},n+1)
      p[0:n-1]=struc_modifyobject[0:n-1]
      cur_patient=n
      struc_modifyobject=p
    ENDELSE
    ;; initialize the elements for THIS struc_modifyobject
    struc_modifyobject[cur_patient].modify_step = 0
    struc_modifyobject[cur_patient].current_step_index = 0
    struc_modifyobject[cur_patient].exist_num[0] = ncount
    struc_modifyobject[cur_patient].hide_status = 1
    struc_modifyobject[cur_patient].hide_status[0, 0:(ncount-1)] = 0
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    FOR i=0,ncount-1 DO BEGIN
      oPolygon = model_lists[i]->get(pos=0)
      oPolygon->getProperty,POLYGONS=POLYGONS
      oPolygon->getProperty,DATA=DATA
      ;oPolygon->getProperty,COLOR=COLOR
      COLOR = color_rgb[i+1, *]
      struc_modifyobject[cur_patient].p_ori_verts_list[i] = PTR_NEW(DATA)
      struc_modifyobject[cur_patient].p_maskindex_list[i] = temp_maskind_list[i]
      struc_modifyobject[cur_patient].countmask_list[i]=temp_countmask_list[i]
      ;define a new polygon
      vert_colors = DATA
      vert_colors[0,*] = color_rgb[i+1, 0]
      vert_colors[1,*] = color_rgb[i+1, 1]
      vert_colors[2,*] = color_rgb[i+1, 2]
      oPoly = OBJ_NEW('IDLgrPolygon', DATA, POLYGON=POLYGONS,  NAME=STRING(i), STYLE=2, $
        SHADING=1, vert_colors=vert_colors, $
        XCOORD_CONV=sel_one.xcoord_conv, YCOORD_CONV=sel_one.ycoord_conv, ZCOORD_CONV=sel_one.zcoord_conv)
      struc_modifyobject[cur_patient].polygon_list[i] = PTR_NEW(oPoly)
    ENDFOR
    
    d_surfview1,sel_one, *(state.vol_HU_cube), temp_mask, temp_mask, state.vol_HU_cube_resolution, struc_modifyobject, cur_patient, state.icon_dir, additionalDir = *(state.additionalDir)
    ;
    
  ENDIF
    widget_control, ev.top, set_uvalue=state, /no_copy

  end

;;GJ, add to refresh the plot
;  'MYREFRESH':begin
;  widget_control, ev.top, get_uvalue=state, /no_copy
;  if state.selected ne obj_new() then begin
;    oPolygon =state.selected->get(pos=0)
;    oPolygon->setProperty,POLYGONS=Outconn1
;    oPolygon->setProperty,DATA=smoothedOutverts
;    widget_control, ev.top, set_uvalue=state, /no_copy
;  endif else begin
;    widget_control, state.text, set_value="No current selection!"
;  endelse
;
;end

  'MYRESET': begin
        widget_control, ev.top, get_uvalue=state, /no_copy
        ;闁兼儳鍢茶ぐ鍣坥del_lists 闁告瑥锕弳杈ㄦ償閿燂拷       
         state.oCurrentView->GetProperty, uvalue=view_uval
         count = n_elements( $
          *(state.model_lists[view_uval.num]) $
          )
          model_lists=*(state.model_lists[view_uval.num])
        state.oModelMan->SetTarget, obj_new()

        FOR i=0,count-1 DO BEGIN
;          d_objworld2Event,{ $
;          top: ev.top, $
;            handler: 0L, $
;            id: wDel $
;            }
            state.oCurrentTopModel->remove,model_lists[i]
           IF (Obj_Valid(model_lists[i])) then obj_destroy,model_lists[i]
        ENDFOR
        
        tmp1 = obj_new('IDLgrModel',name='Surface') 
        tmp2 = obj_new('IDLgrModel',name='Green Light')
        ptr_free, state.model_lists[view_uval.num]
        state.model_lists[view_uval.num] = ptr_new([tmp1,tmp2, $
          obj_new() $ ; Placeholder NULL at end of each list.
          ])
        

        state.selected = state.oCurrentTopModel
        ; state.selected = tmp_obj->Get(position=1)
     
        state.oModelMan->SetTarget, state.selected
        
        IF PTR_VALID((state.Out)[0]) THEN BEGIN
          myGetInfo_one,state,vol_HU_cube,myindex_min,myindex_max, myRGB,mySTR,cc
          len=size(state.Out,/DIMENSIONS)
          ; RGBtable=[[151,173,172],[104,82,83],[173,137,118],[82,118,137],[255,150,128],[0,105,127],[0,34,40],[255,221,215],[255,94,72],[0,161,183]]
          cc=[0,1]
          if len ne 0 then begin
            ma=max((*(state.Out)[0]).vert)
            mi=min((*(state.Out)[0]).vert)
            bias=(ma-mi)/2
            rl=mi+bias
            rr=ma+bias
            cc=[(-rl)/(rr-rl),1/(rr-rl)]
            xma = max(((*(state.Out)[0]).vert)[0,*])
            xmi = min(((*(state.Out)[0]).vert)[0,*])
            xmid = 0.5*(xma+xmi)*cc[1]
            yma = max(((*(state.Out)[0]).vert)[1,*])
            ymi = min(((*(state.Out)[0]).vert)[1,*])
            ymid = 0.5*(yma+ymi)*cc[1]
            zma = max(((*(state.Out)[0]).vert)[2,*])
            zmi = min(((*(state.Out)[0]).vert)[2,*])
            zmid = 0.5*(zma+zmi)*cc[1]
            FOR j=n_elements((state.Out))-1, 0, -1 DO BEGIN
              oPoly =OBJ_NEW('IDLgrPolygon',COLOR=myRGB(*,j),(*(state.Out)(j)).vert,POLYGON=(*(state.Out)(j)).conn,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
              mymodel=OBJ_NEW('IDLgrModel')
              mymodel->Add, oPoly
              mymodel->SetProperty,name=(*(state.Out)(j)).name
              setSelected,state,myModel
              myModelAdd,state,myModel
            ENDFOR     
          
        
         
          
            ;Make the new object be the current selection...
            ;
            oModel=mymodel
            state.oCurrentView->GetProperty, uvalue=view_uval
            state.selected = oModel
            g = oModel->get(pos=0)
            if (obj_isa(g,'IDLgrText')) then begin
              rect = state.win->gettextdimensions(g)
            end
        
            state.oModelMan->SetTarget, state.selected
            state.selected->GetProperty, name=s
            ; state.selected->GetProperty, uvalue=s
            str = "Current selection:" + s
            widget_control, state.text, set_value=str
            widget_control, state.wModelDeleteButton, sensitive=1
            widget_control, state.wModifyButton, sensitive=1
            widget_control, state.wModifyAllButton, sensitive=0
;            widget_control, state.wAddChildButton, sensitive=1
            widget_control, state.wUnselectButton, sensitive=1
            widget_control, state.wModelModeRadio, sensitive=1
            widget_control, state.wSaveButton, sensitive=([1,0])[lmgr(/demo)]
            widget_control, state.wModelSelectButton, $
              sensitive= $
              n_elements(*(state.model_lists[view_uval.num])) gt 2
            demo_draw, state.win, state.scene, debug=state.debug
          endif
        ENDIF ELSE BEGIN
          myWorldAdd_one,state
          ;GJ, 20180906
;          print, '(PTR_VALID(state.Info_tables)) = ', (PTR_VALID(state.Info_tables))
;          print, 'PTR_VALID(state.vol_HU_cube) = ', PTR_VALID(state.vol_HU_cube)
          IF PTR_VALID(state.vol_HU_cube) THEN BEGIN
            widget_control, state.wSplitFive, sensitive=1
            widget_control, state.wModifyAllButton, sensitive=0
          ENDIF ELSE BEGIN
            widget_control, state.wModifyButton, sensitive=0
            widget_control, state.wModifyAllButton, sensitive=0
          ENDELSE
          
          ;GJ, 20180906 end
        ENDELSE
               
        widget_control, ev.top, set_uvalue=state, /no_copy
   
    end
    
    'SPLITFIVE':begin
      widget_control, ev.top, get_uvalue=state, /no_copy
      
      if ((state.selected ne obj_new()) AND $
        (state.selected ne state.oCurrentTopModel)) then begin
        if (state.cur_tool eq 0) then begin
          state.oModelMan->SetTarget, obj_new()
          state.selected->GetProperty, parent=p
          p->remove, state.selected
          state.selected->getProperty,uvalue=obj_uval
          IF (N_Elements(obj_uval) gt 0 && Obj_Valid(obj_uval[0])) then obj_destroy,obj_uval
          obj_destroy, state.selected
          state.oCurrentView->GetProperty, uvalue=view_uval
          indx = where( $
            obj_valid(*(state.model_lists[view_uval.num])), $
            count $
            )
          myWorldAdd_five,state
          widget_control, state.wSplitFive, sensitive=0
          widget_control, state.wModifyAllButton, sensitive=1
        endif
      endif
    
      widget_control, ev.top, set_uvalue=state, /no_copy
    end

  'MYBACK':begin

end
'SAVE' : begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  if ((state.selected NE obj_new()) $
    and (state.selected NE state.oCurrentTopModel) $
    and (state.cur_tool NE 1)) then begin
    file = dialog_pickfile(/write, $
      file='untitled.sav', $
      group=ev.top, $
      filter='*.sav' $
      )
    if (file NE '') then begin
      ;
      ;               Isolate tmp_obj from the tree.
      ;
      state.selected->GetProperty, parent=parent
      parent->remove, state.selected
      state.oModelMan->SetTarget, obj_new()
      tmp_obj = state.selected
      ;
      ;               Save it.
      ;
      save, tmp_obj, filename=file
      ;
      ;               Repair the tree.
      ;
      parent->add, state.selected
      state.oModelMan->SetTarget, state.selected
    end
  end
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'MODELSELECT': begin ; Select next object.
  widget_control, ev.top, get_uvalue=state, /no_copy
  wDraw = state.wDraw
  widget_control, ev.top, set_uvalue=state, /no_copy
  d_objworld2Event, { $
    id: wDraw, $
    top: ev.top, $
    handler: 0L, $
    type: 0, $; Button down
    press: 1, $ ; Left Mouse-button.
    x: -2, $
    y: -2  $
  }
end
'VIEWSELECT': begin ; Select next view.
  widget_control, ev.top, get_uvalue=state, /no_copy
  wDraw = state.wDraw
  widget_control, ev.top, set_uvalue=state, /no_copy
  d_objworld2Event, { $
    id: wDraw, $
    top: ev.top, $
    handler: 0L, $
    type: 0, $; Button down
    press: 1, $ ; Left Mouse-button.
    x: -2, $
    y: -2  $
  }
end

'UNSELECT': begin
  widget_control, /hourglass
  widget_control, ev.top, get_uvalue=state, /no_copy
  state.selected = state.oCurrentTopModel

  widget_control, state.wModelDeleteButton, sensitive=0
  widget_control, state.wModifyButton, sensitive=0
  widget_control, state.wModifyAllButton, sensitive=1
 ; widget_control, state.wAddChildButton, sensitive=0
  widget_control, state.wUnselectButton, sensitive=0
  widget_control, state.wModelModeRadio, sensitive=0
  widget_control, state.wModelSelectButton, sensitive=1
  widget_control, state.wSaveButton, sensitive=0

  widget_control, state.text, set_value="No current selection"
  state.oModelMan->SetTarget, obj_new()
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end

'TOOL': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  tool = d_objworld2ToggleState(state.wToolButton)
  case 1 of
    (state.cur_tool eq tool): ; Do Nothing...
    tool eq 1: begin ; View Manipulator tool selected.
      state.oModelMan->SetTarget, obj_new()
      state.oViewMan->SetTarget, state.oCurrentView, state.win
      state.selected = state.oCurrentView
      state.selected->GetProperty, uvalue=view_uvalue
      widget_control, state.text, $
        set_value="Current selection:" + view_uvalue.name
      demo_draw, state.win, state.scene, debug=state.debug

      widget_control, state.wModelModeRadio, sensitive=0
      widget_control, state.wViewControlBase,  map=1
      widget_control, state.wModelControlBase, map=0
;      widget_control, state.wLoadButton, sensitive=0
      widget_control, state.wLoadDBButton, sensitive=0
      widget_control, state.wSaveButton, sensitive=0

      state.cur_tool = 1
    end
    tool eq 0: begin ; Model Manipulator tool selected.
      state.oViewMan->SetTarget, obj_new(), state.win

      wDraw = state.wDraw
      state.cur_tool = 0
      widget_control, ev.top, set_uvalue=state, /no_copy
      d_objworld2Event, { $
        id: wDraw, $
        top: ev.top, $
        handler: 0L, $
        type: 0, $; Button down
        press: 1, $ ; Left Mouse-button.
        x: -1, $
        y: -1  $
      }
      widget_control, ev.top, get_uvalue=state, /no_copy
      widget_control, state.wViewControlBase,  map=0
      widget_control, state.wModelControlBase, map=1
      state.oCurrentView->GetProperty, uvalue=view_uval
      num = n_elements( $
        *(state.model_lists[view_uval.num]) $
        )
      widget_control, $
        state.wModelSelectButton, $
        sensitive=([0,1])[num gt 2]
;      widget_control, state.wLoadButton, sensitive=1
         widget_control, state.wLoadDBButton, sensitive=1
      widget_control, $
        state.wSaveButton, $
        sensitive=([1,0])[lmgr(/demo)]
    end
  endcase
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'MODELMODE': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  d_objworld2NewMode, state, ev.value
  print, 'modelmode = ', ev.value 
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'MODELMODE0': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  ;1 translate, 0 rotate, 2 scale
  d_objworld2NewMode, state, 0
  print, 'modelmode = ', 0
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'MODELMODE1': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
    ;1 translate, 0 rotate, 2 scale
  d_objworld2NewMode, state, 1
  print, 'modelmode = ', 1
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'MODELMODE2': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
    ;1 translate, 0 rotate, 2 scale
  d_objworld2NewMode, state, 2
  print, 'modelmode = ', 2
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end

'ADDVIEW': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  state.win->GetProperty,dim=wdim
  state.oCurrentView = d_objworld2MakeView( $
    wdim[0], $
    wdim[1], $
    {name:'ObjView ' + strcompress(state.highest_view_count), $
    num: state.highest_view_count}, state.icon_dir $
    )
  state.model_lists[state.highest_view_count] = ptr_new(obj_new())
  state.oTrackballMB1->Reset, wdim/2.0, wdim[0]/2.0, mouse=4
  state.oTrackballMB2->Reset, wdim/2.0, wdim[0]/2.0, mouse=2
  state.view_count = state.view_count + 1
  state.highest_view_count = state.highest_view_count + 1
  state.scene->add, state.oCurrentView
  state.oViewMan->SetTarget, state.oCurrentView, state.win
  state.selected = state.oCurrentView
  state.oCurrentView->GetProperty, uvalue=view_uvalue
  widget_control, state.wViewDeleteButton, $
    sensitive=([0,1])[view_uvalue.num ne 0]
  widget_control, state.wViewSelectButton, sensitive=1
  d_objworld2GetViewObjs, state.selected, w,b,t
  state.oWorldRotModel = w
  state.oBasePlatePolygon = b
  state.oCurrentTopModel = t
  str = "Current selection:" + view_uvalue.name
  widget_control, state.text, set_value=str
  demo_draw, state.win, state.scene, debug=state.debug
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'ADD': begin
  
; The original code for adding one STL object, canceled by THL,2020.10.9

;  widget_control, /hourglass
;  widget_control, ev.top, get_uvalue=state, /no_copy
;  if state.oBasePlatePolygon ne obj_new() then begin
;    numObj = state.numObj
;    oModel = d_objworld2MakeGroupObj((where(state.addgroup_subjects eq uval[1]))[0], state.theFont, numObj)
;    state.numObj = numObj
;    IF ~OBJ_VALID(oModel) THEN BEGIN
;      widget_control, ev.top, set_uvalue=state, /no_copy
;      return
;    ENDIF
;    d_objworld2Add, state, oModel, isaPolygon
;    ;GJ, 2019/6/5, after loading the object, hiding baseplate with company icons
;    state.obaseplate->SetProperty, hide = 1
;    
;    IF isaPolygon THEN widget_control, state.wModifyButton,  sensitive = 1 ELSE widget_control, state.wModifyButton,  sensitive = 0
;    widget_control, state.wModifyAllButton,  sensitive = 0
;;    widget_control, state.wModelModeRadio, get_value=mode
;    widget_control, state.wRotateRadio, SET_BUTTON = 1
;    state.oModelMan->SetProperty, mode=1
;    ;
;    ;               Determine how many things there thee are to select.
;    ;
;    if state.cur_tool eq 0 then begin
;      state.oCurrentView->GetProperty, uval=view_uval
;      num_selectables = n_elements( $
;        *(state.model_lists[view_uval.num]) $
;        ) - 1 ; Last item on list is obj_new().
;    end $
;    else $
;      num_selectables = state.view_count
;    print, 'number of selectable objects: ', num_selectables-2
;    num_omodels = num_selectables-2
;    arrange_nine = DBLARR(8, 3)
;    d_nine = 0.8
;    arrange_nine[0,*] = [d_nine, -d_nine, 0]
;    arrange_nine[1,*] = [0, -d_nine, 0]
;    arrange_nine[2,*] = [-d_nine, -d_nine, 0]
;    arrange_nine[3,*] = [d_nine, 0, 0]
;    arrange_nine[4,*] = [-d_nine, 0, 0]
;    arrange_nine[5,*] = [d_nine, d_nine, 0]
;    arrange_nine[6,*] = [0, d_nine, 0]
;    arrange_nine[7,*] = [-d_nine, d_nine, 0]
;    IF num_omodels GE 2 THEN BEGIN
;      (*(state.model_lists[view_uval.num]))[1] -> translate,arrange_nine[(num_omodels-2) MOD 8, 0], arrange_nine[(num_omodels-2) MOD 8, 1], arrange_nine[(num_omodels-2) MOD 8, 2]
;    ENDIF
;    demo_draw, state.win, state.scene, debug=state.debug
;  end
;  widget_control, ev.top, set_uvalue=state, /no_copy
  
  ;Adding one STL object,and generate three another copies in different directions automatically
  ;Note that this operation will override the existed STL objects
  ;THL,2020.10.9
  widget_control, /hourglass
  widget_control, ev.top, get_uvalue=state, /no_copy
  if state.oBasePlatePolygon ne obj_new() then begin
    ;GJ, 2019/6/5, after loading the object, hiding baseplate with company icons
    state.obaseplate->SetProperty, hide = 1
    demo_draw, state.win, state.scene, debug=state.debug
  
    if state.cur_tool eq 0 then begin
      state.oCurrentView->GetProperty, uval=view_uval
      num_selectables = n_elements( $
        *(state.model_lists[view_uval.num]) $
        ) - 1 ; Last item on list is obj_new().
      print,"here"
    end $
    else $
      num_selectables = state.view_count
    print, 'number of selectable objects: ', num_selectables-2
    num_omodels = num_selectables-2
  
  
    IF num_omodels gt 3 THEN BEGIN
      for i =0,3  do begin
        (*(state.model_lists[view_uval.num]))[i]->GetProperty, parent=p
        p->remove, (*(state.model_lists[view_uval.num]))[i]
        ;IF obj_valid((*(state.model_lists[view_uval.num]))[i,0]) THEN obj_destroy,(*(state.model_lists[view_uval.num]))[i]
        ;obj_destroy, (*(state.model_lists[view_uval.num]))[i+1]
        state.view_count =state.view_count - 1
      endfor
    ENDIF
    for i =0,3  do begin
      numObj = state.numObj
      oModel = d_objworld2MakeGroupObj((where(state.addgroup_subjects eq uval[1]))[0], state.theFont, numObj)
      state.numObj = numObj
      IF ~OBJ_VALID(oModel) THEN BEGIN
        widget_control, ev.top, set_uvalue=state, /no_copy
        return
      ENDIF
      d_objworld2Add, state, oModel, isaPolygon
    endfor
    translate_op = DBLARR(1, 3)
    translate_op[*] = [0.245 ,0.7 ,0]
    rotate_op = DBLARR(1, 4)
    rotate_op[*] = [-90,-180,-270,0]
    state.oCurrentView->GetProperty, uval=view_uval
    for i =0,3  do begin
      (*(state.model_lists[view_uval.num]))[3-i] ->translate,translate_op[0],translate_op[1],translate_op[2]
      (*(state.model_lists[view_uval.num]))[3-i] ->rotate,[0,0,1],rotate_op[i]
    endfor
    IF isaPolygon THEN widget_control, state.wModifyButton,  sensitive = 1 ELSE widget_control, state.wModifyButton,  sensitive = 0
    widget_control, state.wModifyAllButton,  sensitive = 0
    ;    widget_control, state.wModelModeRadio, get_value=mode
    widget_control, state.wRotateRadio, SET_BUTTON = 1
    state.oModelMan->SetProperty, mode=1
    ;
    ;               Determine how many things there thee are to select.
    ;
    demo_draw, state.win, state.scene, debug=state.debug
  end
  widget_control, ev.top, set_uvalue=state, /no_copy
  
end
'ADDCHILD': begin
  widget_control, /hourglass
  widget_control, ev.top, get_uvalue=state, /no_copy
  if state.oBasePlatePolygon ne obj_new() then begin
    d_objworld2Add, $
      state, $
      d_objworld2MakeObj( $
      (where(state.addable_subjects eq uval[1]))[0], $
      state.theFont $
      ), isaPolygon, $
      /as_child
   IF isaPolygon THEN widget_control, state.wModifyButton,  sensitive = 1 ELSE widget_control, state.wModifyButton,  sensitive = 0
    widget_control, state.wModifyAllButton,  sensitive = 0
;    widget_control, state.wModelModeRadio, get_value=mode
    widget_control, state.wRotateRadio, SET_BUTTON = 1
    state.oModelMan->SetProperty, mode=1
  end
  widget_control, ev.top, set_uvalue=state, /no_copy
end

'DEL': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  
  ;GJ, 2019/6/5, after loading the object, company pictures hide
; state.oCurrentView->SetProperty, hide = 1

  if ((state.selected ne obj_new()) AND $
    (state.selected ne state.oCurrentTopModel)) then begin
    if (state.cur_tool eq 1) then begin
      state.selected->GetProperty, uvalue=uvalue
      if (uvalue.num ne 0) then begin ; cannot delete first one
        state.oViewMan->SetTarget, obj_new(), state.win
        state.selected->GetProperty, parent=p
        p->remove, state.selected
        state.selected->getProperty,uvalue=uval
        IF obj_valid(uval[0]) THEN obj_destroy,uval
        obj_destroy, state.selected
        state.view_count =state.view_count - 1

        widget_control, state.wViewSelectButton, $
          sensitive=([0,1])[state.view_count gt 1]
       
        ;
        ;               Determine how many things there thee are to select.
        ;
        if state.cur_tool eq 0 then begin
          state.oCurrentView->GetProperty, uval=view_uval
          num_selectables = n_elements( $
            *(state.model_lists[view_uval.num]) $
            ) - 1 ; Last item on list is obj_new().
        end $
        else $
          num_selectables = state.view_count
        print, 'number of selectable objects: ', num_selectables-2

        ; select next view.
        wDraw = state.wDraw
        widget_control, ev.top, set_uvalue=state, /no_copy
        d_objworld2Event, { $
          id: wDraw, $
          top: ev.top, $
          handler: 0l, $
          type: 0, $  ; Button down
          press: 1, $ ; Left Mouse-button.
          x: -2, $
          y: -2  $
        }
        widget_control, ev.top, get_uvalue=state, /no_copy
      end $
      else begin
        widget_control, state.text, $
          set_value="Cannot delete initial view"
      end
    end $
    else begin ; Current tool is Model Manipulator
      state.oModelMan->SetTarget, obj_new()
      state.selected->GetProperty, parent=p
      IF OBJ_VALID(p) THEN BEGIN
        p->remove, state.selected
        state.selected->getProperty,uvalue=obj_uval
        IF (N_Elements(obj_uval) gt 0 && Obj_Valid(obj_uval[0])) then obj_destroy,obj_uval
        obj_destroy, state.selected
      ENDIF
      
      state.oCurrentView->GetProperty, uvalue=view_uval
      indx = where( $
        obj_valid(*(state.model_lists[view_uval.num])), $
        count $
        )
      if indx[0] eq -1 then begin
        *(state.model_lists[view_uval.num]) = obj_new()
        state.selected = state.oCurrentTopModel
        str = "No current selection"
        widget_control, state.text, set_value=str
        widget_control, state.wModelDeleteButton, sensitive=0
        widget_control, state.wModifyButton, sensitive=0
        widget_control, state.wModifyAllButton, sensitive=1
;        widget_control, state.wAddChildButton, sensitive=0
        widget_control, state.wUnselectButton, sensitive=0
        widget_control, state.wModelSelectButton, sensitive=0
        widget_control, state.wSaveButton, sensitive=0
        widget_control, state.wModelModeRadio, sensitive=0
        demo_draw, state.win, state.scene, debug=state.debug
      end $
      else begin
        *(state.model_lists[view_uval.num]) = [ $
          (*(state.model_lists[view_uval.num])) $
          [indx], $
          obj_new() $
          ]
        ;
        ;                   Select something.
        ;
        
        ;
        ;               Determine how many things there thee are to select.
        ;
        if state.cur_tool eq 0 then begin
          state.oCurrentView->GetProperty, uval=view_uval
          num_selectables = n_elements( $
            *(state.model_lists[view_uval.num]) $
            ) - 1 ; Last item on list is obj_new().
        end $
        else $
          num_selectables = state.view_count
        print, 'number of selectable objects: ', num_selectables-2
        wDraw = state.wDraw
        widget_control, ev.top, set_uvalue=state, /no_copy
        d_objworld2Event, { $
          id: wDraw, $
          top: ev.top, $
          handler: 0L, $
          type: 0, $  ; Button down
          press: 1, $ ; Left Mouse-button.
          x: -1, $
          y: -1  $
        }
        return
      end
    end
  end
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'DRAGQLOW' : begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  state.dragq = 0
  widget_control, state.wDragQLow,    sensitive=0
  widget_control, state.wDragQMedium, sensitive=1
  widget_control, state.wDragQHigh,   sensitive=1
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'DRAGQMEDIUM' : begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  state.dragq = 1
  widget_control, state.wDragQLow,    sensitive=1
  widget_control, state.wDragQMedium, sensitive=0
  widget_control, state.wDragQHigh,   sensitive=1
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'DRAGQHIGH' : begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  state.dragq = 2
  widget_control, state.wDragQLow,    sensitive=1
  widget_control, state.wDragQMedium, sensitive=1
  widget_control, state.wDragQHigh,   sensitive=0
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'GRID' : begin
  widget_control, /hourglass
  widget_control, ev.top, get_uvalue=state, /no_copy
  if (OBJ_VALID(state.oCurrentView)) then begin
    if (OBJ_VALID(state.oBasePlatePolygon)) then begin
      state.oBasePlatePolygon->SetProperty, $
        hide=1-d_objworld2ToggleState(state.wGridButton)
      demo_draw, state.win, state.scene, debug=state.debug
    end
  end
  widget_control, ev.top, set_uvalue=state, /no_copy
end

;GJ 2019/2/6, draw base event to rotate the whole continuously
'TIMER': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  print, 'rotating'
  FOR i=0, 180 DO BEGIN
    state.oWorldRotModel->rotate, [0,1,0], 2
    demo_draw, state.win, state.scene, debug=state.debug
    wait, 0.01  
  ENDFOR
  widget_control, ev.top, set_uvalue=state, /no_copy
end

'DRAW': begin
  widget_control, ev.top, get_uvalue=state, /no_copy

;  print, 'ev.type = ', ev.type
;  print, 'ev.press = ', ev.press

;if(n_elements(Outconn1) ne 0 and n_elements(smoothedOutverts) ne 0)then begin
;    if state.selected ne obj_new() then begin
;      
;      
;      tempModel =state.selected->get(pos=0)
;      oPolygon =tempModel->get(pos=0)
;      help, oPolygon, /struc
;      oPolygon->getProperty,POLYGONS=Outconn_Old
;      oPolygon->getProperty,DATA=Outverts_Old      
;      
;      IF N_ELEMENTS(Outverts_old) LT N_ELEMENTS(smoothedOutverts) THEN BEGIN
;        
;;        IF PTR_VALID((state.vol_HU_cube)) THEN BEGIN
;;          vol_HU_cube = (*(state.vol_HU_cube) GE threshold_min_value) * (*(state.vol_HU_cube) LE threshold_max_value) * (threshold_min_value+1+1000) -1000 ; * (vol_HU_cube LT threshold_max_value)
;;          SHADE_VOLUME, vol_HU_cube,threshold_min_value, Outverts, Outconn1
;;          smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1);, ITERATIONS=200)   ;smooth the segmented volume; GJ 2017/12/22
;;        ENDIF
;        oPolygon->setProperty, DATA = smoothedOutverts
;        oPolygon->setProperty, POLYGONS = Outconn1
;      ENDIF ELSE BEGIN
;        oPolygon->setProperty, POLYGONS = Outconn1
;        oPolygon->setProperty, DATA = smoothedOutverts
;      ENDELSE
;
;      state.win->draw,state.oCurrentView
;    endif else begin
;      widget_control, state.text, set_value="No current selection!"
;    endelse
;    Outconn1=[]
;    smoothedOutverts=[]
;endif
  
  
  
  
  ; Expose.
  if (ev.type EQ 4) then $
    demo_draw, state.win, state.scene, debug=state.debug
  ; Handle trackball updates.
  if state.oTrackballMB2->Update(ev, transform=qmat) then begin
    state.oWorldRotModel->GetProperty, transform=t
    mt = t # qmat
    state.oWorldRotModel->setproperty,transform=mt
    demo_draw, state.win, state.scene, debug=state.debug
  end
  have_mb1_transform = $
    state.oTrackballMB1->Update(ev, transform=mb1_transform)

  ; Handle other events (selection, etc.) ...
  case ev.type of
    0 : begin ; Button press
      case 1 of
        (ev.press EQ 1) AND (state.cur_tool EQ 1): begin
          widget_control, /hourglass
          if ev.x eq -2 then begin
            picked = state.scene->Get()
          end $
          else begin
            state.oViewMan->SetProperty,hide=1
            picked = state.win->select(state.scene, [ev.x,ev.y])
            state.oViewMan->SetProperty, hide=0
          end
          if obj_valid(picked[0]) then begin
            state.selected = picked[0]
            state.selected->GetProperty,uvalue=view_uvalue
            widget_control, $
              state.wViewDeleteButton, $
              sensitive=([1,0])[view_uvalue.num eq 0]
            str = "Current selection:" + view_uvalue.name

            ; "pop" the view
            state.scene->remove, picked[0]
            state.scene->add, picked[0]
          end $
          else begin
            state.selected = obj_new()
            str = "No current selection"
          end

          ; point the oViewMan at the node...
          state.oViewMan->GetProperty, target=manip
          if (manip ne state.selected) then begin
            state.oViewMan->SetTarget,obj_new(),state.win
          end
          if state.selected ne obj_new() then begin
            state.oCurrentView = state.selected
            state.oCurrentView->GetProperty,dim=dim,loc=loc
            state.oTrackballMB1->Reset, loc + dim/2.0, $
              dim[0]/2.0, $
              mouse=4
            state.oTrackballMB2->Reset, loc + dim/2.0, $
              dim[0]/2.0, $
              mouse=2
            d_objworld2GetViewObjs,state.selected,w,b,t
            state.oWorldRotModel = w
            state.oBasePlatePolygon = b
            state.oCurrentTopModel = t
            state.oViewMan->SetTarget,state.selected,state.win
          end
          widget_control, state.text, set_value=str
          state.win->draw,state.scene
          demo_draw, state.win, state.scene, debug=state.debug
        end
        ev.press EQ 1: begin
          widget_control, /hourglass
          if ev.x lt 0 then begin
            state.oCurrentView->GetProperty, uvalue=view_uval
            if (n_elements(*(state.model_lists[view_uval.num])) gt 1) then begin
              state.model_cycle_pos = state.model_cycle_pos + ([0,1])[abs(ev.x) - 1]
              ; Last item on list is obj_new()
              state.model_cycle_pos = state.model_cycle_pos $
                mod (n_elements(*(state.model_lists[view_uval.num])) - 1)
              picked = (state.model_cycle_pos ge 0) ? $
                ((*(state.model_lists[view_uval.num]))[state.model_cycle_pos])->get() : obj_new()
            endif else begin
              picked = obj_new()
            endelse
          endif else begin
            state.oModelMan->setproperty,hide=1
            picked = state.win->select( $
              state.oCurrentView,[ev.x,ev.y] $
              )
            state.oModelMan->setproperty,hide=0
          endelse
          if obj_valid(picked[0]) then begin
            if (picked[0] EQ state.oBasePlatePolygon) $
              then begin
              state.selected = state.oCurrentTopModel
              str = "No current selection"

              widget_control, state.wModelDeleteButton, $
                sensitive=0
              widget_control, state.wModifyButton, $
                sensitive=0
              widget_control, state.wModifyAllButton, $
                sensitive=0
 ;             widget_control, state.wAddChildButton, $
;                sensitive=0
              widget_control, state.wUnselectButton, $
                sensitive=0
              widget_control, state.wSaveButton, $
                sensitive=0
            end $
            else begin
              ;HACK - for the sphere objs, we use the
              ;IDLgrModel::getproperty directly
              if (obj_isa(picked[0],'IDLgrModel') EQ 1) $
                then begin
                picked[0]->IDLgrModel::getproperty,parent=p
              end $
              else begin
                picked[0]->GetProperty, parent=p
              end

              if (state.selected EQ p) then begin
                state.oModelMan->GetProperty, mode=mode
                d_objworld2NewMode, $
                  state, $
                  (mode + 1) mod 3
                IF (mode + 1) mod 3 EQ 0 THEN WIDGET_CONTROL,  state.wTranslateRadio, /SET_BUTTON
                IF (mode + 1) mod 3 EQ 1 THEN WIDGET_CONTROL,  state.wRotateRadio, /SET_BUTTON
                IF (mode + 1) mod 3 EQ 2 THEN WIDGET_CONTROL,  state.wScaleRadio, /SET_BUTTON
              end

              state.oCurrentView->GetProperty, $
                uvalue=view_uval
              state.model_cycle_pos = $
                where(*(state.model_lists[view_uval.num]) eq p)
              state.selected = p
              state.selected->GetProperty, name=s
              str = "Current selection:" + s
              widget_control, state.wModelDeleteButton, $
                sensitive=1
              widget_control, state.wModifyButton, $
                sensitive=1
              widget_control, state.wModifyAllButton, $
                sensitive=0
   ;           widget_control, state.wAddChildButton, $
    ;            sensitive=1
              widget_control, state.wUnselectButton, $
                sensitive=1
              widget_control, state.wModelModeRadio, $
                sensitive=1
              widget_control, state.wSaveButton, $
                sensitive=([1,0])[lmgr(/demo)]

              state.oCurrentView->GetProperty, $
                uvalue=view_uval
              if n_elements( $
                *(state.model_lists[view_uval.num]) $
                ) le 2 then widget_control, $
                state.wModelSelectButton, $
                sensitive=0
            end
          end $
          else begin
            state.selected = state.oCurrentTopModel
            str = "No current selection"

            widget_control, state.wModelDeleteButton, $
              sensitive=0
            widget_control, state.wModifyButton, sensitive=0
            widget_control, state.wModifyAllButton, sensitive=1
;            widget_control, state.wAddChildButton, sensitive=0
            widget_control, state.wUnselectButton, sensitive=0
            widget_control, state.wModelModeRadio, sensitive=0

            ; try to change the current view...
            if ev.x ge 0 then begin
              state.oViewMan->setproperty,hide=1
              picked = state.win->select(state.scene,[ev.x,ev.y])
              state.oViewMan->setproperty, hide=0
              si = size(picked)
              if (si[0] ne 0) then begin
                if (picked[0] ne state.oCurrentView) then begin
                  state.oCurrentView = picked[0]
                  state.oCurrentView->GetProperty,dim=dim,loc=loc
                  state.oTrackballMB1->Reset, loc + dim/2.0, $
                    dim[0]/2.0, $
                    mouse=4
                  state.oTrackballMB2->Reset, loc + dim/2.0, $
                    dim[0]/2.0, $
                    mouse=2
                  d_objworld2GetViewObjs, $
                    state.oCurrentView, $
                    w,$
                    b,$
                    t
                  state.oWorldRotModel = w
                  state.oBasePlatePolygon = b
                  state.oCurrentTopModel = t
                  state.selected = state.oCurrentTopModel
                  str = "New view selected"
                  widget_control, $
                    state.wModelDeleteButton, $
                    sensitive=0
                  widget_control, $
                    state.wModifyButton, $
                    sensitive=0
                  widget_control, $
                    state.wModifyAllButton, $
                    sensitive=1
       ;           widget_control, $
        ;            state.wAddChildButton, $
          ;          sensitive=0

                  ; pop it
                  state.scene->remove, state.oCurrentView
                  state.scene->add, state.oCurrentView
                end
              end
            end

            state.oCurrentView->GetProperty, uvalue=view_uval
            if n_elements( $
              *(state.model_lists[view_uval.num]) $
              ) gt 1 then widget_control, $
              state.wModelSelectButton, $
              sensitive=1

          end

          ; point the oModelMan at the node...
          state.oModelMan->GetProperty,target=manip
          if (manip ne state.selected) then begin
            state.oModelMan->SetTarget,obj_new()
          end
          if ((state.selected ne state.oCurrentTopModel) and $
            (state.selected ne obj_new())) then begin
            state.oModelMan->SetTarget,state.selected
          end

          widget_control, state.text, set_value=str
          demo_draw, state.win, state.scene, debug=state.debug

        end
        ev.press EQ 2: begin
          state.win->setproperty, QUALITY=state.dragq
          widget_control,state.wDraw,/draw_motion
        end
        (ev.press EQ 4) and (state.cur_tool EQ 1): begin
          if (state.selected ne obj_new()) then begin
            state.oViewMan->MouseDown,[ev.x,ev.y],state.win
            state.btndown = 1b
            state.win->setproperty, QUALITY=state.dragq
            widget_control,state.wDraw,/draw_motion
            demo_draw, $
              state.win, $
              state.scene, $
              debug=state.debug
          end
        end
        ev.press EQ 4: begin
          state.win->SetProperty, QUALITY=state.dragq
          widget_control, state.wDraw, /draw_motion
          state.btndown = 1b
          if ((state.selected ne state.oCurrentTopModel) and $
            (state.selected ne obj_new())) then begin
            state.oModelMan->MouseDown, $
              [ev.x,ev.y], $
              state.win
          end
        end
        else: print, 'Ouch!'
      endcase
    end

    2: begin ; Button motion.
      if state.btndown eq 1b then begin
        case 1 of
          state.cur_tool EQ 1: begin
            state.oViewMan->MouseTrack, [ev.x,ev.y], state.win
            demo_draw, $
              state.win, $
              state.scene, $
              debug=state.debug
            state.oCurrentView->GetProperty,dim=dim,loc=loc
            state.oTrackballMB1->Reset, loc + dim/2.0, $
              dim[0]/2.0, $
              mouse=4
            state.oTrackballMB2->Reset, loc + dim/2.0, $
              dim[0]/2.0, $
              mouse=2
          end
          (state.selected ne state.oCurrentTopModel) and $
            (state.selected ne obj_new()): begin
            state.oModelMan->MouseTrack, [ev.x,ev.y], $
              state.win
            demo_draw, $
              state.win, $
              state.scene, $
              debug=state.debug
          end
          else: begin
            ; Rotate.
            if have_mb1_transform then begin
              state.oWorldRotModel->GetProperty, $
                transform=t
              state.oWorldRotModel->SetProperty, $
                transform=t # mb1_transform
              demo_draw, $
                state.win, $
                state.scene, $
                debug=state.debug
            end
          end
        endcase
      end
    end

    1: begin ; Button release.
      if state.btndown eq 1b then begin
        case 1 of
          state.cur_tool EQ 1: $
            state.oViewMan->MouseUp,[ev.x,ev.y],state.win
          (state.selected ne state.oCurrentTopModel) and $
            (state.selected ne obj_new()): $
            state.oModelMan->MouseUp, [ev.x,ev.y], state.win
          else:
        endcase
      end
      state.btndown = 0b
      state.win->setproperty, QUALITY=2
      demo_draw, state.win, state.scene, debug=state.debug
      widget_control,state.wDraw,draw_motion=0
    end
    else:
  endcase
  widget_control, ev.top, set_uvalue=state, /no_copy
end
'HOTKEY': begin
  widget_control, ev.top, get_uvalue=state, /no_copy
  case strupcase(ev.ch) of
    ' ': begin ; "Next" select.
      ;
      ;               Determine how many things there thee are to select.
      ;
      if state.cur_tool eq 0 then begin
        state.oCurrentView->GetProperty, uval=view_uval
        num_selectables = n_elements( $
          *(state.model_lists[view_uval.num]) $
          ) - 1 ; Last item on list is obj_new().
      end $
      else $
        num_selectables = state.view_count
      ;
      ;               Select something.
      ;
      case 1 of
        num_selectables gt 1: begin
          wDraw = state.wDraw
          widget_control, ev.top, set_uvalue=state, /no_copy
          d_objworld2Event, { $
            id: wDraw, $
            top: ev.top, $
            handler: 0L, $
            type: 0, $ ; Button down
            press: 1, $ ; Left Mouse-button.
            x: -2, $
            y: -2  $
          }
          widget_control, ev.top, get_uvalue=state, /no_copy
        end
        num_selectables eq 1: begin
          if state.cur_tool eq 0 and $
            state.selected eq state.oCurrentTopModel then begin
            wDraw = state.wDraw
            widget_control, ev.top, set_uvalue=state, /no_copy
            d_objworld2Event, { $
              id: wDraw, $
              top: ev.top, $
              handler: 0L, $
              type: 0, $ ; Button down
              press: 1, $ ; Left Mouse-button.
              x: -1, $
              y: -1  $
            }
            widget_control, ev.top, get_uvalue=state, /no_copy
          end
        end
        else:
      endcase
      ;
    end
    'S': begin ; Scale.
      if state.cur_tool eq 0 then begin ; (0=Model Manipulator)
        ;widget_control, ev.top, get_uvalue=state, /no_copy
        ;0 translate, 1 rotate, 2 scale
        WIDGET_CONTROL,  state.wScaleRadio, /SET_BUTTON
        d_objworld2NewMode, state, 2
        print, 'modelmode = ', 2
        demo_draw, state.win, state.scene, debug=state.debug
        widget_control, ev.top, set_uvalue=state, /no_copy
      end
    end
    'R': begin ; Rotate.
      if state.cur_tool eq 0 then begin ; (0=Model Manipulator)
          ;widget_control, ev.top, get_uvalue=state, /no_copy
          ;0 translate, 1 rotate, 2 scale
          WIDGET_CONTROL,  state.wRotateRadio, /SET_BUTTON
          d_objworld2NewMode, state, 1
          print, 'modelmode = ', 1
          demo_draw, state.win, state.scene, debug=state.debug
          widget_control, ev.top, set_uvalue=state, /no_copy
      end
    end
    'T': begin ; Translate.
      if state.cur_tool eq 0 then begin ; (0=Model Manipulator)
          ;widget_control, ev.top, get_uvalue=state, /no_copy
          ;0 translate, 1 rotate, 2 scale
          WIDGET_CONTROL,  state.wTranslateRadio, /SET_BUTTON
          d_objworld2NewMode, state, 0
          print, 'modelmode = ', 0
          demo_draw, state.win, state.scene, debug=state.debug
          widget_control, ev.top, set_uvalue=state, /no_copy
      end
    end
    'U': begin ; Unselect
      if state.cur_tool eq 0 then begin
        wUnselectButton = state.wUnselectButton
        widget_control, ev.top, set_uvalue=state, /no_copy
        d_objworld2Event, { $
          top: ev.top, $
          handler: 0L, $
          id: wUnselectButton $
        }
        widget_control, ev.top, get_uvalue=state, /no_copy
      end
    end
    'D': begin ; Delete
      wDel = state.wModelDeleteButton
      widget_control, ev.top, set_uvalue=state, /no_copy
      d_objworld2Event, { $
        top: ev.top, $
        handler: 0L, $
        id: wDel $
      }
      widget_control, ev.top, get_uvalue=state, /no_copy
    end
    'V': begin ; Toggle Manipulate Views mode
      wToolButton = state.wToolButton
      widget_control, ev.top, set_uvalue=state, /no_copy
      d_objworld2Event, { $
        top: ev.top, $
        handler: 0L, $
        id: wToolButton $
      }
      widget_control, ev.top, get_uvalue=state, /no_copy
    end
    else:
  endcase
  widget_control, ev.top, set_uvalue=state, /no_copy
end
endcase
widget_control, ev.top, get_uvalue=state, /no_copy
if xregistered('demo_tour') eq 0 AND N_ELEMENTS(state) then begin
  widget_control, state.wHotKeyReceptor, /input_focus
end
widget_control, ev.top, set_uvalue=state, /no_copy
end


PRO myimage_read

  afn = 'E:\idlProject\demo\S0004_130952_Recon3\'
  ; afn = DIALOG_PICKFILE(/MUST_EXIST, TITLE="dicom files", /DIRECTORY)
  IF STRLEN(afn) NE 0 THEN BEGIn
    a = FILE_SEARCH(afn, '*.dcm', count=nrfile)
    IF nrfile NE 0 THEN BEGIN
      obj = OBJ_NEW('IDLffDICOMex', a[0], /NO_PIXEL_DATA)
      obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
      OBJ_DESTROY, obj
      vol = DBLARR(rows, cols, nrfile)
      voxelsize = DBLARR(nrfile, 5) * 0.
      FOR i=0L, nrfile-1L DO BEGIN
        obj = OBJ_NEW('IDLffDICOMex', a[i])
        result=QUERY_DICOM(a[i], info)
        IF result NE 0 THEN BEGIN
          obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
          order=1
          vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
          IF i eq 0 then begin  
            IF (obj->QueryValue('0018,0088')) EQ 2 THEN spacing = obj->GetValue('0018,0088') ELSE spacing = 0.
            IF (obj->QueryValue('0018,0050')) EQ 2 THEN sliceTh = obj->GetValue('0018,0050') ELSE sliceTh = 0.
            IF (obj->QueryValue('0028,0030')) EQ 2 THEN pixelSp = obj->GetValue('0028,0030') ELSE pixelSp = 0.
            IF sliceTh GT 0. THEN BEGIN
              voxelsize[i, 0] = pixelSp[0]
              voxelsize[i, 1] = pixelSp[1]
              voxelsize[i, 2] = sliceTh + spacing
              voxelsize[i, 3] = rows
              voxelsize[i, 4] = cols
            ENDIF
            
            ;convert grey pixel to HU
            ;rescale intercept
            IF (obj->QueryValue('0028,1052')) EQ 2 THEN rescale_intercept = obj->GetValue('0028,1052') ELSE rescale_intercept = 0.
            ;rescale slope
            IF (obj->QueryValue('0028,1053')) EQ 2 THEN rescale_slope = obj->GetValue('0028,1053') ELSE rescale_slope = 0.
          endif
          OBJ_DESTROY, obj
          

          ;smooth the image and store to 3D
          vol[*,*,i] = vPixels
        ENDIF
      ENDFOR
      vol_HU = vol * DOUBLE(rescale_slope) + DOUBLE(rescale_intercept)
      vol_HU_cube=CONGRID(vol_HU,  FLOOR(voxelsize[0, 3]*voxelsize[0, 0]/voxelsize[0, 2]), FLOOR(voxelsize[0, 4]*voxelsize[0, 1]/voxelsize[0, 2]), nrfile)



      save,vol_HU_cube,filename='img_info.sav'
      ;xobjview_XDU, mymodel,BACKGROUND = [0, 0, 0]

      ;      oWindow = OBJ_NEW('IDLgrWindow', dimension = [600, 600])
      ;      oView = OBJ_NEW('IDLgrView', viewPlane_Rect = [0,0,200,200])
      ;      oModel = OBJ_NEW('IDLgrModel')
      ;      oView->Add, oModel
      ;      oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], Outverts, POLYGON=Outconn)
      ;      oModel->Add, oPoly1
      ;      oModel->Rotate, [1,0,0], 90
      ;      oWindow->Draw, oView

      ;vol_bone = vol_HU_cube*0.
      ;vol_bone[outverts[0,*], outverts[1,*], outverts[2,*]] =255.
      ;      iimage, vol_bone[*,*,20]
      ;ivolume, vol_bone
      ; return,mymodel
      ;obj_destroy, /all
    ENDIF
  ENDIF
END


PRO setSelected,state,model
  state.selected = model
  ; state.selected = tmp_obj->Get(position=1)
  state.oModelMan->SetTarget, state.selected
  state.selected->getproperty,name=s
  str = "current selection:" + s
  widget_control, state.text, set_value=str
end

;PRO myWorldShow,state,myStructs
;  len=size(Out,/DIMENSIONS)
;  RGBtable=[[151,173,172],[104,82,83],[173,137,118],[82,118,137],[255,150,128],[0,105,127],[0,34,40],[255,221,215],[255,94,72],[0,161,183]]
;  cc=[0,1]
;  if len ne 0 then begin
;    ma=max(myStructs(0).verts)
;    mi=min(myStructs(0).verts)
;    bias=(ma-mi)/2
;    rl=mi+bias
;    rr=ma+bias
;    cc=[(-rl)/(rr-rl),1/(rr-rl)]
;    FOR j=0,4 DO BEGIN
;      oPoly =OBJ_NEW('IDLgrPolygon',COLOR=myRGB(*,j),myStructs(j).verts,POLYGON=myStructs(j).conn,XCOORD_CONV=cc,YCOORD_CONV=cc,ZCOORD_CONV=cc)
;      mymodel=OBJ_NEW('IDLgrModel')
;      mymodel->Add, oPoly
;      mymodel->SetProperty,name=(*Out(j)).name
;      setSelected,state,myModel
;      d_objworld2Add,state,myModel,as_child=0
;    ENDFOR
;    setSelected,state,myModel
;
;  endif
;
;end

;閻犱礁澧介悿鍡橈紣濠婂棗顥涢悹瀣暢婢瑰﹪寮剁捄銊ф惣
PRO myGetInfo_five,state,vol_HU_cube, HU_min_value,HU_max_value, color,str,range
    IF(PTR_VALID(state.vol_HU_cube)) THEN vol_HU_cube=*(state.vol_HU_cube)
;    IF(PTR_VALID(state.vol_HU)) THEN vol_HU=*(state.vol_HU)  
;    IF(~PTR_VALID(state.Info_tables)) THEN BEGIN
      ;replace 5 tissues with 1 bone tissue
      sizeVHc = SIZE(vol_HU_cube)
      mean_image = IMAGE_THRESHOLD(BYTSCL(vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
      image_statistics, vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
      bone_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
      IF STRCMP(state.p_modal, 'CT') THEN bone_value = bone_value * 0.5
;      print, 'patient age', state.patient_age
;      IF (state.patient_age LT 16) AND (state.patient_age GT 0) THEN bone_value = 140 ELSE bone_value = 250
      max_bone = MAX(vol_HU_cube)
      color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
      Info_tables={HU_min_value:[bone_value,bone_value,bone_value,bone_value,bone_value], HU_max_value:[max_bone,max_bone,max_bone,max_bone,max_bone], color:TRANSPOSE(color_rgb[1:5, *]), $
;        color:[[254,206,180],[216,50,46],[230,197,183],[215,220,221],[255,255,255]],$
        str:['threshold 1','threshold 2','threshold 3','threshold 4','threshold 5'],$
        range:[0,1]$
;        Info_tables={HU_min_value:[250],$
;          HU_max_value:[2000],$
;          color:[[255,255,255]],$
;          str:['bone'],$
;          range:[0,1]$
      }
      state.Info_tables=PTR_NEW(Info_tables)
;    ENDIF 
    myAlpha=[0.9,0.5,0.8,0.99,1]  
    Info_tables=*(state.Info_tables)
    HU_min_value=Info_tables.HU_min_value
    HU_max_value=Info_tables.HU_max_value
    color=Info_tables.color
    str=Info_tables.str
    range=Info_tables.range
END

PRO myGetInfo_one,state,vol_HU_cube, HU_min_value,HU_max_value, color,str,range
  IF(PTR_VALID(state.vol_HU_cube)) THEN vol_HU_cube=*(state.vol_HU_cube)
  ;    IF(PTR_VALID(state.vol_HU)) THEN vol_HU=*(state.vol_HU)
  IF(~PTR_VALID(state.Info_tables)) THEN BEGIN
    sizeVHc = SIZE(vol_HU_cube)
    mean_image = IMAGE_THRESHOLD(BYTSCL(vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
    image_statistics, vol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
    thres_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5 ;GJ 2019/5/19;-100.
    IF STRCMP(state.p_modal, 'CT') THEN thres_value = thres_value * 0.5
;    IF (state.patient_age LT 16) AND (state.patient_age GT 0) THEN bone_value = 140 ELSE bone_value = 250
    Info_tables={HU_min_value:[thres_value],$
      HU_max_value:MAX(vol_HU_cube),$
      color:[[255,255,255]],$
      str:['threshold'],$
      range:[0,1]$
    }
    state.Info_tables=PTR_NEW(Info_tables)
  ENDIF
  ;myAlpha=[0.9,0.5,0.8,0.99,1]
  Info_tables=*(state.Info_tables)
  HU_min_value=Info_tables.HU_min_value
  HU_max_value=Info_tables.HU_max_value
  color=Info_tables.color
  str=Info_tables.str
  range=Info_tables.range
END


PRO myModelAdd,state,omodel
  ;
  state.oCurrentTopModel->add, oModel
  state.oCurrentView->GetProperty, uvalue=view_uval
  *(state.model_lists[view_uval.num]) = $
    [oModel, *(state.model_lists[view_uval.num])]
  state.model_cycle_pos = 0 
END




PRO  myWorldAdd_five,state, segmentation=segmentation
 
 
  myGetInfo_five,state,vol_HU_cube,myindex_min,myindex_max, myRGB,mySTR,cc
  if(N_ELEMENTS(vol_HU_cube) EQ 0) THEN RETURN  

  ;濞达綀娉曢弫銈咁嚕閿熺晫绠ョ紒鐘愁殰椤︾敻鎮堕崱妤佺闁绘鎷�
  side = 15
  strucElem = FLTARR(side, side, side) *0. + 1.
 segmentation = 1
  IF KEYWORD_SET(segmentation) THEN BEGIN
 
    Out=ptrArr(n_elements(myindex_min))
    threshold_min_value = myindex_min(0)
    threshold_max_value = myindex_max(0)
    growmaxfive,vol_HU_cube,threshold_min_value,threshold_max_value,b,n_grows
    if(n_grows eq 0) then begin
      print,'not found'
    endif else begin
      state.growmaxfive_mask = PTR_NEW(b)
      num_objs = n_grows
      FOR i=1, num_objs DO BEGIN
        temp_mask = (b EQ i)
        SHADE_VOLUME, temp_mask,0.9, v, Outconn1
        smoothedOutverts = MESH_SMOOTH(v, Outconn1, ITERATIONS=50)
        print, 'Point numbers: ', N_ELEMENTS(smoothedOutverts)
        struct={name:mySTR(i-1),vert:smoothedOutverts,conn:Outconn1}
        Out(i-1)=ptr_new(struct)
      ENDFOR
    endelse
  ENDIF ELSE BEGIN
    Out=ptrArr(1)
    threshold_min_value = myindex_min(0)
    threshold_max_value = myindex_max(0)
    vol_HU_cube_mask = (vol_HU_cube GE threshold_min_value) * (vol_HU_cube LE threshold_max_value)
    state.vol_HU_cube_mask = PTR_NEW(vol_HU_cube_mask)
    state.vol_HU_cube_mask_modify = PTR_NEW(vol_HU_cube_mask*0 + 1)
    
    btm_value = MIN([MIN(vol_HU_cube), -1000])
    vol_HU_cube =  vol_HU_cube_mask * (threshold_min_value+1-btm_value) +btm_value ; * (vol_HU_cube LT threshold_max_value)
    
    ;delete the variable
    vol_HU_cube_mask = !NULL

    vol_HU_cube[0, *, *] = vol_HU_cube[0, *, *]*0. +btm_value
    vol_HU_cube[(size(vol_HU_cube))[1]-1, *, *] = vol_HU_cube[(size(vol_HU_cube))[1]-1, *, *]*0. +btm_value
    vol_HU_cube[*, 0, *] = vol_HU_cube[*, 0, *]*0. +btm_value
    vol_HU_cube[*, (size(vol_HU_cube))[2]-1, *] = vol_HU_cube[*, (size(vol_HU_cube))[2]-1, *]*0. +btm_value
    vol_HU_cube[*, *, 0] = vol_HU_cube[*, *, 0]*0. +btm_value
    vol_HU_cube[*, *, (size(vol_HU_cube))[3]-1] = vol_HU_cube[*, *, (size(vol_HU_cube))[3]-1]*0. +btm_value
    SHADE_VOLUME, vol_HU_cube,threshold_min_value, Outverts, Outconn1
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1);, ITERATIONS=200)   ;smooth the segmented volume; GJ 2017/12/22

    print, 'Point numbers: ', N_ELEMENTS(smoothedOutverts)
    struct={name:mySTR(0),vert:smoothedOutverts,conn:Outconn1};,vert_max:smoothedOutverts_max,conn_max:Outconn_max};
    Out(0)=ptr_new(struct)
    
    ;delete the variables
    vol_HU_cube = !NULL
    Outverts = !NULL
    Outconn1 = !NULL
    smoothedOutverts =!NULL
    
  ENDELSE
   
  len=size(Out,/DIMENSIONS)
 ; RGBtable=[[151,173,172],[104,82,83],[173,137,118],[82,118,137],[255,150,128],[0,105,127],[0,34,40],[255,221,215],[255,94,72],[0,161,183]]
  cc=[0,1]
  if len ne 0 then begin
    ma=max((*Out[0]).vert)
    mi=min((*Out[0]).vert)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max(((*Out[0]).vert)[0,*])
    xmi = min(((*Out[0]).vert)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max(((*Out[0]).vert)[1,*])
    ymi = min(((*Out[0]).vert)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max(((*Out[0]).vert)[2,*])
    zmi = min(((*Out[0]).vert)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]
    FOR j=n_elements(Out)-1, 0, -1 DO BEGIN
      oPoly =OBJ_NEW('IDLgrPolygon',COLOR=myRGB(*,j),(*Out(j)).vert,POLYGON=(*Out(j)).conn,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
      mymodel=OBJ_NEW('IDLgrModel')
      mymodel->Add, oPoly
      mymodel->SetProperty,name=(*Out(j)).name
      myModel->Rotate, [-1,0,0], 90
      myModel->Rotate, [0,0,1], 160
      setSelected,state,myModel
      myModelAdd,state,myModel
    ENDFOR     
  
    ;Make the new object be the current selection...
    ;
    oModel=mymodel
    state.oCurrentView->GetProperty, uvalue=view_uval
    state.selected = oModel
    g = oModel->get(pos=0)
    if (obj_isa(g,'IDLgrText')) then begin
      rect = state.win->gettextdimensions(g)
    end

    state.oModelMan->SetTarget, state.selected
    state.selected->GetProperty, name=s
    ; state.selected->GetProperty, uvalue=s
    str = "Current selection:" + s
    widget_control, state.text, set_value=str
    widget_control, state.wModelDeleteButton, sensitive=1
    widget_control, state.wModifyButton, sensitive=1
    widget_control, state.wModifyAllButton, sensitive=0
 ;   widget_control, state.wAddChildButton, sensitive=1
    widget_control, state.wUnselectButton, sensitive=1
    widget_control, state.wModelModeRadio, sensitive=1
    widget_control, state.wSaveButton, sensitive=([1,0])[lmgr(/demo)]
    widget_control, state.wModelSelectButton, $
      sensitive= $
      n_elements(*(state.model_lists[view_uval.num])) gt 2
    demo_draw, state.win, state.scene, debug=state.debug
  endif

  state.out = Out
;  FOR i=0,n_elements(Out)-1 DO BEGIN 
;    if ptr_valid(Out[i]) then ptr_free,Out
;  ENDFOR
    

END

PRO  myWorldAdd_one,state


  myGetInfo_one,state,vol_HU_cube,myindex_min,myindex_max, myRGB,mySTR,cc
  if(N_ELEMENTS(vol_HU_cube) EQ 0) THEN RETURN

  ;濞达綀娉曢弫銈咁嚕閿熺晫绠ョ紒鐘愁殰椤︾敻鎮堕崱妤佺闁绘鎷�
  side = 15
  strucElem = FLTARR(side, side, side) *0. + 1.

  Out=ptrArr(n_elements(myindex_min))
  for i=0,n_elements(myindex_min)-1 do begin
    threshold_min_value = myindex_min(i)
    threshold_max_value = myindex_max(i)

    ;reduce dimesion
    ;      dims = size(vol_HU_cube)
    ;      vol_HU_cube = CONGRID(vol_HU_cube, dims[1]/2, dims[2]/2, dims[3]/2)
    vol_HU_cube_mask = (vol_HU_cube GE threshold_min_value) * (vol_HU_cube LE threshold_max_value)
    state.vol_HU_cube_mask = PTR_NEW(vol_HU_cube_mask)
    state.vol_HU_cube_mask_modify = PTR_NEW(vol_HU_cube_mask*0 + 1)

    btm_value = MIN([MIN(vol_HU_cube), -1000])
    
    vol_HU_cube =  vol_HU_cube_mask * (threshold_min_value+1-btm_value) +btm_value ; * (vol_HU_cube LT threshold_max_value)

    ;
    vol_HU_cube[0, *, *] = vol_HU_cube[0, *, *]*0. +btm_value
    vol_HU_cube[(size(vol_HU_cube))[1]-1, *, *] = vol_HU_cube[(size(vol_HU_cube))[1]-1, *, *]*0. +btm_value
    vol_HU_cube[*, 0, *] = vol_HU_cube[*, 0, *]*0. +btm_value
    vol_HU_cube[*, (size(vol_HU_cube))[2]-1, *] = vol_HU_cube[*, (size(vol_HU_cube))[2]-1, *]*0. +btm_value
    vol_HU_cube[*, *, 0] = vol_HU_cube[*, *, 0]*0. +btm_value
    vol_HU_cube[*, *, (size(vol_HU_cube))[3]-1] = vol_HU_cube[*, *, (size(vol_HU_cube))[3]-1]*0. +btm_value
    SHADE_VOLUME, vol_HU_cube,threshold_min_value, Outverts, Outconn1
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1);, ITERATIONS=200)   ;smooth the segmented volume; GJ 2017/12/22

    print, 'Point numbers: ', N_ELEMENTS(smoothedOutverts)
    struct={name:mySTR(i),vert:smoothedOutverts,conn:Outconn1};,vert_max:smoothedOutverts_max,conn_max:Outconn_max};
    Out(i)=ptr_new(struct)
  endfor



  len=size(Out,/DIMENSIONS)
  ; RGBtable=[[151,173,172],[104,82,83],[173,137,118],[82,118,137],[255,150,128],[0,105,127],[0,34,40],[255,221,215],[255,94,72],[0,161,183]]
  cc=[0,1]
  if len ne 0 then begin
    ma=max((*Out[0]).vert)
    mi=min((*Out[0]).vert)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max(((*Out[0]).vert)[0,*])
    xmi = min(((*Out[0]).vert)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max(((*Out[0]).vert)[1,*])
    ymi = min(((*Out[0]).vert)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max(((*Out[0]).vert)[2,*])
    zmi = min(((*Out[0]).vert)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]
    FOR j=0,N_ELEMENTS(Out)-1 DO BEGIN
      color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
      oPoly =OBJ_NEW('IDLgrPolygon',COLOR=TRANSPOSE(color_rgb[j+1, *]),(*Out(j)).vert,POLYGON=(*Out(j)).conn,XCOORD_CONV=[-xmid, cc[1]],YCOORD_CONV=[-ymid, cc[1]],ZCOORD_CONV=[-zmid, cc[1]])
      ;      oPoly_max =OBJ_NEW('IDLgrPolygon',COLOR=[0,0,0],(*Out(j)).vert_max,POLYGON=(*Out(j)).conn_max,XCOORD_CONV=cc,YCOORD_CONV=cc,ZCOORD_CONV=cc)
      mymodel=OBJ_NEW('IDLgrModel')
      mymodel->Add, oPoly
      ;      mymodel->Add, oPoly_max
      mymodel->SetProperty,name=(*Out(j)).name
      myModel->Rotate, [-1,0,0], 90
      myModel->Rotate, [0,0,1], 160
      setSelected,state,myModel
      myModelAdd,state,myModel
    ENDFOR




    ;Make the new object be the current selection...
    ;
    oModel=mymodel
    state.oCurrentView->GetProperty, uvalue=view_uval
    state.selected = oModel
    g = oModel->get(pos=0)
    if (obj_isa(g,'IDLgrText')) then begin
      rect = state.win->gettextdimensions(g)
    end

    state.oModelMan->SetTarget, state.selected
    state.selected->GetProperty, name=s
    ; state.selected->GetProperty, uvalue=s
    str = "Current selection:" + s
    widget_control, state.text, set_value=str
    widget_control, state.wModelDeleteButton, sensitive=1
    widget_control, state.wModifyButton, sensitive=1
  ;  widget_control, state.wAddChildButton, sensitive=1
    widget_control, state.wUnselectButton, sensitive=1
    widget_control, state.wModelModeRadio, sensitive=1
    widget_control, state.wSaveButton, sensitive=([1,0])[lmgr(/demo)]
    widget_control, state.wModelSelectButton, $
      sensitive= $
      n_elements(*(state.model_lists[view_uval.num])) gt 2
    demo_draw, state.win, state.scene, debug=state.debug
  endif


  ;闂佹彃锕ラ弬涓盪T
  FOR i=0,n_elements(Out)-1 DO BEGIN
    if ptr_valid(Out[i]) then ptr_free,Out
  ENDFOR


END


PRO myWorldExpell


end



;-----------------------------------------------------------------
;+
; :Description:
;    Describe the procedure.
;
;
;
; :Keywords:
;    record_to_filename
;    group
;    debug
;    apptlb
;
; :Author: Administrator
;-
pro d_objworld21, $
  record_to_filename=record_to_filename, $
  group=group, $               ; IN: (opt) group identifier
  debug=debug, $               ; IN: (opt) debug flag
  apptlb=apptlb                ; OUT: (opt) TLB of this application

;  !PATH = Expand_Path('+.lib\coyoteprograms\coyote\') + ';' + !PATH
;  !PATH = Expand_Path('+C:\Program Files\Exelis\IDL83\lib\') + ';' + !PATH
;  !PATH = Expand_Path('+C:\Users\vit\Desktop\software_20180124c\lib\coyoteprograms\coyote\') + ';' + !PATH
  language_file = './lib/icon_pics/language_file.txt'
  ;whether this file exist or not
  IF FILE_TEST(language_file) THEN read_language_file,language_file,icon_dir ELSE icon_dir = './lib/icon_pics/Chinese/'
  widget_control, /hourglass
  ;
  ;Check the validity of the group identifier.
  ;
  ngroup = n_elements(group)
  if (ngroup ne 0) then begin
    check = widget_info(group, /valid_id)
    if (check ne 1) then begin
      print,'Error: the group identifier is not valid.'
      print,'Returning to the main application.'
      return
    end
    groupbase = group
  end $
  else $
    groupbase = 0l
  ;
  ;Get the screen size.
  ;
  device, get_screen_size = screensize
  ;
  ;Set up dimensions of the drawing (viewing) area.
  ;
  xdim = screensize[0]*0.6
  ydim = xdim*0.8
  ;
  ;Get the current color vectors to restore
  ;when this application is exited.
  ;
  tvlct, savedr, savedg, savedb, /get
  ;
  ;Build color table from color vectors
  ;
  colortable = [[savedr],[savedg],[savedb]]
  ;
  ;Get the data size
  sz = size(u)
  ;
  ;
  req_ysize = screensize[1]*0.9;950L
  ;Create widgets.
  ;
  if (n_elements(group) eq 0) then begin
    IF req_ysize GT screensize[1] THEN BEGIN
      wTopBase = widget_base( $
      /column, $
      title="AIMIS3D: AI-based Medical Image Segmentation for 3D Printing and Naked Eye 3D Visualization", $
      bitmap=icon_dir+'BASE1.bmp', $
      xpad=0, $
      ypad=0, $
      /tlb_kill_request_events, ysize =  req_ysize, y_scroll_size  = screensize[1], $
;      tlb_frame_attr=1, $
      mbar=barbase $
      )
    ENDIF ELSE BEGIN
      wTopBase = widget_base( $
        /column, $
        title="AIMIS3D: AI-based Medical Image Segmentation for 3D Printing and Naked Eye 3D Visualization", $
        bitmap=icon_dir+'BASE1.bmp', $
        xpad=0, $
        ypad=0, $
        /tlb_kill_request_events, ysize =  req_ysize, $
        ;      tlb_frame_attr=1, $
        mbar=barbase $
        )
    ENDELSE
  end $
  else begin
    IF req_ysize GT screensize[1] THEN BEGIN
      wTopBase = widget_base( $
      /column, $
      title="AIMIS3D: AI-based Medical Image Segmentation for 3D Printing and Naked Eye 3D Visualization", $
      bitmap=icon_dir+'BASE1.bmp', $
      xpad=0, $
      ypad=0, $
      /tlb_kill_request_events, $
      group_leader=group, ysize =  req_ysize, y_scroll_size  = screensize[1], $
;      tlb_frame_attr=1, $
      mbar=barbase $
      )
    ENDIF ELSE BEGIN
      wTopBase = widget_base( $
        /column, $
        title="AIMIS3D: AI-based Medical Image Segmentation for 3D Printing and Naked Eye 3D Visualization", $
        bitmap=icon_dir+'BASE1.bmp', $
        xpad=0, $
        ypad=0, $
        /tlb_kill_request_events, $
        group_leader=group, ysize =  req_ysize, $
        ;      tlb_frame_attr=1, $
        mbar=barbase $
        )
    ENDELSE
  end
  apptlb = wTopBase ; Return parameter.
  widget_control, wTopBase, set_uname='d_objworld2:tlb'
  ;
  ;Create the menu bar.
  ;
  wFileButton = widget_button(barbase, value='New', /menu)
;  wLoadButton = widget_button( $
;    wFileButton, $
;    value="New", $
;    uval='LOAD' $
;    )
  wLoadDBButton=widget_button($
   wFileButton,$
   value="Open...",$
   uval='DBLOAD'$
   )
  wSortButton=widget_button($
   wFileButton,$
   value="Sort...",$
   uval='SORT'$
   )
  wCrawlerButton=widget_button($
   wFileButton,$
   value="Crawler...",$
   uval='CRAWLER'$
   )
  wSaveButton = widget_button( $
    wFileButton, $
    value="Save selection", $
    uval='SAVE' $
    )
  wPrintButton = widget_button( $
    wFileButton, $
    value="Print", $
    uval='PRINT' $
    )
  wVRMLButton = widget_button( $
    wFileButton, $
    value="Export VRML", $
    uval='VRML' $
    )
  void = widget_button( $
    wFileButton, $
    value='Quit', $
    /separator, $
    uname='d_objworld2:quit', $
    uvalue='QUIT' $
    )   
  ;
  ;
  ;
  wMPIButton = widget_button(barbase, value='MPI_Simulation', /menu)
  wSinExcitationButton = widget_button( $
    wMPIButton, $
    value="Sine Excitation", $
    uval='SinExcitation' $
    )
  
  ;@GJ, 2023/12/22, cheybshev polynomial 2nd kind   
  wSinExcitationCheybshevButton = widget_button( $
    wMPIButton, $
    value="Focal Modulation", $
    uval='Sin_Excitation_Chebyshev_MPI' $
    )
    
  wPulseExcitationButton = widget_button( $
    wMPIButton, $
    value="Pulsed Excitation", $
    uval='PulsedExcitation' $
    )
  wRelaxButton = widget_button( $
    wMPIButton, $
    value="Relaxation Kernel", $
    uval='RelaxationKernel' $
    )
  
  wAndersonButton = widget_button( $
    wMPIButton, $
    value="Anderson Coil", $
    uval='Anderson' $
    )
    
  wLangevinButton = widget_button( $
    wMPIButton, $
    value="Langevin", $
    uval='Langevin' $
    )

  wSignal_H_Button = widget_button( $
    wMPIButton, $
    value="Signal vs H", $
    uval='Signal_H' $
    )
  wMPI_FBP_Button = widget_button( $
    wMPIButton, $
    value="MPI FBP Simulation", $
    uval='MPI_FBP_Sim' $
    )
  wMPS_Data_Analyze_Button = widget_button( $
    wMPIButton, $
    value="MPS Data Analyze", $
    uval='MPS_Data_Analyze' $
    )
  
  wTwin_DCM_Analyze_Button = widget_button( $
    wMPIButton, $
    value="Twin DCM Analyze", $
    uval='Twin_DCM_Analyze' $
    )
    
  ;
  ;Options menu.
  ;
  wOptionsButton = widget_button(barbase, /menu, value="Options")
  wDragQ = widget_button(wOptionsButton, /menu, value="Drag Quality")
  wDragQLow = widget_button( $
    wDragQ, $
    value='Low', $
    uval='DRAGQLOW' $
    )
  wDragQMedium = widget_button( $
    wDragQ, $
    value='Medium', $
    uval='DRAGQMEDIUM' $
    )
  wDragQHigh = widget_button( $
    wDragQ, $
    value='High', $
    uval='DRAGQHIGH' $
    )
  wGridButton = widget_button( $
    wOptionsButton, $
    value="Show Background Logo (on )", $
    uname='d_objworld2:togglegrid', $
    uval='GRID' $
    )
  void = widget_button(wOptionsButton, value="Anti-Alias" ,uval='AA')
  wToolButton = widget_button( $
    wOptionsButton, $
    value="Manipulate Views (off)", $
    uval='TOOL' $
    )
  wClipboardButton = widget_button( $
    wOptionsButton, $
    value="Copy to Clipboard", $
    uval='CLIPBOARD' $
    )
   
  ;Help menu.
  ;
  wHelpButton = widget_button(barbase, /menu, value="Help")
  wWebsite = widget_button(wHelpButton, value="MPI+ Website...", uval='WEBSITE')
  wAbout = widget_button(wHelpButton, value="About MPI+...", uval='AboutAIMIS3D')
   
   
  if (lmgr(/demo)) then begin
    widget_control, wPrintButton, sensitive=0
    widget_control, wClipboardButton, sensitive=0
    widget_control, wVRMLButton, sensitive=0
    widget_control, wSaveButton, sensitive=0
  end

  wTopRowBase = widget_base(wTopBase,/row,/frame)

  addable_subjects = [ $
    'stl File', $
    'Sphere', $
    'Cube', $
    'Tetrahedron', $
    'Cone', $
    'Green Light', $
    'Surface', $
    'Image', $
    'White Light', $
    '3D Text', $
    'Plot', $
    'Ribbon Plot', $
    'Bar Plot', $
    'Textured Surface' $
    ]

  wGuiBase = widget_base(wTopRowBase, /column )

  wStackerBase = widget_base(wGuiBase, xpad=0, ypad=0)

  wViewControlBase = widget_base(wStackerBase, $
    xpad=0, $
    ypad=0, $
    /column, $
    MAP=0)
  void = widget_button( $
    wViewControlBase, $
    value='Add View', $
    uname='d_objworld2:addview', $
    uvalue='ADDVIEW' $
    )
  wViewDeleteButton = widget_button( $
    wViewControlBase, $
    value="Delete", $
    uname='d_objworld2:viewdelete', $
    uval='DEL' $
    )
  widget_control, wViewDeleteButton, sensitive=0

  wViewSelectButton = widget_button( $
    wViewControlBase, $
    value='Select', $
    uname='d_objworld2:viewselect', $
    uval='VIEWSELECT' $
    )
  widget_control, wViewSelectButton, sensitive=0
  wModelControlBase = widget_base(wStackerBase, $
    xpad=0, $
    ypad=0, $
    /column $
    )
  
  ;GJ 2019/2/7
  language_subjects = [ $
    'Chinese', $
    'English', $
    'Arabic', $
    'German', $
    'Japanese', $
    'Korea', $
    'Russia', $
    'Vietinam' $
    ]

  wLanguageButton = widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'language.bmp', /BITMAP,$
    /menu $
    )
  ;this is for adding different language
  ;GJ 2019/2/7
  ;  for i=0, n_elements(language_subjects)-1 do begin
  for i=0, 4 do begin
    void = widget_button(wLanguageButton, $
      VALUE=icon_dir+language_subjects[i]+'.bmp', /BITMAP,$
      uname='d_objworld2:language' + language_subjects[i], $
      uvalue=['LANGUAGE', language_subjects[i]] $
      )
  endfor
  
    
  wOpenButton=widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'open.bmp', /BITMAP, $
    /menu)

  open_subjects = [ $
    'DBLOAD', $
    'CDU']

  wOpen_listButton = LONARR(N_ELEMENTS(open_subjects))
  for i=0, n_elements(open_subjects)-1 do begin
    wOpen_listButton[i] = widget_button(wOpenButton, $
      VALUE=icon_dir+open_subjects[i]+'.bmp', /BITMAP,$
      uname='d_objworld2:open' + open_subjects[i], $
      uvalue=open_subjects[i])
  endfor

    
  wMyAddFolderButton=widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'addfolder.bmp', /BITMAP,$
    uname='d_objworld2: wAaddFolder',$
    uval='ADDFOLDER' $
    )
    
    
  wSplitFive=widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'splitfive.bmp', /BITMAP, $
    uname='d_objworld2:SPLITFIVE',$
    uval='SPLITFIVE' $
    )

  wAddButton = widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'add.bmp', /BITMAP,$
    /menu $
    )
  ;this is for adding different objects, such as body, plants, animals, terracotta, etc
  ;GJ 2019/2/7
  addgroup_subjects = [ $
    'terracotta', $
    'animal', $
    'plant', $
    'mermaide', $
    'humanbody', $
    'computerparts', $
    'key_ring', $
    'readstl' $
    ]

  wAdd_listButton = LONARR(N_ELEMENTS(addgroup_subjects))
  for i=0, n_elements(addgroup_subjects)-1 do begin
    wAdd_listButton[i] = widget_button(wAddButton, $
      VALUE=icon_dir+addgroup_subjects[i]+'.bmp', /BITMAP,$
      uname='d_objworld2:add' + addgroup_subjects[i], $
      uvalue=['ADD', addgroup_subjects[i]] $
      )
  endfor
;  for i=0,n_elements(addable_subjects)-1 do begin
;    void = widget_button(wAddButton, $
;      value=addable_subjects[i], $
;      uname='d_objworld2:add' + addable_subjects[i], $
;      uvalue=['ADD', addable_subjects[i]] $
;      )
;  end
;  wAddChildButton = widget_button( $
;    wModelControlBase, $
;    VALUE=icon_dir+'addchild.bmp', /BITMAP, $
;    /menu $
;    )
;  wAddChild_listButton = LONARR(N_ELEMENTS(addable_subjects))
;  for i=0,n_elements(addable_subjects)-1 do begin
;    wAddChild_listButton[i] = widget_button(wAddChildButton, $
;      value=addable_subjects[i], $
;      uname='d_objworld2:addchild' + addable_subjects[i], $
;      uvalue=['ADDCHILD', addable_subjects[i]] $
;      )
;  end
  wModelDeleteButton = widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'delete.bmp', /BITMAP, $
    uname='d_objworld2:delete', $
    uval='DEL' $
    )
  wModelSelectButton = widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'select.bmp', /BITMAP, $
    uname='d_objworld2:modelselect', $
    uvalue='MODELSELECT' $
    )
  wUnselectButton = widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'unselect.bmp', /BITMAP,$
    uname='d_objworld2:unselect', $
    uvalue='UNSELECT' $
    )
  wModifyButton=widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'modifyselected.bmp', /BITMAP,$
    uname='d_objworld2: wModify',$
    uval='MODIFY' $
    )
  
  wModifyAllButton=widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'modifyall.bmp', /BITMAP,$
    uname='d_objworld2:wModifyAll',$
    uval='MODIFYALL' $
    )
 
;  widget_control, wAddChildButton, sensitive = 0
  widget_control, wModifyButton,  sensitive = 0
  widget_control, wModifyAllButton,  sensitive = 0
  ;GJ, add refresh button
;  wMyRefreshButton=widget_button( $
;    wModelControlBase, $
;    value='Refresh',$
;    uname='d_objworld2: wMyRefresh',$
;    uval='MYREFRESH' $
;    ) 
;    
    
  wMyResetButton=widget_button( $
    wModelControlBase, $
    VALUE=icon_dir+'reset.bmp', /BITMAP,$
    uname='d_objworld2: wMyReset',$
    uval='MYRESET' $
    )

  ;  wModelModeRadio = cw_bgroup( $
  ;    wModelControlBase, $
  ;    ['Translate', 'Rotate', 'Scale'], $
  ;    /exclusive, $
  ;    /no_release, $
  ;    set_value=1, $
  ;    uvalue='MODELMODE' $
  ;    )
  ;
  ;  widget_control, wModelModeRadio, $
  ;    set_uname='d_objworld2:modelmoderadio', sensitive = 0


  _tools = STRUPCASE(['Rotate', 'Translate', 'Scale'])
  if N_ELEMENTS(_tools) gt 1 then begin

    ; Create the exclusive button toolbar.
    wModelModeRadio = WIDGET_BASE(wModelControlBase, /ALIGN_CENTER, /COLUMN, $
      /EXCLUSIVE, SPACE=0, /TOOLBAR)

    for i=0,N_ELEMENTS(_tools)-1 do begin
      case STRUPCASE(_tools[i]) of
        'ROTATE': begin
          wRotateRadio = WIDGET_BUTTON( $
            wModelModeRadio, $
            VALUE=FILEPATH('rotate.bmp', ROOT_DIR=icon_dir), $
            XSIZE =  120, /BITMAP, $
            TOOLTIP='Rotate', $
            UNAME='d_objworld2:rotatemode', $
            UVALUE='MODELMODE1' $
            )
        end
        'TRANSLATE': begin
          wTranslateRadio = WIDGET_BUTTON( $
            wModelModeRadio, $
            VALUE=FILEPATH('translate.bmp', ROOT_DIR=icon_dir), $
            XSIZE =  120, /BITMAP, $
            ;                  bitmap='./lib/icon_pics/Circularmaskpoly.bmp', $
            TOOLTIP='Translate', $
            UNAME='d_objworld2:translatemode', $
            UVALUE='MODELMODE0' $
            )
        end
        'SCALE': begin
          wScaleRadio = WIDGET_BUTTON( $
            wModelModeRadio, $
            VALUE=FILEPATH('scale.bmp', ROOT_DIR=icon_dir), $
            XSIZE =  120, /BITMAP, $
            TOOLTIP='Scale', $
            UNAME='d_objworld2:scalemode', $
            UVALUE='MODELMODE2' $
            )
        end
      endcase
    endfor


    ;  WIDGET_CONTROL,  WIDGET_INFO(wModelModeRadio, /CHILD), /SET_BUTTON
    WIDGET_CONTROL,  wRotateRadio, /SET_BUTTON

    widget_control, wModelModeRadio, $
      set_uname='d_objworld2:modelmoderadio', sensitive = 0


  endif




  ;GJ 2018/11/2, dicom
  wDCMbase = WIDGET_BASE(wModelControlBase, XSIZE=200, YSIZE=150, FRAME=1, $
    /COLUMN)

  wDCMlabel = WIDGET_LABEL(wDCMbase, value='DICOM Crawler')

  wDCMsort = WIDGET_BUTTON(wDCMbase, $
    value=icon_dir+'dcmsort.bmp', /bitmap, $
    UVALUE='SORT', $
    UNAME='d_objworld2:DCMsort')

  wDCMsum = WIDGET_BUTTON(wDCMbase, $
    value=icon_dir+'dcmsum.bmp', /bitmap, $
    UVALUE='DCMSUM', $
    UNAME='d_objworld2:DCMsum')

  wDCMone = WIDGET_BUTTON(wDCMbase, $
    value=icon_dir+'dcmone.bmp', /bitmap, $
    UVALUE='DCMONE', $
    UNAME='d_objworld2:DCMone')

  wDCMpict = WIDGET_BUTTON(wDCMbase, $
    value=icon_dir+'dcmpict.bmp', /bitmap, $
    UVALUE='DCMPICT', $
    UNAME='d_objworld2:DCMpict')
  ;GJ 2018/11/2, dicom

  ;THL, 2019/11/12,Media control
  wMediaCollectbase = WIDGET_BASE(wModelControlBase, XSIZE=200, YSIZE=80, FRAME=1, $
    /COLUMN)
  ;THL, 2019/11/12,gesture control
  wMediaCollectorlabel = WIDGET_LABEL(wMediaCollectbase, value='Media Control')
  wGestureControlButton = WIDGET_BUTTON(wMediaCollectbase, $
    value=icon_dir+'gesture_control.bmp', /bitmap, $
    UVALUE='GESTURECONTROL', $
    UNAME='d_objworld2:Gesturecontrol')

  ;THL, 2020/06/16,voice control
  wVoiceControlButton = WIDGET_BUTTON(wMediaCollectbase, $
    value=icon_dir+'voice_control.bmp', /bitmap, $
    UVALUE='VOICECONTROL', $
    UNAME='d_objworld2:Voicecontrol')
  ;THL, 2019/11/12,Media control

  wStackerBase = widget_base(wTopRowBase, xpad=0, ypad=0, uval='TIMER', uname='d_objworld2:timer')
  wDraw = widget_draw(wStackerBase, $
    xsize=xdim, $
    ysize=ydim, $
    /button_ev, $
    uval='DRAW', $
    retain=0, $
    /expose_ev, $
    uname='d_objworld2:draw', $
    graphics_level=2 $
    )
  wHotKeyReceptor = widget_text(wStackerBase, $
    /all_events, $
    uvalue='HOTKEY', $
    uname='d_objworld2:hotkey' $
    )
  ;
  ;Status readout widget.
  ;
  wGuiBase2 = widget_base(wTopBase, /col)
  wText = widget_text(wGuiBase2,value="object selection",XSIZE = 170)
  wText1 = widget_text(wGuiBase2,value="folder selection",XSIZE = 170, YSIZE=6, /SCROLL);,/dynamic_resize)
  ;
  ;Realize the base widget.
  ;
  widget_control, wTopBase, /realize
  ;
  ;Add demo tips widgets.
  ;
;  sText = demo_gettips( $
;    demo_filepath( $
;    'objworld2.tip', $
;    subdir=['examples','demo', 'demotext'] $
;    ), $
;    wTopBase, $
;    widget_base(wTopBase, map=0, /row) $
;    )
  ;
  ;Get the window id of the drawable.
  ;
  widget_control, wdraw, get_value=win
  ;
  ;Build the scene.
  ;
  scene=obj_new('IDLgrScene')
  oCurrentView = d_objworld2makeview(xdim, ydim, {name:'ObjView', num:0}, icon_dir, obaseplate)
  oCurrentView->getproperty, dim=dim, loc=loc
  scene->add, oCurrentView
  d_objworld2getviewobjs, $
    oCurrentView, $
    oWorldRotModel, $
    oBasePlatePolygon, $
    oCurrentTopModel
  ;
  ;Make a font for the demo.
  ;
  thefont = objarr(4)
  thefont[0] = obj_new('IDLgrFont','times',size=30)
  thefont[1] = obj_new('IDLgrFont','hershey*3',size=9)
  thefont[2] = obj_new('IDLgrFont','helvetica',size=40)
  thefont[3] = obj_new('IDLgrFont','helvetica',size=12)
  ;
  if n_elements(record_to_filename) eq 0 then $
    record_to_filename = ''
  ;
  ;Save state.
  ;
  ;GJ
  ;vol_HU_cube=[!NULL]
  
  
  state = { $
    icon_dir:icon_dir, $
    oTrackballMB1: obj_new('trackball', $
    (loc + dim/2.0), $
    dim[0] / 2.0, $
    mouse=4b $
    ), $
    oTrackballMB2: obj_new('trackball', $
    (loc + dim/2.0), $
    dim[0] / 2.0, $
    mouse=2b $
    ), $
    btndown: 0b,              $
    thefont: thefont,         $
    pt0: fltarr(3),           $
    pt1: fltarr(3),           $
    wDraw: wDraw,             $
    wToolButton: wToolButton,     $
    oWorldRotModel: oWorldRotModel,     $
    oBasePlatePolygon: oBasePlatePolygon, $
    oCurrentView: oCurrentView,       $
    obaseplate: obaseplate,       $;GJ 2019; for hiding the baseplate after loading object
    oModelMan : obj_new('IDLexModelManip', $
    translate=[1,1,1], $
    selector_color=[255,255,255], $
    manipulator_color=[255, 60, 60] $
    ), $
    oViewMan : obj_new('IDLexViewManip', $
    color=[255, 0, 0] $
    ), $
    addgroup_subjects: addgroup_subjects, $
    addable_subjects: addable_subjects, $
    language_subjects: language_subjects, $; GJ 2019/2/7 for international language
    text: wtext,              $
    text1: wtext1,            $;GJ 2020/5/19 for file folder selection
    win: win,                 $
    oCurrentTopModel: oCurrentTopModel, $
    cur_tool: 0,              $
    selected: oCurrentTopModel, $
    scene: scene,             $
    view_count: 1,            $
    highest_view_count: 1, $ ; "High water mark" for view count.
    dragq: 1,                 $
    groupbase: groupbase,     $
    model_lists: ptrarr(50), $ ; One list for each view.
    Out: ptrarr(5), $  ;GJ, 2018/6/30 one list for each polygon
    numObj: 0, $    ;GJ 2019/2/26 for labeling number of object
    wDCMsort: wDCMsort, $ ;GJ, 2019/2/7
    wDCMsum: wDCMsum, $
    wDCMone: wDCMone, $
    wDCMpict: wDCMpict, $; GJ, 2019/2/7
    wGestureControlButton: wGestureControlButton,$  ; THL,2019/11/12
    wVoiceControlButton: wVoiceControlButton,$  ; THL,2020/06/16
    wViewControlBase: wViewControlBase, $
    wViewSelectButton: wViewSelectButton, $
    wModelControlBase: wModelControlBase, $
    wRotateRadio: wRotateRadio, $;GJ, 2019/1/1 
    wTranslateRadio: wTranslateRadio, $;GJ, 2019/1/1
    wScaleRadio: wScaleRadio, $;GJ, 2019/1/1
    wViewDeleteButton: wViewDeleteButton, $
    wModelDeleteButton: wModelDeleteButton, $
    wAddButton: wAddButton, $;GJ 2019/2/7
    wAdd_listButton: wAdd_listButton, $;GJ 2019/2/7
 ;   wAddChildButton: wAddChildButton, $
  ;  wAddChild_listButton: wAddChild_listButton, $; GJ 2019/2/7
    wModelSelectButton: wModelSelectButton, $
    wModelModeRadio: wModelModeRadio, $
    wUnselectButton: wUnselectButton, $
    wOpenButton: wOpenButton,$
    open_subjects: open_subjects, $;GJ 2019/3/23
    wOpen_listButton: wOpen_listButton, $;GJ 2019/2/13
    wSplitFive: wSplitFive, $
    wModifyButton: wModifyButton,$
    wModifyAllButton: wModifyAllButton, $
   ; wMyRefreshButton: wMyRefreshButton, $
    wDragQLow: wDragQLow, $
    wDragQMedium: wDragQMedium, $
    wDragQHigh: wDragQHigh, $
    wGridButton: wGridButton, $
    wHotKeyReceptor: wHotKeyReceptor, $
;    wLoadButton: wLoadButton, $
    wLoadDBButton:wLoadDBButton,$
    wMyAddFolderButton:wMyAddFolderButton,$
    wSortButton:wSortButton,$ ;2018/6/3, sorting dicom images
    wCrawlerButton:wCrawlerButton,$ ;2018/10/27, crawling dicom images
    ;wDBdialog:wDBdialog,$
    wTopBase:wTopBase,$
   ; wDBdialogButtonOk:wDBdialogButtonOk,$
    ;wwDBdialogButtonCancel:wDBdialogButtonCancel,$
    wSaveButton:wSaveButton, $
    wMyResetButton:wMyResetButton,$
    model_cycle_pos: 1, $
    record_to_filename: record_to_filename, $
    debug: keyword_set(debug), $
    colortable: colortable,    $
    vol_HU_cube:PTR_NEW(), $
    vol_HU_cube_mask:PTR_NEW(), $           ;GJ, 2018/5/18
    vol_HU_cube_mask_modify:PTR_NEW(), $           ;GJ, 2018/5/26 ;
    growmaxfive_mask:PTR_NEW(), $                 ;GJ, 2018/6/20
    vol_HU_cube_resolution:0., $
    patient_age:40., $        ;GJ, 2018/8/3
    additionalDir:PTR_NEW(), $  ;GJ 2018/11/12
    p_modal:'', $   ;GJ 2019/05/22
    Info_tables:PTR_NEW(), $
    database:PTR_NEW(), $
    vol_HU_sample:PTR_NEW() $
  }
;  INITDatabase,state
  ;
  ;Restore a sample scene containing a surface and a light
  ;
  ;restore, $
  ;   demo_filepath('objw_surf.sav', $
  ;        subdir=['examples','demo','demodata'] $
  ;        ), $
  ;    /relaxed_structure_assignment
  ;; set up and add names to new objects
  tmp1 = obj_new('IDLgrModel',name='Surface')
  ;surf = tmp_obj->Get(position=0)
  ;tmp_obj->remove,surf
  ;tmp1->add,surf
  tmp2 = obj_new('IDLgrModel',name='None')
  ;gl = tmp_obj->Get(position=0)
  ;tmp_obj->remove,gl
  ;tmp2->add,gl
  ;obj_destroy,tmp_obj
  ;myimage_read
  ;restore,'img_info.sav'

;123

  ;tmp1->translate, 0, 0, .001, /premultiply ; Lift off of baseplate.
  ;tmp2->translate, 0, 0, .001, /premultiply ; Lift off of baseplate.
  ; tmp_obj->translate, 0, 0, .001, /premultiply ; Lift off of baseplate.
  ;
  ;Add tmp_obj to the current tree.
  ;
  ;state.selected->add, tmp1
  ;state.selected->add, tmp2
  ; state.selected->add, tmp_obj
  ;
  ;Add our restored objects to array of selectable objects.
  ;
  state.model_lists[0] = ptr_new([tmp1,tmp2, $
    obj_new() $ ; Placeholder NULL at end of each list.
    ])
  ; state.model_lists[0] = ptr_new([ $
  ;     tmp_obj, $
  ;     tmp_obj->get(position=1), $
  ;     obj_new() $ ; Placeholder NULL at end of each list.
  ;     ])
  ;
  ;Target the Green Light.
  ;
  state.selected = tmp2
  ; state.selected = tmp_obj->Get(position=1)
  state.oModelMan->SetTarget, state.selected
  state.selected->getproperty,name=s
  str = "current selection:" + s
  widget_control, state.text, set_value=str
  ;
  ;Add some rotation on the whole thing.
  ;
;  state.oWorldRotModel->rotate, [-1,0,0], 90
;  state.oWorldRotModel->rotate, [0,1,0], 20
  state.oWorldRotModel->rotate, [0,0,1], -160
  ;



;  restore,'E:\idlProject\demo\img_info.sav'
;  state.vol_HU_cube= ptr_new(vol_HU_cube)
;  myWorldAdd,state



  widget_control, wSplitFive, sensitive=0
  widget_control, wDragQMedium, sensitive=0
  widget_control, wMyAddFolderButton, sensitive=0
  widget_control, wTopBase, set_uvalue=state, /no_copy



  xmanager, 'd_objworld2', wTopBase, $
    event_handler='d_objworld2event', $
    /no_block, $
    cleanup='d_objworld2cleanup'

end

; docformat = 'rst'

;+
; Open an url in the default web browser. On Windows and Mac this is easy. On
; UNIX platforms, the first time this routine is called it will ask for the
; location of your preferred web browser and save this location in
; APP_USER_DIR.
;
; :Params:
;    url : in, required, type=string
;       url to goto in the default web browser
;
; :Requires:
;    IDL 6.1
;
; :Author:
;    Michael Galloy, 2006
;-
pro mg_open_url, url
  compile_opt strictarr

  ; launch the default web browser with the url, unfortunately, this is
  ; platform dependent
  case !version.os_family of
    'Windows' : spawn, 'start ' + url, /hide, /nowait
    else : begin
      ; Mac OS X has a nice way of doing this...
      if (!version.os_name eq 'Mac OS X') then begin
        spawn, 'open ' + url
        return
      endif

      ; ...but the other UNIX platforms don't
      app_readme_text = $
        ['This is the configuration directory for MG_OPEN_URL ', $
        'routine. It is used to save the location of the default ', $
        'web browser between MG_OPEN_URL invocations on UNIX ', $
        'platforms.', $
        '', $
        'It is safe to remove this directory, as it', $
        'will be recreated on demand. Note that all', $
        'settings (e.g. smoke injection depth, juicitron', $
        'angle, etc.) will revert to their default settings.']

      prefdir = app_user_dir('mg', $
        'Michael Galloy', $
        'default-browser', $
        'Default browser location', $
        app_readme_text, 1)
      preffile = filepath('default-browser', root=prefdir)

      if (file_test(preffile)) then begin
        openr, lun, preffile, /get_lun
        browser = ''
        readf, lun, browser
        free_lun, lun
        spawn, browser + ' ' + url
      endif else begin
        browser_location = dialog_pickfile()
        openw, lun, preffile, /get_lun
        printf, lun, browser_location
        free_lun, lun
        msg = ['Your browser location has been stored in:', '', $
          '    ' + preffile, '']
        ok = dialog_message(msg, /info)
        spawn, browser_location + ' ' + url
      endelse
    end
  endcase
end



;$Id: //depot/Release/ENVI51_IDL83/idl/idldir/lib/utilities/xobjview_XDU.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  FILE:
;       xobjview_XDU.pro
;
;  CALLING SEQUENCE:
;       xobjview_XDU, oObj
;
;  INPUTS:
;       Argument oObj is a reference to an atomic graphics object, a
;       reference to an IDLgrModel, or an array of such references.
;
;  INPUT KEYWORDS:
;       XSIZE           Pixels.  Size of drawable.
;       YSIZE           Pixels.  Size of drawable.
;       STATIONARY      Non-moveable objects (e.g. lights) to be viewed.
;       GROUP           Group leader widget.
;       BLOCK           If set, block IDL command line.
;       JUST_REG        If set, register and return (XMANAGER)
;       MODAL           If set, block other IDL widgets from receiving
;                           events.
;       SCALE           Size factor for initial view.  Default: 1/sqrt(3).
;       TITLE           String to appear in xobjview_XDU's title bar.
;       TEST            If set, don't require oObj arg.  Draw a blue
;                           sinusoidal surface instead.
;       DOUBLE_VIEW     Set this keyword for high precision (i.e.
;                           "double precision") graphics.
;       BACKGROUND      Three-element [r,g,b] color vector.
;       REFRESH         Set this keyword to a Top-level base of an xobjview_XDU
;                           to force a redraw of that xobjview_XDU's graphics.
;       XOFFSET         Pixels.  Horzontal postion of xobjview_XDU on the screen.
;       YOFFSET         Pixels.  Vertical postion of xobjview_XDU on the screen.
;       USE_INSTANCING  Set this keyword to enable "instancing."  xobjview_XDU
;                           will use IDL's CREATE_INSTANCE functionality when
;                           drawing graphics, and will use DRAW_INSTANCE if
;                           and when the current graphic is hidden and then
;                           exposed by the user.  See IDLLgrWindow::Draw
;                           documentation for definitions of CREATE_INSTANCE
;                           and DRAW_INSTANCE.
;       DEBUG           Set this keyword to disable error catching and
;                           recovery.  With this keyword set, execution
;                           of xobjview_XDU and/or its supporting routines
;                           will stop where and if an error occurs (unless
;                           other routines with error handling are overriding
;                           this behavior).  This stop-at-error behavior
;                           is often helpful for troubleshooting.
;
;  OUTPUT KEYWORD:
;       TLB             Top-level base widget.
;
;  PURPOSE:
;       Provide a quick way to see graphics objects on the screen.
;
;  CATEGORY:
;       Object Graphics.
;
;  REFERENCE: IDL Reference Guide, IDL User's Guide
;
;  A note about calling xobjview_XDU with the MODAL keyword set:
;       To be modal, xobjview_XDU does not require that its caller specify
;       a group leader.  This is unlike other IDL widget procedures
;       such as XLOADCT, which, to be modal, do require that their
;       caller specify a group leader.  In those other procedures,
;       the requirement exists to encourage the caller to create
;       a modal widget that will be well-behaved with respect to
;       Layering and Iconizing.  (See WIDGET_BASE in the IDL Reference
;       Guide for explanations of "layering" and "iconizing".)   For
;       the same reason the requirement exists in those procedures,
;       supplying an appropriate group leader when invoking
;       xobjview_XDU, /MODAL is good programming practice.  Sometimes,
;       however, it is desirable to invoke xobjview_XDU, /MODAL in a
;       program that uses no other IDL widgets. For this reason,
;       XOBJVIEW allows itself to be invoked MODAL with no group leader.
;
;
;  EXAMPLE 1:
;       IDL> oObj = obj_new('IDLgrSurface', dist(30))
;       IDL> xobjview_XDU, oObj, /bloc
;       IDL> obj_destroy, oObj
;
;  EXAMPLE 2:
;       oObj = obj_new('IDLgrSurface', dist(30))
;       xobjview_XDU, oObj, tlb=tlb
;       ...
;       oObj->SetProperty, color=[255, 0, 0]
;       xobjview_XDU, refresh=tlb
;       ...
;       obj_destroy, oObj
;
;  MODIFICATION HISTORY:
;       9/1999  PCS   - created.
;       11/1999 PCS   - added Update method and related changes.
;       01/2000 PCS   - added TITLE keyword and functionality.
;       07/2000 DMS   - added REFRESH, TLB and BACKGROUND keywords.
;-
;
;--------------------------------------------------------------------
;--------------------------------------------------------------------
pro xobjview_XDU_cleanup, wID
  compile_opt hidden

  widget_control, wID, get_uvalue=pState

  if (*pState).group_leader_is_fabricated then begin
    if widget_info(*(*pState).pGroupLeader, /valid_id) then begin
      widget_control, *(*pState).pGroupLeader, /destroy
    endif
  endif

  obj_destroy, (*pState).oObjviewWid
  obj_destroy, (*pState).oTestSurface
  ptr_free, (*pState).pGroupLeader
  ptr_free, pState
end
;--------------------------------------------------------------------
pro xobjview_XDU_event, event
  compile_opt hidden
  ;
  ;Handle resize events.
  ;
  if tag_names(event, /structure_name) eq 'WIDGET_BASE' then begin
    widget_control, event.top, get_uvalue=pState
    widget_control, /hourglass
    pad = 4 ; Estimate.
    xsize = event.x - pad
    ysize = event.y - pad
    (*pState).oObjviewWid->SetProperty, xsize=MEAN([xsize,ysize]), ysize=MEAN([xsize,ysize])
  endif
end
;--------------------------------------------------------------------
pro xobjview_XDU, $
  oObj, $           ; IN
  xsize=xsize, $    ; IN: (opt) pixels.  Size of drawable.
  ysize=ysize, $    ; IN: (opt) pixels.  Size of drawable.
  stationary=oStationary, $
  group=group_leader, $   ; IN: (opt) group leader widget
  block=block, $
  just_reg=just_reg, $
  modal=modal, $
  scale=scale, $    ; IN: (opt) Scale size of initial view.
  test=test, $      ; IN: (opt) Bool.  Don't require oObj arg.
  title=title, $
  xoffset=xoffset, $; IN: (opt)
  yoffset=yoffset, $; IN: (opt)
  tlb=tlb, $        ;OUT: (opt), widget id of top level base
  refresh=refresh_tlb, $  ; IN: (opt) top-level base of existing xobjview_XDU.
  double_view=double_view, $  ; IN: (opt) set for high precision view.
  background=background, $    ; IN: (opt) view background color.
  renderer=renderer, $        ; IN: (opt) 1 = IDL's software renderer.
  use_instancing=use_instancing, $ ; IN: (opt) For expose events.
  debug=debug       ; IN: (opt)

  compile_opt hidden
  on_error, 2 ; Return to caller on error.

  catch, error_status
  if error_status ne 0 then begin
    catch, /cancel
    ;
    ;   Clean up.
    ;
    IF N_ELEMENTS(oTestSurface) THEN Obj_Destroy, oTestSurface
    if keyword_set(group_leader_is_fabricated) then begin
      widget_control, group_leader, /destroy
      ptr_free, ptr_new(group_leader, /no_copy)
    endif
    IF N_ELEMENTS(oWindow) THEN Obj_Destroy, oWindow
    IF N_ELEMENTS(oScene) THEN Obj_Destroy, oScene
    if n_elements(oOriginalObj) gt 0 then $
      oObj = oOriginalObj
    if keyword_set(do_destroy_view) then begin
      Obj_Destroy, oCurrentView
    endif
    ;
    ;   Re-throw the error.
    ;
    message, !error_state.msg + ' ' + !error_state.sys_msg
  endif

  if keyword_set(debug) then begin
    on_error, 0
    catch, /cancel
  endif

  if n_elements(refresh_tlb) gt 0 then begin
    if not widget_info(refresh_tlb, /valid) then $
      message, 'Specified widget for REFRESH is invalid.', /noname

    widget_control, refresh_tlb, get_uvalue = pState
    if obj_valid((*pState).oObjviewWid) then $
      (*pState).oObjviewWid->Draw, /hourglass
    return
  end

  if keyword_set(test) then begin
    if n_elements(oObj) gt 0 then $
      oOriginalObj = oObj
    oObj = obj_new('IDLgrSurface', $
      beselj(shift(dist(40), 20, 20) / 2,0) * 20, $
      color=[60, 60, 255], $
      style=2, $
      shading=1, $
      name='Test Surface'$
      )
    oTestSurface = oObj
  endif $
  else begin
    oTestSurface = obj_new()
  endelse

  if n_elements(oObj) eq 0 then begin
    if arg_present(oObj) then $
      message, 'Argument is undefined.', /noname $
    else $
      message, 'requires an argument.', /noname
  endif

  if size(oObj[0], /tname) ne 'OBJREF' then $
    message, 'Argument must be of object reference type.', /noname

  if n_elements(uniq(oObj, sort(oObj))) ne n_elements(oObj) then $
    message, 'Array argument must contain unique values.', /noname

  if obj_isa(oObj[0], 'IDLgrView') then $
    message, 'Note: View arguments are an undocumented feature.', /inform

  if n_elements(oStationary) gt 0 then begin
    if size(oStationary, /tname) ne 'OBJREF' then $
      message, $
      'Keyword STATIONARY must be of object reference type.', $
      /noname
  end

  if n_elements(scale) eq 3 then begin
    ;
    ;   The restriction imposed by the following message is not technically
    ;   necessary, but it serves to simplify xobjview_XDU's command interface.
    ;   It leaves one way of doing independent x, y, z scaling, instead of
    ;   two ways that would interact with each other.
    ;
    message, 'SCALE must be a 1-element value.  (If you desire ' + $
      'different ammounts of scale in x, y & z, use an IDLgrModel.)'
  end

  if n_elements(group_leader) ne 0 then begin
    if not widget_info(group_leader, /valid_id) then begin
      message, 'Specified Group Leader is not valid.', /noname
    endif
  endif $
  else begin
    if keyword_set(modal) then begin
      ;
      ;       Modal widgets require a group leader.  A group leader was not
      ;       specified, so fabricate an invisible one.
      ;
      group_leader = widget_base(map=0)
      group_leader_is_fabricated = 1b
    endif
  endelse
  ;
  ;Create widgets.
  ;
  if keyword_set(modal) then begin
    tlb = widget_base( $
      /column, $
      xpad=0, $
      ypad=0, $
      xoffset=xoffset, $
      yoffset=yoffset, $
      title=n_elements(title) eq 0 ? 'xobjview_XDU' : title, $
      /tlb_size_events, $
      /modal, $
      group_leader=group_leader $
      )
  endif $
  else begin
    tlb = widget_base( $
      /column, $
      xpad=0, $
      ypad=0, $
      xoffset=xoffset, $
      yoffset=yoffset, $
      title=n_elements(title) eq 0 ? 'xobjview_XDU' : title, $
      /tlb_size_events, $
      group_leader=group_leader, $
      mbar=mbar $
      )
  endelse

  oObjviewWid = obj_new('IDLexObjviewWid', $
    tlb, $
    oObj, $
    menu_parent=mbar, $
    scale=scale, $
    draw_xsize=xsize, $
    draw_ysize=ysize, $
    background=background, $
    stationary=oStationary, $
    double_view=double_view, $
    renderer=renderer, $
    use_instancing=use_instancing, $
    /include_refresh_button, $
    /include_full_reset_button, $
    debug=debug $
    )
  if not obj_valid(oObjviewWid) then begin
    message, !error_state.msg + ' ' + !error_state.sys_msg, /noname
  end

  widget_control, tlb, set_uvalue=ptr_new({ $
    oObjviewWid: oObjviewWid, $
    oTestSurface: obj_valid(oTestSurface) ? oTestSurface : obj_new(), $
    pGroupLeader: ptr_new(group_leader), $
    group_leader_is_fabricated: keyword_set(group_leader_is_fabricated) $
  })
  widget_control, tlb, /realize
  xmanager, $
    "xobjview_XDU", $
    tlb, $
    just_reg=keyword_set(just_reg),$
    no_block=keyword_set(block) eq 0, $
    cleanup='xobjview_XDU_cleanup'

  if keyword_set(group_leader_is_fabricated) then begin
    ;
    ;   Leave GROUP_LEADER parameter like we found it: undefined.
    ;
    ptr_free, ptr_new(group_leader, /no_copy)
  endif
  ;
  ;Leave oObj argument unchanged.
  ;
  if n_elements(oOriginalObj) gt 0 then $
    oObj = oOriginalObj

end

;----------------------------------------------------------------------------
;
;  Purpose:  Handle the event.
;
pro d_surfview1Event, $
  sEvent         ; IN: event structure

COMMON GJPOLY, smoothedOutverts, Outconn1, threshold_min_value, threshold_max_value

  WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
  demo_record, sEvent, filename=xd_sState.record_to_filename
  WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

   ;;;;;;;;;;HJX;;;;;;;;;;;;;;;;;;;;; 
 
  uName = WIDGET_INFO(sEvent .id,/uName)
  IF uName EQ 'table1' THEN BEGIN
    ;
    IF sEvent .type EQ 4 THEN BEGIN
      IF sEvent .sel_left NE -1 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
   ;     xd_sState.DATAXYZDISTANCE=string(xd_sState.DATAXYZDISTANCE,format='(f4.2)')
        xd_sState.mDATAXYZDISTANCE = STRING(xd_sState.dataXYZdistance, $
          FORMAT='(F6.2)')
          IF sEvent.SEL_RIGHT LE 1  then begin
        xd_sState.arr2(4*sEvent.SEL_TOP+sEvent.SEL_RIGHT) =xd_sState.mDATAXYZDISTANCE
        print,xd_sState.mDATAXYZDISTANCE
        if xd_sState.arr2(4*sEvent.SEL_TOP+sEvent.SEL_RIGHT) Ge 4 THEN  begin
                xd_sState.arr2(4*sEvent.SEL_TOP+sEvent.SEL_RIGHT+2) ='Yes'  
                ;WIDGET_CONTROL,sEvent .id, set_value =xd_sState.arr2,background_color=[[255,255,255],[255,255,255],[0,255,255],[255,0,255]] 
        endif     
        if xd_sState.arr2(4*sEvent.SEL_TOP+sEvent.SEL_RIGHT) Lt 4 THEN  begin
          xd_sState.arr2(4*sEvent.SEL_TOP+sEvent.SEL_RIGHT+2) ='No'
          ;WIDGET_CONTROL,sEvent .id, set_value =xd_sState.arr2,background_color=[[255,255,255],[255,255,255],[0,255,255],[255,0,255]] 
        endif
          endif      
        WIDGET_CONTROL,sEvent .id, set_value =xd_sState.arr2
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
      ENDIF
    ENDIF
  ENDIF
IF uName NE 'table1' THEN BEGIN  
  
;;;;;;;;;;HJX;;;;;;;;;;;;;;;;;;;;;
    
    

  ;  Quit the application using the close box.
  ;
  if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ $
    'WIDGET_KILL_REQUEST') then begin
    WIDGET_CONTROL, sEvent.top, /DESTROY
    RETURN
  endif

  WIDGET_CONTROL, sEvent.id, GET_UVAL=uval

  case uval of
        
        'ColorIndex': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          WIDGET_CONTROL, xd_sState.wColorIndex, GET_VALUE=ColorIndex_value
          
          ;whether selected is true
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          xd_sState.selected->getProperty, data=verts
          IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;assign color to the selected object
          color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
          ;define a new polygon
          xd_sState.selected->GetProperty, DATA=DATA, hide=hideStatus
          vert_colors = DATA
          vert_colors[0,*] = color_rgb[ColorIndex_value, 0]
          vert_colors[1,*] = color_rgb[ColorIndex_value, 1]
          vert_colors[2,*] = color_rgb[ColorIndex_value, 2]
          xd_sState.selected->setProperty, vert_colors = vert_colors
          WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
          DEVICE, DECOMPOSED=1
          WSET, wColorIndexID
          ERASE, ColorIndexBGR[ColorIndex_value]
          IF hideStatus EQ 1 THEN BEGIN
            plots, [0,40], [0,30], color=[0,0,0], /device
            plots, [0,40], [30,0], color=[0,0,0], /device
          ENDIF
          ;          IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
          xd_sState.tvol_white.SetProperty, color=vert_colors[*,0];, FILL_COLOR=fill_color
          
          ;Hide the thickness colorbar
          xd_sState.colorbarmodel->setProperty, Hide = 1
          
          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ortho2D_plot, xd_sState
          ;draw_3DView, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        
        'ObjectIndex': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          print, 'get_value = ', ObjectIndex_value
          WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')

          ;plot thw colorIndex_draw window
          WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
          DEVICE, DECOMPOSED=1
          WSET, wColorIndexID
          color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            widget_control, xd_sState.wColorIndex, sensitive=0
            
            widget_control, xd_sState.wSmooth, sensitive=0
            widget_control, xd_sState.wSmooth_text, sensitive=0
            ;reset the smooth level value
            WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=0.
            WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(0., format='(f5.1)')
                        
            widget_control, xd_sState.wtransSlider, sensitive=0
            widget_control, xd_sState.wTrans_text, sensitive=0
            ;reset the alpha-channel transparent value
            WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=100.
            WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(100., format='(f5.1)')
            
            widget_control, xd_sState.w3DSegd, sensitive=0
            widget_control, xd_sState.w3Drg, sensitive=0
            widget_control, xd_sState.wthickness, sensitive=0
            widget_control, xd_sState.wpoint, sensitive=0
            widget_control, xd_sState.wwire, sensitive=0
            widget_control, xd_sState.wfill, sensitive=0
;            widget_control, xd_sState.wexpstl, sensitive=0
;            widget_control, xd_sState.whide, sensitive=0
            widget_control, xd_sState.wsinglecenter, sensitive=0
            
            ERASE, ColorIndexBGR[LONG(ObjectIndex_value)]
            xd_sState.selected=OBJ_NEW()
            xd_sState.ogreenaxis1->SetProperty, HIDE=1
            xd_sState.oyellowaxis1->SetProperty, HIDE=1
            xd_sState.oBlueAxis1->SetProperty, HIDE=1
            xd_sState.oS0->SetProperty, HIDE=1
            xd_sState.oS1->SetProperty, HIDE=1
            xd_sState.oS2->SetProperty, HIDE=1
            xd_sState.ckr1->setProperty,hide=1
            xd_sState.ckb->setProperty,hide=1
            xd_sState.ckg->setProperty, hide=1
            xd_sState.ckr->setProperty, hide=1
            xd_sState.ckag->setProperty, hide=1
            xd_sState.oLine1->setProperty,hide=1
            xd_sState.oLine2->setProperty,hide=1
            xd_sState.oLine3->setProperty,hide=1
            xd_sState.tvol_white.SetProperty,hide=1
          ENDIF ELSE BEGIN
            widget_control, xd_sState.wColorIndex, sensitive=1
            widget_control, xd_sState.wSmooth, sensitive=1
            widget_control, xd_sState.wSmooth_text, sensitive=1
            widget_control, xd_sState.wtransSlider, sensitive=1
            widget_control, xd_sState.wTrans_text, sensitive=1
            widget_control, xd_sState.w3DSegd, sensitive=1
            widget_control, xd_sState.w3Drg, sensitive=1
            widget_control, xd_sState.wthickness, sensitive=1
            widget_control, xd_sState.wpoint, sensitive=1
            widget_control, xd_sState.wwire, sensitive=1
            widget_control, xd_sState.wfill, sensitive=1
            widget_control, xd_sState.wexpstl, sensitive=1
;            widget_control, xd_sState.whide, sensitive=1
            widget_control, xd_sState.wsinglecenter, sensitive=1
            
            ERASE, ColorIndexBGR[LONG(ObjectIndex_value)]
            print, 'color = ', ColorIndexBGR[LONG(ObjectIndex_value)]
            xd_sState.selected=xd_sState.oPoly[ObjectIndex_value-1]
            xd_sState.selected->GetProperty, DATA = DATA, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel, hide=hideStatus, THICK=thickness
            
            ;set the alpha-channel transparent value
            WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
            WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')
            
            ;thick was used to record smoothness
            smoothness = (thickness-FLOOR(thickness))*1000
            WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=smoothness
            WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothness, format='(f5.1)')
            
            xd_sState.tvol_white.SetProperty,hide=0
            count = xd_sState.struc_modifyobject.countmask_list[LONG(ObjectIndex_value)-1]
            xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
            Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
;            IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
            IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
                        
            dataxyz = REFORM(DATA[*, N_ELEMENTS(DATA(0,*))/2])

            ;Enable color selection
            WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1

            ;the object is multi-color, we don't need to do anything
            ;if the object is single color, we will find right color index
            IF STDDEV(vert_colors[0,*]) LT 1 THEN BEGIN
;              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
              FOR m=0, N_ELEMENTS(ColorIndex)-1 DO BEGIN
                IF vert_colors[0,0] EQ color_rgb[m,0] AND vert_colors[1,0] EQ color_rgb[m,1] AND vert_colors[2,0] EQ color_rgb[m,2] THEN colorindexS = m
              ENDFOR
              IF N_ELEMENTS(colorindexS) EQ 0 THEN colorindexS = 0
              WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=colorindexS
              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              ERASE, ColorIndexBGR[colorindexS]
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDIF ELSE BEGIN
              
              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              ;use multiple lines to label multiple colro
              plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
              plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
              plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDELSE
            
            ;GJ 2019/2/22, summarize all plot pickdata into a program
            PLOT_pickdata, dataxyz, xd_sState
          ENDELSE

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
          ;;draw_3DView, xd_sState
          ;GJ 2019/5/19; disable the draw_3DView for object index

          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        'ObjectIndex_text': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          WIDGET_CONTROL, xd_sState.wObjectIndex_text, GET_VALUE=temp
          ObjectIndex_value = LONG(temp)
          n_Poly = xd_sState.struc_modifyobject.exist_num[xd_sState.struc_modifyobject.modify_Step]
          IF ObjectIndex_value GT n_Poly THEN ObjectIndex_value = n_Poly
          IF ObjectIndex_value LT 0 THEN ObjectIndex_value = 0
          print, 'get_value = ', ObjectIndex_value
          WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
          WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value
;          IF LONG(ObjectIndex_value) GT 0 THEN BEGIN
;            widget_control, xd_sState.wColorIndex, sensitive=1
;            WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=LONG(ObjectIndex_value)
;          ENDIF ELSE BEGIN
;            widget_control, xd_sState.wColorIndex, sensitive=0
;            WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=LONG(ObjectIndex_value)
;          ENDELSE
;          WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=LONG(ObjectIndex_value)
          ;plot thw colorIndex_draw window
          WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
          DEVICE, DECOMPOSED=1
          WSET, wColorIndexID
          color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            widget_control, xd_sState.wColorIndex, sensitive=0
            widget_control, xd_sState.wSmooth, sensitive=0
            widget_control, xd_sState.wSmooth_text, sensitive=0
            widget_control, xd_sState.wtransSlider, sensitive=0
            widget_control, xd_sState.wTrans_text, sensitive=0
            widget_control, xd_sState.w3DSegd, sensitive=0
            widget_control, xd_sState.w3Drg, sensitive=0
            widget_control, xd_sState.wthickness, sensitive=0
            widget_control, xd_sState.wpoint, sensitive=0
            widget_control, xd_sState.wwire, sensitive=0
            widget_control, xd_sState.wfill, sensitive=0
;            widget_control, xd_sState.wexpstl, sensitive=0
;            widget_control, xd_sState.whide, sensitive=0
            widget_control, xd_sState.wsinglecenter, sensitive=0
                        
            ERASE, ColorIndexBGR[LONG(ObjectIndex_value)]
            xd_sState.selected=OBJ_NEW()
            xd_sState.ogreenaxis1->SetProperty, HIDE=1
            xd_sState.oyellowaxis1->SetProperty, HIDE=1
            xd_sState.oBlueAxis1->SetProperty, HIDE=1
            xd_sState.oS0->SetProperty, HIDE=1
            xd_sState.oS1->SetProperty, HIDE=1
            xd_sState.oS2->SetProperty, HIDE=1
            xd_sState.ckr1->setProperty,hide=1
            xd_sState.ckb->setProperty,hide=1
            xd_sState.ckg->setProperty, hide=1
            xd_sState.ckr->setProperty, hide=1
            xd_sState.ckag->setProperty, hide=1
            xd_sState.oLine1->setProperty,hide=1
            xd_sState.oLine2->setProperty,hide=1
            xd_sState.oLine3->setProperty,hide=1
            xd_sState.tvol_white.SetProperty,hide=1
          ENDIF ELSE BEGIN
            widget_control, xd_sState.wColorIndex, sensitive=1
            widget_control, xd_sState.wSmooth, sensitive=1
            widget_control, xd_sState.wSmooth_text, sensitive=1
            widget_control, xd_sState.wtransSlider, sensitive=1
            widget_control, xd_sState.wTrans_text, sensitive=1
            widget_control, xd_sState.w3DSegd, sensitive=1
            widget_control, xd_sState.w3Drg, sensitive=1
            widget_control, xd_sState.wthickness, sensitive=1
            widget_control, xd_sState.wpoint, sensitive=1
            widget_control, xd_sState.wwire, sensitive=1
            widget_control, xd_sState.wfill, sensitive=1
            widget_control, xd_sState.wexpstl, sensitive=1
;            widget_control, xd_sState.whide, sensitive=1
            widget_control, xd_sState.wsinglecenter, sensitive=1
            
            ERASE, ColorIndexBGR[LONG(ObjectIndex_value)]
            print, 'color = ', ColorIndexBGR[LONG(ObjectIndex_value)]
            xd_sState.selected=xd_sState.oPoly[ObjectIndex_value-1]
            xd_sState.selected->GetProperty, DATA = DATA, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel, hide=hideStatus, THICK=thickness

            ;set the alpha-channel transparent value
            WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
            WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

            ;thick was used to record smoothness
            smoothness = (thickness-FLOOR(thickness))*1000
            WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=smoothness
            WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothness, format='(f5.1)')
            
            dataxyz = REFORM(DATA[*, N_ELEMENTS(DATA(0,*))/2])

            xd_sState.tvol_white.SetProperty,hide=0
            count = xd_sState.struc_modifyobject.countmask_list[LONG(ObjectIndex_value)-1]
            xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
            Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
;            IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
            IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
            
            ;Enable color selection
            WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1

            ;the object is multi-color, we don't need to do anything
            ;if the object is single color, we will find right color index
            IF STDDEV(vert_colors[0,*]) LT 1 THEN BEGIN
;              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
              FOR m=0, N_ELEMENTS(ColorIndex)-1 DO BEGIN
                IF vert_colors[0,0] EQ color_rgb[m,0] AND vert_colors[1,0] EQ color_rgb[m,1] AND vert_colors[2,0] EQ color_rgb[m,2] THEN colorindexS = m
              ENDFOR
              IF N_ELEMENTS(colorindexS) EQ 0 THEN colorindexS = 0
              WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=colorindexS
              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              ERASE, ColorIndexBGR[colorindexS]
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDIF ELSE BEGIN
              
              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              ;use multiple lines to label multiple colro
              plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
              plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
              plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDELSE
            
            ;GJ 2019/2/22, summarize all plot pickdata into a program
            PLOT_pickdata, dataxyz, xd_sState
          ENDELSE

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
          ;;draw_3DView, xd_sState
          ;GJ, 2019/5/19;disable the draw_3DView
          
          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        ;this is for colortable display, doing nothing
        'ColorIndexlabel2': begin
        end
        
        'singlecenter': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

          ;whether selected is true
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN 
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;center the selected object
          xd_sState.selected->getProperty, data=verts
          IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          Adjust_whole_polygons, verts, xd_sState

          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        'wholecenter': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          ;adjust the center of all 3D objects
          n_Poly = xd_sState.struc_modifyobject.exist_num[xd_sState.struc_modifyobject.modify_Step]
          FOR i = 0, n_Poly-1 DO BEGIN
            xd_sstate.opoly[i]->getProperty, data=verts
            IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
              WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
              RETURN
            ENDIF
            IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
          ENDFOR
          ;center the selected object
          Adjust_whole_polygons, whole_verts, xd_sState

          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
                
        'singlehide': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          if (sEvent.type EQ 0) then begin
            ;whether selected is true
            IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
              WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
              RETURN
            ENDIF

            ;whether selected is a real idlgrpolygon
            IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
              WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
              RETURN
            ENDIF

            xd_sState.selected->getProperty, data=verts, NAME=name_poly, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel, hide=hideStatus, THICK=thickness
            
            IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
              WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
              RETURN
            ENDIF
            
            ;set the alpha-channel transparent value
            WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
            WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')
            
            ;thick was used to record smoothness
            smoothness = (thickness-FLOOR(thickness))*1000
            WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=smoothness
            WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothness, format='(f5.1)')

            ;update modifyStep and current_step_index
            
            modify_step = xd_sState.struc_modifyobject.current_step_index + 1
            xd_sState.struc_modifyobject.modify_step =  modify_step
            xd_sState.struc_modifyobject.current_step_index = modify_step
            xd_sState.struc_modifyobject.exist_num[modify_step] = xd_sState.struc_modifyobject.exist_num[modify_step-1]
            widget_control, xd_sState.wUndo, SENSITIVE = 1
            widget_control, xd_sState.wRedo, SENSITIVE = 0
            ;hideStatus = xd_sState.struc_modifyobject.hide_status[modifyStep, name_poly]

            IF hideStatus EQ 0 THEN BEGIN
              xd_sState.selected->setProperty, hide = 1
              xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
              xd_sState.struc_modifyobject.hide_status[modify_step, name_poly] = 1
              
              ;hide the total volume
              xd_sState.tvol_white.SetProperty,hide=1
              
              ;disable these options
              widget_control, xd_sState.wColorIndex, sensitive=0
              widget_control, xd_sState.wSmooth, sensitive=0
              widget_control, xd_sState.wSmooth_text, sensitive=0
              widget_control, xd_sState.wtransSlider, sensitive=0
              widget_control, xd_sState.wTrans_text, sensitive=0
              widget_control, xd_sState.w3DSegd, sensitive=0
              widget_control, xd_sState.w3Drg, sensitive=0
              widget_control, xd_sState.wthickness, sensitive=0
              widget_control, xd_sState.wpoint, sensitive=0
              widget_control, xd_sState.wwire, sensitive=0
              widget_control, xd_sState.wfill, sensitive=0
;             widget_control, xd_sState.wexpstl, sensitive=0
;             widget_control, xd_sState.whide, sensitive=0
              widget_control, xd_sState.wsinglecenter, sensitive=0
            ENDIF ELSE BEGIN
              xd_sState.selected->setProperty, hide = 0
              xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
              xd_sState.struc_modifyobject.hide_status[modify_step, name_poly] = 0
              
              ;hide the total volume
              xd_sState.tvol_white.SetProperty,hide=0
              
              ;enable these options
              widget_control, xd_sState.wColorIndex, sensitive=1
              widget_control, xd_sState.wSmooth, sensitive=1
              widget_control, xd_sState.wSmooth_text, sensitive=1
              widget_control, xd_sState.wtransSlider, sensitive=1
              widget_control, xd_sState.wTrans_text, sensitive=1
              widget_control, xd_sState.w3DSegd, sensitive=1
              widget_control, xd_sState.w3Drg, sensitive=1
              widget_control, xd_sState.wthickness, sensitive=1
              widget_control, xd_sState.wpoint, sensitive=1
              widget_control, xd_sState.wwire, sensitive=1
              widget_control, xd_sState.wfill, sensitive=1
              widget_control, xd_sState.wexpstl, sensitive=1
              ;            widget_control, xd_sState.whide, sensitive=0
              widget_control, xd_sState.wsinglecenter, sensitive=1
            ENDELSE
            
            ;the object is multi-color, we don't need to do anything
            ;if the object is single color, we will find right color index
            xd_sState.selected->getProperty, hide = hideStatus
            IF STDDEV(vert_colors[0,*]) LT 1 THEN BEGIN
              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
              FOR m=0, N_ELEMENTS(ColorIndex)-1 DO BEGIN
                IF vert_colors[0,0] EQ color_rgb[m,0] AND vert_colors[1,0] EQ color_rgb[m,1] AND vert_colors[2,0] EQ color_rgb[m,2] THEN colorindexS = m
              ENDFOR
              IF N_ELEMENTS(colorindexS) EQ 0 THEN colorindexS = 0
              WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=colorindexS
              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              ERASE, ColorIndexBGR[colorindexS]
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDIF ELSE BEGIN

              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
              WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
              ERASE, ColorIndexBGR[LONG(ObjectIndex_value)]
              ;use multiple lines to label multiple colro
              plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
              plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
              plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDELSE
            
            ;redraw and display 3D naked eye
            xd_sState.oWindow->Draw, xd_sState.oView
            ;draw_3DView, xd_sState
            ortho2D_plot, xd_sState
          endif  
       
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end

        'pointstyle': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          
          ;whether selected is true
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;whether selected is a real idlgrpolygon
          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN 
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;use point to draw polygon
          xd_sState.selected->setProperty, style = 0
          
          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        'wirestyle': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          
          ;whether selected is true
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;whether selected is a real idlgrpolygon
          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN 
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;use wire to draw polygon
          xd_sState.selected->setProperty, style = 1
          
          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end

        'fillstyle': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          
          ;whether selected is true
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;whether selected is a real idlgrpolygon
          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN 
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;use point to draw polygon
          xd_sState.selected->setProperty, style = 2
          
          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        'expstl': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

          ;whether selected is true
          print, 'OBJ_VALID(xd_sState.selected) = ', OBJ_VALID(xd_sState.selected)
          print, "obj_isa(xd_sState.selected,'IDLgrPolygon') = ", obj_isa(xd_sState.selected,'IDLgrPolygon')
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            ;save all shown IDLgrPolygons
            curr_modify_Step = xd_sState.struc_modifyobject.current_step_index
            n_Poly = xd_sState.struc_modifyobject.exist_num[curr_modify_Step]
            hide_status = xd_sState.struc_modifyobject.hide_status[curr_modify_step, *]
            count = 0
            FOR i = 0, n_Poly-1 DO BEGIN
              IF hide_status[i] EQ 0 THEN BEGIN
                xd_sState.opoly[i]->getProperty, data=verts1, POLYGONS=conn1
                IF count EQ 0 THEN BEGIN
                  Verts = Verts1
                  Conn = Conn1
                  count = count + 1
                ENDIF ELSE BEGIN
                  Result = MESH_MERGE(Verts, Conn, Verts1, Conn1)
                  count = count + 1
                ENDELSE
              ENDIF
            ENDFOR

            outputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
              TITLE='Select stl File', FILTER='*.stl')
            IF STRLEN(outputFile) GT 0 AND STRCMP(STRMID(outputFile, STRLEN(outputFile)-4, 4), '.stl') NE 1 THEN outputFile = outputFile + '.stl'

            ;use mesh_decimate to save stl file
            IF STRLEN(outputFile) GT 0 THEN BEGIN
              numberVertices = MESH_DECIMATE(Verts, Conn, p2, VERTICES = v2)
              writestl,outputfile,v2*xd_sState.vol_HU_cube_resolution,p2
            ENDIF
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;whether selected is a real idlgrpolygon
          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
           ;save all shown IDLgrPolygons
            curr_modify_Step = xd_sState.struc_modifyobject.current_step_index
            n_Poly = xd_sState.struc_modifyobject.exist_num[curr_modify_Step]
            hide_status = xd_sState.struc_modifyobject.hide_status[curr_modify_step, *]
            count = 0
            FOR i = 0, n_Poly-1 DO BEGIN
              IF hide_status[i] EQ 0 THEN BEGIN
                xd_sState.opoly[i]->getProperty, data=verts1, POLYGONS=conn1
                IF count EQ 0 THEN BEGIN
                  Verts = Verts1
                  Conn = Conn1
                  count = count + 1
                ENDIF ELSE BEGIN
                  Result = MESH_MERGE(Verts, Conn, Verts1, Conn1)
                  count = count + 1
                ENDELSE
              ENDIF
            ENDFOR

            outputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
              TITLE='Select stl File', FILTER='*.stl')
            IF STRLEN(outputFile) GT 0 AND STRCMP(STRMID(outputFile, STRLEN(outputFile)-4, 4), '.stl') NE 1 THEN outputFile = outputFile + '.stl'

            ;use mesh_decimate to save stl file
            IF STRLEN(outputFile) GT 0 THEN BEGIN
              numberVertices = MESH_DECIMATE(Verts, Conn, p2, VERTICES = v2)
              writestl,outputfile,v2*xd_sState.vol_HU_cube_resolution,p2
            ENDIF
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;pick vert colors based on a new matrix ;GJ 2019/2/13         
          xd_sState.selected->getProperty, data=verts, POLYGONS=Outconn
          IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
            ;save all shown IDLgrPolygons
            curr_modify_Step = xd_sState.struc_modifyobject.current_step_index
            n_Poly = xd_sState.struc_modifyobject.exist_num[curr_modify_Step]
            hide_status = xd_sState.struc_modifyobject.hide_status[curr_modify_step, *]
            count = 0
            FOR i = 0, n_Poly-1 DO BEGIN
              IF hide_status[i] EQ 0 THEN BEGIN
                xd_sState.opoly[i]->getProperty, data=verts1, POLYGONS=conn1
                IF count EQ 0 THEN BEGIN
                  Verts = Verts1
                  Conn = Conn1
                  count = count + 1
                ENDIF ELSE BEGIN
                  Result = MESH_MERGE(Verts, Conn, Verts1, Conn1)
                  count = count + 1
                ENDELSE
              ENDIF
            ENDFOR

            outputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
              TITLE='Select stl File', FILTER='*.stl')
            IF STRLEN(outputFile) GT 0 AND STRCMP(STRMID(outputFile, STRLEN(outputFile)-4, 4), '.stl') NE 1 THEN outputFile = outputFile + '.stl'

            ;use mesh_decimate to save stl file
            IF STRLEN(outputFile) GT 0 THEN BEGIN
              numberVertices = MESH_DECIMATE(Verts, Conn, p2, VERTICES = v2)
              writestl,outputfile,v2*xd_sState.vol_HU_cube_resolution,p2
            ENDIF
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          outputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
            TITLE='Select stl File', FILTER='*.stl')
          IF STRLEN(outputFile) GT 0 AND STRCMP(STRMID(outputFile, STRLEN(outputFile)-4, 4), '.stl') NE 1 THEN outputFile = outputFile + '.stl'
          
          ;use mesh_decimate to save stl file
          IF STRLEN(outputFile) GT 0 THEN BEGIN
            numberVertices = MESH_DECIMATE(verts, Outconn, p2, VERTICES = v2)
            writestl,outputfile,v2*xd_sState.vol_HU_cube_resolution,p2
          ENDIF

          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, /HOURGLASS
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        'eccentricity': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;whether selected is a real idlgrpolygon
          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;pick vert colors based on a new matrix ;GJ 2019/2/13
          ;modify this algorithm, GJ 2020/5/22
          xd_sState.selected->getProperty, data=verts, POLYGONS=conn
          IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF

          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE

          ;this is the fastest calculation by smooth
          ;print, "additional dir: ", xd_sState.additionalDir
          ;central_loc = CENTRAL_LOC_CAL(xd_sState.additionalDir[0])
          vert_colors = PICK_VERT_COLORS_by_hough_segmentation(mask_vol_cube, verts, conn, rgb_table)


;          vert_colors = PICK_VERT_COLORS_by_eccentricity(mask_vol_cube, verts, rgb_table)
          ;GJ 2019/5/10
          ;vert_colors = PICK_VERT_COLORS_by_region_grow(mask_vol_cube, verts, rgb_table)
          xd_sState.selected->setProperty, vert_colors = vert_colors
          xd_sState.colorbarmodel->setProperty, Hide = 0
          WIDGET_CONTROL, /HOURGLASS

          ;Add color lines to label color
          color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
          WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
          DEVICE, DECOMPOSED=1
          WSET, wColorIndexID
          ;use multiple lines to label multiple colro
          plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
          plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
          plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device

          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
        
        'thickness': begin
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
          
          IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;whether selected is a real idlgrpolygon
          IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN 
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
 
          ;pick vert colors based on a new matrix ;GJ 2019/2/13
          ;modify this algorithm, GJ 2020/5/22         
          xd_sState.selected->getProperty, data=verts, POLYGONS=conn
          IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
            WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
            RETURN
          ENDIF
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;this is the fastest calculation by smooth
;          print, "additional dir: ", xd_sState.additionalDir
;          central_loc = CENTRAL_LOC_CAL(xd_sState.additionalDir[0])
;          vert_colors = PICK_VERT_COLORS_by_thickness_smooth(mask_vol_cube, central_loc, verts, conn, rgb_table)
          
          ;this is the old one
          vert_colors = PICK_VERT_COLORS_by_thickness(mask_vol_cube, verts, rgb_table)
          ;GJ 2020/5/23
          ;vert_colors = PICK_VERT_COLORS_by_vessel_thickness(verts, conn, rgb_table)
          ;GJ 2019/5/10
          ;vert_colors = PICK_VERT_COLORS_by_region_grow(mask_vol_cube, verts, rgb_table)
          xd_sState.selected->setProperty, vert_colors = vert_colors
          xd_sState.colorbarmodel->setProperty, Hide = 0
          WIDGET_CONTROL, /HOURGLASS
          
          ;Add color lines to label color
          color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
          WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
          DEVICE, DECOMPOSED=1
          WSET, wColorIndexID
          ;use multiple lines to label multiple colro
          plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
          plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
          plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device
      
          ;redraw and display 3D naked eye
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        end
    
;    ;;;;GJ, 2018/12/19
;    'xscaleslider':begin
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;
;      widget_control, xd_sState.x_scale_slider, get_value = xscale
;      widget_control, xd_sState.y_scale_slider, get_value = yscale
;      widget_control, xd_sState.z_scale_slider, get_value = zscale
;      widget_control, xd_sState.x_move_slider, get_value = xmove
;      widget_control, xd_sState.y_move_slider, get_value = ymove
;      widget_control, xd_sState.z_move_slider, get_value = zmove
;  
;      print, 'x scale = ', xscale
;      threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
;      volume_temp = xd_sState.red_mask_vol_cube
;      Transform_Volume, volume_temp, new_vol, Scale=[xscale, 1, 1]
;      SHADE_VOLUME, new_vol,threshold_min_value, Outverts, Outconn1
;
;      IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
;        index = WHERE(new_vol GT threshold_min_value, count)
; ;       IF count NE 0 THEN xd_sState.vol_red = count * (xd_sState.vol_HU_cube_resolution^3)/1000.
; ;       print, 'volume (red) = ', xd_sState.vol_red, ' cm3'
; ;       Stringred = STRING(xd_sState.vol_red,FORMAT='("vol(red) =", F9.3, " cm3")')
; ;       xd_sState.tvol_red.SetProperty,strings=stringred
;
;        smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1, ITERATIONS=400)
;        ;plot segmentation result
;        xd_sState.oPoly1->getProperty, DATA = Outverts_old
;        ;define to the rotation center
;        IF N_ELEMENTS(Outverts_old) LT N_ELEMENTS(smoothedOutverts) THEN BEGIN
;          xd_sState.oPoly1->setProperty, DATA = smoothedOutverts
;          xd_sState.oPoly1->setProperty, POLYGONS = Outconn1
;        ENDIF ELSE BEGIN
;          xd_sState.oPoly1->setProperty, POLYGONS = Outconn1
;          xd_sState.oPoly1->setProperty, DATA = smoothedOutverts
;        ENDELSE
;      
;        widget_control, xd_sState.x_scale_slider, set_value = 1
;        xd_sState.oPoly1->setProperty, ALPHA_CHANNEL = 0.5
;        xd_sState.oWindow->Draw, xd_sState.oView
;        ortho2D_plot, xd_sState
;        ;draw_3DView, xd_sState
;      ENDIF
;
;      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;    end
;    
;    
;    'yscaleslider':begin
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;  
;      widget_control, xd_sState.x_scale_slider, get_value = xscale
;      widget_control, xd_sState.y_scale_slider, get_value = yscale
;      widget_control, xd_sState.z_scale_slider, get_value = zscale
;      widget_control, xd_sState.x_move_slider, get_value = xmove
;      widget_control, xd_sState.y_move_slider, get_value = ymove
;      widget_control, xd_sState.z_move_slider, get_value = zmove
;  
;      print, 'x move = ', xmove
;  
;      xd_sState.oWindow->Draw, xd_sState.oView
;      ortho2D_plot, xd_sState
;      ;draw_3DView, xd_sState
;      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;
;
;    end
;
;
;    'zscaleslider':begin
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;
;      widget_control, xd_sState.x_scale_slider, get_value = xscale
;      widget_control, xd_sState.y_scale_slider, get_value = yscale
;      widget_control, xd_sState.z_scale_slider, get_value = zscale
;      widget_control, xd_sState.x_move_slider, get_value = xmove
;      widget_control, xd_sState.y_move_slider, get_value = ymove
;      widget_control, xd_sState.z_move_slider, get_value = zmove
;  
;      print, 'x move = ', xmove
;      
;      xd_sState.oWindow->Draw, xd_sState.oView
;      ortho2D_plot, xd_sState
;      ;draw_3DView, xd_sState
;      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;
;
;    end


;    'xmoveslider':begin
;       WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;
;
;       widget_control, xd_sState.x_scale_slider, get_value = xscale
;       widget_control, xd_sState.y_scale_slider, get_value = yscale
;       widget_control, xd_sState.z_scale_slider, get_value = zscale
;       widget_control, xd_sState.x_move_slider, get_value = xmove
;       widget_control, xd_sState.y_move_slider, get_value = ymove
;       widget_control, xd_sState.z_move_slider, get_value = zmove
;       
;       print, 'x move = ', xmove
;       xd_sState.oPoly1->getProperty, DATA = Outverts_old
;       Outverts_new = Outverts_old
;       Outverts_new[0,*] = Outverts_old[0,*] + xmove
;       xd_sState.oPoly1->setProperty, DATA = Outverts_new
;
;       widget_control, xd_sState.x_move_slider, set_value = 0
;
;       xd_sState.oPoly1->setProperty, ALPHA_CHANNEL = 0.5
;       xd_sState.oWindow->Draw, xd_sState.oView
;       ortho2D_plot, xd_sState
;       ;draw_3DView, xd_sState
;       WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;
;
;    end
;
;
;    'ymoveslider':begin
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;
;      widget_control, xd_sState.x_scale_slider, get_value = xscale
;      widget_control, xd_sState.y_scale_slider, get_value = yscale
;      widget_control, xd_sState.z_scale_slider, get_value = zscale
;      widget_control, xd_sState.x_move_slider, get_value = xmove
;      widget_control, xd_sState.y_move_slider, get_value = ymove
;      widget_control, xd_sState.z_move_slider, get_value = zmove
;  
;      print, 'y move = ', ymove
;      threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
;      xd_sState.oPoly1->getProperty, DATA = Outverts
;      minx=MIN(Outverts[0,*])
;      maxx=MAX(Outverts[0,*])
;      miny=MIN(Outverts[1,*])
;      maxy=MAX(Outverts[1,*])
;      minz=MIN(Outverts[2,*])
;      maxz=MAX(Outverts[2,*])
;      volume_temp = xd_sState.red_mask_vol_cube
;      volume_temp[minx:maxx, miny:maxy, minz:maxz] = threshold_min_value+200.
;      Transform_Volume, volume_temp, new_vol, Translate=[0, ymove, 0]
;      SHADE_VOLUME, new_vol,threshold_min_value, Outverts, Outconn1
;
;      IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
;        index = WHERE(new_vol GT threshold_min_value, count)
; ;       IF count NE 0 THEN xd_sState.vol_red = count * (xd_sState.vol_HU_cube_resolution^3)/1000.
; ;       print, 'volume (red) = ', xd_sState.vol_red, ' cm3'
; ;       Stringred = STRING(xd_sState.vol_red,FORMAT='("vol(red) =", F9.3, " cm3")')
; ;       xd_sState.tvol_red.SetProperty,strings=stringred
;
;        smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1, ITERATIONS=400)
;        ;plot segmentation result
;        xd_sState.oPoly1->getProperty, DATA = Outverts_old
;        ;define to the rotation center
;        IF N_ELEMENTS(Outverts_old) LT N_ELEMENTS(smoothedOutverts) THEN BEGIN
;          xd_sState.oPoly1->setProperty, DATA = smoothedOutverts
;          xd_sState.oPoly1->setProperty, POLYGONS = Outconn1
;        ENDIF ELSE BEGIN
;          xd_sState.oPoly1->setProperty, POLYGONS = Outconn1
;          xd_sState.oPoly1->setProperty, DATA = smoothedOutverts
;        ENDELSE
;
;        xd_sState.oPoly1->setProperty, ALPHA_CHANNEL = 0.5
;        widget_control, xd_sState.y_move_slider, set_value = 0
;        xd_sState.oWindow->Draw, xd_sState.oView
;        ortho2D_plot, xd_sState
;        ;draw_3DView, xd_sState
;      ENDIF
;  
;      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;    end
;
;
;    'zmoveslider':begin
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;
;      widget_control, xd_sState.x_scale_slider, get_value = xscale
;      widget_control, xd_sState.y_scale_slider, get_value = yscale
;      widget_control, xd_sState.z_scale_slider, get_value = zscale
;      widget_control, xd_sState.x_move_slider, get_value = xmove
;      widget_control, xd_sState.y_move_slider, get_value = ymove
;      widget_control, xd_sState.z_move_slider, get_value = zmove
;  
;      print, 'x move = ', xmove
;  
;      xd_sState.oWindow->Draw, xd_sState.oView
;      ortho2D_plot, xd_sState
;      ;draw_3DView, xd_sState
;      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;
;
;    end
;    
;    
;    'xyzscalemove':begin
;      ;;zy
;      
;      
;    end
    
    'OPOLY1DEL':begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      
      xd_sState.oPoly1->getProperty, DATA = Outverts_old
      xd_sState.oPoly1->setProperty, DATA = Outverts_old*0.0001
;      xd_sState.vol_red = 0.
;      print, 'volume (red) = ', xd_sState.vol_red, ' cm3'
;      Stringred = STRING(xd_sState.vol_red,FORMAT='("vol(red) =", F9.3, " cm3")')
;      xd_sState.tvol_red.SetProperty,strings=stringred
      
;      WIDGET_CONTROL, xd_sState.wSmooth1, sensitive = 0
;      WIDGET_CONTROL, xd_sState.wtransSlider1, sensitive = 0
;      WIDGET_CONTROL, xd_sState.wloadctSlider1, sensitive = 0
      
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    ;;;;GJ, 2018/12/19

    'Save' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
       txtfile= DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', FILTER='*.xls',TITLE="Choose txtfile saving data.")
       IF STRLEN(txtfile) NE 0 THEN BEGIN
         arr  = xd_sState.arr2
         openw,lun,txtfile,/Get_lun
         for i=0,4 do printf,lun,strtrim(arr[0,i],2)+string(9b)+strtrim(arr[1,i],2)+string(9b)+strtrim(arr[2,i],2)+string(9b)+strtrim(arr[3,i],2)+string(9b)
         free_lun,lun
       ENDIF
       WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
     end

    'Saveall': begin
     WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
     
     date_time=JULDAY()
     filename = STRING(date_time, format = '(C(CMOI2.2, "_", CDI2.2, "_", CYI))')
      file = dialog_pickfile(PATH='./test_data/S0002_162200/', /write, $
        file=filename+'.sav', $
        filter='*.sav' $
        )
      if (file NE '') then begin
        temp_new_mask_vol_cube = xd_sState.new_mask_vol_cube
        temp_mask_vol_cube = xd_sState.mask_vol_cube
        temp_mask_vol_cube_modify = xd_sState.mask_vol_cube_modify
        temp_resolution = xd_sState.vol_HU_cube_resolution
        xd_sstate.opoly[0]->getProperty, POLYGONS = Outconn1
        xd_sstate.opoly[0]->getProperty, DATA = smoothedOutverts
        SAVE, temp_new_mask_vol_cube, temp_mask_vol_cube, temp_mask_vol_cube_modify, temp_resolution, Outconn1, smoothedOutverts, FILENAME = file
      endif
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'Rotate360': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      redvalue = replicate(245,256,1)
      greenvalue = BINDGEN(256)
      bluevalue = reverse(greenvalue)
      FOR i=0, 180 DO BEGIN
        ;GJ 2019/1/28
        axis = [0,1,0]
        angle = 2
        xd_sState.oScalingModel->Rotate, axis, angle
        wloadctslidervalue = FLOOR(i*10.*256./180.) MOD 256
        xd_sstate.opoly[0]->setProperty, color = [redvalue[wloadctslidervalue],greenvalue[wloadctslidervalue],bluevalue[wloadctslidervalue]]
        xd_sState.oWindow->Draw, xd_sState.oView
        ;end GJ 2019/1/28
        draw_3DView, xd_sState
        wait, 0.02
      ENDFOR
      
      ortho2D_plot, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'Restoreall': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      file = dialog_pickfile(PATH='./test_data/S0002_162200/', filter='*.sav')
      if (file NE '') then begin
        RESTORE, FILENAME = file;, RESTORED_OBJECTS = xd_sState
        xd_sState.new_mask_vol_cube = temp_new_mask_vol_cube
        xd_sState.mask_vol_cube = temp_mask_vol_cube
        xd_sState.mask_vol_cube_modify = temp_mask_vol_cube_modify
        xd_sState.vol_HU_cube_resolution = temp_resolution
        xd_sstate.opoly[0]->setProperty, POLYGONS = Outconn1
        xd_sstate.opoly[0]->setProperty, DATA = smoothedOutverts
        verts = smoothedOutverts
        xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
        xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
        xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
        verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
          [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
          [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
          [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
          [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
          [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
          [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
          [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
        xd_sState.oPoly_border->setProperty, DATA = verts_border

        ;define to the rotation center
        ma=max(smoothedOutverts)
        mi=min(smoothedOutverts)
        bias=(ma-mi)/2
        rl=mi+bias
        rr=ma+bias
        cc=[(-rl)/(rr-rl),1/(rr-rl)]
        xma = max((smoothedOutverts)[0,*])
        xmi = min((smoothedOutverts)[0,*])
        xmid = 0.5*(xma+xmi)*cc[1]
        yma = max((smoothedOutverts)[1,*])
        ymi = min((smoothedOutverts)[1,*])
        ymid = 0.5*(yma+ymi)*cc[1]
        zma = max((smoothedOutverts)[2,*])
        zmi = min((smoothedOutverts)[2,*])
        zmid = 0.5*(zma+zmi)*cc[1]

        xd_sstate.opoly[0]->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oPoly1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oPoly2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oPoly_border->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oS0->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oS1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oS2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
        xd_sState.ogreenaxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oyellowaxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oblueAxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
        xd_sState.oLine1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oLine2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        xd_sState.oLine3->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
        
        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        xd_sState.oPoly_border->SetProperty, STYLE = 0
        xd_sState.oWindow->Draw, xd_sState.oView
        ortho2D_plot, xd_sState
        ;draw_3DView, xd_sState
      endif
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    
    ;  Animate the surface accordingly to scenario number 0.
    ;
    'ANIMATE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      d_surfviewAnimateSurface, xd_sState, 0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ; of ANIMATE
    
     
   ;;;;;;;;hxn1-start;;;;;;;;
    'rdstl': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ; Allow the user to select a DICOM file.
      inputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
        TITLE='Select stl File', FILTER='*.stl', /MUST_EXIST)
;      inputfile = 'E:\Research\MED_visualization\6_spineCases\Li  guizhen\LI GUIZHEN.stl'
      IF STRLEN(inputFile) GT 0 THEN BEGIN
        readstl,inputfile,DATA_0,POLYGONS_0
        POLYGONS_0[0,*] = POLYGONS_0[0,*]-MEAN(POLYGONS_0[0,*])
        POLYGONS_0[1,*] = POLYGONS_0[1,*]-MEAN(POLYGONS_0[1,*])
        POLYGONS_0[2,*] = POLYGONS_0[2,*]-MEAN(POLYGONS_0[2,*])
        result = MESH_DECIMATE(DATA_0, POLYGONS_0, Outconn1, vertices=smoothedOutverts)
        print, 'finish loading'
        Add_new_by_readstl, xd_sState, Outconn1, smoothedOutverts
        
        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        xd_sState.oWindow->Draw, xd_sState.oView
        ;draw_3DView, xd_sState

        ;replot the 3 draw windows
        ortho2D_plot, xd_sState
      ENDIF
      
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end
    
    'UNDO': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS

      IF xd_sState.struc_modifyobject.current_step_index GE 1 AND xd_sState.struc_modifyobject.current_step_index LE xd_sState.struc_modifyobject.modify_step THEN BEGIN
        xd_sState.struc_modifyobject.current_step_index -= 1
        current_step_index = xd_sState.struc_modifyobject.current_step_index
        IF current_step_index LE 0 THEN BEGIN
          widget_control, xd_sState.wUndo, SENSITIVE = 0
          widget_control, xd_sState.wRedo, SENSITIVE = 1
        ENDIF ELSE BEGIN
          widget_control, xd_sState.wUndo,  SENSITIVE = 1
          widget_control, xd_sState.wRedo,  SENSITIVE = 1
        ENDELSE
        hide_status_diff = ABS(xd_sState.struc_modifyobject.hide_status[current_step_index, *]-xd_sState.struc_modifyobject.hide_status[current_step_index+1, *])
;        print, 'current hide_status = ', xd_sState.struc_modifyobject.hide_status[current_step_index, *]
;        print, 'old hide_status = ', xd_sState.struc_modifyobject.hide_status[current_step_index+1, *]
        max_n_Poly = MAX(xd_sState.struc_modifyobject.exist_num);define the maximum total number of objects
;        print, 'max_n_Poly = ', max_n_Poly
        
        ;setting the hide status of each oPoly
        FOR i=0, max_n_Poly-1 DO BEGIN
          IF hide_status_diff[0,i] EQ 1 THEN BEGIN
            xd_sState.opoly[i]->getProperty, hide = oldhide
            xd_sState.opoly[i]->setProperty, hide = 1-oldhide
          ENDIF
        ENDFOR
        
        
        ;updating the hiding status of the displayed object
        WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
        IF hide_status_diff[0, ObjectIndex_value-1] EQ 1 THEN BEGIN
          xd_sState.opoly[ObjectIndex_value-1]->getProperty, hide = newhide
          IF newhide EQ 1 THEN BEGIN
            ;disable these options
            widget_control, xd_sState.wColorIndex, sensitive=0
            widget_control, xd_sState.wSmooth, sensitive=0
            widget_control, xd_sState.wSmooth_text, sensitive=0
            widget_control, xd_sState.wtransSlider, sensitive=0
            widget_control, xd_sState.wTrans_text, sensitive=0
            widget_control, xd_sState.w3DSegd, sensitive=0
            widget_control, xd_sState.w3Drg, sensitive=0
            widget_control, xd_sState.wthickness, sensitive=0
            widget_control, xd_sState.wpoint, sensitive=0
            widget_control, xd_sState.wwire, sensitive=0
            widget_control, xd_sState.wfill, sensitive=0
;            widget_control, xd_sState.wexpstl, sensitive=0
            ;            widget_control, xd_sState.whide, sensitive=0
            widget_control, xd_sState.wsinglecenter, sensitive=0

            ;change the icon to hidecancel
            ;              WIDGET_CONTROL, xd_sState.wHide, /BITMAP, SET_VALUE=xd_sState.icon_dir+'singlehidecancel.bmp'

            ;plot a cross to label hiding
            WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
            DEVICE, DECOMPOSED=1
            WSET, wColorIndexID
            plots, [0,40], [0,30], color=[0,0,0], /device
            plots, [0,40], [30,0], color=[0,0,0], /device
          ENDIF ELSE BEGIN
            ;enable these options
            widget_control, xd_sState.wColorIndex, sensitive=1
            widget_control, xd_sState.wSmooth, sensitive=1
            widget_control, xd_sState.wSmooth_text, sensitive=1
            widget_control, xd_sState.wtransSlider, sensitive=1
            widget_control, xd_sState.wTrans_text, sensitive=1
            widget_control, xd_sState.w3DSegd, sensitive=1
            widget_control, xd_sState.w3Drg, sensitive=1
            widget_control, xd_sState.wthickness, sensitive=1
            widget_control, xd_sState.wpoint, sensitive=1
            widget_control, xd_sState.wwire, sensitive=1
            widget_control, xd_sState.wfill, sensitive=1
            widget_control, xd_sState.wexpstl, sensitive=1
            ;            widget_control, xd_sState.whide, sensitive=0
            widget_control, xd_sState.wsinglecenter, sensitive=1
            ;change the icon to hidecancel
            ;              WIDGET_CONTROL, xd_sState.wHide, /BITMAP, SET_VALUE=xd_sState.icon_dir+'singlehide.bmp'

            ;remove the cross to label cancel hiding
            color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
            WIDGET_CONTROL, xd_sState.wColorIndex, GET_VALUE=ColorIndex_value
            WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
            DEVICE, DECOMPOSED=1
            WSET, wColorIndexID
            ERASE, ColorIndexBGR[ColorIndex_value]
          ENDELSE
        ENDIF
        
        ;hide the total volume
        xd_sState.tvol_white.SetProperty,hide=1
        
        ;replotting whole stuff
        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        xd_sState.oWindow->Draw, xd_sState.oView
        ;draw_3DView, xd_sState

        ;replot the 3 draw windows
        ortho2D_plot, xd_sState
      ENDIF
      
      
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

    end
    
    'REDO': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS

      IF xd_sState.struc_modifyobject.current_step_index GE 0 AND xd_sState.struc_modifyobject.current_step_index LE xd_sState.struc_modifyobject.modify_step-1 THEN BEGIN
        xd_sState.struc_modifyobject.current_step_index += 1
        current_step_index = xd_sState.struc_modifyobject.current_step_index
        IF current_step_index GE xd_sState.struc_modifyobject.modify_step THEN BEGIN
          widget_control, xd_sState.wRedo, SENSITIVE = 0
          widget_control, xd_sState.wUndo, SENSITIVE = 1
        ENDIF ELSE BEGIN
          widget_control, xd_sState.wRedo, SENSITIVE = 1
          widget_control, xd_sState.wUndo, SENSITIVE = 1
        ENDELSE
        hide_status_diff = ABS(xd_sState.struc_modifyobject.hide_status[current_step_index, *]-xd_sState.struc_modifyobject.hide_status[current_step_index-1, *])
;        print, 'current hide_status = ', xd_sState.struc_modifyobject.hide_status[current_step_index, *]
;        print, 'old hide_status = ', xd_sState.struc_modifyobject.hide_status[current_step_index-1, *]
        max_n_Poly = MAX(xd_sState.struc_modifyobject.exist_num);define the maximum total number of objects
;        print, 'max_n_Poly = ', max_n_Poly
        ;setting the hide status of each oPoly
        FOR i=0, max_n_Poly-1 DO BEGIN
          IF hide_status_diff[0,i] EQ 1 THEN BEGIN
            xd_sState.opoly[i]->getProperty, hide = oldhide
            xd_sState.opoly[i]->setProperty, hide = 1-oldhide
          ENDIF
        ENDFOR
        
        ;updating the hiding status of the displayed object
        WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
        IF hide_status_diff[0, ObjectIndex_value-1] EQ 1 THEN BEGIN
          xd_sState.opoly[ObjectIndex_value-1]->getProperty, hide = newhide
          IF newhide EQ 1 THEN BEGIN
            ;disable these options
            widget_control, xd_sState.wColorIndex, sensitive=0
            widget_control, xd_sState.wSmooth, sensitive=0
            widget_control, xd_sState.wSmooth_text, sensitive=0
            widget_control, xd_sState.wtransSlider, sensitive=0
            widget_control, xd_sState.wTrans_text, sensitive=0
            widget_control, xd_sState.w3DSegd, sensitive=0
            widget_control, xd_sState.w3Drg, sensitive=0
            widget_control, xd_sState.wthickness, sensitive=0
            widget_control, xd_sState.wpoint, sensitive=0
            widget_control, xd_sState.wwire, sensitive=0
            widget_control, xd_sState.wfill, sensitive=0
;            widget_control, xd_sState.wexpstl, sensitive=0
            ;            widget_control, xd_sState.whide, sensitive=0
            widget_control, xd_sState.wsinglecenter, sensitive=0

            ;change the icon to hidecancel
            ;              WIDGET_CONTROL, xd_sState.wHide, /BITMAP, SET_VALUE=xd_sState.icon_dir+'singlehidecancel.bmp'

            ;plot a cross to label hiding
            WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
            DEVICE, DECOMPOSED=1
            WSET, wColorIndexID
            plots, [0,40], [0,30], color=[0,0,0], /device
            plots, [0,40], [30,0], color=[0,0,0], /device
          ENDIF ELSE BEGIN
            ;enable these options
            widget_control, xd_sState.wColorIndex, sensitive=1
            widget_control, xd_sState.wSmooth, sensitive=1
            widget_control, xd_sState.wSmooth_text, sensitive=1
            widget_control, xd_sState.wtransSlider, sensitive=1
            widget_control, xd_sState.wTrans_text, sensitive=1
            widget_control, xd_sState.w3DSegd, sensitive=1
            widget_control, xd_sState.w3Drg, sensitive=1
            widget_control, xd_sState.wthickness, sensitive=1
            widget_control, xd_sState.wpoint, sensitive=1
            widget_control, xd_sState.wwire, sensitive=1
            widget_control, xd_sState.wfill, sensitive=1
            widget_control, xd_sState.wexpstl, sensitive=1
            ;            widget_control, xd_sState.whide, sensitive=0
            widget_control, xd_sState.wsinglecenter, sensitive=1
            ;change the icon to hidecancel
            ;              WIDGET_CONTROL, xd_sState.wHide, /BITMAP, SET_VALUE=xd_sState.icon_dir+'singlehide.bmp'

            ;remove the cross to label cancel hiding
            color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
            WIDGET_CONTROL, xd_sState.wColorIndex, GET_VALUE=ColorIndex_value
            WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
            DEVICE, DECOMPOSED=1
            WSET, wColorIndexID
            ERASE, ColorIndexBGR[ColorIndex_value]
          ENDELSE
        ENDIF
        
        ;hide the total volume
        xd_sState.tvol_white.SetProperty,hide=1
        
        ;replotting whole stuff        
        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        xd_sState.oWindow->Draw, xd_sState.oView
        ;draw_3DView, xd_sState

        ;replot the 3 draw windows
        ortho2D_plot, xd_sState
      ENDIF

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end

    
    'wrstl': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      outputFile = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/', $
        TITLE='Select stl File', FILTER='*.stl')
      ;      outputfile= 'C:\Users\vit\4.stl'
      IF STRCMP(STRMID(outputFile, STRLEN(outputFile)-4, 4), '.stl') NE 1 THEN outputFile = outputFile + '.stl'
      
      IF STRLEN(outputFile) GT 0 THEN BEGIN
        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
        widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
        widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
        widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1

        ;      vol_HU_cube = xd_sState.vol_HU_cube
        ;      vol_HU_cube_small = vol_HU_cube[ xvalue:xvalue1,yvalue: yvalue1, zvalue: zvalue1]
        ;
        ;      SHADE_VOLUME, vol_HU_cube_small,xd_sState.initThresholdValue_min, Outverts, Outconn1
        ;      smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1);, ITERATIONS=200)   ;smooth the segmented volume; GJ 2017/12/22

        ;;;;;20180417;;;;;;;
        xd_sState.oPoly[0]->getProperty, POLYGONS = Outconn1
        xd_sState.oPoly[0]->getProperty, DATA = smoothedOutverts
        ;;;;;20180417;;;;;;;

        numberVertices = MESH_DECIMATE(smoothedOutverts, Outconn1, p2, VERTICES = v2)



        ;       Out=ptrArr(n_elements(myindex_min))
        ;       for i=0,n_elements(myindex_min)-1 do begin
        ;         threshold_min_value = myindex_min(i)
        ;

        ;      xd_sstate.opoly[0]->getProperty, DATA = vers
        ;      xd_sstate.opoly[0]->getProperty, POLYGONS= conn
        print, 'vol_HU_cube_resolution = ', xd_sState.vol_HU_cube_resolution
        writestl,outputfile,v2*xd_sState.vol_HU_cube_resolution,p2
      ENDIF
      
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end
     '3DSEGD': begin
       WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
       
       ;whether selected is true
       IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
         WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
         RETURN
       ENDIF

       IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
         WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
         RETURN
       ENDIF

       ;make sure verts and object are selected.
       xd_sState.selected->getProperty, data=verts
       IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
         WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
         RETURN
       ENDIF

       ;Do screen-capture of draw3 window
       xd_sState.oPoly_border->SetProperty, STYLE = 2
       xd_sState.oWindow->Draw, xd_sState.oView
       pick_border = xd_sState.oWindow->PickData(xd_sState.oView, $
         xd_sState.oPoly_border, [xd_sState.draw3_size/2., xd_sState.draw3_size/2.], dataxyz_border)
       dist_cent = 3; xd_sState.draw3_size/10.

       if (pick_border ne 0) then begin
         print, 'border dataxyz = ', dataxyz_border
         IF ABS(dataxyz_border[0]-xd_sState.oPoly_min[0]) LT dist_cent OR ABS(dataxyz_border[0]-xd_sState.oPoly_max[0]) LT dist_cent OR $
           ABS(dataxyz_border[0]-xd_sState.oPoly_min[0]+xd_sState.oPoly_size[0]) LT dist_cent OR ABS(dataxyz_border[0]-xd_sState.oPoly_max[0]-xd_sState.oPoly_size[0]) LT dist_cent THEN BEGIN
           print, 'on yz plane'
           xd_sState.direction = 'yz'
           verts_border = [[MIN(verts[0,*]), MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*]), MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*]), MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MIN(verts[0,*]), MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MIN(verts[0,*]), MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*]), MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*]), MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]], $
             [MIN(verts[0,*]), MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]]]
           xd_sState.oPoly_border->setProperty, DATA = verts_border
         ENDIF
         IF ABS(dataxyz_border[1]-xd_sState.oPoly_min[1]) LT dist_cent OR ABS(dataxyz_border[1]-xd_sState.oPoly_max[1]) LT dist_cent OR $
           ABS(dataxyz_border[1]-xd_sState.oPoly_min[1]+xd_sState.oPoly_size[1]) LT dist_cent OR ABS(dataxyz_border[1]-xd_sState.oPoly_max[1]-xd_sState.oPoly_size[1]) LT dist_cent THEN BEGIN
           print, 'on xz plane'
           xd_sState.direction = 'xz'
           verts_border = [[MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MIN(verts[1,*]), MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MIN(verts[1,*]), MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MAX(verts[1,*]), MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MAX(verts[1,*]), MIN(verts[2,*])-2.*xd_sState.oPoly_size[2]], $
             [MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MIN(verts[1,*]), MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MIN(verts[1,*]), MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MAX(verts[1,*]), MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]], $
             [MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MAX(verts[1,*]), MAX(verts[2,*])+2.*xd_sState.oPoly_size[2]]]
           xd_sState.oPoly_border->setProperty, DATA = verts_border
         ENDIF
         IF ABS(dataxyz_border[2]-xd_sState.oPoly_min[2]) LT dist_cent OR ABS(dataxyz_border[2]-xd_sState.oPoly_max[2]) LT dist_cent OR $
           ABS(dataxyz_border[2]-xd_sState.oPoly_min[2]+xd_sState.oPoly_size[2]) LT dist_cent OR ABS(dataxyz_border[2]-xd_sState.oPoly_max[2]-xd_sState.oPoly_size[2]) LT dist_cent THEN BEGIN
           print, 'on xy plane'
           xd_sState.direction = 'xy'
           verts_border = [[MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MIN(verts[2,*])], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MIN(verts[2,*])], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MIN(verts[2,*])], $
             [MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MIN(verts[2,*])], $
             [MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MAX(verts[2,*])], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MIN(verts[1,*])-2.*xd_sState.oPoly_size[1], MAX(verts[2,*])], $
             [MAX(verts[0,*])+2.*xd_sState.oPoly_size[0], MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MAX(verts[2,*])], $
             [MIN(verts[0,*])-2.*xd_sState.oPoly_size[0], MAX(verts[1,*])+2.*xd_sState.oPoly_size[1], MAX(verts[2,*])]]
           xd_sState.oPoly_border->setProperty, DATA = verts_border
         ENDIF
         widget_control, /hourglass
         xd_sState.oWindow->Draw, xd_sState.oView
 ;        ;draw_3DView, xd_sState
       endif
       xd_sState.oPoly_border->SetProperty, STYLE = 0
       xd_sState.oWindow->Draw, xd_sState.oView
       
       odraw_img = xd_sState.oWindow->READ()
       odraw_img->getProperty, DATA = draw_img
       
       xd_sState.oPoly_border->SetProperty, STYLE = 2
       xd_sState.oWindow->Draw, xd_sState.oView
       ;      print, 'draw img prop: ', SIZE(draw_img, /DIMENSIONS)
       ;      WRITE_BMP, 'test.bmp', draw_img
       XROI_01, draw_img, TITLE = 'ROI selection', REGIONS_OUT = seg_ROIout, X_SCROLL_SIZE=xd_sState.draw3_size, Y_SCROLL_SIZE=xd_sState.draw3_size, /BLOCK
       ;;previous ROI is drawn
       vol_dims = SIZE(xd_sState.vol_HU_cube, /DIMENSIONS)
       screenSize = xd_sState.draw3_size
       IF OBJ_VALID(seg_ROIout[N_ELEMENTS(seg_ROIout)-1]) THEN BEGIN
         seg_ROIout[N_ELEMENTS(seg_ROIout)-1]->GetProperty, DATA=roiData, N_Verts = nverts, ROI_XRANGE=xr, ROI_YRANGE=yr
         print, 'roi data = ', roiData
         roi_3D_data_a = DBLARR(2, nverts)
         location_data_a = DBLARR(nverts)
         FOR j=0, nverts-1 DO BEGIN
           pick_border = xd_sState.oWindow->PickData(xd_sState.oView, $
             xd_sState.oPoly_border, roiData[0:1, j], dataxyz_border)
           print, 'dataxyz_border = ', dataxyz_border
           IF STRCMP(xd_sState.direction, 'xy') THEN BEGIN
             roi_3D_data_a[*,j] = [dataxyz_border[0], dataxyz_border[1]]
             location_data_a[j] = dataxyz_border[2]
           ENDIF
           IF STRCMP(xd_sState.direction, 'yz') THEN BEGIN
             roi_3D_data_a[*,j] = [dataxyz_border[1], dataxyz_border[2]]
             location_data_a[j] = dataxyz_border[0]
           ENDIF
           IF STRCMP(xd_sState.direction, 'xz') THEN BEGIN
             roi_3D_data_a[*,j] = [dataxyz_border[0], dataxyz_border[2]]
             location_data_a[j] = dataxyz_border[1]
           ENDIF
         ENDFOR

         xd_sState.selected->getProperty, DATA = verts
         xd_sState.oPoly_border->getProperty, DATA = verts_border
         temp_verts = verts
         temp_verts_border = verts_border
         location_a = MEDIAN(location_data_a)
         IF STRCMP(xd_sState.direction, 'xy') THEN BEGIN
           IF ABS(location_a-MIN(verts[2,*])) LT ABS(location_a-MAX(verts[2,*])) THEN BEGIN
             temp_verts[2,*] = verts[2,*] + (MAX(verts[2,*]) - MIN(verts[2,*]))
             temp_verts_border[2,*] = verts_border[2,*] + (MAX(verts[2,*]) - MIN(verts[2,*]))
           ENDIF ELSE BEGIN
             temp_verts[2,*] = verts[2,*] - (MAX(verts[2,*]) - MIN(verts[2,*]))
             temp_verts_border[2,*] = verts_border[2,*] - (MAX(verts[2,*]) - MIN(verts[2,*]))
           ENDELSE
         ENDIF
         IF STRCMP(xd_sState.direction, 'yz') THEN BEGIN
           IF ABS(location_a-MIN(verts[0,*])) LT ABS(location_a-MAX(verts[0,*])) THEN BEGIN
             temp_verts[0,*] = verts[0,*] + (MAX(verts[0,*]) - MIN(verts[0,*]))
             temp_verts_border[0,*] = verts_border[0,*] + (MAX(verts[0,*]) - MIN(verts[0,*]))
           ENDIF ELSE BEGIN
             temp_verts[0,*] = verts[0,*] - (MAX(verts[0,*]) - MIN(verts[0,*]))
             temp_verts_border[0,*] = verts_border[0,*] - (MAX(verts[0,*]) - MIN(verts[0,*]))
           ENDELSE
         ENDIF
         IF STRCMP(xd_sState.direction, 'xz') THEN BEGIN
           IF ABS(location_a-MIN(verts[1,*])) LT ABS(location_a-MAX(verts[1,*])) THEN BEGIN
             temp_verts[1,*] = verts[1,*] + (MAX(verts[1,*]) - MIN(verts[1,*]))
             temp_verts_border[1,*] = verts_border[1,*] + (MAX(verts[1,*]) - MIN(verts[1,*]))
           ENDIF ELSE BEGIN
             temp_verts[1,*] = verts[1,*] - (MAX(verts[1,*]) - MIN(verts[1,*]))
             temp_verts_border[1,*] = verts_border[1,*] - (MAX(verts[1,*]) - MIN(verts[1,*]))
           ENDELSE
         ENDIF
         xd_sState.selected->setProperty, DATA = temp_verts
         xd_sState.oPoly_border->setProperty, DATA = temp_verts_border
         xd_sState.oWindow->Draw, xd_sState.oView
         ;plot another contour on the back wall
         roi_3D_data_b = DBLARR(2, nverts)
         location_data_b = DBLARR(nverts)
         FOR j=0, nverts-1 DO BEGIN
           pick_border = xd_sState.oWindow->PickData(xd_sState.oView, $
             xd_sState.oPoly_border, roiData[0:1, j], dataxyz_border)
           print, 'new dataxyz_border = ', dataxyz_border
           IF STRCMP(xd_sState.direction, 'xy') THEN BEGIN
             roi_3D_data_b[*,j] = [dataxyz_border[0], dataxyz_border[1]]
             location_data_b[j] = dataxyz_border[2]
           ENDIF
           IF STRCMP(xd_sState.direction, 'yz') THEN BEGIN
             roi_3D_data_b[*,j] = [dataxyz_border[1], dataxyz_border[2]]
             location_data_b[j] = dataxyz_border[0]
           ENDIF
           IF STRCMP(xd_sState.direction, 'xz') THEN BEGIN
             roi_3D_data_b[*,j] = [dataxyz_border[0], dataxyz_border[2]]
             location_data_b[j] = dataxyz_border[1]
           ENDIF
         ENDFOR

         xd_sState.selected->setProperty, DATA = verts
         xd_sState.oPoly_border->setProperty, DATA = verts_border
         xd_sState.oPoly_border->setProperty, STYLE = 0
         xd_sState.oWindow->Draw, xd_sState.oView
  ;       ;draw_3DView, xd_sState
         location_b = MEDIAN(location_data_b)
         roi_3D_data = DBLARR(2, nverts)

         IF STRCMP(xd_sState.direction, 'xy') THEN BEGIN
           FOR imind=0, vol_dims[2]-1 DO BEGIN
             location = imind
             roi_3D_data[0,*] = (roi_3D_data_b[0,*]-roi_3D_data_a[0,*])/(location_b-location_a)*(location-location_a)+roi_3D_data_a[0,*]
             roi_3D_data[1,*] = (roi_3D_data_b[1,*]-roi_3D_data_a[1,*])/(location_b-location_a)*(location-location_a)+roi_3D_data_a[1,*]
             ROI = OBJ_NEW('IDLgrROI', roi_3D_data)
             maskResult = ROI -> ComputeMask(DIMENSIONS = [vol_dims[0], vol_dims[1]])
             xd_sState.mask_vol_cube_modify[*, *, imind] = (CONGRID(maskResult, vol_dims[0], vol_dims[1]) GT 0); (REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0
             OBJ_DESTROY, ROI
           ENDFOR
         ENDIF


         IF STRCMP(xd_sState.direction, 'yz') THEN BEGIN
           FOR imind=0, vol_dims[0]-1 DO BEGIN
             location = imind
             roi_3D_data[0,*] = (roi_3D_data_b[0,*]-roi_3D_data_a[0,*])/(location_b-location_a)*(location-location_a)+roi_3D_data_a[0,*]
             roi_3D_data[1,*] = (roi_3D_data_b[1,*]-roi_3D_data_a[1,*])/(location_b-location_a)*(location-location_a)+roi_3D_data_a[1,*]
             ROI = OBJ_NEW('IDLgrROI', roi_3D_data)
             maskResult = ROI -> ComputeMask(DIMENSIONS = [vol_dims[1], vol_dims[2]])
             xd_sState.mask_vol_cube_modify[imind, *, *] = (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
             OBJ_DESTROY, ROI
           ENDFOR
         ENDIF

         IF STRCMP(xd_sState.direction, 'xz') THEN BEGIN
           FOR imind=0, vol_dims[1]-1 DO BEGIN
             location = imind
             roi_3D_data[0,*] = (roi_3D_data_b[0,*]-roi_3D_data_a[0,*])/(location_b-location_a)*(location-location_a)+roi_3D_data_a[0,*]
             roi_3D_data[1,*] = (roi_3D_data_b[1,*]-roi_3D_data_a[1,*])/(location_b-location_a)*(location-location_a)+roi_3D_data_a[1,*]
             ROI = OBJ_NEW('IDLgrROI', roi_3D_data)
             maskResult = ROI -> ComputeMask(DIMENSIONS = [vol_dims[0], vol_dims[2]])
             xd_sState.mask_vol_cube_modify[*, imind, *] = (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)
             OBJ_DESTROY, ROI
           ENDFOR
         ENDIF

;         xd_sState.mask_vol_cube_modify = 1-xd_sState.mask_vol_cube_modify
         
         ;define the mask
         WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
         IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
           mask_vol_cube = xd_sState.new_mask_vol_cube
           ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
           ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
         ENDIF ELSE BEGIN
           mask_vol_cube = xd_sState.new_mask_vol_cube*0
           mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
         ENDELSE
         
         modify_roi = 8 ; remove the ROI region
;         Mask_segmentation, xd_sState, mask_vol_cube
         lasso_segmentation_3D, xd_sState, mask_vol_cube
       ENDIF

       WIDGET_CONTROL, sEvent.top, /HOURGLASS
       xd_sState.oPoly_border->SetProperty, STYLE = 0
       xd_sState.oWindow->Draw, xd_sState.oView
       ;draw_3DView, xd_sState

       ;replot the 3 draw windows
       ortho2D_plot, xd_sState

       WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
     end

    '3DRG': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      
       dataxyz = xd_sState.seed_location
       IF ABS(dataxyz[0]) GT 0.1 AND ABS(dataxyz[1]) GT 0.1 AND ABS(dataxyz[2]) GT 0.1 THEN BEGIN
         ;define the mask
         WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
         IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
           mask_vol_cube = xd_sState.new_mask_vol_cube
           ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
           ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
         ENDIF ELSE BEGIN
           mask_vol_cube = xd_sState.new_mask_vol_cube*0
           mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
         ENDELSE
         
         growpixel8, mask_vol_cube, dataxyz[0], dataxyz[1], dataxyz[2], 0.9, b, notb
         sizeb=size(b)
         IF(sizeb[0] NE 0) THEN BEGIN
           xd_sState.mask_vol_cube_modify = (b GT 0.9)
           modify_roi = 20
           lasso_segmentation_3D, xd_sState, mask_vol_cube

           WIDGET_CONTROL, sEvent.top, /HOURGLASS
           xd_sState.oWindow->Draw, xd_sState.oView
           ;draw_3DView, xd_sState
  
           ;replot the 3 draw windows
           ortho2D_plot, xd_sState
         ENDIF
       ENDIF

       WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    ;GJ 2020/3/22
    'AIVesselSeg': begin ;AI-based vessel segmentation
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS

      Add_new_by_AIVesselSegmentation, xd_sState

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of AI-based vessl segmentation
    
    ;vessl wall contrast enhancement segmentation
    'VWCESeg': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      
       dataxyz = xd_sState.seed_location
       IF ABS(dataxyz[0]) GT 0.1 AND ABS(dataxyz[1]) GT 0.1 AND ABS(dataxyz[2]) GT 0.1 THEN BEGIN
         ;define the mask
         WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
         IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
           mask_vol_cube = xd_sState.new_mask_vol_cube
           ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
           ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
         ENDIF ELSE BEGIN
           mask_vol_cube = xd_sState.new_mask_vol_cube*0
           mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
         ENDELSE
         
         
         Add_new_by_VWCESegmentation, xd_sState, mask_vol_cube
         
         WIDGET_CONTROL, sEvent.top, /HOURGLASS
         xd_sState.oWindow->Draw, xd_sState.oView
         ;draw_3DView, xd_sState

         ;replot the 3 draw windows
         ortho2D_plot, xd_sState
       ENDIF

       WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end       ;  of vessl wall contrast enhancement segmentation
    ;GJ 2019/8/16
    ;'Facepic': begin
    'FacePicReg': begin
      ;    '3DRG': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
      IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
        mask_vol_cube = xd_sState.new_mask_vol_cube
        ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
        ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
      ENDIF ELSE BEGIN
        mask_vol_cube = xd_sState.new_mask_vol_cube*0
        mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
;        img = PROJECT_VOL(mask_vol_cube*255., 64, 64, 256, /AVG_INTENSITY)
        img = DBLARR(256, 256)*0.
        FOR i=0, 255 DO BEGIN
          FOR j=0, 255 DO BEGIN
            FOR k = 0, 255 DO BEGIN
              ;extract the face surface
               IF mask_vol_cube[i, k, j] GT 0 THEN break
               img[i,j] = img[i,j] + 1
            ENDFOR
          ENDFOR
        ENDFOR
;        iimage, img
        ;make the image minimum as 0.
        img_temp = BYTSCL(256.-img)
        min_img = MIN(img_temp)
;        print, 'minimum img =', min_img
        img = img_temp - min_img
        index = WHERE(img GT 1, count)
        IF count NE 0 THEN BEGIN
          index2d = ARRAY_INDICES(img, index)
          xdim = max(index2d[0,*]) - min(index2d[0,*]) + 1
          ydim = max(index2d[1,*]) - min(index2d[1,*]) + 1
          xarray = DBLARR(xdim, ydim)
          yarray = DBLARR(xdim, ydim)
          zarray = DBLARR(xdim, ydim)
          FOR i=0, xdim-1 DO BEGIN
            FOR j=0, ydim-1 DO BEGIN
              xarray[i, j] = min(index2d[0,*]) + i
              yarray[i, j] = min(index2d[1,*]) + j
              zarray[i, j] = img[min(index2d[0,*]) + i, min(index2d[1,*]) + j]
            ENDFOR
          ENDFOR
        ENDIF
        
        WRITE_JPEG, 'D:\AIMIS_3D\CBCTfaceKelly1.jpg', BYTSCL(HIST_EQUAL(zarray)) 
        img_ratio = CEIL(MAX([xdim, ydim])/512.)
        IF img_ratio EQ 1 THEN img_ratio = MAX([xdim, ydim])/512.
        WINDOW, 0, XSIZE = xdim/img_ratio, YSIZE = ydim/img_ratio, TITLE = 'CT Face Projection Image: Please mark left eye'
        TV, BYTSCL(HIST_EQUAL(CONGRID(zarray, xdim/img_ratio, ydim/img_ratio)))
        WSET, 0
        CURSOR, xi1, yi1, /DEVICE
        print, 'xi1, yi1 = ', xi1, yi1
        WDELETE, 0
        WAIT, 2 ;wait 2 sec

        WINDOW, 0, XSIZE = xdim/img_ratio, YSIZE = ydim/img_ratio, TITLE = 'CT Face Projection Image: Please mark right eye'
        TV, BYTSCL(HIST_EQUAL(CONGRID(zarray, xdim/img_ratio, ydim/img_ratio)))
        WSET, 0
        CURSOR, xi2, yi2, /DEVICE
        print, 'xi2, yi2 = ', xi2, yi2
        WDELETE, 0
        WAIT, 2 ;wait 2 sec
        
        eye1_CT = [xi1*img_ratio, yi1*img_ratio]
        eye2_CT = [xi2*img_ratio, yi2*img_ratio]
        eyemiddle_CT = (eye1_CT + eye2_CT)/2.
        eyedistance_CT = SQRT((eye1_CT[0]-eye2_CT[0])^2 + (eye1_CT[1]-eye2_CT[1])^2)
;        xroi, zarray
                        
        ;imagefile='D:\AIMIS_3D\kelly2.jpg';lucy2.jpg';kelly.jpg';lucy2.jpg';kelly.jpg'; gj.jpg
        imagefile = DIALOG_PICKFILE(PATH='D:/AIMIS_3D/', $
          TITLE='Select face photo File', FILTER='*.jpg', /MUST_EXIST)
        ;      inputfile = 'E:\Research\MED_visualization\6_spineCases\Li  guizhen\LI GUIZHEN.stl'
        IF STRLEN(imagefile) GT 0 THEN BEGIN
          read_jpeg, imagefile, idata1
          image=idata1
          image_size = SIZE(idata1, /dim)
          image_ratio = CEIL(MAX(image_size)/512.)
          IF image_ratio EQ 1 THEN image_ratio = MAX(image_size)/512.
          WINDOW, 0, XSIZE = image_size[1]/image_ratio, YSIZE = image_size[2]/image_ratio, TITLE = 'Face Photo Image: Please mark left eye'
          TV, CONGRID(idata1, 3, image_size[1]/image_ratio, image_size[2]/image_ratio), TRUE=1
          WSET, 0
          CURSOR, xi10, yi10, /DEVICE
          WDELETE, 0
          WAIT, 2 ;wait 2 sec

          WINDOW, 0, XSIZE = image_size[1]/image_ratio, YSIZE = image_size[2]/image_ratio, TITLE = 'Face Photo Image: Please mark right eye'
          TV, CONGRID(idata1, 3, image_size[1]/image_ratio, image_size[2]/image_ratio), TRUE=1
          WSET, 0
          CURSOR, xi20, yi20, /DEVICE
          WDELETE, 0

          ;GJ 2019/8/18
          ;semi-automatically adjust the images
          eye1 = [xi10*image_ratio, yi10*image_ratio]
          eye2 = [xi20*image_ratio, yi20*image_ratio]
          eyemiddle = (eye1 + eye2)/2.
          eyedistance = SQRT((eye1[0]-eye2[0])^2 + (eye1[1]-eye2[1])^2)
          ratio = eyedistance/eyedistance_CT
          xmin = (eyemiddle[0]-ratio*eyemiddle_CT[0])
          IF xmin LT 0 THEN xmin = 0
          xmax = (eyemiddle[0]+ratio*(xdim-eyemiddle_CT[0]))
          IF xmax GT image_size[1] THEN xmax = image_size[1]-1
          ymin = (eyemiddle[1]-ratio*eyemiddle_CT[1])
          IF ymin LT 0 THEN ymin = 0
          ymax = (eyemiddle[1]+ratio*(ydim-eyemiddle_CT[1]))
          IF ymax GT image_size[2] THEN ymax = image_size[2]-1
          new_image = image[*, xmin:xmax, ymin:ymax]

          xd_sState.oImage1->SetProperty, DATA=new_image, INTERLEAVE=0
          xd_sState.oSurface->SetProperty, DATAZ=GAUSS_SMOOTH(zarray, 2), DATAX=xarray, DATAY=yarray, TEXTURE_MAP=xd_sState.oImage1, HIDE=0   
          ;smooth zarray using gauss_smooth
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState

          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

          ;        ;GJ, 2019/8/16 displaying the result
          ;        d_flythru, SURFACE=img, image=image
        ENDIF
        
      ENDELSE
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'SpineSegmentation1': begin
;    '3DRG': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      dataxyz = xd_sState.seed_location
      IF ABS(dataxyz[0]) GT 0.1 AND ABS(dataxyz[1]) GT 0.1 AND ABS(dataxyz[2]) GT 0.1 THEN BEGIN
        ;define the mask
        WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
        IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
          mask_vol_cube = xd_sState.new_mask_vol_cube
          ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
          ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
        ENDIF ELSE BEGIN
          mask_vol_cube = xd_sState.new_mask_vol_cube*0
          mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
        ENDELSE

        Add_new_by_region_grow_spine, xd_sState, mask_vol_cube

        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        xd_sState.oWindow->Draw, xd_sState.oView
        ;draw_3DView, xd_sState

        ;replot the 3 draw windows
        ortho2D_plot, xd_sState

        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
      ENDIF

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end


    'hfslicexslider': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ;  (*pState).rVolume->SetProperty, DATA0 = (*pState).originalvolumedata
      ;   widget_control, (*pState).hSlice_x_slider, get_value = xvalue
      ;         widget_control, (*pState).hSlice_y_slider, get_value = yvalue
      ;         widget_control, xd_sState.hSlice_z_slider, get_value = zvalue


      if  (xd_sState.flag_x eq 1)then begin

        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
        widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
        widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
        widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1
        
        xd_sState.opoly[0]->SetProperty, $
          ; CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)]]
          CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]
       endif else begin


       if ( xd_sState.hSlice_on_status  eq  1)then begin
;          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
;          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
;          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
;
;          xd_sstate.opoly[0]->SetProperty, $
;            CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
;        endif  else begin
          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue

          xd_sState.opoly[0]->SetProperty, $
            CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]

          if xd_sState.flag_other then begin
            xd_sState.opoly1->SetProperty, $
              CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
            xd_sState.opoly2->SetProperty, $
              CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
          endif
        endif

      endelse


      ;mask_segmentation,xd_sState
      ;updatepoly,xd_sState
      ;xd_sState.orotationmodel->Add,xd_sState.opoly
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end

    'hfslicexslider1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ;  xd_sState.rVolume->SetProperty, DATA0 = xd_sState.originalvolumedata
      ;   widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;         widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;         widget_control, xd_sState.hSlice_z_slider, get_value = zvalue


      if  (xd_sState.flag_x eq 1)then begin

        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
        widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
        widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
        widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1

        xd_sState.opoly[0]->SetProperty, $
          ; CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)]]
          CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]
      endif
      ;updatepoly,xd_sState
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end



    'hfsliceyslider': begin

      ;            xd_sState.rVolume->SetProperty, DATA0 = xd_sState.originalvolumedata
      ;         widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;   widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;         widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      if  (xd_sState.flag_x eq 1)then begin

        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
        widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
        widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
        widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1
        
        xd_sState.opoly[0]->SetProperty, $
          ; CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)]]
          CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]
      endif else begin

        if ( xd_sState.hSlice_on_status  eq  1)then begin
;          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
;          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
;          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
;
;          xd_sstate.opoly[0]->SetProperty, $
;            CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
;        endif  else begin
          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue

          xd_sState.opoly[0]->SetProperty, $
            CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]


          if xd_sState.flag_other then begin
            xd_sState.opoly1->SetProperty, $
              CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
            xd_sState.opoly2->SetProperty, $
              CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
          endif
        endif
      endelse
      ;xd_sState.rWindow->Draw, xd_sState.rView
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end

    'hfsliceyslider1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ;  xd_sState.rVolume->SetProperty, DATA0 = xd_sState.originalvolumedata
      ;   widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;         widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;         widget_control, xd_sState.hSlice_z_slider, get_value = zvalue


      if  (xd_sState.flag_x eq 1)then begin

      widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
      widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
      widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1

     xd_sState.opoly[0]->SetProperty, $
        ; CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)]]
        CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]
endif

      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end





    'hfslicezslider': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ;            xd_sState.rVolume->SetProperty, DATA0 = xd_sState.originalvolumedata
      ;         widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;         widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;         widget_control, xd_sState.hSlice_z_slider, get_value = zvalue



      if  (xd_sState.flag_x eq 1)then begin

        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
        widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
        widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
        widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1

         xd_sState.opoly[0]->SetProperty, $
          ; CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)]]
          CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]

      endif else begin

        if ( xd_sState.hSlice_on_status  eq  1)then begin
;          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
;          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
;          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
;
;          xd_sstate.opoly[0]->SetProperty, $
;            CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
;        endif  else begin
          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue

          xd_sState.opoly[0]->SetProperty, $
            CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]

          if xd_sState.flag_other then begin
            xd_sState.opoly1->SetProperty, $
              CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
            xd_sState.opoly2->SetProperty, $
              CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
          endif
        endif

      endelse
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY

    end

    'hfslicezslider1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ;  xd_sState.rVolume->SetProperty, DATA0 = xd_sState.originalvolumedata
      ;   widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;         widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;         widget_control, xd_sState.hSlice_z_slider, get_value = zvalue


      if  (xd_sState.flag_x eq 1)then begin

      widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
      widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
      widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1

       xd_sState.opoly[0]->SetProperty, $
        ; CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)]]
        CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]
 
      endif  
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end






;    'hfsliceon': begin
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;      xd_sState.hSlice_on_status =  sEvent .select
;
;      if xd_sState.hSlice_on_status then begin
;
;
;
;        widget_control, xd_sState.hSlice_x_slider, SENSITIVE = 1
;        widget_control, xd_sState.hSlice_y_slider, SENSITIVE = 1
;        widget_control, xd_sState.hSlice_z_slider, SENSITIVE = 1
;        widget_control, xd_sState.hSlice_on_x_button, SENSITIVE = 1
;        widget_control, xd_sState.hSlice_on_y_button, SENSITIVE = 1
;        widget_control, xd_sState.hSlice_on_z_button, SENSITIVE = 1
;        widget_control, xd_sState.hSlice_save,     SENSITIVE = 1
;        widget_control, xd_sState.hSlice_three,     SENSITIVE = 1
;        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
;        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
;        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
;        widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
;        widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
;        widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1
;
;        if(xd_sState.flag_x eq 1) then begin
;          widget_control, xd_sState.hSlice_x_slider1, SENSITIVE = 1
;          widget_control, xd_sState.hSlice_y_slider1, SENSITIVE = 1
;          widget_control, xd_sState.hSlice_z_slider1, SENSITIVE = 1
;          widget_control, xd_sState.hSlice_save,     SENSITIVE = 0
;
;        endif
;
;        xd_sstate.opoly[0]->SetProperty, $
;          CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
;
;        ;xd_sState.oWindow->Draw, xd_sState.oView
;
;      endif else begin
;
;        widget_control, xd_sState.hSlice_x_slider,    SENSITIVE = 0
;        widget_control, xd_sState.hSlice_y_slider,    SENSITIVE = 0
;        widget_control, xd_sState.hSlice_z_slider,    SENSITIVE = 0
;        widget_control, xd_sState.hSlice_x_slider1,    SENSITIVE = 0
;        widget_control, xd_sState.hSlice_y_slider1,    SENSITIVE = 0
;        widget_control, xd_sState.hSlice_z_slider1,    SENSITIVE = 0
;        widget_control, xd_sState.hSlice_on_x_button, SENSITIVE = 0
;        widget_control, xd_sState.hSlice_on_y_button, SENSITIVE = 0
;        widget_control, xd_sState.hSlice_on_z_button, SENSITIVE = 0
;        widget_control, xd_sState.hSlice_save,        SENSITIVE = 0
;        widget_control, xd_sState.hSlice_three,        SENSITIVE = 0
;        widget_control, xd_sState.hSlice_on_x_button, SET_BUTTON = 0
;        widget_control, xd_sState.hSlice_on_y_button, SET_BUTTON = 0
;        widget_control, xd_sState.hSlice_on_z_button, SET_BUTTON = 0
;
;        xd_sState.hSlice_on_x_button_status = 0
;        xd_sState.hSlice_on_x_button_status = 0
;        xd_sState.hSlice_on_x_button_status = 0
;        xd_sState.flag_select=0
;        xd_sState.flag_other =0
;        ;xd_sState. flag_x=0
;        ;          widget_control, xd_sState.hSlice_x_slider, set_value = 1
;        ;          widget_control, xd_sState.hSlice_y_slider, set_value = 1
;        ;          widget_control, xd_sState.hSlice_z_slider, set_value = 1
;        ;          widget_control, xd_sState.hSlice_x_slider1, set_value = 1
;        ;          widget_control, xd_sState.hSlice_y_slider1, set_value = 1
;        ;          widget_control, xd_sState.hSlice_z_slider1, set_value = 1
;
;        ;             tmpdd = intarr(xd_sState.volumexlength,xd_sState.volumeylength,xd_sState.volumezlength)
;        ;                tmpdd = xd_sState.originalvolumedata
;        ;                xd_sState.rVolume->SetProperty,DATA = tmpdd
;
;        xd_sstate.opoly[0]->SetProperty, CLIP_PLANES = 0
;        xd_sState.opoly1->SetProperty, CLIP_PLANES=[1, 1, 1, -3]
;        ;xd_sState.rWindow->Draw, xd_sState.rView
;
;      endelse
;      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
;    end

'hfsliceon': begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
        xd_sState.hSlice_on_status =  sEvent .select
       
        if xd_sState.hSlice_on_status then begin
         

  
  xd_sState.flag_other=1
          widget_control, xd_sState.hSlice_x_slider, SENSITIVE = 1
          widget_control, xd_sState.hSlice_y_slider, SENSITIVE = 1
          widget_control, xd_sState.hSlice_z_slider, SENSITIVE = 1
          widget_control, xd_sState.hSlice_x_slider1, SENSITIVE = 0
          widget_control, xd_sState.hSlice_y_slider1, SENSITIVE =0
          widget_control, xd_sState.hSlice_z_slider1, SENSITIVE =0
          
          widget_control, xd_sState.hSlice_on_x_button, SENSITIVE = 1
          widget_control, xd_sState.hSlice_on_y_button, SENSITIVE = 1
         ; widget_control, xd_sState.hSlice_on_z_button, SENSITIVE = 1
          widget_control, xd_sState. hSlice_output,     SENSITIVE = 0
          widget_control, xd_sState.hSlice_three,     SENSITIVE = 1
          widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
          widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
          widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
          
          
         xd_sState.opoly1->SetProperty, $
         CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
         xd_sState.opoly2->SetProperty, $
         CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
         xd_sstate.opoly[0]->SetProperty, $
         CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
        
        endif
        xd_sState.oWindow->Draw, xd_sState.oView
        ortho2D_plot, xd_sState
        ;draw_3DView, xd_sState
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
      end

      'hthree': begin
        WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
        xd_sState.flag_x = sevent.select
       
        if xd_sState.flag_x then begin
           xd_sState.flag_other=0
           xd_sState.opoly1->SetProperty, $
                      CLIP_PLANES=[1, 1, 1, -3]
          
                  xd_sState.flag_select =0
                    
                      xd_sstate.opoly[0]->SetProperty, $
                        ALPHA_CHANNEL=1
          
                      
         widget_control, xd_sState.hSlice_x_slider, SENSITIVE = 1
          widget_control, xd_sState.hSlice_y_slider, SENSITIVE = 1
          widget_control, xd_sState.hSlice_z_slider, SENSITIVE = 1
          widget_control, xd_sState.hSlice_x_slider1, SENSITIVE = 1
          widget_control, xd_sState.hSlice_y_slider1, SENSITIVE = 1
          widget_control, xd_sState.hSlice_z_slider1, SENSITIVE = 1
          
          widget_control, xd_sState.hSlice_on_x_button, SENSITIVE = 1
          widget_control, xd_sState.hSlice_on_y_button, SENSITIVE = 0
        ;  widget_control, xd_sState.hSlice_on_z_button, SENSITIVE = 1
         widget_control, xd_sState. hSlice_output,     SENSITIVE = 1
          widget_control, xd_sState.hSlice_three,     SENSITIVE = 1
        
          


        endif
        xd_sState.oWindow->Draw, xd_sState.oView
        ortho2D_plot, xd_sState
        ;draw_3DView, xd_sState
        WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
      end

      ;    'hthree': begin
      ;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      ;      xd_sState.flag_x = sevent.select
      ;
      ;      if xd_sState.flag_x then begin
      ;        widget_control, xd_sState.hSlice_save,     SENSITIVE = 0
      ;        widget_control, xd_sState.hSlice_x_slider1, SENSITIVE = 1
      ;        widget_control, xd_sState.hSlice_y_slider1, SENSITIVE = 1
      ;        widget_control, xd_sState.hSlice_z_slider1, SENSITIVE = 1
      ;      endif else begin
      ;
      ;        widget_control, xd_sState.hSlice_save,     SENSITIVE =1
      ;        widget_control, xd_sState.hSlice_x_slider1, SENSITIVE = 0
      ;        widget_control, xd_sState.hSlice_y_slider1, SENSITIVE = 0
      ;        widget_control, xd_sState.hSlice_z_slider1, SENSITIVE = 0
      ;
      ;
      ;      endelse
      ;
      ;      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
      ;    end
    'hfsliceonxbutton': begin

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      vol_HU_cube_size=size(xd_sState.vol_HU_cube)
      widget_control, xd_sState.hSlice_x_slider, set_value =1
      widget_control, xd_sState.hSlice_y_slider, set_value =1
      widget_control, xd_sState.hSlice_z_slider, set_value =1

      widget_control, xd_sState.hSlice_x_slider1, set_value = vol_HU_cube_size[1]-1
      widget_control, xd_sState.hSlice_y_slider1, set_value = vol_HU_cube_size[2]-1
      widget_control, xd_sState.hSlice_z_slider1, set_value = vol_HU_cube_size[3]-1
      
      widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
      widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
      widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1
      
      widget_control, xd_sState.wThresholdSlider_min, set_value = xd_sState.initThresholdValue_min
      widget_control, xd_sState.wThresholdSlider_max, set_value = xd_sState.initThresholdValue_max
      
      xd_sstate.opoly[0]->SetProperty, $
        CLIP_PLANES=[[-xvalue, -1, -1, (xvalue*xvalue + 2)],[-1,-yvalue,-1, (yvalue*yvalue + 2)],[-1,-1,-zvalue, (zvalue*zvalue + 2)],[xvalue1, 1, 1, -(xvalue1*xvalue1 + 2)],[1,yvalue1,1, -(yvalue1*yvalue1 + 2)],[1,1,zvalue1, -(zvalue1*zvalue1 + 2)]]
      
     
     xd_sState.mask_vol_cube = xd_sState.ori_mask_vol_cube
     xd_sState.mask_vol_cube_modify = xd_sState.ori_mask_vol_cube
     xd_sState.new_mask_vol_cube = xd_sState.ori_mask_vol_cube
     index = WHERE(xd_sState.new_mask_vol_cube, count)
     IF count NE 0 THEN xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)

     
     xd_sstate.opoly[0]->getProperty, DATA = Outverts_old
     ;define to the rotation center
     IF N_ELEMENTS(Outverts_old) LT N_ELEMENTS(XD_SSTATE.sel_one.verts) THEN BEGIN
       xd_sstate.opoly[0]->setProperty, DATA =XD_SSTATE.sel_one.verts
       xd_sstate.opoly[0]->setProperty, POLYGONS =  XD_SSTATE.sel_one.conn
     ENDIF ELSE BEGIN
       xd_sstate.opoly[0]->setProperty, POLYGONS =  XD_SSTATE.sel_one.conn
       xd_sstate.opoly[0]->setProperty, DATA = XD_SSTATE.sel_one.verts
     ENDELSE
    
      verts = XD_SSTATE.sel_one.verts
      xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
      xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
      xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
      verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
        [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
        [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
        [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
        [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
        [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
        [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
        [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
      xd_sState.oPoly_border->setProperty, DATA = verts_border

      ;define to the rotation center
      ma=max(XD_SSTATE.sel_one.verts)
      mi=min(XD_SSTATE.sel_one.verts)
      bias=(ma-mi)/2
      rl=mi+bias
      rr=ma+bias
      cc=[(-rl)/(rr-rl),1/(rr-rl)]
      xma = max((XD_SSTATE.sel_one.verts)[0,*])
      xmi = min((XD_SSTATE.sel_one.verts)[0,*])
      xmid = 0.5*(xma+xmi)*cc[1]
      yma = max((XD_SSTATE.sel_one.verts)[1,*])
      ymi = min((XD_SSTATE.sel_one.verts)[1,*])
      ymid = 0.5*(yma+ymi)*cc[1]
      zma = max((XD_SSTATE.sel_one.verts)[2,*])
      zmi = min((XD_SSTATE.sel_one.verts)[2,*])
      zmid = 0.5*(zma+zmi)*cc[1]

      xd_sstate.opoly[0]->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oPoly1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oPoly2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oPoly_border->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oS0->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oS1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oS2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
      xd_sState.ogreenaxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oyellowaxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oblueAxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
      xd_sState.oLine1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oLine2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      xd_sState.oLine3->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
      s = obj_new('orb',radius=0.1,shading=1)
      s->getproperty, POLYGONS=conn
      s->getproperty, DATA=verts
      xd_sState.dataXYZ0=fltarr(3)
      xd_sState.dataXYZ1=fltarr(3)
      xd_sState.dataXYZ2=fltarr(3)
      xd_sState.COUNT123=0
      xd_sState.dataXYZind=0


      xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
      xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
      xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
      xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
      xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
      xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
      xd_sState.ogreenaxis1->setProperty,  hide=1
      xd_sState.oyellowaxis1->setProperty, hide=1
      xd_sState.oblueAxis1->setProperty,  hide=1
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713

      xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
      xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
      xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
      
      xd_sState.ckr1->setProperty,hide=1
      xd_sState.ckg->setProperty,hide=1
      xd_sState.ckr->setProperty, hide=1
      xd_sState.ckb->setProperty, hide=1
      xd_sState.ckag->setProperty,hide=1

      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY

    end

;    'hfsliceonybutton': begin
;
;      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
;      xd_sState.flag_select = sevent.select
;      ;xd_sState.flag_other = sevent.select
;
;      if xd_sState.flag_select then begin
;widget_control,xd_sState. hSlice_on_x_button, SENSITIVE = 0
;        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
;        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
;        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
;        xd_sState.opoly1->SetProperty, $
;          CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
;      endif else begin
;widget_control,xd_sState. hSlice_on_x_button, SENSITIVE = 1
;        xd_sState.opoly1->SetProperty, $
;          CLIP_PLANES=[1, 1, 1, -3]
;
;
;      endelse
;
;      xd_sState.oWindow->Draw, xd_sState.oView
;
;      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
;
;    end
    'hfsliceonybutton': begin

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      xd_sState.flag_select = sevent.select
      ;xd_sState.flag_other = sevent.select

      if xd_sState.flag_select then begin
        ;widget_control,xd_sState. hSlice_on_x_button, SENSITIVE = 0
        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
        xd_sState.opoly1->SetProperty, $
          ALPHA_CHANNEL=0.3
          xd_sstate.opoly[0]->SetProperty, $
          ALPHA_CHANNEL=1
      endif else begin
       ; widget_control,xd_sState. hSlice_on_x_button, SENSITIVE = 1
         xd_sstate.opoly[0]->SetProperty, $
          ALPHA_CHANNEL=0.3
          xd_sState.opoly1->SetProperty, $
          ALPHA_CHANNEL=1

      endelse

      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState

      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY

    end





    'hfslicesave': begin

      WIDGET_CONTROL, sEvent.top, gET_UVALUE=xd_sState, /NO_COPY

      ;20180318hxn

      ;
      ;        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      ;        print,'xvalue',xvalue
      ;        print,'yvalue',yvalue
      ;        print,'zvalue',zvalue
      ;        print,'dvalue',-(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)
      ;xd_sstate.opoly[0]->getProperty,data=v1
      ;xd_sstate.opoly[0]->getProperty,polygons=p1
      ;opoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255, 127, 127], v1, POLYGON=p1,/NO_COPY)
      ;        opoly1->SetProperty, $
      ;          CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
      ;rTopModel1 = OBJ_NEW('IDLgrModel')
      ;rTopModel1->Add,opoly1
      ;        ;(*pState).rWindow->Draw, (*pState).rView
      ;        xobjview_XDU,rTopModel1, BACKGROUND = [0, 0, 0]

      ;hxn20180324hxn
      ;        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      ;        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      ;        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      ;        xd_sstate.opoly[0]->getProperty,data=v1
      ;        xd_sstate.opoly[0]->getProperty,polygons=p1
      ;        print,'sava'
      ;opoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255, 127, 127], v1, POLYGON=p1,/NO_COPY)
      ;
      ;xd_sState.opoly = OBJ_NEW('IDLgrPolygon',COLOR = [255, 127, 127], v1, POLYGON=p1,CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)],/NO_COPY)
      ;opoly1->SetProperty, CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
      ;xd_sState.oRotationModel1->Add,opoly1
      ;xd_sState.oWindow->Draw, xd_sState.oView

      ;20180327hxn
      widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      xd_sState.opoly1->SetProperty, CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
      xd_sState.opoly2->SetProperty, CLIP_PLANES=[xvalue, yvalue, zvalue, -(xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]
      xd_sstate.opoly[0]->SetProperty, CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue*xvalue + yvalue*yvalue + zvalue*zvalue)]


      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState


      WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end

    'hfsliceoutput': begin
       WIDGET_CONTROL, sEvent.top, gET_UVALUE=xd_sState, /NO_COPY

      widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
      widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
      widget_control, xd_sState.hSlice_z_slider, get_value = zvalue
      widget_control, xd_sState.hSlice_x_slider1, get_value = xvalue1
      widget_control, xd_sState.hSlice_y_slider1, get_value = yvalue1
      widget_control, xd_sState.hSlice_z_slider1, get_value = zvalue1

            
            xd_sState.clippoint=[xvalue,xvalue1,yvalue,yvalue1,zvalue,zvalue1]
            print,xd_sState.clippoint
            modify_roi = 1
            Mask_segmentation, xd_sState, modify_roi
            ortho2D_plot, xd_sState

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
            s = obj_new('orb',radius=0.1,shading=1)
            s->getproperty, POLYGONS=conn
            s->getproperty, DATA=verts
            xd_sState.dataXYZ0=fltarr(3)
            xd_sState.dataXYZ1=fltarr(3)
            xd_sState.dataXYZ2=fltarr(3)
            xd_sState.COUNT123=0
            xd_sState.dataXYZind=0


            xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
            xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
            xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
            xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
            xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
            xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
            xd_sState.ogreenaxis1->setProperty,  hide=1
            xd_sState.oyellowaxis1->setProperty, hide=1
            xd_sState.oblueAxis1->setProperty,  hide=1
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713

            xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]

            xd_sState.ckr1->setProperty,hide=1
            xd_sState.ckg->setProperty,hide=1
            xd_sState.ckr->setProperty, hide=1
            xd_sState.ckb->setProperty, hide=1
            xd_sState.ckag->setProperty,hide=1
            
            
            xd_sState.oWindow->Draw, xd_sState.oView
            ortho2D_plot, xd_sState
            ;draw_3DView, xd_sState
            WIDGET_CONTROL, sEvent.top, sET_UVALUE=xd_sState, /NO_COPY
    end


    ;;;;;;;hxn1-end;;;;;;;
   
    
    
    
    ;  婵炴潙顑夐崳鍝勵湤鎼达紕澶勯柡宥囨嚀閺勫倻锟介敓锟�    ;
    'ZGG' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      IF xd_sState.dataXYZdistance GT 4.0 THEN BEGIN
        ZGGString = STRING(xd_sState.dataXYZdistance, $
          FORMAT='("ZGG Thickness = ", F6.2, " mm!")')
        result_zGG = DIALOG_MESSAGE(ZGGString + ' Please perform pedical screw palcement!', /INFORMATION, TITLE='Suggestions');;DIALOG_MESSAGE('閻犲洨鍏橀敓浠嬫偨閵婏富妯嬬�顔芥尰閻楀绱旈锟藉珶', /INFORMATION, TITLE='閻犲洤锕ラ弻鍥ь嚈妤︽鍞�)
        WIDGET_CONTROL, xd_sState.wZGG, SENSITIVE = 0
        WIDGET_CONTROL, xd_sState.wCKcomplete, SENSITIVE = 0
        WIDGET_CONTROL, xd_sState.wCKcompleteNot, SENSITIVE = 0
        WIDGET_CONTROL, xd_sState.wZB, SENSITIVE = 0
      ENDIF ELSE BEGIN
        WIDGET_CONTROL, xd_sState.wZGG, SENSITIVE = 0
        WIDGET_CONTROL, xd_sState.wCKcomplete, SENSITIVE = 1
        WIDGET_CONTROL, xd_sState.wCKcompleteNot, SENSITIVE = 1
        WIDGET_CONTROL, xd_sState.wZB, SENSITIVE = 0
      ENDELSE
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ; of ZGG

    ;  濞撴皜鍕垫▼閻庣懓鏈弳锟�
    ;
    'CKcomplete' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      result_CKcomplete = DIALOG_MESSAGE('CK is complete! Please perform lateral mass screw palcement!', /INFORMATION, TITLE='Suggestions');;DIALOG_MESSAGE('閻犲洨鍏橀敓浠嬫偨閵娿倖娅犻柛褎顨熼悿鍡涙煢閿燂拷 /INFORMATION, TITLE='閻犲洤锕ラ弻鍥ь嚈妤︽鍞�)
      WIDGET_CONTROL, xd_sState.wZGG, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wCKcomplete, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wCKcompleteNot, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wZB, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ; of CKcomplete
    
    ;  濞撴皜鍕垫▼閻庣懓鏈弳锟�
    ;
    'CKcompleteNot' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wZGG, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wCKcomplete, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wCKcompleteNot, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wZB, SENSITIVE = 1
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ; of CKcompleteNot
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    'colorbutton' : begin
    end
    
    'colorbutton1' : begin
    end
    ;;;;;;;;;;;;;;;;;;;;;;;;;;hjx1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    
    
    
    'ZB' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      ZBString = STRING(xd_sState.dataXYZdistance, $
        FORMAT='("ZB Thickness = ", F6.2, " mm!")')
      result_CKcomplete = DIALOG_MESSAGE(ZBString, /INFORMATION, TITLE='Suggestions');;DIALOG_MESSAGE('閻犲洨鍏橀敓浠嬫偨閵婏富妯嬮柡澶庢硶閻ゅ棝鏌﹂敓锟�/INFORMATION, TITLE='閻犲洤锕ラ弻鍥ь嚈妤︽鍞�)
      WIDGET_CONTROL, xd_sState.wZGG, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wCKcomplete, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wCKcompleteNot, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wZB, SENSITIVE = 0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ; of ZB
    
    ;  Change the font.
    ;
    'FONT' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      case sEvent.index of
        0: xd_sState.oFont->SetProperty, NAME='Helvetica'
        1: xd_sState.oFont->SetProperty, NAME='Times'
        2: begin
          if (!VERSION.OS_FAMILY eq "Windows") then $
            xd_sState.oFont->SetProperty, NAME='Courier New' $
          else xd_sState.oFont->SetProperty, NAME='Courier'
        end    ; of case 2
      endcase
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ; of FONT

    ;  Set the shading to flat.
    ;
    'FLAT' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oSurface->SetProperty, SHADING=0
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wFlatButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wGouraudButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of FLAT

    ;  Set the shading to gouraud.
    ;
    'GOURAUD' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oSurface->SetProperty, SHADING=1
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wFlatButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wGouraudButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of GOURAUD

    ;  Set the style to point.
    ;
    'POINT' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, SENSITIVE=0
      xd_sState.oSurface->SetProperty, STYLE=0
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of POINT

    ;  Set the style to wire.
    ;
    'WIRE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, SENSITIVE=0
      xd_sState.oSurface->SetProperty, STYLE=1
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of POINT

    ;  Set the style to solid.
    ;
    'SOLID' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, /SENSITIVE

      ;  Make linestyle solid, because it will show when dragging
      ;  at low quality.
      ;
      j = d_surfviewToggleOffOn(xd_sState.wLineStyleButton)
      case j of
        0: xd_sState.oSurface->SetProperty, LINESTYLE=0 ; solid
        1: if d_surfviewToggleOffOn(xd_sState.wLineStyleButton) ne 0 $
          then $
          PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
        -1: PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
      endcase

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oSurface->SetProperty, STYLE=2
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of SOLID

    ;  Set the style to ruled xz.
    ;
    'RULEDXZ' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, SENSITIVE=0
      xd_sState.oSurface->SetProperty, STYLE=3
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of RULEDXZ

    ;  Set the style to ruled yz.
    ;
    'RULEDYZ' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, SENSITIVE=0
      xd_sState.oSurface->SetProperty, STYLE=4
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of RULEDYZ

    ;  Set the style to lego wire.
    ;
    'LEGOWIRE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, SENSITIVE=0
      xd_sState.oSurface->SetProperty, STYLE=5
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of LEGOWIRE

    ;  Set the style to lego solid.
    ;
    'LEGOSOLID' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      ;  Make linestyle solid, because it will show when dragging
      ;  at low quality.
      ;
      j = d_surfviewToggleOffOn(xd_sState.wLineStyleButton)
      case j of
        0: xd_sState.oSurface->SetProperty, LINESTYLE=0 ; solid
        1: if d_surfviewToggleOffOn(xd_sState.wLineStyleButton) ne 0 $
          then $
          PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
        -1: PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
      endcase
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTracingButton, SENSITIVE=0
      xd_sState.oSurface->SetProperty, STYLE=6
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wPointButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSolidButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledXZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wRuledYZButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoWireButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wLegoSolidButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wHiddenButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wLineStyleButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wShadingButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of LEGOSOLID

    ; set hole size
    ;
    'HOLESIZE': begin   ; this handles event from either slider
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wHoleSlider, GET_VALUE=hole_size
      xd_sState.initHoleValue = hole_size
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of HOLESIZE


   ;
   ;
   'FUSIONMODE': begin   ; this handles event from fusion
     WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
     WIDGET_CONTROL, sEvent.top, /HOURGLASS
     widget_control, xd_sState.wFusionModeRadio, get_value=fusionmode
     print, 'fusion mode:', fusionmode
     xd_sState.fusionmode = fusionmode
     
     WIDGET_CONTROL, sEvent.top, /HOURGLASS
     xd_sState.oWindow->Draw, xd_sState.oView
     ;draw_3DView, xd_sState

     IF fusionmode EQ 1 THEN widget_control, xd_sState.wfusionViewBase, map=1 ELSE widget_control, xd_sState.wfusionViewBase, map=0

     ;replot the 3 draw windows
     ortho2D_plot, xd_sState

     WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'XFLIPMODE': begin   ; this handles event from x flip
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      widget_control, xd_sState.wxFlipModeRadio, get_value=xflipmode
      print, 'x flip mode:', xflipmode
      delta_xflip = xflipmode - xd_sState.xflipmode
      xd_sState.xflipmode = xflipmode
      IF ABS(delta_xflip) GT 0.001 THEN BEGIN
        xd_sState.vol_HU_cube_a1 = REVERSE(xd_sState.vol_HU_cube_a1, 1)
        xd_sState.svol_HU_cube_a1 = CONGRID(xd_sState.vol_HU_cube_a1, xd_sState.imageSize, xd_sState.imageSize, xd_sState.imageSize)
      ENDIF
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'YFLIPMODE': begin   ; this handles event from y flip
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      widget_control, xd_sState.wyFlipModeRadio, get_value=yflipmode
      print, 'y flip mode:', yflipmode
      delta_yflip = yflipmode - xd_sState.yflipmode
      xd_sState.yflipmode = yflipmode
      IF ABS(delta_yflip) GT 0.001 THEN BEGIN
        xd_sState.vol_HU_cube_a1 = REVERSE(xd_sState.vol_HU_cube_a1, 2)
        xd_sState.svol_HU_cube_a1 = CONGRID(xd_sState.vol_HU_cube_a1, xd_sState.imageSize, xd_sState.imageSize, xd_sState.imageSize)
      ENDIF
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'ZFLIPMODE': begin   ; this handles event from x flip
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      widget_control, xd_sState.wzFlipModeRadio, get_value=zflipmode
      print, 'z flip mode:', zflipmode
      delta_zflip = zflipmode - xd_sState.zflipmode
      xd_sState.zflipmode = zflipmode
      IF ABS(delta_zflip) GT 0.001 THEN BEGIN
        xd_sState.vol_HU_cube_a1 = REVERSE(xd_sState.vol_HU_cube_a1, 3)
        xd_sState.svol_HU_cube_a1 = CONGRID(xd_sState.vol_HU_cube_a1, xd_sState.imageSize, xd_sState.imageSize, xd_sState.imageSize)
      ENDIF
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'xmoveslider': begin   ; this handles event from x move
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      widget_control, xd_sState.wx_move_slider, get_value=x_move
      print, 'x move:', x_move
      delta_x_move = x_move - xd_sState.x_move
      xd_sState.vol_HU_cube_a1 = SHIFT(xd_sState.vol_HU_cube_a1, delta_x_move, 0, 0)
      xd_sState.svol_HU_cube_a1 = CONGRID(xd_sState.vol_HU_cube_a1, xd_sState.imageSize, xd_sState.imageSize, xd_sState.imageSize)
      xd_sState.x_move = x_move
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'ymoveslider': begin   ; this handles event from y move
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      widget_control, xd_sState.wy_move_slider, get_value=y_move
      print, 'y move:', y_move
       delta_y_move = y_move - xd_sState.y_move
      xd_sState.vol_HU_cube_a1 = SHIFT(xd_sState.vol_HU_cube_a1, 0, delta_y_move, 0)
      xd_sState.svol_HU_cube_a1 = CONGRID(xd_sState.vol_HU_cube_a1, xd_sState.imageSize, xd_sState.imageSize, xd_sState.imageSize)
      xd_sState.y_move = y_move
     WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'zmoveslider': begin   ; this handles event from z move
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      widget_control, xd_sState.wz_move_slider, get_value=z_move
      print, 'z move:', z_move
      delta_z_move = z_move - xd_sState.z_move
      xd_sState.vol_HU_cube_a1 = SHIFT(xd_sState.vol_HU_cube_a1, 0, 0, delta_z_move)
      xd_sState.svol_HU_cube_a1 = CONGRID(xd_sState.vol_HU_cube_a1, xd_sState.imageSize, xd_sState.imageSize, xd_sState.imageSize)
      xd_sState.z_move = z_move
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'xyzflipmove_save': begin   ; this handles event from saving xyx flip and move results
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'xyzflipmove_apply': begin   ; this handles event from applying xyx flip and move results
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    ;  Set threshold value.
    ;
    'THRESHOLD_min': begin   ; this handles event from either slider
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdSlider_min, GET_VALUE=threshold_min_value
      WIDGET_CONTROL, xd_sState.wThresholdText_min, SET_VALUE=string(format='(I5)', threshold_min_value)
;      xd_sState.initThresholdValue_min = threshold_min_value
      threshold_max_value = xd_sState.initThresholdValue_max
      xd_sState.sel_one.HU_min_value = threshold_min_value
      IF threshold_min_value GT threshold_max_value THEN BEGIN
        temp = threshold_min_value
        threshold_min_value = threshold_max_value
        threshold_max_value = temp
      ENDIF
      
      Add_new_by_thresholding, xd_sState
      
;      ;mask segmentation
;      modify_roi = 24
;      Mask_segmentation, xd_sState, modify_roi
;      ;GJ 2018/12/22 testing HJX's code
;;      xd_sstate.opoly[0]->getProperty, DATA = smoothedOutverts
;;      change_color,smoothedOutverts, xd_sState.new_mask_vol_cube
;     
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_min

    ;  Set threshold value.
    ;
    'THRESHOLD_min_text': begin   ; this handles event from either text
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdText_min, GET_VALUE=temp
      threshold_min_value = LONG(temp)
      IF threshold_min_value GT xd_sState.maxValue THEN threshold_min_value = xd_sState.maxValue
      IF threshold_min_value LT xd_sState.minValue THEN threshold_min_value = xd_sState.minValue
      WIDGET_CONTROL, xd_sState.wThresholdText_min, SET_VALUE=string(format='(I5)', threshold_min_value)
      WIDGET_CONTROL, xd_sState.wThresholdSlider_min, SET_VALUE=threshold_min_value
      ;      xd_sState.initThresholdValue_min = threshold_min_value
      threshold_max_value = xd_sState.initThresholdValue_max
      xd_sState.sel_one.HU_min_value = threshold_min_value
      IF threshold_min_value GT threshold_max_value THEN BEGIN
        temp = threshold_min_value
        threshold_min_value = threshold_max_value
        threshold_max_value = temp
      ENDIF


      Add_new_by_thresholding, xd_sState

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_min


    ;  Set threshold value.
    ;
    'THRESHOLD_max': begin   ; this handles event from either slider
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdSlider_max, GET_VALUE=threshold_max_value
      WIDGET_CONTROL, xd_sState.wThresholdText_max, SET_VALUE=string(format='(I5)', threshold_max_value)
;      xd_sState.initThresholdValue_max = threshold_max_value
      threshold_min_value = xd_sState.initThresholdValue_min
      xd_sState.sel_one.HU_max_value = threshold_max_value
      IF threshold_min_value GT threshold_max_value THEN BEGIN
        temp = threshold_min_value
        threshold_min_value = threshold_max_value
        threshold_max_value = temp
      ENDIF
      
      Add_new_by_thresholding, xd_sState
     
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState
      
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_max
    
    ;  Set threshold value.
    ;
    'THRESHOLD_max_text': begin   ; this handles event from either text
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdText_max, GET_VALUE=temp
      threshold_max_value = LONG(temp)
      IF threshold_max_value GT xd_sState.maxValue THEN threshold_max_value = xd_sState.maxValue
      IF threshold_max_value LT xd_sState.minValue THEN threshold_max_value = xd_sState.minValue
      WIDGET_CONTROL, xd_sState.wThresholdText_max, SET_VALUE=string(format='(I5)', threshold_max_value)
      WIDGET_CONTROL, xd_sState.wThresholdSlider_max, SET_VALUE=threshold_max_value
      threshold_min_value = xd_sState.initThresholdValue_min
      xd_sState.sel_one.HU_max_value = threshold_max_value
      IF threshold_min_value GT threshold_max_value THEN BEGIN
        temp = threshold_min_value
        threshold_min_value = threshold_max_value
        threshold_max_value = temp
      ENDIF
      
      Add_new_by_thresholding, xd_sState
     
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState
      
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_max
    
    ;  Set threshold value.
    ;
    'THRESHOLD_min_a1': begin   ; this handles event from either slider
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdSlider_min_a1, GET_VALUE=threshold_min_value_a1
      WIDGET_CONTROL, xd_sState.wThresholdText_min_a1, SET_VALUE=string(format='(I5)', threshold_min_value_a1)
      ;      xd_sState.initThresholdValue_min = threshold_min_value
      threshold_max_value_a1 = xd_sState.initThresholdValue_max_a1
;      xd_sState.sel_one.HU_min_value = threshold_min_value
      IF threshold_min_value_a1 GT threshold_max_value_a1 THEN BEGIN
        temp = threshold_min_value_a1
        threshold_min_value_a1 = threshold_max_value_a1
        threshold_max_value_a1 = temp
      ENDIF
      
      xd_sState.threshold_min_value_a1 = threshold_min_value_a1
      xd_sState.threshold_max_value_a1 = threshold_max_value_a1
      Add_new_by_thresholding_a1, xd_sState

      ;      ;mask segmentation
      ;      modify_roi = 24
      ;      Mask_segmentation, xd_sState, modify_roi
      ;      ;GJ 2018/12/22 testing HJX's code
      ;;      xd_sstate.opoly[0]->getProperty, DATA = smoothedOutverts
      ;;      change_color,smoothedOutverts, xd_sState.new_mask_vol_cube
      ;
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_min

    ;  Set threshold value.
    ;
    'THRESHOLD_min_text_a1': begin   ; this handles event from either text
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdText_min_a1, GET_VALUE=temp
      threshold_min_value_a1 = LONG(temp)
      IF threshold_min_value_a1 GT xd_sState.maxValue_a1 THEN threshold_min_value_a1 = xd_sState.maxValue_a1
      IF threshold_min_value_a1 LT xd_sState.minValue_a1 THEN threshold_min_value_a1 = xd_sState.minValue_a1
      WIDGET_CONTROL, xd_sState.wThresholdText_min_a1, SET_VALUE=string(format='(I5)', threshold_min_value_a1)
      WIDGET_CONTROL, xd_sState.wThresholdSlider_min_a1, SET_VALUE=threshold_min_value_a1
      ;      xd_sState.initThresholdValue_min = threshold_min_value
      threshold_max_value_a1 = xd_sState.initThresholdValue_max_a1
;      xd_sState.sel_one.HU_min_value = threshold_min_value
      IF threshold_min_value_a1 GT threshold_max_value_a1 THEN BEGIN
        temp = threshold_min_value_a1
        threshold_min_value_a1 = threshold_max_value_a1
        threshold_max_value_a1 = temp
      ENDIF

      xd_sState.threshold_min_value_a1 = threshold_min_value_a1
      xd_sState.threshold_max_value_a1 = threshold_max_value_a1
      
      Add_new_by_thresholding_a1, xd_sState

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_min


    ;  Set threshold value.
    ;
    'THRESHOLD_max_a1': begin   ; this handles event from either slider
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdSlider_max_a1, GET_VALUE=threshold_max_value_a1
      WIDGET_CONTROL, xd_sState.wThresholdText_max_a1, SET_VALUE=string(format='(I5)', threshold_max_value_a1)
      ;      xd_sState.initThresholdValue_max = threshold_max_value
      threshold_min_value_a1 = xd_sState.initThresholdValue_min_a1
;      xd_sState.sel_one.HU_max_value = threshold_max_value
      IF threshold_min_value_a1 GT threshold_max_value_a1 THEN BEGIN
        temp = threshold_min_value_a1
        threshold_min_value_a1 = threshold_max_value_a1
        threshold_max_value_a1 = temp
      ENDIF
      
      xd_sState.threshold_min_value_a1 = threshold_min_value_a1
      xd_sState.threshold_max_value_a1 = threshold_max_value_a1
      
      Add_new_by_thresholding_a1, xd_sState

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_max

    ;  Set threshold value.
    ;
    'THRESHOLD_max_text_a1': begin   ; this handles event from either text
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wThresholdText_max, GET_VALUE=temp
      threshold_max_value_a1 = LONG(temp)
      IF threshold_max_value_a1 GT xd_sState.maxValue_a1 THEN threshold_max_value_a1 = xd_sState.maxValue_a1
      IF threshold_max_value_a1 LT xd_sState.minValue_a1 THEN threshold_max_value = xd_sState.minValue_a1
      WIDGET_CONTROL, xd_sState.wThresholdText_max_a1, SET_VALUE=string(format='(I5)', threshold_max_value_a1)
      WIDGET_CONTROL, xd_sState.wThresholdSlider_max_a1, SET_VALUE=threshold_max_value_a1
      threshold_min_value_a1 = xd_sState.initThresholdValue_min_a1
;      xd_sState.sel_one.HU_max_value = threshold_max_value
      IF threshold_min_value_a1 GT threshold_max_value_a1 THEN BEGIN
        temp = threshold_min_value_a1
        threshold_min_value_a1 = threshold_max_value_a1
        threshold_max_value_a1 = temp_a1
      ENDIF
      
      
      xd_sState.threshold_min_value_a1 = threshold_min_value_a1
      xd_sState.threshold_max_value_a1 = threshold_max_value_a1
      Add_new_by_thresholding_a1, xd_sState

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of THRESHOLD_max
    
    
    'wtransSlider': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wtransslider, GET_VALUE=temp
      wtransslidervalue=float(temp)/100.0
      xd_sState.wtransslidervalue = wtransslidervalue
      WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(wtransslidervalue*100, format='(f5.1)')

      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;whether there is enough verts
      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF
      
      ;modify the object's transparency
      xd_sState.selected->setProperty, ALPHA_CHANNEL=wtransslidervalue
      xd_sState.selected->setProperty, data=verts

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

    end

    'Trans_text': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTrans_text, GET_VALUE=temp
      wtransslidervalue=float(temp)/100.0
      IF wtransslidervalue GT 1. THEN wtransslidervalue = 1.
      IF wtransslidervalue LT 0. THEN wtransslidervalue = 0.
      
      WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=wtransslidervalue*100
      WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(wtransslidervalue*100, format='(f5.1)')
      
      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;modify the object's transparency
      xd_sState.selected->setProperty, ALPHA_CHANNEL=wtransslidervalue
      xd_sState.selected->setProperty, data=verts

      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
    
    'Smooth': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wSmooth, GET_VALUE=smoothvalue
      WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothvalue, format='(f5.1)')

      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      xd_sState.selected->getProperty, data=verts, NAME=name_poly, POLYGONS=conn
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      PRINT, 'smooth value ', smoothvalue
      IF smoothvalue EQ 0 THEN BEGIN
        smoothedOutverts = *(xd_sState.struc_modifyobject.p_ori_verts_list[LONG(name_poly)])
      ENDIF ELSE BEGIN
        smoothedOutverts = MESH_SMOOTH(*(xd_sState.struc_modifyobject.p_ori_verts_list[LONG(name_poly)]), $
          conn, /FIXED_EDGE_VERTICES, ITERATIONS=smoothvalue*100, LAMBDA=0.5);01);5);0.01)
;          ;GJ 2020/5/21, adding lambda, i.e. from 0.05 to 0.5
;        numberVertices = MESH_DECIMATE(smoothedOutverts, conn, $
;            decimatedConn, VERTICES = decimatedVerts, PERCENT_VERTICES = 20)
;        print, 'number of vertices = ', numberVertices
;        smoothedOutverts = decimatedVerts
;        conn = decimatedConn
      ENDELSE
      ;assign selected object with smoothed verts,use thick property to record smoothnews
      ;xd_sState.selected->setProperty, DATA = smoothedOutverts, thick = smoothvalue/1000.
      xd_sState.selected->setProperty, DATA = smoothedOutverts, POLYGONS=conn, thick = smoothvalue/1000.

      smoothvalue = smoothvalue/ 400.0
      smoothvalue = smoothvalue*100.0

;      SmoothString = 'Smooth : ' + STRING(smoothvalue, $
;        FORMAT='(f5.1)') + ' %'
;      print,SmoothString
;      WIDGET_CONTROL, xd_sState.wSmoothLabel, $
;        SET_VALUE=SmoothString


      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

    end
    
    'Smooth_text': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wSmooth_text, GET_VALUE=temp
      smoothvalue = LONG(temp)
      IF smoothvalue GT 400 THEN smoothvalue = 400
      IF smoothvalue LT 0 THEN smoothvalue = 0
      WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothvalue, format='(f5.1)')
      WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=smoothvalue

      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      xd_sState.selected->getProperty, data=verts, NAME=name_poly, POLYGONS=conn
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      PRINT, 'smooth value ', smoothvalue
      IF smoothvalue EQ 0 THEN BEGIN
        smoothedOutverts = *(xd_sState.struc_modifyobject.p_ori_verts_list[LONG(name_poly)])
      ENDIF ELSE BEGIN
        smoothedOutverts = MESH_SMOOTH(*(xd_sState.struc_modifyobject.p_ori_verts_list[LONG(name_poly)]), $
          conn, ITERATIONS=smoothvalue)
      ENDELSE
      ;assign selected object with smoothed verts; use thick property to record smooth value
      xd_sState.selected->setProperty, DATA = smoothedOutverts, thick = smoothvalue/1000.
      
      smoothvalue = smoothvalue/ 400.0
      smoothvalue = smoothvalue*100.0

      ;      SmoothString = 'Smooth : ' + STRING(smoothvalue, $
      ;        FORMAT='(f5.1)') + ' %'
      ;      print,SmoothString
      ;      WIDGET_CONTROL, xd_sState.wSmoothLabel, $
      ;        SET_VALUE=SmoothString


      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ;draw_3DView, xd_sState

      ;replot the 3 draw windows
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
      
    end
   
    ;  Set scaling
    ;
    'SCALING': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale

      scale = ABS(0.75 + FLOAT(scale) / 100.0)
      scalep = scale*100.0
      scalingString = 'Scaling : ' + STRING(scalep, $
        FORMAT='(f5.1)') + ' %'
;      WIDGET_CONTROL, xd_sState.wScalingLabel, $
;        SET_VALUE=scalingString

;      xd_sState.oRuler1model->getproperty,transform=transform_matrix
;      print, 'transform_matrix = ', transform_matrix
;      print, 'scale = ', scale
;      print, 'scalep = ', scalep
;      print, '(xd_sState.oPoly_size[0]*xd_sState.vol_HU_cube_resolution) = ', (xd_sState.oPoly_size[0]*xd_sState.vol_HU_cube_resolution)
      transform = [[scale, 0, 0, 0.0], [0, scale, 0, 0.0], $
        [0, 0, scale, 0.0], [0, 0, 0, 1]]
      xd_sState.oScalingModel->SetProperty, TRANSFORM=transform
      xd_sState.oTraceScalingModel->SetProperty, TRANSFORM=transform
      scale_ruler = scale*2.3*280./(xd_sState.oPoly_size[0]*xd_sState.vol_HU_cube_resolution)
      xd_sState.oRuler1model->setproperty,transform=[[scale_ruler, 0, 0, 0.0],[0, 1, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
      xd_sState.oRuler2model->setproperty,transform=[[1, 0, 0, 0.0],[0, scale_ruler, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
      xd_sState.oRuler1Model->scale,0.008,0.05,1
      xd_sState.oRuler1Model->translate,-0.2,0.4,0.45
      xd_sState.oRuler2Model->scale,0.05,0.008,1
      xd_sState.oRuler2Model->translate,0.4,-0.2,0.45

      

;      
      ;GJ, 2018/6/8
;      xd_sState.oRotationModel->GetProperty, TRANSFORM=xdTransformRotation
;      xd_sState.oScalingModel->GetProperty, TRANSFORM=xdTransformScaling
;      print, 'xdTransformRotation = ', xdTransformRotation
;      print, 'xdTransformScaling = ', xdTransformScaling
      
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of SCALING

    'CONSTRAINT': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      case sEvent.index of
        0: begin ; unconstrained rotation
          xd_sState.oRotationModel->SetProperty, $
            CONSTRAIN=0
          xd_sState.oTraceRotationModel->SetProperty, $
            CONSTRAIN=0
        end
        1: begin
          xd_sState.oRotationModel->SetProperty, $
            AXIS=0, /CONSTRAIN
          xd_sState.oTraceRotationModel->SetProperty, $
            AXIS=0, /CONSTRAIN
        end
        2: begin
          xd_sState.oRotationModel->SetProperty, $
            AXIS=1, /CONSTRAIN
          xd_sState.oTraceRotationModel->SetProperty, $
            AXIS=1, /CONSTRAIN
        end
        3: begin
          xd_sState.oRotationModel->SetProperty, $
            AXIS=2, /CONSTRAIN
          xd_sState.oTraceRotationModel->SetProperty, $
            AXIS=2, /CONSTRAIN
        end
        4: begin
          xd_sState.oRotationModel->SetProperty, $
            AXIS=0, CONSTRAIN=2
          xd_sState.oTraceRotationModel->SetProperty, $
            AXIS=0, CONSTRAIN=2
        end
        5: begin
          xd_sState.oRotationModel->SetProperty, $
            AXIS=1, CONSTRAIN=2
          xd_sState.oTraceRotationModel->SetProperty, $
            AXIS=1, CONSTRAIN=2
        end
        6: begin
          xd_sState.oRotationModel->SetProperty, $
            AXIS=2, CONSTRAIN=2
          xd_sState.oTraceRotationModel->SetProperty, $
            AXIS=2, CONSTRAIN=2
        end
      endcase

      xd_sState.oRotationModel->GetProperty, CONSTRAIN=constrain
      xd_sState.oRedAxis->GetProperty, HIDE=red_axis_hidden

      if (constrain eq 2) and red_axis_hidden then begin
        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
        xd_sState.oRedAxis->SetProperty, HIDE=0
        xd_sState.oGreenAxis->SetProperty, HIDE=0
        xd_sState.oBlueAxis->SetProperty, HIDE=0
        xd_sState.oWindow->Draw, xd_sState.oView
        ortho2D_plot, xd_sState
        ;draw_3DView, xd_sState
        WIDGET_CONTROL, xd_sState.wHideAxes, /BITMAP, SET_VALUE=xd_sState.icon_dir+'hideaxes.bmp'
        WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      end
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end

    'HOTKEY' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      case STRUPCASE(sEvent.ch) of
        'U': begin ; unconstrained rotation
          xd_sState.oRotationModel->SetProperty, CONSTRAIN=0
          xd_sState.oTraceRotationModel->SetProperty, CONSTRAIN=0
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=0
          redvalue = replicate(245,256,1)
          greenvalue = BINDGEN(256)
          bluevalue = reverse(greenvalue)
          ;;;;;;;;;;;;;;;hjx1226

          FOR i=0, 180 DO BEGIN
            ;GJ 2019/1/28
            axis = [0,1,0]
            angle = 2
            xd_sState.oScalingModel->Rotate, axis, angle
            wloadctslidervalue = FLOOR(i*10.*256./180.) MOD 256
            xd_sstate.opoly[0]->setProperty, color = [redvalue[wloadctslidervalue],greenvalue[wloadctslidervalue],bluevalue[wloadctslidervalue]]
            xd_sState.oWindow->Draw, xd_sState.oView
            ortho2D_plot, xd_sState
            draw_3DView, xd_sState
            wait, 0.02
            ;end GJ 2019/1/28
          ENDFOR
        end
        'X': begin
          xd_sState.oRotationModel->SetProperty, AXIS=0, $
            /CONSTRAIN
          xd_sState.oTraceRotationModel->SetProperty, AXIS=0, $
            /CONSTRAIN
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=1
        end
        'Y': begin
          xd_sState.oRotationModel->SetProperty, AXIS=1, $
            /CONSTRAIN
          xd_sState.oTraceRotationModel->SetProperty, AXIS=1, $
            /CONSTRAIN
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=2
        end
        'Z': begin
          xd_sState.oRotationModel->SetProperty, AXIS=2, $
            /CONSTRAIN
          xd_sState.oTraceRotationModel->SetProperty, AXIS=2, $
            /CONSTRAIN
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=3
        end
        'R': begin
          xd_sState.oRotationModel->SetProperty, AXIS=0, $
            CONSTRAIN=2
          xd_sState.oTraceRotationModel->SetProperty, AXIS=0, $
            CONSTRAIN=2
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=4
        end
        'G': begin
          xd_sState.oRotationModel->SetProperty, AXIS=1, $
            CONSTRAIN=2
          xd_sState.oTraceRotationModel->SetProperty, AXIS=1, $
            CONSTRAIN=2
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=5
        end
        'B': begin
          xd_sState.oRotationModel->SetProperty, AXIS=2, $
            CONSTRAIN=2
          xd_sState.oTraceRotationModel->SetProperty, AXIS=2, $
            CONSTRAIN=2
;          WIDGET_CONTROL, xd_sState.wConstraintsDroplist, $
;            SET_DROPLIST_SELECT=6
        end
        'H': begin ; Toggle hide-show rotation axes
          id = xd_sState.wHideAxes
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
          d_surfview1Event, {id: id, top: sEvent.top, handler:0l}
          WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
        end
        else:
      endcase

      xd_sState.oRotationModel->GetProperty, CONSTRAIN=constrain
      xd_sState.oRedAxis->GetProperty, HIDE=red_axis_hidden

      if (constrain eq 2) $
        and red_axis_hidden $
        and STRUPCASE(sEvent.ch) ne 'H' $
        then begin
        WIDGET_CONTROL, sEvent.top, /HOURGLASS
        WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
        xd_sState.oRedAxis->SetProperty, HIDE=0
        xd_sState.oGreenAxis->SetProperty, HIDE=0
        xd_sState.oBlueAxis->SetProperty, HIDE=0
        xd_sState.oWindow->Draw, xd_sState.oView
        ortho2D_plot, xd_sState
        ;draw_3DView, xd_sState
        WIDGET_CONTROL, xd_sState.wHideAxes, /BITMAP, SET_VALUE=xd_sState.icon_dir+'hideaxes.bmp'
;          SET_VALUE='Hide Axes'
        WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      end

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end

    ;  Hide the Rotation Axes
    ;
    'HIDEAXES': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oRedAxis->GetProperty, HIDE=red_axis_hidden
      xd_sState.oRedAxis->SetProperty, HIDE=1-red_axis_hidden
      xd_sState.oGreenAxis->SetProperty, HIDE=1-red_axis_hidden
      xd_sState.oBlueAxis->SetProperty, HIDE=1-red_axis_hidden
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wHideAxes, /BITMAP, SET_VALUE=xd_sState.icon_dir+(['Show','Hide'])[red_axis_hidden]+'Axes.bmp'
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end

    ;  Reset the initial orientation of the surface
    ;
    'RESETTRANSFORM': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oRotationModel->SetProperty, $
        TRANSFORM=xd_sState.initTransformRotation
      xd_sState.oTraceRotationModel->SetProperty, $
        TRANSFORM=xd_sState.initTransformRotation
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of RESETTRANSFORM

    ;  Toggle off and on the texture mapping.
    ;
    'TEXTURE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      j = d_surfviewToggleOffOn(sEvent.id)
      case j of

        0: begin
          xd_sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
        end    ;   of  0

        1: begin
          xd_sState.oSurface->SetProperty, $
            TEXTURE_MAP=xd_sState.oTextureImage

          ;  Turn off vertex coloring.
          ;
          if d_surfviewToggleOffOn(xd_sState.wVertexButton) eq 1 $
            then $
            void = d_surfviewToggleOffOn(xd_sState.wVertexButton)
          xd_sState.oSurface->SetProperty, VERT_COLORS=0

          if d_surfviewToggleOffOn(xd_sState.wTracingMaskButton) eq 1 $
            then $
            void = d_surfviewToggleOffOn( $
            xd_sState.wTracingMaskButton $
            )
          ;  Turn off tracing mode.
          ;
          if d_surfviewToggleOffOn(xd_sState.wTracingModeButton) eq 1 $
            then $
            void = d_surfviewToggleOffOn( $
            xd_sState.wTracingModeButton $
            )
          xd_sState.tracingMode = 0
          xd_sState.oTracePolyline->SetProperty, /HIDE

          if not xd_sState.tracingMode then $
            WIDGET_CONTROL, xd_sState.wStyleButton, SENSITIVE=1

          xd_sState.oWindow->Draw, xd_sState.oView
          ;ortho2D_plot, xd_sState
          ;draw_3DView, xd_sState
          demo_putTips, xd_sState, ['infor','instr'], [11,12], /LABEL
        end    ;   of  1

        -1: PRINT, 'Error in d_surfviewToggleOffOn/texture.'

      endcase
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of TEXTURE

    ;  Toggle off and on the vertex colors.
    ;
    'VERTEXCOLOR' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      j = d_surfviewToggleOffOn(sEvent.id)
      case j of
        0: xd_sState.oSurface->SetProperty, VERT_COLORS=0
        1: begin

          ;  Turn off texture mapping.
          ;
          if d_surfviewToggleOffOn(xd_sState.wTextureButton) eq 1 $
            then $
            void = d_surfviewToggleOffOn(xd_sState.wTextureButton) $
          else begin
            xd_sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
            xd_sState.oSurface->GetProperty, STYLE=style
            WIDGET_CONTROL, xd_sState.wTracingButton, $
              SENSITIVE=([0,1])[style eq 2];Solid can be traced
          end

          ;  Turn on vertex colors.
          ;
          xd_sState.oSurface->SetProperty, $
            VERT_COLORS=xd_sState.vertexColors
        end
        -1: PRINT, 'Error in d_surfviewToggleOffOn/vertexcolors.'
      endcase
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of VERTEXCOLOR

    ;;;;;;;;;;;;;;;;;;;;;hjx1;;;;;;;;;;;;;;;;;;;;;;;;;;

    'white' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, /HOURGLASS

      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;center the selected object
      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      xd_sState.selected->getProperty, vert_colors = old_vert_colors
      old_vert_colors[0,*] = 255
      old_vert_colors[1,*] = 255
      old_vert_colors[2,*] = 255
      xd_sState.selected->setProperty, vert_colors = old_vert_colors
      WIDGET_CONTROL, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of COLORs


    'blue' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, /HOURGLASS

      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;center the selected object
      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      xd_sState.selected->getProperty, vert_colors = old_vert_colors
      old_vert_colors[0,*] = 176
      old_vert_colors[1,*] = 224
      old_vert_colors[2,*] = 230
      xd_sState.selected->setProperty, vert_colors = old_vert_colors
      WIDGET_CONTROL, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of COLORs

    'pink' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, /HOURGLASS

      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;center the selected object
      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      xd_sState.selected->getProperty, vert_colors = old_vert_colors
      old_vert_colors[0,*] = 255
      old_vert_colors[1,*] = 192
      old_vert_colors[2,*] = 203
      xd_sState.selected->setProperty, vert_colors = old_vert_colors
      WIDGET_CONTROL, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end


    'salmon' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, /HOURGLASS
      
      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;center the selected object
      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF
      
        xd_sState.selected->getProperty, vert_colors = old_vert_colors
        old_vert_colors[0,*] = 255
        old_vert_colors[1,*] = 140
        old_vert_colors[2,*] = 105
        xd_sState.selected->setProperty, vert_colors = old_vert_colors
      WIDGET_CONTROL, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end

    'wloadctSlider' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.wloadctslider, GET_VALUE=wloadctslidervalue
      
      ;whether selected is true
      IF ~OBJ_VALID(xd_sState.selected) THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      IF ~obj_isa(xd_sState.selected,'IDLgrPolygon') THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

      ;center the selected object
      xd_sState.selected->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF

;      redvalue =  BINDGEN(256)
;      greenvalue =  replicate(180,256,1)
;      bluevalue = reverse(redvalue)

        ;;;;;;;;;;;;;;;hjx1226
        redvalue = replicate(245,256,1)
        greenvalue = BINDGEN(256)
        bluevalue = reverse(greenvalue)
        ;;;;;;;;;;;;;;;hjx1226

      xd_sState.wloadctslidervalue = wloadctslidervalue
      xd_sState.selected->getProperty, vert_colors = old_vert_colors
      old_vert_colors[0,*] = redvalue[wloadctslidervalue]
      old_vert_colors[1,*] = greenvalue[wloadctslidervalue]
      old_vert_colors[2,*] = bluevalue[wloadctslidervalue]
      xd_sState.selected->setProperty, vert_colors = old_vert_colors
      WIDGET_CONTROL, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end
   
    ;  Toggle off and on the hidden points and lines.
    ;
    'HIDE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      j = d_surfviewToggleOffOn(sEvent.id)
      if (j eq -1) then begin
        PRINT, 'Error in d_surfviewToggleOffOn/hide.'
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      endif
      xd_sState.oSurface->SetProperty, HIDDEN_LINES=j
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of HIDE

    ;  Toggle between solid and dash linestyles.
    ;
    'LINESTYLE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      j = d_surfviewToggleOffOn(sEvent.id)
      case j of
        0: xd_sState.oSurface->SetProperty, LINESTYLE=0 ; solid
        1: xd_sState.oSurface->SetProperty, LINESTYLE=4
        -1: PRINT, 'Error in d_surfviewToggleOffOn/linestyle.'
      endcase
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of LINESTYLE

    ;  Show no skirt.
    ;
    'SKIRTNONE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oSurface->SetProperty, SHOW_SKIRT = 0
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wSkirtNoneButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wSkirt10Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt20Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt30Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of SKIRTNONE

    ;  Set skirt to -0.5.
    ;
    'SKIRT10' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oSurface->GetProperty, ZRANGE=zrange
      xd_sState.oSurface->SetProperty, SKIRT=zrange[0], /SHOW_SKIRT
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wSkirtNoneButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt10Button, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wSkirt20Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt30Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of SKIR10

    ;  Set skirt to 0.0.
    ;
    'SKIRT20' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oSurface->GetProperty, ZRANGE=zrange
      xd_sState.oSurface->SetProperty, $
        SKIRT=(zrange[1] - zrange[0]) / 2. + zrange[0], $
        /SHOW_SKIRT
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wSkirtNoneButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt10Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt20Button, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wSkirt30Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of SKIR20

    ;  Set skirt to 0.5.
    ;
    'SKIRT30' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=0
      xd_sState.oSurface->GetProperty, ZRANGE=zrange
      xd_sState.oSurface->SetProperty, SKIRT=zrange[1], /SHOW_SKIRT
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, xd_sState.wSkirtNoneButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt10Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt20Button, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wSkirt30Button, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wTopBase, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of SKIR30

    ;  Set drag quality to low.
    ;
    'LOW' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      xd_sState.dragq = 0
      WIDGET_CONTROL, xd_sState.wLowButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wMediumButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHighButton, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of LOW

    ;  Set drag quality to medium.
    ;
    'MEDIUM' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      xd_sState.dragq = 1
      WIDGET_CONTROL, xd_sState.wLowButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wMediumButton, SENSITIVE=0
      WIDGET_CONTROL, xd_sState.wHighButton, SENSITIVE=1
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of MEDIUM

    ;  Set drag quality to high.
    ;
    'HIGH' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      xd_sState.dragq = 2
      WIDGET_CONTROL, xd_sState.wLowButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wMediumButton, SENSITIVE=1
      WIDGET_CONTROL, xd_sState.wHighButton, SENSITIVE=0
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end        ;  of HIGH

    ;  Draw only the surface within the tracing
    ;  contour (trace mask on),
    ;  or draw the whole surface (trace mask off).
    ;  Toggle between these 2 options.
    ;
    'TRACING_MASK' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS

      case d_surfviewToggleOffOn(sEvent.id) of
        0: begin
          xd_sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
          if xd_sState.tracingMode eq 0 then begin
            WIDGET_CONTROL, xd_sState.wStyleButton, SENSITIVE=1
            WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=1
            demo_putTips, $
              xd_sState, $
              ['infor','instr'], $
              [11,12], $
              /LABEL
          end
        end
        1: begin
          xd_sState.oTracePolyline->GetProperty, $
            DATA=data, POLYLINE=connectivityList
          xd_sState.oTracingMask->GetProperty, DATA=idata
          siz = SIZE(idata)
          mask = BYTARR(siz[1], siz[2])
          idata = BYTARR(siz[1], siz[2], 4) + 255b
          if (connectivityList[0] ge 3) then begin
            ;  Get the number of current points minus 1
            ;
            nCurrentPointsM1 = connectivityList[0] - 1
            x = data[0, 0:nCurrentPointsM1] * 10.0
            y = data[1, 0:nCurrentPointsM1] * 10.0
            fill = POLYFILLV(x, y, siz[1], siz[2])
            mask[*, *] = 0
            if fill[0] gt -1 then $
              mask[fill] = 255
          endif else begin
            mask[*, *] = 255
          endelse
          idata[*, *, 3] = mask
          xd_sState.oTracingMask->SetProperty, DATA=idata
          xd_sState.oSurface->SetProperty, $
            TEXTURE_MAP=xd_sState.oTracingMask

          if d_surfviewToggleOffOn(xd_sState.wTextureButton) eq 1 $
            then $
            void = d_surfviewToggleOffOn(xd_sState.wTextureButton)

          WIDGET_CONTROL, sEvent.top, /HOURGLASS

          xd_sState.oWindow->Draw, xd_sState.oView
          ortho2D_plot, xd_sState
          ;draw_3DView, xd_sState
          WIDGET_CONTROL, xd_sState.wStyleButton, SENSITIVE=0
          if xd_sState.tracingMode eq 1 then begin
            demo_putTips, $
              xd_sState, $
              ['reset','erase'], $
              [11,12], $
              /LABEL
          end else begin
            demo_putTips, $
              xd_sState, $
              ['mask2','disp2'], $
              [11,12], $
              /LABEL
          end
        end          ; of 1
        -1: PRINT, 'Error in d_surfviewToggleOffOn/tracingMask.'
      endcase

      xd_sState.oWindow->Draw, xd_sState.oView
      ortho2D_plot, xd_sState
      ;draw_3DView, xd_sState
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end         ; of TRACING_MASK

    ;  Enable (on) or disable (off) the tracing mode.
    ;  Toggle between these 2 options.
    ;
    'TRACING_MODE' : begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      j = d_surfviewToggleOffOn(sEvent.id)
      case j of
        -1: begin
          PRINT, 'Error in d_surfviewToggleOffOn/tracingMode.'
          WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
          RETURN
        end
        1: begin
          xd_sState.tracingMode = 1
          xd_sState.oTracePolyline->GetProperty, $
            POLYLINE=connectivityList
          connectivityList[0] = 0
          connectivityList[1] = -1
          xd_sState.oTracePolyline->SetProperty, $
            POLYLINE=connectivityList, $
            HIDE=0
          if d_surfviewToggleOffOn(xd_sState.wTextureButton) eq 1 $
            then $
            void = d_surfviewToggleOffOn(xd_sState.wTextureButton)
          xd_sState.oSurface->SetProperty, TEXTURE_MAP=OBJ_NEW()
          xd_sState.oWindow->Draw, xd_sState.oView
          ortho2D_plot, xd_sState
          ;draw_3DView, xd_sState

          WIDGET_CONTROL, xd_sState.wTracingMaskButton, SENSITIVE=0

          WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=0
;          WIDGET_CONTROL, xd_sState.wAnimateButton, SENSITIVE=0
          WIDGET_CONTROL, xd_sState.wStyleButton, SENSITIVE=0
          demo_putTips, xd_sState, ['regio','right'], [11,12], /LABEL
        end
        0: begin
          xd_sState.tracingMode = 0
          xd_sState.oTracePolyline->SetProperty, /HIDE
          xd_sState.oWindow->Draw, xd_sState.oView
          ortho2D_plot, xd_sState
          ;draw_3DView, xd_sState
          xd_sState.oSurface->GetProperty, TEXTURE_MAP=texture_map
          if OBJ_VALID(texture_map) then begin
            demo_putTips, $
              xd_sState, $
              ['mask2','disp2'], $
              [11,12], $
              /LABEL
          endif else begin
            WIDGET_CONTROL, xd_sState.wStyleButton, SENSITIVE=1
            demo_putTips, $
              xd_sState, $
              ['infor','instr'], $
              [11,12], $
              /LABEL
          end
          WIDGET_CONTROL, xd_sState.wTextureButton, SENSITIVE=1
;          WIDGET_CONTROL, xd_sState.wAnimateButton, SENSITIVE=1

        end
      endcase
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end         ; of TRACING_MODE


    ;
    'DRAW0' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN
;      print, 'sEvent.press = ', sEvent.press
;      print, 'sEvent.type = ', sEvent.type
      
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
     
      
      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor
      
      if (sEvent.type EQ 1) then begin
        

        ;
        if (sEvent.release EQ 4) then begin
            ;right button press.
            location3D = xd_sState.location3D
            imageSize = xd_sState.imageSize * 3
            vol_HU_cube = xd_sState.vol_HU_cube
            
            ;calculate the additional image's threshold
            sizeVHc = SIZE(xd_sState.svol_HU_cube)
            mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
            image_statistics, xd_sState.svol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
            IF tf[0] GT 0.00001 THEN BEGIN
              threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
            ENDIF ELSE BEGIN
              threshold_max_value = img_max
            ENDELSE
            threshold_min_value = img_min
            threshold = [threshold_min_value, threshold_max_value]
           
            xyImage = DBLARR(imageSize, imageSize)
            temp_xyImage = CONGRID(vol_HU_cube[*, *, location3D[2]*(SIZE(vol_HU_cube))[3]/xd_sState.imageSize], imageSize, imageSize, 1)
            xyImage[*, *] = temp_xyImage[*, *, 0]
            dims = SIZE(xyImage, /DIMENSIONS)
            vol_dims = SIZE(vol_HU_cube, /DIMENSIONS)
            
            ;define the mask
            WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
            IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
              mask_vol_cube = xd_sState.new_mask_vol_cube
              ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
              ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
            ENDIF ELSE BEGIN
              mask_vol_cube = xd_sState.new_mask_vol_cube*0
              mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
            ENDELSE
            
            ;draw new ROI
;            print, 'xd_sState.location3D = ', xd_sState.location3D
;            print, 'xd_sState.imageSize = ', xd_sState.imageSize
            ;GJ 2020/4/11, displaying the image for better contrast
            XD_XROI, (REVERSE(xyImage, 2)), xd_sState.svol_HU_cube, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, ObjectIndex_value=ObjectIndex_value, ROI_GEOMETRY = ROIgeom, direction = 'xy',  title='xy plane ROI', REGIONS_OUT = xyROIout, BLOCK=1
            IF OBJ_VALID(xyROIout[N_ELEMENTS(xyROIout)-1]) THEN BEGIN
              result_question = DIALOG_MESSAGE('Start new 3D reconstruction?', /Question, TITLE='Continue...')
              IF result_question EQ 'Yes' THEN BEGIN
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252 
                s = obj_new('orb',radius=0.1,shading=1)
                s->getproperty, POLYGONS=conn
                s->getproperty, DATA=verts
                xd_sState.dataXYZ0=fltarr(3)
                xd_sState.dataXYZ1=fltarr(3)
                xd_sState.dataXYZ2=fltarr(3)
                xd_sState.COUNT123=0
                xd_sState.dataXYZind=0
                
                
                xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
                xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
                xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
                xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
                xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
                xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
                xd_sState.ogreenaxis1->setProperty,  hide=1
                xd_sState.oyellowaxis1->setProperty, hide=1
                xd_sState.oblueAxis1->setProperty,  hide=1
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
                
                xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                
                xd_sState.ckr1->setProperty,hide=1
                xd_sState.ckg->setProperty,hide=1
                xd_sState.ckr->setProperty, hide=1
                xd_sState.ckb->setProperty, hide=1
                xd_sState.ckag->setProperty,hide=1
    
                ;define a new mask
                maskResult = xyROIout[N_ELEMENTS(xyROIout)-1] -> ComputeMask(DIMENSIONS = dims)
                xy_mask_ROI = ((REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0)
                xd_sState.xy_mask_ROI = xy_mask_ROI
                xyROIout[N_ELEMENTS(xyROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
                print, 'ContOrInde = ', ContOrInde
                print, 'ac = ', ac
                modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;12=addnew; 4=addlayer, 8=removelayer
                print, 'modify roi index = ', modify_roi
                imageSize = xd_sState.imageSize * 3
                Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, xyROIout, modify_roi, ContOrInde, imageSize, direction = 'xy'
                ;not adding red or green objects
                ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi
    
                ;replace old xyROIout in xd_sState
                IF OBJ_VALID(xyROIout[N_ELEMENTS(xyROIout)-1]) THEN xd_sState.xyROIout = xyROIout[N_ELEMENTS(xyROIout)-1]
                ;draw_3DView, xd_sState
              ENDIF
            ENDIF
            
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView


            ;replot the 3 draw windows
            ortho2D_plot, xd_sState
           
        endif else begin
          
            
          imageSize = xd_sState.imageSize
          
          location3D = xd_sState.location3D
          
          location3D[0] = sEvent.x
          
          location3D[1] = imageSize - sEvent.y
          
          xd_sState.location3D = location3D
          
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;     
          vol_HU_cube = xd_sState.vol_HU_cube
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          if (sEvent.release EQ 1) then begin

            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D
            
            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0          
            
  ;;;;;;;;;;;;红绿蓝三点依次显示         
            scale = 0.75 + FLOAT(scale) / 100.0
            s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
            s->getproperty, POLYGONS=conn
            s->getproperty, DATA=verts
            
            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            IF alpha0 EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha1 EQ 0 THEN BEGIN
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              ENDIF ELSE BEGIN
                IF alpha2 EQ 0 THEN BEGIN
                  xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                ENDIF
              ENDELSE
            ENDELSE
            
            IF xd_sState.dataXYZind EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, DATA = verts
              xd_sState.oS0->setProperty, POLYGONS = conn
               xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ0 = dataxyz
              xd_sState.dataXYZind = 1
            ENDIF else begin
              if xd_sState.dataXYZind EQ 1 THEN BEGIN
                xd_sState.oS1->setProperty, DATA = verts
                xd_sState.oS1->setProperty, POLYGONS = conn
                 xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ1 = dataxyz
                xd_sState.dataXYZind = 2
              ENDIF    else begin
                xd_sState.oS2->setProperty, POLYGONS = conn
                xd_sState.oS2->setProperty, DATA = verts
                 xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ2 = dataxyz
                xd_sState.dataXYZind = 0
              ENDelse
            endelse

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
;;;;;;;;;;;;;;;只有两个点时
            datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
            datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
            datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
            if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
              xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

              statusString = STRING(xd_sState.dataXYZdistance, $
                FORMAT='("Thickness = ", F6.2, " mm")')


              demo_putTips, xd_sState, 'locat', 11, /LABEL
              demo_putTips, xd_sState, statusString, 12

              xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
              ;  endif

            endif else begin
;;;;;;;;;;;;;;;出现三个点时
              if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

                if(xd_sState.COUNT123  eq 0) THEN begin

                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.COUNT123=1

                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                  result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'
                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+','+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckg->setProperty,hide=1
                  xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 1) THEN begin

                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                    xd_sState.COUNT123=2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'


                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')
                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckb->setProperty,hide=1
                    xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                  endif else begin
                    if(xd_sState.COUNT123  eq 2) THEN begin
                      xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                      xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                      xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                      xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                      xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                      xd_sState.COUNT123=0
                      print,'xd_sState.COUNT123=2'
                      print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                      print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                      xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                      xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                      xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                      result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                      angle=acos(result)*!radeg
                      print,'2 lines angle:',angle,'degree'

                      print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                      statusString = STRING(xd_sState.dataXYZdistance, $
                        FORMAT='("Distance = ", F6.2, " mm")')

                      statusString1 = STRING(angle, $
                        FORMAT='("Angle = ", F6.2, " degree")')

                      statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                        FORMAT='("Distance = ", F6.2, " mm")')


                      xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                      demo_putTips, xd_sState, 'locat', 11, /LABEL
                      demo_putTips, xd_sState, statusString+statusString1, 12
                      xd_sState.ckr1->setProperty,hide=1
                      xd_sState.ckr->setProperty,hide=1
                      xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                      xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                      xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                    endif
                  endelse
                endelse
              endif
            endelse
 
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;draw_3DView, xd_sState
          endif
          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;    
          ortho2D_plot, xd_sState
         
          WIDGET_CONTROL, xd_sState.DrawSlider[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider[2], SET_VALUE = location3D[1]
          
        endelse
      endif

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

    end        ;  of DRAW0

    'DRAW1' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      
      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor


       if (sEvent.type EQ 1) then begin
        if (sEvent.release EQ 4) then begin
          ;right button press.
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube = xd_sState.vol_HU_cube
          
          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          
          yzImage = DBLARR(imageSize, imageSize)
          temp_yzImage = CONGRID(vol_HU_cube[location3D[0]*(SIZE(vol_HU_cube))[1]/xd_sState.imageSize, *, *], 1, imageSize, imageSize)
          yzImage[*, *] = temp_yzImage[0, *, *]
          dims = SIZE(yzImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube, /DIMENSIONS)
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, yzImage, xd_sState.svol_HU_cube, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, ObjectIndex_value=ObjectIndex_value, direction = 'yz',  title='yz plane ROI', REGIONS_OUT = yzROIout, BLOCK=1
          IF OBJ_VALID(yzROIout[N_ELEMENTS(yzROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start new 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252 
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0
              
              
              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
                        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              
              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1

              ;define a new mask
              maskResult = yzROIout[N_ELEMENTS(yzROIout)-1] -> ComputeMask(DIMENSIONS = dims)
              yz_mask_ROI = (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
              xd_sState.yz_mask_ROI = yz_mask_ROI
              yzROIout[N_ELEMENTS(yzROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
              modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
              print, 'modify roi index = ', modify_roi
              imageSize = xd_sState.imageSize * 3
              Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, yzROIout, modify_roi, ContOrInde, imageSize, direction = 'yz'
              ;not adding red or green objects
              ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi
  
              ;replace old yzROIout in xd_sState
              IF OBJ_VALID(yzROIout[N_ELEMENTS(yzROIout)-1]) THEN xd_sState.yzROIout = yzROIout[N_ELEMENTS(yzROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF
          
          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
          ;replot the 3 draw windows
          ortho2D_plot, xd_sState
  
         endif else begin
        
            ;      print, sEvent.x, sEvent.y
            imageSize = xd_sState.imageSize
            location3D = xd_sState.location3D
            location3D[1] = sEvent.x
            location3D[2] = sEvent.y
            xd_sState.location3D = location3D
  
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
             WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
            vol_HU_cube = xd_sState.vol_HU_cube

            if (sEvent.release EQ 1) then begin

            
              dataxyz = DBLARR(3)
              dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube))[1]))/xd_sState.imageSize
              dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube))[2]))/xd_sState.imageSize
              dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube))[3]))/xd_sState.imageSize
              xd_sState.location3D = location3D

;              scale = 0.75 + FLOAT(scale) / 100.0

;              s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
;              s->getproperty, POLYGONS=conn
;              s->getproperty, DATA=verts
;              ;          print,'dataxyz=',dataxyz,',conn=',conn,',verts=',verts
;              xd_sState.oS3->setProperty, DATA = verts
;              xd_sState.oS3->setProperty, POLYGONS = conn
;              xd_sState.oS3->setProperty, ALPHA_CHANNEL = 1
              

              s1 = OBJ_NEW('IDLgrPolyline', $
                [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
              s2 = OBJ_NEW('IDLgrPolyline', $
                [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
              s3 = OBJ_NEW('IDLgrPolyline', $
                [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


              s1->getproperty, POLYLINES=conn1
              s1->getproperty, DATA=verts1

              s2->getproperty, POLYLINES=conn2
              s2->getproperty, DATA=verts2

              s3->getproperty, POLYLINES=conn3
              s3->getproperty, DATA=verts3

              xd_sState.ogreenaxis1->setProperty, DATA = verts1
              xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
              xd_sState.ogreenaxis1->SetProperty, HIDE=0

              xd_sState.oyellowaxis1->setProperty, DATA = verts2
              xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
              xd_sState.oyellowaxis1->SetProperty, HIDE=0

              xd_sState.oBlueAxis1->setProperty, DATA = verts3
              xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
              xd_sState.oBlueAxis1->SetProperty, HIDE=0
              
              
              WIDGET_CONTROL, sEvent.top, /HOURGLASS
              xd_sState.oWindow->Draw, xd_sState.oView
              ortho2D_plot, xd_sState
              ;;draw_3DView, xd_sState
            endif
   ;;;;;;;;;;;;;;;;;;;;三个点依次显示         
            scale = 0.75 + FLOAT(scale) / 100.0
            s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
            s->getproperty, POLYGONS=conn
            s->getproperty, DATA=verts

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            IF alpha0 EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha1 EQ 0 THEN BEGIN
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              ENDIF ELSE BEGIN
                IF alpha2 EQ 0 THEN BEGIN
                  xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                ENDIF
              ENDELSE
            ENDELSE

            IF xd_sState.dataXYZind EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, DATA = verts
              xd_sState.oS0->setProperty, POLYGONS = conn
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ0 = dataxyz
              xd_sState.dataXYZind = 1
         
            ENDIF else begin
              if xd_sState.dataXYZind EQ 1 THEN BEGIN
                xd_sState.oS1->setProperty, DATA = verts
                xd_sState.oS1->setProperty, POLYGONS = conn
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ1 = dataxyz
                xd_sState.dataXYZind = 2

              ENDIF    else begin
                xd_sState.oS2->setProperty, POLYGONS = conn
                xd_sState.oS2->setProperty, DATA = verts
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ2 = dataxyz
                xd_sState.dataXYZind = 0

              ENDelse
            endelse

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
            ;;;;;;;;;;;;;;;只有两个点时
            datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
            datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
            datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
            if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
              xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
              xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

              statusString = STRING(xd_sState.dataXYZdistance, $
                FORMAT='("Distance = ", F6.2, " mm")')


              demo_putTips, xd_sState, 'locat', 11, /LABEL
              demo_putTips, xd_sState, statusString, 12

              xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
              ;  endif

            endif else begin
              ;;;;;;;;;;;;;;;出现三个点时
              if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

                if(xd_sState.COUNT123  eq 0) THEN begin

                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.COUNT123=1

                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                  result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'
                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+','+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckg->setProperty,hide=1
                  xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 1) THEN begin

                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                    xd_sState.COUNT123=2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'


                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')
                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckb->setProperty,hide=1
                    xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                  endif else begin
                    if(xd_sState.COUNT123  eq 2) THEN begin
                      xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                      xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                      xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                      xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                      xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                      xd_sState.COUNT123=0
                      print,'xd_sState.COUNT123=2'
                      print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                      print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                      xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                      xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                      xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                      result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                      angle=acos(result)*!radeg
                      print,'2 lines angle:',angle,'degree'

                      print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                      statusString = STRING(xd_sState.dataXYZdistance, $
                        FORMAT='("Distance = ", F6.2, " mm")')

                      statusString1 = STRING(angle, $
                        FORMAT='("Angle = ", F6.2, " degree")')

                      statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                        FORMAT='("Distance = ", F6.2, " mm")')


                      xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                      demo_putTips, xd_sState, 'locat', 11, /LABEL
                      demo_putTips, xd_sState, statusString+statusString1, 12
                      xd_sState.ckr1->setProperty,hide=1
                      xd_sState.ckr->setProperty,hide=1
                      xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                      xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                      xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                    endif
                  endelse
                endelse
              endif
            endelse
            xd_sState.oWindow->Draw, xd_sState.oView
            ;draw_3DView, xd_sState
            
            
            ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  
  
            ;ortho2D_plot, xd_sState

  
            ;
            WIDGET_CONTROL, xd_sState.DrawSlider[0], SET_VALUE = location3D[2]
            WIDGET_CONTROL, xd_sState.DrawSlider[1], SET_VALUE = location3D[0]
            WIDGET_CONTROL, xd_sState.DrawSlider[2], SET_VALUE = location3D[1]
  
          endelse
        endif

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end        ;  of DRAW1

    ; 
    ;
    'DRAW2' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      
      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor


       if (sEvent.type EQ 1) then begin
        if (sEvent.release EQ 4) then begin
          ;Left button press
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube = xd_sState.vol_HU_cube
          
          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          
          ;
          xzImage = DBLARR(imageSize, imageSize)
          temp_xzImage = CONGRID(vol_HU_cube[*, location3D[1]*(SIZE(vol_HU_cube))[2]/xd_sState.imageSize, *], imageSize, 1, imageSize)
          xzImage[*, *] = temp_xzImage[*, 0, *]
          dims = SIZE(xzImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube, /DIMENSIONS)

          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE

          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, xzImage, xd_sState.svol_HU_cube, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, ObjectIndex_value=ObjectIndex_value, direction = 'xz',  title='xz plane ROI', REGIONS_OUT = xzROIout, BLOCK=1
          IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252 
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0
              
              
              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              
              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1
            
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN BEGIN
                ;define a new mask
                maskResult = xzROIout[N_ELEMENTS(xzROIout)-1] -> ComputeMask(DIMENSIONS = dims)
                xz_mask_ROI = (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)
                xd_sState.xz_mask_ROI = xz_mask_ROI
                xzROIout[N_ELEMENTS(xzROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
                modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
                print, 'modify roi index = ', modify_roi
                screenSize = xd_sState.imageSize * 3
                Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, xzROIout, modify_roi, ContOrInde, screenSize, direction = 'xz'
                ;not adding red or green objects
                ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi
              ENDIF
            ENDIF
            
            ;replace old yzROIout in xd_sState
            IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN xd_sState.xzROIout = xzROIout[N_ELEMENTS(xzROIout)-1]
            ;draw_3DView, xd_sState
          ENDIF
          
          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
         ;replot the 3 draw windows
          ortho2D_plot, xd_sState

         endif else begin

          ;      print, sEvent.x, sEvent.y
          imageSize = xd_sState.imageSize
          location3D = xd_sState.location3D
          location3D[0] = sEvent.x
          location3D[2] = sEvent.y
          xd_sState.location3D = location3D
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          vol_HU_cube = xd_sState.vol_HU_cube

          if (sEvent.release EQ 1) then begin


            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

 ;           scale = 0.75 + FLOAT(scale) / 100.0
 ;;;;;;;;;;;;;;;;;;;显示点上的三个轴           
            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0
            
            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;;draw_3DView, xd_sState
          endif
  ;;;;;;;;;;;;;;红绿蓝三个点依次显示        
          scale = 0.75 + FLOAT(scale) / 100.0
          s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
          s->getproperty, POLYGONS=conn
          s->getproperty, DATA=verts

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          IF alpha0 EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
          ENDIF ELSE BEGIN
            IF alpha1 EQ 0 THEN BEGIN
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha2 EQ 0 THEN BEGIN
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              ENDIF
            ENDELSE
          ENDELSE

          IF xd_sState.dataXYZind EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, DATA = verts
            xd_sState.oS0->setProperty, POLYGONS = conn
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            xd_sState.dataXYZ0 = dataxyz
            xd_sState.dataXYZind = 1
          ENDIF else begin
            if xd_sState.dataXYZind EQ 1 THEN BEGIN
              xd_sState.oS1->setProperty, DATA = verts
              xd_sState.oS1->setProperty, POLYGONS = conn
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ1 = dataxyz
              xd_sState.dataXYZind = 2
            ENDIF    else begin
              xd_sState.oS2->setProperty, POLYGONS = conn
              xd_sState.oS2->setProperty, DATA = verts
              xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ2 = dataxyz
              xd_sState.dataXYZind = 0
            ENDelse
          endelse

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString, 12

          ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
          ;;;;;;;;;;;;;;;只有两个点时
          datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
          datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
          datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
          if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
            xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
            xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

            statusString = STRING(xd_sState.dataXYZdistance, $
              FORMAT='("Distance = ", F6.2, " mm")')


            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
            ;  endif

          endif else begin
            ;;;;;;;;;;;;;;;出现三个点时
            if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

              if(xd_sState.COUNT123  eq 0) THEN begin

                xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                xd_sState.COUNT123=1

                xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                angle=acos(result)*!radeg
                print,'2 lines angle:',angle,'degree'
                print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                statusString = STRING(xd_sState.dataXYZdistance, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                statusString1 = STRING(angle, $
                  FORMAT='("Angle = ", F6.2, " degree")')
                demo_putTips, xd_sState, 'locat', 11, /LABEL
                demo_putTips, xd_sState, statusString+','+statusString1, 12
                xd_sState.ckr1->setProperty,hide=1
                xd_sState.ckg->setProperty,hide=1
                xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


              endif else begin
                if(xd_sState.COUNT123  eq 1) THEN begin

                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                  xd_sState.COUNT123=2
                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                  result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'


                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckb->setProperty,hide=1
                  xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 2) THEN begin
                    xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                    xd_sState.COUNT123=0
                    print,'xd_sState.COUNT123=2'
                    print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                    print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'

                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')

                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')

                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')


                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckr->setProperty,hide=1
                    xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                  endif
                endelse
              endelse
            endif
          endelse
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          
          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ortho2D_plot, xd_sState
    
          ;
          WIDGET_CONTROL, xd_sState.DrawSlider[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider[2], SET_VALUE = location3D[1]

        endelse
      endif


      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY



    end        ;  of DRAW2

    ;
    'DRAWSLIDER0': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider[0], GET_VALUE=location_value
      
      xd_sState.location3D[2] = location_value - 1
      
      ortho2D_plot, xd_sState
      
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
      
      
    end       ; of DRAWSLIDER2 

    ;
    'DRAWSLIDER1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider[1], GET_VALUE=location_value

      xd_sState.location3D[0] = location_value - 1
      
      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2

    
    'DRAWSLIDER2': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider[2], GET_VALUE=location_value

      xd_sState.location3D[1] = location_value - 1

      ortho2D_plot, xd_sState
      
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2


    ;
    'DRAW0_A1' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN
;      print, 'sEvent.press = ', sEvent.press
;      print, 'sEvent.type = ', sEvent.type

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY


      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a1[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor

      if (sEvent.type EQ 1) then begin


        ;
        if (sEvent.release EQ 4) then begin
          ;Right button press.
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube_a1 = xd_sState.svol_HU_cube_a1 ;GJ 2020/4/18, fix the bug
          
          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube_a1)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          print, 'threshold = ', threshold
          
          xyImage = DBLARR(imageSize, imageSize)
          temp_xyImage = CONGRID(vol_HU_cube_a1[*, *, location3D[2]*(SIZE(vol_HU_cube_a1))[3]/xd_sState.imageSize], imageSize, imageSize, 1)
          xyImage[*, *] = temp_xyImage[*, *, 0]
          dims = SIZE(xyImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube_a1, /DIMENSIONS)

          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, (REVERSE(xyImage, 2)), xd_sState.svol_HU_cube_a1, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, direction = 'xy',  title='xy plane ROI', REGIONS_OUT = xyROIout, BLOCK=1
          IF OBJ_VALID(xyROIout[N_ELEMENTS(xyROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start new 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0


              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713

              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]

              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1

              ;define a new mask
              maskResult = xyROIout[N_ELEMENTS(xyROIout)-1] -> ComputeMask(DIMENSIONS = dims)
              xy_mask_ROI = ((REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0)
              xd_sState.xy_mask_ROI = xy_mask_ROI
              xyROIout[N_ELEMENTS(xyROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
              print, 'ac = ', ac
              modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
              print, 'modify roi index = ', modify_roi
              imageSize = xd_sState.imageSize * 3
              Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, xyROIout, modify_roi, ContOrInde, imageSize, direction = 'xy'
              ;not adding red or green objects
              ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi

              ;replace old xyROIout in xd_sState
              IF OBJ_VALID(xyROIout[N_ELEMENTS(xyROIout)-1]) THEN xd_sState.xyROIout = xyROIout[N_ELEMENTS(xyROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
;          ;draw_3DView, xd_sState

          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

        endif else begin


          imageSize = xd_sState.imageSize

          location3D = xd_sState.location3D

          location3D[0] = sEvent.x

          location3D[1] = imageSize - sEvent.y

          xd_sState.location3D = location3D

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          vol_HU_cube = xd_sState.vol_HU_cube
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          if (sEvent.release EQ 1) then begin

            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube_a1))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube_a1))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube_a1))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0

            ;;;;;;;;;;;;红绿蓝三点依次显示
            scale = 0.75 + FLOAT(scale) / 100.0
            s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
            s->getproperty, POLYGONS=conn
            s->getproperty, DATA=verts

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            IF alpha0 EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha1 EQ 0 THEN BEGIN
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              ENDIF ELSE BEGIN
                IF alpha2 EQ 0 THEN BEGIN
                  xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                ENDIF
              ENDELSE
            ENDELSE

            IF xd_sState.dataXYZind EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, DATA = verts
              xd_sState.oS0->setProperty, POLYGONS = conn
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ0 = dataxyz
              xd_sState.dataXYZind = 1
            ENDIF else begin
              if xd_sState.dataXYZind EQ 1 THEN BEGIN
                xd_sState.oS1->setProperty, DATA = verts
                xd_sState.oS1->setProperty, POLYGONS = conn
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ1 = dataxyz
                xd_sState.dataXYZind = 2
              ENDIF    else begin
                xd_sState.oS2->setProperty, POLYGONS = conn
                xd_sState.oS2->setProperty, DATA = verts
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ2 = dataxyz
                xd_sState.dataXYZind = 0
              ENDelse
            endelse

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
            ;;;;;;;;;;;;;;;只有两个点时
            datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
            datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
            datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
            if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
              xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

              statusString = STRING(xd_sState.dataXYZdistance, $
                FORMAT='("Distance = ", F6.2, " mm")')


              demo_putTips, xd_sState, 'locat', 11, /LABEL
              demo_putTips, xd_sState, statusString, 12

              xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
              ;  endif

            endif else begin
              ;;;;;;;;;;;;;;;出现三个点时
              if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

                if(xd_sState.COUNT123  eq 0) THEN begin

                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.COUNT123=1

                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                  result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'
                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+','+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckg->setProperty,hide=1
                  xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 1) THEN begin

                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                    xd_sState.COUNT123=2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'


                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')
                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckb->setProperty,hide=1
                    xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                  endif else begin
                    if(xd_sState.COUNT123  eq 2) THEN begin
                      xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                      xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                      xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                      xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                      xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                      xd_sState.COUNT123=0
                      print,'xd_sState.COUNT123=2'
                      print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                      print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                      xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                      xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                      xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                      result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                      angle=acos(result)*!radeg
                      print,'2 lines angle:',angle,'degree'

                      print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                      statusString = STRING(xd_sState.dataXYZdistance, $
                        FORMAT='("Distance = ", F6.2, " mm")')

                      statusString1 = STRING(angle, $
                        FORMAT='("Angle = ", F6.2, " degree")')

                      statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                        FORMAT='("Distance = ", F6.2, " mm")')


                      xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                      demo_putTips, xd_sState, 'locat', 11, /LABEL
                      demo_putTips, xd_sState, statusString+statusString1, 12
                      xd_sState.ckr1->setProperty,hide=1
                      xd_sState.ckr->setProperty,hide=1
                      xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                      xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                      xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                    endif
                  endelse
                endelse
              endif
            endelse

            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;draw_3DView, xd_sState
          endif
          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;
          ortho2D_plot, xd_sState

          WIDGET_CONTROL, xd_sState.DrawSlider_a1[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[2], SET_VALUE = location3D[1]

        endelse
      endif

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

    end        ;  of DRAW0

    'DRAW1_A1' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a1[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor


      if (sEvent.type EQ 1) then begin
        if (sEvent.release EQ 4) then begin
          ;Right button press.
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube_a1 = xd_sState.svol_HU_cube_a1

          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube_a1)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          print, 'threshold = ', threshold
          
          yzImage = DBLARR(imageSize, imageSize)
          temp_yzImage = CONGRID(vol_HU_cube_a1[location3D[0]*(SIZE(vol_HU_cube_a1))[1]/xd_sState.imageSize, *, *], 1, imageSize, imageSize)
          yzImage[*, *] = temp_yzImage[0, *, *]
          dims = SIZE(yzImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube_a1, /DIMENSIONS)
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, yzImage, xd_sState.svol_HU_cube_a1, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, direction = 'yz',  title='yz plane ROI', REGIONS_OUT = yzROIout, BLOCK=1
          IF OBJ_VALID(yzROIout[N_ELEMENTS(yzROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start new 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0


              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]

              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1

              ;define a new mask
              maskResult = yzROIout[N_ELEMENTS(yzROIout)-1] -> ComputeMask(DIMENSIONS = dims)
              yz_mask_ROI = (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
              xd_sState.yz_mask_ROI = yz_mask_ROI
              yzROIout[N_ELEMENTS(yzROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
              modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
              print, 'modify roi index = ', modify_roi
              imageSize = xd_sState.imageSize * 3
              Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, yzROIout, modify_roi, ContOrInde, imageSize, direction = 'yz'
              ;not adding red or green objects
              ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi

              ;replace old yzROIout in xd_sState
              IF OBJ_VALID(yzROIout[N_ELEMENTS(yzROIout)-1]) THEN xd_sState.yzROIout = yzROIout[N_ELEMENTS(yzROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
          ;replot the 3 draw windows
          ortho2D_plot, xd_sState
        endif else begin

          ;      print, sEvent.x, sEvent.y
          imageSize = xd_sState.imageSize
          location3D = xd_sState.location3D
          location3D[1] = sEvent.x
          location3D[2] = sEvent.y
          xd_sState.location3D = location3D

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          vol_HU_cube_a1 = xd_sState.vol_HU_cube_a1

          if (sEvent.release EQ 1) then begin


            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube_a1))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube_a1))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube_a1))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

            ;              scale = 0.75 + FLOAT(scale) / 100.0

            ;              s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
            ;              s->getproperty, POLYGONS=conn
            ;              s->getproperty, DATA=verts
            ;              ;          print,'dataxyz=',dataxyz,',conn=',conn,',verts=',verts
            ;              xd_sState.oS3->setProperty, DATA = verts
            ;              xd_sState.oS3->setProperty, POLYGONS = conn
            ;              xd_sState.oS3->setProperty, ALPHA_CHANNEL = 1


            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0


            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;;draw_3DView, xd_sState
          endif
          ;;;;;;;;;;;;;;;;;;;;三个点依次显示
          scale = 0.75 + FLOAT(scale) / 100.0
          s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
          s->getproperty, POLYGONS=conn
          s->getproperty, DATA=verts

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          IF alpha0 EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
          ENDIF ELSE BEGIN
            IF alpha1 EQ 0 THEN BEGIN
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha2 EQ 0 THEN BEGIN
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              ENDIF
            ENDELSE
          ENDELSE

          IF xd_sState.dataXYZind EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, DATA = verts
            xd_sState.oS0->setProperty, POLYGONS = conn
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            xd_sState.dataXYZ0 = dataxyz
            xd_sState.dataXYZind = 1

          ENDIF else begin
            if xd_sState.dataXYZind EQ 1 THEN BEGIN
              xd_sState.oS1->setProperty, DATA = verts
              xd_sState.oS1->setProperty, POLYGONS = conn
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ1 = dataxyz
              xd_sState.dataXYZind = 2

            ENDIF    else begin
              xd_sState.oS2->setProperty, POLYGONS = conn
              xd_sState.oS2->setProperty, DATA = verts
              xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ2 = dataxyz
              xd_sState.dataXYZind = 0

            ENDelse
          endelse

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString, 12

          ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
          ;;;;;;;;;;;;;;;只有两个点时
          datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
          datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
          datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
          if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
            xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
            xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

            statusString = STRING(xd_sState.dataXYZdistance, $
              FORMAT='("Distance = ", F6.2, " mm")')


            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
            ;  endif

          endif else begin
            ;;;;;;;;;;;;;;;出现三个点时
            if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

              if(xd_sState.COUNT123  eq 0) THEN begin

                xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                xd_sState.COUNT123=1

                xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                angle=acos(result)*!radeg
                print,'2 lines angle:',angle,'degree'
                print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                statusString = STRING(xd_sState.dataXYZdistance, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                statusString1 = STRING(angle, $
                  FORMAT='("Angle = ", F6.2, " degree")')
                demo_putTips, xd_sState, 'locat', 11, /LABEL
                demo_putTips, xd_sState, statusString+','+statusString1, 12
                xd_sState.ckr1->setProperty,hide=1
                xd_sState.ckg->setProperty,hide=1
                xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


              endif else begin
                if(xd_sState.COUNT123  eq 1) THEN begin

                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                  xd_sState.COUNT123=2
                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                  result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'


                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckb->setProperty,hide=1
                  xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 2) THEN begin
                    xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                    xd_sState.COUNT123=0
                    print,'xd_sState.COUNT123=2'
                    print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                    print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'

                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')

                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')

                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')


                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckr->setProperty,hide=1
                    xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                  endif
                endelse
              endelse
            endif
          endelse
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState

          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



          ortho2D_plot, xd_sState

          ;
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[2], SET_VALUE = location3D[1]

        endelse
      endif

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end        ;  of DRAW1

    ;
    ;
    'DRAW2_A1' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a1[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor


      if (sEvent.type EQ 1) then begin
        if (sEvent.release EQ 4) then begin
          ;Right button press
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube_a1 = xd_sState.svol_HU_cube_a1   ;GJ 2020/4/18, fix the bug
          ;
          xzImage = DBLARR(imageSize, imageSize)
          temp_xzImage = CONGRID(vol_HU_cube_a1[*, location3D[1]*(SIZE(vol_HU_cube_a1))[2]/xd_sState.imageSize, *], imageSize, 1, imageSize)
          xzImage[*, *] = temp_xzImage[*, 0, *]
          dims = SIZE(xzImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube_a1, /DIMENSIONS)

          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube_a1)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, xzImage, xd_sState.svol_HU_cube_a1, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, direction = 'xz',  title='xz plane ROI', REGIONS_OUT = xzROIout, BLOCK=1
          IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0
  
  
              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
  
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
  
              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1
  
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN BEGIN
                ;define a new mask
                maskResult = xzROIout[N_ELEMENTS(xzROIout)-1] -> ComputeMask(DIMENSIONS = dims)
                xz_mask_ROI = (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)
                xd_sState.xz_mask_ROI = xz_mask_ROI
                xzROIout[N_ELEMENTS(xzROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
                modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
                print, 'modify roi index = ', modify_roi
                screenSize = xd_sState.imageSize * 3
                Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, xzROIout, modify_roi, ContOrInde, screenSize, direction = 'xz'
                ;not adding red or green objects
                ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi
              ENDIF
  
              ;replace old yzROIout in xd_sState
              IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN xd_sState.xzROIout = xzROIout[N_ELEMENTS(xzROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView

          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

        endif else begin

          ;      print, sEvent.x, sEvent.y
          imageSize = xd_sState.imageSize
          location3D = xd_sState.location3D
          location3D[0] = sEvent.x
          location3D[2] = sEvent.y
          xd_sState.location3D = location3D

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          vol_HU_cube_a1 = xd_sState.vol_HU_cube_a1

          if (sEvent.release EQ 1) then begin
          ;Left button

            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube_a1))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube_a1))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube_a1))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

;            scale = 0.75 + FLOAT(scale) / 100.0
            ;;;;;;;;;;;;;;;;;;;显示点上的三个轴
            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0

            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;;draw_3DView, xd_sState
          endif
          ;;;;;;;;;;;;;;红绿蓝三个点依次显示
          scale = 0.75 + FLOAT(scale) / 100.0
          s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
          s->getproperty, POLYGONS=conn
          s->getproperty, DATA=verts

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          IF alpha0 EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
          ENDIF ELSE BEGIN
            IF alpha1 EQ 0 THEN BEGIN
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha2 EQ 0 THEN BEGIN
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              ENDIF
            ENDELSE
          ENDELSE

          IF xd_sState.dataXYZind EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, DATA = verts
            xd_sState.oS0->setProperty, POLYGONS = conn
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            xd_sState.dataXYZ0 = dataxyz
            xd_sState.dataXYZind = 1
          ENDIF else begin
            if xd_sState.dataXYZind EQ 1 THEN BEGIN
              xd_sState.oS1->setProperty, DATA = verts
              xd_sState.oS1->setProperty, POLYGONS = conn
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ1 = dataxyz
              xd_sState.dataXYZind = 2
            ENDIF    else begin
              xd_sState.oS2->setProperty, POLYGONS = conn
              xd_sState.oS2->setProperty, DATA = verts
              xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ2 = dataxyz
              xd_sState.dataXYZind = 0
            ENDelse
          endelse

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString, 12

          ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
          ;;;;;;;;;;;;;;;只有两个点时
          datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
          datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
          datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
          if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
            xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
            xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

            statusString = STRING(xd_sState.dataXYZdistance, $
              FORMAT='("Distance = ", F6.2, " mm")')


            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
            ;  endif

          endif else begin
            ;;;;;;;;;;;;;;;出现三个点时
            if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

              if(xd_sState.COUNT123  eq 0) THEN begin

                xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                xd_sState.COUNT123=1

                xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                angle=acos(result)*!radeg
                print,'2 lines angle:',angle,'degree'
                print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                statusString = STRING(xd_sState.dataXYZdistance, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                statusString1 = STRING(angle, $
                  FORMAT='("Angle = ", F6.2, " degree")')
                demo_putTips, xd_sState, 'locat', 11, /LABEL
                demo_putTips, xd_sState, statusString+','+statusString1, 12
                xd_sState.ckr1->setProperty,hide=1
                xd_sState.ckg->setProperty,hide=1
                xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


              endif else begin
                if(xd_sState.COUNT123  eq 1) THEN begin

                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                  xd_sState.COUNT123=2
                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                  result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'


                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckb->setProperty,hide=1
                  xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 2) THEN begin
                    xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                    xd_sState.COUNT123=0
                    print,'xd_sState.COUNT123=2'
                    print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                    print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'

                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')

                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')

                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')


                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckr->setProperty,hide=1
                    xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                  endif
                endelse
              endelse
            endif
          endelse
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState

          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ortho2D_plot, xd_sState

          ;
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider_a1[2], SET_VALUE = location3D[1]

        endelse
      endif


      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY



    end        ;  of DRAW2

    ;
    'DRAWSLIDER0_A1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider_a1[0], GET_VALUE=location_value

      xd_sState.location3D[2] = location_value - 1

      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2

    ;
    'DRAWSLIDER1_A1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider_a1[1], GET_VALUE=location_value

      xd_sState.location3D[0] = location_value - 1

      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2


    'DRAWSLIDER2_A1': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider_a1[2], GET_VALUE=location_value

      xd_sState.location3D[1] = location_value - 1

      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2
    
    
    
    ;
    'DRAW0_A2' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN
;      print, 'sEvent.press = ', sEvent.press
;      print, 'sEvent.type = ', sEvent.type

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY


      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a2[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor

      if (sEvent.type EQ 1) then begin


        ;
        if (sEvent.release EQ 4) then begin
          ;Right button press.
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube_a2 = xd_sState.svol_HU_cube_a2  ;GJ 2020/4/18, fix the bug
          
          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube_a2)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          
          xyImage = DBLARR(imageSize, imageSize)
          temp_xyImage = CONGRID(vol_HU_cube_a2[*, *, location3D[2]*(SIZE(vol_HU_cube_a2))[3]/xd_sState.imageSize], imageSize, imageSize, 1)
          xyImage[*, *] = temp_xyImage[*, *, 0]
          dims = SIZE(xyImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube_a2, /DIMENSIONS)
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, (REVERSE(xyImage, 2)), xd_sState.svol_HU_cube_a2, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, direction = 'xy',  title='xy plane ROI', REGIONS_OUT = xyROIout, BLOCK=1
          IF OBJ_VALID(xyROIout[N_ELEMENTS(xyROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start new 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0


              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713

              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]

              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1

              ;define a new mask
              maskResult = xyROIout[N_ELEMENTS(xyROIout)-1] -> ComputeMask(DIMENSIONS = dims)
              xy_mask_ROI = ((REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0)
              xd_sState.xy_mask_ROI = xy_mask_ROI
              xyROIout[N_ELEMENTS(xyROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
              print, 'ac = ', ac
              modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
              print, 'modify roi index = ', modify_roi
              imageSize = xd_sState.imageSize * 3
              Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, xyROIout, modify_roi, ContOrInde, imageSize, direction = 'xy'
              ;not adding red or green objects
              ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi

              ;replace old xyROIout in xd_sState
              IF OBJ_VALID(xyROIout[N_ELEMENTS(xyROIout)-1]) THEN xd_sState.xyROIout = xyROIout[N_ELEMENTS(xyROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView

          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

        endif else begin


          imageSize = xd_sState.imageSize

          location3D = xd_sState.location3D

          location3D[0] = sEvent.x

          location3D[1] = imageSize - sEvent.y

          xd_sState.location3D = location3D

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          vol_HU_cube = xd_sState.vol_HU_cube
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          if (sEvent.release EQ 1) then begin

            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0

            ;;;;;;;;;;;;红绿蓝三点依次显示
            scale = 0.75 + FLOAT(scale) / 100.0
            s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
            s->getproperty, POLYGONS=conn
            s->getproperty, DATA=verts

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            IF alpha0 EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha1 EQ 0 THEN BEGIN
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              ENDIF ELSE BEGIN
                IF alpha2 EQ 0 THEN BEGIN
                  xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                ENDIF
              ENDELSE
            ENDELSE

            IF xd_sState.dataXYZind EQ 0 THEN BEGIN
              xd_sState.oS0->setProperty, DATA = verts
              xd_sState.oS0->setProperty, POLYGONS = conn
              xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ0 = dataxyz
              xd_sState.dataXYZind = 1
            ENDIF else begin
              if xd_sState.dataXYZind EQ 1 THEN BEGIN
                xd_sState.oS1->setProperty, DATA = verts
                xd_sState.oS1->setProperty, POLYGONS = conn
                xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ1 = dataxyz
                xd_sState.dataXYZind = 2
              ENDIF    else begin
                xd_sState.oS2->setProperty, POLYGONS = conn
                xd_sState.oS2->setProperty, DATA = verts
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
                xd_sState.dataXYZ2 = dataxyz
                xd_sState.dataXYZind = 0
              ENDelse
            endelse

            xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
            xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
            xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
            ;;;;;;;;;;;;;;;只有两个点时
            datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
            datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
            datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
            if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
              xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

              statusString = STRING(xd_sState.dataXYZdistance, $
                FORMAT='("Distance = ", F6.2, " mm")')


              demo_putTips, xd_sState, 'locat', 11, /LABEL
              demo_putTips, xd_sState, statusString, 12

              xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
              ;  endif

            endif else begin
              ;;;;;;;;;;;;;;;出现三个点时
              if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

                if(xd_sState.COUNT123  eq 0) THEN begin

                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.COUNT123=1

                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                  result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'
                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+','+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckg->setProperty,hide=1
                  xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 1) THEN begin

                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                    xd_sState.COUNT123=2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'


                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')
                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')
                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckb->setProperty,hide=1
                    xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                  endif else begin
                    if(xd_sState.COUNT123  eq 2) THEN begin
                      xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                      xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                      xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                      xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                      xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                      xd_sState.COUNT123=0
                      print,'xd_sState.COUNT123=2'
                      print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                      print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                      xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                      xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                      xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                      result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                      angle=acos(result)*!radeg
                      print,'2 lines angle:',angle,'degree'

                      print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                      statusString = STRING(xd_sState.dataXYZdistance, $
                        FORMAT='("Distance = ", F6.2, " mm")')

                      statusString1 = STRING(angle, $
                        FORMAT='("Angle = ", F6.2, " degree")')

                      statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                        FORMAT='("Distance = ", F6.2, " mm")')


                      xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                      demo_putTips, xd_sState, 'locat', 11, /LABEL
                      demo_putTips, xd_sState, statusString+statusString1, 12
                      xd_sState.ckr1->setProperty,hide=1
                      xd_sState.ckr->setProperty,hide=1
                      xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                      xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                      xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                    endif
                  endelse
                endelse
              endif
            endelse

            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;draw_3DView, xd_sState
          endif
          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;
          ortho2D_plot, xd_sState

          WIDGET_CONTROL, xd_sState.DrawSlider_a2[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[2], SET_VALUE = location3D[1]

        endelse
      endif

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY

    end        ;  of DRAW0

    'DRAW1_A2' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a2[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor


      if (sEvent.type EQ 1) then begin
        if (sEvent.release EQ 4) then begin
          ;Right button press.
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube_a2 = xd_sState.svol_HU_cube_a2 ;GJ 2020/4/18, fix the bug

          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube_a2)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]

          yzImage = DBLARR(imageSize, imageSize)
          temp_yzImage = CONGRID(vol_HU_cube_a2[location3D[0]*(SIZE(vol_HU_cube_a2))[1]/xd_sState.imageSize, *, *], 1, imageSize, imageSize)
          yzImage[*, *] = temp_yzImage[0, *, *]
          dims = SIZE(yzImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube_a2, /DIMENSIONS)
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, yzImage, xd_sState.svol_HU_cube_a2, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, direction = 'yz',  title='yz plane ROI', REGIONS_OUT = yzROIout, BLOCK=1
          IF OBJ_VALID(yzROIout[N_ELEMENTS(yzROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start new 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0


              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]

              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1

              ;define a new mask
              maskResult = yzROIout[N_ELEMENTS(yzROIout)-1] -> ComputeMask(DIMENSIONS = dims)
              yz_mask_ROI = (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
              xd_sState.yz_mask_ROI = yz_mask_ROI
              yzROIout[N_ELEMENTS(yzROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
              modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
              print, 'modify roi index = ', modify_roi
              imageSize = xd_sState.imageSize * 3
              Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, yzROIout, modify_roi, ContOrInde, imageSize, direction = 'yz'
              ;not adding red or green objects
              ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi

              ;replace old yzROIout in xd_sState
              IF OBJ_VALID(yzROIout[N_ELEMENTS(yzROIout)-1]) THEN xd_sState.yzROIout = yzROIout[N_ELEMENTS(yzROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView
          ;replot the 3 draw windows
          ortho2D_plot, xd_sState
        endif else begin

          ;      print, sEvent.x, sEvent.y
          imageSize = xd_sState.imageSize
          location3D = xd_sState.location3D
          location3D[1] = sEvent.x
          location3D[2] = sEvent.y
          xd_sState.location3D = location3D

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          vol_HU_cube = xd_sState.vol_HU_cube

          if (sEvent.release EQ 1) then begin


            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

            ;              scale = 0.75 + FLOAT(scale) / 100.0

            ;              s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
            ;              s->getproperty, POLYGONS=conn
            ;              s->getproperty, DATA=verts
            ;              ;          print,'dataxyz=',dataxyz,',conn=',conn,',verts=',verts
            ;              xd_sState.oS3->setProperty, DATA = verts
            ;              xd_sState.oS3->setProperty, POLYGONS = conn
            ;              xd_sState.oS3->setProperty, ALPHA_CHANNEL = 1


            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0


            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;;draw_3DView, xd_sState
          endif
          ;;;;;;;;;;;;;;;;;;;;三个点依次显示
          scale = 0.75 + FLOAT(scale) / 100.0
          s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
          s->getproperty, POLYGONS=conn
          s->getproperty, DATA=verts

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          IF alpha0 EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
          ENDIF ELSE BEGIN
            IF alpha1 EQ 0 THEN BEGIN
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha2 EQ 0 THEN BEGIN
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              ENDIF
            ENDELSE
          ENDELSE

          IF xd_sState.dataXYZind EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, DATA = verts
            xd_sState.oS0->setProperty, POLYGONS = conn
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            xd_sState.dataXYZ0 = dataxyz
            xd_sState.dataXYZind = 1

          ENDIF else begin
            if xd_sState.dataXYZind EQ 1 THEN BEGIN
              xd_sState.oS1->setProperty, DATA = verts
              xd_sState.oS1->setProperty, POLYGONS = conn
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ1 = dataxyz
              xd_sState.dataXYZind = 2

            ENDIF    else begin
              xd_sState.oS2->setProperty, POLYGONS = conn
              xd_sState.oS2->setProperty, DATA = verts
              xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ2 = dataxyz
              xd_sState.dataXYZind = 0

            ENDelse
          endelse

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString, 12

          ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
          ;;;;;;;;;;;;;;;只有两个点时
          datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
          datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
          datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
          if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
            xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
            xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

            statusString = STRING(xd_sState.dataXYZdistance, $
              FORMAT='("Distance = ", F6.2, " mm")')


            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
            ;  endif

          endif else begin
            ;;;;;;;;;;;;;;;出现三个点时
            if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

              if(xd_sState.COUNT123  eq 0) THEN begin

                xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                xd_sState.COUNT123=1

                xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                angle=acos(result)*!radeg
                print,'2 lines angle:',angle,'degree'
                print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                statusString = STRING(xd_sState.dataXYZdistance, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                statusString1 = STRING(angle, $
                  FORMAT='("Angle = ", F6.2, " degree")')
                demo_putTips, xd_sState, 'locat', 11, /LABEL
                demo_putTips, xd_sState, statusString+','+statusString1, 12
                xd_sState.ckr1->setProperty,hide=1
                xd_sState.ckg->setProperty,hide=1
                xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


              endif else begin
                if(xd_sState.COUNT123  eq 1) THEN begin

                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                  xd_sState.COUNT123=2
                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                  result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'


                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckb->setProperty,hide=1
                  xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 2) THEN begin
                    xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                    xd_sState.COUNT123=0
                    print,'xd_sState.COUNT123=2'
                    print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                    print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'

                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')

                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')

                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')


                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckr->setProperty,hide=1
                    xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                  endif
                endelse
              endelse
            endif
          endelse
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState

          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



          ortho2D_plot, xd_sState

          ;
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[2], SET_VALUE = location3D[1]

        endelse
      endif

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end        ;  of DRAW1

    ;
    ;
    'DRAW2_A2' : begin

      ;  Return if the button  was not the left mouse button or
      ;  the composite view does not exist or the
      ;  sinogram image does not exist.
      ;
      if (sEvent.press NE 0) then RETURN

      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a2[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor


      if (sEvent.type EQ 1) then begin
        if (sEvent.release EQ 4) then begin
          ;Right button press
          location3D = xd_sState.location3D
          imageSize = xd_sState.imageSize * 3
          vol_HU_cube_a2 = xd_sState.svol_HU_cube_a2  ;GJ 2020/4/18, fix the bug
          
          ;calculate the additional image's threshold
          sizeVHc = SIZE(xd_sState.svol_HU_cube_a2)
          mean_image = IMAGE_THRESHOLD(BYTSCL(xd_sState.svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
          image_statistics, xd_sState.svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
          IF tf[0] GT 0.00001 THEN BEGIN
            threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
          ENDIF ELSE BEGIN
            threshold_max_value = img_max
          ENDELSE
          threshold_min_value = img_min
          threshold = [threshold_min_value, threshold_max_value]
          
          xzImage = DBLARR(imageSize, imageSize)
          temp_xzImage = CONGRID(vol_HU_cube_a2[*, location3D[1]*(SIZE(vol_HU_cube_a2))[2]/xd_sState.imageSize, *], imageSize, 1, imageSize)
          xzImage[*, *] = temp_xzImage[*, 0, *]
          dims = SIZE(xzImage, /DIMENSIONS)
          vol_dims = SIZE(vol_HU_cube_a2, /DIMENSIONS)
          
          ;define the mask
          WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
          IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube
            ;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
            ;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
          ENDIF ELSE BEGIN
            mask_vol_cube = xd_sState.new_mask_vol_cube*0
            mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
          ENDELSE
          
          ;draw new ROI
          ;GJ 2020/4/11, displaying the image for better contrast
          XD_XROI, xzImage, xd_sState.svol_HU_cube_a2, threshold, mask_vol_cube, xd_sState.location3D, xd_sState.imageSize, xd_sState.icon_dir, direction = 'xz',  title='xz plane', REGIONS_OUT = xzROIout, BLOCK=1
          IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN BEGIN
            result_question = DIALOG_MESSAGE('Start 3D reconstruction?', /Question, TITLE='Continue...')
            IF result_question EQ 'Yes' THEN BEGIN
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              s = obj_new('orb',radius=0.1,shading=1)
              s->getproperty, POLYGONS=conn
              s->getproperty, DATA=verts
              xd_sState.dataXYZ0=fltarr(3)
              xd_sState.dataXYZ1=fltarr(3)
              xd_sState.dataXYZ2=fltarr(3)
              xd_sState.COUNT123=0
              xd_sState.dataXYZind=0
  
  
              xd_sState.oS0->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS0->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS1->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS1->setProperty, POLYGONS = fltarr(n_elements(conn))
              xd_sState.oS2->setProperty, DATA = fltarr(3,n_elements(verts)/3)
              xd_sState.oS2->setProperty, POLYGONS = fltarr(n_elements(conn))
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
              xd_sState.ogreenaxis1->setProperty,  hide=1
              xd_sState.oyellowaxis1->setProperty, hide=1
              xd_sState.oblueAxis1->setProperty,  hide=1
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
  
              xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
              xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
  
              xd_sState.ckr1->setProperty,hide=1
              xd_sState.ckg->setProperty,hide=1
              xd_sState.ckr->setProperty, hide=1
              xd_sState.ckb->setProperty, hide=1
              xd_sState.ckag->setProperty,hide=1
  
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
              IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN BEGIN
                ;define a new mask
                maskResult = xzROIout[N_ELEMENTS(xzROIout)-1] -> ComputeMask(DIMENSIONS = dims)
                xz_mask_ROI = (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)
                xd_sState.xz_mask_ROI = xz_mask_ROI
                xzROIout[N_ELEMENTS(xzROIout)-1] -> getProperty, ALPHA_CHANNEL=ac, LINESTYLE = ContOrInde
                modify_roi = CEIL((ac * 10000. - FLOOR(ac * 10000.)) * 100)   ;1=border, 2=add, 3=remove
                print, 'modify roi index = ', modify_roi
                screenSize = xd_sState.imageSize * 3
                Generate_Mask_ROI_modify, xd_sState, mask_vol_cube, xzROIout, modify_roi, ContOrInde, screenSize, direction = 'xz'
                ;not adding red or green objects
                ;IF ABS(modify_roi) LE 10 THEN Mask_segmentation, xd_sState, modify_roi
              ENDIF
  
              ;replace old yzROIout in xd_sState
              IF OBJ_VALID(xzROIout[N_ELEMENTS(xzROIout)-1]) THEN xd_sState.xzROIout = xzROIout[N_ELEMENTS(xzROIout)-1]
              ;draw_3DView, xd_sState
            ENDIF
          ENDIF

          WIDGET_CONTROL, sEvent.top, /HOURGLASS
          xd_sState.oWindow->Draw, xd_sState.oView

          ;replot the 3 draw windows
          ortho2D_plot, xd_sState

        endif else begin

          ;      print, sEvent.x, sEvent.y
          imageSize = xd_sState.imageSize
          location3D = xd_sState.location3D
          location3D[0] = sEvent.x
          location3D[2] = sEvent.y
          xd_sState.location3D = location3D

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          vol_HU_cube_a2 = xd_sState.vol_HU_cube_a2

          if (sEvent.release EQ 1) then begin


            dataxyz = DBLARR(3)
            dataxyz[0] = (location3D[0]*((SIZE(xd_sState.vol_HU_cube_a2))[1]))/xd_sState.imageSize
            dataxyz[1] = (location3D[1]*((SIZE(xd_sState.vol_HU_cube_a2))[2]))/xd_sState.imageSize
            dataxyz[2] = (location3D[2]*((SIZE(xd_sState.vol_HU_cube_a2))[3]))/xd_sState.imageSize
            xd_sState.location3D = location3D

;            scale = 0.75 + FLOAT(scale) / 100.0
            ;;;;;;;;;;;;;;;;;;;显示点上的三个轴
            s1 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
            s2 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
            s3 = OBJ_NEW('IDLgrPolyline', $
              [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])


            s1->getproperty, POLYLINES=conn1
            s1->getproperty, DATA=verts1

            s2->getproperty, POLYLINES=conn2
            s2->getproperty, DATA=verts2

            s3->getproperty, POLYLINES=conn3
            s3->getproperty, DATA=verts3

            xd_sState.ogreenaxis1->setProperty, DATA = verts1
            xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
            xd_sState.ogreenaxis1->SetProperty, HIDE=0

            xd_sState.oyellowaxis1->setProperty, DATA = verts2
            xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
            xd_sState.oyellowaxis1->SetProperty, HIDE=0

            xd_sState.oBlueAxis1->setProperty, DATA = verts3
            xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
            xd_sState.oBlueAxis1->SetProperty, HIDE=0

            WIDGET_CONTROL, sEvent.top, /HOURGLASS
            xd_sState.oWindow->Draw, xd_sState.oView
            ;;draw_3DView, xd_sState
          endif
          ;;;;;;;;;;;;;;红绿蓝三个点依次显示
          scale = 0.75 + FLOAT(scale) / 100.0
          s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
          s->getproperty, POLYGONS=conn
          s->getproperty, DATA=verts

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          IF alpha0 EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
          ENDIF ELSE BEGIN
            IF alpha1 EQ 0 THEN BEGIN
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
            ENDIF ELSE BEGIN
              IF alpha2 EQ 0 THEN BEGIN
                xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              ENDIF
            ENDELSE
          ENDELSE

          IF xd_sState.dataXYZind EQ 0 THEN BEGIN
            xd_sState.oS0->setProperty, DATA = verts
            xd_sState.oS0->setProperty, POLYGONS = conn
            xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1
            xd_sState.dataXYZ0 = dataxyz
            xd_sState.dataXYZind = 1
          ENDIF else begin
            if xd_sState.dataXYZind EQ 1 THEN BEGIN
              xd_sState.oS1->setProperty, DATA = verts
              xd_sState.oS1->setProperty, POLYGONS = conn
              xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ1 = dataxyz
              xd_sState.dataXYZind = 2
            ENDIF    else begin
              xd_sState.oS2->setProperty, POLYGONS = conn
              xd_sState.oS2->setProperty, DATA = verts
              xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1
              xd_sState.dataXYZ2 = dataxyz
              xd_sState.dataXYZind = 0
            ENDelse
          endelse

          xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
          xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
          xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString, 12

          ;;;;;;;;;;;;;;;;;;;;;;;红绿蓝三点连线，并测量距离和角度
          ;;;;;;;;;;;;;;;只有两个点时
          datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
          datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
          datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
          if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
            xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1
            xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

            statusString = STRING(xd_sState.dataXYZdistance, $
              FORMAT='("Distance = ", F6.2, " mm")')


            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString, 12

            xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
            ;  endif

          endif else begin
            ;;;;;;;;;;;;;;;出现三个点时
            if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN

              if(xd_sState.COUNT123  eq 0) THEN begin

                xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0.
                xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                xd_sState.COUNT123=1

                xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
                result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                angle=acos(result)*!radeg
                print,'2 lines angle:',angle,'degree'
                print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                statusString = STRING(xd_sState.dataXYZdistance, $
                  FORMAT='("Distance = ", F6.2, " mm")')
                xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
                statusString1 = STRING(angle, $
                  FORMAT='("Angle = ", F6.2, " degree")')
                demo_putTips, xd_sState, 'locat', 11, /LABEL
                demo_putTips, xd_sState, statusString+','+statusString1, 12
                xd_sState.ckr1->setProperty,hide=1
                xd_sState.ckg->setProperty,hide=1
                xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
                xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
                xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0


              endif else begin
                if(xd_sState.COUNT123  eq 1) THEN begin

                  xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0.
                  xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3   ;dotted line
                  xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
                  xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                  xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
                  xd_sState.COUNT123=2
                  xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

                  xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
                  xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                  result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
                  angle=acos(result)*!radeg
                  print,'2 lines angle:',angle,'degree'


                  print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                  statusString = STRING(xd_sState.dataXYZdistance, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  statusString1 = STRING(angle, $
                    FORMAT='("Angle = ", F6.2, " degree")')
                  statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                    FORMAT='("Distance = ", F6.2, " mm")')
                  xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                  demo_putTips, xd_sState, 'locat', 11, /LABEL
                  demo_putTips, xd_sState, statusString+statusString1, 12
                  xd_sState.ckr1->setProperty,hide=1
                  xd_sState.ckb->setProperty,hide=1
                  xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
                  xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
                  xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0


                endif else begin
                  if(xd_sState.COUNT123  eq 2) THEN begin
                    xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
                    xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0  ;dotted line
                    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3

                    xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0
                    xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
                    xd_sState.COUNT123=0
                    print,'xd_sState.COUNT123=2'
                    print,'xd_sState.dataXYZ0',xd_sState.dataXYZ0
                    print,'xd_sState.dataXYZ2',xd_sState.dataXYZ2
                    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution

                    xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
                    xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

                    result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
                    angle=acos(result)*!radeg
                    print,'2 lines angle:',angle,'degree'

                    print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
                    statusString = STRING(xd_sState.dataXYZdistance, $
                      FORMAT='("Distance = ", F6.2, " mm")')

                    statusString1 = STRING(angle, $
                      FORMAT='("Angle = ", F6.2, " degree")')

                    statusString_old = STRING(xd_sState.dataXYZdistance_old, $
                      FORMAT='("Distance = ", F6.2, " mm")')


                    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

                    demo_putTips, xd_sState, 'locat', 11, /LABEL
                    demo_putTips, xd_sState, statusString+statusString1, 12
                    xd_sState.ckr1->setProperty,hide=1
                    xd_sState.ckr->setProperty,hide=1
                    xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
                    xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
                    xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0

                  endif
                endelse
              endelse
            endif
          endelse
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState

          ;   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          ortho2D_plot, xd_sState

          ;
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[0], SET_VALUE = location3D[2]
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[1], SET_VALUE = location3D[0]
          WIDGET_CONTROL, xd_sState.DrawSlider_a2[2], SET_VALUE = location3D[1]

        endelse
      endif


      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY



    end        ;  of DRAW2

    ;
    'DRAWSLIDER0_A2': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider_A2[0], GET_VALUE=location_value

      xd_sState.location3D[2] = location_value - 1

      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2

    ;
    'DRAWSLIDER1_A2': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider_a2[1], GET_VALUE=location_value

      xd_sState.location3D[0] = location_value - 1

      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2


    'DRAWSLIDER2_A2': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
      WIDGET_CONTROL, sEvent.top, /HOURGLASS
      WIDGET_CONTROL, xd_sState.DrawSlider_a2[2], GET_VALUE=location_value

      xd_sState.location3D[1] = location_value - 1

      ortho2D_plot, xd_sState

      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY


    end       ; of DRAWSLIDER2


    ;  Handle the event that occurs in the drawing (viewing) area.
    ;  These are : expose and mouse button(press, motion, release).
    ;
    'DRAW3': begin
      WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY

      ;  Expose.
      ;
      if (sEvent.type eq 4) then begin
        xd_sState.oWindow->draw, xd_sState.oView
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      endif

      ;  User rotation.
      ;
      void = xd_sState.oTraceRotationModel->Update(sEvent, mouse=4)
      
      if xd_sState.oRotationModel->Update(sEvent) then xd_sState.oWindow->Draw, xd_sState.oView
      if xd_sState.oRotationModel1->Update(sEvent) then xd_sState.oWindow->Draw, xd_sState.oView 

      if (sEvent.type EQ 0) then begin
        ;Left button press.
        if (sEvent.press EQ 1) then begin
          pick = xd_sState.oWindow->PickData(xd_sState.oView, $
            xd_sState.oPoly[0], [sEvent.x, sEvent.y], dataxyz)
          picked = xd_sState.oWindow->select(xd_sState.oView, [sEvent.x,sEvent.y])
          IF ~OBJ_VALID(picked[0]) THEN BEGIN
            IF picked[0] EQ -1 THEN BEGIN
              WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
              RETURN
            ENDIF
          ENDIF
          xd_sState.selected = picked[0]
          print, 'picked = ', xd_sState.selected
          print, 'mouse location = ', [sEvent.x, sEvent.y]
          ;GJ 2019/8/16
          print, '3d location = ', dataxyz
          xd_sState.seed_location = dataxyz
          if (pick ne 0) then begin
            xd_sState.selected->GetProperty, NAME=name_poly, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel, hide=hideStatus, THICK=thickness
            IF N_ELEMENTS(vert_colors) EQ 1 THEN BEGIN
              WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
              RETURN
            ENDIF
            
            ;set the alpha-channel transparent value
            WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
            WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')
            
            ;thick was used to record smoothness
            smoothness = (thickness-FLOOR(thickness))*1000
            WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=smoothness
            WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothness, format='(f5.1)')
            
            ;set the object index correctly
            WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=LONG(name_poly)+1
            WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(LONG(name_poly)+1, format='(I3)')
            
            ;Enable color selection
            widget_control, xd_sState.wColorIndex, sensitive=1
            widget_control, xd_sState.wSmooth, sensitive=1
            widget_control, xd_sState.wSmooth_text, sensitive=1
            widget_control, xd_sState.wtransSlider, sensitive=1
            widget_control, xd_sState.wTrans_text, sensitive=1
            widget_control, xd_sState.w3DSegd, sensitive=1
            widget_control, xd_sState.w3Drg, sensitive=1
            widget_control, xd_sState.wthickness, sensitive=1
            widget_control, xd_sState.wpoint, sensitive=1
            widget_control, xd_sState.wwire, sensitive=1
            widget_control, xd_sState.wfill, sensitive=1
            widget_control, xd_sState.wexpstl, sensitive=1
;            widget_control, xd_sState.whide, sensitive=1
            widget_control, xd_sState.wsinglecenter, sensitive=1
            
            ;the object is multi-color, we don't need to do anything
            ;if the object is single color, we will find right color index
            IF STDDEV(vert_colors[0,*]) LT 1 THEN BEGIN
              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
              FOR m=0, N_ELEMENTS(ColorIndex)-1 DO BEGIN
                IF vert_colors[0,0] EQ color_rgb[m,0] AND vert_colors[1,0] EQ color_rgb[m,1] AND vert_colors[2,0] EQ color_rgb[m,2] THEN colorindexS = m
              ENDFOR
              IF N_ELEMENTS(colorindexS) EQ 0 THEN colorindexS = 0
              WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=colorindexS
              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              ERASE, ColorIndexBGR[colorindexS]
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDIF ELSE BEGIN

              WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
              DEVICE, DECOMPOSED=1
              WSET, wColorIndexID
              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
              ERASE, ColorIndexBGR[LONG(name_poly)+1]
              ;use multiple lines to label multiple colro
              plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
              plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
              plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device
              IF hideStatus EQ 1 THEN BEGIN
                plots, [0,40], [0,30], color=[0,0,0], /device
                plots, [0,40], [30,0], color=[0,0,0], /device
              ENDIF
            ENDELSE
            
            xd_sState.tvol_white.SetProperty,hide=0
            count = xd_sState.struc_modifyobject.countmask_list[LONG(name_poly)]
            xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
            Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
;            IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
            print, 'vert_colors[*,0] = ', vert_colors[*,0]
            IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
            
              ;GJ 2019/2/22, summarize all plot pickdata into a program        
            PLOT_pickdata, dataxyz, xd_sState                                   
;            wZGG_sensi = WIDGET_INFO(xd_sState.wZGG, /sensitive)
;            wCKcomplete_sensi = WIDGET_INFO(xd_sState.wCKcomplete, /sensitive)
;            wCKcompleteNot_sensi = WIDGET_INFO(xd_sState.wCKcompleteNot, /sensitive)
;            wZB_sensi = WIDGET_INFO(xd_sState.wZB, /sensitive)
;            IF wZGG_sensi EQ 0 AND wCKcomplete_sensi EQ 0 AND wCKcompleteNot_sensi EQ 0 AND wZB_sensi EQ 0 THEN WIDGET_CONTROL, xd_sState.wZGG, SENSITIVE = 1
;            IF wZGG_sensi EQ 0 AND wCKcomplete_sensi EQ 1 AND wCKcompleteNot_sensi EQ 1 AND wZB_sensi EQ 0 THEN WIDGET_CONTROL, xd_sState.wZB, SENSITIVE = 1
            xd_sState.oWindow->Draw, xd_sState.oView
               
            if (xd_sState.tracingMode eq 1) then begin
              d_surfviewAddTracePoint, xd_sState.oTracePolyline, $
                dataxyz, check
              xd_sState.firstPoint = dataxyz
              if (check eq -1) then begin
                PRINT, $
                 'Error in d_surfviewAddTracePoint/Draw.'
              endif
              xd_sState.oWindow->Draw, xd_sState.oView
            endif
            xd_sState.btndown = 1b
            WIDGET_CONTROL, xd_sState.wDraw, /DRAW_MOTION
          endif else begin
            demo_putTips, xd_sState, "Data point:In background", 12
          endelse
        endif else begin
          ;not left mouse button
          xd_sState.btndown = 4b
          xd_sState.oWindow->SetProperty, QUALITY=xd_sState.dragq
          WIDGET_CONTROL, xd_sState.wDraw, /DRAW_MOTION
        endelse
      endif    ; ev.typ EQ 0

      ;  Button motion and tracing mode
      if ((sEvent.type eq 2) and (xd_sState.btndown eq 1b)) then begin
        pick = xd_sState.oWindow->PickData(xd_sState.oView, $
          xd_sState.oPoly[0], [sEvent.x, sEvent.y], dataxyz)
        if (pick ne 0) then begin
          statusString = STRING(dataxyz[0], $
            dataxyz[1],dataxyz[2], $
            FORMAT='("X=", F6.2,' + $
            ' ", Y=",F6.2,", Z=",F6.2)')
          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString, 12
         
          if (xd_sState.tracingMode ne 0) then begin
            d_surfviewAddTracePoint, xd_sState.oTracePolyline, $
              dataxyz, check
            if (check eq -1) then begin
              PRINT,'Error in d_surfviewAddTracePoint/Draw.'
            endif
            xd_sState.oWindow->Draw, xd_sState.oView
          endif
        endif else begin
          demo_putTips, xd_sState, "Data point:In background", 12
        endelse
      endif

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;zxh
      ;wheel event, zhang xianghuai
      if (sEvent.type eq 7) then begin
        print, 'wheel test1'
        if (sEvent.clicks lt 0) then begin
          
          print, 'wheel test2'
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          scale = ABS(FLOAT(SCALE)/100 - FLOAT(10)/ 100.0) ;ABS is used to avoid negative scale, GJ 2019/2/7
          scalep = scale*100.0
          scalingString = 'Scaling : ' + STRING(scalep, $
            FORMAT='(f5.1)') + ' %'
;          WIDGET_CONTROL, xd_sState.wScalingLabel, $
;            SET_VALUE=scalingString
          WIDGET_CONTROL, xd_sState.wScalingSlider, SET_VALUE=scalep
          transform = [[scale, 0, 0, 0.0], [0, scale, 0, 0.0], $
            [0, 0, scale, 0.0], [0, 0, 0, 1]]
          xd_sState.oScalingModel->SetProperty, TRANSFORM=transform
          xd_sState.oRuler1model->setproperty,transform=[[scale, 0, 0, 0.0],[0, 1, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
          xd_sState.oRuler2model->setproperty,transform=[[1, 0, 0, 0.0],[0, scale, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
          xd_sState.oRuler1Model->scale,0.008,0.05,1
          xd_sState.oRuler1Model->translate,-0.2,0.4,0.45
          xd_sState.oRuler2Model->scale,0.05,0.008,1
          xd_sState.oRuler2Model->translate,0.4,-0.2,0.45
          xd_sState.oTraceScalingModel->SetProperty, TRANSFORM=transform
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
        endif else begin
          WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
          scale = ABS(FLOAT(SCALE)/100 + FLOAT(10)/ 100.0)
          scalep = scale*100.0
          scalingString = 'Scaling : ' + STRING(scalep, $
            FORMAT='(f5.1)') + ' %'
;          WIDGET_CONTROL, xd_sState.wScalingLabel, $
;            SET_VALUE=scalingString
          WIDGET_CONTROL, xd_sState.wScalingSlider, SET_VALUE=scalep
          transform = [[scale, 0, 0, 0.0], [0, scale, 0, 0.0], $
            [0, 0, scale, 0.0], [0, 0, 0, 1]]
          xd_sState.oScalingModel->SetProperty, TRANSFORM=transform
          xd_sState.oRuler1model->setproperty,transform=[[scale, 0, 0, 0.0],[0, 1, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
          xd_sState.oRuler2model->setproperty,transform=[[1, 0, 0, 0.0],[0, scale, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
          xd_sState.oRuler1Model->scale,0.008,0.05,1
          xd_sState.oRuler1Model->translate,-0.2,0.4,0.45
          xd_sState.oRuler2Model->scale,0.05,0.008,1
          xd_sState.oRuler2Model->translate,0.4,-0.2,0.45
          xd_sState.oTraceScalingModel->SetProperty, TRANSFORM=transform
          xd_sState.oWindow->Draw, xd_sState.oView
          ;draw_3DView, xd_sState
          ortho2D_plot, xd_sState
        endelse
      endif
      
      ;Button release.
      if (sEvent.type eq 1) then begin
        ;draw_3DView, xd_sState
        ortho2D_plot, xd_sState
        if (xd_sState.btndown EQ 4b) then begin
          xd_sState.oWindow->SetProperty, QUALITY=2 ; High
          if xd_sState.dragq ne 2 then begin
            xd_sstate.opoly[0]->GetProperty, STYLE=style
            if style eq 5 then $ ; Lego
              widget_control, /hourglass
            if style eq 6 then $ ; Filled Lego
              widget_control, /hourglass
            xd_sState.oWindow->Draw, xd_sState.oView
          endif
        endif else if ((xd_sState.btndown EQ 1b) $
          AND (xd_sState.tracingmode EQ 1) ) then begin
          d_surfviewAddTracePoint, xd_sState.oTracePolyline, $
            xd_sState.firstPoint, check
          if (check eq -1) then begin
            PRINT,'Error in d_surfviewAddTracePoint/Draw.'
          endif
          xd_sState.oWindow->Draw, xd_sState.oView
          xd_sState.tracingMode = 1
          xd_sState.oTracePolyline->GetProperty, POLYLINE=polyline
          if polyline[0] ge 3 then begin
            WIDGET_CONTROL, xd_sState.wTracingMaskButton, SENSITIVE=1
            demo_putTips, xd_sState, $
              ['mask1','displ'], [11,12], /LABEL
            if d_surfviewToggleOffOn(xd_sState.wTracingMaskButton) $
              eq 1 $
              then $
              void = d_surfviewToggleOffOn( $
              xd_sState.wTracingMaskButton $
              )
          endif
        endif
        xd_sState.btndown = 1b
        WIDGET_CONTROL, xd_sState.wDraw, DRAW_MOTION = 0
      endif
      
      WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
    end   ; of DRAW

    ;  Quit the application.
    ;
    'QUIT' : begin
      WIDGET_CONTROL, sEvent.top, /DESTROY
      RETURN
    end   ; of QUIT

  endcase   ; of case uval of
ENDIF

  if XREGISTERED('demo_tour') eq 0 then begin
    WIDGET_CONTROL, sEvent.top, GET_UVALUE=xd_sState, /NO_COPY
    WIDGET_CONTROL, xd_sState.wHotKeyReceptor, /INPUT_FOCUS
    WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
  end
end     ; of event handler


;----------------------------------------------------------------------------
;
;  Purpose:  Restore the previous color table and
;            destroy the top objects.
;
pro d_surfview1Cleanup, $
  wTopBase          ;  IN: top level base ID.
  WIDGET_CONTROL, wTopBase, GET_UVALUE=xd_sState, /NO_COPY

  ;Clean up heap variables.
  ;
  for i=0,n_tags(xd_sState)-1 do begin
    case size((xd_sState).(i), /TNAME) of
      'POINTER': $
        ptr_free, (xd_sState).(i)
      'OBJREF': $
        obj_destroy, (xd_sState).(i)
      else:
    endcase
  end

  ;  Silently flush any accumulated math errors.
  ;
  void = check_math()

  ;  Restore math error behavior.
  ;
  !except = xd_sState.orig_except

  ;  Restore the color table
  ;
  TVLCT, xd_sState.colorTable

  ;  Map the group leader base if it exists.
  ;
  if (WIDGET_INFO(xd_sState.groupBase, /VALID_ID)) then $
    WIDGET_CONTROL, xd_sState.groupBase, /MAP

end   ;  of d_surfview1Cleanup



;----------------------------------------------------------------------------
;
;  Purpose:  Display a surface.
;
pro d_surfview1,sel_one,vol_HU_cube, vol_HU_cube_mask, vol_HU_cube_mask_modify, vol_HU_cube_resolution, struc_modifyobject, cur_patient, icon_dir, additionalDir = additionalDir, $  ;sel_one ---->閻炴稏鍔庨妵姘辨偖椤愶讣鎷烽柟灏佹櫊閸庢挳宕氶崱娆愮暠濞ｅ洠鍓濇导鍛存晬鐏炶棄寰斿ù锝嗘尵濠�拷MODIFY'闁告繂绉寸花鍙夌鐎ｂ晜顐藉ù鐙呯悼閻栵拷
  ;    ALT_FUNC=ALT_FUNC, $       ; IN: (opt) Alternative function : sine dist
  TRANSPARENT=transparent, $ ; IN: (opt) Transparent across a plane
  GROUP=group, $             ; IN: (opt) group identifier
  DEBUG=debug, $             ; IN: (opt) debug mode
  RECORD_TO_FILENAME=record_to_filename, $
  APPTLB = appTLB            ; OUT: (opt) TLB of this application


  ; Check the validity of the group identifier
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

 
; lqrname, Outverts, Outconn
 
 
 
 
  ;  Set up dimensions of the drawing (viewing) area.
  ;
  device, GET_SCREEN_SIZE=scr
  xdim = scr[0]*0.6
  ydim = xdim*0.8
  
  ;print, 'xdim = ', xdim

  ;oPoly = OBJ_NEW('IDLgrPolygon', sel_one.verts, POLYGON=sel_one.conn, STYLE=2, $
  ;  SHADING=1, vert_colors=vert_colors, $
  ;  XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV)
  
  ;add a cube with the size of verts ;GJ, 2018/6/09
  oPoly_size = [MAX(sel_one.verts[0,*])-MIN(sel_one.verts[0,*]), MAX(sel_one.verts[1,*])-MIN(sel_one.verts[1,*]), MAX(sel_one.verts[2,*])-MIN(sel_one.verts[2,*])]
  print, 'oPoly_size = ', oPoly_size
  print, 'oPoly in mm = ', oPoly_size*vol_HU_cube_resolution
  ;struc_modifyobject[cur_patient].hide_Status = 1
  ;print, 'hide status = ', struc_modifyobject[cur_patient].hide_Status
  print, 'struc_modifyobject.modify_step = ', struc_modifyobject[cur_patient].modify_step
  print, 'struc_modifyobject.exist_list = ', struc_modifyobject[cur_patient].exist_num[struc_modifyobject[cur_patient].modify_step]
  print, 'struc_modifyobject.countmask_list = ', struc_modifyobject[cur_patient].countmask_list[0:10]
;  print, 'struc_modifyobject.verts_min_list = ', struc_modifyobject[cur_patient].verts_min_list[0:10, *]
;  print, 'struc_modifyobject.verts_max_list = ', struc_modifyobject[cur_patient].verts_max_list[0:10, *]
  oPoly_min = [MIN(sel_one.verts[0,*]), MIN(sel_one.verts[1,*]), MIN(sel_one.verts[2,*])]
  oPoly_max = [MAX(sel_one.verts[0,*]), MAX(sel_one.verts[1,*]), MAX(sel_one.verts[2,*])]
  verts = [[MIN(sel_one.verts[0,*]), MIN(sel_one.verts[1,*]), MIN(sel_one.verts[2,*])], $
           [MAX(sel_one.verts[0,*]), MIN(sel_one.verts[1,*]), MIN(sel_one.verts[2,*])], $
           [MAX(sel_one.verts[0,*]), MAX(sel_one.verts[1,*]), MIN(sel_one.verts[2,*])], $
           [MIN(sel_one.verts[0,*]), MAX(sel_one.verts[1,*]), MIN(sel_one.verts[2,*])], $
           [MIN(sel_one.verts[0,*]), MIN(sel_one.verts[1,*]), MAX(sel_one.verts[2,*])], $
           [MAX(sel_one.verts[0,*]), MIN(sel_one.verts[1,*]), MAX(sel_one.verts[2,*])], $
           [MAX(sel_one.verts[0,*]), MAX(sel_one.verts[1,*]), MAX(sel_one.verts[2,*])], $
           [MIN(sel_one.verts[0,*]), MAX(sel_one.verts[1,*]), MAX(sel_one.verts[2,*])]]
 ;          verts = [[-0.1,-0.1,-0.1],[0.1,-0.1,-0.1],[0.1,0.1,-0.1], $
 ;          [-0.1,0.1,-0.1], $
 ;          [-0.1,-0.1,0.1],[0.1,-0.1,0.1],[0.1,0.1,0.1],$
 ;          [-0.1,0.1,0.1]]
  conn = [[4,3,2,1,0],[4,4,5,6,7],[4,0,1,5,4],[4,1,2,6,5], [4,2,3,7,6],[4,3,0,4,7]]
  oPoly_border = obj_new('IDLgrPolygon',verts,poly=conn,color=[255,255,255], STYLE=0, $
    SHADING=1, ALPHA_CHANNEL = 0.01, $
    XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV)
    ;end GJ

  oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR =[255,127,127], sel_one.verts, POLYGON=sel_one.conn, STYLE=2, $
    SHADING=1, $
    XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV)

  ;  oPoly->SetProperty,CLIP_PLANES=[-70, -108, -0, (70*70 +108*108 + 0)]
  oPoly1->SetProperty,CLIP_PLANES=[1, 1, 1, -(3)]


  oPoly2 = OBJ_NEW('IDLgrPolygon', COLOR =[1,255,25], sel_one.verts, POLYGON=sel_one.conn, STYLE=2, $
    SHADING=1, $
    XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV)

  ;  oPoly->SetProperty,CLIP_PLANES=[-70, -108, -0, (70*70 +108*108 + 0)]
  oPoly2->SetProperty,CLIP_PLANES=[1, 1, 1, -(3)]

  oRotationModel = OBJ_NEW('IDLexRotator',  [xdim/2.0, ydim/2.0], ydim/2.0)
  ;/////////////////
  oRotationModel1 = OBJ_NEW('IDLexRotator', [xdim/2.0, ydim/2.0], ydim/2.0)
;  oPoly->getProperty,  XCOORD_CONV= XCOORD_CONV
;  print,' XCOORD_CONV11', XCOORD_CONV
  oRotationModel1->Add, oPoly1
  oRotationModel1->Add, oPoly2
  
  ;GJ 2019/8/17; add a surface
  ;
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
  ;GJ 2019/10/10, for general purpose
  imagefile=icon_dir+'XidianUniv.jpg';kelly.jpg';lucy2.jpg';kelly.jpg'; gj.jpg
  read_jpeg, imagefile, idata1

  idata1_0 = CONGRID(idata1, 3, 512, 512)
  FOR i=0,511 DO BEGIN
    FOR j=0, 511 DO BEGIN
      idata[i, j, *] = idata1_0[*, i, j]
    ENDFOR
  ENDFOR

  oImage1 = OBJ_NEW('IDLgrImage', idata, INTERLEAVE=2)

  ;  Create the surface object.
  ;
  oSurface = OBJ_NEW('IDLgrSurface', z, $
    STYLE=2, $
    SHADING=1, $
    HIDE=1, $
    COLOR=[255,255,255], $
    TEXTURE_MAP=oImage1, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV)
  ;GJ 2019/8/17
  oRotationModel->Add, oSurface
    
  ;/////////////////
  n_Polys = struc_modifyobject[cur_patient].exist_num[struc_modifyobject[cur_patient].modify_step]
  ;make the total number of 30 objects
  oPoly = MAKE_ARRAY(30, /OBJ)
  
  ;Change the draw order, draw the smallest one first, draw the largest one finally
  ;so the transparent effect for the large object could be guaranted
  FOR i= n_Polys-1, 0, -1 DO BEGIN
    oPoly[i] = *(struc_modifyobject[cur_patient].polygon_list[i])
    oRotationModel->Add, oPoly[i]
    oPoly[i]->getProperty, NAME=name_poly
    print, 'poly_name = ', LONG(name_poly)
    IF i EQ 0 THEN BEGIN
      oPoly[i]->getProperty, vert_colors=vert_colors
      print, 'vert_colors = ', vert_colors[*,0]
    ENDIF
  ENDFOR
  ;GJ, 2018/6/09
  oRotationModel->Add, oPoly_border
  ;end GJ 2019/2/21
      
  v2=sel_one.verts
  ma=max(v2)
  mi=min(v2)
  bias=(ma-mi)/2
  rl=mi+bias
  rr=ma+bias
  cc=[(-rl)/(rr-rl),1/(rr-rl)]
  xma = max((v2)[0,*])
  xmi = min((v2)[0,*])
  xmid = 0.5*(xma+xmi)*cc[1]
  yma = max((v2)[1,*])
  ymi = min((v2)[1,*])
  ymid = 0.5*(yma+ymi)*cc[1]
  zma = max((v2)[2,*])
  zmi = min((v2)[2,*])
  zmid = 0.5*(zma+zmi)*cc[1]
  vol_HU_cube_size=SIZE(vol_HU_cube)
;  print,'xma',xma
;  print,'yma',yma
;  print,'zma',zma
;  help,vol_HU_cube
  ;;;;;;;;;;;;;;hxn-end;;;;;;;;;;;;

  ;    ;update sel_one
  ;    print, 'N_ELEMENTS(xd_sState.sel_one.verts) = ', N_ELEMENTS(*(xd_sState.sel_one.p_verts))
  ;    print, 'N_ELEMENTS(xd_sState.sel_one.conn) = ', N_ELEMENTS(*(xd_sState.sel_one.p_conn))
  sel_one.p_verts = PTR_NEW(sel_one.verts)
  sel_one.p_conn = PTR_NEW(sel_one.conn)


  ;  Get the current color vectors to restore
  ;  when this application is exited.
  TVLCT, savedR, savedG, savedB, /GET

  ; Build color table from color vectors
  ;
  colorTable = [[savedR],[savedG],[savedB]]

  ;  Load an initial file .
  ;
  file = 'abnorm.dat'
  demo_getdata, NewImage, FILENAME=file, /TWO_DIM
  z=NewImage[*,*,7]
  z=smooth(z, 3, /EDGE_TRUNCATE)
  siz = SIZE(z)
  z=REBIN(z, siz[1]/2, siz[2]/2)

  req_ysize = 1000L
  ;  Create widgets.
  ;
  if (N_ELEMENTS(group) EQ 0) then begin
    IF req_ysize GT scr[1] THEN BEGIN
      wTopBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, $
        /TLB_KILL_REQUEST_EVENTS, $
        MBAR=barBase, $
        bitmap=icon_dir+'BASE1.bmp', $
        UNAME='d_surfview1:tlb', $
        TITLE="MPI+: Modifying...", ysize = req_ysize, y_scroll_size = scr[1])
    ENDIF ELSE BEGIn
      wTopBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, $
        /TLB_KILL_REQUEST_EVENTS, $
        bitmap=icon_dir+'BASE1.bmp', $
        MBAR=barBase, $
        UNAME='d_surfview1:tlb', $
        TITLE="MPI+: Modifying...", ysize = req_ysize)
    ENDELSE
    
  endif else begin
    IF req_ysize LT scr[1] THEN BEGIN
      wTopBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, $
        /TLB_KILL_REQUEST_EVENTS, $
        bitmap=icon_dir+'BASE1.bmp', $
        MBAR=barBase, $
        GROUP_LEADER=group, $
        UNAME='d_surfview1:tlb', $
        TITLE="MPI+: Modifying...", ysize = req_ysize, y_scroll_size = scr[1])
    ENDIF ELSE BEGIN
      wTopBase = WIDGET_BASE(/COLUMN, XPAD=0, YPAD=0, $
        /TLB_KILL_REQUEST_EVENTS, $
        bitmap=icon_dir+'BASE1.bmp', $
        MBAR=barBase, $
        GROUP_LEADER=group, $
        UNAME='d_surfview1:tlb', $
        TITLE="MPI+: Modifying...", ysize = req_ysize)
    ENDELSE
    
  endelse

  ;  Create the menu bar. It contains the file,
  ;  edit, and help buttons.
  ;
  wFileButton = WIDGET_BUTTON(barBase, VALUE = 'File', /MENU)
;  wRdstl=WIDGET_BUTTON(wFileButton, value='Load stl file', $
;    UVALUE='rdstl', $
;    UNAME='d_surfview1:rdstl')
  wWrstl=WIDGET_BUTTON(wFileButton, value='save stl file', $
    UVALUE='wrstl', $
    UNAME='d_surfview1:wrstl')
  wsave = WIDGET_BUTTON(wFileButton, value='Save ZGG meas', $
    UVALUE='Save', $
    UNAME='d_surfview1:Save')
  wSaveAll = widget_button( wFileButton, value="Save all", $
    uval='Saveall', $
    UNAME='d_surfview1:Saveall')
  wRestoreAll = widget_button( wFileButton, value="Restore all", $
    uval='Restoreall', $
    UNAME='d_surfview1:Restoreall')
  wRotate360 = widget_button( wFileButton, value="Rotate 360", $
    uval='Rotate360', $
    UNAME='d_surfview1:Rotate360')
  wQuitButton = WIDGET_BUTTON(wFileButton, VALUE='Quit', $
    UNAME='d_surfview1:quit', $
    UVALUE='QUIT')

  ;  Create the menu bar item Edit
  ;  that has the shade and style options.
  ;
  wOptionButton = WIDGET_BUTTON(barBase, VALUE = 'Surface Options', /MENU)

  ;  Select the plot shading.
  ;
  wShadingButton = WIDGET_BUTTON(wOptionButton, $
    VALUE='Shading', UVALUE='SHADING', MENU=1)

  wFlatButton = WIDGET_BUTTON(wShadingButton, $
    VALUE='Flat', UVALUE='FLAT', UNAME='d_surfview1:flat')

  wGouraudButton = WIDGET_BUTTON(wShadingButton, $
    VALUE='Gouraud', UVALUE='GOURAUD', $
    UNAME='d_surfview:gouraud')

  ;  Select the plot style.
  ;
  wStyleButton = WIDGET_BUTTON(wOptionButton, $
    VALUE='Style', UVALUE='STYLE', /MENU)

  wPointButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Point', UVALUE='POINT', UNAME='d_surfview1:point')

  wWireButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Wire', UVALUE='WIRE', UNAME='d_surfview1:wire')

  wSolidButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Solid', UVALUE='SOLID', UNAME='d_surfview1:solid')

  wRuledXZButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Ruled XZ', UVALUE='RULEDXZ', $
    UNAME='d_surfview1:ruledxz')

  wRuledYZButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Ruled YZ', UVALUE='RULEDYZ', $
    UNAME='d_surfview1:ruledyz')

  wLegoWireButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Lego Wire', UVALUE='LEGOWIRE', $
    UNAME='d_surfview1:legowire')

  wLegoSolidButton = WIDGET_BUTTON(wStyleButton, $
    VALUE='Lego Solid', UVALUE='LEGOSOLID', $
    UNAME='d_surfview1:legosolid')

  ;  Select the skirt value.
  ;
  wSkirtButton = WIDGET_BUTTON(wOptionButton, $
    VALUE='Skirt', UVALUE='SKIRT', /MENU)

  ; Remove the skirt
  ;
  wSkirtNoneButton = WIDGET_BUTTON(wSkirtButton, $
    VALUE='None', UVALUE='SKIRTNONE', $
    UNAME='d_surfview1:skirtnone')

  ; Skirt10 is  -0.5.
  ;
  wSkirt10Button = WIDGET_BUTTON(wSkirtButton, $
    VALUE='z=6', UVALUE='SKIRT10', $
    UNAME='d_surfview1:skirt6')

  ; Skirt20 is  0.0
  ;
  wSkirt20Button = WIDGET_BUTTON(wSkirtButton, $
    VALUE='z=102', UVALUE='SKIRT20', $
    UNAME='d_surfview1:skirt102')

  ; Skirt30 is  0.5
  ;
  wSkirt30Button = WIDGET_BUTTON(wSkirtButton, $
    VALUE='z=198', UVALUE='SKIRT30', $
    UNAME='d_surfview1:skirt198')

  ;  Set up the drag quality. Low is wire, medium is
  ;  polygons, high is smoothed polygons.
  ;
  wDragButton = Widget_Button(wOptionButton, $
    VALUE="Drag Quality", UVALUE='DRAGQ', /MENU)

  wLowButton = WIDGET_BUTTON(wDragButton, $
    VALUE='Low', UVALUE='LOW', UNAME='d_surfview1:lowdrag')

  wMediumButton = WIDGET_BUTTON(wDragButton, $
    VALUE='Medium', UVALUE='MEDIUM', $
    UNAME='d_surfview1:meddrag')

  wHighButton = WIDGET_BUTTON(wDragButton, $
    VALUE='High', UVALUE='HIGH', $
    UNAME='d_surfview1:highdrag')

  ;  Allows to trace a contour on the surface, and
  ;  then to display only the surface contained within
  ;  that contour.
  ;
  wTracingButton = WIDGET_BUTTON(wOptionButton, $
    VALUE='Tracing', UVALUE='TRACING', /MENU)

  wTracingModeButton = WIDGET_BUTTON(wTracingButton, $
    VALUE='Activate Outlining', $
    UVALUE='TRACING_MODE', $
    UNAME='d_surfview1:tracing_mode')

  wTracingMaskButton = WIDGET_BUTTON(wTracingButton, $
    VALUE='Activate Trace Mask', UVALUE='TRACING_MASK', $
    UNAME='d_surfview1:tracing_mask')

  ;  Toggle between showing or not showing
  ;  the hidden points and lines.
  ;
  wHiddenButton = WIDGET_BUTTON(wOptionButton, $
    VALUE="Activate Hiding", UVALUE='HIDE', $
    UNAME='d_surfview1:hidden')
    
;  ;;;;;;;;;;;;;;;;;;;;;hjx1;;;;;;;;;;;;;;;;;;;;;;;;;;
;  wCOLButton = WIDGET_BUTTON(wOptionButton, $
;    VALUE="changecolors", UVALUE='changecolors',/menu)
;
;  wpinkButton = WIDGET_BUTTON(wCOLButton, $
;    VALUE="pink", UVALUE='pink', $
;    UNAME='d_surfview1:pink')
;
;  wwhiteButton = WIDGET_BUTTON(wCOLButton, $
;    VALUE="white", UVALUE='white', $
;    UNAME='d_surfview1:white')
;
;  wsalmonButton = WIDGET_BUTTON(wCOLButton, $
;    VALUE="salmon", UVALUE='salmon', $
;    UNAME='d_surfview1:salmon')
;
;  wblueButton = WIDGET_BUTTON(wCOLButton, $
;    VALUE="blue", UVALUE='blue', $
;    UNAME='d_surfview1:blue')
;
;
;  ;;;;;;;;;;;;;;;;;;;;;hjx1;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;  Toggle between showing or not showing
  ;  the vertices in colors.
  ;
  wVertexButton = WIDGET_BUTTON(wOptionButton, $
    VALUE="Activate Vertex Colors", UVALUE='VERTEXCOLOR', $
    UNAME='d_surfview1:vertcolor')

  ;  Toggle between showing or not showing
  ;  the texture mapping.
  ;
  wTextureButton = WIDGET_BUTTON(wOptionButton, $
    VALUE="Deactivate Texture Map", UVALUE='TEXTURE', $
    UNAME='d_surfview1:texture')

  ;  Toggle between a solid or dash line style.
  ;
  wLineStyleButton = WIDGET_BUTTON(wOptionButton, $
    VALUE="Activate Line Style", UVALUE='LINESTYLE', $
    UNAME='d_surfview1:linestyle')

  ;  Create the menu bar item Edit
  ;  that has the shade and style options.
  ;
;  wViewButton = WIDGET_BUTTON(barBase, VALUE = 'View', /MENU)
;
;  wAnimateButton = WIDGET_BUTTON(wViewButton, $
;    VALUE="Animate", UVALUE='ANIMATE', $
;    UNAME='d_surfview1:animate')
;
;  wResetButton = WIDGET_BUTTON(wViewButton, $
;    VALUE="Reset Orientation", UVALUE='RESETTRANSFORM', $
;    UNAME='d_surfview1:reset_transfrm')


  ;  Create a sub base of the top base (wBase).
  ;
  wSubBase = WIDGET_BASE(wTopBase, COLUMN=2)

  ;  Create the left Base that contains the functionality buttons.
  ;  Here the only button is to animate the object.
  ;
  wLeftbase = WIDGET_BASE(wSubBase, $ ;/BASE_ALIGN_CENTER, $
    COLUMN=1)

  ;define hole size, GJ
  minHoleValue = 1
  maxHoleValue = 9
  initHoleValue = 5
;  wHoleLabel = WIDGET_LABEL(wLeftBase, $
;    VALUE='Hole Size:')
;
;  wHoleSlider = WIDGET_SLIDER(wLeftBase, $
;    MINIMUM=minHoleValue, $
;    MAXIMUM=maxHoleValue, VALUE=initHoleValue, $
;    UVALUE='HOLESIZE', $
;    UNAME='d_surfview1:holesize_slider')
   ;end changing hole size

  ;define threshold value, GJ
  initThresholdValue_min = sel_one.HU_min_value
  initThresholdValue_max = sel_one.HU_max_value
  minValue = MIN(vol_HU_cube)+3;MIN([MIN(vol_HU_cube)-1, initThresholdValue_min-500]); -299
  maxValue = MAX([MAX(vol_HU_cube)-5, initThresholdValue_max]);MAX([MAX(vol_HU_cube)+1, initThresholdValue_max+500]) ; 3000
  wThreshold_min_base = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
;  wThresholdLabel = WIDGET_LABEL(wThreshold_min_base, $
;    VALUE='Threshold_min:')
    
  wThresholdminLabel = WIDGET_BUTTON(wThreshold_min_base, /ALIGN_CENTER, $
    value=icon_dir+'thresholdmin.bmp', /bitmap, $
    UVALUE='thresholdminlabel', $
    UNAME='d_surfview1:thresholdminlabel')
  widget_control, wThresholdminLabel, sensitive=0
  
  IF minValue GT initThresholdValue_min THEN minValue = initThresholdValue_min-0.1
  wThresholdSlider_min = WIDGET_SLIDER(wThreshold_min_base, $
    MINIMUM=minValue, $
    MAXIMUM=maxValue, VALUE=initThresholdValue_min, $
    SCR_XSIZE=256, UVALUE='THRESHOLD_min', $
    UNAME='d_surfview1:threshold_min_slider')
  
  wThresholdText_min = WIDGET_TEXT(wThreshold_min_base, VALUE=string(format='(I5)', initThresholdValue_min), $
    UVALUE='THRESHOLD_min_text', YOFFSET = 5, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:threshold_min_text')

  wThreshold_max_base = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
;  wThresholdLabel = WIDGET_LABEL(wThreshold_max_base, $
;    VALUE='Threshold_max:')
 
  wThresholdmaxLabel = WIDGET_BUTTON(wThreshold_max_base, /ALIGN_CENTER, $
    value=icon_dir+'thresholdmax.bmp', /bitmap, $
    UVALUE='thresholdmaxlabel', $
    UNAME='d_surfview1:thresholdmaxlabel')
  widget_control, wThresholdmaxLabel, sensitive=0

  wThresholdSlider_max = WIDGET_SLIDER(wThreshold_max_base, $
    MINIMUM=minValue, $
    MAXIMUM=maxValue, VALUE=initThresholdValue_max, $
    SCR_XSIZE=256, UVALUE='THRESHOLD_max', $
    UNAME='d_surfview1:threshold_max_slider')
    
  wThresholdText_max = WIDGET_TEXT(wThreshold_max_base, VALUE=string(format='(I5)', initThresholdValue_max), $
    UVALUE='THRESHOLD_max_text', YOFFSET = 200, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:threshold_max_text')
  
  

  wScalingBase = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)

  percent = 100
  scalingString = 'Scaling : ' + STRING(percent, $
    FORMAT='(f5.1)') + ' %'
;  wScalingLabel = WIDGET_LABEL(wScalingBase, $
;    VALUE=scalingString)
  wScalingLabel = WIDGET_BUTTON(wScalingBase, /ALIGN_CENTER, $
    value=icon_dir+'scale.bmp', /bitmap, $
    UVALUE='scalinglabel', $
    UNAME='d_surfview1:scalinglabel')
  widget_control, wScalingLabel, sensitive=0    

  wScalingSlider = WIDGET_SLIDER(wScalingBase, $
    MINIMUM=0, $
    MAXIMUM=425, VALUE=25, SCR_XSIZE=300, $
;    /SUPPRESS_VALUE, $
    UVALUE='SCALING', $
    UNAME='d_surfview1:zoom')

 
  wObjectIndex_base = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  wObjectIndexLabel = WIDGET_BUTTON(wObjectIndex_base, /ALIGN_CENTER, $
    value=icon_dir+'objectindex.bmp', /bitmap, $
    UVALUE='ObjectIndexlabel', $
    UNAME='d_surfview1:ObjectIndexlabel')
  widget_control, wObjectIndexLabel, sensitive=0
  wObjectIndex = WIDGET_SLIDER(wObjectIndex_base, $
    MINIMUM=0, $
    MAXIMUM=n_Polys, VALUE=1, $ ;GJ 2020/5/25 change value to 1
    SCR_XSIZE=256, UVALUE='ObjectIndex', $
    UNAME='d_surfview1:ObjectIndex')
  wObjectIndex_text = WIDGET_TEXT(wObjectIndex_base, VALUE=string(format='(I3)', 1), $ ;GJ 2020/5/25 change value to 1
    UVALUE='ObjectIndex_text', YOFFSET = 200, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:ObjectIndex_text')
    
  wColorIndex_base = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
;  wColorIndexLabel1 = WIDGET_BUTTON(wColorIndex_base, /ALIGN_CENTER, $
;    value=icon_dir+'colorindex.bmp', /bitmap, $
;    UVALUE='ColorIndexlabel1', $
;    UNAME='d_surfview1:ColorIndexlabel1')
;  widget_control, wColorIndexLabel1, sensitive=0
  wColorIndexSliderBase = WIDGET_BASE(wColorIndex_base, /COL)
  wColorIndexLabel2 = WIDGET_button(wColorIndexSliderBase, /ALIGN_CENTER, SCR_XSIZE=286, $
    VALUE=icon_dir+'loadct.bmp', /BITMAP, $
    UVALUE='ColorIndexlabel2', $
    UNAME='d_surfview1:ColorIndexlabel2')

  wColorIndex = WIDGET_SLIDER(wColorIndexSliderBase, /ALIGN_CENTER, $
    MINIMUM=0, $
    MAXIMUM=30, VALUE=1, /SUPPRESS_VALUE, $
    SCR_XSIZE=316, UVALUE='ColorIndex', $
    UNAME='d_surfview1:ColorIndex')
  widget_control, wColorIndex, sensitive=1
;   
  wColorIndex_draw = WIDGET_DRAW(wColorIndex_base, /BUTTON_EVENTS, $
    UVALUE='singlehide', YOFFSET = 200, SCR_XSIZE=40, SCR_YSIZE=30, $
    UNAME='d_surfview1:singlehide')
  
  wSmooth_base = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  percent = 12.5
  SmoothString = 'Smooth : ' + STRING(percent, $
    FORMAT='(f5.1)') + ' %'

  wSmoothLabel = WIDGET_BUTTON(wSmooth_base, /ALIGN_CENTER, $
    value=icon_dir+'smooth.bmp', /bitmap, $
    UVALUE='Smoothlabel', $
    UNAME='d_surfview1:Smoothlabel')
  widget_control, wSmoothLabel, sensitive=0
  
  wSmooth= WIDGET_SLIDER(wSmooth_base, $
    MINIMUM=0, $
    MAXIMUM=400, VALUE=0, $
    SCR_XSIZE=266, $
    UVALUE='Smooth', $
    UNAME='d_surfview1:smooth_slider')
  widget_control, wSmooth, sensitive=1
    
  wSmooth_text = WIDGET_TEXT(wSmooth_base, VALUE=string(format='(f5.1)', 0), $
    UVALUE='Smooth_text', YOFFSET = 200, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:Smooth_text')
  widget_control, wSmooth_text, sensitive=1
  
  WTRANSSLIDERVALUE=100
  wTrans_base = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  percent = 100
  transString = 'transparency : ' + STRING(percent, $
    FORMAT='(f5.1)') + ' %'
;  wtransLabel = WIDGET_LABEL(wTrans_base, $
;    VALUE=transString)
  wTransLabel = WIDGET_BUTTON(wTrans_base, /ALIGN_CENTER, $
    value=icon_dir+'transparent.bmp', /bitmap, $
    UVALUE='Translabel', $
    UNAME='d_surfview1:Translabel')
  widget_control, wtransLabel, sensitive=0

  wtransSlider = WIDGET_SLIDER(wTrans_base, $
    MINIMUM=0, $
    MAXIMUM=100, VALUE=100, $
    SCR_XSIZE=266, $
    UVALUE='wtransSlider', $
    UNAME='d_surfview1:wtransSlider')
  widget_control, wtransSlider, sensitive=1
  
  wTrans_text = WIDGET_TEXT(wTrans_base, VALUE=string(format='(f5.1)', 100), $
    UVALUE='Trans_text', YOFFSET = 200, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:Trans_text')
  widget_control, wTrans_text, sensitive=1
    
  ;seg model
  wseg = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  w3DSegd = WIDGET_BUTTON(wseg,  /ALIGN_CENTER, $
    value=icon_dir+'3DSEGD.bmp', /bitmap, $
    UVALUE='3DSEGD', $
    UNAME='d_surfview1:3dsegd')
  widget_control, w3DSegd, sensitive=1
  w3Drg = WIDGET_BUTTON(wseg,  /ALIGN_CENTER, $
    value=icon_dir+'3DRG.bmp', /bitmap, $
    UVALUE='3DRG', $
    UNAME='d_surfview1:3drg')
 widget_control, w3Drg, sensitive=1
  wthickness = WIDGET_BUTTON(wseg,  /ALIGN_CENTER, $
    value=icon_dir+'thickness.bmp', /bitmap, $
    UVALUE='thickness', $
    UNAME='d_surfview1:thickness')
  widget_control, wthickness, sensitive=1
  WIDGET_CONTROL,wseg,/realize
  
  ;style model
  wstylebase = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  wpoint = WIDGET_BUTTON(wstylebase,  /ALIGN_CENTER, $
    value=icon_dir+'pointstyle.bmp', /bitmap, $
    UVALUE='pointstyle', $
    UNAME='d_surfview1:pointstyle')
  widget_control, wpoint, sensitive=1
  wwire = WIDGET_BUTTON(wstylebase,  /ALIGN_CENTER, $
    value=icon_dir+'wirestyle.bmp', /bitmap, $
    UVALUE='wirestyle', $
    UNAME='d_surfview1:wirestyle')
  widget_control, wwire, sensitive=1
  wfill = WIDGET_BUTTON(wstylebase,  /ALIGN_CENTER, $
    value=icon_dir+'fillstyle.bmp', /bitmap, $
    UVALUE='fillstyle', $
    UNAME='d_surfview1:fillstyle')
  widget_control, wfill, sensitive=1
  WIDGET_CONTROL,wstylebase,/realize
  
;  ;stl model
;  wstl = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
;;  wrdstl = WIDGET_BUTTON(wstl,  /ALIGN_CENTER, $
;;    value=icon_dir+'readstl.bmp', /bitmap, $
;;    UVALUE='rdstl', $
;;    UNAME='d_surfview1:rdstl')
;  wexpstl = WIDGET_BUTTON(wstl,  /ALIGN_CENTER, $
;    value=icon_dir+'expstl.bmp', /bitmap, $
;    UVALUE='expstl', $
;    UNAME='d_surfview1:expstl')
;  widget_control, wexpstl, sensitive=0
;  WIDGET_CONTROL,wstl,/realize
  
  ;hiding model
  whidingbase = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
;  whide = WIDGET_BUTTON(whidingbase,  /ALIGN_CENTER, $
;    value=icon_dir+'singlehide.bmp', /bitmap, $
;    UVALUE='singlehide', $
;    UNAME='d_surfview1:singlehide')
;  widget_control, whide, sensitive=0
;  whidecancel = WIDGET_BUTTON(whidingbase,  /ALIGN_CENTER, $
;    value=icon_dir+'singlehidecancel.bmp', /bitmap, $
;    UVALUE='singlehidecancel', $
;    UNAME='d_surfview1:singlehidecancel')
;  widget_control, whidecancel, sensitive=0
    
  wsinglecenter = WIDGET_BUTTON(whidingbase,  /ALIGN_CENTER, $
    value=icon_dir+'singlecenter.bmp', /bitmap, $
    UVALUE='singlecenter', $
    UNAME='d_surfview1:singlecenter')
  widget_control, wsinglecenter, sensitive=1
  wwholecenter = WIDGET_BUTTON(whidingbase,  /ALIGN_CENTER, $
    value=icon_dir+'wholecenter.bmp', /bitmap, $
    UVALUE='wholecenter', $
    UNAME='d_surfview1:wholecenter')
    
  wHideAxes = WIDGET_BUTTON(whidingbase, /ALIGN_CENTER, $
    value=icon_dir+'showaxes.bmp', /bitmap, $
    UVALUE='HIDEAXES', $
    UNAME='d_surfview1:hide_axes')
  WIDGET_CONTROL,whidingbase,/realize
  
  ;undo redo model
  wdobase = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  wUndo = Widget_Button(wdobase,  $
    UNAME='d_surfview1:undo' ,  /ALIGN_CENTER, $
    value=icon_dir+'undo.bmp', /bitmap, $
    UVALUE='UNDO')
  wRedo = Widget_Button(wdobase,  $
    UNAME='d_surfview1:redo' ,  /ALIGN_CENTER, $
    value=icon_dir+'redo.bmp', /bitmap, $
    UVALUE='REDO')
  wexpstl = WIDGET_BUTTON(wdobase,  /ALIGN_CENTER, $
    value=icon_dir+'expstl.bmp', /bitmap, $
    UVALUE='expstl', $
    UNAME='d_surfview1:expstl')
;  widget_control, wexpstl, sensitive=0
  WIDGET_CONTROL,wdobase,/realize
  
  ;seg model
  wseg_face = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  w3Drg_face = WIDGET_BUTTON(wseg_face,  /ALIGN_CENTER, $
    value='FacePicReg', $
    UVALUE='FacePicReg', $
    UNAME='d_surfview1:FacePicReg')
  widget_control, w3Drg_face, sensitive=1
  
  wAIVesselSeg = WIDGET_BUTTON(wseg_face,  /ALIGN_CENTER, $
    value='AIVesselSeg', $
    UVALUE='AIVesselSeg', $
    UNAME='d_surfview1:AIVesselSeg')
  widget_control, wAIVesselSeg, sensitive=1

  wVesselEcc = WIDGET_BUTTON(wseg_face,  /ALIGN_CENTER, $
    value='VesselRec', $
    UVALUE='eccentricity', $
    UNAME='d_surfview1:eccentricity')
  widget_control, wVesselEcc, sensitive=1
  
  wVWCESeg = WIDGET_BUTTON(wseg_face,  /ALIGN_CENTER, $
    value='VWCESeg', $
    UVALUE='VWCESeg', $
    UNAME='d_surfview1:VWCESeg')
  widget_control, wVWCESeg, sensitive=1

  WIDGET_CONTROL,wseg_face,/realize
  
  ;define threshold value, GJ
  vol_HU_cube_a1 = vol_HU_cube
  initThresholdValue_min_a1 = initThresholdValue_min
  initThresholdValue_max_a1 = initThresholdValue_max

  minValue_a1 = MIN(vol_HU_cube)+3;MIN([MIN(vol_HU_cube)-1, initThresholdValue_min-500]); -299
  maxValue_a1 = MAX([MAX(vol_HU_cube)-5, initThresholdValue_max_a1]);MAX([MAX(vol_HU_cube)+1, initThresholdValue_max+500]) ; 3000
  wThreshold_min_base_a1 = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)

  wThresholdminLabel_a1 = WIDGET_BUTTON(wThreshold_min_base_a1, /ALIGN_CENTER, $
    value=icon_dir+'thresholdmin.bmp', /bitmap, $
    UVALUE='thresholdminlabel_a1', $
    UNAME='d_surfview1:thresholdminlabel_a1')
  widget_control, wThresholdminLabel_a1, sensitive=0

  IF minValue_a1 GT initThresholdValue_min_a1 THEN minValue_a1 = initThresholdValue_min_a1-0.1
  wThresholdSlider_min_a1 = WIDGET_SLIDER(wThreshold_min_base_a1, $
    MINIMUM=minValue_a1, $
    MAXIMUM=maxValue_a1, VALUE=initThresholdValue_min_a1, $
    SCR_XSIZE=256, UVALUE='THRESHOLD_min_a1', $
    UNAME='d_surfview1:threshold_min_slider_a1')

  wThresholdText_min_a1 = WIDGET_TEXT(wThreshold_min_base_a1, VALUE=string(format='(I5)', initThresholdValue_min_a1), $
    UVALUE='THRESHOLD_min_text_a1', YOFFSET = 5, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:threshold_min_text_a1')

  wThreshold_max_base_a1 = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)

  wThresholdmaxLabel_a1 = WIDGET_BUTTON(wThreshold_max_base_a1, /ALIGN_CENTER, $
    value=icon_dir+'thresholdmax.bmp', /bitmap, $
    UVALUE='thresholdmaxlabel_a1', $
    UNAME='d_surfview1:thresholdmaxlabel_a1')
  widget_control, wThresholdmaxLabel_a1, sensitive=0

  wThresholdSlider_max_a1 = WIDGET_SLIDER(wThreshold_max_base_a1, $
    MINIMUM=minValue_a1, $
    MAXIMUM=maxValue_a1, VALUE=initThresholdValue_max_a1, $
    SCR_XSIZE=256, UVALUE='THRESHOLD_max_a1', $
    UNAME='d_surfview1:threshold_max_slider_a1')

  wThresholdText_max_a1 = WIDGET_TEXT(wThreshold_max_base_a1, VALUE=string(format='(I5)', initThresholdValue_max_a1), $
    UVALUE='THRESHOLD_max_text_a1', YOFFSET = 200, SCR_XSIZE=40, SCR_YSIZE=30, $
    /EDITABLE, UNAME='d_surfview1:threshold_max_text_a1')

  IF N_ELEMENTS(additionalDir) GE 2 THEN BEGIN
    widget_control, wThreshold_min_base_a1, map=1
    widget_control, wThreshold_max_base_a1, map=1
  ENDIF ELSE BEGIN
    widget_control, wThreshold_min_base_a1, map=0
    widget_control, wThreshold_max_base_a1, map=0
  ENDELSE
  
  ;define a fusion and flip and move base
  wFusionFlipControlBase = WIDGET_BASE(wLeftBase, FRAME=1, /ROW)
  
  ;Fusion mode
  wFusionControlBase = WIDGET_BASE(wFusionFlipControlBase, FRAME=1, /ROW)
;  wFusionLabel = WIDGET_LABEL(wFusionControlBase, VALUE='MPI-CT Fusion')
  wFusionModeRadio = cw_bgroup( $
    wFusionControlBase, $
    ['No', 'Fusion'], $
    /exclusive, $
    /no_release, $
    set_value=0, $
    uvalue='FUSIONMODE' $
    )
  widget_control, wFusionModeRadio, set_uname='d_surfview1:fusionmoderadio', sensitive=1

  ;x flip mode
  wxFlipControlBase = WIDGET_BASE(wFusionFlipControlBase, FRAME=1, /ROW)
  wxFlipModeRadio = cw_bgroup( $
    wxFlipControlBase, $
    ['No Flip', 'x Flip'], $
    /exclusive, $
    /no_release, $
    set_value=0, $
    uvalue='XFLIPMODE' $
    )
  widget_control, wxFlipModeRadio, set_uname='d_surfview1:xflipmoderadio', sensitive=1

  ;y flip mode
  wyFlipControlBase = WIDGET_BASE(wFusionFlipControlBase, FRAME=1, /ROW)
  wyFlipModeRadio = cw_bgroup( $
    wyFlipControlBase, $
    ['No Flip', 'y Flip'], $
    /exclusive, $
    /no_release, $
    set_value=0, $
    uvalue='YFLIPMODE' $
    )
  widget_control, wyFlipModeRadio, set_uname='d_surfview1:yflipmoderadio', sensitive=1

  ;z flip mode
  wzFlipControlBase = WIDGET_BASE(wFusionFlipControlBase, FRAME=1, /ROW)
  wzFlipModeRadio = cw_bgroup( $
    wzFlipControlBase, $
    ['No Flip', 'z Flip'], $
    /exclusive, $
    /no_release, $
    set_value=0, $
    uvalue='ZFLIPMODE' $
    )
  widget_control, wzFlipModeRadio, set_uname='d_surfview1:zflipmoderadio', sensitive=1

  ;@GJ, show or hide the fusion or flip mode
  IF N_ELEMENTS(additionalDir) GE 2 THEN widget_control, wFusionFlipControlBase, map=1 ELSE widget_control, wFusionFlipControlBase, map=0

  
  ;;;;;;;;hxn1-start;;;;;;;;;;;
  hSlice_base1 = Widget_Base(wLeftBase,  $
    UNAME='hSlice_base1' ,FRAME=1   $
    ,SCR_XSIZE=169, SCR_YSIZE=140, SPACE=3, XPAD=3, YPAD=3)

  x_move_label = Widget_Label(hSlice_base1,  $
    UNAME='x_move_label' ,XOFFSET=115 ,YOFFSET=2  $
    ,SCR_XSIZE=18 ,SCR_YSIZE=25 ,/ALIGN_CENTER ,VALUE='XMOVE')

  y_move_label = Widget_Label(hSlice_base1,  $
    UNAME='y_move_label' ,XOFFSET=115 ,YOFFSET=29  $
    ,SCR_XSIZE=18 ,SCR_YSIZE=25 ,/ALIGN_CENTER ,VALUE='YMOVE')

  z_move_label = Widget_Label(hSlice_base1,  $
    UNAME='z_move_label' ,XOFFSET=115 ,YOFFSET=62  $
    ,SCR_XSIZE=19 ,SCR_YSIZE=25 ,/ALIGN_CENTER ,VALUE='ZMOVE')


  wx_move_slider = Widget_Slider(hSlice_base1,  $
    UNAME='x_move_slider' ,XOFFSET=150 ,YOFFSET=2  $
    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=-(vol_HU_cube_size[1]-1)/2 ,MAXIMUM=(vol_HU_cube_size[1]-1)/2 ,VALUE=0 ,$
    UVALUE='xmoveslider')

  wy_move_slider = Widget_Slider(hSlice_base1,  $
    UNAME='y_move_slider' ,XOFFSET=150 ,YOFFSET=29  $
    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=-(vol_HU_cube_size[2]-1)/2 ,MAXIMUM=(vol_HU_cube_size[2]-1)/2 ,VALUE=0 ,$
    UVALUE='ymoveslider')

  wz_move_slider = Widget_Slider(hSlice_base1,  $
    UNAME='z_move_slider' ,XOFFSET=150 ,YOFFSET=62  $
    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=-(vol_HU_cube_size[3]-1)/2 ,MAXIMUM=(vol_HU_cube_size[3]-1)/2 ,VALUE=0 ,$
    UVALUE='zmoveslider')

    wxyzflipmove_output = Widget_Button(hSlice_base1,  $
    UNAME='xyzflipmove_save' ,XOFFSET=40 ,YOFFSET=95 $
    ,SCR_XSIZE=62 ,SCR_YSIZE=20 ,/ALIGN_CENTER ,value=icon_dir+'analyse.bmp', /bitmap, $
    UVALUE='xyzflipmove_save')
    
    wxyzflipmove_apply = Widget_Button(hSlice_base1,  $
    UNAME='xyzflipmove_apply' ,XOFFSET=150 ,YOFFSET=95 $
    ,SCR_XSIZE=62 ,SCR_YSIZE=20 ,/ALIGN_CENTER ,value=icon_dir+'analyse.bmp', /bitmap, $
    UVALUE='xyzflipmove_apply')
    
    WIDGET_CONTROL, wx_move_slider, sensitive = 1
    WIDGET_CONTROL, wy_move_slider, sensitive = 1
    WIDGET_CONTROL, wz_move_slider, sensitive = 1
    WIDGET_CONTROL, wxyzflipmove_output, sensitive = 1
    WIDGET_CONTROL, wxyzflipmove_apply, sensitive = 1
    
    IF N_ELEMENTS(additionalDir) GE 2 THEN widget_control, hSlice_base1, map=1 ELSE widget_control, hSlice_base1, map=0
    
;
;    oPoly1_DeleteButton= Widget_Button(hSlice_base1,  $
;    UNAME='opoly1delete' ,XOFFSET=140 ,YOFFSET=2  $
;    ,SCR_XSIZE=62 ,SCR_YSIZE=20 ,/ALIGN_CENTER ,value=icon_dir+'delete.bmp', /bitmap, $
;    UVALUE='OPOLY1DEL')
;    
;;    wloadctbutton = WIDGET_button(hSlice_base1, XOFFSET=10, YOFFSET=135, SCR_XSIZE=300, $
;;    VALUE='./lib/icon_pics/loadct.bmp', /BITMAP,/DYNAMIC_RESIZE,UVALUE='colorbutton1', $
;;    UNAME='d_surfview1:colorbutton1')
;;    wloadctSlider = WIDGET_SLIDER(hSlice_base1, XOFFSET=10, YOFFSET=155, $
;;    MINIMUM=minValue_color, SCR_XSIZE=300, $
;;    MAXIMUM=maxValue_color, VALUE=wloadctSlidervalue, /SUPPRESS_VALUE, $
;;    UVALUE='wloadctSlider1', $
;;    UNAME='d_surfview1:wloadctSlider1')
;    
;    WIDGET_CONTROL, x_scale_slider, sensitive = 0
;    WIDGET_CONTROL, y_scale_slider, sensitive = 0
;    WIDGET_CONTROL, z_scale_slider, sensitive = 0
;    WIDGET_CONTROL, x_move_slider, sensitive = 0
;    WIDGET_CONTROL, y_move_slider, sensitive = 0
;    WIDGET_CONTROL, z_move_slider, sensitive = 0
;    WIDGET_CONTROL, xyzscalemove_output, sensitive = 0
;
;  ;;;;;gj20181219
;
; 
;  ;;;;;;;;hxn1-start;;;;;;;;;;;
;  hSlice_base = Widget_Base(wLeftBase,  $
;    UNAME='hSlice_base' ,FRAME=1   $
;    ,SCR_XSIZE=169 ,SCR_YSIZE=165 ,SPACE=3 ,XPAD=3 ,YPAD=3)
;
;  ;  Slice_label = Widget_Label(hSlice_base,  $
;  ;    UNAME='Slice_label' ,XOFFSET=12 ,YOFFSET=3  $
;  ;    ,SCR_XSIZE=150 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='Slice Operation')
;  ;
;  ;  hSlice_base_line = Widget_Base(hSlice_base,  $
;  ;    UNAME='hSlice_base_line' ,FRAME=1 ,XOFFSET=6  $
;  ;    ,YOFFSET=25 ,SCR_XSIZE=146 ,SCR_YSIZE=2 ,SPACE=3 ,XPAD=3 ,YPAD=3)
;
;  exbase = WIDGET_BASE(hSlice_base,XOFFSET=2 ,YOFFSET=5,/EXCLUSIVE,ROW=1)
;  hSlice_on = WIDGET_BUTTON(exbase, XOFFSET=2 ,YOFFSET=5,value ='Split',UVALUE='hfsliceon')
;
;  ;  hSlice_on = Widget_Button(Slice_on_base,  $
;  ;    UNAME='Slice_on' ,YOFFSET=3 ,SCR_XSIZE=66  $
;  ;    ,SCR_YSIZE=20 ,/ALIGN_LEFT ,VALUE='On/Off' ,$
;  ;    UVALUE='hfsliceon')
;  hSlice_three = WIDGET_BUTTON(exbase, XOFFSET=30 ,YOFFSET=5,value ='clip',UVALUE='hthree')
;
;  ;  hSlice_on_y_button = Widget_Button(exbase,  XOFFSET=130 ,YOFFSET=73,VALUE='Reverse' , UVALUE='hfsliceonybutton')
;
;  ;  hSlice_on_y_button = Widget_Button(Slice_on_y_base,  $
;  ;    UNAME='Slice_on_y_button' ,/ALIGN_LEFT ,VALUE='Reverse' ,$
;  ;    UVALUE='hfsliceonybutton')
;
;  ;  hSlice_three = Widget_Button(Slice_on_three_base,  $
;  ;    UNAME='Slice_three',XOFFSET=125 ,YOFFSET=165  $
;  ;    ,/ALIGN_LEFT ,VALUE='3_plane' ,$
;  ;    UVALUE='hthree')
;
;  hSlice_output = Widget_Button(hSlice_base,  $
;    UNAME='Slice_output' ,XOFFSET=40 ,YOFFSET=140  $
;    ,SCR_XSIZE=62 ,SCR_YSIZE=20 ,/ALIGN_CENTER ,value=icon_dir+'analyse.bmp', /bitmap, $
;    UVALUE='hfsliceoutput')
;
;
;  Slice_x_label = Widget_Label(hSlice_base,  $
;    UNAME='Slice_x_label' ,XOFFSET=2 ,YOFFSET=46  $
;    ,SCR_XSIZE=18 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='X')
;
;  Slice_y_label = Widget_Label(hSlice_base,  $
;    UNAME='Slice_y_label' ,XOFFSET=2 ,YOFFSET=80  $
;    ,SCR_XSIZE=18 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='Y')
;
;  Slice_z_label = Widget_Label(hSlice_base,  $
;    UNAME='Slice_z_label' ,XOFFSET=2 ,YOFFSET=113  $
;    ,SCR_XSIZE=19 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='Z')
;
;  hSlice_x_slider = Widget_Slider(hSlice_base,  $
;    UNAME='Slice_x_slider' ,XOFFSET=20 ,YOFFSET=29  $
;    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=1 ,MAXIMUM=vol_HU_cube_size[1]-1 ,VALUE=1 ,$
;    UVALUE='hfslicexslider')
;
;  hSlice_y_slider = Widget_Slider(hSlice_base,  $
;    UNAME='Slice_y_slider' ,XOFFSET=20 ,YOFFSET=62  $
;    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=1 ,MAXIMUM=vol_HU_cube_size[2]-1,VALUE=1 ,$
;    UVALUE='hfsliceyslider')
;
;  hSlice_z_slider = Widget_Slider(hSlice_base,  $
;    UNAME='Slice_z_slider' ,XOFFSET=20 ,YOFFSET=95  $
;    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=1 ,MAXIMUM=vol_HU_cube_size[3]-1 ,VALUE=1 ,$
;    UVALUE='hfslicezslider')
;
;
;  Slice_x_label1 = Widget_Label(hSlice_base,  $
;    UNAME='Slice_x_label1' ,XOFFSET=123 ,YOFFSET=46  $
;    ,SCR_XSIZE=18 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='-X')
;
;  Slice_y_label1 = Widget_Label(hSlice_base,  $
;    UNAME='Slice_y_label1' ,XOFFSET=123 ,YOFFSET=80  $
;    ,SCR_XSIZE=18 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='-Y')
;
;  Slice_z_label1 = Widget_Label(hSlice_base,  $
;    UNAME='Slice_z_label1' ,XOFFSET=123 ,YOFFSET=113  $
;    ,SCR_XSIZE=19 ,SCR_YSIZE=18 ,/ALIGN_CENTER ,VALUE='-Z')
;
;
;  hSlice_x_slider1 = Widget_Slider(hSlice_base,  $
;    UNAME='Slice_x_slider1' ,XOFFSET=150 ,YOFFSET=29  $
;    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=1 ,MAXIMUM=vol_HU_cube_size[1]-1 ,VALUE=vol_HU_cube_size[1]-1 ,$
;    UVALUE='hfslicexslider1')
;
;  hSlice_y_slider1 = Widget_Slider(hSlice_base,  $
;    UNAME='Slice_y_slider1' ,XOFFSET=150 ,YOFFSET=62  $
;    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=0 ,MAXIMUM=vol_HU_cube_size[2]-1 ,VALUE=vol_HU_cube_size[2]-1,$
;    UVALUE='hfsliceyslider1')
;
;  hSlice_z_slider1 = Widget_Slider(hSlice_base,  $
;    UNAME='Slice_z_slider1' ,XOFFSET=150 ,YOFFSET=95  $
;    ,SCR_XSIZE=100 ,SCR_YSIZE=32 ,MINIMUM=0 ,MAXIMUM=vol_HU_cube_size[3]-1 ,VALUE=vol_HU_cube_size[3]-1 ,$
;    UVALUE='hfslicezslider1')
;
;
;
;  Slice_on_label = Widget_Label(hSlice_base,  $
;    UNAME='Slice_on_label' ,XOFFSET=165,YOFFSET=13  $
;    ,/ALIGN_LEFT ,VALUE='Transparency')
;
;  ;  Slice_on_label = Widget_Label(hSlice_base,  $
;  ;    UNAME='Slice_on_label' ,XOFFSET=160 ,YOFFSET=110  $
;  ;    ,/ALIGN_LEFT ,VALUE='Reset')
;
;  ;  Slice_on_x_base = Widget_Base(hSlice_base,  $
;  ;    UNAME='Slice_on_x_base' ,XOFFSET=130 ,YOFFSET=103  $
;  ;    ,SCR_XSIZE=18 ,SCR_YSIZE=22 ,COLUMN=1 ,/NONEXCLUSIVE)
;  ;
;  ;  hSlice_on_x_button = Widget_Button(Slice_on_x_base,  $
;  ;    UNAME='Slice_on_x_button' ,/ALIGN_LEFT ,VALUE='' ,$
;  ;    UVALUE='hfsliceonxbutton')
;
;  Slice_on_y_base = Widget_Base(hSlice_base,  $
;    UNAME='Slice_on_y_base' ,XOFFSET=140 ,YOFFSET=5  $
;    ,SCR_XSIZE=18 ,SCR_YSIZE=22 ,ROW=1 ,/NONEXCLUSIVE)
;
;  hSlice_on_y_button = Widget_Button(Slice_on_y_base,  $
;    UNAME='Slice_on_y_button' ,/ALIGN_LEFT ,VALUE='Reverse' ,$
;    UVALUE='hfsliceonybutton')
;
;  WIDGET_CONTROL, hSlice_on_y_button, sensitive = 0
;;  WIDGET_CONTROL, hSlice_on_x_button, sensitive = 1
;  WIDGET_CONTROL, hSlice_output, sensitive = 0
;  WIDGET_CONTROL, hSlice_x_slider, sensitive = 0
;  WIDGET_CONTROL, hSlice_x_slider1, sensitive = 0
;  WIDGET_CONTROL, hSlice_y_slider, sensitive = 0
;  WIDGET_CONTROL, hSlice_y_slider1, sensitive = 0
;  WIDGET_CONTROL, hSlice_z_slider, sensitive = 0
;  WIDGET_CONTROL, hSlice_z_slider1, sensitive = 0
;
;
;wJZbase = WIDGET_BASE(wLeftBase, FRAME=1, $
;  /COLUMN)
;
;wJZlabel = WIDGET_LABEL(wJZbase, value='JZ analysis')
;;
;;wJZsegmentation = WIDGET_BUTTON(wJZbase, $
;;  value='JZ segmentation', $
;;  UVALUE='JZsegmentation', $
;;  UNAME='d_surfview1:jzsegmentation')
;;
;;
;
;  wZGG = WIDGET_BUTTON(wJZbase, $
;    value=icon_dir+'zgg.bmp', /bitmap, $
;    UVALUE='ZGG', $
;    UNAME='d_surfview1:zgg')
;
;;  Dir_imgckcomplete= FILE_WHICH('ckcomplete.bmp',/INCLUDE_CURRENT_DIR)
;;  Dir_imgckcompleteNot= FILE_WHICH('ckcompleteNot.bmp',/INCLUDE_CURRENT_DIR)
;
;  wCKbase = WIDGET_BASE(wJZbase, $
;    /ROW)
;  
;
;   wCKcomplete = WIDGET_BUTTON(wCKbase, $
;    value=icon_dir+'ckcomplete.bmp', /bitmap, $
;    UVALUE='CKcomplete', $
;    UNAME='d_surfview1:CKcomplete')
;    
;
;  wCKcompleteNot = WIDGET_BUTTON(wCKbase, $
;    value=icon_dir+'ckcompleteNot.bmp', /bitmap, $
;    UVALUE='CKcompleteNot', $
;    UNAME='d_surfview1:CKcompleteNot')
;  
;;  Dir_zb= FILE_WHICH('zbthickness.bmp',/INCLUDE_CURRENT_DIR)
;
;  wZB = WIDGET_BUTTON(wJZbase, $
;    value=icon_dir+'zbthickness.bmp', /bitmap, $
;    UVALUE='ZB', $
;    UNAME='d_surfview1:zb')
;  
;  WIDGET_CONTROL, wZGG, sensitive = 0
;  WIDGET_CONTROL, wCKcomplete, sensitive = 0
;  WIDGET_CONTROL, wCKcompleteNot, sensitive = 0
;  WIDGET_CONTROL, wZB, sensitive = 0
;
;  ;;;;;;;;;;HJX;;;;;;;;;;;;;;;;;;;;;
;  tlb = WIDGET_BASE(wLeftBase,title='measure ZGG', FRAME=1, $
;    uname='measure', /COLUMN)
;
;  wJZlabel = WIDGET_LABEL(tlb, value='ZGG measurement')
;
;  table1 = WIDGET_TABLE(tlb, $
;    uName  = 'table1',uv  = 'table1', $
;    /Editable,$
;    column_labels =['left ZGG','right ZGG','left ZGG','right ZGG'],$
;    row_labels =['3','4','5','6','7'],$
;    ;;;;;;;;;;HJX;;;;;;;;;;;;;;;;;;;;;
;    xsize = 4, ysize = 5, $
;    ;;;;;;;;;;HJX;;;;;;;;;;;;;;;;;;;;;
;    /All_Events,x_scroll_size = 2, y_scroll_size = 5, /scroll,/RESIZEABLE_COLUMNS)
;



;    ;;;;;;;;;;;;;;;;;;;;hjx1212
;
;    wdistance = WIDGET_BUTTON(tlb, $
;      value='vessel_distance', $
;      UVALUE='vessel_distance', $
;      UNAME='d_surfview1:vessel_distance')
;
;    ;;;;;;;;;;;;;;;;;;;;hjx1212


;  tlb2= WIDGET_BASE(tlb,title='Save', /ALIGN_CENTER,$
;    uname='Save')
;  wsave = WIDGET_BUTTON(tlb2,  /ALIGN_CENTER, $
;    value='Save ZGG meas', $
;    UVALUE='Save', $
;    UNAME='d_surfview1:Save')
;  WIDGET_CONTROL,tlb,/realize
;  WIDGET_CONTROL,tlb2,/realize



  max_dim = MAX((SIZE(vol_HU_cube))[1:3], max_ind)
  
  ;mask_vol_cube = vol_HU_cube * 0. -1000.
  xy_mask = BYTARR(max_dim, max_dim) * 0 + 1
  xy_mask_ROI = xy_mask
  xy_mask_clip = xy_mask
  yz_mask = BYTARR(max_dim, max_dim) * 0 + 1
  yz_mask_ROI = yz_mask
  yz_mask_clip = yz_mask
  xz_mask = BYTARR(max_dim, max_dim) * 0 + 1
  xz_mask_ROI = xz_mask
  xz_mask_clip = xz_mask
  
  ;


  ;  Initialize sizes.
  ;
  if (N_ELEMENTS(ImageSizeIn) LE 0) then begin
    DEVICE, GET_SCREEN = x
    if (x[0] GT 1024) then ImageSize=x[0]/5. $; 192*2. $
    else if x[0] ge 800 then ImageSize=128*2. $
    else ImageSize = 64*2.
  endif else begin
    imageSize = imageSizeIn
  endelse

  ;  Initialize working parameters and arrays.
  ;
  draw = LONARR(4)
  DrawSlider = LONARR(3)
  window = LONARR(6)
  ocolors = INTARR(4)             ; Overlay colors

    simageSize = FLOOR(imageSize*0.7)
  location3D = LONARR(3)*0 + FLOOR(simageSize)/2
  svol_HU_cube = CONGRID(vol_HU_cube, simageSize, simageSize, simageSize) 
  ;  Create the right Base that has the drawing area.
  ;
  wRightbase = WIDGET_BASE(wSubBase, /ROW)
  wRow1Base = WIDGET_BASE(wRightBase)
  wRow2Base = WIDGET_BASE(wRightBase, /COLUMN)
   
  temp = LONARR(3)
  for i=0,2 do begin
    temp[i] = WIDGET_BASE(wRow2Base, /COLUMN)
    Draw[i] = WIDGET_DRAW(temp[i], /BUTTON, $
      RETAIN=0, /WHEEL_EVENTS, $
      XSIZE=simageSize, YSIZE=simageSize)
    DrawSlider[i] = WIDGET_SLIDER(temp[i], $
      MINIMUM=1, /DRAG, $
      MAXIMUM=simageSize, VALUE=simageSize/2)
  endfor
  WIDGET_CONTROL, Draw[0], SET_UVALUE = 'DRAW0', SET_UNAME='d_surfview1:draw0'
  WIDGET_CONTROL, DrawSlider[0], SET_UVALUE = 'DRAWSLIDER0', SET_UNAME='d_surfview1:drawSlider0'
  WIDGET_CONTROL, Draw[1], SET_UVALUE = 'DRAW1', SET_UNAME='d_surfview1:draw1'
  WIDGET_CONTROL, DrawSlider[1], SET_UVALUE = 'DRAWSLIDER1', SET_UNAME='d_surfview1:drawSlider1'
  WIDGET_CONTROL, Draw[2], SET_UVALUE = 'DRAW2', SET_UNAME='d_surfview1:draw2'
  WIDGET_CONTROL, DrawSlider[2], SET_UVALUE = 'DRAWSLIDER2', SET_UNAME='d_surfview1:drawSlider2'
  
  temp_a1 = LONARR(3)
  draw_a1 = LONARR(3)
  DrawSlider_a1 = LONARR(3)
  temp_a2 = LONARR(3)
  draw_a2 = LONARR(3)
  DrawSlider_a2 = LONARR(3)
  vol_HU_cube_a1 = vol_HU_cube * 0.
  vol_HU_cube_a2 = vol_HU_cube * 0.
  svol_HU_cube_a1 = svol_HU_cube * 0.
  svol_HU_cube_a2 = svol_HU_cube * 0.
;  print, 'folders = ', additionalDir
;  additionalDir = ['C:\Temp\Research\Med_visualization\4330003_LIU_YING_20180907_HR_MRA\S0201_165711_s3DI_MC_HR\','C:\Temp\Research\Med_visualization\4330003_LIU_YING_20180907_HR_MRA\S0401_170837_3D_T1WIVesselWall\', 'C:\Temp\Research\Med_visualization\4330003_LIU_YING_20180907_HR_MRA\S1001_173516_3D_T1WIVesselWall\']
  vol_HU_cube_resolution1=vol_HU_cube_resolution
  IF KEYWORD_SET(additionalDir) THEN BEGIN
    IF N_ELEMENTS(additionalDir) GE 2 THEN BEGIN
      image_read, additionalDir[0], vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal = p_modal0
      image_read, additionalDir[1], vol_HU_cube_a1, vol_HU_cube_resolution1, patient_age1, direction1, origin_loc1, p_modal = p_modal1

      ;@GJ, 2023/5/16, do image alignment
      IF STRCMP(p_modal0, 'CT', 2, /FOLD_CASE) AND STRCMP(p_modal1, 'MPI', 3, /FOLD_CASE) THEN BEGIN
        vol_HU_cube_a1 = REVERSE(vol_HU_cube_a1, 2)
        x_shift = FIX((SIZE(vol_HU_cube, /DIMENSIONS))[0]/20.+(SIZE(vol_HU_cube, /DIMENSIONS))[1]/572.*26.)
        y_shift = FIX((SIZE(vol_HU_cube, /DIMENSIONS))[1]/10.)
        z_shift = FIX((SIZE(vol_HU_cube, /DIMENSIONS))[2]/100.*92)
        vol_HU_cube_a1 = SHIFT(vol_HU_cube_a1, -x_shift, -y_shift, z_shift)
      ENDIF
      ;@GJ, 2023/5/16, do image alignment
      IF STRCMP(p_modal0, 'MPI', 3, /FOLD_CASE) AND STRCMP(p_modal1, 'CT', 2, /FOLD_CASE) THEN BEGIN
        vol_HU_cube_a1 = REVERSE(vol_HU_cube_a1, 2)
        x_shift = FIX((SIZE(vol_HU_cube, /DIMENSIONS))[0]/20.+(SIZE(vol_HU_cube, /DIMENSIONS))[1]/572.*26.)
        y_shift = FIX((SIZE(vol_HU_cube, /DIMENSIONS))[1]/10.)
        z_shift = FIX((SIZE(vol_HU_cube, /DIMENSIONS))[2]/100.*92.)
        vol_HU_cube_a1 = SHIFT(vol_HU_cube_a1, x_shift, y_shift, -z_shift)
      ENDIF
      
      temp=volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, vol_HU_cube_a1, vol_HU_cube_resolution1, direction1, origin_loc1)
      ;@GJ, 2023/5/16, correcting the bug
      vol_HU_cube_a1 = temp
      IF MAX(temp)-MIN(temp) GT 0.001 THEN BEGIN
        ;this is MRI or data with correct fusion
        svol_HU_cube_a1 = CONGRID(temp, simageSize, simageSize, simageSize)
      ENDIF ELSE BEGIN
        ;this is PET/CT data, there is no need to do fusion
        ;GJ 2020/4/17
        svol_HU_cube_a1 = CONGRID(vol_HU_cube1, simageSize, simageSize, simageSize)
      ENDELSE
      
     ;calculate the additional image's threshold
      sizeVHc = SIZE(svol_HU_cube_a1)
      mean_image = IMAGE_THRESHOLD(BYTSCL(svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
      image_statistics, svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
      ;define threshold value, GJ
      initThresholdValue_min_a1 = MIN(svol_HU_cube_a1)
      initThresholdValue_max_a1 = MAX(svol_HU_cube_a1)
      minValue_a1 = MIN(vol_HU_cube_a1)+3;MIN([MIN(vol_HU_cube)-1, initThresholdValue_min-500]); -299
      maxValue_a1 = MAX([MAX(vol_HU_cube_a1)-5, initThresholdValue_max_a1]);MAX([MAX(vol_HU_cube)+1, initThresholdValue_max+500]) ; 3000
      IF minValue_a1 GT initThresholdValue_min_a1 THEN minValue_a1 = initThresholdValue_min_a1-0.1
      
      WIDGET_CONTROL, wThresholdSlider_min_a1, SET_VALUE=string(format='(I5)', initThresholdValue_min_a1), SET_SLIDER_MAX=maxValue_a1, SET_SLIDER_MIN=minValue_a1
      WIDGET_CONTROL, wThresholdText_min_a1, SET_VALUE=string(format='(I5)', initThresholdValue_min_a1)
      WIDGET_CONTROL, wThresholdSlider_max_a1, SET_VALUE=string(format='(I5)', initThresholdValue_max_a1), SET_SLIDER_MAX=maxValue_a1, SET_SLIDER_MIN=minValue_a1
      WIDGET_CONTROL, wThresholdText_max_a1, SET_VALUE=string(format='(I5)', initThresholdValue_max_a1)

      wRow3Base = WIDGET_BASE(wRightBase, /COLUMN)

      for i=0,2 do begin
        temp_a1[i] = WIDGET_BASE(wRow3Base, /COLUMN)
        Draw_a1[i] = WIDGET_DRAW(temp_a1[i], /BUTTON, $
          RETAIN=0, /WHEEL_EVENTS, $
          XSIZE=simageSize, YSIZE=simageSize)
        DrawSlider_a1[i] = WIDGET_SLIDER(temp_a1[i], $
          MINIMUM=1, /DRAG, $
          MAXIMUM=simageSize, VALUE=simageSize/2)
      endfor
      WIDGET_CONTROL, Draw_a1[0], SET_UVALUE = 'DRAW0_A1', SET_UNAME='d_surfview1:draw0_A1'
      WIDGET_CONTROL, DrawSlider_a1[0], SET_UVALUE = 'DRAWSLIDER0_A1', SET_UNAME='d_surfview1:drawSlider0_A1'
      WIDGET_CONTROL, Draw_a1[1], SET_UVALUE = 'DRAW1_A1', SET_UNAME='d_surfview1:draw1_A1'
      WIDGET_CONTROL, DrawSlider_a1[1], SET_UVALUE = 'DRAWSLIDER1_A1', SET_UNAME='d_surfview1:drawSlider1_A1'
      WIDGET_CONTROL, Draw_a1[2], SET_UVALUE = 'DRAW2_A1', SET_UNAME='d_surfview1:draw2_A1'
      WIDGET_CONTROL, DrawSlider_a1[2], SET_UVALUE = 'DRAWSLIDER2_A1', SET_UNAME='d_surfview1:drawSlider2_A1'
    
      IF N_ELEMENTS(additionalDir) GE 3 THEN BEGIN
        image_read, additionalDir[2], vol_HU_cube_a2, vol_HU_cube_resolution2, patient_age2, direction2, origin_loc2, p_modal = p_modal

        temp=volume_fusion(vol_HU_cube, vol_HU_cube_resolution, direction, origin_loc, vol_HU_cube_a2, vol_HU_cube_resolution2, direction2, origin_loc2)
        vol_HU_cube_a2 = temp
        IF MAX(temp)-MIN(temp) GT 0.001 THEN BEGIN
          ;this is MRI or data with correct fusion
          svol_HU_cube_a2 = CONGRID(temp, simageSize, simageSize, simageSize)
        ENDIF ELSE BEGIN
          ;this is PET/CT data, there is no need to do fusion
          ;GJ 2020/4/17
          svol_HU_cube_a2 = CONGRID(vol_HU_cube2, simageSize, simageSize, simageSize)
        ENDELSE

        ;calculate the additional image's threshold GJ 2020/4/19
        sizeVHc = SIZE(svol_HU_cube_a2)
        mean_image = IMAGE_THRESHOLD(BYTSCL(svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
        image_statistics, svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
        initThresholdValue_min_a2 = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
        
        wRow4Base = WIDGET_BASE(wRightBase, /COLUMN)

        for i=0,2 do begin
          temp_a2[i] = WIDGET_BASE(wRow4Base, /COLUMN)
          Draw_a2[i] = WIDGET_DRAW(temp_a2[i], /BUTTON, $
            RETAIN=0, /WHEEL_EVENTS, $
            XSIZE=simageSize, YSIZE=simageSize)
          DrawSlider_a2[i] = WIDGET_SLIDER(temp_a2[i], $
            MINIMUM=1, /DRAG, $
            MAXIMUM=simageSize, VALUE=simageSize/2)
        endfor
        WIDGET_CONTROL, Draw_a2[0], SET_UVALUE = 'DRAW0_A2', SET_UNAME='d_surfview1:draw0_A2'
        WIDGET_CONTROL, DrawSlider_a2[0], SET_UVALUE = 'DRAWSLIDER0_A2', SET_UNAME='d_surfview1:drawSlider0_A2'
        WIDGET_CONTROL, Draw_a2[1], SET_UVALUE = 'DRAW1_A2', SET_UNAME='d_surfview1:draw1_A2'
        WIDGET_CONTROL, DrawSlider_a2[1], SET_UVALUE = 'DRAWSLIDER1_A2', SET_UNAME='d_surfview1:drawSlider1_A2'
        WIDGET_CONTROL, Draw_a2[2], SET_UVALUE = 'DRAW2_A2', SET_UNAME='d_surfview1:draw2_A2'
        WIDGET_CONTROL, DrawSlider_a2[2], SET_UVALUE = 'DRAWSLIDER2_A2', SET_UNAME='d_surfview1:drawSlider2_A2'
      ENDIF
    ENDIF  
  ENDIF
  
  
  Draw[3] = WIDGET_DRAW(wRow1Base, $
    GRAPHICS_LEVEL=2, $
    XSIZE=imageSize*2.4, YSIZE=imageSize*2.4, /BUTTON_EVENTS, $
    UVALUE='DRAW3', RETAIN=0, /EXPOSE_EVENTS, /WHEEL_EVENTS, $
    UNAME='d_surfview1:draw3')
  wDraw = Draw[3]
  

  ;@GJ, 2023/5/16, display the orth 3D fusion  
  ;wRow5Base = WIDGET_BASE(wRightBase)
  wfusionViewBase =  WIDGET_BASE(wRightBase, XSIZE=imageSize*3.2, YSIZE=imageSize*2.2, /COLUMN)
  temp_draw3dView1 = WIDGET_BASE(wfusionViewBase, XSIZE=imageSize*3., YSIZE=imageSize, /ROW)
  temp_draw3dView2 = WIDGET_BASE(wfusionViewBase, XSIZE=imageSize*3., YSIZE=imageSize, /ROW)
  temp_3dView = LONARR(6)
  draw_3dView = LONARR(6)
  for i=0,2 do begin
    temp_3dView[i] = WIDGET_BASE(temp_draw3dView1, XSIZE=imageSize, YSIZE=imageSize)
    Draw_3dView[i] = WIDGET_DRAW(temp_3dView[i], /BUTTON, $
      RETAIN=0, /WHEEL_EVENTS, $
      XSIZE=imageSize, YSIZE=imageSize)
  endfor
  for i=3,5 do begin
    temp_3dView[i] = WIDGET_BASE(temp_draw3dView2, XSIZE=imageSize, YSIZE=imageSize)
    Draw_3dView[i] = WIDGET_DRAW(temp_3dView[i], /BUTTON, $
      RETAIN=0, /WHEEL_EVENTS, $
      XSIZE=imageSize, YSIZE=imageSize)
  endfor
  WIDGET_CONTROL, Draw_3dView[0], SET_UVALUE = 'DRAW3_3DView0', SET_UNAME='d_surfview1:draw3_3dView0'
  WIDGET_CONTROL, Draw_3dView[1], SET_UVALUE = 'DRAW3_3DView1', SET_UNAME='d_surfview1:draw3_3dView1'
  WIDGET_CONTROL, Draw_3dView[2], SET_UVALUE = 'DRAW3_3DView2', SET_UNAME='d_surfview1:draw3_3dView2'
  WIDGET_CONTROL, Draw_3dView[3], SET_UVALUE = 'DRAW3_3DView3', SET_UNAME='d_surfview1:draw3_3dView3'
  WIDGET_CONTROL, Draw_3dView[4], SET_UVALUE = 'DRAW3_3DView4', SET_UNAME='d_surfview1:draw3_3dView4'
  WIDGET_CONTROL, Draw_3dView[5], SET_UVALUE = 'DRAW3_3DView5', SET_UNAME='d_surfview1:draw3_3dView5'
  
  IF N_ELEMENTS(additionalDir) GE 2 THEN widget_control, wfusionViewBase, map=1 ELSE widget_control, wfusionViewBase, map=0
  
  
  wDraw_3DView = Draw_3DView
  ;GJ 2019/1/29
  
  wHotKeyReceptor = WIDGET_TEXT(wRow1Base, $
    /ALL_EVENTS, $
    UVALUE='HOTKEY', $
    UNAME='d_surfview1:hotkey')

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
  sText = demo_gettips('./lib/source_code/surfview1.tip', wTopBase, wStatusBase)

  ;  Grab the window id of the drawable.
  ;
  WIDGET_CONTROL, wDraw, GET_VALUE=oWindow
  
  WIDGET_CONTROL, wDraw_3DView[0], GET_VALUE=oWindow_3DView ;GJ 2019/1/29 for naked eye 3D view

  WIDGET_CONTROL, wTopBase, SENSITIVE=0


  bias = -0.5
  ;  Compute viewplane rectangle based on aspect ratio.
  ;
  aspect = 1. ;float(xdim)/float(ydim)
  if (aspect > 1) then $
    myview = [(1.0-aspect)/2.0+bias, 0.0+bias, aspect, 1.0] $
  else $
    myview = [0.0+bias, (1.0-(1.0/aspect))/2.0+bias, 1.0, (1.0/aspect)]

  ;  Create view 
  ;  
  ;  
  ;  
  ;  object.
  ;
  oView = OBJ_NEW('IDLgrView', PROJECTION=2, EYE=3, $
    ZCLIP=[1.5,-1.5], VIEW=myview, COLOR=[0, 0, 0]) ;

  ;  Make the text location to be centered.
  ;
  textLocation = [myview[0]+0.5*myview[2], myview[1]+0.5*myview[3]]

  ;  Create and display the PLEASE WAIT text.
  ;
  oFont = OBJ_NEW('IDLgrFont', 'Helvetica', SIZE=18)
  oText = OBJ_NEW('IDLgrText', $
    'MPI+ Starting up  Please wait...', $
    ALIGN=0.5, $
    LOCATION=textLocation, $
    COLOR=[255,255,0], FONT=oFont)


  ruler1 = obj_new('IDLgrAXis',0,range=[0,40],thick=1,$
    color=[255,255,0],major=5,minor=9,TICKDIR = 0,extend=1)
  ruler2 = obj_new('IDLgrAXis',1,range=[0,40],thick=1,$
    color=[255,255,0],major=5,minor=9,TICKDIR = 0,extend=1)
 
  ckr1=obj_new('IDLgrText')

  
  ckr1.SetProperty,color=[255,102,102],locations=[-0.45,-0.41,0],strings='Distance = ',hide=1;,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]

  ckr=obj_new('IDLgrText')
  ckr.SetProperty,color=[255,102,102],locations=[-0.45,-0.41,0],strings='Distance = ',hide=1;,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]
  
  ckg=obj_new('IDLgrText')
  ckg.SetProperty,color=[153,255,0],locations=[-0.45,-0.435,0],strings='Distance = ',hide=1;,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]
  
  ckb=obj_new('IDLgrText')
  ckb.SetProperty,color=[0,153,255],locations=[-0.45,-0.46,0],strings='Distance = ',hide=1;,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]
  
  ckag=obj_new('IDLgrText')
  ckag.SetProperty,locations=[-0.45,-0.485,0],strings='Angle = ',hide=1;,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]

  ;;;;;;;;;;;;;;;;;;;;;;;hjx0912

  vol_white=0.0
;  vol_red=0.0
;  vol_green=0.0
  
  Stringwhite = STRING(vol_white,FORMAT='("Volume =", F9.3, " cm3")')
;  Stringred = STRING(vol_red,FORMAT='("vol(red) =", F9.3, " cm3")')
;  Stringgreen = STRING(vol_green,FORMAT='("vol(green) =", F9.3, " cm3")')

  tvol_white=obj_new('IDLgrText')
  tvol_white->SetProperty,color=[0,0,0],locations=[0.16,-0.485,0],strings=stringwhite,hide=1;,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]
;  tvol_red=obj_new('IDLgrText')
;  tvol_red.SetProperty,color=[255,0,0],locations=[0.16,-0.435,0],strings=stringred,hide=0,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]

;  tvol_green=obj_new('IDLgrText')
;  tvol_green.SetProperty,color=[0,100,0],locations=[0.16,-0.46,0],strings=stringgreen,hide=0,/FILL_BACKGROUND,FILL_COLOR=[255,255,255]


  ;;;;;;;;;;;;;;;;;;;;;;;hjx0912


  ;  Create model.
  ;
  oTopModel = OBJ_NEW('IDLgrModel')
  oRuler1Model = OBJ_NEW('IDLgrMOdel')
  oRuler2Model = OBJ_NEW('IDLgrMOdel')
  obottomModel = OBJ_NEW('IDLgrModel')
  oScalingModel = OBJ_NEW('IDLgrModel')
  

  oTopModel->Add, oScalingModel



  obottommodel->Add, ckr1
  obottommodel->Add, ckr
  obottommodel->Add, ckg
  obottommodel->Add, ckb
  obottommodel->Add,ckag
  oRuler1Model->Add,ruler1
  oRuler2Model->Add,ruler2
  ;;;;;;;;;;;;;;;;;;;;;;;hjx0912
;  obottommodel->Add,tvol_red
  obottommodel->Add,tvol_white
;  obottommodel->Add,tvol_green
  ;;;;;;;;;;;;;;;;;;;;;;;hjx0912
  
  oRuler1Model->scale,0.008,0.05,1
  oRuler1Model->translate,-0.2,0.4,0.45
  oRuler2Model->scale,0.05,0.008,1
  oRuler2Model->translate,0.4,-0.2,0.45

  ;  To avoid conflict with the surface, trace outlines will be
  ;  offset slightly.
  ;
  oTraceScalingModel = OBJ_NEW('IDLgrModel')
  oTraceRotationModel = OBJ_NEW('IDLexRotator', $
    [xdim/2.0, ydim/2.0], ydim/2.0)
  oTraceOffset = OBJ_NEW('IDLgrModel')

  oTopModel->Add, oTraceOffset
  oTraceOffset->Add, oTraceScalingModel
  oTraceScalingModel->Add, oTraceRotationModel

  oTraceOffset->Translate, 0, 0, .01 ; Value determined empirically.

  ;  Place the model in the view.
  ;
  oView->Add, oTopModel
  oView->Add,obottommodel
  oView->Add,oRuler1Model
  oView->Add,oRuler2Model
  
 
  ;end hjx
;  
;  oTopModel = OBJ_NEW('IDLgrModel')
;  oScalingModel = OBJ_NEW('IDLgrModel')
;  oRotationModel = OBJ_NEW('IDLexRotator', $
;    [xdim/2.0, ydim/2.0], ydim/2.0)
;
;  oTopModel->Add, oScalingModel
;  oScalingModel->Add, oRotationModel
;
;  ;  To avoid conflict with the surface, trace outlines will be
;  ;  offset slightly.
;  ;
;  oTraceScalingModel = OBJ_NEW('IDLgrModel')
;  oTraceRotationModel = OBJ_NEW('IDLexRotator', $
;    [xdim/2.0, ydim/2.0], ydim/2.0)
;  oTraceOffset = OBJ_NEW('IDLgrModel')
;
;  oTopModel->Add, oTraceOffset
;  oTraceOffset->Add, oTraceScalingModel
;  oTraceScalingModel->Add, oTraceRotationModel
;
;  oTraceOffset->Translate, 0, 0, .01 ; Value determined empirically.
;
;  ;  Place the model in the view.
;  ;
;  oView->Add, oTopModel

  ;  Scale the top model to fit the viewing area.
  ;
  sct = 0.6
  oTopModel->Scale, sct, sct, sct

  oTopModel->Add, oText
  
;  ;GJ 2019/1/28
;  axis = [0,1,0]
;  angle = 90
;  oScalingModel->Rotate, axis, angle
;  ;end GJ 2019/1/28

  ;  Draw the starting up screen.
  ;
  owindow.setproperty,graphics_tree=oview
  oWindow->Draw, oView
  
  ;婵＄偛楠哥槐锕傚冀閻熺増瀚叉俊鐐插濠㈡﹢宕㈠顒�唺
  ;;;HJX 5305-5321
  ;婵＄偛楠哥槐锕傚冀閻熺増瀚叉俊鐐插濠㈡﹢宕㈠顒�唺
  dataXYZ0 = DBLARR(3)*0.
  dataXYZ1 = DBLARR(3)*0.
  dataXYZ2 = DBLARR(3)*0.
  
  ;;;;;;;;;;;;;;hjx05252
  datasum0=0.
  datasum1=0.
  datasum2=0.
  ;;;;;;;;;;;;;;hjx05252
  
   arr1=fltarr(4,5)
   arr2=string(arr1,format='(f4.2)')
   arr2=reform(arr2,4,5)
   MDATAXYZDISTANCE=string('1',format='(f4.2)')
   
   dataXYZdistance_old = DISTANCE_MEASURE([[dataXYZ1], [dataXYZ2]]) * vol_HU_cube_resolution
   dataXYZdistance = DISTANCE_MEASURE([[dataXYZ1], [dataXYZ2]]) * vol_HU_cube_resolution
   dataXYZdistance1 = DISTANCE_MEASURE([[dataXYZ1], [dataXYZ2]]) * vol_HU_cube_resolution
   dataXYZdistance2 = DISTANCE_MEASURE([[dataXYZ1], [dataXYZ2]]) * vol_HU_cube_resolution
   
   dataXYZind = 0
   COUNT123=0
    

    
  ;  Compute coordinate conversion to normalize.
  ;
  sz = SIZE(z)
  maxx = sz[1] - 1
  maxy = sz[2] - 1
  maxz = MAX(z, MIN=minz)
  xs = [0+bias, 1.0/maxx]
  ys = [0+bias, 1.0/maxy]
  minz2 = minz - 1
  maxz2 = maxz + 1
  zs = [-minz2/(maxz2-minz2)+bias, 1.0/(maxz2-minz2)]

  ;  For height-fields, use the following vertex colors.
  ;
  vertexColors = BYTARR(3, sz[1]*sz[2], /NOZERO)
  cbins= $
    [[255,   0, 0],$
    [255,  85, 0],$
    [255, 170, 0],$
    [255, 255, 0],$
    [170, 255, 0],$
    [85,  255, 0],$
    [0,   255, 0]]

  zi = ROUND(z/float(maxz) * 6.0)
  vertexColors[*, *] = cbins[*, zi]

  oPalette = OBJ_NEW('IDLgrPalette')
  oPalette->LOADCT, 13

  oTextureImage = OBJ_NEW('IDLgrImage', $
    BYTSCL(REBIN(z,256,256)), PALETTE=oPalette)

  ;  Create the tracing objects.
  ;
  workData = BYTARR((size(z))[1]*10, (size(z))[2]*10, 4) + 255b
  oTracingMask = OBJ_NEW('IDLgrImage', workData, INTERLEAVE=2)
  tracingData = FLTARR(3, 1024)
  tracingConnectivityList = LONARR(1024)
  tracingConnectivityList[0] = 0
  tracingConnectivityList[1] = -1

  ;  Create polyline object in the same space as the surface.
  ;
  oTracePolyline = OBJ_NEW('IDLgrPolyline', $
    tracingData, POLYLINES=tracingConnectivityList, $
    COLOR=[255, 0, 0], $
    XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs, $
    THICK=3)
  oTraceRotationModel->Add, oTracePolyline

  
  ;婵烇綀顕ф慨鐐村閻樻彃寮抽柣銊ュ煢oly
;;;;;;;;;hxn-start;;;;;;;;;;
  ;oPoly = OBJ_NEW('IDLgrPolygon', COLOR =sel_one.color, sel_one.verts, POLYGON=sel_one.conn, STYLE=2, $
   ;   SHADING=1, $
    ;  XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV)
 ; oRotationModel->Add, oPoly
  ;;;;;;;;;hxn-end;;;;;;;;;;

  
  s = obj_new('orb',radius=0.1,shading=1)
  s->getproperty, POLYGONS=conn
  s->getproperty, DATA=verts
  
  oS0 = obj_new('IDLgrPolygon',verts,poly=conn,color=[255,0,0],$
    shading=0, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, ALPHA_CHANNEL = 0.)
  oRotationModel->Add, oS0
  oS1 = obj_new('IDLgrPolygon',verts,poly=conn,color=[0,255,0],$
    shading=0, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, ALPHA_CHANNEL = 0.)
  oRotationModel->Add, oS1

  oS2 = obj_new('IDLgrPolygon',verts,poly=conn,color=[0,0,255],$
    shading=0, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, ALPHA_CHANNEL = 0.)
  oRotationModel->Add, oS2


  oLine1 = obj_new('IDLgrPolyline', [[dataxyz0], [dataxyz1]], color = [255, 0, 0], $
    shading=0, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, ALPHA_CHANNEL = 0.)
  oRotationModel->Add, oLine1


  oLine2 = obj_new('IDLgrPolyline', [[dataxyz1], [dataxyz2]], color = [0, 255, 0], $
    shading=0, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, ALPHA_CHANNEL = 0.)
  oRotationModel->Add, oLine2


  oLine3 = obj_new('IDLgrPolyline', [[dataxyz2], [dataxyz0]], color = [0, 0, 255], $
    shading=0, XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, ALPHA_CHANNEL = 0.)
  oRotationModel->Add, oLine3

  
  ;Create axes objects
  ;
  stick = [-.5, .5] * 1.5
  oRedAxis = OBJ_NEW('IDLgrPolyline', $
    stick, [0,0], [0,0], COLOR=[255,0,0], /HIDE)
  oGreenAxis = OBJ_NEW('IDLgrPolyline', $
    [0,0], stick, [0,0], COLOR=[0,255,0], /HIDE)
  oBlueAxis = OBJ_NEW('IDLgrPolyline', $
    [0,0], [0,0], stick, COLOR=[0,0,255], /HIDE)
 ;;;;;;;;;;;;;;;;;;;;;;;;hjx0713

 ogreenaxis1 = OBJ_NEW('IDLgrPolyline', $
  XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, COLOR=[0,255,0], /HIDE)
 oyellowaxis1 = OBJ_NEW('IDLgrPolyline', $
   XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, COLOR=[255,255,0], /HIDE)
 oBlueAxis1 = OBJ_NEW('IDLgrPolyline', $
   XCOORD_CONV=sel_one.XCOORD_CONV, YCOORD_CONV=sel_one.YCOORD_CONV, ZCOORD_CONV=sel_one.ZCOORD_CONV, COLOR=[0,0,255], /HIDE)
 oRotationModel->Add, ogreenaxis1
 oRotationModel->Add, oyellowaxis1
 oRotationModel->Add, oBlueAxis1
 ;;;;;;;;;;;;;;;;;;;;;;;hjx0713
    
    
  oRotationModel->Add, oRedAxis
  oRotationModel->Add, oGreenAxis
  oRotationModel->Add, oBlueAxis

  ;  Create a light.
  ;
;  oSunLight = OBJ_NEW('IDLgrLight', LOCATION=[1.5, 0, 1], $
;    TYPE=1, INTENSITY=0.5)
;  otopModel->Add, oSunLight
;  oSunLight = OBJ_NEW('IDLgrLight', TYPE=0, $
;    INTENSITY=0.75)
;  otopModel->Add, oSunLight
  
  ;婵烇綀顕ф慨鐐寸▔閵堝嫰鍤嬮柣蹇ｅ灠閸橈拷
  oTopModel->add, obj_new('IDLgrLight', $
    loc=[2,2,5], $
    type=2, $ ; Directional (parallel rays).
    color=[255,255,255], $
    intensity=.5 $
    )
  oTopModel->add, obj_new('IDLgrLight', $
    type=0, $ ; Ambient.
    intensity=.5, $
    color=[255,255,255] $
    )

  ; Rotate to standard view for first draw.
  ;
     ;;;;;;hxn-start;;;;;
 oRotationModel->Rotate,      [1, 0, 0], -90
  oRotationModel1->Rotate,      [1, 0, 0], -90
  oTraceRotationModel->Rotate, [1, 0, 0], -90
  oRotationModel->Rotate, [0, 1, 0], 30
  oRotationModel->Rotate, [1, 0, 0], 30
  oRotationModel1->Rotate, [0, 1, 0], 30
  oRotationModel1->Rotate, [1, 0, 0], 30
  
  oTraceRotationModel->Rotate, [0, 1, 0], 30
  oTraceRotationModel->Rotate, [1, 0, 0], 30
  
  oScalingModel->Add, oRotationModel
  oScalingModel->Add, oRotationModel1
   ;;;;;;hxn-end;;;;;


  oRotationModel->GetProperty, TRANSFORM=initTransformRotation
  oScalingModel->GetProperty, TRANSFORM=initTransformScaling
  
  ;GJ 2019/2/24, add color bar for thickness measurment
  colorbarmodel = OBJ_NEW('IDLgrModel')
  mytitle = OBJ_NEW('IDLgrText', 'Thickness')
  barDims = [0.02, 0.2]
  LOADCT, 13, RGB_TABLE = rgb_table
  redValues = rgb_table[*,0]
  greenValues = rgb_table[*,1]
  blueValues = rgb_table[*,2]
  mycolorbar = OBJ_NEW('IDLgrColorbar', redValues, $
    greenValues, blueValues, TITLE=mytitle, $
    DIMENSIONS=barDims);, /SHOW_AXIS, /SHOW_OUTLINE)
  colorbarmodel->Add, mycolorbar
  oView->Add, colorbarmodel
  ; Center the colorbar in the window.
  ; Note that you must use the ComputeDimensions method to
  ; get the dimensions of the colorbar.
  barPlusTextDims = mycolorbar->ComputeDimensions(oWindow)
  colorbarmodel->Translate, -barDims[0]-(barPlusTextDims[0]*20.), $
    -barDims[1]+(barPlusTextDims[1]*2.), 0
  colorbarmodel->setProperty, Hide = 1

  if not keyword_set(record_to_filename) then $
    record_to_filename = ''

  ;  Create the xd_sState
  ;
  xd_sState = { $
    icon_dir:icon_dir,$
    ;;;HJX ;;;;;;;;;;;;;;;;;;;;;;;;;
    arr1:arr1,$
      arr2:arr2,$
      mdataXYZdistance:mdataXYZdistance, $
   ;;;HJX ;;;;;;;;;;;;;;;;;;;;;;;;;
    colorTable: colorTable, $
    center: xdim/2., $                           ; Center of view
    radius: ydim/2, $                            ; Radius
    xSizeData: sz[1], $                          ; x data dimension
    ySizeData: sz[2], $                          ; x data dimension
    btndown: 0b, $                               ; Botton down flag
    pt0: FLTARR(3), $                            ; Initial point
    pt1: FLTARR(3), $                            ; Final point
    dragq: 2, $                                  ; Drag quality
    firstPoint: FLTARR(3), $
    
    wDraw: wDraw, $                              ; Widget draw
    wDraw_3DView: wDraw_3DView, $               ;GJ 2019/1/29 ; for naked eye 3D view
    Draw_3DView: Draw_3DView, $                 ;GJ 2019/1/29 for naked eye 3D view
;    imageSize_3DView: imageSize_3DView, $       ;GJ 2019/1/29 for naked eye 3D view
    Draw: Draw, $
;    Draw_3d: Draw_3d, $
    Draw3_size: imageSize*2.4, $                             
    DrawSlider: DrawSlider, $
    additionalDir: additionalDir, $                
    draw_a1: draw_a1, $
    DrawSlider_a1: DrawSlider_a1, $
    draw_a2: draw_a2, $
    DrawSlider_a2: DrawSlider_a2, $
    vol_HU_cube_a1: vol_HU_cube_a1, $   ;Original is bright blood, a1 is for black blood images
    vol_HU_cube_a2: vol_HU_cube_a2, $   ;Original is bright blood, a1 is for black blood images with contrast agent
    svol_HU_cube_a1: svol_HU_cube_a1, $   ;;Original is bright blood, a1 is for black blood images
    svol_HU_cube_a2: svol_HU_cube_a2, $   ;Original is bright blood, a1 is for black blood images with contrast agent
    window: window, $                        
    imageSize: simageSize, $                  
    location3D_previous: location3D, $
    location3D: location3D, $            
    oTopModel: oTopModel, $                      ; Top model
    oRuler1Model:oRuler1Model ,$
    oRuler2Model:oRuler2Model  ,$            ;ruler model
    oScalingModel: oScalingModel, $
    oTraceScalingModel: oTraceScalingModel, $
    oRotationModel: oRotationModel, $
    oRotationModel1: oRotationModel1, $          ;;;;;;;hxn;;;;;
    oTraceRotationModel: oTraceRotationModel, $
    oImage1: oImage1, $     ;GJ, 2019/8/17
    oSurface: oSurface, $                        ; Surface object
    ;Outverts: sel_one.verts, $                              ;Polygon object GJ, 2017.6.16
    ;Outconn: sel_one.conn, $                              ;Polygon object GJ, 2017.6.16
  
    sel_one:sel_one, $                    
  
    ;lqrname: lqrname, $
    ;myindex_m: myindex_m, $
   
    struc_modifyobject: struc_modifyobject, $; 2019/2/20 
    cur_patient: cur_patient, $ ;2019/2/20
    oPoly: oPoly, $                              ;Polygon object GJ, 2017.6.16, now is an array, GJ 2019/2/21
    oPoly_size: oPoly_size, $                    ;Polygon size GJ, 2018.6.9
    oPoly_min: oPoly_min, $
    oPoly_max: oPoly_max, $
    oPoly_border: oPoly_border, $                ;Polygon border GJ, 2018.6.9
      
    oPoly1: oPoly1, $                              ;;;;;;;;;hxn;;;;;;;;;;;
    oPoly2: oPoly2, $                              ;;;;;;;;;hxn;;;;;;;;;;;
    clippoint:intarr(6)*0., $          ;;;;;;;;;hxn0702;;;;;;;;;;;
    oS0: oS0, $                                  
    oS1: oS1, $                                 
    ckr1: ckr1, $  
    ckr: ckr, $                                       ;text object
    ckg: ckg, $    
    ckb: ckb, $       
    ckag: ckag, $  
    ;;;;;;;;;;;;;;;;;;;;;;;hjx0912
;    tvol_red:tvol_red, $                               ;text object
;    tvol_green:tvol_green, $
    tvol_white:tvol_white, $


;    vol_red:vol_red, $                               ;number
;    vol_green:vol_green, $
    vol_white:vol_white, $

    ;;;;;;;;;;;;;;;;;;;;;;;hjx0912
    dataXYZdistance_old:dataXYZdistance_old, $
    
;    wloadctSlider:wloadctSlider, $
;    wloadctSlidervalue:wloadctSlidervalue, $
    
;    wloadctSlider1:wloadctSlider1, $
;    wloadctSlidervalue1:wloadctSlidervalue1, $
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;hjx2;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtransLabel:wtransLabel, $
    wtransSlider:wtransSlider, $
    wtransSlidervalue:wtransSlidervalue, $
    wTrans_text: wTrans_text, $
    
;    wtransLabel1:wtransLabel1, $
;    wtransSlider1:wtransSlider1, $
;    wtransSlidervalue1:wtransSlidervalue1, $
    ;;;;;;;;;;;;;;;;;;;;;;;;;;hjx2;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    wCOLButton: wCOLButton, $              ; Vertex coloring option
;    wblueButton: wblueButton, $
;    wpinkButton: wpinkButton, $
;    wsalmonButton: wsalmonButton, $
;    wwhiteButton: wwhiteButton, $


    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ;;;HJX  5526-5534     
    COUNT123:COUNT123, $                              
      oS2: oS2, $    
    oLine1: oLine1, $                                    
    oLine2: oLine2, $                                    
    oLine3: oLine3, $                               
    dataXYZ2: dataXYZ2, $                             
    dataXYZdistance1:dataXYZdistance1, $             
    dataXYZdistance2:dataXYZdistance2, $
    vol_HU_cube: vol_HU_cube, $   
    svol_HU_cube: svol_HU_cube, $
    vol_HU_cube_resolution: vol_HU_cube_resolution, $
    vol_HU_cube_resolution1: vol_HU_cube_resolution1, $
    mask_vol_cube: vol_HU_cube_mask, $; *0+1, $
    mask_vol_cube_modify: vol_HU_cube_mask_modify, $
    ori_mask_vol_cube: vol_HU_cube_mask, $
    new_mask_vol_cube: vol_HU_cube_mask, $    ;GJ 2018/5/26, for new result
    distance_vol_HU_cube: vol_HU_cube, $      ;GJ 2018/12/18, for distance map
    xyROIout: OBJ_NEW('IDLgrROI'), $  
    yzROIout: OBJ_NEW('IDLgrROI'), $  
    xzROIout: OBJ_NEW('IDLgrROI'), $  
    xy_mask: xy_mask, $
    xz_mask: xz_mask, $
    yz_mask: yz_mask, $
    xy_mask_ROI: xy_mask_ROI, $
    xz_mask_ROI: xz_mask_ROI, $
    yz_mask_ROI: yz_mask_ROI, $
    xy_mask_clip: xy_mask_clip, $
    xz_mask_clip: xz_mask_clip, $
    yz_mask_clip: yz_mask_clip, $
    dataXYZ0: dataXYZ0, $                
    dataXYZ1: dataXYZ1, $               
    dataXYZdistance:dataXYZdistance, $    
    dataXYZind: dataXYZind, $            
    initHoleValue: initHoleValue, $
    initThresholdValue_min: initThresholdValue_min, $      ;initial threshold value
    initThresholdValue_max: initThresholdValue_max, $      ;initial threshold value
    initThresholdValue_min_a1: initThresholdValue_min_a1, $      ;initial threshold value
    initThresholdValue_max_a1: initThresholdValue_max_a1, $      ;initial threshold value
    threshold_min_value_a1:initThresholdValue_min_a1, $
    threshold_max_value_a1:initThresholdValue_max_a1, $
    oRedAxis: oRedAxis, $
    oGreenAxis: oGreenAxis, $
    oBlueAxis: oBlueAxis, $
    ;;;;;;;;;;;;;;;;;;;hjx0713
    ogreenaxis1: ogreenaxis1, $
    oyellowaxis1: oyellowaxis1, $
    oBlueAxis1: oBlueAxis1, $ 
    ;;;;;;;;;;;;;;;;;;hjx0713
    
    initTransformScaling: initTransformScaling, $
    initTransformRotation: initTransformRotation, $
    oView: oView, $                              ; Main view object
    oFont: oFont, $                              ; Font object
    oText: oText, $                              ; Text object
    vertexColors: vertexColors, $                ; Vertex colors(RGB)
    oTracePolyline: oTracePolyline, $            ; Trace object
    oTracingMask: oTracingMask, $                ; Tracing mask object
    oTextureImage: oTextureImage, $              ; Texture image object
    oPalette: oPalette, $                        ; Palette for image
    tracingMode: 0, $                            ; 0=off,1=on
    oWindow: oWindow, $                          ; Window object
    oWindow_3DView: oWindow_3DView, $             ;GJ 2019/1/29 for naked eye 3D view
    colorbarmodel: colorbarmodel, $       ;GJ 2019/2/24 for thickness colorbar
    wTopBase : wTopBase, $                       ; Top level base
    wShadingButton : wShadingButton, $
    wFlatButton : wFlatButton, $                 ; Shading options
    wGouraudButton : wGouraudButton, $
    wStyleButton: wStyleButton, $
    wPointButton : wPointButton, $             ; Styles options
    wWireButton : wWireButton, $
    wSolidButton : wSolidButton, $
    wRuledXZButton : wRuledXZButton, $
    wRuledYZButton : wRuledYZButton, $
    wLegoWireButton : wLegoWireButton, $
    wLegoSolidButton : wLegoSolidButton, $
    wSkirtNoneButton : wSkirtNoneButton, $       ; Skirt options
    wSkirt10Button : wSkirt10Button, $
    wSkirt20Button : wSkirt20Button, $
    wSkirt30Button : wSkirt30Button, $
    wLowButton : wLowButton, $                   ; Drag quality options
    wMediumButton : wMediumButton, $
    wHighButton : wHighButton, $
    wSmooth: wSmooth, $
    wSmoothLabel:wSmoothLabel, $
    wSmooth_text:wSmooth_text, $
    wFusionModeRadio:wFusionModeRadio, $
    wxFlipModeRadio:wxFlipModeRadio, $
    wyFlipModeRadio:wyFlipModeRadio, $
    wzFlipModeRadio:wzFlipModeRadio, $
    wx_move_slider:wx_move_slider, $
    wy_move_slider:wy_move_slider, $
    wz_move_slider:wz_move_slider, $
    wxyzflipmove_output:wxyzflipmove_output, $
    wxyzflipmove_apply:wxyzflipmove_apply, $
    fusionmode:0, $
    xFlipMode:0, $
    yFlipMode:0, $
    zFlipMode:0, $
    x_Move:0, $
    y_Move:0, $
    z_Move:0, $
  
;    wSmooth1: wSmooth1, $
;    wSmoothLabel1:wSmoothLabel1, $
    
    wThresholdSlider_min: wThresholdSlider_min, $                    ; Sliders IDs
    wThresholdSlider_max: wThresholdSlider_max, $                    ; Sliders IDs
    wThresholdText_min: wThresholdText_min, $
    wThresholdText_max: wThresholdText_max, $                 ;GJ 2018/12/22
    minValue:minValue, $
    maxValue:maxValue, $
    ;    wHoleSlider: wHoleSlider, $                              ;Slider hole size ID
    wThresholdmaxLabel: wThresholdmaxLabel, $                      ; Slider Threshold ID
    wThresholdminLabel: wThresholdminLabel, $                      ; Slider Threshold ID
    
    wThresholdSlider_min_a1: wThresholdSlider_min_a1, $                    ; Sliders IDs
    wThresholdSlider_max_a1: wThresholdSlider_max_a1, $                    ; Sliders IDs
    wThresholdText_min_a1: wThresholdText_min_a1, $
    wThresholdText_max_a1: wThresholdText_max_a1, $                 ;GJ 2018/12/22
    minValue_a1:minValue_a1, $
    maxValue_a1:maxValue_a1, $
    ;    wHoleSlider: wHoleSlider, $                              ;Slider hole size ID
    wThresholdmaxLabel_a1: wThresholdmaxLabel_a1, $                      ; Slider Threshold ID
    wThresholdminLabel_a1: wThresholdminLabel_a1, $                      ; Slider Threshold ID
    wfusionViewBase:wfusionViewBase, $
;    wHoleLabel: wHoleLabel, $
    wScalingSlider: wScalingSlider, $
    wScalingLabel: wScalingLabel, $
    wObjectIndex: wObjectIndex, $   ;GJ 2019/2/21, object index ID
    wObjectIndex_text: wObjectIndex_text, $    ;showing object indes
    wColorIndex: wColorIndex, $   ;showing color index
    wColorIndex_draw: wColorIndex_draw, $       ;showing color index color
    ;wColorIndexID: wColorIndexID, $         ;color index window ID
    wTracingButton : wTracingButton, $           ; Tracing options
    wTracingModeButton : wTracingModeButton, $
    wTracingMaskButton : wTracingMaskButton, $
    wHiddenButton : wHiddenButton, $             ; Hidden option
    wHideAxes: wHideAxes, $
;    wAnimateButton: wAnimateButton, $
    wHotKeyReceptor: wHotKeyReceptor, $
;    wConstraintsDroplist: wConstraintsDroplist, $
;    wJZsegmentation: wJZsegmentation, $
;    wZGG: wZGG, $
;    wCKcomplete: wCKcomplete, $
;    wCKcompleteNot: wCKcompleteNot, $
;    wZB: wZB, $
    wTextureButton : wTextureButton, $           ; Texture mapping option
    wVertexButton: wVertexButton, $              ; Vertex coloring option
    sText: sText, $                              ; Text tips structure
    wLineStyleButton : wLineStyleButton, $       ; Linestyle option
    orig_except : !except, $
    record_to_filename: record_to_filename, $
    debug: keyword_set(debug), $
    groupBase: groupBase, $                       ; Base of Group Leader
    
    ;;;;;;;;hxn-start;;;;;;;;
    lmb_scale:0,$ 
    sc:FLTARR(3), $
;    wrdstl:wrdstl, $
    wwrstl:wwrstl, $
    wexpstl:wexpstl, $
;    w3DSeg:w3DSeg, $
    w3DSegd:w3DSegd, $
    w3Drg:w3Drg, $
    w3Drg_face:w3Drg_face, $
    wAIVesselSeg:wAIVesselSeg, $; GJ 2020/3/22
    wVesselEcc:wVesselEcc, $; GJ 2020/3/23
    wVWCESeg:wVWCESeg, $; GJ 2020/3/22
;    w3Drgd:w3Drgd, $
    wUndo:wUndo, $
    wRedo:wRedo, $
    wthickness:wthickness, $ ;GJ 2019/2/16
    ;GJ 2019/2/26
    wpoint: wpoint, $
    wwire: wwire, $
    wfill: wfill, $
;    whide: whide, $  ;for single hide
;    whidecancel: whidecancel, $
    wsinglecenter: wsinglecenter, $
    ;GJ 2019/2/26
    selected: OBJ_NEW('IDLgrPolygon'), $ ;GJ 2019/2/17
    ShowList:1, $ ;GJ 2019/2/17
    modify_obj_list: OBJARR(6), $
    modify_mask_list: PTRARR(6), $
    modify_index: 0, $
    modify_steps: 0, $
    seed_location:DBLARR(3)*0., $
    direction:'yz', $
    ; rVolume: rVolume,rVolume1: rVolume1,   $
    ;rScaleToys: rScaleToys,rScaleToys1: rScaleToys1, rModel:rModel, rView: rView, rWindow: rWindow, $

    ;originalvolumedata: originalvolumedata,    $
    ;    volumexlength: i[1],            $
    ;    volumeylength: i[2],            $
    ;    volumezlength: i[3],            $
;
;    hSlice_base: hSlice_base,                      $
;    hSlice_x_slider: hSlice_x_slider,              $
;    hSlice_y_slider: hSlice_y_slider,              $
;    hSlice_z_slider: hSlice_z_slider,              $
;
;    hSlice_x_slider1: hSlice_x_slider1,              $
;    hSlice_y_slider1: hSlice_y_slider1,              $
;    hSlice_z_slider1: hSlice_z_slider1,              $
;
;
;;    hSlice_on_x_button: hSlice_on_x_button,        $
;    hSlice_on_y_button: hSlice_on_y_button,        $
;;    hSlice_on_z_button: hSlice_on_z_button,        $
;    hSlice_on_x_button_status: 0,                           $
;    hSlice_on_y_button_status: 0,                           $
;    hSlice_on_z_button_status: 0, $
    flag_x:0, $
    flag_other:1 , $
    flag_select:0 , $
;    ;clip_plane_cur:clip_plane_cur,$
;    hSlice_on: hSlice_on,                          $
;    hSlice_on_status: 0,                                    $ ; 1=yes, 0=no
;    ;hSlice_save: hSlice_save  ,                     $
;    hSlice_output: hSlice_output  ,   $
;    hSlice_three:hSlice_three, $


    ;;;;;;;;hxn-end;;;;;;;
    
    ;;GJ 20181219
      red_mask_vol_cube:vol_HU_cube_mask*0 $
;      hSlice_base1:hSlice_base1, $
;      xyzscalemove_output:xyzscalemove_output, $
;      x_scale_label:x_scale_label, $
;      y_scale_label:y_scale_label, $
;      z_scale_label:z_scale_label, $
;      x_scale_slider:x_scale_slider, $
;      y_scale_slider:y_scale_slider, $
;      z_scale_slider:z_scale_slider, $
;      x_move_label:x_move_label, $
;      y_move_label:y_move_label, $
;      z_move_label:z_move_label, $
;      x_move_slider:x_move_slider, $
;      y_move_slider:y_move_slider, $
;      z_move_slider:z_move_slider $
    
    ;;GJ 20181219 end


  }

  ;
  ;
  ;
  WIDGET_CONTROL, wColorIndex_draw, GET_VALUE=wColorIndexID
;  DEVICE, DECOMPOSED=0
;;  WSET, wColorIndexID
;;  ; Create a set of R, G, and B colormap vectors:
;  R = BINDGEN(256)
;  G = BINDGEN(256)
;  B = BINDGEN(256)
;
;  ; Load these vectors into the color table:
;  TVLCT, R, G, B
  DEVICE, DECOMPOSED=1
  WSET, wColorIndexID
  color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
  ERASE, ColorIndexBGR[0];'000000'x;, nothing selected, and use black color
  
  ;update the modify list
  xd_sState.modify_obj_list[0] = OBJ_NEW('IDLgrPolygon', DATA = sel_one.verts, POLYGONS=sel_one.conn)
  ;    xobjview_XDU, xd_sState.modify_obj_list[xd_sState.modify_index MOD 5]
  xd_sState.modify_mask_list[0] = PTR_NEW(xd_sState.new_mask_vol_cube)
  widget_control, xd_sState.wUndo, SENSITIVE = 0
  widget_control, xd_sState.wRedo, SENSITIVE = 0
  
  ;update the volume of 3d regions
  index = WHERE(xd_sState.new_mask_vol_cube, count)
  IF count NE 0 THEN xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000.
  print, "Volume = ", xd_sState.vol_white, ' cm3'
  Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
  xd_sState.tvol_white.SetProperty,strings=stringwhite
  
  ;adjust rulers
  ;[MAX(sel_one.verts[0,*])-MIN(sel_one.verts[0,*])
  scale = 2.3*280./(oPoly_size[0]*xd_sState.vol_HU_cube_resolution)
  xd_sState.oRuler1model->setproperty,transform=[[scale, 0, 0, 0.0],[0, 1, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
  xd_sState.oRuler2model->setproperty,transform=[[1, 0, 0, 0.0],[0, scale, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
  xd_sState.oRuler1Model->scale,0.008,0.05,1
  xd_sState.oRuler1Model->translate,-0.2,0.4,0.45
  xd_sState.oRuler2Model->scale,0.05,0.008,1
  xd_sState.oRuler2Model->translate,0.4,-0.2,0.45
  
  ;GJ 2020/5/25
  ObjectIndex_value = 1
  ERASE, ColorIndexBGR[LONG(ObjectIndex_value)]
  print, 'color = ', ColorIndexBGR[LONG(ObjectIndex_value)]
  xd_sState.selected=xd_sState.oPoly[ObjectIndex_value-1]
  xd_sState.selected->GetProperty, DATA = DATA, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel, hide=hideStatus, THICK=thickness
  xd_sState.tvol_white.SetProperty,hide=0
  count = xd_sState.struc_modifyobject.countmask_list[LONG(ObjectIndex_value)-1]
  xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
  Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
  ;            IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
  IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color

  dataxyz = REFORM(DATA[*, N_ELEMENTS(DATA(0,*))/2])
  ;the object is multi-color, we don't need to do anything
  ;if the object is single color, we will find right color index
  IF STDDEV(vert_colors[0,*]) LT 1 THEN BEGIN
    ;              color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
    FOR m=0, N_ELEMENTS(ColorIndex)-1 DO BEGIN
      IF vert_colors[0,0] EQ color_rgb[m,0] AND vert_colors[1,0] EQ color_rgb[m,1] AND vert_colors[2,0] EQ color_rgb[m,2] THEN colorindexS = m
    ENDFOR
    IF N_ELEMENTS(colorindexS) EQ 0 THEN colorindexS = 0
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=colorindexS
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[colorindexS]
    IF hideStatus EQ 1 THEN BEGIN
      plots, [0,40], [0,30], color=[0,0,0], /device
      plots, [0,40], [30,0], color=[0,0,0], /device
    ENDIF
  ENDIF ELSE BEGIN

    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ;use multiple lines to label multiple colro
    plots, [0,40], [21,21], color=[109,255,109], THICK=3, /device
    plots, [0,40], [14,14], color=[255,255,255]*255, THICK=3, /device
    plots, [0,40], [7,7], color=[128,0,128], THICK=3, /device
    IF hideStatus EQ 1 THEN BEGIN
      plots, [0,40], [0,30], color=[0,0,0], /device
      plots, [0,40], [30,0], color=[0,0,0], /device
    ENDIF
  ENDELSE
  ;GJ 2020/5/25----end of selection
  
  ;  Unless we are debugging, accumulate any math warnings silently.
  ;
  !except = ([0, 2])[keyword_set(debug)]

  ;  Desensitize the defaults buttons.
  ;
  WIDGET_CONTROL, wGouraudButton, SENSITIVE=0
  WIDGET_CONTROL, wSolidButton, SENSITIVE=0
  WIDGET_CONTROL, wSkirtNoneButton, SENSITIVE=0
  WIDGET_CONTROL, wHighButton, SENSITIVE=0
  WIDGET_CONTROL, wHiddenButton, SENSITIVE=0
  WIDGET_CONTROL, wLineStyleButton, SENSITIVE=0
  WIDGET_CONTROL, wTracingMaskButton, SENSITIVE=0

;  help, *(xd_sState.sel_one.p_conn)
;  help, *(xd_sState.sel_one.p_verts)
;  smoothedOutverts = MESH_SMOOTH(*(xd_sState.sel_one.p_verts), *(xd_sState.sel_one.p_conn), ITERATIONS=50)

;  PRINT, 'p_conn; p_verts'


  WIDGET_CONTROL, wTopBase, /HOURGLASS

  otopModel->Remove, oText
  oWindow->Draw, oView
  
  
  WIDGET_CONTROL, wTopBase, SENSITIVE=1
  
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  ;    LOADCT, 0


  topClr = !D.TABLE_SIZE - 6
  ;  Load the grey scale colr table.
  ;
  LOADCT, 0, /SILENT, NCOLORS=topClr+1

  ;  Allocate working colors : red, green, yellow, blue, white.
  ;
  TVLCT, [255,0,255,0,0],[0,255,255,0,255],[0,0,0,255,255], topClr+1
  
    ; 闁活枎鎷风紓浣侯棞琚欓柛鎾寸墪濞达拷
    ;  Get the windows (drawing areas) identifiers.
    ;
    for i=0,2 do begin
      WIDGET_CONTROL, draw[i], GET_VALUE=j
      window[i] = j
    endfor
    
    mask_vol_cube = vol_HU_cube_mask
    
    imageSize = simageSize

    xyImage = DBLARR(imageSize, imageSize)
    temp_xyImage = svol_HU_cube[*, *, (imageSize/2)]
    xyImage[*, *] = temp_xyImage[*, *, 0]
    WSET, window[0]
    xyMask = CONGRID(mask_vol_cube[*, *, (SIZE(mask_vol_cube))[3]/2], imageSize, imageSize, 1)
    TV, REVERSE(BYTSCL(xyImage, TOP = topClr, MAX=initThresholdValue_min, MIN=MIN([-300, initThresholdValue_min-300])), 2)
    CONTOUR, REVERSE(xyMask,2), LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
    FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
      line = [LINDGEN(PathInfo(I).N), 0]
      oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
      DRAW_ROI, oROI, COLOR =topClr+1, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
      OBJ_DESTROY, oROI
    ENDFOR
    plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
    plots, [0,imageSize], [imageSize - location3D[1],imageSize - location3D[1]], COLOR=topClr+3, /DEVICE
    

    yzImage = DBLARR(imageSize, imageSize)
    temp_yzImage = svol_HU_cube[(imageSize/2), *, *]
    yzImage[*, *] = temp_yzImage[0, *, *]
    WSET, window[1]
    yzMask = CONGRID(mask_vol_cube[(SIZE(mask_vol_cube))[1]/2, *, *], 1, imageSize, imageSize)
    TV, BYTSCL(yzImage, TOP = topClr, MAX=initThresholdValue_min, MIN=MIN([-300, initThresholdValue_min-300]))
    CONTOUR, yzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
    FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
      line = [LINDGEN(PathInfo(I).N), 0]
      oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
      DRAW_ROI, oROI, COLOR =topClr+2, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
      OBJ_DESTROY, oROI
    ENDFOR
    plots, [location3D[1],location3D[1]], [0,imageSize], COLOR=topClr+3, /DEVICE
    plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE
    
    xzImage = DBLARR(imageSize, imageSize)
    xzMask = DBLARR(imageSize, imageSize)
    temp_xzImage = svol_HU_cube[*, (imageSize/2), *]
    xzImage[*, *] = temp_xzImage[*, 0, *]
    WSET, window[2]
    temp_xzMask = CONGRID(mask_vol_cube[*, (SIZE(mask_vol_cube))[2]/2, *], imageSize, 1, imageSize)
    xzMask[*, *] = temp_xzMask[*, 0, *]
    TV, BYTSCL(xzImage, TOP = topClr, MAX=initThresholdValue_min, MIN=MIN([-300, initThresholdValue_min-300]))
    CONTOUR, xzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
    FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
      line = [LINDGEN(PathInfo(I).N), 0]
      oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
      DRAW_ROI, oROI, COLOR =topClr+3, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
      OBJ_DESTROY, oROI
    ENDFOR
   plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
    plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE
    
    
    IF KEYWORD_SET(additionalDir) THEN BEGIN
      IF N_ELEMENTS(additionalDir) GE 2 THEN BEGIN
        for i=0,2 do begin
          WIDGET_CONTROL, draw_a1[i], GET_VALUE=j
          window[i] = j
        endfor

        xyImage = DBLARR(imageSize, imageSize)
        temp_xyImage = svol_HU_cube_a1[*, *, (imageSize/2)]
        xyImage[*, *] = temp_xyImage[*, *, 0]
        WSET, window[0]
        TV, HIST_EQUAL(REVERSE(BYTSCL(xyImage, TOP = topClr, MAX=initThresholdValue_min_a1, MIN=MIN([-300, initThresholdValue_min_a1-300])), 2))
        CONTOUR, REVERSE(xyMask,2), LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
        FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
          line = [LINDGEN(PathInfo(I).N), 0]
          oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
          DRAW_ROI, oROI, COLOR =topClr+1, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
          OBJ_DESTROY, oROI
        ENDFOR
        plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
        plots, [0,imageSize], [imageSize - location3D[1],imageSize - location3D[1]], COLOR=topClr+3, /DEVICE


        yzImage = DBLARR(imageSize, imageSize)
        temp_yzImage = svol_HU_cube_a1[(imageSize/2), *, *]
        yzImage[*, *] = temp_yzImage[0, *, *]
        WSET, window[1]
        TV, BYTSCL(yzImage, TOP = topClr, MAX=initThresholdValue_min_a1, MIN=MIN([-300, initThresholdValue_min_a1-300]))
        CONTOUR, yzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
        FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
          line = [LINDGEN(PathInfo(I).N), 0]
          oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
          DRAW_ROI, oROI, COLOR =topClr+2, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
          OBJ_DESTROY, oROI
        ENDFOR
        plots, [location3D[1],location3D[1]], [0,imageSize], COLOR=topClr+3, /DEVICE
        plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

        xzImage = DBLARR(imageSize, imageSize)
        temp_xzImage = svol_HU_cube_a1[*, (imageSize/2), *]
        xzImage[*, *] = temp_xzImage[*, 0, *]
        WSET, window[2]
        TV, BYTSCL(xzImage, TOP = topClr, MAX=initThresholdValue_min_a1, MIN=MIN([-300, initThresholdValue_min_a1-300]))
        CONTOUR, xzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
        FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
          line = [LINDGEN(PathInfo(I).N), 0]
          oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
          DRAW_ROI, oROI, COLOR =topClr+3, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
          OBJ_DESTROY, oROI
        ENDFOR
        plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
        plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE
        
        IF N_ELEMENTS(additionalDir) GE 3 THEN BEGIN
          for i=0,2 do begin
            WIDGET_CONTROL, draw_a2[i], GET_VALUE=j
            window[i] = j
          endfor

          xyImage = DBLARR(imageSize, imageSize)
          temp_xyImage = svol_HU_cube_a2[*, *, (imageSize/2)]
          xyImage[*, *] = temp_xyImage[*, *, 0]
          WSET, window[0]
          TV, REVERSE(BYTSCL(xyImage, TOP = topClr, MAX=initThresholdValue_min_a2, MIN=MIN([-300, initThresholdValue_min_a2-300])), 2)
          CONTOUR, REVERSE(xyMask,2), LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
          FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
            line = [LINDGEN(PathInfo(I).N), 0]
            oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
            DRAW_ROI, oROI, COLOR =topClr+1, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
            OBJ_DESTROY, oROI
          ENDFOR
          plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
          plots, [0,imageSize], [imageSize - location3D[1],imageSize - location3D[1]], COLOR=topClr+3, /DEVICE


          yzImage = DBLARR(imageSize, imageSize)
          temp_yzImage = svol_HU_cube_a2[(imageSize/2), *, *]
          yzImage[*, *] = temp_yzImage[0, *, *]
          WSET, window[1]
          TV, BYTSCL(yzImage, TOP = topClr, MAX=initThresholdValue_min_a2, MIN=MIN([-300, initThresholdValue_min_a2-300]))
          CONTOUR, yzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
          FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
            line = [LINDGEN(PathInfo(I).N), 0]
            oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
            DRAW_ROI, oROI, COLOR =topClr+2, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
            OBJ_DESTROY, oROI
          ENDFOR
          plots, [location3D[1],location3D[1]], [0,imageSize], COLOR=topClr+3, /DEVICE
          plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

          xzImage = DBLARR(imageSize, imageSize)
          temp_xzImage = svol_HU_cube_a2[*, (imageSize/2), *]
          xzImage[*, *] = temp_xzImage[*, 0, *]
          WSET, window[2]
          TV, BYTSCL(xzImage, TOP = topClr, MAX=initThresholdValue_min_a2, MIN=MIN([-300, initThresholdValue_min_a2-300]))
          CONTOUR, xzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS
          FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
            line = [LINDGEN(PathInfo(I).N), 0]
            oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
            DRAW_ROI, oROI, COLOR =topClr+3, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
            OBJ_DESTROY, oROI
          ENDFOR
          plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
          plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE
          
        ENDIF
        
      ENDIF
    ENDIF
     
  ;2019/2/3, GJ, plotting the naked 3D view window
  ;draw_3DView, xd_sState
    
  WIDGET_CONTROL, wTopBase, SET_UVALUE=xd_sState, /NO_COPY
  XMANAGER, 'd_surfview1', wTopBase, EVENT_HANDLER='d_surfview1Event', $
    CLEANUP='d_surfview1Cleanup', /NO_BLOCK

  WIDGET_CONTROL, wHotKeyReceptor, /INPUT_FOCUS

end

;; The program to calculate the thickess , jia.11@osu.edu
;1






PRO MorphThinMidline_2D, binaryImg, midline_vertices

  ; Prepare the display device and load grayscale color
  ; table.
  ;DEVICE, DECOMPOSED = 0, RETAIN = 2
  ;LOADCT, 0

  ;; Load an image.
  ;;file = FILEPATH('pollens.jpg', $
  ;;   SUBDIRECTORY = ['examples', 'demo', 'demodata'])
  ;tifai='K:\Knee_3T\2007_01_23_knee_volunteer_04_billy\Volume Measurements\cartilage\cartilage_SL042.tif'
  ;result=read_tiff(tifai, channels=0)
  ;img=result
  ;READ_JPEG, file, img, /GRAYSCALE

  ; Get the image size, prepare a display window and
  ; display the image.
  dims = SIZE(binaryImg, /DIMENSIONS)
  ;WINDOW, 0, XSIZE = 2*dims[0], YSIZE = 2*dims[1], $
  ;   TITLE = 'Original, Binary and Thinned Images'
  ;TVSCL, img, 0

  ;; Generate a binary image by thresholding.
  ;binaryImg = img GE 140
  ;TVSCL, binaryImg, 1

  ; Prepare hit and miss structures for thinning.
  h0 = [[0b, 0, 0], [0, 1, 0], [1, 1, 1]]
  m0 = [[1b, 1, 1], [0, 0, 0], [0, 0, 0]]
  h1 = [[0b, 0, 0], [1, 1, 0], [1, 1, 0]]
  m1 = [[0b, 1, 1], [0, 0, 1], [0, 0, 0]]
  h2 = [[1b, 0, 0], [1, 1, 0], [1, 0, 0]]
  m2 = [[0b, 0, 1], [0, 0, 1], [0, 0, 1]]
  h3 = [[1b, 1, 0], [1, 1, 0], [0, 0, 0]]
  m3 = [[0b, 0, 0], [0, 0, 1], [0, 1, 1]]
  h4 = [[1b, 1, 1], [0, 1, 0], [0, 0, 0]]
  m4 = [[0b, 0, 0], [0, 0, 0], [1, 1, 1]]
  h5 = [[0b, 1, 1], [0, 1, 1], [0, 0, 0]]
  m5 = [[0b, 0, 0], [1, 0, 0], [1, 1, 0]]
  h6 = [[0b, 0, 1], [0, 1, 1], [0, 0, 1]]
  m6 = [[1b, 0, 0], [1, 0, 0], [1, 0, 0]]
  h7 = [[0b, 0, 0], [0, 1, 1], [0, 1, 1]]
  m7 = [[1b, 1, 0], [1, 0, 0], [0, 0, 0]]

  ; Iterate until the thinned image is identical to
  ; the input image for a given iteration.
  bCont = 1b
  iIter = 1
  thinImg = binaryImg
  WHILE bCont EQ 1b DO BEGIN
    PRINT,'Iteration: ', iIter
    inputImg = thinImg

    ; Perform the thinning using the first pair
    ; of structure elements.
    thinImg = MORPH_THIN(inputImg, h0, m0)

    ; Perform the thinning operation using the
    ; remaining structural element pairs.
    thinImg = MORPH_THIN(thinImg, h1, m1)
    thinImg = MORPH_THIN(thinImg, h2, m2)
    thinImg = MORPH_THIN(thinImg, h3, m3)
    thinImg = MORPH_THIN(thinImg, h4, m4)
    thinImg = MORPH_THIN(thinImg, h5, m5)
    thinImg = MORPH_THIN(thinImg, h6, m6)
    thinImg = MORPH_THIN(thinImg, h7, m7)

    ; Display the results of thinning and wait a second for
    ; display purposes.
    ;   TVSCL, thinImg, 2
    ;   WAIT, 1

    ; Test the condition and increment the loop.
    bCont = MAX(inputImg - thinImg)
    iIter = iIter + 1

    ; End WHILE loop statements.
  ENDWHILE


  ;;;
  ;;;; the start of removing samll bars in midline

  ;Use Count to get the number of nonzero elements
  index=WHERE(thinImg, count)

  ; Only scbscript the array if it's safe
  IF count EQ 0 THEN  GOTO, convertVertices

  ;;;;;;do another search find the end points with incorrect midlines, and remove them
  ;;;;;;do another search find the end points
  FOR i=0, 10 DO BEGIN
    ; Show inverse of final result.
    TVSCL, 1 - thinImg
    ;WAIT, 1

    Calculate_connectivity, thinImg, dims, L, index, Midline_pixels
    ;calculate number of end point
    indEP = WHERE(L[*,1] LE 1, nEP)
    ;calculate number of branching point
    indBP = WHERE(L[*,1] GE 3, nBP)
    IF (nEP NE 0) AND (nBP NE 0) THEN thinImg[L[indEP]]=0.
  ENDFOR

  ;; Show inverse of final result.
  ;TVSCL, 1 - thinImg
  ;Calculate_connectivity, thinImg, dims, L, index, Midline_pixels
  ;;calculate number of end point
  ;indEP = WHERE(L[*,1] LE 1, nEP)
  ;;calculate number of branching point
  ;indBP = WHERE(L[*,1] GE 3, nBP)
  ;WHILE (nEP NE 0) AND (nBP NE 0) DO BEGIN
  ;        thinImg[L[indEP]]=0.
  ;        Calculate_connectivity, thinImg, dims, L, index, Midline_pixels
  ;        ;calculate number of end point
  ;        indEP = WHERE(L[*,1] LE 1, nEP)
  ;        ;calculate number of branching point
  ;        indBP = WHERE(L[*,1] GE 3, nBP)
  ;ENDWHILE
  ;window, 2
  ;TVSCL, 1 - thinImg

  ;;;;;;;do another search find the end points, and grow them based on old direction and
  ; Calculate_connectivity, thinImg, dims, L, index, Midline_pixels
  ;;calculate number of end point
  ;indEP = WHERE(L[*,1] EQ 1, nEP)
  ;IF (nEP NE 0) THEN BEGIN
  ;    FOR i=0L, nEP-1L DO BEGIN
  ;       pixel_1 = ARRAY_INDICES(thinImg, L[indEP[i], 0])
  ;       pixel_2 = ARRAY_INDICES(thinImg, L[indEP[i], 2])
  ;       pixel_2_Ind = WHERE(L[*,0] EQ L[indEP[i], 2])
  ;       IF L[pixel_2_Ind, 2] EQ L[indEP[i], 0] THEN BEGIN
  ;            pixel_3 = ARRAY_INDICES(thinImg, L[pixel_2_Ind, 3])
  ;       ENDIF ELSE BEGIN
  ;            pixel_3 = ARRAY_INDICES(thinImg, L[pixel_2_Ind, 2])
  ;       ENDELSE
  ;       deltaPixel_13 = pixel_1 - pixel_3
  ;       deltaPixel_12 = pixel_1 - pixel_2
  ;       newPixel_13 = pixel_1 + deltaPixel_13
  ;       newPixel_12 = pixel_1 + deltaPixel_12
  ;       WHILE binaryImg[newPixel_13[0], newPixel_13[1]] NE 0 DO BEGIN
  ;            thinImg[newPixel_13[0], newPixel_13[1]] = 1.
  ;            thinImg[newPixel_12[0], newPixel_12[1]] = 1.
  ;            midline_pixels = [[midline_pixels], [newPixel_12]]
  ;            midline_pixels = [[midline_pixels], [newPixel_13]]
  ;            newPixel_13 = newPixel_13 + deltaPixel_13
  ;            newPixel_12 = newPixel_13 + deltaPixel_12
  ;       ENDWHILE
  ;    ENDFOR
  ;ENDIF
  window, 3
  TVSCL, thinImg+binaryImg
  ;iimage, thinImg+binaryImg
  ;WAIT, 1
  ;
  Calculate_connectivity, thinImg, dims, L, index, Midline_pixels
  print, 'before removing the branching...', L

  ;calculate number of branching point
  indBP = WHERE(L[*,1] GE 3, nBP)

  ; Only scbscript the array if it's safe
  ; there is no branching point, exit
  IF nBP EQ 0 THEN GOTO, findlesion

  ;calculate the branching distances and store each path
  FOR j=0L, nBP-1L DO BEGIN
    numPath=L[indBP[j],1]
    Path_array=PTRARR(numPath, /ALLOCATE_HEAP)
    Path_length=LONARR(numPath)*0L

    FOR k=0L, numPath-1L DO BEGIN
      prevPixel=L[indBP[j],0]
      currPixel=L[indBP[j],k+2]
      Path_array[k]=PTR_NEW(currPixel)
      Path_length[k]=1L
      ;find the current pixel and find nearest pixels
      currInd = WHERE(L[*,0] EQ currPixel)
      WHILE (L[currInd, 1] EQ 2) DO BEGIN
        IF prevPixel EQ L[currInd, 2] THEN nextPixel= L[currInd, 3] ELSE nextPixel= L[currInd, 2]
        *Path_array[k]=[*Path_array[k], nextPixel]
        Path_length[k]=Path_length[k] + 1L
        prevPixel=currPixel
        currPixel=nextPixel
        currInd = WHERE(L[*,0] EQ currPixel[0])
      ENDWHILE
      IF (L[currInd, 1] NE 1) THEN Path_length[k] = 0 ;;IF EQ 0, end point is vertices; IF NE 0, end point is bifurcation point
    ENDFOR

    indZero=WHERE(Path_length EQ 0, countZero)
    IF countZero NE 0 THEN Path_length[indZero] = max(Path_length) + 1000L
    ;find the path with shortest length
    minlength = MIN(Path_length, minSubscript)
    ;assign the pixels in the shortest length as 0, assuming they are short bars, not the main midline of cartilage
    thinImg[*Path_array[minSubscript]]=BYTE(0)
    PTR_FREE, Path_array
  ENDFOR

  ;;;
  ;;;
  ;;;;; the end of removing samll bars in midline

  FINDLESION: print, 1
  convertVertices: index=WHERE(thinimg EQ 1, count)
  index_2D = ARRAY_INDICES(thinimg, index)
  midline_vertices = index_2D

  ; Show inverse of final result.
  ;TVSCL, 1 - thinImg
  ;WAIT, 1

  ;iimage, 1 - thinImg

  ;print, 1
END

;;闁跨喐鏋婚幏绌搃a.11@osu.edu, July 29, 2008, this is the key part of the program
;;;thinImg; the mask with only line segments
;;;L: is the connectivity of each point on line segments
PRO Calculate_connectivity, thinImg, dims, L, index, Midline_pixels

  ;;;;
  ;;;;
  ;Use Count to get the number of nonzero elements
  index=WHERE(thinImg, count)

  ; Only scbscript the array if it's safe
  IF count EQ 0 THEN  return

  ;Define a array to store the topology of the line
  L = LINDGEN(count, 6) ; [vertix index,  # nearst pixels, pixel a, pixel b, pixel c, pixel d]
  L = L * 0. - 1.
  L[*,0] = index

  ;Find out the 2D indices of midline pixels
  Midline_pixels=ARRAY_INDICES(thinImg, index)

  ;find the number of nereat pixels
  FOR i=0L, count-1L DO BEGIN
    nNP=0
    Pixel_0=Midline_pixels[0,i]-1
    Pixel_1=Midline_pixels[1,i]
    IF thinImg[Pixel_0, Pixel_1] EQ 1 THEN BEGIN
      nNP = nNP + 1
      L[i, nNP+1] = Pixel_1 * dims[0] + Pixel_0
    ENDIF
    Pixel_0=Midline_pixels[0,i]+1
    Pixel_1=Midline_pixels[1,i]
    IF thinImg[Pixel_0, Pixel_1] EQ 1 THEN BEGIN
      nNP = nNP + 1
      L[i, nNP+1] = Pixel_1 * dims[0] + Pixel_0
    ENDIF
    Pixel_0=Midline_pixels[0,i]
    Pixel_1=Midline_pixels[1,i]-1
    IF thinImg[Pixel_0, Pixel_1] EQ 1 THEN BEGIN
      nNP = nNP + 1
      L[i, nNP+1] = Pixel_1 * dims[0] + Pixel_0
    ENDIF
    Pixel_0=Midline_pixels[0,i]
    Pixel_1=Midline_pixels[1,i]+1
    IF thinImg[Pixel_0, Pixel_1] EQ 1 THEN BEGIN
      nNP = nNP + 1
      L[i, nNP+1] = Pixel_1 * dims[0] + Pixel_0
    ENDIF

    ;Also assign number of Nearest Pixel (nNP) to L
    ;nNP = 1, end point
    ;nNP = 2, connecting point
    ;nNP =3 or 4, branching point
    L[i, 1] = nNP
  ENDFOR

END


PRO ortho2D_plot, xd_sState
  
  ;@GJ, 2023/5/16, generating fused image
  ;@GJ, 2023/5/18, adding the option of fusion
  IF xd_sState.fusionmode EQ 1 THEN image_fusion_MPI_CT, xd_sState
  
  
  ;闁告帗绻傞‖濠囧礌閺嶎剦鍟庣紓鍐挎嫹
  location3D = xd_sState.location3D
  ;GJ 2019/6/1, negative location3D values may give error, bug was fixed.
  IF location3D[0] LT 0 THEN location3D[0] = 0
  IF location3D[1] LT 0 THEN location3D[1] = 0
  IF location3D[2] LT 0 THEN location3D[2] = 0
  
  ;update the location3D_previous
  IF TOTAL(ABS(location3D - xd_sState.location3D_previous)) GT 0.01 THEN BEGIN
    xd_sState.location3D_previous = location3D
  ENDIF
    
  ;print,'2D_plot location3D',location3D
  imageSize = xd_sState.imageSize
  ;print,'imagesize',imagesize
  IF MAX(location3D) GE imageSize THEN RETURN
  threshold_min_value = xd_sState.sel_one.HU_min_value
  threshold_max_value = xd_sState.sel_one.HU_max_value
  
  IF threshold_min_value GT threshold_max_value THEN BEGIN
    temp = threshold_min_value
    threshold_min_value = threshold_max_value
    threshold_max_value = temp
  ENDIF
  
  ;make sure the images are not too dark
  IF threshold_max_value GT ABS(threshold_min_value)+500 THEN threshold_max_value = ABS(threshold_min_value)+500
  threshold_min_value =threshold_min_value/3.-500.
  
  vol_HU_cube = xd_sState.vol_HU_cube
  WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=ObjectIndex_value
  IF LONG(ObjectIndex_value) EQ 0 THEN BEGIN
    mask_vol_cube = xd_sState.new_mask_vol_cube
;    print, 'MAX(mask_vol_cube) = ', MAX(mask_vol_cube)
;    print, 'MIN(mask_vol_cube) = ', MIN(mask_vol_cube)
  ENDIF ELSE BEGIN
    mask_vol_cube = xd_sState.new_mask_vol_cube*0
    mask_vol_cube[*(xd_sState.struc_modifyobject.p_maskindex_list[LONG(ObjectIndex_value)-1])] = 1
  ENDELSE
  svol_HU_cube = xd_sState.svol_HU_cube

  ; 闁活枎鎷风紓浣侯棞琚欓柛鎾寸墪濞达拷
  DEVICE, DECOMPOSED = 0, RETAIN = 2
  ;    LOADCT, 0


  topClr = !D.TABLE_SIZE - 6
  ;  Load the grey scale colr table.
  ;
  LOADCT, 0, /SILENT, NCOLORS=topClr+1

  ;  Allocate working colors : red, green, yellow, blue, white.
  ;
  TVLCT, [255,0,255,0,0],[0,255,255,0,255],[0,0,0,255,255], topClr+1

  ;  Get the windows (drawing areas) identifiers.
  ;
  for i=0,2 do begin
    WIDGET_CONTROL, xd_sState.draw[i], GET_VALUE=j
    xd_sState.window[i] = j
  endfor
;for dilation
radius = 2
strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius


xyImage = DBLARR(imageSize, imageSize)
temp_xyImage = svol_HU_cube[*, *, location3D[2]]
xyImage[*, *] = temp_xyImage[*, *, 0]
WSET, xd_sState.window[0]
xyMask = CONGRID(mask_vol_cube[*, *, location3D[2]*(SIZE(vol_HU_cube))[3]/imageSize], imageSize, imageSize, 1)
xyMask = DILATE(xyMask, strucElem)
;xyMask = BYTSCL(xyImage, TOP = 512)  * (xyImage GE threshold_min_value) * (xyImage LE threshold_max_value)
TV, REVERSE(BYTSCL(xyImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value), 2)
CONTOUR, REVERSE(xyMask,2), LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

;delete the variable
temp_xyImage = !NULL
xyImage = !NULL
xyMask = !NULL

FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
  line = [LINDGEN(PathInfo(I).N), 0]
  oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
  DRAW_ROI, oROI, COLOR =topClr+1, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
  OBJ_DESTROY, oROI
ENDFOR
plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
plots, [0,imageSize], [imageSize - location3D[1],imageSize - location3D[1]], COLOR=topClr+3, /DEVICE

yzImage = DBLARR(imageSize, imageSize)
temp_yzImage = svol_HU_cube[location3D[0], *, *]
yzImage[*, *] = temp_yzImage[0, *, *]
WSET, xd_sState.window[1]
yzMask = CONGRID(mask_vol_cube[location3D[0]*(SIZE(vol_HU_cube))[1]/imageSize, *, *], 1, imageSize, imageSize)
yzMask = DILATE(REFORM(yzMask), strucElem)
;yzMask = BYTSCL(yzImage, TOP = 512) * (yzImage GE threshold_min_value) * (yzImage LE threshold_max_value)
TV, BYTSCL(yzImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value)
CONTOUR, yzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

;delete the variable
temp_yzImage = !NULL
yzImage = !NULL
yzMask = !NULL

FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
  line = [LINDGEN(PathInfo(I).N), 0]
  oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
  DRAW_ROI, oROI, COLOR =topClr+2, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
  OBJ_DESTROY, oROI
ENDFOR
plots, [location3D[1],location3D[1]], [0,imageSize], COLOR=topClr+3, /DEVICE
plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

xzImage = DBLARR(imageSize, imageSize)
temp_xzImage = svol_HU_cube[*, location3D[1], *]
xzImage[*, *] = temp_xzImage[*, 0, *]
WSET, xd_sState.window[2]
xzMask = DBLARR(imageSize, imageSize)
temp_xzMask = CONGRID(mask_vol_cube[*, location3D[1]*(SIZE(vol_HU_cube))[2]/imageSize, *], imageSize, 1, imageSize)
xzMask[*, *] = temp_xzMask[*, 0, *]
xzMask = DILATE(REFORM(xzMask), strucElem)
;xzMask = BYTSCL(xzImage, TOP = 512) * (xzImage GE threshold_min_value) * (xzImage LE threshold_max_value)
TV, BYTSCL(xzImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value)
CONTOUR, xzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

;delete the variable
temp_xzImage = !NULL
xzImage = !NULL
xzMask = !NULL

FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
  line = [LINDGEN(PathInfo(I).N), 0]
  oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
  DRAW_ROI, oROI, COLOR =topClr+3, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
  OBJ_DESTROY, oROI
ENDFOR
plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

;delete the variables
vol_HU_cube = !NULL
svol_HU_cube = !NULL

;LOADCT, 12
;IF OBJ_VALID(xd_sState.xyROIout) THEN BEGIN
;  xd_sState.xyROIout->getProperty, DATA = xyROIverts
;  IF N_ELEMENTS(xyROIverts) GE 3 THEN BEGIN
;    ;overplot ROI region
;    WSET, xd_sState.window[0]
;    xd_sState.xyROIout->setProperty, DATA = xyROIverts/3
;    DRAW_ROI, xd_sState.xyROIout, /LINE_FILL, COLOR = 80, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;    xd_sState.xyROIout->setProperty, DATA = xyROIverts
;  ENDIF
;ENDIF
;
;IF OBJ_VALID(xd_sState.yzROIout) THEN BEGIN
;  xd_sState.yzROIout->getProperty, DATA = yzROIverts
;  IF N_ELEMENTS(yzROIverts) GE 3 THEN BEGIN
;    ;overplot ROI region
;    WSET, xd_sState.window[1]
;    xd_sState.yzROIout->setProperty, DATA = yzROIverts/3
;    DRAW_ROI, xd_sState.yzROIout, /LINE_FILL, COLOR = 42, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;    xd_sState.yzROIout->setProperty, DATA = yzROIverts
;  ENDIF
;ENDIF
;
;IF OBJ_VALID(xd_sState.xzROIout) THEN BEGIN
;  xd_sState.xzROIout->getProperty, DATA = xzROIverts
;  IF N_ELEMENTS(xzROIverts) GE 3 THEN BEGIN
;    ;overplot ROI region
;    WSET, xd_sState.window[2]
;    xd_sState.xzROIout->setProperty, DATA = xzROIverts/3
;    DRAW_ROI, xd_sState.xzROIout, /LINE_FILL, COLOR = 120, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;    xd_sState.xzROIout->setProperty, DATA = xzROIverts
;  ENDIF
;ENDIF

;print, 'ortho2D_plot done!'

IF KEYWORD_SET(xd_sState.additionalDir) THEN BEGIN
  IF N_ELEMENTS(xd_sState.additionalDir) GE 2 THEN BEGIN
      vol_HU_cube_a1 = xd_sState.vol_HU_cube_a1
      svol_HU_cube_a1 = xd_sState.svol_HU_cube_a1
      
        ;calculate the additional image's threshold
        sizeVHc = SIZE(svol_HU_cube_a1)
        mean_image = IMAGE_THRESHOLD(BYTSCL(svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
        image_statistics, svol_HU_cube_a1[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
        threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
        threshold_min_value = img_min
              
      ; 闁活枎鎷风紓浣侯棞琚欓柛鎾寸墪濞达拷

      topClr = !D.TABLE_SIZE - 6
      ;  Load the grey scale colr table.
      ;
      LOADCT, 0, /SILENT, NCOLORS=topClr+1

      ;  Allocate working colors : red, green, yellow, blue, white.
      ;
      TVLCT, [255,0,255,0,0],[0,255,255,0,255],[0,0,0,255,255], topClr+1
      ;  Get the windows (drawing areas) identifiers.
      ;
      for i=0,2 do begin
        WIDGET_CONTROL, xd_sState.draw_a1[i], GET_VALUE=j
        xd_sState.window[i] = j
      endfor

      xyImage = DBLARR(imageSize, imageSize)
      temp_xyImage = svol_HU_cube_a1[*, *, location3D[2]]
      xyImage[*, *] = temp_xyImage[*, *, 0]
      WSET, xd_sState.window[0]
      xyMask = CONGRID(mask_vol_cube[*, *, location3D[2]*(SIZE(vol_HU_cube_a1))[3]/imageSize], imageSize, imageSize, 1)
      xyMask = DILATE(xyMask, strucElem)
      ;xyMask = BYTSCL(xyImage, TOP = 512)  * (xyImage GE threshold_min_value) * (xyImage LE threshold_max_value)
      TV, REVERSE(BYTSCL(xyImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value), 2)
      CONTOUR, REVERSE(xyMask,2), LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

      ;delete the variable
      temp_xyImage = !NULL
      xyImage = !NULL
      xyMask = !NULL

      FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
        line = [LINDGEN(PathInfo(I).N), 0]
        oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
        DRAW_ROI, oROI, COLOR =topClr+1, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
        OBJ_DESTROY, oROI
      ENDFOR
      plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
      plots, [0,imageSize], [imageSize - location3D[1],imageSize - location3D[1]], COLOR=topClr+3, /DEVICE

      yzImage = DBLARR(imageSize, imageSize)
      temp_yzImage = svol_HU_cube_a1[location3D[0], *, *]
      yzImage[*, *] = temp_yzImage[0, *, *]
      
      WSET, xd_sState.window[1]
      yzMask = CONGRID(mask_vol_cube[location3D[0]*(SIZE(vol_HU_cube_a1))[1]/imageSize, *, *], 1, imageSize, imageSize)
      yzMask = DILATE(REFORM(yzMask), strucElem)
      ;yzMask = BYTSCL(yzImage, TOP = 512) * (yzImage GE threshold_min_value) * (yzImage LE threshold_max_value)
      TV, BYTSCL(yzImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value)
      CONTOUR, yzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

      ;delete the variable
      temp_yzImage = !NULL
      yzImage = !NULL
      yzMask = !NULL

      FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
        line = [LINDGEN(PathInfo(I).N), 0]
        oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
        DRAW_ROI, oROI, COLOR =topClr+2, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
        OBJ_DESTROY, oROI
      ENDFOR
      plots, [location3D[1],location3D[1]], [0,imageSize], COLOR=topClr+3, /DEVICE
      plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

      xzImage = DBLARR(imageSize, imageSize)
      temp_xzImage = svol_HU_cube_a1[*, location3D[1], *]
      xzImage[*, *] = temp_xzImage[*, 0, *]
      WSET, xd_sState.window[2]
      xzMask = DBLARR(imageSize, imageSize)
      temp_xzMask = CONGRID(mask_vol_cube[*, location3D[1]*(SIZE(vol_HU_cube_a1))[2]/imageSize, *], imageSize, 1, imageSize)
      xzMask[*, *] = temp_xzMask[*, 0, *]
      xzMask = DILATE(REFORM(xzMask), strucElem)
      ;xzMask = BYTSCL(xzImage, TOP = 512) * (xzImage GE threshold_min_value) * (xzImage LE threshold_max_value)
      TV, BYTSCL(xzImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value)
      CONTOUR, xzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

      ;delete the variable
      temp_xzImage = !NULL
      xzImage = !NULL
      xzMask = !NULL

      FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
        line = [LINDGEN(PathInfo(I).N), 0]
        oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
        DRAW_ROI, oROI, COLOR =topClr+3, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
        OBJ_DESTROY, oROI
      ENDFOR
      plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
      plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

      ;delete the variables
      vol_HU_cube_a1 = !NULL
      svol_HU_cube_a1 = !NULL
      
      

;      LOADCT, 12
;      IF OBJ_VALID(xd_sState.xyROIout) THEN BEGIN
;        xd_sState.xyROIout->getProperty, DATA = xyROIverts
;        IF N_ELEMENTS(xyROIverts) GE 3 THEN BEGIN
;          ;overplot ROI region
;          WSET, xd_sState.window[0]
;          xd_sState.xyROIout->setProperty, DATA = xyROIverts/3
;          DRAW_ROI, xd_sState.xyROIout, /LINE_FILL, COLOR = 80, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;          xd_sState.xyROIout->setProperty, DATA = xyROIverts
;        ENDIF
;      ENDIF
;
;      IF OBJ_VALID(xd_sState.yzROIout) THEN BEGIN
;        xd_sState.yzROIout->getProperty, DATA = yzROIverts
;        IF N_ELEMENTS(yzROIverts) GE 3 THEN BEGIN
;          ;overplot ROI region
;          WSET, xd_sState.window[1]
;          xd_sState.yzROIout->setProperty, DATA = yzROIverts/3
;          DRAW_ROI, xd_sState.yzROIout, /LINE_FILL, COLOR = 42, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;          xd_sState.yzROIout->setProperty, DATA = yzROIverts
;        ENDIF
;      ENDIF
;
;      IF OBJ_VALID(xd_sState.xzROIout) THEN BEGIN
;        xd_sState.xzROIout->getProperty, DATA = xzROIverts
;        IF N_ELEMENTS(xzROIverts) GE 3 THEN BEGIN
;          ;overplot ROI region
;          WSET, xd_sState.window[2]
;          xd_sState.xzROIout->setProperty, DATA = xzROIverts/3
;          DRAW_ROI, xd_sState.xzROIout, /LINE_FILL, COLOR = 120, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;          xd_sState.xzROIout->setProperty, DATA = xzROIverts
;        ENDIF
;      ENDIF
      
      IF N_ELEMENTS(xd_sState.additionalDir) GE 3 THEN BEGIN
        vol_HU_cube_a2 = xd_sState.vol_HU_cube_a2
        svol_HU_cube_a2 = xd_sState.svol_HU_cube_a2
  
        ;calculate the additional image's threshold
        sizeVHc = SIZE(svol_HU_cube_a2)
        mean_image = IMAGE_THRESHOLD(BYTSCL(svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2]), THRESHOLD=tf, /MAXENTROPY)
        image_statistics, svol_HU_cube_a2[sizeVHc[1]/4:sizeVHc[1]*3/4,sizeVHc[2]/4:sizeVHc[2]*3/4,sizeVHc[3]/2], MAX=img_max, MIN=img_min
        threshold_max_value = (tf[0]*(img_max-img_min)/256.+img_min);*0.5;GJ, 2019/5/19;-100.
        threshold_min_value = img_min
        
        ; 闁活枎鎷风紓浣侯棞琚欓柛鎾寸墪濞达拷

        topClr = !D.TABLE_SIZE - 6
        ;  Load the grey scale colr table.
        ;
        LOADCT, 0, /SILENT, NCOLORS=topClr+1

        ;  Allocate working colors : red, green, yellow, blue, white.
        ;
        TVLCT, [255,0,255,0,0],[0,255,255,0,255],[0,0,0,255,255], topClr+1
        ;  Get the windows (drawing areas) identifiers.
        ;
        for i=0,2 do begin
          WIDGET_CONTROL, xd_sState.draw_a2[i], GET_VALUE=j
          xd_sState.window[i] = j
        endfor



        xyImage = DBLARR(imageSize, imageSize)
        temp_xyImage = svol_HU_cube_a2[*, *, location3D[2]]
        xyImage[*, *] = temp_xyImage[*, *, 0]
        WSET, xd_sState.window[0]
        xyMask = CONGRID(mask_vol_cube[*, *, location3D[2]*(SIZE(vol_HU_cube_a2))[3]/imageSize], imageSize, imageSize, 1)
        xyMask = DILATE(xyMask, strucElem)
        ;xyMask = BYTSCL(xyImage, TOP = 512)  * (xyImage GE threshold_min_value) * (xyImage LE threshold_max_value)
        TV, REVERSE(BYTSCL(xyImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value), 2)
        CONTOUR, REVERSE(xyMask,2), LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

        ;delete the variable
        xyImage2 = xyImage
        temp_xyImage = !NULL
        xyImage = !NULL
        xyMask = !NULL

        FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
          line = [LINDGEN(PathInfo(I).N), 0]
          oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
          DRAW_ROI, oROI, COLOR =topClr+1, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
          OBJ_DESTROY, oROI
        ENDFOR
        plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
        plots, [0,imageSize], [imageSize - location3D[1],imageSize - location3D[1]], COLOR=topClr+3, /DEVICE

        yzImage = DBLARR(imageSize, imageSize)
        temp_yzImage = svol_HU_cube_a2[location3D[0], *, *]
        yzImage[*, *] = temp_yzImage[0, *, *]
        WSET, xd_sState.window[1]
        yzMask = CONGRID(mask_vol_cube[location3D[0]*(SIZE(vol_HU_cube_a2))[1]/imageSize, *, *], 1, imageSize, imageSize)
        yzMask = DILATE(REFORM(yzMask), strucElem)
        ;yzMask = BYTSCL(yzImage, TOP = 512) * (yzImage GE threshold_min_value) * (yzImage LE threshold_max_value)
        TV, BYTSCL(yzImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value)
        CONTOUR, yzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

        ;delete the variable
        yzImage2 = yzImage
        temp_yzImage = !NULL
        yzImage = !NULL
        yzMask = !NULL

        FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
          line = [LINDGEN(PathInfo(I).N), 0]
          oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
          DRAW_ROI, oROI, COLOR =topClr+2, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
          OBJ_DESTROY, oROI
        ENDFOR
        plots, [location3D[1],location3D[1]], [0,imageSize], COLOR=topClr+3, /DEVICE
        plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

        xzImage = DBLARR(imageSize, imageSize)
        temp_xzImage = svol_HU_cube_a2[*, location3D[1], *]
        xzImage[*, *] = temp_xzImage[*, 0, *]
        WSET, xd_sState.window[2]
        xzMask = DBLARR(imageSize, imageSize)
        temp_xzMask = CONGRID(mask_vol_cube[*, location3D[1]*(SIZE(vol_HU_cube_a2))[2]/imageSize, *], imageSize, 1, imageSize)
        xzMask[*, *] = temp_xzMask[*, 0, *]
        xzMask = DILATE(REFORM(xzMask), strucElem)
        ;xzMask = BYTSCL(xzImage, TOP = 512) * (xzImage GE threshold_min_value) * (xzImage LE threshold_max_value)
        TV, BYTSCL(xzImage, TOP = topClr, MAX=threshold_max_value, MIN=threshold_min_value)
        CONTOUR, xzMask, LEVEL = 1, XMARGIN = [0, 0], YMARGIN = [0, 0], /NOERASE, PATH_INFO = pathInfo, PATH_XY = pathXY, XSTYLE = 5, YSTYLE = 5, /PATH_DATA_COORDS

        ;delete the variable
        xzImage2 = xzImage
        temp_xzImage = !NULL
        xzImage = !NULL
        xzMask = !NULL

        FOR I = 0,(N_ELEMENTS(PathInfo) - 1 ) DO BEGIN
          line = [LINDGEN(PathInfo(I).N), 0]
          oROI = OBJ_NEW('IDLanROI', (pathXY(*, pathInfo(I).OFFSET + line))[0, *], (pathXY(*, pathInfo(I).OFFSET + line))[1, *])
          DRAW_ROI, oROI, COLOR =topClr+3, /LINE_FILL, SPACING=0.005, ORIENTATION=45, /DEVICE
          OBJ_DESTROY, oROI
        ENDFOR
        plots, [location3D[0],location3D[0]], [0,imageSize], COLOR=topClr+2, /DEVICE
        plots, [0,imageSize], [location3D[2],location3D[2]], COLOR=topClr+4, /DEVICE

        ;delete the variables
        vol_HU_cube_a2 = !NULL
        svol_HU_cube_a2 = !NULL

;        LOADCT, 12
;        IF OBJ_VALID(xd_sState.xyROIout) THEN BEGIN
;          xd_sState.xyROIout->getProperty, DATA = xyROIverts
;          IF N_ELEMENTS(xyROIverts) GE 3 THEN BEGIN
;            ;overplot ROI region
;            WSET, xd_sState.window[0]
;            xd_sState.xyROIout->setProperty, DATA = xyROIverts/3
;            DRAW_ROI, xd_sState.xyROIout, /LINE_FILL, COLOR = 80, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;            xd_sState.xyROIout->setProperty, DATA = xyROIverts
;          ENDIF
;        ENDIF
;
;        IF OBJ_VALID(xd_sState.yzROIout) THEN BEGIN
;          xd_sState.yzROIout->getProperty, DATA = yzROIverts
;          IF N_ELEMENTS(yzROIverts) GE 3 THEN BEGIN
;            ;overplot ROI region
;            WSET, xd_sState.window[1]
;            xd_sState.yzROIout->setProperty, DATA = yzROIverts/3
;            DRAW_ROI, xd_sState.yzROIout, /LINE_FILL, COLOR = 42, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;            xd_sState.yzROIout->setProperty, DATA = yzROIverts
;          ENDIF
;        ENDIF
;
;        IF OBJ_VALID(xd_sState.xzROIout) THEN BEGIN
;          xd_sState.xzROIout->getProperty, DATA = xzROIverts
;          IF N_ELEMENTS(xzROIverts) GE 3 THEN BEGIN
;            ;overplot ROI region
;            WSET, xd_sState.window[2]
;            xd_sState.xzROIout->setProperty, DATA = xzROIverts/3
;            DRAW_ROI, xd_sState.xzROIout, /LINE_FILL, COLOR = 120, SPACING = 0.1, ORIENTATION = 315, /DEVICE
;            xd_sState.xzROIout->setProperty, DATA = xzROIverts
;          ENDIF
;        ENDIF
      ENDIF
    ENDIF
  ENDIF

mask_vol_cube = !NULL

end

;@GJ, 2023/5/16, generating the fusion map of MPI and CT
PRO image_fusion_MPI_CT, xd_sState
  
  location3D = xd_sState.location3D
  ;GJ 2019/6/1, negative location3D values may give error, bug was fixed.
  IF location3D[0] LT 0 THEN location3D[0] = 0
  IF location3D[1] LT 0 THEN location3D[1] = 0
  IF location3D[2] LT 0 THEN location3D[2] = 0
  ;print,'2D_plot location3D',location3D
  imageSize = xd_sState.imageSize
  
  vol_HU_cube = xd_sState.vol_HU_cube
  vol_HU_cube_a1 = xd_sState.vol_HU_cube_a1
  vol_size = SIZE(vol_HU_cube, /DIMENSIONS)
  vol_a1_size = SIZE(vol_HU_cube_a1, /DIMENSIONS)
  
  FOR i=0, 5 DO BEGIN
    WIDGET_CONTROL, xd_sState.Draw_3DView[i], GET_VALUE=j
    xd_sState.window[i] = j
  ENDFOR
  
  ;calculate the cross-sectional image
  xyImage0 = REFORM(vol_HU_cube[*, *, FIX(location3D[2]*1.0/imageSize*vol_size[2])])
  xyImage1 = REFORM(vol_HU_cube_a1[*, *, FIX(location3D[2]*1.0/imageSize*vol_a1_size[2])])
  yzImage0 = REFORM(vol_HU_cube[FIX(location3D[0]*1.0/imageSize*vol_size[0]), *, *])
  yzImage1 = REFORM(vol_HU_cube_a1[FIX(location3D[0]*1.0/imageSize*vol_a1_size[0]), *, *])
  xzImage0 = REFORM(vol_HU_cube[*, FIX(location3D[1]*1.0/imageSize*vol_size[1]), *])
  xzImage1 = REFORM(vol_HU_cube_a1[*, FIX(location3D[1]*1.0/imageSize*vol_a1_size[1]), *])

  ;@GJ, 2023/5/18, calculating the uniform scaling factor
  Max_xy = MAX(xyImage1, MIN=Min_xy)
  th_result = IMAGE_THRESHOLD(xyImage1, threshold=th_xy, /MAXENTROPY)
  Max_yz = MAX(yzImage1, MIN=Min_yz)
  th_result = IMAGE_THRESHOLD(yzImage1, threshold=th_yz, /MAXENTROPY)
  Max_xz = MAX(xzImage1, MIN=Min_xz)
  th_result = IMAGE_THRESHOLD(xzImage1, threshold=th_xz, /MAXENTROPY)
  Max_overall = MEAN([Max_xy, Max_yz, Max_xz])
  Min_overall = MEAN([Min_xy, Min_yz, Min_xz])
  
  ;@GJ, 2023/5/18, presetting the color table
  max_ct=255;200;max_ct=2000
  min_ct=MIN([th_xy, th_yz, th_xz]);150
  major=FIX((max_ct-min_ct)/10.)
  
  ;fusion xy  
  foregroundImage = BYTSCL(xyImage1, MAX=Max_overall, MIN=Min_overall)
  foregroundImage = BYTSCL(foregroundImage, MAX=max_ct, MIN=min_ct)
;  WSET, xd_sState.window[0]
;  cgImage, xyImage0, CTIndex=0
;  cgImage, foregroundImage, CTIndex=3, Transparent=50, Missing_Value=0.5
  

  closeYes=1
  dir_temp = xd_sState.additionaldir[0]
  pos_dir = STRPOS(STRMID(dir_temp, 0, STRLEN(dir_temp)-1), '\', /REVERSE_SEARCH)
  systm = STRING(SYSTIME(/JULIAN, /UTC), FORMAT='(i12)')
  cur_dir = STRMID(dir_temp, 0, pos_dir+1) + STRTRIM(systm,1) + '_'
  
  filename0=cur_dir+'fusion_c3_xy.tiff';'unet.jpg'
  Image_Blend, BYTSCL(xyImage0), foregroundImage, COLORTABLE=3, blendTitle='MPI-CT', major, max_ct, min_ct, filename0, closeYes, 'MPI-CT'
  image = READ_TIFF(filename0, R, G, B)
  WSET, xd_sState.window[0]
  cgImage, REVERSE(image, 3)
  
  filename3=cur_dir+'fusion_c4_xy.tiff';'unet.jpg'
  Image_Blend, BYTSCL(xyImage0), foregroundImage, COLORTABLE=4, blendTitle='MPI-CT', major, max_ct, min_ct, filename3, closeYes, 'MPI-CT'
  image = READ_TIFF(filename3, R, G, B)
  WSET, xd_sState.window[3]
  cgImage, REVERSE(image, 3)
    
  ;fusion yz
  foregroundImage = BYTSCL(yzImage1, MAX=Max_overall, MIN=Min_overall)
  foregroundImage = BYTSCL(foregroundImage, MAX=max_ct, MIN=min_ct)
;  iimage, foregroundImage
;  foregroundImage([WHERE(foregroundImage LT min_ct, /NULL)]) = 0
  filename1=cur_dir+'fusion_c3_yz.tiff';'unet.jpg'
  Image_Blend, BYTSCL(yzImage0), foregroundImage, COLORTABLE=3, blendTitle='MPI-CT', major, max_ct, min_ct, filename1, closeYes, 'MPI-CT'
  image = READ_TIFF(filename1, R, G, B)
  WSET, xd_sState.window[1]
  cgImage, REVERSE(image, 3)
  
  filename4=cur_dir+'fusion_c4_yz.tiff';'unet.jpg'
  Image_Blend, BYTSCL(yzImage0), foregroundImage, COLORTABLE=4, blendTitle='MPI-CT', major, max_ct, min_ct, filename4, closeYes, 'MPI-CT'
  image = READ_TIFF(filename4, R, G, B)
  WSET, xd_sState.window[4]
  cgImage, REVERSE(image, 3)

  ;fusion xz
  foregroundImage = BYTSCL(xzImage1, MAX=Max_overall, MIN=Min_overall)
  foregroundImage = BYTSCL(foregroundImage, MAX=max_ct, MIN=min_ct)
  filename2=cur_dir+'fusion_c3_xz.tiff';'unet.jpg'
  Image_Blend, BYTSCL(xzImage0), foregroundImage, COLORTABLE=3, blendTitle='MPI-CT', major, max_ct, min_ct, filename2, closeYes, 'MPI-CT'
  image = READ_TIFF(filename2, R, G, B)
  WSET, xd_sState.window[2]
  cgImage, REVERSE(image, 3)
    
  filename5=cur_dir+'fusion_c4_xz.tiff';'unet.jpg'
  Image_Blend, BYTSCL(xzImage0), foregroundImage, COLORTABLE=4, blendTitle='MPI-CT', major, max_ct, min_ct, filename5, closeYes, 'MPI-CT'
  image = READ_TIFF(filename5, R, G, B)
  WSET, xd_sState.window[5]
  cgImage, REVERSE(image, 3)

END

PRO Generate_Mask_ROI_modify, xd_sState, old_mask_vol_cube, ROIout, modify_roi, ContOrInde, screenSize, direction = direction

n_ROIs = N_ELEMENTS(ROIout)
Loc_array = INTARR(n_ROIs)
Area_array = DBLARR(n_ROIs)

Image = DBLARR(screenSize, screenSize)
dims = SIZE(Image, /DIMENSIONS)
vol_dims = SIZE(xd_sState.vol_HU_cube, /DIMENSIONS)
mask_ROI = BYTARR(vol_dims[0], vol_dims[1], N_ELEMENTS(ROIout))

;initialize the modifiy mask
xd_sState.mask_vol_cube_modify = xd_sState.mask_vol_cube_modify * 0

FOR i=0, n_ROIs-1 DO BEGIN
  ROIout[i] -> getProperty, ALPHA_CHANNEL=ac
  print, 'ac = ', ac
  result = ROIout[i]->ComputeGeometry(AREA=area, CENTROID=cent)
  Area_array[i] = area
  maskResult = ROIout[i] -> ComputeMask(DIMENSIONS = [screenSize, screenSize])
  IF STRCMP(direction, 'xy') THEN BEGIN
    Loc_array[i] = FLOOR(FLOOR(ac * 10000.)* (vol_dims[2]) / xd_sState.imageSize)
    xd_sState.mask_vol_cube_modify[*,*, Loc_array[i]] = xd_sState.mask_vol_cube_modify[*,*, Loc_array[i]] OR ((REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0)
  ENDIF
  IF STRCMP(direction, 'yz') THEN BEGIN
    Loc_array[i] = FLOOR(FLOOR(ac * 10000.)* (vol_dims[0]) / xd_sState.imageSize)
    xd_sState.mask_vol_cube_modify[Loc_array[i],*,*] = xd_sState.mask_vol_cube_modify[Loc_array[i],*,*] OR (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
  ENDIF
  IF STRCMP(direction, 'xz') THEN BEGIN
    Loc_array[i] = FLOOR(FLOOR(ac * 10000.)* (vol_dims[1]) / xd_sState.imageSize)
    xd_sState.mask_vol_cube_modify[*,Loc_array[i],*] = xd_sState.mask_vol_cube_modify[*,Loc_array[i],*] OR (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)
  ENDIF
ENDFOR

  print, 'loc_array = ', loc_array
;multiple ROIs
IF n_ROIs GT 1 THEN BEGIN
  sort_ind = SORT(Loc_array)
  Loc_array = Loc_array[sort_ind]
  ROIout = ROIout[sort_ind]
  Area_array = Area_array[sort_ind]

  FOR j=0, n_ROIs-2 DO BEGIN
    FOR k=j+1, n_ROIs-1 DO BEGIN
      Loc_diff = Loc_array[k] - Loc_array[j]
      IF (Loc_diff GT 0) AND (ContOrInde EQ 1) THEN BEGIN
        ROI_a = ROIout[j]
        ROI_b = ROIout[k]
        loc_a = Loc_array[j]
        loc_b = Loc_array[k]
        temp=ROI_pair_to_mask(xd_sState.vol_HU_cube, xd_sState.mask_vol_cube_modify, ROI_a, ROI_b, loc_a, Loc_b, direction, screenSize)
        xd_sState.mask_vol_cube_modify = temp
        k = n_ROIs
      ENDIF
    ENDFOR
  ENDFOR
ENDIF

;iimage, xd_sState.mask_vol_cube_modify[*,*,128]
print, 'modify_roi = ', modify_roi
;make new polygon, add new polygon
IF ABS(modify_roi-12) LE 1 THEN BEGIN
  threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  IF threshold_min_value LT 0 THEN threshold_min_value = 0
  mask_vol_cube = xd_sState.mask_vol_cube_modify * (threshold_min_value+1+1000) -1000
ENDIF

;for add and remove
IF ABS(modify_roi) LT 10 THEN BEGIN
  ;segmentation
  threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  threshold_max_value = MAX([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  print, 'threshold_min_value = ', threshold_min_value
  print, 'threshold_max_value = ', threshold_max_value
  thresh_ave = 0.5 * (threshold_min_value + threshold_max_value)
  ;mask the 3D volume
  mask_vol_cube = DOUBLE(xd_sState.mask_vol_cube) * 0. ;zero the mask

  FOR imind = 0, N_ELEMENTS(xd_sState.vol_HU_cube[0, 0, *])-1 DO BEGIN
    ;add
    IF ABS(modify_roi-4) LE 1 THEN BEGIN
      test_mask = (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0) OR (old_mask_vol_cube[*, *, imind] GT 0)
      mask_vol_cube[*, *, imind] = test_mask*thresh_ave

      ;delete the variable
      test_mask = !NULL
    ENDIF

    ;remove selected region and 3d segmentation 3dseg
    IF ABS(modify_roi-8) LE 1 THEN BEGIN
      test_mask = (1-(xd_sState.mask_vol_cube_modify[*, *, imind] GT 0)) AND (old_mask_vol_cube[*, *, imind] GT 0)
      mask_vol_cube[*, *, imind] = test_mask*thresh_ave

      ;delete the variable
      test_mask = !NULL
    ENDIF
  ENDFOR
ENDIF

;iimage, mask_vol_cube[*,*,128]

SHADE_VOLUME, mask_vol_cube,threshold_min_value, Outverts, Outconn1
  
IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
  smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
  verts = smoothedOutverts
  xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
  xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
  xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
  verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
    [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
    [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
    [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
    [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
    [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
    [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
    [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
  xd_sState.oPoly_border->setProperty, DATA = verts_border

  ;define to the rotation center
  ma=max(smoothedOutverts)
  mi=min(smoothedOutverts)
  bias=(ma-mi)/2
  rl=mi+bias
  rr=ma+bias
  cc=[(-rl)/(rr-rl),1/(rr-rl)]
  xma = max((smoothedOutverts)[0,*])
  xmi = min((smoothedOutverts)[0,*])
  xmid = 0.5*(xma+xmi)*cc[1]
  yma = max((smoothedOutverts)[1,*])
  ymi = min((smoothedOutverts)[1,*])
  ymid = 0.5*(yma+ymi)*cc[1]
  zma = max((smoothedOutverts)[2,*])
  zmi = min((smoothedOutverts)[2,*])
  zmid = 0.5*(zma+zmi)*cc[1]
  
  ;hide old object
  WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=old_ObjectIndex_value
  
  ;Update modify_Step
  modify_Step = xd_sState.struc_modifyobject.current_step_index +1
  xd_sState.struc_modifyobject.modify_step = modify_Step
  xd_sState.struc_modifyobject.current_step_index = modify_Step
  widget_control, xd_sState.wUndo, SENSITIVE = 1

  ;update exist number
  exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
  xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
  n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

  ObjectIndex_value = LONG(n_Poly)
  print, 'get_new_object_value = ', ObjectIndex_value
  WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
  WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
  WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

  xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
  xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0
  ;If modifying current object, hide old one after modification
  ;If adding new object, we don't need to hid old ones. GJ, 2019/4/19, gjia@xidian.edu.cn
  IF (old_ObjectIndex_value GE 1) AND (ABS(modify_roi) LT 10) THEN BEGIN
    print, 'hiding object =', old_ObjectIndex_value-1
    xd_sState.struc_modifyobject.hide_status[modify_Step, old_ObjectIndex_value-1] = 1
    xd_sState.oPoly[old_ObjectIndex_value-1]->SetProperty, Hide = 1
  ENDIF
  
  ;define a new polygon with color
  color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
  vert_colors = smoothedOutverts
  vert_colors[0,*] = color_rgb[n_Poly, 0]
  vert_colors[1,*] = color_rgb[n_Poly, 1]
  vert_colors[2,*] = color_rgb[n_Poly, 2]

  ;define mask_vol_cube
  index = WHERE(mask_vol_cube GT threshold_min_value, count)
  IF count GT 1 THEN BEGIN
    xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
    xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

    xd_sState.tvol_white.SetProperty,hide=0
    xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
    Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
    ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
    IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
  ENDIF

  alpha_channel = 1.0
  oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
    SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
  xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
  xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
  xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
  xd_sState.selected = xd_sState.oPoly[n_Poly-1]
  xd_sState.oRotationModel->Add, oPoly
  FOR i = 0, n_Poly-1 DO BEGIN
    xd_sState.opoly[i]->getProperty, data=verts
    IF N_ELEMENTS(verts) EQ 0 THEN RETURN
    IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
  ENDFOR
  ;center the selected object
  Adjust_whole_polygons, whole_verts, xd_sState



  dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
  ;GJ 2019/2/22, summarize all plot pickdata into a program
  PLOT_pickdata, dataxyz, xd_sState


  ;set the alpha-channel transparent value
  WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
  WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

  ;Enable color selection
  WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
  WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
  WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
  DEVICE, DECOMPOSED=1
  WSET, wColorIndexID
  ERASE, ColorIndexBGR[n_Poly]

  ;enable these options
  ;    widget_control, xd_sState.wColorIndex, sensitive=1
  widget_control, xd_sState.wSmooth, sensitive=1
  widget_control, xd_sState.wSmooth_text, sensitive=1
  widget_control, xd_sState.wtransSlider, sensitive=1
  widget_control, xd_sState.wTrans_text, sensitive=1
  widget_control, xd_sState.w3DSegd, sensitive=1
  widget_control, xd_sState.w3Drg, sensitive=1
  widget_control, xd_sState.wthickness, sensitive=1
  widget_control, xd_sState.wpoint, sensitive=1
  widget_control, xd_sState.wwire, sensitive=1
  widget_control, xd_sState.wfill, sensitive=1
  widget_control, xd_sState.wexpstl, sensitive=1
  ;            widget_control, xd_sState.whide, sensitive=0
  widget_control, xd_sState.wsinglecenter, sensitive=1

  widget_control, xd_sState.wUndo, SENSITIVE = 1
  ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

ENDIF ELSE BEGIN
  infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
ENDELSE

END

FUNCTION ROI_pair_to_mask, vol_HU_cube, mask_vol_cube_modify, ROI_a, ROI_b, loc_a, Loc_b, direction, screenSize

ROI_a->GetProperty, DATA=roiData_a, N_Verts = nverts_a, ROI_XRANGE=xr_a, ROI_YRANGE=yr_a
xrange_a = xr_a[1] - xr_a[0]
yrange_a = yr_a[1] - yr_a[0]
ROI_b->GetProperty, DATA=roiData_b, N_Verts = nverts_b, ROI_XRANGE=xr_b, ROI_YRANGE=yr_b
xrange_b = xr_b[1] - xr_b[0]
yrange_b = yr_b[1] - yr_b[0]
resultG_a = ROI_a->ComputeGeometry(AREA=area_a, CENTROID=cent_a, SPATIAL_OFFSET = spao_a, SPATIAL_SCALE= spac_a)
resultG_b = ROI_b->ComputeGeometry(AREA=area_b, CENTROID=cent_b, SPATIAL_OFFSET = spao_b, SPATIAL_SCALE= spac_b)

nframes = DOUBLE(ABS(loc_b - loc_a))
cent_dif = (cent_b - cent_a) / nframes
;a to b
for i=1, nframes-1 do begin
  temp_cent = cent_a + cent_dif * DOUBLE(i)
  ROI_a->Translate, temp_cent[0] - cent_a[0], temp_cent[1] - cent_a[1]
  ROI_b->Translate, temp_cent[0] - cent_b[0], temp_cent[1] - cent_b[1]
  
  t_a = DOUBLE(i) / (nframes)
  Sx_a = 1. + t_a * (xrange_b - xrange_a)/xrange_a
  Sy_a = 1. + t_a * (yrange_b - yrange_a)/yrange_a
  resultGt_a = ROI_a->ComputeGeometry(AREA=tarea_a, CENTROID=tcent_a)
  ROI_a->Translate, -tcent_a
  ROI_a->Scale, Sx_a, Sy_a
  ROI_a->Translate, tcent_a

  t_b = DOUBLE(nframes-i) / (nframes)
  Sx_b = 1. + t_b * (xrange_a - xrange_b)/xrange_b
  Sy_b = 1. + t_b * (yrange_a - yrange_b)/yrange_b
  resultGt_b = ROI_b->ComputeGeometry(AREA=tarea_b, CENTROID=tcent_b)
  ROI_b->Translate, -tcent_b
  ROI_b->Scale, Sx_b, Sy_b
  ROI_b->Translate, tcent_b
  
  vol_dims = SIZE(vol_HU_cube, /DIMENSIONS)
  maskResult = ROI_a -> ComputeMask(DIMENSIONS = [screenSize, screenSize])
  IF STRCMP(direction, 'xy') THEN maskResult_a = ((REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0)
  IF STRCMP(direction, 'yz') THEN maskResult_a = (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
  IF STRCMP(direction, 'xz') THEN maskResult_a = (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)

  maskResult = ROI_b -> ComputeMask(DIMENSIONS = [screenSize, screenSize])
  IF STRCMP(direction, 'xy') THEN maskResult_b = ((REVERSE(CONGRID(maskResult, vol_dims[0], vol_dims[1]), 2)) GT 0)
  IF STRCMP(direction, 'yz') THEN maskResult_b = (CONGRID(maskResult, vol_dims[1], vol_dims[2]) GT 0)
  IF STRCMP(direction, 'xz') THEN maskResult_b = (CONGRID(maskResult, vol_dims[0], vol_dims[2]) GT 0)
  
  maskResult = maskResult_a OR maskResult_b

  ;  iimage, maskResult
  IF STRCMP(direction, 'xy') THEN mask_vol_cube_modify[*, *, (loc_a+i)] = mask_vol_cube_modify[*, *, (loc_a+i)] OR maskResult[*,*]
  IF STRCMP(direction, 'yz') THEN mask_vol_cube_modify[(loc_a+i), *, *] = mask_vol_cube_modify[(loc_a+i), *, *] OR maskResult[*,*]
  IF STRCMP(direction, 'xz') THEN mask_vol_cube_modify[*, (loc_a+i), *] = mask_vol_cube_modify[*, (loc_a+i), *] OR maskResult[*,*]
  
  ;return
  resultGt_a = ROI_a->ComputeGeometry(AREA=tarea_a, CENTROID=tcent_a)
  ROI_a->Translate, -tcent_a
  ROI_a->Scale, 1./Sx_a, 1./Sy_a
  ROI_a->Translate, tcent_a

  resultGt_b = ROI_b->ComputeGeometry(AREA=tarea_b, CENTROID=tcent_b)
  ROI_b->Translate, -tcent_b
  ROI_b->Scale, 1./Sx_b, 1./Sy_b
  ROI_b->Translate, tcent_b

  ROI_a->Translate, -(temp_cent[0] - cent_a[0]), -(temp_cent[1] - cent_a[1])
  ROI_b->Translate, -(temp_cent[0] - cent_b[0]), -(temp_cent[1] - cent_b[1])
  
endfor

RETURN, mask_vol_cube_modify

END

;GJ 2019/5/31, negative threshold may affect segmentation, bug was fixed
;GJ 2019/3/3
PRO lasso_segmentation_3D, xd_sState, old_mask_vol_cube
  threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  threshold_max_value = MAX([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  print, 'threshold_min_value = ', threshold_min_value
  print, 'threshold_max_value = ', threshold_max_value
  thresh_ave = 0.5 * (threshold_min_value + threshold_max_value)
  ;mask the 3D volume
  mask_vol_cube_1 = DOUBLE(xd_sState.mask_vol_cube) * 0. ;zero the mask
  mask_vol_cube_2 = DOUBLE(xd_sState.mask_vol_cube) * 0. ;zero the mask

  FOR imind = 0, N_ELEMENTS(xd_sState.vol_HU_cube[0, 0, *])-1 DO BEGIN
    test_mask_1 = (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0) AND (old_mask_vol_cube[*, *, imind] GT 0)
    mask_vol_cube_1[*, *, imind] = test_mask_1*thresh_ave

    test_mask_2 = (1-(xd_sState.mask_vol_cube_modify[*, *, imind] GT 0)) AND (old_mask_vol_cube[*, *, imind] GT 0)
    mask_vol_cube_2[*, *, imind] = test_mask_2*thresh_ave
  ENDFOR

  SHADE_VOLUME, mask_vol_cube_1,0.9, Outverts, Outconn1

  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;hide old object
    WIDGET_CONTROL, xd_sState.wObjectIndex, GET_VALUE=old_ObjectIndex_value

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_new_object_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0
    IF old_ObjectIndex_value GE 1 THEN BEGIN
      print, 'hiding object =', old_ObjectIndex_value-1
      xd_sState.struc_modifyobject.hide_status[modify_Step, old_ObjectIndex_value-1] = 1
      xd_sState.oPoly[old_ObjectIndex_value-1]->SetProperty, Hide = 1
    ENDIF

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(mask_vol_cube_1 GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
    ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
;    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
;    RETURN
  ENDELSE

  SHADE_VOLUME, mask_vol_cube_2,0.9, Outverts, Outconn1

  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_new_object_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(mask_vol_cube_2 GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

  ENDIF ELSE BEGIN
;    infowarning = DIALOG_MESSAGE('No volume generated! Please reselect!', /ERROR)
;    RETURN
  ENDELSE


END

;GJ 2020/3/22
;generate 3 MIPs, run u-net, load results, generate new vessel results
PRO Add_new_by_AIVesselSegmentation, xd_sState
  ;generate 3 direction MIP
  size_vol = SIZE(xd_sState.vol_HU_cube, /dimensions)
  mip_slice2d = DBLARR(size_vol[0], size_vol[1])

  ;xy plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR j=0, size_vol[1]-1 DO BEGIN
      mip_slice2d[i, j] = MAX(xd_sState.vol_HU_cube[i, j, *])
    ENDFOR
  ENDFOR

  dcm_dir = DIALOG_PICKFILE(PATH=dcm_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select folder to save MIP:')
  FILE_MKDIR, dcm_dir+'MIP_Unet\'
  WRITE_JPEG,  dcm_dir+'MIP_Unet\000.jpeg', BYTSCL(CONGRID(mip_slice2d, 512, 512)), QUALITY=200

  ;yz plane
  FOR j=0, size_vol[1]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[j, k] = MAX(xd_sState.vol_HU_cube[*, j, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\001.jpeg', BYTSCL(CONGRID(mip_slice2d, 512, 512)), QUALITY=200

  ;yz plane
  FOR i=0, size_vol[0]-1 DO BEGIN
    FOR k=0, size_vol[2]-1 DO BEGIN
      mip_slice2d[i, k] = MAX(xd_sState.vol_HU_cube[i, *, k])
    ENDFOR
  ENDFOR
  WRITE_JPEG,  dcm_dir+'MIP_Unet\002.jpeg', BYTSCL(CONGRID(mip_slice2d, 512, 512)), QUALITY=200

  ;call U-net to segment vessel
  ;ZXY
  ;Pleas put code here
  ;

  ;check whether the U-net results are ready
  ;dcm_dir+'Unet'+STRING(m, format='(I3.3)')+'.png'
  a = FILE_SEARCH(dcm_dir, 'Unet*.png', count=nrfile)
  IF nrfile EQ 0 THEN BEGIN
    infowarning = DIALOG_MESSAGE('No U-net result found!', /ERROR)
    RETURN
  ENDIF

  ;if the u-net results are ready, process
  unet_reconstruction_3directions, dcm_dir, vol_HU_cube_mask, Outverts, Outconn1

  ;plot the AI segmwented vessel results
  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    ;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(vol_HU_cube_mask GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState


    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
    ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
  ENDELSE


END


;GJ 2020/3/22
;generate 3 MIPs, run u-net, load results, generate new vessel results
PRO Add_new_by_VWCESegmentation, xd_sState, vessel_mask_vol
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
  
  ;select BBC directory
  IF N_ELEMENTS(BBC_dir) EQ 0 THEN BEGIN
    ;BBC_dir = 'C:\D_drive\AIMIS_3D\images\MRI_BB\2647120_GUO_CAI_QIN_20181009_HR_MRA\S1401_144602_3D_T1WITRA+C\'
    BBC_dir = DIALOG_PICKFILE(PATH=TOF_dir, /MUST_EXIST, /DIRECTORY, TITLE='Select BB+C dcm folder:')
    filearr=file_search(BBC_dir,'*.dcm',count=num)
    IF num EQ 0 THEN RETURN
  ENDIF
  
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
  SHADE_VOLUME, VW_BBC_vol_HU_cube_mask,0.5, Outverts0, Outconn1
  Outverts = MESH_SMOOTH(Outverts0, Outconn1, ITERATIONS=50)
  vol_HU_cube_mask=VW_BBC_vol_HU_cube_mask
  
  ;plot the AI segmwented vessel results
  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    ;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(vol_HU_cube_mask GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState


    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
    ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
  ENDELSE

  
END

;GJ 2023/5/18, adding the registered volume
PRO Add_new_by_thresholding_a1, xd_sState
  
  vol_HU_cube_mask_a1 = (xd_sState.vol_HU_cube_a1 GE xd_sState.threshold_min_value_a1) * (xd_sState.vol_HU_cube_a1 LE xd_sState.threshold_max_value_a1)
  ;  mask_vol_cube = vol_HU_cube_mask * thresh_ave + (1-vol_HU_cube_mask)*(+btm_value-thresh_ave)

  ;delete the variable
  ;  vol_HU_cube_mask = !NULL
  ;
  ;  print, 'max(mask_vol_cube) = ', max(mask_vol_cube)
  ;  print, 'min(mask_vol_cube) = ', min(mask_vol_cube)

  SHADE_VOLUME, vol_HU_cube_mask_a1, 0.9, Outverts, Outconn1
  ;  SHADE_VOLUME, mask_vol_cube,threshold_min_value, Outverts, Outconn1

  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    ;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(vol_HU_cube_mask_a1 GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution1^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
    ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
  ENDELSE

END


;GJ 2019/3/2
PRO Add_new_by_thresholding, xd_sState
  ;segmentation
  threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  threshold_max_value = MAX([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  print, 'threshold_min_value = ', threshold_min_value
  print, 'threshold_max_value = ', threshold_max_value
;  thresh_ave = 0.5 * (threshold_min_value + threshold_max_value)
;  btm_value = MIN([MIN(xd_sState.vol_HU_cube), -1000])
;  print, 'btm_value', btm_value
  ;mask the 3D volume
;  mask_vol_cube = DOUBLE(xd_sState.mask_vol_cube) * 0. ;zero the mask
  
  vol_HU_cube_mask = (xd_sState.vol_HU_cube GE threshold_min_value) * (xd_sState.vol_HU_cube LE threshold_max_value)
;  mask_vol_cube = vol_HU_cube_mask * thresh_ave + (1-vol_HU_cube_mask)*(+btm_value-thresh_ave)

  ;delete the variable
;  vol_HU_cube_mask = !NULL
;  
;  print, 'max(mask_vol_cube) = ', max(mask_vol_cube)
;  print, 'min(mask_vol_cube) = ', min(mask_vol_cube)
  
  SHADE_VOLUME, vol_HU_cube_mask, 0.9, Outverts, Outconn1
;  SHADE_VOLUME, mask_vol_cube,threshold_min_value, Outverts, Outconn1
  
  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    ;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(vol_HU_cube_mask GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
  ENDELSE

END

;GJ 2019/3/2
PRO Add_new_by_readstl, xd_sState, Outconn1, smoothedOutverts
  IF N_ELEMENTS(smoothedOutverts) GT 10 THEN BEGIN
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    ;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    vol_HU_cube_mask = BYTARR(256, 256, 256)*0.
    FOR ii=0, (SIZE(smoothedOutverts, /DIM))[1]-1 DO vol_HU_cube_mask[smoothedOutverts[0, ii], smoothedOutverts[1, ii], smoothedOutverts[2, ii]] = 1
    ; Create the shape operator:
    S = REPLICATE(1, 10, 10, 10)
    ; "Dilate" operator:
    vol_HU_cube_mask = DILATE(vol_HU_cube_mask, S)
    
    index = WHERE(vol_HU_cube_mask EQ 1, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
    ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
  ENDELSE

END

;name: mask_segmentation
;input: 
PRO Mask_segmentation, xd_sState, modify_roi

  
  ;segmentation
  threshold_min_value = MIN([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  threshold_max_value = MAX([xd_sState.sel_one.HU_min_value, xd_sState.sel_one.HU_max_value])
  print, 'threshold_min_value = ', threshold_min_value
  print, 'threshold_max_value = ', threshold_max_value
  thresh_ave = 0.5 * (threshold_min_value + threshold_max_value)
  ;mask the 3D volume
  mask_vol_cube = DOUBLE(xd_sState.mask_vol_cube) * 0. ;zero the mask
  
  btm_value = MIN([MIN(xd_sState.vol_HU_cube), -1000])
  
;  xd_sState.xy_mask = xd_sState.xy_mask_ROI * xd_sState.xy_mask_clip
;  xy_index = WHERE(xd_sState.xy_mask GT 0, xy_count)
;  xd_sState.yz_mask = xd_sState.yz_mask_ROI * xd_sState.yz_mask_clip
;  yz_index = WHERE(xd_sState.yz_mask GT 0, yz_count)
;  xd_sState.xz_mask = xd_sState.xz_mask_ROI * xd_sState.xz_mask_clip
;  xz_index = WHERE(xd_sState.xz_mask GT 0, xz_count)
  
  IF N_ELEMENTs(modify_roi) THEN BEGIN
    
    print,'modify_roi = ', modify_roi
    FOR imind = 0, N_ELEMENTS(xd_sState.vol_HU_cube[0, 0, *])-1 DO BEGIN
      ;border
      IF ABS(modify_roi-1) LE 1 THEN BEGIN
        test_mask = (xd_sState.mask_vol_cube[*, *, imind] GT 0) AND (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0)
        mask_vol_cube[*, *, imind] = xd_sState.vol_HU_cube[*, *, imind] * test_mask + (1-test_mask)*btm_value;(-1000) ;  + (test_mask *1000 -1000)
        
        ;delete the variable
        test_mask = !NULL
      ENDIF
      
      ;add
      IF ABS(modify_roi-4) LE 1 THEN BEGIN
        test_mask = (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0) AND (1-(xd_sState.new_mask_vol_cube[*, *, imind] GT 0))
        ;mask_vol_cube[*, *, imind] = xd_sState.vol_HU_cube[*, *, imind] * (xd_sState.new_mask_vol_cube[*, *, imind] GT 0) + test_mask*thresh_ave
        mask_vol_cube[*, *, imind] = ((xd_sState.new_mask_vol_cube[*, *, imind] GT 0) + test_mask)*thresh_ave
        
        ;delete the variable
        test_mask = !NULL
      ENDIF

      ;remove selected region and 3d segmentation 3dseg
      IF ABS(modify_roi-8) LE 1 THEN BEGIN
        mask_vol_cube[*, *, imind] = (xd_sState.new_mask_vol_cube[*, *, imind] GT 0)*thresh_ave + (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0)*(btm_value-thresh_ave) + (xd_sState.new_mask_vol_cube[*, *, imind] LE 0)*(btm_value-thresh_ave)
      ENDIF
      
      ;add red object
      IF ABS(modify_roi-12) LE 1 THEN BEGIN
        test_mask = (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0) AND (1-(xd_sState.new_mask_vol_cube[*, *, imind] GT 0))
        ;mask_vol_cube[*, *, imind] = xd_sState.vol_HU_cube[*, *, imind] * (xd_sState.new_mask_vol_cube[*, *, imind] GT 0) + test_mask*thresh_ave
        mask_vol_cube[*, *, imind] = ((xd_sState.new_mask_vol_cube[*, *, imind] GT 0) + test_mask)*thresh_ave
        
        ;delete the variable
        test_mask = !NULL
      ENDIF
      
      ;add green object
      IF ABS(modify_roi-16) LE 1 THEN BEGIN
        test_mask = (xd_sState.mask_vol_cube_modify[*, *, imind] GT 0) AND (1-(xd_sState.new_mask_vol_cube[*, *, imind] GT 0))
        ;mask_vol_cube[*, *, imind] = xd_sState.vol_HU_cube[*, *, imind] * (xd_sState.new_mask_vol_cube[*, *, imind] GT 0) + test_mask*thresh_ave
        mask_vol_cube[*, *, imind] = ((xd_sState.new_mask_vol_cube[*, *, imind] GT 0) + test_mask)*thresh_ave
        
        ;delete the variable
        test_mask = !NULL
      ENDIF
     
      ;3d region grow 3drg
      IF ABS(modify_roi-20) LE 1 THEN BEGIN
        mask_vol_cube[*, *, imind] = (xd_sState.new_mask_vol_cube[*, *, imind] GT 0)*thresh_ave + (xd_sState.new_mask_vol_cube[*, *, imind] LE 0)*(btm_value-thresh_ave)
      ENDIF
      
     
     ENDFOR
    
    IF MAX(xd_sState.clippoint) GT 0 THEN BEGIN
      print, 'xd_sState.clippoint = ', xd_sState.clippoint
      print, 'size(mask_vol_cube) = ', size(mask_vol_cube)
      mask_vol_cube[0:(xd_sState.clippoint[0]-1), *, *] = mask_vol_cube[0:(xd_sState.clippoint[0]-1), *, *]*0. +btm_value
      mask_vol_cube[(xd_sState.clippoint[1]-3):(size(mask_vol_cube))[1]-1, *, *] = mask_vol_cube[(xd_sState.clippoint[1]-3):(size(mask_vol_cube))[1]-1, *, *]*0. +btm_value
      mask_vol_cube[*, 0:(xd_sState.clippoint[2]-1), *] = mask_vol_cube[*, 0:(xd_sState.clippoint[2]-1), *]*0. +btm_value
      mask_vol_cube[*, (xd_sState.clippoint[3]-3):(size(mask_vol_cube))[2]-1, *] = mask_vol_cube[*, (xd_sState.clippoint[3]-3):(size(mask_vol_cube))[2]-1, *]*0. +btm_value
      mask_vol_cube[*, *, 0:(xd_sState.clippoint[4]-1)] = mask_vol_cube[*, *, 0:(xd_sState.clippoint[4]-1)]*0. +btm_value
      mask_vol_cube[*, *, (xd_sState.clippoint[5]-3):(size(mask_vol_cube))[3]-1] = mask_vol_cube[*, *, (xd_sState.clippoint[5]-3):(size(mask_vol_cube))[3]-1]*0. +btm_value
    ENDIF ELSE BEGIN
      mask_vol_cube[0, *, *] = mask_vol_cube[0, *, *]*0. +btm_value
      mask_vol_cube[(size(mask_vol_cube))[1]-1, *, *] = mask_vol_cube[(size(mask_vol_cube))[1]-1, *, *]*0. +btm_value
      mask_vol_cube[*, 0, *] = mask_vol_cube[*, 0, *]*0. +btm_value
      mask_vol_cube[*, (size(mask_vol_cube))[2]-1, *] = mask_vol_cube[*, (size(mask_vol_cube))[2]-1, *]*0. +btm_value
      mask_vol_cube[*, *, 0] = mask_vol_cube[*, *, 0]*0. +btm_value
      mask_vol_cube[*, *, (size(mask_vol_cube))[3]-1] = mask_vol_cube[*, *, (size(mask_vol_cube))[3]-1]*0. +btm_value
    ENDELSE
    
    ;threshold modifying
    IF ABS(modify_roi-24) LE 1 THEN BEGIN
      vol_HU_cube_mask = (xd_sState.vol_HU_cube GE threshold_min_value) * (xd_sState.vol_HU_cube LE threshold_max_value)
;      mask_vol_cube = vol_HU_cube_mask * xd_sState.mask_vol_cube * thresh_ave + (1-vol_HU_cube_mask * xd_sState.mask_vol_cube)*(+btm_value-thresh_ave)
      mask_vol_cube = vol_HU_cube_mask * thresh_ave + (1-vol_HU_cube_mask)*(+btm_value-thresh_ave)
      
      ;delete the variable
      vol_HU_cube_mask = !NULL
;      mask_vol_cube = xd_sState.vol_HU_cube * xd_sState.mask_vol_cube + (1-xd_sState.mask_vol_cube)*(-1000)
;      xd_sState.oPoly1->getProperty, DATA = Outverts_old
;      xd_sState.oPoly1->setProperty, DATA = Outverts_old*0.0001
;      xd_sState.vol_red = 0.
;      print, 'volume (red) = ', xd_sState.vol_red, ' cm3'
;      Stringred = STRING(xd_sState.vol_red,FORMAT='("vol(red) =", F9.3, " cm3")')
;      xd_sState.tvol_red.SetProperty,strings=stringred
    ENDIF
   
    SHADE_VOLUME, mask_vol_cube,threshold_min_value, Outverts, Outconn1

    ;change mask back to 0 and 1
    xd_sState.new_mask_vol_cube = (mask_vol_cube GT threshold_min_value) AND (mask_vol_cube LT threshold_max_value)
    
    ;not for region growing/threshold modifying
    IF modify_roi LT 19 THEN BEGIN
      ;GJ 2018/12/21, keeping removing part
      mask_vol_oPoly1 = xd_sState.mask_vol_cube * (1-xd_sState.new_mask_vol_cube); * xd_sState.mask_vol_cube
;      iimage, mask_vol_oPoly1[*,*,128] 
      SHADE_VOLUME, mask_vol_oPoly1,0.9, Outverts1, Outconn11
      
      IF N_ELEMENTS(Outverts1) GT 10 THEN BEGIN
        index = WHERE(mask_vol_oPoly1, count)
;        IF count NE 0 THEN xd_sState.vol_red = count * (xd_sState.vol_HU_cube_resolution^3)/1000.
;        print, 'volume (red) = ', xd_sState.vol_red, ' cm3'
;        Stringred = STRING(xd_sState.vol_red,FORMAT='("vol(red) =", F9.3, " cm3")')
;        xd_sState.tvol_red.SetProperty,strings=stringred

        smoothedOutverts1 = MESH_SMOOTH(Outverts1, Outconn11, ITERATIONS=50)
        ;smoothedOutverts1 = Outverts1
        ;plot segmentation result
        xd_sState.oPoly1->getProperty, DATA = Outverts_old
        ;define to the rotation center
        IF N_ELEMENTS(Outverts_old) LT N_ELEMENTS(smoothedOutverts1) THEN BEGIN
          xd_sState.oPoly1->setProperty, DATA = smoothedOutverts1
          xd_sState.oPoly1->setProperty, POLYGONS = Outconn11
        ENDIF ELSE BEGIN
          xd_sState.oPoly1->setProperty, POLYGONS = Outconn11
          xd_sState.oPoly1->setProperty, DATA = smoothedOutverts1
        ENDELSE

        xd_sstate.opoly[0]->getProperty, XCOORD_CONV=xc, YCOORD_CONV=yc, ZCOORD_CONV=zc
        xd_sState.oPoly1->setProperty, XCOORD_CONV=xc, YCOORD_CONV=yc, ZCOORD_CONV=zc
        xd_sState.oPoly2->setProperty, XCOORD_CONV=xc, YCOORD_CONV=yc, ZCOORD_CONV=zc

        widget_control, xd_sState.hSlice_x_slider, get_value = xvalue
        widget_control, xd_sState.hSlice_y_slider, get_value = yvalue
        widget_control, xd_sState.hSlice_z_slider, get_value = zvalue

        xd_sState.opoly1->SetProperty, $
          CLIP_PLANES=[-xvalue, -yvalue, -zvalue, (xvalue^2 + yvalue^2 + zvalue^2)];, ALPHA_CHANNEL=0.5


        WIDGET_CONTROL, xd_sState.x_scale_slider, sensitive = 1
        WIDGET_CONTROL, xd_sState.y_scale_slider, sensitive = 1
        WIDGET_CONTROL, xd_sState.z_scale_slider, sensitive = 1
        WIDGET_CONTROL, xd_sState.x_move_slider, sensitive = 1
        WIDGET_CONTROL, xd_sState.y_move_slider, sensitive = 1
        WIDGET_CONTROL, xd_sState.z_move_slider, sensitive = 1
        WIDGET_CONTROL, xd_sState.xyzscalemove_output, sensitive = 1
        
;        WIDGET_CONTROL, xd_sState.wSmooth1, sensitive = 1
;        WIDGET_CONTROL, xd_sState.wtransSlider1, sensitive = 1
;        WIDGET_CONTROL, xd_sState.wloadctSlider1, sensitive = 1
        
;        oPoly1 = OBJ_NEW('IDLgrPolygon', COLOR = [255, 127, 127], smoothedOutverts1, POLYGON=Outconn11)
;        mymodel1 = OBJ_NEW('IDLgrModel')
;        mymodel1->Add, oPoly1
;        ;
;        xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]
        
      ENDIF
;      xd_sState.mask_vol_cube = xd_sState.new_mask_vol_cube
    ENDIF

    ;IF ABS(modify_roi-24) GT 1 THEN 
    xd_sState.mask_vol_cube = xd_sState.new_mask_vol_cube
    
    
    index = WHERE(xd_sState.new_mask_vol_cube, count)
    IF count NE 0 THEN BEGIN
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000.
;      print, 'volume = ', xd_sState.vol_white, ' cm3'
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      xd_sState.tvol_white.SetProperty,strings=stringwhite

      index3d =  ARRAY_INDICES(xd_sState.new_mask_vol_cube, index)
      center = FLTARR(3)
      center[0] = MEAN(index3d[0,*])
      center[1] = MEAN(index3d[1,*])
      center[2] = MEAN(index3d[2,*])
      ;      center = findcenter(xd_sState.new_mask_vol_cube)
      print, 'white volume center1 = ', center
      x_range=[MIN(index3d[0,*]), MAX(index3d[0,*])]
      print, 'white volume x range = ', x_range
      y_range=[MIN(index3d[1,*]), MAX(index3d[1,*])]
      print, 'white volume y range = ', y_range
      z_range=[MIN(index3d[2,*]), MAX(index3d[2,*])]
      print, 'white volume z range = ', z_range
    ENDIF

  ENDIF ELSE BEGIN
    
    print,'222'
    mask_vol_cube = xd_sState.vol_HU_cube
    ;reduce the dimension of vol_HU_cube
;    dims = size(xd_sState.vol_HU_cube)
;    vol_HU_cube_sm = CONGRID(xd_sState.vol_HU_cube, dims[1]/2, dims[2]/2, dims[3]/2)
    vol_HU_cube_sm = xd_sState.vol_HU_cube
    vol_HU_cube_sm = (vol_HU_cube_sm GE threshold_min_value) * (vol_HU_cube_sm LE threshold_max_value)

    vol_HU_cube_sm =  vol_HU_cube_sm * (threshold_min_value+1-btm_value) +btm_value ; * (vol_HU_cube_sm LT threshold_max_value)

    IF MAX(xd_sState.clippoint) GT 0 THEN BEGIN
      vol_HU_cube_sm[0:(xd_sState.clippoint[0]-1), *, *] = vol_HU_cube_sm[0:(xd_sState.clippoint[0]-1), *, *]*0. +btm_value
      vol_HU_cube_sm[(xd_sState.clippoint[1]-3):(size(vol_HU_cube_sm))[1]-1, *, *] = vol_HU_cube_sm[(xd_sState.clippoint[1]-3):(size(vol_HU_cube_sm))[1]-1, *, *]*0. +btm_value
      vol_HU_cube_sm[*, 0:(xd_sState.clippoint[2]-1), *] = vol_HU_cube_sm[*, 0:(xd_sState.clippoint[2]-1), *]*0. +btm_value
      vol_HU_cube_sm[*, (xd_sState.clippoint[3]-3):(size(vol_HU_cube_sm))[2]-1, *] = vol_HU_cube_sm[*, (xd_sState.clippoint[3]-3):(size(vol_HU_cube_sm))[2]-1, *]*0. +btm_value
      vol_HU_cube_sm[*, *, 0:(xd_sState.clippoint[4]-1)] = vol_HU_cube_sm[*, *, 0:(xd_sState.clippoint[4]-1)]*0. +btm_value
      vol_HU_cube_sm[*, *, (xd_sState.clippoint[5]-3):(size(vol_HU_cube_sm))[3]-1] = vol_HU_cube_sm[*, *, (xd_sState.clippoint[5]-3):(size(vol_HU_cube_sm))[3]-1]*0. +btm_value
    ENDIF ELSE BEGIN
      vol_HU_cube_sm[0, *, *] = vol_HU_cube_sm[0, *, *]*0. +btm_value
      vol_HU_cube_sm[(size(vol_HU_cube_sm))[1]-1, *, *] = vol_HU_cube_sm[(size(vol_HU_cube_sm))[1]-1, *, *]*0. +btm_value
      vol_HU_cube_sm[*, 0, *] = vol_HU_cube_sm[*, 0, *]*0. +btm_value
      vol_HU_cube_sm[*, (size(vol_HU_cube_sm))[2]-1, *] = vol_HU_cube_sm[*, (size(vol_HU_cube_sm))[2]-1, *]*0. +btm_value
      vol_HU_cube_sm[*, *, 0] = vol_HU_cube_sm[*, *, 0]*0. +btm_value
      vol_HU_cube_sm[*, *, (size(vol_HU_cube_sm))[3]-1] = vol_HU_cube_sm[*, *, (size(vol_HU_cube_sm))[3]-1]*0. +btm_value
    ENDELSE
    SHADE_VOLUME, vol_HU_cube_sm,threshold_min_value, Outverts, Outconn1
    
    ;delete the variable
    vol_HU_cube_sm = !NULL
    
    ;not for threshold modifying
    IF ABS(modify_roi-24) GT 1 THEN xd_sState.mask_vol_cube = mask_vol_cube
    
    xd_sState.new_mask_vol_cube = (mask_vol_cube GT threshold_min_value) AND (mask_vol_cube LT threshold_max_value); mask_vol_cube GT -50000.
    
    ;delete the variable
    mask_vol_cube = !NULL
    
    index = WHERE(xd_sState.new_mask_vol_cube, count)
    IF count NE 0 THEN BEGIN
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000.
;      print, 'volume (white) = ', xd_sState.vol_white, ' cm3'
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      xd_sState.tvol_white.SetProperty,strings=stringwhite
      
      index3d =  ARRAY_INDICES(xd_sState.new_mask_vol_cube, index)
      center = FLTARR(3)
      center[0] = MEAN(index3d[0,*])
      center[1] = MEAN(index3d[1,*])
      center[2] = MEAN(index3d[2,*])
      ;      center = findcenter(xd_sState.new_mask_vol_cube)
      print, 'white volume center1 = ', center
      x_range=[MIN(index3d[0,*]), MAX(index3d[0,*])]
      print, 'white volume x range = ', x_range
      y_range=[MIN(index3d[1,*]), MAX(index3d[1,*])]
      print, 'white volume y range = ', y_range
      z_range=[MIN(index3d[2,*]), MAX(index3d[2,*])]
      print, 'white volume z range = ', z_range
    ENDIF

  ENDELSE
  
  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]
    
    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1
    
    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value
    
;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0
    
    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]
    
    ;define mask_vol_cube
    index = WHERE(mask_vol_cube GT threshold_min_value, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF
    
    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=xd_sState, /NO_COPY
        RETURN
      ENDIF
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState

    
    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')
    
    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]
    
    ;enable these options
;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1
  
    IF N_ELEMENTS(modify_roi) THEN BEGIN
      widget_control, xd_sState.wUndo, SENSITIVE = 1
      widget_control, xd_sState.wRedo, SENSITIVE = 1
    ENDIF
    
 ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
 ENDELSE
 
END

;adjust whole drawing
PRO Adjust_whole_polygons, smoothedOutverts, xd_sState
  ;define to the rotation center
  ma=max(smoothedOutverts)
  mi=min(smoothedOutverts)
  bias=(ma-mi)/2
  rl=mi+bias
  rr=ma+bias
  cc=[(-rl)/(rr-rl),1/(rr-rl)]
  xma = max((smoothedOutverts)[0,*])
  xmi = min((smoothedOutverts)[0,*])
  xmid = 0.5*(xma+xmi)*cc[1]
  yma = max((smoothedOutverts)[1,*])
  ymi = min((smoothedOutverts)[1,*])
  ymid = 0.5*(yma+ymi)*cc[1]
  zma = max((smoothedOutverts)[2,*])
  zmi = min((smoothedOutverts)[2,*])
  zmid = 0.5*(zma+zmi)*cc[1]

  print, 'cc = ', cc
  
  ;adjust the center of all 3D objects
  n_Poly = xd_sState.struc_modifyobject.exist_num[xd_sState.struc_modifyobject.modify_Step]
  FOR i = 0, n_Poly-1 DO BEGIN
    xd_sstate.opoly[i]->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  ENDFOR
  
  xd_sState.oPoly1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oPoly2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oPoly_border->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oS0->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oS1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oS2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713
  xd_sState.ogreenaxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oyellowaxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oblueAxis1->setProperty,  XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx0713

  xd_sState.oLine1->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oLine2->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]
  xd_sState.oLine3->setProperty, XCOORD_CONV=[-xmid, cc[1]], YCOORD_CONV=[-ymid, cc[1]], ZCOORD_CONV=[-zmid, cc[1]]

  ;adjust rulers
  ;[MAX(sel_one.verts[0,*])-MIN(sel_one.verts[0,*])
  scale = 2.3*280./((xma-xmi)*xd_sState.vol_HU_cube_resolution)
  xd_sState.oRuler1model->setproperty,transform=[[scale, 0, 0, 0.0],[0, 1, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
  xd_sState.oRuler2model->setproperty,transform=[[1, 0, 0, 0.0],[0, scale, 0, 0.0],[0, 0, 1, 0.0],[0, 0, 0, 1.0]]
  xd_sState.oRuler1Model->scale,0.008,0.05,1
  xd_sState.oRuler1Model->translate,-0.2,0.4,0.45
  xd_sState.oRuler2Model->scale,0.05,0.008,1
  xd_sState.oRuler2Model->translate,0.4,-0.2,0.45
  ;    ;update sel_one
  xd_sState.sel_one.xcoord_conv = [-xmid, cc[1]]
  xd_sState.sel_one.ycoord_conv = [-ymid, cc[1]]
  xd_sState.sel_one.zcoord_conv = [-zmid, cc[1]]

;  smoothvalue = 12.5
;  WIDGET_CONTROL, xd_sState.wSmooth, SET_VALUE=smoothvalue
;  SmoothString = 'Smooth : ' + STRING(smoothvalue, $
;    FORMAT='(f5.1)') + ' %'
;  print,SmoothString
;  WIDGET_CONTROL, xd_sState.wSmooth_text, SET_VALUE=STRING(smoothvalue, format='(f5.1)')

END


;濠⒀呭仜婵偤寮甸敓鑺ユ緰闁汇劌瀚伴鍥ㄥ緞閺夋垵琚遍幖杈炬嫹
;input: image,濞存粌鐬煎ǎ顕�炊閸撗冾暬, binary type
;output: newimage,濠⒀呭仜鐢倝宕ユ惔锝嗙暠濞存粌鐬煎ǎ顕�炊閸撗冾暬, binary type
FUNCTION ThickIncrease, binaryImg

  radius = 7
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
  binarImg = MORPH_OPEN(binaryImg, strucElem)

  
  ; Get the image size, prepare a display window and
  ; display the image.
  dims = SIZE(binaryImg, /DIMENSIONS)
  ;WINDOW, 0, XSIZE = 2*dims[0], YSIZE = 2*dims[1], $
  ;   TITLE = 'Original, Binary and Thinned Images'
  ;TVSCL, img, 0

  ;; Generate a binary image by thresholding.
  ;binaryImg = img GE 140
  ;TVSCL, binaryImg, 1

; Prepare hit and miss structures for thinning.
h0 = [[0b,0,0], $
  [0,1,0], $
  [1,1,1]]
m0 = [[1b,1,1], $
  [0,0,0], $
  [0,0,0]]
h1 = [[0b,0,0], $
  [1,1,0], $
  [1,1,0]]
m1 = [[0b,1,1], $
  [0,0,1], $
  [0,0,0]]
h2 = [[1b,0,0], $
  [1,1,0], $
  [1,0,0]]
m2 = [[0b,0,1], $
  [0,0,1], $
  [0,0,1]]
h3 = [[1b,1,0], $
  [1,1,0], $
  [0,0,0]]
m3 = [[0b,0,0], $
  [0,0,1], $
  [0,1,1]]
h4 = [[1b,1,1], $
  [0,1,0], $
  [0,0,0]]
m4 = [[0b,0,0], $
  [0,0,0], $
  [1,1,1]]
h5 = [[0b,1,1], $
  [0,1,1], $
  [0,0,0]]
m5 = [[0b,0,0], $
  [1,0,0], $
  [1,1,0]]
h6 = [[0b,0,1], $
  [0,1,1], $
  [0,0,1]]
m6 = [[1b,0,0], $
  [1,0,0], $
  [1,0,0]]
h7 = [[0b,0,0], $
  [0,1,1], $
  [0,1,1]]
m7 = [[1b,1,0], $
  [1,0,0], $
  [0,0,0]]
  
; Iterate until the thinned image is identical to
; the input image for a given iteration.
bCont = 1b
iIter = 1
thinImg = binaryImg
thickImg = binaryImg
WHILE bCont EQ 1b DO BEGIN
   PRINT,'Iteration: ', iIter
   inputImg = thinImg

   ; Perform the thinning using the first pair
   ; of structure elements.
   thinImg = MORPH_THIN(inputImg, h0, m0)

   ; Perform the thinning operation using the
   ; remaining structural element pairs.
   thinImg = MORPH_THIN(thinImg, h1, m1)
   thinImg = MORPH_THIN(thinImg, h2, m2)
   thinImg = MORPH_THIN(thinImg, h3, m3)
   thinImg = MORPH_THIN(thinImg, h4, m4)
   thinImg = MORPH_THIN(thinImg, h5, m5)
   thinImg = MORPH_THIN(thinImg, h6, m6)
   thinImg = MORPH_THIN(thinImg, h7, m7)

   ; Display the results of thinning and wait a second for
   ; display purposes.
;   TVSCL, thinImg, 2
;   WAIT, 1

   ; 閻犱焦婢樼紞宥夊储濮橆剙顔�
   thickImg = thickImg + inputImg
   
   ; Test the condition and increment the loop.
   bCont = MAX(inputImg - thinImg)
   iIter = iIter + 1

; End WHILE loop statements.
ENDWHILE
  ;;;
  ;;;; the start of removing samll bars in midline
  
;  iimage, binaryImg
;  iimage, thinImg
  ;Use Count to get the number of nonzero elements
  index=WHERE(thinImg, count)
  
  radius = 2
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
;  PRINT, strucElem
  dilateImg = DILATE(thinImg, strucElem)
  
  thickerImg = dilateImg OR binaryImg
;  return, thickerImg
  return, binaryImg

END

;GJ 2019/01/29
;generating 4 images (front, back, left, right) and packing together for naked eye 3D view
;
PRO draw_3DView, xd_sState

  imagesize_3DView = xd_sState.imageSize_3DView
  odraw_img = xd_sState.oWindow->READ()
  odraw_img->getProperty, DATA = draw_img_tmp
  draw_img = CUT_CORNER(draw_img_tmp)
  ;iimage, draw_img
  front_img_t = CONGRID(draw_img, 3, 0.5*imagesize_3DView,  0.5*imagesize_3DView, /INTERP)
  front_img = REVERSE(front_img_t, 3, /OVERWRITE)

  ;GJ 2019/1/28
  axis = [0,1,0]
  angle = 90
  xd_sState.oScalingModel->Rotate, axis, angle
  ;end GJ 2019/1/28
  xd_sState.oWindow->Draw, xd_sState.oView
  odraw_img = xd_sState.oWindow->READ()
  odraw_img->getProperty, DATA = draw_img_tmp
  draw_img = CUT_CORNER(draw_img_tmp)
  left_img_t = CONGRID(draw_img, 3, 0.5*imagesize_3DView,  0.5*imagesize_3DView, /INTERP)
  left_img = left_img_t
  left_img[0,*,*] = ROTATE(REFORM(left_img_t[0,*,*]), 3)
  left_img[1,*,*] = ROTATE(REFORM(left_img_t[1,*,*]), 3)
  left_img[2,*,*] = ROTATE(REFORM(left_img_t[2,*,*]), 3)
  left_img = REVERSE(left_img, 2, /OVERWRITE)

  ;GJ 2019/1/28
  axis = [0,1,0]
  angle = 90
  xd_sState.oScalingModel->Rotate, axis, angle
  ;end GJ 2019/1/28
  xd_sState.oWindow->Draw, xd_sState.oView
  odraw_img = xd_sState.oWindow->READ()
  odraw_img->getProperty, DATA = draw_img_tmp
  draw_img = CUT_CORNER(draw_img_tmp)
  back_img = CONGRID(draw_img, 3, 0.5*imagesize_3DView,  0.5*imagesize_3DView, /INTERP)
;  back_img = REVERSE(back_img_t, 3, /OVERWRITE)

  ;GJ 2019/1/28
  axis = [0,1,0]
  angle = 90
  xd_sState.oScalingModel->Rotate, axis, angle
  ;end GJ 2019/1/28
  xd_sState.oWindow->Draw, xd_sState.oView
  odraw_img = xd_sState.oWindow->READ()
  odraw_img->getProperty, DATA = draw_img_tmp
  draw_img = CUT_CORNER(draw_img_tmp)
  right_img_t = CONGRID(draw_img, 3, 0.5*imagesize_3DView,  0.5*imagesize_3DView, /INTERP)
  right_img = right_img_t
  right_img[0,*,*] = ROTATE(REFORM(right_img_t[0,*,*]), 1)
  right_img[1,*,*] = ROTATE(REFORM(right_img_t[1,*,*]), 1)
  right_img[2,*,*] = ROTATE(REFORM(right_img_t[2,*,*]), 1)
  right_img = REVERSE(right_img, 2, /OVERWRITE)
  right_img = REVERSE(right_img, 3, /OVERWRITE)
  
  ;check img size
  img_size = SIZE(left_img, /dim)
  big_img_l = BYTARR(3, imagesize_3DView, imagesize_3DView)*0
  big_img_l[*, 0:(img_size[1]-1), FLOOR(imagesize_3DView/4):(FLOOR(imagesize_3DView/4)+img_size[2]-1)]=left_img
  big_img_r = BYTARR(3, imagesize_3DView, imagesize_3DView)*0
  big_img_r[*, FLOOR(imagesize_3DView/2):(FLOOR(imagesize_3DView/2)+img_size[1]-1), FLOOR(imagesize_3DView/4):(FLOOR(imagesize_3DView/4)+img_size[2]-1)]=right_img
  big_img_f = BYTARR(3, imagesize_3DView, imagesize_3DView)*0
  big_img_f[*, FLOOR(imagesize_3DView/4):(FLOOR(imagesize_3DView/4)+img_size[1]-1), 0:(img_size[2]-1)]=front_img
  big_img_b = BYTARR(3, imagesize_3DView, imagesize_3DView)*0
  big_img_b[*, FLOOR(imagesize_3DView/4):(FLOOR(imagesize_3DView/4)+img_size[1]-1), FLOOR(imagesize_3DView/2):(FLOOR(imagesize_3DView/2)+img_size[2]-1)]=back_img
  big_img = big_img_l + big_img_r + big_img_f + big_img_b
 
  WSET, xd_sState.oWindow_3DView
  tv, big_img, TRUE=1
;  tv, front_img, TRUE=1, 3
;  tv, left_img, TRUE=1, 2
;  tv, back_img, TRUE=1, 0
;  tv, right_img, TRUE=1, 1


  ;GJ 2019/1/28
  axis = [0,1,0]
  angle = 90
  xd_sState.oScalingModel->Rotate, axis, angle
  ;end GJ 2019/1/28
  xd_sState.oWindow->Draw, xd_sState.oView

  ;GJ 2019/1/28

END

;pick vert colors based on a new matrix ;GJ 2019/2/13
;vert_colors = PICK_VERT_COLORS_thickness(vol_HU_cube, sel_one.verts, rgb_table)
FUNCTION PICK_VERT_COLORS_by_thickness, vol_cube_mask, verts, rgb_table

  IF N_ELEMENTS(verts) EQ 0 THEN RETURN, 0

  ;define a color table
  ;redvalue = BYTSCL(RANDOMN(seed, 256))
  ;redvalue = replicate(245,256,1)
  ;greenvalue = BINDGEN(256)
  ;bluevalue = reverse(greenvalue)
  LOADCT, 13, RGB_TABLE = rgb_table

  mask_dim = SIZE(vol_cube_mask, /DIM)
  min_thick = 20

  radius = 2
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius

  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='pink', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start

  ;scale vol cube to 0-255
  temp_vol_cube_mask_xy = BYTSCL(vol_cube_mask)
  temp_vol_thin_xy = BYTSCL(vol_cube_mask)*0.
  FOR i = 0, mask_dim[2]-1 DO BEGIN
    count=(i+1.)/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;i = mask_dim[2]/2
    mask_img = DILATE(temp_vol_cube_mask_xy[*,*,i], strucElem)
    result = LABEL_REGION(mask_img, /ALL_NEIGHBORS)
    thinImg = mask_img*0.
    FOR j = 1, MAX(result) DO BEGIN
      temp_mask_img = mask_img*0
      index = WHERE(result EQ j, count)
      IF count NE 0 THEN temp_mask_img[index] = 1
      thick = image_thinning(temp_mask_img)
      ;    print, 'thick = ', thick
      mask_img[index] = thick
      IF thick LT min_thick THEN min_thick = thick
    ENDFOR
    temp_vol_cube_mask_xy[*,*,i] = mask_img
    temp_vol_thin_xy[*,*,i] = DILATE(thinImg, strucElem)
  ENDFOR

  ;scale vol cube to 0-255
  temp_vol_cube_mask_yz = BYTSCL(vol_cube_mask)
  temp_vol_thin_yz = BYTSCL(vol_cube_mask)*0.
  FOR i = 0, mask_dim[0]-1 DO BEGIN
    count=(i+1.+mask_dim[2])/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;i = mask_dim[2]/2
    mask_img = DILATE(REFORM(temp_vol_cube_mask_yz[i,*,*]), strucElem)
    thinImg = mask_img*0.
    result = LABEL_REGION(mask_img, /ALL_NEIGHBORS)
    FOR j = 1, MAX(result) DO BEGIN
      temp_mask_img = mask_img*0
      index = WHERE(result EQ j, count)
      IF count NE 0 THEN temp_mask_img[index] = 1
      thick = image_thinning(temp_mask_img)
      mask_img[index] = thick
      IF thick LT min_thick THEN min_thick = thick
    ENDFOR
    temp_vol_cube_mask_yz[i,*,*] = mask_img
    temp_vol_thin_yz[i,*,*] = DILATE(thinImg, strucElem)
  ENDFOR

  ;scale vol cube to 0-255
  temp_vol_cube_mask_xz = BYTSCL(vol_cube_mask)
  temp_vol_thin_xz = BYTSCL(vol_cube_mask)*0.
  FOR i = 0, mask_dim[1]-1 DO BEGIN
    count=(i+1.+mask_dim[2]+mask_dim[1])/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;i = mask_dim[2]/2
    mask_img = DILATE(REFORM(temp_vol_cube_mask_xz[*,i,*]), strucElem)
    thinImg = mask_img*0.
    result = LABEL_REGION(mask_img, /ALL_NEIGHBORS)
    FOR j = 1, MAX(result) DO BEGIN
      temp_mask_img = mask_img*0
      index = WHERE(result EQ j, count)
      IF count NE 0 THEN temp_mask_img[index] = 1
      thick = image_thinning(temp_mask_img)
      mask_img[index] = thick
      IF thick LT min_thick THEN min_thick = thick
    ENDFOR
    temp_vol_cube_mask_xz[*,i,*] = mask_img
    temp_vol_thin_xz[*,i,*] = DILATE(thinImg, strucElem)
  ENDFOR

  ;destroy the progress bar
  progressbar -> Destroy
 
; ;GJ, 2020/06/24, calculate the center line
;  temp_vol_thin = temp_vol_thin_xy * temp_vol_thin_yz * temp_vol_thin_xz
;  SHADE_VOLUME, temp_vol_thin,0.9, v, p
;  smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
;  oPoly = OBJ_NEW('IDLgrPolygon', COLOR = [255, 0, 0], smoothedVertices , POLYGON=p)
;  mymodel = OBJ_NEW('IDLgrModel')
;  mymodel->Add, oPoly
;  xobjview_XDU, mymodel, BACKGROUND = [0, 0, 0]
  
  
  max_thick = MAX([MAX(temp_vol_cube_mask_xy), MAX(temp_vol_cube_mask_yz), MAX(temp_vol_cube_mask_xz)])
  print, 'max_thick = ', max_thick
  print, 'min_thick = ', min_thick
  temp_vol_cube_xy = BYTSCL(temp_vol_cube_mask_xy, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, temp_vol_cube_xy[*,*, mask_dim[2]/2]
  temp_vol_cube_yz = BYTSCL(temp_vol_cube_mask_yz, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, REFORM(temp_vol_cube_yz[mask_dim[0]/2,*,*])
  temp_vol_cube_xz = BYTSCL(temp_vol_cube_mask_xz, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, REFORM(temp_vol_cube_xz[*,mask_dim[1]/2,*])

  vert_colors = verts
  FOR i=0L, N_ELEMENTS(verts[0,*])-1L DO BEGIN
    vert_color_ind_xy = temp_vol_cube_xy[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind_yz = temp_vol_cube_yz[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind_xz = temp_vol_cube_xz[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind = MEAN([vert_color_ind_xy, vert_color_ind_xz, vert_color_ind_yz])
    ;Fractional Anisotropy
    L1 = vert_color_ind_xy*1.
    L2 = vert_color_ind_xz*1.
    L3 = vert_color_ind_yz*1.
    max_L = MAX([L1,L2,L3])
    min_L = MIN([L1,L2,L3])
    mid_L = L1+L2+L3-max_L-min_L
    IF L1+L2+L3 GT 0 THEN BEGIN
      eccentricity = ABS(max_L-mid_L)/(max_L+mid_L)     ;SQRT((L1-L2)^2 + (L2-L3)^2 +(L3-L1)^2)/2./SQRT(L1^2+L2^2+L3^2)
      vert_color_eccentricity = eccentricity*255.
    ENDIF ELSE BEGIN
      eccentricity = 0.
      vert_color_eccentricity = eccentricity*255.
    ENDELSE

    vert_colors[*,i] = TRANSPOSE(rgb_table[255-vert_color_ind, *])
    ;  vert_colors[*,i] = TRANSPOSE(rgb_table[vert_color_eccentricity, *])


    ;  vert_colors[0,i] = redvalue[vert_color_ind]
    ;  vert_colors[1,i] = greenvalue[vert_color_ind]
    ;  vert_colors[2,i] = bluevalue[vert_color_ind]
  ENDFOR

  RETURN, vert_colors
END

;pick vert colors based on a new matrix ;GJ 2019/2/13
;vert_colors = PICK_VERT_COLORS(vol_HU_cube, sel_one.verts)
;计算血管形态; calculate the shape and morphology of vessel
;modified to assign different vessel branches; GJ 2020/7/5
FUNCTION PICK_VERT_COLORS_by_hough_segmentation, vol_cube_mask, verts, conn, rgb_table

  IF N_ELEMENTS(verts) EQ 0 THEN RETURN, 0

  ;calculate the central point of original vessels
  ;remove single straight lines, only keep the cross points
  ves_ind = WHERE(vol_cube_mask EQ 1, count_ves)
  IF count_ves GT 10 THEN BEGIN
    ind3d = ARRAY_INDICES(vol_cube_mask, ves_ind)
    ves_center = [MEAN(ind3d[0,*]), MEAN(ind3d[1,*]), MEAN(ind3d[2,*])]
  ENDIF; ELSE BEGIN
  ;  ves_center = central_loc
  ;ENDELSE

  ;smoothvalue = 10000
  ;smo_verts = MESH_SMOOTH(verts, conn, ITERATIONS=smoothvalue, LAMBDA=0.01)
  ;numberVertices = MESH_DECIMATE(smo_verts, conn, $
  ;  decimatedConn, VERTICES = decimatedVerts, PERCENT_VERTICES = 20)

  ;find the iterations without further changes
  lambda_value = 0.5
  ;lambda_value = 0.01
  iterations_step = 1000
  iterations = 0L
  old_verts = verts
  old_vol_cube_mask = vol_cube_mask*0.
  old_vol_cube_mask[old_verts[0,*], old_verts[1,*], old_verts[2,*]] = 1
  temp_distance_ind_list=DBLARR(50, N_ELEMENTS(verts[0,*]))*0.
  i=0
  REPEAT BEGIN
    new_verts = MESH_SMOOTH(old_verts, conn, ITERATIONS=iterations_step, LAMBDA=lambda_value)
    IF i LT 50 THEN temp_distance_ind_list[i,*] = SQRT((new_verts[0,*]-verts[0,*])^2 + (new_verts[1,*]-verts[1,*])^2 + (new_verts[2,*]-verts[2,*])^2)

    new_vol_cube_mask = vol_cube_mask*0.
    new_vol_cube_mask[new_verts[0,*], new_verts[1,*], new_verts[2,*]] = 1

    mask_diff = new_vol_cube_mask - old_vol_cube_mask
    diff_new_index = WHERE(mask_diff EQ 1, diff_new_count)
    diff_old_index = WHERE(mask_diff EQ -1, diff_old_count)

    old_verts = new_verts
    old_vol_cube_mask = new_vol_cube_mask
    iterations = iterations + iterations_step
    print, 'iterations =', iterations
    print, 'diff_old_count =', diff_old_count, ', diff_new_count =', diff_new_count
    i++
  ENDREP UNTIL diff_new_count LT 50;10;50;200;50
  temp_distance_ind_list_size=i-1
  smo_verts = new_verts
;  ;plot the oOply0
;  oPoly0 = OBJ_NEW('IDLgrPolygon',COLOR = [255,255,255], new_verts, POLYGON=conn, STYLE=1)
;  mymodel0 = OBJ_NEW('IDLgrModel')
;  mymodel0->Add, oPoly0
;  xobjview_XDU, mymodel0, BACKGROUND = [0, 0, 0]
  
  ;GJ 2020/7/7 test the image cube
  ;iimage, REFORM(vol_cube_mask[*,127,*])
  
  new_verts_range = [MAX(new_verts[0,*])-MIN(new_verts[0,*]), MAX(new_verts[1,*])-MIN(new_verts[1,*]), MAX(new_verts[2,*])-MIN(new_verts[2,*])]
  new_verts_order = SORT(new_verts_range)
  ;only assume yz plane
  new_verts_order = [1, 2, 0]
  direc0 = MIN(new_verts_order[1:2])
  direc1 = MAX(new_verts_order[1:2])
  ;GJ, 2020/7/5
  size_vol = SIZE(new_vol_cube_mask, /dimensions)
  mip_slice2d = DBLARR(size_vol[direc0], size_vol[direc1])*0.
  IF new_verts_order[0] EQ 2 THEN BEGIN
    ;xy plane
    FOR i=0, size_vol[0]-1 DO BEGIN
      FOR j=0, size_vol[1]-1 DO BEGIN
        mip_slice2d[i, j] = MAX(new_vol_cube_mask[i, j, *])
      ENDFOR
    ENDFOR
    ;iimage, mip_slice2dxy
  ENDIF

  IF new_verts_order[0] EQ 1 THEN BEGIN
    ;xy plane
    FOR i=0, size_vol[0]-1 DO BEGIN
      FOR k=0, size_vol[2]-1 DO BEGIN
        mip_slice2d[i, k] = MAX(new_vol_cube_mask[i, *, k])
      ENDFOR
    ENDFOR
    ;iimage, mip_slice2dxz
  ENDIF
  
  IF new_verts_order[0] EQ 0 THEN BEGIN
    ;yz plane
    FOR j=0, size_vol[1]-1 DO BEGIN
      FOR k=0, size_vol[2]-1 DO BEGIN
        mip_slice2d[j, k] = MAX(new_vol_cube_mask[*, j, k])
      ENDFOR
    ENDFOR
    ;iimage, mip_slice2dyz
  ENDIF
  
  ;segment the vessle by hough transform
  vlines_hough_transformj, mip_slice2d, vline_info, back_maskend
  
  ;load the color table 13
  LOADCT, 13, RGB_TABLE = rgb_table
  ;analyze the vessel segment
  vert_colors = verts*0.
  ;temp_vline_ind=DBLARR(N_ELEMENTS(verts[0,*]))*0.
  back_maskend_vline = BYTSCL(back_maskend)
;  iimage, back_maskend_vline
  FOR m=0, N_ELEMENTS(vert_colors[1,*])-1 DO BEGIN
    IF new_verts_order[0] EQ 2 THEN temp_vline_ind = back_maskend_vline[new_verts[0,m], new_verts[1,m]]
    IF new_verts_order[0] EQ 1 THEN temp_vline_ind = back_maskend_vline[new_verts[0,m], new_verts[2,m]]
    IF new_verts_order[0] EQ 0 THEN temp_vline_ind = back_maskend_vline[new_verts[1,m], new_verts[2,m]]
    vert_colors[*,m] = TRANSPOSE(255-rgb_table[temp_vline_ind, *])
  ENDFOR
  
  RETURN, vert_colors
END


;pick vert colors based on a new matrix ;GJ 2019/2/13
;vert_colors = PICK_VERT_COLORS(vol_HU_cube, sel_one.verts)
;计算血管形态; calculate the shape and morphology of vessel
FUNCTION PICK_VERT_COLORS_by_thickness_smooth_old, vol_cube_mask, central_loc, verts, conn, rgb_table

IF N_ELEMENTS(verts) EQ 0 THEN RETURN, 0

;calculate the central point of original vessels
;remove single straight lines, only keep the cross points
ves_ind = WHERE(vol_cube_mask EQ 1, count_ves)
IF count_ves GT 10 THEN BEGIN
  ind3d = ARRAY_INDICES(vol_cube_mask, ves_ind)
  ves_center = [MEAN(ind3d[0,*]), MEAN(ind3d[1,*]), MEAN(ind3d[2,*])]
ENDIF ELSE BEGIN
  ves_center = central_loc
ENDELSE

;smoothvalue = 10000
;smo_verts = MESH_SMOOTH(verts, conn, ITERATIONS=smoothvalue, LAMBDA=0.01)
;numberVertices = MESH_DECIMATE(smo_verts, conn, $
;  decimatedConn, VERTICES = decimatedVerts, PERCENT_VERTICES = 20)

;find the iterations without further changes
lambda_value = 0.5
;lambda_value = 0.01
iterations_step = 1000
iterations = 0L
old_verts = verts
old_vol_cube_mask = vol_cube_mask*0.
old_vol_cube_mask[old_verts[0,*], old_verts[1,*], old_verts[2,*]] = 1
temp_distance_ind_list=DBLARR(50, N_ELEMENTS(verts[0,*]))*0.
i=0
REPEAT BEGIN
  new_verts = MESH_SMOOTH(old_verts, conn, ITERATIONS=iterations_step, LAMBDA=lambda_value)
  IF i LT 50 THEN temp_distance_ind_list[i,*] = SQRT((new_verts[0,*]-verts[0,*])^2 + (new_verts[1,*]-verts[1,*])^2 + (new_verts[2,*]-verts[2,*])^2)
  
  new_vol_cube_mask = vol_cube_mask*0.
  new_vol_cube_mask[new_verts[0,*], new_verts[1,*], new_verts[2,*]] = 1
  
  mask_diff = new_vol_cube_mask - old_vol_cube_mask
  diff_new_index = WHERE(mask_diff EQ 1, diff_new_count)
  diff_old_index = WHERE(mask_diff EQ -1, diff_old_count)
  
  old_verts = new_verts
  old_vol_cube_mask = new_vol_cube_mask
  iterations = iterations + iterations_step
  print, 'iterations =', iterations
  print, 'diff_old_count =', diff_old_count, ', diff_new_count =', diff_new_count
  i++
ENDREP UNTIL diff_new_count LT 50;10;50;200;50
temp_distance_ind_list_size=i-1
smo_verts = new_verts
;plot the oOply0
oPoly0 = OBJ_NEW('IDLgrPolygon',COLOR = [255,255,255], new_verts, POLYGON=conn, STYLE=1)
mymodel0 = OBJ_NEW('IDLgrModel')
mymodel0->Add, oPoly0
xobjview_XDU, mymodel0, BACKGROUND = [0, 0, 0]
;omodel = OBJ_NEW('IDLgrModel')
;omodel->Add, oPoly0
;oview = OBJ_NEW('IDLgrView',COLOR = [0,0,0])
;oview->Add, omodel
;owindow = OBJ_NEW('IDLgrWindow', dimension = [600, 600])
;owindow->Draw, oview
;odraw_img = owindow->READ()
;odraw_img->getProperty, DATA = draw_img_tmp
;; Write a PNG file to the temporary directory
;; Note the use of the TRUE keyword to TVRD
;filename = FILEPATH('test.png', ROOT_DIR='C:\D_drive\AIMIS_3D\images\')
;WRITE_PNG, filename, draw_img_tmp
;PRINT, 'File written to ', filename

;image thinning to the image
;do this to all three directions
mask_dim = SIZE(new_vol_cube_mask, /DIMENSIONS)
FOR k = 0, mask_dim[2]-1 DO BEGIN
  mask_img = image_thinning_img(new_vol_cube_mask[*,*,k])
  new_vol_cube_mask[*,*,k] = mask_img
ENDFOR
FOR j = 0, mask_dim[1]-1 DO BEGIN
  mask_img = image_thinning_img(REFORM(new_vol_cube_mask[*,j,*]))
  new_vol_cube_mask[*,j,*] = mask_img
ENDFOR
FOR i = 0, mask_dim[0]-1 DO BEGIN
  mask_img = image_thinning_img(REFORM(new_vol_cube_mask[i,*,*]))
  new_vol_cube_mask[i,*,*] = mask_img
ENDFOR

;remove single straight lines, only keep the cross points
old_vol_cube_mask = new_vol_cube_mask
smo_vessel_ind = WHERE(old_vol_cube_mask EQ 1, count_vessel)
IF count_vessel GT 10 THEN BEGIN
  ind3d = ARRAY_INDICES(old_vol_cube_mask, smo_vessel_ind)
  FOR i=0L, count_vessel-1L DO BEGIN
    small_mat = old_vol_cube_mask[(ind3d[0,i]-1):(ind3d[0,i]+1), (ind3d[1,i]-1):(ind3d[1,i]+1), (ind3d[2,i]-1):(ind3d[2,i]+1)]
    vind = WHERE(small_mat EQ 1, count)
    ;optimized 7 because 8 is too large
    IF count LT 6 THEN BEGIN
      old_vol_cube_mask[ind3d[0,i], ind3d[1,i], ind3d[2,i]]=0
      print, 'Less tan 6 i = ', i
    ENDIF ELSE BEGIN
      print, 'count = ', count
    ENDELSE
    
;    IF count GT 8 THEN BEGIN
;      old_vol_cube_mask[ind3d[0,i], ind3d[1,i], ind3d[2,i]]=0
;      print, 'Greater than 10 i = ', i
;    ENDIF
  ENDFOR
ENDIF

;dilate the cross points
radius = 1
strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
mask_dim = SIZE(old_vol_cube_mask, /DIMENSIONS)
FOR i = 0L, mask_dim[0]-1 DO BEGIN
  mask_img = DILATE(REFORM(old_vol_cube_mask[i,*,*]), strucElem)
  old_vol_cube_mask[i,*,*] = mask_img
ENDFOR
FOR j = 0L, mask_dim[1]-1 DO BEGIN
  mask_img = DILATE(REFORM(old_vol_cube_mask[*,j,*]), strucElem)
  old_vol_cube_mask[*,j,*] = mask_img
ENDFOR
FOR k = 0L, mask_dim[2]-1 DO BEGIN
  mask_img = DILATE(old_vol_cube_mask[*,*,k], strucElem)
  old_vol_cube_mask[*,*,k] = mask_img
ENDFOR
oPoly01 = OBJ_NEW('IDLgrPolygon',COLOR = [0,255,0], new_verts, POLYGON=conn, STYLE=1, alpha_channel=0.6)
SHADE_VOLUME, old_vol_cube_mask, 0.9, v, p
smoothedVertices = MESH_SMOOTH(v, p, ITERATIONS=200)
oPoly1 = OBJ_NEW('IDLgrPolygon',COLOR = [255,127,127], smoothedVertices, POLYGON=p, STYLE=1)
mymodel1 = OBJ_NEW('IDLgrModel')
mymodel1->Add, oPoly01
mymodel1->Add, oPoly1
xobjview_XDU, mymodel1, BACKGROUND = [0, 0, 0]

;pick the cross points
labelImg = LABEL_REGION(old_vol_cube_mask, $
  ALL_NEIGHBORS=allNeighbors, $
  ULONG=ulong)
labelmax=max(labelImg)
h = HISTOGRAM(labelImg, REVERSE_INDICES=r)
h[0]=0
print,'max', max(h)
print,'labelmax = ', labelmax
;the distance to x direction
x_dist_arr = DBLARR(labelmax)
ves_dist_arr = DBLARR(labelmax)
cross_ind3d_mean = DBLARR(3, labelmax)
;iplot, h
FOR i=1, labelmax DO BEGIN
  vessel_cross_ind = WHERE(labelImg EQ i, count_vessel_cross)
  cross_ind3d = ARRAY_INDICES(labelImg, vessel_cross_ind)
  ;print, i, ' region=: ', MEAN(cross_ind3d[0,*])-central_loc[0], ',', MEAN(cross_ind3d[1,*])-central_loc[1], ',',MEAN(cross_ind3d[2,*])-central_loc[2]
  print, i, ' region=: ', MEAN(cross_ind3d[0,*]), ',', MEAN(cross_ind3d[1,*]), ',',MEAN(cross_ind3d[2,*])
  ;record the mean of vessel cross points
  cross_ind3d_mean[*,i-1] = [MEAN(cross_ind3d[0,*]), MEAN(cross_ind3d[1,*]), MEAN(cross_ind3d[2,*])]
  ;record the distance to x direction
  x_dist_arr[i-1] = ABS(MEAN(cross_ind3d[0,*])-central_loc[0]);;;;;;ABS(MEAN(cross_ind3d[0,*])-ves_center[0]);
  ves_dist_arr[i-1] = SQRT((MEAN(cross_ind3d[0,*])-ves_center[0])^2 + (MEAN(cross_ind3d[1,*])-ves_center[1])^2 + (MEAN(cross_ind3d[2,*])-ves_center[2])^2)
ENDFOR


;find the shortest distance to x central pixel
sorted_x_dist_arr_ind = SORT(x_dist_arr)
print, 'sorted x distance array:', x_dist_arr[sorted_x_dist_arr_ind]
IF x_dist_arr[sorted_x_dist_arr_ind[1]] LT 2.*x_dist_arr[sorted_x_dist_arr_ind[0]] THEN BEGIN
  temp = [ves_dist_arr[sorted_x_dist_arr_ind[0]], ves_dist_arr[sorted_x_dist_arr_ind[1]], ves_dist_arr[sorted_x_dist_arr_ind[2]]]
  print, 'sorted vessel distance array:', temp
  min_temp = MIN(temp, min_ind)
  n = sorted_x_dist_arr_ind[min_ind]
ENDIF ELSE BEGIN
  n=sorted_x_dist_arr_ind[0]
ENDELSE
print, 'region index = ', n+1
print, 'point loc = ', cross_ind3d_mean[*,n]

;  ;check whether there are 2 to 3 branches
;  s=3
;  ;left + right , at least 1
;  temp_left = new_vol_cube_mask[cross_ind3d_mean[0,n]-s, cross_ind3d_mean[1,n]-s:cross_ind3d_mean[1,n]+s, cross_ind3d_mean[2,n]-s:cross_ind3d_mean[2,n]+s]
;  ind_temp_left = WHERE(temp_left GE 1, count_left)
;  temp_right = new_vol_cube_mask[cross_ind3d_mean[0,n]+s, cross_ind3d_mean[1,n]-s:cross_ind3d_mean[1,n]+s, cross_ind3d_mean[2,n]-s:cross_ind3d_mean[2,n]+s]
;  ind_temp_right = WHERE(temp_right GE 1, count_right)
;  ;top, at least 1
;  temp_top = new_vol_cube_mask[cross_ind3d_mean[0,n]-s:cross_ind3d_mean[0,n]+s, cross_ind3d_mean[1,n]-s:cross_ind3d_mean[1,n]+s, cross_ind3d_mean[2,n]+s]
;  ind_temp_top = WHERE(temp_top GE 1, count_top)
;  ;bottom must be 0
;  temp_bottom = new_vol_cube_mask[cross_ind3d_mean[0,n]-s:cross_ind3d_mean[0,n]+s, cross_ind3d_mean[1,n]-s:cross_ind3d_mean[1,n]+s, cross_ind3d_mean[2,n]-s]
;  ind_temp_bottom = WHERE(temp_bottom GE 1, count_bottom)
;  print, 'count l,r,t,b = ', count_left, count_right, count_top, count_bottom
  
  ;check this is the right cross-point
  ;2020/6/2
;  IF (count_left+count_right GT 0 AND (count_top GT 0 OR count_bottom EQ 0)) THEN BEGIN
    temp_matrix = new_vol_cube_mask[cross_ind3d_mean[0,n]-30:cross_ind3d_mean[0,n]+30, 0:cross_ind3d_mean[1,n]+20, *]
    temp_image = DBLARR(N_ELEMENTS(temp_matrix[0,*,0]), N_ELEMENTS(temp_matrix[0,0,*]))*0.
    FOR j=0, N_ELEMENTS(temp_image[*,0])-1 DO BEGIN
      FOR k=0, N_ELEMENTS(temp_image[0,*])-1 DO BEGIN
        temp_image[j,k]=MAX(temp_matrix[*,j,k])
      ENDFOR
    ENDFOR
    
    window, 2
    tvscl, temp_image
    
    
    DEVICE, DECOMPOSED = 0, RETAIN = 2
    LOADCT, 0
    dims_temp = SIZE(temp_image, /DIMENSIONS)
    bigger_temp_image = DBLARR(dims_temp[0]+100, dims_temp[1]+100)*0.
    bigger_temp_image[50:49+dims_temp[0],50:49+dims_temp[1]] = temp_image
    A_image = ROT(bigger_temp_image, 45);ATE(temp_image, 3)
    pixel_ind = WHERE(A_image GT 0, count)
    IF count GT 5 THEN BEGIN
      pixel_ind2d = ARRAY_INDICES(A_image, pixel_ind)
      width_0 = MAX(pixel_ind2d[0,*])-MIN(pixel_ind2d[0,*])+1
      width_1 = MAX(pixel_ind2d[1,*])-MIN(pixel_ind2d[1,*])+1
      diff = ABS(width_1 - width_0)
      IF width_1 GT width_0 THEN BEGIN
        B_image = A_image[MIN(pixel_ind2d[0,*])-2-diff/2:MAX(pixel_ind2d[0,*])+2+diff/2,MIN(pixel_ind2d[1,*])-2:MAX(pixel_ind2d[1,*])+2]
      ENDIF ELSE BEGIN
        B_image = A_image[MIN(pixel_ind2d[0,*])-2:MAX(pixel_ind2d[0,*])+2,MIN(pixel_ind2d[1,*])-2-diff/2:MAX(pixel_ind2d[1,*])+2+diff/2]
      ENDELSE
      dims = SIZE(B_image, /DIMENSIONS)
      print, 'B image size = ', dims
      WINDOW, 0, XSIZE = 3*dims[0], YSIZE = 3*dims[1], $
        TITLE='Brain Vessel Display'
      ;tvscl, temp_image, 0
      tvscl, MAX(B_image)-B_image, 1
    ENDIF
;     tvscl, ROTATE(temp_image, 1), 1
;  ENDIF
;ENDFOR


;SHADE_VOLUME, new_vol_cube_mask, 0.9, v0, p0
;smoothedVertices0 = MESH_SMOOTH(v0, p0, ITERATIONS=200)

;this is from mesh_smooth 2020/6/7
;smoothedVertices0 = new_verts
;p0 = conn
;oPoly0 = OBJ_NEW('IDLgrPolygon',COLOR = [0,255,0], smoothedVertices0, POLYGON=p0, STYLE=1, alpha_channel=0.6)



;define a color table
;redvalue = BYTSCL(RANDOMN(seed, 256))
;redvalue = replicate(245,256,1)
;greenvalue = BINDGEN(256)
;bluevalue = reverse(greenvalue)
LOADCT, 13, RGB_TABLE = rgb_table
;temp_distance_ind=REFORM(verts[0,*]*0.)
;temp_distance_ind=SQRT((smo_verts[0,*]-verts[0,*])^2 + (smo_verts[1,*]-verts[1,*])^2 + (smo_verts[2,*]-verts[2,*])^2)

;analyze the distance
temp_distance_ind=DBLARR(N_ELEMENTS(verts[0,*]))*0.
window,11
FOR m=0, N_ELEMENTS((new_verts[0,*]))-1 DO BEGIN
  IF m EQ 0 THEN BEGIN
    plot, temp_distance_ind_list[0:(temp_distance_ind_list_size-1),m]
  ENDIF
  IF (m MOD 2000) EQ 0 THEN BEGIN
    oplot, temp_distance_ind_list[0:(temp_distance_ind_list_size-1),m]
  ENDIF
  ;determine the distance
  xa=FLOOR((temp_distance_ind_list_size-1)/2)
  ya=temp_distance_ind_list[xa,m]
  xb=temp_distance_ind_list_size-1
  yb=temp_distance_ind_list[xb,m]
  y0=yb-xb*(yb-ya)/(xb-xa)
  IF y0 LE 0 THEN temp_distance_ind[m]=temp_distance_ind_list[0,m] ELSE temp_distance_ind[m]=y0
ENDFOR

;average the distance measurements
distance_ind=temp_distance_ind*0.
smo_vert_x = REFORM(smo_verts[0,*])
smo_vert_y = REFORM(smo_verts[1,*])
smo_vert_z = REFORM(smo_verts[2,*])
FOR i=0L, N_ELEMENTS(smo_verts[0,*])-1L DO BEGIN
  range = 0.5;1.;0.5
  index_list = WHERE(smo_vert_x GT smo_verts[0,i]-range AND smo_vert_x LT smo_verts[0,i]+range $
    AND smo_vert_y GT smo_verts[1,i]-range AND smo_vert_y LT smo_verts[1,i]+range $
    AND smo_vert_z GT smo_verts[2,i]-range AND smo_vert_z LT smo_verts[2,i]+range, count)
;  print, "count = ", count
  IF count GT 3 THEN BEGIN
;    distance =  STDEV(temp_distance_ind[index_list])/MEAN(temp_distance_ind[index_list])
    distance =  MEAN(temp_distance_ind[index_list])
    distance_ind[index_list] = distance
  ENDIF
ENDFOR


vert_color_ind = BYTSCL(distance_ind)

vert_colors = verts
FOR i=0L, N_ELEMENTS(verts[0,*])-1L DO BEGIN
  vert_colors[*,i] = TRANSPOSE(rgb_table[255-vert_color_ind[i], *])
ENDFOR

RETURN, vert_colors
END

;calculate the pixel with 0
;central_loc = CENTRAL_LOC_CAL(xd_sState.additionalDir[0])
FUNCTION CENTRAL_LOC_CAL, afn
  filearr=file_search(afn,'*.dcm',count=num)
  IF num EQ 0 THEN BEGIN
    RETURN, ROUND([(SIZE[0]-1.)/2.0,(SIZE[1]-1.)/2.0,(SIZE[2]-1.)/2.0])
  ENDIF
  image_read, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal
  size = SIZE(vol_HU_cube, /dimension)
;  center_loc = origin_loc + TRANSPOSE((direction[*,0]*(SIZE[0]-1.)/2.0 + direction[*,1]*(SIZE[1]-1.)/2. + direction[*,2]*(SIZE[2]-1.)/2.)) * vol_HU_cube_resolution
  center_loc = -REFORM(ROUND(INVERT(TRANSPOSE(direction))##origin_loc/vol_HU_cube_resolution))
;  IF MIN(center_loc) LT 0 OR MIN(size-center_loc) LT 0 THEN BEGIN
;    RETURN, ROUND([(SIZE[0]-1.)/2.0,(SIZE[1]-1.)/2.0,(SIZE[2]-1.)/2.0])
;  ENDIF
  return, center_loc
END

;pick vert colors based on vessel thickness ;GJ 2020/5/23
FUNCTION PICK_VERT_COLORS_by_vessel_thickness, verts, conn, rgb_table

  IF N_ELEMENTS(verts) EQ 0 THEN RETURN, 0
  
  ;smooth using many iterations
  smoothvalue = 1000
  smo_verts = MESH_SMOOTH(verts, conn, ITERATIONS=smoothvalue, /FIXED_EDGE_VERTICES)

  ;calculate distance between original verts and smoothed verts
  distance_ind=REFORM(verts[0,*]*0.)
  distance_ind=SQRT((smo_verts[0,*]-verts[0,*])^2 + (smo_verts[1,*]-verts[1,*])^2 + (smo_verts[2,*]-verts[2,*])^2)
  
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='yellow', Text='0%', Title='Calculating vessel thickness...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start
  
  ;calculate thickness by averaging the distance based on each individual pixel
  thickness_ind=REFORM(verts[0,*]*0.)
  smo_vert_x = REFORM(smo_verts[0,*])
  smo_vert_y = REFORM(smo_verts[1,*])
  smo_vert_z = REFORM(smo_verts[2,*])
  step_inc=2
  FOR i=FLOOR(MIN(smo_vert_x)), CEIL(MAX(smo_vert_x)), step_inc DO BEGIN
    count=(i-FLOOR(MIN(smo_vert_x))+1.)/FLOOR(MAX(smo_vert_x)-MIN(smo_vert_x))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    FOR j=FLOOR(MIN(smo_vert_y)), CEIL(MAX(smo_vert_y)), step_inc DO BEGIN
      FOR k=FLOOR(MIN(smo_vert_z)), CEIL(MAX(smo_vert_z)), step_inc DO BEGIN
        index_list = WHERE(smo_vert_x GE i AND smo_vert_x LT i+step_inc $
              AND smo_vert_y GE j AND smo_vert_y LT j+step_inc $
              AND smo_vert_z GE k AND smo_vert_z LT k+step_inc, count)
        IF count GT 0 THEN BEGIN
          ;average_thickness = STDDEV(distance_ind[index_list])/MEDIAN(distance_ind[index_list])
          average_thickness = MEDIAN(distance_ind[index_list])
          thickness_ind[index_list] = average_thickness
        ENDIF
      ENDFOR
    ENDFOR
  ENDFOR
  
  ;destroy the progress bar
  progressbar -> Destroy
  
  ;assigned calculated vessel thickness as verts color index
  vert_color_ind = BYTSCL(thickness_ind)

  ;assign color index using color table
  LOADCT, 13, RGB_TABLE = rgb_table
  vert_colors = verts
  FOR i=0L, N_ELEMENTS(verts[0,*])-1L DO BEGIN
    vert_colors[*,i] = TRANSPOSE(rgb_table[255-vert_color_ind[i], *])
  ENDFOR

  RETURN, vert_colors
END


;pick vert colors based on a new matrix ;GJ 2019/2/13
;vert_colors = PICK_VERT_COLORS(vol_HU_cube, sel_one.verts)
FUNCTION PICK_VERT_COLORS_by_eccentricity, vol_cube_mask, verts, rgb_table

  IF N_ELEMENTS(verts) EQ 0 THEN RETURN, 0

  ;define a color table
  ;redvalue = BYTSCL(RANDOMN(seed, 256))
  ;redvalue = replicate(245,256,1)
  ;greenvalue = BINDGEN(256)
  ;bluevalue = reverse(greenvalue)
  LOADCT, 13, RGB_TABLE = rgb_table

  mask_dim = SIZE(vol_cube_mask, /DIM)
  min_thick = 20

  radius = 2
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius

  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='pink', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start

  ;scale vol cube to 0-255
  temp_vol_cube_mask_xy = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[2]-1 DO BEGIN
    count=(i+1.)/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;i = mask_dim[2]/2
    mask_img = DILATE(temp_vol_cube_mask_xy[*,*,i], strucElem)
    result = LABEL_REGION(mask_img, /ALL_NEIGHBORS)
    FOR j = 1, MAX(result) DO BEGIN
      temp_mask_img = mask_img*0
      index = WHERE(result EQ j, count)
      IF count NE 0 THEN temp_mask_img[index] = 1
      thick = image_thinning(temp_mask_img)
      ;    print, 'thick = ', thick
      mask_img[index] = thick
      IF thick LT min_thick THEN min_thick = thick
    ENDFOR
    temp_vol_cube_mask_xy[*,*,i] = mask_img
  ENDFOR

  ;scale vol cube to 0-255
  temp_vol_cube_mask_yz = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[0]-1 DO BEGIN
    count=(i+1.+mask_dim[2])/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;i = mask_dim[2]/2
    mask_img = DILATE(REFORM(temp_vol_cube_mask_yz[i,*,*]), strucElem)
    result = LABEL_REGION(mask_img, /ALL_NEIGHBORS)
    FOR j = 1, MAX(result) DO BEGIN
      temp_mask_img = mask_img*0
      index = WHERE(result EQ j, count)
      IF count NE 0 THEN temp_mask_img[index] = 1
      thick = image_thinning(temp_mask_img)
      mask_img[index] = thick
      IF thick LT min_thick THEN min_thick = thick
    ENDFOR
    temp_vol_cube_mask_yz[i,*,*] = mask_img
  ENDFOR

  ;scale vol cube to 0-255
  temp_vol_cube_mask_xz = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[1]-1 DO BEGIN
    count=(i+1.+mask_dim[2]+mask_dim[1])/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    ;i = mask_dim[2]/2
    mask_img = DILATE(REFORM(temp_vol_cube_mask_xz[*,i,*]), strucElem)
    result = LABEL_REGION(mask_img, /ALL_NEIGHBORS)
    FOR j = 1, MAX(result) DO BEGIN
      temp_mask_img = mask_img*0
      index = WHERE(result EQ j, count)
      IF count NE 0 THEN temp_mask_img[index] = 1
      thick = image_thinning(temp_mask_img)
      mask_img[index] = thick
      IF thick LT min_thick THEN min_thick = thick
    ENDFOR
    temp_vol_cube_mask_xz[*,i,*] = mask_img
  ENDFOR

  ;destroy the progress bar
  progressbar -> Destroy

  max_thick = MAX([MAX(temp_vol_cube_mask_xy), MAX(temp_vol_cube_mask_yz), MAX(temp_vol_cube_mask_xz)])
  print, 'max_thick = ', max_thick
  print, 'min_thick = ', min_thick
  temp_vol_cube_xy = BYTSCL(temp_vol_cube_mask_xy, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, temp_vol_cube_xy[*,*, mask_dim[2]/2]
  temp_vol_cube_yz = BYTSCL(temp_vol_cube_mask_yz, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, REFORM(temp_vol_cube_yz[mask_dim[0]/2,*,*])
  temp_vol_cube_xz = BYTSCL(temp_vol_cube_mask_xz, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, REFORM(temp_vol_cube_xz[*,mask_dim[1]/2,*])

  vert_colors = verts
  FOR i=0L, N_ELEMENTS(verts[0,*])-1L DO BEGIN
    vert_color_ind_xy = temp_vol_cube_xy[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind_yz = temp_vol_cube_yz[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind_xz = temp_vol_cube_xz[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind = MEAN([vert_color_ind_xy, vert_color_ind_xz, vert_color_ind_yz])
    ;Fractional Anisotropy
    L1 = vert_color_ind_xy*1.
    L2 = vert_color_ind_xz*1.
    L3 = vert_color_ind_yz*1.
    max_L = MAX([L1,L2,L3])
    min_L = MIN([L1,L2,L3])
    mid_L = L1+L2+L3-max_L-min_L
    IF L1+L2+L3 GT 0 THEN BEGIN
      eccentricity = ABS(max_L-mid_L)/(max_L+mid_L)     ;SQRT((L1-L2)^2 + (L2-L3)^2 +(L3-L1)^2)/2./SQRT(L1^2+L2^2+L3^2)
      vert_color_eccentricity = eccentricity*255.
    ENDIF ELSE BEGIN
      eccentricity = 0.
      vert_color_eccentricity = eccentricity*255.
    ENDELSE

    ;vert_colors[*,i] = TRANSPOSE(rgb_table[255-vert_color_ind, *])
    vert_colors[*,i] = TRANSPOSE(rgb_table[vert_color_eccentricity, *])+100.


    ;  vert_colors[0,i] = redvalue[vert_color_ind]
    ;  vert_colors[1,i] = greenvalue[vert_color_ind]
    ;  vert_colors[2,i] = bluevalue[vert_color_ind]
  ENDFOR

  RETURN, vert_colors
END


;GJ 2019/3/2
PRO Add_new_by_region_grow_spine, xd_sState, vol_cube_mask
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start
  
  mask_dim = SIZE(vol_cube_mask, /DIM)
  ;define a temporary
  temp_vol_cube_mask_xy = vol_cube_mask*0.
  ;  atemp_vol_cube_mask_xy = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[2]-1 DO BEGIN
    count=(i+1.)/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

    mask_img = 1-vol_cube_mask[*,*,i]
    original_mask_img = 1- mask_img
    temp_mask_img = mask_img * 0
    index_mask_img = WHERE(mask_img EQ 0, ncount)
    IF ncount GT 0 THEN BEGIN
      index_2d = ARRAY_INDICES(mask_img, index_mask_img)

      labelImg = LABEL_REGION(mask_img, /ALL_NEIGHBORS)

      h = HISTOGRAM(labelImg)

      hsize=size(h)
      ;    print,'hszie',hsize[1]
      IF hsize[1] GT 2 THEN BEGIN
        for k=2,hsize[1]-1 do begin
          IF h[k] GT ncount/6 AND h[k] LT ncount THEN BEGIN
            pos=where(labelImg EQ k, npos)
            pos2da = ARRAY_INDICES(labelImg, pos)
            pos2d = [MEAN(pos2da[0,*]), MEAN(pos2da[1,*])]
            IF pos2d[0] GT MEAN(index_2d[0,*])-STDDEV(index_2d[0,*]) AND pos2d[0] LT MEAN(index_2d[0,*])+STDDEV(index_2d[0,*]) AND pos2d[1] GT MEAN(index_2d[1,*])-STDDEV(index_2d[1,*]) AND pos2d[1] LT MEAN(index_2d[1,*])+STDDEV(index_2d[1,*]) THEN BEGIN
              temp_mask_img[pos]=h[k]
              n_angles = 36
              spine_thickness = INTARR(n_angles)*0
              temp_spine_thickness = spine_thickness
              radius = 0.
              
              ;define spine orientation
              ;index_spine = WHERE(original_mask_img GE 1, ncount_spine)
              ;index2d_spine = ARRAY_INDICES(original_mask_img, index_spine)
              max_value = MAX(pos2da[0,*], max_ind)
              min_value = MIN(pos2da[0,*], min_ind)
              initial_angle =  ATAN(1.*(pos2da[1,max_ind]-pos2da[1,min_ind])/(pos2da[0,max_ind]-pos2da[0,min_ind]))
              REPEAT BEGIN
                radius = radius + 1
                spine_thickness = temp_spine_thickness
                temp_spine_thickness = spine_thickness * 0
                FOR l = 0, n_angles-1 DO FOR m = 0, radius-1 DO BEGIN
                  index_x = pos2d[0]+m * COS(l*10./180.*!PI + initial_angle)
                  index_y = pos2d[1]+m * SIN(l*10./180.*!PI + initial_angle)
                  IF index_x LT N_ELEMENTS(original_mask_img[*,0]) AND index_y LT N_ELEMENTS(original_mask_img[0,*]) THEN BEGIN
                    temp_spine_thickness[l] = temp_spine_thickness[l] + original_mask_img[index_x, index_y]
                  ENDIF
                ENDFOR
              ENDREP UNTIL ((MAX(temp_spine_thickness) GT 0) AND MAX(ABS(spine_thickness-temp_spine_thickness)) EQ 0)
              print, 'finish with spine orientation angle = ', initial_angle/!PI*180.
              print, 'vertebral pedicle thickness 1 = ', MIN(spine_thickness[19:21])
              print, 'vertebral pedicle thickness 2 = ', MIN(spine_thickness[33:35])
              print, 'radius = ', radius-1
            ENDIF
          ENDIF
        endfor
      ENDIF
      temp_vol_cube_mask_xy[*,*,i] = temp_mask_img
    ENDIF
  ENDFOR
   
  ;define a temporary
  temp_vol_cube_mask_xy_rot = vol_cube_mask*0.
  TRANSFORM_VOLUME, vol_cube_mask, vol_cube_mask_rot, Rotation=[-15, 0, 0]
  ;  atemp_vol_cube_mask_xy = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[2]-1 DO BEGIN
    count=(i+mask_dim[0]+1.)/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

    mask_img = 1-vol_cube_mask_rot[*,*,i]
    original_mask_img = 1- mask_img
    temp_mask_img = mask_img * 0
    index_mask_img = WHERE(mask_img EQ 0, ncount)
    IF ncount GT 0 THEN BEGIN
      index_2d = ARRAY_INDICES(mask_img, index_mask_img)

      labelImg = LABEL_REGION(mask_img, /ALL_NEIGHBORS)

      h = HISTOGRAM(labelImg)

      hsize=size(h)
      ;    print,'hszie',hsize[1]
      IF hsize[1] GT 2 THEN BEGIN
        for k=2,hsize[1]-1 do begin
          IF h[k] GT ncount/6 AND h[k] LT ncount THEN BEGIN
            pos=where(labelImg EQ k, npos)
            pos2da = ARRAY_INDICES(labelImg, pos)
            pos2d = [MEAN(pos2da[0,*]), MEAN(pos2da[1,*])]
            IF pos2d[0] GT MEAN(index_2d[0,*])-STDDEV(index_2d[0,*]) AND pos2d[0] LT MEAN(index_2d[0,*])+STDDEV(index_2d[0,*]) AND pos2d[1] GT MEAN(index_2d[1,*])-STDDEV(index_2d[1,*]) AND pos2d[1] LT MEAN(index_2d[1,*])+STDDEV(index_2d[1,*]) THEN BEGIN
              temp_mask_img[pos]=h[k]
              n_angles = 36
              spine_thickness = INTARR(n_angles)*0
              temp_spine_thickness = spine_thickness
              radius = 0.
              
              ;define spine orientation
              ;index_spine = WHERE(original_mask_img GE 1, ncount_spine)
              ;index2d_spine = ARRAY_INDICES(original_mask_img, index_spine)
              max_value = MAX(pos2da[0,*], max_ind)
              min_value = MIN(pos2da[0,*], min_ind)
              initial_angle =  ATAN(1.*(pos2da[1,max_ind]-pos2da[1,min_ind])/(pos2da[0,max_ind]-pos2da[0,min_ind]))
              REPEAT BEGIN
                radius = radius + 1
                spine_thickness = temp_spine_thickness
                temp_spine_thickness = spine_thickness * 0
                FOR l = 0, n_angles-1 DO FOR m = 0, radius-1 DO BEGIN
                  index_x = pos2d[0]+m * COS(l*10./180.*!PI + initial_angle)
                  index_y = pos2d[1]+m * SIN(l*10./180.*!PI + initial_angle)
                  IF index_x LT N_ELEMENTS(original_mask_img[*,0]) AND index_y LT N_ELEMENTS(original_mask_img[0,*]) THEN BEGIN
                    temp_spine_thickness[l] = temp_spine_thickness[l] + original_mask_img[index_x, index_y]
                  ENDIF
                ENDFOR
              ENDREP UNTIL ((MAX(temp_spine_thickness) GT 0) AND MAX(ABS(spine_thickness-temp_spine_thickness)) EQ 0)
              print, 'finish with spine orientation angle = ', initial_angle/!PI*180.
              print, 'vertebral pedicle thickness 1 = ', MIN(spine_thickness[19:21])
              print, 'vertebral pedicle thickness 2 = ', MIN(spine_thickness[33:35])
              print, 'radius = ', radius-1
            ENDIF
          ENDIF
        endfor
      ENDIF
      
      temp_vol_cube_mask_xy_rot[*,*,i] = temp_mask_img
    ENDIF
  ENDFOR
  TRANSFORM_VOLUME, temp_vol_cube_mask_xy_rot, new_temp_vol_cube_mask_xy_rot, Rotation=[15, 0, 0]
 
  ;define a temporary
  temp_vol_cube_mask_xy_rot2 = vol_cube_mask*0.
  TRANSFORM_VOLUME, vol_cube_mask, vol_cube_mask_rot2, Rotation=[-30, 0, 0]
  ;  atemp_vol_cube_mask_xy = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[2]-1 DO BEGIN
    count=(i+mask_dim[0]*2+1.)/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'

    mask_img = 1-vol_cube_mask_rot2[*,*,i]
    original_mask_img = 1- mask_img
    temp_mask_img = mask_img * 0
    index_mask_img = WHERE(mask_img EQ 0, ncount)
    IF ncount GT 0 THEN BEGIN
      index_2d = ARRAY_INDICES(mask_img, index_mask_img)

      labelImg = LABEL_REGION(mask_img, /ALL_NEIGHBORS)

      h = HISTOGRAM(labelImg)

      hsize=size(h)
      ;    print,'hszie',hsize[1]
      IF hsize[1] GT 2 THEN BEGIN
        for k=2,hsize[1]-1 do begin
          IF h[k] GT ncount/6 AND h[k] LT ncount THEN BEGIN
            pos=where(labelImg EQ k, npos)
            pos2da = ARRAY_INDICES(labelImg, pos)
            pos2d = [MEAN(pos2da[0,*]), MEAN(pos2da[1,*])]
            IF pos2d[0] GT MEAN(index_2d[0,*])-STDDEV(index_2d[0,*]) AND pos2d[0] LT MEAN(index_2d[0,*])+STDDEV(index_2d[0,*]) AND pos2d[1] GT MEAN(index_2d[1,*])-STDDEV(index_2d[1,*]) AND pos2d[1] LT MEAN(index_2d[1,*])+STDDEV(index_2d[1,*]) THEN BEGIN
              temp_mask_img[pos]=h[k]
              n_angles = 36
              spine_thickness = INTARR(n_angles)*0
              temp_spine_thickness = spine_thickness
              radius = 0.
              
              ;define spine orientation
              ;index_spine = WHERE(original_mask_img GE 1, ncount_spine)
              ;index2d_spine = ARRAY_INDICES(original_mask_img, index_spine)
              max_value = MAX(pos2da[0,*], max_ind)
              min_value = MIN(pos2da[0,*], min_ind)
              initial_angle =  ATAN(1.*(pos2da[1,max_ind]-pos2da[1,min_ind])/(pos2da[0,max_ind]-pos2da[0,min_ind]))
              REPEAT BEGIN
                radius = radius + 1
                spine_thickness = temp_spine_thickness
                temp_spine_thickness = spine_thickness * 0
                FOR l = 0, n_angles-1 DO FOR m = 0, radius-1 DO BEGIN
                  index_x = pos2d[0]+m * COS(l*10./180.*!PI + initial_angle)
                  index_y = pos2d[1]+m * SIN(l*10./180.*!PI + initial_angle)
                  IF index_x LT N_ELEMENTS(original_mask_img[*,0]) AND index_y LT N_ELEMENTS(original_mask_img[0,*]) THEN BEGIN
                    temp_spine_thickness[l] = temp_spine_thickness[l] + original_mask_img[index_x, index_y]
                  ENDIF
                ENDFOR
              ENDREP UNTIL ((MAX(temp_spine_thickness) GT 0) AND MAX(ABS(spine_thickness-temp_spine_thickness)) EQ 0)
              print, 'finish with spine orientation angle = ', initial_angle/!PI*180.
              print, 'vertebral pedicle thickness 1 = ', MIN(spine_thickness[19:21])
              print, 'vertebral pedicle thickness 2 = ', MIN(spine_thickness[33:35])
              print, 'radius = ', radius-1
            ENDIF
          ENDIF
        endfor
      ENDIF

      temp_vol_cube_mask_xy_rot2[*,*,i] = temp_mask_img
    ENDIF
  ENDFOR
  TRANSFORM_VOLUME, temp_vol_cube_mask_xy_rot2, new_temp_vol_cube_mask_xy_rot2, Rotation=[30, 0, 0]

  
  ;destroy the progress bar
  progressbar -> Destroy
  
  vol_HU_cube_mask = (temp_vol_cube_mask_xy + new_temp_vol_cube_mask_xy_rot+new_temp_vol_cube_mask_xy_rot2) GT 5
  SHADE_VOLUME, vol_HU_cube_mask, 0.9, Outverts, Outconn1
  ;  SHADE_VOLUME, mask_vol_cube,threshold_min_value, Outverts, Outconn1

  IF N_ELEMENTS(Outverts) GT 10 THEN BEGIN
    smoothedOutverts = MESH_SMOOTH(Outverts, Outconn1)
    verts = smoothedOutverts
    xd_sState.oPoly_size = [MAX(verts[0,*])-MIN(verts[0,*]), MAX(verts[1,*])-MIN(verts[1,*]), MAX(verts[2,*])-MIN(verts[2,*])]
    xd_sState.oPoly_min = [MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])]
    xd_sState.oPoly_max = [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]
    verts_border = [[MIN(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MIN(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MIN(verts[2,*])], $
      [MIN(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MIN(verts[1,*]), MAX(verts[2,*])], $
      [MAX(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])], $
      [MIN(verts[0,*]), MAX(verts[1,*]), MAX(verts[2,*])]]
    xd_sState.oPoly_border->setProperty, DATA = verts_border

    ;define to the rotation center
    ma=max(smoothedOutverts)
    mi=min(smoothedOutverts)
    bias=(ma-mi)/2
    rl=mi+bias
    rr=ma+bias
    cc=[(-rl)/(rr-rl),1/(rr-rl)]
    xma = max((smoothedOutverts)[0,*])
    xmi = min((smoothedOutverts)[0,*])
    xmid = 0.5*(xma+xmi)*cc[1]
    yma = max((smoothedOutverts)[1,*])
    ymi = min((smoothedOutverts)[1,*])
    ymid = 0.5*(yma+ymi)*cc[1]
    zma = max((smoothedOutverts)[2,*])
    zmi = min((smoothedOutverts)[2,*])
    zmid = 0.5*(zma+zmi)*cc[1]

    ;Update modify_Step
    modify_Step = xd_sState.struc_modifyobject.current_step_index +1
    xd_sState.struc_modifyobject.modify_step = modify_Step
    xd_sState.struc_modifyobject.current_step_index = modify_Step
    widget_control, xd_sState.wUndo, SENSITIVE = 1

    ;update exist number
    exist_num = xd_sState.struc_modifyobject.exist_num[modify_Step-1]
    xd_sState.struc_modifyobject.exist_num[modify_Step] = exist_num+1
    n_Poly = xd_sState.struc_modifyobject.exist_num[modify_Step]

    ObjectIndex_value = LONG(n_Poly)
    print, 'get_value = ', ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex_text, SET_VALUE=STRING(ObjectIndex_value, format='(I3)')
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_SLIDER_MAX=ObjectIndex_value
    WIDGET_CONTROL, xd_sState.wObjectIndex, SET_VALUE=ObjectIndex_value

    ;    xd_sState.struc_modifyobject.hide_status[modify_Step+1, 0:(n_Poly-1)] = 1
    xd_sState.struc_modifyobject.hide_status[modify_step, *] = xd_sState.struc_modifyobject.hide_status[modify_step-1, *]
    xd_sState.struc_modifyobject.hide_status[modify_Step, n_Poly-1] = 0

    ;define a new polygon with color
    color_rgb = COLORTABLE_rgb_3D(colorIndex, ColorIndexBGR)
    vert_colors = smoothedOutverts
    vert_colors[0,*] = color_rgb[n_Poly, 0]
    vert_colors[1,*] = color_rgb[n_Poly, 1]
    vert_colors[2,*] = color_rgb[n_Poly, 2]

    ;define mask_vol_cube
    index = WHERE(vol_HU_cube_mask GT 0.9, count)
    IF count GT 1 THEN BEGIN
      xd_sState.struc_modifyobject.p_maskindex_list[n_Poly-1] = PTR_NEW(index)
      xd_sState.struc_modifyobject.countmask_list[n_Poly-1]=count

      xd_sState.tvol_white.SetProperty,hide=0
      xd_sState.vol_white = count * (xd_sState.vol_HU_cube_resolution^3)/1000. ;in cm3
      Stringwhite = STRING(xd_sState.vol_white,FORMAT='("Volume =", F9.3, " cm3")')
      ;      IF TOTAL(vert_colors[*,0]) LT 128*3. THEN fill_color=[255,255,255] ELSE fill_color=[0,0,0]
      IF TOTAL(vert_colors[*,0]) GT 1 THEN xd_sState.tvol_white.SetProperty,strings=stringwhite, color=vert_colors[*,0];, FILL_COLOR=fill_color
    ENDIF

    alpha_channel = 1.0
    oPoly = OBJ_NEW('IDLgrPolygon', smoothedOutverts,  NAME=STRING(n_Poly-1), POLYGON=Outconn1, STYLE=2, $
      SHADING=1, vert_colors=vert_colors, ALPHA_CHANNEL=alpha_channel)
    xd_sState.struc_modifyobject.polygon_list[n_Poly-1] = PTR_NEW(oPoly)
    xd_sState.struc_modifyobject.p_ori_verts_list[n_Poly-1] = PTR_NEW(smoothedOutverts)
    xd_sState.oPoly[n_Poly-1] = *(xd_sState.struc_modifyobject.polygon_list[n_Poly-1])
    xd_sState.selected = xd_sState.oPoly[n_Poly-1]
    xd_sState.oRotationModel->Add, oPoly
    FOR i = 0, n_Poly-1 DO BEGIN
      xd_sState.opoly[i]->getProperty, data=verts
      IF N_ELEMENTS(verts) EQ 0 THEN RETURN
      IF i EQ 0 THEN whole_verts = verts ELSE whole_verts = [[whole_verts], [verts]]
    ENDFOR
    ;center the selected object
    Adjust_whole_polygons, whole_verts, xd_sState



    dataxyz = REFORM(smoothedOutverts[*, N_ELEMENTS(verts(0,*))/2])
    ;GJ 2019/2/22, summarize all plot pickdata into a program
    PLOT_pickdata, dataxyz, xd_sState


    ;set the alpha-channel transparent value
    WIDGET_CONTROL, xd_sState.wtransslider, SET_VALUE=alpha_channel*100
    WIDGET_CONTROL, xd_sState.wTrans_text, SET_VALUE=STRING(alpha_channel*100, format='(f5.1)')

    ;Enable color selection
    WIDGET_CONTROL, xd_sState.wColorIndex, sensitive=1
    WIDGET_CONTROL, xd_sState.wColorIndex, SET_VALUE=n_Poly
    WIDGET_CONTROL, xd_sState.wColorIndex_draw, GET_VALUE=wColorIndexID
    DEVICE, DECOMPOSED=1
    WSET, wColorIndexID
    ERASE, ColorIndexBGR[n_Poly]

    ;enable these options
    ;    widget_control, xd_sState.wColorIndex, sensitive=1
    widget_control, xd_sState.wSmooth, sensitive=1
    widget_control, xd_sState.wSmooth_text, sensitive=1
    widget_control, xd_sState.wtransSlider, sensitive=1
    widget_control, xd_sState.wTrans_text, sensitive=1
    widget_control, xd_sState.w3DSegd, sensitive=1
    widget_control, xd_sState.w3Drg, sensitive=1
    widget_control, xd_sState.wthickness, sensitive=1
    widget_control, xd_sState.wpoint, sensitive=1
    widget_control, xd_sState.wwire, sensitive=1
    widget_control, xd_sState.wfill, sensitive=1
    widget_control, xd_sState.wexpstl, sensitive=1
    ;            widget_control, xd_sState.whide, sensitive=0
    widget_control, xd_sState.wsinglecenter, sensitive=1

    widget_control, xd_sState.wUndo, SENSITIVE = 1
    ;    widget_control, xd_sState.wRedo, SENSITIVE = 1

  ENDIF ELSE BEGIN
    infowarning = DIALOG_MESSAGE('No volume selection! Please reselect!', /ERROR)
  ENDELSE

END

;pick vert colors based on a new matrix ;GJ 2019/2/13
;vert_colors = PICK_VERT_COLORS(vol_HU_cube, sel_one.verts)
FUNCTION PICK_VERT_COLORS_by_region_grow, vol_cube_mask, verts, rgb_table

  IF N_ELEMENTS(verts) EQ 0 THEN RETURN, 0

  ;define a color table
  ;redvalue = BYTSCL(RANDOMN(seed, 256))
  ;redvalue = replicate(245,256,1)
  ;greenvalue = BINDGEN(256)
  ;bluevalue = reverse(greenvalue)
  LOADCT, 13, RGB_TABLE = rgb_table

  mask_dim = SIZE(vol_cube_mask, /DIM)
  radius = 2
  strucElem = SHIFT(DIST(2*radius+1), radius, radius) LE radius
  
  ;      ; Create the progress bar.
  progressbar = Obj_New('progressbar', Color='green', Text='0%', Title='Calculating...', /NOCANCEL, XOFFSET=100, XSIZE=345, YOFFSET=45)
  ;      ; Place the progress bar on the display.
  progressbar -> Start

  ;scale vol cube to 0-255
  temp_vol_cube_mask_xy = vol_cube_mask*0.
;  atemp_vol_cube_mask_xy = BYTSCL(vol_cube_mask)
  FOR i = 0, mask_dim[2]-1 DO BEGIN
    count=(i+1.)/(TOTAL(mask_dim))*100.0
    progressbar -> Update, count, Text=STRTRIM(STRING(FLOOR(count)),2) + '%'
    
    mask_img = 1-vol_cube_mask[*,*,i]
    index_mask_img = WHERE(mask_img EQ 0, ncount)
    IF ncount GT 0 THEN BEGIN
      index_2d = ARRAY_INDICES(mask_img, index_mask_img)

      labelImg = LABEL_REGION(mask_img, /ALL_NEIGHBORS)

      h = HISTOGRAM(labelImg)

      hsize=size(h)
      ;    print,'hszie',hsize[1]
      for k=1,hsize[1]-1 do begin
        IF h[k] GT ncount/4 AND h[k] LT ncount THEN BEGIN
          pos=where(labelImg EQ k, npos)
          pos2d = ARRAY_INDICES(labelImg, pos[npos/2])
          IF pos2d[0] GT MIN(index_2d[0,*]) AND pos2d[0] LT MAX(index_2d[0,*]) AND pos[1] GT MAX(index_2d[0,*]) AND pos2d[1] LT MAX(index_2d[1,*]) THEN BEGIN
            labelImg[pos]=h[k]
          ENDIF
        ENDIF
      endfor

      temp_vol_cube_mask_xy[*,*,i] = labelImg
    ENDIF
  ENDFOR

  ;destroy the progress bar
  progressbar -> Destroy

  max_thick = MAX(temp_vol_cube_mask_xy)
  min_thick = MIN(temp_vol_cube_mask_xy)
;  print, 'max_thick = ', max_thick
;  print, 'min_thick = ', min_thick
  temp_vol_cube_xy = BYTSCL(temp_vol_cube_mask_xy, MIN=min_thick, MAX=max_thick, TOP=255)
  ;iimage, temp_vol_cube_xy[*,*, mask_dim[2]/2]

  vert_colors = verts
  FOR i=0L, N_ELEMENTS(verts[0,*])-1L DO BEGIN
    vert_color_ind_xy = temp_vol_cube_xy[FLOOR(verts[0,i]), FLOOR(verts[1,i]), FLOOR(verts[2,i])]
    vert_color_ind = vert_color_ind_xy
    vert_colors[*,i] = TRANSPOSE(rgb_table[255-vert_color_ind, *])
    ;  vert_colors[0,i] = redvalue[vert_color_ind]
    ;  vert_colors[1,i] = greenvalue[vert_color_ind]
    ;  vert_colors[2,i] = bluevalue[vert_color_ind]
  ENDFOR

  RETURN, vert_colors
END

FUNCTION image_thinning_img, Img

  binaryImg = img GE 1
  h0 = [[0b,0,0], $
    [0,1,0], $
    [1,1,1]]
  m0 = [[1b,1,1], $
    [0,0,0], $
    [0,0,0]]
  h1 = [[0b,0,0], $
    [1,1,0], $
    [1,1,0]]
  m1 = [[0b,1,1], $
    [0,0,1], $
    [0,0,0]]
  h2 = [[1b,0,0], $
    [1,1,0], $
    [1,0,0]]
  m2 = [[0b,0,1], $
    [0,0,1], $
    [0,0,1]]
  h3 = [[1b,1,0], $
    [1,1,0], $
    [0,0,0]]
  m3 = [[0b,0,0], $
    [0,0,1], $
    [0,1,1]]
  h4 = [[1b,1,1], $
    [0,1,0], $
    [0,0,0]]
  m4 = [[0b,0,0], $
    [0,0,0], $
    [1,1,1]]
  h5 = [[0b,1,1], $
    [0,1,1], $
    [0,0,0]]
  m5 = [[0b,0,0], $
    [1,0,0], $
    [1,1,0]]
  h6 = [[0b,0,1], $
    [0,1,1], $
    [0,0,1]]
  m6 = [[1b,0,0], $
    [1,0,0], $
    [1,0,0]]
  h7 = [[0b,0,0], $
    [0,1,1], $
    [0,1,1]]
  m7 = [[1b,1,0], $
    [1,0,0], $
    [0,0,0]]

  bCont = 1b
  iIter = 0
  thinImg = binaryImg

  WHILE bCont EQ 1b DO BEGIN & $
    ;    PRINT,'Iteration: ', iIter & $
    inputImg = thinImg & $
    thinImg = MORPH_THIN(inputImg, h0, m0) & $
    thinImg = MORPH_THIN(thinImg, h1, m1) & $
    thinImg = MORPH_THIN(thinImg, h2, m2) & $
    thinImg = MORPH_THIN(thinImg, h3, m3) & $
    thinImg = MORPH_THIN(thinImg, h4, m4) & $
    thinImg = MORPH_THIN(thinImg, h5, m5) & $
    thinImg = MORPH_THIN(thinImg, h6, m6) & $
    thinImg = MORPH_THIN(thinImg, h7, m7) & $
    ;    TVSCL, thinImg, 2 & $
    ;    WAIT, 1 & $
    bCont = MAX(inputImg - thinImg) & $
    iIter = iIter + 1 & $
  ENDWHILE

RETURN, thinImg
END

FUNCTION image_thinning, Img, thinImg=thinImg

  binaryImg = img GE 1
  h0 = [[0b,0,0], $
    [0,1,0], $
    [1,1,1]]
  m0 = [[1b,1,1], $
    [0,0,0], $
    [0,0,0]]
  h1 = [[0b,0,0], $
    [1,1,0], $
    [1,1,0]]
  m1 = [[0b,1,1], $
    [0,0,1], $
    [0,0,0]]
  h2 = [[1b,0,0], $
    [1,1,0], $
    [1,0,0]]
  m2 = [[0b,0,1], $
    [0,0,1], $
    [0,0,1]]
  h3 = [[1b,1,0], $
    [1,1,0], $
    [0,0,0]]
  m3 = [[0b,0,0], $
    [0,0,1], $
    [0,1,1]]
  h4 = [[1b,1,1], $
    [0,1,0], $
    [0,0,0]]
  m4 = [[0b,0,0], $
    [0,0,0], $
    [1,1,1]]
  h5 = [[0b,1,1], $
    [0,1,1], $
    [0,0,0]]
  m5 = [[0b,0,0], $
    [1,0,0], $
    [1,1,0]]
  h6 = [[0b,0,1], $
    [0,1,1], $
    [0,0,1]]
  m6 = [[1b,0,0], $
    [1,0,0], $
    [1,0,0]]
  h7 = [[0b,0,0], $
    [0,1,1], $
    [0,1,1]]
  m7 = [[1b,1,0], $
    [1,0,0], $
    [0,0,0]]
    
  bCont = 1b
  iIter = 0
  thinImg = binaryImg
  
  WHILE bCont EQ 1b DO BEGIN & $
;    PRINT,'Iteration: ', iIter & $
    inputImg = thinImg & $
    thinImg = MORPH_THIN(inputImg, h0, m0) & $
    thinImg = MORPH_THIN(thinImg, h1, m1) & $
    thinImg = MORPH_THIN(thinImg, h2, m2) & $
    thinImg = MORPH_THIN(thinImg, h3, m3) & $
    thinImg = MORPH_THIN(thinImg, h4, m4) & $
    thinImg = MORPH_THIN(thinImg, h5, m5) & $
    thinImg = MORPH_THIN(thinImg, h6, m6) & $
    thinImg = MORPH_THIN(thinImg, h7, m7) & $
;    TVSCL, thinImg, 2 & $
;    WAIT, 1 & $
    bCont = MAX(inputImg - thinImg) & $
    iIter = iIter + 1 & $
  ENDWHILE

RETURN, iIter
END

;cut corner for naked eye 3D view
FUNCTION CUT_CORNER, draw_img
  draw_img = REVERSE(draw_img, 3)
  draw_img_size = SIZE(draw_img, /dim)
  FOR i=0, draw_img_size[1]/2-1 DO BEGIN
    draw_img[*, 0:(draw_img_size[1]/2-1-i), draw_img_size[1]-1-i] = 0.
    draw_img[*, (draw_img_size[1]/2+i):(draw_img_size[1]-1), draw_img_size[1]-1-i] = 0.
  ENDFOR

new_img = REVERSE(draw_img, 3)
RETURN, new_img
END

;add a point to the 3D object
;2019/2/22
PRO PLOT_pickdata, dataxyz, xd_sState
   statusString = STRING(dataxyz[0], $
    dataxyz[1],dataxyz[2], $
    FORMAT='("X=", F6.2,' + $
    ' ", Y=",F6.2,", Z=",F6.2)')
  demo_putTips, xd_sState, 'locat', 11, /LABEL
  demo_putTips, xd_sState, statusString, 12

  location3D = DBLARR(3)
  ;GJ 2019/5/21, fix the bug of mask display after using xd_xroi
  location3D[0] = (dataxyz[0]/((SIZE(xd_sState.vol_HU_cube))[1]))*xd_sState.imageSize + 1
  location3D[1] = (dataxyz[1]/((SIZE(xd_sState.vol_HU_cube))[2]))*xd_sState.imageSize + 1
  location3D[2] = (dataxyz[2]/((SIZE(xd_sState.vol_HU_cube))[3]))*xd_sState.imageSize + 1
  IF location3D[0] GT xd_sState.imageSize THEN location3D[0] = xd_sState.imageSize
  IF location3D[1] GT xd_sState.imageSize THEN location3D[1] = xd_sState.imageSize
  IF location3D[2] GT xd_sState.imageSize THEN location3D[2] = xd_sState.imageSize
  xd_sState.location3D = location3D

  WIDGET_CONTROL, xd_sState.DrawSlider[0], SET_VALUE = location3D[2]
  WIDGET_CONTROL, xd_sState.DrawSlider[1], SET_VALUE = location3D[0]
  WIDGET_CONTROL, xd_sState.DrawSlider[2], SET_VALUE = location3D[1]

  WIDGET_CONTROL, xd_sState.wScalingSlider, GET_VALUE=scale
  scale = 0.75 + FLOAT(scale) / 100.0

  s1 = OBJ_NEW('IDLgrPolyline', $
    [dataxyz[0]+10,dataxyz[0]-10], [dataxyz[1],dataxyz[1]], [dataxyz[2],dataxyz[2]], COLOR=[255,0,0])
  s2 = OBJ_NEW('IDLgrPolyline', $
    [dataxyz[0],dataxyz[0]], [dataxyz[1]+10, dataxyz[1]-10], [dataxyz[2],dataxyz[2]], COLOR=[0,255,0])
  s3 = OBJ_NEW('IDLgrPolyline', $
    [dataxyz[0],dataxyz[0]], [dataxyz[1],dataxyz[1]], [dataxyz[2]+10, dataxyz[2]-10] , COLOR=[0,0,255])

  s1->getproperty, POLYLINES=conn1
  s1->getproperty, DATA=verts1

  s2->getproperty, POLYLINES=conn2
  s2->getproperty, DATA=verts2

  s3->getproperty, POLYLINES=conn3
  s3->getproperty, DATA=verts3

  xd_sState.ogreenaxis1->setProperty, DATA = verts1
  xd_sState.ogreenaxis1->setProperty, POLYLINES = conn1
  xd_sState.ogreenaxis1->SetProperty, HIDE=0

  xd_sState.oyellowaxis1->setProperty, DATA = verts2
  xd_sState.oyellowaxis1->setProperty, POLYLINES = conn2
  xd_sState.oyellowaxis1->SetProperty, HIDE=0

  xd_sState.oBlueAxis1->setProperty, DATA = verts3
  xd_sState.oBlueAxis1->setProperty, POLYLINES = conn3
  xd_sState.oBlueAxis1->SetProperty, HIDE=0

  s = obj_new('orb',color=[255,0,0],radius=3./scale, POS = dataxyz)
  s->getproperty, POLYGONS=conn
  s->getproperty, DATA=verts

  xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
  xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
  xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2

  IF alpha0 EQ 0 THEN BEGIN
    xd_sState.oS0->setProperty, ALPHA_CHANNEL = 1, hide=0
  ENDIF ELSE BEGIN
    IF alpha1 EQ 0 THEN BEGIN
      xd_sState.oS1->setProperty, ALPHA_CHANNEL = 1, hide=0
    ENDIF ELSE BEGIN
      IF alpha2 EQ 0 THEN BEGIN
        xd_sState.oS2->setProperty, ALPHA_CHANNEL = 1, hide=0
      ENDIF
    ENDELSE
  ENDELSE
  IF xd_sState.dataXYZind EQ 0 THEN BEGIN
    xd_sState.oS0->setProperty, DATA = verts
    xd_sState.oS0->setProperty, POLYGONS = conn
    xd_sState.dataXYZ0 = dataxyz
    xd_sState.dataXYZind = 1
    xd_sState.oS0->setProperty, hide=0
  ENDIF ELSE BEGIN
    IF xd_sState.dataXYZind EQ 1 THEN BEGIN
      xd_sState.oS1->setProperty, DATA = verts
      xd_sState.oS1->setProperty, POLYGONS = conn
      xd_sState.dataXYZ1 = dataxyz
      xd_sState.dataXYZind = 2
      xd_sState.oS1->setProperty, hide=0
    ENDIF ELSE BEGIN
      xd_sState.oS2->setProperty, POLYGONS = conn
      xd_sState.oS2->setProperty, DATA = verts
      xd_sState.dataXYZ2 = dataxyz
      xd_sState.dataXYZind = 0
      xd_sState.oS2->setProperty, hide=0
    ENDELSE
  ENDELSE
  xd_sState.oS0->getProperty, ALPHA_CHANNEL = alpha0
  xd_sState.oS1->getProperty, ALPHA_CHANNEL = alpha1
  xd_sState.oS2->getProperty, ALPHA_CHANNEL = alpha2
  demo_putTips, xd_sState, 'locat', 11, /LABEL
  demo_putTips, xd_sState, statusString, 12
  datasum0=xd_sState.dataXYZ0[0]+xd_sState.dataXYZ0[1]+xd_sState.dataXYZ0[2]
  datasum1=xd_sState.dataXYZ1[0]+xd_sState.dataXYZ1[1]+xd_sState.dataXYZ1[2]
  datasum2=xd_sState.dataXYZ2[0]+xd_sState.dataXYZ2[1]+xd_sState.dataXYZ2[2]
  if (datasum0 NE 0) AND (datasum1 NE 0) and (datasum2 eq 0) then begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;hjx05252
    xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1, hide=0
    xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
    xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
    xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

    statusString = STRING(xd_sState.dataXYZdistance, $
      FORMAT='("Distance = ", F6.2, " mm")')
    demo_putTips, xd_sState, 'locat', 11, /LABEL
    demo_putTips, xd_sState, statusString, 12

    xd_sState.ckr1->setProperty, STRINGS= statusString,hide=0
  endif else begin
    if (datasum0 NE 0) AND (datasum1 NE 0) AND (datasum2 NE 0) THEN BEGIN
      if(xd_sState.COUNT123  eq 0) THEN begin
        xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 0., hide=0
        xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=3, hide=0  ;dotted line
        xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
        xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0, hide=0
        xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
        xd_sState.COUNT123=1

        xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
        xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
        xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution
        result=(xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1+xd_sState.dataXYZdistance*xd_sState.dataXYZdistance-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
        angle=acos(result)*!radeg
        print,'2 lines angle:',angle,'degree'
        print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
        statusString_old = STRING(xd_sState.dataXYZdistance_old, $
          FORMAT='("Distance = ", F6.2, " mm")')
        statusString = STRING(xd_sState.dataXYZdistance, $
          FORMAT='("Distance = ", F6.2, " mm")')
        xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
        statusString1 = STRING(angle, $
          FORMAT='("Angle = ", F6.2, " degree")')
        demo_putTips, xd_sState, 'locat', 11, /LABEL
        demo_putTips, xd_sState, statusString+','+statusString1, 12
        xd_sState.ckr1->setProperty,hide=1
        xd_sState.ckg->setProperty,hide=1
        xd_sState.ckr->setProperty, STRINGS= statusString_old,hide=0
        xd_sState.ckb->setProperty, STRINGS= statusString,hide=0
        xd_sState.ckag->setProperty,color=[153,255,0], STRINGS= statusString1,hide=0
      endif else begin
        if(xd_sState.COUNT123  eq 1) THEN begin
          xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 0., hide=0
          xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=3, hide=0   ;dotted line
          xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]
          xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 1,linestyle=0, hide=0
          xd_sState.oLine1->setProperty, DATA = [[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]
          xd_sState.COUNT123=2
          xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ2], [xd_sState.dataXYZ0]]) * xd_sState.vol_HU_cube_resolution

          xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
          xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

          result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2-xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1)/(2*xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance)
          angle=acos(result)*!radeg
          print,'2 lines angle:',angle,'degree'
          print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
          statusString = STRING(xd_sState.dataXYZdistance, $
            FORMAT='("Distance = ", F6.2, " mm")')
          statusString1 = STRING(angle, $
            FORMAT='("Angle = ", F6.2, " degree")')
          statusString_old = STRING(xd_sState.dataXYZdistance_old, $
            FORMAT='("Distance = ", F6.2, " mm")')
          xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance

          demo_putTips, xd_sState, 'locat', 11, /LABEL
          demo_putTips, xd_sState, statusString+statusString1, 12
          xd_sState.ckr1->setProperty,hide=1
          xd_sState.ckb->setProperty,hide=1
          xd_sState.ckg->setProperty, STRINGS= statusString_old,hide=0
          xd_sState.ckr->setProperty, STRINGS= statusString,hide=0
          xd_sState.ckag->setProperty,color=[0,153,255], STRINGS= statusString1,hide=0
        endif else begin
          if (xd_sState.COUNT123  eq 2) THEN begin
            xd_sState.oLine3->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]
            xd_sState.oLine3->setProperty, ALPHA_CHANNEL = 1,linestyle=0, hide=0  ;dotted line
            xd_sState.oLine1->setProperty, ALPHA_CHANNEL = 0.,linestyle=3, hide=0

            xd_sState.oLine2->setProperty, ALPHA_CHANNEL = 1,linestyle=0, hide=0
            xd_sState.oLine2->setProperty, DATA = [[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]
            xd_sState.COUNT123=0
            xd_sState.dataXYZdistance = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ1]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance1 = DISTANCE_MEASURE([[xd_sState.dataXYZ0], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution
            xd_sState.dataXYZdistance2 = DISTANCE_MEASURE([[xd_sState.dataXYZ1], [xd_sState.dataXYZ2]]) * xd_sState.vol_HU_cube_resolution

            result=(xd_sState.dataXYZdistance*xd_sState.dataXYZdistance+xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance1-xd_sState.dataXYZdistance2*xd_sState.dataXYZdistance2)/(2*xd_sState.dataXYZdistance1*xd_sState.dataXYZdistance)
            angle=acos(result)*!radeg
            print,'2 lines angle:',angle,'degree'
            print, '2 points distance: ', xd_sState.dataXYZdistance,'mm'
            statusString = STRING(xd_sState.dataXYZdistance, $
              FORMAT='("Distance = ", F6.2, " mm")')
            statusString1 = STRING(angle, $
              FORMAT='("Angle = ", F6.2, " degree")')
            statusString_old = STRING(xd_sState.dataXYZdistance_old, $
              FORMAT='("Distance = ", F6.2, " mm")')
            xd_sState.dataXYZdistance_old=xd_sState.dataXYZdistance
            demo_putTips, xd_sState, 'locat', 11, /LABEL
            demo_putTips, xd_sState, statusString+statusString1, 12
            xd_sState.ckr1->setProperty,hide=1
            xd_sState.ckr->setProperty,hide=1
            xd_sState.ckb->setProperty, STRINGS= statusString_old,hide=0
            xd_sState.ckg->setProperty, STRINGS= statusString,hide=0
            xd_sState.ckag->setProperty,color=[255,102,102], STRINGS= statusString1,hide=0
          endif
        endelse
      endelse
    endif
  endelse

END





;need to define 30 colors for 3D IDLgrPolygon objects
;color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
FUNCTION COLORTABLE_rgb_3D, ColorIndex, ColorIndexBGR
  color_rgb = BYTARR(31, 3)
;  color_rgb[0:6, *] = [[254,206,180],[255,134,1],[255,255,1],[1,255,25],[1,1,255],[92,1,255],[255,1,243]]
  color_rgb = [['00'x, '00'x, '00'x], $;0
               ['FF'x, '66'x, '66'x], $;1
               ['FF'x, 'FF'x, '00'x], $;2
               ['33'x, '99'x, 'CC'x], $;3
               ['00'x, '99'x, '33'x], $;4
               ['FF'x, 'FF'x, 'CC'x], $;5
               ['FF'x, '99'x, '33'x], $;6
               ['00'x, '99'x, 'CC'x], $;7
               ['FF'x, '66'x, '00'x], $;10
               ['FF'x, 'FF'x, '66'x], $;11
               ['00'x, '99'x, '66'x], $;12
               ['CC'x, '00'x, '66'x], $;13
               ['00'x, '99'x, '99'x], $;14
               ['FF'x, 'CC'x, '33'x], $;15
               ['99'x, 'CC'x, '00'x], $;16
               ['CC'x, '99'x, '99'x], $;17
               ['CC'x, 'CC'x, '99'x], $;18
               ['FF'x, '99'x, 'CC'x], $;19
               ['FF'x, 'FF'x, '33'x], $;20
               ['FF'x, '99'x, '00'x], $;21
               ['CC'x, '99'x, '00'x], $;22
               ['FF'x, 'FF'x, '99'x], $;23
               ['66'x, '99'x, '99'x], $;24
               ['CC'x, '66'x, '00'x], $;25
               ['99'x, '99'x, '99'x], $;26
               ['CC'x, 'CC'x, '33'x], $;27
               ['CC'x, '66'x, '33'x], $;28
               ['FF'x, 'CC'x, '99'x], $;29
               ['CC'x, '66'x, 'FF'x], $;30
               ['CC'x, 'CC'x, 'CC'x], $;8
               ['FF'x, 'FF'x, 'FF'x]];9
 color_rgb = TRANSPOSE(color_rgb)
 ColorIndex = ['000000'x, $
               'FF6666'x, 'FFFF00'x, '3399CC'x, '009933'x, 'FFFFCC'x, $
               'FF9933'x, '0099CC'x, 'FF6600'x, $
               'FFFF66'x, '009966'x, 'CC0066'x, '009999'x, 'FFCC33'x, $
               '99CC00'x, 'CC9999'x, 'CCCC99'x, 'FF99CC'x, 'FFFF33'x, $
               'FF9900'x, 'CC9900'x, 'FFFF99'x, '669999'x, 'CC6600'x, $
               '999999'x, 'CCCC33'x, 'CC6633'x, 'FFCC99'x, 'CC66FF'x, 'CCCCCC'x, 'FFFFFF'x]
 ColorIndexBGR = ['000000'x, $
               '6666FF'x, '00FFFF'x, 'CC9933'x, '339900'x, 'CCFFFF'x, $
               '3399FF'x, 'CC9900'x, '0066FF'x, $
               '66FFFF'x, '669900'x, '6600CC'x, '999900'x, '33CCFF'x, $
               '00CC99'x, '9999CC'x, 'CC99CC'x, 'CC99FF'x, '33FFFF'x, $
               '0099FF'x, '0099CC'x, '99FFFF'x, '999966'x, '0066CC'x, $
               '999999'x, '33CCCC'x, '3366CC'x, '99CCFF'x, 'FF66CC'x, 'CCCCCC'x, 'FFFFFF'x]
 RETURN, color_rgb
END

PRO test_plot_colortable_precolors
  wTlb = WIDGET_BASE(/ROW)
  ;temp = LONARR(31)
  draw = LONARR(31)
  for i=0,30 do begin
    ;temp[i] = WIDGET_BASE(wRow2Base, /COLUMN)
    Draw[i] = WIDGET_DRAW(wTlb, XSIZE=30, YSIZE=120)
  endfor
  WIDGET_CONTROL, wTlb, /REALIZE
  
  color_rgb = COLORTABLE_rgb_3D(ColorIndex, ColorIndexBGR)
  DEVICE, DECOMPOSED=1
  for i=0,30 do begin
    WIDGET_CONTROL, Draw[i], GET_VALUE=winID
    WSET, winID
    ERASE, ColorIndexBGR[i]
  endfor
  XMANAGER, 'test_plot_colortable_precolors', wTlb, /NO_BLOCK
  odraw_img = wTlb->READ()
  odraw_img->getProperty, DATA = draw_img
  help, draw_img
END