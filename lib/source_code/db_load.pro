;把变量转换成 字符串变量
FUNCTION TYPE_CONVERSION,value,type
type=STRUPCASE(type)
CASE type of
  'VARCHAR':BEGIN
          rvalue=string(value)
     END
  'CHAR': BEGIN
        rvalue=string(value)
    END
  else : BEGIN
          rvalue=value
    END 
ENDCASE
  return,rvalue
END


;根据source，user_id, password连接数据库
;连接数据库成功，返回1，连接不成功，返回0
PRO CONNECT,database
 ;connect database
  if(~obj_valid((*database).obj) ) then $
    (*database).obj=OBJ_NEW('IDLdbDatabase')
 
  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
 ;   user_id='gjia' ; 'root'
 ;   password='ltpltp123' ; 'root'
 ;   (*database).obj->CONNECT,DATASOURCE=(*database).source,USER_ID=user_id,PASSWORD=password
    OBJ_DESTROY,(*database).obj
    (*database).isconnect=0
    RETURN
  ENDIF
(*database).obj->CONNECT,DATASOURCE=(*database).source,USER_ID=(*database).user_id,PASSWORD=(*database).password
(*database).obj->GETPROPERTY,IS_CONNECTED = connectStat
IF connectStat EQ 0 THEN BEGIN
  OBJ_DESTROY,(*database).obj
  (*database).isconnect=0
  RETURN
ENDIF
(*database).isconnect=1
END

;初始化数据库，将datasource, user_id, password做成列表，放在一起
PRO INITDatabase,state
; initialize the database
;Run isolated
;Connect test
ENABLE_DATABASE=0; 开启和关闭数据库功能

account_N = 4

datasource = STRARR(account_N)
user_id = STRARR(account_N)
password = STRARR(account_N)

;gjia
datasource[3]='IDL_vis'
user_id[3]='gjia' ; 'root'
password[3]='ltpltp123' ; 'root'

;xwt
datasource[2]='IDL_vis'
user_id[2]='root' ; 'root'
password[2]='root' ; 'root'

;huangxunan
datasource[1]='IDL_vis'
user_id[1]='root' ; 'root'
password[1]='010201' ; 'root'


;haojiaxue
datasource[0]='root'
user_id[0]='root' ; 'root'
password[0]='0315' ; 'root'


  if(N_elements(state) eq 0) then state={database:ptr_new()}

  IF PTR_VALID(state.database) NE 1 THEN BEGIN
   IF DB_EXISTS() EQ 0 THEN BEGIN
;disable this message! 2020/5/19
;    msg = DIALOG_MESSAGE('We do not have ODBC! Please install it.',/Error)
    RETURN
   ENDIF
 ENDIF
    
  ;create database object and set datasource
  oDatabase=OBJ_NEW('IDLdbDatabase')
  
  ; Establish error handler. When errors occur, the index of the
  ; error is returned in the variable Error_status:
  CATCH, Error_status
  ;This statement begins the error handler:
  IF Error_status NE 0 THEN BEGIN
    PRINT, 'Data source driver is not available!'
    CATCH, /CANCEL
    IF N_ELEMENTS(i) THEN BEGIN
      IF i LT account_N-2 THEN i=i+1
      print, i
    ENDIF
    goto, recal
  ENDIF
  ;check the validity of the Excel database
  sources =oDatabase.GETDATASOURCES()

  found = 0
  FOR i=0, account_N-1 DO BEGIN
recal:    IF N_ELEMENTS(sources) EQ 0 THEN RETURN ELSE index = WHERE(sources.DATASOURCE EQ datasource[i],count)
    print, 'i = ', i
    IF count EQ 1 THEN BEGIN
      found = 1
      state.database=ptr_new({obj:oDatabase,source:datasource[i],user_id:user_id[i],password:password[i],isconnect:0})
      print, *(state.database)
      connect,state.database
      if((*state.database).isconnect eq 1) then begin
        found = 2
        ;close
        (*state.database).obj->cleanup
        break
      endif
    ENDIF
  ENDFOR
  
  IF found LT 2 THEN BEGIN
;disable this message 2020/5/19
;      msg = DIALOG_MESSAGE('The system does not have ODBC or MySQL installed!',/Error)
      OBJ_DESTROY,oDatabase
      PTR_FREE, state.database
      RETURN
  ENDIF
  
  RETURN  
  
END

;执行SQL script得到病例列表
;to excute sql script from an file,and return the recorde(if exists) by the oRecordset
PRO EXCUTE_SQL_script,database,sql_file_path,oRecordset=oRecordset
  if ~FILE_TEST(sql_file_path) then BEGIN
    msg=DIALOG_MESSAGE('the file: '+sql_file_path+' do not exists!',/Error)
    RecorSet=obj_new();
    RETURN
  end
  Openr,sql_file,sql_file_path,/get_lun
  if sql_file eq -1 then begin
    msg=DIALOG_MESSAGE('file error!',/Error)
    RecorSet=obj_new();
    RETURN
  end
  sql=''
  temp=''
  if (*database).isconnect eq 0 then begin
    msg=DIALOG_MESSAGE('database has not been connected!',/Error)
    RecorSet=obj_new();
    FREE_LUN,sql_file
    RETURN
  endif
  WHILE(~EOF(sql_file)) DO BEGIN
    READF,sql_file,temp
    ;if(temp ne '') then $
     ; oRecordset = OBJ_NEW('IDLdbRecordset',(*database).obj,SQL=temp)
     IF TEMP NE '' THEN BEGIN
       P=STRPOS(TEMP,'#')
       IF P EQ -1 THEN BEGIN
          sql+=temp
       ENDIF ELSE BEGIN 
          SQL+=STRMID(TEMP,0,P)
       ENDELSE
     ENDIF
  ENDWHILE
  ;sql=STRMID(sql,0,STRLEN(sql)-1)
  pos=0
  len=strlen(sql)
  ;excute
  sql=STRTRIM(sql,2)
 ; sql=STRUPCASE(sql) 
  sql_c=STRSPLIT(SQL,';', /EXTRACT)
  FOR i=0,n_elements(sql_c)-1 DO BEGIN
      temp=sql_c[i]+';'
      ;TEMP=STRJOIN(STRSPLIT(TEMP,STRING(10),/EXTRACT),' ')
      p=STRPOS(temp,' ')
      IF P NE -1 THEN BEGIN
        subtemp=STRUPCASE(STRMID(TEMP,0,P))
        IF STRPOS(subtemp,'SELECT',P,/REVERSE_SEARCH) NE -1 ||  STRPOS(subtemp,'SHOW',P,/REVERSE_SEARCH) NE -1 THEN BEGIN
           oRecordset=OBJ_NEW('IDLdbRecordset',(*database).obj,SQL=temp);  
        ENDIF ELSE BEGIN $
            (*database).obj->ExecuteSQL,temp
        ENDELSE
      END
    ENDFOR  
;  WHILE(pos LT len-1) DO BEGIN
;    npos=STRPOS(sql,';',pos)
;    if npos EQ -1 THEN $
;      npos=len-1
;    temp=STRMID(sql,pos,npos-pos+1)
;    temp=STRTRIM(temp,2)
;    IF (TEMP NE '' OR TEMP NE string(10B))  THEN $
;        oRecordset = OBJ_NEW('IDLdbRecordset',(*database).obj,SQL=temp)
;     pos=npos+1
;  ENDWHILE
 
    
   FREE_LUN,sql_file
END


;得到病例记录，从表格中，根据鼠标相应的位置，选取相应的病例
FUNCTION GET_RECORD,oRecordSet
  ;labels=['a','b','c','d','e','f']
  ;d1=[123,234]
  ;data=[d1,d1,d1,d1]
  IF N_ELEMENTS(oRecordSet) eq 0 THEN BEGIN
    RETURN,[]
  ENDIF
  oRecordset->GETPROPERTY,field_info = fieldinfo
  nfield=n_elements(fieldinfo)
  labels=STRLOWCASE(fieldinfo.FIELD_NAME)
  colwidth=intarr(nfield)
  type=fieldinfo.TYPE_NAME
  data=[]
  max_strlen=max([max(strlen(labels)),10])
  IF oRecordset->MOVECURSOR(/first) EQ 1 THEN BEGIN
    d=[]
    FOR j=0,nfield-1 DO BEGIN
      Value =TYPE_CONVERSION(oRecordset->GETFIELD(j),TYPE(j))
      d=create_struct(d,labels(j),Value)
      colwidth[j]=max([colwidth[j],strlen(strtrim(string(value),2)),strlen(labels[j])])
    ENDFOR
    data=[data,d]
    WHILE oRecordset->MOVECURSOR(/next) DO BEGIN
      d=[]
      FOR j=0,nfield-1 DO BEGIN
        Value =TYPE_CONVERSION(oRecordset->GETFIELD(j),TYPE(j))
        d=create_struct(d,labels(j),Value)
        colwidth[j]=max([colwidth[j],strlen(strtrim(string(value),2)),strlen(labels[j])])
      ENDFOR
      data=[data,d]
    ENDWHILE
  ENDIF else begin
    d=[]
    FOR j=0,nfield-1 DO BEGIN
      Value =''
      d=create_struct(d,labels(j),Value)
      colwidth[j]=max([colwidth[j],strlen(strtrim(string(value),2)),strlen(labels[j])])
    ENDFOR
    data=[data,d]
  Endelse
    result={data:data,colwidth:colwidth} 
  RETURN,result
END








;根据选中的病例，把详细信息展示出来
PRO SET_TABLE_INFO,wDBtable,oRecordSet

 result=GET_RECORD(oRecordSet)
 if n_elements(result) GE 1 then begin
  labels=STRLOWCASE(tag_names(result.data))
  col_widths= MAKE_ARRAY(1, n_elements(labels), /INTEGER, VALUE =110)
  WIDGET_CONTROL, wDBtable,COLUMN_LABELS=labels
  WIDGET_CONTROL,wDBtable,/UNITS
  WIDGET_CONTROL,wDBtable,COLUMN_WIDTHS=col_widths
  WIDGET_CONTROL,wDBtable,xsize=n_elements(labels)  
  WIDGET_CONTROL,wDBtable,ysize=n_elements(result.data)
  WIDGET_CONTROL, wDBtable,set_value=result.data
  endif else begin
  endelse
END


















;db_load的窗口设置
PRO DB_DIALOG_EVENT,ev
WIDGET_CONTROL,ev.id,get_uvalue=uval
case uval[0] of
  'OK' :begin
    WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
    (*(info.button_info)).cancel=0;
    select=WIDGET_INFO(info.table,/TABLE_SELECT)
    if(select[0] ne select[2]) then begin
        WIDGET_CONTROL,info.table,get_value=value,/no_copy
        connect,info.database
        sql='use test;'
        (*info.database).obj->ExecuteSQL,sql
        t=select[1]
        b=select[3];   
        if n_elements(value[t].id) and value[t].id ne '' then begin
          sql='set global max_allowed_packet=268435456;'
          (*info.database).obj->ExecuteSQL,sql
          
          sql='call id_search('+STRTRIM(string(value[t].id),2)+');'
          oRecordset=OBJ_NEW('IDLdbRecordset',(*info.database).obj,sql=sql);
          if(obj_valid(oRecordset)) then begin
            oRecordset->GETPROPERTY,field_info = fieldinfo
            type=fieldinfo.TYPE_NAME
            dcm_dir =TYPE_CONVERSION(oRecordset->GETFIELD(10),TYPE(10))
            WHILE (((I = STRPOS(dcm_dir, '*'))) NE -1) DO STRPUT, dcm_dir, '\', I
            print, 'dcm_dir = ', dcm_dir
            image_read, dcm_dir, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal
            IF N_ELEMENTS(vol_HU_cube) GT 0 THEN BEGIN
              (*(info.button_info))=create_struct((*(info.button_info)),'data',vol_HU_cube, 'resolution', vol_HU_cube_resolution, 'age', patient_age, 'path', dcm_dir)
            ENDIF ELSE BEGIN
              msg = DIALOG_MESSAGE('image reading failed',/Error)
              (*info.database).obj->cleanup
              WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
              WIDGET_CONTROL,ev.top,/destroy
              return
            ENDELSE
          endif
          (*info.database).obj->cleanup
          WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
          WIDGET_CONTROL,ev.top,/destroy
        endif else begin
           (*info.database).obj->cleanup
           
        endelse
    endif 
     ; WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
    end
    
    'NEW' :begin
      WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
      (*(info.button_info)).cancel=0;
      state_temp = info.state
      DB_ADD,state=state_temp
      IF PTR_VALID(state_temp.vol_HU_cube) THEN BEGIN
        info.state = state_temp
        (*(info.button_info))=create_struct((*(info.button_info)),'data',*(state_temp.vol_HU_cube), 'resolution', state_temp.vol_HU_cube_resolution, 'age', state_temp.patient_age, 'path', state_temp.additionalDir)
        print, 'resolution', state_temp.vol_HU_cube_resolution
      ENDIF ELSE BEGIN
        (*(info.button_info)).cancel=1;
      ENDELSE
      WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
      WIDGET_CONTROL,ev.top,/destroy
    end
    
   'CANCEL':begin
   WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
   (*(info.button_info)).cancel=1;
   WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
   WIDGET_CONTROL,ev.top,/destroy
    end
   'TABLE':begin
          WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
            IF ev.type EQ 4 THEN BEGIN
                IF ev.sel_left NE -1 THEN BEGIN
              
                ;widget_info,wDBtable,
               
                select=WIDGET_INFO(info.table,/TABLE_SELECT)
                WIDGET_CONTROL,info.table,get_value=value  
                table_view=WIDGET_INFO(info.table,/TABLE_VIEW) 
                    if select[1] eq select[3] then begin
                      col_widths = WIDGET_INFO(info.table, /COLUMN_WIDTHS)
                     ; scol_size=WIDGET_INFO(INFO.table,/y_scroll_size)
                      select[0]=0
                      select[2]=n_elements(tag_names(value))-1;
                      widget_control,info.table,SET_TABLE_SELECT=select
                      widget_control,info.table,SET_TABLE_VIEW=table_view
                     ; widget_control,info.table,set_y_scrol_size=scol_size
                    endif
                ENDIF
           ENDIF
         WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
    end
   'SEARCH' :BEGIN
        WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
        ;get text value
        WIDGET_CONTROL,info.menu.wName,get_value=str_name
        ;get selected value 
        ;modaling
        select=WIDGET_INFO(info.menu.wmodal,/droplist_select)
        WIDGET_CONTROL,info.menu.wmodal,get_value= value
        str_modal=value[select];
        ;date1 year
        select=WIDGET_INFO(info.menu.wDate1.y,/droplist_select)
        WIDGET_CONTROL,info.menu.wDate1.y,get_value= value
        str_y1=value[select];
        ;date1 month
        select=WIDGET_INFO(info.menu.wDate1.m,/droplist_select)
        WIDGET_CONTROL,info.menu.wDate1.m,get_value= value
        str_m1=value[select];
        ;date1 day
        select=WIDGET_INFO(info.menu.wDate1.d,/droplist_select)
        WIDGET_CONTROL,info.menu.wDate1.d,get_value= value
        str_d1=value[select];
        
        ;date2 year
        select=WIDGET_INFO(info.menu.wDate2.y,/droplist_select)
        WIDGET_CONTROL,info.menu.wDate2.y,get_value= value
        str_y2=value[select];
        ;date2 month
        select=WIDGET_INFO(info.menu.wDate2.m,/droplist_select)
        WIDGET_CONTROL,info.menu.wDate2.m,get_value= value
        str_m2=value[select];
        ;date2 day
        select=WIDGET_INFO(info.menu.wDate2.d,/droplist_select)
        WIDGET_CONTROL,info.menu.wDate2.d,get_value= value
        str_d2=value[select];
        
        
        ;check_valid,
        
        ;
        ;search order by {name,modal,time} if eq ''
       
        if str_y1 ne '' then  begin
        d1=string(str_y1,str_m1,str_d1,format='(i04,i02,i02)')
        endif else begin
          d1=str_y1
        endelse
        
        if str_y2 ne '' then begin
          d2=string(str_y2,str_m2,str_d2,format='(i04,i02,i02)')
        endif else begin
          d2=str_y2
        endelse
        
        connect,info.database
        if strpos(d1,'*') ne -1 then d1=''
        if strpos(d2,'*') ne -1 then d2=''
        sql='use test;'
        (*info.database).obj->ExecuteSQL,sql
        sql="call search_info('"+ $
        STRJOIN(STRSPLIT(str_name, /EXTRACT,"'"),"\'") + $
        "','"+STRJOIN(STRSPLIT(str_modal, /EXTRACT,"'"),"\'")+ $
        "','"+ STRJOIN(STRSPLIT(d1, /EXTRACT,"'"),"\'")+ $
        "','"+ STRJOIN(STRSPLIT(d2, /EXTRACT,"'"),"\'")+ $
        "');"    
        oRecordset=OBJ_NEW('IDLdbRecordset',(*info.database).obj,sql=sql); 
        SET_TABLE_INFO,info.table,oRecordSet
        (*info.database).obj->cleanup
        WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
    END
    'DELETE':BEGIN
        WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
        
        select=WIDGET_INFO(info.table,/TABLE_SELECT) 
        if(select[0] ne select[2]) then begin
          WIDGET_CONTROL,info.table,get_value=value,/no_copy
          connect,info.database
          sql='use test;'
          (*info.database).obj->ExecuteSQL,sql      
          t=select[1]
          b=select[3];
          for i=t,b do begin
            if n_elements(value[i].id) and value[i].id ne '' then begin
              sql='call delete_info('+string(value[i].id)+');'
              (*info.database).obj->ExecuteSQL,sql
            endif
          endfor
        
          if n_elements(value[t].id) and value[t].id ne '' then begin
              ;delete value[t]~value[b]
            if t ne 0 then begin
              tvalue=value[0:t-1]
            endif else begin
              tvalue=[]
            endelse
            if b ne n_elements(value)-1 then begin
              bvalue=value[b+1:n_elements(value)-1]
            endif else begin
             bvalue=[];
            endelse
            value=[tvalue,bvalue]
          endif 
          (*info.database).obj->cleanup
          WIDGET_CONTROL,info.table,ysize=n_elements(value)
          WIDGET_CONTROL,info.table,set_value=value,/no_copy
        endif         
        WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
      END
   'DETAILS': BEGIN
         WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
         select=WIDGET_INFO(info.table,/TABLE_SELECT) 
         t=select[1]
         l=select[0]
         r=select[2]
         b=select[3]
         if t ne -1 and t eq b and t ne r then begin
         WIDGET_CONTROL,info.table,get_value=value,/no_copy
         if n_elements(value[t].id) and value[t].id ne '' then $
         DB_ADD,state=info,case_id=value[t].id   
         WIDGET_CONTROL,info.table,set_value=value,/no_copy
         endif
         ;          
         ;          d_objworld2Event,{ $
         ;          top: ev.top, $
         ;            handler: 0L, $
         ;            id: wDel $
         ;            }
         wButtonSearch=info.wButtonSearch
         WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy  
         DB_DIALOG_EVENT,{ $
          top: ev.top, $
          handler: 0L, $
          id:wButtonSearch $
         }
         
        
    END   
      
      
   else: begin
    
    end 
endcase


END

;db_load的主程序
PRO DB_LOAD,state
;load DATA by database

;THE RETURN INFO OF DIALOG
button_info=ptr_New({cancel:1})
;create database dialog
SCR_SIZE={x:1300,y:400}
MARGIN=1
MENUBORDER={x:SCR_SIZE.x,y:50}
MENUSIZE={x:MENUBORDER.x-MARGIN*2,y:MENUBORDER.y-MARGIN*2}
TABLEBORDER={x:SCR_SIZE.x,y:SCR_SIZE.y-MENUBORDER.y}
TABLESIZE={x:TABLEBORDER.x-MARGIN*2,y:TABLEBORDER.y-MARGIN*2-9}


SIZE={scr:scr_size,ma:margin,mb:MENUBORDER,me:MENUSIZE,tb:TABLEBORDER,ta:TABLESIZE}


IF N_Elements(state) NE 0 THEN BEGIN
wDBdialog= WIDGET_BASE(GROUP_LEADER =state.wTopBase, $
  xsize =SCR_SIZE.x,ysize =SCR_SIZE.y, $
  xOffset =400, $
  yOffset =200, $
  TLB_FRAME_ATTR = 1, $
  /Floating,$
  /Base_Align_Center,$
  /modal,$
  /column,$
  title ='Database loading')
ENDIF ELSE BEGIN
  RETURN
;  wDBdialog= WIDGET_BASE($
;    xsize =SCR_SIZE.x,ysize =SCR_SIZE.y, $
;    xOffset =400, $
;    yOffset =200, $
;    TLB_FRAME_ATTR = 1,$
;    /Base_Align_Center,$
;     /column,$
;    title ='search')
ENDELSE
wDBmenuBorder=WIDGET_BASE(wDBdialog,$
  ;/GRID_LAYOUT,$
  ;ROW=3,$
  xsize=size.mb.x,$
  ysize=size.mb.y,$
  /COLUMN ,$
  /FRAME $
  )
 wDBmenu=WIDGET_BASE($
  wDBmenuBorder,$
  /ALIGN_CENTER,$
  /BASE_ALIGN_CENTER,$
  /ROW,$
  xsize=size.me.x,ysize=size.me.y)
  
  
  wSubform122=WIDGET_BASE(wDBmenu,/COLUMN,/frame)
  wDBbuttonGroup2=WIDGET_BASE(wSubForm122,/ROW,/ALIGN_CENTER,/BASE_ALIGN_CENTER)

  wDBdialogButtonNew=WIDGET_BUTTON($
    wDBbuttonGroup2,$
    /align_center,$
    VALUE=state.icon_dir+'newcase.bmp', /BITMAP,$
    uname='New',$
    uval='NEW'$
    )

  
 wSubForm11=WIDGET_BASE(wDBmenu,/row,/tab_mode)
 wlebase1=WIDGET_BASE(wSubForm11,/ROW,/ALIGN_CENTER,/frame)   
 
 wlebase11=WIDGET_BASE(wlebase1,/row)
 label_name=WIDGET_LABEL(wlebase11,value='name:',uname='name')   
 edit_name=WIDGET_TEXT(wlebase11,uname='ediname',uval='NAME',/EDITABLE,xsize=10)
 
 str_modal=['','CT','MR','DSA','MRI']
 wlebase12=WIDGET_BASE(wlebase1,/row)     
 label_modaling=WIDGET_LABEL(wlebase12,value='modality:')
 drop_modaling=WIDGET_Droplist(wlebase12,uval='MODALING',value=str_modal)


 wlebase2=WIDGET_BASE(wSubForm11,/row,/ALIGN_CENTER,/frame)
 label3=WIDGET_LABEL(wlebase2,value='time:')
 wlebase21=WIDGET_BASE(wlebase2,/row)
 str_y=[[''],[STRTRIM(string(indgen(1,16)-8+2015),2)]]
 str_m=[[''],[STRTRIM(string(indgen(1,12)+1),2)]]
 str_d=[[''],[STRTRIM(string(indgen(1,31)+1),2)]]
 drop1_y=WIDGET_Droplist(wlebase21,uval='DROP',value=str_y)
 drop1_m=WIDGET_Droplist(wlebase21,title='-',uval='DROP',value=str_m)
 drop1_d=WIDGET_Droplist(wlebase21,title='-',uval='DROP',value=str_d)
 
 label4=WIDGET_LABEL(wlebase2,value='to')
 wlebase22=WIDGET_BASE(wlebase2,/row)
 drop2_y=WIDGET_Droplist(wlebase22,uval='DROP',value=[[''],[STRTRIM(string(indgen(1,16)-8+2015),2)]])
 drop2_m=WIDGET_Droplist(wlebase22,title='-',uval='DROP',value=[[''],[STRTRIM(string(indgen(1,12)+1),2)]])
 drop2_d=WIDGET_Droplist(wlebase22,title='-',uval='DROP',value=[[''],[STRTRIM(string(indgen(1,31)+1),2)]])

wDate1={y:drop1_y,m:drop1_m,d:drop1_d}
sDate={y:str_Y,m:str_m,d:str_d}
wDate2={y:drop2_y,m:drop2_m,d:drop2_d}
wMenu={wName:edit_name,wModal:drop_modaling,wDate1:wDate1,wDate2:wDate2}


wSubform12=WIDGET_BASE(wDBmenu,/COLUMN,/frame)
wDBbuttonGroup=WIDGET_BASE(wSubForm12,/ROW,/ALIGN_CENTER,/BASE_ALIGN_CENTER)

wDBdialogButtonSearch=WIDGET_BUTTON($
  wDBbuttonGroup,$
  /align_center,$
  VALUE=state.icon_dir+'search.bmp', /BITMAP,$
  uname='search',$
  uval='SEARCH'$
  )
; wDBdialogButtonAdd=WIDGET_BUTTON($
;   wDBbuttonGroup,$
;   /align_center,$
;   value='Add',$
;   uname='Add',$
;   uval='ADD'$
;  )
  wDBdialogButtonUpdate=WIDGET_BUTTON($
    wDBbuttonGroup,$
    /align_center,$
    VALUE=state.icon_dir+'details.bmp', /BITMAP,$
    uname='DETAILS',$
    uval='DETAILS'$
  )
  wDBdialogButtonDelete=WIDGET_BUTTON($
    wDBbuttonGroup,$
    /align_center,$
    VALUE=state.icon_dir+'deletecase.bmp', /BITMAP,$
    uname='delete',$
    uval='DELETE'$
  )
  wDBdialogButtonLoad=WIDGET_BUTTON($
    wDBbuttonGroup,$
    /align_center,$
    VALUE=state.icon_dir+'load.bmp', /BITMAP,$
    uname='load',$
    uval='OK'$
    )
 

 ;widget_control, wDBdialogButtonLoad, sensitive=1
  
;wDBdialogButtonOk=WIDGET_BUTTON($
;  wDBbuttonGroup,$
;  /align_center,$
;  value='ok',$
;  uname='ok',$
;  uval='OK'$
;  )
;wDBdialogButtonCancel=WIDGET_BUTTON($
;  wDBbuttonGroup,$
;  /align_center,$
;  value='cancel',$
;  uname='cancel',$
;  uval='CANCEL'$
; )
   
wDBtableBorder=WIDGET_BASE(wDBdialog,$
  xsize=size.tb.x,$
 ysize=size.tb.y,$
 /COLUMN)
wDBtable=WIDGET_TABLE(wDBtableBorder,$
  uval='TABLE',$
 scr_xsize=size.ta.x,scr_ysize=size.ta.y,$
  /RESIZEABLE_COLUMNS,$
  /ALIGN_CENTER,$
  /ROW_MAJOR,$
  /FRAME,$
  /NO_ROW_HEADERS);,$
widget_control,wDBtable,/ALL_TABLE_EVENTS

  
INITDatabase,state  ;if error happens please disable this

IF ~PTR_VALID(state.database) THEN BEGIN
reselect:  first_file = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/',/MUST_EXIST, filter=['*.dcm', '*.nii'], TITLE='Select DICOM Patient File')
  afn = file_dirname(first_file)
  IF (STRLEN(afn) EQ 0) OR (STRLEN(first_file) EQ 0) THEN RETURN
  filearr=file_search(afn,'*.{nii,dcm}',count=num)
  IF num EQ 0 THEN BEGIN
    dcm_info = DIALOG_MESSAGE('Please sort images first!', /INFORMATION, TITLE='Suggestions')
    Separate_to_series, afn, afn_des
    first_file = DIALOG_PICKFILE(PATH=afn_des,/MUST_EXIST, filter=['*.dcm', '*.nii'], TITLE='Select DICOM Patient File')
    afn = file_dirname(first_file)
    IF (STRLEN(afn) EQ 0) OR (STRLEN(first_file) EQ 0) THEN RETURN
    filearr=file_search(afn,'*.{nii,dcm}',count=num)
    IF num NE 0 THEN BEGIN
      image_read, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal
      ;确保vol_HU_iso不为空
      IF((size(vol_HU_cube))[0] GT 2) THEN BEGIN
        ptr_free,state.vol_HU_cube
        state.vol_HU_cube = ptr_new(vol_HU_cube)
        state.vol_HU_cube_resolution = vol_HU_cube_resolution
        state.patient_age = patient_age
        state.additionalDir = ptr_new(afn)
        print, 'patient age 2 = ', state.patient_age
        print, 'database resolution', state.vol_HU_cube_resolution
      ENDIF
    ENDIF
  ENDIF ELSE BEGIN
    ;GJ 2020/7/5; removing the bugs about too few images
    IF num GT 0 THEN BEGIN
      image_read, afn, vol_HU_cube, vol_HU_cube_resolution, patient_age, direction, origin_loc, p_modal=p_modal
    ENDIF ELSE BEGIN
      dcm_info = DIALOG_MESSAGE('Too few images! Plese select another folder!', /INFORMATION, TITLE='Suggestions')
      goto, reselect
    ENDELSE
    ;确保vol_HU_iso不为空
    ;@GJ, 2024/10/11, select the sav file from large bore MPI
    afn_sav = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X8mm\'
    sav_file = DIALOG_PICKFILE(PATH=afn_sav,/MUST_EXIST, filter=['*.sav'], TITLE='Select LY MPI Sav File')
    restore, sav_file
;    IF N_ELEMENTS(save_file) GT 0 THEN BEGIN
;      restore, sav_file;'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X8mm\vol_HU_cube.sav'
;      ;gaussian blur may be added
;      ;@GJ, 2024/10/7
;    END
;    restore, 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X10mm\vol_HU_cube.sav'
    IF((size(vol_HU_cube))[0] GT 2) THEN BEGIN
      ptr_free,state.vol_HU_cube
      state.vol_HU_cube = ptr_new(BYTSCL(vol_HU_cube))
      state.vol_HU_cube_resolution = vol_HU_cube_resolution
      state.patient_age = patient_age
      state.additionalDir = ptr_new(afn)
      print, 'patient age 3 = ', state.patient_age
      print, 'database resolution', state.vol_HU_cube_resolution
    ENDIF
  ENDELSE
  RETURN
ENDIF

connect,state.database
if((*state.database).isconnect ne 1) then begin
    msg = DIALOG_MESSAGE('connect failed',/Error)
    return 
end
sql='use test;'
(*state.database).obj->ExecuteSQL,sql
sql="call search_info('','','','');"
oRecordset=OBJ_NEW('IDLdbRecordset',(*state.database).obj,sql=sql);

if(obj_valid(oRecordset)) then begin
  SET_TABLE_INFO,wDBtable,oRecordSet
endif 
;
(*state.database).obj->cleanup
;



info={wTopBase:wDBdialog,button_info:button_info,database:state.database,table:wDBtable,menu:wMenu,wButtonSearch:wDBdialogButtonSearch,state:state}


WIDGET_CONTROL,wDBdialog,/REALIZE
WIDGET_CONTROL,wDBdialog,set_uvalue=info,/no_copy
xmanager ,'db_load',wDBdialog,$
  Event_Handler='DB_DIALOG_EVENT'

if(*button_info).cancel eq 0 AND n_tags(*button_info) GT 2 then begin
  state.vol_HU_cube =ptr_new((*button_info).data)
  state.vol_HU_cube_resolution = (*button_info).resolution
  state.patient_age = (*button_info).age
  IF PTR_VALID((*button_info).path) THEN BEGIN
    IF STRLEN(*(*button_info).path) GT 1 THEN state.additionalDir = ptr_new(*((*button_info).path))
  ENDIF ELSE BEGIN
    IF STRLEN((*button_info).path) GT 1 THEN BEGIN
      state.additionalDir = ptr_new((*button_info).path)
    ENDIF ELSE BEGIN
      state.additionalDir = ptr_new((*(*button_info).form_info).Path)
    ENDELSE
  ENDELSE
endif
print ,(*button_info).cancel 
ptr_free,button_info


END


;添加新的病例，事件处理程序
PRO DB_ADD_EVENT,ev
WIDGET_CONTROL,ev.id,get_uvalue=uval
case uval[0] of
  'OK' :begin
  WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
 ; WIDGET_CONTROL,info.wPath,get_value=afn
  (*(info.button_info)).cancel=0;
 ; (*(info.button_info)).path=afn;
  
  select=WIDGET_INFO(info.wGen,/droplist_select)
  WIDGET_CONTROL,info.wGen,get_value= value
  str_Gender=STRTRIM(value[select],2);
  
  select=WIDGET_INFO(info.wModal,/droplist_select)
  WIDGET_CONTROL,info.wModal,get_value= value
  str_Modal=STRTRIM(value[select],2);
  ;date1 year
  select=WIDGET_INFO(info.wBirth.y,/droplist_select)
  WIDGET_CONTROL,info.wBirth.y,get_value= value
  str_y=STRTRIM(value[select],2);
  ;date1 month
  select=WIDGET_INFO(info.wBirth.m,/droplist_select)
  WIDGET_CONTROL,info.wBirth.m,get_value= value
  str_m=STRTRIM(value[select],2);
  ;date1 day
  select=WIDGET_INFO(info.wBirth.d,/droplist_select)
  WIDGET_CONTROL,info.wBirth.d,get_value= value
  str_d=STRTRIM(value[select],2);


  if str_y ne '' then  begin
    str_Birth=string(str_y,str_m,str_d,format='(i04,i02,i02)')
  endif else begin
    str_Birth=str_y
  endelse
  if strpos(str_Birth,'*') ne -1 then str_Birth=''
  
  WIDGET_CONTROL,info.wName,get_value=str_Name
  WIDGET_CONTROL,info.wHospital,get_value=str_Hospital
  WIDGET_CONTROL,info.wDoctor,get_value=str_Doctor
  WIDGET_CONTROL,info.wDescription,get_value=str_Description
  WIDGET_CONTROL,info.wNote,get_value=str_Note
  WIDGET_CONTROL,info.wPath,get_value=str_Path

  
  
  
  (*(info.button_info)).form_info=ptr_new({Name:STRJOIN(STRSPLIT(str_Name, /EXTRACT,"'"),"\'") $
  ,Gender:STRJOIN(STRSPLIT(str_Gender, /EXTRACT,"'"),"\'")  $
  ,Birth:STRJOIN(STRSPLIT(str_Birth, /EXTRACT,"'"),"\'")  $
  ,Modal:STRJOIN(STRSPLIT(str_Modal, /EXTRACT,"'"),"\'")  $
  ,Hospital:STRJOIN(STRSPLIT(str_Hospital, /EXTRACT,"'"),"\'") $
  ,Doctor:STRJOIN(STRSPLIT(str_Doctor, /EXTRACT,"'"),"\'")  $
  ,Description:STRJOIN(STRSPLIT(STRJOIN(str_Description,"\n"), /EXTRACT,"'"),"\'")  $
  ,Note:STRJOIN(STRSPLIT(STRJOIN(str_Note,"\n"), /EXTRACT,"'"),"\'") $
  ,Path:STRJOIN(STRSPLIT(str_Path, /EXTRACT,"'"),"\'")  $
  })

  WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
  WIDGET_CONTROL,ev.top,/destroy
end
'CANCEL':begin
      WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
      (*(info.button_info)).cancel=1;
      WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
      WIDGET_CONTROL,ev.top,/destroy
      end
  'LOADDATAPATH': BEGIN
        WIDGET_CONTROL,ev.top,get_uvalue=info,/no_copy
        
        afn = DIALOG_PICKFILE(PATH='./test_data/S0002_162200/',/MUST_EXIST, TITLE="file choose", /DIRECTORY)
        IF STRLEN(afn) EQ 0 THEN RETURN
        filearr=file_search(afn,'*.{nii,dcm}',count=num)
        IF num EQ 0 THEN BEGIN
          dcm_info = DIALOG_MESSAGE('Please sort images first!', /INFORMATION, TITLE='Suggestions')
          Separate_to_series, afn, afn_des
          afn = DIALOG_PICKFILE(PATH=afn_des,/MUST_EXIST, TITLE="file choose", /DIRECTORY)
        ENDIF
        
        filearr=file_search(afn,'*.{nii,dcm}',count=num)
        print,filearr[num/2]
      
        oImg = OBJ_NEW('IDLffDicom', filearr[num/2])
        p_name = STRMID(STRING(*(oImg->GetValue('0010'x,'0010'x,/NO_COPY))[0]), 0, 10)
        p_id = STRMID(STRING(*(oImg->GetValue('0010'x,'0020'x,/NO_COPY))[0]), 18)
        p_gender = STRING(*(oImg->GetValue('0010'x,'0040'x,/NO_COPY))[0])
        IF STRCMP(p_gender, 'M') THEN p_gender='male'
        p_dob = STRING(*(oImg->GetValue('0010'x,'0030'x,/NO_COPY))[0])
        p_age = STRING(*(oImg->GetValue('0010'x,'1010'x,/NO_COPY))[0])
        p_modal = STRING(*(oImg->GetValue('0008'x,'0060'x,/NO_COPY))[0])
        p_hospital = STRMID(STRING(*(oImg->GetValue('0008'x,'0080'x,/NO_COPY))[0]), 0, 20)
        p_doctor = STRMID(STRING(*(oImg->GetValue('0008'x,'0090'x,/NO_COPY))[0]), 0, 10)
        OBJ_DESTROY, oImg
        ;p_name = oImg->GetValue('0010,0010')
        ;p_id = oImg->getvalue('0010,0020')
        ;p_gender = oImg->GetValue('0010,0040')
        ;p_age = oImg->GetValue('0010,0030')
      ;  help,p_age
      ;  p_modal = oImg->GetValue('0008,0060')
      ;  p_hospital = oImg->GetValue('0008,0080')
        out = long(STRSPLIT(p_age,/EXTRACT))
      
        p_year = strMid(p_dob,0,4)
        p_month = strMid(p_dob,4,2)
        p_day = strMid(P_dob,6,2)
        print,p_name
        WIDGET_CONTROL,info.wGen,set_value = p_gender
        WIDGET_CONTROL,info.wName,set_value = p_name
        WIDGET_CONTROL,info.wBirth.Y,set_value = p_year
        WIDGET_CONTROL,info.wBirth.M,set_value = p_month
        WIDGET_CONTROL,info.wBirth.D,set_value = p_day
        WIDGET_CONTROL,info.wHospital,set_value = p_hospital
        WIDGET_CONTROL,info.wDoctor,set_value = p_doctor
        WIDGET_CONTROL,info.wModal,set_value = p_modal
        WIDGET_CONTROL,info.wPath,set_value=afn
      
        WIDGET_CONTROL,ev.top,set_uvalue=info,/no_copy
      
        END

ELSE: BEGIN


END
ENDCASE

END


;添加新的病例的主程序
;for ADD and update
PRO DB_ADD,state=state,case_id=case_id
button_info=ptr_New({cancel:1,form_info:ptr_new(),case_id:-1})

SCR_SIZE={x:400,y:400}
MARGIN=2
n_row=14
BORDER={x:SCR_SIZE.x-MARGIN*2,y:SCR_SIZE.y-MARGIN*2} 
subSize={x:BORDER.x,y:(SCR_SIZE.y-MARGIN*2)/n_row-MARGIN}

wDBDialogAdd=WIDGET_BASE(GROUP_LEADER=state.wTopBase,$
  /modal,$
  /column,$
  xsize=SCR_SIZE.x,$
  ysize=SCR_SIZE.y,$
  title='ADD')
wBorder=WIDGET_BASE(wDBDialogAdd,/align_center,/column,xsize=BORDER.x,ysize=BORDER.Y,/frame,/tab_mode)
wSubBase0=WIDGET_BASE(wBorder,xpad=0,/align_center,/column,xsize=subSize.x,ysize=subSize.y,/frame)
label_name=WIDGET_LABEL(wSubBase0,/align_center,value='information acquisition')

wSubBase8=WIDGET_BASE(wBorder,xpad=0,/align_center,/row,xsize=subSize.x,ysize=subSize.y,/frame)
wSubBase81=WIDGET_BASE(wSubBase8,/align_center,/column,xsize=subSize.x/6+margin*3,ysize=subSize.y)
wSubBase82=WIDGET_BASE(wSubBase8,/align_center,/row,xsize=subSize.x/6*5,ysize=subSize.y)
wSubBase821=WIDGET_BASE(wSubBase82,/align_center,/row,xsize=subSize.x/6*5,ysize=subSize.y)

label_data=WIDGET_LABEL(wSubBase81,/align_center,value='choose data:')
edit_data=WIDGET_TEXT(wSubBase821,/align_center,uname='note',uval='note',value = '',/EDITABLE,ysize=subSize.y,scr_xsize=(subSize.x/6*5-margin*2)*0.6,scr_ysize=subSize.y)
button_data=WIDGET_BUTTON(wSubBase821,/align_center,value='path',uval='LOADDATAPATH')

wSubBase1=WIDGET_BASE(wBorder,xpad=0,/row,xsize=subSize.x,ysize=subSize.y,/frame)      
  wSubBase11=WIDGET_BASE(wSubBase1,/align_center,/row,xsize=subSize.x/2,ysize=subSize.y)
  wSubBase12=WIDGET_BASE(wSubBase1,/align_center,/row,xsize=subSize.x/2,ysize=subSize.y)

wSubBase111=WIDGET_BASE(wSubBase11,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*1/3)
wSubBase112=WIDGET_BASE(wSubBase11,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*2/3)
wSubBase121=WIDGET_BASE(wSubBase12,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*1/3)
wSubBase122=WIDGET_BASE(wSubBase12,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*2/3)


label_name=WIDGET_LABEL(wSubBase111,/align_center,value='name:')
edit_name=WIDGET_TEXT(wSubBase112,/align_center,uname='ediname',uval='NAME',/EDITABLE,scr_ysize=subSize.y)
label_gender=WIDGET_LABEL(wSubBase121,/align_center,value='gender:')
drop_gender=WIDGET_Droplist(wSubBase122,/align_center,uval='gender',value=['               ','male','female'])

; 



wSubBase2=WIDGET_BASE(wBorder,xpad=0,/row,xsize=subSize.x,ysize=subSize.y,/frame)
wSubBase21=WIDGET_BASE(wSubBase2,/align_center,/row,xsize=subSize.x/2,ysize=subSize.y)
wSubBase22=WIDGET_BASE(wSubBase2,/align_center,/row,xsize=subSize.x/2,ysize=subSize.y)

wSubBase211=WIDGET_BASE(wSubBase21,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*1/3)
wSubBase212=WIDGET_BASE(wSubBase21,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*2/3)
wSubBase221=WIDGET_BASE(wSubBase22,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*1/3)
wSubBase222=WIDGET_BASE(wSubBase22,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*2/3)



label_birth=WIDGET_LABEL(wSubBase211,/align_center,value='birth:')
drop2_y=WIDGET_Droplist(wSubBase212,uval='BIRTH',value=[[''],[STRTRIM(string(indgen(1,100)+2030-100),2)]])
drop2_m=WIDGET_Droplist(wSubBase212,title='',uval='BIRTH',value=[[''],[STRTRIM(string(indgen(1,12)+1),2)]])
drop2_d=WIDGET_Droplist(wSubBase212,title='',uval='BIRTH',value=[[''],[STRTRIM(string(indgen(1,31)+1),2)]])


str_modal=['               ','CT','MR','DSA','MRI']
label_modaling=WIDGET_LABEL(wSubBase221,/align_center,value='modaling:')
drop_modaling=WIDGET_Droplist(wSubBase222,/align_center,uval='MODALING',value=str_modal)




wSubBase3=WIDGET_BASE(wBorder,xpad=0,/row,xsize=subSize.x,ysize=subSize.y,/frame)
wSubBase31=WIDGET_BASE(wSubBase3,/align_center,/row,xsize=subSize.x/2,scr_ysize=subSize.y)
wSubBase32=WIDGET_BASE(wSubBase3,/align_center,/row,xsize=subSize.x/2,scr_ysize=subSize.y)

wSubBase311=WIDGET_BASE(wSubBase31,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*1/3)
wSubBase312=WIDGET_BASE(wSubBase31,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*2/3)
wSubBase321=WIDGET_BASE(wSubBase32,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*1/3)
wSubBase322=WIDGET_BASE(wSubBase32,/row,/align_center,ysize=subSize.y,xsize=subSize.x/2*2/3)


label_hospital=WIDGET_LABEL(wSubBase311,/align_center,value='hospital:')
edit_hospital=WIDGET_TEXT(wSubBase312,/align_center,uname='edihos',uval='',/EDITABLE,scr_ysize=subSize.y)
label_doctor=WIDGET_LABEL(wSubBase321,/align_center,value='doctor:')
edit_doctor=WIDGET_TEXT(wSubBase322,/align_center,uname='edidoctor',uval='',/EDITABLE,scr_ysize=subSize.y)





;wSubBase4=WIDGET_BASE(wBorder,xpad=0,/align_center,/row,xsize=subSize.x,ysize=subSize.y*0.125)



wSubBase5=WIDGET_BASE(wBorder,xpad=0,/align_center,/row,xsize=subSize.x,ysize=subSize.y*6,/frame)
wSubBase51=WIDGET_BASE(wSubBase5,/align_center,/column,xsize=subSize.x/6+margin*4,ysize=subSize.y*6)
wSubBase52=WIDGET_BASE(wSubBase5,/align_center,/row,xsize=subSize.x/6*5,ysize=subSize.y*6)
label_description=WIDGET_LABEL(wSubBase51,/align_center,value='description:')
edit_description=WIDGET_TEXT(wSubBase52,/align_center,uname='description',uval='description',/EDITABLE,ysize=subSize.y*4,scr_xsize=subSize.x/6*5-margin*2,scr_ysize=subSize.y*6)



wSubBase6=WIDGET_BASE(wBorder,xpad=0,/align_center,/row,xsize=subSize.x,ysize=subSize.y*2,/frame)
wSubBase61=WIDGET_BASE(wSubBase6,/align_center,/column,xsize=subSize.x/6+margin*4,ysize=subSize.y*2)
wSubBase62=WIDGET_BASE(wSubBase6,/align_center,/row,xsize=subSize.x/6*5,ysize=subSize.y*2)
label_description=WIDGET_LABEL(wSubBase61,/align_center,value='note:')
edit_note=WIDGET_TEXT(wSubBase62,/align_center,uname='note',uval='note',/EDITABLE,ysize=subSize.y*4,scr_xsize=subSize.x/6*5-margin*2,scr_ysize=subSize.y*2)


wSubBase9=WIDGET_BASE(wBorder,xpad=0,/align_center,/row,xsize=subSize.x,ysize=subSize.y)
wSubBase91=WIDGET_BASE(wSubBase9,/align_center,/column,xsize=subSize.x/8*3,ysize=subSize.y)
wSubBase92=WIDGET_BASE(wSubBase9,/align_center,/row,xsize=subSize.x/4*1,ysize=subSize.y)
wDBbuttonGroup=WIDGET_BASE(wSubBase92,column=2,/align_center)
wDBAddButtonOk=WIDGET_BUTTON($
  wDBbuttonGroup,$
  /align_center,$
  value='ok',$
  uname='ok',$
  uval='OK'$
  ) 
wDBAddButtonCancel=WIDGET_BUTTON($
  wDBbuttonGroup,$
  /align_center,$
  value='CANCEL',$
  uname='CANCEL',$
  uval='CANCEL'$
  )





info={button_info:button_info $
  ,wName:edit_name $
  ,wGen:drop_gender $
  ,wBirth:{y:drop2_y,m:drop2_m,d:drop2_d}$
  ,wModal:drop_modaling $
  ,wHospital:edit_hospital $
  ,wDoctor:edit_doctor $
  ,wDescription:edit_description $
  ,wNote:edit_note $
  ,wPath:edit_data $
  ,wPathButton:button_data $
  }
if n_elements(case_id) then begin
    connect,state.database
    sql='use test;'
    (*state.database).obj->ExecuteSQL,sql
    sql='call id_search('+STRTRIM(string(case_id),2)+');'
    oRecordset=OBJ_NEW('IDLdbRecordset',(*state.database).obj,sql=sql);
    
    if(obj_valid(oRecordset)) then begin
      result=get_record(oRecordset)
      data=result.data
      ;individual info
      WIDGET_CONTROL,info.wName,set_value=data.NAME
      
      WIDGET_CONTROL,info.wGen,get_value=gen_list
      WIDGET_CONTROL,info.wGen,SET_DROPLIST_SELECT=WHERE(STRTRIM(string(gen_list),2) EQ data.GENDER)
      
      ;modaling
      WIDGET_CONTROL,info.wModal,get_value=modal_list
      WIDGET_CONTROL,info.wModal,SET_DROPLIST_SELECT=WHERE(STRTRIM(string(modal_list),2) EQ data.MODALING)
      
      ;date  
      WIDGET_CONTROL,info.wBirth.y,get_value=y_list
      WIDGET_CONTROL,info.wBirth.m,get_value=m_list
      WIDGET_CONTROL,info.wBirth.d,get_value=d_list
     
      IF STRPOS(data.y, '0') EQ 0 THEN data.y = STRMID(data.y, 1, 1)
      WIDGET_CONTROL,info.wBirth.y,SET_DROPLIST_SELECT=WHERE(STRTRIM(string(y_list),2) EQ data.y)
      IF STRPOS(data.m, '0') EQ 0 THEN data.m = STRMID(data.m, 1, 1)
      WIDGET_CONTROL,info.wBirth.m,SET_DROPLIST_SELECT=WHERE(STRTRIM(string(m_list),2) EQ data.m)
      IF STRPOS(data.d, '0') EQ 0 THEN data.d = STRMID(data.d, 1, 1)
      WIDGET_CONTROL,info.wBirth.d,SET_DROPLIST_SELECT=WHERE(STRTRIM(string(d_list),2) EQ data.d)
      
      ;hospital info
      WIDGET_CONTROL,info.wHospital,set_value=data.HOSPITAL  
      WIDGET_CONTROL,info.wDoctor,set_value=data.Doctor  
      
      ;Description and note
      WIDGET_CONTROL,info.wDescription,set_value=STRSPLIT(data.Description,/EXTRACT,string(10b))
      WIDGET_CONTROL,info.wNote,set_value=STRSPLIT(data.Note,/EXTRACT,string(10b))
      
      WIDGET_CONTROL,info.wPath,SENSITIVE=0
      WIDGET_CONTROL,info.wPathButton,SENSITIVE=0
        
    endif
    (*state.database).obj->cleanup
  
endif
  

WIDGET_CONTROL,wDBDialogAdd,/REALIZE
WIDGET_CONTROL,wDBDialogAdd,set_uvalue=info,/no_copy
xmanager ,'db_add',wDBDialogAdd,$
  Event_Handler='DB_ADD_EVENT'


if((*button_info).cancel ne 1) then begin
 if n_elements(case_id) then begin
   if case_id ne '' then begin
     form=*((*button_info).form_info)
     connect,state.database
     sql='use test;'
    (*state.database).obj->ExecuteSQL,sql
    path_temp = form.Path
    WHILE (((I = STRPOS(path_temp, '\'))) NE -1) DO STRPUT, path_temp, '*', I
    sql="call  update_info"               + $
      " ( "+STRTRIM(string(case_id),2)    + $
      " ,'"+form.Name          + $
      "','"+form.Gender        + $
      "','"+form.Birth         + $
      "', "+"null"             + $
      " ,'"+form.Modal         + $
      "','"+form.Hospital      + $
      "','"+form.Doctor        + $
      "','"+form.Description   + $
      "','"+form.Note          + $
      "','"+path_temp          + $
      "');"
    (*state.database).obj->ExecuteSQL,sql
    (*state.database).obj->cleanup
   endif 
 endif else begin
      image_read, (*(*button_info).form_info).Path, vol_HU_cube, vol_HU_cube_resolution, patient_age, p_modal = p_modal;, vol_HU_iso, vol_HU_iso_db, vol_HU_iso_resolution
      ;确保vol_HU_iso不为空
      IF((size(vol_HU_cube))[0] GT 2) THEN BEGIN
          ptr_free,state.vol_HU_cube
          state.vol_HU_cube = ptr_new(vol_HU_cube)
          state.vol_HU_cube_resolution = vol_HU_cube_resolution
          state.patient_age = patient_age
          state.additionalDir = ptr_new((*(*button_info).form_info).Path) ;GJ, 2018/11/18
          print, 'patient age 4 = ', state.patient_age
          print, 'database resolution', state.vol_HU_cube_resolution
        form=*((*button_info).form_info)
        connect,state.database
        sql='use test;'
        (*state.database).obj->ExecuteSQL,sql
        path_temp = form.Path
        WHILE (((I = STRPOS(path_temp, '\'))) NE -1) DO STRPUT, path_temp, '*', I
        sql="call  insert_info('"+form.Name + $
          "','"+form.Gender        + $
          "','"+form.Birth         + $
          "', "+"null"             + $
          " , "+"null"             + $
          " ,'"+form.Modal         + $
          "','"+form.Hospital      + $
          "', "+"now()"            + $
          " ,'"+form.Doctor        + $
          "','"+form.Description   + $
          "','"+form.Note          + $
          "','"+path_temp          + $
          "');"     
        oRecordset=OBJ_NEW('IDLdbRecordset',(*state.database).obj,sql=sql);
        IF oRecordset->MOVECURSOR(/first) EQ 1 THEN BEGIN
             oRecordset->GETPROPERTY,field_info = fieldinfo
             type=fieldinfo.TYPE_NAME
             id=TYPE_CONVERSION(oRecordset->GETFIELD(0),TYPE(0))
             ;(*state.database).obj->cleanup
             ;connect,state.database
             ttable='data_table'
             oRecordset=OBJ_NEW('IDLdbRecordset',(*state.database).obj,table=ttable);
             if (obj_valid(oRecordset)) then begin
                sql='set global max_allowed_packet=268435456;'
                 (*state.database).obj->ExecuteSQL,sql
                 s=size(vol_HU_cube)
                 ;only save a y-z mid-section image
                 oRecordset->AddRecord,id,BYTSCL(DIST(100)),100,100,vol_HU_cube[*,s[2]/2,*],s[1],1,s[3]
                ;(*state.database).obj->ExecuteSQL,'call insert_data()'
             endif 
        ENDIF
        (*state.database).obj->cleanup
      ENDIF    
   endelse
endif
ptr_free,button_info
END


