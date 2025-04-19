
;
;
PRO USING_ODBC_EXCEL
  filename = FILE_DIRNAME(ROUTINE_FILEPATH('Using_ODBC_EXCEL'))+'\SatelliteInformation.xlsx'

  IF DB_EXISTS() EQ 0 THEN BEGIN
    msg = DIALOG_MESSAGE('ODBC!',/Error)
    RETURN
  ENDIF

  oDatabase = OBJ_NEW('IDLdbDatabase')

  sources = oDatabase->GETDATASOURCES()
  index = WHERE(sources.DATASOURCE EQ 'Excel Files',count)
  IF count EQ 0 THEN BEGIN
    msg = DIALOG_MESSAGE('ODBC',/Error)
    OBJ_DESTROY,oDatabase
    RETURN
  ENDIF

  IF ~FILE_TEST(filename) THEN BEGIN
    msg = DIALOG_MESSAGE('ODBC',/Error)
    OBJ_DESTROY,oDatabase
    RETURN
  ENDIF

  oDatabase->CONNECT,DATASOURCE='Excel Files;DBQ='+filename

  oDatabase->GETPROPERTY,IS_CONNECTED = connectStat
  IF connectStat EQ 0 THEN BEGIN
    msg = DIALOG_MESSAGE('ODBC',/Error)
    OBJ_DESTROY,oDatabase
    RETURN
  ENDIF
  

  tables = oDatabase->GETTABLES()
  nTables = N_ELEMENTS(tables)
  FOR i=0, nTables-1 DO BEGIN
    ; 
    tname = '[' + tables[i].NAME + ']'
    PRINT, 'table name', tname
    
    oRecordset = OBJ_NEW('IDLdbRecordset',oDatabase,table=tname);, SQL=sqlstr)
    ; 
    oRecordset->GETPROPERTY,field_info = fieldinfo
    NFileds = N_ELEMENTS(fieldinfo)
    ;
    IF oRecordset->MOVECURSOR(/first) THEN BEGIN
      FOR j=0, NFileds-1 DO BEGIN
        Value = oRecordset->GETFIELD(j)
        PRINT, 'Talbe: ' + (fieldinfo.TABLE_NAME)[j] + ', ' + $
          'Filed Name: ' + (fieldinfo.FIELD_NAME)[j] + ', ' + $
          'Value: ', Value
      ENDFOR
      WHILE oRecordset->MOVECURSOR(/next) DO BEGIN
        FOR j=0, NFileds-1 DO BEGIN
          Value = oRecordset->GETFIELD(j)
          PRINT, 'Talbe: ' + (fieldinfo.TABLE_NAME)[j] + ', ' + $
            'Filed Name: ' + (fieldinfo.FIELD_NAME)[j] + ', ' + $
            'Value: ', Value
        ENDFOR
      ENDWHILE
    ENDIF ELSE BEGIN ; 
      msg = DIALOG_MESSAGE('ODBC',/Infor)
      OBJ_DESTROY, oRecordset
    ENDELSE
    ;
    OBJ_DESTROY, oRecordset
  ENDFOR
  OBJ_DESTROY,oDatabase
END