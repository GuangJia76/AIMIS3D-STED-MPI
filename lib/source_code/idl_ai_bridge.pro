;
;+
;Name: search_filename.pro
;
;PURPOSE: Search a filename in a filename array
;
;MODULES CONTAINED:
;
;CATEGORY:
;
;CALLING SEQUENCE: status=search_filename(filenamearray, filename, flag, tError)
;
;OUTPUTS: flag is 1 for filename in filenamearray, 0 for filename not in filenamearray.
;     status is 0 for successfully finish search, 1 for fail to search.
;
;COMMON BLOCKS:
;
;SIDE EFFECTS:
;
;
;
FUNCTION search_filename, filenamearray, filename, flag, tError

  tError=''
  ;; first turn off math error output
  except_old=!except
  !except=0 ;; never report math errors immediately after the statement which caused them, copied from slice_batch.pro

  ON_IOERROR, error

  ;;BEGIN OF CATCH ERROR HANDLER
  CATCH, Error_status
  IF Error_status NE 0 THEN BEGIN
    ;;Handle the error
    tError=!ERR_STRING
    CATCH, /CANCEL;; needed??
    GOTO, rerror
  ENDIF
  ;; END OF CATCH ERROR HANDLER

  IF N_ELEMENTS(filenamearray) LE 0 THEN BEGIN
    tError="No filename array can be used to serach in?"
    GOTO, rerror
  ENDIF

  ;flag used to show whether filename in filenamearray or not
  flag=0
  FOR i=0L,N_ELEMENTS(filenamearray)-1 DO BEGIN
    IF STRCMP(filenamearray[i], filename, /FOLD_CASE) THEN BREAK
  ENDFOR

  IF i LT N_ELEMENTS(filenamearray) THEN flag=1

  ;; return with status OK
  RETURN,0
  IF 0 EQ 1 THEN BEGIN
    ;; return with error
    error: tError=!ERROR_STATE.MSG
    rerror:PRINT, tError
    RETURN, 1
  ENDIF
END



;
;+
;Name: compare_filenamearrays_find_diff.pro
;
;PURPOSE: Compare two filename arrays and find out one element in second array which is not in the first array
;
;MODULES CONTAINED:
;
;CATEGORY:
;
;CALLING SEQUENCE: status=search_filename(filenamearray1, filenamearray2, index_fnarray2, tError)
;
;OUTPUTS: index_fnarray2 is the index with which the element in filenamearray2 is not included in filenamearray1.
;     If all elements in filenamearray2 are included in filenamearray1, index_fnarray2 is equal to -1.
;     status is 0 for successfully finish search, 1 for fail to search.
;
; ?? I don't understand that. Is this the index of the first element
; in array 1 which has no match in array 2??
;
;
;
;COMMON BLOCKS:
;
;SIDE EFFECTS:
;
;
;
FUNCTION compare_filenamearrays_find_diff, filenamearray1, filenamearray2, index_fnarray2, tError

  tError=''
  ON_IOERROR, error

  ;;BEGIN OF CATCH ERROR HANDLER
  CATCH, Error_status
  IF Error_status NE 0 THEN BEGIN
    ;;Handle the error
    tError=!ERR_STRING
    CATCH, /CANCEL;; needed??
    GOTO, rerror
  ENDIF
  ;; END OF CATCH ERROR HANDLER

  ;initialization
  index_fnarray2=-1

  IF (N_ELEMENTS(filenamearray1) LE 0) OR (N_ELEMENTS(filenamearray2) LE 0) THEN BEGIN
    tError="No filename arrays can be used to compare?"
    GOTO, rerror
  ENDIF

  FOR i=0L, N_ELEMENTS(filenamearray2)-1 DO BEGIN
    status=search_filename(filenamearray1, filenamearray2[i], flag, tError)
    IF status THEN GOTO, rerror
    IF flag EQ 0 THEN BREAK
  ENDFOR

  IF i LT N_ELEMENTS(filenamearray2) THEN index_fnarray2=i

  ;; return with status OK
  RETURN,0

  ;; return with error
  error: tError=!ERROR_STATE.MSG
  rerror:PRINT, tError
  RETURN, 1
END



;
;+
;Name: process_batches.pro
;
;PURPOSE: Process every batch in a directory, when processing, you can add a new .batch to this directory and this function
;         will process it.  You can also delete a .batch you don't want to be processed as long as you know it has not been processed
;
;MODULES CONTAINED:
;
;CATEGORY:
;
;CALLING SEQUENCE: status=process_batches(tError)
;
;OUTPUTS: status is 0 for successfully finish all batches, 1 for fail to process batches.
;
;
;COMMON BLOCKS:
;
;SIDE EFFECTS:
;
;MODIFICATION HISTORY:
;��,GJ, Feb. 28, 2003, Creation;
;��,GJ, March. 30, 2003, Make it neat and divide it into several functions
;
;
;
;
FUNCTION process_batches,aif_parent, tError
  tError=''
  ;; first turn off math error output
  except_old=!except
  !except=0 ;; never report math errors immediately after the statement which caused them, copied from slice_batch.pro

  ON_IOERROR, error

  ;;BEGIN OF CATCH ERROR HANDLER
  CATCH, Error_status
  IF Error_status NE 0 THEN BEGIN
    ;;Handle the error
    tError=!ERR_STRING
    CATCH, /CANCEL;; needed??
    GOTO, rerror
  ENDIF
  ;; END OF CATCH ERROR HANDLER

  ;Select the directory where the .batch files are saved
  batch_dir=DIALOG_PICKFILE(/DIRECTORY, /MUST_EXIST, $
    TITLE='select batch directory')

  IF STRLEN(batch_dir) EQ 0 THEN BEGIN
    tError='No batch processed'
    GOTO, finish
  ENDIF
  ;Always search .batch file in a directory
  WHILE (1 NE 0) DO BEGIN

    ;; find out all .batch files in this directory, put this
    ;; inside the loop because someone may add or delete any .batch files
    fdecomp, batch_dir, disk, dir, name, qual, version
    fcomp, batch_dir_all_file, disk, dir, "*", "batch", version
    batch_fn_list=FILE_SEARCH(batch_dir_all_file, COUNT=number_batch_fn,FOLD_CASE=1)
    IF number_batch_fn EQ 0 THEN GOTO, finish

    ;; Compare whether every file in batch_fn_list is
    ;; inside the batch_fn_list_processed
    IF number_batch_fn NE 0 THEN BEGIN
      number_batch_fn_processed=N_ELEMENTS(batch_fn_list_processed)
      ;;No batch file was processed, just take the first one in batch_fn_list
      IF number_batch_fn_processed EQ 0 THEN BEGIN
        batch_fn_list_processed=[batch_fn_list[0]]
        number_batch_fn_processed=1
        ;;also create a .log file, which save the infomation
        fcomp, batch_log_fn, disk, dir, "batch", "log", version
        OPENW,unit,batch_log_fn, /APPEND,/get_lun
        FREE_LUN, unit
      ENDIF ELSE BEGIN
        status=compare_filenamearrays_find_diff(batch_fn_list_processed, batch_fn_list, diff_index, tError)
        IF status THEN GOTO, rerror
        IF diff_index EQ -1 THEN GOTO, finish
        batch_fn_list_processed=[batch_fn_list_processed,$
          batch_fn_list[diff_index]]
        number_batch_fn_processed=number_batch_fn_processed+1L
      ENDELSE
    ENDIF

    ;;
    ;;Copy the newly added .batch file to a temp file and do
    ;;slice_batch, in case someone delete the newly added .batch file
    batch_filename=batch_fn_list_processed[number_batch_fn_processed-1]
    fdecomp, batch_filename, disk, dir, name, qual, version
    fcomp, copy_batch_fn, disk, dir, "On_Processing_"+name, qual, version
    file_copy, batch_filename, copy_batch_fn
    status=read_batch(batch_filename,dyn_fn, roi_fn, model, output_dir, tError)
    IF status THEN BEGIN
      GOTO, error
    ENDIF

    ;;Start processing the batch, now write info to the .log file
    OPENU,unit,batch_log_fn, /APPEND,/get_lun
    PRINTF, unit, name+"."+qual+" started at " + systime()+"!"
    FREE_LUN, unit
    ;; ERROR HANDLER !!!!


    ;;Process the batch
    fdecomp,roi_fn, disk, dir, name, qual, version
    IF STRCMP(qual, "roi", /FOLD_CASE) THEN BEGIN
;      s=slice_batch_ds(dyn_fn,model, output_dir,ROI_ds=roi_fn,tError)
    ENDIF ELSE BEGIN
;      s=slice_batch(aif_parent, dyn_fn,model, output_dir,MASK_FN=roi_fn,tError)
    ENDELSE
    ;;s=1
    IF s THEN BEGIN
      GOTO, error
    ENDIF

    ;;  WAIT, 10 ;; for GUANGs testing

    ;;Finish processing the batch, now write info to the .log file
    OPENU,unit,batch_log_fn, /APPEND,/get_lun
    PRINTF, unit, name+"."+qual+" finished at " + systime()+"!"
    FREE_LUN, unit
    file_delete, copy_batch_fn

    ;; and an error handler
    IF 0 EQ 1 THEN BEGIN
      error: tError='No batch processed'
      rerror: OPENU,unit,batch_log_fn, /APPEND,/get_lun
      PRINTF, unit, name+"."+qual+" " + terror +" "+ systime()+"? Fail!"
      FREE_LUN, unit
      IF STRLEN(FILE_SEARCH(copy_batch_fn)) EQ 0 THEN BEGIN
        terror=copy_batch_fn+"does not exist!"
      ENDIF ELSE BEGIN
        file_delete, copy_batch_fn
      ENDELSE
    ENDIF

  ENDWHILE

  finish:
  RETURN, (0) ;; successfull completion

END

;
;Created by GJ, 2021/2/16
;
;name: idl_ai_bridge
;function: prepare images for python AI algorithm, load output mask into IDL
PRO idl_ai_bridge


END