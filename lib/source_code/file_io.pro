;+
; NAME: file_io.pro
;
;
;
; PURPOSE: Collection of routines for image input/output for
;          dynamik - image postprocessing program
;
; MODULES CONTAINED:
;
;
;
; CATEGORY:
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
; MODIFICATION HISTORY:
;  Sep. 5, 2002 Creation (K. Baudendistel, M. Knopp)
;-

;; check for a directory and create it, if it does not exist
;;
FUNCTION dircreate,path,tError
 tError=''
 IF FILE_TEST(path) EQ 1 THEN BEGIN
     IF FILE_TEST(path,/DIRECTORY) NE 1 THEN BEGIN
         tError='file <' + path + '> exists and is NOT a directory'
         return, 1
     ENDIF ELSE BEGIN
         return, 0;; EXIT SUCCESS
     ENDELSE
 ENDIF ELSE BEGIN ;; try to create it
     FILE_MKDIR,path ;; as idl help does not mention any return status, check for error
     IF FILE_TEST(path,/DIRECTORY) NE 1 THEN BEGIN
         tError='file <' + path + '> could not be created'
         return, 1
     ENDIF ELSE BEGIN
         return, 0;; EXIT SUCCESS
     ENDELSE
ENDELSE
END





;COMMON FILE_INPUT, FN_ARR_IN, PATH_IN, REF_FN_IN, EXT_IN
; THE FOLLOWING IS A ROUTINE FROM
; http://idlastro.gsfc.nasa.gov/ftp/v53/
function gettok,st,char
;+
; NAME:
;   GETTOK
; PURPOSE:
;   Retrieve the first part of a (vector) string up to a specified character
; EXPLANATION:
;   GET TOKen - Retrieve first part of string until the character char
;   is encountered.    IDL V5.3 or later!
;
; CALLING SEQUENCE:
;   token = gettok( st, char )
;
; INPUT:
;   char - character separating tokens, scalar string
;
; INPUT-OUTPUT:
;   st - string to get token from (on output token is removed),
;            scalar or vector
;
; OUTPUT:
;   token - extracted string value is returned, same dimensions as st
;
; EXAMPLE:
;   If ST is ['abc=999','x=3.4234'] then gettok(ST,'=') would return
;   ['abc','x'] and ST would be left as ['999','3.4234']
;
; PROCEDURE CALLS:
;       REPCHR
; HISTORY
;   version 1  by D. Lindler APR,86
;   Remove leading blanks    W. Landsman (from JKF)    Aug. 1991
;   Converted to IDL V5.0   W. Landsman   September 1997
;       V5.3 version, accept vector input   W. Landsman February 2000
;       Slightly faster implementation  W. Landsman   February 2001
;-
;----------------------------------------------------------------------
  On_error,2                           ;Return to caller

  st = strtrim(st,1)              ;Remove leading blanks

; if char is a blank treat tabs as blanks

  tab = string(9b)
  if max(strpos(st,tab)) GE 0 then st = repchr(st,tab,' ')

  token = st

; find character in string

  pos = strpos(st,char)
  test = pos EQ -1
  bad = where(test, Nbad)
  if Nbad GT 0 then st[bad] = ''

; extract token

 good = where(1b-test, Ngood)
 if Ngood GT 0 then begin
    stg = st[good]
    pos = reform( pos[good], 1, Ngood )
    token[good] = strmid(stg,0,pos)
    st[good] = strmid(stg,pos+1)
 endif

;  Return the result.

 return,token
 end








; THE FOLLOWING IS A ROUTINE FROM
; file '/net/www/deutsch/idl/idllib/astron/contrib/landsman/allpro/fdecomp.pro'
pro fdecomp, filename, disk, dir, name, qual, version, OSfamily = osfamily
;+
; NAME:
;   FDECOMP
; PURPOSE:
;   Routine to decompose a file name for any operating system
;
; CALLING SEQUENCE:
;   FDECOMP, filename, disk, dir, name, qual, version, [OSFamily = ]
;
; INPUT:
;   filename - string file name, scalar
;
; OUTPUTS:
;   All the output parameters are scalar strings
;   disk - disk name, always '' on a Unix machine, scalar string
;   dir - directory name, scalar string
;   name - file name, scalar string
;   qual - qualifier, set equal to the characters beyond the last "."
;   version - version number, always '' on a non-VMS machine, scalar string
;
; OPTIONAL INPUT KEYWORD:
;   OSFamily - one of the four scalar strings specifying the operating
;     system:  'vms','windows','MacOS' or 'unix'.    If not supplied,
;     then !VERSION.OS_FAMILY is used to determine the OS.
; EXAMPLES:
;   Consider the following file names
;
;   Unix:    file = '/rsi/idl40/avg.pro'
;   VMS:     file = '$1$dua5:[rsi.idl40]avg.pro;3
;   Mac:     file = 'Macintosh HD:Programs:avg.pro'
;   Windows: file =  'd:\rsi\idl40\avg.pro'
;
;   then IDL> FDECOMP,  file, disk, dir, name, qual, version
;   will return the following
;
;       Disk             Dir          Name        Qual     Version
;   Unix:      ''            '/rsi/idl40/'  'avg'       'pro'       ''
;   VMS:     '$1$dua5'       '[RSI.IDL40]'  'avg'       'pro'       '3'
;   Mac:     'Macintosh HD'  ':Programs:'   'avg'       'pro'       ''
;   Windows:    'd:'         \rsi\idl40\    'avg'       'pro'       ''
;
; NOTES:
;   (1) All tokens are removed between
;     1) name and qual  (i.e period is removed)
;     2) qual and ver   (i.e. VMS semicolon is removed)
;   (2) On VMS the filenames "MOTD" and "MOTD." are distinguished by the
;       fact that qual = '' for the former and qual = ' ' for the latter.
;
; ROUTINES CALLED:
;   Function GETTOK()
; HISTORY
;   version 1  D. Lindler  Oct 1986
;   Include VMS DECNET machine name in disk    W. Landsman  HSTX  Feb. 94
;   Converted to Mac IDL, I. Freedman HSTX March 1994
;
;   Converted to IDL V5.0   W. Landsman   September 1997
;-
;   Converted to IDL V5.0   W. Landsman   September 1997
;--------------------------------------------------------
;
  On_error,2                            ;Return to caller

  if N_params() LT 2 then begin
     print, 'Syntax - FDECOMP, filename, disk, [dir, name, qual, ver ] '
     return
  endif

; Find out what machine you're on, and take appropriate action.
 if not keyword_set(OSFAMILY) then osfamily = !VERSION.OS_FAMILY

 case OSFAMILY of

  "MacOS": begin

; disk name is all characters up to the first colon
; directory is string of folders
; file name+qualifier is all characters after the last colon
; version   is null string

  st = filename
  if strpos(st,':') GE 0 then disk = gettok(st,':')  else disk = ''

     dir = ':' & tok = ''
     REPEAT BEGIN
    oldtok = tok
    tok = gettok(st,':')
    dir = dir + oldtok + ':'
     ENDREP UNTIL tok EQ ''

       dir = strmid(dir,1,strpos(dir,oldtok)-1)

     fname = oldtok & qual = ''
     pos = strpos(fname,'.')
     if pos GE 0 then begin
       name = gettok(fname,'.')
       qual   = fname
     endif

    version = ''

    end

 "vms":  begin                     ; begin VMS version

    st = filename

; get disk

    nodepos = strpos(st,'::')          ; Node name included in directory?
    if nodepos GE 0 then begin
    disk = strmid(st,0,nodepos+2)
    st = strmid(st,nodepos+2, 999 )
    endif else disk = ''
    if strpos(st,':') GE 0 then disk = disk + gettok(st,':') + ':' else $
                                disk = disk + ''

; get dir

    if strpos( st, ']' ) GE 0 then dir = gettok( st, ']' ) + ']' else dir=''
    if strpos( st, ']' ) GE 0 then dir = dir + gettok( st, ']' ) + ']'

; get name

    sv_name = st
    name = gettok(st,'.')

; get qualifier

    if (name + '.') EQ sv_name then qual = ' ' else $
    qual = gettok(st,';')

; get version

    version = st

  end   ;  end VMS version

 "Windows": begin

     st = filename
     pos = strpos( st, ':')         ; DOS diskdrive (i.e. c:)
     if (pos gt 0) then disk = gettok(st,':') + ':' else disk=''
; ====== INSERTED BY KB, SEP 2002: DEALING WITH WIN FOR WORKGROUP FILENAMES ========
     pos = strpos( st, '\\')
     if (pos NE -1) THEN BEGIN; filename starts with a '\\'
     cp=strmid(st,2) ; computer + pathname without leading backslashes
     disk='\\'+gettok(cp,'\')
     st='\'+cp ; correct for gettok removing the backslash from the path
     ENDIF
; ================ END OF INSERTION ================================================

;  Search the path name (i.e. \dos\idl\) and locate all backslashes

     lpos = -1  ; directory position path (i.e. \dos\idl\)
     pos = -1
     repeat begin
    pos = strpos(st, '\',pos+1)
    if (pos GE 0) then lpos = pos
     endrep until pos lt 0

     ;  Parse off the directory path

     if lpos ge 0 then begin
    dir = strmid(st, 0, lpos+1)
    len = strlen(st)
    if lpos eq (len-1) then $
       st = '' else st = strmid(st,lpos+1,len-lpos-1)
     endif else dir=''

; get Windows name and qualifier (extension)...qual is optional

     lpos=-1
     repeat begin
    pos = strpos(st,'.',pos+1)
        if (pos ge 0) then lpos = pos
    endrep until pos lt 0

    ; Parse name and qual (if a qual was found )

     if lpos ge 0 then begin
    len = strlen(st)
    name = strmid(st,0,lpos)
        qual = strmid(st,lpos+1,len-lpos-1)
     endif else begin
    name = st
    qual = ''
     endelse

     version = ''     ; no version numbers in Windows
     end

 ELSE: begin

    st = filename

; get disk

    disk = ''

; get dir

    lpos = -1
    pos = -1
    repeat begin
        pos = strpos(st, '/', pos+1)
        if (pos GE 0) then lpos = pos
    endrep until pos LT 0

    if lpos GE 0 then begin
        dir = strmid(st, 0, lpos+1)
        len = strlen(st)
        if lpos eq (len-1) then st = '' else $
                                    st = strmid(st,lpos+1,len-lpos-1)
    endif else dir = ''

; get name and qual

    pos = -1
    lpos = -1
    repeat begin
             pos = strpos(st,'.',pos+1)
             if (pos GE 0) then lpos = pos
    endrep until pos LT 0

    if lpos GE 0 then begin
             len = strlen(st)
             name = strmid(st,0,lpos)
             qual = strmid(st,lpos+1,len-lpos-1)
     endif else begin
         name = st
         qual = ''
     endelse

    version = ''

 end

ENDCASE        ; end OTHER version


  return
  end

;;; AND THE INVERSE FUNCTION WRITTEN BY KB
;; COMPOSE A FILENAME FROM THE THINGS RETURNED BY THE CODE ABOVE
pro fcomp, filename, disk, dir, name, qual, version, OSfamily = osfamily
;+
; NAME:
;   FCOMP
; PURPOSE:
;   Routine to compose a file name for any operating system
;
; CALLING SEQUENCE:
;   FCOMP, filename, disk, dir, name, qual, version, [OSFamily = ]
;
;
; INPUTS:
;   All the input parameters are scalar strings
;   disk - disk name, always '' on a Unix machine, scalar string
;   dir - directory name, scalar string
;   name - file name, scalar string
;   qual - qualifier, set equal to the characters beyond the last "."
;   version - version number, always '' on a non-VMS machine, scalar string
;
; OUTPUT:
;   filename - string file name, scalar
;
; OPTIONAL INPUT KEYWORD:
;   OSFamily - one of the four scalar strings specifying the operating
;     system:  'vms','windows','MacOS' or 'unix'.    If not supplied,
;     then !VERSION.OS_FAMILY is used to determine the OS.
; EXAMPLES:
; << STILL TO BE MODIFIED >>
; HISTORY
;--------------------------------------------------------
;
  On_error,2                            ;Return to caller

  if N_params() LT 2 then begin
     print, 'Syntax - FCOMP, filename, disk, [dir, name, qual, ver ] '
     return
  endif

; Find out what machine you're on, and take appropriate action.
 if not keyword_set(OSFAMILY) then osfamily = !VERSION.OS_FAMILY

 case OSFAMILY of

  "MacOS": begin

    end

 "vms":  begin                     ; begin VMS version


  end   ;  end VMS version

 "Windows": begin

;  put the things together
    filename=disk+dir+name+'.'+qual
 end ;   end Windows version

 "unix": begin

;  put the things together
     filename=disk+dir+name+'.'+qual

     end



 ELSE: begin


 end

ENDCASE        ; end OTHER version

RETURN
END


;+
; NAME:FUNCTION cat_fpaths,core,add
;
;
;
; PURPOSE: OS-system dependent creation of pathnames by concatenating
;          an old valid path with a new one. The old path must not be
;          empty
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
PRO cat_fpaths,result,core,add,OSFAMILY=osfamily
; Find out what machine you're on, and take appropriate action.
result='';; empty string
if not keyword_set(OSFAMILY) then osfamily = !VERSION.OS_FAMILY

case OSFAMILY of

    "MacOS": begin

    end

    "vms":  begin               ; begin VMS version


    end                         ;  end VMS version

    "Windows": begin
;  put the things together
        filename=core + '\' + add + '\'
        ;; windows does not like two or more backslashes
        ;; in the middle of a path, so replace two backslashes
        ;; by one. Do not start at the beginning, as this may
        ;; be a WfW share filename



        f   = STREGEX(STRMID(filename,2),'\\{2}') ;; search for '\\'

        REPEAT BEGIN
            IF f NE -1 THEN BEGIN ;; eliminate the '\\'
                result = strmid(filename,0,f+2)+strmid(filename,f+3)
                filename=result
            ENDIF
            f   = STREGEX(STRMID(filename,2),'\\{2}')
        ENDREP UNTIL (f eq -1)
    END                         ;   end Windows version

    "unix": begin
;; unix is simple: (/ are redundant, so 2 or more in succession don't matter)
        result = core + '/' + add + '/'
    end

    ELSE: begin
    end

ENDCASE        ; end OTHER version

RETURN
END



;+
; NAME:FUNCTION cat_fpfn,path, file, qual
;
;
;
; PURPOSE: OS-system dependent creation of filenames by adding path,
; file, and qual
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
;;;            patient[0].sl_pos = (*(patient[0].pMorphInfo))[patient[0].cur_slice,0,0].struc_image_info.sl_pos
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
PRO cat_fpfn,result,path,file,qual,OSFAMILY=osfamily
; Find out what machine you're on, and take appropriate action.
result='';; empty string
if not keyword_set(OSFAMILY) then osfamily = !VERSION.OS_FAMILY

case OSFAMILY of

    "MacOS": begin

    end

    "vms":  begin               ; begin VMS version


    end                         ;  end VMS version

    "Windows": begin
        result = path + file + '.' + qual
    END                         ;   end Windows version
    "unix": begin
        result = path + file + '.' + qual
    end                         ;   end unix version
    ELSE: begin
    end

ENDCASE        ; end OTHER version

RETURN
END








;+
; NAME: fn_create_numa3_internal
;
;
;
; PURPOSE: create a filename in SIEMENS NUMARIS3/internal
;          file format
;          onr-stu-image.ima e.g. 123-2-23
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;  status = create_fn_numa3_internal(fn_out,           $
;           onr,firststudy,firstimage,slicesperstudy,  $
;           timepoint,slice)
;
;
;
;
; INPUTS:
;        onr:             SIEMENS NUMA3 ordernumber
;        firstudy:        first SIEMENS NUMA3 studynumber of series
;        firstimage:      first SIEMENS NUMA3 imagenumber of series
;        slicesperstudy:  number of slices per study
;        timepoint       :timeframe# (0..#of timepoints -1)
;        slice           :slice#     (0..#of slices-1)
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;        RETURNS ERROR CODE as returnvalue
;
;        OUTPUT-Arguments:
;        fn_out:         SIEMENS NUMA3 filename (without extension)
;                        AND without PATH!!!
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


FUNCTION create_fn_numa3_internal,fn_out,$
                                  onr,firststudy,firstimage,imasperstudy,$
                                  timepoint,slice

study=firststudy+timepoint
image=firstimage+imasperstudy*timepoint+slice

fn_out=string(onr,study,image,FORMAT='(%"%d-%d-%d")')

RETURN, (0)
END



;+
; NAME: fn_create_hst
;
;
;
; PURPOSE: create a filename in handnummer slice timepoint format
;          e.g. h0230103
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;  status = create_hst(fn_out,           $
;                      handnr,timepoint,slice)
;
; INPUTS:
;        handnrr         :patient identification number
;        timepoint       :timeframe# (0..#of timepoints -1)
;        slice           :slice#     (0..#of slices-1)
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;        RETURNS ERROR CODE as returnvalue
;
;        OUTPUT-Arguments:
;        fn_out:         filename (without extension)
;                        AND without PATH!!!
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


FUNCTION init_fn_hst, handnr,timepoint, slice

fn_out=STRING(handnr,timepoint,slice,FORMAT='(%"h%03d%02d%02d")')

RETURN, (0)
END

;


;+
; NAME: init_fn_array_numa3_internal
;
;
;
; PURPOSE: INITIALISE FILENAME ARRAY FOR NUMARIS3/INTERNAL IMAGES
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
;           FILE_INPUT: DECLARES AND PRESETS fn_arr_in
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
FUNCTION init_fn_array_numa3_internal, nrslices, nrtimepoints, $
                                       onr,firststudy,firstimage,slicesperstudy

COMMON FILE_INPUT, FN_ARR_IN, PATH_IN, REF_FN_IN

fn_arr_in=STRARR(nrtimepoints,nrslices)

FOR time=0L,nrtimepoints-1L DO BEGIN
    FOR slice=0L,nrslices-1L DO BEGIN
        status = create_fn_numa3_internal(fn_out,                        $
                                 onr,firststudy,firstimage,slicesperstudy,  $
                                 timepoint,slice)
        fn_arr_in[time,nrslices]=fn_out
    ENDFOR
ENDFOR

RETURN,  (0)
END



;+
; NAME: init_fn_array_hst
;
;
;
; PURPOSE: INITIALISE FILENAME ARRAY (hst FORMAT)
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
;           FILE_INPUT: DECLARES AND PRESETS fn_arr_in
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
FUNCTION init_fn_array_hst, nrslices, nrtimepoints, $
                            handnr

COMMON FILE_INPUT, FN_ARR_IN, PATH_IN, REF_FN_IN

fn_arr_in=STRARR(nrtimepoints,nrslices)

FOR time=0L,nrtimepoints-1L DO BEGIN
    FOR slice=0L,nrslices-1L DO BEGIN
        status = create_fn_hst(fn_out,              $
                               handnr,time,slice)
        fn_arr_in[time,nrslices]=fn_out
    ENDFOR
ENDFOR

RETURN,  (0)
END

;+
; NAME: init_fn_array_nih
;
;
;
; PURPOSE: INITIALISE FILENAME ARRAY (nih PACS FORMAT)
;          JUST A QUICK HACK, NEEDS TO BE REPLACED BY SOMETHING SOUND
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
;           FILE_INPUT: DECLARES AND PRESETS fn_arr_in
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

FUNCTION init_fn_array_nih, nrslices, nrtimepoints

COMMON FILE_INPUT, FN_ARR_IN, PATH_IN, REF_FN_IN

FN_ARR_IN=STRARR(nrtimepoints,nrslices)

FOR time=0L,nrtimepoints-1L DO BEGIN
    FOR slice=0L,nrslices-1L DO BEGIN
        fn_arr_in[time,slice]=string(time*nrslices+slice, $
                                     FORMAT='(%"i%d")')
    ENDFOR
ENDFOR

RETURN,  (0)
END



;+
; NAME: write_batch_data.pro
;
;
;
; PURPOSE: write the contents of the filename array to a ASCII file
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
FUNCTION write_batch_data, fn,          $
                           inf_duration,$ ;;inf. duration in seconds
                           rows,        $ ;; output matrix size
                           cols,        $ ;; output matrix size,
                           slices,      $ ;; number of slices
                           timepoints , $ ;; number of timepoints
                           type,        $ ;; file data type
                           fnarr,       $ ;; filename[slices,timepts]
                           reltime,     $ ;; acq. time[slices,timepts]
                           pIDString,   $ ;; ��GJ, 03, 2004, patient ID string
                           tError,      $ ;; Error message
                           Slo=Slo,     $ ;; Signal level offset[slices,timepts]
                           Noi=Noi        ;; Background noise level[slices,timepts]

;; clear Error status
tError=''
ON_IOERROR,ioerror

; check fnarr for 2 dimensions (slice, time)
IF (size(fnarr))[0] NE 2 THEN BEGIN
 ;;; ISSUE ERROR MESSAGE AND EXIT WITH ERROR
    tError='write_batch.pro: fnarr is not of dimension 2!!'
    GOTO,error
ENDIF

; ERROR HANDLING IS MISSING
OPENW,unit,fn,/get_lun

; write a header
printf,unit,"##batchfile v1.0"
printf,unit,inf_duration,        FORMAT='(%"%f")'
printf,unit,slices,timepoints,   FORMAT='(%"%d %d")'
printf,unit,rows,cols,           FORMAT='(%"%d %d")'
printf,unit,type,                FORMAT='(%"%s")'

; write the data
FOR slice=0L,slices-1 DO BEGIN
    FOR time=0L,timepoints-1 DO BEGIN
        printf,unit,slice,time,reltime[slice,time],fnarr[slice,time],$
          FORMAT='(%"%d %d %f %s")'
    ENDFOR
ENDFOR

printf,unit,pIDString,           FORMAT='(%"%s")'

IF KEYWORD_SET(Slo) AND KEYWORD_SET(Noi) THEN BEGIN
    ; write the data
    FOR slice=0L,slices-1 DO BEGIN
        FOR time=0L,timepoints-1 DO BEGIN
            printf,unit,slice,time,reltime[slice,time],Slo[slice,time],Noi[slice,time],$
            FORMAT='(%"%d %d %f %f %f")'
        ENDFOR
    ENDFOR
ENDIF

FREE_LUN,unit
status=0
RETURN,  (0)
IF 0 EQ 1 THEN BEGIN
ioerror: tError=!ERROR_STATE.MSG
error: RETURN,(1)
ENDIF
END


;+
; NAME: read_batch_data.pro
;
;
;
; PURPOSE: read a data file produced by write_batch_data
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
FUNCTION read_batch_data,  fn,          $
                           infusion_dur,$ ;; inf. duration in seconds
                           rows,        $ ;; output matrix size
                           cols,        $ ;; output matrix size,
                           slices,      $ ;; number of slices
                           timepoints,  $ ;; number of timepoints
                           type,        $ ;; file data type
                           fnarr,       $ ;; filename[slices,timepts]
                           reltime,     $ ;; acq. time[slices,timepts]
                           pIDString,   $ ;; ��GJ, 03, 2004, patient ID string
                           tError,      $
                           Slo,         $ ;; ��GJ, July 08, 2004, Signal level offset
                           Noi            ;; ��GJ, July 08, 2004, Background noise level
;; clear Error status
tError=''
ON_IOERROR,ioerror

OPENR,unit,fn,/get_lun

; predeclare all strings
input_path=''
type=''
output_path=''
output_ext=''
version=''
;; and because IDL initialises all variables as floats
;; also predefine all things which need to stay integer
rows=0L
cols=0L
slices=0L
timepoints=0L

;; read a header
readf,unit, version,FORMAT='(%"%s")'
readf,unit, infusion_dur,FORMAT='(%"%f")'
readf,unit, slices, timepoints,FORMAT='(%"%d %d")'
readf,unit, rows, cols,FORMAT='(%"%d %d")'
readf,unit, type,FORMAT='(%"%s")'

;; generate arrays
reltime=dblarr(slices,timepoints)
fnarr=strarr(slices,timepoints)

slc=0L
tm=0L
r=0.D0
fi=''
;; read the slice data
FOR s=0L,slices-1 DO BEGIN
    FOR t=0L,timepoints-1 DO BEGIN
        readf,unit,slc,tm,r,fi, FORMAT='(%"%d %d %f %s")'
        reltime[slc,tm]=r
        ;; KB: despite of the format string above, the first character
        ;; of fi is the " " used as a field separator (think this is
        ;; an IDL bug, too) - so call strtrim to get rid of leading spaces
        fnarr[slc,tm]=strtrim(fi,1)
        ;;        print,' read filename <',fnarr[slc,tm],'>'
    ENDFOR
ENDFOR

pIDString=''
IF EOF(unit) EQ 0 THEN BEGIN
    readf,unit, pIDString
ENDIF

IF EOF(unit) EQ 0 THEN BEGIN
    ;; generate arrays
    Slo=dblarr(slices,timepoints)
    Noi=dblarr(slices,timepoints)
    slc=0L
    tm=0L
    slo_temp=0.D0
    noi_temp=0.D0

    FOR s=0L,slices-1 DO BEGIN
        FOR t=0L,timepoints-1 DO BEGIN
            readf,unit,slc,tm,r,slo_temp,noi_temp,FORMAT='(%"%d %d %f %f %f")'
            Slo[slc,tm]=slo_temp
            Noi[slc,tm]=noi_temp
            ;reltime[slc,tm]=r
            ;; KB: despite of the format string above, the first character
            ;; of fi is the " " used as a field separator (think this is
            ;; an IDL bug, too) - so call strtrim to get rid of leading spaces
            ;fnarr[slc,tm]=strtrim(fi,1)
            ;;        print,' read filename <',fnarr[slc,tm],'>'
        ENDFOR
    ENDFOR
ENDIF

FREE_LUN,unit
status=0

RETURN,  (0)
ioerror: tError=!ERROR_STATE.MSG
error: RETURN,(1)

END


;; write a dat file
;; major modifications in writing process
;; compared to Ulf
pro write_dat_file,filename,mean,times,comment,tau
Ndata=N_ELEMENTS(mean)

openw,unit,filename,/get_lun

;; row 1-3: comment, number of data points and duration of infusion
printf,unit,comment,FORMAT='(%"%s")'
printf,unit,Ndata,FORMAT='(%"%d")'
printf,unit,tau,FORMAT='(%"%f")'

;; write the data points
FOR i=0,Ndata-1 DO BEGIN
a=mean[i]
b=times[i]
printf,unit,a,b,FORMAT='(%"%f %f")'
ENDFOR

;; close file and END execution
FREE_LUN, unit
END

;; read a dat file
;; major modifications in reading process (read ASCII, not buffer and
;; process later)
pro read_dat_file,filename,mean,times,comment,Ndata,tau
text='' &Ndata = 1L &tau = 0.D0

openr,unit,filename,/get_lun

;; row 1-3: comment, number of data points and duration of infusion
readf,unit,text,FORMAT='(%"%s")'
comment=text
readf,unit,Ndata,FORMAT='(%"%d")'
readf,unit,tau,FORMAT='(%"%f")'

;; create mean and times
mean=dindgen(Ndata)
times=mean

;; read in the data points
FOR i=0L,Ndata-1 DO BEGIN
readf,unit,a,b,FORMAT='(%"%f %f")'
mean[i]=a
times[i]=b
ENDFOR

;; close file and END execution
FREE_LUN, unit
END




;; write a dat file
;; major modifications in writing process
;; compared to Ulf
pro write_language_file,filename,icon_dir

  openw,unit,filename,/get_lun

  ;; row 1-3: comment, number of data points and duration of infusion
  printf,unit,icon_dir,FORMAT='(%"%s")'
  ;; close file and END execution
  FREE_LUN, unit
END

;; read a dat file
;; major modifications in reading process (read ASCII, not buffer and
;; process later)
pro read_language_file,filename,icon_dir
  text=''
  
  openr,unit,filename,/get_lun

  ;; row 1-3: comment, number of data points and duration of infusion
  readf,unit,text,FORMAT='(%"%s")'
  icon_dir=text
  ;; close file and END execution
  FREE_LUN, unit
END
















;+
; NAME:get_dicom3_hdr,fn,header
;
;
;
; PURPOSE: return a dicom3 header structure from an existing dicom file
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
FUNCTION get_dicom3_hdr,fn,header
header=0;; initialize header
obj = OBJ_NEW( 'IDLffDICOMex', fn)
status = obj->Read(fn)
;;ERROR_HANDLING             ;
if (status eq 0) THEN BEGIN
    print,FORMAT='(%"get_dicom3_hdr: ERROR on reading file <%s>: IDLffDICOM::read() return status is %d")',fn,status
    return,1
ENDIF

;; get the group
refs = obj->GetReference()

count=0
FOR i = 0L, N_ELEMENTS(refs)-1 DO BEGIN
    IF PTR_VALID(obj->GetValue(REFERENCE=refs[i])) THEN BEGIN
    count=count+1
    ENDIF
ENDFOR
print,FORMAT='(%"get_dicom3_hdr: Number of non-<NULL> ptrs: %d (of %d tags)")',count,N_ELEMENTS(refs)

;; create a dicom header structure of appropriate size
header=replicate({struc_dicom},count)

;; IDL 5.x is not able of assigning a VR for (0x0028,0x0120)
;; the interesting thing is, that there are tags like that,
;; whose content can be dumped using ->DumpElements, but
;; on using GetValue for this tag a NULL pointer will be returned

;; This means that in case of safety only tags with valid pointers
;; will be handled - so count the number of valid tags before
;;                   proceeding



;;
    ;; not even UK - it just returns '  ' (2 blanks) - Thanks, RSI
    ;; which is either US or SS, so deal with this tag by hand
;;


;; copy tags into the structure
ccount=0
FOR i = 0L, N_ELEMENTS(refs)-1 DO BEGIN
    IF PTR_VALID(obj->GetValue(REFERENCE=refs[i])) THEN BEGIN
    ;; (KB) I'm no fan of IDL pointers - so derefence it in two steps
    ;; Seems so awkward in comparison to C
    ;; create a copy of the pointer - to be freed later (must NOT use
    ;;                                /NO_COPY!!)
    ;; The construction below fixes that for strings the pointer is
    ;; just copied and deallocated on destroying the IDLffDICOM object
    ;; ::this is the syntax needed for accessing the object
    ;; - Oh how I love plain C!!
    ;;    c=*((header[i].pVal)[0])
    ;; help, header[i].pVal
        val = *((obj->GetValue(REFERENCE=refs[i]))[0])
        header[ccount].pVal =PTR_NEW(val)

        header[ccount].group = obj->GetGroup(REFERENCE=refs[i])
        header[ccount].element = obj->GetElement(REFERENCE=refs[i])
        header[ccount].vr = obj->GetVR(REFERENCE=refs[i])
        header[ccount].length = obj->GetLength(REFERENCE=refs[i])
        ccount = ccount +1
    ENDIF
ENDFOR
OBJ_DESTROY, obj

RETURN,0
END













; BEAUTIFIED VERSION OF ULFs extract_dicom
FUNCTION extract_dcm2,header,type,first,second,TAG=tag,START=start

; FUNCTION to extract ACR-NEMA header information.
; If the TAG-keyword is set, the parameters first and second are the
; group and element numbers of the desired variable.
; ELSE (DEFAULT) first is the starting byte and second the number of bytes.

; Autor Ulf Hoffmann            Datum: 25.10.1996

on_error,2                      ; return to caller
within_header = 1

IF KEYWORD_SET(TAG) THEN BEGIN  ; first=group, second=element
    tag={group:0,$
         element:0,$
         length:0L}
    i = 0
    s=size(header)
    sh=s(n_elements(s)-1)

    WHILE (tag.group NE first) OR (tag.element NE second) $
      AND within_header DO BEGIN

        tag.group   = FIX(header,i+0)
        tag.element = FIX(header,i+2)
        tag.length  = LONG(header,i+4)

        i = i + 8 + tag.length
        if i ge sh THEN within_header = 0
    ENDWHILE
    start = i - tag.length
    ende  = i - 1
ENDIF $
ELSE   BEGIN
    start = first
    ende  = second
ENDELSE

if not within_header THEN BEGIN
    message, 'TAG not found' ,/traceback,/informational
    return,-1
ENDIF

CASE type of
    2 :return,FIX   (header,start)
    3 :return,LONG  (header,start)
    4 :return,FLOAT (header,start)
    5 :return,DOUBLE(header,start)
    7 :return,STRING(header(start:ende))
ENDCASE

end                             ; of extract_dcm2

;****************************************************************************

FUNCTION get_struc_image_info_dcm3, obj
struc_image_info ={struc_image_info}; defined in struc_image_info__define.pro

; Version: 09/06/02

;Image Matrix
IF (obj->QueryValue('0028,0010')) EQ 2 THEN struc_image_info.rows = obj->GetValue('0028,0010')

IF (obj->QueryValue('0028,0011')) EQ 2 THEN struc_image_info.columns = obj->GetValue('0028,0011')

;; pat info
IF (obj->QueryValue('0010,0010')) EQ 2 THEN struc_image_info.pat_name = obj->GetValue('0010,0010')

;; pat id
IF (obj->QueryValue('0010,0020')) EQ 2 THEN struc_image_info.pat_id = obj->GetValue('0010,0020')

; Slice location
IF (obj->QueryValue('0020,1041')) EQ 2 THEN struc_image_info.sl_pos = obj->GetValue('0020,1041') ELSE struc_image_info.sl_pos =0.

; Acquisition time
IF (obj->QueryValue('0008,0032')) EQ 2 THEN BEGIN
    time = obj->GetValue('0008,0032')
ENDIF ELSE BEGIN
     IF (obj->QueryValue('0008,0033')) EQ 2 THEN time = obj->GetValue('0008,0033')
ENDELSE

if time(0) ne -1 then begin
    time = STRTRIM(time,1) ;; fixed: REMOVE trailing blanks
    h = fix(strmid(time,0,2))
    m = fix(strmid(time,2,2))
    s = fix(strmid(time,4,2))
    f = double(strmid(time,6))
endif else begin
    h=0 & m=0 & s= 0 & f=0
    ; Now attempt to get time from the trigger time

    ; read in the trigger time
    ;value = obj->GetValue('0018'x,'1060'x, /NO_COPY)
    ;t=0.0
    ;if(ptr_valid(value[0])) then t=*value[0]
    ;if(t gt 0.0) then struc_image_info.time = t/1000.0 $ ; time in second
    ;else struc_image_info.time=0             ;

endelse

; Acquisition date

year = 1800
mon  = 1
day  = 1

; Acquisition time
IF (obj->QueryValue('0008,0022')) EQ 2 THEN BEGIN
    date = obj->GetValue('0008,0022')
ENDIF ELSE BEGIN
     IF (obj->QueryValue('0008,0023')) EQ 2 THEN date = obj->GetValue('0008,0023')
ENDELSE

if date ne -1 then begin
    date = STRTRIM(date,1) ;; fixed: REMOVE trailing blanks
    year = fix(strmid(date,0,4))
    mon  = fix(strmid(date,4,2))
    day  = fix(strmid(date,6,2))
end

; Conversion to JULDAY
if (f NE 0) THEN PRINT,'FRACTION=',f
struc_image_info.time=JULDAY(FLOAT(mon),FLOAT(day),FLOAT(year),FLOAT(h),FLOAT(m),s+f)

;��GJ, 04, 2004, add the trigger time
IF (obj->QueryValue('0018,1060')) EQ 2 THEN struc_image_info.triggertime = obj->GetValue('0018,1060') ELSE struc_image_info.triggertime =0.

;��GJ, 04, 2004, add the reference start time
IF (obj->QueryValue('0054,1300')) EQ 2 THEN struc_image_info.frameStarttime = DOUBLE(obj->GetValue('0054,1300'))/1000. ELSE struc_image_info.frameStarttime =0.

; Sequence name
IF (obj->QueryValue('0018,0020')) EQ 2 THEN struc_image_info.seq_name = (obj->GetValue('0018,0020'))[0] ELSE struc_image_info.seq_name =''

;Pixel spacing:
IF (obj->QueryValue('0028,0030')) EQ 2 THEN struc_image_info.pixel_spacing = (obj->GetValue('0028,0030'))[0] ELSE struc_image_info.pixel_spacing =0.

; Slice Thickness
IF (obj->QueryValue('0018,0050')) EQ 2 THEN struc_image_info.sl_thick = obj->GetValue('0018,0050') ELSE struc_image_info.sl_thick =0.

;��GJ, add the echo time of each image 03, 2004
; Echo time
IF (obj->QueryValue('0018,0081')) EQ 2 THEN struc_image_info.echotime = obj->GetValue('0018,0081') ELSE struc_image_info.echotime =0.

;��GJ, add the flip angle of each image Aug 24, 2004
; flip angle
IF (obj->QueryValue('0018,1314')) EQ 2 THEN struc_image_info.flangle = obj->GetValue('0018,1314') ELSE struc_image_info.flangle =0.

debug=0
IF debug NE 0 THEN BEGIN
print,'   We are in the DICOM Format'
print,'   Acquisition time:  ',h,' h ',m,' min ',s,' sec ',f
print,'   Acquisition date:  ',day,'. ',mon,'. ',year
print,'   Sequence File Name:  ',struc_image_info.seq_name
print,'   Matrix            :',struc_image_info.rows,' x ',struc_image_info.columns
print,'   Pixel size:      ',struc_image_info.pixel_spacing
print,'   Slice thickness: ',struc_image_info.sl_thick
ENDIF

RETURN, struc_image_info
END




FUNCTION get_messung_struct_numa3, header
struc_image_info ={struc_image_info}; defined in struc_image_info__define.pro

;Image Matrix - G28.Pre.Rows
s=fix(header,4994)
byteorder,s,/SWAP_IF_LITTLE_ENDIAN,/SSWAP
struc_image_info.rows = s

;; G28.Pre.Columns
s=fix(header,4996)
byteorder,s,/SWAP_IF_LITTLE_ENDIAN,/SSWAP
struc_image_info.columns = s

;; pat info - G10.Pat.PatientName
struc_image_info.pat_name = string(header(768:768+25))
;; G10.Pat.PatientId
struc_image_info.pat_id = string(header(795:795+11))

;; Slice location - was header.G51.Txt.SlicePosition (string(header(5805:5805+8)))
;; switched to G21.Rel1.CM.ImageDistance
d= double(header,3816)
byteorder,d,/SWAP_IF_LITTLE_ENDIAN,/L64SWAP
struc_image_info.sl_pos=d

;; Acquisition time (h,m,s,f as consecutive longwords)
t = long(header,52,4)
byteorder,t, /SWAP_IF_LITTLE_ENDIAN,/LSWAP
;; cast the type to double to produce a double Juldate
t=double(t)

;; Acquisition date (y,m,d as  consecutive longwords)
y = long(header,12,3)
byteorder,y, /SWAP_IF_LITTLE_ENDIAN,/LSWAP
;; cast the type to double to produce a double Juldate
y=double(y)

; Conversion to JULDAY
struc_image_info.time=JULDAY(y[1],y[2],y[0],t[0],t[1],t[2]+t[3]/1.D3)

;; Sequence name
struc_image_info.seq_name = string(header[3009:3009+64])

;Pixel size:
;; I hope using only one direction will be OK
x_pix = double(header,5000)
byteorder,x_pix,/SWAP_IF_LITTLE_ENDIAN,/L64SWAP
struc_image_info.pixel_spacing = x_pix

; Slice Thickness
d = DOUBLE( header,1544)
byteorder,d,/SWAP_IF_LITTLE_ENDIAN,/L64SWAP
struc_image_info.sl_thick= d


debug=0
IF debug NE 0 THEN BEGIN
print,'   We are in the SIEMENS VISION Format'
n=N_TAGS(struc_image_info)
FOR i=0,n-1 DO BEGIN
    PRINT,FORMAT='(%"%20s: <%s>")',  (TAG_NAMES(struc_image_info))[i],struc_image_info.(i)
ENDFOR

;;    print,'   Sequence File Name:  ',struc_image_info.seq_name
;;    print,'   Matrix            :',struc_image_info.rows,' x ',struc_image_info.columns
;;    print,'   Pixel size:      ',struc_image_info.pixel_spacing
;;    print,'   Slice thickness: ',struc_image_info.sl_thick
ENDIF


RETURN, struc_image_info
END


;; this may also serve as a
;; template for file handling
FUNCTION write_ampk21_cluster_text,fn,a_buf,k21_buf,tError
tError='' ;; clear error status
ON_IOERROR,error

count = size(k21_buf, /N_ELEMENTS)

IF (count EQ 0) THEN RETURN,0

OPENW, lun, fn, /GET_LUN

PRINTF, lun, FORMAT = '("Amp   kep")'
FOR i=0L, count-1 DO BEGIN
    printf, lun, a_buf[i], k21_buf[i], FORMAT='((f7.3,x), f7.3)'
ENDFOR

FREE_LUN, lun

RETURN,0

error:
tError=!ERROR_STATE.MSG
RETURN,1
END

;; this may also serve as a
;; template for file handling
FUNCTION write_ampkep_kel_cluster_text,fn,mask_index, a_buf,kep_buf, kel_buf, tError
tError='' ;; clear error status
ON_IOERROR,error

count = size(kep_buf, /N_ELEMENTS)

IF (count EQ 0) THEN RETURN,0

OPENW, lun, fn, /GET_LUN

PRINTF, lun, FORMAT = '("Mask_index", %"\t","Amp", %"\t", "kep", %"\t", "kel")'
FOR i=0L, count-1 DO BEGIN
    printf, lun, mask_index[i], a_buf[i], kep_buf[i], kel_buf[i], FORMAT='((f7.0,x),%"\t",(f7.3,x),%"\t", (f7.3,x), %"\t", (f7.3,x))'
ENDFOR

FREE_LUN, lun

RETURN,0

error:
tError=!ERROR_STATE.MSG
RETURN,1
END



;; this may also serve as a
;; template for file handling
FUNCTION read_ampkep_kel_cluster_text,fn,mask_index, a_buf,kep_buf, kel_buf, tError
tError='' ;; clear error status
ON_IOERROR,error


OPENR, lun, fn, /GET_LUN

title=''
mask_index=0.D
a_buf=0.D
kep_buf=0.D
kel_buf=0.D

READF, lun, title, FORMAT='(%"%s")'
WHILE ~ EOF(lun) DO BEGIN
    mask_index_temp=0.D
    temp_a=0.D
    temp_kep=0.D
    temp_kel=0.D
    READF, lun, mask_index_temp, temp_a, temp_kep, temp_kel,FORMAT='(4F0)'
    mask_index = [[mask_index], mask_index_temp]
    a_buf=[[a_buf], temp_a]
    kep_buf=[[kep_buf], temp_kep]
    kel_buf=[[kel_buf], temp_kel]
ENDWHILE

FREE_LUN, lun

mask_index=mask_index[1:N_ELEMENTS(mask_index)-1]
a_buf=a_buf[1:N_ELEMENTS(a_buf)-1]
kep_buf=kep_buf[1:N_ELEMENTS(kep_buf)-1]
kel_buf=kel_buf[1:N_ELEMENTS(kel_buf)-1]

RETURN,0
count = size(kep_buf, /N_ELEMENTS)

IF (count EQ 0) THEN RETURN,0

error:
tError=!ERROR_STATE.MSG
RETURN,1
END






FUNCTION get_messung_struct_dcm2, header
messung ={pat_name:'',pat_id:'',sl_pos:0.0,time:0.0}

messung.pat_name = extract_dcm2(header,7,'10'XL,'10'XL,/TAG)
messung.pat_id   = extract_dcm2(header,7,'10'XL,'20'XL,/TAG)
                                ; LOCATION (RETARDED IN DICOM3.0!!)
messung.sl_pos   = float(extract_dcm2(header,7,'20'XL,'50'XL,/TAG))

                                ; acquisition time
time =  extract_dcm2(header,7,'8'XL,'32'XL,/TAG)
if time(0) ne -1 THEN BEGIN
    h = fix(strmid(time,0,2))
    m = fix(strmid(time,3,2))
    s = fix(strmid(time,6,2))
        f = fix(strmid(time,9,3))
     ENDIF ELSE BEGIN
        h=0 & m=0 & s= 0 & f=0
     endELSE

     ;Acquisition date
     date = extract_dcm2(header,7,'8'XL,'22'XL,/TAG)
     if date ne -1 THEN BEGIN
       year = fix(strmid(date,0,4))
       mon  = fix(strmid(date,5,2))
       day  = fix(strmid(date,8,2))
     end ELSE BEGIN
       year = 0
       mon  = 0
       day  = 0
     end

     ; Sequence name
     seq_file_name = extract_dcm2(header,7,'18'XL,'20'XL,/TAG)

     ;Pixel size:
     pix = extract_dcm2(header,7,'28'XL,'30'XL,/TAG)

     ; Slice Thickness
     slthk = extract_dcm2(header,7,'18'XL,'50'XL,/TAG)

RETURN, messung
END

;+
; NAME: read_mri_data.pro
;
; PURPOSE: general routine for reading of an MRI image
;          returns pixel data and struc_image_info header-structure
;
; CATEGORY:
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

FUNCTION read_mri_data, fn, pixeldata, header_struc, type, tError,$
                        NUMA2=numa2, NUMA3=numa3, DICOM3=dicom3, HEADER=header
; read in a mri data file and try to make a correct guess of its
; data format

; check for DICOM3.0: IF IDL can read it, it has to be dicom
; unfortunately, idl opens dicom2.0 too, and issues an error message
; instead of properly returning

; initialise output data type

ON_IOERROR,ioerror
type='undefined'

;; as IDL has some problems with ACR-NEMA1.0 ("DICOM2.0") ( ACR/NEMA Standards Publication No. 300-1988 )
;; first try to detect an ACR-NEMA1.0 file before using IDL for
;; opening the file

;; unfortunately this forces us to open the file twice - thanX IDL
OPENR,unit,fn,/get_lun
fsize=(fstat(unit)).size

 ;; read the first 0x2b bytes of the file and close it
IF fsize EQ 0 THEN BEGIN
    tError=STRING(FORMAT='(%"File <%s> has size  %d bytes")',fn,fsize)
    FREE_LUN, UNIT
    goto,rerror
ENDIF

IF fsize GT '2b'x THEN fsize='2b'x
data = bytarr(fsize)
READU, unit, data
FREE_LUN, UNIT

;; look for tag (0008,0010) -  "RecognitionCode" which should be "ACR-NEMA 1.0"
;; for the good old SIEMENS MAGNETOM SP
;; I don't want to create potential problems with non Siemens
;; ACR-NEMA 1.0 appl. (I don't want to process them)
;; SO LETS RELY ON SIEMENS' FIXED HEADER AND ONE OF ITS
;; ELEMENTS WHICH I DEEM TO BE CHARACTERISTIC - bytes 0x20-0x2b
;; WHICH CONTAIN THE Regcognition code

IF fsize GT '2b'x THEN BEGIN
    IF STRCMP(data['20'x:'2b'x],'ACR-NEMA 1.0',12) THEN BEGIN
        ;; heureka, it's ACR-NEMA

;; the tag struc was not defined, so define it now
;; someone will remedy the ACR-NEMA file format problems anyways
        tag={tag_struc,$
             group:0,$
             element:0,$
             length:0L}

        row_tag = {tag_struc,'28'XL,'10'XL,2L}
        col_tag = {tag_struc,'28'XL,'11'XL,2L}
        data_tag = {tag_struc,'7fe0'XL,'10'XL,0L}
        type='DCM2'
                                ; separate header and pixel data
        tag.group = FIX(data,0)

        ;;KB: What is the proper initialisation for tag.element??
        i=0
        WHILE (tag.group NE data_tag.group) OR $
          (tag.element NE data_tag.element) DO BEGIN

            tag.group   = FIX(data,i+0)
            tag.element = FIX(data,i+2)
            tag.length  = LONG(data,i+4) ;,index=[0,1,2]
            if      comp_struc(tag,row_tag) THEN rows    = fix(data,i+8) $
            ELSE if comp_struc(tag,col_tag) THEN columns = fix(data,i+8)
            i = i + 8 + tag.length
        ENDWHILE
        start = i - tag.length
        imageformat = extract_dcm2(data,7,'8'XL,'10'XL,/TAG)
;;        print, 'File <', fn, '> is a Numaris 2 file!! and',fsize,' bytes long'
                                ;header kopieren
        header = bytarr(start)
        header = data(0:start-1)
        pixeldata = intarr(7,13) ;something weird for testing as we don't extract pixeldata yet
        goto, done ;; don't spend time on beautifying useless IF statements ...
    ENDIF
ENDIF

;; try to open it as a DICOM3.0  file - rely on IDLs crappy code

;;��jia.11@osu.edu, IDLffDICOM is replaced by IDLffDICOMex
obj = OBJ_NEW( 'IDLffDICOMex', fn)
;var = obj->Read(fn)             ;

;query (jun) needs to added here.
var=1

IF (var EQ 1) THEN BEGIN
;;    print, 'File <', fn, '> is a dicom file!!' ;
    type='DCM3'

    ; get the pixeldata
    obj->GetProperty, BITS_ALLOCATED = vBA, ROWS=rows,COLUMNS=cols, SAMPLES_PER_PIXEL=samples
    order=1
    vPixels = obj->GetPixelData(ORDER=order, COUNT=cnt)
    ; select monochrome image
    ;select the last sample data (need to be added here).
    IF samples gt 1 THEN BEGIN
        pixeldata=vPixels
;        a=dialog_message('multiple sample per pixel, add codes here, line 2129, jia.11@osu.edu. file_io.pro')
    ENDIF ELSE BEGIN
        pixeldata=vPixels
    ENDELSE

;;    help, pixeldata
    header_struc = get_struc_image_info_dcm3(obj);
    OBJ_DESTROY, obj

ENDIF ELSE BEGIN ; No DCM3.0 - check for SIEMENS NUMARIS XX header
    OBJ_DESTROY, obj   ; get rid of the possibly existing dicom object
;; ERROR HANDLING IS MISSING
    header=bytarr(6144)         ; length of a NUMARIS3/internal header

    OPENR,unit,fn,/get_lun
    fsize=(fstat(unit)).size

    IF (fsize LE 6144) THEN BEGIN
        ;; ERROR HANDLING - THIS CANNOT BE A VALID NUMARIS3 FILE
        tError=STRING(FORMAT='(%"File <%s> is not a NUMARIS3 file. Size: %d bytes")',fn,fsize)
        FREE_LUN, UNIT
        goto,rerror
    ENDIF

                                ; read the file and close it
    data = bytarr(fsize)
    READU, unit, data
    FREE_LUN, UNIT

   ;; now check for a NUMARIS3 header
    IF STRCMP(data[6*16:6*16+6],'SIEMENS',7) THEN BEGIN
        type='NUMA3'
        print, 'File <', fn, '> is a Numaris 3 file!! and',fsize,' bytes long'

        ;; image size: 128x128 or 256x256??
        if fsize eq 530432 then begin
            columns = 512       ; 256
            rows    = 512       ; 256
        end else $
          if fsize eq 137216 then begin
            columns = 256
            rows    = 256
            end else $
              if fsize eq 38912 then begin
                columns = 128
                rows    = 128
            end else begin
                tError=STRING(FORMAT='(%"Wrong size for VISION NUMARIS3 image: File has only %d bytes")',fsize)
                goto, rerror
        end

        pixeldata=FIX(data,6144L,columns,rows)
        byteorder,pixeldata,/SWAP_IF_LITTLE_ENDIAN
        header_struc = get_messung_struct_numa3(data[0:6143]) ;

    ENDIF ELSE BEGIN        ; is it a pseudo-NUMARIS3 parameter image?
        IF ARRAY_EQUAL(data[0:512*12-1],0*bytarr(512*12)) THEN BEGIN

            PRINT,'(KB) read parameter image, <NUMA3> format'
                      ; 128^2, 256^2 or 512^2??
            if fsize eq 530432 then begin
                columns = 512   ; 256
                rows    = 512   ; 256
            end else $
              if fsize eq 137216 then begin
                columns = 256
                rows    = 256
            end else $
              if fsize eq 38912 then begin
                columns = 128
                rows    = 128
            end else begin
                tError=STRING(FORMAT='(%"Wrong size for VISION NUMARIS3 image: File has only %d bytes")',fsize)
                goto, rerror
            end
        ENDIF
    ENDELSE
ENDELSE
;; insert file name in header structure
header_struc.file_name=fn
;done: ;;print,header_struc
done:RETURN,0

ioerror:tError=!ERROR_STATE.MSG
rerror:RETURN,1
END




;+
; NAME: assign_alternate_morphology
;
;
;
; PURPOSE: Find read in slices with same position as TimeSeries images
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:p=assign_alternate_morphology(patient,filelist)
;
;
;
; INPUTS:
; pMorphInfo: Pointer to structure containing TimeSeries info
; filelist  : List of files to be processed
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
; Pointer array, dim.=#of slices, found slices with valid pointers,
;                other with NULL pointers
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
; Passing the whole patient structure is overkill and a bit slow .....
;
;
; RESTRICTIONS: IF more than one images match with the slice, the
; last processed one will be taken, which may not be the best match.

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
; Dec. 6, 20002 creation baudendistel.4@osu.edu
;-
FUNCTION assign_alternate_morphology,patient,fnames,p,tError
;; create a pointer array with #of slices
ON_IOERROR,ioerror
IF fnames[0] NE '' THEN BEGIN ;; 7d-safe
    p=PTRARR(patient.n_slices)

    ;; diff. in slice thickness
    dt=.5 * (*(*(patient[0].pMorphInfo)).pSAimage_info)[0,0].sl_thick
    ;; process the files
    FOR i=0L,N_ELEMENTS(fnames)-1 DO BEGIN
        IF read_mri_data(fnames[i],pixel,h,t,tError) THEN goto,ioerror
        ;; check, whether the locations and matrix sizes match
        ;; this is only a crude check and should be more elaborated
        ;; by comparison of ALL relevant features
        FOR sl=0L,patient.n_slices-1 DO BEGIN
            pos=(*(*(patient[cur_pat].pMorphInfo)).pSAimage_info)[sl,0].sl_pos
            IF (patient.ima_rows EQ h.rows) AND $
              (patient.ima_rows EQ h.columns) THEN BEGIN

                ;; now check for the slice position
                IF (ABS(pos-h.sl_pos) LT dt) THEN BEGIN
                    IF PTR_VALID(p[sl]) THEN PTR_FREE,p[sl]
                    ;; just insert it
                    print,FORMAT='(%"assign_alternate_morphology: FOUND MATCH <%s> for slice (%d). %f matched by %f")',$
                      fnames[i],sl,pos,h.sl_pos
                    p[sl]=PTR_NEW(pixel,/NO_COPY)
                ENDIF
            ENDIF
        ENDFOR
    ENDFOR
ENDIF

RETURN,0

ioerror:tError=!ERROR_STATE.MSG
FOR i=0,N_ELEMENTS(p)-1 DO BEGIN
    IF PTR_VALID(p[i]) THEN PTR_FREE,p[i]
ENDFOR
RETURN,1
END



;
;

;; Is s1 < s2 (TRUE=1, FALSE=0)
;; compare up to minimum length of both names

FUNCTION compare_fn,s1,l1,s2,l2
b=0
len =min([l1,l2])

; compare up to minimum length of both names
b= strmid(s1,0,len) LE strmid(s2,0,len)

;fix problem with  "i2", "i100" (-> "i2" < "i1")
IF (b EQ 0) AND (l1 LT l2) THEN RETURN, 1

; both equal up to len characters
IF (b EQ 1) AND (l1 GT l2) THEN RETURN, 0


RETURN, B
END









;+
; NAME: sillysort
;
;
;
; PURPOSE: sort filenames in numeric order (i1,i2, i10, i11, ..i99, i100, ...)
;
; It may look silly to split the extension from the files and to sort
; the result and then resort the original array using an index, but
; this even works if the file extensions differ in small/capital letters
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE: sillysort,inputarray,outputarray
;                   (works with inputarray=outputarray)
;
;
;
; INPUTS: inputarray: Array of filenames
;
;
;
; OPTIONAL INPUTS: -
;
;
;
; KEYWORD PARAMETERS: -
;
;
;
; OUTPUTS: outputarray: Array of sorted filenames
;
;
;
; OPTIONAL OUTPUTS: -
;
;
;
; COMMON BLOCKS: -
;
;
;
; SIDE EFFECTS: -
;
;
;
; RESTRICTIONS: -
;
;
;
; PROCEDURE: -
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;  Sep, 10, 2002 Creation  (K. Baudendistel)
;-
PRO sillysort,fnarr_in, fnarr_out

n=N_ELEMENTS(fnarr_in)
; check

IF (n LT 1) THEN BEGIN
    print,'sillysort: only ',n,' filenames given -> ERROR'
    print, 5/0
ENDIF


; create an empty array for the file names (without extensions)
f=strarr(n)

; array for the indices (used for resorting original input array)
ind = replicate(n,n); force idl to create an array from the datatype of n (could be 16, 32, or 64bit)

; array for the string lengths
l = lonarr(n)

; remove path and extension

FOR i=0UL, n-1 DO BEGIN
 fdecomp,fnarr_in[i],dr,dir,file,q,v
 f[i]=file
 ind[i]=i
 l[i]=STRLEN(file)
ENDFOR

; General sorting idea: Look for a smaller element in the array.
; If there is one, swap these two

FOR i=0UL, n-1 DO BEGIN         ; source position
    REPEAT BEGIN
        restart=0
        FOR j=i, n-1 DO BEGIN   ; filename to be compared against
            IF i EQ j THEN CONTINUE ; would compare against itself


            b=compare_fn(f[i],l[i],f[j],l[j])

            ; SWAP??
            IF (b EQ 0) AND (i LT j) THEN BEGIN
                ; swap files indexed by i and j
                h=f[i]    & f[i]  =f[j] & f[j]=h
                h=l[i]    & l[i]  =l[j] & l[j]=h
                h=ind[i]  & ind[i]=ind[j]    & ind[j]=h
                restart=1
                BREAK
            ENDIF

        ENDFOR
    ENDREP UNTIL (restart EQ 0)
ENDFOR

; now we have resorted the indices - reorder original input array
fnarr_out=fnarr_in; declare output array
tmp=fnarr_in;       copy in case of in-place argument

fnarr_out=tmp[ind]; resort it
END

FUNCTION browse_directory, ref_fn, fn_files, tError
; The initial file was picked and is passed by ref_fn.
; All filenames in the same directory and extension are stored in fn_files

 fdecomp,ref_fn, disk, dir, name, qual, version

; print,disk, dir, name, qual, version,FORMAT='(%"<%s+%s+%s+%s+%s>\n")'

 fcomp,f,disk,dir,'*',qual,version

 IF STRMID(f, 0, 1, /reverse_offset) EQ '.' THEN f=STRMID(f, 0, STRLEN(f)-1)

 print,'browse_directory: browsing for files matching pattern <',f,'>'

 fn_files=file_search(f)

 IF fn_files[0] NE '' THEN BEGIN
     tError=''
     return, (0)
 ENDIF ELSE BEGIN
     return, (1)
 ENDELSE

END

;;��GJ, August 05, 2003, fix the bug that the sort is not based on slice location. (Houston data were obtained as image # 1,2,3,4, with same slice locations)
;;��GJ, March 20, 2003, fix the bug that different series are not allowed to have different number of images
FUNCTION sort_DCM_by_ser_and_ima_number, fnarr_in, fnarr_out, numSlices, tError
;;Sort DICOM3.0 images according their series number and image number
;;Call method:
;;file_io.pro, FUNCTION read_all_images_to_matrix, ref_fn, sl_slab, Imat, Imat_inf,tError
;;replace sillysort to SORT_BY_IMAGE_NUMBER
;; we hope that this gives us the image filenames in the order:
;; Timept1 (slice1, slice2, ..., sliceN)
;; ...
;; TimeptM (slice1, slice2, ..., sliceN)
;;
;; if not, then we'll need something else
;;
;; Images are first sorted by the SERIES number
;; And within the SERIES by their IMAGE number

n=N_ELEMENTS(fnarr_in)
; check number of input files

IF (n LT 1) THEN BEGIN
    tError=STRING(FORMAT='(%"ERROR: sort_by_image_number: can not sort %d filenames")',number)
    GOTO,rerror
ENDIF

;PRINT, 'Before (BEFORE SORTING):'
;PRINT, fnarr_in

;; create structure for sorting

file_info_array=REPLICATE({file_index:0UL,       $
                           series_number:0UL,    $
                           image_number:0UL,     $
                           slice_location:0.D0,  $
                           image_position:DBLARR(3)*0.D0,  $
                           slice_thickness:0.D0}, $
                           n)
ON_IOERROR, error

;; for catching errors with DICOM input try CATCH here
;; (I hate IDL for the fact that the CATCH code has to
;; come prior to possible errors instead of having implemented
;; a try-catch loop) - so this has to jump to the error handler
;; at the END of the subroutine

;;BEGIN OF CATCH ERROR HANDLER
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    ;;Handle the error
    tError=!ERR_STRING
    CATCH, /CANCEL;; needed??
    GOTO, rerror
ENDIF
;; END OF CATCH ERROR HANDLER

FOR i=0UL, n-1 DO BEGIN
    file_info_array[i].file_index=i

    obj = OBJ_NEW( 'IDLffDICOMex', fnarr_in[i], /NO_PIXEL_DATA)

    ; Series Number
    IF (obj->QueryValue('0020,0011')) EQ 2 THEN file_info_array[i].series_number = obj->GetValue('0020,0011') ELSE file_info_array[i].series_number =0

    ; Image Number
    IF (obj->QueryValue('0020,0013')) EQ 2 THEN file_info_array[i].image_number = obj->GetValue('0020,0013') ELSE file_info_array[i].image_number =0
    ;for multiple image number
    ;value = obj->GetValue('0020'x,'0013'x, /NO_COPY)
    ;nele=N_ELEMENTS(value)
    ;tempvalue=DBLARR(nele)
    ;FOR j=0, nele-1 DO tempvalue[j]=*value[j]
    ;file_info_array[i].image_number= MAX(tempvalue)

    ; Slice position
    IF (obj->QueryValue('0020,1041')) EQ 2 THEN file_info_array[i].slice_location = obj->GetValue('0020,1041') ELSE file_info_array[i].slice_location =0

    ; Image position
    IF (obj->QueryValue('0020,0032')) EQ 2 THEN file_info_array[i].image_position = obj->GetValue('0020,0032')

    ; Slice thickness
    IF (obj->QueryValue('0018,0050')) EQ 2 THEN file_info_array[i].slice_thickness = obj->GetValue('0018,0050') ELSE file_info_array[i].slice_thickness =0


    OBJ_DESTROY, obj

ENDFOR


;; check for different slice thicknesses
max_th=MAX(file_info_array.slice_thickness, min=min_th)

IF max_th NE min_th THEN BEGIN
    tError=STRING(FORMAT='(%"Error:Data contains differing slice thicknesses from %f mm to %f mm!")',min_th,max_th)
        goto,rerror
ENDIF

min_series=MIN(file_info_array.series_number)
file_info_array.image_number=(file_info_array.series_number-min_series)*2000.+file_info_array.image_number
file_info_array=file_info_array(SORT(file_info_array.image_Number))
tempIndex=0L
k=0L
FOR i=0L, n-2 DO BEGIN
    IF file_info_array[i].image_Number NE file_info_array[i+1].image_Number THEN BEGIN
        tempIndex=[tempIndex, i+1]
        k=k+1
    ENDIF
ENDFOR
newfilearray=file_info_array[tempIndex]
file_info_array=newfilearray

;;Basic sorting based on image_number
numSlices=N_ELEMENTS(UNIQ(file_info_array.slice_location, SORT(file_info_array.slice_location)))
numPositions=DBLARR(3)
FOR j=0, 2 DO numPositions[j]=N_ELEMENTS(UNIQ(file_info_array.image_position[j], SORT(file_info_array.image_position[j])))
IF numSlices NE MAX(numPositions, maxind) THEN BEGIN
    file_info_array.slice_location=file_info_array.image_position[maxind]
    numSlices=N_ELEMENTS(UNIQ(file_info_array.slice_location, SORT(file_info_array.slice_location)))
ENDIF
numTimepoints=N_ELEMENTS(file_info_array.slice_location)/numSlices
PRINT, '# of slices: ', numSlices, ', # of timepoints: ', numTimepoints
IF numSlices*numTimepoints NE N_ELEMENTS(file_info_array) THEN BEGIN
    tError='Incorrect number of images'
    GOTO, rerror
ENDIF

UniqSL=(file_info_array.slice_location)(UNIQ(file_info_array.slice_location, SORT(file_info_array.slice_location)))
IF uniqSL[0] NE file_info_array[0].slice_location THEN UniqSL=REVERSE(UniqSL)
nSL=N_ELEMENTS(UniqSL)
IF nSL NE n THEN BEGIN
    SLinfoarray=file_info_array
    FOR i=0, nSL-1 DO BEGIN
        indexk=0
        FOR j=0, n-1 DO BEGIN
            IF ABS(UniqSL[i]-file_info_array[j].slice_location) LT 0.00001 THEN BEGIN
              SLinfoarray[indexk*nSL+i]=file_info_array[j]
              indexk=indexk+1
            ENDIF
        ENDFOR
    ENDFOR
ENDIF
file_info_array=SLinfoarray

; now we have resorted the indices - reorder original input array

fnarr_out=fnarr_in; declare output array
tmp=fnarr_in;       copy in case of in-place argument

fnarr_out=tmp[file_info_array.file_index]; resort it

;; exit with status OK
RETURN,0
;; return with error
error: tError=!ERROR_STATE.MSG
rerror: PRINT, tError & RETURN, 1
END

;;##################################TESTROUTINE+##########################################
;;PRO test_sort_DCM_by_ser_and_ima_number
;;filename='D:\images\Image_test_DCE_MRI\browser\kentucky\1eew6.dcm'
;;filename='D:\images\Image_test_DCE_MRI\browser\prostate\1.2.840.113619.2.5.1762527380.1977.1040568491.580.dcm'
;;filename = dialog_pickfile(filter=['*.dcm'],/MUST_EXIST,TITLE='select first image file to be processed')
;;dummy=browse_directory(filename,fn_files,status)
;;IF (N_ELEMENTS(fn_files) NE 0) THEN BEGIN
;;  status=sort_DCM_by_ser_and_ima_number(fn_files, fn_files, sl_slab,tError)
;;  PRINT, 'status: ', status
;;ENDIF
;;END
;;##################################TESTROUTINE+##########################################


;;;; unfortunately the sorting is too slow using the sillysort for
;;;; SIEMENS images with internal naming. Because of that the following "crutch" is added

pro vision_internal_format_crutch, ref_fn, fn_files
;; had big fun with regexps: Try
;; g='-4--3-4-5-65-6---43-4-3-23-3-3-3-232-' AND
;; l=strsplit(g,"\-[0-9]+\-",/REGEX,LENGTH=m) & for i=0,N_ELEMENTS(l)-1 DO BEGIN & print,strmid(g,l[i],m[i]) & ENDFOR
;; this will aslo return the "-" of the "---" group, DESPITE THE FACT THAT
;; the regexp specifies at least one number between two "-" as
;; matching pattern -AAAAAAAARRRRRRRRRGGGGGGGGHHHHHHHHHHHH!!!!
;; so this is the best approx. I could manage (just 1h of
;; wasting time, ThanX Creaso)

;; look at the first file again:
;; if it does not contain a hyphen it may be a ULF-renamed file
;; if the return of strsplit is the complete string, then most
;; likely the search pattern is not present - but there is no
;; docu about retvals in this case

;; the ref-filename should contain the number of the first image to be processed
fdecomp, ref_fn, disk, dir,name, qual, version
l=strsplit(name,'-',LENGTH=m)
IF m[0] NE strlen(name) THEN BEGIN ;; THEN IT IS A VALID INTERNAL VISION FORMAT
    i=N_ELEMENTS(l)-1
    l1=l[i] & m1=m[i]
    first_nr=long(strmid(name,l1,m1))

    ;; get the last file to be processed
    lfn=''
    IF !version.os EQ "linux" THEN BEGIN
        lfn='/home/klaus/newdata/b072_bopta/2308-5-244.ima'
        lfn='/home/klaus/GUANG_DOUBTFUL/b072_bopta/2308-9-500.ima'
;;        lfn='/home/klaus/newdata/b072_bopta/2308-5-181.ima'
    ENDIF ELSE BEGIN
        REPEAT BEGIN
            cat_fpaths,p,disk,dir
            lfn=DIALOG_PICKFILE(/MUST_EXIST,PATH=p,TITLE='select last image to be processed',FILTER='*.ima')
        ENDREP UNTIL (lfn NE '')
    ENDELSE


    ;; and get the image number of the last input file
    fdecomp, lfn, disk, dir,name, qual, version
    l=strsplit(name,'-',LENGTH=m)
    IF m[0] NE strlen(name) THEN BEGIN
        ;; get the last image
        fdecomp, lfn, disk, dir,name, qual, version
        l=strsplit(lfn,'-',LENGTH=m)
        i=N_ELEMENTS(l)-1
        l1=l[i] & m1=m[i]
        last_nr=long(strmid(lfn,l1,m1))

        n=lonarr(N_ELEMENTS(fn_files))
        ;; take the input array and get the last "-" position
        FOR f=0L,N_ELEMENTS(fn_files)-1 DO BEGIN
            l=strsplit(fn_files[f],'-',LENGTH=m)
            i=N_ELEMENTS(l)-1
            l1=l[i] & m1=m[i]
            n[f]=long(strmid(fn_files[f],l1,m1))
        ENDFOR


        ;; sort by image number (last number string)
        s=SORT(n)

        ;; check against ambiguity here

        IF N_ELEMENTS(s) NE N_ELEMENTS(UNIQ(SORT(n))) THEN BEGIN
            ;; ERROR HANDLER
            STOP
        ENDIF

        ;; find the positions of the first and last image to be processed
        first=where(n[s] EQ first_nr)

        IF (first[0] EQ -1) THEN BEGIN
            ;; ERROR HANDLER
            STOP
        ENDIF

        last=where(n[s] EQ last_nr)
        IF (last[0] EQ -1) THEN BEGIN
            ;; ERROR HANDLER
            STOP
        ENDIF

        ;; resort the file names
        fn_files=fn_files[s]

        ;; and get rid of the things outside the range to be processed
        tmp=fn_files[first:last]
        fn_files=tmp
    ENDIF ELSE BEGIN
        ;; error handler missing
        STOP
    ENDELSE
ENDIF ELSE BEGIN
    ;; ULF special format (Urgh)
    ;; the following is complete crap, but at leasts it works ....
    ;; I try to do everything without opening the files, which would
    ;; be a lot easier but more time consuming
    ;; this really is sorted in timepoints/slices, so rearrange it
    ;; into slices/timepoints

    ;; get the rightmost four digits of the filename

    ;; all filenames have the same length - so getting indices from
    ;;                                      ref. is OK
    fdecomp, ref_fn, disk, dir,name, qual, version

    ;; now take the last four characters from the right
    pattern=strmid(name,3,4,/REVERSE_OFFSET)

    ;; find out where the "." of the qualifier starts
    q = STRPOS(ref_fn,qual,/REVERSE_OFFSET,/REVERSE_SEARCH)

    ;; mysteriously q is counted from the beginning
    ;; what is good for me

    ;; and get the slice/timepoint string  location in the absolute path
    sp=q-5
    tp=q-3

    n=N_ELEMENTS(fn_files)
    ;; and now sort the images so consecutive slice positions will be
    ;; in consecutive files

    h=long(STRMID(fn_files[0:n-1],tp,2))*1000+long(STRMID(fn_files[0:n-1],sp,2))

    s=SORT(h)

    ;; check against ambiguity here

    IF N_ELEMENTS(s) NE N_ELEMENTS(UNIQ(SORT(h))) THEN BEGIN
        ;; ERROR HANDLER
        STOP
    ENDIF

    ;; resort the file names
    fn_files=fn_files[s]
ENDELSE
END

;; look for images consecutively until same pos occurs
FUNCTION fancycount_nr_slices,fn_files,sl_slab,tError
;; Read first slice ( ...again, I know)
num_images=N_ELEMENTS(fn_files)
status=read_mri_data(fn_files[0],pixeldata,struc_image_info, type, tError)
;; Error handling
IF status THEN BEGIN
    GOTO, error_ret
ENDIF
IF sl_slab eq 0 THEN BEGIN
    ;;read until we encounter the same slice position again
    ;; because consecutive images represent successive slice positions
    ;; at least we assume this until we have our BROWSER for
    ;; getting DICOM input
    test=double(struc_image_info.sl_pos)
    i=0
    REPEAT BEGIN
        i=i+1
        status=read_mri_data(fn_files[i],pixeldata,struc_image_info_curr, type, tError)
        IF status THEN BEGIN
            GOTO, error_ret
        ENDIF
        sl_pos = double(struc_image_info_curr.sl_pos)
    ENDREP UNTIL ((ABS(sl_pos - test) LE .1D-4) OR (i EQ N_ELEMENTS(fn_files)-1))

    sl_slab=i

    ;; added because of problems with Hee's testdata
    IF (i EQ  N_ELEMENTS(fn_files)-1) AND (struc_image_info_curr.sl_pos NE test) THEN BEGIN
        ;; in this case the same slice position did never occur again
        sl_slab = N_ELEMENTS(fn_files)
    ENDIF
ENDIF
RETURN, (0)
error_ret: return, (1)
END


;; get all 2-komp. parametric model  image input
;; look for all images matching <digit><digit>parext based on
;; ref. file with <digit><digit>_amp or <digit><digit>_AMP

FUNCTION get_2comp_param_image_filenames, ref_fn, found, tError,PAREXT=parext
IF NOT (KEYWORD_SET(parext)) THEN parext='_amp'
searchext='_amp'
;; get a list of _amp.dcm, _kep.dcm, _t.dcm files

;; separate the input filename from the name_part <digit><digit>_amp.dcm
fdecomp, ref_fn, disk, dir, name, qual, version

;; first: check for the file extension
;;   -not yet, because we may have not a file extension
;; be case-insensitive
q = STRPOS(STRUPCASE(name),STRUPCASE(searchext),/REVERSE_OFFSET,/REVERSE_SEARCH)

;; mysteriously q is counted from the beginning
;; what is good for me

IF q EQ (STRLEN(name)-STRLEN(searchext)) THEN BEGIN
    n=STRMID(name,0,q-2)+'??'+parext
    fcomp, fn, disk, dir, n, qual, version

    ;; and glob it ...
    found=findfile(fn)
    IF found[0] NE '' THEN BEGIN
        RETURN,(0)
    ENDIF
ENDIF
tError='No match for pattern <digit><digit>'+parext
found=''
error: RETURN, (1)
END


;;�� YXY: Check if the reference images and the DCE study images in Bruno Morgan model
;;          were arranged in the same order.
FUNCTION check_reference_images, struc_imat_inf_ref, struc_imat_inf, tError

    row_ref = struc_imat_inf_ref.n_rows
    row_study = struc_imat_inf.n_rows
    col_ref = struc_imat_inf_ref.n_cols
    col_study = struc_imat_inf.n_cols
    IF (row_ref-row_study GT 0.01) OR (col_ref-col_study GT 0.01) THEN BEGIN
        tError = 'The reference images should have the same size with the study images!'
        GOTO, error_ret
    ENDIF

    slice_ref = struc_imat_inf_ref.n_slices
    slice_study = struc_imat_inf.n_slices
    IF (slice_ref - slice_study GT 0.01) THEN BEGIN
        tError = 'The reference images should have the same number of slices with the study images!'
        GOTO, error_ret
    ENDIF

    loc_ref = DBLARR(slice_ref)
    loc_study = DBLARR(slice_ref)
    FOR sl = 0, slice_ref-1 DO BEGIN
        loc_ref[sl] = (*(struc_imat_inf_ref.pSAimage_info))[sl, 0].sl_pos
        loc_study[sl] = (*(struc_imat_inf.pSAimage_info))[sl, 0].sl_pos
    ENDFOR
    dloc = ABS(loc_ref - loc_study)
    IF MAX(dloc) GT 1.0e-3 THEN BEGIN
        tError = 'The reference images should have the same slice locations with the study images!'
        print, 'Reference Image slice locations:', loc_ref
        print, 'Study Image slice locations:', loc_study
        GOTO, error_ret
    ENDIF

    RETURN, (0)

    error_ret: RETURN, (1)

END


;; and sort the images in the correct way
;; read the images and compare them to the patient's list
FUNCTION read_and_sort_param_images, filelist, pMorphInfo, pParamVol,tError,PAREXT=parext

IF filelist[0] EQ '' THEN BEGIN
    PRINT,'read_and_sort_param_images: filelist is empty!!'
    tError='read_and_sort_param_images: filelist is empty!!'
    GOTO,error_ret
ENDIF

;; NEW BY KB - SWITCH FOR ANALYZE FORMAT
CASE ((*pMorphInfo).type) OF
    "ANALYZE":BEGIN
        ;; open a multi-dimensional analyze file

        STOP
    END
    ELSE: BEGIN
        FOR i=0L,N_ELEMENTS(filelist)-1 DO BEGIN
            ;; read in the file header - as it is dicom use idls
            ;;                           function for that
            obj = OBJ_NEW( 'IDLffDICOM' )
            status = obj->Read(filelist[i])
            IF (status NE 1) THEN BEGIN ;; ERROR HANDLING
                print,FORMAT='(%"Error reading dicom file: <%s>")',filelist[i]
                stop
            ENDIF
            ;; now look to which image it belongs
            ;; Assume we are sure about patient, just compare
            ;; slice positions (still to be defined what has to
            ;; be compared!!! - could be relocated to function compare(ima_struc,dcm3header))

            FOR sl=0L,(*pMorphInfo).n_slices-1 DO BEGIN ;; slices

;;        print,FORMAT='(%"processing file <%s>")',filelist[sl]

                loc = (*((*pMorphInfo).pSAimage_info))[sl,0].sl_pos
                ;; slice location
                value = obj->GetValue('0020'x,'1041'x, /NO_COPY)
                slloc = *value[0] ;
;;        print,FORMAT='(%"AMP: read sl=%d Compare ref %f with %f")',sl,slloc,loc
                IF (ABS(slloc-loc) LT 1.0e-3) THEN BEGIN ;; (rounding error <1.0e-3 seems reasonable)

                    ;; rows
                    value = obj->GetValue('0028'x,'0010'x, /NO_COPY)
                    rows = *value[0] ;

                    ;; columns
                    value = obj->GetValue('0028'x,'0011'x, /NO_COPY)
                    columns = *value[0] ;

                    ;; pixeldata
                    value = obj->GetValue('7fe0'x,'0010'x, /NO_COPY)
                    pixeldata = *value[0]

                    ;; OK I know that dimensions should be checked ....
                    CASE (parext) OF
                        '_amp':BEGIN
                            print,FORMAT='(%"Found slice #%d: file <%s>")',sl,filelist[sl]
                            (*pParamVol)[0:columns-1,0:rows-1,sl,0]=FLOAT(pixeldata)/400.
                        END
                        '_kep':BEGIN
                            print,FORMAT='(%"Found slice #%d: file <%s>")',sl,filelist[sl]
                            (*pParamVol)[0:columns-1,0:rows-1,sl,1]=FLOAT(pixeldata)/20.
                        END
                        '_kel':BEGIN
                            print,FORMAT='(%"Found slice #%d: file <%s>")',sl,filelist[sl]
                            (*pParamVol)[0:columns-1,0:rows-1,sl,2]=(FLOAT(pixeldata-2048))/2000.
                        END
                        '_t':BEGIN
                            print,FORMAT='(%"Found slice #%d: file <%s>")',sl,filelist[sl]
                            (*pParamVol)[0:columns-1,0:rows-1,sl,3]=FLOAT(pixeldata)/500.
                        END
                        '_auc':BEGIN
                            print,FORMAT='(%"Found slice #%d: file <%s>")',sl,filelist[sl]
                            (*pParamVol)[0:columns-1,0:rows-1,sl,4]=FLOAT(pixeldata)/1000.
                        END
                        '_slope':BEGIN
                            print,FORMAT='(%"Found slice #%d: file <%s>")',sl,filelist[sl]
                            (*pParamVol)[0:columns-1,0:rows-1,sl,5]=FLOAT(pixeldata)/100.
                        END
                        ELSE: BEGIN
                            print,FORMAT='(%"Unhandled switch <%s>")',parext
                            STOP
                        ENDELSE
                    ENDCASE
                    break ; we have found a match so assume we are done
                ENDIF
            ENDFOR
            OBJ_DESTROY, obj
        ENDFOR ;;  FOR i=0L,N_ELEMENTS(found) DO BEGIN
    ENDELSE
ENDCASE
RETURN, (0)

error: tError=!ERROR_STATE.MSG
error_ret:RETURN, (1)
END

;; This could be enlarged to read other images like analyze format, too....
FUNCTION read_all_images_from_analyze_file, ref_fn, sl_slab, Imat, Imat_inf,tError
ON_IOERROR,error
;; now the ref_fn must be the name of the .img file
;; check using query_analyze
IF NOT query_analyze(ref_fn,info) THEN GOTO, ERROR
;; SUCCESSFULL QUERY -> read the file

data=read_analyze(ref_fn,  info)
IF N_ELEMENTS(data) EQ 1 THEN BEGIN
    tError='FILE NOT FOUND'
    GOTO, error_ret
ENDIF
;; create a struc_image_info - structure
;; and preset it from analyze header -
;; this will be the universal image header for all the
;; images in this analyze file
image_info={struc_image_info}
image_info.file_name=ref_fn
image_info.pat_name =info.name
image_info.pat_id   =info.scan_id
image_info.sl_pos=0.D0
image_info.time=0.D0 ;; julian date of first image
image_info.seq_name='unknown'
image_info.pixel_spacing=0.D0
image_info.sl_thick =info.thickness
image_info.reltime  =0.D0 ;; in seconds
image_info.rows     =info.dimensions[0]
image_info.columns  =info.dimensions[1]

sl_slab=N_ELEMENTS(data[0,0,*])
time_point=1
IF (size(data))[0] EQ 4 THEN BEGIN
    time_point=N_ELEMENTS(data[0,0,0,*])
ENDIF

;; and try to preset the image structure with the timepoints
;; extracted here
; create the array holding all images
IMAT=make_array(image_info.columns,image_info.rows,sl_slab,$
                time_point,/UINT)

;; create a 3D-structure-array with the image info
;; as kind of a meta-header
imat_inf={struc_imat_inf}

;; now create the structure info for the individual images ....
;; replicate as many times as needed
ASimage_info=replicate({struc_image_info}, $
                       sl_slab*time_point)

;; and change the indices to the desired 3D array
ASimage_info=reform(ASimage_info, sl_slab, time_point, /OVERWRITE)


;; Read images into the matrix

FOR index_time_point=0L, time_point-1 DO BEGIN
    FOR index_sl_slab=0L, sl_slab-1 DO BEGIN
        ;; copy to return arguments
        imat[0:image_info.columns-1,0:image_info.rows-1,$
             index_sl_slab, index_time_point]=$
          data[0:image_info.columns-1,0:image_info.rows-1,$
               index_sl_slab, index_time_point]

        ;; copy our info structur
        ASimage_info[index_sl_slab,index_time_point]=image_info
        ;; as there is no time info in the analyze header, initialise
        ;; reltime[] as 0 seconds
        ASimage_info[index_sl_slab,index_time_point].reltime= 0.0

        ;; and as slice position we take the index of the slice in the volume
        ASimage_info[index_sl_slab,index_time_point].sl_pos = index_sl_slab
    ENDFOR
ENDFOR
imat_inf.type        = 'ANALYZE'
imat_inf.fn          = ref_fn
imat_inf.n_rows      = info.dimensions[0]
imat_inf.n_cols      = info.dimensions[1]
imat_inf.n_slices    = sl_slab
imat_inf.n_timepoints= time_point
;; and insert image info in imat_inf structure
imat_inf.pSAimage_info = PTR_NEW(ASimage_info,/NO_COPY)
;; and carry on the analyze header in our meta header
imat_inf.pAnalyzehdr=PTR_NEW(info,/NO_COPY)
RETURN,0;; NO ERROR
error: tError=!ERROR_STATE.MSG
error_ret:RETURN, (1)
END

;; This could be enlarged to read other images like analyze format, too....
FUNCTION read_all_images_from_dyn_file, ref_fn, sl_slab, Imat, Imat_inf,infusion_dur,tError
ON_IOERROR,error
;; .dyn-file : read the .dyn file
s=read_batch_data(  ref_fn,      $
                    infusion_dur,$ ;; inf. duration in seconds
                    rows,        $ ;; output matrix size
                    cols,        $ ;; output matrix size,
                    sl_slab,     $ ;; number of slices
                    time_point,  $ ;; number of timepoints
                    type,        $ ;; file extension
                    fnarr,       $ ;; filename[slices,timepts]
                    reltime,     $ ;; acq. time[slices,timepts]
                    pIDString,   $ ;; ��GJ,patient ID string
                    tError,      $
                    Slo,         $ ;; ��GJ, Signal level offset
                    Noi)           ;; ��GJ, Background noise level

IF S THEN BEGIN
    print, tError
    goto, error
ENDIF
;; bring the filenames in a useful order
;; the expected order is TIME1(slice1,..,sliceN) TIME2(slice1,...,sliceN)
fn_files=strarr(sl_slab*time_point)

;; and a simple loop, because otherwise I have to figure out the
;; mapping of the arrays
FOR t=0L,time_point-1 DO BEGIN
    FOR sl=0L,sl_slab-1 DO BEGIN
        fn_files[t*sl_slab+sl]=fnarr[sl,t]
    ENDFOR
ENDFOR
;;
CASE type OF
    "ANALYZE":BEGIN
        ;; all timepoints have to refer to the same analyze file
        FOR i=1L,N_ELEMENTS(fn_files)-1 DO BEGIN
            IF fn_files[i] NE fn_files[0] THEN BEGIN
                tError='different analyze input files specified'
                GOTO, error_ret
            ENDIF
        ENDFOR
        ;; read the analyze file
        s=read_all_images_from_analyze_file(fn_files[0], sl_slab, Imat, Imat_inf,tError)
        IF s THEN BEGIN
            GOTO, error_ret
        ENDIF

        ;; and modify reltime according to the specified times
        (*imat_inf.pSAimage_info)[0:sl_slab-1,0:time_point-1].reltime=reltime[0:sl_slab-1,0:time_point-1]
    END
    ELSE:BEGIN
        ;; open the first image and get its file type
        status=read_mri_data(fn_files[0],pixeldata,struc_image_info, type, tError)
        IF status THEN BEGIN
            GOTO, error_ret
        ENDIF
        ;; OK, I know this is not good -but let us take the first real
        ;; image name as reference for the filenames in the dynamik-program
        ref_fn=fn_files[0]

        ;; and try to preset the image structure with the timepoints
        ;; extracted here
        ;; create the array holding all images
        IMAT=make_array(struc_image_info.columns,struc_image_info.rows,sl_slab,$
                        time_point,/UINT)


        ;; create a kind of a meta-header
        imat_inf={struc_imat_inf}

        ;; now create the structure info for the individual images ....
        ;; replicate as many times as needed
        ASimage_info=replicate({struc_image_info}, $
                               sl_slab*time_point)

        ;; and change the indices to the desired 3D array
        ASimage_info=reform(ASimage_info, sl_slab, time_point, /OVERWRITE)

        ;; Read images into the matrix
        FOR index_time_point=0L, time_point-1 DO BEGIN
            FOR index_sl_slab=0L, sl_slab-1 DO BEGIN
                index=(index_time_point*(sl_slab))+index_sl_slab
                status=read_mri_data(fn_files[index],pixeldata,struc_image_info, type)
                IF status THEN BEGIN
                    GOTO, error_ret
                ENDIF
                ;; copy to return arguments
                imat[0:struc_image_info.columns-1,0:struc_image_info.rows-1,$
                     index_sl_slab, index_time_point]=pixeldata

                ASimage_info[index_sl_slab,index_time_point]=struc_image_info

                ;; reltime has to be replaced by the already determined time
                ;; offsets provided by the .dyn file (in seconds)
                ASimage_info[index_sl_slab,index_time_point].reltime=$
                  reltime[index_sl_slab,index_time_point]
                IF N_ELEMENTS(Slo) NE 0 AND N_ELEMENTS(Noi) NE 0 THEN BEGIN
                    ASimage_info[index_sl_slab,index_time_point].Slo=$
                      Slo[index_sl_slab,index_time_point]
                    ASimage_info[index_sl_slab,index_time_point].Noi=$
                      Noi[index_sl_slab,index_time_point]
                ENDIF
            ENDFOR
        ENDFOR

        imat_inf.type          = type
        imat_inf.n_rows        = ASimage_info[0,0].rows
        imat_inf.n_cols        = ASimage_info[0,0].columns
        imat_inf.n_slices      = sl_slab
        imat_inf.n_timepoints  = time_point
        ;; and insert image info in imat_inf structure
        imat_inf.pSAimage_info = PTR_NEW(ASimage_info,/NO_COPY)
        ;;��GJ
        imat_inf.pIDString=pIDString
    ENDELSE
ENDCASE

RETURN,0;; NO ERROR
error: tError=!ERROR_STATE.MSG
error_ret: RETURN, (1)
END


;; get all morphology image input
FUNCTION read_all_images_to_matrix, ref_fn, sl_slab, Imat, Imat_inf,tError,DYN_FILE=dyn_file
ON_IOERROR,error
;; ref_fn - first image file
;; sl_slab - number of images per slab, if eq 0 then determine


;; INSERT FEB2003: IF INPUT IS A .dyn-file, THEN BEHAVE DIFFERENTLY
;; reset the dyn_file. I know, this is not the best way to do it, but
;; a convenient shortcut for the time being (KB)
dyn_file=[0,0];; is initialized as NOT a structure anymore, is changed to a structure, if a .dyn is read

fdecomp,ref_fn, disk, dir, name, qual, version
;; some modification for implementing analyze and dyn format.
;; For single images, never trust an extension - as we have
;; some SIEMENS .ima -files which are in reality DICOM3.0
;; but have the wrong qualifier

CASE STRUPCASE(qual) OF
    'DYN': BEGIN ;; our own proprietary format
        s=read_all_images_from_dyn_file(ref_fn, sl_slab, Imat, Imat_inf,infusion_dur,tError)
        IF s THEN BEGIN
            GOTO, error_ret
        ENDIF

        dyn_file={DYN_MODIFIES_PARAS,dyn_file:1.0,infusion_dur:infusion_dur}
    END
    'IMG': BEGIN ;; analyze data cubes
        s=read_all_images_from_analyze_file(ref_fn, sl_slab, Imat, Imat_inf,tError)
        IF s THEN BEGIN
            GOTO, error_ret
        ENDIF

    END
    ELSE: BEGIN ;; ordinary images
        ;;
        ;; it based on rep. of the same slice position in input files

        ;; get all filename candidates
        dummy=browse_directory(ref_fn,fn_files,tError)

        ;; error check
        IF fn_files[0] EQ '' THEN BEGIN
            tError='NO MATCH FOR PATTERN'+ref_fn
            goto, error_ret
        ENDIF

        ;; open the first image and get its file type
        status=read_mri_data(fn_files[0],pixeldata,struc_image_info, type, tError)
        IF status THEN BEGIN
            goto, error_ret
        ENDIF


        ;; unfortunately reading all images in all supported formats
        ;; in the correct
        ;; order is endlessly complicated
        ;; so deal with it using a CASE statement


        ;;           I M P O R T A N T
        ;; at the end of the case statement fn_arr has to contain the
        ;; filenames in the order timept1: slice1, slice2, ... timept2:
        ;; slice1, slice2 ......

        ;; num_images has to contain the total number of images
        ;; sl_slab has to contain the number of slices

        CASE (type) OF
            'NUMA3': BEGIN
                ;; we need to call the vision internal format crutch
                ;; to sort out the images

                vision_internal_format_crutch, ref_fn, fn_files
                num_images=N_ELEMENTS(fn_files)
                s=fancycount_nr_slices(fn_files,sl_slab,tError)
                IF status THEN BEGIN
                    GOTO, error_ret
                ENDIF
            END

;;        'NUMA2': BEGIN
;;            STILL TO BE IMPLEMENTED. SO FOR NOW USE ERROR HANDLING BELOW
;;        END

            'DCM3': BEGIN
                ;;           sillysort,fn_files,fn_files
                sl_slab=0
                IF sort_DCM_by_ser_and_ima_number(fn_files, fn_files, sl_slab, tError) THEN GOTO,error_ret
                s=fancycount_nr_slices(fn_files,sl_slab,tError)
                ;IF s THEN BEGIN
                 ;   GOTO, error_ret
                ;ENDIF
                num_images=N_ELEMENTS(fn_files)
            END
            ELSE: BEGIN
                tError=STRING(FORMAT='(%"Unhandled data format type <%s>")',type)
                goto, error_ret
            ENDELSE
        ENDCASE

        ;; Determine the number of time points
        time_point=num_images/sl_slab

        print,time_point,sl_slab, $
          FORMAT='(%"read_all_images_to_matrix: we seem to have %d slices and %d time points")'

        IF num_images MOD sl_slab NE 0 THEN BEGIN
            tError='Number of images is inconsistent'
            goto, error_ret
        ENDIF
        ;; create the array holding all images
        IMAT=make_array(struc_image_info.columns,struc_image_info.rows,sl_slab,$
                        time_point,/UINT)

        ;; create a kind of a meta-header
        imat_inf={struc_imat_inf}

        ;; now create the structure info for the individual images ....
        ;; replicate as many times as needed
        ASimage_info=replicate({struc_image_info}, $
                               sl_slab*time_point)

        ;; and change the indices to the desired 3D array
        ASimage_info=reform(ASimage_info, sl_slab, time_point, /OVERWRITE)

        ;; Read images into the matrix
        FOR index_time_point=0L, time_point-1 do begin
            FOR index_sl_slab=0L, sl_slab-1 do begin
                index=(index_time_point*(sl_slab))+index_sl_slab

                ;; and write what we are reading
                print,FORMAT='(%"timepoint: %04d slice:%04d <%s>")',$
                  index_time_point+1,index_sl_slab+1,fn_files[index]

                status=read_mri_data(fn_files[index],pixeldata,struc_image_info, type, tError)
                ;; Error handling
                IF status THEN BEGIN
                    GOTO, error_ret
                ENDIF


                ;; copy to return arguments
                imat[0:struc_image_info.columns-1,0:struc_image_info.rows-1,$
                     index_sl_slab, index_time_point]=pixeldata
;;              help, ((*imat_inf.PSAIMAGE_INFO)[0,2]),/STRUC

                ASimage_info[index_sl_slab,index_time_point]=struc_image_info

                ;; added by KB: I want the time differences in seconds to the
                ;; first image of the same slice

                ASimage_info[index_sl_slab,index_time_point].reltime= 86400D0 * $
                  (struc_image_info.time -ASimage_info[index_sl_slab,0].time)
                ;; I AM SORRY, BUT THE TIMES ARE JUST CRAP FOR GE-images.
            ENDFOR
        ENDFOR

        ;��GJ, April 7, 2004, add the method to get reltime from trigger time and image acquisiton time
        IF ASimage_info[0,1].reltime-ASimage_info[0,0].reltime LT 0.1 THEN BEGIN
            IF ASimage_info[0,1].triggertime-ASimage_info[0,0].triggertime GT 0.1 THEN BEGIN
                FOR index_time_point=1L, time_point-1 do begin
                    IF ASimage_info[0,index_time_point].reltime-ASimage_info[0,(index_time_point-1)].reltime LT 0.1 THEN BEGIN
                        delta_time=ASimage_info[0,index_time_point].triggertime-ASimage_info[0,(index_time_point-1)].triggertime
                        ASimage_info[*,index_time_point].reltime=ASimage_info[*,index_time_point-1].reltime+delta_time
                    ENDIF ELSE BEGIN
                        ;keep the old reltime, need to do nothing
                    ENDELSE
                ENDFOR
            ENDIF
        ENDIF


        ;��GJ, April 22, 2005, add the method to get reltime from reference start time
        IF ASimage_info[0,1].reltime-ASimage_info[0,0].reltime LT 0.1 THEN BEGIN
            IF ASimage_info[0,1].frameStarttime-ASimage_info[0,0].frameStarttime GT 0.1 THEN BEGIN
                FOR index_time_point=1L, time_point-1 do begin
                    IF ASimage_info[0,index_time_point].reltime-ASimage_info[0,(index_time_point-1)].reltime LT 0.1 THEN BEGIN
                        delta_time=ASimage_info[0,index_time_point].frameStarttime-ASimage_info[0,(index_time_point-1)].frameStarttime
                        ASimage_info[*,index_time_point].reltime=ASimage_info[*,index_time_point-1].reltime+delta_time
                    ENDIF ELSE BEGIN
                        ;keep the old reltime, need to do nothing
                    ENDELSE
                ENDFOR
            ENDIF
        ENDIF


;
;  ��GJ, April 12, 2005 to Kishore to read co-registered png images, without sorting them by assuming image number is file name
;        filename=DIALOG_PICKFILE(TITLE='png new data')
;        ;; get all filename candidates
;        IF STRLEN(filename) NE 0 THEN BEGIN
;            dummy=browse_directory(filename,files,tError)
;            FOR i=0, N_ELEMENTS(files)-1 DO BEGIN
;              result=READ_PNG(files[i])
;              imat[0:struc_image_info.columns-1,0:struc_image_info.rows-1,$
;                     0, i]=REVERSE(result, 2)
;            ENDFOR
;        ENDIF
        imat_inf.type          = type
        imat_inf.n_rows        = ASimage_info[0,0].rows
        imat_inf.n_cols        = ASimage_info[0,0].columns
        imat_inf.n_slices      = sl_slab
        imat_inf.n_timepoints  = time_point
        ;; and insert image info in imat_inf structure
        imat_inf.pSAimage_info = PTR_NEW(ASimage_info,/NO_COPY)
    ENDELSE
ENDCASE
RETURN, (0)

error: tError=!ERROR_STATE.MSG
error_ret:RETURN, (1)
END

;; For DICOM output:
;; took Jacks original win dicom writer and modified it a bit
;;
;; and he in turn seems to have borrowed from the internet
;;
;convert a value, val, that is num bytes long, into
;a series of ordered bytes
function getbytes, val, num
    ret=BYTARR(num)
    offset=0
;work in big endian ONLY
    ;val=swap_endian(val)
;if (!version.arch eq 'x86') then begin REPLACED with little endian test
;thanks to David Fanning :)
    little_endian = (BYTE(1, 0, 1))[0]
    if (little_endian) then begin
       byteorder,val,/SWAP_IF_BIG_ENDIAN
    endif else begin
       byteorder,val,/SWAP_IF_LITTLE_ENDIAN
    endelse
    for i=0,(num-1) do begin
       tmpres=BYTE(ISHFT(val, offset) AND 255)
       ret[i]=tmpres
;;   ret[i]=tmpres[0];; makes no sense, but no error message
       offset=offset-8
    endfor

    return, ret
end

;;
;; ROUTINES FOR DEALING WITH DICOM OUTPUT START HERE
;; - TO BE REPLACED AS SOON AS IDL SUPPORTS WRITING DICOM FILES
;;

;; all routines suppose that the header is a list of ascending
;; group/elements, as is DICOM3.0 standard at time of writing

;; remove a header tag
FUNCTION header_tag_rm
return,0
END

;; add a header tag
FUNCTION header_tag_add
return,0
END

;generate any tag
function generate_anytag, group, element, length, data, STR = str

    pad=BYTE(0)

;check to see if string type is set - if it is, change
;padding byte to a space
    IF KEYWORD_SET(STR) then begin
       pad=BYTE(STRING(' '))
    endif

    rs=[getbytes(group,2),getbytes(element,2)]



;correct to even length if necessary
    bs=BYTE(data)
    nl=n_elements(bs)

;; KBs idea: if someone passes -1 as length, then determine length here
        IF length eq -1 THEN BEGIN
            length = nl
            if ((length mod 2) ne 0) then length=length+1 ;
        ENDIF ELSE BEGIN
            if ((length mod 2) ne 0) then length=length+1 ;
        ENDELSE

;   if ((nl mod 2) ne 0) then begin
    if (nl lt length) then begin
       pads = bytarr(length-nl);
       pads[*] = pads[*] + pad
;     bs=[bs,pads]
       if (pad[0] eq 32) then begin
         bs=[pads,bs]
       endif else bs=[bs,pads]
;     nl=nl+1
    end

;size of field
;   rs=[rs,getbytes(nl,2)]
;padding, since nl is a short integer
;   rs=[rs,[0,0]]
    rs = [rs, getbytes(length, 4)]

;string itself
    rs=[rs,bs]

    return, rs
end

;generate explict VR tag
function generate_exp_VRtag, group, element, vr, length, data, STR = str

    pad=BYTE(0)

;check to see if string type is set - if it is, change
;padding byte to a space
    IF KEYWORD_SET(STR) then begin
       pad=BYTE(STRING(' '))
    endif

    rs=[getbytes(group,2),getbytes(element,2)]

    vrb = byte(vr)
    rs=[rs,vrb];

;correct to even length if necessary
    if ((length mod 2) ne 0) then length=length+1;

    bs=BYTE(data)
    nl=n_elements(bs)
    if (nl lt length) then begin
       pads = bytarr(length-nl);
       pads[*] = pads[*] + pad
       bs=[bs,pads]
    end

;size of field
;   rs=[rs,getbytes(nl,2)]
;padding, since nl is a short integer
;   rs=[rs,[0,0]]
;   rs = [rs, getbytes(length, 4)]

    if((vr eq 'OB') or (vr eq 'OW') or (vr eq 'SQ') or (vr eq 'UN') or (vr eq 'UT')) then begin
       tmpi = bytarr(2);
       rs = [rs, tmpi]; reserved field
       ;size of field
       rs=[rs,getbytes(length,4)]
    endif else begin
       ;size of field
       nl = fix(length)
       rs=[rs,getbytes(nl,2)]
    endelse

;string itself
    rs=[rs,bs]

    return, rs
end


;generate empty tag
function generate_empty_tag, group, element

    rs=[getbytes(group,2),getbytes(element,2)]

;size of field
    nl = 0L;
    rs=[rs,getbytes(nl,4)]

    return, rs
end


function generate_VRtag, group, element, VR, length, data

;   print, group, element, VR, length, data
    CASE VR of

       'AE': begin
;Application Entity - normal string tag, truncated to 16 bytes
         dval=STRMID(data,0,16)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'AS': begin
;Age String - should be nnnX, where X={D,W,M,Y} (days, weeks, months, years),
;truncated to 4 bytes
         dval=STRMID(data,0,4)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'AT': begin
;Attribute tag - should be a pair of unsigned integers representing a data
;element tag eg. ['0018'x,'00FF'x]
         dval=[getbytes(UINT(data[0]),2), getbytes(UINT(data[1]),2)]
         rs=generate_anytag(group,element,length,dval)
       end

       'CS': begin
;Code string - 64 byte string
         dval=STRMID(data,0,64)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'DA': begin
;Date string - 8 bytes fixed, formay yyyymmdd, or 10 bytes fixed
;yyyy.mm.dd, which is compatible with versions prior dicom v3.0 -
;so thats what will be used
         dval=STRMID(data,0,10)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'DS': begin
;Decimal string - convert an float into a string, and store
;16 bytes maximum
         dval=STRTRIM(STRING(data),1)
         dval=STRMID(dval,0,64)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'DT': begin
;Date/time string - 26 byte maximum string of format:
;YYYMMDDGGMMSS.FFFFFF
         dval=STRMID(data,0,26)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'FL': begin
;Floating point single - 4 byte fp single val
;storing as LITTLE ENDIAN - needs to be checked!!!!
         dval=FLOAT(data)
;explicily cast to bytes
         dvaltmp=BYTE(dval,0,4)
;fix byteorder (little-endian) if necessary; word length is 16 bits (?is this correct?)
         dval=[getbytes(UINT(dvaltmp[0:1],0,1),2), $
               getbytes(UINT(dvaltmp[2:3],0,1),2)]
         rs=generate_anytag(group,element,length,dval)
       end

       'FD': begin
;Floating point double - 8 byte fp double val
;storing as LITTLE ENDIAN - needs to be checked!!!!
         dval=DOUBLE(data)
;explicily cast to bytes
         dvaltmp=BYTE(dval,0,8)
;fix byteorder (little-endian) if necessary; word length is 16 bits (?is this correct?)
         dval=[getbytes(UINT(dvaltmp[0:1],0,1),2), $
               getbytes(UINT(dvaltmp[2:3],0,1),2), $
               getbytes(UINT(dvaltmp[4:5],0,1),2), $
               getbytes(UINT(dvaltmp[6:7],0,1),2)]
         rs=generate_anytag(group,element,length,dval)
       end

       'IS': begin
;Decimal string - convert an int into a string, and store
;12 bytes maximum
         dval=STRTRIM(STRING(data),1)
         dval=STRMID(dval,0,12)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'LO': begin
;long string - IDL doesn't care about this one, 64 bytes max
         dval=STRMID(data,0,64)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'LT': begin
;long text - IDL doesn't care about this one too much, 10240 bytes max
         dval=STRMID(data,0,10240)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'OB': begin
;other byte string - padded by 00H
         dval=data
         rs=generate_anytag(group,element,length,dval)
       end

       'OW': begin
;other word string - padded by 00H. not sure if this is working
         dval=data
         rs=generate_anytag(group,element,length,dval)
       end

;     'PN': begin
;person name - not supported! (yet?)
;      print, 'PN currently unsupported!'
;      rs=BYTE(0)
;     end

       'SH': begin
;short string - 16 bytes max
         dval=STRMID(data,0,16)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'SL': begin
;signed long
         dval=getbytes(LONG(data),4)
         rs=generate_anytag(group,element,length,dval)
       end

;     'SQ': begin
;sequence of items - not supported!
;      print, 'SQ currently unsupported!'
;      rs=BYTE(0)
;     end

       'SS': begin
;signed short
         dval=getbytes(FIX(data),2)
         rs=generate_anytag(group,element,length,dval)
       end

       'ST': begin
;short text - 1024 bytes max
         dval=STRMID(data,0,1024)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'TM': begin
;time - of format hhmmss.frac, 16 bytes maximum
         dval=STRMID(data,0,16)
         rs=generate_anytag(group,element,length,dval,/STR)
       end

       'UI': begin
;unique identifier, 64 bytes maximum
         dval=STRMID(data,0,64)
         rs=generate_anytag(group,element,length,dval)
       end

       'UL': begin
;unsigned long
         dval = getbytes(ULONG(data[0]),4)
         for i=1, size(data, /n_elements)-1 do begin
          dval=[dval, getbytes(ULONG(data[i]),4)];
         endfor
         rs=generate_anytag(group,element,length,dval)
       end

       'UN': begin
;unknown - do whatever you please with this one
         dval=data
         rs=generate_anytag(group,element,length,dval)
       end

       'US': begin
;unsigned short
         dval = getbytes(UINT(data[0]),2)
         for i=1, size(data, /n_elements)-1 do begin
          dval=[dval, getbytes(UINT(data[i]),2)];
         endfor
         rs=generate_anytag(group,element,length,dval)
       end

       'UT': begin
;unlimited text; could be huge!
         dval=data
         rs=generate_anytag(group,element,length,dval)
       end
       else: begin
         dval=data
         rs=generate_anytag(group,element,length,dval)
       end
    ENDCASE
    return, BYTE(rs)
end

;generate pixel tag
function generate_pixeltag, group, element, val
    return, BYTE([getbytes(group,2),getbytes(element,2), getbytes(val,4)])
end


; write a tag into a buffer (byte array) and return it
function write_buff_tag, group, element, vr, length, value

    if(group ne '0002'x) then begin
       if(group ne '7fe0'x or element ne '0010'x) then begin
         if(length ne 0) then begin
;;          if(PTR_VALID(value)) then begin
          if((PTR_VALID(value))[0]) then begin;; array tags
              ret = generate_VRtag(group,element,vr,length,*value)
          endif else begin
              tmpvalue = bytarr(length);
              ret = generate_VRtag(group,element,vr,length,tmpvalue)
          endelse
         endif else begin
          ret = generate_empty_tag(group,element);
         endelse
       endif else begin  ; write the pixel data block
         ret = generate_pixeltag('7FE0'x,'0010'x,length)
         ret = [ret, *value]
       endelse
    endif else begin   ;write explicit VR for meta head
       if(length ne 0) then begin
         if(element eq '0010'x) then begin ; transfer syntax
          tmpvalue = '1.2.840.10008.1.2'; implicit little
          ret = generate_exp_VRtag(group,element,vr,length,tmpvalue)
         endif else begin
          if(PTR_VALID(value)) then begin
              ret = generate_exp_VRtag(group,element,vr,length,*value)
          endif else begin
              tmpvalue = bytarr(length);
              ret = generate_exp_VRtag(group,element,vr,length,tmpvalue)
          endelse
         endelse
       endif
    endelse

    return, ret;
end


; write a tag into a file
pro write_tag, U, group, element, vr, length, value

    if(group ne '0002'x) then begin
       if(group ne '7fe0'x or element ne '0010'x) then begin
         if(length ne 0) then begin
          if(PTR_VALID(value)) then begin
              WRITEU, U, generate_VRtag(group,element,vr,length,*value)
          endif else begin
              tmpvalue = bytarr(length);
              WRITEU, U, generate_VRtag(group,element,vr,length,tmpvalue)
          endelse
         endif else begin
          WRITEU, U, generate_empty_tag(group,element);
         endelse
       endif else begin  ; write the pixel data block
         WRITEU, U, generate_pixeltag('7FE0'x,'0010'x,length)
         WRITEU, U, *value
       endelse
    endif else begin   ;write explicit VR for meta head
       if(length ne 0) then begin
         if(element eq '0010'x) then begin ; transfer syntax
          tmpvalue = '1.2.840.10008.1.2'; implicit little
          WRITEU, U, generate_exp_VRtag(group,element,vr,length,tmpvalue)
         endif else begin
          if(PTR_VALID(value)) then begin
              WRITEU, U, generate_exp_VRtag(group,element,vr,length,*value)
          endif else begin
              tmpvalue = bytarr(length);
              WRITEU, U, generate_exp_VRtag(group,element,vr,length,tmpvalue)
          endelse
         endelse
       endif
    endelse
end


;; JACK's win_dicom_writer - it now has to take its input out of the
;;  header structure array insted of an existing dicom image - otherwise we will
;; needs some cleaning - e.g. should this thing modify any header tags?????
;;  or just write out valid headers

pro win_dicom_writer_data_kb_backup,  dcm_result_fn, new_image, dcm_header
;; and relicts from this routine
 series_num=1
 image_num=244
;; end of relicts

    dim = size(new_image, /dimensions);
    num_col = fix(dim[0]);
    num_row = fix(dim[1]);
    n_ele = size(new_image, /n_elements);

    itype = size(new_image, /type);
    if(itype eq 1) then begin
       bitsa=8
       bitss=8
       highbit=7
       n_buf = n_ele
    end;
    if(itype eq 2) then begin
       bitsa=16
       bitss=16
       highbit=15
       n_buf = n_ele * 2
    end;

    filename = dcm_result_fn
        OPENW, U, filename,/GET_LUN

        print, FORMAT='(%"Open for writing <%s> ")', filename

    ; if metaheader exist, should include metaheader
    group = dcm_header[0].group; test the first tag
    if(group[0] eq '0002'x) then begin
       tmpvalue = bytarr(128);
       WRITEU, U, tmpvalue;
       tmpvalue = 'DICM';
       WRITEU, U, tmpvalue;
    endif

    ; get system time string
    t = BIN_DATE(SYSTIME());
    YY = strtrim(t(0),2)
    MM = strtrim(t(1),2)
    if(strlen(MM) eq 1) then MM = '0'+MM;
    DD = strtrim(t(2),2)
    if(strlen(DD) eq 1) then DD = '0'+DD;
    HH = strtrim(t(3),2)
    if(strlen(HH) eq 1) then HH = '0'+HH;
    NN = strtrim(t(4),2)
    if(strlen(NN) eq 1) then NN = '0'+NN;
    SS = strtrim(t(5),2)
    if(strlen(SS) eq 1) then SS = '0'+SS;

    date = YY+MM+DD;
    time = HH+NN+SS;

    ; init the setting
    set_imaget = 0;
    set_cdate = 0;
    set_ctime = 0;
    set_idate = 0;
    set_itime = 0;
    set_sdate = 0;
    set_stime = 0;
    set_sdevice = 0;
    set_sversion = 0;
    set_snum = 0;
    set_inum = 0;
    set_sample = 0;
    set_photo = 0;
    set_bitsa = 0;
    set_bitss = 0;
    set_highbit = 0;
    set_pixelr = 0;
    set_small = 0;
    set_large = 0;
    set_width = 0;
    set_center = 0;
    set_inter = 0;
    set_slope = 0;
    set_stype = 0;
    set_plane = 0;
    set_sdes = 0;

    ; start to fill the tags
    last_group = 0
    group_buf = bytarr(2)

;; minimal invasive ::obj->GetReference()
        FOR i = 0L, N_ELEMENTS(dcm_header)-1 DO BEGIN
            group   =dcm_header[i].group
            element =dcm_header[i].element
            vr      =dcm_header[i].vr
            length  =dcm_header[i].length
;;    value   =*((dcm_header[i].pVal)[0])
            value   =dcm_header[i].pVal

            com_tag = long(group[0])*'10000'xL+long(element[0])
       cur_buf = bytarr(2)

       case com_tag of
                    '0008000'x: begin
                    end
;      '00080008'x: begin
;          end
;;       '00080016'x: begin
                    '7fe00000'x:begin
                    end
                    '7fe00010'x:begin
                    end
                    else: begin
                        if(element[0] ne 0) then begin
                            cur_buf=[cur_buf,write_buff_tag(group[0], element[0], vr[0], length[0], value[0])] ;
                        endif
                    end
                endcase

                size_buf = size(cur_buf,/n_elements)
                if(size_buf gt 2) then begin
                    if(size(group_buf,/n_elements) lt 4) then begin
                        group_buf=cur_buf[2:size_buf-1]
                    endif else begin
                        group_buf = [group_buf, cur_buf[2:size_buf-1]]
                    endelse
                endif

                if(element[0] ne 0) then last_group=group[0]
;;                print,FORMAT='(%"length of group_buf: %d")',N_ELEMENTS(group_buf)
            ENDFOR
;;            print,FORMAT='(%"win_dicom_writer_data_kb: WRITING ON Index: %5d (%4.4x %4.4x) ")',i,group,element

;; THE FOLLOWING TRIES TO WRITE SORT OF A PREAMBLE WHICH CAUSES IDL TO FAIL
;;            WRITEU, U, generate_VRtag(last_group, '0000'x, 'UL', 4, size(group_buf, /n_elements))

            ;; write the group buffer
            WRITEU, U, group_buf



                ;; seems to finish writing header on first
                ;; element with value 0

;;   if(element[0] eq 0 and last_group ne 0) then begin
;;                    print,FORMAT='(%"win_dicom_writer_data_kb: WRITING ON Index: %5d (%4.4x %4.4x) ")',i,group,element
;;       WRITEU, U, generate_VRtag(last_group, '0000'x, 'UL', 4, size(group_buf, /n_elements))
;;       WRITEU, U, group_buf
;;       group_buf = bytarr(2)
;;   endif
;;  Endfor

    WRITEU, U, generate_VRtag('7fe0'x,'0000'x,'UL',4,n_buf+8)
    WRITEU, U, generate_pixeltag('7FE0'x,'0010'x,n_buf)
    WRITEU, U, new_image

    CLOSE, U
    FREE_LUN, U

        print, FORMAT='(%"Successfully wrote <%s>with n_buf=%d ")', filename,n_buf
end

;; save fit parameters
FUNCTION read_fitpar,fn,comment,a,sa,tError
;; preset return status
text='' & tError='' & ret=0
a=DBLARR(4)
sa=DBLARR(4)

;; model is necessary, so we can check and take action on the fitparas

OPENR,unit,fn,/GET_LUN

readf,unit,text,FORMAT='(%"%s")'
comment=text

readf,unit,text,FORMAT='(%"%s")'
FOR i=0, 3 DO BEGIN
    readf,unit,tempa,tempb,FORMAT='(%"%10.4f %10.4f")'
    a[i]=tempa
    sa[i]=tempb
ENDFOR

FREE_LUN,unit
;; error handling
RETURN, (0)
END


;; save fit parameters
FUNCTION save_fitpar,fn,comment,a,sa,model,tError
;; preset return status
tError='' & ret=0

;; model is necessary, so we can check and take action on the fitparas

ON_IOERROR,bad_name

CASE model.name OF
    "TWO_KOMP_MODEL":BEGIN
        ;; convert tlag in seconds
        IF N_ELEMENTS(a) GE 4 THEN BEGIN
            a[3]  =  a[3] * 60. ; sec
            sa[3] = sa[3] * 60. ; sec
        ENDIF

        OPENW,unit,fn,/GET_LUN
        status=FSTAT(unit)
        IF comment NE '' THEN printf,unit,comment
        IF N_ELEMENTS(a) EQ 4 THEN BEGIN
            printf,unit,'!amp sd_amp, kep sd_kep kel sd_kel tlag sd_tlag'
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[0],sa[0])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[1],sa[1])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[2],sa[2])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[3],sa[3])
        ENDIF
        IF N_ELEMENTS(a) EQ 6 THEN BEGIN
            printf,unit,'!amp sd_amp, kep sd_kep kel sd_kel tlag sd_tlag'
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[0],sa[0])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[1],sa[1])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[2],sa[2])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[3],sa[3])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[4],sa[4])
            printf,unit,STRING(FORMAT='(%"%10.4f %10.4f")',a[5],sa[5])
        ENDIF
        FREE_LUN,unit
        PRINT, FORMAT='(%"Successfully wrote:<%s>")', fn
    END
    ELSE: BEGIN
        PRINT,FORMAT='(%"Unhandled model name <%s>")',model.name
        STOP
    END
ENDCASE
;; error handling
RETURN, (0)
bad_name:  tError=!ERROR_STATE.MSG
RETURN, (1)
END



;; save the intensity of the pixels inside the histogram parameters
FUNCTION save_intensity_histogram,fn,data,tError
;; preset return status
tError='' & ret=0

;; model is necessary, so we can check and take action on the fitparas

ON_IOERROR,bad_name

OPENW,unit,fn,/GET_LUN
FOR i=0L,N_ELEMENTS(data)-1 DO BEGIN
    printf,unit,STRING(FORMAT='(%"%10.1f")',data[i])
ENDFOR
FREE_LUN,unit
PRINT, FORMAT='(%"Successfully wrote:<%s>")', fn
;; error handling
RETURN, (0)
bad_name:  tError=!ERROR_STATE.MSG
RETURN, (1)
END



;; and some things from shortnam.pro
;; save_bmp is separated into four different modules ....
;; and I am not going to read bmps back right now
;; "BMP of Actual Image':"
;; save a bitmap of graphics window win_id
;; so we are able to keep it that short ....
PRO kb_save_bmp,filename,win_id,widgets
;;ON_IOERROR,bad_bmp_name

;; images are displayed from  0:bottom to top  or  1:top to bottom^

!order=0 ;;
; true-color-Bild vom Bildschirm lesen:
bmp_ima = tvrd(true=1)
!order=1

;; and write true color
WRITE_BMP,filename,bmp_ima

;; some place for error handling:
; IF 1 NE 0 THEN BEGIN
;     bad_bmp_name:
;     widget_control,widgets.sanduhr,SET_VALUE = 'Bad filename'
;     ON_IOERROR,NULL
;     print,string(7B),string(7B) ;bell
;     message,!ERR_STRING,/NONAME,/IOERROR,/CONTINUE
;     filename = get_text(TITLE="filename:",DEFAULT=filename)
;     IF filename EQ '' THEN save = 0 ; Abbruch
; ENDIF


END



;; this routine is placed here because of the need of having the
;; dicom-routines compiled before accessing this stuff

;; insert new entry in appropriate place in list
;; as the tag list is sorted in ascending order, this is easy
FUNCTION insert_tag_in_taglist,taglist,grp,elm,vr,value,CREATE=create,GET_LENGTH=get_length
;+
; NAME:insert_tag_in_taglist
;
;
;
; PURPOSE: insert dicom tags with elements into an existing dicom
; header structure
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


;; I (KB) would wish this completely different, but there is no way
;; to check whether an (only potentially) existing variable is a structure of a desired
;; kind - taglist should be checked and created, if no valid
;;        struc_dicom is passed


;; the bad thing is, that we have to create a correct length for the
;; TAG .... so I modified generate_VRtag to return the length of the
;; item if l is set to a value < 0. If a valid length is to be passed,
;; it has to be > 0....

    l=-1
    buffer=generate_VRtag( grp, elm, VR, l, value)

;;    print,FORMAT='(%"(get_length: %4.4x,%4.4x)  ",$)',grp,elm
;;    FOR i=0L,N_ELEMENTS(buffer)-1 DO BEGIN
;;        print,FORMAT='(%"%x",$)',buffer[i]
;;    END
;;print,' Length:',l ;; and print newline
;; ENDIF
if KEYWORD_SET(create) THEN BEGIN ;; -1 means var does not exist, and I won't create empty or less-than-two elements taglists
    taglist=replicate({struc_dicom},1)
    created=1
ENDIF ELSE BEGIN
    created=0
    ;;   goto, created
ENDELSE

;;if elm EQ 0 THEN stop

IF created EQ 1 THEN BEGIN
    ;; create first tag
    place=0
    taglist[place].group=grp
    taglist[place].element=elm
    taglist[place].vr=vr
    taglist[place].length=strlen(value);; KB invalid
    taglist[place].length=l
    taglist[place].pVal=PTR_NEW(value)
ENDIF ELSE BEGIN ;; check for create or replace
    place = WHERE ((taglist[*].group eq grp) AND (taglist[*].element eq elm))
    IF place[0] NE -1 THEN BEGIN
        ;; replace this tag
;;        print, place,grp,elm,FORMAT='(%"replaced element #%d (%04x,%04x)")'
        taglist[place].group=grp
        taglist[place].element=elm
        taglist[place].vr=vr
        taglist[place].length=strlen(value);; KB invalid
;;        IF KEYWORD_SET(GET_LENGTH)
        taglist[place].length=l
        ;; avoid memory leak by freeing the old pointer
        IF PTR_VALID(taglist[place].pVal) THEN PTR_FREE,taglist[place].pVal
        taglist[place].pVal=PTR_NEW(value)
        goto, done
    ENDIF ELSE BEGIN
        ;; create this tag
        ;; first find position of previous group tag
        place = WHERE ((taglist[*].group eq grp) AND (taglist[*].element lt elm))
        IF place[0] NE -1 THEN BEGIN
            place = place[N_ELEMENTS(place)-1] ;; first lesser elm
        ENDIF ELSE BEGIN
            place = WHERE (taglist[*].group lt grp)
            IF place[0] EQ -1 THEN BEGIN
                place=-1 ;; generate new first entry for the list
            ENDIF ELSE BEGIN
                place = place[N_ELEMENTS(place)-1];; last lesser elm
            ENDELSE
        ENDELSE
    ENDELSE
    ;; now insert element in location place + 1
;;    print, grp,elm,FORMAT='(%"inserted element (%04x,%04x)")'

    ;; get some properties of the old taglist
    ntags =  N_ELEMENTS(taglist)

    ;; first create taglist with one more tag
    newtaglist =replicate(taglist[0],ntags +1 )

    ;; copy the old contents to the new places
    IF place NE -1 THEN BEGIN
        newtaglist[0:place]=taglist[0:place]
    ENDIF
    IF place+2 LE ntags THEN BEGIN
        newtaglist[place+2:ntags]=taglist[place+1:ntags-1]
    ENDIF

    ;; and insert the new tag
    newtaglist[place+1].group=grp
    newtaglist[place+1].element=elm
    newtaglist[place+1].vr=vr
    newtaglist[place+1].length=l                     ;; KB invalid
    newtaglist[place+1].pVal=PTR_NEW(value)

    ;; replace the old taglist with the new one
    taglist = newtaglist & newtaglist = ''
ENDELSE


;;done: print,FORMAT='(%"Taglist has %d elements!")',N_ELEMENTS(taglist)
done:return, taglist
END


PRO free_taglist,taglist
;+
; NAME: free_taglist,taglist
;
;
;
; PURPOSE: free pointer to elements in a taglist and destroy the
; variable
FOR i=0L,N_ELEMENTS(taglist)-1 DO BEGIN
    IF PTR_VALID(taglist[i].pVal) THEN PTR_FREE,taglist[i].pVal
ENDFOR

;; and destroy the variable by handing its contents to a pointer
a=PTR_NEW(taglist,/NO_COPY)
;; and destroy the pointer
PTR_FREE,a
END



FUNCTION fake_dicom_header
;+
; NAME: fake_dicom_header
;
;
;
; PURPOSE: create a dicom header with the minimum number of entries
; necessary for postprocessing the images with IDL
;
; The tags are taken from the dicom_writer example available from the INTERNET
;
; This routine is especially necessary for creating DICOM output from
; proprietary input format e.g. SIEMENS VISION
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
; oct 16, 2002 baudendistel.4@osu.edu -creation
;-

; DICOM tags - feel free to add more!
; DICOM tags - feel free to add more!
;dummy fill-in variables.
random= '123456'
;SOP class set to MR - see Annex A in PS 3.6-2001
SOPClass = '1.2.840.10008.5.1.4.1.1.20'
SOPInstance = '1.2.840.10008.5.1.4.1.1.20.1'
StudyInstanceUID = SOPInstance + random
SeriesInstanceUID = StudyInstanceUID
RelFrameOfReferenceUID = StudyInstanceUID
SeriesInstanceUID = SeriesInstanceUID + '.1'
RelFrameOfReferenceUID = RelFrameOfReferenceUID + '.2'

StudyID = 0
Seriesnum = 0
Acqnum = 0
Imagenum = 0

rows=256
cols=256
bpp=2;; BYTES per pixel!!

thickness=1.0
spacing='1.0\\1.0'

StudyID = 0
Seriesnum = 0
Acqnum = 0
Imagenum = 0

patient='unknown'
physician='unknown'
patient_id='unknown'

;0008 tags

;MR type
    t=insert_tag_in_taglist(t,'0008'x,'0008'x,'CS','ORIGINAL\\PRIMARY\\OTHER',/CREATE)
;Instance date
    t=insert_tag_in_taglist(t,'0008'x,'0012'x,'DA','2002.01.14')
;Instance time
    t=insert_tag_in_taglist(t,'0008'x,'0013'x,'TM','150101')
;SOP class
    t=insert_tag_in_taglist(t,'0008'x,'0016'x,'UI',SOPClass)
;SOP instance
    t=insert_tag_in_taglist(t,'0008'x,'0018'x,'UI',SOPInstance)
;Modality
    t=insert_tag_in_taglist(t,'0008'x,'0060'x,'CS','MR')
;Manufacturer
    t=insert_tag_in_taglist(t,'0008'x,'0070'x,'LO','GE')
;Study Physicians Name
    t=insert_tag_in_taglist(t,'0008'x,'0090'x,'UI',physician)

;0010 tags

;Patient name
    t=insert_tag_in_taglist(t,'0010'x,'0010'x,'UI',patient)
;Patient ID
    t=insert_tag_in_taglist(t,'0010'x,'0020'x,'LO',patient_id)
;Patient birth date
    t=insert_tag_in_taglist(t,'0010'x,'0030'x,'DA','20020114')
;Patient sex
    t=insert_tag_in_taglist(t,'0010'x,'0040'x,'CS','M')

;0018 tags

;Acquisition type
    t=insert_tag_in_taglist(t,'0018'x,'0023'x,'CS','2D')
;Slice thickness
    t=insert_tag_in_taglist(t,'0018'x,'0050'x,'DS',thickness)

;0020 tags

;Study instance
    t=insert_tag_in_taglist(t,'0020'x,'000D'x,'UI',StudyInstanceUID)
;Series instance UID
    t=insert_tag_in_taglist(t,'0020'x,'000E'x,'UI',SeriesInstanceUID)
;StudyID
    t=insert_tag_in_taglist(t,'0020'x,'0010'x,'IS',StudyID)
;Series number
    t=insert_tag_in_taglist(t,'0020'x,'0011'x,'IS',seriesnum)
;Acquisition number
    t=insert_tag_in_taglist(t,'0020'x,'0012'x,'IS',acqnum)
;Image number
    t=insert_tag_in_taglist(t,'0020'x,'0013'x,'IS',imagenum)

;0028 tags

;samples per pixel
    t=insert_tag_in_taglist(t,'0028'x,'0002'x,'US',1)
;Photometric interpretation
    t=insert_tag_in_taglist(t,'0028'x,'0004'x,'CS','MONOCHROME2')
;Rows in image
    t=insert_tag_in_taglist(t,'0028'x,'0010'x,'US',rows)
;Columns in image
    t=insert_tag_in_taglist(t,'0028'x,'0011'x,'US',cols)
;pixel spacing
    t=insert_tag_in_taglist(t,'0028'x,'0030'x,'DS',spacing)
;bits allocated per sample
    t=insert_tag_in_taglist(t,'0028'x,'0100'x,'US',bpp*8)
;bits stored per sample
    t=insert_tag_in_taglist(t,'0028'x,'0101'x,'US',bpp*8)
;high bit
    t=insert_tag_in_taglist(t,'0028'x,'0102'x,'US',(bpp*8)-1)
;pixel representation
    t=insert_tag_in_taglist(t,'0028'x,'0103'x,'US','0000'x)
;; and an empty tag to close the whole thing off
;;        t=insert_tag_in_taglist(t,'0032'x,'0000'x,'US','0000'x)

return, (t)
END



;; remove tag index w from taglist
FUNCTION remove_taglist_entry,struc_tag,w
n=N_ELEMENTS(struc_tag)

;; free the value pointer
IF PTR_VALID(*(struc_tag[w].pVal[0])) THEN BEGIN
    PTR_FREE, struc_tag[w].pVal[0]
ENDIF

IF w NE n-1 THEN BEGIN
    t=[struc_tag[0:w-1],struc_tag[w+1:n-1]]
ENDIF ELSE BEGIN
    t=struc_tag[0:w-1]
ENDELSE
RETURN, t
END

;; and something to strip an existing header structure of problem tags
;; remove everything which JACK removed originally
FUNCTION strip_off_problem_tags,struc_tag
;; the list of header tags which JACK was removing in win_dicom_writer_data
list =['00080008'x,'00080012'x,'00080013'x,'00080016'x,'00080018'x, $
       '00080023'x,'00080033'x,'00181012'x,'00181014'x,'00181016'x, $
       '00181020'x,'0020000e'x,'00200011'x,'00200013'x,'00280011'x, $
       '00280100'x,'00280101'x,'00280102'x,'7fe00000'x,'7fe00010'x]

;; for self-consisteny eliminate possible double tags
list=list[uniq[list]]


;; as we assume that header is correctly sorted in ascending gr/el
;; just go backwards and eliminate entries
FOR i = N_ELEMENTS(struc_tag)-1 , 0 ,-1 DO BEGIN
    group   =struc_tag[i].group
    element =struc_tag[i].element
    com_tag= long(group)*'10000'xL+long(element)

    w=where[com_tag EQ list]
    if w[0] NE -1 THEN BEGIN
        dcm_header=remove_taglist_entry(struc_tag,w)
        struc_tag=dcm_header
    ENDIF
ENDFOR
RETURN, (struc_tag)
END

;; a routine for faking headers containing goodies from our
;; header structure
FUNCTION preset_header_from_info_struc,hs
;; generate the header completely from scratch
header=fake_dicom_header()

;; add the slice location ('0020'x,'1041'x)
header=insert_tag_in_taglist(header,'0020'xL,'1041'xL,'DS',$
                             STRING(hs.sl_pos))

;; add the patient name ('0010'x,'0010'x)
header=insert_tag_in_taglist(header,'0010'xL,'0010'xL,'PN',$
                             hs.pat_name)

;; add the patient ID ('0010'x,'0020'x)
header=insert_tag_in_taglist(header,'0010'xL,'0020'xL,'LO',$
                             hs.pat_id)

;; add the pixel spacing
;; if I interpret things right, then pixel spacing should be a
;; two dimensional array....
;; so if it is a 1d-array, simply double the first component
;; if it is 2d, then build one single string of it....
ps=hs.pixel_spacing
if N_ELEMENTS(ps) EQ 1 THEN BEGIN
    h=ps[0]
    ps=STRCOMPRESS(STRING(FORMAT='(%"%10.3f\\%10.3f")',h,h),/REMOVE_ALL)
ENDIF ELSE BEGIN
    h=STRCOMPRESS(STRING(FORMAT='(%"%10.3f\\%10.3f")',ps[0],ps[1]),/REMOVE_ALL)
    ps=h
ENDELSE
header=insert_tag_in_taglist(header,'0028'xL,'0030'xL,'DS',$
                             ps)

;; the slice thickness
header=insert_tag_in_taglist(header,'0018'xL,'0050'xL,'DS',$
                             hs.sl_thick)

;; and the following are copied from JACK
header=insert_tag_in_taglist(header,'0008'x, '0008'x, 'CS', 'DERIVED\SECONDARY\SCREEN SAVE\IP')


buff = 'IDL ' + !VERSION.RELEASE
header=insert_tag_in_taglist(header,'0018'x, '1016'x, 'LO',  buff)

buff = 'Dynamic Enhancement 1.1'
header=insert_tag_in_taglist(header,'0018'x, '1020'x, 'LO',  buff)
RETURN, header
END



;; and a routine for writing true color bmps as DICOM images
;; more or less copied from JACK

;; the STOP should be replaced by proper Error Handling (will you ever
;; do this, Klaus??)-started, but not complete
PRO dicom_writer_bmp,filename,truecolorbmp,dcm_header
tError=''
ON_IOERROR,error
;; just do it the way we want it in future
OPENW, U, filename,/GET_LUN

print, FORMAT='(%"Open for writing <%s> ")', filename


;; test for format RGBxWIDTHxHEIGHT
s=size(truecolorbmp)
IF s[0] NE 3 THEN BEGIN
    tError='Wrong number of dimensions: '+STRING(s[0])+'INSTEAD OF 3'
    goto, error
ENDIF
IF s[1] NE 3 THEN BEGIN
    tError='Wrong number of dimensions: '+STRING(s[0])+'INSTEAD OF 3'
    goto, error
END
cols=s[2]
rows=s[3]

;; and set the necessary output entries for DATA TYPE bitmap

;;Rows in image
t=insert_tag_in_taglist(dcm_header,'0028'x,'0010'x,'US',rows)

;;Columns in image
t=insert_tag_in_taglist(dcm_header,'0028'x,'0011'x,'US',cols)

buff = 'DERIVED\SECONDARY\SCREEN SAVE\IP'
dcm_header=insert_tag_in_taglist(dcm_header,'0008'x, '0008'x, 'CS', buff)

buff= 'Screen save '
dcm_header=insert_tag_in_taglist(dcm_header,'0008'x, '103e'x, 'LO',  buff)

buff = 'IDL ' + !VERSION.RELEASE
dcm_header=insert_tag_in_taglist(dcm_header,'0018'x, '1016'x, 'LO',  buff)

buff = 'Dynamic Enhancement 1.1'
dcm_header=insert_tag_in_taglist(dcm_header,'0018'x, '1020'x, 'LO',  buff)

;; samples per pixel
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0002'x, 'US', 3)

;; Photometric interpretation
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0004'x, 'CS', 'RGB ')

;; Planar configuration
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0006'x, 'US', 0)

;; bits allocated
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0100'x, 'US', 8)

;; bits stored
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0101'x, 'US', 8)

;; high bit
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0102'x, 'US', 7)

;; pixel representation
dcm_header=insert_tag_in_taglist(dcm_header,'0028'x, '0103'x, 'US', 0)

;0020 tags
 ;;the unique image identifier - just take the juldate
 uid=STRCOMPRESS(STRING(FORMAT='(%"%30.10f")',1.e6*SYSTIME(/JULIAN)),/REMOVE_ALL)
 StudyInstanceUID='1.2.840.10008.5.1.4.1.1.20.1'+uid
 SeriesInstanceUID = StudyInstanceUID


 RelFrameOfReferenceUID = StudyInstanceUID
 SeriesInstanceUID = SeriesInstanceUID + '.1'
 RelFrameOfReferenceUID = RelFrameOfReferenceUID + '.2'

;Study instance
 dcm_header=insert_tag_in_taglist(dcm_header,'0020'x,'000D'x,'UI',StudyInstanceUID)

;Series instance UID
 dcm_header=insert_tag_in_taglist(dcm_header,'0020'x,'000E'x,'UI',SeriesInstanceUID)

;; Smallest image pixel value ... skipped till 28.0x106-1052


;;
                                ; if metaheader exist, should include metaheader
group = dcm_header[0].group     ; test the first tag

if(group[0] eq '0002'x) then begin
    tmpvalue = bytarr(128)      ;
    WRITEU, U, tmpvalue         ;
    tmpvalue = 'DICM'           ;
    WRITEU, U, tmpvalue         ;
endif

;; process the true color bitmap

                                ; get system time string
t = BIN_DATE(SYSTIME())         ;
YY = strtrim(t(0),2)
MM = strtrim(t(1),2)
if(strlen(MM) eq 1) then MM = '0'+MM ;
DD = strtrim(t(2),2)
if(strlen(DD) eq 1) then DD = '0'+DD ;
HH = strtrim(t(3),2)
if(strlen(HH) eq 1) then HH = '0'+HH ;
NN = strtrim(t(4),2)
if(strlen(NN) eq 1) then NN = '0'+NN ;
SS = strtrim(t(5),2)
if(strlen(SS) eq 1) then SS = '0'+SS ;

date = YY+MM+DD                 ;
time = HH+NN+SS                 ;

                                ; start to fill the tags
last_group = 0
group_buf = bytarr(2)

;; minimal invasive ::obj->GetReference()
FOR i = 0L, N_ELEMENTS(dcm_header)-1 DO BEGIN
    group   =dcm_header[i].group
    element =dcm_header[i].element
    vr      =dcm_header[i].vr
    length  =dcm_header[i].length
    value   =dcm_header[i].pVal

    com_tag = long(group[0])*'10000'xL+long(element[0])
    cur_buf = bytarr(2)

    case com_tag of
        '0008000'x: begin
        end
        '7fe00000'x:begin
        end
        '7fe00010'x:begin
        end
        else: begin
            if(element[0] ne 0) then begin
                cur_buf=[cur_buf,write_buff_tag(group[0], element[0], vr[0], length[0], value[0])] ;
            endif
        end
    endcase

    size_buf = size(cur_buf,/n_elements)
    if(size_buf gt 2) then begin
        if(size(group_buf,/n_elements) lt 4) then begin
            group_buf=cur_buf[2:size_buf-1]
        endif else begin
            group_buf = [group_buf, cur_buf[2:size_buf-1]]
        endelse
    endif

    if(element[0] ne 0) then last_group=group[0]
;;                print,FORMAT='(%"length of group_buf: %d")',N_ELEMENTS(group_buf)
ENDFOR
print,FORMAT='(%"dicom_writer_bmp: WRITING ON Index: %5d (%4.4x %4.4x) ")',i,group,element

;; write the group buffer
WRITEU, U, group_buf

n_bytes=N_ELEMENTS(truecolorbmp)
WRITEU, U, generate_VRtag('7fe0'x,'0000'x,'UL',4,n_bytes+8)
WRITEU, U, generate_pixeltag('7FE0'x,'0010'x,n_bytes)
WRITEU, U, truecolorbmp

FREE_LUN, U

print, FORMAT='(%"Successfully wrote:<%s>")', filename
tError=''
RETURN
error:  tError=!ERROR_STATE.MSG
error2: RETURN
END




;; JACK's win_dicom_writer - it now has to take its input out of the
;;  header structure array insted of an existing dicom image
;;  and I removed the stripping and creation of some tags, because
;;  the header should be prepared prior to that by another routine

;; Nov 01, 2002: KB:: I decided that this routine has to set image
;; data type, number of rows and columns, etc. based on the pixeldata
;; to avoid failure due to mismatching entries in the prepared headers

pro win_dicom_writer_data_kb,  dcm_result_fn, new_image, dcm_header,tError
;; clear error status, if present
IF N_PARAMS() EQ 4 THEN tError=''
ON_IOERROR, ioerror

 ;;determine type of image and set bpp
imtype=size(new_image,/type)

 CASE (imtype) OF
     1:  bpp=1L                  ;byte
     2:  bpp=2L                  ;int
     3:  bpp=4L                  ;long int
     12: bpp=2L                  ;unsigned int
     13: bpp=4L                  ;unsigned long int
     ELSE: bpp=0L
 ENDCASE

 IF (bpp EQ 0L) THEN BEGIN
     PRINT, 'Only integer type images (byte/int/long) supported at this time, sorry!'
     RETURN
 ENDIF

;; get the number of rows and cols
 s=SIZE(new_image)
 IF s[0] NE 2 THEN BEGIN
     PRINT, 'Only 2-dimensional data supported at this time, sorry!'
     RETURN
 ENDIF
 rows=s[2] & hrows = rows
 cols=s[1] & hcols = cols

 byteswapped_image=new_image

 ;; I dont understand the purpose of the following lines
 ; if (bpp ge 2) then begin
;      little_endian = (BYTE(1, 0, 1))[0]
;      if (little_endian ne 1) then begin
;          print, 'Big endian architecture detected'
;          byteorder,byteswapped_image,/SWAP_IF_BIG_ENDIAN
;      endif else begin
;          print, 'Little endian architecture detected'
;          byteorder,byteswapped_image,/SWAP_IF_LITTLE_ENDIAN
;      endelse
;  endif

;; I think swapping this way works better
;; image will be little endian as little endian is transfer syntax byteorder
byteorder, byteswapped_image, /SWAP_IF_BIG_ENDIAN


 ;; manipulate the pixel - dependend header entries
 ;;0028 tags

 ;;samples per pixel
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0002'x,'US',1)

 ;;Photometric interpretation
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0004'x,'CS','MONOCHROME2')

 ;;Rows in image
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0010'x,'US',hrows)

 ;;Columns in image
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0011'x,'US',hcols)

 ;;bits allocated per sample
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0100'x,'US',bpp*8)

 ;;bits stored per sample
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0101'x,'US',bpp*8)

 ;;high bit
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0102'x,'US',(bpp*8)-1)

 ;;pixel representation
 dcm_header=insert_tag_in_taglist(dcm_header,'0028'x,'0103'x,'US','0000'x)

;0020 tags
 ;;the unique image identifier - just take the juldate
 uid=STRCOMPRESS(STRING(FORMAT='(%"%30.10f")',1.e6*SYSTIME(/JULIAN)),/REMOVE_ALL)
 StudyInstanceUID='1.2.840.10008.5.1.4.1.1.20.1'+uid
 SeriesInstanceUID = StudyInstanceUID

 RelFrameOfReferenceUID = StudyInstanceUID
 SeriesInstanceUID = SeriesInstanceUID + '.1'
 RelFrameOfReferenceUID = RelFrameOfReferenceUID + '.2'

;Study instance
 dcm_header=insert_tag_in_taglist(dcm_header,'0020'x,'000D'x,'UI',StudyInstanceUID)

;Series instance UID
 dcm_header=insert_tag_in_taglist(dcm_header,'0020'x,'000E'x,'UI',SeriesInstanceUID)

 filename = dcm_result_fn
 OPENW, U, filename,/GET_LUN

 print, FORMAT='(%"Open for writing <%s> ")', filename

                                ; if metaheader exist, should include metaheader
 group = dcm_header[0].group    ; test the first tag
 if(group[0] eq '0002'x) then begin
     tmpvalue = bytarr(128)     ;
     WRITEU, U, tmpvalue        ;
     tmpvalue = 'DICM'          ;
     WRITEU, U, tmpvalue        ;
 endif


                                ; start to fill the tags
 last_group = 0
 group_buf = bytarr(2)

;; minimal invasive ::obj->GetReference()
 FOR i = 0L, N_ELEMENTS(dcm_header)-1 DO BEGIN
     group   =dcm_header[i].group
     element =dcm_header[i].element
     vr      =dcm_header[i].vr
     length  =dcm_header[i].length
;;    value   =*((dcm_header[i].pVal)[0])
     value   =dcm_header[i].pVal

     com_tag = long(group[0])*'10000'xL+long(element[0])
     cur_buf = bytarr(2)

     case com_tag of
         '0008000'x: begin
         end
;      '00080008'x: begin
;          end
;;       '00080016'x: begin
         '7fe00000'x:begin
         end
         '7fe00010'x:begin
         end
         else: begin
             if(element[0] ne 0) then begin
                 cur_buf=[cur_buf,write_buff_tag(group[0], element[0], vr[0], length[0], value[0])] ;
             endif
         end
     endcase

     size_buf = size(cur_buf,/n_elements)
     if(size_buf gt 2) then begin
         if(size(group_buf,/n_elements) lt 4) then begin
             group_buf=cur_buf[2:size_buf-1]
         endif else begin
             group_buf = [group_buf, cur_buf[2:size_buf-1]]
         endelse
     endif

     if(element[0] ne 0) then last_group=group[0]
;;                print,FORMAT='(%"length of group_buf: %d")',N_ELEMENTS(group_buf)
 ENDFOR
;; print,FORMAT='(%"win_dicom_writer_data_kb: WRITING ON Index: %5d (%4.4x %4.4x) ")',i,group,element

;; THE FOLLOWING TRIES TO WRITE SORT OF A PREAMBLE WHICH CAUSES IDL TO FAIL
;;            WRITEU, U, generate_VRtag(last_group, '0000'x, 'UL', 4, size(group_buf, /n_elements))

 ;; write the group buffer
 WRITEU, U, group_buf

 ;; seems to finish writing header on first
 ;; element with value 0

;;   if(element[0] eq 0 and last_group ne 0) then begin
;;                    print,FORMAT='(%"win_dicom_writer_data_kb: WRITING ON Index: %5d (%4.4x %4.4x) ")',i,group,element
;;       WRITEU, U, generate_VRtag(last_group, '0000'x, 'UL', 4, size(group_buf, /n_elements))
;;       WRITEU, U, group_buf
;;       group_buf = bytarr(2)
;;   endif
;;  Endfor

 imsize=LONG(rows)*LONG(cols)* LONG(bpp)
 WRITEU, U, generate_VRtag('7fe0'x,'0000'x,'UL',4,imsize+8)
 WRITEU, U, generate_pixeltag('7FE0'x,'0010'x,imsize)
 WRITEU, U, byteswapped_image

 FREE_LUN, U

 print, FORMAT='(%"Successfully wrote <%s> with imsize=%d ")', filename,imsize
 ;; or am I allowed to use "return" here?
 IF 1 EQ 0 THEN BEGIN
     ioerror: IF N_PARAMS() EQ 4 THEN BEGIN
         tError=!ERROR_STATE.MSG
     ENDIF ELSE BEGIN
         print,"win_dicom_writer_data_kb: Error",!ERROR_STATE.MSG
     ENDELSE
 ENDIF
END

;@GJ, gjia@xidian.edu.cn, 2023/3/27, adding the information to MPI dcm
PRO insert_meta_info_MPI_dcm
  afn='C:\D_drive\MPI_Tianjie\LY_MPICT\mouse\mouse\001_liposomem47h3DMPI_2022-10-18_18-42-19\tomographic_combined_image\'
  afn = DIALOG_PICKFILE(PATH=afn,/MUST_EXIST, /DIRECTORY, TITLE='Select Series folder:')

  
 a = FILE_SEARCH(afn, '*.dcm', count=nrfile)
 
 
 afn_des='C:\D_drive\MPI_Tianjie\LY_MPICT\mouse\MPI\'
 afn_des = DIALOG_PICKFILE(PATH=afn,/MUST_EXIST, /DIRECTORY, TITLE='Select saving folder:')
 FILE_MKDIR, afn_des+'MPI\'
 afn_des = afn_des+'MPI\'
 ;fn_MPI='C:\D_drive\liposomem47h3DMPI_0021.dcm'
 ;fn_MPI_new = 'C:\D_drive\liposomem47h3DMPI_0021new5.dcm'
 
 fn_CT='C:\D_drive\00007_IM7.dcm'
 data_1 = bytarr(158)
 OPENR,unit,fn_CT,/get_lun
 READU, unit, data_1
 FREE_LUN, UNIT
 ;@GJ, 2023/3/27, the length including all 0002 elements
 ;data_1[140]=182
 
 FOR i=0, nrfile-1 DO BEGIN
   fn_MPI = a[i]
   print, 'i = ', i+1, 'out of ', nrfile
   data_2 = read_binary(fn_MPI, DATA_START=132)
   
   ; Use Count to get the number of nonzero elements:
   index = WHERE(data_2 EQ 8, count)
   ; Only subscript the array if it is safe:
   ;define the File Meta Group Element Length 0002,0000
   IF count NE 0 THEN BEGIN
    FOR j=0, count-1 DO BEGIN
      IF data_2[index[j]+1] EQ 0 THEN BEGIN
        data_1[140]=index[j]+14
        break
      ENDIF
    ENDFOR
   ENDIF

   ;Call procedure in Klaus' code
   fdecomp, fn_MPI, disk, dir, name, qual, version
   fn_MPI_new = afn_des+name+'.'+qual
   OPENW, U, fn_MPI_new,/GET_LUN
   WRITEU, U, data_1
   WRITEU, U, data_2
   FREE_LUN, U
   
   obj = OBJ_NEW('IDLffDICOM', fn_MPI_new)

   ; Get the row & column size of the image(s):
   rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
   cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
   OBJ_DESTROY, obj

   print, 'rows: ', rows
   print, 'cols: ', cols
 ENDFOR

; test=read_dicom(fn_MPI_new)
; iimage, test
; 
; obj = OBJ_NEW('IDLffDICOM', fn_MPI_new)
;
; ; Get the row & column size of the image(s):
; rows = *(obj->GetValue('0028'x,'0010'x,/NO_COPY))[0]
; cols = *(obj->GetValue('0028'x,'0011'x,/NO_COPY))[0]
; OBJ_DESTROY, obj
; 
; print, 'rows: ', rows
; print, 'cols: ', cols
END

PRO insert_meta_info_single_MPI_dcm

  fn_MPI='C:\D_drive\liposomem47h3DMPI_0021.dcm'
  fn_CT='C:\D_drive\00007_IM7.dcm'

  fn_MPI_new = 'C:\D_drive\liposomem47h3DMPI_0021new2.dcm'

  data_1 = bytarr(158)
  OPENR,unit,fn_CT,/get_lun
  READU, unit, data_1
  FREE_LUN, UNIT
  ;@GJ, 2023/3/27, the length including all 0002 elements
  ;data_1[140]=182

  data_2 = read_binary(fn_MPI, DATA_START=132)
  
  ; Use Count to get the number of nonzero elements:
  index = WHERE(data_2 EQ 8, count)
  ; Only subscript the array if it is safe:
  ;define the File Meta Group Element Length 0002,0000
  IF count NE 0 THEN BEGIN
    FOR j=0, count-1 DO BEGIN
      IF data_2[index[j]+1] EQ 0 THEN BEGIN
        data_1[140]=index[j]
        break
      ENDIF
    ENDFOR
  ENDIF
  
  OPENW, U, fn_MPI_new,/GET_LUN
  WRITEU, U, data_1
  WRITEU, U, data_2
  FREE_LUN, U

  test=read_dicom(fn_MPI_new)
  iimage, test

  obj = OBJ_NEW('IDLffDICOM', fn_MPI_new)

  ; Get the row & column size of the image(s):
  rows = *(obj->GetValue('0028'x,'0010'x))[0]
  cols = *(obj->GetValue('0028'x,'0011'x))[0]
  OBJ_DESTROY, obj

  print, 'rows: ', rows
  print, 'cols: ', cols

END


;; ====================================================================================
;; read_nifti
;;
;; ARGUMENTS
;; fname        Nifti filename (either .nii or a pair of .hdr and .img)
;;
;;
;; OPTIONS
;; header      to get a structure with the Nifti header on return
;;
;; just_header reads the header information only
;;
;; show        show header information on screen
;;
;; compress    reads gnuzip Nifti files
;;
;; flip_ver    performs a top/bottom flip on the image/volume
;;             ==> Note this is the default
;;
;; flip_hor    performs a horizontal flip on the volume (e.g. for HyperView)
;;
;;
;; requirements:
;;    - needs:  IDLffNifti__define.pro
;;
;; Author: Juergen Dammers, J.Dammers@fz-juelich.de
;; ====================================================================================
FUNCTION read_nifti, fname, header=hdr, just_header=just_header,show=show, compress=compress, $
  flip_ver=flip_ver, flip_hor=flip_hor


  ;; ------------------------------------------------------------------
  ;; check input
  ;; ------------------------------------------------------------------
  if (file_test(fname) ne 1) then return,-1
  if (strlen(fname)-3 eq strpos(fname,'.gz')) then compress=1  ; check for compressed file


  ;; ------------------------------------------------------------------
  ;; Create IDLffNifti Object
  ;; ------------------------------------------------------------------
  oNifti = obj_new('IDLffNifti')



  ;; ------------------------------------------------------------------
  ;; read header file
  ;; ------------------------------------------------------------------
  hdr = oNifti->ReadHeaderFile(fname, compress=compress)
  ;; show header ?
  if keyword_set(show) then begin
    print, '>>> filename:'
    print, fname
    oNifti->OutputHeader, hdr,-1
  endif
  ;; return ?
  if keyword_set(just_header) then begin
    obj_destroy, oNifti
    return, hdr
  endif


  ;; ------------------------------------------------------------------
  ;; read data and header
  ;; obj->ReadFile, [f] [,_EXTRA=extra])
  ;; Parameters:
  ;; fname       filename  (hdr, nii)
  ;; ------------------------------------------------------------------
  oNifti->ReadFile, fname, compress=compress
  vol = oNifti->GetData()
  ;; Note for IDL we need to flip top/bottom ==> this is the default
  if not keyword_set(flip_ver) then vol = reverse(vol,2,/overwrite)

  ;; check if we need a left right flip (needed for HyperView)
  if keyword_set(flip_hor) then begin
    case hdr.dim[0] of
      2:  vol = reverse(vol,1,/overwrite)
      3:  vol = reverse(vol,3,/overwrite)
      else: stop
    endcase
  endif

  ;; ------------------------------------------------------------------
  ;; free pointer
  ;; ------------------------------------------------------------------
  obj_destroy, oNifti


  return, vol

END

;; ====================================================================================
;;
;;  save_nifti
;;
;; ARGUMENTS
;;   fnifti:    output filename of nifti header
;;   vol   :    3D volume to save (optional)
;;
;; OPTIONS
;;   header:  use this header as a header template
;;            (otherwise defaults are used as stored in IDLffNifti__define.pro)
;;
;;   just_header: store only header file, no volume will be stored
;;       - dimension: dimension of the volume (used togehter with just_header)
;;       - byte_type: idl byte type of the data (used togehter with just_header)
;;
;;   show:    show header info on screen after saving the data
;;
;;   dimension:  user defined image/volume dimension
;;   pixdim:  dimension of the voxel in mm
;;            default is pixdim = [1.0,1.0,1.0]
;;
;;   origin:  set new origin
;;
;;   vox_offset:  set new voxel offset
;;   quatern_bcd: set quaterions b,c,d
;;   smatrix    : instead of quatern_bcd you can provide the full 3x4 srow matrix
;;                example:
;;                smatrix = [ $
;;                          [ [-1, 0, 0] * pixdim[0], origin[0] ], $ ; left - right  orientation (radiological)
;;                          [ [ 0, 0,-1] * pixdim[1], origin[1] ], $ ; inferior - superior orientation  (from bottom to top)
;;                          [ [ 0, 1, 0] * pixdim[2], origin[2] ]]   ; posterior - anterior
;;
;;   noreverse:  do not flip top-bottom image orientation
;;               ==> the default is to flip the image because of IDL's image origin
;;
;; requirements:
;;    - needs:  IDLffNifti__define.pro
;;
;; Author: Juergen Dammers, J.Dammers@fz-juelich.de
;; ====================================================================================
FUNCTION save_nifti, fnifti, vol, header=hdr, just_header=just_header, dimension=dimension, show=show, $
  pixdim=pixdim, origin=origin, vox_offset=vox_offset,quatern_bcd=quatern_bcd, $
  noreverse=noreverse, smatrix=smatrix


  ;; ------------------------------------------------------------------
  ;; check input
  ;; ------------------------------------------------------------------
  sz  = size(vol)
  sz_dim = size(vol,/dim)
  if (sz_dim[0] ne 0) then ndim = n_elements(sz_dim) else ndim=0
  if (total(sz) lt 3) then just_header=1
  ext   = strmid(fnifti,strlen(fnifti)-3,3)
  if (ext eq 'nii') then begin
    if not keyword_set(vox_offset) then vox_offset=352
  endif else vox_offset=0
  if (ndim eq 4) then nvol = sz_dim[3] else nvol=1


  ;; ------------------------------------------------------------------
  ;; create Nifti Object
  ;; ------------------------------------------------------------------
  oNifti = obj_new('IDLffNifti')



  ;; ------------------------------------------------------------------
  ;; pass input header
  ;; ------------------------------------------------------------------
  if keyword_set(hdr) then begin
    header = hdr
    oNifti -> SetProperty,hdr
  endif else header = oNifti -> GetProperty()


  ;; ------------------------------------------------------------------
  ;; check keywords
  ;; ------------------------------------------------------------------
  ;; dimension
  if keyword_set(dimension) then dim=dimension else begin
    if (ndim gt 1) then begin
      dim         = make_array(8,/integer,value=1)
      dim[0]      = ndim
      dim[1:ndim] = sz_dim[0:ndim-1]
      dim[4]      = nvol
    endif else dim = header.dim
  endelse
  dim_ok = onifti->isnum(dim,/Positive) * ndim
  if (dim_ok lt ndim) then begin
    tmp = dialog_message('ERROR (save_nifti): Invalid dimension: must be 3, 4 or 5.')
    return,-1
  endif
  oNifti -> SetDim, dim

  ;; origin
  if keyword_set(origin) then begin
    orig_ok = (onifti->isnum(origin)) * n_elements(origin)
    origin  =  float(origin)
    if (orig_ok eq 3) then oNifti -> SetOrigin, origin
  endif

  ;; pixdim
  if keyword_set(pixdim) then begin
    ok = (onifti->isnum(pixdim)) * n_elements(pixdim)
    pixdim  =  float(pixdim)
    pixdim_array = header.pixdim
    pixdim_array[1:3] = pixdim
    if (ok eq 3) then oNifti -> SetPixdim, pixdim_array
  endif

  ;; sform matrix
  if keyword_set(smatrix) then begin
  ;@GJ, 2023/5/19, error may happen when saving nifti
    oNifti->CalcQuatern,smatrix=smatrix
  endif

  ;; quatern_bcd
  if keyword_set(quatern_bcd) then oNifti->SetQuatern_bcd, float(quatern_bcd)

  ;; extension
  if (ext eq 'nii') and (header.vox_offset eq 0) then oNifti -> SetVoxoffset,vox_offset



  ;; ------------------------------------------------------------------
  ;; pass data
  ;; ------------------------------------------------------------------
  if not keyword_set(just_header) then begin
    ;; check if we need to flip top -> bottom
    if not keyword_set(noreverse) then volume = reverse(vol,2) else volume=vol
    oNifti -> SetData, volume
    header =  oNifti -> GetProperty()
  endif



  ;; ------------------------------------------------------------------
  ;; write data
  ;; and get header on return
  ;; ------------------------------------------------------------------
  oNifti -> WriteFile, fnifti, just_header=just_header
  header = oNifti -> GetProperty()


  ;; ------------------------------------------------------------------
  ;; show header info
  ;; ------------------------------------------------------------------
  if keyword_set(show) then begin
    print, '>>> filename:'
    print, fnifti
    oNifti->OutputHeader, header,-1
  endif

  obj_destroy, oNifti

  return,header
END
