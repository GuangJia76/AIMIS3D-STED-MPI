; $Id: //depot/Release/ENVI51_IDL83/idl/idldir/examples/demo/demosrc/demo_filepath.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;

FUNCTION demo_filepath, filename, ROOT_DIR=root_dir, SUBDIRECTORY=subdirectory

;+
; NAME:
;     DEMO_FILEPATH
;
; PURPOSE:
;     Given the name of a file, DEMO_FILEPATH
;     returns the fully-qualified path to use in
;     opening the file. Unlike the standard FILEPATH routine,
;     DEMO_FILEPATH looks first in the local directory for the
;     file.  If the file is found there, DEMO_FILEPATH returns
;     the local name.  This is useful if the demo routine is to
;     be run with local data files, allowing a user to supply a
;     file without having to place it in the IDL distribution.
;
; CATEGORY:
;     Demo System, File Management.
;
; CALLING SEQUENCE:
;     Result = DEMO_FILEPATH('filename')
;
; INPUTS:
;     filename:   The name of the file to be opened. No device
;        or directory information should be included.
;
; KEYWORDS:
;     ROOT_DIR: The name of the directory from which the resulting path
;        should be based. If not present, and if the supplied filename is
;        not found in the current directory, the value of !DIR is used.
;
;        If the input filename is found in the current directory, then
;        keyword ROOT_DIR is ignored.
;
;     SUBDIRECTORY: The name of the subdirectory in which the file
;        should be found. This variable can be either a scalar
;        string or a string array with the name of each level of
;        subdirectory depth represented as an element of the array.
;
;        If the input filename is found in the current directory, then
;        keyword SUBDIRECTORY is ignored.
;
; OUTPUTS:
;     The fully-qualified file path is returned.
;
; COMMON BLOCKS:
;     None.
;
; EXAMPLE:
;     To get a path to the file cities.dat in the "examples/demo/demodata"
;     subdirectory of the IDL directory, enter:
;
;     path = DEMO_FILEPATH("cities.dat", $
;                          SUBDIRECTORY = ["examples", "demo", "demodata"])
;
;     The variable "path" contains a string that is the fully-qualified file
;     path for the file cities.dat.  Note that if the specified filename is
;     found in the current directory the full path to it will be returned.
;
; MODIFICATION HISTORY:
;     January, 1998, ACY and PCS, New for 5.1 demo system.
;
;-

if filename ne '' then $
   void = FILE_SEARCH(filename, count=count) $
else $
   count = 0

if (count gt 0) then begin
   cd, current=current
   case !version.os_family of
      'Windows': begin
         len = strlen(current)
;
;        Trim trailing backslash.
;
         if strpos(current, '\') eq (len - 1) then begin
            current = strmid(current, 0, len - 1)
            end
;
         end
      else:
      endcase
   return, filepath(filename, root_dir=current, subdirectory='')
   end $
else begin
   return, filepath(filename, root_dir=root_dir, subdirectory=subdirectory)
   end

end
