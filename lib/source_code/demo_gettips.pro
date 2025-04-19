; $Id: //depot/Release/ENVI51_IDL83/idl/idldir/examples/demo/demosrc/demo_gettips.pro#1 $
;
;  Copyright (c) 1997-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;----------------------------------------------------------
;
;    PURPOSE  Read a file with a specified format (see below) and
;             return the tip structure.
;
;             File format :
;
;             ncolumns nrows
;             id1 | text1
;             id2 | text2
;              ..    ..
;             idN | textN
;
;             ncolumns: number of columns for tip display
;             nrows   : number of rows of tip display
;             idX     : a label used to refer to a tip
;             textX   : the text of the corresponding tip
;             |       : the separator string. This can be specified
;                       with the SEPARATOR keyword.
;
;             The text associated with the first (ncolumns*nrows) tags
;             will be displayed in the text widgets at application startup.
;
;    RETURNS  A structure containing the tip information, including
;             the ids of the text widgets, number of columns and rows,
;             the tags and the text strings of the tips and the text
;             array currently loaded into the text widgets.
;

function demo_gettips, $
    filename, $             ; IN: string containing the tip file name
    wTopBase, $             ; IN: top level base
    wTipBase, $             ; IN: base which will hold the tip widgets
    SEPARATOR = separator   ; IN: (opt) single character used as separator.

    if (N_ELEMENTS(separator) EQ 0) then separator = '|'
    OPENR, lun, filename, /GET_LUN

    nrows=0 & ncols = 0
    READF, lun,  ncols, nrows

    tags = strarr(100)
    text = strarr(100)

    i=0
    while not eof(lun) do begin     ;Each message
        str = ''
        READF, lun,  str
        str = STRTOK(str, separator, /EXTRACT, /PRESERVE_NULL)
        tags[i] = str[0]
        text[i] = str[1]
        i = i + 1
    endwhile
    ntags = i
    FREE_LUN, lun

    ; Construct and set the initial values of the text widgets here.
    wTextId = lonarr(ncols)
    contents = strarr(ncols, nrows)
    for i=0, (ntags < ncols * nrows)-1 do $
      contents[i / nrows, i mod nrows] = text[i]

    for i = 0, ncols-1  do $
        wTextId[i] = WIDGET_TEXT(wTipBase, $
            XSIZE=43, $                  ; based on trial and error
                                         ; on 640x480 laptop screen
            YSIZE=nrows)

    for i=0, ncols-1 do widget_control, wtextid[i], set_val=contents[i,*]

    widget_control, wTipBase, /realize

    ;  resize text widgets if they can be larger
    tlbGeometry = WIDGET_INFO(wTopBase, /GEOMETRY)
    tlbXSize = tlbGeometry.xsize + (2 * tlbGeometry.margin)
    newTipXsize = tlbXsize/ncols - 5
    tipGeometry = WIDGET_INFO(wTextId[0], /GEOMETRY)
    if (tipGeometry.scr_xsize LT newTipXsize) then $
       for i=0, ncols-1 do $
          WIDGET_CONTROL, wTextId[i], SCR_XSIZE = newTipXsize

    widget_control, wTipBase, /map

    Result = { wTextId: wTextId, $ ;Widget ID of the text widgets
           wTipBase: 0L, $        ;Widget id of top level base. needed???
           ncols: ncols, $        ;# of columns of tips on screen
           nrows: nrows, $        ;# of rows of tips
           tags: tags[0:ntags-1], $ ;The tags from the tip file
           text: text[0:ntags-1], $ ;Text associated with each tag
           contents: contents     $ ;Current contents of text widget
        }

    RETURN, result

end

