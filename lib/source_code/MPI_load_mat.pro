;+
; NAME:
;
;   LOAD_MAT
;
; PURPOSE:
;
;   Read MATLAB MAT-files in IDL (see README for more information).
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
;   PRO load_mat, <filename>, <path>, STORE_LEVEL=store_level, $
;                 VERBOSE=verbose, DEBUG=debug
;
; MODIFICATION HISTORY:
;
;   See changelog.
;
; COPYRIGHT:
;
;   Copyright (C) 2009 Gordon Farquharson <gordonfarquharson@gmail.com>
;
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;      
;   You should have received a copy of the GNU General Public License
;   along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;-

FUNCTION size_of_data_type, data_symbol

    SWITCH data_symbol OF
        'miINT8'   :
        'miUINT8'  :
        'miUTF8'   : return, 1
        'miINT16'  :
        'miUINT16' :
        'miUTF16'  : return, 2
        'miINT32'  :
        'miUINT32' :
        'miUTF32'  :
        'miSINGLE' : return, 4
        'miINT64'  :
        'miUINT64' :
        'miDOUBLE' : return, 8
    ENDSWITCH

END

PRO skip_padding_bytes, lun, DEBUG=debug

    ;; All data elements are aligned on a 64 bit boundary. Calculate
    ;; how many padding bytes exist, and advance the file pointer
    ;; appropriately.

    point_lun, -lun, position

    IF (position MOD 8) NE 0 THEN BEGIN

        number_OF_padding_bytes = 8 - (position MOD 8)

        IF keyword_set(debug) THEN $
            print, 'Skipping ', number_of_padding_bytes, ' bytes'

        padding_bytes = bytarr(number_of_padding_bytes)
        readu, lun, padding_bytes

    ENDIF

END

PRO read_int8_data, lun, element_tag, data

    ;; FIXME: Not sure how to represent signed 8-bit data in IDL.

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = bytarr(number_of_elements)
    readu, lun, data

END

PRO read_uint8_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = bytarr(number_of_elements)
    readu, lun, data

END

PRO read_int16_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = intarr(number_of_elements)
    readu, lun, data

END

PRO read_uint16_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = uintarr(number_of_elements)
    readu, lun, data

END

PRO read_int32_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = lonarr(number_of_elements)
    readu, lun, data

END

PRO read_uint32_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = ulonarr(number_of_elements)
    readu, lun, data

END

PRO read_single_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = fltarr(number_of_elements)
    readu, lun, data

END

PRO read_double_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = dblarr(number_of_elements)
    readu, lun, data

END

PRO read_int64_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = lon64arr(number_of_elements)
    readu, lun, data

END

PRO read_uint64_data, lun, element_tag, data

    number_of_elements = element_tag.number_of_bytes / $
                         size_of_data_type(element_tag.data_symbol)
    data = ulon64arr(number_of_elements)
    readu, lun, data

END

PRO read_utf8_data, lun, element_tag, data

    data = bytarr(element_tag.number_of_bytes)
    readu, lun, data

END

PRO read_utf16_data, lun, element_tag, data

    data = intarr(element_tag.number_of_bytes)
    readu, lun, data

END

PRO read_utf32_data, lun, element_tag, data

    data = lonarr(element_tag.number_of_bytes)
    readu, lun, data

END

PRO skip_unknown_data, lun, element_tag

    data_bytes = bytarr(element_tag.number_of_bytes)
    readu, lun, data_bytes

END

FUNCTION element_tag_struct

    return, { mat_v5_element_tag, $
              data_type             : 0UL, $
              data_type_description : '', $
              data_symbol           : '', $
              number_of_bytes       : 0UL, $
              small_element_format  : 0B $
            }

END

PRO read_element_tag, lun, element_struct, DEBUG=debug

    data_type = 0UL
    number_of_bytes = 0UL

    readu, lun, data_type

    IF (data_type AND 'FFFF0000'XUL) EQ 0UL THEN BEGIN

        readu, lun, number_of_bytes

        element_struct.data_type = data_type
        element_struct.number_of_bytes = number_of_bytes
        element_struct.small_element_format = 0B

    ENDIF ELSE BEGIN

        ;; Small data element format

        element_struct.number_of_bytes = $
            ishft(data_type AND 'FFFF0000'XUL, -16)
        element_struct.data_type = data_type AND '0000FFFF'XUL
        element_struct.small_element_format = 1B

    ENDELSE

    data_type_description = ''
    data_symbol = ''

    CASE element_struct.data_type OF
        1  : BEGIN
            data_type_description = '8 bit, signed'
            data_symbol = 'miINT8'
        END
        2  : BEGIN
            data_type_description = '8 bit, unsigned'
            data_symbol = 'miUINT8'
        END
        3  : BEGIN
            data_type_description = '16 bit, signed'
            data_symbol = 'miINT16'
        END
        4  : BEGIN
            data_type_description = '16 bit, unsigned'
            data_symbol = 'miUINT16'
        END
        5  : BEGIN
            data_type_description = '32 bit, signed'
            data_symbol = 'miINT32'
        END
        6  : BEGIN
            data_type_description = '32 bit, unsigned'
            data_symbol = 'miUINT32'
        END
        7  : BEGIN
            data_type_description = 'IEEE 754 single format'
            data_symbol = 'miSINGLE'
        END
        8  : BEGIN
            data_type_description = 'Reserved (8)'
            data_symbol = ''
        END
        9  : BEGIN
            data_type_description = 'IEEE 754 double format'
            data_symbol = 'miDOUBLE'
        END
        10 : BEGIN
            data_type_description = 'Reserved'
            data_symbol = ''
        END
        11 : BEGIN
            data_type_description = 'Reserved'
            data_symbol = ''
        END
        12 : BEGIN
            data_type_description = '64 bit, signed'
            data_symbol = 'miINT64'
        END
        13 : BEGIN
            data_type_description = '64 bit, unsigned'
            data_symbol = 'miUINT64'
        END
        14 : BEGIN
            data_type_description = 'MATLAB array'
            data_symbol = 'miMATRIX'
        END
        15 : BEGIN
            data_type_description = 'Compressed data'
            data_symbol = 'miCOMPRESSED'
        END
        16 : BEGIN
            data_type_description = 'Unicode UTF-8 encoded character data'
            data_symbol = 'miUTF8'
        END
        17 : BEGIN
            data_type_description = 'Unicode UTF-16 encoded character data'
            data_symbol = 'miUTF16'
        END
        18 : BEGIN
            data_type_description = 'Unicode UTF-32 encoded character data'
            data_symbol = 'miUTF32'
        END
    ENDCASE

    element_struct.data_type_description = data_type_description
    element_struct.data_symbol = data_symbol

    IF element_struct.small_element_format THEN $
        small_element_text = 'True' $
    ELSE $
        small_element_text = 'False'

    IF keyword_set(DEBUG) THEN BEGIN
        print, 'Data type       : ', element_struct.data_type_description
        print, 'Data symbol     : ', element_struct.data_symbol
        print, 'Number of bytes : ', element_struct.number_of_bytes
        print, 'Small element   : ', small_element_text
    ENDIF

END

FUNCTION subelement_array_flags_struct

    return, { mat_v5_subelement_array_flags, $
              flag_word_1 : 0UL, $
              flag_word_2 : 0UL, $
              complex     : 0B, $
              global      : 0B, $
              logical     : 0B, $
              class       : 0B, $
              class_description : '', $
              class_symbol      : '' $
            }

END

PRO read_subelement_array_flags, lun, subelement_tag, subelement_struct, $
                                 DEBUG=debug

    flags1 = 0UL
    flags2 = 0UL

    readu, lun, flags1, flags2

    subelement_struct.flag_word_1 = flags1
    subelement_struct.flag_word_2 = flags2
    subelement_struct.complex = flags1 AND '00000800'XL
    subelement_struct.global = flags1 AND '00000400'XL
    subelement_struct.logical = flags1 AND '00000200'XL
    subelement_struct.class = flags1 AND '000000FF'XL

    IF keyword_set(debug) THEN BEGIN
        print, 'Complex           : ', subelement_struct.complex
        print, 'Global            : ', subelement_struct.global
        print, 'Logical           : ', subelement_struct.logical
    ENDIF

    class_description = ''
    class_symbol = ''

    CASE subelement_struct.class OF
        1 : BEGIN
            class_description = 'Cell array'
            class_symbol = 'mxCELL_CLASS'
        END
        2 : BEGIN
            class_description = 'Structure'
            class_symbol = 'mxSTRUCT_CLASS'
        END
        3 : BEGIN
            class_description = 'Object'
            class_symbol = 'mxOBJECT_CLASS'
        END
        4 : BEGIN
            class_description = 'Character array'
            class_symbol = 'mxCHAR_CLASS'
        END
        5 : BEGIN
            class_description = 'Sparse array'
            class_symbol = 'mxSPARSE_CLASS'
        END
        6 : BEGIN
            class_description = 'Double precision array'
            class_symbol = 'mxDOUBLE_CLASS'
        END
        7 : BEGIN
            class_description = 'Single precision array'
            class_symbol = 'mxSINGLE_CLASS'
        END
        8 : BEGIN
            class_description = '8-bit, signed integer'
            class_symbol = 'mxINT8_CLASS'
        END
        9 : BEGIN
            class_description = '8-bit, unsigned integer'
            class_symbol = 'mxUINT8_CLASS'
        END
        10 : BEGIN
            class_description = '16-bit, signed integer'
            class_symbol = 'mxINT16_CLASS'
        END
        11 : BEGIN
            class_description = '16-bit, unsigned integer'
            class_symbol = 'mxUINT16_CLASS'
        END
        12 : BEGIN
            class_description = '32-bit, signed integer'
            class_symbol = 'mxINT32_CLASS'
        END
        13 : BEGIN
            class_description = '32-bit, unsigned integer'
            class_symbol = 'mxUINT32_CLASS'
        END
    ENDCASE

    subelement_struct.class_description = class_description
    subelement_struct.class_symbol = class_symbol

    IF keyword_set(debug) THEN BEGIN
        print, 'Class description : ', subelement_struct.class_description
        print, 'Class symbol      : ', subelement_struct.class_symbol
    ENDIF

END

FUNCTION subelement_dimensions_array_struct

    ;; I think that IDL allows a maximum of 8 dimensions.

    return, { mat_v5_subelement_dimensions_array, $
              number_of_dimensions : 0L, $
              dimensions           : lonarr(8) $
            }

END

PRO read_subelement_dimensions_array, lun, subelement_tag, subelement_struct, $
                                      DEBUG=debug

    number_of_dimensions = subelement_tag.number_of_bytes / $
                           size_of_data_type(subelement_tag.data_symbol)

    subelement_struct.number_of_dimensions = number_of_dimensions
    
    ;; I don't know if this case statement is necessary. The
    ;; documentation is not clear on whether the dimensions array type
    ;; is always miINT32.

    dimensions = lonarr(number_of_dimensions)

    CASE size_of_data_type(subelement_tag.data_symbol) OF
        1 : BEGIN
            dimension = 0B
            FOR i = 0, number_of_dimensions-1 DO BEGIN
                readu, lun, dimension
                dimensions[i] = dimension
            ENDFOR
        END
        2 : BEGIN
            dimension = 0
            FOR i = 0, number_of_dimensions-1 DO BEGIN
                readu, lun, dimension
                dimensions[i] = dimension
            ENDFOR
        END
        4 : BEGIN
            dimension = 0L
            FOR i = 0, number_of_dimensions-1 DO BEGIN
                readu, lun, dimension
                dimensions[i] = dimension
            ENDFOR
        END
    ENDCASE

    subelement_struct.dimensions = dimensions

    IF keyword_set(debug) THEN BEGIN
        print, 'Number of dimensions : ', subelement_struct.number_of_dimensions
        print, 'Dimensions           : ', subelement_struct.dimensions
    ENDIF

    skip_padding_bytes, lun, DEBUG=debug

END

PRO read_subelement_array_name, lun, subelement_tag, array_name, DEBUG=debug

    ;; Assume that data type is always miINT8.

    array_name_bytes = bytarr(subelement_tag.number_of_bytes)
    readu, lun, array_name_bytes
    array_name = string(array_name_bytes)

    IF keyword_set(debug) THEN print, 'Array name : ', array_name

    skip_padding_bytes, lun, DEBUG=debug

END

PRO read_element_data, lun, element_tag, data, DEBUG=debug

    data_recognized = 1

    SWITCH element_tag.data_symbol OF

        'miINT8'       : BEGIN
            read_int8_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUINT8'      : BEGIN
            read_uint8_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miINT16'      : BEGIN
            read_int16_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUINT16'     : BEGIN
            read_uint16_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miINT32'      : BEGIN
            read_int32_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUINT32'     : BEGIN
            read_uint32_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miSINGLE'     : BEGIN
            read_single_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miDOUBLE'     : BEGIN
            read_double_data, lun, element_tag, data
            BREAK
        END

        'miINT64'      : BEGIN
            read_int64_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUINT64'     : BEGIN
            read_uint64_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miMATRIX'     :
        'miCOMPRESSED' : BEGIN
            print, '*** ', element_tag.data_symbol, ' NOT IMPLEMENTED ***'
            skip_unknown_data, lun, element_tag
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUTF8'       : BEGIN
            read_utf8_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUTF16'      : BEGIN
            read_utf8_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END

        'miUTF32'      : BEGIN
            read_utf8_data, lun, element_tag, data
            skip_padding_bytes, lun, DEBUG=debug
            BREAK
        END
        
        ELSE           : BEGIN
            skip_unknown_data, lun, element_tag
            skip_padding_bytes, lun, DEBUG=debug
            data_recognized = 0
        END

    ENDSWITCH

    skip_padding_bytes, lun, DEBUG=debug

    IF keyword_set(debug) THEN BEGIN
        IF data_recognized THEN $
            print, 'Data : ', data $
        ELSE $
            print, 'UNKNOWN DATA ELEMENT'
    ENDIF

END

FUNCTION format_array_element_data, data, array_flags, dimensions_array

    dimensions = $
        dimensions_array.dimensions[0:dimensions_array.number_of_dimensions-1]

    ;; Prevent 1x1 arrays from being created.

    IF size(data, /N_ELEMENTS) NE 1 THEN $
        _data = reform(reform(data, dimensions)) $
    ELSE $
        _data = data[0]

    CASE array_flags.class_symbol OF

        'mxCELL_CLASS'   : BEGIN
            print, '*** Formatting ', array_flags.class_symbol, $
                   ' not supported ***'
        END

        'mxSTRUCT_CLASS' : BEGIN
            print, '*** Formatting ', array_flags.class_symbol, $
                   ' not supported ***'
        END

        'mxOBJECT_CLASS' : BEGIN
            print, '*** Formatting ', array_flags.class_symbol, $
                   ' not supported ***'
        END

        'mxCHAR_CLASS'   : BEGIN
            data = string(_data)
        END

        'mxSPARSE_CLASS' : BEGIN
            print, '*** Formatting ', array_flags.class_symbol, $
                   ' not supported ***'
        END

        'mxDOUBLE_CLASS' : BEGIN
            IF array_flags.complex THEN $
                data = dcomplex(_data) $
            ELSE $
                data = double(_data)
        END

        'mxSINGLE_CLASS' : BEGIN
            IF array_flags.complex THEN $
                data = complex(_data) $
            ELSE $
                data = float(_data)
        END

        'mxINT8_CLASS'   : BEGIN
            print, '*** Formatting ', array_flags.class_symbol, $
                   ' not supported ***'
        END

        'mx_UINT8_CLASS' : BEGIN
            data = byte(_data)
        END

        'mxINT16_CLASS'  : BEGIN
            data = fix(_data, TYPE=2)
        END

        'mxUINT16_CLASS' : BEGIN
            data = uint(_data)
        END

        'mxINT32_CLASS'  : BEGIN
            data = long(_data)
        END

        'mxUINT32_CLASS' : BEGIN
            data = ulong(_data)
        END

    ENDCASE

    ;; Todo: transpose the data ?

    return, data

END

PRO load_mat, filename, path, data, STORE_LEVEL=store_level, $
              VERBOSE=verbose, DEBUG=debug

    header = { mat_v5_header, $
               description: "", $
               subsys_data_offset: 0ULL, $
               version: 0U, $
               endian_indicator: "" $
             }


    file = filepath(filename, ROOT_DIR=path)

    file_information = file_info(file)

    IF file_information.exists EQ 0 THEN BEGIN
        print, "File does not exist (", file, ")"
        return
    ENDIF

    IF file_information.directory EQ 1 THEN BEGIN
        print, "File is a directory (", file, ")"
        return
    ENDIF

    openr, lun, file, /GET_LUN

    ;; By default, create the variables on the $MAIN$ level

    IF NOT keyword_set(store_level) THEN store_level = 1

    IF keyword_set(debug) THEN BEGIN
        print
        print, '* Header'
    ENDIF

    ;; Todo: put this header code into a procedure

    description = bytarr(116)
    subsys_data_offset = 0ULL
    version = 0U
    endian_indicator = 0

    readu, lun, description
    readu, lun, subsys_data_offset
    readu, lun, version
    readu, lun, endian_indicator

    header.description = string(description)
    header.subsys_data_offset = subsys_data_offset
    header.version = version
    header.endian_indicator = $
        string(byte(ISHFT(endian_indicator AND 'FF00'XS, -8))) + $
        string(byte(endian_indicator AND '00FF'XS))

    IF keyword_set(DEBUG) THEN BEGIN
        print, 'Description           : ', header.description
        print, 'Subsystem data offset : ', $
               header.subsys_data_offset, FORMAT='(A,Z016)'
        print, 'Header version        : ', header.version, FORMAT='(A,Z04)'
        print, 'Endian                : ', header.endian_indicator
    ENDIF

    ;; Todo: must implement endian swapping

    data = 0
    data_element_number = 0

    WHILE NOT(eof(lun)) DO BEGIN

        IF keyword_set(debug) THEN BEGIN
            print, '=========================================================='
            print, '* Data Element ', data_element_number++
        ENDIF

        element_tag = element_tag_struct()
        read_element_tag, lun, element_tag, DEBUG=debug

        SWITCH element_tag.data_symbol OF

            'miMATRIX' : BEGIN

                ;; Array flags subelement

                IF keyword_set(debug) THEN BEGIN
                    print
                    print, '* Array flags subelement tag'
                ENDIF

                array_flags_tag = element_tag_struct()
                read_element_tag, lun, array_flags_tag, DEBUG=debug

                IF keyword_set(debug) THEN BEGIN
                    print, '* Array flags subelement data'
                ENDIF

                array_flags = subelement_array_flags_struct()
                read_subelement_array_flags, lun, array_flags_tag, array_flags, $
                                             DEBUG=debug

                ;; Dimensions array subelement

                IF keyword_set(debug) THEN BEGIN
                    print
                    print, '* Dimensions array subelement tag'
                ENDIF

                dimensions_array_tag = element_tag_struct()
                read_element_tag, lun, dimensions_array_tag, DEBUG=debug

                IF keyword_set(debug) THEN BEGIN
                    print, '* Dimensions array subelement data'
                ENDIF

                dimensions_array = subelement_dimensions_array_struct()
                read_subelement_dimensions_array, lun, dimensions_array_tag, $
                                                  dimensions_array, DEBUG=debug

                ;; Array name subelement

                IF keyword_set(debug) THEN BEGIN
                    print
                    print, '* Array name subelement tag'
                ENDIF

                array_name_tag = element_tag_struct()
                read_element_tag, lun, array_name_tag, DEBUG=debug

                IF keyword_set(debug) THEN BEGIN
                    print, '* Array name subelement data'
                ENDIF

                array_name = ''
                read_subelement_array_name, lun, array_name_tag, array_name, $
                                            DEBUG=debug

                IF keyword_set(verbose) THEN print, array_name

                ;; Real part (pr) subelement

                IF keyword_set(debug) THEN BEGIN
                    print
                    print, '* Real part (pr) subelement tag'
                ENDIF

                real_part_tag = element_tag_struct()
                read_element_tag, lun, real_part_tag, DEBUG=debug

                IF keyword_set(debug) THEN BEGIN
                    print, '* Real part (pr) subelement data'
                ENDIF

                read_element_data, lun, real_part_tag, real_data, DEBUG=debug

                data = real_data

                IF array_flags.complex THEN BEGIN

                    ;; Imaginary part (pi) subelement

                    IF keyword_set(debug) THEN BEGIN
                        print
                        print, '* Imaginary part (pi) subelement tag'
                    ENDIF

                    imag_part_tag = element_tag_struct()
                    read_element_tag, lun, imag_part_tag, DEBUG=debug

                    IF keyword_set(debug) THEN BEGIN
                        print, '* Imaginary part (pi) subelement data'
                    ENDIF

                    read_element_data, lun, imag_part_tag, imag_data, $
                                       DEBUG=debug

                    data = complex(real_data, imag_data)

                ENDIF

                data = format_array_element_data(data, array_flags, $
                                                 dimensions_array)
                
            END

        ENDSWITCH

        ;; Create a variable on the main level using the undocumented
        ;; IDL routine ROUTINE_NAMES. This only works for IDL 5.3 and
        ;; higher.

        foo = routine_names(array_name, data, STORE=store_level)

        IF keyword_set(debug) THEN BEGIN
            point_lun, -lun, current_file_position
            print, 'Current file position : ', current_file_position, $
                   FORMAT='(A, Z08)'
        ENDIF

    ENDWHILE

    close, lun
    free_lun, lun

END

;@GJ, 2024/7/21, processing the experimental data from Zhongwei Bian
PRO STED_MPI_gaussian_donut_BUAA_phantom

  ;  filename='test0.mat'
  ;  filename='test1.mat'
  ;  filename='test2.mat'
  ;  filename='test10.mat'
  ;  filename='test20.mat'
;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\2_10mT_xz_matdata\'
;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\2_15mT_xz_matdata\'
;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\x_axis_cali\'
;  path = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\', /MUST_EXIST, TITLE="mat files", /DIRECTORY)
  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\BUAA_PSF\'
;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\20241128\1_2_NWU\'
;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\20241128\pr3_perimag\'
  IF N_ELEMENTS(path) EQ 0 THEN RETURN
  IF STRLEN(path) EQ 0 THEN RETURN
  read_BZW_mat_file, path, n_lines, n_pxls, abs_3rd, real_3rd, imag_3rd, phase_3rd, abs_2nd, real_2nd, imag_2nd, phase_2nd, abs_3rd_fov, phase_3rd_fov, real_3rd_fov, imag_3rd_fov

  ;@GJ, 2024/12/2, do the angle calculation
  ;scatter plot the phase
  iplot, real_3rd, imag_3rd, xtitle='Real', ytitle='Imag', color='blue', SYM_INDEX=1, LINESTYLE=6, title='PSF Scatter Plot', /NO_SAVEPROMPT
  real_1d = REFORM(real_3rd, n_lines*n_pxls)
  imag_1d = REFORM(imag_3rd, n_lines*n_pxls)
  result = POLY_FIT(real_1d, imag_1d, 6, yfit=imag_1d_fit)
  iplot, real_1d, imag_1d_fit, COLOR='blue', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
  
  ;rotate to get the donut-shaped focal spot
  max_real = MAX(ABS(real_1d), maxid)
  angle = 180./!PI*ATAN(imag_1d_fit[maxid], real_1d[maxid], /phase)

  ;@GJ, 2024/8/23, plot the rotated
  phi_g = angle / 180. * !PI - 0.5 * !PI
  new_3rd_fov_g = -real_3rd_fov * sin(phi_g) + imag_3rd_fov * cos(phi_g)
  new_3rd_fov_imag_g = -real_3rd_fov * sin(phi_g) + imag_3rd_fov * cos(phi_g)
  new_3rd_fov_real_g = imag_3rd_fov * sin(phi_g) + real_3rd_fov * cos(phi_g)
  iplot, new_3rd_fov_real_g, new_3rd_fov_imag_g, COLOR='green', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  real_1d_g = imag_1d_fit * sin(phi_g) + real_1d * cos(phi_g)
  imag_1d_g = -real_1d * sin(phi_g) + imag_1d_fit * cos(phi_g)
  iplot, real_1d_g, imag_1d_g, COLOR='green', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
  max_real_g = MAX(ABS(real_1d_g), max_real_id_g)
  ;  max_imag_g = MAX(ABS(imag_1d_g), max_imag_id_g)
  
  phi_d = angle / 180. * !PI
  new_3rd_fov_d = -real_3rd_fov * sin(phi_d) + imag_3rd_fov * cos(phi_d)
  new_3rd_fov_imag_d = -real_3rd_fov * sin(phi_d) + imag_3rd_fov * cos(phi_d)
  new_3rd_fov_real_d = imag_3rd_fov * sin(phi_d) + real_3rd_fov * cos(phi_d)
  iplot, new_3rd_fov_real_d, new_3rd_fov_imag_d, COLOR='red', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  real_1d_d = imag_1d_fit * sin(phi_d) + real_1d * cos(phi_d)
  imag_1d_d = -real_1d * sin(phi_d) + imag_1d_fit * cos(phi_d)
  iplot, real_1d_d, imag_1d_d, COLOR='red', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
;  max_real_d = MAX(ABS(real_1d_d), max_real_id_d)
  max_imag_d = MAX(ABS(imag_1d_d), max_imag_id_d)
  
  ;calculate dog factor
  dog_factor = ABS(imag_1d_g[max_real_id_g]) / max_imag_d

  ;scatter plot the data
  iplot, new_3rd_fov_real_d, new_3rd_fov_imag_d*dog_factor, COLOR='pink', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  iplot, real_1d_d, imag_1d_d*dog_factor, COLOR='pink', SYM_INDEX=0, THICK=2, LINESTYLE=0, /OVERPLOT
;  iplot, real_1d_d, imag_1d_g - imag_1d_d*dog_factor, xtitle='Real', ytitle='Imag', color='red', title='STED-MPI'
  
  ;@GJ, 2025/4/3, CEO, curve of excitation offset
  sort_ind = SORT(real_1d_d)
  iplot, new_3rd_fov_real_d, new_3rd_fov_imag_d, xrange=[-MAX(abs(real_1d_d))*0.2, MAX(abs(real_1d_d))*1.1], yrange=[-MAX(abs(real_1d_d))*0.2, MAX(abs(real_1d_d))*1.1], COLOR='blue', SYM_INDEX=1, LINESTYLE=6, xtitle='G3R', ytitle='G3I', title='CEO (Curve of Excitation Offset)', /NO_SAVEPROMPT
  iplot, real_1d_d[sort_ind], imag_1d_d[sort_ind], COLOR='red', SYM_INDEX=0, THICK=3, LINESTYLE=0, /OVERPLOT
  
  ;calculate the STED image
  new_3rd_fov_imag_STED = new_3rd_fov_imag_g - new_3rd_fov_imag_d*dog_factor
  
  ;@GJ, 2024/12/5, plot the PSFs of Gaussian, donut, and STED
  max_center = MAX(new_3rd_fov_imag_STED, max_index)
  center_2d = ARRAY_INDICES(new_3rd_fov_imag_STED, max_index)
  offset_field_fov_imag = new_3rd_fov_imag_STED * 0.
  gradient_x = 1.7; mT/mm
  gradient_y = 1.7; mT/mm
  fov = 30; mm
  delta_offset_field_x = gradient_x * fov / n_pxls
  delta_offset_field_y = gradient_y * fov / n_pxls
  FOR i=0, n_pxls-1 DO BEGIN
    FOR j=0, n_pxls-1 DO BEGIN
      IF i NE center_2d[0] AND SIGNUM(j-center_2d[1]) THEN BEGIN
        offset_field_fov_imag[i, j] = SIGNUM(i-center_2d[0])*SIGNUM(j-center_2d[1])*SQRT((i-center_2d[0])^2*delta_offset_field_x^2 + (j-center_2d[1])^2*delta_offset_field_y^2)
      ENDIF ELSE BEGIN
        offset_field_fov_imag[i, j] = SQRT((i-center_2d[0])^2*delta_offset_field_x^2 + (j-center_2d[1])^2*delta_offset_field_y^2)
      ENDELSE
    ENDFOR
  ENDFOR
  iplot, offset_field_fov_imag, new_3rd_fov_imag_g, xtitle='Offset Field [mT]', ytitle='I3', color='red', SYM_INDEX=3, LINESTYLE=6, title='STED-MPI PSFs', /NO_SAVEPROMPT
  iplot, offset_field_fov_imag, new_3rd_fov_imag_d * dog_factor, color='green', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  iplot, offset_field_fov_imag, new_3rd_fov_imag_STED, color='blue', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  offset_field_fov_imag_1d = REFORM(offset_field_fov_imag, n_pxls*n_pxls)
  new_3rd_fov_imag_g_1d = REFORM(new_3rd_fov_imag_g, n_pxls*n_pxls)
  new_3rd_fov_imag_sted_1d = REFORM(new_3rd_fov_imag_sted, n_pxls*n_pxls)
  sort_offset = SORT(offset_field_fov_imag_1d)
  offset_field_fov_imag_1d = offset_field_fov_imag_1d[sort_offset]
  new_3rd_fov_imag_g_1d = new_3rd_fov_imag_g_1d[sort_offset]
  new_3rd_fov_imag_sted_1d = new_3rd_fov_imag_sted_1d[sort_offset]
  yfit_psf_g = GAUSSFIT(offset_field_fov_imag_1d, new_3rd_fov_imag_g_1d, coeff_g, NTERMS=4)
;  print, 'Gaussian PSF Fit Result: ', coeff_g[0:4-1]
  print, 'Gaussian FWHM: ', 2*SQRT(2*ALOG(2))*coeff_g[2], ' mT, ', 2*SQRT(2*ALOG(2))*coeff_g[2] / gradient_x, ' mm'
  iplot, offset_field_fov_imag_1d, yfit_psf_g, LINSTYLE=0, thick=4, color='red', /OVERPLOT
  yfit_psf_STED = GAUSSFIT(offset_field_fov_imag_1d, new_3rd_fov_imag_sted_1d, coeff_STED, NTERMS=4)
;  print, 'STED PSF Fit Result: ', coeff_STED[0:4-1]
  print, 'STED FWHM: ', 2*SQRT(2*ALOG(2))*coeff_STED[2], ' mT, ', 2*SQRT(2*ALOG(2))*coeff_STED[2] / gradient_x, ' mm'
  iplot, offset_field_fov_imag_1d, yfit_psf_sted, LINSTYLE=0, thick=4, color='blue', /OVERPLOT
  yfit_psf_d = yfit_psf_g - yfit_psf_STED
  iplot, offset_field_fov_imag_1d, yfit_psf_d, LINSTYLE=0, thick=4, color='green', /OVERPLOT
  
  ;@GJ, 2025/3/25, evaluate the Radon 1d RL filter
  NRHO = 132.;parameters.N_points; 237.;137;;137
  NTHETA = 90.;parameters.N_angles;27;1. * CEIL(!PI * NRHO/2.);127.;181.;27.;181
  new_3rd_fov_psf_sted = new_3rd_fov_imag_STED
  sinogram_STED = RADON(new_3rd_fov_psf_sted, RHO=rho, THETA=theta, NRHO=NRHO, NTHETA=NTHETA)
;  filtered_resultSTED_decon = sinogram_STED * 0.
;  Filter_type = 'Ramp'
;  FOR i=0, NTHETA-1 DO BEGIN
;    blurred = reform(sinogram_STED[i, *])
;    restored = blurred;rl_deconv1d(blurred, REFORM(new_3rd_fov_imag_STED[*, center_2d[1]]), iterations=100)
;    filtered_resultSTED_decon[i, *] = Filter_FFT(restored, filter_type)
;  ENDFOR
;  backproject_STED_decon = RADON(filtered_resultSTED_decon, /BACKPROJECT, RHO=rho, THETA=theta)
  ;display the images
  
  iimage, BYTSCL(new_3rd_fov_imag_g), view_title='Gaussian', RGB_TABLE=3, VIEW_GRID=[3,1], DIMENSIONS=[n_pxls*3., n_pxls], WINDOW_TITLE='STED-PSF', /NO_SAVEPROMPT
  iimage, BYTSCL(new_3rd_fov_imag_d), view_title='Donut', RGB_TABLE=8, /VIEW_NEXT
  iimage, BYTSCL(new_3rd_fov_imag_STED), view_title='STED', RGB_TABLE=1, /VIEW_NEXT
;  iimage, BYTSCL(backproject_STED_decon), view_title='STED Decon', /VIEW_NEXT
 
  ;@GJ, 2025/1/13, STED filter kernel
;  STED_filter_FFT = FFT(new_3rd_fov_imag_STED)
;  iplot, ABS(STED_filter_FFT), title='STED kernel'
  ;@GJ, 2025/1/13
  ;evaluate the linearity
  
;  ;@GJ, do gaussian curve fitting
;  x_offset_field = (FINDGEN(n_pxls) - center_2d[0]) *  delta_offset_field_x
;  x_psf_g = REFORM(new_3rd_fov_imag_g[*, center_2d[1]])
;  x_psf_d = REFORM(new_3rd_fov_imag_d[*, center_2d[1]])
;  x_psf_STED = REFORM(new_3rd_fov_imag_STED[*, center_2d[1]])
;  iplot,x_offset_field, x_psf_g, xtitle='Offset Field [mT]', ytitle='PSF', color='green', SYM_INDEX=3, LINESTYLE=0, title='PSF', /NO_SAVEPROMPT
;  iplot, x_offset_field, x_psf_d/MAX(x_psf_d)*MAX(x_psf_g), color='red', SYM_INDEX=3, LINESTYLE=0, /OVERPLOT
;  iplot, x_offset_field, x_psf_STED/MAX(x_psf_STED)*MAX(x_psf_g), color='blue', SYM_INDEX=3, LINESTYLE=0, /OVERPLOT
    
  psf_sted_filename = path+'psf_sted.tif'
  write_tiff, psf_sted_filename, BYTSCL(new_3rd_fov_imag_g - new_3rd_fov_imag_d*dog_factor), /long
  psf_sted_filename = path+'psf_sted.png'
  write_png, psf_sted_filename, BYTSCL(new_3rd_fov_imag_g - new_3rd_fov_imag_d*dog_factor)
  
  path = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\', /MUST_EXIST, TITLE="mat files", /DIRECTORY)
  ;path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\BUAA_B2\'
  IF N_ELEMENTS(path) EQ 0 THEN RETURN
  IF STRLEN(path) EQ 0 THEN RETURN
  read_BZW_mat_file, path, n_lines, n_pxls, abs_3rd, real_3rd, imag_3rd, phase_3rd, abs_2nd, real_2nd, imag_2nd, phase_2nd, abs_3rd_fov, phase_3rd_fov, real_3rd_fov, imag_3rd_fov
  iplot, real_3rd, imag_3rd, xtitle='Real', ytitle='Imag', color='blue', SYM_INDEX=1, LINESTYLE=6, title='Image Scatter Plot', /NO_SAVEPROMPT

  ;@GJ, calculate the border of the LEO
  half_range_index = WHERE(abs_3rd_fov GT 0.5*MAX(abs_3rd_fov), half_count)
  min_phase = MIN(phase_3rd_fov[half_range_index])
  max_phase = MAX(phase_3rd_fov[half_range_index])
  angle_image = min_phase
  
  ;Rotate to gaussian
  phi_g = angle_image / 180. * !PI - 0.5 * !PI
  new_3rd_fov_g = -real_3rd_fov * sin(phi_g) + imag_3rd_fov * cos(phi_g)
  new_3rd_fov_imag_g = -real_3rd_fov * sin(phi_g) + imag_3rd_fov * cos(phi_g)
  new_3rd_fov_real_g = imag_3rd_fov * sin(phi_g) + real_3rd_fov * cos(phi_g)
  iplot, new_3rd_fov_real_g, new_3rd_fov_imag_g, COLOR='green', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  
  ;Rotate the donut
  phi_d = angle_image / 180. * !PI
  new_3rd_fov_d = -real_3rd_fov * sin(phi_d) + imag_3rd_fov * cos(phi_d)
  new_3rd_fov_imag_d = -real_3rd_fov * sin(phi_d) + imag_3rd_fov * cos(phi_d)
  new_3rd_fov_real_d = imag_3rd_fov * sin(phi_d) + real_3rd_fov * cos(phi_d)
  iplot, new_3rd_fov_real_d, new_3rd_fov_imag_d, COLOR='red', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  iplot, new_3rd_fov_real_d, new_3rd_fov_imag_d*dog_factor, COLOR='pink', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT
  
  ;@GJ, 2025/4/3, CEO, curve of excitation offset
;  sort_ind = SORT(real_1d_d)
  iplot, new_3rd_fov_real_d, new_3rd_fov_imag_d, xrange=[-MAX(abs(new_3rd_fov_real_d))*0.2, MAX(abs(new_3rd_fov_real_d))*1.1], yrange=[-MAX(abs(new_3rd_fov_real_d))*0.2, MAX(abs(new_3rd_fov_real_d))*1.1], COLOR='blue', SYM_INDEX=1, LINESTYLE=6, xtitle='G3R', ytitle='G3I', title='CEO (Curve of Excitation Offset)', /NO_SAVEPROMPT
;  iplot, real_1d_d[sort_ind], imag_1d_d[sort_ind], COLOR='red', SYM_INDEX=0, THICK=3, LINESTYLE=0, /OVERPLOT
  
  ;@GJ, 2024/12/4, do STED subtraction
  new_3rd_fov_imag_sted = new_3rd_fov_imag_g - new_3rd_fov_imag_d*dog_factor
  ;@GJ, 2025/3/25, evaluate the Radon 1d RL filter
  sinogram_fov_imag_STED = RADON(new_3rd_fov_imag_STED, RHO=rho, THETA=theta, NRHO=NRHO, NTHETA=NTHETA)
  filtered_resultSTED_decon_fov_imag = sinogram_fov_imag_STED * 0.
  Filter_type = 'Ramp'
  PSF_kernel = MEAN(sinogram_STED, dimension=1, /double)
  x_fit = FINDGEN(NRHO)-center_2d[0]/n_pxls*NRHO
  ;yfit_psf_kernel = GAUSSFIT(x_fit, PSF_kernel, coeff_g, NTERMS=4)
  yfit_psf_kernel = PSF_kernel
  yfit_psf_kernel = yfit_psf_kernel - MIN(yfit_psf_kernel)
  ;@GJ, 2025/2/3, do the deconvolution of G3STED
;  temp_psf_FFT = FFT(yfit_psf_kernel)
;  k = MAX(ABS(temp_psf_FFT))/2.
;  wiener_filter = CONJ(temp_psf_FFT) / (ABS(temp_psf_FFT)^2 + k)
;  thre_decon = MAX(ABS(temp_psf_FFT))/5.
;  psf_FFT = temp_psf_FFT * (ABS(temp_psf_FFT) GT thre_decon) + (ABS(temp_psf_FFT) LT thre_decon)

;  iplot, PSF_kernel, title='PSF_kernel', /NO_SAVEPROMPT
  FOR i=0, NTHETA-1 DO BEGIN
    image_1d = reform(sinogram_fov_imag_STED[i, *])
    N = N_ELEMENTS(image_1d)
    T = 0.1
    X = (FINDGEN((N - 1)/2) + 1)
    is_N_even = (N MOD 2) EQ 0
    if (is_N_even) then $
      freqs = [0.0, X, N/2, -N/2 + X]/(N*T) $
    else $
      freqs = [0.0, X, -(N/2 + 1) + X]/(N*T)
      
    P_FFT = FFT(image_1d)
    filtered_P_FFT = P_FFT
    filter = abs(freqs) / abs(FFT(yfit_psf_kernel)); * cos(!PI * freqs / (N - 1))
    filter = filter * (abs(FFT(yfit_psf_kernel)) GE 0.1*MAX(abs(FFT(yfit_psf_kernel))))
;    filter = abs(freqs) / psf_FFT
    filtered_P_FFT[1:*] = P_FFT[1:*] * filter[1:*]
    abs_filtered_P_FFT = (FFT(filtered_P_FFT, /INVERSE))
    filtered_resultSTED_decon_fov_imag[i, *] = abs_filtered_P_FFT
    IF i EQ 0 THEN BEGIN
      iplot, image_1d, /NO_SAVEPROMPT
      iplot, abs_filtered_P_FFT/max(abs_filtered_P_FFT)*max(image_1d), color='red', /overplot
;      iplot, yfit_psf_kernel/MAX(ABS(yfit_psf_kernel))*MAX(ABS(sinogram_fov_imag_STED)), color='red', /overplot
;;      iplot, rl_deconv1d(blurred_fov_imag, PSF_kernel_temp, iterations=100), color='blue', /overplot
;      iplot, convol(restored_fov_imag, yfit_psf_kernel, /center, /edge_zero, /normalize), color='green', /overplot
    ENDIF
  ENDFOR
  backproject_STED_decon_fov_imag = RADON(filtered_resultSTED_decon_fov_imag, /BACKPROJECT, RHO=rho, THETA=theta)
  fbp_size = size(backproject_STED_decon_fov_imag, /dim)
  STED_decon_fov_imag = backproject_STED_decon_fov_imag[fbp_size[0]/2-n_pxls/2:fbp_size[0]/2+n_pxls/2-1, fbp_size[1]/2-n_pxls/2:fbp_size[1]/2+n_pxls/2-1]
;  iplot, PSF_kernel, title='PSF_kernel replot'
  ; Filter the degraded image with the Wiener filter
;  powerClean = ABS(FFT(new_3rd_fov_imag_sted, /CENTER))^2
;  powerNoise = ABS(FFT(new_3rd_fov_psf_sted, /CENTER))^2
;  imageFiltered = WIENER_FILTER(new_3rd_fov_imag_sted, new_3rd_fov_psf_sted, powerClean, powerNoise)
  iimage, BYTSCL(new_3rd_fov_imag_g), view_title='Gaussian', RGB_TABLE=3, VIEW_GRID=[4,1], DIMENSIONS=[n_pxls*4., n_pxls], WINDOW_TITLE='STED-Image', /NO_SAVEPROMPT
  iimage, BYTSCL(new_3rd_fov_imag_d), view_title='Donut', RGB_TABLE=8, /VIEW_NEXT
  iimage, BYTSCL(new_3rd_fov_imag_g - new_3rd_fov_imag_d*dog_factor), RGB_TABLE=1, view_title='STED', /VIEW_NEXT
  iimage, BYTSCL(STED_decon_fov_imag), view_title='STED Decon', RGB_TABLE=64, /VIEW_NEXT

  image_gaussian_filename = path+'image_gaussian.tif'
  write_tiff, image_gaussian_filename, BYTSCL(new_3rd_fov_imag_g), /long
  image_donut_filename = path+'image_donut.tif'
  write_tiff, image_donut_filename, BYTSCL(new_3rd_fov_imag_d), /long
  image_sted_filename = path+'image_sted.tif'
  write_tiff, image_sted_filename, BYTSCL(new_3rd_fov_imag_sted), /long
  image_sted_filename = path+'image_sted_decon.tif'
  write_tiff, image_sted_filename, BYTSCL(STED_decon_fov_imag), /long
  
  image_gaussian_filename = path+'image_gaussian.png'
  write_png, image_gaussian_filename, BYTSCL(new_3rd_fov_imag_g)
  image_donut_filename = path+'image_donut.png'
  write_png, image_donut_filename, BYTSCL(new_3rd_fov_imag_d)
  image_sted_filename = path+'image_sted.png'
  write_png, image_sted_filename, BYTSCL(new_3rd_fov_imag_sted)
  image_sted_filename = path+'image_sted_decon.png'
  write_png, image_sted_filename, BYTSCL(STED_decon_fov_imag)
END

;@GJ, 2024/12/4, display the reconstructed by RL
PRO STED_MPI_Richardson_Lucy_Display
  path = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\', /MUST_EXIST, TITLE="mat files", /DIRECTORY)
  image_gaussian_filename = path+'image_gaussian.tif'
  image_gaussian = read_tiff(image_gaussian_filename)
  image_donut_filename = path+'image_donut.tif'
  image_donut = read_tiff(image_donut_filename)
  image_sted_filename = path+'image_sted.tif'
  image_sted = read_tiff(image_sted_filename)
  image_RL_filename = path+'image_RL.tif'
  image_RL = read_tiff(image_RL_filename)

  iimage, BYTSCL(image_gaussian), view_title='Gaussian', VIEW_GRID=[4,1], DIMENSIONS=[(size(image_RL, /dim))[0]*4., (size(image_RL, /dim))[1]], WINDOW_TITLE='STED-Decon', /NO_SAVEPROMPT
  iimage, BYTSCL(image_donut), view_title='Donut', /VIEW_NEXT
  iimage, BYTSCL(image_sted), view_title='STED', /VIEW_NEXT
  iimage, BYTSCL(image_RL), view_title='Deconvolution', /VIEW_NEXT

END

;@GJ, 2024/12/3, load the file
PRO read_BZW_mat_file, path, n_lines, n_pxls, abs_3rd, real_3rd, imag_3rd, phase_3rd, abs_2nd, real_2nd, imag_2nd, phase_2nd, abs_3rd_fov, phase_3rd_fov, real_3rd_fov, imag_3rd_fov
   
  IF N_ELEMENTS(path) EQ 0 THEN RETURN
  IF STRLEN(path) EQ 0 THEN RETURN

  file_array = file_search(path, 'test*.mat', count=n_lines_temp)
  n_lines = n_lines_temp; 31.
  n_pxls = 50. ;number of pixels

  abs_3rd = DBLARR(n_lines, n_pxls) * 0.
  real_3rd = DBLARR(n_lines, n_pxls) * 0.
  imag_3rd = DBLARR(n_lines, n_pxls) * 0.
  phase_3rd = DBLARR(n_lines, n_pxls) * 0.

  abs_2nd = DBLARR(n_lines, n_pxls) * 0.
  real_2nd = DBLARR(n_lines, n_pxls) * 0.
  imag_2nd = DBLARR(n_lines, n_pxls) * 0.
  phase_2nd = DBLARR(n_lines, n_pxls) * 0.

  ;save the reslults
  n_harmonics = 10
  abs_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  real_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  imag_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  phase_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  amp_phase_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.

  n_3rd_mp = 11
  abs_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  real_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  imag_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  phase_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  amp_phase_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.

  ;@GJ, 2024/7/22, fov is 30 mm
  fov = 30. ;mm
  N = 50000. ;total time elements
  T = 0.4 ;us
  t_array = FINDGEN(N) * T ;us
  x_array = FINDGEN(n_pxls)/(n_pxls-1)*fov - fov/2. ; fov = 30 mm
  ; N is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X = (FINDGEN((N - 1)/2) + 1)
  is_N_even = (N MOD 2) EQ 0
  if (is_N_even) then $
    freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
  else $
    freq = [0.0, X, -(N/2 + 1) + X]/(N*T)

  freq_base = 0.025 ;MHz
  ;@GJ, 2024/12/6, standardize the starting point of the excitation
  N_period = 1. / freq_base / T
  ;calculate the shift elements
  filename = 'test0.mat'
  load_mat, filename, path, signal_data, STORE_LEVEL=1
  temp_signal = REFORM(signal_data[0:N-1, n_pxls/2])
  one_period_signal = temp_signal[0:N_period+N_period/10]
  max_half = MAX(one_period_signal, max_index)
  print, 'shift: ', max_index
  temp_signal = SHIFT(temp_signal, -max_index)
;  iplot, temp_signal[0:1000], title='Shifted Time Signal'

  FOR k=0, n_lines-1 DO BEGIN
    filename = 'test'+STRING(k, format='(I0)')+'.mat'
    load_mat, filename, path, signal_data, STORE_LEVEL=1

    line_abs = DBLARR(n_pxls, N) * 0.
    line_real = DBLARR(n_pxls, N) * 0.
    line_imag = DBLARR(n_pxls, N) * 0.
    line_phase = DBLARR(n_pxls, N) * 0.

    FOR i=0, n_pxls-1 DO BEGIN
      temp_signal = SHIFT(REFORM(signal_data[0:N-1, i]), -max_index)
      temp_signal_fft = FFT(temp_signal)
      line_abs[i, *] = ABS(temp_signal_fft)
      line_real[i, *] = Real_part(temp_signal_fft)
      line_imag[i, *] = Imaginary(temp_signal_fft)
;      line_imag[i, *] = -Imaginary(temp_signal_fft)
      ;@GJ, 2024/1/5, correcting the phase calculation
      FOR j=0, N-1 DO BEGIN
        IF line_abs[i, j] GT 1e-6 THEN BEGIN
          line_phase[i, j] = 180./!PI*ATAN(temp_signal_fft[j], /phase)
          IF line_phase[i, j] LT 0 THEN line_phase[i, j] += 360.
          ;          line_phase[i, j] = 180./!PI*ATAN(Imaginary(temp_signal_fft[j]), Real_part(temp_signal_fft[j]))
          ;          IF line_phase[i, j] GT 180. THEN line_phase[i, j] = 360. - line_phase[i, j]
          ;          IF line_phase[i, j] LT 0. THEN line_phase[i, j] = 360. + line_phase[i, j]
        ENDIF
      ENDFOR
      IF k EQ FLOOR(n_lines/2) THEN BEGIN
        IF i EQ FLOOR(n_pxls/2) THEN BEGIN
;          iplot, freq/freq_base, temp_signal_fft, xrange=[0., 10.], color='red', xtitle='Harmonics', title='Signal FFT', /NO_SAVEPROMPT ELSE iplot, freq/freq_base, temp_signal_fft, color=350000*(i+1), /overplot
          ind_1st = freq_base*1./freq[1]
          ind_1st_neg = N - ind_1st
          temp_signal_fft[ind_1st] *= 0. ;temp_signal_fft[3.*ind_1st]
          temp_signal_fft[ind_1st_neg] *= 0.; temp_signal_fft[3.*ind_1st]
          filtered_signal = FFT(temp_signal_fft, /inverse)
          temp_signal_period = DBLARR(N_period) * 0.
          FOR m=0, N/N_period-1 DO temp_signal_period += temp_signal[m*N_period:(m+1)*N_period-1]
          filtered_signal_period = DBLARR(N_period) * 0.
          FOR m=0, N/N_period-1 DO filtered_signal_period += filtered_signal[m*N_period:(m+1)*N_period-1]
;          iplot, temp_signal_period/MAX(temp_signal_period)
;          iplot, filtered_signal_period/MAX(filtered_signal_period), color='red', /overplot
;          iplot, temp_signal[0:1000]
;          iplot, filtered_signal[0:1000], color='red', /overplot
          
          ;process in a different
          temp_signal_period_fft = FFT(temp_signal_period)
          temp_signal_period_fft[1] *= 0.
          temp_signal_period_fft[99] *= 0.
          filtered_signal_period = FFT(temp_signal_period_fft, /inverse)
;          iplot, temp_signal_period/MAX(temp_signal_period)
;          iplot, filtered_signal_period/MAX(filtered_signal_period), color='red', /overplot

        ENDIF
        ;    IF i EQ 0 THEN iplot, freq/freq_base, REFORM(line_phase[i, *]), xrange=[0., 10.], xtitle='Harmonics', title='Phase', /NO_SAVEPROMPT ELSE iplot, freq/freq_base, REFORM(line_phase[i, *]), /overplot
      ENDIF
    ENDFOR

    ind_3rd = freq_base*3./freq[1]
    abs_3rd[k, *] = REFORM(line_abs[*, ind_3rd])
    real_3rd[k, *] = REFORM(line_real[*, ind_3rd])
    imag_3rd[k, *] = REFORM(line_imag[*, ind_3rd])
    phase_3rd[k, *] = REFORM(line_phase[*, ind_3rd])

    ind_2nd = freq_base*2./freq[1]
    abs_2nd[k, *] = REFORM(line_abs[*, ind_2nd])
    real_2nd[k, *] = REFORM(line_real[*, ind_2nd])
    imag_2nd[k, *] = REFORM(line_imag[*, ind_2nd])
    phase_2nd[k, *] = REFORM(line_phase[*, ind_2nd])

    ;@GJ, 2024/8/11, calculate the harmonics
    FOR l=0, n_harmonics-1 DO BEGIN
      ind_har = freq_base*(l+1.)/(freq[1]-freq[0])
      abs_array[l, k, *] = REFORM(line_abs[*, ind_har])
      real_array[l, k, *] = REFORM(line_real[*, ind_har])
      imag_array[l, k, *] = REFORM(line_imag[*, ind_har])
      phase_array[l, k, *] = REFORM(line_phase[*, ind_har])
      amp_phase_array[l, k, *] = REFORM(line_abs[*, ind_har]) * REFORM(line_phase[*, ind_har])
    ENDFOR

    ;@GJ, 2024/8/11, calculate the 2.5-3.5 Harmonics
    FOR m=0, n_3rd_mp-1 DO BEGIN
      ind_har = freq_base*(2.5+m/10.)/(freq[1]-freq[0])
      abs_3rd_mp_array[m, k, *] = REFORM(line_abs[*, ind_har])
      real_3rd_mp_array[m, k, *] = REFORM(line_real[*, ind_har])
      imag_3rd_mp_array[m, k, *] = REFORM(line_imag[*, ind_har])
      phase_3rd_mp_array[m, k, *] = REFORM(line_phase[*, ind_har])
      amp_phase_3rd_mp_array[m, k, *] = REFORM(line_abs[*, ind_har]) * REFORM(line_phase[*, ind_har])
    ENDFOR
  ENDFOR
  
  ;Fingerprinting
;  max_id = MAX(REFORM(abs_array[1, *, *]), max_index)
;  max_ind_2d = ARRAY_INDICES(REFORM(abs_array[1, *, *]), max_index)
;  iplot, REFORM(real_array[1, max_ind_2d[0], *]), REFORM(imag_array[1, max_ind_2d[0], *])
  color_array = ['red', 'blue', 'green', 'purple', 'peru', 'orchid', 'magenta', 'violet', 'navy', 'teal', 'sienna']
  max_abs = MAX(abs_array[1:*, *, *])
  flag = 0;
  FOR i=1, 9 DO BEGIN
    IF MAX(abs_array[i, *, *])/max_abs GT 0.02 THEN BEGIN
      ;@GJ, calculate the border of the LEO
      half_range_index = WHERE(REFORM(abs_array[i, *, *]) GT 0.5*MAX(abs_array[i, *, *]), half_count)
      min_phase = MIN((REFORM(abs_array[i, *, *]))[half_range_index])
      max_phase = MAX((REFORM(abs_array[i, *, *]))[half_range_index])
      angle_image = min_phase
      
      ;Rotate to gaussian
      phi_d = angle_image / 180. * !PI
      new_imag_d = -REFORM(real_array[i, *, *]) * sin(phi_d) + REFORM(imag_array[i, *, *]) * cos(phi_d)
      new_real_d = REFORM(imag_array[i, *, *]) * sin(phi_d) + REFORM(real_array[i, *, *]) * cos(phi_d)

      sign = (-1)^(i)
      IF flag EQ 0 THEN BEGIN
;        iplot, sign*new_real_d/MAX(abs_array[i, *, *])*max_abs, sign*new_imag_d/MAX(abs_array[i, *, *])*max_abs, color=color_array[i-1], SYM_INDEX=3, THICK=2, LINESTYLE=6, xtitle='Real', ytitle='Imag', title='Fingerprint', /NO_SAVEPROMPT
      ENDIF ELSE BEGIN
;        iplot, sign*new_real_d/MAX(abs_array[i, *, *])*max_abs, sign*new_imag_d/MAX(abs_array[i, *, *])*max_abs, color=color_array[i-1], SYM_INDEX=3, THICK=2, LINESTYLE=6, /overplot
      ENDELSE
      flag = 1
    ENDIF
  ENDFOR

  abs_3rd_fov = CONGRID(abs_3rd, n_pxls, n_pxls)
  phase_3rd_fov = CONGRID(phase_3rd, n_pxls, n_pxls)
  real_3rd_fov = CONGRID(real_3rd, n_pxls, n_pxls)
  imag_3rd_fov = CONGRID(imag_3rd, n_pxls, n_pxls)
;  iimage, abs_3rd_fov, view_title='3rd Amp', VIEW_GRID=[5,7], DIMENSIONS=[n_pxls*6., n_pxls*8.], WINDOW_TITLE='3rd Harmonic', /NO_SAVEPROMPT
;  iimage, real_3rd_fov, view_title='3rd Real', /VIEW_NEXT
;  iimage, imag_3rd_fov, view_title='3rd Imag', /VIEW_NEXT
;  iimage, BYTSCL(phase_3rd_fov), view_title='3rd Phase', /VIEW_NEXT
;  iimage, 255.-BYTSCL(phase_3rd_fov), view_title='3rd -Phase', /VIEW_NEXT
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (100.+i*10.)/180.
    new_3rd_fov = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    new_3rd_fov_imag = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    new_3rd_fov_real = imag_3rd_fov * sin(phi) + real_3rd_fov * cos(phi)
;    iimage, BYTSCL(new_3rd_fov), view_title='3rd New', /VIEW_NEXT
  ENDFOR
  
  ;@GJ, 2024/7/24, complete image
  gap = 3.
  whole_harmonics = DBLARR(n_harmonics*(n_pxls+gap), 5.*(n_pxls+gap)) * 0.
  FOR i=0, n_harmonics-1 DO BEGIN
    ;@GJ, 2024/8/7, phase is minus
    temp_phase = CONGRID(REFORM(-phase_array[i, *, *]), n_pxls, n_pxls)
    max_p = MAX(temp_phase[n_pxls/3.:n_pxls*2./3., n_pxls/3.:n_pxls*2./3.], min=min_p)
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 0:n_pxls-1] = BYTSCL(temp_phase, max=max_p, min=min_p)
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, n_pxls+gap:n_pxls+gap+n_pxls-1] = BYTSCL(CONGRID(REFORM(imag_array[i, *, *]), n_pxls, n_pxls))
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 2.*(n_pxls+gap):2.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(real_array[i, *, *]), n_pxls, n_pxls))
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 3.*(n_pxls+gap):3.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(abs_array[i, *, *]), n_pxls, n_pxls))
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 4.*(n_pxls+gap):4.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(amp_phase_array[i, *, *]), n_pxls, n_pxls))
  ENDFOR

  ;@GJ, 2024/7/24, save the results
  whole_harmonics_file = path + 'complete_harmonics.jpg'
  WRITE_JPEG, whole_harmonics_file, whole_harmonics

  ;@GJ, 2024/8/11, 2.5 to 3.5th harmonic
  whole_3rd_mp_harmonics = DBLARR(n_3rd_mp*(n_pxls+gap), 5.*(n_pxls+gap)) * 0.
  FOR i=0, n_3rd_mp-1 DO BEGIN
    ;@GJ, 2024/8/7, phase is minus
    temp_phase = CONGRID(REFORM(-phase_3rd_mp_array[i, *, *]), n_pxls, n_pxls)
    max_p = MAX(temp_phase[n_pxls/3.:n_pxls*2./3., n_pxls/3.:n_pxls*2./3.], min=min_p)
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 0:n_pxls-1] = BYTSCL(temp_phase, max=max_p, min=min_p)
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, n_pxls+gap:n_pxls+gap+n_pxls-1] = BYTSCL(CONGRID(REFORM(imag_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 2.*(n_pxls+gap):2.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(real_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 3.*(n_pxls+gap):3.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(abs_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 4.*(n_pxls+gap):4.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(amp_phase_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
  ENDFOR

  ;@GJ, 2024/7/24, save the results
  whole_3rd_mp_harmonics_file = path + 'complete_3rd_harmonics.jpg'
  WRITE_JPEG, whole_3rd_mp_harmonics_file, whole_3rd_mp_harmonics

END

;@GJ, 2024/7/21, processing the experimental data from Zhongwei Bian
PRO MPI_load_mat_LEO_Bian_Zhongwei
  
;  filename='test0.mat'
;  filename='test1.mat'
;  filename='test2.mat'
;  filename='test10.mat'
;  filename='test20.mat'
  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\2_10mT_xz_matdata\'
  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\2_15mT_xz_matdata\'
  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\x_axis_cali\'
  path = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\', /MUST_EXIST, TITLE="mat files", /DIRECTORY)
  IF N_ELEMENTS(path) EQ 0 THEN RETURN
  
  file_array = file_search(path, 'test*.mat', count=n_lines_temp) 
  n_lines = n_lines_temp; 31.
  n_pxls = 50. ;number of pixels
  
  abs_3rd = DBLARR(n_lines, n_pxls) * 0.
  real_3rd = DBLARR(n_lines, n_pxls) * 0.
  imag_3rd = DBLARR(n_lines, n_pxls) * 0.
  phase_3rd = DBLARR(n_lines, n_pxls) * 0.
  
  abs_2nd = DBLARR(n_lines, n_pxls) * 0.
  real_2nd = DBLARR(n_lines, n_pxls) * 0.
  imag_2nd = DBLARR(n_lines, n_pxls) * 0.
  phase_2nd = DBLARR(n_lines, n_pxls) * 0.
  
  ;save the reslults
  n_harmonics = 10
  abs_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  real_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  imag_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  phase_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  amp_phase_array = DBLARR(n_harmonics, n_lines, n_pxls) * 0.
  
  n_3rd_mp = 11
  abs_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  real_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  imag_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  phase_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  amp_phase_3rd_mp_array = DBLARR(n_3rd_mp, n_lines, n_pxls) * 0.
  
  ;@GJ, 2024/7/22, fov is 30 mm
  fov = 30. ;mm
  
  FOR k=0, n_lines-1 DO BEGIN
    filename = 'test'+STRING(k, format='(I0)')+'.mat'
    load_mat, filename, path, signal_data, STORE_LEVEL=1
    
    
    N = 50000. ;total time elements
    line_abs = DBLARR(n_pxls, N) * 0.
    line_real = DBLARR(n_pxls, N) * 0.
    line_imag = DBLARR(n_pxls, N) * 0.
    line_phase = DBLARR(n_pxls, N) * 0.

    T = 0.4 ;us
    t_array = FINDGEN(N) * T ;us
    x_array = FINDGEN(n_pxls)/(n_pxls-1)*fov - fov/2. ; fov = 30 mm
    ; N is an integer giving the number of elements in a particular dimension
    ; T is a floating-point number giving the sampling interval
    X = (FINDGEN((N - 1)/2) + 1)
    is_N_even = (N MOD 2) EQ 0
    if (is_N_even) then $
      freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
    else $
      freq = [0.0, X, -(N/2 + 1) + X]/(N*T)

    freq_base = 0.025 ;MHz
    FOR i=0, n_pxls-1 DO BEGIN
      temp_signal = REFORM(signal_data[0:N-1, i])
      temp_signal_fft = FFT(temp_signal)
      line_abs[i, *] = ABS(temp_signal_fft)
      line_real[i, *] = Real_part(temp_signal_fft)
      line_imag[i, *] = Imaginary(temp_signal_fft)
      ;@GJ, 2024/1/5, correcting the phase calculation
      FOR j=0, N-1 DO BEGIN
        IF line_abs[i, j] GT 1e-6 THEN BEGIN
          line_phase[i, j] = 180./!PI*ATAN(temp_signal_fft[j], /phase)
          IF line_phase[i, j] LT 0 THEN line_phase[i, j] += 360.
;          line_phase[i, j] = 180./!PI*ATAN(Imaginary(temp_signal_fft[j]), Real_part(temp_signal_fft[j]))
;          IF line_phase[i, j] GT 180. THEN line_phase[i, j] = 360. - line_phase[i, j]
;          IF line_phase[i, j] LT 0. THEN line_phase[i, j] = 360. + line_phase[i, j]
        ENDIF
      ENDFOR
      IF k EQ FLOOR(n_lines/2) THEN BEGIN
        IF i EQ 0 THEN iplot, freq/freq_base, temp_signal_fft, xrange=[0., 10.], color='red', xtitle='Harmonics', title='Signal FFT', /NO_SAVEPROMPT ELSE iplot, freq/freq_base, temp_signal_fft, /overplot
        ;    IF i EQ 0 THEN iplot, freq/freq_base, REFORM(line_phase[i, *]), xrange=[0., 10.], xtitle='Harmonics', title='Phase', /NO_SAVEPROMPT ELSE iplot, freq/freq_base, REFORM(line_phase[i, *]), /overplot
      ENDIF
    ENDFOR

    ind_3rd = freq_base*3./freq[1]
    abs_3rd[k, *] = REFORM(line_abs[*, ind_3rd])
    real_3rd[k, *] = REFORM(line_real[*, ind_3rd])
    imag_3rd[k, *] = REFORM(line_imag[*, ind_3rd])
    phase_3rd[k, *] = REFORM(line_phase[*, ind_3rd])
    
    ind_2nd = freq_base*2./freq[1]
    abs_2nd[k, *] = REFORM(line_abs[*, ind_2nd])
    real_2nd[k, *] = REFORM(line_real[*, ind_2nd])
    imag_2nd[k, *] = REFORM(line_imag[*, ind_2nd])
    phase_2nd[k, *] = REFORM(line_phase[*, ind_2nd])
    
    ;@GJ, 2024/8/11, calculate the harmonics
    FOR l=0, n_harmonics-1 DO BEGIN
      ind_har = freq_base*(l+1.)/(freq[1]-freq[0])
      abs_array[l, k, *] = REFORM(line_abs[*, ind_har])
      real_array[l, k, *] = REFORM(line_real[*, ind_har])
      imag_array[l, k, *] = REFORM(line_imag[*, ind_har])
      phase_array[l, k, *] = REFORM(line_phase[*, ind_har])
      amp_phase_array[l, k, *] = REFORM(line_abs[*, ind_har]) * REFORM(line_phase[*, ind_har])
    ENDFOR
    
    ;@GJ, 2024/8/11, calculate the 2.5-3.5 Harmonics
    FOR m=0, n_3rd_mp-1 DO BEGIN
      ind_har = freq_base*(2.5+m/10.)/(freq[1]-freq[0])
      abs_3rd_mp_array[m, k, *] = REFORM(line_abs[*, ind_har])
      real_3rd_mp_array[m, k, *] = REFORM(line_real[*, ind_har])
      imag_3rd_mp_array[m, k, *] = REFORM(line_imag[*, ind_har])
      phase_3rd_mp_array[m, k, *] = REFORM(line_phase[*, ind_har])
      amp_phase_3rd_mp_array[m, k, *] = REFORM(line_abs[*, ind_har]) * REFORM(line_phase[*, ind_har])
    ENDFOR

;    iplot, x_array, line_phase[*, ind_3rd], xtitle='Pixels in mm', ytitle='Phase', title='3rd Harmonic Phase', /NO_SAVEPROMPT
;    iplot, x_array, line_abs[*, ind_3rd], xtitle='Pixels in mm', ytitle='Amp/Real/Imag', title='3rd Harmonic Amp'
;    iplot, x_array, line_real[*, ind_3rd], color='red', /overplot;xtitle='Pixels in mm', ytitle='Real', title='3rd Harmonic Phase'
;    iplot, x_array, line_imag[*, ind_3rd], color='blue', /overplot;xtitle='Pixels in mm', ytitle='Imaginary', title='3rd Harmonic Amp'
  ENDFOR
  
  
;  iplot, abs_3rd[*, 19], xtitle='Pixels in mm', title='3rd Harmonics', /NO_SAVEPROMPT
  ;iplot, phase_3rd[*, 19], color='red', xtitle='Pixels in mm', ytitle='degree', title='3rd Harmonics', /NO_SAVEPROMPT;, /overplot
  ;iplot, real_3rd[*, 19]/MAX(real_3rd[*, 19])*200., color='green', /overplot
  ;iplot, imag_3rd[*, 19]/MAX(imag_3rd[*, 19])*200., color='blue', /overplot

  ;@GJ, 2024/7/21, 6th harmonics
  ;iplot, phase_array[5, *, 19], color='red', xtitle='Pixels in mm', ytitle='degree', title='6th Harmonics', /NO_SAVEPROMPT;, /overplot
  ;iplot, real_array[5, *, 19]/MAX(real_array[5, *, 19])*200., color='green', /overplot
  ;iplot, imag_array[5, *, 19]/MAX(imag_array[5, *, 19])*200., color='blue', /overplot
 
  ;@GJ, 2024/7/21, 6th
  ;iplot, REFORM(real_array[5, *, 19]), REFORM(imag_array[5, *, 19]), color='red', sym_index=1, sym_color='blue', xtitle='real',  ytitle='imaginary', title='phase 6th', /NO_SAVEPROMPT
  ;iplot, REFORM(real_3rd[*, 19]), REFORM(imag_3rd[*, 19]), color='red', sym_index=1, sym_color='blue', xtitle='real',  ytitle='imaginary', title='phase 3rd', /NO_SAVEPROMPT
  
  ;@GJ, 2024/7/22, plot the 3rd harmonic phase, ;@GJ, 2024/7/22, fov is 30 mm
  ;iimage, CONGRID(REFORM(imag_array[2, *, *]), n_pxls, n_pxls), title='3rd Imaginary'
  
  abs_3rd_fov = CONGRID(abs_3rd, n_pxls, n_pxls)
  phase_3rd_fov = CONGRID(phase_3rd, n_pxls, n_pxls)
  real_3rd_fov = CONGRID(real_3rd, n_pxls, n_pxls)
  imag_3rd_fov = CONGRID(imag_3rd, n_pxls, n_pxls)
  iimage, abs_3rd_fov, view_title='3rd Amp', VIEW_GRID=[5,7], DIMENSIONS=[n_pxls*6., n_pxls*8.], WINDOW_TITLE='3rd Harmonic', /NO_SAVEPROMPT
  iimage, real_3rd_fov, view_title='3rd Real', /VIEW_NEXT
  iimage, imag_3rd_fov, view_title='3rd Imag', /VIEW_NEXT
  iimage, BYTSCL(phase_3rd_fov), view_title='3rd Phase', /VIEW_NEXT
  iimage, 255.-BYTSCL(phase_3rd_fov), view_title='3rd -Phase', /VIEW_NEXT
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (100.+i*10.)/180.
    new_3rd_fov = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    new_3rd_fov_imag = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    new_3rd_fov_real = imag_3rd_fov * sin(phi) + real_3rd_fov * cos(phi)
;    iplot, new_3rd_fov_real, new_3rd_fov_imag, color='blue', /overplot
    iimage, BYTSCL(new_3rd_fov), view_title='3rd New', /VIEW_NEXT
  ENDFOR
   
  i=12;
  phi = -!PI * i/180.
  new_image_imag = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
  new_image_real = imag_3rd_fov * sin(phi) + real_3rd_fov * cos(phi)
  new_image_abs = SQRT(new_image_imag^2 + new_image_real^2)
  iimage, new_image_imag, title='3rd Imag'
  n=26
  mean_imag = DBLARR(n)
  mean_real = DBLARR(n)
  linear_ratio = DBLARR(n)
  FOR line = 21, 21+n-1 DO BEGIN; No DC offset
    mean_imag[line-21] = MEAN(new_image_imag[line,*])
    mean_real[line-21] = MEAN(new_image_real[line,*])
    fit_result = LINFIT(new_image_real, new_image_imag, CHISQR=fit_chisqr)
    linear_ratio[line-21] = fit_chisqr
    max_imag = MAX(ABS(new_image_imag))
    max_real = MAX(ABS(new_image_real))
    wait, 2
    IF line EQ 21 THEN BEGIN
      iplot, new_image_real[line,*], new_image_imag[line,*], xtitle='Real', ytitle='Imag', title='Scatter Plot', yrange=[-max_imag, max_imag], SYM_INDEX=3, LINESTYLE=0, color='red',  /NO_SAVEPROMPT
    ENDIF ELSE BEGIN
      iplot, new_image_real[line,*], new_image_imag[line,*], xrange=[-max_real, max_real], yrange=[-max_imag, max_imag], SYM_INDEX=2, LINESTYLE=0, color=350000*(line+1), /overplot
    ENDELSE
  ENDFOR
  
  ;@GJ, do subtraction with and without DC offset
  line_a_ind = 23
  line_b_ind = 29
  line_a_imag = REFORM(new_image_imag[line_a_ind, *])
  line_a_real = REFORM(new_image_real[line_a_ind, *])
  line_b_imag = REFORM(new_image_imag[line_b_ind, *])
  line_b_real = REFORM(new_image_real[line_b_ind, *])
  new_image_abs = SQRT(new_image_imag^2 + new_image_real^2)
  iplot, line_b_real, line_b_imag, xtitle='Real', ytitle='Imag', title='Scatter Plot', SYM_INDEX=2, LINESTYLE=0, color='cyan', /NO_SAVEPROMPT
  ratio = MAX(REFORM(new_image_abs[line_a_ind, *])) / MAX(REFORM(new_image_abs[line_b_ind, *]))
  phi = !PI * 2./180.
  new_line_b_imag = -line_b_real * ratio * sin(phi) + line_b_imag * ratio * cos(phi)
  new_line_b_real = line_b_imag * ratio * sin(phi) + line_b_real * ratio * cos(phi)
;  iplot, line_b_real*ratio, line_b_imag*ratio, SYM_INDEX=2, LINESTYLE=0, color=350000*(line_b_ind+3), /overplot
  iplot, new_line_b_real, new_line_b_imag, SYM_INDEX=2, LINESTYLE=0, color='green', /overplot
  iplot, line_a_real, line_a_imag, SYM_INDEX=2, LINESTYLE=0, color='red', /overplot
  diff_imag = line_a_imag - new_line_b_imag
  diff_real = line_a_real - new_line_b_real
  iplot, diff_real, diff_imag, SYM_INDEX=2, LINESTYLE=0, color='blue', /overplot
  iplot, x_array, line_a_imag, xtitle='Pixel #/mm', color='red', /NO_SAVEPROMPT
  iplot, x_array, new_line_b_imag, color='green', /overplot
  iplot, x_array, diff_imag, color='blue', /overplot
  
  ;@GJ, 2024/9/1, modify the result
  FOR i=0,15 DO BEGIN
    new_image_imag[i, *] = diff_imag
  ENDFOR
  iimage, new_image_imag, title='Modified 3rd Imag'
  
  ;@GJ, 2024/8/31, plot the phase under different DC offset
  iplot, mean_real, mean_imag, xtitle='Real', ytitle='Imag', yerror=linear_ratio/MAX(linear_ratio), title='Mean Plot vs DC Offset', /NO_SAVEPROMPT
   
;    i = 6
;    phi = -!PI * (100.+i*10.)/180.
;    new_3rd_1_imag = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
;    new_3rd_1_real = imag_3rd_fov * sin(phi) + real_3rd_fov * cos(phi)
;    iplot, new_3rd_1_real, new_3rd_1_imag, color='red'
;    i = 7
;    phi = -!PI * (100.+i*10.)/180.
;    new_3rd_2_imag = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
;    new_3rd_2_real = imag_3rd_fov * sin(phi) + real_3rd_fov * cos(phi)
;    iplot, new_3rd_2_real, new_3rd_2_imag, color='blue', /overplot
    
  iimage, abs_3rd_fov, view_title='3rd Amp', VIEW_GRID=[5,7], DIMENSIONS=[n_pxls*6., n_pxls*8.], WINDOW_TITLE='3rd Harmonic', /NO_SAVEPROMPT
  iimage, real_3rd_fov, view_title='3rd Real', /VIEW_NEXT
  iimage, imag_3rd_fov, view_title='3rd Imag', /VIEW_NEXT
  iimage, BYTSCL(phase_3rd_fov), view_title='3rd Phase', /VIEW_NEXT
  iimage, 255.-BYTSCL(phase_3rd_fov), view_title='3rd -Phase', /VIEW_NEXT
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * i/180.;-!PI * (310.+i)/180.
    new_3rd_fov = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    new_3rd_fov_imag = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    new_3rd_fov_real = imag_3rd_fov * sin(phi) + real_3rd_fov * cos(phi)
    iimage, BYTSCL(new_3rd_fov), view_title='3rd New', /VIEW_NEXT
  ENDFOR
 ; iplot, phase_3rd_fov-360., ytitle='Phase', title='3rd Phase Angle'+STRING(MEAN(phase_3rd_fov-360.))+'deg', /NO_SAVEPROMPT
;  phase_3rd_temp = CONGRID(phase_3rd, n_pxls, n_pxls)
;  max_p = MAX(phase_3rd_temp[n_pxls/3.:n_pxls*2./3., n_pxls/3.:n_pxls*2./3.], min=min_p)
;  iimage, BYTSCL(phase_3rd_temp, max=max_p, min=min_p), title='3rd Phase'
  delta_line = fov / n_pxls
  ; n_pxls is an integer giving the number of elements in a particular dimension
  ; T is a floating-point number giving the sampling interval
  X_line = (FINDGEN((n_pxls - 1)/2) + 1)
  is_n_pxls_even = (n_pxls MOD 2) EQ 0
  N21 = n_pxls/2 + 1
  if (is_N_even) then $
    freq_line = [0.0, X_line, n_pxls/2, -n_pxls/2 + X_line]/(n_pxls*delta_line) $
  else $
    freq_line = [0.0, X_line, -(n_pxls/2 + 1) + X_line]/(n_pxls*delta_line)
  
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (100.+i*10.)/180.
    new_3rd_fov = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    y_profile = REFORM(new_3rd_fov[24,*])
    IF i EQ 1 THEN iplot, SHIFT(freq_line, -N21), SHIFT(ABS(FFT(y_profile)), -N21), COLOR='red', SYM_INDEX=3, xtitle='Spatial Freq [/mm]', title='Vertical Profile FFT' ELSE iplot, SHIFT(freq_line, -N21), SHIFT(ABS(FFT(y_profile)), -N21), color=35000*(i+1), SYM_INDEX=3, /overplot
    print, (100.+i*10.), ' degree'
;    wait, 1
  ENDFOR
  
  FOR i = 1, 30 DO BEGIN
    phi = -!PI * (310.+i)/180.
    new_3rd_fov = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
    y_profile = REFORM(new_3rd_fov[24,*])
    IF i EQ 1 THEN iplot, SHIFT(freq_line, -N21), SHIFT(ABS(FFT(y_profile)), -N21), COLOR='red', SYM_INDEX=3, xtitle='Spatial Freq [/mm]', title='Vertical Profile FFT' ELSE iplot, SHIFT(freq_line, -N21), SHIFT(ABS(FFT(y_profile)), -N21), color=35000*(i+1), SYM_INDEX=3, /overplot
    print, (310.+i), ' degree'
;    wait, 1
  ENDFOR
  
  ;@GJ, 2024/8/28, plot the vertical line profile
  phi=-!PI*320./180.
  new_3rd_fov = -real_3rd_fov * sin(phi) + imag_3rd_fov * cos(phi)
  y_profile = REFORM(new_3rd_fov[24,*])
  iplot, SHIFT(freq_line, -N21), SHIFT(ABS(FFT(y_profile)), -N21), COLOR='blue', SYM_INDEX=3, xtitle='Spatial Freq [/mm]', title='Vertical Profile FFT', /NO_SAVEPROMPT
  iimage, new_3rd_fov, title='3rd Harmonic'
  
  ;scatter plot the phase
  iplot, real_3rd, imag_3rd, xtitle='Real', ytitle='Imag', SYM_INDEX=1, LINESTYLE=6, title='Scatter Plot', /NO_SAVEPROMPT
  ;@GJ, 2024/8/23, plot the 2nd
  iplot, REFORM(real_array[1, *, *]), REFORM(imag_array[1, *, *]), xtitle='Real', ytitle='Imag', COLOR='red', SYM_INDEX=3, LINESTYLE=6, /OVERPLOT

  ;@GJ, 2024/8/14, get the image
  ind_filtered = WHERE(abs_3rd_fov GT MAX(abs_3rd_fov)*0.20, n_ind, COMPLEMENT=ind_unf, NCOMPLEMENT=n_ind_unf)
  phase_3rd_fov[ind_filtered] *= -1
  phase_3rd_fov[ind_unf] = MIN(phase_3rd_fov[ind_filtered])
 ; iimage, phase_3rd_fov, title='3rd Phase Filtered', /NO_SAVEPROMPT
  phase_filtered_filename = path+'phase_filtered.tif'
  write_tiff, phase_filtered_filename, BYTSCL(phase_3rd_fov), /long

  ;@GJ, 2024/7/24, complete image
  gap = 3.
  whole_harmonics = DBLARR(n_harmonics*(n_pxls+gap), 5.*(n_pxls+gap)) * 0.
  FOR i=0, n_harmonics-1 DO BEGIN
    ;@GJ, 2024/8/7, phase is minus
    temp_phase = CONGRID(REFORM(-phase_array[i, *, *]), n_pxls, n_pxls)
    max_p = MAX(temp_phase[n_pxls/3.:n_pxls*2./3., n_pxls/3.:n_pxls*2./3.], min=min_p)
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 0:n_pxls-1] = BYTSCL(temp_phase, max=max_p, min=min_p)
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, n_pxls+gap:n_pxls+gap+n_pxls-1] = BYTSCL(CONGRID(REFORM(imag_array[i, *, *]), n_pxls, n_pxls))
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 2.*(n_pxls+gap):2.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(real_array[i, *, *]), n_pxls, n_pxls))
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 3.*(n_pxls+gap):3.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(abs_array[i, *, *]), n_pxls, n_pxls))
    whole_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 4.*(n_pxls+gap):4.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(amp_phase_array[i, *, *]), n_pxls, n_pxls))
  ENDFOR
  
  ;@GJ, 2024/7/24, save the results
  whole_harmonics_file = path + 'complete_harmonics.jpg'
  WRITE_JPEG, whole_harmonics_file, whole_harmonics
  
  ;@GJ, 2024/8/11, 2.5 to 3.5th harmonic
  whole_3rd_mp_harmonics = DBLARR(n_3rd_mp*(n_pxls+gap), 5.*(n_pxls+gap)) * 0.
  FOR i=0, n_3rd_mp-1 DO BEGIN
    ;@GJ, 2024/8/7, phase is minus
    temp_phase = CONGRID(REFORM(-phase_3rd_mp_array[i, *, *]), n_pxls, n_pxls)
    max_p = MAX(temp_phase[n_pxls/3.:n_pxls*2./3., n_pxls/3.:n_pxls*2./3.], min=min_p)
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 0:n_pxls-1] = BYTSCL(temp_phase, max=max_p, min=min_p)
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, n_pxls+gap:n_pxls+gap+n_pxls-1] = BYTSCL(CONGRID(REFORM(imag_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 2.*(n_pxls+gap):2.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(real_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 3.*(n_pxls+gap):3.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(abs_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
    whole_3rd_mp_harmonics[i*(n_pxls+gap):i*(n_pxls+gap)+n_pxls-1, 4.*(n_pxls+gap):4.*(n_pxls+gap)+n_pxls-1] = BYTSCL(CONGRID(REFORM(amp_phase_3rd_mp_array[i, *, *]), n_pxls, n_pxls))
  ENDFOR

  ;@GJ, 2024/7/24, save the results
  whole_3rd_mp_harmonics_file = path + 'complete_3rd_harmonics.jpg'
  WRITE_JPEG, whole_3rd_mp_harmonics_file, whole_3rd_mp_harmonics

;   iplot, REFORM(real_3rd[*, 19]), REFORM(imag_3rd[*, 19]), color='red', sym_index=1, sym_color='blue', xtitle='real',  ytitle='imaginary', title='phase'
;  
;  iimage, BYTSCL(abs_2nd), title='2nd Amp'
;  iimage, BYTSCL(real_2nd), title='2nd Real'
;  iimage, BYTSCL(imag_2nd), title='2nd Imag'
;  iimage, BYTSCL(phase_2nd), title='2nd Phase'
   
END


;@GJ, 2024/9/25, Lingying Primate-Size MPI signal processing based on line of excitation offset
PRO primate_MPI_LEO_v01_Hui_Hui

;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\2_10mT_xz_matdata\'
;  path ='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\2_15mT_xz_matdata\'
;  path = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X8mm\'
  path = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X10mm\'
;  path = DIALOG_PICKFILE(PATH='C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\zhongweidata\', /MUST_EXIST, TITLE="mat files", /DIRECTORY)
;  IF N_ELEMENTS(path) EQ 0 THEN RETURN

  temp_file_array = file_search(path, '*.txt', count=n_slices)
  filename_list = LONG(FILE_BASENAME(temp_file_array, '.txt'))
  file_array = temp_file_array[SORT(filename_list)]
  n_tp = 1000000.
  sample_rate = 1000000.;Hz
  delta_t = 1. / sample_rate
  excitation_freq = 5000.; Hz
  n_time_period_per_excitation = sample_rate / excitation_freq
  n_cycle_total = n_tp / n_time_period_per_excitation
  y_freq = 10.; Hz
  n_time_half_y_cycle = sample_rate / y_freq / 2.
  y_n_sample = 20.
  n_cycle_per_y = 1./y_freq/2./y_n_sample * excitation_freq
  n_time_point_per_y = n_cycle_per_y * n_time_period_per_excitation
  y_total = n_tp / n_time_point_per_y
  signal_3rd_harmonic_complex_3d = DCOMPLEXARR(n_slices, y_n_sample, y_n_sample)
  signal_3rd_harmonic_3d_real = DBLARR(n_slices, y_n_sample, y_n_sample)
  signal_3rd_harmonic_3d_imag = DBLARR(n_slices, y_n_sample, y_n_sample)
  signal_3rd_harmonic_3d_phase = DBLARR(n_slices, y_n_sample, y_n_sample)
  
  harmonic_map = DCOMPLEXARR(n_slices, n_time_point_per_y/2) * 0.
  harmonic_50_map = DCOMPLEXARR(n_slices, 50) * 0.
  ;@GJ, 2024/8/11, 2.5 to 3.5th harmonic
  gap = 0.;3.
  whole_3rd_harmonic = DBLARR(n_slices *(y_n_sample+gap), 4.*(y_n_sample+gap)) * 0.
  FOR k=0, n_slices-1 DO BEGIN; k=0, n_slices-1
;  FOR k=n_slices-3, n_slices-3 DO BEGIN
    OPENR, lun, file_array[k], /get_lun
    y_current_ini = 0.
    x_excitation_ini = 0.
    x_voltage_ini = 0.
    readf, lun, y_current_ini, x_excitation_ini, x_voltage_ini, FORMAT='(%"%f%f%f")'
    FREE_LUN, lun
    
    signal_3rd_harmonic = COMPLEXARR(y_total) * 0.
    y_current_startup_array = DBLARR(y_total) * 0.
    y_current_step_array = (ASIN(FINDGEN(21)/10. - 1.) / !PI + 0.5) * n_time_half_y_cycle
    OPENR, lun, file_array[k], /get_lun
    i = 0.
    FOR j=0, y_total-1 DO BEGIN
      num_in_pixel = 0
      REPEAT BEGIN
        y_current_temp = 0.
        x_excitation_temp = 0.
        x_voltage_temp = 0.
        readf, lun, y_current_temp, x_excitation_temp, x_voltage_temp, FORMAT='(%"%f%f%f")'
        IF num_in_pixel EQ 0 THEN BEGIN
          y_current_array = [y_current_temp]
          x_excitation_array = [x_excitation_temp]
          x_voltage_array = [x_voltage_temp]
          y_current_startup_array[j] = y_current_temp
        ENDIF ELSE BEGIN
          y_current_array = [y_current_array, y_current_temp]
          x_excitation_array = [x_excitation_array, x_excitation_temp]
          x_voltage_array = [x_voltage_array, x_voltage_temp]
        ENDELSE
        i++
        num_in_pixel++
        count_index = FLOOR(j/20) * n_time_half_y_cycle + y_current_step_array[j mod 20 + 1]
      ENDREP UNTIL (i GE count_index) ;@GJ, 2024/9/29, read the time points
      
      ;@GJ, to find the singular point
;      IF j EQ 339 THEN BEGIN
;        print, 'test'
;      ENDIF
;      
;      print, 'i of this pixels: ', i
      freq = 1. / (delta_t * num_in_pixel)
      base_index = excitation_freq / freq
      ;@GJ, plot the result
      IF k EQ 0 THEN BEGIN
        IF j EQ y_n_sample/2-1 THEN iplot, FINDGEN(num_in_pixel)* delta_t, x_excitation_array, xtitle='t [s]', ytitle='x drive', xrange=[0, num_in_pixel* delta_t], title='x drive 1/4 cycle'
        IF j EQ y_n_sample-1 THEN iplot, FINDGEN(num_in_pixel)* delta_t, x_excitation_array, xtitle='t [s]', ytitle='x drive', xrange=[0, num_in_pixel* delta_t], title='x drive half cycle'
      ENDIF
;      IF j EQ 0 THEN iplot, FINDGEN(num_in_pixel)/base_index, ABS((FFT(x_voltage_array))), xtitle='Freq (# of Base)', ytitle='Amp', title='FFT', /NO_SAVEPROMPT ELSE iplot, FINDGEN(num_in_pixel)/base_index, ABS((FFT(x_voltage_array))), /overplot
;      signal_3rd_harmonic[j] = (FFT(x_voltage_array))[3. * base_index]
      ;@GJ, 2024/9/30, calculate the 2.9-3.1rd harmonic
      short_3rd = (FFT(x_voltage_array))[2.9*base_index:3.1*base_index]
      real_3rd = MEAN(REAL_PART(short_3rd))
      imag_3rd = MEAN(IMAGINARY(short_3rd))
      signal_3rd_harmonic[j] = DCOMPLEX(real_3rd, imag_3rd)
      ;@GJ, 2024/9/27, calculate the total frequency components
      harmonic_map[k, *] += (FFT(x_voltage_array))[0:n_time_point_per_y/2-1]
      harmonic_50_map[k, *] += (FFT(x_voltage_array))[0:49]
    ENDFOR
    
    FREE_LUN, lun
    print, '# of file: ', k+1
;    IF k EQ FLOOR(n_slices/2) THEN iplot, y_current_startup_array, SYM_INDEX=2, LINESTYLE=0, xtitletitle = 'y current', /NO_SAVEPROMPT
    ;@GJ, 2024/9/30, delete the singular points
    mean_signal_3rd_harmonic = MEAN(ABS(signal_3rd_harmonic))
    stddev_signal_3rd_harmonic = STDDEV(ABS(signal_3rd_harmonic))
    FOR m = 0, y_total-1 DO BEGIN
      IF ABS(signal_3rd_harmonic[m]) GT mean_signal_3rd_harmonic+3.*stddev_signal_3rd_harmonic THEN signal_3rd_harmonic[m] *= 0.
    ENDFOR
    
    ;@GJ, 2024/9/25, get the harmonic signal
    signal_3rd_harmonic_real = REFORM(REAL_PART(signal_3rd_harmonic), y_n_sample, y_n_sample)
    signal_3rd_harmonic_imag = REFORM(IMAGINARY(signal_3rd_harmonic), y_n_sample, y_n_sample)
    
    ;@GJ, 2024/10/1, reverse odd line due to scanning sequence
    FOR i_odd=1, y_n_sample-1, 2 DO BEGIN
      signal_3rd_harmonic_real[*,i_odd] = REVERSE(signal_3rd_harmonic_real[*,i_odd])
      signal_3rd_harmonic_imag[*,i_odd] = REVERSE(signal_3rd_harmonic_imag[*,i_odd])
    ENDFOR
    signal_3rd_harmonic_3d_real[k, *, *] = signal_3rd_harmonic_real
    signal_3rd_harmonic_3d_imag[k, *, *] = signal_3rd_harmonic_imag
    signal_3rd_harmonic_complex = DCOMPLEX(signal_3rd_harmonic_real, signal_3rd_harmonic_imag)
    signal_3rd_harmonic_complex_3d[k, *, *] = signal_3rd_harmonic_complex
    
    ;@GJ, 2024/10/1, do phase calculation
    signal_3rd_harmonic_phase = 180./!PI*ATAN(signal_3rd_harmonic_complex, /phase)
    FOR ip=0, y_n_sample-1 DO BEGIN
      FOR jp=0, y_n_sample-1 DO BEGIN
        IF signal_3rd_harmonic_phase[ip, jp] LT 0 THEN signal_3rd_harmonic_phase[ip, jp] += 360.
      ENDFOR
    ENDFOR
    signal_3rd_harmonic_3d_phase[k, *, *] = signal_3rd_harmonic_phase
    
    IF k EQ FLOOR(n_slices/2) THEN BEGIN
      iimage, CONGRID(ABS(signal_3rd_harmonic_complex), 256, 256), title='3rd Harmonic Amp'
      iimage, BYTSCL(CONGRID(signal_3rd_harmonic_phase, 256, 256)), title='3rd Harmonic Phase'
    ENDIF
    
    whole_3rd_harmonic[k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1, 0:y_n_sample-1] = BYTSCL(signal_3rd_harmonic_phase)
    whole_3rd_harmonic[k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1, y_n_sample+gap:y_n_sample+gap+y_n_sample-1] = BYTSCL(signal_3rd_harmonic_imag)
    whole_3rd_harmonic[k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1, 2.*(y_n_sample+gap):2.*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(signal_3rd_harmonic_real)
    whole_3rd_harmonic[k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1, 3.*(y_n_sample+gap):3.*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(ABS(signal_3rd_harmonic_complex))
  ENDFOR
  
  FOR m_x=0, y_n_sample-1 DO BEGIN
    IF m_x EQ 0 THEN BEGIN
      iplot, REAL_PART(signal_3rd_harmonic_complex_3d[*, *, m_x]), IMAGINARY(signal_3rd_harmonic_complex_3d[*, *, m_x]), xtitle='Real', ytitle='Imag', title='Scatter Plot z_y', SYM_INDEX=3, LINESTYLE=6, color='red',  /NO_SAVEPROMPT
    ENDIF ELSE BEGIN
      iplot, REAL_PART(signal_3rd_harmonic_complex_3d[*, *, m_x]), IMAGINARY(signal_3rd_harmonic_complex_3d[*, *, m_x]), SYM_INDEX=2, LINESTYLE=6, color='red', /overplot
    ENDELSE
  ENDFOR
  iimage, BYTSCL(whole_3rd_harmonic), title='3rd harmonic map xy'
  whole_3rd_harmonics_file = path + '3rd_harmonics_xy.jpg'
  WRITE_JPEG, whole_3rd_harmonics_file, whole_3rd_harmonic
  
  gap = 2.
  whole_3rd_harmonic_zx = DBLARR(n_slices *(y_n_sample+gap), 4.*(y_n_sample+gap)) * 0.
  whole_3rd_harmonic_zy = DBLARR(n_slices *(y_n_sample+gap), 4.*(y_n_sample+gap)) * 0.
  FOR k=0, y_n_sample-1 DO BEGIN
    whole_3rd_harmonic_zx[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, 0:y_n_sample-1] = REFORM(BYTSCL(signal_3rd_harmonic_3d_phase[*, k, *]))
    whole_3rd_harmonic_zx[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, y_n_sample+gap:y_n_sample+gap+y_n_sample-1] = BYTSCL(signal_3rd_harmonic_3d_imag[*, k, *])
    whole_3rd_harmonic_zx[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, 2.*(y_n_sample+gap):2.*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(signal_3rd_harmonic_3d_real[*, k, *])
    whole_3rd_harmonic_zx[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, 3.*(y_n_sample+gap):3.*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(ABS(signal_3rd_harmonic_complex_3d[*, k, *]))

    whole_3rd_harmonic_zy[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, 0:y_n_sample-1] = REFORM(BYTSCL(signal_3rd_harmonic_3d_phase[*, *, k]))
    whole_3rd_harmonic_zy[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, y_n_sample+gap:y_n_sample+gap+y_n_sample-1] = BYTSCL(signal_3rd_harmonic_3d_imag[*, *, k])
    whole_3rd_harmonic_zy[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, 2.*(y_n_sample+gap):2.*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(signal_3rd_harmonic_3d_real[*, *, k])
    whole_3rd_harmonic_zy[k*(y_n_sample+gap):k*(y_n_sample+gap)+n_slices-1, 3.*(y_n_sample+gap):3.*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(ABS(signal_3rd_harmonic_complex_3d[*, *, k]))

  ENDFOR  

  iimage, BYTSCL(whole_3rd_harmonic_zx), title='3rd harmonic map zx'
  whole_3rd_harmonics_file = path + '3rd_harmonics_zx.jpg'
  WRITE_JPEG, whole_3rd_harmonics_file, whole_3rd_harmonic_zx
  
  iimage, BYTSCL(whole_3rd_harmonic_zy), title='3rd harmonic map zy'
  whole_3rd_harmonics_file = path + '3rd_harmonics_zy.jpg'
  WRITE_JPEG, whole_3rd_harmonics_file, whole_3rd_harmonic_zy
  
;  iimage, BYTSCL(CONGRID(ABS(harmonic_map), 256, 500)), title='harmonic map' 
  whole_harmonics_file = path + 'whole_harmonic.jpg'
  WRITE_JPEG, whole_harmonics_file, CONGRID(BYTSCL(ABS(harmonic_map)), 256, 500)  
;  
;  iimage, BYTSCL(CONGRID(ABS(harmonic_50_map), 256, 500)), title='harmonic 50 map'
  whole_harmonics_file = path + 'whole_harmonic50.jpg'
  WRITE_JPEG, whole_harmonics_file, CONGRID(BYTSCL(ABS(harmonic_50_map)), 256, 500)
  
  ;@GJ, 2024/10/9, save the real and imag parts as jpg for denoise
  n_i = 7
  n_j = 3
  temp_image = DBLARR(3, n_i * y_n_sample, n_j * y_n_sample) * 0.
  FOR i=0, n_i-1 DO BEGIN
    FOR j=0, n_j-1 DO BEGIN
      temp_image[0, i*y_n_sample:i*y_n_sample+y_n_sample-1, j*y_n_sample:j*y_n_sample+y_n_sample-1] = REFORM(signal_3rd_harmonic_3d_real[i*j, *, *])
      temp_image[1, i*y_n_sample:i*y_n_sample+y_n_sample-1, j*y_n_sample:j*y_n_sample+y_n_sample-1] = REFORM(signal_3rd_harmonic_3d_imag[i*j, *, *])
      temp_image[2, i*y_n_sample:i*y_n_sample+y_n_sample-1, j*y_n_sample:j*y_n_sample+y_n_sample-1] = REFORM(ABS(signal_3rd_harmonic_complex_3d[i*j, *, *]))
    ENDFOR
  ENDFOR
  
  iimage, BYTSCL(temp_image), title='image for N2N'
  min_real = MIN(temp_image[0, *, *])
  max_real = MAX(temp_image[0, *, *])
  min_imag = MIN(temp_image[1, *, *])
  max_imag = MAX(temp_image[1, *, *])
  top = 255.
  harmonic_3_file = path + 'harmonic3_N2N.jpg'
  WRITE_JPEG, harmonic_3_file, BYTSCL(temp_image), /true
  
  ;@GJ, 2024/10/10, load the N2N denoise data
  fn_n2n = path + 'N2N_b.jpg'
  denoise_image_n2n = read_image(fn_n2n)
  denoise_image_n2n_real = DOUBLE(CONGRID(REFORM(denoise_image_n2n[0,*,*]), n_i * y_n_sample, n_j * y_n_sample))
  denoise_image_n2n_imag = DOUBLE(CONGRID(REFORM(denoise_image_n2n[1,*,*]), n_i * y_n_sample, n_j * y_n_sample))
  FOR i=0, n_i-1 DO BEGIN
    FOR j=0, n_j-1 DO BEGIN
      signal_3rd_harmonic_3d_real[i*j, *, *] = min_real + (max_real - min_real) / (top + 1.) * denoise_image_n2n_real[i*y_n_sample:i*y_n_sample+y_n_sample-1, j*y_n_sample:j*y_n_sample+y_n_sample-1]
      signal_3rd_harmonic_3d_imag[i*j, *, *] = min_imag + (max_imag - min_imag) / (top + 1.) * denoise_image_n2n_imag[i*y_n_sample:i*y_n_sample+y_n_sample-1, j*y_n_sample:j*y_n_sample+y_n_sample-1]
    ENDFOR
  ENDFOR
  
  ;@GJ, 2024/10/10, plot the z-y plane
  FOR i_x=0, y_n_sample-1 DO BEGIN
    IF i_x EQ 0 THEN BEGIN
      iplot, REFORM(signal_3rd_harmonic_3d_real[*, *, i_x]), REFORM(signal_3rd_harmonic_3d_imag[*, *, i_x]), xtitle='Real', ytitle='Imag', title='Scatter Plot N2N (z_y)', xrange=[min_real, max_real], yrange=[min_imag, max_imag], SYM_INDEX=2, LINESTYLE=6, color='red', /NO_SAVEPROMPT
    ENDIF ELSE BEGIN
      iplot, REFORM(signal_3rd_harmonic_3d_real[*, *, i_x]), REFORM(signal_3rd_harmonic_3d_imag[*, *, i_x]), SYM_INDEX=2, LINESTYLE=6, color=350000*(i_x+1), /overplot
    ENDELSE
  ENDFOR
  
  ;do angle correction
  mean_temp = MEAN(signal_3rd_harmonic_complex_3d)
  current_angle = 180./!PI*ATAN(mean_temp, /phase)
  IF current_angle LT 0 THEN current_angle += 360.
  rot_angle = -(90. - current_angle) / 180. * !PI
  signal_3rd_harmonic_3d_imag = -REAL_PART(signal_3rd_harmonic_complex_3d) * sin(rot_angle) + IMAGINARY(signal_3rd_harmonic_complex_3d) * cos(rot_angle)
  signal_3rd_harmonic_3d_real = IMAGINARY(signal_3rd_harmonic_complex_3d) * sin(rot_angle) + REAL_PART(signal_3rd_harmonic_complex_3d) * cos(rot_angle)
  signal_3rd_harmonic_complex_3d = DCOMPLEX(signal_3rd_harmonic_3d_real, signal_3rd_harmonic_3d_imag)
  iplot, signal_3rd_harmonic_3d_real, signal_3rd_harmonic_3d_imag, xtitle='Real', ytitle='Imag', title='Scatter Plot after angle corr', SYM_INDEX=3, LINESTYLE=6, color='red',  /NO_SAVEPROMPT
  ;@GJ, 2024/10/6, save the data for AIMIS3D view
  vol_HU_cube = BYTSCL(CONGRID(signal_3rd_harmonic_3d_imag, 256, 256, 256, /CENTER, /INTERP), MIN=0.)
  vol_HU_cube_resolution = 40. / 255.
  vol_ori_filename = path + 'vol_HU_cubeOri_N2N.sav'
  save, vol_HU_cube, vol_HU_cube_resolution, filename = vol_ori_filename
  
  ;@GJ, 2024/10/2, filter the low amplitude part
  gap = 2.
  rot_whole_3rd_harmonic_zy = DBLARR(39 *(n_slices+gap), y_n_sample*(y_n_sample+gap)) * 0.
  rot_whole_3rd_harmonic_MIP = DBLARR(39 *(n_slices+gap), 3*(y_n_sample+gap)) * 0.
  ;the original image
  FOR k=0, y_n_sample-1 DO BEGIN
    ;@GJ, 2024/10/2, original abs image
    rot_whole_3rd_harmonic_zy[0:n_slices-1, k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1] = REFORM(BYTSCL(ABS(signal_3rd_harmonic_complex_3d[*, *, k])))
    ;@GJ, 2024/10/2, original imag image
    rot_whole_3rd_harmonic_zy[(n_slices+gap):(n_slices+gap)+n_slices-1, k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1] = REFORM(BYTSCL(signal_3rd_harmonic_3d_imag[*, *, k]))   
  ENDFOR
  ;calculate the MIP
  FOR j=0, y_n_sample-1 DO BEGIN
    FOR k=0, y_n_sample-1 DO BEGIN
      rot_whole_3rd_harmonic_MIP[j, k] = MAX(ABS(signal_3rd_harmonic_complex_3d[*, j, k]))
      rot_whole_3rd_harmonic_MIP[(n_slices+gap)+j, k] = MAX(signal_3rd_harmonic_3d_imag[*, j, k])
    ENDFOR
  ENDFOR
  temp = rot_whole_3rd_harmonic_MIP[0:y_n_sample-1, 0:y_n_sample-1]
  rot_whole_3rd_harmonic_MIP[0:y_n_sample-1, 0:y_n_sample-1] = BYTSCL(temp)
  temp = rot_whole_3rd_harmonic_MIP[(n_slices+gap):(n_slices+gap)+y_n_sample-1, 0:y_n_sample-1]
  rot_whole_3rd_harmonic_MIP[(n_slices+gap):(n_slices+gap)+y_n_sample-1, 0:y_n_sample-1] = BYTSCL(temp)
  ;calculate the MIP
  FOR i=0, n_slices-1 DO BEGIN
    FOR k=0, y_n_sample-1 DO BEGIN
      rot_whole_3rd_harmonic_MIP[i, (y_n_sample+gap)+k] = MAX(ABS(signal_3rd_harmonic_complex_3d[i, *, k]))
      rot_whole_3rd_harmonic_MIP[(n_slices+gap)+i, (y_n_sample+gap)+k] = MAX(signal_3rd_harmonic_3d_imag[i, *, k])
    ENDFOR
  ENDFOR
  temp = rot_whole_3rd_harmonic_MIP[0:y_n_sample-1, (y_n_sample+gap):(y_n_sample+gap)+y_n_sample-1]
  rot_whole_3rd_harmonic_MIP[0:y_n_sample-1, (y_n_sample+gap):(y_n_sample+gap)+y_n_sample-1] = BYTSCL(temp)
  temp = rot_whole_3rd_harmonic_MIP[(n_slices+gap):(n_slices+gap)+y_n_sample-1, (y_n_sample+gap):(y_n_sample+gap)+y_n_sample-1]
  rot_whole_3rd_harmonic_MIP[(n_slices+gap):(n_slices+gap)+y_n_sample-1, (y_n_sample+gap):(y_n_sample+gap)+y_n_sample-1] = BYTSCL(temp)
  ;calculate the MIP
  FOR i=0, n_slices-1 DO BEGIN
    FOR j=0, y_n_sample-1 DO BEGIN
      rot_whole_3rd_harmonic_MIP[i, 2*(y_n_sample+gap)+j] = MAX(ABS(signal_3rd_harmonic_complex_3d[i, j, *]))
      rot_whole_3rd_harmonic_MIP[(n_slices+gap)+i, 2*(y_n_sample+gap)+j] = MAX(signal_3rd_harmonic_3d_imag[i, j, *])
    ENDFOR
  ENDFOR
  temp = rot_whole_3rd_harmonic_MIP[0:y_n_sample-1, 2*(y_n_sample+gap):2*(y_n_sample+gap)+y_n_sample-1]
  rot_whole_3rd_harmonic_MIP[0:y_n_sample-1, 2*(y_n_sample+gap):2*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(temp)
  temp = rot_whole_3rd_harmonic_MIP[(n_slices+gap):(n_slices+gap)+y_n_sample-1, 2*(y_n_sample+gap):2*(y_n_sample+gap)+y_n_sample-1]
  rot_whole_3rd_harmonic_MIP[(n_slices+gap):(n_slices+gap)+y_n_sample-1, 2*(y_n_sample+gap):2*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(temp)
  
  ;@GJ, 2024/10/5, calculate the threshold
  c1 = REFORM(signal_3rd_harmonic_3d_real, n_slices * y_n_sample * y_n_sample)
  c2 = REFORM(signal_3rd_harmonic_3d_imag, n_slices * y_n_sample * y_n_sample)
;  array = TRANSPOSE([[c1],[c2]])
;  ; Compute cluster weights, using three clusters:
;  weights = CLUST_WTS(array, N_CLUSTERS = 3)
;  ; Compute the classification of each sample:
;  result = CLUSTER(array, weights, N_CLUSTERS = 3)
;  ; Plot each cluster using a different symbol:
;  IPLOT, array[*, WHERE(result eq 0)], $
;    LINESTYLE = 6, color='red', SYM_INDEX = 2
;  IPLOT, array[*, WHERE(result eq 1)], /OVERPLOT, $
;    LINESTYLE = 6, color='blue', SYM_INDEX = 4
;  IPLOT, array[*, WHERE(result eq 2)], /OVERPLOT, $
;    LINESTYLE = 6, SYM_INDEX = 1
  pdf = HISTOGRAM(c2, NBINS=50, LOCATIONS=bin_ind)
  iplot, bin_ind, pdf, title='imag histogram', xtitle='imag', ytitle='#'
  max_pdf = MAX(pdf, max_pdf_ind)
  min_pdf_10p = MIN(ABS(pdf[max_pdf_ind:*]-0.1*max_pdf), min_10p_ind)
  threshold =  bin_ind[max_pdf_ind + min_10p_ind]
  index_low_imag = WHERE(signal_3rd_harmonic_3d_imag LT threshold, count)
  IF count NE 0 THEN BEGIN
    signal_3rd_harmonic_3d_phase[index_low_imag] *= 0.
    signal_3rd_harmonic_3d_real[index_low_imag] *= 0.
    signal_3rd_harmonic_3d_imag[index_low_imag] *= 0.
    signal_3rd_harmonic_complex_3d[index_low_imag] *= 0.
  ENDIF
  
  ;@GJ, 2024/10/2, scatter plot of the real and imag
  FOR m_x=0, y_n_sample-1 DO BEGIN
    IF m_x EQ 0 THEN BEGIN
      iplot, REAL_PART(signal_3rd_harmonic_complex_3d[*, *, m_x]), IMAGINARY(signal_3rd_harmonic_complex_3d[*, *, m_x]), xtitle='Real', ytitle='Imag', title='Scatter Plot z_y after filter', SYM_INDEX=3, LINESTYLE=6, color='blue',  /NO_SAVEPROMPT
    ENDIF ELSE BEGIN
      iplot, REAL_PART(signal_3rd_harmonic_complex_3d[*, *, m_x]), IMAGINARY(signal_3rd_harmonic_complex_3d[*, *, m_x]), SYM_INDEX=2, LINESTYLE=6, color='blue', /overplot
    ENDELSE
  ENDFOR
  
  IF k EQ FLOOR(n_slices/2) THEN BEGIN
    iimage, CONGRID(ABS(signal_3rd_harmonic_complex_3d), 256, 256), title='3rd Harmonic Amp'
    iimage, BYTSCL(CONGRID(signal_3rd_harmonic_phase, 256, 256)), title='3rd Harmonic Phase'
  ENDIF
  
  ;@GJ, 2024/10/2, do rotation after filter
  FOR i_rot = 0, 36 DO BEGIN
    phi = -!PI * i_rot * 10. / 180.
    new_signal_3rd_harmonic_3d_imag = -signal_3rd_harmonic_3d_real * sin(phi) + signal_3rd_harmonic_3d_imag * cos(phi)
    new_signal_3rd_harmonic_3d_real = signal_3rd_harmonic_3d_imag * sin(phi) + signal_3rd_harmonic_3d_real * cos(phi)
    FOR k=0, y_n_sample-1 DO BEGIN
      rot_whole_3rd_harmonic_zy[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+n_slices-1, k*(y_n_sample+gap):k*(y_n_sample+gap)+y_n_sample-1] = REFORM(BYTSCL(new_signal_3rd_harmonic_3d_imag[*, *, k]))
    ENDFOR
    ;calculate the MIP
    FOR j=0, y_n_sample-1 DO BEGIN
      FOR k=0, y_n_sample-1 DO BEGIN
        rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap)+j, k] = MAX(new_signal_3rd_harmonic_3d_imag[*, j, k])
      ENDFOR
    ENDFOR
    temp = rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+y_n_sample-1, 0:y_n_sample-1]
    rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+y_n_sample-1, 0:y_n_sample-1] = BYTSCL(temp)
    ;calculate the MIP
    FOR i=0, n_slices-1 DO BEGIN
      FOR k=0, y_n_sample-1 DO BEGIN
        rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap)+i, (y_n_sample+gap)+k] = MAX(new_signal_3rd_harmonic_3d_imag[i, *, k])
      ENDFOR
    ENDFOR
    temp = rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+n_slices-1, (y_n_sample+gap):(y_n_sample+gap)+y_n_sample-1]
    rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+n_slices-1, (y_n_sample+gap):(y_n_sample+gap)+y_n_sample-1] = BYTSCL(temp)
    ;calculate the MIP
    FOR i=0, n_slices-1 DO BEGIN
      FOR j=0, y_n_sample-1 DO BEGIN
        rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap)+i, 2*(y_n_sample+gap)+j] = MAX(new_signal_3rd_harmonic_3d_imag[i, j, *])
      ENDFOR
    ENDFOR
    temp = rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+n_slices-1, 2*(y_n_sample+gap):2*(y_n_sample+gap)+y_n_sample-1]
    rot_whole_3rd_harmonic_MIP[(i_rot+2)*(n_slices+gap):(i_rot+2)*(n_slices+gap)+n_slices-1, 2*(y_n_sample+gap):2*(y_n_sample+gap)+y_n_sample-1] = BYTSCL(temp)

    IF i_rot EQ 9 THEN BEGIN
      FOR m_x=0, y_n_sample-1 DO BEGIN
        IF m_x EQ 0 THEN BEGIN
          iplot, new_signal_3rd_harmonic_3d_real[*, *, m_x], new_signal_3rd_harmonic_3d_imag[*, *, m_x], xtitle='Real', ytitle='Imag', title='Scatter Plot z_y after 90 degree', SYM_INDEX=3, LINESTYLE=6, color='green', /NO_SAVEPROMPT
        ENDIF ELSE BEGIN
          iplot, new_signal_3rd_harmonic_3d_real[*, *, m_x], new_signal_3rd_harmonic_3d_imag[*, *, m_x], SYM_INDEX=2, LINESTYLE=6, color='green', /overplot
        ENDELSE
      ENDFOR
      ;@GJ, 2024/10/6, save the data for AIMIS3D view
      vol_HU_cube = BYTSCL(CONGRID(new_signal_3rd_harmonic_3d_imag, 256, 256, 256, /CENTER, /INTERP), MIN=0.)
      vol_HU_cube_resolution = 40. / 255.
      vol_filename = path + 'vol_HU_cubeRot_N2N.sav'
      save, vol_HU_cube, vol_HU_cube_resolution, filename = vol_filename
    ENDIF
  ENDFOR
  
  iimage, BYTSCL(rot_whole_3rd_harmonic_zy), title='Rotation'
  rot_harmonics_file = path + 'rot_3d_zy.jpg'
  WRITE_JPEG, rot_harmonics_file, BYTSCL(rot_whole_3rd_harmonic_zy)
  
  iimage, BYTSCL(rot_whole_3rd_harmonic_MIP[0:13.*(n_slices+gap)+n_slices-1, *]), title='Rotation MIP'
  rot_harmonics_MIP_file = path + 'rot_MIP.jpg'
  WRITE_JPEG, rot_harmonics_MIP_file, BYTSCL(rot_whole_3rd_harmonic_MIP[0:13.*(n_slices+gap)+n_slices-1, *])

END

;@GJ, 2024/10/9, compare the images
PRO read_denoise_image_x10mm
  fn_ori = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X10mm\harmonic3_N2N.jpg'
  image_ori = read_image(fn_ori)
  iplot, REFORM(image_ori[0,*,*]), REFORM(image_ori[1,*,*]), xtitle='Real', ytitle='Imag', title='Scatter Plot original', SYM_INDEX=3, LINESTYLE=6, xrange=[0, 256], yrange=[0, 256], color='red', /NO_SAVEPROMPT
  
  fn_n2n = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X10mm\N2N_b.jpg'
  image_n2n = read_image(fn_n2n)
  iplot, REFORM(image_n2n[0,*,*]), REFORM(image_n2n[1,*,*]), xtitle='Real', ytitle='Imag', title='Scatter Plot N2N', SYM_INDEX=3, LINESTYLE=6, xrange=[0, 256], yrange=[0, 256], color='green', /NO_SAVEPROMPT
  
  fn_dncnn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X10mm\dnCNN_b.jpg'
  image_dncnn = read_image(fn_dncnn)
  iplot, REFORM(image_dncnn[0,*,*]), REFORM(image_dncnn[1,*,*]), xtitle='Real', ytitle='Imag', title='Scatter Plot dnCNN', SYM_INDEX=3, LINESTYLE=6, xrange=[0, 256], yrange=[0, 256], color='blue', /NO_SAVEPROMPT

END

;@GJ, 2024/10/9, compare the images
PRO read_denoise_image_x8mm
  fn_ori = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X8mm\harmonic3_N2N.jpg'
  image_ori = read_image(fn_ori)
  iplot, REFORM(image_ori[0,*,*]), REFORM(image_ori[1,*,*]), xtitle='Real', ytitle='Imag', title='Scatter Plot original', SYM_INDEX=3, LINESTYLE=6, xrange=[0, 256], yrange=[0, 256], color='red', /NO_SAVEPROMPT

  fn_n2n = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X8mm\N2N_b.jpg'
  image_n2n = read_image(fn_n2n)
  iplot, REFORM(image_n2n[0,*,*]), REFORM(image_n2n[1,*,*]), xtitle='Real', ytitle='Imag', title='Scatter Plot N2N', SYM_INDEX=3, LINESTYLE=6, xrange=[0, 256], yrange=[0, 256], color='green', /NO_SAVEPROMPT

  fn_dncnn = 'C:\D_drive\MPI_Tianjie\RESI_Donut_MPI\HuiKunming\20221026-X8mm\dnCNN_b.jpg'
  image_dncnn = read_image(fn_dncnn)
  iplot, REFORM(image_dncnn[0,*,*]), REFORM(image_dncnn[1,*,*]), xtitle='Real', ytitle='Imag', title='Scatter Plot dnCNN', SYM_INDEX=3, LINESTYLE=6, xrange=[0, 256], yrange=[0, 256], color='blue', /NO_SAVEPROMPT

END

;+
; :AUTHOR: Cao zhigang
; :Copyright:CAS-NIGLAS
; :email:zhigang_niglas@163.com
; :blog:blog.sina.com.cn/ahnucao
;-


PRO PRINCOMP,IN_DATA = in_data,LOADINGS = loadings,SCORES = scores,LATENT = latent

  ; Principal component analysis (PCA) on data
  ;
  ; IN_DATA: n*p matrix as input,n-observations, p-variables.
  ; LOADINGS is a p-by-p matrix
  ; LATENT: a vector containing the eigenvalues of the covariance matrix of IN_DATA
  ; SCORES: the principal component scores

  ; Attention: the structure of array in IDL is: col * row, must transpose the matrix according to actual situation.
  ;
  ;
  ; To determin the deminsions of in_data

  IF N_ELEMENTS(SIZE(in_data,/dimension)) NE 2 THEN BEGIN
    PRINT,'Input data must be 2-d form.'
    RETURN
  ENDIF
  ;
  ; Get the col and row
  ;
  dims = SIZE(in_data,/dimensions)
  col = dims[0]
  row = dims[1]
  ;------------------------PCA----------------------------------------------
  ;
  ;1. Normalize the data,i.e., to minus the mean of every column

  avg = FLTARR(col)
  avg = MEAN(in_data,dimension = 2); Get the mean values of every column

  nor_data = FLTARR(col,row)
  ; Normalize
  FOR i=0,col-1 DO BEGIN
    nor_data[i,*] = in_data[i,*] - avg[i]
  ENDFOR

  ;2 Get covariance matrix
  cov = IMSL_COVARIANCES(TRANSPOSE(nor_data))

  ;3 Calc eigvalues and eig-vector by IMSL Advanced Math and Statistics
  eig = IMSL_EIG(cov,vector = eig_matrix)
  ;
  scores = nor_data##FLOAT(eig_matrix)
  ;
  latent = FLOAT(eig)
  loadings =FLOAT(eig_matrix)
END

PRO TEST_PRINCOMP
  ;
  x1 = [2.5,0.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1]
  x2 = [2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,0.9]
  x = [[x1],[x2]]
  PRINCOMP,in_data = TRANSPOSE(x),loadings = loadings,scores= scores,latent= latent
  HELP,loadings,scores,latent
  PRINT,loadings,STRING(13b)
  PRINT,scores,STRING(13b)
  PRINT,latent,STRING(13b)
END