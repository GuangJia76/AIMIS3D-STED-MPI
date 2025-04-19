.RESET_SESSION
.RESET_SESSION

!PATH = Expand_Path('+C:\D_drive\AIMIS_3D\coyote\') + ';' + !PATH

restore,'./lib/source_code/idl_progress_bar.sav'
.compile ./lib/source_code/fsc_color.pro
.compile ./lib/source_code/progressbar__define.pro
.compile ./lib/source_code/rdstl.pro
.compile ./lib/source_code/separate.pro
.compile ./lib/source_code/db_load.pro
.compile ./lib/source_code/using_odbc_excel.pro
.compile ./lib/source_code/xd_xroi.pro
.compile ./lib/source_code/xroi_01.pro
.compile ./lib/source_code/d_surfview.pro
.compile ./lib/source_code/idlexmodelmanip__define.pro
.compile ./lib/source_code/idlexviewmanip__define.pro
.compile ./lib/source_code/idlexrotator__define.pro
.compile ./lib/source_code/idlexobjviewwid__define.pro
.compile ./lib/source_code/idlexwidget__define.pro
.compile ./lib/source_code/idlexinscribingview__define.pro
.compile ./lib/source_code/idlexobjview__define.pro
.compile ./lib/source_code/idlexviewgroup__define.pro
.compile ./lib/source_code/idlexvolview__define.pro
.compile ./lib/source_code/idlexvolviewwid__define.pro
.compile ./lib/source_code/idlgrcolorbar__define.pro
.compile ./lib/source_code/orb__define.pro
.compile ./lib/source_code/struc_modifyobject__define.pro
.compile ./lib/source_code/test_struc_vlineobject.pro
.compile ./lib/source_code/idlffnifti__define.pro
.compile ./lib/source_code/cgDefCharSize.pro
.compile ./lib/source_code/cgQuery.pro
.compile ./lib/source_code/WindowAvailable.pro
.compile ./lib/source_code/cgErrorMsg.pro
.compile ./lib/source_code/cgPercentiles.pro
.compile ./lib/source_code/cgScaleVector.pro
.compile ./lib/source_code/file_io.pro
.compile ./lib/source_code/separate_to_frames.pro
.compile ./lib/source_code/crawler.pro
.compile ./lib/source_code/copy_to_one_direct.pro
.compile ./lib/source_code/convert_to_bmp_images.pro
.compile ./lib/source_code/grow.pro
.compile ./lib/source_code/image_blend.pro
.compile ./lib/source_code/d_objworld21.pro
.compile ./lib/source_code/load_dcm_to_stl.pro
.compile ./lib/source_code/idl_ai_bridge.pro
.compile ./lib/source_code/MPI_01.pro
.compile ./lib/source_code/MPI_load_mat.pro
.compile ./lib/source_code/mpi_filter.pro
.compile ./lib/source_code/sin_excitation.pro
.compile ./lib/source_code/Sin_Excitation_Chebyshev_MPI.pro
.compile ./lib/source_code/pulse_excitation.pro
.compile ./lib/source_code/hk_Anderson_coil.pro
.compile ./lib/source_code/dataFileRenaming_widgets.pro
.compile ./lib/source_code/AIMIS3D.pro

.compile ./lib/source_code/fsc_color.pro
.compile ./lib/source_code/progressbar__define.pro
.compile ./lib/source_code/rdstl.pro
.compile ./lib/source_code/separate.pro
.compile ./lib/source_code/db_load.pro
.compile ./lib/source_code/using_odbc_excel.pro
.compile ./lib/source_code/xd_xroi.pro
.compile ./lib/source_code/xroi_01.pro
.compile ./lib/source_code/d_surfview.pro
.compile ./lib/source_code/idlexmodelmanip__define.pro
.compile ./lib/source_code/idlexviewmanip__define.pro
.compile ./lib/source_code/idlexrotator__define.pro
.compile ./lib/source_code/idlexobjviewwid__define.pro
.compile ./lib/source_code/idlexwidget__define.pro
.compile ./lib/source_code/idlexinscribingview__define.pro
.compile ./lib/source_code/idlexobjview__define.pro
.compile ./lib/source_code/idlexviewgroup__define.pro
.compile ./lib/source_code/idlexvolview__define.pro
.compile ./lib/source_code/idlexvolviewwid__define.pro
.compile ./lib/source_code/idlgrcolorbar__define.pro
.compile ./lib/source_code/orb__define.pro
.compile ./lib/source_code/struc_modifyobject__define.pro
.compile ./lib/source_code/test_struc_vlineobject.pro
.compile ./lib/source_code/idlffnifti__define.pro
.compile ./lib/source_code/cgDefCharSize.pro
.compile ./lib/source_code/cgQuery.pro
.compile ./lib/source_code/WindowAvailable.pro
.compile ./lib/source_code/cgErrorMsg.pro
.compile ./lib/source_code/cgPercentiles.pro
.compile ./lib/source_code/cgScaleVector.pro
.compile ./lib/source_code/file_io.pro
.compile ./lib/source_code/separate_to_frames.pro
.compile ./lib/source_code/crawler.pro
.compile ./lib/source_code/copy_to_one_direct.pro
.compile ./lib/source_code/convert_to_bmp_images.pro
.compile ./lib/source_code/grow.pro
.compile ./lib/source_code/image_blend.pro
.compile ./lib/source_code/d_objworld21.pro
.compile ./lib/source_code/load_dcm_to_stl.pro
.compile ./lib/source_code/idl_ai_bridge.pro
.compile ./lib/source_code/MPI_01.pro
.compile ./lib/source_code/MPI_load_mat.pro
.compile ./lib/source_code/mpi_filter.pro
.compile ./lib/source_code/sin_excitation.pro
.compile ./lib/source_code/Sin_Excitation_Chebyshev_MPI.pro
.compile ./lib/source_code/pulse_excitation.pro
.compile ./lib/source_code/hk_Anderson_coil.pro
.compile ./lib/source_code/dataFileRenaming_widgets.pro
.compile ./lib/source_code/AIMIS3D.pro

RESOLVE_ALL

;IDLITRESOLVEITOOLS
SAVE, /ROUTINES, FILENAME='AIMIS3D.sav'