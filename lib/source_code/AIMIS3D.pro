;
;  Copyright (c) 2017-2025, Guang Jia @ Xidian University. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
;  NAME:
;       AIMIS3D
;
;  CALLING SEQUENCE:
;       AIMIS3D
;
;  PURPOSE:
;       AI-based Medical Imaging Segmentation for 3D Printing and Naked Eye 3D Visualization
;       AIMIS3D
;
;  ARGUMENTS:
;       NONE
;
;  KEYWORDS:
;       NONE
;
;  MODIFICATION HISTORY:
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/6/3, Original
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/6/25, Fix a bug in sorting and add catch error to DB_LOAD; If error exists about database, please comment this line (INITDatabase,state) in DB_LOAD
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/7/23, Right mouse button to pick point in 3D model and 2D cross-sectional plane images
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/8/3, 1) Add smooth levels with a slider, 2) Let user pick whether to continue to split to 5 objects
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/8/10, 包含如下更新：1。模型还原后，取消自动分块，让用户选择是否分离块；2。第二个界面3D模型窗口实现鼠标中键滚轮---缩放；3。CT3视平面图上选择点与模型上选点实现同步；4。可以添加两个额外模型，画出肿瘤或血管；
;                                                                                     5。添加smooth重建模型的平滑度滚动条；6。在掩模操作过程中，取消了关闭窗口自动重建功能，提供了窗口，让用户选择是否重建。
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/9/6, Fix the bug asking split into 5 objects.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/9/7, Fix the bug of 3D segmentation when using negative min threshold value.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/9/7, Fix the bug of smoothing.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/9/7, Fix the bug of 'MR' modality adding to database for reading MR images.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/9/7, Fix the bug of ruler display when using 'scaling'.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/9/7, Fix the bug of Smooth values.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/11/1, Release storage and varirables.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/11/1, Load stl from main program.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/11/2, Load stl speed was dramatically improved. Now very fast. DCM analysis buttons added to main program.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/11/12, Add folders for image fusion and comparison.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/12/11, If no patient images were selected, the program do nothing instead of error message.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/12/25, Fix stl file saving bug without multipling resolution.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2018/12/29, Compile all necessary files (IDLex...pro) for saving stl in order to use xd_xwt_3D.sav runtime version.    
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/1/2, Add painting brush function to xd_xroi.pro.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/1/4, Replace old colortable by warm colors due to black background.            
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/1/29, Fix the bug of sorting DCM images.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/1/29, 3D view needs front, back, left, right images.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/3, Hot key rts accompanis with button set.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/4, Bug about rotating and naked 3D visualizaiton was fixed. Bug about main window and modify window switching was fixted.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/7, Obj_isa is added to check for modifying; negative scaling using wheel event is disabled.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/8, International languages are being added for selection; Add and Addchild buttons are changed to pictures.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/9, Multiple objects are organized and displayed without overlap.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/10, Left and right mouse button functions are switched.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/11, Speed for drawing circular mask ROI was updated.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/13, (1) Bug about getting image threshold for segmentation was fixed. 
;                                                                                     (2) if there are two many files, only pick 1/3 of the files.
;                                                                                     (3) Either choose database or choose CD to load images.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/14, Circular mask ROI with dynamic diameter changes wwas improved.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/16, Vessel thickness measurement was added with appropriate color and color bar.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/18, For multiple object modification, a structure (struc_modifyobject__define) was defined.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/19, As a 1st step, the structure (struc_modifyobject__define) works in the main interface.       
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/20, As a 2nd step, the structure (struc_modifyobject__define) works in the modify interface.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/22, 3D object could be picked by using a widget_slider with dynamic color change and free color selection (from 30 colors).
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/2/25, Change the draw order, draw the smallest one first, draw the largest one finally
;                                                                                     so the transparent effect for the large object could be guaranted
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/3/21, International version, including Chinease, English, Japanese, languages were added.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/10, A bug about right-mouse-button selecting objects was fixed.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/16, A bug about black-screen appearance in lasso-based 3d segmentation was fixed.
;                                                                                     When clicking redo or undo in modify window, volume display in lower-right corner was hidden.  
;                                                                                     When clicking color box for hide or show single object in modify window, volume display in lower-right corner was hidden or shown.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/17, Color-table for object coloring was adjusted because 8 and 9 are both grey, which was moved to the end.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/19, Circular ROI for painting in xd_xroi.pro was improved by fixing several bugs.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/21, Modifying stl file was significantly improved.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/26, Smoothness level was recorded using thick property of IDLgrPolygon and updated when selecting different object.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/4/27, Circular mask drawing was modified with a faster speed. Bugs about removing layer with circular mask ROI were fixed.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/5/8,  A bug about saving all shown objects as a stl file when object index is 0 was fixed.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/5/19, Some bugs in xd_xroi were fixed and the speed of drawing circular mask ROI was improved.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/5/19, Mask could be displayed after drawing a mask due to location3D mismatch.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/5/31, Negative threshold may affect segmentation in PRO lasso_segmentation_3D, bug was fixed
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/6/1, Negative location3D cause program error in PRO ortho2D_plot, bug was fixed
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/6/5, 1. Software name was changed to AIMIS3D.
;                                                                                    2. After loading the object, hiding baseplate with company icons
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/6/7, 1. Help menu including "AIMIS3D website" and "About AIMIS3D" was added.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/6/14, 1. A bug about loading CD was fixed. 2. A title was added as Progressbar title.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/6/23, Arabic language verison was added into the program.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/8/16, Adding face picture to CT data. 
;                                                                                     disabling draw_3DView by replacing "draw_3DView, xd_sState" by ";draw_3DView, xd_sState" 
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/8/21, Registering face photo to CT reconstructed object.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2019/3/8, Based on AI (U-net algorithm), Brain vessel and vessel wall segmentation and 3D visualization.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/15, VW_CE_reconstruction function, modified to use CE instead of signal difference
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/19, VW_CE_reconstruction, if there is no U-net result, use threshold to do analysis
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/22, AI-based vessel segmentation was added.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/22, Vessel wall enhancement segmentation was added.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/23, Vessel eccentricity measurement was added.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/28, Replace IDLffDICOMex by IDLffDICOM for generating exe file
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/3/28, make_rt, 'aimis3d',  'C:\D_drive\AIMIS_3D\Program', savefile='C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20200328a\AIMIS3D.sav', /overwrite, /DICOMEX, /WIN64, /DATAMINER
;                                                                                     copy 'aimis3d.ini' and 'splash.bmp' to 'CC:\D_drive\AIMIS_3D\Program\aimis3d\'
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/4/13, Fix a bug in dicom_crawler, differentiating MR or non-MR images.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/4/18, Fix a bug in volume_fusion, the probram changed ori_vol_HU_cube_1's value. Now fix the bug.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/4/18, PET/CT data could be loaded and displayed.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/4/18, Fix the bug of 2nd column 3-view diaplay.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/4/19, Fix the bug of window & level setting of 2nd & 3rd column 3-view display.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/5/19, Disable the message about ODBC database missing.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/5/19, A bug about 3-view (last one) ball size is not scaled right was fixted.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/5/23, Vessl thickness was calculated based on smoothing for morph analysis and fast speed.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/5/26, Vessl cross points are calculated when calculating vessel thickness (PICK_VERT_COLORS_by_thickness_smooth).
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/6/24, Generate 3D mask based on stl file for vessel analysis with smooth surface and correct resolution.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/7/05, Ask to reselect dicom image folder if there are less than 5 images. 2) combing vessel thinning and Hough transform together.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/7/07, Brain vessel recongnization was set up, which still needs optimization.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/7/07, make_rt, 'aimis3d',  'C:\D_drive\AIMIS_3D\Program', savefile='C:\D_drive\AIMIS_3D\Software_versions\AIMIS3D_20200707c\AIMIS3D.sav', /overwrite, /DICOMEX, /WIN64, /DATAMINER
;                                                                                     copy 'aimis3d.ini' and 'splash.bmp' to 'CC:\D_drive\AIMIS_3D\Program\aimis3d\'
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/7/12, brain vessels become straight lines using mesh_smooth and separated into line segments, L/R A1 and A2 were segmented.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2020/10/17, 1) Four-view display of loading stl object; 2). Load real size dicom image, instead of icon image. 
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/02/13, Fix a bug of loading CT images with two matrix per dicom file (128*128 and 512*512) in load_dcm_to_stl.pro file.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/02/19, Fix a bug of loading and sorting dicom file in load_dcm_to_stl.pro file. Sort images based on image numer instead of file name.
;                                                                                      Prostate images ADC and T2 could not be displayed appropriately. (C:\D_drive\AIMIS_3D\images\Prostate_Patient778\T2W)
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/02/21, Fix a bug of calculating center position as original_loc in load_dcm_to_stl.pro file for the volume_fusion procedure.
;       ;                                                                              Correct for same orientation of image series with same resolution for registration. 
;       ;                                                                              (1) need to adjust based on resolution; (2) need to adjust based on orientation, axial, sag, or coronnal etc.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/02/22, Fix a bug of resolution adjustment in the volume_fusion procedure.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/02/23, Rewrite the volume_fusion procedure with 3 steps of fusing two volume cubes based on orientation, center point and resolution.
;       ;                                                                              Reference for Euler angle: https://www.pianshen.com/article/48361722399/
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/03/8, Rewrite the procedure, unet_reconstruction_3directions, and works well for vessel 3D restruction based on 3-direction back projection.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/06/30, MPI_01.pro was added to do simulation and analysis for Magnetic Particle Imaging or Tomography.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/07/23, MPI_01.pro, magnetic field was changed from sin function to -cos function in MPI.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/07/28, MPI_filter.pro was added to show the effect of relaxation kernel. MPI simulation was added as a file menu.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/08/24, MPI_01.pro Magnetic particle tomography uisng MRI gradients.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/08/30, MPI_01.pro Using 3rd harmonic and MRI gradient encoding for Magnetic Particle Tomography.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/08/30, MPI_01.pro Relaxation time is calculated and added to signal change.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2021/11/08, MPI_01.pro Fix the bug of non-zero mean signal value for relaxation convolution in non_adiabatic function.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/04/13, MPI_01.pro Add UPEN analysis input and output files.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/04/16, MPI_01.pro Generate multi-color digital phantom with RGB representing AUC, Neel, and Brownian relaxation times respectively
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/05/09, MPI_01.pro Automatically calculate flat ratio and correctly read data from Tian's Lab
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/05/17, MPI_01.pro Correct parameters output in batch_pulsed_fitting procedure
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/06/01,  MPI_01.pro Close the plot objects in batch_pulsed_fitting procedure and correct the output_filename
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/06/08,  Add pro to rename and analyze MRT dat file
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/06/16,  MPI_01.pro Add M vs time curves with and without curve fitting
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/06/16,  dataFileRenaming_widgets.por Add name modification history log to the end of the file 
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/06/17,  dataFileRenaming_widgets.por Display current date when running the program
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2022/06/18,  MPI_01.pro fix the bug for negative M values and dataFileRenaming_widgets.por update the date and particle/buffer/pulse name when loading the old file name
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/01/18,  MPI_01.pro adding the gradient-based rotational drift for spatial encoding in MPI.
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/04/04,  file_io.pro adding the nifti functions (read_nifti and save_nifti) and idlffnifti__define.pro (https://github.com/jdammers/IDL-read-write-Nifti/)
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/05/17,  image_fusion_MPI_CT was added to d_objworld21.pro to generate MPI-CT fused images
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/05/18,  Bugs were fixed in image_fusion_MPI_CT was added to d_objworld21.pro to generate MPI-CT fused images
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/05/18,  DICOM Crawler, "Separate_to_frames.pro", Bugs in "browse_separate" were fixed to modify MPI DICOM format by adding "0002, 00000 File Meta Element Group Len"
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/05/18,  Controlling the fusion widgets' visibility in d_surface1 function of d_objworld21.pro
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/08/04,  Resolution enhancement based on random sampling and twin patter recognization in MPI_01.pro
;       Guang Jia, Xidian University, gjia@xidian.edu.cn, +86-13609118250, 2023/09/16,  Adding dark-spot RESI imaging to improve MPI image resolution in sin_excitation.pro
;       Guang Jia, Xidian University, gjoa@xidian.edu.cn, +86-13609118250, 2023/12/10,  Fully modified to include the offset DC and Gradient fields for focal modulation MPI (FM-MPI) in sin_excitation.pro
;       Guang Jia, Xidian University, gjoa@xidian.edu.cn, +86-13609118250, 2024/04/05,  STED-MPI was fully modified in Sin_Excitation_Chebyshev_MPI.pro. Two bugs about zero tau and wrong harmonic map display for large AC field were fixed.
;       Guang Jia, Xidian University, gjoa@xidian.edu.cn, +86-13609118250, 2024/04/23,  Dog factor optimization and deconvolution was added in Sin_Excitation_Chebyshev_MPI.pro. The pSNR and SSIM can be displayed on reconstructed images.
;       Guang Jia, Xidian University, gjoa@xidian.edu.cn, +86-13609118250, 2024/05/03,  Batch calculation of donut threshold at different excitation fields.
;       Guang Jia, Xidian University, gjoa@xidian.edu.cn, +86-13609118250, 2024/09/12,  Fix the bug of mouse selection of A_DC and A_DC_1; Add a text widget to select any harmonics, e.g. 3.1 or 2.7
;       Guang Jia, Xidian University, gjoa@xidian.edu.cn, +86-13609118250, 2025/02/08,  Update donut and phase modulation using STED-MPI imaging method in sin_excitation.pro and Sin_Excitation_Chebyshev_MPI.pro, double-check the SSIM, pSNR, and MSE accuracy with the help of AI doubao

;
;       ;;
;
;;       
;       ;       ;       ;-

;---------------------------------------------------------------------

PRO AIMIS3D
  
  test_fn = './lib/icon_pics/MPI_test.dcm'
  
  IF (FILE_INFO(test_fn)).EXISTS EQ 0 THEN BEGIN
    msg = ['MPI+ demo version expired. Please contact 779574877@qq.com or wechat 13609118250 for license purchase!']
    ok = dialog_message(msg, /info)
    RETURN
  ENDIF ELSE BEGIN
    IF QUERY_DICOM(test_fn) EQ 1 THEN BEGIN
      obj = OBJ_NEW('IDLffDICOM', test_fn)
      p_modal_temp = STRING(*(obj->GetValue('0008'x,'0060'x,/NO_COPY))[0])
      p_modal = STRCOMPRESS(p_modal_temp, /REMOVE_ALL)
      IF STRCMP(p_modal, 'MPI', 3, /FOLD_CASE) THEN BEGIN
        print, 'licensed user'
      ENDIF ELSE BEGIN
        msg = ['MPI+ demo version expired. Please contact 779574877@qq.com or wechat 13609118250 for license purchase!']
        ok = dialog_message(msg, /info)
        RETURN
      ENDELSE
    ENDIF ELSE BEGIN
      ;limt the use times, demo version
      data_1 = read_binary(test_fn, DATA_START=0,  DATA_DIMS=132)
      data_2 = read_binary(test_fn, DATA_START=132)

      index = WHERE(data_1 EQ 0B, count)
      IF count EQ 0 THEN BEGIN
        msg = ['MPI+ demo version expired. Please contact 779574877@qq.com or wechat 13609118250 for license purchase!']
        ok = dialog_message(msg, /info)
        RETURN
      ENDIF

      j=0
      FOR i=0, 131 DO BEGIN
        IF data_1[i] EQ 0B AND j LT 10 THEN BEGIN
          data_1[i] = 10B
          j++
        ENDIF
      ENDFOR

      OPENW, U, test_fn,/GET_LUN
      WRITEU, U, data_1
      WRITEU, U, data_2
      FREE_LUN, U
    ENDELSE
  ENDELSE
  
  
  ;call object world program
  d_objworld21

END
