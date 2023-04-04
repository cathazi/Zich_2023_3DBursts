# Zich_2023_3DBursts

This repository contains the scripts and software to run the data analysis published in:

Spatiotemporal organization of human sensorimotor beta burst activity. Catharina Zich, Andrew J Quinn, James J Bonaiuto, George O'Neill, Lydia C Mardell, Nick S Ward, Sven Bestmann (2023)

paper: https://doi.org/10.7554/eLife.80160

data: https://doi.org/10.17605/OSF.IO/9CQHB



Requirements

The analysis requires the following software running on a Unix-Type operating system.

- MatLab 2022a or greater (may run on earlier versions but not tested)
- Matlab Toolboxes:
    - Statistics and Machine Learning Toolbox
    - Image Processing Toolbox
    - Computer Vision Toolbox

- Additional toolboxes & functions:
    - OSL2 MatLab Toolbox (https://ohba-analysis.github.io/osl-docs/)
    - SPM12 MatLab Toolbox (http://www.fil.ion.ucl.ac.uk/spm/software/download/)
    - MarsBaR region of interest toolbox for SPM (https://marsbar-toolbox.github.io/index.html)
    
    - Circular Statistics Toolbox (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
    - Exact geodesic for triangular meshes (https://uk.mathworks.com/matlabcentral/fileexchange/18168-exact-geodesic-for-triangular-meshes)
    - inhull.m (https://uk.mathworks.com/matlabcentral/fileexchange/10226-inhull)
    - pointbiserial.m (https://uk.mathworks.com/matlabcentral/fileexchange/11222-point-biserial-correlation)
    - TriangleRayIntersection.m (https://uk.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection)


Getting started

Download or clone this repository to your computer
Ensure that you have a recent MatLab with access to the Signal Processing Toolbox and Wavelet Toolbox.
Edit the file paths in the top of hmm_0_initialise to point to the location of these toolboxes on your computer
run hmm_0_initialise in MatLab, if this returns without error then you are good to go. If you see warnings then some dependencies may be missing, follow the instructions in the warning message.
Work through the hmm tutorial scripts

Contents
tss_0_initialise runs the initial setup and configuration for these analyses

tss_1_main_3DBurst.m preprocessin & main 3D Burst analysis
tss_2_get_distances.m get distances for Speed calucaltion
tss_3_main_TW.m quantify direction and speed per 3D Burst

