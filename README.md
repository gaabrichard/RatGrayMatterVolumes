#  <img src="assets/logo.png" alt="logo" width="50"/>   RatGrayMatterVolumes            


**RatGrayMatterVolumes** Provides the code and notebooks used for the region-based analysis of gray matter volumes in the rat.


## Processing pipeline            

T2* brain images were converted from Varian .fdf format to NIfTI using the script *fdf2nii.py*.
Typical use : *fdf2nii.py /data/fsems_01.img /data/nii/fsems_01 --scale 10*
Where a scale factor of 10 is applied to the voxel size to facilitate processing using SPM8 tools developed for humans.

Segmentation of gray and white matter was performed in SPM8 (MATLAB) using the SIGMA tissue probability maps.

Regional volumes were calculated using the SIGMA atlas and the *calculate_stats.ipynb* notebook relying on ANTsPy.

To perform the same analysis using the WHS atlas, the SIGMA atlas was registered to the WHS atlas using *register_SIGMA_to_new_atlas.ipynb* (ANTsPy). The linear transformation was applied to the tissue segmentations.
Regional volumes for the WHS atlas were also calculated using *calculate_stats.ipynb*.
   
## Test dataset

The atlases and sample tissue maps (in SIGMA space) are provided in the *data* folder to test the pipeline.



