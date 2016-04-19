for sub in 20120307.6833 20120712.3934

do

# Freesurfer
# Folder structure: /usr/local/freesurfer/subjects/ID#/mri/orig/001.mgz
# recon-all -s ID# -all

# ExploreDTI
# .nii are in /data/LAS
# .bval & .bvec are in /data/ORIGINAL
# 2bet folder in ORIGINAL
# 61grads folder in DWIs
# matlab --> set path to /usr/local/bin/exploreDTI --> MainExploreDTI
# Plugins --> Flip/permute dimension(s) of...
## Permute dimensions: 1 2 3
## Flip dimensions: 0 1 0
## Multiple
## .nii: /data/LAS
## output: /data/ORIGINAL

# Script
#cd /data/ORIGINAL
#mv *.nii 2bet/.
#cd 2bet
#bet ${sub}.nii ../${sub}_brain.nii -m -R -f 0.3
#fslmaths ${sub}.nii -mas ../${sub}_brain_mask.nii.gz ../${sub}.nii
#cd ../
#gzip -d ${sub}.nii.gz
#/data/scripts/reorder_dki.sh /data/ORIGINAL /data/REORDERED ${sub}

# ExploreDTI
# Plugins --> Convert --> *.bval/*.bvec
## /data/REORDERED
# Calculate DTI --> Convert raw data to 'DTI *.mat'
## Permute gradient components: y x z
## Flip sign of gradient: x y -z
## B-value: 1000
## Voxel size: 2 2 2
## Number of non-DW images: 7
## Number of DW images: 61
## Matrix size: 128 128 64
## Start converting
## Nifti path: /data/REORDERED/sub.nii
## Text path: /data/REORDERED/sub.txt
## Output path: /data/REORDERED

# Script
#/usr/local/freesurfer/bin/mri_convert /data/freesurfer/${sub}/mri/brain.mgz /data/T1_1mm/${sub}.T1.nii --out_orientation RAS

# ExploreDTI
# Plugins --> Resample 3D/4D *.nii file(s)
## Resolution: 2 2 2
## Interpolator: 2
## Multiple
## Directory: REORDERED
## Folder: T1_1mm
# Settings --> SM/EC/EPI Correction -->
## Also register to other data? --> Yes, to do the EPI correction (non-rigid) --> Suffix: _T1_resampled.nii
## Also register to other data? --> Define image type for registration --> FA
## Registration details for SM/EC correction --> Regularize for low SNR data
## Registration details for SM/EC correction --> Interpolation method --> Cubic spine
# In MATLAB command window:
## MatlabPath = getenv('LD_LIBRARY_PATH');
## MatlabPath = [MatlabPath, ':/usr/local/bin/exploreDTI/Source/MD_cor_E/linux64'];
## setenv('LD_LIBRARY_PATH', MatlabPath)
# Plugins --> Correct for subject motion & EC/EPI distortions --> Multiple
## Uncorrected: /data/REORDERED
## Output: /data/MC_EPI
# Plugins --> Export stuff to *.nii file(s)
## Multiple
## Select:
### DWIs with B0(s) ('_DWIs.nii')
### Gradient direcions ('_grads.txt')
## Folder: /data/MC_EPI
## Output: /data/DWIs

# Terminal
# mv /data/DWIs/*.txt /data/DWIs/61grads
# /data/scripts/convertGrads.sh /data/DWIs/61grads /data/DWIs

# Script
#mkdir -p /data/tracula/${sub}/orig
#cp /data/DWIs/${sub}_MD_C_native_DWIs.nii /data/tracula/${sub}/orig
#cp /data/DWIs/${sub}_MD_C_native_grads68.txt /data/tracula/${sub}/orig
#cd /data/tracula/${sub}/orig
#/usr/local/freesurfer/bin/mri_convert ${sub}_MD_C_native_DWIs.nii ${sub}_MD_C_native_DWIs.nii.gz

# GEDIT
# replace sub in dmrirc.template

# TMUX
# trac-all -c dmrirc.sub -prep
# bedpostx /data/tracula/sub/dmri
# trac-all -c dmrirc.sub -path

done
