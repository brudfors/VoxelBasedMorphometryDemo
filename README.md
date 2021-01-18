
# VoxelBasedMorphometryDemo

MATLAB script demonstrating Voxel Based Morphometry (VBM) using SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/). This script is based on John Ashburners VBM Tutorial (https://www.fil.ion.ucl.ac.uk/~john/misc/VBMclass10.pdf). Inputs are assumed a bunch of nifti images (e.g., T1w MRIs), and covariates, these subjects' age and sex. 

The script is run as a pipeline of SPM batch jobs that does:

1. Creates normalised, modulated GM, WM and CSF segmentations, pre-fixed `mwc[1-3]*`, in the same folder as the input images.
2. Compute total intercranial volume (TIV), to correct for brain volume in statistical testing.
3. Gaussian smoothing of GM segmentations, pre-fixed `smwc1*`, in the same folder as the input images
4. Defines a statistical model.
5. Fits the statistical model.

Results can be viewed via the SPM GUI (run `spm` in the MATLAB command window) by pressing the `Results` button, and selecting the generated `SPM.mat` file.
