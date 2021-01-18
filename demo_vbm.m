% This script does VBM on a bunch of T1w MRIs using the SPM batch
% interface. It was based on John Ashburners VBM Tutorial:
%   https://www.fil.ion.ucl.ac.uk/~john/misc/VBMclass10.pdf
%__________________________________________________________________________

clear; clc;

% Parameters
N0 = Inf;                     % Inf uses all images
dir_t1w = '/pth/to/niis';     % data directory (N0 niftis in this directory will be selected)
dir_res = '/pth/to/results';  % results directory
fwhm = 10;                    % amount of smoothing

% Age and sex are covariates, as we here do not have these values, we just 
% pick some random numbers. When specifying them, ensure that they are in
% the same order as the input files read on line 24.
sex = rand(N0,1);
age = rand(N0,1);

%% 1. Unified segmentation
% Normalises+segments the input images, which will be written prefixed 
% mwc[1-3]* (modulated warped GM, WM and CSF) in the same folder as the 
% input images.
files = spm_select('FPList',dir_t1w,'^.*\.nii$');
files = files(1:min(N0,size(files,1)),:);
N     = size(files,1);

% Ensure correct number of covariates
sex = sex(1:N);
age = age(1:N);

matlabbatch = {};
matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(files);
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/home/mbrud/Code/matlab/spm/trunk/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/home/mbrud/Code/matlab/spm/trunk/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/home/mbrud/Code/matlab/spm/trunk/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/home/mbrud/Code/matlab/spm/trunk/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/home/mbrud/Code/matlab/spm/trunk/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/home/mbrud/Code/matlab/spm/trunk/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];                                                                                                                             
spm_jobman('run',matlabbatch);

%% 2. Compute total intercranial volume (TIV)
% Corrects for brain volume in statistical testing.
files_mwc    = cell(3,1);
files_mwc{1} = spm_select('FPList',dir_t1w,'^mwc1.*\.nii$'); files_mwc{1} = files_mwc{1}(1:N,:);
files_mwc{2} = spm_select('FPList',dir_t1w,'^mwc2.*\.nii$'); files_mwc{2} = files_mwc{2}(1:N,:);
files_mwc{3} = spm_select('FPList',dir_t1w,'^mwc3.*\.nii$'); files_mwc{3} = files_mwc{3}(1:N,:);
tiv = zeros(N,1);
vx = 1.5;  % SPM atlas has this voxel size (isotropic)
for n=1:N
    tiv_n = 0;
    for k=1:numel(files_mwc)
        Nii = nifti(files_mwc{k}(n,:));
        tiv_n = tiv_n + sum(Nii.dat(:) > 0.5);
    end
    tiv(n) = vx^3*tiv_n;
end

%% 3. Smooth normalised segmentations
files = spm_select('FPList',dir_t1w,'^mwc1.*\.nii$');
files = files(1:N,:);

matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.data = cellstr(files);
matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);

%% 4. Define statistical model
files = spm_select('FPList',dir_t1w,'^smwc1.*\.nii$');
files = files(1:N,:);

matlabbatch = {};
matlabbatch{1}.spm.stats.factorial_design.dir = {dir_res};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = cellstr(files);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).c = sex;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).cname = 'Sex';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).c = age;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).cname = 'Age';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_user.global_uval = tiv;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',matlabbatch);

%% 5. Fit model
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(dir_res,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
