%% compute volunteer GM smoothed image
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_GM_fn.nii'), 'C:\mritempfiles\MNI_GM_fn.nii')
    spm('Defaults','PET');
    spm_jobman('initcfg');
%     matlabbatch{1}.spm.spatial.smooth.data = {char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal\', p, '\MRF_VBM\MNI_GM_prob.nii'))};
    matlabbatch{1}.spm.spatial.smooth.data = {'C:\mritempfiles\MNI_GM_fn.nii'};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch, cell(0, 1));
    copyfile('C:\mritempfiles\sMNI_GM_fn.nii', strcat(MRF_path,'\',p,'\MRF_VBM\sMNI_GM_fn.nii'))
end

%% compute patient GM smoothed image
subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    p
    copyfile(strcat(MRF_path,'\',p,'\MRF_VBM\MNI_GM_fn.nii'), 'C:\mritempfiles\MNI_GM_fn.nii')
    spm('Defaults','PET');
    spm_jobman('initcfg');
%     matlabbatch{1}.spm.spatial.smooth.data = {char(strcat('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients\', p, '\MRF_VBM\MNI_GM_prob.nii'))};
    matlabbatch{1}.spm.spatial.smooth.data = {char(strcat('C:\mritempfiles\MNI_GM_fn.nii'))};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch, cell(0, 6));
    copyfile('C:\mritempfiles\sMNI_GM_fn.nii', strcat(MRF_path,'\',p,'\MRF_VBM\sMNI_GM_fn.nii'))
end
%% generate mean extension and std images
subjID=[(studyID)];
mExt = zeros(182,218,182);
sdExt = zeros(182,218,182);
n = 1;
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
%     f1 = 'sMNI_GM_prob.nii';
    f1 = 'sMNI_GM_fn.nii';
    a = niftiread(f1);
    oma = mExt;
    mExt = (mExt.*(n-1) + a) ./ n;
    if n > 1
        sdExt = sqrt(((n-2).*sdExt.^2  +  (a-mExt).*(a-oma))./(n-1));
    end
    n = n + 1
end
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
% save('MRFextension_mean_and_sd.mat','mExt','sdExt');
save('extension_mean_and_sd.mat','mExt','sdExt');

copyfile extension_mean.nii MRFextension_mean.nii
file = load_untouch_nii('extension_mean.nii');
file.img = mExt;
save_untouch_nii(file, 'extension_mean.nii');

copyfile extension_std.nii MRFextension_std.nii
file = load_untouch_nii('extension_std.nii');
file.img = sdExt;
save_untouch_nii(file, 'extension_std.nii');

copyfile('Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\extension_std.nii', 'C:\mritempfiles\extension_std.nii')
spm('Defaults','PET');
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.smooth.data = {char(strcat('Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\extension_std.nii'))};
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run', matlabbatch, cell(0, 1));
copyfile('C:\mritempfiles\sextension_std.nii', 'Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sextension_std.nii')

%% patient extension maps
subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
load extension_mean_and_sd.mat
sdext = double(load_untouch_nii('sextension_std.nii').img);

for p = subjID
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    p
%     copyfile MNI_T1.nii MNI_MRFextenz.nii
%     t1wz = load_untouch_nii('MNI_MRFextenz.nii');
%     GMprob = single(load_untouch_nii('MNI_GM_prob.nii').img).*atlasi;
%     WMprob = single(load_untouch_nii('MNI_WM_prob.nii').img).*atlasi;
%     t1wz.img = (double(load_untouch_nii('sMNI_GM_prob.nii').img)-mExt)./sdext;

    copyfile MNI_T1w.nii MNI_extenz.nii
    t1wz = load_untouch_nii('MNI_extenz.nii');
    GMprob = single(load_untouch_nii('MNI_GM_fn.nii').img).*atlasi;
    WMprob = single(load_untouch_nii('MNI_WM_fn.nii').img).*atlasi;
    t1wz.img = (double(load_untouch_nii('sMNI_GM_fn.nii').img)-mExt)./sdext;

    t1wz.img((GMprob<0.05)&(WMprob<0.05)) = 0;
    t1wz.img(ventmask == 1) = 0;
%     t1wz.img(sdT1w < 15) = 0;
%     save_untouch_nii(t1wz, 'MNI_MRFextenz.nii');
    save_untouch_nii(t1wz, 'MNI_extenz.nii');
end