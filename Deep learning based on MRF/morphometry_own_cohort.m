subjID=[(studyIDs)];
CSFthresh = 0.1
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;

for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    copyfile MNI_GM_fn.nii MNI_binary.nii
    binaryimage = load_untouch_nii('MNI_binary.nii');
    GMprob = single(load_untouch_nii('MNI_GM_fn.nii').img).*atlasi;
    WMprob = single(load_untouch_nii('MNI_WM_fn.nii').img).*atlasi;
    CSFprob = single(load_untouch_nii('MNI_CSF_fn.nii').img).*atlasi;
    BrainProb = GMprob + WMprob;

    upperthresh = 0.75;
    lowerthresh = 0.25;
    bina = zeros(182,218,182);
    bina((GMprob>=lowerthresh)&(GMprob<=upperthresh)&(BrainProb>=0.9)) = 1;
    bina((WMprob>=lowerthresh)&(WMprob<=upperthresh)&(BrainProb>=0.9)) = 1;
    
    kern = ones(5,5,5);
    bina = convn(bina,kern,'same');
    bina(ventmask==1) = 0;
    bina(CSFprob>=CSFthresh) = 0;
    binaryimage.img = bina;
    save_untouch_nii(binaryimage, 'MNI_binary.nii');
    
end
%% patient junction z score
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
subjID=[(studyID)];
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
% load MRFjunction_mean_and_sd.mat
% mT1wjunc = mT1w;
% mmrfT1junc = mT1;
% mmrfT2junc = mT2;
% sdT1wjunc = sdT1w;
% sdmrfT1junc = sdT1;
copyfile junction_mean.nii MRFjunction_mean.nii
copyfile junction_mean_std.nii MRFjunction_mean_std.nii
mT1wjunc = single(load_untouch_nii('junction_mean.nii').img);
sdT1wjunc = double(load_untouch_nii('junction_mean_std.nii').img);
sdT1wjunc(isnan(sdT1wjunc))=0;

for p = subjID
    p
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile MNI_T1z.nii MNI_juncz.nii
  
    t1wz = load_untouch_nii('MNI_juncz.nii');
    t1wzi = double(t1wz.img);


    GMprob = single(load_untouch_nii('MNI_GM_fn.nii').img).*atlasi;
    WMprob = single(load_untouch_nii('MNI_WM_fn.nii').img).*atlasi;
    BrainProb = GMprob + WMprob;
    upperthresh = 0.75;
    lowerthresh = 0.25;
    bina = zeros(182,218,182);
    bina((GMprob>=lowerthresh)&(GMprob<=upperthresh)&(BrainProb>=0.90)) = 1;
    bina((WMprob>=lowerthresh)&(WMprob<=upperthresh)&(BrainProb>=0.90)) = 1;

    t1wz.img = (double(load_untouch_nii('MNI_binary.nii').img)-mT1wjunc)./sdT1wjunc;
    t1wz.img((GMprob<0.05)&(WMprob<0.05)) = 0;
    t1wz.img(ventmask == 1) = 0;
%     t1wz.img(sdT1w < 15) = 0;
    save_untouch_nii(t1wz, 'MNI_juncz.nii');

end