subjID=["PXX_XXXXX"];
subjIDROI = subjID;
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile T1overT2.nii T1oT2_nmCC.nii
    original = single(load_untouch_nii('T1overT2.nii').img);
    normalized = load_untouch_nii('T1oT2_nmCC.nii');
    ccmask = single(load_untouch_nii('antcorpus_mask_coreg_binary.nii').img);
    normalized.img = original/mean(original(ccmask>0))

    save_untouch_nii(normalized, 'T1oT2_nmCC.nii')

end

%% Generate biasfield corrected T1overT2 
subjID=["PXX_XXXXX"];
subjIDROI = subjID;

subjID = ["P42_13702" "P72_14364" "P33_13509"...
    "P61_14203"]; 
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
%     copyfile T1overT2.nii T1oT2_bico.nii
%     T1bc = single(load_untouch_nii('mT1.nii').img);
%     T2bc = single(load_untouch_nii('mT2im.nii').img);
%     t1ot2 = load_untouch_nii('T1oT2_bico.nii');
%     t1ot2.img = T1bc./T2bc;
%     save_untouch_nii(t1ot2, 'T1oT2_bico.nii');
%     
    copyfile T1overT2.nii T1oT2_nmCC_bico.nii
    original = single(load_untouch_nii('T1oT2_bico.nii').img);
    normalized = load_untouch_nii('T1oT2_nmCC_bico.nii');
    ccmask = single(load_untouch_nii('antcorpus_mask_coreg_binary.nii').img);
    normalized.img = original/mean(original(ccmask>0))
    save_untouch_nii(normalized, 'T1oT2_nmCC_bico.nii');

end
%% normal control normalize
subjID=["VXX_XXXXX"];

MRF_path='Z:\Imaging\Multimodal\Myelin\Normal';
for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile T1overT2.nii T1oT2_nmCC.nii
    original = single(load_untouch_nii('T1overT2.nii').img);
    normalized = load_untouch_nii('T1oT2_nmCC.nii');
    ccmask = single(load_untouch_nii('antcorpus_mask_coreg_binary.nii').img);
    normalized.img = original/mean(original(ccmask>0))
    save_untouch_nii(normalized, 'T1oT2_nmCC.nii')

end
%% file copy ops
subjID=["VXX_XXXXX"];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Normal\';
local_path = 'C:\myelin\';
for p = subjID
%     path = strcat(myelin_path,'\',subjID(p));
%     mkdir(strcat(local_path, p))
    source = strcat(MRF_path,'\', p, '\VBM_newANN_32\T1_dog.nii')
    dest = strcat(myelin_path, p, '\T1.nii')
    copyfile(source, dest)
%     source = strcat(myelin_path,'\', p, '\T2im.nii')
%     dest = strcat(local_path, p, '\T2im.nii')
%     copyfile(source, dest)
end