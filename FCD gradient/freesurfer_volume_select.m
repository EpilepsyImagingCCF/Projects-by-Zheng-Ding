% make cortical mask in freesurfer space

subjID=["PXX_XXXXX"];
subjIDROI = subjID; 
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';


for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile labels.nii cortical_mask.nii
    label = load_untouch_nii('labels.nii');
    labeli = single(label.img);
    
    roi = load_untouch_nii('ROI_T1w.nii');
    roii = single(roi.img);
    roii(roii>0) = 1;
    [roiblob, nBlobs] = bwlabeln(roii);
    centroid = regionprops3(roiblob, "Centroid");
    c = round(centroid.Centroid);
    cx = c(2);
    cy = c(1);
    cz = c(3);
    
%     if cx > 120 
%         labeli((labeli == 2)|((labeli>1000)&(labeli<2000))) = 3000;
%         labeli(labeli < 3000) = 0;
%         labeli(labeli == 3000) = 1;
%     else
        labeli((labeli == 41)|(labeli>2000)) = 3000;
        labeli(labeli < 3000) = 0;
        labeli(labeli == 3000) = 1;
%     end
    
    corm = load_untouch_nii('cortical_mask.nii');
    corm.img = labeli;
    save_untouch_nii(corm, 'cortical_mask.nii')

end
% now go coregister the masks

%% Get corpus callosum

subjID=["PXX_XXXXX"];
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';


for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile labels.nii antcorpus_mask.nii
    label = load_untouch_nii('labels.nii');
    labeli = single(label.img);  
    labeli(labeli ~= 255) = 0;
    labeli(labeli == 255) = 1;
    
    corm = load_untouch_nii('antcorpus_mask.nii');
    corm.img = labeli;
    save_untouch_nii(corm, 'antcorpus_mask.nii')

end

% now go coregister the masks on linux
%% modify coregistered antcorpus_mask
subjID=["PXX_XXXXX"];
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile antcorpus_mask_coreg.nii antcorpus_mask_coreg_binary.nii 
    label = load_untouch_nii('antcorpus_mask_coreg_binary.nii');
    labeli = single(label.img);
    labeli(labeli<1) = 0;
%     labeli(labeli>0) = 1;
   
    label.img = labeli;
    save_untouch_nii(label, 'antcorpus_mask_coreg_binary.nii')

end
%% special case for P83
subjID=["studyXXXXXII"];
subjIDROI = ["P83_XXXXX"]; 
MRF_path='T:\Imaging\Multimodal\Myelin\Patients';

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile labels.nii cortical_mask.nii
    label = load_untouch_nii('labels.nii');
    labeli = single(label.img);
    
    roi = load_untouch_nii('ROI_T1w.nii');
    roii = single(roi.img);
    roii(roii>0) = 1;
    [roiblob, nBlobs] = bwlabeln(roii);
    centroid = regionprops3(roiblob, "Centroid");
    c = round(centroid.Centroid);
    cx = c(2);
    cy = c(1);
    cz = c(3);
    

    labeli((labeli == 41)|(labeli>2000)) = 3000;
    labeli(labeli < 3000) = 0;
    labeli(labeli == 3000) = 1;

    
    corm = load_untouch_nii('cortical_mask.nii');
    corm.img = labeli;
    save_untouch_nii(corm, 'cortical_mask.nii')

end
% now go coregister the masks
%% modify coregistered cortical mask
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';
subjID=["PXX_XXXXX"];
for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile cortical_mask_coreg.nii cortical_mask_coreg_binary.nii 
    label = load_untouch_nii('cortical_mask_coreg_binary.nii');
    labeli = single(label.img);
    labeli(labeli<0.2) = 0;
    labeli(labeli>0) = 1;
   
    label.img = labeli;
    save_untouch_nii(label, 'cortical_mask_coreg_binary.nii')

end

%% Checking the label
subjID=["PXX_XXXXX"];
subjIDROI = ["PXX_XXXXX"]; 
MRF_path='T:\Imaging\Multimodal\Myelin\Patients';

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    label = load_untouch_nii('cortical_mask_coreg_binary.nii');
    labeli = single(label.img);
    checklabel = regionprops3(labeli,'all')
end
