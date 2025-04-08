%% create new sMNI images for mean of T1 and T2
subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;

t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;
CSFthresh = 0.05;
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
%     copyfile sMNI_T1.nii sMNI_T1_nocsf.nii
%     copyfile sMNI_T2.nii sMNI_T2_nocsf.nii
    T1i = single(load_untouch_nii('MNI_T1.nii').img).*atlasi;
    T2i = single(load_untouch_nii('MNI_T2.nii').img).*atlasi;
    csfi = single(load_untouch_nii('MNI_CSF_prob.nii').img);

    kernelSize = 3; 
    [filteredImage1, filteredImage2] = smooth3Dnonzero(T1i, T2i, kernelSize, csfi);
%     sigma = 6/(2*sqrt(2*log(2)));
%     
%     % Create a 3D Gaussian kernel
%     [x, y, z] = meshgrid(-floor(kernelSize/2):floor(kernelSize/2));
%     gaussianKernel = exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2));
%     gaussianKernel = gaussianKernel / sum(gaussianKernel(:));
%     filteredImage = zeros(size(T1i));
%     
%     % Apply Gaussian filtering
%     for i = 1:size(T1i, 1)
%         for j = 1:size(T1i, 2)
%             for k = 1:size(T1i, 3)
%                 % Define the neighborhood
%                 rowRange = max(1, i-floor(kernelSize/2)):min(size(T1i, 1), i+floor(kernelSize/2));
%                 colRange = max(1, j-floor(kernelSize/2)):min(size(T1i, 2), j+floor(kernelSize/2));
%                 depthRange = max(1, k-floor(kernelSize/2)):min(size(T1i, 3), k+floor(kernelSize/2));
%     
%                 % Extract the neighborhood
%                 neighborhood1 = T1i(rowRange, colRange, depthRange);
%                 neighborhood2 = T2i(rowRange, colRange, depthRange);
%                 csfhood = csfi(rowRange, colRange, depthRange);
% 
%                 mask1 = (neighborhood1 < t1max)&(neighborhood1 > t1min)&(csfhood<CSFthresh);
%                 mask2 = (neighborhood2 < t2max)&(neighborhood2 > t2min)&(csfhood<CSFthresh);
% %                 mask1 = (csfhood<CSFthresh);
% %                 mask2 = (csfhood<CSFthresh);
%                 filteredNeighborhood1 = neighborhood1 .* mask1;
%                 filteredNeighborhood2 = neighborhood2 .* mask2;
% 
%                 gaussianKernelAdjusted = gaussianKernel(1:length(rowRange), 1:length(colRange), 1:length(depthRange));
%                 gaussianKernelAdjusted1 = gaussianKernelAdjusted .* mask1;
%                 gaussianKernelAdjusted2 = gaussianKernelAdjusted .* mask2;
%                 kernelSum1 = sum(gaussianKernelAdjusted1(:));
%                 kernelSum2 = sum(gaussianKernelAdjusted2(:));
%     
%                 if kernelSum1 > 0
%                     gaussianKernelNormalized1 = gaussianKernelAdjusted1 / kernelSum1;
%                     filteredValue1 = sum(sum(sum(filteredNeighborhood1 .* gaussianKernelNormalized1)));
%                 else
%                     filteredValue1 = T1i(i, j, k);
%                 end
% 
%                 if kernelSum2 > 0
%                     gaussianKernelNormalized2 = gaussianKernelAdjusted2 / kernelSum2;
%                     filteredValue2 = sum(sum(sum(filteredNeighborhood2 .* gaussianKernelNormalized2)));
%                 else
%                     filteredValue2 = T2i(i, j, k);
%                 end
%     
%                 filteredImage1(i, j, k) = filteredValue1;
%                 filteredImage2(i, j, k) = filteredValue2;
% 
%             end
%         end
%     end
    
    newim1 = load_untouch_nii('sMNI_T1_nocsf.nii');
    newim1.img = filteredImage1;
    save_untouch_nii(newim1, 'sMNI_T1_nocsf.nii');

    newim2 = load_untouch_nii('sMNI_T2_nocsf.nii');
    newim2.img = filteredImage2;
    save_untouch_nii(newim2, 'sMNI_T2_nocsf.nii');

end

%% patient new smoothed T1 and T2 maps
subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;

t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;
CSFthresh = 0.05;
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    T1i = single(load_untouch_nii('MNI_T1.nii').img);
    T2i = single(load_untouch_nii('MNI_T2.nii').img);
    GMi = single(load_untouch_nii('MNI_GM_prob.nii').img);
    WMi = single(load_untouch_nii('MNI_WM_prob.nii').img);
    csfi = single(load_untouch_nii('MNI_CSF_prob.nii').img);
    copyfile sMNI_T1.nii sMNI_T1_nocsf.nii
    copyfile sMNI_T2.nii sMNI_T2_nocsf.nii
    kernelSize = 3; 
    [filteredImage1, filteredImage2] = smooth3Dnonzero(T1i, T2i, kernelSize, csfi);
%     sigma = 6/(2*sqrt(2*log(2)));
%     
%     % Create a 3D Gaussian kernel
%     [x, y, z] = meshgrid(-floor(kernelSize/2):floor(kernelSize/2));
%     gaussianKernel = exp(-(x.^2 + y.^2 + z.^2)/(2*sigma^2));
%     gaussianKernel = gaussianKernel / sum(gaussianKernel(:));
%     filteredImage = zeros(size(T1i));
%     
%     % Apply Gaussian filtering
%     for i = 1:size(T1i, 1)
%         for j = 1:size(T1i, 2)
%             for k = 1:size(T1i, 3)
%                 % Define the neighborhood
%                 rowRange = max(1, i-floor(kernelSize/2)):min(size(T1i, 1), i+floor(kernelSize/2));
%                 colRange = max(1, j-floor(kernelSize/2)):min(size(T1i, 2), j+floor(kernelSize/2));
%                 depthRange = max(1, k-floor(kernelSize/2)):min(size(T1i, 3), k+floor(kernelSize/2));
%     
%                 % Extract the neighborhood
%                 neighborhood1 = T1i(rowRange, colRange, depthRange);
%                 neighborhood2 = T2i(rowRange, colRange, depthRange);
%                 csfhood = csfi(rowRange, colRange, depthRange);
% 
%                 mask1 = (neighborhood1 < t1max)&(neighborhood1 > t1min)&(csfhood<CSFthresh);
%                 mask2 = (neighborhood2 < t2max)&(neighborhood2 > t2min)&(csfhood<CSFthresh);
%                 filteredNeighborhood1 = neighborhood1 .* mask1;
%                 filteredNeighborhood2 = neighborhood2 .* mask2;
% 
%                 gaussianKernelAdjusted = gaussianKernel(1:length(rowRange), 1:length(colRange), 1:length(depthRange));
%                 gaussianKernelAdjusted1 = gaussianKernelAdjusted .* mask1;
%                 gaussianKernelAdjusted2 = gaussianKernelAdjusted .* mask2;
%                 kernelSum1 = sum(gaussianKernelAdjusted1(:));
%                 kernelSum2 = sum(gaussianKernelAdjusted2(:));
%     
%                 if kernelSum1 > 0
%                     gaussianKernelNormalized1 = gaussianKernelAdjusted1 / kernelSum1;
%                     filteredValue1 = sum(sum(sum(filteredNeighborhood1 .* gaussianKernelNormalized1)));
%                 else
%                     filteredValue1 = T1i(i, j, k);
%                 end
% 
%                 if kernelSum2 > 0
%                     gaussianKernelNormalized2 = gaussianKernelAdjusted2 / kernelSum2;
%                     filteredValue2 = sum(sum(sum(filteredNeighborhood2 .* gaussianKernelNormalized2)));
%                 else
%                     filteredValue2 = T2i(i, j, k);
%                 end
%     
%                 % Store the filtered value
%                 filteredImage1(i, j, k) = filteredValue1;
%                 filteredImage2(i, j, k) = filteredValue2;
% 
%             end
%         end
%     end
    
    newim1 = load_untouch_nii('sMNI_T1_nocsf.nii');
    newim1.img = filteredImage1;
    save_untouch_nii(newim1, 'sMNI_T1_nocsf.nii');

    newim2 = load_untouch_nii('sMNI_T2_nocsf.nii');
    newim2.img = filteredImage2;
    save_untouch_nii(newim2, 'sMNI_T2_nocsf.nii');

end

%% Calculate HC mean and std T1 and T2 
% % subjID=[(studyID)];
% MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
% cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
% atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
% atlasi(atlasi>0) = 1;
% cd('Z:\Imaging\Multimodal\MRF\Peter');
% ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
% atlasi(ventmask==1) = 0;
% 
% t1max = 1759; 
% t2max = 106;
% t1min = 743;
% t2min = 29;
% CSFthresh = 0.2;
% for p = subjID
%     p
%     path = strcat(MRF_path,'\',p,'\MRF_VBM');
%     cd(path)
% %     copyfile MNI_GM_prob.nii MNI_MRFbinary.nii
%     binaryimage = load_untouch_nii('MNI_MRFbinary.nii');
%     mpg = double(load_untouch_nii('MNI_T1w.nii').img).*atlasi;
%     GMprob = single(load_untouch_nii('MNI_GM_prob.nii').img).*atlasi;
%     WMprob = single(load_untouch_nii('MNI_WM_prob.nii').img).*atlasi;
%     BrainProb = GMprob + WMprob;
% 
% %     lowerthresh = mean(mpg(GMprob>=0.5)) + 0.5*std(mpg(GMprob>=0.5));
% %     upperthresh = mean(mpg(WMprob>=0.5)) - 0.5*std(mpg(WMprob>=0.5));
% %     mpg((mpg>=lowerthresh)&(mpg<=upperthresh)) = 1;
% %     mpg(mpg<1) = 0;
% %     bina = mpg;
% %     bina((BrainProb<0.99)) = 0;
% 
% %     upperthresh = mean(GMprob(GMprob>0)) + 0.5*std(GMprob(GMprob>0));
% %     lowerthresh = mean(WMprob(WMprob>0)) - 0.5*std(WMprob(WMprob>0));
%     upperthresh = 0.75;
%     lowerthresh = 0.25;
%     bina = zeros(182,218,182);
%     bina((GMprob>=lowerthresh)&(GMprob<=upperthresh)&(BrainProb>=0.90)) = 1;
%     bina((WMprob>=lowerthresh)&(WMprob<=upperthresh)&(BrainProb>=0.90)) = 1;
%     
%     T1i = single(load_untouch_nii('sMNI_T1_nocsf.nii').img).*bina;
%     T2i = single(load_untouch_nii('sMNI_T2_nocsf.nii').img).*bina;
%     T1i(T1i>t1max) = 0;
%     T2i(T2i>t2max) = 0;
% 
%     kern = ones(5,5,5);
%     bina = convn(bina,kern,'same');
%     bina(ventmask==1) = 0;
%     binaryimage.img = bina;
%     save_untouch_nii(binaryimage, 'MNI_MRFbinary.nii');
%     
% %     kern2 = ones(5,5,5);
% %     T1bina = convn(T1i,kern2,'same');
%     T1bina = T1i;
%     T1bina(ventmask==1) = 0;
% %     copyfile MNI_GM_prob.nii MNI_MRFbinary_T1.nii
%     binaryimage = load_untouch_nii('MNI_MRFbinary_T1.nii');
%     binaryimage.img = T1bina;
%     save_untouch_nii(binaryimage, 'MNI_MRFbinary_T1.nii');
% 
% %     T2bina = convn(T2i,kern2,'same');
%     T2bina = T2i;
%     T2bina(ventmask==1) = 0;
% %     copyfile MNI_GM_prob.nii MNI_MRFbinary_T2.nii
%     binaryimage = load_untouch_nii('MNI_MRFbinary_T2.nii');
%     binaryimage.img = T2bina;
%     save_untouch_nii(binaryimage, 'MNI_MRFbinary_T2.nii');
% end
% calculate binary image with motion mask

subjID=[(studyID)];
mT1w = zeros(182,218,182);
sdT1w = zeros(182,218,182);
mT1 = zeros(182,218,182);
sdT1 = zeros(182,218,182);
mT2 = zeros(182,218,182);
sdT2 = zeros(182,218,182);
n_matrix = zeros(182,218,182);
n_matrix_p = zeros(182,218,182);
n_matrix_p2 = zeros(182,218,182);

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';

for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    p
    GMprob = single(load_untouch_nii('MNI_GM_prob.nii').img).*atlasi;
    WMprob = single(load_untouch_nii('MNI_WM_prob.nii').img).*atlasi;
    BrainProb = GMprob + WMprob;

    f1 = 'MNI_MRFbinary.nii';
%     f2 = 'sMNI_T1_nocsf.nii';
%     f3 = 'sMNI_T2_nocsf.nii';
    f2 = 'MNI_T1.nii';
    f3 = 'MNI_T2.nii';
    f4 = 'MNI_motion_mask.nii';

    a = niftiread(f1);
    b = niftiread(f2);
    c = niftiread(f3);
    d = niftiread(f4);
    d(b==0) = 0;
%     d(b==0) = 0;
%     d(c==0) = 0;

    n_matrix_p2 = n_matrix_p;
    n_matrix_p = n_matrix;
    n_matrix = n_matrix+d;

    oma = mT1w;
    omb = mT1;
    omc = mT2;

    mT1w = (mT1w.*(n_matrix_p) + a.*d) ./ n_matrix;
    mT1 = (mT1.*(n_matrix_p) + b.*d) ./ n_matrix;
    mT2 = (mT2.*(n_matrix_p) + c.*d) ./ n_matrix;

    if (p~="V04_12578")
        sdT1w(d>0) = sqrt(((n_matrix_p2(d>0)).*sdT1w(d>0).^2  +  (a(d>0)-mT1w(d>0)).*(a(d>0)-oma(d>0)))./(n_matrix_p(d>0)));
        sdT1(d>0) = sqrt(((n_matrix_p2(d>0)).*sdT1(d>0).^2  +  (b(d>0)-mT1(d>0)).*(b(d>0)-omb(d>0)))./(n_matrix_p(d>0)));
        sdT2(d>0) = sqrt(((n_matrix_p2(d>0)).*sdT2(d>0).^2  +  (c(d>0)-mT2(d>0)).*(c(d>0)-omc(d>0)))./(n_matrix_p(d>0)));
    end
end
%
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
% save('MRFjunction_mean_and_sd.mat','mT1w','mT1', 'mT2', 'sdT1w' ,'sdT1','sdT2');

% file = load_untouch_nii('junction_mean.nii');
% file.img = zeros(182,218,182);
% file.img = mT1w;
% save_untouch_nii(file, 'junction_mean.nii');

% file = load_untouch_nii('junction_mean_std.nii');
% file.img = zeros(182,218,182);
% file.img = sdT1w;
% save_untouch_nii(file, 'junction_mean_std.nii');

% was junction_T1_mean.nii and junction_T2_mean.nii for submission to Neurology
% , but really it is just T1 mean and T2 mean

file = load_untouch_nii('T1_mean_nosmooth.nii');
file.img = zeros(182,218,182);
file.img = mT1;
save_untouch_nii(file, 'T1_mean_nosmooth.nii');

file = load_untouch_nii('T1_std_nosmooth.nii');
file.img = zeros(182,218,182);
file.img = sdT1;
save_untouch_nii(file, 'T1_std_nosmooth.nii');

file = load_untouch_nii('T2_mean_nosmooth.nii');
file.img = zeros(182,218,182);
file.img = mT2;
save_untouch_nii(file, 'T2_mean_nosmooth.nii');

file = load_untouch_nii('T2_std_nosmooth.nii');
file.img = zeros(182,218,182);
file.img = sdT2;
save_untouch_nii(file, 'T2_std_nosmooth.nii');
%
% copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_mean_std.nii C:\MRFzmaps\junction_mean_std.nii
% copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_T1_mean.nii C:\MRFzmaps\junction_T1_mean.nii
% copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_T2_mean.nii C:\MRFzmaps\junction_T2_mean.nii
% copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_T1_std.nii C:\MRFzmaps\junction_T1_std.nii
% copyfile Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\junction_T2_std.nii C:\MRFzmaps\junction_T2_std.nii
% 
% spm('Defaults','PET');
% spm_jobman('initcfg');
% matlabbatch{1}.spm.spatial.smooth.data = {
%     char(strcat('C:\MRFzmaps\junction_mean_std.nii'))
%     char(strcat('C:\MRFzmaps\junction_T1_mean.nii'))
%     char(strcat('C:\MRFzmaps\junction_T2_mean.nii'))
%     char(strcat('C:\MRFzmaps\junction_T1_std.nii'))
%     char(strcat('C:\MRFzmaps\junction_T2_std.nii'))
%     };
% matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
% matlabbatch{1}.spm.spatial.smooth.dtype = 0;
% matlabbatch{1}.spm.spatial.smooth.im = 0;
% matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% 
% spm_jobman('run', matlabbatch, cell(0, 8));
% 
% copyfile C:\MRFzmaps\sjunction_mean_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sjunction_mean_std.nii
% copyfile C:\MRFzmaps\sjunction_T1_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sjunction_T1_mean.nii 
% copyfile C:\MRFzmaps\sjunction_T2_mean.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sjunction_T2_mean.nii 
% copyfile C:\MRFzmaps\sjunction_T1_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sjunction_T1_std.nii 
% copyfile C:\MRFzmaps\sjunction_T2_std.nii Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps\sjunction_T2_std.nii 

%% patient junction maps

CSFthresh = 0.1
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;
subjID=[(studyID)];
for p = subjID
    p
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    copyfile MNI_GM_prob.nii MNI_MRFbinary.nii
    binaryimage = load_untouch_nii('MNI_MRFbinary.nii');
%     mpg = single(load_untouch_nii('MNI_T1w.nii').img).*atlasi;
    GMprob = single(load_untouch_nii('MNI_GM_prob.nii').img).*atlasi;
    WMprob = single(load_untouch_nii('MNI_WM_prob.nii').img).*atlasi;
    CSFprob = single(load_untouch_nii('MNI_CSF_prob.nii').img).*atlasi;
    BrainProb = GMprob + WMprob;

%     lowerthresh = mean(mpg(GMprob>=0.5)) + 0.5*std(mpg(GMprob>=0.5));
%     upperthresh = mean(mpg(WMprob>=0.5)) - 0.5*std(mpg(WMprob>=0.5));
%     mpg((mpg>=lowerthresh)&(mpg<=upperthresh)) = 1;
%     mpg(mpg<1) = 0;
%     bina = mpg;
%     bina((BrainProb<0.99)) = 0;

%     upperthresh = mean(GMprob(GMprob>0)) + 0.5*std(GMprob(GMprob>0));
%     lowerthresh = mean(WMprob(WMprob>0)) - 0.5*std(WMprob(WMprob>0));
    upperthresh = 0.75;
    lowerthresh = 0.25;
    bina = zeros(182,218,182);
    bina((GMprob>=lowerthresh)&(GMprob<=upperthresh)&(BrainProb>=0.9)) = 1;
    bina((WMprob>=lowerthresh)&(WMprob<=upperthresh)&(BrainProb>=0.9)) = 1;
    
    T1i = single(load_untouch_nii('sMNI_T1_nocsf.nii').img).*bina;
    T2i = single(load_untouch_nii('sMNI_T2_nocsf.nii').img).*bina;
    T1i(T1i>t1max) = 0;
    T2i(T2i>t2max) = 0;

    kern = ones(5,5,5);
    bina = convn(bina,kern,'same');
    bina(ventmask==1) = 0;
    bina(CSFprob>=CSFthresh) = 0;
    binaryimage.img = bina;
    save_untouch_nii(binaryimage, 'MNI_MRFbinary.nii');
    
%     T1bina = convn(T1i,kern2,'same');
%     T1bina = T1i;
%     T1bina(ventmask==1) = 0;
%     T1bina(CSFprob>=CSFthresh) = 0;
%     copyfile MNI_GM_prob.nii MNI_MRFbinary_T1.nii
%     binaryimage = load_untouch_nii('MNI_MRFbinary_T1.nii');
%     binaryimage.img = T1bina;
%     save_untouch_nii(binaryimage, 'MNI_MRFbinary_T1.nii');

%     T2bina = convn(T2i,kern2,'same');
%     T2bina = T2i;
%     T2bina(ventmask==1) = 0;
%     T1bina(CSFprob>=CSFthresh) = 0;
%     copyfile MNI_GM_prob.nii MNI_MRFbinary_T2.nii
%     binaryimage = load_untouch_nii('MNI_MRFbinary_T2.nii');
%     binaryimage.img = T2bina;
%     save_untouch_nii(binaryimage, 'MNI_MRFbinary_T2.nii');

end
%% patient junction z score and T1 T2z
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;

subjID=[(studyID)];
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
% load MRFjunction_mean_and_sd.mat
% mT1wjunc = mT1w;
% mmrfT1junc = mT1;
% mmrfT2junc = mT2;
% sdT1wjunc = sdT1w;
% sdmrfT1junc = sdT1;
% mT1wjunc = single(load_untouch_nii('junction_mean.nii').img);
% mmrfT1junc = single(load_untouch_nii('junction_T1_mean.nii').img);
% sdT1wjunc = double(load_untouch_nii('junction_mean_std.nii').img);
% sdmrfT1junc = double(load_untouch_nii('sjunction_T1_std.nii').img);
% sdmrfT2junc = double(load_untouch_nii('sjunction_T2_std.nii').img);


% mmrfT1junc = double(load_untouch_nii('T1_mean_nocsf.nii').img);
% mmrfT2junc = double(load_untouch_nii('sT2_mean.nii').img);
% sdmrfT1junc = double(load_untouch_nii('sT1_std.nii').img);
% sdmrfT2junc = double(load_untouch_nii('sT2_std.nii').img);
mmrfT1junc = double(load_untouch_nii('T1_mean_nocsf.nii').img);
mmrfT2junc = double(load_untouch_nii('T2_mean_nocsf.nii').img);
sdmrfT1junc = double(load_untouch_nii('T1_std_nocsf.nii').img);
sdmrfT2junc = double(load_untouch_nii('T2_std_nocsf.nii').img);

% !!!!!!!usually use this section please
% mmrfT1junc = double(load_untouch_nii('junction_T1_mean.nii').img);
% mmrfT1junc(isnan(mmrfT1junc))=0;
% mmrfT2junc = double(load_untouch_nii('junction_T2_mean.nii').img);
% mmrfT2junc(isnan(mmrfT2junc))=0;
% sdmrfT1junc = double(load_untouch_nii('junction_T1_std.nii').img);
% sdmrfT1junc(isnan(sdmrfT1junc))=0;
% sdmrfT2junc = double(load_untouch_nii('junction_T2_std.nii').img);
% sdmrfT2junc(isnan(sdmrfT2junc))=0;
% csfi = double(load_untouch_nii('CSF_mean.nii').img);
% [sdmrfT1junc, sdmrfT2junc] = smooth3Dnonzero_std(double(load_untouch_nii('junction_T1_mean.nii').img), double(load_untouch_nii('junction_T2_mean.nii').img), sdmrfT1junc, sdmrfT2junc, 3, csfi);

% newim1 = load_untouch_nii('sjunction_T1_mean.nii');
% newim1.img = mmrfT1junc;
% save_untouch_nii(newim1, 'sjunction_T1_mean.nii');
% 
% newim2 = load_untouch_nii('sjunction_T2_mean.nii');
% newim2.img = mmrfT2junc;
% save_untouch_nii(newim2, 'sjunction_T2_mean.nii');
% 
% newim1 = load_untouch_nii('sjunction_T1_std.nii');
% newim1.img = sdmrfT1junc;
% save_untouch_nii(newim1, 'sjunction_T1_std.nii');
% 
% newim2 = load_untouch_nii('sjunction_T2_std.nii');
% newim2.img = sdmrfT2junc;
% save_untouch_nii(newim2, 'sjunction_T2_std.nii');
%
for p = subjID
    p
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
%     copyfile MNI_T1.nii MNI_MRFjuncz.nii
%     copyfile MNI_T1.nii MNI_MRFT1juncz.nii
%     copyfile MNI_T1.nii MNI_MRFT2juncz.nii
   
%     t1wz = load_untouch_nii('MNI_MRFjuncz.nii');
%     t1wzi = double(t1wz.img);
%     t1z = load_untouch_nii('MNI_MRFT1juncz.nii');
%     t1zi = double(t1z.img);
%     t2z = load_untouch_nii('MNI_MRFT2juncz.nii');
%     t2zi = double(t2z.img);


%     copyfile MNI_GM_prob.nii MNI_T1z.nii
%     copyfile MNI_GM_prob.nii MNI_T2z.nii
%     t1z = load_untouch_nii('MNI_T1z.nii');
%     t1zi = double(t1z.img);
%     t2z = load_untouch_nii('MNI_T2z.nii');
%     t2zi = double(t2z.img);
    copyfile MNI_GM_prob.nii MNI_T1z_nosmth.nii
    copyfile MNI_GM_prob.nii MNI_T2z_nosmth.nii
    t1z = load_untouch_nii('MNI_T1z_nosmth.nii');
    t1zi = double(t1z.img);
    t2z = load_untouch_nii('MNI_T2z_nosmth.nii');
    t2zi = double(t2z.img);

    GMprob = single(load_untouch_nii('MNI_GM_prob.nii').img).*atlasi;
    WMprob = single(load_untouch_nii('MNI_WM_prob.nii').img).*atlasi;
    BrainProb = GMprob + WMprob;
    upperthresh = 0.75;
    lowerthresh = 0.25;
    bina = zeros(182,218,182);
    bina((GMprob>=lowerthresh)&(GMprob<=upperthresh)&(BrainProb>=0.90)) = 1;
    bina((WMprob>=lowerthresh)&(WMprob<=upperthresh)&(BrainProb>=0.90)) = 1;

%     t1wz.img = (double(load_untouch_nii('MNI_MRFbinary.nii').img)-mT1wjunc)./sdT1wjunc;
%     t1z.img = (double(load_untouch_nii('sMNI_T1_nocsf.nii').img)-mmrfT1junc)./sdmrfT1junc.*bina;
%     t2z.img = (double(load_untouch_nii('sMNI_T2_nocsf.nii').img)-mmrfT2junc)./sdmrfT2junc.*bina;
% !!!!! use this
%     t1wz.img = (double(load_untouch_nii('MNI_MRFbinary.nii').img)-mT1wjunc)./sdT1wjunc;
%     t1z.img = (double(load_untouch_nii('sMNI_T1_nocsf.nii').img)-mmrfT1junc)./sdmrfT1junc;
%     t2z.img = (double(load_untouch_nii('sMNI_T2_nocsf.nii').img)-mmrfT2junc)./sdmrfT2junc;


    t1z.img = (double(load_untouch_nii('MNI_T1.nii').img)-mmrfT1junc)./sdmrfT1junc;
    t2z.img = (double(load_untouch_nii('MNI_T2.nii').img)-mmrfT2junc)./sdmrfT2junc;


%     t1wz.img((GMprob<0.05)&(WMprob<0.05)) = 0;
%     t1wz.img(ventmask == 1) = 0;

    t1z.img((GMprob<0.05)&(WMprob<0.05)) = 0;
    t1z.img(t1z.img>10) = 0;
    t1z.img(t1z.img<-10) = 0;
    t1z.img(ventmask == 1) = 0;

    t2z.img((GMprob<0.05)&(WMprob<0.05)) = 0;
    t2z.img(t2z.img>10) = 0;
    t2z.img(t2z.img<-10) = 0;
    t2z.img(ventmask == 1) = 0;
%     t2z.img(sdmrfT2junc < 15) = 0;
%     save_untouch_nii(t1wz, 'MNI_MRFjuncz.nii');
%     save_untouch_nii(t1z, 'MNI_MRFT1juncz.nii');
%     save_untouch_nii(t2z, 'MNI_MRFT2juncz.nii');
%     save_untouch_nii(t1wz, 'MNI_MRFjuncz.nii');
    save_untouch_nii(t1z, 'MNI_T1z_nosmth.nii');
    save_untouch_nii(t2z, 'MNI_T2z_nosmth.nii');

end

%% threshold junction map and look at some metrics
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
subjID=[(studyID)];

foundlistjunc = [];
NumFPlistjunc = [];  
foundlistT1 = [];
NumFPlistT1 = []; 
foundlistT2 = [];
NumFPlistT2 = []; 
foundlistext = [];
% NumFPlist = [];
thresh = 2.5;
% threshper = 99.9;
voxelthresh = 5;
for p = subjID
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    p
    cd(path)
    t1zscore = single(load_untouch_nii('MNI_T1z.nii').img);
    t1zscore(isnan(t1zscore)) = 0; 
    t1zscore(t1zscore<thresh) = 0;
    t1zscore(t1zscore>=8) = 0;
    t1zscore(t1zscore>=thresh) = 1;
    

%     pixelValues = t1zscore(:);
%     thresholdValue = prctile(pixelValues, threshper);
%     t1zscore = t1zscore> thresholdValue;
%     copyfile MNI_MRFT1juncz.nii MNI_T1z_binary.nii

    file = load_untouch_nii('MNI_T1z_binary.nii');
    file.img = medfilt3(t1zscore);
    save_untouch_nii(file, 'MNI_T1z_binary.nii');

    roii = single(load_untouch_nii('MNI_ROI_final.nii').img);
    roii(isnan(roii)) = 0; 
    
    overlap = file.img & roii;
    
    if sum(overlap(:)) > voxelthresh
        foundlistT1 = vertcat(foundlistT1, 1);
    else
        foundlistT1 = vertcat(foundlistT1, 0);
    end
%     NumFP = numberOfBlobs - size(trueind,1);
%     NumFPlistjunc = vertcat(NumFPlistjunc, NumFP);   
    

    junczscore = single(load_untouch_nii('MNI_MRFjuncz.nii').img);
    junczscore(isnan(junczscore)) = 0; 
    junczscore(junczscore<thresh) = 0;
    junczscore(junczscore>=8) = 0;
    junczscore(junczscore>=thresh) = 1;

%     pixelValues = junczscore(:);
%     thresholdValue = prctile(pixelValues, threshper);
%     junczscore = junczscore > thresholdValue;
%     copyfile MNI_MRFjuncz.nii MNI_juncz_binary.nii
%     file = load_untouch_nii('MNI_juncz_binary.nii');
% %     file.img = junczscore;
%     file.img = medfilt3(junczscore);
%     save_untouch_nii(file, 'MNI_juncz_binary.nii')
    
    overlap = file.img & roii;
    if sum(overlap(:)) > voxelthresh
        foundlistjunc = vertcat(foundlistjunc, 1);
    else
        foundlistjunc = vertcat(foundlistjunc, 0);
    end 

%     junczscore = single(load_untouch_nii('MNI_MRFT2juncz.nii').img);
    junczscore = single(load_untouch_nii('MNI_T2z.nii').img);
    junczscore(isnan(junczscore)) = 0; 
    junczscore(junczscore<thresh) = 0;
    junczscore(junczscore>10) = 0;
    junczscore(junczscore>=thresh) = 1;

%     pixelValues = junczscore(:);
%     thresholdValue = prctile(pixelValues, threshper);
%     junczscore = junczscore > thresholdValue;
%     copyfile MNI_MRFT2juncz.nii MNI_T2z_binary.nii

    file = load_untouch_nii('MNI_T2z_binary.nii');
%     file.img = junczscore;
    file.img = medfilt3(junczscore);
    save_untouch_nii(file, 'MNI_T2z_binary.nii')


    overlap = file.img & roii;
    overlap(isnan(overlap)) = 0; 
    if sum(overlap(:)) > voxelthresh
        foundlistT2 = vertcat(foundlistT2, 1);
    else
        foundlistT2 = vertcat(foundlistT2, 0);
    end 

    extzscore = single(load_untouch_nii('MNI_MRFextenz.nii').img);
    extzscore(isnan(extzscore)) = 0;
    extzscore(extzscore<thresh) = 0;
    extzscore(extzscore>8) = 0;
    extzscore(extzscore>=thresh) = 1;

%     pixelValues = extzscore(:);
%     thresholdValue = prctile(pixelValues, threshper);
%     extzscore = extzscore > thresholdValue;
 
    overlap = extzscore & roii;
    overlap(isnan(overlap)) = 0; 
    if sum(overlap(:)) > voxelthresh
        foundlistext = vertcat(foundlistext, 1);
    else
        foundlistext = vertcat(foundlistext, 0);
    end 

end
either = ((foundlistjunc == 1) | (foundlistT1 == 1) | (foundlistT2 == 1) | (foundlistext == 1));
summary = horzcat(foundlistjunc, foundlistT1,  foundlistT2, foundlistext, either);

%% clipping
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
subjID=[(studyID)];
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps

for p = subjID
    p
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
   
    t1wz = load_untouch_nii('MNI_MRFjuncz.nii');
    t1wzi = double(t1wz.img);
    t1wzi(t1wzi>10) = 10;
    t1wzi(t1wzi<-10) = -10;
    t1wz.img = t1wzi;
    save_untouch_nii(t1wz, 'MNI_MRFjuncz.nii');


    extz = load_untouch_nii('MNI_MRFextenz.nii');
    extzi = double(extz.img);
    extzi(extzi>10) = 10;
    extzi(extzi<-10) = -10;
    extz.img = extzi;
    save_untouch_nii(extz, 'MNI_MRFextenz.nii');

end

%% clipping
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
subjID=[(studyID)];

for p = subjID
    p
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    copyfile MNI_WM_prob.nii MNI_CSF_prob.nii
    t1wz = load_untouch_nii('MNI_T1w.nii');
    t1wzi = double(t1wz.img);
    t1wzi(t1wzi>0) = 1;

    GMz = load_untouch_nii('MNI_GM_prob.nii');
    GMzi = double(GMz.img);
    WMz = load_untouch_nii('MNI_WM_prob.nii');
    WMzi = double(WMz.img);

    extz = load_untouch_nii('MNI_CSF_prob.nii');
    extzi = double(extz.img);
    extzi = t1wzi - GMzi - WMzi;
    t1wzi(t1wzi<0) = 0;
    extz.img = extzi;
    save_untouch_nii(extz, 'MNI_CSF_prob.nii');

end