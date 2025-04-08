subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
mT1 = zeros(182,218,182);
sdT1 = mT1;
mT2 = mT1;
sdT2 = mT1;
mGM = mT1;
sdGM = mT1;
mWM = mT1;
sdWM = mT1;
mCSF = mT1;
sdCSF = mT1;
mjunc = mT1;
sdjunc = mT1;

% t1max = 2000; 
% t2max = 180;
% t1min = 600;
% t2min = 20;
t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;
% n = 1;
n_matrix = zeros(182,218,182);
n_matrix_p = zeros(182,218,182);
n_matrix_p2 = zeros(182,218,182);
for p = subjID
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path);

%            f1 = char(append('n_syN_T1_', strcat(p), '_brain.nii'));
%            f2 = char(append('n_syN_T2_', strcat(p), '_brain.nii'));
%            f3 = char(append('n_syN_GM_', strcat(p), '_brain.nii'));
%            f4 = char(append('n_syN_WM_', strcat(p), '_brain.nii'));
    f1 = 'MNI_T1.nii';
    f2 = 'MNI_T2.nii';
    f3 = 'MNI_GM_prob.nii';
    f4 = 'MNI_WM_prob.nii';
    f5 = 'MNI_CSF_prob.nii';
    
    a = niftiread(f1);
    b = niftiread(f2);
    c = niftiread(f3);
    d = niftiread(f4);
    e = niftiread(f5);

    f6 = 'MNI_motion_mask.nii';
    f7 = 'MNI_MRFbinary.nii';

    f = niftiread(f6);
    g = niftiread(f7);

%     e(a==0) = 0;
%     e((a>t1max)|(a<t1min)|(b>t2max)|(b<t2min)) = 0;
%     e(g>0.05) = 0;
%     d(b==0) = 0;
%     d(c==0) = 0;

    n_matrix_p2 = n_matrix_p;
    n_matrix_p = n_matrix;
    n_matrix = n_matrix+f;

    oma = mT1;
    omb = mT2;
    omc = mGM;
    omd = mWM;
    ome = mCSF;
    omg = mjunc;

    mT1 = (mT1.*(n_matrix_p) + a.*f) ./ n_matrix;
    mT2 = (mT2.*(n_matrix_p) + b.*f) ./ n_matrix;
    mGM = (mGM.*(n_matrix_p) + c.*f) ./ n_matrix;
    mWM = (mWM.*(n_matrix_p) + d.*f) ./ n_matrix;
    mCSF = (mCSF.*(n_matrix_p) + e.*f) ./ n_matrix;
    mjunc = (mjunc.*(n_matrix_p) + g.*f) ./ n_matrix;

    if (p~="V04_12578")
        sdT1 = sqrt(((n_matrix_p2).*sdT1.^2  +  (a-mT1).*(a-oma))./(n_matrix_p));
        sdT2 = sqrt(((n_matrix_p2).*sdT2.^2  +  (b-mT2).*(b-omb))./(n_matrix_p));
        sdGM = sqrt(((n_matrix_p2).*sdGM.^2  +  (c-mGM).*(c-omc))./(n_matrix_p));
        sdWM = sqrt(((n_matrix_p2).*sdWM.^2  +  (d-mWM).*(d-omd))./(n_matrix_p));
        sdCSF = sqrt(((n_matrix_p2).*sdCSF.^2  +  (e-mCSF).*(e-ome))./(n_matrix_p));
        sdjunc = sqrt(((n_matrix_p2).*sdjunc.^2  +  (g-mjunc).*(g-omg))./(n_matrix_p));
    end
end

cd Z:\Imaging\Multimodal\MRF\Peter
save('MRFvalues_withcsf.mat','mT1','mT2','sdT1','sdT2', 'mGM', 'sdGM', 'mWM', 'sdWM', 'mCSF', 'sdCSF');

cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps

file = load_untouch_nii('T1_mean.nii');
file.img = mT1.*atlasi;
save_untouch_nii(file, 'T1_mean.nii');

file = load_untouch_nii('T2_mean.nii');
file.img = mT2.*atlasi;
save_untouch_nii(file, 'T2_mean.nii');

file = load_untouch_nii('T1_std.nii');
file.img = sdT1.*atlasi;
save_untouch_nii(file, 'T1_std.nii');

file = load_untouch_nii('T2_std.nii');
file.img = sdT2.*atlasi;
save_untouch_nii(file, 'T2_std.nii');

file = load_untouch_nii('GM_mean.nii');
file.img = mGM.*atlasi;
save_untouch_nii(file, 'GM_mean.nii');

file = load_untouch_nii('WM_mean.nii');
file.img = mWM.*atlasi;
save_untouch_nii(file, 'WM_mean.nii');

file = load_untouch_nii('CSF_mean.nii');
file.img = mCSF.*atlasi;
save_untouch_nii(file, 'CSF_mean.nii');

file = load_untouch_nii('GM_std.nii');
file.img = sdGM.*atlasi;
save_untouch_nii(file, 'GM_std.nii');

file = load_untouch_nii('WM_std.nii');
file.img = sdWM.*atlasi;
save_untouch_nii(file, 'WM_std.nii');

file = load_untouch_nii('CSF_std.nii');
file.img = sdCSF.*atlasi;
save_untouch_nii(file, 'CSF_std.nii');

file = load_untouch_nii('junction_mean.nii');
file.img = mjunc.*atlasi;
save_untouch_nii(file, 'junction_std.nii');

file = load_untouch_nii('junction_std.nii');
file.img = sdjunc.*atlasi;
save_untouch_nii(file, 'junction_std.nii');

% cd Z:\Imaging\Multimodal\MRF\Peter
% save('meanMRFvalues_nocsf.mat','mT1','mT2','sdT1','sdT2','n_matrix');
%
% cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps
% file = load_untouch_nii('T1_mean_nocsf.nii');
% file.img = zeros(182,218,182);
% file.img = mT1;
% save_untouch_nii(file, 'T1_mean_nocsf.nii');
% 
% file = load_untouch_nii('T2_mean_nocsf.nii');
% file.img = zeros(182,218,182);
% file.img = mT2;
% save_untouch_nii(file, 'T2_mean_nocsf.nii');
% 
% file = load_untouch_nii('T1_std_nocsf.nii');
% file.img = zeros(182,218,182);
% file.img = sdT1;
% save_untouch_nii(file, 'T1_std_nocsf.nii');
% 
% file = load_untouch_nii('T2_std_nocsf.nii');
% file.img = zeros(182,218,182);
% file.img = sdT2;
% save_untouch_nii(file, 'T2_std_nocsf.nii');
%% includes CSF
subjID=[(studyID)];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
mT1 = zeros(182,218,182);
sdT1 = mT1;
mT2 = mT1;
sdT2 = mT1;
mGM = mT1;
sdGM = mT1;
mWM = mT1;
sdWM = mT1;
mCSF = mT1;
sdCSF = mT1;


% t1max = 2000; 
% t2max = 180;
% t1min = 600;
% t2min = 20;
t1max = 1759;
t2max = 106;
t1min = 743;
t2min = 29;

n = 1;
n_matrix = zeros(182,218,182);
for p = subjID
    path = strcat(MRF_path,'\',p);
    cd(path);
    mkdir MRFz
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path);

%            f1 = char(append('n_syN_T1_', strcat(p), '_brain.nii'));
%            f2 = char(append('n_syN_T2_', strcat(p), '_brain.nii'));
%            f3 = char(append('n_syN_GM_', strcat(p), '_brain.nii'));
%            f4 = char(append('n_syN_WM_', strcat(p), '_brain.nii'));
       f1 = 'MNI_T1.nii';
       f2 = 'MNI_T2.nii';
       f3 = 'MNI_GM_prob.nii';
       f4 = 'MNI_WM_prob.nii';
       f5 = 'MNI_CSF_prob.nii';

       a = niftiread(f1);
       b = niftiread(f2);
       c = niftiread(f3);
       d = niftiread(f4);
       e = niftiread(f5);

       oma = mT1;
       omb = mT2;
       omc = mGM;
       omd = mWM;
       ome = mCSF;

       mT1 = (mT1.*(n-1) + a) ./ n;
       mT2 = (mT2.*(n-1) + b) ./ n;
       mGM = (mGM.*(n-1) + c) ./ n;
       mWM = (mWM.*(n-1) + d) ./ n;
       mCSF = (mCSF.*(n-1) + e) ./ n;

       if n > 1 
%             sdT1 = sqrt(((n-2)*sdT1.^2+(n-1)/n*(mT1.^2 + a.^2 - 2.*mT1.*a))./(n-1));
%             sdT2 = sqrt(((n-2)*sdT2.^2+(n-1)/n*(mT2.^2 + b.^2 - 2.*mT2.*b))./(n-1));
            sdT1 = sqrt(((n-2).*sdT1.^2  +  (a-mT1).*(a-oma))./(n-1));
            sdT2 = sqrt(((n-2).*sdT2.^2  +  (b-mT2).*(b-omb))./(n-1));
            sdGM = sqrt(((n-2).*sdGM.^2  +  (c-mGM).*(c-omc))./(n-1));
            sdWM = sqrt(((n-2).*sdWM.^2  +  (d-mWM).*(d-omd))./(n-1));
            sdCSF = sqrt(((n-2).*sdCSF.^2  +  (e-mCSF).*(e-ome))./(n-1));
       end
%             sdT1 = sqrt(((n-2).*sdT1.^2  +  (a-mT1).*(a-oma))./(n-1));
%             sdT2 = sqrt(((n-2).*sdT2.^2  +  (b-mT2).*(b-omb))./(n-1));

%        sdGM(good2) = sqrt( ((n_matrix(good2)-2).*sdGM(good2).^2  +  (c(good2)-mGM(good2)).*(c(good2)-omc(good2)) )./(n_matrix(good2)-1));
%        sdWM(good2) = sqrt( ((n_matrix(good2)-2).*sdWM(good2).^2  +  (d(good2)-mWM(good2)).*(d(good2)-omd(good2)) )./(n_matrix(good2)-1));     
       n = n + 1
end
n_matrix = ones(182,218,182)*length(subjID);
cd Z:\Imaging\Multimodal\MRF\Peter
save('MRFvalues_withcsf.mat','mT1','mT2','sdT1','sdT2', 'mGM', 'sdGM', 'mWM', 'sdWM', 'mCSF', 'sdCSF');
cd Z:\Imaging\Multimodal\MRF\Peter\MRFzmaps

file = load_untouch_nii('T1_mean.nii');
file.img = mT1;
save_untouch_nii(file, 'T1_mean.nii');

file = load_untouch_nii('T2_mean.nii');
file.img = mT2;
save_untouch_nii(file, 'T2_mean.nii');

file = load_untouch_nii('T1_std.nii');
file.img = sdT1;
save_untouch_nii(file, 'T1_std.nii');

file = load_untouch_nii('T2_std.nii');
file.img = sdT2;
save_untouch_nii(file, 'T2_std.nii');

file = load_untouch_nii('GM_mean.nii');
file.img = mGM;
save_untouch_nii(file, 'GM_mean.nii');

file = load_untouch_nii('WM_mean.nii');
file.img = mWM;
save_untouch_nii(file, 'WM_mean.nii');

file = load_untouch_nii('CSF_mean.nii');
file.img = mCSF;
save_untouch_nii(file, 'CSF_mean.nii');

file = load_untouch_nii('GM_std.nii');
file.img = sdGM;
save_untouch_nii(file, 'GM_std.nii');

file = load_untouch_nii('WM_std.nii');
file.img = sdWM;
save_untouch_nii(file, 'WM_std.nii');

file = load_untouch_nii('CSF_std.nii');
file.img = sdCSF;
save_untouch_nii(file, 'CSF_std.nii');
%% Linear regression based MRF mean and var
clear all
subjID=[(studyID)];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Normal';
mT1oT2_NMCC = zeros(182,218,182,length(subjID));
mT1 = zeros(182,218,182,length(subjID));
mT2 = zeros(182,218,182,length(subjID));

b0T1oT2 = zeros(182,218,182);
b0T1 = zeros(182,218,182);
b0T2 = zeros(182,218,182);

b1T1oT2 = zeros(182,218,182);
b1T1 = zeros(182,218,182);
b1T2 = zeros(182,218,182);

b2T1oT2 = zeros(182,218,182);
b2T1 = zeros(182,218,182);
b2T2 = zeros(182,218,182);

RMSET1oT2 = zeros(182,218,182);
RMSET1 = zeros(182,218,182);
RMSET2 = zeros(182,218,182);

cd('Z:\Imaging\Multimodal\MRF\Peter')
clinifo = readmatrix('HCinfo.csv', 'Range','1:69');
age = clinifo(2:length(subjID)+1,1);
sex = clinifo(2:length(subjID)+1,3) - 1;

Y1 = zeros(length(subjID),1);
Y2 = zeros(length(subjID),1);
for p = 1:length(subjID)
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    cd(path);
    f1 = 'MNI_T1.nii';
    f2 = 'MNI_T2.nii';
%     f3 = 'MNI_T1oT2_nmCC.nii'
    a = niftiread(f1);
    b = niftiread(f2);
end
for i = 1:182
    i
    datetime
    for j = 1:218
        for k = 1:182
            Y1 = mT1(i,j,k,:);
            Y2 = mT2(i,j,k,:);
%             Y3 = mT1oT2_NMCC(i,j,k,:);
            mdl1 = fitlm([age, sex], Y1);
            mdl2 = fitlm([age, sex], Y2);
%             mdl3 = fitlm([age, sex], Y3);

%             b0T1oT2(i,j,k) = mdl3.Coefficients.Estimate(1);
%             b1T1oT2(i,j,k) = mdl3.Coefficients.Estimate(2);
%             b2T1oT2(i,j,k) = mdl3.Coefficients.Estimate(3);

            b0T1(i,j,k) = mdl1.Coefficients.Estimate(1);
            b1T1(i,j,k) = mdl1.Coefficients.Estimate(2);
            b2T1(i,j,k) = mdl1.Coefficients.Estimate(3);

            b0T2(i,j,k) = mdl2.Coefficients.Estimate(1);
            b1T2(i,j,k) = mdl2.Coefficients.Estimate(2);
            b2T2(i,j,k) = mdl2.Coefficients.Estimate(3);  
        end
    end
end

% cd T:\Imaging\Multimodal\MRF\Peter
% save('MRFvalues_linreg.mat','mT1','mT2','sdT1','sdT2');
% save('MRFvalues_linreg.mat', 'b0T1', 'b1T1', 'b2T1', ...
%                              'b0T2', 'b1T2', 'b2T2');
cd('Z:\Imaging\Multimodal\MRF\Peter')
save('MRFvalues_linreg_T1.mat', 'b0T1', 'b1T1', 'b2T1', 'RMSET1');

%% Linear regression based MRF mean and separate T1 and T2
clear all
subjID=[(studyID)];

% mT1oT2_NMCC = zeros(182,218,182,length(subjID));
mT1 = zeros(182,218,182,length(subjID));
mT2 = zeros(182,218,182,length(subjID));

% b0T1oT2 = zeros(182,218,182);
b0T1 = zeros(182,218,182);
b0T2 = zeros(182,218,182);

% b1T1oT2 = zeros(182,218,182);
b1T1 = zeros(182,218,182);
b1T2 = zeros(182,218,182);

% b2T1oT2 = zeros(182,218,182);
b2T1 = zeros(182,218,182);
b2T2 = zeros(182,218,182);

% RMSET1oT2 = zeros(182,218,182);
RMSET1 = zeros(182,218,182);
RMSET2 = zeros(182,218,182);

cd('Z:\Imaging\Multimodal\MRF\Peter')
clinifo = readmatrix('HCinfo.csv', 'Range','1:68');
age = clinifo(2:length(subjID)+1,1);
sex = clinifo(2:length(subjID)+1,3) - 1;

% Y1 = zeros(length(subjID),1);
% Y2 = zeros(length(subjID),1);
for p = 1:length(subjID)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    cd(path);
    f1 = 'MNI_T1.nii';
    f2 = 'MNI_T2.nii';
%     f1 = 'sMNI_T1_nocsf.nii';
%     f2 = 'sMNI_T2_nocsf.nii';
%     f3 = 'MNI_T1oT2_nmCC.nii'
    a = niftiread(f1);
    b = niftiread(f2);
    mT1(:,:,:,p) = a;
    mT2(:,:,:,p) = b;
end

for i = 1:182
    i
    for j = 1:218
        for k = 1:182
            Y1 = reshape(mT1(i,j,k,:), [], 1);
            Y2 = reshape(mT2(i,j,k,:), [], 1);
%             Y3 = mT1oT2_NMCC(i,j,k,:);
            if (sum(Y1) ~= 0)&(mean(Y1)<2900)
                mdl1 = fitlm([age, sex], Y1);
                mdl2 = fitlm([age, sex], Y2);
    
                b0T1(i,j,k) = mdl1.Coefficients.Estimate(1);
                b1T1(i,j,k) = mdl1.Coefficients.Estimate(2);
                b2T1(i,j,k) = mdl1.Coefficients.Estimate(3);
                RMSET1(i,j,k) = mdl1.RMSE;
                b0T2(i,j,k) = mdl2.Coefficients.Estimate(1);
                b1T2(i,j,k) = mdl2.Coefficients.Estimate(2);
                b2T2(i,j,k) = mdl2.Coefficients.Estimate(3);  
                RMSET2(i,j,k) = mdl2.RMSE;
            end
        end
    end
end

cd('Z:\Imaging\Multimodal\MRF\Peter')
% save('MRFvalues_linreg.mat', 'b0T1', 'b1T1', 'b2T1', ...
%                              'b0T2', 'b1T2', 'b2T2');
save('MRFvalues_linreg_complete.mat', 'b0T1', 'b1T1', 'b2T1', 'RMSET1',...
                                'b0T2', 'b1T2', 'b2T2', 'RMSET2');

%% patient predict and score
subjID=[(studyID)];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
cd('Z:\Imaging\Multimodal\MRF\Recon_MRF_3T');
atlasi = single(load_untouch_nii('MNI152_T1_1mm_brain.nii').img);
atlasi(atlasi>0) = 1;
cd('Z:\Imaging\Multimodal\MRF\Peter');
ventmask= single(load_untouch_nii('ventricular_mask.nii').img);
atlasi(ventmask==1) = 0;
cd('Z:\Imaging\Multimodal\MRF\Peter')
load MRFvalues_linreg_complete.mat
% load MRFvalues_linreg_smooth.mat
clinifo = readmatrix('Ptinfo.xlsx', 'Range','1:44');
age = clinifo(2:length(subjID)+1,1);
sex = clinifo(2:length(subjID)+1,2) - 1;

for p = 1:length(subjID)
    p
    path = strcat(MRF_path,'\',subjID(p),'\MRF_VBM');
    cd(path);
    copyfile MNI_T1z.nii MNI_T1LMz.nii
    copyfile MNI_T1z.nii MNI_T2LMz.nii

    f1 = 'MNI_T1.nii';
    f2 = 'MNI_T2.nii';
%     f1 = 'sMNI_T1_nocsf.nii';
%     f2 = 'sMNI_T2_nocsf.nii';

    a = niftiread(f1);
    b = niftiread(f2);

    t1 = b0T1 + age(p)*b1T1 + sex(p)*b2T1;
    t2 = b0T2 + age(p)*b1T2 + sex(p)*b2T2;

    zt1 = (a - t1)./RMSET1;
    zt2 = (b - t2)./RMSET2;
    zt1(isnan(zt1))=0;
    zt1(isinf(zt1))=0;
%     zt1(zt1>7) = 0;
    zt1(zt1>10) = 10;
    zt1(zt1<-10) = -10;
    zt1(a>1800) = 0;
    zt1(b>110) = 0;

    zt2(isnan(zt2))=0;
    zt2(isinf(zt2))=0;
%     zt2(zt2>7) = 0;
%     zt2(zt2<3) = 0;
    zt2(zt2>10) = 10;
    zt2(zt2<-10) = -10;
    zt2(a>1800) = 0;
    zt2(b>110) = 0;
%     zt1(1:10,1:10,90)
%     zt1(90:100,90:100,90)

%     zt2(isnan(zt2))=0;

%     zt2(isinf(zt2))=0;

    imt1 = load_untouch_nii('MNI_T1LMz.nii');
    imt1.img = double(zt1).*atlasi;
    save_untouch_nii(imt1, 'MNI_T1LMz.nii');

    imt2 = load_untouch_nii('MNI_T2LMz.nii');
    imt2.img = double(zt2).*atlasi;
    save_untouch_nii(imt2, 'MNI_T2LMz.nii');

end
%%
% MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
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

    file = load_untouch_nii('MNI_T1LMz.nii');
    file.img(file.img<3) = 0;
    file.img(file.img>8) = 0;
    t1zscore = file.img;
    t1zscore(t1zscore>0) = 1;
    
    cc = bwconncomp(t1zscore, 6); % Using 26-connectivity for 3D % Step 2: Measure the size of each component and filter 
    componentSizes = cellfun(@numel, cc.PixelIdxList); 
    minVoxelSize = 5; 
    componentsToKeep = find(componentSizes >= minVoxelSize); 
    filteredImg = false(size(t1zscore)); 
    for i = 1:length(componentsToKeep)
        filteredImg(cc.PixelIdxList{componentsToKeep(i)}) = true; 
    end 
    t1zscore = filteredImg;
    
    %     t1zscore = medfilt3(t1zscore);
    file = load_untouch_nii('MNI_T1z_binary.nii');
    file.img = t1zscore;

    save_untouch_nii(file, 'MNI_T1z_binary.nii');

    roii = single(load_untouch_nii('MNI_ROI_final.nii').img);
    roii(isnan(roii)) = 0; 
    overlap = file.img & roii;
    
    if sum(overlap(:)) > voxelthresh
        foundlistT1 = vertcat(foundlistT1, 1);
    else
        foundlistT1 = vertcat(foundlistT1, 0);
    end


end
either = (foundlistT1 == 1);
% summary = horzcat(foundlistjunc, foundlistT1,  foundlistT2, foundlistext, either);