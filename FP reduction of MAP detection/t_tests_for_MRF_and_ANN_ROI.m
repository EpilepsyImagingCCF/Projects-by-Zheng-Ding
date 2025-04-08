%% all promising patients and separate into GM and WM
subjID=["PXX_XXXXX"];
subtypes = ["mMCD" "mMCD" "IIB" "mMCD"...
    "IIB" "IIB" "mMCD" "mMCD" "IIA"...
    "IIA" "IIA" "IIB" "IIA" "IIB"...
    "mMCD" "mMCD" "MOGHE" "IIB" "MOGHE"...
    "IIA" "IIB" "IIB" "IIB" "IIB"...
    "IIA" "IIA" "mMCD" "IIA" "IIB"...
    "IIB" "mMCD" "IIA" "IIB" "mMCD"...
    "mMCD" "IIA" "mMCD" "IIB" "mMCD"];

% MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
addpath('/home/dingz/beegfs/code_and_scripts');
MRF_path='/home/dingz/beegfs/all_patients';

t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;

for volthresh = 1
volthresh
xtallW = [];
ytallW = [];
xfallW = [];
yfallW = [];
xtallG = [];
ytallG = [];
xfallG = [];
yfallG = [];
xtallWV = [];
ytallWV = [];
xfallWV = [];
yfallWV = [];
xtallGV = [];
ytallGV = [];
xfallGV = [];
yfallGV = [];

for p = subjID
    
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    path
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    if p == "P31_13375"
        anni(:,242:end,1:57) = 0;
    end
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    mask(mask>0.95) = 1;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img).*mask;
    GM(GM<0.6) = 0;
    GM(GM>=0.6) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img).*mask;
    WM(WM<0.6) = 0;
    WM(WM>=0.6) = 1;

    t1s = load_untouch_nii('MNI_T1.nii');
    t1si = single(t1s.img);
    t2s = load_untouch_nii('MNI_T2.nii');
    t2si = single(t2s.img);

    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>volthresh)];
    numberOfBlobs = size(allinds,1);

    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    if p == "P83_XXXXX"
        roii(:,269:end,:) = 0;
    end

    [roiblob, nBlobs] = bwlabeln(roii);
    roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
    roiinds = roiMeasurements.VoxelIdxList{1};
    
    for i = 1:numberOfBlobs
        xtGM = [];
        ytGM = [];
        xfGM = [];
        yfGM = [];
        xtWM = [];
        ytWM = [];
        xfWM = [];
        yfWM = [];
        A = blobMeasurements.VoxelIdxList{i};
        if sum(ismember(roiinds,A))>0
            xtGM = vertcat(xtGM, t1si(allinds{i}).*GM(allinds{i}));
            ytGM = vertcat(ytGM, t2si(allinds{i}).*GM(allinds{i}));
            xtGMc = xtGM;
            ytGMc = ytGM;
            xtGM((xtGMc > t1max)|(xtGMc< t1min)|(ytGMc > t2max)|(ytGMc< t2min)) = []; 
            ytGM((ytGMc > t2max)|(ytGMc< t2min)|(xtGMc > t1max)|(xtGMc< t1min)) = [];

            xtWM = vertcat(xtWM, t1si(allinds{i}).*WM(allinds{i}));
            ytWM = vertcat(ytWM, t2si(allinds{i}).*WM(allinds{i}));
            xtWMc = xtWM;
            ytWMc = ytWM;
            xtWM((xtWMc > t1max)|(xtWMc< t1min)|(ytWMc > t2max)|(ytWMc< t2min)) = []; 
            ytWM((ytWMc > t2max)|(ytWMc< t2min)|(xtWMc > t1max)|(xtWMc< t1min)) = [];

        else
            xfGM = vertcat(xfGM, t1si(allinds{i}).*GM(allinds{i}));
            yfGM = vertcat(yfGM, t2si(allinds{i}).*GM(allinds{i}));
            xfGMc = xfGM;
            yfGMc = yfGM;
            xfGM((xfGMc > t1max)|(xfGMc< t1min)|(yfGMc > t2max)|(yfGMc< t2min)) = []; 
            yfGM((yfGMc > t2max)|(yfGMc< t2min)|(xfGMc > t1max)|(xfGMc< t1min)) = [];
            
            xfWM = vertcat(xfWM, t1si(allinds{i}).*WM(allinds{i}));
            yfWM = vertcat(yfWM, t2si(allinds{i}).*WM(allinds{i}));
            xfWMc = xfWM;
            yfWMc = yfWM;
            xfWM((xfWMc > t1max)|(xfWMc< t1min)|(yfWMc > t2max)|(yfWMc< t2min)) = []; 
            yfWM((yfWMc > t2max)|(yfWMc< t2min)|(xfWMc > t1max)|(xfWMc< t1min)) = [];
        end
        xtallG = vertcat(xtallG, mean(xtGM));
        ytallG = vertcat(ytallG, mean(ytGM));
        xfallG = vertcat(xfallG, mean(xfGM));
        yfallG = vertcat(yfallG, mean(yfGM));
        xtallGV = vertcat(xtallGV, length(xtGM));
        ytallGV = vertcat(ytallGV, length(ytGM));
        xfallGV = vertcat(xfallGV, length(xfGM));
        yfallGV = vertcat(yfallGV, length(yfGM));
        
        xtallW = vertcat(xtallW, mean(xtWM));
        ytallW = vertcat(ytallW, mean(ytWM));
        xfallW = vertcat(xfallW, mean(xfWM));
        yfallW = vertcat(yfallW, mean(yfWM));
        xtallWV = vertcat(xtallWV, length(xtWM));
        ytallWV = vertcat(ytallWV, length(ytWM));
        xfallWV = vertcat(xfallWV, length(xfWM));
        yfallWV = vertcat(yfallWV, length(yfWM));
    end 
end

% calculate mean and std 
xtallW = rmmissing(xtallW);
ytallW = rmmissing(ytallW);
xfallW = rmmissing(xfallW);
yfallW = rmmissing(yfallW);
xtallG = rmmissing(xtallG);
ytallG = rmmissing(ytallG);
xfallG = rmmissing(xfallG);
yfallG = rmmissing(yfallG);
xtallWV = nonzeros(xtallWV);
ytallWV = nonzeros(ytallWV);
xfallWV = nonzeros(xfallWV);
yfallWV = nonzeros(yfallWV);
xtallGV = nonzeros(xtallGV);
ytallGV = nonzeros(ytallGV);
xfallGV = nonzeros(xfallGV);
yfallGV = nonzeros(yfallGV);
%white matter
meanwt1t = round(mean(xtallW));
stdwt1t = round(std(xtallW));
meanwt2t = round(mean(ytallW),1);
stdwt2t = round(std(ytallW),1);
meanwt1f = round(mean(xfallW));
stdwt1f = round(std(xfallW));
meanwt2f = round(mean(yfallW),1);
stdwt2f = round(std(yfallW),1);
[h1,p1w] = ttest2(xtallW, xfallW,'Vartype','unequal');
[h2,p2w] = ttest2(ytallW, yfallW,'Vartype','unequal');
%gray matter
meangt1t = round(mean(xtallG));
stdgt1t = round(std(xtallG));
meangt2t = round(mean(ytallG),1);
stdgt2t = round(std(ytallG),1);
meangt1f = round(mean(xfallG));
stdgt1f = round(std(xfallG));
meangt2f = round(mean(yfallG),1);
stdgt2f = round(std(yfallG),1);
[h1,p1g] = ttest2(xtallG, xfallG, 'Vartype','unequal');
[h2,p2g] = ttest2(ytallG, yfallG, 'Vartype','unequal');
% make tables
GM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meangt1t),string(char(177)),string(stdgt1t)); strcat(string(meangt1f),string(char(177)),string(stdgt1f));string(p1g)];
T2 = [strcat(string(meangt2t),string(char(177)),string(stdgt2t)); strcat(string(meangt2f),string(char(177)),string(stdgt2f));string(p2g)];
clust = [length(xtallGV);length(xfallGV);length(xtallGV)+length(xfallGV)]
GMtab1 = table(GM,T1,T2,clust)

WM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meanwt1t),string(char(177)),string(stdwt1t)); strcat(string(meanwt1f),string(char(177)),string(stdwt1f));string(p1w)];
T2 = [strcat(string(meanwt2t),string(char(177)),string(stdwt2t)); strcat(string(meanwt2f),string(char(177)),string(stdwt2f));string(p2w)];
clust = [length(xtallWV);length(xfallWV);length(xtallWV)+length(xfallWV)]
WMtab1 = table(WM,T1,T2,clust)

filename = 'ttests.xlsx';
writetable(GMtab1,filename,'Sheet',1)
writetable(WMtab1,filename,'Sheet',2)

% All promising patients and separate into GM and WM normalized

xtallW = [];
ytallW = [];
xfallW = [];
yfallW = [];
xtallG = [];
ytallG = [];
xfallG = [];
yfallG = [];

for p = subjID
    path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;
    mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
    mask(mask>0.95) = 1;
    mask(mask<0.95) = 0;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img).*mask;
    GM(GM<0.6) = 0;
    GM(GM>=0.6) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img).*mask;
    WM(WM<0.6) = 0;
    WM(WM>=0.6) = 1;
    t1s = load_untouch_nii('MNI_T1.nii');
    t1si = single(t1s.img);
    t2s = load_untouch_nii('MNI_T2.nii');
    t2si = single(t2s.img);
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>volthresh)];
    numberOfBlobs = size(allinds,1);
    
    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    [roiblob, nBlobs] = bwlabeln(roii);
    roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
    roiinds = roiMeasurements.VoxelIdxList{1};
    
%     region = zeros(182, 218, 182);
%     region(allinds{i}) = 1; 
%     [mT1, sdT1, mT2, sdT2] = MRFclusterwisezscore(region);
    
    for i = 1:numberOfBlobs
        xtGM = [];
        ytGM = [];
        xfGM = [];
        yfGM = [];
        xtWM = [];
        ytWM = [];
        xfWM = [];
        yfWM = [];
        A = blobMeasurements.VoxelIdxList{i};
        region = zeros(182, 218, 182);
        region(allinds{i}) = 1; 
        [mT1, sdT1, mT2, sdT2, mGMT1, sdGMT1, mWMT1, sdWMT1, mGMT2, sdGMT2, mWMT2, sdWMT2] = MRFclusterwisezscore(region);

        if sum(ismember(roiinds,A))>0
            xtGM = vertcat(xtGM, t1si(allinds{i}).*GM(allinds{i}));
            ytGM = vertcat(ytGM, t2si(allinds{i}).*GM(allinds{i}));
            xtGMc = xtGM;
            ytGMc = ytGM;
            xtGM((xtGMc > t1max)|(xtGMc< t1min)|(ytGMc > t2max)|(ytGMc< t2min)) = []; 
            ytGM((ytGMc > t2max)|(ytGMc< t2min)|(xtGMc > t1max)|(xtGMc< t1min)) = [];    
            xtallG = vertcat(xtallG, (mean(xtGM) - mT1)/sdT1);
            ytallG = vertcat(ytallG, (mean(ytGM) - mT2)/sdT2);
            xtWM = vertcat(xtWM, t1si(allinds{i}).*WM(allinds{i}));
            ytWM = vertcat(ytWM, t2si(allinds{i}).*WM(allinds{i}));
            xtWMc = xtWM;
            ytWMc = ytWM;
            xtWM((xtWMc > t1max)|(xtWMc< t1min)|(ytWMc > t2max)|(ytWMc< t2min)) = []; 
            ytWM((ytWMc > t2max)|(ytWMc< t2min)|(xtWMc > t1max)|(xtWMc< t1min)) = [];       
            xtallW = vertcat(xtallW, (mean(xtWM) - mT1)/sdT1);
            ytallW = vertcat(ytallW, (mean(ytWM) - mT2)/sdT2);
        else
            xfGM = vertcat(xfGM, t1si(allinds{i}).*GM(allinds{i}));
            yfGM = vertcat(yfGM, t2si(allinds{i}).*GM(allinds{i}));
            xfGMc = xfGM;
            yfGMc = yfGM;
            xfGM((xfGMc > t1max)|(xfGMc< t1min)|(yfGMc > t2max)|(yfGMc< t2min)) = []; 
            yfGM((yfGMc > t2max)|(yfGMc< t2min)|(xfGMc > t1max)|(xfGMc< t1min)) = [];
            xfallG = vertcat(xfallG, (mean(xfGM) - mT1)/sdT1);
            yfallG = vertcat(yfallG, (mean(yfGM) - mT1)/sdT2);
            xfWM = vertcat(xfWM, t1si(allinds{i}).*WM(allinds{i}));
            yfWM = vertcat(yfWM, t2si(allinds{i}).*WM(allinds{i}));
            xfWMc = xfWM;
            yfWMc = yfWM;
            xfWM((xfWMc > t1max)|(xfWMc< t1min)|(yfWMc > t2max)|(yfWMc< t2min)) = []; 
            yfWM((yfWMc > t2max)|(yfWMc< t2min)|(xfWMc > t1max)|(xfWMc< t1min)) = [];
            xfallW = vertcat(xfallW, (mean(xfWM) - mT1)/sdT1);
            yfallW = vertcat(yfallW, (mean(yfWM) - mT2)/sdT2);
        end
    end
end


% calculate mean and std 
%white matter
xtallW = rmmissing(xtallW);
ytallW = rmmissing(ytallW);
xfallW = rmmissing(xfallW);
yfallW = rmmissing(yfallW);
xtallG = rmmissing(xtallG);
ytallG = rmmissing(ytallG);
xfallG = rmmissing(xfallG);
yfallG = rmmissing(yfallG);
meanwt1t = round(mean(xtallW),3);
stdwt1t = round(std(xtallW),3);
meanwt2t = round(mean(ytallW),3);
stdwt2t = round(std(ytallW),3);
meanwt1f = round(mean(xfallW),3);
stdwt1f = round(std(xfallW),3);
meanwt2f = round(mean(yfallW),3);
stdwt2f = round(std(yfallW),3);
[h1,p1w] = ttest2(xtallW, xfallW,'Vartype','unequal');
[h2,p2w] = ttest2(ytallW, yfallW,'Vartype','unequal');
%gray matter
meangt1t = round(mean(xtallG),3);
stdgt1t = round(std(xtallG),3);
meangt2t = round(mean(ytallG),3);
stdgt2t = round(std(ytallG),3);
meangt1f = round(mean(xfallG),3);
stdgt1f = round(std(xfallG),3);
meangt2f = round(mean(yfallG),3);
stdgt2f = round(std(yfallG),3);
[h1,p1g] = ttest2(xtallG, xfallG,'Vartype','unequal');
[h2,p2g] = ttest2(ytallG, yfallG,'Vartype','unequal');
% make tables
GM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meangt1t),string(char(177)),string(stdgt1t)); strcat(string(meangt1f),string(char(177)),string(stdgt1f));string(p1g)];
T2 = [strcat(string(meangt2t),string(char(177)),string(stdgt2t)); strcat(string(meangt2f),string(char(177)),string(stdgt2f));string(p2g)];
GMtab3 = table(GM,T1,T2)

WM = ["TP";"FP";"p-value"];
T1 = [strcat(string(meanwt1t),string(char(177)),string(stdwt1t)); strcat(string(meanwt1f),string(char(177)),string(stdwt1f));string(p1w)];
T2 = [strcat(string(meanwt2t),string(char(177)),string(stdwt2t)); strcat(string(meanwt2f),string(char(177)),string(stdwt2f));string(p2w)];
WMtab3 = table(WM,T1,T2)
end
filename = 'ttests.xlsx';
writetable(GMtab3,filename,'Sheet',3)
writetable(WMtab3,filename,'Sheet',4)