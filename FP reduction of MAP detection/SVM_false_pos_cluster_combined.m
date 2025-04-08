cd T:\Imaging\Multimodal\MRF\Peter
% load meanMRFvalues.mat
% load MRFstdev_voxel.mat
load MRFvalues_withcsf.mat

subjID=["PXX_XXXXX"];
subtypes = ["IIB" "mMCD" "IIB" "mMCD"...
    "IIB" "IIB" "mMCD" "mMCD" "IIA"...
    "IIA" "IIA" "IIB" "IIA" "IIB"...
    "mMCD" "mMCD" "MOGHE" "IIB" "MOGHE"...
    "IIA" "IIB" "IIB" "IIB" "IIB"...
    "IIA" "IIA" "mMCD" "IIA" "IIB"...
    "IIB" "mMCD" "IIA" "IIB" "mMCD"...
    "IIA" "IIB" "mMCD"];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
dict = readtable('atlas.xlsx');
atlas = load_untouch_nii('Talairach-labels-1mm.nii');
atlasi = single(atlas.img);

for volthresh = 1
xallW = [];
yallW = [];
xallG = [];
yallG = [];
xallWn = [];
yallWn = [];
xallGn = [];
yallGn = [];
xallWV = [];
yallWV = [];
xallGV = [];
yallGV = [];
clusttp =[];
ptids = [];
sutys = [];
lobes = [];
% gyri = [];
centers = [];

hastps = [];
t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;
% t1max = 1700; 
% t2max = 100;
% t1min = 750;
% t2min = 30;

cursor = 0;
for p = subjID
    p
    cursor = cursor + 1;
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;

    if p == "P31_XXXXX"
        anni(:,242:end,1:57) = 0;
    end
    mask = single(load_untouch_nii('MNI_GM.nii').img) + single(load_untouch_nii('MNI_WM.nii').img);
    mask(mask>0.95) = 1;
%     mask(mask<0.95) = 0;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img).*mask;
    GM(GM<0.8) = 0;
    GM(GM>=0.8) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img).*mask;
    WM(WM<0.8) = 0;
    WM(WM>=0.8) = 1;
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList', 'Centroid');
	allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>volthresh)];
    allcs = round(blobMeasurements.Centroid);
    allcs = allcs(blobMeasurements.Volume>volthresh,:);
    numberOfBlobs = size(allinds,1);


    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    if p == "P83_XXXXX"
        roii(:,269:end,:) = 0;
    end
    [roiblob, nBlobs] = bwlabeln(roii);
    roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
    roiinds = roiMeasurements.VoxelIdxList{1};
    t1s = load_untouch_nii('MNI_T1.nii');
    t1si = single(t1s.img);
    t2s = load_untouch_nii('MNI_T2.nii');
    t2si = single(t2s.img);

    for i = 1:numberOfBlobs
        xGM = [];
        yGM = [];
        xWM = [];
        yWM = [];
        xGMn = [];
        yGMn = [];
        xWMn = [];
        yWMn = [];
        xGMnsd = [];
        yGMnsd = [];
        xWMnsd = [];
        yWMnsd = [];
        nsGM =[];
        nsWM=[];

        A = blobMeasurements.VoxelIdxList{i};
        if sum(ismember(roiinds,A))>0    
            xGM = vertcat(xGM, t1si(allinds{i}).*GM(allinds{i}));
            yGM = vertcat(yGM, t2si(allinds{i}).*GM(allinds{i}));
            xGMn = vertcat(xGMn, mT1(allinds{i}).*GM(allinds{i}));
            yGMn = vertcat(yGMn, mT2(allinds{i}).*GM(allinds{i}));
            xGMnsd = vertcat(xGMnsd, sdT1(allinds{i}).*GM(allinds{i}));
            yGMnsd = vertcat(yGMnsd, sdT2(allinds{i}).*GM(allinds{i}));
            nsGM = vertcat(nsGM, n_matrix(allinds{i}).*GM(allinds{i}));
            xGMc = xGM;
            yGMc = yGM;
            xGM((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];
            xGMn((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMn((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            xGMnsd((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMnsd((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            nsGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];

            xWM = vertcat(xWM, t1si(allinds{i}).*WM(allinds{i}));
            yWM = vertcat(yWM, t2si(allinds{i}).*WM(allinds{i}));
            xWMn = vertcat(xWMn, mT1(allinds{i}).*WM(allinds{i}));
            yWMn = vertcat(yWMn, mT2(allinds{i}).*WM(allinds{i}));
            xWMnsd = vertcat(xWMnsd, sdT1(allinds{i}).*WM(allinds{i}));
            yWMnsd = vertcat(yWMnsd, sdT2(allinds{i}).*WM(allinds{i}));
            nsWM = vertcat(nsWM, n_matrix(allinds{i}).*WM(allinds{i}));
            xWMc = xWM;
            yWMc = yWM;
            xWM((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMn((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMn((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMnsd((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMnsd((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            nsWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];

            if (isempty(xGM) == 0 & isempty(yGM) == 0)|(isempty(xWM) == 0 & isempty(yWM) == 0)
                clusttp = vertcat(clusttp, 1);
                ptids = vertcat(ptids, p);
                sutys = vertcat(sutys, subtypes(cursor));

                cx = allcs(i,2);
                cy = allcs(i,1);
                cz = allcs(i,3);
                labnum = atlasi(cx,cy,cz);
                centers = vertcat(centers, [cx,cy,cz]);

                if labnum > 0
                    brainregion = dict(labnum,2);
                end
                parsed = split(string(brainregion.Var2),".");
                if length(parsed)>1
                    lobes = vertcat(lobes, parsed(2));
%                     gyri = vertcat(gyri, parsed(3));
                else
                    lobes = vertcat(lobes, parsed(1));
%                     gyri = vertcat(gyri, parsed(1));
                end

                xallG = vertcat(xallG, mean(xGM));
                yallG = vertcat(yallG, mean(yGM));
                xallGV = vertcat(xallGV, length(xGM));
                yallGV = vertcat(yallGV, length(yGM));
                xallW = vertcat(xallW, mean(xWM));
                yallW = vertcat(yallW, mean(yWM));
                xallWV = vertcat(xallWV, length(xWM));
                yallWV = vertcat(yallWV, length(yWM));

                ROIxGMsd = sqrt(sum(xGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIyGMsd = sqrt(sum(yGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIxWMsd = sqrt(sum(xWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
                ROIyWMsd = sqrt(sum(yWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
%                 xallGn = vertcat(xallGn, (mean(xGM) - mean(xGMn))/ROIxGMsd);
%                 yallGn = vertcat(yallGn, (mean(yGM) - mean(yGMn))/ROIyGMsd);
%                 xallWn = vertcat(xallWn, (mean(xWM) - mean(xWMn))/ROIxWMsd);
%                 yallWn = vertcat(yallWn, (mean(yWM) - mean(yWMn))/ROIyWMsd);

%                 xallGn = vertcat(xallGn, (mean(xGM) - sum(xGMn.*nsGM)/sum(nsGM))/ROIxGMsd);
%                 yallGn = vertcat(yallGn, (mean(yGM) - sum(yGMn.*nsGM)/sum(nsGM))/ROIyGMsd);
%                 xallWn = vertcat(xallWn, (mean(xWM) - sum(xWMn.*nsWM)/sum(nsWM))/ROIxWMsd);
%                 yallWn = vertcat(yallWn, (mean(yWM) - sum(yWMn.*nsWM)/sum(nsWM))/ROIyWMsd);
                
                xallGn = vertcat(xallGn, mean(xGM)/(sum(xGMn.*nsGM)/sum(nsGM)));
                yallGn = vertcat(yallGn, mean(yGM)/(sum(yGMn.*nsGM)/sum(nsGM)));
                xallWn = vertcat(xallWn, mean(xWM)/(sum(xWMn.*nsWM)/sum(nsWM)));
                yallWn = vertcat(yallWn, mean(yWM)/(sum(yWMn.*nsWM)/sum(nsWM)));
                        
            end
        else
            xGM = vertcat(xGM, t1si(allinds{i}).*GM(allinds{i}));
            yGM = vertcat(yGM, t2si(allinds{i}).*GM(allinds{i}));
            xGMn = vertcat(xGMn, mT1(allinds{i}).*GM(allinds{i}));
            yGMn = vertcat(yGMn, mT2(allinds{i}).*GM(allinds{i}));
            xGMnsd = vertcat(xGMnsd, sdT1(allinds{i}).*GM(allinds{i}));
            yGMnsd = vertcat(yGMnsd, sdT2(allinds{i}).*GM(allinds{i}));
            nsGM = vertcat(nsGM, n_matrix(allinds{i}).*GM(allinds{i}));
            xGMc = xGM;
            yGMc = yGM;
            xGM((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];
            xGMn((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMn((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            xGMnsd((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMnsd((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            nsGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];

            xWM = vertcat(xWM, t1si(allinds{i}).*WM(allinds{i}));
            yWM = vertcat(yWM, t2si(allinds{i}).*WM(allinds{i}));
            xWMn = vertcat(xWMn, mT1(allinds{i}).*WM(allinds{i}));
            yWMn = vertcat(yWMn, mT2(allinds{i}).*WM(allinds{i}));
            xWMnsd = vertcat(xWMnsd, sdT1(allinds{i}).*WM(allinds{i}));
            yWMnsd = vertcat(yWMnsd, sdT2(allinds{i}).*WM(allinds{i}));
            nsWM = vertcat(nsWM, n_matrix(allinds{i}).*WM(allinds{i}));
            xWMc = xWM;
            yWMc = yWM;
            xWM((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMn((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMn((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMnsd((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMnsd((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            nsWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];

            if (isempty(xGM) == 0 & isempty(yGM) == 0)|(isempty(xWM) == 0 & isempty(yWM) == 0)
                clusttp = vertcat(clusttp, 0);
                ptids = vertcat(ptids, p);
                sutys = vertcat(sutys, subtypes(cursor));
                
                cx = allcs(i,2);
                cy = allcs(i,1);
                cz = allcs(i,3);
                labnum = atlasi(cx,cy,cz);
                centers = vertcat(centers, [cx,cy,cz]);
                if labnum > 0
                    brainregion = dict(labnum,2);
                end
                parsed = split(string(brainregion.Var2),".");
                if length(parsed)>1
                    lobes = vertcat(lobes, parsed(2));
%                     gyri = vertcat(gyri, parsed(3));
                else
                    lobes = vertcat(lobes, parsed(1));
%                     gyri = vertcat(gyri, parsed(1));
                end
                xallG = vertcat(xallG, mean(xGM));
                yallG = vertcat(yallG, mean(yGM));
                xallGV = vertcat(xallGV, length(xGM));
                yallGV = vertcat(yallGV, length(yGM));
                xallW = vertcat(xallW, mean(xWM));
                yallW = vertcat(yallW, mean(yWM));
                xallWV = vertcat(xallWV, length(xWM));
                yallWV = vertcat(yallWV, length(yWM));

                ROIxGMsd = sqrt(sum(xGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIyGMsd = sqrt(sum(yGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIxWMsd = sqrt(sum(xWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
                ROIyWMsd = sqrt(sum(yWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
%                 xallGn = vertcat(xallGn, (mean(xGM) - sum(xGMn.*nsGM)/sum(nsGM))/ROIxGMsd);
%                 yallGn = vertcat(yallGn, (mean(yGM) - sum(yGMn.*nsGM)/sum(nsGM))/ROIyGMsd);
%                 xallWn = vertcat(xallWn, (mean(xWM) - sum(xWMn.*nsWM)/sum(nsWM))/ROIxWMsd);
%                 yallWn = vertcat(yallWn, (mean(yWM) - sum(yWMn.*nsWM)/sum(nsWM))/ROIyWMsd);

                xallGn = vertcat(xallGn, mean(xGM)/(sum(xGMn.*nsGM)/sum(nsGM)));
                yallGn = vertcat(yallGn, mean(yGM)/(sum(yGMn.*nsGM)/sum(nsGM)));
                xallWn = vertcat(xallWn, mean(xWM)/(sum(xWMn.*nsWM)/sum(nsWM)));
                yallWn = vertcat(yallWn, mean(yWM)/(sum(yWMn.*nsWM)/sum(nsWM)));
            end
        end
    end


%     figure()
%     scatter(xtWM, ytWM, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     hold on
%     scatter(xfWM, yfWM, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
%     title(p,"WM")
%     xlabel('T1(ms)');
%     ylabel('T2(ms)');
end

% figure()
% scatter(xtallW, ytallW, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% hold on
% scatter(xfallW, yfallW, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% title('all patients WM')
% xlabel('T1(ms)');
% ylabel('T2(ms)');
% xlim([0 2000]);
% ylim([0 400]);

% figure()
% scatter(xtallG, ytallG, 6, 'r', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% hold on
% scatter(xfallG, yfallG, 6, 'b', 'filled','MarkerFaceAlpha',.35,'MarkerEdgeAlpha',.35)
% title('all patients GM')
% xlabel('T1(ms)');
% ylabel('T2(ms)');
% xlim([0 2000]);
% ylim([0 400]);

xallG(isnan(xallG)) = 0;
xallW(isnan(xallW)) = 0;
yallG(isnan(yallG)) = 0;
yallW(isnan(yallW)) = 0;
ad = horzcat(xallG,xallW,yallG,yallW,yallGV, yallWV);  % combine 
% ad(isnan(ad)) = 0;
% adw = horzcat(xallG.*xallGV./sum(xallGV),xallW.*xallWV./sum(xallWV),yallG.*yallGV./sum(yallGV),yallW.*yallWV./sum(yallWV));
% adnov = horzcat(xallG,xallW,yallG,yallW);
% adnovG = horzcat(xallG,yallG);
% adnovW = horzcat(xallW,yallW);

xallGn(isnan(xallGn)) = 0;
xallWn(isnan(xallWn)) = 0;
yallGn(isnan(yallGn)) = 0;
yallWn(isnan(yallWn)) = 0;
adn = horzcat(xallGn,xallWn, yallGn,yallWn,yallGV, yallWV);  % combine 

adn(isinf(adn)|isnan(adn)) = 0;
ad(isinf(ad)|isnan(ad)) = 0;
% normed =horzcat(adn, ptids, sutys, clusttp);
% writematrix(normed, 'norm_trainingatpoint0.xlsx')
% unnormed =horzcat(ad, ptids, sutys, clusttp);
% writematrix(unnormed, 'unnorm_trainingatpoint0.xlsx')
% lesinfo =horzcat(ptids, sutys, lobes, centers,yallGV,yallWV,clusttp);
% writematrix(lesinfo, 'lesinfo_trainingatpoint0.xlsx')

normed =horzcat(adn, ptids, sutys, clusttp);
writematrix(normed, 'norm_trainingatpoint0old.xlsx')
unnormed =horzcat(ad, ptids, sutys, clusttp);
writematrix(unnormed, 'unnorm_trainingatpoint0old.xlsx')
lesinfo =horzcat(ptids, sutys, lobes, centers,yallGV,yallWV,clusttp);
writematrix(lesinfo, 'lesinfo_trainingatpoint0old.xlsx')

% normed =horzcat(adn, ptids, sutys, clusttp);
% writematrix(normed, 'norm_trainingatpoint0_fn.xlsx')
% unnormed =horzcat(ad, ptids, sutys, clusttp);
% writematrix(unnormed, 'unnorm_trainingatpoint0_fn.xlsx')
% lesinfo =horzcat(ptids, sutys, lobes, centers,yallGV,yallWV,clusttp);
% writematrix(lesinfo, 'lesinfo_trainingatpoint0_fn.xlsx')

end
%% tsne all
cd C:\Users\eegrvw\Downloads
clusttp = readmatrix("response_15.csv");
adn = readmatrix("preds_norm_15.csv");

% Y1 = tsne(ad(:,1:4), 'Perplexity', 30);
% Y2 = tsne(ad(:,1:6), 'Perplexity', 30);
Y3 = tsne(adn(:,1:4), 'Perplexity', 40);
Y4 = tsne(adn(:,1:6), 'Perplexity', 40);
tcl = tiledlayout(2,2);
lcn = 'southwest';
labels = string(clusttp);
labels(clusttp == 1) = "TP clusters";
labels(clusttp == 0) = "FP clusters";
nexttile(tcl)
% gscatter(Y1(:,1),Y1(:,2),labels,'rb');
title('MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
% gscatter(Y2(:,1),Y2(:,2),labels,'rb');
title('MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y3(:,1),Y3(:,2),labels,'rb');
title('Normalized MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y4(:,1),Y4(:,2),labels,'rb');
title('Normalized MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
%%
GMtrain = horzcat(adnovG((clusttp == 1),1), adnovG((clusttp == 1),2));
GMtrain = smote(GMtrain, 2, 5);
GMtrain = smote(GMtrain, 2, 5);
WMtrain = horzcat(adnovW((clusttp == 1),1), adnovW((clusttp == 1),2));
WMtrain = smote(WMtrain, 2, 5);
WMtrain = smote(WMtrain, 2, 5);
GMcontrol = horzcat(adnovG((clusttp == 0),1), adnovG((clusttp == 0),2));
WMcontrol = horzcat(adnovW((clusttp == 0),1), adnovW((clusttp == 0),2));
adnovG = vertcat(GMcontrol,GMtrain);
adnovW = vertcat(WMcontrol,WMtrain);
Glbl = vertcat(zeros(size(GMcontrol,1),1),ones(size(GMtrain,1),1));
Wlbl = vertcat(zeros(size(WMcontrol,1),1),ones(size(WMtrain,1),1));

%% do SVM and plot
X = adnovG;
lbl = categorical(clusttp,[0 1],{'normal' 'lesion'});
% mdl = fitcsvm(X, lbl, 'OptimizeHyperparameters','auto', ...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
%     'expected-improvement-plus'));
mdl = fitcsvm(X, lbl);
sv = mdl.SupportVectors;
CV = crossval(mdl);
kfoldLoss(CV)
figure
gscatter(X(:,1),X(:,2),lbl)
hold on
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)


set(gcf,'Visible','on')
scatter(X(:,1),X(:,2),10,wlbl,'filled','MarkerEdgeColor','r' )
xlabel('T1');ylabel('T2');
% hold on
% plot3(sv(:,1),sv(:,2),sv(:,3),'ko','MarkerSize',10)
% hold on
% numGrid = 100;
% [x1Grid,x2Grid,x3Grid] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),numGrid),...
%     linspace(min(X(:,2)),max(X(:,2)),numGrid),linspace(min(X(:,3)),max(X(:,3)),numGrid));
% xList = [x1Grid(:),x2Grid(:),x3Grid(:)];
% [~,scores] = predict(mdlw,xList);
% [faces,verts] = isosurface(x1Grid,x2Grid,x3Grid, reshape(scores(:,2),size(x1Grid)),0);
% p=patch('Vertices', verts, 'Faces', faces, 'FaceColor','k','edgecolor', 'none', 'FaceAlpha', 0.5);
% p.FaceColor = 'red';
% grid on; box on

CV = crossval(mdl);
kfoldLoss(CV)

% mdlSVM = fitPosterior(mdlw);
% [~,score_svm] = resubPredict(mdlSVM);
% [Xsvm,Ysvm,~,AUCsvm] = perfcurve(clusttpw,score_svm(:,2),1);
% figure()
% plot(Xsvm,Ysvm)

%%
Xg = horzcat(gt1,gt2,xallGV);
% X = horzcat(wt1,wt2);
glbl = categorical(clusttpg,[0 1],{'normal' 'lesion'});
% mdlw = fitcsvm(X, wlbl, 'OptimizeHyperparameters','auto', ...
%     'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
%     'expected-improvement-plus'));
mdlg = fitcsvm(Xg, glbl);
sv = mdlg.SupportVectors;
figure
set(gcf,'Visible','on')
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),10,glbl,'filled','MarkerEdgeColor','r' );
xlabel('T1');ylabel('T2');zlabel('volume');

CVg = crossval(mdlg);
kfoldLoss(CVg)
%%
x = adnovG;
r = Glbl'; %first row lesion
r2 = -(r-1); %second row normal
r = vertcat(r,r2);
% net = patternnet(4);
% net.layers{1}.transferFcn = 'tansig';
% net.layers{2}.transferFcn = 'softmax';
%% organize data for non-normalized
poscases = ad((clusttp==1),:);
posres = ones(size(poscases,1),1);
negcases = ad((clusttp==0),:);
negres = zeros(size(negcases,1),1);
poscases = repmat(poscases,14,1);
posres = repmat(posres,14,1);
allcases = vertcat(poscases, negcases);
allres = vertcat(posres, negres);
indices = randperm(size(allcases,1));
shcases = allcases(indices,:);
shres = allres(indices,:);
oa = 0;
x = shcases;
r = shres;
r = r';
r2 = -(r-1);
r = vertcat(r,r2);


fold = 5;
for i = 1:fold-1
    if i == 1
        x = shcases(round(2/fold*size(shcases,1)):end,:);
        r = shcases(round(2/fold*size(shcases,1)):end,:);
    else
        x = vertcat(shcases(1:round((i-1)/fold*size(negcases,1)),:),...
            shcases(round((i+1)/fold*size(negcases,1)):end,:));
        r = vertcat(shres(1:round((i-1)/fold*size(negcases,1)),:),...
            shres(round((i+1)/fold*size(negcases,1)):end,:));
    end
    r = r';
    r2 = -(r-1);
    r = vertcat(r,r2);
    
    xv = shcases(round(i/fold*size(shcases,1)):round((i+1)/fold*size(shcases,1)),:);
    rv = shres(round(i/fold*size(shcases,1)):round((i+1)/fold*size(shcases,1)),:);
    rv = rv';
    rv2 = -(rv-1);
    rv = vertcat(rv,rv2);
    net = patternnet([6, 6, 6]);
    net.layers{1}.transferFcn = 'tansig';
    net.layers{2}.transferFcn = 'tansig';
    net.layers{3}.transferFcn = 'tansig';
    net.layers{4}.transferFcn = 'softmax';

%     net = patternnet([8, 8, 8]);
%     net.layers{1}.transferFcn = 'poslin';
%     net.layers{2}.transferFcn = 'poslin';
%     net.layers{3}.transferFcn = 'poslin';
%     net.layers{4}.transferFcn = 'softmax';
    
    net.trainFcn = 'trainscg';
    net.trainParam.epochs = 500;
    net.trainParam.max_fail = 20;
    net.divideParam.trainRatio = 100/100;
    net.divideParam.valRatio = 0/100;
    net.divideParam.testRatio = 0/100;
    
    net = train(net,x',r);
    net.trainFcn = 'trainscg';
    y = net(xv');
    y(y>=0.5) = 1;
    y(y<0.5) = 0;
%     size(y)
    accuracy = sum(y(1,:)==rv(1,:))/size(y(1,:),2)
%     perf = perform(net,y,r)
    oa = oa + accuracy;
end
oa = oa/(fold-1)
%% organize data for normalized
poscases = adn((clusttp==1),:);
posres = ones(size(poscases,1),1);
negcases = adn((clusttp==0),:);
negres = zeros(size(negcases,1),1);
poscases = repmat(poscases,14,1);
posres = repmat(posres,14,1);
allcases = vertcat(poscases, negcases);
allres = vertcat(posres, negres);
indices = randperm(size(allcases,1));
shcases = allcases(indices,:);
shres = allres(indices,:);
oa = 0;
x = shcases;
r = shres;
r = r';
r2 = -(r-1);
r = vertcat(r,r2);
%%
x = ad;
r = clusttp'; %first row lesion
r2 = -(r-1); %second row normal
r = vertcat(r,r2);
net = patternnet([6 6]);
net.layers{1}.transferFcn = 'tansig';
net.layers{2}.transferFcn = 'tansig';
net.layers{3}.transferFcn = 'softmax';

net.trainFcn = 'trainscg';
net.trainParam.epochs = 500;
net.trainParam.max_fail = 10;
net.divideParam.trainRatio = 80/100;
net.divideParam.valRatio = 00/100;
net.divideParam.testRatio =20/100;

net = train(net,x',r);
net.trainFcn = 'trainscg';
y = net(x');
perf = perform(net,y,r)

%% tsne none zero
ad0s = ad(:,5).*ad(:,6);
adn0s = adn(:,5).*adn(:,6);
adn0 = ad(ad0s~=0,:);
adnn0 = adn(adn0s~=0,:);
Y1 = tsne(adn0(:,1:4), 'Perplexity', 20);
Y2 = tsne(adn0, 'Perplexity', 20);
Y3 = tsne(adnn0(:,1:4), 'Perplexity', 20);
Y4 = tsne(adnn0, 'Perplexity', 20);
clusttpn0 = clusttp(ad0s~=0,:);

tcl = tiledlayout(2,2);
lcn = 'southwest';
labels = string(clusttpn0);
labels(clusttpn0 == 1) = "TP clusters";
labels(clusttpn0 == 0) = "FP clusters";
nexttile(tcl)
gscatter(Y1(:,1),Y1(:,2),labels);
title('MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y2(:,1),Y2(:,2),labels);
title('MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y3(:,1),Y3(:,2),labels);
title('Normalized MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y4(:,1),Y4(:,2),labels);
title('Normalized MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;

%% tsne WM only

adn0 = ad(ad(:,6)~=0,[2,4,6]);
adnn0 = adn(adn(:,6)~=0,[2,4,6]);
Y1 = tsne(adn0(:,1:3), 'Perplexity', 50);
Y2 = tsne(adn0, 'Perplexity', 50);
Y3 = tsne(adnn0(:,1:3), 'Perplexity', 50);
Y4 = tsne(adnn0, 'Perplexity', 50);
clusttpn0 = clusttp(ad(:,6)~=0,:);

tcl = tiledlayout(2,2);
lcn = 'southwest';
labels = string(clusttpn0);
labels(clusttpn0 == 1) = "TP clusters";
labels(clusttpn0 == 0) = "FP clusters";
nexttile(tcl)
gscatter(Y1(:,1),Y1(:,2),labels);
title('MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y2(:,1),Y2(:,2),labels);
title('MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y3(:,1),Y3(:,2),labels);
title('Normalized MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y4(:,1),Y4(:,2),labels);
title('Normalized MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
%% tsne divide volume
ad0s = ad(:,5).*ad(:,6);
adn0s = adn(:,5).*adn(:,6);
adn0 = ad(ad0s~=0,:);
adnn0 = adn(adn0s~=0,:);
% adn0 = adn0(:,1:4);
% adnn0 = adnn0(:,1:4);
adn0(:,1) = adn0(:,1)./adn0(:,5);
adn0(:,2) = adn0(:,2)./adn0(:,6);
adn0(:,3) = adn0(:,3)./adn0(:,5);
adn0(:,4) = adn0(:,4)./adn0(:,6);
adnn0(:,1) = adnn0(:,1)./adn0(:,5);
adnn0(:,2) = adnn0(:,2)./adn0(:,6);
adnn0(:,3) = adnn0(:,3)./adn0(:,5);
adnn0(:,4) = adnn0(:,4)./adn0(:,6);

Y1 = tsne(adn0(:,1:4), 'Perplexity', 20);
Y2 = tsne(adn0, 'Perplexity', 20);
Y3 = tsne(adnn0(:,1:4), 'Perplexity', 20);
Y4 = tsne(adnn0, 'Perplexity', 20);
clusttpn0 = clusttp(ad0s~=0,:);

tcl = tiledlayout(2,2);
lcn = 'southwest';
labels = string(clusttpn0);
labels(clusttpn0 == 1) = "TP clusters";
labels(clusttpn0 == 0) = "FP clusters";
nexttile(tcl)
gscatter(Y1(:,1),Y1(:,2),labels);
title('MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y2(:,1),Y2(:,2),labels);
title('MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y3(:,1),Y3(:,2),labels);
title('Normalized MRF data', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
nexttile(tcl)
gscatter(Y4(:,1),Y4(:,2),labels);
title('Normalized MRF data and Volume', 'FontSize', 20)
legend('Location',lcn)
ax=gca;
ax.FontSize = 15;
%% z score based
cd T:\Imaging\Multimodal\MRF\Peter
% load meanMRFvalues.mat
% load MRFstdev_voxel.mat

subjID=["PXX_XXXXX"];
subtypes = ["IIB" "mMCD" "mMCD"...
    "IIB" "IIB" "mMCD" "mMCD" "IIA"...
    "IIA" "IIA" "IIB" "IIA" "IIB"...
    "mMCD" "mMCD" "MOGHE" "IIB" "MOGHE"...
    "IIA" "IIB"...
    "IIA" "IIA" "mMCD" "IIA" "IIB"...
    "IIB" "mMCD" "IIA" "IIB" "mMCD"...
    "IIA" "IIB" "mMCD"];

MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
dict = readtable('atlas.xlsx');
atlas = load_untouch_nii('Talairach-labels-1mm.nii');
atlasi = single(atlas.img);

for volthresh = 1
xallW = [];
yallW = [];
xallG = [];
yallG = [];
xallWn = [];
yallWn = [];
xallGn = [];
yallGn = [];
xallWV = [];
yallWV = [];
xallGV = [];
yallGV = [];
clusttp =[];
ptids = [];
sutys = [];
lobes = [];
% gyri = [];
centers = [];

hastps = [];
t1max = 1759; 
t2max = 106;
t1min = 743;
t2min = 29;
% t1max = 1700; 
% t2max = 100;
% t1min = 750;
% t2min = 30;

cursor = 0;
for p = subjID
    p
    cursor = cursor + 1;
    path = strcat(MRF_path,'\',p,'\MRF_VBM');
    cd(path)
    ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
    anni = single(ann.img);
    anni(anni>0) = 1;

    if p == "P31_XXXXX"
        anni(:,242:end,1:57) = 0;
    end
    mask = single(load_untouch_nii('MNI_GM.nii').img) + single(load_untouch_nii('MNI_WM.nii').img);
    mask(mask>0.95) = 1;
%     mask(mask<0.95) = 0;
    anni = anni.*mask;
    GM = single(load_untouch_nii('MNI_GM.nii').img).*mask;
    GM(GM<0.5) = 0;
    GM(GM>=0.5) = 1;
    WM = single(load_untouch_nii('MNI_WM.nii').img).*mask;
    WM(WM<0.5) = 0;
    WM(WM>=0.5) = 1;
    
    % GM compute and show graphs
    [labeledImage, numberOfBlobs] = bwlabeln(anni);
	blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList', 'Centroid');
	allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>volthresh)];
    allcs = round(blobMeasurements.Centroid);
    allcs = allcs(blobMeasurements.Volume>volthresh,:);
    numberOfBlobs = size(allinds,1);


    roi = load_untouch_nii('MNI_ROI_final.nii');
    roii = single(roi.img);
    if p == "P83_XXXXX"
        roii(:,269:end,:) = 0;
    end
    [roiblob, nBlobs] = bwlabeln(roii);
    roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
    roiinds = roiMeasurements.VoxelIdxList{1};
    t1s = load_untouch_nii('MNI_T1.nii');
    t1si = single(t1s.img);
    t2s = load_untouch_nii('MNI_T2.nii');
    t2si = single(t2s.img);

    for i = 1:numberOfBlobs
        xGM = [];
        yGM = [];
        xWM = [];
        yWM = [];
        xGMn = [];
        yGMn = [];
        xWMn = [];
        yWMn = [];
        xGMnsd = [];
        yGMnsd = [];
        xWMnsd = [];
        yWMnsd = [];
        nsGM =[];
        nsWM=[];

        A = blobMeasurements.VoxelIdxList{i};
        if sum(ismember(roiinds,A))>0    
            xGM = vertcat(xGM, t1si(allinds{i}).*GM(allinds{i}));
            yGM = vertcat(yGM, t2si(allinds{i}).*GM(allinds{i}));
            xGMn = vertcat(xGMn, mT1(allinds{i}).*GM(allinds{i}));
            yGMn = vertcat(yGMn, mT2(allinds{i}).*GM(allinds{i}));
            xGMnsd = vertcat(xGMnsd, sdT1(allinds{i}).*GM(allinds{i}));
            yGMnsd = vertcat(yGMnsd, sdT2(allinds{i}).*GM(allinds{i}));
            nsGM = vertcat(nsGM, n_matrix(allinds{i}).*GM(allinds{i}));
            xGMc = xGM;
            yGMc = yGM;
            xGM((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];
            xGMn((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMn((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            xGMnsd((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMnsd((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            nsGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];

            xWM = vertcat(xWM, t1si(allinds{i}).*WM(allinds{i}));
            yWM = vertcat(yWM, t2si(allinds{i}).*WM(allinds{i}));
            xWMn = vertcat(xWMn, mT1(allinds{i}).*WM(allinds{i}));
            yWMn = vertcat(yWMn, mT2(allinds{i}).*WM(allinds{i}));
            xWMnsd = vertcat(xWMnsd, sdT1(allinds{i}).*WM(allinds{i}));
            yWMnsd = vertcat(yWMnsd, sdT2(allinds{i}).*WM(allinds{i}));
            nsWM = vertcat(nsWM, n_matrix(allinds{i}).*WM(allinds{i}));
            xWMc = xWM;
            yWMc = yWM;
            xWM((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMn((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMn((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMnsd((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMnsd((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            nsWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];

            if (isempty(xGM) == 0 & isempty(yGM) == 0)|(isempty(xWM) == 0 & isempty(yWM) == 0)
                clusttp = vertcat(clusttp, 1);
                ptids = vertcat(ptids, p);
                sutys = vertcat(sutys, subtypes(cursor));

                cx = allcs(i,2);
                cy = allcs(i,1);
                cz = allcs(i,3);
                labnum = atlasi(cx,cy,cz);
                centers = vertcat(centers, [cx,cy,cz]);

                if labnum > 0
                    brainregion = dict(labnum,2);
                end
                parsed = split(string(brainregion.Var2),".");
                if length(parsed)>1
                    lobes = vertcat(lobes, parsed(2));
                else
                    lobes = vertcat(lobes, parsed(1));
                end

                xallG = vertcat(xallG, mean(xGM));
                yallG = vertcat(yallG, mean(yGM));
                xallGV = vertcat(xallGV, length(xGM));
                yallGV = vertcat(yallGV, length(yGM));
                xallW = vertcat(xallW, mean(xWM));
                yallW = vertcat(yallW, mean(yWM));
                xallWV = vertcat(xallWV, length(xWM));
                yallWV = vertcat(yallWV, length(yWM));

                ROIxGMsd = sqrt(sum(xGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIyGMsd = sqrt(sum(yGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIxWMsd = sqrt(sum(xWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
                ROIyWMsd = sqrt(sum(yWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
                xallGn = vertcat(xallGn, (mean(xGM) - mean(xGMn))/ROIxGMsd);
                yallGn = vertcat(yallGn, (mean(yGM) - mean(yGMn))/ROIyGMsd);
                xallWn = vertcat(xallWn, (mean(xWM) - mean(xWMn))/ROIxWMsd);
                yallWn = vertcat(yallWn, (mean(yWM) - mean(yWMn))/ROIyWMsd);

%                 xallGn = vertcat(xallGn, (mean(xGM) - sum(xGMn.*nsGM)/sum(nsGM))/ROIxGMsd);
%                 yallGn = vertcat(yallGn, (mean(yGM) - sum(yGMn.*nsGM)/sum(nsGM))/ROIyGMsd);
%                 xallWn = vertcat(xallWn, (mean(xWM) - sum(xWMn.*nsWM)/sum(nsWM))/ROIxWMsd);
%                 yallWn = vertcat(yallWn, (mean(yWM) - sum(yWMn.*nsWM)/sum(nsWM))/ROIyWMsd);
                
%                 xallGn = vertcat(xallGn, mean(xGM)/(sum(xGMn.*nsGM)/sum(nsGM)));
%                 yallGn = vertcat(yallGn, mean(yGM)/(sum(yGMn.*nsGM)/sum(nsGM)));
%                 xallWn = vertcat(xallWn, mean(xWM)/(sum(xWMn.*nsWM)/sum(nsWM)));
%                 yallWn = vertcat(yallWn, mean(yWM)/(sum(yWMn.*nsWM)/sum(nsWM)));
                        
            end
        else
            xGM = vertcat(xGM, t1si(allinds{i}).*GM(allinds{i}));
            yGM = vertcat(yGM, t2si(allinds{i}).*GM(allinds{i}));
            xGMn = vertcat(xGMn, mT1(allinds{i}).*GM(allinds{i}));
            yGMn = vertcat(yGMn, mT2(allinds{i}).*GM(allinds{i}));
            xGMnsd = vertcat(xGMnsd, sdT1(allinds{i}).*GM(allinds{i}));
            yGMnsd = vertcat(yGMnsd, sdT2(allinds{i}).*GM(allinds{i}));
            nsGM = vertcat(nsGM, n_matrix(allinds{i}).*GM(allinds{i}));
            xGMc = xGM;
            yGMc = yGM;
            xGM((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];
            xGMn((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMn((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            xGMnsd((xGMc > t1max)|(xGMc< t1min)|(yGMc > t2max)|(yGMc< t2min)) = []; 
            yGMnsd((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = []; 
            nsGM((yGMc > t2max)|(yGMc< t2min)|(xGMc > t1max)|(xGMc< t1min)) = [];

            xWM = vertcat(xWM, t1si(allinds{i}).*WM(allinds{i}));
            yWM = vertcat(yWM, t2si(allinds{i}).*WM(allinds{i}));
            xWMn = vertcat(xWMn, mT1(allinds{i}).*WM(allinds{i}));
            yWMn = vertcat(yWMn, mT2(allinds{i}).*WM(allinds{i}));
            xWMnsd = vertcat(xWMnsd, sdT1(allinds{i}).*WM(allinds{i}));
            yWMnsd = vertcat(yWMnsd, sdT2(allinds{i}).*WM(allinds{i}));
            nsWM = vertcat(nsWM, n_matrix(allinds{i}).*WM(allinds{i}));
            xWMc = xWM;
            yWMc = yWM;
            xWM((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMn((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMn((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            xWMnsd((xWMc > t1max)|(xWMc < t1min)|(yWMc > t2max)|(yWMc< t2min)) = []; 
            yWMnsd((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];
            nsWM((yWMc > t2max)|(yWMc < t2min)|(xWMc > t1max)|(xWMc< t1min)) = [];

            if (isempty(xGM) == 0 & isempty(yGM) == 0)|(isempty(xWM) == 0 & isempty(yWM) == 0)
                clusttp = vertcat(clusttp, 0);
                ptids = vertcat(ptids, p);
                sutys = vertcat(sutys, subtypes(cursor));
                
                cx = allcs(i,2);
                cy = allcs(i,1);
                cz = allcs(i,3);
                labnum = atlasi(cx,cy,cz);
                centers = vertcat(centers, [cx,cy,cz]);
                if labnum > 0
                    brainregion = dict(labnum,2);
                end
                parsed = split(string(brainregion.Var2),".");
                if length(parsed)>1
                    lobes = vertcat(lobes, parsed(2));
%                     gyri = vertcat(gyri, parsed(3));
                else
                    lobes = vertcat(lobes, parsed(1));
%                     gyri = vertcat(gyri, parsed(1));
                end
                xallG = vertcat(xallG, mean(xGM));
                yallG = vertcat(yallG, mean(yGM));
                xallGV = vertcat(xallGV, length(xGM));
                yallGV = vertcat(yallGV, length(yGM));
                xallW = vertcat(xallW, mean(xWM));
                yallW = vertcat(yallW, mean(yWM));
                xallWV = vertcat(xallWV, length(xWM));
                yallWV = vertcat(yallWV, length(yWM));

                ROIxGMsd = sqrt(sum(xGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIyGMsd = sqrt(sum(yGMnsd.^2.*(nsGM-1))/sum((nsGM-1)));
                ROIxWMsd = sqrt(sum(xWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
                ROIyWMsd = sqrt(sum(yWMnsd.^2.*(nsWM-1))/sum((nsWM-1)));
                xallGn = vertcat(xallGn, (mean(xGM) - mean(xGMn))/ROIxGMsd);
                yallGn = vertcat(yallGn, (mean(yGM) - mean(yGMn))/ROIyGMsd);
                xallWn = vertcat(xallWn, (mean(xWM) - mean(xWMn))/ROIxWMsd);
                yallWn = vertcat(yallWn, (mean(yWM) - mean(yWMn))/ROIyWMsd);
%                 xallGn = vertcat(xallGn, (mean(xGM) - sum(xGMn.*nsGM)/sum(nsGM))/ROIxGMsd);
%                 yallGn = vertcat(yallGn, (mean(yGM) - sum(yGMn.*nsGM)/sum(nsGM))/ROIyGMsd);
%                 xallWn = vertcat(xallWn, (mean(xWM) - sum(xWMn.*nsWM)/sum(nsWM))/ROIxWMsd);
%                 yallWn = vertcat(yallWn, (mean(yWM) - sum(yWMn.*nsWM)/sum(nsWM))/ROIyWMsd);

            end
        end
    end


end

xallG(isnan(xallG)) = 0;
xallW(isnan(xallW)) = 0;
yallG(isnan(yallG)) = 0;
yallW(isnan(yallW)) = 0;
ad = horzcat(xallG,xallW,yallG,yallW,yallGV, yallWV);  % combine 

xallGn(isnan(xallGn)) = 0;
xallWn(isnan(xallWn)) = 0;
yallGn(isnan(yallGn)) = 0;
yallWn(isnan(yallWn)) = 0;
adn = horzcat(xallGn,xallWn, yallGn,yallWn,yallGV, yallWV);  % combine 

adn(isinf(adn)|isnan(adn)) = 0;
ad(isinf(ad)|isnan(ad)) = 0;

normed =horzcat(adn, ptids, sutys, clusttp);
% writematrix(normed, 'norm_trainingatpoint0old_Bronly.xlsx')
% unnormed =horzcat(ad, ptids, sutys, clusttp);
% writematrix(unnormed, 'unnorm_trainingatpoint0old_Bronly.xlsx')
% lesinfo =horzcat(ptids, sutys, lobes, centers,yallGV,yallWV,clusttp);
% writematrix(lesinfo, 'lesinfo_trainingatpoint0old_Bronly.xlsx')
cd Z:\Imaging\Multimodal\MRF\Peter
writematrix(normed, 'norm_trainingatpoint0old_excRadio.xlsx')
unnormed =horzcat(ad, ptids, sutys, clusttp);
writematrix(unnormed, 'unnorm_trainingatpoint0old_excRadio.xlsx')
lesinfo =horzcat(ptids, sutys, lobes, centers,yallGV,yallWV,clusttp);
writematrix(lesinfo, 'lesinfo_trainingatpoint0old_excRadio.xlsx')

end