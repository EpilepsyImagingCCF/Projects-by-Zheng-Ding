%% all promising patientscluster level
subjID=["PXX_XXXXX"];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
tplst = [];
fplst = [];
% for thres = 0.1:0.1:0.9
for thres = 0
    thres
    tpth = 0;
    fpth = 0;
    tnth = 0;
    fnth = 0;
    for p = subjID
        path = strcat(MRF_path,'\',p,'\MRF_VBM');
        cd(path)
        ann = load_untouch_nii('MNI_lesionprob_ANN_fn.nii');
        anni = single(ann.img);
        anni(anni>=thres) = 1;
        anni(anni<thres) = 0;
        mask = single(load_untouch_nii('MNI_GM_fn.nii').img) + single(load_untouch_nii('MNI_WM_fn.nii').img);
        mask(mask>=0.99) = 1;
        mask(mask<0.99) = 0;
        anni = anni.*mask;
        

        [labeledImage, numberOfBlobs] = bwlabeln(anni);
	    blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	    allinds = [blobMeasurements.VoxelIdxList];
        numberOfBlobs = size(allinds,1);
    
        roi = load_untouch_nii('MNI_ROI_final.nii');
        roii = single(roi.img);
        [roiblob, nBlobs] = bwlabeln(roii);
        roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
        roiinds = roiMeasurements.VoxelIdxList{1};
        
        trueind = [];
        for i = 1:numberOfBlobs
            A = blobMeasurements.VoxelIdxList{i};
            if sum(ismember(roiinds,A))>0
                trueind = horzcat(trueind, i);
            end
        end
    
        for i = 1:numberOfBlobs
            if ismember(i,trueind) == 1
                tpth = tpth + 1;
            else
                fpth = fpth + 1;
            end
        end
    
        ann = load_untouch_nii('MNI_lesionprob_ANN_fn.nii');
        anni = single(ann.img);
        anni(anni>=thres) = 0;
        anni((anni<thres)&(anni>0)) = 1;
        mask = single(load_untouch_nii('MNI_Brain_Mask.nii').img);
        mask(mask>0.95) = 1;
        anni = anni.*mask;
    
        [labeledImage, numberOfBlobs] = bwlabeln(anni);
	    blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	    allinds = [blobMeasurements.VoxelIdxList];
        numberOfBlobs = size(allinds,1);
    
        roi = load_untouch_nii('MNI_ROI_final.nii');
        roii = single(roi.img);
        [roiblob, nBlobs] = bwlabeln(roii);
        roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
        roiinds = roiMeasurements.VoxelIdxList{1};
        
        trueind = [];
        for i = 1:numberOfBlobs
            A = blobMeasurements.VoxelIdxList{i};
            if sum(ismember(roiinds,A))>0
                trueind = horzcat(trueind, i);
            end
        end
    
        for i = 1:numberOfBlobs
            if ismember(i,trueind) == 1
                fnth = fnth + 1;
            else
                tnth = tnth + 1;
            end
        end
    end
    [tpth, fpth, fnth, tnth]
    
    acc = (tpth+tnth)/(tpth + fpth +fnth + tnth)*100;
    tpr = (tpth)/(tpth + fnth)*100;
    fpr = (fpth)/(tnth + fpth)*100;
    fnr = (fnth)/(tpth + fnth)*100;
    tnr = (tnth)/(tnth + fpth)*100;
    [acc, tpr, fpr, fnr, tnr]

    tplst = vertcat(tplst, tpr);
    fplst = vertcat(fplst, fpr);
end
figure()
plot(fplst, tplst)

%% Patient level
subjID=["PXX_XXXXX"];
subtypes = ["IIB" "MOGHE" "IIB" "mMCD"...
    "IIB" "IIB" "mMCD" "mMCD" "IIA"...
    "IIA" "IIA" "IIB" "IIA" "IIB"...
    "mMCD" "mMCD" "MOGHE" "IIB" "MOGHE"...
    "IIA" "IIB" "IIB" "IIB" "IIB"...
    "IIA" "IIA" "mMCD" "IIA" "IIB"...
    "IIB" "mMCD" "IIA" "IIB" "mMCD"...
    "IIA" "IIB" "mMCD"];


MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

num_pt_clust_lst = [];
tplst = [];
fplst_0 = [];
fplst_1 = [];
fplst_2 = [];
fplst_3 = [];
fplst_4 = [];
fnlst = [];
tnlst = [];
fpthlst = [];

% for thres = 0.1:0.1:0.9
for thres = 0.001
    num_pt_clust = 0;
    num_pt_tp = 0;
    num_pt_fp0 = 0;
    num_pt_fp1 = 0;
    num_pt_fp2 = 0;
    num_pt_fp3 = 0;
    num_pt_fp4 = 0;
    num_pt_fn = 0;
    num_pt_tn = 0;
    sub_cursor = 1;
    keys = {'IIA'; 'IIB'; 'mMCD'; 'MOGHE'};
    FP0 = [0; 0; 0; 0];
    FP1 = [0; 0; 0; 0];
    FP2 = [0; 0; 0; 0];
    FP3 = [0; 0; 0; 0];
    FP4 = [0; 0; 0; 0];
    T = table(FP0, FP1, FP2, FP3, FP4, 'RowNames', keys);
    totalfp = 0;

    for p = subjID
        fpth = 0;
        tnth = 0;
        fnth = 0;


        path = strcat(MRF_path,'\',p,'\MRF_VBM');
        cd(path)
        ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
        anni = single(ann.img);
        anni(anni>=thres) = 1;
        anni(anni<thres) = 0;
        mask = single(load_untouch_nii('MNI_GM.nii').img) + single(load_untouch_nii('MNI_WM.nii').img);
        mask(mask>=0.95) = 1;
        mask(mask<0.95) = 0;
        anni = anni.*mask;
        detected = 0;

        % above thresh
        [labeledImage, numberOfBlobs] = bwlabeln(anni);

	    blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	    allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>1)];
        numberOfBlobs = size(allinds,1);
    
        roi = load_untouch_nii('MNI_ROI_final.nii');
        roii = single(roi.img);
        [roiblob, nBlobs] = bwlabeln(roii);
        roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
        roiinds = roiMeasurements.VoxelIdxList{1};
        
        trueind = [];
        for i = 1:numberOfBlobs
            A = blobMeasurements.VoxelIdxList{i};
            B = double(ismember(roiinds,A));
            if sum(B)>0
                trueind = horzcat(trueind, i);
            end
        end

        if isempty(trueind) ~= 1
            num_pt_clust = num_pt_clust + 1;
        end

        for i = 1:numberOfBlobs
            if ismember(i,trueind) == 1
                if detected == 0
                    detected = 1;
                end
            else
                fpth = fpth + 1;
            end
        end

        % below thresh
%         ann = load_untouch_nii('MNI_lesionprob_ANN.nii');
%         anni = single(ann.img);
%         anni(anni>=thres) = 0;
%         anni((anni<thres)&(anni>0)) = 1;
%         mask = single(load_untouch_nii('MNI_GM.nii').img) + single(load_untouch_nii('MNI_WM.nii').img);
%         mask(mask>=0.95) = 1;
%         mask(mask<0.95) = 0;
%         anni = anni.*mask;
% 
%         [labeledImage, numberOfBlobs] = bwlabeln(anni);
% 	    blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
% 	    allinds = [blobMeasurements.VoxelIdxList];
%         numberOfBlobs = size(allinds,1);
%     
%         roi = load_untouch_nii('MNI_ROI_final.nii');
%         roii = single(roi.img);
%         [roiblob, nBlobs] = bwlabeln(roii);
%         roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
%         roiinds = roiMeasurements.VoxelIdxList{1};
%         
%         trueind = [];
%         for i = 1:numberOfBlobs
%             A = blobMeasurements.VoxelIdxList{i};
%             if sum(ismember(roiinds,A))>0
%                 trueind = horzcat(trueind, i);
%             end
%         end
%     
%         for i = 1:numberOfBlobs
%             if ismember(i,trueind) == 1
%                 fnth = fnth + 1;
%             else
%                 tnth = tnth + 1;
%             end
%         end
    
        if detected == 1
            num_pt_tp = num_pt_tp + 1;
        end

%         if tnth == numberOfBlobs
%             num_pt_tn = num_pt_tn + 1;
%         end

%         if fnth > 0
%             num_pt_fn = num_pt_fn + 1;
%         end
        totalfp = totalfp + fpth;
        fpthlst = horzcat(fpthlst, fpth)
        if detected == 1
            switch fpth
                case 0
                    num_pt_fp0 = num_pt_fp0+1;
                    T{char(subtypes(sub_cursor)),1} = T{char(subtypes(sub_cursor)),1} + 1;
                case 1
                    num_pt_fp1 = num_pt_fp1+1;
                    T{char(subtypes(sub_cursor)),2} = T{char(subtypes(sub_cursor)),2} + 1;
                case 2
                    num_pt_fp2 = num_pt_fp2+1;
                    T{char(subtypes(sub_cursor)),3} = T{char(subtypes(sub_cursor)),3} + 1;
                case 3
                    num_pt_fp3 = num_pt_fp3+1;
                    T{char(subtypes(sub_cursor)),4} = T{char(subtypes(sub_cursor)),4} + 1;
                otherwise
                    num_pt_fp4 = num_pt_fp4+1;
                    T{char(subtypes(sub_cursor)),5} = T{char(subtypes(sub_cursor)),5} + 1;
        end
        
        end
    sub_cursor = sub_cursor + 1;
    end
    num_pt_clust_lst = vertcat(num_pt_clust_lst, num_pt_clust);
    tplst = vertcat(tplst, num_pt_tp);
    fplst_0 = vertcat(fplst_0, num_pt_fp0);
    fplst_1 = vertcat(fplst_1, num_pt_fp1);
    fplst_2 = vertcat(fplst_2, num_pt_fp2);
    fplst_3 = vertcat(fplst_3, num_pt_fp3);
    fplst_4 = vertcat(fplst_4, num_pt_fp4);
    fnlst = vertcat(tplst, num_pt_fn);
    tnlst = vertcat(tplst, num_pt_tn);
    thres
    T
    totalfp/37
    std(fpthlst)
end

%% ANN threshold just checking overlap
subjID=["PXX_XXXXX"];
MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';

overlaplst = [];
noollst =[];
num_pt_tp = [];

for thres = 0.5

    for p = subjID
        path = strcat(MRF_path,'\',p,'\MRF_VBM');
        cd(path)
        ann = load_untouch_nii('MNI_lesionprob_ANN_fn.nii');
        anni = single(ann.img);
        anni(anni>=thres) = 1;
        anni(anni<thres) = 0;
%         mask = single(load_untouch_nii('MNI_GM_fn.nii').img) + single(load_untouch_nii('MNI_WM_fn.nii').img);
%         mask(mask>=0.95) = 1;
%         mask(mask<0.95) = 0;
%         anni = anni.*mask;
        detected = 0;

        % above thresh
        [labeledImage, numberOfBlobs] = bwlabeln(anni);

	    blobMeasurements = regionprops3(labeledImage, 'Volume', 'VoxelIdxList', 'VoxelList');
	    allinds = [blobMeasurements.VoxelIdxList(blobMeasurements.Volume>1)];
        numberOfBlobs = size(allinds,1);
    
        roi = load_untouch_nii('MNI_ROI_final.nii');
        roii = single(roi.img);
        [roiblob, nBlobs] = bwlabeln(roii);
        roiMeasurements = regionprops3(roiblob, 'VoxelIdxList');
        roiinds = roiMeasurements.VoxelIdxList{1};
        
        trueind = [];
        for i = 1:numberOfBlobs
            A = blobMeasurements.VoxelIdxList{i};
            B = double(ismember(roiinds,A));
            if sum(B)>0
                trueind = horzcat(trueind, i);
            end
        end


        for i = 1:numberOfBlobs
            if ismember(i,trueind) == 1
                if detected == 0
                    detected = 1;
                end
            end
        end

    
        if detected == 1
            num_pt_tp = num_pt_tp + 1;
            overlaplst = horzcat(overlaplst, p);
        else
            noollst = horzcat(noollst, p);
        end

    end
end

overlaplst
noollst
num_pt_tp