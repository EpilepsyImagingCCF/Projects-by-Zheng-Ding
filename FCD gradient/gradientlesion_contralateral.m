%% separate into 30% 60% ect.
nhood = zeros(3,3,3);
nhood(2,2,1) = 1;
nhood(1,2,2) = 1;
nhood(2,1,2) = 1;
nhood(2,2,2) = 1;
nhood(2,3,2) = 1;
nhood(3,2,2) = 1;
nhood(2,2,3) = 1;

subjID = ["PXX_XXXXX"];
subjIDROI = subjID;
types = ["IIB" "IIB" "IIA"...
    "IIA" "IIB" "IIB"...
    "IIB" "IIB" "IIA" "IIB" "IIB" "IIB" ...
    "mMCD" "mMCD" ...
    "mMCD" "mMCD"...
    "mMCD" "mMCD" "mMCD"];
temporal_exc = ["P33_XXXXX" "P82_XXXXX" "P61_XXXXX"]


MRF_path='Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
myelin_path = 'Z:\Imaging\Multimodal\Myelin\Patients';
ROI_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
sz = [10 8];
varTypes = ["string","string","double","double","double","double","double","double"];
varNames = ["Subject","type","GM T1/T2","WM T1/T2","GM T1","WM T1","GM T2","WM T2"];
temps = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
centroids = [];

n2b=0;
n2a=0;
nmcd=0;

% dms = nthroot([0 30 60 100 200 300 400 500 600 800 1000]./100,3)*100;
% dms0 = nthroot([30 60 100 200 300 400 500 600 800 1000]./100,3)*100;
dms = nthroot([0 30 60 100 200 300 400 500 600 800]./100,3)*100;
dms0 = nthroot([30 60 100 200 300 400 500 600 800]./100,3)*100;
% dms = nthroot([0 25 50 100 200 330 480 640 830]./100,3)*100;
% dms0 = nthroot([25 50 100 200 330 480 640 830]./100,3)*100;
IIBGM = zeros(size(dms0'));
IIBGMt1 = zeros(size(dms0'));
IIBGMt2 = zeros(size(dms0'));
IIBWM= zeros(size(dms0'));
IIBWMt1 = zeros(size(dms0'));
IIBWMt2 = zeros(size(dms0'));
IIAGM = zeros(size(dms0'));
IIAGMt1 = zeros(size(dms0'));
IIAGMt2 = zeros(size(dms0'));
IIAWM= zeros(size(dms0'));
IIAWMt1 = zeros(size(dms0'));
IIAWMt2 = zeros(size(dms0'));
MCDGM = zeros(size(dms0'));
MCDGMt1 = zeros(size(dms0'));
MCDGMt2 = zeros(size(dms0'));
MCDWM= zeros(size(dms0'));
MCDWMt1 = zeros(size(dms0'));
MCDWMt2 = zeros(size(dms0'));

IIBGMn = zeros(size(dms0'));
IIBGMt1n = zeros(size(dms0'));
IIBGMt2n = zeros(size(dms0'));
IIBWMn= zeros(size(dms0'));
IIBWMt1n = zeros(size(dms0'));
IIBWMt2n = zeros(size(dms0'));
IIAGMn = zeros(size(dms0'));
IIAGMt1n = zeros(size(dms0'));
IIAGMt2n = zeros(size(dms0'));
IIAWMn = zeros(size(dms0'));
IIAWMt1n = zeros(size(dms0'));
IIAWMt2n = zeros(size(dms0'));
MCDGMn = zeros(size(dms0'));
MCDGMt1n = zeros(size(dms0'));
MCDGMt2n = zeros(size(dms0'));
MCDWMn = zeros(size(dms0'));
MCDWMt1n = zeros(size(dms0'));
MCDWMt2n = zeros(size(dms0'));


IIBGMse = zeros(size(dms0'));
IIBGMt1se = zeros(size(dms0'));
IIBGMt2se = zeros(size(dms0'));
IIBWMse = zeros(size(dms0'));
IIBWMt1se = zeros(size(dms0'));
IIBWMt2se = zeros(size(dms0'));
IIAGMse = zeros(size(dms0'));
IIAGMt1se = zeros(size(dms0'));
IIAGMt2se = zeros(size(dms0'));
IIAWMse= zeros(size(dms0'));
IIAWMt1se = zeros(size(dms0'));
IIAWMt2se = zeros(size(dms0'));
MCDGMse = zeros(size(dms0'));
MCDGMt1se = zeros(size(dms0'));
MCDGMt2se = zeros(size(dms0'));
MCDWMse = zeros(size(dms0'));
MCDWMt1se = zeros(size(dms0'));
MCDWMt2se = zeros(size(dms0'));
gmgs_tsz2b = zeros(size(dms0'));
wmgs_tsz2b = zeros(size(dms0'));
gmgs_tsz2a = zeros(size(dms0'));
wmgs_tsz2a = zeros(size(dms0'));
gmgs_tszmcd = zeros(size(dms0'));
wmgs_tszmcd = zeros(size(dms0'));

IIBGMsen = zeros(size(dms0'));
IIBGMt1sen = zeros(size(dms0'));
IIBGMt2sen = zeros(size(dms0'));
IIBWMsen = zeros(size(dms0'));
IIBWMt1sen = zeros(size(dms0'));
IIBWMt2sen = zeros(size(dms0'));
IIAGMsen = zeros(size(dms0'));
IIAGMt1sen = zeros(size(dms0'));
IIAGMt2sen = zeros(size(dms0'));
IIAWMsen = zeros(size(dms0'));
IIAWMt1sen = zeros(size(dms0'));
IIAWMt2sen = zeros(size(dms0'));
MCDGMsen = zeros(size(dms0'));
MCDGMt1sen = zeros(size(dms0'));
MCDGMt2sen = zeros(size(dms0'));
MCDWMsen = zeros(size(dms0'));
MCDWMt1sen = zeros(size(dms0'));
MCDWMt2sen = zeros(size(dms0'));
gmgs_tsz2bn = zeros(size(dms0'));
wmgs_tsz2bn = zeros(size(dms0'));
gmgs_tsz2an = zeros(size(dms0'));
wmgs_tsz2an = zeros(size(dms0'));
gmgs_tszmcdn = zeros(size(dms0'));
wmgs_tszmcdn = zeros(size(dms0'));


n2b=1;
n2a=1;
nmcd=1;
for p = 1:length(subjID)
    subjID(p)
    path = strcat(myelin_path,'\',subjID(p));
    
    cd(path)
%     copyfile MNI_T1oT2.nii masks
%     cd(strcat(myelin_path,'\',subjID(p),'\masks'))
%     t1ot2 = load_untouch_nii('MNI_T1oT2.nii');
%     image = single(t1ot2.img);
%     image = image.*0;
%     t1ot2.img = image;
%     save_untouch_nii(t1ot2, 'MNI_T1oT2.nii');
    
    cd(path)

    t1ot2 = load_untouch_nii('MNI_T1oT2_nmCC.nii');
    image = single(t1ot2.img);

%     t1ot2 = load_untouch_nii('MNI_T1oT2_nmCC_bico.nii');
%     image = single(t1ot2.img);

    corm = load_untouch_nii('MNI_cortical_mask.nii');
    cormi = single(corm.img);
    cormi(cormi > 0) = 1;

    cd(strcat(MRF_path,'\',subjID(p),'\MRF_VBM')); 
    t1 = load_untouch_nii('MNI_T1.nii');
    t1i = single(t1.img);
    t2 = load_untouch_nii('MNI_T2.nii');
    t2i = single(t2.img);
    roi = load_untouch_nii('MNI_ROI.nii');
    roii = single(roi.img);

    % threshold lesion label
    roii(roii >= 0.5) = 1;
    roii(roii < 0.5) = 0;
    roiprop = regionprops3(roii, 'Volume', 'VoxelIdxList', 'VoxelList');
    totalvol = roiprop.Volume;
    roiie = roii;

    gm = load_untouch_nii('MNI_GM_fn.nii');
    gmi = single(gm.img);
    gmic = gmi;
    wm = load_untouch_nii('MNI_WM_fn.nii');
    wmi = single(wm.img);
    wmic = wmi;
    brain = gmic + wmic;
    gmi((gmic>=0.95)&(brain>=0.95)) = 1;
    wmi((wmic>=0.5)&(brain>=0.95)) = 1;
    gmi((gmic>=0.95)&(brain<0.95)) = 0;
    wmi((wmic>=0.5)&(brain<0.95)) = 0;
    gmi(gmic<0.95) = 0;
    wmi(wmic<0.5) = 0;
    gmif = gmi.*flip(cormi,1);
    wmif = wmi.*flip(cormi,1);
    gmi = gmi.*cormi;
    wmi = wmi.*cormi;

 
    dmeang = [];
    dmeanw = [];
    dseg = [];
    dsew = [];
    dstdg = [];
    dstdw = [];
    dmeangt1 = [];
    dstdgt1 = [];
    dsegt1 = [];
    dmeanwt1 = [];
    dstdwt1 = [];
    dsewt1 = [];
    dmeangt2 = [];
    dstdgt2 = [];
    dsegt2 = [];
    dmeanwt2 = [];
    dstdwt2 = [];
    dsewt2 = [];

    dmeangn = [];
    dmeanwn = [];
    dsegn = [];
    dsewn = [];
    dstdgn = [];
    dstdwn = [];
    dmeangt1n = [];
    dstdgt1n = [];
    dsegt1n = [];
    dmeanwt1n = [];
    dstdwt1n = [];
    dsewt1n = [];
    dmeangt2n = [];
    dstdgt2n = [];
    dsegt2n = [];
    dmeanwt2n = [];
    dstdwt2n = [];
    dsewt2n = [];

    gpoints = [];
    wpoints = [];
    gpointst1 = [];
    wpointst1 = [];
    gpointst2 = [];
    wpointst2 = [];

    gpointsn = [];
    wpointsn = [];
    gpointst1n = [];
    wpointst1n = [];
    gpointst2n = [];
    wpointst2n = [];
    
    gr = [];
    wr = [];
    grt1 = [];
    wrt1 = [];
    grt2 = [];
    wrt2 = [];

    grn = [];
    wrn = [];
    grt1n = [];
    wrt1n = [];
    grt2n = [];
    wrt2n = [];

    percs = [];

% keep eroding the ROI until no voxel left
    dis = 0;
%     cd(strcat(myelin_path,'\',subjID(p),'\masks'))
%     copyfile('MNI_T1oT2.nii', 'the1mask.nii')
%     copyfile('MNI_T1oT2.nii', 'percentagemask.nii')

    while sum(roiie(:)) > 1
        previous = roiie;
        roiie = imerode(roiie, nhood); 
        remain = regionprops3(roiie, 'Volume', 'VoxelIdxList', 'VoxelList');
        perc = nthroot(remain.Volume/totalvol,3)*100;
        percs = vertcat(percs, perc);
        area = previous - roiie;
        dis = dis + 1;
%         if subjID(p) == "study14129"
%             area = gradientLesionZexclude(subjID(p),area);
%         end
%         if subjID(p) == "study14516II"
%             area = gradientLesionZexclude(subjID(p),area);
%         end

%         gms = gmi.*area*1;
%         wms = wmi.*area*2;
%         maskg = load_untouch_nii('the1mask.nii');
%         maskg.img = maskg.img + gms + wms;
%         save_untouch_nii(maskg, 'the1mask.nii')   
%         maskg = load_untouch_nii('percentagemask.nii');
%         maskg.img = maskg.img + gms + wms;
%         save_untouch_nii(maskg, 'percentagemask.nii')
        
        gmall = image.*gmi.*area;
        gmall = gmall((gmall > 0)&(gmall < 8));   % gm T1oT2 points within range of 1 voxel dilation
        gpoints = vertcat(gpoints, gmall);        % add to list of points
        gr = vertcat(gr, ones(size(gmall))*dis);  % record the distance indices of the group of points
        gmean = mean(gmall);                      % mean of the group of points
        gstd = std(gmall);                        % std of group of points
        gse = std(gmall)/sqrt(size(gmall, 1));    % se of group of points

        gmalln = image.*gmif.*flip(area,1);
        gmalln = gmalln((gmalln > 0)&(gmalln < 8));
        gpointsn = vertcat(gpointsn, gmalln);
        grn = vertcat(grn, ones(size(gmalln))*dis);
        gmeann = mean(gmalln);
        gstdn = std(gmalln);
        gsen = std(gmalln)/sqrt(size(gmalln, 1));

        wmall = image.*wmi.*area;
        wmall = wmall((wmall > 0)&(wmall < 8));
        wpoints = vertcat(wpoints, wmall);
        wr = vertcat(wr, ones(size(wmall))*dis);
        wmean = mean(wmall);
        wstd = std(wmall);
        wse = std(wmall)/sqrt(size(wmall, 1));

        wmalln = image.*wmif.*flip(area,1);
        wmalln = wmalln((wmalln > 0)&(wmalln < 8));
        wpointsn = vertcat(wpointsn, wmalln);
        wrn = vertcat(wrn, ones(size(wmalln))*dis);
        wmeann = mean(wmalln);
        wstdn = std(wmalln);
        wsen = std(wmalln)/sqrt(size(wmalln, 1));

        dmeang = vertcat(dmeang, gmean);        % record the mean in a list
        dstdg = vertcat(dstdg, gstd);           % record the standard dev in a list
        dseg = vertcat(dseg, gse);              % record the standard error in a list
        dmeanw = vertcat(dmeanw, wmean);            
        dstdw = vertcat(dstdw, wstd);
        dsew = vertcat(dsew, wse);

        dmeangn = vertcat(dmeangn, gmeann);
        dstdgn = vertcat(dstdgn, gstdn);
        dsegn = vertcat(dsegn, gsen);
        dmeanwn = vertcat(dmeanwn, wmeann);
        dstdwn = vertcat(dstdwn, wstdn);
        dsewn = vertcat(dsewn, wsen);

        gmallt1 = t1i.*gmi.*area;
        gmallt1 = gmallt1((gmallt1 > 600)&(gmallt1 < 1600));
        gpointst1 = vertcat(gpointst1, gmallt1);
        grt1 = vertcat(grt1, ones(size(gmallt1))*dis);
        gmeant1 = mean(gmallt1);
        gstdt1 = std(gmallt1);
        gset1 = std(gmallt1)/sqrt(size(gmallt1, 1));

        gmallt1n = t1i.*gmif.*flip(area,1);
        gmallt1n = gmallt1n((gmallt1n > 600)&(gmallt1n < 1600));
        gpointst1n = vertcat(gpointst1n, gmallt1n);
        grt1n = vertcat(grt1n, ones(size(gmallt1n))*dis);
        gmeant1n = mean(gmallt1n);
        gstdt1n = std(gmallt1n);
        gset1n = std(gmallt1n)/sqrt(size(gmallt1n, 1));

        wmallt1 = t1i.*wmi.*area;
        wmallt1 = wmallt1((wmallt1 > 600)&(wmallt1 < 1600));
        wpointst1 = vertcat(wpointst1, wmallt1);
        wrt1 = vertcat(wrt1, ones(size(wmallt1))*dis);
        wmeant1 = mean(wmallt1);
        wstdt1 = std(wmallt1);
        wset1 = std(wmallt1)/sqrt(size(wmallt1, 1));

        wmallt1n = t1i.*wmif.*flip(area,1);
        wmallt1n = wmallt1n((wmallt1n > 600)&(wmallt1n < 1600));
        wpointst1n = vertcat(wpointst1n, wmallt1n);
        wrt1n = vertcat(wrt1n, ones(size(wmallt1n))*dis);
        wmeant1n = mean(wmallt1n);
        wstdt1n = std(wmallt1n);
        wset1n = std(wmallt1n)/sqrt(size(wmallt1n, 1));

        dmeangt1 = vertcat(dmeangt1, gmeant1);
        dstdgt1 = vertcat(dstdgt1, gstdt1);
        dsegt1 = vertcat(dsegt1, gset1);
        dmeanwt1 = vertcat(dmeanwt1, wmeant1);
        dstdwt1 = vertcat(dstdwt1, wstdt1);
        dsewt1 = vertcat(dsewt1, wset1);

        dmeangt1n = vertcat(dmeangt1n, gmeant1n);
        dstdgt1n = vertcat(dstdgt1n, gstdt1n);
        dsegt1n = vertcat(dsegt1n, gset1n);
        dmeanwt1n = vertcat(dmeanwt1n, wmeant1n);
        dstdwt1n = vertcat(dstdwt1n, wstdt1n);
        dsewt1n = vertcat(dsewt1n, wset1n);

        gmallt2 = t2i.*gmi.*area;
        gmallt2 = gmallt2((gmallt2 > 20)&(gmallt2 < 60));
        gpointst2 = vertcat(gpointst2, gmallt2);
        grt2 = vertcat(grt2, ones(size(gmallt2))*dis);
        gmeant2 = mean(gmallt2);
        gstdt2 = std(gmallt2);
        gset2 = std(gmallt2)/sqrt(size(gmallt2, 1));

        gmallt2n = t2i.*gmif.*flip(area,1);
        gmallt2n = gmallt2n((gmallt2n > 20)&(gmallt2n < 60));
        gpointst2n = vertcat(gpointst2n, gmallt2n);
        grt2n = vertcat(grt2n, ones(size(gmallt2n))*dis);
        gmeant2n = mean(gmallt2n);
        gstdt2n = std(gmallt2n);
        gset2n = std(gmallt2n)/sqrt(size(gmallt2n, 1));

        wmallt2 = t2i.*wmi.*area;
        wmallt2 = wmallt2((wmallt2 > 20)&(wmallt2 < 60));
        wpointst2 = vertcat(wpointst2, wmallt2);
        wrt2 = vertcat(wrt2, ones(size(wmallt2))*dis);
        wmeant2 = mean(wmallt2);
        wstdt2 = std(wmallt2);
        wset2 = std(wmallt2)/sqrt(size(wmallt2, 1));

        wmallt2n = t2i.*wmif.*flip(area,1);
        wmallt2n = wmallt2n((wmallt2n > 20)&(wmallt2n < 60));
        wpointst2n = vertcat(wpointst2n, wmallt2n);
        wrt2n = vertcat(wrt2n, ones(size(wmallt2n))*dis);
        wmeant2n = mean(wmallt2n);
        wstdt2n = std(wmallt2n);
        wset2n = std(wmallt2n)/sqrt(size(wmallt2n, 1));

        dmeangt2 = vertcat(dmeangt2, gmeant2);
        dstdgt2 = vertcat(dstdgt2, gstdt2);
        dsegt2 = vertcat(dsegt2, gset2);
        dmeanwt2 = vertcat(dmeanwt2, wmeant2);
        dstdwt2 = vertcat(dstdwt2, wstdt2);
        dsewt2 = vertcat(dsewt2, wset2);

        dmeangt2n = vertcat(dmeangt2n, gmeant2n);
        dstdgt2n = vertcat(dstdgt2n, gstdt2n);
        dsegt2n = vertcat(dsegt2n, gset2n);
        dmeanwt2n = vertcat(dmeanwt2n, wmeant2n);
        dstdwt2n = vertcat(dstdwt2n, wstdt2n);
        dsewt2n = vertcat(dsewt2n, wset2n);        

%         cd(strcat(myelin_path,'\',subjID(p),'\masks'))
%         copyfile('MNI_T1oT2.nii', strcat('g', string(-dis),'.nii'))
%         copyfile('MNI_T1oT2.nii', strcat('r', string(-dis),'.nii'))
%         
%         maskg = load_untouch_nii(char(strcat('g', string(-dis),'.nii')))
%         maskg.img = gmif.*flip(area,1);
%         save_untouch_nii(maskg, char(strcat('g', string(-dis),'.nii')))
%         
%         maskw = load_untouch_nii(char(strcat('w', string(-dis),'.nii')))
%         maskw.img = wmif.*flip(area,1);
%         save_untouch_nii(maskw, char(strcat('w', string(-dis),'.nii')))
    end
% center of lesion
%     roiie = imerode(roiie, nhood);
    area = roiie;
    dis = dis + 1;

%     gms = gmi.*area*1;
%     wms = wmi.*area*2;
%     maskg = load_untouch_nii('the1mask.nii');
%     maskg.img = maskg.img + gms + wms;
%     save_untouch_nii(maskg, 'the1mask.nii')
%     maskg = load_untouch_nii('percentagemask.nii');
%     maskg.img = maskg.img + gms + wms;
%     save_untouch_nii(maskg, 'percentagemask.nii')

    gmall = image.*gmi.*area;
    gmall = gmall((gmall > 0)&(gmall < 8));
    gpoints = vertcat(gpoints, gmall);
    gr = vertcat(gr, ones(size(gmall))*dis);
    gmean = mean(gmall);
    gstd = std(gmall);
    gse = std(gmall)/sqrt(size(gmall, 1));

    gmalln = image.*gmif.*flip(area,1);
    gmalln = gmalln((gmalln > 0)&(gmalln < 8));
    gpointsn = vertcat(gpointsn, gmalln);
    grn = vertcat(grn, ones(size(gmalln))*dis);
    gmeann = mean(gmalln);
    gstdn = std(gmalln);
    gsen = std(gmalln)/sqrt(size(gmalln, 1));

    wmall = image.*wmi.*area;
    wmall = wmall((wmall > 0)&(wmall < 8));
    wpoints = vertcat(wpoints, wmall);
    wr = vertcat(wr, ones(size(wmall))*dis);
    wmean = mean(wmall);
    wstd = std(wmall);
    wse = std(wmall)/sqrt(size(wmall, 1));

    wmalln = image.*wmif.*flip(area,1);
    wmalln = wmalln((wmalln > 0)&(wmalln < 8));
    wpointsn = vertcat(wpointsn, wmalln);
    wrn = vertcat(wrn, ones(size(wmalln))*dis);
    wmeann = mean(wmalln);
    wstdn = std(wmalln);
    wsen = std(wmalln)/sqrt(size(wmalln, 1));

    dmeang = vertcat(dmeang, gmean);
    dstdg = vertcat(dstdg, gstd);
    dseg = vertcat(dseg, gse);
    dmeanw = vertcat(dmeanw, wmean);
    dstdw = vertcat(dstdw, wstd);
    dsew = vertcat(dsew, wse);

    dmeangn = vertcat(dmeangn, gmeann);
    dstdgn = vertcat(dstdgn, gstdn);
    dsegn = vertcat(dsegn, gsen);
    dmeanwn = vertcat(dmeanwn, wmeann);
    dstdwn = vertcat(dstdwn, wstdn);
    dsewn = vertcat(dsewn, wsen);

    gmallt1 = t1i.*gmi.*area;
    gmallt1 = gmallt1((gmallt1 > 600)&(gmallt1 < 1600));
    gpointst1 = vertcat(gpointst1, gmallt1);
    grt1 = vertcat(grt1, ones(size(gmallt1))*dis);
    gmeant1 = mean(gmallt1);
    gstdt1 = std(gmallt1);
    gset1 = std(gmallt1)/sqrt(size(gmallt1, 1));

    gmallt1n = t1i.*gmif.*flip(area,1);
    gmallt1n = gmallt1n((gmallt1n > 600)&(gmallt1n < 1600));
    gpointst1n = vertcat(gpointst1n, gmallt1n);
    grt1n = vertcat(grt1n, ones(size(gmallt1n))*dis);
    gmeant1n = mean(gmallt1n);
    gstdt1n = std(gmallt1n);
    gset1n = std(gmallt1n)/sqrt(size(gmallt1n, 1));

    wmallt1 = t1i.*wmi.*area;
    wmallt1 = wmallt1((wmallt1 > 600)&(wmallt1 < 1600));
    wpointst1 = vertcat(wpointst1, wmallt1);
    wrt1 = vertcat(wrt1, ones(size(wmallt1))*dis);
    wmeant1 = mean(wmallt1);
    wstdt1 = std(wmallt1);
    wset1 = std(wmallt1)/sqrt(size(wmallt1, 1));

    wmallt1n = t1i.*wmif.*flip(area,1);
    wmallt1n = wmallt1n((wmallt1n > 600)&(wmallt1n < 1600));
    wpointst1n = vertcat(wpointst1n, wmallt1n);
    wrt1n = vertcat(wrt1n, ones(size(wmallt1n))*dis);
    wmeant1n = mean(wmallt1n);
    wstdt1n = std(wmallt1n);
    wset1n = std(wmallt1n)/sqrt(size(wmallt1n, 1));

    dmeangt1 = vertcat(dmeangt1, gmeant1);
    dstdgt1 = vertcat(dstdgt1, gstdt1);
    dsegt1 = vertcat(dsegt1, gset1);
    dmeanwt1 = vertcat(dmeanwt1, wmeant1);
    dstdwt1 = vertcat(dstdwt1, wstdt1);
    dsewt1 = vertcat(dsewt1, wset1);

    dmeangt1n = vertcat(dmeangt1n, gmeant1n);
    dstdgt1n = vertcat(dstdgt1n, gstdt1n);
    dsegt1n = vertcat(dsegt1n, gset1n);
    dmeanwt1n = vertcat(dmeanwt1n, wmeant1n);
    dstdwt1n = vertcat(dstdwt1n, wstdt1n);
    dsewt1n = vertcat(dsewt1n, wset1n);

    gmallt2 = t2i.*gmi.*area;
    gmallt2 = gmallt2((gmallt2 > 20)&(gmallt2 < 60));
    gpointst2 = vertcat(gpointst2, gmallt2);
    grt2 = vertcat(grt2, ones(size(gmallt2))*dis);
    gmeant2 = mean(gmallt2);
    gstdt2 = std(gmallt2);
    gset2 = std(gmallt2)/sqrt(size(gmallt2, 1));

    gmallt2n = t2i.*gmif.*flip(area,1);
    gmallt2n = gmallt2n((gmallt2n > 20)&(gmallt2n < 60));
    gpointst2n = vertcat(gpointst2n, gmallt2n);
    grt2n = vertcat(grt2n, ones(size(gmallt2n))*dis);
    gmeant2n = mean(gmallt2n);
    gstdt2n = std(gmallt2n);
    gset2n = std(gmallt2n)/sqrt(size(gmallt2n, 1));

    wmallt2 = t2i.*wmi.*area;
    wmallt2 = wmallt2((wmallt2 > 20)&(wmallt2 < 60));
    wpointst2 = vertcat(wpointst2, wmallt2);
    wrt2 = vertcat(wrt2, ones(size(wmallt2))*dis);
    wmeant2 = mean(wmallt2);
    wstdt2 = std(wmallt2);
    wset2 = std(wmallt2)/sqrt(size(wmallt2, 1));

    wmallt2n = t2i.*wmif.*flip(area,1);
    wmallt2n = wmallt2n((wmallt2n > 20)&(wmallt2n < 60));
    wpointst2n = vertcat(wpointst2n, wmallt2n);
    wrt2n = vertcat(wrt2n, ones(size(wmallt2n))*dis);
    wmeant2n = mean(wmallt2n);
    wstdt2n = std(wmallt2n);
    wset2n = std(wmallt2n)/sqrt(size(wmallt2n, 1));

    dmeangt2 = vertcat(dmeangt2, gmeant2);
    dstdgt2 = vertcat(dstdgt2, gstdt2);
    dsegt2 = vertcat(dsegt2, gset2);
    dmeanwt2 = vertcat(dmeanwt2, wmeant2);
    dstdwt2 = vertcat(dstdwt2, wstdt2);
    dsewt2 = vertcat(dsewt2, wset2);
    
    dmeangt2n = vertcat(dmeangt2n, gmeant2n);
    dstdgt2n = vertcat(dstdgt2n, gstdt2n);
    dsegt2n = vertcat(dsegt2n, gset2n);
    dmeanwt2n = vertcat(dmeanwt2n, wmeant2n);
    dstdwt2n = vertcat(dstdwt2n, wstdt2n);
    dsewt2n = vertcat(dsewt2n, wset2n);
   
%     copyfile('MNI_T1oT2.nii', strcat('g', string(-dis),'.nii'))
%     copyfile('MNI_T1oT2.nii', strcat('w', string(-dis),'.nii'))
% 
%     maskg = load_untouch_nii(char(strcat('g', string(-dis),'.nii')))
%     maskg.img = gmif.*flip(area,1);
%     save_untouch_nii(maskg, char(strcat('g', string(-dis),'.nii')))
% 
%     maskw = load_untouch_nii(char(strcat('w', string(-dis),'.nii')))
%     maskw.img = wmif.*flip(area,1);
%     save_untouch_nii(maskw, char(strcat('w', string(-dis),'.nii')))

% flip
    dmeang = flip(dmeang);
    dstdg = flip(dstdg);
    dseg = flip(dseg);
    dmeanw = flip(dmeanw);
    dstdw = flip(dstdw);
    dsew = flip(dsew);
    dmeangt1 = flip(dmeangt1);
    dstdgt1 = flip(dstdgt1);
    dsegt1 = flip(dsegt1);
    dmeanwt1 = flip(dmeanwt1);
    dstdwt1 = flip(dstdwt1);
    dsewt1 = flip(dsewt1);
    dmeangt2 = flip(dmeangt2);
    dstdgt2 = flip(dstdgt2);
    dsegt2 = flip(dsegt2);
    dmeanwt2 = flip(dmeanwt2);
    dstdwt2 = flip(dstdwt2);
    dsewt2 = flip(dsewt2);
    gpoints = flip(gpoints);
    wpoints = flip(wpoints);
    gpointst1 = flip(gpointst1);
    wpointst1 = flip(wpointst1);
    gpointst2 = flip(gpointst2);
    wpointst2 = flip(wpointst2);
    gr = flip(gr);
    gr = abs(gr-dis)+1;
    wr = flip(wr);
    wr = abs(wr-dis)+1;
    grt1 = flip(grt1);
    grt1 = abs(grt1-dis)+1;
    wrt1 = flip(wrt1);
    wrt1 = abs(wrt1-dis)+1;
    grt2 = flip(grt2);
    grt2 = abs(grt2-dis)+1;
    wrt2 = flip(wrt2);
    wrt2 = abs(wrt2-dis)+1;

    dmeangn = flip(dmeangn);
    dstdgn = flip(dstdgn);
    dsegn = flip(dsegn);
    dmeanwn = flip(dmeanwn);
    dstdwn = flip(dstdwn);
    dsewn = flip(dsewn);
    dmeangt1n = flip(dmeangt1n);
    dstdgt1n = flip(dstdgt1n);
    dsegt1n = flip(dsegt1n);
    dmeanwt1n = flip(dmeanwt1n);
    dstdwt1n = flip(dstdwt1n);
    dsewt1n = flip(dsewt1n);
    dmeangt2n = flip(dmeangt2n);
    dstdgt2n = flip(dstdgt2n);
    dsegt2n = flip(dsegt2n);
    dmeanwt2n = flip(dmeanwt2n);
    dstdwt2n = flip(dstdwt2n);
    dsewt2n = flip(dsewt2n);
    gpointsn = flip(gpointsn);
    wpointsn = flip(wpointsn);
    gpointst1n = flip(gpointst1n);
    wpointst1n = flip(wpointst1n);
    gpointst2n = flip(gpointst2n);
    wpointst2n = flip(wpointst2n);
    grn = flip(grn);
    grn = abs(grn-dis)+1;
    wrn = flip(wrn);
    wrn = abs(wrn-dis)+1;
    grt1n = flip(grt1n);
    grt1n = abs(grt1n-dis)+1;
    wrt1n = flip(wrt1n);
    wrt1n = abs(wrt1n-dis)+1;
    grt2n = flip(grt2n);
    grt2n = abs(grt2n-dis)+1;
    wrt2n = flip(wrt2n);
    wrt2n = abs(wrt2n-dis)+1;

    percs = flip(percs);

    rea = 25;
    rea1 = 24;
%     if subjID(p) == "P109_16008"
%         rea = 50;
%         rea1 = 49;
%     end
    roiie = roii;
    for r = 1:rea
        previous = roiie;
        remain = regionprops3(roiie, 'Volume', 'VoxelIdxList', 'VoxelList');
        perc = nthroot(remain.Volume/totalvol,3)*100;
        percs = vertcat(percs, perc);
        roiie = gradientLesionDilate_MNI(subjID(p), roiie);
  
        area = roiie - previous;
%         area = gradientLesionZexclude_MNI(subjID(p), area);
        dis = dis + 1;
        
%         gms = gmi.*area*3;
%         wms = wmi.*area*4;
% 
%         maskg = load_untouch_nii('the1mask.nii');
%         maskg.img = maskg.img + gms + wms;
%         save_untouch_nii(maskg, 'the1mask.nii')    
% 
%         if perc <= 250
%             maskg = load_untouch_nii('percentagemask.nii');
%             maskg.img = maskg.img + gms + wms;
%             save_untouch_nii(maskg, 'percentagemask.nii')
%         end

        gmall = image.*gmi.*area;
        gmall = gmall((gmall > 0)&(gmall < 8));
        gpoints = vertcat(gpoints, gmall);
        gr = vertcat(gr, ones(size(gmall))*dis);
        gmean = mean(gmall);
        gstd = std(gmall);
        gse = std(gmall)/sqrt(size(gmall, 1));

        gmalln = image.*gmif.*flip(area,1);
        gmalln = gmalln((gmalln > 0)&(gmalln < 8));
        gpointsn = vertcat(gpointsn, gmalln);
        grn = vertcat(grn, ones(size(gmalln))*dis);
        gmeann = mean(gmalln);
        gstdn = std(gmalln);
        gsen = std(gmalln)/sqrt(size(gmalln, 1));

        wmall = image.*wmi.*area;
        wmall = wmall((wmall > 0)&(wmall < 8));
        wpoints = vertcat(wpoints, wmall);
        wr = vertcat(wr, ones(size(wmall))*dis);
        wmean = mean(wmall);
        wstd = std(wmall);
        wse = std(wmall)/sqrt(size(wmall, 1));

        wmalln = image.*wmif.*flip(area,1);
        wmalln = wmalln((wmalln > 0)&(wmalln < 8));
        wpointsn = vertcat(wpointsn, wmalln);
        wrn = vertcat(wrn, ones(size(wmalln))*dis);
        wmeann = mean(wmalln);
        wstdn = std(wmalln);
        wsen = std(wmalln)/sqrt(size(wmalln, 1));

        dmeang = vertcat(dmeang, gmean);
        dstdg = vertcat(dstdg, gstd);
        dseg = vertcat(dseg, gse);
        dmeanw = vertcat(dmeanw, wmean);
        dstdw = vertcat(dstdw, wstd);
        dsew = vertcat(dsew, wse);

        dmeangn = vertcat(dmeangn, gmeann);
        dstdgn = vertcat(dstdgn, gstdn);
        dsegn = vertcat(dsegn, gsen);
        dmeanwn = vertcat(dmeanwn, wmeann);
        dstdwn = vertcat(dstdwn, wstdn);
        dsewn = vertcat(dsewn, wsen);

        gmallt1 = t1i.*gmi.*area;
        gmallt1 = gmallt1((gmallt1 > 600)&(gmallt1 < 1600));
        gpointst1 = vertcat(gpointst1, gmallt1);
        grt1 = vertcat(grt1, ones(size(gmallt1))*dis);
        gmeant1 = mean(gmallt1);
        gstdt1 = std(gmallt1);
        gset1 = std(gmallt1)/sqrt(size(gmallt1, 1));

        gmallt1n = t1i.*gmif.*flip(area,1);
        gmallt1n = gmallt1n((gmallt1n > 600)&(gmallt1n < 1600));
        gpointst1n = vertcat(gpointst1n, gmallt1n);
        grt1n = vertcat(grt1n, ones(size(gmallt1n))*dis);
        gmeant1n = mean(gmallt1n);
        gstdt1n = std(gmallt1n);
        gset1n = std(gmallt1n)/sqrt(size(gmallt1n, 1));

        wmallt1 = t1i.*wmi.*area;
        wmallt1 = wmallt1((wmallt1 > 600)&(wmallt1 < 1600));
        wpointst1 = vertcat(wpointst1, wmallt1);
        wrt1 = vertcat(wrt1, ones(size(wmallt1))*dis);
        wmeant1 = mean(wmallt1);
        wstdt1 = std(wmallt1);
        wset1 = std(wmallt1)/sqrt(size(wmallt1, 1));

        wmallt1n = t1i.*wmif.*flip(area,1);
        wmallt1n = wmallt1n((wmallt1n > 600)&(wmallt1n < 1600));
        wpointst1n = vertcat(wpointst1n, wmallt1n);
        wrt1n = vertcat(wrt1n, ones(size(wmallt1n))*dis);
        wmeant1n = mean(wmallt1n);
        wstdt1n = std(wmallt1n);
        wset1n = std(wmallt1n)/sqrt(size(wmallt1n, 1));

        dmeangt1 = vertcat(dmeangt1, gmeant1);
        dstdgt1 = vertcat(dstdgt1, gstdt1);
        dsegt1 = vertcat(dsegt1, gset1);
        dmeanwt1 = vertcat(dmeanwt1, wmeant1);
        dstdwt1 = vertcat(dstdwt1, wstdt1);
        dsewt1 = vertcat(dsewt1, wset1);

        dmeangt1n = vertcat(dmeangt1n, gmeant1n);
        dstdgt1n = vertcat(dstdgt1n, gstdt1n);
        dsegt1n = vertcat(dsegt1n, gset1n);
        dmeanwt1n = vertcat(dmeanwt1n, wmeant1n);
        dstdwt1n = vertcat(dstdwt1n, wstdt1n);
        dsewt1n = vertcat(dsewt1n, wset1n);

        gmallt2 = t2i.*gmi.*area;
        gmallt2 = gmallt2((gmallt2 > 20)&(gmallt2 < 60));
        gpointst2 = vertcat(gpointst2, gmallt2);
        grt2 = vertcat(grt2, ones(size(gmallt2))*dis);
        gmeant2 = mean(gmallt2);
        gstdt2 = std(gmallt2);
        gset2 = std(gmallt2)/sqrt(size(gmallt2, 1));

        gmallt2n = t2i.*gmif.*flip(area,1);
        gmallt2n = gmallt2n((gmallt2n > 20)&(gmallt2n < 60));
        gpointst2n = vertcat(gpointst2n, gmallt2n);
        grt2n = vertcat(grt2n, ones(size(gmallt2n))*dis);
        gmeant2n = mean(gmallt2n);
        gstdt2n = std(gmallt2n);
        gset2n = std(gmallt2n)/sqrt(size(gmallt2n, 1));

        wmallt2 = t2i.*wmi.*area;
        wmallt2 = wmallt2((wmallt2 > 20)&(wmallt2 < 60));
        wpointst2 = vertcat(wpointst2, wmallt2);
        wrt2 = vertcat(wrt2, ones(size(wmallt2))*dis);
        wmeant2 = mean(wmallt2);
        wstdt2 = std(wmallt2);
        wset2 = std(wmallt2)/sqrt(size(wmallt2, 1));

        wmallt2n = t2i.*wmif.*flip(area,1);
        wmallt2n = wmallt2n((wmallt2n > 20)&(wmallt2n < 60));
        wpointst2n = vertcat(wpointst2n, wmallt2n);
        wrt2n = vertcat(wrt2n, ones(size(wmallt2n))*dis);
        wmeant2n = mean(wmallt2n);
        wstdt2n = std(wmallt2n);
        wset2n = std(wmallt2n)/sqrt(size(wmallt2n, 1));

        dmeangt2 = vertcat(dmeangt2, gmeant2);
        dstdgt2 = vertcat(dstdgt2, gstdt2);
        dsegt2 = vertcat(dsegt2, gset2);
        dmeanwt2 = vertcat(dmeanwt2, wmeant2);
        dstdwt2 = vertcat(dstdwt2, wstdt2);
        dsewt2 = vertcat(dsewt2, wset2);

        dmeangt2n = vertcat(dmeangt2n, gmeant2n);
        dstdgt2n = vertcat(dstdgt2n, gstdt2n);
        dsegt2n = vertcat(dsegt2n, gset2n);
        dmeanwt2n = vertcat(dmeanwt2n, wmeant2n);
        dstdwt2n = vertcat(dstdwt2n, wstdt2n);
        dsewt2n = vertcat(dsewt2n, wset2n);

%         cd(strcat(myelin_path,'\',subjID(p),'\masks'))
%         copyfile('MNI_T1oT2.nii', char(strcat('g', string(dis),'.nii')))
%         copyfile('MNI_T1oT2.nii', char(strcat('w', string(dis),'.nii')))
%         
%         maskg = load_untouch_nii(char(strcat('g', string(dis),'.nii')))
%         maskg.img = gmif.*flip(area,1);
%         save_untouch_nii(maskg, char(strcat('g', string(dis),'.nii')))
%         
%         maskw = load_untouch_nii(char(strcat('w', string(dis),'.nii')))
%         maskw.img = wmif.*flip(area,1);
%         save_untouch_nii(maskw, char(strcat('w', string(dis),'.nii')))    
    end
%     dvarg = dstdg./dmeang;
%     dvarw = dstdw./dmeanw;
    gmgs = [];
    wmgs = [];
    gt1gs = [];
    wt1gs = [];
    gt2gs = [];
    wt2gs = [];
    gmgs_mn = [];
    wmgs_mn = [];
    gt1gs_mn = [];
    wt1gs_mn = [];
    gt2gs_mn = [];
    wt2gs_mn = [];
    gmgs_se = [];
    wmgs_se = [];
    gt1gs_se = [];
    wt1gs_se = [];
    gt2gs_se = [];
    wt2gs_se = [];
    gmgs_sz = [];
    wmgs_sz = [];
    gt1gs_sz = [];
    wt1gs_sz = [];
    gt2gs_sz = [];
    wt2gs_sz = [];

    gmgsn = [];
    wmgsn = [];
    gt1gsn = [];
    wt1gsn = [];
    gt2gsn = [];
    wt2gsn = [];
    gmgs_mnn = [];
    wmgs_mnn = [];
    gt1gs_mnn = [];
    wt1gs_mnn = [];
    gt2gs_mnn = [];
    wt2gs_mnn = [];
    gmgs_sen = [];
    wmgs_sen = [];
    gt1gs_sen = [];
    wt1gs_sen = [];
    gt2gs_sen = [];
    wt2gs_sen = [];
    gmgs_szn = [];
    wmgs_szn = [];
    gt1gs_szn = [];
    wt1gs_szn = [];
    gt2gs_szn = [];
    wt2gs_szn = [];

    percs = percs(1:end);
    dislist = linspace(1, dis-2, dis-2);

    for i = 1:size(dms,2)-1
        dises = dislist((percs>dms(i))&(percs<=dms(i+1)));
        gmg = gpoints((gr>=dises(1))&(gr<=dises(end)));
        wmg = wpoints((wr>=dises(1))&(wr<=dises(end)));
        gt1g = gpointst1((grt1>=dises(1))&(grt1<=dises(end)));
        wt1g = wpointst1((wrt1>=dises(1))&(wrt1<=dises(end)));
        gt2g = gpointst2((grt2>=dises(1))&(grt2<=dises(end)));
        wt2g = wpointst2((wrt2>=dises(1))&(wrt2<=dises(end)));
        
        gmgs_mn = vertcat(gmgs_mn, mean(gmg));
        wmgs_mn = vertcat(wmgs_mn, mean(wmg));
        gt1gs_mn = vertcat(gt1gs_mn, mean(gt1g));
        wt1gs_mn = vertcat(wt1gs_mn, mean(wt1g));
        gt2gs_mn = vertcat(gt2gs_mn, mean(gt2g));
        wt2gs_mn = vertcat(wt2gs_mn, mean(wt2g));

        gmgs_se = vertcat(gmgs_se, std(gmg));
        wmgs_se = vertcat(wmgs_se, std(wmg));
        gt1gs_se = vertcat(gt1gs_se, std(gt1g));
        wt1gs_se = vertcat(wt1gs_se, std(wt1g));
        gt2gs_se = vertcat(gt2gs_se, std(gt2g));
        wt2gs_se = vertcat(wt2gs_se, std(wt2g));
        
        gmgs_sz = vertcat(gmgs_sz, size(gmg,1));
        wmgs_sz = vertcat(wmgs_sz, size(wmg,1));
        gt1gs_sz = vertcat(gt1gs_sz, size(gt1g,1));
        wt1gs_sz = vertcat(wt1gs_sz, size(wt1g,1));
        gt2gs_sz = vertcat(gt2gs_sz, size(gt2g,1));
        wt2gs_sz = vertcat(wt2gs_sz, size(wt2g,1));
     
        gmgn = gpointsn((grn>=dises(1))&(grn<=dises(end)));
        wmgn = wpointsn((wrn>=dises(1))&(wrn<=dises(end)));
        gt1gn = gpointst1n((grt1n>=dises(1))&(grt1n<=dises(end)));
        wt1gn = wpointst1n((wrt1n>=dises(1))&(wrt1n<=dises(end)));
        gt2gn = gpointst2n((grt2n>=dises(1))&(grt2n<=dises(end)));
        wt2gn = wpointst2n((wrt2n>=dises(1))&(wrt2n<=dises(end)));
        
        gmgs_mnn = vertcat(gmgs_mnn, mean(gmgn));
        wmgs_mnn = vertcat(wmgs_mnn, mean(wmgn));
        gt1gs_mnn = vertcat(gt1gs_mnn, mean(gt1gn));
        wt1gs_mnn = vertcat(wt1gs_mnn, mean(wt1gn));
        gt2gs_mnn = vertcat(gt2gs_mnn, mean(gt2gn));
        wt2gs_mnn = vertcat(wt2gs_mnn, mean(wt2gn));

        gmgs_sen = vertcat(gmgs_sen, std(gmgn));
        wmgs_sen = vertcat(wmgs_sen, std(wmgn));
        gt1gs_sen = vertcat(gt1gs_sen, std(gt1gn));
        wt1gs_sen = vertcat(wt1gs_sen, std(wt1gn));
        gt2gs_sen = vertcat(gt2gs_sen, std(gt2gn));
        wt2gs_sen = vertcat(wt2gs_sen, std(wt2gn));
        
        gmgs_szn = vertcat(gmgs_szn, size(gmgn,1));
        wmgs_szn = vertcat(wmgs_szn, size(wmgn,1));
        gt1gs_szn = vertcat(gt1gs_szn, size(gt1gn,1));
        wt1gs_szn = vertcat(wt1gs_szn, size(wt1gn,1));
        gt2gs_szn = vertcat(gt2gs_szn, size(gt2gn,1));
        wt2gs_szn = vertcat(wt2gs_szn, size(wt2gn,1));
    end
% switch lesize(p)
%     case "S"
%         gmgs_sz(7:end) = 0;
%         wmgs_sz(7:end) = 0;
% end   
% switch types(p) 
%     case "IIB"
%         IIBGM= horzcat(IIBGM, gmgs_mn/(gmgs_mn(end)));
%         IIBGMt1 = horzcat(IIBGMt1, gt1gs_mn/(gt1gs_mn(end)));
%         IIBGMt2 = horzcat(IIBGMt2, gt2gs_mn/(gt2gs_mn(end)));
%         IIBWM= horzcat(IIBWM, wmgs_mn/(wmgs_mn(end)));
%         IIBWMt1 = horzcat(IIBWMt1, wt1gs_mn/(wt1gs_mn(end)));
%         IIBWMt2 = horzcat(IIBWMt2, wt2gs_mn/(wt2gs_mn(end)));
%     case "IIA"
%         IIAGM= horzcat(IIAGM, gmgs_mn/(gmgs_mn(end)));
%         IIAGMt1 = horzcat(IIAGMt1, gt1gs_mn/(gt1gs_mn(end)));
%         IIAGMt2 = horzcat(IIAGMt2, gt2gs_mn/(gt2gs_mn(end)));
%         IIAWM= horzcat(IIAWM, wmgs_mn/(wmgs_mn(end)));
%         IIAWMt1 = horzcat(IIAWMt1, wt1gs_mn/(wt1gs_mn(end)));
%         IIAWMt2 = horzcat(IIAWMt2, wt2gs_mn/(wt2gs_mn(end)));
%     case "MCD"
%         MCDGM= horzcat(MCDGM, gmgs_mn/(gmgs_mn(end)));
%         MCDGMt1 = horzcat(MCDGMt1, gt1gs_mn/(gt1gs_mn(end)));
%         MCDGMt2 = horzcat(MCDGMt2, gt2gs_mn/(gt2gs_mn(end)));
%         MCDWM= horzcat(MCDWM, wmgs_mn/(wmgs_mn(end)));
%         MCDWMt1 = horzcat(MCDWMt1, wt1gs_mn/(wt1gs_mn(end)));
%         MCDWMt2 = horzcat(MCDWMt2, wt2gs_mn/(wt2gs_mn(end)));
% end
% by max-min
% sqrt(((n1-1)*s1*s1 + (n2-1)*s2*s2 + n1 * n2 / (n1 + n2) * (m1*m1 + m2*m2 - 2 * m1 * m2)) / (n1 + n2 -1))

% separate plots 
    figure()
    subplot(2,3,1)
    shadedErrorBar(dms(2:5),gmgs_mn(1:4), gmgs_se(1:4)./sqrt(gmgs_sz(1:4)),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),gmgs_mn(4:end), gmgs_se(4:end)./sqrt(gmgs_sz(4:end)),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(2:5),gmgs_mnn(1:4), gmgs_sen(1:4)./sqrt(gmgs_szn(1:4)),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),gmgs_mnn(4:end), gmgs_sen(4:end)./sqrt(gmgs_szn(4:end)),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
    title('GM - by ROI dilation and erosion', 'FontSize', 20)
    subtitle(strcat(subjID(p)," ", types(p)), 'FontSize', 18)
    xlabel('Percentage Distance from Lesion Center', 'FontSize', 20)
    ylabel('T1/T2 value', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;

    subplot(2,3,4)
%     errorbar(dms(2:4),wmgs_mn(1:3),wmgs_se(1:3)./sqrt(wmgs_sz(1:3)), "Color", "#A30234")
%     hold on
%     errorbar(dms(5:end),wmgs_mn(4:end),wmgs_se(4:end)./sqrt(wmgs_sz(4:end)), "Color", "#ce8080")
%     hold on
%     errorbar(dms(2:4),wmgs_mnn(1:3),wmgs_sen(1:3)./sqrt(wmgs_szn(1:3)), "Color", "#002157")
%     hold on
%     errorbar(dms(5:end),wmgs_mnn(4:end),wmgs_sen(4:end)./sqrt(wmgs_szn(4:end)), "Color", "#0076c0")
    shadedErrorBar(dms(2:5),wmgs_mn(1:4), wmgs_se(1:4)./sqrt(wmgs_sz(1:4)),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),wmgs_mn(4:end), wmgs_se(4:end)./sqrt(wmgs_sz(4:end)),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(2:5),wmgs_mnn(1:4), wmgs_sen(1:4)./sqrt(wmgs_szn(1:4)),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),wmgs_mnn(4:end), wmgs_sen(4:end)./sqrt(wmgs_szn(4:end)),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
    title('WM - by ROI dilation and erosion', 'FontSize', 20)
    subtitle(strcat(subjID(p)," ", types(p)), 'FontSize', 18)
    xlabel('Percentage Distance from Lesion Center', 'FontSize', 20)
    ylabel('T1/T2 value', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;

%     subplot(2,4,4)
%     scatter([1:1:dis],dvarg)
%     title('GM CV')
%     xlabel('distance from center (0.94mm)')
%     ylabel('T1/T2 CV')
% 
%     subplot(2,4,8)
%     scatter([1:1:dis],dvarw)
%     title('WM CV')
%     xlabel('distance from center (0.94mm)')
%     ylabel('T1/T2 CV')

    subplot(2,3,2)
%     errorbar(dms(2:4), gt1gs_mn(1:3), gt1gs_se(1:3)./sqrt(gt1gs_sz(1:3)), "Color", "#A30234")
%     hold on
%     errorbar(dms(5:end), gt1gs_mn(4:end), gt1gs_se(4:end)./sqrt(gt1gs_sz(4:end)), "Color", "#ce8080")    
%     hold on
%     errorbar(dms(2:4), gt1gs_mnn(1:3), gt1gs_sen(1:3)./sqrt(gt1gs_szn(1:3)), "Color", "#002157")
%     hold on
%     errorbar(dms(5:end), gt1gs_mnn(4:end), gt1gs_sen(4:end)./sqrt(gt1gs_szn(4:end)), "Color", "#0076c0")
    shadedErrorBar(dms(2:5),gt1gs_mn(1:4), gt1gs_se(1:4)./sqrt(gt1gs_sz(1:4)),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),gt1gs_mn(4:end), gt1gs_se(4:end)./sqrt(gt1gs_sz(4:end)),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(2:5),gt1gs_mnn(1:4), gt1gs_sen(1:4)./sqrt(gt1gs_szn(1:4)),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),gt1gs_mnn(4:end), gt1gs_sen(4:end)./sqrt(gt1gs_szn(4:end)),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
    title('GM T1', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('Percentage Distance from Lesion Center', 'FontSize', 20)
    ylabel('T1 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;

    subplot(2,3,5)
%     errorbar(dms(2:4),wt1gs_mn(1:3),wt1gs_se(1:3)./sqrt(wt1gs_sz(1:3)), "Color", "#A30234")
%     hold on
%     errorbar(dms(5:end),wt1gs_mn(4:end),wt1gs_se(4:end)./sqrt(wt1gs_sz(4:end)), "Color", "#ce8080")
%     hold on
%     errorbar(dms(2:4),wt1gs_mnn(1:3),wt1gs_sen(1:3)./sqrt(wt1gs_szn(1:3)), "Color", "#002157")
%     hold on
%     errorbar(dms(5:end),wt1gs_mnn(4:end),wt1gs_sen(4:end)./sqrt(wt1gs_szn(4:end)), "Color", "#0076c0")
    shadedErrorBar(dms(2:5),wt1gs_mn(1:4), wt1gs_se(1:4)./sqrt(wt1gs_sz(1:4)),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),wt1gs_mn(4:end), wt1gs_se(4:end)./sqrt(wt1gs_sz(4:end)),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(2:5),wt1gs_mnn(1:4), wt1gs_sen(1:4)./sqrt(wt1gs_szn(1:4)),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),wt1gs_mnn(4:end), wt1gs_sen(4:end)./sqrt(wt1gs_szn(4:end)),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
    title('WM T1', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('Percentage Distance from Lesion Center', 'FontSize', 20)
    ylabel('T1 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;

    subplot(2,3,3)
%     errorbar(dms(2:4),gt2gs_mn(1:3),gt2gs_se(1:3)./sqrt(gt2gs_sz(1:3)), "Color", "#A30234")
%     hold on
%     errorbar(dms(5:end),gt2gs_mn(4:end),gt2gs_se(4:end)./sqrt(gt2gs_sz(4:end)), "Color", "#ce8080")
%     hold on
%     errorbar(dms(2:4),gt2gs_mnn(1:3),gt2gs_sen(1:3)./sqrt(gt2gs_szn(1:3)), "Color", "#002157")
%     hold on
%     errorbar(dms(5:end),gt2gs_mnn(4:end),gt2gs_sen(4:end)./sqrt(gt2gs_szn(4:end)), "Color", "#0076c0")
    shadedErrorBar(dms(2:5),gt2gs_mn(1:4), gt2gs_se(1:4)./sqrt(gt2gs_sz(1:4)),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),gt2gs_mn(4:end), gt2gs_se(4:end)./sqrt(gt2gs_sz(4:end)),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(2:5),gt2gs_mnn(1:4), gt2gs_sen(1:4)./sqrt(gt2gs_szn(1:4)),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),gt2gs_mnn(4:end), gt2gs_sen(4:end)./sqrt(gt2gs_szn(4:end)),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
    title('GM T2', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('Percentage Distance from Lesion Center', 'FontSize', 20)
    ylabel('T2 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;

    subplot(2,3,6)
%     errorbar(dms(2:4),wt2gs_mn(1:3),wt2gs_se(1:3)./sqrt(wt2gs_sz(1:3)), "Color", "#A30234")
%     hold on
%     errorbar(dms(5:end),wt2gs_mn(4:end),wt2gs_se(4:end)./sqrt(wt2gs_sz(4:end)), "Color", "#ce8080")
%     hold on
%     errorbar(dms(2:4),wt2gs_mnn(1:3),wt2gs_sen(1:3)./sqrt(wt2gs_szn(1:3)), "Color", "#002157")
%     hold on
%     errorbar(dms(5:end),wt2gs_mnn(4:end),wt2gs_sen(4:end)./sqrt(wt2gs_szn(4:end)), "Color", "#0076c0")
    shadedErrorBar(dms(2:5),wt2gs_mn(1:4), wt2gs_se(1:4)./sqrt(wt2gs_sz(1:4)),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),wt2gs_mn(4:end), wt2gs_se(4:end)./sqrt(wt2gs_sz(4:end)),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(2:5),wt2gs_mnn(1:4), wt2gs_sen(1:4)./sqrt(wt2gs_szn(1:4)),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
    hold on
    shadedErrorBar(dms(5:end),wt2gs_mnn(4:end), wt2gs_sen(4:end)./sqrt(wt2gs_szn(4:end)),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
    title('WM T2', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('Percentage Distance from Lesion Center', 'FontSize', 20)
    ylabel('T2 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
   

    switch types(p) 
        case "IIB"
            n2b = n2b+1;
            oma = IIBGM;
            omb = IIBGMt1;
            omc = IIBGMt2;
            omd = IIBWM;
            ome = IIBWMt1;
            omf = IIBWMt2;
            IIBGM = (IIBGM.*gmgs_tsz2b + (gmgs_mn).*gmgs_sz)./(gmgs_tsz2b + gmgs_sz);
            IIBGMt1 = (IIBGMt1.*gmgs_tsz2b + (gt1gs_mn).*gmgs_sz)./(gmgs_tsz2b + gmgs_sz);
            IIBGMt2 = (IIBGMt2.*gmgs_tsz2b + (gt2gs_mn).*gmgs_sz)./(gmgs_tsz2b + gmgs_sz);
            IIBWM= (IIBWM.*wmgs_tsz2b + (wmgs_mn).*wmgs_sz)./(wmgs_tsz2b + wmgs_sz);
            IIBWMt1 = (IIBWMt1.*wmgs_tsz2b + (wt1gs_mn).*wmgs_sz)./(wmgs_tsz2b + wmgs_sz);
            IIBWMt2 = (IIBWMt2.*wmgs_tsz2b + (wt2gs_mn).*wmgs_sz)./(wmgs_tsz2b + wmgs_sz);
            IIBGMse = sqrt(((gmgs_tsz2b-1).*IIBGMse.^2 + (gmgs_sz-1).*gmgs_se.^2 + gmgs_tsz2b.*gmgs_sz./(gmgs_tsz2b+gmgs_sz).*(oma.^2+IIBGM.^2-2*oma.*IIBGM))./(gmgs_tsz2b+gmgs_sz-1));
            IIBGMt1se = sqrt(((gmgs_tsz2b-1).*IIBGMt1se.^2 + (gmgs_sz-1).*gt1gs_se.^2 + gmgs_tsz2b.*gmgs_sz./(gmgs_tsz2b+gmgs_sz).*(omb.^2+IIBGMt1.^2-2*omb.*IIBGMt1))./(gmgs_tsz2b+gmgs_sz-1));
            IIBGMt2se = sqrt(((gmgs_tsz2b-1).*IIBGMt2se.^2 + (gmgs_sz-1).*gt2gs_se.^2 + gmgs_tsz2b.*gmgs_sz./(gmgs_tsz2b+gmgs_sz).*(omc.^2+IIBGMt2.^2-2*omc.*IIBGMt2))./(gmgs_tsz2b+gmgs_sz-1));
            IIBWMse = sqrt(((wmgs_tsz2b-1).*IIBWMse.^2 + (wmgs_sz-1).*wmgs_se.^2 + wmgs_tsz2b.*wmgs_sz./(wmgs_tsz2b+wmgs_sz).*(omd.^2+IIBWM.^2-2*omd.*IIBWM))./(wmgs_tsz2b+wmgs_sz-1));
            IIBWMt1se = sqrt(((wmgs_tsz2b-1).*IIBWMt1se.^2 + (wmgs_sz-1).*wt1gs_se.^2 + wmgs_tsz2b.*wmgs_sz./(wmgs_tsz2b+wmgs_sz).*(ome.^2+IIBWMt1.^2-2*ome.*IIBWMt1))./(wmgs_tsz2b+wmgs_sz-1));
            IIBWMt2se = sqrt(((wmgs_tsz2b-1).*IIBWMt2se.^2 + (wmgs_sz-1).*wt2gs_se.^2 + wmgs_tsz2b.*wmgs_sz./(wmgs_tsz2b+wmgs_sz).*(omf.^2+IIBWMt2.^2-2*omf.*IIBWMt2))./(wmgs_tsz2b+wmgs_sz-1));
            gmgs_tsz2b = gmgs_tsz2b + gmgs_sz;
            wmgs_tsz2b = wmgs_tsz2b + wmgs_sz;
            
            

            oman = IIBGMn;
            ombn = IIBGMt1n;
            omcn = IIBGMt2n;
            omdn = IIBWMn;
            omen = IIBWMt1n;
            omfn = IIBWMt2n;
            IIBGMn = (IIBGMn.*gmgs_tsz2bn + (gmgs_mnn).*gmgs_szn)./(gmgs_tsz2bn + gmgs_szn);
            IIBGMt1n = (IIBGMt1n.*gmgs_tsz2bn + (gt1gs_mnn).*gmgs_szn)./(gmgs_tsz2bn + gmgs_szn);
            IIBGMt2n = (IIBGMt2n.*gmgs_tsz2bn + (gt2gs_mnn).*gmgs_szn)./(gmgs_tsz2bn + gmgs_szn);
            IIBWMn= (IIBWMn.*wmgs_tsz2bn + (wmgs_mnn).*wmgs_szn)./(wmgs_tsz2bn + wmgs_szn);
            IIBWMt1n = (IIBWMt1n.*wmgs_tsz2bn + (wt1gs_mnn).*wmgs_szn)./(wmgs_tsz2bn + wmgs_szn);
            IIBWMt2n = (IIBWMt2n.*wmgs_tsz2bn + (wt2gs_mnn).*wmgs_szn)./(wmgs_tsz2bn + wmgs_szn);
            IIBGMsen = sqrt(((gmgs_tsz2bn-1).*IIBGMsen.^2 + (gmgs_szn-1).*gmgs_sen.^2 + gmgs_tsz2bn.*gmgs_szn./(gmgs_tsz2bn+gmgs_szn).*(oman.^2+IIBGMn.^2-2*oman.*IIBGMn))./(gmgs_tsz2bn+gmgs_szn-1));
            IIBGMt1sen = sqrt(((gmgs_tsz2bn-1).*IIBGMt1sen.^2 + (gmgs_szn-1).*gt1gs_sen.^2 + gmgs_tsz2bn.*gmgs_szn./(gmgs_tsz2bn+gmgs_szn).*(ombn.^2+IIBGMt1n.^2-2*ombn.*IIBGMt1n))./(gmgs_tsz2bn+gmgs_szn-1));
            IIBGMt2sen = sqrt(((gmgs_tsz2bn-1).*IIBGMt2sen.^2 + (gmgs_szn-1).*gt2gs_sen.^2 + gmgs_tsz2bn.*gmgs_szn./(gmgs_tsz2bn+gmgs_szn).*(omcn.^2+IIBGMt2n.^2-2*omcn.*IIBGMt2n))./(gmgs_tsz2bn+gmgs_szn-1));
            IIBWMsen = sqrt(((wmgs_tsz2bn-1).*IIBWMsen.^2 + (wmgs_szn-1).*wmgs_sen.^2 + wmgs_tsz2bn.*wmgs_szn./(wmgs_tsz2bn+wmgs_szn).*(omdn.^2+IIBWMn.^2-2*omdn.*IIBWMn))./(wmgs_tsz2bn+wmgs_szn-1));
            IIBWMt1sen = sqrt(((wmgs_tsz2bn-1).*IIBWMt1sen.^2 + (wmgs_szn-1).*wt1gs_sen.^2 + wmgs_tsz2bn.*wmgs_szn./(wmgs_tsz2bn+wmgs_szn).*(omen.^2+IIBWMt1n.^2-2*omen.*IIBWMt1n))./(wmgs_tsz2bn+wmgs_szn-1));
            IIBWMt2sen = sqrt(((wmgs_tsz2bn-1).*IIBWMt2sen.^2 + (wmgs_szn-1).*wt2gs_sen.^2 + wmgs_tsz2bn.*wmgs_szn./(wmgs_tsz2bn+wmgs_szn).*(omfn.^2+IIBWMt2n.^2-2*omfn.*IIBWMt2n))./(wmgs_tsz2bn+wmgs_szn-1));
            gmgs_tsz2bn = gmgs_tsz2bn + gmgs_szn;
            wmgs_tsz2bn = wmgs_tsz2bn + wmgs_szn;
    
        case "IIA"
            n2a = n2a+1;
            oma = IIAGM;
            omb = IIAGMt1;
            omc = IIAGMt2;
            omd = IIAWM;
            ome = IIAWMt1;
            omf = IIAWMt2;
            IIAGM = (IIAGM.*gmgs_tsz2a + (gmgs_mn).*gmgs_sz)./(gmgs_tsz2a + gmgs_sz);
            IIAGMt1 = (IIAGMt1.*gmgs_tsz2a + (gt1gs_mn).*gmgs_sz)./(gmgs_tsz2a + gmgs_sz);
            IIAGMt2 = (IIAGMt2.*gmgs_tsz2a + (gt2gs_mn).*gmgs_sz)./(gmgs_tsz2a + gmgs_sz);
            IIAWM= (IIAWM.*wmgs_tsz2a + (wmgs_mn).*wmgs_sz)./(wmgs_tsz2a + wmgs_sz);
            IIAWMt1 = (IIAWMt1.*wmgs_tsz2a + (wt1gs_mn).*wmgs_sz)./(wmgs_tsz2a + wmgs_sz);
            IIAWMt2 = (IIAWMt2.*wmgs_tsz2a + (wt2gs_mn).*wmgs_sz)./(wmgs_tsz2a + wmgs_sz);
            IIAGMse = sqrt(((gmgs_tsz2a-1).*IIAGMse.^2 + (gmgs_sz-1).*gmgs_se.^2 + gmgs_tsz2a.*gmgs_sz./(gmgs_tsz2a+gmgs_sz).*(oma.^2+IIAGM.^2-2*oma.*IIAGM))./(gmgs_tsz2a+gmgs_sz-1));
            IIAGMt1se = sqrt(((gmgs_tsz2a-1).*IIAGMt1se.^2 + (gmgs_sz-1).*gt1gs_se.^2 + gmgs_tsz2a.*gmgs_sz./(gmgs_tsz2a+gmgs_sz).*(omb.^2+IIAGMt1.^2-2*omb.*IIAGMt1))./(gmgs_tsz2a+gmgs_sz-1));
            IIAGMt2se = sqrt(((gmgs_tsz2a-1).*IIAGMt2se.^2 + (gmgs_sz-1).*gt2gs_se.^2 + gmgs_tsz2a.*gmgs_sz./(gmgs_tsz2a+gmgs_sz).*(omc.^2+IIAGMt2.^2-2*omc.*IIAGMt2))./(gmgs_tsz2a+gmgs_sz-1));
            IIAWMse = sqrt(((wmgs_tsz2a-1).*IIAWMse.^2 + (wmgs_sz-1).*wmgs_se.^2 + wmgs_tsz2a.*wmgs_sz./(wmgs_tsz2a+wmgs_sz).*(omd.^2+IIAWM.^2-2*omd.*IIAWM))./(wmgs_tsz2a+wmgs_sz-1));
            IIAWMt1se = sqrt(((wmgs_tsz2a-1).*IIAWMt1se.^2 + (wmgs_sz-1).*wt1gs_se.^2 + wmgs_tsz2a.*wmgs_sz./(wmgs_tsz2a+wmgs_sz).*(ome.^2+IIAWMt1.^2-2*ome.*IIAWMt1))./(wmgs_tsz2a+wmgs_sz-1));
            IIAWMt2se = sqrt(((wmgs_tsz2a-1).*IIAWMt2se.^2 + (wmgs_sz-1).*wt2gs_se.^2 + wmgs_tsz2a.*wmgs_sz./(wmgs_tsz2a+wmgs_sz).*(omf.^2+IIAWMt2.^2-2*omf.*IIAWMt2))./(wmgs_tsz2a+wmgs_sz-1));
            gmgs_tsz2a = gmgs_tsz2a + gmgs_sz;
            wmgs_tsz2a = wmgs_tsz2a + wmgs_sz;
            
            oman = IIAGMn;
            ombn = IIAGMt1n;
            omcn = IIAGMt2n;
            omdn = IIAWMn;
            omen = IIAWMt1n;
            omfn = IIAWMt2n;
            IIAGMn = (IIAGMn.*gmgs_tsz2an + (gmgs_mnn).*gmgs_szn)./(gmgs_tsz2an + gmgs_szn);
            IIAGMt1n = (IIAGMt1n.*gmgs_tsz2an + (gt1gs_mnn).*gmgs_szn)./(gmgs_tsz2an + gmgs_szn);
            IIAGMt2n = (IIAGMt2n.*gmgs_tsz2an + (gt2gs_mnn).*gmgs_szn)./(gmgs_tsz2an + gmgs_szn);
            IIAWMn= (IIAWMn.*wmgs_tsz2an + (wmgs_mnn).*wmgs_szn)./(wmgs_tsz2an + wmgs_szn);
            IIAWMt1n = (IIAWMt1n.*wmgs_tsz2an + (wt1gs_mnn).*wmgs_szn)./(wmgs_tsz2an + wmgs_szn);
            IIAWMt2n = (IIAWMt2n.*wmgs_tsz2an + (wt2gs_mnn).*wmgs_szn)./(wmgs_tsz2an + wmgs_szn);
            IIAGMsen = sqrt(((gmgs_tsz2an-1).*IIAGMsen.^2 + (gmgs_szn-1).*gmgs_sen.^2 + gmgs_tsz2an.*gmgs_szn./(gmgs_tsz2an+gmgs_szn).*(oman.^2+IIAGMn.^2-2*oman.*IIAGMn))./(gmgs_tsz2an+gmgs_szn-1));
            IIAGMt1sen = sqrt(((gmgs_tsz2an-1).*IIAGMt1sen.^2 + (gmgs_szn-1).*gt1gs_sen.^2 + gmgs_tsz2an.*gmgs_szn./(gmgs_tsz2an+gmgs_szn).*(ombn.^2+IIAGMt1n.^2-2*ombn.*IIAGMt1n))./(gmgs_tsz2an+gmgs_szn-1));
            IIAGMt2sen = sqrt(((gmgs_tsz2an-1).*IIAGMt2sen.^2 + (gmgs_szn-1).*gt2gs_sen.^2 + gmgs_tsz2an.*gmgs_szn./(gmgs_tsz2an+gmgs_szn).*(omcn.^2+IIAGMt2n.^2-2*omcn.*IIAGMt2n))./(gmgs_tsz2an+gmgs_szn-1));
            IIAWMsen = sqrt(((wmgs_tsz2an-1).*IIAWMsen.^2 + (wmgs_szn-1).*wmgs_sen.^2 + wmgs_tsz2an.*wmgs_szn./(wmgs_tsz2an+wmgs_szn).*(omdn.^2+IIAWMn.^2-2*omdn.*IIAWMn))./(wmgs_tsz2an+wmgs_szn-1));
            IIAWMt1sen = sqrt(((wmgs_tsz2an-1).*IIAWMt1sen.^2 + (wmgs_szn-1).*wt1gs_sen.^2 + wmgs_tsz2an.*wmgs_szn./(wmgs_tsz2an+wmgs_szn).*(omen.^2+IIAWMt1n.^2-2*omen.*IIAWMt1n))./(wmgs_tsz2an+wmgs_szn-1));
            IIAWMt2sen = sqrt(((wmgs_tsz2an-1).*IIAWMt2sen.^2 + (wmgs_szn-1).*wt2gs_sen.^2 + wmgs_tsz2an.*wmgs_szn./(wmgs_tsz2an+wmgs_szn).*(omfn.^2+IIAWMt2n.^2-2*omfn.*IIAWMt2n))./(wmgs_tsz2an+wmgs_szn-1));
            gmgs_tsz2an = gmgs_tsz2an + gmgs_szn;
            wmgs_tsz2an = wmgs_tsz2an + wmgs_szn;
        case "mMCD"
            if ~ismember(subjID(p), temporal_exc)
                nmcd = nmcd+1;
                oma = MCDGM;
                omb = MCDGMt1;
                omc = MCDGMt2;
                omd = MCDWM;
                ome = MCDWMt1;
                omf = MCDWMt2;
                MCDGM = (MCDGM.*gmgs_tszmcd + (gmgs_mn).*gmgs_sz)./(gmgs_tszmcd + gmgs_sz);
                MCDGMt1 = (MCDGMt1.*gmgs_tszmcd + (gt1gs_mn).*gmgs_sz)./(gmgs_tszmcd + gmgs_sz);
                MCDGMt2 = (MCDGMt2.*gmgs_tszmcd + (gt2gs_mn).*gmgs_sz)./(gmgs_tszmcd + gmgs_sz);
                MCDWM= (MCDWM.*wmgs_tszmcd + (wmgs_mn).*wmgs_sz)./(wmgs_tszmcd + wmgs_sz);
                MCDWMt1 = (MCDWMt1.*wmgs_tszmcd + (wt1gs_mn).*wmgs_sz)./(wmgs_tszmcd + wmgs_sz);
                MCDWMt2 = (MCDWMt2.*wmgs_tszmcd + (wt2gs_mn).*wmgs_sz)./(wmgs_tszmcd + wmgs_sz);
                MCDGMse = sqrt(((gmgs_tszmcd-1).*MCDGMse.^2 + (gmgs_sz-1).*gmgs_se.^2 + gmgs_tszmcd.*gmgs_sz./(gmgs_tszmcd+gmgs_sz).*(oma.^2+MCDGM.^2-2*oma.*MCDGM))./(gmgs_tszmcd+gmgs_sz-1));
                MCDGMt1se = sqrt(((gmgs_tszmcd-1).*MCDGMt1se.^2 + (gmgs_sz-1).*gt1gs_se.^2 + gmgs_tszmcd.*gmgs_sz./(gmgs_tszmcd+gmgs_sz).*(omb.^2+MCDGMt1.^2-2*omb.*MCDGMt1))./(gmgs_tszmcd+gmgs_sz-1));
                MCDGMt2se = sqrt(((gmgs_tszmcd-1).*MCDGMt2se.^2 + (gmgs_sz-1).*gt2gs_se.^2 + gmgs_tszmcd.*gmgs_sz./(gmgs_tszmcd+gmgs_sz).*(omc.^2+MCDGMt2.^2-2*omc.*MCDGMt2))./(gmgs_tszmcd+gmgs_sz-1));
                MCDWMse = sqrt(((wmgs_tszmcd-1).*MCDWMse.^2 + (wmgs_sz-1).*wmgs_se.^2 + wmgs_tszmcd.*wmgs_sz./(wmgs_tszmcd+wmgs_sz).*(omd.^2+MCDWM.^2-2*omd.*MCDWM))./(wmgs_tszmcd+wmgs_sz-1));
                MCDWMt1se = sqrt(((wmgs_tszmcd-1).*MCDWMt1se.^2 + (wmgs_sz-1).*wt1gs_se.^2 + wmgs_tszmcd.*wmgs_sz./(wmgs_tszmcd+wmgs_sz).*(ome.^2+MCDWMt1.^2-2*ome.*MCDWMt1))./(wmgs_tszmcd+wmgs_sz-1));
                MCDWMt2se = sqrt(((wmgs_tszmcd-1).*MCDWMt2se.^2 + (wmgs_sz-1).*wt2gs_se.^2 + wmgs_tszmcd.*wmgs_sz./(wmgs_tszmcd+wmgs_sz).*(omf.^2+MCDWMt2.^2-2*omf.*MCDWMt2))./(wmgs_tszmcd+wmgs_sz-1));
                gmgs_tszmcd = gmgs_tszmcd + gmgs_sz;
                wmgs_tszmcd = wmgs_tszmcd + wmgs_sz;
    
                oman = MCDGMn;
                ombn = MCDGMt1n;
                omcn = MCDGMt2n;
                omdn = MCDWMn;
                omen = MCDWMt1n;
                omfn = MCDWMt2n;
                MCDGMn = (MCDGMn.*gmgs_tszmcdn + (gmgs_mnn).*gmgs_szn)./(gmgs_tszmcdn + gmgs_szn);
                MCDGMt1n = (MCDGMt1n.*gmgs_tszmcdn + (gt1gs_mnn).*gmgs_szn)./(gmgs_tszmcdn + gmgs_szn);
                MCDGMt2n = (MCDGMt2n.*gmgs_tszmcdn + (gt2gs_mnn).*gmgs_szn)./(gmgs_tszmcdn + gmgs_szn);
                MCDWMn= (MCDWMn.*wmgs_tszmcdn + (wmgs_mnn).*wmgs_szn)./(wmgs_tszmcdn + wmgs_szn);
                MCDWMt1n = (MCDWMt1n.*wmgs_tszmcdn + (wt1gs_mnn).*wmgs_szn)./(wmgs_tszmcdn + wmgs_szn);
                MCDWMt2n = (MCDWMt2n.*wmgs_tszmcdn + (wt2gs_mnn).*wmgs_szn)./(wmgs_tszmcdn + wmgs_szn);
                MCDGMsen = sqrt(((gmgs_tszmcdn-1).*MCDGMsen.^2 + (gmgs_szn-1).*gmgs_sen.^2 + gmgs_tszmcdn.*gmgs_szn./(gmgs_tszmcdn+gmgs_szn).*(oman.^2+MCDGMn.^2-2*oman.*MCDGMn))./(gmgs_tszmcdn+gmgs_szn-1));
                MCDGMt1sen = sqrt(((gmgs_tszmcdn-1).*MCDGMt1sen.^2 + (gmgs_szn-1).*gt1gs_sen.^2 + gmgs_tszmcdn.*gmgs_szn./(gmgs_tszmcdn+gmgs_szn).*(ombn.^2+MCDGMt1n.^2-2*ombn.*MCDGMt1n))./(gmgs_tszmcdn+gmgs_szn-1));
                MCDGMt2sen = sqrt(((gmgs_tszmcdn-1).*MCDGMt2sen.^2 + (gmgs_szn-1).*gt2gs_sen.^2 + gmgs_tszmcdn.*gmgs_szn./(gmgs_tszmcdn+gmgs_szn).*(omcn.^2+MCDGMt2n.^2-2*omcn.*MCDGMt2n))./(gmgs_tszmcdn+gmgs_szn-1));
                MCDWMsen = sqrt(((wmgs_tszmcdn-1).*MCDWMsen.^2 + (wmgs_szn-1).*wmgs_sen.^2 + wmgs_tszmcdn.*wmgs_szn./(wmgs_tszmcdn+wmgs_szn).*(omdn.^2+MCDWMn.^2-2*omdn.*MCDWMn))./(wmgs_tszmcdn+wmgs_szn-1));
                MCDWMt1sen = sqrt(((wmgs_tszmcdn-1).*MCDWMt1sen.^2 + (wmgs_szn-1).*wt1gs_sen.^2 + wmgs_tszmcdn.*wmgs_szn./(wmgs_tszmcdn+wmgs_szn).*(omen.^2+MCDWMt1n.^2-2*omen.*MCDWMt1n))./(wmgs_tszmcdn+wmgs_szn-1));
                MCDWMt2sen = sqrt(((wmgs_tszmcdn-1).*MCDWMt2sen.^2 + (wmgs_szn-1).*wt2gs_sen.^2 + wmgs_tszmcdn.*wmgs_szn./(wmgs_tszmcdn+wmgs_szn).*(omfn.^2+MCDWMt2n.^2-2*omfn.*MCDWMt2n))./(wmgs_tszmcdn+wmgs_szn-1));
                gmgs_tszmcdn = gmgs_tszmcdn + gmgs_szn;
                wmgs_tszmcdn = wmgs_tszmcdn + wmgs_szn;
            end
end
end

% final plot for group response
IIBGMse = IIBGMse./sqrt(gmgs_tsz2b);
IIBGMt1se = IIBGMt1se./sqrt(gmgs_tsz2b);
IIBGMt2se = IIBGMt2se./sqrt(gmgs_tsz2b);
IIBWMse = IIBWMse./sqrt(gmgs_tsz2b);
IIBWMt1se = IIBWMt1se./sqrt(gmgs_tsz2b);
IIBWMt2se = IIBWMt2se./sqrt(gmgs_tsz2b);
IIAGMse = IIAGMse./sqrt(gmgs_tsz2a);
IIAGMt1se = IIAGMt1se./sqrt(gmgs_tsz2a);
IIAGMt2se = IIAGMt2se./sqrt(gmgs_tsz2a);
IIAWMse = IIAWMse./sqrt(gmgs_tsz2a);
IIAWMt1se = IIAWMt1se./sqrt(gmgs_tsz2a);
IIAWMt2se = IIAWMt2se./sqrt(gmgs_tsz2a);
MCDGMse = MCDGMse./sqrt(gmgs_tszmcd);
MCDGMt1se = MCDGMt1se./sqrt(gmgs_tszmcd);
MCDGMt2se = MCDGMt2se./sqrt(gmgs_tszmcd);
MCDWMse = MCDWMse./sqrt(gmgs_tszmcd);
MCDWMt1se = MCDWMt1se./sqrt(gmgs_tszmcd);
MCDWMt2se = MCDWMt2se./sqrt(gmgs_tszmcd);

IIBGMsen = IIBGMsen./sqrt(gmgs_tsz2bn);
IIBGMt1sen = IIBGMt1sen./sqrt(gmgs_tsz2bn);
IIBGMt2sen = IIBGMt2sen./sqrt(gmgs_tsz2bn);
IIBWMsen = IIBWMsen./sqrt(gmgs_tsz2bn);
IIBWMt1sen = IIBWMt1sen./sqrt(gmgs_tsz2bn);
IIBWMt2sen = IIBWMt2sen./sqrt(gmgs_tsz2bn);
IIAGMsen = IIAGMsen./sqrt(gmgs_tsz2an);
IIAGMt1sen = IIAGMt1sen./sqrt(gmgs_tsz2an);
IIAGMt2sen = IIAGMt2sen./sqrt(gmgs_tsz2an);
IIAWMsen = IIAWMsen./sqrt(gmgs_tsz2an);
IIAWMt1sen = IIAWMt1sen./sqrt(gmgs_tsz2an);
IIAWMt2sen = IIAWMt2sen./sqrt(gmgs_tsz2an);
MCDGMsen = MCDGMsen./sqrt(gmgs_tszmcdn);
MCDGMt1sen = MCDGMt1sen./sqrt(gmgs_tszmcdn);
MCDGMt2sen = MCDGMt2sen./sqrt(gmgs_tszmcdn);
MCDWMsen = MCDWMsen./sqrt(gmgs_tszmcdn);
MCDWMt1sen = MCDWMt1sen./sqrt(gmgs_tszmcdn);
MCDWMt2sen = MCDWMt2sen./sqrt(gmgs_tszmcdn);

%% shaded error bar plot
% IIB group
figure()
subplot(2,3,1)
% errorbar(dms(2:5),IIBGM(1:4),IIBGMse(1:4),"Color", "#A30234")
% errorbar(dms(5:end),IIBGM(4:end),IIBGMse(4:end), "Color", "#ce8080")
% errorbar(dms(2:5),IIBGMn(1:4),IIBGMsen(1:4),"Color", "#002157")
% errorbar(dms(5:end),IIBGMn(4:end),IIBGMsen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIBGM(1:4), IIBGMse(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBGM(4:end), IIBGMse(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIBGMn(1:4), IIBGMsen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBGMn(4:end), IIBGMsen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('GM T1w/T2w', 'FontSize', 20)
subtitle('Type IIB', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
% yticks(linspace(0.2,0.9,15))
ylim([0.2 0.6])
ax=gca;
ax.FontSize = 15;

subplot(2,3,4)
% errorbar(dms(2:5),IIBWM(1:4),IIBWMse(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIBWM(4:end),IIBWMse(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIBWMn(1:4),IIBWMsen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIBWMn(4:end),IIBWMsen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIBWM(1:4), IIBWMse(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBWM(4:end), IIBWMse(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIBWMn(1:4), IIBWMsen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBWMn(4:end), IIBWMsen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('WM T1w/T2w', 'FontSize', 20)
subtitle('Type IIB', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
ylim([0.5 1])
ax=gca;
ax.FontSize = 15;

subplot(2,3,2)
% errorbar(dms(2:5),IIBGMt1(1:4),IIBGMt1se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIBGMt1(4:end),IIBGMt1se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIBGMt1n(1:4),IIBGMt1sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIBGMt1n(4:end),IIBGMt1sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIBGMt1(1:4), IIBGMt1se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBGMt1(4:end), IIBGMt1se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIBGMt1n(1:4), IIBGMt1sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBGMt1n(4:end), IIBGMt1sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('GM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
ylim([1000 1300])
ax=gca;
ax.FontSize = 15;

subplot(2,3,5)
% errorbar(dms(2:5),IIBWMt1(1:4),IIBWMt1se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIBWMt1(4:end),IIBWMt1se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIBWMt1n(1:4),IIBWMt1sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIBWMt1n(4:end),IIBWMt1sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIBWMt1(1:4), IIBWMt1se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBWMt1(4:end), IIBWMt1se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIBWMt1n(1:4), IIBWMt1sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBWMt1n(4:end), IIBWMt1sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('WM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
ylim([750 1200])
ax=gca;
ax.FontSize = 15;

subplot(2,3,3)
% errorbar(dms(2:5),IIBGMt2(1:4),IIBGMt2se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIBGMt2(4:end),IIBGMt2se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIBGMt2n(1:4),IIBGMt2sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIBGMt2n(4:end),IIBGMt2sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIBGMt2(1:4), IIBGMt2se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBGMt2(4:end), IIBGMt2se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIBGMt2n(1:4), IIBGMt2sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBGMt2n(4:end), IIBGMt2sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('GM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
ylim([40 60])
ax=gca;
ax.FontSize = 15;

subplot(2,3,6)
% errorbar(dms(2:5),IIBWMt2(1:4),IIBWMt2se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIBWMt2(4:end),IIBWMt2se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIBWMt2n(1:4),IIBWMt2sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIBWMt2n(4:end),IIBWMt2sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIBWMt2(1:4), IIBWMt2se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBWMt2(4:end), IIBWMt2se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIBWMt2n(1:4), IIBWMt2sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIBWMt2n(4:end), IIBWMt2sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('WM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
ylim([40 50])
ax=gca;
ax.FontSize = 15;

% IIA group
figure()
subplot(2,3,1)
% errorbar(dms(2:5),IIAGM(1:4),IIAGMse(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIAGM(4:end),IIAGMse(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIAGMn(1:4),IIAGMsen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIAGMn(4:end),IIAGMsen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIAGM(1:4), IIAGMse(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAGM(4:end), IIAGMse(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIAGMn(1:4), IIAGMsen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAGMn(4:end), IIAGMsen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('GM T1/T2', 'FontSize', 20)
subtitle('Type IIA', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
ylim([0.2 0.6])
ax=gca;
ax.FontSize = 15;

subplot(2,3,4)
% errorbar(dms(2:5),IIAWM(1:4),IIAWMse(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIAWM(4:end),IIAWMse(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIAWMn(1:4),IIAWMsen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIAWMn(4:end),IIAWMsen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIAWM(1:4), IIAWMse(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAWM(4:end), IIAWMse(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIAWMn(1:4), IIAWMsen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAWMn(4:end), IIAWMsen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('WM T1w/T2w', 'FontSize', 20)
subtitle('Type IIA', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
ylim([1000 1300])
ax=gca;
ax.FontSize = 15;

subplot(2,3,2)
% errorbar(dms(2:5),IIAGMt1(1:4),IIAGMt1se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIAGMt1(4:end),IIAGMt1se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIAGMt1n(1:4),IIAGMt1sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIAGMt1n(4:end),IIAGMt1sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIAGMt1(1:4), IIAGMt1se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAGMt1(4:end), IIAGMt1se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIAGMt1n(1:4), IIAGMt1sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAGMt1n(4:end), IIAGMt1sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('GM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,5)
% errorbar(dms(2:5),IIAWMt1(1:4),IIAWMt1se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIAWMt1(4:end),IIAWMt1se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIAWMt1n(1:4),IIAWMt1sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIAWMt1n(4:end),IIAWMt1sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIAWMt1(1:4), IIAWMt1se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAWMt1(4:end), IIAWMt1se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIAWMt1n(1:4), IIAWMt1sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAWMt1n(4:end), IIAWMt1sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('WM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,3)
% errorbar(dms(2:5),IIAGMt2(1:4),IIAGMt2se(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),IIAGMt2(4:end),IIAGMt2se(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),IIAGMt2n(1:4),IIAGMt2sen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),IIAGMt2n(4:end),IIAGMt2sen(4:end), "Color", "#0076c0")
shadedErrorBar(dms(2:5),IIAGMt2(1:4), IIAGMt2se(1:4),'lineprops',{'-','Color', '#A30234'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAGMt2(4:end), IIAGMt2se(4:end),'lineprops',{'-','Color', '#ce8080'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(2:5),IIAGMt2n(1:4), IIAGMt2sen(1:4),'lineprops',{'-','Color', '#002157'},'patchSaturation',0.33)
hold on
shadedErrorBar(dms(5:end),IIAGMt2n(4:end), IIAGMt2sen(4:end),'lineprops',{'-','Color', '#0076c0'},'patchSaturation',0.33)
title('GM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,6)
errorbar(dms(2:5),IIAGMt2(1:4),IIAGMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAGMt2(4:end),IIAGMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAGMt2n(1:4),IIAGMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAGMt2n(4:end),IIAGMt2sen(4:end), "Color", "#0076c0")
title('WM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

% MCD group
figure()
subplot(2,3,1)
% errorbar(dms(2:5),MCDGM(1:4),MCDGMse(1:4),"Color", "#A30234")
% hold on
% errorbar(dms(5:end),MCDGM(4:end),MCDGMse(4:end), "Color", "#ce8080")
% hold on
% errorbar(dms(2:5),MCDGMn(1:4),MCDGMsen(1:4),"Color", "#002157")
% hold on
% errorbar(dms(5:end),MCDGMn(4:end),MCDGMsen(4:end), "Color", "#0076c0")
title('GM T1/T2', 'FontSize', 20)
subtitle('mMCD', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,4)
errorbar(dms(2:5),MCDWM(1:4),MCDWMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDWM(4:end),MCDWMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDWMn(1:4),MCDWMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDWMn(4:end),MCDWMsen(4:end), "Color", "#0076c0")
title('WM T1/T2', 'FontSize', 20)
subtitle('mMCD', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,2)
errorbar(dms(2:5),MCDGMt1(1:4),MCDGMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDGMt1(4:end),MCDGMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDGMt1n(1:4),MCDGMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDGMt1n(4:end),MCDGMt1sen(4:end), "Color", "#0076c0")
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,5)
errorbar(dms(2:5),MCDWMt1(1:4),MCDWMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDWMt1(4:end),MCDWMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDWMt1n(1:4),MCDWMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDWMt1n(4:end),MCDWMt1sen(4:end), "Color", "#0076c0")
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,3)
errorbar(dms(2:5),MCDGMt2(1:4),MCDGMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDGMt2(4:end),MCDGMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDGMt2n(1:4),MCDGMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDGMt2n(4:end),MCDGMt2sen(4:end), "Color", "#0076c0")
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,6)
errorbar(dms(2:5),MCDWMt2(1:4),MCDWMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDWMt2(4:end),MCDWMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDWMt2n(1:4),MCDWMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDWMt2n(4:end),MCDWMt2sen(4:end), "Color", "#0076c0")
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

%% individually normalize
figure()
tcl = tiledlayout(2,3);

ax1 = nexttile(tcl)
plot(dms(2:end),IIBGM(1:end)-(IIBGMn.*gmgs_tsz2bn)./(gmgs_tsz2bn)+1,'r-o')
hold on
plot(dms(2:end),IIAGM(1:end)-(IIAGMn.*gmgs_tsz2an)./(gmgs_tsz2an)+1,'b-*')
hold on
plot(dms(2:end),MCDGM(1:end)-(MCDGMn.*gmgs_tszmcdn)./(gmgs_tszmcdn)+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
% hold on
% plot(dms, 1-f1.a*exp(f1.b*dms))
% hold on
% plot(dms, 1-f2.a*exp(f2.b*dms))
% hold on
% plot(dms, 1-f3.a*exp(f3.b*dms))
title('GM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
% yticks(ax1, linspace(0.9,1,11))
ylim(ax1, [0.8, 1.05])
xticks(ax1, linspace(60,220,9))
xlim([50 210])
ax=gca;
ax.FontSize = 15;


ax2 = nexttile(tcl)
plot(dms(2:end),IIBGMt1(1:end)-(IIBGMt1n)+1,'r-o')
hold on
plot(dms(2:end),IIAGMt1(1:end)-(IIAGMt1n)+1,'b-*')
hold on
plot(dms(2:end),MCDGMt1(1:end)-(MCDGMt1n)+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))) ,'g-')
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(ax2, linspace(0.9,1,11))
xticks(ax2, linspace(60,220,9))
xlim([50 210])
ax=gca;
ax.FontSize = 15;

ax3=nexttile(tcl)
plot(dms(2:end),IIBGMt2(1:end)-(IIBGMt2n)+1,'r-o')
hold on
plot(dms(2:end),IIAGMt2(1:end)-(IIAGMt2n)+1,'b-*')
hold on
plot(dms(2:end),MCDGMt2(1:end)-(MCDGMt2n)+1,'k-+')
hold on 
plot(dms(2:end),zeros(size(dms(2:end))),'g-')
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(ax3, linspace(0.9,1,11))
yticks(linspace(-2,6,5))
ylim([-2 6])
xticks(ax3, linspace(60,220,9))
xlim([50 210])
ax=gca;
ax.FontSize = 15;

ax4=nexttile(tcl)
plot(dms(2:end),IIBWM(1:end)-(IIBWMn.*wmgs_tsz2bn)./(wmgs_tsz2bn)+1,'r-o')
hold on
plot(dms(2:end),IIAWM(1:end)-(IIAWMn.*wmgs_tsz2an)./(wmgs_tsz2an)+1,'b-*')
hold on
plot(dms(2:end),MCDWM(1:end)-(MCDWMn.*wmgs_tszmcdn)./(wmgs_tszmcdn)+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax4, linspace(0.5,1,11))
xticks(ax4, linspace(60,220,9))
xlim([50 210])
ax=gca;
ax.FontSize = 15;

ax5=nexttile(tcl)
plot(dms(2:end),IIBWMt1(1:end)-(IIBWMt1n)+1,'r-o')
hold on
plot(dms(2:end),IIAWMt1(1:end)-(IIAWMt1n)+1,'b-*')
hold on
plot(dms(2:end),MCDWMt1(1:end)-(MCDWMt1n)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
xticks(ax5, linspace(60,220,9))
xlim([50 210])
ax=gca;
ax.FontSize = 15;


ax6=nexttile(tcl)
l1 = plot(dms(2:end),IIBWMt2(1:end)-(IIBWMt2n)+1,'r-o','DisplayName','IIB')
hold on
l3 = plot(dms(2:end),IIAWMt2(1:end)-(IIAWMt2n)+1,'b-*','DisplayName','IIA')
hold on
l5 = plot(dms(2:end),MCDWMt2(1:end)-(MCDWMt2n)+1,'k-+','DisplayName','mMCD')
hold on 
l7 = plot(dms(2:end), zeros(size(dms(2:end))), 'g-', 'DisplayName','Healthy Control')
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(-2,20,12))
ylim([-2 18])
xticks(ax6, linspace(60,220,9))
xlim([50 210])
ax=gca;
ax.FontSize = 15;

hL = legend([l1,l3,l5,l7]); 
hL.Layout.Tile = 'East';
%% individually normalize without mMCD
figure()
tcl = tiledlayout(2,3);

ax1 = nexttile(tcl)
% plot(dms(2:end),IIBGM(1:end)-(IIBGMn.*gmgs_tsz2bn)./(gmgs_tsz2bn)+1,'r-o')
% hold on
% plot(dms(2:end),IIAGM(1:end)-(IIAGMn.*gmgs_tsz2an)./(gmgs_tsz2an)+1,'b-*')
% hold on
% plot(dms(2:end),ones(size(dms(2:end))),'g-')
errorbar(dms(2:end),IIBGM(1:end)-IIBGMn+1, IIBGMse, 'r-')
hold on
errorbar(dms(2:end),IIAGM(1:end)-IIAGMn+1, IIAGMse, 'b-')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
% hold on
% plot(dms, 1-f1.a*exp(f1.b*dms))
% hold on
% plot(dms, 1-f2.a*exp(f2.b*dms))
% hold on
% plot(dms, 1-f3.a*exp(f3.b*dms))
title('GM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
ylim(ax1, [0.8, 1.05])
ax=gca;
ax.FontSize = 15;

ax2 = nexttile(tcl)
% plot(dms(2:end),IIBGMt1(1:end)-(IIBGMt1n),'r-o')
% hold on
% plot(dms(2:end),IIAGMt1(1:end)-(IIAGMt1n),'b-*')
% hold on
% plot(dms(2:end),zeros(size(dms(2:end))) ,'g-')
errorbar(dms(2:end),IIBGMt1(1:end)-IIBGMt1n, IIBGMt1se, 'r-')
hold on
errorbar(dms(2:end),IIAGMt1(1:end)-IIAGMt1n, IIAGMt1se, 'b-')
hold on
plot(dms(2:end),zeros(size(dms(2:end))),'g-')
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
ax=gca;
ax.FontSize = 15;

ax3=nexttile(tcl)
% plot(dms(2:end),IIBGMt2(1:end)-(IIBGMt2n),'r-o')
% hold on
% plot(dms(2:end),IIAGMt2(1:end)-(IIAGMt2n),'b-*')
% hold on
% plot(dms(2:end),zeros(size(dms(2:end))),'g-')
errorbar(dms(2:end),IIBGMt2(1:end)-IIBGMt2n, IIBGMt2se, 'r-')
hold on
errorbar(dms(2:end),IIAGMt2(1:end)-IIAGMt2n, IIAGMt2se, 'b-')
hold on
plot(dms(2:end),zeros(size(dms(2:end))),'g-')
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
ylim(ax3, [-1, 15])
% yticks(ax3, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax4=nexttile(tcl)
% plot(dms(2:end),IIBWM(1:end)-(IIBWMn.*wmgs_tsz2bn)./(wmgs_tsz2bn)+1,'r-o')
% hold on
% plot(dms(2:end),IIAWM(1:end)-(IIAWMn.*wmgs_tsz2an)./(wmgs_tsz2an)+1,'b-*')
% hold on
% plot(dms(2:end),ones(size(dms(2:end))),'g-')
errorbar(dms(2:end),IIBWM(1:end)-IIBWMn+1, IIBWMse, 'r-')
hold on
errorbar(dms(2:end),IIAWM(1:end)-IIAWMn+1, IIAWMse, 'b-')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
ylim(ax4, [0.8, 1.05])
ax=gca;
ax.FontSize = 15;

ax5=nexttile(tcl)
% plot(dms(2:end),IIBWMt1(1:end)-(IIBWMt1n),'r-o')
% hold on
% plot(dms(2:end),IIAWMt1(1:end)-(IIAWMt1n),'b-*')
% hold on
% plot(dms(2:end),zeros(size(dms(2:end))),'g-')
errorbar(dms(2:end),IIBWMt1(1:end)-IIBWMt1n, IIBWMt1se, 'r-')
hold on
errorbar(dms(2:end),IIAWMt1(1:end)-IIAWMt1n, IIAWMt1se, 'b-')
hold on
plot(dms(2:end),zeros(size(dms(2:end))),'g-')
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;


ax6=nexttile(tcl)
l1 = errorbar(dms(2:end),IIBWMt2(1:end)-(IIBWMt2n), IIBWMt2se, 'r-','DisplayName','IIB ROI')
hold on
l3 = errorbar(dms(2:end),IIAWMt2(1:end)-(IIAWMt2n), IIAWMt2se, 'b-','DisplayName','IIA ROI')
hold on
l7 = plot(dms(2:end), zeros(size(dms(2:end))), 'g-', 'DisplayName','Homotopic Control ROI')
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
ylim(ax6, [-1, 15])
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

hL = legend([l1,l3,l7]); 
hL.Layout.Tile = 'East';
%% plot as one

figure()
tcl = tiledlayout(2,3);

nexttile(tcl)
plot(dms(2:4),IIBGM(1:3),'r-o')
hold on
plot(dms(4:end),IIBGM(3:end),'r-o')
hold on
plot(dms(2:4),IIAGM(1:3),'b-*')
hold on
plot(dms(4:end),IIAGM(3:end),'b-*')
hold on
plot(dms(2:4),MCDGM(1:3),'k-+')
hold on
plot(dms(4:end),MCDGM(3:end),'k-+')
hold on 
plot(dms(2:end),(IIBGMn + IIAGMn + MCDGMn)/3,'g-')
% hold on
% plot(dms, 1-f1.a*exp(f1.b*dms))
% hold on
% plot(dms, 1-f2.a*exp(f2.b*dms))
% hold on
% plot(dms, 1-f3.a*exp(f3.b*dms))
title('GM T1w/T2w', 'FontSize', 20)
subtitle('Type IIB', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

nexttile(tcl)
plot(dms(2:4),IIBGMt1(1:3),'r-o')
hold on
plot(dms(4:end),IIBGMt1(3:end),'r-o')
hold on
plot(dms(2:4),IIAGMt1(1:3),'b-*')
hold on
plot(dms(4:end),IIAGMt1(3:end),'b-*')
hold on
plot(dms(2:4),MCDGMt1(1:3),'k-+')
hold on
plot(dms(4:end),MCDGMt1(3:end),'k-+')
hold on 
plot(dms(2:end),(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3,'g-')
title('GM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

nexttile(tcl)
plot(dms(2:4),IIBGMt2(1:3),'r-o')
hold on
plot(dms(4:end),IIBGMt2(3:end),'r-o')
hold on
plot(dms(2:4),IIAGMt2(1:3),'b-*')
hold on
plot(dms(4:end),IIAGMt2(3:end),'b-*')
hold on
plot(dms(2:4),MCDGMt2(1:3),'k-+')
hold on
plot(dms(4:end),MCDGMt2(3:end),'k-+')
hold on 
plot(dms(2:end),(IIBGMt2n + IIAGMt2n + MCDGMt2n)/3,'g-')
title('GM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

nexttile(tcl)
plot(dms(2:4),IIBWM(1:3),'r-o')
hold on
plot(dms(4:end),IIBWM(3:end),'r-o')
hold on
plot(dms(2:4),IIAWM(1:3),'b-*')
hold on
plot(dms(4:end),IIAWM(3:end),'b-*')
hold on
plot(dms(2:4),MCDWM(1:3),'k-+')
hold on
plot(dms(4:end),MCDWM(3:end),'k-+')
hold on 
plot(dms(2:end),(IIBWMn + IIAWMn + MCDWMn)/3,'g-')
title('WM T1w/T2w', 'FontSize', 20)
subtitle('Type IIB', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

nexttile(tcl)
plot(dms(2:4),IIBWMt1(1:3),'r-o')
hold on
plot(dms(4:end),IIBWMt1(3:end),'r-o')
hold on
plot(dms(2:4),IIAWMt1(1:3),'b-*')
hold on
plot(dms(4:end),IIAWMt1(3:end),'b-*')
hold on
plot(dms(2:4),MCDWMt1(1:3),'k-+')
hold on
plot(dms(4:end),MCDWMt1(3:end),'k-+')
hold on 
plot(dms(2:end),(IIBWMt1n + IIAWMt1n + MCDWMt1n)/3,'g-')
title('WM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;


nexttile(tcl)
l1 = plot(dms(2:4),IIBWMt2(1:3),'r-o','DisplayName','IIB ROI')
hold on
l2 = plot(dms(4:end),IIBWMt2(3:end),'r-o','DisplayName','IIB Healthy')
hold on
l3 = plot(dms(2:4),IIAWMt2(1:3),'b-*','DisplayName','IIA ROI')
hold on
l4 = plot(dms(4:end),IIAWMt2(3:end),'b-*','DisplayName','IIA Healthy')
hold on
l5 = plot(dms(2:4),MCDWMt2(1:3),'k-+','DisplayName','mMCD ROI')
hold on
l6 = plot(dms(4:end),MCDWMt2(3:end),'k-+','DisplayName','mMCD Healthy')
hold on 
l7 = plot(dms(2:end),(IIBWMt2n + IIAWMt2n + MCDWMt2n)/3, 'g-', 'DisplayName','Healthy Control')
title('WM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

hL = legend([l1,l2,l3,l4,l5,l6,l7]); 
hL.Layout.Tile = 'East';
%% plot as one after normalization

% n2a = ones(size(dms(2:end)))*n2a;
% n2b = ones(size(dms(2:end)))*n2b;
% nmcd = ones(size(dms(2:end)))*nmcd;

figure()
tcl = tiledlayout(2,3);

ax1 = nexttile(tcl)
plot(dms(2:end),IIBGM(1:end)-(IIBGMn.*gmgs_tsz2bn + IIAGMn.*gmgs_tsz2an + MCDGMn.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'r-o')
hold on
plot(dms(2:end),IIAGM(1:end)-(IIBGMn.*gmgs_tsz2bn + IIAGMn.*gmgs_tsz2an + MCDGMn.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'b-*')
hold on
plot(dms(2:end),MCDGM(1:end)-(IIBGMn.*gmgs_tsz2bn + IIAGMn.*gmgs_tsz2an + MCDGMn.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
% hold on
% plot(dms, 1-f1.a*exp(f1.b*dms))
% hold on
% plot(dms, 1-f2.a*exp(f2.b*dms))
% hold on
% plot(dms, 1-f3.a*exp(f3.b*dms))
title('GM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax1, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;


ax2 = nexttile(tcl)
plot(dms(2:end),IIBGMt1(1:end)-(IIBGMt1n.*gmgs_tsz2bn + IIAGMt1n.*gmgs_tsz2an + MCDGMt1n.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'r-o')
hold on
plot(dms(2:end),IIAGMt1(1:end)-(IIBGMt1n.*gmgs_tsz2bn + IIAGMt1n.*gmgs_tsz2an + MCDGMt1n.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'b-*')
hold on
plot(dms(2:end),MCDGMt1(1:end)-(IIBGMt1n.*gmgs_tsz2bn + IIAGMt1n.*gmgs_tsz2an + MCDGMt1n.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'k-+')
hold on 
plot(dms(2:end),(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3-(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3+1,'g-')
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(ax2, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax3=nexttile(tcl)
plot(dms(2:end),IIBGMt2(1:end)-(IIBGMt2n.*gmgs_tsz2bn + IIAGMt2n.*gmgs_tsz2an + MCDGMt2n.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'r-o')
hold on
plot(dms(2:end),IIAGMt2(1:end)-(IIBGMt2n.*gmgs_tsz2bn + IIAGMt2n.*gmgs_tsz2an + MCDGMt2n.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'b-*')
hold on
plot(dms(2:end),MCDGMt2(1:end)-(IIBGMt2n.*gmgs_tsz2bn + IIAGMt2n.*gmgs_tsz2an + MCDGMt2n.*gmgs_tszmcdn)./(gmgs_tsz2bn+gmgs_tsz2an+gmgs_tszmcdn)+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(ax3, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax4=nexttile(tcl)
plot(dms(2:end),IIBWM(1:end)-(IIBWMn.*wmgs_tsz2bn + IIAWMn.*wmgs_tsz2an + MCDWMn.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'r-o')
hold on
plot(dms(2:end),IIAWM(1:end)-(IIBWMn.*wmgs_tsz2bn + IIAWMn.*wmgs_tsz2an + MCDWMn.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'b-*')
hold on
plot(dms(2:end),MCDWM(1:end)-(IIBWMn.*wmgs_tsz2bn + IIAWMn.*wmgs_tsz2an + MCDWMn.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax4, linspace(0.5,1,11))
ax=gca;
ax.FontSize = 15;

ax5=nexttile(tcl)
plot(dms(2:end),IIBWMt1(1:end)-(IIBWMt1n.*wmgs_tsz2bn + IIAWMt1n.*wmgs_tsz2an + MCDWMt1n.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'r-o')
hold on
plot(dms(2:end),IIAWMt1(1:end)-(IIBWMt1n.*wmgs_tsz2bn + IIAWMt1n.*wmgs_tsz2an + MCDWMt1n.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'b-*')
hold on
plot(dms(2:end),MCDWMt1(1:end)-(IIBWMt1n.*wmgs_tsz2bn + IIAWMt1n.*wmgs_tsz2an + MCDWMt1n.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;


ax6=nexttile(tcl)
l1 = plot(dms(2:end),IIBWMt2(1:end)-(IIBWMt2n.*wmgs_tsz2bn + IIAWMt2n.*wmgs_tsz2an + MCDWMt2n.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'r-o','DisplayName','IIB ROI')
hold on
l3 = plot(dms(2:end),IIAWMt2(1:end)-(IIBWMt2n.*wmgs_tsz2bn + IIAWMt2n.*wmgs_tsz2an + MCDWMt2n.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'b-*','DisplayName','IIA ROI')
hold on
l5 = plot(dms(2:end),MCDWMt2(1:end)-(IIBWMt2n.*wmgs_tsz2bn + IIAWMt2n.*wmgs_tsz2an + MCDWMt2n.*wmgs_tszmcdn)./(wmgs_tsz2bn+wmgs_tsz2an+wmgs_tszmcdn)+1,'k-+','DisplayName','mMCD ROI')
hold on 
l7 = plot(dms(2:end), ones(size(dms(2:end))), 'g-', 'DisplayName','Healthy Control')
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

hL = legend([l1,l3,l5,l7]); 
hL.Layout.Tile = 'East';

%% plot as one after normalization by a weighted average

% n2a = ones(size(dms(2:end)))*n2a;
% n2b = ones(size(dms(2:end)))*n2b;
% nmcd = ones(size(dms(2:end)))*nmcd;

figure()
tcl = tiledlayout(2,3);

ax1 = nexttile(tcl)
plot(dms(2:end),IIBGM(1:end)-(IIBGMn*7 + IIAGMn*3)./10+1,'r-o')
hold on
plot(dms(2:end),IIAGM(1:end)-(IIBGMn*7 + IIAGMn*3)./10+1,'b-*')
hold on
plot(dms(2:end),MCDGM(1:end)-(IIBGMn*7 + IIAGMn*3)./10+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
% hold on
% plot(dms, 1-f1.a*exp(f1.b*dms))
% hold on
% plot(dms, 1-f2.a*exp(f2.b*dms))
% hold on
% plot(dms, 1-f3.a*exp(f3.b*dms))
title('GM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax1, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;


ax2 = nexttile(tcl)
plot(dms(2:end),IIBGMt1(1:end)-(IIBGMt1n*7 + IIAGMt1n*3)./10+1,'r-o')
hold on
plot(dms(2:end),IIAGMt1(1:end)-(IIBGMt1n*7 + IIAGMt1n*3)./10+1,'b-*')
hold on
plot(dms(2:end),MCDGMt1(1:end)-(IIBGMt1n*7 + IIAGMt1n*3)./10+1,'k-+')
hold on 
plot(dms(2:end),(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3-(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3+1,'g-')
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(ax2, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax3=nexttile(tcl)
plot(dms(2:end),IIBGMt2(1:end)-(IIBGMt2n*7 + IIAGMt2n*3)./10+1,'r-o')
hold on
plot(dms(2:end),IIAGMt2(1:end)-(IIBGMt2n*7 + IIAGMt2n*3)./10+1,'b-*')
hold on
plot(dms(2:end),MCDGMt2(1:end)-(IIBGMt2n*7 + IIAGMt2n*3)./10+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(ax3, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax4=nexttile(tcl)
plot(dms(2:end),IIBWM(1:end)-(IIBWMn*7 + IIAWMn*3)./(10)+1,'r-o')
hold on
plot(dms(2:end),IIAWM(1:end)-(IIBWMn*7 + IIAWMn*3)./(10)+1,'b-*')
hold on
plot(dms(2:end),MCDWM(1:end)-(IIBWMn*7 + IIAWMn*3)./(10)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax4, linspace(0.5,1,11))
ax=gca;
ax.FontSize = 15;

ax5=nexttile(tcl)
plot(dms(2:end),IIBWMt1(1:end)-(IIBWMt1n*7 + IIAWMt1n*3)./(10)+1,'r-o')
hold on
plot(dms(2:end),IIAWMt1(1:end)-(IIBWMt1n*7 + IIAWMt1n*3)./(10)+1,'b-*')
hold on
plot(dms(2:end),MCDWMt1(1:end)-(IIBWMt1n*7 + IIAWMt1n*3)./(10)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;


ax6=nexttile(tcl)
l1 = plot(dms(2:end),IIBWMt2(1:end)-(IIBWMt2n*7 + IIAWMt2n*3)./(10)+1,'r-o','DisplayName','IIB ROI')
hold on
l3 = plot(dms(2:end),IIAWMt2(1:end)-(IIBWMt2n*7 + IIAWMt2n*3)./(10)+1,'b-*','DisplayName','IIA ROI')
hold on
l5 = plot(dms(2:end),MCDWMt2(1:end)-(IIBWMt2n*7 + IIAWMt2n*3)./(10)+1,'k-+','DisplayName','mMCD ROI')
hold on 
l7 = plot(dms(2:end), ones(size(dms(2:end))), 'g-', 'DisplayName','Healthy Control')
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

hL = legend([l1,l3,l5,l7]); 
hL.Layout.Tile = 'East';




%% plot as one after normalization by a simple average

% n2a = ones(size(dms(2:end)))*n2a;
% n2b = ones(size(dms(2:end)))*n2b;
% nmcd = ones(size(dms(2:end)))*nmcd;

figure()
tcl = tiledlayout(2,3);

ax1 = nexttile(tcl)
plot(dms(2:end),IIBGM(1:end)-(IIBGMn + IIAGMn + MCDGMn)./3+1,'r-o')
hold on
plot(dms(2:end),IIAGM(1:end)-(IIBGMn + IIAGMn + MCDGMn)./3+1,'b-*')
hold on
plot(dms(2:end),MCDGM(1:end)-(IIBGMn + IIAGMn + MCDGMn)./3+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
% hold on
% plot(dms, 1-f1.a*exp(f1.b*dms))
% hold on
% plot(dms, 1-f2.a*exp(f2.b*dms))
% hold on
% plot(dms, 1-f3.a*exp(f3.b*dms))
title('GM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax1, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;


ax2 = nexttile(tcl)
plot(dms(2:end),IIBGMt1(1:end)-(IIBGMt1n + IIAGMt1n + MCDGMt1n)./3+1,'r-o')
hold on
plot(dms(2:end),IIAGMt1(1:end)-(IIBGMt1n + IIAGMt1n + MCDGMt1n)./3+1,'b-*')
hold on
plot(dms(2:end),MCDGMt1(1:end)-(IIBGMt1n + IIAGMt1n + MCDGMt1n)./3+1,'k-+')
hold on 
plot(dms(2:end),(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3-(IIBGMt1n + IIAGMt1n + MCDGMt1n)/3+1,'g-')
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(ax2, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax3=nexttile(tcl)
plot(dms(2:end),IIBGMt2(1:end)-(IIBGMt2n + IIAGMt2n + MCDGMt2n)./3+1,'r-o')
hold on
plot(dms(2:end),IIAGMt2(1:end)-(IIBGMt2n + IIAGMt2n + MCDGMt2n)./3+1,'b-*')
hold on
plot(dms(2:end),MCDGMt2(1:end)-(IIBGMt2n + IIAGMt2n + MCDGMt2n)./3+1,'k-+')
hold on 
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(ax3, linspace(0.9,1,11))
ax=gca;
ax.FontSize = 15;

ax4=nexttile(tcl)
plot(dms(2:end),IIBWM(1:end)-(IIBWMn + IIAWMn + MCDWMn)./(3)+1,'r-o')
hold on
plot(dms(2:end),IIAWM(1:end)-(IIBWMn + IIAWMn + MCDWMn)./(3)+1,'b-*')
hold on
plot(dms(2:end),MCDWM(1:end)-(IIBWMn + IIAWMn + MCDWMn)./(3)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1w/T2w', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1w/T2w value', 'FontSize', 20)
yticks(ax4, linspace(0.5,1,11))
ax=gca;
ax.FontSize = 15;

ax5=nexttile(tcl)
plot(dms(2:end),IIBWMt1(1:end)-(IIBWMt1n + IIAWMt1n + MCDWMt1n)./(3)+1,'r-o')
hold on
plot(dms(2:end),IIAWMt1(1:end)-(IIBWMt1n + IIAWMt1n + MCDWMt1n)./(3)+1,'b-*')
hold on
plot(dms(2:end),MCDWMt1(1:end)-(IIBWMt1n + IIAWMt1n + MCDWMt1n)./(3)+1,'k-+')
hold on
plot(dms(2:end),ones(size(dms(2:end))),'g-')
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;


ax6=nexttile(tcl)
l1 = plot(dms(2:end),IIBWMt2(1:end)-(IIBWMt2n + IIAWMt2n + MCDWMt2n)./(3)+1,'r-o','DisplayName','IIB ROI')
hold on
l3 = plot(dms(2:end),IIAWMt2(1:end)-(IIBWMt2n + IIAWMt2n + MCDWMt2n)./(3)+1,'b-*','DisplayName','IIA ROI')
hold on
l5 = plot(dms(2:end),MCDWMt2(1:end)-(IIBWMt2n + IIAWMt2n + MCDWMt2n)./(3)+1,'k-+','DisplayName','mMCD ROI')
hold on 
l7 = plot(dms(2:end), ones(size(dms(2:end))), 'g-', 'DisplayName','Healthy Control')
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
% yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

hL = legend([l1,l3,l5,l7]); 
hL.Layout.Tile = 'East';
%% error bar plot
% IIB group
figure()
subplot(2,3,1)
errorbar(dms(2:5),IIBGM(1:4),IIBGMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIBGM(4:end),IIBGMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIBGMn(1:4),IIBGMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIBGMn(4:end),IIBGMsen(4:end), "Color", "#0076c0")
title('GM T1/T2', 'FontSize', 20)
subtitle('Type IIB', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0.2,0.9,15))
ax=gca;
ax.FontSize = 15;

subplot(2,3,4)
errorbar(dms(2:5),IIBWM(1:4),IIBWMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIBWM(4:end),IIBWMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIBWMn(1:4),IIBWMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIBWMn(4:end),IIBWMsen(4:end), "Color", "#0076c0")

title('WM T1/T2', 'FontSize', 20)
subtitle('Type IIB', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0.2,0.9,15))
ax=gca;
ax.FontSize = 15;

subplot(2,3,2)
errorbar(dms(2:5),IIBGMt1(1:4),IIBGMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIBGMt1(4:end),IIBGMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIBGMt1n(1:4),IIBGMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIBGMt1n(4:end),IIBGMt1sen(4:end), "Color", "#0076c0")

title('GM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,5)
errorbar(dms(2:5),IIBWMt1(1:4),IIBWMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIBWMt1(4:end),IIBWMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIBWMt1n(1:4),IIBWMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIBWMt1n(4:end),IIBWMt1sen(4:end), "Color", "#0076c0")
title('WM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,3)
errorbar(dms(2:5),IIBGMt2(1:4),IIBGMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIBGMt2(4:end),IIBGMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIBGMt2n(1:4),IIBGMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIBGMt2n(4:end),IIBGMt2sen(4:end), "Color", "#0076c0")
title('GM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,6)
errorbar(dms(2:5),IIBWMt2(1:4),IIBWMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIBWMt2(4:end),IIBWMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIBWMt2n(1:4),IIBWMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIBWMt2n(4:end),IIBWMt2sen(4:end), "Color", "#0076c0")
title('WM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

% IIA group
figure()
subplot(2,3,1)
errorbar(dms(2:5),IIAGM(1:4),IIAGMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAGM(4:end),IIAGMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAGMn(1:4),IIAGMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAGMn(4:end),IIAGMsen(4:end), "Color", "#0076c0")

title('GM T1/T2', 'FontSize', 20)
subtitle('Type IIA', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,4)
errorbar(dms(2:5),IIAWM(1:4),IIAWMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAWM(4:end),IIAWMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAWMn(1:4),IIAWMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAWMn(4:end),IIAWMsen(4:end), "Color", "#0076c0")
title('WM T1/T2', 'FontSize', 20)
subtitle('Type IIA', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,2)
errorbar(dms(2:5),IIAGMt1(1:4),IIAGMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAGMt1(4:end),IIAGMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAGMt1n(1:4),IIAGMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAGMt1n(4:end),IIAGMt1sen(4:end), "Color", "#0076c0")

title('GM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,5)
errorbar(dms(2:5),IIAWMt1(1:4),IIAWMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAWMt1(4:end),IIAWMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAWMt1n(1:4),IIAWMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAWMt1n(4:end),IIAWMt1sen(4:end), "Color", "#0076c0")
title('WM T1', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,3)
errorbar(dms(2:5),IIAGMt2(1:4),IIAGMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAGMt2(4:end),IIAGMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAGMt2n(1:4),IIAGMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAGMt2n(4:end),IIAGMt2sen(4:end), "Color", "#0076c0")
title('GM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,6)
errorbar(dms(2:5),IIAWMt2(1:4),IIAWMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),IIAWMt2(4:end),IIAWMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),IIAWMt2n(1:4),IIAWMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),IIAWMt2n(4:end),IIAWMt2sen(4:end), "Color", "#0076c0")
title('WM T2', 'FontSize', 20)
subtitle(subjID(p), 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

% MCD group
figure()
subplot(2,3,1)
errorbar(dms(2:5),MCDGM(1:4),MCDGMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDGM(4:end),MCDGMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDGMn(1:4),MCDGMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDGMn(4:end),MCDGMsen(4:end), "Color", "#0076c0")
title('GM T1/T2', 'FontSize', 20)
subtitle('mMCD', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,4)
errorbar(dms(2:5),MCDWM(1:4),MCDWMse(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDWM(4:end),MCDWMse(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDWMn(1:4),MCDWMsen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDWMn(4:end),MCDWMsen(4:end), "Color", "#0076c0")
title('WM T1/T2', 'FontSize', 20)
subtitle('mMCD', 'FontSize', 18)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1/T2 value', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,2)
errorbar(dms(2:5),MCDGMt1(1:4),MCDGMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDGMt1(4:end),MCDGMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDGMt1n(1:4),MCDGMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDGMt1n(4:end),MCDGMt1sen(4:end), "Color", "#0076c0")
title('GM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,5)
errorbar(dms(2:5),MCDWMt1(1:4),MCDWMt1se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDWMt1(4:end),MCDWMt1se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDWMt1n(1:4),MCDWMt1sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDWMt1n(4:end),MCDWMt1sen(4:end), "Color", "#0076c0")
title('WM T1', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T1 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,3)
errorbar(dms(2:5),MCDGMt2(1:4),MCDGMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDGMt2(4:end),MCDGMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDGMt2n(1:4),MCDGMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDGMt2n(4:end),MCDGMt2sen(4:end), "Color", "#0076c0")
title('GM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;

subplot(2,3,6)
errorbar(dms(2:5),MCDWMt2(1:4),MCDWMt2se(1:4),"Color", "#A30234")
hold on
errorbar(dms(5:end),MCDWMt2(4:end),MCDWMt2se(4:end), "Color", "#ce8080")
hold on
errorbar(dms(2:5),MCDWMt2n(1:4),MCDWMt2sen(1:4),"Color", "#002157")
hold on
errorbar(dms(5:end),MCDWMt2n(4:end),MCDWMt2sen(4:end), "Color", "#0076c0")
title('WM T2', 'FontSize', 20)
xlabel('Percentage Distance from Center of ROI', 'FontSize', 20)
ylabel('T2 (ms)', 'FontSize', 20)
yticks(linspace(0,1,6))
ax=gca;
ax.FontSize = 15;
