%%  ROI shape dilation\
nhood = zeros(5,5,3);
% nhood(2,3,1) = 1;
% nhood(3,2:4,1) = 1;
% nhood(4,3,1) = 1;
nhood(3,3,1) = 1;
nhood(1,3,2) = 1;
nhood(2,2:4,2) = 1;
nhood(3,1:5,2) = 1;
nhood(4,2:4,2) = 1;
nhood(5,3,2) = 1;
% nhood(2,3,3) = 1;
% nhood(3,2:4,3) = 1;
% nhood(4,3,3) = 1;
nhood(3,3,3) = 1;

subjID = ["PXX_XXXXX"]; 
subjIDROI = subjID; 
types = ["IIB" "IIB" "IIA" "IIA" "MCD"...
    "mMCD" "IIA" "IIB" "IIB" "IIA" "IIB"...
    "IIB" "IIB" "MCD" "MCD"...
    "mMCD" "mMCD" "mMCD" "IIB" "mMCD" "IIB"];

MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';
ROI_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
sz = [10 8];
varTypes = ["string","string","double","double","double","double","double","double"];
varNames = ["Subject","type","GM T1/T2","WM T1/T2","GM T1","WM T1","GM T2","WM T2"];
temps = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
centroids = [];

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
    cd(path)
    t1ot2 = load_untouch_nii('T1oT2_nmCC.nii');
    image = single(t1ot2.img);
%     t1ot2 = load_untouch_nii('T1oT2_bico.nii');
%     image = single(t1ot2.img);
%     t1ot2 = load_untouch_nii('T1oT2_nmCC_bico.nii');
%     image = single(t1ot2.img);

    t1 = load_untouch_nii('MRF_T1.nii');
    t1i = single(t1.img);
    t2 = load_untouch_nii('MRF_T2.nii');
    t2i = single(t2.img);
    roi = load_untouch_nii('ROI_T1w.nii');
    roii = single(roi.img);
    % threshold lesion label
    roii(roii >= 0.9) = 1;
    roii(roii < 0.9) = 0;
    roiie = roii;


    corm = load_untouch_nii('cortical_mask_coreg_binary.nii');
    cormi = single(corm.img);
    gm = load_untouch_nii('T1_brain_pve_1.nii');
    gmi = single(gm.img);
    gmic = gmi;
    wm = load_untouch_nii('T1_brain_pve_2.nii');
    wmi = single(wm.img);
    wmic = wmi;
    brain = gmic + wmic;
    gmi((gmic>=0.5)&(brain>=0.95)) = 1;
    wmi((wmic>=0.5)&(brain>=0.95)) = 1;
    gmi((gmic>=0.5)&(brain<0.95)) = 0;
    wmi((wmic>=0.5)&(brain<0.95)) = 0;
    gmi(gmic<0.5) = 0;
    wmi(wmic<0.5) = 0;
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

    gpoints = [];
    wpoints = [];
    gpointst1 = [];
    wpointst1 = [];
    gpointst2 = [];
    wpointst2 = [];
    gr = [];
    wr = [];
    grt1 = [];
    wrt1 = [];
    grt2 = [];
    wrt2 = [];
    % keep eroding the ROI until no voxel left
    dis = 0;
    
    while sum(roiie(:)) > 10
        sum(roiie(:))
        previous = roiie;
        roiie = imerode(roiie, nhood); 
        area = previous - roiie;
        dis = dis + 1;
        
        gmall = image.*gmi.*area;
        gmall = gmall((gmall > 0)&(gmall < 8));
        gpoints = vertcat(gpoints, gmall);
        gr = vertcat(gr, ones(size(gmall))*dis);
        gmean = mean(gmall);
        gstd = std(gmall);
        gse = std(gmall)/sqrt(size(gmall, 1));

        wmall = image.*wmi.*area;
        wmall = wmall((wmall > 0)&(wmall < 8));
        wpoints = vertcat(wpoints, wmall);
        wr = vertcat(wr, ones(size(wmall))*dis);
        wmean = mean(wmall);
        wstd = std(wmall);
        wse = std(wmall)/sqrt(size(wmall, 1));

        dmeang = vertcat(dmeang, gmean);
        dstdg = vertcat(dstdg, gstd);
        dseg = vertcat(dseg, gse);
        dmeanw = vertcat(dmeanw, wmean);
        dstdw = vertcat(dstdw, wstd);
        dsew = vertcat(dsew, wse);

        gmallt1 = t1i.*gmi.*area;
        gmallt1 = gmallt1((gmallt1 > 600)&(gmallt1 < 1600));
        gpointst1 = vertcat(gpointst1, gmallt1);
        grt1 = vertcat(grt1, ones(size(gmallt1))*dis);
        gmeant1 = mean(gmallt1);
        gstdt1 = std(gmallt1);
        gset1 = std(gmallt1)/sqrt(size(gmallt1, 1));

        wmallt1 = t1i.*wmi.*area;
        wmallt1 = wmallt1((wmallt1 > 600)&(wmallt1 < 1600));
        wpointst1 = vertcat(wpointst1, wmallt1);
        wrt1 = vertcat(wrt1, ones(size(wmallt1))*dis);
        wmeant1 = mean(wmallt1);
        wstdt1 = std(wmallt1);
        wset1 = std(wmallt1)/sqrt(size(wmallt1, 1));

        dmeangt1 = vertcat(dmeangt1, gmeant1);
        dstdgt1 = vertcat(dstdgt1, gstdt1);
        dsegt1 = vertcat(dsegt1, gset1);
        dmeanwt1 = vertcat(dmeanwt1, wmeant1);
        dstdwt1 = vertcat(dstdwt1, wstdt1);
        dsewt1 = vertcat(dsewt1, wset1);

        gmallt2 = t2i.*gmi.*area;
        gmallt2 = gmallt2((gmallt2 > 20)&(gmallt2 < 60));
        gpointst2 = vertcat(gpointst2, gmallt2);
        grt2 = vertcat(grt2, ones(size(gmallt2))*dis);
        gmeant2 = mean(gmallt2);
        gstdt2 = std(gmallt2);
        gset2 = std(gmallt2)/sqrt(size(gmallt2, 1));

        wmallt2 = t2i.*wmi.*area;
        wmallt2 = wmallt2((wmallt2 > 20)&(wmallt2 < 60));
        wpointst2 = vertcat(wpointst2, wmallt2);
        wrt2 = vertcat(wrt2, ones(size(wmallt2))*dis);
        wmeant2 = mean(wmallt2);
        wstdt2 = std(wmallt2);
        wset2 = std(wmallt2)/sqrt(size(wmallt2, 1));

        dmeangt2 = vertcat(dmeangt2, gmeant2);
        dstdgt2 = vertcat(dstdgt2, gstdt2);
        dsegt2 = vertcat(dsegt2, gset2);
        dmeanwt2 = vertcat(dmeanwt2, wmeant2);
        dstdwt2 = vertcat(dstdwt2, wstdt2);
        dsewt2 = vertcat(dsewt2, wset2);
    end
% center of lesion
%     roiie = imerode(roiie, nhood);
    area = roiie;
    dis = dis + 1;

    gmall = image.*gmi.*area;
    gmall = gmall((gmall > 0)&(gmall < 8));
    gpoints = vertcat(gpoints, gmall);
    gr = vertcat(gr, ones(size(gmall))*dis);
    gmean = mean(gmall);
    gstd = std(gmall);
    gse = std(gmall)/sqrt(size(gmall, 1));

    wmall = image.*wmi.*area;
    wmall = wmall((wmall > 0)&(wmall < 8));
    wpoints = vertcat(wpoints, wmall);
    wr = vertcat(wr, ones(size(wmall))*dis);
    wmean = mean(wmall);
    wstd = std(wmall);
    wse = std(wmall)/sqrt(size(wmall, 1));

    dmeang = vertcat(dmeang, gmean);
    dstdg = vertcat(dstdg, gstd);
    dseg = vertcat(dseg, gse);
    dmeanw = vertcat(dmeanw, wmean);
    dstdw = vertcat(dstdw, wstd);
    dsew = vertcat(dsew, wse);

    gmallt1 = t1i.*gmi.*area;
    gmallt1 = gmallt1((gmallt1 > 600)&(gmallt1 < 1600));
    gpointst1 = vertcat(gpointst1, gmallt1);
    grt1 = vertcat(grt1, ones(size(gmallt1))*dis);
    gmeant1 = mean(gmallt1);
    gstdt1 = std(gmallt1);
    gset1 = std(gmallt1)/sqrt(size(gmallt1, 1));

    wmallt1 = t1i.*wmi.*area;
    wmallt1 = wmallt1((wmallt1 > 600)&(wmallt1 < 1600));
    wpointst1 = vertcat(wpointst1, wmallt1);
    wrt1 = vertcat(wrt1, ones(size(wmallt1))*dis);
    wmeant1 = mean(wmallt1);
    wstdt1 = std(wmallt1);
    wset1 = std(wmallt1)/sqrt(size(wmallt1, 1));

    dmeangt1 = vertcat(dmeangt1, gmeant1);
    dstdgt1 = vertcat(dstdgt1, gstdt1);
    dsegt1 = vertcat(dsegt1, gset1);
    dmeanwt1 = vertcat(dmeanwt1, wmeant1);
    dstdwt1 = vertcat(dstdwt1, wstdt1);
    dsewt1 = vertcat(dsewt1, wset1);

    gmallt2 = t2i.*gmi.*area;
    gmallt2 = gmallt2((gmallt2 > 20)&(gmallt2 < 60));
    gpointst2 = vertcat(gpointst2, gmallt2);
    grt2 = vertcat(grt2, ones(size(gmallt2))*dis);
    gmeant2 = mean(gmallt2);
    gstdt2 = std(gmallt2);
    gset2 = std(gmallt2)/sqrt(size(gmallt2, 1));

    wmallt2 = t2i.*wmi.*area;
    wmallt2 = wmallt2((wmallt2 > 20)&(wmallt2 < 60));
    wpointst2 = vertcat(wpointst2, wmallt2);
    wrt2 = vertcat(wrt2, ones(size(wmallt2))*dis);
    wmeant2 = mean(wmallt2);
    wstdt2 = std(wmallt2);
    wset2 = std(wmallt2)/sqrt(size(wmallt2, 1));

    dmeangt2 = vertcat(dmeangt2, gmeant2);
    dstdgt2 = vertcat(dstdgt2, gstdt2);
    dsegt2 = vertcat(dsegt2, gset2);
    dmeanwt2 = vertcat(dmeanwt2, wmeant2);
    dstdwt2 = vertcat(dstdwt2, wstdt2);
    dsewt2 = vertcat(dsewt2, wset2);

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
    gr = abs(gr-dis);
    wr = flip(wr);
    wr = abs(wr-dis);
    grt1 = flip(grt1);
    grt1 = abs(grt1-dis);
    wrt1 = flip(wrt1);
    wrt1 = abs(wrt1-dis);
    grt2 = flip(grt2);
    grt2 = abs(grt2-dis);
    wrt2 = flip(wrt2);
    wrt2 = abs(wrt2-dis);

    rea = 20;
    rea1 = 19;
    roiie = roii;
    for r = 1:rea
        previous = roiie;
        roiie = imdilate(roiie, nhood);

        area = roiie - previous;
        dis = dis + 1;
        
        gmall = image.*gmi.*area;
        gmall = gmall((gmall > 0)&(gmall < 8));
        gpoints = vertcat(gpoints, gmall);
        gr = vertcat(gr, ones(size(gmall))*dis);
        gmean = mean(gmall);
        gstd = std(gmall);
        gse = std(gmall)/sqrt(size(gmall, 1));

        wmall = image.*wmi.*area;
        wmall = wmall((wmall > 0)&(wmall < 8));
        wpoints = vertcat(wpoints, wmall);
        wr = vertcat(wr, ones(size(wmall))*dis);
        wmean = mean(wmall);
        wstd = std(wmall);
        wse = std(wmall)/sqrt(size(wmall, 1));

        dmeang = vertcat(dmeang, gmean);
        dstdg = vertcat(dstdg, gstd);
        dseg = vertcat(dseg, gse);
        dmeanw = vertcat(dmeanw, wmean);
        dstdw = vertcat(dstdw, wstd);
        dsew = vertcat(dsew, wse);

        gmallt1 = t1i.*gmi.*area;
        gmallt1 = gmallt1((gmallt1 > 600)&(gmallt1 < 1600));
        gpointst1 = vertcat(gpointst1, gmallt1);
        grt1 = vertcat(grt1, ones(size(gmallt1))*dis);
        gmeant1 = mean(gmallt1);
        gstdt1 = std(gmallt1);
        gset1 = std(gmallt1)/sqrt(size(gmallt1, 1));

        wmallt1 = t1i.*wmi.*area;
        wmallt1 = wmallt1((wmallt1 > 600)&(wmallt1 < 1600));
        wpointst1 = vertcat(wpointst1, wmallt1);
        wrt1 = vertcat(wrt1, ones(size(wmallt1))*dis);
        wmeant1 = mean(wmallt1);
        wstdt1 = std(wmallt1);
        wset1 = std(wmallt1)/sqrt(size(wmallt1, 1));

        dmeangt1 = vertcat(dmeangt1, gmeant1);
        dstdgt1 = vertcat(dstdgt1, gstdt1);
        dsegt1 = vertcat(dsegt1, gset1);
        dmeanwt1 = vertcat(dmeanwt1, wmeant1);
        dstdwt1 = vertcat(dstdwt1, wstdt1);
        dsewt1 = vertcat(dsewt1, wset1);

        gmallt2 = t2i.*gmi.*area;
        gmallt2 = gmallt2((gmallt2 > 20)&(gmallt2 < 60));
        gpointst2 = vertcat(gpointst2, gmallt2);
        grt2 = vertcat(grt2, ones(size(gmallt2))*dis);
        gmeant2 = mean(gmallt2);
        gstdt2 = std(gmallt2);
        gset2 = std(gmallt2)/sqrt(size(gmallt2, 1));

        wmallt2 = t2i.*wmi.*area;
        wmallt2 = wmallt2((wmallt2 > 20)&(wmallt2 < 60));
        wpointst2 = vertcat(wpointst2, wmallt2);
        wrt2 = vertcat(wrt2, ones(size(wmallt2))*dis);
        wmeant2 = mean(wmallt2);
        wstdt2 = std(wmallt2);
        wset2 = std(wmallt2)/sqrt(size(wmallt2, 1));

        dmeangt2 = vertcat(dmeangt2, gmeant2);
        dstdgt2 = vertcat(dstdgt2, gstdt2);
        dsegt2 = vertcat(dsegt2, gset2);
        dmeanwt2 = vertcat(dmeanwt2, wmeant2);
        dstdwt2 = vertcat(dstdwt2, wstdt2);
        dsewt2 = vertcat(dsewt2, wset2);
    
    
    end
    dvarg = dstdg./dmeang;
    dvarw = dstdw./dmeanw;
    
    fir3g = gpoints((gr>0)&(gr<4));
    fir3w = wpoints((wr>0)&(wr<4));
    fir3gt1 = gpointst1((grt1>0)&(grt1<4));
    fir3wt1 = wpointst1((wrt1>0)&(wrt1<4));
    fir3gt2 = gpointst2((grt2>0)&(grt2<4));
    fir3wt2 = wpointst2((wrt2>0)&(wrt2<4));
    
    las3g = gpoints((gr>dis-3)&(gr==dis));
    las3w = wpoints((wr>dis-3)&(wr==dis));
    las3gt1 = gpointst1((grt1>dis-3)&(grt1==dis));
    las3wt1 = wpointst1((wrt1>dis-3)&(wrt1==dis));
    las3gt2 = gpointst2((grt2>dis-3)&(grt2==dis));
    las3wt2 = wpointst2((wrt2>dis-3)&(wrt2==dis));
    
    [hg,pg] = ttest2(fir3g, las3g, 'Vartype', 'unequal');
    [hw,pw] = ttest2(fir3w, las3w, 'Vartype', 'unequal');
    [hgt1,pgt1] = ttest2(fir3gt1, las3gt1, 'Vartype', 'unequal');
    [hwt1,pwt1] = ttest2(fir3wt1, las3wt1, 'Vartype', 'unequal');
    [hgt2,pgt2] = ttest2(fir3gt2, las3gt2, 'Vartype', 'unequal');
    [hwt2,pwt2] = ttest2(fir3wt2, las3wt2, 'Vartype', 'unequal');
    
    temps(p,:) = {subjID(p), types(p), hg, hw, hgt1, hwt1, hgt2, hwt2};
    
% separate plots
    figure()
    subplot(2,3,1)
    errorbar([1:1:dis-rea],dmeang(1:1:dis-rea),dseg(1:1:dis-rea))
    hold on
    errorbar([dis-rea1:1:dis],dmeang(dis-rea1:1:dis),dseg(dis-rea1:1:dis))
    title('GM T1w/T2w', 'FontSize', 20)
    subtitle(strcat(subjID(p)," ", types(p)), 'FontSize', 18)
    xlabel('distance from center (0.94mm)', 'FontSize', 20)
    ylabel('T1w/T2w value', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
    xticks(0:2:dis)
    ylim([0.2 0.6])

    subplot(2,3,4)
    errorbar([1:1:dis-rea],dmeanw(1:1:dis-rea),dsew(1:1:dis-rea))
    hold on
    errorbar([dis-rea1:1:dis],dmeanw(dis-rea1:1:dis),dsew(dis-rea1:1:dis))
    title('WM T1w/T2w', 'FontSize', 20)
    subtitle(strcat(subjID(p)," ", types(p)), 'FontSize', 18)
    xlabel('distance from center (0.94mm)', 'FontSize', 20)
    ylabel('T1w/T2w value', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
    xticks(0:2:dis)
    ylim([0.4 0.8])

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
    errorbar([1:1:dis-rea],dmeangt1(1:1:dis-rea),dsegt1(1:1:dis-rea))
    hold on
    errorbar([dis-rea1:1:dis],dmeangt1(dis-rea1:1:dis),dsegt1(dis-rea1:1:dis))
    title('GM T1', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('distance from center (0.94mm)', 'FontSize', 20)
    ylabel('T1 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
    xticks(0:2:dis)
    ylim([1200 1600])

    subplot(2,3,5)
    errorbar([1:1:dis-rea],dmeanwt1(1:1:dis-rea),dsewt1(1:1:dis-rea))
    hold on
    errorbar([dis-rea1:1:dis],dmeanwt1(dis-rea1:1:dis),dsewt1(dis-rea1:1:dis))
    title('WM T1', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('distance from center (0.94mm)', 'FontSize', 20)
    ylabel('T1 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
    xticks(0:2:dis)
    

    subplot(2,3,3)
    errorbar([1:1:dis-rea],dmeangt2(1:1:dis-rea),dsegt2(1:1:dis-rea))
    hold on
    errorbar([dis-rea1:1:dis],dmeangt2(dis-rea1:1:dis),dsegt2(dis-rea1:1:dis))
    title('GM T2', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('distance from center (0.94mm)', 'FontSize', 20)
    ylabel('T2 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
    xticks(0:2:dis)
    ylim([48 68])

    subplot(2,3,6)
    errorbar([1:1:dis-rea],dmeanwt2(1:1:dis-rea),dsewt2(1:1:dis-rea))
    hold on
    errorbar([dis-rea1:1:dis],dmeanwt2(dis-rea1:1:dis),dsewt2(dis-rea1:1:dis))
    title('WM T2', 'FontSize', 20)
    subtitle(subjID(p), 'FontSize', 18)
    xlabel('distance from center (0.94mm)', 'FontSize', 20)
    ylabel('T2 (ms)', 'FontSize', 20)
    ax=gca;
    ax.FontSize = 15;
    xticks(0:2:dis)
    legend('within lesion','outside lesion');
%     figure()
%     boxplot(gpoints, gr)
%     boxplot(wpoints, wr)
    
end

%% Generate Lesion Masks
nhood = zeros(5,5,3);
% nhood(2,3,1) = 1;
% nhood(3,2:4,1) = 1;
% nhood(4,3,1) = 1;
nhood(3,3,1) = 1;
nhood(1,3,2) = 1;
nhood(2,2:4,2) = 1;
nhood(3,1:5,2) = 1;
nhood(4,2:4,2) = 1;
nhood(5,3,2) = 1;
% nhood(2,3,3) = 1;
% nhood(3,2:4,3) = 1;
% nhood(4,3,3) = 1;
nhood(3,3,3) = 1;

subjID=["study13375" "study13702" "study13890" "study14364" "study13509"...
    "study14203" "study14218" "study14439" "study14473" "study14516" "study14705"...
    "study13923" "study14129" "study14287" "study14494"];
subjIDROI = ["P31_13375" "P42_13702" "P50_13890" "P72_14364"...
    "P33_13509" "P61_14203" "P63_14218" "P76_14439" "P81_14473" "P83_14516"...
    "P85_14705" "P52_13923" "P57_14129" "P68_14287" "P82_14494" ...
    "P102_15731" "P91_15140" "P100_15582" "P101_15675" "P89_14878" "P46_13728"]; 
types = ["IIB" "IIB" "IIA" "IIA"...
    "MCD" "MCD" "IIA" "IIB" "IIB" "IIA"...
    "IIB" "IIB" "IIB" "MCD" "MCD"];

subjID=["study13375" "study13702" "study13890" "study14364" "study13509"];
subjIDROI = ["P31_13375" "P42_13702" "P50_13890" "P72_14364"...
    "P33_13509"]; 
types = ["IIB" "IIB" "IIA" "IIA"...
    "MCD"];


MRF_path='T:\Imaging\Multimodal\Myelin\Patients';
ROI_path = 'T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
centroids = [];

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
    cd(path)
    mkdir masks
    copyfile 'T1overT2.nii' masks
    cd(strcat(MRF_path,'\',subjID(p),'\masks'))
    t1ot2 = load_untouch_nii('T1overT2.nii');
    image = single(t1ot2.img);
    image = image.*0;
    t1ot2.img = image;
    save_untouch_nii(t1ot2, 'T1overT2.nii');
    cd(path)
    
    roi = load_untouch_nii('ROI_T1w.nii');
    roii = single(roi.img);

    roii(roii >= 0.9) = 1;
    roii(roii < 0.9) = 0;
    roiie = roii;

    corm = load_untouch_nii('cortical_mask_coreg_binary.nii');
    cormi = single(corm.img);
    gm = load_untouch_nii('T1_brain_pve_1.nii');
    gmi = single(gm.img);
    gmic = gmi;
    wm = load_untouch_nii('T1_brain_pve_2.nii');
    wmi = single(wm.img);
    wmic = wmi;
    brain = gmic + wmic;
    gmi((gmic>=0.5)&(brain>=0.95)) = 1;
    wmi((wmic>=0.5)&(brain>=0.95)) = 1;
    gmi((gmic>=0.5)&(brain<0.95)) = 0;
    wmi((wmic>=0.5)&(brain<0.95)) = 0;
    gmi(gmic<0.5) = 0;
    wmi(wmic<0.5) = 0;

    gmi = gmi.*cormi;
    wmi = wmi.*cormi;

    % keep eroding the ROI until no voxel left
    dis = 0;
    while sum(roiie(:)) > 1
        previous = roiie;
        roiie = imerode(roiie, nhood);
        area = previous - roiie;
        dis = dis + 1;
        
        gmall = gmi.*area*1;
        wmall = wmi.*area*2;
%         figure()
%         histogram(gmall,50)
%         figure()
%         histogram(wmall,50)

        cd(strcat(path,'\masks'))
        copyfile('T1overT2.nii', strcat('g', string(-dis),'.nii'))
        copyfile('T1overT2.nii', strcat('w', string(-dis),'.nii'))
        
        maskg = load_untouch_nii(char(strcat('g', string(-dis),'.nii')))
        maskg.img = gmall;
        save_untouch_nii(maskg, char(strcat('g', string(-dis),'.nii')))
        
        maskw = load_untouch_nii(char(strcat('w', string(-dis),'.nii')))
        maskw.img = wmall;
        save_untouch_nii(maskw, char(strcat('w', string(-dis),'.nii')))

        
    end

  
    roiie = roii;
    for r = 1:20
        previous = roiie;
        roiie = imdilate(roiie, nhood);
        area = roiie - previous;
        dis = dis + 1;
        
        gmall = gmi.*area*3;
        wmall = wmi.*area*4;
    
        cd(strcat(MRF_path,'\',subjID(p),'\masks'))
        copyfile('T1overT2.nii', char(strcat('g', string(dis),'.nii')))
        copyfile('T1overT2.nii', char(strcat('w', string(dis),'.nii')))
        
        maskg = load_untouch_nii(char(strcat('g', string(dis),'.nii')))
        maskg.img = gmall;
        save_untouch_nii(maskg, char(strcat('g', string(dis),'.nii')))
        
        maskw = load_untouch_nii(char(strcat('w', string(dis),'.nii')))
        maskw.img = wmall;
        save_untouch_nii(maskw, char(strcat('w', string(dis),'.nii')))    
    end

end
%% single ROI mask

nhood = zeros(5,5,3);
% nhood(2,3,1) = 1;
% nhood(3,2:4,1) = 1;
% nhood(4,3,1) = 1;
nhood(3,3,1) = 1;
nhood(1,3,2) = 1;
nhood(2,2:4,2) = 1;
nhood(3,1:5,2) = 1;
nhood(4,2:4,2) = 1;
nhood(5,3,2) = 1;
% nhood(2,3,3) = 1;
% nhood(3,2:4,3) = 1;
% nhood(4,3,3) = 1;
nhood(3,3,3) = 1;

subjID=["study13375" "study13702" "study13890" "study14364" "study13509"...
    "study14203" "study14218" "study14439" "study14473" "study14516" "study14705"...
    "study13923" "study14129" "study14287" "study14494"];
subjIDROI = ["P31_13375" "P42_13702" "P50_13890" "P72_14364"...
    "P33_13509" "P61_14203" "P63_14218" "P76_14439" "P81_14473" "P83_14516"...
    "P85_14705" "P52_13923" "P57_14129" "P68_14287" "P82_14494"]; 
types = ["IIB" "IIB" "IIA" "IIA"...
    "MCD" "MCD" "IIA" "IIB" "IIB" "IIA"...
    "IIB" "IIB" "IIB" "MCD" "MCD"];


subjID=["study14203" "study14218" "study14439" "study14473" "study14516" "study14705"...
    "study13923" "study14129" "study14287" "study14494"];
subjIDROI = ["P61_14203" "P63_14218" "P76_14439" "P81_14473" "P83_14516"...
    "P85_14705" "P52_13923" "P57_14129" "P68_14287" "P82_14494"]; 
types = ["MCD" "IIA" "IIB" "IIB" "IIA"...
    "IIB" "IIB" "IIB" "MCD" "MCD"];

MRF_path='T:\Imaging\Multimodal\Myelin\Patients';
ROI_path = 'T:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';
centroids = [];

for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
    cd(path)
    mkdir masks
    copyfile 'T1overT2.nii' masks
    cd(strcat(MRF_path,'\',subjID(p),'\masks'))
    t1ot2 = load_untouch_nii('T1overT2.nii');
    image = single(t1ot2.img);
    image = image.*0;
    t1ot2.img = image;
    save_untouch_nii(t1ot2, 'T1overT2.nii');
    cd(path)
    
    roi = load_untouch_nii('ROI_T1w.nii');
    roii = single(roi.img);

    roii(roii >= 0.9) = 1;
    roii(roii < 0.9) = 0;
    roiie = roii;

    corm = load_untouch_nii('cortical_mask_coreg_binary.nii');
    cormi = single(corm.img);
    gm = load_untouch_nii('T1_brain_pve_1.nii');
    gmi = single(gm.img);
    gmic = gmi;
    wm = load_untouch_nii('T1_brain_pve_2.nii');
    wmi = single(wm.img);
    wmic = wmi;
    csf = load_untouch_nii('T1_brain_pve_0.nii');
    csfi = single(csf.img);
    gmi((gmic>=0.5)&(csfi<=0.1)) = 1;
    wmi((wmic>=0.5)&(csfi<=0.1)) = 1;
    gmi((gmic>=0.5)&(csfi>0.1)) = 0;
    wmi((wmic>=0.5)&(csfi>0.1)) = 0;
    gmi(gmic<0.5) = 0;
    wmi(wmic<0.5) = 0;

    gmi = gmi.*cormi;
    wmi = wmi.*cormi;

    % keep eroding the ROI until no voxel left
    cd(strcat(MRF_path,'\',subjID(p),'\masks'))
    copyfile('T1overT2.nii', 'the1mask.nii')
   
    dis = 0;
    while sum(roiie(:)) > 1
        previous = roiie;
        roiie = imerode(roiie, nhood);
        area = previous - roiie;
        dis = dis + 1;
        
        gms = gmi.*area*1;
        wms = wmi.*area*2;
        
        maskg = load_untouch_nii('the1mask.nii');
        maskg.img = maskg.img + gms + wms;
        save_untouch_nii(maskg, 'the1mask.nii')   
    end

  
    roiie = roii;
    for r = 1:20
        previous = roiie;
        roiie = imdilate(roiie, nhood);
        area = roiie - previous;
        dis = dis + 1;
        
        gms = gmi.*area*3;
        wms = wmi.*area*4;
        maskg = load_untouch_nii('the1mask.nii');
        maskg.img = maskg.img + gms + wms;
        save_untouch_nii(maskg, 'the1mask.nii')    
    end

end
%%
subjID=["PXX_XXXXX"];
subjIDROI = subjID;  
MRF_path='Z:\Imaging\Multimodal\Myelin\Patients';
ROI_path = 'Z:\Imaging\Multimodal\MRF\Recon_MRF_3T\Patients';


for p = 1:length(subjID)
    subjID(p)
    path = strcat(MRF_path,'\',subjID(p));
%     path = strcat(MRF_path,'/',p,'/MRF_VBM');
    cd(path)
    
    t1ot2 = load_untouch_nii('T1overT2.nii');
    image = single(t1ot2.img);
    roi = load_untouch_nii('ROI_T1w.nii');
    roii = single(roi.img);
    roii(roii>0) = 1;
    [roiblob, nBlobs] = bwlabeln(roii);
    centroid = regionprops3(roiblob, "Centroid");
    c = round(centroid.Centroid);
    cx = c(2);
    cy = c(1);
    cz = c(3);

    t1mrf = load_untouch_nii('MRF_T1.nii');
    t1mrfi = single(t1mrf.img);
    t2mrf = load_untouch_nii('MRF_T2.nii');
    t2mrfi = single(t2mrf.img);
    gm = load_untouch_nii('T1_brain_pve_1.nii');
    gmi = single(gm.img);
    gmi(gmi>=0.8) = 1;
    gmi(gmi<=0.8) = 0;
    t1mrfig = t1mrfi.*gmi;
    t2mrfig = t2mrfi.*gmi;
    wm = load_untouch_nii('T1_brain_pve_2.nii');
    wmi = single(wm.img);
    wmi(wmi>=0.8) = 1;
    wmi(wmi<=0.8) = 0;
    t1mrfiw = t1mrfi.*wmi;
    t2mrfiw = t2mrfi.*wmi;
    gmi = image.*gmi;
    wmi = image.*wmi;
    ginroi = gmi.*roii;
    winroi = wmi.*roii;
    meanROIg = mean(ginroi(ginroi>0));
    meanROIw = mean(winroi(winroi>0));
   
    dmeang = [];
    dmeanw = [];
    dstdg = [];
    dstdw = [];
    gpoints = [];
    wpoints = [];
    gpointsx = [];
    wpointsx = [];
    dmeangnorm = [];
    dstdgnorm = [];
    dmeanwnorm = [];
    dstdwnorm = [];
    dmeangt1 = [];
    dstdgt1 = [];
    dmeanwt1 = [];
    dstdwt1 = [];
    dmeangnormt1 = [];
    dstdgnormt1 = [];
    dmeanwnormt1 = [];
    dstdwnormt1 = [];
    dmeangt2 = [];
    dstdgt2 = [];
    dmeanwt2 = [];
    dstdwt2 = [];
    dmeangnormt2 = [];
    dstdgnormt2 = [];
    dmeanwnormt2 = [];
    dstdwnormt2 = [];

    for r = 1:2:37
        x = -38:1:38;
        y = -38:1:38;
        z = -38:2:38;
        [x y z] = meshgrid(x,y,z);
        sphere = x.^ 2 + y.^ 2 + z.^2;
        region = 1*logical((sphere <= r ^ 3));
        if r > 1
            region = 1*logical((sphere <= r ^ 3)) - 1*logical((sphere <= (r-2) ^ 3));
        end
        gmall = gmi(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19);
        gmean = mean(gmall((gmall > 0)&(gmall < 15)));
        gstd = std(gmall((gmall > 0)&(gmall < 15)))/sqrt(size(gmall,1));
        gpoints = vertcat(gpoints, gmall((gmall > 0)&(gmall < 15)));
        gpointsx = vertcat(gpointsx, r*ones(size(gmall((gmall > 0)&(gmall < 15)))));
        
        gmnorm = gmi(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*abs(roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19)-1);
        gnmean = mean(gmnorm((gmnorm > 0)&(gmnorm < 15)));
        gnstd = std(gmnorm((gmnorm > 0)&(gmnorm < 15)))/sqrt(size(gmnorm,1));

        wmall = wmi(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19);
        wmean = mean(wmall((wmall > 0)&(wmall < 15)));
        wstd = std(wmall((wmall > 0)&(wmall < 15)))/sqrt(size(wmall,1));
        wpoints = vertcat(wpoints, wmall((wmall > 0)&(wmall < 15)));
        wpointsx = vertcat(wpointsx, r*ones(size(wmall((wmall > 0)&(wmall < 15)))));

        wmnorm = wmi(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*abs(roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19)-1);
        wnmean = mean(wmnorm((wmnorm > 0)&(wmnorm < 15)));
        wnstd = std(wmnorm((wmnorm > 0)&(wmnorm < 15)))/sqrt(size(wmnorm,1));

        dmeang = vertcat(dmeang, gmean);
        dstdg = vertcat(dstdg, gstd);
        dmeanw = vertcat(dmeanw, wmean);
        dstdw = vertcat(dstdw, wstd);

        dmeangnorm = vertcat(dmeangnorm, gnmean);
        dstdgnorm = vertcat(dstdgnorm, gnstd);
        dmeanwnorm = vertcat(dmeanwnorm, wnmean);
        dstdwnorm = vertcat(dstdwnorm, wnstd);

        % MRF T1 
        gmallt1 = t1mrfig(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19);
        gmeant1 = mean(gmallt1((gmallt1 > 600)&(gmallt1 < 1600)));
        gstdt1 = std(gmallt1((gmallt1 > 600)&(gmallt1 < 1600)))/sqrt(size(gmallt1,1));
        
        gmnormt1 = t1mrfig(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*abs(roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19)-1);
        gnmeant1 = mean(gmnormt1((gmnormt1 > 600)&(gmnormt1 < 1600)));
        gnstdt1 = std(gmnormt1((gmnormt1 > 600)&(gmnormt1 < 1600)))/sqrt(size(gmnormt1,1));

        wmallt1 = t1mrfiw(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19);
        wmeant1 = mean(wmallt1((wmallt1 > 600)&(wmallt1 < 1600)));
        wstdt1 = std(wmallt1((wmallt1 > 600)&(wmallt1 < 1600)))/sqrt(size(wmallt1,1));

        wmnormt1 = t1mrfiw(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*abs(roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19)-1);
        wnmeant1 = mean(wmnormt1((wmnormt1 > 600)&(wmnormt1 < 1600)));
        wnstdt1 = std(wmnormt1((wmnormt1 > 600)&(wmnormt1 < 1600)))/sqrt(size(wmnormt1,1));

        dmeangt1 = vertcat(dmeangt1, gmeant1);
        dstdgt1 = vertcat(dstdgt1, gstdt1);
        dmeanwt1 = vertcat(dmeanwt1, wmeant1);
        dstdwt1 = vertcat(dstdwt1, wstdt1);

        dmeangnormt1 = vertcat(dmeangnormt1, gnmeant1);
        dstdgnormt1 = vertcat(dstdgnormt1, gnstdt1);
        dmeanwnormt1 = vertcat(dmeanwnormt1, wnmeant1);
        dstdwnormt1 = vertcat(dstdwnormt1, wnstdt1);

        % MRF T2
        gmallt2 = t2mrfig(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19);
        gmeant2 = mean(gmallt2((gmallt2 > 20)&(gmallt2 < 60)));
        gstdt2 = std(gmallt2((gmallt2 > 20)&(gmallt2 < 60)))/sqrt(size(gmallt2,1));
        
        gmnormt2 = t2mrfig(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*abs(roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19)-1);
        gnmeant2 = mean(gmnormt2((gmnormt2 > 20)&(gmnormt2 < 60)));
        gnstdt2 = std(gmnormt2((gmnormt2 > 20)&(gmnormt2 < 60)))/sqrt(size(gmnormt2,1));

        wmallt2 = t2mrfiw(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19);
        wmeant2 = mean(wmallt2((wmallt2 > 20)&(wmallt2 < 60)));
        wstdt2 = std(wmallt2((wmallt2 > 20)&(wmallt2 < 60)))/sqrt(size(wmallt2,1));

        wmnormt2 = t2mrfiw(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19).*region.*abs(roii(cx-38:cx+38,cy-38:cy+38,cz-19:cz+19)-1);
        wnmeant2 = mean(wmnormt2((wmnormt2 > 20)&(wmnormt2 < 60)));
        wnstdt2 = std(wmnormt2((wmnormt2 > 20)&(wmnormt2 < 60)))/sqrt(size(wmnormt2,1));

        dmeangt2 = vertcat(dmeangt2, gmeant2);
        dstdgt2 = vertcat(dstdgt2, gstdt2);
        dmeanwt2 = vertcat(dmeanwt2, wmeant2);
        dstdwt2 = vertcat(dstdwt2, wstdt2);

        dmeangnormt2 = vertcat(dmeangnormt2, gnmeant2);
        dstdgnormt2 = vertcat(dstdgnormt2, gnstdt2);
        dmeanwnormt2 = vertcat(dmeanwnormt2, wnmeant2);
        dstdwnormt2 = vertcat(dstdwnormt2, wnstdt2);

    end
        dvarg = dstdg./dmeang;
        dvarw = dstdw./dmeanw;
    
    dmeang(7:9) = NaN;
    dstdg(7:9) = NaN;
    dmeangt1(7:9) = NaN;
    dstdgt1(7:9) = NaN;
    dmeangt2(7:9) = NaN;
    dstdgt2(7:9) = NaN;
    overlaps = isnan(dmeang) + isnan(dmeangnorm);
    inds = find(overlaps == 0);
    dmeangnorm(inds) = NaN;
    dstdgnorm(inds) = NaN;
    dmeangnormt1(inds) = NaN;
    dstdgnormt1(inds) = NaN;
    dmeangnormt2(inds) = NaN;
    dstdgnormt2(inds) = NaN;
    
    dmeanw(7:9) = NaN;
    dstdw(7:9) = NaN;
    dmeanwt1(7:9) = NaN;
    dstdwt1(7:9) = NaN;
    dmeanwt2(7:9) = NaN;
    dstdwt2(7:9) = NaN;
    overlaps = isnan(dmeanw) + isnan(dmeanwnorm);
    inds = find(overlaps == 0);
    dmeanwnorm(inds) = NaN;
    dstdwnorm(inds) = NaN;
    dmeanwnormt1(inds) = NaN;
    dstdwnormt1(inds) = NaN;
    dmeanwnormt2(inds) = NaN;
    dstdwnormt2(inds) = NaN;

    figure()
    subplot(2,4,1)
%     scatter([1:1:25],dmeang)
    errorbar([1:1:19],dmeang,dstdg)
    hold on
    errorbar([1:1:19],dmeangnorm,dstdgnorm)
%     scatter(gpointsx,gpoints)
    title('GM mean & SD')
    subtitle(strcat(subjIDROI(p), 'mean in ROI = ', string(meanROIg)))
    xlabel('distance from center (0.94mm)')
    ylabel('T1/T2 value')

    subplot(2,4,5)
    errorbar([1:1:19],dmeanw,dstdw)
    hold on
    errorbar([1:1:19],dmeanwnorm,dstdwnorm)
    title('WM  mean & SD')
    subtitle(strcat(subjIDROI(p), 'mean in ROI = ', string(meanROIw)))
    xlabel('distance from center (0.94mm)')
    ylabel('T1/T2 value')

    subplot(2,4,2)
    scatter([1:1:19],dvarg)
    title('GM CV')
    subtitle(strcat('mean in ROI = ', string(meanROIg)))
    xlabel('distance from center (0.94mm)')
    ylabel('T1/T2 CV')

    subplot(2,4,6)
    scatter([1:1:19],dvarw)
    title('WM CV')
    subtitle(strcat('mean in ROI = ', string(meanROIw)))
    xlabel('distance from center (0.94mm)')
    ylabel('T1/T2 CV')

    subplot(2,4,3)
    errorbar([1:1:19],dmeangt1,dstdgt1)
    hold on
    errorbar([1:1:19],dmeangnormt1,dstdgnormt1)
    title('GM T1')
    subtitle(strcat(subjIDROI(p), 'mean in ROI = ', string(meanROIg)))
    xlabel('distance from center (0.94mm)')
    ylabel('T1 value')

    subplot(2,4,7)
    errorbar([1:1:19],dmeanwt1,dstdwt1)
    hold on
    errorbar([1:1:19],dmeanwnormt1,dstdwnormt1)
    title('WM T1')
    subtitle(strcat(subjIDROI(p), 'mean in ROI = ', string(meanROIw)))
    xlabel('distance from center (0.94mm)')
    ylabel('T1 value')

    subplot(2,4,4)
    errorbar([1:1:19],dmeangt2,dstdgt2)
    hold on
    errorbar([1:1:19],dmeangnormt2,dstdgnormt2)
    title('GM T2')
    subtitle(strcat(subjIDROI(p), 'mean in ROI = ', string(meanROIg)))
    xlabel('distance from center (0.94mm)')
    ylabel('T2 value')

    subplot(2,4,8)
    errorbar([1:1:19],dmeanwt2,dstdwt2)
    hold on
    errorbar([1:1:19],dmeanwnormt2,dstdwnormt2)
    title('WM T2')
    subtitle(strcat(subjIDROI(p), 'mean in ROI = ', string(meanROIw)))
    xlabel('distance from center (0.94mm)')
    ylabel('T2 value')

end