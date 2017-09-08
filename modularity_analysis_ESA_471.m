% This script messes with modularity values/ setups as used by Andric &
% Hasson 2015 and Cocchi et al., 2015.
% uses a BaseLine Modularity defined with the code commented out below.

clear
do_makevals = 1;
do_theshold_matrices = 0;
do_create_mod = 1;
do_WMD = 1;
do_degree = 1;
do_PC = 1;
do_WD = 1;
do_BMD_Module_specific = 1;
do_BMD_not_specific = 1;
do_controllability = 1;
do_stats = 1;

thr = 0.4;
beta = 1;
modular_type= 5; % 1 = baselinemod, 2 = baseline fmri (BCT algorithm), 3 = baseline fmri (BLUE algorithm), 4 = dti (BCT algorithm), 5 = dti (BLUE algorithm), 6 = a priori DTI mod
gamma = 1;

subjects = {'002_TP', '003_RL', '004_RS', '005_TWS', '006_GH', '007_VD', '009_CY', '010_MU', '011_LD', '012_SD', '013_RH', '014_JM', '015_BM'};
% load /Users/swd4/Dropbox/duke/TMSEnhance/Analysis/BrainNet_files/Jun18_lower_gamma/baselinemod_noCB.mat
% load /Users/swd4/Dropbox/Duke/TMSEnhance/Connectomes/dtiMOD_mod.mat
gender = [1	1	2	1	1	2	2	2	2	1	2	2];
age = [65.58	66.76	63.75	63.75	61.55	62.96	65.59	72.73	65.61	68.60	74.31	69.52];

if do_makevals == 1
    % set shit up for avg matrix, pull degree while we're there
    for i = 1:12
        baselineHitfile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s__BASE_ESA.txt', subjects{i});
        exciteHitfile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s__5Hz_ESA.txt', subjects{i});
        inhibHitfile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s__1Hz_ESA.txt', subjects{i});
        dtifile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s_HOAsp_connectome.csv', subjects{i});
        baseHit = dlmread(baselineHitfile); baseHit(isnan(baseHit)) = 1;
        exciteHit = dlmread(exciteHitfile);
        inhibHit = dlmread(inhibHitfile);
        dti = dlmread(dtifile);
        
        if do_theshold_matrices == 1
            %             base(abs(base)<thr) = NaN;
            %             excite(abs(excite)<thr) = NaN;
            %             inhib(abs(inhib)<thr) = NaN;
            %             dti(abs(dti)<thr) = NaN;
            baseHit = ((baseHit+1)/2)^beta;
            exciteHit = ((exciteHit+1)/2)^beta;
            inhibHit = ((inhibHit+1)/2)^beta;
            dti = ((dti+1)/2)^beta;
        end
        
        base_all(i, :,:) = baseHit;
        excite_all(i,:,:) = exciteHit;
        inhib_all(i,:,:) = inhibHit;
        dti_all(i,:,:) = dti;
    end
    
    % create the averages
    base_avg = squeeze(nanmean(base_all, 1));
    excite_avg = squeeze(nanmean(excite_all, 1));
    inhib_avg = squeeze(nanmean(inhib_all, 1));
    dti_avg = squeeze(nanmean(dti_all, 1));
    
    if do_create_mod == 1
        %% create modularity from group avg; first get some iterations
        for modcounter = 1:100
            %             [modularity.values.baseline(modcounter,:), modularity.Q.baseline(modcounter)] = community_louvain(base_avg,0.7);
            %             [modularity.values.excite(modcounter,:), modularity.Q.excite(modcounter)] = community_louvain(excite_avg,0.7);
            %             [modularity.values.inhib(modcounter,:), modularity.Q.inhib(modcounter)] = community_louvain(inhib_avg,0.7);
            [modularity.values.dti(modcounter,:), modularity.Q.dti(modcounter)] = community_louvain(dti_avg, gamma);
        end
        
        % find the best modularity out of the 100 tries
        %         [valueB, indexB]= max(modularity.Q.baseline');
        %         [valueE, indexE]= max(modularity.Q.excite');
        %         [valueI, indexI]= max(modularity.Q.inhib');
        [valueD, indexD]= max(modularity.Q.dti');
    end
    
    if modular_type == 1
        Ci = baselinemod;
    elseif modular_type == 2
        Ci = modularity.values.baseline(indexB,:)';
    elseif modular_type == 3
        [Ci_base,Q_base,B_base] = Mod_BLUE_Oct5(base_avg, 1,1);
        Ci = Ci_base;
    elseif modular_type == 4
        Ci = modularity.values.dti(indexD,:);
    elseif modular_type == 5
        [Ci_dti,Q_dti,B_dti] = Mod_BLUE_Oct5(dti_avg, 1,1);
        Ci = Ci_dti;
    elseif modular_type == 6
        Ci = Ci;
    end
    
    clear base excite inhib dti i modcounter
    
    for i = 1:12
        % load everything again
        baselineHitfile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s__BASE_ESA.txt', subjects{i});
        exciteHitfile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s__5Hz_ESA.txt', subjects{i});
        inhibHitfile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s__1Hz_ESA.txt', subjects{i});
        dtifile = sprintf('/Users/swd4/Dropbox/duke/TMSEnhance/Connectomes471/%s_HOAsp_connectome.csv', subjects{i});
        
        baseHit = dlmread(baselineHitfile); baseHit(isnan(baseHit)) = 1;
        exciteHit = dlmread(exciteHitfile);
        inhibHit = dlmread(inhibHitfile);
        dti = dlmread(dtifile);
        
        if do_theshold_matrices == 1
            %             base(abs(base)<thr) = NaN;
            %             excite(abs(excite)<thr) = NaN;
            %             inhib(abs(inhib)<thr) = NaN;
            %             dti(abs(dti)<thr) = NaN;
            baseHit = ((baseHit+1)/2)^beta;
            exciteHit = ((exciteHit+1)/2)^beta;
            inhibHit = ((inhibHit+1)/2)^beta;
            dti = ((dti+1)/2)^beta;
        end
        
        if do_degree == 1
            %% degree
            degree.baseline(i,:) = strengths_und(baseHit);
            degree.inhib(i,:) = strengths_und(inhibHit);
            degree.excite(i,:) = strengths_und(exciteHit);
        end
        
        if do_controllability == 1
            %blah
        end
        
        if do_WD == 1
            WD.dti(i,:,:) = distance_wei(dti);
        end
        
        if do_PC == 1
            %% participation coefficient
            PC.baseline(i,:) = participation_coef(baseHit, Ci);
            PC.excite(i,:) = participation_coef(exciteHit, Ci);
            PC.inhib(i,:) = participation_coef(inhibHit, Ci);
        end
        
        if do_WMD == 1
            %% within-module degree (WMD), using the best baseline modularity decomposition
            WMD.baseline(i,:) = module_degree_zscore(baseHit, Ci);
            WMD.excite(i,:) = module_degree_zscore(exciteHit, Ci);
            WMD.inhib(i,:) = module_degree_zscore(inhibHit, Ci);
        end
        
        if do_BMD_Module_specific == 1
            %% Extra-module degree
            BMD_S.base = NaN(max(Ci),max(Ci),471);
            BMD_S.inhib = NaN(max(Ci),max(Ci),471);
            BMD_S.excite = NaN(max(Ci),max(Ci),471);
            for xx = 1:max(Ci)
                for zz = 1:max(Ci)
                    Koi.base=nansum(baseHit(Ci==xx,Ci==zz),2);
                    BMD_S.base(xx, zz,Ci==xx)=(Koi.base-mean(Koi.base))./std(Koi.base);
                    Koi.inhib=nansum(inhibHit(Ci==xx,Ci==zz),2);
                    BMD_S.inhib(xx, zz,Ci==xx)=(Koi.inhib-mean(Koi.inhib))./std(Koi.inhib);
                    Koi.excite=nansum(exciteHit(Ci==xx,Ci==zz),2);
                    BMD_S.excite(xx, zz,Ci==xx)=(Koi.excite-mean(Koi.excite))./std(Koi.excite);
                end
            end
            
            % this gives a subject x ROI x ROI x module matrix
            BMD_ROI.base(i,:,:,:) = BMD_S.base; % baseline
            BMD_MOD_S.base = nanmean(BMD_S.base, 3);
            BMD_MOD.base(i,:,:) = BMD_MOD_S.base;
            BMD_ROI.inhib(i,:,:,:) = BMD_S.inhib; % inhib
            BMD_MOD_S.inhib = nanmean(BMD_S.inhib, 3);
            BMD_MOD.inhib(i,:,:) = BMD_MOD_S.inhib;
            BMD_ROI.excite(i,:,:,:) = BMD_S.excite; % excite
            BMD_MOD_S.excite = nanmean(BMD_S.excite, 3);
            BMD_MOD.excite(i,:,:) = BMD_MOD_S.excite;
        end
        
        if do_BMD_not_specific == 1
            %% BMD, collapsing all extra-module connections
            BMD_S2.baseline = NaN(max(Ci),471);
            BMD_S2.inhib = NaN(max(Ci),471);
            BMD_S2.excite = NaN(max(Ci),471);
            for xx = 1:max(Ci)
                % hits
                Koi2.base=nansum(baseHit(Ci==xx,Ci~=xx),2);
                BMD_S2.baseline(xx, Ci==xx)=(Koi2.base-mean(Koi2.base))./std(Koi2.base);
                Koi2.inhib=nansum(inhibHit(Ci==xx,Ci~=xx),2);
                BMD_S2.inhib(xx, Ci==xx)=(Koi2.inhib-mean(Koi2.inhib))./std(Koi2.inhib);
                Koi2.excite=nansum(exciteHit(Ci==xx,Ci~=xx),2);
                BMD_S2.excite(xx, Ci==xx)=(Koi2.excite-mean(Koi2.excite))./std(Koi2.excite);
            end
            
            % this gives a subject x ROI x module matrix
            BMD_ROI2.baseline(i,:) = nanmean(BMD_S2.baseline, 1)'; % baseline
            BMD_ROI2.inhib(i,:) = nanmean(BMD_S2.inhib, 1)'; % inhib
            BMD_ROI2.excite(i,:) = nanmean(BMD_S2.excite, 1)'; % excite
        end
    end
end

% save modularity_values.mat SNSC1 SNSC2 SNSC3 WMD modularity degree
clear x y i

if do_stats == 1
    %% stats - module level
    %     %% ExMod Module Values
    %     for xx = 1:max(Ci)
    %         for yy = 1:max(Ci)
    %             [~,~,~,stats] = ttest(BMD_MOD.base(:,xx,yy));
    %             BMD_MODstats.base(xx,yy) = stats.tstat;
    %             [~,~,~,stats] = ttest(BMD_MOD.inhib(:,xx,yy));
    %             BMD_MODstats.inhib(xx,yy) = stats.tstat;
    %             [~,~,~,stats] = ttest(BMD_MOD.excite(:,xx,yy));
    %             BMD_MODstats.excite(xx,yy) = stats.tstat;
    %             [~,~,~,stats] = ttest(BMD_MOD.inhib(:,xx,yy), BMD_MOD.base(:,xx,yy));
    %             BMD_MODstats.HvM(xx,yy) = stats.tstat;
    %             [~,~,~,stats] = ttest(BMD_MOD.excite(:,xx,yy), BMD_MOD.base(:,xx,yy));
    %             BMD_MODstats.EvB(xx,yy) = stats.tstat;
    %             [~,~,~,stats] = ttest(BMD_MOD.excite(:,xx,yy), BMD_MOD.inhib(:,xx,yy));
    %             BMD_MODstats.EvI(xx,yy) = stats.tstat;
    %         end
    %     end
    
    %% stats - node level
    %     load modularity_values.mat
    for i = 1:471
        %%%% ROI-specific values
        %%%% degree
        [~,~,~,stats] = ttest(degree.excite(:,i));
        degreestats.excite(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(degree.inhib(:,i));
        degreestats.inhib(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(degree.baseline(:,i));
        degreestats.baseline(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(degree.excite(:,i), degree.baseline(:,i));
        degreestats.exciteVbaseline(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(degree.inhib(:,i), degree.baseline(:,i));
        degreestats.inhibVbaseline(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(degree.excite(:,i), degree.inhib(:,i));
        degreestats.exciteVinhib(i,1) = stats.tstat;
        
        %%%% PC
        [~,~,~,stats] = ttest(PC.baseline(:,i));
        PCstats.baseline(i) = stats.tstat;
        [~,~,~,stats] = ttest(PC.excite(:,i));
        PCstats.excite(i) = stats.tstat;
        [~,~,~,stats] = ttest(PC.inhib(:,i));
        PCstats.inhib(i) = stats.tstat;
        [~,~,~,stats] = ttest(PC.excite(:,i), PC.baseline(:,i));
        PCstats.exciteVbaseline(i) = stats.tstat;
        [~,~,~,stats] = ttest(PC.inhib(:,i), PC.baseline(:,i));
        PCstats.inhibVbaseline(i) = stats.tstat;
        [~,~,~,stats] = ttest(PC.excite(:,i), PC.inhib(:,i));
        PCstats.exciteVinhib(i) = stats.tstat;
        
        %%% WMD
        [~,~,~,stats] = ttest(WMD.baseline(:,i));
        WMDstats.baseline(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(WMD.excite(:,i));
        WMDstats.excite(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(WMD.inhib(:,i));
        WMDstats.inhib(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(WMD.excite(:,i), WMD.baseline(:,i));
        WMDstats.exciteVbaseline(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(WMD.inhib(:,i), WMD.baseline(:,i));
        WMDstats.inhibVbaseline(i,1) = stats.tstat;
        [~,~,~,stats] = ttest(WMD.excite(:,i), WMD.inhib(:,i));
        WMDstats.exciteVinhib(i,1) = stats.tstat;
        
        %%% BMD ROI values - split by module
        for xx = 1:max(Ci)
            for yy = 1:max(Ci)
                if isnan(BMD_ROI.base(1,xx,yy,i))
                    BMD_ROIstats.base(xx,yy,i) = NaN;
                else
                    [~,~,~,stats] = ttest(BMD_ROI.base(:,xx,yy,i));
                    BMD_ROIstats.base(xx,yy,i) = stats.tstat;
                    [~,~,~,stats] = ttest(BMD_ROI.inhib(:,xx,yy,i));
                    BMD_ROIstats.inhib(xx,yy,i) = stats.tstat;
                    [~,~,~,stats] = ttest(BMD_ROI.excite(:,xx,yy,i));
                    BMD_ROIstats.excite(xx,yy,i) = stats.tstat;
                    [~,~,~,stats] = ttest(BMD_ROI.inhib(:,xx,yy,i), BMD_ROI.base(:,xx,yy,i));
                    BMD_ROIstats.IvB(xx,yy,i) = stats.tstat;
                    [~,~,~,stats] = ttest(BMD_ROI.excite(:,xx,yy,i), BMD_ROI.base(:,xx,yy,i));
                    BMD_ROIstats.EvB(xx,yy,i) = stats.tstat;
                end
            end
        end
        
        %%% BMD ROI values - averaged over all extra-module connections
        [~,~,~,stats] = ttest(BMD_ROI2.baseline(:,i));
        BMD_ROI2stats.base(i,1) = stats.tstat';
        [~,~,~,stats] = ttest(BMD_ROI2.inhib(:,i));
        BMD_ROI2stats.inhib(i,1) = stats.tstat';
        [~,~,~,stats] = ttest(BMD_ROI2.excite(:,i));
        BMD_ROI2stats.excite(i,1) = stats.tstat';
        [~,~,~,stats] = ttest(BMD_ROI2.inhib(:,i), BMD_ROI2.baseline(:,i));
        BMD_ROI2stats.IvB(i,1) = stats.tstat';
        [~,~,~,stats] = ttest(BMD_ROI2.excite(:,i), BMD_ROI2.baseline(:,i));
        BMD_ROI2stats.EvB(i,1) = stats.tstat';
        [~,~,~,stats] = ttest(BMD_ROI2.excite(:,i), BMD_ROI2.inhib(:,i));
        BMD_ROI2stats.EvI(i,1) = stats.tstat';
        
        
    end
    %     save modularity_stats.mat WMDstats SNSCstats degreestats
    
end
% big_WMD = [reshape(WMD.baseline, 1080,1); reshape(WMD.baselineM, 1080,1); reshape(WMD.excite, 1080,1); reshape(WMD.exciteM, 1080,1); reshape(WMD.inhib, 1080,1); reshape(WMD.inhibM, 1080,1) ];
% big_mBMD = [reshape(BMD_ROI2.baseline, 1080,1); reshape(BMD_ROI2.baselineM, 1080,1); reshape(BMD_ROI2.excite, 1080,1); reshape(BMD_ROI2.exciteM, 1080,1); reshape(BMD_ROI2.inhib, 1080,1); reshape(BMD_ROI2.inhibM, 1080,1) ];

disp('fin')