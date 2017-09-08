%clear all
% TMS Enhance network comparison
subjects = {'002_TP' '003_RL' '004_RS' '005_TWS' '006_GH' '007_VD' '008_JS' '009_CY' '010_MU' '011_LD' '012_SD' '013_RH' '014_JM'};

% Module ROIs (9 Modules, 8 is left PFC)
Module_ROIs = {1,3,6,9,35,41,42,64,65,66,93,94,112,135,150,168,181,188,196,200,209,215,217,236,239,242,247,252,293,302,306,314,335,367,372,374,392,399,403,404,415,425,431,433,443,447,451,464,466,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];18,34,46,49,60,63,67,70,71,76,78,88,102,103,107,117,120,121,129,142,144,148,182,185,194,214,225,227,229,243,263,266,277,283,288,294,308,316,323,330,332,338,343,355,356,371,411,422,436,446,461,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];4,8,11,19,23,30,32,39,43,59,61,69,72,79,81,97,99,110,115,122,128,133,145,152,155,156,157,162,165,167,169,171,193,218,220,222,238,255,257,264,274,276,282,286,298,322,324,339,340,376,377,384,394,407,420,421,439,453,457,471,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];12,13,16,28,37,45,62,82,83,85,87,101,108,109,113,116,125,126,130,136,140,141,153,158,163,183,192,195,201,208,211,241,249,259,260,270,291,296,299,337,350,362,383,385,386,387,388,397,398,412,452,460,465,468,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];2,14,20,22,89,114,139,147,170,172,179,180,184,198,205,206,226,237,240,256,290,313,318,325,336,344,349,357,366,368,369,381,382,390,395,396,405,406,408,409,413,416,417,430,437,441,444,459,462,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];5,7,15,50,54,57,73,75,84,98,119,124,134,138,149,151,175,178,186,187,189,204,213,216,219,228,245,246,254,273,278,287,304,309,310,311,315,321,345,346,351,360,361,373,400,410,419,424,434,450,456,463,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];17,27,29,33,36,38,40,44,47,68,74,77,96,100,104,106,118,123,132,160,176,177,234,258,261,272,280,295,300,307,317,319,327,342,347,358,363,370,380,389,393,401,402,414,423,426,432,458,469,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];10,24,25,31,51,52,53,56,58,86,90,92,105,146,159,161,164,174,190,202,203,207,210,212,221,224,230,231,232,235,244,248,250,251,253,265,268,269,271,275,279,281,284,285,289,292,301,303,305,312,320,326,328,329,333,334,348,352,353,359,364,375,378,379,418,427,428,429,435,440,442,445,449,455,467,470;21,26,48,55,80,91,95,111,127,131,137,143,154,166,173,191,197,199,223,233,262,267,297,331,341,354,365,391,438,448,454,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};

% Subsequent Memory ROIs - WMD 1Hz, WMD 5Hz, BMD 1Hz, BMD 5Hz
% SME_ROIs =
% {1,21,137,149,213,280,338,402,411,450,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[];1,21,30,58,72,95,104,115,122,153,154,189,200,215,218,220,222,231,242,267,287,325,326,360,365,371,381,409,421,427,448,459;13,78,95,98,163,189,190,213,218,251,252,348,361,432,446,449,456,463,[],[],[],[],[],[],[],[],[],[],[],[],[],[];21,27,41,42,113,118,125,135,171,176,180,181,182,202,205,206,228,259,288,299,322,360,371,394,432,[],[],[],[],[],[],[]}; % t > 2.1
SME_ROIs = {1,21,72,115,122,153,154,200,267,365,421,427;149,280,411,[],[],[],[],[],[],[],[],[];176,228,288,[],[],[],[],[],[],[],[],[];78,190,213,252,361,446,[],[],[],[],[],[]};% t > 2.5


adjustment = 2; % 1 for inverse, 2 for sigmoid decay function
regionwise = 0; % 1 for one region, 2 for SME ROIs, 3 for Module ROIs, 0 for whole-brain vector

for i = 1:13
    % read in files
    baselineHitfile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s__BASE_HIT.txt', subjects{i});
    exciteHitfile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s__5Hz_HIT.txt', subjects{i});
    inhibHitfile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s__1Hz_HIT.txt', subjects{i});
    baselineESAfile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s__BASE_ESA.txt', subjects{i});
    exciteESAfile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s__5Hz_ESA.txt', subjects{i});
    inhibESAfile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s__1Hz_ESA.txt', subjects{i});
    
    dtifile = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471/%s_HOAsp_connectome.csv', subjects{i});
    dtifile2 = sprintf('/Users/simonwdavis/Dropbox/duke/TMSEnhance/Connectomes471_DSI/%s_matrix_gfa.mat', subjects{i});
    dti = dlmread(dtifile);
    
    baseHit = dlmread(baselineHitfile); baseHit(isnan(baseHit)) = 1;
    exciteHit = dlmread(exciteHitfile);
    inhibHit = dlmread(inhibHitfile);
    baseESA = dlmread(baselineESAfile);
    exciteESA = dlmread(exciteESAfile);
    inhibESA = dlmread(inhibESAfile);
%     pre_dti = load(dtifile2);
%     dti2 = pre_dti.connectivity;
    
%     if size(dti1,1)==size(dti2,1)
%         penis = ~(eye(471));
%         dti1 = dti1.*penis; dti2=dti2.*penis;
%         nums(i,1)=nnz(dti1); nums(i,2)=nnz(dti2);
%         dti1(dti1>0)=1; dti2(dti2>0)=1;
%         dti1(dti1==0) = NaN; dti2(dti2==0) = NaN;
%         
%         [cock,dick] = corrcoef(reshape(dti1,471*471,1), reshape(dti2,471*471,1), 'rows', 'pairwise');
%         dicks(i) = dick(2,1);
%     else
%         fprintf('no go on subject %s\n', subjects{i});
%     end
%     
    if adjustment == 1
        % inverse adjustment
        dti_adj = 1./dti; dti_adj(dti_adj==Inf) = NaN;
        %     hist(reshape(dti_adj, 471*471,1), 100)
        
    elseif adjustment == 2
        % sigmoidal adjustment
        dti_scaled = mat2gray(dti); dti_scaled(dti_scaled==0) = NaN;
        [Y,I]=sort(reshape(dti_scaled,471*471,1));
        bval = round(0.95*sum(~isnan(Y))); b = Y(bval);
        dti_adj =  1 ./ (1+exp(12*(dti_scaled-b)));
        adjusted_dti(i,:,:) = dti_adj;
        %     hist(reshape(dti_adj, 471*471,1), 100)
    end
    
    % perform the weighted distance calculation on the adjusted DTI matrix
    [distance(i,:,:), numpaths(i,:,:)] = distance_wei(dti_adj);
    
    for jj = 1:25
        cock = squeeze(distance(i,:,:)); pussy = squeeze(numpaths(i,:,:));
        pussy(pussy>jj)=NaN; pussy(pussy>0)=1; pussy(pussy==0)=NaN;
        dildo = cock.*pussy;
        
        if regionwise == 1
            for kk = 58 % do just one ROI, treat as one vector
                % interesting regions: 190 391 58
                RvalB = corrcoef(dildo(kk,:), baseESA(kk,:), 'rows', 'pairwise');
                CrossCors.base(i,jj) = RvalB(2,1);
                RvalE = corrcoef(dildo(kk,:), exciteESA(kk,:), 'rows', 'pairwise');
                CrossCors.excite(i,jj) = RvalE(2,1);
                RvalI = corrcoef(dildo(kk,:), inhibESA(kk,:), 'rows', 'pairwise');
                CrossCors.inhib(i,jj) = RvalI(2,1);
            end
            clear cock pussy dildo
        elseif regionwise == 2  %%% SME ROIs, treat as one vector
            ROIs_5Hz = cell2mat(SME_ROIs(1,:));
            ROIs_1Hz = cell2mat(SME_ROIs(4,:));  
            structMatrix_1Hz = dildo(ROIs_1Hz,:); structMatrix_5Hz = dildo(ROIs_5Hz,:);  %if useing SME ROIs, make sure it matches the stimulation type
            exciteESA2=exciteESA(ROIs_5Hz,:);inhibESA2=inhibESA(ROIs_1Hz,:); % baseESA2=baseESA(ROIs,:);
            
%             RvalB = corrcoef(reshape(dildo2,size(ROIs,2)*471,1), reshape(baseESA2,size(ROIs,2)*471,1), 'rows', 'pairwise');
%             CrossCors.base(i,jj) = RvalB(2,1);
            RvalE = corrcoef(reshape(structMatrix_5Hz,size(ROIs_5Hz,2)*471,1), reshape(exciteESA2,size(ROIs_5Hz,2)*471,1), 'rows', 'pairwise');
            CrossCors.excite(i,jj) = RvalE(2,1);
            RvalI = corrcoef(reshape(structMatrix_1Hz,size(ROIs_1Hz,2)*471,1), reshape(inhibESA2,size(ROIs_1Hz,2)*471,1), 'rows', 'pairwise');
            CrossCors.inhib(i,jj) = RvalI(2,1);
        elseif regionwise == 3 %%% Module ROIs
            ROIs = cell2mat(Module_ROIs(4,:));
            structMatrix = dildo(ROIs,:); 
            baseESA2=baseESA(ROIs,:);exciteESA2=exciteESA(ROIs,:);inhibESA2=inhibESA(ROIs,:);
            RvalB = corrcoef(reshape(structMatrix,size(ROIs,2)*471,1), reshape(baseESA2,size(ROIs,2)*471,1), 'rows', 'pairwise');
            CrossCors.base(i,jj) = RvalB(2,1);
            RvalE = corrcoef(reshape(structMatrix,size(ROIs,2)*471,1), reshape(exciteESA2,size(ROIs,2)*471,1), 'rows', 'pairwise');
            CrossCors.excite(i,jj) = RvalE(2,1);
            RvalI = corrcoef(reshape(structMatrix,size(ROIs,2)*471,1), reshape(inhibESA2,size(ROIs,2)*471,1), 'rows', 'pairwise');
            CrossCors.inhib(i,jj) = RvalI(2,1);
        elseif regionwise == 0  % whole brain vector
            RvalB = corrcoef(reshape(dildo,471*471,1), reshape(baseESA,471*471,1), 'rows', 'pairwise');
            CrossCors.base(i,jj) = RvalB(2,1);
            RvalE = corrcoef(reshape(dildo,471*471,1), reshape(exciteESA,471*471,1), 'rows', 'pairwise');
            CrossCors.excite(i,jj) = RvalE(2,1);
            RvalI = corrcoef(reshape(dildo,471*471,1), reshape(inhibESA,471*471,1), 'rows', 'pairwise');
            CrossCors.inhib(i,jj) = RvalI(2,1);
            clear cock pussy dildo
        end
    end
    fprintf('done with subject %d\n', i);
end
