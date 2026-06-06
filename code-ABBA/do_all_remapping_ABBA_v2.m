%  ========= Remapping =========
clear
clc
Session='ABBA';
repetitions=1;
for icon = 1:2
    if icon == 1
        if length(Session)==3
            directories_AAA_WT
        elseif length(Session)==4
            directories_ABBA_WT_all
        end
        mark = 'WT';
    else
        if length(Session)==3
            directories_AAA_APP
        elseif length(Session)==4
            directories_ABBA_APP_all
        end
        mark = 'APP';
    end
    file_in = ['remapping_',Session,'.mat'];
    
    for ns = 1:isession 
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(file_in)
            fprintf(['=====day %.3g do not have stack cell data=====' '\n'], ns)
            continue
        end
        load(file_in)
        savefile_day=strcat('remapping_4measures_',Session,'_cellID'); % 在rate overlap内加入了筛选过后的cell id
         
        cutoff = 1; % cutoff by peak rate (peak rate >= 1Hz)
        nonzero = 0; % cutoff by firing rate (firing rate >= 0Hz)
        Stack = {stack_A1,stack_B2,stack_B3,stack_A4};
        RC = getactivebins(Stack,nonzero);
        Peak = [peakRate_A1,peakRate_B2,peakRate_B3,peakRate_A4];
        KC = get4session_keepcell_wxl(Peak,cutoff,nonzero); 
        RC = RC(KC);
        Stack_A1 = stack_A1(:,:,KC);
        Stack_B2 = stack_B2(:,:,KC);
        Stack_B3 = stack_B3(:,:,KC);
        Stack_A4 = stack_A4(:,:,KC);
        
        % A1-B2
        peak1=peakRate_A1;
        peak2=peakRate_B2;
        stack1=Stack_A1;
        stack2=Stack_B2;
        RateOverlap = RM_RateOverlap_cz_lzd(peak1,peak2,repetitions,cutoff,nonzero,KC);% 放电率重叠比；
        RateDI = RM_RateDI_cz(peak1,peak2,repetitions,cutoff,nonzero,KC);% 考察神经元峰值放电率在两个session中的变化趋势
        [PFCorr,stack1_active,stack2_active] = RM_SpatialCorr_cz_v3(stack1,stack2,repetitions,cutoff,nonzero,RC);
        PVCorr = RM_PVCorr_cz(stack1_active,stack2_active,repetitions,cutoff,nonzero);
        PVCrossCorr = RM_PVCrossCorr_cz(stack1_active,stack2_active,cutoff,nonzero,1); 
        RateOverlap_ABBA_A1B2=RateOverlap; % 放电率重叠比；
        RateDI_ABBA_A1B2=RateDI;% 神经元峰值放电率在两个session中的变化趋势
        PFCorr_ABBA_A1B2=PFCorr;% 同一个cell在不同session下的相关性
        PVCorr_ABBA_A1B2=PVCorr;% cell集群在不同session下的相关性
        PVCrossCorr_ABBA_A1B2=PVCrossCorr;
        
        % B2-B3
        peak1=peakRate_B2;
        peak2=peakRate_B3;
        stack1=Stack_B2;
        stack2=Stack_B3;
        RateOverlap = RM_RateOverlap_cz_lzd(peak1,peak2,repetitions,cutoff,nonzero,KC);
        RateDI = RM_RateDI_cz(peak1,peak2,repetitions,cutoff,nonzero,KC);
        [PFCorr,stack1_active,stack2_active] = RM_SpatialCorr_cz_v3(stack1,stack2,repetitions,cutoff,nonzero,RC);
        PVCorr = RM_PVCorr_cz(stack1_active,stack2_active,repetitions,cutoff,nonzero);
        PVCrossCorr = RM_PVCrossCorr_cz(stack1_active,stack2_active,cutoff,nonzero,1); 
        RateOverlap_ABBA_B2B3=RateOverlap;
        RateDI_ABBA_B2B3=RateDI;
        PFCorr_ABBA_B2B3=PFCorr;
        PVCorr_ABBA_B2B3=PVCorr;
        PVCrossCorr_ABBA_B2B3=PVCrossCorr;
        
        % B3-A4
        peak1=peakRate_B3;
        peak2=peakRate_A4;
        stack1=Stack_B3;
        stack2=Stack_A4;
        RateOverlap = RM_RateOverlap_cz_lzd(peak1,peak2,repetitions,cutoff,nonzero,KC);
        RateDI = RM_RateDI_cz(peak1,peak2,repetitions,cutoff,nonzero,KC);
        [PFCorr,stack1_active,stack2_active] = RM_SpatialCorr_cz_v3(stack1,stack2,repetitions,cutoff,nonzero,RC);
        PVCorr = RM_PVCorr_cz(stack1_active,stack2_active,repetitions,cutoff,nonzero);
        PVCrossCorr = RM_PVCrossCorr_cz(stack1_active,stack2_active,cutoff,nonzero,1); 
        RateOverlap_ABBA_B3A4=RateOverlap;
        RateDI_ABBA_B3A4=RateDI;
        PFCorr_ABBA_B3A4=PFCorr;
        PVCorr_ABBA_B3A4=PVCorr;
        PVCrossCorr_ABBA_B3A4=PVCrossCorr;
        
        % A1-A4
        peak1=peakRate_A1;
        peak2=peakRate_A4;
        stack1=Stack_A1;
        stack2=Stack_A4;
        RateOverlap = RM_RateOverlap_cz_lzd(peak1,peak2,repetitions,cutoff,nonzero,KC);
        RateDI = RM_RateDI_cz(peak1,peak2,repetitions,cutoff,nonzero,KC);
        [PFCorr,stack1_active,stack2_active] = RM_SpatialCorr_cz_v3(stack1,stack2,repetitions,cutoff,nonzero,RC);
        PVCorr = RM_PVCorr_cz(stack1_active,stack2_active,repetitions,cutoff,nonzero);
        PVCrossCorr = RM_PVCrossCorr_cz(stack1_active,stack2_active,cutoff,nonzero,1); 
        RateOverlap_ABBA_A1A4=RateOverlap;
        RateDI_ABBA_A1A4=RateDI;
        PFCorr_ABBA_A1A4=PFCorr;
        PVCorr_ABBA_A1A4=PVCorr;
        PVCrossCorr_ABBA_A1A4=PVCrossCorr;
        
        save(savefile_day,'RateOverlap_ABBA_A1B2','RateDI_ABBA_A1B2','PFCorr_ABBA_A1B2','PVCorr_ABBA_A1B2','PVCrossCorr_ABBA_A1B2',...
            'RateOverlap_ABBA_B2B3','RateDI_ABBA_B2B3','PFCorr_ABBA_B2B3','PVCorr_ABBA_B2B3','PVCrossCorr_ABBA_B2B3',...
            'RateOverlap_ABBA_B3A4','RateDI_ABBA_B3A4','PFCorr_ABBA_B3A4','PVCorr_ABBA_B3A4','PVCrossCorr_ABBA_B3A4',...
            'RateOverlap_ABBA_A1A4','RateDI_ABBA_A1A4','PFCorr_ABBA_A1A4','PVCorr_ABBA_A1A4','PVCrossCorr_ABBA_A1A4','KC');
        fprintf(['=====day %.3g done=====' '\n'], ns)
    end
end

%% stack all data
clear
clc
Session='ABBA';
for icon = 1:2
    if icon == 1
        if length(Session)==3
            directories_AAA_WT
        elseif length(Session)==4
            directories_ABBA_WT %
        end
        mark = 'WT';
    else
        if length(Session)==3
            directories_AAA_APP
        elseif length(Session)==4
            directories_ABBA_APP_v2 
        end
        mark = 'APP';
    end
    PVC=[];RateOverlap=cell(isession,4);PFCorr=cell(isession,4);Ncell=[];PFsize_PFCorr={};
    RateID=cell(isession,4);
    PFsize_PFCorr_all=[];PFCorr_all=[];KC_all=[];
    for ns = 1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        file_in = ['remapping_4measures_',Session,'_cellID.mat'];
        if ~exist(file_in,'file')
             PVC(ns,1:4)=nan;
            continue
        end
        load(file_in,'PVCorr_ABBA_A1B2','PVCorr_ABBA_B3A4','PVCorr_ABBA_B2B3','PVCorr_ABBA_A1A4',...
            'RateOverlap_ABBA_A1B2','RateOverlap_ABBA_B3A4','RateOverlap_ABBA_B2B3','RateOverlap_ABBA_A1A4',...
            'PFCorr_ABBA_A1B2','PFCorr_ABBA_B3A4','PFCorr_ABBA_B2B3','PFCorr_ABBA_A1A4','KC',...
            'RateDI_ABBA_A1B2')
        load('Cell_feature.mat','pfFeature');
        load(['remapping_',Session,'.mat']);
        Peak = [peakRate_A1,peakRate_B2,peakRate_B3,peakRate_A4];
        cutoff = 1; % cutoff by peak rate (peak rate >= 1Hz)
        nonzero = 0; % cutoff by firing rate (firing rate >= 0Hz)
        KC = get4session_keepcell_wxl(Peak,cutoff,nonzero); % 只要有一个session 放电率大于cutoff就留下
        
        KC_all=[KC_all;[KC,ones(length(KC),1)*ns]];
        for nseg=1:4
            temp = [pfFeature{nseg}.placefieldsize];
            PFsize_PFCorr{ns,1}(:,nseg) = temp(KC)';
        end
        PFsize_PFCorr_all = [PFsize_PFCorr_all;PFsize_PFCorr{ns,1}];
        
        PVC(ns,1)=nanmean( nanmean(PVCorr_ABBA_A1B2.PVCorr_matrix)); %A1B1 (实验顺序是：A1,B1,B2,A2)
        PVC(ns,2)=nanmean( nanmean(PVCorr_ABBA_B3A4.PVCorr_matrix)); %A2B2
        PVC(ns,3)=nanmean( nanmean(PVCorr_ABBA_B2B3.PVCorr_matrix)); %B1B2
        PVC(ns,4)=nanmean( nanmean(PVCorr_ABBA_A1A4.PVCorr_matrix)); %A1A2
        RateOverlap{ns,1} = [RateOverlap{ns,1},RateOverlap_ABBA_A1B2.RateRatio];
        RateOverlap{ns,2} = [RateOverlap{ns,2},RateOverlap_ABBA_B3A4.RateRatio];
        RateOverlap{ns,3} = [RateOverlap{ns,3},RateOverlap_ABBA_B2B3.RateRatio];
        RateOverlap{ns,4} = [RateOverlap{ns,4},RateOverlap_ABBA_A1A4.RateRatio];
        RateID{ns,1} = RateOverlap_ABBA_A1B2.Cell_id;
        RateID{ns,2} = RateOverlap_ABBA_B3A4.Cell_id;
        RateID{ns,3} = RateOverlap_ABBA_B2B3.Cell_id;
        RateID{ns,4} = RateOverlap_ABBA_A1A4.Cell_id;
        PFCorr{ns,1} = [PFCorr{ns,1},PFCorr_ABBA_A1B2.PFCorrActive];
        PFCorr{ns,2} = [PFCorr{ns,2},PFCorr_ABBA_B3A4.PFCorrActive];
        PFCorr{ns,3} = [PFCorr{ns,3},PFCorr_ABBA_B2B3.PFCorrActive];
        PFCorr{ns,4} = [PFCorr{ns,4},PFCorr_ABBA_A1A4.PFCorrActive];
        PFCorr0=[PFCorr_ABBA_A1B2.PFCorrActive',PFCorr_ABBA_B3A4.PFCorrActive',PFCorr_ABBA_B2B3.PFCorrActive',PFCorr_ABBA_A1A4.PFCorrActive'];
        PFCorr0=[PFCorr0,day_Rat(ns)*ones(size(PFCorr0,1),1)];
        PFCorr_all = [PFCorr_all;PFCorr0];        
        Ncell(ns,1) = length(RateOverlap_ABBA_A1B2.RateRatio);
    end
    ind_nan = find(Ncell<2); % cell小于2个的天要去掉PVC，
    PVC(ind_nan,:)=nan;
    if icon==1
        PVC_WT=PVC;
        RateOverlap_WT = RateOverlap;
        RateID_WT = RateID;
        PFCorr_WT = PFCorr;
        PFCorr_all_WT = PFCorr_all;
        day_Rat_WT = day_Rat;
        Ind_Rat_WT=Ind_Rat;
        Ncell_WT=Ncell;
        PFsize_PFCorr_all_WT = PFsize_PFCorr_all;
        PFsize_PFCorr_WT = PFsize_PFCorr;
        KC_WT=KC_all;
    else
        PVC_APP=PVC;
        RateOverlap_APP = RateOverlap;
        RateID_APP = RateID;
        PFCorr_APP = PFCorr;
        PFCorr_all_APP = PFCorr_all;
        day_Rat_APP =day_Rat;
        Ind_Rat_APP=Ind_Rat;
        Ncell_APP=Ncell;
        PFsize_PFCorr_all_APP = PFsize_PFCorr_all;
        PFsize_PFCorr_APP = PFsize_PFCorr;
        KC_APP=KC_all;
    end
end
save('H:\group data_APP_WT\remapping\ABBA_day_all.mat','PVC_WT','RateOverlap_WT','PFCorr_WT','day_Rat_WT','Ind_Rat_WT',...
    'PVC_APP','RateOverlap_APP','PFCorr_APP','day_Rat_APP','Ind_Rat_APP','Ncell_APP','Ncell_WT',...
    'PFsize_PFCorr_all_WT','PFsize_PFCorr_WT','PFsize_PFCorr_all_APP','PFsize_PFCorr_APP','PFCorr_all_WT','PFCorr_all_APP',...
    'KC_APP','KC_WT','RateID_APP','RateID_WT');

%===========================================================%
close all
n_day = 4;
title_ses={'A1-B1','A2-B2'};
for ses=1:2
PC1_stack=[];PC2_stack=[];
mean1=[];mean2=[];std1=[];std2=[];
for ngroup=1: n_day %
    ind1 = find(day_Rat_WT==ngroup);
    ind2 = find(day_Rat_APP==ngroup);
    
    PC1{ngroup} =  PVC_WT(ind1,ses);
    PC2{ngroup} =  PVC_APP(ind2,ses);
    PC1_stack = [PC1_stack;[PC1{ngroup},ones(length(PC1{ngroup}),1)*ngroup]];
    PC2_stack = [PC2_stack;[PC2{ngroup},ones(length(PC2{ngroup}),1)*ngroup]];
    mean1(ngroup,1) =  nanmean(PC1{ngroup});
    mean2(ngroup,1) =  nanmean(PC2{ngroup});
    std1(ngroup,1) = nanstd(PC1{ngroup})/sqrt(length(PC1{ngroup}));
    std2(ngroup,1) = nanstd(PC2{ngroup})/sqrt(length(PC2{ngroup}));
end
color_WT='b';
color_APP='r';
title12={'WT','APP'};
label0='Population vector correlation';
title0='Population vector correlation';
figure('Units','normalized','Position',[0.2 0.2 0.3 0.3]);
hold on
h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
set(gca,'xtick',[1:n_day]);
xlim([0.5,n_day+0.5])
lhd = legend([h1 h2],title12,'Location','best');
set(lhd,'Box', 'off')
title(title_ses{ses});
set(gca,'fontsize',15);
ylabel(label0);
ylim([0,1]);
end
%%
load('H:\group data_APP_WT\remapping\ABBA_day_all.mat')
%% PVC : bar chart
close all
n_day = 6;
title_ses={'A1-B1','A2-B2'};

PC1_stack=[];PC2_stack=[];
mean1=[];mean2=[];std1=[];std2=[];
PC1={};PC2={};
for ngroup=1: n_day %天数 % size(RateOverlap_ABBA_WT,1)
    ind1 = find(day_Rat_WT==ngroup);
    ind2 = find(day_Rat_APP==ngroup);
    
    PC1{ngroup} =  nanmean([PVC_WT(ind1,1),PVC_WT(ind1,2)],2); 
    PC2{ngroup} =  nanmean([PVC_APP(ind2,1),PVC_APP(ind2,2)],2);
    
    PC1{ngroup}(isnan(PC1{ngroup}))=[];PC2{ngroup}(isnan(PC2{ngroup}))=[];
    PC1_stack = [PC1_stack;[PC1{ngroup},ones(length(PC1{ngroup}),1)*ngroup]];
    PC2_stack = [PC2_stack;[PC2{ngroup},ones(length(PC2{ngroup}),1)*ngroup]];
    mean1(ngroup,1) =  nanmean(PC1{ngroup});
    mean2(ngroup,1) =  nanmean(PC2{ngroup});
    std1(ngroup,1) = nanstd(PC1{ngroup})/sqrt(length(PC1{ngroup}));
    std2(ngroup,1) = nanstd(PC2{ngroup})/sqrt(length(PC2{ngroup}));
end
color_WT='b';
color_APP='r';
title12={'WT','APP'};
label0='Population vector correlation';
title0='Population vector correlation';
figure('Units','normalized','Position',[0.2 0.2 0.3 0.3]);
hold on
h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
set(gca,'xtick',[1:n_day]);
% set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
xlim([0.5,n_day+0.5])
lhd = legend([h1 h2],title12,'Location','best');
set(lhd,'Box', 'off')
title('种群向量相关性');
set(gca,'fontsize',15);
ylabel(label0);
ylim([0,1]);


% 区分A1B1,A2B2
PC1_rat_WT={};PC1_rat_APP={};
mean_A=[];sd_A=[];mean_B=[];sd_B=[];SPSS=[];SPSS_wt=[];SPSS_app=[];
for day=1:n_day
    ind1 = find(day_Rat_WT==day);
    temp_A1B1_day = PVC_WT(ind1,1);
    temp_A2B2_day = PVC_WT(ind1,2);
    temp_A1B1_day(isnan(temp_A1B1_day))=[];
    temp_A2B2_day(isnan(temp_A2B2_day))=[];
    PC1_rat_WT{day,1} =  temp_A1B1_day;
    PC1_rat_WT{day,2} =  temp_A2B2_day;
    
    % 重复测量的方差分析
    SPSS_wt = [SPSS_wt;[ones(length(temp_A1B1_day),1),temp_A1B1_day,temp_A2B2_day]];
    
    mean_A(day,1) = nanmean(temp_A1B1_day);
    sd_A(day,1) = nanstd(temp_A1B1_day)/sqrt(length(temp_A1B1_day));
    mean_A(day,2) = nanmean(temp_A2B2_day);
    sd_A(day,2) = nanstd(temp_A2B2_day)/sqrt(length(temp_A2B2_day));
    %====================================================%    
    ind2 = find(day_Rat_APP==day);
    temp_A1B1_day = PVC_APP(ind1,1);
    temp_A2B2_day = PVC_APP(ind1,2);
    temp_A1B1_day(isnan(temp_A1B1_day))=[];
    temp_A2B2_day(isnan(temp_A2B2_day))=[];
    PC1_rat_APP{day,1} =  temp_A1B1_day;
    PC1_rat_APP{day,2} =  temp_A2B2_day;
    
    SPSS_app = [SPSS_app;[ones(length(temp_A1B1_day),1)*2,temp_A1B1_day,temp_A2B2_day]];
    mean_B(day,1) = nanmean(temp_A1B1_day);
    sd_B(day,1) = nanstd(temp_A1B1_day)/sqrt(length(temp_A1B1_day));
    mean_B(day,2) = nanmean(temp_A2B2_day);
    sd_B(day,2) = nanstd(temp_A2B2_day)/sqrt(length(temp_A2B2_day));
end
SPSS=[SPSS_wt;SPSS_app];
figure('Units','normalized','Position',[0.2 0.2 0.2 0.3]);
hold on
% WT-A1B1
h1 = errorbar([1:length(mean_A(:,1))],mean_A(:,1),sd_A(:,1),'s-','LineWidth',1,'Color',color_WT);
% WT-A2B2
h2 = errorbar([1:length(mean_A(:,2))],mean_A(:,2),sd_A(:,2),'s--','LineWidth',1,'Color',color_WT);
% APP-A1B1
h3 = errorbar([1:length(mean_B(:,1))],mean_B(:,1),sd_B(:,1),'s-','LineWidth',1,'Color',color_APP);
% APP-A2B2
h4 = errorbar([1:length(mean_B(:,2))],mean_B(:,2),sd_B(:,2),'s--','LineWidth',1,'Color',color_APP);

set(gca,'xtick',[1:n_day]);
xlim([0.5,n_day+0.5])
lhd = legend([h1 h2 h3 h4],{'WT-A1B1','WT-A2B2','APP-A1B1','APP-A2B2'},'Location','best');
set(lhd,'Box', 'off')
title('A-B');
set(gca,'fontsize',15);
label0 = 'PVC';
ylabel(label0);
ylim([0,1]);
title('种群向量相关性')
xlabel('Day');

%% Spatial correlation
PC1_rat_WT={};PC1_rat_APP={};
mean_A=[];sd_A=[];mean_B=[];sd_B=[];SPSS=[];
for day=1:n_day
    ind = find(day_Rat_WT==day);
    temp_day=[];
    for i = 1:length(ind)
%         temp = [PFCorr_WT{ind(i),1},PFCorr_WT{ind(i),2}]'; % A1B1 和 A2B2(
%         temp = [PFCorr_WT{ind(i),1}]'; 
        temp = nanmean([PFCorr_WT{ind(i),1};PFCorr_WT{ind(i),2}]',2);
        temp(isnan(temp))=[];
        temp_day = [temp_day;temp];
    end
    PC1_rat_WT{day} =  temp_day;
    SPSS = [SPSS;[ones(length(temp_day),1),ones(length(temp_day),1)*day,temp_day]];
    
    mean_A(day) = nanmean(temp_day);
    sd_A(day) = nanstd(temp_day)/sqrt(length(temp_day));
    %====================================================%
    ind = find(day_Rat_APP==day);
    temp_day=[];
    for i = 1:length(ind)
        temp = nanmean([PFCorr_APP{ind(i),1};PFCorr_APP{ind(i),2}]',2);        
        temp(isnan(temp))=[];
        temp_day = [temp_day;temp];
    end
    PC1_rat_APP{day} =  temp_day;
    SPSS = [SPSS;[ones(length(temp_day),1)*2,ones(length(temp_day),1)*day,temp_day]];
    mean_B(day) = nanmean(temp_day);
    sd_B(day) = nanstd(temp_day)/sqrt(length(temp_day));
end
figure('Units','normalized','Position',[0.2 0.2 0.3 0.3]);
hold on
h1 = errorbar([1:length(mean_A)],mean_A,sd_A,'s-','LineWidth',2,'Color',color_WT);
h2 = errorbar([1:length(mean_B)],mean_B,sd_B,'s-','LineWidth',2,'Color',color_APP);
set(gca,'xtick',[1:n_day]);
% set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
xlim([0.5,n_day+0.5])
lhd = legend([h1 h2],title12,'Location','best');
set(lhd,'Box', 'off')
title('A-B');
set(gca,'fontsize',15);
label0 = 'Spatial correlation';
ylabel(label0);
ylim([0,1]);
ylim([-0.1,0.6]);
title('空间相关性')
xlabel('Day');

% 区分A1B1,A2B2
PC1_rat_WT={};PC1_rat_APP={};
mean_A=[];sd_A=[];mean_B=[];sd_B=[];SPSS=[];SPSS_wt=[];
for day=1:n_day
    ind = find(day_Rat_WT==day);
    temp_A1B1_day=[]; temp_A2B2_day=[];
    for i = 1:length(ind)
        temp_A1B1 = PFCorr_WT{ind(i),1}';
        temp_A2B2 = PFCorr_WT{ind(i),2}';
        temp_A1B1_day = [temp_A1B1_day;temp_A1B1];
        temp_A2B2_day = [temp_A2B2_day;temp_A2B2];   
    end
    temp_A1B1_day(isnan(temp_A1B1_day))=[];
    temp_A2B2_day(isnan(temp_A2B2_day))=[];
    
    PC1_rat_WT{day,1} =  temp_A1B1_day;
    PC1_rat_WT{day,2} =  temp_A2B2_day;
    
    % 重复测量的方差分析
    SPSS_wt = [SPSS_wt;[ones(length(temp_A1B1_day),1),temp_A1B1_day,temp_A2B2_day]];
    
    mean_A(day,1) = nanmean(temp_A1B1_day);
    sd_A(day,1) = nanstd(temp_A1B1_day)/sqrt(length(temp_A1B1_day));
    mean_A(day,2) = nanmean(temp_A2B2_day);
    sd_A(day,2) = nanstd(temp_A2B2_day)/sqrt(length(temp_A2B2_day));
    %====================================================%
    ind = find(day_Rat_APP==day);
    temp_day=[];
    for i = 1:length(ind)
        temp_A1B1 = PFCorr_APP{ind(i),1}';
        temp_A2B2 = PFCorr_APP{ind(i),2}';
        temp_A1B1_day = [temp_A1B1_day;temp_A1B1];
        temp_A2B2_day = [temp_A2B2_day;temp_A2B2];   
    end
    temp_A1B1_day(isnan(temp_A1B1_day))=[];
    temp_A2B2_day(isnan(temp_A2B2_day))=[];
    PC1_rat_APP{day,1} =  temp_A1B1_day;
    PC1_rat_APP{day,2} =  temp_A2B2_day;
%     SPSS = [SPSS;[ones(length(temp_day),1)*2,ones(length(temp_day),1)*day,temp_day]];
    mean_B(day,1) = nanmean(temp_A1B1_day);
    sd_B(day,1) = nanstd(temp_A1B1_day)/sqrt(length(temp_A1B1_day));
    mean_B(day,2) = nanmean(temp_A2B2_day);
    sd_B(day,2) = nanstd(temp_A2B2_day)/sqrt(length(temp_A2B2_day));
end
figure('Units','normalized','Position',[0.2 0.2 0.2 0.3]);
hold on
% WT-A1B1
h1 = errorbar([1:length(mean_A(:,1))],mean_A(:,1),sd_A(:,1),'s-','LineWidth',1,'Color',color_WT);
% WT-A2B2
h2 = errorbar([1:length(mean_A(:,2))],mean_A(:,2),sd_A(:,2),'s--','LineWidth',1,'Color',color_WT);
% APP-A1B1
h3 = errorbar([1:length(mean_B(:,1))],mean_B(:,1),sd_B(:,1),'s-','LineWidth',1,'Color',color_APP);
% APP-A2B2
h4 = errorbar([1:length(mean_B(:,2))],mean_B(:,2),sd_B(:,2),'s--','LineWidth',1,'Color',color_APP);

set(gca,'xtick',[1:n_day]);
% set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
xlim([0.5,n_day+0.5])
lhd = legend([h1 h2 h3 h4],{'WT-A1B1','WT-A2B2','APP-A1B1','APP-A2B2'},'Location','best');
set(lhd,'Box', 'off')
title('A-B');
set(gca,'fontsize',15);
label0 = 'Spatial correlation';
ylabel(label0);
ylim([0,1]);
ylim([0.1,0.4]);
title('空间相关性')
xlabel('Day');

%% 跨天的rate overlap; 
close all
PC1_rat_WT={};PC1_rat_APP={};
mean_A=[];sd_A=[];mean_B=[];sd_B=[];SPSS=[];
for day=1:n_day
    ind = find(day_Rat_WT==day);
    temp_day=[];
    for i = 1:length(ind)
        temp = nanmean([RateOverlap_WT{ind(i),1};RateOverlap_WT{ind(i),2}],2);
        temp(isnan(temp))=[];
        temp_day = [temp_day;temp];
    end
    PC1_rat_WT{day} =  temp_day;
    SPSS = [SPSS;[ones(length(temp_day),1),ones(length(temp_day),1)*day,temp_day]];
    
    mean_A(day) = nanmean(temp_day);
    sd_A(day) = nanstd(temp_day)/sqrt(length(temp_day));
    %====================================================%
    ind = find(day_Rat_APP==day);
    temp_day=[];
    for i = 1:length(ind)              
        temp = nanmean([RateOverlap_APP{ind(i),1};RateOverlap_APP{ind(i),2}],2);        
        temp(isnan(temp))=[];
        temp_day = [temp_day;temp];
    end
    PC1_rat_APP{day} =  temp_day;
    SPSS = [SPSS;[ones(length(temp_day),1)*2,ones(length(temp_day),1)*day,temp_day]];
    mean_B(day) = nanmean(temp_day);
    sd_B(day) = nanstd(temp_day)/sqrt(length(temp_day));
end
figure('Units','normalized','Position',[0.2 0.2 0.3 0.3]);
hold on
h1 = errorbar([1:length(mean_A)],mean_A,sd_A,'s-','LineWidth',2,'Color',color_WT);
h2 = errorbar([1:length(mean_B)],mean_B,sd_B,'s-','LineWidth',2,'Color',color_APP);
set(gca,'xtick',[1:n_day]);
xlim([0.5,n_day+0.5])
lhd = legend([h1 h2],title12,'Location','best');
set(lhd,'Box', 'off')
title('A-B');
set(gca,'fontsize',15);
label0 = 'Rate overlap';
ylabel(label0);
ylim([0.3,0.8]);
title('放电率重叠比')
xlabel('Day');
