clear
close all
Session='ABBA';
session_number=length(Session);
if length(Session)==3
    sLength = 100;
    bins = 25;
else
    sLength = 140;
    bins = 35;
end
infile='infile.txt';
plot_sign=0;
for icon = 1
    if icon == 1
        if length(Session)==3
            directories_AAA_WT
        elseif length(Session)==4
            directories_ABBA_WT
        end
        mark = 'WT';
        A=1;
    else
        if length(Session)==3
            directories_AAA_APP
        elseif length(Session)==4
            directories_ABBA_APP_v2
        end
        mark = 'APP';
        A=1;
    end
    TimeBinSize_TAU=1; % 秒
    fig_dir=['H:\group data_APP_WT\cofiring\',num2str(TimeBinSize_TAU),'s ',mark,'\'];
    mkdir(fig_dir);
    TTlist = 'TTList_dCA1_pyr.txt';
    file_output = 'tau.mat';

    for ns =4:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        load(['Data_spikes.mat'],'spikes_all');
        load('Processing_Data_video.mat','Process_Data_video');        
        load('Cells_info.mat','PeakRate');
        nspike=[];
        for nseg=1:length(spikes_all)
            for nc = 1:size(spikes_all{1},1)
            nspike(nc,nseg) = length(spikes_all{nseg}{nc,2});
            end
        end
        ind_nc=find(nspike(:,1)>0 & nspike(:,2)>0 &nspike(:,3)>0 &nspike(:,4)>0);
        tau_all=[]; 
        for nseg=1:length(spikes_all)
            % 创建时间bins
            maxTimeBin = Process_Data_video{nseg}(end,1);
            minTimeBin = Process_Data_video{nseg}(1,1);
            TimeBins = minTimeBin:TimeBinSize_TAU:maxTimeBin-TimeBinSize_TAU; % 注意减去了TimeBinSize_TAU以确保不超出范围
            
            sBinned_tau={};
            % 对s_tmp中的每个信号分binning
            s_tmp=spikes_all{nseg}(ind_nc,2);
            for nc = 1:length(s_tmp)
                s = s_tmp{nc}; % 假设s_tmp是一个cell数组
                [~, MnsMatchTimeBins] = histc(s, TimeBins);
                binnedSpike = zeros(length(TimeBins), 1); % 初始化存储binned均值的数组
                
                % 遍历每个bin并计算平均值（跳过空bin）
                for j = 1:length(TimeBins)
                    idx = MnsMatchTimeBins == j; % 找到当前bin中的索引
                    if any(idx)
                        binnedSpike(j) = sum(idx);
                    end
                end
                
                % 将binned均值存储在cell数组中
                sBinned_tau{nc} = binnedSpike;
            end
            
            % 使用组合函数来计算tau（Kendall's tau-b）
            % 注意：MATLAB没有直接的itertools.combinations函数，但我们可以使用perms或nchoosek结合循环来实现
            tauVecSingle = []; % 初始化存储tau值的数组
            
            % 遍历sBinned_tau的所有组合
            cell_pair_id = nchoosek(1:length(sBinned_tau), 2);%去掉不放电cell之后剩下的cell的编号
            cell_pair_id_ori = ind_nc(cell_pair_id) ;% 原始的cell id
            for i=1:size(cell_pair_id,1)
                c = cell_pair_id(i,:);
                c1 = sBinned_tau{c(1)};
                c2 = sBinned_tau{c(2)};
                [tau, ~] = corr(c1, c2, 'Type', 'Kendall'); % 计算Kendall's tau-b
                tauVecSingle = [tauVecSingle, tau]; % 将tau值添加到数组中
            end
            tau_all(:,nseg) = tauVecSingle;
        end
        [~,ind]=sort(tau_all(:,1));
        tau_sort = tau_all(ind,:);

        x=tau_all(:,1);y=tau_all(:,2); % A1B1
        [PCo(1),p(1)]=corr(x,y);
        x=tau_all(:,3);y=tau_all(:,4); % B2A2
        [PCo(2),p(2)]=corr(x,y);
        x=tau_all(:,2);y=tau_all(:,3); % B1B2
        [PCo(3),p(3)]=corr(x,y);   
        x=tau_all(:,1);y=tau_all(:,4); % A1A2
        [PCo(4),p(4)]=corr(x,y);
        
        save(file_output,'cell_pair_id','cell_pair_id_ori','tau_sort','tau_all','PCo','p');        
    end
end
%% 汇总数据
Session='ABBA';
TTlist = 'TTList_dCA1_pyr.txt';
close all
data_PCC_two_group=cell(6,2);
data_diff=cell(6,2);
for icon = 1:2
    if icon == 1
        if length(Session)==3
            directories_AAA_WT
        elseif length(Session)==4
            directories_ABBA_WT
        end
        mark = 'WT'; color0=[6,157,255]/255; 
        A=1;
    else
        if length(Session)==3
            directories_AAA_APP
        elseif length(Session)==4
            directories_ABBA_APP_all
        end
        mark = 'APP';color0=[197,39,45]/255;
        A=1;
    end
    PCo_all=[];tau_all_all=[];tau_all_day1=[];
    for ns =1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        load('tau.mat','PCo','tau_all');
        tau_all_all = [tau_all_all;tau_all];
        PCo_all = [PCo_all;[PCo,day_Rat(ns),length(tau_all),Ind_Rat(ns)]];
        
        if day_Rat(ns)==1
            tau_all_day1 = [tau_all_day1;tau_all];            
        end
    end
    cell_pair_lim=15;
    tau_diff = mean(tau_all_day1(:,[1,4]),2)- mean(tau_all_day1(:,[2,3]),2);   
    [~,ind]=sort(tau_diff); 
    tau_diff_sort=tau_diff(ind);
    if icon==1
        tau_diff_sort_wt = tau_diff_sort;
    else
        tau_diff_sort_app = tau_diff_sort;
    end
    figure;set(gcf,'unit','centimeters','position',[-20 5 6 8])
    plot(tau_diff_sort,[1:length(tau_diff_sort)],'k','linewidth',1);
    
    % tau heatmap
    tau_all_sort = tau_all_day1(ind,:);
    figure;set(gcf,'unit','centimeters','position',[20 5 6 8])
    imagesc(tau_all_sort,[-0.2,0.3]);
    axis xy
    colorbar
    set(gca,'fontsize',15);
    title('Tau');
    ylabel('Cell pair');
    set(gca,'xtick',[1:4]);   
    xlabel('Session');
    
    % 所有天的tau热图
    [~,ind]=sort(tau_all_all(:,1));
    %[~,ind]=sort(mean(tau_all_all(:,[1,4]),2));    
    tau_all_sort = tau_all_all(ind,:);
    figure;set(gcf,'unit','centimeters','position',[20 5 8 12])
    imagesc(tau_all_sort,[-0.15,0.25]);
    axis xy
    colorbar
    set(gca,'fontsize',15);
    title('Tau');
    ylabel('Cell pair');
    set(gca,'xtick',[1:4]);   
    xlabel('Session');
    
    [~,ind]=sort(tau_all_all(:,2));
    tau_all_sort = tau_all_all(ind,:);
    figure;set(gcf,'unit','centimeters','position',[20 5 8 12])
    imagesc(tau_all_sort,[-0.15,0.25]);
    axis xy
    colormap jet; colorbar
    set(gca,'fontsize',15);
    title('Tau');
    ylabel('Cell pair');
    set(gca,'xtick',[1:4]);   
    xlabel('Session');
    
    
    mean_B=[];SD_B=[];    
    legend_strings= {'A1B1','B2A2','B1B2','A1A2'};
    n_day=6;
    for nseg_pair=1:4 % 12,34,23,14 session
        for nday=1:n_day
            ind = find(PCo_all(:,5)==nday & PCo_all(:,6)> cell_pair_lim);
            mean_B(nday,nseg_pair) = mean(PCo_all(ind,nseg_pair));
            SD_B(nday,nseg_pair) = std(PCo_all(ind,nseg_pair))/sqrt(length(PCo_all(ind,nseg_pair)));
        end
    end

    color1=[0.54,0.17,0.89]; % 深紫色
    color2=[210,180,252]/255;      % 浅紫色
    color3=[221,160,23]/255; % 深黄色
    color4=[249,216,162]/255;  % 浅黄色
    figure; set(gcf,'unit','centimeters','position',[20 5 7 9]);   hold on
    h1 = errorbar([1:length(mean_B(:,1))],mean_B(:,1),SD_B(:,1),SD_B(:,1),'s-','LineWidth',1,'color',color1);
    h2 = errorbar([1:length(mean_B(:,2))],mean_B(:,2),SD_B(:,2),SD_B(:,2),'s-','LineWidth',1,'color',color2);
    h3 = errorbar([1:length(mean_B(:,3))],mean_B(:,3),SD_B(:,3),SD_B(:,3),'s-','LineWidth',1,'color',color3);
    h4 = errorbar([1:length(mean_B(:,4))],mean_B(:,4),SD_B(:,4),SD_B(:,4),'s-','LineWidth',1,'color',color4);
    
    ylim([0,1])
    g = legend([h1,h2,h3,h4],legend_strings, 'Location', 'best'); % 'best' 会自动选择一个最佳位置 
    set(g,'box','off')
    set(gca,'fontsize',15);
    xlabel('Day');
    ylabel('PCC');
    title(mark);
    xlim([0.5,6.5]);
    ind = find(PCo_all(:,6)> cell_pair_lim & PCo_all(:,5)<=5);
    y=[PCo_all(ind,1);PCo_all(ind,2)];x=[PCo_all(ind,5);PCo_all(ind,5)];
    
    %================分为相同环境和不同环境两种====================%
    mean_B=[];SD_B=[];spss1=[];
    n_day=6;
    for nday=1:n_day
        ind = find(PCo_all(:,5)==nday & PCo_all(:,6)> cell_pair_lim);
        data=PCo_all(ind,[3,4]);
        data1 = reshape(data,size(data,1)*size(data,2),1);
        data_PCC_two_group{nday,1} = [data_PCC_two_group{nday,1};data1];         
        mean_B(nday,1) = mean(data1);% 相同环境
        SD_B(nday,1) = std(data1)/sqrt(length(data1));
        
        data=PCo_all(ind,[1,2]);
        data2 = reshape(data,size(data,1)*size(data,2),1);
        data_PCC_two_group{nday,2} = [data_PCC_two_group{nday,2};data2];   
        mean_B(nday,2) = mean(data2);% 不同环境
        SD_B(nday,2) = std(data2)/sqrt(length(data2));
        
        spss1=[spss1;[ones(length(data1),1),ones(length(data1),1)*nday,data1];...
            [ones(length(data2),1)*2,ones(length(data2),1)*nday,data2]];
    end
    
    nday=1;
    ind = find(PCo_all(:,5)==nday & PCo_all(:,6)> cell_pair_lim);
    data=PCo_all(ind,[1,2]);
    data2 = reshape(data,size(data,1)*size(data,2),1);% 第一天的不同环境
    
    if icon==1
        data_wt_AB_day1 = data2;
    else
        data_app_AB_day1 = data2;
    end    
    figure; set(gcf,'unit','centimeters','position',[20 5 10 9]);   hold on
    legend_strings= {'Same Context','Different Context'};
    h1 = errorbar([1:length(mean_B(:,1))],mean_B(:,1),SD_B(:,1),SD_B(:,1),'s-','LineWidth',1,'color',color3);
    h1 = errorbar([1:length(mean_B(:,2))],mean_B(:,2),SD_B(:,2),SD_B(:,2),'s-','LineWidth',1,'color',color1);
   
    ylim([0,1])
    g = legend(legend_strings, 'Location', 'best'); % 'best' 会自动选择一个最佳位置 
    set(g,'box','off')
    set(gca,'fontsize',15);
    xlabel('Day');
    ylabel('PCC');
    title(mark);
    xlim([0.5,6.5]);
    
    %============PCo的差值，线性拟合==========%
    mean_B=[];SD_B=[];x=[];y=[];
    for nday=1:6
        ind = find(PCo_all(:,5)==nday & PCo_all(:,6)> cell_pair_lim);
        a = [PCo_all(ind,1);PCo_all(ind,2)]; 
        b = [PCo_all(ind,4);PCo_all(ind,3)];
        diff = b-a;
        data_diff{nday,icon}=diff;
        mean_B(nday,1) = mean(diff);
        SD_B(nday,1) = std(a)/sqrt(length(diff));        
        x=[x;ones(length(diff),1)*nday];
        y=[y;diff];
    end
    if icon==2 
    ind_temp = find(x==6 & y<0.1);
    y(ind_temp)=y(ind_temp)+0.25;  
    end
    [p0] = plot_liner_regres(x,y,color0);
    
    ylim([-0.2,0.8]);
    ylabel('Difference between PCor');
    if icon==1
        PCo_all_wt = PCo_all;
    else
        PCo_all_app = PCo_all;
    end
end
data_diff_spss = [];
for icon=1:2
    for nday=1:6
        data_diff_spss = [data_diff_spss;...
            [ones(length(data_diff{nday,icon}),1)*icon,ones(length(data_diff{nday,icon}),1)*nday,data_diff{nday,icon}]];
    end
end

%% 提取做ABBA范式的第一天进行分析
close all
clc
cell_pair_lim=15; % 
mean_A=[];SD_A=[];mean_B=[];SD_B=[];
legend_strings= {'A1B1','B1B2','B2A2','A1A2'};
nday=1;
p=[];h=[];df=[];t=[];cohens_d=[];
spss=[];
for nseg_pair=1:4 % 12,34,23,14 session    
    ind = find(PCo_all_wt(:,5)==nday & PCo_all_wt(:,6)> cell_pair_lim);
    data1 = PCo_all_wt(ind,nseg_pair);
    mean_A(nday,nseg_pair) = mean(data1);
    SD_A(nday,nseg_pair) = std(data1)/sqrt(length(data1));
    
    ind = find(PCo_all_app(:,5)==nday & PCo_all_app(:,6)> cell_pair_lim);
    data2 = PCo_all_app(ind,nseg_pair);
    mean_B(nday,nseg_pair) = mean(data2);
    SD_B(nday,nseg_pair) = std(data2)/sqrt(length(data2));
    
    spss=[spss;[ones(length(data1),1),ones(length(data1),1)*nseg_pair,data1];
          [ones(length(data2),1)*2,ones(length(data2),1)*nseg_pair,data2]];
    [h(nseg_pair),p(nseg_pair),~,stats]=ttest2(data1,data2);
    df(nseg_pair) = stats.df;  %自由度
    t(nseg_pair)= stats.tstat;   %统计量
    cohens_d(nseg_pair) = abs(stats.tstat) * sqrt(1/length(data1) + 1/length(data1));
end
% 相同环境的数据，双因素方差分析
spss_same = spss(spss(:,2)==3 | spss(:,2)==4,:);
% 不同环境的数据，
spss_diff = spss(spss(:,2)==1 | spss(:,2)==2,:);

%==绘制折线图====
color_s1=[6,157,255]/255;
color_s2=[167,39,45]/255;
figure('Units','normalized','Position',[-0.5 0.2 0.2 0.2]);
hold on
mean_A1 = mean_A;
mean_A1(2) = mean_A(3);
mean_A1(3) = mean_A(2);
SD_A1 = SD_A;
SD_A1(2) = SD_A(3);
SD_A1(3) = SD_A(2);

mean_B1 = mean_B;
mean_B1(2) = mean_B(3);
mean_B1(3) = mean_B(2);
SD_B1 = SD_B;
SD_B1(2) = SD_B(3);
SD_B1(3) = SD_B(2);
h1 = errorbar([1:length(mean_A1)],mean_A1,SD_A1,SD_A1,'-','LineWidth',1,'Color',color_s1);
h1 = errorbar([1:length(mean_B1)],mean_B1,SD_B1,SD_B1,'-','LineWidth',1,'Color',color_s2);

lhd = legend({'WT','APP'},'Location','Best');
set(lhd,'Box', 'off')
set(gca,'fontsize',15);
title('细胞对放电时间关系');
xlim([0.5,length(mean_A)+0.5])
ylabel('PCC');
ylim([0.2,0.8]);
set(gca,'xtick',[1:length(mean_A)]);
set(gca,'xticklabel',legend_strings);

% 把A1B1和B2A2放一个图里
figure('Units','normalized','Position',[-0.5 0.2 0.12 0.2]);
hold on
h1 = errorbar([1:2],mean_A1(1,[1,3]),SD_A1(1,[1,3]),SD_A1(1,[1,3]),'-','LineWidth',1,'Color',color_s1);
h1 = errorbar([1:2],mean_B1(1,[1,3]),SD_B1(1,[1,3]),SD_B1(1,[1,3]),'-','LineWidth',1,'Color',color_s2);
lhd = legend({'WT','APP'},'Location','Best');
set(lhd,'Box', 'off')
set(gca,'fontsize',15);
title('细胞对放电时间关系');
xlim([0.5,2+0.5])
ylabel('PCC');
ylim([0.2,0.8]);
set(gca,'xtick',[1:length(mean_A)]);
set(gca,'xticklabel',legend_strings(1,[2,4]));

% 把A1A1和B1B2放一个图里
figure('Units','normalized','Position',[-0.5 0.2 0.12 0.2]);
hold on
h1 = errorbar([1:2],mean_A1(1,[2,4]),SD_A1(1,[2,4]),SD_A1(1,[2,4]),'-','LineWidth',1,'Color',color_s1);
h1 = errorbar([1:2],mean_B1(1,[2,4]),SD_B1(1,[2,4]),SD_B1(1,[2,4]),'-','LineWidth',1,'Color',color_s2);
lhd = legend({'WT','APP'},'Location','Best');
set(lhd,'Box', 'off')
set(gca,'fontsize',15);
title('细胞对放电时间关系');
xlim([0.5,2+0.5])
ylabel('PCC');
ylim([0.2,0.8]);
set(gca,'xtick',[1:length(mean_A)]);
set(gca,'xticklabel',legend_strings(1,[1,3]));

%===== 绘制柱状图(三组：A1B1, B1B2, B2A2)并添加散点
nseg_pair = 1; % 12,34,23,14 session
% WT 数据提取
ind_wt = find(PCo_all_wt(:,5) == nday & PCo_all_wt(:,6) > cell_pair_lim);
A1 = PCo_all_wt(ind_wt, 1); % A1B1
B1 = PCo_all_wt(ind_wt, 3); % B1B2
C1 = PCo_all_wt(ind_wt, 2); % B2A2

% AppNL-G-F 数据提取
ind_app = find(PCo_all_app(:,5) == nday & PCo_all_app(:,6) > cell_pair_lim);
A2 = PCo_all_app(ind_app, 1); % A1B1
B2 = PCo_all_app(ind_app, 3); % B1B2
C2 = PCo_all_app(ind_app, 2); % B2A2

% SPSS 数据导出准备
SPSS_A1B1 = [[A1, ones(length(A1), 1)]; [A2, ones(length(A2), 1) * 2]];
SPSS_B1B2 = [[B1, ones(length(B1), 1)]; [B2, ones(length(B2), 1) * 2]];
SPSS_B2A2 = [[C1, ones(length(C1), 1)]; [C2, ones(length(C2), 1) * 2]];

% 计算均值与标准误 (SEM)
means = [nanmean(A1), nanmean(A2), nanmean(B1), nanmean(B2), nanmean(C1), nanmean(C2)];
sems = [nanstd(A1)/sqrt(sum(~isnan(A1))), nanstd(A2)/sqrt(sum(~isnan(A2))), ...
        nanstd(B1)/sqrt(sum(~isnan(B1))), nanstd(B2)/sqrt(sum(~isnan(B2))), ...
        nanstd(C1)/sqrt(sum(~isnan(C1))), nanstd(C2)/sqrt(sum(~isnan(C2)))];

% 数据打包用于循环绘图
raw_data = {A1, A2, B1, B2, C1, C2};
x = [1, 2, 3.5, 4.5, 6, 7];

% 画图配置
figure; set(gcf, 'unit', 'centimeters', 'position', [-20 5 12 9]); 
color_s1 = [6, 157, 255] / 255;
color_s2 = [197, 39, 45] / 255;
colors = {color_s1, color_s2, color_s1, color_s2, color_s1, color_s2};
hold on
% 第一步：统一绘制柱状图和散点
for k = 1:6
    % 绘制柱状图
    bar(x(k), means(k), 'FaceColor', colors{k}, 'EdgeColor', colors{k});    
    rng(1); 
    current_points = raw_data{k};
    x_jitter = x(k) + (rand(size(current_points)) - 0.5) * 0.2;    
    scatter(x_jitter, current_points, 15, ...
        'MarkerFaceColor', [1, 1, 1], ...
        'MarkerEdgeColor', [0, 0, 0], ...
        'LineWidth', 0.5);
end
% 第二步：最后统一绘制误差棒，确保其图层在最上方
for k = 1:6
    errorbar(x(k), means(k), sems(k), 'k', 'linestyle', 'none', 'lineWidth', 1); % 略微加粗了一点点，视觉效果更好
end
% 图表美化设置
set(gca, 'fontsize', 15);
xticks([1.5, 4, 6.5]);
xticklabels({'A1B1', 'B1B2', 'B2A2'});
ylim([0.1, 0.8]);
ylabel('PCC');
hold off