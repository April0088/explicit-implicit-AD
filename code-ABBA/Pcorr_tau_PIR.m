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
            directories_ABBA_APP_all
        end
        mark = 'APP';
        A=2;
    end
    TimeBinSize_TAU=1; % 秒
    fig_dir=['H:\group data_APP_WT\cofiring\',num2str(TimeBinSize_TAU),'s ',mark,'\'];
    mkdir(fig_dir);
    TTlist = 'TTList_dCA1_pyr.txt';
    file_output = 'PIR_tau.mat';

    for ns =[17,21,22] % A:isession
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
        load('position_independent_rate','r_obs_all','r_pti_all');
        
        % 去掉放电率完全为0的cell（因为放电率为0时，计算得到的相关性是nan）
        nspike=[]; 
        for nseg=1:length(spikes_all)
            for nc = 1:size(spikes_all{1},1)
                nspike(nc,nseg) = nansum(r_pti_all{nseg}(nc,:));
            end
        end
        ind_nc=find(nspike(:,1)~=0 & nspike(:,2)~=0 &nspike(:,3)~=0 &nspike(:,4)~=0);
        tau_all=[]; 
        for nseg=1:length(spikes_all)
    
            r_pti=r_pti_all{nseg}(ind_nc,:);
   
            % 使用组合函数来计算tau（Kendall's tau-b）
            % 注意：MATLAB没有直接的itertools.combinations函数，但我们可以使用perms或nchoosek结合循环来实现
            tauVecSingle = []; % 初始化存储tau值的数组
            
            % 遍历sBinned_tau的所有组合
            cell_pair_id = nchoosek(1:length(ind_nc), 2);%去掉不放电cell之后剩下的cell的编号
            cell_pair_id_ori = ind_nc(cell_pair_id) ;% 原始的cell id
            for i=1:size(cell_pair_id,1)
                c = cell_pair_id(i,:);
                c1 = r_pti(c(1),:)';
                c2 = r_pti(c(2),:)';
                c1(isnan(c1))=0;
                c2(isnan(c2))=0;
                [tau, ~] = corr(c1, c2, 'Type', 'Kendall'); % 计算Kendall's tau-b
                tauVecSingle = [tauVecSingle, tau]; % 将tau值添加到数组中
            end
            tau_all(:,nseg) = tauVecSingle;
        end
        [~,ind]=sort(tau_all(:,1));
        tau_sort = tau_all(ind,:);

        x=tau_all(:,1);y=tau_all(:,2); % 1A-1B
        [PCo(1),p(1)]=corr(x,y);
        x=tau_all(:,3);y=tau_all(:,4); % 2A-2B
        [PCo(2),p(2)]=corr(x,y);
        x=tau_all(:,2);y=tau_all(:,3); % 1B-2B
        [PCo(3),p(3)]=corr(x,y);   
        x=tau_all(:,1);y=tau_all(:,4); % 1A-2A
        [PCo(4),p(4)]=corr(x,y);
        
        save(file_output,'cell_pair_id','cell_pair_id_ori','tau_sort','tau_all','PCo','p');
        fprintf('第 %d 天计算结束...\n', ns);
    end
end

%% 汇总数据 并画统计图
Session='ABBA';
 TTlist = 'TTList_dCA1_pyr.txt';
close all
for icon = 1:2
    if icon == 1
        if length(Session)==3
            directories_AAA_WT
        elseif length(Session)==4
            directories_ABBA_WT
        end
        mark = 'WT'; color0=[6,157,255]/255;
    else
        if length(Session)==3
            directories_AAA_APP
        elseif length(Session)==4
            directories_ABBA_APP_all
        end
        mark = 'APP';color0=[197,39,45]/255;
    end
    PCo_all=[];tau_all_all=[];tau_all_day1=[];
    tau_all_day=cell(1,8);
    for ns =1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        load('PIR_tau.mat','PCo','tau_all');
        tau_all_all = [tau_all_all;tau_all];
        
        %顺序为12,34,23,14，第几天，几个cell pairs
        PCo_all = [PCo_all;[PCo,day_Rat(ns),length(tau_all)]];
        
        if day_Rat(ns)==1
            tau_all_day1 = [tau_all_day1;tau_all];            
        end
        
        nday = day_Rat(ns);
        tau_all_day{nday} = [tau_all_day{nday};tau_all];
    end
    if icon==1
        tau_all_day_wt=tau_all_day;
    else
        tau_all_day_app=tau_all_day;
    end
    cell_pair_lim=15;    
  
    %================分为相同环境和不同环境两种====================%
    mean_B=[];SD_B=[];
    figure; set(gcf,'unit','centimeters','position',[20 5 10 9]);   hold on
    legend_strings= {'Same Context','Different Context'};
    n_day=6;
    for nday=1:n_day
        ind = find(PCo_all(:,5)==nday & PCo_all(:,6)> cell_pair_lim);
        %data1 = mean(PCo_all(ind,[3,4]),2);
        data=PCo_all(ind,[3,4]);
        data1 = reshape(data,size(data,1)*size(data,2),1);
        mean_B(nday,1) = mean(data1);% 相同环境
        SD_B(nday,1) = std(data1)/sqrt(length(data1));
        
        data=PCo_all(ind,[1,2]);
        data2 = reshape(data,size(data,1)*size(data,2),1);
        data_AB_day{nday} = data2;
        mean_B(nday,2) = mean(data2);% 不同环境
        SD_B(nday,2) = std(data2)/sqrt(length(data2));
    end
    if icon==1
        PCo_all_wt = PCo_all;
    else
        PCo_all_app = PCo_all;
    end
    
    h1 = errorbar([1:length(mean_B(:,1))],mean_B(:,1),SD_B(:,1),SD_B(:,1),'s-','LineWidth',1.5,'color',color2);
    h1 = errorbar([1:length(mean_B(:,2))],mean_B(:,2),SD_B(:,2),SD_B(:,2),'s-','LineWidth',1.5,'color',color1);
   
    ylim([0,1])
    g = legend(legend_strings, 'Location', 'best'); % 'best' 会自动选择一个最佳位置 
    set(g,'box','off')
    set(gca,'fontsize',15);
    xlabel('Day');
    ylabel('PCo(r)');
    title(mark);
    xlim([0.5,6.5]);
end
%%
file_name=['H:\group data_APP_WT\PCC_PIR.mat'];
save(file_name,'PCo_all_wt','PCo_all_app');

%% plot bar chart
load('H:\group data_APP_WT\PCC_PIR.mat');
close all
color1=[6,157,255]/255;
color2=[197,39,45]/255;

% 观察两组大鼠的B2A2的结果
figure('Units','normalized','Position',[-0.4 0.2 0.08 0.25]);
nbar=0;
x=[1,2];
ind = find(PCo_all_wt(:,5)<=6 & PCo_all_wt(:,6)> cell_pair_lim);
A = PCo_all_wt(ind,2);
ind = find(PCo_all_app(:,5)<=6 & PCo_all_app(:,6)> cell_pair_lim);
B = PCo_all_app(ind,2);
[H,P,~,stats]=ttest2(A,B)
% 计算 Cohen's d 并直接打印所有结果
d = (nanmean(A) - nanmean(B)) / sqrt(((sum(~isnan(A))-1)*nanstd(A)^2 + (sum(~isnan(B))-1)*nanstd(B)^2) / stats.df);
fprintf('df = %d, t = %.3f, p = %.4e, Cohen''s d = %.3f\n', stats.df, stats.tstat, P, d);
spss_data=[[ones(length(A),1),A];[ones(length(B),1)*2,B]];

mean_A =nanmean(A);
mean_B = nanmean(B);
sd_A = nanstd(A)/sqrt(length(A));
sd_B = nanstd(B)/sqrt(length(B));
hold on
nbar=nbar+1;
bar(x(nbar),mean_A,'FaceColor',color1,'EdgeColor', color1);
Err1 = errorbar(x(nbar),mean_A,sd_A,'k', 'linestyle', 'none', 'lineWidth', 1);
nbar=nbar+1;
bar(x(nbar),mean_B,'FaceColor',color2,'EdgeColor', color2);
Err2 = errorbar(x(nbar),mean_B,sd_B,'k', 'linestyle', 'none', 'lineWidth', 1);
set(gca,'fontsize',15);
ylim([0,0.8]);
title('B2A2');

 raw_data={A,B};
for i=1:2
    current_points = raw_data{i};
    x_jitter=x(i) + (rand(size(current_points)) - 0.5) * 0.5;
    scatter(x_jitter, current_points, 10, ...
        'MarkerFaceColor', [0.7, 0.7, 0.7], ...
        'MarkerEdgeColor', [0.7, 0.7, 0.7], ...
        'MarkerFaceAlpha', 0.5, ...          % 填充透明度（0为完全透明，1为不透明）
        'MarkerEdgeAlpha', 0.5, ...          % 边缘透明度
        'LineWidth', 0.5);
end
uistack(Err1, 'top');
uistack(Err2, 'top');

