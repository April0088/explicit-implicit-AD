% 计算位置域重合面积，筛选位置域部分重合的cell  
% 然后计算它们的时序CCG  然后计算CCG是否由theta调制，得到theta_modulation_index
clc
% clear all
plot_sign=1; %是否画每一个细胞对的theta调制图
area_lim=0; 
for icon = 1:2
    if icon == 1
        directories_ABBA_WT; 
        mark = 'WT';
    else
        directories_ABBA_APP_v2; 
        mark = 'APP';
    end
    file_input = 'Cells_info.mat';
    file_input2 = 'Data_spikes.mat'; % spike的时间信息
    file_input3 = 'Processing_Data_video.mat'; % 行为学信息
    file_input4 = 'Cell_feature.mat';
    file_out='space_overlap_CCG'; % 存在空间上位置域重合的cell的CCG
    for ns = 1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        load(file_input,'fieldcomp','Ncell','fieldProp','RateMap');
        load(file_input2,'spikes_all');
        load(file_input3,'Process_Data_video','vfs');
        load(file_input4,'pfFeature');
        if Ncell<10
           continue
        end
        % 0. 位置域质心的位置
        COM_x=cell(1,length(fieldcomp));
        COM_y=cell(1,length(fieldcomp));
        for nseg=1:length(fieldcomp)
            for nc=1:Ncell
                if length(pfFeature{nseg}(nc).COM_position)==1
                    COM_x{nseg} = [COM_x{nseg};nan];
                    COM_y{nseg} = [COM_y{nseg};nan];
                else
                    COM_x{nseg} = [COM_x{nseg};pfFeature{nseg}(nc).COM_position(1)-70]; 
                    COM_y{nseg} = [COM_y{nseg};pfFeature{nseg}(nc).COM_position(2)-70];
                end                
            end
        end
        
        % 1. 找到位置域有重合的细胞对
        overlap_cell_id={};overlap_ratio={};
        for nseg=1:length(fieldcomp)
            map = fieldcomp{nseg};
            overlapMatrix = calculateOverlapArea(map); % 计算细胞对之间的重合面积
            
            [overlap_cell_id0,overlap_ratio0 ]= findOverlappingPairs(overlapMatrix,area_lim);
            overlap_cell_id{nseg}=overlap_cell_id0; % 重合比例大于10%的细胞对id
            overlap_ratio{nseg}=overlap_ratio0; % 所有放过电的细胞对的重合比例和实际重合面积。(面积单位是bin数，一个bin边长4cm)
        end
        
        % 2. 计算上述细胞对的spike之间的时间差,并得到直方图
        CCG_nseg={};theta_index_nseg={};
        r_theta_index_nseg={};
        theta_index_db_nseg={};
        r_theta_index_db_nseg={};
        
        for nseg=1:length(fieldcomp)            
            overlap_cell_id0 = overlap_cell_id{nseg};
            spk_ts=spikes_all{nseg}(:,2);
            spk_xy = spikes_all{nseg}(:,3);
            max_lag=0.5;%spike之间最大的时间差范围，即只计算0.5s内的
            comx = COM_x{nseg};
            comy = COM_y{nseg};
            
            [spike_time_diffs,spk_n] = SpikeTimeDiffs(spk_ts, spk_xy, overlap_cell_id0, max_lag,comx,comy); % 得到spike之间的时间差,以及spike的数量
            % spike_time的两列，对应着两个方向（C1到C2，C2到C1）
            bin_width=0.005;
            % bin_width=0.011;% 画例子用的，wt1 的第11行
            [histogram_data, bin_centers] = getSmoothedHistogramData(spike_time_diffs, bin_width, max_lag); % 高斯平滑之后的CCG
            
            % 计算细胞对被theta调制的程度
            
            theta_index=[];% theta频段最大功率作为index
            r_theta_index=[];%相对的theta功率，用theta除以总功率
            theta_index_db=[];% 将功率谱转换为dB单位，然后theta频段最大功率作为index
            r_theta_index_db=[];% 相对的，衡量theta功率比基线高出多少分贝
             
            for n=1:size(histogram_data,1) % 例子是wt组第一天的n=11
                if plot_sign==1
                    ffa=figure('position',[-350,300,350,300]);
                end
                for i=1:size(histogram_data,2)
                    FS=length(histogram_data{n,i}); % CCG的采样频率，由于我这里的CCG时间是1s的，其bin长度刚好是采样频率
                    [theta_index(n,i),r_theta_index(n,i),theta_index_db(n,i),r_theta_index_db(n,i)]...
                        =CCG_theta_modu_index(histogram_data{n,i},bin_centers,bin_width,FS,plot_sign);
                    ylim([0,0.1])
                    if plot_sign==1
                        figure(ffa);hold on
                        plot(bin_centers * 1000,histogram_data{n,i},'linewidth',1)
                        xlim([-250,250]);
                        set(gca,'fontsize',15);
                        xlabel('Time lag(ms)');
                        ylabel('CCG');
                    end
                end
                close all
            end
            spk_n_nseg{nseg}=spk_n;
            CCG_nseg{nseg}=histogram_data;
            theta_index_nseg{nseg}=theta_index;
            r_theta_index_nseg{nseg}=r_theta_index;
            theta_index_db_nseg{nseg}=theta_index_db;
            r_theta_index_db_nseg{nseg}=r_theta_index_db;            
        end
        
        save(file_out,'overlap_cell_id','overlap_ratio','CCG_nseg','overlap_cell_id','bin_centers',...
            'spk_n_nseg','theta_index_nseg','r_theta_index_nseg','theta_index_db_nseg','r_theta_index_db_nseg');        
        fprintf('完成会话 %d \n', ns);
    end
end

close all
figure;hold on
plot(bin_centers*1000,histogram_data{26,1});
plot(bin_centers*1000,histogram_data{26,2});
xlim([-200,200]);
legend({'Cell 1 to Cell 2','Cell 2 to Cell 1'});
xlabel('Time lag(ms)');
ylabel('Cross-correlogram');
set(gca,'fontsize',15);

figure;hold on
plot(bin_centers*1000,histogram_data{1,1});
plot(bin_centers*1000,histogram_data{1,2});
xlim([-200,200]);
legend({'Cell 1 to Cell 2','Cell 2 to Cell 1'});
xlabel('Time lag(ms)');
ylabel('Cross-correlogram');
set(gca,'fontsize',15);
%% 统计theta index的大小
close all
for icon = 1:2
    if icon == 1
        directories_ABBA_WT; % 假设此函数设置了WT组的路径
        mark = 'WT';
    else
        directories_ABBA_APP_v2; % 假设此函数设置了APP组的路径
        mark = 'APP';
    end
    file_tau_cell_pair = ['cell_pair_sorded.mat'];
    file_input = 'space_overlap_CCG.mat'; % 存在空间上位置域重合的cell的CCG
    theta_index_all=cell(4,1);
    r_theta_index_all=cell(4,1);
    theta_index_db_all=cell(4,1);
    r_theta_index_db_all=cell(4,1);
    
    for ns = 1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(file_input)
            continue
        end
        load(file_input,'spk_n_nseg','theta_index_nseg','r_theta_index_nseg','theta_index_db_nseg','r_theta_index_db_nseg','overlap_ratio');
        for nseg=1:4
            ind = find(overlap_ratio{nseg}(:,3)>=0.1 | overlap_ratio{nseg}(:,4)>=0.1);
            theta_index_all{nseg} = [theta_index_all{nseg};...
                [theta_index_nseg{nseg}(ind,:),ones(length(ind),1)*day_Rat(ns)]];         
            r_theta_index_all{nseg} = [r_theta_index_all{nseg};...
                [r_theta_index_nseg{nseg}(ind,:),ones(length(ind),1)*day_Rat(ns)]];            
            theta_index_db_all{nseg} = [theta_index_db_all{nseg};...
                [theta_index_db_nseg{nseg}(ind,:),ones(length(ind),1)*day_Rat(ns)]];  
            r_theta_index_db_all{nseg} = [r_theta_index_db_all{nseg};...
                [r_theta_index_db_nseg{nseg}(ind,:),ones(length(ind),1)*day_Rat(ns)]]; 
            
        end
    end
    
    if icon==1
        theta_index_wt=theta_index_all;
        r_theta_index_wt=r_theta_index_all;
        theta_index_db_wt=theta_index_db_all;
        r_theta_index_db_wt=r_theta_index_db_all;        
    else
        theta_index_app=theta_index_all;    
        r_theta_index_app=r_theta_index_all;
        theta_index_db_app=theta_index_db_all;
        r_theta_index_db_app=r_theta_index_db_all;
    end
end
%% 统计theta index值在两组之间的差异

close all
n_day=6;
for nseg=1:4
    mean1=[];mean2=[];std1=[];std2=[];SPSS=[];
    for ngroup=1:n_day
        ind1 = find(theta_index_wt{nseg}(:,3)==ngroup);
        ind2 = find(theta_index_app{nseg}(:,3)==ngroup);
        A = nanmean(theta_index_wt{nseg}(ind1,[2]),2);
        B = nanmean(theta_index_app{nseg}(ind2,[2]),2);
        mean1(ngroup,1) =  mean(A);
        mean2(ngroup,1) =  mean(B);
        std1(ngroup,1) = std(A)/sqrt(length(A));
        std2(ngroup,1) = std(B)/sqrt(length(B));
        SPSS=[SPSS;[ones(length(A),1),ones(length(A),1)*ngroup,A];...
            [ones(length(B),1)*2,ones(length(B),1)*ngroup,B]];
    end
    
    color_WT='b';
    color_APP='r';
    title12={'WT','APP'};
    figure('Units','normalized','Position',[0.2 0.2 0.15 0.3]);
    hold on
    h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
    h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
    % set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
    xlim([0.5,n_day+0.5])
    lhd = legend([h1 h2],title12,'Location','best');
    set(lhd,'Box', 'off')
    title(['S',num2str(nseg)]);
    set(gca,'fontsize',15);
    ylabel('Theta index');
    ylim([0,0.6]);
    xlabel('Day');
end


% 相对的theta功率，用theta除以总功率
for nseg=1:4
    mean1=[];mean2=[];std1=[];std2=[];SPSS=[];
    for ngroup=1:n_day
        ind1 = find(r_theta_index_wt{nseg}(:,3)==ngroup);
        ind2 = find(r_theta_index_app{nseg}(:,3)==ngroup);
        A = nanmean(r_theta_index_wt{nseg}(ind1,[2]),2);
        B = nanmean(r_theta_index_app{nseg}(ind2,[2]),2);
        A(isnan(A))=[];
        B(isnan(B))=[];
        mean1(ngroup,1) =  nanmean(A);
        mean2(ngroup,1) =  nanmean(B);
        std1(ngroup,1) = std(A)/sqrt(length(A));
        std2(ngroup,1) = std(B)/sqrt(length(B));
        SPSS=[SPSS;[ones(length(A),1),ones(length(A),1)*ngroup,A];...
            [ones(length(B),1)*2,ones(length(B),1)*ngroup,B]];
    end
    
    color_WT='b';
    color_APP='r';
    title12={'WT','APP'};
    figure('Units','normalized','Position',[0.2 0.2 0.15 0.3]);
    hold on
    h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
    h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
    % set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
    xlim([0.5,n_day+0.5])
    lhd = legend([h1 h2],title12,'Location','best');
    set(lhd,'Box', 'off')
    title(['S',num2str(nseg)]);
    set(gca,'fontsize',15);
    ylabel('Theta index');
    ylim([0.1,0.2]);
    xlabel('Day');
end

% 将功率谱转换为dB单位，然后theta频段最大功率作为index           
for nseg=1:4
    mean1=[];mean2=[];std1=[];std2=[];SPSS=[];
    for ngroup=1:n_day
        ind1 = find(theta_index_db_wt{nseg}(:,3)==ngroup);
        ind2 = find(theta_index_db_app{nseg}(:,3)==ngroup);
        A = nanmean(theta_index_db_wt{nseg}(ind1,[2]),2);
        B = nanmean(theta_index_db_app{nseg}(ind2,[2]),2);
        A(isinf(A))=[];
        B(isinf(B))=[];
        mean1(ngroup,1) =  nanmean(A);
        mean2(ngroup,1) =  nanmean(B);
        std1(ngroup,1) = std(A)/sqrt(length(A));
        std2(ngroup,1) = std(B)/sqrt(length(B));
        SPSS=[SPSS;[ones(length(A),1),ones(length(A),1)*ngroup,A];...
            [ones(length(B),1)*2,ones(length(B),1)*ngroup,B]];
    end
    
    color_WT='b';
    color_APP='r';
    title12={'WT','APP'};
    figure('Units','normalized','Position',[0.2 0.2 0.15 0.3]);
    hold on
    h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
    h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
    % set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
    xlim([0.5,n_day+0.5])
    lhd = legend([h1 h2],title12,'Location','best');
    set(lhd,'Box', 'off')
    title(['S',num2str(nseg)]);
    set(gca,'fontsize',15);
    ylabel('Theta index');
    ylim([-25,-10]);
    xlabel('Day');
end


% 相对的，衡量theta功率比基线高出多少分贝      
for nseg=1:4
    mean1=[];mean2=[];std1=[];std2=[];SPSS=[];
    for ngroup=1:n_day
        ind1 = find(r_theta_index_db_wt{nseg}(:,3)==ngroup);
        ind2 = find(r_theta_index_db_app{nseg}(:,3)==ngroup);
        A = nanmean(r_theta_index_db_wt{nseg}(ind1,[2]),2);
        B = nanmean(r_theta_index_db_app{nseg}(ind2,[2]),2);
        A(isnan(A))=[];
        B(isnan(B))=[];
        mean1(ngroup,1) =  nanmean(A);
        mean2(ngroup,1) =  nanmean(B);
        std1(ngroup,1) = std(A)/sqrt(length(A));
        std2(ngroup,1) = std(B)/sqrt(length(B));
        SPSS=[SPSS;[ones(length(A),1),ones(length(A),1)*ngroup,A];...
            [ones(length(B),1)*2,ones(length(B),1)*ngroup,B]];
    end
    
    color_WT='b';
    color_APP='r';
    title12={'WT','APP'};
    figure('Units','normalized','Position',[0.2 0.2 0.15 0.3]);
    hold on
    h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
    h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
    % set(gca,'xticklabel',{'A1-B2','B2-B3','B3-A4','A1-A4'});
    xlim([0.5,n_day+0.5])
    lhd = legend([h1 h2],title12,'Location','best');
    set(lhd,'Box', 'off')
    title(['S',num2str(nseg)]);
    set(gca,'fontsize',15);
    ylabel('Theta index');
    ylim([14,18]);
    xlabel('Day');
end

% 不考虑session，不考虑跨天
data1=[];data2=[];
for nseg=1:4
    data1 = [data1;r_theta_index_wt{nseg}];
    data2 = [data2;r_theta_index_app{nseg}];
end
A = nanmean(data1(:,[1,2]),2);
B = nanmean(data2(:,[1,2]),2);
[h,p,~,stats]=ttest2(A,B);
df = stats.df; % 自由度
t = stats.tstat;   %统计量
cohens_d = abs(stats.tstat) * sqrt(1/length(A) + 1/length(B));

mean_A =nanmean(A);
mean_B = nanmean(B);
sd_A = nanstd(A)/sqrt(length(A));
sd_B = nanstd(B)/sqrt(length(B));

color1=[6,157,255]/255;
color2=[197,39,45]/255;
x = [1,2];
figure
hold on
bar(x(1),mean_A,'FaceColor',color1,'EdgeColor', color1);
errorbar(1,mean_A,sd_A,'k', 'linestyle', 'none', 'lineWidth', 2);
bar(x(2),mean_B,'FaceColor',color2,'EdgeColor', color2);
errorbar(2,mean_B,sd_B,'k', 'linestyle', 'none', 'lineWidth', 2);
set(gca,'fontsize',15);
xticks([1,2])
xticklabels({'WT','App'})
ylabel('Normalized power');

%% 以下都是子函数
%%  计算细胞对之间的重合面积
function overlapMatrix = calculateOverlapArea(placefields)
    % 计算细胞的总数
    nCells = length(placefields);
    
    % 初始化一个 nCells x nCells 的矩阵来存储结果，用0填充
    overlapMatrix = zeros(nCells, nCells);
    
    % 使用双重循环遍历所有唯一的细胞对
    % 外层循环从第1个细胞到倒数第2个
    for i = 1:nCells
        % 内层循环从 i+1 开始，以避免重复计算 (i,j) 和 (j,i)，以及自身比较
        for j = (i + 1):nCells
            
            % 从元胞数组中提取两个细胞的位置野矩阵
            field_A = placefields{i};
            field_B = placefields{j};
            
            % 核心步骤：使用逻辑“与”(&)运算找到重合区域
            % 只有在field_A和field_B中对应元素都为1时，结果矩阵的元素才为1
            overlap_map = field_A & field_B;
            
            % 计算重合面积，即对结果矩阵的所有元素求和
            overlap_area = sum(overlap_map, 'all');
            
            % 将结果存入对称矩阵的相应位置
            overlapMatrix(i, j) = overlap_area;
            overlapMatrix(j, i) = overlap_area; % 因为是对称的
        end
    end
    
    % (可选) 计算对角线元素，即每个细胞自身位置野的面积
    for i = 1:nCells
        overlapMatrix(i, i) = sum(placefields{i}, 'all');
    end
    
end
%% 寻找重合面积大于25%的细胞对
function [overlap_cell_id,overlap_ratio] = findOverlappingPairs(overlapMatrix,threshold)
    if nargin < 2
        threshold = 0.25; % 25%
    end
    
    % 获取细胞总数
    nCells = size(overlapMatrix, 1);
    
    % 从对角线提取每个细胞自身的总面积
    totalAreas = diag(overlapMatrix);
    
    % 初始化一个空的数组，用于存储符合条件的细胞对索引
    overlap_cell_id = [];
    overlap_ratio = []; %细胞对重合面积比
    
    % 同样使用高效的循环遍历所有唯一的细胞对
    for i = 1:nCells
        for j = (i + 1):nCells
            
            % 提取细胞i和j的重合面积，以及它们各自的总面积
            overlap_area = overlapMatrix(i, j);
            area_i = totalAreas(i);
            area_j = totalAreas(j);
            
            % --- 核心判断逻辑 ---
            % 检查以确保总面积不为0，避免除零错误
            if area_i > 0 && area_j > 0
                
                % 计算重合面积分别占各自总面积的百分比
                ratio_i = overlap_area / area_i;
                ratio_j = overlap_area / area_j;
                overlap_ratio = [overlap_ratio; [i,j,ratio_i,ratio_j,overlap_area]]; 
                % 如果任何一个比例大于设定的阈值，则认为它们是重合的细胞对
                if ratio_i >= threshold || ratio_j >= threshold
                    % 将这对细胞的索引 [i, j] 添加到结果列表中
                    overlap_cell_id = [overlap_cell_id; i, j];
                end
            end
        end
    end
end

%% 计算上述细胞对的spike之间的时间差
function [spike_time_diffs,spk_n] = SpikeTimeDiffs(spikes_ts, spikes_xy, overlap_cell_id, max_lag,comx,comy)
    % 如果用户没有提供max_lag，设置一个默认值
    if nargin < 4
        max_lag = 0.5; % 默认计算前后500ms的时间差
    end
    
    % 获取重合细胞对的数量
    num_pairs = size(overlap_cell_id, 1);
    
    % 初始化一个元胞数组来存储每一对的结果
    spike_time_diffs = cell(num_pairs, 2);
    spk_n=[]; % 存在相邻spike时的spike的个数，注意，并不是两个细胞本身的是spike个数
    
    % 遍历每一对重合细胞
    for p = 1:num_pairs
        % 获取这对细胞的索引(第一列的记为细胞A，第二列的记为细胞B)
        cell_idx_A = overlap_cell_id(p, 1);
        cell_idx_B = overlap_cell_id(p, 2);
        
        % 获取这对细胞的位置野中心
        center_A = [comx(cell_idx_A), comy(cell_idx_A)];
        center_B = [comx(cell_idx_B), comy(cell_idx_B)];
        
        % 计算“场地参考向量”(从A指向B)
        V_field = center_B - center_A;
        
        % 提取这两个细胞的放电时间戳和位置
        spikes_A = spikes_ts{cell_idx_A}(:);
        spikes_B = spikes_ts{cell_idx_B}(:);
        spk_A_x = spikes_xy{cell_idx_A}(:,1);
        spk_A_y = spikes_xy{cell_idx_A}(:,2);
        spk_B_x = spikes_xy{cell_idx_B}(:,1);
        spk_B_y = spikes_xy{cell_idx_B}(:,2);
        
        
        % 初始化2个空的向量，用于存储这对细胞的时间差
        diffs_A_to_B = []; diffs_B_to_A=[];
        
        % 检查是否有spike，避免空数组错误
        if isempty(spikes_A) || isempty(spikes_B)
            spike_time_diffs{p} = []; % 如果任一细胞没有放电，则时间差为空
            continue; % 继续下一次循环
        end
        
        % 核心算法：以细胞A的每个spike为“参考点”
        for i = 1:length(spikes_A)
            ref_spike_time = spikes_A(i);
            
            % 寻找细胞B中，落在参考点前后 max_lag 窗口内的所有棘波
            % 这种筛选可以极大提升效率，避免不必要的计算
            ind_near = find(spikes_B > (ref_spike_time - max_lag) & ...
                       spikes_B < (ref_spike_time + max_lag));
            nearby_spikes = spikes_B(ind_near);
            
            % 如果找到了附近的spike
            if ~isempty(nearby_spikes)
                % 计算时间差 (B的时间 - A的时间)
                diffs = nearby_spikes - ref_spike_time;
                
                % 计算每对spike发生时动物的移动方向
                A_xy = [spk_A_x(i),spk_A_y(i)];
                B_xy = [spk_B_x(ind_near),spk_B_y(ind_near)];
                V_animal = B_xy - A_xy;
                
                % 用点积判断方向(两个向量角度查小于90°，可以认为动物是从A到B，大于90°则是从B到A)
                direction_score = V_animal * V_field'; 
                
                % direction_score大于0说明是从A位置域到B位置域，相反是从B到A
                ind_AB = find(direction_score>0);
                ind_BA = find(direction_score<0);
                
                % 将这次计算出的时间差追加到总列表中
                diffs_A_to_B = [diffs_A_to_B; diffs(ind_AB)];
                diffs_B_to_A = [diffs_B_to_A; diffs(ind_BA)];
            end
        end
        
        % 将这对细胞的所有时间差存储到输出的元胞数组中
        spike_time_diffs{p,1} = diffs_A_to_B;
        spike_time_diffs{p,2} = diffs_B_to_A;
        spk_n(p,1)=length(diffs_A_to_B);
        spk_n(p,2)=length(diffs_B_to_A);
    end
end

%% 计算并平滑spike之间的时间差直方图
function [smoothed_counts, bin_centers] = getSmoothedHistogramData(spike_time_diffs, bin_width_seconds, max_lag_seconds, smoothing_width_seconds)
    % 设置默认参数
    if nargin < 2 || isempty(bin_width_seconds)
        bin_width_seconds = 0.005; % 默认bin宽度为5ms
    end
    if nargin < 3 || isempty(max_lag_seconds)
        max_lag_seconds = 0.5; % 默认最大时间差为500ms
    end
    if nargin < 4 || isempty(smoothing_width_seconds)
        smoothing_width_seconds = 0.05; % 默认高斯平滑核宽度为50ms
    end

    % 获取细胞对的数量
    num_pairs = length(spike_time_diffs);
    
    % --- 新增：计算平滑窗口大小 (单位: bins) ---
    % round函数确保我们得到一个整数
    window_size_bins = round(smoothing_width_seconds / bin_width_seconds);

    % 初始化一个元胞数组来存储最终的平滑结果
    smoothed_counts = cell(num_pairs, 1);

    % 定义直方图的边界
    hist_edges = -max_lag_seconds:bin_width_seconds:max_lag_seconds;
    
    % 计算bin的中心点
    bin_centers = hist_edges(1:end-1) + (bin_width_seconds / 2);
    
    for i=1:size(spike_time_diffs,2)
        % 遍历每一对细胞的时间差数据
        for p = 1:num_pairs
            time_lags = spike_time_diffs{p,i};
            
            if isempty(time_lags)
                % 如果没有数据，返回全零向量
                counts = zeros(1, length(bin_centers));
            else
                % 1. 计算原始的直方图计数值
                counts = histcounts(time_lags, hist_edges);
            end
            
            % --- 新增：对计数值进行高斯平滑 ---
            % 2. 使用 smoothdata 函数
            if sum(counts) > 0 % 只在有数据时进行平滑
                smoothed_counts{p,i} = smoothdata(counts, 'gaussian', window_size_bins);
            else
                smoothed_counts{p,i} = counts; % 如果全是0，则无需平滑
            end
        end
    end
end