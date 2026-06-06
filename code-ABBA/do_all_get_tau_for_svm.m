% 在每个session内，以2分钟为一个bin，计算tau，用于后续的svm分类
clear
close all
fig_sign=0;
for icon = 1
    if icon == 1
        directories_ABBA_WT
        mark = 'WT';
       
    else
        directories_ABBA_APP_all
        mark = 'APP';
    end
    TimeBinSize_TAU=1; % 秒
    file_out=['tau_for_SVM.mat']; 
    fig_dir=['H:\group data_APP_WT\cofiring\SVM\'];
    mkdir(fig_dir);
    TTlist = 'TTList_dCA1_pyr.txt';
    daily_results = struct( 'best_C', [], 'best_KernelScale', []); % 交叉验证SVM的参数
    
    for ns = 1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        fprintf('处理会话 %d\n', ns);
        load(['Data_spikes.mat'],'spikes_all');
        load('Processing_Data_video.mat','Process_Data_video');
        load('Cells_info.mat','PeakRate');
        if  length(spikes_all{1})<15
            test_accuracy_day = [nan,nan,nan,nan];
            n_cellpair0_day = [nan,nan,nan,nan];
            save(file_out,'test_accuracy_day','n_cellpair0_day')
           continue 
        end
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
            TimeBin1 = [minTimeBin:120:maxTimeBin];
            if length(TimeBin1)==10
                TimeBin1=[TimeBin1,maxTimeBin];
            end
            Time_start = TimeBin1(1:10);
            Time_end = TimeBin1(2:11);
            tau_session=[];
            for ind_min =1:length(Time_start) % 在每个一分钟内求tau
                TimeBins = Time_start(ind_min):TimeBinSize_TAU:Time_end(ind_min)-TimeBinSize_TAU; % 注意减去了TimeBinSize_TAU以确保不超出范围
                sBinned_tau={};
                % 对s_tmp中的每个信号分binning
                s_tmp=spikes_all{nseg}(ind_nc,2);
                for nc = 1:length(s_tmp)
                    s = s_tmp{nc}; % spike的时间
                    [~, MnsMatchTimeBins] = histc(s, TimeBins);
                    binnedSpike = zeros(length(TimeBins), 1); % 初始化存储binned均值的数组
                    
                    % 遍历每个bin并计算平均值（跳过空bin）
                    for j = 1:length(TimeBins)
                        idx = MnsMatchTimeBins == j; % 找到当前bin中的索引
                        if any(idx)
                            binnedSpike(j) = sum(idx); % spike的个数
                        end
                    end
                    % 将binned均值存储在cell数组中
                    sBinned_tau{nc} = binnedSpike;
                end
                
                % 使用组合函数来计算tau（Kendall's tau-b）
                % 注意：MATLAB没有直接的itertools.combinations函数，但我们可以使用perms或nchoosek结合循环来实现
                tauVecSingle = []; % 初始化存储tau值的数组
                % 遍历两两细胞对的所有组合
                cell_pair_id = nchoosek(1:length(sBinned_tau), 2); %去掉不放电cell之后剩下的cell的编号
                cell_pair_id_ori = ind_nc(cell_pair_id) ;% 原始的cell id
                for i=1:size(cell_pair_id,1)
                    c = cell_pair_id(i,:);
                    c1 = sBinned_tau{c(1)};
                    c2 = sBinned_tau{c(2)};
                    [tau, ~] = corr(c1, c2, 'Type', 'Kendall'); % 计算Kendall's tau-b
                    tauVecSingle = [tauVecSingle, tau]; % 将tau值添加到数组中
                end
                tau_session(:,ind_min) = tauVecSingle';% 一个session内的结果 每一列是2分钟，每一行是一个细胞对
            end
            tau_all = [tau_all,tau_session];
        end
        save(file_out,'tau_all','cell_pair_id','cell_pair_id_ori');
    end
end