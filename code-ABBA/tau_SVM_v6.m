% =========================================================================
% 把A1B1作为训练集，B2A2作为测试集
% 交叉参数测试
% 只sort训练集
% =========================================================================
% --- 设置 ---
clear
close all
load('H:\group data_APP_WT\cofiring\SVM\SVM_para.mat', 'daily_results');
for icon = 1:2
    if icon == 1
        directories_ABBA_WT; % 假设此函数设置了WT组的路径
        mark = 'WT';
    else
        directories_ABBA_APP_all; % 假设此函数设置了APP组的路径
        mark = 'APP';
    end
    
    TTlist = 'TTList_dCA1_pyr.txt';
    file_out_suffix = '_v6.mat'; % 使用新的文件名以防覆盖
    file_out = ['cofiring_SVM' file_out_suffix];
    for ns = 1:isession
        path_ns = path{ns};
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        % --- 数据加载和基本检查 ---
        tau_all=[];
        load('tau_for_SVM.mat', 'tau_all');
        if length(tau_all)<2
            test_accuracy_day = [nan;nan;nan;nan];
            n_cellpair0_day = [nan;nan;nan;nan];
            save(file_out,'test_accuracy_day','n_cellpair0_day')
            fprintf('跳过会话 %d: 不存在tau\n', ns);
            continue;
        end
        
        % --- NaN数据预处理--
        tau_imputed = local_preprocess_tau_data(tau_all);
        if isempty(tau_imputed)
            fprintf('跳过会话 %d: 预处理后无有效数据\n', ns);
            continue;
        end
        
        % --- 特征数量循环 ---
        percent_ForSVM_all = [1.0]; % 1.0代表全部特征
        test_accuracy_day = nan(length(percent_ForSVM_all), length(daily_results));
        n_cellpair0_day = nan(length(percent_ForSVM_all), 1);
        B2_test_accuracy_day = nan(length(percent_ForSVM_all), length(daily_results));
        A2_test_accuracy_day = nan(length(percent_ForSVM_all), length(daily_results));

        for p = 1:length(percent_ForSVM_all)
            percent_ForSVM = percent_ForSVM_all(p);                
            % --- 1. 定义固定的训练集和测试集 (基于时间) ---
            % 训练集: 前半段 (A1, B1)
            train_indices_time = 1:20; 
            % 测试集: 后半段 (B2, A2) - 注意，在原始tau矩阵中，A2在B2之后
            % 所以我们需要拼接 B2 (21-30) 和 A2 (31-40)
            test_indices_time = 21:40;
            tau_train_raw = tau_imputed(:, train_indices_time);
            tau_test_raw = tau_imputed(:, test_indices_time);            
            y_train = [ones(10, 1); 2*ones(10, 1)]; % A1=1, B1=2
            y_test = [2*ones(10, 1); ones(10, 1)];  % B2=2, A2=1

            % --- 2. 【关键】只在训练集上进行特征选择 ---
            a_train = mean(tau_train_raw(:, y_train==1), 2);
            b_train = mean(tau_train_raw(:, y_train==2), 2);
            [~, sorted_indices_train] = sort(a_train - b_train, 'descend');
            
            if percent_ForSVM < 1.0
                n_over = round(size(tau_imputed, 1) * percent_ForSVM / 2);
                if n_over < 1; continue; end
                selected_indices = [sorted_indices_train(1:n_over); sorted_indices_train(end-n_over+1:end)];
            else
                selected_indices = (1:size(tau_imputed, 1))';
            end
            
            % --- 3. 创建最终的训练和测试矩阵 ---
            X_train = tau_train_raw(selected_indices, :)';
            X_test = tau_test_raw(selected_indices, :)';
            
            % 4. 使用所有WT超参数进行测试
            accuracies=[];
            B2_accuracies=[];
            A2_accuracies=[];
            for n = 1:length(daily_results)
                best_C = daily_results(p, n).best_C;
                best_KS = daily_results(p, n).best_KernelScale;                
                model = fitcsvm(X_train, y_train, 'KernelFunction', 'rbf','Standardize', true, ...
                    'BoxConstraint', best_C, 'KernelScale', best_KS);                
                y_pred = predict(model, X_test);
                accuracies(1, n) = sum(y_pred == y_test) / length(y_test);                
                B2_accuracies(1,n) = sum(y_pred(1:10) == y_test(1:10)) / length(y_test(1:10));
                A2_accuracies(1,n) = sum(y_pred(11:20) == y_test(11:20)) / length(y_test(11:20));
            end            
            % 计算当前百分比下的平均准确率
            test_accuracy_day(p, :) = accuracies;
            B2_test_accuracy_day(p, :) = B2_accuracies;
            A2_test_accuracy_day(p, :) = A2_accuracies;
            n_cellpair0_day(p, 1) = size(X_train, 2); % 记录实际使用的特征数
        end
        save(file_out, 'test_accuracy_day', 'n_cellpair0_day','B2_test_accuracy_day','A2_test_accuracy_day');
    end
end

% --- 辅助函数 ---
function tau_imputed = local_preprocess_tau_data(tau_all)
    % 此函数封装了NaN处理代码
    max_nan_ratio = 0.3;
    nan_ratio_per_row = sum(isnan(tau_all), 2) / size(tau_all, 2);
    rows_to_keep = nan_ratio_per_row <= max_nan_ratio;
    tau_filtered = tau_all(rows_to_keep, :);
    
    if isempty(tau_filtered)
        tau_imputed = [];
        return;
    end

    session_boundaries = [1, 10; 11, 20; 21, 30; 31, 40];
    tau_imputed = tau_filtered;
    for i = 1:size(tau_imputed, 1)
        for j = 1:size(session_boundaries, 1)
            start_col = session_boundaries(j, 1);
            end_col = session_boundaries(j, 2);
            segment = tau_imputed(i, start_col:end_col);
            nan_in_segment = isnan(segment);
            if any(nan_in_segment)
                seg_mean = mean(segment, 'omitnan');
                if isnan(seg_mean); seg_mean = 0; end
                global_indices = find(nan_in_segment) + start_col - 1;
                tau_imputed(i, global_indices) = seg_mean;
            end
        end
    end
end
