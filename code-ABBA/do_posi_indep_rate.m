% 计算空间无关的放电率（position-tuning independent rate）
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
plot_sign=1;

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
            %             directories_ABBA_APP
            directories_ABBA_APP_all
        end
        mark = 'APP';
        A=1;
    end

    TTlist = 'TTList_dCA1_pyr.txt';
    file_input = 'ratemap_no_mooth.mat';
    file_out='position_independent_rate';
    for ns =1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        load(file_input,'Ncell','PeakRate','RateMap','mapAxis','timeMaps','sessions','F','pathCoord','timeStamps');
        load(['Data_spikes.mat'],'spikes_all');
        r_obs_all={};
        r_exp_all={};
        r_pti_all={};
        PTI_rate_map_all={};
        obs_pti_cor=[];spatial_corr=[];
        for nseg=1:4
            %for nc=1:Ncell 
            for nc=11 %例子，wt第一天 cell 11    
                rateMap = RateMap{nseg,nc};
                x = pathCoord{nseg,1};
                y = pathCoord{nseg,2};
                t = pathCoord{nseg,3};
                timeMap = timeMaps{nseg,1};
                vfs=25;
                spike_times = spikes_all{nseg}{nc,2};
                [r_obs, r_exp, r_pti, binEdges] = compute_pti( ...
                    rateMap, x, y, t, spike_times, mapAxis, vfs,plot_sign);
                
                PTI_rate_map = compute_pti_ratemap(timeMap, ...
                    x, y, t, r_pti, binEdges, mapAxis, vfs);
                
                r_obs_all{nseg}(nc,:)=r_obs;
                r_exp_all{nseg}(nc,:)=r_exp;
                r_pti_all{nseg}(nc,:)=r_pti;
                PTI_rate_map_all{nc,nseg}=PTI_rate_map;
                
                r_pti(isnan(r_pti))=0;
                obs_pti_cor(nc,nseg) = corr(r_pti',r_obs');
                
                % 计算原始空间ratemap和PIR的ratemap的相关性
                [r_nonnan,c_nonnan] = find(~isnan(rateMap));
                [r_nonnan_pti,c_nonnan_pti] = find(~isnan(PTI_rate_map));
                RC_nonnan = intersect([r_nonnan_pti,c_nonnan_pti],[r_nonnan,c_nonnan],'rows');
                r = RC_nonnan(:,1);
                c = RC_nonnan(:,2);
                % NaN firing rate in unvisited bins
                x=diag(rateMap(r,c),0);  % x=X(r,c,i); returns matrix    从ratemap中提取出对应r、c的值
                y=diag(PTI_rate_map(r,c),0);
                spatial_corr(nc,nseg) = corr(x,y);
            
                
                if plot_sign ==1
                    figure('Units','normalized','Position',[-0.4 0.2 0.28 0.2]);
                    subplot(1,2,1);
                    rateMap_sm = smooth2a(rateMap,2);
                    h=imagesc(rateMap_sm,[0,max(max(rateMap))/1.5]); colormap jet
                    set (h,'alphadata',~isnan (rateMap))
                    subplot(1,2,2);
                    h=imagesc(PTI_rate_map,[0,max(max(rateMap))/1.5]); colormap jet
                    set (h,'alphadata',~isnan (rateMap))
                    
                end
            end
        end
        
        save(file_out,'r_obs_all','r_exp_all','r_pti_all','PTI_rate_map_all','obs_pti_cor');
        fprintf('第 %d 天计算结束...\n', ns);
    end
end