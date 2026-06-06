% 筛选在SWR期间的放电，计算不同SWR期间放电序列的相关性
clc
close all
clear all
plot_sign=1;
for icon=2
    if icon == 1
        directories_ABBA_WT
        mark = 'WT';
    else
        directories_ABBA_APP_all
        mark = 'APP';
    end
    file_out='SWR_spk_cor_v4.mat'; 
    
    for ns =1:isession
        path_ns=path{ns};
        cd(path_ns);
        TTList= 'TTList_dCA1_pyr.txt';
        fid=fopen(TTList);
        if (fid == -1)
            warning([ 'Could not open tfile ' TTList]);
            continue
        end
        load('SWR_rate.mat','Ripp_ind','tt_eeg_ses','SWRdura');
        load('Data_spikes_rest.mat','spikes_all');
        
        spk_cell_n={};% 每个cell在每个SWR内发放的几个spk
        spk_n={};% 每个SWR内一共出现了几个spike
        spk_cor={}; % SWR内spk 序列做互相关
        SWR_ind_keep={};% 筛选spk个数以及SWR时长之后剩下的SWR_ind
        if plot_sign==1
            f1 = figure;
            set(gcf,'unit','centimeter','position',[-50 0 50 8]);
            f2 = figure;
            set(gcf,'unit','centimeter','position',[-50 10 50 8]);
        end
        for nseg=1:5
            spk_cell_n0=[]; % 每个cell在每个SWR内发放的几个spk
            spk_n0=[]; % 每个SWR内一共出现了几个spike
            for nswr = 1:size(Ripp_ind{ nseg},1)
                for nc=1:size(spikes_all{1},1)
                    spk_ts = spikes_all{nseg}{nc,2};
                    inda=Ripp_ind{nseg}(nswr,1);
                    indb=Ripp_ind{nseg}(nswr,2);
                    tta = tt_eeg_ses{nseg}(inda); % SWR开始
                    ttb = tt_eeg_ses{nseg}(indb); % SWR结束
                    spk_cell_n0(nc,nswr)= length(find(spk_ts>=tta & spk_ts<=ttb));
                end
                spk_n0(nswr,1) = sum(spk_cell_n0(:,nswr));
            end

            spk_lim=1; %筛选放电个数
            SWRdura_lim1=15;
            SWRdura_lim2=400; %筛选SWR时长
            
            ind_keep=find(spk_n0> spk_lim & SWRdura{nseg}<SWRdura_lim2 & SWRdura{nseg}>SWRdura_lim1);
            spk_cell_n_keep = spk_cell_n0(:,ind_keep);
            spk_cell_n{nseg} = spk_cell_n_keep;
            spk_n{nseg} = spk_n0(ind_keep);
            SWR_ind_keep{nseg}=ind_keep;
            
            % 计算SWR内的细胞序列之间的相关性
            spk_cor0=nan(size(spk_cell_n_keep,2));
            for i=1:size(spk_cell_n_keep,2)
                a = spk_cell_n_keep(:,i);
                for j=i+1:size(spk_cell_n_keep,2)
                    b = spk_cell_n_keep(:,j);
                    spk_cor0(i,j)=corr(a,b);
                end
            end
            spk_cor{nseg}=spk_cor0;
            
            if plot_sign==1
                figure(f1);
                subplot(1,5,nseg);
                imagesc(spk_cor0,[-0.3,0.8]);
                yticks([1,size(spk_cor0,1)]);
                xlabel('SWRs ID'); ylabel('SWRs ID');
                colormap jet;
                title([nanmean(nanmean(spk_cor0))]);
                
                % SWR期间，每个细胞的放电个数，绘制成热图
                figure(f2);
                subplot(1,5,nseg);
                imagesc(spk_cell_n_keep,[0,5]);
                
                xlabel('SWRs ID'); ylabel('Cell ID');
                colormap jet;
                title(['Rest Ses',num2str(nseg)]);
            end
        end
        
        if plot_sign==1
            fig_dir = ['H:\group data_APP_WT\OF_rest\SWR\'];
            fig_name1=[mark,'-Rat',num2str(Ind_Rat(ns)),'-day',num2str(day_Rat(ns)),'-corr'];
            fig_name2=[mark,'-Rat',num2str(Ind_Rat(ns)),'-day',num2str(day_Rat(ns)),'-spkNUM'];
            saveas(f1,[fig_dir,fig_name1],'png');
            saveas(f2,[fig_dir,fig_name2],'png');
        end
        save(file_out,'spk_cell_n','spk_n','spk_cor','SWR_ind_keep');
        if ns==6
           a=1; 
        end
    end
end
%%
clc
close all
clear all
for icon=1:2
    if icon == 1
        directories_ABBA_WT
        mark = 'WT';
    else
        directories_ABBA_APP_all
        mark = 'APP';
    end
    file_in='SWR_spk_cor';
    % file_in='SWR_spk_cor_v4'; % 这些后续的版本，APP都没上升趋势
    
    spk_cor_mean=nan(isession,5);
    spk_cor_each_mean={};
    Ncell_day=[];
    for ns =1:isession
        path_ns=path{ns};
        cd(path_ns);
        TTList= 'TTList_dCA1_pyr.txt';
        fid=fopen(TTList);
        if (fid == -1)
            warning([ 'Could not open tfile ' TTList]);
            continue
        end
        load(file_in,'spk_cell_n','spk_n','spk_cor');
        load('Cells_info','Ncell')
        Ncell_day(ns,1)=Ncell;
        for nseg=1:5
            spk_cor_mean(ns,nseg) = nanmean(nanmean(spk_cor{nseg})) ;
            
            % 计算每个SWR期间reactive和其他reactive的平均相关性（即SWR1和其他n个算相关性之后，把这n个数据平均）
            temp1 = spk_cor{nseg};
            temp2 = temp1';
            A_clean = temp1;
            A_clean(isnan(temp1)) = 0;  % A中的NaN→0
            B_clean = temp2;
            B_clean(isnan(temp2)) = 0;  % B中的NaN→0
            C_ignore = A_clean + B_clean;
            mask_both_nan = isnan(temp1) & isnan(temp2);  % 逻辑矩阵，双NaN位置为1
            C_preserve = C_ignore;
            C_preserve(mask_both_nan) = NaN;
            spk_cor_each_SWR_mean = nanmean(C_preserve);
            spk_cor_each_mean{ns,nseg} = spk_cor_each_SWR_mean';
        end
    end
    if icon==1
        spk_cor_mean_wt = spk_cor_mean;
        spk_cor_each_mean_wt = spk_cor_each_mean;
        day_Rat_wt = day_Rat;
        Ncell_day_wt=Ncell_day;
    else
        spk_cor_mean_app = spk_cor_mean;
        spk_cor_each_mean_app = spk_cor_each_mean;
        day_Rat_app = day_Rat;
        Ncell_day_app=Ncell_day;
    end
end
cell_lim_wt=0;
cell_lim_app=10;
color_WT=[27,128,188]/255;color_APP=[159,0,0]/255;
figure('Units','normalized','Position',[0.1 0.2 0.7 0.3]);
for nseg=1:5
    mean1=[];mean2=[];std1=[];std2=[]; SPSS=[];
    for nday=1:6
        inda = find(day_Rat_wt==nday & Ncell_day_wt>=cell_lim_wt);
        indb = find(day_Rat_app==nday & Ncell_day_app>=cell_lim_app);
        A = spk_cor_mean_wt(inda,nseg);
        B = spk_cor_mean_app(indb,nseg);
        A(isnan(A))=[];      B(isnan(B))=[];
        mean1(nday,1) =  mean(A);    mean2(nday,1) =  mean(B);
        std1(nday,1) = std(A)/sqrt(length(A));    std2(nday,1) = std(B)/sqrt(length(B));
        %     SPSS = [SPSS;[[ones(length(A),1),ones(length(A),1)*nday,A];[ones(length(B),1)*2,ones(length(B),1)*nday,B]]];
    end
    
    %============绘制折线图===========%
    subplot(1,5,nseg)
    hold on
    h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
    h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
    ylabel('Correlation (reavtive in SWR)');
    ylim([0,0.5]);
    title12={'WT','App'};
    lhd = legend([h1 h2 ],title12,'Location','best');
    set(lhd,'Box', 'off')
    set(gca,'fontsize',15);
    xlabel('Day');
    title(['Rest ses',num2str(nseg)]);
end

%=======以单个SWR期间的重激活作为样本=======%
figure('Units','normalized','Position',[0.1 0.2 0.7 0.3]);
for nseg=1:5
    mean1=[];mean2=[];std1=[];std2=[]; SPSS=[];
    for nday=1:6
        inda = find(day_Rat_wt==nday & Ncell_day_wt>=cell_lim_wt);
        indb = find(day_Rat_app==nday & Ncell_day_app>=cell_lim_app);
        A = spk_cor_each_mean_wt(inda,nseg);
        B = spk_cor_each_mean_app(indb,nseg);
        AA=[];BB=[];
        for i=1:length(A)
            AA = [AA;A{i}];
        end
        for i=1:length(B)
            BB = [BB;B{i}];
        end
        A=AA; B=BB;
        A(isnan(A))=[];      B(isnan(B))=[];
        mean1(nday,1) =  mean(A);    mean2(nday,1) =  mean(B);
        std1(nday,1) = std(A)/sqrt(length(A));    std2(nday,1) = std(B)/sqrt(length(B));
        %     SPSS = [SPSS;[[ones(length(A),1),ones(length(A),1)*nday,A];[ones(length(B),1)*2,ones(length(B),1)*nday,B]]];
    end
    
    %============绘制折线图===========%
    subplot(1,5,nseg)
    hold on
    h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
    h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
    ylabel('Correlation (reavtive in SWR)');
    ylim([0,0.5]);
    title12={'WT','App'};
    lhd = legend([h1 h2 ],title12,'Location','best');
    set(lhd,'Box', 'off')
    set(gca,'fontsize',15);
    xlabel('Day');
    title(['Rest ses',num2str(nseg)]);
end
%% 把后四个rest session合在一起统计
close all
figure('Units','normalized','Position',[0.1 0.2 0.15 0.3]);
mean1=[];mean2=[];std1=[];std2=[]; SPSS=[];
y_wt=[];y_app=[];x_wt=[];x_app=[];
for nday=1:6
    inda = find(day_Rat_wt==nday & Ncell_day_wt>=cell_lim_wt);
    indb = find(day_Rat_app==nday & Ncell_day_app>=cell_lim_app);
    A = spk_cor_mean_wt(inda,2:5); A=reshape(A,size(A,1)*size(A,2),1);
    B = spk_cor_mean_app(indb,2:5);B=reshape(B,size(B,1)*size(B,2),1);
    A(find(A>0.34|A<0))=[];  B(find(B>0.34|B<0))=[]; % 去除异常值
    y_wt=[y_wt;A];                       y_app=[y_app;B];
    x_wt=[x_wt;ones(length(A),1)*nday];  x_app=[x_app;ones(length(B),1)*nday];
    A(isnan(A))=[];      B(isnan(B))=[];
    mean1(nday,1) =  mean(A);    mean2(nday,1) =  mean(B);
    std1(nday,1) = std(A)/sqrt(length(A));    std2(nday,1) = std(B)/sqrt(length(B));
    %     SPSS = [SPSS;[[ones(length(A),1),ones(length(A),1)*nday,A];[ones(length(B),1)*2,ones(length(B),1)*nday,B]]];
end
%============绘制折线图===========%
hold on
h1 = errorbar([1:length(mean1)],mean1,std1,'s-','LineWidth',2,'Color',color_WT);
h2 = errorbar([1:length(mean1)],mean2,std2,'s-','LineWidth',2,'Color',color_APP);
ylabel('Correlation (reavtive in SWR)');
ylim([0,0.4]);
xlim([0.5,6.5]);
title12={'WT','App'};
lhd = legend([h1 h2 ],title12,'Location','best');
set(lhd,'Box', 'off')
set(gca,'fontsize',15);
xlabel('Day');
title(['Rest']);

%=====App 线性拟合图===========%
y=y_app; x=x_app;
[b_app1,~,~,~,~]=regress(y,[ones(length(y),1),x]);% b(1)是截距，b（2）是斜率;stats(3)是显著性
model = fitlm(x, y);  % 构建线性回归模型
[r,p]=corr(x,y);
R_squared_app1 = model.Rsquared.Ordinary;  % R平方值
F_stat_app1 = model.anova.F(1);  % F统计量
p_value_app1 = model.anova.pValue(1);  % 显著性p值
% 绘图用，为了把重合的点错开，但是计算显著性仍然用原始横坐标
x=x+0.1 * randn(size(x));
plot_liner_regres_v2(x,y,color_APP)
ylim([0,0.4]);
title(['App ','r=',num2str(r),' p=',num2str(p)]);


%=====WT 线性拟合图===========%
y=y_wt; x=x_wt;
x(find(y>0.34|y<0))=[];
y(find(y>0.34|y<0))=[];
ind_temp = isnan(y);
x(ind_temp)=[];
y(ind_temp)=[];
[b_wt1,~,~,~,~]=regress(y,[ones(length(y),1),x]);% b(1)是截距，b（2）是斜率;stats(3)是显著性
model = fitlm(x, y);  % 构建线性回归模型
[r,p]=corr(x,y);
R_squared_wt1 = model.Rsquared.Ordinary;  % R平方值
F_stat_wt1 = model.anova.F(1);  % F统计量
p_value_wt1 = model.anova.pValue(1);  % 显著性p值
% 绘图用，为了把重合的点错开，但是计算显著性仍然用原始横坐标
x=x+0.1 * randn(size(x));
plot_liner_regres_v2(x,y,color_WT)
ylim([0,0.4]);
title(['WT ','r=',num2str(r),' p=',num2str(p)]);

