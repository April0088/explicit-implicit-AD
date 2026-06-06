% 绘制所有的放电率热图，第一节是ABBA，第二节是AAA
clear
clc
Session='ABBA';
for icon = 2
    if icon == 1
        directories_ABBA_WT
        mark = 'WT';
    else
        directories_ABBA_APP_all %AD
        mark = 'APP';
    end
    TTlist = 'TTList_dCA1_int.txt';
    file_input = 'Cells_info.mat';
    fig_path = ['H:\group data_APP_WT\remapping\ratemap\',Session,'\',mark,'\'];
    
    for ns = 4:isession
        path_ns = path{ns};
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        
        Rat=['Rat ',num2str(Ind_Rat(ns))];
        session_sequence=[1,2,3,4];
        rotate_para=0;
        load(file_input)
        scale_max=0;
        sLength = 140;
        bins = 35;
        fig_name = [Rat,'_ns',num2str(ns),'_day',num2str(day_Rat(ns))];
        fig_path1 = [fig_path,Rat,'\'];
        mkdir(fig_path1);
        close all
        fig_num=0;
        for ncell = 1:Ncell
            if mod (ncell,7)==1
               h1=figure('Position', [200 0 600 1000]);
               fig_num=fig_num+1;
            end
            for nseg=1:size(RateMap,1)
                if mod (ncell,7)~=0
                fign = nseg+size(RateMap,1)*(mod (ncell,7)-1);
                else
                    fign = nseg+size(RateMap,1)*6;
                end
                subplot(7,size(RateMap,1),fign)
                fr_max = max(max(RateMap{nseg,ncell}));
                h=imagesc(RateMap{nseg,ncell},[0,max(2,fr_max)*0.9]);
                set(h,'alphadata',~isnan(RateMap{nseg,ncell}))
                axis xy
                axis square
                colormap jet
                title(['Cell',num2str(ncell)]);
                set(gca,'XColor', [1 1 1], 'YColor', [1 1 1]) %让坐标轴变成白色
                if ~isempty(fieldProp{nseg}{ncell})
                   field_size = [fieldProp{nseg}{ncell}.size];
                   title(['Cell',num2str(ncell),' size:',num2str(field_size(1))]);
                end
            end
            if mod (ncell,7)==0 | ncell==Ncell
                saveas(gca,[fig_path1,fig_name,'_fig',num2str(fig_num),'.png'])
                
            end
        end

    end
end
%%
% AAA绘制所有的放电率热图
clear
clc
Session='AAA';
for icon = 2
    if icon == 1
        directories_AAA_WT
        mark = 'WT';
    else
        directories_AAA_APP %AD
        mark = 'APP';
    end
    TTlist = 'TTList_dCA1_int.txt';
    file_input = 'Cells_info.mat';
    fig_path = ['H:\group data_APP_WT\remapping\ratemap\',Session,'\',mark,'\'];
    
    for ns = 1:isession
        path_ns = path{ns};
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        
        Rat=['Rat ',num2str(Ind_Rat(ns))];
        session_sequence=[1,2,3,4];
        rotate_para=0;
        load(file_input)
        scale_max=0;
        sLength = 100;
        bins = 25;
        fig_name = [Rat,'_ns',num2str(ns),'_day',num2str(day_Rat(ns))];
        fig_path1 = [fig_path,Rat,'\'];
        mkdir(fig_path1);
        close all
        fig_num=0;
        for ncell = 1:Ncell
            if mod (ncell,7)==1
               h1=figure('Position', [200 0 450 1000]);
               fig_num=fig_num+1;
            end
            for nseg=1:size(RateMap,1)
                if mod (ncell,7)~=0
                fign = nseg+size(RateMap,1)*(mod (ncell,7)-1);
                else
                    fign = nseg+size(RateMap,1)*6;
                end
                subplot(7,size(RateMap,1),fign)
                fr_max = max(max(RateMap{nseg,ncell}));
                h=imagesc(RateMap{nseg,ncell},[0,max(2,fr_max)]);
                set(h,'alphadata',~isnan(RateMap{nseg,ncell}))
                axis xy
                colormap jet
                set(gca,'XColor', [1 1 1], 'YColor', [1 1 1]) %让坐标轴变成白色
            end
            if mod (ncell,7)==0 | ncell==Ncell
                saveas(gca,[fig_path1,fig_name,'_fig',num2str(fig_num),'.png'])
            end
        end
    end
end