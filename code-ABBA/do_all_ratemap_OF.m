% 计算旷场内的位置域信息，包括个数，大小，峰值放电率等
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

for icon =1:2
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
    fig_dir=['H:\group data_OF\cell\',mark,'_ratemap_',Session,'\'];
    mkdir(fig_dir);
    TTlist = 'TTList_dCA1_pyr.txt';
    file_output = 'Cells_info.mat';
    % limit the posang_ontrack with running speed > vel_threshold
    vel_threshold = 5; % cm/s  %IMPORTANT: change to 0 if do not want to limit
    
    for ns =1:isession
        path_ns = path{ns};
        if isempty(path_ns)
            continue
        end
        cd(path_ns);
        if ~exist(TTlist)
            continue
        end
        
        rotate_para=0;
        [RateMap,timeMaps,PeakRate,mapAxis,sessions,F,pathCoord,timeStamps,Ncell] ...
            = equalPlot_RM_step1_rotate_cz_wxl(TTlist,25*60,rotate_para,session_number,sLength,bins);
        %             [RateMap,timeMaps,PeakRate,mapAxis,sessions,F,pathCoord,timeStamps] ...
        %                 = equalPlot_RM_step1_rotate_cz(infile,20*60,rotate_para,session_number,sLength,bins);
        save(file_output,'TTlist','path_ns','Ncell','PeakRate','RateMap','mapAxis','timeMaps','sessions','F','pathCoord','timeStamps');
        
        %% calculate the place fields
        p.minNumBins = 5;
        % Bins with rate at p.fieldTreshold * peak rate and higher will be
        % considered as part of a place field
        % Actually, these parameters are not used- see peakAll variable at approximately line 429, which is
        % set to 1 Hz
        p.fieldTreshold = 0.2;
        % Lowest field rate in Hz.
        p.lowestFieldRate = 1;
        peak_all = 0.5;
        
        % place field for single laps
        nFields = cell(size(RateMap,1),1);
        fieldProp = cell(size(RateMap,1),1);
        fieldcomp = cell(size(RateMap,1),1);
        mkdir([path_ns,'placefiled']);
        for nseg = 1:size(RateMap,1)
            for nc = 1:Ncell
                ratemap_nc = RateMap{nseg,nc};
                binWidth = sLength/bins;
                fig_name=[fig_dir,'Rat ',num2str(Ind_Rat(ns)),'-Day ',num2str(ns),'-Cell ',num2str(nc),'-Session ',num2str(nseg)];
                [nFields_nc,fieldProp_nc,fieldcomp_nc] = placefield_2D_WN_v3(ratemap_nc,p,binWidth,peak_all,plot_sign,fig_name);
                nFields{nseg}(1,nc) = nFields_nc;
                fieldProp{nseg}{1,nc} = fieldProp_nc;
                fieldcomp{nseg}{1,nc} = fieldcomp_nc;
                close all
            end
        end
        
        save(file_output,'p','peak_all',...
            'nFields','fieldProp','fieldcomp','-append');
        clear RateMap timeMaps fieldProp
        %% Return
        cd ../
    end
end