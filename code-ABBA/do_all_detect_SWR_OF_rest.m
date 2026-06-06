%% 计算每天SWR发生率； 把所有数据汇总； 用汇总后的数据统计

clc
close all
clear all
%SWRs相关参数
Fs=4:2:100;
sampfreq = 2000;
gap =1;% tov detect SWRs
plot_sign=1;
for icon=1:2
    if icon == 1
        directories_ABBA_WT
        mark = 'WT';
    else
        directories_ABBA_APP_all
        mark = 'APP';
    end
    file_path=['H:\group data_APP_WT\OF_rest\SWR\',mark,'\'];
    mkdir(file_path);
    file_outall=[file_path,'SWRsrate_replayrate_',mark,'.mat'];
    file_out=['SWR_rate.mat'];
    for ns =1:isession
        path_ns=path{ns};
        cd(path_ns);
        disp(strcat('Start calculating the SWR events in ',path_ns));
        load('Data_eventTS.mat','Ts_sleep');
        %取SWRs的时间点
        SWRnum=[];SWRrate=[];SWRdura={};Ripp_ind={};tt_eeg_ses={};
        
        Csclist = CSClist_CA1{ns};
        nCSC= length(Csclist);
        all_removeeg = {};
        %所有remove的eeg点
        EEG = cell(0);
        time_EEG=cell(0);
        for ncsc = 1:nCSC
            filnam = strcat('CSC',num2str(Csclist(ncsc)),'.ncs');
            if(exist(filnam)~=2)
                continue
            end
            % read EEG
            [sample,tt, ~,~] = loadCSC_new_cz_allrecording(filnam);
            EEG{ncsc,1}=sample;
            time_EEG{ncsc}=tt;
            % kick out time points
            [removinds,~] = fixeeg_cz(sample,gap);  % find out-of-range EEG
            removeeg=unique(removinds);
            all_removeeg{ncsc}=removeeg;
        end
        choseeeg = 1;
        time_EEG_chosen = time_EEG{choseeeg};%挑选其中一导的时间
        
        for ses=1:5
            EEG_seg={};
            ind = find(time_EEG_chosen>Ts_sleep(ses,1) & time_EEG_chosen<Ts_sleep(ses,2));
            tt_seg=time_EEG_chosen(ind);
            for ncsc = 1:nCSC
                EEG_seg{ncsc,1} = EEG{ncsc}(ind);
            end
            
            [RippleOnsetIndex, RippleOffsetIndex] = DetectRipples_v4(EEG_seg,sampfreq);%SWR时长大于15ms(0.015s),std=3 对所有导找SWRs
            time_ripples={};
            for i = 1:length(RippleOnsetIndex)
                time_ripples{i} = tt_seg(RippleOnsetIndex(i):1:RippleOffsetIndex(i));
            end
            %compute SWR rate of rest time
            SWRnum(ses) = length(time_ripples);
            SWRrate(ses) = SWRnum(ses)./(tt_seg(end)-tt_seg(1));
            SWRdura{ses} = (RippleOffsetIndex-RippleOnsetIndex+1)/2; % SWR的持续时间，单位毫秒
            Ripp_ind{ses}(:,1) = RippleOnsetIndex;
            Ripp_ind{ses}(:,2) = RippleOffsetIndex;
            tt_eeg_ses{ses} = tt_seg;
            
            %=====画SWR图=====%
            if plot_sign
                bin_max = max(SWRdura{ses})*2; %最长的SWR的bin数
                num_swr=length(RippleOnsetIndex);
                num_seq_row = 3; % 一幅图5列seq
                num_seq_col_max = 3;%一幅图最多3行
                num_seq_col = ceil(num_swr/num_seq_row); % 一共多少行
                num_fig = ceil(num_seq_col/num_seq_col_max); % 一共画几张图
                if num_seq_col>num_seq_col_max
                    num_seq_col=num_seq_col_max;
                end
                num_figure=0;
                for nswr=1:length(RippleOnsetIndex)
                    nseq_fig = mod(nswr,num_seq_row*num_seq_col_max); %在图中的第几个seq
                    if nseq_fig ==0;
                        nseq_fig=num_seq_row*num_seq_col_max;
                    end
                    if nseq_fig==1
                        num_figure= num_figure+1;
                        figure(num_figure)
                        set(gcf,'unit','centimeter','position',[-50 0 15*num_seq_row 10*num_seq_col]);
                    end
                    time_length = 1000;% 最多画多少毫秒的长度
                    temp_length = floor((time_length - SWRdura{ses}(nswr))/2)*2;%前后延长多少个bin
                    if (RippleOnsetIndex(nswr)-temp_length)<0
                        temp_length = RippleOnsetIndex(nswr)-1;
                    end
                    if RippleOffsetIndex(nswr)+temp_length>length(tt_seg)
                        temp_length = length(tt_seg)-RippleOffsetIndex(nswr);
                    end
                    subplot(num_seq_col,num_seq_row,nseq_fig); hold on
                    eeg_ind0 = [RippleOnsetIndex(nswr)-temp_length:RippleOffsetIndex(nswr)+temp_length];
                    eeg_ind = [RippleOnsetIndex(nswr):RippleOffsetIndex(nswr)];
                    
                    for ee=1:length(EEG_seg)
                        plot(tt_seg(eeg_ind0),EEG_seg{ee}(eeg_ind0)+ee*600,'k');
                        plot(tt_seg(eeg_ind),EEG_seg{ee}(eeg_ind)+ee*600,'r')
                    end
                    xlim([tt_seg(eeg_ind0(1)),tt_seg(eeg_ind0(end))]);
                    xticks([tt_seg(eeg_ind0(1)),tt_seg(eeg_ind0(end))-0.05]);
                    set(gca,'xticklabel',[0,time_length]);
                    xlabel('Time(ms)');
                    ylim([-1000,length(EEG)*750])
                    
                    if  nseq_fig== num_seq_row*num_seq_col_max || nswr==num_swr
                        fig_name = ['Rat',num2str(Ind_Rat(ns)),'-day',num2str(ns),'-session',num2str(ses),'-fig',num2str(num_figure)];
                        fig_dir = ['H:\group data_APP_WT\OF_rest\SWR\',mark,'\'];
                        saveas(gcf,[fig_dir,fig_name,'.png']);
                        close all
                    end
                end
            end
        end
        cd(path_ns);
        save(file_out,'SWRnum','SWRrate','SWRdura','Ripp_ind','tt_eeg_ses')
    end
end
