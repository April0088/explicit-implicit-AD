clear

vfs = 25; % sampling frequenct of the video

%% load event file for time stamps
load('Data_eventTS.mat');
load('Data_video.mat');
savefile=('Processing_Data_video.mat');
ndirs = size(Ts_begin,1);
Process_Data_video={};
%% 将矿场之外的反光点线性差值
data_video_cha = data_video;

if ndirs == 3  % AAA session
    posx = data_video(:,2);
    ind1 = find(abs(posx)>50);
    cha_part = find_continue_part(ind1);
    for ii =1:size(cha_part,2)
        num = length(cha_part{ii});
        if num>1
            icha1 = cha_part{ii}(1)-1;
            icha2 = cha_part{ii}(num)+1;
        else
            icha1 = cha_part{ii}(1)-1;
            icha2 = cha_part{ii}(1)+1;
        end
        if icha2 > length(posx) || icha1 <1
            pre_x2 = ones(1,length(cha_part{ii})+2)*NaN;
        else
            pre_x2 =interp1([icha1,icha2],[posx(icha1),posx(icha2)],icha1:icha2,'linear');
        end
        posx(icha1+1:icha2-1) = pre_x2(2:num+1);
    end
    posy = data_video(:,3);
    ind2 = find(abs(posy)>50);
    cha_part = find_continue_part(ind2);
    for ii =1:size(cha_part,2)
        num = length(cha_part{ii});
        if num>1
            icha1 = cha_part{ii}(1)-1;
            icha2 = cha_part{ii}(num)+1;
        else
            icha1 = cha_part{ii}(1)-1;
            icha2 = cha_part{ii}(1)+1;
        end
        if icha2 > length(posy) || icha1 <1
            pre_y2 = ones(1,length(cha_part{ii})+2)*NaN;
        else
            pre_y2 =interp1([icha1,icha2],[posy(icha1),posy(icha2)],icha1:icha2,'linear');
        end
        posy(icha1+1:icha2-1) = pre_y2(2:num+1);
    end
    data_video_cha(:,2)= posx;
    data_video_cha(:,3)= posy;
else % ABBA
    for i = 1 : 4
        if i == 1 || i ==4
            posx = data_video{i}(:,2);
            ind1 = find(abs(posx)>50);
            cha_part = find_continue_part(ind1);
            for ii =1:size(cha_part,2)
                num = length(cha_part{ii});
                if num>1
                    icha1 = cha_part{ii}(1)-1;
                    icha2 = cha_part{ii}(num)+1;
                else
                    icha1 = cha_part{ii}(1)-1;
                    icha2 = cha_part{ii}(1)+1;
                end
                if icha2 > length(posx) || icha1 <1
                    pre_x2 = ones(1,length(cha_part{ii})+2)*NaN;
                else
                    pre_x2 =interp1([icha1,icha2],[posx(icha1),posx(icha2)],icha1:icha2,'linear');
                end
                posx(icha1+1:icha2-1) = pre_x2(2:num+1);
            end
            posy = data_video{i}(:,3);
            ind2 = find(abs(posy)>50);
            cha_part = find_continue_part(ind2);
            for ii =1:size(cha_part,2)
                num = length(cha_part{ii});
                if num>1
                    icha1 = cha_part{ii}(1)-1;
                    icha2 = cha_part{ii}(num)+1;
                else
                    icha1 = cha_part{ii}(1)-1;
                    icha2 = cha_part{ii}(1)+1;
                end
                if icha2 > length(posy) || icha1 <1
                    pre_y2 = ones(1,length(cha_part{ii})+2)*NaN;
                else
                    pre_y2 =interp1([icha1,icha2],[posy(icha1),posy(icha2)],icha1:icha2,'linear');
                end
                posy(icha1+1:icha2-1) = pre_y2(2:num+1);
            end
            data_video_cha{i}(:,2)= posx;
            data_video_cha{i}(:,3)= posy;
        else
            posx = data_video{i}(:,2);
            posy = data_video{i}(:,3);
            dist = data_video{i}(:,4);
            ind = find(dist>64); % 假设B session为原型矿场，假定圆的半径
            cha_part = find_continue_part(ind);
            for ii =1:size(cha_part,2)
                num = length(cha_part{ii});
                if num>1
                    icha1 = cha_part{ii}(1)-1;
                    icha2 = cha_part{ii}(num)+1;
                else
                    icha1 = cha_part{ii}(1)-1;
                    icha2 = cha_part{ii}(1)+1;
                end
                if icha2 > length(posy) || icha1 <1
                    pre_x2 = ones(1,length(cha_part{ii})+2)*NaN;
                    pre_y2 = ones(1,length(cha_part{ii})+2)*NaN;
                else
                    pre_x2 =interp1([icha1,icha2],[posx(icha1),posx(icha2)],icha1:icha2,'linear');
                    pre_y2 =interp1([icha1,icha2],[posy(icha1),posy(icha2)],icha1:icha2,'linear');
                end
                posx(icha1+1:icha2-1) = pre_x2(2:num+1);
                posy(icha1+1:icha2-1) = pre_y2(2:num+1);
            end
            data_video_cha{i}(:,2)= posx;
            data_video_cha{i}(:,3)= posy;
        end
    end
end

%% 平滑轨迹
if ndirs == 3
    for i=1:ndirs              % AAA session
        [value,ind(1)]=min(abs(data_video(:,1) - Ts_begin(i,1)));
        [value,ind(2)]=min(abs(data_video(:,1) - Ts_begin(i,2)));
        Process_Data_video{i}=data_video_cha(ind(1):ind(2),:);
        figure
        plot(Process_Data_video{i}(:,2),Process_Data_video{i}(:,3))
        %         axis([-50,50,-50,50])
        
        posx=Process_Data_video{i}(:,2);
        posy=Process_Data_video{i}(:,3);
        % Moving window mean filter
        for cc = 8:length(posx)-7
            posx(cc) = nanmean(posx(cc-7:cc+7));
            posy(cc) = nanmean(posy(cc-7:cc+7));
        end
        Process_Data_video{i}(:,2)=posx;
        Process_Data_video{i}(:,3)=posy;
        figure
        plot(Process_Data_video{i}(:,2),Process_Data_video{i}(:,3))
        %         axis([-50,50,-50,50])
    end
else
    for i=1:ndirs              % ABBA session
        [value,ind(1)]=min(abs(data_video{i}(:,1) - Ts_begin(i,1)));
        [value,ind(2)]=min(abs(data_video{i}(:,1) - Ts_begin(i,2)));
        Process_Data_video{i}=data_video_cha{i}(ind(1):ind(2),:);
        figure
        plot(Process_Data_video{i}(:,2),Process_Data_video{i}(:,3))
        %         axis([-50,50,-50,50])
        
        posx=Process_Data_video{i}(:,2);
        posy=Process_Data_video{i}(:,3);
        % Moving window mean filter
        for cc = 8:length(posx)-7
            posx(cc) = nanmean(posx(cc-7:cc+7));
            posy(cc) = nanmean(posy(cc-7:cc+7));
        end
        Process_Data_video{i}(:,2)=posx;
        Process_Data_video{i}(:,3)=posy;
        figure
        plot(Process_Data_video{i}(:,2),Process_Data_video{i}(:,3))
        %         axis([-50,50,-50,50])
    end
end

save(savefile,'Process_Data_video','vfs','scale','centre');
