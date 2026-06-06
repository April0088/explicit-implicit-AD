% to load video data in the session of open field
% edited by Guo M and Zheng C, 2021/5/16

clear

vfs = 25; % sampling frequenct of the video
x_length = 100 ; % real size of the field
y_length = 100 ; % real size of the field
cir_r = 0.6370 ; %

%% load event file for time stamps
load('Data_eventTS.mat');
ndirs = size(Ts_begin,1);

%% load video data for running sessions
fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 1; % Extracted X
fieldSelection(3) = 1; % Extracted Y
fieldSelection(4) = 0; % Extracted Angel
fieldSelection(5) = 0;  % Targets
fieldSelection(6) = 0; % Points

% Do we return header 1 = Yes, 0 = No.
extractHeader = 0;
% 5 different extraction modes, see help file for Nlx2MatVT
extractMode = 1; % Extract all data
[t, x, y] = Nlx2MatVT('VT1.nvt',fieldSelection,extractHeader,extractMode);
t = t ./ 1000000;
% find zero values and set to nan
ind = find(x==0 & y==0);
x(ind) = nan;
y(ind) = nan;
% combine behavior data in running sessions
ind_begin = {};
for nd = 1:ndirs
    [~,ind1] = min(abs(t-Ts_begin(nd,1)));
    [~,ind2] = min(abs(t-Ts_begin(nd,2)));
    x_nd = x(ind1:ind2);
    y_nd = y(ind1:ind2);
    t_nd = t(ind1:ind2);
    
    % Threshold for how far a rat can move (100cm/s)
    % assume the pixel to cm scale is about 0.3
    threshold = (120/vfs)/0.3;
    [x_nd,y_nd] = removeBadTracking_cz(x_nd,y_nd,threshold);
    
    % find zero values and do interpolation
    ind = find((x_nd>0 & y_nd>0)==1 & ~isnan(x_nd) & ~isnan(y_nd));
    x_nd_fix=interp1(t_nd(ind),x_nd(ind),t_nd,'linear');
    y_nd_fix=interp1(t_nd(ind),y_nd(ind),t_nd,'linear');
    
    % save the fixed behavior data into the original data
    x(ind1:ind2) = x_nd_fix;
    y(ind1:ind2) = y_nd_fix;
    
    ind_begin{nd} = [ind1:ind2];
end
centre = [];
centre2 = [];
%% find scale x scale y center from session 1, and apply to all other sessions
for nd = 1
    if nd == 2 || nd ==4
        % plot1: original data
        x_all = x(ind_begin{nd});
        y_all = y(ind_begin{nd});
        if isempty(centre)
            figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
            scatter(x_all,y_all,'+')
            axis square
            title('Point out the diagonal corners without the box')
            [xc,yc] = ginput(2); % point to the corners without the box
            centre = [mean(xc),mean(yc)];
            % plot2: centered data
            x_center = x_all-centre(1);
            y_center = y_all-centre(2);
            figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
            scatter(x_all,y_all,'+')
            axis square
            title('Point out the corners on x-axis')
            
            % plot3: rotate with the angle
            [xc,yc] = ginput(2);
            l = sqrt((xc(1)-xc(2)).^2+(yc(1)-yc(2)).^2);
            pos = (yc(2)-yc(1))./l;    % get the cos
            ang_rotate = asin(pos);
            disp(ang_rotate);
            x_rot = x_center.*cos(ang_rotate)+y_center.*sin(ang_rotate);  % clockwise rotation
            y_rot = -x_center.*sin(ang_rotate)+y_center.*cos(ang_rotate);
            figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
            scatter(x_rot,y_rot,'+')
            axis square
            title('Point out the inner-border of the field')
            
            % plot4: re-scale the field, and find pixels/cm scale
            [xc,yc] = ginput(4); % point to the corners without the box
            scale_x = x_length./(max(xc)-min(xc));
            scale_y = y_length./(max(yc)-min(yc));
            x2 = x_rot*scale_x;
            y2 = y_rot*scale_y;
        else
            % session1和session4用一样的标准来归一化
            x_center = x_all-centre(1);
            y_center = y_all-centre(2);
            x_rot = x_center.*cos(ang_rotate)+y_center.*sin(ang_rotate);  % clockwise rotation
            y_rot = -x_center.*sin(ang_rotate)+y_center.*cos(ang_rotate);
            x2 = x_rot*scale_x;
            y2 = y_rot*scale_y;
        end
        
        figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
        hold on
        plot(x2,y2,'.')
        plot([-x_length,x_length],[0,0],'r--')
        plot([0,0],[-y_length,y_length],'r--')
        plot([-x_length/2,x_length/2],[-y_length/2,-y_length/2],'r')
        plot([-x_length/2,x_length/2],[y_length/2,y_length/2],'r')
        plot([-x_length/2,-x_length/2],[-y_length/2,y_length/2],'r')
        plot([x_length/2,x_length/2],[-y_length/2,y_length/2],'r')
        hold off
        xlim([-x_length/2-20,x_length/2+20])
        ylim([-y_length/2-20,y_length/2+20])
        axis square
        x0 = -x_length/2:x_length/4:x_length/2;
        y0 = -y_length/2:y_length/4:y_length/2;
        set(gca,'FontSize',20);
        set(gca, 'XTick',x0);
        set(gca, 'YTick',x0);
        xlabel('Position x (cm)')
        ylabel('Position y (cm)')
        title('Video data for whole session')
        set(gcf,'PaperPositionMode','auto')
        close all
        
        %%%=========================save the data ======================%%%
        timelimit=120;   % ==== set the highest speed, unit is cm/s ====
        vel0=speed2D(x2,y2,t); %velocity in cm/s
        vel0(vel0>=timelimit) = 0.5*(vel0(circshift((vel0>=timelimit),-3)) + vel0(circshift((vel0>=timelimit),3)));
        t_nd = t(ind_begin{nd});
        % distance to centre
        dist_all = [];
        for i = 1 : length(x2)
            a = [x2(i),y2(i)];
            b = [0,0];
            dist = norm(a-b);
            dist_all = [dist_all,dist];
        end
        % Save data per session
        data_video{nd} = nan(length(t_nd),6);
        data_video{nd}(:,1)=t_nd';% Timesample
        data_video{nd}(:,2)=x2';% x*scale-center
        data_video{nd}(:,3)=y2';% y*scale-center
        data_video{nd}(:,4)=dist_all'; % distance to center
        data_video{nd}(:,5)=vel0;% velocity
        data_video{nd}(:,6)=nan;% angle velocity

    else
        % session2和session3用一样的标准来归一化
        x_all = x(ind_begin{nd});
        y_all = y(ind_begin{nd});
        if isempty(centre2)
            figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
            scatter(x_all,y_all,'+')
            axis square
            title('Point out the quadrant vertex')
            [xc,yc] = ginput(4); % point to the corners without the box
            centre2 = [mean(xc),mean(yc)];
        end
        
        % plot2: centered data
        x_center = x_all-centre2(1);
        y_center = y_all-centre2(2);
        
        % plot3: rotate with the angle get from session1
        x_rot = x_center.*cos(ang_rotate)+y_center.*sin(ang_rotate);  % clockwise rotation
        y_rot = -x_center.*sin(ang_rotate)+y_center.*cos(ang_rotate);
        
        % plot4: re-scale the field, and find pixels/cm scale
        x2 = x_rot*scale_x;
        y2 = y_rot*scale_y;
        figure('Units','normalized','Position',[0.1 0.1 0.6 0.8]);
        hold on
        plot(x2,y2,'.')
        plot([-x_length,x_length],[0,0],'r--')
        plot([0,0],[-y_length,y_length],'r--')
        plot([-x_length/2,x_length/2],[-y_length/2,-y_length/2],'r')
        plot([-x_length/2,x_length/2],[y_length/2,y_length/2],'r')
        plot([-x_length/2,-x_length/2],[-y_length/2,y_length/2],'r')
        plot([x_length/2,x_length/2],[-y_length/2,y_length/2],'r')
        hold off
        xlim([-x_length/2-20,x_length/2+20])
        ylim([-y_length/2-20,y_length/2+20])
        axis square
        x0 = -x_length/2:x_length/4:x_length/2;
        y0 = -y_length/2:y_length/4:y_length/2;
        set(gca,'FontSize',20);
        set(gca, 'XTick',x0);
        set(gca, 'YTick',x0);
        xlabel('Position x (cm)')
        ylabel('Position y (cm)')
        title('Video data for whole session')
        set(gcf,'PaperPositionMode','auto')
        close all
        
        %%%=========================save the data ======================%%%
        timelimit=120;   % ==== set the highest speed, unit is cm/s ====
        vel0=speed2D(x2,y2,t); %velocity in cm/s
        vel0(vel0>=timelimit) = 0.5*(vel0(circshift((vel0>=timelimit),-3)) + vel0(circshift((vel0>=timelimit),3)));
        t_nd = t(ind_begin{nd});
        % distance to centre
        dist_all = [];
        for i = 1 : length(x2)
            a = [x2(i),y2(i)];
            b = [0,0];
            dist = norm(a-b);
            dist_all = [dist_all,dist];
        end
        % Save data per session
        data_video{nd} = nan(length(t_nd),6);
        data_video{nd}(:,1)=t_nd';% Timesample
        data_video{nd}(:,2)=x2';% x*scale-center
        data_video{nd}(:,3)=y2';% y*scale-center
        data_video{nd}(:,4)=dist_all'; % distance to center
        data_video{nd}(:,5)=vel0;% velocity
        data_video{nd}(:,6)=nan;% angle velocity
        
    end
end

scale=[scale_x,scale_y];% scale of x and y
save('Data_video.mat','data_video','vfs','scale','ang_rotate','centre','centre2','ind_begin');
