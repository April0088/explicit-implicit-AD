% placefield identifies the placefields in the firing map. It returns the
% number of placefields and the location of the peak within each
% placefield.
% original file: placefield_circ_CZ
%
%
%
% map           Rate map
% pTreshold     Field treshold
% pBins         Minimum number of bins in a field
% mapAxis       The map axis
%
function [nFields,fieldProp,fieldcomp] = placefield_2D_WN_v3(map,p,binWidth,peak_all,plot_sign,fig_name)

% Counter for the number of fields
nFields = 0;
% Field properties will be stored in this struct array
fieldProp = [];
% Allocate memory to the arrays
[N_x,N_y]= size(map);
% Array that contain the bins of the map this algorithm has visited
% Field component
fieldcomp = zeros(N_x,N_y);
visited = zeros(N_x,N_y);
nanInd = isnan(map);
visited(nanInd) = 1;
visited2 = visited;
All_area_ind = ~isnan(map);
% fill Nan with 0.001
map_im = map;
map_im(isnan(map_im))=0.001;

% find peak
peak_max = max(max(map_im));
% Check if peak rate is high enough
if peak_max < p.lowestFieldRate
    return;
end

if p.fieldTreshold * peak_max < peak_all
    peak_all = p.fieldTreshold * peak_max;
end
visited2(map < peak_all) = 1;

%% find multiple peak
n=0;
Niob_ind= {};
% 找到合适的连通域
FR_threshold = peak_max*0.2; % 位置域的放电率阈值
n=n+1;
BW(:,:,n)= logical(map_im>FR_threshold);% 二值化的ratemap
% Label connected components in BW image
[L(:,:,n),Niob(n)] = bwlabel(BW(:,:,n),4); % defult 8-connected iobects

ind1 = [];%储存大的iob
ind2 = [];%储存小的iob
for iob = 1:Niob(n)
    iobInd = find(L(:,:,n)==iob);
    iobAreabin = length(iobInd);%计算第iob个区域的由多少个bin组成
    if iobAreabin >= p.minNumBins
        ind1 = [ind1,iob];
    else
        ind2 = [ind2,iob];
    end
end
% 位置域范围太小的不要
if isempty(ind1)
    return;
end

Niob_ind{1,n} = ind1;
Niob_ind{2,n} = ind2;
iobNum(1,n) = length(ind1);% 大的个数
iobNum(2,n) = length(ind2);% 小的个数


iobNum_max_ind = find(iobNum(1,:)==max(iobNum(1,:)));
indx = iobNum_max_ind(1); % 合适的连通域的ind

map_BW1 = BW(:,:,1);% 用来大于阈值的BW
map_BW2 = BW(:,:,indx);% 用来找峰的合适的BW

L1 = L(:,:,1);
if ~isempty(Niob_ind{2,1})
    for i = Niob_ind{2,1}
        L1(L1==i)=0;
    end
end

visited2 = double((~logical(L1)) | logical(visited2));
ntemp = 0;
for iob = Niob_ind{1,1}
    ntemp = ntemp+1;
    % Ind = find(L1==iob);
    AllPeak1(ntemp) = max(map_im(L1==iob)); % iobective peak
end

if isempty(AllPeak1)
    return
end
[p1sort,pI1]= sort(AllPeak1,'descend');
L2 = L(:,:,indx);

for iob =1: Niob(indx)
    %     Ind = find(L2==iob);
    AllPeak2(iob) = max(map_im(L2==iob)); % iobective peak
end
[p2sort,pI2]= sort(AllPeak2,'descend');

AllPeak = union(AllPeak1,AllPeak2,'stable');
[Allpeak, pI]= sort(AllPeak,'descend');

nPeak = length(Allpeak);
ntemp = 0;
% 按峰值大小，排序component
for iob = 1:nPeak
    if iob<=Niob(indx) && length( find(L2==pI2(iob)) )>=p.minNumBins
        ntemp = ntemp +1;
        fieldcomp(L2==pI2(iob))=ntemp;
        peakInd(ntemp) = iob;
    else if ~isempty(find(p1sort==Allpeak(iob)));
            ntemp = ntemp +1;
            pItemp = find(p1sort==Allpeak(iob));
            fieldcomp(L1==pI1(pItemp))=ntemp;
            peakInd(ntemp) = iob
        end
    end
end

nPeak = ntemp;
visited3 = double((~logical(L2))|logical(visited2));

if length(find(visited3==0))< length(find(visited2==0))
    visited3 = ~xor(visited3,visited2);
    px = nan(1,nPeak);
    py = nan(1,nPeak);
    for ipeak = 1:nPeak
        [px(ipeak),py(ipeak),~] = find(map_im==Allpeak(peakInd( ipeak)));% 所有峰值的坐标
    end
    
    % visited3 = double((~logical(L2))|logical(visited2));
    % comparesize = length(visited3)-length(visited2);
    
    [vx,vy,~] = find(visited2==0);% 所有还没有遍历的点的坐标
    nv_num = length(vx);
    pX = repmat(px,nv_num,1);
    pY = repmat(py,nv_num,1);
    vX = repmat(vx,1,nPeak);
    vY = repmat(vy,1,nPeak);
    dist = sqrt( (vX-pX).^2+(vY-pY).^2 );
    [~,belongto] = min(dist,[],2);
    for i = 1:nv_num
        fieldcomp(vx(i),vy(i)) = belongto(i);
    end
end

%% calculate property
x=[];y=[]; % 质心的横纵坐标
for iob = 1:nPeak
    iobInd = find(fieldcomp==iob);
    iobPeak = max(map_im(iobInd)); % iobective peak
    iobAreabin = length(map_im(iobInd));%计算第iob个区域的由多少个bin组成
    if iobAreabin < p.minNumBins % iob太小换下一个iob
        continue;
    end
    [iobPeaky,iobPeakx] = find(map_im==iobPeak);% iobective peak bin
    iobPeakX = iobPeakx * binWidth;iobPeakY = iobPeaky * binWidth; % iobective peak bin real position
    
    sum_x=0;    sum_y=0;    areaFR=0; %初始化
    for i=1:N_x
        for j=1:N_y
            if fieldcomp(i,j)==iob
                sum_x=sum_x+j * map(i,j);  %计算第Ｋ区域的横坐标总和
                sum_y=sum_y+i * map(i,j);  %计算第Ｋ区域的纵坐标总和
                areaFR=areaFR+map(i,j);    %计算该区域的放电率总和
            end
        end
    end
    
    iobCOMx=sum_x/areaFR;  %计算第Ｋ区域的质心横坐标
    iobCOMy=sum_y/areaFR;%计算第Ｋ区域的质心纵坐标
    iobAveFR = areaFR/iobAreabin;%计算平均放电率
    
    iobCOMX = iobCOMx * binWidth;% COM的实际x坐标
    iobCOMY = iobCOMy * binWidth;% COM的实际y坐标
    
    iobArea = iobAreabin * binWidth^2; % 实际面积
    
    % X or Y mean the real position(cm), and x or y unit is bin
    fieldProp = [fieldProp; ...
        struct('XY_COM',[iobCOMX,iobCOMY],...
        'peakRate',iobPeak,'XY_peak',[iobPeakX,iobPeakY],...
        'avgRate',iobAveFR,'size',iobArea,...
        'xy_COM',[iobCOMx,iobCOMy],... % com在第几个bin
        'xy_peak',[iobPeakx,iobPeaky], 'sizebin',iobAreabin)];
    nFields = nFields + 1;
    x(iob) = iobCOMx;
    y(iob) = iobCOMy;
end
%===========画图=================%
if plot_sign==1
    figure
    set(gcf,'Units','normalized','Position',[0.2 0.2 0.6 0.5]);
    subplot(1,2,1)
    h1=imagesc(fieldcomp);
    set(h1,'alphadata',~isnan(map));
    axis xy
    axis xy;axis square; axis off
    subplot(1,2,2)
    h2=imagesc(map);
    set(h2,'alphadata',~isnan(map));
    axis xy
    colormap(jet)
    title(['where the COM    ','PF:',num2str(round(peak_max,2)),'Hz'])
    axis xy;axis square;axis off
    hold on;
    for iob=1:nPeak
        subplot(1,2,2)
        plot(x(iob),y(iob),'k*', 'MarkerSize', 15, 'LineWidth', 2)
    end
    saveas(gca,[fig_name,'.png'])
    hold off
end
end






