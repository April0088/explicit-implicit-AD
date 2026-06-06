% to load event timestamps
% edited by Guo M and Zheng C, 2021/5/28

clear
clc

file = 'Events.nev';

fieldSelection(1) = 1; % Timestamps
fieldSelection(2) = 0; % Event IDs
fieldSelection(3) = 0; % TTLs
fieldSelection(4) = 0; % Extras
fieldSelection(5) = 1;  % Event Strings

% Do we return header 1 = Yes, 0 = No.
extractHeader = 1;
% 5 different extraction modes, see help file for Nlx2MatEV
extractMode = 1; % Extract all data
modeArray = 1;
[TimeStamps, EventStrings, Header] = Nlx2MatEV(file, fieldSelection, extractHeader,extractMode, modeArray );
TimeStamps = TimeStamps ./ 1000000;
open EventStrings
%% match timestamps with each event

% settings: NEED to check!!!
% Ind_begin = [4,5;8,9;12,13];
% Ind_sleep = [2 3;6,7;10,11;14,15];
Ind_begin = [4,5;8,9;12,13;16,17];
Ind_sleep = [2,3;6,7;10,11;14,15;18,19];
% Ind_begin = [4,5;8,9;12,13;16,17]+1;
% Ind_sleep = [2,3;7,8;11,12;15,16;19,20];
% Ind_begin = [2,3;6,7];
% Ind_sleep = [4,5;8,9];

duration_begin = 20*60;  % 20mins
duration_sleep = 10*60;  % 10mins

% get time stamps for each session
Ts_begin = TimeStamps(Ind_begin);
Ts_sleep = TimeStamps(Ind_sleep);

% fix the time stamps based on the duration
dt = diff(Ts_begin');
ind = find(dt>duration_begin);
Ts_begin(ind,2) = Ts_begin(ind,1)+duration_begin;
% rat016 20220112-ABBA
% Ts_begin(2,1) = Ts_begin(2,2)-duration_begin;
% Ts_begin(ind,2) = Ts_begin(ind,1)+duration_begin;

dt = diff(Ts_sleep');
ind = find(dt>duration_sleep);
Ts_sleep(ind,2) = Ts_sleep(ind,1)+duration_sleep;

save('Data_eventTS.mat','Ts_begin','Ts_sleep','duration_begin','duration_sleep');



