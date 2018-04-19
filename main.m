% TODO
% - use libsvm as SVM method
% - SVM window 5 - 10 sec
% - REMEBER TO CHANGE pasDIR to 'real' data path
% direct translation between MAC and windows
clear; clc; close all

if strfind(computer,'PC')==1
    pasDir = dir(fullfile('..\Data\*')); % !!!!!~~CHANGE TO CORRECT PATH~~!!!!!
    dDir = pasDir(1).folder;
elseif strfind(computer,'MAC')==1
    pasDir = dir(fullfile('/Data/*')); % !!!!!~~CHANGE TO CORRECT PATH~~!!!!!
    dDir = pwd;
end

%% import all data to struct -> TestSubject
for j = 1:length(pasDir)
    if j == 1 || j == 2
        % do nothing for folder '.' and '..'
    else
        % indexing for TestSubject
        k = j-2;
        % initialize dirs for each file.
        BVPdir = fullfile(dDir,pasDir(j).name,'BVP.csv');
        BVPraw = load(BVPdir);
        ACCdir = fullfile(dDir,pasDir(j).name,'ACC.csv');
        ACCraw = load(ACCdir);
        EDAdir = fullfile(dDir,pasDir(j).name,'EDA.csv');
        EDAraw = load(EDAdir);
        TAGSdir = fullfile(dDir,pasDir(j).name,'tags.csv');
        TAGSraw = load(TAGSdir);
        
        % Load meta data; test ID, initial Time, time instances for tags, normalized
        % value of the time instances for tags.
        TestSubject{k}.ID = pasDir(j).name;
        TestSubject{k}.meta.iniTime = BVPraw(1);
        TestSubject{k}.meta.tags = TAGSraw;
        TestSubject{k}.meta.tagsNorm = TAGSraw-BVPraw(1);
        
        % Load ACC, BVP and EDA data
        TestSubject{k}.ACC.fs = ACCraw(2,1);
        TestSubject{k}.ACC.data = ACCraw(3:end,:);
        TestSubject{k}.BVP.fs = BVPraw(2);
        TestSubject{k}.BVP.data = BVPraw(3:end);
        TestSubject{k}.EDA.fs = EDAraw(2);
        TestSubject{k}.EDA.data = EDAraw(3:end);
    end
end

%% calculate windows for signals


%% EDA - peak count and slope

%% BVP - HR calc and HRV feature extraction
ite = 1; % remove later

dataBVP = TestSubject{ite}.BVP.data;
fsBVP = TestSubject{ite}.BVP.fs;
load('FIRfilt.mat')

d = designfilt('bandpass', 'PassbandFrequency',0.15,'StopbandFrequency',0.2, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple');

% butter_order = 4;
% butter_cut = [.5,15];
% contruct bandpass butterworth filter
filtfilt(filterCoe,dataBVP)
[b,a] = butter(butter_order,butter_cut/(2*fsBVP),'bandpass');

% filtfilt(filterCoe

[~,locPeakBVP] = findpeaks(-dataBVP);

% for i = 1:N
%     % array structure of data
%     ecg{i} = ecg_data(i,~isnan(ecg_data(i,:)));
%     % apply filter to data
%     ecg_prepros_filt{i} = filtfilt(b,a,ecg{i});
%     timeAx = 1:length(ecg_prepros_filt{i}); % data time axis
% end

close all
figure;
subplot(2,1,1); plot(dataBVP); title('Pure signal')
% subplot(2,1,2); plot(detrendDataBVP); title('Detrend signal')

%% failure

k = 1:fsBVP;

N = length(X);
M = nan(fsBVP,N);
for k = 1:fsBVP
%     W_k = 2*k;
    for i = k+2:N-k+1
        if detrendDataBVP(i-1)>detrendDataBVP(i-k-1) && detrendDataBVP(i-1)>detrendDataBVP(i+k-1)
            M(k,i) = 1;
        else
            M(k,i) = 0;
        end
    end
end

%% SVM - 2-way classifier