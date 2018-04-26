% TODO
% - use libsvm as SVM method
% - SVM window 5 - 10 sec
% - REMEBER TO CHANGE pasDIR to 'real' data path
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
        TestSubject{k}.meta.iniTime = num2str(BVPraw(1));
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
ite = 1; % remove later

dataEDA = TestSubject{ite}.EDA.data;
fsEDA = 2*TestSubject{ite}.EDA.fs; %Samplingsrate of EDA (Hz)

 %% preprocessing
 %lowpass butterworth filter of order 6 on the EDA signal, 
 %[0.05 1.5](range energy of SCR) , cutoff frequency would be 1.5
 Wn = 1.5/(fsEDA/2);     %Cutoff frequency normalized
 [b,a] = butter(6, Wn,'low');    % 6th- order butterworth  band-pass filter
 
filter_EDA = filtfilt(b,a,dataEDA); 

figure;
subplot(2,1,1); plot(dataEDA); title('Pure signal')
subplot(2,1,2); plot(filter_EDA); title('Filtered')

%% Detection of peaks
detrendDataFilterEDA = detrend(filter_EDA);
[valPeakEDA,locPeakEDA] = findpeaks(-detrendDataFilterEDA); %locPeakEDA gives all peaks in the signal

length(locPeakEDA) 
timeEDA = 1:length(detrendDataFilterEDA);

close all
figure
% subplot(2,1,1)
plot(timeEDA,detrendDataFilterEDA,'g', timeEDA(locPeakEDA),-valPeakEDA,'or')
% subplot(2,1,2)
% plot(timeBVP,filtDataBVP,'g',timeBVP(locPeakBVP),valPeakBVP,'or')

% figure;
% plot(detrendDataFilterEDA); title('Detrend signal')

%% BVP - HR calc and HRV feature extraction
ite = 2; % remove later

dataBVP = TestSubject{ite}.BVP.data;
fsBVP = TestSubject{ite}.BVP.fs;
% filterdesigner > bandpass, FIR(equiripple), minimum order, Density
% Factor(20), (hz, 64, 0, 0.5, 15, 15,5), (dB, 60, 1, 80)
% load('FIRfilt.mat')

bpFirFilt = designfilt('bandpassfir', ... % Matteo 0 doesnt work?
       'StopbandFrequency1',eps, 'PassbandFrequency1',.5, ...
       'PassbandFrequency2',15, 'StopbandFrequency2',15.5, ...
       'StopbandAttenuation1',60, 'PassbandRipple',1, 'StopbandAttenuation2',80, ...
       'DesignMethod','Equiripple', 'SampleRate',fsBVP);
% fvtool(bpFirFilt) %conform the designed filter
filtDataBVP = filtfilt(bpFirFilt,dataBVP);
timeBVP = 1:length(filtDataBVP);

[valPeakBVP,locPeakBVP] = PPG2PEAK(filtDataBVP, 20, 0.75);
% [valPeakBVP,locPeakBVP] = findpeaks(-filtDataBVP,'MinPeakHeight',10);
% 
% close all
% figure
% plot(timeBVP,filtDataBVP,'g',timeBVP(locPeakBVP),-valPeakBVP,'or')
% 
% diffStats = quantile(diff(locPeakBVP),[0.25 0.50 0.75]);
% 
% figure; plot(diff(locPeakBVP),'o')
% 
% failedPeak = diff(locPeakBVP) < floor(diffStats(1)*0.75);
% locSuccesPeak = find(failedPeak ==0);
% locPeakBVP2 = locPeakBVP(locSuccesPeak);
% valPeakBVP2 = -valPeakBVP(locSuccesPeak);
close all;
[valPeakBVPtest,locPeakBVPtest] = findpeaks(-filtDataBVP,'MinPeakHeight',10);

% plot peak result
figure
subplot(2,1,1)
plot(timeBVP,filtDataBVP,'g',timeBVP(locPeakBVPtest),-valPeakBVPtest,'or')
subplot(2,1,2)
plot(timeBVP,filtDataBVP,'g',timeBVP(locPeakBVP),valPeakBVP,'or')
% startT = TestSubject{1}.meta.iniTime;
% endT = (length(dataBVP)-1)*fsBVP+startT;
% time = [startT:fsBVP:endT];

peakDataBVP = full(sparse(1,locPeakBVP,1,1,length(filtDataBVP)));

f_resample = 8;
[HRV, qrs_loc_time, HRV_resample, qrs_loc_time_resample]=get_HRV(peakDataBVP, f_resample, fsBVP);

%% SVM - 2-way classifier