% TODO
% - use libsvm as SVM method
% - SVM window 5 - 10 sec
% - REMEBER TO CHANGE pasDIR to 'real' data path
clear; clc; close all

if strfind(computer,'PC')==1
    pasDir = dir(fullfile('D:\Data2\*')); % !!!!!~~CHANGE TO CORRECT PATH~~!!!!!
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
        TestSubject{k}.ID = str2num(pasDir(j).name);
        TestSubject{k}.meta.iniTime = BVPraw(1);
        TestSubject{k}.meta.tags = [BVPraw(1);TAGSraw];
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
ite = 4;
% disp({str2num(TestSubject{ite}.ID), str2num(TestSubject{ite}.meta.iniTime)})

for i = 1:length(TestSubject)

dataBVP{ite} = TestSubject{ite}.BVP.data;
fsBVP{ite} = TestSubject{ite}.BVP.fs;
T1{ite} = TestSubject{ite}.meta.iniTime;
tStep{ite} = 1/fsBVP{ite};
T2{ite} = (length(dataBVP{ite})-1)*tStep{ite}+T1{ite};
tAx{ite} = T1{ite}:tStep{ite}:T2{ite};
tags{ite} = TestSubject{ite}.meta.tags;
[~,locTagBVP{ite}] = min(abs(tAx{ite}-tags{ite}),[],2);

if TestSubject{ite}.ID==459892
    locTagBVP{ite}(4) = locTagBVP{ite}(3);
    locTagBVP{ite}(3) = locTagBVP{ite}(4)-5*60*fsBVP;
elseif TestSubject{ite}.ID==459273 || TestSubject{ite}.ID==460540 || TestSubject{ite}.ID==463522
    locTagBVP{ite}(3:4) = locTagBVP{ite}(2:3);
    locTagBVP{ite}(2) = 5*60*fsBVP{ite};
else 
    [~,locTagBVP{ite}] = min(abs(tAx{ite}-tags{ite}),[],2);
end

end

close all
figure;plot(tAx{ite},dataBVP{ite},'g')
hold on 
plot(tAx{ite}(locTagBVP{ite}),dataBVP{ite}(locTagBVP{ite})./dataBVP{ite}(locTagBVP{ite}),'x','MarkerSize',5,'LineWidth',3)

% startT = TestSubject{1}.meta.iniTime;
% endT = (length(dataBVP)-1)*fsBVP+startT;
% time = [startT:fsBVP:endT];

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
% ite = 2; % remove later
for ite = 1:length(pasDir)-2
dataBVP{ite} = TestSubject{ite}.BVP.data;
fsBVP{ite} = TestSubject{ite}.BVP.fs;
% filterdesigner > bandpass, FIR(equiripple), minimum order, Density
% Factor(20), (hz, 64, 0, 0.5, 15, 15,5), (dB, 60, 1, 80)
% load('FIRfilt.mat')

bpFirFilt{ite} = designfilt('bandpassfir', ... % Matteo 0 doesnt work?
       'StopbandFrequency1',eps, 'PassbandFrequency1',.5, ...
       'PassbandFrequency2',15, 'StopbandFrequency2',15.5, ...
       'StopbandAttenuation1',60, 'PassbandRipple',1, 'StopbandAttenuation2',80, ...
       'DesignMethod','Equiripple', 'SampleRate',fsBVP{ite});
% fvtool(bpFirFilt) %conform the designed filter
filtDataBVP{ite} = filtfilt(bpFirFilt{ite},dataBVP{ite});
timeBVP = 1:length(filtDataBVP{ite});

[valPeakBVP{ite},locPeakBVP{ite}] = PPG2PEAK(filtDataBVP{ite}, 20, 0.75);
[valPeakBVPtest{ite},locPeakBVPtest{ite}] = findpeaks(-filtDataBVP{ite},'MinPeakHeight',20);
% plot peak result
figure(ite)
subplot(2,1,1)
plot(timeBVP,filtDataBVP{ite},'g',timeBVP(locPeakBVPtest{ite}),-valPeakBVPtest{ite},'or', 'LineWidth',2)
subplot(2,1,2)
plot(timeBVP,filtDataBVP{ite},'g',timeBVP(locPeakBVP{ite}),valPeakBVP{ite},'or', 'LineWidth',2)


peakDataBVP{ite} = full(sparse(1,locPeakBVP{ite},1,1,length(filtDataBVP{ite})));

f_resample = 8;

[HRV{ite}, qrs_loc_time{ite}, HRV_resample{ite}, qrs_loc_time_resample{ite}]=get_HRV(peakDataBVP{ite}, f_resample, fsBVP{ite});
end

%% HRV feature extraction
N = (length(pasDir)-2);
f_resample=8; % declare interpolation resample rate, ???
msInt = 50/1000; % NN50 time interval margin
count=1;

for i = 1:N
% time feature extraction
HRV_feature(i).mean_RR = mean(HRV{i}); % mean of NN intervals
mean_RR(count) = mean(HRV{i});
HRV_feature(i).SD_NN = std(HRV{i}); % standard deviation of NN intervals [SDNN]
SD_NN(count) = std(HRV{i});
HRV_feature(i).RMS_NN = rms(diff(HRV{i})); % RMS of difference between adjacent NN intervals [RMSNN]
RMS_NN(count) = rms(diff(HRV{i}));
HRV_feature(i).SDSD = std(diff(HRV{i})); % STD of difference between adjacent NN intervals [SDSD]
SDSD(count) = std(diff(HRV{i}));
HRV_feature(i).NN50 = sum(diff(HRV{i}) > msInt); % NN intervals that differ by interval margin [NN50]
NN50(count) = sum(diff(HRV{i}) > msInt);
HRV_feature(i).pNN50 = HRV_feature(i).NN50/length(HRV{i})*100; % percentage of NN50 intervals in signal [pNN50]
pNN50(count) = HRV_feature(i).NN50/length(HRV{i})*100;

%  frequency feature extraction
HRV_freq{i} = fft(HRV_resample{i});
freq_axis{i} = length(HRV_resample{i});
%         P2{i,j} = abs(HRV_freq{i,j}/freq_axis{i,j}); % what is this Erik?
f{i} = f_resample*(0:(freq_axis{i}/2))/freq_axis{i}; % need better name?
%         HRV_fft{i,j} = P2{i,j}(1:ceil(freq_axis{i,j}/2)+1);
HRV_power{i} = (abs(HRV_freq{i}).^2)/freq_axis{i};
% total power
HRV_feature(i).power = sum(HRV_power{i});
power(count) = sum(HRV_power{i});

% calculating VLF area
HRV_vlf{i} = HRV_power{i}(1:find(f{i}<0.04,1,'last'));
f_vlf{i} = f{i}(find(f{i}<0.04));
HRV_feature(i).area_vlf = trapz(f_vlf{i},HRV_vlf{i});
area_vlf(count) = trapz(f_vlf{i},HRV_vlf{i});

% calculation LF area
HRV_lf{i} = HRV_power{i}(find(f{i}>0.05,1):find(f{i}<0.15,1,'last'));
f_lf{i} = f{i}(find(f{i}>0.05,1):find(f{i}<0.15,1,'last'));
area_lf{i} = trapz(f_lf{i},HRV_lf{i});
HRV_feature(i).LF_norm = area_lf{i}/(HRV_feature(i).power-HRV_feature(i).area_vlf);
LF_norm(count) = area_lf{i}/(HRV_feature(i).power-HRV_feature(i).area_vlf);

% calculating HF area
HRV_hf{i} = HRV_power{i}(find(f{i}>0.16,1):find(f{i}<0.40,1,'last'));
f_hf{i} = f{i}(find(f{i}>0.16,1):find(f{i}<0.40,1,'last'));
area_hf{i} = trapz(f_hf{i},HRV_hf{i});
HRV_feature(i).HF_norm = area_hf{i}/(HRV_feature(i).power-HRV_feature(i).area_vlf);
HF_norm(count) = area_hf{i}/(HRV_feature(i).power-HRV_feature(i).area_vlf);

% calculating LF/HF ratio
HRV_feature(i).LF_HF_ratio = HRV_feature(i).LF_norm/HRV_feature(i).HF_norm;
LF_HF_ratio(count) = HRV_feature(i).LF_norm/HRV_feature(i).HF_norm;
count=1+count;
end
%% SVM - 2-way classifier
SVM_feature = [mean_RR;SD_NN;RMS_NN;SDSD;NN50;pNN50;power;area_vlf;LF_norm;HF_norm;LF_HF_ratio]';
lengthArray =[1:(length(pasDir)-2)*size(peakDataBVP,2)];
SVM_tag = classes(lengthArray);
t = templateSVM('Standardize',1,'KernelFunction','rbf');
% t = templateSVM('Standardize',1);

SVMModel = fitcecoc(SVM_feature,SVM_tag,'Learners',t,'ObservationsIn','rows');
SVM_fit = resubPredict(SVMModel);
[SVM_groupLabelsCM,SVM_grpOrder] = confusionmat(SVM_tag, SVM_fit);
SVM_confusion = SVM_groupLabelsCM';
SVM_ResubErr = resubLoss(SVMModel);
SVM_xVal = crossval(SVMModel, 'kfold',5);
SVM_xValErr = kfoldLoss(SVM_xVal);

%% Ideas


[HRV, qrs_loc_time, HRV_resample, qrs_loc_time_resample]=get_HRV(peakDataBVP, f_resample, fsBVP);


% % time axis can be calculated by the .csv timeelement and 1/fs