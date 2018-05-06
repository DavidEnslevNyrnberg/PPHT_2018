% TODO
% - use libsvm as SVM method
% - SVM window 15 sec
% Use real tags instead of "stupid"
clear; clc; close all

doPlot = 1;
% selecting all test subject directories for either Windows or MAC
if strfind(computer,'PC')==1
    testDir = dir(fullfile('D:\31566_Personal_portable_health_technologies\7week\Data\')); % !!!!!~~CHANGE TO CORRECT PATH~~!!!!!
    dDir = testDir(1).folder;
elseif strfind(computer,'MAC')==1
    testDir = dir(fullfile('/Data/*')); % !!!!!~~CHANGE TO CORRECT PATH~~!!!!!
    dDir = pwd;
end
% github created 'nothing.txt' here we remove it as a path
for i = 1:length(testDir)
    if strcmp(testDir(i).name,'nothing.txt')==1
        testDir(i) = [];
    end
end

%% import all data to struct -> TestSubject
missTagSession = []; % initialize for hard coded test subjects with missing tags
for j = 1:length(testDir)
    if j == 1 || j == 2
        % do nothing for folder '.' and '..'
    else
        % indexing for TestSubject
        k = j-2;
        % initialize dirs for each file.
        BVPdir = fullfile(dDir,testDir(j).name,'BVP.csv');
        BVPraw = load(BVPdir);
        ACCdir = fullfile(dDir,testDir(j).name,'ACC.csv');
        ACCraw = load(ACCdir);
        EDAdir = fullfile(dDir,testDir(j).name,'EDA.csv');
        EDAraw = load(EDAdir);
        TAGSdir = fullfile(dDir,testDir(j).name,'tags.csv');
        TAGSraw = load(TAGSdir);
        
        % Load meta data; test ID, initial Time, time instances for tags, normalized
        % value of the time instances for tags.
        TestSubject{k}.ID = str2num(testDir(j).name);
        % if ID is hardcodded for missing tag. Find row number
        if ismember(TestSubject{k}.ID,[459892, 460540, 459273, 463522])==1
            missTagSession=[missTagSession,k];
        end
        TestSubject{k}.meta.iniTime = BVPraw(1);
        TestSubject{k}.meta.tags = [BVPraw(1);TAGSraw];
        TestSubject{k}.meta.tagStupid = [BVPraw(1),BVPraw(1)+5*60,BVPraw(1)+12*60,BVPraw(1)+17*60];
        
        
        % Load ACC, BVP and EDA data
        TestSubject{k}.ACC.fs = ACCraw(2,1);
        TestSubject{k}.ACC.data = ACCraw(3:end,:);
        TestSubject{k}.BVP.fs = BVPraw(2);
        TestSubject{k}.BVP.data = BVPraw(3:end);
        TestSubject{k}.EDA.fs = EDAraw(2);
        TestSubject{k}.EDA.data = EDAraw(3:end);
    end
end

%% calculate tag Matrix for signal windows
% disp({str2num(TestSubject{ite}.ID), str2num(TestSubject{ite}.meta.iniTime)})

for ite = 1:length(TestSubject)
%     clear T1 tStep T2
%     dataBVP = TestSubject{ite}.BVP.data;
%     fsBVP = TestSubject{ite}.BVP.fs;
    T1 = TestSubject{ite}.meta.iniTime;
    tStep = 1/TestSubject{ite}.BVP.fs;
    T2 = (length(TestSubject{ite}.BVP.data)-1)*tStep+T1;
    TestSubject{ite}.BVP.timeBVPax = T1:tStep:T2;
    [~,TestSubject{ite}.BVP.tagLocBVP] = maxk(ismember(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.meta.tagStupid),4);
    TestSubject{ite}.BVP.tagBoolBVP = ismember(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.meta.tagStupid);

%     dataEDA = TestSubject{ite}.EDA.data;
%     fsEDA = TestSubject{ite}.EDA.fs;
    T1 = TestSubject{ite}.meta.iniTime;
    tStep = 1/TestSubject{ite}.EDA.fs;
    T2 = (length(TestSubject{ite}.EDA.data)-1)*tStep+T1;
    TestSubject{ite}.EDA.timeEDAax = T1:tStep:T2;
    [~,TestSubject{ite}.EDA.tagLocEDA] = maxk(ismember(TestSubject{ite}.EDA.timeEDAax,TestSubject{ite}.meta.tagStupid),4);
    TestSubject{ite}.EDA.tagBoolEDA = ismember(TestSubject{ite}.EDA.timeEDAax,TestSubject{ite}.meta.tagStupid);
end

if doPlot == 1
    ite = 24;
    close all
    figure; 
    subplot(2,1,1)
    plot(TestSubject{ite}.EDA.timeEDAax,TestSubject{ite}.EDA.data,'r', 'LineWidth',0.75);
    hold on;
    plot(TestSubject{ite}.EDA.timeEDAax(TestSubject{ite}.EDA.tagBoolEDA),TestSubject{ite}.EDA.data(TestSubject{ite}.EDA.tagBoolEDA),'s', 'MarkerSize',4, 'LineWidth',4)
    hold off
    legend('EDA signal', 'tags for test transition', 'Location','southeast')
    title('Raw EDA signal with tag marks')
    xlabel('Time'); ylabel('EDA')
    subplot(2,1,2)
%     figure
    plot(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.BVP.data,'r')
    hold on
    plot(TestSubject{ite}.BVP.timeBVPax(TestSubject{ite}.BVP.tagBoolBVP),TestSubject{ite}.BVP.data(TestSubject{ite}.BVP.tagBoolBVP),'x','MarkerSize',5,'LineWidth',3)
    hold off
    legend('PPG signal', 'tags for test transition')
    title('Raw PPG signal with tag marks')
    xlabel('Time'); ylabel('PPG')
%     plot(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.BVP.data,'g'); hold on; plot(TestSubject{ite}.BVP.timeBVPax(TestSubject{ite}.BVP.tagBoolBVP),TestSubject{ite}.BVP.data(TestSubject{ite}.BVP.tagBoolBVP),'x','MarkerSize',5,'LineWidth',3)
end

% datetime(TestSubject{ite}.meta.tags,'TimeZone','local','ConvertFrom','posixtime')

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
timeCut = 30*64; % remove 30sec before and after 'windowBVP/EDA'.
f_resample = 8; % amount of interpolation steps for HRV_resample

for ite = 1:length(TestSubject)
%     dataBVP{ite} = TestSubject{ite}.BVP.data;
%     fsBVP{ite} = TestSubject{ite}.BVP.fs;
%     T1{ite} = TestSubject{ite}.meta.iniTime;
%     tStep{ite} = 1/fsBVP{ite};
%     T2{ite} = (length(dataBVP{ite})-1)*tStep{ite}+T1{ite};
%     timeBVPax{ite} = T1{ite}:tStep{ite}:T2{ite};
%     timeBVP2{ite} = 1:length(timeBVPax{ite});
%     windowBVP{ite} = ismember(timeBVPax{ite},TestSubject{ite}.meta.tagStupid);
    % filterdesigner > bandpass, FIR(equiripple), minimum order, Density
    % Factor(20), (hz, 64, 0, 0.5, 15, 15,5), (dB, 60, 1, 80)
    % load('FIRfilt.mat')
    bpFirFilt{ite} = designfilt('bandpassfir', ... % Matteo 0 doesnt work?
           'StopbandFrequency1',eps, 'PassbandFrequency1',.5, ...
           'PassbandFrequency2',15, 'StopbandFrequency2',15.5, ...
           'StopbandAttenuation1',60, 'PassbandRipple',1, 'StopbandAttenuation2',80, ...
           'DesignMethod','Equiripple', 'SampleRate',TestSubject{ite}.BVP.fs);
    % fvtool(bpFirFilt) %conform the designed filter
    filtDataBVP{ite} = filtfilt(bpFirFilt{ite},TestSubject{ite}.BVP.data);

    [valPeakBVP{ite},locPeakBVP{ite}] = PPG2PEAK(filtDataBVP{ite}, 20, 0.75, 0);
    [valPeakBVPtest{ite},locPeakBVPtest{ite}] = findpeaks(-filtDataBVP{ite},'MinPeakHeight',20);
    
    if doPlot == 1
        % plot peak result
        figure(ite)
        subplot(2,1,1)
        plot(TestSubject{ite}.BVP.timeBVPax,filtDataBVP{ite},'g'); hold on
        plot(TestSubject{ite}.BVP.timeBVPax(locPeakBVPtest{ite}),-valPeakBVPtest{ite},'or', 'LineWidth',2)
        subplot(2,1,2)
        plot(TestSubject{ite}.BVP.timeBVPax,filtDataBVP{ite},'g'); hold on
        plot(TestSubject{ite}.BVP.timeBVPax(locPeakBVP{ite}),valPeakBVP{ite},'or', 'LineWidth',2)
    end
    peakDataBVP{ite} = full(sparse(1,locPeakBVP{ite},1,1,length(filtDataBVP{ite})));
    
    % creating windows
    windowLoc1 = timeCut+TestSubject{ite}.BVP.tagLocBVP(1):TestSubject{ite}.BVP.tagLocBVP(2)-timeCut;
    windowLoc2 = timeCut+TestSubject{ite}.BVP.tagLocBVP(2):TestSubject{ite}.BVP.tagLocBVP(3)-timeCut;
    windowLoc3 = timeCut+TestSubject{ite}.BVP.tagLocBVP(3):TestSubject{ite}.BVP.tagLocBVP(4)-timeCut;
    
    TestSubject{ite}.BVP.windowLoc = [windowLoc1,windowLoc2,windowLoc3];
    TestSubject{ite}.BVP.windowData1 = filtDataBVP{ite}(windowLoc1)';
    TestSubject{ite}.BVP.windowData2 = filtDataBVP{ite}(windowLoc2)';
    TestSubject{ite}.BVP.windowData3 = filtDataBVP{ite}(windowLoc3)';
    TestSubject{ite}.BVP.windowPeak1 = peakDataBVP{ite}(windowLoc1);
    TestSubject{ite}.BVP.windowPeak2 = peakDataBVP{ite}(windowLoc2);
    TestSubject{ite}.BVP.windowPeak3 = peakDataBVP{ite}(windowLoc3);
    
    if doPlot == 1
        figure
        plot(TestSubject{ite}.BVP.timeBVPax(windowLoc1), TestSubject{ite}.BVP.windowData1, 'b');
        hold on
        plot(TestSubject{ite}.BVP.timeBVPax(logical(TestSubject{ite}.BVP.windowPeak1)), TestSubject{ite}.BVP.windowData1(logical(TestSubject{ite}.BVP.windowPeak1)), 'ro', 'MarkerSize',4);
        hold off
    end
    
    [HRV{ite,1}, qrs_loc_time{ite,1}, HRV_resample{ite,1}, qrs_loc_time_resample{ite,1}]=get_HRV(TestSubject{ite}.BVP.windowPeak1, f_resample, TestSubject{ite}.BVP.fs);
    [HRV{ite,2}, qrs_loc_time{ite,2}, HRV_resample{ite,2}, qrs_loc_time_resample{ite,2}]=get_HRV(TestSubject{ite}.BVP.windowPeak2, f_resample, TestSubject{ite}.BVP.fs);
    [HRV{ite,3}, qrs_loc_time{ite,3}, HRV_resample{ite,3}, qrs_loc_time_resample{ite,3}]=get_HRV(TestSubject{ite}.BVP.windowPeak3, f_resample, TestSubject{ite}.BVP.fs);
end

%% HRV feature extraction
N = length(TestSubject);
f_resample=8; % declare interpolation resample rate, ???
msInt = 50/1000; % NN50 time interval margin
count=1;

for i = 1:1
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
lengthArray =[1:(length(testDir)-2)*size(peakDataBVP,2)];
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

%15sec int idea
[sample_start: 15: sample_end]

[HRV, qrs_loc_time, HRV_resample, qrs_loc_time_resample]=get_HRV(peakDataBVP, f_resample, fsBVP);

% startT = TestSubject{1}.meta.iniTime;
% endT = (length(dataBVP)-1)*fsBVP+startT;
% time = [startT:fsBVP:endT];
[HRV, qrs_loc_time, HRV_resample, qrs_loc_time_resample]=get_HRV(peakDataBVP, f_resample, fsBVP);

% % time axis can be calculated by the .csv timeelement and 1/fs