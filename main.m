% TODO
% - use libsvm as SVM method
% - SVM window 15 sec
% Use real tags instead of "stupid"
clear; clc; close all

doPlot = 0;
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
missTag2 = []; % initialize for hard coded test subjects with missing tags
missTag3 = []; % its known that 3 miss the stroopStart tag and 1 the stroopEnd tag
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
        % if ID is hardcodded for missing tag. Find row number and store
        if ismember(TestSubject{k}.ID,[460540, 459273, 463522])==1
            missTag3=[missTag3,k];
            % tagCalc is tags calculated from stroop-start and-end for each
            % testsubject
            TestSubject{k}.meta.tagCalc = [TAGSraw(2)-4.5*60; TAGSraw(2); TAGSraw(2)+7*60; TAGSraw(2)+11.5*60];
        elseif ismember(TestSubject{k}.ID,[459892])==1
            missTag2=[missTag2,k];
            TestSubject{k}.meta.tagCalc = [TAGSraw(3)-11.5*60 ;TAGSraw(3)-7*60; TAGSraw(3) ;TAGSraw(3)+4.5*60];
        else
            TestSubject{k}.meta.tagCalc = [TAGSraw(2)-4.5*60; TAGSraw(2); TAGSraw(3); TAGSraw(3)+4.5*60];
        end
        TestSubject{k}.meta.iniTime = BVPraw(1);
        TestSubject{k}.meta.tags = [BVPraw(1);TAGSraw];
        TestSubject{k}.meta.tagStupid = [BVPraw(1),BVPraw(1)+5*60,BVPraw(1)+12*60,BVPraw(1)+17*60];
        if ismember(TestSubject{k}.ID,[460857])==1 % hard code end data point for sub 17min test
            TestSubject{k}.meta.tagStupid(4) = BVPraw(1)+17*60-20;
        end
        
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
for ite = 1:length(TestSubject)
    % calc time axis for BVP
    T1 = TestSubject{ite}.meta.iniTime;
    tStep = 1/TestSubject{ite}.BVP.fs;
    T2 = (length(TestSubject{ite}.BVP.data)-1)*tStep+T1;
    TestSubject{ite}.BVP.timeBVPax = T1:tStep:T2;
    % find tag location for BVP
    [~,TestSubject{ite}.BVP.tagLocBVP] = maxk(ismember(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.meta.tagStupid),4);
    TestSubject{ite}.BVP.tagBoolBVP = ismember(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.meta.tagStupid);
    
    % calc time axis for EDA
    T1 = TestSubject{ite}.meta.iniTime;
    tStep = 1/TestSubject{ite}.EDA.fs;
    T2 = (length(TestSubject{ite}.EDA.data)-1)*tStep+T1;
    TestSubject{ite}.EDA.timeEDAax = T1:tStep:T2;
    % find tag location for BVP
    [~,TestSubject{ite}.EDA.tagLocEDA] = maxk(ismember(TestSubject{ite}.EDA.timeEDAax,TestSubject{ite}.meta.tagStupid),4);
    TestSubject{ite}.EDA.tagBoolEDA = ismember(TestSubject{ite}.EDA.timeEDAax,TestSubject{ite}.meta.tagStupid);
end

% plot the raw data with time tags.
if doPlot == 1
    ite = 19;
    close all
    figure; 
    subplot(2,1,1)
    plot(TestSubject{ite}.EDA.timeEDAax,TestSubject{ite}.EDA.data,'r', 'LineWidth',0.75);
    hold on;
    plot(TestSubject{ite}.EDA.timeEDAax(TestSubject{ite}.EDA.tagBoolEDA),TestSubject{ite}.EDA.data(TestSubject{ite}.EDA.tagBoolEDA),'s', 'MarkerSize',4, 'LineWidth',4)
    hold off
    legend('EDA signal', 'tags for test transition', 'Location','northwest')
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

% datetime(TestSubject{ite}.meta.tags, 'TimeZone','local', 'ConvertFrom','posixtime')

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
    % filterdesigner > bandpass, FIR(equiripple), minimum order, Density
    % Factor(20), (hz, 64, 0, 0.5, 15, 15,5), (dB, 60, 1, 80)
    % load('FIRfilt.mat')
    bpFirFilt{ite} = designfilt('bandpassfir', ... % actually a LP-filter since
           'StopbandFrequency1',eps, 'PassbandFrequency1',.5, ... % StopbandFrequency1 = 0 = eps
           'PassbandFrequency2',15, 'StopbandFrequency2',15.5, ...
           'StopbandAttenuation1',60, 'PassbandRipple',1, 'StopbandAttenuation2',80, ...
           'DesignMethod','Equiripple', 'SampleRate',TestSubject{ite}.BVP.fs);
    % fvtool(bpFirFilt) %conform the designed filter
    % forward and backward filtering with the LP-filter of the BVP signal
    filtDataBVP{ite} = filtfilt(bpFirFilt{ite},TestSubject{ite}.BVP.data);
    
    % find the BVP peaks in the filtered PPG signal.
    [valPeakBVP{ite},locPeakBVP{ite}] = PPG2PEAK(filtDataBVP{ite}, 20, 0.75, 0);
    [valPeakBVPtest{ite},locPeakBVPtest{ite}] = findpeaks(-filtDataBVP{ite},'MinPeakHeight',20);
    % asigns the location of locPeakBVP to filtDataBVP's sequence size
    peakDataBVP{ite} = full(sparse(1,locPeakBVP{ite},1,1,length(filtDataBVP{ite})));
    
    if doPlot == 1
        % plot detected peaks of peaks for findpeaks and PPG2PEAK.
        figure(ite)
        subplot(2,1,1)
        plot(TestSubject{ite}.BVP.timeBVPax,filtDataBVP{ite},'g'); hold on
        plot(TestSubject{ite}.BVP.timeBVPax(locPeakBVPtest{ite}),-valPeakBVPtest{ite},'or', 'LineWidth',2)
        subplot(2,1,2)
        plot(TestSubject{ite}.BVP.timeBVPax,filtDataBVP{ite},'g'); hold on
        plot(TestSubject{ite}.BVP.timeBVPax(locPeakBVP{ite}),valPeakBVP{ite},'or', 'LineWidth',2)
    end
    
    % calculating windows with time offsets for - preStroop - Stroop - postStroop
    windowLoc1 = timeCut+TestSubject{ite}.BVP.tagLocBVP(1):TestSubject{ite}.BVP.tagLocBVP(2)-timeCut;
    windowLoc2 = timeCut+TestSubject{ite}.BVP.tagLocBVP(2):TestSubject{ite}.BVP.tagLocBVP(3)-timeCut;
    windowLoc3 = timeCut+TestSubject{ite}.BVP.tagLocBVP(3):TestSubject{ite}.BVP.tagLocBVP(4)-timeCut;
    
    
    TestSubject{ite}.BVP.windowLoc = [windowLoc1,windowLoc2,windowLoc3];
    % 15s segment for Window1
    TestSubject{ite}.BVP.windowData1 = filtDataBVP{ite}(windowLoc1);
    TestSubject{ite}.BVP.windowPeak1 = peakDataBVP{ite}(windowLoc1);
%     segmentWindow1 = [1: TestSubject{ite}.BVP.fs*15: (length(TestSubject{ite}.BVP.windowData1))]
    
    % 15s segment for Window2
    TestSubject{ite}.BVP.windowData2 = filtDataBVP{ite}(windowLoc2);
    TestSubject{ite}.BVP.windowData3 = filtDataBVP{ite}(windowLoc3);
    % 15s segment for Window3
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
% NN intervals is normal peak to peak intervals
N = length(TestSubject);
f_resample=8; % declare interpolation resample rate, ???
msInt = 50/1000; % NN50 time interval margin
count=1;

for i = 1:size(HRV,1)
    for j = 1:size(HRV,2)
        % Time features for HRV
        % Mean of NN intervals
        TestSubject{i}.HRV_features(j).mean = mean(HRV{i,j});
        % Standard deviation of NN intervals [SDNN]
        TestSubject{i}.HRV_features(j).sd = std(HRV{i,j});
        % RMS of difference between adjacent NN intervals [RMSNN]
        TestSubject{i}.HRV_features(j).RMS_NN = rms(diff(HRV{i,j}));
        % STD of difference between adjacent NN intervals [SDSD]
        TestSubject{i}.HRV_features(j).SDSD = std(diff(HRV{i,j}));
        % NN intervals that differ by interval margin [NN50]
        TestSubject{i}.HRV_features(j).NN50 = sum(diff(HRV{i}) > msInt);
        % Percentage of NN50 intervals in signal [pNN50]
        TestSubject{i}.HRV_features(j).pNN50 = TestSubject{i}.HRV_features(j).NN50/length(HRV{i})*100;
        
        % Frequency freature for HRV
        % Calculate frequency- and power-spectrum
        HRV_freq{i,j} = fft(HRV_resample{i,j});
        freq_axis{i,j} = linspace(0,length(HRV_freq{i,j})/TestSubject{i}.BVP.fs,length(HRV_freq{i,j}));
        HRV_power{i,j} = (abs(HRV_freq{i,j}).^2)/freq_axis{i,j};
        % total power
        TestSubject{i}.HRV_feature(i,j).power = sum(HRV_power{i,j});
        % calculate VLF for freq <= 0.04Hz
%         HRV_vlf{i,j} = HRV_power{i,j}(1:find(freq_axis{i,j}<=0.04,1,'last'));
%         freq_vlf{i,j} = freq_axis{i,j}(find(freq_axis{i,j}<=0.04));
%         TestSubject{i}.HRV_feature(i,j).area_vlf = trapz(freq_vlf{i,j},HRV_vlf{i,j});
%         calculate LF for freq = ]0.04 : 0.15 Hz]
%         HRV_LF{i,j} = HRV_power{i,j}(find(freq_axis{i,j}>0.04,1):find(freq_axis{i,j}<=0.15,1,'last'));
%         freq_LF{i,j} = freq_axis{i,j}(find(freq_axis{i,j}>0.04,1):find(freq_axis{i,j}<=0.15,1,'last'));
%         area_LF{i,j} = trapz(freq_LF{i,j},HRV_LF{i,j});
%         TestSubject{i}.HRV_feature(i,j).LF_norm = area_LF{i,j}/(TestSubject{i}.HRV_feature(i,j).power-HRV_feature(i,j).area_vlf);
%         calculate HF for freq = ]0.15 : 0.4 Hz]

    end
end

%%
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
freq_vlf{i} = f{i}(find(f{i}<0.04));
HRV_feature(i).area_vlf = trapz(freq_vlf{i},HRV_vlf{i});
area_vlf(count) = trapz(freq_vlf{i},HRV_vlf{i});

% calculation LF area
HRV_LF{i} = HRV_power{i}(find(f{i}>0.05,1):find(f{i}<0.15,1,'last'));
freq_LF{i} = f{i}(find(f{i}>0.05,1):find(f{i}<0.15,1,'last'));
area_LF{i} = trapz(freq_LF{i},HRV_LF{i});
HRV_feature(i).LF_norm = area_LF{i}/(HRV_feature(i).power-HRV_feature(i).area_vlf);
LF_norm(count) = area_LF{i}/(HRV_feature(i).power-HRV_feature(i).area_vlf);

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