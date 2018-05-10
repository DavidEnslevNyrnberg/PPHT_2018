% TODO
% - use libsvm as SVM method
% - problem if 15s window have 0 or 1 peak in get_HRV
% - Use real tags instead of "stupid"
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
for path = 1:length(testDir)
    if path == 1 || path == 2
        % do nothing for folder '.' and '..'
    else
        % indexing for TestSubject
        pas = path-2;
        % initialize dirs for each file.
        BVPdir = fullfile(dDir,testDir(path).name,'BVP.csv');
        BVPraw = load(BVPdir);
        ACCdir = fullfile(dDir,testDir(path).name,'ACC.csv');
        ACCraw = load(ACCdir);
        EDAdir = fullfile(dDir,testDir(path).name,'EDA.csv');
        EDAraw = load(EDAdir);
        TAGSdir = fullfile(dDir,testDir(path).name,'tags.csv');
        TAGSraw = load(TAGSdir);
        
        % Load meta data; test ID, initial Time, time instances for tags, normalized
        % value of the time instances for tags.
        TestSubject{pas}.ID = str2num(testDir(path).name);
        % if ID is hardcodded for missing tag. Find row number and store
        if ismember(TestSubject{pas}.ID,[460540, 459273, 463522])==1
            missTag3=[missTag3,pas];
            % tagCalc is tags calculated from stroop-start and-end for each
            % testsubject
            TestSubject{pas}.meta.tagCalc = [TAGSraw(2)-4.5*60; TAGSraw(2); TAGSraw(2)+7*60; TAGSraw(2)+11.5*60];
        elseif ismember(TestSubject{pas}.ID,[459892])==1
            missTag2=[missTag2,pas];
            TestSubject{pas}.meta.tagCalc = [TAGSraw(3)-11.5*60 ;TAGSraw(3)-7*60; TAGSraw(3) ;TAGSraw(3)+4.5*60];
        else
            TestSubject{pas}.meta.tagCalc = [TAGSraw(2)-4.5*60; TAGSraw(2); TAGSraw(3); TAGSraw(3)+4.5*60];
        end
        TestSubject{pas}.meta.iniTime = BVPraw(1);
        TestSubject{pas}.meta.tags = [BVPraw(1);TAGSraw];
        TestSubject{pas}.meta.tagStupid = [BVPraw(1),BVPraw(1)+5*60,BVPraw(1)+12*60,BVPraw(1)+17*60];
        if ismember(TestSubject{pas}.ID,[460857])==1 % hard code end data point for sub 17min test
            TestSubject{pas}.meta.tagStupid(4) = BVPraw(1)+17*60-20;
        end
        
        % Load ACC, BVP and EDA data
        TestSubject{pas}.ACC.fs = ACCraw(2,1);
        TestSubject{pas}.ACC.data = ACCraw(3:end,:);
        TestSubject{pas}.BVP.fs = BVPraw(2);
        TestSubject{pas}.BVP.data = BVPraw(3:end);
        TestSubject{pas}.EDA.fs = EDAraw(2);
        TestSubject{pas}.EDA.data = EDAraw(3:end);
    end
end

%% calculate tag Matrix for signal windows
for pas = 1:length(TestSubject)
    % calc time axis for BVP
    T1 = TestSubject{pas}.meta.iniTime;
    tStep = 1/TestSubject{pas}.BVP.fs;
    T2 = (length(TestSubject{pas}.BVP.data)-1)*tStep+T1;
    TestSubject{pas}.BVP.timeBVPax = T1:tStep:T2;
    % find tag location for BVP
    [~,TestSubject{pas}.BVP.tagLocBVP] = maxk(ismember(TestSubject{pas}.BVP.timeBVPax,TestSubject{pas}.meta.tagStupid),4);
    TestSubject{pas}.BVP.tagBoolBVP = ismember(TestSubject{pas}.BVP.timeBVPax,TestSubject{pas}.meta.tagStupid);
    
    % calc time axis for EDA
    T1 = TestSubject{pas}.meta.iniTime;
    tStep = 1/TestSubject{pas}.EDA.fs;
    T2 = (length(TestSubject{pas}.EDA.data)-1)*tStep+T1;
    TestSubject{pas}.EDA.timeEDAax = T1:tStep:T2;
    % find tag location for BVP
    [~,TestSubject{pas}.EDA.tagLocEDA] = maxk(ismember(TestSubject{pas}.EDA.timeEDAax,TestSubject{pas}.meta.tagStupid),4);
    TestSubject{pas}.EDA.tagBoolEDA = ismember(TestSubject{pas}.EDA.timeEDAax,TestSubject{pas}.meta.tagStupid);
end

% plot the raw data with time tags.
if doPlot == 1
    pas = 19;
    close all
    figure
    subplot(2,1,1)
    timeAxEDA = datetime(TestSubject{pas}.EDA.timeEDAax,'ConvertFrom','posixtime','TimeZone','local');
    plot(timeAxEDA,TestSubject{pas}.EDA.data,'r', 'LineWidth',0.75);
    hold on;
    plot(timeAxEDA(TestSubject{pas}.EDA.tagBoolEDA),TestSubject{pas}.EDA.data(TestSubject{pas}.EDA.tagBoolEDA),'s', 'MarkerSize',4, 'LineWidth',4)
    hold off
    legend('EDA signal', 'tags for test transition', 'Location','northwest')
    title('Raw EDA signal with tag marks')
    xlabel('Time'); ylabel('EDA')
    subplot(2,1,2)
%     figure
    timeAxEDA = datetime(TestSubject{pas}.BVP.timeBVPax,'ConvertFrom','posixtime','TimeZone','local');
    plot(timeAxEDA,TestSubject{pas}.BVP.data,'r')
    hold on
    plot(timeAxEDA(TestSubject{pas}.BVP.tagBoolBVP),TestSubject{pas}.BVP.data(TestSubject{pas}.BVP.tagBoolBVP),'x','MarkerSize',5,'LineWidth',3)
    hold off
    legend('PPG signal', 'tags for test transition')
    title('Raw PPG signal with tag marks')
    xlabel('Time'); ylabel('PPG')
%     plot(TestSubject{ite}.BVP.timeBVPax,TestSubject{ite}.BVP.data,'g'); hold on; plot(TestSubject{ite}.BVP.timeBVPax(TestSubject{ite}.BVP.tagBoolBVP),TestSubject{ite}.BVP.data(TestSubject{ite}.BVP.tagBoolBVP),'x','MarkerSize',5,'LineWidth',3)
end

% datetime(TestSubject{ite}.meta.tags, 'TimeZone','local', 'ConvertFrom','posixtime')

%% EDA - peak count and slope
pas = 4; % remove later

dataEDA = TestSubject{pas}.EDA.data;
fsEDA = 2*TestSubject{pas}.EDA.fs; %Samplingsrate of EDA (Hz)

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
[valPeakEDA,locPeakEDA] = findpeaks(detrendDataFilterEDA); %locPeakEDA gives all peaks in the signal

length(locPeakEDA) 
timeEDA = 1:length(detrendDataFilterEDA);

close all
figure
% subplot(2,1,1)
plot(timeEDA,detrendDataFilterEDA,'g', timeEDA(locPeakEDA),valPeakEDA,'or')
% subplot(2,1,2)
% plot(timeBVP,filtDataBVP,'g',timeBVP(locPeakBVP),valPeakBVP,'or')

% figure;
% plot(detrendDataFilterEDA); title('Detrend signal')

%% BVP - HR calc and HRV feature extraction
timeCut = 30*64; % remove 30sec before and after 'windowBVP/EDA'.
f_resample = 8; % amount of interpolation steps for HRV_resample

for pas = 1:length(TestSubject)
    % filterdesigner > bandpass, FIR(equiripple), minimum order, Density
    % Factor(20), (hz, 64, 0, 0.5, 15, 15,5), (dB, 60, 1, 80)
    % load('FIRfilt.mat')
    bpFirFilt{pas} = designfilt('bandpassfir', ... % actually a LP-filter since
           'StopbandFrequency1',eps, 'PassbandFrequency1',.5, ... % StopbandFrequency1 = 0 = eps
           'PassbandFrequency2',15, 'StopbandFrequency2',15.5, ...
           'StopbandAttenuation1',60, 'PassbandRipple',1, 'StopbandAttenuation2',80, ...
           'DesignMethod','Equiripple', 'SampleRate',TestSubject{pas}.BVP.fs);
%     fvtool(bpFirFilt) %conform the designed filter
    % forward and backward filtering with the LP-filter of the BVP signal
    filtDataBVP{pas} = filtfilt(bpFirFilt{pas},TestSubject{pas}.BVP.data);
    
    % find the BVP peaks in the filtered PPG signal.
    [valPeakBVP{pas},locPeakBVP{pas}] = PPG2PEAK(filtDataBVP{pas}, 20, 0.75, 0);
    [valPeakBVPtest{pas},locPeakBVPtest{pas}] = findpeaks(-filtDataBVP{pas},'MinPeakHeight',20);
    % asigns the location of locPeakBVP to filtDataBVP's sequence size
    peakDataBVP{pas} = full(sparse(1,locPeakBVP{pas},1,1,length(filtDataBVP{pas})));
    
    if doPlot == 1
        % plot detected peaks of peaks for findpeaks and PPG2PEAK.
        timeAx = datetime(TestSubject{pas}.BVP.timeBVPax,'ConvertFrom','posixtime','TimeZone','local');
        figure(pas)
        subplot(2,1,1)
        plot(timeAx,filtDataBVP{pas},'g'); hold on
        plot(timeAx(locPeakBVPtest{pas}),-valPeakBVPtest{pas},'or', 'LineWidth',2)
        subplot(2,1,2)
        plot(timeAx,filtDataBVP{pas},'g'); hold on
        plot(timeAx(locPeakBVP{pas}),valPeakBVP{pas},'or', 'LineWidth',2)
    end
    
    % calculating windows with time offsets for - preStroop - Stroop - postStroop
    windowLoc1 = timeCut+TestSubject{pas}.BVP.tagLocBVP(1):TestSubject{pas}.BVP.tagLocBVP(2)-timeCut;
    windowLoc2 = timeCut+TestSubject{pas}.BVP.tagLocBVP(2):TestSubject{pas}.BVP.tagLocBVP(3)-timeCut;
    windowLoc3 = timeCut+TestSubject{pas}.BVP.tagLocBVP(3):TestSubject{pas}.BVP.tagLocBVP(4)-timeCut;
    % 15s segments for preStroop
    windowData1 = filtDataBVP{pas}(windowLoc1);
    windowPeak1 = peakDataBVP{pas}(windowLoc1);
    windowLengthStart1 = [1: TestSubject{pas}.BVP.fs*15: (length(windowData1))];
    windowLengthStart1(end) = [];
    windowLengthEnd1 = [1: TestSubject{pas}.BVP.fs*15: (length(windowData1))]-1;
    windowLengthEnd1(1) = [];
    for i = 1:length(windowLengthStart1)
        TestSubject{pas}.BVP.segmentWindowPeak1{i} = windowPeak1(windowLengthStart1(i):windowLengthEnd1(i));
    end
    % 15s segments for Stroop
    windowData2 = filtDataBVP{pas}(windowLoc2);
    windowPeak2 = peakDataBVP{pas}(windowLoc2);
    windowLengthStart2 = [1: TestSubject{pas}.BVP.fs*15: (length(windowData2))];
    windowLengthStart2(end) = [];
    windowLengthEnd2 = [1: TestSubject{pas}.BVP.fs*15: (length(windowData2))]-1;
    windowLengthEnd2(1) = [];
    for i = 1:length(windowLengthStart2)
        TestSubject{pas}.BVP.segmentWindowPeak2{i} = windowPeak2(windowLengthStart2(i):windowLengthEnd2(i));
    end
    % 15s segments for postStroop
    windowData3 = filtDataBVP{pas}(windowLoc3);
    windowPeak3 = peakDataBVP{pas}(windowLoc3);
    windowLengthStart3 = [1: TestSubject{pas}.BVP.fs*15: (length(windowData3))];
    windowLengthStart3(end) = [];
    windowLengthEnd3 = [1: TestSubject{pas}.BVP.fs*15: (length(windowData3))]-1;
    windowLengthEnd3(1) = [];
    for i = 1:length(windowLengthStart3)
        TestSubject{pas}.BVP.segmentWindowPeak3{i} = windowPeak3(windowLengthStart3(i):windowLengthEnd3(i));
    end
    
    if doPlot == 1
        figure
        plot(TestSubject{pas}.BVP.timeBVPax(windowLoc1), TestSubject{pas}.BVP.windowData1, 'b');
        hold on
        plot(TestSubject{pas}.BVP.timeBVPax(logical(TestSubject{pas}.BVP.windowPeak1)), TestSubject{pas}.BVP.windowData1(logical(TestSubject{pas}.BVP.windowPeak1)), 'ro', 'MarkerSize',4);
        hold off
    end
    for i = 1:length(windowLengthStart1)
        [HRV.pre{pas,i}, qrs_loc_time.pre{pas,i}, HRV_resample.pre{pas,i}, qrs_loc_time_resample.pre{pas,i}]=get_HRV(TestSubject{pas}.BVP.segmentWindowPeak1{i}, f_resample, TestSubject{pas}.BVP.fs, [pas,i]);
    end
    for i = 1:length(windowLengthStart2)
        [HRV.Stroop{pas,i}, qrs_loc_time.Stroop{pas,i}, HRV_resample.Stroop{pas,i}, qrs_loc_time_resample.Stroop{pas,i}]=get_HRV(TestSubject{pas}.BVP.segmentWindowPeak2{i}, f_resample, TestSubject{pas}.BVP.fs, [pas,i]);
    end
    for i = 1:length(windowLengthStart3)
        [HRV.post{pas,i}, qrs_loc_time.post{pas,i}, HRV_resample.post{pas,i}, qrs_loc_time_resample.post{pas,i}]=get_HRV(TestSubject{pas}.BVP.segmentWindowPeak3{i}, f_resample, TestSubject{pas}.BVP.fs, [pas,i]);
    end
%     [HRV{pas,2}, qrs_loc_time{pas,2}, HRV_resample{pas,2}, qrs_loc_time_resample{pas,2}]=get_HRV(TestSubject{pas}.BVP.windowPeak2, f_resample, TestSubject{pas}.BVP.fs);
%     [HRV{pas,3}, qrs_loc_time{pas,3}, HRV_resample{pas,3}, qrs_loc_time_resample{pas,3}]=get_HRV(TestSubject{pas}.BVP.windowPeak3, f_resample, TestSubject{pas}.BVP.fs);
end

%% HRV feature extraction
% NN intervals is normal peak to peak intervals
N = length(TestSubject);
f_resample=8; % declare interpolation resample rate, ???
msInt = 50/1000; % NN50 time interval margin

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
        freq_axis{i,j} = linspace(0,length(HRV_freq{i,j})/TestSubject{i}.BVP.fs/f_resample,length(HRV_freq{i,j}));
        HRV_power{i,j} = (abs(HRV_freq{i,j}).^2)./length(freq_axis{i,j});
        % total power
        TestSubject{i}.HRV_features(j).power = sum(HRV_power{i,j});
        % calculate VLF for freq <= 0.04Hz
        HRV_vlf{i,j} = HRV_power{i,j}(1:find(freq_axis{i,j}<=0.04,1,'last'));
        freq_vlf{i,j} = freq_axis{i,j}(find(freq_axis{i,j}<=0.04));
        TestSubject{i}.HRV_features(j).area_vlf = trapz(freq_vlf{i,j},HRV_vlf{i,j});
        % calculate LF for freq = ]0.04 : 0.15 Hz]
        HRV_LF{i,j} = HRV_power{i,j}(find(freq_axis{i,j}>0.04,1):find(freq_axis{i,j}<=0.15,1,'last'));
        freq_LF{i,j} = freq_axis{i,j}(find(freq_axis{i,j}>0.04,1):find(freq_axis{i,j}<=0.15,1,'last'));
        area_LF{i,j} = trapz(freq_LF{i,j},HRV_LF{i,j});
        TestSubject{i}.HRV_features(j).LF_norm = area_LF{i,j}/(TestSubject{i}.HRV_features(j).power-TestSubject{i}.HRV_features(j).area_vlf);
        % calculate HF for freq = ]0.15 : 0.4 Hz]
        HRV_HF{i,j} = HRV_power{i,j}(find(freq_axis{i,j}>0.15,1):find(freq_axis{i,j}<=0.4,1,'last'));
        freq_HF{i,j} = freq_axis{i,j}(find(freq_axis{i,j}>0.15,1):find(freq_axis{i,j}<=0.4,1,'last'));
        area_HF{i,j} = trapz(freq_HF{i,j},HRV_HF{i,j});
        TestSubject{i}.HRV_features(j).HF_norm = area_HF{i,j}/(TestSubject{i}.HRV_features(j).power-TestSubject{i}.HRV_features(j).area_vlf);
        % calculating LF/HF ratio
        TestSubject{i}.HRV_features(j).LF_HF_ratio = TestSubject{i}.HRV_features(j).LF_norm/TestSubject{i}.HRV_features(j).HF_norm;
    end
end

%% LibSVM lib 2-way classifier - 11 features
featureOneTest = cell2mat(struct2cell(TestSubject{1}.HRV_features(1)));

TestSubject{i}.HRV_features
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