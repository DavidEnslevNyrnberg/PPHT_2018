% TODO
% - use libsvm as SVM method
% - SVM window 5 - 10 sec
% - REMEBER TO CHANGE pasDIR to 'real' data path
clear; clc; close all

wDir = pwd;
pasDir = dir(fullfile('..\Data\*')); % !!!!!~~CHANGE TO CORRECT PATH~~!!!!!

%% import all data to struct -> TestSubject
for j = 1:length(pasDir)
    if j == 1 || j == 2
        % do nothing for folder '.' and '..'
    else
        % indexing for TestSubject
        k = j-2;
        % initialize dirs for each file.
        BVPdir = fullfile([pasDir(j).folder,'\',pasDir(j).name,'\BVP.csv']);
        BVPraw = load(BVPdir);
        ACCdir = fullfile([pasDir(j).folder,'\',pasDir(j).name,'\ACC.csv']);
        ACCraw = load(ACCdir);
        EDAdir = fullfile([pasDir(j).folder,'\',pasDir(j).name,'\EDA.csv']);
        EDAraw = load(EDAdir);
        TAGSdir = fullfile([pasDir(j).folder,'\',pasDir(j).name,'\tags.csv']);
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

%% SVM - 2-way classifier