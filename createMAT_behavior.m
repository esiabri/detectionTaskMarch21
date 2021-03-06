% this m-file create the MAT file for storing the behavioral data for one
% animal

% Enter the mouse number through the pup-up window
prompt = {'Enter the mouse number:','Mouse Initial Weight:'}; 
titleBox = 'Input';
dims = [1 35]; 
dialogBoxInputs = inputdlg(prompt,titleBox,dims);

mouseNumber = dialogBoxInputs{1};
mouseInitialWeight = str2double(dialogBoxInputs{2});

% Select the folder for saving the file, this is the main root of the
% external drive on the behavior computer

defaultPath = 'D:\';
selpath = uigetdir(defaultPath);

dataFileAdd = strcat(selpath,'\',mouseNumber,'.mat');

Date = [];
weightYesterday = [];
addWaterYesterday = [];
trainingStage = [];
hitCount = {};
stimCount = {};
firstLickDist = {};
totalTodayReward = [];
FA_count = [];
extendedStimCount = [];


save(dataFileAdd,'Date','weightYesterday','addWaterYesterday','trainingStage','hitCount',...
    'stimCount','firstLickDist','totalTodayReward',...
    'mouseInitialWeight','FA_count','extendedStimCount');

% Date = [Date, string(date)];
% weightYesterday = [weightYesterday, 27];
% 
% save(dataFileAdd,'Date','weight','-append');




% clear all
% % load the saved file
% defaultPath = 'D:\DATA';
% [fileName,pathToFile] = uigetfile(defaultPath);
% dataFileAdd = strcat(pathToFile,fileName);
% 
% load(dataFileAdd)
% 
% newVar = [];
% save(dataFileAdd,'newVar','-append')