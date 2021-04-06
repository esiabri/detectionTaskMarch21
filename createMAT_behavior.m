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
hitCountReinforceStatic = [];
stimCountReinforceStatic = [];
hitCountDetectionStatic = [];
stimCountDetectionStatic = [];
hitCountReinforceDrifting = [];
stimCountReinforceDrifting = [];
hitCountDetectionDrifting = [];
stimCountDetectionDrifting = [];
hitCountContrastLearningStatic = [];
stimCountContrastLearningStatic = [];
hitCountContrastLearningDrifting = [];
stimCountContrastLearningDrifting = [];
firstLickDistReinforceStatic = {};
firstLickDistDetectionStatic = {};
firstLickDistReinforceDrifting = {};
firstLickDistDetectionDrifting = {};
firstLickDistContrastLearningStatic = {};
firstLickDistContrastLearningDrifting = {};
totalTodayReward = [];
FA_countReinforceStatic = [];
FA_countDetectionStatic = [];
FA_countReinforceDrifting = [];
FA_countDetectionDrifting = [];
extendedStimCountReinforceStatic = [];
extendedStimCountDetectionStatic = [];
extendedStimCountReinforceDrifting = [];
extendedStimCountDetectionDrifting = [];


save(dataFileAdd,'Date','weightYesterday','addWaterYesterday','trainingStage','hitCountReinforceStatic',...
    'stimCountReinforceStatic','hitCountDetectionStatic','stimCountDetectionStatic',...
    'hitCountReinforceDrifting','stimCountReinforceDrifting','hitCountDetectionDrifting','stimCountDetectionDrifting',...
    'hitCountContrastLearningStatic','stimCountContrastLearningStatic','hitCountContrastLearningDrifting','stimCountContrastLearningDrifting',...
    'firstLickDistReinforceStatic',...
    'firstLickDistDetectionStatic','firstLickDistReinforceDrifting','firstLickDistDetectionDrifting',...
    'firstLickDistContrastLearningStatic','firstLickDistContrastLearningDrifting','totalTodayReward',...
    'mouseInitialWeight','FA_countReinforceStatic','FA_countDetectionStatic','FA_countReinforceDrifting',...
    'FA_countDetectionDrifting','extendedStimCountReinforceStatic','extendedStimCountDetectionStatic',...
    'extendedStimCountReinforceDrifting','extendedStimCountDetectionDrifting');

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