clear all
close all 

stageOfTraining = 'Detection';

  
Screen('Preference', 'SkipSyncTests', 1);

% Detection Task. The final code for detection.
% GOAL = The mouse should lick the spout during the stimulus presentation to get the reward.
% STRUCTURE = The animal will be presented with visual stimulus trials for
% 1.5 seconds. For the trial to start, animal should not lick the spout for 5 seconds in half of the trials
% and for 6 seconds in the other half of the trials. If animal lick the
% spout during the last second of 6-second trials we count that trial as a
% FA

localDirectory = 'D:\recordedData\Behavioral\Ehsan';
googleDriveLocalFolderAddress = 'G:\My Drive\behavioralRestults';
networkDirectory = 'Z:\recordedData\Behavioral\Ehsan'; %this is the address that the copy will be attempted at the end of the
% of the session and also will be synced during the night with the hard
% drive (in case the network is not valid)

networkValid = 0;
if exist(networkDirectory, 'dir')
   networkValid = 1;   
end  
    
% entering session Information  
prompt = {'Mouse Number:','Mouse Weight (Yesterday):','Additional Water (Yesterday):','Box Number:','Mouse Session Number:','Reward Vol:','Trial Numbers (even):'}; 
titleBox = 'Input';
dims = [1 40]; 
dialogBoxInputs = inputdlg(prompt,titleBox,dims);

mouseNumber = dialogBoxInputs{1};
MouseWeight = dialogBoxInputs{2};
addWater = dialogBoxInputs{3};
boxNumber = dialogBoxInputs{4};
sessionNumber = dialogBoxInputs{5};
earnedRewardVol = str2num(dialogBoxInputs{6}); %Reward volume per trial
totalTrialNo = str2num(dialogBoxInputs{7}); %Total#of trials
 
% Making the folder for saving the data data folder name
dataFolderName = 'Mouse' + string(mouseNumber) + '_' + datestr(date,'mm-dd-yyyy') + '_Session' + sessionNumber + '_' + mfilename;

mkdir(localDirectory,dataFolderName);
dataFolderAdd = string(localDirectory) + '\' + dataFolderName;

slackAddFile = 'D:\webhookAdd.txt';
fileID = fopen(slackAddFile,'r');
slackNotifAddEhsan = fscanf(fileID,'%c');
fclose(fileID);
% slackNotifAddEhsan = 'https://hooks.slack.com/services/T8R4YBRS8/B01SHM907M2/WCxxVFN52Q7UNt9D8X5SqgU';

%% Initialization of the required daq card sessions
niDevName = 'Dev1';

%Analog input session for recording the following signals:
signalsRecordingSession = daq.createSession('ni');  
%1- output of the right lick sensor AI1
%2- copy of the step motor command AI2
%3- output to the speaker AI3
%4- trial tags (sent through the daq card) AI5
%5- photodiode signal (sensed on the screens) AI0
%6- copy of the sound sent to the speaker AI4 

addAnalogInputChannel(signalsRecordingSession,niDevName,0,'Voltage');
addAnalogInputChannel(signalsRecordingSession,niDevName,1,'Voltage');
addAnalogInputChannel(signalsRecordingSession,niDevName,2,'Voltage');
% addAnalogInputChannel(signalsRecordingSession,niDevName,3,'Voltage');
% addAnalogInputChannel(signalsRecordingSession,niDevName,4,'Voltage');
addAnalogInputChannel(signalsRecordingSession,niDevName,5,'Voltage');
% addAnalogInputChannel(signalsRecordingSession,niDevName,6,'Voltage');
% addAnalogInputChannel(signalsRecordingSession,niDevName,7,'Voltage'); 
% addAnalogInputChannel(signalsRecordingSession,niDevName,15,'Voltage');

%TAGGING TRIALS. Digital output session for tagging the trials

trialDigitalTagSession = daq.createSession('ni');

stimTagPortLine = 'port0/line2';
addDigitalChannel(trialDigitalTagSession,niDevName,stimTagPortLine,'OutputOnly');

trialDigitalTagSession.outputSingleScan(0);

% MONITORING LICK SENSOR. Digital input session for monitoring the spout
spoutSession = daq.createSession('ni');

sensorCopyPortLine1 = 'port0/line1';
addDigitalChannel(spoutSession,niDevName,sensorCopyPortLine1,'InputOnly');


% DELIVERING REWARD WITH STEP MOTOR. Digital Output session for right reward control
rewardStepMotorCtl1 = daq.createSession('ni');
rewardPortLine1 = 'port0/line0';

rewardStepMotorCtl1.addDigitalChannel(niDevName,rewardPortLine1,'OutputOnly');

rewardStepMotorEnable = daq.createSession('ni');
motorEnablePortLine = 'port0/line5';
rewardStepMotorEnable.addDigitalChannel(niDevName,motorEnablePortLine,'OutputOnly');
rewardStepMotorEnable.outputSingleScan(1);


%Configuring the session for recording analog inputs
signalsRecordingSession.Rate = 2e3;
for chNo=1:size(signalsRecordingSession.Channels,2)
    signalsRecordingSession.Channels(1,chNo).TerminalConfig = 'SingleEnded';
end
signalsRecordingSession.IsNotifyWhenDataAvailableExceedsAuto = 0;
signalsRecordingSession.IsContinuous = true;
inputSavingDur = 1; %based on the warning, 0.05 seconds is the minimum saving time that is possible, (higher interval to less affect the timings in the code!)
signalsRecordingSession.NotifyWhenDataAvailableExceeds = floor(signalsRecordingSession.Rate * inputSavingDur);

%recording data through the listener, we will also analyze the input data
%for detecting the licks by sending a copy of the lick sensor to a digital
%input
binFile = dataFolderAdd + '\' + 'synchedNI-CardInputs.bin';
fid1 = fopen(binFile,'w');
lh = signalsRecordingSession.addlistener('DataAvailable',@(src, event)logData(src, event, fid1));

%enable the pump
rewardStepMotorEnable.outputSingleScan(0);
%% ---------------------------Configuration of PTB--------------------------

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

%Screen('Preference', 'SkipSyncTests', 1);
% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
screens = Screen('Screens');


% selecting the screen for the stimulus presentation: mirroring two
% identical screens show them with one number here.
screenNumber = 2;%max(screens);%

window = Screen('OpenWindow', screenNumber);

% load C:\Users\Stimulation\Documents\MatlabScripts\AsusGammaTable23April2019SophiePhotometer
% Screen('LoadNormalizedGammaTable', window, gammaTable*[1 1 1]);


% Define black and white (white will be 1 and black 0). This is because
% in general luminace values are defined between 0 and 1 with 255 steps in
% between. All values in Psychtoolbox are defined between 0 and 1
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Do a simple calculation to calculate the luminance value for grey. This
% will be half the luminace values for white
gray=(white+black)/2;
photoDiodeGray1 = gray/5;  %during gray screen
photoDiodeGray2 = gray/2;  %during stimulus
% if round(gray)==white
%     gray=black;
% end

% Taking the absolute value of the difference between white and gray will
% help keep the grating consistent regardless of whether the CLUT color
% code for white is less or greater than the CLUT color code for black.
absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);

% Open an on screen window using PsychImaging and color it white. (or the
% default stimulus with white under the photo-diode area)

[scrWidthPix, scrHeightPix]=Screen('WindowSize',screenNumber);
% Measure the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', window);

% set the priority level of this thread as the highest
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);


%Defines the patch of screen sitting in front of the photodiaode in order
%ot send the change stim signal before downsampling anything
% patchRect = [(scrWidthPix - (scrWidthPix/25)) 0 scrWidthPix (scrHeightPix/15)];
% patchRect = [(scrWidthPix*2/3 - (scrWidthPix/3/25)) 0 scrWidthPix*2/3 (scrHeightPix/15)];
patchRect = [scrWidthPix-(scrWidthPix/10) 0 (scrWidthPix) (scrHeightPix/6)];


%--------------------
%% Visual Stim


StimDuration = 1.5; %in sec
StimFrames = round(StimDuration/ifi);

%Asus Screen Size: width (20.9235 inches, 53.1456 cm) height (11.7694
%inches, 29.8944 cm), pixel density: 91.76 pixels/inch
pixelDensityCM = (1280+720)/(15.41+9.05); %pixels/cm, this is approximate because the pixel density don't match for width and height

% Dimension of the region where will draw the Gabor in pixels
% Dimension of the region where will draw the Gabor in pixels
gaborDimCM_Height = 9.5; 
gaborDimPixHeight = floor(gaborDimCM_Height*pixelDensityCM); %windowRect(4) / 2;

gaborDimCM_Width = 14; 
gaborDimPixWidth = floor(gaborDimCM_Width*pixelDensityCM); %windowRect(4) / 2;

% Sigma of Gaussian
gaussianSigma = 0; %5
sigma = gaborDimCM_Height / gaussianSigma;

% Obvious Parameters
orientationPreferred = 0;
orientationNonPreferred = 90;
contrast = 100; %100
contrastMock = 0; %0
aspectRatio = 2.0;
phase = 0;

temporalFreq = 2;
degPerSec = 360 * temporalFreq;
degPerFrame =  degPerSec * ifi;
phase = 0;

% Numer of frames to wait before re-drawing
waitframes = 1;

% Spatial Frequency (Cycles Per Pixel)
% One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
spatialFrequency = 0.3; %cycles/cm

%freq = numCycles / gaborDimPix;
freq = spatialFrequency / pixelDensityCM; %cycles/pixel



% Build a procedural gabor texture (Note: to get a "standard" Gabor patch
% we set a grey background offset, disable normalisation, and set a
% pre-contrast multiplier of 0.5.
% For full details see:
% https://groups.yahoo.com/neo/groups/psychtoolbox/conversations/topics/9174
backgroundOffset = [0.5 0.5 0.5 0.0];
disableNorm = 1;
preContrastMultiplier = 0.5;
gabortex = CreateProceduralGabor(window, gaborDimPixWidth, gaborDimPixHeight, [],...
    backgroundOffset, disableNorm, preContrastMultiplier);

% gabortexMirror = CreateProceduralGabor(window, gaborDimPix, gaborDimPix, [],...
%     backgroundOffset, disableNorm, preContrastMultiplier);

% Randomise the phase of the Gabors and make a properties matrix.
propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];
propertiesMatmock = [phase, freq, sigma, contrastMock, aspectRatio, 0, 0, 0];

%--------------------
heightOffsetInCM = 0;
heightOffset = floor(heightOffsetInCM*pixelDensityCM);
widthOffset = 0;

stimHeightOffset = [0,heightOffset,0,gaborDimPixHeight+heightOffset];
righImageHorzPos = [0,0,gaborDimPixWidth,0];

%% ---------------------Camera and NI-Card Recording 

% Black is the default value under the photodiode and on the screen, Flip the default screen
% with White under the photo-diode
Screen('FillRect', window, gray);
Screen('FillRect',window, black, patchRect);
Screen('Flip', window);


signalsRecordingSession.startBackground();
disp('Start recording...')

Screen('FillRect', window, gray);
Screen('FillRect',window, black, patchRect);
Screen('Flip', window);


%--------------------------------------------------------------------------
%TRIAL PARAMETERS
afterStimGrayTime = 1; %in sec
afterStimGrayFrames = round(afterStimGrayTime/ifi);

 
stimRewardDelay = 0; %delay between the stim presentation and sensor monitoring start time

% Reward Volume:
rewardVol = earnedRewardVol; %in microL
syringeVol = 5;


%COUNTERS
StimCounter = 0;
hitCounter = 0;
missedCounter = 0;
FA_counter = 0;
CR_counter = 0;
extendedITT_trialsCounter = 0;
waitPeriodLickCounter = 0;
waitPeriodTotalTime = 0;
trialWaitPeriodLickCounter = 0;


earnedRewardVolTotal = 0;

%Vectors
allPhases = [];
stimPresTime = [];
TrialLick =[];
allRandNoLickTimes = [];
allWaitPeriodLickNo = [];

rewardCounter = 1;
rewardCompRate = 1.2;

minLickContact = 0.05;

noLickDurBeforeStim = 2; % in sec

%% stimVector
% here we want to build a stim vector with 5 or 6 s ITI trials
% 5s ITI: normal trials, encoded with 1
% 6s ITI: extended trials, encoded with 0
normalTrials = ones(1,totalTrialNo/2);
extendedTrials = zeros(1,totalTrialNo/2);

stimVector = [normalTrials extendedTrials];
stimVector = stimVector(randperm(length(stimVector)));

%% Free reward at start
deliverRewardMarch21(earnedRewardVol,syringeVol,rewardStepMotorCtl1); 

%% Session Trials

disp('');
disp('press any key to start the session');
KbWait;
startRecTime = GetSecs();

for trialNo=1:totalTrialNo
% Water will be delivered if the mice stop licking for 10 consecutive trials    
    
    
%     phase = 360*rand;
    phase = 0;
    allPhases = [allPhases phase];
    propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];

    % ESC session if necessary 
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
    if (find(keyCode) == 27)  %Press escape to manually finish the loop
        break;
    end
    
    
       
      
    if (find(keyCode) == 82)  %Press r to show the stim and give a reward
            Screen('FillRect', window, gray);
            Screen('FillRect',window, white, patchRect);
            Screen('DrawTextures', window, gabortex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
            kPsychDontDoRotation, propertiesMat');
        
        
            trialDigitalTagSession.outputSingleScan(1);
            [vblStim StimulusOnsetTime FlipTimestampStim MissedStim BeamposStim] = Screen('Flip', window, cuePresTime + (1 - 0.5) * ifi);

            vbl = vblStim;
            rewardDelivered = 0;

            for driftFrames=1:StimFrames-1

                % calculating the phase of the grating in this frame based on the
                % temporal frequency
                phase = phase + degPerFrame;
                propertiesMat(:, 1) = phase';

                % loading the next frame
                Screen('FillRect', window, gray);
                Screen('FillRect',window, white, patchRect);
                Screen('DrawTextures', window, gabortex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
                        kPsychDontDoRotation, propertiesMat');

                % Flip to the screen
                vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                if rewardDelivered == 0
                    if (GetSecs - vblStim) > stimRewardDelay
                        deliverRewardMarch21(earnedRewardVol,syringeVol,rewardStepMotorCtl1);
                        rewardDelivered = 1;
                    end 
                end


            end
            
            Screen('FillRect', window, gray);
            Screen('FillRect',window, black, patchRect);
            vblAfterStimGrayTime = Screen('Flip', window, vblStim + (StimFrames - 0.5) * ifi);
            trialDigitalTagSession.outputSingleScan(0);
            
            StimCounter = StimCounter + 1;
        
    end
    
    % Present gray screen
    Screen('FillRect', window, gray);
    Screen('FillRect',window, white, patchRect);
    
                
%no lick time is 5 seconds at least in all trials 
    noLickDurBeforeStim = 5;
    extendedNoLickDurBeforeStim = 1;

    % checking if this is an extended ITI trials for another second or not
    if stimVector(trialNo) == 1
        normalLickDurBeforeStim = 1;
    else
        normalLickDurBeforeStim = 0;
    end
% Check that the animal is not licking

    FA_flag = 0;

    normalWaitTime = 0;
    waitStartTime = GetSecs;
    while (1)
        scanStartTime = GetSecs;
        [lickFlag, relDetectionTime] = detectLickMarch21(noLickDurBeforeStim, spoutSession);
        if ~lickFlag
            cuePresTime = scanStartTime + relDetectionTime;
            waitPeriodTotalTime = waitPeriodTotalTime + (GetSecs - waitStartTime);
            normalWaitTime = 1;
        else
            % to correct for the multiple transitions in the output of the
            % lick sensor and counting the real number of the licks, after
            % lick detection we wait here until the lick sensor go back to
            % zero and then if the duration of the lick sensor output being
            % 1 is more than minLickContact, we count this transition in
            % the output of the lick sensor as a lick
            while (inputSingleScan(spoutSession))
                ;
            end
            if (GetSecs - (scanStartTime + relDetectionTime)) > minLickContact
                trialWaitPeriodLickCounter = trialWaitPeriodLickCounter + 1;
            end
        end
        if normalWaitTime
            if normalLickDurBeforeStim
                break;
            else
                scanStartTime = GetSecs;
                [lickFlag, relDetectionTime] = detectLickMarch21(extendedNoLickDurBeforeStim, spoutSession);
                if ~lickFlag
                    cuePresTime = scanStartTime + relDetectionTime;
                    waitPeriodTotalTime = waitPeriodTotalTime + (GetSecs - waitStartTime);
                    break;
                else
                    normalWaitTime = 0;
                    FA_flag = 1;
                end
            end
        end
    end
        

    StimCounter = StimCounter + 1;
        
        
% this is the real stimulus
    Screen('DrawTextures', window, gabortex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
        kPsychDontDoRotation, propertiesMat');
% 
%         Screen('DrawTextures', window, gabortex, [], righMirrorImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
%             kPsychDontDoRotation, propertiesMat');
    trialDigitalTagSession.outputSingleScan(1);
    [vblStim StimulusOnsetTime FlipTimestampStim MissedStim BeamposStim] = Screen('Flip', window, cuePresTime + (1 - 0.5) * ifi);

    % wait 0.5 seconds to give the reward
    
    vbl = vblStim;
    
    % In each iteration of this loop, next frame of this drifting stimulus
    % is presented
    lickCheck = 1;
    for driftFrames=1:StimFrames-1

        % calculating the phase of the grating in this frame based on the
        % temporal frequency
        phase = phase + degPerFrame;
        propertiesMat(:, 1) = phase';

        % loading the next frame
        Screen('FillRect', window, gray);
        Screen('FillRect',window, white, patchRect);
        Screen('DrawTextures', window, gabortex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
                kPsychDontDoRotation, propertiesMat');

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        if lickCheck == 1
            [lickFlag, relDetectionTime] = detectLickMarch21(waitframes*ifi, spoutSession);
            if lickFlag == 1
                deliverRewardMarch21(earnedRewardVol,syringeVol,rewardStepMotorCtl1);
                lickCheck = 0;
            end
        end

    end
    
    Screen('FillRect', window, gray);
    Screen('FillRect',window, black, patchRect);
    vblAfterStimGrayTime = Screen('Flip', window, vblStim + (StimFrames - 0.5) * ifi);
    trialDigitalTagSession.outputSingleScan(0);
 
    

%     trialDigitalTagSession.outputSingleScan(0);
 
            
    if lickFlag == 1 %if the MOUSE LICKS
        %Change Counters
        hitCounter = hitCounter + 1;
        rewardedTrial = 1;

%         rewardStepMotorEnable.outputSingleScan(0);
%         deliverReward(earnedRewardVol,syringeVol,rewardStepMotorCtl1,rewardStepMotorEnable);
        earnedRewardVolTotal = earnedRewardVolTotal + earnedRewardVol;
        deliverRewardFlag = 1;

        disp(['Hit: ', num2str(hitCounter), ' / ', num2str(StimCounter)]);
        disp(['Out of Stim Licks:', num2str(trialWaitPeriodLickCounter)]);
        disp(['passed time: ',num2str(floor((GetSecs()-startRecTime)/60)), ' Minutes']);

%         vblAfterStimGrayTime = Screen('Flip', window, vblStim + (StimFrames - 0.5) * ifi);
%         trialDigitalTagSession.outputSingleScan(0);


    else %If the MOUSE DOESN'T LICK
        %Change Counters
        missedCounter = missedCounter + 1;
        rewardedTrial = 0;

        disp(['Missed: ', num2str(missedCounter), ' / ', num2str(StimCounter)]);
        disp(['Out of Stim Licks:', num2str(trialWaitPeriodLickCounter)]);
        disp(['passed time: ',num2str(floor((GetSecs()-startRecTime)/60)), ' Minutes']);

%         vblAfterStimGrayTime = Screen('Flip', window, vblStim + (StimFrames - 0.5) * ifi);
%         trialDigitalTagSession.outputSingleScan(0);


    end
    
    if normalLickDurBeforeStim == 0
        extendedITT_trialsCounter = extendedITT_trialsCounter + 1;
        if FA_flag 
            FA_counter = FA_counter + 1;
            disp(' ');
            disp(['FA: ', num2str(FA_counter), ' / ', num2str(extendedITT_trialsCounter)]);
            
        else
            CR_counter = CR_counter + 1;
            disp(' ');
            disp(['CR: ', num2str(CR_counter), ' / ', num2str(extendedITT_trialsCounter)]);
        end
        
        
    end
        
    disp(' ');
    
    allPhases = [allPhases, phase];
    stimPresTime = [stimPresTime vblStim];
    allRandNoLickTimes = [allRandNoLickTimes noLickDurBeforeStim];
    allWaitPeriodLickNo = [allWaitPeriodLickNo trialWaitPeriodLickCounter];

    trialWaitPeriodLickCounter = 0;
    % ITI wait and counting the number of the licks
    ITI_lickSensorScanIntervals = 0.2;
    while((GetSecs - vblAfterStimGrayTime) < afterStimGrayTime)
        scanStartTime = GetSecs;
        [lickFlag, relDetectionTime] = detectLickMarch21(ITI_lickSensorScanIntervals, spoutSession);
        if lickFlag
            % to correct for the multiple transitions in the output of the
            % lick sensor and counting the real number of the licks, after
            % lick detection we wait here until the lick sensor go back to
            % zero and then if the duration of the lick sensor output being
            % 1 is more than minLickContact, we count this transition in
            % the output of the lick sensor as a lick
            while (inputSingleScan(spoutSession))
                ;
            end
            if (GetSecs - (scanStartTime + relDetectionTime)) > minLickContact
                trialWaitPeriodLickCounter = trialWaitPeriodLickCounter + 1;
            end
        end
    end
    
    waitPeriodTotalTime = waitPeriodTotalTime + afterStimGrayTime;
    
    waitPeriodLickCounter = waitPeriodLickCounter + trialWaitPeriodLickCounter;
    
    
end  
%%
FinalTime = num2str(floor((GetSecs()-startRecTime)/60));
sessionEndTime = now;

% notification to the experimenter to enable trigger in camera setting to stop frame recording 

SendSlackNotification(slackNotifAddEhsan,strcat('box ',num2str(boxNumber),' is done!'));


%Stop the recording
disp('Saving the session...')

%disable the pump

rewardStepMotorEnable.outputSingleScan(1);
pause(inputSavingDur) %To be sure that the whole session is recorded!

signalsRecordingSession.stop()

sca;
% Close the audio device
% PsychPortAudio('Close', pahandle);

delete(lh);
fclose(fid1);



%% Saving the variables of the session and the code file in the recording
%directory
save(dataFolderAdd + '\' + 'workspaceVariables');
copyfile(string(mfilename('fullpath')) + '.m', dataFolderAdd);

disp(' ');
disp(['Hit: ', num2str(hitCounter), ' / ', num2str(StimCounter)]);
disp(['Miss: ', num2str(missedCounter), ' / ', num2str(StimCounter)]);
disp(['FA: ', num2str(FA_counter), ' / ', num2str(extendedITT_trialsCounter)]);
disp(['CR: ', num2str(CR_counter), ' / ', num2str(extendedITT_trialsCounter)]);
disp(' ');
disp(['Total Reward: ', num2str(earnedRewardVolTotal),'uL in ', num2str(floor((GetSecs()-startRecTime)/60)), ' minutes']);
disp(' ');
disp(['lick rate during wait period: ', num2str(waitPeriodLickCounter/waitPeriodTotalTime), ' lick/sec'])
disp(' ');
disp(['Finish Time: ', datestr(sessionEndTime,'HH:MM:SS.FFF')])


%% reading the recorded analog data file
fid2 = fopen(binFile,'r');
% testData = fread(fid2,'double');
[data,count] = fread(fid2,[5,inf],'double');
fclose(fid2);

% figure()
t = data(1,:);
ch = data(2:5,:);

save(dataFolderAdd + '\' + 'workspaceVariables');
copyfile(string(mfilename('fullpath')) + '.m', dataFolderAdd);

figure()
plot(t, ch);


%% lick response based on the stim start detected by photodiode sensor

lickSensor = ch(2,:);
%getting rid of the noise
lickSensor(lickSensor<2)=0;
lickSensor(lickSensor>2)=1;

lickSensorDiff = [0 lickSensor(2:end)-lickSensor(1:end-1)];
lickTimes = t(find(lickSensorDiff==1));

photoDiodeSig = ch(1,:);

% plot(movMeanPhotoDiodeSignal);

movMeanPhotoDiodeSignal = photoDiodeSig;
% movMeanPhotoDiodeSignal = movmean(photoDiodeSig,100);
% movMeanPhotoDiodeSignal = movmean(movMeanPhotoDiodeSignal,100);
% movMeanPhotoDiodeSignal = movmean(movMeanPhotoDiodeSignal,50);
% movMeanPhotoDiodeSignal = movmean(movMeanPhotoDiodeSignal,50);

% finding the peaks of the bimodal distribution of the values in the
% photodiode signal to determine the cut level and get the start time of
% the visual stimulus
[idx,Centers] = kmeans(movMeanPhotoDiodeSignal',2);

cutLevel = (Centers(1) + Centers(2))/2;
% bar(histCenters,histCounts);
%getting rid of the noise
movMeanPhotoDiodeSignal(movMeanPhotoDiodeSignal<cutLevel)=0;
movMeanPhotoDiodeSignal(movMeanPhotoDiodeSignal>cutLevel)=1;

stimPhotodiodeDiff = [0 movMeanPhotoDiodeSignal(2:end)-movMeanPhotoDiodeSignal(1:end-1)];
stimStartPhotodiode = t(find(stimPhotodiodeDiff==1));
stimStopPhotodiode = t(find(stimPhotodiodeDiff==-1));

% if data starts with high photodiode level, getting rid of the first stop
if stimStopPhotodiode(1) < stimStartPhotodiode(1)
    stimStopPhotodiode = stimStopPhotodiode(2:end);
end

% if data ends with high level photodiode, getting rid of the last start
if length(stimStopPhotodiode) ~= length(stimStartPhotodiode)
    stimStartPhotodiode = stimStartPhotodiode(1:end-1);
end

% calculate the duration of each detected stim
if length(stimStopPhotodiode) == length(stimStartPhotodiode)
    stimDur = stimStopPhotodiode - stimStartPhotodiode;
else
    disp('photo diode signal can not be interepreted to calc the stim start time')
end
    
% dropping the fluctutations due to over hearing of the other channels by
% checking the length of the stimDur
minimumExpectedStimDur = 1;
stimStartPhotodiode = stimStartPhotodiode(stimDur>minimumExpectedStimDur);

if length(stimStartPhotodiode) ~= StimCounter
    disp('the number of stimuli detected on the photodiode is not matching the presented stimuli');
end


windowBeforeStimStart = 3;
windowAfterStimStart = 5;

stimLockedLickTimesAllTrials = [];
firstLickAfterStimAllTrials = [];

maxTimeToLookForTheFirstLick = windowAfterStimStart;

for stimCounter=1:size(stimStartPhotodiode,2)
    
    stimLockedLickTimes = lickTimes(lickTimes>(stimStartPhotodiode(stimCounter)-windowBeforeStimStart) & lickTimes<(stimStartPhotodiode(stimCounter)+windowAfterStimStart)) - stimStartPhotodiode(stimCounter);
    
    firstLickAfterStim = lickTimes(lickTimes>stimStartPhotodiode(stimCounter) & lickTimes<(stimStartPhotodiode(stimCounter)+windowAfterStimStart)) - stimStartPhotodiode(stimCounter);
    
    firstLickAfterStim = firstLickAfterStim(firstLickAfterStim<maxTimeToLookForTheFirstLick);
    if size(firstLickAfterStim)
        firstLickAfterStim = firstLickAfterStim(1);
    else
        firstLickAfterStim = double.empty(1,0);
    end
        
    stimLockedLickTimesAllTrials = [stimLockedLickTimesAllTrials stimLockedLickTimes];
    firstLickAfterStimAllTrials = [firstLickAfterStimAllTrials firstLickAfterStim];
end

histStep = 0.05;
windowStartBeforeStim = windowBeforeStimStart;
histBins = (-windowStartBeforeStim:histStep:windowAfterStimStart)+histStep/2;

figure()
hist(stimLockedLickTimesAllTrials,histBins)
title(strcat('distribution for the time of the all licks'));

histStep = 0.05;
windowStartBeforeStim = 1;
windowAfterStimStart = maxTimeToLookForTheFirstLick;
histBins = (-windowStartBeforeStim:histStep:windowAfterStimStart)+histStep/2;

colorVec = ["#4569D3","#97D345","#D3C845","#D38145","#D34550"];

figure()
[histCounts,histCenters] = hist(firstLickAfterStimAllTrials,histBins);

bar(histCenters,histCounts,1,'FaceColor',colorVec(1),'EdgeColor','none');
title(strcat('distribution for the time of the first lick'));

xlabel('time from stim start (seconds)');

yl = ylim;
text(0.5*maxTimeToLookForTheFirstLick,yl(2)*0.75,['median delay: ',num2str(round(1000*median(firstLickAfterStimAllTrials))),' ms'])

hold on
plot([0 StimDuration],[yl(2)*0.9 yl(2)*0.9],'k','linewidth',2)
%% updating the database for the animal

% load the data base for the animal
defaultPath = 'D:\animalsMatFiles';

animalMatDataFileAdd = strcat(defaultPath,'\',num2str(mouseNumber),'.mat');
load(animalMatDataFileAdd)

Date = [Date, string(date)];
weightYesterday = [weightYesterday, string(MouseWeight)];
addWaterYesterday = [addWaterYesterday, string(addWater)];
trainingStage = [trainingStage, string(stageOfTraining)];
hitCount{end+1} = hitCounter;  % add empty vectors for the days without data for this item
stimCount{end+1} = StimCounter;
firstLickDist{end+1} = firstLickAfterStimAllTrials;
totalTodayReward = [totalTodayReward, earnedRewardVolTotal];
FA_count = [FA_count, FA_counter];
extendedStimCount = [extendedStimCount, extendedITT_trialsCounter];

save(animalMatDataFileAdd,'Date','weightYesterday','addWaterYesterday','trainingStage','hitCount',...
    'stimCount','firstLickDist','totalTodayReward',...
    'mouseInitialWeight','FA_count','extendedStimCount','-append');

%% determine the day number of detection

sessionID_detection = find(trainingStage == "Detection");

sessionID_detectionToInclude = []; %to exclude the extra sessions that should be ignored in

for sessionCounter=1:length(sessionID_detection)
    
    sessionInd = sessionID_detection(sessionCounter);
    
    % check to see if there is any other sessions been recorded in this day
    
    sessionDateTemp = Date(sessionInd);
    
    sessionsIndInThisDay = find(Date == sessionDateTemp);
    
    if sessionInd == sessionsIndInThisDay(end) % just include the last session from this day
        sessionID_detectionToInclude = [sessionID_detectionToInclude sessionInd];
    end
end
    
todayDayNoOfDetection = length(sessionID_detectionToInclude);

%% distribution of all licks

histStep = 0.05;
windowStartBeforeStim = windowBeforeStimStart;
histBins = (-windowStartBeforeStim:histStep:windowAfterStimStart)+histStep/2;

h = figure();
hist(stimLockedLickTimesAllTrials,histBins)
% title(strcat('distribution for the time of the all licks'));

xlabel('time from stim start (seconds)');

histStep = 0.05;
windowStartBeforeStim = 1;
windowAfterStimStart = maxTimeToLookForTheFirstLick;
histBins = (-windowStartBeforeStim:histStep:windowAfterStimStart)+histStep/2;


sgtitle({strcat('#',mouseNumber, ' Day', num2str(todayDayNoOfDetection)) ,string(date),...
    'Distribution of All Lick Times across Trials'})

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

fileName = strcat('allLicksDist_Day', num2str(todayDayNoOfDetection));
fileAdd = strcat(googleDriveLocalFolderAddress,'\',mouseNumber,'\',fileName);
print(h,fileAdd,'-dpdf','-r0')

%% first lick distribution for different contrasts
h = figure();
[histCounts,histCenters] = hist(firstLickAfterStimAllTrials,histBins);

colorVec = ["#4569D3","#97D345","#D3C845","#D38145","#D34550"];
 
bar(histCenters,histCounts,1,'FaceColor',colorVec(1),'EdgeColor','none');
% title(strcat('distribution for the time of the first lick'));

xlabel('time from stim start (seconds)');

yl = ylim;
text(0.5*maxTimeToLookForTheFirstLick,yl(2)*0.75,['median delay: ',num2str(round(1000*median(firstLickAfterStimAllTrials))),' ms'])

hold on
plot([0 StimDuration],[yl(2)*0.9 yl(2)*0.9],'k','linewidth',2)

sgtitle({strcat('#',mouseNumber, ' Day', num2str(todayDayNoOfDetection)) ,string(date),...
    'Distribution of First Lick Times across Trials'})

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

fileName = strcat('firstLicksDist_Day', num2str(todayDayNoOfDetection));
fileAdd = strcat(googleDriveLocalFolderAddress,'\',mouseNumber,'\',fileName);
print(h,fileAdd,'-dpdf','-r0')


%% generate the plots and replace on the google drive

% performance change across days 

h = figure('Position', [50 50 600 1000]);

    
tempContrastHitRate = [];

for sessionCounter=1:todayDayNoOfDetection
    sessionInd = sessionID_detectionToInclude(sessionCounter);
    tempContrastHitRate(end+1) = (hitCount{sessionInd})/(stimCount{sessionInd})*100;
end

p = plot(tempContrastHitRate);
p.Color = colorVec(1);
p.LineStyle = '--';
p.Marker = 'o';



xlabel('day');


yl = ylim;

%     text(0.5,(yl(1)+yl(2))/2,['Contrast',num2str(ContrastCounter)])

xlim([0, todayDayNoOfDetection+1])

xticks(1:todayDayNoOfDetection)

ylabel('Hit Rate (%)');
    
hold on

tempFA_Rate = [];

for sessionCounter=1:todayDayNoOfDetection
        sessionInd = sessionID_detectionToInclude(sessionCounter);
        tempFA_Rate(end+1) = (FA_count(sessionInd))/(extendedStimCount(sessionInd))*100;
end
p = plot(tempFA_Rate);
p.Color = 'k';
p.LineStyle = '--';
p.Marker = 'o';

legend('Hit','FA')

sgtitle({strcat('#',mouseNumber), ' Performance Across Days'})

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

fileName = 'PerformanceAcrossDays_Detection';
fileAdd = strcat(googleDriveLocalFolderAddress,'\',mouseNumber,'\',fileName);
print(h,fileAdd,'-dpdf','-r0')

%% response delay change across days for different contrasts

h = figure('Position', [50 50 600 1000]);
hold on
tempDelayMedian = [];
tempDelaySEM = [];

for sessionCounter=1:todayDayNoOfDetection
    sessionInd = sessionID_detectionToInclude(sessionCounter);

    firstLickTimes = firstLickDist{sessionInd};
    firstLickTimes = firstLickTimes(firstLickTimes<2.5);

    tempDelayMedian(end+1) = median(firstLickTimes);
    tempDelaySEM(end+1) = std(firstLickTimes)/sqrt(length(firstLickTimes));

end

p = plot(tempDelayMedian*1000);
%     p = errorbar(tempDelayMean*1000,tempDelaySEM*1000);
p.Color = colorVec(1);
p.LineStyle = '--';
p.Marker = 'o';



xlabel('day');


ylabel('delay(ms)');

yl = ylim;

%     text(0.5,(yl(1)+yl(2))/2,['Contrast',num2str(ContrastCounter)])

xlim([0, todayDayNoOfDetection+1])

xticks(1:todayDayNoOfDetection)

ylabel('ms');

sgtitle({strcat('#',mouseNumber), ' Response Delay Across Days'})

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

fileName = 'ResponseDelayAcrossDays_Detection';
fileAdd = strcat(googleDriveLocalFolderAddress,'\',mouseNumber,'\',fileName);
print(h,fileAdd,'-dpdf','-r0')


%% Saving the variables of the session and the code file in the recording
%directory
save(dataFolderAdd + '\' + 'workspaceVariables');