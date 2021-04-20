Screen('Preference', 'SkipSyncTests', 1);

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
%% Gabor information
%--------------------

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
phaseLine = 0;

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
% gabortex = CreateProceduralGabor(window, gaborDimPixWidth, gaborDimPixHeight, [],...
%     backgroundOffset, disableNorm, preContrastMultiplier);

preContrastMultiplier = 0.5;

gratingTex = CreateProceduralSquareWaveGrating(window, gaborDimPixWidth, gaborDimPixHeight, backgroundOffset, [],preContrastMultiplier);



% gabortexMirror = CreateProceduralGabor(window, gaborDimPix, gaborDimPix, [],...
%     backgroundOffset, disableNorm, preContrastMultiplier);

heightOffsetInCM = 0;
heightOffset = floor(heightOffsetInCM*pixelDensityCM);
widthOffset = 0;

stimHeightOffset = [0,heightOffset,0,gaborDimPixHeight+heightOffset];
righImageHorzPos = [0,0,gaborDimPixWidth,0];

StimDuration = 1.5; %in sec
StimFrames = round(StimDuration/ifi);
% % Randomise the phase of the Gabors and make a properties matrix.
% propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];
% propertiesMatmock = [phase, freq, sigma, contrastMock, aspectRatio, 0, 0, 0];

% contrastVector = [0.001, 0.01, 0.1, 1];

contrast = 0.8;

% propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];
propertiesMat = [phase, freq, contrast, 0];

Screen('FillRect', window, gray);
Screen('FillRect',window, white, patchRect);
% Screen('DrawTextures', window, gabortex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
% kPsychDontDoRotation, propertiesMat');
Screen('DrawTextures', window, gratingTex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
            kPsychUseTextureMatrixForRotation, propertiesMat');
vblStim = Screen('Flip', window);

vbl = vblStim;

while (~KbCheck)
        % calculating the phase of the grating in this frame based on the
        % temporal frequency
        phaseLine = phaseLine + degPerFrame;
        propertiesMat(:, 1) = phaseLine';

        % loading the next frame
        Screen('FillRect', window, gray);
        Screen('FillRect',window, black, patchRect);
%         Screen('DrawTextures', window, gabortex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
%                 kPsychDontDoRotation, propertiesMat');
        Screen('DrawTextures', window, gratingTex, [], righImageHorzPos+stimHeightOffset, orientationPreferred, [], [], [], [],...
            kPsychUseTextureMatrixForRotation, propertiesMat');

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
end 

KbWait;
sca;
