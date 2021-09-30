% error('more things need to be saved. Do that b4 running');   
subjID = input('Enter subject ID: \n','s'); %e.g. 001 002 010 099
WM_run = input('Enter run number: \n'); %e.g. 1 2 3

environment = 'behav';
version = 'V1';

fileID = [subjID,'_',num2str(WM_run)];

%% Subj ID & Paths

fileName = [fileID, '_WMdistAnt_',version];
fileDir = fullfile(pwd,'data',subjID);
if ~exist(fileDir,'dir')
    mkdir(fileDir);
end
outFile = fullfile(fileDir,[fileName, '.mat']);
if exist(outFile)>0
    error('Log file already exists! Select different run or subject number')
end

PsychDefaultSetup(2); %converts color from 0-255 to float 0-1

whichScreen = 0; %Screen ID.

% REMOVE THIS WHEN ACTUALLY TESTING!!!! 
Screen('Preference', 'SkipSyncTests', 1)


%--------------------------------------------------------------------------
% Add functions to path
addpath(genpath('Functions'));

% Reset rng
rng('default');
rng('shuffle');

%--------------------------------------------------------------------------
%% Define environment and load appropriate settings
background = 102; % dark gray
background = background/255;

mriPulse = 0; % Flag to determine if we are going to wait for an MRI pulse before starting
% Gets turned to 1 if environment is MRI by default

% Monitor settings for each environment
% SreenSize - either one value (diagonal) OR two values (x y) - ALL IN CM!
switch environment
    case 'scanner' 
        error('Change the parameters')
        params.distFromScreen = 144;
        params.ScreenSize = [69.84 39.29]; % MANUALLY ENTER THE SCREEN SIZE in cm
        delete(instrfind)
        mriPulse = 1; % Flag for waiting for mriPulse
    case 'behav'
        params.distFromScreen = 50;
        params.ScreenSize = 39.624; %MANUALLY ENTER THE SCREEN SIZE in cm
end

%% Prepare

% Open window
% [window, rect] = PsychImaging('OpenWindow', whichScreen, background);
[window, rect] = PsychImaging('OpenWindow', whichScreen, background, [], 32, 2, [], 6,  []);


HideCursor(window);

%% Stimulus Size Calculation:

if length(params.ScreenSize) == 1
    Ratio.pixelsAcrossScreen = sqrt(rect(3)^2+rect(4)^2);
    Ratio.ppCM = (Ratio.pixelsAcrossScreen/params.ScreenSize); %pixels per cm. %2.54 cm/inch
elseif length(params.ScreenSize) == 2
    Ratio.ppCM_x = rect(3)/params.ScreenSize(1);
    Ratio.ppCM_y = rect(4)/params.ScreenSize(2);
    if round(Ratio.ppCM_x, 1)== round(Ratio.ppCM_y,1)
        Ratio.ppCM = round(Ratio.ppCM_x);
    else
        error('Please make sure the entered ScreenSize values (in cm) are correct OR decrease the stringency of the round function (by changing N in round(x, N))')
    end
end

%% 

% define parameters for stimulus

% define parameters related to task

task.cues = [100, 0, 50];
task.cuesColor = ['blue','orange','green'];
task.distractorOrientations = [30, 15, 5, -5, -15, -30];
task.numUniqueTrials = length(task.cues)*length(task.distractorOrientations); 






params.numBlocksPer = 3;

params.numTrials = 32; % trials per block

params.scale=4;
params.luminance=0.8;
params.max_con=1;
params.min_con=0.1;
params.contrasts=logspace(log10(params.min_con),log10(params.max_con),12); % number of possible contrasts
params.mean_con=params.contrasts(4);

params.contrast_probe = logspace(log10(params.min_con),log10(params.max_con),12); %number of possible options for an answer


% define the spatial geometry of the stimulus
params.n_circ=4; % number of circles/2
params.n_ang=12; % number of sectors per circle/2
params.r_scaling=2; % length scaling of sectors.
params.circle_param = 8; %how many circles are included in one item  
params.angle_param = 24; % how many thetas are included in one item % CHANGED TO 2
params.numSegments = (2*params.n_circ/params.circle_param)*(2*params.n_ang/params.angle_param);

params.stimrange=255; %221 orig
params.visAngleMin = 1; %enter in degrees will convert to radians later.
params.visAngleMax = 10; 

PalRGB=repmat(linspace(params.luminance/2-params.luminance*params.max_con/2,params.luminance/2+params.luminance*params.max_con/2,params.stimrange)',1,3);
bgcol=ceil(params.stimrange/2);
PalRGB(params.stimrange+1,:)=params.luminance;
conrange=(params.stimrange-1)/2;
PalContrasts=round(params.contrasts/params.max_con*conrange);
PalContrasts = PalContrasts';
PalContrasts_probe = round(params.contrast_probe/params.max_con*conrange);

params.mean_con=round(params.mean_con/params.max_con*conrange);

% define timing parameters
time.TR = 1.5;
time.scannerWaitTime = time.TR*1; 
time.stimulus = 0.8;
time.ISI = 0.5; 
time.cue = 0.8;
time.delay = 5; 
time.probe = 0.5; 
time.answer = 5; 
time.feedback = 1;
time.ITI = 3.5;
time.delayOneStim = time.delay+2*time.ISI+time.stimulus+time.cue; 

[ch_board, indic, ch_board_parts, ch_board_min]=random_checker_sin(params.n_circ,params.n_ang,params.r_scaling,params.scale, params.circle_param, params.angle_param, rect, params.visAngleMin, params.visAngleMax, params.distFromScreen, Ratio.ppCM, params.numSegments);

%% Estimate flip interval and present text on screen
Screen('TextFont',window, 'Arial');
Screen('TextSize',window, 20);
% Fill screen with mid-grey background
Screen('FillRect', window, background);

% Measure flip interval for accurate stimulus timings!
Priority(1); % Set priority
Screen('DrawText',window,'Estimating monitor flip interval...', 100, 100);
Screen('DrawText',window,'(This may take up to 20s)', 100, 120);
Screen('Flip',window);
halfifi = Screen('GetFlipInterval', window, 100, 0.00005, 20);
halfifi = halfifi/2;
Priority(0); % reset priority level

Screen('DrawText',window,'Generating stimuli...', 100, 100);
Screen('DrawText',window,'(This may take up to a couple of minutes)', 100, 120);
Screen('Flip',window);