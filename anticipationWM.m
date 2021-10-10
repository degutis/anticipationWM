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
Screen('Preference', 'SkipSyncTests', 0)


%--------------------------------------------------------------------------
% Add functions to path
addpath(genpath('Functions'));

% Reset rng
rng('default');
rng('shuffle');

%--------------------------------------------------------------------------
%% Define environment and load appropriate settings
%background = 102; % dark gray
%background = background/255;
[background,~,~] = calibrate_lum2(); 

mriPulse = 0; % Flag to determine if we are going to wait for an MRI pulse before starting
% Gets turned to 1 if environment is MRI by default

% Monitor settings for each environment
% SreenSize - either one value (diagonal) OR two values (x y) - ALL IN CM!
switch environment
    case 'scanner' 
        error('Change the parameters')
        stim.distFromScreen = 144;
        stim.ScreenSize = [69.84 39.29]; % MANUALLY ENTER THE SCREEN SIZE in cm
        delete(instrfind)
        mriPulse = 1; % Flag for waiting for mriPulse
    case 'behav'
        stim.distFromScreen = 60;
        stim.ScreenSize = 45; %MANUALLY ENTER THE SCREEN SIZE in cm
    case 'workPC'
        stim.distFromScreen = 50;
        stim.ScreenSize = 39.624; %MANUALLY ENTER THE SCREEN SIZE in cm
end

%% Prepare

% Open window
% [window, rect] = PsychImaging('OpenWindow', whichScreen, background);
%[window, rect] = PsychImaging('OpenWindow', whichScreen, background, [], 32, 2, [], 6,  []);
% HideCursor(window);

[window, rect] = Screen('OpenWindow',whichScreen);
%%%%%HideCursor(window);

%% Stimulus Size Calculation:

if length(stim.ScreenSize) == 1
    Ratio.pixelsAcrossScreen = sqrt(rect(3)^2+rect(4)^2);
    Ratio.ppCM = (Ratio.pixelsAcrossScreen/stim.ScreenSize); %pixels per cm. %2.54 cm/inch
elseif length(stim.ScreenSize) == 2
    Ratio.ppCM_x = rect(3)/stim.ScreenSize(1);
    Ratio.ppCM_y = rect(4)/stim.ScreenSize(2);
    if round(Ratio.ppCM_x, 1)== round(Ratio.ppCM_y,1)
        Ratio.ppCM = round(Ratio.ppCM_x);
    else
        error('Please make sure the entered ScreenSize values (in cm) are correct OR decrease the stringency of the round function (by changing N in round(x, N))')
    end
end

%% Task Parameters

task.cues = [100, 0, 50];
%task.cuesColor = [{'blue'},{'orange'},{'green'}];
%task.cuesColorRGB = {[0, 0, 255] [255,165,0], [0,128,0]};
task.distractorOrientations = [30, 15, 5, -5, -15, -30];
task.numUniqueTrials = length(task.cues)*length(task.distractorOrientations); 
task.trialScalar = 4;


task.numUniqueTrials = task.numUniqueTrials*task.trialScalar; %increase num of trials

task.numOrientationBins = 6;
task.memoryOrientations = [randi([5,35],task.numUniqueTrials,1),randi([36,65],task.numUniqueTrials,1),randi([66,85],task.numUniqueTrials,1),randi([95,125],task.numUniqueTrials,1),randi([126,155],task.numUniqueTrials,1),randi([156,175],task.numUniqueTrials,1)];

task.cuesFull = repelem(task.cues,task.numOrientationBins,length(task.distractorOrientations)*task.trialScalar)';
task.distractorOrientationsFull = repmat(task.distractorOrientations,task.numOrientationBins,length(task.cues)*task.trialScalar)';

task.numBlocks = 6;
task.numTrialsBlock=(size(task.memoryOrientations,1)*size(task.memoryOrientations,2))*2/task.numBlocks; %the extra 2 is to balance the 50 trials

task.trialMatrix = [task.memoryOrientations(:),task.cuesFull(:),task.distractorOrientationsFull(:)];
task.trialMatrix(task.trialMatrix(:,2)==0,3)=0;
task.matrixBalance = task.trialMatrix;
task.matrixBalance(task.matrixBalance(:,2)==50,3)=0; 


BlockA = [];
BlockB = [];
for blockcue = 1:length(task.trialMatrix)/task.numBlocks
    
    currentBlockA = Shuffle(1:task.numBlocks)';
    currentBlockB = Shuffle(1:task.numBlocks)';
    
    BlockA = [BlockA;currentBlockA];
    BlockB = [BlockB;currentBlockB];
end

task.trialMatrix = [task.trialMatrix,BlockA];
task.matrixBalance = [task.matrixBalance,BlockB];

task.trialMatrixFinal = zeros(task.numTrialsBlock,size(task.trialMatrix,2),task.numBlocks);

for blockk = 1:task.numBlocks
    currentblock_matrix = [task.trialMatrix(task.trialMatrix(:,4)==blockk,:);task.matrixBalance(task.matrixBalance(:,4)==blockk,:)];
    currentblock_matrix = currentblock_matrix(randperm(size(currentblock_matrix, 1)),:);
    task.trialMatrixFinal(:,:,blockk) = currentblock_matrix; 
end

%% Time settings

time.cue = .5; 
time.ISI = .2; 
time.MemoryGabor = .5; % Time to present Gabor
time.DelayLength = 3;
time.DelayDist = time.DelayLength/4;
time.Distractor = time.DelayLength/2;
time.DistOn = 0.25;
time.DistOff = 0.25; 
time.DistFreq = 1/time.DistOn;
time.DistPrezT = time.Distractor*time.DistFreq/2;
time.scannerWaitTime = 2;
time.answer = 5;
time.ITI = 1.5;
%% Gabor parameters

stim.contrast = 0.5;
stim.stimSize = 16; % Diameter in degrees
stim.sf = 1; % cpd
stim.phases = linspace(0, 2*pi, 13); %12 possible phases
stim.phases(1)=[];
stim.phaseMemory = Shuffle(repmat(stim.phases, 1,task.numTrialsBlock*task.numBlocks/length(stim.phases)));
stim.phaseMemory = reshape(stim.phaseMemory, task.numTrialsBlock,task.numBlocks);

stim.phaseDistractor = repmat(stim.phases, 1,(task.numTrialsBlock/2)*task.numBlocks*time.DistPrezT/length(stim.phases));
stim.phaseDistractor = reshape(stim.phaseDistractor, task.numTrialsBlock/2,time.DistPrezT, task.numBlocks);

% cue, memory, delay, distractor*on/off [stim.phaseMemory*2], delay2,
% probe, iti. 
% per TRIAL
stim.numPresentations = 6+size(stim.phaseDistractor,2)*2;
stim.initialProbePos = randi([0,179],task.numTrialsBlock*task.numBlocks,1);
stim.initialProbePos = Shuffle(stim.initialProbePos(:));
stim.initialProbePos = reshape(stim.initialProbePos,[],task.numBlocks);



%% Estimate flip interval and generate Gabors

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

stimTex = zeros(length(stim.phases),1); % Preallocate

for b = 1:length(stim.phases)
    
    width = degrees2pixels(stim.stimSize, stim.distFromScreen, Ratio.ppCM);
    inner_degree = 1;
    size_inner_most_circle = degrees2pixels(inner_degree, stim.distFromScreen, Ratio.ppCM);
    start_linear_decay_in_degree = .5; 
    start_linear_decay = degrees2pixels(start_linear_decay_in_degree, stim.distFromScreen, Ratio.ppCM);

    nCycles = stim.stimSize*stim.sf; % number of cycles in a stimulus

    % compute the pixels per grating period
    pixelsPerGratingPeriod = width / nCycles;

    spatialFrequency = 1 / pixelsPerGratingPeriod; % How many periods/cycles are there in a pixel?
    radiansPerPixel = spatialFrequency * (2 * pi); % = (periods per pixel) * (2 pi radians per period)

    [background, Lmin_rgb, Lmax_rgb] = calibrate_lum2(stim.contrast);

    halfWidthOfGrid = width / 2;
    widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  
    % widthArray is used in creating the meshgrid.
    % Creates a sinusoidal grating, where the period of the sinusoid is 
    % approximately equal to "pixelsPerGratingPeriod" pixels.
    % Note that each entry of gratingMatrix varies between minus one and
    % one; -1 <= gratingMatrix(x0, y0)  <= 1

    % Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
    % the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
    % and y = y(x0, y0) corresponds to the y-coordinate of element "i"
    [x_mesh y_mesh] = meshgrid(widthArray);

    % Creates a sinusoidal grating, where the period of the sinusoid is 
    % approximately equal to "pixelsPerGratingPeriod" pixels.
    % Note that each entry of gratingMatrix varies between minus one and
    % one; -1 <= gratingMatrix(x0, y0)  <= 1

    % the grating is oriented vertically unless otherwise specified.

    stimulusMatrix = sin(radiansPerPixel * x_mesh + stim.phases(b));
    % Because the difference between Lmax_rgb and background is not equal to
    % the difference between Lmin_rgb and background (the differences are equal
    % in luminance, but not in rgb), correct the range of the positive values
    % of stimulusMatrix.

    rgb_range_up = Lmax_rgb - background;
    rgb_range_down = background - Lmin_rgb;
    
    rgb_range_up_down_ratio = rgb_range_up/rgb_range_down;
    stimulusMatrix(stimulusMatrix > 0) = stimulusMatrix(stimulusMatrix > 0) * rgb_range_up_down_ratio;

    % Make a fading annulus, to use as a mask.
    annulusMatrix = makeLinearMaskCircleAnn(width+1,width+1,size_inner_most_circle,start_linear_decay,width/2);
    stimulusMatrix = stimulusMatrix .* annulusMatrix;

    %Make the grating
    gaborPatch = background + rgb_range_down * stimulusMatrix;

    stimTex(b) = Screen('MakeTexture', window, gaborPatch);
end

Screen('PreloadTextures',window,stimTex);

% Generate probe 

lineTex = gaborPatch;
lineTex(lineTex~=background)=max(max(gaborPatch)); 
size_line = degrees2pixels(0.1, stim.distFromScreen, Ratio.ppCM);
size_line = round(size_line);
lineTex(:,1:size(gaborPatch,2)/2-size_line)=background;
lineTex(:,size(gaborPatch,2)/2+size_line:end)=background;
lineTex(size(gaborPatch,2)/2:size(gaborPatch,2)/2+size_inner_most_circle,:)=background;
lineTex(size(gaborPatch,2)/2-size_inner_most_circle:size(gaborPatch,2)/2,:)=background;

probeTex = Screen('MakeTexture', window, lineTex);
Screen('PreloadTextures',window,probeTex);

%% Wait for mouse click to begin
clicks = 0;
Screen('FillRect',window,background);
instructions = ['Instructions row 1!\n\nInstructions row two\n\nRun length: X minutes\n\n'];
Screen('DrawText',window,'[to Experimenter] waiting for mouse click...', 100, 100);
DrawFormattedText(window, instructions, 'center', 'center');
Screen('Flip',window);
while clicks < 1
    clicks = GetClicks(window);
end

% should make a table

data = table(NaN(task.numTrialsBlock*task.numBlocks,1),NaN(task.numTrialsBlock*task.numBlocks,1),reshape(task.trialMatrixFinal(:,1,:),[],1), reshape(task.trialMatrixFinal(:,2,:),[],1),reshape(task.trialMatrixFinal(:,3,:),[],1),reshape(task.trialMatrixFinal(:,4,:),[],1));
data.Properties.VariableNames = {'Answer','RT','StimulusOr','Trial_type','DistractorOr','Block'};

presentation_timings = table(NaN(stim.numPresentations*task.numTrialsBlock*task.numBlocks,1),NaN(stim.numPresentations*task.numTrialsBlock*task.numBlocks,1),strings(stim.numPresentations*task.numTrialsBlock*task.numBlocks,1),strings(stim.numPresentations*task.numTrialsBlock*task.numBlocks,1),strings(stim.numPresentations*task.numTrialsBlock*task.numBlocks,1),repmat(repelem([1:task.numTrialsBlock],1,stim.numPresentations)',task.numBlocks,1),repelem([1:task.numBlocks],1,stim.numPresentations*task.numTrialsBlock)',repelem([1:task.numBlocks*task.numTrialsBlock],1,stim.numPresentations)');
presentation_timings.Properties.VariableNames = {'TimeRealFlip','TimeExpected','Stim_type','Trial_type','DistractorOr','TrialNumBlock','Block','TrialsTotal'};


expectedTime = time.scannerWaitTime;
draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
[~,startTime] = Screen('Flip',window); % Store time experiment started
Screen('Flip', window, startTime + expectedTime - halfifi);

expectedKeys = {'a','l','space'};


KbQueueCreate;
KbQueueStart;

for block = 1:task.numBlocks
    if block~=1
        % Function to give the subject a short break ;)
        DrawFormattedText(window,['Feel free to take a short break\n Click the mouse to continue\n'],'center', 'center');
        Screen('Flip',window);

        clicks = 0;
        % Wait for click
        while clicks < 1
            clicks = GetClicks(window);
        end

        % Draw fixation and wait for a couple of seconds
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,stim.distFromScreen,Ratio.ppCM);    
        Screen('DrawingFinished', window);
        % Reset timings since we may have waited a while
        startTime = Screen('Flip', window);
        expectedTime = 0;

        expectedTime = expectedTime + time.scannerWaitTime;
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,stim.distFromScreen,Ratio.ppCM);    
        Screen('DrawingFinished', window);
        Screen('Flip', window, startTime + expectedTime - halfifi);
    end
    
    distractorTrial=0;

    for trialNum = 1:task.numTrialsBlock
        
        presentation_timing_index = (block-1)*(task.numTrialsBlock*stim.numPresentations)+(trialNum-1)*stim.numPresentations;
        
        % Colored cue
        %draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],task.cuesColorRGB{task.cues==task.trialMatrixFinal(trialNum,2,block)},task.cuesColorRGB{task.cues==task.trialMatrixFinal(trialNum,2,block)},stim.distFromScreen,Ratio.ppCM);
        %Screen('DrawingFinished', window);
 
        % Percent number cue
        DrawFormattedText(window,[num2str(task.trialMatrixFinal(trialNum,2,block)),'%'],'center', 'center');
        Screen('DrawingFinished', window);
        [~, presentation_timings.TimeRealFlip(presentation_timing_index+1)] = Screen('Flip', window, startTime + expectedTime - halfifi);
        presentation_timings.TimeRealFlip(presentation_timing_index+1) = presentation_timings.TimeRealFlip(presentation_timing_index+1)-startTime;
        presentation_timings.TimeExpected(presentation_timing_index+1) = expectedTime;
        presentation_timings.Stim_type(presentation_timing_index+1)= 'Cue';
        presentation_timings.Trial_type(presentation_timing_index+1) = num2str(task.trialMatrixFinal(trialNum,2,block));
        presentation_timings.DistractorOr(presentation_timing_index+1) = num2str(task.trialMatrixFinal(trialNum,3,block));
        expectedTime = expectedTime+time.cue; 
        
        % Present memory stimulus
        Screen('DrawTexture', window, stimTex(stim.phases==stim.phaseMemory(trialNum,block)), [],[],task.trialMatrixFinal(trialNum,1,block));
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
        Screen('DrawingFinished', window);
        [~, presentation_timings.TimeRealFlip(presentation_timing_index+2)] = Screen('Flip', window, startTime + expectedTime - halfifi);
        presentation_timings.TimeRealFlip(presentation_timing_index+2) = presentation_timings.TimeRealFlip(presentation_timing_index+2)-startTime;
        presentation_timings.TimeExpected(presentation_timing_index+2) = expectedTime;
        presentation_timings.Stim_type(presentation_timing_index+2)= 'MemoryStim';
        presentation_timings.DistractorOr(presentation_timing_index+2) = num2str(task.trialMatrixFinal(trialNum,3,block));
        presentation_timings.Trial_type(presentation_timing_index+2) = num2str(task.trialMatrixFinal(trialNum,2,block));   
        expectedTime = expectedTime+time.MemoryGabor; 


        if task.trialMatrixFinal(trialNum,3,block)==0
            % Delay extended
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
            Screen('DrawingFinished', window);        
            [~, presentation_timings.TimeRealFlip(presentation_timing_index+3)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings.TimeRealFlip(presentation_timing_index+3) = presentation_timings.TimeRealFlip(presentation_timing_index+3)-startTime;
            presentation_timings.TimeExpected(presentation_timing_index+3) = expectedTime;
            presentation_timings.Stim_type(presentation_timing_index+3)= 'DelayExtended';
            presentation_timings.DistractorOr(presentation_timing_index+3) = num2str(task.trialMatrixFinal(trialNum,3,block));
            presentation_timings.Trial_type(presentation_timing_index+3) = num2str(task.trialMatrixFinal(trialNum,2,block));   
            expectedTime = expectedTime+time.DelayLength; 
        
        else
            distractorTrial=distractorTrial+1; 
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
            Screen('DrawingFinished', window);
            [~, presentation_timings.TimeRealFlip(presentation_timing_index+3)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings.TimeRealFlip(presentation_timing_index+3) = presentation_timings.TimeRealFlip(presentation_timing_index+3)-startTime;
            presentation_timings.TimeExpected(presentation_timing_index+3) = expectedTime;
            presentation_timings.Stim_type(presentation_timing_index+3)= 'Delay1';
            presentation_timings.DistractorOr(presentation_timing_index+3) = num2str(task.trialMatrixFinal(trialNum,3,block));
            presentation_timings.Trial_type(presentation_timing_index+3) = num2str(task.trialMatrixFinal(trialNum,2,block));   
            expectedTime = expectedTime+time.DelayDist;
 
            for numDist = 1:time.DistPrezT
                Screen('DrawTexture', window, stimTex(stim.phases==stim.phaseDistractor(distractorTrial,numDist,block)), [],[],task.trialMatrixFinal(trialNum,1,block)+task.trialMatrixFinal(trialNum,3,block));
                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
                Screen('DrawingFinished', window);
                [~, presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+1)] = Screen('Flip', window, startTime + expectedTime - halfifi);
                presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+1) = presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+1)-startTime;
                presentation_timings.TimeExpected(presentation_timing_index+3+(numDist-1)*2+1) = expectedTime;
                presentation_timings.Stim_type(presentation_timing_index+3+(numDist-1)*2+1)= 'DistractorOn';
                presentation_timings.DistractorOr(presentation_timing_index+3+(numDist-1)*2+1) = num2str(task.trialMatrixFinal(trialNum,3,block));
                presentation_timings.Trial_type(presentation_timing_index+3+(numDist-1)*2+1) = num2str(task.trialMatrixFinal(trialNum,2,block));   
                expectedTime = expectedTime+time.DistOn; 

                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
                Screen('DrawingFinished', window);
                [~, presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+2)] = Screen('Flip', window, startTime + expectedTime - halfifi);
                presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+2) = presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+2)-startTime;
                presentation_timings.TimeExpected(presentation_timing_index+3+(numDist-1)*2+2) = expectedTime;
                presentation_timings.Stim_type(presentation_timing_index+3+(numDist-1)*2+2)= 'DistractorOff';
                presentation_timings.DistractorOr(presentation_timing_index+3+(numDist-1)*2+2) = num2str(task.trialMatrixFinal(trialNum,3,block));
                presentation_timings.Trial_type(presentation_timing_index+3+(numDist-1)*2+2) = num2str(task.trialMatrixFinal(trialNum,2,block));   
                expectedTime = expectedTime+time.DistOff;
            end
            
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
            Screen('DrawingFinished', window);
            [~, presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+3)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+3) = presentation_timings.TimeRealFlip(presentation_timing_index+3+(numDist-1)*2+3)-startTime;
            presentation_timings.TimeExpected(presentation_timing_index+3+(numDist-1)*2+3) = expectedTime;
            presentation_timings.Stim_type(presentation_timing_index+3+(numDist-1)*2+3)= 'Delay2';
            presentation_timings.DistractorOr(presentation_timing_index+3+(numDist-1)*2+3) = num2str(task.trialMatrixFinal(trialNum,3,block));
            presentation_timings.Trial_type(presentation_timing_index+3+(numDist-1)*2+3) = num2str(task.trialMatrixFinal(trialNum,2,block));   
            expectedTime = expectedTime+time.DelayDist;
        end
        
        % Present probe
        positionProbe = stim.initialProbePos(trialNum,block);
        Screen('DrawTexture', window, probeTex, [],[],positionProbe);
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
        Screen('DrawingFinished', window);
        [~, presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations-1)] = Screen('Flip', window, startTime + expectedTime - halfifi);
        presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations-1) = presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations-1)-startTime;
        presentation_timings.TimeExpected(presentation_timing_index+stim.numPresentations-1) = expectedTime;
        presentation_timings.Stim_type(presentation_timing_index+stim.numPresentations-1)= 'ProbeStim';
        presentation_timings.DistractorOr(presentation_timing_index+stim.numPresentations-1) = num2str(task.trialMatrixFinal(trialNum,3,block));
        presentation_timings.Trial_type(presentation_timing_index+stim.numPresentations-1) = num2str(task.trialMatrixFinal(trialNum,2,block));         

        response = [];
        close = 0;
        KbQueueFlush;

        while strcmp(response,'m')==0
            if GetSecs > startTime+expectedTime+time.answer-halfifi
                response = 'm'; %m for miss
                break
            elseif close==1
                Screen('DrawTexture', window, probeTex, [],[], positionProbe);
                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,stim.distFromScreen,Ratio.ppCM);
                Screen('DrawingFinished', window);
                %[~,data.RT((block-1)*task.numTrialsBlock+trialNum)] = Screen('Flip', window, startTime + expectedTime - halfifi);
                Screen('Flip', window, startTime + expectedTime - halfifi);
                
                % RT measured by subtracting the startTime and the onset of
                % the probe. 
                data.RT((block-1)*task.numTrialsBlock+trialNum) = RT_current-startTime-presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations-1);
                data.Answer((block-1)*task.numTrialsBlock+trialNum) = positionProbe;
                % MEASURE RT HERE
            else
                response = [];
                %[~,keyCode]=KbQueueCheck;
                [keyIsDown,~,keyCode,~] = KbCheck;
                if keyCode == 0
                    continue
                else
                    response = KbName(keyCode);
                    try
                        response = response{1};
                    catch
                        response;
                    end
                    % Check if the key pressed was one we expected
%                    rep = find(strcmp(response{1},expectedKeys) == 1, 1);
                    % If the response was not one of the keys we were expecting, then make
                    % response empty again and wait for another response.
%                    if isempty(rep)
%                        response = [];
%                    end
                    if strcmp(response,'a')
                        while keyIsDown && strcmp(response,'a')
                            if GetSecs > startTime+expectedTime+time.answer-halfifi
                                response = 'm'; %m for miss
                                break
                            elseif strcmp(response,'space')
                                close=1;
                                RT_current = GetSecs;
                                break                                
                            elseif positionProbe == 1
                                positionProbe = 180;
                            else
                                positionProbe = positionProbe-1;
                            end
                            Screen('DrawTexture', window, probeTex, [],[], positionProbe);
                            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
                            Screen('DrawingFinished', window);
                            Screen('Flip', window, startTime + expectedTime - halfifi);
                            keyIsDown=0;
                            [keyIsDown,~,keyCode,~] = KbCheck;
                            response = KbName(keyCode);
                            try
                                response = response{1};
                            catch
                                response;
                            end
                        end
                    elseif strcmp(response,'l')
                        while keyIsDown && strcmp(response,'l')
                            if GetSecs > startTime+expectedTime+time.answer-halfifi
                                response = 'm'; %m for miss
                                break
                            elseif strcmp(response,'space')
                                close=1;
                                RT_current = GetSecs;                                
                                break    
                            elseif positionProbe == 180
                                positionProbe = 1;
                            else
                                positionProbe = positionProbe+1;
                            end
                            Screen('DrawTexture', window, probeTex, [],[], positionProbe);
                            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
                            Screen('DrawingFinished', window);
                            Screen('Flip', window, startTime + expectedTime - halfifi);
                            keyIsDown=0;                        
                            [keyIsDown,~,keyCode,~] = KbCheck;
                            response = KbName(keyCode); 
                            try
                                response = response{1};
                            catch
                                response;
                            end                        
                        end
                    elseif strcmp(response,'space')
                        close=1;
                        RT_current = GetSecs;
                        %break
                    end
                end
            end
%            Screen('DrawTexture', window, probeTex, [],[], positionProbe);
%            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
%            Screen('DrawingFinished', window);
%            Screen('Flip', window, startTime + expectedTime - halfifi);
%        end
        end
        expectedTime = expectedTime+time.answer;            
        
        % ITI
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,stim.distFromScreen,Ratio.ppCM);
        Screen('DrawingFinished', window);        
        [~, presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations)] = Screen('Flip', window, startTime + expectedTime - halfifi);
        presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations) = presentation_timings.TimeRealFlip(presentation_timing_index+stim.numPresentations)-startTime;        
        presentation_timings.TimeExpected(presentation_timing_index+stim.numPresentations) = expectedTime;
        presentation_timings.Stim_type(presentation_timing_index+stim.numPresentations)= 'ITI';
        presentation_timings.DistractorOr(presentation_timing_index+stim.numPresentations) = num2str(task.trialMatrixFinal(trialNum,3,block));
        presentation_timings.Trial_type(presentation_timing_index+stim.numPresentations) = num2str(task.trialMatrixFinal(trialNum,2,block));   
        expectedTime = expectedTime+time.ITI;               
    end
end

sca;
save(outFile, 'data', 'stim', 'Ratio','presentation_timings', 'task', 'time');
writetable(data,fullfile(fileDir,['Data_',fileID,'.csv']),'Delimiter',',')  
writetable(data,fullfile(fileDir,['PresentationTimings_',fileID,'.csv']),'Delimiter',',')  





