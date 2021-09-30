subjID = input('Enter subject ID: \n','s'); %e.g. 001 002 010 099
WM_run = input('Enter run number: \n'); %e.g. 1 2 3

environment = 'behav';

fileID = [subjID,'_',num2str(WM_run)];

%% Subj ID & Paths

fileName = [fileID, '_WM_OneStim'];
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
params.scale=4;
params.luminance=0.8;
params.max_con=1;
params.min_con=0.1;
params.contrasts=logspace(log10(params.min_con),log10(params.max_con),7);
% params.mean_con=logspace(log10(params.contrasts(2)),log10(params.contrasts(3)),3);
params.mean_con=params.contrasts(4);

params.contrast_probe = logspace(log10(params.min_con),log10(params.max_con),13);

% define the spatial geometry of the stimulus
params.n_circ=4; % number of circles/2
params.n_ang=12; % number of sectors per circle/2
params.r_scaling=2; % length scaling of sectors.
params.circle_param = 8; %how many circles are included in one item  
params.angle_param = 6; % how many thetas are included in one item
params.numSegments = (2*params.n_circ/params.circle_param)*(2*params.n_ang/params.angle_param);

params.stimrange=255; %221 orig
params.visAngleMin = 1; %enter in degrees will convert to radians later.
params.visAngleMax = 10; 

% define parameters related to task
params.numStimTrial = 2; 
params.numTrials = 20; 
params.nStim = 10;

params.phase = repmat([1;-1],1,params.numTrials);
params.phase = Shuffle(params.phase);

% MIGHT WANT TO TAKE A LOOK AT IMPROVING SIGNAL WITH
                                  % A GOOD WAY OF SWITCHING BETWEEN CUE 1
                                  % AND 2
% select the probes at random without replacement from the x number of possible parts
% the number of areas probed and their occurances should be the same for both cue 1 and cue 2
% shall we probe each part once?  
% I could also have the same number of draws from circle 1 and circle 2


% define timing parameters
time.TR = 1.5;
time.scannerWaitTime = time.TR*1; 
time.stimulus = 1.2;
time.ISI = 0.5; 
time.cue = 1;
time.delay = 5; 
time.probe = 0.5; 
time.answer = 5; 
time.feedback = 1;
time.ITI = 5;


PalRGB=repmat(linspace(params.luminance/2-params.luminance*params.max_con/2,params.luminance/2+params.luminance*params.max_con/2,params.stimrange)',1,3);
bgcol=ceil(params.stimrange/2);
PalRGB(params.stimrange+1,:)=params.luminance;
conrange=(params.stimrange-1)/2;
PalContrasts=round(params.contrasts/params.max_con*conrange);
PalContrasts = PalContrasts';
PalContrasts_probe = round(params.contrast_probe/params.max_con*conrange);

params.mean_con=round(params.mean_con/params.max_con*conrange);

[ch_board, indic, ch_board_parts, ch_board_min]=random_checker_sin(params.n_circ,params.n_ang,params.r_scaling,params.scale, params.circle_param, params.angle_param, rect, params.visAngleMin, params.visAngleMax, params.distFromScreen, Ratio.ppCM);

try
    params.probe = repmat(randsample(1:length(ch_board_parts),params.numTrials/2,false),2,1);  
catch % if cannot sample without replacement then sample with replacement
    params.probe = repmat(randsample(1:length(ch_board_parts),params.numTrials/2,true),2,1);  
end

params.probe(2,:) = Shuffle(params.probe(2,:));
cue1 = 1;
cue2 = 1;
params.probeFull = [];
for asg = 1:length(params.cue)
    if params.cue(asg) ==1
        params.probeFull = [params.probeFull, params.probe(1,cue1)];
        cue1 = cue1+1;
    elseif params.cue(asg) ==2
        params.probeFull = [params.probeFull, params.probe(2,cue2)];
        cue2 = cue2+1;
    end 
end

contrasts_stim = [0:1/(length(params.contrasts)-1):1];

Stimulus_texture_1 = zeros(1,params.numTrials);
Stimulus_texture_2 = zeros(1,params.numTrials);
Stimulus_texture_probe = zeros(1,length(indic));

Stimulus_matrix_full_1 = {};
Stimulus_matrix_full_2 = {};
Contrasts_matrix_stim1 = {};
Contrasts_matrix_stim2 = {};


%Each stimulus has the same TOTAL contrast level
%The number of different individual contrast level items differs across
%trials 

% should each stimulus have the same overall luminance on each given trial?
% probably yes... 

% should the delta also be controled? StimA(x)-StimA(y) = Stim(B)(x)-StimB(y)?

combos = nchoosek(contrasts_stim,length(indic));
% Unique values
[~,idxu,idxc] = unique(sum(combos,2));
combos_all = {};
for c = 1:max(idxc)
    combos_all{c}=combos(idxc==c,:);
end

idx = cellfun(@(x) size(x,1)>1,combos_all);
combos_all = combos_all(idx);

combos_allocations = cellfun(@(x) nchoosek(1:size(x,1),2),combos_all,'UniformOutput',false);

contrast_vector = cellfun(@(x) sum(size(x,1)), combos_allocations);
contrast_vector = repelem(1:length(contrast_vector),contrast_vector); 

contrast_vector_i = [];
for u = unique(contrast_vector)
    [~,~,count] = unique(contrast_vector);
    contrast_vector_i = [contrast_vector_i, 1:sum(count==u)];
end
contrast_vector = [contrast_vector;contrast_vector_i];
contrast_vector = contrast_vector(:,randperm(size(contrast_vector,2))); 
select_contrasts = contrast_vector(:,1:params.numTrials); 


for stim = 1:params.numTrials
   
    trialStim_choice = select_contrasts(:,stim);
    
    contrasts_full1 = combos_all{trialStim_choice(1)}(combos_allocations{trialStim_choice(1)}(trialStim_choice(2),1),:);
    contrasts_full2 = combos_all{trialStim_choice(1)}(combos_allocations{trialStim_choice(1)}(trialStim_choice(2),2),:);
    
    contrasts_full1 = contrasts_full1(randperm(length(contrasts_full1)));
    Contrasts_matrix_stim1{stim} = contrasts_full1; 
    
    
    while 1 % make sure that contrast of stim1(A) ~= stim2(A).
        contrasts_full2 = contrasts_full2(randperm(length(contrasts_full2)));
        if contrasts_full2 ~= contrasts_full1
            Contrasts_matrix_stim2{stim} = contrasts_full2; 
            break
        end
    end
    Contrasts_matrix_stim2{stim} = contrasts_full2; 
    
    contrasts_full1 = contrasts_full1';
    contrasts_full2 = contrasts_full2'; 
    
    
    con_im_1=ones(rect(3),rect(4));
    con_im_2=ones(rect(3),rect(4));
    
    
    for ll=1:params.numSegments
        indic_now = indic{ll};
        for rc = 1:length(indic_now)
            con_im_1(indic_now(rc,2),indic_now(rc,1))=round(contrasts_full1(ll)*(length(params.contrasts)-1))+1;
            con_im_2(indic_now(rc,2),indic_now(rc,1))=round(contrasts_full2(ll)*(length(params.contrasts)-1))+1;

        end 
    end
    
    

    con_im_1 = con_im_1';
    con_im_2 = con_im_2';


    Stimulus_matrix_1=PalContrasts(con_im_1).*ch_board;
    Stimulus_matrix_1(Stimulus_matrix_1~=0) = Stimulus_matrix_1(Stimulus_matrix_1~=0)+params.stimrange/2;
    Stimulus_matrix_1 = Stimulus_matrix_1/params.stimrange;
    if params.phase(1,stim)==-1
        Stimulus_matrix_1(Stimulus_matrix_1~=0) = abs(Stimulus_matrix_1(Stimulus_matrix_1~=0)-1);
    end
    Stimulus_matrix_1(Stimulus_matrix_1==0) = Stimulus_matrix_1(Stimulus_matrix_1==0)+background;
    Stimulus_matrix_full_1{stim} = Stimulus_matrix_1; 
    
    Stimulus_matrix_2=PalContrasts(con_im_2).*ch_board;
    Stimulus_matrix_2(Stimulus_matrix_2~=0) = Stimulus_matrix_2(Stimulus_matrix_2~=0)+params.stimrange/2;
    Stimulus_matrix_2 = Stimulus_matrix_2/params.stimrange;
    if params.phase(2,stim)==-1
        Stimulus_matrix_2(Stimulus_matrix_2~=0) = abs(Stimulus_matrix_2(Stimulus_matrix_2~=0)-1);
    end
    Stimulus_matrix_2(Stimulus_matrix_2==0) = Stimulus_matrix_2(Stimulus_matrix_2==0)+background;
    Stimulus_matrix_full_2{stim} = Stimulus_matrix_2; 
    
    Stimulus_texture_1(stim)=Screen('MakeTexture', window, Stimulus_matrix_1);
    Stimulus_texture_2(stim)=Screen('MakeTexture', window, Stimulus_matrix_2);
    
end

% Assignment of correct responses
cue1 = 1;
cue2 = 1;
params.answerFull = [];
for asg2 = 1:length(params.cue)
    if params.cue(asg2) ==1
        params.answerFull = [params.answerFull, Contrasts_matrix_stim1{asg2}(params.probe(1,cue1))];
        cue1 = cue1+1;
    elseif params.cue(asg2) ==2
        params.answerFull = [params.answerFull, Contrasts_matrix_stim2{asg2}(params.probe(2,cue2))];
        cue2 = cue2+1;
    end 
end

params.answerFull = round(params.answerFull*3)+1;
params.answerFull = PalContrasts(params.answerFull);



Probe_stim_full = {};
Probe_stim_texture = zeros(1,length(PalContrasts_probe)); 
for posContrast = 1:length(PalContrasts_probe)
    Probe_stim=PalContrasts_probe(posContrast).*ch_board_min;
    Probe_stim(Probe_stim~=0) = Probe_stim(Probe_stim~=0)+params.stimrange/2;
    Probe_stim = Probe_stim/params.stimrange;
    Probe_stim(Probe_stim==0) = Probe_stim(Probe_stim==0)+background;
    Probe_stim_full{posContrast} = Probe_stim;
    Probe_stim_texture(posContrast)=Screen('MakeTexture', window, Probe_stim);
end

for part = 1:length(indic)
    ch_board_parts_now = ch_board_parts{part};
    ch_board_parts_now(ch_board_parts_now==0) = background; 
    ch_board_parts_now(ch_board_parts_now==1) = 0;
    Stimulus_texture_probe(part) = Screen('MakeTexture',window,ch_board_parts_now); 
end

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


%% Wait for mouse click to begin
clicks = 0;
Screen('FillRect',window,background);
instructions = 'Instructions row 1!\n\nInstructions row two\n\nRun length: X minutes';
Screen('DrawText',window,'[to Experimenter] waiting for mouse click...', 100, 100);
DrawFormattedText(window, instructions, 'center', 'center');
Screen('Flip',window);
while clicks < 1
    clicks = GetClicks(window);
end


%% Wait to detect first scanner pulse before starting experiment
if mriPulse == 1
    error('This part needs to be implemented  - depends on the scanner')
    Screen('DrawText',window,'[to Experimenter] waiting for mri pulse...', 100, 100);
    Screen('Flip',window);
end

%% MAIN TRIAL LOOP

% Begin presenting stimuli
% Start with fixation cross
expectedTime = time.scannerWaitTime;
draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
[~,startTime] = Screen('Flip',window); % Store time experiment started
Screen('Flip', window, startTime + expectedTime - halfifi);
presentation_timings = {};
presentation_timings.timeExpected = zeros(params.numTrials*params.nStim,1);
presentation_timings.timeActual = zeros(params.numTrials*params.nStim,1);
presentation_timings.type = strings(params.numTrials*params.nStim,1);
presentation_timings.trialType = strings(params.numTrials*params.nStim,1);

expectedKeys = {'a','l','space'};


% loop over blocks
KbQueueCreate;
KbQueueStart;

data.answer = zeros(params.numTrials,1)'; 
data.correctAnswer = params.answerFull';
data.numClick = zeros(params.numTrials,1)';
data.cue = params.cue;


for trialNum = 1:params.numTrials
    % Present stimulus 1
    Screen('DrawTexture', window, Stimulus_texture_1(trialNum));
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+1).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+1).timeActual = presentation_timings((trialNum-1)*params.nStim+1).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+1).type = 'Stim1';
    presentation_timings((trialNum-1)*params.nStim+1).trialType = num2str(params.cue(trialNum));
    presentation_timings((trialNum-1)*params.nStim+1).timeExpected = expectedTime;
    expectedTime = expectedTime+time.stimulus; 

    
    % ISI
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+2).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+2).timeActual = presentation_timings((trialNum-1)*params.nStim+2).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+2).type = 'ISI';
    presentation_timings((trialNum-1)*params.nStim+2).trialType = num2str(params.cue(trialNum));   
    presentation_timings((trialNum-1)*params.nStim+2).timeExpected = expectedTime;
    expectedTime = expectedTime+time.ISI; 

    % Present stimulus 2
    Screen('DrawTexture', window, Stimulus_texture_2(trialNum));
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+3).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+3).timeActual = presentation_timings((trialNum-1)*params.nStim+3).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+3).type = 'Stim2';  
    presentation_timings((trialNum-1)*params.nStim+3).trialType = num2str(params.cue(trialNum));       
    presentation_timings((trialNum-1)*params.nStim+3).timeExpected = expectedTime;
    expectedTime = expectedTime+time.stimulus; 
    
    % ISI
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+4).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+4).timeActual = presentation_timings((trialNum-1)*params.nStim+4).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+4).type = 'ISI';
    presentation_timings((trialNum-1)*params.nStim+4).trialType = num2str(params.cue(trialNum));      
    presentation_timings((trialNum-1)*params.nStim+4).timeExpected = expectedTime;
    expectedTime = expectedTime+time.ISI; 
    
    % Cue
    DrawFormattedText(window, num2str(params.cue(trialNum)) , 'center', 'center');
    [~, presentation_timings((trialNum-1)*params.nStim+5).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+5).timeActual = presentation_timings((trialNum-1)*params.nStim+5).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+5).type = 'Cue';
    presentation_timings((trialNum-1)*params.nStim+5).trialType = num2str(params.cue(trialNum));         
    presentation_timings((trialNum-1)*params.nStim+5).timeExpected = expectedTime;
    expectedTime = expectedTime+time.cue; 
    
    % Delay
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+6).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+6).timeActual = presentation_timings((trialNum-1)*params.nStim+6).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+6).type = 'Delay';
    presentation_timings((trialNum-1)*params.nStim+6).trialType = num2str(params.cue(trialNum));             
    presentation_timings((trialNum-1)*params.nStim+6).timeExpected = expectedTime;
    expectedTime = expectedTime+time.delay; 
    
    % Probe
    Screen('DrawTexture', window, Stimulus_texture_probe(params.probeFull(trialNum))); 
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+7).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+7).timeActual = presentation_timings((trialNum-1)*params.nStim+7).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+7).type = 'Probe'; 
    presentation_timings((trialNum-1)*params.nStim+7).trialType = num2str(params.cue(trialNum));                 
    presentation_timings((trialNum-1)*params.nStim+7).timeExpected = expectedTime;
    expectedTime = expectedTime+time.probe; 
    
    close=0;
    % Answer   
    contrast_number = randperm(length(Probe_stim_texture),1);
    Screen('DrawTexture', window, Probe_stim_texture(contrast_number)); 
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+8).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+8).timeActual = presentation_timings((trialNum-1)*params.nStim+8).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+8).type = 'Answer'; 
    presentation_timings((trialNum-1)*params.nStim+8).trialType = num2str(params.cue(trialNum));                    
    presentation_timings((trialNum-1)*params.nStim+8).timeExpected = expectedTime;
    
    %  
    response = [];
    close = 0;
    KbQueueFlush;
    
    while strcmp(response,'m')==0 
        response = [];
        if GetSecs > startTime+expectedTime+time.answer-halfifi
            response = 'm'; %m for miss
            break
        else
            response = [];
            loophere = 1; 
            [~,keyCode]=KbQueueCheck;            
            if keyCode == 0
                continue
            else
                response = KbName(keyCode);
                % Check if the key pressed was one we expected
                rep = find(strcmp(response,expectedKeys) == 1, 1);
                % If the response was not one of the keys we were expecting, then make
                % response empty again and wait for another response. 
                if isempty(rep)
                    response = [];
                end
            end
            
            if response == 'a'
                if contrast_number == 1
                    contrast_number = length(Probe_stim_texture);
                else
                    contrast_number = contrast_number-1;
                end
                data.numClick(trialNum)=data.numClick(trialNum)+1; 
            elseif response == 'l'
                if contrast_number == length(Probe_stim_texture)
                    contrast_number = 1;
                else
                    contrast_number = contrast_number+1;
                end
                data.numClick(trialNum)=data.numClick(trialNum)+1;
            elseif strcmp(response,'space')
                close=1;
                break
            end
            
%             rt = secs - startTime;
        end
        Screen('DrawTexture', window, Probe_stim_texture(contrast_number)); 
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);    
        Screen('Flip', window, startTime + expectedTime - halfifi);       
    end
    if close==1
        Screen('DrawTexture', window, Probe_stim_texture(contrast_number));
        data.answer(trialNum) = PalContrasts_probe(contrast_number); 
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,params.distFromScreen,Ratio.ppCM);    
        Screen('Flip', window, startTime + expectedTime - halfifi);      
    end
    expectedTime = expectedTime+time.answer; 
    
    % Feedback 
    if strcmp(response,'m')
        displayFeedback = 'Miss';
        data.answer(trialNum) = 999; 
    elseif data.answer(trialNum) == params.answerFull(trialNum)
        displayFeedback = 'Correct';
    elseif data.answer(trialNum) ~= params.answerFull(trialNum)
        displayFeedback = 'Incorrect';
    end
    
    DrawFormattedText(window, displayFeedback , 'center', 'center');
    [~, presentation_timings((trialNum-1)*params.nStim+9).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+9).timeActual = presentation_timings((trialNum-1)*params.nStim+9).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+9).type = 'Feedback';
    presentation_timings((trialNum-1)*params.nStim+9).trialType = num2str(params.cue(trialNum));                       
    presentation_timings((trialNum-1)*params.nStim+9).timeExpected = expectedTime;
    expectedTime = expectedTime+time.feedback; 
    
    draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
    [~, presentation_timings((trialNum-1)*params.nStim+10).timeActual] = Screen('Flip', window, startTime + expectedTime - halfifi);
    presentation_timings((trialNum-1)*params.nStim+10).timeActual = presentation_timings((trialNum-1)*params.nStim+10).timeActual - startTime;
    presentation_timings((trialNum-1)*params.nStim+10).type = 'ITI';
    presentation_timings((trialNum-1)*params.nStim+10).trialType = num2str(params.cue(trialNum));                          
    presentation_timings((trialNum-1)*params.nStim+10).timeExpected = expectedTime;
    expectedTime = expectedTime+time.ITI;   
        
end

sca;
save(outFile, 'data', 'params', 'time', 'presentation_timings');




% Function to give the subject a short break ;)
% %     function takeBreak
% %         DrawFormattedText(window,'Feel free to take a short break\n Click the mouse to continue','center', 'center');
% %         Screen('Flip',window);
% %         
% %         clicks = 0;
% %         % Wait for click
% %         while clicks < 1
% %             clicks = GetClicks(window);
% %         end
% %                         
% %         % Draw fixation and wait for a couple of seconds
% %         draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
% %         Screen('DrawingFinished', window);
% %         % Reset timings since we may have waited a while
% %         startTime = Screen('Flip', window);
% %         expectedTime = 0;
% %           
% %         expectedTime = expectedTime + 2;
% %         draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
% %         Screen('DrawingFinished', window);
% %         Screen('Flip', window, startTime + expectedTime - halfifi);
% %     end

