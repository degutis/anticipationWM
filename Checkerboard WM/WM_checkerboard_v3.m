% error('more things need to be saved. Do that b4 running');   
subjID = input('Enter subject ID: \n','s'); %e.g. 001 002 010 099
WM_run = input('Enter run number: \n'); %e.g. 1 2 3

environment = 'behav';

fileID = [subjID,'_',num2str(WM_run)];

%% Subj ID & Paths

fileName = [fileID, '_WM_retro'];
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
params.numBlocksDiff = 2;
params.numBlocksPer = 3;
params.Blocks=[];
for k=1:params.numBlocksPer
    params.Blocks=[params.Blocks,randperm(params.numBlocksDiff)];
end

% params.Blocks = params.Blocks(randperm(length(params.Blocks))); 
params.Blocks_names = {'Two stimuli','One stimulus'}; 

params.numTrials = 32; % trials per block
params.nStim = [10,6]; %10 stim in block"1" and 6 in block"2" 

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
params.angle_param = 6; % how many thetas are included in one item
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


params.phase = {}; 
params.cue = cell(1,length(params.Blocks)); 

for bl0 = 1:length(params.Blocks)
    if params.Blocks(bl0) == 1
        params.phase{bl0} = repmat([1;-1],1,params.numTrials);
        params.phase{bl0} = Shuffle(params.phase{bl0});
        params.cue{bl0} = repmat([1,2],1,params.numTrials/2); 
        params.cue{bl0} = Shuffle(params.cue{bl0}); 

    elseif params.Blocks(bl0) == 2
        params.phase{bl0} = repmat([1,-1],1,params.numTrials/2);
        params.phase{bl0} = params.phase{bl0}(randperm(length(params.phase{bl0}))); 
    end
end


params.numSegments


[ch_board, indic, ch_board_parts, ch_board_min]=random_checker_sin(params.n_circ,params.n_ang,params.r_scaling,params.scale, params.circle_param, params.angle_param, rect, params.visAngleMin, params.visAngleMax, params.distFromScreen, Ratio.ppCM);

params.probe = {};
for bl = 1:length(params.Blocks)
    if params.Blocks(bl)==1
        try
            params.probe{bl} = repmat(randsample(1:length(ch_board_parts),params.numTrials/2,false),2,1);  
        catch % if cannot sample without replacement then sample with replacement
            params.probe{bl} = repmat(1:length(ch_board_parts),2,floor(params.numTrials/length(ch_board_parts)/2))
            params.probe{bl} = [params.probe{bl},repmat(randsample(1:length(ch_board_parts),params.numTrials/2-size(params.probe{bl},2),false),2,1)];
        end
        params.probe{bl} = params.probe{bl}(:,randperm(length(params.probe{bl}(1,:))));
        params.probe{bl}(2,:) = params.probe{bl}(2,randperm(length(params.probe{bl}(2,:))));
        
    elseif params.Blocks(bl)==2
        try
            params.probe{bl} = repmat(randsample(1:length(ch_board_parts),params.numTrials,false),1,1);  
        catch % if cannot sample without replacement then sample with replacement
            params.probe{bl} = repmat(1:length(ch_board_parts),1,floor(params.numTrials/length(ch_board_parts)))
            params.probe{bl} = [params.probe{bl},(randsample(1:length(ch_board_parts),params.numTrials-size(params.probe{bl},2),false))];
        end        
        params.probe{bl} = params.probe{bl}(:,randperm(length(params.probe{bl}(1,:))));
    end
end
        

params.probeFull = cell(1,length(params.Blocks));
for bl2 = 1:length(params.Blocks)
    cue1 = 1;
    cue2 = 1;
    if params.Blocks(bl2) ==1
        for asg = 1:length(params.cue{bl2})
            if params.cue{bl2}(asg) ==1
                params.probeFull{bl2} = [params.probeFull{bl2}, params.probe{bl2}(1,cue1)];
                cue1 = cue1+1;
            elseif params.cue{bl2}(asg) ==2
                params.probeFull{bl2} = [params.probeFull{bl2}, params.probe{bl2}(2,cue2)];
                cue2 = cue2+1;
            end
        end
    elseif params.Blocks(bl2) ==2
        params.probeFull{bl2} = params.probe{bl2};
    end
end

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


%% Generate stimuli 
contrasts_stim = [0:1/(length(params.contrasts)-1):1];

Stimulus_texture_1 = cell(1,length(params.Blocks));
Stimulus_texture_2 = cell(1,length(params.Blocks));

% Stimulus_matrix_full_1 = cell(1,cell(1:params.numTrials));
% Stimulus_matrix_full_2 = cell(1,length(params.Blocks));
Contrasts_matrix_stim1 = cell(1,length(params.Blocks));
Contrasts_matrix_stim2 = cell(1,length(params.Blocks));


%Each stimulus has the same TOTAL contrast level
%The number of different individual contrast level items differs across
%trials 

% should each stimulus have the same overall luminance on each given trial?
% probably yes... 

% should the delta also be controled? StimA(x)-StimA(y) = Stim(B)(x)-StimB(y)?
if length(indic)<length(contrasts_stim)
    combos = nchoosek(contrasts_stim,length(indic));
elseif length(indic)>length(contrasts_stim)
    combos = nchoosek(repelem(contrasts_stim,2),length(indic));

%     combos = [];
%     for conextra = 1:length(contrasts_stim)
%         combos = [combos;[contrasts_stim,contrasts_stim(conextra)]];
%     end
end
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

select_contrasts = cell(1,length(params.Blocks));
for bl3 = 1:length(params.Blocks)
    
    contrast_vector_randperm = contrast_vector(:,randperm(size(contrast_vector,2))); 
    select_contrasts{bl3} = contrast_vector_randperm(:,1:params.numTrials); 

    for stim = 1:params.numTrials

        trialStim_choice = select_contrasts{bl3}(:,stim);

        contrasts_full1 = combos_all{trialStim_choice(1)}(combos_allocations{trialStim_choice(1)}(trialStim_choice(2),1),:);
        contrasts_full2 = combos_all{trialStim_choice(1)}(combos_allocations{trialStim_choice(1)}(trialStim_choice(2),2),:);

        contrasts_full1 = contrasts_full1(randperm(length(contrasts_full1)));
        Contrasts_matrix_stim1{bl3}{stim} = contrasts_full1; 

        if params.Blocks(bl3)==1
            while 1 % make sure that contrast of stim1(A) ~= stim2(A).
                contrasts_full2 = contrasts_full2(randperm(length(contrasts_full2)));
                if contrasts_full2 ~= contrasts_full1
                    Contrasts_matrix_stim2{bl3}{stim} = contrasts_full2; 
                    break
                end
            end
%             Contrasts_matrix_stim2{bl}{stim} = contrasts_full2; 
        end

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
        if params.phase{bl3}(1,stim)==-1
            Stimulus_matrix_1(Stimulus_matrix_1~=0) = abs(Stimulus_matrix_1(Stimulus_matrix_1~=0)-1);
        end
        Stimulus_matrix_1(Stimulus_matrix_1==0) = Stimulus_matrix_1(Stimulus_matrix_1==0)+background;
%         Stimulus_matrix_full_1{} = Stimulus_matrix_1; 

        if params.Blocks(bl3)==1

            Stimulus_matrix_2=PalContrasts(con_im_2).*ch_board;
            Stimulus_matrix_2(Stimulus_matrix_2~=0) = Stimulus_matrix_2(Stimulus_matrix_2~=0)+params.stimrange/2;
            Stimulus_matrix_2 = Stimulus_matrix_2/params.stimrange;
            if params.phase{bl3}(2,stim)==-1
                Stimulus_matrix_2(Stimulus_matrix_2~=0) = abs(Stimulus_matrix_2(Stimulus_matrix_2~=0)-1);
            end
            Stimulus_matrix_2(Stimulus_matrix_2==0) = Stimulus_matrix_2(Stimulus_matrix_2==0)+background;
%             Stimulus_matrix_full_2{} = Stimulus_matrix_2; 
        end
        Stimulus_texture_1{bl3}(stim)=Screen('MakeTexture', window, Stimulus_matrix_1);
        if params.Blocks(bl3)==1
            Stimulus_texture_2{bl3}(stim)=Screen('MakeTexture', window, Stimulus_matrix_2);
        end
    end
end


params.answerFull = cell(1,length(params.Blocks));
   
% Assignment of correct responses
for bl4 = 1:length(params.Blocks)

    cue1 = 1;
    cue2 = 1;
    if params.Blocks(bl4)==1
        for asg2 = 1:length(params.cue{bl4})
            if params.cue{bl4}(asg2) ==1
                params.answerFull{bl4} = [params.answerFull{bl4}, Contrasts_matrix_stim1{bl4}{asg2}(params.probe{bl4}(1,cue1))];
                cue1 = cue1+1;
            elseif params.cue{bl4}(asg2) ==2
                params.answerFull{bl4} = [params.answerFull{bl4}, Contrasts_matrix_stim2{bl4}{asg2}(params.probe{bl4}(2,cue2))];
                cue2 = cue2+1;
            end 
        end

    elseif params.Blocks(bl4) ==2
        for asg4 = 1:params.numTrials
            params.answerFull{bl4} = [params.answerFull{bl4},Contrasts_matrix_stim1{bl4}{asg4}(params.probe{bl4}(asg4))];
        end
    end        
    params.answerFull{bl4} = round(params.answerFull{bl4}*(length(params.contrasts)-1))+1;
    params.answerFull{bl4} = PalContrasts(params.answerFull{bl4});
    
end


% Generate all possible contrasts for the answer stimuli
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

% Generate the probe stimuli
Stimulus_texture_probe = zeros(1,length(indic));
for part = 1:length(indic)
    ch_board_parts_now = ch_board_parts{part};
    ch_board_parts_now(ch_board_parts_now==0) = background; 
    ch_board_parts_now(ch_board_parts_now==1) = 0;
    Stimulus_texture_probe(part) = Screen('MakeTexture',window,ch_board_parts_now); 
end



%% Wait for mouse click to begin
clicks = 0;
Screen('FillRect',window,background);
instructions = ['Instructions row 1!\n\nInstructions row two\n\nRun length: X minutes\n\n First block is: ',params.Blocks_names{params.Blocks(1)}];
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

data = cell(1,length(params.Blocks));
presentation_timings = cell(1,length(params.Blocks));

for blockz = 1:length(params.Blocks)
        data{blockz}.answer = zeros(params.numTrials,1)'; 
        data{blockz}.correctAnswer = params.answerFull{blockz}';
        data{blockz}.numClick = zeros(params.numTrials,1)';
        data{blockz}.cue = params.cue{blockz};
        if params.Blocks(blockz)==1
            presentation_timings{blockz}.timeExpected = zeros(params.numTrials*params.nStim(1),1);
            presentation_timings{blockz}.timeActual = zeros(params.numTrials*params.nStim(1),1);
            presentation_timings{blockz}.type = strings(params.numTrials*params.nStim(1),1);
            presentation_timings{blockz}.trialType = strings(params.numTrials*params.nStim(1),1);
        elseif params.Blocks(blockz)==2
            presentation_timings{blockz}.timeExpected = zeros(params.numTrials*params.nStim(2),1);
            presentation_timings{blockz}.timeActual = zeros(params.numTrials*params.nStim(2),1);
            presentation_timings{blockz}.type = strings(params.numTrials*params.nStim(2),1);
            presentation_timings{blockz}.trialType = strings(params.numTrials*params.nStim(2),1);
        end
end


% Begin presenting stimuli
% Start with fixation cross
expectedTime = time.scannerWaitTime;
draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
[~,startTime] = Screen('Flip',window); % Store time experiment started
Screen('Flip', window, startTime + expectedTime - halfifi);

expectedKeys = {'a','l','space'};


KbQueueCreate;
KbQueueStart;


for block = 1:length(params.Blocks)
    if block~=1
        % Function to give the subject a short break ;)
        DrawFormattedText(window,['Feel free to take a short break\n Click the mouse to continue\n Next block is: ',params.Blocks_names{params.Blocks(block)}],'center', 'center');
        Screen('Flip',window);

        clicks = 0;
        % Wait for click
        while clicks < 1
            clicks = GetClicks(window);
        end

        % Draw fixation and wait for a couple of seconds
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,params.distFromScreen,Ratio.ppCM);    
        Screen('DrawingFinished', window);
        % Reset timings since we may have waited a while
        startTime = Screen('Flip', window);
        expectedTime = 0;

        expectedTime = expectedTime + time.scannerWaitTime;
        draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,params.distFromScreen,Ratio.ppCM);    
        Screen('DrawingFinished', window);
        Screen('Flip', window, startTime + expectedTime - halfifi);
    end
    if params.Blocks(block) == 1
        for trialNum = 1:params.numTrials
            % Present stimulus 1
            Screen('DrawTexture', window, Stimulus_texture_1{block}(trialNum));
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+1)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+1) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+1) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+1) = 'Stim1';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+1) = num2str(params.cue{block}(trialNum));
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+1) = expectedTime;
            expectedTime = expectedTime+time.stimulus; 


            % ISI
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+2)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+2) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+2) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+2)= 'ISI';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+2) = num2str(params.cue{block}(trialNum));   
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+2) = expectedTime;
            expectedTime = expectedTime+time.ISI; 

            % Present stimulus 2
            Screen('DrawTexture', window, Stimulus_texture_2{block}(trialNum));
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+3)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+3) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+3) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+3) = 'Stim2';  
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+3) = num2str(params.cue{block}(trialNum));       
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+3) = expectedTime;
            expectedTime = expectedTime+time.stimulus; 

            % ISI
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+4)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+4) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+4) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+4) = 'ISI';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+4) = num2str(params.cue{block}(trialNum));      
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+4) = expectedTime;
            expectedTime = expectedTime+time.ISI; 

            % Cue
            DrawFormattedText(window, num2str(params.cue{block}(trialNum)) , 'center', 'center');
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+5)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+5) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+5) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+5) = 'Cue';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+5) = num2str(params.cue{block}(trialNum));         
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+5) = expectedTime;
            expectedTime = expectedTime+time.cue; 

            % Delay
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+6)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+6) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+6) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+6) = 'Delay';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+6) = num2str(params.cue{block}(trialNum));             
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+6) = expectedTime;
            expectedTime = expectedTime+time.delay; 

            % Probe
            Screen('DrawTexture', window, Stimulus_texture_probe(params.probeFull{block}(trialNum))); 
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+7)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+7) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+7) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+7) = 'Probe'; 
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+7) = num2str(params.cue{block}(trialNum));                 
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+7) = expectedTime;
            expectedTime = expectedTime+time.probe; 

            % Answer   
            contrast_number = randperm(length(Probe_stim_texture),1);
            Screen('DrawTexture', window, Probe_stim_texture(contrast_number)); 
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+8)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+8) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+8) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+8) = 'Answer'; 
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+8) = num2str(params.cue{block}(trialNum));                    
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+8) = expectedTime;


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
                        data{block}.numClick(trialNum)=data{block}.numClick(trialNum)+1;
                    elseif response == 'l'
                        if contrast_number == length(Probe_stim_texture)
                            contrast_number = 1;
                        else
                            contrast_number = contrast_number+1;
                        end
                        data{block}.numClick(trialNum)=data{block}.numClick(trialNum)+1;
                    elseif strcmp(response,'space')
                        close=1;
                        break
                    end
                end
                Screen('DrawTexture', window, Probe_stim_texture(contrast_number));
                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
                Screen('Flip', window, startTime + expectedTime - halfifi);
            end
            if close==1
                Screen('DrawTexture', window, Probe_stim_texture(contrast_number));
                data{block}.answer(trialNum) = PalContrasts_probe(contrast_number);
                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,params.distFromScreen,Ratio.ppCM);
                Screen('Flip', window, startTime + expectedTime - halfifi);
            end
            expectedTime = expectedTime+time.answer;            


            % Feedback 
            if strcmp(response,'m')
                displayFeedback = 'Miss';
                data{block}.answer(trialNum) = 999; 
            elseif data{block}.answer(trialNum) == params.answerFull{block}(trialNum)
                displayFeedback = 'Correct';
            elseif data{block}.answer(trialNum) ~= params.answerFull{block}(trialNum)
                displayFeedback = 'Incorrect';
            end

            DrawFormattedText(window, displayFeedback , 'center', 'center');
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+9)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+9) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+9) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+9) = 'Feedback';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+9) = num2str(params.cue{block}(trialNum));                       
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+9) = expectedTime;
            expectedTime = expectedTime+time.feedback; 

            % ITI
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+10)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+10) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(1)+10) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(1)+10) = 'ITI';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(1)+10) = num2str(params.cue{block}(trialNum));                          
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(1)+10) = expectedTime;
            expectedTime = expectedTime+time.ITI;   
        end

    elseif params.Blocks(block) == 2
        for trialNum = 1:params.numTrials
            % Present stimulus 1
            Screen('DrawTexture', window, Stimulus_texture_1{block}(trialNum));
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+1)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+1) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+1) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(2)+1) = 'Stim1';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(2)+1) = 'oneStim';
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(2)+1) = expectedTime;
            expectedTime = expectedTime+time.stimulus; 


            % Delay
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+2)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+2) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+2) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(2)+2) = 'Delay';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(2)+2) = 'oneStim';             
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(2)+2) = expectedTime;
            expectedTime = expectedTime+time.delayOneStim; 

            % Probe
            Screen('DrawTexture', window, Stimulus_texture_probe(params.probeFull{block}(trialNum))); 
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+3)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+3) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+3) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(2)+3) = 'Probe'; 
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(2)+3) = 'oneStim';                 
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(2)+3) = expectedTime;
            expectedTime = expectedTime+time.probe; 

            % Answer   
            contrast_number = randperm(length(Probe_stim_texture),1);
            Screen('DrawTexture', window, Probe_stim_texture(contrast_number)); 
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+4)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+4) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+4) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(2)+4) = 'Answer'; 
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(2)+4) = 'oneStim';                    
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(2)+4) = expectedTime;
            
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
                        data{block}.numClick(trialNum)=data{block}.numClick(trialNum)+1;
                    elseif response == 'l'
                        if contrast_number == length(Probe_stim_texture)
                            contrast_number = 1;
                        else
                            contrast_number = contrast_number+1;
                        end
                        data{block}.numClick(trialNum)=data{block}.numClick(trialNum)+1;
                    elseif strcmp(response,'space')
                        close=1;
                        break
                    end
                end
                Screen('DrawTexture', window, Probe_stim_texture(contrast_number));
                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
                Screen('Flip', window, startTime + expectedTime - halfifi);
            end
            if close==1
                Screen('DrawTexture', window, Probe_stim_texture(contrast_number));
                data{block}.answer(trialNum) = PalContrasts_probe(contrast_number);
                draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],255,255,params.distFromScreen,Ratio.ppCM);
                Screen('Flip', window, startTime + expectedTime - halfifi);
            end
            expectedTime = expectedTime+time.answer;

            % Feedback 
            if strcmp(response,'m')
                displayFeedback = 'Miss';
                data{block}.answer(trialNum) = 999; 
            elseif data{block}.answer(trialNum) == params.answerFull{block}(trialNum)
                displayFeedback = 'Correct';
            elseif data{block}.answer(trialNum) ~= params.answerFull{block}(trialNum)
                displayFeedback = 'Incorrect';
            end

            DrawFormattedText(window, displayFeedback , 'center', 'center');
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+5)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+5) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+5) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(2)+5) = 'Feedback';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(2)+5) = 'oneStim';                      
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(2)+5) = expectedTime;
            expectedTime = expectedTime+time.feedback; 

            % ITI
            draw_fixationdot(window,background, rect(3), rect(4),[0.25 0.05],0,0,params.distFromScreen,Ratio.ppCM);
            [~, presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+6)] = Screen('Flip', window, startTime + expectedTime - halfifi);
            presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+6) = presentation_timings{block}.timeActual((trialNum-1)*params.nStim(2)+6) - startTime;
            presentation_timings{block}.type((trialNum-1)*params.nStim(2)+6) = 'ITI';
            presentation_timings{block}.trialType((trialNum-1)*params.nStim(2)+6) = 'oneStim';                          
            presentation_timings{block}.timeExpected((trialNum-1)*params.nStim(2)+6) = expectedTime;
            expectedTime = expectedTime+time.ITI;   
        end        
    end
end

sca;
save(outFile, 'data', 'params', 'Ratio','time', 'presentation_timings', 'select_contrasts','Contrasts_matrix_stim1','Contrasts_matrix_stim2');





