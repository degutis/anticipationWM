function Psychophysics_V6
% Run psychophysics working memory experiment to establish 75% correct
% threshold for orientation dscrimination
%
% Present two reference stimuli (45 and 135 degrees), cue one to be
% remembered, after delay period present test stimulus with orientation
% near cued reference stimulus, discriminate relative to cued reference
%
% Use QUEST to measure 75% correct threshold over 100 trials
%
% 10/11/2016 SJDL wrote it.
%
% 03/11/2016 V3 edit = using a different function to generate stimuli that
% smoothes the edges of the stimulus.
%
% V4 edit - stimulus orientation now jittered by +/- 3 degrees
%
% V5 edit - gives feedback after each trial
% V6 edit --> added a fixation cross instead of a dot. 
TR = 3.2; % TR will determine stimulus timings

%--------------------------------------------------------------------------
% Add functions to path
addpath('.\..\Functions');
addpath(genpath('.\..\Toolbox'));
% Add java memory cleaner
javaaddpath(which('MatlabGarbageCollector.jar'))

% Reset rng
rng('default');
rng('shuffle');

%--------------------------------------------------------------------------
% Define environment and load appropriate luminance settings
scanner = 'prisma';
environment = 'behav';
background = calibrate_lum(1, environment,scanner);

bitsiFlag = 0;
% Monitor settings for each environment
switch environment
    case 'work_station'
        distFromScreen = 60;
        pixelsPerCm = 33.7;
    case 'dummy'
        distFromScreen = 80;
        pixelsPerCm = 26.7;
    case 'mri'
        distFromScreen = 90;
        pixelsPerCm = 27.74;
    case 'behav'
        distFromScreen = 60;
        pixelsPerCm = 26.3;
        bitsiFlag = 1; % Will use button box for responses
        bitsiResponse = Bitsi_Scanner('com1');
        bitsi_key1 = 97; % Key code
        bitsi_key2 = 98;
        
end

%--------------------------------------------------------------------------
% Get experiment information and set up file for saving
% Experiment information
subjID = input('Enter subject ID\n');
practice = input('Is this a practice session? Enter 1 for yes\n');
if practice == 1
    myFile = sprintf('WM_psychophysics_S%d_Practice',subjID);
else
    SessionID = input('Enter session number\n');
    myFile = sprintf('WM_psychophysics_S%d_Session%d',subjID,SessionID);
end
[~,myDir] = uiputfile(myFile,'Choose file directory');
outFile = [myDir, myFile, '.mat'];

%--------------------------------------------------------------------------
% Set up key experiment parameters
whichScreen = 0; %Screen ID

params = struct();
params.stimSize = 16; % Diameter in degrees
params.refOrient = [45, 135]; % reference stimuli are 45 and 135 degrees
params.sf = 1; % cpd
params.contrast = .5;
params.refPres = .2; % Reference stimulus presentation time
params.ISI = .45; % Interval between reference stimuli
params.cueTime = .8; % Time to present cue for
params.trialLength = TR * 4;
params.delay = params.trialLength - (params.refPres + params.ISI + params.refPres + params.ISI + params.cueTime); % Delay period
params.testPres = .5; % test stimulus presentation time
params.response = 2; % response window
params.feedback = .3; % Feedback presentation
params.ITI = 2.5 - (params.response + params.feedback); % Time between trials
params.phases = linspace(0, 2*pi, 11); % 10 possible stimulus phases, will be selected at random each trial
params.phases(1) = [];
params.dotSize = 4; % Size of fixation dot in pixels
if practice == 1
    params.numTrials = 20;
else
    params.numTrials = 120; % Number of trials we will do, must be divisible by 20
end

%--------------------------------------------------------------------------
% Set up key QUEST staircase parameters
if practice == 1
    tGuess = 10; % Initial threshold estimate
    tGuessSd = 9.9; % Standard deviation for estimate (be generous!)
    range = 20; % Largest difference between max and min intensity that internal table can store
else
    % USED PRACTICE SESSION RESULTS FOR VARIABLES BELOW
    tGuess = 5; % Initial threshold estimate
    tGuessSd = 4.9; % Standard deviation for estimate (be generous!)
    range = 10; % Largest difference between max and min intensity that internal table can store
end
pThreshold = 0.75; % 75% correct threshold
beta = 3.5; % steepnest of psychometric function, 3.5 is typical
delta = 0.01; % fraction of trials observer presses blindly, 0.01 is typical
gamma = 0.5; % fraction of trials that will generate correct response when intensity == -inf (.5 = chance)
grain = 0.01; % Step size of internal table

% Initialise staircase
params.q = QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range);
params.q.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

%--------------------------------------------------------------------------
% Open window
[window, rect] = Screen('OpenWindow',whichScreen);
HideCursor(window);
% Get resolution of screen
x = rect(3);
y = rect(4);
% Change font
Screen('TextFont',window, 'Arial');
Screen('TextSize',window, 20);
% Fill screen with mid-grey background
Screen('FillRect', window, background);

%--------------------------------------------------------------------------
% Generate stimulus textures
Screen('DrawText',window,'Generating textures...', 100, 100);
Screen('Flip',window);

% Will need a different texture for each stimulus phase
stimTex = zeros(1,length(params.phases)); % Preallocate
for i = 1:length(params.phases)
    stim = makeGaborStimulus(params.contrast,params.stimSize,params.phases(i),params.sf,1,distFromScreen,environment,scanner, whichScreen);
    stimTex(i) = Screen('MakeTexture', window, stim);
end

% Preload textures into video memory
Screen('PreloadTextures',window,stimTex);

%--------------------------------------------------------------------------
% Preallocate some variables for storing information
correct = zeros(1, params.numTrials); % Will store if response was correct
testOrientation = correct; % Will store orientation of test stimulus (relative to reference)
if bitsiFlag == 1
    response = correct;
else
    response = cell(1,params.numTrials); % Will restore key presses
end
rt = correct; % Will store reaction times
buttonPr = correct; % Store button presses if using button box

%--------------------------------------------------------------------------
% Determine random variables for this session

% On half trials 45 degree ref stimulus presented first (1), 135 degree first
% for other half (2), in pseudorandom order (no more than 3 of the same in
% a row)
refOrder = Create_pseudo_randomised_list([1, 2], params.numTrials/2);
% Jitters for each reference stimulis
jitters1 = Create_pseudo_randomised_list([3 -3], params.numTrials/2);
jitters2 = Create_pseudo_randomised_list([3 -3], params.numTrials/2);
% Going to pseudorandomise cue so that on 45 degree cued on 50% of trials
% and 135 degree cued on other 50%. This works by taking trials where 45
% degree is presented first (defined in refOrder) and cueing first stimulus
% on half of those trials (25% of total) and the second for the other half,
% then doing the same for trials where 135 degree was presented first.
cue = zeros(size(refOrder)); % First preallocate
cueTemp = Create_pseudo_randomised_list([1, 2], params.numTrials/4); % Create list defining which stimulus will be cued
cue(refOrder==1) = cueTemp; % Sub list into trials where 45 degree presented first
cueTemp = Create_pseudo_randomised_list([1, 2], params.numTrials/4); % Create new list defining which stimulus will be cued
cue(refOrder==2) = cueTemp; % Sub list into trials where 135 degrees was presented first

% And for stimulus phase - 3 stimuli per trial so need 3 phases for each
% trial
phases = zeros(3,params.numTrials);
for i = 1:size(phases,1)
    phases(i,:) = Create_pseudo_randomised_list(1:length(params.phases), params.numTrials/length(params.phases));
end
% Test stimulus should be clockwise of reference (1) for half trials and
% anti clockwise (-1) for other half
testDirection = Create_pseudo_randomised_list([-1 1],params.numTrials/2);

%--------------------------------------------------------------------------
% Measure flip interval for accurate stimulus timings!
Priority(1); % Set priority
Screen('DrawText',window,'Estimating monitor flip interval...', 100, 100);
Screen('DrawText',window,'(This may take up to 20s)', 100, 120);
Screen('Flip',window);
halfifi = Screen('GetFlipInterval', window, 100, 0.00005, 20);
halfifi = halfifi/2;
Priority(0); % reset priority level

%--------------------------------------------------------------------------
% Wait for mouse click to begin
clicks = 0;
Screen('DrawText',window,'Waiting for mouse click...', 100, 100);
Screen('Flip',window);
while clicks < 1
    clicks = GetClicks(window);
end

%--------------------------------------------------------------------------
% MAIN TRIAL LOOP
Priority(1); % Set priority

% Begin presenting stimuli
% Start with fixation cross
draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
startTime = Screen('Flip',window); % Store time experiment started
expectedTime = 0; % This will record what the expected event duration should be

tookBreak = 0; % Will get set to 1 when we take a break

for trialNum = 1:params.numTrials
    
    % wait a couple of seconds on first trial and after break
    if trialNum == 1 || tookBreak == 1
        expectedTime = expectedTime + 2;
        tookBreak = 0; % Reset
    else 
        % Draw first reference stimulus after ITI
        expectedTime = expectedTime + params.ITI;
    end
    
    % Jitters set separately for each trial depending on which reference
    % gets presented first, to ensure equal jitters for each stimulus
    switch refOrder(trialNum)
        case 1
            Screen('DrawTexture', window, stimTex(phases(1,trialNum)), [], [], params.refOrient(refOrder(trialNum))+jitters1(trialNum));
            firstRefJitter = jitters1(trialNum); % Need to store what jitter was so that test stimulus is correct
        case 2
            Screen('DrawTexture', window, stimTex(phases(1,trialNum)), [], [], params.refOrient(refOrder(trialNum))+jitters2(trialNum));
            firstRefJitter = jitters2(trialNum);
    end
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Return to fixation after reference presentation time
    expectedTime = expectedTime + params.refPres;
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Draw second reference stimulus after ISI
    expectedTime = expectedTime + params.ISI;
    switch refOrder(trialNum)
        case 1
            Screen('DrawTexture', window, stimTex(phases(2,trialNum)), [], [], params.refOrient(2)+jitters2(trialNum));
            secondRefJitter = jitters2(trialNum);
        case 2
            Screen('DrawTexture', window, stimTex(phases(2,trialNum)), [], [], params.refOrient(1)+jitters1(trialNum));
            secondRefJitter = jitters1(trialNum);
    end
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Return to fixation after reference presentation time
    expectedTime = expectedTime + params.refPres;
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Present cue after ISI
    expectedTime = expectedTime + params.ISI;
    DrawFormattedText(window, num2str(cue(trialNum)), 'center', 'center');
    Screen('DrawingFinished',window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
        
    % Return to fixation after cue period
    expectedTime = expectedTime + params.cueTime;
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Present test stimulus after delay period
    expectedTime = expectedTime + params.delay;
    % Get test orientation for this trial from quest
    testOrientation(trialNum) = QuestQuantile(params.q);
    % Draw the stimulus
    switch cue(trialNum)
        case 1 % If cued first stimulus then add orientation to that stimulus (stored in refOrder)
            switch testDirection(trialNum) % Find out if going clockwise or anticlockwise
                case 1
                    Screen('DrawTexture', window, stimTex(phases(3,trialNum)), [], [], params.refOrient(refOrder(trialNum)) + testOrientation(trialNum) + firstRefJitter);
                case -1
                    Screen('DrawTexture', window, stimTex(phases(3,trialNum)), [], [], params.refOrient(refOrder(trialNum)) - testOrientation(trialNum) + firstRefJitter);
            end
        case 2 % If cued second stimulus thenn add oreintation to second stimulus
            switch refOrder(trialNum) % Add it to second stimulus
                case 1
                    switch testDirection(trialNum) % Find out if going clockwise or anticlockwise
                        case 1
                            Screen('DrawTexture', window, stimTex(phases(2,trialNum)), [], [], params.refOrient(2) + testOrientation(trialNum) + secondRefJitter);
                        case -1
                            Screen('DrawTexture', window, stimTex(phases(2,trialNum)), [], [], params.refOrient(2) - testOrientation(trialNum) + secondRefJitter);
                    end
                case 2
                    switch testDirection(trialNum) % Find out if going clockwise or anticlockwise
                        case 1
                            Screen('DrawTexture', window, stimTex(phases(2,trialNum)), [], [], params.refOrient(1) + testOrientation(trialNum) + secondRefJitter);
                        case -1
                            Screen('DrawTexture', window, stimTex(phases(2,trialNum)), [], [], params.refOrient(1) - testOrientation(trialNum) + secondRefJitter);
                    end
            end
    end
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Fixation after test presentation time
    expectedTime = expectedTime + params.testPres;
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    if bitsiFlag == 1 % If using button box
        bitsiResponse.clearResponses(); % Clear bitsi lines so we read the right button press
        currState = bitsiResponse.serobj.BytesAvailable;
        stateChanged = currState;
        currTime = GetSecs;
        
        while stateChanged == currState;
            stateChanged = bitsiResponse.serobj.BytesAvailable;
            % Break out of loop if out of time
            if GetSecs > currTime + (params.response - .1)
                break
            end
        end
        % calculate reaction time
        rt(trialNum) = GetSecs - currTime;
        % Read response
        if stateChanged > currState % If button was pressed
            buttonPr(trialNum) = fread(bitsiResponse.serobj,1);
            switch buttonPr(trialNum)
                case bitsi_key1
                    response(trialNum) = 1; % Anti clockwise
                case bitsi_key2
                    response(trialNum) = 2; % Clockwise
            end
        end
        % Determine if response was correct
        % If orientatation was more clockiwise (larger) and pressed button 2, correct
        switch response(trialNum)
            case 2
                switch testDirection(trialNum)
                    case 1
                        correct(trialNum) = 1;
                    case -1
                        correct(trialNum) = 0;
                end
            case 1
                switch testDirection(trialNum)
                    case -1
                        correct(trialNum) = 1;
                    case 1
                        correct(trialNum) = 0;
                end
        end
        
        % Throw to workspace
        assignin('base', 'correct', correct);
        
    else % Otherwise use keyboard
        
        % Wait for response - wait for a maximum of ITI - 100ms for response,
        % the -100ms is a buffer to ensure we reach next trial on time
        [response{trialNum}, rt(trialNum)] = read_response_and_rt_window(params.response-.1);
        
        % Determine if response was correct
        % If orientatation was more clockiwise (larger) and pressed x, correct
        switch response{trialNum}
            case 'm'
                correct(trialNum) = 0; % Mark misses as incorrect
            case 'x'
                switch testDirection(trialNum)
                    case 1
                        correct(trialNum) = 1;
                    case -1
                        correct(trialNum) = 0;
                end
                % If orientation was anti clockwise (smaller) and pressed z,
                % correct
            case 'z'
                switch testDirection(trialNum)
                    case -1
                        correct(trialNum) = 1;
                    case 1
                        correct(trialNum) = 0;
                end
            case 'q'
                save_and_quit;
                return
        end
    end
    
    % Safe quit mechanism (hold q to quit)
    [keyPr,~,key,~] = KbCheck;
    if keyPr == 1 && KbName(key) == 'q'
        save_and_quit;
        return
    end
    
    % Update quest
    params.q = QuestUpdate(params.q,testOrientation(trialNum),correct(trialNum));
    
    % Give feedback
    expectedTime = expectedTime + params.response;
    switch correct(trialNum)
        case 1
            text = 'Correct';
        case 0
            text = 'Wrong';
    end
    DrawFormattedText(window, text, 'center', 'center');
    Screen('DrawingFinished', window);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % Flip back after feedback
    expectedTime = expectedTime + params.feedback;
    draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
    Screen('Flip', window, startTime + expectedTime - halfifi);
    
    % END OF TRIAL LOOP
    
    % Take a break every 10 minutes (31 trials)
    if mod(trialNum, 31) == 0
        takeBreak;
        tookBreak = 1; % To track that we just took a break
    end
end

Priority(0); % Reset priority

% Display end message
Screen('DrawText',window,'End of experiment, thank you for participating!', 100, 100);
Screen('DrawingFinished',window);
Screen('Flip',window);

pause(2);

% -----------------------------------------------------------------
% call function to save results, close window and clean up
save_and_quit;

    function save_and_quit
        
        % Caculate threshold
        threshold = QuestMean(params.q);
        thresh_sd = QuestSd(params.q);
        
        % First dump everything to workspace just in case something goes
        % wrong
        assignin('base', 'correct', correct);
        assignin('base', 'threshold', threshold);
        assignin('base', 'thresh_sd', thresh_sd);
        
        % Save out results
        try
            save(outFile, 'params', 'refOrder', 'phases', 'cue', 'testDirection', 'response', 'correct', 'testOrientation', 'rt', 'threshold', 'thresh_sd');
        catch
            disp('Could not save out_file - may not exist');
        end
        
        % Close window
        Screen('Close');
        Screen('Closeall');
        
        % Close serial port
        if bitsiFlag == 1
            bitsiResponse.close();
            disp('Serial port closed');
        end
        
        % Clean up 
        try
            jheapcl; % clean the java heap space
        catch
            disp('Could not clean java heap space');
        end
        
        disp('Quit safely');
    end

% Function to give the subject a short break ;)
    function takeBreak
        DrawFormattedText(window,'Feel free to take a short break\n Click the mouse to continue','center', 'center');
        Screen('Flip',window);
        
        clicks = 0;
        % Wait for click
        while clicks < 1
            clicks = GetClicks(window);
        end
                        
        % Draw fixation and wait for a couple of seconds
        draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
        Screen('DrawingFinished', window);
        % Reset timings since we may have waited a while
        startTime = Screen('Flip', window);
        expectedTime = 0;
          
        expectedTime = expectedTime + 2;
        draw_fixationdot(window,background, x, y,[0.24 0.06],0,0,distFromScreen);
        Screen('DrawingFinished', window);
        Screen('Flip', window, startTime + expectedTime - halfifi);
    end
end