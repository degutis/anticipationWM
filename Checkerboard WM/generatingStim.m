%% Generate stimuli 

%% Some pre-specifications 
params.contrasts = 1:7;
params.numBlocksDiff = 2;
params.numBlocksPer = 5;
params.numTrials = 20;
params.Blocks = repelem(1:params.numBlocksDiff,params.numBlocksPer);
params.Blocks = params.Blocks(randperm(length(params.Blocks))); 
indic = 1:4; 

%% Code
