% data = data(1:12);
% params.Blocks = params.Blocks(1:12);
contrasts_nonLog = params.contrast_probe;
contrasts_nonLog = 0:1/(length(contrasts_nonLog)-1):1; % converting into a 0-1 linear space. 
contrasts_nonLogAnswers = contrasts_nonLog;
% assuming that the log difference between contrast(min+1)-contrast(min) ==
% contrast(max) - contrast(max-1)
% while the total difference between contrast(max)-contrast(min) = 1; 
conrange=(params.stimrange-1)/2;
PalContrasts_probe = round(params.contrast_probe/params.max_con*conrange);
palcon = {'13','19','27','40',' 59','87','127'};

responseDiff = {};
correctAnswer = {};
answerChosen = {};

for i = 1:length(data)
    for ii = 1:length(contrasts_nonLog)
        data{i}.answer(data{i}.answer == PalContrasts_probe(ii)) = contrasts_nonLog(ii); 
        data{i}.correctAnswer(data{i}.correctAnswer == PalContrasts_probe(ii)) = contrasts_nonLog(ii);
    end 
    
    responseDiff{i} = data{i}.correctAnswer - data{i}.answer;  
    correctAnswer{i} = data{i}.correctAnswer; 
    answerChosen{i} = data{i}.answer; 
end

One_stim = {};
Two_stim = {};

One_stim.total = responseDiff(params.Blocks==2);
One_stim.correctAnswer = correctAnswer(params.Blocks==2);
One_stim.answer = answerChosen(params.Blocks==2);
One_stim.meanBlocks = [];
One_stim.seBlock = [];
% One_stim.stdBlock = [];
One_stim.viableTrials = [];

Two_stim.total = responseDiff(params.Blocks==1);
Two_stim.correctAnswer = correctAnswer(params.Blocks==1);
Two_stim.answer = answerChosen(params.Blocks==1);
Two_stim.cues = data(params.Blocks==1); 
Two_stim.meanBlocks = [];
Two_stim.seBlock = [];
Two_stim.stdBlock = [];
Two_stim.viableTrials = [];
Two_stim.viableCues = [];
Two_stim.viableTrialsBlock = {};
Two_stim.viableCuesBlock = {};

One_stim.viableAnswersCorrect = [];
One_stim.viableAnswersChosen = [];

Two_stim.viableAnswersCorrect = [];
Two_stim.viableAnswersChosen = [];

for mb = 1:length(Two_stim.total)
    One_stim.viableTrials = [One_stim.viableTrials, abs(One_stim.total{mb}(abs(One_stim.total{mb})<1.1))];
    One_stim.viableAnswersCorrect = [One_stim.viableAnswersCorrect, abs(One_stim.correctAnswer{mb}(abs(One_stim.total{mb})<1.1))];
    One_stim.viableAnswersChosen = [One_stim.viableAnswersChosen, abs(One_stim.answer{mb}(abs(One_stim.total{mb})<1.1))];

    One_stim.meanBlocks=[One_stim.meanBlocks,mean(abs(One_stim.total{mb}(abs(One_stim.total{mb})<1.1)))];
    One_stim.seBlock = [One_stim.seBlock, std(abs(One_stim.total{mb}(abs(One_stim.total{mb})<1.1)))/sqrt(length((One_stim.total{mb}(abs(One_stim.total{mb})<1.1))))];
%     One_stim.stdBlock = [One_stim.stdBlock, std(abs(One_stim.total{mb}(abs(One_stim.total{mb})<1.1)))];
    
    Two_stim.viableTrials = [Two_stim.viableTrials, abs(Two_stim.total{mb}(abs(Two_stim.total{mb})<1.1))];
    Two_stim.viableAnswersCorrect = [Two_stim.viableAnswersCorrect, abs(Two_stim.correctAnswer{mb}(abs(Two_stim.total{mb})<1.1))];
    Two_stim.viableAnswersChosen = [Two_stim.viableAnswersChosen, abs(Two_stim.answer{mb}(abs(Two_stim.total{mb})<1.1))];

    Two_stim.viableCues = [Two_stim.viableCues, Two_stim.cues{mb}.cue(abs(Two_stim.total{mb})<1.1)];
    
    Two_stim.viableTrialsBlock{mb} = abs(Two_stim.total{mb}(abs(Two_stim.total{mb})<1.1));
    Two_stim.viableCuesBlock{mb} = Two_stim.cues{mb}.cue(abs(Two_stim.total{mb})<1.1); 
end    
 
answersProbes = [[One_stim.viableAnswersCorrect,Two_stim.viableAnswersCorrect];[One_stim.viableAnswersChosen,Two_stim.viableAnswersChosen]];
% answersProbes = [One_stim.viableAnswersCorrect;One_stim.viableAnswersChosen];
% answersProbes = [Two_stim.viableAnswersCorrect;Two_stim.viableAnswersChosen];

matrixOnes = zeros(length(contrasts_nonLog),length(contrasts_nonLogAnswers));

for conTrue = 1:length(contrasts_nonLogAnswers)
   for conAnswer = 1:length(contrasts_nonLog)
       for trial = 1:length(answersProbes)
           if answersProbes(1,trial)==contrasts_nonLogAnswers(conTrue) && answersProbes(2,trial)==contrasts_nonLog(conAnswer)
               matrixOnes(conAnswer,conTrue) = matrixOnes(conAnswer,conTrue)+1; 
           end
       end
   end
end

fig4 = figure(4)
imagesc(matrixOnes)
colorbar
set(gca,'xticklabel',round(contrasts_nonLogAnswers,2))
set(gca, 'YTick', 1:13,'yticklabel',round(contrasts_nonLog,2))
set(gca, 'XTick', 1:13,'xticklabel',round(contrasts_nonLog,2))

xlabel('Conditions') 
ylabel('Answers') 
title('Across Tasks')



% Mean across blocks

One_stim_mean = mean(One_stim.viableTrials);
One_stim_se = std(One_stim.viableTrials)/sqrt(length(One_stim.viableTrials));

Two_stim_Cue_1_mean = mean(Two_stim.viableTrials(Two_stim.viableCues==1)); %n=60
Two_stim_Cue_2_mean = mean(Two_stim.viableTrials(Two_stim.viableCues==2)); %n=59

Two_stim_Cue_1_se = std(Two_stim.viableTrials(Two_stim.viableCues==1))/sqrt(length(Two_stim.viableTrials(Two_stim.viableCues==1))); 
Two_stim_Cue_2_se = std(Two_stim.viableTrials(Two_stim.viableCues==2))/sqrt(length(Two_stim.viableTrials(Two_stim.viableCues==2)));


X = {'One stimulus','Retrocue: cue 1','Retrocue: cue 2'};
% X = reordercats(X,{'One stimulus','Retrocue: cue 1','Retrocue: cue 2'});
Y = [One_stim_mean,Two_stim_Cue_1_mean,Two_stim_Cue_2_mean];
fig1=figure(1)
b1 = bar(Y,'grouped')
set(gca,'xticklabel',X)
hold on
er = errorbar([One_stim_mean,Two_stim_Cue_1_mean,Two_stim_Cue_2_mean],[One_stim_se,Two_stim_Cue_1_se,Two_stim_Cue_2_se]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xlabel('Conditions') 
ylabel('Distance from correct answer (range from 0-1)') 

%% Take a look at differences for each luminance level when probed
One_stim.viableTrialsIndex = abs([One_stim.total{:}])<1.1;
Two_stim.viableTrialsIndex = abs([Two_stim.total{:}])<1.1;

one_stim_correct = [correctAnswer{params.Blocks==2}];
one_stim_correct = one_stim_correct(One_stim.viableTrialsIndex); 

two_stim_correct = [correctAnswer{params.Blocks==1}];
two_stim_correct = two_stim_correct(Two_stim.viableTrialsIndex); 

meanContrastOne = [];
meanContrastTwo = [];
seContrastOne = [];
seContrastTwo = [];
numContrastOne = [];
numContrastTwo = [];

meanContrastBoth = [];
seContrastBoth = [];
numContrastsBoth = [];

for lum = 1:length(contrasts_nonLogAnswers)
    ix_correct_one = one_stim_correct==contrasts_nonLogAnswers(lum);
    ix_correct_two = two_stim_correct==contrasts_nonLogAnswers(lum);
    meanContrastOne = [meanContrastOne,mean(One_stim.viableTrials(ix_correct_one))];
    meanContrastTwo = [meanContrastTwo,mean(Two_stim.viableTrials(ix_correct_two))];
    seContrastOne = [seContrastOne, std(One_stim.viableTrials(ix_correct_one))/sqrt(length(One_stim.viableTrials(ix_correct_one)))];
    seContrastTwo = [seContrastTwo, std(Two_stim.viableTrials(ix_correct_two))/sqrt(length(Two_stim.viableTrials(ix_correct_two)))];
    numContrastOne = [numContrastOne,length(One_stim.viableTrials(ix_correct_one))];
    numContrastTwo = [numContrastTwo,length(Two_stim.viableTrials(ix_correct_two))];
    
    meanContrastBoth = [meanContrastBoth,mean([One_stim.viableTrials(ix_correct_one),Two_stim.viableTrials(ix_correct_two)])];
    seContrastBoth = [seContrastBoth,std([One_stim.viableTrials(ix_correct_one),Two_stim.viableTrials(ix_correct_two)])/sqrt(length([One_stim.viableTrials(ix_correct_one),Two_stim.viableTrials(ix_correct_two)]))];
    numContrastsBoth = [numContrastsBoth,length([One_stim.viableTrials(ix_correct_one),Two_stim.viableTrials(ix_correct_two)])];
end
Co = strsplit(num2str(round(contrasts_nonLogAnswers,2)));
Co2 = palcon;
fig3 = figure(3)
subplot(3,2,1)
bar(meanContrastOne)
set(gca,'xticklabel',Co)
hold on
er = errorbar(meanContrastOne,seContrastOne);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xlabel('Contrast') 
ylabel('Distance/error') 
ylim([0 0.6])
title('One stimulus task contrast error')

subplot(3,2,2)
bar(numContrastOne)
set(gca,'xticklabel',Co)
xlabel('Contrast') 
ylabel('Number of probes') 
ylim([0 10])
title('One stimulus task: number of probes')

subplot(3,2,3)
bar(meanContrastTwo)
set(gca,'xticklabel',Co)
hold on
er = errorbar(meanContrastTwo,seContrastTwo);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xlabel('Contrast') 
ylabel('Distance/error') 
ylim([0 0.6])
title('Retrocue task contrast error')

subplot(3,2,4)
bar(numContrastTwo)
set(gca,'xticklabel',Co)
xlabel('Contrast') 
ylabel('Number of probes') 
ylim([0 10])
title('Retrocue task: number of probes')

subplot(3,2,5)
bar(meanContrastBoth)
set(gca,'xticklabel',Co2)
hold on
er = errorbar(meanContrastBoth,seContrastBoth);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
xlabel('Contrast') 
ylabel('Distance/error')
ylim([0 0.6])
title('Across tasks contrast error')

subplot(3,2,6)
bar(numContrastsBoth)
set(gca,'xticklabel',Co2)
xlabel('Contrast') 
ylabel('Number of probes') 
title('Across tasks: number of probes')



%% Differences across blocks
Two_stim_cue1_meanPB = [];
Two_stim_cue1_meanPBse = [];

Two_stim_cue2_meanPB = [];
Two_stim_cue2_meanPBse = [];


for ccb = 1:length(Two_stim.viableTrialsBlock)
    Two_stim_cue1_meanPB = [Two_stim_cue1_meanPB,mean(Two_stim.viableTrialsBlock{ccb}(Two_stim.viableCuesBlock{ccb}==1))];
    Two_stim_cue1_meanPBse = [Two_stim_cue1_meanPBse,std(Two_stim.viableTrialsBlock{ccb}(Two_stim.viableCuesBlock{ccb}==1))/sqrt(length(Two_stim.viableTrialsBlock{ccb}(Two_stim.viableCuesBlock{ccb}==1)))];

    Two_stim_cue2_meanPB = [Two_stim_cue2_meanPB,mean(Two_stim.viableTrialsBlock{ccb}(Two_stim.viableCuesBlock{ccb}==2))];
    Two_stim_cue2_meanPBse = [Two_stim_cue2_meanPBse,std(Two_stim.viableTrialsBlock{ccb}(Two_stim.viableCuesBlock{ccb}==2))/sqrt(length(Two_stim.viableTrialsBlock{ccb}(Two_stim.viableCuesBlock{ccb}==2)))];
end

y = [One_stim.meanBlocks;Two_stim_cue1_meanPB;Two_stim_cue2_meanPB]';
erb = [One_stim.seBlock,Two_stim_cue1_meanPBse,Two_stim_cue2_meanPBse];

fig2=figure(2)
bar(y)
legend(X)
xlabel('Blocks') 
ylabel('Distance from correct answer (range from 0-1)') 

% error('Dont overwrite');
saveas(fig1,'KW_MeanBlocks.png')
saveas(fig2,'KW_Blocks.png')
saveas(fig3,'KW_acrossContrast.png')
saveas(fig4,'KW_histogram.png')

