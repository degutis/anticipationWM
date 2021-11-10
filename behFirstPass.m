% Analysis of behavioral data

open 'data/001_KM/001_KM_1_WMdistAnt_V1.mat'

data.error = data.Answer - data.StimulusOr;
data.error = (data.error+90)-floor((data.error+90)/180)*180 - 90;


%% Between different trial types

fig = figure(1);
subplot(1,2,1)
violinplot(struct('ZERO',abs(data.error(data.Trial_type==0)),'FIFTY',abs(data.error(data.Trial_type==50)),'HUNDERED',abs(data.error(data.Trial_type==100))))
ylabel('Absolute angular error (degrees)')
subplot(1,2,2)
violinplot(struct('ZERO',abs(data.RT(data.Trial_type==0)),'FIFTY',abs(data.RT(data.Trial_type==50)),'HUNDERED',abs(data.RT(data.Trial_type==100))))
ylabel('Reaction time (seconds)')

unique0 = unique(data.error(data.Trial_type==0));
unique0 = unique0(~isnan(unique0));
[a,b]=hist(data.error(data.Trial_type==0),unique0)



[p,tbl,stats] = anova1([abs(data.error(data.Trial_type==0));abs(data.error(data.Trial_type==50));abs(data.error(data.Trial_type==100))],[repmat({'zero'},length(abs(data.error(data.Trial_type==0))),1);repmat({'fifty'},length(abs(data.error(data.Trial_type==50))),1);repmat({'hundered'},length(abs(data.error(data.Trial_type==100))),1)])


%% Distractor vs no distractor

fig8 = figure(8);
subplot(1,2,1)
violinplot(struct('noDist',abs(data.error(data.DistractorOr==0)),'Dist',abs(data.error(data.DistractorOr~=0))))
[h2,p2,ci2,stats2]=ttest2(abs(data.error(data.DistractorOr==0)),abs(data.error(data.DistractorOr~=0)));
subplot(1,2,2)
violinplot(struct('noDist',abs(data.RT(data.DistractorOr==0)),'Dist',abs(data.RT(data.DistractorOr~=0))))


fig9 = figure(9);
subplot(1,2,1)
violinplot(struct('noDist50',abs(data.error(data.DistractorOr==0 & data.Trial_type==50)),'noDist0',abs(data.error(data.Trial_type==0))))
[~,pval,~]=ttest2(abs(data.error(data.DistractorOr==0 & data.Trial_type==50)),abs(data.error(data.Trial_type==0)));
title(['No Distractor trials: 50% vs 0% blocks. Pval: ',num2str(pval)])
ylim([0 60])
subplot(1,2,2)
violinplot(struct('Dist50',abs(data.error(data.DistractorOr~=0 & data.Trial_type==50)),'Dist100',abs(data.error(data.Trial_type==100))))
ylim([0 60])
[~,pval,~]=ttest2(abs(data.error(data.DistractorOr~=0 & data.Trial_type==50)),abs(data.error(data.Trial_type==100)));
title(['Distractor trials: 50% vs 100% blocks. Pval: ',num2str(pval)])



%% Distractor orientations
fig10 = figure(10)
scatter(data.DistractorOr(data.DistractorOr~=0),data.error(data.DistractorOr~=0))
hold on
plot(-80:80, repmat(0,length(-80:80),1))
ylim([-60 60])
xlabel('Distractor Orientation')
ylabel('Error')

fig11 = figure(11)
scatter(abs(data.DistractorOr(data.DistractorOr~=0)),data.error(data.DistractorOr~=0))


data.normalizedBias = (data.error./data.DistractorOr)*100;

fig12 = figure(12)
scatter(data.DistractorOr(data.DistractorOr~=0),data.normalizedBias(data.DistractorOr~=0))
hold on
plot(-80:80, repmat(0,length(-80:80),1))
%ylim([-60 60])
xlabel('Distractor Orientation')
ylabel('Error')
