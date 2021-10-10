% Analysis of behavioral data

open 'data/001_KM/001_KM_1_WMdistAnt_V1.mat'

data.Answer = deg2rad(mod(data.Answer,180));
data.StimulusOr = deg2rad(mod(data.StimulusOr,180)); 

data.error = data.Answer - data.StimulusOr;
data.error = (a+90)-floor((a+90)/180)*180 - 90;


%% Between different trial types

fig = figure(1);
violinplot([abs(data.error(data.Trial_type==0)),abs(data.error(data.Trial_type==50)),abs(data.error(data.Trial_type==100))],{'0%','50%','100%'})


fig2 = figure(2);

count(data.error(data.Trial_type==0))


unique0 = unique(data.error(data.Trial_type==0));
unique0 = unique0(~isnan(unique0));


[a,b]=hist(data.error(data.Trial_type==0),unique0)

bar(b,a)

f=fit(b,a','gauss')
plot(f,b',a)


unique(data.error(data.Trial_type==0)))


f = fit(x.',y.','gauss2')


[p,tbl,stats] = anova1([abs(data.error(data.Trial_type==0)),abs(data.error(data.Trial_type==50)),abs(data.error(data.Trial_type==100))])


%% Distractor vs no distractor

fig2 = figure(4);
violinplot([abs(data.error(data.DistractorOr==0)),abs(data.error(data.DistractorOr~=0))],{'No Distractor','Distractor'})

[h2,p2,ci2,stats2]=ttest2(abs(data.error(data.DistractorOr==0)),abs(data.error(data.DistractorOr~=0)));


%% Distractor orientations
fig3 = figure(8)
violinplot([data.error(data.DistractorOr==-30),data.error(data.DistractorOr==-15),data.error(data.DistractorOr==-5),data.error(data.DistractorOr==5),data.error(data.DistractorOr==15),data.error(data.DistractorOr==30)],{'-30','-15','-5','5','15','30'})

[p3,tbl3,stats3] = anova1([data.error(data.DistractorOr==-30),data.error(data.DistractorOr==-15),data.error(data.DistractorOr==-5),data.error(data.DistractorOr==5),data.error(data.DistractorOr==15),data.error(data.DistractorOr==30)])


data.normalizedBias = (data.error./data.DistractorOr)*100;

fig4 = figure(12)
violinplot([data.normalizedBias(data.DistractorOr==-30 ),data.normalizedBias(data.DistractorOr==-15),data.normalizedBias(data.DistractorOr==-5),data.normalizedBias(data.DistractorOr==5),data.normalizedBias(data.DistractorOr==15),data.normalizedBias(data.DistractorOr==30)],{'-30','-15','-5','5','15','30'})


fig4 = figure(16)
violinplot([data.normalizedBias(data.DistractorOr==-30 | data.DistractorOr==30),data.normalizedBias(data.DistractorOr==-15 | data.DistractorOr==15),data.normalizedBias(data.DistractorOr==-5 | data.DistractorOr==5)],{'abs(30)','abs(15)','abs(5)'})

ttest(data.normalizedBias(data.DistractorOr==-30 | data.DistractorOr==30))
ttest(data.normalizedBias(data.DistractorOr==-15 | data.DistractorOr==15))
ttest(data.normalizedBias(data.DistractorOr==-5 | data.DistractorOr==5))


