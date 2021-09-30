scale=4;
luminance=0.8;
max_con=1;
min_con=0.1;
contrasts=logspace(log10(min_con),log10(max_con),4);
mean_con=logspace(log10(contrasts(2)),log10(contrasts(3)),3);
mean_con=mean_con(2);
% define the spatial geometry of the stimulus
n_circ=4; % number of circles/2
n_ang=12; % number of sectors per circle/2
r_scaling=2; % length scaling of sectors.
circle_param = 4; %how many circles are included in one item  
angle_param = 4; % how many thetas are included in one item
stimrange=221;%221 orig
rect = [0 0 1920 1080];

PalRGB=repmat(linspace(luminance/2-luminance*max_con/2,luminance/2+luminance*max_con/2,stimrange)',1,3);
bgcol=ceil(stimrange/2);
PalRGB(stimrange+1,:)=luminance;
conrange=(stimrange-1)/2;
PalContrasts=round(contrasts/max_con*conrange);
% mean_con=round(mean_con/max_con*conrange);

[ch_board, indic]=random_checker_sin(n_circ,n_ang,r_scaling,scale, circle_param, angle_param, rect);
% ch_board = ch_board'; 

con_im=ones(rect(3),rect(4));
hh=zeros(rect(3),rect(4));

num_segments = length(indic);
repetition = 1;
contrasts2 = [0,0.33333,0.66667,1];
contrasts_full = repmat(contrasts2,1,num_segments/(repetition*length(contrasts2)));
contrasts_full = contrasts_full(randperm(length(contrasts_full)));
contrasts_full = repelem(contrasts_full,repetition); 

t_4 = contrasts_full';

for ll=1:num_segments
    indic_now = indic{ll};
    for rc = 1:length(indic_now)
        con_im(indic_now(rc,2),indic_now(rc,1))=round(t_4(ll)*3)+1;
    end 
end
con_im = con_im';
PalContrasts = PalContrasts';

hh=bgcol+PalContrasts(con_im).*ch_board;

figure(2)
imagesc(hh)