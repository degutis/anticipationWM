function stim_circcheck_pal2(subid,triggermode)

sinus_pattern=1;
scale=4;
port_num=1; % needs to be set to 1 at the scanner
screen_num=2;
dummy_Vols=3;
num_images=4;
num_runs=1;
num_pres=num_runs*num_images;
luminance=0.8;
max_con=0.1;
min_con=0.01;
contrasts=logspace(log10(min_con),log10(max_con),4);
mean_con=logspace(log10(contrasts(2)),log10(contrasts(3)),3);
mean_con=mean_con(2);
% define the spatial geometry of the stimulus
n_circ=4; % number of circles
n_ang=12; % number of sectors per circle
r_scaling=2; % length scaling of sectors.

stimrange=221;
PalRGB=repmat(linspace(luminance/2-luminance*max_con/2,luminance/2+luminance*max_con/2,stimrange)',1,3);
bgcol=ceil(stimrange/2);
white=stimrange;
white_landolt=stimrange+1;
black=1;
PalRGB(stimrange+1,:)=luminance;
conrange=(stimrange-1)/2;
PalContrasts=round(contrasts/max_con*conrange);
mean_con=round(mean_con/max_con*conrange);

stim_order=randperm(num_pres);
stim_num=1:4;

curdir=pwd;

clock_int=round(clock);
date_str=sprintf('%02.0f%02.0f%02.0f',clock_int(1)-2000,clock_int(2),clock_int(3));
time_str=sprintf('%02.0f%02.0f',clock_int(4),clock_int(5));

% different timing parameters.
t_gray=3000; % 30000
TR=1500;
t_dummy=dummy_Vols*TR;
t_gray=t_dummy;

tot_image=5000; % total time for flashing one image
t_total=num_images*num_runs*2*tot_image+t_dummy+t_gray; %(Vols+dummy_Vols)*TR;
t_fix_change=800; % frequency of fixation changes 

tchange=zeros(10000,1); n=1; %vector to save timing of images.
resolution=1024*768; % number of pixels per image.

%statistik von fixation überprüfen!!!!
fixation=fixation_random([5 7],t_total/t_fix_change); %t/t_image/2+1);

%% save all the data ign the stuct res
% information about the person, experiment, date etc.
res.info.subj=subid;
res.info.paradigm='Checkerring, local contrast variation';
res.info.date=date_str;
res.info.time=time_str;
res.info.m_file=mfilename;
res.info.stim_file=sprintf('stim%.0f_%.0f',n_circ,n_ang);

%% load the stimulus 
eval(sprintf('load stim%.0f_%.0f',n_circ,n_ang));;
if sinus_pattern
    [ch_board]=random_checker_sin(n_circ,n_ang,r_scaling,scale)';
end

% information about the parameters of the experiment
res.para.fixations=fixation;
res.para.luminance=luminance;
res.para.timage=zeros(num_images,1);
res.para.t_fix_change=t_fix_change;
res.para.t_total=t_total;
res.para.con_values=contrasts;
res.para.con_valuesPAL=PalContrasts;
res.para.tgray=t_gray;
res.para.n_circ=n_circ;
res.para.n_ang=n_ang;
res.para.r_scaling=r_scaling;
res.para.scale=scale;
res.para.contrasts=t_4(:,50+(1:num_images));

% configure cogent for the sepcific run
% and 
if ~isstr(subid)
    subid=num2str(subid);
end
c=clock;
filename=['./logs/CircRun_' subid '_' sprintf('%02d',c(1)) '_' sprintf('%02d',c(2)) '_' int2str(c(3)) '_' sprintf('%02d',c(4)) '_' sprintf('%02d',c(5)) '_' sprintf('%02d',fix(c(6)))];
config_display(screen_num, 3, [0 0 0], [luminance luminance luminance], 'Arial', 24, 5 ,8)
%config_keyboard;
config_keyboard(1000,5,'nonexclusive')

config_serial(port_num,19200,0,0,8);
config_log( [filename '.log']);
save([filename '.mat'],'res');

% start the cogent application and read images into the memory.
start_cogent;

cgcoltab(1,PalRGB);
cgnewpal;

cgflip(bgcol);
cgflip(bgcol);

num_segments=48;
con_im=ones(1024,768);
hh=zeros(1024,768);
disp('%%% loading images ...');
for i=1:num_images%length(files)
    for ll=1:num_segments; 
        con_im(indic{ll})=i;
    end, 
    hh(:,:)=bgcol+PalContrasts(con_im).*ch_board;
    cgloadarray(20+(2*i)-1,1024,768,round(reshape(hh,resolution,1))',PalRGB,0);
    cgloadarray(20+2*i,1024,768,white+1-round(reshape(hh,resolution,1))',PalRGB,0);
end
    for ll=1:num_segments; 
        con_im(indic{ll})=mean_con; 
    end
hh(:,:)=bgcol+con_im.*ch_board(:,:);
cgloadarray(1,1024,768,round(reshape(hh,resolution,1))',PalRGB,0);
cgloadarray(2,1024,768,white+1-round(reshape(hh,resolution,1))',PalRGB,0);
   
clear hh con_im
clear picture

%%%%%%%%%%%%%%%%%%landolt c's zeichnen%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
landolt_c=[-5 5; 5 5; 5 -5; -5 -5];
draw_landold_c(landolt_c,bgcol,white_landolt,PalRGB);

cgsetsprite(0);
cgflip(bgcol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cgtext('READY...',0,0);
cgflip(bgcol);

if triggermode==1
    clearserialbytes(port_num);
    waitserialbyte(port_num,inf,53);
else
    waitkeydown(inf);
end;

show_fix=0;
clearserialbytes(port_num);
t0=time;
disp('Experiment running ...');
res.para.tstart=t0;
start_next=0;
show_fix=0;
z=0;
f=0;
img_counter=0;
flicker=1;
while time<(t0+t_dummy)
    tflip=time;
    z_old=z;
    z=ceil((time-t0)/t_fix_change);
        flicker=not(flicker);
            if flicker
        cgdrawsprite(1,0,0);
    else
        cgdrawsprite(2,0,0);
            end
    if z_old~=z
        if show_fix
            show_fix=0;
        else
            f=f+1;
            show_fix=1;
        end
    end
    if show_fix
        cgdrawsprite(fixation(f),0,0);
    else
        cgdrawsprite(9,0,0);
    end
    
    waituntil(tflip+100);
    cgflip(bgcol);
end
t_start_images=time;
while time<t0+t_total-t_gray
    img_counter=img_counter+1;
    disp(sprintf('Showing image number : %03.0f',img_counter)); 
    start_next=t_start_images+img_counter*tot_image*2;
   
    res.para.timage(img_counter)=time-t0;
    while time<(start_next-tot_image)
        tflip=time;
        flicker=not(flicker);
            if flicker
        cgdrawsprite(20+stim_num(img_counter)*2-1,0,0);
    else
        cgdrawsprite(20+stim_num(img_counter)*2,0,0);
            end

        z_old=z;
        z=ceil((time-t0)/t_fix_change);
        if z_old~=z
            if show_fix
                show_fix=0;
            else
                f=f+1;
                show_fix=1;
            end
        end
        if show_fix
            cgdrawsprite(fixation(f),0,0);
        else
            cgdrawsprite(9,0,0);
        end

        tchange(n)=time;
        n=n+1;
        waituntil(tflip+100);
        cgflip(bgcol);
    end
    while time<start_next
        tflip=time;

        z_old=z;
        z=ceil((time-t0)/t_fix_change);
          flicker=not(flicker);
            if flicker
        cgdrawsprite(1,0,0);
    else
        cgdrawsprite(2,0,0);
            end
        if z_old~=z
            if show_fix
                show_fix=0;
            else
                f=f+1;
                show_fix=1;
            end
        end
        if show_fix
            cgdrawsprite(fixation(f),0,0);
        else
            cgdrawsprite(9,0,0);
        end

        tchange(n)=time;
        n=n+1;
        waituntil(tflip+100);
        cgflip(bgcol);
    end
    res.timeimg=time;

end;
t_end_images=time;
while time<(t0+t_total)
    tflip=time;
    z_old=z;
    z=ceil((time-t0)/t_fix_change);
      flicker=not(flicker);
            if flicker
        cgdrawsprite(1,0,0);
    else
        cgdrawsprite(2,0,0);
            end
    if z_old~=z
        if show_fix
            show_fix=0;
        else
            f=f+1;
            show_fix=1;
        end
    end
    if show_fix
        cgdrawsprite(fixation(f),0,0);
    else
        cgdrawsprite(9,0,0);
    end
        
    
    waituntil(tflip+100);
    cgflip(bgcol);
end

t_end=time;
res.para.tend=t_end;
readserialbytes(port_num);
logserialbytes(port_num);
[k_value,k_time,n]=getserialbytes(port_num);
res.results.serial_values=k_value;
res.results.serial_times=k_time-t0;;
res.para.t_end_images=t_end_images;
res.para.t_start_images=t_start_images;

save([filename '.mat'],'res');

stop_cogent;
end


%% END OF ACTUAL FUNCTION %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Funktionsteil%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filename=my_config_cogent(subid)


end


function draw_landold_c(landolt_c,bgcol,white,PalRGB)

cgcoltab(1,PalRGB);
cgnewpal;
sprite_size=20;
cgmakesprite(5,sprite_size,sprite_size,bgcol);
cgsetsprite(5);
cgpencol(bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white);
i=[1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[4 3 2 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(5,'n');

cgmakesprite(6,sprite_size,sprite_size,bgcol);
cgsetsprite(6);
cgpencol(bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white);
i=[2 3 4 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[1 4 3 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(6,'n');

cgmakesprite(7,sprite_size,sprite_size,bgcol);
cgsetsprite(7);
cgpencol(bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white);
i=[3 4 1 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(7,'n');

cgmakesprite(8,sprite_size,sprite_size,bgcol);
cgsetsprite(8);
cgpencol(bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white);
i=[4 1 2 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(8,'n');

cgmakesprite(9,sprite_size,sprite_size,bgcol);
cgsetsprite(9);
cgpencol(bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white);
i=[4 1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(9,'n');

end

function r=fixation_random(A,n)
B=randperm(length(A));
r(1)=A(B(1));
r(2)=A(B(2));
B=randperm(length(A));
r(3)=A(B(1));
i=3;
while i<n
    B=randperm(length(A));
    i=i+1;
    r(i)=A(B(1));
end
end