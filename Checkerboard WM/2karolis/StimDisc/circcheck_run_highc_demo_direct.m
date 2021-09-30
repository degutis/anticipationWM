% function circcheck_run_highc_demo_direct(subid, run_ind, triggermode)
%
% To start the ciruclar checkerboard stimulation with Landolt C task, execute
% the following two steps:
% 
% Add COGENT (\Toolbox) to path 
% http://www.vislab.ucl.ac.uk/cogent.php  
% e.g. 
%   addpath('C:\Cogent2000v1.32\Toolbox')
% 
% Start a demo by calling
%     circcheck_run_highc_demo_direct('SBJNAME',1,0)
% where 
%   'SBJNAME': is the ID of the subject (can be string or number)
%   1: That the first run should be shown (other numbers would show also
%      show one run but with different stimuli).
%   0: use keyboard input (1 would mean serial port input; if the serial
%      port should be used, set the portnum variable in this script).
% This will show one run which takes 5:30 mins and shows 100 images.
%
% If multiple runs should be shown, call 
%   circcheck_run_highc_demo_direct('SBJNAME',1,0)
%   circcheck_run_highc_demo_direct('SBJNAME',2,0)
%     etc. currently up to 
%   circcheck_run_highc_demo_direct('SBJNAME',10,0)
% 
% Startup might take a bit. Log files will be created in \logs. Press any key
% to start after loading.
%
% If you like to check whether images are shown correctly already in the
% beginning, set "show_images_during_startup = 1" below.
% 
% The script runs correctly if you a fast flickering ciruclar checkerboard 
% (10Hz, currently hardcoded) with a Landolt C task in the middle. Flickering 
% exchanegs white to black and vice versa. 
% 
% In the beginning, the checkboard has low contrast. Then 100 flickering
% checkerboards with different highlighted sectors (medium and high contrast)
% are shown, and finally the uniform low contrast checkboard again. 
%
% If the script does not flicker:
%   - reduce the number of images (e.g. set num_images=5 below) 
%   - or use a computer with a better graphics card (more GPU memory)
% 
% The script has been tested with versions Cogent2000v1.32 and 2000v1.29 on 
% Windows XP, Matlab2010b and Matlab 2015a and 2011b with Cogent2000v1.32
% Windows 8.
% 
% On old computers, that are too slow to show the stimulation and still 
% support a 256 color mode, you can try 
%   circcheck_run_highc_demo_pal('SBJNAME',1,0)
% instead. NOTE: This demo lacks much of the logging capacities of the direct
% demo. If the palette mode stimulation should be used, the logging needs
% to be adapted.
% 
% Courtesy to Jakob Heinzle and John-Dylan Haynes.
% Original script by Jakob Heinzle; tested and adapted by Kai, 2015/08/03

function circcheck_run_highc_demo_direct(subid, run_ind, triggermode)

res = []; % init result struct

sinus_pattern=1;
scale=4;
useSerial = 0; % 1 if serial input is to be used
res.port_num=1; % only important if serial is used
screen_num=0;
dummy_Vols=10;
num_images=100;  %CHANGE
luminance=0.8;
max_con=1;
min_con=0.1;
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

curdir=pwd;

clock_int=round(clock);
date_str=sprintf('%02.0f%02.0f%02.0f',clock_int(1)-2000,clock_int(2),clock_int(3));
time_str=sprintf('%02.0f%02.0f',clock_int(4),clock_int(5));

% different timing parameters.
t_gray=3000; % 30000
TR=1500;
t_dummy=dummy_Vols*TR;
t_gray=t_dummy;

tot_image=3000; % total time for flashing one image
t_total=num_images*tot_image+t_dummy+t_gray; %(Vols+dummy_Vols)*TR;
t_fix_change=800; % frequency of fixation changes

tchange=zeros(10000,1); n=1; %vector to save timing of images.
resolution=1024*768; % number of pixels per image.

%statistik von fixation überprüfen!!!!
fixation=fixation_random([5 7],t_total/t_fix_change); %t/t_image/2+1);

% specify whether images should be loaded already during startup
show_images_during_startup = 1; % 1: show, 0: just load (blank screen)

%% save all the data ign the stuct res
% information about the person, experiment, date etc.
res.info.subj=subid;
res.info.paradigm='Checkerring, local contrast variation';
res.info.date=date_str;
res.info.time=time_str;
res.info.m_file=mfilename;
res.info.stim_file=sprintf('stim%.0f_%.0f',n_circ,n_ang);

%% load the stimulus
eval(sprintf('load stim%.0f_%.0f',n_circ,n_ang));
load randomvect;
rvec=randomvect{run_ind};
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
res.para.randomvect=rvec;

% configure cogent for the sepcific run
% and
if ~isstr(subid)
    subid=num2str(subid);
end
c=clock;
filename=['./logs/CircRun_' subid '_' sprintf('%02d',c(1)) '_' sprintf('%02d',c(2)) '_' int2str(c(3)) '_' sprintf('%02d',c(4)) '_' sprintf('%02d',c(5)) '_' sprintf('%02d',fix(c(6)))];
config_display(screen_num, 3, [0 0 0], [luminance luminance luminance], 'Arial', 24, 5 ,0)
%config_keyboard;
config_keyboard(1000,5,'nonexclusive')

if useSerial
    config_serial(res.port_num,19200,0,0,8);
end
config_log( [filename '.log']);
save([filename '.mat'],'res');

% start the cogent application and read images into the memory.
start_cogent;


cgflip(bgcol/256,bgcol/256,bgcol/256);
cgflip(bgcol/256,bgcol/256,bgcol/256);

num_segments=48;
con_im=ones(1024,768);
hh=zeros(1024,768);
disp('%%% loading images ...');
for i=1:num_images%length(files)
    logstring(sprintf('Loading image %i/%i', i, num_images))
    for ll=1:num_segments;
        con_im(indic{rvec(ll)})=round(t_4(ll,i+50)*3)+1;
    end,
    hh(:,:)=bgcol+PalContrasts(con_im).*ch_board;
    cgloadarray(10+(2*i)-1,1024,768,repmat(round(reshape(hh,resolution,1))/256',1,3));
    cgloadarray(10+2*i,1024,768,repmat((white+1-round(reshape(hh,resolution,1)))/256',1,3));
    
    if show_images_during_startup
        % show all images, nothing to log here (already in GPU memory)
        
        % prepare and show normal image
        cgdrawsprite(10+(2*i)-1,0,0); 
        cgflip;
        
        % prepare and show inverted image
        cgdrawsprite(10+(2*i),0,0); 
        cgflip;
    end
end

for ll=1:num_segments;
    con_im(indic{ll})=mean_con;
end
hh(:,:)=bgcol+con_im.*ch_board(:,:);
cgloadarray(1,1024,768,repmat(round(reshape(hh,resolution,1))/256',1,3));
cgloadarray(2,1024,768,repmat((white+1-round(reshape(hh,resolution,1)))/256',1,3));

clear hh con_im
clear picture

%%%%%%%%%%%%%%%%%%landolt c's zeichnen%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_landolt_c=3;
landolt_c=[-size_landolt_c size_landolt_c; size_landolt_c size_landolt_c; size_landolt_c -size_landolt_c; -size_landolt_c -size_landolt_c];
draw_landold_c(landolt_c,bgcol/256,white_landolt/256);

cgsetsprite(0);
cgflip(bgcol/256,bgcol/256,bgcol/256);

% Define landolt C code (for log file)
landolt_open(1:9) = '?'; % default: here only square, left and right should be shown, so we only define these below
landolt_open(9) = 'S'; % square
landolt_open(5) = 'L'; % left
landolt_open(7) = 'R'; % left

% init result
res.all_sprites_t_bckgrndsprt_imgnr_ldlt = []; % save t background_sprite, imagenr, landolt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cgtext('READY...',0,0);
cgflip(bgcol/256,bgcol/256,bgcol/256);

if triggermode==1
    display('Waiting for first serial scanner pulse (define details in code below)')
    clearserialbytes(res.port_num);
    waitserialbyte(res.port_num,inf,53);
else
    display('Waiting for any key to start (if this should be a scanner pulse, define below)')
    waitkeydown(inf);
end;

show_fix=0;
if useSerial
    clearserialbytes(res.port_num);
else
    clearkeys
end
t0=time;
res.t0=t0;
disp('Experiment running ...');
res.para.tstart=t0;
start_next=0;
show_fix=0;
z=0;
f=0;
img_counter=0;
flicker=1;

% Show standard lowc checkboard and landolt task
last_landolt = 0; % init
while time<(t0+t_dummy)
    tflip=time;
    z_old=z;
    z=ceil((time-t0)/t_fix_change);
    flicker=not(flicker);
    if flicker
        background_sprt = 1; % first standard background (low contrast checkboard, no highlighted areas)
    else
        cgdrawsprite(2,0,0); % inverted first standard background (same like above, but black is white and v.v.)
        background_sprt = 2;
    end
    cgdrawsprite(background_sprt,0,0);
    
    if z_old~=z
        if show_fix
            show_fix=0;
        else
            f=f+1;
            show_fix=1;
        end
    end
    if show_fix
        landolt = fixation(f); % fixation = 5: Landolt C left open, = 7: right open; 
    else
        landolt = 9;% 9: Landolt C closed (square)
    end
    cgdrawsprite(landolt,0,0); 
    
    waituntil(tflip+100);
    t = cgflip(bgcol/256,bgcol/256,bgcol/256);
    % only log if something changed
    if last_landolt ~= landolt
        res = loginput(res, useSerial); % read input from last run
        logstring(['Landolt C: ' int2str(landolt) ': ' landolt_open(landolt) '. Standard lowc checkboard background']); % log new run
        last_landolt = landolt;
    end
    % save image onset and details
    res.all_sprites_t_bckgrndsprt_imgnr_ldlt(end+1, :) = [t, background_sprt, img_counter, landolt];
end

% Show Landolt + Background images (circular checkerboard)
last_img_counter = 0; % init
t_start_images=time;
while time<t0+t_total-t_gray
    img_counter=img_counter+1;
    logstring(sprintf('Preparing image number : %03.0f',img_counter));
    start_next=t_start_images+img_counter*tot_image;
    
    res.para.timage(img_counter)=time-t0;
    while time<start_next
        tflip=time;
        flicker=not(flicker);
        if flicker
            background_sprt = 10+img_counter*2-1;
        else
            background_sprt = 10+img_counter*2;
        end
        cgdrawsprite(background_sprt,0,0);        
        
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
            landolt = fixation(f); % fixation = 5: Landolt C left open, = 7: right open
        else
            landolt = 9; % fixation = 9: Landolt C closed (square)
        end
        cgdrawsprite(landolt,0,0);
        
        tchange(n)=time;
        n=n+1;
        waituntil(tflip+100);
        t = cgflip(bgcol/265,bgcol/256,bgcol/256);
        % only log if something changed
        if last_landolt ~= landolt || last_img_counter ~= img_counter
            res = loginput(res, useSerial); % read input from last run
            logstring(['Landolt C: ' int2str(landolt) ': ' landolt_open(landolt) '. Image: ' int2str(img_counter)]); % log new run
            last_landolt = landolt;
            last_img_counter = img_counter;
        end
        % save image onset and details
        res.all_sprites_t_bckgrndsprt_imgnr_ldlt(end+1, :) = [t, background_sprt, img_counter, landolt];
    end
    res.timeimg=time;
end;

% Show Landolt with standard lowc checkerboard
t_end_images=time;
img_counter = -1; % marker for last image
while time<(t0+t_total)
    tflip=time;
    z_old=z;
    z=ceil((time-t0)/t_fix_change);
    flicker=not(flicker);
    if flicker
        background_sprt = 1;
    else
        background_sprt = 2;
    end
    cgdrawsprite(background_sprt,0,0);
    
    if z_old~=z
        if show_fix
            show_fix=0;
        else
            f=f+1;
            show_fix=1;
        end
    end
    if show_fix
        landolt = fixation(f); % fixation = 5: Landolt C left open, = 7: right open
    else
        landolt = 9; % fixation = 9: Landolt C closed (square)
    end
    cgdrawsprite(landolt,0,0);
    
    waituntil(tflip+100);
    % only log if something changed
    cgflip(bgcol/256,bgcol/256,bgcol/256);
    if last_landolt ~= landolt
        res = loginput(res, useSerial); % read input from last run
        logstring(['Landolt C: ' int2str(landolt) ': ' landolt_open(landolt) '. background']); % log new run
        last_landolt = landolt;
    end
    % save image onset and details
    res.all_sprites_t_bckgrndsprt_imgnr_ldlt(end+1, :) = [t, background_sprt, img_counter, landolt];
end

t_end=time;
res.para.tend=t_end;
% add the input to log (either keys or serial input)
res = loginput(res, useSerial);
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

% function to read input
function res = loginput(res, useSerial)
    if useSerial % !!! read out of button presses only works with serial in this scirpt.
        readserialbytes(res.port_num);
        logserialbytes(res.port_num);
        [k_value,k_time,n]=getserialbytes(res.port_num);
    else
        readkeys
        logkeys
        [k_value,k_time,n]=getkeydown;
    end

    if ~isfield(res, 'results')
        % init
        res.results.input_values = [];
        res.results.input_times = [];
    end

    if ~isequal(size(k_value), size(k_time))
        error('Check size')
    end
    
    res.results.input_values=[res.results.input_values; k_value(:)];
    res.results.input_times=[res.results.input_times; k_time-res.t0(:)];
end

function filename=my_config_cogent(subid)


end


function draw_landold_c(landolt_c,bgcol,white)

sprite_size=20;
cgmakesprite(5,sprite_size,sprite_size,bgcol,bgcol,bgcol);
cgsetsprite(5);
cgpencol(bgcol,bgcol,bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white,white,white);
i=[1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[4 3 2 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(5,'n');

cgmakesprite(6,sprite_size,sprite_size,bgcol,bgcol,bgcol);
cgsetsprite(6);
cgpencol(bgcol,bgcol,bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white,white,white);
i=[2 3 4 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[1 4 3 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(6,'n');

cgmakesprite(7,sprite_size,sprite_size,bgcol,bgcol,bgcol);
cgsetsprite(7);
cgpencol(bgcol,bgcol,bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white,white,white);
i=[3 4 1 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(7,'n');

cgmakesprite(8,sprite_size,sprite_size,bgcol,bgcol,bgcol);
cgsetsprite(8);
cgpencol(bgcol,bgcol,bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white,white,white);
i=[4 1 2 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
%cgtrncol(8,'n');

cgmakesprite(9,sprite_size,sprite_size,bgcol,bgcol,bgcol);
cgsetsprite(9);
cgpencol(bgcol,bgcol,bgcol);
%cgellipse(0,0,30,30,'f');
cgpencol(white,white,white);
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