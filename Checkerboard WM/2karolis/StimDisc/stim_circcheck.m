function stim_images(subid,triggermode)


images2use=randperm(100);
port_num=1;
dummy_Vols=5;
num_images=10;
luminance=0.8;
bgcol=luminance/2;

curdir=pwd;

clock_int=round(clock);
date_str=sprintf('%02.0f%02.0f%02.0f',clock_int(1)-2000,clock_int(2),clock_int(3));
time_str=sprintf('%02.0f%02.0f',clock_int(4),clock_int(5));

% different timing parameters.
t_gray=3000; % 30000
TR=1500;
t_dummy=dummy_Vols*TR+t_gray;

rise_interval=100;
fade=[0:10 ones(1,5)*10 10:-1:0]+10;
num_fade=length(fade);

image_time=rise_interval*num_fade;
tot_image=3000; % total time for flashing one image
t_total=num_images*tot_image+t_dummy+t_gray;%(Vols+dummy_Vols)*TR;
t_fix_change=800;

tchange=zeros(10000,1);n=1; %vector to save timing of images.
resolution=1024*768; % number of pixels per image.

%statistik von fixation überprüfen!!!!
fixation=fixation_random([165 167],t_total/t_fix_change);%t/t_image/2+1);

%% save all the data in the stuct res
% information about the person, experiment, date etc.
res.info.subj=subid;
res.info.paradigm='Checkerring, local contrast variation';
res.info.date=date_str;
res.info.time=time_str;

% information about the parameters of the experiment
res.para.fixations=fixation;
res.para.luminance=luminance;
res.para.timage=zeros(num_images,1);
res.para.t_fix_change=t_fix_change;
res.para.num_fade=num_fade;
res.para.rise_interval=rise_interval;
res.para.imag_order=images2use;
res.para.image_time=image_time;
res.para.t_total=t_total;

% configure cogent for the sepcific run
% and 
if ~isstr(subid)
    subid=num2str(subid);
end
c=clock;
filename=['./logs/Images_' subid '_' sprintf('%02d',c(1)) '_' sprintf('%02d',c(2)) '_' int2str(c(3)) '_' sprintf('%02d',c(4)) '_' sprintf('%02d',c(5)) '_' sprintf('%02d',fix(c(6)))];
config_display( 0, 3, [0 0 0], [luminance luminance luminance], 'Arial', 24, 5 ,32)
%config_keyboard;
config_keyboard(1000,5,'nonexclusive')

config_serial(port_num,19200,0,0,8);
config_log( [filename '.log']);
save([filename '.mat'],'res');

% start the cogent application and read images into the memory.
start_cogent;

cgflip(bgcol , bgcol , bgcol);
cgflip(bgcol , bgcol , bgcol);

load stim4_12;
num_segments=48;
con_im=zeros(768,1024);
hh=zeros(768,1024);
disp('%%% loading images ...');
for i=1:num_images%length(files)
    for ll=1:num_segments; 
        con_im(indic{ll})=t_4(ll,i+150); 
    end, 
    hh(:,:)=con_im(:,:).*ch_board(:,:);
    hh(:)=(hh(:)+1)/2;
    IM{i}=reshape(hh,resolution,1);
end
clear hh con_im
full_img_on=zeros(resolution,3);
full_img_off=zeros(resolution,3);
clear picture

%%%%%%%%%%%%%%%%%%landolt c's zeichnen%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
landolt_c=[-5 5; 5 5; 5 -5; -5 -5];
draw_landold_c(landolt_c,luminance);

cgsetsprite(0);
cgflip(bgcol , bgcol , bgcol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cgtext('READY...',0,0);
cgflip(bgcol , bgcol , bgcol);

if triggermode==1
    clearserialbytes(port_num);
    waitserialbyte(port_num,inf,53);
else
    waitkeydown(inf);
end;

show_fix=0;
clearserialbytes(1);
t0=time;
disp('Experiment running ...');
res.para.tstart=t0;
start_next=0;
show_fix=0;
z=0;
f=0;
img_counter=0;
while time<(t0+t_dummy)
    tflip=time;
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
        cgdrawsprite(169,0,0);
    end
    cgflip(bgcol , bgcol , bgcol);
    waituntil(tflip+100);
end
flicker=1;
t_start_images=time;
while time<t0+t_total-t_gray
    img_counter=img_counter+1;
    disp(sprintf('Showing image number : %03.0f',img_counter)); 
    start_next=t_start_images+img_counter*tot_image;

    full_img_on(:,:)=repmat(IM{img_counter},1,3);
    cgloadarray(11,1024,768,luminance*full_img_on);
   
    full_img_off(:,:)=1-repmat(IM{img_counter},1,3);
    cgloadarray(12,1024,768,luminance*full_img_off);
   
    res.para.timage(img_counter)=time-t0;
    while time<start_next
        tflip=time;
        flicker=xor(flicker,1)
            if flicker
        cgdrawsprite(11,0,0);
    else
        cgdrawsprite(12,0,0);
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
            cgdrawsprite(169,0,0);
        end

        tchange(n)=time;
        n=n+1;
        waituntil(tflip+100);
        cgflip(bgcol,bgcol,bgcol);
    end
    res.timeimg=time;

end;
t_end_images=time;
while time<(t0+t_total)
    tflip=time;
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
        cgdrawsprite(169,0,0);
    end
        
    cgflip(bgcol , bgcol , bgcol);
    waituntil(tflip+100);
end

t_end=time;
res.para.tend=t_end;
readserialbytes(1);
[k_value,k_time,n]=getserialbytes(1);
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

function filename=my_config_cogent_notebook(subid)
if ~isstr(subid)
    subid=num2str(subid);
end
c=clock;
filename=['./logs/NaturalScenes_' subid '_' int2str(c(3)) '_' int2str(c(2)) '_' int2str(c(1)) '_' int2str(c(4)) '_' int2str(c(5)) '_' int2str(c(6))];

config_display( 1, 3, [0 0 0], [luminance,luminance,luminance], 'Arial', 24, 4 ,32)
%config_keyboard;
config_keyboard(1000,5,'nonexclusive')

%config_serial(1,19200,0,0,8);
config_log( [filename '.log']);

end


function draw_landold_c(landolt_c,luminance)
bg_color=0;
bgcol=luminance/2;
sprite_size=40;
cgmakesprite(165,sprite_size,sprite_size,bg_color,bg_color,bg_color);
cgsetsprite(165);
cgpencol(bgcol , bgcol , bgcol);
cgellipse(0,0,30,30,'f');
cgpencol(luminance,luminance,luminance);
i=[1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[4 3 2 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
cgtrncol(165,'n');

cgmakesprite(166,sprite_size,sprite_size,bg_color,bg_color,bg_color);
cgsetsprite(166);
cgpencol(bgcol , bgcol , bgcol);
cgellipse(0,0,30,30,'f');
cgpencol(luminance,luminance,luminance);
i=[2 3 4 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[1 4 3 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
cgtrncol(166,'n');

cgmakesprite(167,sprite_size,sprite_size,bg_color,bg_color,bg_color);
cgsetsprite(167);
cgpencol(bgcol , bgcol , bgcol);
cgellipse(0,0,30,30,'f');
cgpencol(luminance,luminance,luminance);
i=[3 4 1 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
cgtrncol(167,'n');

cgmakesprite(168,sprite_size,sprite_size,bg_color,bg_color,bg_color);
cgsetsprite(168);
cgpencol(bgcol , bgcol , bgcol);
cgellipse(0,0,30,30,'f');
cgpencol(luminance,luminance,luminance);
i=[4 1 2 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
cgtrncol(168,'n');

cgmakesprite(169,sprite_size,sprite_size,bg_color,bg_color,bg_color);
cgsetsprite(169);
cgpencol(bgcol , bgcol , bgcol);
cgellipse(0,0,30,30,'f');
cgpencol(luminance,luminance,luminance);
i=[4 1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
cgtrncol(169,'n');

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
    while A(B(1))==r(i) && r(i)==r(i-1) && r(i)==r(i-2)
        B=randperm(length(A));
    end
    i=i+1;
    r(i)=A(B(1));
end
end