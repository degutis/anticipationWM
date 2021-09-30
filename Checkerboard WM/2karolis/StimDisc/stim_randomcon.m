function rot_localizer(subid,triggermode)

contrast=0.35;
offset=0;
bgcol=contrast/2+offset;

port_num=1;

landolt_c=[-5 5; 5 5; 5 -5; -5 -5];
t_total=300000;  
t_image=1000;
t_dummy=2000*5;
time_cycle=30000;
steps=30;
tpp=time_cycle/steps; %time per picture
fixation=fixation_random([65 67],(t_total+t_dummy)/t_image/2+1);
hf=1/60/2;

if ~isstr(subid)
    subid=num2str(subid);
end
c=clock;
filename=['./logs/Random_Con_' subid '_' sprintf('%02d',c(1)) '_' sprintf('%02d',c(2)) '_' int2str(c(3)) '_' sprintf('%02d',c(4)) '_' sprintf('%02d',c(5)) '_' sprintf('%02d',fix(c(6)))];

config_display( 2, 3, [0 0 0], [contrast+offset contrast+offset contrast+offset], 'Arial', 24, 4 ,32)
%config_keyboard;
config_keyboard(1000,5,'nonexclusive')

config_serial(1,19200,0,0,8);
config_log( [filename '.log']);
start_cogent;
cgalign('c','c');

%%%%%%%%%%%%%%%%%%%bilder laden%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=0;
for i=1:steps
    z=z+1;
    a=double(imread(['./localizer/h' num2str(i) '.BMP']))'/255;
    cgloadarray(z,1024,768,repmat(a(:),1,3)*contrast+offset);
    z=z+1;
    a=double(imread(['./localizer/v' num2str(i) '.BMP']))'/255;
    cgloadarray(z,1024,768,repmat(a(:),1,3)*contrast+offset);
end
frames=z;

%%%%%%%%%%%%%%%%%%landolt c's zeichnen%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cgmakesprite(65,20,20,bgcol,bgcol,bgcol);
cgsetsprite(65);
i=[1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[4 3 2 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end

cgmakesprite(66,20,20,bgcol,bgcol,bgcol);
cgsetsprite(66);
i=[2 3 4 1];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[1 4 3 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end


cgmakesprite(67,20,20,bgcol,bgcol,bgcol);
cgsetsprite(67);
i=[3 4 1 2];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end

cgmakesprite(68,20,20,bgcol,bgcol,bgcol);
cgsetsprite(68);
i=[4 1 2 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end

cgmakesprite(69,20,20,bgcol,bgcol,bgcol);
cgsetsprite(69);
i=[4 1 2 3 4];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end
i=[3 2 1 4 3];
for j=1:length(i)-1
    cgdraw(landolt_c(i(j),1),landolt_c(i(j),2),landolt_c(i(j+1),1),landolt_c(i(j+1),2));
end

cgsetsprite(0);
cgflip(bgcol,bgcol,bgcol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cgtext('READY...',0,0);
cgflip(bgcol,bgcol,bgcol);
if triggermode==1
     clearserialbytes(port_num);
     waitserialbyte(port_num,inf,53);
else
    waitkeydown(inf);
end;
clearserialbytes(1);
t0=time;

show_fix=0;
z=0;
f=0;  
i=0;
flicker=0;    
while time<(t0+t_dummy)
    i=i+1;
    z_old=z;
    z=ceil((time-t0)/t_image);
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
            cgdrawsprite(69,0,0);
    end    
    cgflip(bgcol,bgcol,bgcol);
    waituntil(t0+i*100);
end

% show_fix=0;
% z=0;
% f=0;  
% i=0;
% flicker=0;    
while time<(t0+t_dummy+t_total)
    i=i+1;
    flicker=xor(flicker,1);    
    y=mod(time-t0-t_dummy,time_cycle);
    y=floor(y/tpp)*2+1+flicker;
    cgdrawsprite(y,0,0);
    z_old=z;
    z=ceil((time-t0)/t_image);
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
            cgdrawsprite(69,0,0);
    end    
    cgflip(bgcol,bgcol,bgcol);
    waituntil(t0+i*100);
end
    
    
t_end=time;
readserialbytes(1);
k_value=0;
k_time=0;
n=0;
k_time_corrected=0;
[k_value,k_time,n]=getserialbytes(1);
k_time_corrected=k_time-t0;

save([filename '.mat'],'k_value','k_time','k_time_corrected','t0','t_end','fixation');

stop_cogent;
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