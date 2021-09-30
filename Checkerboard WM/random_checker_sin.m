function [ch_board, indic, ch_board_parts, ch_board_min]=random_checker_sin(n_circ,n_ang,r_scaling,scale, circle_param, angle_param, rect, visAngleMin, visAngleMax, distFromScreen, ppCM, numParts)

% define some parameters.


r_min=round(distFromScreen*tan(degtorad(visAngleMin))*ppCM); %original 20
r_max=round(distFromScreen*tan(degtorad(visAngleMax))*ppCM); %original 400
r_scaling=2;
% circle_param = 4; 
% angle_param = 4; 

invert=0; 
if numParts==2
    invert=1;
end

imw=rect(3);
imh=rect(4);

if invert==1
    imw=rect(4);
    imh=rect(3);
end


[XX,YY]=meshgrid(1:imw,1:imh);
[TH,RR]=cart2pol(XX-imw/2,YY-imh/2);
TH=TH+pi;
%create the checkerboard

ru=linspace(r_min^(1/r_scaling),r_max^(1/r_scaling),n_circ*2+1);
r_values=ru.^r_scaling;
th_values=(0:1/(n_ang*2):1)*2*pi;
th_values(end)=2*pi+0.01;

percent_thickness = 0.01;
r_values_nill = [];
for r_nill = 1+circle_param:circle_param:length(r_values)-circle_param
    r_values_nill = [r_values_nill,r_values(r_nill)-(r_values(r_nill)*percent_thickness+10), r_values(r_nill)+(r_values(r_nill)*percent_thickness+10)];
end

th_values_nill = [];
for t_nill = 1+angle_param:angle_param:length(th_values)
%     th_values_nill = [th_values_nill,th_values(t_nill)-th_values(t_nill)*percent_thickness, th_values(t_nill)+th_values(t_nill)*percent_thickness];
    th_values_nill = [th_values_nill,th_values(t_nill)-0.1, th_values(t_nill)+0.1];
%    th_values_nill = th_values_nill-pi/2; 

end


VV=sin(TH*n_ang);
UU=RR.^(1/r_scaling);
UU=UU-r_min^(1/r_scaling);
UU=sin(UU/(r_max^(1/r_scaling)-r_min^(1/r_scaling))*n_circ*2*pi);

UU_diff=sin(UU/(r_max^(1/r_scaling)-r_min^(1/r_scaling))*2*pi);

ch_board_min = UU_diff.*VV;
ii = find(0.9*r_min<RR);
ch_board_min(ii) = 0;
ii = find(0.2*r_min>RR);
ch_board_min(ii) =0; 
ch_board_min=sign(ch_board_min).*abs(ch_board_min).^(1/scale);

ii=find(RR<r_min);
UU(ii)=0;
ii=find(RR>r_max);
UU(ii)=0;

% UU_inverse = zeros(size(UU)); 

for r_nil2 = 1:2:length(r_values_nill)
    ii = find(r_values_nill(r_nil2)<RR & RR<r_values_nill(r_nil2+1));
    UU(ii)=0;
%     UU_inverse(ii)=1;
end

for t_nil2 = 1:2:length(th_values_nill)
    if t_nil2==length(th_values_nill)-1
        ii_1 = find(th_values_nill(t_nil2)<TH);
        ii_2 = find(TH<th_values_nill(t_nil2+1)-2*pi);   
        initial_val = th_values_nill(t_nil2+1)-2*pi; 
        ii = [ii_1;ii_2];
    else
        ii = find(th_values_nill(t_nil2)<TH & TH<th_values_nill(t_nil2+1));
    end
    UU(ii)=0;
%     UU_inverse(ii)=1;   
end

% ii=find(RR>r_max);
% UU_inverse(ii)=0;
% 
% ii=find(RR<r_min);
% UU_inverse(ii)=0;
% 
% ii=find(r_max-0.01*r_max<RR);
% UU_inverse(ii)=1; 
% 
% ii= find(r_max+0.01*r_max<RR);
% UU_inverse(ii)=0;

ch_board=UU.*VV;

if invert==1
    ch_board = rot90(ch_board);
end

ch_board=sign(ch_board).*abs(ch_board).^(1/scale);



indic = {};
ch_board_parts = {};
r_values_nill2 = [r_min, r_values_nill, r_max];
th_values_nill2 = [initial_val, th_values_nill];
for circ = 1:2:length(r_values_nill2)
    for ang = 1:2:length(th_values_nill)
        UU_new = UU; 
        ii=find(RR<r_values_nill2(circ));
        UU_new(ii)=0;
                
        ii=find(RR>r_values_nill2(circ+1));
        UU_new(ii)=0;
        
        UU_new_inverse = UU_new;
        ii=find(0.90*r_values_nill2(circ)<RR & RR<r_values_nill2(circ));
        UU_new_inverse(ii) = 5; 

        ii=find(RR>r_values_nill2(circ+1)& 1.1*r_values_nill2(circ+1)>RR);
        UU_new_inverse(ii) = 5; 
        
       
        % angles
        ii=find(TH<th_values_nill2(ang));
        UU_new(ii)=0;
        UU_new_inverse(ii)=0;
                
        ii=find(TH>th_values_nill2(ang+1));
        UU_new(ii)=0;    
        UU_new_inverse(ii)=0;
        
        ii=find(TH>th_values_nill2(ang)-0.5*th_values_nill2(1) & TH<th_values_nill2(ang));
        UU_new_inverse(ii) =5; 
        
        ii=find(TH<th_values_nill2(ang+1)+0.5*th_values_nill2(1) & TH>th_values_nill2(ang+1));
        UU_new_inverse(ii) =5; 
        
        ii=find(0.9*r_values_nill2(circ)>RR);
        UU_new_inverse(ii) = 0; 

        ii=find(1.1*r_values_nill2(circ+1)<RR);
        UU_new_inverse(ii) = 0; 
        
        UU_new_inverse(UU_new_inverse~=5) = 0; 
        UU_new_inverse(UU_new_inverse==5) = 1;

        
             
        ch_board_new=UU_new.*VV;
        ch_board_new=sign(ch_board_new).*abs(ch_board_new).^(1/scale);
        if invert==1
            ch_board_new = rot90(ch_board_new);
            UU_new_inverse = rot90(UU_new_inverse);
        end
%         figure
%         imagesc(UU_new_inverse)
        [indic_r, indic_c] = find(ch_board_new);
        ch_board_parts = [ch_board_parts, UU_new_inverse];
        indic = [indic,[indic_r,indic_c]];
    end
end

end