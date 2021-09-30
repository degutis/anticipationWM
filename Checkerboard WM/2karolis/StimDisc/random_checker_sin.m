function [ch_board, indic]=random_checker_sin(n_circ,n_ang,r_scaling,scale);

% define some parameters.


n_fields=n_ang*n_circ;
r_min=20; %original 20
r_max=400; %original 400
r_scaling=2;
circle_param = 4; 


imw=1024;
imh=768;
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

angle_param = 4; 
th_values_nill = [];
for t_nill = 1+angle_param:angle_param:length(th_values)
%     th_values_nill = [th_values_nill,th_values(t_nill)-th_values(t_nill)*percent_thickness, th_values(t_nill)+th_values(t_nill)*percent_thickness];
    th_values_nill = [th_values_nill,th_values(t_nill)-0.1, th_values(t_nill)+0.1];

end

VV=sin(TH*n_ang);
UU=RR.^(1/r_scaling);
UU=UU-r_min^(1/r_scaling);
UU=sin(UU/(r_max^(1/r_scaling)-r_min^(1/r_scaling))*4*2*pi);
ii=find(RR<r_min);
UU(ii)=0;
ii=find(RR>r_max);
UU(ii)=0;

for r_nil2 = 1:2:length(r_values_nill)
    ii = find(r_values_nill(r_nil2)<RR & RR<r_values_nill(r_nil2+1));
    UU(ii)=0;
end

for t_nil2 = 1:2:length(th_values_nill)
    if t_nil2==length(th_values_nill)-1
        ii_1 = find(th_values_nill(t_nil2)<TH);
        ii_2 = find(TH<th_values_nill(t_nil2+1)-2*pi);   
        ii = [ii_1;ii_2];
    else
        ii = find(th_values_nill(t_nil2)<TH & TH<th_values_nill(t_nil2+1));
    end
    UU(ii)=0;
end


ch_board=UU.*VV;
ch_board=sign(ch_board).*abs(ch_board).^(1/scale);

indic = {};
r_values_nill2 = [r_min, r_values_nill, r_max];
th_values_nill2 = [0, th_values_nill];
for circ = 1:2:length(r_values_nill2)
    for ang = 1:2:length(th_values_nill)
        UU_new = UU; 
        ii=find(RR<r_values_nill2(circ));
        UU_new(ii)=0;
        ii=find(RR>r_values_nill2(circ+1));
        UU_new(ii)=0;
        
        ii=find(TH<th_values_nill2(ang));
        UU_new(ii)=0;
        ii=find(TH>th_values_nill2(ang+1));
        UU_new(ii)=0;        
        

        ch_board_new=UU_new.*VV;
        ch_board_new=sign(ch_board_new).*abs(ch_board_new).^(1/scale);
%         figure
%         imagesc(ch_board_new)
        [indic_r, indic_c] = find(ch_board_new); 
        indic = [indic,[indic_r,indic_c]];
    end
end

end