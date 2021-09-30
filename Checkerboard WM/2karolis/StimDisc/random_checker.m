function [ch_board,indic]=random_checker(n_circ,n_ang,r_scaling);

% define some parameters.

n_fields=n_ang*n_circ;

r_min=20;
r_max=400;
% if nargin<3
r_scaling=2;
% end

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

ch_board=zeros(size(TH));
val=1;
indic_all = {};
for k=1:(n_circ*2) 
    for l=1:(n_ang*2) 
        ii=find(((TH>=th_values(l)-eps))&(TH<(th_values(l+1)+eps))&(RR<(r_values(k+1)+eps))&(RR>=(r_values(k))-eps)); 
        ch_board(ii)=(mod(val,2)-0.5)*2;
        indic_all = [indic_all,ii];
        val=val+1;
    end 
    val=val+1;
end

circles = [];

for circ = 1:(n_circ)
    circles = [circles,repmat(repelem(1:n_ang,2),1,2)+(circ-1)*n_ang];
end

indic = {};

for index =1:n_circ*n_ang
    ix = circles==index;
    indics_now = indic_all(ix);
    new_indic = [];
    for iix = 1:4
        new_indic = [new_indic;indics_now{iix}];
    end
    indic = [indic,new_indic];
end
        

end