% this script was used to take data from compute shape and unroll the
% particle  in polar coordinates to get rid of double values functios, but
% it was getting too hard to make that work, so i'll eliminate the fourier
% analysis all together. 

function[tj_A, Rj_A] = unrollParticle(mP,X)

y_bar = mP(2);
x_bar = mP(1);

d_mP = mP-X;  % find the distance from the center of mass to each perimeter point
R_j = sqrt(sum(d_mP.^2,2));
theta_j = atan(((X(:,2)-y_bar)./(X(:,1)-x_bar)));  % find the angle from the center of mass to each perimeter point4


dy = (X(:,2)-y_bar);
dx = (X(:,1)-x_bar);
s = dy./dx;

el = [];
el(:,1) = dy>0 & dx>0;
el(:,2) = dy>0 & dx<0;
el(:,3) = dy<0 & dx<0;
el(:,4) = dy<0 & dx>0;


theta_j = atan(s);

Rj_A = [];
tj_A = [];
pPlus = [0,pi,pi,2*pi];

for t_i = 1:4

    tj = [];
    is = [];
    Rj = [];
    [tj,is] = sort(theta_j(logical(el(:,t_i))),'ascend');
    %                 tj = theta_j(logical(el(:,t_i)))
    tj = tj+pPlus(t_i);
    Rj = R_j(logical(el(:,t_i)));
    Rj = Rj(is);
    Rj_A = [Rj_A;Rj];
    tj_A = [tj_A;tj];
    
end
