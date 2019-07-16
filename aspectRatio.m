function [ AR, R_tMaj  ] = aspectRatio( tj_A, Rj_A, mP)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here


y_bar = mP(2);
x_bar = mP(1);
% figure(1)
% subplot(1,2,2)
% plot(tj_A,Rj_A,'-o')

%%% step through each angle and find a second angle that is closest to being
%%% 180 degrees away. When you've come close to the end for the second then
%%% stop

t_Diff = abs(tj_A - (tj_A(1)+pi));
elD_init = find(min(t_Diff)==t_Diff);
elD_init  = elD_init (1);

for ct = 1:elD_init-1

    t_Diff = abs(tj_A - (tj_A(ct)+pi));
    elD = find(min(t_Diff)==t_Diff);
    elD = elD(1);
    R_t(ct) = Rj_A(ct)+Rj_A(elD);              %store the sum of the radius for each combination of angles
    elD_A(ct) = elD;                         %store the element that is 180 degrees aways
end

mrt = max(R_t); % find the maximum radius within the array
el1M = find(mrt==R_t);
el1M = el1M(1);
el2 = elD_A(el1M);   % go back and retrieve the second point where the greatest distance is a maximum
%     radDiff = tj_A(el2)-tj_A(el1M) % here'e a check that they are are 180 apart

if tj_A(el1M)>(pi/2)

    t_Diff2 = abs(tj_A - (tj_A(el1M)-(pi/2)));  % using the second point, find a point that is 90 degrees away
    elpl90 = find(min(t_Diff2)==t_Diff2);
    elpl90 = elpl90(1);
    t_Diff3 = abs(tj_A - (tj_A(elpl90)+(pi)));  % using the second point, find a point that is 90 degrees away
    elmn90 = find(min(t_Diff3)==t_Diff3);
    elmn90 = elmn90(1);

else

    t_Diff2 = abs(tj_A - (tj_A(el1M)+(pi/2)));  % using the second point, find a point that is 90 degrees away
    elpl90 = find(min(t_Diff2)==t_Diff2);
    elpl90 = elpl90(1);
    t_Diff3 = abs(tj_A - (tj_A(elpl90)+(pi)));  % using the second point, find a point that is 90 degrees away
    elmn90 = find(min(t_Diff3)==t_Diff3);
    elmn90 = elmn90(1);

end

R_tMin = Rj_A(elpl90)+Rj_A(elmn90);
R_tMaj = mrt;

AR = R_tMaj/R_tMin;
% 
% figure(1)
% subplot(2,2,1)
elA = [el1M, el2];

xE = (Rj_A(elA).*cos(tj_A(elA)))+x_bar;
yE = (Rj_A(elA).*sin(tj_A(elA)))+y_bar;
% plot(xE,yE,'k-s','markerfacecolor','r','linewidth',2,'markersize',10)

elA = [elmn90, elpl90];
xE = (Rj_A(elA).*cos(tj_A(elA)))+x_bar;
yE = (Rj_A(elA).*sin(tj_A(elA)))+y_bar;
% plot(xE,yE,'k-s','markerfacecolor','y','linewidth',2,'markersize',10)

% axis equal
% 
% subplot(2,2,2)
% hold on
% plot(tj_A([el1M, el2]),Rj_A([el1M, el2]),'ks','markerfacecolor','r','linewidth',2,'markersize',10)
% hold on
% plot(tj_A([elmn90, elpl90]),Rj_A([elmn90, elpl90]),'ks','markerfacecolor','y','linewidth',2,'markersize',10)
%     if ct >= 1000
%         ['The el1 is ' num2str(el1M) ' and its angle is ' num2str(tj_A(el1M)*180/pi) ' and its radius is ' num2str(Rj_A(el1M)) ]
%         ['The el1 is ' num2str(el2) ' and its angle is ' num2str(tj_A(el2)*180/pi) ' and its radius is ' num2str(Rj_A(el2)) ]
%         ['The sum of the radii are ' num2str(Rj_A(el1M)+Rj_A(el2))  ]        
%         ['The el1 is ' num2str(elmn90) ' and its angle is ' num2str(tj_A(elmn90)*180/pi) ' and its radius is ' num2str(Rj_A(elmn90)) ]
%         ['The el1 is ' num2str(elpl90) ' and its angle is ' num2str(tj_A(elpl90)*180/pi) ' and its radius is ' num2str(Rj_A(elpl90)) ]
%        ['The sum of the small radii are ' num2str(Rj_A(elmn90)+Rj_A(elpl90))  ]        
% 
%  keyboard
%     end    
% pause
% end
% AR
% keyboard

end

