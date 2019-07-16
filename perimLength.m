function [ lgth ] = perimLength(x,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%     msfc = 25/38;
msfc = 1
    for i=1:length(x)-1;
        d(i) = sqrt((x(i+1)-x(i)).^2+(y(i+1)-y(i)).^2);
    end
    lgth = sum(d)*msfc;


end

