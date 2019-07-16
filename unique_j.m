function [ X ] = unique_j( X )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

x = X(:,1);
y = X(:,2);
ct = 1;
ctA = [];
for i = 1:length(x)-1
    xi = x(i);
    yi = y(i);
    
    for j = i+1:length(x)
       xj = x(j);
       yj = y(j);
       if (abs(xj - xi))<1e-2 & (abs(yj - yi))<1e-2
           ctA(ct) = j;
           ct = ct+1 ;
       end
    end
end

X(ctA,:) = [];