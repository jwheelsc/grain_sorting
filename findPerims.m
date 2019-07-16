function [  ptclP , cvxhP, arrN, arrNe ] = findPerims( rc )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    x = (rc(:,2));
    y = (rc(:,1));
    xR = max(x)-min(x);
    yR = max(y)-min(y);
    xN = x-min(x)+1;
    yN = y-min(y)+1;
    arrN = [xN,yN];

    newI = zeros(int32(xR+4),int32(yR+4));

    xOS = [-0.5,0.5,0.5,-0.5];
    yOS = [-0.5,-0.5,0.5,0.5];
    newAx = zeros(length(xN*4),1);
    newAy = zeros(length(xN*4),1);
    for t = 1:length(x)
        newAx([1:4]+(4*(t-1))) = xN(t)+xOS;
        newAy([1:4]+(4*(t-1))) = yN(t)+yOS;
    end

    xNe = newAx+0.5;
    yNe = newAy+0.5;
    arrNe = [xNe,yNe];

    for t = 1:length(xNe)
        newI(int32(xNe(t)+1),int32(yNe(t)+1)) = 1;
    end

    rBP = min(find(newI(:,2)==1));
    contour = bwtraceboundary(newI,[rBP 2],'W',8,Inf,'counterclockwise');
%             figure
%             plot(xN+1.5,yN+1.5,'+')
%             hold on
%             plot(xNe+1,yNe+1,'mo')
%             hold on
%             plot(contour(:,1),contour(:,2),'k-o')
    [kC] = convhull(double(xNe),double(yNe));
%             plot(xNe(kC)+1,yNe(kC)+1,'r-o')
     cvxhP = [xNe(kC),yNe(kC)];
     ptclP = [contour(:,1)-1,contour(:,2)-1];


end

