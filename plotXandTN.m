function [] = plotXandTN( x, tn, Nc, t0, tf )
%PLOTXANDTN Summary of this function goes here
%   Detailed explanation goes here

subplot(2,1,1);
plot(1:size(tn,2),tn);
subplot(2,1,2);

yDisp = repmat([1,1,0,0],size(x));
xDisp = zeros(1, size(x,2)*2);
xDisp(1:2:end) = 24*3600/Nc - x;
xDisp(2:2:end) = x;
xDisp = cumsum(xDisp);
xDisp = xDisp - x(1);
xDisp = repelem2014(xDisp,repmat(2, size(xDisp)));
xDisp = circshift(xDisp, [0 -1]);
xDisp(end) = tf;
%xDisp(1:2:end) = xDisp(1:2:end) + 0.01; % Add offset, so we have different values
plot(xDisp,yDisp);
axis([t0, tf, -0.5, 1.5]);

end

