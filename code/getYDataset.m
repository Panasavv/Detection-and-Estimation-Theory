function y = getYDataset(a, discardSampleSize, sampleSize, et, a0, y0)
%GETYS Summary of this function goes here
%   Detailed explanation goes here
%the discard sample size must be greater than p

p = length(a);

totalSampleSize = sampleSize + discardSampleSize;

%initialize y array
y = zeros(totalSampleSize, 1);

% y0
% y1 = a0 + a1*y0; 
% y2 = a0 + a1*y1 + a2*y0
% y3 = a0 + a1*y2 + a2*y1 + a3*y0 <--last "special" becauses it uses y0
%first p y values are special (y0 is not used at all)
for yi=1:p
    %disp(['yi=' num2str(yi)]);
    y(yi) = a0;
    for ai = 1:yi - 1
       y(yi) = y(yi) + a(ai)*y(yi-ai);
    end
    %max a index is yi (see above)
    y(yi) = y(yi) + a(yi)*y0;

    %apply a "noise" value that follows normal distribution (et~N(0,errorSigma))
    y(yi) = y(yi) + et(yi);
end

%get the rest of y values
for yi=p+1:totalSampleSize
    y(yi) = a0;
    for ai = 1:p
       y(yi) = y(yi) + a(ai)*y(yi-ai);
    end
    
    %apply a "noise" value that follows normal distribution (et~N(0,errorSigma))
    y(yi) = y(yi) + et(yi);
end
end

