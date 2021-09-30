function [aEstimate,et_p,yModel] = getAEstimates(y,p,discardSampleSize,sampleSize,et,a0,y0)
    %GETAESTIMATES Summary of this function goes here
    %   Detailed explanation goes here
    %find the p-dependent a coefficients based on the least-squared approach
    %https://textbooks.math.gatech.edu/ila/least-squares.html
    
    %https://siva82kb.github.io/2018/09/least-square-estimation-of-ar-models-and-whitening-part-i

    %get the useful y values (after the discarded samples)
    %the intercept is a given constant in this case (a0)
    totalSampleSize = sampleSize + discardSampleSize;
    b = y(discardSampleSize+1:end)-a0;
    A = zeros(sampleSize,p);
    for i=discardSampleSize+1:totalSampleSize
        irow = i-discardSampleSize;
        for j = 1:p
            A(irow,j)=y(i-j);
        end
    end
    %find the a coefficients (using the least squares method)
    aEstimate = ((A'*A)\A')*b;

    %find the y values based on the estimated model
    yModel = getYDataset(aEstimate,discardSampleSize,sampleSize,et,a0,y0);

    %find the new p-dependent-error estimates
    et_p = yModel - y;

    %get only the errors for the main sample
    %et_p = et_p(discardSampleSize+1:end);
end

