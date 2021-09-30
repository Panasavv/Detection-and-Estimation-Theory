function [AIC,FPE,SIC,HQC,BIC,sp2] = getCriteria(y,p,discardSampleSize,sampleSize,et,a0,y0)
%   getCriteria Summary of this function goes here
%   Detailed explanation goes here

    %get modeled y based on recalculation of a (p-dependent)
    [~,et_p, yModel] = getAEstimates(y,p,discardSampleSize,sampleSize,et,a0,y0);

    %omit discarded data
    yModel = yModel(discardSampleSize+1:end);
    et_p = et_p(discardSampleSize+1:end);

    T =  sampleSize;
    sp2 = sum(et_p(p:end).^2)/(T-p-1);

    %https://sccn.ucsd.edu/wiki/Chapter_3.5._Model_order_selection
    %note that the sign for the AIC in the paper is incorrect
    AIC =  2*T*log(sp2)+ 2*p;

    FPE = sp2 / (T - p) * (T + p);
    SIC = log(sp2)+ p * log(T) / T;
    HQC = log(sp2)+ 2 / T * p * log(log(T));
    BIC = (T - p) * log(T * sp2/(T - p))+T*(1+log(sqrt(2*pi)))+p*log((sum(yModel.^2)-T*sp2)/p);
end

