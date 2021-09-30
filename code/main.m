function main

clear;clc

y0 = 0;
a0 = 0; %intercept term
errorSigma = 0.0001;

%the discard sample size must be greater than p
discardSampleSize = 100;

sampleSizes = [25 50 100 200 400 800 1600];
replications = 1000;

pMax = 20;
ps = [2 3 5];

for ip = 1:length(ps)
    p = ps(ip);
    disp(['p = ' num2str(p)])

    %initialize stats counters
    sampleSizesCount = length(sampleSizes);
    AICSuccesses = zeros(sampleSizesCount,1); 
    SICSuccesses = zeros(sampleSizesCount,1);
    FPESuccesses = zeros(sampleSizesCount,1);
    HQCSuccesses = zeros(sampleSizesCount,1);
    BICSuccesses = zeros(sampleSizesCount,1);

    AICUnderEstimates = zeros(sampleSizesCount,1); 
    SICUnderEstimates  = zeros(sampleSizesCount,1);
    FPEUnderEstimates  = zeros(sampleSizesCount,1);
    HQCUnderEstimates  = zeros(sampleSizesCount,1);
    BICUnderEstimates  = zeros(sampleSizesCount,1);

    AICOverEstimates = zeros(sampleSizesCount,1); 
    SICOverEstimates  = zeros(sampleSizesCount,1);
    FPEOverEstimates  = zeros(sampleSizesCount,1);
    HQCOverEstimates  = zeros(sampleSizesCount,1);
    BICOverEstimates  = zeros(sampleSizesCount,1);

    %for each sample size (S) run 1000 replications
    for iSampleSize = 1:sampleSizesCount
        sampleSize = sampleSizes(iSampleSize);
        disp(['Running for sample size ' num2str(sampleSize) '...'])
        for r = 1:replications

            %get a coeficients in the (-1,1) range so that
            %|a1 + a2 + a3+ ...| < 1 in order to enause stationary AR process
            a = -1 + 2 * rand(p,1);
            while (sum(abs(a)))>1
                a = -1 + 2 * rand(p,1);
            end
           

            %calculate "noise" values: "the same simulated series" should be modeled so
            %we keep the same noise values for each t
            et = normrnd(0,errorSigma, sampleSize + discardSampleSize, 1);

            %calculate y values using the original a parameters
            y = getYDataset(a,discardSampleSize,sampleSize,et,a0,y0);

            %run single case scenarios
            pOptimum = getOptimumPCriteria(y,discardSampleSize,sampleSize,et,a0,y0,pMax,false);
            %increment in the case of successful prediction
            if pOptimum.AIC == p
                AICSuccesses(iSampleSize) = AICSuccesses(iSampleSize) +  1;
            elseif pOptimum.AIC < p
                AICUnderEstimates(iSampleSize) = AICUnderEstimates(iSampleSize) +  1;
            else
                AICOverEstimates(iSampleSize) = AICOverEstimates(iSampleSize) +  1;
            end

            if pOptimum.SIC == p
                SICSuccesses(iSampleSize) = SICSuccesses(iSampleSize) +  1;
            elseif pOptimum.SIC < p
                SICUnderEstimates(iSampleSize) = SICUnderEstimates(iSampleSize) + 1;
            else
                SICOverEstimates(iSampleSize) = SICOverEstimates(iSampleSize) + 1;
            end

            if pOptimum.FPE == p
                FPESuccesses(iSampleSize) = FPESuccesses(iSampleSize) +  1;
            elseif pOptimum.FPE < p
                FPEUnderEstimates(iSampleSize) = FPEUnderEstimates(iSampleSize) +1;
            else
                FPEOverEstimates(iSampleSize) = FPEOverEstimates(iSampleSize) +1;
            end

            if pOptimum.HQC == p
                HQCSuccesses(iSampleSize) = HQCSuccesses(iSampleSize) + 1;
            elseif pOptimum.HQC < p
                HQCUnderEstimates(iSampleSize) = HQCUnderEstimates(iSampleSize) + 1;
            else
                HQCOverEstimates(iSampleSize) = HQCOverEstimates(iSampleSize) + 1;
            end

            if pOptimum.BIC == p
                BICSuccesses(iSampleSize) = BICSuccesses(iSampleSize) + 1;
            elseif pOptimum.BIC < p
                BICUnderEstimates(iSampleSize) = BICUnderEstimates(iSampleSize)+1;
            else
                BICOverEstimates(iSampleSize) = BICOverEstimates(iSampleSize)+1;
            end

        end
        
    end
    
    %convert frequency to relative frequency
    AICSuccesses = AICSuccesses / replications;
    SICSuccesses = SICSuccesses / replications;
    FPESuccesses = FPESuccesses / replications;
    HQCSuccesses = HQCSuccesses / replications;
    BICSuccesses = BICSuccesses / replications;
    AICUnderEstimates = AICUnderEstimates/replications; 
    SICUnderEstimates = SICUnderEstimates/replications;
    FPEUnderEstimates = FPEUnderEstimates/replications;
    HQCUnderEstimates = HQCUnderEstimates/replications;
    BICUnderEstimates = BICUnderEstimates/replications;
    AICOverEstimates = AICOverEstimates/replications; 
    SICOverEstimates = SICOverEstimates/replications;
    FPEOverEstimates = FPEOverEstimates/replications;
    HQCOverEstimates = HQCOverEstimates/replications;
    BICOverEstimates = BICOverEstimates/replications;

    figure(1)
    subplot(length(ps),3,3*(ip-1)+1)
    hold off
    semilogx(sampleSizes,AICSuccesses, 'b')
    hold on
    semilogx(sampleSizes,SICSuccesses,'r')
    hold on
    semilogx(sampleSizes,FPESuccesses,'g')
    hold on
    semilogx(sampleSizes,HQCSuccesses,'k')
    hold on
    semilogx(sampleSizes,BICSuccesses,'m')
    legend('AIC','SIC','FPE','HQC','BIC','Location','southeast')
    title(['Successes (p = ' num2str(p) ')'])
    xlabel('Sample size')
    ylabel('Relative frequency')
    ylim([0 1])

    subplot(length(ps),3,3*(ip-1)+2)
    hold off
    semilogx(sampleSizes,AICUnderEstimates, 'b')
    hold on
    semilogx(sampleSizes,SICUnderEstimates,'r')
    hold on
    semilogx(sampleSizes,FPEUnderEstimates,'g')
    hold on
    semilogx(sampleSizes,HQCUnderEstimates,'k')
    hold on
    semilogx(sampleSizes,BICUnderEstimates,'m')
    legend('AIC','SIC','FPE','HQC','BIC','Location','northeast')
    title(['Underestimations (p = ' num2str(p) ')'])
    xlabel('Sample size')
    ylabel('Relative frequency')
    ylim([0 1])
    
    subplot(length(ps),3,3*(ip-1)+3)
    hold off
    semilogx(sampleSizes,AICOverEstimates, 'b')
    hold on
    semilogx(sampleSizes,SICOverEstimates,'r')
    hold on
    semilogx(sampleSizes,FPEOverEstimates,'g')
    hold on
    semilogx(sampleSizes,HQCOverEstimates,'k')
    hold on
    semilogx(sampleSizes,BICOverEstimates,'m')
    legend('AIC','SIC','FPE','HQC','BIC','Location','northeast')
    title(['Overestimations (p = ' num2str(p) ')'])
    xlabel('Sample size')
    ylabel('Relative frequency')
    ylim([0 1])
    
disp([AICSuccesses SICSuccesses FPESuccesses HQCSuccesses BICSuccesses])
end

% function y = incrementCounters(p, calculatedP,counters)
%     y = counters;
%     if p<calculatedP
%         y(1)=y(1)+1;
%     elseif p==calculatedP
%         y(2)=y(2)+1;
%     else
%         y(3)=y(3)+1;
%     end
% end