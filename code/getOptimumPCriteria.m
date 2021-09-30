function optimumPs = getOptimumPCriteria(y,discardSampleSize,sampleSize,et,a0,y0,pMax,showPlots)
%runSingleCaseScenarios Summary of this function goes here
%initialize our arrays

AIC = zeros(pMax,1); FPE = zeros(pMax,1); SIC= zeros(pMax,1); 
HQC= zeros(pMax,1); BIC= zeros(pMax,1); sp2= zeros(pMax,1); 

for p = 1:pMax
    [AIC(p),FPE(p),SIC(p),HQC(p),BIC(p),sp2(p)] =  getCriteria(y,p,discardSampleSize,sampleSize,et,a0,y0);
end
[~,I] = min(AIC); optimumPs.AIC = I;
[~,I] = min(FPE); optimumPs.FPE = I;
[~,I] = min(SIC); optimumPs.SIC = I;
[~,I] = min(HQC); optimumPs.HQC = I;
[~,I] = min(BIC); optimumPs.BIC = I;
[~,I] = min(sp2); optimumPs.sp2 = I;

if ~showPlots
    return
end

ps = 1:pMax;
subplot(2,3,1)
plot(ps,AIC,'b');set(gca,'XTick',1:20);grid on; title('AIC (Akaike information)')
subplot(2,3,2)
plot(ps,FPE,'r');set(gca,'XTick',1:20);grid on; title('FPE (Final prediction)')
subplot(2,3,3)
plot(ps,SIC,'g');set(gca,'XTick',1:20);grid on; title('SIC (Schwarz information)')
subplot(2,3,4)
plot(ps,HQC,'k');set(gca,'XTick',1:20);grid on; title('HQC (Hannan-Quinn)')
subplot(2,3,5)
plot(ps,BIC,'k');set(gca,'XTick',1:20);grid on; title('BIC (Bayesian information)')
subplot(2,3,6)
plot(ps,sp2,'k');set(gca,'XTick',1:20);grid on; title('(Just) ó_p^2')

end

