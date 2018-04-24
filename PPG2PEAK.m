function [p2pBVPval, p2pBVPloc] = PPG2PEAK(signal, minPeakHeight, peakCut)

[valPeakBVP,locPeakBVP] = findpeaks(-signal,'MinPeakHeight',minPeakHeight);

diffStats = quantile(diff(locPeakBVP),[0.25]);
fprintf('output for lower quantile: %.1f\n',diffStats)

failedPeak = diff(locPeakBVP) < floor(diffStats*peakCut);
locSuccesPeak = find(failedPeak ==0);
locPeakBVP2 = locPeakBVP(locSuccesPeak);
valPeakBVP2 = -valPeakBVP(locSuccesPeak);

p2pBVPval = valPeakBVP2;
p2pBVPloc = locPeakBVP2;

end