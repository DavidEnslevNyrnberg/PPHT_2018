function [p2pBVPval, p2pBVPloc] = PPG2PEAK(signal, minPeakHeight, peakCut)

[valPeakBVP,locPeakBVP] = findpeaks(-signal, 'MinPeakHeight', minPeakHeight);

diffStats = quantile(diff(locPeakBVP),[0.25, 0.5]);
fprintf('Meadian: %.1f\n', diffStats(1))
fprintf('Lower quantile: %.1f\n', diffStats(2))

failedPeak = diff(locPeakBVP) < floor(diffStats(1)*peakCut);
locSuccesPeak = find(failedPeak==0);
locPeakBVP2 = locPeakBVP(locSuccesPeak);
valPeakBVP2 = -valPeakBVP(locSuccesPeak);

p2pBVPval = valPeakBVP2;
p2pBVPloc = locPeakBVP2;

end