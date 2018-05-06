% Calculates a BVP, peak to peak (p2p) value and location from a PPG signal
% Author: David Enslev Nyrnberg
% Inputs:
%  - signal: time line of filtered PPG signal 
%  - minPeakHeight: threashhold for peaks to avoid BVP notches
%  - peakCut: percent scale for size of p2p distance
%  - printOut: (1) to enable function prints
% Outputs:
% - p2pBVPval: sequence of peak values
% - p2pBVPloc: locations of detected peaks on signal time line

function [p2pBVPval, p2pBVPloc] = PPG2PEAK(signal, minPeakHeight, peakCut, printOut)

% find all data points who are higher or equal to neighboring values with a
% threshhold to find relevandt values.
[valPeakBVP,locPeakBVP] = findpeaks(-signal, 'MinPeakHeight', minPeakHeight);

% 25% and 50% quantiles for the distance of peaks
diffStats = quantile(diff(locPeakBVP),[0.25, 0.5]);
if printOut == 1
    fprintf('Median: %.1f\n', diffStats(1))
    fprintf('Lower quantile: %.1f\n', diffStats(2))
end

% calculate peaks, which should be remove, for too low distance to next point.
failedPeak = diff(locPeakBVP) < floor(diffStats(1)*peakCut);
locSuccesPeak = find(failedPeak==0);
locPeakBVP2 = locPeakBVP(locSuccesPeak);
valPeakBVP2 = -valPeakBVP(locSuccesPeak);

% asign correct detected values and locations to output.
p2pBVPval = valPeakBVP2;
p2pBVPloc = locPeakBVP2;

end