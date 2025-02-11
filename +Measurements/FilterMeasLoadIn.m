function [observedMeas] = FilterMeasLoadIn(yHist)
% I found myself copying and pasting this into several filters so I created
% a function to do it

% --- Measurements ---
% noisey measurements - these are from sensor
rangeMeas    = yHist.Range;
rangeDotMeas = yHist.RangeRate;

% Observed measurements (from sensor)
observedMeas = [rangeMeas, rangeDotMeas];
end
