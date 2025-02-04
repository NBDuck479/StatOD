function [statNumOb, measDelta] = StationObs(observedMeas, computedMeas, visibilityMask, i)
% Made a function to paste into each filter

% --- Observation Deviation --- 
    % [range from stations] [range dot from stations]
    OC = observedMeas(i,:) - computedMeas(i,:);
    
    % range first row; range rate second
    measDelta = [OC(1:3); OC(4:6)];
    
    % --- Observation State Matrix
    % Determine which statoin made the ob and do calculations just for that
    stationOb = isnan(visibilityMask(i,:));
    
    % index where it is NOT Nan is the station number that made the ob
    statNumOb = find(stationOb == 0);
    
    % put error if mulipltle observation sat the same time!
    assert(length(statNumOb) == 1 || isempty(statNumOb), 'Multiple Station Obs at same time!!!')
end