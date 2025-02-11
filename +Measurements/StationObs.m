function [statNumOb] = StationObs(visibilityMask, i)
% Determines which station made the observation

    % --- Observation State Matrix
    % Determine which statoin made the ob and do calculations just for that
    stationOb = isnan(visibilityMask(i,:));
    
    % index where it is NOT Nan is the station number that made the ob
    statNumOb = find(stationOb == 0);
    
    % put error if mulipltle observation sat the same time!
    assert(length(statNumOb) == 1 || isempty(statNumOb), 'Multiple Station Obs at same time!!!')
end