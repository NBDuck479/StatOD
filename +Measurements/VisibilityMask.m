function [visibilityMask, viewingAngles, statNumOb, obTime] = VisibilityMask(stationPos, scPos, elevAngle, currTime, deltaT)
% This computes the visibilty mask for a given ground station via the
% elevation angle of the spacecraft relative to the ground station.
%
% --- Inputs ---
% stationPos:       {n x 3} cell for each ground station position
% scPos:            [3 x n] Pos of Single Spacecraft
% elevAngle:        [1 x 1] Elevation Angle Limit for Ground Station
% tSpan:            [1 x n] The time to run this calculation
%
% Units: Deg, km, km/s
%
% Assumption: EVERYTHING INPUT IS IN ECI

[~, numStations] = size(stationPos);

% Initialize - set these as empty for EKF
viewingAngles = [];
statNumOb = [];
obTime = [];


count = 1; 

for i = 1:length(currTime)
    
    for j = 1:numStations
        
        % grab position of single station at time step
        singleStatPos = stationPos{i,j};
        
        % Distance between ground station and spacecraft
        rho = scPos(:,i) - singleStatPos;
        
        % Local Vertical Unit Vector
        ehat = singleStatPos / norm(singleStatPos);
        
        % elevation angle of spacecraft relative to ground station
        theta = asind(dot(rho,ehat) / norm(rho));
        
        % Enforce elevation angle for visibility mask
        if theta > elevAngle
            % Spacecraft is high enough above horizon to be visible
            visibilityMask(i,j) = 1;
            
            % Store angles ground could see
            viewingAngles(i,j) = theta;
            
            % which station made the observation 
            % BREAKS DOWN IF MULTIPLE 
            statNumOb(i,j) = j; 

            % time of observation 
            obTime(i,j) = currTime(i);
             
            
        else
            % Spacecraft is not high enough to be visible
            visibilityMask(i,j) = NaN;
            
            % no viewing angles
            viewingAngles(i,j) = NaN;
            
            statNumOb(i,j) = j; 
            
            % time of observation 
            obTime(i,j) = NaN;
            
        end
        
    end
end