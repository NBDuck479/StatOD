function [stationECI] = StationStateECI(stations, theta0, Re, tSpan)
% This function computes the state of each station
%
% --- Inputs ---
% INPUT IS ECEF
% stationsLatLong:  struct with lat and long for each station
% theta0:           Initial angle between ECI and ECEF
% Re:               Radius of the Earth
%
% Units - Deg, km, km/s
%
% --- Outputs ---
% OUTPUT WILL BE ECI

% Number of stations
numStations = length(fields(stations)); % Get the number of stations dynamically

% Earth rotation rate
earthRotRad = (2*pi) / (24*60*60); % rad/sec
earthRotDeg = 360 / (24*60*60); % deg/sec

% Predefine the station states (in ECI) structure
stationECI = cell(numStations, length(tSpan)); 


for t = 1:length(tSpan)
    
    for i = 1:numStations
        
        % the station lat and long (dynamically accessing fields)
        lat = stations.(['stat' num2str(i)]).lat;
        long = stations.(['stat' num2str(i)]).long;
        
        % convert from spherical to cartesian coordinates (ECEF)
        [x, y, z] = sph2cart(deg2rad(long), deg2rad(lat), Re);  
        
        % Station position in ECEF frame
        stationPosECEF = [x; y; z];
        
        % --- Compute station velocity
        % Project station position onto XY plane
        stationXYProj = [stationPosECEF(1:2); 0];
        
        % Get linear velocity of station
        stationVelMag = earthRotRad * norm(stationXYProj);
        
        % Velocity unit vector (cross product for direction)
        statVelUnitVec = cross([0; 0; 1], stationXYProj / (norm([0; 0; 1]) * norm(stationXYProj)));
        
        % Station velocity vector
        stationVelECEF = stationVelMag * statVelUnitVec;
        
        % --- Transform from ECEF to ECI
        % Rotation matrix for simple transformation
        Rz = @(Theta) [cosd(Theta), -sind(Theta), 0; sind(Theta), cosd(Theta), 0; 0, 0, 1];
        
        % Earth rotation angle at time t
        thetaCurrent = tSpan(t) * earthRotDeg + theta0;
        
        % Compute the station's ECI coordinates (apply rotation to ECEF)
        statPosECI = Rz(thetaCurrent) * stationPosECEF;
        
        statVelECi = Rz(thetaCurrent) * stationVelECEF;
        
        stationECI{i, t} = [statPosECI; statVelECi];
         
    end
    
end

stationECI = stationECI';

end