function [stationECI] = StationStateECIproj1(stations, omegaEarth, Re, obTimes)


for t = 1:length(obTimes)
    
    % set current time 
    currTime = obTimes(t);
    
    % calcualte how far Earth has rotated
    theta = rad2deg(omegaEarth * currTime); 
    
    for i = 1:3
        % Station position in ECEF frame
        stationPosECEF = stations.(['stat' num2str(i) 'pos']);
        
        % --- Compute station velocity
        % Project station position onto XY plane
        stationXYProj = [stationPosECEF(1:2); 0];
        
        % Get linear velocity of station
        stationVelMag = omegaEarth * norm(stationXYProj);
        
        % Velocity unit vector (cross product for direction)
        statVelUnitVec = cross([0; 0; 1], stationXYProj / (norm([0; 0; 1]) * norm(stationXYProj)));
        
        % Station velocity vector
        stationVelECEF = stationVelMag * statVelUnitVec;
        
        % --- Transform from ECEF to ECI
        % Rotation matrix for simple transformation
        Rz = @(Theta) [cosd(Theta), -sind(Theta), 0; sind(Theta), cosd(Theta), 0; 0, 0, 1];
        
        % Compute the station's ECI coordinates (apply rotation to ECEF)
        statPosECI = Rz(theta) * stationPosECEF;
        
        statVelECi = Rz(theta) * stationVelECEF;
        
        stationECI{i, t} = [statPosECI; statVelECi];
         
    end
    
end

stationECI = stationECI';

end