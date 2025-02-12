classdef OrbConst
    properties (Constant)
        % Define your constants here
        % SpeedOfLight = 299792458; % m/s
        muEarth     = 398600.4415; % km^3/s^2
        muSun       = 1.32712428*10^11; 
        EarthRadius = 6378.13; %km
        J2          = 0.00108248; 
        J3          = -0.0000025323; 
        
        DragRho0 = 0.3614; % kg/km^3 
        r0Drag   = 7000 + EarthRadius; % km
        DragH    = 88.667; % km
    end
end