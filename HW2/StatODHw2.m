%% StatOD Homework 2
% setup add to path 
addpath('C:\Users\user\OneDrive - UCB-O365\Desktop\Spring 2025\StatOD\repo');

%% Setting up Problem 

% Set Constants
mu = Const.OrbConst.muEarth;
J2 = Const.OrbConst.J2; 
Re = Const.OrbConst.EarthRadius;

% Given Orbit
SMA   = 10000;  % km 
eccen = 0.001; 
inc   = 40; % deg
RAAN  = 80; % deg
AOP   = 40; % deg
TA0   = 0; % deg

% Convert from Orbital Elemenets to Cartesian to get initial state vector
[r0,v0] = Utility.OrbCart(SMA,eccen,inc,RAAN,AOP,TA0,Const.OrbConst.muEarth);

% Period of orbit
period = 2*pi*sqrt(SMA^3/Const.OrbConst.muEarth);

% propagate for 15 orbits
t = 0:10:15*period;

% initial state vector
Y0 = [r0;v0]; 

% --- Reference Orbit Propagation --- 
odeoptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[T,Y] = ode45(@Dynamics.NumericJ2Prop, t, Y0, odeoptions, Const.OrbConst.muEarth, J2, Re);

% reference trajectory
refPos = Y(:,1:3);
refVel = Y(:,4:6);





%% Setting Up Measurements

% load in measurements from HW1 4a
MeasStruct = load('MeasHW1prob4a.mat');

% extract measurements from struct
Range     = MeasStruct.Meas4aHW1.RangeMeas;
RangeRate = MeasStruct.Meas4aHW1.RangeRateMeas;
ElevAng   = MeasStruct.Meas4aHW1.ElevAng;

% keep the observations without noise for reference
yHistRef.Range     = Range;
yHistRef.RangeRate = RangeRate; 

% add noise to these measurements
sigmaRange     = 0.001; % km
sigmaRangeRate = 0.000001; % km/s

% create sensor measurement noise covariance matrox
R = [sigmaRange^2 0; 0  sigmaRangeRate^2]; 

% seed randn for repeatability 
seed = 1; 
rng(seed);

noisyRange     = Range + sigmaRange * randn(size(Range));
noisyRangeRate = RangeRate + sigmaRangeRate * randn(size(Range));

% measurements put together in struct
yHist.Range     = noisyRange;
yHist.RangeRate = noisyRangeRate;

% Create A priori covariance matrix
varPos = 1^2; 
varVel = 0.001^2; 

P0 = blkdiag(varPos, varPos, varPos, varVel, varVel, varVel);

% --- State Devaition 
% deltaPos = 1 * randn(3,1);
% deltaVel = 0.001 * randn(3,1);

 deltaPos = [0.539596402821098; -0.689119254306932; -0.34029836838684];
 deltaVel = [0.00104246260737242; -0.000846607337481823; -0.000738866839881026];
 
 % initial state deviaiotn vector
 deltax0 = [deltaPos; deltaVel];
 
 
%% Linear Kalman Filter testing

% staion lat and long
stations.stat1.lat = -35.398333;
stations.stat1.long = 148.981944;

stations.stat2.lat = 40.427222;
stations.stat2.long = 355.749444;

stations.stat3.lat = 35.247164;
stations.stat3.long = 243.205;

% initial rotation of ECEF wrt ECI
theta0 = 122; 
Re = Const.OrbConst.EarthRadius;

% calculate the ECI for each station
[stationECI] = Measurements.StationStateECI(stations, theta0, Re, t);

% get just the stations position in ECI
for i = 1:length(t)
    for j = 1:3
        stationPosECI{i,j} = stationECI{i,j}(1:3);
        stationVelECI{i,j} = stationECI{i,j}(4:6);
    end
end

% Get Visibility Mask
[visibilityMask, viewingAngles] = Measurements.VisibilityMask(stationPosECI, refPos', 10, t);

%---- PASS THIS INTO FUNC Set Initial Conditions
r0   = [-3515.49032703351; 8390.71631024339; 4127.62735255368];
v0   = [-4.35767632217815; -3.35657913876455; 3.1118929278699];
Phi0 = eye(6,6);

% Stack all initial conditions into vector
IC = [r0; v0; reshape(Phi0, [36,1])];

% initial perturbation set earlier
pert = deltax0; 

[xhist, measResHist] = Filters.LinearKF(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, t, mu, J2, Re)