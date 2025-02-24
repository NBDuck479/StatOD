%% Homework 3 - StatOD
addpath('C:\Users\user\OneDrive - UCB-O365\Desktop\Spring 2025\StatOD\repo')
%% Create new truth data using J2 and J3 dynamics
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
tOverall = 10:10:15*period;

% initial state vector
Y0 = [r0;v0; reshape(eye(6), [6*6,1])]; 

% Set constant values
mu = Const.OrbConst.muEarth;
J2 = Const.OrbConst.J2;
J3 = Const.OrbConst.J3;
Re = Const.OrbConst.EarthRadius;

% propagate the initial state vector
odeoptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[T,Y] = ode45(@Dynamics.DynamicsA_J2_J3, tOverall, Y0, odeoptions, mu, J2 , J3, Re);

% reference trajectory
refPos = Y(:,1:3);
refVel = Y(:,4:6);

%% propagate the station positions
% staion lat and long
stations.stat1.lat = -35.398333;
stations.stat1.long = 148.981944;

stations.stat2.lat = 40.427222;
stations.stat2.long = 355.749444;

stations.stat3.lat = 35.247164;
stations.stat3.long = 243.205;

% initial rotation of ECEF wrt ECI
theta0 = 122; 

% calculate the ECI for each station
[stationECI] = Measurements.StationStateECI(stations, theta0, Re, tOverall);

% get just the stations position in ECI
for i = 1:length(tOverall)
    for j = 1:3
        stationPosECI{i,j} = stationECI{i,j}(1:3);
        stationVelECI{i,j} = stationECI{i,j}(4:6);
    end
end
%% this gives visibility mask 
[visibilityMask, viewingAngles, statNumOb, obTime] = Measurements.VisibilityMask(stationPosECI, refPos', 10, tOverall);

% spacecraft reference state
spacecraftState = [refPos'; refVel'];

%% Measurements

% for i = 1:length(tOverall)
%     for j = 1:3
%         Htilde{i,j} = Measurements.HtildeSC(spacecraftState(:,i), stationECI{i,j}, 3) * visibilityMask(i,j);
%     end
% end

% get measurements for the spacecraft from each station

for i = 1:length(tOverall)
    for j = 1:3    
        range{i,j}        = spacecraftState(1:3,i) - stationPosECI{i,j};
        rangeNorm(i,j)    = norm(range{i,j});
        rangeRate{i,j}    = dot(range{i,j}, spacecraftState(4:6,i) - stationVelECI{i,j}) / rangeNorm(i,j);
    end
end

% find when the obs occur 
% obTimeInd = find(obTime > 0);

% detemrine which station made the ob
% stationObNum = stationObNum(stationObNum ~= 0);

% keep the observations without noise for reference!!
yHistRef.Range     = rangeNorm;
yHistRef.RangeRate = cell2mat(rangeRate); 

%% Add noise to these measurements
sigmaRange     = 0.001; % km
sigmaRangeRate = 0.000001; % km/s

% create sensor measurement noise covariance matrox
R = [sigmaRange^2 0; 0  sigmaRangeRate^2]; 

% seed randn for repeatability 
seed = 1; 
rng(seed);

% Adding Noise!
noisyRange     = yHistRef.Range + sigmaRange * randn(size(yHistRef.Range));
noisyRangeRate = yHistRef.RangeRate + sigmaRangeRate * randn(size(yHistRef.RangeRate));

% measurements put together in struct - sensor measurements!!
yHist.Range     = noisyRange;
yHist.RangeRate = noisyRangeRate;
yHist.obTime    = obTime(obTime ~=0 )'; 
yHist.statNo    = statNumOb'; 

% I also want the time of each of my observations!
stat1vis = find(visibilityMask(:,1) == 1);
stat2vis = find(visibilityMask(:,2) == 1);
stat3vis = find(visibilityMask(:,3) == 1);

% No overlapping measurements! 
intersect(stat1vis, stat2vis);
intersect(stat1vis, stat3vis);
intersect(stat2vis, stat3vis);

rmmissing(visibilityMask)

% Create A priori covariance matrix
varPos = 1^2; 
varVel = .001^2; 

% build initial estimation error covarinace matrix
P0 = blkdiag(varPos, varPos, varPos, varVel, varVel, varVel);


%% Set up parameters for filters
NumStates   = 6; 
pert        = zeros(6,1);
MeasFlag    = 3; % process both range and range rate

% Set up state noise comp
sigmaXdotdot = 10^-5; 
sigmaYdotdot = 10^-5;
sigmaZdotdot = 10^-5; 

covSNC = blkdiag(sigmaXdotdot^2, sigmaYdotdot^2, sigmaZdotdot^2); 


%% LKF
% LKF Input struct 
LKFinputs.IC                = Y0; 
LKFinputs.NumStates         = NumStates; 
LKFinputs.pert              = pert; 
LKFinputs.P0                = P0; 
LKFinputs.R                 = R; 
LKFinputs.mu                = mu; 
LKFinputs.J2                = J2;
LKFinputs.yHist             = yHist; 
LKFinputs.yHistRef          = yHistRef;
LKFinputs.stationECI        = stationECI;
LKFinputs.visibilityMask    = visibilityMask;
LKFinputs.tOverall          = tOverall; 
LKFinputs.Re                = Re; 
% LKFinputs.omegaEarth        = omegaEarth;
% LKFinputs.Area = Area; LKFinputs.Mass = Mass; LKFinputs.DragH = DragH; LKFinputs.r0Drag =  r0Drag; LKFinputs.DragRho0 = DragRho0;
%    ----- LKFinputs.obsHist           = obsHist; % whole struct for observation histories
LKFinputs.MeasFlag          = MeasFlag; % 1 for only range, 2 for range rate, 3 for all
LKFinputs.Q                 = covSNC; 


[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = Filters.LinearKF(LKFinputs);
fig = 1; 
Utility.plotLKFEstError(covPlus, NumStates, xhist, fig, obTime, tOverall)
Utility.LKFResidPlotter(measDeltaHist, measResHist, fig, tOverall)



%% EKF 

% make an input struct for EKF also 
EKFinputs.X0         = Y0; 
EKFinputs.P0         = P0;
EKFinputs.R          = R; 
EKFinputs.yHist      = yHist;
EKFinputs.tVec       = tOverall; 
EKFinputs.stationECI = stationECI;
EKFinputs.mu         = mu; 
EKFinputs.J2         = J2; 
EKFinputs.Re         = Re;
EKFinputs.MeasFlag   = MeasFlag; 
EKFinputs.Q          = covSNC; 

[XrefHist,PplusHist,measDeltaHist] = Filters.ExtendedKF(EKFinputs);

plot(measDeltaHist(1,:), 'o')

Utility.plotEKFEstError(PplusHist, XrefHist, refPos, refVel, fig, tOverall)
