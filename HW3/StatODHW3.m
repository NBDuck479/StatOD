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
tOverall = 0:10:15*period;

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
TruePos = Y(:,1:3);
TrueVel = Y(:,4:6);

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
[visibilityMask, viewingAngles, statNumOb, obTime] = Measurements.VisibilityMask(stationPosECI, TruePos', 10, tOverall);

% spacecraft reference state
spacecraftState = [TruePos'; TrueVel'];

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

% import HW1a measurements
% HW1Meas = load('MeasHW1prob4a.mat');

% HW1Range = HW1Meas.MeasStruct.Meas4aHW1.RangeMeas;

% The reference measurements should only be propagated with what matches
% the filter dynamics! (I guess IDK honestly but let's try it)

[T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tOverall, Y0, odeoptions, mu, J2, Re, 0);

% reference states 
% reference trajectory
refPos = TrajNom(:,1:3);
refVel = TrajNom(:,4:6);

spacecraftStateRef = [refPos'; refVel'];

for i = 1:length(tOverall)
    for j = 1:3    
        rangeRef{i,j}        = spacecraftStateRef(1:3,i) - stationPosECI{i,j};
        rangeNormRef(i,j)    = norm(rangeRef{i,j});
        rangeRateRef{i,j}    = dot(rangeRef{i,j}, spacecraftStateRef(4:6,i) - stationVelECI{i,j}) / rangeNormRef(i,j);
    end
end

yHistRef.Range     = rangeNormRef;
yHistRef.RangeRate = cell2mat(rangeRateRef);

%% Add noise to these measurements
sigmaRange     = 0.001; % km
sigmaRangeRate = 0.000001; % km/s

% create sensor measurement noise covariance matrox
R = [sigmaRange^2 0; 0  sigmaRangeRate^2]; 

% seed randn for repeatability 
seed = 1; 
rng(seed);

% Adding Noise!
noisyRange     = rangeNorm + sigmaRange * randn(size(yHistRef.Range));
noisyRangeRate = cell2mat(rangeRate) + sigmaRangeRate * randn(size(yHistRef.RangeRate));

% measurements put together in struct - sensor measurements!!
yHist.Range     = noisyRange;
yHist.RangeRate = noisyRangeRate;
yHist.obTime    = obTime(obTime(2:end) ~=0 )'; 
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
pert        = [0.539596402821098; -0.689119254306932; -0.34029836838684; 0.00104246260737242; -0.000846607337481823; -0.000738866839881026]% ; 0;0;0];
MeasFlag    = 3; % process both range and range rate

% Set up state noise comp
sigmaXdotdot = 10^-2 * 0.001 % 10^-4; 
sigmaYdotdot = 10^-2 * 0.001% 10^-4;
sigmaZdotdot = 10^-2 * 0.001% 10^-4; 

covSNC = blkdiag(sigmaXdotdot^2, sigmaYdotdot^2, sigmaZdotdot^2); 


%% LKF
% LKF Input struct 
LKFinputs.IC                = Y0; %[Y0(1:6); 0;0;0; reshape(eye(9,9), [9^2, 1])]; % Y0; 
LKFinputs.NumStates         = NumStates; 
LKFinputs.pert              = pert; 
LKFinputs.P0                = P0; %[P0 zeros(6,3); zeros(3,6) blkdiag(0.000001^2, 0.000001^2, 0.000001^2)]; %P0; 
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
LKFinputs.Qframe            = 'RSW';
LKFinputs.DMC               = 0; % 1 for DMC and 0 for none
LKFinputs.tau               = [period/30; period/30; period/30]; % time constants for DMC
LKFinputs.sigmaDMC          = [0.0001; 0.0001; 0.0001]; % uncertainty for DMC accel error states


[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus, refTrajStatesHist] = Filters.LinearKF(LKFinputs);
fig = 1; 
[threeSigmaStates, fig] = Utility.plotLKFEstError(covPlus, NumStates, xhist, fig, obTime, tOverall);
Utility.LKFResidPlotter(measDeltaHist, measResHist, fig, tOverall)

% est error RMS
posX_RMS = sqrt(sum(xhist(1,:).^2)) / length(xhist(1,:));
posY_RMS = sqrt(sum(xhist(2,:).^2)) / length(xhist(2,:));
posZ_RMS = sqrt(sum(xhist(3,:).^2)) / length(xhist(3,:));

velX_RMS = sqrt(sum(xhist(4,:).^2)) / length(xhist(4,:));
velY_RMS = sqrt(sum(xhist(5,:).^2)) / length(xhist(5,:));
velZ_RMS = sqrt(sum(xhist(6,:).^2)) / length(xhist(6,:));

% Save the RMS state values for each run 
posVelRMSerror(:,count) = [posX_RMS;posY_RMS;posZ_RMS; velX_RMS;velY_RMS;velZ_RMS];

count = count +1; 

% ---- post fit residuals ----- 
% Sigma = 10^-2 m/ss
RhoSigma(1)    = 0.000334071;
RhoDotSigma(1) = 3.89549e-07;

% sigma = 10^-4
RhoSigma(2)    = 0.000146496 ;
RhoDotSigma(2) = 2.43586e-07 ;

% sigma = 10^-6
RhoSigma(3)    = 0.000151267 ;
RhoDotSigma(3) = 2.46564e-07 ;

% Sigma = 10^-8 
RhoSigma(4)    = 0.000775564 ;
RhoDotSigma(4) = 7.58683e-07 ;

% Sigma = 10^-10
RhoSigma(5)    = 0.000869199 ;
RhoDotSigma(5) = 8.68993e-07 ;

% Sigma = 10^-12
RhoSigma(6)    = 0.000869235 ;
RhoDotSigma(6) = 8.69034e-07 ;

% Sigma = 10^-14
RhoSigma(7)    = 0.000869235 ;
RhoDotSigma(7) = 8.69034e-07 ;

% Plot the RMS resutls with Sigma
figure(fig)
plot(RhoSigma, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('LKF Post-Fit Range Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 

figure(fig)
plot(RhoDotSigma, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('LKF Post-Fit Range Rate Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 

figure(fig)

% Subplot 1 (RMS Position Error in X Direction)
subplot(3,1,1)
plot(posVelRMSerror(1,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('X Position Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 2 (RMS Position Error in Y Direction)
subplot(3,1,2)
plot(posVelRMSerror(2,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('Y Position Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 3 (RMS Position Error in Z Direction)
subplot(3,1,3)
plot(posVelRMSerror(3,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
title('Z Position Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
grid on;

% Set x-tick labels only for the third subplot (Z Direction)
subplot(3,1,3)  % Make sure we are still working with the third subplot
xticks(1:length(posVelRMSerror(3,:)));  % Set the x-tick positions
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', ...
             '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Add xlabel only to the last subplot (common x-label for all subplots)
xlabel('Sigma', 'Interpreter', 'latex', 'FontSize', 12);
sgtitle('LKF Position RMS')

fig = fig + 1;


figure(fig)

% Subplot 1 (RMS velocity Error in X Direction)
subplot(3,1,1)
plot(posVelRMSerror(4,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('X Velocity Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 2 (RMS Position Error in Y Direction)
subplot(3,1,2)
plot(posVelRMSerror(5,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('Y Velocity Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 3 (RMS Position Error in Z Direction)
subplot(3,1,3)
plot(posVelRMSerror(6,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
title('Z Velocity Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
grid on;

% Set x-tick labels only for the third subplot (Z Direction)
subplot(3,1,3)  % Make sure we are still working with the third subplot
xticks(1:length(posVelRMSerror(3,:)));  % Set the x-tick positions
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', ...
             '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Add xlabel only to the last subplot (common x-label for all subplots)
xlabel('Sigma', 'Interpreter', 'latex', 'FontSize', 12);

sgtitle('LKF Velocity RMS')

fig = fig + 1; 




%% EKF 

% Set up state noise comp
sigmaXdotdot = 10^-4 * 0.001 % 10^-4; 
sigmaYdotdot = 10^-4 * 0.001% 10^-4;
sigmaZdotdot = 10^-4 * 0.001% 10^-4; 

covSNC = blkdiag(sigmaXdotdot^2, sigmaYdotdot^2, sigmaZdotdot^2); 

% make an input struct for EKF also 
EKFinputs.X0         = Y0; %[Y0(1:6); 0;0;0; reshape(eye(9,9), [9^2, 1])]; %[Y0(1:6); reshape(eye(6,6), [6^2, 1])];   
EKFinputs.NumStates  = 6; 
EKFinputs.P0         = P0; %[P0 zeros(6,3); zeros(3,6) blkdiag(0.000001^2, 0.000001^2, 0.000001^2)];% P0;  
EKFinputs.R          = R; 
EKFinputs.yHist      = yHist;
EKFinputs.tVec       = tOverall; 
EKFinputs.stationECI = stationECI;
EKFinputs.mu         = mu; 
EKFinputs.J2         = J2; 
EKFinputs.Re         = Re;
EKFinputs.MeasFlag   = MeasFlag; 
EKFinputs.Q          = covSNC; 
EKFinputs.Qframe     = 'RSW';
EKFinputs.DMC        = 0; % 1 for DMC, 0 for no DMC
EKFinputs.tau        = [period/30; period/30; period/30]; % time constants for DMC
EKFinputs.sigmaDMC   = [0.000001; 0.000001; 0.000001]; % uncertainty for DMC accel error states

[XrefHist,PplusHist,measDeltaHist, post_fit_meas] = Filters.ExtendedKF(EKFinputs);

plot(measDeltaHist(1,:), 'o')

[fig, EKFesterr] = Utility.plotEKFEstError(PplusHist, XrefHist, TruePos, TrueVel, fig, tOverall);

[fig] = Utility.EKFResidPlotter(measDeltaHist, post_fit_meas, fig, tOverall)

% est error RMS
posX_RMS = sqrt(sum(EKFesterr(1,:).^2)) / length(EKFesterr(1,:));
posY_RMS = sqrt(sum(EKFesterr(2,:).^2)) / length(EKFesterr(2,:));
posZ_RMS = sqrt(sum(EKFesterr(3,:).^2)) / length(EKFesterr(3,:));

velX_RMS = sqrt(sum(EKFesterr(4,:).^2)) / length(EKFesterr(4,:));
velY_RMS = sqrt(sum(EKFesterr(5,:).^2)) / length(EKFesterr(5,:));
velZ_RMS = sqrt(sum(EKFesterr(6,:).^2)) / length(EKFesterr(6,:));

% Save the RMS state values for each run 
posVelRMSerrorEKF(:,count) = [posX_RMS;posY_RMS;posZ_RMS; velX_RMS;velY_RMS;velZ_RMS];

count = count +1; 


% RMS for each sigma
% Sigma = 10^-2 m/ss
RhoSigmaEKF(1)    = 5.45935e-06 
RhoDotSigmaEKF(1) = 2.72358e-11 

% Sigma = 10^-4 
RhoSigmaEKF(2)    = 5.467e-06    
RhoDotSigmaEKF(2) = 2.49738e-09

% Sigma = 10^-6
RhoSigmaEKF(3) = 1.22208e-05
RhoDotSigmaEKF(3) = 1.48643e-08

% Sigma = 10 ^-8
RhoSigmaEKF(4)    = 0.00035854
RhoDotSigmaEKF(4) = 4.39504e-07

% Sigma = 10^-10
RhoSigmaEKF(5)    = 0.00042944
RhoDotSigmaEKF(5) = 5.77915e-07

% Sigma = 10 ^-12
RhoSigmaEKF(6)    = 0.000429485
RhoDotSigmaEKF(6) = 5.77993e-07 

% Sigma = 10 ^-14
RhoSigmaEKF(7) = 0.000429485
RhoDotSigmaEKF(7) = 5.77993e-07




% Plot the RMS resutls with Sigma
figure(fig)
plot(RhoSigmaEKF, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('EKF Post-Fit Range Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 

figure(fig)
plot(RhoDotSigmaEKF, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('EKF Post-Fit Range Rate Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Range RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 

figure(fig)

% Subplot 1 (RMS Position Error in X Direction)
subplot(3,1,1)
plot(posVelRMSerrorEKF(1,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('X Position Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 2 (RMS Position Error in Y Direction)
subplot(3,1,2)
plot(posVelRMSerrorEKF(2,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('Y Position Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 3 (RMS Position Error in Z Direction)
subplot(3,1,3)
plot(posVelRMSerrorEKF(3,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
title('Z Position Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
grid on;

% Set x-tick labels only for the third subplot (Z Direction)
subplot(3,1,3)  % Make sure we are still working with the third subplot
xticks(1:length(posVelRMSerror(3,:)));  % Set the x-tick positions
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', ...
             '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Add xlabel only to the last subplot (common x-label for all subplots)
xlabel('Sigma', 'Interpreter', 'latex', 'FontSize', 12);
sgtitle('EKF Position RMS')

fig = fig + 1;


figure(fig)

% Subplot 1 (RMS velocity Error in X Direction)
subplot(3,1,1)
plot(posVelRMSerrorEKF(4,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('X Velocity Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 2 (RMS Position Error in Y Direction)
subplot(3,1,2)
plot(posVelRMSerrorEKF(5,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('Y Velocity Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
grid on;

% Subplot 3 (RMS Position Error in Z Direction)
subplot(3,1,3)
plot(posVelRMSerrorEKF(6,:), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
title('Z Velocity Error', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS Error', 'Interpreter', 'latex', 'FontSize', 12);  % Assuming units in meters
grid on;

% Set x-tick labels only for the third subplot (Z Direction)
subplot(3,1,3)  % Make sure we are still working with the third subplot
xticks(1:length(posVelRMSerror(3,:)));  % Set the x-tick positions
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', ...
             '$\sigma = 10^{-10}$', '$\sigma = 10^{-12}$', '$\sigma = 10^{-14}$'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

% Add xlabel only to the last subplot (common x-label for all subplots)
xlabel('Sigma', 'Interpreter', 'latex', 'FontSize', 12);

sgtitle('EKF Velocity RMS')

fig = fig + 1; 



%% DMC Added to Filters ----- 

NumStates   = 9; 
pert        = [0.539596402821098; -0.689119254306932; -0.34029836838684; 0.00104246260737242; -0.000846607337481823; -0.000738866839881026; 0;0;0];
MeasFlag    = 3; % process both range and range rate

% Set up state noise comp
sigmaXdotdot = 0; %10^-2 * 0.001 % 10^-4; 
sigmaYdotdot = 0; %10^-2 * 0.001% 10^-4;
sigmaZdotdot = 0; %10^-2 * 0.001% 10^-4; 

covSNC = blkdiag(sigmaXdotdot^2, sigmaYdotdot^2, sigmaZdotdot^2); 


% --- LKF ----
% LKF Input struct 
LKFinputs.IC                = [Y0(1:6); 0;0;0; reshape(eye(9,9), [9^2, 1])]; % Y0; 
LKFinputs.NumStates         = NumStates; 
LKFinputs.pert              = pert; 
LKFinputs.P0                = [P0 zeros(6,3); zeros(3,6) blkdiag((1e-07)^2, (1e-07)^2, (1e-07)^2)]; %P0; 
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
LKFinputs.Qframe            = 'RSW';
LKFinputs.DMC               = 1; % 1 for DMC and 0 for none
LKFinputs.tau               = [period/10; period/10; period/10]; % time constants for DMC
LKFinputs.sigmaDMC          = [10^-4 * 0.001; 10^-4 * 0.001; 10^-4 * 0.001]; % uncertainty for DMC accel error states


[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus, refTrajStatesHist] = Filters.LinearKF(LKFinputs);
fig = 1; 
[threeSigmaStates, fig] = Utility.plotLKFEstError(covPlus, NumStates, xhist, fig, obTime, tOverall);
Utility.LKFResidPlotter(measDeltaHist, measResHist, fig, tOverall)

% save the RMS for post fit measurement residuals 

% sigma = 10^-2
RhoRMSLkf(1)    = 0.000154153 
RhoDotRMSLkf(1) = 2.59395e-07 

% sigma = 10^-4
RhoRMSLkf(2)    = 0.000151867
RhoDotRMSLkf(2) = 2.62457e-07 

% sigma = 10^-6
RhoRMSLkf(3)    = 0.00391277
RhoDotRMSLkf(3) = 3.7665e-06

% sigma = 10^-8
RhoRMSLkf(4)    = 0.0001561
RhoDotRMSLkf(4) = 2.48532e-07



% Plot the RMS resutls with Sigma
figure(fig)
plot(1:1:4,RhoRMSLkf, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('LKF Post-Fit Range Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$'});
xticks([1,2,3,4])
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 

figure(fig)
plot(RhoDotRMSLkf, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('LKF Post-Fit Range Rate Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$'});
xticks([1,2,3,4])
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 



% ---- EKF ---- 

% Set up state noise comp
sigmaXdotdot = 0%10^-4 * 0.001 % 10^-4; 
sigmaYdotdot = 0%10^-4 * 0.001% 10^-4;
sigmaZdotdot = 0%10^-4 * 0.001% 10^-4; 

covSNC = blkdiag(sigmaXdotdot^2, sigmaYdotdot^2, sigmaZdotdot^2); 

% make an input struct for EKF also 
EKFinputs.X0         = [Y0(1:6); 0;0;0; reshape(eye(9,9), [9^2, 1])]; %[Y0(1:6); reshape(eye(6,6), [6^2, 1])];   
EKFinputs.NumStates  = 9; 
EKFinputs.P0         = [P0 zeros(6,3); zeros(3,6) blkdiag((1e-09)^2, (1e-09)^2, (1e-09)^2)];% P0;  
EKFinputs.R          = R; 
EKFinputs.yHist      = yHist;
EKFinputs.tVec       = tOverall; 
EKFinputs.stationECI = stationECI;
EKFinputs.mu         = mu; 
EKFinputs.J2         = J2; 
EKFinputs.Re         = Re;
EKFinputs.MeasFlag   = MeasFlag; 
EKFinputs.Q          = covSNC; 
EKFinputs.Qframe     = 'ECI';
EKFinputs.DMC        = 1; % 1 for DMC, 0 for no DMC
EKFinputs.tau        = [period/10; period/10; period/10]; % time constants for DMC
EKFinputs.sigmaDMC   = [10^-6 * 0.001; 10^-6 * 0.001; 10^-6 * 0.001]; % uncertainty for DMC accel error states

[XrefHist,PplusHist,measDeltaHist, post_fit_meas] = Filters.ExtendedKF(EKFinputs);

plot(measDeltaHist(1,:), 'o')

[fig, EKFesterr] = Utility.plotEKFEstError(PplusHist, XrefHist, TruePos, TrueVel, fig, tOverall);

[fig] = Utility.EKFResidPlotter(measDeltaHist, post_fit_meas, fig, tOverall)


% save the RMS for post fit measurement residuals 

% sigma = 10^-2
RhoRMSEkf(1)    = 4.96294e-07
RhoDotRMSEkf(1) =  6.16482e-11

% sigma = 10^-4
RhoRMSEkf(2)   = 0.000462683
RhoDotRMSEkf(2) = 1.57715e-07

% sigma = 10^-6
RhoRMSEkf(3)    = 1.2779e-05
RhoDotRMSEkf(3) = 1.08528e-08

% sigma = 10^-8
RhoRMSEkf(4)    = 1.77333e-05
RhoDotRMSEkf(4) = 1.88103e-08

%sigma = 10^-10
RhoRMSEkf(5)   = 0.000662839
RhoDotRMSEkf(5) = 7.47025e-07


% Plot the RMS resutls with Sigma
figure(fig)
plot(RhoRMSEkf, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('EKF Post-Fit Range Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', '$\sigma = 10^{-10}$'});
xticks([1,2,3,4, 5])
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 

figure(fig)
plot(RhoDotRMSEkf, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('EKF Post-Fit Range Rate Residuals vs $\sigma$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('RMS', 'Interpreter', 'latex', 'FontSize', 12);
xticklabels({'$\sigma = 10^{-2}$', '$\sigma = 10^{-4}$', '$\sigma = 10^{-6}$', '$\sigma = 10^{-8}$', '$\sigma = 10^{-10}$'});
xticks([1,2,3,4, 5])
set(gca, 'TickLabelInterpreter', 'latex');
fig = fig + 1; 
