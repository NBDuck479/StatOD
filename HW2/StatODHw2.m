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
% Integrate Trajectory
[T, Y] = ode45(@Dynamics.NumericJ2Prop, t, [Y0; reshape(eye(6,6), [36,1])], odeoptions, mu, J2, Re);

% reference trajectory
refPos = Y(:,1:3);
refVel = Y(:,4:6);





%% Setting Up Measurements

% load in measurements from HW1 4a
MeasStruct = load('MeasHW1prob4a.mat');

% extract measurements from struct
Range     = MeasStruct.MeasStruct.Meas4aHW1.RangeMeas;
RangeRate = MeasStruct.MeasStruct.Meas4aHW1.RangeRateMeas;
ElevAng   = MeasStruct.MeasStruct.Meas4aHW1.ElevAng;

% keep the observations without noise for reference!!
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

% NO NOISE IN NOISY MEASUREMENTS FOR DEBUGGING!
noisyRange     = Range + sigmaRange * randn(size(Range));
noisyRangeRate = RangeRate + sigmaRangeRate * randn(size(Range));

% measurements put together in struct - sensor measurements!!
yHist.Range     = noisyRange;
yHist.RangeRate = noisyRangeRate;

% Create A priori covariance matrix
varPos = 1^2; 
varVel = .001^2; 

% build initial estimation error covarinace matrix
P0 = blkdiag(varPos, varPos, varPos, varVel, varVel, varVel);

% --- State Devaition 
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

% plot station position to check
statPos = cell2mat(stationPosECI)';
plot(statPos(:,3:3:end)')

statVel = cell2mat(stationVelECI)';
plot(statVel(:,3:3:end)')

% Get Visibility Mask
[visibilityMask, viewingAngles] = Measurements.VisibilityMask(stationPosECI, refPos', 10, t);

%---- PASS THIS INTO FUNC Set Initial Conditions
r0   = refPos(1,:);
v0   = refVel(1,:);
Phi0 = eye(6,6);

% Stack all initial conditions into vector
IC = [r0'; v0'; reshape(Phi0, [36,1])];

% initial perturbation set earlier
pert = deltax0; 
%% Filter testing
% No noisy measurements!
% testing with NO A PRIOIR STATE ERROR

[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = Filters.LinearKF(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, t, mu, J2, Re);

covPlusmat = cell2mat(covPlus);
% get each state uncertainty
stateUncert1 = covPlusmat(1, 1:6:end);
stateUncert2 = covPlusmat(2, 2:6:end);
stateUncert3 = covPlusmat(3, 3:6:end);
stateUncert4 = covPlusmat(4, 4:6:end);
stateUncert5 = covPlusmat(5, 5:6:end);
stateUncert6 = covPlusmat(6, 6:6:end);

threeSigmaState1 = 3*sqrt(stateUncert1);
threeSigmaState2 = 3*sqrt(stateUncert2);
threeSigmaState3 = 3*sqrt(stateUncert3);
threeSigmaState4 = 3*sqrt(stateUncert4);
threeSigmaState5 = 3*sqrt(stateUncert5);
threeSigmaState6 = 3*sqrt(stateUncert6);


% RMS state error
% Initialize a variable to accumulate the squared errors
squared_errors = 0;

rangeMeasRes = measResHist(1,:);
%% Loop over each time step to calculate the squared error
for i = 1:length(xhist)
    % Compute the state error at time t
    if isnan(rangeMeasRes(i))
        % do nothing
    else 
        % add to error
    state_error = rangeMeasRes(i);
    
    % Square the norm of the state error and accumulate it
    squared_errors = squared_errors + state_error.^2;
    end
end

% Compute the RMS state error
rms_state_error = sqrt(squared_errors / length(xhist));

% Display the RMS state error
disp(rms_state_error);

disp(['RMS State Error: ', num2str(rms_state_error)]);


% plot the estimation error
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState1), threeSigmaState1, 'r--', 1:length(threeSigmaState1), -threeSigmaState1, 'r--');

hold on
plot(xhist(1,:)')
grid on
ylabel('X position [km]')

subplot(3,1,2)
plot(1:length(threeSigmaState2), threeSigmaState2, 'r--', 1:length(threeSigmaState2), -threeSigmaState2, 'r--');

hold on
plot(xhist(2,:)')
grid on
ylabel('Y position [km]')

subplot(3,1,3)
plot(1:length(threeSigmaState3), threeSigmaState3, 'r--', 1:length(threeSigmaState3), -threeSigmaState3, 'r--');

hold on
plot(xhist(3,:)')
grid on
ylabel('Z position [km]')

sgtitle('LKF Estimation Error Position')
fig = fig+1;

% velcoity estimation error 

figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState4), threeSigmaState4, 'r--', 1:length(threeSigmaState4), -threeSigmaState4, 'r--');

hold on
plot(xhist(4,:)')
grid on
ylabel('X velocity [km/s]')

subplot(3,1,2)
plot(1:length(threeSigmaState5), threeSigmaState5, 'r--', 1:length(threeSigmaState5), -threeSigmaState5, 'r--');

hold on
plot(xhist(5,:)')
grid on
ylabel('Z velocity [km/s]')

subplot(3,1,3)
plot(1:length(threeSigmaState6), threeSigmaState6, 'r--', 1:length(threeSigmaState6), -threeSigmaState6, 'r--');

hold on
plot(xhist(6,:)')
grid on
ylabel('Z velocity [km/s]')

sgtitle('LKF Estimation Error Velocity')
fig = fig+1;




subplot(3,1,1)
plot(statNumObHist', 'o')
title('statNumOb')

subplot(3,1,2)
plot(cell2mat(measDeltaHist)', 'o')
title('measDeltaHist')

subplot(3,1,3)
plot(cell2mat(estimatedDeviationOb)', 'o');
title('estimatedDeviationOb');
ylim([-1 1])


% plot the post fit residuals 
figure(fig)
subplot(2,1,1)
plot(measResHist(1,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range [km]')

subplot(2,1,2)
plot(measResHist(2,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range Rate [km/s]')

sgtitle('LKF Post-Fit Residuals')
fig = fig+1; 


cell2mat(measDeltaHist)
figure(11); plot(xhist'); title('xhist');ylim([-1 1])
figure(12); plot(measResHist', 'o'); title('measResHist');
figure(13); plot(cell2mat(measDeltaHist)', 'o'); title('measDeltaHist');
figure(14); plot(cell2mat(estimatedDeviationOb)', 'o'); title('estimatedDeviationOb'); ylim([-1 1])


%% EKF Testing

% init total state
totState = [r0';v0']; %+ pert;

X0 = [totState; reshape(eye(6,6), [36,1])];
[XrefHist,PplusHist,measDeltaHist] = Filters.ExtendedKF(X0, P0, R, yHist, t, stationECI, mu, J2, Re);

% plot function for EKF estimation error
Utility.plotEKFEstError(PplusHist, XrefHist, refPos, refVel, fig)


figure(fig)
subplot(2,1,1)
plot(measDeltaHist(1,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range [km]')


subplot(2,1,2)
plot(measDeltaHist(2,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range Rate [km/s]')

sgtitle('EKF Post-Fit Residuals')
fig = fig + 1; 

%% Batch filter testing


[xhat0est, P0est, phiHist, measDeltaHist, Htilde] = Filters.BatchFilter(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, t, mu, J2, Re);

% --- get the post fit residuals ---
for k = 1:length(t)/2-1
    xhatEst(:,k) = reshape(phiHist(k,:), [6,6]) * xhat0est;
    
    % predicted measurement
    predMeas(:,k) = Htilde{k}*xhatEst(:,k);
    
    postFit(:,k) = measDeltaHist(:,k) - predMeas(:,k);
    
    % covariance
    Pprop{k} = reshape(phiHist(k,:), [6,6]) * P0est * reshape(phiHist(k,:), [6,6])';
end



% 2 sigma covariance plots
PpropHistMat = abs(cell2mat(Pprop));

% get each state uncertainty
stateUncert1 = PpropHistMat(1, 1:6:end);
stateUncert2 = PpropHistMat(2, 2:6:end);
stateUncert3 = PpropHistMat(3, 3:6:end);
stateUncert4 = PpropHistMat(4, 4:6:end);
stateUncert5 = PpropHistMat(5, 5:6:end);
stateUncert6 = PpropHistMat(6, 6:6:end);

threeSigmaState1 = 3*sqrt(stateUncert1);
threeSigmaState2 = 3*sqrt(stateUncert2);
threeSigmaState3 = 3*sqrt(stateUncert3);
threeSigmaState4 = 3*sqrt(stateUncert4);
threeSigmaState5 = 3*sqrt(stateUncert5);
threeSigmaState6 = 3*sqrt(stateUncert6);


% position error for Batch
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState1), threeSigmaState1, 'r--', 1:length(threeSigmaState1), -threeSigmaState1, 'r--');

hold on
plot(xhatEst(1,:)')
grid on
ylabel('X position [km]')

subplot(3,1,2)
plot(1:length(threeSigmaState2), threeSigmaState2, 'r--', 1:length(threeSigmaState2), -threeSigmaState2, 'r--');

hold on
plot(xhatEst(2,:)')
grid on
ylabel('Y position [km]')

subplot(3,1,3)
plot(1:length(threeSigmaState3), threeSigmaState3, 'r--', 1:length(threeSigmaState3), -threeSigmaState3, 'r--');

hold on
plot(xhatEst(3,:)')
grid on
ylabel('Z position [km]')

sgtitle('Batch Estimation Error Position')
fig = fig+1;


% Velocity error for batch
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState4), threeSigmaState4, 'r--', 1:length(threeSigmaState4), -threeSigmaState4, 'r--');
hold on
plot(xhatEst(4,:)')

grid on
ylabel('X Velocity [km/s]')

subplot(3,1,2)
plot(1:length(threeSigmaState5), threeSigmaState5, 'r--', 1:length(threeSigmaState5), -threeSigmaState5, 'r--');
hold on
plot(xhatEst(5,:)')

grid on
ylabel('Y Velocity [km/s]')

subplot(3,1,3)
plot(1:length(threeSigmaState6), threeSigmaState6, 'r--', 1:length(threeSigmaState6), -threeSigmaState6, 'r--');
hold on
plot(xhatEst(6,:)')

grid on
ylabel('Z Velocity [km/s]')

sgtitle('Batch Estimation Error Velocity')
fig = fig+1;

plot(xhatEst')

plot(predMeas(1,:)', 'o')
plot(predMeas(2,:)', 'o')

figure(fig)
subplot(2,1,1)
plot(postFit(1,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range [km]')

subplot(2,1,2)
plot(postFit(2,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range Rate [km/s]')

sgtitle('Batch Filter Post-Fit Residuals')

% get the 3 sigma for each state form propagated covariance
Ppropmat = cell2mat(Pprop);

state1uncert = Ppropmat(1, 1:6:end);

ThreeSigmaState1 = 3*sqrt(abs(state1uncert));





%% Create new data with J3
% Set Constants
mu = Const.OrbConst.muEarth;
J2 = Const.OrbConst.J2; 
Re = Const.OrbConst.EarthRadius;
J3 = Const.OrbConst.J3;

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
% Integrate Trajectory
[T, Y_J3] = ode45(@Dynamics.DynamicsA_J2_J3, t, [Y0; reshape(eye(6,6), [36,1])], odeoptions, mu, J2, J3, Re);

% reference trajectory
refPosJ3 = Y_J3(:,1:3);
refVelJ3 = Y_J3(:,4:6);

% plot difference between J2 and J3
plot(refPosJ3 - refPos)


%% create measurements with J2 + J3
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


%% this gives visibility mask 
[visibilityMask, viewingAngles] = Measurements.VisibilityMask(stationPosECI, refPosJ3', 10, t);

% spacecraft reference state
spacecraftState = [refPosJ3'; refVelJ3'];

%% linearized H matrix from reference state
for i = 1:length(t)
    for j = 1:3
        Htilde{i,j} = Measurements.HtildeSC(spacecraftState(:,i), stationECI{i,j}) * visibilityMask(i,j);
    end
end

% get measurements for the spacecraft from each station

for i = 1:length(t)
    for j = 1:3
        range{i,j}        = spacecraftState(1:3,i) - stationPosECI{i,j} * visibilityMask(i,j);
        rangeNormJ3(i,j)    = norm(range{i,j});
        rangeRateJ3{i,j}    = dot(range{i,j}, spacecraftState(4:6,i) - stationVelECI{i,j}) / rangeNormJ3(i,j);
    end
end

% test the J3 measurements
figure(fig)
plot(rangeNormJ3)
fig = fig + 1; 

figure(fig)
plot(cell2mat(rangeRateJ3))

fig = fig +1; 


% load in J2 measurements!
MeasStruct = load('MeasHW1prob4a.mat');

% extract measurements from struct
RangeJ2     = MeasStruct.MeasStruct.Meas4aHW1.RangeMeas;
RangeRateJ2 = MeasStruct.MeasStruct.Meas4aHW1.RangeRateMeas;
ElevAngJ2   = MeasStruct.MeasStruct.Meas4aHW1.ElevAng;


% plot differences
figure(fig)
plot(rangeNormJ3 - RangeJ2)
grid on
xlabel('Time [sec]')
ylabel('Delta Range [km]')
legend('station 1', 'station 2', 'station 3')
title('Difference in Range Between J2 and J3')
fig = fig +1; 

figure(fig)
plot(cell2mat(rangeRateJ3) - RangeRateJ2)
xlabel('Time [sec]')
grid on
ylabel('Delta Range Rate [km/s]')
legend('station 1', 'station 2', 'station 3')
title('Difference in Range Rate Between J2 and J3')
fig = fig +1; 

% plot difference between J2 and J3
figure(fig)
plot(refPosJ3 - refPos)
grid on
ylabel('Spacecraft Position [km]')
xlabel('Time')
title('Difference in spacecraft Position J2 and J3')
fig = fig +1; 


figure(fig)
plot(refVelJ3 - refVel)
grid on
ylabel('Spacecraft Velocity [km/s]')
xlabel('Time')
title('Difference in spacecraft Velocity J2 and J3')
fig = fig +1; 

%% Use J3 data for filter to process but J2 dynamics 
% NO NOISE IN NOISY MEASUREMENTS FOR DEBUGGING!
sigmaRange     = 0.001; % km
sigmaRangeRate = 0.000001; % km/s

noisyRangeJ3     = rangeNormJ3 + sigmaRange * randn(size(rangeNormJ3));
noisyRangeRateJ3 = cell2mat(rangeRateJ3) + sigmaRangeRate * randn(size(cell2mat(rangeRateJ3)));

% measurements put together in struct - sensor measurements!!
yHist.Range     = noisyRangeJ3;
yHist.RangeRate = noisyRangeRateJ3;


% initial conditoin
IC = [refPosJ3(1,:)'; refVelJ3(1,:)'; reshape(eye(6,6), [36,1])]

% filters 
[xhat0est, P0est, phiHist, measDeltaHist, Htilde] = Filters.BatchFilter(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, t, mu, J2, Re);

% --- get the post fit residuals ---
for k = 1:length(t)-1
    xhatEst(:,k) = reshape(phiHist(k,:), [6,6]) * xhat0est;
    
    % predicted measurement
    predMeas(:,k) = Htilde{k}*xhatEst(:,k);
    
    postFit(:,k) = measDeltaHist(:,k) - predMeas(:,k);
    
    % covariance
    Pprop{k} = reshape(phiHist(k,:), [6,6]) * P0est * reshape(phiHist(k,:), [6,6])';
end



% 2 sigma covariance plots
PpropHistMat = abs(cell2mat(Pprop));

% get each state uncertainty
stateUncert1 = PpropHistMat(1, 1:6:end);
stateUncert2 = PpropHistMat(2, 2:6:end);
stateUncert3 = PpropHistMat(3, 3:6:end);
stateUncert4 = PpropHistMat(4, 4:6:end);
stateUncert5 = PpropHistMat(5, 5:6:end);
stateUncert6 = PpropHistMat(6, 6:6:end);

threeSigmaState1 = 3*sqrt(stateUncert1);
threeSigmaState2 = 3*sqrt(stateUncert2);
threeSigmaState3 = 3*sqrt(stateUncert3);
threeSigmaState4 = 3*sqrt(stateUncert4);
threeSigmaState5 = 3*sqrt(stateUncert5);
threeSigmaState6 = 3*sqrt(stateUncert6);


% position error for Batch
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState1), threeSigmaState1, 'r--', 1:length(threeSigmaState1), -threeSigmaState1, 'r--');

hold on
plot(xhatEst(1,:)')
grid on
ylabel('X position [km]')

subplot(3,1,2)
plot(1:length(threeSigmaState2), threeSigmaState2, 'r--', 1:length(threeSigmaState2), -threeSigmaState2, 'r--');

hold on
plot(xhatEst(2,:)')
grid on
ylabel('Y position [km]')

subplot(3,1,3)
plot(1:length(threeSigmaState3), threeSigmaState3, 'r--', 1:length(threeSigmaState3), -threeSigmaState3, 'r--');

hold on
plot(xhatEst(3,:)')
grid on
ylabel('Z position [km]')

sgtitle('Batch Estimation Error Position')
fig = fig+1;


% Velocity error for batch
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState4), threeSigmaState4, 'r--', 1:length(threeSigmaState4), -threeSigmaState4, 'r--');
hold on
plot(xhatEst(4,:)')

grid on
ylabel('X Velocity [km/s]')

subplot(3,1,2)
plot(1:length(threeSigmaState5), threeSigmaState5, 'r--', 1:length(threeSigmaState5), -threeSigmaState5, 'r--');
hold on
plot(xhatEst(5,:)')

grid on
ylabel('Y Velocity [km/s]')

subplot(3,1,3)
plot(1:length(threeSigmaState6), threeSigmaState6, 'r--', 1:length(threeSigmaState6), -threeSigmaState6, 'r--');
hold on
plot(xhatEst(6,:)')

grid on
ylabel('Z Velocity [km/s]')

sgtitle('Batch Estimation Error Velocity')
fig = fig+1;

[XrefHist,PplusHist,measDeltaHist] = Filters.ExtendedKF(X0, P0, R, yHist, t, stationECI, mu, J2, Re);

% 2 sigma covariance plots
PplusHistMat = cell2mat(PplusHist);

% get each state uncertainty
stateUncert1 = PplusHistMat(1, 1:6:end);
stateUncert2 = PplusHistMat(2, 2:6:end);
stateUncert3 = PplusHistMat(3, 3:6:end);
stateUncert4 = PplusHistMat(4, 4:6:end);
stateUncert5 = PplusHistMat(5, 5:6:end);
stateUncert6 = PplusHistMat(6, 6:6:end);

threeSigmaState1 = 3*sqrt(stateUncert1);
threeSigmaState2 = 3*sqrt(stateUncert2);
threeSigmaState3 = 3*sqrt(stateUncert3);
threeSigmaState4 = 3*sqrt(stateUncert4);
threeSigmaState5 = 3*sqrt(stateUncert5);
threeSigmaState6 = 3*sqrt(stateUncert6);
% estimation error
EKFesterr = XrefHist - [refPos(1:14928, :), refVel(1:14928,:)]';

% position error for EKF
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState1), threeSigmaState1, 'r--', 1:length(threeSigmaState1), -threeSigmaState1, 'r--');
ylim([-0.03 0.03])
hold on
plot(EKFesterr(1,:)')
grid on
ylabel('X position [km]')

subplot(3,1,2)
plot(1:length(threeSigmaState2), threeSigmaState2, 'r--', 1:length(threeSigmaState2), -threeSigmaState2, 'r--');
ylim([-0.03 0.03])
hold on
plot(EKFesterr(2,:)')
grid on
ylabel('Y position [km]')

subplot(3,1,3)
plot(1:length(threeSigmaState3), threeSigmaState3, 'r--', 1:length(threeSigmaState3), -threeSigmaState3, 'r--');
ylim([-0.03 0.03])
hold on
plot(EKFesterr(3,:)')
grid on
ylabel('Z position [km]')

sgtitle('EKF Estimation Error Position')
fig = fig+1;

% Velocity error for EKF
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState4), threeSigmaState4, 'r--', 1:length(threeSigmaState4), -threeSigmaState4, 'r--');
hold on
plot(EKFesterr(4,:)')
ylim([-0.00001 0.00001])
grid on
ylabel('X Velocity [km/s]')

subplot(3,1,2)
plot(1:length(threeSigmaState5), threeSigmaState5, 'r--', 1:length(threeSigmaState5), -threeSigmaState5, 'r--');
hold on
plot(EKFesterr(5,:)')
ylim([-0.00001 0.00001])
grid on
ylabel('Y Velocity [km/s]')

subplot(3,1,3)
plot(1:length(threeSigmaState6), threeSigmaState6, 'r--', 1:length(threeSigmaState6), -threeSigmaState6, 'r--');
hold on
plot(EKFesterr(6,:)')
ylim([-0.00001 0.00001])
grid on
ylabel('Z Velocity [km/s]')

sgtitle('EKF Estimation Error Velocity')
fig = fig+1;


figure(fig)
subplot(2,1,1)
plot(measDeltaHist(1,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range [km]')


subplot(2,1,2)
plot(measDeltaHist(2,:)', 'o')
grid on
xlabel('seconds')
ylabel('Range Rate [km/s]')

sgtitle('EKF Post-Fit Residuals')
fig = fig + 1; 


[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = Filters.LinearKF(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, t, mu, J2, Re);
covPlusmat = cell2mat(covPlus);
% get each state uncertainty
stateUncert1 = covPlusmat(1, 1:6:end);
stateUncert2 = covPlusmat(2, 2:6:end);
stateUncert3 = covPlusmat(3, 3:6:end);
stateUncert4 = covPlusmat(4, 4:6:end);
stateUncert5 = covPlusmat(5, 5:6:end);
stateUncert6 = covPlusmat(6, 6:6:end);

threeSigmaState1 = 3*sqrt(stateUncert1);
threeSigmaState2 = 3*sqrt(stateUncert2);
threeSigmaState3 = 3*sqrt(stateUncert3);
threeSigmaState4 = 3*sqrt(stateUncert4);
threeSigmaState5 = 3*sqrt(stateUncert5);
threeSigmaState6 = 3*sqrt(stateUncert6);

% plot the estimation error
figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState1), threeSigmaState1, 'r--', 1:length(threeSigmaState1), -threeSigmaState1, 'r--');

hold on
plot(xhist(1,:)')
grid on
ylabel('X position [km]')

subplot(3,1,2)
plot(1:length(threeSigmaState2), threeSigmaState2, 'r--', 1:length(threeSigmaState2), -threeSigmaState2, 'r--');

hold on
plot(xhist(2,:)')
grid on
ylabel('Y position [km]')

subplot(3,1,3)
plot(1:length(threeSigmaState3), threeSigmaState3, 'r--', 1:length(threeSigmaState3), -threeSigmaState3, 'r--');

hold on
plot(xhist(3,:)')
grid on
ylabel('Z position [km]')

sgtitle('LKF Estimation Error Position')
fig = fig+1;

% velcoity estimation error 

figure(fig)
subplot(3,1,1)
plot(1:length(threeSigmaState4), threeSigmaState4, 'r--', 1:length(threeSigmaState4), -threeSigmaState4, 'r--');

hold on
plot(xhist(4,:)')
grid on
ylabel('X velocity [km/s]')

subplot(3,1,2)
plot(1:length(threeSigmaState5), threeSigmaState5, 'r--', 1:length(threeSigmaState5), -threeSigmaState5, 'r--');

hold on
plot(xhist(5,:)')
grid on
ylabel('Z velocity [km/s]')

subplot(3,1,3)
plot(1:length(threeSigmaState6), threeSigmaState6, 'r--', 1:length(threeSigmaState6), -threeSigmaState6, 'r--');

hold on
plot(xhist(6,:)')
grid on
ylabel('Z velocity [km/s]')

sgtitle('LKF Estimation Error Velocity')
fig = fig+1;