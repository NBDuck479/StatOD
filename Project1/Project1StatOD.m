%% StatOD project 1 
%% define the constants 
Re       = 6378.13;
DragRho0 = 3.6140e-22; % kg/km^3 
r0Drag   = 700 + Re; % km
DragH    = 88.667;  % km
Area     = 3*10^-6; % km^2
Mass     = 970;     % kg
Cd       = 2.0; 
mu       = 398600.4415;
J2       = 0.00108248;
omegaEarth = 7.2921158553*10^-5; % rad
fig = 1; 

% import the given observations
obsHistGiven = load('projectGivenObservations.txt');

% get the observation time
obsHist.time = obsHistGiven(:,1);

% get what station made ob
obsHist.statNo = obsHistGiven(:,2);
% Loop through each element of statNo to check and set the station number
for i = 1:length(obsHist.statNo)
    if obsHist.statNo(i) == 101
        obsHist.statNo(i) = 1;  % Set station number 101 to 1
    elseif obsHist.statNo(i) == 337
        obsHist.statNo(i) = 2;  % Set station number 337 to 2
    elseif obsHist.statNo(i) == 394
        obsHist.statNo(i) = 3;  % Set station number 394 to 3
    end
end

% get range observation
obsHist.Range = obsHistGiven(:,3) * 0.001; % get into km

% get range rate observations
obshist.RangeRate = obsHistGiven(:,4) * 0.001; % get into km/s

% structure observation history appriopriately
yHist.Range = obsHist.Range;
yHist.RangeRate = obshist.RangeRate;

% initial estimation error covariance
mm2km2 = .001 * 0.001;

P0 = diag([1*10^6 * mm2km2, 1*10^6 *mm2km2, 1*10^6 *mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1^10^20 * (0.001*0.001*0.001)^2,...
    1*10^6, 1*10^6,...
    1*10^-10 * mm2km2, 1*10^-10 * mm2km2, 1*10^-10 * mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2]);

% Orbit Intial Conditions
Y0 = [757.7; 5222.607; 4851.5; 2.21321; 4.67834; -5.37130]; 

%% propagate the station positions

% station 101
stations.stat1pos = [-5127.510; -3794.160; 0];

% station 337
stations.stat2pos = [3860.910; 3238.490; 3898.094];

% station 394
stations.stat3pos = [549.505; -1380.872; 6182.197];

% --- overall time
tOverall = 0:18339;

% get station position
[stationECI] = Measurements.StationStateECIproj1(stations, omegaEarth, Re, obsHist.time); % obsHist.time

% get just the stations position in ECI
for i = 1:length(obsHist.time) % obsHist.time
    for j = 1:3
        stationPosECI{i,j} = stationECI{i,j}(1:3);
        stationVelECI{i,j} = stationECI{i,j}(4:6);
    end
end

stationPos = cell2mat(stationPosECI);
stationVel = cell2mat(stationVelECI);


%% Derive the equations of motion for satellite under J2 and Drag 

% provided intial state
Y0(1:6) = [757.700; 5222.607; 4851.500; 2.21321; 4.67834; -5.37130];

% initial mu
Y0(7) = 3.986004415*10^14 * 0.001 * 0.001 *0.001;
Y0(8) = 1.082626925638815*10^-3; 
Y0(9) = 2; 

% initial station states
Y0(10:12) = stations.stat1pos;
Y0(13:15) = stations.stat2pos;
Y0(16:18) = stations.stat3pos;

% initial STM 
Y0(19:342) = reshape(eye(18,18), [324,1]);

% integrate state
odeoptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[T,Y] = ode45(@Dynamics.Numeric_J2_Drag_Prop, obsHist.time, Y0, odeoptions, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0);

% spacecraft reference position
scRefpos = Y(:,1:3);
scRefvel = Y(:,4:6);

stat1posRef = Y(:,10:12);
stat2posRef = Y(:,13:15);
plot(stat2posRef(:,1))


%% Simulate reference Measurements for spacecraft on trajectory

% determine visibillity for each station
[visibilityMask, viewingAngles] = Measurements.VisibilityMask(stationPosECI, scRefpos', -5, obsHist.time);

% manually set visibilty mask because IDk what they are for each GS
visibilityMask(242:244,3) = [NaN; NaN; NaN];
visibilityMask(240:241,2) = [NaN; NaN];
visibilityMask(359:362,3) = [NaN; NaN; NaN; NaN];
visibilityMask(353:358,2) = [NaN; NaN; NaN; NaN; NaN; NaN];
% computed measurements for what is expected
for i = 1:length(obsHist.time)
    for j = 1:3
        range{i,j}        = scRefpos(i,:)' - stationPosECI{i,j} * visibilityMask(i,j);
        rangeNorm(i,j)    = norm(range{i,j});
        rangeRate{i,j}    = dot(range{i,j}, scRefvel(i,:)' - stationVelECI{i,j}) / rangeNorm(i,j);
    end
end

% save off these reference measurements for later use
yHistRef.Range =  rangeNorm;
yHistRef.RangeRate = cell2mat(rangeRate);
%% Set up the Batch Filter
% based on info in project 
sigma_rhoCM = 1; % cm
sigma_rhoKM = sigma_rhoCM * 0.01 * 0.001; % km

sigma_rhoDotmm = 1;
sigma_rhoDotkm = sigma_rhoDotmm * 10^-6; % km/s

% change from st dev to variance
rhoVar = sqrt(sigma_rhoKM);
rhoDotVar = sqrt(sigma_rhoDotkm);

% construct R based on information given
R = [rhoVar 0; 0 rhoDotVar];

% create perturbation
pert = zeros(18,1);

% set convergence tol
critConv = 10^-3; 

% Total number in state vector
NumStates = 18;

% flag for range (1), range rate (2), all (3)
MeasFlag = 3; 

[xhat0est, P0est, phiHist, resid_pfHist, Htilde, preFit_res] = Filters.BatchFilter(Y0, NumStates, pert, P0, R, critConv, yHist, yHistRef, stationECI, visibilityMask, obsHist.time, Y0(7), Y0(8), Re, omegaEarth, Area, Mass, DragRho0, r0Drag, DragH, obsHist.time, obsHist.statNo, MeasFlag);

% --- propagate the estimated Xhat0 in time
% propagate with ode45 as IC + xhat0est
[~,YBatchX0est] = ode45(@Dynamics.Numeric_J2_Drag_Prop, tOverall, [Y0(1:18) + xhat0est; reshape(eye(18,18), [18^2,1])], odeoptions, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0);
[~,YNominal] = ode45(@Dynamics.Numeric_J2_Drag_Prop, tOverall, [Y0(1:18); reshape(eye(18,18), [18^2,1])], odeoptions, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0);

% subtract out nominal state 
xhat0Hist = YBatchX0est - YNominal;

% use STM from ode45 pert prop to get the cov propagation 
batchSTM = YBatchX0est(:,19:end);

for i = 1:length(batchSTM)
   
    P0histBatch{i} = reshape(batchSTM(i,:), [18,18]) * P0est * reshape(batchSTM(i,:), [18,18])'; 
end

Utility.BatchEstErrPlotter(P0histBatch, 18, xhat0Hist, fig, obsHist.time , tOverall)

% wanting to plot residuals on overall time scale
Utility.BatchResidPlotter(preFit_res, resid_pfHist, obsHist, fig);

finalposCov = P0histBatch{end}(1:3,1:3);

% covariance error ellipse at final time
fig = Utility.CovUncertEllipse(finalposCov, scRefpos, xhat0Hist, fig);

[Tbatch, TrajNom] = ode45(@Dynamics.Numeric_J2_Drag_Prop, tOverall, [pert; reshape(eye(18,18), [18^2, 1])], odeoptions, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0);

%% --- get the state deviation throughout time ---
for k = 1:length(tOverall)
    xhatEst(:,k) = reshape(TrajNom(k,18+1:end), [18,18]) * xhat0est;
    
end



%% Set up linearized Kalman Filter

% LKF Input struct 
LKFinputs.IC                = Y0; 
LKFinputs.NumStates         = NumStates; 
LKFinputs.pert              = pert; 
LKFinputs.P0                = P0; 
LKFinputs.R                 = R; 
LKFinputs.yHist             = yHist; 
LKFinputs.yHistRef          = yHistRef;
LKFinputs.stationECI        = stationECI;
LKFinputs.visibilityMask    = visibilityMask;
LKFinputs.tOverall          = tOverall; 
LKFinputs.Re                = Re; 
LKFinputs.omegaEarth        = omegaEarth;
LKFinputs.Area = Area; LKFinputs.Mass = Mass; LKFinputs.DragH = DragH; LKFinputs.r0Drag =  r0Drag; LKFinputs.DragRho0 = DragRho0;
LKFinputs.obsHist           = obsHist; % whole struct for observation histories
LKFinputs.MeasFlag          = MeasFlag; % 1 for only range, 2 for range rate, 3 for all


[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = Filters.LinearKF(LKFinputs);
fig = 1; 
[threeSigmaStates] = Utility.plotLKFEstError(covPlus, NumStates, xhist, fig, obsHist.time, tOverall);

measDeltaHistMat = cell2mat(measDeltaHist);

Utility.LKFhistogramPlotter(measDeltaHistMat)

% plot lkf residuals
Utility.LKFResidPlotter(measDeltaHist, measResHist, obsHist, fig, tOverall)

RMS_Range = sqrt(sum(rmmissing(measResHist(1,:)).^2))/length(rmmissing(measResHist(1,:)));
RMS_RangeRate = sqrt(sum(rmmissing(measResHist(2,:)).^2))/length(rmmissing(measResHist(2,:)));
        
measDeltaHist = cell2mat(measDeltaHist);
RMS_Range_pf = sqrt(sum(rmmissing(measDeltaHist(1,:)).^2)) / length(rmmissing(measDeltaHist(1,:)));
RMS_Rangerate_pf = sqrt(sum(rmmissing(measDeltaHist(2,:)).^2)) / length(rmmissing(measDeltaHist(2,:)));

%% --- final covaraince ellipse plotting
% Pull the final covaraince matrix
ckfFinalCov = covPlus{end};

% extract position final covaraince
finalposCov = ckfFinalCov(1:3,1:3);

fig = Utility.CovUncertEllipse(finalposCov, scRefpos, xhist, fig);

%% Relative trength of measurements

% only prcess the range measurements in batch
measFlag = 1; % 1 is range, 2 is range rate, 3 is both

[xhat0est, P0est, phiHist, resid_pfHist, Htilde, preFit_res] = Filters.BatchFilter(Y0, NumStates, pert, P0, R, critConv, yHist, yHistRef, stationECI, visibilityMask, obsHist.time, Y0(7), Y0(8), Re, omegaEarth, Area, Mass, DragRho0, r0Drag, DragH, obsHist.time, obsHist.statNo, measFlag);

finalposCov = P0histBatch{end}(1:3,1:3);

% Plot the estimation error and final covariance ellipse
fig = Utility.CovUncertEllipse(finalposCov, scRefpos, xhat0Hist, fig);

% estimation error
Utility.BatchEstErrPlotter(P0histBatch, 18, xhat0Hist, fig, obsHist.time , tOverall)


measFlag = 2; % 1 is range, 2 is range rate, 3 is both

[xhat0est, P0est, phiHist, resid_pfHist, Htilde, preFit_res] = Filters.BatchFilter(Y0, NumStates, pert, P0, R, critConv, yHist, yHistRef, stationECI, visibilityMask, obsHist.time, Y0(7), Y0(8), Re, omegaEarth, Area, Mass, DragRho0, r0Drag, DragH, obsHist.time, obsHist.statNo, measFlag);

finalposCov = P0histBatch{end}(1:3,1:3);

% Plot the estimation error and final covariance ellipse
fig = Utility.CovUncertEllipse(finalposCov, scRefpos, xhat0Hist, fig);

% estimation error
Utility.BatchEstErrPlotter(P0histBatch, 18, xhat0Hist, fig, obsHist.time , tOverall);

%% Changing fixed stations

% -- No stations fixed
P0 = diag([1*10^6 * mm2km2, 1*10^6 *mm2km2, 1*10^6 *mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1^10^20 * (0.001*0.001*0.001)^2,...
    1*10^6, 1*10^6,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2]);

[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = Filters.LinearKF(Y0, NumStates, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, tOverall, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0, obsHist.time, obsHist.statNo, MeasFlag);

[threeSigmaStates] = Utility.plotLKFEstError(covPlus, NumStates, xhist, fig, obsHist.time, tOverall);


% -- Different station fixed 
P0 = diag([1*10^6 * mm2km2, 1*10^6 *mm2km2, 1*10^6 *mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1^10^20 * (0.001*0.001*0.001)^2,...
    1*10^6, 1*10^6,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1*10^6 * mm2km2, 1*10^6 * mm2km2, 1*10^6 * mm2km2,...
    1*10^-10 * mm2km2, 1*10^-10 * mm2km2, 1*10^-10 * mm2km2]);

[xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = Filters.LinearKF(Y0, NumStates, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, tOverall, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0, obsHist.time, obsHist.statNo, MeasFlag);

[threeSigmaStates] = Utility.plotLKFEstError(covPlus, NumStates, xhist, fig, obsHist.time, tOverall);


%% Extra Elements

% make a priori covaraince more realistic
P0 = diag([1 * mm2km2, 1 *mm2km2, 1 *mm2km2,...
    0.5 * mm2km2, 0.5 * mm2km2, 0.5 * mm2km2,...
    1 * (0.001*0.001*0.001)^2,...
    1*10^-6, 1*10^-6,...
    1 * mm2km2, 1 * mm2km2, 1* mm2km2,...
    1 * mm2km2, 1 * mm2km2, 1 * mm2km2,...
    1 * mm2km2, 1 * mm2km2, 1 * mm2km2]);


measFlag = 3; % 1 is range, 2 is range rate, 3 is both

[xhat0est, P0est, phiHist, resid_pfHist, Htilde, preFit_res] = Filters.BatchFilter(Y0, NumStates, pert, P0, R, critConv, yHist, yHistRef, stationECI, visibilityMask, obsHist.time, Y0(7), Y0(8), Re, omegaEarth, Area, Mass, DragRho0, r0Drag, DragH, obsHist.time, obsHist.statNo, measFlag);

finalposCov = P0histBatch{end}(1:3,1:3);

% estimation error
Utility.BatchEstErrPlotter(P0histBatch, 18, xhat0Hist, fig, obsHist.time , tOverall)