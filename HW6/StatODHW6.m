%% StatOD Homework 6 
addpath('C:\Users\user\OneDrive - UCB-O365\Desktop\Spring 2025\StatOD\repo');
mu = Const.OrbConst.muEarth;
Re = Const.OrbConst.EarthRadius; 
J2 = Const.OrbConst.J2;

DragRho0 = 3.6140e-22; % kg/km^3 
r0Drag   = 700 + Re; % km
DragH    = 88.667;  % km
Area     = 3*10^-6; % km^2
Mass     = 970;     % kg

omegaEarth = 7.2921158553*10^-5; % rad

% Load in observations from previous Assignment
% yHistLoad = load('yHist.mat');
% yHist = yHistLoad.yHist; 

 yHistLoadJ3 = load('yHistJ3.mat');
 yHist = yHistLoadJ3.yHist; 

% from previous homework to make sensor observations
sigmaRange     = 0.001; % km
sigmaRangeRate = 0.000001; % km/s

% create sensor measurement noise covariance matrox
R = [sigmaRange^2 0; 0  sigmaRangeRate^2]; 

% --- Dynamically unpack previous homework struct
HomeworkDataStruct = load('HomeworkData.mat');

% Get all field names in the struct
fields = fieldnames(HomeworkDataStruct.HomeworkData);

% Loop through each field and access its value dynamically
for i = 1:length(fields)
    % Get the field name
    fieldName = fields{i};
    
    % Access the value of the field using dynamic field referencing
    fieldValue = HomeworkDataStruct.HomeworkData.(fieldName);
    
    eval([fieldName ' = fieldValue']);
    
end


%% Implement the UKF

% Set up state noise comp
sigmaXdotdot = 10^-2 * 0.001; % 10^-4; 
sigmaYdotdot = 10^-2 * 0.001; % 10^-4;
sigmaZdotdot = 10^-2 * 0.001; % 10^-4; 

covSNC = blkdiag(sigmaXdotdot^2, sigmaYdotdot^2, sigmaZdotdot^2); 

% SRIF input struct
UKFinputs.IC                = Y0; %[Y0(1:6); 0;0;0; reshape(eye(9,9), [9^2, 1])]; % Y0; 
UKFinputs.NumStates         = NumStates; 
UKFinputs.P0                = 10^4*P0; %[P0 zeros(6,3); zeros(3,6) blkdiag(0.000001^2, 0.000001^2, 0.000001^2)]; %P0; 
UKFinputs.X0                = Y0(1:6) + pert;
UKFinputs.mu                = mu;  
UKFinputs.J2                = J2;
UKFinputs.yHist             = yHist; 
UKFinputs.stationECI        = stationECI;
UKFinputs.visibilityMask    = visibilityMask;
UKFinputs.tOverall          = tOverall; 
UKFinputs.Re                = Re; 
UKFinputs.MeasFlag          = 3; % 1 for only range, 2 for range rate, 3 for all
UKFinputs.a                 = 10^-4; 
UKFinputs.B                 = 2; 
UKFinputs.R                 = R; 
UKFinputs.Q                 = covSNC; 
UKFinputs.Qframe            = 'ECI';
% UKFinputs.processNoise      = 1; 

[Xhist, Phist, OminusC] = Filters.UnscentedKF(UKFinputs)

% integrate the nominal trqjectory
odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
% [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tOverall, Y0, odeOptions, mu, J2, Re);
[T, TrajNom] = ode45(@Dynamics.DynamicsA_J2_J3, tOverall, Y0, odeOptions, mu, J2 , Const.OrbConst.J3, Re);

 Utility.plotUKFEstError(Phist, Xhist, TrajNom(2:end,1:6), 1, tOverall, NumStates)

% Estimation Error 
estErr = Xhist(1:3,:) - TrajNom(2:end,1:3)';

OCrange = OminusC(1,:);
OCrangeRate = OminusC(2,:);

OCrange(OCrange == 0) = NaN;
OCrangeRate(OCrangeRate == 0) = NaN; 

% plot residuals 
figure(100)
plot(OCrangeRate, 'o')
grid on
title('UKF O Minus C Range Rate')
ylabel('Residual')


