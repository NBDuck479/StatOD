%% StatOD Homework 5 
addpath('C:\Users\user\OneDrive - UCB-O365\Desktop\Spring 2025\StatOD\repo')
mu = Const.OrbConst.muEarth;
Re = Const.OrbConst.EarthRadius; 
J2 = Const.OrbConst.J2;


% Import Sensor Measurements from previous homeworks
yHistLoad = load('yHist.mat');
yHist = yHistLoad.yHist; 

% from previous homework to make sensor observations
sigmaRange     = 0.001; % km
sigmaRangeRate = 0.000001; % km/s

% create sensor measurement noise covariance matrox
R = [sigmaRange^2 0; 0  sigmaRangeRate^2]; 

% Going to use V for whitening - pass into filter and whiten on the spot?
V = chol(R);


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

%% ---  Implement the Square Root Information Filter

% SRIF input struct
SRIFinputs.IC                = Y0; %[Y0(1:6); 0;0;0; reshape(eye(9,9), [9^2, 1])]; % Y0; 
SRIFinputs.NumStates         = NumStates; 
SRIFinputs.pert              = zeros(6,1); 
SRIFinputs.P0                = P0; %[P0 zeros(6,3); zeros(3,6) blkdiag(0.000001^2, 0.000001^2, 0.000001^2)]; %P0; 
SRIFinputs.V                 = V; % Whitening  
SRIFinputs.mu                = mu;  
SRIFinputs.J2                = J2;
SRIFinputs.yHist             = yHist; 
SRIFinputs.stationECI        = stationECI;
SRIFinputs.visibilityMask    = visibilityMask;
SRIFinputs.tOverall          = tOverall; 
SRIFinputs.Re                = Re; 
SRIFinputs.MeasFlag          = 3; % 1 for only range, 2 for range rate, 3 for all
% SRIFinputs.Q                 = covSNC; 
% SRIFinputs.Qframe            = 'RSW';
% SRIFinputs.DMC               = 0; % 1 for DMC and 0 for none
% SRIFinputs.tau               = [period/30; period/30; period/30]; % time constants for DMC
% SRIFinputs.sigmaDMC          = [0.0001; 0.0001; 0.0001]; % uncertainty for DMC accel error states

% Run the SRIF
[xhat, P] = Filters.SRIF(SRIFinputs)
