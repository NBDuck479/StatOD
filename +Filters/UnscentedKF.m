function [Xhist, Phist, OminusC] = UnscentedKF(UKFinputs)

% This implements a UKF for state esimtation also known as a sigma point
% filter.
% It does a really good job with estimating the second moment (covariance)
% of the distribution. It doesn't just linearly propagate the nonlinear
% coaviarnce, it actually tries to use the sigma points to reconstruct the
% distribution.

% --- Dynamically unpack the input struct for UKF

% Get all field names in the struct
fields = fieldnames(UKFinputs);

% Loop through each field and access its value dynamically
for i = 1:length(fields)
    % Get the field name
    fieldName = fields{i};
    
    % Access the value of the field using dynamic field referencing
    fieldValue = UKFinputs.(fieldName);
    
    eval([fieldName ' = fieldValue']);
    
end


% --- Measurements ---
% noisey measurements - these are from sensor
rangeMeas    = yHist.Range;
rangeDotMeas = yHist.RangeRate;

% Observed measurements (from sensor)
observedMeasAll = [rangeMeas, rangeDotMeas];


% --- Compute Weights for Sigma Points ---
% Weight Parameters
L       = NumStates;
kappa   = 3 - L;
lambda  = a^2*(L+kappa)-L;
gamma   = sqrt(L+lambda);

% Weight for Mean and Covariance
Wm(1) = lambda / (L+lambda);
Wc(1) = (lambda / (L+lambda)) + (1 - a^2 + B);

% add weights for other sigma points
for j = 2:2*L+1
    Wm(j) = 1 / (2*(L+lambda));
    Wc(j) = 1 / (2*(L+lambda));
end

% weights for mean must be 1
assert(.999 < sum(Wm) < 1.0001, 'Weights of mean sigma points must be 1')

% --- Initialize outside of Filter
% set filter initial conditions
Xprev = X0;
sqrtPprev = sqrtm(P0);
timePrev = 0;

% Set time for last observation
obTimePrev = [];

% Set filter to loop over number of observations
for i = 2:length(tOverall)
    
    % --- Compute Sigma Points
    % Use previous time step results
    chiPrev = [Xprev(1:NumStates), Xprev(1:NumStates)+gamma*sqrtPprev, Xprev(1:NumStates)-gamma*sqrtPprev];
    
    % --- Propagate Sigma Points ---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    for SigPt = 1:2*L+1
        [T, propSigmaPt] = ode45(@Dynamics.NumericJ2Prop, [timePrev,tOverall(i)], chiPrev(:,SigPt), odeOptions, mu, J2, Re);
        
        % Save the propagated Sigma points
        chiMinus(SigPt,:) = propSigmaPt(end,:);
    end
    
    % --- Time update ---
    weightedSigmaPt = zeros(NumStates,1);
    Xprev           = zeros(NumStates,1);
    CovNoPNloop     = zeros(NumStates,NumStates);
    CovNoPN         = zeros(NumStates,NumStates);
    
    % For loop State Time Update
    for q = 1:2*L+1
        % Loop over each sigma point q for State
        weightedSigmaPt = Wm(q) * chiMinus(q,:);
        
        % Summation
        Xprev = weightedSigmaPt' + Xprev;
    end
    
    for w = 1:2*L+1
        % Loop each sigma point for covariance
        CovNoPNloop = Wc(w) * (chiMinus(w,:)' - Xprev) * (chiMinus(w,:)' - Xprev)';
        
        % Summation
        CovNoPN = CovNoPN + CovNoPNloop;
    end
    
    % State Noise Compensation
    if tOverall(i) - obTimePrev < 15
        [Q_s, Gamma] = Dynamics.StateNoiseComp(tOverall(i) - timePrev, Q, 0, Qframe);
    else
        Q_s = zeros(NumStates, NumStates);
    end
    
    % Covariance Update with process noise
    Pprev = Q_s + CovNoPN;
    
    % ENSURE MEASUREMENT HAPPENS AT THIS TIME!!!
    % determine if time aligns with observation
    if ismember(tOverall(i), yHist.obTime)
        % --- Recompute Sigma Points after Propagation
        ChiMinusProp = [Xprev, Xprev+gamma*sqrtm(Pprev), Xprev-gamma*sqrtm(Pprev)];
        
        % Loop over each element of tOverall and find the matching indices
        % Find indices where tOverall(i) matches elements in yHist.obTime
        [row, col] = find(yHist.obTime == tOverall(i));
        
        % the row is the index that matches with time
        obInd = row;
        % the column is the station number
        statNumOb = col;
        
        % --- Compute the Measurements at Time
        [computedMeas] = Measurements.ComputeFilterMeasurements(stationECI, i, statNumOb, ChiMinusProp);
        
        % Initialize before loop
        meanPredMeas = zeros(2,1);
        
        % Mean Predicted Measurement
        for SP = 1:2*L+1
            meanPredMeas = Wm(SP) * computedMeas(:,SP) + meanPredMeas;
        end
        
        Pyy = zeros(2,2);
        Pxy = zeros(6,2);
        % --- Innovation and Cross Covariances
        for innov = 1:2*L+1
            % Innovation Covariance
            Pyy = Wc(innov) * (computedMeas(:,innov) - meanPredMeas) * (computedMeas(:,innov) - meanPredMeas)' + Pyy;
            
            % Cross Covairnace
            Pxy = Wc(innov) * (ChiMinusProp(:,innov) - Xprev) * (computedMeas(:,innov) - meanPredMeas)' + Pxy;
        end
        % Innovation Covariance
        Pyy = R + Pyy;
        
        % --- Actual Sensor observation for this time
        ObsMeasRange    = yHist.Range(obInd, statNumOb);
        ObMeasRangeRate = yHist.RangeRate(obInd, statNumOb);
        
        ObMeasurement = [ObsMeasRange; ObMeasRangeRate];
        
        % --- Compute Kalman Gain
        Kk = Pxy * inv(Pyy);
        
        % --- Measurement Update
        Xhat = Xprev + Kk * (ObMeasurement - meanPredMeas);
        
        P = Pprev - Kk*Pyy*Kk';
        
        % save off residual 
        OminusC(:,i) = ObMeasurement - meanPredMeas;
        
        % set previous observation time 
        obTimePrev = tOverall(i);
        
    else
        % No observation at this time
        Xhat = Xprev;
        P = Pprev;
    end
    
    % reset everything for next go around
    Xprev = Xhat;
    Pprev = P;
    sqrtPprev = sqrtm(Pprev);
    timePrev = tOverall(i);
    
    % Save off Histories
    Xhist(:,i-1) = Xhat;
    Phist{i-1} = P;
    
    
end