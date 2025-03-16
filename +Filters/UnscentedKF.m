function [] = UnscentedKF(UKFinputs)
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
Wc(1) = lambda / (L+lambda) + (1 - a^2 + B);

% add weights for other sigma points
for j = 2:2*L+1
    Wm(j) = 1 / (2*(L+lambda));
    Wc(j) = 1 / (2*(L+lambda)); 
end

% weights for mean must be 1
assert(.999 < sum(Wm) < 1.0001, 'Weights of mean sigma points must be 1')

% Set filter to loop over number of observations
for i = 1:length(tOverall)
    
    if i == 1
        % set filter initial conditions
        XrefPrev = X0;
        sqrtPprev = chol(P0);
        timePrev = 0;
        
        % Set time for last observation
        obTimePrev = [];
        
        % get initial phi
        phi = reshape(XrefPrev(NumStates+1:end), [NumStates,NumStates]);
        
        % initial reference position
        refPos = XrefPrev(1:3);
        refVel = XrefPrev(4:6);
        
    else
        % --- Integrate Reference Traj ---
        % Set integrator options
        odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
        [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [timePrev,tVec(i)], XrefPrev, odeOptions, mu, J2, Re);
    end
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(end,1:NumStates);
    
    % reference traj position
    refPos = refTrajStates(1:3);
    refVel = refTrajStates(4:6);
    
    % Extract the Integrated STM - this maps previous time to current time
    phi = reshape(TrajNom(end,NumStates+1:end), [NumStates,NumStates]);
    
    % --- Compute Sigma Points
    chiPrev = [XrefPrev, XrefPrev+gamma*sqrtPprev, XrefPrev+gamma*sqrtPprev];
    
    % --- Propagate Sigma Points ---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    for SigPt = 1:2*L+1
        [T, propSigmaPt] = ode45(@Dynamics.NumericJ2Prop, [timePrev,tVec(i)], chiPrev(SigPt), odeOptions, mu, J2, Re);
        
        % Save the propagated Sigma points 
        chiMinus(SigPt,:) = propSigmaPt(end,:); 
    end
    
    % --- Time update ---
    weightedSigmaPt = [];
    Xprev           = [];
    CovNoPNloop         = [];
    
    % For loop State Time Update
    for q = 1:2*L+1
        % Loop over each sigma point q for State
        weightedSigmaPt = Wm(q) * chiMinus(q);
        
        % Summation
        Xprev = weightedSigmaPt + Xprev;
        
        % Loop each sigma point for covariance
        CovNoPNloop = Wc(q) * (chiMinus(q) - Xprev) * (chiMinus(q) - Xprev)';
        
        % Summation
        CovNoPN = CovNoPN + CovNoPNloop;
    end
    
    % Covariance Update with process noise
    Pprev = Q_s + CovNoPN; 

    % --- Recompute Sigma Points after Propagation
    ChiMinusProp = [Xprev, Xprev+gamma*chol(Pprev), Xprev-gamma*chol(Pprev)]; 
    
    
    % Computed Measurements (predicted) for each sigma point
    
        % --- Calculate Computed measurement
        rangeMeasComp     = refPos - stationPosECI';
        rangeNormComp     = norm(rangeMeasComp);
        rangeRateComp     = dot(rangeMeasComp, refVel - stationVelECI') / rangeNormComp;
    
    
    
end