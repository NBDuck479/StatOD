function [] = LinearKF(IC, pert, P0, R, yHist, yHistRef, yHistTimeTag, stationECI, tVec, mu)

%%%%%%%%%%%%%% INPUTS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC:           [6 x 1] Initial Total State condition that we'll propagate with
% pert:         [6 x 1] Initial Perturbation State 
% P0:           [6 x 6] Initial uncertainty of IC
% R:            [6 x 6] Assumed constant measurement noise for filter
% yHist:        [Struct] History of measurements for all stations starting
%                   at t0+dt
% yHistTimeTag
% tVec:         [1 x n] time history of observations
% mu:           [1 x 1] Scalar of planet gravitational parameter
% stationECI:   [] History of all station states in ECI 

%%%%%%%%%%%%% OUTPUTS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pertHist:      [6 x n] History of perturbation state estimate
% covPostHist:   [] History of final, posterior state
%                uncertainty
% covPriorHist:  [] History of pre-measurement state uncertainty
% measResidHist: [] History of measurement residuals/innovations

% Set Constants
mu = Const.OrbConst.muEarth;
J2 = Const.OrbConst.J2; 
Re = Const.OrbConst.EarthRadius;

%--- Integrate ref traj & STM before running filter --- 

% integrate for the different observation times
tVec = 0:yHistTimeTag;

% Set integrator options
odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);

% Integrate Trajectory 
[T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tVec, IC, odeOptions, mu, J2, Re);

% Extract the reference trajectory states
refTrajStates = TrajNom(:,1:6);

% Extract the Integrated STM - STM is by Row!
phi = TrajNom(:,7:end);

% --- Measurements --- 
% noisey measurements - these are from sensor
rangeMeas    = yHist.Range;
rangeDotMeas = yHist.RangeRate;

% Observed measurements (from sensor)
observedMeas = [rangeMeas, rangeDotMeas];

% reference measurements - these are perfect No Noise!
refRange     = yHistRef.Range;
refRangeRate = yHistRef.RangeRate;

% Computed Measurements (reference)
computedMeas = [refRange; refRangeRate];

% --- Implement LKF Algorithm 
% Naming convention
% Prev: Value from previous time step
% Minus: A Priori Value
% Plus: A Posteriori Value
% 
% set filter initial conditions
xhatPrev = pert;
pPrev = P0;
timePrev = 0;

% Set filter to loop over number of observations
for i = length(yhist)

    % Get the time tag for this Observation
    time = yHistTimeTag(i) - timePrev;
    
    % Setup phi for run
    if i == 1
        % For first go phi is IC 
        STM = reshape(phi(i,:), [6,6]);
        
    else
        % After first go phi multiplies itself 
        STM = reshape(phi(i,:), [6,6])*reshape(phi(i-1, :), [6,6]);
        
    end
    
    
    % --- Time Update Step 
    xMinus = STM * xhatPrev; 
    pMinus = STM * pPrev * STM'; 
    
    % -- Observation Deviation
    measResidual = observedMeas(:,i) - computedMeas(:,i);
    
    % --- Observation State Matrix
    
    
    % Set time for next iteration
    timePrev = time; 
end