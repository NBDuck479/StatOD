function [] = LinearKF(IC, pert, P0, R, yHist, yHistTimeTag, tVec, mu)

%%%%%%%%%%%%%% INPUTS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IC:           [6 x 1] Initial Total State condition that we'll propagate with
% pert:         [6 x 1] Initial Perturbation State 
% P0:           [6 x 6] Initial uncertainty of IC
% R:            [6 x 6] Assumed constant measurement noise for filter
% yHist:        [SIZE x SIZE] History of measurements for all stations starting
%                   at t0+dt
% yHistTimeTag
% tVec:         [1 x n] time history of observations
% mu:           [1 x 1] Scalar of planet gravitational parameter

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

%---- PASS THIS INTO FUNC Set Initial Conditions
r0   = [-3515.49032703351; 8390.71631024339; 4127.62735255368];
v0   = [-4.35767632217815; -3.35657913876455; 3.1118929278699];
Phi0 = eye(6,6);

% integrate for the different observation times
tVec = 0:yHistTimeTag;

% Stack all initial conditions into vector
IC = [r0; v0; reshape(Phi0, [36,1])];

% Set integrator options
odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);

% Integrate Trajectory 
[T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tVec, IC, odeOptions, mu, J2, Re);

% Extract the reference trajectory states
refTrajStates = TrajNom(:,1:6);

% Extract the Integrated STM - STM is by Row!
phi = TrajNom(:,7:end);


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
    
    % -- Computed Observation
    Measurements.VisibilityMask(stationPos, scPos, elevAngle, tSpan)
    
    
    % Set time for next iteration
    timePrev = time; 
end