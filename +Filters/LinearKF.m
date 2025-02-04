function [xhist, measResHist] = LinearKF(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, tVec, mu, J2, Re)

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

%--- Integrate ref traj & STM before running filter ---


% Set integrator options
odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);

% Integrate Trajectory
[T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tVec, IC, odeOptions, mu, J2, Re);

% Extract the reference trajectory states
refTrajStates = TrajNom(:,1:6);

% Extract the Integrated STM - STM is by Row!
phi = TrajNom(:,7:end);

% Function to help sort measurements
[computedMeas, observedMeas] = Measurements.FilterMeasLoadIn(yHist);

% --- Implement LKF Algorithm ---
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
for i = 1:length(tVec)
    
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
    
    % function to determine which filter observed and calc measurement
    % delta
    [statNumOb, measDelta] = Measurements.StationObs(observedMeas, computedMeas, visibilityMask, i);
    
    if ~isempty(statNumOb)
        % each column is station
        Htilde{i} = Measurements.HtildeSC(refTrajStates(i,:)', stationECI{i,statNumOb});
        
        % --- Kalman Gain
        Kk = pMinus * Htilde{i}' * inv(Htilde{i} * pMinus * Htilde{i}' + R);
        
        % -- Measurement Update
        measRes = (measDelta(:,statNumOb) - Htilde{i}*xMinus);
        xhatPlus = xMinus + Kk * measRes;
        Pplus = (eye(6,6) - Kk*Htilde{i}) * pPrev * (eye(6,6) - Kk*Htilde{i})' + Kk*R*Kk';
    else
        % no station observation - simply propagating the state w/o meas
        xhatPlus = xMinus;
        Pplus    = pMinus;
        
    end
    
    % save off histories
    xhist(:, i) = xhatPlus;
    measResHist(:,i) = measRes;
    
    % update variables for next go around
    xhatPrev = xhatPlus;
    pPrev    = Pplus;
    
end

end