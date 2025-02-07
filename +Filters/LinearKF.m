function [xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = LinearKF(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, tVec, mu, J2, Re)

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


% Function to help sort measurements
[computedMeas, observedMeas] = Measurements.FilterMeasLoadIn(yHist, yHistRef);

% --- Implement LKF Algorithm ---
% Naming convention
% Prev: Value from previous time step
% Minus: A Priori Value
% Plus: A Posteriori Value
%
% set filter initial conditions
xhatPrev = pert;
pPrev = P0;
refState = IC;
timePrev = 0;

% Set filter to loop over number of observations
for i = 1:length(tVec)-1
        
    %--- Integrate ref traj & STM between each time step---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    % Integrate Trajectory
    [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [tVec(i):tVec(i+1)], refState, odeOptions, mu, J2, Re);
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(1:end,1:6);
    
    % Extract the Integrated STM - STM is by Row!
    phi = TrajNom(:,7:end);
    
    
    
    
    
    
    
    if i == 1
        % just grab STM from propagated ref traj
        STM = reshape(phi(i,:), [6,6]);
    else
        % Want phi(ti, ti-1) on each iteration
        
        % invert previous STM
        prevSTMinv = inv(reshape(phi(i, :), [6,6]));
        
        % new STM maps from prev time step to current
        STM = reshape(phi(i+1,:), [6,6]) * prevSTMinv;
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
        Kk = pMinus * Htilde{i}' / (Htilde{i} * pMinus * Htilde{i}' + R);
        
        % -- Measurement Update
        measRes = measDelta(:,statNumOb) - Htilde{i}*xMinus;
        xhatPlus = xMinus + Kk * measRes;
        Pplus = (eye(6,6) - Kk*Htilde{i}*STM) * pPrev * (eye(6,6) - Kk*Htilde{i}*STM)' + Kk*R*Kk';
        % Pplus = (eye(6,6) - Kk*Htilde{i}) * pMinus;
        
        % CHECKING P IS GETTING SMALLER
        % assert(trace(Pplus) < trace(pMinus), 'Trace of P is not getting smaller');
        
        % debugging saving this
        estimatedDeviationOb{i} = Htilde{i}*xMinus;
    else
        % no station observation - simply propagating the state w/o meas
        xhatPlus = xMinus;
        Pplus    = pMinus;
        
        % debugging saving this
        estimatedDeviationOb{i} = [NaN; NaN];
        measResHist(:,i) = [NaN; NaN];
        
    end
    
    % save off histories
    xhist(:, i) = xhatPlus;
    measResHist(:,i) = measRes;
    covMins{i} = pMinus;
    covPlus{i} = Pplus;
    measDeltaHist{i} = measDelta;
    if isempty(statNumOb)
        statNumOb = 0;
    end
    statNumObHist(i) = statNumOb;
    
    % update variables for next go around
    xhatPrev = xhatPlus;
    pPrev    = Pplus;
    
end

end