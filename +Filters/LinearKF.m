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


% Function to help sort measurements - these are noisy measurements of ref
[observedMeas] = Measurements.FilterMeasLoadIn(yHist);

% --- Implement LKF Algorithm ---
% Naming convention
% Prev: Value from previous time step
% Minus: A Priori Value
% Plus: A Posteriori Value
%
% set filter initial conditions


% Set filter to loop over number of observations
for i = 1:length(tVec)
    
    if i == 1
        % first loop set up filter
        xhatPrev = pert;
        pPrev = P0;
        TrajNom = IC'; 
        timePrev = 0;
        PhiTotal = reshape(TrajNom(end,7:end), [6,6]);
    else
        
    %--- Integrate ref traj & STM between each time step---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    % Integrate Trajectory
    [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [timePrev:tVec(i)], refState, odeOptions, mu, J2, Re);
    
    % set previous time
    timePrev = tVec(i);
    
    end
    
    timePrev = tVec(i);
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(end,1:6);
    
    % Extract the Integrated STM (maps previous to current time)
    phi = TrajNom(end,7:end);
    
    % reshapre the STM
    STM = reshape(phi, [6,6]);
    
    PhiTotal = STM*PhiTotal; 
    
    % --- Time Update Step
    xMinus = STM * xhatPrev; 
    pMinus = STM * pPrev * STM';
    
    % function to determine which station observed the spacecraft
    [statNumOb] = Measurements.StationObs(visibilityMask, i);
    
    if ~isempty(statNumOb)
        % each column is station
        Htilde{i} = Measurements.HtildeSC(refTrajStates', stationECI{i,statNumOb});
        
        % Computed measurement for filter estimated state
     %   Xcomp = refTrajStates' + xMinus;
        
        % calcualte the computed measurements
        refRangeMeas     = yHistRef.Range(i,statNumOb);
        refRangeRateMeas = yHistRef.RangeRate(i,statNumOb);
        
        compMeas = [refRangeMeas; refRangeRateMeas];        
        
        % pre-fit measurement residuals
        measDelta = rmmissing(observedMeas(i,:))' - compMeas;
        
        % --- Kalman Gain
        Kk = pMinus * Htilde{i}' / (Htilde{i} * pMinus * Htilde{i}' + R);
        
        % -- Measurement Update
        measRes = measDelta - Htilde{i}*xMinus;
        xhatPlus = xMinus + Kk * measRes;
        Pplus = (eye(6,6) - Kk*Htilde{i}) * pMinus * (eye(6,6) - Kk*Htilde{i})' + Kk*R*Kk';
      %  Pplus = (eye(6,6) - Kk*Htilde{i}) * pMinus;
        
        % CHECKING P IS GETTING SMALLER
        assert(trace(Pplus) < trace(pMinus), 'Trace of P is not getting smaller');
        
        % debugging saving this
        estimatedDeviationOb{i} = Htilde{i}*xMinus;
        measDeltaHist{i} = measDelta;
    else
        % no station observation - simply propagating the state w/o meas
        xhatPlus = xMinus;
        Pplus    = pMinus;
        
        % debugging saving this
        estimatedDeviationOb{i} = [NaN; NaN];
        measResHist(:,i) = [NaN; NaN];
        measDeltaHist{i} = [NaN; NaN];
        measRes = [NaN; NaN];
        
    end
    
    % save off histories
    xhist(:, i) = xhatPlus;
    measResHist(:,i) = measRes;
    covMins{i} = pMinus;
    covPlus{i} = Pplus;
    if isempty(statNumOb)
        statNumOb = 0;
    end
    statNumObHist(i) = statNumOb;
    
    % update variables for next go around
    xhatPrev = xhatPlus;
    pPrev    = Pplus;
    refState = [refTrajStates'; reshape(eye(6,6), [36,1])]; 
    
end

end