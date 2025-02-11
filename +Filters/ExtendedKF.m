function [XrefHist,PplusHist,measDeltaHist] = ExtendedKF(X0, P0, R, yHist, tVec, stationECI, mu, J2, Re)

%%%%%%%%%%% INPUTS: %%%%%%%%%%
% X0:           [6 x 1] Initial Total State we'll propagate with + STM
% P0:           [6 x 6] Initial uncertainty of initial condition
% Rkf:          [2 x 2] Assumed constant measurement noise (rho,rhoDot)
% yHist:        [Struct] History of measurements for all stations starting
%                   at t0+dt
% tVec:         [1 x n] time history of observations
% traj:         The truth trajectory, only here for debugging reasons
% stationECI:   [] History of all station states in ECI

%%%%%%%%%%% OUTPUTS: %%%%%%%%%
% xHist:         [6 x n] History of Total State Estimate


% --- Measurements ---
% noisey measurements - these are from sensor
rangeMeas    = yHist.Range;
rangeDotMeas = yHist.RangeRate;

% Observed measurements (from sensor)
observedMeasAll = [rangeMeas, rangeDotMeas];

% --- Implement EKF Algorithm ---
% Naming convention
% Prev: Value from previous time step
% Minus: A Priori Value
% Plus: A Posteriori Value
%
% set filter initial conditions
XrefPrev = X0;
pPrev = P0;

% Set filter to loop over number of OBSERVATIONS
for i = 1:length(tVec)-1
    
    % --- Integrate Reference Traj ---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    % Integrate Trajectory
    [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [tVec(i),tVec(i+1)], XrefPrev, odeOptions, mu, J2, Re);
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(end,1:6);
    
    % reference traj position
    refPos = refTrajStates(1:3);
    refVel = refTrajStates(4:6);
    
    % Extract the Integrated STM - this maps previous time to current time
    phi = reshape(TrajNom(end,7:end), [6,6]);
    
    % --- Time Update ---
    pMinus = phi * pPrev * phi';
    
    % --- Computed Measurements --- 
    % each station state
    for j = 1:3
        stationPosECI{i,j} = stationECI{i,j}(1:3);
        stationVelECI{i,j} = stationECI{i,j}(4:6);
    end
    
    % determine station visibility
    [visibilityMask, ~] = Measurements.VisibilityMask(stationPosECI(i,:), refPos', 10, tVec(i));
    
    % Determines which station number made the ob
    statNumOb = find(isnan(visibilityMask) == 0);
    
    if ~isempty(statNumOb)
        % observation was made by a station
        
        if length(statNumOb) > 1
            % if multiple obs then just take one (just for now)
            statNumOb = statNumOb(1);
        else
            % Nothing
        end
        
        % Computed measurement from the station!
        rangeMeasComp     = refPos - cell2mat(stationPosECI(i,statNumOb))';
        rangeNormComp     = norm(rangeMeasComp);
        rangeRateComp     = dot(rangeMeasComp, refVel - cell2mat(stationVelECI(i,statNumOb))') / rangeNormComp;
        
        % stack computed measurements
        computedMeas = [rangeNormComp; rangeRateComp];
        
        % stack observed measurements
        observedMeas = rmmissing(observedMeasAll(i,:));
        
        % Measurement delta: Observed - Computed
        if isempty(observedMeas)
            % If there is no observed measurement bu filter thinks there
            % should be - just go along with and see what happens
            measDelta = computedMeas;
        else
            measDelta = observedMeas' - computedMeas;
        end
        
        % Compute Htilde
        Htilde{i} = Measurements.HtildeSC(refTrajStates', stationECI{i,statNumOb});
        
        % Kalman Gain
        Kk = pMinus*Htilde{i}' / (Htilde{i}*pMinus*Htilde{i}' + R);
        
        xhat = Kk * measDelta;
        
        % reference state + deviation
        Xplus = refTrajStates' + xhat;
        
        % covariance update
        Pplus = (eye(6,6) - Kk*Htilde{i})*pMinus;
        
        assert(trace(Pplus) < trace(pMinus), 'Covariance not decrease!');
        
    else
        % no station made an observation - simply propagate state w/o meas
        Xplus = refTrajStates';
        Pplus = pMinus;
        measDelta = [NaN; NaN];
        
    end
    
    
    % save histories
    XrefHist(:,i) = Xplus;
    PplusHist{i} = Pplus; 
    measDeltaHist(:,i) = measDelta; 
    
    % reset everything for next iteration
    XrefPrev = [Xplus; reshape(eye(6,6), [36,1])];
    pPrev    = Pplus;
    
end

end
