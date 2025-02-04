function [] = ExtendedKF(X0, P0, Rkf, yHist, tVec, traj, stationECI, mu, J2, Re)

%%%%%%%%%%% INPUTS: %%%%%%%%%%
% X0:           [6 x 1] Initial Total State we'll propagate with
% P0:           [6 x 6] Initial uncertainty of initial condition
% Rkf:          [2 x 2] Assumed constant measurement noise (rho,rhoDot)
% yHist:        [Struct] History of measurements for all stations starting
%                   at t0+dt
% tVec:         [1 x n] time history of observations
% traj:         The truth trajectory, only here for debugging reasons
% stationECI:   [] History of all station states in ECI

%%%%%%%%%%% OUTPUTS: %%%%%%%%%
% xHist:         [6 x n] History of Total State Estimate


% function to help sort measurements
[computedMeas, observedMeas] = Measurement.FilterMeasLoadIn(yHist);

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
for i = 1:length(tVec)
    
    % get change in time
    dt = tVec(k+1) - tVec(k);
    
    % --- Integrate Reference Traj ---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    % Integrate Trajectory
    [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tVec, XrefPrev, odeOptions, mu, J2, Re);
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(:,1:6);
    
    % Extract the Integrated STM - STM is by Row!
    phi = TrajNom(:,7:end);
    
    % --- Time Update ---
    pMinus = phi * pPrev * phi';
    
    % function to determine which filter observed and calc measurement
    % delta
    [statNumOb, measDelta] = Measurements.StationObs(observedMeas, computedMeas, visibilityMask, i);
    
    % Check if station made an observation
    if ~isempty(statNumOb)
        % Compute Htilde
        Htilde{i} = Measurements.HtildeSC(refTrajStates(i,:)', stationECI{i,statNumOb});
        
        % Kalman Gain
        Kk = pMinus*Htilde{i}' * inv(Htilde{i}*pMinus*Htilde{i}' + R);
        
        xhat = Kk * measDelta;
        
        % reference state + deviation
        Xplus = refTrajStates(i,:) + xhat;
        
        % covariance update
        Pplus = (eye(6,6) - Kk*Htilde{i}) * pPrev * (eye(6,6) - Kk*Htilde{i})' + Kk*R*Kk';
        
    else
        % no station made an observation - simply propagate state w/o meas
        Xplus = refTrajStates(i,:);
        Pplus = pMinus;
        
    end
    
    
    % save histories
    XrefHist(:,i) = Xplus;
    
    
    % reset everything for next iteration
    XrefPrev = Xplus;
    pPrev    = Pplus;
    
end

end
