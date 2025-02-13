function [xhat0est, P0est, phiHist, measDeltaHist, Htilde] = BatchFilter(IC, pert, P0, R, yHist, yHistRef, stationECI, visibilityMask, tVec, mu, J2, Re)
% Batch filter to process spaceceraft observation measurements
%
%%%%%%%%%% INPUTS %%%%%%%%%%
% IC        Initial full state and STM

% Function to help sort measurements - these are noisy measurements of ref
[observedMeas] = Measurements.FilterMeasLoadIn(yHist);

% set starting condition for filter
Xprev   = IC(1:6);
STMprev = IC(7:end);

% Build first intitial condition
XrefPrev = [Xprev; STMprev];

% set a priori estimates
delta = inv(P0);
N = inv(P0)\pert;

% Integrate reference trajectory and STM
% Set integrator options
odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);

% Integrate Trajectory for time step
[T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, tVec, XrefPrev, odeOptions, mu, J2, Re);

% Process each observation
for i = 1:length(tVec)-1
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(i+1,1:6);
    
    % Extract the Integrated STM - STM is by Row!
    phi = TrajNom(i+1,7:end);
    
    STM = reshape(phi, [6,6]);
    
    % keep track of STM mapped back to time zero
    
    % --- Accumulate current observation
    % function to determine which filter observed and calc measurement
    % delta
    [statNumOb] = Measurements.StationObs(visibilityMask, i);
    
    % Check if station made an observation
    if ~isempty(statNumOb)
        % Compute Htilde
        Htilde{i} = Measurements.HtildeSC(refTrajStates', stationECI{i,statNumOb});
        
        % calcualte the computed measurements
        refRangeMeas     = yHistRef.Range(i,statNumOb);
        refRangeRateMeas = yHistRef.RangeRate(i,statNumOb);
        
        compMeas = [refRangeMeas; refRangeRateMeas];
        
        % measurement delta
        measDelta = rmmissing(observedMeas(i,:))' - compMeas;
        
        % propagate Htilde
        H = Htilde{i} * STM;
        
        % update delta and N
        delta = delta + H'*inv(R)*H;
        N = N + H'*inv(R)*measDelta;
    else
        % no ob given
        Htilde{i} = [NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN];
    end
    
    if i < length(tVec)/2-2
        % Go back and read in the next observation
        
    else
        % Finished processing!
        
        % --- Solve Normal Equations
        
        % invert the information matrix, delta
        deltaChol = chol(delta);
        
        % solve for initial pert condition
        xhat0est = deltaChol\(deltaChol'\N);
        
        % check the convergence of xhat0
        
        % covariance
        P0est = inv(deltaChol);
        
        
    end % end of obs measurements
    
    % save off some histories 
    phiHist(i, :) = phi; 
    measDeltaHist(:,i) = measDelta;
    
end