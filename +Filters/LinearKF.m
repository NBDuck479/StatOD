function [xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covPlus] = LinearKF(LKFinputs)

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

% --- Dynamically unpack the input struct for LKF

% Get all field names in the struct
fields = fieldnames(LKFinputs);

% Loop through each field and access its value dynamically
for i = 1:length(fields)
    % Get the field name
    fieldName = fields{i};
    
    % Access the value of the field using dynamic field referencing
    fieldValue = LKFinputs.(fieldName);
    
    eval([fieldName ' = fieldValue']);
    
end


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
for i = 1:length(tOverall)
    
    if i == 1
        % first loop set up filter
        xhatPrev = pert;
        pPrev = P0;
        TrajNom = IC'; 
        timePrev = 0;
        PhiTotal = reshape(TrajNom(end,NumStates+1:end), [NumStates,NumStates]);
    else
        
    %--- Integrate ref traj & STM between each time step---
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    % Integrate Trajectory
    [T, TrajNom] = ode45(@Dynamics.Numeric_J2_Drag_Prop, [timePrev:tOverall(i)], refState, odeOptions, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0);
    
    % propagated stations in ECEF so rotate those states into ECI!
    
    % set previous time
    timePrev = tOverall(i);
    
    end
    
    timePrev = tOverall(i);

    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(end,1:NumStates);
    
    % Extract the Integrated STM (maps previous to current time)
    phi = TrajNom(end,NumStates+1:end);
    
    % reshapre the STM
    STM = reshape(phi, [NumStates,NumStates]);
    
    PhiTotal = STM*PhiTotal; 
    
    % --- Time Update Step
    xMinus = STM * xhatPrev; 
    pMinus = STM * pPrev * STM';
    
    % determine if time aligns with observation
    if ismember(tOverall(i), obsHist.time)
        
        % get index of observaiton
        obInd = find(obsHist.time == tOverall(i));
        
        % get station number 
        statNumOb = obsHist.statNo(obInd); 
        
        % each column is station
        Htilde{i} = Measurements.HtildeSCProj1(refTrajStates', stationECI{obInd,statNumOb}, statNumOb, MeasFlag);
        
        % Computed measurement for filter estimated state
     %   Xcomp = refTrajStates' + xMinus;
        
        % calcualte the computed measurements
        refRangeMeas     = yHistRef.Range(obInd,statNumOb);
        refRangeRateMeas = yHistRef.RangeRate(obInd,statNumOb);
        
        compMeas = [refRangeMeas; refRangeRateMeas];        
        
        % pre-fit measurement residuals
        measDelta = rmmissing(observedMeas(obInd,:))' - compMeas;
        
        % --- Kalman Gain
        Kk = pMinus * Htilde{i}' / (Htilde{i} * pMinus * Htilde{i}' + R);
        
        % -- Measurement Update
        measRes = measDelta - Htilde{i}*xMinus;
        xhatPlus = xMinus + Kk * measRes;
        Pplus = (eye(NumStates,NumStates) - Kk*Htilde{i}) * pMinus * (eye(NumStates,NumStates) - Kk*Htilde{i})' + Kk*R*Kk';
      %  Pplus = (eye(6,6) - Kk*Htilde{i}) * pMinus;
        
        % CHECKING P IS GETTING SMALLER
    %    assert(trace(Pplus) < trace(pMinus), 'Trace of P is not getting smaller');
        
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
    refState = [refTrajStates'; reshape(eye(NumStates,NumStates), [NumStates^2,1])]; 
    
end

end