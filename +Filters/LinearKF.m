function [xhist, measResHist, measDeltaHist, estimatedDeviationOb, statNumObHist, covMins, covPlus, refTrajStatesHist, STMhist, PhiTotal] = LinearKF(LKFinputs)

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
        
        % initialize as empty in case no station makes imeediate ob
        statNumOb = [];
        
        % initialize time of first observation
        prevObTime = [];
    else
        
        %--- Integrate ref traj & STM between each time step---
        % Set integrator options
        odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
        
        if DMC == 1
            % integrate using DMC
            [T, TrajNom] = ode45(@Dynamics.NumericJ2PropDMC, [timePrev:tOverall(i)], refState, odeOptions, mu, J2, Re, tau);
        else
            % integrate without DMC
            % Integrate Trajectory
            [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [timePrev:tOverall(i)], refState, odeOptions, mu, J2, Re, DMC);
            
            %      [T,TrajNom] = ode45(@Dynamics.DynamicsA_J2_J3, [timePrev:tOverall(i)], refState, odeOptions, mu, J2 , -2.5323e-06, Re);
            
        end
    end
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(end,1:NumStates);
    
    % Extract the Integrated STM (maps previous to current time)
    phi = TrajNom(end,NumStates+1:end);
    
    % reshapre the STM
    STM = reshape(phi, [NumStates,NumStates]);
    
    PhiTotal = STM*PhiTotal;
    
    if DMC == 1
        % include DMC to time update
        
        % Determine how long it's been since last observation
        if tOverall(i) - prevObTime > 10
            Q = zeros(NumStates, NumStates);
        else
            % calculate Q
            Q = Dynamics.J2_DMC_Q_Matrix(tOverall(i) - timePrev, tau, sigmaDMC);
        end
        
        % --- Time Update
        pMinus = STM * pPrev * STM' + Q;
        
    else
        % Add SNC instead of DMC
        
        % Determine if time gap is too large for added SNC
        if tOverall(i) - prevObTime > 10
            GammaQGamma = zeros(NumStates, NumStates);
        else
            % time is small enough to add SNC
            [GammaQGamma, ~] = Dynamics.StateNoiseComp(tOverall(i) - timePrev, Q, refTrajStates, Qframe);
        end
        
        % SNC gets added to covariance update
        pMinus = STM * pPrev * STM' + GammaQGamma;
        
    end
    
    
    % --- Time Update Step
    xMinus = STM * xhatPrev;
    
    % determine if time aligns with observation
    if ismember(tOverall(i), yHist.obTime)
        
        % - find indices of where the times match
        
        % Loop over each element of tOverall and find the matching indices
        % Find indices where tOverall(i) matches elements in yHist.obTime
        [row, col] = find(yHist.obTime == tOverall(i));
        
        % the row is the index that matches with time
        obInd = row;
        % the column is the station number
        statNumOb = col;
        
        % each column is station
        Htilde{i} = Measurements.HtildeSC(refTrajStates', stationECI{i,statNumOb}, MeasFlag);
        
        if DMC == 1
            Htilde{i} = [Htilde{i}, zeros(2,3)];
        end
        
        % Computed measurement for filter estimated state
        %   Xcomp = refTrajStates' + xMinus;
        
        % calcualte the computed measurements
        refRangeMeas     = yHistRef.Range(i,statNumOb);
        refRangeRateMeas = yHistRef.RangeRate(i,statNumOb);
        
        compMeas = [refRangeMeas; refRangeRateMeas];
        
        % measurements from all stations at time
        fullObsMeas = observedMeas(i,:)';
        
        % range and range rate at time
        ObsMeasRange     = fullObsMeas(1:3);
        ObsMeasRangerate = fullObsMeas(4:6);
        
        % pre-fit measurement residuals
        measDelta = [ObsMeasRange(statNumOb); ObsMeasRangerate(statNumOb)] - compMeas;
        
        % --- Kalman Gain
        Kk = pMinus * Htilde{i}' / (Htilde{i} * pMinus * Htilde{i}' + R);
        
        % -- Measurement Update
        measRes = measDelta - Htilde{i}*xMinus;
        xhatPlus = xMinus + Kk * measRes;
        Pplus = (eye(NumStates,NumStates) - Kk*Htilde{i}) * pMinus * (eye(NumStates,NumStates) - Kk*Htilde{i})' + Kk*R*Kk';
        
        % CHECKING P IS GETTING SMALLER
        assert(trace(Pplus) < trace(pMinus), 'Trace of P is not getting smaller');
        
        % save this obsevration time as the previous observation time
        prevObTime = yHist.obTime(obInd);
        
        % debugging saving this
        estimatedDeviationOb{i} = Htilde{i}*xMinus;
        measDeltaHist(1:2,i) = measDelta;
    else
        % no station observation - simply propagating the state w/o meas
        xhatPlus = xMinus;
        Pplus    = pMinus;
        
        % debugging saving this
        estimatedDeviationOb{i} = [NaN; NaN];
        measResHist(:,i) = [NaN; NaN];
        measDeltaHist(1:2,i) = [NaN; NaN];
        measRes = [NaN; NaN];
        
    end
    
    % save off histories
    xhist(:, i) = xhatPlus;
    measResHist(:,i) = measRes;
    covMins{i} = pMinus;
    covPlus{i} = Pplus;
    STMhist{i} = STM;
    if isempty(statNumOb)
        statNumOb = 0;
    end
    statNumObHist(i) = statNumOb;
    
    % update variables for next go around
    xhatPrev = xhatPlus;
    pPrev    = Pplus;
    refState = [refTrajStates'; reshape(eye(NumStates,NumStates), [NumStates^2,1])];
    timePrev = tOverall(i);
    
    refTrajStatesHist(i,:) = refTrajStates;
    
end

% pre-fit residual RMS
measDeltaRMS_Rho    = rmmissing(measDeltaHist(1,:));
measDeltaRMS_RhoDot = rmmissing(measDeltaHist(2,:));

RMS_pre_Rho    = sqrt(sum(measDeltaRMS_Rho.^2)) / length(measDeltaRMS_Rho);
RMS_pre_RhoDot = sqrt(sum(measDeltaRMS_RhoDot.^2)) / length(measDeltaRMS_RhoDot);

% post-fit residual RMS
measResRMS_Rho    = rmmissing(measResHist(1,:));
measResRMS_RhoDot = rmmissing(measResHist(2,:));

RMS_post_Rho    = sqrt(sum(measResRMS_Rho.^2)) / length(measResRMS_Rho);
RMS_post_RhoDot = sqrt(sum(measResRMS_RhoDot.^2)) / length(measResRMS_RhoDot);

% Display the RMS residual results for the LKF
% print out what is happening
fprintf('--- Residual RMS info: ---\n');
fprintf('Pre-fit residual RMS values\n');
fprintf('Rho    = %g km\n', RMS_pre_Rho );
fprintf('RhoDot = %g km\n', RMS_pre_RhoDot );
fprintf('Post-fit residual RMS values\n');
fprintf('Rho    = %g km\n', RMS_post_Rho );
fprintf('RhoDot = %g km\n', RMS_post_RhoDot );

fprintf('\n---------------------------------------------\n\n');

end