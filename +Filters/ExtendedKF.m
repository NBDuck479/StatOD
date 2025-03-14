function [XrefHist,PplusHist,measDeltaHist, post_fit_meas] = ExtendedKF(EKFinputs)

%%%%%%%%%%% INPUTS: %%%%%%%%%%
% X0:           [6 x 1] Initial Total State we'll propagate with + STM
% P0:           [6 x 6] Initial uncertainty of initial condition
% Rkf:          [2 x 2] Assumed constant measurement noise (rho,rhoDot)
% yHist:        [Struct] History of measurements for all stations starting
%                   at t0+dt
% tVec:         [1 x n] time of whole simulation
% traj:         The truth trajectory, only here for debugging reasons
% stationECI:   [] History of all station states in ECI

%%%%%%%%%%% OUTPUTS: %%%%%%%%%
% xHist:         [6 x n] History of Total State Estimate


% --- Dynamically unpack the input struct for LKF

% Get all field names in the struct
fields = fieldnames(EKFinputs);

% Loop through each field and access its value dynamically
for i = 1:length(fields)
    % Get the field name
    fieldName = fields{i};
    
    % Access the value of the field using dynamic field referencing
    fieldValue = EKFinputs.(fieldName);
    
    eval([fieldName ' = fieldValue']);
    
end

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

% Set filter to loop 
for i = 1:length(tVec)
    
    if i == 1
        % set filter initial conditions
        XrefPrev = X0;
        pPrev = P0;
        timePrev = 0;
        
        % Set time for last observation
        obTimePrev = [];
        
        % get initial phi 
        phi = reshape(XrefPrev(NumStates+1:end), [NumStates,NumStates]);
        
        % initial reference position
        refPos = XrefPrev(1:3);
        refVel = XrefPrev(4:6);
        
        if DMC == 1
            refTrajStates = [refPos', refVel', X0(7:9)'];
        else
            refTrajStates = [refPos', refVel'];
        end
        
    else
        
        % --- Integrate Reference Traj ---
        % Set integrator options
        odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
        
        % Choose integration with or without DMC
        if DMC == 1
            % add the accel error to state   
            
            % integrate with DMC
            [T, TrajNom] = ode45(@Dynamics.NumericJ2PropDMC, [timePrev,tVec(i)], XrefPrev, odeOptions, mu, J2, Re, tau);

            assert(NumStates == 9, 'wrong number of states for DMC')
        else
            % integrate without DMC 
            % Integrate Trajectory
            [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [timePrev,tVec(i)], XrefPrev, odeOptions, mu, J2, Re);
        %    [T, TrajNom] = ode45(@Dynamics.DynamicsA_J2_J3, [timePrev,tVec(i)], XrefPrev, odeOptions, mu, J2, -2.5323e-06, Re);
            assert(NumStates == 6, 'wrong number of states for no DMC')
        end
        
        % Extract the reference trajectory states
        refTrajStates = TrajNom(end,1:NumStates);
        
        % reference traj position
        refPos = refTrajStates(1:3);
        refVel = refTrajStates(4:6);
        
        % Extract the Integrated STM - this maps previous time to current time
        phi = reshape(TrajNom(end,NumStates+1:end), [NumStates,NumStates]);
        
    end
    
    % --- Time update --- 
    if DMC == 1
        % include DMC
        % calculate Q 
        Q = Dynamics.J2_DMC_Q_Matrix(tVec(i) - timePrev, tau, sigmaDMC);
        
        % --- Time Update
        pMinus = phi * pPrev * phi' + Q;
        
    else
        % Add SNC insteaf of DMC
        if tVec(i) - obTimePrev < 15
            % Add SNC because gap is small enough
            
            % State Noise Compensation
            GammaQGamma = Dynamics.StateNoiseComp(tVec(i) - timePrev, Q, refTrajStates, Qframe);
        else
            % gap too big, don't add SNC
            GammaQGamma = zeros(NumStates, NumStates);
        end
        
        % --- Time Update ---
        pMinus = phi * pPrev * phi' + GammaQGamma;
        
    end
    
    % --- Computed Measurements --- 
    % each station state
    for j = 1:3
        stationPosECI{i,j} = stationECI{i,j}(1:3);
        stationVelECI{i,j} = stationECI{i,j}(4:6);
    end
    
    % determine station visibility
    [visibilityMask, viewingAngles, ~, obTime] = Measurements.VisibilityMask(stationPosECI(i,:), refPos', 10, tVec(i), 10);
    
    if ~isempty(obTime)
        % observation was made by a station
        statNumOb = find(visibilityMask == 1);
        
        if length(statNumOb) > 1
            % if multiple obs then just take one (just for now)
            multiObs = length(statNumOb);
            warning('multiple observations at same time!')
        else
            % Nothing
            multiObs = 1; 
        end
        
        for q = 1%:multiObs % ignore multiple obs for now
            
        % Computed measurement from the station!
        rangeMeasComp     = refPos - cell2mat(stationPosECI(i,statNumOb(q)))';
        rangeNormComp     = norm(rangeMeasComp);
        rangeRateComp     = dot(rangeMeasComp, refVel - cell2mat(stationVelECI(i,statNumOb(q)))') / rangeNormComp;
        
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
            
            % measurements from all stations at time
            fullObsMeas = observedMeas;
            
            % range and range rate at time
            ObsMeasRange     = fullObsMeas(1:3);
            ObsMeasRangerate = fullObsMeas(4:6);
            
            measDelta = [ObsMeasRange(statNumOb(q)); ObsMeasRangerate(statNumOb(q))] - computedMeas;
            
        end

        if measDelta(1) > 10
        %    measDelta = [0;0]; 
        end
        
        % Compute Htilde
        Htilde{i} = Measurements.HtildeSC(refTrajStates', stationECI{i,statNumOb(q)}, MeasFlag);
        
        if DMC == 1
            % insert zeros to pad Htilde with DMC
            Htilde{i} = [Htilde{i}, zeros(2,3)];
        end
            
        % Kalman Gain
        Kk = pMinus*Htilde{i}' / (Htilde{i}*pMinus*Htilde{i}' + R);
        
        xhat = Kk * measDelta;
        
        % reference state + deviation
        Xplus = refTrajStates' + xhat;
        
        % Save off post fit residual - not needed for Filter computations
        % though
        post_fit_meas(1:2,i) = measDelta - Htilde{i} * xhat;
        
        % covariance update
        Pplus = (eye(NumStates,NumStates) - Kk*Htilde{i}) * pMinus * (eye(NumStates,NumStates) - Kk*Htilde{i})' + Kk*R*Kk';
        
        assert(trace(Pplus) < trace(pMinus), 'Covariance not decrease!');
        
        % set previous observation time as last time an ob occured
        obTimePrev = obTime;
        
        end
        
    else
        % no station made an observation - simply propagate state w/o meas
        Xplus = refTrajStates';
        Pplus = pMinus;
        measDelta = [NaN; NaN];
        post_fit_meas(1:2,i) = [NaN; NaN];
        
    end
    
    
    % save histories
    XrefHist(:,i) = Xplus;
    PplusHist{i} = Pplus; 
    measDeltaHist(:,i) = measDelta; 
    
    % reset everything for next iteration
    XrefPrev = [Xplus; reshape(eye(NumStates,NumStates), [NumStates^2,1])];
    pPrev    = Pplus;
    timePrev = tVec(i);
    
end

% print out final RMS values

% pre-fit residual RMS
measDeltaRMS_Rho    = rmmissing(measDeltaHist(1,:));
measDeltaRMS_RhoDot = rmmissing(measDeltaHist(2,:));

RMS_pre_Rho    = sqrt(sum(measDeltaRMS_Rho.^2)) / length(measDeltaRMS_Rho);
RMS_pre_RhoDot = sqrt(sum(measDeltaRMS_RhoDot.^2)) / length(measDeltaRMS_RhoDot);

% post-fit residual RMS 
measResRMS_Rho    = rmmissing(post_fit_meas(1,:));
measResRMS_RhoDot = rmmissing(post_fit_meas(2,:));

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
