function [xhat0est, P0est, phiHist, resid_pfHist, Htilde, preFit_res] = BatchFilter(Batchinputs)
% Batch filter to process spaceceraft observation measurements
%
%%%%%%%%%% INPUTS %%%%%%%%%%
% IC        Initial full state and STM



% Get all field names in the struct
fields = fieldnames(Batchinputs);

% Loop through each field and access its value dynamically
for i = 1:length(fields)
    % Get the field name
    fieldName = fields{i};
    
    % Access the value of the field using dynamic field referencing
    fieldValue = Batchinputs.(fieldName);
    
    eval([fieldName ' = fieldValue']);
    
end


% Function to help sort measurements - these are noisy measurements of ref
[observedMeas] = Measurements.FilterMeasLoadIn(yHist);

if MeasFlag == 1
    observedMeas = observedMeas(:,1);
    R = R(1,1);
elseif MeasFlag == 2
    observedMeas = observedMeas(:,2);
    R = R(2,2);
else
    % process all measurements
end

% set starting condition for filter
Xprev   = IC(1:NumStates);
STMprev = IC(NumStates+1:end);

% Build first intitial condition
XrefPrev = [Xprev; STMprev];

% iteration counter
j = 0;

% intial start for overall while loop
% xhatMag = 10*critConv;

% overall convergence loop
while j < iterations
    
    % update iteration counter
    j = j + 1;
    
    % set a priori estimates
    deltaR = inv(chol(P0));
    
    delta = deltaR * deltaR';
    
    % solve for initial N
    N = delta * pert;
    
    % Integrate reference trajectory and STM
    % Set integrator options
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    
    % Integrate Trajectory for time step
   % [T, TrajNom] = ode45(@Dynamics.Numeric_J2_Drag_Prop, tVec, XrefPrev, odeOptions, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0);
    [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [tOverall(1):tOverall(end)], XrefPrev, odeOptions, mu, J2, Re);
    
    % save off phiHist
    phiHist = TrajNom(:,NumStates+1:end);
    
    % Process each observation
    for i = 1:length(tOverall)
        
        % Extract the reference trajectory states
        refTrajStates = TrajNom(i,1:NumStates);
        
        % Extract the Integrated STM - STM is by Row!
        phi = TrajNom(i,NumStates+1:end);
        
        STM = reshape(phi, [NumStates,NumStates]);
        
        % keep track of STM mapped back to time zero
        
        % --- Accumulate current observation
        % function to determine which filter observed and calc measurement
        % delta
%        [statNumOb] = Measurements.StationObs(visibilityMask, i);
        
        % Check if station made an observation
        if ismember(tOverall(i), yHist.obTime)
            
            % Find indices where tOverall(i) matches elements in yHist.obTime
            [row, col] = find(yHist.obTime == tOverall(i));
            
            % the row is the index that matches with time
            obInd = row;
            
            % the column is the station number
            statNumOb = col; 
            
            % Compute Htilde
            Htilde{i} = Measurements.HtildeSC(refTrajStates', stationECI{i,statNumOb}, MeasFlag);

            if MeasFlag == 1
                refRangeMeas = yHistRef.Range(obInd,statNumOb);
                compMeas = refRangeMeas;
                 
            elseif MeasFlag == 2
                refRangeRateMeas = yHistRef.RangeRate(obInd,statNumOb);
                compMeas = refRangeRateMeas;
                
            elseif MeasFlag == 3
                % calcualte the computed measurements
                refRangeMeas     = yHistRef.Range(obInd,statNumOb);
                refRangeRateMeas = yHistRef.RangeRate(obInd,statNumOb);
                
                compMeas = [refRangeMeas; refRangeRateMeas];
                
            end
            
            % measurements from all stations at time
            fullObsMeas = observedMeas(i,:)';
            
            % range and range rate at time
            ObsMeasRange     = fullObsMeas(1:3);
            ObsMeasRangerate = fullObsMeas(4:6);
            
            % measurement delta
            measDelta(:,i) = [ObsMeasRange(statNumOb); ObsMeasRangerate(statNumOb)] - compMeas;
            
            % propagate Htilde
            H = Htilde{i} * STM;
            
            % update delta and N
            delta = delta + H'*inv(R)*H;
            N = N + H'*inv(R)*measDelta(:,i);
        else
            % no ob given
            
        end
        
    end
    
    
    % Finished processing!
    
    % --- Solve Normal Equations
    
    % invert the information matrix, delta
    deltaCholR = chol(delta);
    
    deltaCholRinv = inv(deltaCholR);
    
    P0est = deltaCholRinv*deltaCholRinv';
    
    % solve for initial pert condition
    xhat0est = P0est*N;
    
    % magnitude of estimated x
    xhatMag = norm(xhat0est);
    
    % update and re-run the batch again!
    XrefPrev(1:NumStates) = XrefPrev(1:NumStates) + xhat0est;
    
    % Linearized post-fits Only Loop over Observation Times!
    for k = 1:obInd
        % get STM at each step
        STM = reshape(TrajNom(k, NumStates+1:end), [NumStates, NumStates]);
        
        % prop linearized H
        if isempty(Htilde{k}) && measDelta(2,k) < 0.000001
            % no ob during time - Maybe just NaN???
            H = zeros(2, NumStates);
        else
            % linearize as normal with ob
            H = Htilde{k} * STM;
        end
        resid_pf(:,k) = measDelta(:,k) - H*xhat0est;
    end
    
    
    if MeasFlag == 1
        RMS_Rho = sqrt(sum(measDelta(1,:).^2))/length(measDelta(1,:));
        RMS_RhoDot = NaN;
        
        RMS_Rho_pf = sqrt(sum(resid_pf(1,:).^2)) / length(resid_pf(1,:));
        RMS_RhoDot_pf = NaN;
    elseif MeasFlag == 2
        RMS_Rho = NaN;
        RMS_RhoDot = sqrt(sum(measDelta(1,:).^2))/length(measDelta(1,:));
        
        RMS_RhoDot_pf = sqrt(sum(resid_pf(1,:).^2)) / length(resid_pf(1,:));
        RMS_Rho_pf = NaN;
    else
        RMS_Rho = sqrt(sum(measDelta(1,:).^2))/length(measDelta(1,:));
        RMS_RhoDot = sqrt(sum(measDelta(2,:).^2))/length(measDelta(2,:));
        
        RMS_Rho_pf = sqrt(sum(resid_pf(1,:).^2)) / length(resid_pf(1,:));
        RMS_RhoDot_pf = sqrt(sum(resid_pf(2,:).^2)) / length(resid_pf(2,:));
    end
    
    
    % print out what is happening
    fprintf('--- Iteration %d info: ---\n', j );
    fprintf('Pre-fit residual RMS values\n');
    fprintf('Rho    = %g km\n', RMS_Rho );
    fprintf('RhoDot = %g km\n', RMS_RhoDot );
    fprintf('Post-fit residual RMS values\n');
    fprintf('Rho    = %g km\n', RMS_Rho_pf );
    fprintf('RhoDot = %g km\n', RMS_RhoDot_pf );

    fprintf('\n---------------------------------------------\n\n');

    
end % end of obs measurements


% save off some histories
resid_pfHist = resid_pf;

preFit_res = measDelta; 

end
