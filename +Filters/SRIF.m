function [xhat, P, eps, epsNonWhite] = SRIF(SRIFinputs)
% This function implements a square root information filter 
% Similar structure to LKF
% 
% Reformulate information equation
% Guarantees PosDef matrix
% Decreases Condition number of manipulated matrices (good numerics!)
%
%%%%%%%%%%% Inputs %%%%%%%%% 
%
% P0 and x0 to start 
%
%
%%%%%%%%%%% Outputs %%%%%%%%
%
%
%
% --- Dynamically unpack the input struct for SRIF

% Get all field names in the struct
fields = fieldnames(SRIFinputs);

% Loop through each field and access its value dynamically
for i = 1:length(fields)
    % Get the field name
    fieldName = fields{i};
    
    % Access the value of the field using dynamic field referencing
    fieldValue = SRIFinputs.(fieldName);
    
    eval([fieldName ' = fieldValue']);
    
end

% load in observed sensor measurements 
[observedMeas] = Measurements.FilterMeasLoadIn(yHist);

% if process noise then compute square root 
if processNoise == 1
    % chol of Q 
    cholQ = chol(Q);
    Ru = inv(cholQ);
    
end

% --- SRIF algorithm --- 
% Notation
% Prev: Value from previous time step
% Minus: A Priori Value
% Plus: A Posteriori Value

for k = 1:length(tOverall)
    
    if k ==1 
        % initialize for t0 
        xhatPrev = pert; 
        TrajNom = IC'; 
        timePrev = 0; 
        % - Given P0 and x0 to start 
        cholP0 = chol(P0); 
        Rprev = inv(cholP0);
        
        bPrev = Rprev * xhatPrev; 
    else
        % integrate reference trajectory at each time step
        odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
        [T, TrajNom] = ode45(@Dynamics.NumericJ2Prop, [timePrev:tOverall(k)], refState, odeOptions, mu, J2, Re);
    end
    
    % Extract the reference trajectory states
    refTrajStates = TrajNom(end,1:NumStates);
    
    % Reference Traj position and velocity
    refPos = refTrajStates(1:3);
    refVel = refTrajStates(4:6);
    
    % Extract the Integrated STM (maps previous to current time)
    phi = TrajNom(end,NumStates+1:end);
    
    % reshape the STM
    STM = reshape(phi, [NumStates,NumStates]);
    
    % --- Time Update --- 
    RbarMinus = Rprev * inv(STM); 
    xhatMinus = STM * xhatPrev;
        
    if processNoise == 1 && tOverall(k)-timePrev <20
        % ---- Process Noise Time Update
        [GammaQGamma, Gamma] = Dynamics.StateNoiseComp(tOverall(k)-timePrev, Q, refTrajStates, Qframe);
        
        updateMat = [Ru zeros(3,NumStates) zeros(3,1); -RbarMinus*Gamma RbarMinus bPrev];
        updTrans = Transformation.HouseHolder(updateMat);
        
        % pull out the updates
        Rbar = updTrans(4:end, 4:end-1);
        bBar = updTrans(4:end, end);
        
        % extras for process noise
        %         RbarUk = updTrans(1:3, 1:NumStates);
        %         RbarUxk = updTrans(1:3, NumStates+1:2*NumStates+1);
        %         bTildeUk = updTrans(1:3, end);
    else
        
        % Householder for time update
         updateMat = [RbarMinus bPrev];
         updTrans = Transformation.HouseHolder(updateMat);
         
         % Time updated R and b
         Rbar = updTrans(1:NumStates, 1:NumStates);
        bBar = updTrans(1:NumStates, end);

    end
        
    % Determine if time aligns with an observation 
    if ismember(tOverall(k), yHist.obTime)
        
        % Loop over each element of tOverall and find the matching indices
        % Find indices where tOverall(i) matches elements in yHist.obTime
        [row, col] = find(yHist.obTime == tOverall(k));
        
        % the row is the index that matches with time
        obInd = row;
        
        % the column is the station number
        statNumOb = col;
        
        % station state
        stationState = stationECI{k,statNumOb};
        stationPosECI = stationState(1:3);
        stationVelECI = stationState(4:6);
        
        % each column is station
        Htilde{k} = Measurements.HtildeSC(refTrajStates', stationState, MeasFlag);
        
        % White the sensed matrix H 
        HtildeWhite = inv(V) * Htilde{k};
        
        % --- Calculate Computed measurement
        rangeMeasComp     = refPos - stationPosECI';
        rangeNormComp     = norm(rangeMeasComp);
        rangeRateComp     = dot(rangeMeasComp, refVel - stationVelECI') / rangeNormComp;
        
        % stack computed measurements
        computedMeas = [rangeNormComp; rangeRateComp];
        
        % take observation measurements from time 
        ObsMeas = observedMeas(k,:)';
        
        % range and range rate at time for all stations!
        ObsMeasRange     = ObsMeas(1:3);
        ObsMeasRangerate = ObsMeas(4:6);
        
        % whitened measurement
     %   yWhite = inv(V) * [ObsMeasRange(statNumOb); ObsMeasRangerate(statNumOb)];
        
        % Observed - Commputed 
        yWhite = inv(V)*([ObsMeasRange(statNumOb); ObsMeasRangerate(statNumOb)] - computedMeas);
        
        % --- Measurement Update
        dataMatrix = [Rbar  bBar; HtildeWhite yWhite];
        
        % Householder Transform: TA = [Rk bk; 0 e]
        TA = Transformation.HouseHolder(dataMatrix);
        
        % take out sections of transformed data matrix
        Rk = TA(1:NumStates, 1:NumStates);
        bk = TA(1:NumStates, NumStates+1);
        
        % can backwards substitute for xhat state estimate
        xhat(:,k) = Utility.BackwardsSubsitution(Rk, bk);
        
        % post error - should be [0,1]
        eps(:,k) = yWhite - HtildeWhite * xhat(:,k);
        
        epsNonWhite(:,k) = V *(yWhite - HtildeWhite * xhat(:,k));
        
        
        % reset everything to go back for next iteration 
        Rprev = Rk; 
        bPrev = bk; 
        xhatPrev = xhat(:,k);
        % can save off est error covariance too 
        P{k} = inv(Rk) * inv(Rk)';
        
    else
        % simply time update if no observation
        xhatPrev = xhatMinus;
        Rprev = RbarMinus;
        bPrev = bBar;
        P{k} = inv(Rk) * inv(Rk)';
        xhat(:,k) = xhatMinus;
    end
    % Update for next go around
    timePrev = tOverall(k);
    refState = [refTrajStates'; reshape(eye(NumStates,NumStates), [NumStates^2,1])];
end
        
        
end