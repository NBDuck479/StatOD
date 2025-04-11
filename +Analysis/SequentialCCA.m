function [] = SequentialCCA(Pxx0, C, Ctrue, Pcc0)
% This function does a sequential Consider Covariance Analysis 
%
% This is covariance analysis so the state estiamte is not done. Only the
% covariance estimate

mu = Const.OrbConst.muEarth;
J2 = Const.OrbConst.J2;
J3 = Const.OrbConst.J3;
Re = Const.OrbConst.EarthRadius;

% little c 
c = C - Ctrue; 

% Initialize Consider Analysis Filter
for i = 1:length(tOverall)
    
    if i == 1
        % --- Initialize filter 
        % Need spacecraft state 
        [HtildeXk0] = Measurements.HtildeSC(scState, statState, MeasFlag);
        
        % initial kalman gain
        K0 = Pxx0 * HtildeXk0' * inv(HtildeXk0 * Pxx0 * HtildeXk0' + R);
        
        % J3 doesn't have effect on spacecraft or station state
        HtildeCk0 = [0 0];
    
        % Initial sensitivity 
        S = -K0*HtildeCk0; 
        
        % estimated state covariance 
        P0plus = (eye - K0 * HtildeXk0) * P0minus;
        
        % estimated consider covaraince 
        PC0 = P0plus + S * Pcc0 * S';
        
        % --- Propagate to next time step 
        odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
        thetaInt = ode45(@Dynamics.ConsiderAnalysisJ3, [timePrev,tOverall(i)], Y, odeOptions, theta);
        
        % Also need to output phi 
        
        
        % STM 
        
        Pminus = 
    end
    
    %% --- Measurement sensitivity Matrices
    
    % -- partials wrt state
    [HtildeXk] = Measurements.HtildeSC(scState, statState, MeasFlag);
    %         range{i,j}        = spacecraftState(1:3,i) - stationPosECI{i,j};
    %         rangeNorm(i,j)    = norm(range{i,j});
    %         rangeRate{i,j}    = dot(range{i,j}, spacecraftState(4:6,i) - stationVelECI{i,j}) / rangeNorm(i,j);
    
    % -- partials wrt consider parameters
    % J3 doesn't have effect on spacecraft or station state
    HtildeCk = [0 0];
    
    % Integrate Theta 
    odeOptions = odeset('AbsTol',1e-12,'RelTol', 1e-12);
    thetaInt = ode45(@Dynamics.ConsiderAnalysisJ3, [timePrev,tOverall(i)], Y, odeOptions, theta); 
    
    % A priori sensitivity
    Sbar_k = STM * Sprev + thetaInt; 
    
    % DO I GRAB STATE ESTIMATES FROM FILTER TO USE? 
    % ---- Time Update of state Deviation 
    xbar_ck = xbar + Sbar_k * cbar; 
    
end