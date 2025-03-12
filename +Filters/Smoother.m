function [smoothed_xhat, smoothed_P] = Smoother (xhat_k_k, P_apriori, P_aposteriori, phi, PhiTotal, ProcessNoise)
%%%%%%%%%% Inputs %%%%%%%%%%
%
% xhat_k_k      [6 x K time steps] Filter state estiamte at time k
% P_k_k         [NumStates x NumStates] Filiter Covariance at time k
% P_k_kplus     [NumStates x NumStates] Time update of Covariance matrix
% phi           [NumStates x NumStates] STM from k -> k+1
% Gamma         [NumStates x NumStates] maps process noise to covariance
%
%
%     Initially xhat_k_k is xhat_l_l and things are worked backward
%
%     REMEMBER 'plus' is like the PREVIOUS time step - WORKING BACKWARDS
%
%%%%%%%%% Outputs %%%%%%%%%
%
%
if ProcessNoise == 1
    for i = 14929:-1:1
        % pick this length so it goes back all the way to start
        
        % Schur Identity
        cholAprior = chol(P_apriori{i});
        cholApriorInv = inv(cholAprior);
        P_aprioriINV = cholApriorInv * cholApriorInv'; 
        
        
        Sk = P_aposteriori{i} * phi{i}' * P_aprioriINV; %inv(P_apriori{i});
        
        if i == 14929
            % I dont have a smoothed state/covriance to start with 
            xhat_l_kplus = xhat_k_k(:,i);
            P_l_kplus = P_aposteriori{i};
        end
    
        % Smoothed State Estimate
        xhat_l_k = xhat_k_k(:,i) + Sk * (xhat_l_kplus - phi{i} * xhat_k_k(:,i));
        
        if isnan(xhat_l_k(1))
            hejr = 0;
        end
        
        % Smoothed Covariance
        P_l_k = P_aposteriori{i} + Sk * (P_l_kplus - P_apriori{i}) * Sk';
        
        % go back for another step so reset
        P_l_kplus = P_l_k;
        
        xhat_l_kplus = xhat_l_k;
        
    end
    
else
    % no process noise
    % If no process noise then simply prop backwards!
    xhat_l_k = inv(PhiTotal) * xhat_k_k(:,end);
    
    P_l_k = inv(PhiTotal) * P_apriori{end} * inv(PhiTotal)';
end

% final time k smoothed estimates
smoothed_xhat = xhat_l_kplus;
smoothed_P    = P_l_kplus;