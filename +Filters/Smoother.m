function [smoothed_xhat, smoothed_P] = Smoother (xhat_k_k, P_k_k, P_k_kplus, phi, Gamma, TimeStepsBack)
%%%%%%%%%% Inputs %%%%%%%%%%
%
% xhat_k_k      [NumStates x 1] Filter state estiamte at time k
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

for i = 1:TimeStepsBack
    % Schur Identity
    Sk = P_k_k * phi' * inv(P_k_kplus);
    
    % Smoothed State Estimate
    xhat_l_k = xhat_k_k + Sk * (xhat_l_kplus - phi * xhat_k_k);
    
    % Smoothed Covariance 
    P_l_k = P_k_k + Sk * (P_l_kplus - P_k_kplus) * Sk'; 
    
    % go back another step so reset 
    P_k_kplus = P_l_k; 
    
    xhat_l_kplus = xhat_l_k;
    
end

% final time k smoothed estimates
smoothed_xhat = xhat_l_k;
smoothed_P    = P_l_k;