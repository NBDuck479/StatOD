function [HtildeSC] = HtildeSC(scState, statState)
% computes the linearized H matrix 
% scState - spacecraft state
% statState - station state
% use km and km/s please!

% break states into components
R = scState(1:3);
V = scState(4:6);

Rs = statState(1:3);
Vs = statState(4:6);

% compute range and range rate
rho    = norm(R - Rs); 
rhoDot = dot((R - Rs),(V - Vs)) / rho;

% construct Htilde matrix
HtildeSC = [(R-Rs)'/rho, zeros(1,3); ((V-Vs)'*rho-rhoDot*(R-Rs)') / (rho^2), (R-Rs)'/rho];
end