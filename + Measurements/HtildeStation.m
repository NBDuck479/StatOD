function [HtildeStation] = HtildeStation(scState, statState)
% computes the linearized H matrix - station! 
% scState - spacecraft state
% statState - station state
% use km and km/s please!

% break states into components
R = scState(1:3);
V = scState(4:6);

Rs = statState(1:3);
Vs = statState(4:6);

% compute range
rho    = norm(R - Rs); 
rhoDot = dot((R - Rs),(V - Vs)) / rho;

% construct Htilde matrix - Station
HtildeStation = [(-R+Rs)'/rho; (rho*(Vs-V)'-rhoDot*(Rs-R)')/(rho^2)];
end