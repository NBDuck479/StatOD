function [GammaQGamma, Gamma] = StateNoiseComp(deltaT, Q, scState, Qframe)
% Computes the State Noise compenstation for the Filter
% This assumes the state has 6 parameters
%
%%%%% Inputs %%%%%
% deltaT        [1 x 1] time step size
% Q             [3 x 3] uncertainty in ACCELERATION for each direction
% scState       [6 x 1] state of spacecraft in ECI
% Qframe        [ char ] Frame of input Q matrix
%
%%%%% Output %%%%
% GammaQGamma
%
%

% Gamma from tk to tk+1
Gamma = deltaT * [deltaT/2 * eye(3); eye(3)];

if strcmp(Qframe, 'ECI')
    
    % output
    GammaQGamma = Gamma * Q * Gamma';
    
elseif strcmp(Qframe, 'RSW')
    % Rotation form ECI to RSW frame
    Rhat = -scState(1:3) / norm(scState(1:3));
    What = -cross(scState(1:3), scState(4:6)) / norm(cross(scState(1:3), scState(4:6)));
    Shat = cross(What, Rhat);
    
    % Construct frame transformation
    T_ECI2RSW = [Rhat', What', Shat'];
    
    GammaQGamma = Gamma * T_ECI2RSW' * Q * T_ECI2RSW * Gamma'; 
    
else
    warning('Wrong Q frame!!!')
end