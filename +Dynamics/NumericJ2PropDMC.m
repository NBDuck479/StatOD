function [ydot] = NumericJ2PropDMC(t, Y, mu, J2, Re, tau)
% This function can propagate the state and STM with Dynamic Model Compensation due to the effects of 2
% body and J2 gravitaional dynamics. 
%
%%%%%%%% Inputs %%%%%%%%
% t:    [1 x n] time array
% Y:    [9+81 x 1] Current State [r; v; w] + STM [9^2 x 1]
% mu:   [1 x 1] Gravitional parameter
% J2:   [1 x 1] J2 coefficient
% Re:   [1 x 1] Earth Radius [km]
%
%%%%%%%% Output %%%%%%% 
% ydot:  [9+81 x 1] output state and STM 
% Q:     [9 x 9] covariance of process noise
%
%
%-------------- What is happening -------- 
% Add deterministic acceleration to filter state to estimate values [NOT white noise!]
%
%           X    = [r; v; w]
%           Xdot = [v; vdot + w; wdot]
%
% Notice how the acceleration [w] is added to velocity state and integrated
%
%------------------------------------------

% Number of state
stateLength = 9; 

%--- assign states
% position
x = Y(1);
y = Y(2);
z = Y(3);

% velocity
vx = Y(4);
vy = Y(5);
vz = Y(6);

%--- construct A matrix to be evaluated for each new state
A = zeros(9,9);

% velocity states map to themselves
A(1:3, 4:6) = eye(3,3);

% Accelx partials
A(4,1) = -(2*mu*(x^2 + y^2 + z^2)^3 - 6*mu*x^2*(x^2 + y^2 + z^2)^2 + 3*J2*Re^2*mu*(x^2 + y^2 + z^2)^2 + 105*J2*Re^2*mu*x^2*z^2 - 15*J2*Re^2*mu*x^2*(x^2 + y^2 + z^2) - 15*J2*Re^2*mu*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(4,2) = (6*mu*x*y*(x^2 + y^2 + z^2)^2 - 105*J2*Re^2*mu*x*y*z^2 + 15*J2*Re^2*mu*x*y*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(4,3) = (6*mu*x*z*(x^2 + y^2 + z^2)^2 - 105*J2*Re^2*mu*x*z^3 + 45*J2*Re^2*mu*x*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(4,4:6) = zeros(1,3);

% Accely partials
A(5,1) = (6*mu*x*y*(x^2 + y^2 + z^2)^2 - 105*J2*Re^2*mu*x*y*z^2 + 15*J2*Re^2*mu*x*y*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(5,2) = -(2*mu*(x^2 + y^2 + z^2)^3 - 6*mu*y^2*(x^2 + y^2 + z^2)^2 + 3*J2*Re^2*mu*(x^2 + y^2 + z^2)^2 + 105*J2*Re^2*mu*y^2*z^2 - 15*J2*Re^2*mu*y^2*(x^2 + y^2 + z^2) - 15*J2*Re^2*mu*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(5,3) =(6*mu*y*z*(x^2 + y^2 + z^2)^2 - 105*J2*Re^2*mu*y*z^3 + 45*J2*Re^2*mu*y*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(5,4:6) = zeros(1,3);

% Accelz partials
A(6,1) = (6*mu*x*z*(x^2 + y^2 + z^2)^2 - 105*J2*Re^2*mu*x*z^3 + 45*J2*Re^2*mu*x*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(6,2) = (6*mu*y*z*(x^2 + y^2 + z^2)^2 - 105*J2*Re^2*mu*y*z^3 + 45*J2*Re^2*mu*y*z*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(6,3) = -(2*mu*(x^2 + y^2 + z^2)^3 - 6*mu*z^2*(x^2 + y^2 + z^2)^2 + 9*J2*Re^2*mu*(x^2 + y^2 + z^2)^2 + 105*J2*Re^2*mu*z^4 - 90*J2*Re^2*mu*z^2*(x^2 + y^2 + z^2))/(2*(x^2 + y^2 + z^2)^(9/2));

A(6,4:6) = zeros(1,3);


% --- propagate the state
%%%% THIS IS CORRECT??? %%%%%

% % convert tau into noise mapping matrix
 B = blkdiag(inv(tau(1)), inv(tau(2)), inv(tau(3)));

% velocity maps to itself 
ydot(1:3,1) = Y(4:6,1);

% compute r
r = sqrt(x^2 + y^2 + z^2);

% accerlation due to J2 perturbation

apertx = -((mu*x)/(r^3))*(1-J2*(3/2)*(Re/r)^2*(5*(z/r)^2-1));
aperty = -((mu*y)/(r^3))*(1-J2*(3/2)*(Re/r)^2*(5*(z/r)^2-1));
apertz = -((mu*z)/(r^3))*(1-J2*(3/2)*(Re/r)^2*(5*(z/r)^2-3));

% Accelerations plus the acceleration error state!
ydot(4) = apertx + Y(7);
ydot(5) = aperty + Y(8);
ydot(6) = apertz + Y(9);

% propagate the error states in time
wdot = -B * Y(7:9);

ydot(7) = wdot(1);
ydot(8) = wdot(2);
ydot(9) = wdot(3);


% propagate the STM
AprimeSTM = [zeros(3,3), zeros(3,3), zeros(3,3); zeros(3,3), zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3), -B];

% maybe add the Aprime for error states to the non DMC A matrix? 
Astm = AprimeSTM + A;

phidot = Astm * reshape(Y(stateLength+1:length(Y)), [stateLength, stateLength]);

ydot(stateLength+1:length(Y)) = reshape(phidot, [stateLength^2, 1]);


end
