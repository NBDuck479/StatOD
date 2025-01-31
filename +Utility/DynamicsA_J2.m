function [Amat] = NumericProp_J2(state)
% This function produces the A matrix (a jacobian of the accelerations)
% experienced by a spacecraft in orbit

% start with the gravity potential
% declare symbolic variables for diff
syms x y z vx vy vz

a = 6378;
J2 = 0.00108248;
mu = 398600.4415;

r = sqrt(x^2 + y^2 + z^2);

U = mu/r - mu*a^2*J2*(3*z^2-r^2)/(2*r^5);

% take partial of U to get acceleration each component
accelxyz = jacobian(U, [x y z]);

% take the partial of the acceleration to build up STM
stateVec = [x y z vx vy vz];
Aaccel = jacobian(accelxyz, stateVec);

% final A matrix
Amat = [zeros(3,3), eye(3,3), zeros(3,3); ...
    double(subs(Aaccel,stateVec,state)) ; ...
    zeros(2,9)];

end