function [Amat] = DynamicsA_J2_J3(state)
% This function produces the A matrix (a jacobian of the accelerations)
% experienced by a spacecraft in orbit

% start with the gravity potential
% declare symbolic variables for diff
syms x y z vx vy vz mu J2 J3

a = 6378;

r = sqrt(x^2 + y^2 + z^2);

U = mu/r - mu*a^2*J2*(3*z^2-r^2)/(2*r^5) - mu*a^3*J3*(5*z^3-3*z*r^2)/(2*r^7);

% take partial of U to get acceleration each component
accelxyz = jacobian(U, [x y z]);

% %--- here the acceleratoin vector can be rotated to align with ECI
% T_XYZ_xyz = [cosd(alpha) -sind(alpha) 0; sind(alpha) cosd(alpha) 0; 0 0 1];
% 
% % ---

% take the partial of the acceleration to build up STM
stateVec = [x y z vx vy vz mu J2 J3];
Aaccel = jacobian(accelxyz, stateVec);

% final A matrix
Amat = [zeros(3,3), eye(3,3), zeros(3,3); ...
    double(subs(Aaccel,stateVec,state)) ; ...
    zeros(3,9)];

end