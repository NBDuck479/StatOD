function [Aaccel] = J2propJacobian()
% this function is used to pre-compute the jacobian needed for the
% numerical integration using J2.

syms x y z vx vy vz

a = 6378;
J2 = 0.00108248;
mu = Const.OrbConst.muEarth;

r = sqrt(x^2 + y^2 + z^2);

U = mu/r - mu*a^2*J2*(3*z^2-r^2)/(2*r^5);

% take partial of U to get acceleration each component
accelxyz = jacobian(U, [x y z]);

% take the partial of the acceleration to build up STM
stateVec = [x y z vx vy vz];
Aaccel = jacobian(accelxyz, stateVec);
end