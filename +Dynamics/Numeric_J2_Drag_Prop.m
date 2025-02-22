function [ydot] = Numeric_J2_Drag_Prop(t, Y, Re, omegaEarth, Area, Mass, DragH, r0Drag, DragRho0)
% This function can propagate the state and STM due to the effects of 2
% body and J2 gravitaional dynamics with air drag on the spacecraft.
% It will also propagate the whole state provided in project 1
%
%%%%%%%%%%% Inputs %%%%%%%%%%
% t:   [n x 1] time array
% Y:   [342 x 1] Current State and STM
% Re:  Earth Radius [km]
%
% NOTE: Stations input in ECEF and output in ECEF!!!
% Spacecraft state input in ECI and output in ECI!!!
%
%%%%%%%%%% Outputs %%%%%%%%%
% ydot:  [342 x 1] Propagated State and STM

%--- assign states
% position
x = Y(1);
y = Y(2);
z = Y(3);

% velocity
vx = Y(4);
vy = Y(5);
vz = Y(6);

% gravity
mu = Y(7);
J2 = Y(8);

% drag
Cd = Y(9);

% --- Building A matrix derivatives
    % syms x y z vx vy vz mu J2 Cd r0Drag DragH omegaEarth Area Mass DragRho0 Re

% compute r
r = sqrt(x^2 + y^2 + z^2);

% Realtive velocity for drag
Vr = [vx; vy; vz] - cross([0;0;omegaEarth], [x;y;z]);

% --- Atmo Drag Model
rho = DragRho0 * exp(-(r-r0Drag)/DragH);

% --- Acceleration due to Drag
accelDrag = -0.5 * (Cd*Area/Mass) * rho * norm(Vr) * Vr;

% --- Accerlation due to J2 perturbation
aJ2x = -((mu*x)/(r^3))*(1-J2*(3/2)*(Re/r)^2*(5*(z/r)^2-1));
aJ2y = -((mu*y)/(r^3))*(1-J2*(3/2)*(Re/r)^2*(5*(z/r)^2-1));
aJ2z = -((mu*z)/(r^3))*(1-J2*(3/2)*(Re/r)^2*(5*(z/r)^2-3));

% Map velocity states to themselves
ydot(1:3,1) = [vx; vy; vz];

% How the accelerations map to output with J2 plus Drag
accelJ2Drag = [aJ2x + accelDrag(1); aJ2y + accelDrag(2); aJ2z + accelDrag(3)];

ydot(4:6,1) = accelJ2Drag;

% the rate of change of gravity coeff and Cd is zero
ydot(7:9) = zeros(3,1);

% stations assumed in ECEF and therefore do not change
% % Rotation matrix for simple transformation
% Rz = @(Theta) [cosd(Theta), -sind(Theta), 0; sind(Theta), cosd(Theta), 0; 0, 0, 1];
 
 ydot(10:12) = cross([0;0;omegaEarth], Y(10:12));
 
 ydot(13:15) = cross([0;0;omegaEarth], Y(13:15));
 
 ydot(16:18) = cross([0;0;omegaEarth], Y(16:18));

% ydot(10:18) = zeros(9,1);

    % Take jacobian of acceleration with J2 Drag WRT all states
    % stateVec = [x y z vx vy vz mu J2 Cd];
    % accelJ2DragWRTStates = jacobian(accelJ2Drag, stateVec);

% going to propagate the STM too!
%--- construct A matrix to be evaluated for each new state

% Initialize A with zeros to make easier to fill in
A = zeros(18,18);

% velocity states map to themselves
A(1:3, 4:6) = eye(3,3);

% Accelx partials
A(4,1) = (mu*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x^2*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) - (mu*x*((15*J2*Re^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (3*J2*Re^2*x*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) + (Area*Cd*DragRho0*omegaEarth*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vy - omegaEarth*x)*sign(vy - omegaEarth*x)*(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2)) + (Area*Cd*DragRho0*x*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vx + omegaEarth*y)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2));

A(4,2) = (Area*Cd*DragRho0*y*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vx + omegaEarth*y)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2)) - (3*mu*x*y*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) - (Area*Cd*DragRho0*omegaEarth*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass) - (Area*Cd*DragRho0*omegaEarth*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vx + omegaEarth*y)*sign(vx + omegaEarth*y)*(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2)) - (mu*x*((15*J2*Re^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (3*J2*Re^2*y*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2);

A(4,3) = (mu*x*((3*J2*Re^2*((10*z)/(x^2 + y^2 + z^2) - (10*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (3*J2*Re^2*z*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x*z*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) + (Area*Cd*DragRho0*z*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vx + omegaEarth*y)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2));

A(4,4) = - (Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass) - (Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vx + omegaEarth*y)*sign(vx + omegaEarth*y)*(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(4,5) = -(Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vy - omegaEarth*x)*sign(vy - omegaEarth*x)*(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(4,6) = -(Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vz)*sign(vz)*(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(4,7) = (x*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2);

A(4,8) = (3*Re^2*mu*x*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(5/2));

A(4,9) = -(Area*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vx + omegaEarth*y)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass);

% Accely partials
A(5,1) = (Area*Cd*DragRho0*omegaEarth*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass) - (3*mu*x*y*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) - (mu*y*((15*J2*Re^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (3*J2*Re^2*x*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) + (Area*Cd*DragRho0*omegaEarth*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vy - omegaEarth*x)*sign(vy - omegaEarth*x)*(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2)) + (Area*Cd*DragRho0*x*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vy - omegaEarth*x)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2));

A(5,2) = (mu*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*y^2*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) - (mu*y*((15*J2*Re^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (3*J2*Re^2*y*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (Area*Cd*DragRho0*omegaEarth*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vx + omegaEarth*y)*sign(vx + omegaEarth*y)*(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2)) + (Area*Cd*DragRho0*y*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vy - omegaEarth*x)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2));

A(5,3) = (mu*y*((3*J2*Re^2*((10*z)/(x^2 + y^2 + z^2) - (10*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (3*J2*Re^2*z*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*y*z*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) + (Area*Cd*DragRho0*z*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vy - omegaEarth*x)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2));

A(5,4) = -(Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vx + omegaEarth*y)*sign(vx + omegaEarth*y)*(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(5,5) = - (Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass) - (Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vy - omegaEarth*x)*sign(vy - omegaEarth*x)*(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(5,6) = -(Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vz)*sign(vz)*(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(5,7) = (y*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 1)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2);

A(5,8) = (3*Re^2*mu*y*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(5/2));

A(5,9) = -(Area*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(vy - omegaEarth*x)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass);

% Accelz partials
A(6,1) = (Area*Cd*DragRho0*vz*x*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2)) - (3*mu*x*z*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 3)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) - (mu*z*((15*J2*Re^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (3*J2*Re^2*x*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) + (Area*Cd*DragRho0*omegaEarth*vz*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vy - omegaEarth*x)*sign(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(6,2) = (Area*Cd*DragRho0*vz*y*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2)) - (3*mu*y*z*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 3)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) - (mu*z*((15*J2*Re^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (3*J2*Re^2*y*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (Area*Cd*DragRho0*omegaEarth*vz*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vx + omegaEarth*y)*sign(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(6,3) = (mu*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 3)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2) + (mu*z*((3*J2*Re^2*((10*z)/(x^2 + y^2 + z^2) - (10*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (3*J2*Re^2*z*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*z^2*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 3)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2) + (Area*Cd*DragRho0*vz*z*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*DragH*Mass*(x^2 + y^2 + z^2)^(1/2));

A(6,4) = -(Area*Cd*DragRho0*vz*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vx + omegaEarth*y)*sign(vx + omegaEarth*y))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(6,5) = -(Area*Cd*DragRho0*vz*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vy - omegaEarth*x)*sign(vy - omegaEarth*x))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(6,6) = - (Area*Cd*DragRho0*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass) - (Area*Cd*DragRho0*vz*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*abs(vz)*sign(vz))/(2*Mass*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2));

A(6,7) = (z*((3*J2*((5*z^2)/(x^2 + y^2 + z^2) - 3)*Re^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2);

A(6,8) = (3*Re^2*mu*z*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(2*(x^2 + y^2 + z^2)^(5/2));

A(6,9) = -(Area*DragRho0*vz*exp((r0Drag - (x^2 + y^2 + z^2)^(1/2))/DragH)*(abs(vy - omegaEarth*x)^2 + abs(vx + omegaEarth*y)^2 + abs(vz)^2)^(1/2))/(2*Mass);

% mu partials
A(7,1:9) = zeros(1,9);

% J2 partials
A(8,1:9) = zeros(1,9);

% Cd partials
A(9,1:9) = zeros(1,9);

% Partials of the stations
A(10:12,10:12) = [0,-omegaEarth,0;omegaEarth,0, 0; 0, 0, 0]; 

A(13:15,13:15) = [0,-omegaEarth,0;omegaEarth,0, 0; 0, 0, 0];

A(16:18,16:18) = [0,-omegaEarth,0;omegaEarth,0, 0; 0, 0, 0];

% Phi is the end of the Y column vector
phiCol = Y(18+1:end);

% reshape to be matrix
phi = reshape(phiCol, [18,18]);

% STM propagation
phiDot = A * phi;

% The state is the first part and phi is the second
ydot(18+1:length(Y)) = reshape(phiDot, [18^2, 1]);

% end

end
