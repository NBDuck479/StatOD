function [ydot] = DynamicsA_J2_J3(t, Y, mu, J2, J3, Re)
% velocity maps to itself - state prop
ydot(1:3,1) = Y(4:6,1);

% assign states
x = Y(1);
y = Y(2);
z = Y(3);

% compute r
r = sqrt(x^2 + y^2 + z^2);

% Accel partials
apertx = (J2*Re^2*mu*x)/(x^2 + y^2 + z^2)^(5/2) - (mu*x)/(x^2 + y^2 + z^2)^(3/2) + (3*J3*Re^3*mu*x*z)/(x^2 + y^2 + z^2)^(7/2) - (5*J2*Re^2*mu*x*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)) - (7*J3*Re^3*mu*x*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(9/2));
aperty = (J2*Re^2*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (mu*y)/(x^2 + y^2 + z^2)^(3/2) + (3*J3*Re^3*mu*y*z)/(x^2 + y^2 + z^2)^(7/2) - (5*J2*Re^2*mu*y*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)) - (7*J3*Re^3*mu*y*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(9/2));
apertz = (J3*Re^3*mu*(3*x^2 + 3*y^2 - 6*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)) - (mu*z)/(x^2 + y^2 + z^2)^(3/2) - (2*J2*Re^2*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (5*J2*Re^2*mu*z*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(7/2)) - (7*J3*Re^3*mu*z^2*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(9/2));


% Accel map to output
ydot(4) = apertx;
ydot(5) = aperty;
ydot(6) = apertz;


% % take paertials to get equations
% syms x y z vx vy vz mu J2 J3 Re
% 
% aperts = [apertx, aperty, apertz];
% 
% Apartials = simplify(jacobian(aperts, [x y z vx vy vz]))


% ---- Construct A Dynamics ----- 

% Allocate A with all zeros
A = zeros(6,6);

% velocity states map to themselves
A(1:3,4:6) = eye(3,3);

% partials WRT X 
A(4,1) = (mu*(4*x^8 + 10*x^6*y^2 + 10*x^6*z^2 + 12*J2*Re^2*x^6 + 6*x^4*y^4 + 12*x^4*y^2*z^2 + 21*J2*Re^2*x^4*y^2 + 6*x^4*z^4 - 69*J2*Re^2*x^4*z^2 + 90*J3*Re^3*x^4*z - 2*x^2*y^6 - 6*x^2*y^4*z^2 + 6*J2*Re^2*x^2*y^4 - 6*x^2*y^2*z^4 - 63*J2*Re^2*x^2*y^2*z^2 + 75*J3*Re^3*x^2*y^2*z - 2*x^2*z^6 - 69*J2*Re^2*x^2*z^4 - 205*J3*Re^3*x^2*z^3 - 2*y^8 - 8*y^6*z^2 - 3*J2*Re^2*y^6 - 12*y^4*z^4 + 6*J2*Re^2*y^4*z^2 - 15*J3*Re^3*y^4*z - 8*y^2*z^6 + 21*J2*Re^2*y^2*z^4 + 5*J3*Re^3*y^2*z^3 - 2*z^8 + 12*J2*Re^2*z^6 + 20*J3*Re^3*z^5))/(2*(x^2 + y^2 + z^2)^(11/2));

A(4,2) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) - (10*J2*Re^2*mu*x*y)/(x^2 + y^2 + z^2)^(7/2) + (35*J2*Re^2*mu*x*y*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(9/2)) - (42*J3*Re^3*mu*x*y*z)/(x^2 + y^2 + z^2)^(9/2) + (63*J3*Re^3*mu*x*y*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(11/2));

A(4,3) = (3*mu*x*(2*x^6*z + 6*x^4*y^2*z + 6*x^4*z^3 + 15*J2*Re^2*x^4*z - 5*J3*Re^3*x^4 + 6*x^2*y^4*z + 12*x^2*y^2*z^3 + 30*J2*Re^2*x^2*y^2*z - 10*J3*Re^3*x^2*y^2 + 6*x^2*z^5 - 5*J2*Re^2*x^2*z^3 + 60*J3*Re^3*x^2*z^2 + 2*y^6*z + 6*y^4*z^3 + 15*J2*Re^2*y^4*z - 5*J3*Re^3*y^4 + 6*y^2*z^5 - 5*J2*Re^2*y^2*z^3 + 60*J3*Re^3*y^2*z^2 + 2*z^7 - 20*J2*Re^2*z^5 - 40*J3*Re^3*z^4))/(2*(x^2 + y^2 + z^2)^(11/2));

A(4,4:6) = zeros(1,3);

% partials WRT Y
A(5,1) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) - (10*J2*Re^2*mu*x*y)/(x^2 + y^2 + z^2)^(7/2) + (35*J2*Re^2*mu*x*y*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(9/2)) - (42*J3*Re^3*mu*x*y*z)/(x^2 + y^2 + z^2)^(9/2) + (63*J3*Re^3*mu*x*y*z*(3*x^2 + 3*y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^(11/2));

A(5,2) = (mu*(- 2*x^8 - 2*x^6*y^2 - 8*x^6*z^2 - 3*J2*Re^2*x^6 + 6*x^4*y^4 - 6*x^4*y^2*z^2 + 6*J2*Re^2*x^4*y^2 - 12*x^4*z^4 + 6*J2*Re^2*x^4*z^2 - 15*J3*Re^3*x^4*z + 10*x^2*y^6 + 12*x^2*y^4*z^2 + 21*J2*Re^2*x^2*y^4 - 6*x^2*y^2*z^4 - 63*J2*Re^2*x^2*y^2*z^2 + 75*J3*Re^3*x^2*y^2*z - 8*x^2*z^6 + 21*J2*Re^2*x^2*z^4 + 5*J3*Re^3*x^2*z^3 + 4*y^8 + 10*y^6*z^2 + 12*J2*Re^2*y^6 + 6*y^4*z^4 - 69*J2*Re^2*y^4*z^2 + 90*J3*Re^3*y^4*z - 2*y^2*z^6 - 69*J2*Re^2*y^2*z^4 - 205*J3*Re^3*y^2*z^3 - 2*z^8 + 12*J2*Re^2*z^6 + 20*J3*Re^3*z^5))/(2*(x^2 + y^2 + z^2)^(11/2));

A(5,3) = (3*mu*y*(2*x^6*z + 6*x^4*y^2*z + 6*x^4*z^3 + 15*J2*Re^2*x^4*z - 5*J3*Re^3*x^4 + 6*x^2*y^4*z + 12*x^2*y^2*z^3 + 30*J2*Re^2*x^2*y^2*z - 10*J3*Re^3*x^2*y^2 + 6*x^2*z^5 - 5*J2*Re^2*x^2*z^3 + 60*J3*Re^3*x^2*z^2 + 2*y^6*z + 6*y^4*z^3 + 15*J2*Re^2*y^4*z - 5*J3*Re^3*y^4 + 6*y^2*z^5 - 5*J2*Re^2*y^2*z^3 + 60*J3*Re^3*y^2*z^2 + 2*z^7 - 20*J2*Re^2*z^5 - 40*J3*Re^3*z^4))/(2*(x^2 + y^2 + z^2)^(11/2));

A(5,4:6) = zeros(1,3);

% partials WRT Z
A(6,1) = (3*mu*x*(2*x^6*z + 6*x^4*y^2*z + 6*x^4*z^3 + 15*J2*Re^2*x^4*z - 5*J3*Re^3*x^4 + 6*x^2*y^4*z + 12*x^2*y^2*z^3 + 30*J2*Re^2*x^2*y^2*z - 10*J3*Re^3*x^2*y^2 + 6*x^2*z^5 - 5*J2*Re^2*x^2*z^3 + 60*J3*Re^3*x^2*z^2 + 2*y^6*z + 6*y^4*z^3 + 15*J2*Re^2*y^4*z - 5*J3*Re^3*y^4 + 6*y^2*z^5 - 5*J2*Re^2*y^2*z^3 + 60*J3*Re^3*y^2*z^2 + 2*z^7 - 20*J2*Re^2*z^5 - 40*J3*Re^3*z^4))/(2*(x^2 + y^2 + z^2)^(11/2));

A(6,2) = (3*mu*y*(2*x^6*z + 6*x^4*y^2*z + 6*x^4*z^3 + 15*J2*Re^2*x^4*z - 5*J3*Re^3*x^4 + 6*x^2*y^4*z + 12*x^2*y^2*z^3 + 30*J2*Re^2*x^2*y^2*z - 10*J3*Re^3*x^2*y^2 + 6*x^2*z^5 - 5*J2*Re^2*x^2*z^3 + 60*J3*Re^3*x^2*z^2 + 2*y^6*z + 6*y^4*z^3 + 15*J2*Re^2*y^4*z - 5*J3*Re^3*y^4 + 6*y^2*z^5 - 5*J2*Re^2*y^2*z^3 + 60*J3*Re^3*y^2*z^2 + 2*z^7 - 20*J2*Re^2*z^5 - 40*J3*Re^3*z^4))/(2*(x^2 + y^2 + z^2)^(11/2));

A(6,3) = -(mu*(2*x^8 + 8*x^6*y^2 + 2*x^6*z^2 + 9*J2*Re^2*x^6 + 12*x^4*y^4 + 6*x^4*y^2*z^2 + 27*J2*Re^2*x^4*y^2 - 6*x^4*z^4 - 63*J2*Re^2*x^4*z^2 + 75*J3*Re^3*x^4*z + 8*x^2*y^6 + 6*x^2*y^4*z^2 + 27*J2*Re^2*x^2*y^4 - 12*x^2*y^2*z^4 - 126*J2*Re^2*x^2*y^2*z^2 + 150*J3*Re^3*x^2*y^2*z - 10*x^2*z^6 - 48*J2*Re^2*x^2*z^4 - 200*J3*Re^3*x^2*z^3 + 2*y^8 + 2*y^6*z^2 + 9*J2*Re^2*y^6 - 6*y^4*z^4 - 63*J2*Re^2*y^4*z^2 + 75*J3*Re^3*y^4*z - 10*y^2*z^6 - 48*J2*Re^2*y^2*z^4 - 200*J3*Re^3*y^2*z^3 - 4*z^8 + 24*J2*Re^2*z^6 + 40*J3*Re^3*z^5))/(2*(x^2 + y^2 + z^2)^(11/2));

A(6,4:6) = zeros(1,3); 

% ---- Propagate STM ----

% STM
phiCol = Y(7:end);

% reshapre to be matrix
phi = reshape(phiCol, [6,6]);

% STM propagation
phiDot = A * phi;

% STM is second part of state dot output

ydot(7:length(Y)) = reshape(phiDot, [36, 1]);


end