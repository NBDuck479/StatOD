function [ydot] = Numeric2bodyProp(t, Y, mu, Re)

% velocity maps to itself
ydot(1:3,1) = Y(4:6,1);

% assign states
x = Y(1);
y = Y(2);
z = Y(3);

% compute r
r = sqrt(x^2 + y^2 + z^2);

% accelration partials
apertx = -x*mu / (r^3);
aperty = -y*mu / (r^3);
apertz = -z*mu / (r^3);

% How the accelerations map to output
ydot(4) = apertx;
ydot(5) = aperty;
ydot(6) = apertz;

% Construct A for STM propagation
A = zeros(6,6);

A(1:3,4:6) = eye(3,3);

A(4,1) = -(mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);

A(4,2) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);

A(4,3) = (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);

A(5,1) = (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2);

A(5,2) = -(mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2);

A(5,3) = (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);

A(6,1) = (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2);

A(6,2) = (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2);

A(6,3) = -(mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2);


% Phi is the end of the Y column vector
phiCol = Y(6+1:end);

% reshape to be matrix
phi = reshape(phiCol, [6,6]);

% STM propagation
phiDot = A * phi;

% The state is the first part and phi is the second
ydot(6+1:length(Y)) = reshape(phiDot, [6^2, 1]);
end
