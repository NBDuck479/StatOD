function [Q] = J2_DMC_Q_Matrix(deltaT, tau, sigmaDMC)
%
% deltaT:    [1 x 1] time interval 
% tau:       [3 x 1] DMC time constant 
% sigmaDMC:  [3 x 1] uncertainty in unmodeled acceleration
%

% get Beta from tau
Beta = [inv(tau(1)), inv(tau(2)), inv(tau(3))];

%%%%%%% Construct the Q matrix %%%%%%%%% 
% position component wrt itself
Qrr = @(sigmaDMC, Beta, deltaT) sigmaDMC^2 * ( 1/(3*Beta^2)*deltaT^3 - 1/(Beta^3)*deltaT^2 + 1/(Beta^4)*deltaT - 2/(Beta^4)*exp(-Beta*deltaT)*deltaT + 1/(2*Beta^5)*(1 -exp(-2*Beta*deltaT)));

% position wrt same component velocity
Qrv = @(sigmaDMC, Beta, deltaT) sigmaDMC^2 * (1/(2*Beta^2)*deltaT^2 - 1/(Beta^3)*deltaT + 1/(Beta^3)*exp(-Beta*deltaT)*deltaT + 1/(Beta^4)*(1-exp(-Beta*deltaT)) - 1/(2*Beta^4)*(1-exp(-2*Beta*deltaT)));

Qvr = Qrv; 

% position wrt same error state 
Qrw = @(sigmaDMC, Beta, deltaT) sigmaDMC^2 * (1/(2*Beta^3)*(1-exp(-2*Beta*deltaT)) - 1/(Beta^2)*exp(-Beta*deltaT)*deltaT);

Qwr = Qrw; 

% velocity wrt same velocity state
Qvv = @(sigmaDMC, Beta, deltaT) sigmaDMC^2 * (1/(Beta^2)*deltaT - 2/(Beta^3)*(1-exp(-Beta*deltaT)) + 1/(2*Beta^3)*(1-exp(-2*Beta*deltaT)));

% velocity wrt same error state
Qvw = @(sigmaDMC, Beta, deltaT) sigmaDMC^2 * (1/(2*Beta^2)*(1+exp(-2*Beta*deltaT)) - 1/(Beta^2)*exp(-Beta*deltaT));

Qwv = Qvw; 

% error state wrt same error component error state
Qww = @(sigmaDMC, Beta, deltaT) sigmaDMC^2/(2*Beta) * (1 - exp(-2*Beta*deltaT)); 


% --- Initialize Q as 9 x 9
Q = zeros(9,9); 

% position x component
Q(1,1) = Qrr(sigmaDMC(1), Beta(1), deltaT);

Q(1,4) = Qrv(sigmaDMC(1), Beta(1), deltaT);

Q(1,7) = Qrw(sigmaDMC(1), Beta(1), deltaT); 

% position y component
Q(2,2) = Qrr(sigmaDMC(2), Beta(2), deltaT);

Q(2,5) = Qrv(sigmaDMC(2), Beta(2), deltaT);

Q(2,8) = Qrw(sigmaDMC(2), Beta(2), deltaT);

% position z component
Q(3,3) = Qrr(sigmaDMC(3), Beta(3), deltaT);

Q(3,6) = Qrv(sigmaDMC(3), Beta(3), deltaT);

Q(3,9) = Qrw(sigmaDMC(3), Beta(3), deltaT);

% velocity x component
Q(4,1) = Qvr(sigmaDMC(1), Beta(1), deltaT); 

Q(4,4) = Qvv(sigmaDMC(1), Beta(1), deltaT);

Q(4,7) = Qvw(sigmaDMC(1), Beta(1), deltaT);

% velocity y component
Q(5,2) = Qvr(sigmaDMC(2), Beta(2), deltaT);

Q(5,5) = Qvv(sigmaDMC(2), Beta(2), deltaT);

Q(5,8) = Qvw(sigmaDMC(2), Beta(2), deltaT);

% velocity z component
Q(6,3) = Qvr(sigmaDMC(3), Beta(3), deltaT);

Q(6,6) = Qvv(sigmaDMC(3), Beta(3), deltaT);

Q(6,9) = Qvw(sigmaDMC(3), Beta(3), deltaT);

% accel error x component
Q(7,1) = Qwr(sigmaDMC(1), Beta(1), deltaT);

Q(7,4) = Qwv(sigmaDMC(1), Beta(1), deltaT);

Q(7,7) = Qww(sigmaDMC(1), Beta(1), deltaT);

% accel error y component
Q(8,2) = Qwr(sigmaDMC(2), Beta(2), deltaT);

Q(8,5) = Qwv(sigmaDMC(2), Beta(2), deltaT);

Q(8,8) = Qww(sigmaDMC(2), Beta(2), deltaT);

% accel error z component
Q(9,3) = Qwr(sigmaDMC(3), Beta(3), deltaT);

Q(9,6) = Qwv(sigmaDMC(3), Beta(3), deltaT);

Q(9,9) = Qww(sigmaDMC(3), Beta(3), deltaT);