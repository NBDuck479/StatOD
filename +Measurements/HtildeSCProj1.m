function [Htilde] = HtildeSCProj1(entireState, statState, statNo, MeasFlag)
% computes the linearized H matrix the mu, J2, Cd as added states
%
%%%%%%%%%%% Inputs %%%%%%%%%%
% scState       [6 x 1] position and velocity of spacecraft
% statState     [6 x 3] position and velocity (row) of each station (col) 
% statNo        [1 x 1] what station made the observation
% 
%%%%%%%%%%% Outputs %%%%%%%%%
% Htilde        [2 x 18] range and range rate for each station 
% 
%%%% Units: km, km/s, deg

% Spacecraft States
x = entireState(1);
y = entireState(2);
z = entireState(3);   
vx = entireState(4);
vy = entireState(5);
vz = entireState(6);

% station states
% - station1
if statNo == 1
    xS1  =  statState(1);
    yS1  =  statState(2);
    zS1  =  statState(3);
    vxS1 =  statState(4);
    vyS1 =  statState(5);
    vzS1 =  statState(6);
    
elseif statNo == 2
    % - station2
    xS2  =  statState(1);
    yS2  =  statState(2);
    zS2  =  statState(3);
    vxS2 =  statState(4);
    vyS2 =  statState(5);
    vzS2 =  statState(6);
    
elseif statNo == 3
    % - station3
    xS3  =  statState(1);
    yS3  =  statState(2);
    zS3  =  statState(3);
    vxS3 =  statState(4);
    vyS3 =  statState(5);
    vzS3 =  statState(6);
end

%% Process to Compute Htilde
    % syms x y z vx vy vz mu J2 Cd xS1 yS1 zS1 xS2 yS2 zS2 xS3 yS3 zS3 vxS3 vyS3 vzS3 

    % measurement rho
    % rho = norm([x; y; z] - [xS3; yS3; zS3]);

    % measurement range rate
    % rhoDot = dot([x; y; z] - [xS3; yS3; zS3], [vx; vy; vz] - [vxS3; vyS3; vzS3]) / rho;

    % construct state vector [pos vel mu J2 Cd Stat1 Stat2 Stat3]
    % stateVec = [x y z vx vy vz mu J2 Cd xS1 yS1 zS1 xS2 yS2 zS2 xS3 yS3 zS3];


    % take partials of rho WRT state 
    % rhoWRTstate = jacobian(rho, stateVec);

    % take partials of rhoDot WRT state
    % rhoDotWRTstate = jacobian(rhoDot, stateVec);



%% --- Construct Htilde matrix
% Going to have gates based on what station made obs

if statNo == 1
    % have only partials WRT to station 101
    
    % - partials of rho WRT entire state
    Htilde(1,1) = (abs(x - xS1)*sign(x - xS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(1,2) = (abs(y - yS1)*sign(y - yS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(1,3) = (abs(z - zS1)*sign(z - zS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(1,4:6) = zeros(1,3); % partials are zero wrt velocity states
    
    Htilde(1,7:9) = zeros(1,3); % partials are zero wrt grav & Cd
    
    Htilde(1,10) = -(abs(x - xS1)*sign(x - xS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(1,11) = -(abs(y - yS1)*sign(y - yS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(1,12) = -(abs(z - zS1)*sign(z - zS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);

    Htilde(1,13:18) = zeros(1,6); % partials are zero wrt other stations

    
    % - partials of rhoDot WRT entire state
    Htilde(2,1) = (vx - vxS1)/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2) - (abs(x - xS1)*sign(x - xS1)*((conj(x) - conj(xS1))*(vx - vxS1) + (conj(y) - conj(yS1))*(vy - vyS1) + (conj(z) - conj(zS1))*(vz - vzS1)))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(3/2);

    Htilde(2,2) = (vy - vyS1)/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2) - (abs(y - yS1)*sign(y - yS1)*((conj(x) - conj(xS1))*(vx - vxS1) + (conj(y) - conj(yS1))*(vy - vyS1) + (conj(z) - conj(zS1))*(vz - vzS1)))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(3/2);

    Htilde(2,3) = (vz - vzS1)/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2) - (abs(z - zS1)*sign(z - zS1)*((conj(x) - conj(xS1))*(vx - vxS1) + (conj(y) - conj(yS1))*(vy - vyS1) + (conj(z) - conj(zS1))*(vz - vzS1)))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(3/2);

    Htilde(2,4) = (conj(x) - conj(xS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(2,5) = (conj(y) - conj(yS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(2,6) = (conj(z) - conj(zS1))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(2,7:9) = zeros(1,3); % partials are zero wrt grav & Cd
    
    Htilde(2,10) = (abs(x - xS1)*sign(x - xS1)*((conj(x) - conj(xS1))*(vx - vxS1) + (conj(y) - conj(yS1))*(vy - vyS1) + (conj(z) - conj(zS1))*(vz - vzS1)))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(3/2) - (vx - vxS1)/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);

    Htilde(2,11) = (abs(y - yS1)*sign(y - yS1)*((conj(x) - conj(xS1))*(vx - vxS1) + (conj(y) - conj(yS1))*(vy - vyS1) + (conj(z) - conj(zS1))*(vz - vzS1)))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(3/2) - (vy - vyS1)/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);

    Htilde(2,12) = (abs(z - zS1)*sign(z - zS1)*((conj(x) - conj(xS1))*(vx - vxS1) + (conj(y) - conj(yS1))*(vy - vyS1) + (conj(z) - conj(zS1))*(vz - vzS1)))/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(3/2) - (vz - vzS1)/(abs(x - xS1)^2 + abs(y - yS1)^2 + abs(z - zS1)^2)^(1/2);
    
    Htilde(2,13:18) = zeros(1,6); % partials are zero wrt other stations

    
elseif statNo == 2
    % have only partials WRT to station 101
    
    % - partials of rho WRT entire state
    Htilde(1,1) = (abs(x - xS2)*sign(x - xS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);
    
    Htilde(1,2) = (abs(y - yS2)*sign(y - yS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);
    
    Htilde(1,3) = (abs(z - zS2)*sign(z - zS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);
    
    Htilde(1,4:6) = zeros(1,3); % zeros wrt velocity
    
    Htilde(1,7:9) = zeros(1,3); % zero wrt grav & Cd
    
    Htilde(1,10:12) = zeros(1,3); % zeros wrt station1
    
    Htilde(1,13) = -(abs(x - xS2)*sign(x - xS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);
    
    Htilde(1,14) = -(abs(y - yS2)*sign(y - yS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);
    
    Htilde(1,15) = -(abs(z - zS2)*sign(z - zS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);
    
    Htilde(1,16:18) = zeros(1,3); % zeros wrt station 3
    
    % - partials of rhoDot WRT entire state
    Htilde(2,1) = (vx - vxS2)/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2) - (abs(x - xS2)*sign(x - xS2)*((conj(x) - conj(xS2))*(vx - vxS2) + (conj(y) - conj(yS2))*(vy - vyS2) + (conj(z) - conj(zS2))*(vz - vzS2)))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(3/2);
    
    Htilde(2,2) = (vy - vyS2)/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2) - (abs(y - yS2)*sign(y - yS2)*((conj(x) - conj(xS2))*(vx - vxS2) + (conj(y) - conj(yS2))*(vy - vyS2) + (conj(z) - conj(zS2))*(vz - vzS2)))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(3/2);

    Htilde(2,3) = (vz - vzS2)/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2) - (abs(z - zS2)*sign(z - zS2)*((conj(x) - conj(xS2))*(vx - vxS2) + (conj(y) - conj(yS2))*(vy - vyS2) + (conj(z) - conj(zS2))*(vz - vzS2)))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(3/2);

    Htilde(2,4) = (conj(x) - conj(xS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);

    Htilde(2,5) = (conj(y) - conj(yS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);

    Htilde(2,6) = (conj(z) - conj(zS2))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);

    Htilde(2,7:9) = zeros(1,3); %zero WRT grav & Cd
    
    Htilde(2,10:12) = zeros(1,3); % zeros WRT station 1
    
    Htilde(2,13) = (abs(x - xS2)*sign(x - xS2)*((conj(x) - conj(xS2))*(vx - vxS2) + (conj(y) - conj(yS2))*(vy - vyS2) + (conj(z) - conj(zS2))*(vz - vzS2)))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(3/2) - (vx - vxS2)/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);

    Htilde(2,14) = (abs(y - yS2)*sign(y - yS2)*((conj(x) - conj(xS2))*(vx - vxS2) + (conj(y) - conj(yS2))*(vy - vyS2) + (conj(z) - conj(zS2))*(vz - vzS2)))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(3/2) - (vy - vyS2)/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);

    Htilde(2,15) = (abs(z - zS2)*sign(z - zS2)*((conj(x) - conj(xS2))*(vx - vxS2) + (conj(y) - conj(yS2))*(vy - vyS2) + (conj(z) - conj(zS2))*(vz - vzS2)))/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(3/2) - (vz - vzS2)/(abs(x - xS2)^2 + abs(y - yS2)^2 + abs(z - zS2)^2)^(1/2);

    Htilde(2,16:18) = zeros(1,3); % zeros WRT station 3
    
    
elseif statNo == 3
    % have only partials WRT to station 101
    
    % - partials of rho WRT entire state
    Htilde(1,1) = (abs(x - xS3)*sign(x - xS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(1,2) = (abs(y - yS3)*sign(y - yS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(1,3) = (abs(z - zS3)*sign(z - zS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(1,4:6) = zeros(1,3); % zeros wrt velo
    
    Htilde(1,7:9) = zeros(1,3); % zeros wrt grav & Cd
    
    Htilde(1,10:12) = zeros(1,3); % zeros wrt station 1
    
    Htilde(1,13:15) = zeros(1,3); % zeros wrt staiotn 2
    
    Htilde(1,16) = -(abs(x - xS3)*sign(x - xS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);

    Htilde(1,17) = -(abs(y - yS3)*sign(y - yS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);

    Htilde(1,18) = -(abs(z - zS3)*sign(z - zS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);

    
    % - partials of RhoDot WRT entire state
    Htilde(2,1) = (vx - vxS3)/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2) - (abs(x - xS3)*sign(x - xS3)*((conj(x) - conj(xS3))*(vx - vxS3) + (conj(y) - conj(yS3))*(vy - vyS3) + (conj(z) - conj(zS3))*(vz - vzS3)))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(3/2);

    Htilde(2,2) = (vy - vyS3)/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2) - (abs(y - yS3)*sign(y - yS3)*((conj(x) - conj(xS3))*(vx - vxS3) + (conj(y) - conj(yS3))*(vy - vyS3) + (conj(z) - conj(zS3))*(vz - vzS3)))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(3/2);
    
    Htilde(2,3) = (vz - vzS3)/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2) - (abs(z - zS3)*sign(z - zS3)*((conj(x) - conj(xS3))*(vx - vxS3) + (conj(y) - conj(yS3))*(vy - vyS3) + (conj(z) - conj(zS3))*(vz - vzS3)))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(3/2);

    Htilde(2,4) = (conj(x) - conj(xS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(2,5) = (conj(y) - conj(yS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(2,6) = (conj(z) - conj(zS3))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(2,7:9) = zeros(1,3); % zeros wrt grav & Cd
    
    Htilde(2,10:12) = zeros(1,3); % zeros wrt station 1
    
    Htilde(2,13:15) = zeros(1,3); % zeros wrt staion 2
    
    Htilde(2,16) = (abs(x - xS3)*sign(x - xS3)*((conj(x) - conj(xS3))*(vx - vxS3) + (conj(y) - conj(yS3))*(vy - vyS3) + (conj(z) - conj(zS3))*(vz - vzS3)))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(3/2) - (vx - vxS3)/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);

    Htilde(2,17) = (abs(y - yS3)*sign(y - yS3)*((conj(x) - conj(xS3))*(vx - vxS3) + (conj(y) - conj(yS3))*(vy - vyS3) + (conj(z) - conj(zS3))*(vz - vzS3)))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(3/2) - (vy - vyS3)/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    
    Htilde(2,18) = (abs(z - zS3)*sign(z - zS3)*((conj(x) - conj(xS3))*(vx - vxS3) + (conj(y) - conj(yS3))*(vy - vyS3) + (conj(z) - conj(zS3))*(vz - vzS3)))/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(3/2) - (vz - vzS3)/(abs(x - xS3)^2 + abs(y - yS3)^2 + abs(z - zS3)^2)^(1/2);
    

end


% If theres a measurement flag then just save wanted measurements
if MeasFlag == 1
    Htilde = Htilde(1,:);
elseif MeasFlag == 2
    Htilde = Htilde(2,:);
elseif MeasFlag == 3
    % no need to change
end


end