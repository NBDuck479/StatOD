function [threeSigmaStates, fig] = plotLKFEstError(covPlus, NumStates, xhist, fig, obTime, simuTime)
% Convert the cell array to a matrix
covPlusmat = cell2mat(covPlus);

% Initialize a cell array to store the threeSigmaState values for each state
threeSigmaStates = cell(NumStates, 1);

% Loop through each state to extract the uncertainty and calculate 3-sigma values
for i = 1:NumStates
    % Extract the uncertainty for the i-th state
    stateUncertainty(i,:) = covPlusmat(i, i:NumStates:end);
    
    % Calculate the 3-sigma value (3 * sqrt(uncertainty)) for the whole vector
    threeSigmaStates{i} = 3 * sqrt(stateUncertainty(i,:));
end

% trace of position and velocity covariance
for i = 1:length(stateUncertainty)
    
    sumTracePos(i) = sum(stateUncertainty(1:3,i));
    
    sumTracevel(i) = sum(stateUncertainty(4:6,i));
    
end


% Create hours array based on the observation times
[overlapTimeIdn, ~] = ismember(simuTime, obTime);
hours = simuTime / 60 / 60;

% Initialize arrays for residuals
fillsumTracePos = NaN(1, 18341);
obTime(1) = 1;


figure(fig)

% plot the trace of the position over time
semilogy(hours, sumTracePos, 'o')
xlabel('Time [hours]')
ylabel('Trace of Position')
title('Trace of Spacecraft Position States')
grid on

fig = fig+1;

figure(fig)

% plot the trace of the position over time
semilogy(hours, sumTracevel, 'o')
xlabel('Time [hours]')
ylabel('Trace of Velocity')
title('Trace of Spacecraft Velocity States')
grid on

fig = fig+1;
% Plotting the 3-sigma estimation error
figure(fig)
subplot(3,1,1)
plot(hours, threeSigmaStates{1}, 'r--', ...
    hours, -threeSigmaStates{1}, 'r--');
hold on
plot(hours,xhist(1,:))
grid on
ylabel('X [km]');

subplot(3,1,2)
plot(hours, threeSigmaStates{2}, 'r--', ...
    hours, -threeSigmaStates{2}, 'r--');
hold on
plot(hours, xhist(2,:))
grid on
ylabel('Y [km]');

subplot(3,1,3)
plot(hours, threeSigmaStates{3}, 'r--', ...
    hours, -threeSigmaStates{3}, 'r--');
hold on
plot(hours, xhist(3,:))
grid on
ylabel('Z [km]');
xlabel('Time [hours]')

sgtitle('LKF Spacecraft Position Estimation Error')

fig = fig+1;

% Second figure: States 4 to 6
figure(fig)
subplot(3,1,1)
plot(hours, threeSigmaStates{4}, 'r--', ...
    hours, -threeSigmaStates{4}, 'r--');
hold on
plot(hours, xhist(4,:))
grid on
ylabel('$\dot{x}$ [km/s]', 'Interpreter', 'latex');

subplot(3,1,2)
plot(hours, threeSigmaStates{5}, 'r--', ...
    hours, -threeSigmaStates{5}, 'r--');
hold on
plot(hours, xhist(5,:))
grid on
ylabel('$\dot{y}$ [km/s]', 'Interpreter', 'latex');

subplot(3,1,3)
plot(hours, threeSigmaStates{6}, 'r--', ...
    hours, -threeSigmaStates{6}, 'r--');
hold on
plot(hours, xhist(6,:))
grid on
ylabel('$\dot{z}$ [km/s]', 'Interpreter', 'latex');
xlabel('Time [hours]')
fig = fig+1;

sgtitle('LKF Spacecraft Velocity Estimation Error')

% if length(xhist(:,1)) > 6
    % Third figure: States 7 to 9
    figure(fig);
    subplot(3,1,1)
    plot(hours, threeSigmaStates{7}, 'r--', ...
        hours, -threeSigmaStates{7}, 'r--');
    hold on
    plot(hours,xhist(7,:))
    grid on
    ylabel('$w_x$', 'Interpreter', 'latex');
    
    subplot(3,1,2)
    plot(hours, threeSigmaStates{8}, 'r--', ...
        hours, -threeSigmaStates{8}, 'r--');
    hold on
    plot(hours,xhist(8,:))
    grid on
    ylabel('$w_y$', 'Interpreter', 'latex');
    
    subplot(3,1,3)
    plot(hours, threeSigmaStates{9}, 'r--', ...
        hours, -threeSigmaStates{9}, 'r--');
    hold on
    plot(hours, xhist(9,:))
    grid on
    ylabel('$w_z$', 'Interpreter', 'latex');
    xlabel('Time [hours]')
    
    sgtitle('LKF acceleration error states', 'Interpreter', 'latex')
    
    fig = fig+1;
%     
%     % Fourth figure: States 10 to 12
%     figure(fig);
%     subplot(3,1,1)
%     plot(hours, threeSigmaStates{10}, 'r--', ...
%         hours, -threeSigmaStates{10}, 'r--');
%     hold on
%     plot(hours,xhist(10,:))
%     grid on
%     ylabel('X [km]');
%     
%     subplot(3,1,2)
%     plot(hours, threeSigmaStates{11}, 'r--', ...
%         hours, -threeSigmaStates{11}, 'r--');
%     hold on
%     plot(hours,xhist(11,:))
%     grid on
%     ylabel('Y [km]');
%     
%     subplot(3,1,3)
%     plot(hours, threeSigmaStates{12}, 'r--', ...
%         hours, -threeSigmaStates{12}, 'r--');
%     hold on
%     plot(hours, xhist(12,:))
%     grid on
%     ylabel('Z [km]');
%     xlabel('Time [hours]')
%     
%     sgtitle('Station 101 Position Estimation Error')
%     
%     fig = fig+1;
%     
%     % Fifth figure: States 13 to 15
%     figure(fig)
%     subplot(3,1,1)
%     plot(hours, threeSigmaStates{13}, 'r--', ...
%         hours, -threeSigmaStates{13}, 'r--');
%     hold on
%     plot(hours,xhist(13,:))
%     grid on
%     
%     ylabel('X [km]');
%     
%     subplot(3,1,2)
%     plot(hours, threeSigmaStates{14}, 'r--', ...
%         hours, -threeSigmaStates{14}, 'r--');
%     hold on
%     plot(hours,xhist(14,:))
%     grid on
%     
%     ylabel('Y [km]');
%     
%     subplot(3,1,3)
%     plot(hours, threeSigmaStates{15}, 'r--', ...
%         hours, -threeSigmaStates{15}, 'r--');
%     hold on
%     plot(hours,xhist(15,:))
%     grid on
%     ylabel('Z [km]');
%     xlabel('Time [hours]')
%     
%     sgtitle('Station 337 Position Estimation Error')
%     
%     fig = fig+1;
%     
%     % Sixth figure: States 16 to 18
%     figure(fig)
%     subplot(3,1,1)
%     plot(hours, threeSigmaStates{16}, 'r--', ...
%         hours, -threeSigmaStates{16}, 'r--');
%     hold on
%     plot(hours,xhist(16,:))
%     grid on
%     ylabel('X [km]');
%     
%     subplot(3,1,2)
%     plot(hours, threeSigmaStates{17}, 'r--', ...
%         hours, -threeSigmaStates{17}, 'r--');
%     hold on
%     plot(hours,xhist(17,:))
%     grid on
%     ylabel('Y [km]');
%     
%     subplot(3,1,3)
%     plot(hours, threeSigmaStates{18}, 'r--', ...
%         hours, -threeSigmaStates{18}, 'r--');
%     hold on
%     plot(hours,xhist(18,:))
%     grid on
%     ylabel('Z [km]');
%     xlabel('Time [hours]')
%     
%     sgtitle('Station 394 Position Estimation Error')
    
    
    % irgnore

end