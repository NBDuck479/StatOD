function [threeSigmaStates, fig] = SRIFPlotter(covPlus, NumStates, xhist, fig, simuTime)
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
    hours = [1:1:length(simuTime)] / 60 / 60;

    % Initialize arrays for residuals
    fillsumTracePos = NaN(1, 18341);
obTime(1) = 1; 

figure(fig)

% plot the trace of the position over time
semilogy(hours, sumTracePos, 'o')
xlabel('Time [hours]')
ylabel('Tr(Pos)')
title('Trace of Spacecraft Position States')
grid on

fig = fig+1;

figure(fig)

% plot the trace of the position over time
semilogy(hours, sumTracevel, 'o')
xlabel('Time [hours]')
ylabel('Tr(Vel)')
title('Trace of Spacecraft Velocity States')
grid on
fig = fig + 1;

% Plotting the 3-sigma estimation error

% Plotting the 3-sigma estimation error
figure(fig)
subplot(3,1,1)
plot(hours, threeSigmaStates{1}, 'r-', 'LineWidth', 2);
hold on
plot(hours, -threeSigmaStates{1}, 'r-', 'LineWidth', 2);
plot(hours, xhist(1,:), '-b', 'LineWidth', 2)
grid on
ylabel('X [km]');

subplot(3,1,2)
plot(hours, threeSigmaStates{2}, 'r-', 'LineWidth', 2);
hold on
plot(hours, -threeSigmaStates{2}, 'r-', 'LineWidth', 2);
plot(hours, xhist(2,:), '-b', 'LineWidth', 2)
grid on
ylabel('Y [km]');

subplot(3,1,3)
plot(hours, threeSigmaStates{3}, 'r-', 'LineWidth', 2);
hold on
plot(hours, -threeSigmaStates{3}, 'r-', 'LineWidth', 2);
plot(hours, xhist(3,:), '-b', 'LineWidth', 2)
grid on
ylabel('Z [km]');
xlabel('Time [hours]')

sgtitle('SRIF Spacecraft Position Estimation Error')

fig = fig + 1;

% Second figure: States 4 to 6
figure(fig)
subplot(3,1,1)
plot(hours, threeSigmaStates{4}, 'r-', 'LineWidth', 2);
hold on
plot(hours, -threeSigmaStates{4}, 'r-', 'LineWidth', 2);
plot(hours, xhist(4,:), '-b', 'LineWidth', 2)
grid on
ylabel('$\dot{x}$ [km/s]', 'Interpreter', 'latex');

subplot(3,1,2)
plot(hours, threeSigmaStates{5}, 'r-', 'LineWidth', 2);
hold on
plot(hours, -threeSigmaStates{5}, 'r-', 'LineWidth', 2);
plot(hours, xhist(5,:), '-b', 'LineWidth', 2)
grid on
ylabel('$\dot{y}$ [km/s]', 'Interpreter', 'latex');

subplot(3,1,3)
plot(hours, threeSigmaStates{6}, 'r-', 'LineWidth', 2);
hold on
plot(hours, -threeSigmaStates{6}, 'r-', 'LineWidth', 2);
plot(hours, xhist(6,:), '-b', 'LineWidth', 2)
grid on
ylabel('$\dot{z}$ [km/s]', 'Interpreter', 'latex');
xlabel('Time [hours]')
fig = fig + 1;

sgtitle('SRIF Spacecraft Velocity Estimation Error')


end