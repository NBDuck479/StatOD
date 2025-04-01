function [fig, UKFesterr] = plotUKFEstError(PplusHist, XrefHist, refStates, fig, tOverall, NumStates)
    % Convert PplusHist from cell to matrix
    PplusHistMat = cell2mat(PplusHist);

    % Get each state uncertainty
    stateUncert = cell(NumStates, 1);
    for i = 1:NumStates
        stateUncert{i} = PplusHistMat(i, i:NumStates:end);
    end

    % Calculate 3 sigma bounds for each state
    threeSigmaState = cell(NumStates, 1);
    for i = 1:NumStates
        threeSigmaState{i} = 3 * sqrt(stateUncert{i});
    end

    % Compute estimation error
    UKFesterr = XrefHist(1:NumStates,:) - refStates'; %[refPos(1:length(tOverall)-1, :), refVel(1:length(tOverall)-1,:)]';
   

    % Plot position error for UKF
    figure(fig)
    for i = 1:3
        subplot(3, 1, i)
        
        % Plot 3-sigma bounds
        plot(1:length(threeSigmaState{i}), threeSigmaState{i}, 'r-', 1:length(threeSigmaState{i}), -threeSigmaState{i}, 'r-', 'LineWidth', 2);
        hold on
        
        % Plot estimation error
        plot(UKFesterr(i, :), '-b', 'LineWidth', 2)
        
        grid on
        if i == 1
            ylabel('Position X [km]')
        elseif i == 2
            ylabel('Position Y [km]')
        else
            ylabel('Position Z [km]')
        end
    end
    sgtitle('UKF Estimation Error Position')
    fig = fig + 1;

    % Plot velocity error for UKF
    figure(fig)
    for i = 4:6
        subplot(3, 1, i - 3)
        
        % Plot 3-sigma bounds
        plot(1:length(threeSigmaState{i}), threeSigmaState{i}, 'r-', 1:length(threeSigmaState{i}), -threeSigmaState{i}, 'r-', 'LineWidth', 2);
        hold on
        
        % Plot estimation error
        plot(UKFesterr(i, :),'-b', 'LineWidth', 2)
        
        grid on
        if i == 4
            ylabel('Velocity X [km/s]')
        elseif i == 5
            ylabel('Velocity Y [km/s]')
        else
            ylabel('Velocity Z [km/s]')
        end
    end
    sgtitle('UKF Estimation Error Velocity')
    fig = fig + 1;
    
    if NumStates > 6
    % Plot Mu and J2
    subplot(2,1,1)
    % Plot 3-sigma bounds
    plot(1:length(threeSigmaState{7}), threeSigmaState{7}, 'r-', 1:length(threeSigmaState{7}), -threeSigmaState{7}, 'r-', 'LineWidth', 2);
    hold on
    % Plot estimation error
    plot(UKFesterr(7, :),'-b', 'LineWidth', 2)
    grid on 
    ylabel('mu')
    
    subplot(2,1,2)
    plot(1:length(threeSigmaState{8}), threeSigmaState{8}, 'r-', 1:length(threeSigmaState{8}), -threeSigmaState{8}, 'r-', 'LineWidth', 2);
    hold on
    % Plot estimation error
    plot(UKFesterr(8, :),'-b', 'LineWidth', 2)
    grid on
    ylabel('J2')
    
    sgtitle('mu and J2')
    end
end
