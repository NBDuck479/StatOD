function [fig, UKFesterr] = plotUKFEstError(PplusHist, XrefHist, refPos, refVel, fig, tOverall)
    % Convert PplusHist from cell to matrix
    PplusHistMat = cell2mat(PplusHist);

    % Get each state uncertainty
    stateUncert = cell(6, 1);
    for i = 1:6
        stateUncert{i} = PplusHistMat(i, i:6:end);
    end

    % Calculate 3 sigma bounds for each state
    threeSigmaState = cell(6, 1);
    for i = 1:6
        threeSigmaState{i} = 3 * sqrt(stateUncert{i});
    end

    % Compute estimation error
    UKFesterr = XrefHist(1:6,:) - [refPos(1:length(tOverall)-1, :), refVel(1:length(tOverall)-1,:)]';
   

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
end
