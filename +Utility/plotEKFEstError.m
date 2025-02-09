function plotEKFEstError(PplusHist, XrefHist, refPos, refVel, fig)
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
    EKFesterr = XrefHist - [refPos(1:14928, :), refVel(1:14928,:)]';

    % Plot position error for EKF
    figure(fig)
    for i = 1:3
        subplot(3, 1, i)
        
        % Plot 3-sigma bounds
        plot(1:length(threeSigmaState{i}), threeSigmaState{i}, 'r--', 1:length(threeSigmaState{i}), -threeSigmaState{i}, 'r--');
        hold on
        
        % Plot estimation error
        plot(EKFesterr(i, :))
        
        % Dynamically set ylim based on the data range
        yMin = min([min(threeSigmaState{i}), min(EKFesterr(i,:))]);
        yMax = max([max(threeSigmaState{i}), max(EKFesterr(i,:))]);
        ylim([yMin, yMax])
        
        grid on
        ylabel(['Position ' char('X' + (i-1)) ' [km]'])
    end
    sgtitle('EKF Estimation Error Position')
    fig = fig + 1;

    % Plot velocity error for EKF
    figure(fig)
    for i = 4:6
        subplot(3, 1, i - 3)
        
        % Plot 3-sigma bounds
        plot(1:length(threeSigmaState{i}), threeSigmaState{i}, 'r--', 1:length(threeSigmaState{i}), -threeSigmaState{i}, 'r--');
        hold on
        
        % Plot estimation error
        plot(EKFesterr(i, :))
        
        % Dynamically set ylim based on the data range
        yMin = min([min(threeSigmaState{i}), min(EKFesterr(i,:))]);
        yMax = max([max(threeSigmaState{i}), max(EKFesterr(i,:))]);
        ylim([yMin, yMax])
        
        grid on
        ylabel(['Velocity ' char('X' + (i-1)) ' [km/s]'])
    end
    sgtitle('EKF Estimation Error Velocity')
    fig = fig + 1;
end