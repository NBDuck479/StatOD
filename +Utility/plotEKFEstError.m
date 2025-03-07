function [fig, EKFesterr] = plotEKFEstError(PplusHist, XrefHist, refPos, refVel, fig, tOverall)
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
    EKFesterr = XrefHist(1:6,:) - [refPos(1:length(tOverall), :), refVel(1:length(tOverall),:)]';
    
    % show estimate DMC error states
   accelError1 = XrefHist(7,:);
   accelError2 = XrefHist(8,:);
   accelError3 = XrefHist(9,:);
   
   figure(fig)
   subplot(3,1,1)
   plot(1:length(threeSigmaState{7}), threeSigmaState{7}, 'r--', 1:length(threeSigmaState{7}), -threeSigmaState{7}, 'r--')
   hold on
   plot(EKFesterr(1,:))
   grid on
   ylabel('$w_x$', 'Interpreter', 'latex');
   
   subplot(3,1,2)
   plot(1:length(threeSigmaState{8}), threeSigmaState{8}, 'r--', 1:length(threeSigmaState{8}), -threeSigmaState{8}, 'r--')
   hold on
   plot(EKFesterr(2,:))
   grid on
   ylabel('$w_y$', 'Interpreter', 'latex');
   
   subplot(3,1,3)
   plot(1:length(threeSigmaState{9}), threeSigmaState{9}, 'r--', 1:length(threeSigmaState{9}), -threeSigmaState{9}, 'r--')
   hold on
   plot(EKFesterr(3,:))
   grid on
   ylabel('$w_z$', 'Interpreter', 'latex');
   
   sgtitle('EKF acceleration error states')
   

    % Plot position error for EKF
    figure(fig)
    for i = 1:3
        subplot(3, 1, i)
        
        % Plot 3-sigma bounds
        plot(1:length(threeSigmaState{i}), threeSigmaState{i}, 'r--', 1:length(threeSigmaState{i}), -threeSigmaState{i}, 'r--');
        hold on
        
        % Plot estimation error
        plot(EKFesterr(i, :))
        
        grid on
        if i == 1
            ylabel('Position X [km]')
        elseif i == 2
            ylabel('Position Y [km]')
        else
            ylabel('Position Z [km]')
        end
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
        
        grid on
        if i == 4
            ylabel('Velocity X [km/s]')
        elseif i == 5
            ylabel('Velocity Y [km/s]')
        else
            ylabel('Velocity Z [km/s]')
        end
    end
    sgtitle('EKF Estimation Error Velocity')
    fig = fig + 1;
end
