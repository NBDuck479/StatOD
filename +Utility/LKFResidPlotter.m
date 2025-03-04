function LKFResidPlotter(preFit_res, resid_pfHist, fig, simuTime)
    % LKFResidPlotter - Plots pre-fit and post-fit residuals for a given dataset
    %
    % Inputs:
    %   preFit_res    - Pre-fit residuals (2 x N matrix, where N is the number of observations)
    %   resid_pfHist  - Post-fit residuals (2 x N matrix)
    %   obsHist       - Observation history structure (with .time field)
    %   fig           - Figure number (used to plot)

    % Create hours array based on the observation times
  %  [overlapTimeIdn, ~] = ismember(simuTime, obsHist.time);
    hours = [simuTime] / 60 / 60;
    tLength = length(simuTime);

    % convert cell to double
    preFit_res =preFit_res;
    
    % Fill residual arrays
    fillpreFit_res_rho = preFit_res(1, :);
    fillpreFit_res_rhoDot = preFit_res(2, :);
    fillpostFit_res_rho = resid_pfHist(1, :);
    fillpostFit_res_rhoDot = resid_pfHist(2, :);

    % Plot Pre-Fit residuals
    figure(fig)
    subplot(2,1,1)
    plot( hours, fillpreFit_res_rho', 'o')
    ylabel('Range [km]')
    xlabel('Time [hours]')
    grid on

    subplot(2,1,2)
    plot( hours, fillpreFit_res_rhoDot', 'o')
    ylabel('Range Rate [km/sec]')
    xlabel('Time [hours]')
    sgtitle('Pre-Fit LKF Residuals')
    grid on

    fig = fig + 1;

    % Plot Post-Fit residuals
    figure(fig)
    subplot(2,1,1)
    plot( hours, fillpostFit_res_rho', 'o')
    ylabel('Range [km]')
    xlabel('Time [hours]')
    ylim([-1*10^-3, 1*10^-3])
    grid on

    subplot(2,1,2)
    plot( hours, fillpostFit_res_rhoDot', 'o')
    ylabel('Range Rate [km/sec]')
    xlabel('Time [hours]')
    ylim([-1*10^-5 1*10^-5])
    sgtitle('Post-Fit LKF Residuals')
    grid on

    fig = fig + 1;
    
    
    % --- plot residual histograms
    figure(fig)
    subplot(2,1,1)
    histogram(fillpreFit_res_rho)
    grid on
    ylabel('Count')
    title('Range')
    
    subplot(2,1,2)
    histogram(fillpreFit_res_rhoDot)
    grid on
    ylabel('Count')
    xlabel('Residual')
    title('Range Rate')
    
    sgtitle('LKF Pre-Fit Residuals Histogram')
    
     fig = fig + 1;

     
    figure(fig)
    subplot(2,1,1)
    histogram(fillpostFit_res_rho, 'BinWidth', 0.0002)
    xlim([-0.001 0.001])
    grid on
    ylabel('Count')
    title('Range')
    
    subplot(2,1,2)
    histogram(fillpostFit_res_rhoDot, 'BinWidth', 0.000002)
    grid on
    xlim([-0.00001 0.00001])
    ylabel('Count')
    xlabel('Residual')
    title('Range Rate')
    
    sgtitle('LKF Post-Fit Residuals Histogram')
    
     fig = fig + 1;

end
