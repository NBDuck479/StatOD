function BatchResidPlotter(preFit_res, resid_pfHist, obsHist, fig)
    % BatchResidPlotter - Plots pre-fit and post-fit residuals for a given dataset
    %
    % Inputs:
    %   preFit_res    - Pre-fit residuals (2 x N matrix, where N is the number of observations)
    %   resid_pfHist  - Post-fit residuals (2 x N matrix)
    %   obsHist       - Observation history structure (with .time field)
    %   fig           - Figure number (used to plot)

    % Create hours array based on the observation times
    [overlapTimeIdn, ~] = ismember(0:1:18340, obsHist.time);
    hours = [0:1:18340] / 60 / 60;

    % Initialize arrays for residuals
    fillpreFit_res_rho = NaN(1, 18340);
    fillpreFit_res_rhoDot = NaN(1, 18340);
    fillpostFit_res_rho = NaN(1, 18340);
    fillpostFit_res_rhoDot = NaN(1, 18340);

    % Fill residual arrays
    fillpreFit_res_rho(overlapTimeIdn) = preFit_res(1, :);
    fillpreFit_res_rhoDot(overlapTimeIdn) = preFit_res(2, :);
    fillpostFit_res_rho(overlapTimeIdn) = resid_pfHist(1, :);
    fillpostFit_res_rhoDot(overlapTimeIdn) = resid_pfHist(2, :);

    % Plot Pre-Fit residuals
    figure(fig)
    subplot(2,1,1)
    plot(hours, fillpreFit_res_rho', 'o')
    ylabel('Range [km]')
    xlabel('Time [hours]')
    grid on

    subplot(2,1,2)
    plot(hours, fillpreFit_res_rhoDot', 'o')
    ylabel('Range Rate [km/sec]')
    xlabel('Time [hours]')
    sgtitle('Pre-Fit Batch Residuals')
    grid on

    fig = fig + 1;

    % Plot Post-Fit residuals
    figure(fig)
    subplot(2,1,1)
    plot(hours, fillpostFit_res_rho', 'o')
    ylabel('Range [km]')
    xlabel('Time [hours]')
    grid on

    subplot(2,1,2)
    plot(hours, fillpostFit_res_rhoDot', 'o')
    ylabel('Range Rate [km/sec]')
    xlabel('Time [hours]')
    sgtitle('Post-Fit Batch Residuals')
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
    
    sgtitle('Batch Pre-Fit Residuals Histogram')
    
     fig = fig + 1;

     
    figure(fig)
    subplot(2,1,1)
    histogram(fillpostFit_res_rho)
    grid on
    ylabel('Count')
    title('Range')
    
    subplot(2,1,2)
    histogram(fillpostFit_res_rhoDot)
    grid on
    ylabel('Count')
    xlabel('Residual')
    title('Range Rate')
    
    sgtitle('Batch Post-Fit Residuals Histogram')
    
     fig = fig + 1;

end
