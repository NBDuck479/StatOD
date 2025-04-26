function CornerPlots(YmonteCat, Ynom, timeIdx, time)
    % Corner Plots Function
    % YmonteCat is expected to be a matrix where each row represents a variable

    % Extract data (every 6th column)
    X = YmonteCat(timeIdx, 1:6:end);
    Y = YmonteCat(timeIdx, 2:6:end);
    Z = YmonteCat(timeIdx, 3:6:end);
    
    VX = YmonteCat(timeIdx, 4:6:end);
    VY = YmonteCat(timeIdx, 5:6:end);
    VZ = YmonteCat(timeIdx, 6:6:end);
    
    % grab nominal states at time 
    nomStates = Ynom(timeIdx, :); 
    
    % Create Position Corner Plots
    figure(1)
    
    % Create subplots
    subplot(3, 3, 1)
    plotHistogram(X, 'X Distribution');
    
    subplot(3, 3, 4)
    plotScatterWithCovariance(X, Y, 'X vs Y', nomStates(1:2));
    
    subplot(3, 3, 5)
    plotHistogram(Y, 'Y Distribution');
    
    subplot(3, 3, 7)
    plotScatterWithCovariance(X, Z, 'X vs Z', nomStates([1,3]));
    
    subplot(3, 3, 8)
    plotScatterWithCovariance(Y, Z, 'Y vs Z', nomStates([2,3]));
    
    subplot(3, 3, 9)
    plotHistogram(Z, 'Z Distribution');
    
    % Overall title for all subplots
sgtitle(sprintf('Position Corner Plots at %.2f hours', time(timeIdx)/(60*60)));
    
    % Now Velocity Corner Plots 
     figure(2);

    % VX Distribution
    subplot(3, 3, 1);
    plotHistogram(VX, 'VX Distribution');

    % VX vs VY
    subplot(3, 3, 4);
    plotScatterWithCovariance(VX, VY, 'VX vs VY', nomStates(4:5));

    % VY Distribution
    subplot(3, 3, 5);
    plotHistogram(VY, 'VY Distribution');

    % VX vs VZ
    subplot(3, 3, 7);
    plotScatterWithCovariance(VX, VZ, 'VX vs VZ', nomStates([4,6]));

    % VY vs VZ
    subplot(3, 3, 8);
    plotScatterWithCovariance(VY, VZ, 'VY vs VZ', nomStates([5,6]));

    % VZ Distribution
    subplot(3, 3, 9);
    plotHistogram(VZ, 'VZ Distribution');

sgtitle(sprintf('Velocity Corner Plots at %.2f hours', time(timeIdx)/(60*60)));

end

function plotHistogram(data, titleText)
    % Helper function to plot histogram
    histogram(data)
    title(titleText);
    mu = mean(data); 
    sigma = std(data); 
    muText = sprintf('\\mu = %.2f', mu); % Format mean value
    sigmaText = sprintf('\\sigma = %.2f', sigma);
    statsText = sprintf('%s\n%s', muText, sigmaText); % Multiline text
    
    % Get current axis limits
    xLimits = xlim;
    yLimits = ylim;
    
    % Set position slightly inset from the top-right corner
    xPos = xLimits(2) - 0.1 * range(xLimits);
    yPos = yLimits(2) - 0.1 * range(yLimits);
    
    % Display the text with a box
    text(xPos, yPos, statsText, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'Margin', 5, ...
        'FontWeight', 'bold');
end

function plotScatterWithCovariance(X, Y, titleText, nomStates)
    % Helper function to plot scatter and covariance ellipse
    plot(X, Y, '.', 'MarkerSize', 0.1);
    hold on;
    
    % plot nom state
plot(nomStates(1), nomStates(2), 'o', ...
    'Color', [1, 0.5, 0], ...        % Orange RGB
    'MarkerFaceColor', [1, 0.5, 0], ... % Fill the dot
    'MarkerSize', 4);    

    % Plot mean point
    plot(mean(X), mean(Y), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Plot covariance ellipse
    plotCovarianceEllipse(X, Y);
    
    title(titleText);
    hold off;
end

function plotCovarianceEllipse(X, Y)
    % Helper function to plot 3-sigma covariance ellipse
    % Compute covariance matrix
    C = cov(X', Y');
    
    % Generate unit circle
    theta = linspace(0, 2*pi, 100);
    circle = [cos(theta); sin(theta)];
    
    % SVD for ellipse orientation
    [U, S, ~] = svd(C);
    
    % Scale and rotate the unit circle
    ellipse = 3 * U * sqrt(S) * circle;
    
    % Shift ellipse to the mean position
    ellipse = ellipse + [mean(X); mean(Y)];
    
    % Plot the ellipse
    plot(ellipse(1, :), ellipse(2, :), 'r-', 'LineWidth', 1);
end
