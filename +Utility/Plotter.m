function Plotter(x, y, plotTitle, xLabel, yLabel)
    % createPlot: A function to generate a plot from input data.
    % Input:
    %   x - vector of x values
    %   y - vector of y values
    %   plotTitle - title of the plot (string)
    %   xLabel - label for the x-axis (string)
    %   yLabel - label for the y-axis (string)
    
    % Check if x and y have the same length
    if length(x) ~= length(y)
        error('x and y must have the same length');
    end
    
    % Create the plot
    figure;  % Open a new figure window
    plot(x, y, 'b-', 'LineWidth', 2);  % Plot x vs y in blue with line width 2
    
    % Customize the plot
    title(plotTitle);  % Set the plot title
    xlabel(xLabel);    % Set the x-axis label
    ylabel(yLabel);    % Set the y-axis label
    grid on;           % Enable grid
    axis tight;        % Adjust axes to fit the data
end
