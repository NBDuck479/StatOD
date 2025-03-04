function LKFhistogramPlotter(measResHist)
% LKFhistogramPlotter generates histograms for each row in measResHist.

% Number of rows in measResHist
numRows = size(measResHist, 1);

% Create a figure for the histograms
figure;

% Loop over each row
for i = 1:numRows
    % Extract the data for the current row
    data = measResHist(i, :);
    
    % Create a subplot for each row's histogram
    subplot(numRows, 1, i);
    
    % For the first subplot, set the bin width to 0.0001
    if i == 1
        histogram(data);
        % Add labels and title
        title('Range');
        ylabel('Count');
        grid on;
    else
        histogram(data);  % Default bin width for other subplots
        % Add labels and title
        title('Range Rate');
        xlabel('Residuals');
        ylabel('Count');
        grid on;
        
        sgtitle('LKF Post-Fit Residuals Histogram')
    end
    
end
end
