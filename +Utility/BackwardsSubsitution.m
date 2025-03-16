function [x] = BackwardsSubsitution(U, b)
% Function solves backward subsitution problem 
% can be used instead of matrix inverse! Assumes input is upper triangular

    n = length(b);      % Number of equations
    x = zeros(n, 1);    % Initialize solution vector
    
    % Perform backward substitution
    for i = n:-1:1
        % Calculate the i-th variable
        x(i) = (b(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end