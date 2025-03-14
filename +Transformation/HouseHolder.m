function [TA] = HouseHolder(A)
%
% This function computes the householder transformation for a given matrix.
% The goal is to upper triangularize an (n + m) x (n + 1) matrix 
%
%%%%%%%%% Inputs %%%%%%%%
% 
%   
%   A = [Rbar bbar; H y]
%


% Ensure given matrix A is correct size
[numRow, numCol] = size(A); 


% Implement algorithm 
for k = 1:numCol-1
    
    % set sigma value
    sigma = sign(A(k,k)) * sqrt(sum(A(k:numRow, k).^2));
    
    % calculate uk
    u(k) = A(k,k) + sigma; 
    
    % set value of A(k,k) to be sigma
    A(k,k) = -sigma; 
    
    % get ui value 
    u(k+1:numRow) = A(k+1:numRow, k);
    
    % set Beta 
    Beta = 1/(sigma*u(k));
    
    % nested loop 
    for j = k+1:numCol
        
        % get Gamma 
        Gamma = Beta * sum(u(k:numRow)'.*A(k:numRow,j));
        
        % fill in values for A
        A(k:numRow,j) = A(k:numRow,j) - Gamma * u(k:numRow)';
    end
    
    % set ik value for A
    A(k+1:numRow,k) = 0; % probably zeros
    
end
    
% Transformed A is output
TA = A;
