function [matrix] = construct_laplace_matrix(N)

% this is the function that generate the discrete laplacian matrix for a

% grid of size N, with dirichlet boundary conditons at both extremities of the grid
matrix = zeros(N,N);
matrix(1,1) =1;
matrix(N,N) =1;
for i=2:N-1    
    matrix(i,i-1:i+1) = [-1,2,-1]*N^2;
    %% add the coefficients
end
end