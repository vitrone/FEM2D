function [xi, quadW] = fem2d_Gdata(degree)

% Gauss quadrature data on the reference triangular element (K_t).
% xi    - Gauss quadrature nodes (in K_t)
% quadW - associated weights
%
% Reference triangle (K_t):
%
%                 xi_2
%                 ^  
%     (-1,1)      |
%       +         +--> xi_1 
%       |\
%       | \
%       |  \
%       |   \
%       +----+
% (-1,-1)     (1,-1)
% 
% 

numnodes = rule_full_size ( degree );
vert = [-1, -1;
         1, -1;
        -1,  1];
[ xi, quadW ] = triasymq( degree, vert(1, :), vert(2, :), vert(3, :),  numnodes);

